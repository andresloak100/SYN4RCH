"""Curation command line entry point.

    python3 -m syn4rch.curate.run units
    python3 -m syn4rch.curate.run gate PROPOSALS.jsonl

``units`` (stages 4-5 of docs/CURATION_PROTOCOL.md)
    Reads ``data/interim/candidates.jsonl`` and ``data/raw/fulltext/*.xml``,
    splits every retrieved paper into passages, and writes one row per
    localized (population, method, outcome) co-occurrence to
    ``data/interim/extraction_units.jsonl``. Those rows are *candidates for a
    reader*, not records; see :mod:`syn4rch.curate.units`.

``gate`` (stages 7-8)
    Reads a proposed-records JSONL, runs :func:`syn4rch.schema.validate.validate`
    and the verbatim gate :func:`syn4rch.curate.verify.gate`, writes admitted
    records to ``data/curated/records.jsonl`` and every rejection to
    ``reports/curation_rejects.md``. Per DECISIONS D-017 this is the only
    writer of the curated file.

Exit status:

    0  nothing was rejected and every input line was readable
    1  at least one proposal was rejected, or a line could not be read
    2  usage error (missing input file, bad arguments)

A non-zero ``gate`` status is the normal outcome of a reader that proposed
something unsupported, and ``reports/curation_rejects.md`` is the point of it.
Rejections are counted and published rather than suppressed
(``docs/CURATION_PROTOCOL.md`` 4).

**Every number this CLI prints is counted, never estimated.** Where a quantity
is not known (how many papers exist, how many records a corpus could yield) it
is not printed.

Python 3.10 compatible, standard library only. No network.
"""

from __future__ import annotations

import argparse
import datetime
import glob
import json
import os
import sys
from typing import Dict, List, Optional, Tuple

from ..schema.io import JsonlError, loads_record, write_jsonl
from ..schema.record import SCHEMA_VERSION, TierARecord
from ..schema.validate import Severity, validate
from .jats import ParsedPaper, parse_jats_file
from .passages import split_passages
from .units import paper_term_summary, units_from_paper
from .verify import CODES, NORMALISATION, Rejection, gate as gate_records, index_papers

__all__ = ["main", "build_parser", "cmd_units", "cmd_gate"]


# --------------------------------------------------------------------------- #
# paths
# --------------------------------------------------------------------------- #


def _repo_root() -> str:
    try:
        from ..harvest.cache import repo_root

        return repo_root()
    except Exception:
        return os.path.abspath(os.getcwd())


def _default(*parts: str) -> str:
    return os.path.join(_repo_root(), *parts)


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# --------------------------------------------------------------------------- #
# reading inputs
# --------------------------------------------------------------------------- #


def read_candidates(path: str) -> Tuple[List[dict], int]:
    """Read ``candidates.jsonl`` tolerantly. Returns ``(rows, unreadable)``.

    ``candidates.jsonl`` is screening output (plain JSON objects), not Tier A
    records, so it is read with :mod:`json` and not through
    :mod:`syn4rch.schema.io`. A missing file is not an error: the harvest may
    not have written it yet, and the full-text directory is the authority on
    what text exists.
    """
    rows: List[dict] = []
    unreadable = 0
    if not os.path.isfile(path):
        return rows, unreadable
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except ValueError:
                unreadable += 1
                continue
            if isinstance(row, dict):
                rows.append(row)
            else:
                unreadable += 1
    return rows, unreadable


def fulltext_files(directory: str) -> List[str]:
    """Every ``*.xml`` in the full-text directory, sorted.

    ``data/raw/fulltext/`` is append-only and may be being written right now,
    so the listing is a snapshot and a file in it may be incomplete. That is
    handled by :func:`syn4rch.curate.jats.parse_jats_file`, which never
    raises; nothing here retries or waits.
    """
    if not os.path.isdir(directory):
        return []
    return sorted(glob.glob(os.path.join(directory, "*.xml")))


def _candidate_index(rows: List[dict]) -> Dict[str, dict]:
    index: Dict[str, dict] = {}
    for row in rows:
        pmcid = row.get("pmcid")
        if pmcid:
            index.setdefault(str(pmcid).strip().upper(), row)
    return index


def _study_id_for(paper: ParsedPaper, candidate: Optional[dict]) -> Optional[str]:
    """``study_id`` for a paper: the candidate's if it has one, else the
    paper's own ids in the order SCHEMA.md 3.1 prefers (``doi:`` then
    ``pmid:``). Never invented."""
    if candidate:
        existing = candidate.get("study_id")
        if existing:
            return str(existing)
        doi = candidate.get("doi")
        if doi:
            return "doi:%s" % str(doi).lower()
        pmid = candidate.get("pmid")
        if pmid:
            return "pmid:%s" % pmid
    if paper.doi:
        return "doi:%s" % paper.doi
    if paper.pmid:
        return "pmid:%s" % paper.pmid
    if paper.pmcid:
        return "pmcid:%s" % paper.pmcid
    return None


# --------------------------------------------------------------------------- #
# units
# --------------------------------------------------------------------------- #


def cmd_units(args) -> int:
    candidates, bad_candidate_lines = read_candidates(args.candidates)
    by_pmcid = _candidate_index(candidates)
    paths = fulltext_files(args.fulltext_dir)

    counts = {
        "candidate_rows_read": len(candidates),
        "candidate_lines_unreadable": bad_candidate_lines,
        "candidates_with_pmcid": len(by_pmcid),
        "fulltext_files_found": len(paths),
        "papers_matched_to_a_candidate": 0,
        "papers_not_in_candidates": 0,
        "papers_with_parse_warnings": 0,
        "papers_with_no_readable_passages": 0,
        "papers_without_abstract": 0,
        "papers_with_no_study_id": 0,
        "passages_total": 0,
        "papers_with_at_least_one_unit": 0,
        "papers_mentioning_all_three_but_yielding_no_unit": 0,
        "units_written": 0,
        "units_dropped_by_per_passage_cap": 0,
    }
    by_section: Dict[str, int] = {}
    rows: List[dict] = []

    for path in paths:
        paper = parse_jats_file(path)
        pmcid = (paper.pmcid or "").upper()
        candidate = by_pmcid.get(pmcid) if pmcid else None
        if candidate is not None:
            counts["papers_matched_to_a_candidate"] += 1
        else:
            counts["papers_not_in_candidates"] += 1
        if paper.parse_warnings:
            counts["papers_with_parse_warnings"] += 1
        if any("no abstract found" in w for w in paper.parse_warnings):
            counts["papers_without_abstract"] += 1
        passages = split_passages(paper)
        counts["passages_total"] += len(passages)
        if not passages:
            counts["papers_with_no_readable_passages"] += 1
            continue

        study_id = _study_id_for(paper, candidate)
        if study_id is None:
            counts["papers_with_no_study_id"] += 1

        uncapped = units_from_paper(paper, study_id=study_id, passages=passages)
        units = (
            units_from_paper(
                paper,
                study_id=study_id,
                passages=passages,
                max_per_passage=args.max_per_passage,
            )
            if args.max_per_passage > 0
            else uncapped
        )
        counts["units_dropped_by_per_passage_cap"] += len(uncapped) - len(units)
        if units:
            counts["papers_with_at_least_one_unit"] += 1
        else:
            mentions = paper_term_summary(paper, passages=passages)
            if mentions["populations"] and mentions["methods"] and mentions["outcomes"]:
                counts["papers_mentioning_all_three_but_yielding_no_unit"] += 1
        for unit in units:
            by_section[unit.section_label] = by_section.get(unit.section_label, 0) + 1
            rows.append(unit.to_dict())

    rows.sort(
        key=lambda row: (
            -int(row["score"]),
            str(row["study_id"] or ""),
            row["locator"],
            row["candidate_population_id"],
            row["candidate_method_id"],
            row["candidate_outcome_id"],
        )
    )
    counts["units_written"] = len(rows)

    if not args.dry_run:
        parent = os.path.dirname(os.path.abspath(args.out))
        if parent:
            os.makedirs(parent, exist_ok=True)
        with open(args.out, "w", encoding="utf-8", newline="\n") as handle:
            for row in rows:
                handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True))
                handle.write("\n")

    print("stage 4-5: passage segmentation and unit building")
    print("  candidates:   %s" % args.candidates)
    print("  full text:    %s" % args.fulltext_dir)
    print("  units out:    %s%s" % (args.out, " (not written: --dry-run)" if args.dry_run else ""))
    print("")
    for key in (
        "candidate_rows_read",
        "candidate_lines_unreadable",
        "candidates_with_pmcid",
        "fulltext_files_found",
        "papers_matched_to_a_candidate",
        "papers_not_in_candidates",
        "papers_with_parse_warnings",
        "papers_with_no_readable_passages",
        "papers_without_abstract",
        "papers_with_no_study_id",
        "passages_total",
        "papers_with_at_least_one_unit",
        "papers_mentioning_all_three_but_yielding_no_unit",
        "units_dropped_by_per_passage_cap",
        "units_written",
    ):
        print("  %-48s %d" % (key, counts[key]))
    if by_section:
        print("")
        print("  units by section label (counted):")
        for label, count in sorted(by_section.items(), key=lambda kv: (-kv[1], kv[0])):
            print("    %-20s %d" % (label, count))
    print("")
    print("  A unit is a place worth reading, not a record. Expect a low")
    print("  units-to-records yield; that is the protocol working.")
    return 0


# --------------------------------------------------------------------------- #
# gate
# --------------------------------------------------------------------------- #


def read_proposals(
    path: str, expect_version: Optional[str] = SCHEMA_VERSION
) -> Tuple[List[TierARecord], List[Tuple[int, str]]]:
    """Read proposed records. Returns ``(records, [(line_number, error)])``.

    A line that cannot be read is collected rather than raised on, so that one
    malformed proposal does not stop the other thousand from being gated. Every
    such line is reported and makes the run exit non-zero.
    """
    records: List[TierARecord] = []
    errors: List[Tuple[int, str]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                records.append(
                    loads_record(
                        line,
                        expect_version=expect_version,
                        path=path,
                        line_number=number,
                    )
                )
            except JsonlError as exc:
                errors.append((number, str(exc)))
            except Exception as exc:  # pragma: no cover - defensive
                errors.append((number, "%s: %s" % (type(exc).__name__, exc)))
    return records, errors


def load_papers(directory: str) -> Tuple[List[ParsedPaper], Dict[str, int]]:
    """Parse every full-text file in ``directory``. Never raises."""
    papers: List[ParsedPaper] = []
    stats = {"files": 0, "with_parse_warnings": 0, "with_no_passages": 0}
    for path in fulltext_files(directory):
        stats["files"] += 1
        paper = parse_jats_file(path)
        if paper.parse_warnings:
            stats["with_parse_warnings"] += 1
        if paper.is_empty:
            stats["with_no_passages"] += 1
        papers.append(paper)
    return papers, stats


def render_rejects_md(
    rejections: List[Rejection],
    admitted: List[TierARecord],
    *,
    proposals_path: str,
    fulltext_dir: str,
    paper_stats: Dict[str, int],
    read_errors: List[Tuple[int, str]],
    strict_locator_passage: bool,
    note_counts: Dict[str, int],
    started: str,
    excerpt_chars: int = 400,
) -> str:
    """The rejection report. Every number in it is counted from the inputs."""
    total = len(admitted) + len(rejections)
    L: List[str] = []
    L.append("# Curation gate: rejections")
    L.append("")
    L.append("Generated %s by `python3 -m syn4rch.curate.run gate`." % started)
    L.append("")
    L.append("| input | value |")
    L.append("|---|---|")
    L.append("| proposals | `%s` |" % proposals_path)
    L.append("| retrieved full text | `%s` |" % fulltext_dir)
    L.append("| full-text files parsed | %d |" % paper_stats.get("files", 0))
    L.append(
        "| of those, with parse warnings | %d |"
        % paper_stats.get("with_parse_warnings", 0)
    )
    L.append(
        "| of those, yielding no passages | %d |"
        % paper_stats.get("with_no_passages", 0)
    )
    L.append("| proposals read | %d |" % total)
    L.append("| proposal lines unreadable | %d |" % len(read_errors))
    L.append("| admitted | %d |" % len(admitted))
    L.append("| rejected | %d |" % len(rejections))
    L.append(
        "| locator-passage agreement | %s |"
        % ("blocking (default)" if strict_locator_passage else "reported only (--lenient-passage)")
    )
    L.append("")
    L.append("## The check")
    L.append("")
    L.append("`support_quote` must occur as a substring of the paper's retrieved")
    L.append("text under exactly one normalisation: %s" % NORMALISATION)
    L.append("")
    L.append("`source_locator` must name a passage that exists in that paper.")
    L.append("")
    L.append("A failed proposal is rejected and recorded. It is never repaired,")
    L.append("re-quoted or nudged into passing (docs/CURATION_PROTOCOL.md 4).")
    L.append("")

    L.append("## Rejections by reason")
    L.append("")
    if rejections:
        by_code: Dict[str, int] = {}
        for rejection in rejections:
            for code in rejection.codes:
                by_code[code] = by_code.get(code, 0) + 1
        L.append("| code | rejections | meaning |")
        L.append("|---|---|---|")
        for code, count in sorted(by_code.items(), key=lambda kv: (-kv[1], kv[0])):
            meaning = CODES.get(code, (True, "schema validation finding"))[1]
            L.append("| `%s` | %d | %s |" % (code, count, meaning.replace("|", "/")))
    else:
        L.append("None. Every proposal read was admitted.")
    L.append("")

    if note_counts:
        L.append("## Non-blocking notes on records that were admitted")
        L.append("")
        L.append("These did not block admission. They are published because a")
        L.append("silently-held caveat is not a caveat.")
        L.append("")
        L.append("| code | records | meaning |")
        L.append("|---|---|---|")
        for code, count in sorted(note_counts.items(), key=lambda kv: (-kv[1], kv[0])):
            meaning = CODES.get(code, (False, ""))[1]
            L.append("| `%s` | %d | %s |" % (code, count, meaning.replace("|", "/")))
        L.append("")

    if read_errors:
        L.append("## Unreadable proposal lines")
        L.append("")
        for number, message in read_errors:
            L.append("* line %d: %s" % (number, message))
        L.append("")

    L.append("## Every rejection")
    L.append("")
    if not rejections:
        L.append("None.")
        L.append("")
    for position, rejection in enumerate(rejections, start=1):
        record = rejection.record
        L.append(
            "### %d. `%s` - %s"
            % (position, rejection.record_id or "(no record_id)", rejection.code)
        )
        L.append("")
        L.append("| field | value |")
        L.append("|---|---|")
        L.append("| study_id | `%s` |" % (record.study_id or "null"))
        L.append("| source_locator | `%s` |" % (record.source_locator or "null"))
        L.append("| population_raw | %s |" % (record.population_raw or "null"))
        L.append("| outcome_raw | %s |" % (record.outcome_raw or "null"))
        L.append(
            "| direction | %s |"
            % (record.direction.value if record.direction else "null")
        )
        L.append(
            "| extraction_method | %s |"
            % (record.extraction_method.value if record.extraction_method else "null")
        )
        verification = rejection.verification
        if verification is not None:
            L.append("| paper matched by | `%s` |" % (verification.paper_key or "none"))
            L.append("| locator resolution | %s |" % verification.locator_resolution)
        L.append("")
        for code, reason in zip(rejection.codes, rejection.reasons):
            L.append("* **%s** %s" % (code, reason))
        if record.support_quote:
            L.append("")
            L.append("Quote as proposed (not found; shown verbatim, not repaired):")
            L.append("")
            L.append("> %s" % str(record.support_quote)[:excerpt_chars].replace("\n", " "))
        if verification is not None and verification.located_passages:
            L.append("")
            L.append(
                "Locator resolved to: %s"
                % ", ".join("`%s`" % loc for loc in verification.located_passages[:5])
            )
        if verification is not None and verification.quote_found_in:
            L.append("")
            L.append(
                "The quote does occur in: %s"
                % ", ".join("`%s`" % loc for loc in verification.quote_found_in[:5])
            )
        L.append("")
    return "\n".join(L) + "\n"


def cmd_gate(args) -> int:
    started = _now()
    if not os.path.isfile(args.proposals):
        sys.stderr.write("error: no such file: %s\n" % args.proposals)
        return 2

    records, read_errors = read_proposals(
        args.proposals, expect_version=None if args.any_version else SCHEMA_VERSION
    )
    papers, paper_stats = load_papers(args.fulltext_dir)
    index = index_papers(papers)

    admitted, rejections = gate_records(
        records,
        index,
        strict_locator_passage=(not args.lenient_passage),
    )

    # Non-blocking notes on admitted records, recomputed from the same gate so
    # the report cannot disagree with the verdict.
    from .verify import verify_support

    note_counts: Dict[str, int] = {}
    warning_counts: Dict[str, int] = {}
    for record in admitted:
        paper_index, _key = _lookup(record, index)
        verification = verify_support(
            record, paper_index, strict_locator_passage=(not args.lenient_passage)
        )
        for finding in verification.notes:
            note_counts[finding.code] = note_counts.get(finding.code, 0) + 1
        for finding in validate(record):
            if finding.severity is Severity.WARNING:
                warning_counts[finding.code] = warning_counts.get(finding.code, 0) + 1

    report = render_rejects_md(
        rejections,
        admitted,
        proposals_path=args.proposals,
        fulltext_dir=args.fulltext_dir,
        paper_stats=paper_stats,
        read_errors=read_errors,
        strict_locator_passage=(not args.lenient_passage),
        note_counts=note_counts,
        started=started,
    )

    if not args.dry_run:
        parent = os.path.dirname(os.path.abspath(args.report))
        if parent:
            os.makedirs(parent, exist_ok=True)
        with open(args.report, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(report)
        write_jsonl(args.out, admitted)

    print("stage 7-8: verbatim gate and schema validation")
    print("  proposals:    %s" % args.proposals)
    print("  full text:    %s" % args.fulltext_dir)
    print("  curated out:  %s%s" % (args.out, " (not written: --dry-run)" if args.dry_run else ""))
    print("  reject report:%s%s" % (" " + args.report, " (not written: --dry-run)" if args.dry_run else ""))
    print("")
    print("  %-48s %d" % ("full_text_files_parsed", paper_stats.get("files", 0)))
    print("  %-48s %d" % ("full_text_files_with_parse_warnings", paper_stats.get("with_parse_warnings", 0)))
    print("  %-48s %d" % ("full_text_files_with_no_passages", paper_stats.get("with_no_passages", 0)))
    print("  %-48s %d" % ("papers_indexed_by_identifier", len(index)))
    print("  %-48s %d" % ("proposals_read", len(records)))
    print("  %-48s %d" % ("proposal_lines_unreadable", len(read_errors)))
    print("  %-48s %d" % ("admitted", len(admitted)))
    print("  %-48s %d" % ("rejected", len(rejections)))
    if rejections:
        print("")
        print("  rejections by code (counted):")
        by_code: Dict[str, int] = {}
        for rejection in rejections:
            for code in rejection.codes:
                by_code[code] = by_code.get(code, 0) + 1
        for code, count in sorted(by_code.items(), key=lambda kv: (-kv[1], kv[0])):
            print("    %-44s %d" % (code, count))
    if note_counts:
        print("")
        print("  non-blocking notes on admitted records (counted):")
        for code, count in sorted(note_counts.items(), key=lambda kv: (-kv[1], kv[0])):
            print("    %-44s %d" % (code, count))
    if warning_counts:
        print("")
        print("  schema warnings on admitted records (counted):")
        for code, count in sorted(warning_counts.items(), key=lambda kv: (-kv[1], kv[0])):
            print("    %-44s %d" % (code, count))
    for number, message in read_errors:
        sys.stderr.write("unreadable proposal line %d: %s\n" % (number, message))
    return 1 if (rejections or read_errors) else 0


def _lookup(record: TierARecord, index):
    from .verify import find_paper

    return find_paper(record, index)


# --------------------------------------------------------------------------- #
# argument parsing
# --------------------------------------------------------------------------- #


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python3 -m syn4rch.curate.run",
        description=(
            "Turn retrieved open-access full text into extraction unit "
            "candidates, and gate proposed Tier A records against the text "
            "they claim to come from."
        ),
    )
    sub = parser.add_subparsers(dest="command")

    units = sub.add_parser(
        "units",
        help="segment retrieved full text and write extraction unit candidates",
        description=(
            "Stages 4-5. Units are places worth reading, not records "
            "(docs/CURATION_PROTOCOL.md 2)."
        ),
    )
    units.add_argument(
        "--candidates",
        default=None,
        metavar="PATH",
        help="screening output (default: data/interim/candidates.jsonl)",
    )
    units.add_argument(
        "--fulltext-dir",
        default=None,
        metavar="DIR",
        help="retrieved JATS XML (default: data/raw/fulltext)",
    )
    units.add_argument(
        "--out",
        default=None,
        metavar="PATH",
        help="where to write units (default: data/interim/extraction_units.jsonl)",
    )
    units.add_argument(
        "--max-per-passage",
        type=int,
        default=0,
        metavar="N",
        help="keep at most N highest-scoring units per passage (0 = no cap); "
        "the number dropped is reported",
    )
    units.add_argument(
        "--dry-run",
        action="store_true",
        help="count and print, write nothing",
    )
    units.set_defaults(func=cmd_units)

    gate = sub.add_parser(
        "gate",
        help="verify proposed records against the retrieved text and admit or reject",
        description=(
            "Stages 7-8. The verdict is mechanical and final: a failed "
            "proposal is rejected and recorded, never repaired."
        ),
    )
    gate.add_argument("proposals", help="path to a proposed-records JSONL")
    gate.add_argument(
        "--fulltext-dir",
        default=None,
        metavar="DIR",
        help="retrieved JATS XML (default: data/raw/fulltext)",
    )
    gate.add_argument(
        "--out",
        default=None,
        metavar="PATH",
        help="where admitted records go (default: data/curated/records.jsonl)",
    )
    gate.add_argument(
        "--report",
        default=None,
        metavar="PATH",
        help="rejection report (default: reports/curation_rejects.md)",
    )
    gate.add_argument(
        "--lenient-passage",
        action="store_true",
        help="do NOT reject a verbatim quote that occurs elsewhere in the paper "
        "than the passage its locator names; report it only. Off by default "
        "(DECISIONS.md D-018): use only for hand-curated records whose "
        "paragraph count may differ from this parser's",
    )
    gate.add_argument(
        "--any-version",
        action="store_true",
        help="read proposals whose schema_version differs from this code's",
    )
    gate.add_argument(
        "--dry-run",
        action="store_true",
        help="count and print, write nothing",
    )
    gate.set_defaults(func=cmd_gate)
    return parser


def _apply_defaults(args) -> None:
    if args.command == "units":
        if args.candidates is None:
            args.candidates = _default("data", "interim", "candidates.jsonl")
        if args.fulltext_dir is None:
            args.fulltext_dir = _default("data", "raw", "fulltext")
        if args.out is None:
            args.out = _default("data", "interim", "extraction_units.jsonl")
    elif args.command == "gate":
        if args.fulltext_dir is None:
            args.fulltext_dir = _default("data", "raw", "fulltext")
        if args.out is None:
            args.out = _default("data", "curated", "records.jsonl")
        if args.report is None:
            args.report = _default("reports", "curation_rejects.md")


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        return 2
    _apply_defaults(args)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
