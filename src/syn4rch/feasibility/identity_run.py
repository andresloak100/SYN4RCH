"""Perturbation identity command line entry point.

    python3 -m syn4rch.feasibility.identity_run
    python3 -m syn4rch.feasibility.identity_run \\
        --records data/curated/records.jsonl \\
        --out reports/perturbation_identity.md \\
        [--json reports/perturbation_identity.json] [--scope usable]

Defaults read ``data/curated/records.jsonl`` and write
``reports/perturbation_identity.md``.

**Scope.** ``--scope all`` (the default) analyses every record read.
``--scope usable`` restricts to the records that count as usable observations
under :mod:`syn4rch.feasibility.coverage`, which is the subset a model would
actually be fitted on. Both are reported in the written Markdown, because the
two answer different questions: the first is a property of the corpus, the
second a property of the trainable dataset. Neither is a default that hides
the other.

**A missing or empty input is not an error**, and neither is a finding that an
entire perturbation class is unlearnable. Both exit ``0``: SCHEMA-conformant
evidence that the corpus cannot support cross-study prediction for a class of
records is this module's success case, and an exit code that failed on it
would train whoever runs it to stop running it (the same reasoning as
docs/DECISIONS.md D-012).

Exit status:

    0  a report was produced (whatever it concludes)
    1  the input exists but a line could not be read as a Tier A record
    2  usage error (unwritable output path, bad arguments)

This CLI never writes to ``data/``: it reads the curated file and writes only
its own report.

Python 3.10 compatible, standard library only.
"""

import argparse
import json
import os
import sys
from typing import List, Optional

from ..schema.io import JsonlError, iter_jsonl
from ..schema.record import SCHEMA_VERSION
from .coverage import partition_records
from .identity import build_report, render_markdown, top_line

__all__ = ["main", "build_parser", "DEFAULT_RECORDS", "DEFAULT_REPORT", "SCOPES"]

#: Path, relative to the repository root, of the dataset (SCHEMA.md 4).
DEFAULT_RECORDS = os.path.join("data", "curated", "records.jsonl")

#: Path, relative to the repository root, of the Markdown report. Tracked in
#: version control per docs/DECISIONS.md D-016.
DEFAULT_REPORT = os.path.join("reports", "perturbation_identity.md")

#: Which record set the headline analysis runs over.
SCOPES = ("all", "usable")

_SCOPE_LABELS = {
    "all": "all records read",
    "usable": (
        "usable observations only (syn4rch.feasibility.coverage, "
        "schema v0.2.0 predicate)"
    ),
}


def _repo_root() -> str:
    """Best guess at the repository root: ``src/syn4rch/feasibility/`` + 3 up.

    Falls back to the working directory when this file is not laid out that
    way (e.g. an installed wheel), so the defaults still resolve to something
    a caller can see.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    candidate = os.path.dirname(os.path.dirname(os.path.dirname(here)))
    if os.path.isdir(os.path.join(candidate, "docs")) and os.path.isdir(
        os.path.join(candidate, "data")
    ):
        return candidate
    return os.getcwd()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python3 -m syn4rch.feasibility.identity_run",
        description=(
            "Measure, per declared grouping, what fraction of records have a "
            "perturbation identity that occurs in exactly one group and are "
            "therefore unlearnable by construction (SCHEMA.md 5b). An entire "
            "class being unlearnable is a valid and expected conclusion."
        ),
    )
    parser.add_argument(
        "--records",
        default=None,
        metavar="PATH",
        help="input JSONL (default: <repo>/%s). A missing or empty file is "
        "reported as an empty corpus, not an error." % DEFAULT_RECORDS,
    )
    parser.add_argument(
        "--out",
        default=None,
        metavar="PATH",
        help="Markdown report to write (default: <repo>/%s)" % DEFAULT_REPORT,
    )
    parser.add_argument(
        "--json",
        dest="json_path",
        default=None,
        metavar="PATH",
        help="also write the full measurement as JSON",
    )
    parser.add_argument(
        "--scope",
        choices=SCOPES,
        default="all",
        help="which record set the headline analysis runs over: 'all' records "
        "read (default) or only the records usable as observations. The "
        "written report contains both; this chooses which one leads.",
    )
    parser.add_argument(
        "--stdout",
        action="store_true",
        help="write the Markdown report to stdout instead of a file",
    )
    parser.add_argument(
        "--any-version",
        action="store_true",
        help="read records whose schema_version differs from this code's "
        "(migration inspection only; off by default)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="suppress the progress and verdict lines on stdout; the report "
        "file is still written",
    )
    return parser


def _load(path: str, any_version: bool):
    """Return ``(records, load_errors, exists)``; never raises on a bad line."""
    if not os.path.isfile(path):
        return [], [], False
    records = []
    errors: List[str] = []
    try:
        for record in iter_jsonl(
            path, expect_version=None if any_version else SCHEMA_VERSION
        ):
            records.append(record)
    except JsonlError as exc:
        errors.append(str(exc))
    return records, errors, True


def _display_path(path: str, root: str) -> str:
    """The input path as the report should name it.

    Relative to the repository root when it is inside it, so a committed
    report does not carry one machine's home directory, and absolute
    otherwise so it still says where the numbers came from.
    """
    try:
        relative = os.path.relpath(os.path.abspath(path), root)
    except ValueError:  # different drive on Windows
        return path
    if relative.startswith(os.pardir):
        return path
    return relative.replace(os.sep, "/")


def _write_text(path: str, text: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def _render_both(primary, secondary, secondary_scope: str) -> str:
    """The primary report, with the other scope's headline table appended.

    The second scope is not a separate document: a reader who sees one
    unseen-identity fraction should see the other in the same file, or the
    choice of scope becomes a way to pick the more comfortable number.
    """
    text = render_markdown(primary)
    lines = [
        text.rstrip(),
        "",
        "## The other scope, for comparison",
        "",
        "The analysis above ran over **%s**. The same measurement over "
        "**%s** is below, so that the choice of scope cannot be a choice of "
        "number." % (primary.scope, secondary_scope),
        "",
        "| grouping | records | placed | single-group | unseen-identity fraction |",
        "|---|---|---|---|---|",
    ]
    for coverage in secondary.coverages:
        value = coverage.fraction_single_group
        lines.append(
            "| `%s` | %d | %d | %d | %s |"
            % (
                coverage.grouping,
                secondary.n_records,
                coverage.n_records_placed,
                coverage.n_records_single_group,
                "not measurable" if value is None else "%.1f%%" % (100.0 * value),
            )
        )
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    root = _repo_root()
    records_path = args.records or os.path.join(root, DEFAULT_RECORDS)
    out_path = args.out or os.path.join(root, DEFAULT_REPORT)

    records, load_errors, exists = _load(records_path, args.any_version)
    usable, _, _ = partition_records(records)
    source = _display_path(records_path, root)

    scopes = {
        "all": records,
        "usable": usable,
    }
    primary_name = args.scope
    secondary_name = "usable" if primary_name == "all" else "all"

    primary = build_report(
        scopes[primary_name],
        source_path=source,
        source_exists=exists,
        load_errors=load_errors,
        scope=_SCOPE_LABELS[primary_name],
    )
    secondary = build_report(
        scopes[secondary_name],
        source_path=source,
        source_exists=exists,
        load_errors=load_errors,
        scope=_SCOPE_LABELS[secondary_name],
        generated_at=primary.generated_at,
    )
    text = _render_both(primary, secondary, _SCOPE_LABELS[secondary_name])

    if args.stdout:
        sys.stdout.write(text)
    else:
        try:
            _write_text(out_path, text)
        except OSError as exc:
            sys.stderr.write("error: cannot write %s: %s\n" % (out_path, exc))
            return 2
        if not args.quiet:
            sys.stdout.write("wrote %s\n" % out_path)

    if args.json_path:
        payload = {
            "primary_scope": primary_name,
            "primary": primary.to_dict(),
            "secondary_scope": secondary_name,
            "secondary": secondary.to_dict(),
        }
        try:
            parent = os.path.dirname(os.path.abspath(args.json_path))
            if parent:
                os.makedirs(parent, exist_ok=True)
            with open(args.json_path, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2, ensure_ascii=False)
                handle.write("\n")
        except OSError as exc:
            sys.stderr.write("error: cannot write %s: %s\n" % (args.json_path, exc))
            return 2
        if not args.quiet:
            sys.stdout.write("wrote %s\n" % args.json_path)

    if not args.quiet and not args.stdout:
        sys.stdout.write(top_line(primary) + "\n")

    for err in load_errors:
        sys.stderr.write("load error: %s\n" % err)
    return 1 if load_errors else 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
