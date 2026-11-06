"""Fill the ``background_genotype`` dimension from a declared table (D-060).

The genotype dimension arrived with schema v0.4.0 as a *declared* field, null on
all 508 records. This pass fills it, and it is a **fill** rather than a
correction in the precise sense D-039 set: it writes where the value is ``None``
and refuses where it is not. Overwriting a curated value is what
:mod:`syn4rch.curate.correct_conditions` exists for, under a different and
narrower set of locks.

Why a third tool rather than extending one of the two
-----------------------------------------------------
:mod:`syn4rch.curate.fill_conditions` is the v0.2.0 -> v0.3.0 pass and is pinned
to that migration so its published report still recomputes (D-046). It also
reads papers for the four dimensions a Methods section states about a cohort,
and a background genotype is not one of those: deciding whether an allele is the
record's own perturbation or the substrate it was applied to is a per-record
judgement no extractor can make (D-060). ``correct_conditions`` refuses a
dimension with no evidence entry, which is every record here. So this is its own
pass, with its own table and its own refusals.

What it will not do
-------------------
**It decides nothing.** Every value is a row of
``data/curated/background_genotype_assignments.json``, carrying the record, the
value, the verbatim quote, the locator, the source and the reasoning. A record
absent from the table keeps ``None``, which D-060 defines as "the curator has
not addressed genotype here". The table is tracked and reviewable in a diff; it
is not inferred at runtime.

**It validates the grammar rather than trusting it.** A value is one of the
three sentinels or one or more ``gene:effect:scope`` segments whose gene resolves
in ``data/vocab/targets.json``, whose effect and scope are declared enums, and
which are deduplicated and sorted. An unresolvable gene is refused, not stored
with a ``raw:`` prefix: the vocabulary was extended for precisely the genes this
pass needs, so an unknown one means the table and the vocabulary disagree and a
curator has to say which is wrong.

**It attaches evidence, and refuses to attach the wrong kind.** Every value
except ``not_reported`` is a positive claim and gets a
:class:`~syn4rch.schema.record.ConditionEvidence` entry, so it cannot earn W017.
``not_reported`` asserts an absence and an absence cannot be quoted, so it gets
no entry — A009 and A011 enforce that in both directions, and this pass refuses
a ``not_reported`` row that carries a quote rather than silently dropping it.

**It re-reads the papers.** Every quote is checked against the passage its
locator names, the locator-strict check from the fill pass, not the weaker
"occurs somewhere in the paper". Any failure blocks the whole write.

**It re-identifies nothing, and proves it.** ``background_genotype`` is not a
``RECORD_ID_FIELDS`` component (D-060), so no ``record_id`` may move. The pass
asserts that per record and refuses if one does, exactly as
:mod:`syn4rch.curate.migrate_0_4_0` does, because neither has an id map to offer.

Python 3.10 compatible, standard library only.
"""

import argparse
import datetime
import json
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from ..schema.enums import ConditionSource
from ..schema.io import read_jsonl
from ..schema.record import (
    RECORD_ID_FIELDS,
    SCHEMA_VERSION,
    ConditionEvidence,
    TierARecord,
)
from ..schema.validate import Severity, validate
from .fill_conditions import check_evidence_locators
from .migrate import _atomic_write_text, _records_text

__all__ = [
    "DIMENSION",
    "EFFECTS",
    "OUTCOMES",
    "SCOPES",
    "SENTINELS",
    "Assignment",
    "FillOutcome",
    "FillResult",
    "GrammarError",
    "IdentityMovedError",
    "apply_assignments",
    "build_parser",
    "default_table_path",
    "fill_rule",
    "main",
    "parse_value",
    "plan_assignments",
    "read_table",
    "render_report_md",
]

DIMENSION = "background_genotype"

#: Values that are not segments. ``null`` is the absence of a row, not a value.
SENTINELS: Tuple[str, ...] = ("not_reported", "wild_type", "mixed")

#: The one sentinel that asserts an absence, and so carries no evidence.
NO_EVIDENCE: Tuple[str, ...] = ("not_reported",)

EFFECTS: Tuple[str, ...] = ("loss_of_function", "gain_of_function", "unclear")
SCOPES: Tuple[str, ...] = ("global", "region_restricted", "cell_type_restricted")

#: Every outcome a row can have, as data, so the report cannot grow a state the
#: tests do not know about.
OUTCOMES: Tuple[str, ...] = (
    "filled",
    "already_filled",
    "refused_record_not_found",
    "refused_value_already_set",
    "refused_bad_grammar",
    "refused_unknown_gene",
    "refused_missing_quote",
    "refused_quote_on_an_absence",
    "refused_duplicate_row",
)


class GrammarError(ValueError):
    """A declared value is not admissible for this dimension."""


class IdentityMovedError(RuntimeError):
    """A ``record_id`` moved, which this pass has no map to explain."""


def fill_rule(decision: str) -> str:
    """The ``rule`` written onto an evidence entry.

    Names the decision that authorised the reading rather than claiming an
    extractor produced it, the same convention
    :mod:`syn4rch.curate.correct_conditions` uses.
    """
    return "read_under_the_arm_rule_%s" % decision.lower().replace("-", "")


def parse_value(value: str, known_genes: Sequence[str]) -> List[Tuple[str, str, str]]:
    """Validate a declared value; return its segments, empty for a sentinel.

    Raises :class:`GrammarError` with the reason. The grammar is checked rather
    than pattern-matched loosely, because a value that is nearly right is the
    one that enters a key and is never looked at again.
    """
    if not isinstance(value, str) or not value.strip():
        raise GrammarError("value is empty")
    value = value.strip()
    if value in SENTINELS:
        return []
    if value == "null":
        raise GrammarError(
            "'null' is the absence of a row, not a value: omit the record from "
            "the table instead"
        )
    segments: List[Tuple[str, str, str]] = []
    seen = set()
    for part in value.split("+"):
        bits = part.split(":")
        if len(bits) != 3:
            raise GrammarError(
                "segment %r is not gene:effect:scope" % part
            )
        gene, effect, scope = (b.strip() for b in bits)
        if gene != gene.lower() or not gene:
            raise GrammarError("gene %r must be lowercase and non-empty" % gene)
        if gene not in known_genes:
            raise GrammarError(
                "gene %r does not resolve in data/vocab/targets.json. The "
                "vocabulary was extended for the genes this pass needs, so an "
                "unknown one means the table and the vocabulary disagree."
                % gene
            )
        if effect not in EFFECTS:
            raise GrammarError(
                "effect %r is not one of %s" % (effect, ", ".join(EFFECTS))
            )
        if scope not in SCOPES:
            raise GrammarError(
                "scope %r is not one of %s" % (scope, ", ".join(SCOPES))
            )
        key = (gene, effect, scope)
        if key in seen:
            raise GrammarError("segment %r appears twice" % part)
        seen.add(key)
        segments.append(key)
    canonical = "+".join(sorted("%s:%s:%s" % s for s in segments))
    if canonical != value:
        raise GrammarError(
            "segments must be sorted and deduplicated: %r should be %r"
            % (value, canonical)
        )
    return segments


@dataclass(frozen=True)
class Assignment:
    """One declared value for one record."""

    record_id: str
    study_id: str
    value: str
    quote: str = ""
    locator: str = ""
    source: str = "passage"
    confidence: Optional[str] = None
    reasoning: str = ""
    rule: Optional[str] = None

    @property
    def needs_evidence(self) -> bool:
        return self.value not in NO_EVIDENCE

    def to_dict(self) -> Dict[str, object]:
        return {
            "record_id": self.record_id,
            "study_id": self.study_id,
            "value": self.value,
            "quote": self.quote,
            "locator": self.locator,
            "source": self.source,
            "confidence": self.confidence,
            "reasoning": self.reasoning,
        }


@dataclass
class FillOutcome:
    """What happened to one assignment, and why."""

    assignment: Assignment
    outcome: str
    detail: Optional[str] = None
    stored_before: Optional[str] = None

    @property
    def refused(self) -> bool:
        return self.outcome.startswith("refused_")

    def to_dict(self) -> Dict[str, object]:
        return {
            "assignment": self.assignment.to_dict(),
            "outcome": self.outcome,
            "detail": self.detail,
            "stored_before": self.stored_before,
        }


@dataclass
class FillResult:
    """Every outcome, the rewritten records, and the audit of the pass."""

    records: List[TierARecord] = field(default_factory=list)
    outcomes: List[FillOutcome] = field(default_factory=list)
    ids_unchanged: int = 0
    evidence_added: int = 0
    locators_checked: int = 0
    locator_failures: List[dict] = field(default_factory=list)
    hard_findings_before: Dict[str, int] = field(default_factory=dict)
    hard_findings_after: Dict[str, int] = field(default_factory=dict)

    @property
    def filled(self) -> List[FillOutcome]:
        return [o for o in self.outcomes if o.outcome == "filled"]

    @property
    def refused(self) -> List[FillOutcome]:
        return [o for o in self.outcomes if o.refused]

    def counts(self) -> Dict[str, int]:
        out = {name: 0 for name in OUTCOMES}
        for item in self.outcomes:
            out[item.outcome] = out.get(item.outcome, 0) + 1
        return out

    def to_dict(self) -> Dict[str, object]:
        return {
            "dimension": DIMENSION,
            "schema_version": SCHEMA_VERSION,
            "n_records": len(self.records),
            "counts": self.counts(),
            "n_filled": len(self.filled),
            "n_refused": len(self.refused),
            "ids_unchanged": self.ids_unchanged,
            "record_id_map_required": False,
            "evidence_entries_added": self.evidence_added,
            "locators_checked": self.locators_checked,
            "locator_failures": list(self.locator_failures),
            "hard_findings_before": dict(sorted(self.hard_findings_before.items())),
            "hard_findings_after": dict(sorted(self.hard_findings_after.items())),
            "outcomes": [o.to_dict() for o in self.outcomes],
        }


def read_table(path: str) -> Tuple[List[Assignment], Dict[str, object]]:
    """Read the declared assignments table."""
    with open(path, "r", encoding="utf-8") as handle:
        doc = json.load(handle)
    meta = doc.get("_meta") or {}
    rows = doc.get("assignments") or []
    out: List[Assignment] = []
    for row in rows:
        out.append(
            Assignment(
                record_id=str(row.get("record_id") or ""),
                study_id=str(row.get("study_id") or ""),
                value=str(row.get("value") or ""),
                quote=str(row.get("quote") or ""),
                locator=str(row.get("locator") or ""),
                source=str(row.get("source") or "passage"),
                confidence=row.get("confidence"),
                reasoning=str(row.get("reasoning") or ""),
                rule=row.get("rule"),
            )
        )
    return out, meta


def _known_genes(path: Optional[str] = None) -> List[str]:
    if path is None:
        path = os.path.join(_repo_root(), "data", "vocab", "targets.json")
    with open(path, "r", encoding="utf-8") as handle:
        return sorted(json.load(handle).get("entries") or {})


def plan_assignments(
    records: Sequence[TierARecord],
    assignments: Sequence[Assignment],
    known_genes: Sequence[str],
) -> List[FillOutcome]:
    """Decide every assignment's outcome without changing anything."""
    index = {r.record_id: r for r in records if r.record_id}
    seen: Dict[str, int] = {}
    out: List[FillOutcome] = []
    for item in assignments:
        seen[item.record_id] = seen.get(item.record_id, 0) + 1
        if seen[item.record_id] > 1:
            out.append(
                FillOutcome(item, "refused_duplicate_row",
                            "this table already assigns %s" % item.record_id)
            )
            continue
        record = index.get(item.record_id)
        if record is None:
            out.append(FillOutcome(item, "refused_record_not_found"))
            continue
        stored = getattr(record, DIMENSION, None)
        stored = getattr(stored, "value", stored)
        try:
            parse_value(item.value, known_genes)
        except GrammarError as exc:
            kind = (
                "refused_unknown_gene"
                if "does not resolve" in str(exc)
                else "refused_bad_grammar"
            )
            out.append(FillOutcome(item, kind, str(exc), stored))
            continue
        if stored == item.value:
            out.append(FillOutcome(item, "already_filled", None, stored))
            continue
        if stored is not None:
            out.append(
                FillOutcome(
                    item,
                    "refused_value_already_set",
                    "stored %r; this is a FILL and will not overwrite a curated "
                    "value (D-039). Use correct_conditions for that." % stored,
                    stored,
                )
            )
            continue
        if item.needs_evidence and not item.quote.strip():
            out.append(
                FillOutcome(
                    item,
                    "refused_missing_quote",
                    "%r is a positive claim and needs a verbatim quote" % item.value,
                    stored,
                )
            )
            continue
        if not item.needs_evidence and item.quote.strip():
            out.append(
                FillOutcome(
                    item,
                    "refused_quote_on_an_absence",
                    "%r asserts an absence and an absence cannot be quoted "
                    "(A009, A011); drop the quote or change the value"
                    % item.value,
                    stored,
                )
            )
            continue
        out.append(FillOutcome(item, "filled", None, stored))
    return out


def apply_assignments(
    records: Sequence[TierARecord],
    assignments: Sequence[Assignment],
    *,
    known_genes: Optional[Sequence[str]] = None,
    fulltext_dir: Optional[str] = None,
    decision: str = "D-062",
) -> FillResult:
    """Fill every applicable assignment, attach evidence, and re-check quotes.

    ``records`` is not mutated.
    """
    genes = list(known_genes) if known_genes is not None else _known_genes()
    before = [TierARecord.from_dict(r.to_dict()) for r in records]
    result = FillResult()
    result.hard_findings_before = _hard(before)
    outcomes = plan_assignments(before, assignments, genes)
    by_id = {r.record_id: r for r in before if r.record_id}

    for item in outcomes:
        if item.outcome != "filled":
            continue
        a = item.assignment
        record = by_id[a.record_id]
        data = record.to_dict()
        data[DIMENSION] = a.value
        if a.needs_evidence:
            entries = list(data.get("condition_evidence") or [])
            entries = [
                e for e in entries
                if (e.get("dimension") if isinstance(e, dict) else e.dimension)
                != DIMENSION
            ]
            entries.append(
                ConditionEvidence(
                    dimension=DIMENSION,
                    value=a.value,
                    source=ConditionSource(a.source),
                    locator=a.locator,
                    quote=a.quote,
                    rule=a.rule or fill_rule(decision),
                ).to_dict()
            )
            data["condition_evidence"] = entries
            result.evidence_added += 1
        rebuilt = TierARecord.from_dict(data)
        if rebuilt.record_id != record.record_id:
            raise IdentityMovedError(
                "record_id moved from %r to %r. %s is not a component of "
                "RECORD_ID_FIELDS (%s), so no identity may change and this "
                "pass has no id map to resolve one that did."
                % (record.record_id, rebuilt.record_id, DIMENSION,
                   ", ".join(RECORD_ID_FIELDS))
            )
        by_id[a.record_id] = rebuilt

    for record in before:
        current = by_id.get(record.record_id, record)
        result.ids_unchanged += int(current.record_id == record.record_id)
        result.records.append(current)

    result.outcomes = outcomes
    result.hard_findings_after = _hard(result.records)
    if fulltext_dir is not None:
        checked, failures = check_evidence_locators(result.records, fulltext_dir)
        result.locators_checked = checked
        result.locator_failures = failures
    return result


def _hard(records: Sequence[TierARecord]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for record in records:
        for finding in validate(record):
            if finding.severity is Severity.HARD:
                out[finding.code] = out.get(finding.code, 0) + 1
    return out


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #

def _table(header: Sequence[str], rows: Sequence[Sequence[object]]) -> List[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    for row in rows:
        out.append("| " + " | ".join(str(c) for c in row) + " |")
    return out


def _dist(records: Sequence[TierARecord]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for record in records:
        value = getattr(record, DIMENSION, None)
        value = getattr(value, "value", value)
        key = "null" if value is None else (
            value if value in SENTINELS else "segment(s)"
        )
        out[key] = out.get(key, 0) + 1
    return dict(sorted(out.items()))


def render_report_md(
    result: FillResult,
    before: Sequence[TierARecord],
    *,
    generated: str,
    records_path: str,
    table_path: str,
) -> str:
    data = result.to_dict()
    lines: List[str] = []
    add = lines.append
    add("# Background genotype fill (D-060, result in D-062)")
    add("")
    add(
        "Schema v0.4.0 declared `background_genotype` and left it null on all "
        "508 records. This pass fills it from a declared table, reviewable in a "
        "diff, produced by eight readers — one per study carrying an "
        "adversarially confirmed background genotype — each deciding every "
        "record of its study under the arm rule."
    )
    add("")
    add("## Run")
    add("")
    lines.extend(
        _table(
            ["field", "value"],
            [
                ["generated (UTC)", generated],
                ["input", "`%s`" % records_path],
                ["assignments table", "`%s`" % table_path],
                ["records", data["n_records"]],
                ["**filled**", "**%s**" % data["n_filled"]],
                ["refused", data["n_refused"]],
                ["evidence entries added", data["evidence_entries_added"]],
                ["quotes re-checked at their locator", data["locators_checked"]],
                ["locator failures", len(data["locator_failures"])],
                ["**record ids unchanged**", "**%s of %s**" % (data["ids_unchanged"], data["n_records"])],
                ["identity map required", "no"],
            ],
        )
    )
    add("")
    counts = data["counts"]
    shown = ", ".join("`%s` %d" % (k, v) for k, v in counts.items() if v)
    add("Outcomes: " + (shown or "none") + ".")
    add("")
    add("## The dimension, before and after")
    add("")
    first, second = _dist(before), _dist(result.records)
    keys = sorted(set(first) | set(second))
    lines.extend(
        _table(
            ["value", "before", "after", "delta"],
            [
                [
                    "`%s`" % k,
                    first.get(k, 0),
                    second.get(k, 0),
                    "%+d" % (second.get(k, 0) - first.get(k, 0)),
                ]
                for k in keys
            ],
        )
    )
    add("")
    add(
        "`null` is not a gap to be closed. D-060 defines it as \"the curator has "
        "not addressed genotype here\", and `condition_key` pools it with "
        "`not_reported` into `unknown`. A record absent from the table keeps it."
    )
    add("")
    add("## What `not_reported` and `wild_type` mean here, and why they differ")
    add("")
    add(
        "`not_reported` asserts that the paper was read for this record and "
        "establishes no background manipulation. It carries **no evidence "
        "entry**, because an absence cannot be quoted (A009, A011), and this "
        "pass refuses a `not_reported` row that arrives with a quote rather "
        "than dropping it silently."
    )
    add("")
    add(
        "`wild_type` is the opposite: a *positive, quoted* establishment that "
        "the animals carry nothing beyond the recorded perturbation. It does "
        "**not** pool with `unknown`, so an unstated background and an "
        "established wild-type background never share a claim cell (D-044). A "
        "background strain is not sufficient for it — a strain says what the "
        "animals are, not that nothing further was done to them."
    )
    add("")
    add(
        "`mixed` is a positive value for a record whose own quote reports a "
        "joint result across two or more genotypes. It is a fact about the "
        "record rather than an absence, and the record must **not** be split "
        "into one per genotype: deriving two observations from a sentence "
        "reporting one is what D-044 refuses in general."
    )
    add("")

    if result.refused:
        add("## Refusals")
        add("")
        lines.extend(
            _table(
                ["record", "value", "outcome", "why"],
                [
                    [
                        "`%s`" % o.assignment.record_id,
                        "`%s`" % o.assignment.value,
                        o.outcome,
                        (o.detail or "")[:160],
                    ]
                    for o in result.refused[:20]
                ],
            )
        )
        add("")

    if result.locator_failures:
        add("## Quotes that do not verify at their locator")
        add("")
        add(
            "Every one of these blocks the write. A quote carrying the wrong "
            "address sends an auditor to the wrong paragraph, which is worse "
            "than no quote because it looks sound."
        )
        add("")
        for failure in result.locator_failures[:20]:
            add(
                "* `%s` %s at `%s`: %s"
                % (
                    failure.get("record_id"),
                    failure.get("dimension"),
                    failure.get("locator"),
                    failure.get("why"),
                )
            )
        add("")

    add("## Hard findings, before and after")
    add("")
    lines.extend(
        _table(
            ["when", "hard findings"],
            [
                ["before", data["hard_findings_before"] or "none"],
                ["after", data["hard_findings_after"] or "none"],
            ],
        )
    )
    add("")
    add("## What this pass does not do")
    add("")
    add(
        "* It decides nothing. Every value is a row of the declared table, "
        "carrying its own quote, locator and reasoning."
    )
    add(
        "* It does not overwrite. A record whose `background_genotype` is "
        "already set is refused, not updated — that is "
        "`correct_conditions`'s job and D-039's rule."
    )
    add(
        "* It re-identifies nothing. `background_genotype` is not a "
        "`RECORD_ID_FIELDS` component (D-060), the pass checks that per record, "
        "and it refuses rather than writing an id map it does not have."
    )
    add("")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _repo_root() -> str:
    here = os.path.abspath(__file__)
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(here))))


def _default(*parts: str) -> str:
    return os.path.join(_repo_root(), *parts)


def default_table_path() -> str:
    return _default("data", "curated", "background_genotype_assignments.json")


def _now() -> str:
    return (
        datetime.datetime.now(datetime.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python3 -m syn4rch.curate.fill_background_genotype",
        description=(
            "Fill background_genotype from a declared, reviewable table. Fills "
            "only where the value is null and refuses to overwrite (D-039)."
        ),
    )
    parser.add_argument("records", nargs="?", help="curated JSONL to rewrite")
    parser.add_argument("--table", default=None, help="assignments JSON")
    parser.add_argument("--fulltext-dir", default=None)
    parser.add_argument("--backup", default=None)
    parser.add_argument("--report", default=None)
    parser.add_argument("--json", dest="json_path", default=None)
    parser.add_argument("--decision", default="D-062")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--allow-refusals",
        action="store_true",
        help=(
            "write even if a row was refused. Off by default: a refusal means "
            "the table and the dataset disagree, and writing anyway hides it"
        ),
    )
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    records_path = args.records or _default("data", "curated", "records.jsonl")
    table_path = args.table or default_table_path()
    fulltext_dir = args.fulltext_dir or _default("data", "raw", "fulltext")

    try:
        records = list(read_jsonl(records_path))
    except Exception as exc:
        sys.stderr.write("error: could not read %s: %s\n" % (records_path, exc))
        return 2
    if not records:
        sys.stderr.write("error: %s holds no records\n" % records_path)
        return 2
    try:
        assignments, meta = read_table(table_path)
    except Exception as exc:
        sys.stderr.write("error: could not read %s: %s\n" % (table_path, exc))
        return 2
    if not assignments:
        sys.stderr.write("error: %s declares no assignments\n" % table_path)
        return 2

    try:
        result = apply_assignments(
            records,
            assignments,
            fulltext_dir=fulltext_dir,
            decision=args.decision,
        )
    except IdentityMovedError as exc:
        sys.stderr.write("error: %s\n" % exc)
        return 1

    counts = result.counts()
    if not args.quiet:
        sys.stdout.write(
            "filled %d | already %d | refused %d | evidence +%d | ids unchanged %d\n"
            % (
                counts.get("filled", 0),
                counts.get("already_filled", 0),
                len(result.refused),
                result.evidence_added,
                result.ids_unchanged,
            )
        )
    if result.locator_failures:
        sys.stderr.write(
            "error: %d quote(s) do not verify at their locator; nothing "
            "written\n" % len(result.locator_failures)
        )
        for failure in result.locator_failures[:10]:
            sys.stderr.write("  %s\n" % failure)
        return 1
    if result.refused and not args.allow_refusals:
        sys.stderr.write(
            "error: %d row(s) refused; nothing written. Pass --allow-refusals "
            "only once each refusal is understood:\n" % len(result.refused)
        )
        for item in result.refused[:10]:
            sys.stderr.write(
                "  %s  %s  %s\n"
                % (item.assignment.record_id, item.outcome, item.detail or "")
            )
        return 1
    introduced = {
        code: n
        for code, n in result.hard_findings_after.items()
        if n > result.hard_findings_before.get(code, 0)
    }
    if introduced:
        sys.stderr.write(
            "error: the fill would introduce hard finding(s) %s; nothing "
            "written\n" % introduced
        )
        return 1

    generated = _now()
    relative = os.path.relpath(os.path.abspath(records_path), _repo_root())
    table_rel = os.path.relpath(os.path.abspath(table_path), _repo_root())
    report = render_report_md(
        result,
        records,
        generated=generated,
        records_path=relative,
        table_path=table_rel,
    )
    report_path = args.report or _default(
        "reports", "background_genotype_fill.md"
    )
    if args.dry_run:
        if not args.quiet:
            sys.stdout.write("dry run: nothing written\n")
        return 0

    backup_path = args.backup or (records_path + ".pre-background-genotype.bak")
    _atomic_write_text(backup_path, _records_text(records))
    _atomic_write_text(records_path, _records_text(result.records))
    _atomic_write_text(report_path, report)
    if args.json_path:
        _atomic_write_text(
            args.json_path, json.dumps(result.to_dict(), indent=2) + "\n"
        )
    if not args.quiet:
        sys.stdout.write(
            "wrote %s\nbackup %s\nreport %s\n"
            % (records_path, backup_path, report_path)
        )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
