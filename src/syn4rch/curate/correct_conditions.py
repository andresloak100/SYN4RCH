"""Change a condition value that was curated wrongly, under a named decision.

``syn4rch.curate.fill_conditions`` fills a condition dimension that is empty and
**refuses to overwrite one that is not** (D-039). That refusal is correct and it
leaves a gap: D-053 fixed a convention — a food-intake measurement taken after a
fast is ``refed``, because ``challenge`` is the state *at the moment of
measurement* and an animal whose intake is being measured has food — and found
six records stored the other way, plus one whose ``challenge`` came from a
comparator rather than from the animals measured. Nothing in the pipeline could
change them.

This module is the corrector, and it is deliberately not a general editor.

What it will not do
-------------------
**It does not decide anything.** Every change is a :class:`Correction` written
out in :data:`CORRECTIONS` — the record, the dimension, the value that must be
there now, the value to put there, the fragment of the stored quote that must
match, and the decision entry that authorises it. A correction this table does
not contain cannot be applied, and the table is reviewed in a diff rather than
inferred at runtime.

**It is addressed to the evidence, not only to the record.** A correction names
``quote_contains``, a fragment that must occur in the stored evidence quote.
A ``record_id`` is a hash of twelve fields, so it survives re-curation of the
quote; if the sentence behind a value is ever replaced, the correction refuses
rather than applying itself to different evidence. Both locks have to open.

**It re-reads the paper.** After the corrections are applied, every condition
quote in the whole dataset is re-checked against the passage its locator names
(the locator-strict check from the fill, not the weaker "occurs somewhere in the
paper"). A correction that would leave the audit trail misdirected is refused
before anything is written.

**It never invents a rule.** The corrected evidence entry keeps its quote, its
locator and its source, because the same sentence is still the evidence. Only
the value changes, and ``rule`` becomes
:data:`CORRECTION_RULE` + the decision id, so the entry says a convention
decided it rather than claiming an extractor read it that way.

**A removal drops the evidence.** Where ``to_value`` is ``None`` the dimension
goes back to ``null`` and its evidence entry is deleted. Leaving the quote
attached to a null dimension would be an entry asserting a value that is no
longer stored, and the validator refuses that in both directions (A009, A011).

**It is idempotent and it re-identifies.** ``challenge`` is a
``RECORD_ID_FIELDS`` component, so a correction changes the record's identity:
the id is recomputed, a complete id map is written, ``conflicts_with``
references are rewritten through it, and a collision refuses the run. Applying
the same corrections twice reports ``already_applied`` and writes nothing.

Python 3.10 compatible, standard library only.
"""

import argparse
import datetime
import json
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from ..schema.io import read_jsonl
from ..schema.record import (
    CONDITION_DIMENSIONS,
    RECORD_ID_FIELDS,
    SCHEMA_VERSION,
    TierARecord,
)
from ..schema.validate import Severity, validate
from .fill_conditions import check_evidence_locators
from .migrate import _atomic_write_text, _records_text, remap_conflicts

__all__ = [
    "CORRECTIONS",
    "CORRECTION_RULE",
    "Correction",
    "CorrectionOutcome",
    "CorrectionResult",
    "OUTCOMES",
    "apply_corrections",
    "build_parser",
    "main",
    "plan_corrections",
    "render_report_md",
]

#: Prefix of the ``rule`` written onto a corrected evidence entry. The decision
#: id is appended, so the entry names the convention that moved it and a reader
#: can find the argument rather than only the fact.
CORRECTION_RULE = "corrected_by_decision"

#: Every outcome a correction can have, as data, so the report cannot grow a
#: state the tests do not know about.
OUTCOMES: Tuple[str, ...] = (
    "applied",
    "already_applied",
    "refused_record_not_found",
    "refused_dimension_not_a_condition",
    "refused_value_is_not_the_declared_from",
    "refused_no_evidence_entry",
    "refused_quote_does_not_match",
    "refused_duplicate_correction",
)


@dataclass(frozen=True)
class Correction:
    """One declared change to one condition dimension of one record.

    ``from_value`` is not optional and not a wildcard. A corrector that would
    overwrite whatever it found is a corrector that cannot be reviewed, and the
    whole reason this module exists separately from the fill is that
    overwriting a curated value is a thing to do deliberately and narrowly.
    """

    record_id: str
    dimension: str
    from_value: Optional[str]
    to_value: Optional[str]
    quote_contains: str
    reason: str
    decision: str

    @property
    def removes(self) -> bool:
        return self.to_value is None

    def to_dict(self) -> Dict[str, object]:
        return {
            "record_id": self.record_id,
            "dimension": self.dimension,
            "from": self.from_value,
            "to": self.to_value,
            "quote_contains": self.quote_contains,
            "reason": self.reason,
            "decision": self.decision,
        }


#: The corrections D-053 identified, one per record.
#:
#: Six are the same reading mistake: ``challenge`` is the state at the moment of
#: measurement (SCHEMA.md 3.3b), and a food-intake measurement after a fast is
#: therefore ``refed`` however the paper words it — the animal is eating. All six
#: are ``food_intake``, which is the convention biting exactly where it should:
#: for any other readout "fasted" is still the state the animal is in.
#:
#: The seventh is a different error. ``c1b98d981dba36bb`` took its ``challenge``
#: from "comparable to fasted mice", where "fasted mice" is the *comparator* the
#: effect is measured against and not the state of the animals measured.
CORRECTIONS: Tuple[Correction, ...] = (
    Correction(
        record_id="0397b7caf40ad06b",
        dimension="challenge",
        from_value="fasted",
        to_value="refed",
        quote_contains="in fasted viral POMC-ChR2 mice resulted in a significant reduction in food intake",
        reason=(
            "food intake is the readout, so the animals had food at the moment "
            "of measurement; 'in fasted mice' describes the state they were in "
            "before it was returned"
        ),
        decision="D-053",
    ),
    Correction(
        record_id="087960c4021bf6f2",
        dimension="challenge",
        from_value="fasted",
        to_value="refed",
        quote_contains="food intake following an overnight fast",
        reason="deprive, return food, measure intake: the measurement is of the refeeding",
        decision="D-053",
    ),
    Correction(
        record_id="41fca07890c238c8",
        dimension="challenge",
        from_value="fasted",
        to_value="refed",
        quote_contains="hyperphagia",
        reason=(
            "fasting-induced hyperphagia is eating, which cannot be measured in "
            "an animal that still has no food"
        ),
        decision="D-053",
    ),
    Correction(
        record_id="fcea91b0669ad33f",
        dimension="challenge",
        from_value="fasted",
        to_value="refed",
        quote_contains="food intake after a fast",
        reason="deprive, return food, measure intake: the measurement is of the refeeding",
        decision="D-053",
    ),
    Correction(
        record_id="93c2e8f0c13a30ed",
        dimension="challenge",
        from_value="fasted",
        to_value="refed",
        quote_contains="suppressed feeding after a fast",
        reason="feeding after a fast is the refed state, not the fasted one",
        decision="D-053",
    ),
    Correction(
        record_id="fdb78695221a531e",
        dimension="challenge",
        from_value="fasted",
        to_value="refed",
        quote_contains="in fasted SF1-hM3Dq mice suppressed feeding",
        reason=(
            "'in fasted mice' names the state before food was returned; the "
            "readout is the feeding that followed"
        ),
        decision="D-053",
    ),
    Correction(
        record_id="c1b98d981dba36bb",
        dimension="challenge",
        from_value="fasted",
        to_value=None,
        quote_contains="comparable to fasted mice",
        reason=(
            "'fasted mice' is the comparator the effect is measured against, "
            "not the state of the animals measured; the quote supports no "
            "challenge for this record"
        ),
        decision="D-053",
    ),
)


@dataclass
class CorrectionOutcome:
    """What happened to one correction, and why."""

    correction: Correction
    outcome: str
    detail: Optional[str] = None
    old_record_id: Optional[str] = None
    new_record_id: Optional[str] = None
    stored_value: Optional[str] = None

    @property
    def refused(self) -> bool:
        return self.outcome.startswith("refused_")

    def to_dict(self) -> Dict[str, object]:
        return {
            "correction": self.correction.to_dict(),
            "outcome": self.outcome,
            "detail": self.detail,
            "old_record_id": self.old_record_id,
            "new_record_id": self.new_record_id,
            "stored_value": self.stored_value,
        }


@dataclass
class CorrectionResult:
    """Every outcome, the rewritten records, and the identity map."""

    records: List[TierARecord] = field(default_factory=list)
    outcomes: List[CorrectionOutcome] = field(default_factory=list)
    id_map: Dict[str, str] = field(default_factory=dict)
    unmapped_conflicts: List[Tuple[str, str]] = field(default_factory=list)
    locators_checked: int = 0
    locator_failures: List[dict] = field(default_factory=list)

    @property
    def applied(self) -> List[CorrectionOutcome]:
        return [o for o in self.outcomes if o.outcome == "applied"]

    @property
    def refused(self) -> List[CorrectionOutcome]:
        return [o for o in self.outcomes if o.refused]

    @property
    def ids_changed(self) -> int:
        return sum(1 for old, new in self.id_map.items() if old != new)

    def counts(self) -> Dict[str, int]:
        out = {name: 0 for name in OUTCOMES}
        for item in self.outcomes:
            out[item.outcome] = out.get(item.outcome, 0) + 1
        return out


def _value_of(record: TierARecord, dimension: str) -> Optional[str]:
    raw = getattr(record, dimension, None)
    return getattr(raw, "value", raw)


def _evidence_for(record: TierARecord, dimension: str):
    for entry in record.condition_evidence or []:
        if entry.dimension == dimension:
            return entry
    return None


def plan_corrections(
    records: Sequence[TierARecord],
    corrections: Sequence[Correction] = CORRECTIONS,
) -> List[CorrectionOutcome]:
    """Decide the outcome of every correction without changing anything.

    Separated from :func:`apply_corrections` so that ``--dry-run`` runs the
    same decisions the real pass runs, rather than a second implementation of
    them that could disagree.
    """
    index: Dict[str, TierARecord] = {}
    for record in records:
        if record.record_id:
            index[record.record_id] = record

    seen: Dict[Tuple[str, str], int] = {}
    outcomes: List[CorrectionOutcome] = []
    for correction in corrections:
        key = (correction.record_id, correction.dimension)
        seen[key] = seen.get(key, 0) + 1
        if seen[key] > 1:
            outcomes.append(
                CorrectionOutcome(
                    correction=correction,
                    outcome="refused_duplicate_correction",
                    detail=(
                        "this table already corrects %s on %s; two corrections "
                        "of one dimension cannot both be the intended value"
                        % (correction.dimension, correction.record_id)
                    ),
                )
            )
            continue
        if correction.dimension not in CONDITION_DIMENSIONS:
            outcomes.append(
                CorrectionOutcome(
                    correction=correction,
                    outcome="refused_dimension_not_a_condition",
                    detail="declared condition dimensions are %s"
                    % ", ".join(CONDITION_DIMENSIONS),
                )
            )
            continue
        record = index.get(correction.record_id)
        if record is None:
            outcomes.append(
                CorrectionOutcome(
                    correction=correction,
                    outcome="refused_record_not_found",
                    detail=(
                        "no record with this id; a correction is addressed to "
                        "an identity, and this one may already have been "
                        "applied and re-identified"
                    ),
                )
            )
            continue

        stored = _value_of(record, correction.dimension)
        if stored == correction.to_value:
            outcomes.append(
                CorrectionOutcome(
                    correction=correction,
                    outcome="already_applied",
                    detail="the dimension already holds the corrected value",
                    old_record_id=record.record_id,
                    new_record_id=record.record_id,
                    stored_value=stored,
                )
            )
            continue
        if stored != correction.from_value:
            outcomes.append(
                CorrectionOutcome(
                    correction=correction,
                    outcome="refused_value_is_not_the_declared_from",
                    detail=(
                        "stored %r, correction declares from=%r. Refusing: a "
                        "value this table did not expect was curated by "
                        "someone who saw something this correction does not "
                        "account for"
                        % (stored, correction.from_value)
                    ),
                    old_record_id=record.record_id,
                    stored_value=stored,
                )
            )
            continue

        entry = _evidence_for(record, correction.dimension)
        if entry is None:
            outcomes.append(
                CorrectionOutcome(
                    correction=correction,
                    outcome="refused_no_evidence_entry",
                    detail=(
                        "the dimension holds a value with no evidence entry, so "
                        "there is no quote to check this correction against"
                    ),
                    old_record_id=record.record_id,
                    stored_value=stored,
                )
            )
            continue
        if correction.quote_contains not in (entry.quote or ""):
            outcomes.append(
                CorrectionOutcome(
                    correction=correction,
                    outcome="refused_quote_does_not_match",
                    detail=(
                        "the stored quote does not contain %r, so the sentence "
                        "behind this value is not the one the correction was "
                        "written about"
                        % correction.quote_contains
                    ),
                    old_record_id=record.record_id,
                    stored_value=stored,
                )
            )
            continue

        outcomes.append(
            CorrectionOutcome(
                correction=correction,
                outcome="applied",
                old_record_id=record.record_id,
                stored_value=stored,
            )
        )
    return outcomes


class CorrectionCollision(ValueError):
    """Two records would share a ``record_id`` after correction."""


def apply_corrections(
    records: Sequence[TierARecord],
    corrections: Sequence[Correction] = CORRECTIONS,
    *,
    fulltext_dir: Optional[str] = None,
) -> CorrectionResult:
    """Apply every applicable correction, re-identify, and re-check the quotes.

    ``records`` is not mutated: the result carries fresh records built from
    ``to_dict``/``from_dict``, so a caller that wants the before and after can
    hold both.
    """
    before = [TierARecord.from_dict(r.to_dict()) for r in records]
    outcomes = plan_corrections(before, corrections)
    by_id = {r.record_id: r for r in before if r.record_id}

    for item in outcomes:
        if item.outcome != "applied":
            continue
        record = by_id[item.correction.record_id]
        dimension = item.correction.dimension
        to_value = item.correction.to_value
        if to_value is None:
            setattr(record, dimension, None)
            record.condition_evidence = [
                e for e in (record.condition_evidence or []) if e.dimension != dimension
            ]
        else:
            setattr(record, dimension, to_value)
            for entry in record.condition_evidence or []:
                if entry.dimension != dimension:
                    continue
                # The quote, the locator and the source are untouched: the same
                # sentence is still the evidence, and only how it is read has
                # changed. The rule says which decision read it.
                entry.value = to_value
                entry.rule = "%s_%s" % (
                    CORRECTION_RULE,
                    item.correction.decision.lower().replace("-", ""),
                )
        # Re-coerce so the enum fields hold enums rather than the strings just
        # assigned, by the same path that reading a file takes.
        rebuilt = TierARecord.from_dict(record.to_dict())
        rebuilt.record_id = rebuilt.compute_id()
        item.new_record_id = rebuilt.record_id
        by_id[item.correction.record_id] = rebuilt

    result = CorrectionResult()
    id_map: Dict[str, str] = {}
    for record in before:
        old = record.record_id or ""
        current = by_id.get(old, record)
        id_map[old] = current.record_id or old
        result.records.append(current)

    taken: Dict[str, str] = {}
    for record in result.records:
        rid = record.record_id or ""
        if rid in taken:
            raise CorrectionCollision(
                "record_id %s would be shared by two records after correction "
                "(%s and this one). RECORD_ID_FIELDS are %s; two records "
                "differing only in a corrected condition collapse into one, "
                "and silently keeping one of them would delete an observation."
                % (rid, taken[rid], ", ".join(RECORD_ID_FIELDS))
            )
        taken[rid] = rid

    result.outcomes = outcomes
    result.id_map = id_map
    result.unmapped_conflicts = remap_conflicts(result.records, id_map)
    if fulltext_dir is not None:
        checked, failures = check_evidence_locators(result.records, fulltext_dir)
        result.locators_checked = checked
        result.locator_failures = failures
    return result


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #

def _table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> List[str]:
    out = ["| " + " | ".join(header) + " |"]
    out.append("|" + "|".join("---" for _ in header) + "|")
    for row in rows:
        out.append("| " + " | ".join(str(c) for c in row) + " |")
    return out


def _dist(records: Sequence[TierARecord], dimension: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for record in records:
        value = _value_of(record, dimension)
        key = "null" if value is None else str(value)
        out[key] = out.get(key, 0) + 1
    return dict(sorted(out.items()))


def render_report_md(
    result: CorrectionResult,
    before: Sequence[TierARecord],
    *,
    generated: str,
    records_path: str,
) -> str:
    lines: List[str] = []
    add = lines.append
    add("# Condition corrections (D-053)")
    add("")
    add(
        "`fill_conditions` fills an empty condition dimension and refuses to "
        "overwrite one that is not empty (D-039). That refusal is right, and it "
        "left seven records this pass corrects: six storing `fasted` for a "
        "fast-refeed intake test, and one whose `challenge` came from a "
        "comparator rather than from the animals measured."
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
                ["schema version", SCHEMA_VERSION],
                ["records", str(len(before))],
                ["corrections declared", str(len(CORRECTIONS))],
                ["applied", str(len(result.applied))],
                ["refused", str(len(result.refused))],
                ["record ids changed", str(result.ids_changed)],
                ["condition quotes re-checked at their locator", str(result.locators_checked)],
                ["locator failures", str(len(result.locator_failures))],
            ],
        )
    )
    add("")
    counts = result.counts()
    add("Outcomes: " + ", ".join(
        "`%s` %d" % (name, counts[name]) for name in OUTCOMES if counts.get(name)
    ) + ".")
    add("")

    add("## What changed")
    add("")
    rows = []
    for item in result.outcomes:
        correction = item.correction
        rows.append(
            [
                "`%s`" % (item.old_record_id or correction.record_id),
                "`%s`" % (item.new_record_id or "-"),
                correction.dimension,
                "`%s`" % correction.from_value,
                "`%s`" % ("null" if correction.removes else correction.to_value),
                item.outcome,
            ]
        )
    lines.extend(_table(["old id", "new id", "dimension", "from", "to", "outcome"], rows))
    add("")
    add("### Why, one by one")
    add("")
    for item in result.outcomes:
        add(
            "* `%s` — %s"
            % (item.old_record_id or item.correction.record_id, item.correction.reason)
        )
    add("")

    add("## The `challenge` distribution, before and after")
    add("")
    first = _dist(before, "challenge")
    second = _dist(result.records, "challenge")
    keys = sorted(set(first) | set(second))
    lines.extend(
        _table(
            ["value", "before", "after", "delta"],
            [
                [
                    "`%s`" % key,
                    first.get(key, 0),
                    second.get(key, 0),
                    "%+d" % (second.get(key, 0) - first.get(key, 0)),
                ]
                for key in keys
            ],
        )
    )
    add("")

    if result.locator_failures:
        add("## Condition quotes that no longer verify at their locator")
        add("")
        add(
            "Every one of these blocks the run. A corrected value whose quote "
            "cannot be found at the address it names leaves an audit trail that "
            "misdirects, which is worse than none because it looks sound."
        )
        add("")
        for failure in result.locator_failures:
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

    if result.unmapped_conflicts:
        add("## `conflicts_with` references not in this dataset")
        add("")
        add(
            "Left exactly as they were. The id map is complete over the input, "
            "so a reference it cannot resolve was already dangling; dropping it "
            "would turn a broken pointer into an absence of conflict."
        )
        add("")
        for record_id, reference in result.unmapped_conflicts:
            add("* `%s` -> `%s`" % (record_id, reference))
        add("")

    add("## What this pass does not do")
    add("")
    add(
        "* It applies no judgement of its own. Every change is a row of "
        "`CORRECTIONS` in `syn4rch.curate.correct_conditions`, naming the "
        "record, the value that had to be there, the value to put there, the "
        "quote fragment that had to match and the decision that authorises it."
    )
    add(
        "* It does not re-read the papers for new conditions. That is "
        "`fill_conditions`, which only fills what is empty."
    )
    add(
        "* It does not touch any other field. `challenge` is a "
        "`RECORD_ID_FIELDS` component, so the identity moves with it and the "
        "id map is how every external reference follows."
    )
    add("")
    return "\n".join(lines) + "\n"


def _id_map_document(
    result: CorrectionResult, *, records_path: str, generated: str
) -> dict:
    return {
        "_meta": {
            "schema_version": SCHEMA_VERSION,
            "generated": generated,
            "generated_by": "python3 -m syn4rch.curate.correct_conditions run",
            "source": records_path,
            "decisions": sorted({c.decision for c in CORRECTIONS}),
            "records": len(result.records),
            "ids_changed": result.ids_changed,
            "record_id_fields": list(RECORD_ID_FIELDS),
            "how_to_read": (
                "map[old_id] = new_id. Complete over the input: every id is a "
                "key exactly once, and maps to itself where the record was not "
                "corrected. Follow it to resolve any reference written before "
                "this pass; the dataset's own conflicts_with lists were "
                "rewritten in place."
            ),
        },
        "map": dict(sorted(result.id_map.items())),
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _repo_root() -> str:
    return os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    )


def _default(*parts: str) -> str:
    return os.path.join(_repo_root(), *parts)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python3 -m syn4rch.curate.correct_conditions",
        description=(
            "Change a condition value that was curated wrongly, under a named "
            "decision. Every change is declared in CORRECTIONS; this tool "
            "decides nothing."
        ),
    )
    sub = parser.add_subparsers(dest="command", required=True)

    show = sub.add_parser("show", help="print the declared corrections")
    show.add_argument("--json", action="store_true", help="as JSON")

    run = sub.add_parser("run", help="apply them and rewrite the dataset")
    run.add_argument("records", nargs="?", help="curated JSONL to rewrite")
    run.add_argument("--fulltext-dir", default=None, help="directory of PMC*.xml")
    run.add_argument("--backup", default=None)
    run.add_argument("--id-map", default=None)
    run.add_argument("--report", default=None)
    run.add_argument("--json", dest="json_path", default=None)
    run.add_argument("--dry-run", action="store_true")
    run.add_argument(
        "--allow-refusals",
        action="store_true",
        help=(
            "write even if some correction was refused. Off by default: a "
            "refused correction means the dataset is not what the table "
            "expected, and writing anyway hides that"
        ),
    )
    return parser


def _now() -> str:
    return (
        datetime.datetime.now(datetime.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _cmd_show(args) -> int:
    if args.json:
        sys.stdout.write(
            json.dumps([c.to_dict() for c in CORRECTIONS], indent=2) + "\n"
        )
        return 0
    for correction in CORRECTIONS:
        sys.stdout.write(
            "%s  %s: %s -> %s  [%s]\n    %s\n"
            % (
                correction.record_id,
                correction.dimension,
                correction.from_value,
                "null" if correction.removes else correction.to_value,
                correction.decision,
                correction.reason,
            )
        )
    return 0


def _cmd_run(args) -> int:
    records_path = args.records or _default("data", "curated", "records.jsonl")
    fulltext_dir = args.fulltext_dir or _default("data", "raw", "fulltext")
    try:
        before = list(read_jsonl(records_path))
    except Exception as exc:  # pragma: no cover - surfaced to the operator
        sys.stderr.write("error: could not read %s: %s\n" % (records_path, exc))
        return 2
    if not before:
        sys.stderr.write("error: %s holds no records\n" % records_path)
        return 2

    try:
        result = apply_corrections(before, fulltext_dir=fulltext_dir)
    except CorrectionCollision as exc:
        sys.stderr.write("error: %s\n" % exc)
        return 1

    counts = result.counts()
    for item in result.outcomes:
        sys.stdout.write(
            "%-16s %-10s %s -> %s  %s\n"
            % (
                item.old_record_id or item.correction.record_id,
                item.outcome,
                item.correction.from_value,
                "null" if item.correction.removes else item.correction.to_value,
                item.new_record_id or "",
            )
        )
    sys.stdout.write(
        "applied %d, already applied %d, refused %d; ids changed %d\n"
        % (
            counts.get("applied", 0),
            counts.get("already_applied", 0),
            len(result.refused),
            result.ids_changed,
        )
    )
    if result.locator_failures:
        sys.stderr.write(
            "error: %d condition quote(s) no longer verify at their locator; "
            "nothing written\n" % len(result.locator_failures)
        )
        for failure in result.locator_failures:
            sys.stderr.write("  %s\n" % failure)
        return 1
    if result.refused and not args.allow_refusals:
        sys.stderr.write(
            "error: %d correction(s) refused; nothing written. Pass "
            "--allow-refusals only once each refusal is understood\n"
            % len(result.refused)
        )
        return 1

    hard = [
        f
        for record in result.records
        for f in validate(record)
        if f.severity is Severity.HARD
    ]
    if hard:
        sys.stderr.write(
            "error: %d hard schema violation(s) after correction; nothing "
            "written\n" % len(hard)
        )
        for finding in hard[:10]:
            sys.stderr.write("  %s\n" % finding)
        return 1

    generated = _now()
    relative = os.path.relpath(os.path.abspath(records_path), _repo_root())
    report = render_report_md(
        result, before, generated=generated, records_path=relative
    )
    report_path = args.report or _default(
        "reports", "condition_corrections_d053.md"
    )
    id_map_path = args.id_map or _default(
        "data", "curated", "record_id_map_d053.json"
    )

    if args.dry_run:
        sys.stdout.write("dry run: nothing written\n")
        return 0

    backup_path = args.backup or (records_path + ".pre-d053.bak")
    _atomic_write_text(backup_path, _records_text(before))
    _atomic_write_text(records_path, _records_text(result.records))
    _atomic_write_text(
        id_map_path,
        json.dumps(
            _id_map_document(result, records_path=relative, generated=generated),
            indent=2,
        )
        + "\n",
    )
    _atomic_write_text(report_path, report)
    if args.json_path:
        _atomic_write_text(
            args.json_path,
            json.dumps(
                {
                    "generated": generated,
                    "counts": counts,
                    "ids_changed": result.ids_changed,
                    "locators_checked": result.locators_checked,
                    "outcomes": [o.to_dict() for o in result.outcomes],
                },
                indent=2,
            )
            + "\n",
        )
    sys.stdout.write(
        "wrote %s\nbackup %s\nid map %s\nreport %s\n"
        % (records_path, backup_path, id_map_path, report_path)
    )
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "show":
        return _cmd_show(args)
    return _cmd_run(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
