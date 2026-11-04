"""Fill the condition dimension from the papers' own text (v0.2.0 -> v0.3.0).

This is the second writer to touch ``data/curated/records.jsonl`` after the
curation pass that created it, and the second migration in the project's
history (D-031 established that a migration is a writer like any other and
must leave a map behind). What it does is narrow:

* read ``sex``, ``diet`` and ``age_group`` out of each paper's Methods, but
  **only for the dimensions that paper holds constant** across every
  experiment it reports;
* read ``challenge`` out of each record's own passage, because the acute
  metabolic state is a property of the measurement and no Methods section says
  which cohort an individual figure panel came from;
* attach a :class:`~syn4rch.schema.record.ConditionEvidence` entry -- verbatim
  sentence plus locator -- to every value it fills or confirms;
* recompute ``record_id`` and leave a complete old-to-new map, because all four
  condition fields are hash components (D-026).

What it refuses to do matters as much. It never writes a value the paper does
not state, it never overwrites a curated value with a paper-level one, and it
never resolves a disagreement between the two: a stored ``male`` against a
Methods section that describes both sexes is reported and left alone. The
dataset's job is to preserve that, not to pick.

The four outcomes per (record, dimension), which the report counts separately:

``filled``
    the field was ``null`` and the paper states a value. ``record_id`` changes.
``corrected``
    the field said ``not_reported`` and the paper states a value. The old value
    was asserted against the record's own passage by an extraction that could
    not see the Methods, so at paper scope it was simply wrong. ``record_id``
    changes.
``confirmed``
    the field already held this value. Only evidence is attached, so the
    ``record_id`` does not change. This is how a value curated before v0.3.0
    stops being an unevidenced assertion (W017).
``conflicted``
    the field holds a different positive value. Nothing is applied.

Python 3.10, standard library only (D-001).
"""

import argparse
import json
import os
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ..schema.record import (
    CONDITION_DIMENSIONS,
    RECORD_ID_FIELDS,
    ConditionEvidence,
    TierARecord,
)
from ..schema.validate import hard_findings, validate
from .conditions import PaperConditions, read_paper_conditions, read_record_conditions
from .jats import parse_jats_file
from .passages import build_index

# ``read_records`` and ``remap_conflicts`` are the migration reader and the
# conflicts_with rewriter, and ``_atomic_write_text`` is shared rather than
# reimplemented on purpose: two implementations of an atomic replace is how a
# project ends up with one that is not atomic and a half-written dataset.
from .migrate import (  # noqa: F401
    _atomic_write_text,
    _records_text,
    remap_conflicts,
)

__all__ = [
    "FROM_VERSION",
    "ACCEPTED_INPUT_VERSIONS",
    "TO_VERSION",
    "FILLED",
    "CORRECTED",
    "CONFIRMED",
    "CONFLICTED",
    "DimensionOutcome",
    "RecordPlan",
    "FillResult",
    "plan_record",
    "fill_records",
    "read_v0_2_0_records",
    "RETIRED_FIELDS",
    "fill_stats",
    "check_evidence_locators",
    "verify_plan",
    "REFUSED_NOT_AT_LOCATOR",
    "render_report_md",
    "build_parser",
    "main",
]

#: The version this pass reads.
FROM_VERSION = "0.2.0"
#: The version it writes.
TO_VERSION = "0.3.0"

#: Input versions this pass accepts. ``0.2.0`` is the original migration;
#: ``0.3.0`` is a **top-up**, which is what a later curation pass needs: new
#: records arrive with whatever conditions their own passage stated, and the
#: Methods of their papers have never been read. A top-up fills only what is
#: missing, confirms what agrees, and refuses to overwrite a curated value,
#: because those are the same three rules the first run followed.
#:
#: Running it twice over the same file is close to a no-op by construction:
#: everything it filled the first time comes back `confirmed`, which attaches
#: no second evidence entry and changes no ``record_id``.
ACCEPTED_INPUT_VERSIONS = ("0.2.0", "0.3.0")

FILLED = "filled"
CORRECTED = "corrected"
CONFIRMED = "confirmed"
CONFLICTED = "conflicted"

#: The one condition value that is an assertion about the paper rather than a
#: reading of it, and so is the one value a wider-scope reading may replace.
NOT_REPORTED = "not_reported"

#: Dimensions this pass may take from a paper's Methods. ``challenge`` is
#: deliberately absent: an acute state is a property of a measurement, and a
#: Methods section that describes a fast and a glucose tolerance test does not
#: say which one any given figure panel reports (see
#: ``conditions.refused_housing_statement_is_not_a_challenge``).
METHODS_DIMENSIONS = ("sex", "diet", "age_group")

#: Dimensions this pass may take from a record's own passage: all four. A
#: reading of the passage the record came from is never worse than a reading of
#: the paper's Methods, because it describes this measurement rather than the
#: study, so the passage route is open to every dimension and wins any it
#: reaches. In practice ``challenge`` is the one it usually reaches, since a
#: Results sentence rarely restates the animals' sex.
PASSAGE_DIMENSIONS = CONDITION_DIMENSIONS

#: Fields v0.3.0 retires (D-037). A v0.2.0 file carries them, and
#: ``TierARecord.from_dict`` refuses an unknown key by design, so this pass
#: strips exactly these two and nothing else. ``strict=False`` would have read
#: the file as well and would also have swallowed a misspelled field name,
#: which is how a column quietly stops being curated.
RETIRED_FIELDS = ("condition_source", "condition_locator")


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))


def _default(*parts: str) -> str:
    return os.path.join(_repo_root(), *parts)


def read_v0_2_0_records(
    path: str,
) -> Tuple[List[TierARecord], List[Tuple[int, str]], Dict[str, int]]:
    """Read a curated JSONL, retiring exactly the two fields v0.3.0 removes.

    Accepts an input at any of :data:`ACCEPTED_INPUT_VERSIONS`. A v0.3.0 file
    carries neither retired field, so the strip is a no-op on it.

    Returns ``(records, errors, versions_seen)``, like
    :func:`syn4rch.curate.migrate.read_records`.

    A retired field holding a **non-null** value is an error, not something to
    drop: both were null on all 362 records of the first corpus, and that is
    the only reason removing them was free. If one ever held a value, this pass
    would be deleting curated data, and it refuses instead.
    """
    records: List[TierARecord] = []
    errors: List[Tuple[int, str]] = []
    versions: Dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                data = json.loads(line)
            except ValueError as exc:
                errors.append((number, "not JSON: %s" % exc))
                continue
            if not isinstance(data, dict):
                errors.append((number, "not a JSON object"))
                continue
            kept = []
            for name in RETIRED_FIELDS:
                if name in data:
                    value = data.pop(name)
                    if value is not None:
                        kept.append("%s=%r" % (name, value))
            if kept:
                errors.append(
                    (
                        number,
                        "refusing to drop a retired field that holds a value: %s"
                        % ", ".join(kept),
                    )
                )
                continue
            try:
                record = TierARecord.from_dict(data)
            except Exception as exc:
                errors.append((number, "%s: %s" % (type(exc).__name__, exc)))
                continue
            key = str(record.schema_version)
            versions[key] = versions.get(key, 0) + 1
            records.append(record)
    return records, errors, versions


def _value_of(record: TierARecord, dimension: str) -> Optional[str]:
    raw = getattr(record, dimension, None)
    return getattr(raw, "value", raw)


@dataclass
class DimensionOutcome:
    """What happened to one condition dimension on one record."""

    dimension: str
    outcome: str
    stored: Optional[str]
    extracted: Optional[str]
    source: Optional[str] = None
    rule: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "dimension": self.dimension,
            "outcome": self.outcome,
            "stored": self.stored,
            "extracted": self.extracted,
            "source": self.source,
            "rule": self.rule,
        }


@dataclass
class RecordPlan:
    """Everything this pass decided about one record, before anything is written."""

    record_id_in: str
    study_id: Optional[str]
    updates: Dict[str, str] = field(default_factory=dict)
    evidence: List[ConditionEvidence] = field(default_factory=list)
    outcomes: List[DimensionOutcome] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def changes_identity(self) -> bool:
        return any(name in RECORD_ID_FIELDS for name in self.updates)

    def counts(self) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for item in self.outcomes:
            out[item.outcome] = out.get(item.outcome, 0) + 1
        return out


def _as_evidence(extracted, dimension: str) -> ConditionEvidence:
    return ConditionEvidence(
        dimension=dimension,
        value=extracted.value,
        source=extracted.source,
        locator=extracted.locator,
        quote=extracted.quote,
        rule=extracted.rule,
    )


def plan_record(
    record: TierARecord,
    paper: Optional[PaperConditions],
    record_extractions: Optional[Sequence[object]] = None,
) -> RecordPlan:
    """Decide what to do with one record. Nothing is modified.

    Precedence is fixed and is the whole argument of this module: a reading of
    the record's own passage outranks a reading of the paper's Methods, because
    the passage describes this measurement and the Methods describe the study.
    A Methods-level value therefore never overwrites a curated one; where the
    two disagree the disagreement is recorded and the curated value stands.

    ``record_extractions`` exists so the precedence can be tested without
    reaching into the extractor; left ``None`` it reads the record itself.
    """
    plan = RecordPlan(record_id_in=record.record_id, study_id=record.study_id)
    already = {
        e.dimension for e in (record.condition_evidence or []) if e.dimension
    }

    # Record-level first, so that a passage reading wins a dimension outright.
    found: Dict[str, object] = {}
    if record_extractions is None:
        # The extractor's refusals go into the report rather than vanishing: a
        # record whose own wording names two states, or a state the enum cannot
        # express, is a curation question and not a silence.
        refusals: List[str] = []
        record_extractions = read_record_conditions(record, refusals)
        plan.notes.extend(refusals)
    for extracted in record_extractions:
        if extracted.dimension in PASSAGE_DIMENSIONS:
            found.setdefault(extracted.dimension, extracted)
    if paper is not None:
        for dimension in METHODS_DIMENSIONS:
            if dimension in found:
                continue
            constant = paper.constants.get(dimension)
            if constant is not None:
                found[dimension] = constant

    for dimension in CONDITION_DIMENSIONS:
        extracted = found.get(dimension)
        if extracted is None:
            continue
        stored = _value_of(record, dimension)
        value = extracted.value
        if stored == value:
            outcome = CONFIRMED
        elif stored is None:
            outcome = FILLED
        elif stored == NOT_REPORTED:
            outcome = CORRECTED
        else:
            outcome = CONFLICTED

        plan.outcomes.append(
            DimensionOutcome(
                dimension=dimension,
                outcome=outcome,
                stored=stored,
                extracted=value,
                source=extracted.source,
                rule=extracted.rule,
            )
        )
        if outcome == CONFLICTED:
            plan.notes.append(
                "%s: curated %r, paper says %r (%s, %s) -- left as curated"
                % (dimension, stored, value, extracted.source, extracted.locator)
            )
            continue
        if outcome in (FILLED, CORRECTED):
            plan.updates[dimension] = value
        if dimension not in already:
            plan.evidence.append(_as_evidence(extracted, dimension))

    # age_value / age_unit ride along with age_group when the paper states a
    # single scalar age: they are not hash components, and leaving them null
    # beside a filled age_group would make the band unverifiable by arithmetic.
    if (
        paper is not None
        and "age_group" in plan.updates
        and paper.age_value is not None
        and record.age_value is None
    ):
        plan.updates["age_value"] = paper.age_value
        plan.updates["age_unit"] = paper.age_unit
    if (
        paper is not None
        and paper.species is not None
        and record.species is None
        and ("age_group" in plan.updates or "age_group" in already)
    ):
        plan.updates["species"] = paper.species

    return plan


#: Name of the refusal applied when a proposed evidence entry's quote is not
#: present in the passage its locator names.
REFUSED_NOT_AT_LOCATOR = "refused_quote_not_at_its_locator"


def verify_plan(plan: RecordPlan, index) -> int:
    """Drop evidence whose quote is absent from the passage it names.

    Returns the number dropped. **The value goes with the evidence**: where an
    entry was the only support for a fill or a correction, the update is
    withdrawn too, because a value whose evidence was just refused is a value
    with nothing behind it.

    This is what caught the route that read a state out of the record's own
    ``outcome_raw``. ``outcome_raw`` is the paper's own wording of the readout,
    but it is a *field*, not a located quotation: a curator writing
    ``"food intake (fed)"`` was condensing a figure panel's label, and that
    string does not occur in the Results paragraph the record's
    ``source_locator`` names. 8 of 23 such extractions do occur there and are
    kept; the other 15 are refused. Treating a field as a quote was an error in
    the specification of that rule, and the per-entry check is what makes the
    consequence a counted refusal rather than a plausible-looking entry
    pointing at a paragraph that does not contain it.
    """
    from .verify import quote_occurs_in

    if not plan.evidence:
        return 0
    kept: List[ConditionEvidence] = []
    dropped = 0
    for entry in plan.evidence:
        passages, _how = index.resolve(entry.locator) if index is not None else ([], "")
        ok = bool(passages) and any(
            quote_occurs_in(entry.quote, passage.text) for passage in passages
        )
        if ok:
            kept.append(entry)
            continue
        dropped += 1
        plan.notes.append(
            "%s: %s evidence quoted %r at %r, which does not contain it"
            % (
                REFUSED_NOT_AT_LOCATOR,
                entry.dimension,
                entry.quote[:60],
                entry.locator,
            )
        )
        plan.updates.pop(entry.dimension, None)
        plan.outcomes = [o for o in plan.outcomes if o.dimension != entry.dimension]
    plan.evidence = kept
    return dropped


@dataclass
class FillResult:
    records: List[TierARecord] = field(default_factory=list)
    plans: List[RecordPlan] = field(default_factory=list)
    papers: Dict[str, PaperConditions] = field(default_factory=dict)
    id_map: Dict[str, str] = field(default_factory=dict)
    counts: Dict[str, int] = field(default_factory=dict)
    collisions: List[Tuple[str, List[str]]] = field(default_factory=list)
    missing_fulltext: List[str] = field(default_factory=list)
    unmapped_conflicts: List[Tuple[str, str]] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


def _paper_path(record: TierARecord, fulltext_dir: str) -> Optional[str]:
    if not record.pmcid:
        return None
    path = os.path.join(fulltext_dir, "%s.xml" % record.pmcid)
    return path if os.path.exists(path) else None


def fill_records(
    records: Sequence[TierARecord], fulltext_dir: str
) -> FillResult:
    """Apply the pass to a whole dataset, in memory.

    Every record is rewritten to :data:`TO_VERSION` whether or not its values
    changed, because the file is one dataset at one schema version and a mixed
    file is what the version check exists to refuse.
    """
    result = FillResult()

    # One parse per paper, not one per record: 362 records come from 42 papers,
    # and the parsed document is reused for the locator index below.
    papers: Dict[str, Optional[PaperConditions]] = {}
    indexes: Dict[str, object] = {}
    for record in records:
        study = record.study_id
        if study in papers:
            continue
        path = _paper_path(record, fulltext_dir)
        if path is None:
            papers[study] = None
            indexes[record.pmcid or ""] = None
            result.missing_fulltext.append(study)
            continue
        parsed = parse_jats_file(path)
        papers[study] = read_paper_conditions(parsed, study, record.pmcid)
        indexes[record.pmcid or ""] = build_index(parsed)
    result.papers = {k: v for k, v in papers.items() if v is not None}

    out: List[TierARecord] = []
    dropped_total = 0
    for record in records:
        plan = plan_record(record, papers.get(record.study_id))
        dropped_total += verify_plan(plan, indexes.get(record.pmcid or ""))
        result.plans.append(plan)

        data = record.to_dict()
        data.update(plan.updates)
        data["schema_version"] = TO_VERSION
        existing = list(record.condition_evidence or [])
        data["condition_evidence"] = [
            e.to_dict() for e in existing + plan.evidence
        ]
        # Drop the stored id so __post_init__ recomputes it from the new
        # values. An id carried over from before the fill would be a lie about
        # what the record now says.
        data["record_id"] = None
        filled = TierARecord.from_dict(data)
        result.id_map[record.record_id] = filled.record_id
        out.append(filled)

    # Filling a field can only ever make two records more distinguishable --
    # except where two records differed *only* in that one was unfilled. Then
    # the fill merges their identities, which would silently drop a row.
    by_new: Dict[str, List[str]] = {}
    for old, new in result.id_map.items():
        by_new.setdefault(new, []).append(old)
    result.collisions = [
        (new, sorted(olds)) for new, olds in sorted(by_new.items()) if len(olds) > 1
    ]

    result.records = out
    # The map is complete over the input, so an unmapped reference is a
    # dangling pointer that was already there. Reported, never dropped.
    result.unmapped_conflicts = remap_conflicts(out, result.id_map)

    changed = sum(1 for old, new in result.id_map.items() if old != new)
    tally: Dict[str, int] = {
        "evidence_refused_not_at_locator": dropped_total,
        "records_in": len(records),
        "records_out": len(out),
        "papers": len(result.papers),
        "ids_changed": changed,
        "ids_unchanged": len(result.id_map) - changed,
    }
    for plan in result.plans:
        for name, count in plan.counts().items():
            tally[name] = tally.get(name, 0) + count
    result.counts = tally
    return result


# --------------------------------------------------------------------------- #
# statistics
# --------------------------------------------------------------------------- #


def _dist(values: Iterable[Optional[str]]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for value in values:
        key = "null" if value is None else str(value)
        out[key] = out.get(key, 0) + 1
    return dict(sorted(out.items(), key=lambda kv: (-kv[1], kv[0])))


def fill_stats(
    before: Sequence[TierARecord],
    after: Sequence[TierARecord],
    result: FillResult,
    fulltext_dir: Optional[str] = None,
) -> dict:
    """Numbers the report prints. Computed from the records, not from the plans.

    Recomputing from the written records rather than trusting the plan tally is
    the point: a plan that says it filled 80 values and a dataset that holds 70
    is a bug, and only an independent count can see it.
    """
    dims: Dict[str, dict] = {}
    for dimension in CONDITION_DIMENSIONS:
        dims[dimension] = {
            "before": _dist(_value_of(r, dimension) for r in before),
            "after": _dist(_value_of(r, dimension) for r in after),
        }

    outcomes: Dict[str, Dict[str, int]] = {}
    rules: Dict[str, int] = {}
    for plan in result.plans:
        for item in plan.outcomes:
            bucket = outcomes.setdefault(item.dimension, {})
            bucket[item.outcome] = bucket.get(item.outcome, 0) + 1
            if item.outcome != CONFLICTED and item.rule:
                rules[item.rule] = rules.get(item.rule, 0) + 1

    verdicts: Dict[str, Dict[str, int]] = {}
    for paper in result.papers.values():
        for dimension, verdict in paper.verdicts.items():
            bucket = verdicts.setdefault(dimension, {})
            bucket[verdict] = bucket.get(verdict, 0) + 1

    evidenced = 0
    evidence_values = 0
    for record in after:
        entries = record.condition_evidence or []
        if entries:
            evidenced += 1
        evidence_values += len(entries)

    refusals: Dict[str, int] = {}
    for plan in result.plans:
        for note in plan.notes:
            if note.startswith("refused_"):
                name = note.split(":", 1)[0]
                refusals[name] = refusals.get(name, 0) + 1
    for paper in result.papers.values():
        for note in paper.notes:
            if note.startswith("refused_"):
                name = note.split(":", 1)[0]
                refusals[name] = refusals.get(name, 0) + 1

    conflicts = [
        {
            "record_id": plan.record_id_in,
            "study_id": plan.study_id,
            "notes": list(plan.notes),
        }
        for plan in result.plans
        if any(o.outcome == CONFLICTED for o in plan.outcomes)
    ]

    hard_before = _findings(before)
    hard_after = _findings(after)

    locators_checked, locator_failures = (
        check_evidence_locators(after, fulltext_dir)
        if fulltext_dir is not None
        else (0, [])
    )

    return {
        "records": len(after),
        "papers": result.counts.get("papers", 0),
        "dimensions": dims,
        "outcomes": outcomes,
        "rules": dict(sorted(rules.items(), key=lambda kv: (-kv[1], kv[0]))),
        "paper_verdicts": verdicts,
        "refusals": dict(sorted(refusals.items(), key=lambda kv: (-kv[1], kv[0]))),
        "records_with_evidence": evidenced,
        "evidence_entries": evidence_values,
        "ids_changed": result.counts.get("ids_changed", 0),
        "ids_unchanged": result.counts.get("ids_unchanged", 0),
        "collisions": [{"new_id": new, "from": olds} for new, olds in result.collisions],
        "conflicts": conflicts,
        "missing_fulltext": list(result.missing_fulltext),
        "unmapped_conflicts": [
            {"record_id": rid, "reference": ref}
            for rid, ref in result.unmapped_conflicts
        ],
        "locators_checked": locators_checked,
        "locator_failures": locator_failures,
        "hard_findings_before": hard_before,
        "hard_findings_after": hard_after,
    }


def check_evidence_locators(
    records: Sequence[TierARecord], fulltext_dir: str
) -> Tuple[int, List[dict]]:
    """Confirm every condition quote occurs *at the locator it names*.

    Stronger than the gate the extractor applies to itself, which asks only
    whether a quote occurs somewhere in the paper. A quote that is real but
    carries the wrong address sends an auditor to the wrong paragraph, and an
    audit trail that misdirects is worse than none because it looks sound.

    Returns ``(checked, failures)``. A record whose paper is not on disk is not
    counted as checked: an unverifiable claim is not a verified one.
    """
    from .passages import build_index
    from .verify import quote_occurs_in

    indexes: Dict[str, object] = {}
    checked = 0
    failures: List[dict] = []
    for record in records:
        entries = record.condition_evidence or []
        if not entries:
            continue
        key = record.pmcid or ""
        if key not in indexes:
            path = _paper_path(record, fulltext_dir)
            indexes[key] = (
                build_index(parse_jats_file(path)) if path is not None else None
            )
        index = indexes[key]
        if index is None:
            continue
        for entry in entries:
            passages, how = index.resolve(entry.locator)
            if not passages:
                failures.append(
                    {
                        "record_id": record.record_id,
                        "dimension": entry.dimension,
                        "locator": entry.locator,
                        "why": "locator does not resolve in this paper",
                    }
                )
                continue
            checked += 1
            if not any(quote_occurs_in(entry.quote, p.text) for p in passages):
                failures.append(
                    {
                        "record_id": record.record_id,
                        "dimension": entry.dimension,
                        "locator": entry.locator,
                        "resolved_by": how,
                        "why": "quote does not occur in the passage it names",
                    }
                )
    return checked, failures


def _findings(records: Sequence[TierARecord]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for record in records:
        for finding in hard_findings(validate(record)):
            out[finding.code] = out.get(finding.code, 0) + 1
    return dict(sorted(out.items()))


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #


def _table(header: Tuple[str, ...], rows: Sequence[Sequence[object]]) -> List[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    for row in rows:
        out.append("| " + " | ".join(str(c) for c in row) + " |")
    return out


def _inline(dist: Dict[str, int]) -> str:
    if not dist:
        return "none"
    return ", ".join("%s %d" % (k, v) for k, v in dist.items())


def render_report_md(stats: dict, *, generated: str) -> str:
    A: List[str] = []
    add = A.append
    add("# Condition fill, schema v0.2.0 -> v0.3.0")
    add("")
    add("Generated %s by `python3 -m syn4rch.curate.fill_conditions run`." % generated)
    add("")
    add(
        "Read %d records from %d papers. The pass reads `sex`, `diet` and "
        "`age_group` from each paper's Methods, and only for the dimensions "
        "that paper holds constant; it reads `challenge` from each record's own "
        "passage. Nothing is inferred: every value below is backed by a "
        "verbatim sentence and a locator stored on the record."
        % (stats["records"], stats["papers"])
    )
    add("")

    add("## What the papers could and could not settle")
    add("")
    add(
        "A dimension a paper holds constant can be applied to every record from "
        "that paper. A dimension that varies within the paper cannot be applied "
        "at all, because the Methods do not say which cohort an individual "
        "figure panel came from. `varies` is therefore not a failure of this "
        "pass; it is the paper declining to answer the question at the scope "
        "the question was asked."
    )
    add("")
    rows = []
    for dimension in CONDITION_DIMENSIONS:
        verdicts = stats["paper_verdicts"].get(dimension, {})
        rows.append(
            [
                "`%s`" % dimension,
                verdicts.get("constant", 0),
                verdicts.get("varies", 0),
                verdicts.get("absent", 0),
            ]
        )
    A.extend(_table(("dimension", "constant", "varies", "absent"), rows))
    add("")

    add("## What changed, per dimension")
    add("")
    rows = []
    for dimension in CONDITION_DIMENSIONS:
        o = stats["outcomes"].get(dimension, {})
        rows.append(
            [
                "`%s`" % dimension,
                o.get(FILLED, 0),
                o.get(CORRECTED, 0),
                o.get(CONFIRMED, 0),
                o.get(CONFLICTED, 0),
            ]
        )
    A.extend(
        _table(
            ("dimension", "filled", "corrected", "confirmed", "conflicted"), rows
        )
    )
    add("")
    add(
        "`filled` was null. `corrected` said `not_reported`, which an extraction "
        "working from one passage could assert and the Methods then contradict. "
        "`confirmed` already held the value and gained only its evidence, so its "
        "`record_id` did not change. `conflicted` holds a different positive "
        "value and was **left exactly as curated**; those records are listed "
        "below for a human to settle."
    )
    add("")

    add("## Before and after")
    add("")
    for dimension in CONDITION_DIMENSIONS:
        d = stats["dimensions"][dimension]
        add("* `%s`" % dimension)
        add("  * before: %s" % _inline(d["before"]))
        add("  * after:  %s" % _inline(d["after"]))
    add("")

    add("## Which rule read which value")
    add("")
    if stats["rules"]:
        A.extend(
            _table(
                ("rule", "values"),
                [["`%s`" % k, v] for k, v in stats["rules"].items()],
            )
        )
        add("")
        add(
            "This table is the audit handle. A rule that turns out to be wrong "
            "is wrong for every value it produced, and this says how many that "
            "is and therefore how much of the dataset a correction has to touch."
        )
    else:
        add("No rule fired. Nothing was filled.")
    add("")

    add("## Identity")
    add("")
    add(
        "All four condition fields are `record_id` components (D-026), so a "
        "filled value re-identifies its record. %d of %d ids changed; %d did "
        "not, those being the records whose values were confirmed rather than "
        "changed."
        % (
            stats["ids_changed"],
            stats["ids_changed"] + stats["ids_unchanged"],
            stats["ids_unchanged"],
        )
    )
    add("")
    if stats["collisions"]:
        add(
            "**%d identity collisions.** Two records that differed only in that "
            "one had an unfilled condition become one record once the gap is "
            "closed. That is a real duplicate pair surfacing, not a hash "
            "failure, but it is not this pass's business to merge records, so "
            "the run refused:"
            % len(stats["collisions"])
        )
        add("")
        for item in stats["collisions"]:
            add("* `%s` <- %s" % (item["new_id"], ", ".join(item["from"])))
    else:
        add("No two records collapsed onto one id.")
    add("")

    add("## Evidence")
    add("")
    add(
        "%d of %d records now carry at least one condition evidence entry, %d "
        "entries in all. An entry means a sentence was found in the retrieved "
        "text and checked verbatim; it is not the same statement as a value of "
        "`not_reported`, which asserts that no such sentence exists, and the "
        "two are never added together."
        % (stats["records_with_evidence"], stats["records"], stats["evidence_entries"])
    )
    add("")

    if stats["refusals"]:
        add("## What the extractor refused, and why (%d)" % sum(stats["refusals"].values()))
        add("")
        add(
            "A refusal is a reading the rules declined rather than a sentence "
            "they failed to find. Each is named, so the set of records a rule "
            "affects can be re-examined without re-reading the corpus."
        )
        add("")
        A.extend(
            _table(
                ("refusal", "records"),
                [["`%s`" % k, v] for k, v in stats["refusals"].items()],
            )
        )
        add("")

    if stats["conflicts"]:
        add("## Disagreements left unresolved (%d records)" % len(stats["conflicts"]))
        add("")
        add(
            "Each of these holds a curated value that the paper's own text "
            "appears to contradict. Nothing was changed. A human has to read "
            "the paper and decide, and until then the record says what the "
            "original curation said."
        )
        add("")
        for item in stats["conflicts"]:
            add("* `%s` (%s)" % (item["record_id"], item["study_id"]))
            for note in item["notes"]:
                add("  * %s" % note)
        add("")

    if stats["missing_fulltext"]:
        add("## Papers with no retrieved full text (%d)" % len(stats["missing_fulltext"]))
        add("")
        add(
            "No condition was read for these. SCHEMA.md 2 rule 5 forbids "
            "creating or extending a record from a paper whose text was not "
            "retrieved, and that applies to this pass as much as to curation."
        )
        add("")
        for study in stats["missing_fulltext"]:
            add("* %s" % study)
        add("")

    if stats["unmapped_conflicts"]:
        add(
            "## `conflicts_with` references that are not in this dataset (%d)"
            % len(stats["unmapped_conflicts"])
        )
        add("")
        add(
            "These were left exactly as they were. The id map is complete over "
            "the input, so a reference it cannot resolve was already dangling "
            "before this pass; dropping it would turn a broken pointer into an "
            "absence of conflict, which is the more expensive failure."
        )
        add("")
        for item in stats["unmapped_conflicts"]:
            add("* `%s` -> `%s`" % (item["record_id"], item["reference"]))
        add("")

    add("## Does every quote sit where it says it does?")
    add("")
    if stats["locators_checked"] or stats["locator_failures"]:
        add(
            "%d of %d evidence entries were resolved to the passage their "
            "locator names and the quote found there. %d failed."
            % (
                stats["locators_checked"],
                stats["locators_checked"] + len(stats["locator_failures"]),
                len(stats["locator_failures"]),
            )
        )
        add("")
        add(
            "This is a stricter check than the extractor's own gate, which asks "
            "only whether a quote occurs somewhere in the paper. A real quote "
            "carrying the wrong address sends an auditor to the wrong paragraph, "
            "and an audit trail that misdirects is worse than none, because it "
            "looks sound."
        )
        if stats["locator_failures"]:
            add("")
            for item in stats["locator_failures"][:20]:
                add(
                    "* `%s` %s at `%s`: %s"
                    % (
                        item["record_id"],
                        item["dimension"],
                        item["locator"],
                        item["why"],
                    )
                )
    else:
        add("No evidence entries to check.")
    add("")

    add("## Hard findings")
    add("")
    add("* before: %s" % (_inline(stats["hard_findings_before"]) or "none"))
    add("* after:  %s" % (_inline(stats["hard_findings_after"]) or "none"))
    add("")
    add(
        "A pass that introduces a hard finding has written something the schema "
        "forbids and must not be committed."
    )
    add("")
    return "\n".join(A) + "\n"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def _id_map_document(result: FillResult, *, records_path: str, generated: str) -> dict:
    return {
        "_meta": {
            "from_schema_version": FROM_VERSION,
            "to_schema_version": TO_VERSION,
            "generated": generated,
            "generated_by": "python3 -m syn4rch.curate.fill_conditions run",
            "source": os.path.relpath(records_path, _repo_root()),
            "records": result.counts.get("records_in", 0),
            "ids_changed": result.counts.get("ids_changed", 0),
            "ids_unchanged": result.counts.get("ids_unchanged", 0),
            "record_id_fields": list(RECORD_ID_FIELDS),
            "how_to_read": (
                "map[old_id] = new_id, over v0.2.0 ids as they stood after the "
                "0.1.0 -> 0.2.0 migration. Complete over the input: every id is "
                "a key exactly once. An id maps to itself where the fill only "
                "attached evidence, which is not a hash component. Follow this "
                "map to resolve any external reference, including conflicts_with "
                "lists in QC output written before the fill; the dataset's own "
                "conflicts_with lists were rewritten in place."
            ),
        },
        "map": dict(sorted(result.id_map.items())),
    }


def _write_json(path: str, document: dict) -> None:
    _atomic_write_text(path, json.dumps(document, indent=2, sort_keys=False) + "\n")


def cmd_run(args) -> int:
    records_path = args.records
    if not os.path.exists(records_path):
        print("no such file: %s" % records_path, file=sys.stderr)
        return 2

    records, errors, versions = read_v0_2_0_records(records_path)
    if errors:
        print(
            "refusing to run: %d line(s) could not be read" % len(errors),
            file=sys.stderr,
        )
        for number, message in errors[:5]:
            print("  line %d: %s" % (number, message), file=sys.stderr)
        return 1
    if not records:
        print("refusing to run: no records", file=sys.stderr)
        return 1
    if len(set(versions)) != 1:
        print(
            "refusing to run: a dataset is one version, found %s"
            % _inline(versions),
            file=sys.stderr,
        )
        return 1
    input_version = next(iter(versions))
    if input_version not in ACCEPTED_INPUT_VERSIONS:
        print(
            "refusing to run: this pass reads %s, found %s"
            % (" or ".join(ACCEPTED_INPUT_VERSIONS), input_version),
            file=sys.stderr,
        )
        return 1
    top_up = input_version == TO_VERSION

    stale = [r.record_id for r in records if not r.id_is_consistent()]
    if stale:
        print(
            "refusing to run: %d record_id(s) do not match their own fields; "
            "fix the input before re-identifying it" % len(stale),
            file=sys.stderr,
        )
        for rid in stale[:5]:
            print("  %s" % rid, file=sys.stderr)
        return 1

    before = list(records)
    result = fill_records(records, args.fulltext_dir)

    if result.collisions and not args.allow_collisions:
        print(
            "refusing to write: %d identity collision(s); filling a condition "
            "merged records that differed only by that gap"
            % len(result.collisions),
            file=sys.stderr,
        )
        for new, olds in result.collisions[:5]:
            print("  %s <- %s" % (new, ", ".join(olds)), file=sys.stderr)
        return 1

    stats = fill_stats(before, result.records, result, args.fulltext_dir)
    new_hard = {
        code: count
        for code, count in stats["hard_findings_after"].items()
        if count > stats["hard_findings_before"].get(code, 0)
    }
    if new_hard and not args.allow_new_findings:
        print(
            "refusing to write: the fill introduced hard finding(s) %s"
            % _inline(new_hard),
            file=sys.stderr,
        )
        return 1

    if stats["locator_failures"] and not args.allow_locator_failures:
        print(
            "refusing to write: %d evidence quote(s) do not occur at the "
            "locator they name" % len(stats["locator_failures"]),
            file=sys.stderr,
        )
        for item in stats["locator_failures"][:5]:
            print(
                "  %s %s at %r: %s"
                % (item["record_id"], item["dimension"], item["locator"], item["why"]),
                file=sys.stderr,
            )
        return 1

    generated = _now()
    if args.dry_run:
        print("dry run: nothing written")
        _print_summary(stats, result)
        return 0

    stem, ext = os.path.splitext(records_path)
    if args.backup == "%s.%s%s.bak" % (stem, FROM_VERSION, ext):
        args.backup = "%s.%s%s.bak" % (stem, input_version, ext)
    if args.id_map == _default(
        "data", "curated", "record_id_map_%s_to_%s.json" % (FROM_VERSION, TO_VERSION)
    ) and top_up:
        args.id_map = _default(
            "data", "curated", "record_id_map_%s_top_up.json" % TO_VERSION
        )
    if top_up and args.report == _default("reports", "condition_fill_0.3.0.md"):
        args.report = _default("reports", "condition_fill_0.3.0_top_up.md")

    if os.path.exists(args.backup) and not args.force:
        print(
            "refusing to overwrite an existing backup at %s; this pass has "
            "probably already run" % args.backup,
            file=sys.stderr,
        )
        return 1
    with open(records_path, "r", encoding="utf-8") as handle:
        _atomic_write_text(args.backup, handle.read())

    _atomic_write_text(records_path, _records_text(result.records))
    _write_json(
        args.id_map, _id_map_document(result, records_path=records_path, generated=generated)
    )
    _atomic_write_text(args.report, render_report_md(stats, generated=generated))
    if args.json:
        _write_json(args.json, stats)

    print("wrote %s (%d records at %s)" % (records_path, len(result.records), TO_VERSION))
    print("backup %s" % args.backup)
    print("id map %s" % args.id_map)
    print("report %s" % args.report)
    _print_summary(stats, result)
    return 0


def _print_summary(stats: dict, result: FillResult) -> None:
    for dimension in CONDITION_DIMENSIONS:
        o = stats["outcomes"].get(dimension, {})
        print(
            "  %-10s filled %3d  corrected %3d  confirmed %3d  conflicted %3d"
            % (
                dimension,
                o.get(FILLED, 0),
                o.get(CORRECTED, 0),
                o.get(CONFIRMED, 0),
                o.get(CONFLICTED, 0),
            )
        )
    print(
        "  ids changed %d, unchanged %d"
        % (stats["ids_changed"], stats["ids_unchanged"])
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="syn4rch-fill-conditions",
        description=(
            "Fill the condition dimension from the papers' own text, with a "
            "verbatim quote and a locator for every value (v0.2.0 -> v0.3.0)."
        ),
    )
    sub = parser.add_subparsers(dest="command")
    run = sub.add_parser("run", help="read the papers and rewrite the dataset")
    run.add_argument("records", nargs="?", help="curated JSONL to rewrite")
    run.add_argument("--fulltext-dir", default=None, help="directory of PMC*.xml")
    run.add_argument("--backup", default=None)
    run.add_argument("--id-map", default=None)
    run.add_argument("--report", default=None)
    run.add_argument("--json", default=None, help="also write the stats as JSON")
    run.add_argument("--dry-run", action="store_true")
    run.add_argument(
        "--force", action="store_true", help="overwrite an existing backup"
    )
    run.add_argument(
        "--allow-collisions",
        action="store_true",
        help="write even if two records collapse onto one id (they are duplicates)",
    )
    run.add_argument(
        "--allow-locator-failures",
        action="store_true",
        help="write even if a quote does not occur at the locator it names",
    )
    run.add_argument(
        "--allow-new-findings",
        action="store_true",
        help="write even if the fill introduces a hard schema finding",
    )
    run.set_defaults(func=cmd_run)
    return parser


def _apply_defaults(args) -> None:
    if getattr(args, "records", None) is None:
        args.records = _default("data", "curated", "records.jsonl")
    if getattr(args, "fulltext_dir", None) is None:
        args.fulltext_dir = _default("data", "raw", "fulltext")
    if getattr(args, "backup", None) is None:
        # Matches the 0.1.0 -> 0.2.0 migration's naming so the two backups sort
        # together and read as one chain: records.0.1.0.jsonl.bak, then
        # records.0.2.0.jsonl.bak. Both are gitignored (.gitignore, D-016).
        stem, ext = os.path.splitext(args.records)
        args.backup = "%s.%s%s.bak" % (stem, FROM_VERSION, ext)
    if getattr(args, "id_map", None) is None:
        args.id_map = _default(
            "data",
            "curated",
            "record_id_map_%s_to_%s.json" % (FROM_VERSION, TO_VERSION),
        )
    if getattr(args, "report", None) is None:
        args.report = _default("reports", "condition_fill_0.3.0.md")


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "func", None) is None:
        parser.print_help()
        return 2
    _apply_defaults(args)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
