"""Methodological fields a record needs before it can be read as evidence.

Four of them: a control, a group size, a species, and a perturbation valence.
Missing ones are reported, not filled in -- SCHEMA.md 2 rule 4 forbids
defaulting an unreported field.

Renamed at schema v0.2.0, and not only for tidiness. This module used to be
called ``conditions`` and to say that its four fields were "needed before a
record can be compared with another record". v0.2.0 introduced the
*experimental condition dimension* -- ``(sex, diet, age_group, challenge)``,
in ``syn4rch.qc.keys.condition_key`` -- which is the thing that actually decides
whether two records are comparable. Two records agreeing on all four fields
here are still not comparable if one is a fasted male on a high-fat diet and
the other a fed female on chow.

Leaving both meanings of the word "condition" in one codebase was not an
option. This project has twice had to undo a name that meant two things at
once (``docs/DECISIONS.md`` D-019 on ``direction``, D-029 on
``cognitive_decline``), and both times the cost was a record filed backwards.
See D-030.

These checks overlap with the WARNING-level schema findings on purpose: the
schema check answers "may this record be stored?", while this one answers "is
this record readable as evidence?", and produces one per-record summary a
curator can act on.

Python 3.10 compatible, standard library only.
"""

from typing import List

from ..schema.enums import ControlType, Valence
from .findings import QCFinding, QCSeverity

__all__ = ["REQUIRED_COMPARABILITY_FIELDS", "missing_comparability_fields", "find_missing_comparability_fields"]

#: Condition name -> what counts as reported.
REQUIRED_COMPARABILITY_FIELDS = ("control_types", "n", "species", "valence")


def missing_comparability_fields(record) -> List[str]:
    """Names of the required experimental conditions absent from ``record``."""
    missing: List[str] = []

    controls = [
        c for c in record.control_types if c is not ControlType.none_reported
    ]
    if not controls:
        missing.append("control_types")

    if (
        record.n_treatment is None
        and record.n_control is None
        and record.n_animals is None
    ):
        missing.append("n")

    if record.species is None:
        missing.append("species")

    # valence=unclear is a legitimate, expected answer (SCHEMA.md 3.4) and is
    # NOT missing. Only a null valence is missing.
    if record.valence is None:
        missing.append("valence")

    return missing


def find_missing_comparability_fields(records) -> List[QCFinding]:
    out: List[QCFinding] = []
    for record in records:
        missing = missing_comparability_fields(record)
        if not missing:
            continue
        out.append(
            QCFinding(
                code="QC_MISSING_COMPARABILITY_FIELDS",
                severity=QCSeverity.WARNING,
                message=(
                    "Record %s does not report: %s. Left null, not inferred. "
                    "Note this is not the experimental condition dimension "
                    "(sex/diet/age/challenge), which decides comparability."
                    % (record.record_id, ", ".join(missing))
                ),
                record_ids=[record.record_id],
                field_name=missing[0],
                details={"missing": missing},
            )
        )
    return out


def comparability_coverage(records) -> dict:
    """``condition -> number of records reporting it`` (for the report header)."""
    counts = {name: 0 for name in REQUIRED_COMPARABILITY_FIELDS}
    total = 0
    for record in records:
        total += 1
        absent = set(missing_comparability_fields(record))
        for name in REQUIRED_COMPARABILITY_FIELDS:
            if name not in absent:
                counts[name] += 1
    counts["_total_records"] = total
    return counts


def unclear_valence_count(records) -> int:
    """How many records honestly report ``valence=unclear`` (not a defect)."""
    return sum(1 for r in records if r.valence is Valence.unclear)
