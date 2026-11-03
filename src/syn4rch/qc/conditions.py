"""Missing experimental condition detection.

Four conditions are needed before a record can be compared with another
record: a control, a group size, a species, and a perturbation valence.
Missing ones are reported, not filled in -- SCHEMA.md 2 rule 4 forbids
defaulting an unreported field.

These overlap with the WARNING-level schema findings on purpose: the schema
check answers "may this record be stored?", while this check answers "can
this record be used in a comparison?", and produces one per-record summary a
curator can act on.

Python 3.10 compatible, standard library only.
"""

from typing import List

from ..schema.enums import ControlType, Valence
from .findings import QCFinding, QCSeverity

__all__ = ["REQUIRED_CONDITIONS", "missing_conditions", "find_missing_conditions"]

#: Condition name -> what counts as reported.
REQUIRED_CONDITIONS = ("control_types", "n", "species", "valence")


def missing_conditions(record) -> List[str]:
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


def find_missing_conditions(records) -> List[QCFinding]:
    out: List[QCFinding] = []
    for record in records:
        missing = missing_conditions(record)
        if not missing:
            continue
        out.append(
            QCFinding(
                code="QC_MISSING_CONDITIONS",
                severity=QCSeverity.WARNING,
                message=(
                    "Record %s does not report: %s. Left null, not inferred."
                    % (record.record_id, ", ".join(missing))
                ),
                record_ids=[record.record_id],
                field_name=missing[0],
                details={"missing": missing},
            )
        )
    return out


def condition_coverage(records) -> dict:
    """``condition -> number of records reporting it`` (for the report header)."""
    counts = {name: 0 for name in REQUIRED_CONDITIONS}
    total = 0
    for record in records:
        total += 1
        absent = set(missing_conditions(record))
        for name in REQUIRED_CONDITIONS:
            if name not in absent:
                counts[name] += 1
    counts["_total_records"] = total
    return counts


def unclear_valence_count(records) -> int:
    """How many records honestly report ``valence=unclear`` (not a defect)."""
    return sum(1 for r in records if r.valence is Valence.unclear)
