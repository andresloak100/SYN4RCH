"""Unit consistency against the ``outcomes`` controlled vocabulary.

Each ``outcomes`` vocabulary entry declares the unit(s) its readout is
measured in. This check compares a record's ``effect_size_units`` and
``assay_units`` against that declaration.

It is deliberately conservative:

* Dimensionless effect sizes (``percent_change``, ``fold_change``,
  ``cohens_d``) are *not* compared against the outcome's physical unit -- a
  percent change in blood glucose is in percent, not mg/dL. They are only
  checked for being a recognised dimensionless unit.
* Only ``absolute_difference`` and ``auc_difference`` carry the outcome's own
  unit and are compared to it.
* If the vocabulary is unavailable or the entry declares no unit, the check
  reports that fact at ``INFO`` and compares nothing. It never invents an
  expected unit.

The vocabulary is reached only through the agreed interface
``syn4rch.vocab.load(kind) -> {id: entry}`` and is imported lazily, so this
module is importable (and testable) before that package exists.

Python 3.10 compatible, standard library only.
"""

from typing import Callable, Dict, List, Optional

from ..schema.enums import EffectSizeType
from .findings import QCFinding, QCSeverity

__all__ = [
    "DIMENSIONLESS_UNITS",
    "DIMENSIONFUL_EFFECT_TYPES",
    "normalize_unit",
    "unit_variants",
    "load_outcome_units",
    "check_units",
]

#: Accepted spellings for the unit of a dimensionless effect size, by type.
DIMENSIONLESS_UNITS = {
    EffectSizeType.percent_change: frozenset({"%", "percent", "pct", "percentage"}),
    EffectSizeType.fold_change: frozenset({"fold", "x", "ratio", "fold_change"}),
    EffectSizeType.cohens_d: frozenset({"sd", "sd_units", "d", "z"}),
}

#: Effect size types expressed in the outcome's own physical unit.
DIMENSIONFUL_EFFECT_TYPES = frozenset(
    {EffectSizeType.absolute_difference, EffectSizeType.auc_difference}
)

#: Unit spellings treated as the same unit. Mapping is textual only; no
#: numeric conversion is ever performed, and no record is rewritten.
_UNIT_ALIASES = {
    "mg/dl": "mg/dl",
    "mgdl": "mg/dl",
    "mg per dl": "mg/dl",
    "mmol/l": "mmol/l",
    "mm": "mmol/l",
    "ng/ml": "ng/ml",
    "pg/ml": "pg/ml",
    "ul": "uL",
    "µl": "uL",
    "degc": "degC",
    "°c": "degC",
    "c": "degC",
    "kcal/day": "kcal/day",
    "kcal/d": "kcal/day",
    "g": "g",
    "gram": "g",
    "grams": "g",
    "pct": "%",
    "percent": "%",
    "percentage": "%",
    "%": "%",
}


def normalize_unit(unit) -> str:
    """Case-fold, strip whitespace, collapse spaces, apply textual aliases."""
    if unit is None:
        return ""
    text = " ".join(str(unit).split()).strip().casefold()
    text = text.replace(" / ", "/")
    return _UNIT_ALIASES.get(text, text)


def unit_variants(unit):
    """Accepted spellings of one declared unit.

    Vocabulary entries annotate units in prose, e.g. ``"degC (BAT
    temperature)"``, so the parenthetical qualifier is accepted as optional:
    both ``degc (bat temperature)`` and ``degc`` match. Nothing else is
    stripped, and no numeric conversion is implied.
    """
    full = normalize_unit(unit)
    out = {full}
    if "(" in full:
        out.add(normalize_unit(full.split("(", 1)[0]))
    return {v for v in out if v}


def _expected_units(entry) -> List[str]:
    """Pull declared unit(s) out of a vocabulary entry, tolerantly."""
    if not isinstance(entry, dict):
        return []
    raw = entry.get("units", entry.get("unit"))
    if raw is None:
        return []
    if isinstance(raw, (list, tuple, set)):
        return [str(u) for u in raw if u is not None and str(u).strip()]
    if str(raw).strip():
        return [str(raw)]
    return []


def load_outcome_units(loader: Optional[Callable[[str], Dict]] = None):
    """Return ``(units_by_outcome_id, note)``.

    ``units_by_outcome_id`` maps ``outcome_id -> [declared unit, ...]``.
    ``note`` is ``None`` on success, or a human-readable reason the
    vocabulary could not be consulted. ``loader`` defaults to
    ``syn4rch.vocab.load`` and exists so tests can inject a stub.
    """
    if loader is None:
        try:
            from ..vocab import load as loader  # type: ignore
        except Exception as exc:  # ImportError, or vocab files missing
            return {}, "outcomes vocabulary unavailable (%s: %s)" % (
                type(exc).__name__,
                exc,
            )
    try:
        entries = loader("outcomes")
    except Exception as exc:
        return {}, "outcomes vocabulary could not be loaded (%s: %s)" % (
            type(exc).__name__,
            exc,
        )
    if not isinstance(entries, dict):
        return {}, "outcomes vocabulary did not return a mapping"
    return (
        {str(k): _expected_units(v) for k, v in entries.items()},
        None,
    )


def check_units(records, loader: Optional[Callable[[str], Dict]] = None) -> List[QCFinding]:
    """Compare each record's units against its outcome's declared unit(s)."""
    units_by_outcome, note = load_outcome_units(loader)
    if note is not None:
        return [
            QCFinding(
                code="QC_VOCAB_UNAVAILABLE",
                severity=QCSeverity.INFO,
                message=(
                    "Unit consistency not checked: %s. This is a tooling gap, "
                    "not a statement about the records." % note
                ),
                field_name="effect_size_units",
            )
        ]

    out: List[QCFinding] = []
    for record in records:
        outcome_id = record.outcome_id
        if not outcome_id:
            continue  # HARD schema finding elsewhere
        if outcome_id not in units_by_outcome:
            out.append(
                QCFinding(
                    code="QC_OUTCOME_ID_NOT_IN_VOCAB",
                    severity=QCSeverity.WARNING,
                    message=(
                        "outcome_id %r is not in the outcomes vocabulary; its "
                        "units cannot be checked." % outcome_id
                    ),
                    record_ids=[record.record_id],
                    field_name="outcome_id",
                )
            )
            continue

        declared = units_by_outcome[outcome_id]
        expected = sorted({v for u in declared for v in unit_variants(u)})
        effect_units = normalize_unit(record.effect_size_units)
        assay_units = normalize_unit(record.assay_units)
        etype = record.effect_size_type

        if not expected:
            if effect_units or assay_units:
                out.append(
                    QCFinding(
                        code="QC_OUTCOME_UNITS_UNDECLARED",
                        severity=QCSeverity.INFO,
                        message=(
                            "outcome_id %r declares no units in the vocabulary, "
                            "so %r cannot be verified."
                            % (outcome_id, record.effect_size_units or record.assay_units)
                        ),
                        record_ids=[record.record_id],
                        field_name="outcome_id",
                    )
                )
            continue

        if assay_units and assay_units not in expected:
            out.append(
                QCFinding(
                    code="QC_UNIT_MISMATCH",
                    severity=QCSeverity.WARNING,
                    message=(
                        "assay_units=%r for outcome %r, which the vocabulary "
                        "declares as %s."
                        % (record.assay_units, outcome_id, "; ".join(declared))
                    ),
                    record_ids=[record.record_id],
                    field_name="assay_units",
                    details={"expected": list(declared), "found": record.assay_units},
                )
            )

        if not effect_units:
            continue

        if etype in DIMENSIONLESS_UNITS:
            allowed = {normalize_unit(u) for u in DIMENSIONLESS_UNITS[etype]}
            if effect_units not in allowed:
                out.append(
                    QCFinding(
                        code="QC_UNIT_DIMENSIONLESS_MISMATCH",
                        severity=QCSeverity.WARNING,
                        message=(
                            "effect_size_type=%s is dimensionless but "
                            "effect_size_units=%r is not a recognised "
                            "dimensionless unit (%s)."
                            % (
                                etype.value,
                                record.effect_size_units,
                                ", ".join(sorted(allowed)),
                            )
                        ),
                        record_ids=[record.record_id],
                        field_name="effect_size_units",
                        details={"expected": sorted(allowed), "found": record.effect_size_units},
                    )
                )
        elif etype in DIMENSIONFUL_EFFECT_TYPES:
            if effect_units not in expected:
                out.append(
                    QCFinding(
                        code="QC_UNIT_MISMATCH",
                        severity=QCSeverity.WARNING,
                        message=(
                            "effect_size_units=%r with effect_size_type=%s for "
                            "outcome %r, which the vocabulary declares as %s."
                            % (
                                record.effect_size_units,
                                etype.value,
                                outcome_id,
                                "; ".join(declared),
                            )
                        ),
                        record_ids=[record.record_id],
                        field_name="effect_size_units",
                        details={"expected": list(declared), "found": record.effect_size_units},
                    )
                )
        else:
            out.append(
                QCFinding(
                    code="QC_UNIT_UNCHECKABLE",
                    severity=QCSeverity.INFO,
                    message=(
                        "effect_size_units=%r present but effect_size_type is "
                        "null, so the unit cannot be interpreted."
                        % (record.effect_size_units,)
                    ),
                    record_ids=[record.record_id],
                    field_name="effect_size_type",
                )
            )
    return out
