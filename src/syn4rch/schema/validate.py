"""Record validation: the gate into ``data/curated/``.

Two severities, never blurred:

``HARD``
    The record violates SCHEMA.md 2 (anti-fabrication) or is structurally
    unusable. It **must not** enter ``data/curated/``. The QC CLI exits
    non-zero when any HARD finding exists.

``WARNING``
    The record is admissible but weaker evidence, or carries a derived field
    that disagrees with its rule. Warnings are reported, counted and
    down-weighted downstream; they never block admission and never cause a
    record to be edited automatically.

Codes are stable identifiers. Families: ``A`` anti-fabrication (SCHEMA.md 2),
``S`` structural, ``V`` value sanity, ``W`` warnings.

Python 3.10 compatible, standard library only.
"""

import re
from dataclasses import dataclass
from enum import Enum
from typing import Any, List, Optional

from .enums import (
    ASSERTIVE_DIRECTIONS,
    ControlType,
    Direction,
    EvidenceClass,
    RetrievedFrom,
)
from .record import SCHEMA_VERSION, TierARecord
from .rules import assign_age_group, assign_confidence

__all__ = [
    "Severity",
    "ValidationError",
    "RULES",
    "validate",
    "hard_findings",
    "warning_findings",
    "is_admissible",
]

_DIGIT = re.compile(r"\d")


class Severity(str, Enum):
    HARD = "HARD"
    WARNING = "WARNING"

    def __str__(self) -> str:  # pragma: no cover - trivial
        return str(self.value)


@dataclass(frozen=True)
class ValidationError:
    """One finding about one record."""

    code: str
    field: Optional[str]
    message: str
    severity: Severity = Severity.HARD
    record_id: Optional[str] = None

    @property
    def is_hard(self) -> bool:
        return self.severity is Severity.HARD

    def __str__(self) -> str:
        return "[%s] %s (%s): %s" % (
            self.severity.value,
            self.code,
            self.field if self.field else "-",
            self.message,
        )

    def to_dict(self) -> dict:
        return {
            "code": self.code,
            "field": self.field,
            "message": self.message,
            "severity": self.severity.value,
            "record_id": self.record_id,
        }


#: Human-readable catalogue of every code this module can emit.
#: ``code -> (severity, field, what it means)``
RULES = {
    # --- A: SCHEMA.md 2, anti-fabrication (hard) ------------------------
    "A001_MISSING_SUPPORT_QUOTE": (
        Severity.HARD,
        "support_quote",
        "support_quote is mandatory verbatim source text and must be non-empty.",
    ),
    "A002_MISSING_SOURCE_LOCATOR": (
        Severity.HARD,
        "source_locator",
        "source_locator is mandatory and must name where the value came from.",
    ),
    "A003_DIRECTION_WITHOUT_QUOTE": (
        Severity.HARD,
        "direction",
        "A direction other than not_reported must be derivable from a quote.",
    ),
    "A004_EFFECT_SIZE_UNSUPPORTED": (
        Severity.HARD,
        "effect_size_value",
        "A numeric effect size requires a support_quote containing a digit.",
    ),
    "A005_MISSING_OUTCOME_ID": (
        Severity.HARD,
        "outcome_id",
        "outcome_id is required; an unmappable outcome cannot be a Tier A record.",
    ),
    "A006_EVIDENCE_CLASS_NOT_CAUSAL": (
        Severity.HARD,
        "evidence_class",
        "Tier A is causal_interventional by construction; other classes belong "
        "to the knowledge graph, not this table.",
    ),
    "A007_MISSING_EXTRACTION_METHOD": (
        Severity.HARD,
        "extraction_method",
        "extraction_method must honestly record whether a human reviewed this.",
    ),
    # --- S: structural (hard) ------------------------------------------
    "S001_MISSING_STUDY_ID": (
        Severity.HARD,
        "study_id",
        "study_id is required; it anchors provenance and leakage control.",
    ),
    "S002_MISSING_POPULATION_RAW": (
        Severity.HARD,
        "population_raw",
        "population_raw is required (the paper's own wording).",
    ),
    "S003_MISSING_OUTCOME_RAW": (
        Severity.HARD,
        "outcome_raw",
        "outcome_raw is required (the paper's own wording).",
    ),
    "S004_MISSING_DIRECTION": (
        Severity.HARD,
        "direction",
        "direction is required; use not_reported rather than leaving it null.",
    ),
    "S005_RECORD_ID_INCONSISTENT": (
        Severity.HARD,
        "record_id",
        "Stored record_id does not match the hash of its identifying fields.",
    ),
    "S006_SCHEMA_VERSION_MISMATCH": (
        Severity.HARD,
        "schema_version",
        "schema_version is not the version this code implements.",
    ),
    # --- V: value sanity (hard) -----------------------------------------
    "V001_P_VALUE_OUT_OF_RANGE": (
        Severity.HARD,
        "p_value",
        "p_value must lie in [0, 1].",
    ),
    "V002_CI_BOUNDS_INVERTED": (
        Severity.HARD,
        "ci_low",
        "ci_low must not exceed ci_high.",
    ),
    "V003_NEGATIVE_COUNT": (
        Severity.HARD,
        None,
        "Animal counts must be non-negative.",
    ),
    "V004_NEGATIVE_DISPERSION": (
        Severity.HARD,
        "dispersion_value",
        "dispersion_value must be non-negative.",
    ),
    # --- W: warnings -----------------------------------------------------
    "W001_NO_CONTROL_TYPE": (
        Severity.WARNING,
        "control_types",
        "No control type recorded (or only none_reported).",
    ),
    "W002_NO_SAMPLE_SIZE": (
        Severity.WARNING,
        "n_treatment",
        "No group size reported (n_treatment, n_control and n_animals all null).",
    ),
    "W003_POPULATION_ID_NULL": (
        Severity.WARNING,
        "population_id",
        "Population did not resolve to the vocabulary; queue for extension.",
    ),
    "W004_ABSTRACT_ONLY": (
        Severity.WARNING,
        "full_text_available",
        "Abstract-only source: record cannot be re-checked against full text.",
    ),
    "W005_UNREVIEWED_LLM_EXTRACTION": (
        Severity.WARNING,
        "review_status",
        "LLM extraction with no human review yet.",
    ),
    "W006_NO_EFFECT_SIZE": (
        Severity.WARNING,
        "effect_size_value",
        "Direction reported without a quantitative effect size.",
    ),
    "W007_CONFIDENCE_RULE_MISMATCH": (
        Severity.WARNING,
        "confidence",
        "Stored confidence differs from the documented rule's output.",
    ),
    "W008_SIGNIFICANCE_WITHOUT_TEST": (
        Severity.WARNING,
        "significant",
        "significant is set but no statistical_test or p_value is recorded.",
    ),
    "W009_AGE_GROUP_RULE_MISMATCH": (
        Severity.WARNING,
        "age_group",
        "Stored age_group differs from the documented rule's output.",
    ),
    "W010_NO_LAB_IDENTITY": (
        Severity.WARNING,
        "last_author",
        "Missing last_author/affiliation: lab-level leakage control is weakened.",
    ),
    "W011_NO_OUTCOME_TIER": (
        Severity.WARNING,
        "outcome_tier",
        "outcome_tier not set (A physiological vs B aging annotation).",
    ),
    "W012_MISSING_BIBLIOGRAPHY": (
        Severity.WARNING,
        "title",
        "title and/or year missing; record is hard to audit by hand.",
    ),
}


def _blank(value: Any) -> bool:
    """True if a text field is absent or whitespace-only."""
    return value is None or not str(value).strip()


def _err(record, code, message=None, field=None):
    severity, default_field, default_message = RULES[code]
    return ValidationError(
        code=code,
        field=field if field is not None else default_field,
        message=message if message is not None else default_message,
        severity=severity,
        record_id=getattr(record, "record_id", None),
    )


def validate(record: TierARecord) -> List[ValidationError]:
    """Return every finding about ``record``, hard and warning alike.

    Order is stable: anti-fabrication, then structural, then value sanity,
    then warnings. An empty list means the record is admissible with no
    caveats. The record is never modified.
    """
    out: List[ValidationError] = []
    quote = record.support_quote
    quote_blank = _blank(quote)

    # --- A: anti-fabrication (SCHEMA.md 2) -------------------------------
    if quote_blank:
        out.append(_err(record, "A001_MISSING_SUPPORT_QUOTE"))
    if _blank(record.source_locator):
        out.append(_err(record, "A002_MISSING_SOURCE_LOCATOR"))
    if (
        record.direction is not None
        and record.direction in ASSERTIVE_DIRECTIONS
        and quote_blank
    ):
        out.append(
            _err(
                record,
                "A003_DIRECTION_WITHOUT_QUOTE",
                "direction=%s asserts an observed effect but no support_quote "
                "is present to derive it from." % record.direction.value,
            )
        )
    if record.effect_size_value is not None and (
        quote_blank or not _DIGIT.search(str(quote))
    ):
        out.append(
            _err(
                record,
                "A004_EFFECT_SIZE_UNSUPPORTED",
                "effect_size_value=%r but support_quote contains no digit, so "
                "the number is not traceable to the source."
                % (record.effect_size_value,),
            )
        )
    if _blank(record.outcome_id):
        out.append(_err(record, "A005_MISSING_OUTCOME_ID"))
    if record.evidence_class is not EvidenceClass.causal_interventional:
        out.append(
            _err(
                record,
                "A006_EVIDENCE_CLASS_NOT_CAUSAL",
                "evidence_class=%s; Tier A records must be "
                "causal_interventional."
                % (
                    record.evidence_class.value
                    if record.evidence_class is not None
                    else "null"
                ),
            )
        )
    if record.extraction_method is None:
        out.append(_err(record, "A007_MISSING_EXTRACTION_METHOD"))

    # --- S: structural ----------------------------------------------------
    if _blank(record.study_id):
        out.append(_err(record, "S001_MISSING_STUDY_ID"))
    if _blank(record.population_raw):
        out.append(_err(record, "S002_MISSING_POPULATION_RAW"))
    if _blank(record.outcome_raw):
        out.append(_err(record, "S003_MISSING_OUTCOME_RAW"))
    if record.direction is None:
        out.append(_err(record, "S004_MISSING_DIRECTION"))
    if not record.id_is_consistent():
        out.append(
            _err(
                record,
                "S005_RECORD_ID_INCONSISTENT",
                "record_id=%r but identifying fields hash to %r."
                % (record.record_id, record.compute_id()),
            )
        )
    if record.schema_version != SCHEMA_VERSION:
        out.append(
            _err(
                record,
                "S006_SCHEMA_VERSION_MISMATCH",
                "schema_version=%r, this code implements %r."
                % (record.schema_version, SCHEMA_VERSION),
            )
        )

    # --- V: value sanity --------------------------------------------------
    if record.p_value is not None and not (0.0 <= float(record.p_value) <= 1.0):
        out.append(
            _err(record, "V001_P_VALUE_OUT_OF_RANGE", "p_value=%r." % (record.p_value,))
        )
    if (
        record.ci_low is not None
        and record.ci_high is not None
        and float(record.ci_low) > float(record.ci_high)
    ):
        out.append(
            _err(
                record,
                "V002_CI_BOUNDS_INVERTED",
                "ci_low=%r > ci_high=%r." % (record.ci_low, record.ci_high),
            )
        )
    for name in ("n_animals", "n_treatment", "n_control"):
        value = getattr(record, name)
        if value is not None and int(value) < 0:
            out.append(
                _err(
                    record,
                    "V003_NEGATIVE_COUNT",
                    "%s=%r." % (name, value),
                    field=name,
                )
            )
    if record.dispersion_value is not None and float(record.dispersion_value) < 0:
        out.append(
            _err(
                record,
                "V004_NEGATIVE_DISPERSION",
                "dispersion_value=%r." % (record.dispersion_value,),
            )
        )

    # --- W: warnings ------------------------------------------------------
    controls = [c for c in record.control_types if c is not ControlType.none_reported]
    if not controls:
        out.append(_err(record, "W001_NO_CONTROL_TYPE"))
    if (
        record.n_treatment is None
        and record.n_control is None
        and record.n_animals is None
    ):
        out.append(_err(record, "W002_NO_SAMPLE_SIZE"))
    if record.population_id is None:
        out.append(
            _err(
                record,
                "W003_POPULATION_ID_NULL",
                "population_raw=%r did not resolve to the vocabulary."
                % (record.population_raw,),
            )
        )
    if (
        record.full_text_available is not True
        or record.retrieved_from is RetrievedFrom.europepmc_abstract
    ):
        out.append(_err(record, "W004_ABSTRACT_ONLY"))
    if (
        record.extraction_method is not None
        and record.extraction_method.value == "llm"
        and record.review_status.value == "unreviewed"
    ):
        out.append(_err(record, "W005_UNREVIEWED_LLM_EXTRACTION"))
    if (
        record.direction is not None
        and record.direction in ASSERTIVE_DIRECTIONS
        and record.direction is not Direction.no_change
        and record.effect_size_value is None
    ):
        out.append(_err(record, "W006_NO_EFFECT_SIZE"))
    if record.confidence is not None:
        expected = assign_confidence(
            record.full_text_available,
            record.has_actuator_control,
            record.n_treatment is not None or record.n_control is not None,
            record.effect_size_value is not None,
        )
        if record.confidence is not expected:
            out.append(
                _err(
                    record,
                    "W007_CONFIDENCE_RULE_MISMATCH",
                    "confidence=%s but the rule yields %s."
                    % (record.confidence.value, expected.value),
                )
            )
    if record.significant is not None and (
        _blank(record.statistical_test) and record.p_value is None
    ):
        out.append(_err(record, "W008_SIGNIFICANCE_WITHOUT_TEST"))
    if record.age_group is not None:
        expected_age = assign_age_group(
            record.age_value, record.age_unit, record.species
        )
        if record.age_group is not expected_age:
            out.append(
                _err(
                    record,
                    "W009_AGE_GROUP_RULE_MISMATCH",
                    "age_group=%s but the rule yields %s."
                    % (record.age_group.value, expected_age.value),
                )
            )
    if _blank(record.last_author) or _blank(record.last_author_affiliation):
        out.append(_err(record, "W010_NO_LAB_IDENTITY"))
    if record.outcome_tier is None:
        out.append(_err(record, "W011_NO_OUTCOME_TIER"))
    if _blank(record.title) or record.year is None:
        out.append(_err(record, "W012_MISSING_BIBLIOGRAPHY"))

    return out


def hard_findings(findings) -> List[ValidationError]:
    return [f for f in findings if f.severity is Severity.HARD]


def warning_findings(findings) -> List[ValidationError]:
    return [f for f in findings if f.severity is Severity.WARNING]


def is_admissible(record: TierARecord) -> bool:
    """True if ``record`` may enter ``data/curated/`` (no HARD findings)."""
    return not hard_findings(validate(record))
