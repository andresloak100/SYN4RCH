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
    Challenge,
    ConditionSource,
    ControlType,
    Direction,
    EvidenceClass,
    PerturbationClass,
    RetrievedFrom,
)
from .record import (
    CONDITION_DIMENSIONS,
    EVIDENCE_DIMENSIONS,
    EVIDENCE_FOR_CONDITION,
    SCHEMA_VERSION,
    TierARecord,
)
from .rules import (
    assign_age_group,
    assign_confidence,
    perturbation_class_from_method,
    perturbation_polarity_from,
)

__all__ = [
    "Severity",
    "ValidationError",
    "RULES",
    "validate",
    "hard_findings",
    "warning_findings",
    "is_admissible",
    "unevidenced_conditions",
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
    "A008_CONDITION_EVIDENCE_WITHOUT_LOCATOR": (
        Severity.HARD,
        "condition_evidence",
        "A condition evidence entry must name where its quote came from. The "
        "Methods route in particular is only admissible if it can be audited, "
        "which means naming the passage, not merely the section.",
    ),
    "A009_CONDITION_EVIDENCE_WITHOUT_QUOTE": (
        Severity.HARD,
        "condition_evidence",
        "A condition evidence entry exists because a sentence was found. "
        "Without the sentence the entry asserts that evidence exists while "
        "withholding it, which is worse than leaving the dimension null.",
    ),
    "A010_CONDITION_EVIDENCE_MALFORMED": (
        Severity.HARD,
        "condition_evidence",
        "A condition evidence entry must name one of the four condition "
        "dimensions exactly once, and its source must be passage or methods. "
        "not_reported is not a source: an absence cannot be quoted.",
    ),
    "A011_CONDITION_EVIDENCE_VALUE_MISMATCH": (
        Severity.HARD,
        "condition_evidence",
        "A condition evidence entry must agree with the field it describes. An "
        "entry that disagrees with its own field, or that supports a dimension "
        "left null or recorded as not_reported, is provenance pointing at a "
        "value that is not there.",
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
    "S007_PERTURBATION_CLASS_MISMATCH": (
        Severity.HARD,
        "perturbation_class",
        "perturbation_class is a total function of method and admits no "
        "override. A stored value that disagrees with the rule means the "
        "record cannot be keyed reliably, so it is refused rather than "
        "warned about.",
    ),
    "S008_TARGET_ON_NON_MOLECULAR": (
        Severity.HARD,
        "target_raw",
        "A molecular target belongs to a molecular or pharmacological "
        "perturbation. Setting one on a circuit or structural perturbation "
        "corrupts the perturbation key, which is a component of record_id.",
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
    "W013_PERTURBATION_POLARITY_RULE_MISMATCH": (
        Severity.WARNING,
        "perturbation_polarity",
        "Stored perturbation_polarity differs from the documented rule's "
        "output. A warning rather than a refusal: for a locally infused drug "
        "a curator who has read the pharmacology may legitimately know the "
        "sign when valence does not.",
    ),
    "W014_MOLECULAR_TARGET_UNIDENTIFIED": (
        Severity.WARNING,
        "target_id",
        "A molecular perturbation with no identified target pools with every "
        "other molecular perturbation of the same population. A knockout of "
        "Raptor and a knockout of p53 are not one perturbation, and treating "
        "them as one is how a cell comes to look contested.",
    ),
    "W015_CONDITION_UNSPECIFIED": (
        Severity.WARNING,
        "challenge",
        "None of sex, diet, age_group or challenge is known, so this record "
        "can only ever be compared with other condition-blind records, and "
        "any agreement or disagreement found that way may be condition "
        "heterogeneity instead.",
    ),
    "W016_POLARITY_MAY_BE_INVERTED_BY_ALLELE": (
        Severity.WARNING,
        "perturbation_polarity",
        "The targeting text names a dominant-negative or catalytically dead "
        "allele. perturbation_polarity is derived from the method, so "
        "overexpressing such an allele is recorded as a gain of function while "
        "it acts as a loss of function. Resolve by hand; this rule never "
        "rewrites the stored value.",
    ),
    "W017_CONDITION_WITHOUT_EVIDENCE": (
        Severity.WARNING,
        "condition_evidence",
        "A condition dimension holds a positive value but carries no evidence "
        "entry, so a reader cannot check it without re-reading the paper. The "
        "362 records curated before v0.3.0 are all in this state; the warning "
        "exists to count the remaining gap rather than to hide it.",
    ),
}

#: Wordings that mark an allele whose expression *reduces* the target's
#: function. Overexpressing one is ``conditional_overexpression`` by method and
#: therefore ``perturbation_polarity = increase`` by rule, while functionally it
#: is a loss of function -- the D-019 sign error, one level below the one
#: :class:`~syn4rch.schema.enums.PerturbationPolarity` was introduced to fix.
#:
#: This pattern only ever raises a warning. It never rewrites a polarity, for
#: the same reason the outcomes vocabulary's polarity annotations never rewrite
#: a direction (docs/DECISIONS.md D-021): a mechanical guess at a sign is how
#: the inversion gets into the data in the first place. See D-032.
_INVERTING_ALLELE = re.compile(
    r"dominant[-\s]negative"
    r"|\bdominant[-\s]neg\b"
    r"|kinase[-\s]dead"
    r"|deacetylase[-\s]dead"
    r"|catalytically[-\s](?:in|un)?active"
    r"|catalytic(?:ally)?[-\s]dead"
    r"|\bH355Y\b"
    r"|\bK\d+R\b"
    r"|\bD\d+[AN]\b"
    r"|loss[-\s]of[-\s]function\s+mutant"
    r"|inactive\s+mutant",
    re.IGNORECASE,
)

#: Perturbation classes that may carry a molecular target.
_TARGET_BEARING_CLASSES = frozenset(
    {PerturbationClass.molecular, PerturbationClass.pharmacological}
)

#: Condition values that mean "not known", for W015. Kept local rather than
#: imported from :mod:`syn4rch.qc.keys`: ``qc`` depends on ``schema`` and not
#: the other way round, and inverting that would make the schema package
#: unimportable on its own.
_CONDITION_UNKNOWN = frozenset({None, "not_reported"})


def _blank(value: Any) -> bool:
    """True if a text field is absent or whitespace-only."""
    return value is None or not str(value).strip()


def _condition_field_value(record, dimension: str):
    """The record's own value for one condition dimension, as a plain string."""
    raw = getattr(record, dimension, None)
    return getattr(raw, "value", raw)


def unevidenced_conditions(record) -> List[str]:
    """Condition dimensions holding a positive value with no evidence behind it.

    Public because the QC report counts the same population. Two
    implementations of one rule agree until the day one of them is corrected,
    so there is one.

    ``not_reported`` is excluded deliberately: it asserts an absence, and an
    absence has no sentence to quote (see :class:`ConditionEvidence`).

    ``age_group`` is answered for by an ``age_value`` entry rather than one of
    its own, because the band is derived from the age by rule and no paper
    writes the word ``young_adult``
    (:data:`~syn4rch.schema.record.EVIDENCE_FOR_CONDITION`, D-045). An age band
    is therefore auditable exactly when the age it was computed from is quoted.
    """
    evidenced = {
        e.dimension
        for e in (record.condition_evidence or [])
        if isinstance(getattr(e, "dimension", None), str)
    }
    out = []
    for dimension in CONDITION_DIMENSIONS:
        value = _condition_field_value(record, dimension)
        if value in _CONDITION_UNKNOWN:
            continue
        if EVIDENCE_FOR_CONDITION[dimension] not in evidenced:
            out.append(dimension)
    return out


def _same_value(stored: Any, claimed: Any) -> bool:
    """Does an entry's ``value`` name the same thing as the field's value?

    Stringly for enums, numerically for ``age_value``: a paper writes "12 weeks
    old", so a curator may record ``"12"`` where the field holds ``12.0``, and
    refusing that would be an equality bug dressed as an integrity rule.
    """
    if stored == claimed:
        return True
    try:
        return float(stored) == float(claimed)
    except (TypeError, ValueError):
        return False


def _check_condition_evidence(record) -> List["ValidationError"]:
    """A008-A011: every condition evidence entry must be auditable and true.

    The four rules are deliberately separate findings rather than one. A
    missing locator, a missing quote, a malformed entry and an entry that
    disagrees with its own field are four different mistakes with four
    different fixes, and a report that collapses them tells a curator only
    that something is wrong with the condition.

    A010 is the exception and reports alone: an entry naming no recognisable
    dimension cannot be checked for anything else, because there is no field to
    compare it against and no way to know which value it was meant to support.
    So a malformed entry that is *also* missing its quote reports one finding
    rather than two. That is deliberate, and it is the right order of work --
    identify the dimension first, then everything else becomes checkable.
    """
    out: List["ValidationError"] = []
    entries = record.condition_evidence or []
    seen = set()
    for entry in entries:
        dimension = getattr(entry, "dimension", None)
        source = getattr(entry, "source", None)
        label = dimension if isinstance(dimension, str) and dimension else "?"
        malformed = (
            dimension not in EVIDENCE_DIMENSIONS
            or dimension in seen
            or source is None
            or source is ConditionSource.not_reported
        )
        if isinstance(dimension, str):
            seen.add(dimension)
        if malformed:
            out.append(
                _err(
                    record,
                    "A010_CONDITION_EVIDENCE_MALFORMED",
                    "entry for %r with source=%s"
                    % (label, getattr(source, "value", source)),
                )
            )
            continue
        if _blank(getattr(entry, "locator", None)):
            out.append(
                _err(
                    record,
                    "A008_CONDITION_EVIDENCE_WITHOUT_LOCATOR",
                    "%s evidence carries no locator" % label,
                )
            )
        if _blank(getattr(entry, "quote", None)):
            out.append(
                _err(
                    record,
                    "A009_CONDITION_EVIDENCE_WITHOUT_QUOTE",
                    "%s evidence carries no quote" % label,
                )
            )
        stored = _condition_field_value(record, dimension)
        claimed = getattr(entry, "value", None)
        if stored in _CONDITION_UNKNOWN or not _same_value(stored, claimed):
            out.append(
                _err(
                    record,
                    "A011_CONDITION_EVIDENCE_VALUE_MISMATCH",
                    "%s evidence says %r but the field says %r"
                    % (label, claimed, stored),
                )
            )
    return out


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
    out.extend(_check_condition_evidence(record))

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
    expected_class = perturbation_class_from_method(record.method)
    if (
        record.perturbation_class is not None
        and expected_class is not None
        and record.perturbation_class is not expected_class
    ):
        out.append(
            _err(
                record,
                "S007_PERTURBATION_CLASS_MISMATCH",
                "perturbation_class=%s but method=%s implies %s."
                % (
                    record.perturbation_class.value,
                    record.method.value if record.method else "null",
                    expected_class.value,
                ),
            )
        )
    # ``method`` is what a curator read off the paper; ``perturbation_class``
    # is derived from it. Where the two disagree the method wins, so that a
    # wrong stored class cannot also suppress the findings below.
    effective_class = expected_class or record.perturbation_class
    if (
        not _blank(record.target_raw) or not _blank(record.target_id)
    ) and effective_class is not None and effective_class not in _TARGET_BEARING_CLASSES:
        out.append(
            _err(
                record,
                "S008_TARGET_ON_NON_MOLECULAR",
                "target_raw=%r / target_id=%r on a %s perturbation."
                % (record.target_raw, record.target_id, effective_class.value),
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
    if record.perturbation_polarity is not None:
        expected_polarity = perturbation_polarity_from(
            record.method, record.valence
        )
        if (
            expected_polarity is not None
            and record.perturbation_polarity is not expected_polarity
        ):
            out.append(
                _err(
                    record,
                    "W013_PERTURBATION_POLARITY_RULE_MISMATCH",
                    "perturbation_polarity=%s but (method=%s, valence=%s) "
                    "yields %s."
                    % (
                        record.perturbation_polarity.value,
                        record.method.value if record.method else "null",
                        record.valence.value if record.valence else "null",
                        expected_polarity.value,
                    ),
                )
            )
    if effective_class is PerturbationClass.molecular and _blank(record.target_id):
        out.append(
            _err(
                record,
                "W014_MOLECULAR_TARGET_UNIDENTIFIED",
                "method=%s with no target_id; targeting_strategy=%r."
                % (
                    record.method.value if record.method else "null",
                    (record.targeting_strategy or "")[:80],
                ),
            )
        )
    condition_values = [
        record.sex,
        record.diet,
        record.age_group,
        record.challenge,
    ]
    if all(
        (v is None or getattr(v, "value", v) in _CONDITION_UNKNOWN)
        for v in condition_values
    ):
        out.append(_err(record, "W015_CONDITION_UNSPECIFIED"))
    unevidenced = unevidenced_conditions(record)
    if unevidenced:
        out.append(
            _err(
                record,
                "W017_CONDITION_WITHOUT_EVIDENCE",
                "no evidence entry for %s" % ", ".join(unevidenced),
            )
        )
    if effective_class is PerturbationClass.molecular:
        targeting_text = " ".join(
            t for t in (record.targeting_strategy, record.actuator, record.target_raw)
            if t
        )
        allele = _INVERTING_ALLELE.search(targeting_text)
        if allele is not None:
            out.append(
                _err(
                    record,
                    "W016_POLARITY_MAY_BE_INVERTED_BY_ALLELE",
                    "targeting text matches %r, which marks an allele that "
                    "reduces the target's function; stored "
                    "perturbation_polarity is %s and may describe the "
                    "direction of expression rather than of function."
                    % (
                        allele.group(0),
                        record.perturbation_polarity.value
                        if record.perturbation_polarity is not None
                        else "null",
                    ),
                )
            )

    return out


def hard_findings(findings) -> List[ValidationError]:
    return [f for f in findings if f.severity is Severity.HARD]


def warning_findings(findings) -> List[ValidationError]:
    return [f for f in findings if f.severity is Severity.WARNING]


def is_admissible(record: TierARecord) -> bool:
    """True if ``record`` may enter ``data/curated/`` (no HARD findings)."""
    return not hard_findings(validate(record))
