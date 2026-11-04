"""SYN4RCH Tier A schema (v0.3.0): enums, record, derived rules, validation, IO.

``docs/SCHEMA.md`` is the binding contract; this package implements it.
"""

from .enums import (  # noqa: F401
    ACTUATOR_CONTROL_TYPES,
    ASSERTIVE_DIRECTIONS,
    Acuity,
    AgeGroup,
    AgeUnit,
    Challenge,
    ConditionSource,
    Confidence,
    ControlType,
    Diet,
    Direction,
    DispersionType,
    EffectSizeType,
    EvidenceClass,
    ExtractionMethod,
    Laterality,
    Neurotransmitter,
    OutcomeTier,
    PerturbationClass,
    PerturbationMethod,
    PerturbationPolarity,
    RetrievedFrom,
    ReviewStatus,
    Sex,
    Species,
    TimeUnit,
    Valence,
)
from .io import (  # noqa: F401
    JsonlError,
    SchemaVersionError,
    append_jsonl,
    dumps_record,
    iter_jsonl,
    loads_record,
    read_jsonl,
    write_jsonl,
)
from .record import (  # noqa: F401
    CONDITION_DIMENSIONS,
    RECORD_ID_FIELDS,
    SCHEMA_VERSION,
    ConditionEvidence,
    SchemaCoercionError,
    TierARecord,
    compute_record_id,
    normalize_for_hash,
)
from .rules import (  # noqa: F401
    METHOD_CLASS,
    age_in_weeks,
    assign_age_group,
    assign_confidence,
    confidence_inputs_met,
    has_actuator_control_from_controls,
    perturbation_class_from_method,
    perturbation_polarity_from,
)
from .validate import (  # noqa: F401
    RULES,
    Severity,
    ValidationError,
    hard_findings,
    is_admissible,
    validate,
    warning_findings,
)

__all__ = [n for n in dir() if not n.startswith("_")]
