"""The Tier A record: one (population x perturbation x readout x study) observation.

Field-for-field implementation of SCHEMA.md 3.1-3.6 (v0.2.0).

Design note -- why almost every field has a default of ``None``:
``TierARecord`` is a *container*, not a gate. Rule 4 of SCHEMA.md ("if a field
is not reported, it is null") means a half-reported observation must be
representable so that it can be inspected and rejected. Admission to
``data/curated/`` is decided by :func:`syn4rch.schema.validate.validate`, which
is where the anti-fabrication rules live. Constructing an invalid record is
deliberately possible; *storing* one is not.

Python 3.10 compatible, standard library only.
"""

import hashlib
import unicodedata
from dataclasses import dataclass, field, fields
from typing import Any, Dict, List, Optional

from .enums import (
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

__all__ = [
    "SCHEMA_VERSION",
    "RECORD_ID_FIELDS",
    "RECORD_ID_LENGTH",
    "SchemaCoercionError",
    "TierARecord",
    "compute_record_id",
    "normalize_for_hash",
]

#: The only schema version this module reads or writes.
#:
#: 0.2.0 added the condition dimension and split the perturbation axis. Every
#: ``record_id`` in the corpus changed, because the hash input set changed.
#: See docs/DECISIONS.md (D-023 .. D-026).
SCHEMA_VERSION = "0.2.0"

#: Length (hex characters) of ``record_id``. 16 hex chars = 64 bits; with a
#: corpus of ~10^5 records the collision probability is ~3e-10.
RECORD_ID_LENGTH = 16

#: Byte that joins hash components. U+001F (unit separator) is not legal in
#: scientific prose, so no component value can forge a boundary.
_HASH_SEP = "\x1f"

#: Record attributes that feed ``record_id``, in fixed order.
#:
#: SCHEMA.md 3.1 specifies the hash over
#: ``(study_id, population_raw, perturbation, outcome_raw, source_locator)``.
#: ``perturbation`` is a *section* of the schema (3.4), not a field, so an
#: identifying subset is named here instead. See docs/DECISIONS.md (D-002) for
#: the original choice and (D-026) for the v0.2.0 extension.
#:
#: v0.2.0 adds five components, for one reason: two records that differ in any
#: of them are two different observations, and under v0.1.0 they collided.
#: Three id collisions in curation batch 04 were caused exactly this way.
#:
#: * ``target_raw`` -- a knockout of *Raptor* and a knockout of *p53* in the
#:   same neurons are not one perturbation. Free text is admitted here, against
#:   the general rule stated below, because the alternative is a certain
#:   collision rather than a hypothetical re-identification, and a gene symbol
#:   is about as stable as free text gets.
#: * ``sex``, ``diet``, ``age_group``, ``challenge`` -- the condition the
#:   experiment ran under. A male chow-fed cohort and a female high-fat cohort
#:   are two observations and must not share an identity.
#:
#: Still excluded: ``targeting_strategy`` and every other prose field, because
#: re-wording one during curation would silently re-identify the record and
#: orphan every ``conflicts_with`` reference to it. Also excluded:
#: ``perturbation_class`` and ``perturbation_polarity``, which are functions of
#: ``method`` and ``valence`` and so add no discriminating power; and
#: ``time_to_readout_value``, because a float in a hash is a portability
#: hazard and the condition-qualifier convention (D-022) already distinguishes
#: timepoints through ``outcome_raw``.
RECORD_ID_FIELDS = (
    "study_id",
    "population_raw",
    "method",
    "valence",
    "actuator",
    "target_raw",
    "sex",
    "diet",
    "age_group",
    "challenge",
    "outcome_raw",
    "source_locator",
)


class SchemaCoercionError(ValueError):
    """A stored value is not admissible for its field (e.g. unknown enum)."""


def normalize_for_hash(value: Any) -> str:
    """Canonicalise one hash component.

    Steps, in order: ``None`` -> ``""``; enum -> its value; ``str()``; Unicode
    NFKC; collapse all whitespace runs to a single space; strip.

    Case is **preserved**: ``record_id`` identifies a raw string as the paper
    wrote it, and conflating ``AgRP`` with ``agrp`` is a vocabulary question,
    not an identity question. NFKC removes the platform-dependent part (e.g.
    composed vs decomposed accents from different PDF extractors) so the digest
    is stable across machines.
    """
    if value is None:
        return ""
    if hasattr(value, "value") and isinstance(getattr(value, "value"), str):
        value = value.value
    text = unicodedata.normalize("NFKC", str(value))
    return " ".join(text.split())


def compute_record_id(
    study_id: Any,
    population_raw: Any,
    method: Any,
    valence: Any,
    actuator: Any,
    target_raw: Any,
    sex: Any,
    diet: Any,
    age_group: Any,
    challenge: Any,
    outcome_raw: Any,
    source_locator: Any,
) -> str:
    """Deterministic ``record_id`` (SCHEMA.md 3.1, v0.2.0).

    ``sha256`` over the NFKC-normalised components joined by U+001F, truncated
    to :data:`RECORD_ID_LENGTH` hex characters. Pure function of its arguments:
    no clock, no randomness, no dict ordering, no locale. ``PYTHONHASHSEED``
    does not affect it (``hash()`` is not used).

    The parameters are positional and their order is the contract. It matches
    :data:`RECORD_ID_FIELDS` exactly, and a mismatch between the two is caught
    by ``tests/test_record_id.py`` rather than left to be noticed later.
    """
    payload = _HASH_SEP.join(
        normalize_for_hash(v)
        for v in (
            study_id,
            population_raw,
            method,
            valence,
            actuator,
            target_raw,
            sex,
            diet,
            age_group,
            challenge,
            outcome_raw,
            source_locator,
        )
    )
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    return digest[:RECORD_ID_LENGTH]


def _coerce_enum(enum_cls, value, field_name):
    if value is None:
        return None
    if isinstance(value, enum_cls):
        return value
    try:
        return enum_cls(value)
    except ValueError:
        raise SchemaCoercionError(
            "field %r: %r is not a valid %s (allowed: %s)"
            % (
                field_name,
                value,
                enum_cls.__name__,
                ", ".join(m.value for m in enum_cls),
            )
        )


def _coerce_enum_list(enum_cls, value, field_name):
    if value is None:
        return []
    if isinstance(value, (str, bytes)):
        raise SchemaCoercionError(
            "field %r expects a list of %s, got a bare string %r"
            % (field_name, enum_cls.__name__, value)
        )
    return [_coerce_enum(enum_cls, v, field_name) for v in value]


def _coerce_str_list(value, field_name):
    if value is None:
        return []
    if isinstance(value, (str, bytes)):
        raise SchemaCoercionError(
            "field %r expects a list of strings, got a bare string %r"
            % (field_name, value)
        )
    return [str(v) for v in value]


#: Fields holding a single enum: name -> enum class.
_ENUM_FIELDS = {
    "retrieved_from": RetrievedFrom,
    "extraction_method": ExtractionMethod,
    "review_status": ReviewStatus,
    "species": Species,
    "sex": Sex,
    "age_unit": AgeUnit,
    "age_group": AgeGroup,
    "diet": Diet,
    "challenge": Challenge,
    "condition_source": ConditionSource,
    "method": PerturbationMethod,
    "valence": Valence,
    "perturbation_class": PerturbationClass,
    "perturbation_polarity": PerturbationPolarity,
    "acuity": Acuity,
    "duration_unit": TimeUnit,
    "time_to_readout_unit": TimeUnit,
    "laterality": Laterality,
    "outcome_tier": OutcomeTier,
    "direction": Direction,
    "effect_size_type": EffectSizeType,
    "dispersion_type": DispersionType,
    "evidence_class": EvidenceClass,
    "confidence": Confidence,
}

#: Fields holding a list of enums.
_ENUM_LIST_FIELDS = {
    "neurotransmitter": Neurotransmitter,
    "control_types": ControlType,
}

#: Fields holding a list of plain strings.
_STR_LIST_FIELDS = ("marker_genes", "conflicts_with")


@dataclass
class TierARecord:
    """One Tier A observation. See SCHEMA.md 3.1-3.6."""

    # --- 3.1 identity and provenance --------------------------------------
    study_id: Optional[str] = None
    record_id: Optional[str] = None
    pmid: Optional[str] = None
    pmcid: Optional[str] = None
    doi: Optional[str] = None
    title: Optional[str] = None
    journal: Optional[str] = None
    year: Optional[int] = None
    first_author: Optional[str] = None
    last_author: Optional[str] = None
    last_author_affiliation: Optional[str] = None
    retrieved_from: Optional[RetrievedFrom] = None
    full_text_available: Optional[bool] = None
    schema_version: str = SCHEMA_VERSION
    extraction_method: Optional[ExtractionMethod] = None
    extracted_by: Optional[str] = None
    extracted_at: Optional[str] = None
    review_status: ReviewStatus = ReviewStatus.unreviewed

    # --- 3.2 neuronal population ------------------------------------------
    population_raw: Optional[str] = None
    population_id: Optional[str] = None
    marker_genes: List[str] = field(default_factory=list)
    region_raw: Optional[str] = None
    region_id: Optional[str] = None
    neurotransmitter: List[Neurotransmitter] = field(default_factory=list)
    projection_specific: Optional[bool] = None
    projection_target_raw: Optional[str] = None
    projection_target_id: Optional[str] = None

    # --- 3.3 subject -------------------------------------------------------
    species: Optional[Species] = None
    strain: Optional[str] = None
    genetic_background: Optional[str] = None
    sex: Optional[Sex] = None
    age_value: Optional[float] = None
    age_unit: Optional[AgeUnit] = None
    age_group: Optional[AgeGroup] = None
    diet: Optional[Diet] = None
    n_animals: Optional[int] = None
    # condition dimension, v0.2.0 (SCHEMA.md 3.3b, DECISIONS.md D-023/D-024)
    challenge: Optional[Challenge] = None
    challenge_detail: Optional[str] = None
    condition_source: Optional[ConditionSource] = None
    condition_locator: Optional[str] = None

    # --- 3.4 perturbation --------------------------------------------------
    method: Optional[PerturbationMethod] = None
    actuator: Optional[str] = None
    valence: Optional[Valence] = None
    # perturbation axis, v0.2.0 (SCHEMA.md 3.4b, DECISIONS.md D-025)
    perturbation_class: Optional[PerturbationClass] = None
    perturbation_polarity: Optional[PerturbationPolarity] = None
    target_raw: Optional[str] = None
    target_id: Optional[str] = None
    targeting_strategy: Optional[str] = None
    acuity: Optional[Acuity] = None
    duration_value: Optional[float] = None
    duration_unit: Optional[TimeUnit] = None
    stim_frequency_hz: Optional[float] = None
    stim_pulse_width_ms: Optional[float] = None
    stim_duty_cycle: Optional[float] = None
    stim_power_mw: Optional[float] = None
    light_wavelength_nm: Optional[float] = None
    ligand: Optional[str] = None
    ligand_dose_value: Optional[float] = None
    ligand_dose_unit: Optional[str] = None
    ligand_route: Optional[str] = None
    laterality: Optional[Laterality] = None
    time_to_readout_value: Optional[float] = None
    time_to_readout_unit: Optional[TimeUnit] = None

    # --- 3.5 readout and effect -------------------------------------------
    outcome_raw: Optional[str] = None
    outcome_id: Optional[str] = None
    outcome_tier: Optional[OutcomeTier] = None
    tissue_raw: Optional[str] = None
    tissue_id: Optional[str] = None
    assay: Optional[str] = None
    assay_units: Optional[str] = None
    direction: Optional[Direction] = None
    effect_size_value: Optional[float] = None
    effect_size_type: Optional[EffectSizeType] = None
    effect_size_units: Optional[str] = None
    dispersion_type: Optional[DispersionType] = None
    dispersion_value: Optional[float] = None
    ci_low: Optional[float] = None
    ci_high: Optional[float] = None
    n_treatment: Optional[int] = None
    n_control: Optional[int] = None
    statistical_test: Optional[str] = None
    p_value: Optional[float] = None
    significant: Optional[bool] = None
    is_primary_outcome: Optional[bool] = None

    # --- 3.6 controls and evidence grading --------------------------------
    control_types: List[ControlType] = field(default_factory=list)
    has_actuator_control: Optional[bool] = None
    evidence_class: Optional[EvidenceClass] = None
    confidence: Optional[Confidence] = None
    conflicts_with: List[str] = field(default_factory=list)
    support_quote: Optional[str] = None
    source_locator: Optional[str] = None
    extraction_notes: Optional[str] = None

    # ----------------------------------------------------------------------

    def __post_init__(self) -> None:
        for name, enum_cls in _ENUM_FIELDS.items():
            setattr(
                self, name, _coerce_enum(enum_cls, getattr(self, name), name)
            )
        for name, enum_cls in _ENUM_LIST_FIELDS.items():
            setattr(
                self, name, _coerce_enum_list(enum_cls, getattr(self, name), name)
            )
        for name in _STR_LIST_FIELDS:
            setattr(self, name, _coerce_str_list(getattr(self, name), name))
        if not self.record_id:
            self.record_id = self.compute_id()

    def compute_id(self) -> str:
        """``record_id`` this record's identifying fields imply, recomputed."""
        return compute_record_id(
            *(getattr(self, name) for name in RECORD_ID_FIELDS)
        )

    def id_is_consistent(self) -> bool:
        """True if the stored ``record_id`` matches the recomputed one."""
        return self.record_id == self.compute_id()

    # --- serialisation -----------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        """Plain JSON-serialisable dict. Key order follows field declaration."""
        out: Dict[str, Any] = {}
        for f in fields(self):
            value = getattr(self, f.name)
            if isinstance(value, list):
                out[f.name] = [
                    v.value if hasattr(v, "value") else v for v in value
                ]
            elif hasattr(value, "value"):
                out[f.name] = value.value
            else:
                out[f.name] = value
        return out

    @classmethod
    def field_names(cls):
        return tuple(f.name for f in fields(cls))

    @classmethod
    def from_dict(cls, data: Dict[str, Any], strict: bool = True) -> "TierARecord":
        """Build a record from a decoded JSON object.

        Unknown keys raise :class:`SchemaCoercionError` when ``strict`` (the
        default): a silently dropped key is how a field quietly stops being
        curated. Pass ``strict=False`` to ignore them.
        """
        if not isinstance(data, dict):
            raise SchemaCoercionError(
                "expected a JSON object, got %s" % type(data).__name__
            )
        known = set(cls.field_names())
        unknown = sorted(set(data) - known)
        if unknown and strict:
            raise SchemaCoercionError(
                "unknown field(s) %s; schema_version=%r"
                % (", ".join(repr(u) for u in unknown), data.get("schema_version"))
            )
        kwargs = {k: v for k, v in data.items() if k in known}
        return cls(**kwargs)
