"""Pure, inspectable assignment rules.

Four fields in the Tier A schema are *derived*, never judged:

* ``age_group``             -- SCHEMA.md 3.3
* ``confidence``            -- SCHEMA.md 3.6
* ``perturbation_class``    -- SCHEMA.md 3.4b (v0.2.0)
* ``perturbation_polarity`` -- SCHEMA.md 3.4b (v0.2.0)

Both are implemented here as side-effect-free functions of their declared
inputs so that they can be read, tested and audited in isolation. Nothing in
this module reads a record, a file or the network.

Python 3.10 compatible, standard library only.
"""

from typing import Optional

from .enums import (
    ACTUATOR_CONTROL_TYPES,
    AgeGroup,
    AgeUnit,
    Confidence,
    ControlType,
    PerturbationClass,
    PerturbationMethod,
    PerturbationPolarity,
    Species,
    Valence,
)

__all__ = [
    "WEEKS_PER_UNIT",
    "METHOD_CLASS",
    "age_in_weeks",
    "assign_age_group",
    "confidence_inputs_met",
    "assign_confidence",
    "has_actuator_control_from_controls",
    "perturbation_class_from_method",
    "perturbation_polarity_from",
]


#: ``method`` -> ``perturbation_class``. Total over :class:`PerturbationMethod`:
#: a method absent from this table is a programming error, not a default, and
#: :func:`perturbation_class_from_method` raises rather than guessing.
#:
#: The split that matters is the first two rows against the next two. A
#: chemogenetic manipulation changes firing; a conditional knockout changes a
#: gene product and leaves firing an open question. Before v0.2.0 both were
#: described by ``valence`` alone, which is why 185 of the first 362 curated
#: records carried ``valence: unclear`` -- 170 of them molecular -- and why
#: half the dataset had no usable signed axis. See docs/DECISIONS.md (D-025).
METHOD_CLASS = {
    PerturbationMethod.optogenetic: PerturbationClass.circuit_activity,
    PerturbationMethod.chemogenetic: PerturbationClass.circuit_activity,
    PerturbationMethod.electrical_stimulation: PerturbationClass.circuit_activity,
    PerturbationMethod.conditional_knockout: PerturbationClass.molecular,
    PerturbationMethod.conditional_overexpression: PerturbationClass.molecular,
    PerturbationMethod.lesion: PerturbationClass.structural,
    PerturbationMethod.genetic_ablation: PerturbationClass.structural,
    PerturbationMethod.local_pharmacology: PerturbationClass.pharmacological,
    PerturbationMethod.other: PerturbationClass.other,
}


#: Conversion factors to weeks. ``months`` and ``years`` use the mean Gregorian
#: year (365.25 d) rather than 4- or 52-week approximations, so that
#: "3 months" -> 13.04 wk rather than 12 wk. See docs/DECISIONS.md (D-004).
WEEKS_PER_UNIT = {
    AgeUnit.days: 1.0 / 7.0,
    AgeUnit.weeks: 1.0,
    AgeUnit.months: 365.25 / 12.0 / 7.0,
    AgeUnit.years: 365.25 / 7.0,
}


def age_in_weeks(age_value, age_unit) -> Optional[float]:
    """Convert ``(age_value, age_unit)`` to weeks, or ``None`` if unreported.

    Returns ``None`` when either input is missing, when ``age_value`` is not a
    finite non-negative number, or when ``age_unit`` is not a known unit. It
    never guesses a unit.
    """
    if age_value is None or age_unit is None:
        return None
    try:
        unit = AgeUnit(age_unit)
    except ValueError:
        return None
    try:
        value = float(age_value)
    except (TypeError, ValueError):
        return None
    if value != value or value in (float("inf"), float("-inf")) or value < 0:
        return None
    return value * WEEKS_PER_UNIT[unit]


def assign_age_group(age_value, age_unit, species) -> AgeGroup:
    """Assign ``age_group`` by explicit rule (SCHEMA.md 3.3).

    Mouse (``mus_musculus``) boundaries, in weeks, exactly as specified:

    ====================  ===============
    age (weeks)           age_group
    ====================  ===============
    ``< 8``               ``juvenile``
    ``8 <= age < 27``     ``young_adult``
    ``27 <= age < 65``    ``middle_aged``
    ``>= 65``             ``aged``
    ====================  ===============

    SCHEMA.md writes the bands as "8-26", "27-64" and ">=65" weeks, i.e. in
    whole weeks. The half-open form above is the integer-identical reading of
    that table, and additionally fixes the behaviour of non-integer ages in the
    gaps (26.5 wk -> ``young_adult``, 64.5 wk -> ``middle_aged``). This
    resolution is recorded in docs/DECISIONS.md (D-003).

    **These boundaries are mouse-specific and do not transfer.** A 65-week rat
    is not an aged rat, and the mapping is meaningless for humans. Any species
    other than ``mus_musculus`` therefore returns ``not_reported``: this
    function will not extrapolate. Supporting rats or humans requires a
    separate, documented, species-specific rule plus a SCHEMA.md version bump.

    ``not_reported`` is also returned when age is missing or uninterpretable.
    """
    if species is None:
        return AgeGroup.not_reported
    try:
        sp = Species(species)
    except ValueError:
        return AgeGroup.not_reported
    if sp is not Species.mus_musculus:
        return AgeGroup.not_reported

    weeks = age_in_weeks(age_value, age_unit)
    if weeks is None:
        return AgeGroup.not_reported
    if weeks < 8.0:
        return AgeGroup.juvenile
    if weeks < 27.0:
        return AgeGroup.young_adult
    if weeks < 65.0:
        return AgeGroup.middle_aged
    return AgeGroup.aged


def has_actuator_control_from_controls(control_types) -> bool:
    """True if ``control_types`` contains a control that isolates the actuator.

    Helper for curation only. The record's own ``has_actuator_control`` field is
    authoritative and is *not* overwritten by this function.
    """
    if not control_types:
        return False
    out = set()
    for raw in control_types:
        try:
            out.add(ControlType(raw))
        except ValueError:
            continue
    return bool(out & ACTUATOR_CONTROL_TYPES)


def confidence_inputs_met(
    full_text_available,
    has_actuator_control,
    n_reported,
    effect_size_reported,
) -> dict:
    """Normalise the four confidence inputs to strict booleans.

    Exposed separately so a QC report can show *why* a record scored as it did.
    Anything other than ``True`` counts as not met (``None`` is not met).
    """
    return {
        "full_text_available": full_text_available is True,
        "has_actuator_control": has_actuator_control is True,
        "n_reported": n_reported is True,
        "effect_size_reported": effect_size_reported is True,
    }


def assign_confidence(
    full_text_available,
    has_actuator_control,
    n_reported,
    effect_size_reported,
) -> Confidence:
    """Assign ``confidence`` by documented rule (SCHEMA.md 3.6).

    SCHEMA.md names the four inputs but does not give the mapping. The mapping
    below is this project's choice and is fixed here (docs/DECISIONS.md D-005).
    Changing it requires a schema version bump, because it changes stored data.

    Let ``S = has_actuator_control + n_reported + effect_size_reported``
    (0..3), evaluated only on strict ``True``:

    * ``low``      -- ``full_text_available`` is false, **or** ``S <= 1``.
      Rationale: without the full text the record cannot be re-checked at all,
      and a single supporting attribute is not gradeable.
    * ``high``     -- ``full_text_available`` and ``S == 3`` (all four met).
    * ``moderate`` -- otherwise (full text plus exactly 2 of 3).

    The full 16-row truth table is pinned in ``tests/test_rules.py``.
    """
    met = confidence_inputs_met(
        full_text_available, has_actuator_control, n_reported, effect_size_reported
    )
    if not met["full_text_available"]:
        return Confidence.low
    score = (
        int(met["has_actuator_control"])
        + int(met["n_reported"])
        + int(met["effect_size_reported"])
    )
    if score <= 1:
        return Confidence.low
    if score == 3:
        return Confidence.high
    return Confidence.moderate


def perturbation_class_from_method(method) -> Optional[PerturbationClass]:
    """``perturbation_class`` implied by ``method`` (SCHEMA.md 3.4b).

    ``None`` in, ``None`` out: an unreported method does not get a class.
    A method that is not a :class:`PerturbationMethod` value also returns
    ``None`` rather than ``other``, because ``other`` is a curator's explicit
    declaration that the method did not fit the vocabulary, and silently
    writing it over an unrecognised string would hide a vocabulary gap.
    """
    if method is None:
        return None
    try:
        m = PerturbationMethod(method)
    except ValueError:
        return None
    try:
        return METHOD_CLASS[m]
    except KeyError:  # pragma: no cover - guarded by test_rules
        raise KeyError(
            "METHOD_CLASS is not total: %r has no perturbation_class. Add it "
            "deliberately rather than letting a default be inferred." % (m,)
        )


def perturbation_polarity_from(method, valence) -> Optional[PerturbationPolarity]:
    """``perturbation_polarity`` implied by ``(method, valence)``.

    The polarity is the sign of the manipulation **along its own class's
    axis**, and the rule therefore branches on the class first. Read
    :class:`~syn4rch.schema.enums.PerturbationPolarity` before using the
    result; a polarity without its class is meaningless and
    :func:`syn4rch.qc.keys.perturbation_key` exists so that the two cannot be
    separated.

    * ``circuit_activity`` -- the polarity restates ``valence`` on the firing
      axis: ``activation`` -> ``increase``, ``inhibition`` -> ``decrease``,
      ``ablation`` -> ``remove``, ``unclear`` -> ``unclear``. A circuit
      manipulation with no recorded valence stays ``unclear``: the sign is not
      recoverable from the method, since the same virus can excite or silence.
    * ``molecular`` -- the axis is the amount of the target gene product, and
      the method states it outright: ``conditional_knockout`` -> ``decrease``
      (loss of function), ``conditional_overexpression`` -> ``increase`` (gain
      of function). ``valence`` is deliberately ignored here, and this is the
      whole point of the v0.2.0 split: a molecular record's ``valence`` is
      honestly ``unclear`` and must stay so, while its *manipulation* has a
      perfectly well defined sign.
    * ``structural`` -- ``remove``. The population is gone.
    * ``pharmacological`` -- read from ``valence``; a locally infused drug can
      do either and the method does not say which.
    * ``other`` -- ``unclear``.

    ``None`` when ``method`` is missing or unrecognised. ``unclear`` is never
    used as a convenient filler: every path to it is named above.
    """
    cls = perturbation_class_from_method(method)
    if cls is None:
        return None

    if cls is PerturbationClass.molecular:
        m = PerturbationMethod(method)
        if m is PerturbationMethod.conditional_knockout:
            return PerturbationPolarity.decrease
        if m is PerturbationMethod.conditional_overexpression:
            return PerturbationPolarity.increase
        return PerturbationPolarity.unclear  # pragma: no cover

    if cls is PerturbationClass.structural:
        return PerturbationPolarity.remove

    if cls is PerturbationClass.other:
        return PerturbationPolarity.unclear

    # circuit_activity and pharmacological both read the valence.
    if valence is None:
        return PerturbationPolarity.unclear
    try:
        v = Valence(valence)
    except ValueError:
        return PerturbationPolarity.unclear
    return {
        Valence.activation: PerturbationPolarity.increase,
        Valence.inhibition: PerturbationPolarity.decrease,
        Valence.ablation: PerturbationPolarity.remove,
        Valence.unclear: PerturbationPolarity.unclear,
    }[v]
