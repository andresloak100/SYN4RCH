"""Pure, inspectable assignment rules.

Two fields in the Tier A schema are *derived*, never judged:

* ``age_group``  -- SCHEMA.md 3.3
* ``confidence`` -- SCHEMA.md 3.6

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
    Species,
)

__all__ = [
    "WEEKS_PER_UNIT",
    "age_in_weeks",
    "assign_age_group",
    "confidence_inputs_met",
    "assign_confidence",
    "has_actuator_control_from_controls",
]


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
