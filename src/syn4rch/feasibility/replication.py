"""Independent replication and direction agreement, per claim.

Answers feasibility questions 3 and 4. The unit is the *claim*, and what a
claim is changed in schema v0.2.0.

What the claim is now
---------------------
:func:`syn4rch.qc.keys.claim_key` is the key, reused so that "the same claim"
means here what it means to conflict detection. Under v0.2.0 it carries nine
components:

    (population, perturbation_class, perturbation_polarity,
     perturbation_target, outcome_id, sex, diet, age_group, challenge)

``valence`` has left the key. For a circuit manipulation it is restated by
``perturbation_polarity``; for a molecular one it is honestly ``unclear`` and
carries no information (SCHEMA.md 3.4 forbids coercing it). The class, the
polarity and the target always travel together, because a polarity means
different physical things in different classes -- ``decrease`` is less firing
under ``circuit_activity`` and less of a gene product under ``molecular`` --
and :func:`syn4rch.qc.keys.perturbation_key` is what makes keying on polarity
alone awkward on purpose.

**Both keys are reported, side by side.** The strict key is the one this
module groups by: only comparable experiments can replicate. But the
difference between a count under the strict key and the same count under
:func:`syn4rch.qc.keys.claim_key_pooled` -- the exact v0.1.0 key -- *is* the
measurement of how much apparent disagreement in this corpus was condition
heterogeneity rather than laboratories contradicting each other. Reporting
only the new number would make the v0.2.0 change unfalsifiable, so
:func:`compare_claim_keys` computes both censuses and the report prints them
together. A pooled cell is never used to decide whether two records conflict.

Condition specificity travels with every figure
-----------------------------------------------
:func:`syn4rch.qc.keys.condition_key` pools "the curator did not say" with
"the paper does not say" into ``unknown``, because both are the same statement
about the experiment. That pooling is what keeps conflict and replication
detectable at all in a corpus where 153 of 362 records state none of the four
conditions -- and it means a cell can match on *"both unknown"* rather than on
*"both the same"*.

A claim replicated across three studies that each specify sex, diet and age is
therefore a different piece of evidence from one replicated across three
studies that specify nothing, and one number printed for both is overclaiming.
This module makes that structural rather than a footnote:

* every :class:`ReplicationCell` carries ``min_condition_specificity``, the
  *minimum* over its records -- the weakest link, since one unspecified record
  is enough to make the match a match on "unknown";
* an agreement fraction is not a float. It is an :class:`Agreement`, which
  **cannot be constructed without its specificity** (:class:`SpecificityMissing`
  is raised) and whose ``__str__`` always prints it. There is no ``__float__``,
  so the bare number cannot be formatted by accident -- the same device
  :class:`syn4rch.feasibility.splits.GroupedMetric` uses for a grouping label;
* :func:`summarize_replication` reports replication and mean agreement
  *broken down by* minimum specificity, not only pooled.

Independence
------------
Counted at two levels, and they are not interchangeable:

* **study level** -- distinct ``study_id``.
* **lab level** -- distinct ``last_author`` + ``last_author_affiliation``
  (SCHEMA.md 5.2). SCHEMA.md 5 names repeat publication by one group as the
  largest leakage risk in this corpus, so three papers from one lab are
  counted as three studies and **one** lab, and the gate that matters keys on
  the lab count.

Records whose lab cannot be determined are counted in
``n_records_without_lab`` rather than merged into one anonymous lab
(docs/DECISIONS.md D-009): a split that cannot group a record's lab has to say
so.

**Agreement.** The denominator is the number of independent units that state a
direction at all; the numerator is how many of them state the modal direction.

* One observation is **not** 100% agreement. With fewer than two independent
  units the quantity is undefined and reported as ``None``. A cell with one
  study and a cell where four studies all agree are not the same evidence, and
  printing ``1.00`` for both is the single most misleading number this module
  could produce.
* ``no_change`` is a direction like any other and is never dropped (SCHEMA.md
  3.5: dropping null results is a leading source of direction bias). It can be
  the modal direction.
* ``not_reported`` is excluded from the agreement *denominator*, because it
  contradicts nothing (docs/DECISIONS.md D-011) and there is no direction to
  agree with. It is counted and reported, never silently skipped. In practice
  such records are already held out by
  :func:`syn4rch.feasibility.coverage.has_measured_outcome`.
* A study that reports two different directions for one claim is recorded as
  :data:`MIXED`. It stays in the denominator and can never match the modal
  direction, so internal disagreement lowers agreement instead of being
  resolved away.

**Why a Wilson score interval.** An agreement fraction is a binomial
proportion on a very small n, and this corpus's typical replicated cell is 2 or
3 studies that all agree. The textbook Wald interval
``p +- z*sqrt(p(1-p)/n)`` has width exactly zero at ``p = 1``, so it would
report "100% agreement, CI [1.00, 1.00]" from two papers -- certainty
manufactured by arithmetic. The Wilson score interval does not degenerate at
the boundaries, needs nothing but a normal quantile, and
``statistics.NormalDist().inv_cdf`` is in the standard library. The exact
Clopper-Pearson interval would be preferable on principle but requires an
incomplete beta quantile, which the standard library does not provide
(docs/DECISIONS.md D-001 forbids scipy). The interval is still an
approximation on n in the single digits and is reported as a width, never as a
test.

Python 3.10 compatible, standard library only.
"""

import math
import statistics
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from ..qc.conflicts import directions_conflict
from ..qc.keys import (
    CLAIM_KEY_FIELDS,
    CONDITION_FIELDS,
    POOLED_CLAIM_KEY_FIELDS,
    UNKNOWN,
    claim_key,
    claim_key_pooled,
    condition_specificity,
    lab_key,
    norm,
    study_key,
)
from ..schema.enums import ASSERTIVE_DIRECTIONS, Direction

__all__ = [
    "CONFIDENCE_LEVEL",
    "KEY_STRICT",
    "KEY_POOLED",
    "KEY_FIELDS",
    "KEY_FUNCTIONS",
    "MIXED",
    "MIN_UNITS_FOR_AGREEMENT",
    "N_CONDITION_SLOTS",
    "Agreement",
    "KeyCensus",
    "KeyComparison",
    "LOSS_PERTURBATION",
    "LOSS_STATED_DIFFERENCE",
    "LOSS_UNKNOWN_VS_STATED",
    "ReplicationCell",
    "SpecificityMissing",
    "build_replication_cells",
    "compare_claim_keys",
    "summarize_replication",
    "wilson_interval",
]

#: Two-sided confidence level for every interval this module reports.
CONFIDENCE_LEVEL = 0.95

#: ``z`` for :data:`CONFIDENCE_LEVEL`, from the standard library's normal
#: quantile rather than a pasted 1.96, so the level and the constant cannot
#: drift apart.
_Z = statistics.NormalDist().inv_cdf(0.5 + CONFIDENCE_LEVEL / 2.0)

#: Direction recorded for one study (or one lab) that reports more than one
#: direction for the same claim. Not a schema value: a bookkeeping marker for
#: internal disagreement, which is never resolved into a single direction.
MIXED = "mixed"

#: Independent units required before an agreement fraction means anything.
#: With one unit there is nothing to agree with, so agreement is ``None``.
MIN_UNITS_FOR_AGREEMENT = 2

#: Name of the schema v0.2.0 claim key. The one this module groups by.
KEY_STRICT = "strict"

#: Name of the schema v0.1.0 claim key, kept for the side-by-side comparison
#: only (see the module docstring).
KEY_POOLED = "pooled"

#: Component names of each key variant, for labelling every printed cell.
KEY_FIELDS: Dict[str, Tuple[str, ...]] = {
    KEY_STRICT: CLAIM_KEY_FIELDS,
    KEY_POOLED: POOLED_CLAIM_KEY_FIELDS,
}

#: The key function behind each variant. Both come from
#: :mod:`syn4rch.qc.keys`; neither is re-derived here.
KEY_FUNCTIONS = {KEY_STRICT: claim_key, KEY_POOLED: claim_key_pooled}

#: How many condition slots a fully specified record fills. Read from
#: :data:`syn4rch.qc.keys.CONDITION_FIELDS` so that adding a fifth condition to
#: the schema cannot leave a stale "out of 4" in a report.
N_CONDITION_SLOTS = len(CONDITION_FIELDS)


class SpecificityMissing(ValueError):
    """An agreement figure was built without its condition specificity.

    SCHEMA.md 5 makes an unlabelled cross-validated number an invalid result.
    An agreement fraction without the specificity of the cell it came from is
    the same kind of defect one level down: it does not say whether the studies
    matched on the same conditions or on both conditions being unknown. So it
    is an error, not a warning.
    """


def wilson_interval(
    successes: int, trials: int, z: float = _Z
) -> Optional[Tuple[float, float]]:
    """Wilson score interval for ``successes/trials``, clipped to ``[0, 1]``.

    Returns ``None`` when ``trials <= 0``. See the module docstring for why
    Wilson and not Wald or Clopper-Pearson.
    """
    if trials <= 0:
        return None
    if successes < 0 or successes > trials:
        raise ValueError(
            "successes=%r out of range for trials=%r" % (successes, trials)
        )
    p = successes / trials
    denominator = 1.0 + (z * z) / trials
    center = (p + (z * z) / (2.0 * trials)) / denominator
    half = (z / denominator) * math.sqrt(
        (p * (1.0 - p) / trials) + (z * z) / (4.0 * trials * trials)
    )
    return (max(0.0, center - half), min(1.0, center + half))


# --------------------------------------------------------------------------- #
# the agreement figure, which cannot be stated without its specificity
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Agreement:
    """A direction-agreement fraction bound to the evidence it describes.

    Deliberately not a ``float``. ``min_condition_specificity`` is required
    and validated in :meth:`__post_init__`, so the number cannot be built
    without saying how much of the condition was actually known in the cell it
    came from; and there is no ``__float__`` or ``__index__``, so a caller that
    wants the bare number has to reach for :attr:`value` and will be reading
    :attr:`min_condition_specificity` on the same line.

    ``level`` names the independent unit the fraction is over (``studies`` or
    ``labs``). Three studies from one lab and three studies from three labs
    produce the same fraction and are not the same evidence, which is why the
    level is part of the object rather than context the caller remembers.
    """

    value: float
    level: str
    n_agreeing: int
    n_units: int
    min_condition_specificity: int
    modal_direction: Optional[str] = None
    modal_direction_is_tied: bool = False
    ci: Optional[Tuple[float, float]] = None
    n_condition_slots: int = N_CONDITION_SLOTS
    confidence_level: float = CONFIDENCE_LEVEL
    interval_method: str = "wilson_score"

    def __post_init__(self) -> None:
        if self.min_condition_specificity is None:
            raise SpecificityMissing(
                "an agreement fraction may not be built without the minimum "
                "condition specificity of its cell: without it the number "
                "does not say whether the studies matched on the same "
                "conditions or on both being unknown"
            )
        if not 0 <= self.min_condition_specificity <= self.n_condition_slots:
            raise SpecificityMissing(
                "min_condition_specificity=%r is outside 0..%d"
                % (self.min_condition_specificity, self.n_condition_slots)
            )
        if self.level not in ("studies", "labs"):
            raise SpecificityMissing(
                "level must name the independent unit ('studies' or 'labs'), "
                "not %r" % (self.level,)
            )

    @property
    def matched_on_unknown_conditions(self) -> bool:
        """True when at least one record in the cell specifies no condition.

        The cell's records were grouped together partly because their
        conditions were unknown, not because they were the same. The agreement
        is still a measurement; it is a measurement of a weaker thing.
        """
        return self.min_condition_specificity == 0

    @property
    def specificity_label(self) -> str:
        text = "min condition specificity %d/%d" % (
            self.min_condition_specificity,
            self.n_condition_slots,
        )
        if self.matched_on_unknown_conditions:
            return text + " (matched on 'both unknown', not on 'both the same')"
        return text

    def __str__(self) -> str:
        ci = "" if self.ci is None else " %d%% CI [%.2f, %.2f]" % (
            int(round(100 * self.confidence_level)),
            self.ci[0],
            self.ci[1],
        )
        return "%.2f%s [%s=%d of %d, %s]" % (
            self.value,
            ci,
            self.level,
            self.n_agreeing,
            self.n_units,
            self.specificity_label,
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "value": self.value,
            "level": self.level,
            "n_agreeing": self.n_agreeing,
            "n_units": self.n_units,
            "modal_direction": self.modal_direction,
            "modal_direction_is_tied": self.modal_direction_is_tied,
            "ci": list(self.ci) if self.ci else None,
            "min_condition_specificity": self.min_condition_specificity,
            "n_condition_slots": self.n_condition_slots,
            "matched_on_unknown_conditions": self.matched_on_unknown_conditions,
            "confidence_level": self.confidence_level,
            "interval_method": self.interval_method,
        }


def _unit_direction(group: Sequence) -> Optional[str]:
    """The direction one independent unit (study or lab) states for a claim.

    ``None`` means the unit states no direction at all (every record
    ``not_reported`` or null). :data:`MIXED` means it states more than one.
    """
    stated = {
        r.direction.value for r in group if r.direction in ASSERTIVE_DIRECTIONS
    }
    if not stated:
        return None
    if len(stated) == 1:
        return next(iter(stated))
    return MIXED


def _modal(directions: Sequence[str]) -> Tuple[Optional[str], bool, int, int]:
    """``(modal, tied, n_agreeing, n_units)`` for one level.

    ``directions`` is one entry per independent unit; :data:`MIXED` entries
    count in ``n_units`` but can never be modal or agreeing.
    """
    n_units = len(directions)
    counts: Dict[str, int] = {}
    for value in directions:
        if value == MIXED:
            continue
        counts[value] = counts.get(value, 0) + 1
    if not counts:
        return None, False, 0, n_units
    top = max(counts.values())
    winners = sorted(k for k, v in counts.items() if v == top)
    # Deterministic choice; ``tied`` flags that the choice is arbitrary.
    return winners[0], len(winners) > 1, top, n_units


def _agreement(
    directions: Sequence[str], level: str, min_specificity: int
) -> Tuple[Optional[str], bool, Optional[Agreement], int, int]:
    """``(modal, tied, agreement, n_agreeing, n_units)``.

    ``agreement`` is ``None`` when fewer than :data:`MIN_UNITS_FOR_AGREEMENT`
    units state a direction, or when no unit states a single direction. It is
    never a float: see :class:`Agreement`.
    """
    modal, tied, n_agreeing, n_units = _modal(directions)
    if modal is None or n_units < MIN_UNITS_FOR_AGREEMENT:
        return modal, tied, None, n_agreeing, n_units
    return (
        modal,
        tied,
        Agreement(
            value=n_agreeing / n_units,
            level=level,
            n_agreeing=n_agreeing,
            n_units=n_units,
            min_condition_specificity=min_specificity,
            modal_direction=modal,
            modal_direction_is_tied=tied,
            ci=wilson_interval(n_agreeing, n_units),
        ),
        n_agreeing,
        n_units,
    )


# --------------------------------------------------------------------------- #
# one cell
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class ReplicationCell:
    """Replication and direction agreement for one claim.

    ``claim`` is the key tuple and ``claim_fields`` names its components, so a
    cell always knows which key variant it was built under and a printed cell
    can never be read as the other one. :attr:`population` and :attr:`outcome`
    are available under both variants; :attr:`perturbation` and the condition
    accessors exist only under the strict key, and :attr:`valence` only under
    the pooled one.

    ``agreement_studies`` / ``agreement_labs`` are ``None`` when undefined --
    which is the normal state of this corpus, not an error -- and are
    :class:`Agreement` objects otherwise, never bare floats.
    """

    claim: Tuple[str, ...] = ()
    claim_fields: Tuple[str, ...] = CLAIM_KEY_FIELDS
    key_variant: str = KEY_STRICT

    n_records: int = 0
    n_studies: int = 0
    n_labs: int = 0
    n_records_without_lab: int = 0
    studies: Tuple[str, ...] = ()
    labs: Tuple[str, ...] = ()

    #: Minimum :func:`syn4rch.qc.keys.condition_specificity` over the records
    #: in this cell -- the weakest link. One record specifying nothing is
    #: enough to make the whole match a match on "unknown", so the minimum,
    #: not the mean, is what qualifies every figure below.
    #:
    #: Under the **strict** key the condition is part of the key, so every
    #: record in a cell shares it and ``min == max`` by construction. Under the
    #: **pooled** key the condition is absent from the key, so one cell can
    #: hold a fully specified record next to a wholly unspecified one -- which
    #: is exactly the heterogeneity the pooled key hid, and the minimum is what
    #: makes it visible. The field is therefore not redundant: it is the one
    #: number that reads the same way under both keys.
    min_condition_specificity: int = 0
    max_condition_specificity: int = 0
    condition_specificity_counts: Dict[int, int] = field(default_factory=dict)

    direction_counts: Dict[str, int] = field(default_factory=dict)  # per record
    study_directions: Dict[str, str] = field(default_factory=dict)
    lab_directions: Dict[str, str] = field(default_factory=dict)
    n_records_direction_not_reported: int = 0

    modal_direction: Optional[str] = None
    modal_direction_is_tied: bool = False

    n_studies_stating_direction: int = 0
    n_studies_agreeing: int = 0
    agreement_studies: Optional[Agreement] = None

    n_labs_stating_direction: int = 0
    n_labs_agreeing: int = 0
    agreement_labs: Optional[Agreement] = None
    modal_direction_labs: Optional[str] = None
    modal_direction_labs_is_tied: bool = False

    #: QC conflict codes present between distinct studies in this cell, from
    #: :func:`syn4rch.qc.conflicts.directions_conflict` (docs/DECISIONS.md
    #: D-011). Reused rather than re-derived so "disagreement" has one
    #: definition in this repository.
    conflict_codes: Tuple[str, ...] = ()

    #: True when some record states ``increase`` and some states ``decrease``
    #: in this one cell, at record level. The pooled/strict comparison counts
    #: these: 16 of the 18 such cells under the v0.1.0 key were a male cohort
    #: against a female one, or a knockout against an overexpression.
    is_directly_opposed: bool = False

    outcome_tiers: Tuple[str, ...] = ()

    # --- claim components -------------------------------------------------

    @property
    def claim_dict(self) -> Dict[str, str]:
        return dict(zip(self.claim_fields, self.claim))

    @property
    def is_pooled(self) -> bool:
        return self.key_variant == KEY_POOLED

    @property
    def population(self) -> str:
        return self.claim_dict.get("population", "")

    @property
    def outcome(self) -> str:
        return self.claim_dict.get("outcome_id", "")

    @property
    def valence(self) -> Optional[str]:
        """Only the pooled (v0.1.0) key carries a valence. ``None`` otherwise."""
        return self.claim_dict.get("valence")

    @property
    def perturbation(self) -> Optional[Tuple[str, str, str]]:
        """``(class, polarity, target)``, or ``None`` under the pooled key.

        The three are returned together and never individually, for the reason
        :class:`syn4rch.schema.enums.PerturbationPolarity` gives: a polarity
        without its class is meaningless.
        """
        data = self.claim_dict
        if "perturbation_class" not in data:
            return None
        return (
            data["perturbation_class"],
            data["perturbation_polarity"],
            data["perturbation_target"],
        )

    @property
    def perturbation_label(self) -> str:
        """One printable token for the perturbation, class first."""
        parts = self.perturbation
        if parts is None:
            return self.valence or "-"
        return "%s / %s / %s" % parts

    @property
    def condition(self) -> Optional[Tuple[str, ...]]:
        data = self.claim_dict
        if any(name not in data for name in CONDITION_FIELDS):
            return None
        return tuple(data[name] for name in CONDITION_FIELDS)

    @property
    def condition_label(self) -> str:
        """The condition as ``sex=male, diet=unknown, ...``, or ``not keyed``."""
        parts = self.condition
        if parts is None:
            return "not keyed (pooled v0.1.0 key)"
        return ", ".join(
            "%s=%s" % (name, value)
            for name, value in zip(CONDITION_FIELDS, parts)
        )

    @property
    def n_condition_slots_known(self) -> int:
        """How many condition slots the *key* fixes (0 under the pooled key)."""
        parts = self.condition
        if parts is None:
            return 0
        return sum(1 for value in parts if value != UNKNOWN)

    # --- replication ------------------------------------------------------

    @property
    def is_single_observation(self) -> bool:
        return self.n_records <= 1

    @property
    def replicated_across_studies(self) -> bool:
        return self.n_studies >= MIN_UNITS_FOR_AGREEMENT

    @property
    def replicated_across_labs(self) -> bool:
        """Independent replication in the sense SCHEMA.md 5.2 cares about."""
        return self.n_labs >= MIN_UNITS_FOR_AGREEMENT

    @property
    def same_lab_only(self) -> bool:
        """Replicated across studies but all from one lab: discounted."""
        return self.replicated_across_studies and self.n_labs <= 1

    @property
    def agreement_is_defined(self) -> bool:
        return self.agreement_studies is not None

    @property
    def agreement_studies_value(self) -> Optional[float]:
        """The bare fraction, for gate arithmetic only.

        Never for reporting: a printed agreement must carry its level and its
        specificity, which is what :class:`Agreement` is for.
        """
        return None if self.agreement_studies is None else self.agreement_studies.value

    @property
    def agreement_studies_ci(self) -> Optional[Tuple[float, float]]:
        return None if self.agreement_studies is None else self.agreement_studies.ci

    @property
    def agreement_labs_value(self) -> Optional[float]:
        return None if self.agreement_labs is None else self.agreement_labs.value

    @property
    def agreement_labs_ci(self) -> Optional[Tuple[float, float]]:
        return None if self.agreement_labs is None else self.agreement_labs.ci

    @property
    def matched_on_unknown_conditions(self) -> bool:
        """True when a record in this cell specifies none of the conditions."""
        return self.min_condition_specificity == 0

    @property
    def is_contested(self) -> bool:
        """True if distinct studies state incompatible directions."""
        return bool(self.conflict_codes) or self.modal_direction_is_tied

    def to_dict(self) -> Dict[str, Any]:
        return {
            "key_variant": self.key_variant,
            "claim_fields": list(self.claim_fields),
            "claim": list(self.claim),
            "population": self.population,
            "outcome": self.outcome,
            "valence": self.valence,
            "perturbation": (
                list(self.perturbation) if self.perturbation else None
            ),
            "condition": list(self.condition) if self.condition else None,
            "n_records": self.n_records,
            "n_studies": self.n_studies,
            "n_labs": self.n_labs,
            "n_records_without_lab": self.n_records_without_lab,
            "studies": list(self.studies),
            "labs": list(self.labs),
            "min_condition_specificity": self.min_condition_specificity,
            "max_condition_specificity": self.max_condition_specificity,
            "n_condition_slots": N_CONDITION_SLOTS,
            "condition_specificity_counts": {
                str(k): v for k, v in sorted(self.condition_specificity_counts.items())
            },
            "matched_on_unknown_conditions": self.matched_on_unknown_conditions,
            "direction_counts": dict(self.direction_counts or {}),
            "study_directions": dict(self.study_directions or {}),
            "lab_directions": dict(self.lab_directions or {}),
            "n_records_direction_not_reported": self.n_records_direction_not_reported,
            "modal_direction": self.modal_direction,
            "modal_direction_is_tied": self.modal_direction_is_tied,
            "n_studies_stating_direction": self.n_studies_stating_direction,
            "n_studies_agreeing": self.n_studies_agreeing,
            "agreement_studies": (
                self.agreement_studies.to_dict() if self.agreement_studies else None
            ),
            "n_labs_stating_direction": self.n_labs_stating_direction,
            "n_labs_agreeing": self.n_labs_agreeing,
            "agreement_labs": (
                self.agreement_labs.to_dict() if self.agreement_labs else None
            ),
            "modal_direction_labs": self.modal_direction_labs,
            "conflict_codes": list(self.conflict_codes),
            "is_directly_opposed": self.is_directly_opposed,
            "replicated_across_studies": self.replicated_across_studies,
            "replicated_across_labs": self.replicated_across_labs,
            "same_lab_only": self.same_lab_only,
            "agreement_is_defined": self.agreement_is_defined,
            "is_contested": self.is_contested,
            "confidence_level": CONFIDENCE_LEVEL,
            "interval_method": "wilson_score",
            "outcome_tiers": list(self.outcome_tiers),
        }


def _conflict_codes(by_study: Dict[str, List]) -> Tuple[str, ...]:
    """Conflict codes across *distinct studies*, via the QC definition.

    Intra-study tension is excluded here exactly as
    :mod:`syn4rch.qc.conflicts` excludes it; it surfaces instead as a
    :data:`MIXED` study direction.
    """
    codes: set = set()
    studies = sorted(by_study)
    for i in range(len(studies)):
        for j in range(i + 1, len(studies)):
            for left in by_study[studies[i]]:
                for right in by_study[studies[j]]:
                    code = directions_conflict(left.direction, right.direction)
                    if code is not None:
                        codes.add(code)
    return tuple(sorted(codes))


def build_replication_cells(
    records: Iterable, key_variant: str = KEY_STRICT
) -> List[ReplicationCell]:
    """One :class:`ReplicationCell` per claim present in ``records``.

    ``key_variant`` is :data:`KEY_STRICT` by default -- the schema v0.2.0 key,
    which is the only one that may decide whether two records replicate or
    conflict. :data:`KEY_POOLED` reproduces the v0.1.0 grouping and exists for
    :func:`compare_claim_keys`; see the module docstring for why that
    comparison is kept rather than deleted.

    Only claims that occur are built: unlike the coverage matrix, there is no
    dense product here, because replication of a claim nobody made is not a
    quantity. The empty cells are :mod:`syn4rch.feasibility.coverage`'s job.

    Sorted by claim key, so the output is deterministic.
    """
    if key_variant not in KEY_FUNCTIONS:
        raise KeyError(
            "unknown claim key variant %r; declared: %s"
            % (key_variant, ", ".join(sorted(KEY_FUNCTIONS)))
        )
    key_of = KEY_FUNCTIONS[key_variant]
    fields = KEY_FIELDS[key_variant]

    groups: Dict[Tuple[str, ...], List] = {}
    for record in records:
        if not norm(record.outcome_id):
            continue  # ungroupable readout; a HARD schema finding, counted by QC
        groups.setdefault(tuple(key_of(record)), []).append(record)

    return [
        _build_cell(claim, group, key_variant, fields)
        for claim, group in sorted(groups.items())
    ]


def _build_cell(
    claim: Tuple[str, ...],
    group: Sequence,
    key_variant: str,
    fields: Tuple[str, ...],
) -> ReplicationCell:
    by_study: Dict[str, List] = {}
    by_lab: Dict[str, List] = {}
    for record in group:
        by_study.setdefault(study_key(record), []).append(record)
        lab = lab_key(record)
        if lab is not None:
            by_lab.setdefault(lab, []).append(record)

    specificities = [condition_specificity(r) for r in group]
    spec_counts: Dict[int, int] = {}
    for value in specificities:
        spec_counts[value] = spec_counts.get(value, 0) + 1
    # The weakest link, not the mean: one record stating no condition is
    # enough to make the cell a match on "unknown".
    min_specificity = min(specificities) if specificities else 0
    max_specificity = max(specificities) if specificities else 0

    direction_counts: Dict[str, int] = {}
    for record in group:
        key = record.direction.value if record.direction else "null"
        direction_counts[key] = direction_counts.get(key, 0) + 1

    study_directions: Dict[str, str] = {}
    for study, members in sorted(by_study.items()):
        value = _unit_direction(members)
        if value is not None:
            study_directions[study] = value

    lab_directions: Dict[str, str] = {}
    for lab, members in sorted(by_lab.items()):
        value = _unit_direction(members)
        if value is not None:
            lab_directions[" | ".join(lab)] = value

    modal, tied, agreement, n_agreeing, n_units = _agreement(
        list(study_directions.values()), "studies", min_specificity
    )
    (
        lab_modal,
        lab_tied,
        lab_agreement,
        n_lab_agreeing,
        n_lab_units,
    ) = _agreement(list(lab_directions.values()), "labs", min_specificity)

    return ReplicationCell(
        claim=tuple(claim),
        claim_fields=fields,
        key_variant=key_variant,
        n_records=len(group),
        n_studies=len(by_study),
        n_labs=len(by_lab),
        n_records_without_lab=sum(1 for r in group if lab_key(r) is None),
        studies=tuple(sorted(by_study)),
        labs=tuple(sorted(by_lab)),
        min_condition_specificity=min_specificity,
        max_condition_specificity=max_specificity,
        condition_specificity_counts=dict(sorted(spec_counts.items())),
        direction_counts=dict(sorted(direction_counts.items())),
        study_directions=study_directions,
        lab_directions=lab_directions,
        n_records_direction_not_reported=sum(
            1
            for r in group
            if r.direction is None or r.direction is Direction.not_reported
        ),
        modal_direction=modal,
        modal_direction_is_tied=tied,
        n_studies_stating_direction=n_units,
        n_studies_agreeing=n_agreeing,
        agreement_studies=agreement,
        n_labs_stating_direction=n_lab_units,
        n_labs_agreeing=n_lab_agreeing,
        agreement_labs=lab_agreement,
        modal_direction_labs=lab_modal,
        modal_direction_labs_is_tied=lab_tied,
        conflict_codes=_conflict_codes(by_study),
        is_directly_opposed=_directly_opposed(group),
        outcome_tiers=tuple(
            sorted({r.outcome_tier.value for r in group if r.outcome_tier})
        ),
    )


def _directly_opposed(group: Sequence) -> bool:
    """True if some record says ``increase`` and some says ``decrease`` here.

    Record level, not study level, and deliberately so: this is the quantity
    the v0.1.0-to-v0.2.0 comparison is about. A cell holding a male cohort's
    ``increase`` and a female cohort's ``decrease`` was a head-to-head
    contradiction under the pooled key whether or not the two came from
    different studies. :attr:`ReplicationCell.conflict_codes` is the separate,
    stricter, across-studies QC definition and is reported next to it.
    """
    directions = {r.direction for r in group if r.direction is not None}
    return Direction.increase in directions and Direction.decrease in directions


# --------------------------------------------------------------------------- #
# the two keys, side by side
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class KeyCensus:
    """What one claim key sees in a set of records.

    Every count here is of *cells* under that key, except ``n_records``.
    """

    key_variant: str
    fields: Tuple[str, ...]
    n_records: int = 0
    n_cells: int = 0
    n_cells_single_observation: int = 0
    n_cells_multi_study: int = 0
    n_cells_multi_lab: int = 0
    n_cells_directly_opposed: int = 0
    n_cells_with_conflict_codes: int = 0
    n_cells_agreement_defined: int = 0
    n_cells_unanimous_where_defined: int = 0
    #: cells by minimum condition specificity, 0..N_CONDITION_SLOTS
    cells_by_min_condition_specificity: Dict[int, int] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "key_variant": self.key_variant,
            "fields": list(self.fields),
            "n_records": self.n_records,
            "n_cells": self.n_cells,
            "n_cells_single_observation": self.n_cells_single_observation,
            "n_cells_multi_study": self.n_cells_multi_study,
            "n_cells_multi_lab": self.n_cells_multi_lab,
            "n_cells_directly_opposed": self.n_cells_directly_opposed,
            "n_cells_with_conflict_codes": self.n_cells_with_conflict_codes,
            "n_cells_agreement_defined": self.n_cells_agreement_defined,
            "n_cells_unanimous_where_defined": self.n_cells_unanimous_where_defined,
            "cells_by_min_condition_specificity": {
                str(k): v
                for k, v in sorted(self.cells_by_min_condition_specificity.items())
            },
        }


def _census(cells: Sequence[ReplicationCell], n_records: int, variant: str) -> KeyCensus:
    by_spec: Dict[int, int] = {}
    for cell in cells:
        by_spec[cell.min_condition_specificity] = (
            by_spec.get(cell.min_condition_specificity, 0) + 1
        )
    defined = [c for c in cells if c.agreement_is_defined]
    return KeyCensus(
        key_variant=variant,
        fields=KEY_FIELDS[variant],
        n_records=n_records,
        n_cells=len(cells),
        n_cells_single_observation=sum(1 for c in cells if c.is_single_observation),
        n_cells_multi_study=sum(1 for c in cells if c.replicated_across_studies),
        n_cells_multi_lab=sum(1 for c in cells if c.replicated_across_labs),
        n_cells_directly_opposed=sum(1 for c in cells if c.is_directly_opposed),
        n_cells_with_conflict_codes=sum(1 for c in cells if c.conflict_codes),
        n_cells_agreement_defined=len(defined),
        n_cells_unanimous_where_defined=sum(
            1 for c in defined if c.agreement_studies.value == 1.0
        ),
        cells_by_min_condition_specificity=dict(sorted(by_spec.items())),
    )


@dataclass(frozen=True)
class KeyComparison:
    """The v0.1.0 and v0.2.0 claim keys measured on the same records.

    This is the falsification test for the v0.2.0 change. If splitting the key
    by perturbation identity and condition were merely fragmenting the corpus,
    ``n_cells`` would rise and the opposed-cell count would fall by the same
    kind of arithmetic that would follow from *any* finer key. What makes the
    change a finding rather than a tautology is the ratio: how many opposed
    cells were resolved per multi-study cell lost.
    """

    pooled: KeyCensus
    strict: KeyCensus
    #: Pooled cells that the strict key splits into two or more cells.
    n_pooled_cells_split: int = 0
    #: Largest number of strict cells any single pooled cell split into.
    max_strict_cells_per_pooled_cell: int = 0
    #: Pooled cells that were directly opposed and whose strict refinements
    #: are all internally consistent: the apparent contradiction was
    #: heterogeneity, not disagreement.
    n_opposed_cells_resolved_by_split: int = 0
    #: Pooled cells still directly opposed after the split.
    n_opposed_cells_surviving: int = 0
    #: Pooled cells with two or more studies whose strict refinements all have
    #: one study: replication evidence the stricter key gives up.
    n_multi_study_cells_lost: int = 0
    #: Of those, the ones lost because the refinements differ in the
    #: *perturbation* -- a knockout against an overexpression, say. The key is
    #: right to separate these and the evidence was never there.
    n_lost_to_perturbation_split: int = 0
    #: Of those, the ones lost because two refinements state *different*
    #: values for a condition slot: male against female, chow against high
    #: fat. Also a correct separation.
    n_lost_to_stated_condition_difference: int = 0
    #: Of those, the ones lost **only** because one refinement states a
    #: condition the other leaves unknown. This is the stricter key's genuine
    #: over-split: the two experiments may well have run under the same
    #: condition, and nothing in the corpus can say whether they did. See
    #: :attr:`cost_that_cannot_be_adjudicated`.
    n_lost_to_unknown_against_stated: int = 0

    @property
    def records_per_strict_cell(self) -> Optional[float]:
        if not self.strict.n_cells:
            return None
        return self.strict.n_records / self.strict.n_cells

    @property
    def cost_that_cannot_be_adjudicated(self) -> int:
        """Lost multi-study cells where the corpus cannot say who was right.

        :data:`UNKNOWN` pools "the curator did not ask" with "the paper does
        not say", and both are honest. What neither can do is match a stated
        value: a paper that says ``male`` and a paper that says nothing are
        two different keys, so they stop replicating each other even if both
        ran male mice. That is not a defect of the key -- the alternative is
        assuming the unstated condition, which SCHEMA.md 2 rule 4 forbids --
        but it is a cost that is *not* a correction, and it is reported on its
        own line rather than counted as heterogeneity resolved.
        """
        return self.n_lost_to_unknown_against_stated

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pooled": self.pooled.to_dict(),
            "strict": self.strict.to_dict(),
            "n_pooled_cells_split": self.n_pooled_cells_split,
            "max_strict_cells_per_pooled_cell": self.max_strict_cells_per_pooled_cell,
            "n_opposed_cells_resolved_by_split": self.n_opposed_cells_resolved_by_split,
            "n_opposed_cells_surviving": self.n_opposed_cells_surviving,
            "n_multi_study_cells_lost": self.n_multi_study_cells_lost,
            "n_lost_to_perturbation_split": self.n_lost_to_perturbation_split,
            "n_lost_to_stated_condition_difference": (
                self.n_lost_to_stated_condition_difference
            ),
            "n_lost_to_unknown_against_stated": (
                self.n_lost_to_unknown_against_stated
            ),
            "cost_that_cannot_be_adjudicated": self.cost_that_cannot_be_adjudicated,
            "records_per_strict_cell": self.records_per_strict_cell,
        }


#: Why a pooled multi-study cell stopped being multi-study, worst first. The
#: first two are correct separations of work that was never comparable; the
#: third is the stricter key's own over-split and is never counted as one of
#: the first two.
LOSS_PERTURBATION = "perturbation_split"
LOSS_STATED_DIFFERENCE = "stated_condition_difference"
LOSS_UNKNOWN_VS_STATED = "unknown_against_stated"


def _loss_reason(claims: Sequence[Tuple[str, ...]]) -> str:
    """Classify why a set of strict refinements no longer share a cell.

    Priority order, and it matters: a cell whose refinements differ both in
    the perturbation and in a condition slot is attributed to the
    perturbation, because that difference alone would have separated them. The
    ``unknown_against_stated`` bucket is therefore only reached when nothing
    else explains the split, which is what makes it an honest count of the
    stricter key's own cost.
    """
    index = {name: i for i, name in enumerate(CLAIM_KEY_FIELDS)}
    perturbation_slots = [
        index[name]
        for name in (
            "perturbation_class",
            "perturbation_polarity",
            "perturbation_target",
        )
    ]
    if any(len({c[i] for c in claims}) > 1 for i in perturbation_slots):
        return LOSS_PERTURBATION
    for name in CONDITION_FIELDS:
        values = {c[index[name]] for c in claims}
        if len({v for v in values if v != UNKNOWN}) > 1:
            return LOSS_STATED_DIFFERENCE
    return LOSS_UNKNOWN_VS_STATED


def compare_claim_keys(records: Iterable) -> KeyComparison:
    """Measure both claim keys on one set of records.

    Returns the two censuses plus the three numbers that say what the v0.2.0
    key bought and what it cost: how many pooled cells it split, how many
    head-to-head contradictions turned out to be condition or perturbation
    heterogeneity, and how many multi-study cells stopped being multi-study.
    """
    records = list(records)
    keyed = [r for r in records if norm(r.outcome_id)]

    pooled_cells = build_replication_cells(records, KEY_POOLED)
    strict_cells = build_replication_cells(records, KEY_STRICT)

    # Which strict cells each pooled cell contains.
    membership: Dict[Tuple[str, ...], set] = {}
    for record in keyed:
        membership.setdefault(tuple(claim_key_pooled(record)), set()).add(
            tuple(claim_key(record))
        )
    strict_by_claim = {c.claim: c for c in strict_cells}

    split = 0
    max_fanout = 0
    resolved = 0
    surviving = 0
    lost = 0
    losses = {
        LOSS_PERTURBATION: 0,
        LOSS_STATED_DIFFERENCE: 0,
        LOSS_UNKNOWN_VS_STATED: 0,
    }
    for cell in pooled_cells:
        children = sorted(membership.get(cell.claim, ()))
        max_fanout = max(max_fanout, len(children))
        if len(children) > 1:
            split += 1
        if cell.is_directly_opposed:
            if any(strict_by_claim[c].is_directly_opposed for c in children):
                surviving += 1
            else:
                resolved += 1
        if cell.replicated_across_studies and not any(
            strict_by_claim[c].replicated_across_studies for c in children
        ):
            lost += 1
            losses[_loss_reason(children)] += 1

    return KeyComparison(
        pooled=_census(pooled_cells, len(keyed), KEY_POOLED),
        strict=_census(strict_cells, len(keyed), KEY_STRICT),
        n_pooled_cells_split=split,
        max_strict_cells_per_pooled_cell=max_fanout,
        n_opposed_cells_resolved_by_split=resolved,
        n_opposed_cells_surviving=surviving,
        n_multi_study_cells_lost=lost,
        n_lost_to_perturbation_split=losses[LOSS_PERTURBATION],
        n_lost_to_stated_condition_difference=losses[LOSS_STATED_DIFFERENCE],
        n_lost_to_unknown_against_stated=losses[LOSS_UNKNOWN_VS_STATED],
    )


# --------------------------------------------------------------------------- #
# corpus summary
# --------------------------------------------------------------------------- #

def summarize_replication(cells: Sequence[ReplicationCell]) -> Dict[str, Any]:
    """Corpus-level counts over replication cells. Measured, not estimated.

    ``mean_agreement_where_defined`` is ``None`` when no cell has a defined
    agreement, rather than ``0.0`` or ``1.0``. It is accompanied by
    ``mean_agreement_by_min_condition_specificity``, and the report prints the
    breakdown rather than the pooled mean alone: a mean taken across cells
    that matched on "both unknown" and cells that matched on three stated
    conditions describes neither.
    """
    defined = [c for c in cells if c.agreement_is_defined]
    unanimous = [c for c in defined if c.agreement_studies.value == 1.0]
    return {
        "key_variant": cells[0].key_variant if cells else KEY_STRICT,
        "claim_fields": list(cells[0].claim_fields) if cells else list(CLAIM_KEY_FIELDS),
        "n_cells": len(cells),
        "n_cells_single_observation": sum(1 for c in cells if c.is_single_observation),
        "n_cells_replicated_across_studies": sum(
            1 for c in cells if c.replicated_across_studies
        ),
        "n_cells_replicated_across_labs": sum(
            1 for c in cells if c.replicated_across_labs
        ),
        "n_cells_same_lab_only": sum(1 for c in cells if c.same_lab_only),
        "n_cells_agreement_undefined": sum(
            1 for c in cells if not c.agreement_is_defined
        ),
        "n_cells_agreement_defined": len(defined),
        "n_cells_unanimous_where_defined": len(unanimous),
        "n_cells_contested": sum(1 for c in cells if c.is_contested),
        "n_cells_directly_opposed": sum(1 for c in cells if c.is_directly_opposed),
        "n_cells_with_unknown_lab_records": sum(
            1 for c in cells if c.n_records_without_lab > 0
        ),
        "n_cells_matched_on_unknown_conditions": sum(
            1 for c in cells if c.matched_on_unknown_conditions
        ),
        "n_condition_slots": N_CONDITION_SLOTS,
        "cells_by_min_condition_specificity": _by_specificity(cells),
        "replicated_across_labs_by_min_condition_specificity": _by_specificity(
            [c for c in cells if c.replicated_across_labs]
        ),
        "replicated_across_studies_by_min_condition_specificity": _by_specificity(
            [c for c in cells if c.replicated_across_studies]
        ),
        "agreement_defined_by_min_condition_specificity": _by_specificity(defined),
        "mean_agreement_by_min_condition_specificity": _mean_agreement_by_specificity(
            defined
        ),
        "mean_agreement_where_defined": (
            sum(c.agreement_studies.value for c in defined) / len(defined)
            if defined
            else None
        ),
        "modal_direction_counts": _modal_direction_counts(cells),
        "confidence_level": CONFIDENCE_LEVEL,
        "interval_method": "wilson_score",
    }


def _by_specificity(cells: Sequence[ReplicationCell]) -> Dict[str, int]:
    """Cell counts keyed by minimum condition specificity, every slot present.

    Zeros are materialised for the same reason the coverage matrix keeps its
    empty cells: "no replicated claim specifies all four conditions" is the
    finding, and it vanishes if the key is simply absent.
    """
    counts = {str(k): 0 for k in range(N_CONDITION_SLOTS + 1)}
    for cell in cells:
        counts[str(cell.min_condition_specificity)] += 1
    return counts


def _mean_agreement_by_specificity(
    defined: Sequence[ReplicationCell],
) -> Dict[str, Optional[float]]:
    buckets: Dict[str, List[float]] = {
        str(k): [] for k in range(N_CONDITION_SLOTS + 1)
    }
    for cell in defined:
        buckets[str(cell.min_condition_specificity)].append(
            cell.agreement_studies.value
        )
    return {
        # None, never 0.0: no cell at this specificity is not "no agreement".
        key: (sum(values) / len(values) if values else None)
        for key, values in buckets.items()
    }


def _modal_direction_counts(cells: Sequence[ReplicationCell]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for cell in cells:
        key = cell.modal_direction or "undetermined"
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))
