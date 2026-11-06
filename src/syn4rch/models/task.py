"""The direction-of-effect prediction task, and the rules that define its data.

    given (population_id, valence, outcome_id, species, acuity,
           perturbation_class, perturbation_polarity, perturbation_target,
           sex, diet, challenge)
    predict direction in {increase, decrease, no_change}

This is deliberately a three-way *classification* of the sign of an effect, not
a quantitative simulation. ``docs/SCHEMA.md`` 6 refuses to represent
dose-response curves at all (the corpus does not report them consistently), and
refuses any aggregate "rejuvenation score". A sign is the strongest claim the
dataset can support, and whether even that is predictable is the open question
this module exists to answer.

What this module will not do
----------------------------
**It never silently drops a record.** Two of the five ``Direction`` values
cannot be turned into a label of this task:

* ``not_reported`` -- the paper measured something and did not say which way it
  went. Guessing would fabricate the observation.
* ``non_monotonic`` -- a real, recorded response shape with no single sign.
  ``docs/DECISIONS.md`` D-011 already refuses to treat it as contradicting a
  monotonic direction; collapsing it into one here would be the same error.

Both are **excluded and counted**, with the reason attached, in
:class:`ExclusionLedger`. A dataset whose exclusions are not reported is a
dataset whose label distribution cannot be interpreted: dropping
``not_reported`` rows quietly is how a benchmark comes to look better than the
literature it was built from.

**It never trains on Tier B.** ``docs/SCHEMA.md`` 6: "Tier B is an annotation
layer, not a training target." Tier B outcomes (lifespan, frailty, grip
strength) are therefore excluded by default and counted under their own reason,
not merged into some other bucket. ``exclude_tier_b=False`` exists only so a
caller can *measure* how much Tier B data there is; using it to train is a
schema violation, and the ledger records that the switch was flipped.

The v0.2.0 feature tuple
------------------------
Schema v0.2.0 gave every record a *signed manipulation axis*. Before it, 185 of
362 curated records carried ``valence: unclear`` -- 170 of them molecular
manipulations whose effect on firing is honestly unknown -- so half the corpus
had no signed feature at all. ``perturbation_class`` and
``perturbation_polarity`` supply one without coercing ``valence``, which stays
``unclear`` where that is the truth (SCHEMA.md 3.4), and which is *kept* as a
feature here so that the comparison against the previous benchmark is
like-for-like.

Two properties of the new features are enforced by this module rather than
left to a convention:

**A polarity never travels without its class.** ``decrease`` means less firing
under ``circuit_activity`` and less of a gene product under ``molecular``;
those are different physical interventions, and pooling them is the sign error
``docs/DECISIONS.md`` D-019 was written against. So:

* the three components come from :func:`syn4rch.qc.keys.perturbation_key`,
  which returns ``(class, polarity, target)`` as one tuple and is not
  reimplemented here;
* :class:`Features` raises :class:`PolarityWithoutClassError` at construction
  if a polarity is present with no class, so no such feature tuple can exist
  anywhere in the pipeline, not merely at dataset-build time;
* :meth:`Features.cell` raises the same error if a cell key names
  ``perturbation_polarity`` without ``perturbation_class``, so a model cannot
  key on the polarity alone. :data:`PERTURBATION_AXIS` is the composite name
  that expands to all three and is what the baselines use.

**A condition feature distinguishes "unknown" from a value.** ``sex``, ``diet``
and ``challenge`` come from :func:`syn4rch.qc.keys.condition_key`, which pools
``None`` with ``not_reported`` into the single token
:data:`UNKNOWN_CONDITION` and documents why (DECISIONS.md D-023). That
convention is reused, not reinvented. Because pooling hides how much was
actually known, :func:`syn4rch.qc.keys.condition_specificity` travels with every
example and its distribution is carried into
:meth:`TaskDataset.describe`: across the first 362 curated records, 153 know no
condition slot, 161 know one, 43 know two, 5 know three and none knows all
four. A feature set that did not say so would be a feature set that lies about
what the model had to work with.

Grouping
--------
Every grouping this task supports is declared in :data:`GROUPINGS` and the
group keys are computed here, at dataset-build time, from
``syn4rch.qc.keys`` -- the same functions QC uses, so that "same lab" and
"same population" mean one thing across the project. ``docs/SCHEMA.md`` 5
makes an unlabelled cross-validated number an invalid result, so the grouping
registry is the single place a grouping name becomes legitimate, and
``syn4rch.models.evaluate`` refuses any name that is not in it.

Python 3.10 compatible, standard library only.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from ..qc.keys import (
    CONDITION_FIELDS,
    NULL as _KEY_NULL,
    UNKNOWN as _KEY_UNKNOWN,
    condition_key,
    condition_specificity,
    lab_key,
    norm,
    perturbation_key,
    population_key,
    study_key,
)
from ..schema.enums import Direction, OutcomeTier

__all__ = [
    "LABELS",
    "LABELLABLE_DIRECTIONS",
    "UNLABELLABLE_DIRECTIONS",
    "FEATURE_FIELDS",
    "PERTURBATION_AXIS",
    "PERTURBATION_AXIS_FIELDS",
    "CONDITION_FEATURE_FIELDS",
    "DERIVED_CELL_FIELDS",
    "NULL_FEATURE",
    "UNKNOWN_CONDITION",
    "ACCEPTED_SCHEMA_VERSIONS",
    "PolarityWithoutClassError",
    "expand_cell_fields",
    "ExclusionReason",
    "EXCLUSION_REASONS",
    "ExclusionLedger",
    "Features",
    "Example",
    "Grouping",
    "GROUPINGS",
    "DEFAULT_GROUPINGS",
    "DirectionTask",
    "TASK",
    "TaskDataset",
    "build_dataset",
    "load_dataset",
    "default_records_path",
    "repo_root",
]

#: The three labels of this task, in a fixed order. Every probability vector,
#: confusion matrix and per-class table in this package uses this order, so
#: that two reports can be compared column by column.
LABELS: Tuple[str, ...] = ("increase", "decrease", "no_change")

#: ``Direction`` values that map onto :data:`LABELS`.
LABELLABLE_DIRECTIONS: Tuple[Direction, ...] = (
    Direction.increase,
    Direction.decrease,
    Direction.no_change,
)

#: ``Direction`` values that are valid records but not labels of this task.
UNLABELLABLE_DIRECTIONS: Tuple[Direction, ...] = (
    Direction.not_reported,
    Direction.non_monotonic,
)

#: The feature tuple, in the order the task statement gives it. The first five
#: are the schema v0.1.0 features and are kept unchanged -- including
#: ``valence``, which is still a separately meaningful field and whose removal
#: would make the comparison against the previous benchmark incomparable. The
#: last six arrived with schema v0.2.0.
FEATURE_FIELDS: Tuple[str, ...] = (
    "population_id",
    "valence",
    "outcome_id",
    "species",
    "acuity",
    "perturbation_class",
    "perturbation_polarity",
    "perturbation_target",
    "sex",
    "diet",
    "challenge",
)

#: The five features this task had under schema v0.1.0, kept as data so that
#: the before/after comparison in the benchmark report is a statement about
#: declared feature sets rather than about which rows someone remembered.
FEATURE_FIELDS_V1: Tuple[str, ...] = (
    "population_id",
    "valence",
    "outcome_id",
    "species",
    "acuity",
)

#: The composite feature name. It expands to :data:`PERTURBATION_AXIS_FIELDS`,
#: and it exists so that the ordinary way to use the perturbation in a cell key
#: is the way that carries the class with the polarity. See
#: :class:`syn4rch.schema.enums.PerturbationPolarity`: a polarity read without
#: its class is meaningless, because ``decrease`` is less firing under
#: ``circuit_activity`` and less of a gene product under ``molecular``.
PERTURBATION_AXIS = "perturbation_axis"

#: The three components of :func:`syn4rch.qc.keys.perturbation_key`, in order.
PERTURBATION_AXIS_FIELDS: Tuple[str, ...] = (
    "perturbation_class",
    "perturbation_polarity",
    "perturbation_target",
)

#: The condition features of this task. ``age_group`` is deliberately not one:
#: it is a *derived* field (``rules.assign_age_group`` returns ``not_reported``
#: for every non-mouse by design), and the task statement names three condition
#: fields, not four. It is still part of
#: :func:`syn4rch.qc.keys.condition_key`, so ``condition_specificity`` counts
#: four slots while a model may use three; both counts are reported by
#: :meth:`TaskDataset.describe` rather than one being quietly substituted for
#: the other.
CONDITION_FEATURE_FIELDS: Tuple[str, ...] = ("sex", "diet", "challenge")

#: Names a cell key may use that are not themselves features of the task.
DERIVED_CELL_FIELDS: Tuple[str, ...] = ("population_key",)

#: Stand-in for an absent feature value inside a composite cell key. It is a
#: visible token rather than an empty string so that a cell key printed in a
#: report cannot be mistaken for a populated one. Taken from
#: ``syn4rch.qc.keys`` rather than redefined, so that a cell key and a QC key
#: print an absent component the same way.
NULL_FEATURE = _KEY_NULL

#: Token for a condition slot that is not known. ``None`` (the curator did not
#: address it) and ``not_reported`` (the paper was checked and does not say) are
#: different statements about the curation and the same statement about the
#: experiment, so :func:`syn4rch.qc.keys.condition_key` pools them, and this
#: module reuses that single convention rather than inventing a second one. It
#: is a value, distinguishable from every real condition value and from
#: :data:`NULL_FEATURE`; a feature set in which "unknown" was indistinguishable
#: from a stated value would let a model learn from an absence.
UNKNOWN_CONDITION = _KEY_UNKNOWN


class PolarityWithoutClassError(ValueError):
    """A perturbation polarity was about to be used without its class.

    Raised by :class:`Features` at construction and by :meth:`Features.cell`
    when a cell key names ``perturbation_polarity`` and not
    ``perturbation_class``. The pairing is structural here for the same reason
    it is structural in :func:`syn4rch.qc.keys.perturbation_key`: a model that
    keyed on the polarity alone would pool "less firing" with "less of a gene
    product" and reintroduce the sign error DECISIONS.md D-019 was written
    against.
    """


def expand_cell_fields(*fields_: str) -> Tuple[str, ...]:
    """Expand :data:`PERTURBATION_AXIS` into its components, in order.

    Used both by :meth:`Features.cell` and by the baselines' ``describe()``, so
    that the chain a report prints is the chain that was actually keyed on.
    """
    out: List[str] = []
    for name in fields_:
        if name == PERTURBATION_AXIS:
            out.extend(PERTURBATION_AXIS_FIELDS)
        else:
            out.append(name)
    return tuple(out)


# --------------------------------------------------------------------------- #
# exclusions
# --------------------------------------------------------------------------- #

class ExclusionReason:
    """Reason codes for records that cannot become examples of this task.

    Plain string constants rather than an ``Enum``: these codes are written
    into reports and JSON, and are compared against strings in tests.
    """

    direction_absent = "direction_absent"
    direction_not_reported = "direction_not_reported"
    direction_non_monotonic = "direction_non_monotonic"
    outcome_id_absent = "outcome_id_absent"
    outcome_tier_b = "outcome_tier_b"
    population_absent = "population_absent"
    valence_absent = "valence_absent"


#: Reason code -> (one-line explanation, the clause that justifies excluding it).
EXCLUSION_REASONS: Dict[str, Tuple[str, str]] = {
    ExclusionReason.direction_absent: (
        "direction is null: the record states no effect at all",
        "SCHEMA.md 2 rule 4 (a field not reported is null) + 3.5",
    ),
    ExclusionReason.direction_not_reported: (
        "direction is 'not_reported': measured, but the paper does not say "
        "which way it went, so no label exists to predict",
        "SCHEMA.md 3.5",
    ),
    ExclusionReason.direction_non_monotonic: (
        "direction is 'non_monotonic': a recorded response shape with no "
        "single sign; collapsing it to one would be an interpretation",
        "SCHEMA.md 3.5 + DECISIONS.md D-011",
    ),
    ExclusionReason.outcome_id_absent: (
        "outcome_id is null, so the readout cannot be identified and the "
        "feature tuple is incomplete",
        "SCHEMA.md 3.5 (outcome_id required)",
    ),
    ExclusionReason.outcome_tier_b: (
        "outcome is Tier B (aging/lifespan/healthspan), which this project "
        "refuses to use as a training target",
        "SCHEMA.md 6",
    ),
    ExclusionReason.population_absent: (
        "neither population_id nor population_raw is present, so the "
        "population cannot be identified even by the paper's own wording",
        "SCHEMA.md 3.2 + DECISIONS.md D-010",
    ),
    ExclusionReason.valence_absent: (
        "valence is null; note that valence 'unclear' is NOT excluded, it is "
        "a legitimate feature value (SCHEMA.md 3.4 forbids coercing it)",
        "SCHEMA.md 3.4",
    ),
}

#: Order in which reasons are tested. The first match is the record's *primary*
#: reason (what the headline ledger counts); every reason that applied is also
#: counted separately, so a record excluded for two reasons shows up in both
#: tallies and neither is hidden behind the other.
_REASON_ORDER: Tuple[str, ...] = (
    ExclusionReason.direction_absent,
    ExclusionReason.direction_not_reported,
    ExclusionReason.direction_non_monotonic,
    ExclusionReason.outcome_tier_b,
    ExclusionReason.outcome_id_absent,
    ExclusionReason.population_absent,
    ExclusionReason.valence_absent,
)


@dataclass
class ExclusionLedger:
    """Counted, reasoned account of every input record that is not an example.

    ``primary`` partitions the excluded records (each counted exactly once).
    ``co_occurring`` counts every reason that applied to any record, so it can
    exceed ``n_excluded``. ``record_ids`` keeps the identifiers, capped, so a
    curator can go and look.
    """

    n_input: int = 0
    n_examples: int = 0
    primary: Dict[str, int] = field(default_factory=dict)
    co_occurring: Dict[str, int] = field(default_factory=dict)
    record_ids: Dict[str, List[str]] = field(default_factory=dict)
    max_ids_per_reason: int = 25
    tier_b_excluded: bool = True

    @property
    def n_excluded(self) -> int:
        return sum(self.primary.values())

    def add(self, record_id: Optional[str], reasons: Sequence[str]) -> None:
        if not reasons:
            raise ValueError("ExclusionLedger.add called with no reason")
        unknown = [r for r in reasons if r not in EXCLUSION_REASONS]
        if unknown:
            raise ValueError("undeclared exclusion reason(s): %r" % (unknown,))
        primary = reasons[0]
        self.primary[primary] = self.primary.get(primary, 0) + 1
        for reason in reasons:
            self.co_occurring[reason] = self.co_occurring.get(reason, 0) + 1
        bucket = self.record_ids.setdefault(primary, [])
        if record_id and len(bucket) < self.max_ids_per_reason:
            bucket.append(record_id)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_input_records": self.n_input,
            "n_examples": self.n_examples,
            "n_excluded": self.n_excluded,
            "excluded_fraction": (
                None if self.n_input == 0 else self.n_excluded / self.n_input
            ),
            "tier_b_excluded": self.tier_b_excluded,
            "primary_reason_counts": dict(sorted(self.primary.items())),
            "co_occurring_reason_counts": dict(sorted(self.co_occurring.items())),
            "reason_explanations": {
                code: {"why": why, "clause": clause}
                for code, (why, clause) in EXCLUSION_REASONS.items()
            },
            "example_record_ids": {
                code: list(ids) for code, ids in sorted(self.record_ids.items())
            },
        }


# --------------------------------------------------------------------------- #
# features, labels, examples
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Features:
    """The feature tuple of one example.

    ``population_id`` is the *resolved* vocabulary id and is ``None`` when the
    paper's population did not resolve (SCHEMA.md 3.2 forbids force-mapping
    it). ``population_key`` is the grouping key from ``syn4rch.qc.keys``, which
    falls back to ``"raw:" + normalised population_raw`` (DECISIONS.md D-010)
    so two unresolved populations never collide into one cell.

    The two are kept apart on purpose: the frequency baselines key on
    ``population_key`` (an unresolved population is still a population), while
    :class:`~syn4rch.models.baselines.NeurotransmitterSign` needs
    ``population_id`` because a transmitter can only be read out of the
    vocabulary, and must abstain when there is no id to look up.

    The three perturbation components are :func:`syn4rch.qc.keys.perturbation_key`
    unpacked, with its ``<null>`` token turned back into ``None``. ``None``
    means the component is absent; it never means "no effect". The three
    condition fields are :func:`syn4rch.qc.keys.condition_key`'s tokens and are
    therefore never ``None``: an unknown condition reads
    :data:`UNKNOWN_CONDITION`.

    ``condition_specificity`` is carried on the feature tuple but is **not** a
    feature -- :meth:`cell` refuses it. It is here so that every report can say
    how much of the condition was actually known behind a number.
    """

    population_id: Optional[str]
    population_key: str
    valence: Optional[str]
    outcome_id: Optional[str]
    species: Optional[str]
    acuity: Optional[str]
    perturbation_class: Optional[str] = None
    perturbation_polarity: Optional[str] = None
    perturbation_target: Optional[str] = None
    sex: str = UNKNOWN_CONDITION
    diet: str = UNKNOWN_CONDITION
    challenge: str = UNKNOWN_CONDITION
    #: How many of ``condition_key``'s four slots are known (0..4). Reporting
    #: only; not a feature.
    condition_specificity: int = 0

    def __post_init__(self) -> None:
        if self.perturbation_polarity is not None and self.perturbation_class is None:
            raise PolarityWithoutClassError(
                "perturbation_polarity=%r with no perturbation_class: a "
                "polarity is the sign of the manipulation along its own "
                "class's axis and means different physical things in "
                "different classes, so this feature tuple cannot exist "
                "(see syn4rch.schema.enums.PerturbationPolarity, "
                "DECISIONS.md D-019)" % (self.perturbation_polarity,)
            )

    @property
    def perturbation_axis(self) -> Tuple[str, str, str]:
        """``(class, polarity, target)`` with absences shown as ``<null>``.

        The only accessor that returns the polarity, and it cannot return it
        without the class.
        """
        return tuple(  # type: ignore[return-value]
            NULL_FEATURE if value in (None, "") else str(value)
            for value in (
                self.perturbation_class,
                self.perturbation_polarity,
                self.perturbation_target,
            )
        )

    def cell(self, *fields_: str) -> Tuple[str, ...]:
        """A hashable cell key over the named fields, nulls made visible.

        :data:`PERTURBATION_AXIS` expands to its three components. Naming
        ``perturbation_polarity`` without ``perturbation_class`` raises
        :class:`PolarityWithoutClassError`: this is the check that makes
        "a model may not key on polarity alone" structural rather than a
        convention someone remembers.
        """
        names = expand_cell_fields(*fields_)
        if "perturbation_polarity" in names and "perturbation_class" not in names:
            raise PolarityWithoutClassError(
                "cell key %r names perturbation_polarity without "
                "perturbation_class. Use %r, which expands to %s and keeps "
                "them together."
                % (list(fields_), PERTURBATION_AXIS, list(PERTURBATION_AXIS_FIELDS))
            )
        out: List[str] = []
        for name in names:
            if name not in DERIVED_CELL_FIELDS and name not in FEATURE_FIELDS:
                raise KeyError("not a feature of this task: %r" % name)
            value = getattr(self, name)
            out.append(NULL_FEATURE if value in (None, "") else str(value))
        return tuple(out)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "population_id": self.population_id,
            "population_key": self.population_key,
            "valence": self.valence,
            "outcome_id": self.outcome_id,
            "species": self.species,
            "acuity": self.acuity,
            "perturbation_class": self.perturbation_class,
            "perturbation_polarity": self.perturbation_polarity,
            "perturbation_target": self.perturbation_target,
            "sex": self.sex,
            "diet": self.diet,
            "challenge": self.challenge,
            "condition_specificity": self.condition_specificity,
        }


@dataclass(frozen=True)
class Example:
    """One labelled example: features, label, and its group keys.

    The group keys are computed once, here, so that the evaluation harness
    never needs the record back and therefore cannot quietly group by something
    that is not declared in :data:`GROUPINGS`.
    """

    record_id: Optional[str]
    features: Features
    label: str
    group_keys: Mapping[str, Optional[str]]

    def group(self, grouping: str) -> Optional[str]:
        if grouping not in GROUPINGS:
            raise KeyError(
                "undeclared grouping %r; declared: %s"
                % (grouping, ", ".join(sorted(GROUPINGS)))
            )
        return self.group_keys.get(grouping)


def _unnull(value: str) -> Optional[str]:
    """``<null>`` back to ``None``; anything else unchanged."""
    return None if value in (NULL_FEATURE, "") else value


def extract_features(record) -> Features:
    """Feature tuple of a ``TierARecord``. Never invents a value.

    The perturbation components come from
    :func:`syn4rch.qc.keys.perturbation_key` and the condition from
    :func:`syn4rch.qc.keys.condition_key` -- the same functions QC uses, not a
    second implementation of them. Both fall back to
    :mod:`syn4rch.schema.rules` where the v0.2.0 fields are not stored, so a
    record written before the migration yields the same features as the same
    record after it.

    Raises :class:`PolarityWithoutClassError` (from :class:`Features`) if the
    record carries a polarity with no class.
    """
    cls, polarity, target = perturbation_key(record)
    # Read by name, not by unpacking condition_key positionally. This was
    # `sex, diet, _age_group, challenge = condition_key(record)` until D-060
    # added a fifth condition dimension, at which point it raised ValueError
    # and took the whole model layer with it. CONDITION_FEATURE_FIELDS is the
    # declaration of which conditions are features; this follows it.
    _condition = dict(zip(CONDITION_FIELDS, condition_key(record)))
    sex = _condition["sex"]
    diet = _condition["diet"]
    challenge = _condition["challenge"]
    return Features(
        population_id=(record.population_id or None),
        population_key=population_key(record),
        valence=(str(record.valence) if record.valence is not None else None),
        outcome_id=(record.outcome_id or None),
        species=(str(record.species) if record.species is not None else None),
        acuity=(str(record.acuity) if record.acuity is not None else None),
        perturbation_class=_unnull(cls),
        perturbation_polarity=_unnull(polarity),
        perturbation_target=_unnull(target),
        sex=sex,
        diet=diet,
        challenge=challenge,
        condition_specificity=condition_specificity(record),
    )


def extract_label(record) -> Optional[str]:
    """The task label, or ``None`` when the record carries no label.

    ``None`` means "this record cannot be labelled", which is a different
    statement from any of the three labels and must not be collapsed into one.
    """
    direction = record.direction
    if direction is None:
        return None
    if direction in LABELLABLE_DIRECTIONS:
        return str(direction)
    return None


def exclusion_reasons(record, exclude_tier_b: bool = True) -> List[str]:
    """Every reason this record cannot be an example, in :data:`_REASON_ORDER`.

    Empty list means the record is usable. ``valence == 'unclear'`` is usable:
    SCHEMA.md 3.4 calls it "a valid and expected answer" and forbids coercing
    it, so it travels as a feature value and the baselines may abstain on it.
    """
    found: List[str] = []
    direction = record.direction
    if direction is None:
        found.append(ExclusionReason.direction_absent)
    elif direction is Direction.not_reported:
        found.append(ExclusionReason.direction_not_reported)
    elif direction is Direction.non_monotonic:
        found.append(ExclusionReason.direction_non_monotonic)

    if exclude_tier_b and record.outcome_tier is OutcomeTier.B:
        found.append(ExclusionReason.outcome_tier_b)
    if not record.outcome_id:
        found.append(ExclusionReason.outcome_id_absent)
    if not record.population_id and not norm(record.population_raw):
        found.append(ExclusionReason.population_absent)
    if record.valence is None:
        found.append(ExclusionReason.valence_absent)

    order = {code: i for i, code in enumerate(_REASON_ORDER)}
    return sorted(set(found), key=lambda c: order[c])


# --------------------------------------------------------------------------- #
# groupings (SCHEMA.md 5)
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Grouping(object):
    """A declared way to group rows so that a split does not leak.

    ``leakage_prone`` marks a grouping that SCHEMA.md 5 warns against relying
    on. It is kept available because a number computed under it is a useful
    *upper bound* -- and because the honest way to show that row-level
    cross-validation flatters a model is to print it next to the grouped
    number, labelled, rather than to omit it.
    """

    name: str
    description: str
    schema_clause: str
    leakage_prone: bool = False
    may_be_ungroupable: bool = False
    #: The name of the equivalent generator in ``syn4rch.feasibility.splits``,
    #: or ``None`` when that module has no counterpart. Used only to *look*
    #: for a fold generator there; a generator found this way is still
    #: validated against the actual data before it is used.
    feasibility_name: Optional[str] = None


def _lab_group(record) -> Optional[str]:
    # The group *is* the lab key: the normalised last author (D-056). It was
    # "author | affiliation" until D-056, which made this grouping identical to
    # study_id on a corpus where no two papers spell one department the same.
    return lab_key(record)


def _population_outcome_group(record) -> Optional[str]:
    outcome = record.outcome_id
    if not outcome:
        return None
    return "%s x %s" % (population_key(record), norm(outcome))


#: Every grouping the evaluation harness will accept. A name absent from here
#: cannot be used, which is how SCHEMA.md 5's "naming the grouping" is made
#: structural rather than a convention someone remembers.
GROUPINGS: Dict[str, Grouping] = {
    "study_id": Grouping(
        name="study_id",
        description="no study appears in both train and test",
        schema_clause="SCHEMA.md 5 item 1",
        feasibility_name="by_study",
    ),
    "lab": Grouping(
        name="lab",
        description=(
            "no last author appears in both train and test (a lab is a "
            "person, D-056); "
            "records with no lab identity cannot be grouped and are reported, "
            "not quietly pooled into one anonymous lab"
        ),
        schema_clause="SCHEMA.md 5 item 2 + DECISIONS.md D-009",
        may_be_ungroupable=True,
        feasibility_name="by_lab",
    ),
    "population": Grouping(
        name="population",
        description="leave-one-population-out: generalisation to a new population",
        schema_clause="SCHEMA.md 5 item 4",
        feasibility_name="leave_one_population_out",
    ),
    "population_outcome": Grouping(
        name="population_outcome",
        description=(
            "leave-one-(population x outcome)-cell-out: generalisation to a "
            "claim never seen, which is the condition a hypothesis generator "
            "actually operates under"
        ),
        schema_clause="SCHEMA.md 5 item 4",
        may_be_ungroupable=True,
    ),
    "record": Grouping(
        name="record",
        description=(
            "each row is its own group, i.e. ordinary row-level "
            "cross-validation. Leakage-prone by construction: the same study "
            "and the same lab appear on both sides"
        ),
        schema_clause="SCHEMA.md 5 (the thing it forbids reporting unlabelled)",
        leakage_prone=True,
    ),
}

#: Groupings the CLI reports by default, in report order.
DEFAULT_GROUPINGS: Tuple[str, ...] = (
    "study_id",
    "lab",
    "population",
    "population_outcome",
    "record",
)

_GROUP_KEY_FNS = {
    "study_id": lambda r: (study_key(r) or None),
    "lab": _lab_group,
    "population": lambda r: (population_key(r) or None),
    "population_outcome": _population_outcome_group,
    "record": lambda r: (r.record_id or None),
}


def group_keys_for(record) -> Dict[str, Optional[str]]:
    """All declared group keys for one record. ``None`` = cannot be grouped."""
    return {name: _GROUP_KEY_FNS[name](record) for name in GROUPINGS}


# --------------------------------------------------------------------------- #
# the task object
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class DirectionTask:
    """The task definition, as data.

    Exists so that a report, a saved prediction file and a test all quote the
    same statement of what was being predicted, rather than three paraphrases.
    """

    name: str = "direction_of_effect"
    #: Bumped from 0.1.0 when the schema v0.2.0 features were added. A report
    #: and a saved prediction file that disagree about the feature tuple are
    #: two reports about two different tasks, so the version and the statement
    #: below are the single place either is written down.
    version: str = "0.2.0"
    statement: str = (
        "given (population_id, valence, outcome_id, species, acuity, "
        "perturbation_class, perturbation_polarity, perturbation_target, "
        "sex, diet, challenge), "
        "predict direction in {increase, decrease, no_change}"
    )
    labels: Tuple[str, ...] = LABELS
    feature_fields: Tuple[str, ...] = FEATURE_FIELDS
    #: The feature tuple this task had at version 0.1.0, named here so that a
    #: before/after table cannot be assembled from a recollection of it.
    previous_feature_fields: Tuple[str, ...] = FEATURE_FIELDS_V1
    previous_version: str = "0.1.0"
    is_classification: bool = True
    refuses: Tuple[str, ...] = (
        "quantitative effect size prediction (SCHEMA.md 6)",
        "any aggregate 'rejuvenation score' (SCHEMA.md 6)",
        "Tier B outcomes as a training target (SCHEMA.md 6)",
        "a cross-validated number without its grouping (SCHEMA.md 5)",
        "a feature keyed on perturbation_polarity without perturbation_class "
        "(SCHEMA.md 3.4b, DECISIONS.md D-019)",
        "a condition feature in which 'unknown' is indistinguishable from a "
        "stated value (DECISIONS.md D-023)",
    )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "task": self.name,
            "task_version": self.version,
            "statement": self.statement,
            "labels": list(self.labels),
            "feature_fields": list(self.feature_fields),
            "previous_version": self.previous_version,
            "previous_feature_fields": list(self.previous_feature_fields),
            "features_added_in_this_version": [
                name
                for name in self.feature_fields
                if name not in self.previous_feature_fields
            ],
            "refuses": list(self.refuses),
        }


#: The singleton task definition.
TASK = DirectionTask()


# --------------------------------------------------------------------------- #
# dataset
# --------------------------------------------------------------------------- #

@dataclass
class TaskDataset:
    """Labelled examples plus the counted account of what was left out."""

    examples: List[Example] = field(default_factory=list)
    #: The ``TierARecord`` each example came from, aligned index-for-index with
    #: :attr:`examples`. Kept so that a fold generator which takes *records*
    #: -- as ``syn4rch.feasibility.splits`` does -- can be handed the real
    #: records rather than a reconstruction of them, and so that a fold it
    #: returns can be mapped back to example indices without guesswork.
    records: List[Any] = field(default_factory=list)
    ledger: ExclusionLedger = field(default_factory=ExclusionLedger)
    source: Optional[str] = None
    load_errors: List[str] = field(default_factory=list)
    task: DirectionTask = TASK
    #: ``schema_version`` -> how many input records declared it. Carried
    #: because the v0.2.0 fields may be *stored* or *derived* depending on
    #: whether the file has been migrated, and a number computed on a
    #: half-migrated file must not be mistaken for one computed on a migrated
    #: one.
    schema_versions: Dict[str, int] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.examples)

    @property
    def is_empty(self) -> bool:
        return not self.examples

    def labels(self) -> List[str]:
        return [e.label for e in self.examples]

    def features(self) -> List[Features]:
        return [e.features for e in self.examples]

    def label_counts(self) -> Dict[str, int]:
        counts = {label: 0 for label in LABELS}
        for example in self.examples:
            counts[example.label] += 1
        return counts

    def perturbation_axis_summary(self) -> Dict[str, Any]:
        """What the signed manipulation axis actually looks like in this data.

        The point of schema v0.2.0 was that every record should have a signed
        axis where only 177 of 362 did. Whether it does is a measurement, so it
        is counted here rather than asserted: ``n_with_signed_polarity`` counts
        examples whose polarity is one of ``increase``/``decrease``/``remove``,
        and ``n_with_known_target`` counts those whose molecular target is
        identified at all -- without which a ``molecular`` cell pools whatever
        genes happen to share a population and a readout.
        """
        signed = {"increase", "decrease", "remove"}
        class_counts: Dict[str, int] = {}
        axis_counts: Dict[str, int] = {}
        n_signed = 0
        n_target = 0
        for example in self.examples:
            features = example.features
            cls = features.perturbation_class or NULL_FEATURE
            polarity = features.perturbation_polarity or NULL_FEATURE
            class_counts[cls] = class_counts.get(cls, 0) + 1
            pair = "%s / %s" % (cls, polarity)
            axis_counts[pair] = axis_counts.get(pair, 0) + 1
            if features.perturbation_polarity in signed:
                n_signed += 1
            if features.perturbation_target:
                n_target += 1
        return {
            "n_examples": len(self.examples),
            "class_counts": dict(sorted(class_counts.items())),
            "class_polarity_counts": dict(sorted(axis_counts.items())),
            "n_with_signed_polarity": n_signed,
            "n_with_known_target": n_target,
            "polarity_is_never_a_feature_alone": True,
            "note": (
                "a polarity is read with its class or not at all; "
                "Features.cell refuses a key that names the polarity without "
                "the class (SCHEMA.md 3.4b, DECISIONS.md D-019)"
            ),
        }

    def condition_summary(self) -> Dict[str, Any]:
        """How much of the condition was known, over the four slots and the three.

        ``by_slots_known`` is ``syn4rch.qc.keys.condition_specificity``: how
        many of ``(sex, diet, age_group, challenge)`` a record states. A
        conflict or a cell match between two records of specificity 0 matched on
        "both unknown", not on "both the same", and a report that does not
        distinguish the two is overclaiming.

        ``by_feature_slots_known`` is the same count over the three slots this
        task actually uses as features, so that the number next to a model's
        accuracy describes what the model had, not what QC has.
        """
        four: Dict[int, int] = {}
        three: Dict[int, int] = {}
        for example in self.examples:
            features = example.features
            four[features.condition_specificity] = (
                four.get(features.condition_specificity, 0) + 1
            )
            known = sum(
                1
                for name in CONDITION_FEATURE_FIELDS
                if getattr(features, name) != UNKNOWN_CONDITION
            )
            three[known] = three.get(known, 0) + 1
        return {
            "slots": list(CONDITION_FIELDS),
            "feature_slots": list(CONDITION_FEATURE_FIELDS),
            "unknown_token": UNKNOWN_CONDITION,
            "by_slots_known": {str(k): v for k, v in sorted(four.items())},
            "by_feature_slots_known": {str(k): v for k, v in sorted(three.items())},
            "note": (
                "None and not_reported are pooled into one 'unknown' token "
                "(DECISIONS.md D-023), which is why the specificity counts are "
                "reported: 0 slots known means a cell matched on 'both "
                "unknown', not on 'both the same'. age_group is in the "
                "condition key and is not a feature of this task, so both "
                "counts are given."
            ),
        }

    def group_values(self, grouping: str) -> List[Optional[str]]:
        """Group key per example, aligned with :attr:`examples`."""
        if grouping not in GROUPINGS:
            raise KeyError(
                "undeclared grouping %r; declared: %s"
                % (grouping, ", ".join(sorted(GROUPINGS)))
            )
        return [e.group_keys.get(grouping) for e in self.examples]

    def grouping_coverage(self, grouping: str) -> Dict[str, Any]:
        """How many rows and distinct groups a grouping can actually cover."""
        values = self.group_values(grouping)
        groupable = [v for v in values if v is not None]
        sizes: Dict[str, int] = {}
        for value in groupable:
            sizes[value] = sizes.get(value, 0) + 1
        return {
            "grouping": grouping,
            "n_rows_total": len(values),
            "n_rows_groupable": len(groupable),
            "n_rows_ungroupable": len(values) - len(groupable),
            "n_groups": len(sizes),
            "largest_group_rows": max(sizes.values()) if sizes else 0,
            "singleton_groups": sum(1 for n in sizes.values() if n == 1),
            "leakage_prone": GROUPINGS[grouping].leakage_prone,
        }

    def subset(self, indices: Iterable[int]) -> List[Example]:
        return [self.examples[i] for i in indices]

    def record_index(self) -> Dict[str, int]:
        """``record_id -> example index``.

        A record id that occurs twice among the examples would make this map
        ambiguous, so it raises: QC detects exact duplicates
        (``syn4rch.qc.duplicates``) and a modelling layer that silently kept
        one of two identical rows would hide that finding.
        """
        out: Dict[str, int] = {}
        for index, example in enumerate(self.examples):
            rid = example.record_id
            if not rid:
                continue
            if rid in out:
                raise ValueError(
                    "duplicate record_id %r among the examples; run "
                    "syn4rch.qc first" % rid
                )
            out[rid] = index
        return out

    def describe(self) -> Dict[str, Any]:
        return {
            "task": self.task.to_dict(),
            "source": self.source,
            "schema_versions_read": dict(sorted(self.schema_versions.items())),
            "n_examples": len(self.examples),
            "label_counts": self.label_counts(),
            "exclusions": self.ledger.to_dict(),
            "perturbation_axis": self.perturbation_axis_summary(),
            "conditions": self.condition_summary(),
            "groupings": {
                name: self.grouping_coverage(name) for name in GROUPINGS
            },
            "load_errors": list(self.load_errors),
        }


def build_dataset(
    records: Iterable[Any],
    *,
    exclude_tier_b: bool = True,
    source: Optional[str] = None,
) -> TaskDataset:
    """Turn ``TierARecord`` objects into a :class:`TaskDataset`.

    Every input record either becomes an example or appears in the ledger with
    at least one declared reason. There is no third outcome; that invariant is
    asserted below rather than documented, because a counted exclusion that
    silently becomes an uncounted one is exactly the failure this ledger
    exists to prevent.
    """
    ledger = ExclusionLedger(tier_b_excluded=exclude_tier_b)
    examples: List[Example] = []
    kept: List[Any] = []
    versions: Dict[str, int] = {}
    for record in records:
        ledger.n_input += 1
        declared = str(getattr(record, "schema_version", None))
        versions[declared] = versions.get(declared, 0) + 1
        reasons = exclusion_reasons(record, exclude_tier_b=exclude_tier_b)
        label = extract_label(record)
        if reasons or label is None:
            if not reasons:
                # Defensive: a label of None with no reason would be a silent
                # drop. Attribute it rather than lose it.
                reasons = [ExclusionReason.direction_absent]
            ledger.add(getattr(record, "record_id", None), reasons)
            continue
        kept.append(record)
        examples.append(
            Example(
                record_id=getattr(record, "record_id", None),
                features=extract_features(record),
                label=label,
                group_keys=group_keys_for(record),
            )
        )
    ledger.n_examples = len(examples)
    if ledger.n_examples + ledger.n_excluded != ledger.n_input:
        raise AssertionError(
            "exclusion ledger does not balance: %d examples + %d excluded "
            "!= %d input records" % (
                ledger.n_examples, ledger.n_excluded, ledger.n_input
            )
        )
    return TaskDataset(
        examples=examples,
        records=kept,
        ledger=ledger,
        source=source,
        schema_versions=versions,
    )


# --------------------------------------------------------------------------- #
# locating and loading the curated table
# --------------------------------------------------------------------------- #

#: Versions older than the current one whose **feature tuple** is identical to
#: it, and which this task may therefore read beside it.
#:
#: This list is not "versions we tolerate". A version belongs here only when
#: every field in :data:`TaskSpec.statement` means the same thing in it, so
#: that two records from two versions are two observations of one task rather
#: than two tasks silently pooled.
#:
#: * ``0.1.0`` -- the pooled claim key, kept so the v0.2.0 key change stays
#:   falsifiable (D-023).
#: * ``0.2.0`` -- differs from 0.3.0 only in condition *provenance*
#:   (``condition_evidence`` replacing a null scalar pair, D-037). Provenance is
#:   not a feature of this task and not a component of any key, so a 0.2.0
#:   record and a 0.3.0 record are the same observation. Listing it is what
#:   lets the conditioning verdict of D-036 be re-measured against the filled
#:   corpus on one task rather than two.
_FEATURE_EQUIVALENT_VERSIONS: Tuple[str, ...] = ("0.3.0", "0.2.0", "0.1.0")


def _accepted_schema_versions() -> Tuple[str, ...]:
    """The schema versions this task will read, current first.

    ``syn4rch.schema.record.SCHEMA_VERSION`` is read at import time rather than
    hard-coded, so this list follows the schema instead of drifting behind it.
    The rest come from :data:`_FEATURE_EQUIVALENT_VERSIONS`, and belong there on
    the argument recorded with it -- not because they are old, but because the
    task's features mean the same thing in them.
    """
    try:
        from ..schema.record import SCHEMA_VERSION
    except Exception:  # pragma: no cover - schema is stdlib-only
        return _FEATURE_EQUIVALENT_VERSIONS
    return tuple(
        dict.fromkeys((SCHEMA_VERSION,) + _FEATURE_EQUIVALENT_VERSIONS)
    )


#: Schema versions :func:`load_dataset` will read. See its docstring for why
#: reading two versions is admissible *here* and nowhere else.
ACCEPTED_SCHEMA_VERSIONS: Tuple[str, ...] = _accepted_schema_versions()


def repo_root() -> Optional[str]:
    """Directory holding ``data/`` and ``pyproject.toml``, or ``None``.

    Walks up from this file, as ``syn4rch.vocab.core.vocab_dir`` does. Returns
    ``None`` instead of guessing, so a caller in an unusual layout gets a clear
    "I do not know where the dataset is" rather than a path that does not exist.
    """
    import os

    current = os.path.abspath(os.path.dirname(__file__))
    for _ in range(8):
        if os.path.isdir(os.path.join(current, "data")) and os.path.isfile(
            os.path.join(current, "pyproject.toml")
        ):
            return current
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent
    return None


def default_records_path() -> Optional[str]:
    """``<repo>/data/curated/records.jsonl``, whether or not it exists."""
    import os

    root = repo_root()
    if root is None:
        return None
    return os.path.join(root, "data", "curated", "records.jsonl")


def load_dataset(
    path: Optional[str] = None,
    *,
    exclude_tier_b: bool = True,
    any_version: bool = False,
) -> TaskDataset:
    """Read the curated JSONL and build the dataset.

    A missing or empty file is **not** an error: it yields an empty dataset
    whose ``source`` records the path that was looked at. The caller (the CLI,
    the report) is then responsible for saying so. Reading one unparseable line
    stops the read and records the message in ``load_errors``; the records read
    before it are kept, because a partial dataset that says it is partial is
    more useful than an exception.

    Both schema versions in :data:`ACCEPTED_SCHEMA_VERSIONS` are read, and
    which ones were actually found is recorded in
    :attr:`TaskDataset.schema_versions` and printed by the report.
    ``syn4rch.schema.io`` refuses a cross-version load by default, and is right
    to: v0.1.0 field *semantics* must not be applied to v0.2.0 data. This layer
    is the one documented exception, for one reason that is checkable --
    :func:`syn4rch.qc.keys.perturbation_key` and
    :func:`syn4rch.qc.keys.condition_key` derive every v0.2.0 feature this task
    uses from ``method`` and ``valence`` when the stored field is absent, so a
    v0.1.0 record yields the same feature tuple as the same record after
    migration. A version outside the accepted set stops the read with a
    message rather than being coerced.
    """
    import os

    from ..schema.io import JsonlError, iter_jsonl

    if path is None:
        path = default_records_path()
    if path is None:
        dataset = build_dataset([], exclude_tier_b=exclude_tier_b, source=None)
        dataset.load_errors.append(
            "could not locate the repository root, so data/curated/records.jsonl "
            "could not be found; pass an explicit path"
        )
        return dataset

    if not os.path.isfile(path):
        dataset = build_dataset([], exclude_tier_b=exclude_tier_b, source=path)
        dataset.load_errors.append("no such file: %s" % path)
        return dataset

    records: List[Any] = []
    errors: List[str] = []
    try:
        for record in iter_jsonl(path, expect_version=None):
            declared = getattr(record, "schema_version", None)
            if not any_version and declared not in ACCEPTED_SCHEMA_VERSIONS:
                errors.append(
                    "refusing to read schema_version=%r; this task reads %s. "
                    "The read stopped here and the %d record(s) before it were "
                    "kept."
                    % (
                        declared,
                        " or ".join(repr(v) for v in ACCEPTED_SCHEMA_VERSIONS),
                        len(records),
                    )
                )
                break
            records.append(record)
    except JsonlError as exc:
        errors.append(str(exc))

    dataset = build_dataset(records, exclude_tier_b=exclude_tier_b, source=path)
    dataset.load_errors.extend(errors)
    return dataset
