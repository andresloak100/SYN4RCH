"""Eight baselines for direction-of-effect prediction.

These are the **bar to clear**. The point of writing them carefully is that a
sophisticated model is worth adopting only if it beats them under a named
grouping (``docs/SCHEMA.md`` 5), and the commonest way a modelling project
deceives itself is by comparing a learned model against a baseline that was
not given the same chance.

The eight, in increasing order of how much structure they assume:

======================================== ==========================================
``GlobalMajority``                       one label for everything; the floor
``OutcomeMajority``                      the modal direction of each readout
``PopulationOutcomeMajority``            the modal direction of each (population, outcome)
``PerturbationAxisCellMajority``         the axis cell against the floor, nothing between (v0.2.0)
``PerturbationAxisMajority``             + the signed perturbation axis, with backoff (v0.2.0)
``PerturbationAxisConditionMajority``    + the condition the work ran under (v0.2.0)
``NeurotransmitterSign``                 a declared sign rule, abstaining when unsure
``KGPathSign``                           sign propagation over a knowledge-graph edge table
======================================== ==========================================

The first three and ``NeurotransmitterSign`` key on the schema v0.1.0 feature
tuple and are **unchanged**: their cell keys, their backoff order and their tie
rule are exactly what produced the previous benchmark, so the before/after
comparison in the report is like-for-like and not a rerun of a different
measurement. :data:`FEATURE_GENERATION_V1` and
:data:`FEATURE_GENERATION_V2` tag which generation each baseline belongs to, as
data, so that table is assembled from a declaration rather than from a memory
of which row was which.

Two deliberate design choices run through all of them.

**Abstention is a first-class output.** ``predict`` returns ``None`` for "I
decline to answer". A baseline that abstains on a third of the test set and is
right on the rest has said something quite different from one that guesses on
everything, and merging the two into a single accuracy hides it. The evaluation
harness therefore reports coverage alongside accuracy and never scores an
abstention as a wrong answer or as a right one.

**Unavailability is reported, not simulated.** A baseline whose inputs do not
exist reports itself ``UNAVAILABLE`` with the reason. ``KGPathSign`` is in that
state in this repository today: there is no knowledge-graph edge table, and
inventing one so that the row in the table could be filled would be
fabrication of exactly the kind ``docs/SCHEMA.md`` 2 exists to prevent.

Python 3.10 compatible, standard library only (no import here can fail for
want of an optional wheel).
"""

import json
import os
from collections import deque
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .task import (
    CONDITION_FEATURE_FIELDS,
    FEATURE_FIELDS_V1,
    LABELS,
    PERTURBATION_AXIS,
    Example,
    Features,
    expand_cell_fields,
)

__all__ = [
    "Availability",
    "BaselineUnavailable",
    "Baseline",
    "FEATURE_GENERATION_V1",
    "FEATURE_GENERATION_V2",
    "BackoffMajority",
    "GlobalMajority",
    "OutcomeMajority",
    "PopulationOutcomeMajority",
    "PerturbationAxisCellMajority",
    "PerturbationAxisMajority",
    "PerturbationAxisConditionMajority",
    "NeurotransmitterSign",
    "KGPathSign",
    "POPULATION_OUTCOME_BACKOFF",
    "PERTURBATION_AXIS_CELL_BACKOFF",
    "PERTURBATION_AXIS_BACKOFF",
    "PERTURBATION_AXIS_CONDITION_BACKOFF",
    "NEUROTRANSMITTER_VALENCE_RULES",
    "UNSIGNED_TRANSMITTERS",
    "EDGE_TABLE_FORMAT",
    "EdgeTable",
    "all_baseline_factories",
    "baselines_by_feature_generation",
    "build_all",
]


# --------------------------------------------------------------------------- #
# availability
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Availability:
    """Whether a baseline can run at all, and if not, exactly why.

    ``probes`` lists each place that was looked at and what was found, so that
    "unavailable" is a reproducible statement about the repository rather than
    an opinion.
    """

    available: bool
    reason: Optional[str] = None
    probes: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "available": self.available,
            "reason": self.reason,
            "probes": list(self.probes),
        }


AVAILABLE = Availability(available=True)


class BaselineUnavailable(RuntimeError):
    """Raised by :meth:`Baseline.require_available` for callers that want a stop."""


# --------------------------------------------------------------------------- #
# which generation of the feature tuple a baseline uses
# --------------------------------------------------------------------------- #

#: A baseline that keys only on the schema v0.1.0 features. Its number in this
#: benchmark is directly comparable with the previous one.
FEATURE_GENERATION_V1 = "v0.1.0"

#: A baseline that keys on at least one schema v0.2.0 feature.
FEATURE_GENERATION_V2 = "v0.2.0"


# --------------------------------------------------------------------------- #
# base class
# --------------------------------------------------------------------------- #

class Baseline:
    """Common interface: ``fit`` / ``predict`` / ``predict_proba``.

    Subclasses override :meth:`_fit` and :meth:`_predict_one`. ``predict``
    returns one of :data:`~syn4rch.models.task.LABELS` or ``None`` (abstain).
    ``predict_proba`` returns ``None`` when no probability is meaningful for
    that baseline -- which is a statement, not a gap: a hard rule with no
    estimated frequency behind it has no calibrated confidence to report, and
    emitting 1.0 would be a lie the calibration plot would then confirm.
    """

    #: Short stable name used as a column heading and a JSON key.
    name = "baseline"
    #: One-line description for the report.
    description = ""
    #: Whether :meth:`predict_proba` can return anything.
    emits_proba = False
    #: Which generation of the task's feature tuple this baseline keys on;
    #: one of :data:`FEATURE_GENERATION_V1` / :data:`FEATURE_GENERATION_V2`.
    feature_generation = FEATURE_GENERATION_V1
    #: The features this baseline actually reads, declared rather than
    #: inferred. The report's before/after table is built from these two
    #: attributes, so a baseline that gained a feature and did not say so would
    #: show up in the wrong column -- which is why
    #: ``tests/test_models_baselines.py`` checks the declaration against the
    #: cell keys the baseline really uses.
    features_used: Tuple[str, ...] = ()

    def __init__(self) -> None:
        self.availability: Availability = AVAILABLE
        self.fitted: bool = False
        self.n_train: int = 0

    # -- availability ------------------------------------------------------

    @property
    def available(self) -> bool:
        return self.availability.available

    def require_available(self) -> None:
        if not self.available:
            raise BaselineUnavailable(
                "%s is unavailable: %s" % (self.name, self.availability.reason)
            )

    # -- fit ---------------------------------------------------------------

    def fit(self, examples: Sequence[Example]) -> "Baseline":
        """Fit on labelled examples. Always returns ``self``, never raises on
        an empty training set: a baseline fitted on nothing abstains on
        everything, which is the correct behaviour and is reportable."""
        self.n_train = len(examples)
        if self.available:
            self._fit(list(examples))
        self.fitted = True
        return self

    def _fit(self, examples: List[Example]) -> None:  # pragma: no cover - abstract
        raise NotImplementedError

    # -- predict -----------------------------------------------------------

    def predict(self, features: Iterable[Features]) -> List[Optional[str]]:
        if not self.available:
            return [None for _ in features]
        if not self.fitted:
            raise RuntimeError("%s.predict called before fit" % self.name)
        return [self._predict_one(f) for f in features]

    def _predict_one(self, features: Features) -> Optional[str]:  # pragma: no cover
        raise NotImplementedError

    def prediction_level(self, features: Features) -> Optional[str]:
        """Which rung of this baseline's structure answered, or ``None``.

        ``None`` means the question does not apply (the baseline has no
        structure to back off through). A baseline that does have one returns
        the level's name, and :attr:`primary_level` names the most specific
        rung. The evaluation harness uses the pair to report a metric
        restricted to the rows whose **primary cell was populated in
        training** -- which is the only subset on which a cell baseline and the
        global-majority floor can even differ, since on an unpopulated cell
        both fall back to the same global mode and their predictions are
        identical by construction. Pooling those rows into the comparison
        dilutes the difference with rows that carry no information about it.
        """
        return None

    @property
    def primary_level(self) -> Optional[str]:
        """Name of this baseline's most specific rung, or ``None``."""
        return None

    def predict_proba(
        self, features: Iterable[Features]
    ) -> Optional[List[Optional[Dict[str, float]]]]:
        """Per-example ``{label: probability}``, or ``None`` for the whole call
        when this baseline has no meaningful probability. An individual entry
        may be ``None`` where that one example is abstained on."""
        return None

    # -- reporting ---------------------------------------------------------

    def describe(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "availability": self.availability.to_dict(),
            "emits_proba": self.emits_proba,
            "feature_generation": self.feature_generation,
            "features_used": list(self.features_used),
            "n_train": self.n_train,
        }


# --------------------------------------------------------------------------- #
# frequency baselines
# --------------------------------------------------------------------------- #

def _count_labels(examples: Iterable[Example]) -> Dict[str, int]:
    counts = {label: 0 for label in LABELS}
    for example in examples:
        counts[example.label] = counts.get(example.label, 0) + 1
    return counts


def _modal(counts: Mapping[str, int]) -> Optional[str]:
    """Modal label of a count dict, or ``None`` if it is empty.

    Ties break on the order of :data:`~syn4rch.models.task.LABELS`. That is a
    documented, deterministic, *arbitrary* rule: a tie means the data does not
    choose, and a baseline that resolved it randomly would make the benchmark
    irreproducible, while one that resolved it by "whichever was seen first"
    would make it depend on file order.
    """
    best: Optional[str] = None
    best_count = -1
    for label in LABELS:
        count = counts.get(label, 0)
        if count > best_count:
            best, best_count = label, count
    if best_count <= 0:
        return None
    return best


def _proba_from_counts(counts: Mapping[str, int]) -> Optional[Dict[str, float]]:
    total = sum(counts.get(label, 0) for label in LABELS)
    if total <= 0:
        return None
    return {label: counts.get(label, 0) / total for label in LABELS}


class GlobalMajority(Baseline):
    """The modal direction of the training set, for every input.

    The floor. Any model that does not beat this has learned nothing, and on a
    corpus whose label distribution is as skewed as this one is expected to be,
    it is a higher floor than it looks: that is why accuracy is reported next
    to macro-F1, which this baseline scores badly on by construction (it
    scores 0 on two of the three classes).
    """

    name = "global_majority"
    description = "modal direction over the whole training set"
    emits_proba = True
    feature_generation = FEATURE_GENERATION_V1
    features_used = ()

    def __init__(self) -> None:
        super().__init__()
        self.counts: Dict[str, int] = {}
        self.prediction: Optional[str] = None

    def _fit(self, examples: List[Example]) -> None:
        self.counts = _count_labels(examples)
        self.prediction = _modal(self.counts)

    def _predict_one(self, features: Features) -> Optional[str]:
        return self.prediction

    def predict_proba(self, features):
        if not self.available or not self.fitted:
            return None
        proba = _proba_from_counts(self.counts)
        return [proba for _ in features]

    def describe(self) -> Dict[str, Any]:
        out = super().describe()
        out["train_label_counts"] = dict(self.counts)
        out["constant_prediction"] = self.prediction
        return out


class OutcomeMajority(Baseline):
    """The modal direction per ``outcome_id``, backing off to the global mode.

    Tests the hypothesis that direction is mostly a property of the readout --
    that food intake tends to go one way whatever you perturb. If this
    baseline is close to a learned model, the learned model has mostly
    rediscovered outcome marginals.
    """

    name = "outcome_majority"
    description = "modal direction per outcome_id, backing off to global"
    emits_proba = True
    feature_generation = FEATURE_GENERATION_V1
    features_used = ("outcome_id",)

    def __init__(self) -> None:
        super().__init__()
        self.cells: Dict[Tuple[str, ...], Dict[str, int]] = {}
        self.global_counts: Dict[str, int] = {}
        self.backoff_counts: Dict[str, int] = {}

    def _fit(self, examples: List[Example]) -> None:
        self.cells = {}
        for example in examples:
            key = example.features.cell("outcome_id")
            cell = self.cells.setdefault(key, {label: 0 for label in LABELS})
            cell[example.label] = cell.get(example.label, 0) + 1
        self.global_counts = _count_labels(examples)
        self.backoff_counts = {}

    def _counts_for(self, features: Features) -> Tuple[Dict[str, int], str]:
        cell = self.cells.get(features.cell("outcome_id"))
        if cell and sum(cell.values()) > 0:
            return cell, "outcome"
        return self.global_counts, "global"

    def _predict_one(self, features: Features) -> Optional[str]:
        counts, level = self._counts_for(features)
        self.backoff_counts[level] = self.backoff_counts.get(level, 0) + 1
        return _modal(counts)

    def predict_proba(self, features):
        if not self.available or not self.fitted:
            return None
        return [_proba_from_counts(self._counts_for(f)[0]) for f in features]

    def describe(self) -> Dict[str, Any]:
        out = super().describe()
        out["n_cells"] = len(self.cells)
        out["backoff_usage"] = dict(sorted(self.backoff_counts.items()))
        return out


#: The backoff chain of :class:`PopulationOutcomeMajority`, as data so that the
#: report can print it and a test can pin it.
POPULATION_OUTCOME_BACKOFF: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("population_outcome", ("population_key", "outcome_id")),
    ("outcome", ("outcome_id",)),
    ("population", ("population_key",)),
)


class BackoffMajority(Baseline):
    """Modal label of the first populated cell in a declared backoff chain.

    The shared machinery of every cell-frequency baseline here. A subclass
    declares :attr:`chain` -- an ordered sequence of ``(level name, cell
    fields)`` -- and nothing else about prediction. Keeping it in one place
    means the three cell baselines cannot drift apart in their tie rule, their
    "is this cell populated" test or their backoff accounting, and that a new
    cell definition is a line of data rather than a reimplementation whose
    numbers would not be comparable with the others'.

    Every level is tried in order and the first one seen in training wins,
    then the global mode, then abstention (reached only when the training set
    was empty). ``backoff_usage`` counts how often each rung was actually used:
    a baseline named for its primary cell that in fact spent its predictions
    two rungs down is not the model its name describes, and a candidate
    compared against it as though it were would be compared against the wrong
    bar.
    """

    #: ``((level name, cell fields), ...)``, most specific first. A field may
    #: be :data:`~syn4rch.models.task.PERTURBATION_AXIS`, which expands to the
    #: class, the polarity and the target together -- the only way this package
    #: permits a polarity into a cell key.
    chain: Tuple[Tuple[str, Tuple[str, ...]], ...] = ()

    def __init__(self) -> None:
        super().__init__()
        self.levels: Dict[str, Dict[Tuple[str, ...], Dict[str, int]]] = {}
        self.global_counts: Dict[str, int] = {}
        self.backoff_counts: Dict[str, int] = {}

    def _fit(self, examples: List[Example]) -> None:
        self.levels = {name: {} for name, _ in self.chain}
        for example in examples:
            for name, fields_ in self.chain:
                key = example.features.cell(*fields_)
                cell = self.levels[name].setdefault(
                    key, {label: 0 for label in LABELS}
                )
                cell[example.label] = cell.get(example.label, 0) + 1
        self.global_counts = _count_labels(examples)
        self.backoff_counts = {}

    def resolve(self, features: Features) -> Tuple[Dict[str, int], str]:
        """The counts this baseline would use, and the backoff level's name."""
        for name, fields_ in self.chain:
            cell = self.levels.get(name, {}).get(features.cell(*fields_))
            if cell and sum(cell.values()) > 0:
                return cell, name
        if sum(self.global_counts.values()) > 0:
            return self.global_counts, "global"
        return {}, "abstain"

    def _predict_one(self, features: Features) -> Optional[str]:
        counts, level = self.resolve(features)
        self.backoff_counts[level] = self.backoff_counts.get(level, 0) + 1
        return _modal(counts)

    def prediction_level(self, features: Features) -> Optional[str]:
        if not self.available or not self.fitted:
            return None
        return self.resolve(features)[1]

    @property
    def primary_level(self) -> Optional[str]:
        return self.chain[0][0] if self.chain else None

    def predict_proba(self, features):
        if not self.available or not self.fitted:
            return None
        return [_proba_from_counts(self.resolve(f)[0]) for f in features]

    def describe(self) -> Dict[str, Any]:
        out = super().describe()
        out["backoff_chain"] = [
            "(%s)" % ", ".join(expand_cell_fields(*fields_))
            for _, fields_ in self.chain
        ] + ["global", "abstain"]
        out["n_cells_per_level"] = {
            name: len(cells) for name, cells in sorted(self.levels.items())
        }
        out["backoff_usage"] = dict(sorted(self.backoff_counts.items()))
        return out


class PopulationOutcomeMajority(BackoffMajority):
    """Modal direction per (population, outcome) cell, with a declared backoff.

    Backoff chain, in order, each level used only when the level above it was
    never seen in training:

    1. ``(population_key, outcome_id)`` -- the cell itself.
    2. ``(outcome_id)`` -- the readout's marginal.
    3. ``(population_key)`` -- the population's marginal.
    4. the global mode.
    5. abstain (``None``) -- reached only when the training set was empty.

    **Why outcome before population.** ``outcome_id`` is required by
    ``docs/SCHEMA.md`` 3.5 and drawn from a 43-entry vocabulary, while
    ``population_id`` may legitimately be null (3.2) and unresolved
    populations each get their own ``raw:`` key (D-010). The outcome marginal
    is therefore both always defined and far better populated. The commonest
    way an unseen cell arises is a *new population* measured on a *familiar
    readout*, so the readout is the side more likely to carry usable
    information. The opposite order is defensible and is the obvious thing to
    vary; what is not defensible is leaving the order undocumented, because
    then the number in the table cannot be reproduced.

    **Valence is deliberately not in the cell key.** The class predicts per
    (population, outcome) as specified, which means it gives activation and
    inhibition of the same population the same answer. That is a real
    weakness and it is the point: :class:`NeurotransmitterSign` is the
    baseline that crosses valence in, so the pair brackets the question of
    whether direction is carried by the claim's identity or by the sign of the
    perturbation. ``backoff_usage`` reports how often the primary cell was
    actually available, which under a ``population_outcome`` grouping is by
    construction never -- and that number belongs in the report, because a
    baseline silently operating one rung down its own chain would otherwise
    look like a cell-level model.
    """

    name = "population_outcome_majority"
    description = (
        "modal direction per (population, outcome); backoff "
        "cell -> outcome -> population -> global -> abstain"
    )
    emits_proba = True
    feature_generation = FEATURE_GENERATION_V1
    features_used = ("population_key", "outcome_id")
    chain = POPULATION_OUTCOME_BACKOFF


# --------------------------------------------------------------------------- #
# the signed perturbation axis (schema v0.2.0)
# --------------------------------------------------------------------------- #

#: The backoff chain of :class:`PerturbationAxisMajority`. The first two rungs
#: are new in v0.2.0; the last three are :data:`POPULATION_OUTCOME_BACKOFF`
#: unchanged, so this baseline is a strict refinement of
#: :class:`PopulationOutcomeMajority` and the two differ only where the axis
#: cell was populated. ``backoff_usage`` is therefore the measurement of how
#: often the new axis had anything to say.
PERTURBATION_AXIS_BACKOFF: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("population_axis_outcome", ("population_key", PERTURBATION_AXIS, "outcome_id")),
    ("axis_outcome", (PERTURBATION_AXIS, "outcome_id")),
) + POPULATION_OUTCOME_BACKOFF

#: The backoff chain of :class:`PerturbationAxisConditionMajority`. Each
#: conditioned rung is followed immediately by its unconditioned twin, so the
#: chain never has to choose between "this exact condition" and "some
#: condition": it prefers the former and falls through to the latter. The pair
#: structure is what lets ``backoff_usage`` answer whether conditioning helped
#: or merely fragmented -- a conditioned rung that is almost never used is a
#: cell definition that split the data without separating anything.
PERTURBATION_AXIS_CONDITION_BACKOFF: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    (
        "population_axis_outcome_condition",
        ("population_key", PERTURBATION_AXIS, "outcome_id")
        + CONDITION_FEATURE_FIELDS,
    ),
    ("population_axis_outcome", ("population_key", PERTURBATION_AXIS, "outcome_id")),
    (
        "axis_outcome_condition",
        (PERTURBATION_AXIS, "outcome_id") + CONDITION_FEATURE_FIELDS,
    ),
    ("axis_outcome", (PERTURBATION_AXIS, "outcome_id")),
) + POPULATION_OUTCOME_BACKOFF


class PerturbationAxisMajority(BackoffMajority):
    """Modal direction per (population, perturbation axis, outcome) cell.

    The natural cell under schema v0.2.0, and the baseline the schema change
    was made for. Before v0.2.0, 185 of 362 records carried
    ``valence: unclear`` and so had no signed axis at all; a cell keyed on
    valence put all of them in one undifferentiated bucket. Now every record
    has ``perturbation_class`` and ``perturbation_polarity``, so a cell can
    distinguish "activation of this population" from "inhibition of it" and
    "loss of a gene product in it" from "gain of one".

    **The class travels with the polarity by construction.** The cell field is
    :data:`~syn4rch.models.task.PERTURBATION_AXIS`, which expands to
    ``(perturbation_class, perturbation_polarity, perturbation_target)``; the
    polarity cannot be named on its own, because
    :meth:`syn4rch.models.task.Features.cell` raises if it is. ``decrease``
    under ``circuit_activity`` is less firing and under ``molecular`` is less
    of a gene product; a cell that pooled them would be the sign error
    DECISIONS.md D-019 was written against, wearing a v0.2.0 badge.

    **The target is part of the cell, and usually absent.** A ``molecular``
    cell without it pools every gene that happens to share a population and a
    readout -- in this corpus, different genes from different laboratories. The
    cell key therefore includes ``perturbation_target`` so that, where the
    target *is* recorded, two different genes do not occupy one prediction
    cell; and the dataset summary reports how many examples have one, so that a
    number computed where none does cannot be read as though the cell were
    gene-specific.
    """

    name = "perturbation_axis_majority"
    description = (
        "modal direction per (population, perturbation_class, "
        "perturbation_polarity, perturbation_target, outcome); backoff "
        "cell -> (axis, outcome) -> (population, outcome) -> outcome -> "
        "population -> global -> abstain"
    )
    emits_proba = True
    feature_generation = FEATURE_GENERATION_V2
    features_used = (
        "population_key",
        "perturbation_class",
        "perturbation_polarity",
        "perturbation_target",
        "outcome_id",
    )
    chain = PERTURBATION_AXIS_BACKOFF


#: The chain of :class:`PerturbationAxisCellMajority`: one rung, then the
#: global mode. Deliberately *not* :data:`PERTURBATION_AXIS_BACKOFF`.
PERTURBATION_AXIS_CELL_BACKOFF: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("population_axis_outcome", ("population_key", PERTURBATION_AXIS, "outcome_id")),
)


class PerturbationAxisCellMajority(BackoffMajority):
    """The axis cell against the floor, with nothing else in between.

    The *measuring instrument*, as opposed to the best predictor. Its only
    difference from :class:`GlobalMajority` is the cell: where the
    ``(population, perturbation axis, outcome)`` cell was populated in
    training it predicts that cell's mode, and where it was not it predicts the
    global mode, exactly as ``GlobalMajority`` would. A paired comparison
    between the two therefore attributes every discordant row to the axis cell
    and to nothing else.

    :class:`PerturbationAxisMajority` is the better predictor, because the
    outcome and population marginals are real information and it uses them. But
    that makes its margin over ``GlobalMajority`` a mixture of the axis and
    those marginals, and the question this benchmark exists to answer is what
    the *axis* bought. Hence two baselines with the same primary cell: one to
    predict with, one to measure with.

    **On the cell key and the target.** The cell is
    :data:`~syn4rch.models.task.PERTURBATION_AXIS`, i.e. ``(class, polarity,
    target)`` from :func:`syn4rch.qc.keys.perturbation_key`, plus the
    population and the outcome. Every record in the pre-migration curated file
    has a null target, so today this key is exactly ``(population, class,
    polarity, outcome)``. After the v0.2.0 migration stores ``target_raw`` it
    becomes strictly finer -- a knockout of *Raptor* and a knockout of *p53* in
    the same population stop sharing a cell. That is intended, and it means the
    molecular numbers below are expected to move when the migration lands:
    their cells are currently pooling different genes, which is the finding
    that the mirroring diagnostic measures directly.
    """

    name = "perturbation_axis_cell_majority"
    description = (
        "modal direction of the (population, perturbation axis, outcome) cell, "
        "falling straight back to the global mode; the instrument that "
        "isolates the axis cell against the floor"
    )
    emits_proba = True
    feature_generation = FEATURE_GENERATION_V2
    features_used = (
        "population_key",
        "perturbation_class",
        "perturbation_polarity",
        "perturbation_target",
        "outcome_id",
    )
    chain = PERTURBATION_AXIS_CELL_BACKOFF


class PerturbationAxisConditionMajority(BackoffMajority):
    """:class:`PerturbationAxisMajority`, conditioned on sex, diet and challenge.

    Exists to answer one question with a number: does conditioning help, or
    does it merely fragment? Schema v0.2.0 added the condition fields because
    fourteen of eighteen apparent contradictions in the corpus turned out to be
    a male cohort compared against a female one, or chow against high-fat
    (DECISIONS.md D-023) -- which is a real finding about *comparability*. It
    does not follow that conditioning improves *prediction*: a finer cell is a
    smaller cell, and on 361 examples of which 153 state no condition at all,
    the conditioned cell can only be emptier than the cell it refines.

    The condition features are ``syn4rch.qc.keys.condition_key``'s tokens, in
    which an unknown slot reads ``unknown`` -- one token for "the curator did
    not address it" and "the paper does not say", pooled deliberately and
    documented there. Two records that are both ``unknown`` therefore share a
    cell by failing to state a condition, not by sharing one, and
    ``condition_specificity`` in the dataset summary is what stops that from
    being read as agreement.
    """

    name = "perturbation_axis_condition_majority"
    description = (
        "modal direction per (population, perturbation axis, outcome, sex, "
        "diet, challenge); each conditioned rung falls through to its "
        "unconditioned twin"
    )
    emits_proba = True
    feature_generation = FEATURE_GENERATION_V2
    features_used = (
        "population_key",
        "perturbation_class",
        "perturbation_polarity",
        "perturbation_target",
        "outcome_id",
    ) + CONDITION_FEATURE_FIELDS
    chain = PERTURBATION_AXIS_CONDITION_BACKOFF


# --------------------------------------------------------------------------- #
# neurotransmitter sign
# --------------------------------------------------------------------------- #

#: Transmitters this project declines to give a sign to, and why. Written out
#: rather than left as a gap in the rule table, so that an abstention is a
#: recorded decision instead of an oversight.
UNSIGNED_TRANSMITTERS: Dict[str, str] = {
    "dopamine": (
        "net effect depends on receptor subtype at the target (D1-like "
        "excitatory vs D2-like inhibitory); no single sign"
    ),
    "serotonin": (
        "at least 14 receptor subtypes with opposing couplings; no single sign"
    ),
    "histamine": (
        "H1/H2 excitatory, H3 autoinhibitory; net sign is target-dependent"
    ),
    "peptidergic": (
        "a release mode, not a sign; the peptide and its receptor decide"
    ),
    "mixed": (
        "the vocabulary entry says the population has both fast transmitter "
        "phenotypes; that is a statement that the sign is unknown"
    ),
    "unknown": (
        "the vocabulary does not state a transmitter for this population "
        "(SCHEMA.md 3.2: only if the paper or the vocab entry states it)"
    ),
}

#: The explicit rule table: (neurotransmitter, valence) -> predicted direction,
#: or ``None`` to abstain. Written out in full rather than derived from two
#: factor tables, because a rule table that has to be computed is a rule table
#: nobody reads.
#:
#: The two signed rows encode one assumption, stated plainly: that the readout
#: moves with the net excitatory drive the population exerts. For many
#: (population, outcome) pairs that is false -- AgRP activation *raises* food
#: intake despite being GABAergic, because the sign that matters is the sign of
#: the circuit path to the readout, not of the transmitter at its first synapse.
#: This baseline is therefore expected to be wrong in a structured way. That is
#: its value: the gap between it and :class:`PopulationOutcomeMajority` measures
#: how much of the corpus's direction signal a transmitter-level prior can
#: reach, and :class:`KGPathSign` is the baseline that would need a real edge
#: table to do better.
#:
#: ``ablation`` is treated as ``inhibition``: it removes the population's
#: output. Knowingly incomplete -- ablation also removes co-released peptides,
#: can trigger developmental compensation, and in adults often produces a
#: different phenotype from acute silencing.
#:
#: ``unclear`` valence always abstains (SCHEMA.md 3.4 forbids coercing it).
#:
#: ``no_change`` is never predicted by this table. The rule has no notion of
#: "no effect", so it is structurally unable to score on one of the three
#: classes; its macro-F1 is capped below 0.67 no matter how good the sign rule
#: is. That is reported, not corrected.
NEUROTRANSMITTER_VALENCE_RULES: Dict[Tuple[str, str], Optional[str]] = {
    ("GABA", "activation"): "decrease",
    ("GABA", "inhibition"): "increase",
    ("GABA", "ablation"): "increase",
    ("GABA", "unclear"): None,
    ("glutamate", "activation"): "increase",
    ("glutamate", "inhibition"): "decrease",
    ("glutamate", "ablation"): "decrease",
    ("glutamate", "unclear"): None,
    ("dopamine", "activation"): None,
    ("dopamine", "inhibition"): None,
    ("dopamine", "ablation"): None,
    ("dopamine", "unclear"): None,
    ("serotonin", "activation"): None,
    ("serotonin", "inhibition"): None,
    ("serotonin", "ablation"): None,
    ("serotonin", "unclear"): None,
    ("histamine", "activation"): None,
    ("histamine", "inhibition"): None,
    ("histamine", "ablation"): None,
    ("histamine", "unclear"): None,
    ("peptidergic", "activation"): None,
    ("peptidergic", "inhibition"): None,
    ("peptidergic", "ablation"): None,
    ("peptidergic", "unclear"): None,
    ("mixed", "activation"): None,
    ("mixed", "inhibition"): None,
    ("mixed", "ablation"): None,
    ("mixed", "unclear"): None,
    ("unknown", "activation"): None,
    ("unknown", "inhibition"): None,
    ("unknown", "ablation"): None,
    ("unknown", "unclear"): None,
}

#: Abstention reason codes, counted so the report can say *why* it abstained.
ABSTAIN_NO_POPULATION_ID = "no_population_id"
ABSTAIN_POPULATION_NOT_IN_VOCAB = "population_not_in_vocabulary"
ABSTAIN_VOCAB_SILENT = "vocabulary_states_no_transmitter"
ABSTAIN_UNSIGNED_TRANSMITTER = "transmitter_has_no_declared_sign"
ABSTAIN_TRANSMITTERS_DISAGREE = "co_transmitters_imply_opposite_signs"
ABSTAIN_VALENCE = "valence_absent_or_unclear"
ABSTAIN_NO_RULE = "no_rule_for_this_pair"


class NeurotransmitterSign(Baseline):
    """Predicts from the population's neurotransmitter crossed with valence.

    The transmitter is read from ``data/vocab/populations.json`` via
    ``syn4rch.vocab.load`` -- never from the record, and never inferred from a
    marker gene or a population name. **Where the vocabulary does not state a
    transmitter, this baseline abstains.** Concretely it abstains when:

    * ``population_id`` is null (SCHEMA.md 3.2 keeps such records rather than
      force-mapping them, so there is nothing to look up);
    * the id is not an entry in the populations vocabulary;
    * the entry's ``neurotransmitter`` is absent, empty, or ``unknown`` /
      ``mixed`` -- the vocabulary's own ways of saying it does not know;
    * every listed transmitter is one of :data:`UNSIGNED_TRANSMITTERS`;
    * two listed transmitters imply opposite signs;
    * ``valence`` is null or ``unclear``.

    ``fit`` uses no labels: this baseline learns nothing. It is fitted only so
    that it shares an interface with the others, and so that fitting it on a
    fold cannot accidentally be the thing that makes it work.

    ``predict_proba`` returns ``None``. The rule is a hard prior with no
    estimated frequency behind it; a confidence of 1.0 would be a fabricated
    number and the calibration plot would dutifully report it as badly
    calibrated rather than as absent.
    """

    name = "neurotransmitter_sign"
    description = (
        "declared (neurotransmitter x valence) rule table; abstains wherever "
        "the vocabulary does not state a transmitter"
    )
    emits_proba = False
    feature_generation = FEATURE_GENERATION_V1
    features_used = ("population_id", "valence")

    def __init__(self, vocab_loader: Optional[Callable[[str], Mapping[str, Any]]] = None):
        """``vocab_loader`` matches ``syn4rch.vocab.load``: ``load(kind) -> dict``.

        Injectable so a test can supply a vocabulary that is silent about a
        population, and so this baseline does not require the real vocabulary
        files to be present on disk.
        """
        super().__init__()
        self._populations: Dict[str, Any] = {}
        self.abstain_counts: Dict[str, int] = {}
        probes: List[str] = []
        if vocab_loader is None:
            try:
                from ..vocab import load as vocab_load

                vocab_loader = vocab_load
                probes.append("syn4rch.vocab.load: imported")
            except Exception as exc:  # pragma: no cover - vocab is stdlib-only
                self.availability = Availability(
                    available=False,
                    reason="could not import syn4rch.vocab: %s" % exc,
                    probes=tuple(probes + ["syn4rch.vocab.load: import failed"]),
                )
                return
        else:
            probes.append("vocab_loader: supplied by caller")
        try:
            populations = vocab_loader("populations")
        except Exception as exc:
            self.availability = Availability(
                available=False,
                reason=(
                    "the populations vocabulary could not be loaded, so no "
                    "transmitter can be read: %s" % exc
                ),
                probes=tuple(probes + ["load('populations'): raised"]),
            )
            return
        self._populations = dict(populations or {})
        probes.append("load('populations'): %d entries" % len(self._populations))
        self.availability = Availability(available=True, probes=tuple(probes))

    # -- transmitter lookup ------------------------------------------------

    def transmitters(self, population_id: Optional[str]) -> Optional[List[str]]:
        """Transmitters the vocabulary states, or ``None`` if it states none.

        ``None`` and ``[]`` are both "the vocabulary is silent"; they are
        distinguished only so the abstention reason can be specific.
        """
        if not population_id:
            return None
        entry = self._populations.get(population_id)
        if entry is None:
            return None
        values = entry.get("neurotransmitter") if isinstance(entry, Mapping) else None
        if values is None and isinstance(entry, dict):
            values = entry.get("neurotransmitter")
        if not values:
            return []
        if isinstance(values, (str, bytes)):
            return [str(values)]
        return [str(v) for v in values]

    def _abstain(self, reason: str) -> None:
        self.abstain_counts[reason] = self.abstain_counts.get(reason, 0) + 1

    def explain(self, features: Features) -> Dict[str, Any]:
        """The prediction and the full reasoning, for a report or a test."""
        population_id = features.population_id
        if not population_id:
            return {
                "prediction": None,
                "abstain_reason": ABSTAIN_NO_POPULATION_ID,
                "transmitters": None,
            }
        entry = self._populations.get(population_id)
        if entry is None:
            return {
                "prediction": None,
                "abstain_reason": ABSTAIN_POPULATION_NOT_IN_VOCAB,
                "transmitters": None,
            }
        transmitters = self.transmitters(population_id) or []
        if not transmitters:
            return {
                "prediction": None,
                "abstain_reason": ABSTAIN_VOCAB_SILENT,
                "transmitters": [],
            }
        valence = features.valence
        if valence in (None, "", "unclear"):
            return {
                "prediction": None,
                "abstain_reason": ABSTAIN_VALENCE,
                "transmitters": transmitters,
            }
        signed: List[str] = []
        unsigned: List[str] = []
        unruled: List[str] = []
        for transmitter in transmitters:
            key = (transmitter, valence)
            if key not in NEUROTRANSMITTER_VALENCE_RULES:
                unruled.append(transmitter)
                continue
            rule = NEUROTRANSMITTER_VALENCE_RULES[key]
            if rule is None:
                unsigned.append(transmitter)
            else:
                signed.append(rule)
        distinct = sorted(set(signed))
        if not distinct:
            reason = ABSTAIN_UNSIGNED_TRANSMITTER if unsigned else ABSTAIN_NO_RULE
            return {
                "prediction": None,
                "abstain_reason": reason,
                "transmitters": transmitters,
                "unsigned": unsigned,
                "unruled": unruled,
            }
        if len(distinct) > 1:
            return {
                "prediction": None,
                "abstain_reason": ABSTAIN_TRANSMITTERS_DISAGREE,
                "transmitters": transmitters,
                "candidate_predictions": distinct,
            }
        return {
            "prediction": distinct[0],
            "abstain_reason": None,
            "transmitters": transmitters,
            "signed_by": [t for t in transmitters if t not in UNSIGNED_TRANSMITTERS],
        }

    def _fit(self, examples: List[Example]) -> None:
        # Deliberately label-free. Reset only the abstention tally so that the
        # counts reported for a fold describe that fold.
        self.abstain_counts = {}

    def _predict_one(self, features: Features) -> Optional[str]:
        result = self.explain(features)
        if result["prediction"] is None and result["abstain_reason"]:
            self._abstain(result["abstain_reason"])
        return result["prediction"]

    def describe(self) -> Dict[str, Any]:
        out = super().describe()
        out["n_vocab_populations"] = len(self._populations)
        out["abstention_reasons"] = dict(sorted(self.abstain_counts.items()))
        out["signed_rules"] = {
            "%s x %s" % pair: label
            for pair, label in sorted(NEUROTRANSMITTER_VALENCE_RULES.items())
            if label is not None
        }
        out["unsigned_transmitters"] = dict(sorted(UNSIGNED_TRANSMITTERS.items()))
        out["cannot_predict"] = ["no_change"]
        return out


# --------------------------------------------------------------------------- #
# knowledge-graph path sign
# --------------------------------------------------------------------------- #

#: The edge-table format this baseline reads. Declared here so that the format
#: exists as a specification *before* any edge does, and so that whoever builds
#: the knowledge graph is not guessing what this baseline wants. Nothing in this
#: repository produces such a table today.
EDGE_TABLE_FORMAT: Dict[str, Any] = {
    "location": [
        "the edge_table= argument",
        "$SYN4RCH_KG_EDGES (path to a JSON or JSONL file)",
        "syn4rch.kg.edges, exposing EDGES or load_edges()",
        "<repo>/data/inferred/kg_edges.jsonl",
    ],
    "row": {
        "source": "node id: a populations-vocabulary id, or any intermediate node id",
        "target": "node id: an outcomes-vocabulary id, or any intermediate node id",
        "sign": "+1 or -1: whether activity at source raises or lowers target",
        "relation": "free text, e.g. 'projects_to', 'secretes', 'innervates'",
        "evidence_class": "anatomical | functional_association (SCHEMA.md 3.6)",
        "provenance": "where the edge came from; required, never 'inferred by model'",
    },
    "storage_note": (
        "such a table is derived data and belongs in data/inferred/ "
        "(SCHEMA.md 4), where it may never be cited as evidence"
    ),
}


@dataclass
class EdgeTable:
    """A signed directed graph, plus where it came from.

    ``source`` is recorded so that a prediction made from a table can always
    be traced to the table, and so that a report cannot print a KG number
    without saying which graph produced it.
    """

    edges: List[Mapping[str, Any]] = field(default_factory=list)
    source: str = "unspecified"

    def adjacency(self) -> Dict[str, List[Tuple[str, int]]]:
        out: Dict[str, List[Tuple[str, int]]] = {}
        for edge in self.edges:
            src = edge.get("source")
            dst = edge.get("target")
            sign = edge.get("sign")
            if src is None or dst is None or sign not in (1, -1, +1):
                continue
            out.setdefault(str(src), []).append((str(dst), int(sign)))
        return out

    def __len__(self) -> int:
        return len(self.edges)


def _probe_edge_table(
    edge_table: Optional[EdgeTable],
) -> Tuple[Optional[EdgeTable], List[str]]:
    """Look in each declared place for an edge table. Creates nothing."""
    probes: List[str] = []

    if edge_table is not None:
        probes.append(
            "edge_table= argument: supplied (%d edges, source=%r)"
            % (len(edge_table), edge_table.source)
        )
        return edge_table, probes
    probes.append("edge_table= argument: not supplied")

    env_path = os.environ.get("SYN4RCH_KG_EDGES")
    if env_path:
        if os.path.isfile(env_path):
            try:
                rows = _read_edge_file(env_path)
            except Exception as exc:
                probes.append("$SYN4RCH_KG_EDGES=%s: unreadable (%s)" % (env_path, exc))
            else:
                probes.append(
                    "$SYN4RCH_KG_EDGES=%s: %d edges" % (env_path, len(rows))
                )
                return EdgeTable(edges=rows, source=env_path), probes
        else:
            probes.append("$SYN4RCH_KG_EDGES=%s: no such file" % env_path)
    else:
        probes.append("$SYN4RCH_KG_EDGES: not set")

    try:
        from ..kg import edges as kg_edges  # type: ignore
    except Exception as exc:
        probes.append("syn4rch.kg.edges: not importable (%s)" % type(exc).__name__)
    else:
        rows = None
        if hasattr(kg_edges, "load_edges"):
            try:
                rows = list(kg_edges.load_edges())
            except Exception as exc:
                probes.append("syn4rch.kg.edges.load_edges(): raised (%s)" % exc)
        if rows is None and hasattr(kg_edges, "EDGES"):
            rows = list(getattr(kg_edges, "EDGES"))
        if rows is None:
            probes.append(
                "syn4rch.kg.edges: imported but exposes neither EDGES nor load_edges()"
            )
        else:
            probes.append("syn4rch.kg.edges: %d edges" % len(rows))
            return EdgeTable(edges=rows, source="syn4rch.kg.edges"), probes

    from .task import repo_root

    root = repo_root()
    if root is None:
        probes.append("data/inferred/kg_edges.jsonl: repository root not found")
    else:
        path = os.path.join(root, "data", "inferred", "kg_edges.jsonl")
        if os.path.isfile(path):
            try:
                rows = _read_edge_file(path)
            except Exception as exc:
                probes.append("%s: unreadable (%s)" % (path, exc))
            else:
                probes.append("%s: %d edges" % (path, len(rows)))
                return EdgeTable(edges=rows, source=path), probes
        else:
            probes.append("%s: no such file" % path)

    return None, probes


def _read_edge_file(path: str) -> List[Mapping[str, Any]]:
    """Read a JSON array or a JSONL file of edge rows. Never writes."""
    with open(path, "r", encoding="utf-8") as handle:
        text = handle.read()
    stripped = text.lstrip()
    if stripped.startswith("["):
        rows = json.loads(stripped)
        if not isinstance(rows, list):
            raise ValueError("expected a JSON array of edge rows")
        return [row for row in rows if isinstance(row, dict)]
    out: List[Mapping[str, Any]] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if isinstance(row, dict):
            out.append(row)
    return out


class KGPathSign(Baseline):
    """Sign propagation from a population to an outcome over a declared edge table.

    **Status in this repository: UNAVAILABLE.** ``src/syn4rch/kg/`` is empty,
    no ``$SYN4RCH_KG_EDGES`` is set, and ``data/inferred/kg_edges.jsonl`` does
    not exist. :attr:`availability` carries the reason and the list of places
    that were looked at. The baseline is left in the roster rather than deleted
    because its absence is a finding: the strongest mechanistic baseline
    available in principle cannot be run, so no claim in the benchmark report
    rests on it, and the report says so where its row would be.

    It would have been easy to make this row appear by writing a few plausible
    hypothalamic edges. That would be fabrication: an invented edge is an
    anatomical claim, ``docs/SCHEMA.md`` 2 rule 5 forbids creating records for
    sources that were not retrieved, 3.6 gives inferred wiring its own
    evidence class precisely so it cannot be passed off as measured, and 4
    forbids citing anything in ``data/inferred/`` as evidence. A baseline that
    beats the frequency baselines using edges nobody measured is worse than no
    baseline at all, because it would set the bar in the wrong place.

    The algorithm, for when a table exists: breadth-first search from
    ``population_id`` to ``outcome_id``; take the product of edge signs along
    each shortest path; if every shortest path agrees on the product, cross it
    with the valence sign (``activation`` +1, ``inhibition``/``ablation`` -1,
    ``unclear`` abstain) and emit ``increase`` or ``decrease``; if the shortest
    paths disagree in sign, abstain, because a graph that supports both answers
    is not evidence for either. Like :class:`NeurotransmitterSign`, it can
    never predict ``no_change`` and emits no probability.
    """

    name = "kg_path_sign"
    description = (
        "sign propagation over a declared KG edge table; UNAVAILABLE while no "
        "such table exists"
    )
    emits_proba = False
    feature_generation = FEATURE_GENERATION_V1
    features_used = ("population_id", "outcome_id", "valence")

    def __init__(self, edge_table: Optional[EdgeTable] = None, max_depth: int = 4):
        super().__init__()
        self.max_depth = max_depth
        self.table: Optional[EdgeTable] = None
        self._adjacency: Dict[str, List[Tuple[str, int]]] = {}
        self.abstain_counts: Dict[str, int] = {}
        table, probes = _probe_edge_table(edge_table)
        if table is None or len(table) == 0:
            self.availability = Availability(
                available=False,
                reason=(
                    "no knowledge-graph edge table exists in this repository, "
                    "so there is no graph to propagate a sign over. Edges will "
                    "not be invented to make this baseline runnable: an "
                    "invented edge is an unsupported anatomical claim "
                    "(SCHEMA.md 2 rule 5, 3.6, 4). Build "
                    "data/inferred/kg_edges.jsonl in the format given by "
                    "baselines.EDGE_TABLE_FORMAT, or point "
                    "$SYN4RCH_KG_EDGES at one, and this baseline becomes "
                    "available with no code change."
                    + (
                        " A table was found but it is empty."
                        if table is not None
                        else ""
                    )
                ),
                probes=tuple(probes),
            )
            return
        self.table = table
        self._adjacency = table.adjacency()
        self.availability = Availability(available=True, probes=tuple(probes))

    # -- sign propagation --------------------------------------------------

    def path_sign(self, source: Optional[str], target: Optional[str]) -> Optional[int]:
        """Product of edge signs along the shortest paths, or ``None``.

        ``None`` means either "no path within ``max_depth``" or "the shortest
        paths disagree". Both are abstentions and are counted separately.
        """
        if not source or not target or source == target:
            return None
        frontier = deque([(source, 1, 0)])
        seen_at_depth: Dict[Tuple[str, int], set] = {}
        found_depth: Optional[int] = None
        signs: set = set()
        while frontier:
            node, sign, depth = frontier.popleft()
            if depth >= self.max_depth:
                continue
            if found_depth is not None and depth >= found_depth:
                continue
            for neighbour, edge_sign in self._adjacency.get(node, ()):  # type: ignore[arg-type]
                next_sign = sign * edge_sign
                if neighbour == target:
                    if found_depth is None:
                        found_depth = depth + 1
                    if depth + 1 == found_depth:
                        signs.add(next_sign)
                    continue
                key = (neighbour, depth + 1)
                bucket = seen_at_depth.setdefault(key, set())
                if next_sign in bucket:
                    continue
                bucket.add(next_sign)
                frontier.append((neighbour, next_sign, depth + 1))
        if len(signs) != 1:
            return None
        return signs.pop()

    def explain(self, features: Features) -> Dict[str, Any]:
        if not self.available:
            return {"prediction": None, "abstain_reason": "baseline_unavailable"}
        valence = features.valence
        if valence in (None, "", "unclear"):
            return {"prediction": None, "abstain_reason": ABSTAIN_VALENCE}
        valence_sign = 1 if valence == "activation" else -1
        if not features.population_id:
            return {"prediction": None, "abstain_reason": ABSTAIN_NO_POPULATION_ID}
        sign = self.path_sign(features.population_id, features.outcome_id)
        if sign is None:
            return {"prediction": None, "abstain_reason": "no_unambiguous_signed_path"}
        net = sign * valence_sign
        return {
            "prediction": "increase" if net > 0 else "decrease",
            "abstain_reason": None,
            "path_sign": sign,
            "valence_sign": valence_sign,
        }

    def _fit(self, examples: List[Example]) -> None:
        self.abstain_counts = {}

    def _predict_one(self, features: Features) -> Optional[str]:
        result = self.explain(features)
        if result["prediction"] is None and result.get("abstain_reason"):
            self.abstain_counts[result["abstain_reason"]] = (
                self.abstain_counts.get(result["abstain_reason"], 0) + 1
            )
        return result["prediction"]

    def describe(self) -> Dict[str, Any]:
        out = super().describe()
        out["edge_table_source"] = self.table.source if self.table else None
        out["n_edges"] = len(self.table) if self.table else 0
        out["max_depth"] = self.max_depth
        out["expected_format"] = EDGE_TABLE_FORMAT
        out["abstention_reasons"] = dict(sorted(self.abstain_counts.items()))
        out["cannot_predict"] = ["no_change"]
        return out


# --------------------------------------------------------------------------- #
# roster
# --------------------------------------------------------------------------- #

def all_baseline_factories() -> "Dict[str, Callable[[], Baseline]]":
    """Name -> zero-argument constructor, in report order.

    Factories rather than instances: the evaluation harness constructs a fresh
    baseline for every fold, so that nothing fitted on fold *k* can survive
    into fold *k+1*. A shared instance is one of the easier ways to leak a
    training set into a test set without noticing.
    """
    return {
        GlobalMajority.name: GlobalMajority,
        OutcomeMajority.name: OutcomeMajority,
        PopulationOutcomeMajority.name: PopulationOutcomeMajority,
        PerturbationAxisCellMajority.name: PerturbationAxisCellMajority,
        PerturbationAxisMajority.name: PerturbationAxisMajority,
        PerturbationAxisConditionMajority.name: PerturbationAxisConditionMajority,
        NeurotransmitterSign.name: NeurotransmitterSign,
        KGPathSign.name: KGPathSign,
    }


def baselines_by_feature_generation() -> Dict[str, List[str]]:
    """``feature generation -> baseline names``, in report order.

    Read from each class's declared :attr:`Baseline.feature_generation`, so the
    before/after table in the report is assembled from the declaration rather
    than from a hand-written list that could disagree with it.
    """
    out: Dict[str, List[str]] = {
        FEATURE_GENERATION_V1: [],
        FEATURE_GENERATION_V2: [],
    }
    for name, factory in all_baseline_factories().items():
        generation = getattr(factory, "feature_generation", FEATURE_GENERATION_V1)
        out.setdefault(generation, []).append(name)
    return out


def build_all() -> "List[Baseline]":
    """One unfitted instance of each baseline, for availability reporting."""
    return [factory() for factory in all_baseline_factories().values()]
