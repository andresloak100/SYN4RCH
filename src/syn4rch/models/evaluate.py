"""Grouped cross-validation, with the grouping label made structural.

``docs/SCHEMA.md`` 5 ends with a sentence that governs this whole module:

    Reporting a single cross-validated number without naming the grouping is
    considered an invalid result in this project.

A rule like that is usually implemented as a convention and then broken by the
first person in a hurry. Here it is implemented as a type. :class:`Estimate`
and :class:`MetricBundle` both **require** a ``grouping`` at construction, both
validate it against :data:`~syn4rch.models.task.GROUPINGS`, and both put it in
every dict they emit. There is no code path in this package that turns a number
from ``syn4rch.models.metrics`` into JSON or Markdown without passing through
one of them, because the metric functions deliberately return bare floats.
Constructing an unlabelled metric raises :class:`UnlabelledMetricError`.

What is reported, and why separately
------------------------------------
Five things, never blended into one headline number:

* **direction accuracy** -- over examples the baseline answered, with coverage
  printed beside it, because an abstaining baseline's accuracy is a different
  quantity from a guessing one's.
* **macro-F1** -- because this corpus is expected to be dominated by one or two
  directions, and accuracy on a skewed three-class problem is mostly a
  measurement of the skew. Reported over classes present in the fold.
* **ranking** -- one-vs-rest macro ROC AUC, only where a baseline emits
  probabilities that actually vary. A constant scorer gets ``None`` and a
  reason, not 0.5.
* **calibration** -- reliability bins plus the multiclass Brier score. Kept
  apart from accuracy: a baseline can be accurate and badly calibrated, and a
  hypothesis-generation platform is used by ranking and thresholding its
  outputs, so calibration is the property that decides whether its numbers can
  be acted on.
* **uncertainty** -- a percentile bootstrap confidence interval resampled over
  **groups**, never rows (see :func:`syn4rch.models.metrics.bootstrap_group_ci`).

And beside all five, one number that is not a metric at all: the fraction of
test rows whose **perturbation identity** never appeared in training
(``docs/SCHEMA.md`` 5b, via
:func:`syn4rch.models.diagnostics.perturbation_identity_novelty`). Grouping
guarantees that no study crosses a split; it guarantees nothing about whether
the held-out perturbation was learnable from the training side. Where that
fraction is high, the accuracy beside it is substantially a measurement of the
label prior, and SCHEMA.md 5b is explicit that this is a statement about how
the field publishes -- one lab, one gene, one paper -- and not a modelling
failure to be tuned away.

Fold generation
---------------
Fold generators are imported from ``syn4rch.feasibility.splits`` when that
module exists and exposes something this harness can validate; otherwise a
built-in grouped fold generator is used and every emitted bundle says which was
used, in ``fold_source``. Both states are normal: that module is written by
another component and may legitimately be absent, half-written, or export a
different signature. A probe that cannot be *validated* is not used, because a
fold generator that silently leaks a group is worse than no fold generator.

Python 3.10 compatible, standard library only.
"""

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

from . import metrics as M
from ._compat import optional_dependency_report
from .baselines import Baseline, all_baseline_factories
from .diagnostics import dataset_diagnostics, perturbation_identity_novelty
from .task import (
    GROUPINGS,
    LABELS,
    DEFAULT_GROUPINGS,
    NULL_FEATURE,
    Example,
    Grouping,
    TaskDataset,
)

__all__ = [
    "UnlabelledMetricError",
    "ComparisonSlice",
    "Comparison",
    "compare_under_grouping",
    "COMPARISON_SLICES",
    "COMPARISON_REFERENCES",
    "HEADLINE_COMPARISON",
    "MIN_ROWS_FOR_METRIC",
    "MIN_GROUPS_FOR_CI",
    "MIN_GROUPS_FOR_CV",
    "MIN_N_PER_RELIABILITY_BIN",
    "Estimate",
    "MetricBundle",
    "FoldAssignment",
    "FoldProvider",
    "resolve_fold_provider",
    "builtin_grouped_folds",
    "BUILTIN_SOURCE",
    "cross_validate",
    "BenchmarkResult",
    "run_benchmark",
]


# --------------------------------------------------------------------------- #
# the grouping guard
# --------------------------------------------------------------------------- #

class UnlabelledMetricError(ValueError):
    """A metric was about to be emitted without a declared grouping label.

    SCHEMA.md 5 calls such a number an invalid result. Raising here is the
    mechanism that makes it hard to produce one by accident.
    """


def _require_grouping(grouping: Any) -> str:
    if grouping is None or grouping == "":
        raise UnlabelledMetricError(
            "a metric may not be emitted without a grouping label "
            "(SCHEMA.md 5: an unnamed cross-validated number is an invalid "
            "result). Declared groupings: %s" % ", ".join(sorted(GROUPINGS))
        )
    if not isinstance(grouping, str):
        raise UnlabelledMetricError(
            "grouping must be the name of a declared grouping, got %r" % (grouping,)
        )
    if grouping not in GROUPINGS:
        raise UnlabelledMetricError(
            "undeclared grouping %r: a grouping must be registered in "
            "task.GROUPINGS so that its leakage properties are documented. "
            "Declared: %s" % (grouping, ", ".join(sorted(GROUPINGS)))
        )
    return grouping


# --------------------------------------------------------------------------- #
# how small is too small
# --------------------------------------------------------------------------- #

#: Below this many scored rows a metric is reported with an explicit "too few
#: rows to carry meaning" caveat rather than printed bare. Not a hard floor:
#: suppressing the number entirely would hide how little data there is.
MIN_ROWS_FOR_METRIC = 20

#: Below this many independent groups a confidence interval is reported with a
#: caveat: a cluster bootstrap over four labs is resampling four things.
MIN_GROUPS_FOR_CI = 5

#: Grouped cross-validation needs at least two groups to have a test set that
#: is independent of the training set at all.
MIN_GROUPS_FOR_CV = 2

#: Reliability bins holding fewer than this are marked not meaningful.
MIN_N_PER_RELIABILITY_BIN = 10


# --------------------------------------------------------------------------- #
# labelled containers
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Estimate:
    """One number, inseparable from the grouping it was computed under.

    ``value`` of ``None`` means the metric was undefined, and
    ``undefined_reason`` says why. ``caveats`` is non-empty whenever n is too
    small for the number to carry meaning; the renderer prints the caveat
    *next to* the number rather than printing it bare or hiding it.
    """

    name: str
    grouping: str
    value: Optional[float]
    n_rows: int
    n_groups: int
    ci: Optional[Dict[str, Any]] = None
    caveats: Tuple[str, ...] = ()
    undefined_reason: Optional[str] = None
    detail: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        _require_grouping(self.grouping)

    @property
    def meaningful(self) -> bool:
        return self.value is not None and not self.caveats

    def to_dict(self) -> Dict[str, Any]:
        return {
            "metric": self.name,
            # Present in every emitted metric dict, by construction.
            "grouping": self.grouping,
            "grouping_description": GROUPINGS[self.grouping].description,
            "value": self.value,
            "n_rows": self.n_rows,
            "n_groups": self.n_groups,
            "ci": self.ci,
            "caveats": list(self.caveats),
            "undefined_reason": self.undefined_reason,
            "meaningful": self.meaningful,
            "detail": dict(self.detail),
        }

    def render(self, places: int = 3) -> str:
        """Human-readable form that always carries n and any caveat."""
        if self.value is None:
            return "n/a (%s)" % (self.undefined_reason or "undefined")
        text = ("%." + str(places) + "f") % self.value
        if self.ci and self.ci.get("ci_low") is not None:
            text += " [%.3f, %.3f]" % (self.ci["ci_low"], self.ci["ci_high"])
        text += " (n=%d rows / %d %s groups)" % (
            self.n_rows,
            self.n_groups,
            self.grouping,
        )
        if self.caveats:
            text += " -- " + "; ".join(self.caveats)
        return text


def make_estimate(
    name: str,
    grouping: str,
    value: Optional[float],
    n_rows: int,
    n_groups: int,
    ci: Optional[Dict[str, Any]] = None,
    undefined_reason: Optional[str] = None,
    detail: Optional[Dict[str, Any]] = None,
    extra_caveats: Sequence[str] = (),
) -> Estimate:
    """Build an :class:`Estimate`, attaching the small-n caveats automatically.

    Centralising the caveats means a new metric cannot be added that forgets
    them.
    """
    caveats: List[str] = [str(c) for c in extra_caveats]
    if value is not None and n_rows < MIN_ROWS_FOR_METRIC:
        caveats.append(
            "n=%d scored rows is below %d: this number does not carry meaning"
            % (n_rows, MIN_ROWS_FOR_METRIC)
        )
    if value is not None and n_groups < MIN_GROUPS_FOR_CI:
        caveats.append(
            "only %d independent %s group(s): treat any interval as indicative"
            % (n_groups, grouping)
        )
    if ci is not None and ci.get("ci_low") is None and ci.get("note"):
        caveats.append(str(ci["note"]))
    return Estimate(
        name=name,
        grouping=grouping,
        value=value,
        n_rows=n_rows,
        n_groups=n_groups,
        ci=ci,
        caveats=tuple(caveats),
        undefined_reason=undefined_reason,
        detail=dict(detail or {}),
    )


@dataclass
class MetricBundle:
    """Everything measured for one (baseline, grouping) pair.

    ``status`` is one of ``ok``, ``empty_dataset``, ``insufficient_groups`` or
    ``baseline_unavailable``. In the last three there are no numbers at all,
    and the bundle exists to say so: an evaluation that produced nothing must
    produce a statement that it produced nothing, not an absent row.
    """

    baseline: str
    grouping: str
    status: str = "ok"
    status_detail: Optional[str] = None
    fold_source: str = "unset"
    fold_source_note: Optional[str] = None
    n_folds: int = 0
    n_rows_scored: int = 0
    n_rows_total: int = 0
    n_rows_ungroupable: int = 0
    n_rows_not_tested: int = 0
    n_rows_abstained: int = 0
    n_groups: int = 0
    estimates: Dict[str, Estimate] = field(default_factory=dict)
    calibration: Dict[str, Any] = field(default_factory=dict)
    confusion: Dict[str, Dict[str, int]] = field(default_factory=dict)
    per_class: Dict[str, Dict[str, Optional[float]]] = field(default_factory=dict)
    per_fold_accuracy: List[Optional[float]] = field(default_factory=list)
    abstention_reasons: Dict[str, int] = field(default_factory=dict)
    #: How often each rung of a baseline's backoff chain was actually used,
    #: summed over folds. A baseline named for its primary cell that in fact
    #: spent its predictions one rung down is not the model its name
    #: describes, and a candidate compared against it as though it were would
    #: be compared against the wrong bar.
    backoff_usage: Dict[str, int] = field(default_factory=dict)
    #: SCHEMA.md 5b: the fraction of test rows whose perturbation identity
    #: never appeared in training, under *these* folds. Required beside every
    #: cross-validated number in this project, because grouping stops a study
    #: crossing a split but does not establish that the held-out perturbation
    #: was learnable from the training side at all.
    unseen_perturbation_identity: Dict[str, Any] = field(default_factory=dict)
    baseline_detail: Dict[str, Any] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        _require_grouping(self.grouping)

    @property
    def coverage(self) -> Optional[float]:
        total = self.n_rows_scored + self.n_rows_abstained
        if total == 0:
            return None
        return self.n_rows_scored / total

    def add(self, estimate: Estimate) -> None:
        if estimate.grouping != self.grouping:
            raise UnlabelledMetricError(
                "estimate %r was computed under grouping %r but is being "
                "filed under %r; a mislabelled metric is as invalid as an "
                "unlabelled one (SCHEMA.md 5)"
                % (estimate.name, estimate.grouping, self.grouping)
            )
        self.estimates[estimate.name] = estimate

    def to_dict(self) -> Dict[str, Any]:
        grouping = _require_grouping(self.grouping)
        declaration: Grouping = GROUPINGS[grouping]
        return {
            "baseline": self.baseline,
            # Every emitted metric dict carries its grouping, at the top level
            # and again inside each estimate.
            "grouping": grouping,
            "grouping_description": declaration.description,
            "grouping_schema_clause": declaration.schema_clause,
            "grouping_is_leakage_prone": declaration.leakage_prone,
            "status": self.status,
            "status_detail": self.status_detail,
            "fold_source": self.fold_source,
            "fold_source_note": self.fold_source_note,
            "n_folds": self.n_folds,
            "n_rows_total": self.n_rows_total,
            "n_rows_scored": self.n_rows_scored,
            "n_rows_abstained": self.n_rows_abstained,
            "n_rows_ungroupable": self.n_rows_ungroupable,
            "n_rows_not_tested": self.n_rows_not_tested,
            "n_groups": self.n_groups,
            "coverage": self.coverage,
            "metrics": {
                name: estimate.to_dict()
                for name, estimate in sorted(self.estimates.items())
            },
            "calibration": dict(self.calibration),
            "confusion_matrix": self.confusion,
            "per_class": self.per_class,
            "per_fold_accuracy": list(self.per_fold_accuracy),
            "abstention_reasons": dict(sorted(self.abstention_reasons.items())),
            "backoff_usage": dict(sorted(self.backoff_usage.items())),
            "unseen_perturbation_identity": dict(self.unseen_perturbation_identity),
            "baseline_detail": dict(self.baseline_detail),
            "notes": list(self.notes),
        }


# --------------------------------------------------------------------------- #
# fold generation
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class FoldAssignment:
    """Folds for one grouping, plus where they came from.

    ``source`` and ``note`` are per-grouping, not per-run: a provider can
    legitimately serve three groupings from ``syn4rch.feasibility.splits`` and
    the other two from the built-in generator, and a bundle must say which of
    the two produced *its* number.
    """

    folds: Tuple[Tuple[Tuple[int, ...], Tuple[int, ...]], ...]
    source: str
    note: str


def builtin_grouped_folds(
    group_values: Sequence[Optional[str]], n_folds: int = 5
) -> List[Tuple[List[int], List[int]]]:
    """Deterministic grouped k-fold. No group spans train and test.

    Groups are sorted by size descending (ties broken by key, so the result
    does not depend on input order) and assigned greedily to the currently
    smallest fold. This is the standard balanced-partition heuristic; it is not
    optimal, but it is deterministic, which matters more here than balance.

    ``n_folds`` is reduced to the number of groups, so a grouping with three
    groups yields three folds -- i.e. leave-one-group-out -- rather than two
    empty folds and a hidden error.

    Rows whose group key is ``None`` are not placed in any fold. The caller
    reports how many there were (DECISIONS.md D-009).
    """
    if n_folds < 2:
        raise ValueError("n_folds must be >= 2")
    buckets: Dict[str, List[int]] = {}
    for index, value in enumerate(group_values):
        if value is None:
            continue
        buckets.setdefault(value, []).append(index)
    if len(buckets) < MIN_GROUPS_FOR_CV:
        return []
    k = min(n_folds, len(buckets))
    assignment: List[List[int]] = [[] for _ in range(k)]
    ordered = sorted(buckets.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    for _, indices in ordered:
        target = min(range(k), key=lambda i: (len(assignment[i]), i))
        assignment[target].extend(indices)
    folds: List[Tuple[List[int], List[int]]] = []
    for i in range(k):
        test = sorted(assignment[i])
        train = sorted(
            index for j in range(k) if j != i for index in assignment[j]
        )
        if test and train:
            folds.append((train, test))
    return folds


def _as_fold_pairs(raw) -> Optional[List[Tuple[List[int], List[int]]]]:
    """Coerce an index-pair fold object into ``[(train, test), ...]``."""
    try:
        items = list(raw)
    except TypeError:
        return None
    folds: List[Tuple[List[int], List[int]]] = []
    for item in items:
        train = test = None
        if isinstance(item, (tuple, list)) and len(item) == 2:
            train, test = item[0], item[1]
        else:
            for train_attr, test_attr in (
                ("train", "test"),
                ("train_idx", "test_idx"),
                ("train_indices", "test_indices"),
            ):
                if hasattr(item, train_attr) and hasattr(item, test_attr):
                    train = getattr(item, train_attr)
                    test = getattr(item, test_attr)
                    break
        if train is None or test is None:
            return None
        try:
            folds.append((sorted(int(i) for i in train), sorted(int(i) for i in test)))
        except (TypeError, ValueError):
            return None
    return folds or None


def _validate_folds(
    folds: Sequence[Tuple[Sequence[int], Sequence[int]]],
    group_values: Sequence[Optional[str]],
) -> Optional[str]:
    """``None`` if the folds are safe to use, else the reason they are not.

    The checks are exactly the properties this harness depends on, run against
    the **real data** rather than a synthetic probe:

    1. at least two folds with both sides non-empty;
    2. every index in range, and no index in both sides of a fold;
    3. no row tested more than once (a row tested twice is counted twice);
    4. **no group appears on both sides of a fold** -- the leakage property
       that is the entire reason for grouping.

    Failure 4 is a refusal, not a warning. Using a generator that leaks would
    make the benchmark's central claim -- that these numbers are grouped --
    false while still printing the grouping's name beside them, which is worse
    than reporting an ungrouped number honestly.
    """
    usable = [(list(tr), list(te)) for tr, te in folds if len(tr) and len(te)]
    if len(usable) < 2:
        return "fewer than 2 folds have both a train and a test side"
    n = len(group_values)
    tested: List[int] = []
    for train, test in usable:
        if any(i < 0 or i >= n for i in train + test):
            return "a fold contains an out-of-range row index"
        if set(train) & set(test):
            return "a fold's train and test rows overlap"
        leaked = {group_values[i] for i in train} & {group_values[i] for i in test}
        leaked.discard(None)
        if leaked:
            return (
                "group(s) %s appear on both sides of a fold; refused "
                "(SCHEMA.md 5)"
                % ", ".join(sorted(str(g) for g in list(leaked)[:3]))
            )
        tested.extend(test)
    if len(tested) != len(set(tested)):
        return "some rows are tested in more than one fold"
    return None


#: Generic, group-value-style names this harness will also look for in
#: ``syn4rch.feasibility.splits``, for the case where that module grows a
#: different API than the record-based one it has today.
_SPLITS_VALUE_CANDIDATES: Tuple[str, ...] = (
    "grouped_folds",
    "grouped_kfold",
    "group_kfold",
    "grouped_cv_folds",
    "make_folds",
    "generate_folds",
    "kfold_by_group",
    "split_by_group",
)

#: Call shapes tried against a generic candidate, in order.
_CALL_SHAPES: Tuple[str, ...] = (
    "fn(group_values, n_folds=k)",
    "fn(group_values, k)",
    "fn(groups=group_values, n_folds=k)",
    "fn(group_values)",
)


def _invoke(fn: Callable, shape: str, group_values, k):
    if shape == "fn(group_values, n_folds=k)":
        return fn(group_values, n_folds=k)
    if shape == "fn(group_values, k)":
        return fn(group_values, k)
    if shape == "fn(groups=group_values, n_folds=k)":
        return fn(groups=group_values, n_folds=k)
    return fn(group_values)


def _foldset_to_indices(
    foldset: Any, record_index: Dict[str, int]
) -> Optional[List[Tuple[List[int], List[int]]]]:
    """Map a ``syn4rch.feasibility.splits.FoldSet`` onto example indices.

    Record ids that are not among this dataset's examples are dropped: that
    module's folds are built over *all* curated records, including the ones
    this task cannot label (``not_reported``, ``non_monotonic``, Tier B), so a
    fold legitimately names records that are not examples here. Dropping them
    is the correct translation, and the resulting shrinkage is reported as
    rows-not-tested rather than hidden.
    """
    folds_attr = getattr(foldset, "folds", None)
    if folds_attr is None:
        return None
    out: List[Tuple[List[int], List[int]]] = []
    for fold in folds_attr:
        train_ids = getattr(fold, "train_record_ids", None)
        test_ids = getattr(fold, "test_record_ids", None)
        if train_ids is None or test_ids is None:
            return None
        train = sorted(record_index[r] for r in train_ids if r in record_index)
        test = sorted(record_index[r] for r in test_ids if r in record_index)
        if train and test:
            out.append((train, test))
    return out


class FoldProvider(object):
    """Resolves folds for a grouping, preferring ``syn4rch.feasibility.splits``.

    Constructed by :func:`resolve_fold_provider`. Holds whatever it managed to
    find in that module (which may be nothing) and falls back, per grouping,
    to :func:`builtin_grouped_folds`. Every :class:`FoldAssignment` it returns
    names its own source, so a report can state per row which generator
    produced the folds behind that number.
    """

    def __init__(
        self,
        record_api: Optional[Dict[str, Tuple[str, Callable]]] = None,
        value_api: Optional[Tuple[str, Callable, str]] = None,
        resolution_note: str = "",
    ) -> None:
        self.record_api = dict(record_api or {})
        self.value_api = value_api
        self.resolution_note = resolution_note

    @property
    def source(self) -> str:
        parts: List[str] = []
        if self.record_api:
            parts.append(
                "syn4rch.feasibility.splits:{%s}"
                % ",".join(sorted(name for name, _ in self.record_api.values()))
            )
        if self.value_api:
            parts.append("syn4rch.feasibility.splits:%s" % self.value_api[0])
        if not parts:
            return BUILTIN_SOURCE
        parts.append(BUILTIN_SOURCE)
        return " + ".join(parts)

    @property
    def note(self) -> str:
        return self.resolution_note

    def folds_for(
        self, dataset: TaskDataset, grouping: str, n_folds: int
    ) -> FoldAssignment:
        grouping = _require_grouping(grouping)
        group_values = dataset.group_values(grouping)
        rejections: List[str] = []

        entry = self.record_api.get(grouping)
        if entry is not None and dataset.records:
            name, fn = entry
            try:
                record_index = dataset.record_index()
                foldset = fn(dataset.records)
                folds = _foldset_to_indices(foldset, record_index)
            except Exception as exc:
                folds = None
                rejections.append(
                    "%s raised %s: %s" % (name, type(exc).__name__, exc)
                )
            if folds is not None:
                failure = _validate_folds(folds, group_values)
                if failure is None:
                    return FoldAssignment(
                        folds=tuple(
                            (tuple(tr), tuple(te)) for tr, te in folds
                        ),
                        source="syn4rch.feasibility.splits:%s" % name,
                        note=(
                            "folds come from syn4rch.feasibility.splits.%s "
                            "(leave-one-group-out over the curated records), "
                            "validated on this dataset for no train/test group "
                            "overlap" % name
                        ),
                    )
                rejections.append("%s rejected: %s" % (name, failure))
            elif not rejections:
                rejections.append(
                    "%s returned something this harness could not map to "
                    "example indices" % name
                )

        if self.value_api is not None:
            name, fn, shape = self.value_api
            try:
                folds = _as_fold_pairs(_invoke(fn, shape, list(group_values), n_folds))
            except Exception as exc:
                folds = None
                rejections.append(
                    "%s via %s raised %s: %s" % (name, shape, type(exc).__name__, exc)
                )
            if folds is not None:
                failure = _validate_folds(folds, group_values)
                if failure is None:
                    return FoldAssignment(
                        folds=tuple((tuple(tr), tuple(te)) for tr, te in folds),
                        source="syn4rch.feasibility.splits:%s" % name,
                        note=(
                            "folds come from syn4rch.feasibility.splits.%s, "
                            "called as %s, validated on this dataset"
                            % (name, shape)
                        ),
                    )
                rejections.append("%s via %s rejected: %s" % (name, shape, failure))

        folds = builtin_grouped_folds(group_values, n_folds)
        why = (
            "syn4rch.feasibility.splits has no counterpart for the %r grouping"
            % grouping
            if entry is None and self.value_api is None
            else "; ".join(rejections) or "no usable generator was found there"
        )
        return FoldAssignment(
            folds=tuple((tuple(tr), tuple(te)) for tr, te in folds),
            source=BUILTIN_SOURCE,
            note=(
                "folds come from this module's own deterministic grouped "
                "k-fold, because %s. Fold assignment therefore does not match "
                "that module's, so this number is not directly comparable with "
                "one computed there." % why
            ),
        )


#: Source string for the built-in generator.
BUILTIN_SOURCE = "syn4rch.models.evaluate:builtin_grouped_folds"


def resolve_fold_provider(verbose: bool = False) -> FoldProvider:
    """Import fold generators from ``syn4rch.feasibility.splits`` if possible.

    Degrades cleanly through every state that module can legitimately be in
    while it is written concurrently: absent, importable but exporting nothing
    this harness recognises, exporting the record-based ``FoldSet`` API it has
    today, or exporting a group-value API. Nothing is imported from it at
    module-import time, so ``import syn4rch.models.evaluate`` succeeds in all
    four states.

    Anything found is validated against the real data in
    :meth:`FoldProvider.folds_for` before a single number is computed from it.
    """
    try:
        from ..feasibility import splits as splits_module  # type: ignore
    except Exception as exc:
        return FoldProvider(
            resolution_note=(
                "syn4rch.feasibility.splits is not importable (%s: %s), so "
                "every grouping uses this module's own deterministic grouped "
                "k-fold." % (type(exc).__name__, exc)
            )
        )

    record_api: Dict[str, Tuple[str, Callable]] = {}
    for grouping, declaration in GROUPINGS.items():
        name = declaration.feasibility_name
        if not name:
            continue
        fn = getattr(splits_module, name, None)
        if callable(fn):
            record_api[grouping] = (name, fn)

    value_api: Optional[Tuple[str, Callable, str]] = None
    probe = ["g0", "g0", "g0", "g1", "g1", "g1", "g2", "g2", "g2", "g3", "g3", "g3"]
    for name in _SPLITS_VALUE_CANDIDATES:
        fn = getattr(splits_module, name, None)
        if not callable(fn):
            continue
        for shape in _CALL_SHAPES:
            try:
                folds = _as_fold_pairs(_invoke(fn, shape, list(probe), 4))
            except Exception:
                continue
            if folds and _validate_folds(folds, probe) is None:
                value_api = (name, fn, shape)
                break
        if value_api:
            break

    if record_api or value_api:
        found = sorted(name for name, _ in record_api.values())
        note = (
            "syn4rch.feasibility.splits imported; fold generators found: %s%s. "
            "Each is validated on the actual dataset (no train/test group "
            "overlap, each row tested at most once) before any number is "
            "computed from it; a grouping with no counterpart there uses this "
            "module's own grouped k-fold and says so in its own note."
            % (
                ", ".join(found) or "none of the record-based names",
                (
                    " plus the group-value generator %r" % value_api[0]
                    if value_api
                    else ""
                ),
            )
        )
    else:
        note = (
            "syn4rch.feasibility.splits imported but exposed no fold generator "
            "this harness recognised, so every grouping uses this module's own "
            "deterministic grouped k-fold."
        )
    if verbose:  # pragma: no cover - diagnostic
        print(note)
    return FoldProvider(
        record_api=record_api, value_api=value_api, resolution_note=note
    )


# --------------------------------------------------------------------------- #
# cross-validation
# --------------------------------------------------------------------------- #

def _fit_predict_fold(
    factory: Callable[[], Baseline],
    train: Sequence[Example],
    test: Sequence[Example],
) -> Tuple[Baseline, List[Optional[str]], Optional[List[Optional[Dict[str, float]]]]]:
    """A fresh baseline per fold: nothing fitted on fold k survives into k+1."""
    model = factory()
    model.fit(train)
    features = [e.features for e in test]
    predictions = model.predict(features)
    probas = model.predict_proba(features)
    return model, predictions, probas


def cross_validate(
    factory: Callable[[], Baseline],
    dataset: TaskDataset,
    grouping: str,
    *,
    n_folds: int = 5,
    n_bootstrap: int = 1000,
    seed: int = 20240101,
    n_reliability_bins: int = 10,
    fold_provider: Optional[FoldProvider] = None,
) -> MetricBundle:
    """Grouped cross-validation of one baseline under one named grouping.

    Predictions are pooled out-of-fold: every groupable row is predicted
    exactly once, by a model that never saw its group. The metrics are then
    computed on that pooled set, and the bootstrap resamples the *groups* of
    that pooled set. Per-fold accuracies are also reported, because a pooled
    number with a wide fold spread behind it is a different finding from the
    same number with a narrow one.
    """
    grouping = _require_grouping(grouping)
    provider = fold_provider or resolve_fold_provider()
    probe = factory()
    bundle = MetricBundle(
        baseline=probe.name,
        grouping=grouping,
        fold_source=provider.source,
        fold_source_note=provider.note,
        n_rows_total=len(dataset),
        baseline_detail=probe.describe(),
    )

    if not probe.available:
        bundle.status = "baseline_unavailable"
        bundle.status_detail = probe.availability.reason
        bundle.notes.append(
            "No numbers are reported for this baseline because it cannot run. "
            "Its row in the comparison is left empty on purpose."
        )
        return bundle

    if dataset.is_empty:
        bundle.status = "empty_dataset"
        bundle.status_detail = (
            "the dataset contains no labellable examples, so there is nothing "
            "to cross-validate"
        )
        return bundle

    group_values = dataset.group_values(grouping)
    bundle.n_rows_ungroupable = sum(1 for v in group_values if v is None)
    if bundle.n_rows_ungroupable:
        bundle.notes.append(
            "%d of %d rows have no %s group key and were excluded from this "
            "evaluation rather than pooled into one anonymous group "
            "(SCHEMA.md 5, DECISIONS.md D-009)."
            % (bundle.n_rows_ungroupable, len(group_values), grouping)
        )

    distinct = {v for v in group_values if v is not None}
    bundle.n_groups = len(distinct)
    if len(distinct) < MIN_GROUPS_FOR_CV:
        bundle.status = "insufficient_groups"
        bundle.status_detail = (
            "%d distinct %s group(s) among %d groupable rows: grouped "
            "cross-validation needs at least %d, because with one group every "
            "test row shares its group with the training set"
            % (
                len(distinct),
                grouping,
                len(group_values) - bundle.n_rows_ungroupable,
                MIN_GROUPS_FOR_CV,
            )
        )
        return bundle

    assignment = provider.folds_for(dataset, grouping, n_folds)
    folds = assignment.folds
    bundle.fold_source = assignment.source
    bundle.fold_source_note = assignment.note
    if not folds:
        bundle.status = "insufficient_groups"
        bundle.status_detail = (
            "the fold generator (%s) produced no usable folds for %d %s "
            "group(s)" % (assignment.source, len(distinct), grouping)
        )
        return bundle
    bundle.n_folds = len(folds)

    bundle.unseen_perturbation_identity = perturbation_identity_novelty(
        dataset, folds
    )
    unseen_fraction = bundle.unseen_perturbation_identity.get("fraction_unseen")
    if unseen_fraction is not None and unseen_fraction > 0.5:
        bundle.notes.append(
            "%.1f%% of the test rows under this grouping carry a perturbation "
            "identity that never appears in their fold's training set "
            "(SCHEMA.md 5b). For those rows there was nothing specific to have "
            "learned, so the accuracy beside them is substantially a "
            "measurement of the label prior."
            % (100.0 * unseen_fraction)
        )

    pooled_true: List[str] = []
    pooled_pred: List[Optional[str]] = []
    pooled_proba: List[Optional[Dict[str, float]]] = []
    pooled_group: List[Optional[str]] = []
    any_proba = False

    for train_idx, test_idx in folds:
        train = dataset.subset(train_idx)
        test = dataset.subset(test_idx)
        model, predictions, probas = _fit_predict_fold(factory, train, test)
        fold_true = [e.label for e in test]
        bundle.per_fold_accuracy.append(M.accuracy(fold_true, predictions))
        pooled_true.extend(fold_true)
        pooled_pred.extend(predictions)
        pooled_group.extend(group_values[i] for i in test_idx)
        if probas is None:
            pooled_proba.extend([None] * len(test))
        else:
            any_proba = True
            pooled_proba.extend(probas)
        detail = model.describe()
        for reason, count in (detail.get("abstention_reasons") or {}).items():
            bundle.abstention_reasons[reason] = (
                bundle.abstention_reasons.get(reason, 0) + count
            )
        for level, count in (detail.get("backoff_usage") or {}).items():
            bundle.backoff_usage[level] = (
                bundle.backoff_usage.get(level, 0) + count
            )

    bundle.n_rows_not_tested = len(dataset) - len(pooled_pred)
    if bundle.n_rows_not_tested:
        bundle.notes.append(
            "%d of %d examples appear in no fold's test set under this "
            "grouping and this fold generator (%s), so they are not behind any "
            "number in this row."
            % (bundle.n_rows_not_tested, len(dataset), assignment.source)
        )
    bundle.n_rows_abstained = sum(1 for p in pooled_pred if p is None)
    bundle.n_rows_scored = len(pooled_pred) - bundle.n_rows_abstained
    bundle.confusion = M.confusion_matrix(pooled_true, pooled_pred)
    bundle.per_class = M.per_class_prf(pooled_true, pooled_pred)

    groups = M.index_groups(pooled_group)
    n_groups_pooled = len(groups)

    def _accuracy(indices: Sequence[int]) -> Optional[float]:
        return M.accuracy(
            [pooled_true[i] for i in indices], [pooled_pred[i] for i in indices]
        )

    def _macro_f1(indices: Sequence[int]) -> Optional[float]:
        return M.macro_f1(
            [pooled_true[i] for i in indices], [pooled_pred[i] for i in indices]
        )

    accuracy_value = M.accuracy(pooled_true, pooled_pred)
    bundle.add(
        make_estimate(
            "direction_accuracy",
            grouping,
            accuracy_value,
            bundle.n_rows_scored,
            n_groups_pooled,
            ci=M.bootstrap_group_ci(
                _accuracy, groups, n_boot=n_bootstrap, seed=seed
            ),
            undefined_reason=(
                None
                if accuracy_value is not None
                else "the baseline abstained on every test row"
            ),
            detail={
                "scored_rows": bundle.n_rows_scored,
                "abstained_rows": bundle.n_rows_abstained,
                "coverage": bundle.coverage,
                "note": (
                    "computed over answered rows only; abstentions are neither "
                    "right nor wrong, so read this with coverage"
                ),
            },
        )
    )

    macro_value = M.macro_f1(pooled_true, pooled_pred)
    present = [
        label for label in LABELS if bundle.per_class[label]["support"] > 0
    ]
    bundle.add(
        make_estimate(
            "macro_f1",
            grouping,
            macro_value,
            bundle.n_rows_scored,
            n_groups_pooled,
            ci=M.bootstrap_group_ci(
                _macro_f1, groups, n_boot=n_bootstrap, seed=seed + 1
            ),
            undefined_reason=(
                None if macro_value is not None else "no class had any support"
            ),
            detail={
                "classes_present_in_truth": present,
                "classes_absent_from_truth": [
                    label for label in LABELS if label not in present
                ],
                "note": (
                    "averaged over classes present in the truth; a class the "
                    "baseline never predicts contributes 0"
                ),
            },
        )
    )

    # -- ranking, where applicable ----------------------------------------
    if not any_proba:
        bundle.add(
            make_estimate(
                "roc_auc_ovr_macro",
                grouping,
                None,
                bundle.n_rows_scored,
                n_groups_pooled,
                undefined_reason=(
                    "this baseline emits no probabilities, so there is nothing "
                    "to rank"
                ),
            )
        )
    else:
        auc_value, per_class_auc = M.roc_auc_ovr_macro(pooled_true, pooled_proba)
        distinct_vectors = {
            tuple(round(p.get(label, 0.0), 12) for label in LABELS)
            for p in pooled_proba
            if p
        }
        auc_caveats: List[str] = []
        if auc_value is not None and len(distinct_vectors) <= max(len(folds), 1):
            auc_caveats.append(
                "the baseline emitted only %d distinct probability vector(s) "
                "across %d folds, so this ranking reflects which fold a row "
                "landed in rather than any per-example score"
                % (len(distinct_vectors), len(folds))
            )
        constant = {
            label: M.scores_are_constant(
                [p for p in pooled_proba if p], label
            )
            for label in LABELS
        }
        reason = None
        if auc_value is None:
            if all(constant.values()):
                reason = (
                    "every class score is constant across examples, so the "
                    "baseline ranks nothing; an AUC of 0.5 here would be "
                    "arithmetic, not a measurement"
                )
            else:
                reason = (
                    "no class had both positives and negatives with a varying "
                    "score in the pooled test set"
                )

        def _auc(indices: Sequence[int]) -> Optional[float]:
            return M.roc_auc_ovr_macro(
                [pooled_true[i] for i in indices],
                [pooled_proba[i] for i in indices],
            )[0]

        bundle.add(
            make_estimate(
                "roc_auc_ovr_macro",
                grouping,
                auc_value,
                sum(1 for p in pooled_proba if p),
                n_groups_pooled,
                ci=(
                    M.bootstrap_group_ci(
                        _auc, groups, n_boot=n_bootstrap, seed=seed + 2
                    )
                    if auc_value is not None
                    else None
                ),
                undefined_reason=reason,
                extra_caveats=auc_caveats,
                detail={
                    "per_class_auc": per_class_auc,
                    "constant_score_per_class": constant,
                    "n_distinct_score_vectors": len(distinct_vectors),
                    "n_folds": len(folds),
                },
            )
        )

    # -- calibration -------------------------------------------------------
    if not any_proba:
        bundle.calibration = {
            "grouping": grouping,
            "available": False,
            "reason": (
                "this baseline emits no probabilities. That is a property of "
                "the baseline, not a gap in the evaluation: a hard rule with "
                "no estimated frequency behind it has no confidence to "
                "calibrate, and reporting 1.0 would make it look overconfident "
                "rather than silent."
            ),
        }
        bundle.add(
            make_estimate(
                "brier_score",
                grouping,
                None,
                bundle.n_rows_scored,
                n_groups_pooled,
                undefined_reason="this baseline emits no probabilities",
            )
        )
    else:
        bins = M.reliability_bins(
            pooled_true,
            pooled_proba,
            n_bins=n_reliability_bins,
            min_n=MIN_N_PER_RELIABILITY_BIN,
        )
        ece = M.expected_calibration_error(bins)
        n_with_proba = sum(1 for p in pooled_proba if p)
        bundle.calibration = {
            "grouping": grouping,
            "available": True,
            "n_rows_with_probability": n_with_proba,
            "n_groups": n_groups_pooled,
            "n_bins": n_reliability_bins,
            "min_n_per_bin_for_meaning": MIN_N_PER_RELIABILITY_BIN,
            "expected_calibration_error": ece,
            "bins": [b.as_row() for b in bins],
            "note": (
                "bins are on the probability of the predicted class; a bin "
                "holding fewer than %d rows is marked not meaningful and must "
                "not be read as evidence of (mis)calibration in that range"
                % MIN_N_PER_RELIABILITY_BIN
            ),
        }

        def _brier(indices: Sequence[int]) -> Optional[float]:
            return M.brier_score(
                [pooled_true[i] for i in indices],
                [pooled_proba[i] for i in indices],
            )

        brier_value = M.brier_score(pooled_true, pooled_proba)
        bundle.add(
            make_estimate(
                "brier_score",
                grouping,
                brier_value,
                n_with_proba,
                n_groups_pooled,
                ci=M.bootstrap_group_ci(
                    _brier, groups, n_boot=n_bootstrap, seed=seed + 3
                ),
                undefined_reason=(
                    None if brier_value is not None else "no probability vectors"
                ),
                detail={
                    "range": "0 (perfect) to 2 (confidently wrong); lower is better",
                    "expected_calibration_error": ece,
                },
            )
        )
        bundle.add(
            make_estimate(
                "expected_calibration_error",
                grouping,
                ece,
                n_with_proba,
                n_groups_pooled,
                undefined_reason=(None if ece is not None else "no probability vectors"),
                detail={
                    "definition": (
                        "sample-weighted mean |observed accuracy - mean "
                        "confidence| over non-empty reliability bins"
                    )
                },
            )
        )

    return bundle


# --------------------------------------------------------------------------- #
# paired comparison against a reference baseline
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class ComparisonSlice:
    """A paired accuracy difference on one subset of the pooled test rows.

    Four things are reported together, and the reason they are together is
    that on this corpus they can disagree:

    * ``delta`` -- the paired difference in accuracy on this subset;
    * ``ci`` -- a percentile bootstrap of that difference resampling
      **groups**, which respects SCHEMA.md 5's non-independence;
    * ``mcnemar_p`` -- the exact two-sided row-level test on the discordant
      rows, which does **not** respect it, and is reported anyway because it is
      the test a reader will otherwise run in their head;
    * ``n_rows`` / ``n_groups`` -- the size, printed so that a difference
      measured on a hundred rows in twenty studies cannot be read as though it
      were measured on thousands.

    When the clustered interval excludes zero and the row-level test does not,
    the honest conclusion is *suggestive, not established*. The renderer says
    so in those words rather than quoting whichever of the two is more
    flattering.
    """

    name: str
    description: str
    grouping: str
    n_rows: int
    n_groups: int
    accuracy: Optional[float] = None
    reference_accuracy: Optional[float] = None
    delta: Optional[float] = None
    ci: Optional[Dict[str, Any]] = None
    n_only_baseline_right: int = 0
    n_only_reference_right: int = 0
    n_both_right: int = 0
    n_both_wrong: int = 0
    mcnemar_p: Optional[float] = None
    caveats: Tuple[str, ...] = ()

    def __post_init__(self) -> None:
        _require_grouping(self.grouping)

    @property
    def interval_excludes_zero(self) -> Optional[bool]:
        """``True`` only if the **group-clustered** interval excludes 0."""
        if self.ci is None or self.ci.get("ci_low") is None:
            return None
        return self.ci["ci_low"] > 0.0 or self.ci["ci_high"] < 0.0

    @property
    def row_test_significant(self) -> Optional[bool]:
        """``True`` if the row-level exact McNemar test reaches p < 0.05."""
        if self.mcnemar_p is None:
            return None
        return self.mcnemar_p < 0.05

    @property
    def verdict(self) -> str:
        """One of ``no_difference``, ``suggestive``, ``established``, ``undefined``.

        ``established`` requires **both** tests to agree, and is deliberately
        hard to reach: on 116 rows in 20 studies, a clustered interval that
        excludes zero while the row-level test gives p = 0.16 is one piece of
        evidence, not two.
        """
        if self.delta is None:
            return "undefined"
        interval = self.interval_excludes_zero
        row = self.row_test_significant
        if interval is None and row is None:
            return "undefined"
        if interval and row:
            return "established"
        if interval or row:
            return "suggestive"
        return "no_difference"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "slice": self.name,
            "slice_description": self.description,
            "grouping": _require_grouping(self.grouping),
            "n_rows": self.n_rows,
            "n_groups": self.n_groups,
            "accuracy": self.accuracy,
            "reference_accuracy": self.reference_accuracy,
            "delta": self.delta,
            "ci": self.ci,
            "discordant_rows": {
                "only_baseline_right": self.n_only_baseline_right,
                "only_reference_right": self.n_only_reference_right,
                "both_right": self.n_both_right,
                "both_wrong": self.n_both_wrong,
            },
            "mcnemar_p_two_sided_exact": self.mcnemar_p,
            "clustered_interval_excludes_zero": self.interval_excludes_zero,
            "row_level_test_significant": self.row_test_significant,
            "verdict": self.verdict,
            "caveats": list(self.caveats),
        }

    def render(self) -> str:
        if self.delta is None:
            return "n/a (%d rows)" % self.n_rows
        text = "%+.3f" % self.delta
        if self.ci and self.ci.get("ci_low") is not None:
            text += " [%+.3f, %+.3f]" % (self.ci["ci_low"], self.ci["ci_high"])
        text += " (%.3f vs %.3f, n=%d/%dg)" % (
            self.accuracy,
            self.reference_accuracy,
            self.n_rows,
            self.n_groups,
        )
        return text


#: The subsets every comparison is computed on, as data.
#:
#: ``primary_cell_seen`` is the one that matters for attribution. Where the
#: candidate's primary cell was never populated in training, the candidate and
#: the floor fall back to the same global mode and predict identically **by
#: construction**; such rows can contribute nothing to a difference, and
#: including them both dilutes the effect and imports whatever the label prior
#: happens to be on that subset. Both figures are reported, each labelled, so
#: neither can be quoted as the other.
COMPARISON_SLICES: Tuple[Tuple[str, str], ...] = (
    (
        "all_rows",
        "every pooled out-of-fold row the two baselines both answered",
    ),
    (
        "primary_cell_seen",
        "rows whose primary cell was populated in training -- the only rows "
        "on which the two baselines can differ",
    ),
)


@dataclass(frozen=True)
class Comparison:
    """A **paired** accuracy difference between two baselines, with its grouping.

    Two accuracies printed side by side in a table invite the reader to
    subtract them, and the difference of two independent intervals is not an
    interval for the difference. This is the comparison done properly: both
    baselines are run over the *same* folds, the difference is computed on the
    rows **both** of them answered, and the bootstrap resamples groups of that
    paired set, so the correlation between the two baselines (which is high --
    they share a training set and a backoff tail) is kept rather than thrown
    away.

    ``n_only_baseline_right`` and ``n_only_reference_right`` are the discordant
    counts: the rows where exactly one of the two was right. They are the whole
    evidence for a difference, they need no distributional assumption, and on
    this corpus they are usually small enough to be read directly.

    The same difference is computed on several :class:`ComparisonSlice` s: the
    whole pooled test set, the rows whose primary cell was populated in
    training, and -- because a corpus can hide a real effect inside an average
    -- each of those split by ``perturbation_class``. Two halves of a corpus
    that move in opposite directions cancel in the pooled number, and a
    benchmark that reported only the pooled number would record "no effect"
    where the truth is "an effect in one half and the opposite in the other".
    That is exactly what happens here, which is why the stratification is part
    of the measurement rather than a follow-up.

    Like :class:`Estimate`, it cannot be constructed without a declared
    grouping (SCHEMA.md 5).
    """

    grouping: str
    baseline: str
    reference: str
    status: str = "ok"
    status_detail: Optional[str] = None
    n_rows_unpaired: int = 0
    slices: Tuple[ComparisonSlice, ...] = ()
    caveats: Tuple[str, ...] = ()

    def __post_init__(self) -> None:
        _require_grouping(self.grouping)

    def slice_named(self, name: str) -> Optional[ComparisonSlice]:
        for item in self.slices:
            if item.name == name:
                return item
        return None

    @property
    def overall(self) -> Optional[ComparisonSlice]:
        return self.slice_named("all_rows")

    def to_dict(self) -> Dict[str, Any]:
        grouping = _require_grouping(self.grouping)
        return {
            "grouping": grouping,
            "grouping_description": GROUPINGS[grouping].description,
            "baseline": self.baseline,
            "reference": self.reference,
            "status": self.status,
            "status_detail": self.status_detail,
            "n_rows_unpaired": self.n_rows_unpaired,
            "slices": [item.to_dict() for item in self.slices],
            "caveats": list(self.caveats),
        }


def compare_under_grouping(
    factory: Callable[[], Baseline],
    reference_factory: Callable[[], Baseline],
    dataset: TaskDataset,
    grouping: str,
    *,
    n_folds: int = 5,
    n_bootstrap: int = 1000,
    seed: int = 20240101,
    fold_provider: Optional[FoldProvider] = None,
) -> Comparison:
    """Paired out-of-fold accuracy difference, ``factory`` minus ``reference``."""
    grouping = _require_grouping(grouping)
    probe = factory()
    reference_probe = reference_factory()
    base = dict(grouping=grouping, baseline=probe.name, reference=reference_probe.name)

    for model in (probe, reference_probe):
        if not model.available:
            return Comparison(
                status="baseline_unavailable",
                status_detail="%s cannot run: %s"
                % (model.name, model.availability.reason),
                **base,
            )
    if dataset.is_empty:
        return Comparison(
            status="empty_dataset",
            status_detail="no labellable examples, so nothing to compare",
            **base,
        )

    group_values = dataset.group_values(grouping)
    distinct = {v for v in group_values if v is not None}
    if len(distinct) < MIN_GROUPS_FOR_CV:
        return Comparison(
            status="insufficient_groups",
            status_detail="%d distinct %s group(s)" % (len(distinct), grouping),
            **base,
        )

    provider = fold_provider or resolve_fold_provider()
    folds = provider.folds_for(dataset, grouping, n_folds).folds
    if not folds:
        return Comparison(
            status="insufficient_groups",
            status_detail="the fold generator produced no usable folds",
            **base,
        )

    truth: List[str] = []
    mine: List[Optional[str]] = []
    theirs: List[Optional[str]] = []
    pooled_group: List[Optional[str]] = []
    pooled_class: List[str] = []
    primary_seen: List[Optional[bool]] = []
    primary_level = probe.primary_level
    for train_idx, test_idx in folds:
        train = dataset.subset(train_idx)
        test = dataset.subset(test_idx)
        features = [e.features for e in test]
        truth.extend(e.label for e in test)
        model = factory().fit(train)
        mine.extend(model.predict(features))
        theirs.extend(reference_factory().fit(train).predict(features))
        pooled_group.extend(group_values[i] for i in test_idx)
        pooled_class.extend(
            (e.features.perturbation_class or NULL_FEATURE) for e in test
        )
        for item in features:
            level = model.prediction_level(item)
            primary_seen.append(
                None if primary_level is None else (level == primary_level)
            )

    paired = [
        i for i in range(len(truth)) if mine[i] is not None and theirs[i] is not None
    ]
    caveats: List[str] = []
    unpaired = len(truth) - len(paired)
    if unpaired:
        caveats.append(
            "%d of %d test rows are excluded because at least one of the two "
            "baselines abstained there; the difference is computed only where "
            "both answered" % (unpaired, len(truth))
        )
    if not paired:
        return Comparison(
            status="no_paired_rows",
            status_detail=(
                "the two baselines never both answered on the same row, so no "
                "paired difference exists"
            ),
            n_rows_unpaired=unpaired,
            caveats=tuple(caveats),
            **base,
        )
    if primary_level is None:
        caveats.append(
            "%s has no primary cell, so the 'primary_cell_seen' slices are not "
            "computed for this comparison" % probe.name
        )

    correct_mine = {i: 1 if mine[i] == truth[i] else 0 for i in paired}
    correct_theirs = {i: 1 if theirs[i] == truth[i] else 0 for i in paired}

    def _slice(name: str, description: str, rows: Sequence[int]):
        if not rows:
            # A requested slice with no rows is reported as a slice with no
            # rows, not omitted. "No row of this half had its cell populated in
            # training" is a finding -- it is what SCHEMA.md 5b predicts for a
            # perturbation that exists in one study -- and an absent table row
            # would read as an oversight.
            return ComparisonSlice(
                name=name,
                description=description,
                grouping=grouping,
                n_rows=0,
                n_groups=0,
                caveats=(
                    "no paired test row falls in this slice, so the two "
                    "baselines cannot differ here and no difference exists to "
                    "report",
                ),
            )
        groups = M.index_groups([pooled_group[i] for i in rows])
        row_list = list(rows)

        def _delta(positions: Sequence[int]) -> Optional[float]:
            if not positions:
                return None
            n_ = len(positions)
            return (
                sum(correct_mine[row_list[p]] for p in positions) / n_
                - sum(correct_theirs[row_list[p]] for p in positions) / n_
            )

        # index_groups indexes into the slice's own row list, so the bootstrap
        # resamples whole groups *of this slice*.
        n_rows = len(row_list)
        accuracy = sum(correct_mine[i] for i in row_list) / n_rows
        reference_accuracy = sum(correct_theirs[i] for i in row_list) / n_rows
        only_mine = sum(
            1 for i in row_list if correct_mine[i] and not correct_theirs[i]
        )
        only_theirs = sum(
            1 for i in row_list if correct_theirs[i] and not correct_mine[i]
        )
        # Every slice uses the same seed on purpose. Two slices that happen to
        # hold the same rows -- which is what happens when one half of the
        # corpus contributes every row of a subset -- then report the same
        # interval, rather than two Monte Carlo draws of it that a reader would
        # have to reconcile.
        ci = M.bootstrap_group_ci(_delta, groups, n_boot=n_bootstrap, seed=seed)
        slice_caveats: List[str] = []
        if n_rows < MIN_ROWS_FOR_METRIC:
            slice_caveats.append(
                "n=%d paired rows is below %d: this difference does not carry "
                "meaning" % (n_rows, MIN_ROWS_FOR_METRIC)
            )
        if len(groups) < MIN_GROUPS_FOR_CI:
            slice_caveats.append(
                "only %d independent %s group(s): treat the interval as "
                "indicative" % (len(groups), grouping)
            )
        if ci is not None and ci.get("ci_low") is None and ci.get("note"):
            slice_caveats.append(str(ci["note"]))
        return ComparisonSlice(
            name=name,
            description=description,
            grouping=grouping,
            n_rows=n_rows,
            n_groups=len(groups),
            accuracy=accuracy,
            reference_accuracy=reference_accuracy,
            delta=accuracy - reference_accuracy,
            ci=ci,
            n_only_baseline_right=only_mine,
            n_only_reference_right=only_theirs,
            n_both_right=sum(
                1 for i in row_list if correct_mine[i] and correct_theirs[i]
            ),
            n_both_wrong=sum(
                1 for i in row_list if not correct_mine[i] and not correct_theirs[i]
            ),
            mcnemar_p=M.mcnemar_exact_p(only_mine, only_theirs),
            caveats=tuple(slice_caveats),
        )

    seen_rows = [i for i in paired if primary_seen[i]]
    classes = sorted({pooled_class[i] for i in paired})
    requests: List[Tuple[str, str, List[int]]] = [
        ("all_rows", COMPARISON_SLICES[0][1], list(paired)),
    ]
    if primary_level is not None:
        requests.append(("primary_cell_seen", COMPARISON_SLICES[1][1], seen_rows))
    for cls in classes:
        rows = [i for i in paired if pooled_class[i] == cls]
        requests.append(
            (
                cls,
                "perturbation_class = %s, every paired row" % cls,
                rows,
            )
        )
        if primary_level is not None:
            requests.append(
                (
                    "%s__primary_cell_seen" % cls,
                    "perturbation_class = %s, rows whose primary cell was "
                    "populated in training" % cls,
                    [i for i in rows if primary_seen[i]],
                )
            )

    slices = []
    for name, description, rows in requests:
        item = _slice(name, description, rows)
        if item is not None:
            slices.append(item)

    return Comparison(
        status="ok",
        n_rows_unpaired=unpaired,
        slices=tuple(slices),
        caveats=tuple(caveats),
        **base,
    )


#: Which reference each baseline is compared against, as data.
#:
#: Every baseline is compared against ``global_majority`` -- the floor, and the
#: only comparison that answers "is direction predictable at all". The two
#: schema v0.2.0 baselines are *additionally* compared against
#: ``population_outcome_majority``, which is the same cell without the signed
#: axis and is therefore the comparison that isolates what the axis bought.
COMPARISON_REFERENCES: Tuple[Tuple[str, str], ...] = (
    ("outcome_majority", "global_majority"),
    ("population_outcome_majority", "global_majority"),
    ("perturbation_axis_cell_majority", "global_majority"),
    ("perturbation_axis_cell_majority", "population_outcome_majority"),
    ("perturbation_axis_majority", "global_majority"),
    ("perturbation_axis_majority", "population_outcome_majority"),
    ("perturbation_axis_condition_majority", "global_majority"),
    ("perturbation_axis_condition_majority", "perturbation_axis_majority"),
    ("neurotransmitter_sign", "global_majority"),
)

#: The comparison the report's headline is computed from: the axis cell against
#: the floor. Named here so that the headline cannot quietly be recomputed from
#: a different pair than the one the tables show.
HEADLINE_COMPARISON: Tuple[str, str] = (
    "perturbation_axis_cell_majority",
    "global_majority",
)


# --------------------------------------------------------------------------- #
# the whole benchmark
# --------------------------------------------------------------------------- #

@dataclass
class BenchmarkResult:
    """Every (baseline, grouping) bundle, plus the provenance to read them."""

    dataset_description: Dict[str, Any] = field(default_factory=dict)
    bundles: List[MetricBundle] = field(default_factory=list)
    #: Paired comparisons, one per (candidate, reference, grouping). The
    #: before/after statement in the report is read from these and not from
    #: subtracting two cells of the results table.
    comparisons: List[Comparison] = field(default_factory=list)
    baselines: List[Dict[str, Any]] = field(default_factory=list)
    #: Dataset-level diagnostics from :mod:`syn4rch.models.diagnostics`,
    #: computed before any model is fitted. They answer "was there anything
    #: there to learn", which no row of the results table answers.
    diagnostics: Dict[str, Any] = field(default_factory=dict)
    groupings: List[str] = field(default_factory=list)
    fold_source: str = "unset"
    fold_source_note: str = ""
    settings: Dict[str, Any] = field(default_factory=dict)
    environment: Dict[str, Any] = field(default_factory=dict)
    generated_at: Optional[str] = None

    def bundle(self, baseline: str, grouping: str) -> Optional[MetricBundle]:
        for item in self.bundles:
            if item.baseline == baseline and item.grouping == grouping:
                return item
        return None

    def comparison(
        self, baseline: str, reference: str, grouping: str
    ) -> Optional[Comparison]:
        for item in self.comparisons:
            if (
                item.baseline == baseline
                and item.reference == reference
                and item.grouping == grouping
            ):
                return item
        return None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "generated_at": self.generated_at,
            "dataset": self.dataset_description,
            "settings": self.settings,
            "environment": self.environment,
            "fold_source": self.fold_source,
            "fold_source_note": self.fold_source_note,
            "groupings_reported": list(self.groupings),
            "baselines": list(self.baselines),
            "diagnostics": dict(self.diagnostics),
            "results": [b.to_dict() for b in self.bundles],
            "comparisons": [c.to_dict() for c in self.comparisons],
        }


def run_benchmark(
    dataset: TaskDataset,
    *,
    groupings: Sequence[str] = DEFAULT_GROUPINGS,
    baseline_factories: Optional[Dict[str, Callable[[], Baseline]]] = None,
    n_folds: int = 5,
    n_bootstrap: int = 1000,
    seed: int = 20240101,
    n_reliability_bins: int = 10,
    fold_provider: Optional[FoldProvider] = None,
) -> BenchmarkResult:
    """Cross-validate every baseline under every requested grouping."""
    import datetime

    for grouping in groupings:
        _require_grouping(grouping)
    factories = baseline_factories or all_baseline_factories()
    provider = fold_provider or resolve_fold_provider()

    result = BenchmarkResult(
        dataset_description=dataset.describe(),
        diagnostics=dataset_diagnostics(dataset),
        groupings=list(groupings),
        fold_source=provider.source,
        fold_source_note=provider.note,
        settings={
            "n_folds_requested": n_folds,
            "n_bootstrap": n_bootstrap,
            "bootstrap_resampling_unit": "group",
            "seed": seed,
            "n_reliability_bins": n_reliability_bins,
            "min_rows_for_metric": MIN_ROWS_FOR_METRIC,
            "min_groups_for_ci": MIN_GROUPS_FOR_CI,
            "min_groups_for_cv": MIN_GROUPS_FOR_CV,
            "min_n_per_reliability_bin": MIN_N_PER_RELIABILITY_BIN,
        },
        environment=optional_dependency_report(),
        generated_at=datetime.datetime.now(datetime.timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"
        ),
    )
    result.baselines = [factory().describe() for factory in factories.values()]
    for name, factory in factories.items():
        for grouping in groupings:
            result.bundles.append(
                cross_validate(
                    factory,
                    dataset,
                    grouping,
                    n_folds=n_folds,
                    n_bootstrap=n_bootstrap,
                    seed=seed,
                    n_reliability_bins=n_reliability_bins,
                    fold_provider=provider,
                )
            )
    for candidate, reference in COMPARISON_REFERENCES:
        if candidate not in factories or reference not in factories:
            continue
        for grouping in groupings:
            result.comparisons.append(
                compare_under_grouping(
                    factories[candidate],
                    factories[reference],
                    dataset,
                    grouping,
                    n_folds=n_folds,
                    n_bootstrap=n_bootstrap,
                    seed=seed,
                    fold_provider=provider,
                )
            )
    return result
