"""Does the established direction-of-effect result generalise? (D-056)

``docs/DECISIONS.md`` D-055 established one thing: on the circuit half of the
corpus, under leave-one-study-out, on rows whose
``(population, perturbation axis, outcome)`` cell was populated in training, the
axis-cell baseline beats the global mode by +0.208 with a study-clustered
interval that excludes zero *and* a row-level exact McNemar p of 0.004. It ends
by saying the result "rests on a narrow base and the next finding to chase is
whether it generalises".

This module is that chase, and it is deliberately **not** another pass over the
benchmark. ``syn4rch.models.run`` answers "how do the baselines compare"; the
three probes here answer three questions it does not ask, each pre-registered
with a threshold in D-056 before any of this code existed:

**1. Is the signal cross-laboratory?** The same comparison under the ``lab``
grouping. Until D-056 that grouping was ``(last_author, affiliation)``, the
affiliation was free text, and the corpus held 75 studies and 75 "labs" -- so
every leave-one-lab-out number this project published was a leave-one-study-out
number wearing a different label, on a corpus where four laboratories had
published eleven papers between them. The grouping is corrected; this probe
reports the pair, because a key change that moves a result must be shown moving
it (D-050).

**2. Does it need a cell to have been replicated, or merely seen?** For every
test row, how many distinct groups contributed training rows to *its own* cell.
Rows whose cell was reached by one group are a different claim from rows whose
cell was reached by several, and the benchmark pools them. The strata have the
same instrument and the same statistics as the headline slice, so the two are
comparable by construction: :func:`syn4rch.models.evaluate.paired_slice` is
called here, not reimplemented.

**3. Is the verdict rule calibrated at this size?** ``established`` requires a
clustered bootstrap interval excluding zero *and* McNemar p < 0.05. On 154 rows
in 35 clusters that conjunction has never been reached before in this project,
and a rule reached for the first time is a rule worth testing. The permutation
null shuffles labels **within** each study -- preserving every study's own
label distribution and the entire cell structure, destroying only the
association between a cell and a direction *across* studies -- and counts how
often the rule fires on data that cannot support it.

What this module will not do
----------------------------
**It cannot reopen D-048.** That pre-registration named leave-one-study-out,
and its primary was met. Nothing here changes that; a worse number under a
harder split is a statement about the harder split.

**It does not re-estimate anything the benchmark estimates.** Where a number
exists in ``reports/baseline_benchmark.json`` this module reads the same code
path rather than computing a second version of it, because two estimates of one
quantity have to be reconciled and the reconciliation is always the same
discovery: one of them was wrong.

**The permutation null is not a p-value for the real effect.** It is a check on
the *decision rule*, under a null that keeps the clustering. It answers "how
often does this rule fire when it should not", not "is the observed effect
real".

Python 3.10 compatible, standard library only.
"""

import argparse
import dataclasses
import json
import os
import random
import sys
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from .baselines import (
    Baseline,
    GlobalMajority,
    PERTURBATION_AXIS_CELL_BACKOFF,
    PerturbationAxisCellMajority,
)
from .evaluate import (
    ComparisonSlice,
    FoldProvider,
    paired_slice,
    resolve_fold_provider,
)
# The report writer, not a second copy of it: syn4rch.models.run already
# writes atomically and already carries over the replaced file's mode, and a
# reimplementation here would be a second thing to get right.
from .run import _write
from .task import (
    GROUPINGS,
    NULL_FEATURE,
    TaskDataset,
    default_records_path,
    load_dataset,
    repo_root,
)

__all__ = [
    "CANDIDATE",
    "DEPTH_STRATA",
    "GeneralisationResult",
    "HEADLINE_SLICE",
    "PairedRun",
    "PermutationNull",
    "REFERENCE",
    "build_parser",
    "cross_group_probe",
    "default_report_path",
    "depth_probe",
    "main",
    "permutation_probe",
    "render_markdown",
    "run_generalisation",
    "run_paired",
]

#: The candidate and reference of ``evaluate.HEADLINE_COMPARISON``, as
#: factories. Named here so this module cannot be pointed at a different pair
#: than the one D-055 established the result for.
CANDIDATE: Callable[[], Baseline] = PerturbationAxisCellMajority
REFERENCE: Callable[[], Baseline] = GlobalMajority

#: The slice D-055's primary is measured on: circuit rows whose primary cell was
#: populated in training.
HEADLINE_SLICE = "circuit_activity__primary_cell_seen"

#: The perturbation class the primary is restricted to.
CIRCUIT = "circuit_activity"

#: The cell the candidate looks up, from its own backoff chain rather than
#: restated. A copy of this tuple would be a second definition of the cell, and
#: the depth probe would then be measuring the depth of a cell the model does
#: not use.
CELL_FIELDS: Tuple[str, ...] = PERTURBATION_AXIS_CELL_BACKOFF[0][1]

#: Depth strata, as ``(name, lo, hi)`` with ``hi`` inclusive and ``None`` for
#: open. Depth is the number of **distinct groups** of the run's own grouping
#: that contributed training rows to the test row's cell, so under ``lab`` it
#: counts laboratories and under ``study_id`` it counts studies.
DEPTH_STRATA: Tuple[Tuple[str, int, Optional[int]], ...] = (
    ("depth_0", 0, 0),
    ("depth_1", 1, 1),
    ("depth_2_plus", 2, None),
)

#: How often the permutation null may reach ``established`` before the verdict
#: rule is called miscalibrated. D-056 fixed this at 5%, the nominal level the
#: rule's own McNemar threshold claims.
MAX_NULL_ESTABLISHED_RATE = 0.05


@dataclass(frozen=True)
class PairedRun:
    """One pass of the headline comparison over one grouping's folds, row by row.

    Everything a slice can be cut from, kept per pooled out-of-fold row:
    the truth, both predictions, the group, the perturbation class, whether the
    candidate reached its primary cell, and how deep that cell was in training.

    ``cell_depth`` is 0 exactly where ``primary_cell_seen`` is false, and the
    two are computed independently -- one from the baseline's own backoff
    accounting, one from the training set -- so a disagreement between them is
    a bug this class makes visible rather than hides.
    """

    grouping: str
    n_folds: int
    truth: Tuple[str, ...]
    candidate: Tuple[Optional[str], ...]
    reference: Tuple[Optional[str], ...]
    group: Tuple[Optional[str], ...]
    perturbation_class: Tuple[str, ...]
    primary_cell_seen: Tuple[bool, ...]
    cell_depth: Tuple[int, ...]
    paired: Tuple[int, ...]

    @property
    def n_rows(self) -> int:
        return len(self.truth)

    def correctness(self) -> Tuple[Dict[int, int], Dict[int, int]]:
        mine = {
            i: 1 if self.candidate[i] == self.truth[i] else 0 for i in self.paired
        }
        theirs = {
            i: 1 if self.reference[i] == self.truth[i] else 0 for i in self.paired
        }
        return mine, theirs

    def headline_rows(self) -> List[int]:
        """The rows D-055's primary is measured on, under this grouping."""
        return [
            i
            for i in self.paired
            if self.perturbation_class[i] == CIRCUIT and self.primary_cell_seen[i]
        ]

    def depth_disagreements(self) -> int:
        """Rows where the backoff level and the training set disagree.

        Should be 0. Reported rather than asserted, because an assertion here
        would turn a measurement into a crash in the one case where the
        measurement is interesting.
        """
        return sum(
            1
            for i in range(self.n_rows)
            if self.primary_cell_seen[i] != (self.cell_depth[i] > 0)
        )


def _fold_sets(
    dataset: TaskDataset,
    grouping: str,
    n_folds: int,
    fold_provider: Optional[FoldProvider],
) -> Tuple[Tuple[Tuple[int, ...], Tuple[int, ...]], ...]:
    provider = fold_provider or resolve_fold_provider()
    return provider.folds_for(dataset, grouping, n_folds).folds


def run_paired(
    dataset: TaskDataset,
    grouping: str,
    *,
    n_folds: int = 5,
    fold_provider: Optional[FoldProvider] = None,
    folds: Optional[Sequence[Tuple[Sequence[int], Sequence[int]]]] = None,
) -> PairedRun:
    """Fit the headline pair over one grouping's folds and keep every row.

    ``folds`` may be supplied to reuse an already-computed fold set. The folds
    depend only on the grouping values, which no permutation of the labels can
    change, so the permutation probe computes them once rather than 200 times.
    """
    if grouping not in GROUPINGS:
        raise KeyError(
            "undeclared grouping %r; declared: %s"
            % (grouping, ", ".join(sorted(GROUPINGS)))
        )
    group_values = dataset.group_values(grouping)
    if folds is None:
        folds = _fold_sets(dataset, grouping, n_folds, fold_provider)

    truth: List[str] = []
    mine: List[Optional[str]] = []
    theirs: List[Optional[str]] = []
    pooled_group: List[Optional[str]] = []
    pooled_class: List[str] = []
    seen: List[bool] = []
    depth: List[int] = []

    for train_idx, test_idx in folds:
        train = dataset.subset(train_idx)
        groups_per_cell: Dict[Tuple[str, ...], set] = {}
        for position, example in zip(train_idx, train):
            cell = example.features.cell(*CELL_FIELDS)
            value = group_values[position]
            if value is not None:
                groups_per_cell.setdefault(cell, set()).add(value)
            else:
                groups_per_cell.setdefault(cell, set())

        model = CANDIDATE().fit(train)
        floor = REFERENCE().fit(train)
        primary = model.primary_level
        test = dataset.subset(test_idx)
        features = [e.features for e in test]
        truth.extend(e.label for e in test)
        mine.extend(model.predict(features))
        theirs.extend(floor.predict(features))
        pooled_group.extend(group_values[i] for i in test_idx)
        pooled_class.extend(
            (e.features.perturbation_class or NULL_FEATURE) for e in test
        )
        for item in features:
            seen.append(model.prediction_level(item) == primary)
            depth.append(len(groups_per_cell.get(item.cell(*CELL_FIELDS), ())))

    paired = tuple(
        i for i in range(len(truth)) if mine[i] is not None and theirs[i] is not None
    )
    return PairedRun(
        grouping=grouping,
        n_folds=len(folds),
        truth=tuple(truth),
        candidate=tuple(mine),
        reference=tuple(theirs),
        group=tuple(pooled_group),
        perturbation_class=tuple(pooled_class),
        primary_cell_seen=tuple(seen),
        cell_depth=tuple(depth),
        paired=paired,
    )


def _slice_from(
    run: PairedRun,
    name: str,
    description: str,
    rows: Sequence[int],
    *,
    n_bootstrap: int,
    seed: int,
) -> ComparisonSlice:
    mine, theirs = run.correctness()
    return paired_slice(
        name=name,
        description=description,
        grouping=run.grouping,
        rows=rows,
        correct_baseline=mine,
        correct_reference=theirs,
        group_values=run.group,
        n_bootstrap=n_bootstrap,
        seed=seed,
    )


def cross_group_probe(
    dataset: TaskDataset,
    *,
    groupings: Sequence[str] = ("study_id", "lab"),
    n_folds: int = 5,
    n_bootstrap: int = 1000,
    seed: int = 20240101,
    fold_provider: Optional[FoldProvider] = None,
) -> List[Tuple[str, PairedRun, ComparisonSlice]]:
    """The headline slice under each grouping: probe 1 of D-056.

    Returns the pair rather than only the harder number. A grouping change that
    moves a result has to be shown moving it, which means printing the before
    beside the after and not quietly replacing one with the other (D-050).
    """
    out: List[Tuple[str, PairedRun, ComparisonSlice]] = []
    for grouping in groupings:
        run = run_paired(
            dataset, grouping, n_folds=n_folds, fold_provider=fold_provider
        )
        item = _slice_from(
            run,
            HEADLINE_SLICE,
            "circuit rows whose primary cell was populated in training, "
            "leave-one-%s-out" % grouping,
            run.headline_rows(),
            n_bootstrap=n_bootstrap,
            seed=seed,
        )
        out.append((grouping, run, item))
    return out


def depth_probe(
    run: PairedRun,
    *,
    n_bootstrap: int = 1000,
    seed: int = 20240101,
) -> List[ComparisonSlice]:
    """The headline slice cut by training depth of the test row's own cell.

    Probe 2 of D-056. ``depth_0`` is the complement of the headline slice and is
    computed anyway: the candidate and the floor fall back to the same global
    mode there and must be identical **by construction**, so a non-zero delta on
    ``depth_0`` would mean the instrument is not what this project says it is.
    """
    circuit = [i for i in run.paired if run.perturbation_class[i] == CIRCUIT]
    out: List[ComparisonSlice] = []
    for name, low, high in DEPTH_STRATA:
        rows = [
            i
            for i in circuit
            if run.cell_depth[i] >= low
            and (high is None or run.cell_depth[i] <= high)
        ]
        if high is None:
            span = "%d or more" % low
        elif low == high:
            span = "exactly %d" % low
        else:
            span = "%d to %d" % (low, high)
        out.append(
            _slice_from(
                run,
                name,
                "circuit rows whose cell was reached by %s distinct %s "
                "group(s) in training" % (span, run.grouping),
                rows,
                n_bootstrap=n_bootstrap,
                seed=seed,
            )
        )
    return out


def _median(values: Sequence[float]) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2.0


def _quantile(values: Sequence[float], q: float) -> Optional[float]:
    """Nearest-rank quantile. Reuses no interpolation rule it does not need."""
    if not values:
        return None
    ordered = sorted(values)
    index = int(round(q * (len(ordered) - 1)))
    return ordered[max(0, min(index, len(ordered) - 1))]


@dataclass(frozen=True)
class PermutationNull:
    """How often the ``established`` rule fires on labels that cannot support it.

    Probe 3 of D-056. The permutation is **within study**: each study keeps its
    own label multiset and the cell structure is untouched, so what is destroyed
    is only the association between a cell and a direction across studies --
    which is the whole of what the candidate can use.
    """

    grouping: str
    slice_name: str
    n_permutations: int
    observed_delta: Optional[float] = None
    observed_verdict: str = "undefined"
    verdicts: Dict[str, int] = field(default_factory=dict)
    deltas: Tuple[float, ...] = ()
    n_deltas_at_or_above_observed: int = 0
    seconds: Optional[float] = None
    note: str = ""

    @property
    def established_rate(self) -> Optional[float]:
        if not self.n_permutations:
            return None
        return self.verdicts.get("established", 0) / self.n_permutations

    @property
    def calibrated(self) -> Optional[bool]:
        rate = self.established_rate
        if rate is None:
            return None
        return rate < MAX_NULL_ESTABLISHED_RATE

    @property
    def null_centre(self) -> Optional[float]:
        """Median null delta.

        Reported because on this corpus it is **not zero**. Permuting within
        study preserves each study's label multiset, and a cell that one study
        dominates therefore keeps that study's labels -- so cell lookup still
        recovers study-level label structure under the null. A reader who takes
        "delta > 0" as evidence that the cell carries direction information is
        comparing against the wrong origin.
        """
        return _median(self.deltas)

    @property
    def empirical_p(self) -> Optional[float]:
        """``(1 + #{null delta >= observed}) / (1 + n)``.

        The add-one form, so a null distribution that never reached the observed
        delta reports ``1/(n+1)`` rather than 0. A p of exactly 0 from a finite
        permutation set is a claim the set cannot support.
        """
        if not self.n_permutations or self.observed_delta is None:
            return None
        return (1 + self.n_deltas_at_or_above_observed) / (1 + self.n_permutations)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "grouping": self.grouping,
            "slice": self.slice_name,
            "n_permutations": self.n_permutations,
            "observed_delta": self.observed_delta,
            "observed_verdict": self.observed_verdict,
            "verdict_counts": dict(sorted(self.verdicts.items())),
            "established_rate": self.established_rate,
            "max_established_rate": MAX_NULL_ESTABLISHED_RATE,
            "calibrated": self.calibrated,
            "n_deltas_at_or_above_observed": self.n_deltas_at_or_above_observed,
            "empirical_p": self.empirical_p,
            "null_delta_min": min(self.deltas) if self.deltas else None,
            "null_delta_median": self.null_centre,
            "null_delta_mean": (
                sum(self.deltas) / len(self.deltas) if self.deltas else None
            ),
            "null_delta_p95": _quantile(self.deltas, 0.95),
            "null_delta_max": max(self.deltas) if self.deltas else None,
            "seconds": self.seconds,
            "note": self.note,
        }


def _permute_within(dataset: TaskDataset, keys: Sequence[Optional[str]], rng) -> TaskDataset:
    """A copy of ``dataset`` with labels shuffled inside each key's rows."""
    buckets: Dict[Optional[str], List[int]] = {}
    for index, key in enumerate(keys):
        buckets.setdefault(key, []).append(index)
    labels = [e.label for e in dataset.examples]
    for indices in buckets.values():
        pool = [labels[i] for i in indices]
        rng.shuffle(pool)
        for position, index in enumerate(indices):
            labels[index] = pool[position]
    examples = [
        dataclasses.replace(example, label=labels[index])
        for index, example in enumerate(dataset.examples)
    ]
    return dataclasses.replace(dataset, examples=examples)


def permutation_probe(
    dataset: TaskDataset,
    *,
    grouping: str = "study_id",
    n_permutations: int = 200,
    n_folds: int = 5,
    n_bootstrap: int = 1000,
    seed: int = 20240101,
    fold_provider: Optional[FoldProvider] = None,
) -> PermutationNull:
    """Probe 3 of D-056: the null distribution of the verdict rule."""
    folds = _fold_sets(dataset, grouping, n_folds, fold_provider)
    observed_run = run_paired(dataset, grouping, folds=folds)
    observed = _slice_from(
        observed_run,
        HEADLINE_SLICE,
        "observed",
        observed_run.headline_rows(),
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
    study_keys = dataset.group_values("study_id")

    started = time.time()
    verdicts: Dict[str, int] = {}
    deltas: List[float] = []
    at_or_above = 0
    for index in range(n_permutations):
        rng = random.Random(seed + index)
        shuffled = _permute_within(dataset, study_keys, rng)
        run = run_paired(shuffled, grouping, folds=folds)
        # seed + index, not a fixed seed: a bootstrap that drew the same
        # resampling pattern for every permutation would make the null's own
        # Monte Carlo error look smaller than it is.
        item = _slice_from(
            run,
            HEADLINE_SLICE,
            "permuted",
            run.headline_rows(),
            n_bootstrap=n_bootstrap,
            seed=seed + index,
        )
        verdicts[item.verdict] = verdicts.get(item.verdict, 0) + 1
        if item.delta is not None:
            deltas.append(item.delta)
            if observed.delta is not None and item.delta >= observed.delta:
                at_or_above += 1

    return PermutationNull(
        grouping=grouping,
        slice_name=HEADLINE_SLICE,
        n_permutations=n_permutations,
        observed_delta=observed.delta,
        observed_verdict=observed.verdict,
        verdicts=verdicts,
        deltas=tuple(deltas),
        n_deltas_at_or_above_observed=at_or_above,
        seconds=time.time() - started,
        note=(
            "labels permuted within study_id, so every study keeps its own "
            "label distribution and the cell structure is unchanged; what is "
            "destroyed is the association between a cell and a direction "
            "across studies. This is a check on the decision rule, not a "
            "p-value for the observed effect."
        ),
    )


@dataclass
class GeneralisationResult:
    """Everything the three probes measured, plus the provenance to read it."""

    dataset_description: Dict[str, Any] = field(default_factory=dict)
    cross_group: List[Dict[str, Any]] = field(default_factory=list)
    depth: List[Dict[str, Any]] = field(default_factory=list)
    depth_grouping: str = ""
    depth_disagreements: int = 0
    permutation: Optional[Dict[str, Any]] = None
    settings: Dict[str, Any] = field(default_factory=dict)
    generated_at: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "generated_at": self.generated_at,
            "candidate": CANDIDATE().name,
            "reference": REFERENCE().name,
            "slice": HEADLINE_SLICE,
            "cell_fields": list(CELL_FIELDS),
            "dataset": dict(self.dataset_description),
            "settings": dict(self.settings),
            "cross_group": list(self.cross_group),
            "depth": {
                "grouping": self.depth_grouping,
                "strata": list(self.depth),
                "n_rows_where_backoff_and_training_set_disagree": (
                    self.depth_disagreements
                ),
            },
            "permutation": self.permutation,
        }


def run_generalisation(
    dataset: TaskDataset,
    *,
    depth_grouping: str = "lab",
    n_folds: int = 5,
    n_bootstrap: int = 1000,
    seed: int = 20240101,
    n_permutations: int = 200,
    fold_provider: Optional[FoldProvider] = None,
) -> GeneralisationResult:
    """All three probes, once.

    ``depth_grouping`` defaults to ``lab``: the depth question is "how many
    *laboratories* had reached this cell", and asking it under ``study_id``
    would count one laboratory's three papers as three.
    """
    import datetime

    cross = cross_group_probe(
        dataset,
        n_folds=n_folds,
        n_bootstrap=n_bootstrap,
        seed=seed,
        fold_provider=fold_provider,
    )
    runs = {grouping: run for grouping, run, _ in cross}
    depth_run = runs.get(depth_grouping)
    if depth_run is None:
        depth_run = run_paired(
            dataset, depth_grouping, n_folds=n_folds, fold_provider=fold_provider
        )
    strata = depth_probe(depth_run, n_bootstrap=n_bootstrap, seed=seed)

    null = None
    if n_permutations > 0:
        null = permutation_probe(
            dataset,
            grouping="study_id",
            n_permutations=n_permutations,
            n_folds=n_folds,
            n_bootstrap=n_bootstrap,
            seed=seed,
            fold_provider=fold_provider,
        )

    return GeneralisationResult(
        dataset_description={
            "source": dataset.source,
            "n_records_read": dataset.ledger.n_input,
            "n_examples": len(dataset),
            "n_excluded": dataset.ledger.n_excluded,
            "schema_versions": dict(dataset.schema_versions),
            "grouping_coverage": {
                grouping: dataset.grouping_coverage(grouping)
                for grouping in ("study_id", "lab")
            },
        },
        cross_group=[
            {
                "grouping": grouping,
                "n_folds": run.n_folds,
                "n_pooled_rows": run.n_rows,
                "n_paired_rows": len(run.paired),
                "slice": item.to_dict(),
            }
            for grouping, run, item in cross
        ],
        depth=[item.to_dict() for item in strata],
        depth_grouping=depth_grouping,
        depth_disagreements=depth_run.depth_disagreements(),
        permutation=None if null is None else null.to_dict(),
        settings={
            "n_folds_requested": n_folds,
            "n_bootstrap": n_bootstrap,
            "seed": seed,
            "n_permutations": n_permutations,
            "depth_grouping": depth_grouping,
            "max_null_established_rate": MAX_NULL_ESTABLISHED_RATE,
        },
        generated_at=datetime.datetime.now(datetime.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z"),
    )


# --------------------------------------------------------------------------- #
# rendering
# --------------------------------------------------------------------------- #

def _fmt(value: Optional[float], places: int = 3) -> str:
    return "n/a" if value is None else ("%.*f" % (places, value))


def _signed(value: Optional[float]) -> str:
    return "n/a" if value is None else ("%+.3f" % value)


def _interval(item: Dict[str, Any]) -> str:
    ci = item.get("ci") or {}
    if ci.get("ci_low") is None:
        return "n/a"
    return "[%+.3f, %+.3f]" % (ci["ci_low"], ci["ci_high"])


def _p(value: Optional[float]) -> str:
    if value is None:
        return "n/a"
    return "<0.001" if value < 0.001 else "%.3f" % value


def _table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    out = ["| " + " | ".join(header) + " |"]
    out.append("|" + "|".join("---" for _ in header) + "|")
    for row in rows:
        out.append("| " + " | ".join(row) + " |")
    return "\n".join(out)


def render_markdown(result: GeneralisationResult) -> str:
    data = result.to_dict()
    lines: List[str] = []
    lines.append("# SYN4RCH generalisation probes (D-056)")
    lines.append("")
    lines.append(
        "Three questions asked of the result `docs/DECISIONS.md` D-055 "
        "established, each with a threshold fixed in D-056 **before** this "
        "module was written. The comparison is always `%s` against `%s` on "
        "the slice `%s`, the same pair and slice D-055 reports, so nothing "
        "here is a different measurement of a different thing."
        % (data["candidate"], data["reference"], data["slice"])
    )
    lines.append("")
    lines.append(
        "This does **not** reopen D-048. That pre-registration named "
        "leave-one-study-out and its primary was met. A worse number under a "
        "harder split is a statement about the harder split."
    )
    lines.append("")

    lines.append("## Run")
    lines.append("")
    dataset = data["dataset"]
    settings = data["settings"]
    rows = [
        ["generated (UTC)", str(data["generated_at"])],
        ["input", "`%s`" % (dataset.get("source") or "unknown")],
        ["records read", str(dataset.get("n_records_read"))],
        ["labellable examples", str(dataset.get("n_examples"))],
        ["excluded (counted)", str(dataset.get("n_excluded"))],
        ["cell", "`(%s)`" % ", ".join(data["cell_fields"])],
        ["folds requested", str(settings.get("n_folds_requested"))],
        ["bootstrap resamples", str(settings.get("n_bootstrap"))],
        ["permutations", str(settings.get("n_permutations"))],
        ["seed", str(settings.get("seed"))],
    ]
    lines.append(_table(["field", "value"], rows))
    lines.append("")
    coverage = dataset.get("grouping_coverage") or {}
    if coverage:
        lines.append(
            "Group counts, which are the subject of probe 1: "
            + ", ".join(
                "`%s` %d groups (largest %d rows)"
                % (name, item["n_groups"], item["largest_group_rows"])
                for name, item in sorted(coverage.items())
            )
            + "."
        )
        lines.append("")

    lines.append("## Probe 1 - is the signal cross-laboratory?")
    lines.append("")
    lines.append(
        "`lab` was `(last_author, affiliation)` until D-056. The affiliation "
        "is free text, the same department is spelled differently in every "
        "paper, and the corpus therefore held 75 studies and 75 \"labs\" -- so "
        "every leave-one-lab-out number this project published was a "
        "leave-one-study-out number under another name, on a corpus where four "
        "laboratories had published eleven papers between them. Both rows are "
        "printed because a key change that moves a result must be shown moving "
        "it (D-050)."
    )
    lines.append("")
    rows = []
    for item in data["cross_group"]:
        sl = item["slice"]
        rows.append(
            [
                "`%s`" % item["grouping"],
                _signed(sl["delta"]),
                _interval(sl),
                "%s vs %s" % (_fmt(sl["accuracy"]), _fmt(sl["reference_accuracy"])),
                "%d / %d" % (sl["n_rows"], sl["n_groups"]),
                _p(sl["mcnemar_p_two_sided_exact"]),
                "**%s**" % sl["verdict"],
            ]
        )
    lines.append(
        _table(
            [
                "grouping",
                "delta",
                "clustered 95% CI",
                "accuracy",
                "rows / groups",
                "McNemar p",
                "verdict",
            ],
            rows,
        )
    )
    lines.append("")

    lines.append("## Probe 2 - replicated, or merely seen?")
    lines.append("")
    depth = data["depth"]
    lines.append(
        "How many distinct `%s` groups had reached the test row's own cell in "
        "training. `depth_0` is the complement of the headline slice: there the "
        "candidate and the floor both fall back to the global mode, so a "
        "non-zero delta would mean the instrument is not what this project "
        "says it is." % depth["grouping"]
    )
    lines.append("")
    rows = []
    for item in depth["strata"]:
        rows.append(
            [
                "`%s`" % item["slice"],
                _signed(item["delta"]),
                _interval(item),
                "%s vs %s"
                % (_fmt(item["accuracy"]), _fmt(item["reference_accuracy"])),
                "%d / %d" % (item["n_rows"], item["n_groups"]),
                _p(item["mcnemar_p_two_sided_exact"]),
                "**%s**" % item["verdict"],
            ]
        )
    lines.append(
        _table(
            [
                "stratum",
                "delta",
                "clustered 95% CI",
                "accuracy",
                "rows / groups",
                "McNemar p",
                "verdict",
            ],
            rows,
        )
    )
    lines.append("")
    disagreements = depth["n_rows_where_backoff_and_training_set_disagree"]
    lines.append(
        "Rows where the baseline's own backoff level and the training set "
        "disagree about whether the cell was populated: **%d**. The two are "
        "computed independently and must agree; this line is the check."
        % disagreements
    )
    lines.append("")

    null = data.get("permutation")
    lines.append("## Probe 3 - is the `established` rule calibrated here?")
    lines.append("")
    if null is None:
        lines.append("Not run (`--permutations 0`).")
        lines.append("")
    else:
        lines.append(null["note"])
        lines.append("")
        rows = [
            ["permutations", str(null["n_permutations"])],
            ["observed delta", _signed(null["observed_delta"])],
            ["observed verdict", "**%s**" % null["observed_verdict"]],
            [
                "null deltas reaching the observed one",
                "%d of %d" % (
                    null["n_deltas_at_or_above_observed"],
                    null["n_permutations"],
                ),
            ],
            ["empirical p (add-one)", _p(null["empirical_p"])],
            [
                "null delta range",
                "%s to %s"
                % (_signed(null["null_delta_min"]), _signed(null["null_delta_max"])),
            ],
            [
                "null delta median / mean / 95th",
                "%s / %s / %s"
                % (
                    _signed(null["null_delta_median"]),
                    _signed(null["null_delta_mean"]),
                    _signed(null["null_delta_p95"]),
                ),
            ],
            [
                "null reached `established`",
                "%d (%s)"
                % (
                    null["verdict_counts"].get("established", 0),
                    _fmt(null["established_rate"], 4),
                ),
            ],
            [
                "threshold",
                "< %s" % _fmt(null["max_established_rate"], 2),
            ],
            [
                "calibrated",
                {True: "yes", False: "**no**", None: "n/a"}[null["calibrated"]],
            ],
        ]
        lines.append(_table(["field", "value"], rows))
        lines.append("")
        lines.append(
            "Verdicts over the null: "
            + ", ".join(
                "`%s` %d" % (name, count)
                for name, count in sorted(null["verdict_counts"].items())
            )
            + "."
        )
        lines.append("")

    if null is not None and null["null_delta_median"] is not None:
        lines.append(
            "**The null is not centred on zero.** Its median delta is %s. "
            "Permuting within study preserves each study's own label multiset, "
            "so a cell that one study dominates keeps that study's labels and "
            "cell lookup still recovers study-level label structure. The "
            "observed %s has to be read against that origin, not against 0, "
            "and the pre-registered threshold this probe tests is about how "
            "often the *verdict rule* fires -- not about where the null sits."
            % (_signed(null["null_delta_median"]), _signed(null["observed_delta"]))
        )
        lines.append("")

    lines.append("## What this report does not claim")
    lines.append("")
    lines.append(
        "* The permutation null is a check on the **decision rule**, under a "
        "null that keeps the study clustering. It is not a p-value for the "
        "observed effect, and a calibrated rule does not make an effect real."
    )
    lines.append(
        "* Depth is measured in groups that *reached a cell*, not in "
        "independent replications of an experiment. Two papers from one "
        "laboratory are one group here (D-056) and two papers sharing a mouse "
        "cohort are still two, because the schema records no cohort identity "
        "(D-013)."
    )
    lines.append(
        "* Nothing here measures generalisation to an unseen population or to "
        "the molecular half. Both are already answered, negatively and "
        "exactly, in `reports/baseline_benchmark.json`: the candidate's delta "
        "against the floor is 0.000 under leave-one-population-out, and the "
        "molecular subset has no rows whose cell was populated in training."
    )
    lines.append("")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def default_report_path() -> str:
    root = repo_root()
    if root is None:
        return os.path.join("reports", "generalisation.md")
    return os.path.join(root, "reports", "generalisation.md")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python3 -m syn4rch.models.generalise",
        description=(
            "Ask the three D-056 questions of the direction result D-055 "
            "established: is it cross-laboratory, does it need a replicated "
            "cell, and is the 'established' rule calibrated at this size."
        ),
    )
    parser.add_argument(
        "--records",
        dest="records",
        default=None,
        metavar="PATH",
        help="curated JSONL (default: <repo>/data/curated/records.jsonl)",
    )
    parser.add_argument(
        "--out",
        dest="out",
        default=None,
        metavar="PATH",
        help="Markdown report (default: <repo>/reports/generalisation.md)",
    )
    parser.add_argument(
        "--json",
        dest="json_path",
        default=None,
        metavar="PATH",
        help="also write the full result as JSON",
    )
    parser.add_argument(
        "--depth-grouping",
        dest="depth_grouping",
        default="lab",
        metavar="NAME",
        help=(
            "grouping whose groups are counted as a cell's depth "
            "(default: lab -- the depth question is about laboratories)"
        ),
    )
    parser.add_argument("--folds", type=int, default=5, metavar="N")
    parser.add_argument("--bootstrap", type=int, default=1000, metavar="N")
    parser.add_argument(
        "--permutations",
        type=int,
        default=200,
        metavar="N",
        help="permutations for probe 3; 0 skips it",
    )
    parser.add_argument("--seed", type=int, default=20240101, metavar="N")
    parser.add_argument(
        "--any-version",
        action="store_true",
        help="accept records of any schema version",
    )
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.depth_grouping not in GROUPINGS:
        sys.stderr.write(
            "error: undeclared grouping %r. A grouping must be registered in "
            "syn4rch.models.task.GROUPINGS so that its leakage properties are "
            "documented (SCHEMA.md 5). Declared: %s\n"
            % (args.depth_grouping, ", ".join(sorted(GROUPINGS)))
        )
        return 2
    if args.folds < 2:
        sys.stderr.write("error: --folds must be >= 2\n")
        return 2
    if args.bootstrap < 1:
        sys.stderr.write("error: --bootstrap must be >= 1\n")
        return 2
    if args.permutations < 0:
        sys.stderr.write("error: --permutations must be >= 0\n")
        return 2

    dataset = load_dataset(
        args.records or default_records_path(), any_version=args.any_version
    )
    if dataset.is_empty:
        sys.stderr.write(
            "error: no labellable examples, so there is nothing to probe\n"
        )
        return 1

    result = run_generalisation(
        dataset,
        depth_grouping=args.depth_grouping,
        n_folds=args.folds,
        n_bootstrap=args.bootstrap,
        seed=args.seed,
        n_permutations=args.permutations,
    )

    out_path = args.out or default_report_path()
    try:
        _write(out_path, render_markdown(result))
    except OSError as exc:
        sys.stderr.write("error: could not write %s: %s\n" % (out_path, exc))
        return 1
    if args.json_path:
        try:
            _write(
                args.json_path,
                json.dumps(result.to_dict(), indent=2, ensure_ascii=False) + "\n",
            )
        except OSError as exc:
            sys.stderr.write(
                "error: could not write %s: %s\n" % (args.json_path, exc)
            )
            return 1

    if not args.quiet:
        for item in result.cross_group:
            sl = item["slice"]
            sys.stdout.write(
                "%-10s %s  n=%d/%dg  p=%s  %s\n"
                % (
                    item["grouping"],
                    _signed(sl["delta"]),
                    sl["n_rows"],
                    sl["n_groups"],
                    _p(sl["mcnemar_p_two_sided_exact"]),
                    sl["verdict"],
                )
            )
        if result.permutation is not None:
            sys.stdout.write(
                "null established rate: %s (threshold < %s)\n"
                % (
                    _fmt(result.permutation["established_rate"], 4),
                    _fmt(MAX_NULL_ESTABLISHED_RATE, 2),
                )
            )
        sys.stdout.write("wrote %s\n" % out_path)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
