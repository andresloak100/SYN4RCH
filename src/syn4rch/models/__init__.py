"""Baselines and an evaluation harness for direction-of-effect prediction.

    given (population_id, valence, outcome_id, species, acuity)
    predict direction in {increase, decrease, no_change}

This is the only package in the tree permitted to import ``numpy`` or
``scikit-learn``, and it must never *require* them: every import of an optional
library in this package is inside a ``try`` (see :mod:`._compat`), and every
number is computed with the standard library, so ``import syn4rch.models``
succeeds on a bare 3.10 interpreter and a metric computed without numpy is the
same metric.

What this package is for
------------------------
To establish whether direction of effect is predictable **at all** from the
curated table, and to be honest when it is not. The baselines in
:mod:`.baselines` are the bar a sophisticated model has to clear; the harness
in :mod:`.evaluate` is built so that a flattering number is hard to produce by
accident:

* ``docs/SCHEMA.md`` 5 declares a cross-validated number reported without its
  grouping an invalid result, so :class:`.evaluate.Estimate` and
  :class:`.evaluate.MetricBundle` cannot be constructed without a declared
  grouping, and the metric functions in :mod:`.metrics` return bare floats so
  that nothing can be serialised without passing through them.
* Confidence intervals resample **groups**, not rows, because lab-level
  non-independence is what SCHEMA.md 5 names as this corpus's largest risk.
* A cross-validated number is reported beside the fraction of test rows whose
  perturbation identity never appeared in training (``docs/SCHEMA.md`` 5b),
  because grouping does not establish that anything held out was learnable.
* A polarity is never a feature on its own: :mod:`.task` raises if a feature
  tuple or a cell key carries one without its class (``docs/SCHEMA.md`` 3.4b,
  DECISIONS.md D-019).
* Records that cannot be labelled are counted and reported, never dropped.
* Tier B outcomes are excluded from the task (SCHEMA.md 6), and this project
  has no "rejuvenation score".
* A baseline whose inputs do not exist reports itself unavailable rather than
  being made runnable with invented data.

Nothing in this package writes to ``data/``; the CLI writes only its report.
"""

from . import (  # noqa: F401
    _compat,
    baselines,
    diagnostics,
    evaluate,
    metrics,
    report,
    task,
)
from .baselines import (  # noqa: F401
    Availability,
    BackoffMajority,
    Baseline,
    GlobalMajority,
    KGPathSign,
    NeurotransmitterSign,
    OutcomeMajority,
    PerturbationAxisCellMajority,
    PerturbationAxisConditionMajority,
    PerturbationAxisMajority,
    PopulationOutcomeMajority,
    all_baseline_factories,
    baselines_by_feature_generation,
    build_all,
)
from .diagnostics import (  # noqa: F401
    cell_structure,
    dataset_diagnostics,
    molecular_axis_mirroring,
    perturbation_identity_novelty,
)
from .evaluate import (  # noqa: F401
    BenchmarkResult,
    Comparison,
    ComparisonSlice,
    Estimate,
    FoldProvider,
    MetricBundle,
    UnlabelledMetricError,
    compare_under_grouping,
    cross_validate,
    resolve_fold_provider,
    run_benchmark,
)
from .report import headline, render_markdown  # noqa: F401
from .task import (  # noqa: F401
    FEATURE_FIELDS,
    GROUPINGS,
    LABELS,
    PERTURBATION_AXIS,
    TASK,
    UNKNOWN_CONDITION,
    DirectionTask,
    Example,
    ExclusionLedger,
    ExclusionReason,
    Features,
    PolarityWithoutClassError,
    TaskDataset,
    build_dataset,
    load_dataset,
)

__all__ = [n for n in dir() if not n.startswith("_")]
