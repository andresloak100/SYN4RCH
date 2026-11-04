"""Phase 2: can this question be answered from this corpus?

This package measures the curated Tier A table and answers five questions,
each from counted records only:

1. How many neuronal populations have usable intervention data?
2. How many physiological outcomes have been measured?
3. How many neuron-outcome relationships have independent replication?
4. How often do independent studies agree on effect direction?
5. Which outcomes have enough evidence to support predictive modelling?

Three properties of this package are deliberate and load-bearing:

* **"Insufficient" is a result, not a failure.** :mod:`~syn4rch.feasibility.gates`
  can return ``insufficient`` for every task in the corpus, and the report
  will say so in its first sentence. A feasibility module that cannot reach a
  negative conclusion is not measuring anything.
* **Empty cells are kept.** :mod:`~syn4rch.feasibility.coverage` materialises
  the full population x outcome matrix including zeros. Which relationships
  nobody has measured is the scientifically informative half of the answer,
  and it disappears the moment the matrix is stored sparsely and printed as a
  list of what exists.
* **A cross-validated number cannot be produced without its grouping.**
  SCHEMA.md 5 declares an unlabelled cross-validated number an invalid
  result, so :mod:`~syn4rch.feasibility.splits` makes one structurally hard to
  build: metrics attach to a labelled fold set or raise. The same device
  qualifies the two figures schema v0.2.0 added:
  :class:`~syn4rch.feasibility.replication.Agreement` cannot be built without
  the minimum condition specificity of the cell it came from, and
  :class:`~syn4rch.feasibility.identity.IdentityCoverage` cannot be built
  without a grouping declared in :data:`syn4rch.models.task.GROUPINGS`.
* **The before and after of the schema v0.2.0 claim key are both reported.**
  :func:`~syn4rch.feasibility.replication.compare_claim_keys` and
  :func:`~syn4rch.feasibility.coverage.compare_usability` measure both keys on
  the same records. The difference between them is the measurement of how much
  apparent disagreement in this corpus was condition and perturbation
  heterogeneity; printing only the new numbers would make the change
  unfalsifiable.
* **Grouping is not the whole of leakage control.** SCHEMA.md 5b:
  :mod:`~syn4rch.feasibility.identity` measures the fraction of records whose
  perturbation identity occurs in exactly one group and which are therefore
  unlearnable by construction, whatever a fold-averaged accuracy says.

Nothing in this package writes to ``data/``. It reads
``data/curated/records.jsonl`` and writes only its own report.

Python 3.10 compatible, standard library only (docs/DECISIONS.md D-001).
"""

from .coverage import (  # noqa: F401
    EXCLUSION_CODES,
    KEY_VERSION_CURRENT,
    KEY_VERSION_POOLED,
    Cell,
    CoverageMatrix,
    UsabilityComparison,
    build_matrix,
    compare_usability,
    exclusion_note,
    exclusion_reasons,
    has_measured_outcome,
    has_population_identity,
    has_signed_perturbation,
    has_usable_intervention,
    is_usable,
    outcome_key,
    partition_records,
)
from .identity import (  # noqa: F401
    IDENTITY_FIELDS,
    STATUSES,
    TARGET_CONSTITUTIVE_CLASSES,
    ClassIdentityCoverage,
    ClassStructure,
    Identity,
    IdentityCoverage,
    IdentityReport,
    UnlabelledIdentityError,
    analyse,
    class_structure,
    identity_coverage,
    perturbation_identity,
)
from .gates import (  # noqa: F401
    CHANCE_DIRECTION_RATE,
    MIN_DIRECTION_AGREEMENT_FEASIBLE,
    MIN_DIRECTION_AGREEMENT_MARGINAL,
    MIN_LABS_FEASIBLE,
    MIN_RECORDS_FEASIBLE,
    MIN_STUDIES_FEASIBLE,
    MIN_STUDIES_MARGINAL,
    THRESHOLDS,
    VERDICTS,
    PairVerdict,
    TaskVerdict,
    evaluate_cell,
    evaluate_cells,
    evaluate_records,
    pair_verdicts,
    summarize_verdicts,
)
from .replication import (  # noqa: F401
    CONFIDENCE_LEVEL,
    KEY_POOLED,
    KEY_STRICT,
    MIXED,
    N_CONDITION_SLOTS,
    Agreement,
    KeyCensus,
    KeyComparison,
    ReplicationCell,
    SpecificityMissing,
    build_replication_cells,
    compare_claim_keys,
    summarize_replication,
    wilson_interval,
)
from .report import (  # noqa: F401
    FeasibilityReport,
    build_report,
    render_markdown,
    top_line,
)
from .splits import (  # noqa: F401
    GROUPINGS,
    Fold,
    FoldSet,
    GroupedMetric,
    UnlabelledResultError,
    all_fold_sets,
    attach_metric,
    by_lab,
    by_study,
    leave_one_population_out,
    require_grouping,
)

__all__ = [n for n in dir() if not n.startswith("_")]
