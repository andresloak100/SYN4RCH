"""Quality control over the curated Tier A table.

Principles enforced by this package:

* Hard schema violations (SCHEMA.md 2) block admission to ``data/curated/``.
* Warnings inform a curator; they never trigger an automatic edit.
* Conflicts between studies and possible non-independence are *data*: they
  are recorded on both sides and never resolved, averaged or dropped.
* Nothing here writes to ``data/``; the CLI writes only its own report.
"""

from .comparability import (  # noqa: F401
    REQUIRED_COMPARABILITY_FIELDS,
    comparability_coverage,
    find_missing_comparability_fields,
    missing_comparability_fields,
)
from .conflicts import (  # noqa: F401
    annotate_conflicts,
    directions_conflict,
    find_conflicts,
)
from .duplicates import (  # noqa: F401
    affiliation_tokens,
    find_exact_duplicates,
    find_lab_identity_collisions,
    find_locator_collisions,
    find_possible_non_independence,
    run_duplicate_checks,
)
from .findings import QCFinding, QCSeverity  # noqa: F401
from .report import QCReport, render_text, run_qc  # noqa: F401
from .units import check_units, normalize_unit  # noqa: F401

__all__ = [n for n in dir() if not n.startswith("_")]
