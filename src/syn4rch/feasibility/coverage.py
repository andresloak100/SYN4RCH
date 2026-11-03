"""Measured coverage of the curated table: populations, outcomes, and the cells.

Answers feasibility questions 1 and 2 (how many populations carry usable
intervention data, how many outcomes have actually been measured) and builds
the population x outcome matrix the rest of the package reads.

**Why the matrix is dense.** ``build_matrix`` materialises every
``(population, outcome)`` pair, including the ones with zero records, and
reports sparsity. The empty cells are the result: a list of the cells that do
exist answers "what has been studied", while the matrix answers "what has not
been studied", which is what decides whether a prediction task can be posed
at all. Storing the matrix sparsely and printing the non-empty cells would
quietly convert a mostly-unmeasured field into a short table of successes.
Cost, accepted: memory is O(populations x outcomes).

**Why records are excluded, and why nothing is dropped silently.** A record can
be perfectly valid under SCHEMA.md and still not be usable as a modelling
observation -- a ``valence`` of ``unclear`` (SCHEMA.md 3.4 says this is a valid
and expected answer that must not be coerced) states no direction of
perturbation, and a ``direction`` of ``not_reported`` states no outcome.
Such records are excluded from the matrix and **counted by reason** in
:attr:`CoverageMatrix.excluded_counts`, so the exclusions are themselves a
reported statistic rather than an invisible filter. Exclusion here is
exclusion from a feasibility count; it says nothing about admission to
``data/curated/``, which is :mod:`syn4rch.schema.validate`'s decision.

``no_change`` is **not** an exclusion reason. It is a measured outcome like any
other (SCHEMA.md 3.5: dropping null results is forbidden).

Grouping keys are reused from :mod:`syn4rch.qc.keys` so that "same population",
"same study" and "same lab" mean here exactly what they mean in QC. Lab
identity is ``last_author`` + ``last_author_affiliation`` (SCHEMA.md 5.2).

Python 3.10 compatible, standard library only.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from ..qc.keys import lab_key, norm, population_key, study_key
from ..schema.enums import ASSERTIVE_DIRECTIONS, ReviewStatus, Valence

__all__ = [
    "EXCLUSION_CODES",
    "UNRESOLVED_POPULATION_PREFIX",
    "Cell",
    "CoverageMatrix",
    "build_matrix",
    "exclusion_note",
    "exclusion_reasons",
    "has_measured_outcome",
    "has_population_identity",
    "has_usable_intervention",
    "is_usable",
    "outcome_key",
    "partition_records",
]

#: Prefix :func:`syn4rch.qc.keys.population_key` gives an unresolved
#: population (docs/DECISIONS.md D-010). A key equal to the bare prefix means
#: the record names no population at all.
UNRESOLVED_POPULATION_PREFIX = "raw:"

#: Every reason a record can be held out of the feasibility matrix, with the
#: one-line explanation printed in the report. Order is the order checked.
EXCLUSION_CODES: Tuple[Tuple[str, str], ...] = (
    (
        "REJECTED_IN_REVIEW",
        "review_status is 'rejected': a curator ruled it out as evidence",
    ),
    (
        "NO_POPULATION_IDENTITY",
        "neither population_id nor population_raw names a population",
    ),
    (
        "NO_OUTCOME_ID",
        "outcome_id is absent, so the readout cannot be grouped with any other",
    ),
    (
        "NO_PERTURBATION_METHOD",
        "method is absent: the intervention is unspecified",
    ),
    (
        "VALENCE_MISSING",
        "valence is absent: the direction of the intervention is unspecified",
    ),
    (
        "VALENCE_UNCLEAR",
        "valence is 'unclear' -- a valid schema answer (SCHEMA.md 3.4) that "
        "states no direction of perturbation, so it cannot serve as a model "
        "input; kept in the dataset, counted here, never coerced",
    ),
    (
        "DIRECTION_MISSING",
        "direction is absent",
    ),
    (
        "DIRECTION_NOT_REPORTED",
        "direction is 'not_reported': the readout was named but no result "
        "was stated ('no_change' is a result and is NOT excluded)",
    ),
)

_EXCLUSION_NOTES = dict(EXCLUSION_CODES)


def outcome_key(record) -> str:
    """Grouping identity of the readout: the normalised ``outcome_id``.

    SCHEMA.md 3.5 makes ``outcome_id`` required, so unlike populations there is
    no ``raw:`` fallback. A record without one cannot be grouped with any other
    record and is excluded with ``NO_OUTCOME_ID``.
    """
    return norm(record.outcome_id)


def has_population_identity(record) -> bool:
    """True if the record names a population at all (resolved or raw)."""
    key = population_key(record)
    return bool(key) and key != UNRESOLVED_POPULATION_PREFIX


def has_usable_intervention(record) -> bool:
    """True if the perturbation is specified well enough to be a model input.

    Requires a population identity, a ``method``, and a ``valence`` that is
    neither absent nor ``unclear``. This is the predicate behind feasibility
    question 1.
    """
    return (
        has_population_identity(record)
        and record.method is not None
        and record.valence is not None
        and record.valence is not Valence.unclear
    )


def has_measured_outcome(record) -> bool:
    """True if a readout was identified and a result was actually stated.

    ``direction`` must be one of :data:`syn4rch.schema.enums.ASSERTIVE_DIRECTIONS`
    -- which includes ``no_change`` and ``non_monotonic`` and excludes only
    ``not_reported``. This is the predicate behind feasibility question 2.
    """
    return bool(outcome_key(record)) and record.direction in ASSERTIVE_DIRECTIONS


def exclusion_reasons(record) -> List[str]:
    """Codes from :data:`EXCLUSION_CODES` that hold this record out, in order.

    An empty list means the record counts as one usable observation.
    """
    reasons: List[str] = []
    if record.review_status is ReviewStatus.rejected:
        reasons.append("REJECTED_IN_REVIEW")
    if not has_population_identity(record):
        reasons.append("NO_POPULATION_IDENTITY")
    if not outcome_key(record):
        reasons.append("NO_OUTCOME_ID")
    if record.method is None:
        reasons.append("NO_PERTURBATION_METHOD")
    if record.valence is None:
        reasons.append("VALENCE_MISSING")
    elif record.valence is Valence.unclear:
        reasons.append("VALENCE_UNCLEAR")
    if record.direction is None:
        reasons.append("DIRECTION_MISSING")
    elif record.direction not in ASSERTIVE_DIRECTIONS:
        reasons.append("DIRECTION_NOT_REPORTED")
    return reasons


def is_usable(record) -> bool:
    """True if the record contributes one usable observation to the matrix."""
    return not exclusion_reasons(record)


def partition_records(records: Iterable) -> Tuple[List, Dict[str, int], List]:
    """Split into ``(usable, excluded_counts, excluded_records)``.

    ``excluded_counts`` counts *reasons*, so a record failing two predicates is
    counted under both; ``len(excluded_records)`` is the number of records.
    """
    usable: List = []
    excluded: List = []
    counts: Dict[str, int] = {}
    for record in records:
        reasons = exclusion_reasons(record)
        if reasons:
            excluded.append(record)
            for code in reasons:
                counts[code] = counts.get(code, 0) + 1
        else:
            usable.append(record)
    return usable, dict(sorted(counts.items())), excluded


def exclusion_note(code: str) -> str:
    """The printed explanation for an exclusion code."""
    return _EXCLUSION_NOTES.get(code, code)


@dataclass(frozen=True)
class Cell:
    """One ``(population, outcome)`` cell of the coverage matrix.

    A cell with ``n_records == 0`` is a real, deliberately materialised entry:
    it records that nobody has measured this outcome after perturbing this
    population.
    """

    population: str
    outcome: str
    n_records: int = 0
    n_studies: int = 0
    n_labs: int = 0
    n_records_without_lab: int = 0
    studies: Tuple[str, ...] = ()
    labs: Tuple[Tuple[str, str], ...] = ()
    valences: Tuple[str, ...] = ()
    outcome_tiers: Tuple[str, ...] = ()

    @property
    def is_empty(self) -> bool:
        return self.n_records == 0

    @property
    def is_multi_lab(self) -> bool:
        """True if at least two distinct labs contributed (SCHEMA.md 5.2)."""
        return self.n_labs >= 2

    def to_dict(self) -> Dict[str, Any]:
        return {
            "population": self.population,
            "outcome": self.outcome,
            "n_records": self.n_records,
            "n_studies": self.n_studies,
            "n_labs": self.n_labs,
            "n_records_without_lab": self.n_records_without_lab,
            "studies": list(self.studies),
            "labs": [list(lab) for lab in self.labs],
            "valences": list(self.valences),
            "outcome_tiers": list(self.outcome_tiers),
        }


@dataclass
class CoverageMatrix:
    """The dense population x outcome matrix plus its marginals.

    ``cells`` holds one :class:`Cell` for **every** pair in
    ``populations x outcomes``, zeros included.
    """

    populations: Tuple[str, ...] = ()
    outcomes: Tuple[str, ...] = ()
    cells: Dict[Tuple[str, str], Cell] = field(default_factory=dict)
    records_per_population: Dict[str, int] = field(default_factory=dict)
    records_per_outcome: Dict[str, int] = field(default_factory=dict)
    studies_per_population: Dict[str, int] = field(default_factory=dict)
    studies_per_outcome: Dict[str, int] = field(default_factory=dict)
    n_records_read: int = 0
    n_records_usable: int = 0
    n_records_excluded: int = 0
    excluded_counts: Dict[str, int] = field(default_factory=dict)
    n_studies: int = 0
    n_labs: int = 0
    n_records_without_lab: int = 0
    n_populations_unresolved: int = 0
    outcome_tiers: Dict[str, Tuple[str, ...]] = field(default_factory=dict)

    # --- shape ------------------------------------------------------------

    @property
    def n_populations(self) -> int:
        return len(self.populations)

    @property
    def n_outcomes(self) -> int:
        return len(self.outcomes)

    @property
    def n_cells_total(self) -> int:
        return self.n_populations * self.n_outcomes

    @property
    def n_cells_populated(self) -> int:
        return sum(1 for cell in self.cells.values() if not cell.is_empty)

    @property
    def n_cells_empty(self) -> int:
        return self.n_cells_total - self.n_cells_populated

    @property
    def n_cells_multi_lab(self) -> int:
        return sum(1 for cell in self.cells.values() if cell.is_multi_lab)

    @property
    def sparsity(self) -> Optional[float]:
        """Fraction of cells with no records, or ``None`` for an empty matrix.

        ``None`` rather than ``0.0``: a matrix with no rows or no columns has no
        sparsity, and reporting ``0%`` for it would read as perfect coverage.
        """
        if self.n_cells_total == 0:
            return None
        return self.n_cells_empty / self.n_cells_total

    # --- access -----------------------------------------------------------

    def cell(self, population: str, outcome: str) -> Cell:
        """The cell at ``(population, outcome)``; an explicit zero if absent."""
        return self.cells.get(
            (population, outcome), Cell(population=population, outcome=outcome)
        )

    def populated_cells(self) -> List[Cell]:
        """Non-empty cells, most records first then by key (deterministic)."""
        return sorted(
            (c for c in self.cells.values() if not c.is_empty),
            key=lambda c: (-c.n_records, c.population, c.outcome),
        )

    def empty_cells(self) -> List[Cell]:
        return sorted(
            (c for c in self.cells.values() if c.is_empty),
            key=lambda c: (c.population, c.outcome),
        )

    def row(self, population: str) -> List[Cell]:
        return [self.cell(population, outcome) for outcome in self.outcomes]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "populations": list(self.populations),
            "outcomes": list(self.outcomes),
            "n_records_read": self.n_records_read,
            "n_records_usable": self.n_records_usable,
            "n_records_excluded": self.n_records_excluded,
            "excluded_counts": dict(self.excluded_counts),
            "n_studies": self.n_studies,
            "n_labs": self.n_labs,
            "n_records_without_lab": self.n_records_without_lab,
            "n_populations_unresolved": self.n_populations_unresolved,
            "n_cells_total": self.n_cells_total,
            "n_cells_populated": self.n_cells_populated,
            "n_cells_empty": self.n_cells_empty,
            "n_cells_multi_lab": self.n_cells_multi_lab,
            "sparsity": self.sparsity,
            "records_per_population": dict(self.records_per_population),
            "records_per_outcome": dict(self.records_per_outcome),
            "studies_per_population": dict(self.studies_per_population),
            "studies_per_outcome": dict(self.studies_per_outcome),
            # Only the populated cells are serialised: the zeros are implied by
            # populations x outcomes minus these, and writing ~10^4 empty
            # objects to JSON would bury the measured ones.
            "populated_cells": [c.to_dict() for c in self.populated_cells()],
        }


def build_matrix(records: Iterable) -> CoverageMatrix:
    """Count the corpus into a dense population x outcome matrix.

    Rows are populations that have usable intervention data and columns are
    outcomes that have been measured, both as defined by
    :func:`has_usable_intervention` and :func:`has_measured_outcome`; the
    records failing either are counted in ``excluded_counts`` and excluded from
    the axes. Every pair of a row and a column gets a cell, zeros included.
    """
    records = list(records)
    usable, excluded_counts, excluded = partition_records(records)

    matrix = CoverageMatrix(
        n_records_read=len(records),
        n_records_usable=len(usable),
        n_records_excluded=len(excluded),
        excluded_counts=excluded_counts,
    )
    if not usable:
        # No axes, no cells. Deliberately not a 0x0 matrix with a 0% sparsity.
        return matrix

    populations = sorted({population_key(r) for r in usable})
    outcomes = sorted({outcome_key(r) for r in usable})
    matrix.populations = tuple(populations)
    matrix.outcomes = tuple(outcomes)

    grouped: Dict[Tuple[str, str], List] = {}
    per_population: Dict[str, List] = {}
    per_outcome: Dict[str, List] = {}
    tiers: Dict[str, set] = {}
    for record in usable:
        pop = population_key(record)
        out = outcome_key(record)
        grouped.setdefault((pop, out), []).append(record)
        per_population.setdefault(pop, []).append(record)
        per_outcome.setdefault(out, []).append(record)
        if record.outcome_tier is not None:
            tiers.setdefault(out, set()).add(record.outcome_tier.value)

    # The dense product. Explicit zeros are the point (see module docstring).
    cells: Dict[Tuple[str, str], Cell] = {}
    for pop in populations:
        for out in outcomes:
            group = grouped.get((pop, out), [])
            cells[(pop, out)] = _cell_from_group(pop, out, group)
    matrix.cells = cells

    matrix.records_per_population = {
        pop: len(per_population.get(pop, [])) for pop in populations
    }
    matrix.records_per_outcome = {
        out: len(per_outcome.get(out, [])) for out in outcomes
    }
    matrix.studies_per_population = {
        pop: len({study_key(r) for r in per_population.get(pop, []) if r.study_id})
        for pop in populations
    }
    matrix.studies_per_outcome = {
        out: len({study_key(r) for r in per_outcome.get(out, []) if r.study_id})
        for out in outcomes
    }
    matrix.n_studies = len({study_key(r) for r in usable if r.study_id})
    matrix.n_labs = len({lab_key(r) for r in usable if lab_key(r) is not None})
    matrix.n_records_without_lab = sum(1 for r in usable if lab_key(r) is None)
    matrix.n_populations_unresolved = len(
        {p for p in populations if p.startswith(UNRESOLVED_POPULATION_PREFIX)}
    )
    matrix.outcome_tiers = {out: tuple(sorted(tiers.get(out, ()))) for out in outcomes}
    return matrix


def _cell_from_group(population: str, outcome: str, group: Sequence) -> Cell:
    studies = sorted({study_key(r) for r in group if r.study_id})
    labs = sorted({lab_key(r) for r in group if lab_key(r) is not None})
    return Cell(
        population=population,
        outcome=outcome,
        n_records=len(group),
        n_studies=len(studies),
        n_labs=len(labs),
        n_records_without_lab=sum(1 for r in group if lab_key(r) is None),
        studies=tuple(studies),
        labs=tuple(labs),
        valences=tuple(sorted({norm(r.valence) for r in group if r.valence})),
        outcome_tiers=tuple(
            sorted({r.outcome_tier.value for r in group if r.outcome_tier})
        ),
    )
