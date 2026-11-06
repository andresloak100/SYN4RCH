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

What schema v0.2.0 changed here
-------------------------------
Usability used to require a ``valence`` that was neither absent nor
``unclear``. That predicate was correct under v0.1.0, where ``valence`` was
the only signed axis a record had, and it is wrong under v0.2.0.

Measured on the first 362 curated records, it held out 185 of them for
``valence: unclear`` -- and 170 of those 185 are ``conditional_knockout`` or
``conditional_overexpression`` records whose ``valence`` is honestly
``unclear`` (SCHEMA.md 3.4 and 6 both forbid coercing it) but whose
*manipulation* has a perfectly well defined sign: a knockout is ``decrease``
on the amount of the target gene product, an overexpression is ``increase``.
``perturbation_polarity`` is that sign, and :func:`syn4rch.qc.keys.claim_key`
keys on it. So the v0.2.0 predicate reads the polarity, not the valence, and
the 170 molecular records become usable observations.

This is not a loosened threshold. No number in
:mod:`syn4rch.feasibility.gates` moved. It is a corrected predicate: the field
the old one tested stopped being the field that carries the information. The
change admits records rather than excluding them, which is the direction that
needs the most scrutiny, so both predicates stay available --
``key_version="0.1.0"`` reproduces the old count exactly -- and
:func:`compare_usability` reports them side by side so the report can show
what the schema change bought.

**Why records are excluded, and why nothing is dropped silently.** A record can
be perfectly valid under SCHEMA.md and still not be usable as a modelling
observation -- a ``perturbation_polarity`` of ``unclear`` states no direction
of perturbation, and a ``direction`` of ``not_reported`` states no outcome.
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

from ..qc.keys import (
    condition_specificity,
    lab_key,
    norm,
    perturbation_key,
    population_key,
    study_key,
)
from ..schema.enums import ASSERTIVE_DIRECTIONS, ReviewStatus, Valence

__all__ = [
    "EXCLUSION_CODES",
    "KEY_VERSIONS",
    "KEY_VERSION_CURRENT",
    "KEY_VERSION_POOLED",
    "N_CONDITION_SLOTS",
    "UNRESOLVED_POPULATION_PREFIX",
    "Cell",
    "CoverageMatrix",
    "UsabilityComparison",
    "build_matrix",
    "compare_usability",
    "exclusion_note",
    "exclusion_reasons",
    "has_measured_outcome",
    "has_population_identity",
    "has_signed_perturbation",
    "has_usable_intervention",
    "is_usable",
    "outcome_key",
    "partition_records",
    "perturbation_class_of",
    "perturbation_polarity_of",
]

#: Prefix :func:`syn4rch.qc.keys.population_key` gives an unresolved
#: population (docs/DECISIONS.md D-010). A key equal to the bare prefix means
#: the record names no population at all.
UNRESOLVED_POPULATION_PREFIX = "raw:"

#: The schema version whose claim key this module's default predicate matches.
KEY_VERSION_CURRENT = "0.2.0"

#: The previous version's predicate, kept only so the change is falsifiable.
KEY_VERSION_POOLED = "0.1.0"

#: Both, in report order (old first, so a table reads before-then-after).
KEY_VERSIONS: Tuple[str, ...] = (KEY_VERSION_POOLED, KEY_VERSION_CURRENT)

#: Condition slots a fully specified record fills, from
#: :func:`syn4rch.qc.keys.condition_key`.
N_CONDITION_SLOTS = 4

#: Every reason a record can be held out of the feasibility matrix, with the
#: one-line explanation printed in the report. Order is the order checked.
#: Two codes are version-specific and are marked as such in their text: a run
#: reports one of them, never both, and which one it reports is named in the
#: Run table.
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
        "valence is absent: under the v0.1.0 predicate the direction of the "
        "intervention is unspecified (schema v0.1.0 predicate only)",
    ),
    (
        "VALENCE_UNCLEAR",
        "valence is 'unclear' -- a valid schema answer (SCHEMA.md 3.4) that "
        "states no direction of perturbation on the *firing* axis; it is the "
        "correct answer for a molecular perturbation, so under v0.2.0 it is "
        "no longer an exclusion (schema v0.1.0 predicate only)",
    ),
    (
        "NO_PERTURBATION_POLARITY",
        "perturbation_polarity is absent and cannot be derived from method: "
        "the manipulation has no sign on any axis (schema v0.2.0 predicate)",
    ),
    (
        "PERTURBATION_POLARITY_UNCLEAR",
        "perturbation_polarity is 'unclear': the manipulation's sign is not "
        "derivable even along its own class's axis, so it cannot serve as a "
        "model input; kept in the dataset, counted here, never coerced "
        "(schema v0.2.0 predicate)",
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

#: Codes only the v0.1.0 predicate can emit, and codes only v0.2.0 can emit.
#: Published so that a report can say which predicate produced its table.
_V1_ONLY = ("VALENCE_MISSING", "VALENCE_UNCLEAR")
_V2_ONLY = ("NO_PERTURBATION_POLARITY", "PERTURBATION_POLARITY_UNCLEAR")


def outcome_key(record) -> str:
    """Grouping identity of the readout: the normalised ``outcome_id``.

    SCHEMA.md 3.5 makes ``outcome_id`` required, so unlike populations there is
    no ``raw:`` fallback. A record without one cannot be grouped with any other
    record and is excluded with ``NO_OUTCOME_ID``.
    """
    return norm(record.outcome_id)


def perturbation_class_of(record) -> str:
    """The record's perturbation class, via :func:`perturbation_key`.

    Read through the key function rather than the field so that a record
    written before the v0.2.0 migration -- where the class is derived from
    ``method`` -- counts identically to the same record after it.
    """
    return perturbation_key(record)[0]


def perturbation_polarity_of(record) -> str:
    """The record's perturbation polarity, via :func:`perturbation_key`.

    Never read on its own in a report: a polarity without its class means
    different physical things in different classes
    (:class:`syn4rch.schema.enums.PerturbationPolarity`). It is read here only
    to decide whether the manipulation has *a* sign at all.
    """
    return perturbation_key(record)[1]


def has_population_identity(record) -> bool:
    """True if the record names a population at all (resolved or raw)."""
    key = population_key(record)
    return bool(key) and key != UNRESOLVED_POPULATION_PREFIX


def has_signed_perturbation(record) -> bool:
    """True if the manipulation has a sign along its own class's axis.

    The v0.2.0 replacement for "valence is neither absent nor unclear". A
    conditional knockout passes with ``valence: unclear``, because its
    polarity is ``decrease`` on the amount of the target gene product.
    """
    polarity = perturbation_polarity_of(record)
    return bool(polarity) and polarity not in ("unclear", "<null>")


def has_usable_intervention(
    record, key_version: str = KEY_VERSION_CURRENT
) -> bool:
    """True if the perturbation is specified well enough to be a model input.

    Requires a population identity, a ``method``, and a signed manipulation.
    Under ``key_version="0.2.0"`` the sign is ``perturbation_polarity``; under
    ``"0.1.0"`` it is ``valence``, reproducing the pre-migration count exactly
    for the side-by-side comparison. This is the predicate behind feasibility
    question 1.
    """
    _check_key_version(key_version)
    if not has_population_identity(record) or record.method is None:
        return False
    if key_version == KEY_VERSION_POOLED:
        return (
            record.valence is not None and record.valence is not Valence.unclear
        )
    return has_signed_perturbation(record)


def has_measured_outcome(record) -> bool:
    """True if a readout was identified and a result was actually stated.

    ``direction`` must be one of :data:`syn4rch.schema.enums.ASSERTIVE_DIRECTIONS`
    -- which includes ``no_change`` and ``non_monotonic`` and excludes only
    ``not_reported``. This is the predicate behind feasibility question 2.
    """
    return bool(outcome_key(record)) and record.direction in ASSERTIVE_DIRECTIONS


def _check_key_version(key_version: str) -> None:
    if key_version not in KEY_VERSIONS:
        raise KeyError(
            "unknown key_version %r; declared: %s"
            % (key_version, ", ".join(KEY_VERSIONS))
        )


def exclusion_reasons(
    record, key_version: str = KEY_VERSION_CURRENT
) -> List[str]:
    """Codes from :data:`EXCLUSION_CODES` that hold this record out, in order.

    An empty list means the record counts as one usable observation.
    ``key_version`` selects which signed-perturbation predicate is applied;
    see the module docstring.
    """
    _check_key_version(key_version)
    reasons: List[str] = []
    if record.review_status is ReviewStatus.rejected:
        reasons.append("REJECTED_IN_REVIEW")
    if not has_population_identity(record):
        reasons.append("NO_POPULATION_IDENTITY")
    if not outcome_key(record):
        reasons.append("NO_OUTCOME_ID")
    if record.method is None:
        reasons.append("NO_PERTURBATION_METHOD")

    if key_version == KEY_VERSION_POOLED:
        if record.valence is None:
            reasons.append("VALENCE_MISSING")
        elif record.valence is Valence.unclear:
            reasons.append("VALENCE_UNCLEAR")
    else:
        polarity = perturbation_polarity_of(record)
        if not polarity or polarity == "<null>":
            reasons.append("NO_PERTURBATION_POLARITY")
        elif polarity == "unclear":
            reasons.append("PERTURBATION_POLARITY_UNCLEAR")

    if record.direction is None:
        reasons.append("DIRECTION_MISSING")
    elif record.direction not in ASSERTIVE_DIRECTIONS:
        reasons.append("DIRECTION_NOT_REPORTED")
    return reasons


def is_usable(record, key_version: str = KEY_VERSION_CURRENT) -> bool:
    """True if the record contributes one usable observation to the matrix."""
    return not exclusion_reasons(record, key_version=key_version)


def partition_records(
    records: Iterable, key_version: str = KEY_VERSION_CURRENT
) -> Tuple[List, Dict[str, int], List]:
    """Split into ``(usable, excluded_counts, excluded_records)``.

    ``excluded_counts`` counts *reasons*, so a record failing two predicates is
    counted under both; ``len(excluded_records)`` is the number of records.
    """
    usable: List = []
    excluded: List = []
    counts: Dict[str, int] = {}
    for record in records:
        reasons = exclusion_reasons(record, key_version=key_version)
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
class UsabilityComparison:
    """The v0.1.0 and v0.2.0 usability predicates on the same records.

    ``n_admitted_by_v2_only`` is the number the schema change adds. It is
    reported as its own figure, and broken down by perturbation class, because
    a predicate change that admits records is the kind of change that most
    needs to be shown rather than absorbed into a headline count.
    """

    n_records: int = 0
    n_usable_v1: int = 0
    n_usable_v2: int = 0
    n_admitted_by_v2_only: int = 0
    n_admitted_by_v1_only: int = 0
    admitted_by_v2_only_by_class: Dict[str, int] = field(default_factory=dict)
    excluded_counts_v1: Dict[str, int] = field(default_factory=dict)
    excluded_counts_v2: Dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_records": self.n_records,
            "n_usable_v1": self.n_usable_v1,
            "n_usable_v2": self.n_usable_v2,
            "n_admitted_by_v2_only": self.n_admitted_by_v2_only,
            "n_admitted_by_v1_only": self.n_admitted_by_v1_only,
            "admitted_by_v2_only_by_class": dict(self.admitted_by_v2_only_by_class),
            "excluded_counts_v1": dict(self.excluded_counts_v1),
            "excluded_counts_v2": dict(self.excluded_counts_v2),
            "v1_only_codes": list(_V1_ONLY),
            "v2_only_codes": list(_V2_ONLY),
        }


def compare_usability(records: Iterable) -> UsabilityComparison:
    """Apply both usability predicates and report the difference."""
    records = list(records)
    v1, counts1, _ = partition_records(records, KEY_VERSION_POOLED)
    v2, counts2, _ = partition_records(records, KEY_VERSION_CURRENT)
    ids1 = {id(r) for r in v1}
    ids2 = {id(r) for r in v2}
    added = [r for r in records if id(r) in ids2 and id(r) not in ids1]
    by_class: Dict[str, int] = {}
    for record in added:
        key = perturbation_class_of(record)
        by_class[key] = by_class.get(key, 0) + 1
    return UsabilityComparison(
        n_records=len(records),
        n_usable_v1=len(v1),
        n_usable_v2=len(v2),
        n_admitted_by_v2_only=len(added),
        n_admitted_by_v1_only=sum(
            1 for r in records if id(r) in ids1 and id(r) not in ids2
        ),
        admitted_by_v2_only_by_class=dict(sorted(by_class.items())),
        excluded_counts_v1=counts1,
        excluded_counts_v2=counts2,
    )


@dataclass(frozen=True)
class Cell:
    """One ``(population, outcome)`` cell of the coverage matrix.

    A cell with ``n_records == 0`` is a real, deliberately materialised entry:
    it records that nobody has measured this outcome after perturbing this
    population.

    ``min_condition_specificity`` is carried here for the same reason it is
    carried on a replication cell: a cell whose records state none of the four
    conditions is a weaker piece of coverage than one whose records state
    three, and one number for both would say they were the same.
    """

    population: str
    outcome: str
    n_records: int = 0
    n_studies: int = 0
    n_labs: int = 0
    n_records_without_lab: int = 0
    studies: Tuple[str, ...] = ()
    labs: Tuple[str, ...] = ()
    valences: Tuple[str, ...] = ()
    perturbation_classes: Tuple[str, ...] = ()
    perturbation_polarities: Tuple[str, ...] = ()
    min_condition_specificity: Optional[int] = None
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
            "labs": list(self.labs),
            "valences": list(self.valences),
            "perturbation_classes": list(self.perturbation_classes),
            "perturbation_polarities": list(self.perturbation_polarities),
            "min_condition_specificity": self.min_condition_specificity,
            "n_condition_slots": N_CONDITION_SLOTS,
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
    records_per_perturbation_class: Dict[str, int] = field(default_factory=dict)
    condition_specificity_counts: Dict[int, int] = field(default_factory=dict)
    key_version: str = KEY_VERSION_CURRENT
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
    def n_records_condition_unspecified(self) -> int:
        """Usable records stating none of the four conditions."""
        return self.condition_specificity_counts.get(0, 0)

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
            "key_version": self.key_version,
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
            "records_per_perturbation_class": dict(
                self.records_per_perturbation_class
            ),
            "condition_specificity_counts": {
                str(k): v for k, v in sorted(self.condition_specificity_counts.items())
            },
            "n_condition_slots": N_CONDITION_SLOTS,
            # Only the populated cells are serialised: the zeros are implied by
            # populations x outcomes minus these, and writing ~10^4 empty
            # objects to JSON would bury the measured ones.
            "populated_cells": [c.to_dict() for c in self.populated_cells()],
        }


def build_matrix(
    records: Iterable, key_version: str = KEY_VERSION_CURRENT
) -> CoverageMatrix:
    """Count the corpus into a dense population x outcome matrix.

    Rows are populations that have usable intervention data and columns are
    outcomes that have been measured, both as defined by
    :func:`has_usable_intervention` and :func:`has_measured_outcome`; the
    records failing either are counted in ``excluded_counts`` and excluded from
    the axes. Every pair of a row and a column gets a cell, zeros included.

    ``key_version`` selects the usability predicate and is recorded on the
    matrix, so a printed table always says which one produced it.
    """
    _check_key_version(key_version)
    records = list(records)
    usable, excluded_counts, excluded = partition_records(
        records, key_version=key_version
    )

    matrix = CoverageMatrix(
        key_version=key_version,
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
    per_class: Dict[str, int] = {}
    for record in usable:
        key = perturbation_class_of(record)
        per_class[key] = per_class.get(key, 0) + 1
    matrix.records_per_perturbation_class = dict(sorted(per_class.items()))
    spec: Dict[int, int] = {k: 0 for k in range(N_CONDITION_SLOTS + 1)}
    for record in usable:
        spec[condition_specificity(record)] += 1
    matrix.condition_specificity_counts = spec
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
    specificities = [condition_specificity(r) for r in group]
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
        perturbation_classes=tuple(
            sorted({perturbation_class_of(r) for r in group})
        ),
        perturbation_polarities=tuple(
            sorted({perturbation_polarity_of(r) for r in group})
        ),
        # None, not 0, for an empty cell: "nobody measured this" is not the
        # same statement as "the records here state no condition".
        min_condition_specificity=min(specificities) if specificities else None,
        outcome_tiers=tuple(
            sorted({r.outcome_tier.value for r in group if r.outcome_tier})
        ),
    )
