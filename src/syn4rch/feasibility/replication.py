"""Independent replication and direction agreement, per claim.

Answers feasibility questions 3 and 4. The unit is the *claim*
``(population, valence, outcome_id)`` -- :func:`syn4rch.qc.keys.claim_key`,
reused so that "the same claim" means here what it means to conflict
detection. Valence is part of the key because SCHEMA.md 1 forbids aggregating
across valences: pooling an activation result with an inhibition result would
average two opposite-signed interventions into a number describing neither.

Independence is counted at two levels, and they are not interchangeable:

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
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from ..qc.conflicts import directions_conflict
from ..qc.keys import claim_key, lab_key, norm, study_key
from ..schema.enums import ASSERTIVE_DIRECTIONS, Direction

__all__ = [
    "CONFIDENCE_LEVEL",
    "MIXED",
    "MIN_UNITS_FOR_AGREEMENT",
    "ReplicationCell",
    "build_replication_cells",
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


def _agreement(directions: Sequence[str]) -> Tuple[
    Optional[str], bool, Optional[float], Optional[Tuple[float, float]], int, int
]:
    """``(modal, tied, fraction, ci, n_agreeing, n_units)`` for one level.

    ``directions`` is one entry per independent unit; :data:`MIXED` entries
    count in ``n_units`` but can never be modal or agreeing. ``fraction`` and
    ``ci`` are ``None`` when fewer than :data:`MIN_UNITS_FOR_AGREEMENT` units
    state a direction, or when no unit states a single direction.
    """
    n_units = len(directions)
    counts: Dict[str, int] = {}
    for value in directions:
        if value == MIXED:
            continue
        counts[value] = counts.get(value, 0) + 1
    if not counts:
        return None, False, None, None, 0, n_units

    top = max(counts.values())
    winners = sorted(k for k, v in counts.items() if v == top)
    modal = winners[0]  # deterministic; ``tied`` flags that the choice is arbitrary
    tied = len(winners) > 1
    if n_units < MIN_UNITS_FOR_AGREEMENT:
        return modal, tied, None, None, top, n_units
    return modal, tied, top / n_units, wilson_interval(top, n_units), top, n_units


@dataclass(frozen=True)
class ReplicationCell:
    """Replication and direction agreement for one claim.

    ``agreement_studies`` / ``agreement_labs`` are ``None`` when undefined --
    which is the normal state of this corpus, not an error.
    """

    population: str
    valence: str
    outcome: str

    n_records: int = 0
    n_studies: int = 0
    n_labs: int = 0
    n_records_without_lab: int = 0
    studies: Tuple[str, ...] = ()
    labs: Tuple[Tuple[str, str], ...] = ()

    direction_counts: Dict[str, int] = None  # per record
    study_directions: Dict[str, str] = None  # study key -> direction or MIXED
    lab_directions: Dict[str, str] = None  # "author | affiliation" -> direction
    n_records_direction_not_reported: int = 0

    modal_direction: Optional[str] = None
    modal_direction_is_tied: bool = False

    n_studies_stating_direction: int = 0
    n_studies_agreeing: int = 0
    agreement_studies: Optional[float] = None
    agreement_studies_ci: Optional[Tuple[float, float]] = None

    n_labs_stating_direction: int = 0
    n_labs_agreeing: int = 0
    agreement_labs: Optional[float] = None
    agreement_labs_ci: Optional[Tuple[float, float]] = None
    modal_direction_labs: Optional[str] = None
    modal_direction_labs_is_tied: bool = False

    #: QC conflict codes present between distinct studies in this cell, from
    #: :func:`syn4rch.qc.conflicts.directions_conflict` (docs/DECISIONS.md
    #: D-011). Reused rather than re-derived so "disagreement" has one
    #: definition in this repository.
    conflict_codes: Tuple[str, ...] = ()

    outcome_tiers: Tuple[str, ...] = ()

    @property
    def claim(self) -> Tuple[str, str, str]:
        return (self.population, self.valence, self.outcome)

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
    def is_contested(self) -> bool:
        """True if distinct studies state incompatible directions."""
        return bool(self.conflict_codes) or self.modal_direction_is_tied

    def to_dict(self) -> Dict[str, Any]:
        return {
            "population": self.population,
            "valence": self.valence,
            "outcome": self.outcome,
            "n_records": self.n_records,
            "n_studies": self.n_studies,
            "n_labs": self.n_labs,
            "n_records_without_lab": self.n_records_without_lab,
            "studies": list(self.studies),
            "labs": [list(lab) for lab in self.labs],
            "direction_counts": dict(self.direction_counts or {}),
            "study_directions": dict(self.study_directions or {}),
            "lab_directions": dict(self.lab_directions or {}),
            "n_records_direction_not_reported": self.n_records_direction_not_reported,
            "modal_direction": self.modal_direction,
            "modal_direction_is_tied": self.modal_direction_is_tied,
            "n_studies_stating_direction": self.n_studies_stating_direction,
            "n_studies_agreeing": self.n_studies_agreeing,
            "agreement_studies": self.agreement_studies,
            "agreement_studies_ci": (
                list(self.agreement_studies_ci)
                if self.agreement_studies_ci
                else None
            ),
            "n_labs_stating_direction": self.n_labs_stating_direction,
            "n_labs_agreeing": self.n_labs_agreeing,
            "agreement_labs": self.agreement_labs,
            "agreement_labs_ci": (
                list(self.agreement_labs_ci) if self.agreement_labs_ci else None
            ),
            "modal_direction_labs": self.modal_direction_labs,
            "conflict_codes": list(self.conflict_codes),
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


def build_replication_cells(records: Iterable) -> List[ReplicationCell]:
    """One :class:`ReplicationCell` per claim present in ``records``.

    Only claims that occur are built: unlike the coverage matrix, there is no
    dense product here, because replication of a claim nobody made is not a
    quantity. The empty cells are :mod:`syn4rch.feasibility.coverage`'s job.

    Sorted by claim key, so the output is deterministic.
    """
    groups: Dict[Tuple[str, str, str], List] = {}
    for record in records:
        if not norm(record.outcome_id):
            continue  # ungroupable readout; a HARD schema finding, counted by QC
        groups.setdefault(claim_key(record), []).append(record)

    cells: List[ReplicationCell] = []
    for claim, group in sorted(groups.items()):
        cells.append(_build_cell(claim, group))
    return cells


def _build_cell(claim: Tuple[str, str, str], group: Sequence) -> ReplicationCell:
    by_study: Dict[str, List] = {}
    by_lab: Dict[Tuple[str, str], List] = {}
    for record in group:
        by_study.setdefault(study_key(record), []).append(record)
        lab = lab_key(record)
        if lab is not None:
            by_lab.setdefault(lab, []).append(record)

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

    (
        modal,
        tied,
        agreement,
        ci,
        n_agreeing,
        n_units,
    ) = _agreement(list(study_directions.values()))
    (
        lab_modal,
        lab_tied,
        lab_agreement,
        lab_ci,
        n_lab_agreeing,
        n_lab_units,
    ) = _agreement(list(lab_directions.values()))

    return ReplicationCell(
        population=claim[0],
        valence=claim[1],
        outcome=claim[2],
        n_records=len(group),
        n_studies=len(by_study),
        n_labs=len(by_lab),
        n_records_without_lab=sum(1 for r in group if lab_key(r) is None),
        studies=tuple(sorted(by_study)),
        labs=tuple(sorted(by_lab)),
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
        agreement_studies_ci=ci,
        n_labs_stating_direction=n_lab_units,
        n_labs_agreeing=n_lab_agreeing,
        agreement_labs=lab_agreement,
        agreement_labs_ci=lab_ci,
        modal_direction_labs=lab_modal,
        modal_direction_labs_is_tied=lab_tied,
        conflict_codes=_conflict_codes(by_study),
        outcome_tiers=tuple(
            sorted({r.outcome_tier.value for r in group if r.outcome_tier})
        ),
    )


def summarize_replication(cells: Sequence[ReplicationCell]) -> Dict[str, Any]:
    """Corpus-level counts over replication cells. Measured, not estimated.

    ``mean_agreement_where_defined`` is ``None`` when no cell has a defined
    agreement, rather than ``0.0`` or ``1.0``.
    """
    defined = [c for c in cells if c.agreement_is_defined]
    unanimous = [c for c in defined if c.agreement_studies == 1.0]
    return {
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
        "n_cells_with_unknown_lab_records": sum(
            1 for c in cells if c.n_records_without_lab > 0
        ),
        "mean_agreement_where_defined": (
            sum(c.agreement_studies for c in defined) / len(defined)
            if defined
            else None
        ),
        "modal_direction_counts": _modal_direction_counts(cells),
        "confidence_level": CONFIDENCE_LEVEL,
        "interval_method": "wilson_score",
    }


def _modal_direction_counts(cells: Sequence[ReplicationCell]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for cell in cells:
        key = cell.modal_direction or "undetermined"
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))
