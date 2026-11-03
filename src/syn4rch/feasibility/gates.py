"""Promotion gates: is a prediction task modellable at all?

Answers feasibility question 5. Each gate is a named module-level constant
with the reason that number was chosen, because a threshold without a stated
reason is a knob, and a knob gets turned until the answer is the one that was
wanted.

**The unit of a task.** A task is one claim, ``(population, valence,
outcome_id)``, not ``(population, outcome)``. Valence is a model *input*, not a
nuisance to be marginalised: SCHEMA.md 1 forbids aggregating across valences,
and pooling "activating these neurons raises glucose" with "inhibiting them
raises glucose" would score a model on a target built from two opposite
interventions. :func:`pair_verdicts` then rolls the claims up to the
``(population, outcome)`` pairs the Phase 2 question is phrased in: a pair is
as feasible as its best valence, and the rollup names which one.

**Calibration, and why it is pessimistic.** This corpus is expected to consist
overwhelmingly of cells with a single observation. A gate set that would call
such a cell modellable is useless, so every threshold below is placed where a
single paper, and a single lab repeating itself, cannot pass:

* lab-level independence is required for ``feasible``, never study count
  alone;
* agreement at one observation is undefined, and undefined never passes;
* the Wilson lower bound must clear chance, so unanimity among very few
  studies is ``marginal``, not ``feasible``.

The expected output on a real early corpus is ``insufficient`` for nearly every
task. That is a measurement, not a malfunction.

Python 3.10 compatible, standard library only.
"""

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from ..schema.enums import OutcomeTier
from .replication import (
    MIN_UNITS_FOR_AGREEMENT,
    MIXED,
    ReplicationCell,
    build_replication_cells,
)

__all__ = [
    "CHANCE_DIRECTION_RATE",
    "MIN_LABS_FEASIBLE",
    "MIN_LABS_MARGINAL",
    "MIN_RECORDS_FEASIBLE",
    "MIN_STUDIES_FEASIBLE",
    "MIN_STUDIES_MARGINAL",
    "MIN_DIRECTION_AGREEMENT_FEASIBLE",
    "MIN_DIRECTION_AGREEMENT_MARGINAL",
    "REASON_NOTES",
    "THRESHOLDS",
    "VERDICTS",
    "PairVerdict",
    "TaskVerdict",
    "evaluate_cell",
    "evaluate_cells",
    "evaluate_records",
    "pair_verdicts",
    "summarize_verdicts",
]

VERDICT_FEASIBLE = "feasible"
VERDICT_MARGINAL = "marginal"
VERDICT_INSUFFICIENT = "insufficient"

#: The three possible answers, worst to best. ``insufficient`` is the default
#: and the expected one; it is a conclusion, not an error state.
VERDICTS: Tuple[str, ...] = (
    VERDICT_INSUFFICIENT,
    VERDICT_MARGINAL,
    VERDICT_FEASIBLE,
)

# --- thresholds -------------------------------------------------------------

#: Distinct studies needed before the word "replication" applies at all.
#: Two, because that is the definition: one study is one observation and its
#: agreement with itself is not a measurement. Below this a task is
#: ``insufficient`` whatever else it has.
MIN_STUDIES_MARGINAL = 2

#: Distinct studies needed for ``feasible``. Three, not two, because with two
#: studies a single dissent takes agreement to 0.50 -- the agreement estimate
#: has no room to survive one disagreement, so it carries no information about
#: robustness. Three is the smallest n at which a majority direction can
#: absorb one contradicting study and still be a majority.
MIN_STUDIES_FEASIBLE = 3

#: Distinct labs needed for ``marginal``. One: a single lab's work is still
#: evidence and a marginal verdict says exactly that it is not yet independent.
MIN_LABS_MARGINAL = 1

#: Distinct labs needed for ``feasible``. Two. SCHEMA.md 5 item 2 names one
#: group publishing the same circuit repeatedly as the largest leakage risk in
#: this corpus, so a task whose every record shares a last author and
#: affiliation cannot be cross-validated without leaking: by_lab would yield a
#: single fold. Three papers from one lab are three studies and one piece of
#: independent evidence. This is the binding gate, not the study count.
MIN_LABS_FEASIBLE = 2

#: Records needed for ``feasible``. Three, matching MIN_STUDIES_FEASIBLE
#: because each independent study contributes at least one record; it exists as
#: a separate gate only so that a task cannot pass on three studies whose
#: records were all excluded by a later filter. Deliberately not larger: extra
#: records from the same studies are not extra evidence, so raising this would
#: reward verbose papers rather than replication.
MIN_RECORDS_FEASIBLE = 3

#: The trivial baseline an effect-direction classifier has to beat. One third:
#: the modelling target is the three-way choice {increase, decrease,
#: no_change}, so a model that always guesses the commonest direction is right
#: about a third of the time under an even prior. ``no_change`` is in that set
#: because SCHEMA.md 3.5 keeps null results; dropping it would both bias the
#: corpus and flatter the baseline by making the problem binary.
CHANCE_DIRECTION_RATE = 1.0 / 3.0

#: Direction agreement needed for ``marginal``. 0.60: the lowest value that
#: excludes an even split (two studies disagreeing, 0.50) while admitting a
#: 2-of-3 majority (0.667). An even split is not weak evidence for a direction;
#: it is no evidence for either, and the cell is reported as contested instead.
MIN_DIRECTION_AGREEMENT_MARGINAL = 0.60

#: Direction agreement needed for ``feasible``. 0.70: a 2-of-3 majority
#: (0.667) must not qualify, because one dissent in three is a coin-flip away
#: from a tie, while 3-of-4 (0.75) does. Stated as a round number rather than
#: as a fraction of a particular n so that it does not have to be restated
#: every time the corpus grows.
MIN_DIRECTION_AGREEMENT_FEASIBLE = 0.70

#: For ``feasible``, the Wilson 95% lower bound on agreement must exceed
#: :data:`CHANCE_DIRECTION_RATE`. Point agreement of 1.00 from two studies has
#: a lower bound of about 0.34, which barely clears chance and would fail the
#: study and lab gates anyway; from three unanimous studies it is about 0.44.
#: The gate's job is to stop a point estimate of 1.00 on a handful of papers
#: from being read as a settled fact. The interval is an approximation at these
#: sample sizes (see replication.py) and is used as a width, not as a test.
REQUIRE_AGREEMENT_CI_ABOVE_CHANCE = True

#: Tier B outcomes (aging / lifespan / healthspan) are excluded as prediction
#: *targets*. SCHEMA.md 6: Tier B is an annotation layer, not a training
#: target. A Tier B outcome with abundant evidence is still not a task this
#: project will promote, and the verdict says which rule refused it rather
#: than letting it look like a data shortage.
EXCLUDE_TIER_B_AS_TARGET = True

#: Every threshold in one dict, for the report and for JSON output, so a
#: published verdict carries the gates it was judged against.
THRESHOLDS: Dict[str, Any] = {
    "MIN_STUDIES_MARGINAL": MIN_STUDIES_MARGINAL,
    "MIN_STUDIES_FEASIBLE": MIN_STUDIES_FEASIBLE,
    "MIN_LABS_MARGINAL": MIN_LABS_MARGINAL,
    "MIN_LABS_FEASIBLE": MIN_LABS_FEASIBLE,
    "MIN_RECORDS_FEASIBLE": MIN_RECORDS_FEASIBLE,
    "MIN_DIRECTION_AGREEMENT_MARGINAL": MIN_DIRECTION_AGREEMENT_MARGINAL,
    "MIN_DIRECTION_AGREEMENT_FEASIBLE": MIN_DIRECTION_AGREEMENT_FEASIBLE,
    "CHANCE_DIRECTION_RATE": CHANCE_DIRECTION_RATE,
    "REQUIRE_AGREEMENT_CI_ABOVE_CHANCE": REQUIRE_AGREEMENT_CI_ABOVE_CHANCE,
    "EXCLUDE_TIER_B_AS_TARGET": EXCLUDE_TIER_B_AS_TARGET,
}

#: One line per reason code, printed in the report so a verdict is readable
#: without this source file.
REASON_NOTES: Dict[str, str] = {
    "NO_DATA": "no records for this task",
    "SINGLE_OBSERVATION": "exactly one record: nothing to replicate against",
    "NO_INDEPENDENT_REPLICATION": "fewer than %d distinct studies"
    % MIN_STUDIES_MARGINAL,
    "AGREEMENT_UNDEFINED": "fewer than %d studies state a direction, so "
    "agreement is undefined (not 100%%)" % MIN_UNITS_FOR_AGREEMENT,
    "DIRECTION_CONTESTED": "independent studies do not agree on a direction "
    "(agreement below %.2f); reported as contested, not resolved"
    % MIN_DIRECTION_AGREEMENT_MARGINAL,
    "OUTCOME_IS_TIER_B_ANNOTATION": "outcome is Tier B; SCHEMA.md 6 makes "
    "Tier B an annotation layer, not a training target",
    "STUDIES_BELOW_FEASIBLE": "fewer than %d distinct studies"
    % MIN_STUDIES_FEASIBLE,
    "LABS_BELOW_FEASIBLE": "fewer than %d distinct labs: not independent "
    "evidence (SCHEMA.md 5.2)" % MIN_LABS_FEASIBLE,
    "SAME_LAB_ONLY": "replicated across studies but all from one lab; "
    "discounted, and by_lab cross-validation would have one fold",
    "RECORDS_BELOW_FEASIBLE": "fewer than %d records" % MIN_RECORDS_FEASIBLE,
    "AGREEMENT_BELOW_FEASIBLE": "direction agreement below %.2f"
    % MIN_DIRECTION_AGREEMENT_FEASIBLE,
    "AGREEMENT_CI_INCLUDES_CHANCE": "the Wilson %d%% lower bound on agreement "
    "does not clear the %.2f chance rate" % (95, CHANCE_DIRECTION_RATE),
    "MEETS_ALL_FEASIBLE_GATES": "every gate passed",
    "LAB_IDENTITY_MISSING": "some records carry no lab identity and were "
    "counted in no lab (docs/DECISIONS.md D-009); lab-grouped splits must "
    "exclude them",
    "MODAL_DIRECTION_NON_MONOTONIC": "the modal direction is non_monotonic, "
    "which is not one of the three classes the chance rate assumes",
    "DIRECTION_DISAGREEMENT_PRESENT": "distinct studies state incompatible "
    "directions for this claim (QC conflict codes above); the disagreement is "
    "reported, never resolved, averaged or dropped",
    "WITHIN_STUDY_DISAGREEMENT": "at least one study reports more than one "
    "direction for this claim; it counts in the denominator and agrees with "
    "nothing",
}


@dataclass(frozen=True)
class TaskVerdict:
    """The verdict for one ``(population, valence, outcome)`` task.

    ``reasons`` are the gates that decided it; ``notes`` are caveats that did
    not change the verdict but belong next to it.
    """

    population: str
    valence: str
    outcome: str
    verdict: str
    reasons: Tuple[str, ...] = ()
    notes: Tuple[str, ...] = ()
    numbers: Dict[str, Any] = None

    @property
    def task(self) -> Tuple[str, str, str]:
        return (self.population, self.valence, self.outcome)

    @property
    def pair(self) -> Tuple[str, str]:
        return (self.population, self.outcome)

    @property
    def is_feasible(self) -> bool:
        return self.verdict == VERDICT_FEASIBLE

    @property
    def rank(self) -> int:
        """Position in :data:`VERDICTS`; higher is better. For comparisons."""
        return VERDICTS.index(self.verdict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "population": self.population,
            "valence": self.valence,
            "outcome": self.outcome,
            "verdict": self.verdict,
            "reasons": list(self.reasons),
            "notes": list(self.notes),
            "numbers": dict(self.numbers or {}),
            "thresholds": dict(THRESHOLDS),
        }


@dataclass(frozen=True)
class PairVerdict:
    """A ``(population, outcome)`` rollup over the valences of that pair."""

    population: str
    outcome: str
    verdict: str
    best_valence: Optional[str] = None
    verdicts_by_valence: Dict[str, str] = None
    n_studies_max: int = 0
    n_labs_max: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "population": self.population,
            "outcome": self.outcome,
            "verdict": self.verdict,
            "best_valence": self.best_valence,
            "verdicts_by_valence": dict(self.verdicts_by_valence or {}),
            "n_studies_max": self.n_studies_max,
            "n_labs_max": self.n_labs_max,
        }


def _tier_blocks(cell: ReplicationCell) -> bool:
    """True if every outcome tier seen in this cell is Tier B."""
    if not EXCLUDE_TIER_B_AS_TARGET or not cell.outcome_tiers:
        return False
    return set(cell.outcome_tiers) == {OutcomeTier.B.value}


def evaluate_cell(cell: ReplicationCell) -> TaskVerdict:
    """Apply the gates to one replication cell.

    Order matters: blocking conditions are checked first and short-circuit to
    ``insufficient``, so a task with one observation is never described in
    terms of how much agreement it has.
    """
    numbers: Dict[str, Any] = {
        "n_records": cell.n_records,
        "n_studies": cell.n_studies,
        "n_labs": cell.n_labs,
        "n_records_without_lab": cell.n_records_without_lab,
        "n_studies_stating_direction": cell.n_studies_stating_direction,
        "n_studies_agreeing": cell.n_studies_agreeing,
        "modal_direction": cell.modal_direction,
        "modal_direction_is_tied": cell.modal_direction_is_tied,
        "agreement_studies": cell.agreement_studies,
        "agreement_studies_ci": (
            list(cell.agreement_studies_ci) if cell.agreement_studies_ci else None
        ),
        "agreement_labs": cell.agreement_labs,
        "agreement_labs_ci": (
            list(cell.agreement_labs_ci) if cell.agreement_labs_ci else None
        ),
        "conflict_codes": list(cell.conflict_codes),
        "outcome_tiers": list(cell.outcome_tiers),
        "direction_counts": dict(cell.direction_counts or {}),
    }

    notes: List[str] = []
    if cell.n_records_without_lab:
        notes.append("LAB_IDENTITY_MISSING")
    if cell.modal_direction == "non_monotonic":
        notes.append("MODAL_DIRECTION_NON_MONOTONIC")
    if MIXED in set((cell.study_directions or {}).values()):
        notes.append("WITHIN_STUDY_DISAGREEMENT")
    if cell.same_lab_only:
        notes.append("SAME_LAB_ONLY")
    if cell.conflict_codes:
        notes.append("DIRECTION_DISAGREEMENT_PRESENT")

    blocking: List[str] = []
    if _tier_blocks(cell):
        blocking.append("OUTCOME_IS_TIER_B_ANNOTATION")
    if cell.n_records == 0:
        blocking.append("NO_DATA")
    elif cell.n_records == 1:
        blocking.append("SINGLE_OBSERVATION")
    if cell.n_studies < MIN_STUDIES_MARGINAL:
        blocking.append("NO_INDEPENDENT_REPLICATION")
    if cell.agreement_studies is None:
        blocking.append("AGREEMENT_UNDEFINED")
    elif cell.agreement_studies < MIN_DIRECTION_AGREEMENT_MARGINAL:
        blocking.append("DIRECTION_CONTESTED")

    if blocking:
        return TaskVerdict(
            population=cell.population,
            valence=cell.valence,
            outcome=cell.outcome,
            verdict=VERDICT_INSUFFICIENT,
            reasons=tuple(dict.fromkeys(blocking)),
            notes=tuple(dict.fromkeys(notes)),
            numbers=numbers,
        )

    shortfalls: List[str] = []
    if cell.n_studies < MIN_STUDIES_FEASIBLE:
        shortfalls.append("STUDIES_BELOW_FEASIBLE")
    if cell.n_labs < MIN_LABS_FEASIBLE:
        shortfalls.append("LABS_BELOW_FEASIBLE")
        if cell.same_lab_only:
            shortfalls.append("SAME_LAB_ONLY")
    if cell.n_records < MIN_RECORDS_FEASIBLE:
        shortfalls.append("RECORDS_BELOW_FEASIBLE")
    if cell.agreement_studies < MIN_DIRECTION_AGREEMENT_FEASIBLE:
        shortfalls.append("AGREEMENT_BELOW_FEASIBLE")
    if REQUIRE_AGREEMENT_CI_ABOVE_CHANCE:
        lower = (
            cell.agreement_studies_ci[0] if cell.agreement_studies_ci else None
        )
        if lower is None or lower <= CHANCE_DIRECTION_RATE:
            shortfalls.append("AGREEMENT_CI_INCLUDES_CHANCE")

    if shortfalls:
        return TaskVerdict(
            population=cell.population,
            valence=cell.valence,
            outcome=cell.outcome,
            verdict=VERDICT_MARGINAL,
            reasons=tuple(dict.fromkeys(shortfalls)),
            notes=tuple(dict.fromkeys(notes)),
            numbers=numbers,
        )

    return TaskVerdict(
        population=cell.population,
        valence=cell.valence,
        outcome=cell.outcome,
        verdict=VERDICT_FEASIBLE,
        reasons=("MEETS_ALL_FEASIBLE_GATES",),
        notes=tuple(dict.fromkeys(notes)),
        numbers=numbers,
    )


def evaluate_cells(cells: Iterable[ReplicationCell]) -> List[TaskVerdict]:
    """Verdicts for every replication cell, in cell order."""
    return [evaluate_cell(cell) for cell in cells]


def evaluate_records(records: Iterable) -> List[TaskVerdict]:
    """Build replication cells from records and gate every one of them."""
    return evaluate_cells(build_replication_cells(records))


def pair_verdicts(verdicts: Sequence[TaskVerdict]) -> List[PairVerdict]:
    """Roll task verdicts up to ``(population, outcome)`` pairs.

    A pair is as feasible as its best valence, and ``best_valence`` names
    which one. The rollup is for the Phase 2 question's phrasing only; the
    valences are never merged into a single pooled task (SCHEMA.md 1).
    """
    grouped: Dict[Tuple[str, str], List[TaskVerdict]] = {}
    for verdict in verdicts:
        grouped.setdefault(verdict.pair, []).append(verdict)

    out: List[PairVerdict] = []
    for (population, outcome), group in sorted(grouped.items()):
        best = max(group, key=lambda v: (v.rank, v.numbers.get("n_labs", 0)))
        out.append(
            PairVerdict(
                population=population,
                outcome=outcome,
                verdict=best.verdict,
                best_valence=best.valence,
                verdicts_by_valence={
                    v.valence: v.verdict for v in sorted(group, key=lambda v: v.valence)
                },
                n_studies_max=max(v.numbers.get("n_studies", 0) for v in group),
                n_labs_max=max(v.numbers.get("n_labs", 0) for v in group),
            )
        )
    return out


def summarize_verdicts(verdicts: Sequence[TaskVerdict]) -> Dict[str, Any]:
    """Counts by verdict and by reason, plus which outcomes cleared the gates."""
    by_verdict = {name: 0 for name in VERDICTS}
    by_reason: Dict[str, int] = {}
    for verdict in verdicts:
        by_verdict[verdict.verdict] = by_verdict.get(verdict.verdict, 0) + 1
        for reason in verdict.reasons:
            by_reason[reason] = by_reason.get(reason, 0) + 1

    feasible = [v for v in verdicts if v.verdict == VERDICT_FEASIBLE]
    marginal = [v for v in verdicts if v.verdict == VERDICT_MARGINAL]
    return {
        "n_tasks": len(verdicts),
        "by_verdict": by_verdict,
        "by_reason": dict(sorted(by_reason.items())),
        "outcomes_feasible": sorted({v.outcome for v in feasible}),
        "outcomes_marginal_or_better": sorted(
            {v.outcome for v in feasible + marginal}
        ),
        "populations_feasible": sorted({v.population for v in feasible}),
        "any_feasible": bool(feasible),
        "thresholds": dict(THRESHOLDS),
    }
