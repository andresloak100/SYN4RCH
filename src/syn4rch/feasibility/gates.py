"""Promotion gates: is a prediction task modellable at all?

Answers feasibility question 5. Each gate is a named module-level constant
with the reason that number was chosen, because a threshold without a stated
reason is a knob, and a knob gets turned until the answer is the one that was
wanted. **No threshold in this module changed at schema v0.2.0.** What changed
is the unit the thresholds are applied to.

**The unit of a task.** A task is one claim --
:func:`syn4rch.qc.keys.claim_key` -- and under v0.2.0 that is

    (population, perturbation_class, perturbation_polarity,
     perturbation_target, outcome_id, sex, diet, age_group, challenge)

not ``(population, valence, outcome_id)``. The perturbation is a model
*input*, not a nuisance to be marginalised: SCHEMA.md 1 forbids aggregating
across perturbations, and pooling "activating these neurons raises glucose"
with "inhibiting them raises glucose" would score a model on a target built
from two opposite interventions. The condition is in the key for the same
reason one step further out: pooling a fasted male cohort with a fed female
one builds a target out of two experiments that were never comparable.

The arithmetic consequence is unavoidable and is reported rather than
softened: a finer key puts fewer records in each cell, so **fewer cells can
clear a gate that counts studies and labs**. A cell that passed on three
studies under the pooled key may be three cells of one study each under the
strict key, and three cells of one study each is the honest description of
what those papers did. :func:`syn4rch.feasibility.replication.compare_claim_keys`
measures exactly how much replication evidence this costs, and the report
prints it next to the verdicts.

:func:`pair_verdicts` then rolls the claims up to the ``(population,
outcome)`` pairs the Phase 2 question is phrased in: a pair is as feasible as
its best claim, and the rollup names which perturbation and which condition
that was.

**Condition specificity qualifies every verdict.** A verdict is never printed
without the minimum condition specificity of the cell behind it. A cell whose
records state none of the four conditions was grouped partly because those
conditions were unknown, not because they were the same, and
``CONDITIONS_ALL_UNKNOWN`` says so on the verdict. It is a **note, not a
gate**: inventing a specificity threshold would be inventing a knob, and this
module's whole discipline is that a threshold needs a stated reason. The note
travels with the verdict and :func:`summarize_verdicts` breaks the verdict
counts down by specificity, so a reader can see how many `feasible` verdicts
rest on cells that matched on "both unknown".

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

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from ..qc.keys import CLAIM_KEY_FIELDS, NULL
from ..schema.enums import OutcomeTier
from .identity import TARGET_CONSTITUTIVE_CLASSES
from .replication import (
    MIN_UNITS_FOR_AGREEMENT,
    MIXED,
    N_CONDITION_SLOTS,
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

#: There is deliberately **no** minimum condition specificity gate. A cell
#: matched on "both unknown" gets the ``CONDITIONS_ALL_UNKNOWN`` note and is
#: counted separately in :func:`summarize_verdicts`; it is not refused by a
#: number this module would have had to invent. Published as a constant so the
#: absence is a decision on the record rather than an omission.
MIN_CONDITION_SPECIFICITY_GATE = None

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
    "MIN_CONDITION_SPECIFICITY_GATE": MIN_CONDITION_SPECIFICITY_GATE,
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
    "CONDITIONS_ALL_UNKNOWN": "at least one record in this cell states none "
    "of the %d conditions (sex, diet, age_group, challenge), so the records "
    "were grouped partly because those conditions were unknown, not because "
    "they were the same; the replication this verdict rests on matched on "
    "'both unknown'" % N_CONDITION_SLOTS,
    "CONDITIONS_PARTIALLY_UNKNOWN": "the best-specified record in this cell "
    "still leaves condition slots unknown; no claim in this corpus has all "
    "%d" % N_CONDITION_SLOTS,
    "TARGET_UNRESOLVED": "this claim's perturbation class is one whose "
    "identity *is* its target (%s) and the target is null, so the claim is "
    "coarser than the manipulation it describes; a null target on a circuit "
    "or structural claim is correct and does not raise this note. See "
    "`syn4rch.feasibility.identity` and reports/perturbation_identity.md"
    % ", ".join(TARGET_CONSTITUTIVE_CLASSES),
}


@dataclass(frozen=True)
class TaskVerdict:
    """The verdict for one claim-keyed task.

    ``claim`` is the key tuple and ``claim_fields`` names its components, so a
    verdict always carries which key produced it. ``reasons`` are the gates
    that decided it; ``notes`` are caveats that did not change the verdict but
    belong next to it -- including the condition-specificity note, which is
    why a verdict can never be printed without it.
    """

    claim: Tuple[str, ...] = ()
    claim_fields: Tuple[str, ...] = CLAIM_KEY_FIELDS
    verdict: str = VERDICT_INSUFFICIENT
    min_condition_specificity: int = 0
    reasons: Tuple[str, ...] = ()
    notes: Tuple[str, ...] = ()
    numbers: Dict[str, Any] = field(default_factory=dict)
    #: Printable perturbation and condition labels, carried rather than
    #: recomputed so that a serialised verdict is self-describing.
    perturbation_label: str = "-"
    condition_label: str = "-"

    @property
    def claim_dict(self) -> Dict[str, str]:
        return dict(zip(self.claim_fields, self.claim))

    @property
    def population(self) -> str:
        return self.claim_dict.get("population", "")

    @property
    def outcome(self) -> str:
        return self.claim_dict.get("outcome_id", "")

    @property
    def valence(self) -> Optional[str]:
        """Only a pooled (v0.1.0) claim carries a valence. ``None`` otherwise."""
        return self.claim_dict.get("valence")

    @property
    def task(self) -> Tuple[str, ...]:
        return self.claim

    @property
    def pair(self) -> Tuple[str, str]:
        return (self.population, self.outcome)

    @property
    def is_feasible(self) -> bool:
        return self.verdict == VERDICT_FEASIBLE

    @property
    def matched_on_unknown_conditions(self) -> bool:
        return self.min_condition_specificity == 0

    @property
    def rank(self) -> int:
        """Position in :data:`VERDICTS`; higher is better. For comparisons."""
        return VERDICTS.index(self.verdict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "claim_fields": list(self.claim_fields),
            "claim": list(self.claim),
            "population": self.population,
            "outcome": self.outcome,
            "valence": self.valence,
            "perturbation": self.perturbation_label,
            "condition": self.condition_label,
            "verdict": self.verdict,
            "min_condition_specificity": self.min_condition_specificity,
            "n_condition_slots": N_CONDITION_SLOTS,
            "matched_on_unknown_conditions": self.matched_on_unknown_conditions,
            "reasons": list(self.reasons),
            "notes": list(self.notes),
            "numbers": dict(self.numbers or {}),
            "thresholds": dict(THRESHOLDS),
        }


@dataclass(frozen=True)
class PairVerdict:
    """A ``(population, outcome)`` rollup over the claims of that pair.

    ``best_perturbation`` and ``best_condition`` name which claim carried the
    rollup, so a pair-level verdict can always be traced back to the one
    experiment type and condition that earned it. The claims are never merged
    into a single pooled task (SCHEMA.md 1).
    """

    population: str
    outcome: str
    verdict: str
    best_perturbation: Optional[str] = None
    best_condition: Optional[str] = None
    best_min_condition_specificity: Optional[int] = None
    verdicts_by_claim: Dict[str, str] = field(default_factory=dict)
    n_claims: int = 0
    n_studies_max: int = 0
    n_labs_max: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "population": self.population,
            "outcome": self.outcome,
            "verdict": self.verdict,
            "best_perturbation": self.best_perturbation,
            "best_condition": self.best_condition,
            "best_min_condition_specificity": self.best_min_condition_specificity,
            "verdicts_by_claim": dict(self.verdicts_by_claim or {}),
            "n_claims": self.n_claims,
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
    agreement = cell.agreement_studies
    numbers: Dict[str, Any] = {
        "n_records": cell.n_records,
        "n_studies": cell.n_studies,
        "n_labs": cell.n_labs,
        "n_records_without_lab": cell.n_records_without_lab,
        "n_studies_stating_direction": cell.n_studies_stating_direction,
        "n_studies_agreeing": cell.n_studies_agreeing,
        "modal_direction": cell.modal_direction,
        "modal_direction_is_tied": cell.modal_direction_is_tied,
        # The bare fractions appear here for gate arithmetic and for JSON; the
        # specificity below them is not optional and travels in the same dict,
        # which is why the two keys are adjacent.
        "agreement_studies": cell.agreement_studies_value,
        "agreement_studies_ci": (
            list(cell.agreement_studies_ci) if cell.agreement_studies_ci else None
        ),
        "agreement_labs": cell.agreement_labs_value,
        "agreement_labs_ci": (
            list(cell.agreement_labs_ci) if cell.agreement_labs_ci else None
        ),
        "agreement_studies_stated": str(agreement) if agreement else None,
        "min_condition_specificity": cell.min_condition_specificity,
        "max_condition_specificity": cell.max_condition_specificity,
        "n_condition_slots": N_CONDITION_SLOTS,
        "matched_on_unknown_conditions": cell.matched_on_unknown_conditions,
        "perturbation": (
            list(cell.perturbation) if cell.perturbation else None
        ),
        "condition": list(cell.condition) if cell.condition else None,
        "conflict_codes": list(cell.conflict_codes),
        "is_directly_opposed": cell.is_directly_opposed,
        "outcome_tiers": list(cell.outcome_tiers),
        "direction_counts": dict(cell.direction_counts or {}),
    }

    notes: List[str] = []
    if cell.matched_on_unknown_conditions:
        notes.append("CONDITIONS_ALL_UNKNOWN")
    elif cell.max_condition_specificity < N_CONDITION_SLOTS:
        notes.append("CONDITIONS_PARTIALLY_UNKNOWN")
    perturbation = cell.perturbation
    # Only for the classes whose identity *is* the target. A null target on a
    # chemogenetic record is correct, not missing: the manipulation acts on the
    # population, which the claim key carries separately. The list lives in
    # :mod:`syn4rch.feasibility.identity` so the two modules cannot disagree
    # about which classes need one.
    if (
        perturbation is not None
        and perturbation[0] in TARGET_CONSTITUTIVE_CLASSES
        and perturbation[2] == NULL
    ):
        notes.append("TARGET_UNRESOLVED")
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

    agreement_value = cell.agreement_studies_value

    blocking: List[str] = []
    if _tier_blocks(cell):
        blocking.append("OUTCOME_IS_TIER_B_ANNOTATION")
    if cell.n_records == 0:
        blocking.append("NO_DATA")
    elif cell.n_records == 1:
        blocking.append("SINGLE_OBSERVATION")
    if cell.n_studies < MIN_STUDIES_MARGINAL:
        blocking.append("NO_INDEPENDENT_REPLICATION")
    if agreement_value is None:
        blocking.append("AGREEMENT_UNDEFINED")
    elif agreement_value < MIN_DIRECTION_AGREEMENT_MARGINAL:
        blocking.append("DIRECTION_CONTESTED")

    if blocking:
        return _verdict(cell, VERDICT_INSUFFICIENT, blocking, notes, numbers)

    shortfalls: List[str] = []
    if cell.n_studies < MIN_STUDIES_FEASIBLE:
        shortfalls.append("STUDIES_BELOW_FEASIBLE")
    if cell.n_labs < MIN_LABS_FEASIBLE:
        shortfalls.append("LABS_BELOW_FEASIBLE")
        if cell.same_lab_only:
            shortfalls.append("SAME_LAB_ONLY")
    if cell.n_records < MIN_RECORDS_FEASIBLE:
        shortfalls.append("RECORDS_BELOW_FEASIBLE")
    if agreement_value < MIN_DIRECTION_AGREEMENT_FEASIBLE:
        shortfalls.append("AGREEMENT_BELOW_FEASIBLE")
    if REQUIRE_AGREEMENT_CI_ABOVE_CHANCE:
        lower = (
            cell.agreement_studies_ci[0] if cell.agreement_studies_ci else None
        )
        if lower is None or lower <= CHANCE_DIRECTION_RATE:
            shortfalls.append("AGREEMENT_CI_INCLUDES_CHANCE")

    if shortfalls:
        return _verdict(cell, VERDICT_MARGINAL, shortfalls, notes, numbers)

    return _verdict(
        cell, VERDICT_FEASIBLE, ["MEETS_ALL_FEASIBLE_GATES"], notes, numbers
    )


def _verdict(
    cell: ReplicationCell,
    verdict: str,
    reasons: Sequence[str],
    notes: Sequence[str],
    numbers: Dict[str, Any],
) -> TaskVerdict:
    return TaskVerdict(
        claim=cell.claim,
        claim_fields=cell.claim_fields,
        verdict=verdict,
        min_condition_specificity=cell.min_condition_specificity,
        reasons=tuple(dict.fromkeys(reasons)),
        notes=tuple(dict.fromkeys(notes)),
        numbers=numbers,
        perturbation_label=cell.perturbation_label,
        condition_label=cell.condition_label,
    )


def evaluate_cells(cells: Iterable[ReplicationCell]) -> List[TaskVerdict]:
    """Verdicts for every replication cell, in cell order."""
    return [evaluate_cell(cell) for cell in cells]


def evaluate_records(records: Iterable, key_variant: str = "strict") -> List[TaskVerdict]:
    """Build replication cells from records and gate every one of them.

    ``key_variant`` is passed through to
    :func:`syn4rch.feasibility.replication.build_replication_cells`. The
    default is the schema v0.2.0 claim key; the pooled variant is available so
    the before/after comparison can be gated as well as counted, never so that
    a pooled verdict can be reported as the current one.
    """
    return evaluate_cells(build_replication_cells(records, key_variant))


def pair_verdicts(verdicts: Sequence[TaskVerdict]) -> List[PairVerdict]:
    """Roll task verdicts up to ``(population, outcome)`` pairs.

    A pair is as feasible as its best claim, and ``best_perturbation`` /
    ``best_condition`` name which one. The rollup is for the Phase 2
    question's phrasing only; the claims are never merged into a single pooled
    task (SCHEMA.md 1).

    Ties are broken by lab count and then by condition specificity, so that
    between two equally feasible claims the rollup names the better-specified
    one rather than whichever sorted first.
    """
    grouped: Dict[Tuple[str, str], List[TaskVerdict]] = {}
    for verdict in verdicts:
        grouped.setdefault(verdict.pair, []).append(verdict)

    out: List[PairVerdict] = []
    for (population, outcome), group in sorted(grouped.items()):
        best = max(
            group,
            key=lambda v: (
                v.rank,
                v.numbers.get("n_labs", 0),
                v.min_condition_specificity,
            ),
        )
        out.append(
            PairVerdict(
                population=population,
                outcome=outcome,
                verdict=best.verdict,
                best_perturbation=best.perturbation_label,
                best_condition=best.condition_label,
                best_min_condition_specificity=best.min_condition_specificity,
                verdicts_by_claim={
                    "%s | %s" % (v.perturbation_label, v.condition_label): v.verdict
                    for v in sorted(
                        group,
                        key=lambda v: (v.perturbation_label, v.condition_label),
                    )
                },
                n_claims=len(group),
                n_studies_max=max(v.numbers.get("n_studies", 0) for v in group),
                n_labs_max=max(v.numbers.get("n_labs", 0) for v in group),
            )
        )
    return out


def summarize_verdicts(verdicts: Sequence[TaskVerdict]) -> Dict[str, Any]:
    """Counts by verdict and by reason, plus which outcomes cleared the gates.

    Verdict counts are additionally broken down by minimum condition
    specificity. A ``feasible`` count that does not say how many of its cells
    matched on "both unknown" is the pooled number this module exists not to
    print.
    """
    by_verdict = {name: 0 for name in VERDICTS}
    by_reason: Dict[str, int] = {}
    by_note: Dict[str, int] = {}
    for verdict in verdicts:
        by_verdict[verdict.verdict] = by_verdict.get(verdict.verdict, 0) + 1
        for reason in verdict.reasons:
            by_reason[reason] = by_reason.get(reason, 0) + 1
        for note in verdict.notes:
            by_note[note] = by_note.get(note, 0) + 1

    feasible = [v for v in verdicts if v.verdict == VERDICT_FEASIBLE]
    marginal = [v for v in verdicts if v.verdict == VERDICT_MARGINAL]
    return {
        "n_tasks": len(verdicts),
        "claim_fields": list(verdicts[0].claim_fields) if verdicts else list(CLAIM_KEY_FIELDS),
        "by_verdict": by_verdict,
        "by_reason": dict(sorted(by_reason.items())),
        "by_note": dict(sorted(by_note.items())),
        "by_verdict_and_min_condition_specificity": {
            name: _specificity_histogram(
                [v for v in verdicts if v.verdict == name]
            )
            for name in VERDICTS
        },
        "n_feasible_matched_on_unknown_conditions": sum(
            1 for v in feasible if v.matched_on_unknown_conditions
        ),
        "n_marginal_matched_on_unknown_conditions": sum(
            1 for v in marginal if v.matched_on_unknown_conditions
        ),
        "n_condition_slots": N_CONDITION_SLOTS,
        "outcomes_feasible": sorted({v.outcome for v in feasible}),
        "outcomes_marginal_or_better": sorted(
            {v.outcome for v in feasible + marginal}
        ),
        "populations_feasible": sorted({v.population for v in feasible}),
        "perturbations_feasible": sorted({v.perturbation_label for v in feasible}),
        "any_feasible": bool(feasible),
        "thresholds": dict(THRESHOLDS),
    }


def _specificity_histogram(verdicts: Sequence[TaskVerdict]) -> Dict[str, int]:
    counts = {str(k): 0 for k in range(N_CONDITION_SLOTS + 1)}
    for verdict in verdicts:
        counts[str(verdict.min_condition_specificity)] += 1
    return counts
