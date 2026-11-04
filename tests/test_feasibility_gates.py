"""Promotion gates: feasible / marginal / insufficient, and why.

Schema v0.2.0: the gates are unchanged, the unit they apply to is the
nine-component claim key, and no verdict is printed without the minimum
condition specificity of the cell behind it.
"""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.feasibility.gates import (
    CHANCE_DIRECTION_RATE,
    MIN_CONDITION_SPECIFICITY_GATE,
    MIN_DIRECTION_AGREEMENT_FEASIBLE,
    MIN_LABS_FEASIBLE,
    MIN_STUDIES_FEASIBLE,
    REASON_NOTES,
    THRESHOLDS,
    VERDICTS,
    evaluate_records,
    pair_verdicts,
    summarize_verdicts,
)
from syn4rch.feasibility.replication import KEY_POOLED, N_CONDITION_SLOTS


def verdict_of(records):
    verdicts = evaluate_records(records)
    assert len(verdicts) == 1, "expected one task, got %d" % len(verdicts)
    return verdicts[0]


def third_lab_record(**overrides):
    """A third fabricated study from a third fabricated lab."""
    data = dict(
        study_id="synthetic:test-004",
        doi="synthetic:test-004",
        last_author="Third, D.",
        last_author_affiliation="Third Institute of Synthetic Data, Far Away",
        source_locator="Fig S4D (synthetic)",
    )
    data.update(overrides)
    return synthetic.good_record(**data)


class TestEmptyDataset(unittest.TestCase):
    def test_no_records_no_tasks(self):
        self.assertEqual(evaluate_records([]), [])
        self.assertEqual(pair_verdicts([]), [])

    def test_the_summary_of_nothing_says_nothing_is_feasible(self):
        summary = summarize_verdicts([])
        self.assertEqual(summary["n_tasks"], 0)
        self.assertFalse(summary["any_feasible"])
        self.assertEqual(summary["outcomes_feasible"], [])
        self.assertEqual(summary["by_verdict"], {name: 0 for name in VERDICTS})


class TestSingleObservationIsInsufficient(unittest.TestCase):
    def setUp(self):
        self.verdict = verdict_of([synthetic.good_record()])

    def test_verdict(self):
        self.assertEqual(self.verdict.verdict, "insufficient")
        self.assertFalse(self.verdict.is_feasible)

    def test_the_reasons_name_the_single_observation_not_an_agreement_shortfall(self):
        self.assertIn("SINGLE_OBSERVATION", self.verdict.reasons)
        self.assertIn("NO_INDEPENDENT_REPLICATION", self.verdict.reasons)
        self.assertIn("AGREEMENT_UNDEFINED", self.verdict.reasons)
        self.assertNotIn("AGREEMENT_BELOW_FEASIBLE", self.verdict.reasons)

    def test_the_numbers_behind_the_verdict_are_carried(self):
        self.assertEqual(self.verdict.numbers["n_records"], 1)
        self.assertEqual(self.verdict.numbers["n_studies"], 1)
        self.assertEqual(self.verdict.numbers["n_labs"], 1)
        self.assertIsNone(self.verdict.numbers["agreement_studies"])


class TestTwoLabsAgreeingIsMarginalNotFeasible(unittest.TestCase):
    def setUp(self):
        self.verdict = verdict_of(
            [synthetic.good_record(), synthetic.other_lab_record()]
        )

    def test_verdict_is_marginal(self):
        self.assertEqual(self.verdict.verdict, "marginal")

    def test_unanimity_from_two_studies_does_not_clear_the_study_gate(self):
        self.assertEqual(self.verdict.numbers["agreement_studies"], 1.0)
        self.assertEqual(self.verdict.numbers["n_labs"], 2)
        self.assertIn("STUDIES_BELOW_FEASIBLE", self.verdict.reasons)
        self.assertNotIn("LABS_BELOW_FEASIBLE", self.verdict.reasons)


class TestThreeStudiesAcrossThreeLabsCanBeFeasible(unittest.TestCase):
    """The gates must be passable, or they measure nothing."""

    def setUp(self):
        self.verdict = verdict_of(
            [
                synthetic.good_record(),
                synthetic.other_lab_record(),
                third_lab_record(),
            ]
        )

    def test_verdict_is_feasible(self):
        self.assertEqual(self.verdict.verdict, "feasible")
        self.assertEqual(self.verdict.reasons, ("MEETS_ALL_FEASIBLE_GATES",))

    def test_the_gates_it_cleared(self):
        numbers = self.verdict.numbers
        self.assertGreaterEqual(numbers["n_studies"], MIN_STUDIES_FEASIBLE)
        self.assertGreaterEqual(numbers["n_labs"], MIN_LABS_FEASIBLE)
        self.assertGreaterEqual(
            numbers["agreement_studies"], MIN_DIRECTION_AGREEMENT_FEASIBLE
        )
        self.assertGreater(
            numbers["agreement_studies_ci"][0], CHANCE_DIRECTION_RATE
        )

    def test_the_rollup_names_the_perturbation_and_the_condition(self):
        pairs = pair_verdicts([self.verdict])
        self.assertEqual(len(pairs), 1)
        self.assertEqual(pairs[0].verdict, "feasible")
        # Not a valence: valence has left the claim key. The rollup names the
        # perturbation axis and the condition that earned the verdict.
        self.assertEqual(
            pairs[0].best_perturbation, "circuit_activity / increase / <null>"
        )
        self.assertIn("sex=male", pairs[0].best_condition)
        self.assertEqual(
            pairs[0].best_min_condition_specificity, N_CONDITION_SLOTS
        )
        self.assertEqual(pairs[0].n_claims, 1)
        self.assertEqual(pairs[0].n_labs_max, 3)

    def test_the_summary_lists_the_outcome(self):
        summary = summarize_verdicts([self.verdict])
        self.assertTrue(summary["any_feasible"])
        self.assertEqual(summary["outcomes_feasible"], ["syn:outcome_1"])


class TestTwoLabsDisagreeingIsInsufficient(unittest.TestCase):
    def setUp(self):
        self.verdict = verdict_of(
            [
                synthetic.good_record(),
                synthetic.other_lab_record(direction="decrease"),
            ]
        )

    def test_an_even_split_is_not_weak_evidence_it_is_no_evidence(self):
        self.assertEqual(self.verdict.verdict, "insufficient")
        self.assertIn("DIRECTION_CONTESTED", self.verdict.reasons)
        self.assertEqual(self.verdict.numbers["agreement_studies"], 0.5)

    def test_the_disagreement_is_reported_not_resolved(self):
        self.assertEqual(
            self.verdict.numbers["conflict_codes"], ["QC_CONFLICT_OPPOSED"]
        )
        self.assertIn("DIRECTION_DISAGREEMENT_PRESENT", self.verdict.notes)
        self.assertEqual(
            self.verdict.numbers["direction_counts"],
            {"decrease": 1, "increase": 1},
        )

    def test_a_null_result_against_an_effect_is_also_contested(self):
        verdict = verdict_of(
            [
                synthetic.good_record(),
                synthetic.other_lab_record(direction="no_change"),
            ]
        )
        self.assertEqual(verdict.verdict, "insufficient")
        self.assertIn("DIRECTION_CONTESTED", verdict.reasons)
        self.assertEqual(
            verdict.numbers["conflict_codes"], ["QC_CONFLICT_EFFECT_VS_NULL"]
        )

    def test_a_two_of_three_majority_is_marginal_not_insufficient(self):
        verdict = verdict_of(
            [
                synthetic.good_record(),
                synthetic.other_lab_record(),
                third_lab_record(direction="decrease"),
            ]
        )
        self.assertEqual(verdict.verdict, "marginal")
        self.assertAlmostEqual(verdict.numbers["agreement_studies"], 2.0 / 3.0)
        self.assertIn("AGREEMENT_BELOW_FEASIBLE", verdict.reasons)


class TestSameLabReplicationIsDiscounted(unittest.TestCase):
    def setUp(self):
        self.verdict = verdict_of(
            [synthetic.good_record(), synthetic.same_lab_other_study_record()]
        )

    def test_one_lab_can_never_be_feasible(self):
        self.assertEqual(self.verdict.verdict, "marginal")
        self.assertIn("LABS_BELOW_FEASIBLE", self.verdict.reasons)
        self.assertIn("SAME_LAB_ONLY", self.verdict.reasons)
        self.assertEqual(self.verdict.numbers["n_labs"], 1)
        self.assertEqual(self.verdict.numbers["n_studies"], 2)

    def test_three_studies_from_one_lab_are_still_not_feasible(self):
        verdict = verdict_of(
            [
                synthetic.good_record(),
                synthetic.same_lab_other_study_record(),
                synthetic.good_record(
                    study_id="synthetic:test-005",
                    doi="synthetic:test-005",
                    source_locator="Fig S5E (synthetic)",
                ),
            ]
        )
        self.assertEqual(verdict.numbers["n_studies"], 3)
        self.assertEqual(verdict.numbers["n_labs"], 1)
        self.assertEqual(verdict.numbers["agreement_studies"], 1.0)
        self.assertEqual(verdict.verdict, "marginal")
        self.assertIn("LABS_BELOW_FEASIBLE", verdict.reasons)


class TestConditionSpecificityQualifiesEveryVerdict(unittest.TestCase):
    """SCHEMA.md 5b's sibling rule: the number never travels without it."""

    def _unspecified(self, builder, **overrides):
        data = dict(sex=None, diet=None, age_group=None, challenge=None)
        data.update(overrides)
        return builder(**data)

    def test_a_feasible_verdict_on_unknown_conditions_says_so(self):
        verdict = verdict_of(
            [
                self._unspecified(synthetic.good_record),
                self._unspecified(synthetic.other_lab_record),
                self._unspecified(third_lab_record),
            ]
        )
        self.assertEqual(verdict.verdict, "feasible")
        self.assertEqual(verdict.min_condition_specificity, 0)
        self.assertTrue(verdict.matched_on_unknown_conditions)
        self.assertIn("CONDITIONS_ALL_UNKNOWN", verdict.notes)

    def test_it_is_a_note_and_not_a_gate(self):
        # Refusing it would mean inventing a specificity threshold, and this
        # module's discipline is that a threshold needs a stated reason.
        self.assertIsNone(MIN_CONDITION_SPECIFICITY_GATE)
        self.assertIsNone(THRESHOLDS["MIN_CONDITION_SPECIFICITY_GATE"])
        verdict = verdict_of(
            [
                self._unspecified(synthetic.good_record),
                self._unspecified(synthetic.other_lab_record),
                self._unspecified(third_lab_record),
            ]
        )
        self.assertNotIn("CONDITIONS_ALL_UNKNOWN", verdict.reasons)

    def test_a_partially_specified_cell_is_marked_too(self):
        verdict = verdict_of([synthetic.good_record(challenge=None)])
        self.assertEqual(verdict.min_condition_specificity, 3)
        self.assertFalse(verdict.matched_on_unknown_conditions)
        self.assertIn("CONDITIONS_PARTIALLY_UNKNOWN", verdict.notes)
        self.assertNotIn("CONDITIONS_ALL_UNKNOWN", verdict.notes)

    def test_the_specificity_rides_in_the_numbers_dict_next_to_the_fraction(self):
        verdict = verdict_of(
            [synthetic.good_record(), synthetic.other_lab_record()]
        )
        numbers = verdict.numbers
        self.assertEqual(numbers["agreement_studies"], 1.0)
        self.assertEqual(
            numbers["min_condition_specificity"], N_CONDITION_SLOTS
        )
        self.assertEqual(numbers["n_condition_slots"], N_CONDITION_SLOTS)
        self.assertFalse(numbers["matched_on_unknown_conditions"])
        # The stated form carries both, so a report cannot print one alone.
        self.assertIn("min condition specificity", numbers["agreement_studies_stated"])

    def test_the_summary_splits_the_verdict_counts_by_specificity(self):
        specified = [
            synthetic.good_record(),
            synthetic.other_lab_record(),
            third_lab_record(),
        ]
        unspecified = [
            self._unspecified(synthetic.good_record, population_id="syn:pop_beta"),
            self._unspecified(
                synthetic.other_lab_record, population_id="syn:pop_beta"
            ),
            self._unspecified(third_lab_record, population_id="syn:pop_beta"),
        ]
        summary = summarize_verdicts(evaluate_records(specified + unspecified))
        histogram = summary["by_verdict_and_min_condition_specificity"]["feasible"]
        self.assertEqual(histogram["0"], 1)
        self.assertEqual(histogram[str(N_CONDITION_SLOTS)], 1)
        self.assertEqual(summary["n_feasible_matched_on_unknown_conditions"], 1)


class TestTheKeyChangeCostsCellsAndSaysSo(unittest.TestCase):
    """A finer key means fewer cells clear a gate. Reported, not softened."""

    def _two_sexes(self):
        return [
            synthetic.good_record(sex="male"),
            synthetic.other_lab_record(sex="female"),
            third_lab_record(sex="male"),
        ]

    def test_the_pooled_key_promotes_a_cell_the_strict_key_splits(self):
        pooled = evaluate_records(self._two_sexes(), KEY_POOLED)
        strict = evaluate_records(self._two_sexes())
        self.assertEqual(len(pooled), 1)
        self.assertEqual(pooled[0].verdict, "feasible")
        # Three papers, two cohorts: two claims, neither with three studies.
        self.assertEqual(len(strict), 2)
        self.assertEqual(
            sorted(v.verdict for v in strict), ["insufficient", "marginal"]
        )
        self.assertFalse(any(v.is_feasible for v in strict))

    def test_no_threshold_moved_to_make_that_look_better(self):
        self.assertEqual(THRESHOLDS["MIN_STUDIES_MARGINAL"], 2)
        self.assertEqual(THRESHOLDS["MIN_STUDIES_FEASIBLE"], 3)
        self.assertEqual(THRESHOLDS["MIN_LABS_FEASIBLE"], 2)
        self.assertEqual(THRESHOLDS["MIN_RECORDS_FEASIBLE"], 3)
        self.assertEqual(THRESHOLDS["MIN_DIRECTION_AGREEMENT_MARGINAL"], 0.60)
        self.assertEqual(THRESHOLDS["MIN_DIRECTION_AGREEMENT_FEASIBLE"], 0.70)
        self.assertTrue(THRESHOLDS["REQUIRE_AGREEMENT_CI_ABOVE_CHANCE"])
        self.assertTrue(THRESHOLDS["EXCLUDE_TIER_B_AS_TARGET"])

    def test_a_pooled_verdict_knows_it_is_pooled(self):
        pooled = evaluate_records([synthetic.good_record()], KEY_POOLED)[0]
        self.assertEqual(pooled.valence, "activation")
        self.assertEqual(pooled.perturbation_label, "activation")
        self.assertIn("not keyed", pooled.condition_label)

    def test_a_molecular_claim_with_an_unresolved_target_is_noted(self):
        verdict = verdict_of(
            [synthetic.molecular_record(target_raw=None, target_id=None)]
        )
        self.assertIn("TARGET_UNRESOLVED", verdict.notes)

    def test_a_circuit_claim_with_no_target_is_not_noted(self):
        # A null target on a chemogenetic record is correct, not missing.
        verdict = verdict_of([synthetic.good_record()])
        self.assertNotIn("TARGET_UNRESOLVED", verdict.notes)


class TestOtherRefusals(unittest.TestCase):
    def test_a_tier_b_outcome_is_refused_as_a_training_target(self):
        verdict = verdict_of(
            [
                synthetic.good_record(outcome_tier="B"),
                synthetic.other_lab_record(outcome_tier="B"),
                third_lab_record(outcome_tier="B"),
            ]
        )
        self.assertEqual(verdict.verdict, "insufficient")
        self.assertIn("OUTCOME_IS_TIER_B_ANNOTATION", verdict.reasons)

    def test_missing_lab_identity_is_noted_without_changing_the_verdict(self):
        verdict = verdict_of(
            [
                synthetic.good_record(),
                synthetic.other_lab_record(),
                third_lab_record(
                    last_author=None, last_author_affiliation=None
                ),
            ]
        )
        self.assertIn("LAB_IDENTITY_MISSING", verdict.notes)
        self.assertEqual(verdict.numbers["n_records_without_lab"], 1)
        self.assertEqual(verdict.verdict, "feasible")


class TestThresholdsAreDocumentedConstants(unittest.TestCase):
    def test_every_threshold_is_published_with_the_verdict(self):
        verdict = verdict_of([synthetic.good_record()])
        self.assertEqual(verdict.to_dict()["thresholds"], THRESHOLDS)

    def test_every_reason_code_used_has_a_printed_explanation(self):
        records = [
            synthetic.good_record(),
            synthetic.other_lab_record(direction="decrease"),
            synthetic.same_lab_other_study_record(),
            third_lab_record(outcome_id="syn:outcome_2"),
        ]
        for verdict in evaluate_records(records):
            for code in verdict.reasons + verdict.notes:
                self.assertIn(code, REASON_NOTES, code)

    def test_the_gates_are_ordered_consistently(self):
        self.assertLess(
            THRESHOLDS["MIN_STUDIES_MARGINAL"], THRESHOLDS["MIN_STUDIES_FEASIBLE"]
        )
        self.assertLess(
            THRESHOLDS["MIN_LABS_MARGINAL"], THRESHOLDS["MIN_LABS_FEASIBLE"]
        )
        self.assertLess(
            THRESHOLDS["MIN_DIRECTION_AGREEMENT_MARGINAL"],
            THRESHOLDS["MIN_DIRECTION_AGREEMENT_FEASIBLE"],
        )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
