"""Independent replication, direction agreement, and the Wilson interval.

Schema v0.2.0: the claim key is nine components, both keys are reported side
by side, and an agreement fraction cannot be built without the minimum
condition specificity of the cell it came from.
"""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.feasibility.replication import (
    KEY_POOLED,
    KEY_STRICT,
    MIXED,
    N_CONDITION_SLOTS,
    Agreement,
    SpecificityMissing,
    build_replication_cells,
    compare_claim_keys,
    summarize_replication,
    wilson_interval,
)


def cell_of(records):
    cells = build_replication_cells(records)
    assert len(cells) == 1, "expected one claim, got %d" % len(cells)
    return cells[0]


class TestEmptyDataset(unittest.TestCase):
    def test_no_records_no_cells(self):
        self.assertEqual(build_replication_cells([]), [])

    def test_summary_of_nothing_reports_no_mean_agreement(self):
        summary = summarize_replication([])
        self.assertEqual(summary["n_cells"], 0)
        self.assertEqual(summary["n_cells_replicated_across_labs"], 0)
        # Not 0.0 and not 1.0: there is no agreement to average.
        self.assertIsNone(summary["mean_agreement_where_defined"])


class TestWilsonInterval(unittest.TestCase):
    def test_no_trials_no_interval(self):
        self.assertIsNone(wilson_interval(0, 0))

    def test_it_does_not_collapse_at_a_unanimous_two_of_two(self):
        low, high = wilson_interval(2, 2)
        self.assertLess(low, 1.0)
        self.assertGreater(low, 0.0)
        self.assertAlmostEqual(high, 1.0, places=6)
        # The Wald interval would give [1.0, 1.0] here; that is the point.
        self.assertAlmostEqual(low, 0.342, places=2)

    def test_it_narrows_as_evidence_accumulates(self):
        self.assertLess(wilson_interval(2, 2)[0], wilson_interval(3, 3)[0])
        self.assertLess(wilson_interval(3, 3)[0], wilson_interval(8, 8)[0])

    def test_it_stays_inside_the_unit_interval(self):
        for successes, trials in ((0, 1), (0, 3), (1, 3), (5, 5)):
            low, high = wilson_interval(successes, trials)
            self.assertGreaterEqual(low, 0.0)
            self.assertLessEqual(high, 1.0)

    def test_impossible_counts_are_refused(self):
        with self.assertRaises(ValueError):
            wilson_interval(3, 2)


class TestSingleObservation(unittest.TestCase):
    """One observation is not unanimity."""

    def setUp(self):
        self.cell = cell_of([synthetic.good_record()])

    def test_the_cell_exists_and_is_counted(self):
        self.assertEqual(self.cell.n_records, 1)
        self.assertEqual(self.cell.n_studies, 1)
        self.assertEqual(self.cell.n_labs, 1)
        self.assertTrue(self.cell.is_single_observation)

    def test_agreement_is_undefined_not_one(self):
        self.assertIsNone(self.cell.agreement_studies)
        self.assertIsNone(self.cell.agreement_studies_ci)
        self.assertIsNone(self.cell.agreement_labs)
        self.assertFalse(self.cell.agreement_is_defined)
        self.assertNotEqual(self.cell.agreement_studies, 1.0)  # None, not 1.0

    def test_it_is_not_replicated_and_not_contested(self):
        self.assertFalse(self.cell.replicated_across_studies)
        self.assertFalse(self.cell.replicated_across_labs)
        self.assertFalse(self.cell.same_lab_only)
        self.assertFalse(self.cell.is_contested)
        self.assertEqual(self.cell.conflict_codes, ())

    def test_the_modal_direction_is_still_reported(self):
        # The direction is known; only its *agreement* is undefined.
        self.assertEqual(self.cell.modal_direction, "increase")

    def test_a_summary_distinguishes_it_from_a_replicated_cell(self):
        summary = summarize_replication([self.cell])
        self.assertEqual(summary["n_cells_single_observation"], 1)
        self.assertEqual(summary["n_cells_agreement_undefined"], 1)
        self.assertEqual(summary["n_cells_agreement_defined"], 0)
        self.assertIsNone(summary["mean_agreement_where_defined"])


class TestTwoLabsAgreeing(unittest.TestCase):
    def setUp(self):
        self.cell = cell_of(
            [synthetic.fully_specified_record(), synthetic.other_lab_record(background_genotype="wild_type")]
        )

    def test_independence_is_counted_at_both_levels(self):
        self.assertEqual(self.cell.n_records, 2)
        self.assertEqual(self.cell.n_studies, 2)
        self.assertEqual(self.cell.n_labs, 2)
        self.assertTrue(self.cell.replicated_across_studies)
        self.assertTrue(self.cell.replicated_across_labs)
        self.assertFalse(self.cell.same_lab_only)

    def test_agreement_is_defined_and_unanimous(self):
        self.assertTrue(self.cell.agreement_is_defined)
        self.assertEqual(self.cell.agreement_studies.value, 1.0)
        self.assertEqual(self.cell.agreement_labs.value, 1.0)
        self.assertEqual(self.cell.modal_direction, "increase")
        self.assertFalse(self.cell.modal_direction_is_tied)

    def test_the_interval_is_wide_and_does_not_claim_certainty(self):
        low, high = self.cell.agreement_studies_ci
        self.assertLess(low, 0.5)
        self.assertAlmostEqual(high, 1.0, places=6)

    def test_no_conflict_is_recorded(self):
        self.assertEqual(self.cell.conflict_codes, ())
        self.assertFalse(self.cell.is_contested)


class TestTwoLabsDisagreeing(unittest.TestCase):
    def setUp(self):
        self.cell = cell_of(
            [
                synthetic.good_record(),
                synthetic.other_lab_record(direction="decrease"),
            ]
        )

    def test_both_sides_are_kept(self):
        self.assertEqual(self.cell.n_studies, 2)
        self.assertEqual(self.cell.n_labs, 2)
        self.assertEqual(self.cell.direction_counts, {"decrease": 1, "increase": 1})

    def test_agreement_is_half_and_the_mode_is_flagged_as_tied(self):
        self.assertEqual(self.cell.agreement_studies.value, 0.5)
        self.assertTrue(self.cell.modal_direction_is_tied)
        self.assertTrue(self.cell.is_contested)

    def test_the_qc_conflict_definition_is_reused(self):
        self.assertEqual(self.cell.conflict_codes, ("QC_CONFLICT_OPPOSED",))

    def test_an_effect_against_a_null_is_also_a_conflict_and_no_change_is_kept(self):
        cell = cell_of(
            [
                synthetic.good_record(),
                synthetic.other_lab_record(direction="no_change"),
            ]
        )
        self.assertEqual(cell.conflict_codes, ("QC_CONFLICT_EFFECT_VS_NULL",))
        self.assertEqual(cell.direction_counts, {"increase": 1, "no_change": 1})
        self.assertEqual(cell.agreement_studies.value, 0.5)

    def test_no_change_can_be_the_modal_direction(self):
        cell = cell_of(
            [
                synthetic.good_record(direction="no_change"),
                synthetic.other_lab_record(direction="no_change"),
            ]
        )
        self.assertEqual(cell.modal_direction, "no_change")
        self.assertEqual(cell.agreement_studies.value, 1.0)


class TestSameLabReplicationIsDiscounted(unittest.TestCase):
    def setUp(self):
        self.cell = cell_of(
            [synthetic.good_record(), synthetic.same_lab_other_study_record()]
        )

    def test_two_studies_from_one_lab_count_as_one_independent_lab(self):
        self.assertEqual(self.cell.n_studies, 2)
        self.assertEqual(self.cell.n_labs, 1)
        self.assertTrue(self.cell.replicated_across_studies)
        self.assertFalse(self.cell.replicated_across_labs)
        self.assertTrue(self.cell.same_lab_only)

    def test_study_level_agreement_exists_but_lab_level_does_not(self):
        self.assertEqual(self.cell.agreement_studies.value, 1.0)
        self.assertEqual(self.cell.n_labs_stating_direction, 1)
        self.assertIsNone(self.cell.agreement_labs)

    def test_the_summary_reports_it_as_discounted(self):
        summary = summarize_replication([self.cell])
        self.assertEqual(summary["n_cells_replicated_across_studies"], 1)
        self.assertEqual(summary["n_cells_replicated_across_labs"], 0)
        self.assertEqual(summary["n_cells_same_lab_only"], 1)


class TestClaimKeying(unittest.TestCase):
    """What separates two claims under schema v0.2.0.

    Not ``valence``: it has left the key. The perturbation axis
    ``(class, polarity, target)`` is what separates them, and for a circuit
    manipulation the polarity is what the valence used to say. So the *same*
    two records still land in two cells -- for the right reason, and
    ``perturbation_polarity`` is where the separation is visible.
    """

    def test_opposite_perturbation_polarities_are_separate_claims(self):
        cells = build_replication_cells(
            [
                synthetic.good_record(valence="activation"),
                synthetic.other_lab_record(valence="inhibition"),
            ]
        )
        self.assertEqual(len(cells), 2)
        polarities = sorted(cell.perturbation[1] for cell in cells)
        self.assertEqual(polarities, ["decrease", "increase"])
        for cell in cells:
            self.assertEqual(cell.n_studies, 1)
            self.assertIsNone(cell.agreement_studies)
            # valence is simply absent from the strict key.
            self.assertIsNone(cell.valence)

    def test_valence_alone_no_longer_separates_a_claim(self):
        # Two molecular records differing only in valence are ONE claim: a
        # molecular record's valence is honestly unclear and carries no
        # information, so keying on it would have fragmented nothing.
        cells = build_replication_cells(
            [
                synthetic.molecular_record(valence="unclear"),
                synthetic.molecular_record(
                    valence="activation",
                    perturbation_polarity="decrease",
                    study_id="synthetic:test-mol-002",
                    doi="synthetic:test-mol-002",
                    last_author="Other, C.",
                    last_author_affiliation="Second Institute, Elsewhere",
                ),
            ]
        )
        self.assertEqual(len(cells), 1)
        self.assertEqual(cells[0].n_studies, 2)

    def test_the_same_gene_and_a_different_gene_are_different_claims(self):
        cells = build_replication_cells(
            [
                synthetic.molecular_record(),
                synthetic.molecular_record(
                    target_raw="Synthgene2",
                    target_id="syn:synthgene2",
                    study_id="synthetic:test-mol-003",
                    doi="synthetic:test-mol-003",
                ),
            ]
        )
        self.assertEqual(len(cells), 2)

    def test_a_different_condition_is_a_different_claim(self):
        cells = build_replication_cells(
            [
                synthetic.good_record(sex="male"),
                synthetic.other_lab_record(sex="female"),
            ]
        )
        self.assertEqual(len(cells), 2)
        # This is the whole point of v0.2.0: a male cohort and a female cohort
        # are not two observations of one claim, so they cannot contradict.
        for cell in cells:
            self.assertEqual(cell.conflict_codes, ())

    def test_a_record_with_no_outcome_id_cannot_be_grouped(self):
        self.assertEqual(
            build_replication_cells([synthetic.good_record(outcome_id=None)]), []
        )

    def test_a_study_reporting_two_directions_is_mixed_and_never_agrees(self):
        records = [
            synthetic.good_record(source_locator="Fig 1A (synthetic)"),
            synthetic.good_record(
                source_locator="Fig 1B (synthetic)", direction="decrease"
            ),
            synthetic.other_lab_record(),
        ]
        cell = cell_of(records)
        self.assertEqual(cell.n_studies, 2)
        self.assertEqual(
            sorted(cell.study_directions.values()), ["increase", MIXED]
        )
        # The mixed study stays in the denominator and agrees with nothing.
        self.assertEqual(cell.n_studies_stating_direction, 2)
        self.assertEqual(cell.n_studies_agreeing, 1)
        self.assertEqual(cell.agreement_studies.value, 0.5)

    def test_not_reported_is_counted_and_kept_out_of_the_denominator(self):
        cell = cell_of(
            [
                synthetic.good_record(),
                synthetic.other_lab_record(direction="not_reported"),
            ]
        )
        self.assertEqual(cell.n_records, 2)
        self.assertEqual(cell.n_records_direction_not_reported, 1)
        self.assertEqual(cell.n_studies_stating_direction, 1)
        # One study states a direction, so agreement remains undefined.
        self.assertIsNone(cell.agreement_studies)

    def test_an_unknown_key_variant_is_refused(self):
        with self.assertRaises(KeyError):
            build_replication_cells([synthetic.good_record()], "whatever")

    def test_cells_are_deterministic_and_serialisable(self):
        records = [synthetic.fully_specified_record(), synthetic.other_lab_record(background_genotype="wild_type")]
        first = [c.to_dict() for c in build_replication_cells(records)]
        second = [c.to_dict() for c in build_replication_cells(list(reversed(records)))]
        self.assertEqual(first, second)
        self.assertEqual(first[0]["interval_method"], "wilson_score")


class TestAgreementCannotBeStatedWithoutItsSpecificity(unittest.TestCase):
    """The structural part: the number and its qualifier are one object."""

    def test_an_agreement_cannot_be_built_without_a_specificity(self):
        with self.assertRaises(SpecificityMissing):
            Agreement(
                value=1.0,
                level="studies",
                n_agreeing=2,
                n_units=2,
                min_condition_specificity=None,
            )

    def test_a_specificity_outside_the_slot_count_is_refused(self):
        with self.assertRaises(SpecificityMissing):
            Agreement(
                value=1.0,
                level="studies",
                n_agreeing=2,
                n_units=2,
                min_condition_specificity=N_CONDITION_SLOTS + 1,
            )

    def test_an_agreement_must_name_its_independent_unit(self):
        with self.assertRaises(SpecificityMissing):
            Agreement(
                value=1.0,
                level="records",
                n_agreeing=2,
                n_units=2,
                min_condition_specificity=0,
            )

    def test_there_is_no_implicit_conversion_to_a_bare_number(self):
        agreement = Agreement(
            value=1.0,
            level="studies",
            n_agreeing=2,
            n_units=2,
            min_condition_specificity=0,
        )
        with self.assertRaises(TypeError):
            float(agreement)
        with self.assertRaises(TypeError):
            "%f" % agreement

    def test_printing_it_always_carries_the_specificity(self):
        agreement = Agreement(
            value=1.0,
            level="studies",
            n_agreeing=2,
            n_units=2,
            min_condition_specificity=0,
        )
        text = str(agreement)
        self.assertIn("min condition specificity 0/%d" % N_CONDITION_SLOTS, text)
        self.assertIn("both unknown", text)
        self.assertTrue(agreement.matched_on_unknown_conditions)

    def test_a_fully_specified_cell_does_not_claim_to_have_matched_on_unknown(self):
        cell = cell_of([synthetic.fully_specified_record(), synthetic.other_lab_record(background_genotype="wild_type")])
        self.assertEqual(cell.min_condition_specificity, N_CONDITION_SLOTS)
        self.assertFalse(cell.matched_on_unknown_conditions)
        self.assertFalse(cell.agreement_studies.matched_on_unknown_conditions)

    def test_a_strict_cell_cannot_mix_specificities_but_a_pooled_one_can(self):
        """Where the minimum actually bites.

        Under the strict key the condition is *part of the key*, so every
        record in a cell has the same condition and therefore the same
        specificity: ``min == max`` by construction. Under the pooled key the
        condition is not in the key at all, so one cell can hold a fully
        specified record and a wholly unspecified one -- and that is precisely
        the heterogeneity the pooled key hides. The minimum is what makes it
        visible.
        """
        records = [
            synthetic.fully_specified_record(),
            synthetic.other_lab_record(
                sex=None,
                diet=None,
                age_group=None,
                challenge=None,
                background_genotype=None,
            ),
        ]
        strict = build_replication_cells(records, KEY_STRICT)
        self.assertEqual(len(strict), 2)
        for cell in strict:
            self.assertEqual(
                cell.min_condition_specificity, cell.max_condition_specificity
            )

        pooled = build_replication_cells(records, KEY_POOLED)
        self.assertEqual(len(pooled), 1)
        self.assertEqual(pooled[0].n_records, 2)
        self.assertEqual(pooled[0].min_condition_specificity, 0)
        self.assertEqual(pooled[0].max_condition_specificity, N_CONDITION_SLOTS)
        self.assertEqual(
            pooled[0].agreement_studies.min_condition_specificity, 0
        )
        self.assertTrue(
            pooled[0].agreement_studies.matched_on_unknown_conditions
        )

    def test_a_wholly_unspecified_strict_cell_says_it_matched_on_unknown(self):
        cell = cell_of(
            [
                synthetic.good_record(
                    sex=None, diet=None, age_group=None, challenge=None
                ),
                synthetic.other_lab_record(
                    sex=None, diet=None, age_group=None, challenge=None
                ),
            ]
        )
        self.assertEqual(cell.n_studies, 2)
        self.assertEqual(cell.min_condition_specificity, 0)
        self.assertTrue(cell.matched_on_unknown_conditions)
        self.assertTrue(cell.agreement_studies.matched_on_unknown_conditions)
        # It still replicates. The claim is just a weaker claim.
        self.assertTrue(cell.replicated_across_labs)
        self.assertEqual(cell.agreement_studies.value, 1.0)

    def test_the_summary_breaks_replication_down_by_specificity(self):
        specified = cell_of([synthetic.fully_specified_record(), synthetic.other_lab_record(background_genotype="wild_type")])
        unspecified = cell_of(
            [
                synthetic.good_record(
                    population_id="syn:pop_beta",
                    sex=None,
                    diet=None,
                    age_group=None,
                    challenge=None,
                ),
                synthetic.other_lab_record(
                    population_id="syn:pop_beta",
                    sex=None,
                    diet=None,
                    age_group=None,
                    challenge=None,
                ),
            ]
        )
        summary = summarize_replication([specified, unspecified])
        labs = summary["replicated_across_labs_by_min_condition_specificity"]
        self.assertEqual(labs["0"], 1)
        self.assertEqual(labs[str(N_CONDITION_SLOTS)], 1)
        self.assertEqual(summary["n_cells_matched_on_unknown_conditions"], 1)
        # Every slot is materialised, zeros included, so "nothing replicates
        # at specificity 2" is visible rather than absent.
        self.assertEqual(
            sorted(labs), [str(k) for k in range(N_CONDITION_SLOTS + 1)]
        )

    def test_a_specificity_bucket_with_no_cells_has_no_mean_rather_than_zero(self):
        cell = cell_of([synthetic.fully_specified_record(), synthetic.other_lab_record(background_genotype="wild_type")])
        means = summarize_replication([cell])[
            "mean_agreement_by_min_condition_specificity"
        ]
        self.assertEqual(means[str(N_CONDITION_SLOTS)], 1.0)
        self.assertIsNone(means["0"])


class TestBothKeysAreReportedSideBySide(unittest.TestCase):
    """The v0.2.0 change has to stay falsifiable."""

    def _opposed_by_condition(self):
        # Two labs, opposite directions, different sexes: one claim under the
        # pooled key (a contradiction), two claims under the strict key (two
        # incomparable experiments).
        return [
            synthetic.good_record(sex="male", direction="increase"),
            synthetic.other_lab_record(sex="female", direction="decrease"),
        ]

    def test_the_pooled_key_sees_a_contradiction_the_strict_key_does_not(self):
        comparison = compare_claim_keys(self._opposed_by_condition())
        self.assertEqual(comparison.pooled.n_cells, 1)
        self.assertEqual(comparison.strict.n_cells, 2)
        self.assertEqual(comparison.pooled.n_cells_directly_opposed, 1)
        self.assertEqual(comparison.strict.n_cells_directly_opposed, 0)
        self.assertEqual(comparison.n_opposed_cells_resolved_by_split, 1)
        self.assertEqual(comparison.n_opposed_cells_surviving, 0)

    def test_the_cost_is_counted_as_well_as_the_benefit(self):
        comparison = compare_claim_keys(self._opposed_by_condition())
        # The same split that resolved the contradiction cost the one
        # multi-study cell. Both numbers are reported; neither is netted out.
        self.assertEqual(comparison.pooled.n_cells_multi_study, 1)
        self.assertEqual(comparison.strict.n_cells_multi_study, 0)
        self.assertEqual(comparison.n_multi_study_cells_lost, 1)
        self.assertEqual(comparison.n_pooled_cells_split, 1)

    def test_the_three_ways_a_cell_can_be_lost_are_kept_apart(self):
        """The honest part of the cost accounting.

        Two of the three reasons are correct separations of work that was
        never comparable. The third -- one paper states a condition, the other
        leaves it unknown -- is the stricter key's own over-split, and the
        corpus cannot say whether the two experiments actually matched. It is
        never counted as heterogeneity resolved.
        """
        # (a) different perturbations: correct separation.
        perturbation = compare_claim_keys(
            [
                synthetic.molecular_record(),
                synthetic.molecular_record(
                    target_raw="Synthgene2",
                    target_id="syn:synthgene2",
                    study_id="synthetic:mol-2",
                    doi="synthetic:mol-2",
                    last_author="Other, C.",
                    last_author_affiliation="Second Institute, Elsewhere",
                ),
            ]
        )
        self.assertEqual(perturbation.n_multi_study_cells_lost, 1)
        self.assertEqual(perturbation.n_lost_to_perturbation_split, 1)
        self.assertEqual(perturbation.cost_that_cannot_be_adjudicated, 0)

        # (b) male against female: correct separation.
        stated = compare_claim_keys(
            [
                synthetic.good_record(sex="male"),
                synthetic.other_lab_record(sex="female"),
            ]
        )
        self.assertEqual(stated.n_lost_to_stated_condition_difference, 1)
        self.assertEqual(stated.cost_that_cannot_be_adjudicated, 0)

        # (c) male against silence: the corpus cannot say.
        unknown = compare_claim_keys(
            [
                synthetic.good_record(sex="male"),
                synthetic.other_lab_record(sex=None),
            ]
        )
        self.assertEqual(unknown.n_multi_study_cells_lost, 1)
        self.assertEqual(unknown.n_lost_to_unknown_against_stated, 1)
        self.assertEqual(unknown.cost_that_cannot_be_adjudicated, 1)
        self.assertEqual(unknown.n_lost_to_perturbation_split, 0)
        self.assertEqual(unknown.n_lost_to_stated_condition_difference, 0)

    def test_the_loss_buckets_partition_the_losses(self):
        comparison = compare_claim_keys(
            [
                synthetic.good_record(sex="male"),
                synthetic.other_lab_record(sex=None),
                synthetic.good_record(
                    population_id="syn:pop_beta", sex="male"
                ),
                synthetic.other_lab_record(
                    population_id="syn:pop_beta", sex="female"
                ),
            ]
        )
        self.assertEqual(
            comparison.n_lost_to_perturbation_split
            + comparison.n_lost_to_stated_condition_difference
            + comparison.n_lost_to_unknown_against_stated,
            comparison.n_multi_study_cells_lost,
        )

    def test_a_perturbation_difference_outranks_a_condition_difference(self):
        # A cell differing in both would have been split by the perturbation
        # alone, so it is attributed there and the unknown-vs-stated bucket
        # stays an honest count of what only the condition split.
        comparison = compare_claim_keys(
            [
                synthetic.molecular_record(sex="male"),
                synthetic.molecular_record(
                    target_raw="Synthgene2",
                    target_id="syn:synthgene2",
                    sex=None,
                    study_id="synthetic:mol-2",
                    doi="synthetic:mol-2",
                    last_author="Other, C.",
                    last_author_affiliation="Second Institute, Elsewhere",
                ),
            ]
        )
        self.assertEqual(comparison.n_lost_to_perturbation_split, 1)
        self.assertEqual(comparison.n_lost_to_unknown_against_stated, 0)

    def test_a_real_disagreement_survives_the_split(self):
        # Same sex, same perturbation, opposite directions: not heterogeneity.
        comparison = compare_claim_keys(
            [
                synthetic.good_record(direction="increase"),
                synthetic.other_lab_record(direction="decrease"),
            ]
        )
        self.assertEqual(comparison.pooled.n_cells_directly_opposed, 1)
        self.assertEqual(comparison.strict.n_cells_directly_opposed, 1)
        self.assertEqual(comparison.n_opposed_cells_resolved_by_split, 0)
        self.assertEqual(comparison.n_opposed_cells_surviving, 1)

    def test_a_pooled_cell_is_labelled_as_pooled_and_carries_its_valence(self):
        cells = build_replication_cells(
            [synthetic.good_record()], KEY_POOLED
        )
        self.assertEqual(len(cells), 1)
        self.assertTrue(cells[0].is_pooled)
        self.assertEqual(cells[0].valence, "activation")
        self.assertIsNone(cells[0].perturbation)
        self.assertIn("not keyed", cells[0].condition_label)

    def test_a_strict_cell_says_which_key_it_came_from(self):
        cell = cell_of([synthetic.good_record()])
        self.assertEqual(cell.key_variant, KEY_STRICT)
        self.assertFalse(cell.is_pooled)
        self.assertEqual(len(cell.claim), len(cell.claim_fields))
        self.assertEqual(cell.claim_dict["population"], "syn:pop_alpha")
        self.assertEqual(cell.claim_dict["outcome_id"], "syn:outcome_1")

    def test_a_perturbation_is_never_returned_without_its_class(self):
        cell = cell_of([synthetic.molecular_record()])
        self.assertEqual(
            cell.perturbation, ("molecular", "decrease", "syn:synthgene1")
        )
        self.assertIn("molecular", cell.perturbation_label)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
