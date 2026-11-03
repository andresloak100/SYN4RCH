"""Independent replication, direction agreement, and the Wilson interval."""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.feasibility.replication import (
    MIXED,
    build_replication_cells,
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
        self.assertNotEqual(self.cell.agreement_studies, 1.0)

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
            [synthetic.good_record(), synthetic.other_lab_record()]
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
        self.assertEqual(self.cell.agreement_studies, 1.0)
        self.assertEqual(self.cell.agreement_labs, 1.0)
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
        self.assertEqual(self.cell.agreement_studies, 0.5)
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
        self.assertEqual(cell.agreement_studies, 0.5)

    def test_no_change_can_be_the_modal_direction(self):
        cell = cell_of(
            [
                synthetic.good_record(direction="no_change"),
                synthetic.other_lab_record(direction="no_change"),
            ]
        )
        self.assertEqual(cell.modal_direction, "no_change")
        self.assertEqual(cell.agreement_studies, 1.0)


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
        self.assertEqual(self.cell.agreement_studies, 1.0)
        self.assertEqual(self.cell.n_labs_stating_direction, 1)
        self.assertIsNone(self.cell.agreement_labs)

    def test_the_summary_reports_it_as_discounted(self):
        summary = summarize_replication([self.cell])
        self.assertEqual(summary["n_cells_replicated_across_studies"], 1)
        self.assertEqual(summary["n_cells_replicated_across_labs"], 0)
        self.assertEqual(summary["n_cells_same_lab_only"], 1)


class TestClaimKeying(unittest.TestCase):
    def test_valences_are_separate_claims(self):
        cells = build_replication_cells(
            [
                synthetic.good_record(valence="activation"),
                synthetic.other_lab_record(valence="inhibition"),
            ]
        )
        self.assertEqual(len(cells), 2)
        for cell in cells:
            self.assertEqual(cell.n_studies, 1)
            self.assertIsNone(cell.agreement_studies)

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
        self.assertEqual(cell.agreement_studies, 0.5)

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

    def test_cells_are_deterministic_and_serialisable(self):
        records = [synthetic.good_record(), synthetic.other_lab_record()]
        first = [c.to_dict() for c in build_replication_cells(records)]
        second = [c.to_dict() for c in build_replication_cells(list(reversed(records)))]
        self.assertEqual(first, second)
        self.assertEqual(first[0]["interval_method"], "wilson_score")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
