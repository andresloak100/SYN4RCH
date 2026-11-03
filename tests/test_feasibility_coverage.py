"""Coverage counting: axes, explicit zeros, sparsity, and what is held out."""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.feasibility.coverage import (
    build_matrix,
    exclusion_reasons,
    has_measured_outcome,
    has_usable_intervention,
    is_usable,
    partition_records,
)


class TestEmptyDataset(unittest.TestCase):
    def test_no_records_gives_no_axes_and_no_cells(self):
        matrix = build_matrix([])
        self.assertEqual(matrix.n_records_read, 0)
        self.assertEqual(matrix.n_records_usable, 0)
        self.assertEqual(matrix.n_populations, 0)
        self.assertEqual(matrix.n_outcomes, 0)
        self.assertEqual(matrix.n_cells_total, 0)
        self.assertEqual(matrix.cells, {})
        self.assertEqual(matrix.n_studies, 0)
        self.assertEqual(matrix.n_labs, 0)

    def test_sparsity_of_an_empty_matrix_is_undefined_not_zero(self):
        # 0.0 would read as "no cell is empty", i.e. perfect coverage.
        self.assertIsNone(build_matrix([]).sparsity)

    def test_an_absent_cell_is_an_explicit_zero(self):
        cell = build_matrix([]).cell("syn:pop_alpha", "syn:outcome_1")
        self.assertEqual(cell.n_records, 0)
        self.assertTrue(cell.is_empty)


class TestUsability(unittest.TestCase):
    def test_the_baseline_synthetic_record_is_usable(self):
        record = synthetic.good_record()
        self.assertEqual(exclusion_reasons(record), [])
        self.assertTrue(is_usable(record))
        self.assertTrue(has_usable_intervention(record))
        self.assertTrue(has_measured_outcome(record))

    def test_no_change_is_a_measured_outcome_and_is_never_excluded(self):
        record = synthetic.good_record(direction="no_change")
        self.assertEqual(exclusion_reasons(record), [])
        self.assertTrue(has_measured_outcome(record))

    def test_each_exclusion_reason_is_detected(self):
        cases = (
            ({"review_status": "rejected"}, "REJECTED_IN_REVIEW"),
            (
                {"population_id": None, "population_raw": None},
                "NO_POPULATION_IDENTITY",
            ),
            ({"outcome_id": None}, "NO_OUTCOME_ID"),
            ({"method": None}, "NO_PERTURBATION_METHOD"),
            ({"valence": None}, "VALENCE_MISSING"),
            ({"valence": "unclear"}, "VALENCE_UNCLEAR"),
            ({"direction": None}, "DIRECTION_MISSING"),
            ({"direction": "not_reported"}, "DIRECTION_NOT_REPORTED"),
        )
        for override, expected in cases:
            self.assertIn(
                expected, exclusion_reasons(synthetic.good_record(**override)), override
            )

    def test_unclear_valence_is_counted_not_silently_dropped(self):
        records = [synthetic.good_record(), synthetic.good_record(valence="unclear")]
        matrix = build_matrix(records)
        self.assertEqual(matrix.n_records_read, 2)
        self.assertEqual(matrix.n_records_usable, 1)
        self.assertEqual(matrix.n_records_excluded, 1)
        self.assertEqual(matrix.excluded_counts["VALENCE_UNCLEAR"], 1)

    def test_partition_counts_reasons_and_records_separately(self):
        record = synthetic.good_record(valence=None, direction=None)
        usable, counts, excluded = partition_records([record])
        self.assertEqual(usable, [])
        self.assertEqual(len(excluded), 1)
        self.assertEqual(counts["VALENCE_MISSING"], 1)
        self.assertEqual(counts["DIRECTION_MISSING"], 1)


class TestMatrixKeepsExplicitZeros(unittest.TestCase):
    def _two_by_two(self):
        # Two populations x two outcomes, but only the diagonal measured.
        return [
            synthetic.good_record(
                population_id="syn:pop_alpha", outcome_id="syn:outcome_1"
            ),
            synthetic.other_lab_record(
                population_id="syn:pop_beta", outcome_id="syn:outcome_2"
            ),
        ]

    def test_every_pair_gets_a_cell(self):
        matrix = build_matrix(self._two_by_two())
        self.assertEqual(matrix.n_populations, 2)
        self.assertEqual(matrix.n_outcomes, 2)
        self.assertEqual(matrix.n_cells_total, 4)
        self.assertEqual(len(matrix.cells), 4)
        self.assertEqual(matrix.n_cells_populated, 2)
        self.assertEqual(matrix.n_cells_empty, 2)

    def test_the_unmeasured_cells_are_present_and_zero(self):
        matrix = build_matrix(self._two_by_two())
        empty = matrix.empty_cells()
        self.assertEqual(len(empty), 2)
        for cell in empty:
            self.assertEqual(cell.n_records, 0)
            self.assertEqual(cell.n_studies, 0)
            self.assertEqual(cell.n_labs, 0)
        self.assertEqual(
            {(c.population, c.outcome) for c in empty},
            {
                ("syn:pop_alpha", "syn:outcome_2"),
                ("syn:pop_beta", "syn:outcome_1"),
            },
        )

    def test_sparsity_is_reported(self):
        self.assertAlmostEqual(build_matrix(self._two_by_two()).sparsity, 0.5)


class TestCellCountsStudiesAndLabs(unittest.TestCase):
    def test_lab_identity_is_last_author_plus_affiliation(self):
        matrix = build_matrix(
            [synthetic.good_record(), synthetic.other_lab_record()]
        )
        cell = matrix.cell("syn:pop_alpha", "syn:outcome_1")
        self.assertEqual(cell.n_records, 2)
        self.assertEqual(cell.n_studies, 2)
        self.assertEqual(cell.n_labs, 2)
        self.assertTrue(cell.is_multi_lab)

    def test_two_studies_from_one_lab_are_one_lab(self):
        matrix = build_matrix(
            [synthetic.good_record(), synthetic.same_lab_other_study_record()]
        )
        cell = matrix.cell("syn:pop_alpha", "syn:outcome_1")
        self.assertEqual(cell.n_studies, 2)
        self.assertEqual(cell.n_labs, 1)
        self.assertFalse(cell.is_multi_lab)

    def test_same_last_author_at_a_different_affiliation_is_a_different_lab(self):
        matrix = build_matrix(
            [
                synthetic.good_record(),
                synthetic.good_record(
                    study_id="synthetic:test-009",
                    doi="synthetic:test-009",
                    last_author_affiliation="Another Institute, Elsewhere",
                    source_locator="Fig S9A (synthetic)",
                ),
            ]
        )
        self.assertEqual(matrix.cell("syn:pop_alpha", "syn:outcome_1").n_labs, 2)

    def test_records_without_lab_identity_are_counted_not_pooled(self):
        matrix = build_matrix(
            [
                synthetic.good_record(),
                synthetic.good_record(
                    study_id="synthetic:test-010",
                    doi="synthetic:test-010",
                    last_author=None,
                    last_author_affiliation=None,
                    source_locator="Fig S10A (synthetic)",
                ),
            ]
        )
        cell = matrix.cell("syn:pop_alpha", "syn:outcome_1")
        self.assertEqual(cell.n_records, 2)
        self.assertEqual(cell.n_labs, 1)
        self.assertEqual(cell.n_records_without_lab, 1)
        self.assertEqual(matrix.n_records_without_lab, 1)

    def test_unresolved_populations_group_under_the_raw_prefix(self):
        matrix = build_matrix(
            [
                synthetic.good_record(
                    population_id=None, population_raw="unmapped population"
                )
            ]
        )
        self.assertEqual(matrix.populations, ("raw:unmapped population",))
        self.assertEqual(matrix.n_populations_unresolved, 1)

    def test_matrix_serialises(self):
        data = build_matrix([synthetic.good_record()]).to_dict()
        self.assertEqual(data["n_records_usable"], 1)
        self.assertEqual(len(data["populated_cells"]), 1)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
