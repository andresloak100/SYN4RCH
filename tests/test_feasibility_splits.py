"""Leakage-aware folds, and the guard against an unlabelled metric."""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.feasibility.splits import (
    GROUPINGS,
    Fold,
    FoldSet,
    GroupedMetric,
    UnlabelledResultError,
    all_fold_sets,
    attach_metric,
    by_lab,
    by_study,
    leave_one_population_out,
    require_grouping,
)


class TestEmptyDataset(unittest.TestCase):
    def test_every_grouping_returns_an_empty_labelled_fold_set(self):
        for builder in (by_study, by_lab, leave_one_population_out):
            fold_set = builder([])
            self.assertEqual(fold_set.n_folds, 0)
            self.assertEqual(fold_set.n_groups, 0)
            self.assertFalse(fold_set.is_usable)
            self.assertIn(fold_set.grouping, GROUPINGS)

    def test_all_fold_sets_is_keyed_by_grouping(self):
        self.assertEqual(sorted(all_fold_sets([])), sorted(GROUPINGS))


class TestGroupingsSeparateTrainAndTest(unittest.TestCase):
    def setUp(self):
        self.records = [
            synthetic.good_record(),
            synthetic.other_lab_record(),
            synthetic.same_lab_other_study_record(),
        ]
        self.ids = {r.record_id for r in self.records}

    def test_by_study_leaves_one_study_out_per_fold(self):
        fold_set = by_study(self.records)
        self.assertEqual(fold_set.grouping, "by_study")
        self.assertEqual(fold_set.n_folds, 3)
        self.assertTrue(fold_set.is_usable)
        for fold in fold_set.folds:
            self.assertEqual(
                set(fold.train_record_ids) | set(fold.test_record_ids), self.ids
            )
            self.assertFalse(
                set(fold.train_record_ids) & set(fold.test_record_ids)
            )

    def test_by_lab_folds_are_fewer_than_by_study_folds(self):
        # Two of the three synthetic studies share a lab, so lab grouping must
        # merge them: this is the leakage SCHEMA.md 5.2 is about.
        self.assertEqual(by_study(self.records).n_folds, 3)
        self.assertEqual(by_lab(self.records).n_folds, 2)

    def test_by_lab_keeps_one_labs_studies_on_the_same_side(self):
        fold_set = by_lab(self.records)
        same_lab_ids = {
            synthetic.good_record().record_id,
            synthetic.same_lab_other_study_record().record_id,
        }
        for fold in fold_set.folds:
            test = set(fold.test_record_ids)
            self.assertIn(
                test & same_lab_ids,
                (set(), same_lab_ids),
                "one lab's records were split across the boundary",
            )

    def test_leave_one_population_out(self):
        records = [
            synthetic.good_record(population_id="syn:pop_alpha"),
            synthetic.other_lab_record(population_id="syn:pop_beta"),
        ]
        fold_set = leave_one_population_out(records)
        self.assertEqual(fold_set.grouping, "leave_one_population_out")
        self.assertEqual(fold_set.n_folds, 2)
        self.assertEqual(fold_set.groups, ("syn:pop_alpha", "syn:pop_beta"))

    def test_one_group_is_not_cross_validation(self):
        fold_set = leave_one_population_out([synthetic.good_record()])
        self.assertEqual(fold_set.n_folds, 1)
        self.assertTrue(fold_set.folds[0].is_degenerate)
        self.assertFalse(fold_set.is_usable)

    def test_every_fold_set_states_the_shared_cohort_caveat(self):
        for fold_set in all_fold_sets(self.records).values():
            self.assertIn("Shared animal cohort", fold_set.note)


class TestUngroupableRecordsAreExcludedAndNamed(unittest.TestCase):
    def setUp(self):
        self.anonymous = synthetic.good_record(
            study_id="synthetic:test-011",
            doi="synthetic:test-011",
            last_author=None,
            last_author_affiliation=None,
            source_locator="Fig S11A (synthetic)",
        )
        self.records = [
            synthetic.good_record(),
            synthetic.other_lab_record(),
            self.anonymous,
        ]

    def test_a_record_with_no_lab_identity_is_in_no_fold(self):
        fold_set = by_lab(self.records)
        self.assertEqual(
            fold_set.ungroupable_record_ids, (self.anonymous.record_id,)
        )
        self.assertEqual(fold_set.n_ungroupable, 1)
        for fold in fold_set.folds:
            self.assertNotIn(self.anonymous.record_id, fold.train_record_ids)
            self.assertNotIn(self.anonymous.record_id, fold.test_record_ids)

    def test_unknown_labs_are_not_pooled_into_one_anonymous_lab(self):
        second = synthetic.good_record(
            study_id="synthetic:test-012",
            doi="synthetic:test-012",
            last_author=None,
            last_author_affiliation=None,
            source_locator="Fig S12A (synthetic)",
        )
        fold_set = by_lab(self.records + [second])
        self.assertEqual(fold_set.n_groups, 2)  # the two named labs only
        self.assertEqual(fold_set.n_ungroupable, 2)

    def test_the_same_record_is_groupable_by_study(self):
        self.assertEqual(by_study(self.records).n_ungroupable, 0)


class TestTheGroupingLabelGuard(unittest.TestCase):
    """SCHEMA.md 5: an unlabelled cross-validated number is an invalid result."""

    def setUp(self):
        self.fold_set = by_lab(
            [synthetic.good_record(), synthetic.other_lab_record()]
        )

    def test_attaching_to_a_labelled_fold_set_works(self):
        metric = attach_metric(self.fold_set, "direction_accuracy", 0.75)
        self.assertEqual(metric.grouping, "by_lab")
        self.assertEqual(metric.value, 0.75)
        self.assertEqual(metric.n_folds, 2)

    def test_attaching_to_a_bare_list_of_folds_raises(self):
        bare = list(self.fold_set.folds)
        with self.assertRaises(UnlabelledResultError):
            attach_metric(bare, "direction_accuracy", 0.75)
        with self.assertRaises(UnlabelledResultError):
            attach_metric(tuple(bare), "direction_accuracy", 0.75)

    def test_attaching_to_none_raises(self):
        with self.assertRaises(UnlabelledResultError):
            attach_metric(None, "direction_accuracy", 0.75)

    def test_an_object_with_an_empty_grouping_raises(self):
        class Unlabelled:
            grouping = ""
            folds = ()

        with self.assertRaises(UnlabelledResultError):
            attach_metric(Unlabelled(), "direction_accuracy", 0.75)

    def test_an_unrecognised_grouping_raises(self):
        class Mislabelled:
            grouping = "random_80_20"
            folds = ()

        with self.assertRaises(UnlabelledResultError):
            attach_metric(Mislabelled(), "direction_accuracy", 0.75)
        with self.assertRaises(UnlabelledResultError):
            require_grouping(Mislabelled())

    def test_a_label_with_no_folds_raises(self):
        class NoFolds:
            grouping = "by_lab"
            folds = None

        with self.assertRaises(UnlabelledResultError):
            attach_metric(NoFolds(), "direction_accuracy", 0.75)

    def test_a_fold_set_cannot_be_constructed_with_an_unknown_grouping(self):
        with self.assertRaises(UnlabelledResultError):
            FoldSet(grouping="holdout")

    def test_a_grouped_metric_cannot_be_constructed_unlabelled(self):
        with self.assertRaises(UnlabelledResultError):
            GroupedMetric(
                name="accuracy",
                value=0.9,
                grouping="whatever",
                n_folds=2,
                n_usable_folds=2,
                n_groups=2,
            )
        with self.assertRaises(UnlabelledResultError):
            GroupedMetric(
                name="",
                value=0.9,
                grouping="by_lab",
                n_folds=2,
                n_usable_folds=2,
                n_groups=2,
            )

    def test_printing_the_metric_always_names_the_grouping(self):
        text = str(attach_metric(self.fold_set, "direction_accuracy", 0.75))
        self.assertIn("grouping=by_lab", text)
        self.assertIn("direction_accuracy", text)

    def test_the_metric_does_not_silently_become_a_bare_number(self):
        metric = attach_metric(self.fold_set, "direction_accuracy", 0.75)
        self.assertFalse(hasattr(metric, "__float__"))
        with self.assertRaises(TypeError):
            float(metric)

    def test_ungroupable_records_are_named_on_the_metric(self):
        fold_set = by_lab(
            [
                synthetic.good_record(),
                synthetic.other_lab_record(),
                synthetic.good_record(
                    study_id="synthetic:test-013",
                    doi="synthetic:test-013",
                    last_author=None,
                    last_author_affiliation=None,
                    source_locator="Fig S13A (synthetic)",
                ),
            ]
        )
        metric = attach_metric(fold_set, "direction_accuracy", 0.5)
        self.assertEqual(metric.n_ungroupable, 1)
        self.assertIn("ungroupable", str(metric))


class TestFoldShape(unittest.TestCase):
    def test_a_fold_knows_it_is_degenerate(self):
        fold = Fold(index=0, grouping="by_study", held_out="s", test_record_ids=("a",))
        self.assertTrue(fold.is_degenerate)
        self.assertEqual(fold.n_train, 0)
        self.assertEqual(fold.n_test, 1)

    def test_fold_sets_serialise_with_their_grouping(self):
        data = by_study([synthetic.good_record()]).to_dict()
        self.assertEqual(data["grouping"], "by_study")
        self.assertIn("note", data)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
