"""The task definition: features, labels, the exclusion ledger, groupings."""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.models.task import (
    DEFAULT_GROUPINGS,
    EXCLUSION_REASONS,
    GROUPINGS,
    LABELS,
    TASK,
    ExclusionLedger,
    ExclusionReason,
    build_dataset,
    exclusion_reasons,
    extract_features,
    extract_label,
    load_dataset,
)


class TestTaskStatement(unittest.TestCase):
    def test_three_labels_in_a_fixed_order(self):
        self.assertEqual(LABELS, ("increase", "decrease", "no_change"))
        self.assertEqual(TASK.labels, LABELS)

    def test_the_task_refuses_a_rejuvenation_score_and_tier_b(self):
        refused = " ".join(TASK.refuses).lower()
        self.assertIn("rejuvenation score", refused)
        self.assertIn("tier b", refused)
        self.assertIn("grouping", refused)

    def test_features_are_the_five_the_task_names(self):
        self.assertEqual(
            TASK.feature_fields,
            ("population_id", "valence", "outcome_id", "species", "acuity"),
        )


class TestFeatureAndLabelExtraction(unittest.TestCase):
    def test_features_come_from_the_record_unchanged(self):
        record = synthetic.good_record(
            population_id="arc_agrp",
            valence="inhibition",
            outcome_id="food_intake",
            species="mus_musculus",
            acuity="chronic",
        )
        features = extract_features(record)
        self.assertEqual(features.population_id, "arc_agrp")
        self.assertEqual(features.valence, "inhibition")
        self.assertEqual(features.outcome_id, "food_intake")
        self.assertEqual(features.species, "mus_musculus")
        self.assertEqual(features.acuity, "chronic")

    def test_an_unresolved_population_keeps_a_raw_prefixed_key(self):
        record = synthetic.good_record(
            population_id=None, population_raw="POMC-adjacent neurons"
        )
        features = extract_features(record)
        self.assertIsNone(features.population_id)
        self.assertTrue(features.population_key.startswith("raw:"))

    def test_the_three_labellable_directions_become_labels(self):
        for direction in LABELS:
            record = synthetic.good_record(direction=direction)
            self.assertEqual(extract_label(record), direction)

    def test_unlabellable_directions_return_none(self):
        for direction in ("not_reported", "non_monotonic"):
            record = synthetic.good_record(direction=direction)
            self.assertIsNone(extract_label(record))
        self.assertIsNone(extract_label(synthetic.good_record(direction=None)))

    def test_nulls_are_visible_in_a_cell_key(self):
        features = extract_features(synthetic.good_record(acuity=None))
        self.assertEqual(features.cell("acuity"), ("<null>",))


class TestExclusionsAreCountedNotDropped(unittest.TestCase):
    def test_not_reported_and_non_monotonic_are_excluded_with_their_reason(self):
        records = [
            synthetic.good_record(direction="increase", source_locator="Fig 1"),
            synthetic.good_record(direction="not_reported", source_locator="Fig 2"),
            synthetic.good_record(direction="non_monotonic", source_locator="Fig 3"),
        ]
        dataset = build_dataset(records)
        self.assertEqual(len(dataset), 1)
        primary = dataset.ledger.primary
        self.assertEqual(primary[ExclusionReason.direction_not_reported], 1)
        self.assertEqual(primary[ExclusionReason.direction_non_monotonic], 1)
        self.assertEqual(dataset.ledger.n_excluded, 2)
        self.assertEqual(dataset.ledger.n_input, 3)

    def test_the_ledger_always_balances(self):
        records = [
            synthetic.good_record(direction="increase", source_locator="a"),
            synthetic.good_record(direction="not_reported", source_locator="b"),
            synthetic.good_record(outcome_id=None, source_locator="c"),
            synthetic.good_record(valence=None, source_locator="d"),
        ]
        dataset = build_dataset(records)
        self.assertEqual(
            dataset.ledger.n_examples + dataset.ledger.n_excluded,
            dataset.ledger.n_input,
        )

    def test_tier_b_is_excluded_by_default_and_counted(self):
        record = synthetic.good_record(
            outcome_id="median_lifespan", outcome_tier="B"
        )
        dataset = build_dataset([record])
        self.assertEqual(len(dataset), 0)
        self.assertEqual(
            dataset.ledger.primary[ExclusionReason.outcome_tier_b], 1
        )
        self.assertTrue(dataset.ledger.tier_b_excluded)

    def test_tier_b_can_be_measured_but_the_ledger_records_the_switch(self):
        record = synthetic.good_record(
            outcome_id="median_lifespan", outcome_tier="B"
        )
        dataset = build_dataset([record], exclude_tier_b=False)
        self.assertEqual(len(dataset), 1)
        self.assertFalse(dataset.ledger.to_dict()["tier_b_excluded"])

    def test_unclear_valence_is_kept_as_a_feature_value(self):
        record = synthetic.good_record(valence="unclear")
        self.assertEqual(exclusion_reasons(record), [])
        dataset = build_dataset([record])
        self.assertEqual(len(dataset), 1)
        self.assertEqual(dataset.examples[0].features.valence, "unclear")

    def test_every_reason_code_is_declared_with_an_explanation_and_a_clause(self):
        for code, (why, clause) in EXCLUSION_REASONS.items():
            self.assertTrue(why.strip(), code)
            self.assertTrue(clause.strip(), code)

    def test_an_undeclared_reason_cannot_enter_the_ledger(self):
        ledger = ExclusionLedger()
        with self.assertRaises(ValueError):
            ledger.add("rid", ["invented_reason"])
        with self.assertRaises(ValueError):
            ledger.add("rid", [])

    def test_co_occurring_reasons_are_both_counted(self):
        record = synthetic.good_record(direction="not_reported", outcome_id=None)
        dataset = build_dataset([record])
        co = dataset.ledger.co_occurring
        self.assertEqual(co[ExclusionReason.direction_not_reported], 1)
        self.assertEqual(co[ExclusionReason.outcome_id_absent], 1)
        self.assertEqual(dataset.ledger.n_excluded, 1)

    def test_the_report_dict_states_the_excluded_fraction(self):
        dataset = build_dataset(
            [
                synthetic.good_record(direction="increase", source_locator="a"),
                synthetic.good_record(direction="not_reported", source_locator="b"),
            ]
        )
        self.assertAlmostEqual(
            dataset.ledger.to_dict()["excluded_fraction"], 0.5
        )


class TestEmptyDataset(unittest.TestCase):
    def test_an_empty_input_is_an_empty_dataset_not_an_error(self):
        dataset = build_dataset([])
        self.assertTrue(dataset.is_empty)
        self.assertEqual(len(dataset), 0)
        self.assertEqual(dataset.label_counts(), {label: 0 for label in LABELS})
        self.assertEqual(dataset.ledger.n_input, 0)

    def test_describe_works_on_an_empty_dataset(self):
        description = build_dataset([]).describe()
        self.assertEqual(description["n_examples"], 0)
        for name in GROUPINGS:
            self.assertEqual(description["groupings"][name]["n_groups"], 0)

    def test_a_missing_file_is_reported_not_raised(self):
        dataset = load_dataset("/nonexistent/path/records.jsonl")
        self.assertTrue(dataset.is_empty)
        self.assertEqual(len(dataset.load_errors), 1)
        self.assertIn("no such file", dataset.load_errors[0])

    def test_an_empty_file_loads_as_an_empty_dataset(self):
        import os
        import tempfile

        directory = tempfile.mkdtemp()
        path = os.path.join(directory, "records.jsonl")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("")
        dataset = load_dataset(path)
        self.assertTrue(dataset.is_empty)
        self.assertEqual(dataset.load_errors, [])
        self.assertEqual(dataset.source, path)


class TestGroupings(unittest.TestCase):
    def test_every_default_grouping_is_declared(self):
        for name in DEFAULT_GROUPINGS:
            self.assertIn(name, GROUPINGS)
            self.assertTrue(GROUPINGS[name].schema_clause)

    def test_only_the_record_grouping_is_marked_leakage_prone(self):
        prone = {n for n, g in GROUPINGS.items() if g.leakage_prone}
        self.assertEqual(prone, {"record"})

    def test_group_keys_separate_studies_labs_and_populations(self):
        first = synthetic.good_record()
        other = synthetic.other_lab_record()
        dataset = build_dataset([first, other])
        self.assertEqual(len(set(dataset.group_values("study_id"))), 2)
        self.assertEqual(len(set(dataset.group_values("lab"))), 2)
        self.assertEqual(len(set(dataset.group_values("population"))), 1)

    def test_a_record_with_no_lab_identity_cannot_be_grouped_by_lab(self):
        record = synthetic.good_record(
            last_author=None, last_author_affiliation=None
        )
        dataset = build_dataset([record])
        self.assertEqual(dataset.group_values("lab"), [None])
        coverage = dataset.grouping_coverage("lab")
        self.assertEqual(coverage["n_rows_ungroupable"], 1)
        self.assertEqual(coverage["n_groups"], 0)

    def test_an_undeclared_grouping_is_refused(self):
        dataset = build_dataset([synthetic.good_record()])
        with self.assertRaises(KeyError):
            dataset.group_values("by_vibes")
        with self.assertRaises(KeyError):
            dataset.examples[0].group("by_vibes")

    def test_duplicate_record_ids_make_the_index_refuse_rather_than_guess(self):
        record = synthetic.good_record()
        dataset = build_dataset([record, record])
        with self.assertRaises(ValueError):
            dataset.record_index()


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
