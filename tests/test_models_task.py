"""The task definition: features, labels, the exclusion ledger, groupings."""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.models.task import (
    ACCEPTED_SCHEMA_VERSIONS,
    CONDITION_FEATURE_FIELDS,
    DEFAULT_GROUPINGS,
    EXCLUSION_REASONS,
    FEATURE_FIELDS,
    FEATURE_FIELDS_V1,
    GROUPINGS,
    LABELS,
    NULL_FEATURE,
    PERTURBATION_AXIS,
    PERTURBATION_AXIS_FIELDS,
    TASK,
    UNKNOWN_CONDITION,
    ExclusionLedger,
    ExclusionReason,
    Features,
    PolarityWithoutClassError,
    build_dataset,
    exclusion_reasons,
    expand_cell_fields,
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

    def test_features_are_exactly_the_eleven_the_task_names(self):
        self.assertEqual(
            TASK.feature_fields,
            (
                "population_id",
                "valence",
                "outcome_id",
                "species",
                "acuity",
                "perturbation_class",
                "perturbation_polarity",
                "perturbation_target",
                "sex",
                "diet",
                "challenge",
            ),
        )

    def test_the_statement_names_the_feature_tuple_exactly(self):
        # Two reports must never disagree about what was being predicted, so
        # the prose statement and the field list are checked against each
        # other rather than maintained in parallel.
        for name in TASK.feature_fields:
            self.assertIn(name, TASK.statement, name)
        self.assertIn("predict direction in {increase, decrease, no_change}",
                      TASK.statement)
        self.assertEqual(TASK.version, "0.2.0")
        self.assertEqual(TASK.to_dict()["feature_fields"], list(FEATURE_FIELDS))

    def test_valence_is_kept_so_the_previous_benchmark_stays_comparable(self):
        self.assertIn("valence", TASK.feature_fields)
        self.assertEqual(TASK.previous_feature_fields, FEATURE_FIELDS_V1)
        self.assertEqual(
            TASK.to_dict()["features_added_in_this_version"],
            [
                "perturbation_class",
                "perturbation_polarity",
                "perturbation_target",
                "sex",
                "diet",
                "challenge",
            ],
        )

    def test_the_task_refuses_a_polarity_without_its_class(self):
        refused = " ".join(TASK.refuses).lower()
        self.assertIn("perturbation_polarity without perturbation_class", refused)
        self.assertIn("unknown", refused)


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


class TestPolarityNeverTravelsWithoutItsClass(unittest.TestCase):
    """The D-019 guard. A polarity means different physical things in
    different classes, so nothing in this package may key on it alone."""

    def features(self, **overrides):
        return extract_features(synthetic.molecular_record(**overrides))

    def test_a_model_cannot_key_a_cell_on_polarity_alone(self):
        features = self.features()
        with self.assertRaises(PolarityWithoutClassError):
            features.cell("perturbation_polarity")
        with self.assertRaises(PolarityWithoutClassError):
            features.cell("population_key", "perturbation_polarity", "outcome_id")

    def test_the_class_makes_the_same_key_legitimate(self):
        features = self.features()
        self.assertEqual(
            features.cell("perturbation_class", "perturbation_polarity"),
            ("molecular", "decrease"),
        )

    def test_the_composite_expands_to_class_polarity_and_target(self):
        self.assertEqual(
            expand_cell_fields("population_key", PERTURBATION_AXIS, "outcome_id"),
            ("population_key",) + PERTURBATION_AXIS_FIELDS + ("outcome_id",),
        )
        features = self.features()
        self.assertEqual(
            features.cell(PERTURBATION_AXIS),
            ("molecular", "decrease", "syn:synthgene1"),
        )

    def test_a_feature_tuple_with_a_polarity_and_no_class_cannot_exist(self):
        with self.assertRaises(PolarityWithoutClassError):
            Features(
                population_id="p",
                population_key="p",
                valence=None,
                outcome_id="o",
                species=None,
                acuity=None,
                perturbation_class=None,
                perturbation_polarity="decrease",
            )

    def test_a_record_whose_polarity_has_no_class_raises_at_build_time(self):
        # method=None leaves the class underivable; the stored polarity then
        # has nothing to qualify it, and the dataset must refuse rather than
        # quietly carry a meaningless sign.
        record = synthetic.good_record(
            method=None, perturbation_class=None, perturbation_polarity="decrease"
        )
        with self.assertRaises(PolarityWithoutClassError):
            extract_features(record)
        with self.assertRaises(PolarityWithoutClassError):
            build_dataset([record])

    def test_the_same_polarity_in_two_classes_occupies_two_cells(self):
        circuit = extract_features(
            synthetic.good_record(
                method="chemogenetic",
                valence="inhibition",
                perturbation_class="circuit_activity",
                perturbation_polarity="decrease",
            )
        )
        molecular = self.features()
        self.assertEqual(circuit.perturbation_polarity, "decrease")
        self.assertEqual(molecular.perturbation_polarity, "decrease")
        self.assertNotEqual(
            circuit.cell(PERTURBATION_AXIS), molecular.cell(PERTURBATION_AXIS)
        )

    def test_the_axis_components_come_from_qc_keys_not_a_reimplementation(self):
        from syn4rch.qc.keys import perturbation_key

        record = synthetic.molecular_record()
        self.assertEqual(
            extract_features(record).cell(PERTURBATION_AXIS),
            perturbation_key(record),
        )

    def test_a_pre_v0_2_0_record_keys_the_same_as_a_migrated_one(self):
        # qc.keys derives the class and the polarity from method and valence
        # when the stored fields are absent, so the models layer is correct on
        # both sides of the migration.
        stored = synthetic.molecular_record()
        derived = synthetic.molecular_record(
            perturbation_class=None, perturbation_polarity=None
        )
        self.assertEqual(
            extract_features(stored).cell(PERTURBATION_AXIS),
            extract_features(derived).cell(PERTURBATION_AXIS),
        )


class TestConditionFeaturesDistinguishUnknownFromAValue(unittest.TestCase):
    """A condition feature must not make an absence look like a value."""

    def test_none_and_not_reported_pool_into_one_unknown_token(self):
        absent = extract_features(synthetic.good_record(sex=None))
        checked = extract_features(synthetic.good_record(sex="not_reported"))
        self.assertEqual(absent.sex, UNKNOWN_CONDITION)
        self.assertEqual(checked.sex, UNKNOWN_CONDITION)
        self.assertEqual(absent.cell("sex"), checked.cell("sex"))

    def test_unknown_is_distinguishable_from_every_stated_value(self):
        male = extract_features(synthetic.good_record(sex="male"))
        unknown = extract_features(synthetic.good_record(sex=None))
        self.assertEqual(male.sex, "male")
        self.assertNotEqual(male.cell("sex"), unknown.cell("sex"))
        self.assertNotEqual(unknown.sex, NULL_FEATURE)
        self.assertNotIn(UNKNOWN_CONDITION, ("male", "female", "both"))

    def test_the_convention_is_the_one_in_qc_keys(self):
        from syn4rch.qc.keys import UNKNOWN, condition_key

        record = synthetic.good_record(sex=None, diet="high_fat", challenge=None)
        features = extract_features(record)
        sex, diet, _age_group, challenge = condition_key(record)
        self.assertEqual(UNKNOWN_CONDITION, UNKNOWN)
        self.assertEqual((features.sex, features.diet, features.challenge),
                         (sex, diet, challenge))

    def test_condition_specificity_travels_with_the_example(self):
        from syn4rch.qc.keys import condition_specificity

        record = synthetic.good_record(sex="male", diet=None, challenge=None)
        features = extract_features(record)
        self.assertEqual(features.condition_specificity, condition_specificity(record))
        self.assertEqual(features.condition_specificity, 2)  # sex + age_group

    def test_specificity_is_carried_into_the_dataset_summary(self):
        dataset = build_dataset(
            [
                synthetic.good_record(source_locator="a"),
                synthetic.good_record(
                    sex=None, diet=None, age_group=None, challenge=None,
                    source_locator="b",
                ),
            ]
        )
        summary = dataset.describe()["conditions"]
        self.assertEqual(summary["slots"], ["sex", "diet", "age_group", "challenge"])
        self.assertEqual(summary["feature_slots"], list(CONDITION_FEATURE_FIELDS))
        self.assertEqual(summary["by_slots_known"]["0"], 1)
        self.assertEqual(summary["by_slots_known"]["4"], 1)
        self.assertEqual(summary["unknown_token"], UNKNOWN_CONDITION)

    def test_specificity_is_not_a_feature(self):
        features = extract_features(synthetic.good_record())
        with self.assertRaises(KeyError):
            features.cell("condition_specificity")


class TestAxisSummary(unittest.TestCase):
    def test_the_summary_counts_classes_polarities_and_targets(self):
        dataset = build_dataset(
            [
                synthetic.good_record(source_locator="a"),
                synthetic.molecular_record(source_locator="b"),
                synthetic.molecular_record(
                    method="conditional_overexpression",
                    perturbation_polarity="increase",
                    target_raw=None,
                    target_id=None,
                    source_locator="c",
                ),
            ]
        )
        summary = dataset.describe()["perturbation_axis"]
        self.assertEqual(summary["class_counts"],
                         {"circuit_activity": 1, "molecular": 2})
        self.assertEqual(summary["class_polarity_counts"]["molecular / decrease"], 1)
        self.assertEqual(summary["n_with_signed_polarity"], 3)
        self.assertEqual(summary["n_with_known_target"], 1)


class TestSchemaVersionsAreReadAndReported(unittest.TestCase):
    def test_both_accepted_versions_are_declared_current_first(self):
        from syn4rch.schema.record import SCHEMA_VERSION

        self.assertEqual(ACCEPTED_SCHEMA_VERSIONS[0], SCHEMA_VERSION)
        self.assertIn("0.1.0", ACCEPTED_SCHEMA_VERSIONS)

    def test_the_versions_actually_read_are_counted_and_described(self):
        dataset = build_dataset([synthetic.good_record()])
        self.assertEqual(
            dataset.describe()["schema_versions_read"],
            {ACCEPTED_SCHEMA_VERSIONS[0]: 1},
        )

    def test_an_unaccepted_version_stops_the_read_with_a_message(self):
        import json
        import os
        import tempfile

        from syn4rch.schema.io import dumps_record

        directory = tempfile.mkdtemp()
        path = os.path.join(directory, "records.jsonl")
        line = json.loads(dumps_record(synthetic.good_record()))
        line["schema_version"] = "9.9.9"
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(dumps_record(synthetic.good_record()) + "\n")
            handle.write(json.dumps(line) + "\n")
        dataset = load_dataset(path)
        self.assertEqual(len(dataset), 1)
        self.assertEqual(len(dataset.load_errors), 1)
        self.assertIn("9.9.9", dataset.load_errors[0])


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
