"""Missing-condition detection and whole-report assembly."""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.qc import (
    QCSeverity,
    condition_coverage,
    find_missing_conditions,
    missing_conditions,
    render_text,
    run_qc,
)


class TestMissingConditions(unittest.TestCase):
    def test_a_fully_reported_record_is_missing_nothing(self):
        self.assertEqual(missing_conditions(synthetic.good_record()), [])
        self.assertEqual(find_missing_conditions([synthetic.good_record()]), [])

    def test_each_condition_is_detected(self):
        cases = (
            ({"control_types": []}, "control_types"),
            ({"control_types": ["none_reported"]}, "control_types"),
            ({"n_animals": None, "n_treatment": None, "n_control": None}, "n"),
            ({"species": None}, "species"),
            ({"valence": None}, "valence"),
        )
        for override, expected in cases:
            self.assertIn(
                expected, missing_conditions(synthetic.good_record(**override)), override
            )

    def test_partial_n_counts_as_reported(self):
        record = synthetic.good_record(n_treatment=8, n_control=None, n_animals=None)
        self.assertNotIn("n", missing_conditions(record))

    def test_unclear_valence_is_not_missing(self):
        record = synthetic.good_record(valence="unclear")
        self.assertNotIn("valence", missing_conditions(record))

    def test_finding_lists_everything_absent(self):
        record = synthetic.good_record(
            control_types=[],
            species=None,
            valence=None,
            n_animals=None,
            n_treatment=None,
            n_control=None,
        )
        findings = find_missing_conditions([record])
        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0].severity, QCSeverity.WARNING)
        self.assertEqual(
            findings[0].details["missing"],
            ["control_types", "n", "species", "valence"],
        )
        self.assertIn("not inferred", findings[0].message)

    def test_coverage_counts(self):
        records = [
            synthetic.good_record(),
            synthetic.good_record(species=None, source_locator="Fig S1B (synthetic)"),
        ]
        coverage = condition_coverage(records)
        self.assertEqual(coverage["_total_records"], 2)
        self.assertEqual(coverage["species"], 1)
        self.assertEqual(coverage["valence"], 2)


class TestRunQc(unittest.TestCase):
    def setUp(self):
        self.loader = synthetic.stub_vocab_loader()

    def test_clean_records_yield_no_hard_findings(self):
        records = [synthetic.good_record(), synthetic.other_lab_record()]
        report = run_qc(records, vocab_loader=self.loader)
        self.assertEqual(report.hard, [], [str(f) for f in report.hard])
        self.assertTrue(report.ok)
        self.assertEqual(report.n_records, 2)
        self.assertEqual(report.stats["n_studies"], 2)

    def test_hard_violations_make_the_report_not_ok(self):
        records = [synthetic.good_record(support_quote=None)]
        report = run_qc(records, vocab_loader=self.loader)
        self.assertFalse(report.ok)
        self.assertIn(
            "A001_MISSING_SUPPORT_QUOTE", {f.code for f in report.hard}
        )
        self.assertEqual(report.blocking_record_ids, [records[0].record_id])

    def test_conflicts_do_not_make_the_report_not_ok(self):
        up = synthetic.good_record()
        down = synthetic.other_lab_record(
            direction="decrease",
            effect_size_value=-30.0,
            support_quote="Synthetic quote: the readout decreased by 30%.",
        )
        report = run_qc([up, down], vocab_loader=self.loader)
        self.assertTrue(report.ok, "a conflict must not block admission")
        self.assertIn("QC_CONFLICT_OPPOSED", {f.code for f in report.info})
        self.assertEqual(report.stats["n_conflicting_records"], 2)
        self.assertEqual(up.conflicts_with, [down.record_id])

    def test_annotate_false_leaves_conflicts_with_untouched(self):
        up = synthetic.good_record()
        down = synthetic.other_lab_record(
            direction="decrease",
            effect_size_value=-30.0,
            support_quote="Synthetic quote: the readout decreased by 30%.",
        )
        report = run_qc([up, down], annotate=False, vocab_loader=self.loader)
        self.assertEqual(up.conflicts_with, [])
        self.assertIn("QC_CONFLICT_OPPOSED", {f.code for f in report.info})

    def test_every_severity_bucket_is_separate(self):
        records = [
            synthetic.good_record(support_quote=None, effect_size_value=None,
                                  effect_size_type=None, effect_size_units=None,
                                  confidence="moderate"),
            synthetic.other_lab_record(control_types=[], has_actuator_control=False,
                                       confidence="moderate"),
            synthetic.same_lab_other_study_record(),
        ]
        report = run_qc(records, vocab_loader=self.loader)
        self.assertTrue(report.hard)
        self.assertTrue(report.warnings)
        self.assertTrue(report.info)
        severities = {f.severity for f in report.hard}
        self.assertEqual(severities, {QCSeverity.HARD})
        self.assertEqual({f.severity for f in report.warnings}, {QCSeverity.WARNING})
        self.assertEqual({f.severity for f in report.info}, {QCSeverity.INFO})

    def test_report_serialises_to_json_safe_types(self):
        import json

        report = run_qc(
            [synthetic.good_record(), synthetic.other_lab_record(direction="decrease",
             effect_size_value=-30.0,
             support_quote="Synthetic quote: the readout decreased by 30%.")],
            vocab_loader=self.loader,
        )
        text = json.dumps(report.to_dict())
        self.assertIn("counts_by_code", text)
        reloaded = json.loads(text)
        self.assertEqual(reloaded["n_records"], 2)

    def test_text_report_is_readable_and_labels_conflicts_as_data(self):
        up = synthetic.good_record()
        down = synthetic.other_lab_record(
            direction="decrease",
            effect_size_value=-30.0,
            support_quote="Synthetic quote: the readout decreased by 30%.",
        )
        text = render_text(run_qc([up, down], vocab_loader=self.loader))
        self.assertIn("SYN4RCH QC report", text)
        self.assertIn("HARD schema violations", text)
        self.assertIn("these are DATA", text)
        self.assertIn("VERDICT", text)

    def test_vocab_absence_is_reported_but_not_fatal(self):
        def broken(kind):
            raise FileNotFoundError("data/vocab/outcomes.json")

        report = run_qc([synthetic.good_record()], vocab_loader=broken)
        self.assertTrue(report.ok)
        self.assertIn("QC_VOCAB_UNAVAILABLE", {f.code for f in report.info})

    def test_no_record_is_dropped_by_the_qc_pass(self):
        records = [
            synthetic.good_record(),
            synthetic.good_record(),  # exact duplicate
            synthetic.other_lab_record(
                direction="decrease",
                effect_size_value=-30.0,
                support_quote="Synthetic quote: the readout decreased by 30%.",
            ),
        ]
        report = run_qc(records, vocab_loader=self.loader)
        self.assertEqual(report.n_records, 3)
        self.assertEqual(len(records), 3)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
