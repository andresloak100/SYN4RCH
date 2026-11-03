"""Duplicate, locator-collision and non-independence detection."""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.qc import (
    QCSeverity,
    find_exact_duplicates,
    find_locator_collisions,
    find_possible_non_independence,
)


def codes(findings):
    return [f.code for f in findings]


class TestExactDuplicates(unittest.TestCase):
    def test_identical_records_are_flagged_once(self):
        records = [synthetic.good_record(), synthetic.good_record()]
        findings = find_exact_duplicates(records)
        self.assertEqual(codes(findings), ["QC_DUP_EXACT"])
        self.assertEqual(findings[0].severity, QCSeverity.WARNING)
        self.assertEqual(findings[0].details["copies"], 2)

    def test_three_copies_report_the_count(self):
        findings = find_exact_duplicates([synthetic.good_record() for _ in range(3)])
        self.assertEqual(findings[0].details["copies"], 3)

    def test_distinct_records_are_not_flagged(self):
        records = [synthetic.good_record(), synthetic.other_lab_record()]
        self.assertEqual(find_exact_duplicates(records), [])

    def test_same_id_divergent_content_is_a_separate_code(self):
        a = synthetic.good_record()
        b = synthetic.good_record(
            effect_size_value=40.0,
            support_quote="Synthetic quote: a 40% increase was observed.",
        )
        self.assertEqual(a.record_id, b.record_id)  # effect size is not hashed
        findings = find_exact_duplicates([a, b])
        self.assertEqual(codes(findings), ["QC_DUP_ID_CONTENT_DIVERGENT"])

    def test_conflicts_with_does_not_affect_duplicate_detection(self):
        a = synthetic.good_record()
        b = synthetic.good_record()
        b.conflicts_with = ["ffffffffffffffff"]
        self.assertEqual(codes(find_exact_duplicates([a, b])), ["QC_DUP_EXACT"])

    def test_nothing_is_removed(self):
        records = [synthetic.good_record(), synthetic.good_record()]
        find_exact_duplicates(records)
        self.assertEqual(len(records), 2)


class TestLocatorCollisions(unittest.TestCase):
    def test_same_study_same_locator_different_values(self):
        a = synthetic.good_record(outcome_raw="synthetic readout one")
        b = synthetic.good_record(
            outcome_raw="synthetic readout one",
            actuator="synthetic-actuator-2",  # different record_id, same claim+locator
            effect_size_value=40.0,
            support_quote="Synthetic quote: a 40% increase was observed.",
        )
        self.assertNotEqual(a.record_id, b.record_id)
        findings = find_locator_collisions([a, b])
        self.assertEqual(codes(findings), ["QC_LOCATOR_COLLISION"])
        self.assertEqual(sorted(findings[0].record_ids), sorted([a.record_id, b.record_id]))

    def test_different_locators_do_not_collide(self):
        a = synthetic.good_record(source_locator="Fig S1A (synthetic)")
        b = synthetic.good_record(
            source_locator="Fig S1B (synthetic)", actuator="synthetic-actuator-2"
        )
        self.assertEqual(find_locator_collisions([a, b]), [])

    def test_different_studies_sharing_a_locator_string_do_not_collide(self):
        a = synthetic.good_record()
        b = synthetic.other_lab_record(source_locator="Fig S1A (synthetic)")
        self.assertEqual(find_locator_collisions([a, b]), [])

    def test_locator_comparison_ignores_case_and_spacing(self):
        a = synthetic.good_record(source_locator="Fig S1A (synthetic)")
        b = synthetic.good_record(
            source_locator="  fig  s1a   (SYNTHETIC) ",
            actuator="synthetic-actuator-2",
        )
        self.assertEqual(codes(find_locator_collisions([a, b])), ["QC_LOCATOR_COLLISION"])

    def test_exact_duplicates_are_not_double_reported_as_collisions(self):
        records = [synthetic.good_record(), synthetic.good_record()]
        self.assertEqual(find_locator_collisions(records), [])

    def test_records_without_a_locator_are_left_to_schema_validation(self):
        records = [
            synthetic.good_record(source_locator=None),
            synthetic.good_record(source_locator=None, actuator="synthetic-actuator-2"),
        ]
        self.assertEqual(find_locator_collisions(records), [])


class TestPossibleNonIndependence(unittest.TestCase):
    def test_same_lab_two_studies_same_claim_is_flagged_as_info(self):
        a = synthetic.good_record()
        b = synthetic.same_lab_other_study_record()
        findings = find_possible_non_independence([a, b])
        self.assertEqual(codes(findings), ["QC_POSSIBLE_NON_INDEPENDENCE"])
        finding = findings[0]
        self.assertEqual(finding.severity, QCSeverity.INFO)
        self.assertEqual(
            finding.details["studies"], ["synthetic:test-001", "synthetic:test-003"]
        )
        self.assertIn("not an error", finding.message)

    def test_it_is_not_an_error_and_changes_nothing(self):
        a = synthetic.good_record()
        b = synthetic.same_lab_other_study_record()
        before = (a.to_dict(), b.to_dict())
        findings = find_possible_non_independence([a, b])
        self.assertTrue(all(f.severity is QCSeverity.INFO for f in findings))
        self.assertEqual((a.to_dict(), b.to_dict()), before)

    def test_different_labs_are_not_flagged(self):
        findings = find_possible_non_independence(
            [synthetic.good_record(), synthetic.other_lab_record()]
        )
        self.assertEqual(findings, [])

    def test_same_lab_different_direction_is_not_non_independence(self):
        a = synthetic.good_record()
        b = synthetic.same_lab_other_study_record(
            direction="decrease",
            support_quote="Synthetic quote: the readout decreased by 25%.",
        )
        self.assertEqual(find_possible_non_independence([a, b]), [])

    def test_same_lab_same_study_is_not_flagged(self):
        a = synthetic.good_record()
        b = synthetic.good_record(
            actuator="synthetic-actuator-2", source_locator="Fig S1B (synthetic)"
        )
        self.assertEqual(find_possible_non_independence([a, b]), [])

    def test_records_without_lab_identity_are_skipped_not_grouped(self):
        a = synthetic.good_record(last_author=None, last_author_affiliation=None)
        b = synthetic.same_lab_other_study_record(
            last_author=None, last_author_affiliation=None
        )
        self.assertEqual(find_possible_non_independence([a, b]), [])

    def test_same_author_different_affiliation_is_a_different_lab(self):
        a = synthetic.good_record()
        b = synthetic.same_lab_other_study_record(
            last_author_affiliation="A different synthetic institute"
        )
        self.assertEqual(find_possible_non_independence([a, b]), [])

    def test_unresolved_populations_group_by_raw_name_only(self):
        """A null population_id must not group with a resolved one."""
        a = synthetic.good_record(population_id=None)
        b = synthetic.same_lab_other_study_record(population_id="syn:pop_alpha")
        self.assertEqual(find_possible_non_independence([a, b]), [])
        c = synthetic.same_lab_other_study_record(population_id=None)
        self.assertEqual(
            codes(find_possible_non_independence([a, c])),
            ["QC_POSSIBLE_NON_INDEPENDENCE"],
        )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
