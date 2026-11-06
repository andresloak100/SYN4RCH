"""Duplicate, locator-collision and non-independence detection."""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.qc import (
    QCSeverity,
    affiliation_tokens,
    find_exact_duplicates,
    find_lab_identity_collisions,
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

    def test_same_author_different_affiliation_is_the_same_lab(self):
        """D-056. This test asserted the opposite, and the assertion was the defect.

        While the affiliation was part of the lab identity this check could not
        fire on the real corpus at all: four labs had published eleven papers
        and every one of those papers spelled its own department differently,
        so no lab ever appeared in two studies and the detector's silence read
        as a clean result.
        """
        a = synthetic.good_record()
        b = synthetic.same_lab_other_study_record(
            last_author_affiliation="A different synthetic institute"
        )
        self.assertEqual(
            codes(find_possible_non_independence([a, b])),
            ["QC_POSSIBLE_NON_INDEPENDENCE"],
        )

    def test_the_message_names_the_claim_components_it_prints(self):
        """The components are named, not positional.

        They were printed positionally under the wrong names until D-056:
        claim_key stopped carrying `valence` when v0.2.0 added the perturbation
        axis, so the message had been labelling the perturbation class
        "valence", the polarity "outcome" and the target "direction". Nothing
        noticed, because the check never fired.
        """
        found = find_possible_non_independence(
            [
                synthetic.good_record(),
                synthetic.same_lab_other_study_record(
                    last_author_affiliation="A different synthetic institute"
                ),
            ]
        )
        self.assertEqual(len(found), 1)
        message = found[0].message
        self.assertIn("perturbation_class=circuit_activity", message)
        self.assertIn("direction=increase", message)
        self.assertNotIn("valence=", message)
        claim = found[0].details["claim"]
        self.assertEqual(claim["perturbation_class"], "circuit_activity")
        self.assertEqual(found[0].details["lab"], "fixture, b.")

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


class TestLabIdentityCollisions(unittest.TestCase):
    """The homonym check D-056 added when it dropped the affiliation.

    `lab_key` is now the author alone, which can merge two researchers who
    share a name. The merge costs test rows in every leave-one-lab-out split
    and reports two independent results as one lab replicating itself, so it
    is detected rather than assumed away.
    """

    def test_one_author_at_two_unrelated_places_warns(self):
        found = find_lab_identity_collisions(
            [
                synthetic.good_record(
                    last_author_affiliation=(
                        "Department of Physiology, Kyoto Prefectural "
                        "University, Kyoto, Japan"
                    )
                ),
                synthetic.same_lab_other_study_record(
                    last_author_affiliation=(
                        "Rowett Institute, Aberdeen, Scotland"
                    )
                ),
            ]
        )
        self.assertEqual(codes(found), ["QC_POSSIBLE_LAB_HOMONYM"])
        self.assertEqual(found[0].severity, QCSeverity.WARNING)
        self.assertEqual(found[0].field_name, "last_author_affiliation")
        self.assertEqual(found[0].details["lab"], "fixture, b.")

    def test_the_same_department_typed_differently_does_not_warn(self):
        """The four real cases from D-056 are all of this shape."""
        found = find_lab_identity_collisions(
            [
                synthetic.good_record(
                    last_author_affiliation=(
                        "Max Planck Institute for Metabolism Research, "
                        "Department of Neuronal Control of Metabolism, "
                        "Gleueler Str. 50, 50931 Cologne, Germany"
                    )
                ),
                synthetic.same_lab_other_study_record(
                    last_author_affiliation=(
                        "Neuronal Control of Metabolism Group, Max Planck "
                        "Institute for Metabolism Research, Cologne, Germany."
                    )
                ),
            ]
        )
        self.assertEqual(found, [])

    def test_a_shared_generic_word_is_not_a_shared_place(self):
        """Two unrelated labs are both a department at a university somewhere."""
        found = find_lab_identity_collisions(
            [
                synthetic.good_record(
                    last_author_affiliation=(
                        "Department of Medicine, University of Padua, Italy"
                    )
                ),
                synthetic.same_lab_other_study_record(
                    last_author_affiliation=(
                        "Department of Medicine, University of Otago, "
                        "New Zealand"
                    )
                ),
            ]
        )
        self.assertEqual(codes(found), ["QC_POSSIBLE_LAB_HOMONYM"])

    def test_a_blank_affiliation_is_not_evidence_of_two_labs(self):
        found = find_lab_identity_collisions(
            [
                synthetic.good_record(),
                synthetic.same_lab_other_study_record(
                    last_author_affiliation=None
                ),
            ]
        )
        self.assertEqual(found, [])

    def test_one_affiliation_cannot_collide_with_itself(self):
        found = find_lab_identity_collisions(
            [synthetic.good_record(), synthetic.same_lab_other_study_record()]
        )
        self.assertEqual(found, [])

    def test_records_without_a_last_author_are_skipped(self):
        found = find_lab_identity_collisions(
            [
                synthetic.good_record(
                    last_author=None, last_author_affiliation="Place A"
                ),
                synthetic.same_lab_other_study_record(
                    last_author=None, last_author_affiliation="Place B"
                ),
            ]
        )
        self.assertEqual(found, [])

    def test_tokens_drop_generic_words_and_short_ones(self):
        tokens = affiliation_tokens(
            "Department of Neuroscience, Garvan Institute of Medical "
            "Research, NSW 2010, Australia"
        )
        self.assertIn("garvan", tokens)
        self.assertIn("2010", tokens)
        self.assertNotIn("department", tokens)
        self.assertNotIn("institute", tokens)
        self.assertNotIn("australia", tokens)
        self.assertNotIn("nsw", tokens)  # three characters


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
