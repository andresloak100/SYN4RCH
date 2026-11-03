"""Every anti-fabrication rule rejects a crafted bad record and accepts a good one."""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.schema import (
    RULES,
    Severity,
    hard_findings,
    is_admissible,
    validate,
    warning_findings,
)


def codes(findings):
    return {f.code for f in findings}


class TestGoodRecordAccepted(unittest.TestCase):
    def test_fully_reported_record_has_no_findings_at_all(self):
        findings = validate(synthetic.good_record())
        self.assertEqual(
            findings, [], "good record produced findings: %s" % [str(f) for f in findings]
        )
        self.assertTrue(is_admissible(synthetic.good_record()))

    def test_no_change_direction_is_accepted(self):
        """Null results must be storable: dropping them biases direction."""
        record = synthetic.good_record(
            direction="no_change",
            support_quote="Synthetic quote: no difference was observed (n = 8 per group).",
            effect_size_value=None,
            effect_size_type=None,
            effect_size_units=None,
            confidence="moderate",
        )
        self.assertEqual(hard_findings(validate(record)), [])

    def test_unclear_valence_is_accepted(self):
        record = synthetic.good_record(valence="unclear")
        self.assertEqual(hard_findings(validate(record)), [])


class TestAntiFabricationRules(unittest.TestCase):
    """SCHEMA.md 2 -- each violation must be HARD."""

    def assert_hard(self, code, record):
        findings = validate(record)
        self.assertIn(code, codes(hard_findings(findings)), [str(f) for f in findings])
        self.assertFalse(is_admissible(record))

    def test_A001_missing_support_quote(self):
        for value in (None, "", "   ", "\n\t"):
            self.assert_hard(
                "A001_MISSING_SUPPORT_QUOTE",
                synthetic.good_record(
                    support_quote=value, effect_size_value=None, confidence="moderate"
                ),
            )

    def test_A002_missing_source_locator(self):
        for value in (None, "", "  "):
            self.assert_hard(
                "A002_MISSING_SOURCE_LOCATOR", synthetic.good_record(source_locator=value)
            )

    def test_A003_direction_asserted_without_quote(self):
        for direction in ("increase", "decrease", "no_change", "non_monotonic"):
            self.assert_hard(
                "A003_DIRECTION_WITHOUT_QUOTE",
                synthetic.good_record(
                    direction=direction, support_quote=None, effect_size_value=None
                ),
            )

    def test_A003_not_reported_direction_needs_no_direction_quote(self):
        """direction=not_reported asserts nothing, so A003 must not fire.

        A001 still does: support_quote is mandatory for every record.
        """
        record = synthetic.good_record(
            direction="not_reported", support_quote=None, effect_size_value=None
        )
        found = codes(validate(record))
        self.assertIn("A001_MISSING_SUPPORT_QUOTE", found)
        self.assertNotIn("A003_DIRECTION_WITHOUT_QUOTE", found)

    def test_A004_numeric_effect_size_with_digitless_quote(self):
        self.assert_hard(
            "A004_EFFECT_SIZE_UNSUPPORTED",
            synthetic.good_record(
                effect_size_value=25.0,
                support_quote="Synthetic quote: the readout increased markedly.",
            ),
        )

    def test_A004_accepts_a_quote_containing_the_number(self):
        record = synthetic.good_record(
            effect_size_value=25.0,
            support_quote="Synthetic quote: a 25% increase was observed.",
        )
        self.assertNotIn("A004_EFFECT_SIZE_UNSUPPORTED", codes(validate(record)))

    def test_A004_does_not_fire_without_a_numeric_effect_size(self):
        record = synthetic.good_record(
            effect_size_value=None,
            effect_size_type=None,
            effect_size_units=None,
            support_quote="Synthetic quote: the readout increased.",
            confidence="moderate",
        )
        self.assertNotIn("A004_EFFECT_SIZE_UNSUPPORTED", codes(validate(record)))

    def test_A005_missing_outcome_id(self):
        for value in (None, "", " "):
            self.assert_hard(
                "A005_MISSING_OUTCOME_ID", synthetic.good_record(outcome_id=value)
            )

    def test_A006_non_causal_evidence_class(self):
        for value in ("functional_association", "anatomical", None):
            self.assert_hard(
                "A006_EVIDENCE_CLASS_NOT_CAUSAL",
                synthetic.good_record(evidence_class=value),
            )

    def test_A007_missing_extraction_method(self):
        self.assert_hard(
            "A007_MISSING_EXTRACTION_METHOD", synthetic.good_record(extraction_method=None)
        )

    def test_every_anti_fabrication_code_is_exercised_here(self):
        """No A-family rule may be added without a test in this file."""
        with open(__file__, "r", encoding="utf-8") as handle:
            source = handle.read()
        for code in sorted(c for c in RULES if c.startswith("A")):
            self.assertIn(code, source, "no test exercises %s" % code)


class TestStructuralAndValueRules(unittest.TestCase):
    def test_missing_identity_fields_are_hard(self):
        for override, code in (
            ({"study_id": None}, "S001_MISSING_STUDY_ID"),
            ({"population_raw": None}, "S002_MISSING_POPULATION_RAW"),
            ({"outcome_raw": None}, "S003_MISSING_OUTCOME_RAW"),
            ({"direction": None}, "S004_MISSING_DIRECTION"),
            ({"schema_version": "0.2.0"}, "S006_SCHEMA_VERSION_MISMATCH"),
        ):
            findings = hard_findings(validate(synthetic.good_record(**override)))
            self.assertIn(code, codes(findings), override)

    def test_tampered_record_id_is_hard(self):
        record = synthetic.good_record()
        record.record_id = "0000000000000000"
        self.assertIn("S005_RECORD_ID_INCONSISTENT", codes(hard_findings(validate(record))))

    def test_value_sanity(self):
        for override, code in (
            ({"p_value": 1.5}, "V001_P_VALUE_OUT_OF_RANGE"),
            ({"p_value": -0.1}, "V001_P_VALUE_OUT_OF_RANGE"),
            ({"ci_low": 5.0, "ci_high": 1.0}, "V002_CI_BOUNDS_INVERTED"),
            ({"n_treatment": -1}, "V003_NEGATIVE_COUNT"),
            ({"dispersion_value": -2.0}, "V004_NEGATIVE_DISPERSION"),
        ):
            findings = hard_findings(validate(synthetic.good_record(**override)))
            self.assertIn(code, codes(findings), override)


class TestWarningsAreDistinguishable(unittest.TestCase):
    def test_soft_findings_do_not_block_admission(self):
        record = synthetic.good_record(
            control_types=[],
            has_actuator_control=False,
            n_treatment=None,
            n_control=None,
            n_animals=None,
            population_id=None,
            full_text_available=False,
            retrieved_from="europepmc_abstract",
            confidence="low",
        )
        findings = validate(record)
        self.assertEqual(hard_findings(findings), [], [str(f) for f in findings])
        self.assertTrue(is_admissible(record))
        found = codes(warning_findings(findings))
        for code in (
            "W001_NO_CONTROL_TYPE",
            "W002_NO_SAMPLE_SIZE",
            "W003_POPULATION_ID_NULL",
            "W004_ABSTRACT_ONLY",
        ):
            self.assertIn(code, found)

    def test_none_reported_control_counts_as_no_control(self):
        record = synthetic.good_record(
            control_types=["none_reported"], has_actuator_control=False, confidence="moderate"
        )
        self.assertIn("W001_NO_CONTROL_TYPE", codes(validate(record)))

    def test_derived_field_disagreement_is_a_warning_not_an_edit(self):
        record = synthetic.good_record(age_group="aged", confidence="low")
        findings = validate(record)
        self.assertEqual(hard_findings(findings), [])
        self.assertIn("W009_AGE_GROUP_RULE_MISMATCH", codes(findings))
        self.assertIn("W007_CONFIDENCE_RULE_MISMATCH", codes(findings))
        # validate() must not have rewritten anything
        self.assertEqual(record.age_group.value, "aged")
        self.assertEqual(record.confidence.value, "low")

    def test_significance_without_a_test_is_a_warning(self):
        record = synthetic.good_record(statistical_test=None, p_value=None, significant=True)
        findings = validate(record)
        self.assertEqual(hard_findings(findings), [])
        self.assertIn("W008_SIGNIFICANCE_WITHOUT_TEST", codes(findings))

    def test_missing_lab_identity_is_a_warning(self):
        record = synthetic.good_record(last_author=None, last_author_affiliation=None)
        self.assertIn("W010_NO_LAB_IDENTITY", codes(warning_findings(validate(record))))

    def test_unreviewed_llm_extraction_is_a_warning(self):
        record = synthetic.good_record(extraction_method="llm", review_status="unreviewed")
        self.assertIn("W005_UNREVIEWED_LLM_EXTRACTION", codes(validate(record)))

    def test_severity_is_carried_on_every_finding(self):
        record = synthetic.good_record(support_quote=None, control_types=[])
        for finding in validate(record):
            self.assertIn(finding.severity, (Severity.HARD, Severity.WARNING))
            self.assertEqual(finding.severity, RULES[finding.code][0])
            self.assertTrue(finding.code)
            self.assertTrue(finding.message)
            self.assertEqual(finding.is_hard, finding.severity is Severity.HARD)

    def test_finding_carries_field_and_record_id(self):
        record = synthetic.good_record(source_locator=None)
        finding = [f for f in validate(record) if f.code == "A002_MISSING_SOURCE_LOCATOR"][0]
        self.assertEqual(finding.field, "source_locator")
        self.assertEqual(finding.record_id, record.record_id)
        self.assertEqual(finding.to_dict()["severity"], "HARD")


class TestRuleCatalogue(unittest.TestCase):
    def test_every_emitted_code_is_catalogued(self):
        records = [
            synthetic.good_record(),
            synthetic.good_record(support_quote=None, source_locator=None),
            synthetic.good_record(evidence_class=None, outcome_id=None),
            synthetic.good_record(control_types=[], n_animals=None, n_treatment=None,
                                  n_control=None, population_id=None,
                                  full_text_available=False, confidence="low"),
        ]
        for record in records:
            for finding in validate(record):
                self.assertIn(finding.code, RULES)

    def test_catalogue_severities_are_consistent_with_code_family(self):
        for code, (severity, _field, _msg) in RULES.items():
            expected = Severity.WARNING if code.startswith("W") else Severity.HARD
            self.assertIs(severity, expected, code)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
