"""Unit consistency against a stubbed outcomes vocabulary.

The real vocabulary files are owned by another component. This suite must
pass whether or not they exist, so every test injects a stub loader.
"""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.qc import QCSeverity, check_units, normalize_unit


def codes(findings):
    return [f.code for f in findings]


class TestNormalizeUnit(unittest.TestCase):
    def test_case_and_whitespace(self):
        self.assertEqual(normalize_unit("  MG/DL "), "mg/dl")
        self.assertEqual(normalize_unit("mg / dL"), "mg/dl")
        self.assertEqual(normalize_unit(None), "")

    def test_aliases(self):
        self.assertEqual(normalize_unit("percent"), "%")
        self.assertEqual(normalize_unit("pct"), "%")
        self.assertEqual(normalize_unit("°C"), "degC")


class TestUnitVariants(unittest.TestCase):
    def test_parenthetical_qualifier_is_optional(self):
        from syn4rch.qc.units import unit_variants

        self.assertEqual(
            unit_variants("degC (BAT temperature)"),
            # "degc" is alias-folded to the canonical spelling; both sides of
            # a comparison pass through normalize_unit, so they still match.
            {"degc (bat temperature)", "degC"},
        )
        self.assertEqual(normalize_unit("degC"), "degC")
        self.assertEqual(unit_variants("mmHg"), {"mmhg"})

    def test_annotated_vocab_unit_matches_the_bare_unit(self):
        loader = synthetic.stub_vocab_loader(
            {"syn:outcome_1": {"units": ["degC (synthetic probe)"]}}
        )
        record = synthetic.good_record(
            assay_units="degC",
            effect_size_type="absolute_difference",
            effect_size_units="degC (synthetic probe)",
        )
        self.assertEqual(check_units([record], loader=loader), [])


class TestCheckUnits(unittest.TestCase):
    def test_matching_units_pass(self):
        loader = synthetic.stub_vocab_loader()
        self.assertEqual(check_units([synthetic.good_record()], loader=loader), [])

    def test_assay_unit_mismatch_is_a_warning(self):
        loader = synthetic.stub_vocab_loader()
        record = synthetic.good_record(assay_units="furlongs")
        findings = check_units([record], loader=loader)
        self.assertEqual(codes(findings), ["QC_UNIT_MISMATCH"])
        self.assertEqual(findings[0].severity, QCSeverity.WARNING)
        self.assertEqual(findings[0].details["found"], "furlongs")
        self.assertEqual(findings[0].details["expected"], ["syn_unit"])

    def test_dimensionful_effect_size_is_compared_to_the_outcome_unit(self):
        loader = synthetic.stub_vocab_loader()
        bad = synthetic.good_record(
            effect_size_type="absolute_difference", effect_size_units="furlongs"
        )
        self.assertEqual(codes(check_units([bad], loader=loader)), ["QC_UNIT_MISMATCH"])
        good = synthetic.good_record(
            effect_size_type="absolute_difference", effect_size_units="syn_unit"
        )
        self.assertEqual(check_units([good], loader=loader), [])

    def test_percent_change_is_not_compared_to_the_physical_unit(self):
        """A percent change is in percent, not in the readout's own unit."""
        loader = synthetic.stub_vocab_loader()
        record = synthetic.good_record(
            effect_size_type="percent_change", effect_size_units="%"
        )
        self.assertEqual(check_units([record], loader=loader), [])

    def test_dimensionless_type_with_a_physical_unit_is_flagged(self):
        loader = synthetic.stub_vocab_loader()
        record = synthetic.good_record(
            effect_size_type="percent_change", effect_size_units="mg/dL"
        )
        self.assertEqual(
            codes(check_units([record], loader=loader)),
            ["QC_UNIT_DIMENSIONLESS_MISMATCH"],
        )

    def test_units_without_an_effect_size_type_cannot_be_interpreted(self):
        loader = synthetic.stub_vocab_loader()
        record = synthetic.good_record(effect_size_type=None, effect_size_units="syn_unit")
        findings = check_units([record], loader=loader)
        self.assertEqual(codes(findings), ["QC_UNIT_UNCHECKABLE"])
        self.assertEqual(findings[0].severity, QCSeverity.INFO)

    def test_vocab_entry_may_declare_several_units(self):
        loader = synthetic.stub_vocab_loader(
            {"syn:outcome_1": {"units": ["mg/dL", "mmol/L"]}}
        )
        ok = synthetic.good_record(
            assay_units="mmol/L",
            effect_size_type="absolute_difference",
            effect_size_units="mg/dL",
        )
        self.assertEqual(check_units([ok], loader=loader), [])

    def test_unknown_outcome_id_is_a_warning(self):
        loader = synthetic.stub_vocab_loader({"syn:outcome_other": {"units": "syn_unit"}})
        findings = check_units([synthetic.good_record()], loader=loader)
        self.assertEqual(codes(findings), ["QC_OUTCOME_ID_NOT_IN_VOCAB"])
        self.assertEqual(findings[0].severity, QCSeverity.WARNING)

    def test_outcome_without_declared_units_is_info_only(self):
        loader = synthetic.stub_vocab_loader({"syn:outcome_1": {"label": "x"}})
        findings = check_units([synthetic.good_record()], loader=loader)
        self.assertEqual(codes(findings), ["QC_OUTCOME_UNITS_UNDECLARED"])
        self.assertEqual(findings[0].severity, QCSeverity.INFO)

    def test_missing_outcome_id_is_left_to_schema_validation(self):
        loader = synthetic.stub_vocab_loader()
        record = synthetic.good_record(outcome_id=None)
        self.assertEqual(check_units([record], loader=loader), [])

    def test_unavailable_vocabulary_degrades_to_info(self):
        def broken_loader(kind):
            raise FileNotFoundError("data/vocab/outcomes.json")

        findings = check_units([synthetic.good_record()], loader=broken_loader)
        self.assertEqual(codes(findings), ["QC_VOCAB_UNAVAILABLE"])
        self.assertEqual(findings[0].severity, QCSeverity.INFO)
        self.assertIn("tooling gap", findings[0].message)

    def test_loader_returning_a_non_mapping_is_handled(self):
        findings = check_units([synthetic.good_record()], loader=lambda kind: ["nope"])
        self.assertEqual(codes(findings), ["QC_VOCAB_UNAVAILABLE"])

    def test_check_never_modifies_a_record(self):
        loader = synthetic.stub_vocab_loader()
        record = synthetic.good_record(assay_units="furlongs")
        before = record.to_dict()
        check_units([record], loader=loader)
        self.assertEqual(record.to_dict(), before)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
