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
    SCHEMA_VERSION,
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
            # A version that is deliberately not the current one. Spelled as
            # "not SCHEMA_VERSION" rather than as a literal, because the
            # literal that used to sit here became the *correct* version at the
            # v0.2.0 bump and the test silently stopped checking anything.
            ({"schema_version": SCHEMA_VERSION + "-not-this"},
             "S006_SCHEMA_VERSION_MISMATCH"),
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

    def test_an_author_with_no_affiliation_still_has_a_lab_identity(self):
        """W010 tracks lab_key exactly (D-056), and lab_key is the author.

        Seventeen records in the corpus carried W010 while grouping perfectly
        well, because this check asked for the affiliation that D-056 removed
        from the key. The affiliation's absence is still worth a warning -- the
        homonym check needs it -- but it is a different warning.
        """
        record = synthetic.good_record(last_author_affiliation=None)
        found = codes(warning_findings(validate(record)))
        self.assertNotIn("W010_NO_LAB_IDENTITY", found)
        self.assertIn("W018_NO_LAB_AFFILIATION", found)

    def test_a_record_with_no_author_is_not_also_told_about_its_affiliation(self):
        """One statement per record: it has no lab identity at all."""
        record = synthetic.good_record(last_author=None, last_author_affiliation=None)
        found = codes(warning_findings(validate(record)))
        self.assertIn("W010_NO_LAB_IDENTITY", found)
        self.assertNotIn("W018_NO_LAB_AFFILIATION", found)

    def test_the_two_lab_warnings_agree_with_lab_key(self):
        from syn4rch.qc.keys import lab_key

        for overrides in (
            {},
            {"last_author_affiliation": None},
            {"last_author": None, "last_author_affiliation": None},
        ):
            record = synthetic.good_record(**overrides)
            found = codes(warning_findings(validate(record)))
            self.assertEqual(
                lab_key(record) is None,
                "W010_NO_LAB_IDENTITY" in found,
                overrides,
            )

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



class V020ConditionAndPerturbationTest(unittest.TestCase):
    """The rules schema v0.2.0 added (S007, S008, W013..W015).

    A008 was a v0.2.0 rule too, about the record-level
    ``condition_source``/``condition_locator`` pair. v0.3.0 retired that pair
    and restated A008 against ``condition_evidence``, so A008 is exercised in
    :class:`V030ConditionEvidenceTest` below, next to the three hard rules and
    the warning added with it.
    """

    def codes(self, record):
        return [f.code for f in validate(record)]

    def test_a_clean_fixture_triggers_none_of_them(self):
        """Guards against the new rules firing on a well-formed record."""
        new_codes = {
            "S007_PERTURBATION_CLASS_MISMATCH",
            "S008_TARGET_ON_NON_MOLECULAR",
            "W013_PERTURBATION_POLARITY_RULE_MISMATCH",
            "W014_MOLECULAR_TARGET_UNIDENTIFIED",
            "W015_CONDITION_UNSPECIFIED",
        }
        self.assertEqual(set(self.codes(synthetic.good_record())) & new_codes, set())
        self.assertEqual(
            set(self.codes(synthetic.molecular_record())) & new_codes, set()
        )

    def test_class_disagreeing_with_method_is_refused_not_warned(self):
        bad = synthetic.good_record(perturbation_class="molecular")
        findings = [
            f for f in validate(bad)
            if f.code == "S007_PERTURBATION_CLASS_MISMATCH"
        ]
        self.assertEqual(len(findings), 1)
        self.assertTrue(findings[0].is_hard)
        self.assertFalse(is_admissible(bad))

    def test_target_on_a_circuit_perturbation_is_refused(self):
        bad = synthetic.good_record(target_raw="Synthgene1")
        self.assertIn("S008_TARGET_ON_NON_MOLECULAR", self.codes(bad))
        self.assertFalse(is_admissible(bad))

    def test_target_on_a_molecular_perturbation_is_fine(self):
        self.assertNotIn(
            "S008_TARGET_ON_NON_MOLECULAR",
            self.codes(synthetic.molecular_record()),
        )

    def test_a_wrong_stored_class_cannot_suppress_the_target_check(self):
        """Method outranks a derived field that disagrees with it."""
        bad = synthetic.good_record(
            perturbation_class="molecular", target_raw="Synthgene1"
        )
        codes = self.codes(bad)
        self.assertIn("S007_PERTURBATION_CLASS_MISMATCH", codes)
        self.assertIn("S008_TARGET_ON_NON_MOLECULAR", codes)

    def test_polarity_disagreeing_with_the_rule_warns(self):
        r = synthetic.good_record(perturbation_polarity="decrease")
        findings = [
            f for f in validate(r)
            if f.code == "W013_PERTURBATION_POLARITY_RULE_MISMATCH"
        ]
        self.assertEqual(len(findings), 1)
        # A warning, not a refusal: for a locally infused drug a curator who
        # has read the pharmacology may know the sign when valence does not.
        self.assertFalse(findings[0].is_hard)
        self.assertTrue(is_admissible(r))

    def test_molecular_record_without_an_identified_target_warns(self):
        r = synthetic.molecular_record(target_id=None, target_raw=None)
        self.assertIn("W014_MOLECULAR_TARGET_UNIDENTIFIED", self.codes(r))
        self.assertTrue(is_admissible(r))

    def test_condition_blind_record_warns(self):
        blind = synthetic.good_record(
            sex=None, diet=None, age_group=None, challenge=None
        )
        self.assertIn("W015_CONDITION_UNSPECIFIED", self.codes(blind))

    def test_not_reported_everywhere_still_counts_as_condition_blind(self):
        checked = synthetic.good_record(
            sex="not_reported",
            diet="not_reported",
            age_group="not_reported",
            challenge="not_reported",
            age_value=None,
            age_unit=None,
        )
        self.assertIn("W015_CONDITION_UNSPECIFIED", self.codes(checked))

    def test_one_known_condition_clears_the_warning(self):
        r = synthetic.good_record(
            sex="male", diet=None, age_group=None, challenge=None
        )
        self.assertNotIn("W015_CONDITION_UNSPECIFIED", self.codes(r))

    def test_every_new_code_is_declared_in_the_catalogue(self):
        for code in (
            "S007_PERTURBATION_CLASS_MISMATCH",
            "S008_TARGET_ON_NON_MOLECULAR",
            "W013_PERTURBATION_POLARITY_RULE_MISMATCH",
            "W014_MOLECULAR_TARGET_UNIDENTIFIED",
            "W015_CONDITION_UNSPECIFIED",
        ):
            self.assertIn(code, RULES)


class V030ConditionEvidenceTest(unittest.TestCase):
    """The rules schema v0.3.0 added: A008..A011 and W017 (D-037).

    v0.2.0 recorded a condition's provenance in one record-level
    ``condition_source``/``condition_locator`` pair, and A008 asked for a
    locator only when that pair said ``methods``. The pair was null on all 362
    curated records and could not have worked: a record's sex comes from the
    paper's Methods while its acute state, when stated at all, comes from the
    passage the record was read from, and one scalar pair cannot carry two
    provenances. v0.3.0 replaced it with one
    :class:`~syn4rch.schema.record.ConditionEvidence` entry per filled
    dimension, each carrying its own source, locator and quote.
    """

    #: Dimensions a fully specified fixture record fills, in schema order.
    DIMENSIONS = ("sex", "diet", "age_group", "challenge")

    def codes(self, record):
        return [f.code for f in validate(record)]

    def entry(self, **changes):
        """One well-formed ``sex`` evidence entry, with ``changes`` applied.

        Well-formed means: it names a condition dimension, its source is a
        route a condition may legitimately be read by, it carries the sentence
        it was read from and the passage that sentence is in, and its value is
        the value the record's own ``sex`` field holds. Each test below breaks
        exactly one of those.
        """
        entry = {
            "dimension": "sex",
            "value": "male",
            "source": "methods",
            "locator": "Methods para 2 (synthetic)",
            "quote": "Synthetic quote: male mice were used throughout.",
            "rule": "synthetic",
        }
        entry.update(changes)
        return entry

    def assert_hard(self, code, record, note=None):
        findings = validate(record)
        self.assertIn(
            code, codes(hard_findings(findings)), note or [str(f) for f in findings]
        )
        self.assertFalse(is_admissible(record), note)

    def test_the_good_record_triggers_none_of_them(self):
        """The fixture evidences all four of its dimensions, so none fire."""
        new_codes = {
            "A008_CONDITION_EVIDENCE_WITHOUT_LOCATOR",
            "A009_CONDITION_EVIDENCE_WITHOUT_QUOTE",
            "A010_CONDITION_EVIDENCE_MALFORMED",
            "A011_CONDITION_EVIDENCE_VALUE_MISMATCH",
            "W017_CONDITION_WITHOUT_EVIDENCE",
        }
        self.assertEqual(set(self.codes(synthetic.good_record())) & new_codes, set())
        self.assertEqual(
            set(self.codes(synthetic.molecular_record())) & new_codes, set()
        )

    def test_every_condition_evidence_entry_needs_a_locator(self):
        """A008.

        The v0.2.0 rule this replaces asked for a locator only on the Methods
        route, on the reasoning that a record's own passage is already located
        by ``source_locator``. That reasoning does not survive per-dimension
        evidence: the sentence that states the diet is not the sentence the
        effect was read from even when both are in Results, so an entry that
        names no passage cannot be checked whatever its source.
        """
        code = "A008_CONDITION_EVIDENCE_WITHOUT_LOCATOR"
        for locator in (None, "", "   ", "\n\t"):
            record = synthetic.good_record(
                condition_evidence=[self.entry(locator=locator)]
            )
            self.assert_hard(code, record, "locator=%r" % (locator,))
        good = synthetic.good_record(condition_evidence=[self.entry()])
        self.assertNotIn(code, self.codes(good))
        self.assertEqual(hard_findings(validate(good)), [])

    def test_a_passage_sourced_entry_needs_a_locator_too(self):
        """A008, the half of it v0.2.0 exempted.

        Under v0.2.0 ``condition_source=passage`` with no locator was
        admissible. Under v0.3.0 it is not, and this test is the old
        ``test_passage_sourced_condition_needs_no_locator`` inverted on
        purpose rather than deleted: the exemption was the defect.
        """
        code = "A008_CONDITION_EVIDENCE_WITHOUT_LOCATOR"
        bad = synthetic.good_record(
            condition_evidence=[self.entry(source="passage", locator=None)]
        )
        self.assert_hard(code, bad)
        good = synthetic.good_record(
            condition_evidence=[
                self.entry(source="passage", locator="Results para 1 (synthetic)")
            ]
        )
        self.assertNotIn(code, self.codes(good))

    def test_a_condition_evidence_entry_needs_its_quote(self):
        """A009: an entry exists because a sentence was found.

        Withholding the sentence is worse than leaving the dimension null:
        null says nothing is known, while a quoteless entry asserts that
        evidence exists and declines to show it.
        """
        code = "A009_CONDITION_EVIDENCE_WITHOUT_QUOTE"
        for quote in (None, "", "   ", "\n\t"):
            record = synthetic.good_record(
                condition_evidence=[self.entry(quote=quote)]
            )
            self.assert_hard(code, record, "quote=%r" % (quote,))
        self.assertNotIn(code, self.codes(synthetic.good_record()))

    def test_a_malformed_condition_evidence_entry_is_refused(self):
        """A010: an entry must name one of the four dimensions, exactly once.

        ``not_reported`` is not a source. It asserts that the paper was
        checked and says nothing, and an absence has no sentence to quote, so
        an entry claiming it is provenance for nothing.
        """
        code = "A010_CONDITION_EVIDENCE_MALFORMED"
        for changes in (
            {"dimension": "body_weight"},   # not a condition dimension
            {"dimension": "strain"},        # a real field, but not a dimension
            {"dimension": ""},              # no dimension at all
            {"source": None},               # no route
            {"source": "not_reported"},     # not a route
        ):
            record = synthetic.good_record(
                condition_evidence=[self.entry(**changes)]
            )
            self.assert_hard(code, record, repr(changes))
        twice = synthetic.good_record(
            condition_evidence=[
                self.entry(),
                self.entry(locator="Methods para 3 (synthetic)"),
            ]
        )
        self.assert_hard(code, twice, "one dimension, two entries")

    def test_condition_evidence_must_agree_with_the_field_it_describes(self):
        """A011: provenance pointing at a value that is not there.

        Three shapes of the same mistake: the entry contradicts its own field,
        the field is null, or the field says ``not_reported``. The last two
        matter most in practice -- they are what a half-finished curation pass
        leaves behind when a value is edited and its evidence is not.
        """
        code = "A011_CONDITION_EVIDENCE_VALUE_MISMATCH"
        self.assert_hard(
            code,
            synthetic.good_record(condition_evidence=[self.entry(value="female")]),
            "entry says female, field says male",
        )
        self.assert_hard(
            code,
            synthetic.good_record(sex=None, condition_evidence=[self.entry()]),
            "evidence for a dimension left null",
        )
        self.assert_hard(
            code,
            synthetic.good_record(
                sex="not_reported",
                condition_evidence=[self.entry(value="not_reported")],
            ),
            "evidence for an absence",
        )
        self.assertNotIn(code, self.codes(synthetic.good_record()))

    def test_an_unevidenced_condition_is_a_warning_not_a_refusal(self):
        """W017: the gap is counted, not hidden and not blocking.

        Every record curated before v0.3.0 is in this state. Making it hard
        would make the whole corpus inadmissible overnight, and making it
        silent would lose the only measure of how much condition provenance is
        still missing.
        """
        code = "W017_CONDITION_WITHOUT_EVIDENCE"
        record = synthetic.good_record(condition_evidence=[])
        findings = validate(record)
        self.assertEqual(hard_findings(findings), [], [str(f) for f in findings])
        self.assertTrue(is_admissible(record))
        self.assertIn(code, codes(warning_findings(findings)))
        message = [f.message for f in findings if f.code == code][0]
        for dimension in self.DIMENSIONS:
            self.assertIn(dimension, message, dimension)

    def test_a_partly_evidenced_record_is_warned_about_the_rest(self):
        code = "W017_CONDITION_WITHOUT_EVIDENCE"
        record = synthetic.good_record(condition_evidence=[self.entry()])
        message = [f.message for f in validate(record) if f.code == code][0]
        self.assertNotIn("sex", message)
        for dimension in ("diet", "age_group", "challenge"):
            self.assertIn(dimension, message, dimension)

    def test_an_absent_condition_is_not_warned_about(self):
        """W017 asks for evidence of a value, and ``not_reported`` is not one.

        A record that says the paper does not report its sex has nothing to
        quote. Demanding an entry there would push curators to invent one,
        which is the failure mode the rule exists to prevent.
        """
        record = synthetic.good_record(
            sex="not_reported",
            diet="not_reported",
            age_group="not_reported",
            challenge="not_reported",
            age_value=None,
            age_unit=None,
            condition_evidence=[],
        )
        self.assertNotIn("W017_CONDITION_WITHOUT_EVIDENCE", self.codes(record))
        self.assertEqual(hard_findings(validate(record)), [])

    def test_every_new_code_is_declared_in_the_catalogue(self):
        for code in (
            "A008_CONDITION_EVIDENCE_WITHOUT_LOCATOR",
            "A009_CONDITION_EVIDENCE_WITHOUT_QUOTE",
            "A010_CONDITION_EVIDENCE_MALFORMED",
            "A011_CONDITION_EVIDENCE_VALUE_MISMATCH",
            "W017_CONDITION_WITHOUT_EVIDENCE",
        ):
            self.assertIn(code, RULES)

    def test_the_four_hard_rules_are_distinct_findings(self):
        """One record, four defects, four codes: a report must say which.

        A missing locator, a missing quote, a malformed entry and an entry
        that disagrees with its field have four different fixes, and a single
        "something is wrong with the condition" finding would tell a curator
        none of them.
        """
        record = synthetic.good_record(
            condition_evidence=[
                self.entry(dimension="sex", locator=None),
                self.entry(dimension="diet", value="chow", quote=None),
                self.entry(dimension="body_weight"),
                self.entry(dimension="challenge", value="fasted"),
            ]
        )
        self.assertEqual(
            sorted(set(self.codes(record)) & {
                "A008_CONDITION_EVIDENCE_WITHOUT_LOCATOR",
                "A009_CONDITION_EVIDENCE_WITHOUT_QUOTE",
                "A010_CONDITION_EVIDENCE_MALFORMED",
                "A011_CONDITION_EVIDENCE_VALUE_MISMATCH",
            }),
            [
                "A008_CONDITION_EVIDENCE_WITHOUT_LOCATOR",
                "A009_CONDITION_EVIDENCE_WITHOUT_QUOTE",
                "A010_CONDITION_EVIDENCE_MALFORMED",
                "A011_CONDITION_EVIDENCE_VALUE_MISMATCH",
            ],
        )


class W016AlleleInversionTest(unittest.TestCase):
    """D-032: an allele that reduces function, expressed as a gain of function."""

    def codes(self, record):
        return [f.code for f in validate(record)]

    def mol(self, **kw):
        return synthetic.molecular_record(
            method="conditional_overexpression",
            perturbation_polarity="increase",
            **kw
        )

    def test_a_plain_overexpression_is_not_flagged(self):
        self.assertNotIn(
            "W016_POLARITY_MAY_BE_INVERTED_BY_ALLELE", self.codes(self.mol())
        )

    def test_a_dominant_negative_allele_is_flagged(self):
        record = self.mol(
            targeting_strategy="synthetic dominant-negative Synthgene1 knock-in"
        )
        self.assertIn(
            "W016_POLARITY_MAY_BE_INVERTED_BY_ALLELE", self.codes(record)
        )

    def test_the_real_wordings_are_all_caught(self):
        for wording in (
            "dominant negative Synthgene1",
            "dominant-negative allele",
            "kinase-dead Synthgene1",
            "deacetylase dead mutant",
            "catalytically inactive Synthgene1",
            "Rosa26 Synthgene1-H355Y",
            "Synthgene1-K310R knock-in",
            "Synthgene1-D45A substitution",
            "an inactive mutant of Synthgene1",
        ):
            record = self.mol(targeting_strategy="synthetic %s" % wording)
            self.assertIn(
                "W016_POLARITY_MAY_BE_INVERTED_BY_ALLELE",
                self.codes(record),
                wording,
            )

    def test_the_marker_is_read_from_actuator_too(self):
        record = self.mol(
            targeting_strategy="synthetic knock-in",
            actuator="Synthgene1-H355Y",
        )
        self.assertIn(
            "W016_POLARITY_MAY_BE_INVERTED_BY_ALLELE", self.codes(record)
        )

    def test_it_is_a_warning_and_never_rewrites_the_polarity(self):
        record = self.mol(targeting_strategy="synthetic kinase-dead Synthgene1")
        findings = [
            f for f in validate(record)
            if f.code == "W016_POLARITY_MAY_BE_INVERTED_BY_ALLELE"
        ]
        self.assertEqual(len(findings), 1)
        self.assertFalse(findings[0].is_hard)
        self.assertTrue(is_admissible(record))
        # The stored value is reported, not corrected: a mechanical guess at a
        # sign is how an inversion gets into the data (D-021, D-032).
        self.assertEqual(record.perturbation_polarity.value, "increase")
        self.assertIn("increase", findings[0].message)

    def test_a_circuit_perturbation_is_not_checked(self):
        """hM3Dq is not a dominant-negative allele; the rule is molecular-only."""
        record = synthetic.good_record(
            targeting_strategy="synthetic dominant-negative construct"
        )
        self.assertNotIn(
            "W016_POLARITY_MAY_BE_INVERTED_BY_ALLELE", self.codes(record)
        )

    def test_an_ordinary_mutant_word_alone_does_not_fire(self):
        """Too broad a pattern would flag most of the corpus and mean nothing."""
        for innocuous in (
            "floxed Synthgene1 mutant mice",
            "Synthgene1 transgenic line",
            "Synthgene1-GFP fusion",
            "AgRP-Cre x Synthgene1 flox",
        ):
            self.assertNotIn(
                "W016_POLARITY_MAY_BE_INVERTED_BY_ALLELE",
                self.codes(self.mol(targeting_strategy=innocuous)),
                innocuous,
            )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
