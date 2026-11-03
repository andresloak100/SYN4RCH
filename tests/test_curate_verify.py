"""The verbatim gate. Every rejection here is the gate working.

Records are built from tests/synthetic.py (fabricated) and checked against the
hand-written synthetic JATS fixture.
"""

import _path  # noqa: F401
import copy
import os
import unittest

import synthetic
from syn4rch.curate import verify
from syn4rch.curate.jats import parse_jats_file
from syn4rch.curate.verify import (
    gate,
    normalize_for_match,
    quote_occurs_in,
    verify_support,
)

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "jats")
DOI = "10.0000/synthetic.0001"
#: Verbatim from "Results para 1" of the fixture (as flattened).
QUOTE = (
    "Chemogenetic activation of Agrp neurons increased food intake "
    "(n = 8 mice, p < 0.01)"
)


def paper():
    return parse_jats_file(os.path.join(FIX, "synthetic_complete.xml"))


def record(**overrides):
    data = dict(
        study_id="doi:" + DOI,
        doi=DOI,
        support_quote=QUOTE,
        source_locator="Results para 1",
    )
    data.update(overrides)
    return synthetic.good_record(**data)


def codes(result):
    return [f.code for f in result.blocking]


class TestNormalisation(unittest.TestCase):
    def test_whitespace_dashes_quotes_nbsp(self):
        self.assertEqual(
            normalize_for_match("  a  b\n\t c–d “e” ’f "),
            "a b c-d \"e\" 'f",
        )

    def test_case_and_words_preserved(self):
        self.assertEqual(normalize_for_match("Food Intake"), "Food Intake")
        self.assertFalse(quote_occurs_in("food intake", "Food intake rose"))

    def test_empty_quote_never_verifies(self):
        self.assertFalse(quote_occurs_in("", "anything"))
        self.assertFalse(quote_occurs_in("   ", "anything"))

    def test_space_set_covers_unicode_zs(self):
        self.assertTrue(verify.zs_coverage_is_complete())

    def test_no_approximate_matching_imported(self):
        self.assertEqual(verify.forbidden_imports(), [])


class TestVerifySupport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.paper = paper()

    def test_exact_quote_accepted(self):
        result = verify_support(record(), self.paper)
        self.assertTrue(result.ok, result.findings)
        self.assertEqual(result.quote_found_in, ("Results para 1",))
        self.assertEqual(result.located_passages, ("Results para 1",))

    def test_whitespace_only_difference_accepted(self):
        quote = "Chemogenetic   activation of\n Agrp neurons  increased food\tintake"
        self.assertTrue(verify_support(record(support_quote=quote), self.paper).ok)

    def test_unicode_dash_and_quote_folding_accepted_both_ways(self):
        from syn4rch.curate.jats import parse_jats

        src = parse_jats(
            "<article><front><article-meta><article-id pub-id-type='doi'>"
            + DOI + "</article-id></article-meta></front><body><sec><title>Results"
            "</title><p>Intake rose 2\u20133\u00a0fold (\u201cinvented\u201d) in"
            " Agrp\u2011Cre mice\u2019 cages.</p></sec></body></article>"
        )
        ascii_quote = 'Intake rose 2-3 fold ("invented") in Agrp-Cre mice\' cages.'
        self.assertTrue(
            verify_support(record(support_quote=ascii_quote), src).ok
        )
        # Folding does not license a different character: an underscore is not a dash.
        self.assertFalse(
            verify_support(record(support_quote="Intake rose 2_3 fold"), src).ok
        )

    def test_inserted_word_rejected(self):
        quote = QUOTE.replace("activation of", "activation of all")
        result = verify_support(record(support_quote=quote), self.paper)
        self.assertFalse(result.ok)
        self.assertIn("C004_QUOTE_NOT_VERBATIM", codes(result))

    def test_changed_word_rejected(self):
        quote = QUOTE.replace("increased", "decreased")
        result = verify_support(record(support_quote=quote), self.paper)
        self.assertIn("C004_QUOTE_NOT_VERBATIM", codes(result))

    def test_dropped_word_rejected(self):
        quote = QUOTE.replace("food ", "")
        self.assertIn("C004_QUOTE_NOT_VERBATIM",
                      codes(verify_support(record(support_quote=quote), self.paper)))

    def test_different_case_rejected(self):
        for quote in (QUOTE.lower(), QUOTE.upper(),
                      QUOTE.replace("Chemogenetic", "chemogenetic")):
            result = verify_support(record(support_quote=quote), self.paper)
            self.assertIn("C004_QUOTE_NOT_VERBATIM", codes(result), quote)

    def test_quote_stitched_across_paragraphs_rejected(self):
        quote = "was unchanged24h. Blood glucose decreased12"
        self.assertIn("C004_QUOTE_NOT_VERBATIM",
                      codes(verify_support(record(support_quote=quote), self.paper)))

    def test_nonexistent_locator_rejected(self):
        for loc in ("Results para 9", "Fig 7 caption", "Table 3 caption",
                    "Supplementary Figure 2"):
            result = verify_support(record(source_locator=loc), self.paper)
            self.assertFalse(result.ok, loc)
            self.assertIn("C006_LOCATOR_NOT_IN_PAPER", codes(result), loc)

    def test_unparseable_locator_rejected(self):
        result = verify_support(record(source_locator="#4"), self.paper)
        self.assertIn("C005_LOCATOR_UNPARSEABLE", codes(result))

    def test_missing_quote_and_locator_rejected(self):
        result = verify_support(record(support_quote=None, source_locator="  "), self.paper)
        self.assertIn("C002_MISSING_SUPPORT_QUOTE", codes(result))
        self.assertIn("C003_MISSING_SOURCE_LOCATOR", codes(result))

    def test_no_paper_is_rule_2_5(self):
        result = verify_support(record(), None)
        self.assertEqual(codes(result), ["C001_TEXT_NOT_RETRIEVED"])

    def test_empty_paper_is_rule_2_5(self):
        from syn4rch.curate.jats import parse_jats

        result = verify_support(record(), parse_jats("<article"))
        self.assertIn("C001_TEXT_NOT_RETRIEVED", codes(result))

    def test_panel_suffix_resolves_with_a_note(self):
        result = verify_support(
            record(support_quote="(C) Food intake increased after chemogenetic activation",
                   source_locator="Fig 3C caption"),
            self.paper,
        )
        self.assertTrue(result.ok, result.findings)
        self.assertIn("W103_LOCATOR_PANEL_SUFFIX", [f.code for f in result.notes])

    def test_quote_in_other_passage_is_rejected_by_default(self):
        # D-018: a quote must sit in the passage its locator names.
        rec = record(source_locator="Results para 2")
        strict = verify_support(rec, self.paper)
        self.assertFalse(strict.ok)
        self.assertIn("W101_QUOTE_NOT_IN_LOCATED_PASSAGE", codes(strict))
        lenient = verify_support(rec, self.paper, strict_locator_passage=False)
        self.assertTrue(lenient.ok)
        self.assertIn("W101_QUOTE_NOT_IN_LOCATED_PASSAGE", [f.code for f in lenient.notes])

    def test_record_never_modified(self):
        rec = record(support_quote=QUOTE.lower())
        before = copy.deepcopy(rec.to_dict())
        verify_support(rec, self.paper)
        gate([rec], [self.paper])
        self.assertEqual(rec.to_dict(), before)


class TestGate(unittest.TestCase):
    def test_admits_and_rejects_with_machine_readable_reasons(self):
        good = record()
        bad_quote = record(support_quote=QUOTE.replace("increased", "raised"))
        bad_loc = record(source_locator="Results para 40")
        other_paper = record(study_id="doi:10.0000/synthetic.9999", doi="10.0000/synthetic.9999")
        admitted, rejected = gate([good, bad_quote, bad_loc, other_paper], [paper()])
        self.assertEqual(admitted, [good])
        self.assertIs(admitted[0], good)
        by_id = {r.record_id: r for r in rejected}
        self.assertEqual(by_id[bad_quote.record_id].code, "C004_QUOTE_NOT_VERBATIM")
        self.assertEqual(by_id[bad_loc.record_id].code, "C006_LOCATOR_NOT_IN_PAPER")
        self.assertEqual(by_id[other_paper.record_id].code, "C001_TEXT_NOT_RETRIEVED")
        for rejection in rejected:
            self.assertIn(rejection.code, verify.CODES)
            self.assertTrue(rejection.to_dict()["reasons"])

    def test_rejected_record_returned_unrepaired(self):
        bad = record(support_quote=QUOTE.lower())
        _admitted, rejected = gate([bad], [paper()])
        self.assertIs(rejected[0].record, bad)
        self.assertEqual(rejected[0].record.support_quote, QUOTE.lower())

    def test_schema_hard_finding_rejects_even_with_verbatim_quote(self):
        rec = record(evidence_class="anatomical")
        _admitted, rejected = gate([rec], [paper()])
        self.assertEqual(rejected[0].code, "C007_SCHEMA_HARD")
        self.assertTrue(rejected[0].schema_findings)

    def test_no_papers_rejects_everything(self):
        admitted, rejected = gate([record(), record(source_locator="Title")], None)
        self.assertEqual(admitted, [])
        self.assertEqual({r.code for r in rejected}, {"C001_TEXT_NOT_RETRIEVED"})

    def test_lookup_by_pmcid_or_pmid(self):
        rec = record(study_id="pmid:10000001", doi=None, pmid="10000001")
        admitted, _ = gate([rec], [paper()])
        self.assertEqual(admitted, [rec])
        rec = record(study_id="pmcid:PMC0000001", doi=None)
        admitted, _ = gate([rec], {"x": paper()})
        self.assertEqual(admitted, [rec])


if __name__ == "__main__":
    unittest.main()
