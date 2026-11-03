"""Extraction unit candidates: localized, scored by a checkable rule, never records."""

import _path  # noqa: F401
import dataclasses
import os
import unittest

from syn4rch.curate.jats import ParsedPaper, Section, parse_jats_file
from syn4rch.curate.passages import split_passages
from syn4rch.curate.units import (
    SCORE_WEIGHTS,
    ExtractionUnitCandidate,
    match_passage,
    paper_term_summary,
    units_from_paper,
)
from syn4rch.harvest import screen
from syn4rch.schema import TierARecord

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "jats")


def complete():
    return parse_jats_file(os.path.join(FIX, "synthetic_complete.xml"))


def paper_of(*paragraphs):
    return ParsedPaper(
        doi="10.0000/synthetic.unit",
        sections=[Section("Results", "Results", list(paragraphs))],
    )


class TestLocalisation(unittest.TestCase):
    def test_all_three_in_one_passage_yields_units(self):
        units = units_from_paper(
            paper_of("Chemogenetic activation of AgRP neurons increased food intake.")
        )
        self.assertTrue(units)
        self.assertEqual({u.candidate_population_id for u in units}, {"arc_agrp"})
        self.assertEqual({u.candidate_outcome_id for u in units}, {"food_intake"})
        for unit in units:
            self.assertEqual(unit.locator, "Results para 1")

    def test_terms_split_across_passages_yield_nothing(self):
        paper = paper_of(
            "We studied AgRP neurons in this invented work.",
            "Chemogenetic activation increased food intake.",
        )
        self.assertEqual(units_from_paper(paper), [])
        summary = paper_term_summary(paper)
        self.assertEqual(summary["populations"], ["arc_agrp"])
        self.assertIn("food_intake", summary["outcomes"])

    def test_passage_text_is_verbatim(self):
        paper = complete()
        texts = {p.locator: p.text for p in split_passages(paper)}
        for unit in units_from_paper(paper):
            self.assertEqual(unit.passage_text, texts[unit.locator])

    def test_reuses_screen_matcher_and_its_ambiguity_stoplist(self):
        self.assertIn("ARC", screen.AMBIGUOUS_ACRONYMS)
        passage = split_passages(paper_of("Invented sentence naming only the ARC."))[0]
        self.assertEqual(match_passage(passage).regions, ())
        passage = split_passages(paper_of("Invented sentence on the arcuate nucleus."))[0]
        self.assertEqual(match_passage(passage).regions, ("arc",))


class TestScoring(unittest.TestCase):
    def test_score_is_sum_of_documented_weights(self):
        for unit in units_from_paper(complete()):
            for key, value in unit.score_breakdown.items():
                self.assertEqual(SCORE_WEIGHTS[key], value)
            self.assertEqual(unit.score, sum(unit.score_breakdown.values()))

    def test_sorted_descending_and_deterministic(self):
        a = units_from_paper(complete())
        b = units_from_paper(complete())
        self.assertEqual(a, b)
        scores = [u.score for u in a]
        self.assertEqual(scores, sorted(scores, reverse=True))

    def test_results_outranks_introduction(self):
        units = units_from_paper(complete())
        best = {}
        for unit in units:
            best.setdefault(unit.section_label, unit.score)
        self.assertGreater(best["Results"], best["Introduction"])

    def test_cap_per_passage(self):
        paper = complete()
        full = units_from_paper(paper)
        capped = units_from_paper(paper, max_per_passage=1)
        self.assertEqual(len(capped), len({u.locator for u in full}))

    def test_null_result_wording_is_not_ranked_below_positive(self):
        up = units_from_paper(paper_of("Chemogenetic activation of AgRP neurons increased food intake."))
        flat = units_from_paper(paper_of("Chemogenetic activation of AgRP neurons left food intake unchanged."))
        self.assertEqual(max(u.score for u in up), max(u.score for u in flat))


class TestNotARecord(unittest.TestCase):
    def test_type_carries_no_claim_fields(self):
        names = {f.name for f in dataclasses.fields(ExtractionUnitCandidate)}
        for forbidden in ("direction", "valence", "effect_size_value", "confidence",
                          "support_quote", "source_locator", "record_id",
                          "population_id", "outcome_id"):
            self.assertNotIn(forbidden, names)
        self.assertFalse(issubclass(ExtractionUnitCandidate, TierARecord))

    def test_serialised_row_is_stamped(self):
        row = units_from_paper(complete())[0].to_dict()
        self.assertEqual(row["_kind"], "extraction_unit_candidate")
        self.assertIs(row["_not_a_record"], True)
        with self.assertRaises(Exception):
            TierARecord.from_dict(row)

    def test_unit_ids_unique(self):
        units = units_from_paper(complete())
        self.assertEqual(len({u.unit_id for u in units}), len(units))


if __name__ == "__main__":
    unittest.main()
