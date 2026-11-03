"""Passages and the locator grammar: stable, unique, and resolvable."""

import _path  # noqa: F401
import os
import unittest

from syn4rch.curate.jats import Caption, ParsedPaper, Section, parse_jats_file
from syn4rch.curate.passages import (
    KIND_FIGURE_CAPTION,
    build_index,
    parse_locator,
    split_passages,
)

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "jats")


def load(name):
    return parse_jats_file(os.path.join(FIX, name))


class TestLocators(unittest.TestCase):
    def test_complete_fixture_locators(self):
        locators = [p.locator for p in split_passages(load("synthetic_complete.xml"))]
        self.assertEqual(
            locators,
            [
                "Title",
                "Abstract para 1",
                "Introduction para 1",
                "Introduction para 2",
                "Methods para 1",
                "Methods para 2",
                "Results para 1",
                "Results para 2",
                "Discussion para 1",
                "Fig 3 caption",
                "Table 1 caption",
            ],
        )

    def test_stable_across_parses(self):
        a = split_passages(load("synthetic_complete.xml"))
        b = split_passages(load("synthetic_complete.xml"))
        self.assertEqual(a, b)

    def test_unique_within_paper_even_with_collisions(self):
        paper = ParsedPaper(
            sections=[
                Section("Results", "Results", ["one", "two"]),
                Section("Results", "Results", ["three"]),
                Section("Other", None, ["four"]),
            ],
            figures=[
                Caption("figure", "Figure 2", "2", 1, None, "first fig 2"),
                Caption("figure", "Figure 2", "2", 2, None, "second fig 2"),
            ],
        )
        locators = [p.locator for p in split_passages(paper)]
        self.assertEqual(len(locators), len(set(locators)))
        self.assertIn("Results (2) para 1", locators)
        self.assertIn("Section 3 para 1", locators)
        self.assertIn("Fig 2 caption (2)", locators)

    def test_non_canonical_section_uses_its_heading_and_unlabelled_fig_its_order(self):
        passages = split_passages(load("synthetic_no_abstract.xml"))
        by_loc = {p.locator: p for p in passages}
        self.assertIn("Limitations of this invented study para 1", by_loc)
        self.assertTrue(by_loc["Fig 1 caption"].numbered_by_order)
        self.assertIn("Table II caption", by_loc)

    def test_empty_paper_has_no_passages(self):
        self.assertEqual(split_passages(ParsedPaper()), [])

    def test_every_written_locator_parses_and_resolves_to_itself(self):
        for name in ("synthetic_complete.xml", "synthetic_no_abstract.xml",
                     "synthetic_namespaced.xml", "synthetic_truncated.xml"):
            index = build_index(load(name))
            for passage in index.passages:
                self.assertIsNotNone(parse_locator(passage.locator), passage.locator)
                hits, how = index.resolve(passage.locator)
                self.assertIn(passage, hits, (name, passage.locator))
                self.assertEqual(how, "exact")


class TestParseLocator(unittest.TestCase):
    def test_accepted_forms(self):
        cases = {
            "Results para 4": ("section_para", "Results", 4),
            "Methods, paragraph 2": ("section_para", "Methods", 2),
            "Discussion par 1": ("section_para", "Discussion", 1),
            "Title": ("title", "Title", None),
        }
        for text, (kind, name, number) in cases.items():
            ref = parse_locator(text)
            self.assertEqual((ref.kind, ref.name, ref.number), (kind, name, number), text)

    def test_figure_and_table_forms(self):
        ref = parse_locator("Fig 3C caption")
        self.assertEqual((ref.kind, ref.float_number, ref.panel), ("figure", "3", "c"))
        ref = parse_locator("Figure 3")
        self.assertEqual((ref.kind, ref.float_number, ref.panel), ("figure", "3", None))
        ref = parse_locator("Fig. 2B-D legend")
        self.assertEqual((ref.kind, ref.float_number), ("figure", "2"))
        ref = parse_locator("Table II caption")
        self.assertEqual((ref.kind, ref.float_number), ("table", "II"))
        ref = parse_locator("Fig 2 caption (2)")
        self.assertEqual((ref.kind, ref.float_number), ("figure", "2"))

    def test_unparseable(self):
        for text in (None, "", "   ", "page 7 ???", "#4"):
            self.assertIsNone(parse_locator(text), text)


class TestResolve(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index = build_index(load("synthetic_complete.xml"))

    def test_exact(self):
        hits, how = self.index.resolve("Results para 2")
        self.assertEqual(how, "exact")
        self.assertEqual([p.locator for p in hits], ["Results para 2"])

    def test_case_and_spacing_of_locator_are_not_load_bearing(self):
        hits, how = self.index.resolve("  results   PARA 2 ")
        self.assertEqual([p.locator for p in hits], ["Results para 2"])

    def test_panel_suffix(self):
        hits, how = self.index.resolve("Fig 3C caption")
        self.assertEqual(how, "panel_suffix")
        self.assertEqual(hits[0].kind, KIND_FIGURE_CAPTION)

    def test_section_synonym_requires_exact_paragraph_number(self):
        hits, how = self.index.resolve("Materials and Methods para 2")
        self.assertEqual(how, "section_name")
        self.assertEqual(hits[0].locator, "Methods para 2")
        self.assertEqual(self.index.resolve("Materials and Methods para 9")[1], "not_found")

    def test_whole_section(self):
        hits, how = self.index.resolve("Results")
        self.assertEqual(how, "section_whole")
        self.assertEqual(len(hits), 2)

    def test_not_found(self):
        for loc in ("Results para 99", "Fig 9 caption", "Table 4",
                    "Supplementary Figure 2", "Acknowledgements para 1"):
            hits, how = self.index.resolve(loc)
            self.assertEqual(hits, [], loc)
            self.assertEqual(how, "not_found", loc)

    def test_unparseable(self):
        self.assertEqual(self.index.resolve("#4"), ([], "unparseable"))


if __name__ == "__main__":
    unittest.main()
