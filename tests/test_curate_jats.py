"""JATS parsing: structure, inline markup, and never raising on bad input.

All fixtures under tests/fixtures/jats/ are hand-written synthetic XML.
"""

import _path  # noqa: F401
import os
import unittest

from syn4rch.curate.jats import classify_section, parse_jats, parse_jats_file

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "jats")


def fixture(name):
    return os.path.join(FIX, name)


def read(name):
    with open(fixture(name), encoding="utf-8") as handle:
        return handle.read()


class TestCompletePaper(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.paper = parse_jats_file(fixture("synthetic_complete.xml"))

    def test_clean_parse(self):
        self.assertEqual(self.paper.parse_warnings, [])

    def test_identifiers(self):
        self.assertEqual(self.paper.pmcid, "PMC0000001")
        self.assertEqual(self.paper.pmid, "10000001")
        self.assertEqual(self.paper.doi, "10.0000/synthetic.0001")

    def test_bibliographic(self):
        self.assertEqual(
            self.paper.title,
            "Synthetic fixture: chemogenetic activation of Agrp neurons and food intake",
        )
        self.assertEqual(self.paper.journal, "Journal of Synthetic Test Data")
        # pub-date wins over the <history> accepted date (2000).
        self.assertEqual(self.paper.year, 2001)

    def test_authors_and_last_author_affiliation_via_xref(self):
        self.assertEqual(self.paper.first_author, "Firstauthor A. B.")
        self.assertEqual(self.paper.last_author, "Lastauthor D. E.")
        # Resolved through rid="aff2", label "2" dropped.
        self.assertEqual(
            self.paper.last_author_affiliation,
            "Department of Fixtures, University of Nowhere",
        )

    def test_section_order_and_labels(self):
        self.assertEqual(
            [s.label for s in self.paper.sections],
            ["Abstract", "Introduction", "Methods", "Results", "Discussion"],
        )
        self.assertEqual(self.paper.section("Methods").title, "Materials and Methods")

    def test_nested_subsection_paragraphs_are_flattened_in_order(self):
        results = self.paper.section("Results")
        self.assertEqual(len(results.paragraphs), 2)
        self.assertEqual(results.subsection_titles, ["Blood glucose"])
        self.assertEqual(len(self.paper.section("Methods").paragraphs), 2)

    def test_nested_inline_markup(self):
        first = self.paper.section("Results").paragraphs[0]
        self.assertEqual(
            first,
            "Chemogenetic activation of Agrp neurons increased food intake "
            "(n = 8 mice, p < 0.01), and oxygen consumption2 was unchanged24h.",
        )

    def test_xref_text_is_kept_inline(self):
        second = self.paper.section("Results").paragraphs[1]
        self.assertTrue(second.startswith("Blood glucose decreased12 after"))
        self.assertNotIn("\n", second)
        self.assertNotIn("  ", second)

    def test_figure_caption_not_in_running_text(self):
        for para in self.paper.section("Results").paragraphs:
            self.assertNotIn("Invented panel description", para)

    def test_figure_and_table_captions(self):
        self.assertEqual(len(self.paper.figures), 1)
        fig = self.paper.figures[0]
        self.assertEqual(fig.label, "Figure 3")
        self.assertEqual(fig.number, "3")
        self.assertTrue(fig.text.startswith("Figure 3 Food intake after"))
        self.assertEqual(len(self.paper.tables), 1)
        self.assertEqual(self.paper.tables[0].label, "Table 1")
        self.assertEqual(self.paper.tables[0].number, "1")

    def test_back_matter_and_graphical_abstract_excluded(self):
        everything = " ".join(
            p for s in self.paper.sections for p in s.paragraphs
        )
        self.assertNotIn("acknowledgement", everything)
        self.assertNotIn("graphical abstract", everything)
        self.assertNotIn("Invented reference", everything)


class TestNoAbstract(unittest.TestCase):
    def test_no_abstract(self):
        paper = parse_jats_file(fixture("synthetic_no_abstract.xml"))
        self.assertIsNone(paper.section("Abstract"))
        self.assertIn("no abstract found", paper.parse_warnings)
        self.assertEqual(paper.section("Results").paragraphs[0][:20], "Optogenetic stimulat")
        # Single affiliation, no xref: it is the only thing it can be.
        self.assertEqual(paper.last_author_affiliation, "Institute of Synthetic Data, Nowhere")

    def test_unlabelled_figure_and_roman_table(self):
        paper = parse_jats_file(fixture("synthetic_no_abstract.xml"))
        self.assertIsNone(paper.figures[0].number)
        self.assertEqual(paper.figures[0].ordinal, 1)
        self.assertEqual(paper.tables[0].number, "II")

    def test_non_canonical_section_is_other(self):
        paper = parse_jats_file(fixture("synthetic_no_abstract.xml"))
        other = [s for s in paper.sections if s.label == "Other"]
        self.assertEqual(len(other), 1)
        self.assertEqual(other[0].title, "Limitations of this invented study")


class TestBadInputNeverRaises(unittest.TestCase):
    def test_truncated_file_recovers_prefix(self):
        paper = parse_jats_file(fixture("synthetic_truncated.xml"))
        self.assertTrue(paper.parse_warnings)
        self.assertEqual(paper.pmcid, "PMC0000003")
        results = paper.section("Results")
        self.assertIsNotNone(results)
        # The text right before the cut is recovered, unaltered.
        self.assertEqual(
            results.paragraphs,
            [
                "An invented first results paragraph that was written to disk in full.",
                "An invented second results paragraph, cut off right he",
            ],
        )

    def test_malformed_file_recovers_text_before_mismatch(self):
        paper = parse_jats_file(fixture("synthetic_malformed.xml"))
        self.assertTrue(any("not well-formed" in w for w in paper.parse_warnings))
        self.assertEqual(paper.doi, "10.0000/synthetic.0004")
        self.assertEqual(
            paper.section("Results").paragraphs,
            ["An invented results paragraph closed with the wrong tag."],
        )

    def test_every_prefix_of_a_document_parses_without_raising(self):
        text = read("synthetic_complete.xml")
        for cut in range(0, len(text), 37):
            paper = parse_jats(text[:cut])
            if cut < len(text) - 20:
                self.assertTrue(paper.parse_warnings, cut)

    def test_prefix_text_is_always_a_substring_of_the_full_parse(self):
        # Recovery repairs markup, never text.
        full = parse_jats(read("synthetic_complete.xml"))
        full_paras = [p for s in full.sections for p in s.paragraphs]
        text = read("synthetic_complete.xml")
        for cut in range(0, len(text), 53):
            partial = parse_jats(text[:cut])
            for section in partial.sections:
                for para in section.paragraphs:
                    self.assertTrue(
                        any(para.rstrip() in fp for fp in full_paras),
                        (cut, para),
                    )

    def test_degenerate_inputs(self):
        for bad in ("", "   \n", None, "<not xml", "plain text", "<a><b></a>",
                    "<?xml version='1.0'?>", "<!-- only a comment", "\x00\x01"):
            paper = parse_jats(bad)
            self.assertTrue(paper.parse_warnings, repr(bad))
            self.assertTrue(paper.is_empty)

    def test_missing_file(self):
        paper = parse_jats_file(fixture("does_not_exist.xml"))
        self.assertTrue(paper.parse_warnings)
        self.assertTrue(paper.is_empty)


class TestNamespaces(unittest.TestCase):
    def test_namespaced_document_parses(self):
        paper = parse_jats_file(fixture("synthetic_namespaced.xml"))
        self.assertEqual(paper.parse_warnings, [])
        self.assertEqual(paper.pmcid, "PMC0000005")
        self.assertEqual(paper.year, 2005)
        self.assertEqual([s.label for s in paper.sections], ["Abstract", "Results"])
        self.assertEqual(paper.figures[0].label, "Fig. 1")
        self.assertEqual(paper.figures[0].number, "1")


class TestClassifySection(unittest.TestCase):
    def test_mapping(self):
        cases = {
            ("intro", None): "Introduction",
            (None, "Background"): "Introduction",
            ("materials|methods", None): "Methods",
            (None, "Experimental Procedures"): "Methods",
            (None, "3. Results"): "Results",
            (None, "Results and Discussion"): "Results",
            (None, "Conclusions"): "Discussion",
            (None, "Methodological notes"): "Other",
            (None, "Funding"): "Other",
            (None, None): "Other",
        }
        for (sec_type, title), expected in cases.items():
            self.assertEqual(classify_section(sec_type, title), expected, (sec_type, title))


if __name__ == "__main__":
    unittest.main()
