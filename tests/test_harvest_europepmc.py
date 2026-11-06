"""The Europe PMC parser, and the two fields it threw away until D-061.

`search()` has always requested `resultType=core`, the result type that carries
`meshHeadingList` and author affiliations. The parser hardcoded `mesh_terms: []`
and `last_author_affiliation: None`, the second under a comment asserting the
API does not return them. Meanwhile `screen.py` joins mesh into the matched text
and stamps every candidate `screened_on: "title+abstract+mesh"`.

So these tests pin two things: that both fields are now parsed, and that the
shapes Europe PMC actually uses are handled -- the affiliation lives under
`authorList.author[].authorAffiliationDetailsList`, which is not a name a reader
would guess, and a single-element list arrives as a bare object rather than a
list of one.
"""

import unittest

from _path import *  # noqa: F401,F403  (puts src/ on sys.path)

from syn4rch.harvest.europepmc import (
    _last_author_affiliation,
    _mesh_terms,
    _norm_record,
)


def result(**overrides):
    """A `resultType=core` result, shaped as Europe PMC really returns one."""
    data = {
        "pmid": "12345678",
        "pmcid": "pmc7654321",
        "doi": "10.1000/SYNTHETIC.1",
        "title": "Synthetic record",
        "abstractText": "Synthetic abstract.",
        "pubYear": "2024",
        "authorString": "First A, Middle B, Last C.",
        "isOpenAccess": "Y",
        "inEPMC": "Y",
        "meshHeadingList": {
            "meshHeading": [
                {"majorTopic_YN": "N", "descriptorName": "Animals"},
                {
                    "majorTopic_YN": "Y",
                    "descriptorName": "Hypothalamus",
                    "meshQualifierList": {
                        "meshQualifier": [
                            {"abbreviation": "ME", "qualifierName": "metabolism"}
                        ]
                    },
                },
                {"majorTopic_YN": "Y", "descriptorName": "Optogenetics"},
            ]
        },
        "authorList": {
            "author": [
                {
                    "lastName": "First",
                    "authorAffiliationDetailsList": {
                        "authorAffiliation": [{"affiliation": "Institute of First"}]
                    },
                },
                {
                    "lastName": "Last",
                    "authorAffiliationDetailsList": {
                        "authorAffiliation": [{"affiliation": "Institute of Last"}]
                    },
                },
            ]
        },
        "affiliation": "Institute of First",
    }
    data.update(overrides)
    return data


class MeshTermsTest(unittest.TestCase):
    def test_descriptors_are_parsed_in_order(self):
        self.assertEqual(
            _mesh_terms(result()), ["Animals", "Hypothalamus", "Optogenetics"]
        )

    def test_qualifiers_are_dropped(self):
        """The screener matches descriptor names.

        "Hypothalamus/metabolism" would match neither "Hypothalamus" nor
        "metabolism", so a qualified descriptor must still arrive as its plain
        descriptor.
        """
        terms = _mesh_terms(result())
        self.assertIn("Hypothalamus", terms)
        self.assertNotIn("metabolism", " ".join(terms).split())
        for term in terms:
            self.assertNotIn("/", term)

    def test_a_single_heading_arrives_as_an_object_not_a_list(self):
        """Europe PMC collapses one-element lists; the parser must not care."""
        item = result(
            meshHeadingList={"meshHeading": {"descriptorName": "Mice"}}
        )
        self.assertEqual(_mesh_terms(item), ["Mice"])

    def test_duplicates_are_collapsed(self):
        item = result(
            meshHeadingList={
                "meshHeading": [
                    {"descriptorName": "Mice"},
                    {"descriptorName": "Mice"},
                ]
            }
        )
        self.assertEqual(_mesh_terms(item), ["Mice"])

    def test_absent_or_malformed_mesh_yields_an_empty_list(self):
        for item in (
            result(meshHeadingList=None),
            result(meshHeadingList={}),
            result(meshHeadingList={"meshHeading": "not a list"}),
            result(meshHeadingList=[]),
        ):
            self.assertEqual(_mesh_terms(item), [])

    def test_a_heading_with_no_descriptor_is_skipped_not_blank(self):
        item = result(
            meshHeadingList={
                "meshHeading": [{"majorTopic_YN": "N"}, {"descriptorName": "Mice"}]
            }
        )
        self.assertEqual(_mesh_terms(item), ["Mice"])


class LastAuthorAffiliationTest(unittest.TestCase):
    def test_the_last_author_is_preferred_over_the_article_level_field(self):
        """lab_key is about the last author, and `affiliation` is usually the first's."""
        self.assertEqual(
            _last_author_affiliation(result()), "Institute of Last"
        )

    def test_the_article_level_field_is_the_fallback(self):
        item = result(authorList=None)
        self.assertEqual(_last_author_affiliation(item), "Institute of First")

    def test_several_affiliations_on_one_author_are_joined_without_duplicates(self):
        item = result(
            authorList={
                "author": [
                    {
                        "lastName": "Last",
                        "authorAffiliationDetailsList": {
                            "authorAffiliation": [
                                {"affiliation": "Place One"},
                                {"affiliation": "Place Two"},
                                {"affiliation": "Place One"},
                            ]
                        },
                    }
                ]
            }
        )
        self.assertEqual(
            _last_author_affiliation(item), "Place One; Place Two"
        )

    def test_a_single_affiliation_arrives_as_an_object(self):
        item = result(
            authorList={
                "author": {
                    "lastName": "Only",
                    "authorAffiliationDetailsList": {
                        "authorAffiliation": {"affiliation": "Sole Place"}
                    },
                }
            }
        )
        self.assertEqual(_last_author_affiliation(item), "Sole Place")

    def test_none_when_nothing_states_an_affiliation(self):
        item = result(authorList=None, affiliation=None)
        self.assertIsNone(_last_author_affiliation(item))

    def test_a_blank_affiliation_is_none_not_empty_string(self):
        """W010/W018 test for a blank, so "" and None must not diverge here."""
        item = result(authorList=None, affiliation="   ")
        self.assertIsNone(_last_author_affiliation(item))


class NormRecordTest(unittest.TestCase):
    def test_both_fields_reach_the_candidate_record(self):
        """The regression that matters: these were [] and None for ~3,750 candidates."""
        record = _norm_record(result())
        self.assertEqual(
            record["mesh_terms"], ["Animals", "Hypothalamus", "Optogenetics"]
        )
        self.assertEqual(record["last_author_affiliation"], "Institute of Last")

    def test_a_record_with_neither_still_parses(self):
        record = _norm_record(result(meshHeadingList=None, authorList=None, affiliation=None))
        self.assertEqual(record["mesh_terms"], [])
        self.assertIsNone(record["last_author_affiliation"])
        self.assertEqual(record["pmcid"], "PMC7654321")

    def test_the_rest_of_the_record_is_unchanged(self):
        record = _norm_record(result())
        self.assertEqual(record["pmid"], "12345678")
        self.assertEqual(record["doi"], "10.1000/synthetic.1")
        self.assertEqual(record["year"], 2024)
        self.assertEqual(record["last_author"], "Last C.")
        self.assertEqual(record["retrieved_from"], "europepmc_oa")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
