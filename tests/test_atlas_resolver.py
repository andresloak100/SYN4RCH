"""Allen ontology resolver: exact matches only, never forced."""

import _path  # noqa: F401
import unittest

from syn4rch.atlas import ontology


def node(sid, acr, name, children=()):
    return {"id": sid, "acronym": acr, "name": name, "children": list(children)}


GRAPH = {"msg": [node(997, "root", "root", [
    node(1097, "HY", "Hypothalamus", [
        node(223, "ARH", "Arcuate hypothalamic nucleus"),
        node(523, "MPO", "Medial preoptic area"),
        node(515, "MPN", "Medial preoptic nucleus"),
        node(126, "PVp", "Periventricular hypothalamic nucleus, posterior part"),
    ]),
    # Outside HY: must never be hit, even by an exact acronym.
    node(9001, "ARC", "Some cortical thing named ARC"),
])]}


class ResolverTest(unittest.TestCase):
    def resolve(self, entry):
        nodes = ontology.flatten(GRAPH)
        return ontology.resolve_entry(entry, ontology.subtree(nodes, "HY"))

    def test_synonym_acronym_resolves_inside_scope_only(self):
        r = self.resolve({"id": "arc", "label": "Arcuate nucleus of the hypothalamus",
                          "abbreviation": "ARC", "synonyms": ["ARC", "ARH", "arcuate hypothalamic nucleus"]})
        self.assertEqual(r["status"], "resolved")
        self.assertEqual(r["allen_id"], 223)
        self.assertEqual(r["allen_name"], "Arcuate hypothalamic nucleus")

    def test_tier1_wins_and_alternates_are_reported(self):
        r = self.resolve({"id": "mpo", "label": "Medial preoptic area / nucleus", "abbreviation": "MPO",
                          "synonyms": ["medial preoptic area", "medial preoptic nucleus", "MPN"]})
        self.assertEqual(r["allen_id"], 523)
        self.assertEqual([a["allen_id"] for a in r["alternates"]], [515])

    def test_synonym_only_ambiguity_is_not_forced(self):
        r = self.resolve({"id": "x", "label": "Preoptic thing", "abbreviation": "PX",
                          "synonyms": ["medial preoptic area", "medial preoptic nucleus"]})
        self.assertEqual(r["status"], "ambiguous")
        self.assertIsNone(r["allen_id"])

    def test_no_partial_or_prefix_matching(self):
        r = self.resolve({"id": "pvp", "label": "Periventricular hypothalamic nucleus",
                          "abbreviation": "PVa/PVp", "synonyms": ["periventricular nucleus"]})
        self.assertEqual(r["status"], "unresolved")
        self.assertIsNone(r["allen_id"])

    def test_case_sensitive_acronyms(self):
        r = self.resolve({"id": "y", "label": "nothing", "abbreviation": "arh", "synonyms": []})
        self.assertEqual(r["status"], "unresolved")

    def test_scope_includes_hy_itself(self):
        r = self.resolve({"id": "hy", "label": "Hypothalamus", "abbreviation": "HY", "synonyms": []})
        self.assertEqual(r["allen_id"], 1097)

    def test_resolve_regions_requires_scope(self):
        with self.assertRaises(ValueError):
            ontology.resolve_regions({}, {"msg": [node(997, "root", "root")]})


if __name__ == "__main__":
    unittest.main()
