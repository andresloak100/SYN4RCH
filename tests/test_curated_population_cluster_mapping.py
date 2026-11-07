"""Invariants of the committed population -> cluster mapping (D-065).

The mapping is a Tier B annotation built from a 4 GB atlas that is not in the
repository, so these tests cannot recompute it. What they can do is hold the
artefact to the rule it claims to have been built by, and refuse the specific
ways it could quietly go wrong: a tier that the recorded evidence does not
support, a refusal with no reason, a mapping that outlived a vocabulary entry.
"""

import _path  # noqa: F401
import json
import os
import unittest

from syn4rch.atlas.population_clusters import (
    GOOD_RANK,
    MARKER_FLOOR,
    MIN_CONTRAST,
    MIN_REGION_SHARE,
    KEEP_FRACTION,
    OUTCOMES,
    TIERS,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAPPING = os.path.join(ROOT, "data", "curated",
                       "population_cluster_mapping.json")
VOCAB = os.path.join(ROOT, "data", "vocab", "populations.json")


def mapping():
    with open(MAPPING, encoding="utf-8") as handle:
        return json.load(handle)


def vocab():
    with open(VOCAB, encoding="utf-8") as handle:
        return json.load(handle)["entries"]


class ProvenanceTest(unittest.TestCase):
    def test_the_atlas_is_identified_by_doi_and_organism(self):
        atlas = mapping()["atlas"]
        self.assertEqual(atlas["doi"], "10.1038/s42255-022-00657-y")
        self.assertEqual(atlas["organism"], "NCBITaxon:10090")

    def test_the_mapping_declares_itself_tier_B_and_not_a_training_target(self):
        """The standing constraint, written into the artefact rather than only
        into the prose around it."""
        doc = mapping()
        self.assertEqual(doc["tier"], "B")
        self.assertIs(doc["never_a_training_target"], True)

    def test_the_rule_records_that_it_was_declared_after_looking(self):
        """This is not pre-registration and must not be cited as if it were."""
        self.assertIs(mapping()["rule"]["declared_after_inspecting_the_atlas"],
                      True)

    def test_the_recorded_rule_parameters_match_the_code_that_applied_them(self):
        """If a constant moves, the committed artefact is stale, not merely
        inconsistent -- it was produced by different arithmetic."""
        rule = mapping()["rule"]
        self.assertEqual(rule["min_region_share_percent"], MIN_REGION_SHARE)
        self.assertEqual(rule["marker_floor"], MARKER_FLOOR)
        self.assertEqual(rule["keep_fraction"], KEEP_FRACTION)
        self.assertEqual(rule["min_contrast_for_tier_A_or_B"], MIN_CONTRAST)
        self.assertEqual(rule["rank_good_enough_for_tier_B"], GOOD_RANK)


class CoverageOfTheVocabularyTest(unittest.TestCase):
    def test_every_declared_population_has_an_entry(self):
        entries, declared = mapping()["entries"], vocab()
        self.assertEqual(sorted(entries), sorted(declared))

    def test_no_entry_outlives_the_vocabulary_entry_it_describes(self):
        extra = sorted(set(mapping()["entries"]) - set(vocab()))
        self.assertEqual(extra, [], "mappings for populations that no longer exist")

    def test_each_entry_keeps_the_nucleus_its_vocabulary_entry_declares(self):
        declared = vocab()
        for pid, row in mapping()["entries"].items():
            self.assertEqual(row["nucleus"], declared[pid].get("nucleus"), pid)


class EntryShapeTest(unittest.TestCase):
    def test_every_outcome_is_one_the_rule_can_return(self):
        for pid, row in mapping()["entries"].items():
            self.assertIn(row["outcome"], OUTCOMES, pid)

    def test_a_refusal_always_says_why(self):
        for pid, row in mapping()["entries"].items():
            if row["outcome"] != "mapped":
                self.assertTrue(row.get("reason"), pid)
                self.assertNotIn("clusters", row, pid)

    def test_a_mapping_carries_a_level_clusters_and_a_tier(self):
        for pid, row in mapping()["entries"].items():
            if row["outcome"] != "mapped":
                continue
            self.assertIn(row["level"], mapping()["atlas"]["levels"], pid)
            self.assertTrue(row["clusters"], pid)
            self.assertIn(row["tier"], TIERS, pid)

    def test_the_marker_is_one_the_vocabulary_declares(self):
        declared = vocab()
        for pid, row in mapping()["entries"].items():
            if row.get("marker") is None:
                continue
            self.assertIn(row["marker"],
                          declared[pid].get("marker_genes") or [], pid)

    def test_no_cluster_is_claimed_twice_within_one_population(self):
        for pid, row in mapping()["entries"].items():
            cl = row.get("clusters") or []
            self.assertEqual(len(cl), len(set(cl)), pid)

    def test_clusters_are_named_at_the_level_the_entry_declares(self):
        for pid, row in mapping()["entries"].items():
            for cluster in row.get("clusters") or []:
                self.assertTrue(cluster.startswith(row["level"] + "-"),
                                "%s: %s is not a %s cluster"
                                % (pid, cluster, row["level"]))

    def test_the_named_count_cannot_exceed_the_clusters_it_counts(self):
        for pid, row in mapping()["entries"].items():
            if row["outcome"] != "mapped":
                continue
            self.assertEqual(row["total_clusters"], len(row["clusters"]), pid)
            self.assertLessEqual(row["named_clusters"], row["total_clusters"], pid)


class TierConsistencyTest(unittest.TestCase):
    def test_tier_A_requires_every_cluster_name_to_carry_a_marker(self):
        """The one claim tier A makes that a reader will rely on."""
        for pid, row in mapping()["entries"].items():
            if row.get("tier") == "A":
                self.assertEqual(row["named_clusters"], row["total_clusters"], pid)

    def test_tier_A_requires_contrast_at_or_above_the_declared_minimum(self):
        for pid, row in mapping()["entries"].items():
            if row.get("tier") != "A":
                continue
            if row["contrast_is_infinite"]:
                continue
            self.assertGreaterEqual(row["contrast"], MIN_CONTRAST, pid)

    def test_tier_C_is_not_claimed_for_something_with_both_lines_of_support(self):
        for pid, row in mapping()["entries"].items():
            if row.get("tier") != "C":
                continue
            full_names = row["named_clusters"] == row["total_clusters"]
            self.assertFalse(full_names, "%s has full name support but tier C" % pid)

    def test_rank_is_recorded_even_when_it_is_poor(self):
        """D-065: Crh ranks 544th inside its own named cluster. The figure is
        kept precisely so a reader can see that it was not gated on."""
        rows = mapping()["entries"]
        self.assertEqual(rows["pvh_crh"]["tier"], "A")
        self.assertGreater(rows["pvh_crh"]["best_marker_rank"], 200)


class RefusalTest(unittest.TestCase):
    def test_a_cre_line_is_not_a_marker(self):
        """The regularity D-065 records: receptor- and TF-defined populations
        have real Cre lines and not enough transcript to mark a cluster."""
        rows = mapping()["entries"]
        for pid in ("dmh_lepr", "pvh_mc4r", "vmh_sf1"):
            self.assertEqual(rows[pid]["outcome"], "refused_below_floor", pid)
            self.assertTrue(vocab()[pid].get("cre_lines"),
                            "%s is only an interesting refusal because it has "
                            "a Cre line" % pid)

    def test_a_below_floor_refusal_names_the_best_it_found(self):
        for pid, row in mapping()["entries"].items():
            if row["outcome"] == "refused_below_floor":
                self.assertIn("best", row["reason"], pid)
                self.assertIn("cluster", row["reason"], pid)


class NearMissTest(unittest.TestCase):
    """D-070's asymmetry, applied to this artefact.

    A wrongly included cluster sits in `clusters` where a reader will see it.
    A wrongly excluded one would be recorded nowhere, so the near-misses are
    written down -- and held to not being members.
    """

    def test_a_near_miss_is_never_also_a_member(self):
        for pid, row in mapping()["entries"].items():
            members = set(row.get("clusters") or [])
            for nm in row.get("near_miss_excluded") or []:
                self.assertNotIn(nm["cluster"], members, pid)

    def test_every_near_miss_is_below_the_keep_threshold(self):
        """If one were at or above it, the rule would have kept it, and the
        artefact would be inconsistent with the code that built it."""
        for pid, row in mapping()["entries"].items():
            for nm in row.get("near_miss_excluded") or []:
                self.assertLess(nm["fraction_of_best"], KEEP_FRACTION, pid)

    def test_the_note_says_they_are_not_members(self):
        self.assertIn("NOT members", mapping()["near_miss_note"])

    def test_the_pvh_trh_near_miss_that_moves_the_result_is_recorded(self):
        """D-068's sensitivity note turns on this one: a Trh-named, 100%-PVH
        cluster excluded by 0.058, in one of only two rankable populations.
        Adding it takes the result from 2 of 9 rankable to 0 of 9."""
        rows = mapping()["entries"]["pvh_trh"].get("near_miss_excluded") or []
        got = {r["cluster"]: r["fraction_of_best"] for r in rows}
        self.assertIn("C185-12: Ebf1.Trh.GLU-2", got)
        self.assertAlmostEqual(got["C185-12: Ebf1.Trh.GLU-2"], 0.442, places=3)


class NestingTest(unittest.TestCase):
    def test_nesting_is_recorded_so_a_similarity_experiment_can_see_it(self):
        """Two populations sharing clusters have correlated vectors by
        construction; D-059's per-row diagnostics have to be able to tell."""
        doc = mapping()
        self.assertIn("nesting", doc)
        for pair in doc["nesting"]:
            self.assertIn(pair["relation"], ("subset", "superset", "partial"))
            self.assertIn(pair["a"], doc["entries"])
            self.assertIn(pair["b"], doc["entries"])

    def test_every_recorded_nesting_really_shares_cells(self):
        doc = mapping()
        for pair in doc["nesting"]:
            a = set(doc["entries"][pair["a"]].get("leaf_clusters") or [])
            b = set(doc["entries"][pair["b"]].get("leaf_clusters") or [])
            self.assertTrue(a & b, "%s/%s share no cells" % (pair["a"], pair["b"]))

    def test_every_mapped_entry_carries_its_leaf_clusters(self):
        """Without them the overlap test silently degrades to name matching."""
        for pid, row in mapping()["entries"].items():
            if row["outcome"] != "mapped":
                continue
            self.assertTrue(row.get("leaf_clusters"), pid)
            self.assertEqual(row["n_leaf_clusters"], len(row["leaf_clusters"]), pid)
            for leaf in row["leaf_clusters"]:
                self.assertTrue(leaf.startswith("C465-"), "%s: %s" % (pid, leaf))

    def test_the_cross_level_overlap_that_name_matching_missed_is_declared(self):
        doc = mapping()
        listed = {frozenset((p["a"], p["b"])) for p in doc["nesting"]}
        self.assertIn(frozenset(("arc_agrp", "arc_tbx3")), listed)
        a = set(doc["entries"]["arc_agrp"]["clusters"])
        b = set(doc["entries"]["arc_tbx3"]["clusters"])
        self.assertEqual(a & b, set(),
                         "this pair must share no cluster name, or it no longer "
                         "demonstrates what name matching missed")

    def test_no_unrecorded_overlap_exists(self):
        """The invariant that matters: a pair sharing CELLS and not listed
        would be invisible to anything downstream.

        Compared on `leaf_clusters`, because an earlier version of this test
        compared cluster names and passed while ten of eleven overlapping pairs
        went unrecorded -- every pair mapped at two different levels.
        """
        doc = mapping()
        sets = {pid: set(row["leaf_clusters"])
                for pid, row in doc["entries"].items()
                if row["outcome"] == "mapped"}
        listed = {frozenset((p["a"], p["b"])) for p in doc["nesting"]}
        for a in sorted(sets):
            for b in sorted(sets):
                if a >= b or not (sets[a] & sets[b]):
                    continue
                self.assertIn(frozenset((a, b)), listed,
                              "%s and %s share clusters but are not recorded" % (a, b))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
