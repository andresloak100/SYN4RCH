"""The signalling graph: edges are counts of records, and nothing else.

Every assertion here is about a number being a count of rows, a disagreement
staying visible, or an unknown staying unknown.
"""

import _path  # noqa: F401
import json
import os
import tempfile
import unittest

from syn4rch.atlas import signaling

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REAL_RECORDS = os.path.join(ROOT, "data", "curated", "records.jsonl")
POPS = os.path.join(ROOT, "data", "vocab", "populations.json")
OUTCOMES = os.path.join(ROOT, "data", "vocab", "outcomes.json")


def rec(**kw):
    """A record carrying only the fields the derivation reads."""
    r = {
        "record_id": kw.pop("record_id", "r%04d" % rec.n),
        "study_id": "doi:10.1/s1",
        "last_author": "Smith AB",
        "last_author_affiliation": "Some University",
        "population_id": "arc_agrp",
        "population_raw": "ARC AgRP neurons",
        "outcome_id": "food_intake",
        "method": "chemogenetic",
        "valence": "activation",
        "direction": "increase",
        "projection_target_id": None,
    }
    rec.n += 1
    r.update(kw)
    return r


rec.n = 0


class SignalingTest(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        rec.n = 0

    def graph(self, rows, pops=None, outs=None):
        p = os.path.join(self.d, "records.jsonl")
        with open(p, "w", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
        return signaling.build_graph(p, pops or os.path.join(self.d, "no_pops.json"),
                                     outs or os.path.join(self.d, "no_outcomes.json"))

    def edge(self, g, population_id, outcome_id, axis=None):
        hits = [e for e in g["edges"] if e["population_id"] == population_id
                and e["outcome_id"] == outcome_id and (axis is None or e["axis"] == axis)]
        self.assertEqual(len(hits), 1, "expected exactly one edge, got %d" % len(hits))
        return hits[0]

    # --- an edge is always a count of records ------------------------------

    def test_no_records_no_edge(self):
        g = self.graph([])
        self.assertEqual(g["edges"], [])
        self.assertEqual(g["n_records_used"], 0)
        self.assertEqual(g["n_edges"], 0)

    def test_every_edge_carries_at_least_one_record(self):
        g = self.graph([rec(), rec(outcome_id="body_weight", direction="no_change")])
        self.assertEqual(len(g["edges"]), 2)
        for e in g["edges"]:
            self.assertGreaterEqual(e["n_records"], 1)
            self.assertEqual(sum(e["directions"].values()), e["n_records"])
            self.assertGreaterEqual(e["n_studies"], 1)

    def test_a_record_with_no_direction_builds_no_edge(self):
        g = self.graph([rec(direction=None)])
        self.assertEqual(g["edges"], [])
        self.assertEqual(g["excluded"], {signaling.EX_NO_DIRECTION: 1})

    def test_record_ids_are_capped_and_the_cap_is_declared(self):
        g = self.graph([rec() for _ in range(signaling.RECORD_IDS_CAP + 7)])
        e = g["edges"][0]
        self.assertEqual(e["n_records"], signaling.RECORD_IDS_CAP + 7)
        self.assertEqual(len(e["record_ids"]), signaling.RECORD_IDS_CAP)
        self.assertTrue(e["record_ids_truncated"])
        self.assertEqual(g["record_ids_cap"], signaling.RECORD_IDS_CAP)

    # --- the unclear axis stays unclear ------------------------------------

    def test_unclear_axis_is_never_coerced(self):
        g = self.graph([rec(method="conditional_knockout", valence="unclear", direction="decrease"),
                        rec(method="conditional_overexpression", valence="unclear")])
        axes = sorted({e["axis"] for e in g["edges"]})
        self.assertEqual(axes, ["unclear"])
        for e in g["edges"]:
            self.assertEqual(e["axis_class"], "molecular")
            self.assertNotIn(e["axis"], ("activation", "inhibition"))
            self.assertIn("not established", e["axis_phrase"])

    def test_axis_separates_edges_rather_than_pooling_them(self):
        g = self.graph([rec(valence="activation"),
                        rec(valence="inhibition", direction="decrease"),
                        rec(method="conditional_knockout", valence="unclear", direction="no_change")])
        self.assertEqual(len(g["edges"]), 3)
        self.assertEqual(sorted(e["axis"] for e in g["edges"]),
                         ["activation", "inhibition", "unclear"])

    def test_stored_polarity_wins_over_valence_when_it_arrives(self):
        # Written before perturbation_polarity existed; must pick it up with no
        # code change, and must not read valence once it is there.
        g = self.graph([rec(valence="unclear", method="conditional_knockout",
                            perturbation_polarity="decrease",
                            perturbation_class="molecular")])
        e = g["edges"][0]
        self.assertEqual(e["axis"], "decrease")
        self.assertEqual(e["axis_class"], "molecular")
        self.assertEqual(g["axis_source"],
                         "perturbation_polarity when the record carries it, otherwise valence")

    def test_axis_phrase_never_shows_the_bare_polarity_value(self):
        # The polarity trap: axis "increase" and direction "increase" are
        # different statements and must not read the same.
        self.assertEqual(signaling.axis_phrase("increase", "circuit_activity"), "more firing")
        self.assertEqual(signaling.axis_phrase("decrease", "molecular"),
                         "less of the target gene product")
        self.assertNotEqual(signaling.axis_phrase("increase", "molecular"), "increase")
        self.assertIn("unknown", signaling.axis_phrase("wat", "nope").replace("recorded as", "unknown"))

    def test_axis_class_from_method_matches_the_schema_rule(self):
        # The local map mirrors syn4rch.schema.rules; this fails if they drift.
        try:
            from syn4rch.schema.rules import perturbation_class_from_method
        except ImportError:  # pragma: no cover - rule module is optional here
            self.skipTest("syn4rch.schema.rules has no perturbation_class_from_method")
        for method, cls in signaling.METHOD_CLASS.items():
            got = perturbation_class_from_method(method)
            self.assertIsNotNone(got, method)
            self.assertEqual(getattr(got, "value", got), cls, method)

    def test_stored_perturbation_class_wins_over_the_method_map(self):
        g = self.graph([rec(method="chemogenetic", perturbation_class="pharmacological")])
        self.assertEqual(g["edges"][0]["axis_class"], "pharmacological")

    # --- consensus, agreement, contested ----------------------------------

    def test_unanimous_edge(self):
        g = self.graph([rec() for _ in range(5)])
        e = g["edges"][0]
        self.assertEqual(e["directions"], {"increase": 5})
        self.assertEqual(e["consensus"], "increase")
        self.assertEqual(e["agreement"], 1.0)
        self.assertFalse(e["contested"])

    def test_contested_edge_by_the_stated_rule(self):
        # 4 increase, 3 decrease, 3 no_change: agreement 0.4, 10 records.
        rows = ([rec(study_id="doi:10.1/a%d" % i, last_author="A%d" % i) for i in range(4)]
                + [rec(direction="decrease", study_id="doi:10.1/b%d" % i, last_author="B%d" % i) for i in range(3)]
                + [rec(direction="no_change", study_id="doi:10.1/c%d" % i, last_author="C%d" % i) for i in range(3)])
        g = self.graph(rows)
        e = g["edges"][0]
        self.assertEqual(e["directions"], {"decrease": 3, "increase": 4, "no_change": 3})
        self.assertEqual(e["n_records"], 10)
        self.assertEqual(e["consensus"], "increase")
        self.assertEqual(e["agreement"], 0.4)
        self.assertTrue(e["contested"])
        self.assertEqual(g["contested_rule"], signaling.CONTESTED_RULE)

    def test_tie_leaves_no_consensus(self):
        g = self.graph([rec(), rec(study_id="doi:10.1/s2"),
                        rec(direction="decrease"),
                        rec(study_id="doi:10.1/s2", direction="decrease")])
        e = g["edges"][0]
        self.assertEqual(e["directions"], {"decrease": 2, "increase": 2})
        self.assertIsNone(e["consensus"])
        self.assertEqual(e["agreement"], 0.5)
        self.assertTrue(e["contested"])

    def test_a_thin_majority_is_a_consensus_but_a_contested_one(self):
        # consensus is the modal direction, not a claim of a majority; the low
        # agreement is what the page has to show next to it.
        g = self.graph([rec(direction="decrease"), rec(direction="decrease"),
                        rec(), rec(direction="no_change")])
        e = g["edges"][0]
        self.assertEqual(e["consensus"], "decrease")
        self.assertEqual(e["agreement"], 0.5)
        self.assertTrue(e["contested"])

    def test_two_records_are_too_few_to_call_contested(self):
        g = self.graph([rec(), rec(direction="decrease")])
        e = g["edges"][0]
        self.assertEqual(e["agreement"], 0.5)
        self.assertEqual(e["n_records"], 2)
        self.assertFalse(e["contested"])
        self.assertIsNone(e["consensus"])

    def test_high_agreement_is_not_contested_but_a_replicated_minority_is_flagged(self):
        # 10 increase from 2 studies, 2 no_change from 2 other studies:
        # agreement 0.83, so not contested -- but the minority replicated.
        rows = ([rec(study_id="doi:10.1/a") for _ in range(5)]
                + [rec(study_id="doi:10.1/b") for _ in range(5)]
                + [rec(study_id="doi:10.1/c", direction="no_change"),
                   rec(study_id="doi:10.1/d", direction="no_change")])
        e = self.graph(rows)["edges"][0]
        self.assertEqual(e["agreement"], 0.833)
        self.assertFalse(e["contested"])
        self.assertTrue(e["minority_replicated"])

    def test_single_study_minority_is_not_flagged(self):
        rows = [rec(study_id="doi:10.1/a") for _ in range(5)] + [rec(study_id="doi:10.1/a", direction="no_change")]
        e = self.graph(rows)["edges"][0]
        self.assertFalse(e["minority_replicated"])
        self.assertEqual(e["n_studies"], 1)

    # --- study and lab identity -------------------------------------------

    def test_n_labs_does_not_count_a_record_with_no_last_author(self):
        rows = [rec(study_id="doi:10.1/a", last_author="Smith AB"),
                rec(study_id="doi:10.1/b", last_author=None),
                rec(study_id="doi:10.1/c", last_author="")]
        e = self.graph(rows)["edges"][0]
        self.assertEqual(e["n_records"], 3)
        self.assertEqual(e["n_studies"], 3)
        self.assertEqual(e["n_labs"], 1)

    def test_two_unknown_labs_are_not_lumped_into_one(self):
        rows = [rec(study_id="doi:10.1/a", last_author=None),
                rec(study_id="doi:10.1/b", last_author=None)]
        e = self.graph(rows)["edges"][0]
        self.assertEqual(e["n_labs"], 0)
        self.assertEqual(e["n_studies"], 2)

    def test_same_lab_two_studies_counts_one_lab(self):
        rows = [rec(study_id="doi:10.1/a", last_author="Smith AB", last_author_affiliation="U1"),
                rec(study_id="doi:10.1/b", last_author="smith ab", last_author_affiliation="U1")]
        e = self.graph(rows)["edges"][0]
        self.assertEqual((e["n_studies"], e["n_labs"]), (2, 1))

    # --- the partition -----------------------------------------------------

    def test_region_level_records_are_excluded_not_forced_onto_a_population(self):
        g = self.graph([rec(), rec(population_id=None, population_raw="ARC")])
        self.assertEqual(g["n_records_used"], 1)
        self.assertEqual(g["excluded"], {signaling.EX_NO_POP: 1})
        self.assertEqual([p["id"] for p in g["populations"] if p["n_records"]], ["arc_agrp"])

    def test_used_plus_excluded_equals_input(self):
        g = self.graph([rec(), rec(population_id=None), rec(outcome_id=None),
                        rec(direction=None), rec(population_id=None, outcome_id=None)])
        self.assertEqual(g["n_records_in"], 5)
        self.assertEqual(g["n_records_used"] + sum(g["excluded"].values()), g["n_records_in"])

    def test_an_excluded_readout_still_appears_with_its_excluded_count(self):
        g = self.graph([rec(population_id=None, outcome_id="heart_rate")])
        o = [x for x in g["outcomes"] if x["id"] == "heart_rate"][0]
        self.assertEqual((o["n_records"], o["n_records_excluded"]), (0, 1))

    # --- systems are a layout grouping, and total -------------------------

    def test_every_outcome_in_the_real_dataset_has_a_system(self):
        with open(REAL_RECORDS, encoding="utf-8") as fh:
            ids = {json.loads(l)["outcome_id"] for l in fh if l.strip()}
        for oid in sorted(i for i in ids if i):
            self.assertIn(signaling.system_of(oid), signaling.SYSTEM_ORDER, oid)

    def test_the_real_dataset_has_no_unclassified_outcome(self):
        """The protection the hard assert used to give, as a reported number.

        This replaced an `assert` with a declared `unclassified` bucket plus
        this list, because OUTCOME_SYSTEM is a presentational layout aid and a
        missing entry costs an outcome sitting in the wrong place on a page,
        not a build. The gap is still caught here; it just arrives as a number.
        """
        g = signaling.build_graph(REAL_RECORDS, POPS, OUTCOMES)
        self.assertEqual(g["outcomes_without_a_declared_system"], [])

    def test_an_unmapped_outcome_is_reported_rather_than_fatal(self):
        """The behaviour change itself. This used to raise AssertionError.

        It fired twice in two days as the vocabulary grew -- on `nrem_sleep`
        from the sleep split, and in the same shape as `condition_key` going to
        five slots and taking models/task.py:555 with it. A hardcoded map that
        must track a growing vocabulary will keep doing that.
        """
        with open(REAL_RECORDS, encoding="utf-8") as fh:
            rows = [json.loads(l) for l in fh if l.strip()]
        rows[0]["outcome_id"] = "a_readout_the_layout_has_never_heard_of"
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "records.jsonl")
            with open(path, "w", encoding="utf-8") as fh:
                for row in rows:
                    fh.write(json.dumps(row) + "\n")
            g = signaling.build_graph(path, POPS, OUTCOMES)
        self.assertIn(
            "a_readout_the_layout_has_never_heard_of",
            g["outcomes_without_a_declared_system"],
        )
        self.assertTrue(g["edges"], "the graph must still build")

    def test_unclassified_is_a_legal_layout_slot(self):
        self.assertIn(signaling.UNCLASSIFIED_SYSTEM, signaling.SYSTEM_ORDER)
        self.assertEqual(signaling.SYSTEM_ORDER[-1], signaling.UNCLASSIFIED_SYSTEM)

    def test_the_whole_sleep_family_shares_one_system(self):
        """The split must change no edge's classification."""
        family = [
            "sleep", "nrem_sleep", "rem_sleep", "wakefulness",
            "sleep_latency", "sleep_fragmentation",
        ]
        systems = {signaling.system_of(o) for o in family}
        self.assertEqual(systems, {"behaviour"}, systems)

    def test_every_mapped_system_is_in_the_layout_order(self):
        for oid, system in signaling.OUTCOME_SYSTEM.items():
            self.assertIn(system, signaling.SYSTEM_ORDER, oid)
            self.assertIn(system, signaling.SYSTEM_LABEL, oid)

    # --- the real dataset --------------------------------------------------

    def test_real_dataset_is_internally_consistent(self):
        g = signaling.build_graph(REAL_RECORDS, POPS, OUTCOMES)
        self.assertEqual(g["n_records_in"], g["n_records_used"] + sum(g["excluded"].values()))
        self.assertTrue(g["edges"])
        for e in g["edges"]:
            self.assertEqual(sum(e["directions"].values()), e["n_records"])
            self.assertLessEqual(e["n_labs"], e["n_studies"])
            self.assertLessEqual(e["n_studies"], e["n_records"])
            self.assertIn(signaling.system_of(e["outcome_id"]), signaling.SYSTEM_ORDER)
            if e["consensus"] is None:
                self.assertTrue(e["contested"] or e["n_records"] < 3)
        # Every population referenced by an edge is a node, and the zero-record
        # vocabulary populations are listed too so the page can say so.
        ids = {p["id"] for p in g["populations"]}
        for e in g["edges"]:
            self.assertIn(e["population_id"], ids)
        self.assertGreater(len(ids), g["n_populations_with_records"])
        self.assertTrue(all(h for h in g["_honesty"]))
        self.assertEqual(len(g["_honesty"]), len(set(g["_honesty"])))


if __name__ == "__main__":
    unittest.main()
