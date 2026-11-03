"""The five baselines on small synthetic sets with known answers."""

import _path  # noqa: F401
import os
import unittest

import synthetic
from syn4rch.models.baselines import (
    BaselineUnavailable,
    EdgeTable,
    GlobalMajority,
    KGPathSign,
    NEUROTRANSMITTER_VALENCE_RULES,
    NeurotransmitterSign,
    OutcomeMajority,
    PopulationOutcomeMajority,
)
from syn4rch.models.task import Features, build_dataset


def F(pop="syn:p1", outcome="syn:o1", valence="activation", pop_id="__same__"):
    return Features(
        population_id=pop if pop_id == "__same__" else pop_id,
        population_key=pop or "raw:x",
        valence=valence,
        outcome_id=outcome,
        species="mus_musculus",
        acuity="acute",
    )


def examples(spec):
    """spec: list of (population_id, outcome_id, direction)."""
    records = [
        synthetic.good_record(
            population_id=p, outcome_id=o, direction=d, source_locator="Fig %d" % i
        )
        for i, (p, o, d) in enumerate(spec)
    ]
    return build_dataset(records).examples


def stub_populations(entries):
    def load(kind):
        return entries if kind == "populations" else {}
    return load


class TestGlobalMajority(unittest.TestCase):
    def test_predicts_the_mode_everywhere_with_empirical_probabilities(self):
        train = examples(
            [("syn:p1", "syn:o1", "increase")] * 5
            + [("syn:p1", "syn:o1", "decrease")] * 3
        )
        model = GlobalMajority().fit(train)
        self.assertEqual(model.predict([F(), F(outcome="other")]), ["increase"] * 2)
        proba = model.predict_proba([F()])[0]
        self.assertAlmostEqual(proba["increase"], 5 / 8)
        self.assertAlmostEqual(proba["decrease"], 3 / 8)
        self.assertAlmostEqual(proba["no_change"], 0.0)

    def test_ties_break_on_label_order(self):
        train = examples(
            [("syn:p1", "syn:o1", "no_change"), ("syn:p1", "syn:o1", "decrease")]
        )
        self.assertEqual(GlobalMajority().fit(train).predict([F()]), ["decrease"])

    def test_fitted_on_nothing_it_abstains(self):
        model = GlobalMajority().fit([])
        self.assertEqual(model.predict([F()]), [None])
        self.assertEqual(model.predict_proba([F()]), [None])


class TestOutcomeMajority(unittest.TestCase):
    def test_per_outcome_mode_and_global_backoff(self):
        train = examples(
            [("syn:p1", "syn:oA", "decrease")] * 2
            + [("syn:p1", "syn:oA", "increase")]
            + [("syn:p1", "syn:oB", "increase")] * 3
        )
        model = OutcomeMajority().fit(train)
        self.assertEqual(
            model.predict([F(outcome="syn:oA"), F(outcome="syn:oB"), F(outcome="syn:oC")]),
            ["decrease", "increase", "increase"],
        )
        self.assertEqual(model.backoff_counts, {"outcome": 2, "global": 1})


class TestPopulationOutcomeMajority(unittest.TestCase):
    def setUp(self):
        train = examples(
            [("syn:p1", "syn:oA", "decrease")] * 2
            + [("syn:p2", "syn:oA", "increase")] * 3
            + [("syn:p1", "syn:oB", "no_change")] * 4
        )
        self.model = PopulationOutcomeMajority().fit(train)

    def test_each_rung_of_the_backoff_chain(self):
        cases = [
            (F(pop="syn:p1", outcome="syn:oA"), "decrease", "population_outcome"),
            (F(pop="syn:p2", outcome="syn:oB"), "no_change", "outcome"),
            (F(pop="syn:p2", outcome="syn:oZ"), "increase", "population"),
            (F(pop="syn:p9", outcome="syn:oZ"), "no_change", "global"),
        ]
        for features, expected, level in cases:
            counts, used = self.model.resolve(features)
            self.assertEqual(used, level, features)
            self.assertEqual(self.model.predict([features]), [expected], features)

    def test_outcome_is_tried_before_population(self):
        # p2 has a population marginal (increase) and oB an outcome marginal
        # (no_change): the documented chain uses the outcome.
        self.assertEqual(self.model.resolve(F(pop="syn:p2", outcome="syn:oB"))[1], "outcome")

    def test_empty_training_set_abstains(self):
        model = PopulationOutcomeMajority().fit([])
        self.assertEqual(model.predict([F()]), [None])
        self.assertEqual(model.resolve(F())[1], "abstain")

    def test_describe_publishes_the_chain(self):
        chain = self.model.describe()["backoff_chain"]
        self.assertEqual(chain[-2:], ["global", "abstain"])
        self.assertEqual(chain[0], "(population_key, outcome_id)")


class TestNeurotransmitterSign(unittest.TestCase):
    VOCAB = {
        "p_gaba": {"neurotransmitter": ["GABA"]},
        "p_glu": {"neurotransmitter": ["glutamate", "peptidergic"]},
        "p_mixed": {"neurotransmitter": ["mixed"]},
        "p_unknown": {"neurotransmitter": ["unknown"]},
        "p_silent": {"label": "no transmitter field at all"},
        "p_empty": {"neurotransmitter": []},
        "p_conflict": {"neurotransmitter": ["GABA", "glutamate"]},
        "p_da": {"neurotransmitter": ["dopamine"]},
    }

    def model(self):
        return NeurotransmitterSign(vocab_loader=stub_populations(self.VOCAB)).fit([])

    def test_signed_rules(self):
        m = self.model()
        self.assertEqual(
            m.predict([
                F(pop="p_gaba", valence="activation"),
                F(pop="p_gaba", valence="inhibition"),
                F(pop="p_gaba", valence="ablation"),
                F(pop="p_glu", valence="activation"),
                F(pop="p_glu", valence="inhibition"),
            ]),
            ["decrease", "increase", "increase", "increase", "decrease"],
        )

    def test_abstains_when_the_vocabulary_is_silent(self):
        m = self.model()
        for pop in ("p_silent", "p_empty", "p_unknown", "p_mixed"):
            self.assertEqual(m.predict([F(pop=pop)]), [None], pop)
        self.assertEqual(
            m.explain(F(pop="p_silent"))["abstain_reason"],
            "vocabulary_states_no_transmitter",
        )

    def test_abstains_when_the_vocabulary_is_entirely_empty(self):
        m = NeurotransmitterSign(vocab_loader=synthetic.stub_vocab_loader()).fit([])
        self.assertTrue(m.available)
        self.assertEqual(m.predict([F(pop="arc_agrp")]), [None])
        self.assertEqual(
            m.explain(F(pop="arc_agrp"))["abstain_reason"],
            "population_not_in_vocabulary",
        )

    def test_never_reads_the_transmitter_from_the_record(self):
        # good_record carries neurotransmitter=["unknown"]; even a record
        # claiming GABA must not be used when the vocabulary is silent.
        record = synthetic.good_record(population_id="p_silent", neurotransmitter=["GABA"])
        ex = build_dataset([record]).examples
        self.assertEqual(self.model().predict([e.features for e in ex]), [None])

    def test_other_abstentions(self):
        m = self.model()
        self.assertIsNone(m.predict([F(pop="p_gaba", valence="unclear")])[0])
        self.assertIsNone(m.predict([F(pop="p_conflict")])[0])
        self.assertIsNone(m.predict([F(pop="p_da")])[0])
        self.assertIsNone(m.predict([F(pop="raw:x", pop_id=None)])[0])
        self.assertEqual(
            m.explain(F(pop="p_conflict"))["abstain_reason"],
            "co_transmitters_imply_opposite_signs",
        )

    def test_emits_no_probability_and_never_predicts_no_change(self):
        m = self.model()
        self.assertIsNone(m.predict_proba([F(pop="p_gaba")]))
        self.assertNotIn("no_change", set(NEUROTRANSMITTER_VALENCE_RULES.values()))

    def test_a_vocabulary_that_cannot_load_makes_it_unavailable(self):
        def broken(kind):
            raise OSError("vocab missing")
        m = NeurotransmitterSign(vocab_loader=broken)
        self.assertFalse(m.available)
        self.assertIn("vocab missing", m.availability.reason)

    def test_the_real_vocabulary_states_agrp_is_gabaergic(self):
        m = NeurotransmitterSign()
        if not m.available:  # pragma: no cover
            self.skipTest("vocabulary files not present")
        self.assertEqual(m.transmitters("arc_agrp"), ["GABA"])
        self.assertIsNone(m.fit([]).predict([F(pop="avpv_kisspeptin")])[0])


class TestKGPathSign(unittest.TestCase):
    def setUp(self):
        self._env = os.environ.pop("SYN4RCH_KG_EDGES", None)

    def tearDown(self):
        if self._env is not None:
            os.environ["SYN4RCH_KG_EDGES"] = self._env

    def test_reports_itself_unavailable_with_a_reason_and_probes(self):
        m = KGPathSign()
        self.assertFalse(m.available)
        self.assertIn("no knowledge-graph edge table", m.availability.reason)
        self.assertIn("will not be invented", m.availability.reason)
        probes = " ".join(m.availability.probes)
        self.assertIn("SYN4RCH_KG_EDGES", probes)
        self.assertIn("syn4rch.kg.edges", probes)
        self.assertIn("kg_edges.jsonl", probes)

    def test_unavailable_baseline_abstains_and_can_be_made_to_raise(self):
        m = KGPathSign().fit(examples([("syn:p1", "syn:o1", "increase")]))
        self.assertEqual(m.predict([F(), F()]), [None, None])
        self.assertIsNone(m.predict_proba([F()]))
        with self.assertRaises(BaselineUnavailable):
            m.require_available()

    def test_probing_creates_no_file(self):
        from syn4rch.models.task import repo_root
        path = os.path.join(repo_root(), "data", "inferred", "kg_edges.jsonl")
        KGPathSign()
        self.assertFalse(os.path.exists(path))

    def test_an_empty_table_is_still_unavailable(self):
        self.assertFalse(KGPathSign(edge_table=EdgeTable(edges=[], source="synthetic:empty")).available)

    def test_sign_propagation_over_a_supplied_synthetic_table(self):
        # FABRICATED test edges, exercising the algorithm only.
        table = EdgeTable(
            edges=[
                {"source": "syn:p1", "target": "syn:mid", "sign": -1},
                {"source": "syn:mid", "target": "syn:o1", "sign": -1},
                {"source": "syn:p1", "target": "syn:o2", "sign": 1},
                {"source": "syn:p1", "target": "syn:a", "sign": 1},
                {"source": "syn:p1", "target": "syn:b", "sign": -1},
                {"source": "syn:a", "target": "syn:o3", "sign": 1},
                {"source": "syn:b", "target": "syn:o3", "sign": 1},
            ],
            source="synthetic:test-table",
        )
        m = KGPathSign(edge_table=table).fit([])
        self.assertTrue(m.available)
        self.assertEqual(
            m.predict([
                F(outcome="syn:o1"),
                F(outcome="syn:o1", valence="inhibition"),
                F(outcome="syn:o2"),
                F(outcome="syn:o3"),        # shortest paths disagree
                F(outcome="syn:none"),      # no path
            ]),
            ["increase", "decrease", "increase", None, None],
        )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
