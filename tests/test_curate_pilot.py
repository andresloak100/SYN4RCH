"""Pilot selection: reproducible, unbiased within a stratum, honestly labelled."""

import _path  # noqa: F401
import json
import os
import tempfile
import unittest

from syn4rch.curate import pilot


def unit(study, pop, outcome="food_intake", section="Results", locator="Results para 1",
         text="AgRP activation increased food intake.", method="chemo_hm3dq", score=5):
    return {"study_id": study, "candidate_population_id": pop,
            "candidate_method_id": method, "candidate_outcome_id": outcome,
            "section_label": section, "locator": locator, "passage_text": text,
            "score": score, "unit_id": "u-%s-%s-%s" % (study, pop, outcome),
            "evidence_terms": {}}


class EligibilityTest(unittest.TestCase):
    def test_only_result_sections_make_a_paper_eligible(self):
        units = [
            unit("doi:1", "arc_agrp", section="Discussion"),
            unit("doi:1", "arc_agrp", section="Introduction"),
            unit("doi:2", "arc_agrp", section="Results"),
            unit("doi:3", "arc_agrp", section="Methods"),
        ]
        elig = pilot.eligible(units, ["arc_agrp"])
        self.assertEqual(sorted(elig), ["doi:2", "doi:3"])

    def test_population_not_requested_does_not_make_a_paper_eligible(self):
        units = [unit("doi:1", "lha_orexin"), unit("doi:2", "arc_pomc")]
        elig = pilot.eligible(units, ["arc_agrp", "arc_pomc"])
        self.assertEqual(sorted(elig), ["doi:2"])
        self.assertEqual(elig["doi:2"]["covers"], {"arc_pomc"})

    def test_eligible_keeps_all_result_units_of_an_eligible_paper(self):
        # A paper selected for AgRP also carries its other observations: once the
        # passage is being read, anything it reports is admissible.
        units = [unit("doi:1", "arc_agrp"), unit("doi:1", "lha_orexin", outcome="sleep")]
        elig = pilot.eligible(units, ["arc_agrp"])
        self.assertEqual(len(elig["doi:1"]["units"]), 2)


class SelectionTest(unittest.TestCase):
    def _elig(self, n_agrp=10, n_pomc=10, n_both=5):
        units = []
        for i in range(n_agrp):
            units.append(unit("doi:a%d" % i, "arc_agrp"))
        for i in range(n_pomc):
            units.append(unit("doi:p%d" % i, "arc_pomc"))
        for i in range(n_both):
            units.append(unit("doi:b%d" % i, "arc_agrp"))
            units.append(unit("doi:b%d" % i, "arc_pomc"))
        return pilot.eligible(units, ["arc_agrp", "arc_pomc"])

    def test_same_seed_gives_same_selection(self):
        elig = self._elig()
        a = pilot.select(elig, n_papers=12, populations=["arc_agrp", "arc_pomc"], seed="s1")
        b = pilot.select(elig, n_papers=12, populations=["arc_agrp", "arc_pomc"], seed="s1")
        self.assertEqual(a["chosen"], b["chosen"])

    def test_different_seed_gives_a_different_selection(self):
        elig = self._elig()
        a = pilot.select(elig, n_papers=12, populations=["arc_agrp", "arc_pomc"], seed="s1")
        b = pilot.select(elig, n_papers=12, populations=["arc_agrp", "arc_pomc"], seed="s2")
        self.assertNotEqual(a["chosen"], b["chosen"])

    def test_quotas_sum_exactly_to_requested_papers(self):
        elig = self._elig(n_agrp=52, n_pomc=42, n_both=24)
        for n in (1, 7, 13, 50, 99):
            sel = pilot.select(elig, n_papers=n, populations=["arc_agrp", "arc_pomc"], seed="s")
            self.assertEqual(sum(sel["quotas"].values()), n, n)
            self.assertEqual(len(sel["chosen"]), n, n)
            self.assertEqual(len(set(sel["chosen"])), n, n)

    def test_cannot_select_more_papers_than_exist(self):
        elig = self._elig(n_agrp=3, n_pomc=2, n_both=0)
        sel = pilot.select(elig, n_papers=99, populations=["arc_agrp", "arc_pomc"], seed="s")
        self.assertEqual(len(sel["chosen"]), 5)
        self.assertEqual(sel["eligible_total"], 5)

    def test_strata_are_proportional(self):
        elig = self._elig(n_agrp=60, n_pomc=30, n_both=10)
        sel = pilot.select(elig, n_papers=10, populations=["arc_agrp", "arc_pomc"], seed="s")
        self.assertEqual(sel["quotas"]["arc_agrp"], 6)
        self.assertEqual(sel["quotas"]["arc_pomc"], 3)
        self.assertEqual(sel["quotas"]["arc_agrp+arc_pomc"], 1)

    def test_order_is_not_correlated_with_unit_count(self):
        # The richest paper must not be selected preferentially: that would
        # inflate the records-per-paper yield this pilot exists to measure.
        units = []
        for i in range(20):
            for _ in range(i + 1):  # doi:a19 has 20 units, doi:a0 has 1
                units.append(unit("doi:a%d" % i, "arc_agrp", outcome="o%d" % _))
        elig = pilot.eligible(units, ["arc_agrp"])
        sel = pilot.select(elig, n_papers=5, populations=["arc_agrp"], seed="s1")
        richest = "doi:a19"
        picked_across_seeds = [richest in pilot.select(
            elig, n_papers=5, populations=["arc_agrp"], seed="seed%d" % k)["chosen"]
            for k in range(12)]
        self.assertIn(False, picked_across_seeds)
        self.assertEqual(len(sel["chosen"]), 5)


class PackagingTest(unittest.TestCase):
    def test_stub_is_filled_from_harvest_and_labelled_unreviewed(self):
        paper = {"study_id": "doi:1", "pmid": "123", "pmcid": "PMC1", "doi": "10.1/x",
                 "title": "T", "journal": "J", "year": 2020, "first_author": "A",
                 "last_author": "B", "last_author_affiliation": "Aff",
                 "full_text_available": True}
        stub = pilot.build_stub(paper, "0.1.0")
        self.assertEqual(stub["pmid"], "123")
        self.assertEqual(stub["retrieved_from"], "europepmc_oa")
        self.assertEqual(stub["extraction_method"], "llm")
        self.assertEqual(stub["review_status"], "unreviewed")
        self.assertEqual(stub["evidence_class"], "causal_interventional")
        self.assertEqual(stub["schema_version"], "0.1.0")

    def test_stub_without_full_text_is_marked_pubmed_and_false(self):
        stub = pilot.build_stub({"study_id": "doi:1"}, "0.1.0")
        self.assertFalse(stub["full_text_available"])
        self.assertEqual(stub["retrieved_from"], "pubmed")

    def test_passages_are_deduped_and_leads_collected(self):
        units = [unit("doi:1", "arc_agrp", outcome="food_intake"),
                 unit("doi:1", "arc_agrp", outcome="body_weight")]
        elig = pilot.eligible(units, ["arc_agrp"])
        rows = pilot.package(["doi:1"], elig, {"doi:1": {"study_id": "doi:1"}},
                             n_batches=1, schema_version="0.1.0")[0]
        self.assertEqual(rows[0]["n_passages"], 1)
        self.assertEqual(len(rows[0]["passages"][0]["leads"]), 2)
        self.assertIn("not a record", rows[0]["_reminder"])

    def test_batches_are_balanced_and_lose_nothing(self):
        units = [unit("doi:%d" % i, "arc_agrp") for i in range(11)]
        elig = pilot.eligible(units, ["arc_agrp"])
        chosen = sorted(elig)
        batches = pilot.package(chosen, elig, {}, n_batches=5, schema_version="0.1.0")
        self.assertEqual(len(batches), 5)
        self.assertEqual(sum(len(b) for b in batches), 11)
        self.assertLessEqual(max(len(b) for b in batches) - min(len(b) for b in batches), 1)
        flat = {r["study_id"] for b in batches for r in b}
        self.assertEqual(flat, set(chosen))


class CliTest(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def p(self, n):
        return os.path.join(self.d, n)

    def test_missing_input_exits_nonzero(self):
        self.assertEqual(pilot.main(["--units", self.p("no.jsonl"),
                                     "--candidates", self.p("no2.jsonl"),
                                     "--out-dir", self.p("out")]), 1)

    def test_no_eligible_papers_exits_nonzero(self):
        with open(self.p("u.jsonl"), "w") as fh:
            fh.write(json.dumps(unit("doi:1", "lha_orexin")) + "\n")
        with open(self.p("c.jsonl"), "w") as fh:
            fh.write(json.dumps({"study_id": "doi:1"}) + "\n")
        self.assertEqual(pilot.main(["--units", self.p("u.jsonl"),
                                     "--candidates", self.p("c.jsonl"),
                                     "--out-dir", self.p("out"),
                                     "--populations", "arc_agrp"]), 1)

    def test_cli_writes_batches_and_manifest(self):
        with open(self.p("u.jsonl"), "w") as fh:
            for i in range(6):
                fh.write(json.dumps(unit("doi:%d" % i, "arc_agrp")) + "\n")
        with open(self.p("c.jsonl"), "w") as fh:
            for i in range(6):
                fh.write(json.dumps({"study_id": "doi:%d" % i,
                                     "full_text_available": True}) + "\n")
        out = self.p("out")
        self.assertEqual(pilot.main(["--units", self.p("u.jsonl"),
                                     "--candidates", self.p("c.jsonl"),
                                     "--out-dir", out, "--populations", "arc_agrp",
                                     "--papers", "4", "--batches", "2"]), 0)
        man = json.load(open(os.path.join(out, "manifest.json")))
        self.assertEqual(man["papers_selected"], 4)
        self.assertEqual(len(man["selected_study_ids"]), 4)
        self.assertEqual(man["extraction_method"], "llm")
        self.assertEqual(man["review_status"], "unreviewed")
        self.assertIn("sha256(seed|study_id)", man["selection_rule"])
        self.assertTrue(os.path.exists(os.path.join(out, "batch_01.jsonl")))
        self.assertTrue(os.path.exists(os.path.join(out, "batch_02.jsonl")))


if __name__ == "__main__":
    unittest.main()
