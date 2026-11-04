"""Yield assessment: counts are read from units, bounds are never overstated."""

import _path  # noqa: F401
import json
import os
import tempfile
import unittest

from syn4rch.curate import yield_report as y


def unit(study, pop, method, outcome, section="Results"):
    return {"study_id": study, "candidate_population_id": pop,
            "candidate_method_id": method, "candidate_outcome_id": outcome,
            "section_label": section, "_kind": "extraction_unit_candidate",
            "_not_a_record": True}


class YieldTest(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def p(self, name):
        return os.path.join(self.d, name)

    def write(self, name, rows):
        with open(self.p(name), "w", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
        return self.p(name)

    def test_only_result_sections_count_toward_the_ceiling(self):
        units = [
            unit("doi:1", "arc_agrp", "chemo_hm3dq", "food_intake", "Results"),
            unit("doi:1", "arc_agrp", "chemo_hm3dq", "food_intake", "Discussion"),
            unit("doi:1", "arc_agrp", "chemo_hm3dq", "food_intake", "Introduction"),
            unit("doi:1", "arc_agrp", "chemo_hm3dq", "food_intake", "Abstract"),
        ]
        a = y.assess(units, {})
        self.assertEqual(a["n_units"], 4)
        self.assertEqual(a["n_result_units"], 1)
        self.assertEqual(a["n_abstract_units"], 1)
        # Four units describing one claim are one possible record, not four.
        self.assertEqual(a["ceiling_result"], 1)
        self.assertEqual(a["ceiling_all"], 1)

    def test_ceiling_counts_study_times_triple(self):
        units = [
            unit("doi:1", "arc_agrp", "chemo_hm3dq", "food_intake"),
            unit("doi:1", "arc_agrp", "chemo_hm3dq", "body_weight"),
            unit("doi:2", "arc_agrp", "chemo_hm3dq", "food_intake"),
        ]
        a = y.assess(units, {})
        self.assertEqual(a["ceiling_result"], 3)
        self.assertEqual(a["triples_result"], 2)
        self.assertEqual(a["papers_result_unit"], 2)

    def test_same_lab_repeats_do_not_pass_the_gate(self):
        one_lab = {"doi:%d" % i: ("smith", "same place") for i in range(1, 5)}
        units = [unit("doi:%d" % i, "arc_agrp", "chemo_hm3dq", "food_intake")
                 for i in range(1, 5)]
        a = y.assess(units, one_lab)
        cell = a["cells"][("arc_agrp", "food_intake")]
        self.assertEqual(len(cell["studies"]), 4)
        self.assertEqual(len(cell["labs"]), 1)
        self.assertEqual(a["gate_pass"], [])

    def test_three_studies_across_two_labs_passes(self):
        labs = {"doi:1": ("a", "x"), "doi:2": ("b", "y"), "doi:3": ("b", "y")}
        units = [unit("doi:%d" % i, "arc_agrp", "chemo_hm3dq", "food_intake")
                 for i in (1, 2, 3)]
        a = y.assess(units, labs)
        self.assertEqual(a["gate_pass"], [("arc_agrp", "food_intake")])

    def test_single_study_cell_is_counted_as_singleton(self):
        a = y.assess([unit("doi:1", "vmh_sf1", "opto_chr2", "blood_glucose")],
                     {"doi:1": ("a", "x")})
        self.assertEqual(a["singleton"], [("vmh_sf1", "blood_glucose")])
        self.assertEqual(a["replicated"], [])

    def test_missing_lab_identity_is_tracked_not_pooled(self):
        labs = {"doi:1": None, "doi:2": None}
        units = [unit("doi:1", "arc_agrp", "chemo_hm3dq", "food_intake"),
                 unit("doi:2", "arc_agrp", "chemo_hm3dq", "food_intake")]
        a = y.assess(units, labs)
        cell = a["cells"][("arc_agrp", "food_intake")]
        self.assertEqual(len(cell["labs"]), 0)
        self.assertEqual(len(cell["unlabelled_studies"]), 2)
        self.assertEqual(a["studies_without_lab"], 2)

    def test_report_states_bounds_and_renders(self):
        labs = {"doi:1": ("a", "x"), "doi:2": ("b", "y"), "doi:3": ("c", "z")}
        units = [unit("doi:%d" % i, "arc_agrp", "chemo_hm3dq", "food_intake")
                 for i in (1, 2, 3)]
        a = y.assess(units, labs)
        path = self.write("u.jsonl", units)
        text = y.render(a, units_path="u.jsonl", units_ts=y._utc_mtime(path), have_labs=True)
        self.assertIn("upper bound", text)
        self.assertIn("Co-occurrence is not a finding", text)
        self.assertIn("Does not establish", text)
        self.assertIn("`arc_agrp`", text)

    def test_missing_units_file_exits_nonzero(self):
        self.assertEqual(y.main(["--units", self.p("nope.jsonl"),
                                 "--candidates", self.p("nope2.jsonl"),
                                 "--out", self.p("o.md")]), 1)
        self.assertFalse(os.path.exists(self.p("o.md")))

    def test_empty_units_file_exits_nonzero(self):
        path = self.write("empty.jsonl", [])
        self.assertEqual(y.main(["--units", path, "--candidates", self.p("n.jsonl"),
                                 "--out", self.p("o.md")]), 1)

    def test_cli_writes_report_without_lab_index(self):
        path = self.write("u.jsonl", [unit("doi:1", "arc_agrp", "chemo_hm3dq", "food_intake")])
        out = self.p("o.md")
        self.assertEqual(y.main(["--units", path, "--candidates", self.p("absent.jsonl"),
                                 "--out", out]), 0)
        with open(out, encoding="utf-8") as fh:
            text = fh.read()
        self.assertIn("candidates.jsonl absent", text)


if __name__ == "__main__":
    unittest.main()
