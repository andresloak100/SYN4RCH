"""CLI: both subcommands against temp paths only (never data/ or reports/)."""

import _path  # noqa: F401
import io
import json
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout, redirect_stderr

import synthetic
from syn4rch.curate.run import main
from syn4rch.schema.io import read_jsonl, write_jsonl

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "jats")
DOI = "10.0000/synthetic.0001"


class TestCli(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.ft = os.path.join(self.tmp, "fulltext")
        os.makedirs(self.ft)
        for name, pmcid in (("synthetic_complete.xml", "PMC0000001"),
                            ("synthetic_truncated.xml", "PMC0000003"),
                            ("synthetic_malformed.xml", "PMC0000004")):
            shutil.copy(os.path.join(FIX, name), os.path.join(self.ft, pmcid + ".xml"))

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def run_cli(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = main(argv)
        return code, out.getvalue()

    def test_units(self):
        cands = os.path.join(self.tmp, "candidates.jsonl")
        with open(cands, "w") as handle:
            handle.write(json.dumps({"pmcid": "PMC0000001", "study_id": "doi:" + DOI}) + "\n")
            handle.write("{not json\n")
        out_path = os.path.join(self.tmp, "units.jsonl")
        code, out = self.run_cli(["units", "--candidates", cands,
                                  "--fulltext-dir", self.ft, "--out", out_path])
        self.assertEqual(code, 0)
        rows = [json.loads(l) for l in open(out_path)]
        self.assertTrue(rows)
        self.assertTrue(all(r["_not_a_record"] for r in rows))
        self.assertIn("fulltext_files_found                             3", out)
        self.assertIn("candidate_lines_unreadable                       1", out)
        self.assertIn("units_written                                    %d" % len(rows), out)

    def test_gate(self):
        good = synthetic.good_record(
            study_id="doi:" + DOI, doi=DOI, source_locator="Results para 1",
            support_quote="Chemogenetic activation of Agrp neurons increased food intake (n = 8 mice, p < 0.01)")
        bad = synthetic.good_record(
            study_id="doi:" + DOI, doi=DOI, source_locator="Results para 1",
            support_quote="chemogenetic activation of Agrp neurons increased food intake (n = 8 mice, p < 0.01)")
        proposals = os.path.join(self.tmp, "proposals.jsonl")
        write_jsonl(proposals, [good, bad])
        with open(proposals, "a") as handle:
            handle.write("{broken\n")
        out_path = os.path.join(self.tmp, "curated", "records.jsonl")
        report = os.path.join(self.tmp, "rejects.md")
        code, out = self.run_cli(["gate", proposals, "--fulltext-dir", self.ft,
                                  "--out", out_path, "--report", report])
        self.assertEqual(code, 1)
        admitted = read_jsonl(out_path)
        self.assertEqual([r.record_id for r in admitted], [good.record_id])
        text = open(report).read()
        self.assertIn("C004_QUOTE_NOT_VERBATIM", text)
        self.assertIn("| rejected | 1 |", text)
        self.assertIn("| proposal lines unreadable | 1 |", text)

    def test_gate_missing_file_is_usage_error(self):
        code, _ = self.run_cli(["gate", os.path.join(self.tmp, "nope.jsonl"),
                                "--fulltext-dir", self.ft, "--dry-run"])
        self.assertEqual(code, 2)



class FillDerivedTest(unittest.TestCase):
    """The derived perturbation fields are filled at admission, not judged."""

    def rec(self, **kw):
        from syn4rch.schema import TierARecord

        base = dict(
            study_id="synthetic:fd-1",
            population_raw="synthetic population",
            outcome_raw="synthetic readout",
            method="conditional_knockout",
            valence="unclear",
        )
        base.update(kw)
        return TierARecord(**base)

    def test_molecular_knockout_gets_a_signed_axis(self):
        from syn4rch.curate.run import fill_derived

        record = fill_derived(self.rec())
        self.assertEqual(record.perturbation_class.value, "molecular")
        self.assertEqual(record.perturbation_polarity.value, "decrease")
        # valence is untouched: a knockout's effect on firing stays unknown.
        self.assertEqual(record.valence.value, "unclear")

    def test_overexpression_is_the_opposite_sign(self):
        from syn4rch.curate.run import fill_derived

        record = fill_derived(self.rec(method="conditional_overexpression"))
        self.assertEqual(record.perturbation_polarity.value, "increase")

    def test_an_explicit_value_is_never_overwritten(self):
        from syn4rch.curate.run import fill_derived

        record = fill_derived(
            self.rec(method="local_pharmacology", valence=None,
                     perturbation_polarity="decrease")
        )
        self.assertEqual(record.perturbation_polarity.value, "decrease")

    def test_filling_does_not_change_record_id(self):
        """The admission path must never silently re-identify a record."""
        from syn4rch.curate.run import fill_derived

        record = self.rec()
        before = record.record_id
        fill_derived(record)
        self.assertEqual(record.record_id, before)
        self.assertTrue(record.id_is_consistent())

    def test_no_method_means_no_derivation(self):
        from syn4rch.curate.run import fill_derived

        record = fill_derived(self.rec(method=None, valence=None))
        self.assertIsNone(record.perturbation_class)
        self.assertIsNone(record.perturbation_polarity)

    def test_a_filled_record_passes_the_class_check(self):
        from syn4rch.curate.run import fill_derived
        from syn4rch.schema.validate import validate

        from syn4rch.schema.enums import PerturbationMethod

        for method in PerturbationMethod:
            record = fill_derived(self.rec(method=method.value))
            codes = [f.code for f in validate(record)]
            self.assertNotIn("S007_PERTURBATION_CLASS_MISMATCH", codes, method)
            self.assertNotIn(
                "W013_PERTURBATION_POLARITY_RULE_MISMATCH", codes, method
            )


if __name__ == "__main__":
    unittest.main()
