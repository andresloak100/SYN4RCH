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


if __name__ == "__main__":
    unittest.main()
