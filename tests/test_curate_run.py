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
from syn4rch.schema.io import dumps_record, read_jsonl, write_jsonl

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


class GateWillNotSilentlyReplaceADataset(unittest.TestCase):
    """The gate replaces --out. Right for a first corpus, fatal for a second pass."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.out = os.path.join(self.dir, "records.jsonl")
        self.report = os.path.join(self.dir, "rejects.md")
        self.proposals = os.path.join(self.dir, "proposals.jsonl")

    def _write_existing(self, records):
        with open(self.out, "w", encoding="utf-8") as fh:
            for r in records:
                fh.write(dumps_record(r) + "\n")

    def _write_proposals(self, records):
        with open(self.proposals, "w", encoding="utf-8") as fh:
            for r in records:
                data = r.to_dict()
                data.pop("record_id", None)
                fh.write(json.dumps(data) + "\n")

    def _run(self, *extra):
        return main([
            "gate", self.proposals, "--out", self.out, "--report", self.report,
            "--fulltext-dir", self.dir, *extra,
        ])

    def test_a_non_empty_output_is_refused_without_an_explicit_choice(self):
        existing = synthetic.good_record()
        self._write_existing([existing])
        self._write_proposals([synthetic.other_lab_record()])
        code = self._run()
        self.assertEqual(code, 1)
        with open(self.out, encoding="utf-8") as fh:
            kept = [l for l in fh if l.strip()]
        self.assertEqual(len(kept), 1, "the existing dataset must be untouched")

    def test_an_empty_output_needs_no_flag(self):
        open(self.out, "w").close()
        self._write_proposals([])
        self.assertEqual(self._run(), 0)

    def test_dry_run_never_touches_the_output(self):
        """Exit code says what the gate decided; the file must not move either way.

        With no retrieved text the proposal is rejected and the gate exits 1,
        which is its contract. What is at stake here is that the existing
        dataset survives a run that was only ever going to count.
        """
        self._write_existing([synthetic.good_record()])
        self._write_proposals([synthetic.other_lab_record()])
        self._run("--dry-run")
        with open(self.out, encoding="utf-8") as fh:
            self.assertEqual(len([l for l in fh if l.strip()]), 1)

    def test_append_merges_and_refuses_a_record_id_already_present(self):
        existing = synthetic.good_record()
        self._write_existing([existing])
        # the same observation re-proposed: same record_id
        self._write_proposals([synthetic.good_record()])
        code = self._run("--append")
        self.assertEqual(code, 1)
        with open(self.out, encoding="utf-8") as fh:
            self.assertEqual(len([l for l in fh if l.strip()]), 1)


class WithinBatchIdCollisionTest(unittest.TestCase):
    """D-060: two admitted records sharing a record_id must stop the gate.

    `background_genotype` is not a RECORD_ID_FIELDS component, so two records
    differing only in it hash to one id. The check against the existing file
    cannot see a batch that collides with itself.
    """

    def test_two_records_differing_only_in_genotype_share_an_id(self):
        from syn4rch.schema.record import TierARecord

        base = synthetic.good_record().to_dict()
        a = TierARecord.from_dict(
            {**base, "record_id": None, "background_genotype": "wild_type"}
        )
        b = TierARecord.from_dict(
            {
                **base,
                "record_id": None,
                "background_genotype": "npy:loss_of_function:global",
            }
        )
        self.assertEqual(
            a.record_id,
            b.record_id,
            "if this ever differs, genotype has entered the hash and D-060's "
            "no-id-map property is void",
        )
