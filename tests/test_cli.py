"""The QC CLI and the synthetic fixtures it runs on."""

import _path  # noqa: F401
import contextlib
import io
import json
import os
import shutil
import tempfile
import unittest

from syn4rch.qc.run import main
from syn4rch.schema import read_jsonl

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
CLEAN = os.path.join(FIXTURES, "records_clean.jsonl")
BAD = os.path.join(FIXTURES, "records_bad.jsonl")
FOREIGN = os.path.join(FIXTURES, "records_foreign_version.jsonl")


def run_cli(argv):
    out = io.StringIO()
    err = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        status = main(argv)
    return status, out.getvalue(), err.getvalue()


class TestFixturesAreUnmistakablySynthetic(unittest.TestCase):
    def test_every_fixture_record_is_marked_synthetic(self):
        for path in (CLEAN, BAD):
            for record in read_jsonl(path):
                self.assertTrue(
                    (record.study_id or "").startswith("synthetic:"),
                    "%s has study_id=%r" % (path, record.study_id),
                )
                self.assertIsNone(record.pmid)

    def test_fixture_readme_warns_against_use_as_evidence(self):
        with open(os.path.join(FIXTURES, "README.md"), "r", encoding="utf-8") as handle:
            text = handle.read()
        self.assertIn("FABRICATED", text)
        self.assertIn("NEVER EVIDENCE", text)


class TestCliOnCleanFixture(unittest.TestCase):
    def test_exit_zero_and_readable_report(self):
        status, out, _err = run_cli([CLEAN])
        self.assertEqual(status, 0, out)
        self.assertIn("SYN4RCH QC report", out)
        self.assertIn("records read       : 4", out)
        self.assertIn("no hard schema violations", out)

    def test_conflicts_are_reported_without_failing_the_run(self):
        status, out, _err = run_cli([CLEAN])
        self.assertEqual(status, 0)
        self.assertIn("QC_CONFLICT_OPPOSED", out)
        self.assertIn("QC_POSSIBLE_NON_INDEPENDENCE", out)
        self.assertIn("nothing was resolved, averaged or dropped", out)

    def test_json_output_is_written_and_complete(self):
        tmp = tempfile.mkdtemp(prefix="syn4rch-cli-")
        try:
            target = os.path.join(tmp, "nested", "qc.json")
            status, out, _err = run_cli([CLEAN, "--json", target, "--quiet"])
            self.assertEqual(status, 0)
            self.assertEqual(out, "")
            with open(target, "r", encoding="utf-8") as handle:
                data = json.load(handle)
            self.assertEqual(data["n_records"], 4)
            self.assertTrue(data["ok"])
            self.assertIn("QC_CONFLICT_OPPOSED", data["counts_by_code"])
            for finding in data["findings"]:
                self.assertIn(finding["severity"], ("HARD", "WARNING", "INFO"))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_the_cli_does_not_modify_its_input(self):
        with open(CLEAN, "rb") as handle:
            before = handle.read()
        run_cli([CLEAN, "--quiet"])
        with open(CLEAN, "rb") as handle:
            self.assertEqual(handle.read(), before)


class TestCliOnBadFixture(unittest.TestCase):
    def test_exit_nonzero_on_hard_violations(self):
        status, out, _err = run_cli([BAD])
        self.assertEqual(status, 1)
        self.assertIn("hard violation", out)
        for code in (
            "A001_MISSING_SUPPORT_QUOTE",
            "A002_MISSING_SOURCE_LOCATOR",
            "A005_MISSING_OUTCOME_ID",
            "A006_EVIDENCE_CLASS_NOT_CAUSAL",
            "A007_MISSING_EXTRACTION_METHOD",
            "A004_EFFECT_SIZE_UNSUPPORTED",
            # schema v0.3.0: the four condition-evidence rules. One record of
            # the fixture violates each, so the fixture keeps carrying one
            # instance of every anti-fabrication defect the CLI can report.
            "A008_CONDITION_EVIDENCE_WITHOUT_LOCATOR",
            "A009_CONDITION_EVIDENCE_WITHOUT_QUOTE",
            "A010_CONDITION_EVIDENCE_MALFORMED",
            "A011_CONDITION_EVIDENCE_VALUE_MISMATCH",
        ):
            self.assertIn(code, out, out)

    def test_the_bad_fixture_carries_one_record_per_hard_rule(self):
        """The fixture's contract, asserted rather than left to the README.

        A rule added to ``syn4rch.schema.validate`` without a record in
        ``records_bad.jsonl`` is a rule the CLI has never been seen to report.
        Only the A family is required here: the S and V families are reached
        by records that could not be written to a fixture file at all (a
        tampered ``record_id``, a foreign ``schema_version``), and those are
        exercised in ``tests/test_validate.py``.
        """
        from syn4rch.schema import RULES, validate
        from syn4rch.schema.validate import hard_findings

        reported = set()
        for record in read_jsonl(BAD):
            reported |= {f.code for f in hard_findings(validate(record))}
        missing = sorted(
            code
            for code in RULES
            if code.startswith("A") and code not in reported
        )
        self.assertEqual(missing, [], "no record in records_bad.jsonl violates %s" % missing)

    def test_limit_truncates_sections(self):
        status, out, _err = run_cli([BAD, "--limit", "1"])
        self.assertEqual(status, 1)
        self.assertIn("more", out)


class TestCliErrorHandling(unittest.TestCase):
    def test_missing_file_is_a_usage_error(self):
        status, _out, err = run_cli([os.path.join(FIXTURES, "nope.jsonl")])
        self.assertEqual(status, 2)
        self.assertIn("no such file", err)

    def test_foreign_schema_version_fails_the_run(self):
        status, out, _err = run_cli([FOREIGN])
        self.assertEqual(status, 1)
        self.assertIn("LOAD ERRORS", out)
        self.assertIn("refusing to load", out)

    def test_any_version_flag_allows_inspection(self):
        status, out, _err = run_cli([FOREIGN, "--any-version"])
        self.assertIn("S006_SCHEMA_VERSION_MISMATCH", out)
        self.assertEqual(status, 1)

    def test_module_is_runnable_as_a_script(self):
        import subprocess
        import sys

        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        env = dict(os.environ, PYTHONPATH=os.path.join(root, "src"))
        proc = subprocess.run(
            [sys.executable, "-m", "syn4rch.qc.run", CLEAN],
            cwd=root,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr.decode())
        self.assertIn(b"SYN4RCH QC report", proc.stdout)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
