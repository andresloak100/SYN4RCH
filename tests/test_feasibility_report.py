"""The feasibility report and its CLI, including the empty-dataset path."""

import _path  # noqa: F401
import contextlib
import io
import json
import os
import shutil
import tempfile
import unittest

import synthetic
from syn4rch.feasibility.report import (
    SMALL_CORPUS_RECORDS,
    build_report,
    render_markdown,
    top_line,
)
from syn4rch.feasibility.run import DEFAULT_RECORDS, DEFAULT_REPORT, main
from syn4rch.schema import write_jsonl

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
CLEAN = os.path.join(FIXTURES, "records_clean.jsonl")


def run_cli(argv):
    out = io.StringIO()
    err = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        status = main(argv)
    return status, out.getvalue(), err.getvalue()


def third_lab_record(**overrides):
    data = dict(
        study_id="synthetic:test-004",
        doi="synthetic:test-004",
        last_author="Third, D.",
        last_author_affiliation="Third Institute of Synthetic Data, Far Away",
        source_locator="Fig S4D (synthetic)",
    )
    data.update(overrides)
    return synthetic.good_record(**data)


class TestEmptyDataset(unittest.TestCase):
    def test_the_report_is_produced_and_says_so_in_one_sentence(self):
        report = build_report([], source_path="data/curated/records.jsonl")
        line = top_line(report)
        self.assertTrue(report.is_empty)
        self.assertIn("insufficient evidence", line)
        self.assertIn("0 records", line)
        self.assertTrue(line.startswith("**Verdict:"))
        self.assertNotIn("\n", line)  # one sentence, not a paragraph

    def test_a_missing_file_is_described_as_missing_not_as_empty(self):
        report = build_report(
            [], source_path="data/curated/records.jsonl", source_exists=False
        )
        self.assertIn("does not exist yet", top_line(report))

    def test_the_verdict_is_the_first_thing_in_the_markdown(self):
        text = render_markdown(build_report([]))
        lines = [l for l in text.splitlines() if l.strip()]
        self.assertTrue(lines[0].startswith("# SYN4RCH feasibility report"))
        self.assertTrue(lines[1].startswith("**Verdict:"))

    def test_no_degenerate_statistics_are_printed(self):
        text = render_markdown(build_report([]))
        self.assertIn("The five Phase 2 questions", text)
        self.assertIn("arithmetic on nothing", text)
        # A sparsity or an agreement rate over nothing is not printed at all.
        self.assertNotIn("0.0%", text)
        self.assertNotIn("| sparsity", text)
        self.assertNotIn("Mean agreement", text)

    def test_every_section_degrades_to_an_explicit_zero(self):
        report = build_report([])
        self.assertEqual(report.n_records_usable, 0)
        self.assertEqual(report.cells, [])
        self.assertEqual(report.verdicts, [])
        self.assertEqual(report.pairs, [])
        self.assertEqual(report.n_feasible, 0)
        self.assertIsNone(report.coverage.sparsity)

    def test_it_serialises(self):
        data = build_report([]).to_dict()
        self.assertTrue(data["is_empty"])
        self.assertIn("insufficient evidence", data["top_line"])
        json.dumps(data)  # must be JSON-serialisable


class TestSmallCorpusIsDescribedAsSmall(unittest.TestCase):
    def test_a_handful_of_records_does_not_get_a_confident_top_line(self):
        report = build_report(
            [synthetic.good_record(), synthetic.other_lab_record()],
            source_path="x.jsonl",
        )
        line = top_line(report)
        self.assertTrue(report.is_small_corpus)
        self.assertIn("insufficient evidence", line)
        self.assertIn("below the floor of", line)
        self.assertIn(str(SMALL_CORPUS_RECORDS), line)
        self.assertNotIn("\n", line)

    def test_the_tables_are_labelled_a_census(self):
        text = render_markdown(
            build_report([synthetic.good_record(), synthetic.other_lab_record()])
        )
        self.assertIn("census, not as statistics", text)

    def test_sparsity_is_not_printed_for_a_degenerate_matrix(self):
        text = render_markdown(build_report([synthetic.good_record()]))
        self.assertIn("not meaningful at 1 x 1", text)

    def test_explicit_zeros_appear_in_the_grid(self):
        report = build_report(
            [
                synthetic.good_record(
                    population_id="syn:pop_alpha", outcome_id="syn:outcome_1"
                ),
                synthetic.other_lab_record(
                    population_id="syn:pop_beta", outcome_id="syn:outcome_2"
                ),
            ]
        )
        text = render_markdown(report)
        self.assertEqual(report.coverage.n_cells_empty, 2)
        self.assertIn("Record counts per cell (zeros shown)", text)
        self.assertIn("cells with no records | 2", text)
        self.assertIn("| sparsity (fraction of cells empty) | 50.0% |", text)


class TestAFeasibleCorpusIsReportedAsSuch(unittest.TestCase):
    def setUp(self):
        # Enough records to clear the small-corpus floor: ten distinct
        # studies across ten labs, all agreeing, on one claim.
        self.records = [
            synthetic.good_record(
                study_id="synthetic:test-1%02d" % i,
                doi="synthetic:test-1%02d" % i,
                last_author="Author%02d, X." % i,
                last_author_affiliation="Institute %02d" % i,
                source_locator="Fig %dA (synthetic)" % i,
            )
            for i in range(SMALL_CORPUS_RECORDS + 2)
        ]

    def test_the_top_line_reports_the_feasible_count(self):
        report = build_report(self.records)
        self.assertFalse(report.is_small_corpus)
        self.assertEqual(report.n_feasible, 1)
        line = top_line(report)
        self.assertIn("1 of 1 (population, outcome) task(s) are feasible", line)

    def test_the_report_names_the_outcome(self):
        text = render_markdown(build_report(self.records))
        self.assertIn("Feasible for at least one population", text)
        self.assertIn("syn:outcome_1", text)


class TestTheReportShowsBeforeAndAfter(unittest.TestCase):
    """The v0.2.0 change has to be visible in the written report, both ways.

    A reader must be able to see what the stricter key bought (apparent
    contradictions that were heterogeneity) and what it cost (replication
    evidence given up). A report that printed only the new numbers would make
    the change unfalsifiable, which is the reason ``claim_key_pooled`` exists
    at all.
    """

    def _heterogeneous(self):
        # One pooled claim, two strict claims: a male cohort reporting an
        # increase against a female cohort reporting a decrease.
        return [
            synthetic.good_record(sex="male", direction="increase"),
            synthetic.other_lab_record(sex="female", direction="decrease"),
        ]

    def setUp(self):
        self.report = build_report(self._heterogeneous(), source_path="x.jsonl")
        self.text = render_markdown(self.report)

    def test_both_key_censuses_are_computed(self):
        comparison = self.report.key_comparison
        self.assertIsNotNone(comparison)
        self.assertEqual(comparison.pooled.n_cells, 1)
        self.assertEqual(comparison.strict.n_cells, 2)
        self.assertEqual(comparison.pooled.n_cells_directly_opposed, 1)
        self.assertEqual(comparison.strict.n_cells_directly_opposed, 0)

    def test_the_markdown_prints_both_columns(self):
        self.assertIn("pooled key (v0.1.0)", self.text)
        self.assertIn("strict key (v0.2.0)", self.text)
        self.assertIn("What it bought", self.text)
        self.assertIn("What it cost", self.text)

    def test_the_cost_is_stated_and_not_netted_out(self):
        self.assertEqual(self.report.key_comparison.n_multi_study_cells_lost, 1)
        self.assertIn("multi-study", self.text)

    def test_the_usability_change_is_attributed_to_its_class(self):
        report = build_report(
            [synthetic.good_record(), synthetic.molecular_record()]
        )
        comparison = report.usability_comparison
        self.assertEqual(comparison.n_usable_v1, 1)
        self.assertEqual(comparison.n_usable_v2, 2)
        self.assertEqual(
            comparison.admitted_by_v2_only_by_class, {"molecular": 1}
        )
        text = render_markdown(report)
        self.assertIn("The usability predicate also changed", text)
        self.assertIn("corrected predicate", text)

    def test_the_pooled_verdicts_are_labelled_as_not_a_result(self):
        self.assertIn("not a result", self.text)

    def test_every_agreement_table_carries_the_condition_specificity(self):
        report = build_report(
            [synthetic.good_record(), synthetic.other_lab_record()]
        )
        text = render_markdown(report)
        self.assertIn("min cond", text)
        self.assertIn("minimum condition specificity", text)
        self.assertIn("both unknown", text)

    def test_the_report_says_no_claim_specifies_all_four_conditions(self):
        report = build_report(
            [
                synthetic.good_record(challenge=None),
                synthetic.other_lab_record(challenge=None),
            ]
        )
        text = render_markdown(report)
        self.assertIn("No claim in this corpus reaches specificity 5 of 5", text)

    def test_the_splits_section_says_the_key_change_does_not_reach_it(self):
        self.assertIn("claim key change does not reach this table", self.text)

    def test_the_report_points_at_the_identity_analysis(self):
        self.assertIn("perturbation_identity.md", self.text)

    def test_it_serialises_with_both_comparisons(self):
        data = self.report.to_dict()
        self.assertIsNotNone(data["key_comparison"])
        self.assertIsNotNone(data["usability_comparison"])
        self.assertIn("pooled", data["key_comparison"])
        self.assertIn("strict", data["key_comparison"])
        json.dumps(data)


class TestCli(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="syn4rch-feasibility-")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def path(self, name):
        return os.path.join(self.tmp, name)

    def test_a_missing_input_is_not_an_error(self):
        out_path = self.path("nested/feasibility_report.md")
        status, out, err = run_cli(
            ["--records", self.path("absent.jsonl"), "--out", out_path]
        )
        self.assertEqual(status, 0, err)
        self.assertTrue(os.path.isfile(out_path))
        with open(out_path, "r", encoding="utf-8") as handle:
            text = handle.read()
        self.assertIn("does not exist yet", text)
        self.assertIn("insufficient evidence", out)

    def test_an_empty_input_is_not_an_error(self):
        records_path = self.path("records.jsonl")
        open(records_path, "w").close()
        out_path = self.path("report.md")
        status, out, err = run_cli(
            ["--records", records_path, "--out", out_path]
        )
        self.assertEqual(status, 0, err)
        with open(out_path, "r", encoding="utf-8") as handle:
            self.assertIn("is empty (0 records)", handle.read())

    def test_insufficient_evidence_exits_zero(self):
        status, _out, _err = run_cli(
            ["--records", CLEAN, "--out", self.path("report.md")]
        )
        self.assertEqual(status, 0)

    def test_json_output(self):
        status, _out, err = run_cli(
            [
                "--records",
                CLEAN,
                "--out",
                self.path("report.md"),
                "--json",
                self.path("feasibility.json"),
                "--quiet",
            ]
        )
        self.assertEqual(status, 0, err)
        with open(self.path("feasibility.json"), "r", encoding="utf-8") as handle:
            data = json.load(handle)
        self.assertEqual(data["coverage"]["n_records_read"], 4)
        self.assertIn("thresholds", data["verdict_summary"])
        self.assertIn("by_lab", data["fold_sets"])

    def test_stdout_mode_writes_no_file(self):
        out_path = self.path("report.md")
        status, out, _err = run_cli(
            ["--records", CLEAN, "--out", out_path, "--stdout"]
        )
        self.assertEqual(status, 0)
        self.assertFalse(os.path.exists(out_path))
        self.assertIn("# SYN4RCH feasibility report", out)

    def test_quiet_suppresses_stdout_but_still_writes_the_report(self):
        out_path = self.path("report.md")
        status, out, _err = run_cli(
            ["--records", CLEAN, "--out", out_path, "--quiet"]
        )
        self.assertEqual(status, 0)
        self.assertEqual(out, "")
        self.assertTrue(os.path.isfile(out_path))

    def test_an_unreadable_line_exits_one_and_is_reported(self):
        records_path = self.path("broken.jsonl")
        with open(records_path, "w", encoding="utf-8") as handle:
            handle.write("{not json\n")
        out_path = self.path("report.md")
        status, _out, err = run_cli(
            ["--records", records_path, "--out", out_path, "--quiet"]
        )
        self.assertEqual(status, 1)
        self.assertIn("load error", err)
        with open(out_path, "r", encoding="utf-8") as handle:
            text = handle.read()
        self.assertIn("could not be read", text)

    def test_it_does_not_modify_its_input(self):
        records_path = self.path("records.jsonl")
        write_jsonl(records_path, [synthetic.good_record(), synthetic.other_lab_record()])
        with open(records_path, "rb") as handle:
            before = handle.read()
        run_cli(["--records", records_path, "--out", self.path("r.md"), "--quiet"])
        with open(records_path, "rb") as handle:
            self.assertEqual(handle.read(), before)

    def test_it_writes_only_where_it_was_told_to(self):
        out_path = self.path("only/report.md")
        run_cli(
            ["--records", CLEAN, "--out", out_path, "--quiet"]
        )
        written = []
        for root, _dirs, files in os.walk(self.tmp):
            for name in files:
                written.append(os.path.join(root, name))
        self.assertEqual(written, [out_path])

    def test_the_default_paths_are_the_documented_ones(self):
        self.assertEqual(DEFAULT_RECORDS, os.path.join("data", "curated", "records.jsonl"))
        self.assertEqual(DEFAULT_REPORT, os.path.join("reports", "feasibility_report.md"))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
