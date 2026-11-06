"""The three D-056 generalisation probes.

The properties asserted here are the ones that would make the probes lie rather
than merely be wrong: a depth that disagrees with the baseline's own backoff, a
cell key copied instead of read from the baseline, a permutation that changes a
study's label distribution, and a `depth_0` stratum where the instrument is
supposed to be identical to the floor by construction.
"""

import _path  # noqa: F401
import collections
import io
import json
import os
import random
import shutil
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

import synthetic
from syn4rch.models.baselines import PERTURBATION_AXIS_CELL_BACKOFF
from syn4rch.models.generalise import (
    CELL_FIELDS,
    CIRCUIT,
    DEPTH_STRATA,
    HEADLINE_SLICE,
    MAX_NULL_ESTABLISHED_RATE,
    PermutationNull,
    _permute_within,
    build_parser,
    cross_group_probe,
    depth_probe,
    main,
    permutation_probe,
    render_markdown,
    run_generalisation,
    run_paired,
)
from syn4rch.models.task import build_dataset
from syn4rch.schema.io import write_jsonl

OUTCOMES = ("food_intake", "body_weight", "energy_expenditure")
DIRECTIONS = ("increase", "decrease", "no_change")


def corpus(n_studies=12, n_outcomes=3, singletons=True):
    """A corpus with recurring cells, several studies per lab, and singletons.

    Study *k* belongs to ``Author k // 2``, so every lab has two studies and
    the lab grouping is strictly coarser than the study grouping -- the
    property the real corpus has and that the committed code did not express
    until D-056.

    With ``singletons``, each study also contributes one row in a cell of its
    own, so that the `depth_0` stratum is populated: holding that study out
    empties its cell, and the candidate must then fall back to the global mode
    exactly as the floor does.
    """
    rows = []
    for study in range(n_studies):
        for index in range(n_outcomes):
            outcome = OUTCOMES[index % len(OUTCOMES)]
            rows.append(
                synthetic.good_record(
                    study_id="synthetic:study-%02d" % study,
                    doi="synthetic:study-%02d" % study,
                    last_author="Author %d" % (study // 2),
                    last_author_affiliation="Place %d" % study,
                    population_raw="population %d" % (index % 2),
                    population_id="syn:pop_%d" % (index % 2),
                    outcome_id=outcome,
                    outcome_raw=outcome,
                    direction=DIRECTIONS[(study + index) % 3],
                    source_locator="Fig %d%s (synthetic)" % (study, "ABC"[index]),
                    support_quote=(
                        "Synthetic quote: the readout changed in study %d." % study
                    ),
                )
            )
        if singletons:
            rows.append(
                synthetic.good_record(
                    study_id="synthetic:study-%02d" % study,
                    doi="synthetic:study-%02d" % study,
                    last_author="Author %d" % (study // 2),
                    last_author_affiliation="Place %d" % study,
                    population_raw="lone population %d" % study,
                    population_id="syn:lone_%d" % study,
                    outcome_id="body_weight",
                    outcome_raw="body_weight",
                    direction=DIRECTIONS[study % 3],
                    source_locator="Fig %dZ (synthetic)" % study,
                    support_quote=(
                        "Synthetic quote: the lone readout changed in study %d."
                        % study
                    ),
                )
            )
    return rows


def dataset(**kwargs):
    return build_dataset(corpus(**kwargs))


class CellFieldsTest(unittest.TestCase):
    def test_the_cell_is_read_from_the_baselines_own_chain(self):
        """A copied cell key would measure the depth of a cell nothing uses."""
        self.assertEqual(CELL_FIELDS, PERTURBATION_AXIS_CELL_BACKOFF[0][1])


class RunPairedTest(unittest.TestCase):
    def test_every_pooled_row_is_tested_once(self):
        ds = dataset()
        run = run_paired(ds, "study_id")
        self.assertEqual(run.n_rows, len(ds))

    def test_depth_and_the_backoff_level_never_disagree(self):
        """Computed two independent ways; a disagreement is a defect.

        `primary_cell_seen` comes from the fitted baseline's own backoff
        accounting, `cell_depth` from counting groups in the training fold. They
        must say the same thing about whether the cell was populated.
        """
        for grouping in ("study_id", "lab"):
            run = run_paired(dataset(), grouping)
            self.assertEqual(run.depth_disagreements(), 0, grouping)

    def test_depth_zero_is_exactly_the_complement_of_the_headline_slice(self):
        run = run_paired(dataset(), "lab")
        circuit = [i for i in run.paired if run.perturbation_class[i] == CIRCUIT]
        seen = set(run.headline_rows())
        self.assertTrue(seen)
        for i in circuit:
            self.assertEqual(run.cell_depth[i] > 0, i in seen)

    def test_the_lab_grouping_is_coarser_than_the_study_grouping(self):
        ds = dataset()
        self.assertLess(
            ds.grouping_coverage("lab")["n_groups"],
            ds.grouping_coverage("study_id")["n_groups"],
        )

    def test_an_undeclared_grouping_raises_rather_than_defaulting(self):
        with self.assertRaises(KeyError):
            run_paired(dataset(), "by_vibes")

    def test_supplied_folds_are_used_rather_than_regenerated(self):
        ds = dataset()
        run = run_paired(ds, "study_id", folds=[(tuple(range(10)), (10, 11))])
        self.assertEqual(run.n_folds, 1)
        self.assertEqual(run.n_rows, 2)


class CrossGroupProbeTest(unittest.TestCase):
    def test_both_groupings_are_reported_not_only_the_harder_one(self):
        got = cross_group_probe(dataset(), n_bootstrap=50)
        self.assertEqual([grouping for grouping, _, _ in got], ["study_id", "lab"])
        for grouping, _, item in got:
            self.assertEqual(item.grouping, grouping)
            self.assertEqual(item.name, HEADLINE_SLICE)

    def test_the_slice_carries_its_grouping_into_its_dict(self):
        _, _, item = cross_group_probe(
            dataset(), groupings=("lab",), n_bootstrap=50
        )[0]
        self.assertEqual(item.to_dict()["grouping"], "lab")


class DepthProbeTest(unittest.TestCase):
    def test_the_strata_partition_the_circuit_rows(self):
        run = run_paired(dataset(), "lab")
        strata = depth_probe(run, n_bootstrap=50)
        circuit = [i for i in run.paired if run.perturbation_class[i] == CIRCUIT]
        self.assertEqual(sum(item.n_rows for item in strata), len(circuit))
        self.assertEqual([item.name for item in strata], [n for n, _, _ in DEPTH_STRATA])

    def test_depth_zero_has_no_difference_by_construction(self):
        """Where the cell was never populated both baselines are the global mode.

        Not a near-zero tolerance: the two predictions are the same value on
        every such row, so the difference is exactly 0.0. A tolerance here
        would hide the one bug this assertion can catch.
        """
        run = run_paired(dataset(), "lab")
        zero = depth_probe(run, n_bootstrap=50)[0]
        self.assertEqual(zero.name, "depth_0")
        if zero.n_rows:
            self.assertEqual(zero.delta, 0.0)
            self.assertEqual(zero.n_only_baseline_right, 0)
            self.assertEqual(zero.n_only_reference_right, 0)

    def test_depth_zero_is_populated_by_a_cell_only_its_own_study_reached(self):
        run = run_paired(dataset(), "lab")
        zero = depth_probe(run, n_bootstrap=50)[0]
        self.assertGreater(zero.n_rows, 0)

    def test_an_empty_stratum_is_reported_as_empty_not_omitted(self):
        run = run_paired(
            dataset(n_studies=2, n_outcomes=1, singletons=False), "lab"
        )
        strata = {item.name: item for item in depth_probe(run, n_bootstrap=50)}
        self.assertEqual(set(strata), {n for n, _, _ in DEPTH_STRATA})
        for item in strata.values():
            if item.n_rows == 0:
                self.assertIsNone(item.delta)
                self.assertEqual(item.verdict, "undefined")
                self.assertTrue(item.caveats)


class PermuteWithinTest(unittest.TestCase):
    def test_each_group_keeps_its_own_label_multiset(self):
        ds = dataset()
        keys = ds.group_values("study_id")
        before = collections.defaultdict(collections.Counter)
        for key, example in zip(keys, ds.examples):
            before[key][example.label] += 1
        shuffled = _permute_within(ds, keys, random.Random(7))
        after = collections.defaultdict(collections.Counter)
        for key, example in zip(keys, shuffled.examples):
            after[key][example.label] += 1
        self.assertEqual(dict(before), dict(after))

    def test_the_features_and_groups_are_untouched(self):
        ds = dataset()
        shuffled = _permute_within(
            ds, ds.group_values("study_id"), random.Random(7)
        )
        for original, permuted in zip(ds.examples, shuffled.examples):
            self.assertEqual(original.features, permuted.features)
            self.assertEqual(dict(original.group_keys), dict(permuted.group_keys))
            self.assertEqual(original.record_id, permuted.record_id)

    def test_the_original_dataset_is_not_mutated(self):
        ds = dataset()
        before = [e.label for e in ds.examples]
        _permute_within(ds, ds.group_values("study_id"), random.Random(7))
        self.assertEqual([e.label for e in ds.examples], before)

    def test_the_records_stay_aligned_so_folds_do_not_move(self):
        ds = dataset()
        shuffled = _permute_within(
            ds, ds.group_values("study_id"), random.Random(7)
        )
        self.assertEqual(len(shuffled.records), len(ds.records))
        self.assertIs(shuffled.records[0], ds.records[0])


class PermutationProbeTest(unittest.TestCase):
    def test_it_reports_a_verdict_for_every_permutation(self):
        null = permutation_probe(
            dataset(), n_permutations=6, n_bootstrap=50, seed=11
        )
        self.assertEqual(null.n_permutations, 6)
        self.assertEqual(sum(null.verdicts.values()), 6)
        self.assertEqual(null.slice_name, HEADLINE_SLICE)

    def test_the_empirical_p_adds_one_so_it_is_never_zero(self):
        null = permutation_probe(
            dataset(), n_permutations=4, n_bootstrap=50, seed=11
        )
        self.assertIsNotNone(null.empirical_p)
        self.assertGreater(null.empirical_p, 0.0)
        self.assertAlmostEqual(
            null.empirical_p,
            (1 + null.n_deltas_at_or_above_observed) / 5.0,
        )

    def test_the_rate_is_compared_against_the_declared_threshold(self):
        null = PermutationNull(
            grouping="study_id",
            slice_name=HEADLINE_SLICE,
            n_permutations=100,
            verdicts={"established": 2, "no_difference": 98},
        )
        self.assertAlmostEqual(null.established_rate, 0.02)
        self.assertTrue(null.calibrated)
        liberal = PermutationNull(
            grouping="study_id",
            slice_name=HEADLINE_SLICE,
            n_permutations=100,
            verdicts={"established": 20, "no_difference": 80},
        )
        self.assertFalse(liberal.calibrated)
        self.assertLess(MAX_NULL_ESTABLISHED_RATE, 0.2)

    def test_no_permutations_means_no_rate_rather_than_a_rate_of_zero(self):
        null = PermutationNull(
            grouping="study_id", slice_name=HEADLINE_SLICE, n_permutations=0
        )
        self.assertIsNone(null.established_rate)
        self.assertIsNone(null.calibrated)
        self.assertIsNone(null.empirical_p)
        self.assertIsNone(null.null_centre)

    def test_the_null_centre_is_the_median_not_an_assumed_zero(self):
        null = PermutationNull(
            grouping="study_id",
            slice_name=HEADLINE_SLICE,
            n_permutations=3,
            deltas=(0.1, 0.2, 0.9),
        )
        self.assertAlmostEqual(null.null_centre, 0.2)


class RunAndRenderTest(unittest.TestCase):
    def test_the_result_carries_all_three_probes(self):
        result = run_generalisation(
            dataset(), n_bootstrap=50, n_permutations=3, seed=5
        )
        data = result.to_dict()
        self.assertEqual(
            [item["grouping"] for item in data["cross_group"]],
            ["study_id", "lab"],
        )
        self.assertEqual(data["depth"]["grouping"], "lab")
        self.assertEqual(
            data["depth"]["n_rows_where_backoff_and_training_set_disagree"], 0
        )
        self.assertEqual(data["permutation"]["n_permutations"], 3)
        self.assertEqual(data["cell_fields"], list(CELL_FIELDS))

    def test_permutations_zero_skips_probe_three_rather_than_faking_it(self):
        result = run_generalisation(
            dataset(), n_bootstrap=50, n_permutations=0
        )
        self.assertIsNone(result.permutation)
        self.assertIn("Not run", render_markdown(result))

    def test_the_report_names_every_probe_and_its_caveats(self):
        text = render_markdown(
            run_generalisation(dataset(), n_bootstrap=50, n_permutations=3, seed=5)
        )
        for needle in (
            "Probe 1",
            "Probe 2",
            "Probe 3",
            "does **not** reopen D-048",
            "What this report does not claim",
            "leave-one-population-out",
        ):
            self.assertIn(needle, text)

    def test_the_report_states_that_the_null_is_not_assumed_centred(self):
        text = render_markdown(
            run_generalisation(dataset(), n_bootstrap=50, n_permutations=3, seed=5)
        )
        self.assertIn("null is not centred on zero", text)


class CLITest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.records = os.path.join(self.tmp, "records.jsonl")
        write_jsonl(self.records, corpus())

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = main(argv)
        return code, out.getvalue(), err.getvalue()

    def test_it_writes_markdown_and_json(self):
        md = os.path.join(self.tmp, "gen.md")
        js = os.path.join(self.tmp, "gen.json")
        code, out, err = self._run(
            [
                "--records", self.records,
                "--out", md,
                "--json", js,
                "--bootstrap", "50",
                "--permutations", "2",
            ]
        )
        self.assertEqual(code, 0, err)
        self.assertIn("Probe 1", open(md, encoding="utf-8").read())
        data = json.load(open(js, encoding="utf-8"))
        self.assertEqual(len(data["cross_group"]), 2)
        self.assertIn("study_id", out)

    def test_an_undeclared_depth_grouping_is_refused(self):
        code, _, err = self._run(
            ["--records", self.records, "--depth-grouping", "by_vibes"]
        )
        self.assertEqual(code, 2)
        self.assertIn("undeclared grouping", err)

    def test_bad_numeric_arguments_are_refused(self):
        for argv in (["--folds", "1"], ["--bootstrap", "0"], ["--permutations", "-1"]):
            code, _, err = self._run(["--records", self.records] + argv)
            self.assertEqual(code, 2, argv)
            self.assertTrue(err.startswith("error:"), err)

    def test_an_empty_corpus_is_refused_not_reported_as_a_result(self):
        empty = os.path.join(self.tmp, "empty.jsonl")
        open(empty, "w", encoding="utf-8").close()
        code, _, err = self._run(["--records", empty])
        self.assertEqual(code, 1)
        self.assertIn("nothing to probe", err)

    def test_the_parser_declares_the_pre_registered_defaults(self):
        args = build_parser().parse_args([])
        self.assertEqual(args.permutations, 200)
        self.assertEqual(args.depth_grouping, "lab")
        self.assertEqual(args.bootstrap, 1000)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
