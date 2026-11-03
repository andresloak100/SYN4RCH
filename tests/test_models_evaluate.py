"""Metrics, the grouping guard, the group bootstrap, calibration, folds, CLI."""

import _path  # noqa: F401
import os
import sys
import tempfile
import unittest

import synthetic
from syn4rch.models import metrics as M
from syn4rch.models.baselines import Baseline, GlobalMajority, PopulationOutcomeMajority
from syn4rch.models.evaluate import (
    BUILTIN_SOURCE,
    Estimate,
    FoldProvider,
    MetricBundle,
    UnlabelledMetricError,
    builtin_grouped_folds,
    cross_validate,
    make_estimate,
    resolve_fold_provider,
    run_benchmark,
)
from syn4rch.models.report import render_markdown
from syn4rch.models.task import LABELS, build_dataset


def corpus(n=40, n_studies=10, n_labs=7):
    """FABRICATED records spread over synthetic studies, labs and cells."""
    return [
        synthetic.good_record(
            study_id="synthetic:s%d" % (i % n_studies),
            doi="synthetic:s%d" % (i % n_studies),
            last_author="Lab%d" % (i % n_labs),
            last_author_affiliation="Inst%d" % (i % n_labs),
            population_id=["arc_agrp", "vmh_sf1", "arc_pomc", "lha_orexin"][i % 4],
            outcome_id=["food_intake", "body_weight", "core_temperature"][i % 3],
            valence=["activation", "inhibition"][i % 2],
            direction=["increase", "increase", "decrease", "no_change"][i % 4],
            source_locator="Fig %d" % i,
        )
        for i in range(n)
    ]


class TestUnlabelledMetricGuard(unittest.TestCase):
    def test_an_estimate_without_a_grouping_cannot_be_constructed(self):
        for bad in (None, "", "by_vibes", 3):
            with self.assertRaises(UnlabelledMetricError, msg=repr(bad)):
                Estimate(name="acc", grouping=bad, value=0.5, n_rows=10, n_groups=3)

    def test_a_bundle_without_a_grouping_cannot_be_constructed(self):
        with self.assertRaises(UnlabelledMetricError):
            MetricBundle(baseline="b", grouping=None)
        with self.assertRaises(UnlabelledMetricError):
            MetricBundle(baseline="b", grouping="undeclared")

    def test_a_mislabelled_estimate_cannot_be_filed_in_a_bundle(self):
        bundle = MetricBundle(baseline="b", grouping="study_id")
        estimate = make_estimate("acc", "lab", 0.5, 30, 6)
        with self.assertRaises(UnlabelledMetricError):
            bundle.add(estimate)

    def test_every_emitted_metric_dict_carries_its_grouping(self):
        bundles = run_benchmark(build_dataset(corpus()), n_bootstrap=20).bundles
        self.assertTrue(bundles)
        for bundle in bundles:
            d = bundle.to_dict()
            self.assertIn(d["grouping"], ("study_id", "lab", "population", "population_outcome", "record"))
            for metric in d["metrics"].values():
                self.assertEqual(metric["grouping"], d["grouping"])
            if d["calibration"]:
                self.assertEqual(d["calibration"]["grouping"], d["grouping"])

    def test_metric_functions_return_bare_numbers_not_dicts(self):
        y = ["increase", "decrease"]
        self.assertIsInstance(M.accuracy(y, y), float)
        self.assertIsInstance(M.macro_f1(y, y), float)

    def test_cross_validate_refuses_an_undeclared_grouping(self):
        with self.assertRaises(UnlabelledMetricError):
            cross_validate(GlobalMajority, build_dataset(corpus()), "by_vibes")

    def test_small_n_numbers_carry_a_caveat(self):
        e = make_estimate("acc", "study_id", 0.8, 5, 2)
        self.assertFalse(e.meaningful)
        self.assertIn("does not carry meaning", e.render())
        self.assertIn("n=5 rows / 2 study_id groups", e.render())
        self.assertTrue(make_estimate("acc", "study_id", 0.8, 50, 10).meaningful)


class TestMetrics(unittest.TestCase):
    def test_accuracy_ignores_abstentions_and_is_none_when_all_abstain(self):
        self.assertEqual(M.accuracy(["increase", "decrease"], ["increase", None]), 1.0)
        self.assertIsNone(M.accuracy(["increase"], [None]))

    def test_macro_f1_known_value(self):
        y = ["increase"] * 6 + ["decrease"] * 4
        p = ["increase"] * 5 + ["decrease"] * 5
        # increase: P=1, R=5/6 -> 10/11 ; decrease: P=0.8, R=1 -> 8/9
        self.assertAlmostEqual(M.macro_f1(y, p), (10 / 11 + 8 / 9) / 2)

    def test_macro_f1_skips_absent_classes_and_penalises_unpredicted_ones(self):
        self.assertAlmostEqual(M.macro_f1(["increase", "increase"], ["increase", "increase"]), 1.0)
        self.assertAlmostEqual(M.macro_f1(["increase", "decrease"], ["increase", "increase"]), (2 / 3) / 2)

    def test_auc_is_undefined_for_constant_scores(self):
        p = [{"increase": 0.6, "decrease": 0.3, "no_change": 0.1}] * 4
        y = ["increase", "decrease", "increase", "no_change"]
        self.assertIsNone(M.roc_auc_ovr_macro(y, p)[0])

    def test_auc_perfect_ranking(self):
        y = ["increase", "increase", "decrease", "decrease"]
        p = [{"increase": s, "decrease": 1 - s, "no_change": 0.0} for s in (0.9, 0.8, 0.2, 0.1)]
        self.assertEqual(M.roc_auc_ovr(y, p, "increase"), 1.0)


class _Overconfident(Baseline):
    """Says 'increase' with probability 0.95, every time."""

    name = "overconfident_stub"
    emits_proba = True

    def _fit(self, examples):
        pass

    def _predict_one(self, features):
        return "increase"

    def predict_proba(self, features):
        return [{"increase": 0.95, "decrease": 0.03, "no_change": 0.02} for _ in features]


class TestCalibration(unittest.TestCase):
    def test_an_overconfident_predictor_is_shown_to_be_miscalibrated(self):
        y = ["increase", "decrease"] * 20  # right exactly half the time
        p = [{"increase": 0.95, "decrease": 0.03, "no_change": 0.02}] * 40
        bins = M.reliability_bins(y, p, n_bins=10, min_n=10)
        top = [b for b in bins if b.n][0]
        self.assertEqual(top.n, 40)
        self.assertAlmostEqual(top.mean_confidence, 0.95)
        self.assertAlmostEqual(top.observed_accuracy, 0.5)
        self.assertAlmostEqual(top.gap, -0.45)
        self.assertAlmostEqual(M.expected_calibration_error(bins), 0.45)
        # Brier: right rows 0.05^2+.03^2+.02^2, wrong rows .95^2+.97^2+.02^2
        expected = ((0.05**2 + 0.03**2 + 0.02**2) + (0.95**2 + 0.97**2 + 0.02**2)) / 2
        self.assertAlmostEqual(M.brier_score(y, p), expected)

    def test_a_calibrated_predictor_has_near_zero_ece(self):
        y = ["increase"] * 7 + ["decrease"] * 3
        p = [{"increase": 0.7, "decrease": 0.3, "no_change": 0.0}] * 10
        self.assertAlmostEqual(M.expected_calibration_error(M.reliability_bins(y, p)), 0.0)

    def test_thin_bins_are_marked_not_meaningful(self):
        bins = M.reliability_bins(["increase"], [{"increase": 0.9, "decrease": 0.1, "no_change": 0.0}], min_n=10)
        self.assertFalse([b for b in bins if b.n][0].meaningful)

    def test_miscalibration_surfaces_through_the_harness(self):
        records = [
            synthetic.good_record(
                study_id="synthetic:c%d" % i, doi="synthetic:c%d" % i,
                direction=["increase", "decrease"][i % 2], source_locator="F%d" % i,
            )
            for i in range(40)
        ]
        bundle = cross_validate(_Overconfident, build_dataset(records), "study_id", n_bootstrap=20)
        self.assertEqual(bundle.status, "ok")
        self.assertAlmostEqual(bundle.estimates["expected_calibration_error"].value, 0.45)
        self.assertTrue(bundle.calibration["available"])

    def test_a_baseline_without_probabilities_reports_no_calibration(self):
        from syn4rch.models.baselines import NeurotransmitterSign
        bundle = cross_validate(NeurotransmitterSign, build_dataset(corpus()), "study_id", n_bootstrap=10)
        self.assertFalse(bundle.calibration["available"])
        self.assertIsNone(bundle.estimates["brier_score"].value)


class TestGroupBootstrap(unittest.TestCase):
    def test_replicates_contain_whole_groups_only(self):
        groups = [[0, 1], [2, 3], [4, 5, 6]]
        seen = []

        def metric(indices):
            seen.append(list(indices))
            return 0.0

        M.bootstrap_group_ci(metric, groups, n_boot=200, seed=7)
        self.assertEqual(len(seen), 200)
        saw_duplicate_group = False
        for indices in seen:
            for group in groups:
                counts = {i: indices.count(i) for i in group}
                self.assertEqual(len(set(counts.values())), 1, (indices, group))
                if counts[group[0]] > 1:
                    saw_duplicate_group = True
            # exactly len(groups) groups were drawn
            draws = sum(indices.count(g[0]) for g in groups)
            self.assertEqual(draws, len(groups))
        self.assertTrue(saw_duplicate_group, "resampling should be with replacement")

    def test_correlated_groups_give_a_wider_interval_than_rows_would(self):
        # 10 groups of 10: five all-correct, five all-wrong.
        correct = [1.0] * 50 + [0.0] * 50
        groups = [list(range(10 * g, 10 * g + 10)) for g in range(10)]
        rows = [[i] for i in range(100)]

        def acc(indices):
            return sum(correct[i] for i in indices) / len(indices)

        by_group = M.bootstrap_group_ci(acc, groups, n_boot=2000, seed=1)
        by_row = M.bootstrap_group_ci(acc, rows, n_boot=2000, seed=1)
        self.assertEqual(by_group["resampling_unit"], "group")
        self.assertEqual(by_group["n_groups"], 10)
        width_group = by_group["ci_high"] - by_group["ci_low"]
        width_row = by_row["ci_high"] - by_row["ci_low"]
        self.assertGreater(width_group, 2 * width_row)

    def test_one_group_gives_no_interval(self):
        out = M.bootstrap_group_ci(lambda idx: 1.0, [[0, 1, 2]], n_boot=10)
        self.assertIsNone(out["ci_low"])
        self.assertIn("single group", out["note"])

    def test_harness_reports_group_counts_not_row_counts(self):
        bundle = cross_validate(GlobalMajority, build_dataset(corpus()), "lab", n_bootstrap=50)
        ci = bundle.estimates["direction_accuracy"].ci
        self.assertEqual(ci["n_groups"], 7)
        self.assertEqual(ci["n_rows"], 40)
        self.assertEqual(ci["resampling_unit"], "group")


class TestFolds(unittest.TestCase):
    def test_builtin_folds_never_split_a_group(self):
        groups = ["a", "a", "b", "b", "b", "c", "d", "d", None]
        folds = builtin_grouped_folds(groups, 3)
        tested = []
        for train, test in folds:
            self.assertFalse({groups[i] for i in train} & {groups[i] for i in test})
            tested.extend(test)
        self.assertEqual(sorted(tested), list(range(8)))  # the None row is untested

    def test_fewer_groups_than_folds_becomes_leave_one_group_out(self):
        self.assertEqual(len(builtin_grouped_folds(["a", "b", "c"], 10)), 3)
        self.assertEqual(builtin_grouped_folds(["a", "a"], 5), [])

    def test_feasibility_splits_is_used_where_it_has_a_counterpart(self):
        provider = resolve_fold_provider()
        dataset = build_dataset(corpus())
        study = provider.folds_for(dataset, "study_id", 5)
        if not provider.record_api:  # pragma: no cover
            self.skipTest("syn4rch.feasibility.splits not importable here")
        self.assertEqual(study.source, "syn4rch.feasibility.splits:by_study")
        self.assertEqual(len(study.folds), 10)  # leave-one-study-out
        cell = provider.folds_for(dataset, "population_outcome", 5)
        self.assertEqual(cell.source, BUILTIN_SOURCE)
        self.assertIn("no counterpart", cell.note)

    def test_a_leaking_generator_is_refused(self):
        class FakeFold:
            def __init__(self, train, test):
                self.train_record_ids, self.test_record_ids = train, test

        class FakeSet:
            def __init__(self, folds):
                self.folds = folds

        dataset = build_dataset(corpus())
        ids = [e.record_id for e in dataset.examples]

        def leaky(records):  # row-level split labelled as by-study
            return FakeSet([FakeFold(ids[1:], ids[:1]), FakeFold(ids[:1], ids[1:2])])

        provider = FoldProvider(record_api={"study_id": ("by_study", leaky)})
        assignment = provider.folds_for(dataset, "study_id", 5)
        self.assertEqual(assignment.source, BUILTIN_SOURCE)
        self.assertIn("refused", assignment.note)

    def test_absent_feasibility_module_degrades_to_builtin(self):
        import syn4rch.feasibility as parent
        saved = sys.modules.get("syn4rch.feasibility.splits")
        saved_attr = parent.__dict__.pop("splits", None)
        sys.modules["syn4rch.feasibility.splits"] = None  # forces ImportError
        try:
            provider = resolve_fold_provider()
        finally:
            if saved is None:
                del sys.modules["syn4rch.feasibility.splits"]
            else:
                sys.modules["syn4rch.feasibility.splits"] = saved
            if saved_attr is not None:
                parent.splits = saved_attr
        self.assertEqual(provider.source, BUILTIN_SOURCE)
        self.assertIn("not importable", provider.note)
        bundle = cross_validate(GlobalMajority, build_dataset(corpus()), "study_id",
                                n_bootstrap=10, fold_provider=provider)
        self.assertEqual(bundle.status, "ok")
        self.assertEqual(bundle.fold_source, BUILTIN_SOURCE)


class TestHarnessStates(unittest.TestCase):
    def test_empty_dataset_produces_statements_not_numbers(self):
        result = run_benchmark(build_dataset([]), n_bootstrap=10)
        self.assertTrue(result.bundles)
        for bundle in result.bundles:
            self.assertIn(bundle.status, ("empty_dataset", "baseline_unavailable"))
            self.assertEqual(bundle.estimates, {})
        text = render_markdown(result)
        self.assertIn("There is no labellable data", text)
        self.assertNotIn("n=0/0g", text)

    def test_one_group_is_insufficient(self):
        records = [synthetic.good_record(source_locator="F%d" % i) for i in range(5)]
        bundle = cross_validate(GlobalMajority, build_dataset(records), "study_id")
        self.assertEqual(bundle.status, "insufficient_groups")
        self.assertEqual(bundle.estimates, {})

    def test_kg_path_sign_row_is_unavailable_with_no_numbers(self):
        result = run_benchmark(build_dataset(corpus()), n_bootstrap=10)
        rows = [b for b in result.bundles if b.baseline == "kg_path_sign"]
        self.assertTrue(rows)
        for bundle in rows:
            self.assertEqual(bundle.status, "baseline_unavailable")
            self.assertEqual(bundle.estimates, {})

    def test_cell_baseline_backoff_is_reported_under_cell_grouping(self):
        bundle = cross_validate(PopulationOutcomeMajority, build_dataset(corpus()),
                                "population_outcome", n_bootstrap=10)
        self.assertNotIn("population_outcome", bundle.backoff_usage)
        self.assertEqual(sum(bundle.backoff_usage.values()), 40)

    def test_results_are_deterministic(self):
        dataset = build_dataset(corpus())
        a = cross_validate(GlobalMajority, dataset, "lab", n_bootstrap=50).to_dict()
        b = cross_validate(GlobalMajority, dataset, "lab", n_bootstrap=50).to_dict()
        self.assertEqual(a["metrics"], b["metrics"])


class TestCli(unittest.TestCase):
    def test_missing_records_file_writes_a_report_that_says_so(self):
        from syn4rch.models.run import main
        out = os.path.join(tempfile.mkdtemp(), "bench.md")
        rc = main(["--records", "/nonexistent/records.jsonl", "--out", out,
                   "--bootstrap", "10", "--quiet"])
        self.assertEqual(rc, 0)
        text = open(out, encoding="utf-8").read()
        self.assertIn("There is no labellable data", text)
        self.assertIn("no such file", text)
        self.assertIn("UNAVAILABLE", text)

    def test_full_run_writes_markdown_and_json_with_groupings(self):
        import json
        from syn4rch.models.run import main
        from syn4rch.schema.io import write_jsonl
        directory = tempfile.mkdtemp()
        records = os.path.join(directory, "records.jsonl")
        write_jsonl(records, corpus() + [synthetic.good_record(direction="not_reported", source_locator="X")])
        out = os.path.join(directory, "bench.md")
        js = os.path.join(directory, "bench.json")
        rc = main(["--records", records, "--out", out, "--json", js, "--bootstrap", "20", "--quiet"])
        self.assertEqual(rc, 0)
        text = open(out, encoding="utf-8").read()
        for heading in ("`study_id`", "`lab`", "`population`", "`population_outcome`", "`record`"):
            self.assertIn("Grouping: %s" % heading, text)
        self.assertIn("The bar a model has to clear", text)
        self.assertIn("n=40/10g", text)
        data = json.load(open(js, encoding="utf-8"))
        self.assertEqual(data["dataset"]["exclusions"]["n_excluded"], 1)
        for result in data["results"]:
            self.assertIn("grouping", result)

    def test_undeclared_grouping_is_a_usage_error(self):
        from syn4rch.models.run import main
        rc = main(["--groupings", "by_vibes", "--out", os.path.join(tempfile.mkdtemp(), "x.md")])
        self.assertEqual(rc, 2)


class TestOptionalDependencyGuard(unittest.TestCase):
    def test_package_imports_with_numpy_and_sklearn_blocked(self):
        import subprocess
        src = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
        code = (
            "import sys\n"
            "class Block:\n"
            "    def find_spec(self, name, path=None, target=None):\n"
            "        if name.split('.')[0] in ('numpy', 'sklearn'):\n"
            "            raise ImportError('blocked: ' + name)\n"
            "        return None\n"
            "sys.meta_path.insert(0, Block())\n"
            "import syn4rch.models as m\n"
            "assert not m._compat.HAVE_NUMPY and not m._compat.HAVE_SKLEARN\n"
            "print(m.metrics.accuracy(['increase'], ['increase']))\n"
        )
        env = dict(os.environ, PYTHONPATH=src)
        proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout.strip(), "1.0")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
