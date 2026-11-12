"""Independently replay saved completion arithmetic; never fit or rescore geometry."""
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data/interim/axon_completion_v2"
OLD = ROOT / "data/interim/axon_completion_v1"
CHECK = ROOT / "data/interim/axon_completion_v1_root"
OUTPUT = ROOT / "data/interim/axon_completion_v2_root"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def number(value):
    if isinstance(value, dict) and set(value) == {"nonfinite"}:
        return float(value["nonfinite"])
    return float(value)


def close(a, b):
    a, b = float(a), number(b)
    if math.isinf(a) or math.isinf(b):
        assert a == b, (a, b)
    else:
        assert abs(a - b) <= 1e-9, (a, b)


def mean(values):
    return math.fsum(values) / len(values)


def main():
    done = read(BASE / "RUN_COMPLETE.json")  # Refuse an incomplete experiment.
    assert sha(BASE / "results.json") == done["results_sha256"]
    result, forecasts = read(BASE / "results.json"), read(BASE / "forecasts.json")
    seal, started = read(BASE / "FORECASTS_SEALED.json"), read(BASE / "RUN_STARTED.json")
    assert sha(BASE / "forecasts.json") == seal["forecasts_sha256"]
    for path, expected in started["code_sha256"].items():
        assert sha(ROOT / path) == expected, path
    for path, expected_hash in started["input_sha256"].items():
        source = OLD / path.removeprefix("v1/") if path.startswith("v1/") else ROOT / path
        assert sha(source) == expected_hash, path
    ingest = read(ROOT / "data/interim/axon_engine/real/ingestion_log.json")["neurons"]
    expected = {r["neuron"]: r for r in csv.DictReader((CHECK / "structural_cohort.csv").open())}
    assert set(forecasts) == set(expected) and len(forecasts) == 4368
    folds = {f["outer"]: f for f in read(BASE / "folds.json")}
    old_forecasts = read(OLD / "forecasts.json")
    for q, f in forecasts.items():
        assert f["donors"]["completion_selector_v1"] == old_forecasts[q]["donors"]["completion_selector"]
        assert f["completion_sha256"]["completion_selector_v1"] == old_forecasts[q]["completion_sha256"]["completion_selector"]
        for m in old_forecasts[q]["donors"]:
            if m not in {"completion_selector", "switch_random", "permuted_labels"}:
                assert f["donors"][m] == old_forecasts[q]["donors"][m], (q, m)
                assert f["completion_sha256"][m] == old_forecasts[q]["completion_sha256"][m], (q, m)
    for k, info in read(BASE / "fold_info.json").items():
        ids = info["landmark_ids"]
        groups = [ingest[q]["animal"] for q in ids]
        fit = set(folds[int(k)]["fit_animals"])
        assert len(ids) == len(set(ids)) == len(set(groups)) == min(32, len(fit))
        assert set(groups) <= fit
        assert info["fallback_reselected"] is False
    needed = {q: set(f["donors"].values()) for q, f in forecasts.items()}
    scores = defaultdict(dict)
    for source in (OLD / "labels_fit.jsonl", OLD / "labels_score.jsonl", BASE / "labels_missing.jsonl"):
        if not source.exists():
            continue
        for line in source.open():
            row = json.loads(line)
            for pair in row["pairs"]:
                q, d = row["query"], pair["donor"]
                if q not in needed or d not in needed[q]:
                    continue
                if d in scores[q]:
                    assert scores[q][d] == pair, (q, d, "cached rows disagree")
                scores[q][d] = pair
    pair_count = sum(len(v) for v in scores.values())
    assert set(forecasts) <= set(scores)
    methods = set(result["summary"])
    per, cable = defaultdict(list), defaultdict(list)
    bad, inflated, ratios = Counter(), Counter(), defaultdict(list)
    modes = defaultdict(Counter)
    switches, fallback_equal, score_hashes = 0, 0, 0
    for q, f in forecasts.items():
        animal = ingest[q]["animal"]
        fold = folds[f["outer"]]
        assert animal in fold["eval_animals"]
        assert f["obs_sha"] == expected[q]["observation_sha256"]
        assert set(f["donors"]) == methods
        for method, donor in f["donors"].items():
            if donor is not None:
                assert donor in expected and ingest[donor]["animal"] in fold["fit_animals"]
                assert ingest[donor]["animal"] != animal
            row = scores[q][donor]
            if row.get("completion_sha256") is not None:
                assert row["completion_sha256"] == f["completion_sha256"][method]
                score_hashes += 1
            if row["valid"]:
                d, f1 = float(row["distance_um"]), float(row["f1_200"])
                assert math.isfinite(d) and d >= 0 and 0 <= f1 <= 1
            else:
                d, f1 = math.inf, 0.0
                bad[method] += 1
            per[(method, animal)].append((d, f1))
            h, e = row.get("hidden_cable_um"), row.get("extra_cable_um")
            cable[(method, animal)].append(abs(e - h) / h if h and e is not None else math.inf)
            if h:
                ratios[method].append((e or 0.0) / h)
            modes[method].update(row.get("modes", {}))
        if f["why"] == "switched":
            switches += 1
        else:
            assert f["donors"]["completion_selector"] == f["donors"]["adaptive_fallback"]
            assert f["completion_sha256"]["completion_selector"] == f["completion_sha256"]["adaptive_fallback"]
            fallback_equal += 1
    animals = sorted({a for _, a in per})
    assert len(animals) == result["samples"] == 101 and result["neurons"] == 4368
    summary, sample_means, error_means = {}, {}, {}
    max_error = 0.0
    for method in sorted(methods):
        values = np.asarray([[mean([x[j] for x in per[(method, a)]]) for j in (0, 1)] for a in animals])
        sample_means[method] = values
        means = [mean(values[:, j].tolist()) for j in (0, 1)]
        summary[method] = means
        for key, value in zip(("distance_um", "f1_200"), means):
            close(value, result["summary"][method][key])
            if math.isfinite(value):
                max_error = max(max_error, abs(value - number(result["summary"][method][key])))
        assert bad[method] == result["validity"][method]["invalid"]
        err = mean([mean(cable[(method, a)]) for a in animals])
        error_means[method] = err
        close(err, result["validity"][method]["mean_abs_rel_cable_error"])
        share = mean([float(x > 2) for x in ratios[method]]) if ratios[method] else 0.0
        inflated[method] = share
        close(share, result["validity"][method]["inflated_share"])
        assert dict(modes[method]) == result["executed_grafts"][method]
    selector = "completion_selector"
    close(switches / len(forecasts), result["verdict"]["switch_rate"])
    bootstrap = np.random.default_rng(20261110).integers(0, 101, size=(10000, 101))
    required = methods - {selector, "observation_only"}
    gate = bad[selector] == 0 and switches > 0
    for control in sorted(required):
        a, b = sample_means[selector], sample_means[control]
        if bad[control] or not (np.isfinite(a).all() and np.isfinite(b).all()):
            gate = False
            continue
        diff = a - b
        delta = diff.mean(axis=0)
        interval = np.percentile(diff[bootstrap].mean(axis=1), [2.5, 97.5], axis=0)
        original = result["versus"][control]
        close(delta[0], original["distance_delta_um"])
        close(delta[1], original["f1_delta"])
        for j, key in enumerate(("distance_bootstrap95_um", "f1_bootstrap95")):
            for i in (0, 1):
                close(interval[i, j], original[key][i])
        gate &= bool(delta[0] <= -10 and delta[1] >= .010 and interval[1, 0] < 0 and interval[0, 1] > 0)
        gate &= inflated[selector] <= inflated[control] + .02
    gate &= error_means[selector] <= error_means["adaptive_fallback"]
    assert bool(gate) == result["verdict"]["pass"]
    assert sha(ROOT / "data/interim/axon_engine/real/test/metrics.json") == "ae18960e168636b9af9b394b9bc7b20098d2a297de594f3cd8e65d5eff2e351b"
    receipt = {"utc": datetime.now(timezone.utc).isoformat(), "passed": True,
               "kind": "Saved scores, predictions and bootstrap arithmetic; no model fit or geometry rescore",
               "checker_sha256": sha(__file__), "results_sha256": sha(BASE / "results.json"),
               "forecasts_sha256": seal["forecasts_sha256"], "queries": len(forecasts),
               "brain_samples": len(animals), "scored_pairs": pair_count,
               "matched_forecast_score_hashes": score_hashes, "exact_fallback_forecasts": fallback_equal,
               "switched_queries": switches, "max_summary_abs_difference": max_error,
               "all_fold_and_observation_checks_passed": True,
               "fit_only_unique_animal_landmarks": True, "fixed_v1_controls_unchanged": True, "all_frozen_sources_unchanged": True,
               "old_test_unchanged": True, "recomputed_development_gate": bool(gate)}
    OUTPUT.mkdir(parents=True, exist_ok=True)
    path = OUTPUT / "saved_result_check.json"
    path.write_text(json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False) + "\n")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
