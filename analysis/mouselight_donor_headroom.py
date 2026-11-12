"""POST-TEST DEVELOPMENT: truth-assisted capacity of the sealed 16-donor pools.

This is not a predictor or a new test. All MouseLight answers are revealed. The
distance and F1 oracles each select ONE actual attached TRAIN donor and report both
metrics; their independently optimal coordinates are never combined. No fitting,
downloads, source edits or frozen-experiment output writes occur here.

Protocol was recorded before diagnostic computations:
docs/mouselight_donor_headroom_protocol.json. Output is one inspectable JSON record.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "docs/mouselight_donor_headroom_protocol.json"
OUTPUT = ROOT / "reports/mouselight_donor_headroom.json"
EXPERIMENT = ROOT / "data/interim/mouselight_external_test/mouselight_external_geometry_v1"
BASELINE = "E1_local_uniform"
METHODS = ("sealed_E1", "truth_assisted_min_distance", "truth_assisted_max_f1")
TOLERANCE = 1e-9


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def verify_protocol(protocol):
    if protocol["schema"] != "mouselight_donor_headroom_protocol_v1":
        raise ValueError("wrong diagnostic protocol")
    for relative, expected in protocol["source_pins"].items():
        if sha(ROOT / relative) != expected:
            raise ValueError("diagnostic source pin changed: " + relative)
    primary = read(EXPERIMENT / "SCORE_COMPLETE.json")
    if not primary["scored_once"] or primary["metrics_sha256"] != sha(EXPERIMENT / "metrics.json"):
        raise ValueError("primary experiment is not an intact completed test")


def load_adapter():
    name = "mouselight_headroom_frozen_adapter"
    spec = importlib.util.spec_from_file_location(name, ROOT / "analysis/mouselight_external_test.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def select_and_classify(candidates, baseline_id):
    """One attainable donor per oracle; deterministic donor-ID tie rules."""
    candidates = sorted(candidates, key=lambda r: r["donor"])
    baseline = next(r for r in candidates if r["donor"] == baseline_id)
    distance = min(candidates, key=lambda r: (r["distance_um"], r["donor"]))
    f1 = min(candidates, key=lambda r: (-r["f1_200"], r["donor"]))
    joint = [r["donor"] for r in candidates if r["distance_um"] < baseline["distance_um"]
             and r["f1_200"] > baseline["f1_200"]]
    weak = [r["donor"] for r in candidates if r["distance_um"] <= baseline["distance_um"]
            and r["f1_200"] >= baseline["f1_200"] and
            (r["distance_um"] < baseline["distance_um"] or r["f1_200"] > baseline["f1_200"])]
    return {"sealed_E1": baseline, "truth_assisted_min_distance": distance, "truth_assisted_max_f1": f1,
            "joint_improving_donors": joint, "has_strict_both_improvement": bool(joint),
            "weakly_dominating_donors": weak, "E1_is_Pareto_efficient": not weak,
            "E1_has_no_strict_both_improver": not joint}


def sanity_check():
    # The two optima need not be the same donor. Never report the invented pair (5,.9).
    rows = [{"donor": "base", "distance_um": 10., "f1_200": .5},
            {"donor": "distance", "distance_um": 5., "f1_200": .4},
            {"donor": "f1", "distance_um": 11., "f1_200": .9},
            {"donor": "joint", "distance_um": 8., "f1_200": .7}]
    result = select_and_classify(rows, "base")
    assert result["truth_assisted_min_distance"]["f1_200"] == .4
    assert result["truth_assisted_max_f1"]["distance_um"] == 11.
    assert result["joint_improving_donors"] == ["joint"]
    assert not result["E1_is_Pareto_efficient"]
    tied = [{"donor": "b", "distance_um": 10., "f1_200": .5},
            {"donor": "a", "distance_um": 10., "f1_200": .5},
            {"donor": "c", "distance_um": 10., "f1_200": .8}]
    result = select_and_classify(tied, "b")
    assert result["truth_assisted_min_distance"]["donor"] == "a"
    assert not result["has_strict_both_improvement"] and not result["E1_is_Pareto_efficient"]
    return "passed: separate attainable pairs, strict-both versus Pareto distinction, donor-ID ties"


def aggregate(per_query):
    by_animal = {}
    for row in per_query:
        by_animal.setdefault(row["animal"], []).append(row)
    animal = []
    for sample, rows in sorted(by_animal.items()):
        record = {"animal": sample, "neurons": len(rows),
                  "queries_with_strict_both_improver": sum(r["has_strict_both_improvement"] for r in rows),
                  "fraction_with_strict_both_improver": float(np.mean([r["has_strict_both_improvement"] for r in rows])),
                  "E1_Pareto_efficient_count": sum(r["E1_is_Pareto_efficient"] for r in rows),
                  "E1_Pareto_efficient_fraction": float(np.mean([r["E1_is_Pareto_efficient"] for r in rows])),
                  "E1_no_strict_both_improver_count": sum(r["E1_has_no_strict_both_improver"] for r in rows),
                  "methods": {}}
        for method in METHODS:
            record["methods"][method] = {s: float(np.mean([r[method][s] for r in rows]))
                                            for s in ("distance_um", "f1_200")}
        animal.append(record)
    summary = {method: {s: float(np.mean([r["methods"][method][s] for r in animal]))
                        for s in ("distance_um", "f1_200")} for method in METHODS}
    for method in METHODS[1:]:
        summary[method]["distance_delta_vs_E1_um"] = summary[method]["distance_um"] - summary["sealed_E1"]["distance_um"]
        summary[method]["f1_delta_vs_E1"] = summary[method]["f1_200"] - summary["sealed_E1"]["f1_200"]
    n_joint = sum(r["has_strict_both_improvement"] for r in per_query)
    n_pareto = sum(r["E1_is_Pareto_efficient"] for r in per_query)
    counts = {"queries": len(per_query), "animals": len(animal),
              "queries_with_strict_both_improver": n_joint,
              "query_fraction_with_strict_both_improver": n_joint / len(per_query),
              "equal_animal_fraction_with_strict_both_improver": float(np.mean([r["fraction_with_strict_both_improver"] for r in animal])),
              "animals_with_at_least_one_strict_both_improver": sum(r["queries_with_strict_both_improver"] > 0 for r in animal),
              "E1_Pareto_efficient_queries": n_pareto,
              "E1_Pareto_efficient_query_fraction": n_pareto / len(per_query),
              "E1_Pareto_efficient_equal_animal_fraction": float(np.mean([r["E1_Pareto_efficient_fraction"] for r in animal])),
              "E1_no_strict_both_improver_queries": len(per_query) - n_joint,
              "E1_no_strict_both_improver_equal_animal_fraction": 1.0 - float(np.mean([r["fraction_with_strict_both_improver"] for r in animal]))}
    return animal, summary, counts


def main():
    if OUTPUT.exists():
        raise ValueError("diagnostic output already exists; preserve the completed record")
    protocol = read(PROTOCOL)
    verify_protocol(protocol)
    check = sanity_check()
    M = load_adapter()
    library = M.FrozenLibrary(ROOT)
    frozen, queries = M.verify_freeze(ROOT, library)
    E = library.E
    manifest = read(EXPERIMENT / "prediction_manifest.json")
    commit = read(EXPERIMENT / "PREDICTIONS_COMMITTED.json")
    if commit["manifest_sha256"] != sha(EXPERIMENT / "prediction_manifest.json") or manifest["freeze_sha256"] != sha(EXPERIMENT / "FROZEN.json"):
        raise ValueError("sealed prediction/freeze binding differs")
    metrics = read(EXPERIMENT / "metrics.json")
    admission = read(EXPERIMENT / "ADMISSION.json")["queries"]
    acquisition = read(EXPERIMENT / "ACQUISITION.json")["files"]
    calibration = M.calibration_of(ROOT, M.protocol_of(ROOT), queries)
    e1 = {r["query"]: r for r in manifest["predictions"] if r["method"] == BASELINE}
    primary = {r["query"]: r for r in metrics["per_neuron"] if r["method"] == BASELINE}
    expected = {q.neuron for q in queries}
    if len(queries) != 97 or len({q.animal for q in queries}) != 22 or set(e1) != expected or set(primary) != expected:
        raise ValueError("the fixed primary cohort or E1 row membership changed")
    if any(admission[n]["status"] != "included" for n in expected):
        raise ValueError("this diagnostic expected the reported all-97 fixed admission")
    per_query, max_score_error, donor_usage = [], {"distance_um": 0., "f1_200": 0.}, Counter()
    for qi, q in enumerate(queries):
        row = e1[q.neuron]
        ids = row["donor_ids"]
        if len(ids) != 16 or ids != sorted(set(ids)) or not set(ids) <= set(library.ids):
            raise ValueError("sealed E1 pool is not exactly sixteen sorted TRAIN donors")
        if primary[q.neuron]["animal"] != q.animal or row["animal"] != q.animal:
            raise ValueError("query/animal join differs")
        item = acquisition[q.neuron]
        source = M.local_path(ROOT, item["path"])
        if sha(source) != item["sha256"]:
            raise ValueError("revealed source file changed: " + q.neuron)
        truth, log = M.ingest_truth(source, q, calibration, E)
        if truth is None or log["arrays_sha256"] != admission[q.neuron]["arrays_sha256"]:
            raise ValueError("truth does not reproduce previously fixed admission")
        candidates = []
        for donor_id in ids:
            i = library.lib["index"][donor_id]
            attached, details = M.attach(library.donor(i)["arbor"], np.asarray(q.soma), E)
            metric = E.G.geometry_metrics(attached, truth)
            record = {"donor": donor_id, "distance_um": metric["distance_um"], "f1_200": metric["f1_200"], "attachment": details}
            candidates.append(record)
            donor_usage[donor_id] += 1
            if donor_id == row["selected_donor"]:
                sealed_path = M.local_path(EXPERIMENT, row["path"])
                if sha(sealed_path) != row["sha256"] or not E.same_arbor(attached, E.G.load_arbor(sealed_path)):
                    raise ValueError("E1 attachment does not reproduce sealed forecast arrays")
                for score in max_score_error:
                    error = abs(metric[score] - primary[q.neuron][score])
                    max_score_error[score] = max(max_score_error[score], error)
                    if error > TOLERANCE:
                        raise ValueError("E1 score replay differs from primary: " + q.neuron)
        classified = select_and_classify(candidates, row["selected_donor"])
        per_query.append({"query": q.neuron, "animal": q.animal, "neighbour_donor_ids": ids,
                          "candidate_scores": candidates, **classified})
        if qi % 10 == 0 or qi + 1 == len(queries):
            print(f"post-test diagnostic: {qi + 1}/{len(queries)} queries; E1 replay verified", flush=True)
    animal, summary, counts = aggregate(per_query)
    for score in max_score_error:
        if abs(summary["sealed_E1"][score] - metrics["summary"][BASELINE][score]) > TOLERANCE:
            raise ValueError("equal-animal E1 means do not reproduce the primary summary")
    verify_protocol(protocol)
    result = {"schema": "mouselight_donor_headroom_result_v1", "diagnostic_id": protocol["diagnostic_id"],
              "status": "POST_TEST_DEVELOPMENT_TRUTH_ASSISTED_NOT_A_PREDICTOR",
              "completed_utc": datetime.now(timezone.utc).isoformat(), "protocol_sha256": sha(PROTOCOL),
              "code_sha256": sha(Path(__file__)), "source_pins": protocol["source_pins"],
              "sanity_check": check, "E1_max_abs_primary_replay_error": max_score_error,
              "E1_sealed_arrays_reproduced_for_all_queries": True,
              "summary_equal_animal_means": summary, "joint_improvement_and_Pareto": counts,
              "unique_TRAIN_donors_in_fixed_neighbourhoods": len(donor_usage),
              "interpretation": {
                  "capacity": "Truth-assisted donor-selection headroom within the already fixed sixteen-neighbour pools only.",
                  "selection_implication": "A joint improver demonstrates available geometry for that query. It does not demonstrate that permitted inputs can identify it.",
                  "pool_limit": "Queries without a jointly improving donor may lack a joint improvement in this local pool; no claim is made about the full TRAIN library or newly generated shapes.",
                  "honesty": "Both metrics of each actually selected oracle donor are reported. Independently optimal distance and F1 are never combined into a fictitious forecast.",
                  "future_validation": "These MouseLight outcomes are now development data; a changed predictor needs independent untouched animals."},
              "per_animal": animal, "per_query": per_query}
    with OUTPUT.open("x") as f:
        f.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"summary": summary, "counts": counts, "E1_replay": max_score_error}, indent=2), flush=True)


if __name__ == "__main__":
    main()
