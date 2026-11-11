"""Interface tests for analysis/prediction_viewer_export.py. Synthetic fixtures only; no biological data.

The "winning" fixture's candidate is the hidden truth on purpose: only the export interface and its
refusals are tested here, never a method. The fixture declares an upstream fit on its two training
animals, which the referee accepts since its fix of 2026-11-11
(reports/prediction_referee_interop_finding.md).
"""
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if not (ROOT / "analysis/prediction_independent_evaluator.py").exists():
    pytest.skip("the independent referee (another agent's file) is not in this checkout", allow_module_level=True)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


E = load("referee_for_viewer_export_test", "analysis/prediction_independent_evaluator.py")
X = load("viewer_export_tested", "analysis/prediction_viewer_export.py")

COHORT_ROWS = [("d1", "A", "train"), ("d2", "A", "train"), ("d3", "B", "train"), ("d9", "A", "train"), ("v1", "V", "validation"),
               ("q1", "T1", "test"), ("q2", "T2", "test"), ("q3", "T2", "test"), ("qx", "T3", "test")]


def line(soma, axis, length=1000.0):
    n = int(length / 50) + 1
    xyz = np.zeros((n, 3))
    xyz[:, axis] = np.arange(n) * 50.0
    return xyz + np.array(soma, dtype=float), np.arange(-1, n - 1), np.array([1] + [2] * (n - 1))


def locked(tmp, rows=COHORT_ROWS, pin=True):
    manifest = tmp / "cohort.tsv"
    manifest.write_text("neuron\tanimal\tsource_label\tpart\n" + "".join("%s\t%s\tx\t%s\n" % r for r in rows))
    plan = tmp / "plan.md"
    plan.write_text("cohort manifest SHA-256 `%s`\n" % (X.sha(manifest) if pin else "0" * 64))
    return plan, manifest


def bundle(base, *, win=True, matched=True, excluded=("qx",), score=True):
    base.mkdir(parents=True, exist_ok=True)
    donors, train = [], []
    for name, animal, soma, tag, axis in [("d1", "A", [0.0, 0.0, 0.0], "P", 0), ("d2", "A", [0.0, 100.0, 0.0], "P", 1),
                                           ("d3", "B", [200.0, 0.0, 0.0], "Q", 1)]:
        tree, relative = line(soma, axis), "train/%s.npz" % name
        E.save_tree(base / relative, tree)
        row = {"id": name, "animal": animal, "soma_um": soma, "peptide": tag, "arbor_path": relative,
               "arbor_sha256": E.sha(base / relative)}
        train.append(row)
        donors.append({**row, "arbor": tree})
    queries = [{"id": "q1", "animal": "T1", "soma_um": [10.0, 0.0, 0.0], "peptide": "P", "nucleus": "N"},
               {"id": "q2", "animal": "T2", "soma_um": [190.0, 0.0, 0.0], "peptide": None, "nucleus": None},
               {"id": "q3", "animal": "T2", "soma_um": [195.0, 0.0, 0.0], "peptide": "Q", "nucleus": "N"}]
    method = base / "method/protocol.txt"
    method.parent.mkdir()
    method.write_text("synthetic interface fixture\n")
    controls = ["matched_uniform"] if matched else []
    declaration = {
        "geometry_contract": {"units": "um", "axes": ["AP", "DV", "ML"], "frame": "synthetic_same_frame",
                              "admission": "all fixture rows", "preprocessing": "50 um spacing", "clipping": "none"},
        "method_artifacts": {"method/protocol.txt": E.sha(method)},
        "splits": {"train": train, "validation": [{"id": "v1", "animal": "V", "soma_um": [0.0, 0.0, 0.0]}], "test": queries},
        "excluded_ids": list(excluded), "test_excluded_animals": [], "upstream_fit_animals": ["A", "B"],
        "models": {"candidate": "candidate", "copy_baselines": list(E.COPIES), "matched_controls": controls},
        "config": dict(E.CONFIG)}
    E.freeze_declaration(base, declaration)
    copies, rows, truth = E.copying_baselines(donors, queries), [], {}
    for q in queries:
        truth[q["id"]] = line(q["soma_um"], 2)
        for model in ["candidate", *E.COPIES, *controls]:
            source = model if model in E.COPIES else "medoid"
            tree = truth[q["id"]] if (model == "candidate" and win) else copies["predictions"][source][q["id"]]
            relative = "predictions/%s/%s.npz" % (model, q["id"])
            E.save_tree(base / relative, tree)
            row = {"query_id": q["id"], "model": model, "path": relative, "sha256": E.sha(base / relative),
                   "donor_ids": [copies["donor_proof"][source][q["id"]]], "upstream_fit_animals": ["A", "B"]}
            if model == "candidate":
                row.update(hemisphere="Right" if q["id"] == "q1" else "Left", fallback=False, bulky=list(range(2000)))
            rows.append(row)
    E.seal_predictions(base, {"declaration_sha256": E.sha(base / "declaration.json"), "predictions": rows,
                              "upstream_fit_animals": ["A", "B"]})
    if score:
        E.score_once(base, lambda q: truth[q])
    return declaration


def test_passed_test_exports_three_labelled_layers_for_every_scored_neuron(tmp_path):
    bundle(tmp_path / "b")
    plan, manifest = locked(tmp_path)
    out = tmp_path / "out/export.json"
    export = X.export_paths(tmp_path / "b", out, n_initial=12, plan=plan, manifest=manifest)
    assert json.loads(out.read_text()) == json.loads(json.dumps(export))
    assert (tmp_path / "out/export.json.sha256").read_text().strip() == X.sha(out)
    assert export["kind"] == "paths" and export["benchmark_status"]["evidence_level"] == "independent_test"
    assert [n["id"] for n in export["neurons"]] == ["q1", "q2", "q3"]
    metrics = json.loads((tmp_path / "b/metrics.json").read_text())
    summary = metrics["overall"]["summary"]
    strongest = min(E.COPIES, key=lambda m: (summary[m]["distance_um"], -summary[m]["f1_200"], m))
    assert export["roles"]["copy_layer"] == strongest
    for n in export["neurons"]:
        assert set(n["layers"]) == {"measured", "prediction", "copy_baseline"}
        assert n["layers"]["prediction"]["model"] == "candidate" and n["layers"]["copy_baseline"]["model"] == strongest
        assert "not measured" in n["layers"]["prediction"]["status"] and "TRAINING" in n["layers"]["copy_baseline"]["status"]
        for layer in n["layers"].values():
            assert X.sha(tmp_path / "b" / layer["path"]) == layer["sha256"]
        assert set(n["other_controls"]) == set(E.COPIES) - {strongest} | {"matched_uniform"}
        assert set(n["metrics"]) == {"candidate", *E.COPIES, "matched_uniform"}
        assert n["uncertainty"].startswith("uncalibrated")
        assert n["comparison"]["verdict"] == "prediction better on both measures"
        assert n["layers"]["prediction"]["upstream_fit_animals"]["count"] == 2
        assert n["layers"]["prediction"]["provenance"]["fallback"] is False
        assert n["layers"]["prediction"]["provenance_fields_omitted_as_large"] == ["bulky"]
    assert {n["id"]: n["hemisphere"] for n in export["neurons"]}["q1"] == "Right"
    assert export["benchmark_status"]["coverage"] == {"plan_test_neurons": 4, "plan_test_brain_samples": 3, "scored_neurons": 3,
                                                     "scored_brain_samples": 2, "declared_exclusions_in_test": ["qx"]}
    assert export["hashes"]["metrics_sha256"] == X.sha(tmp_path / "b/metrics.json")
    assert export["hashes"]["cohort_manifest_sha256"] == X.sha(manifest)
    assert sum(n["initial_example"] for n in export["neurons"]) == 2  # one per brain sample


def test_losing_candidate_is_refused_and_nothing_is_written(tmp_path):
    bundle(tmp_path / "b", win=False)
    plan, manifest = locked(tmp_path)
    with pytest.raises(X.ExportRefused, match="success condition not met"):
        X.export_paths(tmp_path / "b", tmp_path / "export.json", plan=plan, manifest=manifest)
    assert not (tmp_path / "export.json").exists()


def test_edited_metrics_cannot_turn_a_loss_into_an_export(tmp_path):
    bundle(tmp_path / "b", win=False)
    plan, manifest = locked(tmp_path)
    path = tmp_path / "b/metrics.json"
    metrics = json.loads(path.read_text())
    metrics["overall"]["beats_all_copying"] = metrics["overall"]["beats_all_matched"] = True
    path.write_text(E.canonical(metrics))
    with pytest.raises(X.ExportRefused, match="changed after scoring"):
        X.export_paths(tmp_path / "b", tmp_path / "export.json", plan=plan, manifest=manifest)


def test_win_without_a_matched_control_is_refused(tmp_path):
    bundle(tmp_path / "b", matched=False)
    plan, manifest = locked(tmp_path)
    with pytest.raises(X.ExportRefused, match="matched control"):
        X.export_paths(tmp_path / "b", tmp_path / "export.json", plan=plan, manifest=manifest)


@pytest.mark.parametrize("case", ["development_animal_scored", "test_neuron_dropped", "manifest_not_pinned", "unknown_neuron"])
def test_only_the_whole_locked_test_can_be_exported(tmp_path, case):
    bundle(tmp_path / "b", excluded=() if case == "test_neuron_dropped" else ("qx",))
    rows = list(COHORT_ROWS)
    if case == "development_animal_scored":
        rows = [(n, a, "train" if a == "T1" else p) for n, a, p in rows]
    if case == "unknown_neuron":
        rows = [r for r in rows if r[0] != "q3"]
    plan, manifest = locked(tmp_path, rows, pin=case != "manifest_not_pinned")
    with pytest.raises(X.ExportRefused):
        X.export_paths(tmp_path / "b", tmp_path / "export.json", plan=plan, manifest=manifest)
    assert not (tmp_path / "export.json").exists()


@pytest.mark.parametrize("case", ["forecast", "truth", "unscored", "failed"])
def test_changed_or_incomplete_bundles_are_refused(tmp_path, case):
    bundle(tmp_path / "b", score=case != "unscored")
    plan, manifest = locked(tmp_path)
    if case == "forecast":
        (tmp_path / "b/predictions/candidate/q2.npz").write_bytes(b"changed")
    elif case == "truth":
        (tmp_path / "b/truth/q2.npz").write_bytes(b"changed")
    elif case == "failed":
        (tmp_path / "b/SCORE_FAILED.json").write_text("{}")
    with pytest.raises(X.ExportRefused):
        X.export_paths(tmp_path / "b", tmp_path / "export.json", plan=plan, manifest=manifest)
    assert not (tmp_path / "export.json").exists()


def test_first_examples_depend_on_identifiers_only():
    queries = [{"id": "%s_%d" % (a, i), "animal": a} for a in ("A", "B", "C", "D") for i in range(3)]
    first = X.initial_examples(queries, 3)
    assert len(first) == 3 and len({q.split("_")[0] for q in first}) == 3
    scored = [{**q, "models": {"candidate": {"distance_um": float(i)}}} for i, q in enumerate(reversed(queries))]
    assert X.initial_examples(scored, 3) == first
    assert set(X.initial_examples(queries, 99)) >= set(first) and len(X.initial_examples(queries, 99)) == 4


def test_command_line_refuses_a_bundle_outside_the_locked_cohort(tmp_path, capsys):
    bundle(tmp_path / "b")
    assert X.main(["paths", "--bundle", str(tmp_path / "b"), "--out", str(tmp_path / "export.json")]) == 3
    assert "REFUSED" in capsys.readouterr().out and not (tmp_path / "export.json").exists()


def regional(tmp_path, **status):
    structures = tmp_path / "structures.json"
    structures.write_text(json.dumps({"msg": [{"acronym": "root", "id": 997, "children": [
        {"acronym": "HY", "id": 1097, "children": [{"acronym": "ARH", "id": 223, "children": []}]}]}]}))
    payload = {"task": "regional profile from soma and source label", "scale": "log10(1 + axon length in um)",
               "regions": ["HY", "ARH"],
               "status": {"beats_all_baselines": True, "evidence_level": "exploratory_training_partition",
                          "benchmark": "frozen exploratory benchmark, 101 training brain samples", "protocol_sha256": "a" * 64,
                          "results_sha256": "b" * 64, "model_sha256": "c" * 64, **status},
               "units": [{"unit_id": "n1", "group_id": "A", "predicted": [1.0, 2.0]},
                         {"unit_id": "n2", "group_id": "B", "predicted": [0.5, 0.1], "interval_low": [0.0, 0.0],
                          "interval_high": [1.0, 1.0], "coverage_validated": True}]}
    path = tmp_path / "regional.json"
    path.write_text(json.dumps(payload))
    return path, structures, payload


def test_regional_success_becomes_a_labelled_heat_map_and_never_a_path(tmp_path):
    path, structures, _ = regional(tmp_path)
    export = X.export_regional(path, tmp_path / "regional_export.json", structures=structures)
    assert export["kind"] == "regional" and "neurons" not in export and "layers" not in json.dumps(export)
    assert "Regional prediction" in export["labels"]["heatmap"] and "exploratory" in export["labels"]["heatmap"]
    assert export["regions"] == [{"acronym": "HY", "ccf_structure_id": 1097}, {"acronym": "ARH", "ccf_structure_id": 223}]
    assert [u["uncertainty"].split(":")[0] for u in export["units"]] == ["uncalibrated", "interval with validated coverage"]
    assert export["benchmark_status"]["model_sha256"] == "c" * 64


@pytest.mark.parametrize("status", [{"beats_all_baselines": False}, {"beats_all_baselines": "yes"}, {"evidence_level": "strong"},
                                    {"model_sha256": ""}, {"benchmark": None}])
def test_regional_export_is_refused_without_a_pass_and_its_exact_status(tmp_path, status):
    path, structures, _ = regional(tmp_path, **status)
    with pytest.raises(X.ExportRefused):
        X.export_regional(path, tmp_path / "regional_export.json", structures=structures)
    assert not (tmp_path / "regional_export.json").exists()


def test_regional_export_refuses_unknown_regions_and_ragged_values(tmp_path):
    path, structures, payload = regional(tmp_path)
    payload["regions"] = ["HY", "NOT_A_REGION"]
    path.write_text(json.dumps(payload))
    with pytest.raises(X.ExportRefused, match="not CCF structures"):
        X.export_regional(path, tmp_path / "regional_export.json", structures=structures)
    payload["regions"] = ["HY", "ARH"]
    payload["units"][0]["predicted"] = [1.0]
    path.write_text(json.dumps(payload))
    with pytest.raises(X.ExportRefused, match="do not match"):
        X.export_regional(path, tmp_path / "regional_export.json", structures=structures)
