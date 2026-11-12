"""Synthetic-only unit checks; no biological tables, arbors or model states."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

PATH = Path(__file__).resolve().parents[1] / "analysis/brm_template.py"
SPEC = importlib.util.spec_from_file_location("brm_template_tested", PATH)
M = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)


def route(sign, soma=(0., 0., 0.)):
    points = np.array([[0., 0., 0.], [500. * sign, 0., 0.], [1200. * sign, 100., 0.], [500. * sign, 350., 0.]]) + soma
    raw = M.G.Arbor(points, np.array([-1, 0, 1, 1]), np.array([1, 2, 2, 2]))
    return M.G.resample(raw, 50.)


def fixture():
    donors = []
    for i, (tag, sign, profile) in enumerate((("p", 1, 0.), ("q", -1, 4.))):
        neuron, animal = f"T{i}_001", f"T{i}"
        donors.append({"neuron": neuron, "animal": animal, "full_id": f"{neuron}_{tag}_HY", "tag": tag,
                       "table_soma": np.zeros(3), "swc_soma": np.zeros(3), "profile": np.full(49, profile), "arbor": route(sign)})
    val = [M.Query("V0_001", "V0", "V0_001_p_HY", np.zeros(3), np.zeros(3), "p"),
           M.Query("V1_001", "V1", "V1_001_q_HY", np.zeros(3), np.zeros(3), "q")]
    test = [M.Query("E0_001", "E0", "E0_001_p_HY", np.zeros(3), np.zeros(3), "p"),
            M.Query("E1_001", "E1", "E1_001_q_HY", np.zeros(3), np.zeros(3), "q")]
    continuous = [*M.Geo.R.POSITION, *M.Geo.R.target_fields()]
    train_ids = [d["full_id"] for d in donors]
    model = {"core_sha256": M.Geo.JOINT_SHA, "target_fields": list(M.Geo.R.target_fields()),
             "field_types": {**{f: "continuous" for f in continuous}, M.Geo.R.TAG: "categorical"},
             "training_ids_sha256": hashlib.sha256("\n".join(train_ids).encode()).hexdigest(),
             "field_metadata": {f: {"training_mean": 0., "training_scale": 1.} for f in continuous},
             "blocks": [{"name": "block_0", "fields": [*continuous, M.Geo.R.TAG], "weights": [.5, .5],
                         "means": {f: [0., 0.] if f in M.Geo.R.POSITION else [0., 4.] for f in continuous},
                         "variances": {f: [1., 1.] if f in M.Geo.R.POSITION else [.01, .01] for f in continuous},
                         "tables": {M.Geo.R.TAG: [[.999, .001], [.001, .999]]}}]}
    model["field_metadata"][M.Geo.R.TAG] = {"categories": ["p", "q"]}
    proof = {"name": "synthetic_train_only", "train_ids": train_ids, "train_animals": [d["animal"] for d in donors],
             "test_ids": [q.full_id for q in val + test], "test_animals": [q.animal for q in val + test],
             "model_state_sha256": M._digest(model)}
    regional = M.RegionalState(model, proof)
    identity = lambda q: {"neuron": q.neuron, "animal": q.animal, "full_id": q.full_id}
    gate_proof = {"training_animals": [d["animal"] for d in donors], "validation_animals": [q.animal for q in val],
                  "test_animals": [q.animal for q in test], "training_ids": M.donor_provenance(donors),
                  "validation_ids": list(map(identity, val)), "test_ids": list(map(identity, test)),
                  "regional_states": [{"model": model, "proof": proof}],
                  "state_by_validation_neuron": {q.neuron: proof["model_state_sha256"] for q in val},
                  "state_by_test_animal": {q.animal: proof["model_state_sha256"] for q in test}}
    return donors, regional, val, test, gate_proof


def validation_records(donors, regional, queries, *, reverse=False, same_inputs=False):
    rows = []
    for i, q in enumerate(queries):
        evidence = M.Query(q.neuron, q.animal, q.full_id, q.table_soma, q.swc_soma, "unknown") if same_inputs else q
        truth = route((-1 if i % 2 else 1) * (-1 if reverse else 1))
        metrics = {}
        for policy in M.POLICIES:
            prediction, _ = M.predict_policy(evidence, donors, regional, policy)
            value = M.G.geometry_metrics(prediction, truth)
            metrics[policy] = {key: value[key] for key in ("distance_um", "f1_200")}
        rows.append({"neuron": q.neuron, "animal": q.animal, "full_id": q.full_id, "metrics": metrics})
    return rows


def test_conditional_donor_weights_exact_formula_normalization_zero_mass():
    r = np.array([[1., 0., 0.], [0., 1., 0.], [1., 0., 0.]])
    fit = M.conditional_weights(["A", "A", "B"], r, [1., 0., 0.])
    assert fit["animal_weights"] == pytest.approx([.5, .5, 1.])
    assert fit["component_mass"] == pytest.approx([1.5, .5, 0.])
    assert fit["weights"] == pytest.approx([.3, .1, .6])
    assert fit["component_donor_weights"][:, 2] == pytest.approx([.25, .25, .5])
    assert fit["component_donor_weights"].sum(axis=0) == pytest.approx(np.ones(3))
    assert sum(fit["weights"]) == pytest.approx(1.)


def test_animal_duplication_and_component_permutation_invariance():
    r = np.array([[1., 0.], [0., 1.]])
    original = M.conditional_weights(["A", "B"], r, [.2, .8])
    duplicate = M.conditional_weights(["A", "A", "B"], r[[0, 0, 1]], [.2, .8])
    assert duplicate["weights"][:2].sum() == pytest.approx(original["weights"][0])
    assert duplicate["weights"][2] == pytest.approx(original["weights"][1])
    permuted = M.conditional_weights(["A", "B"], r[:, [1, 0]], [.8, .2])
    assert permuted["weights"] == pytest.approx(original["weights"])


def test_exact_global_baseline_and_input_donor_order_invariance():
    donors, regional, val, test, proof = fixture()
    query = test[1]
    expected = donors[M.G.training_medoid(donors)]["arbor"].shifted(query.swc_soma)
    actual, info = M.predict_policy(query, list(reversed(donors)), regional, M.BASELINE)
    assert np.array_equal(actual.xyz, expected.xyz)
    assert np.array_equal(actual.parent, expected.parent)
    assert not info["fallback"]


def test_exact_geometrically_appropriate_donor_B_recovered_without_query_truth():
    donors, regional, val, test, proof = fixture()
    prediction, info = M.predict_policy(test[1], donors, regional, "conditional_anchor_500")
    assert info["selected_donor"] == "T1_001"
    assert info["weights"][1] > info["weights"][0]
    assert info["executed_policy"] == "conditional_anchor_500"
    assert np.array_equal(prediction.xyz, donors[1]["arbor"].xyz)


def test_hidden_query_counterfactual_and_no_external_data_or_oracle_access(tmp_path, monkeypatch):
    donors, regional, val, test, proof = fixture()
    query = test[1]
    hidden = {"axon": route(-1), "profile": np.full(49, 4.)}
    np.savez(tmp_path / "hidden.npz", xyz=hidden["axon"].xyz, profile=hidden["profile"])
    def forbidden(*args, **kwargs):
        raise AssertionError("Inference tried to read external data or query truth.")
    monkeypatch.setattr(M.Geo.R.pd, "read_csv", forbidden)
    monkeypatch.setattr(np, "loadtxt", forbidden)
    monkeypatch.setattr(np, "load", forbidden)
    monkeypatch.setattr(M.Geo.G, "extract_explicit_axon", forbidden)
    before, info_before = M.predict_policy(query, donors, regional, "conditional_anchor_500")
    hidden["axon"] = route(1, [9000., 9000., 9000.])
    hidden["profile"][:] = 999999.
    np.savez(tmp_path / "hidden.npz", xyz=hidden["axon"].xyz, profile=hidden["profile"])
    after, info_after = M.predict_policy(query, donors, regional, "conditional_anchor_500")
    assert info_before["selected_donor"] == info_after["selected_donor"] == "T1_001"
    assert info_before["weights"] == info_after["weights"]
    assert np.array_equal(before.xyz, after.xyz) and np.array_equal(before.parent, after.parent)
    with pytest.raises(TypeError):
        M.Query(query.neuron, query.animal, query.full_id, query.table_soma, query.swc_soma, query.tag, arbor=hidden["axon"])
    with pytest.raises(TypeError):
        M.predict_policy(query, donors, regional, "conditional_anchor_500", query_profile=hidden["profile"])
    with pytest.raises(ValueError, match="Unknown policy"):
        M.predict_policy(query, donors, regional, "oracle_best_geometry")


def test_compact_attachment_exact_distal_coordinates_edges_topology_and_subdivision():
    a = route(1)
    query_soma = np.array([0., 100., 0.])
    attached, info = M.compact_attachment(a, query_soma, 500.)
    remap = info["canonical_vertex_map"]
    distal = a.path_lengths() >= 500.
    assert np.array_equal(attached.xyz[0], query_soma)
    assert np.array_equal(attached.xyz[np.asarray(remap)[distal]], a.xyz[distal])
    assert attached.branches() == a.branches()
    assert info["edge_ratio_min"] >= .5 and info["edge_ratio_max"] <= 1.5
    assert np.linalg.norm(attached.xyz[1:] - attached.xyz[attached.parent[1:]], axis=1).max() <= 50. + 1e-8
    assert not info["global_injectivity_or_obstacle_claim"]


def test_unsupported_attachment_returns_bitwise_exact_baseline_not_dropped_query():
    donors, regional, val, test, proof = fixture()
    query = M.Query(test[1].neuron, test[1].animal, test[1].full_id, np.array([1000., 0., 0.]), np.array([1000., 0., 0.]), "q")
    baseline, _ = M.predict_policy(query, donors, regional, M.BASELINE)
    for policy in M.POLICIES[1:]:
        actual, info = M.predict_policy(query, donors, regional, policy)
        assert info["fallback"] and info["fallback_reason"] == "soma_shift_exceeds_radius_over_three"
        assert info["executed_policy"] == M.BASELINE
        assert np.array_equal(actual.xyz, baseline.xyz) and np.array_equal(actual.parent, baseline.parent)


def test_informative_validation_selects_conditional_and_improves_independent_test_fixture():
    donors, regional, val, test, proof = fixture()
    gate = M.select_policy(validation_records(donors, regional, val), proof)
    assert gate["selected_policy"] == "conditional_anchor_500"
    assert gate["selected_has_joint_specific_validation_gain"]
    improvements = []
    for i, q in enumerate(test):
        forecast, info = M.predict_gated(q, donors, regional, gate)
        baseline, _ = M.predict_policy(q, donors, regional, M.BASELINE)
        truth = route(-1 if i % 2 else 1)
        fm, bm = M.G.geometry_metrics(forecast, truth), M.G.geometry_metrics(baseline, truth)
        improvements.append((bm["distance_um"] - fm["distance_um"], fm["f1_200"] - bm["f1_200"]))
        assert info["selected_donor"] == f"T{i}_001"
    assert np.mean(improvements, axis=0)[0] > 0 and np.mean(improvements, axis=0)[1] > 0
    json.dumps(gate, allow_nan=False)


@pytest.mark.parametrize("world", ["no_information", "misleading"])
def test_no_information_or_harmful_validation_selects_exact_baseline(world):
    donors, regional, val, test, proof = fixture()
    records = validation_records(donors, regional, val, same_inputs=world == "no_information", reverse=world == "misleading")
    gate = M.select_policy(records, proof)
    assert gate["selected_policy"] == M.BASELINE
    for q in test:
        selected, info = M.predict_gated(q, donors, regional, gate)
        expected, _ = M.predict_policy(q, donors, regional, M.BASELINE)
        assert np.array_equal(selected.xyz, expected.xyz) and np.array_equal(selected.parent, expected.parent)


def test_conditional_must_beat_matching_uniform_not_only_baseline():
    donors, regional, val, test, proof = fixture()
    records = validation_records(donors, regional, val)
    for r in records:
        for policy in M.POLICIES:
            r["metrics"][policy] = {"distance_um": 10., "f1_200": .5}
        r["metrics"]["uniform_anchor_500"] = {"distance_um": 8., "f1_200": .7}
        r["metrics"]["conditional_anchor_500"] = {"distance_um": 7., "f1_200": .6}
    gate = M.select_policy(records, proof)
    assert gate["selected_policy"] == "uniform_anchor_500"
    assert "conditional_anchor_500" not in gate["eligible_nonbaseline_policies"]
    assert not gate["selected_has_joint_specific_validation_gain"]


@pytest.mark.parametrize("bad", ["overlap", "regional_validation_leak", "regional_test_leak", "missing_policy", "missing_validation",
                                "test_metric_record", "unknown_state", "donor_fingerprint", "edited_gate", "edited_state_proof"])
def test_invalid_proofs_and_mutations_raise_instead_of_quiet_fallback(bad):
    donors, regional, val, test, proof = fixture()
    records = validation_records(donors, regional, val)
    gate = M.select_policy(records, proof)
    if bad == "overlap":
        proof["training_animals"].append(proof["validation_animals"][0])
    elif bad in ("regional_validation_leak", "regional_test_leak"):
        animal = "V0" if bad == "regional_validation_leak" else "E0"
        entry = proof["regional_states"][0]
        entry["proof"]["train_ids"].append(f"{animal}_999_p_HY")
        entry["proof"]["train_animals"].append(animal)
        entry["model"]["training_ids_sha256"] = hashlib.sha256("\n".join(entry["proof"]["train_ids"]).encode()).hexdigest()
        identity = M._digest(entry["model"])
        entry["proof"]["model_state_sha256"] = identity
        proof["state_by_validation_neuron"] = {k: identity for k in proof["state_by_validation_neuron"]}
        proof["state_by_test_animal"] = {k: identity for k in proof["state_by_test_animal"]}
    elif bad == "missing_policy":
        records[0]["metrics"].pop("uniform_anchor_500")
    elif bad == "missing_validation":
        records.pop()
    elif bad == "test_metric_record":
        records[0]["neuron"], records[0]["animal"], records[0]["full_id"] = test[0].neuron, test[0].animal, test[0].full_id
    elif bad == "unknown_state":
        proof["state_by_test_animal"]["E0"] = "unlisted"
    elif bad == "donor_fingerprint":
        donors[0]["profile"][0] += 1
        with pytest.raises(ValueError, match="Exact geometry donor"):
            M.predict_gated(test[0], donors, regional, gate)
        return
    elif bad == "edited_gate":
        gate["selected_policy"] = M.BASELINE
        gate["gate_sha256"] = M._digest({k: v for k, v in gate.items() if k != "gate_sha256"})
        with pytest.raises(ValueError, match="not reproducible"):
            M.predict_gated(test[0], donors, regional, gate)
        return
    elif bad == "edited_state_proof":
        changed = copy.deepcopy(regional.proof)
        changed["name"] = "not_the_declared_state"
        with pytest.raises(ValueError, match="complete membership proof"):
            M.predict_gated(test[0], donors, M.RegionalState(regional.model, changed), gate)
        return
    with pytest.raises(ValueError):
        M.select_policy(records, proof)
