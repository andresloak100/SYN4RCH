"""Bounded conditional-template selector, currently for synthetic testing only.

The query supplies soma/tag/identity, never its measured axon or profile.
Conditional donor weights feed a discrete *geometric* medoid decision. Local
attachment preserves distal CCF vertices and has a declared support bound.
The validation gate may select the exact old translated medoid instead; this
is a validation rule, not a guarantee of improvement on future animals.
No model fit or biological data read occurs on import or inside this module.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
import hashlib
import importlib.util
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
GEOMETRY_SHA = "8de317bc4f261e8061b0bccad1dbba537b2a8daf6f6834d94d303645d6611187"
_path = ROOT / "analysis/brm_geometry.py"
if hashlib.sha256(_path.read_bytes()).hexdigest() != GEOMETRY_SHA:
    raise RuntimeError("Pinned geometry dependency changed.")
_spec = importlib.util.spec_from_file_location("brm_template_geometry", _path)
Geo = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = Geo
_spec.loader.exec_module(Geo)
G = Geo.G
Query, RegionalState = Geo.Query, Geo.RegionalState
BASELINE = "baseline"
POLICIES = (BASELINE, "uniform_anchor_500", "conditional_anchor_500",
            "uniform_anchor_1000", "conditional_anchor_1000")
CONFIG = {"policies": list(POLICIES), "lambda_animals": 1., "anchor_radii_um": [500., 1000.],
          "step_um": 50., "maximum_shift_fraction": 1 / 3,
          "edge_ratio_bounds": [.5, 1.5], "pair_cost_frame": "soma_centered",
          "medoid_ties": "ascending_neuron_ID", "selection": "validation_animal_mean_both_metrics",
          "query_outcomes": False, "status": "synthetic_only_not_biologically_validated"}


def _digest(value):
    return hashlib.sha256(Geo.canonical(value).encode()).hexdigest()


def source_provenance():
    return {"code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "geometry_sha256": GEOMETRY_SHA, "geometry_metric_sha256": Geo.CORE_SHA,
            "regional_sha256": Geo.REGIONAL_SHA, "joint_sha256": Geo.JOINT_SHA}


def donor_provenance(donors):
    """Bind the exact training pool, including measured geometry and profiles."""
    result = []
    for d in sorted(donors, key=lambda r: r["neuron"]):
        a = d["arbor"].validate()
        content = {key: d[key] for key in ("neuron", "animal", "full_id", "tag")}
        content.update({"table_soma": np.asarray(d["table_soma"]).tolist(),
                        "swc_soma": np.asarray(d["swc_soma"]).tolist(),
                        "profile": np.asarray(d["profile"]).tolist(),
                        "xyz_um": a.xyz.tolist(), "parent": a.parent.tolist(), "types": a.types.tolist()})
        result.append({"neuron": d["neuron"], "animal": d["animal"], "full_id": d["full_id"],
                       "donor_sha256": _digest(content)})
    return result


def conditional_weights(animals, responsibilities, posterior):
    """Unit-animal prior shrinks every component to the same donor prior."""
    animals, r, pi = np.asarray(animals), np.asarray(responsibilities, float), np.asarray(posterior, float)
    if len(animals) == 0 or animals.ndim != 1 or r.ndim != 2 or r.shape[0] != len(animals) or pi.shape != (r.shape[1],):
        raise ValueError("Donor/responsibility/posterior shapes differ.")
    if any(not isinstance(a, (str, np.str_)) or not a for a in animals):
        raise ValueError("Animal IDs must be nonempty strings.")
    if not np.isfinite(r).all() or not np.isfinite(pi).all() or (r < 0).any() or (pi < 0).any() or not np.allclose(r.sum(axis=1), 1.) or not np.isclose(pi.sum(), 1.):
        raise ValueError("Responsibilities/posterior must be normalized probabilities.")
    counts = Counter(animals)
    a = np.array([1. / counts[x] for x in animals])
    A = len(counts)
    prior = a / A
    masses = (a[:, None] * r).sum(axis=0)
    v = (a[:, None] * r + CONFIG["lambda_animals"] * prior[:, None]) / (masses + CONFIG["lambda_animals"])
    w = v @ pi
    if not np.allclose(v.sum(axis=0), 1., atol=1e-12) or not np.isclose(w.sum(), 1., atol=1e-12) or (w < 0).any():
        raise AssertionError("Conditional donor probabilities failed normalization.")
    return {"weights": w, "animal_weights": a, "animal_prior": prior, "effective_animals": A,
            "component_mass": masses, "component_donor_weights": v}


def pair_costs(donors):
    """Unchanged sampled-curve metric on soma-centred canonical training trees."""
    if not donors:
        raise ValueError("No geometry donors.")
    centered = [d["arbor"].shifted(np.zeros(3)) for d in donors]
    n = len(donors)
    distance, f1_loss = np.zeros((n, n)), np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            result = G.geometry_metrics(centered[i], centered[j])
            distance[i, j] = distance[j, i] = result["distance_um"]
            f1_loss[i, j] = f1_loss[j, i] = 1. - result["f1_200"]
    return distance, f1_loss


def weighted_medoid(distance, weights):
    D, w = np.asarray(distance, float), np.asarray(weights, float)
    if D.ndim != 2 or D.shape != (len(w), len(w)) or len(w) == 0 or not np.isfinite(D).all() or (D < 0).any() or not np.isfinite(w).all() or (w < 0).any() or not np.isclose(w.sum(), 1.):
        raise ValueError("Invalid discrete geometric risk matrix/weights.")
    risk = D @ w
    return int(np.argmin(risk)), risk


class UnsupportedAttachment(ValueError):
    """A valid query outside the declared geometric support, not bad evidence."""

    def __init__(self, reason, details):
        super().__init__(reason)
        self.reason, self.details = reason, details


def compact_attachment(arbor, soma, radius):
    """Geodesic smoothstep; distal vertices remain at exact measured coordinates."""
    if radius not in CONFIG["anchor_radii_um"]:
        raise ValueError("Attachment radius is outside the fixed two-value family.")
    arbor.validate()
    soma = np.asarray(soma, float)
    if soma.shape != (3,) or not np.isfinite(soma).all():
        raise ValueError("Known soma must be a finite three-vector.")
    delta = soma - arbor.xyz[0]
    shift = float(np.linalg.norm(delta))
    details = {"radius_um": float(radius), "soma_shift_um": shift, "maximum_supported_shift_um": radius / 3}
    if shift > radius / 3:
        raise UnsupportedAttachment("soma_shift_exceeds_radius_over_three", details)
    depth = arbor.path_lengths()
    t = np.minimum(depth / radius, 1.)
    taper = 1 - 3 * t**2 + 2 * t**3
    taper[depth >= radius] = 0.  # exact distal identity, not a tiny residual shift
    points = arbor.xyz + taper[:, None] * delta
    points[0] = soma.copy()
    if len(points) > 1:
        old = arbor.xyz[1:] - arbor.xyz[arbor.parent[1:]]
        new = points[1:] - points[arbor.parent[1:]]
        lengths = np.linalg.norm(old, axis=1)
        positive = lengths > 0
        ratio = np.linalg.norm(new[positive], axis=1) / lengths[positive]
        projection = np.einsum("ij,ij->i", old[positive], new[positive])
        details.update({"edge_ratio_min": float(ratio.min()) if len(ratio) else None,
                        "edge_ratio_max": float(ratio.max()) if len(ratio) else None,
                        "zero_length_edges": int((~positive).sum())})
        if not np.isfinite(points).all() or np.any(ratio < .5 - 1e-10) or np.any(ratio > 1.5 + 1e-10) or np.any(projection <= 0):
            raise UnsupportedAttachment("edge_distortion_support_check_failed", details)
        if np.any(np.linalg.norm(new[~positive], axis=1) > 1e-8):
            raise UnsupportedAttachment("zero_length_edge_changed", details)
    if not np.array_equal(points[depth >= radius], arbor.xyz[depth >= radius]):
        raise AssertionError("Distal anatomical coordinates moved.")
    output, remap = Geo.subdivide_preserving_vertices(G.Arbor(points, arbor.parent.copy(), arbor.types.copy()), step=CONFIG["step_um"])
    details.update({"canonical_vertex_map": remap.tolist(), "distal_vertex_count": int((depth >= radius).sum()),
                    "topology_preserved": output.branches() == arbor.branches(),
                    "global_injectivity_or_obstacle_claim": False})
    return output, details


def _training_boundary(query, donors, regional):
    Geo.assert_query_safe(query, regional)
    if not donors or len({d["neuron"] for d in donors}) != len(donors) or len({d["full_id"] for d in donors}) != len(donors):
        raise ValueError("Distinct eligible geometry donors are required.")
    for d in donors:
        if d["animal"] == query.animal or d["neuron"] == query.neuron or d["full_id"] == query.full_id:
            raise ValueError("Query outcome/animal entered geometry donors.")
        if not d["full_id"].startswith(d["neuron"] + "_") or d["animal"] != d["neuron"].split("_", 1)[0]:
            raise ValueError("Donor identity/animal provenance differs.")
        a = d["arbor"].validate()
        if not np.allclose(a.xyz[0], d["swc_soma"], atol=1e-6, rtol=0):
            raise ValueError("Donor geometry has a different known soma.")
        if len(a.xyz) > 1 and np.linalg.norm(a.xyz[1:] - a.xyz[a.parent[1:]], axis=1).max() > 50. + 1e-6:
            raise ValueError("Donor geometry is not canonical50um.")


def predict_policy(query, donors, regional, policy):
    """One deterministic whole-arbor forecast; no query truth is accepted."""
    if policy not in POLICIES:
        raise ValueError("Unknown policy; the fixed family has five members.")
    donors = sorted(donors, key=lambda d: d["neuron"])
    _training_boundary(query, donors, regional)
    baseline_index = G.training_medoid(donors)
    baseline = donors[baseline_index]["arbor"].shifted(query.swc_soma)
    info = {"requested_policy": policy, "executed_policy": BASELINE, "fallback": policy != BASELINE,
            "fallback_reason": None, "baseline_donor": donors[baseline_index]["neuron"],
            "selected_donor": donors[baseline_index]["neuron"],
            "donor_ids": [d["neuron"] for d in donors], "donor_animals": [d["animal"] for d in donors],
            "upstream_state_sha256": regional.proof["model_state_sha256"],
            "source_provenance": source_provenance(), "config": deepcopy(CONFIG)}
    if policy == BASELINE:
        return baseline, info
    conditional = policy.startswith("conditional_")
    radius = float(policy.rsplit("_", 1)[1])
    animals = [d["animal"] for d in donors]
    posterior = Geo.state_posterior(regional.model, np.asarray(query.table_soma)[None], [query.tag])[0]
    responsibilities = Geo.state_posterior(regional.model, np.asarray([d["table_soma"] for d in donors]),
                                          [d["tag"] for d in donors], measured_profiles=np.asarray([d["profile"] for d in donors]))
    fit = conditional_weights(animals, responsibilities, posterior)
    weights = fit["weights"] if conditional else fit["animal_prior"]
    D, F = pair_costs(donors)
    selected, distance_risk = weighted_medoid(D, weights)
    if not conditional:
        selected = baseline_index  # exact existing tie/selection rule
    info.update({"attempted_donor": donors[selected]["neuron"], "radius_um": radius,
                 "guidance": "conditional" if conditional else "uniform",
                 "query_posterior": posterior.tolist(), "donor_responsibilities": responsibilities.tolist(),
                 "weights": weights.tolist(), "animal_prior": fit["animal_prior"].tolist(),
                 "animal_weights": fit["animal_weights"].tolist(), "component_mass": fit["component_mass"].tolist(),
                 "component_donor_weights": fit["component_donor_weights"].tolist(),
                 "pair_distance": D.tolist(), "pair_f1_loss": F.tolist(),
                 "distance_risk": distance_risk.tolist(), "f1_risk": (F @ weights).tolist()})
    try:
        output, attachment = compact_attachment(donors[selected]["arbor"], query.swc_soma, radius)
    except UnsupportedAttachment as unsupported:
        info.update({"fallback_reason": unsupported.reason, "attachment": unsupported.details})
        return baseline, info  # EXACT baseline coordinates/topology, no dropped query
    info.update({"executed_policy": policy, "selected_donor": donors[selected]["neuron"], "fallback": False,
                 "attachment": attachment})
    return output, info


def _identities(records):
    return sorted((r["neuron"], r["animal"], r["full_id"]) for r in records)


def _validate_selection_proof(records, proof):
    keys = {"training_animals", "validation_animals", "test_animals", "training_ids", "validation_ids", "test_ids", "regional_states",
            "state_by_validation_neuron", "state_by_test_animal"}
    if not isinstance(proof, dict) or set(proof) != keys:
        raise ValueError("Selection proof must have the explicit train/validation/test/state schema.")
    sets = []
    for name in ("training_animals", "validation_animals", "test_animals"):
        values = proof[name]
        if not isinstance(values, list) or not values or len(set(values)) != len(values) or any(not isinstance(a, str) or not a for a in values):
            raise ValueError("Every split needs distinct, nonempty animal IDs.")
        sets.append(set(values))
    train, validation, test = sets
    if train & validation or train & test or validation & test:
        raise ValueError("Training, validation and final-test animals overlap.")
    for split, identities in ((train, proof["training_ids"]), (validation, proof["validation_ids"]), (test, proof["test_ids"])):
        if not identities or len({r["neuron"] for r in identities}) != len(identities) or len({r["full_id"] for r in identities}) != len(identities) or {r["animal"] for r in identities} != split:
            raise ValueError("Every split must bind exact distinct neuron identities.")
        for r in identities:
            if not r["full_id"].startswith(r["neuron"] + "_") or r["animal"] != r["neuron"].split("_", 1)[0]:
                raise ValueError("Declared split identity/animal fields differ.")
    for r in proof["training_ids"]:
        value = r.get("donor_sha256")
        if set(r) != {"neuron", "animal", "full_id", "donor_sha256"} or not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
            raise ValueError("Training donor numeric/source fingerprint is required.")
    if not records or len({r["neuron"] for r in records}) != len(records) or len({r["full_id"] for r in records}) != len(records):
        raise ValueError("Distinct complete validation records are required.")
    if _identities(records) != _identities(proof["validation_ids"]) or {r["animal"] for r in records} != validation:
        raise ValueError("Validation identities/animals do not exactly match the declared set.")
    if set(proof["state_by_validation_neuron"]) != {r["neuron"] for r in records} or set(proof["state_by_test_animal"]) != test:
        raise ValueError("Every validation neuron and final-test animal needs a declared regional state.")
    states = {}
    for entry in proof["regional_states"]:
        if set(entry) != {"model", "proof"}:
            raise ValueError("Regional state proof needs numeric model and exact membership.")
        regional = RegionalState(entry["model"], entry["proof"])
        p = regional.proof
        identity = p["model_state_sha256"]
        if identity in states or _digest(regional.model) != identity:
            raise ValueError("Duplicate or changed regional numeric state.")
        ids = p["train_ids"]
        if len(ids) != len(set(ids)) or {i.split("_", 1)[0] for i in ids} != set(p["train_animals"]):
            raise ValueError("Regional training IDs do not prove the listed training animals.")
        if regional.model["training_ids_sha256"] != hashlib.sha256("\n".join(ids).encode()).hexdigest():
            raise ValueError("Regional numeric state does not bind training IDs.")
        if set(p["train_animals"]) & (validation | test):
            raise ValueError("Regional state trained on validation or final-test animal outcomes.")
        states[identity] = regional
    referenced = set(proof["state_by_validation_neuron"].values()) | set(proof["state_by_test_animal"].values())
    if not states or referenced != set(states):
        raise ValueError("Unknown, missing or unused regional state in selection proof.")
    for record in records:
        if not record["full_id"].startswith(record["neuron"] + "_") or record["animal"] != record["neuron"].split("_", 1)[0]:
            raise ValueError("Validation identity/animal fields differ.")
        regional = states[proof["state_by_validation_neuron"][record["neuron"]]]
        query = Query(record["neuron"], record["animal"], record["full_id"], np.zeros(3), np.zeros(3), "proof_only")
        Geo.assert_query_safe(query, regional)  # membership/parameter check; no inference
    for animal, identity in proof["state_by_test_animal"].items():
        if animal not in states[identity].proof["test_animals"]:
            raise ValueError("Declared final-test animal is not excluded by its own regional state.")
    for record in proof["test_ids"]:
        regional = states[proof["state_by_test_animal"][record["animal"]]]
        query = Query(record["neuron"], record["animal"], record["full_id"], np.zeros(3), np.zeros(3), "proof_only")
        Geo.assert_query_safe(query, regional)
    return states


def select_policy(validation_records, proof):
    """Gate uses validation metrics only and retains exact baseline if no gain."""
    records, proof = deepcopy(validation_records), deepcopy(proof)
    _validate_selection_proof(records, proof)
    for r in records:
        if set(r) != {"neuron", "animal", "full_id", "metrics"} or set(r["metrics"]) != set(POLICIES):
            raise ValueError("Every declared validation query needs exactly all five policy metrics.")
        for metrics in r["metrics"].values():
            if set(metrics) != {"distance_um", "f1_200"}:
                raise ValueError("Gate accepts only declared validation distance/F1 metrics.")
            d, f = metrics["distance_um"], metrics["f1_200"]
            if isinstance(d, bool) or isinstance(f, bool) or not math.isfinite(d) or not math.isfinite(f) or d < 0 or not 0 <= f <= 1:
                raise ValueError("Invalid validation metrics.")
    animals = sorted(proof["validation_animals"])
    summary = {}
    per_animal = {}
    for policy in POLICIES:
        per_animal[policy] = [{"animal": a, "neurons": sum(r["animal"] == a for r in records),
                               **{key: float(np.mean([r["metrics"][policy][key] for r in records if r["animal"] == a]))
                                  for key in ("distance_um", "f1_200")}} for a in animals]
        summary[policy] = {key: float(np.mean([r[key] for r in per_animal[policy]])) for key in ("distance_um", "f1_200")}
    better = lambda left, right: summary[left]["distance_um"] < summary[right]["distance_um"] and summary[left]["f1_200"] > summary[right]["f1_200"]
    conditional_benefit = {policy: better(policy, policy.replace("conditional_", "uniform_"))
                           for policy in POLICIES if policy.startswith("conditional_")}
    eligible = [policy for policy in POLICIES if policy != BASELINE and better(policy, BASELINE) and
                (not policy.startswith("conditional_") or conditional_benefit[policy])]
    chosen = min(eligible, key=lambda p: (summary[p]["distance_um"], -summary[p]["f1_200"], POLICIES.index(p))) if eligible else BASELINE
    selected_joint_benefit = chosen.startswith("conditional_") and conditional_benefit[chosen]
    gate = {"schema": "brm_template_gate_v1", "config": deepcopy(CONFIG), "source_provenance": source_provenance(),
            "selected_policy": chosen, "eligible_nonbaseline_policies": eligible,
            "summary": summary, "animal_metrics": per_animal,
            "conditional_beats_matching_uniform": conditional_benefit,
            "selected_has_joint_specific_validation_gain": bool(selected_joint_benefit),
            "selection_reason": "required_validation_comparisons_improved" if eligible else "no_policy_passed_required_validation_comparisons",
            "validation_records": records, "proof": proof, "validation_records_sha256": _digest(records),
            "proof_sha256": _digest(proof), "unseen_test_gain_guaranteed": False}
    gate["gate_sha256"] = _digest(gate)
    return gate


def predict_gated(query, donors, regional, gate):
    """Apply a reproducible validation decision to an independent test animal."""
    if not isinstance(gate, dict) or gate.get("config") != CONFIG or gate.get("source_provenance") != source_provenance():
        raise ValueError("Gate configuration/source provenance differs.")
    saved = deepcopy(gate)
    digest = saved.pop("gate_sha256", None)
    if digest != _digest(saved) or Geo.canonical(select_policy(gate["validation_records"], gate["proof"])) != Geo.canonical(gate):
        raise ValueError("Gate decision/metrics/proof changed or is not reproducible.")
    proof = gate["proof"]
    if (query.neuron, query.animal, query.full_id) not in _identities(proof["test_ids"]):
        raise ValueError("Query is not a declared independent final-test animal.")
    if Geo.canonical(donor_provenance(donors)) != Geo.canonical(sorted(proof["training_ids"], key=lambda r: r["neuron"])):
        raise ValueError("Exact geometry donor identities/data differ from the declared training set.")
    if proof["state_by_test_animal"][query.animal] != regional.proof["model_state_sha256"]:
        raise ValueError("Query uses an unlisted regional state.")
    if _digest(regional.model) != regional.proof["model_state_sha256"]:
        raise ValueError("Query regional numeric state changed.")
    declared_state = next(s for s in proof["regional_states"] if s["proof"]["model_state_sha256"] == regional.proof["model_state_sha256"])
    if Geo.canonical({"model": regional.model, "proof": regional.proof}) != Geo.canonical(declared_state):
        raise ValueError("Query regional numeric state or complete membership proof changed.")
    output, info = predict_policy(query, donors, regional, gate["selected_policy"])
    info["gate_sha256"] = gate["gate_sha256"]
    info["selected_has_joint_specific_validation_gain"] = gate["selected_has_joint_specific_validation_gain"]
    return output, info
