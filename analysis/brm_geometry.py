"""Exploratory joint-regional-to-geometry decoder; no data I/O on import.

Only a query's own frozen regional S1 state may generate its class posterior
and all donor responsibilities. The geometry decoder learns three log-axis
RMS summaries from other animals; it is not joint geometry EM. The resulting
diagonal deformation retains a measured donor's topology and global route.
No hidden query profile or query arbor is an inference argument.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import importlib.util
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = Path("data/interim/brm_geometry")
PROTOCOL = Path("reports/brm_geometry_protocol.md")
REVIEW = Path("reports/brm_geometry_review.md")
COHORT = Path("data/interim/d111_real")
UPSTREAM = Path("data/interim/brm_real")
CORE_SHA = "3b011643b1b224256fdce1ad8a31014e01e644c6a794bd1d49af5e6543090a5e"
REGIONAL_SHA = "5225ea5833b32fee57e735600eda87a9784cb2b1798cd326ee4203e1b6997648"
JOINT_SHA = "ae7f5ec432bcbd2b0eebb65595549bc12158b1b7c047f5ace167c7840c0c937e"
COHORT_DECL_SHA = "a63b321b007d343aa106c8d0b52326f8a6a42e9b5661f6fe66b6b95e2a462175"
COHORT_FROZEN_SHA = "c4f3867f7323c33473389fade0740eefa9cd1a27649efe44eba3677870e3e8ac"
MODELS = ("joint_deformed", "nearest", "medoid", "tag_nearest", "tag_medoid", "guided_copy",
          "guided_global_deformed", "soma")
COPY_MODELS = ("nearest", "medoid", "tag_nearest", "tag_medoid", "guided_copy")
CONFIG = {"frame": "mirrored-left", "coordinate_system": "CCFv3_AP_DV_ML_um",
          "step_um": 50., "rms_floor_um": 50., "shrink_animals": 1., "scale_bounds": [.5, 2.],
          "soma_mismatch_limit_um": 2., "models": list(MODELS), "samples_per_method": 1,
          "regional_evidence_coordinates": "published_TABLE_soma",
          "nearest_coordinates": "published_TABLE_soma", "geometry_origin": "mirrored_SWC_soma",
          "query_profile_evidence": False, "regional_refitting": False,
          "geometry_features": "ln(max(arc_weighted_axis_RMS_about_SWC_soma,50um))",
          "warp": "positive_diagonal_preserve_all_canonical_vertices_subdivide_long_edges",
          "bootstrap_samples": 10000, "seed": 20261110}
REVIEW_SENTINEL = "BRM_GEOMETRY_REVIEW: PASS"


def _load(name, relative, expected_sha):
    path = ROOT / relative
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected_sha:
        raise RuntimeError(f"Pinned dependency changed: {relative}")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


G = _load("brm_geometry_core", "analysis/d110_geometry.py", CORE_SHA)
R = _load("brm_geometry_regional", "analysis/brm_real.py", REGIONAL_SHA)
sha, canonical, read_json, write_json = R.sha, R.canonical, R.read_json, R.write_json
REGIONS = R.REGIONS


@dataclass(frozen=True)
class Query:
    neuron: str
    animal: str
    full_id: str
    table_soma: np.ndarray
    swc_soma: np.ndarray
    tag: str


@dataclass(frozen=True)
class RegionalState:
    model: dict
    proof: dict


def assert_query_safe(query, regional):
    """Fold membership, not a donor's cached OOF state, defines safety."""
    proof, model = regional.proof, regional.model
    if not query.full_id.startswith(query.neuron + "_") or query.animal != query.neuron.split("_", 1)[0]:
        raise ValueError("Query neuron/full-ID/animal join differs.")
    if hashlib.sha256(canonical(model).encode()).hexdigest() != proof["model_state_sha256"]:
        raise ValueError("Regional numeric state differs from its recorded content hash.")
    if np.asarray(query.table_soma).shape != (3,) or np.asarray(query.swc_soma).shape != (3,) or not np.isfinite(query.table_soma).all() or not np.isfinite(query.swc_soma).all():
        raise ValueError("Both known soma inputs must be finite three-vectors.")
    if query.full_id not in proof["test_ids"] or query.animal not in proof["test_animals"]:
        raise ValueError("Query is not in this regional state's own held-out fold.")
    if query.full_id in proof["train_ids"] or query.animal in proof["train_animals"]:
        raise ValueError("Regional model saw the query animal.")
    if any(i.split("_", 1)[0] == query.animal for i in proof["train_ids"]):
        raise ValueError("A regional training ID belongs to the query animal.")
    if model["training_ids_sha256"] != hashlib.sha256("\n".join(proof["train_ids"]).encode()).hexdigest():
        raise ValueError("Regional state does not bind these exact training IDs.")
    if model["core_sha256"] != JOINT_SHA or model["target_fields"] != list(R.target_fields()):
        raise ValueError("Unexpected regional model/target schema.")
    if len(model["blocks"]) != 1 or set(model["blocks"][0]["fields"]) != {*R.POSITION, R.TAG, *R.target_fields()}:
        raise ValueError("Complete regional model must have one known joint block.")
    b = model["blocks"][0]
    k = len(b["weights"])
    labels = model["field_metadata"][R.TAG]["categories"]
    table = np.asarray(b["tables"][R.TAG])
    if table.shape != (k, len(labels)) or not np.isfinite(table).all() or (table <= 0).any() or not np.allclose(table.sum(axis=1), 1.):
        raise ValueError("Regional categorical probabilities are invalid.")
    for field in (*R.POSITION, *R.target_fields()):
        meta = model["field_metadata"][field]
        means, variances = np.asarray(b["means"][field]), np.asarray(b["variances"][field])
        if not math.isfinite(meta["training_mean"]) or not math.isfinite(meta["training_scale"]) or meta["training_scale"] <= 0 or means.shape != (k,) or variances.shape != (k,) or not np.isfinite(means).all() or not np.isfinite(variances).all() or (variances <= 0).any():
            raise ValueError("Regional continuous parameter shape/finiteness differs.")


def state_posterior(model, coordinates, tags, *, measured_profiles=None):
    """Numeric replay: omit all target factors for queries; include for donors."""
    X = np.asarray(coordinates, dtype=float)
    tags = np.asarray(tags, dtype=object)
    if X.ndim != 2 or X.shape[1] != 3 or tags.shape != (len(X),) or not np.isfinite(X).all():
        raise ValueError("Regional evidence must have finite TABLE coordinates and tags.")
    if measured_profiles is not None:
        Y = np.asarray(measured_profiles, dtype=float)
        if Y.shape != (len(X), len(REGIONS)) or not np.isfinite(Y).all() or (Y < 0).any():
            raise ValueError("Measured DONOR profile shape/finiteness differs.")
    b = model["blocks"][0]
    weights = np.asarray(b["weights"], dtype=float)
    if not np.isfinite(weights).all() or (weights <= 0).any() or not np.isclose(weights.sum(), 1.):
        raise ValueError("Invalid regional mixture weights.")
    logits = np.tile(np.log(weights), (len(X), 1))
    values = {field: X[:, j] for j, field in enumerate(R.POSITION)}
    if measured_profiles is not None:
        values.update({field: Y[:, j] for j, field in enumerate(R.target_fields())})
    for field, observed in values.items():
        meta = model["field_metadata"][field]
        value = (observed - meta["training_mean"]) / meta["training_scale"]
        means, variances = np.asarray(b["means"][field]), np.asarray(b["variances"][field])
        if not np.isfinite(value).all() or (np.abs(value) > R.JOINT.MAX_ABS_VALUE).any() or not np.isfinite(variances).all() or (variances <= 0).any():
            raise ValueError("Invalid regional standardized evidence/variance.")
        logits += -.5 * (math.log(2 * math.pi) + np.log(variances) + (value[:, None] - means)**2 / variances)
    labels = model["field_metadata"][R.TAG]["categories"]
    table = np.asarray(b["tables"][R.TAG])
    for tag in set(tags):
        if tag in labels:  # unknown tags are ignored exactly as frozen BRM
            logits[tags == tag] += np.log(table[:, labels.index(tag)])
    shifted = logits - logits.max(axis=1, keepdims=True)
    posterior = np.exp(shifted - logsumexp(shifted, axis=1, keepdims=True))
    if not np.isfinite(posterior).all() or not np.allclose(posterior.sum(axis=1), 1.):
        raise ValueError("Invalid regional posterior.")
    return posterior


def predicted_profile(model, posterior):
    b = model["blocks"][0]
    means = np.column_stack([model["field_metadata"][f]["training_mean"] +
                             model["field_metadata"][f]["training_scale"] * np.asarray(b["means"][f])
                             for f in R.target_fields()])
    return np.maximum(np.asarray(posterior) @ means, 0.)


def axis_features(arbor):
    arbor.validate()
    weights = arbor.weights()
    relative = arbor.xyz - arbor.xyz[0]
    rms = np.sqrt(np.average(relative**2, axis=0, weights=weights))
    return np.log(np.maximum(rms, CONFIG["rms_floor_um"]))


def fit_decoder(features, animals, responsibilities):
    """Three log-axis moments, animal-balanced and shrunk by one animal."""
    H, r, animals = np.asarray(features, float), np.asarray(responsibilities, float), np.asarray(animals)
    if H.shape != (len(animals), 3) or r.ndim != 2 or r.shape[0] != len(animals) or len(animals) < 1:
        raise ValueError("Decoder training array shapes differ.")
    if not np.isfinite(H).all() or not np.isfinite(r).all() or (r < 0).any() or not np.allclose(r.sum(axis=1), 1.):
        raise ValueError("Invalid feature/responsibility values.")
    counts = Counter(animals)
    a = np.array([1. / counts[animal] for animal in animals])
    A = len(counts)
    global_feature = (a[:, None] * H).sum(axis=0) / A
    mass = (a[:, None] * r).sum(axis=0)
    component_means = ((a[:, None] * r).T @ H + CONFIG["shrink_animals"] * global_feature) / (mass[:, None] + CONFIG["shrink_animals"])
    return {"animal_weights": a, "effective_animals": A, "global_feature": global_feature,
            "component_mass": mass, "component_means": component_means}


def subdivide_preserving_vertices(arbor, *, step=50.):
    """Split each edge only; every supplied vertex is retained exactly."""
    arbor.validate()
    xyz, parents, types = [arbor.xyz[0].copy()], [-1], [int(arbor.types[0])]
    remap = np.empty(len(arbor.xyz), dtype=int)
    remap[0] = 0
    for i in range(1, len(arbor.xyz)):
        old_parent = int(arbor.parent[i])
        start, end = arbor.xyz[old_parent], arbor.xyz[i]
        segments = max(1, int(np.ceil(np.linalg.norm(end - start) / step - 1e-12)))
        previous = int(remap[old_parent])
        for j in range(1, segments):
            xyz.append(start + (end - start) * (j / segments))
            parents.append(previous)
            types.append(2)
            previous = len(xyz) - 1
        xyz.append(end.copy())
        parents.append(previous)
        types.append(int(arbor.types[i]))
        remap[i] = len(xyz) - 1
    output = G.Arbor(np.asarray(xyz), np.asarray(parents), np.asarray(types)).validate()
    if not np.array_equal(output.xyz[remap], arbor.xyz) or not np.array_equal(output.types[remap], arbor.types):
        raise AssertionError("Canonical vertices were not preserved.")
    if len(output.xyz) > 1 and np.linalg.norm(output.xyz[1:] - output.xyz[output.parent[1:]], axis=1).max() > step + 1e-6:
        raise AssertionError("Subdivision left an oversized edge.")
    return output, remap


def deform(arbor, anchor, target_feature, donor_feature):
    raw = np.exp(np.asarray(target_feature) - np.asarray(donor_feature))
    if not np.isfinite(raw).all() or raw.shape != (3,) or (raw <= 0).any():
        raise ValueError("Non-finite/nonpositive deformation factors.")
    clipped = np.clip(raw, *CONFIG["scale_bounds"])
    points = (arbor.xyz - arbor.xyz[0]) * clipped + np.asarray(anchor)
    warped, remap = subdivide_preserving_vertices(G.Arbor(points, arbor.parent.copy(), arbor.types.copy()), step=CONFIG["step_um"])
    return warped, {"raw_scale": raw.tolist(), "scale": clipped.tolist(),
                    "clipped_axes": int(np.count_nonzero(raw != clipped)), "canonical_vertex_map": remap.tolist()}


def predict_query(query, donors, regional):
    """Learn only from donors; Query has no measured profile/arbor slot."""
    assert_query_safe(query, regional)
    donors = sorted(donors, key=lambda d: d["neuron"])
    if not donors or len({d["neuron"] for d in donors}) != len(donors):
        raise ValueError("Distinct geometry donors are required.")
    if any(d["neuron"] == query.neuron or d["animal"] == query.animal or d["full_id"] == query.full_id for d in donors):
        raise ValueError("Query neuron/animal entered geometry decoder training.")
    for d in donors:
        d["arbor"].validate()
        if not np.allclose(d["arbor"].xyz[0], d["swc_soma"], atol=1e-6, rtol=0):
            raise ValueError("Donor arbor uses a different SWC soma origin.")
        if len(d["arbor"].xyz) > 1 and np.linalg.norm(d["arbor"].xyz[1:] - d["arbor"].xyz[d["arbor"].parent[1:]], axis=1).max() > CONFIG["step_um"] + 1e-6:
            raise ValueError("Geometry donors must already be canonical at50um.")
    X = np.asarray([d["table_soma"] for d in donors])
    tags = np.asarray([d["tag"] for d in donors])
    profiles = np.asarray([d["profile"] for d in donors])
    posterior = state_posterior(regional.model, np.asarray(query.table_soma)[None], [query.tag])[0]
    predicted = predicted_profile(regional.model, posterior[None])[0]
    donor_r = state_posterior(regional.model, X, tags, measured_profiles=profiles)
    features = np.asarray([axis_features(d["arbor"]) for d in donors])
    decoder = fit_decoder(features, [d["animal"] for d in donors], donor_r)
    expected_feature = posterior @ decoder["component_means"]
    guided = int(np.argmin(((profiles - predicted)**2).sum(axis=1)))
    nearest = int(np.argmin(((X - query.table_soma)**2).sum(axis=1)))
    medoid = G.training_medoid(donors)
    same = [i for i, d in enumerate(donors) if d["tag"] == query.tag]
    tag_nearest = same[int(np.argmin(((X[same] - query.table_soma)**2).sum(axis=1)))] if same else nearest
    tag_medoid = same[G.training_medoid([donors[i] for i in same])] if same else medoid
    selected = {"nearest": nearest, "medoid": medoid, "tag_nearest": tag_nearest,
                "tag_medoid": tag_medoid, "guided_copy": guided}
    predictions, info = {}, {}
    for name, index in selected.items():
        predictions[name] = donors[index]["arbor"].shifted(query.swc_soma)
        info[name] = {"training_template": donors[index]["neuron"], "same_tag_fallback": bool(not same and name.startswith("tag_"))}
    predictions["joint_deformed"], info["joint_deformed"] = deform(donors[guided]["arbor"], query.swc_soma, expected_feature, features[guided])
    predictions["guided_global_deformed"], info["guided_global_deformed"] = deform(donors[guided]["arbor"], query.swc_soma, decoder["global_feature"], features[guided])
    for name in ("joint_deformed", "guided_global_deformed"):
        info[name]["training_template"] = donors[guided]["neuron"]
    predictions["soma"] = G.Arbor(np.asarray(query.swc_soma).reshape(1, 3), np.array([-1]), np.array([1]))
    info["soma"] = {"training_template": None}
    replay = {"neuron": query.neuron, "animal": query.animal, "full_id": query.full_id,
              "upstream_state_sha256": regional.proof["model_state_sha256"], "upstream_fold": regional.proof["name"],
              "query_table_soma": np.asarray(query.table_soma).tolist(), "query_swc_soma": np.asarray(query.swc_soma).tolist(),
              "query_tag": query.tag, "query_posterior": posterior.tolist(), "query_predicted_profile": predicted.tolist(),
              "donor_ids": [d["neuron"] for d in donors], "donor_full_ids": [d["full_id"] for d in donors],
              "donor_animals": [d["animal"] for d in donors], "donor_table_somas": X.tolist(), "donor_tags": tags.tolist(),
              "donor_profiles": profiles.tolist(), "donor_features": features.tolist(), "donor_responsibilities": donor_r.tolist(),
              "animal_weights": decoder["animal_weights"].tolist(), "effective_animals": decoder["effective_animals"],
              "global_feature": decoder["global_feature"].tolist(), "component_mass": decoder["component_mass"].tolist(),
              "component_means": decoder["component_means"].tolist(), "predicted_feature": expected_feature.tolist(),
              "selected_guided_donor": donors[guided]["neuron"], "prediction_info": info}
    return predictions, info, replay


def cohort_mapping(queries, inputs, upstream_declared, manifest, *, root=ROOT):
    """Metadata and frozen prediction replay only; never open SWC/outcomes."""
    mapped = []
    for source in sorted(queries, key=lambda q: q["neuron"]):
        prefix = source["neuron"] + "_"
        hits = [i for i, full in enumerate(inputs.ids) if full.startswith(prefix)]
        if len(hits) != 1:
            raise ValueError("Brain_neuron must join exactly one published full ID.")
        index = hits[0]
        full_id, tag = str(inputs.ids[index]), str(inputs.tags[index])
        if str(inputs.animal[index]) != source["animal"] or full_id.split("_", 1)[0] != source["animal"]:
            raise ValueError("Joined animal ID differs.")
        canonical_tag = {"Orexin": "Wt(orexin)"}.get(source["tag"], source["tag"])
        if canonical_tag != tag:
            raise ValueError("Known source-tag alias does not reproduce table tag.")
        table_soma, swc_soma = inputs.X[index], np.asarray(source["soma_um"], float)
        difference = float(np.linalg.norm(table_soma - swc_soma))
        if difference > CONFIG["soma_mismatch_limit_um"]:
            raise ValueError("TABLE and mirrored SWC soma discrepancy exceeds two micrometres.")
        folds = [f for f in manifest["folds"] if full_id in f["test_ids"]]
        if len(folds) != 1:
            raise ValueError("Query does not belong to exactly one frozen regional test fold.")
        fold = folds[0]
        model = read_json(Path(root) / UPSTREAM / fold["model_state_path"])
        query = Query(source["neuron"], source["animal"], full_id, table_soma, swc_soma, tag)
        assert_query_safe(query, RegionalState(model, fold))
        posterior = state_posterior(model, table_soma[None], [tag])[0]
        profile = predicted_profile(model, posterior[None])[0]
        with np.load(Path(root) / UPSTREAM / fold["path"], allow_pickle=False) as arrays:
            rows = np.flatnonzero(arrays["ids"] == full_id)
            if len(rows) != 1 or not np.allclose(profile, arrays["joint"][rows[0]], atol=1e-10, rtol=1e-10):
                raise ValueError("TABLE-coordinate query replay differs from its frozen regional prediction.")
        mapped.append({**source, "raw_tag": source["tag"], "tag": tag, "full_id": full_id,
                       "source_row": index, "table_soma_um": table_soma.tolist(), "swc_soma_um": swc_soma.tolist(),
                       "table_minus_swc_um": (table_soma - swc_soma).tolist(), "soma_discrepancy_um": difference,
                       "tag_alias_applied": source["tag"] != tag, "upstream_fold": fold["name"],
                       "upstream_state_path": str(UPSTREAM / fold["model_state_path"]),
                       "upstream_state_sha256": fold["model_state_sha256"],
                       "upstream_train_ids": fold["train_ids"], "upstream_train_animals": fold["train_animals"],
                       "upstream_test_ids": fold["test_ids"], "upstream_test_animals": fold["test_animals"],
                       "upstream_predicted_profile": profile.tolist(), "upstream_posterior": posterior.tolist()})
    return mapped


def source_metadata(root=ROOT):
    root = Path(root)
    if sha(root / COHORT / "declaration.json") != COHORT_DECL_SHA or sha(root / COHORT / "frozen.json") != COHORT_FROZEN_SHA:
        raise ValueError("Original frozen D111 cohort metadata changed.")
    cohort = read_json(root / COHORT / "declaration.json")
    cohort_frozen = read_json(root / COHORT / "frozen.json")
    if cohort_frozen["declaration_sha256"] != sha(root / COHORT / "declaration.json"):
        raise ValueError("D111 frozen declaration hash differs.")
    queries = cohort["queries"]
    if len(queries) != 16 or len({q["neuron"] for q in queries}) != 16 or len({q["animal"] for q in queries}) != 14:
        raise ValueError("The fixed explicit-type2 cohort must contain16neurons/14animals.")
    for q in queries:
        if sha(root / q["file"]) != q["sha256"]:
            raise ValueError("A pinned raw SWC changed.")
    upstream_base = root / UPSTREAM
    frozen, declared = R.verify_frozen(upstream_base, root)
    manifest = read_json(upstream_base / "prediction_manifest.json")
    R.validate_prediction_bundle(upstream_base, manifest, declared, frozen)
    complete, audit = read_json(upstream_base / "RUN_COMPLETE.json"), read_json(upstream_base / "independent_audit.json")
    if audit.get("status") != "PASS" or complete["metrics_sha256"] != sha(upstream_base / "metrics.json") or audit["metrics_sha256"] != complete["metrics_sha256"]:
        raise ValueError("Upstream regional saved-output audit/completion is not valid.")
    if complete["prediction_manifest_sha256"] != sha(upstream_base / "prediction_manifest.json") or audit["prediction_manifest_sha256"] != complete["prediction_manifest_sha256"]:
        raise ValueError("Upstream regional output hashes changed.")
    inputs = R.load_inputs(root)
    mapped = cohort_mapping(queries, inputs, declared, manifest, root=root)
    references = [COHORT / "declaration.json", COHORT / "frozen.json",
                  UPSTREAM / "declaration.json", UPSTREAM / "frozen.json", UPSTREAM / "check.json",
                  UPSTREAM / "prediction_manifest.json", UPSTREAM / "PREDICTIONS_COMPLETE",
                  UPSTREAM / "RUN_COMPLETE.json", UPSTREAM / "metrics.json", UPSTREAM / "independent_audit.json",
                  Path("analysis/d110_geometry.py"), Path("analysis/brm_real.py"), Path("analysis/brm_joint.py"),
                  Path(frozen["review_path"]), Path(R.PROTOCOL)]
    for fold in manifest["folds"]:
        references.extend([UPSTREAM / fold["path"], UPSTREAM / fold["model_state_path"]])
    references.extend(Path(rel) for rel in R.INPUTS.values())
    references.extend(Path(q["file"]) for q in queries)
    hashes = {str(path): sha(root / path) for path in sorted(set(references))}
    return {"queries": mapped, "upstream_sha256": hashes,
            "display_case_ids": sorted(q["neuron"] for q in mapped)[:3],
            "folds": [{"animal": animal, "train": [q["neuron"] for q in mapped if q["animal"] != animal],
                       "test": [q["neuron"] for q in mapped if q["animal"] == animal]}
                      for animal in sorted({q["animal"] for q in mapped})]}


def current_declaration(root=ROOT):
    root = Path(root)
    metadata = source_metadata(root)
    return {"config": CONFIG, **metadata, "runtime": R.runtime(),
            "code_sha256": sha(root / "analysis/brm_geometry.py"), "core_sha256": CORE_SHA,
            "tests_sha256": sha(root / "tests/test_brm_geometry.py"),
            "regional_code_sha256": REGIONAL_SHA, "joint_code_sha256": JOINT_SHA,
            "protocol_sha256": sha(root / PROTOCOL),
            "seeds": {name: R.stable_seed(CONFIG["seed"], "brm_geometry_bootstrap", name)
                      for name in MODELS if name != "joint_deformed"}}


def _open_stage(base):
    if (Path(base) / "RUN_STARTED.json").exists():
        raise RuntimeError("Experiment already started; no metadata rewrite or rerun is allowed.")


def prepare(*, base=None, root=ROOT):
    base = Path(base) if base is not None else Path(root) / OUTPUT
    _open_stage(base)
    if (base / "frozen.json").exists():
        raise RuntimeError("Frozen metadata cannot be rewritten.")
    declared = current_declaration(root)
    write_json(base / "declaration.json", declared, preserve_previous=True)
    return declared


def synthetic_checks():
    # Analytic unequal-neuron animal balance and component shrinkage.
    H = np.log(np.array([[50., 100., 200.], [50., 100., 200.], [200., 400., 800.]]))
    r = np.array([[1., 0.], [1., 0.], [0., 1.]])
    fit = fit_decoder(H, ["A", "A", "B"], r)
    if not np.allclose(fit["global_feature"], (H[0] + H[2]) / 2) or not np.allclose(fit["component_mass"], [1., 1.]):
        raise AssertionError("Animal-balanced decoder identity failed.")
    arbor = G.Arbor(np.array([[0., 0., 0.], [40., 0., 0.], [40., 30., 0.], [40., -30., 0.]]),
                    np.array([-1, 0, 1, 1]), np.array([1, 2, 2, 2]))
    warped, info = deform(arbor, [1., 2., 3.], np.log([400., 25., 50.]), np.log([50., 50., 50.]))
    expected = (arbor.xyz - arbor.xyz[0]) * [2., .5, 1.] + [1., 2., 3.]
    if not np.array_equal(warped.xyz[info["canonical_vertex_map"]], expected) or warped.branches() != arbor.branches():
        raise AssertionError("Warp topology/vertex preservation failed.")
    return {"passed": True, "synthetic_only": True,
            "checks": ["animal_balance_and_unit_animal_shrink", "positive_bounded_diagonal_warp",
                       "canonical_vertex_and_branch_preservation"]}


def check(*, base=None, root=ROOT):
    base = Path(base) if base is not None else Path(root) / OUTPUT
    _open_stage(base)
    if (base / "frozen.json").exists():
        raise RuntimeError("Frozen metadata cannot be rewritten.")
    declared = read_json(base / "declaration.json")
    if declared != current_declaration(root):
        raise ValueError("Declaration is stale.")
    result = {**synthetic_checks(), "declaration_sha256": sha(base / "declaration.json"),
              **{key: declared[key] for key in ("code_sha256", "tests_sha256", "core_sha256", "regional_code_sha256", "joint_code_sha256", "protocol_sha256")}}
    write_json(base / "check.json", result, preserve_previous=True)
    return result


def _review_valid(path, declared):
    text = Path(path).read_text()
    if REVIEW_SENTINEL not in text.splitlines():
        raise ValueError("Independent geometry review PASS is missing.")
    if any(declared[key] not in text for key in ("code_sha256", "tests_sha256", "core_sha256", "regional_code_sha256", "joint_code_sha256", "protocol_sha256")):
        raise ValueError("Independent review does not bind every source/protocol hash.")


def freeze(*, base=None, root=ROOT, review=None):
    root = Path(root)
    base = Path(base) if base is not None else root / OUTPUT
    _open_stage(base)
    declared = read_json(base / "declaration.json")
    if declared != current_declaration(root):
        raise ValueError("Declaration is stale.")
    checked = read_json(base / "check.json")
    if checked.get("passed") is not True or checked.get("declaration_sha256") != sha(base / "declaration.json"):
        raise ValueError("Matching synthetic check is missing.")
    path = root / REVIEW if review is None else Path(review)
    if not path.is_absolute():
        path = root / path
    relative = str(path.resolve().relative_to(root.resolve()))
    _review_valid(path, declared)
    frozen = {"declaration_sha256": sha(base / "declaration.json"), "check_sha256": sha(base / "check.json"),
              "review_path": relative, "review_sha256": sha(path), "config": CONFIG,
              **{key: declared[key] for key in ("code_sha256", "tests_sha256", "core_sha256", "regional_code_sha256", "joint_code_sha256", "protocol_sha256")}}
    write_json(base / "frozen.json", frozen, exclusive=True)
    return frozen


def verify_frozen(base, root=ROOT):
    base, root = Path(base), Path(root)
    frozen, declared = read_json(base / "frozen.json"), read_json(base / "declaration.json")
    if frozen.get("config") != CONFIG or frozen.get("declaration_sha256") != sha(base / "declaration.json") or frozen.get("check_sha256") != sha(base / "check.json"):
        raise ValueError("Frozen metadata/check/config integrity differs.")
    if declared != current_declaration(root):
        raise ValueError("Source/input/protocol/runtime changed since declaration.")
    for key in ("code_sha256", "tests_sha256", "core_sha256", "regional_code_sha256", "joint_code_sha256", "protocol_sha256"):
        if frozen.get(key) != declared[key]:
            raise ValueError("Frozen source provenance differs.")
    path = (root / frozen["review_path"]).resolve()
    path.relative_to(root.resolve())
    if sha(path) != frozen["review_sha256"]:
        raise ValueError("Independent review changed.")
    _review_valid(path, declared)
    checked = read_json(base / "check.json")
    if checked.get("passed") is not True or checked.get("declaration_sha256") != sha(base / "declaration.json"):
        raise ValueError("Frozen synthetic check is invalid.")
    return frozen, declared


def query_from_record(record):
    return Query(record["neuron"], record["animal"], record["full_id"], np.asarray(record["table_soma_um"]),
                 np.asarray(record["swc_soma_um"]), record["tag"])


def regional_for_record(record, root=ROOT):
    path = Path(root) / record["upstream_state_path"]
    if sha(path) != record["upstream_state_sha256"]:
        raise ValueError("Query-safe regional state changed.")
    proof = {"name": record["upstream_fold"], "model_state_sha256": record["upstream_state_sha256"],
             "train_ids": record["upstream_train_ids"], "train_animals": record["upstream_train_animals"],
             "test_ids": record["upstream_test_ids"], "test_animals": record["upstream_test_animals"]}
    result = RegionalState(read_json(path), proof)
    assert_query_safe(query_from_record(record), result)
    return result


def _started(base):
    marker = read_json(Path(base) / "RUN_STARTED.json")
    if marker.get("frozen_sha256") != sha(Path(base) / "frozen.json"):
        raise ValueError("Exclusive start marker does not bind the frozen experiment.")


def load_donors(record, declared, base, *, root=ROOT):
    """Read only other-animal geometry/profiles for this query's decoder."""
    root, base = Path(root), Path(base)
    _started(base)
    others = sorted([q for q in declared["queries"] if q["animal"] != record["animal"]], key=lambda q: q["neuron"])
    path = root / R.INPUTS["outcome"]
    if sha(path) != declared["upstream_sha256"][R.INPUTS["outcome"]]:
        raise ValueError("Measured donor regional source changed before reading.")
    table = R.pd.read_csv(path, index_col=0)
    source_ids = read_json(root / UPSTREAM / "declaration.json")["source"]["ids"]
    if list(table.columns) != list(REGIONS) or table.index.duplicated().any() or set(table.index) != set(source_ids):
        raise ValueError("Donor profile IDs/49-column schema differs.")
    # Query rows stay outside the returned training matrix/API.
    Y = np.log10(1. + table.loc[[q["full_id"] for q in others]].to_numpy(float))
    if not np.isfinite(Y).all() or (Y < 0).any():
        raise ValueError("Invalid measured donor regional values.")
    donors = []
    for i, q in enumerate(others):
        raw_path = root / q["file"]
        if sha(raw_path) != q["sha256"]:
            raise ValueError("Raw donor SWC hash changed before extraction.")
        arbor = G.resample(G.extract_explicit_axon(raw_path, q["reflect_to_left"]), CONFIG["step_um"])
        if not np.allclose(arbor.xyz[0], q["swc_soma_um"], rtol=0, atol=1e-6):
            raise ValueError("Mirrored donor soma differs from frozen metadata.")
        donors.append({"neuron": q["neuron"], "animal": q["animal"], "full_id": q["full_id"],
                       "table_soma": np.asarray(q["table_soma_um"]), "swc_soma": np.asarray(q["swc_soma_um"]),
                       "tag": q["tag"], "arbor": arbor, "profile": Y[i]})
    return donors


def _safe(base, relative):
    result = (Path(base) / relative).resolve()
    result.relative_to(Path(base).resolve())
    return result


def validate_prediction_bundle(base, manifest, declared, frozen):
    base = Path(base)
    if canonical(read_json(base / "prediction_manifest.json")) != canonical(manifest):
        raise ValueError("Passed manifest differs from saved canonical manifest.")
    if canonical(read_json(base / "declaration.json")) != canonical(declared) or canonical(read_json(base / "frozen.json")) != canonical(frozen):
        raise ValueError("Passed declaration/frozen records differ from saved records.")
    if (base / "PREDICTIONS_COMPLETE").read_text().strip() != sha(base / "prediction_manifest.json"):
        raise ValueError("Global prediction completion hash differs.")
    if manifest.get("config") != CONFIG or manifest.get("runtime") != declared["runtime"]:
        raise ValueError("Prediction config/runtime differs.")
    if manifest.get("declaration_sha256") != sha(base / "declaration.json") or manifest.get("frozen_sha256") != sha(base / "frozen.json"):
        raise ValueError("Prediction metadata provenance differs.")
    for key in ("code_sha256", "tests_sha256", "core_sha256", "regional_code_sha256", "joint_code_sha256", "protocol_sha256"):
        if manifest.get(key) != declared[key]:
            raise ValueError("Prediction source hash differs.")
    expected = {q["neuron"]: q for q in declared["queries"]}
    records = manifest.get("queries", [])
    if len(records) != len(expected) or len({q.get("neuron") for q in records}) != len(expected) or {q.get("neuron") for q in records} != set(expected):
        raise ValueError("Wrong, missing or duplicate prediction queries.")
    seen_paths = set()
    for q in records:
        source = expected[q["neuron"]]
        for key in ("animal", "full_id", "tag", "table_soma_um", "swc_soma_um", "upstream_fold", "upstream_state_path", "upstream_state_sha256"):
            if q.get(key) != source[key]:
                raise ValueError("Prediction query/state metadata differs.")
        donors = sorted([r for r in declared["queries"] if r["animal"] != q["animal"]], key=lambda r: r["neuron"])
        train_ids, train_animals = [d["neuron"] for d in donors], sorted({d["animal"] for d in donors})
        if q.get("train_ids") != train_ids or q.get("train_animals") != train_animals:
            raise ValueError("Whole-animal geometry source exclusion differs.")
        if q.get("decoder_path") != f"decoders/{q['neuron']}.json":
            raise ValueError("Unexpected decoder artifact path.")
        decoder_path = _safe(base, q["decoder_path"])
        if sha(decoder_path) != q["decoder_sha256"]:
            raise ValueError("Decoder numeric-state hash differs.")
        state = read_json(decoder_path)
        if state.get("neuron") != q["neuron"] or state.get("animal") != q["animal"] or state.get("donor_ids") != train_ids or state.get("donor_animals") != [d["animal"] for d in donors] or state.get("donor_full_ids") != [d["full_id"] for d in donors]:
            raise ValueError("Decoder query/donor provenance differs.")
        if state.get("upstream_state_sha256") != source["upstream_state_sha256"] or state.get("upstream_fold") != source["upstream_fold"]:
            raise ValueError("Decoder did not use its query-safe regional state.")
        if not np.allclose(state["query_predicted_profile"], source["upstream_predicted_profile"], atol=1e-10, rtol=1e-10) or not np.allclose(state["query_posterior"], source["upstream_posterior"], atol=1e-10, rtol=1e-10):
            raise ValueError("Query regional prediction/posterior parity differs.")
        rows = q.get("predictions", [])
        if len(rows) != len(MODELS) or len({p.get("model") for p in rows}) != len(MODELS) or {p.get("model") for p in rows} != set(MODELS):
            raise ValueError("All eight deterministic methods are required exactly once.")
        for p in rows:
            if p.get("sample") != 0 or p.get("seed") is not None or p.get("path") != f"predictions/{q['neuron']}/{p['model']}_0.npz":
                raise ValueError("Unexpected prediction sample/seed/path.")
            path = _safe(base, p["path"])
            if path in seen_paths or sha(path) != p["sha256"]:
                raise ValueError("Duplicate/changed geometry prediction artifact.")
            seen_paths.add(path)
            template = p.get("training_template")
            if p["model"] == "soma":
                if template is not None:
                    raise ValueError("Soma baseline must not have a donor.")
            elif template not in train_ids:
                raise ValueError("A geometry donor was missing or from the query animal.")
            if p["model"] in ("joint_deformed", "guided_global_deformed", "guided_copy") and template != state["selected_guided_donor"]:
                raise ValueError("Guided methods must share the identical selected donor.")
            for key, value in state["prediction_info"][p["model"]].items():
                if p.get(key) != value:
                    raise ValueError("Prediction deformation/donor metadata differs from decoder state.")
            arbor = G.load_arbor(path)
            if not np.allclose(arbor.xyz[0], source["swc_soma_um"], atol=1e-6, rtol=0):
                raise ValueError("Prediction is not anchored at the known query SWC soma.")
            if len(arbor.xyz) != p.get("nodes") or arbor.branches() != p.get("branches") or not math.isclose(arbor.cable(), p["cable_um"], rel_tol=1e-12, abs_tol=1e-8):
                raise ValueError("Prediction geometry metadata differs from saved arrays.")
            if len(arbor.xyz) > 1 and np.linalg.norm(arbor.xyz[1:] - arbor.xyz[arbor.parent[1:]], axis=1).max() > CONFIG["step_um"] + 1e-6:
                raise ValueError("Prediction has an oversized scoring edge.")


def aggregate(queries, declared):
    animals = []
    for animal in sorted({q["animal"] for q in queries}):
        qs = [q for q in queries if q["animal"] == animal]
        animals.append({"animal": animal, "neurons": [q["neuron"] for q in qs],
                        "models": {m: G.mean_metrics([q["models"][m]["mean"] for q in qs]) for m in MODELS}})
    summary = {m: G.mean_metrics([a["models"][m] for a in animals]) for m in MODELS}
    comparison = {}
    for name in MODELS:
        if name == "joint_deformed":
            continue
        distance = np.array([a["models"]["joint_deformed"]["distance_um"] - a["models"][name]["distance_um"] for a in animals])
        f1 = np.array([a["models"]["joint_deformed"]["f1_200"] - a["models"][name]["f1_200"] for a in animals])
        rng = np.random.default_rng(declared["seeds"][name])
        indices = rng.integers(0, len(animals), (CONFIG["bootstrap_samples"], len(animals)))
        comparison[name] = {"distance_delta_um": float(distance.mean()), "f1_delta": float(f1.mean()),
                            "distance_bootstrap95_um": np.quantile(distance[indices].mean(axis=1), [.025, .975]).tolist(),
                            "f1_bootstrap95": np.quantile(f1[indices].mean(axis=1), [.025, .975]).tolist(),
                            "animals_distance_better": int((distance < 0).sum()), "animals_f1_better": int((f1 > 0).sum()),
                            "animals_both_better": int(((distance < 0) & (f1 > 0)).sum()),
                            "animal_differences": [{"animal": a["animal"], "distance_delta_um": float(d), "f1_delta": float(f)}
                                                   for a, d, f in zip(animals, distance, f1)], "descriptive_only": True}
    passes = lambda name: summary["joint_deformed"]["distance_um"] < summary[name]["distance_um"] and summary["joint_deformed"]["f1_200"] > summary[name]["f1_200"]
    copies = all(passes(name) for name in COPY_MODELS)
    global_deformation = passes("guided_global_deformed")
    return animals, summary, {"against": comparison, "beats_all_copying": copies,
                               "beats_global_deformation": global_deformation,
                               "conditional_deformation_useful": copies and global_deformation}


def score_predictions(base, manifest, frozen, declared, *, root=ROOT):
    base, root = Path(base), Path(root)
    validate_prediction_bundle(base, manifest, declared, frozen)
    _started(base)
    expected = {q["neuron"]: q for q in declared["queries"]}
    records = []
    for prediction in manifest["queries"]:
        q = expected[prediction["neuron"]]
        raw_path = root / q["file"]
        if sha(raw_path) != q["sha256"]:
            raise ValueError("Raw query SWC changed immediately before truth extraction.")
        truth = G.resample(G.extract_explicit_axon(raw_path, q["reflect_to_left"]), CONFIG["step_um"])
        if not np.allclose(truth.xyz[0], q["swc_soma_um"], atol=1e-6, rtol=0):
            raise ValueError("Scorer truth soma differs from known query anchor.")
        truth_relative = f"truth/{q['neuron']}.npz"
        truth_path = base / truth_relative
        truth_path.parent.mkdir(exist_ok=True)
        with truth_path.open("xb") as f:
            np.savez_compressed(f, xyz_um=truth.xyz, parent=truth.parent, types=truth.types)
        models = {}
        for p in prediction["predictions"]:
            path = _safe(base, p["path"])
            if sha(path) != p["sha256"]:
                raise ValueError("Saved prediction changed after the global barrier.")
            metrics = G.geometry_metrics(G.load_arbor(path), truth)
            models[p["model"]] = {"samples": [metrics], "mean": metrics}
        records.append({"neuron": q["neuron"], "animal": q["animal"], "full_id": q["full_id"], "tag": q["tag"],
                        "table_soma_um": q["table_soma_um"], "swc_soma_um": q["swc_soma_um"],
                        "truth_path": truth_relative, "truth_sha256": sha(truth_path), "models": models})
    animal_metrics, summary, comparison = aggregate(records, declared)
    result = {"config": CONFIG, "prediction_manifest_sha256": sha(base / "prediction_manifest.json"),
              "queries": records, "animal_metrics": animal_metrics, "summary": summary, "comparison": comparison,
              "deformation": {name: {"clipped_axes": sum(p.get("clipped_axes", 0) for q in manifest["queries"] for p in q["predictions"] if p["model"] == name),
                                      "queries_with_clipping": sum(any(p.get("clipped_axes", 0) > 0 for p in q["predictions"] if p["model"] == name) for q in manifest["queries"])}
                              for name in ("joint_deformed", "guided_global_deformed")},
              "limitations": ["Exploratory reused16neurons/14animals, not independent confirmation.",
                              "Topology and global route are inherited from a measured donor; only three axis scales are learned.",
                              "Regional states are reused without refitting; query-safe state drives query and all donor responsibilities.",
                              "No synapses, anatomical obstacles, cell-specific expression, or calibrated biological uncertainty."]}
    write_json(base / "metrics.json", result, exclusive=True)
    return result


def run(*, base=None, root=ROOT, progress=print):
    root = Path(root)
    base = Path(base) if base is not None else root / OUTPUT
    _open_stage(base)
    frozen, declared = verify_frozen(base, root)
    start = time.monotonic()
    write_json(base / "RUN_STARTED.json", {"frozen_sha256": sha(base / "frozen.json"),
               "started_utc": datetime.now(timezone.utc).isoformat()}, exclusive=True)
    try:
        manifest = {"config": CONFIG, "runtime": declared["runtime"],
                    "declaration_sha256": sha(base / "declaration.json"), "frozen_sha256": sha(base / "frozen.json"),
                    **{key: declared[key] for key in ("code_sha256", "tests_sha256", "core_sha256", "regional_code_sha256", "joint_code_sha256", "protocol_sha256")}, "queries": []}
        for q in declared["queries"]:
            donors = load_donors(q, declared, base, root=root)
            regional = regional_for_record(q, root)
            predictions, info, decoder = predict_query(query_from_record(q), donors, regional)
            if not np.allclose(decoder["query_predicted_profile"], q["upstream_predicted_profile"], atol=1e-10, rtol=1e-10) or not np.allclose(decoder["query_posterior"], q["upstream_posterior"], atol=1e-10, rtol=1e-10):
                raise ValueError("Query TABLE posterior/profile replay differs from frozen regional prediction.")
            decoder_relative = f"decoders/{q['neuron']}.json"
            write_json(base / decoder_relative, decoder, exclusive=True)
            record = {key: q[key] for key in ("neuron", "animal", "full_id", "tag", "table_soma_um", "swc_soma_um", "upstream_fold", "upstream_state_path", "upstream_state_sha256")}
            record.update({"train_ids": sorted(d["neuron"] for d in donors), "train_animals": sorted({d["animal"] for d in donors}),
                           "decoder_path": decoder_relative, "decoder_sha256": sha(base / decoder_relative), "predictions": []})
            for name in MODELS:
                arbor = predictions[name]
                relative = f"predictions/{q['neuron']}/{name}_0.npz"
                path = base / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("xb") as f:
                    np.savez_compressed(f, xyz_um=arbor.xyz, parent=arbor.parent, types=arbor.types)
                record["predictions"].append({"model": name, "sample": 0, "seed": None, "path": relative, "sha256": sha(path),
                                              "nodes": len(arbor.xyz), "cable_um": arbor.cable(), "branches": arbor.branches(), **info[name]})
            manifest["queries"].append(record)
            write_json(base / "prediction_manifest.partial.json", manifest)
            if progress:
                progress(f"{q['neuron']}: all8 prediction trees saved")
        write_json(base / "prediction_manifest.json", manifest, exclusive=True)
        with (base / "PREDICTIONS_COMPLETE").open("x") as f:
            f.write(sha(base / "prediction_manifest.json") + "\n")
        result = score_predictions(base, manifest, frozen, declared, root=root)
        write_json(base / "RUN_COMPLETE.json", {"elapsed_seconds": time.monotonic() - start,
                   "prediction_manifest_sha256": sha(base / "prediction_manifest.json"), "metrics_sha256": sha(base / "metrics.json")}, exclusive=True)
        return result
    except Exception as error:
        write_json(base / "RUN_FAILED.json", {"error_type": type(error).__name__, "error": str(error),
                   "elapsed_seconds": time.monotonic() - start}, exclusive=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "check", "freeze", "run"))
    parser.add_argument("--review", default=None)
    args = parser.parse_args()
    if args.stage == "prepare":
        result = prepare()
        print(canonical({"metadata_only": True, "queries": len(result["queries"]), "declaration_sha256": sha(ROOT / OUTPUT / "declaration.json")}))
    elif args.stage == "check":
        print(canonical(check()))
    elif args.stage == "freeze":
        print(canonical(freeze(review=args.review)))
    else:
        result = run(progress=lambda text: print(text, flush=True))
        print(canonical({"summary": result["summary"], "comparison": result["comparison"]}))


if __name__ == "__main__":
    main()
