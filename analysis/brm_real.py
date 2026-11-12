"""Frozen, one-run regional-profile BRM comparison; no data I/O on import.

prepare/check/freeze are metadata/synthetic only. Only guarded run opens the
regional outcome table. This is prediction of 49 published regional values,
not 3D trajectories, synapses, gene expression, causation, or calibrated
Gaussian physical uncertainty. D109 code is never imported or modified.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import sys
import time
import zlib

for _name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd
import scipy
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = Path("data/interim/brm_real")
PROTOCOL = Path("reports/brm_real_protocol.md")
REVIEW = Path("reports/brm_real_review.md")
CORE = Path("analysis/brm_joint.py")
CORE_SHA = "ae7f5ec432bcbd2b0eebb65595549bc12158b1b7c047f5ace167c7840c0c937e"
INPUTS = {"soma": "data/raw/jiao2025/All_soma_infor_HY_neurons.csv",
          "labels": "data/raw/jiao2025/geno_location_labels.csv",
          "outcome": "data/raw/jiao2025/all_neurons_projection_length_fig1.csv",
          "splits": "data/interim/d109_splits.json"}
INPUT_HASHES = {"soma": "cc01d048bf5a7c9c9f46869fb750f2d67a98ca934903d2e2e71c1a67b1d97cb9",
                "labels": "8260a279767f4d0d2e28b7b426061f1fa0ea17158b65109537738ffbf784c3da",
                "outcome": "2072bc20550607c195be0a42cee9fcbe0e0ab5aaebfaee0a1016b2a43d35e63c",
                "splits": "606dafac82774a77dc873f3aad65a986b21b5f9d55eb5142d0d865f89d972947"}
REGIONS = tuple(("ARH PMv TU TMv MM AHN VMH DMH PVH MPO MPN PH LHA LPO ZI ME MRN PAG SCm SNr VTA DR MY IO NTS PARN IRN "
                 "PB PRNr PCG PA BMA BLA MS BST SI MEA CEA LS ACB LH PVT MD PT AM SPFm HPF Isocortex OLF").split())
POSITION = ("soma_AP", "soma_DV", "soma_ML")
TAG = "source_tag"
MODELS = ("joint", "position_only", "tag_only", "independent_k1", "global_mean", "tag_mean",
          "nearest", "medoid", "tag_nearest", "tag_medoid")
COPY_MODELS = ("nearest", "medoid", "tag_nearest", "tag_medoid")
AVERAGES = ("independent_k1", "global_mean", "tag_mean")
CONFIG = {"seed": 20261110, "k_grid": [1, 2, 4, 8], "inner_folds": 3, "outer_folds": 5,
          "joint_budget": {"n_init": 2, "max_iter": 100, "tol": 1e-5},
          "regions": list(REGIONS), "models": list(MODELS), "target": "log10(1+published_length)",
          "query_inputs": [*POSITION, TAG], "all_targets_hidden_together": True,
          "selection": "full_input_mean_animal_MSE_only_ties_smallest_K",
          "ablations": "same_selected_fit_not_separately_tuned", "nearest_metric": "raw_published_soma_Euclidean",
          "medoid_metric": "exact_animal_weighted_squared_log_Euclidean",
          "bootstrap_samples": 10000, "gaussian_uncertainty_claim": False}
REVIEW_SENTINEL = "BRM_REAL_REVIEW: PASS"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value, *, exclusive=False, preserve_previous=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = canonical(value)
    if preserve_previous and path.exists() and path.read_text() != data:
        archive = path.parent / "archives" / (path.name + "." + sha(path))
        archive.parent.mkdir(exist_ok=True)
        if not archive.exists():
            archive.write_bytes(path.read_bytes())
    with path.open("x" if exclusive else "w") as f:
        f.write(data)


def stable_seed(*parts):
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:4], "big")


def runtime():
    return {"python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "pandas": pd.__version__}


def _core_module():
    path = ROOT / CORE
    if sha(path) != CORE_SHA:
        raise RuntimeError("Unchanged reviewed BRM core hash does not match.")
    spec = importlib.util.spec_from_file_location("brm_joint_real_core", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


JOINT = _core_module()


@dataclass
class Inputs:
    ids: np.ndarray
    animal: np.ndarray
    X: np.ndarray
    tags: np.ndarray

    def subset(self, indices):
        return Inputs(self.ids[indices], self.animal[indices], self.X[indices], self.tags[indices])


def validate_inputs(data):
    n = len(data.ids)
    if n < 1 or len(set(data.ids)) != n or data.X.shape != (n, 3):
        raise ValueError("Unique IDs and three soma columns are required.")
    if data.animal.shape != (n,) or data.tags.shape != (n,) or not np.isfinite(data.X).all():
        raise ValueError("Input shape/finiteness mismatch.")
    if any(not isinstance(x, (str, np.str_)) or not x for values in (data.ids, data.animal, data.tags) for x in values):
        raise ValueError("Identity, animal, and tag metadata must be nonempty strings.")
    animal_tags = {}
    for animal, tag in zip(data.animal, data.tags):
        if animal_tags.setdefault(str(animal), str(tag)) != str(tag):
            raise ValueError("An animal carries more than one source tag.")


def load_inputs(root=ROOT, *, expected=None):
    """Allowlisted input columns only; no outcome-derived type/class fields."""
    root = Path(root)
    s = pd.read_csv(root / INPUTS["soma"], usecols=["file", "geno", "AP", "DV", "ML"],
                    dtype={"file": str, "geno": str})
    g = pd.read_csv(root / INPUTS["labels"], usecols=[0, 1], dtype=str)
    ids = s["file"].to_numpy(str)
    if s["file"].isna().any() or s["geno"].isna().any():
        raise ValueError("Missing source ID/tag.")
    if len(g) != len(ids) or g.iloc[:, 0].duplicated().any() or set(g.iloc[:, 0]) != set(ids):
        raise ValueError("Label and soma IDs do not join one-to-one exactly.")
    label_tag = dict(zip(g.iloc[:, 0], g.iloc[:, 1]))
    tags = s["geno"].to_numpy(str)
    if any(label_tag[i] != tag for i, tag in zip(ids, tags)):
        raise ValueError("Source tag differs between joined tables.")
    parts = [i.split("_", 3) for i in ids]
    if any(len(p) != 4 or p[2] != t for p, t in zip(parts, tags)):
        raise ValueError("Published ID does not agree with source tag.")
    data = Inputs(ids, np.asarray([p[0] for p in parts]), s[["AP", "DV", "ML"]].to_numpy(float), tags)
    validate_inputs(data)
    if expected is None:
        expected = {"neurons": 7180, "animals": 166, "tags": 16, "all_zero": 13}
    if len(ids) != expected["neurons"] or len(set(data.animal)) != expected["animals"] or len(set(tags)) != expected["tags"]:
        raise ValueError("Input cohort count differs from declaration.")
    return data


def validate_folds(data, folds, *, count=None):
    validate_inputs(data)
    n = len(data.ids)
    if count is not None and len(folds) != count:
        raise ValueError("Wrong number of folds.")
    names, tests = set(), []
    universe = set(range(n))
    for f in folds:
        if not isinstance(f.get("name"), str) or f["name"] in names:
            raise ValueError("Fold names must be unique strings.")
        names.add(f["name"])
        tr, te = f["train"], f["test"]
        if not tr or not te or any(type(i) is not int or i < 0 or i >= n for i in tr + te):
            raise ValueError("Invalid fold indices.")
        if len(set(tr)) != len(tr) or len(set(te)) != len(te) or set(tr) & set(te) or set(tr) | set(te) != universe:
            raise ValueError("Fold is not an exact disjoint row partition.")
        if set(data.animal[tr]) & set(data.animal[te]):
            raise ValueError("An animal crosses a fold boundary.")
        tests.extend(te)
    if sorted(tests) != list(range(n)):
        raise ValueError("Each row must be tested exactly once.")


def inner_folds(data, outer_indices, fold_name):
    """Own deterministic metadata-only tag-stratified whole-animal dealing."""
    indices = np.asarray(outer_indices, dtype=int)
    sub = data.subset(indices)
    validate_inputs(sub)
    animals = sorted(set(sub.animal))
    if len(animals) < CONFIG["inner_folds"]:
        raise ValueError("Too few training animals for inner folds.")
    tag_of = dict(zip(sub.animal, sub.tags))
    assigned, offset = {}, 0
    for tag in sorted(set(sub.tags)):
        group = sorted(a for a in animals if tag_of[a] == tag)
        key = "|".join(map(str, ("deal", "inner|S1|" + fold_name, tag)))
        rng = np.random.default_rng([CONFIG["seed"], zlib.crc32(key.encode())])
        for j in rng.permutation(len(group)):
            assigned[group[j]] = offset % CONFIG["inner_folds"]
            offset += 1
    memberships = np.asarray([assigned[a] for a in sub.animal])
    local = [{"name": f"inner{i}", "train": np.flatnonzero(memberships != i).tolist(),
              "test": np.flatnonzero(memberships == i).tolist()} for i in range(CONFIG["inner_folds"])]
    validate_folds(sub, local, count=CONFIG["inner_folds"])
    return [{"name": f["name"], "train": indices[f["train"]].tolist(), "test": indices[f["test"]].tolist()}
            for f in local]


def source_metadata(root=ROOT, *, expected=None, pinned_inputs=None):
    root = Path(root)
    pinned_inputs = INPUT_HASHES if pinned_inputs is None else pinned_inputs
    actual = {key: sha(root / rel) for key, rel in INPUTS.items()}
    if actual != pinned_inputs:
        raise ValueError("Raw input or immutable D109 split hash changed.")
    data = load_inputs(root, expected=expected)
    source = read_json(root / INPUTS["splits"])
    ids_hash = hashlib.sha256("\n".join(data.ids).encode()).hexdigest()
    if source.get("ids_sha256") != ids_hash or source.get("n") != len(data.ids):
        raise ValueError("D109 indices do not bind this exact soma CSV order.")
    if source.get("inputs", {}).get("soma") != actual["soma"] or source.get("inputs", {}).get("labels") != actual["labels"]:
        raise ValueError("D109 splits were declared for other inputs.")
    outer = source["S1"]
    validate_folds(data, outer, count=CONFIG["outer_folds"])
    if any(not set(data.tags[f["test"]]) <= set(data.tags[f["train"]]) for f in outer):
        raise ValueError("A represented-line outer fold has an unseen test tag.")
    inner = {f["name"]: inner_folds(data, f["train"], f["name"]) for f in outer}
    seeds = {"inner": {f["name"]: {str(k): {g["name"]: stable_seed(CONFIG["seed"], f["name"], "inner", g["name"], k)
                         for g in inner[f["name"]]} for k in CONFIG["k_grid"]} for f in outer},
             "outer": {f["name"]: {str(k): stable_seed(CONFIG["seed"], f["name"], "outer", k)
                                    for k in CONFIG["k_grid"]} for f in outer},
             "bootstrap": {name: stable_seed(CONFIG["seed"], "brm_real_bootstrap", name) for name in MODELS if name != "joint"}}
    return data, {"input_sha256": actual, "ids_sha256": ids_hash,
                  "ids": data.ids.tolist(), "animals": data.animal.tolist(), "tags": data.tags.tolist(),
                  "outer_folds": outer, "inner_folds": inner, "seeds": seeds}


def current_declaration(root=ROOT, *, expected=None, pinned_inputs=None):
    root = Path(root)
    _, metadata = source_metadata(root, expected=expected, pinned_inputs=pinned_inputs)
    if sha(root / CORE) != CORE_SHA:
        raise ValueError("Frozen BRM core differs from reviewed implementation.")
    return {"schema": "brm_real_v1", "config": CONFIG, "runtime": runtime(), "source": metadata,
            "code_sha256": sha(root / "analysis/brm_real.py"), "core_sha256": sha(root / CORE),
            "protocol_sha256": sha(root / PROTOCOL),
            "expected": expected or {"neurons": 7180, "animals": 166, "tags": 16, "all_zero": 13}}


def _stage_open(base):
    base = Path(base)
    if (base / "RUN_STARTED.json").exists():
        raise RuntimeError("This experiment has started; metadata cannot be changed or run repeated.")


def prepare(*, base=None, root=ROOT, expected=None, pinned_inputs=None):
    base = Path(base) if base is not None else Path(root) / OUTPUT
    _stage_open(base)
    if (base / "frozen.json").exists():
        raise RuntimeError("Frozen metadata cannot be regenerated.")
    declaration = current_declaration(root, expected=expected, pinned_inputs=pinned_inputs)
    write_json(base / "declaration.json", declaration, preserve_previous=True)
    return declaration


def target_fields(regions=REGIONS):
    return tuple("region_" + r for r in regions)


def fit_joint(train, Y, k, seed, *, regions=REGIONS):
    validate_inputs(train)
    Y = np.asarray(Y, dtype=float)
    if Y.shape != (len(train.ids), len(regions)) or not np.isfinite(Y).all() or (Y < 0).any():
        raise ValueError("Training target matrix must be finite, nonnegative, and training-only.")
    ys = target_fields(regions)
    types = {f: "continuous" for f in (*POSITION, *ys)}
    types[TAG] = "categorical"
    rows = []
    for i in range(len(train.ids)):
        row = {f: float(train.X[i, j]) for j, f in enumerate(POSITION)}
        row[TAG] = str(train.tags[i])
        row.update({f: float(Y[i, j]) for j, f in enumerate(ys)})
        rows.append(row)
    return JOINT.JointModel(types, max_components=k, seed=seed, **CONFIG["joint_budget"]).fit(rows)


def batch_predict(model, query_X, query_tags, *, mode="full", regions=REGIONS):
    """Exact core posterior/means, vectorized; interface has NO query targets."""
    if mode not in ("full", "position", "tag", "none"):
        raise ValueError("Unknown evidence mode.")
    X = np.asarray(query_X, dtype=float)
    tags = np.asarray(query_tags, dtype=object)
    if X.ndim != 2 or X.shape[1] != 3 or tags.shape != (len(X),) or np.isinf(X).any():
        raise ValueError("Query evidence shapes/finite bounds differ.")
    ys = target_fields(regions)
    model._require_fitted()
    prediction = np.empty((len(X), len(ys)))
    done = set()
    for b in model.blocks_:
        logjoint = np.tile(np.log(b.params.weights), (len(X), 1))
        for j, field in enumerate(POSITION):
            if mode not in ("full", "position") or field not in b.fields:
                continue
            meta = model.field_metadata_[field]
            observed = np.isfinite(X[:, j])
            standardized = (X[observed, j] - meta["training_mean"]) / meta["training_scale"]
            if not np.isfinite(standardized).all() or (np.abs(standardized) > JOINT.MAX_ABS_VALUE).any():
                raise ValueError("Query coordinate exceeds core numerical bounds.")
            logjoint[observed] += -.5 * (math.log(2 * math.pi) + np.log(b.params.variances[field]) +
                                        (standardized[:, None] - b.params.means[field])**2 / b.params.variances[field])
        if mode in ("full", "tag") and TAG in b.fields:
            categories = model._category_indices[TAG]
            for raw in set(t for t in tags if not JOINT._missing(t)):
                key = JOINT._category(raw)
                if key in categories:  # unknown category is ignored, exactly as core.infer
                    logjoint[tags == raw] += np.log(b.params.tables[TAG][:, categories[key]])
        shifted = logjoint - logjoint.max(axis=1, keepdims=True)
        weights = np.exp(shifted - logsumexp(shifted, axis=1, keepdims=True))
        for j, field in enumerate(ys):
            if field not in b.fields:
                continue
            meta = model.field_metadata_[field]
            means = meta["training_mean"] + meta["training_scale"] * b.params.means[field]
            prediction[:, j] = weights @ means
            done.add(field)
    if done != set(ys) or not np.isfinite(prediction).all() or (prediction < -1e-10).any():
        raise ValueError("Unsupported or invalid predicted targets.")
    return np.maximum(prediction, 0.)


def animal_weights(animals):
    counts = Counter(animals)
    weights = np.asarray([1. / counts[a] for a in animals])
    return weights / weights.sum()


def baseline_predictions(train, Y, query, *, return_sources=False):
    """Training profiles + query evidence only; exact weighted medoid in O(NR)."""
    validate_inputs(train)
    validate_inputs(query)
    if set(train.ids) & set(query.ids) or set(train.animal) & set(query.animal):
        raise ValueError("Training/query IDs or whole animals overlap.")
    Y = np.asarray(Y, dtype=float)
    if Y.ndim != 2 or Y.shape[0] != len(train.ids) or not np.isfinite(Y).all() or (Y < 0).any():
        raise ValueError("Baseline targets are invalid.")
    weights = animal_weights(train.animal)
    mean = weights @ Y
    # For squared Euclidean loss: argmin_i sum_j w_j ||y_i-y_j||²
    # equals argmin_i ||y_i-sum_j w_j y_j||². Ties use source CSV order.
    medoid = int(np.argmin(((Y - mean)**2).sum(axis=1)))
    nearest = np.empty(len(query.ids), dtype=int)
    for start in range(0, len(query.ids), 128):
        distances = ((query.X[start:start + 128, None] - train.X[None])**2).sum(axis=2)
        nearest[start:start + 128] = np.argmin(distances, axis=1)
    tag_nearest, tag_medoid = nearest.copy(), np.full(len(query.ids), medoid, dtype=int)
    tag_mean = np.repeat(mean[None], len(query.ids), axis=0)
    fallback = np.zeros(len(query.ids), dtype=bool)
    for tag in sorted(set(query.tags)):
        rows, known = np.flatnonzero(query.tags == tag), np.flatnonzero(train.tags == tag)
        if len(known) == 0:
            fallback[rows] = True
            continue
        mu = animal_weights(train.animal[known]) @ Y[known]
        tag_mean[rows] = mu
        tag_medoid[rows] = known[int(np.argmin(((Y[known] - mu)**2).sum(axis=1)))]
        for row in rows:
            tag_nearest[row] = known[int(np.argmin(((train.X[known] - query.X[row])**2).sum(axis=1)))]
    predictions = {"independent_k1": np.repeat(Y.mean(axis=0)[None], len(query.ids), axis=0),
                   "global_mean": np.repeat(mean[None], len(query.ids), axis=0), "tag_mean": tag_mean,
                   "nearest": Y[nearest], "medoid": np.repeat(Y[medoid][None], len(query.ids), axis=0),
                   "tag_nearest": Y[tag_nearest], "tag_medoid": Y[tag_medoid]}
    sources = {"nearest": train.ids[nearest].tolist(), "medoid": [str(train.ids[medoid])] * len(query.ids),
               "tag_nearest": train.ids[tag_nearest].tolist(), "tag_medoid": train.ids[tag_medoid].tolist(),
               "same_tag_fallback_rows": np.flatnonzero(fallback).tolist()}
    return (predictions, sources) if return_sources else predictions


def mean_animal_mse(truth, prediction, animals):
    error = ((np.asarray(truth) - np.asarray(prediction))**2).mean(axis=1)
    return float(np.mean([error[np.asarray(animals) == a].mean() for a in sorted(set(animals))]))


def fit_outer_fold(train, Y_train, query, inner, *, fold_name, regions=REGIONS, progress=None):
    """No query truth argument; inner indices are LOCAL to outer training only."""
    validate_inputs(train)
    validate_inputs(query)
    if set(train.ids) & set(query.ids) or set(train.animal) & set(query.animal):
        raise ValueError("Outer animal/ID leakage.")
    validate_folds(train, inner, count=CONFIG["inner_folds"])
    Y = np.asarray(Y_train, dtype=float)
    selections = []
    for k in CONFIG["k_grid"]:
        animal_errors, fit_records = {}, []
        for split in inner:
            tr, te = np.asarray(split["train"]), np.asarray(split["test"])
            model = fit_joint(train.subset(tr), Y[tr], k,
                              stable_seed(CONFIG["seed"], fold_name, "inner", split["name"], k), regions=regions)
            prediction = batch_predict(model, train.X[te], train.tags[te], regions=regions)
            error = ((Y[te] - prediction)**2).mean(axis=1)
            for animal in sorted(set(train.animal[te])):
                if str(animal) in animal_errors:
                    raise ValueError("Inner animal appears in multiple validation folds.")
                animal_errors[str(animal)] = float(error[train.animal[te] == animal].mean())
            fit_records.append({"inner": split["name"], "fit_seed": model.seed,
                                "train_ids_sha256": hashlib.sha256("\n".join(train.ids[tr]).encode()).hexdigest(),
                                "validation_ids_sha256": hashlib.sha256("\n".join(train.ids[te]).encode()).hexdigest(),
                                "train_animals": sorted(set(train.animal[tr])),
                                "validation_animals": sorted(set(train.animal[te])), "diagnostics": model.diagnostics()})
            if progress:
                progress(f"{fold_name}: K={k} {split['name']} completed")
        if set(animal_errors) != set(train.animal):
            raise ValueError("Inner tuning does not cover every training animal exactly once.")
        selections.append({"k": k, "mean_animal_mse": float(np.mean(list(animal_errors.values()))),
                           "animal_mse": animal_errors, "fits": fit_records})
    selected = min(selections, key=lambda x: (x["mean_animal_mse"], x["k"]))["k"]
    model = fit_joint(train, Y, selected, stable_seed(CONFIG["seed"], fold_name, "outer", selected), regions=regions)
    predictions, sources = baseline_predictions(train, Y, query, return_sources=True)
    predictions.update({"joint": batch_predict(model, query.X, query.tags, regions=regions),
                        "position_only": batch_predict(model, query.X, query.tags, mode="position", regions=regions),
                        "tag_only": batch_predict(model, query.X, query.tags, mode="tag", regions=regions)})
    return predictions, {"selected_k": selected, "selected_at_grid_edge": selected in (min(CONFIG["k_grid"]), max(CONFIG["k_grid"])),
                         "selected_at_grid_low_edge": selected == min(CONFIG["k_grid"]),
                         "selected_at_grid_high_edge": selected == max(CONFIG["k_grid"]),
                         "selection": selections, "outer_fit_diagnostics": model.diagnostics(),
                         "outer_model_state": snapshot_model(model, train.ids, regions=regions), "copy_sources": sources}


def snapshot_model(model, train_ids, *, regions=REGIONS):
    """Selected fit state sufficient to independently recompute all posteriors."""
    return {"core_sha256": CORE_SHA, "fit_seed": model.seed, "field_types": model.field_types,
            "training_ids_sha256": hashlib.sha256("\n".join(train_ids).encode()).hexdigest(),
            "target_fields": list(target_fields(regions)), "field_metadata": model.field_metadata_,
            "training_metadata": model.training_metadata_,
            "blocks": [{"name": b.name, "fields": list(b.fields), "n_components": b.n_components,
                        "weights": b.params.weights.tolist(),
                        "means": {f: a.tolist() for f, a in b.params.means.items()},
                        "variances": {f: a.tolist() for f, a in b.params.variances.items()},
                        "tables": {f: a.tolist() for f, a in b.params.tables.items()}} for b in model.blocks_]}


def profile_metrics(truth, prediction, animals, *, regions=REGIONS):
    truth, prediction, animals = np.asarray(truth), np.asarray(prediction), np.asarray(animals)
    if truth.shape != prediction.shape or truth.shape != (len(animals), len(regions)):
        raise ValueError("Metric shapes differ.")
    if not np.isfinite(truth).all() or not np.isfinite(prediction).all():
        raise ValueError("Metrics require finite matrices.")
    mse, mae = ((truth - prediction)**2).mean(axis=1), np.abs(truth - prediction).mean(axis=1)
    correlations = np.full(len(truth), np.nan)
    true_variable = np.ptp(truth, axis=1) > 0
    predicted_variable = np.ptp(prediction, axis=1) > 0
    for i in np.flatnonzero(true_variable):
        correlations[i] = float(np.corrcoef(truth[i], prediction[i])[0, 1]) if predicted_variable[i] else 0.
    records = []
    for a in sorted(set(animals)):
        mask = animals == a
        valid = mask & np.isfinite(correlations)
        records.append({"animal": str(a), "neurons": int(mask.sum()), "mse": float(mse[mask].mean()),
                        "mae": float(mae[mask].mean()),
                        "profile_r": float(correlations[valid].mean()) if valid.any() else None,
                        "profile_r_valid_neurons": int(valid.sum())})
    r_animals = [r["profile_r"] for r in records if r["profile_r"] is not None]
    return {"summary": {"mse": float(np.mean([r["mse"] for r in records])),
                         "mae": float(np.mean([r["mae"] for r in records])),
                         "profile_r": float(np.mean(r_animals)) if r_animals else None,
                         "profile_r_valid_neurons": int(np.isfinite(correlations).sum()),
                         "profile_r_valid_animals": len(r_animals)}, "animals": records}


def comparisons(metrics):
    joint = metrics["joint"]
    rows = {}
    for name in MODELS:
        if name == "joint":
            continue
        baseline = metrics[name]
        if [r["animal"] for r in joint["animals"]] != [r["animal"] for r in baseline["animals"]]:
            raise ValueError("Unpaired animal comparison.")
        delta = np.asarray([j["mse"] - b["mse"] for j, b in zip(joint["animals"], baseline["animals"])])
        rng = np.random.default_rng(stable_seed(CONFIG["seed"], "brm_real_bootstrap", name))
        resampled = delta[rng.integers(0, len(delta), size=(CONFIG["bootstrap_samples"], len(delta)))].mean(axis=1)
        base_mse, joint_mse = baseline["summary"]["mse"], joint["summary"]["mse"]
        rows[name] = {"mse_delta": float(delta.mean()),
                      "relative_mse_reduction": float(1 - joint_mse / base_mse) if base_mse > 0 else None,
                      "mse_delta_bootstrap95": np.quantile(resampled, [.025, .975]).tolist(),
                      "animals_joint_better": int((delta < 0).sum()), "animals": len(delta),
                      "animal_differences": [{"animal": j["animal"], "joint_mse": j["mse"],
                                               "baseline_mse": b["mse"], "mse_delta": j["mse"] - b["mse"]}
                                              for j, b in zip(joint["animals"], baseline["animals"])],
                      "descriptive_only": True}
    beats_copies = all(joint["summary"]["mse"] < metrics[name]["summary"]["mse"] for name in COPY_MODELS)
    beats_averages = all(joint["summary"]["mse"] < metrics[name]["summary"]["mse"] for name in AVERAGES)
    return {"against": rows, "beats_all_four_copying": beats_copies,
            "beats_all_simple_averages": beats_averages, "overall_useful_prediction": beats_copies and beats_averages}


def synthetic_checks():
    """Fixed small functional invariants; does not touch repository data."""
    rng = np.random.default_rng(721)
    n = 36
    data = Inputs(np.array([f"S{i}" for i in range(n)]), np.array([f"A{i // 3}" for i in range(n)]),
                  rng.normal(size=(n, 3)), np.array(["p" if (i // 3) % 2 else "q" for i in range(n)]))
    regions = ("one", "two", "three")
    Y = np.column_stack([np.maximum(0., data.X[:, 0] + 2), np.full(n, 1.), (data.tags == "p").astype(float) * 3])
    model = fit_joint(data, Y, 2, 381, regions=regions)
    queries_X = np.vstack([data.X[:3], [np.nan, np.nan, np.nan]])
    queries_tag = np.array(["p", "unknown", "q", None], dtype=object)
    for mode in ("full", "position", "tag", "none"):
        batch = batch_predict(model, queries_X, queries_tag, mode=mode, regions=regions)
        for i in range(len(queries_X)):
            obs = {}
            if mode in ("full", "position"):
                obs.update({f: float(queries_X[i, j]) for j, f in enumerate(POSITION) if np.isfinite(queries_X[i, j])})
            if mode in ("full", "tag") and queries_tag[i] is not None:
                obs[TAG] = queries_tag[i]
            public = model.infer(obs)["predictions"]
            if not np.allclose(batch[i], [public[f]["mean"] for f in target_fields(regions)], atol=1e-10, rtol=1e-10):
                raise AssertionError("Batch adapter differs from public core inference.")
    query = Inputs(np.array(["QUERY"]), np.array(["HELDOUT"]), np.zeros((1, 3)), np.array(["absent"]))
    baselines, sources = baseline_predictions(data, Y, query, return_sources=True)
    explicit = np.asarray([np.sum(animal_weights(data.animal) * ((Y - candidate)**2).sum(axis=1)) for candidate in Y])
    expected_medoid = Y[int(np.argmin(explicit))]
    if not np.array_equal(baselines["medoid"][0], expected_medoid):
        raise AssertionError("Linear-time exact medoid identity failed.")
    if sources["same_tag_fallback_rows"] != [0] or not np.array_equal(baselines["tag_nearest"], baselines["nearest"]):
        raise AssertionError("Unknown-tag copy fallback failed.")
    return {"passed": True, "synthetic_only": True,
            "checks": ["public_infer_batch_equivalence_all_evidence_modes", "unknown_and_missing_tag_fallback",
                       "exact_animal_weighted_medoid_identity", "all_targets_hidden_from_batch_interface"]}


def check(*, base=None, root=ROOT, expected=None, pinned_inputs=None):
    base = Path(base) if base is not None else Path(root) / OUTPUT
    _stage_open(base)
    if (base / "frozen.json").exists():
        raise RuntimeError("Frozen metadata cannot be changed.")
    declared = read_json(base / "declaration.json")
    if declared != current_declaration(root, expected=expected, pinned_inputs=pinned_inputs):
        raise ValueError("Declaration is stale.")
    report = {"declaration_sha256": sha(base / "declaration.json"),
              "code_sha256": declared["code_sha256"], "core_sha256": declared["core_sha256"],
              "protocol_sha256": declared["protocol_sha256"], "runtime": runtime(), **synthetic_checks()}
    write_json(base / "check.json", report, preserve_previous=True)
    return report


def _review_valid(review_path, declared):
    text = Path(review_path).read_text()
    if REVIEW_SENTINEL not in text.splitlines():
        raise ValueError("Independent review PASS sentinel is missing.")
    for key in ("code_sha256", "core_sha256", "protocol_sha256"):
        if declared[key] not in text:
            raise ValueError(f"Review does not bind {key}.")


def freeze(*, base=None, root=ROOT, review=None, expected=None, pinned_inputs=None):
    root = Path(root)
    base = Path(base) if base is not None else root / OUTPUT
    _stage_open(base)
    declared = read_json(base / "declaration.json")
    if declared != current_declaration(root, expected=expected, pinned_inputs=pinned_inputs):
        raise ValueError("Declaration is stale.")
    checked = read_json(base / "check.json")
    if not checked.get("passed") or checked.get("declaration_sha256") != sha(base / "declaration.json"):
        raise ValueError("Matching synthetic check is missing.")
    for key in ("code_sha256", "core_sha256", "protocol_sha256"):
        if checked.get(key) != declared[key]:
            raise ValueError("Check source hash differs.")
    review_path = Path(review) if review is not None else root / REVIEW
    if not review_path.is_absolute():
        review_path = root / review_path
    review_rel = str(review_path.resolve().relative_to(root.resolve()))
    _review_valid(review_path, declared)
    frozen = {"schema": "brm_real_frozen_v1", "declaration_sha256": sha(base / "declaration.json"),
              "check_sha256": sha(base / "check.json"), "review_path": review_rel,
              "review_sha256": sha(review_path), "code_sha256": declared["code_sha256"],
              "core_sha256": declared["core_sha256"], "protocol_sha256": declared["protocol_sha256"],
              "input_sha256": declared["source"]["input_sha256"], "config": CONFIG, "runtime": runtime()}
    write_json(base / "frozen.json", frozen, exclusive=True)
    return frozen


def verify_frozen(base, root=ROOT, *, expected=None, pinned_inputs=None):
    base, root = Path(base), Path(root)
    frozen, declared = read_json(base / "frozen.json"), read_json(base / "declaration.json")
    if frozen.get("declaration_sha256") != sha(base / "declaration.json") or frozen.get("check_sha256") != sha(base / "check.json"):
        raise ValueError("Frozen declaration/check integrity mismatch.")
    if declared != current_declaration(root, expected=expected, pinned_inputs=pinned_inputs):
        raise ValueError("Frozen inputs/source/runtime no longer match.")
    for key in ("code_sha256", "core_sha256", "protocol_sha256"):
        if frozen.get(key) != declared[key]:
            raise ValueError("Frozen source hash mismatch.")
    if frozen.get("input_sha256") != declared["source"]["input_sha256"] or frozen.get("config") != CONFIG or frozen.get("runtime") != runtime():
        raise ValueError("Frozen input/config/runtime mismatch.")
    review_path = (root / frozen["review_path"]).resolve()
    review_path.relative_to(root.resolve())
    if sha(review_path) != frozen["review_sha256"]:
        raise ValueError("Frozen review changed.")
    _review_valid(review_path, declared)
    checked = read_json(base / "check.json")
    if checked.get("passed") is not True or checked.get("declaration_sha256") != sha(base / "declaration.json"):
        raise ValueError("Pinned synthetic check is not valid.")
    return frozen, declared


def _read_targets(root, base, frozen, data, *, expected_all_zero, regions=REGIONS):
    """Only guarded run/scorer call this; never prepare/check/freeze."""
    root, base = Path(root), Path(base)
    marker = read_json(base / "RUN_STARTED.json")
    if marker.get("frozen_sha256") != sha(base / "frozen.json"):
        raise ValueError("Run-open marker does not bind the frozen manifest.")
    if sha(root / INPUTS["outcome"]) != frozen["input_sha256"]["outcome"]:
        raise ValueError("Outcome raw hash changed immediately before reading.")
    table = pd.read_csv(root / INPUTS["outcome"], index_col=0)
    if list(table.columns) != list(regions) or len(table) != len(data.ids) or table.index.duplicated().any() or set(table.index) != set(data.ids):
        raise ValueError("Outcome IDs/region columns do not join exactly.")
    raw = table.loc[data.ids.tolist()].to_numpy(float)
    if not np.isfinite(raw).all() or (raw < 0).any():
        raise ValueError("Invalid regional lengths.")
    if int((raw.sum(axis=1) == 0).sum()) != expected_all_zero:
        raise ValueError("Declared all-zero profile count changed; none may be silently removed.")
    return np.log10(1. + raw)


def _safe_artifact(base, relative):
    path = (Path(base) / relative).resolve()
    path.relative_to(Path(base).resolve())
    return path


def validate_prediction_bundle(base, manifest, declared, frozen, *, regions=REGIONS):
    """Global output-complete/integrity barrier BEFORE scorer opens truth."""
    base = Path(base)
    if canonical(read_json(base / "declaration.json")) != canonical(declared) or canonical(read_json(base / "frozen.json")) != canonical(frozen):
        raise ValueError("Passed declaration/frozen records differ from saved canonical records.")
    saved = read_json(base / "prediction_manifest.json")
    if canonical(saved) != canonical(manifest):
        raise ValueError("Passed and saved prediction manifests differ.")
    if (base / "PREDICTIONS_COMPLETE").read_text().strip() != sha(base / "prediction_manifest.json"):
        raise ValueError("Prediction completion marker hash differs.")
    for key in ("code_sha256", "core_sha256", "protocol_sha256"):
        if manifest.get(key) != declared[key]:
            raise ValueError("Prediction source provenance differs.")
    if manifest.get("declaration_sha256") != sha(base / "declaration.json") or manifest.get("frozen_sha256") != sha(base / "frozen.json"):
        raise ValueError("Prediction declaration/frozen provenance differs.")
    if manifest.get("config") != CONFIG or manifest.get("input_sha256") != frozen["input_sha256"]:
        raise ValueError("Prediction input/config provenance differs.")
    expected = {f["name"]: f for f in declared["source"]["outer_folds"]}
    records = manifest.get("folds", [])
    if len(records) != len(expected) or len({f.get("name") for f in records}) != len(expected) or {f.get("name") for f in records} != set(expected):
        raise ValueError("Missing, duplicate or unexpected prediction folds.")
    ids = np.asarray(declared["source"]["ids"])
    animals = np.asarray(declared["source"]["animals"])
    tags = np.asarray(declared["source"]["tags"])
    paths, tested = set(), []
    for record in records:
        fold = expected[record["name"]]
        tr, te = np.asarray(fold["train"]), np.asarray(fold["test"])
        if record.get("train_ids") != ids[tr].tolist() or record.get("test_ids") != ids[te].tolist():
            raise ValueError("Prediction query/training identities differ.")
        if record.get("train_animals") != sorted(set(animals[tr])) or record.get("test_animals") != sorted(set(animals[te])):
            raise ValueError("Prediction animal provenance differs.")
        if record.get("selected_k") not in CONFIG["k_grid"]:
            raise ValueError("Selected K is outside the declared grid.")
        selected = record["selected_k"]
        selection = record.get("selection", [])
        if [s.get("k") for s in selection] != CONFIG["k_grid"]:
            raise ValueError("Incomplete or reordered K-selection grid.")
        for candidate in selection:
            k = candidate["k"]
            inner = declared["source"]["inner_folds"][record["name"]]
            if [s.get("inner") for s in candidate.get("fits", [])] != [f["name"] for f in inner]:
                raise ValueError("Incomplete inner fits.")
            animal_errors = candidate.get("animal_mse", {})
            if set(animal_errors) != set(animals[tr]) or any(not isinstance(v, (int, float)) or not math.isfinite(v) or v < 0 for v in animal_errors.values()):
                raise ValueError("Inner score animal coverage/finiteness differs.")
            if not math.isclose(candidate["mean_animal_mse"], float(np.mean(list(animal_errors.values()))), rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError("Inner score aggregation differs.")
            for fit, split in zip(candidate["fits"], inner):
                itr, ite = np.asarray(split["train"]), np.asarray(split["test"])
                if fit.get("fit_seed") != declared["source"]["seeds"]["inner"][record["name"]][str(k)][split["name"]]:
                    raise ValueError("Inner fit seed differs.")
                if fit.get("train_ids_sha256") != hashlib.sha256("\n".join(ids[itr]).encode()).hexdigest() or fit.get("validation_ids_sha256") != hashlib.sha256("\n".join(ids[ite]).encode()).hexdigest():
                    raise ValueError("Inner training/validation membership differs.")
                if fit.get("train_animals") != sorted(set(animals[itr])) or fit.get("validation_animals") != sorted(set(animals[ite])):
                    raise ValueError("Inner animal membership differs.")
        if selected != min(selection, key=lambda s: (s["mean_animal_mse"], s["k"]))["k"]:
            raise ValueError("Selected K was not the declared inner-score minimizer.")
        if record.get("selected_at_grid_low_edge") != (selected == min(CONFIG["k_grid"])) or record.get("selected_at_grid_high_edge") != (selected == max(CONFIG["k_grid"])) or record.get("selected_at_grid_edge") != (selected in (min(CONFIG["k_grid"]), max(CONFIG["k_grid"]))):
            raise ValueError("K-grid-edge flags differ.")
        state_path = _safe_artifact(base, record["model_state_path"])
        if record["model_state_path"] != f"fits/{record['name']}.json" or sha(state_path) != record["model_state_sha256"]:
            raise ValueError("Selected outer numeric-state integrity differs.")
        state = read_json(state_path)
        expected_fields = {*POSITION, TAG, *target_fields(regions)}
        if state.get("core_sha256") != CORE_SHA or state.get("target_fields") != list(target_fields(regions)) or set(state.get("field_types", {})) != expected_fields:
            raise ValueError("Selected state fields/core differ.")
        if state.get("fit_seed") != declared["source"]["seeds"]["outer"][record["name"]][str(selected)] or state.get("training_ids_sha256") != hashlib.sha256("\n".join(ids[tr]).encode()).hexdigest():
            raise ValueError("Selected outer state seed/membership differs.")
        training = state.get("training_metadata", {})
        if training.get("seed") != state["fit_seed"] or training.get("requested_components") != selected or training.get("informative_rows") != len(tr) or training.get("n_init") != CONFIG["joint_budget"]["n_init"]:
            raise ValueError("Selected outer training metadata differs.")
        outer_diagnostics = record.get("outer_fit_diagnostics", {})
        if outer_diagnostics.get("training") != training or len(outer_diagnostics.get("blocks", [])) != 1 or outer_diagnostics["blocks"][0].get("n_components") != min(selected, len(tr)):
            raise ValueError("Selected outer diagnostics differ from saved state.")
        if len(state.get("blocks", [])) != 1 or set(state["blocks"][0]["fields"]) != expected_fields:
            raise ValueError("Complete real training table must have one connected co-observation block.")
        block = state["blocks"][0]
        k = min(selected, len(tr))
        weights = np.asarray(block["weights"])
        if block.get("n_components") != k or weights.shape != (k,) or not np.isfinite(weights).all() or (weights <= 0).any() or not np.isclose(weights.sum(), 1.):
            raise ValueError("Invalid saved mixture weights.")
        for field in (*POSITION, *target_fields(regions)):
            means, variances = np.asarray(block["means"][field]), np.asarray(block["variances"][field])
            metadata = state["field_metadata"][field]
            if means.shape != (k,) or variances.shape != (k,) or not np.isfinite(means).all() or not np.isfinite(variances).all() or (variances <= 0).any() or not math.isfinite(metadata["training_scale"]) or metadata["training_scale"] <= 0:
                raise ValueError("Invalid saved continuous parameters.")
        table = np.asarray(block["tables"][TAG])
        if table.shape != (k, len(state["field_metadata"][TAG]["categories"])) or not np.isfinite(table).all() or (table <= 0).any() or not np.allclose(table.sum(axis=1), 1.):
            raise ValueError("Invalid saved tag probabilities.")
        path = _safe_artifact(base, record["path"])
        if path in paths or record["path"] != f"predictions/{record['name']}.npz":
            raise ValueError("Duplicate/unexpected prediction path.")
        paths.add(path)
        if sha(path) != record["sha256"]:
            raise ValueError("Saved prediction NPZ checksum differs.")
        sources = record.get("copy_sources", {})
        train_ids = set(ids[tr])
        id_to_index = {str(ids[i]): i for i in tr}
        for name in COPY_MODELS:
            values = sources.get(name, [])
            if len(values) != len(te) or not set(values) <= train_ids:
                raise ValueError("Copy source is missing or from a held-out animal.")
            for q, source in zip(te, values):
                if animals[id_to_index[source]] in set(animals[te]):
                    raise ValueError("A query animal supplied a copying profile.")
                if name.startswith("tag_") and tags[q] in set(tags[tr]) and tags[id_to_index[source]] != tags[q]:
                    raise ValueError("A same-tag copy used a different source tag.")
        with np.load(path, allow_pickle=False) as arrays:
            if set(arrays.files) != {"indices", "ids", *MODELS} or not np.array_equal(arrays["indices"], te) or not np.array_equal(arrays["ids"], ids[te]):
                raise ValueError("Saved prediction row/model completeness differs.")
            for name in MODELS:
                values = arrays[name]
                if values.shape != (len(te), len(regions)) or not np.isfinite(values).all() or (values < 0).any():
                    raise ValueError("Invalid saved regional predictions.")
        tested.extend(te.tolist())
    if sorted(tested) != list(range(len(ids))):
        raise ValueError("Predictions do not cover exact source IDs once.")


def score_predictions(base, manifest, frozen, declared, *, root=ROOT, data=None, regions=REGIONS):
    base, root = Path(base), Path(root)
    validate_prediction_bundle(base, manifest, declared, frozen, regions=regions)
    if data is None:
        data = load_inputs(root, expected=declared["expected"])
    if data.ids.tolist() != declared["source"]["ids"] or data.animal.tolist() != declared["source"]["animals"] or data.tags.tolist() != declared["source"]["tags"]:
        raise ValueError("Scorer input metadata differ from exact declared order.")
    # Truth is opened only after ALL folds/models have passed output integrity.
    truth = _read_targets(root, base, frozen, data, expected_all_zero=declared["expected"]["all_zero"], regions=regions)
    truth_path = base / "truth.npz"
    with truth_path.open("xb") as output:
        np.savez_compressed(output, ids=data.ids, Y_log=truth, regions=np.asarray(regions))
    predictions = {name: np.empty_like(truth) for name in MODELS}
    for record in manifest["folds"]:
        path = _safe_artifact(base, record["path"])
        if sha(path) != record["sha256"]:
            raise ValueError("Prediction changed after the barrier.")
        with np.load(path, allow_pickle=False) as arrays:
            indices = arrays["indices"]
            for name in MODELS:
                predictions[name][indices] = arrays[name]
    primary = {name: profile_metrics(truth, values, data.animal, regions=regions) for name, values in predictions.items()}
    by_fold = {}
    for fold in declared["source"]["outer_folds"]:
        te = np.asarray(fold["test"])
        by_fold[fold["name"]] = {name: profile_metrics(truth[te], values[te], data.animal[te], regions=regions)["summary"]
                                 for name, values in predictions.items()}
    keep = [i for i, r in enumerate(regions) if r != "MY"]
    sensitivity = {name: profile_metrics(truth[:, keep], values[:, keep], data.animal,
                                         regions=tuple(regions[i] for i in keep)) for name, values in predictions.items()}
    result = {"config": CONFIG, "prediction_manifest_sha256": sha(base / "prediction_manifest.json"),
              "truth_path": "truth.npz", "truth_sha256": sha(truth_path),
              "truth_raw_sha256": frozen["input_sha256"]["outcome"], "neurons": len(data.ids),
              "animals": len(set(data.animal)), "all_zero_profiles_preserved": declared["expected"]["all_zero"],
              "models": primary, "fold_metrics": by_fold, "comparison": comparisons(primary),
              "fit_status": {"selected_k": {f["name"]: f["selected_k"] for f in manifest["folds"]},
                             "selected_at_low_edge": sum(f["selected_at_grid_low_edge"] for f in manifest["folds"]),
                             "selected_at_high_edge": sum(f["selected_at_grid_high_edge"] for f in manifest["folds"]),
                             "outer_nonconverged_blocks": sum(not b["converged"] for f in manifest["folds"] for b in f["outer_fit_diagnostics"]["blocks"]),
                             "inner_nonconverged_blocks": sum(not b["converged"] for f in manifest["folds"] for k in f["selection"] for g in k["fits"] for b in g["diagnostics"]["blocks"])},
              "without_MY": {"regions": [regions[i] for i in keep], "models": sensitivity,
                             "comparison": comparisons(sensitivity)},
              "limitations": ["Regional profiles, not 3D paths or synapses.",
                              "Tag is a source-line label confounded with animal/batch, not gene expression or causation.",
                              "Gaussian log-target model used for point means only; no calibrated physical intervals or likelihood claim.",
                              "Conditional ablations share the full-input-selected joint fit and are not separately optimized models.",
                              "Descriptive paired animal bootstrap, not independent confirmatory inference."]}
    write_json(base / "metrics.json", result, exclusive=True)
    return result


def run(*, base=None, root=ROOT, expected=None, pinned_inputs=None, progress=print):
    root = Path(root)
    base = Path(base) if base is not None else root / OUTPUT
    _stage_open(base)
    frozen, declared = verify_frozen(base, root, expected=expected, pinned_inputs=pinned_inputs)
    start = time.monotonic()
    write_json(base / "RUN_STARTED.json", {"frozen_sha256": sha(base / "frozen.json"),
               "started_utc": datetime.now(timezone.utc).isoformat()}, exclusive=True)
    try:
        data = load_inputs(root, expected=declared["expected"])
        all_targets = _read_targets(root, base, frozen, data,
                                    expected_all_zero=declared["expected"]["all_zero"])
        manifest = {"config": CONFIG, "runtime": runtime(), "code_sha256": declared["code_sha256"],
                    "core_sha256": declared["core_sha256"], "protocol_sha256": declared["protocol_sha256"],
                    "declaration_sha256": sha(base / "declaration.json"), "frozen_sha256": sha(base / "frozen.json"),
                    "input_sha256": frozen["input_sha256"], "folds": []}
        for fold in declared["source"]["outer_folds"]:
            tr, te = np.asarray(fold["train"]), np.asarray(fold["test"])
            lookup = {int(value): index for index, value in enumerate(tr)}
            inner = [{"name": f["name"], "train": [lookup[i] for i in f["train"]],
                      "test": [lookup[i] for i in f["test"]]} for f in declared["source"]["inner_folds"][fold["name"]]]
            predictions, record = fit_outer_fold(data.subset(tr), all_targets[tr], data.subset(te), inner,
                                                  fold_name=fold["name"], progress=progress)
            state_relative = f"fits/{fold['name']}.json"
            write_json(base / state_relative, record.pop("outer_model_state"), exclusive=True)
            record.update({"model_state_path": state_relative, "model_state_sha256": sha(base / state_relative)})
            relative = f"predictions/{fold['name']}.npz"
            path = base / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("xb") as output:
                np.savez_compressed(output, indices=te, ids=data.ids[te], **predictions)
            record.update({"name": fold["name"], "train_ids": data.ids[tr].tolist(), "test_ids": data.ids[te].tolist(),
                           "train_animals": sorted(set(data.animal[tr])), "test_animals": sorted(set(data.animal[te])),
                           "path": relative, "sha256": sha(path)})
            manifest["folds"].append(record)
            write_json(base / "prediction_manifest.partial.json", manifest)
            if progress:
                progress(f"{fold['name']}: all10 predictions saved; selectedK={record['selected_k']}")
        del all_targets  # scorer reloads only after the global complete/hash barrier
        write_json(base / "prediction_manifest.json", manifest, exclusive=True)
        with (base / "PREDICTIONS_COMPLETE").open("x") as f:
            f.write(sha(base / "prediction_manifest.json") + "\n")
        result = score_predictions(base, manifest, frozen, declared, root=root, data=data)
        completion = {"elapsed_seconds": time.monotonic() - start, "metrics_sha256": sha(base / "metrics.json"),
                      "prediction_manifest_sha256": sha(base / "prediction_manifest.json")}
        write_json(base / "RUN_COMPLETE.json", completion, exclusive=True)
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
        print(canonical({"stage": "metadata_prepared", "neurons": len(result["source"]["ids"]),
                         "declaration_sha256": sha(ROOT / OUTPUT / "declaration.json")}))
    elif args.stage == "check":
        print(canonical(check()))
    elif args.stage == "freeze":
        print(canonical(freeze(review=args.review)))
    else:
        result = run(progress=lambda message: print(message, flush=True))
        print(canonical({"comparison": result["comparison"],
                         "summary": {name: row["summary"] for name, row in result["models"].items()}}))


if __name__ == "__main__":
    main()
