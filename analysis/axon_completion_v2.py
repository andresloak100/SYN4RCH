"""One error-informed TRAIN-development follow-up to frozen axon-completion v1.

The mask, decoder, folds, pools, labels and fallback choices are unchanged.
The single candidate adds fit-only donor geometry landmarks and a cable-gain
head. No real fitting is permitted until ROOT clears this protocol FINAL and
the new source/protocol/tests are committed. Never opens external traces.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/interim/axon_completion_v2"
OLD = ROOT / "data/interim/axon_completion_v1"
PROTOCOL = ROOT / "reports/axon_completion_v2_protocol.md"
V1_SHA = "27cc42508c7c3872be0a5e7ff1a3e884e6d2a11cdce32db6959ac59302bb31ed"
SELECTOR = "completion_selector"
V1_CONTROL = "completion_selector_v1"
LANDMARKS = 32
_V1 = None
_WORKER_LIB = _WORKER_OUT = None

# Exact immutable completed-v1 inputs. The score-label file is only loaded by
# evaluation AFTER the new forecast seal; hashing a file does not admit labels.
OLD_PINS = {
    "RUN_STARTED.json": "743654b4a54f8adcf3425f6a43268a81ae47956c6b127539a7cfb831abe7fc35",
    "RUN_COMPLETE.json": "bcf7e65184033cee175aad1dc7880cc8e2eb0dd0489e8ace8c0483ed24e67865",
    "index.jsonl": "e0cb708152dbdc1a9dd5881e2721091c3bf8e4d52ecc8dd28528858065862ade",
    "folds.json": "0c0413f50690899cc8f6f72ce996887ee6ff79f038caade86755a13f66c47e8d",
    "fold_info.json": "12d5423e74a1aca607a9a9c4eeae0786fb54846207571ff027803ede3e604916",
    "forecasts.json": "118c94184997f202d1ad4350a16ff0680d3e4fa44e27b26743c2a4fef5f6cfa0",
    "FORECASTS_SEALED.json": "db83d7a2b819eec5d84ceca68172df8975b0b813fdc2653bd806e3e603b64b99",
    "labels_fit.jsonl": "b82a1e323266a65d81cecbcbdc9bba10a16861f926d51e45f5d0d0001f67c2e1",
    "labels_score.jsonl": "46997f1f03a1359635bfed06d7b1be0b7adfc2525f368be266b5e37c4463e3f6",
    "results.json": "bc76687b0bdd99063dca2aa82a089bf2522cec21afdd268c8fed5bd5773ba428",
}


def V1():
    global _V1
    if _V1 is None:
        import hashlib
        path = ROOT / "analysis/axon_completion_v1.py"
        if hashlib.sha256(path.read_bytes()).hexdigest() != V1_SHA:
            raise RuntimeError("Frozen v1 source changed")
        spec = importlib.util.spec_from_file_location("axon_completion_v2_frozen_v1", path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        _V1 = module
    return _V1


def verify_inputs():
    A = V1()
    for name, expected in OLD_PINS.items():
        if A.sha_file(OLD / name) != expected:
            raise RuntimeError(f"Frozen v1 input changed: {name}")
    started = A.read_json(OLD / "RUN_STARTED.json")
    for filename, expected in started["code_sha256"].items():
        if A.sha_file(ROOT / filename) != expected:
            raise RuntimeError(f"Frozen dependency changed: {filename}")
    if A.read_json(OLD / "RUN_COMPLETE.json")["results_sha256"] != OLD_PINS["results.json"]:
        raise RuntimeError("v1 completion receipt disagrees")
    return {**{f"v1/{k}": v for k, v in OLD_PINS.items()}, **started["code_sha256"]}


class GeometryLandmarks:
    """Fitted donor representation. Every D access validates BOTH endpoints."""

    def __init__(self, lib, door, count=LANDMARKS):
        self.lib, self.door = lib, door
        fit = np.flatnonzero(door.mask)
        if len(fit) == 0 or count <= 0:
            raise ValueError("A nonempty fit pool and positive landmark count are required")
        D = np.asarray(door.D(fit, fit), float)
        per_animal = Counter(lib.animal[fit].tolist())
        w = np.array([1.0 / per_animal[a] for a in lib.animal[fit]])
        risk = D @ (w / w.sum())
        first = min(range(len(fit)), key=lambda i: (risk[i], lib.ids[int(fit[i])]))
        chosen, used = [first], {lib.animal[fit[first]]}
        nearest = D[:, first].copy()
        limit = min(int(count), len(per_animal))
        while len(chosen) < limit:
            eligible = [i for i, d in enumerate(fit) if lib.animal[d] not in used]
            selected = min(eligible, key=lambda i: (-nearest[i], lib.ids[int(fit[i])]))
            chosen.append(selected)
            used.add(lib.animal[fit[selected]])
            nearest = np.minimum(nearest, D[:, selected])
        self.indices = fit[chosen]
        self.fit_ids = tuple(lib.ids[int(d)] for d in fit)
        self.ids = tuple(lib.ids[int(d)] for d in self.indices)

    def rows(self, donors):
        donors = np.asarray(donors, int)
        return np.log1p(np.asarray(self.door.D(donors, self.indices), float) / 200.0)


class Featurizer:
    """The unchanged v1 prefix features, plus donor/relative landmark distances."""

    def __init__(self, lib, idx, door):
        self.base = V1().Featurizer(lib, idx, door, prefix=True)
        self.landmarks = GeometryLandmarks(lib, door)
        self.names = None

    def rows(self, local, fallback_donor, *, candidates=None):
        cand = local.cand if candidates is None else np.asarray(candidates, int)
        base = self.base.rows(local, fallback_donor, candidates=cand)
        donor = self.landmarks.rows(cand)
        relative = donor - self.landmarks.rows(np.array([fallback_donor]))
        # Mask the query's own sample: donor-to-own-query landmark distances
        # would be an easy fitting shortcut unavailable to outer-eval queries.
        same_sample = self.base.lib.animal[self.landmarks.indices] == self.base.lib.animal[local.q]
        donor[:, same_sample] = np.nan
        relative[:, same_sample] = np.nan
        names = list(self.base.names) + [f"donor_route_lm_{i:02d}" for i in range(donor.shape[1])]
        names += [f"donor_route_lm_delta_{i:02d}" for i in range(donor.shape[1])]
        V1().V2().check_feature_names(names)
        if self.names is not None and self.names != names:
            raise RuntimeError("Feature schema changed inside a fit")
        self.names = names
        return np.column_stack([base, donor, relative])


def cable_error(label):
    if not label.get("valid", False):
        raise ValueError("Cable targets require valid fitting labels")
    hidden, extra = label.get("hidden_cable_um"), label.get("extra_cable_um")
    if hidden is None or extra is None or not math.isfinite(hidden) or hidden <= 0 or not math.isfinite(extra) or extra < 0:
        raise ValueError("A valid label lacks finite positive hidden cable and nonnegative extra cable")
    return abs(extra - hidden) / hidden


def training_rows(lib, door, feat, locals_, fallback_of):
    """Exactly v1 rows/weights/gd/gf; one additional valid-label cable target."""
    A = V1()
    X, gd, gf, gc, w, skipped = [], [], [], [], [], 0
    counts = Counter(lib.animal[q] for q in locals_)
    for q in sorted(locals_):
        loc, fb = locals_[q], fallback_of[q]
        cand = np.sort(loc.cand)
        lf = door.label(q, fb)
        if not lf["valid"]:
            skipped += len(cand)
            continue
        fb_cable = cable_error(lf)
        rows, keep = feat.rows(loc, fb, candidates=cand), []
        for i, d in enumerate(cand):
            label = door.label(q, int(d))
            if not label["valid"]:
                skipped += 1
                continue
            keep.append(i)
            gd.append(A.logd(lf["distance_um"]) - A.logd(label["distance_um"]))
            gf.append(label["f1_200"] - lf["f1_200"])
            gc.append(fb_cable - cable_error(label))
        X.append(rows[keep])
        w.extend([1.0 / (counts[lib.animal[q]] * max(len(keep), 1))] * len(keep))
    if not w:
        raise RuntimeError("No valid fitting rows")
    w = np.asarray(w); w /= w.mean()
    return np.vstack(X), np.asarray(gd), np.asarray(gf), np.asarray(gc), w, skipped


class ResidualHeads:
    def fit(self, X, gd, gf, gc, w):
        self.models = [V1().make_hgb() for _ in range(3)]
        for model, target in zip(self.models, (gd, gf, gc)):
            model.fit(X, target, sample_weight=w)
        return self

    def predict(self, X):
        return tuple(model.predict(X) for model in self.models)


def residual_choice(candidates, fallback, gd, gf, gc):
    candidates = np.asarray(candidates, int)
    gd, gf, gc = np.asarray(gd, float), np.asarray(gf, float), np.asarray(gc, float)
    if any(x.shape != candidates.shape for x in (gd, gf, gc)):
        raise ValueError("Gain arrays and candidates differ")
    if not all(np.isfinite(x).all() for x in (gd, gf, gc)):
        return int(fallback), "non_finite_prediction"
    ok = (gd > 0) & (gf > 0) & (gc >= 0) & (candidates != fallback)
    if not ok.any():
        return int(fallback), "fallback"
    score = np.where(ok, gd + gf, -np.inf)
    top = int(np.lexsort((candidates, -score))[0])
    return int(candidates[top]), "switched"


def make_controls(lib, loc, fallback, choice, why, permuted_choice, old):
    """Frozen fixed arms, but random switching follows THIS selector's switches."""
    A = V1()
    donors = dict(old["donors"])
    donors[V1_CONTROL] = donors[SELECTOR]
    donors[SELECTOR] = lib.ids[choice]
    donors["permuted_labels"] = lib.ids[permuted_choice]
    if why == "switched":
        others = sorted(int(d) for d in loc.cand if int(d) != fallback)
        pick = np.random.default_rng(A.stable_seed("switch_random", lib.ids[loc.q])).integers(len(others))
        donors["switch_random"] = lib.ids[others[pick]]
    else:
        donors["switch_random"] = lib.ids[fallback]
    return donors


def load_context():
    """Reuse index/known metadata/D only. Never builds a new biological index."""
    A = V1(); lib = A.V2().Library(); cache = A.V2().Cache(lib)
    rows = [json.loads(line) for line in (OLD / "index.jsonl").read_text().splitlines() if line.strip()]
    by_id = {row["id"]: row for row in rows}
    if len(by_id) != len(rows) or set(by_id) != set(lib.ids):
        raise RuntimeError("Frozen index does not cover exactly TRAIN")
    idx = A.Index(lib, [by_id[n] for n in lib.ids])
    plan = A.read_json(OLD / "folds.json")
    if plan != A.fold_plan(sorted(set(lib.animal.tolist()))):
        raise RuntimeError("Frozen five/four fold definitions differ")
    return lib, cache, idx, plan


def locals_for(lib, cache, idx, fold, which):
    A = V1()
    door = A.V2().FitDoor(f"v2_outer{fold['outer']}_{which}", lib, fold["fit_animals"], cache=cache)
    return {q: A.Local(lib, idx, door, q, A.pool_mask(lib, fold["fit_animals"], q))
            for q in A.neurons_of(lib, idx, fold[f"{which}_animals"])}


def assert_local_parity(lib, loc, old):
    """Prove reconstructed observation candidate/control rules are v1's."""
    A = V1()
    for method in A.FALLBACK_SET:
        if lib.ids[loc.choice[method]] != old["donors"][method]:
            raise RuntimeError(f"Frozen local control differs: {method}")
    candidates = {lib.ids[int(d)] for d in loc.cand}
    for method in (SELECTOR, "metadata_only", "switch_random", "permuted_labels"):
        if old["donors"][method] not in candidates:
            raise RuntimeError(f"Frozen {method} is outside the reconstructed candidate union")


def plan_only():
    A = V1(); verify_inputs(); lib, cache, idx, plan = load_context()
    labels = A.Labels(lib); labels.load(OLD / "labels_fit.jsonl")
    old_forecasts = A.read_json(OLD / "forecasts.json")
    count, missing = 0, 0
    per_fold = []
    for fold in plan:
        locs = locals_for(lib, cache, idx, fold, "fit")
        n = sum(len(loc.cand) for loc in locs.values())
        absent = sum(not labels.has(q, int(d)) for q, loc in locs.items() for d in loc.cand)
        count += n; missing += absent
        per_fold.append({"outer": fold["outer"], "fitting_queries": len(locs), "fitting_pairs": n,
                         "missing_fitting_pairs": absent})
        for q, loc in locals_for(lib, cache, idx, fold, "eval").items():
            assert_local_parity(lib, loc, old_forecasts[lib.ids[q]])
    info = {"fitting_row_occurrences": count, "missing_fitting_pair_occurrences": missing,
            "cached_unique_fitting_pairs": len(labels.data), "folds": per_fold,
            "new_eval_pair_upper_bound": 3 * int(idx.eligible.sum()),
            "evaluation_scores_loaded": False, "new_geometry_or_labels_computed": False,
            "frozen_eval_local_parity_queries": len(old_forecasts),
            "runtime_plan": "Five outer fits, three HGB heads plus three permuted heads each. No inner refit/search. "
                            "Recompute zero fitting labels; after seal, only up to three selected new donors/query "
                            "missing from both immutable v1 caches. Decode three selected arms/query once, with bounded donor LRU."}
    OUT.mkdir(parents=True, exist_ok=True); A.write_json(OUT / "plan.json", info)
    if missing:
        raise RuntimeError(f"The frozen cache misses {missing} required fitting pairs")
    return info


def fit_forecasts(lib, cache, idx, plan, labels, out, log):
    A = V1(); old_forecasts = A.read_json(OLD / "forecasts.json"); old_info = A.read_json(OLD / "fold_info.json")
    forecasts, info = {}, {}
    for fold in plan:
        k = fold["outer"]
        checkpoint, receipt = out / f"outer_{k}.json", out / f"outer_{k}_complete.json"
        if receipt.exists():
            if A.sha_file(checkpoint) != A.read_json(receipt)["sha256"]:
                raise RuntimeError("Outer checkpoint changed")
            saved = A.read_json(checkpoint); forecasts.update(saved["forecasts"]); info[str(k)] = saved["info"]
            continue
        fit_locs = locals_for(lib, cache, idx, fold, "fit")
        eval_locs = locals_for(lib, cache, idx, fold, "eval")
        for q, loc in eval_locs.items():
            assert_local_parity(lib, loc, old_forecasts[lib.ids[q]])
        fallback = old_info[str(k)]["fallback"]
        door = A.LabelDoor(f"v2_outer{k}_labels", lib, fold["fit_animals"], labels, cache)
        feat = Featurizer(lib, idx, door)
        fb_of = {q: loc.choice[fallback] for q, loc in fit_locs.items()}
        X, gd, gf, gc, w, skipped = training_rows(lib, door, feat, fit_locs, fb_of)
        model = ResidualHeads().fit(X, gd, gf, gc, w)
        perm = np.random.default_rng(A.stable_seed("permuted", k)).permutation(len(gd))
        permuted = ResidualHeads().fit(X, gd[perm], gf[perm], gc[perm], w)
        probe = next(iter(eval_locs))
        try:
            door.label(probe, int(eval_locs[probe].cand[0]))
        except A.leak_types():
            pass
        else:
            raise RuntimeError("Outer label probe was served")
        ff = {}
        for q, loc in eval_locs.items():
            qid, fb = lib.ids[q], loc.choice[fallback]
            old = old_forecasts[qid]
            if old["outer"] != k or old["obs_sha"] != idx.obs_sha[q] or old["donors"]["adaptive_fallback"] != lib.ids[fb]:
                raise RuntimeError("Frozen query/fold/observation/fallback differs")
            features = feat.rows(loc, fb)
            gains = model.predict(features); choice, why = residual_choice(loc.cand, fb, *gains)
            pchoice, pwhy = residual_choice(loc.cand, fb, *permuted.predict(features))
            donors = make_controls(lib, loc, fb, choice, why, pchoice, old)
            hashes = dict(old["completion_sha256"]); hashes[V1_CONTROL] = hashes[SELECTOR]
            for method in (SELECTOR, "permuted_labels", "switch_random"):
                hashes.pop(method)
            ff[qid] = {"outer": k, "why": why, "permuted_why": pwhy,
                       "metadata_why": old.get("metadata_why"), "v1_why": old["why"],
                       "obs_sha": idx.obs_sha[q], "donors": donors, "completion_sha256": hashes,
                       "candidates": [lib.ids[int(d)] for d in loc.cand]}
        fi = {"fallback": fallback, "metadata_fallback": old_info[str(k)]["metadata_fallback"],
              "train_rows": len(gd), "skipped_invalid_rows": skipped, "features": feat.names,
              "landmark_ids": list(feat.landmarks.ids), "door_refused": dict(door.refused),
              "fallback_reselected": False}
        h = A.write_json(checkpoint, {"forecasts": ff, "info": fi})
        A.write_json(receipt, {"utc": A.utc(), "sha256": h}, exclusive=True)
        forecasts.update(ff); info[str(k)] = fi
        log(f"outer {k} complete: {len(gd)} fitting rows, {len(feat.landmarks.ids)} landmarks")
    return forecasts, info


def init_worker(out):
    global _WORKER_LIB, _WORKER_OUT
    _WORKER_LIB, _WORKER_OUT = V1().V2().Library(), Path(out)
    donor_suffix.cache_clear()


@lru_cache(maxsize=64)
def donor_suffix(nid):
    return V1().DonorSuffix(_WORKER_LIB.load(nid))


def decode_new_job(job):
    A = V1(); qid = job["query"]
    # Acquisition produces only the public prefix; no query suffix is passed on.
    obs = A.C().observe(_WORKER_LIB.load(qid), A.BUDGET_UM)
    if obs.sha256 != job["obs_sha"]:
        raise RuntimeError("Observation differs from frozen v1")
    decoded = {}
    for did in job["donors"]:
        if _WORKER_LIB.animal[_WORKER_LIB.index[did]] == _WORKER_LIB.animal[_WORKER_LIB.index[qid]]:
            raise A.LeakError("Same-animal forecast donor")
        try:
            completion, _ = A.decode(obs, donor_suffix(did))
            decoded[did] = {"completion_sha256": A.completion_sha(completion)}
        except A.CableGuardError as error:
            decoded[did] = {"completion_sha256": None, "error": str(error)}
    return {"job": job["job"], "query": qid, "decoded": decoded}


def evaluation_admitted(out):
    A = V1(); out = Path(out)
    receipt = A.read_json(out / "FORECASTS_SEALED.json")
    if A.sha_file(out / "forecasts.json") != receipt["forecasts_sha256"]:
        raise RuntimeError("Forecast seal differs")
    for filename, key in (("fold_info.json", "fold_info_sha256"), ("RUN_STARTED.json", "run_started_sha256")):
        if A.sha_file(out / filename) != receipt[key]:
            raise RuntimeError(f"Sealed provenance differs: {filename}")
    return receipt


def load_evaluation_labels(lib, out):
    evaluation_admitted(out)
    labels = V1().Labels(lib)
    labels.load(OLD / "labels_fit.jsonl")
    labels.load(OLD / "labels_score.jsonl")
    labels.load(Path(out) / "labels_missing.jsonl")
    return labels


def missing_label_job(job):
    """Only runs after a new forecast seal; the unchanged v1 geometry meaning."""
    evaluation_admitted(_WORKER_OUT)
    A = V1(); qid = job["query"]
    obs, truth = A.C().split_arbor(_WORKER_LIB.load(qid), A.BUDGET_UM)
    if obs.sha256 != job["obs_sha"]:
        raise RuntimeError("Observation differs")
    pairs = []
    for did in job["donors"]:
        if _WORKER_LIB.animal[_WORKER_LIB.index[did]] == _WORKER_LIB.animal[_WORKER_LIB.index[qid]]:
            raise A.LeakError("Same-animal score donor")
        try:
            completion, diag = A.decode(obs, donor_suffix(did))
        except A.CableGuardError as error:
            pairs.append({"donor": did, "valid": False, "status": "failed_cable_guard", "error": str(error)})
            continue
        m = A.C().metrics_to_json(A.C().score_completion(completion, truth, step_um=A.STEP_UM, tolerance_um=A.TOLERANCE_UM))
        pairs.append({"donor": did, "status": m["status"], "valid": bool(m["valid_for_comparison"]),
                      "distance_um": m["distance_um"], "f1_200": m["f1_200"],
                      "extra_cable_um": m["extra_cable_um"], "hidden_cable_um": m["hidden_cable_um"],
                      "modes": diag.get("modes", {}), "merged_duplicate_cable_um": diag.get("merged_duplicate_cable_um", 0.0),
                      "completion_sha256": A.completion_sha(completion)})
    return {"job": job["job"], "query": qid, "pairs": pairs}


def verdict(summary, versus, validity, switch_rate):
    """v1 gate unchanged, with frozen v1 selector added as a required control."""
    A = V1(); result = A.verdict(summary, versus, validity, switch_rate)
    v, val = versus[V1_CONTROL], validity[V1_CONTROL]
    unevaluable = val["invalid"] > 0 or v.get("nonfinite") or v["distance_delta_um"] is None
    beaten = not unevaluable and (v["distance_delta_um"] <= -A.WIN_DISTANCE_UM and v["f1_delta"] >= A.WIN_F1
                                  and v["distance_bootstrap95_um"][1] < 0 and v["f1_bootstrap95"][0] > 0)
    result["margins"][V1_CONTROL] = {"beaten": bool(beaten), "unevaluable": bool(unevaluable)}
    if not beaten:
        result["reasons"].append("not a resolved practical win over frozen completion_selector_v1")
    if validity[SELECTOR]["inflated_share"] > val["inflated_share"] + A.INFLATION_SLACK:
        result["reasons"].append("cable guard: inflated share exceeds frozen completion_selector_v1 plus0.02")
    result["pass"] = not result["reasons"]
    result["sentence"] = "PASS: development completion gate met" if result["pass"] else "FAIL: development completion gate not met"
    result["scope"] = "Error-informed DEVELOPMENT bundle on spent TRAIN; no independent confirmation or promotion"
    return result


def summarise(lib, forecasts, scores, fold_info, out):
    A = V1(); methods = A.METHODS + (V1_CONTROL,); controls = A.CONTROLS + (V1_CONTROL,)
    rows_by, validity, grafts = {}, {}, {}
    for method in methods:
        rows, bad, errors, ratios, modes = {}, Counter(), defaultdict(list), [], Counter()
        merged = 0.0
        for qid, forecast in forecasts.items():
            q = lib.index[qid]; did = forecast["donors"][method]
            label = scores.get(q, None if did is None else lib.index[did])
            expected = forecast["completion_sha256"][method]
            if label.get("completion_sha256") != expected:
                raise RuntimeError(f"{qid}/{method}: cache score is not the sealed emitted geometry")
            rows[q] = label
            if not label["valid"]:
                bad[label["status"]] += 1
            hidden, extra = label.get("hidden_cable_um"), label.get("extra_cable_um")
            if hidden:
                ratios.append((extra or 0.0) / hidden)
            errors[lib.animal[q]].append(abs(extra - hidden) / hidden if hidden and extra is not None else math.inf)
            modes.update(label.get("modes", {})); merged += float(label.get("merged_duplicate_cable_um") or 0.0)
        rows_by[method] = rows; grafts[method] = dict(modes)
        validity[method] = {"invalid": sum(bad.values()), "by_status": dict(bad),
                            "inflated_share": float(np.mean(np.asarray(ratios) > A.INFLATED_RATIO)),
                            "mean_abs_rel_cable_error": float(np.mean([np.mean(v) for v in errors.values()])),
                            "merged_duplicate_cable_um": merged}
    means = {m: A.sample_means(lib, rows) for m, rows in rows_by.items()}
    animals = sorted(means[SELECTOR]); summary = {m: {"distance_um": float(np.mean([means[m][a][0] for a in animals])),
                                                    "f1_200": float(np.mean([means[m][a][1] for a in animals]))} for m in methods}
    versus = {m: A.paired_bootstrap(means[SELECTOR], means[m], animals) for m in controls + A.DIAGNOSTIC}
    switch_rate = float(np.mean([f["why"] == "switched" for f in forecasts.values()]))
    result = {"utc": A.utc(), "neurons": len(forecasts), "samples": len(animals), "summary": summary,
              "validity": validity, "versus": versus, "executed_grafts": grafts, "fold_info": fold_info,
              "why_counts": dict(Counter(f["why"] for f in forecasts.values())),
              "verdict": verdict(summary, versus, validity, switch_rate),
              "scope": "Error-informed development; new representation+cable head is a bundle, no component attribution"}
    A.write_json(out / "results.json", A.finite_json(result)); return result


def code_pins():
    A = V1(); files = ["analysis/axon_completion_v2.py", "tests/test_axon_completion_v2.py", "reports/axon_completion_v2_protocol.md"]
    return {name: A.sha_file(ROOT / name) for name in files}


def run(workers):
    import multiprocessing as mp
    A = V1()
    if workers < 1 or workers > 4 or any(os.environ.get(v) != "1" for v in A.THREAD_VARS):
        raise RuntimeError("Use1–4 workers and all four thread variables1")
    if "STATUS: FINAL" not in PROTOCOL.read_text() or "ROOT_CLEARANCE: APPROVED" not in PROTOCOL.read_text():
        raise RuntimeError("Wait for ROOT's explicit FINAL/fitting clearance")
    owned = ["analysis/axon_completion_v2.py", "tests/test_axon_completion_v2.py", "reports/axon_completion_v2_protocol.md"]
    clean, head = A.git_clean(owned)
    if not clean:
        raise RuntimeError("New owned files must be committed and clean")
    inputs = verify_inputs(); pins = code_pins(); OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "RUN_COMPLETE.json").exists() or (OUT / "RUN_FAILED.json").exists():
        raise RuntimeError("Completed/failed experiment cannot be repeated")
    start = OUT / "RUN_STARTED.json"
    if start.exists():
        old_start = A.read_json(start)
        if old_start["code_sha256"] != pins or old_start["input_sha256"] != inputs:
            raise RuntimeError("Source/protocol/input changes: resume refused")
    else:
        A.write_json(start, {"utc": A.utc(), "git_head": head, "code_sha256": pins, "input_sha256": inputs,
                             "workers": workers, "v1_seed": A.SEED, "landmarks": LANDMARKS}, exclusive=True)
    log = A.Log(OUT / "run.log"); lib, cache, idx, plan = load_context()
    try:
        with mp.get_context("spawn").Pool(workers, initializer=init_worker, initargs=(str(OUT),)) as pool:
            if not (OUT / "FORECASTS_SEALED.json").exists():
                labels = A.Labels(lib); labels.load(OLD / "labels_fit.jsonl")
                forecasts, fold_info = fit_forecasts(lib, cache, idx, plan, labels, OUT, log)
                A.write_json(OUT / "fold_info.json", fold_info)
                jobs = [{"job": f"decode|{q}", "query": q, "obs_sha": f["obs_sha"],
                         "donors": sorted({f["donors"][m] for m in (SELECTOR, "switch_random", "permuted_labels")})}
                        for q, f in sorted(forecasts.items())]
                A.run_jobs(decode_new_job, jobs, OUT / "decoded_new.jsonl", pool, log, "decode new arms")
                decoded = {z["query"]: z["decoded"] for z in map(json.loads, (OUT / "decoded_new.jsonl").read_text().splitlines())}
                for q, f in forecasts.items():
                    for method in (SELECTOR, "switch_random", "permuted_labels"):
                        f["completion_sha256"][method] = decoded[q][f["donors"][method]]["completion_sha256"]
                h = A.write_json(OUT / "forecasts.json", forecasts)
                A.write_json(OUT / "FORECASTS_SEALED.json", {"utc": A.utc(), "forecasts_sha256": h,
                                                           "queries": len(forecasts),
                                                           "fold_info_sha256": A.sha_file(OUT / "fold_info.json"),
                                                           "run_started_sha256": A.sha_file(OUT / "RUN_STARTED.json")}, exclusive=True)
                log(f"FORECASTS_SEALED: {len(forecasts)} queries")
            evaluation_admitted(OUT)
            forecasts = A.read_json(OUT / "forecasts.json"); fold_info = A.read_json(OUT / "fold_info.json")
            scores = load_evaluation_labels(lib, OUT)
            jobs = []
            needed = 0
            for qid, f in sorted(forecasts.items()):
                q = lib.index[qid]
                missing = sorted({d for d in f["donors"].values() if d is not None and not scores.has(q, lib.index[d])})
                if missing:
                    needed += len(missing)
                    jobs.append({"job": f"missing|{qid}", "query": qid, "obs_sha": f["obs_sha"], "donors": missing})
            A.write_json(OUT / "cache_coverage.json", {"new_missing_pairs": needed, "query_jobs": len(jobs),
                                                      "admitted_after_forecast_seal": True})
            A.run_jobs(missing_label_job, jobs, OUT / "labels_missing.jsonl", pool, log, "missing eval labels", store=scores)
        result = summarise(lib, forecasts, scores, fold_info, OUT)
        verify_inputs()
        A.write_json(OUT / "RUN_COMPLETE.json", {"utc": A.utc(), "results_sha256": A.sha_file(OUT / "results.json"),
                                               "forecasts_sha256": A.sha_file(OUT / "forecasts.json")}, exclusive=True)
        log(result["verdict"]["sentence"])
    except A.leak_types() as error:
        A.write_json(OUT / "RUN_FAILED.json", {"utc": A.utc(), "reason": str(error)}, exclusive=True)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("plan", "run")); parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args(argv)
    if args.command == "plan":
        print(json.dumps(plan_only(), indent=2))
    else:
        run(args.workers)


if __name__ == "__main__":
    main()
