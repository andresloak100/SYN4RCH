"""Feasibility for protocol 2 (reports/next_donor_selector_protocol.md): a full-TRAIN donor selector.

TRAINING DATA ONLY. Nothing here predicts, admits, seals, scores or freezes anything of the current
evaluation. It reads:
  * reports/axon_dataset_manifest.tsv, rows with part == "train" only (animal, source label, soma label);
  * data/interim/axon_dataset/local_manifest.tsv and the soma table, through `known_query_inputs`, which skips
    every row whose neuron is not in the requested parts before parsing any of its fields (here: "train" only);
  * data/interim/axon_engine/real/ingestion_log.json (training traces only: table soma, trace soma, status);
  * canonical training arbors under data/interim/axon_engine/real/library/ (read-only, a sample);
  * the 37-trace snapshot data/interim/axon_engine/real_local37/{library,paircache} (complete cache, read-only);
  * the real all-pairs training cache (memory-mapped, read-only) only when its COMPLETE.json exists and its
    D/F/index hashes and every arbor hash verify; every read of it goes through `FitCache`, which refuses any
    neuron outside the fit animals of the current fit call.
No validation or test trace, no regional table at all, nothing under data/raw/jiao2025/traces/.

Subcommands (each writes data/interim/next_donor_selector/feasibility_<subcommand>.json)
  cache-meaning   what the native cache D/F means versus the label the referee would score; per-pair timing of
                  the frozen metric, the referee metric, arbor loading and the attachment ladder
  folds           whole-animal nested folds (outer, inner, final CV), pair counts including the fallback pairs,
                  both-endpoints-fit-only proof, fit-only cache reads
  known-inputs    the known-input loader on the training rows, checked against the ingestion log
  mdd             minimal detectable paired animal difference from exposed development deltas
  all             all of the above plus the compute projection
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/interim/next_donor_selector"
MANIFEST = ROOT / "reports/axon_dataset_manifest.tsv"
REAL = ROOT / "data/interim/axon_engine/real"
LOCAL37 = ROOT / "data/interim/axon_engine/real_local37"
CODEX_RUNS = ROOT / "data/interim/training_donor_selector_development"
SALT = "next_donor_selector|20261111"
OUTER_FOLDS, INNER_FOLDS = 5, 4
CANDIDATE_M = (16, 64)                      # candidate-set sizes the inner CV may choose between (protocol section 6)
ATTACH_LADDER = (500.0, 1000.0)             # the frozen ladder, used only to show what the cache does NOT hold
THREADS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
FINAL_FOLDS = 5                             # the final whole-animal CV on all 101 animals (protocol section 7)
LOCAL_MANIFEST = ROOT / "data/interim/axon_dataset/local_manifest.tsv"
SOMA_TABLE = ROOT / "data/raw/jiao2025/All_soma_infor_HY_neurons.csv"
# pins copied from analysis/axon_engine.py (PLAN_MANIFEST_SHA, SOMA_TABLE_SHA, MIDLINE_ML); checked against its text
PLAN_MANIFEST_SHA = "649f829a592cd3d13bcda5cffdbd049e759568d38bacb144f6a61220610a1ca4"
SOMA_TABLE_SHA = "cc01d048bf5a7c9c9f46869fb750f2d67a98ca934903d2e2e71c1a67b1d97cb9"
MIDLINE_ML = 5700.0
SOMA_COLUMNS_USED = ("file", "neurons", "geno", "AP", "DV", "ML")     # never 'type', 'name' or 'location'
KNOWN_FIELDS = ("id", "animal", "soma_um", "peptide", "nucleus", "hemisphere", "reconstruction_type")
REFEREE_KNOWN_FIELDS = ("id", "animal", "soma_um", "peptide", "nucleus")   # the referee's known(); the rest is dropped
FORBIDDEN_INPUT_NAMES = ("projection_length", "regional")                  # the regional table is never opened


def _module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def frozen_modules():
    """The frozen template (and through it d110_geometry) and the referee's metric module, read-only."""
    T = _module("nds_frozen_template", ROOT / "analysis/brm_template.py")
    R = _module("nds_referee_geometry", ROOT / "analysis/brm_geometry_audit.py")
    return T, T.G, R


def sha_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(x):
    return json.dumps(x, sort_keys=True, separators=(",", ":"), allow_nan=False)


# --------------------------------------------------------------------------------------------
# pure helpers (unit tested)
# --------------------------------------------------------------------------------------------

def assign_folds(animals, count, salt):
    """Whole-animal assignment by a salted hash order; every animal in exactly one fold."""
    ordered = sorted(set(animals), key=lambda a: (hashlib.sha256(f"{SALT}|{salt}|{a}".encode()).hexdigest(), a))
    if len(ordered) < count:
        raise ValueError("fewer animals than folds")
    return {a: i % count for i, a in enumerate(ordered)}


def nested_folds(animals, outer=OUTER_FOLDS, inner=INNER_FOLDS):
    animals = sorted(set(animals))
    o = assign_folds(animals, outer, "outer")
    folds = []
    for k in range(outer):
        fit = [a for a in animals if o[a] != k]
        evaluation = [a for a in animals if o[a] == k]
        i = assign_folds(fit, inner, f"inner|{k}")
        folds.append({"outer": k, "fit_animals": fit, "evaluation_animals": evaluation,
                      "inner": [{"inner": j, "fit_animals": [a for a in fit if i[a] != j],
                                 "evaluation_animals": [a for a in fit if i[a] == j]} for j in range(inner)]})
    return folds


def check_fit_pairs(pairs, animal_of, fit_animals, forbidden_animals):
    """Every training pair has BOTH endpoints in the fit animals, never in a forbidden (evaluation) animal,
    and never the same animal twice. Raises on the first violation; returns the number checked."""
    fit, bad = set(fit_animals), set(forbidden_animals)
    if fit & bad:
        raise ValueError("an animal is both fit-only and evaluation")
    for q, d in pairs:
        aq, ad = animal_of[q], animal_of[d]
        if aq not in fit or ad not in fit:
            raise ValueError(f"pair ({q},{d}) has an endpoint outside the fit animals")
        if aq in bad or ad in bad:
            raise ValueError(f"pair ({q},{d}) touches an evaluation animal")
        if aq == ad:
            raise ValueError(f"pair ({q},{d}) is within one animal")
    return len(pairs)


def candidate_pairs(queries, pool, animal_of, soma, m):
    """Directed (query, donor) pairs: each query with its m nearest POOL donors of other animals by table soma
    (ties by donor ID). For fit labels queries == pool; for evaluation queries the pool is the fit pool only."""
    from scipy.spatial import cKDTree
    pool = sorted(pool)
    xyz = np.array([soma[n] for n in pool], float)
    tree = cKDTree(xyz)
    per_animal = Counter(animal_of[n] for n in pool)
    out = []
    for q in sorted(queries):
        k = min(len(pool), m + per_animal.get(animal_of[q], 0) + 1)
        while True:
            dist, idx = tree.query(np.asarray(soma[q], float), k=k)
            dist, idx = np.atleast_1d(dist), np.atleast_1d(idx)
            keep = [(d, pool[j]) for d, j in zip(dist, idx) if animal_of[pool[j]] != animal_of[q]]
            if len(keep) >= m or k >= len(pool):
                break
            k = min(len(pool), 2 * k)
        keep.sort()
        out += [(q, d) for _, d in keep[:m]]
    return out


def fallback_choice(selected, fallback, predicted):
    """Copying fallback, fixed now: keep the selector's donor only when its predicted distance is lower AND
    its predicted F1 is higher than the predicted values of the fallback donor (the referee 'medoid' of the
    same pool) for this query; otherwise copy the fallback donor. `predicted` maps donor -> (distance, f1)."""
    if selected == fallback:
        return fallback, "selector_equals_fallback"
    ds, fs = predicted[selected]
    df, ff = predicted[fallback]
    if not all(np.isfinite([ds, fs, df, ff])):
        return fallback, "nonfinite_prediction"
    if ds < df and fs > ff:
        return selected, "selector_predicted_to_dominate"
    return fallback, "not_predicted_to_dominate_fallback"


def final_folds(animals, count=FINAL_FOLDS):
    """The final whole-animal CV on all training animals that chooses the configuration of the final fit."""
    animals = sorted(set(animals))
    f = assign_folds(animals, count, "final")
    return [{"final": k, "fit_animals": [a for a in animals if f[a] != k],
             "evaluation_animals": [a for a in animals if f[a] == k]} for k in range(count)]


class FitCache:
    """The only door to the all-pairs training cache during a fit call. It refuses (raises) any neuron whose
    animal is not in this call's fit animals, so no centrality, medoid, kernel, local medoid or label can read an
    evaluation animal (outer or inner). It counts what it served and what it refused."""

    def __init__(self, D, F, ids, animal_of, fit_animals):
        self.D, self.F = D, F
        self.pos = {n: i for i, n in enumerate(ids)}
        self.animal_of = animal_of
        self.fit = frozenset(fit_animals)
        self.reads, self.refused = 0, 0

    def _index(self, names):
        bad = [n for n in names if self.animal_of[n] not in self.fit]
        if bad:
            self.refused += len(bad)
            raise ValueError(f"cache read refused: {len(bad)} neuron(s) outside the fit animals, e.g. {bad[0]}")
        return np.array([self.pos[n] for n in names], dtype=np.intp)

    def block(self, rows, cols, which="D"):
        r, c = self._index(rows), self._index(cols)
        M = self.D if which == "D" else self.F
        out = np.asarray(M[np.ix_(r, c)], dtype=float)
        self.reads += out.size
        return out

    def pairs(self, pairs, which="D"):
        r = self._index([q for q, _ in pairs]); c = self._index([d for _, d in pairs])
        M = self.D if which == "D" else self.F
        out = np.array([M[i, j] for i, j in zip(r, c)], dtype=float)
        self.reads += out.size
        return out


def referee_medoid(pairwise, indices, animals, medoid_fn):
    """The referee's arithmetic (brm_geometry_audit.geometry_medoid): animal-balanced weights, argmin of
    D[ix, ix] @ w. Returns a position within `indices`."""
    return int(medoid_fn(pairwise, list(indices), list(animals)))


def fallback_donors(cache, pool, queries, medoid_fn):
    """Fallback donor of each query: the referee `medoid` of the pool WITHOUT the query's animal (the referee never
    lets a query's animal into its donors). For queries outside the pool's animals this is the pool medoid; for fit
    queries it is the leave-own-animal-out medoid. All reads go through the fit-only cache."""
    pool = sorted(pool)
    Dp = cache.block(pool, pool)
    animals = [cache.animal_of[n] for n in pool]
    memo, out = {}, {}
    for q in sorted(queries):
        a = cache.animal_of[q]
        if a not in memo:
            ix = [i for i, x in enumerate(animals) if x != a]
            memo[a] = pool[ix[referee_medoid(Dp, ix, [animals[i] for i in ix], medoid_fn)]]
        out[q] = memo[a]
    return out


def candidate_sets(pairs, fallback):
    """Candidate set of each query: its M nearest other-animal donors PLUS its fallback donor, so the fallback
    comparison is between two donors the model was trained on (protocol section 6)."""
    sets = {}
    for q, d in pairs:
        sets.setdefault(q, []).append(d)
    for q, f in fallback.items():
        if q in sets and f not in sets[q]:
            sets[q].append(f)
    return sets


def local_medoid(cache, candidates, medoid_fn):
    """Matched control for the chosen M: the referee-arithmetic, animal-balanced medoid of the M nearest donors
    (the candidate set without the added fallback donor), cache D, fit-only."""
    D = cache.block(candidates, candidates)
    return candidates[referee_medoid(D, range(len(candidates)), [cache.animal_of[n] for n in candidates], medoid_fn)]


def stable_seed(*parts):
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:4], "big")


def matched_switch_controls(final_choice, fallback, nearest, candidates, seed=20261111):
    """Matched controls for the switch: on exactly the queries where the combined rule left the fallback, put a
    fixed-rule donor from the same candidate set instead ('switch_nearest': the soma-nearest candidate other than
    the fallback; 'switch_random': a seeded random candidate other than the fallback); elsewhere the fallback.
    Returns (switch_nearest, switch_random, coverage)."""
    sn, sr = {}, {}
    switched = sorted(q for q in final_choice if final_choice[q] != fallback[q])
    for q in final_choice:
        if final_choice[q] == fallback[q]:
            sn[q] = sr[q] = fallback[q]
            continue
        others = [d for d in candidates[q] if d != fallback[q]]
        sn[q] = nearest[q] if nearest[q] != fallback[q] else others[0]
        rng = np.random.default_rng(stable_seed(seed, "switch_random", q))
        sr[q] = sorted(others)[int(rng.integers(len(others)))]
    return sn, sr, {"queries": len(final_choice), "switched": len(switched),
                    "switched_fraction": len(switched) / len(final_choice) if final_choice else float("nan")}


def verdict_sentence(deltas, intervals):
    """deltas: control -> (distance delta um, F1 delta), selector minus control; intervals: control ->
    ((lo, hi) distance, (lo, hi) F1), paired animal bootstrap 95%, descriptive. The decision is the point gate;
    the sentence says whether every interval also excludes 0 in the favourable direction."""
    met = all(dd < 0 and df > 0 for dd, df in deltas.values())
    if not met:
        return "gate not met"
    resolved = all(intervals[c][0][1] < 0 and intervals[c][1][0] > 0 for c in deltas)
    return "gate met and resolved against every control" if resolved else "gate met, not resolved"


def known_query_inputs(parts, manifest=MANIFEST, local_manifest=LOCAL_MANIFEST, soma_table=SOMA_TABLE,
                       pins=True):
    """The ONLY source of a query's inputs (training or held-out): plan manifest (animal, part, source and soma
    labels), local manifest (hemisphere and reconstruction type, cross-checked), soma table restricted to
    SOMA_COLUMNS_USED (table soma, mirrored-left frame). Rows whose neuron is not in `parts` are skipped before any
    of their fields is parsed. The regional table is never opened. Output fields: KNOWN_FIELDS only."""
    parts = set(parts)
    for path in (manifest, local_manifest, soma_table):
        if any(word in Path(path).name for word in FORBIDDEN_INPUT_NAMES):
            raise ValueError(f"{path}: not a permitted input file")
    if pins:
        if sha_file(manifest) != PLAN_MANIFEST_SHA or sha_file(soma_table) != SOMA_TABLE_SHA:
            raise ValueError("manifest or soma table differs from the engine's pin")
    with open(manifest, newline="") as f:
        plan = {r["neuron"]: r for r in csv.DictReader(f, delimiter="\t") if r["part"] in parts}
    with open(local_manifest, newline="") as f:
        local = {r["neuron"]: r for r in csv.DictReader(f, delimiter="\t") if r["neuron"] in plan}
    if set(local) != set(plan):
        raise ValueError("plan and local manifest do not list the same neurons")
    soma = {}
    with open(soma_table, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        col = {name: header.index(name) for name in SOMA_COLUMNS_USED}
        for row in reader:
            n = row[col["neurons"]]
            if n not in plan:
                continue                                    # not requested: nothing else of the row is parsed
            if n in soma:
                raise ValueError(f"duplicate soma row {n}")
            soma[n] = (row[col["file"]], row[col["geno"]], tuple(float(row[col[k]]) for k in ("AP", "DV", "ML")))
    if set(soma) != set(plan):
        raise ValueError("soma table does not list every requested neuron")
    out = {}
    for n in sorted(plan):
        r, l = plan[n], local[n]
        for key in ("animal", "source_label", "soma_label", "part"):
            if l[key] != r[key]:
                raise ValueError(f"local manifest differs from the plan manifest for {n}: {key}")
        file_id, geno, xyz = soma[n]
        if not file_id.startswith(n + "_") or ("Orexin" if geno == "Wt(orexin)" else geno) != r["source_label"]:
            raise ValueError(f"identity or label disagrees for {n}")
        if not all(np.isfinite(xyz)) or xyz[2] > MIDLINE_ML + 1e-6:
            raise ValueError(f"table soma of {n} is not in the mirrored-left frame")
        out[n] = {"id": n, "animal": r["animal"], "soma_um": np.array(xyz, float),
                  "peptide": r["source_label"] or None, "nucleus": r["soma_label"] or None,
                  "hemisphere": l["hemisphere"], "reconstruction_type": l["reconstruction_type"]}
    return out


def referee_row(known):
    """What is handed to the referee: its known() fields only."""
    return {k: known[k] for k in REFEREE_KNOWN_FIELDS}


def mdd(sd, n, alpha=0.05, power=0.8):
    """Paired whole-animal minimal detectable mean difference (normal approximation, two-sided)."""
    from scipy.stats import norm
    return float((norm.ppf(1 - alpha / 2) + norm.ppf(power)) * sd / np.sqrt(n))


def ci_halfwidth(sd, n):
    return float(1.959963984540054 * sd / np.sqrt(n))


# --------------------------------------------------------------------------------------------
# training metadata
# --------------------------------------------------------------------------------------------

def training_metadata():
    with MANIFEST.open() as f:
        rows = {r["neuron"]: r for r in csv.DictReader(f, delimiter="\t") if r["part"] == "train"}
    log = json.loads((REAL / "ingestion_log.json").read_text())["neurons"]
    out = {}
    for n, e in log.items():
        if e.get("part") != "train" or n not in rows:
            raise ValueError(f"{n}: a non-training neuron is in the training ingestion log")
        if e["status"] != "included":
            continue
        r = rows[n]
        if r["animal"] != e["animal"]:
            raise ValueError(f"{n}: animal disagrees between manifest and ingestion log")
        out[n] = {"id": n, "animal": r["animal"], "peptide": r["source_label"] or None, "nucleus": r["soma_label"] or None,
                  "table_soma": np.array(e["table_soma_um"], float), "trace_soma": np.array(e["trace_soma_um"], float),
                  "offset_um": float(e["soma_table_distance_um"]), "arbor_sha256": e["arbor_sha256"]}
    return out, {"manifest_training_rows": len(rows), "ingestion_listed": len(log),
                 "included": len(out), "excluded": len(log) - len(out),
                 "animals": len({v["animal"] for v in out.values()})}


def load_real_cache(meta):
    """The real all-pairs training cache, read-only and memory-mapped, only when COMPLETE.json exists and the D, F
    and index hashes and every arbor hash (against the training ingestion log) verify. None otherwise."""
    done = sorted((REAL / "paircache").glob("*/COMPLETE.json"))
    if not done:
        return None
    cdir_r = done[0].parent
    rec = json.loads((cdir_r / "COMPLETE.json").read_text())
    idx = json.loads((cdir_r / "index.json").read_text())
    if sha_file(cdir_r / "index.json") != rec["index_sha256"]:
        raise ValueError("real cache index differs from its COMPLETE record")
    for name in ("D.npy", "F.npy"):
        if sha_file(cdir_r / name) != rec[name[0] + "_sha256"]:
            raise ValueError(f"real cache {name} differs from its COMPLETE record")
    if {n: h for n, h in zip(idx["ids"], idx["arbor_sha256"])} != {n: meta[n]["arbor_sha256"] for n in meta}:
        raise ValueError("real cache arbors differ from the training ingestion log")
    return {"dir": str(cdir_r.relative_to(ROOT)), "cache_sha256": rec["cache_sha256"], "ids": list(idx["ids"]),
            "pos": {n: i for i, n in enumerate(idx["ids"])},
            "D": np.load(cdir_r / "D.npy", mmap_mode="r"), "F": np.load(cdir_r / "F.npy", mmap_mode="r")}


# --------------------------------------------------------------------------------------------
# 1. cache meaning
# --------------------------------------------------------------------------------------------

def _attach(T, arbor, soma):
    for radius in ATTACH_LADDER:
        try:
            return T.compact_attachment(arbor, soma, radius)[0], f"local_attachment_{int(radius)}"
        except T.UnsupportedAttachment:
            continue
    return arbor.shifted(soma), "translation"


def cache_meaning(max_attach_pairs=60, real_sample_pairs=40, stress_queries=6, stress_donors=5, attach_timing_pairs=20,
                  seed=20261111):
    T, G, R = frozen_modules()
    report = {"real_cache_complete": False}
    real_caches = sorted((REAL / "paircache").glob("*/COMPLETE.json"))
    report["real_cache_complete"] = bool(real_caches)
    report["real_cache_note"] = ("COMPLETE.json present: its D/F hashes are verified and the sampled real pairs below are "
                                 "compared with it (memory-mapped, read-only)"
                                 if real_caches else "no COMPLETE.json under data/interim/axon_engine/real/paircache: "
                                 "the real cache was NOT read; the check below uses the complete 37-trace snapshot cache")
    cdir = next(p.parent for p in sorted((LOCAL37 / "paircache").glob("*/COMPLETE.json")))
    index = json.loads((cdir / "index.json").read_text())
    complete = json.loads((cdir / "COMPLETE.json").read_text())
    for name in ("D.npy", "F.npy"):
        if sha_file(cdir / name) != complete[name[0] + "_sha256"]:
            raise ValueError(f"local37 cache {name} differs from its COMPLETE record")
    D, F = np.load(cdir / "D.npy"), np.load(cdir / "F.npy")
    ids = index["ids"]
    log37 = json.loads((LOCAL37 / "ingestion_log.json").read_text())["neurons"]
    if any(log37[n]["part"] != "train" for n in ids):
        raise ValueError("a non-training neuron in the local37 cache")
    arb = {n: G.load_arbor(LOCAL37 / "library" / f"{n}.npz") for n in ids}
    for n in ids:
        if not np.allclose(arb[n].xyz[0], log37[n]["trace_soma_um"], atol=1e-9):
            raise ValueError(f"{n}: library root is not the trace soma")
    table = {n: np.array(log37[n]["table_soma_um"], float) for n in ids}
    animal = {n: log37[n]["animal"] for n in ids}
    cen = {n: arb[n].shifted(np.zeros(3)) for n in ids}

    def tup(a):
        return (a.xyz, a.parent, a.types)

    rows = []
    for i, q in enumerate(ids):
        for j, d in enumerate(ids):
            if animal[q] == animal[d]:
                continue
            lo, hi = (q, d) if q < d else (d, q)
            a = G.geometry_metrics(cen[lo], cen[hi])                                   # cache definition
            r = R.geometry_metrics(R.shifted(tup(arb[d]), np.zeros(3)), R.shifted(tup(arb[q]), np.zeros(3)))
            lab = R.geometry_metrics(R.shifted(tup(arb[d]), table[q]), tup(arb[q]))    # referee label: translated to table soma
            trc = R.geometry_metrics(R.shifted(tup(arb[d]), arb[q].xyz[0]), tup(arb[q]))  # translated to TRACE soma
            rows.append({"query": q, "donor": d, "D": float(D[i, j]), "F": float(F[i, j]),
                         "cacheformula_D": a["distance_um"], "cacheformula_F": 1 - a["f1_200"],
                         "referee_centred_D": r["distance_um"], "referee_centred_F": 1 - r["f1_200"],
                         "label_D": lab["distance_um"], "label_F": 1 - lab["f1_200"],
                         "trace_D": trc["distance_um"], "trace_F": 1 - trc["f1_200"],
                         "offset_um": float(np.linalg.norm(table[q] - arb[q].xyz[0]))})
    k = len(rows)

    def stat(x, y):
        diff = np.abs(np.array([r[x] for r in rows]) - np.array([r[y] for r in rows]))
        return {"max_abs": float(diff.max()), "mean_abs": float(diff.mean()), "exact_equal": int((diff == 0).sum()),
                "within_1e-9": int((diff <= 1e-9).sum()), "pairs": k}
    agree = {"cache_vs_its_formula_D": stat("D", "cacheformula_D"), "cache_vs_its_formula_F": stat("F", "cacheformula_F"),
             "cache_vs_referee_centred_D": stat("D", "referee_centred_D"), "cache_vs_referee_centred_F": stat("F", "referee_centred_F"),
             "cache_vs_translated_to_trace_soma_D": stat("D", "trace_D"), "cache_vs_translated_to_trace_soma_F": stat("F", "trace_F"),
             "cache_vs_referee_label_table_soma_D": stat("D", "label_D"), "cache_vs_referee_label_table_soma_F": stat("F", "label_F")}
    excess = [abs(r["label_D"] - r["D"]) - r["offset_um"] for r in rows]
    agree["distance_bound_|label-D|<=offset_violations"] = int(sum(e > 1e-9 for e in excess))
    agree["f1_changed_pairs"] = int(sum(abs(r["label_F"] - r["F"]) > 1e-12 for r in rows))
    # does the cache rank donors like the label does? (per query, among its cross-animal donors)
    from scipy.stats import spearmanr
    top_same, rhos, regret = 0, [], []
    for q in ids:
        mine = [r for r in rows if r["query"] == q]
        cD, lD = np.array([r["D"] for r in mine]), np.array([r["label_D"] for r in mine])
        top_same += int(np.argmin(cD) == np.argmin(lD))
        rhos.append(float(spearmanr(cD, lD)[0]))
        regret.append(float(lD[np.argmin(cD)] - lD.min()))
    agree["per_query_top1_same"] = f"{top_same}/{len(ids)}"
    agree["per_query_spearman_min"] = min(rhos)
    agree["per_query_regret_um_max"] = max(regret)
    # attachment (frozen ladder) versus the cache: a different geometry
    rng = np.random.default_rng(seed)
    sample = [rows[i] for i in rng.choice(k, size=min(max_attach_pairs, k), replace=False)]
    att = []
    for r in sample:
        out, how = _attach(T, arb[r["donor"]], table[r["query"]])
        m = R.geometry_metrics(tup(out), tup(arb[r["query"]]))
        att.append({"how": how, "abs_diff_D": abs(m["distance_um"] - r["D"]), "abs_diff_F": abs(1 - m["f1_200"] - r["F"]),
                    "donor_query_soma_um": float(np.linalg.norm(arb[r["donor"]].xyz[0] - table[r["query"]]))})
    by = Counter(a["how"] for a in att)
    agree["attachment_ladder_sample"] = {
        "pairs": len(att), "executed": dict(by),
        **{h: {"max_abs_diff_D_um": max(a["abs_diff_D"] for a in att if a["how"] == h),
               "median_abs_diff_D_um": float(np.median([a["abs_diff_D"] for a in att if a["how"] == h])),
               "max_abs_diff_F": max(a["abs_diff_F"] for a in att if a["how"] == h)} for h in by}}
    report["local37"] = {"cache_dir": str(cdir.relative_to(ROOT)), "cache_sha256": complete["cache_sha256"],
                         "neurons": len(ids), "animals": len(set(animal.values())),
                         "directed_cross_animal_pairs": k, "agreement": agree,
                         "offsets_um": {"max": max(r["offset_um"] for r in rows), "median": float(np.median([r["offset_um"] for r in rows]))}}
    # real library: a seeded random sample of cross-animal training pairs plus a stress sample (the queries with the
    # largest table-to-trace soma offsets, each with its nearest-soma donors), recomputed directly; compared with the
    # real cache only when its COMPLETE.json exists and its hashes verify
    meta, _ = training_metadata()
    rid = sorted(meta)
    picks, seen = [], set()
    while len(picks) < real_sample_pairs:
        q, d = rng.choice(len(rid), 2, replace=False)
        q, d = rid[q], rid[d]
        if meta[q]["animal"] != meta[d]["animal"] and (q, d) not in seen:
            seen.add((q, d)); picks.append((q, d, "random"))
    worst = sorted(rid, key=lambda n: -meta[n]["offset_um"])[:stress_queries]
    for q in worst:
        near = sorted((float(np.linalg.norm(meta[d]["table_soma"] - meta[q]["table_soma"])), d) for d in rid
                      if meta[d]["animal"] != meta[q]["animal"])[:stress_donors]
        picks += [(q, d, "stress_large_offset_nearest_donors") for _, d in near]
    cache = load_real_cache(meta) if real_caches else None
    load, load_times = {}, []
    def arbor(n):
        if n not in load:
            t0 = time.perf_counter()
            a = G.load_arbor(REAL / "library" / f"{n}.npz")
            load_times.append(time.perf_counter() - t0)
            if not np.allclose(a.xyz[0], meta[n]["trace_soma"], atol=1e-9):
                raise ValueError(f"{n}: library root is not the trace soma")
            load[n] = a
        return load[n]
    real, times, rtimes, atimes, atimes_ok, ahow = [], [], [], [], [], Counter()
    for t, (q, d, kind) in enumerate(picks):
        aq, ad = arbor(q), arbor(d)
        lo, hi = (aq, ad) if q < d else (ad, aq)
        t0 = time.perf_counter()
        c = G.geometry_metrics(lo.shifted(np.zeros(3)), hi.shifted(np.zeros(3)))
        times.append(time.perf_counter() - t0)
        t0 = time.perf_counter()
        lab = R.geometry_metrics(R.shifted(tup(ad), meta[q]["table_soma"]), tup(aq))
        rtimes.append(time.perf_counter() - t0)
        if t < attach_timing_pairs or kind != "random":  # cost of the frozen attachment ladder (C4, C5, E1, E2);
            t0 = time.perf_counter()                       # the stress pairs are NEAREST donors, where it can attach
            _, how = _attach(T, ad, meta[q]["table_soma"])
            dt = time.perf_counter() - t0
            atimes.append(dt); ahow[how] += 1
            if how != "translation":
                atimes_ok.append(dt)
        row = {"kind": kind, "dD": abs(lab["distance_um"] - c["distance_um"]), "dF": abs(lab["f1_200"] - c["f1_200"]),
               "offset": meta[q]["offset_um"], "nodes": len(aq.xyz) + len(ad.xyz)}
        if cache is not None:
            i, j = cache["pos"][q], cache["pos"][d]
            row["cache_equals_formula"] = bool(cache["D"][i, j] == c["distance_um"] and cache["F"][i, j] == 1 - c["f1_200"])
            row["cache_vs_label_dD"] = abs(float(cache["D"][i, j]) - lab["distance_um"])
            row["cache_vs_label_dF1"] = abs(1 - float(cache["F"][i, j]) - lab["f1_200"])
        real.append(row)
    summary = {}
    for kind in sorted({r["kind"] for r in real}):
        mine = [r for r in real if r["kind"] == kind]
        summary[kind] = {"pairs": len(mine), "max_abs_diff_D_um": max(r["dD"] for r in mine),
                         "max_abs_diff_F1": max(r["dF"] for r in mine),
                         "pairs_with_F1_changed": sum(r["dF"] > 1e-12 for r in mine),
                         "bound_violations": sum(r["dD"] > r["offset"] + 1e-9 for r in mine),
                         "max_query_offset_um": max(r["offset"] for r in mine)}
        if cache is not None:
            summary[kind].update({"cache_bit_equal_to_formula": f"{sum(r['cache_equals_formula'] for r in mine)}/{len(mine)}",
                                  "cache_vs_label_max_dD_um": max(r["cache_vs_label_dD"] for r in mine),
                                  "cache_vs_label_max_dF1": max(r["cache_vs_label_dF1"] for r in mine)})
    report["real_library_sample"] = {
        "note": "computed directly from data/interim/axon_engine/real/library",
        "cache": None if cache is None else {"dir": cache["dir"], "cache_sha256": cache["cache_sha256"],
                                             "hashes_verified": True, "arbor_hashes_match_ingestion_log": True},
        "by_kind": summary,
        "ms_per_pair_d110_centred_median": 1e3 * float(np.median(times)), "ms_per_pair_d110_centred_mean": 1e3 * float(np.mean(times)),
        "ms_per_pair_referee_label_median": 1e3 * float(np.median(rtimes)), "ms_per_pair_referee_label_mean": 1e3 * float(np.mean(rtimes)),
        "ms_per_arbor_load_mean": 1e3 * float(np.mean(load_times)), "arbors_loaded": len(load_times),
        "attachment_ladder": {"calls": len(atimes), "ms_mean": 1e3 * float(np.mean(atimes)),
                              "ms_max": 1e3 * float(np.max(atimes)), "executed": dict(ahow),
                              "ms_mean_when_attached": 1e3 * float(np.mean(atimes_ok)) if atimes_ok else None,
                              "ms_max_when_attached": 1e3 * float(np.max(atimes_ok)) if atimes_ok else None},
        "timing_note": ("one process, wall-clock per call; machine load moves these by about 10 percent between runs")}
    off = np.array([v["offset_um"] for v in meta.values()])
    report["real_training_table_trace_offsets_um"] = {
        "neurons": int(len(off)), "max": float(off.max()), "median": float(np.median(off)), "p95": float(np.quantile(off, .95)),
        "exactly_zero": int((off == 0).sum()), "below_0.01": int((off < .01).sum()), "above_1": int((off > 1).sum())}
    return report


# --------------------------------------------------------------------------------------------
# 2. folds and pair counts
# --------------------------------------------------------------------------------------------

def fit_configurations(animals):
    """Every fit call of the protocol: (name, fit animals, forbidden animals, evaluation animals). Outer folds,
    their inner folds (forbidden = inner AND outer evaluation animals), the final CV on all animals, and the final
    fit on all animals (no training animal is forbidden; held-out animals are never in the training cache)."""
    configs = []
    for f in nested_folds(animals):
        configs.append((f"outer{f['outer']}", f["fit_animals"], f["evaluation_animals"], f["evaluation_animals"]))
        for x in f["inner"]:
            configs.append((f"outer{f['outer']}_inner{x['inner']}", x["fit_animals"],
                            x["evaluation_animals"] + f["evaluation_animals"], x["evaluation_animals"]))
    for f in final_folds(animals):
        configs.append((f"final{f['final']}", f["fit_animals"], f["evaluation_animals"], f["evaluation_animals"]))
    configs.append(("final_all", sorted(set(animals)), [], []))
    return configs


def folds_report(cache=None):
    meta, counts = training_metadata()
    ids = sorted(meta)
    animal_of = {n: meta[n]["animal"] for n in ids}
    animals = sorted(set(animal_of.values()))
    per_animal = Counter(animal_of.values())
    soma = {n: meta[n]["table_soma"] for n in ids}
    _, _, R = frozen_modules()

    def cross_pairs(group):
        n = sum(per_animal[a] for a in group)
        return n * (n - 1) - sum(per_animal[a] * (per_animal[a] - 1) for a in group)

    union = {m: set() for m in CANDIDATE_M}
    calls, fallback_pairs_total = [], 0
    t_start = time.perf_counter()
    for name, fit_animals, forbidden, evaluation in fit_configurations(animals):
        fit_set = set(fit_animals)
        group = [n for n in ids if animal_of[n] in fit_set]
        evq = [n for n in ids if animal_of[n] in set(evaluation)]
        rec = {"name": name, "fit_animals": len(fit_animals), "fit_neurons": len(group),
               "evaluation_animals": len(evaluation), "evaluation_neurons": len(evq),
               "forbidden_animals": len(forbidden), "full_cross_animal_fit_pairs": cross_pairs(fit_animals)}
        fb_fit = fb_eval = None
        if cache is not None:
            fc = FitCache(cache["D"], cache["F"], cache["ids"], animal_of, fit_animals)
            t0 = time.perf_counter()
            fb_fit = fallback_donors(fc, group, group, R.geometry_medoid)          # leave-own-animal-out medoid
            fb_eval = fallback_donors(fc, group, evq, R.geometry_medoid) if evq else {}
            fb_pairs = sorted(fb_fit.items())
            check_fit_pairs(fb_pairs, animal_of, fit_animals, forbidden)
            fit_labels = fc.pairs(fb_pairs, "D")
            if not np.isfinite(fit_labels).all():
                raise ValueError("a non-finite fallback label")
            refused_probe = 0
            if evq:                                                               # the door must refuse an evaluation neuron
                try:
                    fc.block([evq[0]], [group[0]])
                except ValueError:
                    refused_probe = 1
                if not refused_probe:
                    raise ValueError("FitCache served an evaluation neuron")
            rec.update({"fallback_seconds": time.perf_counter() - t0,
                        "distinct_fallback_donors_fit": len(set(fb_fit.values())),
                        "distinct_fallback_donors_eval": len(set(fb_eval.values())),
                        "fit_cache_reads": fc.reads, "fit_cache_refused_probe": refused_probe})
            fallback_pairs_total += len(fb_pairs) + len(fb_eval)
        for m in CANDIDATE_M:
            pairs = candidate_pairs(group, group, animal_of, soma, m)
            check_fit_pairs(pairs, animal_of, fit_animals, forbidden)
            union[m].update(pairs)
            score = candidate_pairs(evq, group, animal_of, soma, m) if evq else []
            if any(animal_of[d] not in fit_set for _, d in score):
                raise ValueError("a scoring donor outside the fit pool")
            union[m].update(score)
            if fb_fit is not None:
                sets = candidate_sets(pairs, fb_fit)
                added = sum(len(v) for v in sets.values()) - len(pairs)
                rec[f"M{m}_fallback_added_to_candidates"] = added
                union[m].update(fb_fit.items()); union[m].update(fb_eval.items())
        calls.append(rec)
    out = {"counts": counts, "outer_folds": OUTER_FOLDS, "inner_folds": INNER_FOLDS, "final_folds": FINAL_FOLDS,
           "neurons_per_animal": {"min": min(per_animal.values()), "median": float(np.median(list(per_animal.values()))),
                                  "max": max(per_animal.values())},
           "fit_calls": calls,
           "all_directed_cross_animal_pairs": cross_pairs(animals),
           "all_unordered_pairs_cache": len(ids) * (len(ids) - 1) // 2,
           "both_endpoints_fit_only": ("checked for every fit call (5 outer, 20 inner, 5 final-CV, final-all) by "
                                       "check_fit_pairs, label pairs and fallback pairs: 0 violations"),
           "fit_only_cache_reads": ("every cache read of every fit call went through FitCache; a probe read of an "
                                    "evaluation neuron was refused in every call that has evaluation animals"
                                    if cache is not None else "real cache not complete: not exercised"),
           "union_of_candidate_label_pairs": {str(m): len(union[m]) for m in CANDIDATE_M},
           "fallback_pairs_all_calls": fallback_pairs_total,
           "seconds": time.perf_counter() - t_start,
           "note_union": ("every (query, donor) pair any fit or evaluation step of any fit call could read: the M "
                          "nearest other-animal donors within the current pool plus the query's fallback donor")}
    return out


def known_inputs_report():
    """The known-input loader on the TRAINING rows only, checked against the ingestion log's table soma."""
    t0 = time.perf_counter()
    rows = known_query_inputs(("train",))
    meta, _ = training_metadata()
    exact = sum(np.array_equal(rows[n]["soma_um"], meta[n]["table_soma"]) for n in meta)
    labels = sum(rows[n]["peptide"] == meta[n]["peptide"] and rows[n]["nucleus"] == meta[n]["nucleus"] for n in meta)
    engine_text = (ROOT / "analysis/axon_engine.py").read_text()
    pins_match = all(x in engine_text for x in (PLAN_MANIFEST_SHA, SOMA_TABLE_SHA, f"MIDLINE_ML = {MIDLINE_ML}"))
    return {"parts_requested": ["train"], "rows": len(rows), "included_training_neurons_checked": len(meta),
            "table_soma_equal_to_ingestion_log": f"{exact}/{len(meta)}", "labels_equal": f"{labels}/{len(meta)}",
            "fields": list(KNOWN_FIELDS), "referee_fields": list(REFEREE_KNOWN_FIELDS),
            "hemisphere_values": dict(Counter(r["hemisphere"] for r in rows.values())),
            "reconstruction_type_values": dict(Counter(r["reconstruction_type"] for r in rows.values())),
            "pins_equal_to_engine_constants": pins_match, "regional_table_opened": False,
            "seconds": time.perf_counter() - t0}


# --------------------------------------------------------------------------------------------
# 3. minimal detectable difference from exposed development deltas (Codex's 37-neuron runs)
# --------------------------------------------------------------------------------------------

def mdd_report():
    out = {"source": [], "note": ("SDs of per-animal paired differences selector minus control, from Codex's development "
                                  "runs on 37 previously exposed training neurons; planning numbers only, not evidence")}
    sds = {"distance_um": [], "f1_200": []}
    for run in ("run_v1", "run_hgb_v2"):
        p = CODEX_RUNS / run / "results.json"
        r = json.loads(p.read_text())
        out["source"].append({"file": str(p.relative_to(ROOT)), "sha256": sha_file(p)})
        for ctrl in ("nearest", "medoid", "peptide_nearest", "peptide_medoid", "kernel16", "kernel_all"):
            diffs = r["against"][ctrl]["animal_differences"]
            dd = np.array([x["distance_delta_um"] for x in diffs]); ff = np.array([x["f1_delta"] for x in diffs])
            sds["distance_um"].append(float(dd.std(ddof=1))); sds["f1_200"].append(float(ff.std(ddof=1)))
        out["development_animals"] = len(diffs)
    plan = {k: {"median_sd": float(np.median(v)), "max_sd": float(np.max(v))} for k, v in sds.items()}
    out["sd_of_animal_differences"] = plan
    out["mdd_80pct_power_alpha05"] = {
        str(n): {k: {"median_sd": mdd(plan[k]["median_sd"], n), "max_sd": mdd(plan[k]["max_sd"], n)} for k in plan}
        for n in (101, 65, 40, 25)}
    out["ci95_halfwidth"] = {str(n): {k: ci_halfwidth(plan[k]["median_sd"], n) for k in plan} for n in (101, 65, 40, 25)}
    return out


def timing_projection(ms_per_pair, union_pairs, workers=2):
    return {"ms_per_pair": ms_per_pair, "workers": workers,
            **{f"M{m}": {"pairs": n, "processor_hours": n * ms_per_pair / 3.6e6,
                         "wall_hours": n * ms_per_pair / 3.6e6 / workers} for m, n in union_pairs.items()}}


OUTER_ARRAYS_PER_QUERY = 17     # selector+fallback, 4 referee copies, C1, C4, C5, E1, E2, kernel16, kernel_all,
                                # local_medoid_M, switch_nearest, switch_random, permuted-label selector, oracle (descr.)
ATTACH_CALLS_PER_QUERY = 6      # C4 (500), C5 (1000), E1 and E2 ladders (500 then 1000): upper bound
CONFIGURATIONS = 6


def scoring_projection(folds, ms_metric, ms_attach, ms_load, workers=2):
    """Upper bounds of the exact recomputation that scoring needs (every scored number is recomputed by the
    referee metric at the table soma), and of the attachment calls of the attached controls."""
    n = folds["counts"]["included"]
    calls = {c["name"]: c for c in folds["fit_calls"]}
    inner = sum(CONFIGURATIONS * c["evaluation_neurons"] for k, c in calls.items() if "_inner" in k)
    final_cv = sum(CONFIGURATIONS * c["evaluation_neurons"] for k, c in calls.items() if k.startswith("final") and k != "final_all")
    outer = OUTER_ARRAYS_PER_QUERY * n
    pairs = inner + final_cv + outer
    attach = ATTACH_CALLS_PER_QUERY * n
    ph = (pairs * ms_metric + attach * ms_attach + 2 * n * ms_load) / 3.6e6
    return {"inner_cv_scored_pairs_upper": inner, "final_cv_scored_pairs_upper": final_cv,
            "outer_scored_pairs_upper": outer, "all_scored_pairs_upper": pairs, "attachment_calls_upper": attach,
            "ms_per_scored_pair_referee": ms_metric, "ms_per_attachment_call": ms_attach, "ms_per_arbor_load": ms_load,
            "processor_hours": ph, "wall_hours": ph / workers, "workers": workers,
            "not_included": "model fitting (about 160 ridge/HGB fits), bootstrap; not measured here"}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("what", choices=["cache-meaning", "folds", "known-inputs", "mdd", "all"])
    args = p.parse_args(argv)
    bad = [v for v in THREADS if os.environ.get(v) != "1"]
    if bad:
        raise SystemExit(f"set {', '.join(bad)}=1 first")
    t0 = time.time()
    result = {"schema": "next_donor_selector_feasibility_v2", "training_only": True, "held_out_files_opened": 0, "soma_table_heldout_rows_parsed": 0,
              "source_sha256": {str(Path(__file__).relative_to(ROOT)): sha_file(__file__),
                                "analysis/brm_template.py": sha_file(ROOT / "analysis/brm_template.py"),
                                "analysis/d110_geometry.py": sha_file(ROOT / "analysis/d110_geometry.py"),
                                "analysis/brm_geometry_audit.py": sha_file(ROOT / "analysis/brm_geometry_audit.py"),
                                "reports/axon_dataset_manifest.tsv": sha_file(MANIFEST),
                                "data/interim/axon_engine/real/ingestion_log.json": sha_file(REAL / "ingestion_log.json"),
                                "data/interim/axon_dataset/local_manifest.tsv": sha_file(LOCAL_MANIFEST)}}
    if args.what in ("cache-meaning", "all"):
        result["cache_meaning"] = cache_meaning()
    if args.what in ("folds", "all"):
        meta, _ = training_metadata()
        result["folds"] = folds_report(load_real_cache(meta))
    if args.what in ("known-inputs", "all"):
        result["known_inputs"] = known_inputs_report()
    if args.what in ("mdd", "all"):
        result["mdd"] = mdd_report()
    if args.what == "all":
        ms = result["cache_meaning"]["real_library_sample"]["ms_per_pair_d110_centred_mean"]
        result["timing"] = timing_projection(ms, {m: int(v) for m, v in result["folds"]["union_of_candidate_label_pairs"].items()})
        result["timing"]["full_directed_recompute_all_cross_animal_pairs"] = {
            "pairs": result["folds"]["all_directed_cross_animal_pairs"],
            "processor_hours": result["folds"]["all_directed_cross_animal_pairs"] * ms / 3.6e6}
        rs = result["cache_meaning"]["real_library_sample"]
        result["timing"]["scoring"] = scoring_projection(result["folds"], rs["ms_per_pair_referee_label_mean"],
                                                         rs["attachment_ladder"]["ms_max_when_attached"] or rs["attachment_ladder"]["ms_max"],
                                                         rs["ms_per_arbor_load_mean"])
    result["seconds"] = time.time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"feasibility_{args.what}.json"
    path.write_text(json.dumps(result, sort_keys=True, indent=1, allow_nan=False) + "\n")
    try:
        shown = str(path.relative_to(ROOT))
    except ValueError:                      # OUT redirected outside the repository
        shown = str(path)
    print(json.dumps({"output": shown, "sha256": sha_file(path), "seconds": result["seconds"]}))
    return result


if __name__ == "__main__":
    main()
