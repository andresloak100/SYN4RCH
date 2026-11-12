"""Next donor selector, version 2: a learned choice among the 16 local donors, TRAIN-only development experiment.

Protocol: reports/next_donor_selector_v2_protocol.md. Output: data/interim/next_donor_selector_v2/.

Question: from known query inputs only, can a learned donor choice among the same 16 nearest donors that E1 and E2
use beat BOTH exact E1 and exact E2 (and every other control) on BOTH animal-mean distance and animal-mean F1 at 200 um,
with the emitted geometry made by the SAME frozen attachment ladder (500 um, then 1,000 um, else translation)?

TRAINING DATA ONLY. The 25 validation and 40 test brain samples are used (their traces were opened by the completed
E4 evaluation). This module never reads any validation or test row, trace, prediction or truth; it reads:
  * data/interim/axon_engine/real/ingestion_log.json (part == "train", status == "included" rows only)
  * reports/axon_dataset_manifest.tsv and data/interim/axon_dataset/local_manifest.tsv (train rows only, known fields)
  * the canonical training arbors data/interim/axon_engine/real/library/<id>.npz (file hashes checked)
  * the engine's complete all-pairs TRAIN cache D (soma-centred distance; used only for E1/E2 local medoids, the
    copy medoids and donor centrality, always through the fit-only door). The cached F is never read.
  * data/interim/axon_engine/real/develop/ (read-only), only for the parity check of the E1/E2/copy reimplementation.
Every frozen engine, referee, result and raw file is read-only. Nothing here imports the engine.

Subcommands (all output under data/interim/next_donor_selector_v2/):
  parity   reimplementation parity against the engine's committed development forecasts (no new truth needed)
  plan     folds, fit calls, candidate sets and the label-pair union (no truth opened)
  run      the whole experiment, continuously and resumably: labels (fit pairs) -> nested CV -> outer forecasts sealed
           -> outer scoring -> tables. Never repeated after RUN_COMPLETE.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/interim/next_donor_selector_v2"
REAL = ROOT / "data/interim/axon_engine/real"
LIBRARY = REAL / "library"
INGESTION_LOG = REAL / "ingestion_log.json"
PAIRCACHE = REAL / "paircache/85f0f247bacfd142"
DEVELOP = REAL / "develop"
PLAN_MANIFEST = ROOT / "reports/axon_dataset_manifest.tsv"
LOCAL_MANIFEST = ROOT / "data/interim/axon_dataset/local_manifest.tsv"
PROTOCOL = ROOT / "reports/next_donor_selector_v2_protocol.md"

# pins (read-only inputs)
PLAN_MANIFEST_SHA = "649f829a592cd3d13bcda5cffdbd049e759568d38bacb144f6a61220610a1ca4"
INGESTION_LOG_SHA = None          # recorded at run time, not pinned here (the engine re-wrote it once with the same rules)
CACHE_D_SHA = "fd197604091af2a712d9a9a31de1ae81484c7f860c48612ae7d1766b8363f165"
CACHE_IDS_SHA = "a9ea63cd0f407daa3931b21dd1f7b5c04da96d824bd08d79b04d469954d629c8"
TEMPLATE_SHA = "0a3da130c0e0049e8c9cc5df94523442f542f44a6815e82d3ad947a960a8d49a"     # analysis/brm_template.py (frozen)
GEOMETRY_D110_SHA = "3b011643b1b224256fdce1ad8a31014e01e644c6a794bd1d49af5e6543090a5e"  # analysis/d110_geometry.py (frozen)

# the engine's E1/E2 constants (analysis/axon_engine.py CONFIG["engine"]); checked against its source text by `parity`
K = 16
KERNEL_H0_UM = 100.0
KERNEL_RANK = 8
LADDER_UM = (500.0, 1000.0)
ENGINE_CONSTANT_TEXT = ('"k_neighbours": 16, "kernel_h0_um": 100.0, "kernel_rank": 8',
                        '"attachment_ladder_um": [500.0, 1000.0]')

SALT = "next_donor_selector_v2|20261112"
SEED = 20261112
BOOTSTRAP_SEED, BOOTSTRAP_DRAWS = 20261110, 10000
OUTER_FOLDS, INNER_FOLDS = 5, 4
LEARNERS = ("ridge_alpha100", "ridge_alpha1", "hgb")          # tie order: fewer effective parameters first
FALLBACKS = ("E1_local_uniform", "E2_local_kernel")
HGB = dict(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=50)
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
QUERY_KNOWN_FIELDS = ("soma", "peptide", "nucleus", "hemisphere")      # the ONLY query inputs to any prediction
EXCLUDED_QUERY_FIELDS = ("reconstruction_type",)                        # describes the trace, not the cell: never used

SELECTOR = "selector_v2"
CONTROLS = ("E1_local_uniform", "E2_local_kernel", "adaptive_fallback", "nearest", "medoid", "peptide_nearest", "peptide_medoid",
            "local_nearest", "local_label_nearest", "switch_nearest", "switch_random", "permuted_labels")
COPIES = ("nearest", "medoid", "peptide_nearest", "peptide_medoid")       # full fit pool, translated (referee rules)
ATTACHED_METHODS = (SELECTOR, "E1_local_uniform", "E2_local_kernel", "adaptive_fallback", "local_nearest", "local_label_nearest",
                    "switch_nearest", "switch_random", "permuted_labels")
ABLATION = (SELECTOR, "E1_local_uniform", "E2_local_kernel")            # also scored translated (donor-choice only)


# ------------------------------------------------------------------------------------------------------------------
# small helpers
# ------------------------------------------------------------------------------------------------------------------

def utc():
    return datetime.now(timezone.utc).isoformat()


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def arbor_sha(arbor):
    """The engine's content hash of an arbor (analysis/axon_engine.py arbor_sha), reimplemented."""
    h = hashlib.sha256()
    for array, kind in ((arbor.xyz, "<f8"), (arbor.parent, "<i8"), (arbor.types, "<i8")):
        data = np.ascontiguousarray(np.asarray(array).astype(kind, copy=False))
        h.update(str(data.shape).encode())
        h.update(data.tobytes())
    return h.hexdigest()


def write_json(path, value, *, exclusive=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x" if exclusive else "w") as f:
        f.write(json.dumps(value, sort_keys=True, indent=1, allow_nan=False) + "\n")
    return sha_file(path)


def read_json(path):
    return json.loads(Path(path).read_text())


def stable_seed(*parts):
    return int(hashlib.sha256("|".join(str(p) for p in (SEED,) + parts).encode()).hexdigest()[:8], 16)


class LeakError(RuntimeError):
    """A read outside the fit animals of the current fit call, or a held-out row."""


# ------------------------------------------------------------------------------------------------------------------
# frozen geometry (read-only import of the frozen template; it pins its own geometry dependency)
# ------------------------------------------------------------------------------------------------------------------

_T = None


def template():
    global _T
    if _T is None:
        path = ROOT / "analysis/brm_template.py"
        if sha_file(path) != TEMPLATE_SHA or sha_file(ROOT / "analysis/d110_geometry.py") != GEOMETRY_D110_SHA:
            raise RuntimeError("the frozen template or d110 geometry changed")
        spec = importlib.util.spec_from_file_location("next_selector_v2_frozen_template", path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        _T = module
    return _T


def attach(arbor, soma):
    """The engine's ladder (Engine.attach): frozen compact_attachment at 500 um, then 1,000 um, else translation of the
    same arbor. Distal coordinates beyond the radius are preserved exactly by the frozen operator."""
    T = template()
    for radius in LADDER_UM:
        try:
            out, _ = T.compact_attachment(arbor, soma, radius)
        except T.UnsupportedAttachment:
            continue
        return out, f"attach_{int(radius)}"
    return arbor.shifted(soma), "translation"


def emit(arbor, soma, mode):
    if mode == "attached":
        return attach(arbor, soma)
    if mode == "translated":
        return arbor.shifted(soma), "translation"
    raise ValueError(mode)


def metric(prediction, truth):
    m = template().G.geometry_metrics(prediction, truth)
    return float(m["distance_um"]), float(m["f1_200"])


# ------------------------------------------------------------------------------------------------------------------
# the training library: known inputs (query side) and TRAIN arbors (donor side)
# ------------------------------------------------------------------------------------------------------------------

def _tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


class Library:
    """The 4,368 included TRAINING neurons. Known query inputs: table soma (mirrored-left frame, as the engine),
    source label (peptide), soma label (nucleus), hemisphere. Nothing of a validation or test neuron is parsed."""

    def __init__(self):
        if sha_file(PLAN_MANIFEST) != PLAN_MANIFEST_SHA:
            raise RuntimeError("plan manifest changed")
        log = read_json(INGESTION_LOG)["neurons"]
        ids = sorted(n for n, e in log.items() if e["part"] == "train" and e["status"] == "included")
        keep = set(ids)
        plan = {r["neuron"]: r for r in _tsv(PLAN_MANIFEST) if r["neuron"] in keep}
        local = {r["neuron"]: r for r in _tsv(LOCAL_MANIFEST) if r["neuron"] in keep}
        if set(plan) != keep or set(local) != keep:
            raise RuntimeError("known-input rows missing for training neurons")
        self.ids = ids
        self.index = {n: i for i, n in enumerate(ids)}
        self.n = len(ids)
        self.animal = np.array([log[n]["animal"] for n in ids])
        for n in ids:
            if plan[n]["part"] != "train" or plan[n]["animal"] != log[n]["animal"] or local[n]["animal"] != log[n]["animal"]:
                raise RuntimeError(f"{n}: identity/part disagree between sources")
        self.soma = np.array([log[n]["table_soma_um"] for n in ids], dtype=float)
        self.peptide = np.array([plan[n]["source_label"] or "" for n in ids])
        self.nucleus = np.array([plan[n]["soma_label"] or "" for n in ids])
        self.hemisphere = np.array([local[n]["hemisphere"] or "" for n in ids])
        self.file_sha = {n: log[n]["canonical_file_sha256"] for n in ids}
        self.trace_soma = np.array([log[n]["trace_soma_um"] for n in ids], dtype=float)  # DONOR side only (attach shift)
        self.ingestion_log_sha256 = sha_file(INGESTION_LOG)
        self._summary = None

    def path(self, n):
        return LIBRARY / f"{n}.npz"

    def load(self, n, *, verify=True):
        if n not in self.index:
            raise LeakError(f"{n} is not an included training neuron")
        p = self.path(n)
        if verify and sha_file(p) != self.file_sha[n]:
            raise RuntimeError(f"library arbor {n} changed on disk")
        return template().G.load_arbor(p)

    def known(self, i):
        """The query's known inputs, nothing else."""
        return {"soma": self.soma[i].copy(), "peptide": str(self.peptide[i]), "nucleus": str(self.nucleus[i]),
                "hemisphere": str(self.hemisphere[i])}

    def donor_summaries(self):
        """TRAIN morphology summaries of every library arbor (used only on the DONOR side of a pair)."""
        if self._summary is None:
            path = OUT / "donor_summaries.npz"
            if path.exists():
                with np.load(path) as z:
                    if list(z["ids"]) != self.ids:
                        raise RuntimeError("donor summaries belong to another library")
                    self._summary = {k: z[k] for k in z.files if k != "ids"}
            else:
                rows = [summarise(self.load(n)) for n in self.ids]
                self._summary = {k: np.array([r[k] for r in rows], dtype=float) for k in rows[0]}
                OUT.mkdir(parents=True, exist_ok=True)
                np.savez(path, ids=np.array(self.ids), **self._summary)
        return self._summary


SUMMARY_KEYS = ("log_cable", "log_branches", "log_terminals", "rms_ap", "rms_dv", "rms_ml", "centroid_ap",
                "centroid_dv", "centroid_ml", "max_radius", "log_nodes")


def summarise(arbor):
    xyz = arbor.xyz - arbor.xyz[0]
    w = arbor.weights()
    w = w / w.sum()
    centroid = (xyz * w[:, None]).sum(axis=0)
    rms = np.sqrt((xyz ** 2 * w[:, None]).sum(axis=0))
    terminals = int(sum(1 for c in arbor.children()[1:] if not c))
    return {"log_cable": math.log1p(arbor.cable()), "log_branches": math.log1p(arbor.branches()),
            "log_terminals": math.log1p(terminals), "rms_ap": rms[0], "rms_dv": rms[1], "rms_ml": rms[2],
            "centroid_ap": centroid[0], "centroid_dv": centroid[1], "centroid_ml": centroid[2],
            "max_radius": float(np.linalg.norm(xyz, axis=1).max()), "log_nodes": math.log(len(xyz))}


class Cache:
    """The engine's complete TRAIN all-pairs cache: D only (soma-centred distance), verified, in memory."""

    def __init__(self, library):
        complete = read_json(PAIRCACHE / "COMPLETE.json")
        index = read_json(PAIRCACHE / "index.json")
        if complete["D_sha256"] != CACHE_D_SHA or complete["ids_sha256"] != CACHE_IDS_SHA:
            raise RuntimeError("cache pins differ")
        if index["ids"] != library.ids:
            raise RuntimeError("cache ID order differs from the library")
        if sha_file(PAIRCACHE / "D.npy") != CACHE_D_SHA:
            raise RuntimeError("cache D changed")
        self.D = np.load(PAIRCACHE / "D.npy")
        if self.D.shape != (library.n, library.n) or not np.isfinite(self.D).all():
            raise RuntimeError("cache D shape/finite")


class FitDoor:
    """Every read of a fit call goes through here: labels, cache D, centrality. Both endpoints must be neurons of the
    fit animals of this call; the query's animal never equals the donor's. Counts served and refused reads."""

    def __init__(self, name, library, fit_animals, *, labels=None, cache=None):
        self.name, self.lib = name, library
        self.fit_animals = frozenset(fit_animals)
        self.mask = np.isin(library.animal, list(self.fit_animals))
        self._labels, self._cache = labels, cache
        self.served = Counter()
        self.refused = Counter()

    def _check(self, i, what):
        if not self.mask[i]:
            self.refused[what] += 1
            raise LeakError(f"{self.name}: {what} read of {self.lib.ids[i]} outside the fit animals")

    def label(self, q, d):
        self._check(q, "label_query")
        self._check(d, "label_donor")
        if self.lib.animal[q] == self.lib.animal[d]:
            self.refused["label_same_animal"] += 1
            raise LeakError(f"{self.name}: same-animal pair")
        self.served["label"] += 1
        return self._labels.get(q, d)

    def D(self, rows, cols):
        rows, cols = np.asarray(rows), np.asarray(cols)
        for i in np.unique(np.concatenate([rows, cols])):
            self._check(int(i), "cache")
        self.served["cache_block"] += 1
        return self._cache.D[np.ix_(rows, cols)]


# ------------------------------------------------------------------------------------------------------------------
# folds, local candidate sets, E1/E2, copies
# ------------------------------------------------------------------------------------------------------------------

def folds(animals, salt, k):
    order = sorted(set(animals), key=lambda a: (hashlib.sha256(f"{SALT}|{salt}|{a}".encode()).hexdigest(), a))
    return [sorted(order[j::k]) for j in range(k)]


def fold_plan(library):
    animals = sorted(set(library.animal.tolist()))
    outer = folds(animals, "outer", OUTER_FOLDS)
    plan = []
    for k, ev in enumerate(outer):
        fit = sorted(set(animals) - set(ev))
        inner = folds(fit, f"inner|{k}", INNER_FOLDS)
        plan.append({"outer": k, "eval_animals": ev, "fit_animals": fit,
                     "inner": [{"inner": j, "eval_animals": ie, "fit_animals": sorted(set(fit) - set(ie))}
                               for j, ie in enumerate(inner)]})
    return plan


def local_candidates(library, q, pool_mask):
    """The engine's E1-E7 candidate set: the K nearest pool donors by table-soma distance (ties by position = ID),
    then sorted ascending by ID. pool_mask must already exclude the query's own animal."""
    if pool_mask[q] or (pool_mask & (library.animal == library.animal[q])).any():
        raise LeakError("the query's animal is in its own pool")
    pool = np.flatnonzero(pool_mask)
    d2 = ((library.soma[pool] - library.soma[q]) ** 2).sum(axis=1)
    k = min(K, len(pool))
    nearest = np.lexsort((np.arange(len(pool)), d2))[:k]
    index = np.sort(nearest)
    return pool[index], np.sqrt(d2[index])


def local_weights(library, cand, dist):
    animals = library.animal[cand]
    count = Counter(animals.tolist())
    prior = np.array([1.0 / count[a] for a in animals]) / len(count)
    h = max(KERNEL_H0_UM, float(np.sort(dist)[min(KERNEL_RANK, len(cand)) - 1]))
    kernel = np.exp(-dist ** 2 / (2.0 * h * h))
    raw = {"E1_local_uniform": prior, "E2_local_kernel": prior * kernel}
    return {name: w / w.sum() for name, w in raw.items()}, h


def local_choices(library, door, q, pool_mask):
    """Candidates, E1/E2 weights, risks and choices for one query from one pool (cache D through the door)."""
    T = template()
    cand, dist = local_candidates(library, q, pool_mask)
    Dsub = door.D(cand, cand)
    weights, h = local_weights(library, cand, dist)
    choice, risk = {}, {}
    for name, w in weights.items():
        i, r = T.weighted_medoid(Dsub, w)
        choice[name], risk[name] = int(i), r
    return {"q": q, "cand": cand, "dist": dist, "weights": weights, "risk": risk, "choice": choice, "h": h}


def animal_medoid(D, members, animals):
    count = Counter(animals)
    w = np.array([1.0 / count[a] for a in animals])
    w /= w.sum()
    return members[int(np.argmin(D @ w))]


class Copies:
    """The referee's four copying rules on a full pool (copying_baselines arithmetic; D from the verified cache,
    which the referee's own audit proved bit-identical)."""

    def __init__(self, library, door, pool_mask):
        self.lib, self.pool = library, np.flatnonzero(pool_mask)
        Dp = door.D(self.pool, self.pool)
        self.medoid = int(animal_medoid(Dp, self.pool, library.animal[self.pool].tolist()))
        self.label_medoid = {}
        for lab in sorted(set(library.peptide[self.pool].tolist())):
            if not lab:
                continue
            members = np.flatnonzero(library.peptide[self.pool] == lab)
            sub = Dp[np.ix_(members, members)]
            self.label_medoid[lab] = int(self.pool[animal_medoid(sub, members, library.animal[self.pool][members].tolist())])

    def choose(self, q):
        lib, pool = self.lib, self.pool
        if (lib.animal[pool] == lib.animal[q]).any():
            raise LeakError("query animal in the copy pool")
        d2 = ((lib.soma[pool] - lib.soma[q]) ** 2).sum(axis=1)
        nearest = int(pool[int(np.argmin(d2))])
        lab = lib.peptide[q]
        same = np.flatnonzero(lib.peptide[pool] == lab) if lab else np.array([], dtype=int)
        if len(same):
            pn = int(pool[same[int(np.argmin(d2[same]))]])
            pm = self.label_medoid[lab]
        else:
            pn, pm = nearest, self.medoid
        return {"nearest": nearest, "medoid": self.medoid, "peptide_nearest": pn, "peptide_medoid": pm,
                "peptide_fallback": not len(same)}


# ------------------------------------------------------------------------------------------------------------------
# features (known query inputs + TRAIN donor summaries + fit-only centrality + pair/local-structure terms)
# ------------------------------------------------------------------------------------------------------------------

class Featurizer:
    """Fitted on the fit animals of one fit call only: vocabularies, donor centrality. IDs/animals are never features."""

    def __init__(self, library, door):
        lib = self.lib = library
        fit = np.flatnonzero(door.mask)
        self.peptides = sorted(set(lib.peptide[fit].tolist()) - {""})
        self.nuclei = sorted(k for k, c in Counter(lib.nucleus[fit].tolist()).items() if k and c >= 20)
        self.hemis = sorted(set(lib.hemisphere[fit].tolist()) - {""})
        D = door.D(fit, fit)
        count = Counter(lib.animal[fit].tolist())
        w = np.array([1.0 / count[a] for a in lib.animal[fit]])
        w /= w.sum()
        self.centrality = np.full(lib.n, np.nan)
        self.centrality[fit] = D @ w                       # donors are fit neurons only; nan elsewhere (refused below)
        self.summary = lib.donor_summaries()
        self.names = None

    def rows(self, local, known):
        """Feature rows for the K candidates of one query. `known` = the query's known inputs (Library.known)."""
        lib = self.lib
        cand, dist = local["cand"], local["dist"]
        if not np.isfinite(self.centrality[cand]).all():
            raise LeakError("a candidate donor outside the fit animals")
        q_soma = np.asarray(known["soma"], float)
        diff = lib.soma[cand] - q_soma
        shift = np.linalg.norm(lib.trace_soma[cand] - q_soma, axis=1)         # the attachment's soma shift
        e1, e2 = local["choice"]["E1_local_uniform"], local["choice"]["E2_local_kernel"]
        r1, r2 = local["risk"]["E1_local_uniform"], local["risk"]["E2_local_kernel"]
        cols = {}
        for a, name in enumerate(("ap", "dv", "ml")):
            cols[f"q_soma_{name}"] = np.full(len(cand), q_soma[a])
            cols[f"diff_{name}"] = diff[:, a]
            cols[f"absdiff_{name}"] = np.abs(diff[:, a])
        cols["soma_distance"] = dist
        cols["rank"] = np.argsort(np.argsort(dist, kind="stable"), kind="stable").astype(float)
        cols["attach_shift"] = shift
        cols["attach_500_feasible"] = (shift <= LADDER_UM[0] / 3).astype(float)
        cols["attach_1000_feasible"] = (shift <= LADDER_UM[1] / 3).astype(float)
        for p in self.peptides:
            cols[f"q_pep_{p}"] = np.full(len(cand), float(known["peptide"] == p))
        for nuc in self.nuclei:
            cols[f"q_nuc_{nuc}"] = np.full(len(cand), float(known["nucleus"] == nuc))
        for hm in self.hemis:
            cols[f"q_hemi_{hm}"] = np.full(len(cand), float(known["hemisphere"] == hm))
        cols["same_peptide"] = (lib.peptide[cand] == known["peptide"]).astype(float)
        cols["same_nucleus"] = (lib.nucleus[cand] == known["nucleus"]).astype(float)
        cols["same_hemisphere"] = (lib.hemisphere[cand] == known["hemisphere"]).astype(float)
        cols["n_same_peptide_in_set"] = np.full(len(cand), cols["same_peptide"].sum())
        for key in SUMMARY_KEYS:
            cols[f"donor_{key}"] = self.summary[key][cand]
        cols["donor_centrality"] = self.centrality[cand]
        cols["w_E1"] = local["weights"]["E1_local_uniform"]
        cols["w_E2"] = local["weights"]["E2_local_kernel"]
        cols["risk_E1_rel"] = r1 - r1.min()
        cols["risk_E2_rel"] = r2 - r2.min()
        cols["is_E1_choice"] = (np.arange(len(cand)) == e1).astype(float)
        cols["is_E2_choice"] = (np.arange(len(cand)) == e2).astype(float)
        cols["kernel_h"] = np.full(len(cand), local["h"])
        names = list(cols)
        if self.names is None:
            self.names = names
        elif names != self.names:
            raise RuntimeError("feature columns changed within a fit call")
        X = np.column_stack([cols[n] for n in names]).astype(float)
        if not np.isfinite(X).all():
            raise RuntimeError("non-finite feature")
        return X


FORBIDDEN_FEATURE_WORDS = ("id", "animal", "neuron", "sample", "profile", "regional", "endpoint", "truth", "label_")


def check_feature_names(names):
    for n in names:
        low = n.lower()
        if any(w in low.split("_") or low.startswith(w) for w in FORBIDDEN_FEATURE_WORDS):
            raise LeakError(f"feature name {n} looks like an identity or hidden-outcome feature")


# ------------------------------------------------------------------------------------------------------------------
# learners (fixed small family) and the combined rule
# ------------------------------------------------------------------------------------------------------------------

def make_learner(name):
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.linear_model import Ridge
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    if name == "ridge_alpha1":
        return lambda: make_pipeline(StandardScaler(), Ridge(alpha=1.0))
    if name == "ridge_alpha100":
        return lambda: make_pipeline(StandardScaler(), Ridge(alpha=100.0))
    if name == "hgb":
        return lambda: HistGradientBoostingRegressor(random_state=SEED, **HGB)
    raise ValueError(name)


class TwoHeads:
    """Head 1: log1p(distance/200); head 2: F1 at 200 um. Pair weights equalise animals."""

    def __init__(self, name):
        self.name = name
        factory = make_learner(name)
        self.d, self.f = factory(), factory()

    def fit(self, X, yd, yf, w):
        kw = {"sample_weight": w}
        if self.name.startswith("ridge"):
            kw = {"ridge__sample_weight": w}
        self.d.fit(X, yd, **kw)
        self.f.fit(X, yf, **kw)
        return self

    def predict(self, X):
        return self.d.predict(X), self.f.predict(X)


def combined_choice(pred_d, pred_f, fallback_index):
    """Top donor by predicted joint loss (log-distance minus F1; ties by ID position); kept only if it is predicted
    better than the fallback donor on BOTH heads; otherwise the fallback donor (its exact arrays are emitted)."""
    if not (np.isfinite(pred_d).all() and np.isfinite(pred_f).all()):
        return int(fallback_index), "non_finite_prediction"
    score = pred_d - pred_f
    top = int(np.lexsort((np.arange(len(score)), score))[0])
    if top == fallback_index:
        return top, "top_is_fallback"
    if pred_d[top] < pred_d[fallback_index] and pred_f[top] > pred_f[fallback_index]:
        return top, "switched"
    return int(fallback_index), "not_better_on_both_heads"


def joint_loss(distance, f1):
    return math.log1p(distance / 200.0) + (1.0 - f1)


# ------------------------------------------------------------------------------------------------------------------
# label store (exact labels of the emitted geometry at the query's TABLE soma)
# ------------------------------------------------------------------------------------------------------------------

class LabelStore:
    """(query, donor, mode) -> (distance_um, f1_200, executed, arrays_sha256). Loaded from resumable JSONL shards."""

    def __init__(self, library):
        self.lib = library
        self.data = {}

    def add(self, q, d, mode, value):
        key = (q, d, mode)
        if key in self.data and self.data[key][:2] != value[:2]:
            raise RuntimeError(f"label {key} recomputed differently")
        self.data[key] = value

    def get(self, q, d, mode="attached"):
        return self.data[(q, d, mode)]

    def has(self, q, d, mode="attached"):
        return (q, d, mode) in self.data

    def load_jsonl(self, path):
        done = set()
        if not Path(path).exists():
            return done
        with open(path) as f:
            for line in f:
                if not line.strip():
                    continue
                r = json.loads(line)
                q = self.lib.index[r["query"]]
                for d, mode, dist, f1, executed, h in r["pairs"]:
                    self.add(q, self.lib.index[d], mode, (dist, f1, executed, h))
                done.add(r["job"])
        return done


_W_LIB = None


def _worker_init():
    global _W_LIB
    for v in THREAD_VARS:
        os.environ.setdefault(v, "1")
    _W_LIB = Library()
    template()


def label_job(job):
    """One query: load its measured TRAIN arbor (truth) once, emit each donor at the query's TABLE soma, score."""
    lib = _W_LIB
    qid, pairs = job["query"], job["pairs"]
    q = lib.index[qid]
    truth = lib.load(qid)
    soma = lib.soma[q]
    out = []
    cache = {}
    for did, mode in pairs:
        if lib.animal[lib.index[did]] == lib.animal[q]:
            raise LeakError("same-animal pair in a label job")
        if did not in cache:
            cache[did] = lib.load(did)
        arbor, executed = emit(cache[did], soma, mode)
        dist, f1 = metric(arbor, truth)
        out.append([did, mode, dist, f1, executed, arbor_sha(arbor)])
    return {"job": job["job"], "query": qid, "pairs": out}


def run_label_jobs(jobs, path, pool, log, what):
    """Resumable: finished jobs (by job key) are read back; the rest run on the one long-lived worker pool."""
    store_done = set()
    if Path(path).exists():
        with open(path) as f:
            for line in f:
                if line.strip():
                    store_done.add(json.loads(line)["job"])
    todo = [j for j in jobs if j["job"] not in store_done]
    total_pairs = sum(len(j["pairs"]) for j in jobs)
    done_pairs = total_pairs - sum(len(j["pairs"]) for j in todo)
    log(f"{what}: {len(jobs)} query jobs, {total_pairs} pairs; {len(jobs) - len(todo)} jobs already done; {len(todo)} to run")
    start, last, n_done = time.time(), time.time(), 0
    with open(path, "a") as f:
        for result in pool.imap_unordered(label_job, todo, chunksize=2):
            f.write(json.dumps(result, separators=(",", ":")) + "\n")
            f.flush()
            n_done += 1
            done_pairs += len(result["pairs"])
            now = time.time()
            if now - last > 60 or n_done == len(todo):
                rate = (done_pairs - (total_pairs - sum(len(j["pairs"]) for j in todo))) / max(now - start, 1e-9)
                left = total_pairs - done_pairs
                log(f"{what}: {done_pairs}/{total_pairs} pairs ({n_done}/{len(todo)} jobs this session), "
                    f"{rate:.1f} pairs/s, ETA {left / max(rate, 1e-9) / 60:.1f} min")
                last = now
    return {"jobs": len(jobs), "pairs": total_pairs, "resumed_jobs": len(jobs) - len(todo)}


# ------------------------------------------------------------------------------------------------------------------
# the fit calls: candidate sets for every (call, query)
# ------------------------------------------------------------------------------------------------------------------

def pool_mask_for(library, fit_animals, q):
    mask = np.isin(library.animal, list(fit_animals))
    mask &= library.animal != library.animal[q]
    return mask


def call_locals(library, cache, name, pool_animals, queries):
    """Local candidate sets of `queries` drawn from `pool_animals` (minus each query's own animal)."""
    door = FitDoor(name, library, pool_animals, cache=cache)
    base = np.isin(library.animal, list(pool_animals))
    out = {}
    for q in queries:
        out[q] = local_choices(library, door, q, base & (library.animal != library.animal[q]))
    return out, door


def neurons_of(library, animals):
    return [int(i) for i in np.flatnonzero(np.isin(library.animal, list(animals)))]


def build_calls(library, cache, plan):
    """All fit calls: outer fit (queries and pool = outer fit animals), inner fit and inner eval (pool = inner fit
    animals), and outer eval (pool = outer fit animals; scoring only)."""
    calls = {}
    for fold in plan:
        k = fold["outer"]
        fitq = neurons_of(library, fold["fit_animals"])
        calls[f"outer{k}_fit"], _ = call_locals(library, cache, f"outer{k}_fit", fold["fit_animals"], fitq)
        evq = neurons_of(library, fold["eval_animals"])
        calls[f"outer{k}_eval"], _ = call_locals(library, cache, f"outer{k}_eval", fold["fit_animals"], evq)
        for inner in fold["inner"]:
            j = inner["inner"]
            iq = neurons_of(library, inner["fit_animals"])
            calls[f"outer{k}_inner{j}_fit"], _ = call_locals(library, cache, f"outer{k}_inner{j}_fit", inner["fit_animals"], iq)
            eq = neurons_of(library, inner["eval_animals"])
            calls[f"outer{k}_inner{j}_eval"], _ = call_locals(library, cache, f"outer{k}_inner{j}_eval", inner["fit_animals"], eq)
    return calls


def fit_pairs(calls):
    """Union of (query, donor) pairs any FIT step reads (outer fit, inner fit, inner eval). Outer eval excluded."""
    pairs = defaultdict(set)
    for name, locs in calls.items():
        if name.endswith("_eval") and "_inner" not in name:
            continue
        for q, loc in locs.items():
            pairs[q].update(int(d) for d in loc["cand"])
    return pairs


# ------------------------------------------------------------------------------------------------------------------
# one fit call: training rows, fit, predict
# ------------------------------------------------------------------------------------------------------------------

def training_rows(library, door, feat, locs, *, permute_seed=None):
    X, yd, yf, w, owner = [], [], [], [], []
    per_animal = Counter(library.animal[q] for q in locs)
    for q in sorted(locs):
        loc = locs[q]
        X.append(feat.rows(loc, library.known(q)))
        for d in loc["cand"]:
            dist, f1 = door.label(q, int(d))[:2]
            yd.append(math.log1p(dist / 200.0))
            yf.append(f1)
            w.append(1.0 / per_animal[library.animal[q]])
            owner.append(q)
    X = np.vstack(X)
    yd, yf, w = np.array(yd), np.array(yf), np.array(w)
    w = w / w.mean()
    if permute_seed is not None:
        perm = np.random.default_rng(permute_seed).permutation(len(yd))
        yd, yf = yd[perm], yf[perm]
    return X, yd, yf, w


def predict_calls(library, model, feat, locs, fallback):
    out = {}
    for q in sorted(locs):
        loc = locs[q]
        pd_, pf_ = model.predict(feat.rows(loc, library.known(q)))
        fb = loc["choice"][fallback]
        i, why = combined_choice(pd_, pf_, fb)
        out[q] = {"index": i, "why": why, "fallback_index": fb, "pred_d": pd_.tolist(), "pred_f": pf_.tolist()}
    return out


def score_rows(library, door, locs, chosen_index):
    """Per-neuron (animal, distance, f1) for chosen candidate indices, labels through the door."""
    rows = []
    for q in sorted(locs):
        d = int(locs[q]["cand"][chosen_index[q]])
        dist, f1 = door.label(q, d)[:2]
        rows.append((library.animal[q], dist, f1))
    return rows


def animal_mean(rows):
    by = defaultdict(list)
    for a, dist, f1 in rows:
        by[a].append((dist, f1))
    per = {a: (float(np.mean([v[0] for v in vs])), float(np.mean([v[1] for v in vs]))) for a, vs in by.items()}
    return per, float(np.mean([v[0] for v in per.values()])), float(np.mean([v[1] for v in per.values()]))


def objective(rows):
    by = defaultdict(list)
    for a, dist, f1 in rows:
        by[a].append(joint_loss(dist, f1))
    return float(np.mean([np.mean(v) for v in by.values()]))


# ------------------------------------------------------------------------------------------------------------------
# bootstrap and verdict
# ------------------------------------------------------------------------------------------------------------------

def paired_bootstrap(per_a, per_b, animals):
    """Paired animal bootstrap 95% interval of (a - b) on equal-animal means, both measures."""
    a = np.array([per_a[x] for x in animals])
    b = np.array([per_b[x] for x in animals])
    diff = a - b
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    idx = rng.integers(0, len(animals), size=(BOOTSTRAP_DRAWS, len(animals)))
    boot = diff[idx].mean(axis=1)
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    return {"distance_delta_um": float(diff[:, 0].mean()), "f1_delta": float(diff[:, 1].mean()),
            "distance_bootstrap95_um": [float(lo[0]), float(hi[0])], "f1_bootstrap95": [float(lo[1]), float(hi[1])],
            "animals_both_better": int(((diff[:, 0] < 0) & (diff[:, 1] > 0)).sum())}


def verdict(summary, versus, switch_rate):
    sel = summary[SELECTOR]
    beats = {c: sel["distance_um"] < summary[c]["distance_um"] and sel["f1_200"] > summary[c]["f1_200"] for c in CONTROLS}
    # the adaptive (inner-selected, unchanged) fallback is a required control: a pass by adaptive fallback alone,
    # with harmful learned switches, must fail
    main = beats["E1_local_uniform"] and beats["E2_local_kernel"] and beats["adaptive_fallback"]
    gate = all(beats.values()) and switch_rate > 0
    resolved = gate and all(v["distance_bootstrap95_um"][1] < 0 and v["f1_bootstrap95"][0] > 0 for v in versus.values())
    if not gate:
        sentence = "development gate not met"
    elif resolved:
        sentence = "development gate met and resolved against every control"
    else:
        sentence = "development gate met, not resolved"
    return {"beats": beats, "beats_E1_E2_and_adaptive_fallback_on_both": main, "gate_met": gate, "resolved": resolved,
            "switch_rate": switch_rate, "sentence": sentence,
            "scope": ("DEVELOPMENT result on the 101 TRAIN brain samples, which have already informed the design; not an "
                      "unseen test. A development pass would need fresh independent data before any claim.")}


# ------------------------------------------------------------------------------------------------------------------
# parity against the engine's committed development forecasts
# ------------------------------------------------------------------------------------------------------------------

def parity(library, cache, *, n_arrays=None, log=print):
    """E1/E2 choices and arrays, copy donors, query somas and metrics against the engine's develop stage (pool: all
    other training animals, as the engine's development scope)."""
    text = (ROOT / "analysis/axon_engine.py").read_text()
    constants_ok = all(s in text for s in ENGINE_CONSTANT_TEXT)
    manifest = read_json(DEVELOP / "prediction_manifest.json")
    metrics = {r["neuron"]: r["metrics"] for r in read_json(DEVELOP / "metrics.json")["queries"]}
    rec = defaultdict(dict)
    for p in manifest["predictions"]:
        rec[p["query_id"]][p["model"]] = p
    queries = sorted(rec)
    all_animals = sorted(set(library.animal.tolist()))
    door_all = FitDoor("parity", library, all_animals, cache=cache)
    copies_by_animal = {}
    out = {"queries": len(queries), "engine_constants_in_source": constants_ok, "soma_equal": 0,
           "E1_donor_equal": 0, "E2_donor_equal": 0, "copy_donor_equal": 0, "copy_checked": 0,
           "arrays_checked": 0, "arrays_equal": 0, "metric_checked": 0, "metric_bit_equal": 0, "mismatches": []}
    T0 = time.time()
    arr_sample = queries if n_arrays is None else queries[:: max(1, len(queries) // n_arrays)][:n_arrays]
    for qid in queries:
        q = library.index[qid]
        r = rec[qid]
        if np.array_equal(np.array(r["nearest"]["known_soma_um"], float), library.soma[q]):
            out["soma_equal"] += 1
        mask = np.ones(library.n, bool) & (library.animal != library.animal[q])
        loc = local_choices(library, door_all, q, mask)
        for name in FALLBACKS:
            mine = library.ids[int(loc["cand"][loc["choice"][name]])]
            if mine == r[name]["selected_donor"]:
                out[name.split("_")[0] + "_donor_equal"] += 1
            else:
                out["mismatches"].append([qid, name, mine, r[name]["selected_donor"]])
        a = library.animal[q]
        if a not in copies_by_animal:
            copies_by_animal[a] = Copies(library, door_all, mask)
        ch = copies_by_animal[a].choose(q)
        for m in COPIES:
            out["copy_checked"] += 1
            if library.ids[ch[m]] == r[m]["selected_donor"]:
                out["copy_donor_equal"] += 1
            else:
                out["mismatches"].append([qid, m, library.ids[ch[m]], r[m]["selected_donor"]])
        if qid in arr_sample:
            truth = library.load(qid)
            for name in FALLBACKS:
                arbor, _ = attach(library.load(library.ids[int(loc["cand"][loc["choice"][name]])]), library.soma[q])
                out["arrays_checked"] += 1
                out["arrays_equal"] += int(arbor_sha(arbor) == r[name]["arrays_sha256"])
                dist, f1 = metric(arbor, truth)
                out["metric_checked"] += 1
                out["metric_bit_equal"] += int(dist == metrics[qid][name]["distance_um"] and f1 == metrics[qid][name]["f1_200"])
            for m in ("nearest", "medoid"):
                arbor = library.load(library.ids[ch[m]]).shifted(library.soma[q])
                out["arrays_checked"] += 1
                out["arrays_equal"] += int(arbor_sha(arbor) == r[m]["arrays_sha256"])
    out["seconds"] = round(time.time() - T0, 1)
    out["all_equal"] = (constants_ok and out["soma_equal"] == len(queries) and out["E1_donor_equal"] == len(queries)
                        and out["E2_donor_equal"] == len(queries) and out["copy_donor_equal"] == out["copy_checked"]
                        and out["arrays_equal"] == out["arrays_checked"] and out["metric_bit_equal"] == out["metric_checked"])
    out["mismatches"] = out["mismatches"][:50]
    log(f"parity: {json.dumps({k: v for k, v in out.items() if k != 'mismatches'})}")
    return out


# ------------------------------------------------------------------------------------------------------------------
# the experiment
# ------------------------------------------------------------------------------------------------------------------

class Log:
    def __init__(self, path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def __call__(self, message):
        line = f"[{utc()}] {message}"
        print(line, flush=True)
        with self.path.open("a") as f:
            f.write(line + "\n")


def code_pins():
    files = ["analysis/next_donor_selector_v2.py", "reports/next_donor_selector_v2_protocol.md",
             "tests/test_next_donor_selector_v2.py", "analysis/brm_template.py", "analysis/d110_geometry.py",
             "analysis/brm_geometry.py"]
    return {f: sha_file(ROOT / f) for f in files if (ROOT / f).exists()}


def protocol_final():
    text = PROTOCOL.read_text()
    return "STATUS: FINAL" in text and "STATUS: DRAFT" not in text


def git_clean(paths):
    import subprocess
    status = subprocess.run(["git", "status", "--porcelain", "--", *paths], cwd=ROOT, capture_output=True, text=True).stdout
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    return status.strip() == "", head


def run(workers, *, require_committed=True):
    import multiprocessing as mp
    for v in THREAD_VARS:
        if os.environ.get(v) != "1":
            raise RuntimeError(f"set {v}=1")
    OUT.mkdir(parents=True, exist_ok=True)
    log = Log(OUT / "run.log")
    if (OUT / "RUN_COMPLETE.json").exists():
        raise RuntimeError("RUN_COMPLETE exists: this experiment is never repeated")
    if (OUT / "RUN_FAILED.json").exists():
        raise RuntimeError("RUN_FAILED exists: a new attempt needs a new, disclosed protocol version")
    pinned = ["analysis/next_donor_selector_v2.py", "reports/next_donor_selector_v2_protocol.md",
              "tests/test_next_donor_selector_v2.py"]
    clean, head = git_clean(pinned)
    if require_committed and (not clean or not protocol_final()):
        raise RuntimeError("protocol FINAL and code committed and clean are required before the run")
    pins = code_pins()
    started = OUT / "RUN_STARTED.json"
    if started.exists():
        prev = read_json(started)
        if prev["code_sha256"] != pins:
            raise RuntimeError("code or protocol changed since RUN_STARTED; resume refused")
        log(f"resuming run started {prev['utc']}")
    else:
        write_json(started, {"utc": utc(), "git_head": head, "code_sha256": pins, "workers": workers,
                             "salt": SALT, "seed": SEED}, exclusive=True)
        log(f"RUN_STARTED at git {head}, {workers} workers")
    try:
        _run(workers, log, mp)
    except Exception as error:
        if not (OUT / "RUN_COMPLETE.json").exists():
            log(f"run interrupted by {type(error).__name__}: {error} (resumable unless it is a leak or a sealed-stage error)")
            if isinstance(error, LeakError):
                write_json(OUT / "RUN_FAILED.json", {"utc": utc(), "reason": f"{type(error).__name__}: {error}"}, exclusive=True)
        raise


def _run(workers, log, mp):
    T0 = time.time()
    lib = Library()
    cache = Cache(lib)
    lib.donor_summaries()
    log(f"library: {lib.n} training neurons, {len(set(lib.animal.tolist()))} brain samples; cache D verified")

    # parity of the reimplementation (no new truth: develop metrics already exist)
    par_path = OUT / "parity.json"
    if not par_path.exists():
        par = parity(lib, cache, n_arrays=60, log=log)
        write_json(par_path, par)
    par = read_json(par_path)
    if not par["all_equal"]:
        raise RuntimeError("parity with the engine's development forecasts failed; nothing runs")

    plan = fold_plan(lib)
    write_json(OUT / "folds.json", plan)
    calls = build_calls(lib, cache, plan)
    pairs = fit_pairs(calls)
    n_pairs = sum(len(v) for v in pairs.values())
    log(f"plan: {len(calls)} calls; fit-label union {n_pairs} pairs over {len(pairs)} queries")
    for name, locs in calls.items():
        for q, loc in locs.items():
            if (lib.animal[loc["cand"]] == lib.animal[q]).any():
                raise LeakError(f"{name}: same-animal candidate")

    ctx = mp.get_context("spawn")
    with ctx.Pool(workers, initializer=_worker_init) as pool:
        # PHASE 1: exact fit labels (attached emitted geometry at the TABLE soma vs the measured TRAIN arbor)
        jobs = [{"job": f"fit|{lib.ids[q]}|{digest(sorted(lib.ids[d] for d in ds))[:16]}", "query": lib.ids[q],
                 "pairs": [[lib.ids[d], "attached"] for d in sorted(ds)]} for q, ds in sorted(pairs.items())]
        info1 = run_label_jobs(jobs, OUT / "labels_fit.jsonl", pool, log, "fit labels")
        labels = LabelStore(lib)
        labels.load_jsonl(OUT / "labels_fit.jsonl")
        missing = [(q, d) for q, ds in pairs.items() for d in ds if not labels.has(q, d)]
        if missing:
            raise RuntimeError(f"{len(missing)} fit labels missing")
        log(f"fit labels complete: {len(labels.data)}")

        # PHASE 2: nested CV, outer forecasts (no outer-eval truth exists in the store)
        forecasts_path = OUT / "forecasts.json"
        if not (OUT / "FORECASTS_SEALED.json").exists():
            forecasts, selection = nested_cv(lib, cache, plan, calls, labels, log)
            sealed = emit_forecast_arrays(lib, forecasts, pool, log)
            h = write_json(forecasts_path, {"forecasts": sealed, "selection": selection})
            write_json(OUT / "FORECASTS_SEALED.json", {"utc": utc(), "forecasts_sha256": h,
                                                       "queries": len(sealed), "fit_label_pairs": len(labels.data)},
                       exclusive=True)
            log(f"FORECASTS_SEALED: {len(sealed)} outer forecasts, sha {h}")
        seal = read_json(OUT / "FORECASTS_SEALED.json")
        if sha_file(forecasts_path) != seal["forecasts_sha256"]:
            raise RuntimeError("sealed forecasts changed")
        sealed = read_json(forecasts_path)

        # PHASE 3: outer scoring (truth of outer-eval queries opened only now)
        score_jobs = []
        for qid, f in sorted(sealed["forecasts"].items()):
            items = {(d, "attached") for d in f["candidates"]}
            items |= {(f["donors"][m], "translated") for m in COPIES}
            items |= {(f["donors"][m], "translated") for m in ABLATION}
            items |= {(f["donors"][m], "attached") for m in ATTACHED_METHODS}
            score_jobs.append({"job": f"score|{qid}", "query": qid, "pairs": sorted([list(x) for x in items])})
        run_label_jobs(score_jobs, OUT / "labels_score.jsonl", pool, log, "outer scoring")
    scores = LabelStore(lib)
    scores.load_jsonl(OUT / "labels_score.jsonl")
    result = summarise_scores(lib, sealed, scores, log)
    result.update({"fit_label_pairs": info1["pairs"], "seconds_this_session": round(time.time() - T0, 1)})
    h = write_json(OUT / "results.json", result)
    write_tables(result)
    write_json(OUT / "RUN_COMPLETE.json", {"utc": utc(), "results_sha256": h,
                                           "forecasts_sha256": seal["forecasts_sha256"],
                                           "tables_sha256": sha_file(OUT / "tables.md"),
                                           "verdict": result["verdict"]["sentence"]}, exclusive=True)
    log(f"RUN_COMPLETE: {result['verdict']['sentence']}")


def nested_cv(lib, cache, plan, calls, labels, log):
    forecasts, selection = {}, []
    for fold in plan:
        k = fold["outer"]
        outer_fit = fold["fit_animals"]
        # inner CV: fallback (E1 or E2) and learner chosen on inner-eval queries of the OUTER FIT animals only
        scoring_door = FitDoor(f"outer{k}_innerscoring", lib, outer_fit, labels=labels, cache=cache)
        rows_fb = {fb: [] for fb in FALLBACKS}
        rows_cfg = {(lr, fb): [] for lr in LEARNERS for fb in FALLBACKS}
        doors = [scoring_door]
        for inner in fold["inner"]:
            j = inner["inner"]
            door = FitDoor(f"outer{k}_inner{j}", lib, inner["fit_animals"], labels=labels, cache=cache)
            doors.append(door)
            feat = Featurizer(lib, door)
            fit_locs, ev_locs = calls[f"outer{k}_inner{j}_fit"], calls[f"outer{k}_inner{j}_eval"]
            X, yd, yf, w = training_rows(lib, door, feat, fit_locs)
            check_feature_names(feat.names)
            for fb in FALLBACKS:
                rows_fb[fb] += score_rows(lib, scoring_door, ev_locs, {q: l["choice"][fb] for q, l in ev_locs.items()})
            for lr in LEARNERS:
                model = TwoHeads(lr).fit(X, yd, yf, w)
                for fb in FALLBACKS:
                    pred = predict_calls(lib, model, feat, ev_locs, fb)
                    rows_cfg[(lr, fb)] += score_rows(lib, scoring_door, ev_locs, {q: p["index"] for q, p in pred.items()})
        fb_obj = {fb: objective(rows_fb[fb]) for fb in FALLBACKS}
        fallback = min(FALLBACKS, key=lambda fb: (fb_obj[fb], FALLBACKS.index(fb)))
        cfg_obj = {lr: objective(rows_cfg[(lr, fallback)]) for lr in LEARNERS}
        learner = min(LEARNERS, key=lambda lr: (cfg_obj[lr], LEARNERS.index(lr)))
        # outer refit on all outer-fit animals; forecasts for outer-eval animals from known inputs only
        door = FitDoor(f"outer{k}_fit", lib, outer_fit, labels=labels, cache=cache)
        doors.append(door)
        feat = Featurizer(lib, door)
        X, yd, yf, w = training_rows(lib, door, feat, calls[f"outer{k}_fit"])
        check_feature_names(feat.names)
        model = TwoHeads(learner).fit(X, yd, yf, w)
        Xp, ydp, yfp, wp = training_rows(lib, door, feat, calls[f"outer{k}_fit"], permute_seed=stable_seed("permuted", k))
        permuted = TwoHeads(learner).fit(Xp, ydp, yfp, wp)
        ev_locs = calls[f"outer{k}_eval"]
        pred = predict_calls(lib, model, feat, ev_locs, fallback)
        pred_perm = predict_calls(lib, permuted, feat, ev_locs, fallback)
        # a probe read of an outer-eval label must be refused by the outer-fit door
        probe_q = next(iter(ev_locs))
        try:
            door.label(probe_q, int(ev_locs[probe_q]["cand"][0]))
            raise RuntimeError("the outer-fit door served an outer-eval label")
        except LeakError:
            pass
        copies = Copies(lib, door, np.isin(lib.animal, outer_fit))
        for q, loc in ev_locs.items():
            cand = loc["cand"]
            p = pred[q]
            fb_i = loc["choice"][fallback]
            switched = p["why"] == "switched"
            others = [i for i in range(len(cand)) if i != fb_i]
            near = int(np.argmin(loc["dist"]))
            same = [i for i in range(len(cand)) if lib.peptide[cand[i]] == lib.peptide[q] and lib.peptide[q]]
            label_near = min(same, key=lambda i: (loc["dist"][i], i)) if same else near
            sw_near = min(others, key=lambda i: (loc["dist"][i], i)) if switched else fb_i
            sw_rand = others[stable_seed("switch_random", lib.ids[q]) % len(others)] if switched else fb_i
            ch = copies.choose(q)
            idx = {SELECTOR: p["index"], "E1_local_uniform": loc["choice"]["E1_local_uniform"],
                   "E2_local_kernel": loc["choice"]["E2_local_kernel"], "adaptive_fallback": fb_i, "local_nearest": near,
                   "local_label_nearest": label_near, "switch_nearest": sw_near, "switch_random": sw_rand,
                   "permuted_labels": pred_perm[q]["index"]}
            donors = {m: lib.ids[int(cand[i])] for m, i in idx.items()}
            donors.update({m: lib.ids[ch[m]] for m in COPIES})
            forecasts[lib.ids[q]] = {"outer": k, "animal": lib.animal[q], "peptide": lib.peptide[q],
                                     "nucleus": lib.nucleus[q], "candidates": [lib.ids[int(d)] for d in cand],
                                     "fallback": fallback, "learner": learner, "why": p["why"], "switched": switched,
                                     "permuted_switched": pred_perm[q]["why"] == "switched",
                                     "donors": donors, "pred_d": p["pred_d"], "pred_f": p["pred_f"],
                                     "peptide_copy_fallback": ch["peptide_fallback"]}
        refused = sum(sum(d.refused.values()) for d in doors)
        sel = {"outer": k, "fallback_objective": fb_obj, "fallback": fallback, "learner_objective": cfg_obj,
               "learner": learner, "eval_queries": len(ev_locs), "outer_fit_rows": int(len(yd)),
               "switch_rate": float(np.mean([pred[q]["why"] == "switched" for q in ev_locs])),
               "door_reads": {d.name: dict(d.served) for d in doors}, "door_refusals_excluding_probe": refused - 1,
               "features": feat.names}
        if refused != 1:
            raise LeakError(f"outer {k}: {refused - 1} unexpected refused reads")
        selection.append(sel)
        log(f"outer {k}: fallback {fallback} {fb_obj}, learner {learner} {cfg_obj}, switch rate {sel['switch_rate']:.3f}")
    return forecasts, selection


def _emit_job(job):
    lib = _W_LIB
    q = lib.index[job["query"]]
    out = {}
    for m, (did, mode) in job["items"].items():
        arbor, executed = emit(lib.load(did), lib.soma[q], mode)
        out[m] = {"executed": executed, "arrays_sha256": arbor_sha(arbor)}
    return job["query"], out


def emit_forecast_arrays(lib, forecasts, pool, log):
    """Forecast geometry of every method (no truth read): the arrays' hashes are sealed with the donor choices."""
    jobs = []
    for qid, f in sorted(forecasts.items()):
        items = {m: (f["donors"][m], "attached") for m in ATTACHED_METHODS}
        items.update({m: (f["donors"][m], "translated") for m in COPIES})
        jobs.append({"query": qid, "items": items})
    done = 0
    for qid, out in pool.imap_unordered(_emit_job, jobs, chunksize=4):
        forecasts[qid]["arrays"] = out
        done += 1
        if done % 500 == 0:
            log(f"forecast geometry: {done}/{len(jobs)}")
    for qid, f in forecasts.items():
        fbm = f["fallback"]
        if f["donors"]["adaptive_fallback"] != f["donors"][fbm] or f["arrays"]["adaptive_fallback"] != f["arrays"][fbm]:
            raise RuntimeError(f"{qid}: the adaptive-fallback control is not the unchanged inner-selected fallback")
        if not f["switched"] and f["arrays"][SELECTOR] != f["arrays"]["adaptive_fallback"]:
            raise RuntimeError(f"{qid}: a fallback forecast does not carry the fallback's exact arrays")
    return forecasts


def summarise_scores(lib, sealed, scores, log):
    F = sealed["forecasts"]
    methods = (SELECTOR,) + CONTROLS + ("oracle_best_of_16",) + tuple(f"{m}__translated" for m in ABLATION)
    per_neuron = {m: [] for m in methods}
    for qid, f in sorted(F.items()):
        q = lib.index[qid]
        a = f["animal"]
        for m in ATTACHED_METHODS:
            d = lib.index[f["donors"][m]]
            dist, f1, executed, h = scores.get(q, d, "attached")
            if h != f["arrays"][m]["arrays_sha256"]:
                raise RuntimeError(f"{qid} {m}: scored geometry differs from the sealed forecast")
            per_neuron[m].append((a, dist, f1))
        for m in COPIES:
            d = lib.index[f["donors"][m]]
            dist, f1, executed, h = scores.get(q, d, "translated")
            if h != f["arrays"][m]["arrays_sha256"]:
                raise RuntimeError(f"{qid} {m}: scored geometry differs from the sealed forecast")
            per_neuron[m].append((a, dist, f1))
        for m in ABLATION:
            d = lib.index[f["donors"][m]]
            dist, f1 = scores.get(q, d, "translated")[:2]
            per_neuron[f"{m}__translated"].append((a, dist, f1))
        cand = [scores.get(q, lib.index[d], "attached")[:2] for d in f["candidates"]]
        best = min(range(len(cand)), key=lambda i: (cand[i][0], i))
        per_neuron["oracle_best_of_16"].append((a, *cand[best]))
    per_animal, summary = {}, {}
    for m, rows in per_neuron.items():
        per, dist, f1 = animal_mean(rows)
        per_animal[m] = per
        summary[m] = {"distance_um": dist, "f1_200": f1}
    animals = sorted(per_animal[SELECTOR])
    versus = {c: paired_bootstrap(per_animal[SELECTOR], per_animal[c], animals) for c in CONTROLS}
    switch_rate = float(np.mean([f["switched"] for f in F.values()]))
    v = verdict(summary, versus, switch_rate)
    # switched-subset comparison and ablation (fixed, descriptive)
    sw = [qid for qid, f in F.items() if f["switched"]]
    sub = {}
    if sw:
        for m in (SELECTOR, "adaptive_fallback", "E1_local_uniform", "E2_local_kernel", "switch_nearest", "switch_random"):
            rows = [(F[q]["animal"], *scores.get(lib.index[q], lib.index[F[q]["donors"][m]], "attached")[:2]) for q in sw]
            sub[m] = {"neurons": len(rows), "neuron_mean_distance_um": float(np.mean([r[1] for r in rows])),
                      "neuron_mean_f1": float(np.mean([r[2] for r in rows]))}
    ablation = {}
    for m in ABLATION:
        ablation[m] = {"attached": summary[m], "translated": summary[f"{m}__translated"],
                       "attachment_effect_distance_um": summary[m]["distance_um"] - summary[f"{m}__translated"]["distance_um"],
                       "attachment_effect_f1": summary[m]["f1_200"] - summary[f"{m}__translated"]["f1_200"]}
    for c in ("E1_local_uniform", "E2_local_kernel"):
        ablation[f"donor_choice_effect_vs_{c}_translated"] = paired_bootstrap(
            per_animal[f"{SELECTOR}__translated"], per_animal[f"{c}__translated"], animals)
    executed = Counter(F[q]["arrays"][SELECTOR]["executed"] for q in F)
    coverage = defaultdict(lambda: {"neurons": 0, "switched": 0, "animals": set()})
    for qid, f in F.items():
        c = coverage[f["peptide"] or "(none)"]
        c["neurons"] += 1
        c["switched"] += int(f["switched"])
        c["animals"].add(f["animal"])
    coverage = {k: {"neurons": v_["neurons"], "animals": len(v_["animals"]), "switched": v_["switched"],
                    "switch_rate": v_["switched"] / v_["neurons"]} for k, v_ in sorted(coverage.items())}
    label_perf = {}
    for lab in coverage:
        rows_by = {m: [r for r, (qid, f) in zip(per_neuron[m], sorted(F.items())) if (f["peptide"] or "(none)") == lab]
                   for m in (SELECTOR, "adaptive_fallback", "E1_local_uniform", "E2_local_kernel")}
        label_perf[lab] = {m: dict(zip(("distance_um", "f1_200"), animal_mean(rows)[1:])) for m, rows in rows_by.items()}
    return {"schema": "next_donor_selector_v2_results_v1", "label": "DEVELOPMENT (TRAIN nested CV), not an unseen test",
            "queries": len(F), "animals": len(animals), "summary": summary, "versus": versus, "verdict": v,
            "selection": sealed["selection"], "switch_rate": switch_rate,
            "why_counts": dict(Counter(f["why"] for f in F.values())),
            "permuted_switch_rate": float(np.mean([f["permuted_switched"] for f in F.values()])),
            "fallback_by_outer": {s["outer"]: s["fallback"] for s in sealed["selection"]},
            "learner_by_outer": {s["outer"]: s["learner"] for s in sealed["selection"]},
            "selector_emitted_geometry": dict(executed), "switched_subset_neuron_means": sub, "ablation": ablation,
            "per_label_coverage": coverage, "per_label_animal_means": label_perf,
            "per_animal": {m: per_animal[m] for m in methods}}


def write_tables(r):
    s, v = r["summary"], r["verdict"]
    lines = ["# Next donor selector v2: DEVELOPMENT result (TRAIN nested CV)", "",
             f"{r['queries']} TRAIN neurons, {r['animals']} brain samples; neuron means within sample, then equal samples.",
             "This is a development result: TRAIN has already informed the design. Not an unseen test.", "",
             f"**{v['sentence']}.** Beats E1, E2 and the adaptive fallback on both measures: {v['beats_E1_E2_and_adaptive_fallback_on_both']}. "
             f"Switch rate {r['switch_rate']:.3f} (permuted-label control {r['permuted_switch_rate']:.3f}).", "",
             "| method | distance um | F1 at 200 um |", "|---|---:|---:|"]
    for m in (SELECTOR,) + CONTROLS + ("oracle_best_of_16",):
        lines.append(f"| {m} | {s[m]['distance_um']:.1f} | {s[m]['f1_200']:.4f} |")
    lines += ["", "Selector minus control (paired brain-sample bootstrap 95%; negative distance and positive F1 favour the selector):", "",
              "| control | distance um | F1 | samples both better | beats on both |", "|---|---|---|---:|---|"]
    for c, x in r["versus"].items():
        lines.append(f"| {c} | {x['distance_delta_um']:+.1f} ({x['distance_bootstrap95_um'][0]:+.1f} to {x['distance_bootstrap95_um'][1]:+.1f}) | "
                     f"{x['f1_delta']:+.4f} ({x['f1_bootstrap95'][0]:+.4f} to {x['f1_bootstrap95'][1]:+.4f}) | {x['animals_both_better']} | {v['beats'][c]} |")
    lines += ["", "Ablation (fixed): attached vs translated, and donor choice with translation only:", "",
              "| method | attached dist / F1 | translated dist / F1 |", "|---|---|---|"]
    for m in ABLATION:
        a = r["ablation"][m]
        lines.append(f"| {m} | {a['attached']['distance_um']:.1f} / {a['attached']['f1_200']:.4f} | "
                     f"{a['translated']['distance_um']:.1f} / {a['translated']['f1_200']:.4f} |")
    for c in FALLBACKS:
        x = r["ablation"][f"donor_choice_effect_vs_{c}_translated"]
        lines.append(f"\nDonor choice only (both translated), selector minus {c}: {x['distance_delta_um']:+.1f} um "
                     f"({x['distance_bootstrap95_um'][0]:+.1f} to {x['distance_bootstrap95_um'][1]:+.1f}), F1 {x['f1_delta']:+.4f} "
                     f"({x['f1_bootstrap95'][0]:+.4f} to {x['f1_bootstrap95'][1]:+.4f}).")
    lines += ["", f"Inner-CV choices by outer fold: fallback {r['fallback_by_outer']}; learner {r['learner_by_outer']}.",
              f"Combined-rule branches: {r['why_counts']}. Selector emitted geometry: {r['selector_emitted_geometry']}.", "",
              "| source label | neurons | samples | switched | rate |", "|---|---:|---:|---:|---:|"]
    for k, c in r["per_label_coverage"].items():
        lines.append(f"| {k} | {c['neurons']} | {c['animals']} | {c['switched']} | {c['switch_rate']:.3f} |")
    (OUT / "tables.md").write_text("\n".join(lines) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("command", choices=("parity", "plan", "run"))
    p.add_argument("--workers", type=int, default=4)
    args = p.parse_args(argv)
    if args.workers < 1 or args.workers > 4:
        raise SystemExit("1 to 4 workers")
    if args.command == "parity":
        lib = Library()
        out = parity(lib, Cache(lib), n_arrays=60)
        write_json(OUT / "parity_check.json", out)
        return 0 if out["all_equal"] else 1
    if args.command == "plan":
        lib = Library()
        cache = Cache(lib)
        plan = fold_plan(lib)
        t = time.time()
        calls = build_calls(lib, cache, plan)
        pairs = fit_pairs(calls)
        ev = sum(len(v) for n, v in calls.items() if n.endswith("_eval") and "_inner" not in n)
        info = {"folds": [{"outer": f["outer"], "eval_animals": len(f["eval_animals"]), "fit_animals": len(f["fit_animals"]),
                           "inner_eval_animals": [len(i["eval_animals"]) for i in f["inner"]]} for f in plan],
                "calls": len(calls), "fit_label_pairs": sum(len(v) for v in pairs.values()),
                "fit_label_queries": len(pairs), "outer_eval_queries": ev, "seconds": round(time.time() - t, 1)}
        print(json.dumps(info, indent=1))
        write_json(OUT / "plan_check.json", info)
        return 0
    run(args.workers)
    return 0


if __name__ == "__main__":
    sys.exit(main())
