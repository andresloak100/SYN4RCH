"""D-109: predict the regional axon length of held-out hypothalamic neurons (Jiao et al. 2025).

Protocol: reports/d109_protocol.md.  Conventions and sources: reports/d109_jiao_conventions.md.
Input audit: reports/d109_jiao_audit.md.

OUTCOME RULE (absolute)
  The axon-length table (all_neurons_projection_length_fig1.csv) is read by exactly one function,
  `_read_outcome_table`. It is called from two places only:
    * `outcomes_permuted`, which permutes whole outcome rows across neurons with a fixed seed before
      returning anything (used by `check` and `dry`); every line printed in those modes carries
      the word PERMUTED or SYNTHETIC, or CHECK for lines that involve inputs only (split assertions);
    * `outcomes_real`, which exits unless the guard in `run` has passed.
  `splits` and `sources` never open the table. `hashes` (and the hash comparison in `run`) reads its
  bytes to compute the SHA-256 and parses nothing. The outcome-derived columns of the soma tables
  (`type`, `cluster`, `Class`) are never read: the loaders name the columns they take.
  The guard is in-process only: it protects the scripted path. No other code may import this module,
  call `_read_outcome_table`, set `_RUN_OPEN` or read the table before the run (protocol, section 14).

MODES   python3 -I analysis/d109_jiao.py MODE
  splits   write data/interim/d109_splits.json from IDs, soma positions, labels and tags only.
  check    split/leakage assertions, test-invariance of every predictor, and the planted controls
           on SYNTHETIC outcomes (planted tag effect, no tag effect, animal-only leak detector with
           a deliberately leaky positive control). Reduced permutation counts, stated in the output.
  dry      the whole evaluation path of `run` on PERMUTED outcomes. NOT A RESULT.
  hashes   print the SHA-256 of every input, the protocol and the split file (the PINNED block), the
           SHA-256 of this script with the PINNED block blanked, and the library versions.
  sources  verify the local copies of the source documentation against the SHA-256 recorded here;
           `sources --fetch [dir]` fetches them again (small documentation files only).
  run      the pre-registered evaluation, once. Refuses unless PINNED is filled and matches; the
           review file is pinned, carries the line `VERDICT: APPROVED TO RUN` and names the SHA-256 of
           this script (PINNED block blanked), of the protocol and of the split file; this file, the
           protocol, the conventions report and the review are committed and clean; the split file
           reproduces; `check` passed with this script; scikit-learn is the declared version. Writes
           data/interim/d109_run_opened.json before the outcome table is read and refuses if it
           exists, unless D109_RERUN=1 (a repeat writes to new file names and never overwrites).
"""
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import time
import warnings
import zlib

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np  # noqa: E402

# ------------------------------------------------------------------------------- constants
SEED = 20261110
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT_REL = "analysis/d109_jiao.py"
PROTOCOL_REL = "reports/d109_protocol.md"
REVIEW_REL = "reports/d109_review.md"          # written by the independent reviewers, not by this script
INPUTS = {
    "soma": "data/raw/jiao2025/All_soma_infor_HY_neurons.csv",
    "labels": "data/raw/jiao2025/geno_location_labels.csv",
    "outcome": "data/raw/jiao2025/all_neurons_projection_length_fig1.csv",
    "ccf": "data/raw/atlas/structure_graph_1.json",
}
INTERIM = os.path.join(ROOT, "data/interim")
SPLITS = os.path.join(INTERIM, "d109_splits.json")
RUN_MARKER = os.path.join(INTERIM, "d109_run_opened.json")
DOCS = os.path.join(INTERIM, "d109_docs")

# Filled after review from the printout of `hashes`; an empty dict makes `run` refuse.
# Required keys: soma, labels, outcome, ccf, protocol, splits, review.
PINNED = {
}

REGIONS = ("ARH PMv TU TMv MM AHN VMH DMH PVH MPO MPN PH LHA LPO ZI ME MRN PAG SCm SNr VTA DR MY IO NTS PARN IRN "
           "PB PRNr PCG PA BMA BLA MS BST SI MEA CEA LS ACB LH PVT MD PT AM SPFm HPF Isocortex OLF").split()
N_NEURONS, N_ANIMALS, N_TAGS, N_ALLZERO = 7180, 166, 16, 13
FIBRE = "fiber tracts"

K_OUTER, K_INNER = 5, 3
S2_MIN_NEURONS, S2_MIN_ANIMALS = 30, 5        # a nucleus is a hold-out fold only above both
NUCLEUS_MIN = 5                               # back-off threshold of the nucleus-mean baseline
GRID = {
    "global": [None], "nucleus": [None],
    "knn": [5, 10, 20, 40, 80], "knn_tag": [5, 10, 20, 40, 80],
    "ridge": [0.1, 1.0, 10.0, 100.0, 1000.0], "ridge_tag": [0.1, 1.0, 10.0, 100.0, 1000.0],
    "et": [3, 10, 30, 100], "et_tag": [3, 10, 30, 100],       # min_samples_leaf
    "mlp": [1e-4, 1e-2, 1.0], "mlp_tag": [1e-4, 1e-2, 1.0],   # L2 penalty alpha
}
MODELS = list(GRID)
# Under S3 the test tag is never in training, so a with-tag model is its position-only counterpart by the
# fallback; with-tag rows are not fitted or printed there.
MODELS_OF = {"S1": MODELS, "S2": MODELS, "S3": [m for m in MODELS if not m.endswith("_tag")]}
SKLEARN_VERSION = "1.8.0"                     # declared in the protocol; `run` refuses on any other
VERDICT_LINE = "VERDICT: APPROVED TO RUN"     # must be a line of the review file
# The 37 neurons (32 animals) whose SWC reconstructions were downloaded and viewed, and partly scored on 3D
# geometry, by the separate decisions D-110 / D-111 before this protocol was frozen (protocol, section 16).
# `<brain>_<n>` stems from data/raw/jiao2025/swc/MANIFEST.tsv. Used only for one declared sensitivity analysis.
SWC_SAMPLE = ("192710_002 195038_041 195039_009 201221_022 201230_079 201232_040 201232_056 201894_001 201922_073 "
              "202170_016 202457_066 202559_043 202617_025 202618_143 202620_016 202813_009 202813_035 202824_062 "
              "202825_065 202825_072 210007_142 210073_044 210104_007 210105_001 210796_002 210951_004 211266_057 "
              "211353_070 211461_026 211465_137 211595_037 211595_063 211889_096 211890_005 212202_011 212202_026 "
              "212248_071").split()
TAG_SIZE_GROUPS = ((2, 3), (4, 8), (9, 15), (16, 99))     # animals per tag, for the descriptive breakdown
ET_TREES = 100
MLP_HIDDEN, MLP_EPOCHS = (64, 64), 60
PRIMARY = ("et_tag", "et")                    # (with tag, without tag), split S1
MIN_EFFECT = 0.02                             # relative reduction of the per-animal error
N_FLIP = N_BOOT = 10000
N_SHUFFLE, N_SHUFFLE_UNSTRAT = 200, 100       # `run`; `check` and `dry` use the reduced counts below
DRY_SHUFFLE, DRY_SHUFFLE_UNSTRAT = 3, 2
CHECK_SHUFFLE = 39
LABEL = ""                                    # PERMUTED / SYNTHETIC / REAL prefix of every printed line
_RUN_OPEN = False


def say(*a):
    print((LABEL + " " if LABEL else "") + " ".join(str(x) for x in a), flush=True)


def path(rel):
    return os.path.join(ROOT, rel)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rng_for(*key):
    return np.random.default_rng([SEED, zlib.crc32("|".join(str(k) for k in key).encode())])


# ------------------------------------------------------------------------------- inputs (no outcomes)
class Inputs:
    pass


def load_ccf():
    par = {}

    def walk(n, p):
        par[n["acronym"]] = p
        for c in n["children"]:
            walk(c, n["acronym"])
    walk(json.load(open(path(INPUTS["ccf"])))["msg"][0], None)
    return par


def ancestors(a, par):
    """[a, parent, ..., root]"""
    out = []
    while a is not None:
        out.append(a)
        a = par[a]
    return out


def load_inputs():
    import pandas as pd
    s = pd.read_csv(path(INPUTS["soma"]), usecols=["file", "geno", "location", "AP", "DV", "ML"])   # not `type`
    g = pd.read_csv(path(INPUTS["labels"]))
    g.columns = ["id", "geno", "loc16"]
    D = Inputs()
    D.ids = s["file"].to_numpy(str)
    parts = [i.split("_", 3) for i in D.ids]
    assert all(len(p) == 4 for p in parts)
    D.animal = np.array([p[0] for p in parts])
    D.tagname = s["geno"].to_numpy(str)
    D.fine = np.array([p[3] for p in parts])              # the fine label: suffix of the published ID
    D.key50 = s["location"].to_numpy(str)                 # the authors' 50-label column
    D.loc16 = g["loc16"].to_numpy(str)
    D.X = s[["AP", "DV", "ML"]].to_numpy(float)           # as published; micrometres in the CCFv3 10 um frame and
    #                                                       mirrored to one side are INFERRED (conventions report, 2-3)
    D.par = load_ccf()
    # layout checks (inputs only)
    assert len(D.ids) == N_NEURONS == len(set(D.ids)) and (g["id"].to_numpy(str) == D.ids).all()
    assert (g["geno"].to_numpy(str) == D.tagname).all() and all(p[2] == t for p, t in zip(parts, D.tagname))
    D.animals = np.array(sorted(set(D.animal)))
    assert len(D.animals) == N_ANIMALS
    D.animal_tag = {}
    for a, t in zip(D.animal, D.tagname):
        assert D.animal_tag.setdefault(a, t) == t, "an animal with two tags"
    D.tags = np.array(sorted(set(D.tagname)))
    assert len(D.tags) == N_TAGS
    assert all(f in D.par for f in D.fine), "a fine label is not a CCF acronym"
    assert np.isfinite(D.X).all() and (D.X[:, 2] < 5700).all()
    # hold-out key: the fine label, pooled to MBO / ZI when it lies inside them; must equal the authors' column
    D.key = np.array([next((p for p in ("MBO", "ZI") if p in ancestors(f, D.par)), f) for f in D.fine])
    assert (D.key == D.key50).all(), "CCF pooling rule does not reproduce the authors' 50-label column"
    stems = np.array(["_".join(p[:2]) for p in parts])
    assert len(set(stems)) == N_NEURONS and set(SWC_SAMPLE) <= set(stems) and len(set(SWC_SAMPLE)) == 37
    D.swc_animals = sorted(set(x.split("_")[0] for x in SWC_SAMPLE))
    assert len(D.swc_animals) == 32
    return D


# ------------------------------------------------------------------------------- outcomes
def _read_outcome_table(ids):
    """The only reader of the axon-length table. Callers: outcomes_permuted, outcomes_real."""
    import pandas as pd
    t = pd.read_csv(path(INPUTS["outcome"]), index_col=0)
    assert list(t.columns) == REGIONS, "outcome columns differ from the 49 pre-registered regions"
    assert len(t) == N_NEURONS and set(t.index) == set(ids)
    L = t.loc[list(ids)].to_numpy(float)
    assert np.isfinite(L).all() and (L >= 0).all()
    return L


def target(L):
    """log10(1 + axon length in micrometres), 49 regions as published."""
    return np.log10(1.0 + L)


def outcomes_permuted(ids, seed):
    L = _read_outcome_table(ids)
    perm = np.random.default_rng([SEED, 7, seed]).permutation(len(L))
    assert (perm != np.arange(len(L))).mean() > 0.99
    return target(L[perm])            # whole rows moved to other neurons; nothing else leaves this function


def outcomes_real(ids):
    if not _RUN_OPEN:
        sys.exit("refusing: the real outcome table is opened only by a guarded `run`")
    L = _read_outcome_table(ids)
    nz = int((L.sum(1) == 0).sum())
    if nz != N_ALLZERO:
        sys.exit("refusing to score: %d all-zero profiles, the audit recorded %d" % (nz, N_ALLZERO))
    return target(L)


def synthetic(D, kind, seed):
    """SYNTHETIC outcomes built from inputs and a seed only (no table is read).
    tag: position + tag + animal effects; tag_moderate: the same with a weaker tag effect;
    notag: position + animal; animal_only: animal + noise."""
    r = np.random.default_rng([SEED, 11, seed])
    n, R = len(D.ids), len(REGIONS)
    Z = (D.X - D.X.mean(0)) / D.X.std(0)
    F = np.cos(Z @ r.normal(0, 1.0, (3, 24)) + r.uniform(0, 2 * np.pi, 24)) @ r.normal(0, 1, (24, R))
    F = (F - F.mean(0)) / F.std(0)
    s_pos, s_tag, s_an = {"tag": (1, .7, .5), "tag_moderate": (1, .3, .5), "notag": (1, 0, .7), "animal_only": (0, 0, 1)}[kind]
    ai = np.searchsorted(D.animals, D.animal)
    ti = np.searchsorted(D.tags, D.tagname)
    T = r.normal(0, 1, (N_TAGS, R))
    U = r.normal(0, 1, (N_ANIMALS, R))
    return s_pos * F + s_tag * T[ti] + s_an * U[ai] + r.normal(0, 1 if kind != "animal_only" else .5, (n, R))


# ------------------------------------------------------------------------------- splits (inputs only)
def deal(animals, tag_of, k, key):
    """Whole animals to k folds: within each tag (sorted) a seeded shuffle, then round-robin with one
    running offset over tags. Returns {animal: fold}."""
    out, off = {}, 0
    for t in sorted(set(tag_of[a] for a in animals)):
        A = sorted(a for a in animals if tag_of[a] == t)
        A = [A[i] for i in rng_for("deal", key, t).permutation(len(A))]
        for a in A:
            out[a] = off % k
            off += 1
    return out


def strata(D):
    """Animals to 4 spatial strata: AP centroid at or below / above the median over animals, then the same
    for DV inside each half. Used only by the stratified between-animal tag shuffle."""
    c = {a: D.X[D.animal == a].mean(0) for a in D.animals}
    ap = np.median([c[a][0] for a in D.animals])
    half = {a: int(c[a][0] > ap) for a in D.animals}
    st = {}
    for h in (0, 1):
        A = [a for a in D.animals if half[a] == h]
        dv = np.median([c[a][1] for a in A])
        for a in A:
            st[a] = 2 * h + int(c[a][1] > dv)
    return st


def s1_folds(D, tag_of):
    """The S1 rule as a function of the animal -> tag map: whole animals dealt to K_OUTER folds within tags.
    With the real tags this is the S1 split; the tag-shuffle control applies the same rule to shuffled tags."""
    idx = np.arange(len(D.ids))
    f1 = deal(list(D.animals), tag_of, K_OUTER, "S1")
    fold = np.array([f1[a] for a in D.animal])
    return [{"name": "fold%d" % k, "test": idx[fold == k].tolist(), "train": idx[fold != k].tolist()} for k in range(K_OUTER)]


def make_splits(D):
    n = len(D.ids)
    idx = np.arange(n)
    S = {"decision": "D-109", "seed": SEED, "n": n, "ids_sha256": hashlib.sha256("\n".join(D.ids).encode()).hexdigest(),
         "inputs": {k: sha256(path(INPUTS[k])) for k in ("soma", "labels", "ccf")}}
    S["S1"] = s1_folds(D, D.animal_tag)
    S2, skipped = [], []
    for key in sorted(set(D.key)):
        te = D.key == key
        anc = ancestors(key, D.par)
        na = len(set(D.animal[te]))
        if key == "HY" or FIBRE in anc or "HY" not in anc or te.sum() < S2_MIN_NEURONS or na < S2_MIN_ANIMALS:
            skipped.append([key, int(te.sum()), na])
            continue
        purge = np.isin(D.animal, list(set(D.animal[te])))        # every record of every animal with a test neuron
        S2.append({"name": key, "test": idx[te].tolist(), "train": idx[~purge].tolist(), "purged_not_tested": int((purge & ~te).sum())})
    S["S2"], S["S2_not_folds"] = S2, skipped
    S["S3"] = [{"name": t, "test": idx[D.tagname == t].tolist(), "train": idx[D.tagname != t].tolist()} for t in D.tags]
    S["animal_stratum"] = strata(D)
    return S


def dump_splits(S):
    return json.dumps(S, sort_keys=True, separators=(",", ":")) + "\n"


def assert_splits(D, S):
    """Leakage assertions on a split object; returns the list of checks made."""
    done = []
    for name in ("S1", "S2", "S3"):
        for f in S[name]:
            tr, te = np.array(f["train"]), np.array(f["test"])
            assert len(tr) and len(te) and not set(tr) & set(te)
            assert not set(D.animal[tr]) & set(D.animal[te]), "%s %s: an animal on both sides" % (name, f["name"])
            if name == "S1":
                assert len(tr) + len(te) == len(D.ids)
                assert set(D.tagname[te]) <= set(D.tagname[tr]), "S1: a test tag missing from training"
                assert set(D.tagname[tr]) == set(D.tags)
            if name == "S2":
                assert (D.key[te] == f["name"]).all() and (D.key == f["name"]).sum() == len(te)
                assert f["name"] not in set(D.key[tr])
                assert not any(f["name"] in ancestors(x, D.par) for x in set(D.fine[tr])), "S2: a sub-label in training"
            if name == "S3":
                assert (D.tagname[te] == f["name"]).all() and f["name"] not in set(D.tagname[tr])
                assert len(tr) + len(te) == len(D.ids)
        done.append("%s: %d folds, no animal on both sides" % (name, len(S[name])))
    t1 = np.concatenate([f["test"] for f in S["S1"]])
    assert len(t1) == len(set(t1)) == len(D.ids)
    t2 = np.concatenate([f["test"] for f in S["S2"]])
    assert len(t2) == len(set(t2))
    t3 = np.concatenate([f["test"] for f in S["S3"]])
    assert len(t3) == len(set(t3)) == len(D.ids)
    done.append("S1 and S3 test every neuron exactly once; S2 tests %d neurons once" % len(t2))
    return done


def inner_folds(tr, animal, tag, key):
    """K_INNER whole-animal folds inside a training set; `tag` is the tag vector the model is given."""
    A = sorted(set(animal[tr]))
    tag_of = {}
    for a, t in zip(animal[tr], tag[tr]):
        tag_of[a] = t
    d = deal(A, tag_of, K_INNER, "inner|" + key)
    f = np.array([d[a] for a in animal[tr]])
    return [(tr[f != k], tr[f == k]) for k in range(K_INNER)]


# ------------------------------------------------------------------------------- predictors
def onehot(tag_tr, tag_te):
    cats = np.array(sorted(set(tag_tr)))
    return (tag_tr[:, None] == cats).astype(float), (tag_te[:, None] == cats).astype(float)


def poly(Z):
    return np.column_stack([Z, Z ** 2, Z[:, 0] * Z[:, 1], Z[:, 0] * Z[:, 2], Z[:, 1] * Z[:, 2]])


def knn_mean(Xtr, Ytr, Xte, k):
    from scipy.spatial import cKDTree
    k = min(k, len(Xtr))
    _, j = cKDTree(Xtr).query(Xte, k=k)
    return Ytr[j.reshape(len(Xte), k)].mean(1)


def fit_predict(model, hp, Xtr, tag_tr, fine_tr, Ytr, Xte, tag_te, fine_te, par):
    """Everything a predictor may use is in the arguments: training inputs and outcomes, test inputs.
    All preprocessing (means, scales, categories, groups) is computed from the training arguments."""
    mu = Ytr.mean(0)
    if model == "global":
        return np.repeat(mu[None], len(Xte), 0)
    if model == "nucleus":
        out = np.empty((len(Xte), Ytr.shape[1]))
        anc_tr = {f: set(ancestors(f, par)) for f in set(fine_tr)}
        cache = {}
        for f in set(fine_te):
            a = ancestors(f, par)
            chain = [("exact", f)]
            if "HY" in a and f != "HY":
                chain += [("subtree", x) for x in a[1:a.index("HY") + 1]]
            val = mu
            for how, x in chain:
                m = fine_tr == x if how == "exact" else np.array([x in anc_tr[g] for g in fine_tr])
                if m.sum() >= NUCLEUS_MIN:
                    val = Ytr[m].mean(0)
                    break
            cache[f] = val
        for i, f in enumerate(fine_te):
            out[i] = cache[f]
        return out
    if model == "knn":
        return knn_mean(Xtr, Ytr, Xte, hp)
    if model == "knn_tag":
        out = knn_mean(Xtr, Ytr, Xte, hp)                       # fallback: neighbours of any tag
        for t in set(tag_te):
            m = tag_tr == t
            if m.sum() >= hp:                                   # enough same-tag training neurons
                out[tag_te == t] = knn_mean(Xtr[m], Ytr[m], Xte[tag_te == t], hp)
        return out
    if model.endswith("_tag"):
        unseen = ~np.isin(tag_te, tag_tr)
        if unseen.any():                                        # back-off: the position-only model, same setting
            out = np.empty((len(Xte), Ytr.shape[1]))
            out[unseen] = fit_predict(model[:-4], hp, Xtr, tag_tr, fine_tr, Ytr, Xte[unseen], tag_te[unseen], fine_te[unseen], par)
            if (~unseen).any():
                out[~unseen] = fit_predict(model, hp, Xtr, tag_tr, fine_tr, Ytr, Xte[~unseen], tag_te[~unseen], fine_te[~unseen], par)
            return out
    xm, xs = Xtr.mean(0), Xtr.std(0)
    Ztr, Zte = (Xtr - xm) / xs, (Xte - xm) / xs
    Otr, Ote = onehot(tag_tr, tag_te)
    if model in ("ridge", "ridge_tag"):
        Ptr, Pte = poly(Ztr), poly(Zte)
        pm, ps = Ptr.mean(0), Ptr.std(0)
        Ftr, Fte = (Ptr - pm) / ps, (Pte - pm) / ps
        if model == "ridge_tag":
            om = Otr.mean(0)
            Ftr, Fte = np.column_stack([Ftr, Otr - om]), np.column_stack([Fte, Ote - om])
        B = np.linalg.solve(Ftr.T @ Ftr + hp * np.eye(Ftr.shape[1]), Ftr.T @ (Ytr - mu))
        return mu + Fte @ B
    if model in ("et", "et_tag"):
        from sklearn.ensemble import ExtraTreesRegressor
        Ftr, Fte = (Xtr, Xte) if model == "et" else (np.column_stack([Xtr, Otr]), np.column_stack([Xte, Ote]))
        m = ExtraTreesRegressor(n_estimators=ET_TREES, min_samples_leaf=hp, max_features=1.0, n_jobs=-1, random_state=SEED)
        return m.fit(Ftr, Ytr).predict(Fte)
    if model in ("mlp", "mlp_tag"):
        from sklearn.neural_network import MLPRegressor
        Ftr, Fte = (Ztr, Zte) if model == "mlp" else (np.column_stack([Ztr, Otr]), np.column_stack([Zte, Ote]))
        sd = np.maximum(Ytr.std(0), 1e-6)
        m = MLPRegressor(hidden_layer_sizes=MLP_HIDDEN, alpha=hp, max_iter=MLP_EPOCHS, batch_size=128, early_stopping=False,
                         learning_rate_init=1e-3, random_state=SEED, n_iter_no_change=10 ** 6, tol=0.0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return mu + sd * m.fit(Ftr, (Ytr - mu) / sd).predict(Fte)
    raise ValueError(model)


def predict(model, hp, tr, te, D, Y, tag):
    return fit_predict(model, hp, D.X[tr], tag[tr], D.fine[tr], Y[tr], D.X[te], tag[te], D.fine[te], D.par)


def animal_mean(e, idx, animal):
    """Mean of per-neuron values within each animal: (sorted animals, means)."""
    A = np.array(sorted(set(animal[idx])))
    return A, np.array([e[idx][animal[idx] == a].mean() for a in A])


def tune(model, tr, D, Y, tag, key):
    """Nested whole-animal tuning: the value minimising the primary metric inside the training set, i.e.
    the mean over ALL training animals (each is an inner test animal exactly once) of the per-animal mean
    squared error; ties to the first value in GRID order."""
    grid = GRID[model]
    if len(grid) == 1:
        return grid[0]
    full = np.full((len(grid), len(D.ids)), np.nan)
    for itr, ite in inner_folds(tr, D.animal, tag, key):
        for g, hp in enumerate(grid):
            full[g, ite] = ((predict(model, hp, itr, ite, D, Y, tag) - Y[ite]) ** 2).mean(1)
    loss = np.array([animal_mean(full[g], tr, D.animal)[1].mean() for g in range(len(grid))])
    assert np.isfinite(loss).all()
    return grid[int(np.argmin(loss))]


def run_split(D, Y, folds, models, tag, split):
    P = {m: np.full(Y.shape, np.nan) for m in models}
    hps = {m: {} for m in models}
    for f in folds:
        tr, te = np.array(f["train"]), np.array(f["test"])
        for m in models:
            hp = tune(m, tr, D, Y, tag, split + "|" + f["name"])
            P[m][te] = predict(m, hp, tr, te, D, Y, tag)
            hps[m][f["name"]] = hp
    return P, hps


# ------------------------------------------------------------------------------- scoring
def neuron_mse(P, Y, cols=None):
    d = (P - Y) ** 2
    return (d if cols is None else d[:, cols]).mean(1)


def neuron_r(P, Y):
    """Pearson r over the 49 regions; nan when the observed profile is constant, 0 when only the prediction is."""
    Pc, Yc = P - P.mean(1, keepdims=True), Y - Y.mean(1, keepdims=True)
    sp, sy = np.sqrt((Pc ** 2).sum(1)), np.sqrt((Yc ** 2).sum(1))
    r = np.where(sp > 0, (Pc * Yc).sum(1) / np.where(sp * sy > 0, sp * sy, 1), 0.0)
    return np.where(sy > 0, r, np.nan)


def signflip_p(d, key):
    """One-sided (mean > 0) sign-flip p over units, N_FLIP random flips, +1 correction."""
    s = rng_for("flip", key).choice([-1.0, 1.0], (N_FLIP, len(d)))
    return float((1 + ((s * d).mean(1) >= d.mean() - 1e-15).sum()) / (N_FLIP + 1))


def boot_ci(num, den, key):
    """Percentile intervals, bootstrap over units: 95% of mean(num); 95% and 99% of mean(num) / mean(den)."""
    j = rng_for("boot", key).integers(0, len(num), (N_BOOT, len(num)))
    a, b = num[j].mean(1), den[j].mean(1)
    q = lambda v, lo, hi: [float(x) for x in np.percentile(v, [lo, hi])]
    return q(a, 2.5, 97.5), q(a / b, 2.5, 97.5), q(a / b, 0.5, 99.5)


def unit_mean(v, idx, animal, fold_of=None):
    """The reported aggregate: mean over test animals of the per-animal mean (S1); with `fold_of`
    (neuron index -> fold name; S2, S3) the mean over folds of that quantity within each fold."""
    if fold_of is None:
        return float(animal_mean(v, idx, animal)[1].mean())
    names = sorted(set(fold_of[i] for i in idx))
    return float(np.mean([animal_mean(v, np.array([i for i in idx if fold_of[i] == u]), animal)[1].mean() for u in names]))


def compare(e_new, e_ref, idx, animal, key, units=None):
    """Paired comparison `ref minus new` (positive = `new` is better), aggregated per test animal and then,
    when `units` maps neuron index -> fold name, per fold. Returns the statistics dict."""
    if units is None:
        _, a_new = animal_mean(e_new, idx, animal)
        _, a_ref = animal_mean(e_ref, idx, animal)
        unit = "animal"
    else:
        names = sorted(set(units[i] for i in idx))
        a_new = np.array([animal_mean(e_new, np.array([i for i in idx if units[i] == u]), animal)[1].mean() for u in names])
        a_ref = np.array([animal_mean(e_ref, np.array([i for i in idx if units[i] == u]), animal)[1].mean() for u in names])
        unit = "fold"
    d = a_ref - a_new
    ci, ci_rel, ci99_rel = boot_ci(d, a_ref, key)
    return {"unit": unit, "n_units": int(len(d)), "err_new": float(a_new.mean()), "err_ref": float(a_ref.mean()),
            "diff": float(d.mean()), "rel": float(d.mean() / a_ref.mean()), "ci95": ci, "ci95_rel": ci_rel, "ci99_rel": ci99_rel,
            "p_signflip": signflip_p(d, key), "frac_units_better": float((d > 0).mean())}


def shuffled_tags(D, S, b, stratified):
    """Tags permuted BETWEEN animals (each animal's neurons keep one common tag); within the 4 spatial strata
    when `stratified`. Returns the per-neuron tag vector."""
    r = rng_for("shuffle", stratified, b)
    new = {}
    groups = [[a for a in D.animals if S["animal_stratum"][a] == s] for s in range(4)] if stratified else [list(D.animals)]
    for A in groups:
        t = [D.animal_tag[a] for a in A]
        for a, j in zip(A, r.permutation(len(A))):
            new[a] = t[j]
    return np.array([new[a] for a in D.animal])


def primary_stat(D, Y, tag_of, tg):
    """The primary statistic as a function of a tag assignment: S1 folds dealt from `tag_of` (s1_folds), both
    arms tuned and fitted with the per-neuron tag vector `tg`, mean over animals of error(et) - error(et_tag).
    Returns (difference, mean error of the reference arm)."""
    idx = np.arange(len(D.ids))
    P, _ = run_split(D, Y, s1_folds(D, tag_of), list(PRIMARY), tg, "S1")
    a1 = animal_mean(neuron_mse(P[PRIMARY[0]], Y), idx, D.animal)[1]
    a0 = animal_mean(neuron_mse(P[PRIMARY[1]], Y), idx, D.animal)[1]
    return float((a0 - a1).mean()), float(a0.mean())


def tag_shuffle(D, Y, S, obs, n_perm, stratified):
    """Null distribution of the primary statistic under between-animal tag permutation. The WHOLE procedure
    is repeated on each shuffled assignment: S1 folds are dealt again from the shuffled tags (so every tag is
    in every training fold, as in the observed arrangement, and no test animal falls back), and both arms
    are tuned and fitted again."""
    null = []
    for b in range(n_perm):
        tg = shuffled_tags(D, S, b, stratified)
        tag_of = dict(zip(D.animal, tg))
        null.append(primary_stat(D, Y, tag_of, tg)[0])
    null = np.array(null)
    if not n_perm:
        return {"n_perm": 0, "stratified": stratified, "p": float("nan"), "null_mean": None, "null_sd": None, "null_max": None,
                "obs_minus_null_mean_in_sd": None, "null": []}
    sd = float(null.std(ddof=1)) if n_perm > 1 else float("nan")
    return {"n_perm": n_perm, "stratified": stratified, "p": float((1 + (null >= obs - 1e-15).sum()) / (n_perm + 1)),
            "null_mean": float(null.mean()), "null_sd": sd, "null_max": float(null.max()),
            "obs_minus_null_mean_in_sd": float((obs - null.mean()) / sd) if sd > 0 else None, "null": [float(x) for x in null]}


def decide(c, sh, planted_ok):
    """The pre-registered reading of the primary comparison (protocol section 12)."""
    if c["ci99_rel"][1] < MIN_EFFECT:
        return ("NOT SHOWN FOR THIS MODEL: extra-trees did not gain %.0f%% from the tag one-hot (upper end of the 99%% "
                "interval below it); this does not bound the information in the tag" % (100 * MIN_EFFECT))
    if (c["rel"] >= MIN_EFFECT and c["ci95"][0] > 0 and c["p_signflip"] < 0.025 and sh["p"] < 0.05 and planted_ok):
        return "SHOWN: the tag improves prediction for new animals of represented lines (predictive, not causal)"
    return "INSUFFICIENT EVIDENCE"


def versions():
    import pandas
    import scipy
    import sklearn
    return {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
            "pandas": pandas.__version__, "scikit-learn": sklearn.__version__}


def unseen_tag_animals(D, f):
    """Number of test animals of a fold whose tag is absent from its training set (inputs only)."""
    tr, te = np.array(f["train"]), np.array(f["test"])
    seen = set(D.tagname[tr])
    return len(set(a for a, t in zip(D.animal[te], D.tagname[te]) if t not in seen))


def evaluate(D, Y, S, n_perm, n_perm_unstrat, planted_ok, tagfile, marker=None):
    """The whole evaluation. Called by `dry` (PERMUTED outcomes) and by `run`."""
    t0 = time.time()
    assert LABEL in ("PERMUTED", "SYNTHETIC") or (LABEL == "REAL" and _RUN_OPEN)
    out = os.path.join(INTERIM, "d109_%s_results.json" % tagfile)
    out_npz = os.path.join(INTERIM, "d109_%s_predictions.npz" % tagfile)
    if LABEL == "REAL" and (os.path.exists(out) or os.path.exists(out_npz)):
        sys.exit("refusing: %s exists; results of a real run are never overwritten" % os.path.relpath(out, ROOT))
    tag = D.tagname
    res = {"label": LABEL, "splits": {}, "n_perm": n_perm, "versions": versions(), "marker": marker, "grid": GRID}
    say("library versions:", res["versions"])
    preds = {}
    for split in ("S1", "S2", "S3"):
        models = MODELS_OF[split]
        P, hps = run_split(D, Y, S[split], models, tag, split)
        idx = np.concatenate([f["test"] for f in S[split]])
        fold_of = None if split == "S1" else {i: f["name"] for f in S[split] for i in f["test"]}
        E = {m: neuron_mse(P[m], Y) for m in models}
        rows = {}
        for m in models:
            full_r = np.full(len(D.ids), np.nan)
            full_r[idx] = neuron_r(P[m][idx], Y[idx])
            ok = idx[np.isfinite(full_r[idx])]
            vals = list(hps[m].values())
            edge = sum(v in (GRID[m][0], GRID[m][-1]) for v in vals) if len(GRID[m]) > 1 else 0
            rows[m] = {"unit_mse": unit_mean(E[m], idx, D.animal, fold_of), "neuron_mse": float(E[m][idx].mean()),
                       "unit_r": unit_mean(full_r, ok, D.animal, fold_of), "hp": hps[m], "n_folds_at_grid_edge": int(edge)}
        g = rows["global"]["unit_mse"]
        say("%s  test neurons %d, test animals %d, folds %d   [%.0f s]" % (split, len(idx), len(set(D.animal[idx])), len(S[split]), time.time() - t0))
        say("   aggregation: %s" % ("mean over test animals" if split == "S1" else "mean over test animals within each fold, then mean over folds"))
        say("   %-10s %12s %8s %10s   hyperparameters chosen (count over folds); folds at a grid edge" % ("model", "MSE", "skill", "r"))
        for m in models:
            vals = list(rows[m]["hp"].values())
            say("   %-10s %12.4f %8.3f %10.3f   %s; %d" % (m, rows[m]["unit_mse"], 1 - rows[m]["unit_mse"] / g, rows[m]["unit_r"],
                                                         {str(v): vals.count(v) for v in sorted(set(vals), key=str)}, rows[m]["n_folds_at_grid_edge"]))
        comps = {}
        for name, new, ref in COMPARISONS[split]:
            comps[name] = c = compare(E[new], E[ref], idx, D.animal, split + name, fold_of)
            say("   %-34s diff %+.4f (rel %+.2f%%, 95%% CI %+.4f..%+.4f, rel %+.2f%%..%+.2f%%), sign-flip p %.4f, n=%d %ss" % (
                name, c["diff"], 100 * c["rel"], c["ci95"][0], c["ci95"][1], 100 * c["ci95_rel"][0], 100 * c["ci95_rel"][1],
                c["p_signflip"], c["n_units"], c["unit"]))
        per_fold = {}
        if split != "S1":
            for f in S[split]:
                te = np.array(f["test"])
                per_fold[f["name"]] = {"n_neurons": len(te), "n_animals": len(set(D.animal[te])), "n_train": len(f["train"]),
                                       "n_animals_unseen_tag": unseen_tag_animals(D, f),
                                       **{m: float(animal_mean(E[m], te, D.animal)[1].mean()) for m in models}}
            say("   per fold (animal MSE; n = neurons/animals/animals whose tag is not in training): " + "; ".join("%s n=%d/%d/%d global %.3f et %.3f" % (
                k, v["n_neurons"], v["n_animals"], v["n_animals_unseen_tag"], v["global"], v["et"]) for k, v in per_fold.items()))
        res["splits"][split] = {"models": rows, "comparisons": comps, "per_fold": per_fold}
        preds[split] = P
        if split == "S1":
            res["sensitivity"], res["by_tag_size"], res["n_allzero_test"] = s1_extras(D, Y, P, E, idx)
    # primary comparison and its tag-shuffle controls
    c = res["splits"]["S1"]["comparisons"]["PRIMARY tag benefit (et_tag vs et)"]
    for strat, n, key in ((True, n_perm, "tag_shuffle"), (False, n_perm_unstrat, "tag_shuffle_unstratified")):
        res[key] = sh = tag_shuffle(D, Y, S, c["diff"], n, strat)
        say("tag shuffle between animals, %s: %d permutations (folds dealt again, both arms refitted), p %.4f, observed %+.4f, null mean %+.4f, null SD %.4f, "
            "null max %+.4f, observed minus null mean = %s null SD   [%.0f s]" % (
                "within 4 spatial strata" if strat else "unstratified (secondary)", sh["n_perm"], sh["p"], c["diff"], sh["null_mean"], sh["null_sd"],
                sh["null_max"], "%+.2f" % sh["obs_minus_null_mean_in_sd"] if sh["obs_minus_null_mean_in_sd"] is not None else "n/a", time.time() - t0))
    res["decision"] = decide(c, res["tag_shuffle"], planted_ok)
    shown = res["decision"].startswith("SHOWN")
    # Holm over the declared secondary family
    fam = [(s, n) for s in ("S1", "S2", "S3") for n, _, _ in COMPARISONS[s] if n.startswith("H")]
    ps = sorted((res["splits"][s]["comparisons"][n]["p_signflip"], s, n) for s, n in fam)
    adj, run_max = {}, 0.0
    for i, (p, s, n) in enumerate(ps):
        run_max = max(run_max, min(1.0, p * (len(ps) - i)))
        rel = res["splits"][s]["comparisons"][n]["rel"]
        adj[s + " " + n] = {"p": p, "p_holm": run_max, "rel": rel, "claim": bool(run_max < 0.05 and rel >= MIN_EFFECT)}
    gate = {"S1 H1 et vs nucleus mean": "H0", "S1 H2 et vs knn": "H0", "S1 H3 mlp_tag vs mlp": "PRIMARY SHOWN", "S2 H5 et_tag vs et": "PRIMARY SHOWN"}
    h0 = adj["S1 H0 et vs global mean"]["claim"]
    for k, need in gate.items():                 # H1, H2 only if position predicts at all; tag comparisons only if the primary is SHOWN
        adj[k]["claim"] = adj[k]["claim"] and (h0 if need == "H0" else shown)
        adj[k]["gate"] = need
    for k, v in adj.items():
        say("   Holm (family of %d)  %-38s p %.4f  adjusted %.4f  rel %+.2f%%  claim (adjusted p < 0.05 and rel >= %.0f%%%s): %s" % (
            len(ps), k, v["p"], v["p_holm"], 100 * v["rel"], 100 * MIN_EFFECT, ", and " + v["gate"] if "gate" in v else "", "yes" if v["claim"] else "no"))
    res["holm"] = adj
    sh = res["tag_shuffle"]
    say("PRIMARY (S1, et_tag vs et): relative benefit %+.2f%% (95%% CI %+.2f%%..%+.2f%%; 99%% CI %+.2f%%..%+.2f%%), sign-flip p %.4f over %d animals, tag-shuffle p %.4f" % (
        100 * c["rel"], 100 * c["ci95_rel"][0], 100 * c["ci95_rel"][1], 100 * c["ci99_rel"][0], 100 * c["ci99_rel"][1], c["p_signflip"], c["n_units"], sh["p"]))
    say("DECISION RULE OUTPUT: " + res["decision"] + ("" if LABEL == "REAL" else "   (%s DATA: NOT A RESULT)" % LABEL))
    json.dump(res, open(out, "w"), indent=1, sort_keys=True)
    np.savez_compressed(out_npz, **{"%s_%s" % (s, m): preds[s][m].astype(np.float32) for s in preds for m in MODELS_OF[s]})
    say("written", os.path.relpath(out, ROOT), "  total %.0f s" % (time.time() - t0))
    return res


def s1_extras(D, Y, P, E, idx):
    """Declared sensitivity analyses and descriptive breakdowns of the primary comparison (S1). None is decisive."""
    new, ref = PRIMARY
    A, a_new = animal_mean(E[new], idx, D.animal)
    _, a_ref = animal_mean(E[ref], idx, D.animal)
    d = a_ref - a_new
    nz = idx[(Y[idx] != 0).any(1)]
    cols = [i for i, r in enumerate(REGIONS) if r != "MY"]
    top = A[np.argsort(-np.abs(d), kind="stable")[:5]]
    sens = {"no_allzero_test": (compare(E[new], E[ref], nz, D.animal, "sens_nz"), len(nz)),
            "no_MY_column_same_models": (compare(neuron_mse(P[new], Y, cols), neuron_mse(P[ref], Y, cols), idx, D.animal, "sens_my"), len(idx))}
    keep = idx[~np.isin(D.animal[idx], top)]
    sens["without_5_largest_abs_d_animals"] = (compare(E[new], E[ref], keep, D.animal, "sens_top5"), len(keep))
    keep = idx[~np.isin(D.animal[idx], D.swc_animals)]
    sens["without_32_swc_sample_animals"] = (compare(E[new], E[ref], keep, D.animal, "sens_swc"), len(keep))
    for k, (c, n) in sens.items():
        say("   sensitivity %-32s diff %+.4f (rel %+.2f%%, 95%% CI rel %+.2f%%..%+.2f%%), n=%d animals, test neurons %d" % (
            k, c["diff"], 100 * c["rel"], 100 * c["ci95_rel"][0], 100 * c["ci95_rel"][1], c["n_units"], n))
    say("   the 5 animals with the largest |d_a| (animal, tag, d_a): " + "; ".join(
        "%s %s %+.4f" % (a, D.animal_tag[a], d[list(A).index(a)]) for a in top))
    size = {t: sum(1 for a in D.animals if D.animal_tag[a] == t) for t in D.tags}
    groups = {}
    for lo, hi in TAG_SIZE_GROUPS:
        m = np.array([lo <= size[D.animal_tag[a]] <= hi for a in A])
        if m.any():
            groups["%d-%d" % (lo, hi)] = {"n_animals": int(m.sum()), "n_tags": len(set(D.animal_tag[a] for a in A[m])),
                                          "diff": float(d[m].mean()), "rel": float(d[m].mean() / a_ref[m].mean())}
    say("   primary benefit by tag size (animals per tag: animals, tags, rel; descriptive, no per-tag claim below 5 animals): " + "; ".join(
        "%s: %d, %d, %+.2f%%" % (k, v["n_animals"], v["n_tags"], 100 * v["rel"]) for k, v in groups.items()))
    return {k: c for k, (c, n) in sens.items()}, groups, int(len(idx) - len(nz))


# name, new model, reference model. Names starting with "H" form the secondary (Holm) family.
COMPARISONS = {
    "S1": [("PRIMARY tag benefit (et_tag vs et)", "et_tag", "et"),
           ("H0 et vs global mean", "et", "global"),
           ("H1 et vs nucleus mean", "et", "nucleus"),
           ("H2 et vs knn", "et", "knn"),
           ("H3 mlp_tag vs mlp", "mlp_tag", "mlp"),
           ("D knn_tag vs knn", "knn_tag", "knn"),
           ("D ridge_tag vs ridge", "ridge_tag", "ridge"),
           ("D mlp vs nucleus mean", "mlp", "nucleus")],
    "S2": [("H4 et vs global mean", "et", "global"),
           ("H5 et_tag vs et", "et_tag", "et"),
           ("D et vs nucleus back-off mean", "et", "nucleus"),
           ("D mlp vs global mean", "mlp", "global")],
    "S3": [("H6 et vs global mean", "et", "global"),
           ("D et vs nucleus mean", "et", "nucleus"),
           ("D mlp vs global mean", "mlp", "global")],
}


# ------------------------------------------------------------------------------- modes
def load_splits(D):
    if not os.path.exists(SPLITS):
        sys.exit("no split file: run `splits` first")
    txt = open(SPLITS).read()
    if txt != dump_splits(make_splits(D)):
        sys.exit("refusing: data/interim/d109_splits.json is not what `splits` produces from the inputs")
    return json.loads(txt)


def mode_splits():
    D = load_inputs()
    S = make_splits(D)
    for line in assert_splits(D, S):
        print(line)
    open(SPLITS, "w").write(dump_splits(S))
    print("S1 folds (neurons / animals):", [(len(f["test"]), len(set(D.animal[f["test"]]))) for f in S["S1"]])
    print("S2 hold-out nuclei: test neurons / test animals / training neurons / training animals / training tags / purged, not tested / "
          "test animals whose tag is not in training")
    for f in S["S2"]:
        te, tr = np.array(f["test"]), np.array(f["train"])
        print("   %-5s %4d %3d %5d %4d %3d %5d %3d" % (f["name"], len(te), len(set(D.animal[te])), len(tr), len(set(D.animal[tr])),
                                                     len(set(D.tagname[tr])), f["purged_not_tested"], unseen_tag_animals(D, f)))
    print("S2 labels that are not folds (label, neurons, animals):", S["S2_not_folds"])
    print("S3 tags (neurons / animals):", [(f["name"], len(f["test"]), len(set(D.animal[f["test"]]))) for f in S["S3"]])
    st = S["animal_stratum"]
    for s in range(4):
        A = [a for a in D.animals if st[a] == s]
        tg = [D.animal_tag[a] for a in A]
        print("tag-shuffle stratum %d: %d animals, tags %s" % (s, len(A), {t: tg.count(t) for t in sorted(set(tg))}))
    # inputs-only description of the tag / animal-ID confound (animal ID order is the available proxy for batch)
    order = sorted(D.animals, key=lambda a: (len(a), a))
    seq = [D.animal_tag[a] for a in order]
    runs = lambda q: 1 + sum(x != y for x, y in zip(q, q[1:]))
    r = rng_for("id-order")
    sim = [runs([seq[i] for i in r.permutation(len(seq))]) for _ in range(2000)]
    print("tags along the animal-ID order: %d runs of equal tag over %d animals (%d tags; random assignment: mean %.0f, minimum of 2000 draws %d)" % (
        runs(seq), len(seq), N_TAGS, np.mean(sim), min(sim)))
    for t in D.tags:
        A = [a for a in order if D.animal_tag[a] == t]
        print("   %-11s %2d animals, IDs %s to %s, in %d runs" % (t, len(A), A[0], A[-1], sum(1 for i, x in enumerate(seq) if x == t and (i == 0 or seq[i - 1] != t))))
    print("SWC sample (D-110 / D-111): %d neurons in %d animals; S1 sensitivity keeps %d animals" % (
        len(SWC_SAMPLE), len(D.swc_animals), N_ANIMALS - len(D.swc_animals)))
    print("written %s  sha256 %s" % (os.path.relpath(SPLITS, ROOT), sha256(SPLITS)))


def mode_check():
    global LABEL
    D = load_inputs()
    S = load_splits(D)
    ok = []

    def verdict(name, cond, detail=""):
        ok.append(bool(cond))
        say("%s  %s  %s" % ("PASS" if cond else "FAIL", name, detail))
    LABEL = "CHECK"
    say("library versions:", versions())
    for line in assert_splits(D, S):
        say("PASS  split assertion:", line)
    # inner folds: whole animals
    for split in ("S1", "S2", "S3"):
        for f in S[split]:
            for itr, ite in inner_folds(np.array(f["train"]), D.animal, D.tagname, split + "|" + f["name"]):
                assert not set(D.animal[itr]) & set(D.animal[ite]) and set(itr) | set(ite) == set(f["train"])
    say("PASS  inner tuning folds are whole-animal partitions of the training set, for every outer fold")
    # tag shuffle keeps animals whole and the tag multiset fixed within strata
    for strat in (True, False):
        tg = shuffled_tags(D, S, 0, strat)
        assert all(len(set(tg[D.animal == a])) == 1 for a in D.animals)
        if strat:
            for s in range(4):
                A = [a for a in D.animals if S["animal_stratum"][a] == s]
                assert sorted(tg[D.animal == a][0] for a in A) == sorted(D.animal_tag[a] for a in A)
        assert sorted(tg[D.animal == a][0] for a in D.animals) == sorted(D.animal_tag.values())
        for b in range(20):                                        # folds dealt again from shuffled tags: same guarantees as S1
            tgb = shuffled_tags(D, S, b, strat)
            for f in s1_folds(D, dict(zip(D.animal, tgb))):
                tr, te = np.array(f["train"]), np.array(f["test"])
                assert not set(D.animal[tr]) & set(D.animal[te]) and set(tgb[te]) <= set(tgb[tr]) and len(set(tgb[tr])) == N_TAGS
        say("PASS  tag shuffle (stratified=%s): one tag per animal, tag counts preserved; %d of %d animals change tag in permutation 0; "
            "in 20 permutations the re-dealt folds keep animals whole and every tag in every training fold" % (
                strat, sum(tg[D.animal == a][0] != D.animal_tag[a] for a in D.animals), N_ANIMALS))
    assert s1_folds(D, D.animal_tag) == S["S1"]
    say("PASS  the fold rule of the tag shuffle applied to the real tags reproduces S1 exactly")
    LABEL = "SYNTHETIC"
    # test-invariance: a test neuron's prediction does not depend on which other test neurons are present
    Y = synthetic(D, "tag", 0)
    f = S["S1"][0]
    tr, te = np.array(f["train"]), np.array(f["test"])
    worst = 0.0
    for m in MODELS:
        hp = GRID[m][len(GRID[m]) // 2]
        a = predict(m, hp, tr, te, D, Y, D.tagname)
        b = predict(m, hp, tr, te[::7], D, Y, D.tagname)
        worst = max(worst, float(abs(a[::7] - b).max()))
        Y2 = Y.copy()
        Y2[te] = 1e6                                    # test outcomes must not reach any prediction
        c = predict(m, hp, tr, te, D, Y2, D.tagname)
        worst = max(worst, float(abs(a - c).max()))
    verdict("test-invariance of all %d predictors (subset of test inputs; corrupted test outcomes)" % len(MODELS), worst < 1e-8, "max change %.2e" % worst)
    # fallbacks: unseen tag and unseen nucleus
    f = S["S3"][0]
    tr, te = np.array(f["train"]), np.array(f["test"])
    d = max(float(abs(predict(m + "_tag", GRID[m][1], tr, te, D, Y, D.tagname) - predict(m, GRID[m][1], tr, te, D, Y, D.tagname)).max())
            for m in ("knn", "ridge", "et", "mlp"))
    verdict("unseen tag (S3 fold %s): every with-tag model returns exactly its position-only counterpart" % f["name"], d < 1e-8, "max difference %.2e" % d)
    f = S["S2"][0]
    tr, te = np.array(f["train"]), np.array(f["test"])
    p = predict("nucleus", None, tr, te, D, Y, D.tagname)
    verdict("unseen nucleus: nucleus-mean baseline backs off (finite, one value per fine label) on S2 fold %s" % f["name"],
            np.isfinite(p).all() and len(np.unique(p[:, 0])) <= len(set(D.fine[te])))
    # leak detector: outcomes that depend on the animal only
    Y = synthetic(D, "animal_only", 1)
    quick = ["global", "nucleus", "knn", "ridge_tag", "et", "et_tag", "mlp_tag"]
    for split in ("S1", "S2", "S3"):
        P, _ = run_split(D, Y, S[split], quick, D.tagname, split)
        idx = np.concatenate([f["test"] for f in S[split]])
        sk = {m: 1 - animal_mean(neuron_mse(P[m], Y), idx, D.animal)[1].mean() / animal_mean(neuron_mse(P["global"], Y), idx, D.animal)[1].mean() for m in quick}
        verdict("animal-only outcomes, %s: no model beats the global mean by more than 1%%" % split, max(sk.values()) < 0.01,
                "skill " + ", ".join("%s %+.3f" % kv for kv in sk.items()))
    leaky = [{"name": f["name"], "test": f["test"], "train": sorted(set(range(len(D.ids))) - set(f["test"]))} for f in S["S2"][:6]]
    P, _ = run_split(D, Y, leaky, ["global", "knn", "et"], D.tagname, "leaky")
    idx = np.concatenate([f["test"] for f in leaky])
    sk = {m: 1 - animal_mean(neuron_mse(P[m], Y), idx, D.animal)[1].mean() / animal_mean(neuron_mse(P["global"], Y), idx, D.animal)[1].mean() for m in ("knn", "et")}
    verdict("positive control of the detector: S2 WITHOUT the animal purge (first 6 nuclei) breaks the 1% criterion", max(sk.values()) > 0.01,
            "skill " + ", ".join("%s %+.3f" % kv for kv in sk.items()))
    fold = rng_for("leaky-neuron-folds").integers(0, K_OUTER, len(D.ids))        # neurons, not animals, dealt to folds
    leaky = [{"name": "n%d" % k, "test": np.flatnonzero(fold == k).tolist(), "train": np.flatnonzero(fold != k).tolist()} for k in range(K_OUTER)]
    P, _ = run_split(D, Y, leaky, ["global", "knn", "et"], D.tagname, "leaky2")
    idx = np.arange(len(D.ids))
    sk = {m: 1 - animal_mean(neuron_mse(P[m], Y), idx, D.animal)[1].mean() / animal_mean(neuron_mse(P["global"], Y), idx, D.animal)[1].mean() for m in ("knn", "et")}
    verdict("positive control of the detector: neuron-level folds (animals on both sides) break the 1% criterion", max(sk.values()) > 0.01,
            "skill " + ", ".join("%s %+.3f" % kv for kv in sk.items()))
    # planted controls of the primary comparison
    idx = np.arange(len(D.ids))
    out = {}
    for kind, seed in (("tag", 2), ("notag", 3)):
        Y = synthetic(D, kind, seed)
        P, _ = run_split(D, Y, S["S1"], list(PRIMARY), D.tagname, "S1")
        e1, e0 = neuron_mse(P[PRIMARY[0]], Y), neuron_mse(P[PRIMARY[1]], Y)
        c = compare(e1, e0, idx, D.animal, "planted" + kind)
        if kind == "tag":                               # the shuffle's statistic, applied to the real tags, is the observed one
            assert abs(primary_stat(D, Y, D.animal_tag, D.tagname)[0] - c["diff"]) < 1e-12
        sh = tag_shuffle(D, Y, S, c["diff"], CHECK_SHUFFLE, True)
        out[kind] = {"compare": c, "shuffle": {k: v for k, v in sh.items()}, "decision": decide(c, sh, True)}
        say("planted '%s': tag benefit %+.4f (rel %+.2f%%, 95%% CI rel %+.2f%%..%+.2f%%, 99%% CI rel %+.2f%%..%+.2f%%), sign-flip p %.4f, shuffle p %.3f "
            "(%d permutations, smallest possible %.3f; null mean %+.4f, SD %.4f, observed minus null mean %+.2f SD) -> %s" % (
                kind, c["diff"], 100 * c["rel"], 100 * c["ci95_rel"][0], 100 * c["ci95_rel"][1], 100 * c["ci99_rel"][0], 100 * c["ci99_rel"][1],
                c["p_signflip"], sh["p"], CHECK_SHUFFLE, 1 / (CHECK_SHUFFLE + 1), sh["null_mean"], sh["null_sd"], sh["obs_minus_null_mean_in_sd"],
                out[kind]["decision"]))
    # moderate planted tag effect: printed, NOT a gate (shows how much of a real tag effect each model realises)
    Y = synthetic(D, "tag_moderate", 4)
    P, _ = run_split(D, Y, S["S1"], ["et", "et_tag", "ridge", "ridge_tag"], D.tagname, "S1")
    mod = {a: compare(neuron_mse(P[a + "_tag"], Y), neuron_mse(P[a], Y), idx, D.animal, "moderate" + a) for a in ("et", "ridge")}
    out["tag_moderate_not_a_gate"] = mod
    say("planted moderate tag effect (tag SD 0.3 of noise SD; NOT A GATE): benefit realised by et_tag %+.2f%% (95%% CI %+.2f%%..%+.2f%%), by ridge_tag %+.2f%% "
        "(95%% CI %+.2f%%..%+.2f%%); no shuffle run" % (100 * mod["et"]["rel"], 100 * mod["et"]["ci95_rel"][0], 100 * mod["et"]["ci95_rel"][1],
                                                     100 * mod["ridge"]["rel"], 100 * mod["ridge"]["ci95_rel"][0], 100 * mod["ridge"]["ci95_rel"][1]))
    verdict("planted tag effect is detected", out["tag"]["decision"].startswith("SHOWN"))
    verdict("no planted tag effect (position and animal effects only) is not called SHOWN", not out["notag"]["decision"].startswith("SHOWN"))
    # the same pipeline on PERMUTED outcomes, primary pair only (the full path is `dry`)
    LABEL = "PERMUTED"
    Y = outcomes_permuted(D.ids, 0)
    P, _ = run_split(D, Y, S["S1"], ["global"] + list(PRIMARY), D.tagname, "S1")
    c = compare(neuron_mse(P[PRIMARY[0]], Y), neuron_mse(P[PRIMARY[1]], Y), idx, D.animal, "permcheck")
    g = compare(neuron_mse(P["et"], Y), neuron_mse(P["global"], Y), idx, D.animal, "permcheck2")
    say("outcome rows permuted across neurons (seed 0): et vs global rel %+.2f%%; tag benefit rel %+.2f%% (95%% CI %+.2f%%..%+.2f%%), sign-flip p %.3f" % (
        100 * g["rel"], 100 * c["rel"], 100 * c["ci95_rel"][0], 100 * c["ci95_rel"][1], c["p_signflip"]))
    verdict("permuted outcomes: no model skill and no tag benefit", g["rel"] < 0.01 and c["rel"] < MIN_EFFECT and c["ci95"][0] <= 0)
    LABEL = "CHECK"
    json.dump({"planted": out, "all_pass": all(ok), "n_checks": len(ok), "script_sha256": sha256(path(SCRIPT_REL)),
               "script_sha256_pinned_blanked": script_sha256_blanked(), "splits_sha256": sha256(SPLITS), "versions": versions()},
              open(os.path.join(INTERIM, "d109_check_results.json"), "w"), indent=1, sort_keys=True)
    say("ALL CHECKS PASS" if all(ok) else "SOME CHECKS FAIL: %d of %d" % (len(ok) - sum(ok), len(ok)))


def mode_dry():
    global LABEL
    LABEL = "PERMUTED"
    D = load_inputs()
    S = load_splits(D)
    assert_splits(D, S)
    say("DRY RUN: outcome rows permuted across neurons with seed 1. NOT A RESULT. Tag-shuffle permutations reduced to %d / %d (run: %d / %d)." % (
        DRY_SHUFFLE, DRY_SHUFFLE_UNSTRAT, N_SHUFFLE, N_SHUFFLE_UNSTRAT))
    evaluate(D, outcomes_permuted(D.ids, 1), S, DRY_SHUFFLE, DRY_SHUFFLE_UNSTRAT, True, "dry", None)


def current_hashes():
    h = {k: sha256(path(v)) for k, v in INPUTS.items()}      # bytes only; nothing is parsed
    h["protocol"] = sha256(path(PROTOCOL_REL))
    h["splits"] = sha256(SPLITS) if os.path.exists(SPLITS) else None
    h["review"] = sha256(path(REVIEW_REL)) if os.path.exists(path(REVIEW_REL)) else None
    return h


def script_sha256_blanked():
    """SHA-256 of this file with the contents of the PINNED block removed: the hash the review records, so
    that pasting PINNED after review is the only change the run accepts."""
    txt = open(path(SCRIPT_REL), encoding="utf-8").read()
    new, n = re.subn(r"(?ms)^PINNED = \{\n.*?^\}\n", "PINNED = {\n}\n", txt, count=1)
    assert n == 1
    return hashlib.sha256(new.encode("utf-8")).hexdigest()


REVIEW_FIELDS = {"script": r"REVIEWED SCRIPT SHA-256 \(PINNED block blanked\): ([0-9a-f]{64})",
                 "protocol": r"REVIEWED PROTOCOL SHA-256: ([0-9a-f]{64})",
                 "splits": r"REVIEWED SPLITS SHA-256: ([0-9a-f]{64})"}


def read_review():
    """(verdict line present, {field: hash}) from the review file; each line must stand alone."""
    lines = [x.strip() for x in open(path(REVIEW_REL), encoding="utf-8").read().splitlines()]
    got = {}
    for k, pat in REVIEW_FIELDS.items():
        m = [re.fullmatch(pat, x) for x in lines]
        m = [x.group(1) for x in m if x]
        got[k] = m[0] if len(m) == 1 else None
    return lines.count(VERDICT_LINE) == 1, got


def mode_hashes():
    h = current_hashes()
    print("PINNED = {\n" + "".join("    %r: %r,\n" % kv for kv in h.items()) + "}")
    print("script  sha256 %s  git blob %s" % (sha256(path(SCRIPT_REL)), git("hash-object", path(SCRIPT_REL), enforce=False)))
    print("Lines the review file must contain (each on its own line), next to `%s`:" % VERDICT_LINE)
    print("REVIEWED SCRIPT SHA-256 (PINNED block blanked): %s" % script_sha256_blanked())
    print("REVIEWED PROTOCOL SHA-256: %s" % h["protocol"])
    print("REVIEWED SPLITS SHA-256: %s" % h["splits"])
    print("versions:", versions())


def git(*a, enforce=True):
    r = subprocess.run(["git", "-C", ROOT] + list(a), capture_output=True, text=True)
    if r.returncode and enforce:
        sys.exit("refusing: git %s failed: %s" % (" ".join(a), r.stderr.strip()))
    return r.stdout.strip()


def mode_run():
    global LABEL, _RUN_OPEN
    need = ["soma", "labels", "outcome", "ccf", "protocol", "splits", "review"]
    if not PINNED or sorted(PINNED) != sorted(need):
        sys.exit("refusing: PINNED is empty or incomplete. Nothing is scored until the protocol, code and splits are "
                 "reviewed, frozen (`hashes`), pinned and committed. Required keys: %s" % ", ".join(need))
    h = current_hashes()
    bad = [k for k in need if h.get(k) != PINNED[k]]
    if bad:
        sys.exit("refusing: SHA-256 differs from PINNED for: " + ", ".join(bad))
    verdict, rv = read_review()
    if not verdict:
        sys.exit("refusing: %s does not contain exactly one line `%s`" % (REVIEW_REL, VERDICT_LINE))
    now = {"script": script_sha256_blanked(), "protocol": h["protocol"], "splits": h["splits"]}
    bad = [k for k in now if rv[k] != now[k]]
    if bad:
        sys.exit("refusing: the review does not cover this exact %s (its `REVIEWED ... SHA-256` line is missing, repeated or "
                 "different; the script is compared with its PINNED block blanked)" % ", ".join(bad))
    tracked = (SCRIPT_REL, PROTOCOL_REL, REVIEW_REL, "reports/d109_jiao_conventions.md")
    for rel in tracked:
        if git("rev-parse", "HEAD:" + rel) != git("hash-object", path(rel)):
            sys.exit("refusing: %s is not the committed version" % rel)
    if git("status", "--porcelain", "--", *tracked):
        sys.exit("refusing: uncommitted changes in the D-109 files")
    v = versions()
    if v["scikit-learn"] != SKLEARN_VERSION:
        sys.exit("refusing: scikit-learn %s, the protocol declares %s" % (v["scikit-learn"], SKLEARN_VERSION))
    D = load_inputs()
    S = load_splits(D)
    assert_splits(D, S)
    chk = os.path.join(INTERIM, "d109_check_results.json")
    ck = json.load(open(chk)) if os.path.exists(chk) else {}
    if (not ck.get("all_pass") or ck.get("script_sha256_pinned_blanked") != now["script"] or ck.get("splits_sha256") != h["splits"]
            or ck.get("versions") != v):
        sys.exit("refusing: `check` has not passed with this exact script (PINNED block blanked), split file and library versions "
                 "(data/interim/d109_check_results.json)")
    rerun = os.environ.get("D109_RERUN") == "1"
    if os.path.exists(RUN_MARKER) and not rerun:
        sys.exit("refusing: the outcomes were already opened (%s):\n%s\nSet D109_RERUN=1 to repeat (crash or bug only; reported "
                 "as a deviation; written to new file names)." % (os.path.relpath(RUN_MARKER, ROOT), open(RUN_MARKER).read()))
    rec = {"opened_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "commit": git("rev-parse", "HEAD"), "hashes": h,
           "script_sha256": sha256(path(SCRIPT_REL)), "check_results_sha256": sha256(chk), "versions": v,
           "repeat": os.path.exists(RUN_MARKER), "history": []}
    if os.path.exists(RUN_MARKER):
        old = json.load(open(RUN_MARKER))
        rec["history"] = old.get("history", []) + [{k: old[k] for k in ("opened_utc", "commit")}]
    tagfile = "run" if not rec["history"] else "run_repeat%d" % len(rec["history"])
    if os.path.exists(os.path.join(INTERIM, "d109_%s_results.json" % tagfile)):
        sys.exit("refusing: d109_%s_results.json exists; results of a real run are never overwritten" % tagfile)
    json.dump(rec, open(RUN_MARKER, "w"), indent=1)
    LABEL = "REAL"
    say("outcome table opened %s (recorded in %s)%s" % (rec["opened_utc"], os.path.relpath(RUN_MARKER, ROOT), "  THIS IS A REPEAT" if rec["repeat"] else ""))
    say("run marker (copy into the committed results report):", json.dumps(rec, sort_keys=True))
    _RUN_OPEN = True
    evaluate(D, outcomes_real(D.ids), S, N_SHUFFLE, N_SHUFFLE_UNSTRAT, True, tagfile, rec)


# ------------------------------------------------------------------------------- source documentation
ZENODO = "https://zenodo.org/records/14588094/files/Code%20for%20Hypothalamic%20projectome.zip?download=1"
NATURE = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-025-01919-0/MediaObjects/"
PORTAL = "https://mouse.digital-brain.cn/projectome/2/srv/"
PYSWC = "https://raw.githubusercontent.com/txgxxx/pyswcloader/bc64d464976dd4fad88619b64d91acd5999526fd/pyswcloader/"
# local name: (source, SHA-256). "zip:" = member of the Zenodo archive, read by range request.
SOURCES = {
    "zenodo/fig_for_plot_neurons.ipynb": ("zip:Main Figs/fig_for_plot_neurons/fig_for_plot_neurons.ipynb", "f0c5b0ec3a8edc18607a79d8e65db9714e3447b1b21c57506b5538a6d16f83b8"),
    "zenodo/fig5b_5f_ExtendedDataFig8a.ipynb": ("zip:Main Figs/fig5b&5f&Extended Data Fig8a/fig5b&5f&Extended Data Fig8a.ipynb", "e112254b81a48ad0a42b62046f35d5e7a82f69ad45d83cf4f2ffea0ba3798b16"),
    "zenodo/Fig_6a.ipynb": ("zip:Main Figs/fig6a/Fig 6a.ipynb", "0108a2881ec24b7e8636eee5f266da3e6f275380607fb9dc922c36ba4019a357"),
    "zenodo/fig1b_brain_areas.xlsx": ("zip:Main Figs/fig1b/brain_areas.xlsx", "38c22fbac2df9c29136fcd0de2fa7a9198c2ccd4126df10a93b3bc1f4765da3a"),
    "zenodo/brain_area_sort.csv": ("zip:Main Figs/fig1f&Supplementary fig3/brain_area_sort.csv", "96fac9c9079d10cf73e77c00a8f090ea7fd08a9135e0a74ac4b9e642539e7fea"),
    "zenodo/fig_1b.ipynb": ("zip:Main Figs/fig1b/fig 1b.ipynb", "99f9aca8f60943e371f2e6e486bf12502c3d8e68d6c1f0c2765bfac19e35cb09"),
    "zenodo/fig1b_visualization.py": ("zip:Main Figs/fig1b/visualization.py", "77c7e74bb1d347e58e1e2550bd695c7ce0a9cf1051f4b1c43b1e19066058d6fa"),
    "zenodo/fig_1f_Supplementary_fig3.R": ("zip:Main Figs/fig1f&Supplementary fig3/fig_1f&Supplementary_fig3.R", "52be573fe531f4ab208e25d37de3c89822137ce2a7a6e5a3bd1cf6450550ad16"),
    "zenodo/Extended_Data_Fig10a.ipynb": ("zip:Extended Data Fig/Extended Data Fig10a/Extended Data Fig10a.ipynb", "28e8a8edb9fb657cd4e0afb9323fd9287d569bc6194d855a6c87d40651f15e67"),
    "zenodo/fig_5i-j.ipynb": ("zip:Main Figs/fig5i-j/fig 5i-j.ipynb", "f092342bc8ed915b5c285d424b48edb01377d1db388e85a7bc149e4885294e27"),
    "zenodo/Allen_volumne.csv": ("zip:Extended Data Fig/Extended Data Fig2b/Allen_volumne.csv", "a7801789f221d1d87c75893b546a2ac81ce017bcfe1a03923b882014c8ce3df8"),
    "pyswcloader/swc.py": (PYSWC + "swc.py", "5bbd41f890979edfd809dbd91224448fea00df5be5144d0350d10df17ee41c22"),
    "pyswcloader/projection.py": (PYSWC + "projection.py", "814107e57600324f4c11a8f3d0f9de60054dc158c4c4b8bc23951276fb740a80"),
    "pyswcloader/brain.py": (PYSWC + "brain.py", "f1d324a41549e6789845b3cc2908b7d69a6b5fd2b68a02cc63d742699323ce11"),
    "nature/MOESM2_reporting_summary.pdf": (NATURE + "41593_2025_1919_MOESM2_ESM.pdf", "2dd3d2ec841601dd6b0f6b4fd73e4230d8b702fabeb04e01707082334143012d"),
    "nature/MOESM3_supplementary_tables.xlsx": (NATURE + "41593_2025_1919_MOESM3_ESM.xlsx", "893d614b9490c83c3469b176594526716d11055db0461b804b5341ba37b4359b"),
    "portal/hy.neuron.info.json": (PORTAL + "info/mouse/hy/mouse.neuron.info.json", "fdf1fa2eff099bc755e5218dc6437e772f12fe0a934c01a4072d846d9b3f4da3"),
    "portal/pvh_oxt.neuron.info.json": (PORTAL + "info/mouse/pvh_oxt/mouse.neuron.info.json", "61203bb4056330349c260cf0668035a7a9b043c42561113fd8f8b89785d55b86"),
    "portal/pvh_oxt.sampleInfo.json": (PORTAL + "info/mouse/pvh_oxt/sampleInfo/PVH_Oxt.json", "a15fec89f5a8a85062c4e4967e53346509eddc9f15b8fe0e78e12999521562ed"),
}
UA = {"User-Agent": "research-client/0.1"}


class Ranged(io.RawIOBase):
    """A remote file read by HTTP range requests (as in analysis/fetch_jiao.py)."""
    def __init__(s, url):
        import urllib.request
        s.rq = urllib.request
        r = s.rq.urlopen(s.rq.Request(url, headers={**UA, "Range": "bytes=0-0"}))
        s.size, s.url, s.pos = int(r.headers["Content-Range"].split("/")[1]), r.geturl(), 0
    def seekable(s): return True
    def readable(s): return True
    def tell(s): return s.pos
    def seek(s, off, wh=0):
        s.pos = off if wh == 0 else s.pos + off if wh == 1 else s.size + off
        return s.pos
    def read(s, n=-1):
        if n < 0:
            n = s.size - s.pos
        if n == 0 or s.pos >= s.size:
            return b""
        end = min(s.size, s.pos + n) - 1
        b = s.rq.urlopen(s.rq.Request(s.url, headers={**UA, "Range": "bytes=%d-%d" % (s.pos, end)})).read()
        s.pos += len(b)
        return b
    def readinto(s, b):
        d = s.read(len(b))
        b[:len(d)] = d
        return len(d)


def mode_sources():
    fetch = "--fetch" in sys.argv
    out = next((a for a in sys.argv[2:] if not a.startswith("--")), DOCS)
    z, bad = None, 0
    for name, (src, want) in SOURCES.items():
        p = os.path.join(out, name)
        if fetch:
            import urllib.request
            import zipfile
            os.makedirs(os.path.dirname(p), exist_ok=True)
            if src.startswith("zip:"):
                z = z or zipfile.ZipFile(io.BufferedReader(Ranged(ZENODO), 1 << 20))
                data = z.read(src[4:])
            else:
                data = urllib.request.urlopen(urllib.request.Request(src, headers=UA)).read()
            open(p, "wb").write(data)
        got = sha256(p) if os.path.exists(p) else "missing"
        bad += got != want
        print("%-8s %s  %s\n         <- %s" % ("ok" if got == want else "DIFFERS" if got != "missing" else "MISSING", got, os.path.relpath(p, ROOT), src))
    print("%d of %d source documents match the recorded SHA-256" % (len(SOURCES) - bad, len(SOURCES)))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    fn = {"splits": mode_splits, "check": mode_check, "dry": mode_dry, "run": mode_run, "hashes": mode_hashes, "sources": mode_sources}.get(mode)
    if fn is None:
        sys.exit(__doc__)
    fn()
