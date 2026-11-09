"""D-095: does how similar two circuits are predict whether they do the same thing?

Pre-registered in docs/DECISIONS.md D-095. This file was committed before --run was
executed. --run refuses to start unless every data input matches the SHA-256 pinned
below and the code and curated data are unmodified from the commit it prints.

    python3 analysis/d095_pair_test.py --power   # counts, power, size: never the observed statistic
    python3 analysis/d095_pair_test.py --run     # the pre-registered test, run once

The model's inputs, none of which a paper chose (D-091):
  genes         data/interim/hypomap_c465_means.npz from extract_hypomap_means.py: mean
                expression per HypoMap cell type; the 2,000 most variable genes, z-scored
  connectivity  Allen projection density from fetch_connectivity.py. D-093's feature,
                unchanged: every structure outside the hypothalamus subtree, which keeps
                the hypothalamus's ancestor structures and fibre tracts that cross it
Papers enter only as the labels being predicted: each circuit's direction of effect.
"""
import argparse
import collections
import glob
import hashlib
import json
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from syn4rch.atlas.ontology import flatten, subtree  # noqa: E402
from syn4rch.feasibility.coverage import is_usable  # noqa: E402
from syn4rch.qc.keys import CONDITION_FIELDS, UNKNOWN, condition_key, lab_key, perturbation_key  # noqa: E402
from syn4rch.schema.io import read_jsonl  # noqa: E402

SEED = 20261109
N_PERM = 10000
ALPHA = 0.025                    # two-sided, per feature: Bonferroni over the two features
DIRECTIONS = ("increase", "decrease", "no_change")
POLARITIES = ("increase", "decrease", "remove")      # a sign that names one manipulation
CLEAN_BACKGROUND = (None, "not_reported", "wild_type")  # no mutation beyond the perturbation
# The seven laboratories D-087 found under two keys each, merged for the lab sensitivity.
LAB_ALIASES = {"yamanaka a.": "yamanaka a", "de lecea l.": "luis de lecea",
               "fuller pm.": "fuller pm", "sternson sm.": "scott m sternson",
               "huang zl.": "huang zl", "sakurai t.": "takeshi sakurai",
               "atasoy d.": "deniz atasoy"}
GENES_NPZ = "data/interim/hypomap_c465_means.npz"
CONN = "data/raw/atlas/connectivity"
SEPARABILITY_SHARE = 0.5         # a pass losing more than half its size is not separable

# Every data input, by SHA-256. --power prints what it read; --run stops on any difference.
PINNED = {
    "genes_npz": "5b103fecaacc8d814966db2e2f638e120c7de2d907f29e4755ee44af5155a496 (1 file)",
    "hy_injections": "031aa47061ec084737e5ef0128e3027ccf1d9f1dba86864be308352f55ae9fc0 (1 file)",
    "circuit_experiments": "b2ec3ef17fcd589a9b34c3028730ab778a63b8ee16a06161db669b6c2d3e8f8b (1 file)",
    "projections": "aca145fb21d6867db9c526b0a00840c7adbf589f3dab6115e33efab274768692 (43 files)",
    "structure_graph": "b5f0d024d1df09ee18aef15093f55240f4a70b3924ee9689add612dadff80386 (1 file)",
    "records": "1638e5ae81bda8cf53cffe4fa89eaddd7e2673e2613c02bde6a899792a93d85d (1 file)",
    "cluster_mapping": "66b343397db8f928535832b31e48260ce67b8dadfcb83566e9b123e8a6f38a33 (1 file)",
    "populations": "6c870a9203332e83edbda0cbaefe7b918288bb7d62746e892618d064abbe5a26 (1 file)",
}

# power and size: synthetic labels only
POWER_SIMS = 300
POWER_PERM = 2000
SIZE_SIMS = 2000                 # no-link labels, for each design's actual false-positive rate
NN_NULL = 1000                   # D-091's shuffle count, for the nearest-neighbour design
POWER_W = (0.25, 0.5, 1.0)       # latent model: share of the signal similarity explains
POWER_TAU = (0.0, 1.0)           # spread of per-readout dominant directions
THETA = 0.4307                   # N(0,1) thirds: increase / no_change / decrease
INDEP_DELTA = (-0.8, -0.4, 0.2, 0.4, 0.6, 0.8)  # agreement change, least to most similar pair
BASE_AGREE = 0.45
FAMILY_LINKAGE = ("average", "ward")
FAMILY_DIVISOR = (4, 3, 2)       # families = round(circuits / divisor)
FAMILY_PI = (0.5, 0.8, 1.0)      # chance a circuit does what its family does
ALIGN_MIN = 0.3                  # a partition models "function follows similarity" only if
                                 # same-family tracks similarity rank at least this well
TARGET_N = (15, 30, 60)
TARGET_READOUTS = 5
TARGET_PI = 0.8
TARGET_SIMS = 100
TARGET_PERM = 500

Row = collections.namedtuple("Row", "p cell study lab species direction background rid cond")


def val(x):
    return getattr(x, "value", x)


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n else v


def path(rel):
    return os.path.join(ROOT, rel)


def git(*args):
    return subprocess.run(["git", "-C", ROOT] + list(args), capture_output=True, text=True).stdout.strip()


# --- inputs, pinned ---------------------------------------------------------------------

def input_files(genes):
    return {
        "genes_npz": [genes],
        "hy_injections": [path(CONN + "/hy_injections.json")],
        "circuit_experiments": [path(CONN + "/circuit_experiments.json")],
        "projections": sorted(glob.glob(path(CONN + "/proj_*.json"))),
        "structure_graph": [path("data/raw/atlas/structure_graph_1.json")],
        "records": [path("data/curated/records.jsonl")],
        "cluster_mapping": [path("data/curated/population_cluster_mapping.json")],
        "populations": [path("data/vocab/populations.json")],
    }


def input_hashes(genes):
    out = {}
    for name, files in input_files(genes).items():
        h = hashlib.sha256()
        for f in files:
            with open(f, "rb") as fh:
                h.update(fh.read())
        out[name] = "%s (%d file%s)" % (h.hexdigest(), len(files), "" if len(files) == 1 else "s")
    return out


# --- the model's inputs: no paper touches these -----------------------------------------

def atlas_matrix(npz_path):
    """The blind fingerprint of every HypoMap cell type (D-091), and their names."""
    z = np.load(npz_path, allow_pickle=True)
    M, names, genes = z["means"], list(z["clusters"]), z["genes"]
    tech = np.array([str(g).startswith(("mt-", "Rpl", "Rps")) for g in genes])
    var = M.var(axis=0)
    var[tech] = -1
    F = M[:, np.argsort(var)[::-1][:2000]]
    sd = F.std(axis=0)
    sd[sd == 0] = 1
    return (F - F.mean(axis=0)) / sd, names


def gene_features(npz_path):
    """Fingerprints, cell-type membership and D-065 mapping tier per circuit."""
    F, names = atlas_matrix(npz_path)
    idx = {n: i for i, n in enumerate(names)}
    # D-065 locates a circuit in the atlas. Tier B, never a training target: it is used
    # only to find where a paper's circuit sits, so its finding can be checked.
    mp = json.load(open(path("data/curated/population_cluster_mapping.json")))["entries"]
    fp, members, tier = {}, {}, {}
    for e in (mp.values() if isinstance(mp, dict) else mp):
        leaves = [l for l in (e.get("leaf_clusters") or []) if l in idx]
        if leaves:
            fp[e["population_id"]] = unit(F[[idx[l] for l in leaves]].mean(axis=0))
            members[e["population_id"]] = set(leaves)
            tier[e["population_id"]] = e.get("tier")
    return fp, members, tier


def connectivity_features(variant="d093"):
    """D-093's distal feature; the variants exist only to measure what it includes."""
    nodes = flatten(json.load(open(path("data/raw/atlas/structure_graph_1.json"))))
    hy = {n["id"] for n in subtree(nodes, "HY")}
    by_acronym = {n["acronym"]: n for n in nodes}
    ancestors = set(by_acronym["HY"]["path"]) - {by_acronym["HY"]["id"]}
    grey = {n["id"] for n in subtree(nodes, "grey")}
    match = json.load(open(path(CONN + "/circuit_experiments.json")))

    def load(eid):
        rows = json.load(open(path(CONN + "/proj_%d.json" % eid)))
        return {r["structure_id"]: r["projection_density"] for r in rows}

    exps = {e: load(e) for v in match.values() for e in v}
    structs = sorted({s for p in exps.values() for s in p})
    keep = [s for s in structs if s not in hy]
    if variant in ("no_ancestors", "grey_only"):
        keep = [s for s in keep if s not in ancestors]
    if variant == "grey_only":
        keep = [s for s in keep if s in grey]
    fp = {}
    for p, eids in match.items():
        m = np.array([[exps[e].get(s, 0.0) for s in keep] for e in eids])
        fp[p] = unit(np.log1p(m).mean(axis=0))
    return fp, {p: set(v) for p, v in match.items()}


# --- the labels: what the papers found --------------------------------------------------

def lab_of(record):
    k = lab_key(record)
    return LAB_ALIASES.get(k, k)


def make_row(r, cell):
    return Row(r.population_id, cell, r.study_id, lab_of(r), val(r.species), val(r.direction),
               val(r.background_genotype), r.record_id, condition_key(r))


def load_records():
    """Rows on the corrected key, rows on D-091's key, and the exclusion sets.

    The corrected key is (readout, perturbation_key): class, polarity and molecular target
    together, the project's own key for what was done. Keying on valence alone is the
    pooling D-019 and D-025 forbid, and it is what D-091 and D-093 did. A projection-
    specific record is another population (D-050), with no fingerprint here.
    """
    records = list(read_jsonl(path("data/curated/records.jsonl")))
    reported, labs = collections.defaultdict(set), collections.defaultdict(set)
    for r in records:                               # every study and lab naming a circuit
        if r.population_id:
            reported[r.population_id].add(r.study_id)
            if lab_of(r):
                labs[r.population_id].add(lab_of(r))
    usable = [r for r in records if is_usable(r) and r.population_id]
    original = [make_row(r, (val(r.outcome_id), val(r.valence))) for r in usable
                if val(r.outcome_id) and val(r.valence) and val(r.direction)]
    d091 = collections.defaultdict(set)             # D-091's own exclusion set
    for x in original:
        d091[x.p].add(x.study)
    rows = []
    for r in usable:
        if getattr(r, "projection_specific", None):
            continue
        cls, pol, target = perturbation_key(r)
        if not val(r.outcome_id) or pol not in POLARITIES:
            continue
        if val(r.direction) not in DIRECTIONS + ("non_monotonic",):
            continue
        rows.append(make_row(r, (val(r.outcome_id), cls, pol, target)))
    return rows, original, reported, labs, d091, usable


def primary_rows(rows):
    """Wild-type or unstated background, the three directions D-091 scored."""
    return [x for x in rows if x.direction in DIRECTIONS and x.background in CLEAN_BACKGROUND]


def modal(dirs):
    top = collections.Counter(dirs).most_common(2)
    if not top or (len(top) > 1 and top[0][1] == top[1][1]):
        return None                                 # a tie has no defined answer
    return top[0][0]


def rank_matrix(S, pairs, n):
    R = np.full((n, n), np.nan)                     # similarity as a rank over eligible pairs
    sims = np.array([S[a, b] for a, b in pairs])
    ranks = np.empty(len(sims))
    ranks[np.argsort(sims, kind="stable")] = np.arange(1, len(sims) + 1)
    for (a, b), rk in zip(pairs, ranks):
        R[a, b] = R[b, a] = rk
    return R


def build(fp, members, rows, reported, labs=None, skip_pair=None):
    """Every comparison between two circuits on the same readout and manipulation.

    For each pair, records from a study naming both circuits anywhere in the corpus are
    dropped from both sides; with ``labs``, records from a laboratory naming both are too.
    Pairs sharing a cell type or a tracing experiment (D-066), or rejected by
    ``skip_pair``, are not eligible.
    """
    circuits = sorted({x.p for x in rows if x.p in fp})
    by = collections.defaultdict(list)
    for x in rows:
        if x.p in fp:
            by[x.p].append(x)
    pairs, units = [], []
    for a in range(len(circuits)):
        for b in range(a + 1, len(circuits)):
            p, q = circuits[a], circuits[b]
            if members[p] & members[q] or (skip_pair and skip_pair(p, q)):
                continue
            pairs.append((a, b))
            shared = reported[p] & reported[q]
            shared_labs = (labs[p] & labs[q]) if labs is not None else set()
            cp, cq = collections.defaultdict(list), collections.defaultdict(list)
            for side, xs in ((cp, by[p]), (cq, by[q])):
                for x in xs:
                    if x.study not in shared and x.lab not in shared_labs:
                        side[x.cell].append(x.direction)
            for cell in sorted(cp.keys() & cq.keys()):
                dp, dq = modal(cp[cell]), modal(cq[cell])
                if dp and dq:
                    units.append((a, b, cell, dp, dq))
    n = len(circuits)
    S = np.full((n, n), np.nan)
    for a, b in pairs:
        S[a, b] = S[b, a] = float(fp[circuits[a]] @ fp[circuits[b]])
    return circuits, pairs, units, S, rank_matrix(S, pairs, n)


# --- the statistic and its null ---------------------------------------------------------

def permuted_ranks(R, units, n_perm, seed):
    """Mantel null: relabel the circuits of the similarity matrix, keep the findings."""
    rng = np.random.default_rng(seed)
    I = np.array([u[0] for u in units])
    J = np.array([u[1] for u in units])
    P = np.array([rng.permutation(R.shape[0]) for _ in range(n_perm)])
    return R[P[:, I], P[:, J]]                      # (n_perm, n_units); NaN = ineligible pair


class Null:
    def __init__(self, X):
        self.M = ~np.isnan(X)
        self.X0 = np.where(self.M, X, 0.0)
        self.n = self.M.sum(axis=1)
        self.sx = self.X0.sum(axis=1)
        self.sxx = (self.X0 ** 2).sum(axis=1)

    def corr(self, a):
        a = a.astype(float)
        sa = self.M @ a
        saa = self.M @ (a * a)
        sxa = self.X0 @ a
        num = self.n * sxa - self.sx * sa
        den = np.sqrt((self.n * self.sxx - self.sx ** 2) * (self.n * saa - sa ** 2))
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(den > 0, num / den, np.nan)


def corr(x, a):
    a = np.asarray(a, dtype=float)
    if len(a) < 3 or a.std() == 0 or x.std() == 0:
        return float("nan")
    return float(np.corrcoef(x, a)[0, 1])


def two_sided_p(obs, null):
    if np.isnan(obs):
        return float("nan")                         # an undefined statistic cannot pass
    null = null[~np.isnan(null)]
    hi = (1 + (null >= obs).sum()) / (1 + len(null))
    lo = (1 + (null <= obs).sum()) / (1 + len(null))
    return min(1.0, 2 * min(hi, lo))


def evaluate(units, R, n_perm=N_PERM, seed=SEED):
    x = np.array([R[a, b] for a, b, _, _, _ in units])
    agree = np.array([dp == dq for _, _, _, dp, dq in units], dtype=float)
    obs = corr(x, agree)
    null = Null(permuted_ranks(R, units, n_perm, seed)).corr(agree) if units else np.array([])
    return x, agree, obs, null, two_sided_p(obs, null)


# --- power and size: synthetic labels with a known amount of signal ---------------------

def latent_factor(V):
    V = V - V.mean(axis=0)
    V = np.array([unit(v) for v in V])
    lam, U = np.linalg.eigh(V @ V.T)                # centred cosine: a correlation matrix
    return U * np.sqrt(np.clip(lam, 0, None))


def partition(V, method, divisor):
    from scipy.cluster.hierarchy import fcluster, linkage
    k = max(2, int(round(len(V) / divisor)))
    Z = linkage(V, method="ward") if method == "ward" else linkage(V, method="average", metric="cosine")
    return fcluster(Z, t=k, criterion="maxclust")


def family_labels(fam, n_cells, pi, rng):
    """Each family has one direction per readout; a circuit follows it with probability pi."""
    fam_dir = rng.integers(0, 3, size=(fam.max() + 1, n_cells))
    own = rng.random((len(fam), n_cells)) < pi
    return np.where(own, fam_dir[fam], rng.integers(0, 3, size=(len(fam), n_cells)))


def latent_labels(L, n_cells, w, tau, rng):
    n = L.shape[0]
    mu = rng.normal(0, tau, size=n_cells) if tau else np.zeros(n_cells)
    z = (mu[None, :] + np.sqrt(w) * (L @ rng.normal(size=(n, n_cells)))
         + np.sqrt(1 - w) * rng.normal(size=(n, n_cells)))
    return np.where(z > THETA, 2, np.where(z < -THETA, 0, 1))


class PairDesign:
    """The pre-registered pair test, run on synthetic labels over a world of circuits."""

    def __init__(self, circuits, units, R, world, cells, seed):
        wi = {c: i for i, c in enumerate(world)}
        ce = {c: i for i, c in enumerate(cells)}
        self.I = np.array([wi[circuits[u[0]]] for u in units])
        self.J = np.array([wi[circuits[u[1]]] for u in units])
        self.C = np.array([ce[u[2]] for u in units])
        self.x = np.array([R[a, b] for a, b, _, _, _ in units])
        self.null = Null(permuted_ranks(R, units, POWER_PERM, seed))

    def test_agree(self, agree):
        r = corr(self.x, agree)
        return r, two_sided_p(r, self.null.corr(agree))

    def test(self, lab):
        return self.test_agree((lab[self.I, self.C] == lab[self.J, self.C]).astype(float))


def passes(r, p, sign=None):
    """A pass; with ``sign``, only a pass in the modelled direction counts as power."""
    ok = (not np.isnan(p)) and p < ALPHA
    return ok and (sign is None or np.sign(r) == sign)


class NearestDesign:
    """D-091's nearest-neighbour test, as rescore_nearest runs it, on synthetic labels."""

    def __init__(self, fp, members, rows, excl, world, cells, seed=20240101):
        by_cell = collections.defaultdict(list)
        for x in rows:
            if x.p in fp:
                by_cell[(x.p, x.cell)].append(x)
        circuits = sorted({k[0] for k in by_cell})
        cands = {p: [q for q in circuits if q != p and not (members[p] & members[q])]
                 for p in circuits}
        found, avail = [], []
        for (p, cell), xs in sorted(by_cell.items()):
            if modal([x.direction for x in xs]) is None:
                continue
            ok = {q for q in cands[p]
                  if modal([x.direction for x in by_cell.get((q, cell), []) if x.study not in excl[p]])}
            found.append((p, cell))
            avail.append(ok)
        self.total = len(found)
        self.two_plus = sum(len(ok) >= 2 for ok in avail)

        def nearest(order):
            return [next((q for q in order[p] if q in ok), None) for (p, _), ok in zip(found, avail)]

        obs = nearest({p: sorted(cands[p], key=lambda q: -float(fp[p] @ fp[q])) for p in circuits})
        keep = [i for i, q in enumerate(obs) if q is not None]
        rng = np.random.default_rng(seed)
        null = [nearest({p: list(rng.permutation(cands[p])) for p in circuits}) for _ in range(NN_NULL)]
        wi = {c: i for i, c in enumerate(world)}
        ce = {c: i for i, c in enumerate(cells)}
        self.scored = len(keep)
        self.P = np.array([wi[found[i][0]] for i in keep])
        self.C = np.array([ce[found[i][1]] for i in keep])
        self.Q = np.array([wi[obs[i]] for i in keep])
        self.Qn = np.array([[wi[null[s][i]] for i in keep] for s in range(NN_NULL)])

    def p_value(self, lab):
        truth = lab[self.P, self.C]
        acc = (truth == lab[self.Q, self.C]).mean()
        accn = (truth[None, :] == lab[self.Qn, self.C[None, :]]).mean(axis=1)
        return (1 + (accn >= acc).sum()) / (1 + NN_NULL)


def power_tables(fp, members, units, pairs, R, circuits, rescore_rows, d091, seed):
    """Size and power of both designs under three models of how function follows similarity."""
    world = sorted({x.p for x in rescore_rows if x.p in fp} | set(circuits))
    cells = sorted({x.cell for x in rescore_rows if x.p in fp} | {u[2] for u in units})
    V = np.array([fp[c] for c in world])
    pair = PairDesign(circuits, units, R, world, cells, seed + 1)
    near = NearestDesign(fp, members, rescore_rows, d091, world, cells)
    rng = np.random.default_rng(seed + 2)

    def run(gen, sims, sign):
        hp = n05 = n0125 = 0
        for _ in range(sims):
            lab = gen()
            hp += passes(*pair.test(lab), sign=sign)
            pn = near.p_value(lab)
            n05 += pn < 0.05
            n0125 += pn < ALPHA / 2
        return hp / sims, n05 / sims, n0125 / sims

    L = latent_factor(V)
    size = [(tau,) + run(lambda: latent_labels(L, len(cells), 0.0, tau, rng), SIZE_SIMS, None)
            for tau in POWER_TAU]
    latent = []
    for tau in POWER_TAU:
        for w in POWER_W:
            latent.append((tau, w) + run(lambda: latent_labels(L, len(cells), w, tau, rng), POWER_SIMS, 1))
    family = []
    for method in FAMILY_LINKAGE:
        for div in FAMILY_DIVISOR:
            fam = partition(V, method, div)
            sizes = sorted(collections.Counter(fam).values(), reverse=True)
            align = corr(pair.x, (fam[pair.I] == fam[pair.J]).astype(float))
            for pi in FAMILY_PI:
                family.append((method, div, sizes, align, pi)
                              + run(lambda: family_labels(fam, len(cells), pi, rng), POWER_SIMS, 1))
    indep = []
    I = np.array([u[0] for u in units])
    J = np.array([u[1] for u in units])
    allr = np.sort([R[a, b] for a, b in pairs])
    q = np.searchsorted(allr, R[I, J], side="left") / max(1, len(allr) - 1)
    for delta in INDEP_DELTA:
        prob = np.clip(BASE_AGREE + delta * (q - 0.5), 0.01, 0.99)
        hits = sum(passes(*pair.test_agree((rng.random(len(prob)) < prob).astype(float)),
                          sign=np.sign(delta)) for _ in range(POWER_SIMS))
        indep.append((delta, hits / POWER_SIMS))
    return size, latent, family, indep, near


def target_power(pairs_fn, seed, n_sims=TARGET_SIMS):
    """80% family model, every eligible pair compared on TARGET_READOUTS readouts."""
    rng = np.random.default_rng(seed)
    out = {}
    for method in FAMILY_LINKAGE:
        for div in FAMILY_DIVISOR:
            hits, aligns = 0, []
            for _ in range(n_sims):
                Vs, pairs = pairs_fn(rng)
                n = len(Vs)
                S = np.full((n, n), np.nan)
                for a, b in pairs:
                    S[a, b] = S[b, a] = float(Vs[a] @ Vs[b])
                R = rank_matrix(S, pairs, n)
                units = [(a, b, k, None, None) for a, b in pairs for k in range(TARGET_READOUTS)]
                I = np.array([u[0] for u in units])
                J = np.array([u[1] for u in units])
                K = np.array([u[2] for u in units])
                null = Null(permuted_ranks(R, units, TARGET_PERM, int(rng.integers(1 << 31))))
                fam = partition(Vs, method, div)
                aligns.append(corr(R[I, J], (fam[I] == fam[J]).astype(float)))
                lab = family_labels(fam, TARGET_READOUTS, TARGET_PI, rng)
                agree = (lab[I, K] == lab[J, K]).astype(float)
                r = corr(R[I, J], agree)
                hits += passes(r, two_sided_p(r, null.corr(agree)), sign=1)
            out[(method, div)] = (hits / n_sims, float(np.nanmean(aligns)))
    return out


# --- facts about the earlier entries' key, counted ---------------------------------------

def scored_findings(fp, members, rows, excl):
    """Findings the nearest-neighbour procedure scores: a defined truth and a neighbour
    with data. Returned as {(circuit, cell): set of record ids}. Tie status only."""
    by = collections.defaultdict(list)
    for x in rows:
        if x.p in fp:
            by[(x.p, x.cell)].append(x)
    circ = sorted({k[0] for k in by})
    out = {}
    for (p, cell), xs in by.items():
        if modal([x.direction for x in xs]) is None:
            continue
        if any(q != p and not (members[p] & members[q]) and
               modal([x.direction for x in by.get((q, cell), []) if x.study not in excl[p]])
               for q in circ):
            out[(p, cell)] = frozenset(x.rid for x in xs)
    return out


# --- the re-score ------------------------------------------------------------------------

def rescore_nearest(name, fp, members, rows, excl, null_runs=1000, seed=20240101):
    """D-091 / D-093's nearest-neighbour procedure on whichever key ``rows`` carry.

    D-091's exclusion set, candidate exclusion, nearest-with-data rule, tie skipping,
    1,000 shuffles, seed and one-sided add-one p. Printed as a correction: no verdict.
    """
    by_cell = collections.defaultdict(list)
    for x in rows:
        if x.p in fp:
            by_cell[(x.p, x.cell)].append(x)
    circuits = sorted({k[0] for k in by_cell})
    truth = {k: modal([x.direction for x in v]) for k, v in by_cell.items()}
    truth = {k: v for k, v in truth.items() if v}
    cands = {p: [q for q in circuits if q != p and not (members[p] & members[q])] for p in circuits}

    def neighbour(p, cell, order):
        for q in order:
            d = modal([x.direction for x in by_cell.get((q, cell), []) if x.study not in excl[p]])
            if d:
                return d
        return None

    def run(order_for):
        h = n = 0
        for (p, cell), t in truth.items():
            d = neighbour(p, cell, order_for(p))
            if d:
                n += 1
                h += d == t
        return h, n

    order = {p: sorted(cands[p], key=lambda q: -float(fp[p] @ fp[q])) for p in circuits}
    h, n = run(lambda p: order[p])
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(null_runs):
        sh = {p: list(rng.permutation(cands[p])) for p in circuits}
        hh, nn = run(lambda p: sh[p])
        null.append(hh / nn if nn else np.nan)
    null = np.array(null)
    acc = h / n if n else float("nan")
    p_emp = (1 + (null >= acc).sum()) / (1 + len(null))
    print("  %-34s circuits=%-2d findings=%-3d acc=%.3f  null med=%.3f 95th=%.3f  p=%.3f"
          % (name, len(circuits), n, acc, np.nanmedian(null), np.nanpercentile(null, 95), p_emp))


# --- reporting ---------------------------------------------------------------------------

def describe(name, circuits, pairs, units, rows):
    per_pair = collections.Counter((a, b) for a, b, _, _, _ in units)
    per_circuit = collections.Counter()
    for a, b, _, _, _ in units:
        per_circuit[circuits[a]] += 1
        per_circuit[circuits[b]] += 1
    rat = {(x.p, x.cell) for x in rows if x.species == "rattus_norvegicus"}
    n_rat = sum(((circuits[a], c) in rat) or ((circuits[b], c) in rat) for a, b, c, _, _ in units)
    big = sum(circuits[a] in ("lha_orexin", "arc_agrp") or circuits[b] in ("lha_orexin", "arc_agrp")
              for a, b, _, _, _ in units)
    print("%s" % name)
    print("  circuits with findings      : %d, of which %d in a comparison"
          % (len(circuits), len(per_circuit)))
    print("  pairs eligible              : %d" % len(pairs))
    print("  pairs with a comparison     : %d" % len(per_pair))
    print("  comparisons                 : %d   (%d distinct readout x manipulation cells)"
          % (len(units), len({u[2] for u in units})))
    print("  comparisons per pair        : median %d, max %d"
          % (np.median(list(per_pair.values())) if per_pair else 0, max(per_pair.values(), default=0)))
    print("  by manipulation class       : %s" % dict(collections.Counter(u[2][1] for u in units)))
    print("  involving lha_orexin/arc_agrp: %d of %d" % (big, len(units)))
    print("  with a rat record on a side : %d of %d" % (n_rat, len(units)))
    print("  circuit appearances         : %s"
          % ", ".join("%s %d" % kv for kv in per_circuit.most_common()))


def report(label, units, R, decide=False):
    x, agree, obs, null, p = evaluate(units, R)
    if not units:
        print("  %-36s no comparisons" % label)
        return obs, p
    top = x > np.median(sorted({R[a, b] for a, b, _, _, _ in units}))
    verdict = ("  -> " + ("undefined" if np.isnan(p) else "PASS" if p < ALPHA else "not significant")
               if decide else "")
    print("  %-36s n=%-3d r=%+.3f  null 2.5/50/97.5 = %+.3f/%+.3f/%+.3f  p=%.4f%s"
          % (label, len(units), obs, np.nanpercentile(null, 2.5), np.nanmedian(null),
             np.nanpercentile(null, 97.5), p, verdict))
    print("  %-36s agreement: more-similar half %.3f (%d), less-similar half %.3f (%d), overall %.3f"
          % ("", agree[top].mean() if top.any() else float("nan"), top.sum(),
             agree[~top].mean() if (~top).any() else float("nan"), (~top).sum(), agree.mean()))
    return obs, p


def pairs_table(circuits, units, S):
    rows = collections.defaultdict(lambda: [0, 0])
    for a, b, _, dp, dq in units:
        rows[(a, b)][0] += 1
        rows[(a, b)][1] += dp == dq
    for (a, b), (k, ag) in sorted(rows.items(), key=lambda kv: -S[kv[0]]):
        print("    %-15s %-15s sim %+.3f   agree %d of %d" % (circuits[a], circuits[b], S[a, b], ag, k))


def sensitivity_inputs(fp, mem, tier, rows, prim, reported, labs, nucleus, feature):
    """The pre-registered sensitivity rows: name -> (rows, labs, skip_pair). (d) is not here:
    it filters comparisons by their observed directions and exists only in --run."""
    same_nucleus = lambda p, q: nucleus.get(p) == nucleus.get(q)  # noqa: E731
    out = collections.OrderedDict()
    out["a  laboratory excluded"] = (prim, labs, None)
    out["b  same-nucleus pairs out"] = (prim, None, same_nucleus)
    out["c1 without lha_orexin"] = ([x for x in prim if x.p != "lha_orexin"], None, None)
    out["c2 without arc_agrp"] = ([x for x in prim if x.p != "arc_agrp"], None, None)
    out["e  rat records dropped"] = ([x for x in prim if x.species != "rattus_norvegicus"], None, None)
    out["f  mutant and mixed backgrounds in"] = ([x for x in rows if x.direction in DIRECTIONS], None, None)
    if feature == "genes":
        tier_c = {p for p, t in tier.items() if t == "C"}
        out["g  tier-C broad markers out"] = ([x for x in prim if x.p not in tier_c], None, None)
    return out


#: The rows whose loss of sign or size makes a pass "not separable" (D-095).
SEPARABILITY_ROWS = ("a  laboratory excluded", "b  same-nucleus pairs out", "c1 without lha_orexin",
                     "c2 without arc_agrp", "g  tier-C broad markers out")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default=path(GENES_NPZ), help="per-cell-type means (extract_hypomap_means.py)")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--power", action="store_true")
    mode.add_argument("--run", action="store_true")
    args = ap.parse_args()

    hashes = input_hashes(args.genes)
    differs = {k: v for k, v in hashes.items() if PINNED.get(k) != v}
    if args.run:
        if differs:
            sys.exit("refusing to run: inputs differ from the pinned ones: %s" % ", ".join(sorted(differs)))
        dirty = git("status", "--porcelain", "--", "analysis", "src", "data/curated", "data/vocab")
        if dirty:
            sys.exit("refusing to run: uncommitted changes to code or curated data:\n" + dirty)
        print("commit %s; every input matches its pinned SHA-256" % git("rev-parse", "HEAD"))
        print()

    rows, original, reported, labs, d091, usable = load_records()
    prim = primary_rows(rows)
    fg, mg, tier = gene_features(args.genes)
    fc, mc = connectivity_features()
    features = {"genes": (fg, mg), "connectivity": (fc, mc)}
    built = {k: build(fp, mem, prim, reported) for k, (fp, mem) in features.items()}
    nucleus = {p: e.get("nucleus") for p, e in
               json.load(open(path("data/vocab/populations.json")))["entries"].items()}

    if args.power:
        print("inputs (SHA-256):")
        for k, v in hashes.items():
            print("  %-20s %s%s" % (k, v, "" if PINNED.get(k) == v else "   <- not pinned"))
        print()
        mol = sum(perturbation_key(r)[0] == "molecular" for r in usable)
        term = sum(bool(getattr(r, "projection_specific", None)) for r in usable)
        print("circuit records: %d usable; %d molecular; %d projection-specific" % (len(usable), mol, term))
        print("rows: %d somatic with a signed manipulation; %d in the primary analysis"
              % (len(rows), len(prim)))
        print("  set aside, by background and circuit: %s; non_monotonic %d"
              % (", ".join("%s on %s %d" % (b, p, n) for (b, p), n in sorted(collections.Counter(
                  (x.background, x.p) for x in rows if x.background not in CLEAN_BACKGROUND).items())),
                 sum(x.direction == "non_monotonic" for x in rows)))
        known = np.array([[c != UNKNOWN for c in x.cond] for x in prim])
        print("  condition slots known in the primary rows: %s; mean %.2f of %d"
              % (", ".join("%s %d%%" % (f, round(100 * k)) for f, k in zip(CONDITION_FIELDS, known.mean(axis=0))),
                 known.sum(axis=1).mean(), len(CONDITION_FIELDS)))
        mol_cells = collections.defaultdict(set)
        for x in prim:
            if x.cell[1] == "molecular":
                mol_cells[x.cell].add(x.p)
        two = [c for c, ps in mol_cells.items() if len(ps) >= 2]
        print("  molecular cells measured in two or more circuits: %d (%s); compared after the "
              "independence rule: %d"
              % (len(two), ", ".join(sorted({c[3] for c in two})),
                 sum(u[2][1] == "molecular" for k in built for u in built[k][2])))
        for k, (fp, mem) in features.items():
            old = scored_findings(fp, mem, original, d091)
            new = scored_findings(fp, mem, rows, d091)
            unclear = sum(cell[1] == "unclear" for (_, cell) in old)
            term_f = sum(any(getattr(r, "projection_specific", None) for r in usable if r.record_id in s)
                         for s in old.values())
            same = sum(s in set(old.values()) for s in new.values())
            print("  %s: D-091/D-093 key scored %d findings, %d valence-unclear, %d holding a "
                  "projection-specific record; corrected key scores %d, of which %d have exactly the "
                  "record set of an originally scored finding"
                  % (k, len(old), unclear, term_f, len(new), same))
        print()
        fd, _ = connectivity_features("no_ancestors")
        fgr, _ = connectivity_features("grey_only")
        from scipy.stats import spearmanr
        cc = sorted(set(fc) & {x.p for x in prim})
        base = [float(fc[a] @ fc[b]) for i, a in enumerate(cc) for b in cc[i + 1:]]
        for name, fv in (("ancestors of HY removed", fd), ("grey matter only", fgr)):
            other = [float(fv[a] @ fv[b]) for i, a in enumerate(cc) for b in cc[i + 1:]]
            print("connectivity variant, %s: Spearman %.3f against D-093's feature over %d pairs"
                  % (name, spearmanr(base, other).correlation, len(base)))
        print()
        for k, (circuits, pairs, units, S, R) in built.items():
            fp, mem = features[k]
            describe(k, circuits, pairs, units, prim)
            ranked = sorted(collections.Counter((a, b) for a, b, *_ in units).items(), key=lambda kv: -R[kv[0]])
            print("  most similar compared pairs : %s"
                  % "; ".join("%s-%s rank %d/%d (%d)" % (circuits[a], circuits[b], R[a, b], len(pairs), n)
                              for (a, b), n in ranked[:8]))
            sens = sensitivity_inputs(fp, mem, tier, rows, prim, reported, labs, nucleus, k)
            print("  sensitivity comparisons     : %s"
                  % ", ".join("%s %d" % (name.split()[0], len(build(fp, mem, r_, reported, labs=l_, skip_pair=s_)[2]))
                              for name, (r_, l_, s_) in sens.items()))
            size, latent, family, indep, near = power_tables(fp, mem, units, pairs, R, circuits, rows, d091, SEED)
            print("  nearest-neighbour re-score: %d findings, %d scored, %d with 2+ candidate neighbours"
                  % (near.total, near.scored, near.two_plus))
            print("  pair test: two-sided %.3f, power counts passes in the modelled direction only. "
                  "nearest neighbour: one-sided 0.05 (D-091's rule) / %.4f" % (ALPHA, ALPHA / 2))
            print("  size -- no link, %d datasets each:" % SIZE_SIMS)
            for tau, pp, n05, n0125 in size:
                print("    tau %.1f   pair %.3f   nearest %.3f / %.3f" % (tau, pp, n05, n0125))
            print("  latent model -- similarity sets the correlation between circuits:")
            for tau, w, pp, n05, n0125 in latent:
                print("    tau %.1f w %.2f   pair %.2f   nearest %.2f / %.2f" % (tau, w, pp, n05, n0125))
            print("  family model -- a circuit does what its similarity family does, with probability pi:")
            for method, div, sizes, align, pi, pp, n05, n0125 in family:
                print("    %-7s k=n/%d sizes %-24s align %+.2f%s pi %.1f   pair %.2f   nearest %.2f / %.2f"
                      % (method, div, sizes, align, "  " if align >= ALIGN_MIN else " *", pi, pp, n05, n0125))
            print("    (* same-family tracks similarity rank below %.1f on the compared pairs: not a model "
                  "of function following similarity)" % ALIGN_MIN)
            print("  independent comparisons -- agreement changes by delta from least to most similar pair:")
            for delta, pp in indep:
                print("    delta %+.1f   pair %.2f" % (delta, pp))
            print()
        print("data target: family model pi %.1f, every eligible pair compared on the same %d readouts "
              "(power in the modelled direction; align in brackets)" % (TARGET_PI, TARGET_READOUTS))
        for k, (circuits, pairs, units, S, R) in built.items():
            fp, _ = features[k]
            V = np.array([fp[c] for c in circuits])
            res = target_power(lambda rng: (V, pairs), SEED + 10)
            print("  %s, the corpus's own %d circuits (%d comparisons): %s"
                  % (k, len(circuits), len(pairs) * TARGET_READOUTS,
                     "  ".join("%s n/%d %.2f (%+.2f)" % (m, d, v, a) for (m, d), (v, a) in res.items())))
        F, _ = atlas_matrix(args.genes)
        A = F / np.linalg.norm(F, axis=1, keepdims=True)
        for n in TARGET_N:
            def draw(rng, n=n):
                pick = rng.choice(len(A), size=n, replace=False)
                return A[pick], [(a, b) for a in range(n) for b in range(a + 1, n)]
            res = target_power(draw, SEED + 20 + n)
            print("  random HypoMap cell types, %d circuits (%d comparisons): %s"
                  % (n, n * (n - 1) // 2 * TARGET_READOUTS,
                     "  ".join("%s n/%d %.2f (%+.2f)" % (m, d, v, a) for (m, d), (v, a) in res.items())))
        gc = built["genes"][0]
        cg = [float(fg[a] @ fg[b]) for i, a in enumerate(gc) for b in gc[i + 1:]]
        rng = np.random.default_rng(SEED + 30)
        rnd = []
        for _ in range(1000):
            pick = A[rng.choice(len(A), size=len(gc), replace=False)]
            Sx = pick @ pick.T
            rnd.append(Sx[np.triu_indices(len(gc), 1)].mean())
        print("  mean pairwise gene cosine: the corpus's %d circuits %.3f; %d random cell types %.3f"
              % (len(gc), np.mean(cg), len(gc), np.mean(rnd)))
        return

    # ---- the pre-registered test: run once ------------------------------------------
    print("PRIMARY (the only decisions: two-sided p < %.3f)" % ALPHA)
    primary = {}
    for k, (circuits, pairs, units, S, R) in built.items():
        primary[k] = report(k, units, R, decide=True)
    print()

    print("SENSITIVITY (no decision; read by the rule fixed in D-095)")
    sens_r = {}
    for k, (fp, mem) in features.items():
        print(" %s" % k)
        out = {}
        for name, (r_, l_, s_) in sensitivity_inputs(fp, mem, tier, rows, prim, reported, labs,
                                                      nucleus, k).items():
            _, _, units, _, R = build(fp, mem, r_, reported, labs=l_, skip_pair=s_)
            out[name] = report(name, units, R)[0]
        circuits, pairs, units, S, R = built[k]
        report("d  both-no-change dropped", [u for u in units if not (u[3] == u[4] == "no_change")], R)
        sens_r[k] = out
    print()
    for k, (obs, p) in primary.items():
        if np.isnan(p) or p >= ALPHA:
            continue
        lost = [name for name, r in sens_r[k].items() if name in SEPARABILITY_ROWS and
                (np.isnan(r) or np.sign(r) != np.sign(obs) or abs(r) < SEPARABILITY_SHARE * abs(obs))]
        print("  %s passes (%s); %s" % (k, "similar circuits agree more" if obs > 0 else "similar circuits differ more",
                                       ("NOT separable from: " + ", ".join(lost)) if lost
                                       else "keeps its sign and at least half its size in every separability row"))
    print()

    both = sorted(set(fg) & set(fc))
    with_rows = sorted(set(both) & {x.p for x in prim})
    print("SECONDARY (no decision): both features on the %d circuits both place that have findings, "
          "pairs eligible under both" % len(with_rows))
    for k, fp in (("genes", fg), ("connectivity", fc)):
        mem = {p: {("g", x) for x in mg.get(p, ())} | {("c", x) for x in mc.get(p, ())} for p in both}
        _, _, units, _, R = build({p: fp[p] for p in both}, mem, prim, reported)
        report(k, units, R)
    print()

    print("CORRECTION (no decision): D-091 / D-093's nearest-neighbour test on the original and the "
          "corrected key")
    for label, fp, mem in (("genes (D-092)", fg, mg), ("connectivity (D-094)", fc, mc),
                           ("genes, shared circuits", {p: fg[p] for p in both}, mg),
                           ("connectivity, shared circuits", {p: fc[p] for p in both}, mc)):
        rescore_nearest(label + ", original key", fp, mem, original, d091)
        rescore_nearest(label + ", corrected key", fp, mem, rows, d091)
    print()

    for k, (circuits, pairs, units, S, R) in built.items():
        print("every compared pair, by %s similarity:" % k)
        pairs_table(circuits, units, S)
        print()


if __name__ == "__main__":
    main()
