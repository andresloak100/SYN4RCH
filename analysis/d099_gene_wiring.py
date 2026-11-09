"""D-099: can a cell type's molecular identity predict its anatomical connectivity?

    python3 analysis/d099_gene_wiring.py tables        # atlas cluster tables (needs anndata)
    python3 analysis/d099_gene_wiring.py check-d065    # placement reproduces D-065's mapping
    python3 analysis/d099_gene_wiring.py populations   # traced populations -> HypoMap cell types
    python3 analysis/d099_gene_wiring.py fetch         # projection data for every kept experiment
    python3 analysis/d099_gene_wiring.py power         # counts, size and power on synthetic wiring
    python3 analysis/d099_gene_wiring.py run           # the pre-registered test, run once

Both sides are measured and neither comes from a paper: molecular identity from HypoMap
(blind gene fingerprints, D-091), connectivity from the Allen Mouse Brain Connectivity
Atlas's cell-type-specific (Cre) tracing experiments with a hypothalamic injection.

A traced population is (the gene a Cre line drives, the HypoMap region its injection sits
in), placed in HypoMap by the project's own rule, syn4rch.atlas.population_clusters.
assign_population (D-065), with D-065's candidate definitions (`check-d065` verifies that
they reproduce the committed mapping). For a Cre population the driver gene is the only
declared marker, so name support means the driver's symbol in the cluster's name.

`power` uses the real gene fingerprints, pair structure and injection sites with synthetic
wiring only; it never relates genes to the measured wiring.
"""
import collections
import json
import os
import re
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from syn4rch.atlas.population_clusters import Candidate, assign_population  # noqa: E402

GENES_NPZ = os.path.join(ROOT, "data/interim/hypomap_c465_means.npz")
TABLES = os.path.join(ROOT, "data/interim/hypomap_cluster_tables.npz")
POPS = os.path.join(ROOT, "data/interim/d099_populations.json")
CONN = os.path.join(ROOT, "data/raw/atlas/connectivity")
LEVELS = ("C66_named", "C185_named", "C286_named", "C465_named")

# Allen injection structure -> HypoMap's Region_summarized. Only unambiguous matches;
# anything else (PH, SO, SUM, TMv, STN, PSTN, SBPV, ADP) is not placed.
REGION = {
    "ARH": "Arcuate hypothalamic nucleus",
    "VMH": "Ventromedial hypothalamic nucleus",
    "SCH": "Suprachiasmatic nucleus",
    "LHA": "Lateral hypothalamic area",
    "PMv": "(Pre)Mammillary region", "PMd": "(Pre)Mammillary region",
    "MM": "(Pre)Mammillary region", "LM": "(Pre)Mammillary region",
    "MPO": "Medial preoptic area", "MPN": "Medial preoptic area",
    "LPO": "Lateral preoptic area",
    "PVH": "Paraventricular hypothalamic nucleus", "PVHd": "Paraventricular hypothalamic nucleus",
    "AHN": "Anterior hypothalamic nucleus",
    "PVp": "Periventricular hypothalamic nucleus, posterior part",
    "ZI": "Zona incerta",
    "DMH": "Dorsomedial nucleus of the hypothalamus",
    "TU": "Tuberal nucleus",
}

# The gene a Cre line drives, where the line's name is not the gene symbol.
DRIVER_SYNONYMS = {"sf1": "Nr5a1", "vgat": "Slc32a1", "vglut2": "Slc17a6", "cart": "Cartpt"}

# Two populations in one region are not compared if either driver, measured in the other's
# placed cell types, reaches this share of its level in its own: the lines label shared cells.
CROSS_EXPRESSION = 0.25
MIDLINE_Z = 5700.0                       # CCFv3 midline, um; injections are folded onto one side
KERNEL_UM = (250.0, 500.0, 1000.0)       # location covariate: Gaussian proximity scales

ALPHA = 0.05
N_PERM = 10000
SEED = 20261110
SEPARABILITY_SHARE = 0.5


def driver(line, symbols):
    """The longest hyphen-joined prefix of a line's name that is a gene symbol."""
    parts = re.split(r"[-_ ]", line or "")
    if parts and parts[0].lower() in DRIVER_SYNONYMS:
        return DRIVER_SYNONYMS[parts[0].lower()]
    hy = (line or "").split("-")
    for k in range(len(hy), 0, -1):
        cand = "-".join(hy[:k]).split("_")[0]
        if cand.lower() in symbols:
            return symbols[cand.lower()][0]
    return None


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n else v


# --- atlas tables -----------------------------------------------------------------------

def build_tables():
    import anndata as ad
    z = np.load(GENES_NPZ, allow_pickle=True)
    leaves = [str(c) for c in z["clusters"]]
    counts = z["counts"].astype(float)
    a = ad.read_h5ad(os.path.join(ROOT, "data/raw/hypomap/hypomap.h5ad"), backed="r")
    obs = a.obs[list(LEVELS) + ["Region_summarized"]].astype(str)
    regions = sorted(set(obs["Region_summarized"]) - {"NA"})
    li = {c: i for i, c in enumerate(leaves)}
    region_counts = np.zeros((len(leaves), len(regions)))
    ri = {r: j for j, r in enumerate(regions)}
    parent = {lvl: {} for lvl in LEVELS[:-1]}
    for row in obs.itertuples(index=False):
        leaf = getattr(row, "C465_named")
        reg = getattr(row, "Region_summarized")
        if reg in ri:
            region_counts[li[leaf], ri[reg]] += 1
        for lvl in LEVELS[:-1]:
            parent[lvl][leaf] = getattr(row, lvl)
    out = {"leaves": np.array(leaves), "regions": np.array(regions),
           "region_counts": region_counts, "leaf_counts": counts}
    for lvl in LEVELS[:-1]:
        out["parent_" + lvl] = np.array([parent[lvl][c] for c in leaves])
    np.savez_compressed(TABLES, **out)
    print("tables: %d leaves, %d named regions" % (len(leaves), len(regions)))


def load_tables():
    """Both archives read into memory once: an NpzFile decompresses an array on every
    access, and the 465 x 51,216 expression matrix is read in loops."""
    tf = np.load(TABLES, allow_pickle=True)
    zf = np.load(GENES_NPZ, allow_pickle=True)
    t = {k: tf[k] for k in tf.files}
    z = {k: zf[k] for k in zf.files}
    assert [str(c) for c in t["leaves"]] == [str(c) for c in z["clusters"]]
    return t, z


def level_tables(t, z):
    """Per level: cluster names, member leaves, all cells, mean expression, region counts,
    and each cluster's expression sorted descending (for ranks)."""
    leaves = [str(c) for c in t["leaves"]]
    w = t["leaf_counts"]
    out = {}
    for lvl in LEVELS:
        parents = leaves if lvl == "C465_named" else [str(p) for p in t["parent_" + lvl]]
        names = sorted(set(parents))
        members = collections.defaultdict(list)
        for i, p in enumerate(parents):
            members[p].append(i)
        means = np.array([np.average(z["means"][members[n]], axis=0, weights=w[members[n]]) for n in names])
        out[lvl] = {"names": names, "members": members, "means": means,
                    "cells": np.array([w[members[n]].sum() for n in names]),
                    "region_counts": np.array([t["region_counts"][members[n]].sum(axis=0) for n in names]),
                    "sorted_desc": -np.sort(-means, axis=1)}
    return out


def candidates_for(gene_idx, region, lt, regions, markers):
    """D-065's candidates: region share over ALL of a cluster's cells (HypoMap's unlabelled
    cells included in the denominator), the marker's rank among all genes in the cluster,
    and whether the cluster's name carries any declared marker."""
    j = list(regions).index(region)
    pats = [re.compile(r"(^|[^A-Za-z0-9])%s([^A-Za-z0-9]|$)" % re.escape(m), re.I) for m in markers]
    cands = []
    for lvl, d in lt.items():
        for k, name in enumerate(d["names"]):
            share = 100.0 * d["region_counts"][k, j] / d["cells"][k]
            mean = float(d["means"][k, gene_idx])
            row = d["sorted_desc"][k]
            rank = int(np.searchsorted(-row, -mean, side="left")) + 1
            label = name.split(":", 1)[-1]
            cands.append(Candidate(level=lvl, cluster=name, marker_mean=mean, region_share=share,
                                   marker_rank_in_cluster=rank,
                                   name_carries_marker=any(p.search(label) for p in pats)))
    return cands


def check_d065():
    """Placing each single-region D-065 entry from its own marker must reproduce its
    committed level, clusters and contrast."""
    t, z = load_tables()
    lt = level_tables(t, z)
    regions = [str(r) for r in t["regions"]]
    gi = {str(g): i for i, g in enumerate(z["genes"])}
    mp = json.load(open(os.path.join(ROOT, "data/curated/population_cluster_mapping.json")))["entries"]
    same = diff = skipped = 0
    for e in (mp.values() if isinstance(mp, dict) else mp):
        if e.get("outcome") != "mapped" or len(e.get("atlas_regions") or []) != 1 or e.get("marker") not in gi:
            skipped += 1
            continue
        a = assign_population(e["population_id"], e["marker"],
                              candidates_for(gi[e["marker"]], e["atlas_regions"][0], lt, regions,
                                             e.get("declared_markers") or [e["marker"]]))
        ok = (a.level == e["level"] + "_named" or a.level == e["level"]) and \
            [c.split(":")[0] for c in a.clusters] == [c.split(":")[0] for c in e["clusters"]] and \
            round(a.contrast or 0, 2) == round(e["contrast"] or 0, 2)
        same += ok
        diff += not ok
        if not ok:
            print("  differs: %s  committed %s %s %.2f  here %s %s %.2f" % (
                e["population_id"], e["level"], e["clusters"], e["contrast"] or 0,
                a.level, list(a.clusters), a.contrast or 0))
    print("D-065 reproduction: %d identical, %d different, %d not comparable" % (same, diff, skipped))
    return diff == 0


def place_populations():
    t, z = load_tables()
    lt = level_tables(t, z)
    regions = [str(r) for r in t["regions"]]
    genes = [str(g) for g in z["genes"]]
    gi = {g: i for i, g in enumerate(genes)}
    symbols = collections.defaultdict(list)
    for g in genes:
        symbols[g.lower()].append(g)
    rows = json.load(open(os.path.join(CONN, "hy_injections.json")))["msg"]
    pops = collections.OrderedDict()
    dropped = collections.Counter()
    for r in rows:
        line = r.get("transgenic-line")
        if not line:
            dropped["wild type (no cell type)"] += 1
            continue
        region = REGION.get(r["structure-abbrev"])
        if region is None:
            dropped["injection site not in a HypoMap region (%s)" % r["structure-abbrev"]] += 1
            continue
        g = driver(line, symbols)
        if g is None or g not in gi:
            dropped["driver gene not in HypoMap (%s)" % line] += 1
            continue
        key = "%s@%s" % (g, region)
        pops.setdefault(key, {"gene": g, "region": region, "experiments": [], "lines": set()})
        pops[key]["experiments"].append(r["id"])
        pops[key]["lines"].add(line)

    placed = {}
    for key, p in pops.items():
        a = assign_population(key, p["gene"], candidates_for(gi[p["gene"]], p["region"], lt, regions, [p["gene"]]))
        if a.refused:
            dropped["placement refused: %s" % a.outcome] += len(p["experiments"])
            continue
        leaves = sorted({lt["C465_named"]["names"][i] for c in a.clusters
                         for i in lt[a.level]["members"][c]})
        placed[key] = {"gene": p["gene"], "region": p["region"], "experiments": sorted(p["experiments"]),
                       "lines": sorted(p["lines"]), "level": a.level, "clusters": list(a.clusters),
                       "leaf_clusters": leaves, "tier": a.tier, "contrast": a.contrast}
    with open(POPS, "w") as fh:
        json.dump(placed, fh, indent=1)
    per_region = collections.Counter(p["region"] for p in placed.values())
    print("populations placed: %d (%d experiments)" % (len(placed), sum(len(p["experiments"]) for p in placed.values())))
    print("by tier: %s" % dict(collections.Counter(p["tier"] for p in placed.values())))
    print("experiments not used:")
    for k, n in sorted(dropped.items(), key=lambda kv: -kv[1]):
        print("  %3d  %s" % (n, k))
    print("populations per region:")
    for reg, n in per_region.most_common():
        print("  %2d  %s: %s" % (n, reg, ", ".join(sorted(p["gene"] for p in placed.values() if p["region"] == reg))))


def fetch():
    import fetch_connectivity as F
    placed = json.load(open(POPS))
    ids = sorted({e for p in placed.values() for e in p["experiments"]})
    got = 0
    for eid in ids:
        dest = os.path.join(CONN, "proj_%d.json" % eid)
        if os.path.exists(dest):
            got += 1
            continue
        d = F.query("model::ProjectionStructureUnionize,rma::criteria,[section_data_set_id$eq%d],"
                    "[is_injection$eqfalse],[hemisphere_id$eq3],rma::options[num_rows$eqall]" % eid)
        if not d.get("success"):
            print("FAILED %d" % eid)
            continue
        json.dump(d["msg"], open(dest + ".part", "w"))
        os.replace(dest + ".part", dest)
        got += 1
        time.sleep(1.0)
    print("projection data: %d of %d experiments" % (got, len(ids)))


# --- the two fingerprints ---------------------------------------------------------------

def wiring_keep(structs):
    """Grey matter outside the hypothalamus: no HY subtree, none of HY's ancestors (by id:
    the ontology stores a path as acronyms), no fibre tracts or ventricles."""
    from syn4rch.atlas.ontology import flatten, subtree
    nodes = flatten(json.load(open(os.path.join(ROOT, "data/raw/atlas/structure_graph_1.json"))))
    hy = {n["id"] for n in subtree(nodes, "HY")}
    by_ac = {n["acronym"]: n for n in nodes}
    anc = {by_ac[a]["id"] for a in by_ac["HY"]["path"]}
    grey = {n["id"] for n in subtree(nodes, "grey")}
    distal = [s for s in structs if s in grey and s not in hy and s not in anc]
    allt = [s for s in structs if s not in anc]
    assert not (anc & set(distal)) and not (anc & set(allt)) and anc >= {997, 8, 343, 1129}
    return distal, allt


def load_projection(eid):
    rows = json.load(open(os.path.join(CONN, "proj_%d.json" % eid)))
    return {r["structure_id"]: r["projection_density"] for r in rows}


def features(placed, variant="distal"):
    import d095_pair_test as T95
    F, names = T95.atlas_matrix(GENES_NPZ)
    idx = {n: i for i, n in enumerate(names)}
    keys = sorted(placed)
    G = np.array([unit(F[[idx[l] for l in placed[k]["leaf_clusters"]]].mean(axis=0)) for k in keys])
    exps = {e: load_projection(e) for k in keys for e in placed[k]["experiments"]}
    structs = sorted({s for d in exps.values() for s in d})
    distal, allt = wiring_keep(structs)
    keep = distal if variant == "distal" else allt
    W = np.array([unit(np.log1p(np.array([[exps[e].get(s, 0.0) for s in keep]
                                          for e in placed[k]["experiments"]])).mean(axis=0)) for k in keys])
    return keys, G, W, len(keep)


# --- pair structure: who may be compared ---------------------------------------------------

def structure(placed, keys):
    """Region of each population; the 'overlap' relation (shared placed cell types, or a
    driver expressed at >= CROSS_EXPRESSION of its own level in the partner's placed cell
    types); every same-region pair as a candidate; and the pairs eligible as labelled."""
    t, z = load_tables()
    leaves = [str(c) for c in t["leaves"]]
    li = {l: i for i, l in enumerate(leaves)}
    gi = {str(g): i for i, g in enumerate(z["genes"])}
    w = t["leaf_counts"]
    n = len(keys)
    region = np.array([placed[k]["region"] for k in keys])
    own = [[li[l] for l in placed[k]["leaf_clusters"]] for k in keys]

    def level(gene, rows):
        return float(np.average(z["means"][rows, gi[gene]], weights=w[rows]))

    overlap = np.zeros((n, n), dtype=bool)
    cross = np.zeros((n, n))
    for i in range(n):
        gi_own = level(placed[keys[i]]["gene"], own[i])
        for j in range(n):
            if i == j:
                continue
            cross[i, j] = level(placed[keys[i]]["gene"], own[j]) / gi_own if gi_own > 0 else np.inf
    for i in range(n):
        for j in range(n):
            if i != j and (set(own[i]) & set(own[j]) or cross[i, j] >= CROSS_EXPRESSION
                           or cross[j, i] >= CROSS_EXPRESSION):
                overlap[i, j] = True
    cand = [(i, j) for i in range(n) for j in range(i + 1, n) if region[i] == region[j]]
    eligible = [(i, j) for i, j in cand if not overlap[i, j]]
    by_leaf = sum(1 for i, j in cand if set(own[i]) & set(own[j]))
    return region, overlap, cand, eligible, by_leaf


def location(placed, keys, pairs):
    """Per pair: closeness of the injection sites (folded onto one hemisphere), Gaussian
    proximities at KERNEL_UM, and whether the two share an Allen injection structure."""
    rows = {r["id"]: r for r in json.load(open(os.path.join(CONN, "hy_injections.json")))["msg"]}

    def fold(c):
        return np.array([c[0], c[1], abs(c[2] - MIDLINE_Z)], dtype=float)

    pos = np.array([np.mean([fold(rows[e]["injection-coordinates"]) for e in placed[k]["experiments"]], axis=0)
                    for k in keys])
    subs = [frozenset(rows[e]["structure-abbrev"] for e in placed[k]["experiments"]) for k in keys]
    I = np.array([p[0] for p in pairs])
    J = np.array([p[1] for p in pairs])
    d = np.sqrt(((pos[I] - pos[J]) ** 2).sum(axis=1))
    cols = [-d] + [np.exp(-d ** 2 / (2 * s ** 2)) for s in KERNEL_UM] + \
        [np.array([1.0 if subs[i] & subs[j] else 0.0 for i, j in pairs])]
    return np.column_stack(cols), pos


# --- the statistics and their nulls ------------------------------------------------------

def relabellings(region, n_perm, seed, pool=None):
    """Row 0 is the identity; rows 1.. shuffle gene labels within region (and within
    `pool`, if given: populations outside it keep their own label)."""
    rng = np.random.default_rng(seed)
    n = len(region)
    key = region if pool is None else np.array([r if pool[i] else "__fixed_%d" % i for i, r in enumerate(region)])
    groups = [np.where(key == r)[0] for r in sorted(set(key))]
    P = np.empty((n_perm + 1, n), dtype=int)
    P[0] = np.arange(n)
    for k in range(1, n_perm + 1):
        perm = np.arange(n)
        for g in groups:
            if len(g) > 1:
                perm[g] = rng.permutation(g)
        P[k] = perm
    return P


class Partial:
    """Within-region correlation of gene and wiring similarity, after both are centred
    within region and residualised on the location covariates D (none: D = None).
    Exact randomization null: every same-region pair is a candidate, and under each
    relabelling a pair is scored if its relabelled populations do not overlap -- the
    eligible set moves with the gene labels, so every row scores the same number of pairs."""

    def __init__(self, Sg, overlap, cand, region, D, P):
        I = np.array([p[0] for p in cand])
        J = np.array([p[1] for p in cand])
        A, B = P[:, I], P[:, J]
        self.M = (~overlap[A, B]).astype(float)               # (K, pairs)
        regs = sorted(set(region[I]))
        self.Z = np.array([[1.0 if region[i] == r else 0.0 for r in regs] for i in I])
        self.cnt = np.maximum(self.M @ self.Z, 1.0)           # (K, R)
        X = Sg[A, B] * self.M
        Xc = self._centre(X)
        K, npair = X.shape
        if D is None:
            self.Q = np.zeros((K, npair, 0))
        else:
            Dsum = np.einsum("kp,pr,pq->krq", self.M, self.Z, D)
            Dc = (D[None, :, :] - np.einsum("krq,pr->kpq", Dsum / self.cnt[:, :, None], self.Z)) * self.M[:, :, None]
            self.Q = np.zeros((K, npair, D.shape[1]))
            for k in range(K):
                U, s, _ = np.linalg.svd(Dc[k], full_matrices=False)
                keep = s > 1e-9 * max(s.max(), 1e-300)
                self.Q[k, :, :keep.sum()] = U[:, keep]
        self.Xr = self._resid(Xc)

    def _centre(self, X):
        return (X - ((X * self.M) @ self.Z / self.cnt) @ self.Z.T) * self.M

    def _resid(self, Xc):
        return Xc - np.einsum("kpq,kq->kp", self.Q, np.einsum("kpq,kp->kq", self.Q, Xc))

    def corr(self, y):
        Yr = self._resid(self._centre(np.broadcast_to(y, self.M.shape) * self.M))
        num = (self.Xr * Yr).sum(axis=1)
        den = np.sqrt((self.Xr ** 2).sum(axis=1) * (Yr ** 2).sum(axis=1))
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(den > 0, num / den, np.nan)


class Transfer:
    """The transfer rule itself: copy the wiring of the most gene-similar non-overlapping
    population in the same region. Statistic: mean over populations with two or more such
    partners of (wiring cosine to that neighbour) - (mean wiring cosine to all partners);
    tied neighbours are averaged. Null: the same within-region relabellings."""

    def __init__(self, Sg, overlap, region, P):
        n = len(region)
        K = P.shape[0]
        self.Wt = np.zeros((K, n * n))
        for k in range(K):
            p = P[k]
            rows = []
            for i in range(n):
                partners = [j for j in range(n) if j != i and region[j] == region[i] and not overlap[p[i], p[j]]]
                if len(partners) < 2:
                    continue
                g = np.array([Sg[p[i], p[j]] for j in partners])
                best = [partners[t] for t in np.where(np.isclose(g, g.max()))[0]]
                rows.append((i, best, partners))
            for i, best, partners in rows:
                for j in best:
                    self.Wt[k, i * n + j] += 1.0 / len(best) / len(rows)
                for j in partners:
                    self.Wt[k, i * n + j] -= 1.0 / len(partners) / len(rows)

    def gain(self, Sw):
        return self.Wt @ Sw.ravel()


def one_sided_p(r):
    obs, null = r[0], r[1:]
    if np.isnan(obs):
        return float("nan")
    null = null[~np.isnan(null)]
    return (1 + (null >= obs).sum()) / (1 + len(null))


def mantel(Sx, Sy, pairs, n_perm, seed):
    """Correlation of two similarity matrices over `pairs`, null by relabelling all
    populations on the x side (no region constraint)."""
    rng = np.random.default_rng(seed)
    n = Sx.shape[0]
    I = np.array([p[0] for p in pairs])
    J = np.array([p[1] for p in pairs])
    y = Sy[I, J]
    obs = float(np.corrcoef(Sx[I, J], y)[0, 1])
    null = np.empty(n_perm)
    for k in range(n_perm):
        Pm = rng.permutation(n)
        null[k] = np.corrcoef(Sx[Pm[I], Pm[J]], y)[0, 1]
    return obs, (1 + (null >= obs).sum()) / (1 + n_perm)


# --- power and size: synthetic wiring only -------------------------------------------------

POWER_SIMS = 400
POWER_PERM = 1000
POWER_A = (0.0, 0.25, 0.5, 0.75, 1.0)
POSITION_P = (0.5, 1.0)
LATENT = 50


def unitrows(A):
    return A / np.linalg.norm(A, axis=1, keepdims=True)


def gene_factor(Sg):
    """C with C C' = Sg. C @ Z for Gaussian Z has the law of G @ R for a Gaussian random
    projection R, up to scale."""
    lam, U = np.linalg.eigh(Sg)
    return U * np.sqrt(np.clip(lam, 0, None))


def simulate(model, ctx, rng):
    """Synthetic population wiring (unit rows). Every model has a regional component and
    noise. 'gene': + a x placed gene fingerprint through a random projection. 'position':
    + a smooth random field over the real injection sites (500 um), no gene information.
    'cells': + each Cre line's labelled mixture of random cell-type wiring, no gene
    information. 'leaf': wiring generated per cell type from its own genes (strength a) and
    mixed by each line's labelled cells -- a true link, seen through the placements."""
    n, region = ctx["n"], ctx["region"]
    regs = sorted(set(region))
    U = unitrows(rng.normal(size=(len(regs), LATENT)))
    W = U[[regs.index(r) for r in region]] + unitrows(rng.normal(size=(n, LATENT)))
    kind, s = model
    if kind == "gene":
        W = W + s * unitrows(ctx["C"] @ rng.normal(size=(ctx["C"].shape[1], LATENT)))
    elif kind == "position":
        centres = ctx["allpos"][rng.integers(len(ctx["allpos"]), size=60)] + rng.normal(0, 500.0, size=(60, 3))
        V = rng.normal(size=(60, LATENT))
        field = np.array([np.mean([np.exp(-((p - centres) ** 2).sum(axis=1) / (2 * 500.0 ** 2)) @ V
                                   for p in ctx["exp_pos"][i]], axis=0) for i in range(n)])
        W = W + s * unitrows(field)
    elif kind == "cells":
        W = W + s * unitrows(ctx["Wmix"] @ rng.normal(size=(ctx["Wmix"].shape[1], LATENT)))
    elif kind == "leaf":
        leafw = s * unitrows(ctx["Fu"] @ rng.normal(size=(ctx["Fu"].shape[1], LATENT))) + \
            unitrows(rng.normal(size=(ctx["Fu"].shape[0], LATENT)))
        W = W + unitrows(ctx["Wmix"] @ leafw)
    return unitrows(W)


def context(placed, keys, G, region):
    import d095_pair_test as T95
    t, z = load_tables()
    leaves = [str(c) for c in t["leaves"]]
    li = {l: i for i, l in enumerate(leaves)}
    regions = [str(r) for r in t["regions"]]
    gi = {str(g): i for i, g in enumerate(z["genes"])}
    w = t["leaf_counts"]
    Wmix = np.zeros((len(keys), len(leaves)))
    for a, k in enumerate(keys):
        own = [li[l] for l in placed[k]["leaf_clusters"]]
        g = gi[placed[k]["gene"]]
        m_own = np.average(z["means"][own, g], weights=w[own])
        Wmix[a] = t["region_counts"][:, regions.index(placed[k]["region"])] * np.clip(z["means"][:, g] / m_own, 0, 1)
        Wmix[a] /= Wmix[a].sum()
    F, _ = T95.atlas_matrix(GENES_NPZ)
    rows = {r["id"]: r for r in json.load(open(os.path.join(CONN, "hy_injections.json")))["msg"]}
    exp_pos = [[np.array(rows[e]["injection-coordinates"], float) for e in placed[k]["experiments"]] for k in keys]
    return {"n": len(keys), "region": region, "C": gene_factor(G @ G.T), "Wmix": Wmix, "Fu": unitrows(F),
            "exp_pos": exp_pos, "allpos": np.array([p for ps in exp_pos for p in ps])}


def location_strength(y, D, reg):
    """How strongly a pairwise wiring similarity follows the location covariates within
    region: r with closeness, and the share of within-region variance they explain."""
    def centre(v):
        v = np.array(v, float).copy()
        for r in set(reg):
            m = reg == r
            v[m] -= v[m].mean(axis=0)
        return v
    yc, Dc = centre(y), centre(D)
    beta, *_ = np.linalg.lstsq(Dc, yc, rcond=None)
    return float(np.corrcoef(yc, Dc[:, 0])[0, 1]), float(1 - ((yc - Dc @ beta) ** 2).sum() / (yc ** 2).sum())


def describe_and_power():
    placed = json.load(open(POPS))
    keys, G, W_real, n_struct = features(placed)
    region, overlap, cand, eligible, by_leaf = structure(placed, keys)
    D, _ = location(placed, keys, cand)
    Sg = G @ G.T
    print("populations %d; wiring structures (grey matter outside HY, ancestors removed) %d" % (len(keys), n_struct))
    print("same-region pairs %d: %d eligible; %d excluded (%d share a placed cell type, the rest by cross-expression)"
          % (len(cand), len(eligible), len(cand) - len(eligible), by_leaf))
    for r, k in collections.Counter(region[i] for i, _ in eligible).most_common():
        print("  %3d  %s" % (k, r))
    P = relabellings(region, POWER_PERM, SEED + 1)
    prim = Partial(Sg, overlap, cand, region, D, P)
    raw = Partial(Sg, overlap, cand, region, None, P)
    tr = Transfer(Sg, overlap, region, P)
    ctx = context(placed, keys, G, region)
    I = np.array([p[0] for p in cand])
    J = np.array([p[1] for p in cand])
    rng = np.random.default_rng(SEED + 2)
    print("rates over %d synthetic datasets; primary = location-adjusted, one-sided %.2f; transfer = "
          "primary AND nearest-neighbour gain, each one-sided %.2f" % (POWER_SIMS, ALPHA, ALPHA))
    models = [("gene", a) for a in POWER_A] + [("position", s) for s in POSITION_P] + \
             [("cells", 1.0), ("position+cells", 1.0)] + [("leaf", a) for a in (1.0, 2.0, 3.0)]
    for model in models:
        hp = hr = ht = 0
        rs = []
        for _ in range(POWER_SIMS):
            if model[0] == "position+cells":
                W = unitrows(simulate(("position", 1.0), ctx, rng) + simulate(("cells", 1.0), ctx, rng))
            else:
                W = simulate(model, ctx, rng)
            Sw = W @ W.T
            y = Sw[I, J]
            rp = prim.corr(y)
            rs.append(rp[0])
            pp = one_sided_p(rp) < ALPHA
            hp += pp
            hr += one_sided_p(raw.corr(y)) < ALPHA
            ht += pp and one_sided_p(tr.gain(Sw)) < ALPHA
        print("  %-15s %4.2f   mean adjusted r %+.3f   primary %.3f   unadjusted %.3f   transfer %.3f"
              % (model[0], model[1], float(np.nanmean(rs)), hp / POWER_SIMS, hr / POWER_SIMS, ht / POWER_SIMS))

    # Position-only wiring as strong as the real wiring's own dependence on position. Uses
    # the measured wiring and the injection sites only; no gene similarity is involved.
    Ie = np.array([p[0] for p in eligible])
    Je = np.array([p[1] for p in eligible])
    De, _ = location(placed, keys, eligible)
    real_r, real_R2 = location_strength((W_real @ W_real.T)[Ie, Je], De, region[Ie])
    print("real wiring (no genes): within-region r with injection closeness %+.3f; location "
          "covariates explain %.1f%% of within-region wiring similarity" % (real_r, 100 * real_R2))
    best = None
    for s in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        R2 = np.mean([location_strength((lambda Wx: (Wx @ Wx.T)[Ie, Je])(simulate(("position", s), ctx, rng)),
                                        De, region[Ie])[1] for _ in range(100)])
        if best is None or abs(R2 - real_R2) < abs(best[1] - real_R2):
            best = (s, R2)
    hp = ht = 0
    for _ in range(POWER_SIMS):
        Sw = (lambda Wx: Wx @ Wx.T)(simulate(("position", best[0]), ctx, rng))
        pp = one_sided_p(prim.corr(Sw[I, J])) < ALPHA
        hp += pp
        ht += pp and one_sided_p(tr.gain(Sw)) < ALPHA
    print("  position matched to the real data (strength %.1f, explains %.1f%%): false pass, primary %.3f, "
          "transfer %.3f" % (best[0], 100 * best[1], hp / POWER_SIMS, ht / POWER_SIMS))


# --- the run ------------------------------------------------------------------------

def git(*args):
    import subprocess
    return subprocess.run(["git", "-C", ROOT] + list(args), capture_output=True, text=True).stdout.strip()


def input_hashes(placed):
    import hashlib
    files = {"genes_npz": [GENES_NPZ], "cluster_tables": [TABLES], "populations": [POPS],
             "hy_injections": [os.path.join(CONN, "hy_injections.json")],
             "structure_graph": [os.path.join(ROOT, "data/raw/atlas/structure_graph_1.json")],
             "projections": [os.path.join(CONN, "proj_%d.json" % e)
                             for e in sorted({e for p in placed.values() for e in p["experiments"]})]}
    out = {}
    for name, paths in files.items():
        h = hashlib.sha256()
        for f in paths:
            with open(f, "rb") as fh:
                h.update(fh.read())
        out[name] = "%s (%d file%s)" % (h.hexdigest(), len(paths), "" if len(paths) == 1 else "s")
    return out


PINNED = {
    "genes_npz": "5b103fecaacc8d814966db2e2f638e120c7de2d907f29e4755ee44af5155a496 (1 file)",
    "cluster_tables": "8e4a611df5b37abbe820ff59d03f0f561690c8202abdc8d4b6bfd3fba64b6522 (1 file)",
    "populations": "bd8ecebc70b43052d57e37b643b639de6024367a3ef622b0c6279685143636ae (1 file)",
    "hy_injections": "031aa47061ec084737e5ef0128e3027ccf1d9f1dba86864be308352f55ae9fc0 (1 file)",
    "structure_graph": "b5f0d024d1df09ee18aef15093f55240f4a70b3924ee9689add612dadff80386 (1 file)",
    "projections": "a6ef25acf06a4aa11ea7d3895d3f3692c4fda60bf8cb621b40a91252b4bb6e75 (117 files)",
}


def run():
    placed = json.load(open(POPS))
    hashes = input_hashes(placed)
    differs = sorted(k for k, v in hashes.items() if PINNED.get(k) != v)
    if differs:
        sys.exit("refusing to run: inputs differ from the pinned ones: %s" % ", ".join(differs))
    dirty = git("status", "--porcelain", "--", "analysis", "src")
    if dirty:
        sys.exit("refusing to run: uncommitted changes to the code:\n" + dirty)
    print("commit %s; every input matches its pinned SHA-256\n" % git("rev-parse", "HEAD"))

    keys, G, W, n_struct = features(placed)
    region, overlap, cand, eligible, _ = structure(placed, keys)
    D, pos = location(placed, keys, cand)
    Sg, Sw = G @ G.T, W @ W.T
    I = np.array([p[0] for p in cand])
    J = np.array([p[1] for p in cand])
    y = Sw[I, J]
    P = relabellings(region, N_PERM, SEED)

    r = Partial(Sg, overlap, cand, region, D, P).corr(y)
    p = one_sided_p(r)
    g = Transfer(Sg, overlap, region, P).gain(Sw)
    pg = one_sided_p(g)
    passed = not np.isnan(p) and p < ALPHA
    print("DECISION 1 - molecular identity predicts wiring within a region, location partialled out")
    print("  %d eligible pairs: r = %+.3f   null median %+.3f, 95th %+.3f   one-sided p = %.4f  -> %s"
          % (len(eligible), r[0], np.nanmedian(r[1:]), np.nanpercentile(r[1:], 95), p,
             "undefined" if np.isnan(p) else "PASS" if passed else "not significant"))
    print("DECISION 2 - transfer by nearest molecular neighbour beats a random same-region neighbour")
    print("  mean wiring-cosine gain %+.4f   null median %+.4f, 95th %+.4f   one-sided p = %.4f  -> %s"
          % (g[0], np.nanmedian(g[1:]), np.nanpercentile(g[1:], 95), pg,
             "transfer supported" if passed and pg < ALPHA else
             "transfer not shown" if passed else "not tested (decision 1 did not pass)"))
    print()

    print("SECONDARY (no decision)")
    all_pairs = [(i, j) for i in range(len(keys)) for j in range(i + 1, len(keys)) if not overlap[i, j]]
    o, pp = mantel(Sg, Sw, all_pairs, N_PERM, SEED + 1)
    print("  s1 genes, all pairs, location included        r = %+.3f  p = %.4f  (%d pairs)" % (o, pp, len(all_pairs)))
    Dall = -np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(axis=2))
    o, pp = mantel(Dall, Sw, all_pairs, N_PERM, SEED + 2)
    print("  s2 location only, all pairs: closeness/wiring r = %+.3f  p = %.4f" % (o, pp))
    o, pp = mantel(Dall, Sg, all_pairs, N_PERM, SEED + 3)
    print("     closeness/genes                             r = %+.3f  p = %.4f" % (o, pp))
    tier_ab = np.array([placed[k]["tier"] in ("A", "B") for k in keys])
    cand_ab = [(i, j) for i, j in cand if tier_ab[i] and tier_ab[j]]
    if cand_ab:
        Dab, _ = location(placed, keys, cand_ab)
        Pab = relabellings(region, N_PERM, SEED + 4, pool=tier_ab)
        rr = Partial(Sg, overlap, cand_ab, region, Dab, Pab).corr(Sw[[i for i, _ in cand_ab], [j for _, j in cand_ab]])
        print("  s4 tier A/B only (relabelled among A/B)      r = %+.3f  p = %.4f  (%d eligible pairs)"
              % (rr[0], one_sided_p(rr), sum(1 for i, j in cand_ab if not overlap[i, j])))
    _, _, Wa, n_all = features(placed, variant="all")
    Swa = Wa @ Wa.T
    rr = Partial(Sg, overlap, cand, region, D, P).corr(Swa[I, J])
    print("  s5 all targets incl. the hypothalamus         r = %+.3f  p = %.4f  (%d structures)"
          % (rr[0], one_sided_p(rr), n_all))
    rr = Partial(Sg, overlap, cand, region, None, P).corr(y)
    print("  s7 within region, location NOT partialled     r = %+.3f  p = %.4f" % (rr[0], one_sided_p(rr)))
    print("  s6 each region left out of decision 1:")
    s6 = {}
    for reg in sorted(set(region[I])):
        sub = [(i, j) for i, j in cand if region[i] != reg]
        Dsub, _ = location(placed, keys, sub)
        rr = Partial(Sg, overlap, sub, region, Dsub, P[:2001]).corr(Sw[[i for i, _ in sub], [j for _, j in sub]])
        s6[reg] = rr[0]
        print("     without %-52s r = %+.3f  p = %.4f" % (reg, rr[0], one_sided_p(rr)))
    if passed:
        lost = [reg for reg, v in s6.items()
                if np.isnan(v) or np.sign(v) != np.sign(r[0]) or abs(v) < SEPARABILITY_SHARE * abs(r[0])]
        print("  decision 1 passes; %s" % (("NOT separable from: " + "; ".join(lost)) if lost else
                                            "it keeps its sign and at least half its size without any one region"))
    print()
    print("per region (diagnostic): unadjusted correlation over eligible pairs")
    Ie = np.array([p[0] for p in eligible])
    Je = np.array([p[1] for p in eligible])
    for reg in sorted(set(region[Ie])):
        idx = [k for k, (i, j) in enumerate(eligible) if region[i] == reg]
        if len(idx) >= 3:
            print("  %-55s r = %+.3f  (%d pairs)" % (reg, np.corrcoef(Sg[Ie[idx], Je[idx]], Sw[Ie[idx], Je[idx]])[0, 1], len(idx)))
        else:
            print("  %-55s %d pair(s)" % (reg, len(idx)))


if __name__ == "__main__":
    cmd = sys.argv[1:2]
    if cmd == ["tables"]:
        build_tables()
    elif cmd == ["check-d065"]:
        sys.exit(0 if check_d065() else 1)
    elif cmd == ["populations"]:
        place_populations()
    elif cmd == ["fetch"]:
        fetch()
    elif cmd == ["power"]:
        describe_and_power()
    elif cmd == ["run"]:
        run()
    else:
        sys.exit(__doc__)
