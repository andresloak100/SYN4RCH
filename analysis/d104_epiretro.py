"""D-104: does a hypothalamic cell type's gene expression predict where its neurons project,
beyond where the cell type sits? D-099's question, on data where each neuron's cell type and
projection target were read from the same nucleus.

    python3 analysis/fetch_epiretro.py           # inputs (data/raw/epiretro, data/raw/abc)
    python3 analysis/d104_epiretro.py tables     # cell types, positions, fingerprints; no targets shown
    python3 analysis/d104_epiretro.py power      # size and power on synthetic projection data
    python3 analysis/d104_epiretro.py dry        # the whole run on ONE synthetic table: not a result
    python3 analysis/d104_epiretro.py run        # the pre-registered test, run once

(`tables` needs PYTHONPATH=src for the CCF ontology.)

Three measured sources, none from a paper's conclusions:

* Projections: Epi-retro-seq (Zhou et al. 2023, Supplementary Table 2). Each hypothalamic
  nucleus was labelled from one injected target and given an Allen subclass by co-clustering
  its methylome with Allen RNA data. The 1,538 hypothalamic neurons the authors kept.
* Gene expression: the Allen whole-brain 10x taxonomy's cluster means (MapMyCells reference),
  averaged into subclasses weighted by each cluster's cells. Other mice's RNA: these nuclei's
  own RNA was not measured.
* Position: the Allen MERFISH brain C57BL6J-638850, each cell registered to CCFv3 and
  assigned a cluster; a subclass's position is its MERFISH cells inside the hypothalamus.

The unit is an Allen subclass. Its projection profile is how its share of the labelled
neurons rises and falls across targets (log shares, centred over targets and over
subclasses), so a subclass's abundance and a target's labelling depth both drop out.

`tables` and `power` never relate gene similarity to the measured projections. `power` and
`run` read the measured projections without genes to calibrate the synthetic data: the
profiles' size, how far the two sexes agree, and how much of their similarity position
explains. `dry` replaces the measured counts by one synthetic table before anything is
computed from them.
"""
import json
import os
import re
import subprocess
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")          # many tiny matrices: threads only add overhead
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
import d099_gene_wiring as D  # noqa: E402  (relabellings, one_sided_p, mantel)

RAW = {
    "st2": "data/raw/epiretro/st2_cells.xlsx",
    "stats": "data/raw/abc/precomputed_stats.h5",
    "terms": "data/raw/abc/term_counts.csv",
    "genes_10x": "data/raw/abc/genes_10x.csv",
    "genes_merfish": "data/raw/abc/genes_merfish.csv",
    "merfish": "data/raw/abc/merfish_cells.csv",
    "structure_graph": "data/raw/atlas/structure_graph_1.json",
}
TABLES = os.path.join(ROOT, "data/interim/d104_tables.npz")
OUTSIDE_RULE = os.path.join(ROOT, "data/interim/d104_outside_rule.npz")   # analysis/d104_rehearsal.py freeze

MIN_NEURONS = 20                 # labelled neurons per subclass, over all 17 experiments
MIN_PURE = 10                    # the same, counting cleanly clustered neurons only (s15)
MIN_HY_SHARE = 0.5               # share of a subclass's MERFISH cells inside CCF hypothalamus
TARGETS = ("PFC", "MOB", "STR", "TH", "SC", "P", "MY")   # labelled in both sexes
SINGLE_SEX = ("AMY", "PAL", "VTA")                        # one sex only: secondary
PSEUDO = 0.5                     # added to every subclass count in an experiment
N_GENES = 2000                   # most variable genes across the compared subclasses
OVERLAP = 0.25                   # shared methylation clusters at which a pair is not compared
PURE = 0.9                       # a methylation cluster is clean if this share is one subclass
ZONES = ("PVZ", "PVR", "MEZ", "LZ")
NEAR_TIE = ("PH-ant-LHA Otp Bsx Glut", "LZ")     # 472 MERFISH cells in MEZ, 470 in LZ (s14)
MIDLINE_Z_MM = 5.7
KERNEL_UM = (250.0, 500.0, 1000.0)                        # covariates
FIELD_UM = (100.0, 150.0, 250.0, 500.0, 1000.0, 2000.0)  # random-field widths available
GRID_UM = 250.0
SLAB_MM = (5.0, 7.4)             # ARA levels 51-75: the HY dissection, one search's reading
ALPHA = 0.05
HOLM = (0.025, 0.05)             # co-registered with D-105: Holm over the two primary decisions
N_PERM = 10000
SEED = 20261113
# not fingerprint genes: mitochondrial, ribosomal-protein, and sex-linked
GENE_RE = re.compile(r"^(mt-|Rp[ls]\d+[a-z]?\d*(-ps\d+)?$|Rpl36al$|Rps27rt$|Rplp\d|Rpsa$)")
SEX_GENES = ("Xist", "Tsix", "Ddx3y", "Eif2s3y", "Uty", "Kdm5d")


def path(k):
    return os.path.join(ROOT, RAW[k])


def strip(name):
    return re.sub(r"^\d+\s+", "", str(name))


# --- tables ---------------------------------------------------------------------------------

_CACHE = {}


def load_table():
    """Every hypothalamic-source row of Supplementary Table 2 (2,542 nuclei)."""
    if "table" not in _CACHE:
        import pandas as pd
        df = pd.read_excel(path("st2"), sheet_name="TableS2", header=1)
        _CACHE["table"] = df[df["Source"] == "HY"].copy()
    return _CACHE["table"]


def load_neurons():
    """The 1,538 hypothalamic neurons the authors kept."""
    if "neurons" not in _CACHE:
        df = load_table()
        hy = df[df["PassTargetFilter"] == 1].copy()
        assert len(hy) == 1538, len(hy)
        _CACHE["neurons"] = hy
    return _CACHE["neurons"]


def zone_map():
    from syn4rch.atlas.ontology import flatten, subtree
    nodes = flatten(json.load(open(path("structure_graph"))))
    out = {}
    for z in ZONES:
        for n in subtree(nodes, z):
            out[n["acronym"]] = z
    return out


def merfish_cells():
    if "merfish" not in _CACHE:
        import pandas as pd
        cols = ["subclass", "x_ccf", "y_ccf", "z_ccf", "parcellation_division", "parcellation_structure"]
        m = pd.read_csv(path("merfish"), usecols=cols)
        m["sub"] = m["subclass"].map(strip)
        _CACHE["merfish"] = m
    return _CACHE["merfish"]


def subclass_means(names):
    """Mean log2(CPM+1) per subclass: Allen cluster means, weighted by each cluster's cells."""
    import h5py
    import pandas as pd
    f = h5py.File(path("stats"), "r")
    tree = json.loads(f["taxonomy_tree"][()].decode())
    c2r = json.loads(f["cluster_to_row"][()].decode())
    genes = json.loads(f["col_names"][()].decode())
    sub_name = {k: strip(v["name"]) for k, v in tree["name_mapper"]["CCN20230722_SUBC"].items()}
    cells = pd.read_csv(path("terms")).set_index("label")["number_of_cells"]
    sym = pd.read_csv(path("genes_10x")).set_index("gene_identifier")["gene_symbol"]
    out = np.zeros((len(names), len(genes)))
    for a, s in enumerate(names):
        sub = [k for k, v in sub_name.items() if v == s]
        assert len(sub) == 1, (s, sub)
        clusters = [c for sup in tree["CCN20230722_SUBC"][sub[0]] for c in tree["CCN20230722_SUPT"][sup]]
        rows = [c2r[c] for c in clusters]
        w = np.array([cells[c] for c in clusters], float)
        X = np.array([f["sum"][r, :] for r in rows])       # cluster means (n_cells is 1 throughout)
        out[a] = w @ X / w.sum()
    return out, np.array([sym.get(g, g) for g in genes])


def hy_cells(m, s, limit=None, seed=SEED):
    """A subclass's MERFISH cells inside the hypothalamus, folded onto one side, in um
    (columns: anterior-posterior, dorsal-ventral, distance from the midline)."""
    if "hy_xyz" not in _CACHE:                      # one pass over the MERFISH table, not one per call
        h = m[m["parcellation_division"] == "HY"]
        xyz_all = np.column_stack([h["x_ccf"], h["y_ccf"], np.abs(h["z_ccf"] - MIDLINE_Z_MM)]) * 1000.0
        subs = h["sub"].values
        _CACHE["hy_xyz"] = {k: xyz_all[subs == k] for k in np.unique(subs)}
    xyz = _CACHE["hy_xyz"].get(s, np.zeros((0, 3)))
    if limit and len(xyz) > limit:
        xyz = xyz[np.random.default_rng(seed).choice(len(xyz), limit, replace=False)]
    return xyz


def spatial_overlap(m, names, limit=1000):
    """How much two subclasses' hypothalamic cells intermingle at each width s in FIELD_UM:
    the mean Gaussian proximity exp(-d^2 / 4s^2) over cross pairs of cells. This is the
    covariance, between the two subclasses, of a smooth random field of width s averaged
    over each subclass's cells."""
    cells = [hy_cells(m, s, limit) for s in names]
    K = np.zeros((len(FIELD_UM), len(names), len(names)))
    for a in range(len(names)):
        for b in range(a, len(names)):
            d2 = ((cells[a][:, None, :] - cells[b][None, :, :]) ** 2).sum(axis=2)
            for k, sc in enumerate(FIELD_UM):
                K[k, a, b] = K[k, b, a] = np.exp(-d2 / (4 * sc ** 2)).mean()
    return K


def build_tables():
    import pandas as pd
    full = load_table()
    hy = load_neurons()
    vc = hy["Subclass"].value_counts()
    cand = sorted(s for s, n in vc.items() if n >= MIN_NEURONS)
    m = merfish_cells()
    zmap = zone_map()

    exps = sorted(hy["Experiment"].unique())
    exp_target = np.array([hy.loc[hy["Experiment"] == e, "Target"].iloc[0] for e in exps])
    exp_sex = np.array([hy.loc[hy["Experiment"] == e, "Sex"].iloc[0] for e in exps])
    sub, ex, l4col = hy["Subclass"].values, hy["Experiment"].values, hy["L4"].astype(str).values
    counts = np.array([[int(((sub == s) & (ex == e)).sum()) for e in exps] for s in cand])

    # methylation clusters (L4): sharing between subclasses, and the cleanly clustered neurons
    l4 = sorted(set(l4col))
    H = np.array([[float(((sub == s) & (l4col == c)).sum()) for c in l4] for s in cand])
    Hn = H / H.sum(axis=1, keepdims=True)
    share_l4 = np.minimum(Hn[:, None, :], Hn[None, :, :]).sum(axis=2)
    major, clean = {}, {}
    for c in l4:
        v = pd.Series(sub[l4col == c]).value_counts()
        major[c], clean[c] = v.index[0], v.iloc[0] / v.sum() >= PURE
    is_major = np.array([major[c] == s for c, s in zip(l4col, sub)])
    is_pure = is_major & np.array([clean[c] for c in l4col])
    counts_pure = np.array([[int(((sub == s) & (ex == e) & is_pure).sum()) for e in exps] for s in cand])

    hm = m[m["parcellation_division"] == "HY"]
    parcels = sorted(hm["parcellation_structure"].astype(str).unique())
    hy_share, n_hy, centroid, centroid_slab, zone, zone_counts, zone_other, parcel_frac = [], [], [], [], [], [], [], []
    grid = {}
    for a, s in enumerate(cand):
        c = m[m["sub"] == s]
        h = c[c["parcellation_division"] == "HY"]
        hy_share.append(len(h) / len(c))
        n_hy.append(len(h))
        xyz = np.column_stack([h["x_ccf"], h["y_ccf"], np.abs(h["z_ccf"] - MIDLINE_Z_MM)]) * 1000.0
        centroid.append(xyz.mean(axis=0))
        slab = (h["x_ccf"] >= SLAB_MM[0]) & (h["x_ccf"] <= SLAB_MM[1])
        centroid_slab.append(xyz[slab.values].mean(axis=0) if slab.sum() else xyz.mean(axis=0))
        zc = h["parcellation_structure"].map(zmap).value_counts()
        zone.append(zc.index[0] if len(zc) else "none")
        zone_counts.append([int(zc.get(z, 0)) for z in ZONES])
        zone_other.append(len(h) - int(zc.sum()))
        pc = h["parcellation_structure"].astype(str).value_counts()
        parcel_frac.append([pc.get(p, 0) / max(len(h), 1) for p in parcels])
        keys, cnt = np.unique(np.floor(xyz / GRID_UM).astype(int), axis=0, return_counts=True)
        grid[a] = dict(zip(map(tuple, keys), cnt))
    cooc = np.zeros((len(cand), len(cand)))
    for a in range(len(cand)):
        for b in range(len(cand)):
            ka, kb = grid[a], grid[b]
            num = sum(v * kb.get(k, 0) for k, v in ka.items())
            cooc[a, b] = num / np.sqrt(sum(v * v for v in ka.values()) * sum(v * v for v in kb.values()))
    kfield = spatial_overlap(m, cand)
    kemb = np.array([kfield[FIELD_UM.index(sc)] for sc in KERNEL_UM])
    for k in range(len(KERNEL_UM)):
        dg = np.sqrt(np.diag(kemb[k]))
        kemb[k] = kemb[k] / np.outer(dg, dg)            # 1 = identical spatial distributions

    means, symbols = subclass_means(cand)
    panel = sorted(set(pd.read_csv(path("genes_merfish"))["gene_symbol"]) & set(symbols))
    np.savez(TABLES, names=np.array(cand), n_neurons=np.array([vc[s] for s in cand]),
             hy_share=np.array(hy_share), n_hy=np.array(n_hy), centroid=np.array(centroid),
             centroid_slab=np.array(centroid_slab), zone=np.array(zone), zone_counts=np.array(zone_counts),
             zone_other=np.array(zone_other), cooc=cooc, kemb=kemb, kfield=kfield, field_um=np.array(FIELD_UM),
             parcel_names=np.array(parcels), parcel_frac=np.array(parcel_frac),
             share_l4=share_l4, means=means, symbols=symbols, panel=np.array(panel),
             exps=np.array(exps), exp_target=exp_target, exp_sex=exp_sex, counts=counts, counts_pure=counts_pure)

    keep = (np.array([vc[s] for s in cand]) >= MIN_NEURONS) & (np.array(hy_share) >= MIN_HY_SHARE)
    prim = np.isin(exp_target, TARGETS)
    print("passing hypothalamic neurons %d; subclasses with >= %d of them: %d; compared (>= %.0f%% of MERFISH "
          "cells in HY): %d, %d neurons, %d in the %d primary experiments"
          % (len(hy), MIN_NEURONS, len(cand), 100 * MIN_HY_SHARE, keep.sum(), counts[keep].sum(),
             counts[keep][:, prim].sum(), prim.sum()))
    print("experiments (target, sex, neurons of the compared subclasses, astrocyte share of all sorted nuclei):")
    for k, e in enumerate(exps):
        rows = full[full["Experiment"] == e]
        print("  %-12s %-4s %-6s %4d   %3.0f%% of %d" % (e, exp_target[k], exp_sex[k], counts[keep][:, k].sum(),
                                                       100 * (rows["L1 annot"] == "ASC").mean(), len(rows)))
    print("astrocytes among all %d sorted hypothalamic nuclei: %d (%.0f%%)"
          % (len(full), (full["L1 annot"] == "ASC").sum(), 100 * (full["L1 annot"] == "ASC").mean()))
    print("%-36s %7s %7s %5s %10s %8s %5s  %s" % ("subclass", "neurons", "primary", "clean", "MERFISH HY",
                                                  "HY share", "zone", "cells by zone PVZ/PVR/MEZ/LZ/unassigned"))
    for a, s in enumerate(cand):
        print("%-36s %7d %7d %5d %10d %7.0f%% %5s  %s%s" % (
            s, vc[s], counts[a][prim].sum(), counts_pure[a].sum(), n_hy[a], 100 * hy_share[a], zone[a],
            "/".join(str(v) for v in zone_counts[a] + [zone_other[a]]),
            "" if keep[a] else "  (not hypothalamic)"))
    inunit = np.isin(sub, np.array(cand)[keep])
    print("of the compared subclasses' neurons, %.1f%% carry a label other than their methylation cluster's "
          "majority; %d of %d clusters are clean (>= %.0f%% one subclass)"
          % (100 * (~is_major[inunit]).mean(), sum(clean.values()), len(l4), 100 * PURE))
    print("largest share of methylation clusters two compared subclasses have in common: %.2f"
          % np.triu(share_l4[np.ix_(keep, keep)], 1).max())
    print("MERFISH panel genes found in the 10x data: %d; hypothalamic parcels: %d" % (len(panel), len(parcels)))


# --- the measured pieces -------------------------------------------------------------------

def load():
    t = np.load(TABLES, allow_pickle=True)
    return {k: t[k] for k in t.files}


def select(t, min_neurons=MIN_NEURONS, min_hy=MIN_HY_SHARE):
    return np.where((t["n_neurons"] >= min_neurons) & (t["hy_share"] >= min_hy))[0]


def unitrows(A):
    n = np.linalg.norm(A, axis=1, keepdims=True)
    return A / np.where(n > 0, n, 1.0)


def double_centre(A):
    return A - A.mean(axis=1, keepdims=True) - A.mean(axis=0, keepdims=True) + A.mean()


def unit_avg(A):
    return A / np.sqrt((A ** 2).sum() / A.shape[0])


def size_of(L):
    return float(np.sqrt((L ** 2).sum() / L.shape[0]))


def fingerprints(t, idx, genes="variable"):
    """Unit rows of z-scored subclass means: the N_GENES most variable genes across the
    compared subclasses (no mitochondrial, ribosomal-protein or sex-linked genes), or the
    MERFISH panel."""
    M, sym = t["means"][idx], t["symbols"]
    if genes == "panel":
        keep = np.isin(sym, t["panel"])
    else:
        var = M.var(axis=0)
        var[np.array([bool(GENE_RE.match(str(g))) or str(g) in SEX_GENES for g in sym])] = -1
        keep = np.zeros(len(sym), bool)
        keep[np.argsort(var)[::-1][:N_GENES]] = True
    F = M[:, keep]
    sd = F.std(axis=0)
    sd[sd == 0] = 1
    return unitrows((F - F.mean(axis=0)) / sd)


def is_male(sex):
    return np.array([str(x).lower().startswith("m") for x in sex])


def log_shares(counts, exp_target, exp_sex, spec):
    """Per subclass and target: the subclass's share of each experiment's labelled neurons
    (PSEUDO added to every count), averaged over the target's experiments (one per sex, so
    the sexes weigh equally), then logged.

    spec['sex_adjust']: a target labelled in one sex only is moved half the subclass's sex
    difference towards the other sex. The difference is the mean, over the targets labelled
    in both sexes, of (log share in the male experiment - log share in the female one); the
    single-sex column is put on the scale of the two-sex average of shares."""
    S = counts.shape[0]
    q = (counts + PSEUDO) / (counts.sum(axis=0) + PSEUDO * S)
    L = np.log(np.column_stack([q[:, exp_target == tg].mean(axis=1) for tg in spec["targets"]]))
    if spec.get("sex_adjust"):
        male = is_male(exp_sex)
        both = [tg for tg in spec["targets"] if len(set(male[exp_target == tg])) == 2]
        d = np.mean([np.log(q[:, (exp_target == tg) & male].mean(axis=1))
                     - np.log(q[:, (exp_target == tg) & ~male].mean(axis=1)) for tg in both], axis=0)
        for k, tg in enumerate(spec["targets"]):
            if tg not in both:
                # the both-sex columns average shares, so move to that scale: log((1 + e^-+d) / 2)
                L[:, k] += np.log((1 + np.exp(-d if male[exp_target == tg][0] else d)) / 2)
    return L


def profiles(counts, exp_target, exp_sex, spec):
    """Projection profile: log shares centred over targets and over subclasses, unit rows."""
    L = double_centre(log_shares(counts, exp_target, exp_sex, spec))
    return unitrows(L), L


def wcorr(x, y, w):
    w = w / w.sum()
    x = x - (w * x).sum()
    y = y - (w * y).sum()
    den = np.sqrt((w * x * x).sum() * (w * y * y).sum())
    return float((w * x * y).sum() / den) if den > 0 else float("nan")


def sex_agreement(counts, exp_target, exp_sex, targets):
    """How far the two sexes' experiments agree, on the targets labelled in both: the
    correlation between the female-only and the male-only centred log shares, and between the
    two sexes' profile similarities over all pairs of subclasses. None if a sex is missing."""
    male = is_male(exp_sex)
    both = tuple(tg for tg in targets if len(set(male[exp_target == tg])) == 2)
    if len(both) < 2:
        return None
    out = []
    for mm in (False, True):
        use = (male == mm) & np.isin(exp_target, both)
        out.append(profiles(counts[:, use], exp_target[use], exp_sex[use], {"targets": both}))
    (Wf, Lf), (Wm, Lm) = out
    iu = np.triu_indices(counts.shape[0], 1)
    return float(np.corrcoef(Lf.ravel(), Lm.ravel())[0, 1]), float(np.corrcoef((Wf @ Wf.T)[iu], (Wm @ Wm.T)[iu])[0, 1])


def transmitter(name):
    return str(name).split()[-1].split("-")[-1]          # 'Glut' or 'Gaba'


def design(t, idx, n, centroid="centroid", zone=None, same_zone=False, drop_zone=None, drop_unit=None,
           same_transmitter=False):
    """Pairs compared, their 14 covariates, block and weight, the zones, and transfer partners.

    Covariates: centroid distance and its Gaussian proximities (KERNEL_UM); the intermingling
    of the two subclasses' cells at each of those scales; the share of methylation clusters
    they have in common (the route by which mislabelled cells would make a gene-similar pair
    look alike in projection); and the six products of the two subclasses' standardised
    coordinates, c_i[a] c_j[b] + c_i[b] c_j[a], which are what a linear gradient of projection
    across the hypothalamus adds to a pair's similarity.

    n: labelled neurons per subclass in the experiments used; weight 1 / (1/n_i + 1/n_j)."""
    S = len(idx)
    zone = np.array(t["zone"][idx] if zone is None else zone)
    share = t["share_l4"][np.ix_(idx, idx)]
    overlap = share >= OVERLAP
    pairs = [(i, j) for i in range(S) for j in range(i + 1, S) if not overlap[i, j]
             and (not same_zone or zone[i] == zone[j])
             and (drop_zone is None or (zone[i] != drop_zone and zone[j] != drop_zone))
             and (drop_unit is None or drop_unit not in (i, j))]
    I = np.array([p[0] for p in pairs])
    J = np.array([p[1] for p in pairs])
    pos = t[centroid][idx]
    d = np.sqrt(((pos[I] - pos[J]) ** 2).sum(axis=1))
    cols = [-d] + [np.exp(-d ** 2 / (2 * s ** 2)) for s in KERNEL_UM]
    cols += [t["kemb"][k][np.ix_(idx, idx)][I, J] for k in range(len(KERNEL_UM))] + [share[I, J]]
    c = pos - pos.mean(axis=0)
    c = c / c.std(axis=0)
    cols += [c[I, a] * c[J, b] + c[I, b] * c[J, a] for a in range(3) for b in range(a, 3)]
    if same_transmitter:
        tx = np.array([transmitter(s) for s in t["names"][idx]])
        cols += [(tx[I] == tx[J]).astype(float)]
    nn = np.maximum(np.asarray(n, float), 1.0)
    w = 1.0 / (1.0 / nn[I] + 1.0 / nn[J])
    block = np.array(["%s|%s" % tuple(sorted((zone[i], zone[j]))) for i, j in pairs])
    partners = [[j for j in range(S) if j != i and zone[j] == zone[i] and not overlap[i, j]] for i in range(S)]
    return {"n": S, "zone": zone, "overlap": overlap, "pairs": pairs, "I": I, "J": J, "D": np.column_stack(cols),
            "w": w / w.mean(), "block": block, "partners": partners, "pos": pos, "std_pos": c}


# --- the statistics --------------------------------------------------------------------------

class Partial:
    """Weighted correlation of gene and projection similarity over a fixed set of pairs,
    after both are centred within block (the pair's two zones) and residualised on the
    covariates (none: D = None), all by weighted least squares. Rows of P relabel the gene
    side (row 0 the identity). NaN where either side has nothing left after residualising."""

    def __init__(self, Sg, des, P, D="use"):
        I, J = des["I"], des["J"]
        self.sw = np.sqrt(des["w"])
        blocks = sorted(set(des["block"]))
        Z = np.array([[1.0 if b == c else 0.0 for c in blocks] for b in des["block"]])
        if D == "use":
            Z = np.column_stack([Z, des["D"]])
        U, s, _ = np.linalg.svd(Z * self.sw[:, None], full_matrices=False)
        self.Q = U[:, s > 1e-9 * s.max()]
        X = Sg[P[:, I], P[:, J]] * self.sw
        self.Xr = X - (X @ self.Q) @ self.Q.T
        self.xx = (self.Xr ** 2).sum(axis=1)
        self.x_ok = self.xx > 1e-10 * (X ** 2).sum(axis=1)

    def corr(self, y):
        y = y * self.sw
        yr = y - self.Q @ (self.Q.T @ y)
        yy = (yr ** 2).sum()
        ok = self.x_ok & (yy > 1e-10 * (y ** 2).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(ok, (self.Xr @ yr) / np.sqrt(self.xx * yy), np.nan)


class Transfer:
    """The transfer rule: copy the profile of the most gene-similar partner (same zone, not
    sharing methylation clusters). Statistic: mean over subclasses with two or more partners
    of (profile cosine to that neighbour) - (mean cosine to all its partners); ties averaged."""

    def __init__(self, Sg, des, P):
        n = des["n"]
        self.Wt = np.zeros((P.shape[0], n * n))
        rows = [i for i in range(n) if len(des["partners"][i]) >= 2]
        for k in range(P.shape[0]):
            p = P[k]
            for i in rows:
                part = des["partners"][i]
                g = np.array([Sg[p[i], p[j]] for j in part])
                best = [part[x] for x in np.where(np.isclose(g, g.max()))[0]]
                for j in best:
                    self.Wt[k, i * n + j] += 1.0 / len(best) / len(rows)
                for j in part:
                    self.Wt[k, i * n + j] -= 1.0 / len(part) / len(rows)
        self.n_rows = len(rows)

    def gain(self, Sw):
        return self.Wt @ Sw.ravel()


def location_r2(y, des):
    """Share of the block-centred projection similarity the covariates explain (weighted)."""
    sw = np.sqrt(des["w"])
    blocks = sorted(set(des["block"]))
    B = np.array([[1.0 if b == c else 0.0 for c in blocks] for b in des["block"]]) * sw[:, None]
    y = y * sw
    rb = y - B @ np.linalg.lstsq(B, y, rcond=None)[0]
    Z = np.column_stack([B, des["D"] * sw[:, None]])
    rz = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    tot = (rb ** 2).sum()
    return float(1 - (rz ** 2).sum() / tot) if tot > 0 else float("nan")


# --- one variant of the comparison: its data, design and simulation context ---------------------

FIELD_CENTRES = 60
CELLS_PER_SUBCLASS = 300


def factor(K):
    lam, U = np.linalg.eigh(K)
    return U * np.sqrt(np.clip(lam, 0, None))


def case(t, idx, targets=TARGETS, sex=None, sex_adjust=False, genes="variable", centroid="centroid",
         counts="counts", **design_kw):
    """Everything one variant needs. The experiments used are those labelled from `targets`
    (of one sex, if given). A synthetic dataset is built from marginals of the real table
    (each subclass's labelled total, each experiment's total), MERFISH positions and the gene
    similarity, and goes through the same profile code as the measured counts."""
    use = np.isin(t["exp_target"], targets)
    if sex is not None:
        use &= is_male(t["exp_sex"]) == (sex == "male")
    C = t[counts][idx][:, use].astype(float)
    et, es = t["exp_target"][use], t["exp_sex"][use]
    spec = {"targets": tuple(tg for tg in targets if (et == tg).any()), "sex_adjust": sex_adjust}
    n = C.sum(axis=1)
    G = fingerprints(t, idx, genes)
    Sg = G @ G.T
    des = design(t, idx, n, centroid=centroid, **design_kw)
    W, L = profiles(C, et, es, spec)

    m = merfish_cells()
    rng = np.random.default_rng(SEED + 7)
    cells = [hy_cells(m, s, CELLS_PER_SUBCLASS, seed=int(rng.integers(1 << 30))) for s in t["names"][idx]]
    hy = load_neurons()
    names = list(t["names"][idx])
    sub, l4col = hy["Subclass"].values, hy["L4"].astype(str).values
    l4 = sorted(set(l4col[np.isin(sub, names)]))
    A = np.array([[float(((sub == s) & (l4col == c)).sum()) for c in l4] for s in names])
    mix = (A / A.sum(axis=1, keepdims=True)) @ (A / A.sum(axis=0, keepdims=True)).T   # via shared clusters
    ctx = {"n_s": np.maximum(n, PSEUDO), "N_e": C.sum(axis=0).astype(int), "exp_target": et, "exp_sex": es,
           "spec": spec, "cells": cells, "allcells": np.vstack(cells), "C": factor(Sg), "Sg": Sg,
           "gp": {sc: factor(t["kfield"][k][np.ix_(idx, idx)]) for k, sc in enumerate(FIELD_UM)},
           "parcel": t["parcel_frac"][idx], "lin": des["std_pos"], "mix": {"l4": mix_rows(mix)}}
    return {"idx": idx, "names": t["names"][idx], "counts": C, "Sg": Sg, "des": des, "ctx": ctx, "spec": spec,
            "W": W, "L": L, "Sw": W @ W.T, "n": n}


def mix_rows(M):
    """Where a mislabelled cell goes: never back to its own subclass; rows sum to 1 (or 0)."""
    M = np.array(M, float)
    np.fill_diagonal(M, 0.0)
    s = M.sum(axis=1, keepdims=True)
    return M / np.where(s > 0, s, 1.0)


def mix_of(ctx, mix, temp):
    key = mix if mix == "l4" else "gene%.3f" % temp
    if key not in ctx["mix"]:
        ctx["mix"][key] = mix_rows(np.exp(ctx["Sg"] / temp))       # towards gene-similar subclasses
    return ctx["mix"][key]


def component(kind, ctx, rng, T, scale=500.0, axes=None):
    S = len(ctx["n_s"])
    if kind == "noise":
        A = rng.normal(size=(S, T))
    elif kind == "gene":
        A = ctx["C"] @ rng.normal(size=(ctx["C"].shape[1], T))
    elif kind == "position":
        cen = ctx["allcells"][rng.integers(len(ctx["allcells"]), size=FIELD_CENTRES)] + \
            rng.normal(0, scale, size=(FIELD_CENTRES, 3))
        V = rng.normal(size=(FIELD_CENTRES, T))
        A = np.array([(np.exp(-((c[:, None, :] - cen[None, :, :]) ** 2).sum(axis=2) / (2 * scale ** 2)) @ V).mean(axis=0)
                      for c in ctx["cells"]])
    elif kind == "gp":
        F = ctx["gp"][scale]
        A = F @ rng.normal(size=(F.shape[1], T))
    elif kind == "parcel":
        A = ctx["parcel"] @ rng.normal(size=(ctx["parcel"].shape[1], T))
    elif kind == "linear":
        X = ctx["lin"] if axes is None else ctx["lin"][:, list(axes)]
        A = X @ rng.normal(size=(X.shape[1], T))
    else:
        raise ValueError(kind)
    return unit_avg(double_centre(A))


def simulate(ctx, rng, sigma, tau=0.0, pos=0.0, field="position", scale=500.0, axes=None, pos2=0.0,
             field2="gp", scale2=500.0, gene=0.0, confusion=0.0, mix="l4", temp=0.1):
    """Synthetic counts with the real experiment totals. An experiment labelled from target t
    samples subclass s in proportion to (its labelled total) x exp(sigma * Lambda[s, t] +
    tau * e), e a fresh standard normal per experiment and subclass (what makes two experiments
    on one target disagree by more than counting noise). Lambda = noise + pos x (a position
    field) + gene x (a random projection of the gene similarity), rescaled to unit average row
    norm. Position fields, each averaged over a subclass's MERFISH cells: 'position' a sum of 60
    random bumps of width `scale` placed near the cells; 'gp' a stationary Gaussian field of
    width `scale`; 'parcel' a value per CCF structure; 'linear' a gradient (along `axes`:
    0 anterior-posterior, 1 dorsal-ventral, 2 medial-lateral; all three if None). `confusion`
    moves that share of each subclass's cells to other subclasses: along shared methylation
    clusters (mix 'l4') or in proportion to exp(gene cosine / temp) (mix 'gene')."""
    T = len(ctx["spec"]["targets"])
    S = len(ctx["n_s"])
    L = component("noise", ctx, rng, T)
    if pos:
        L = L + pos * component(field, ctx, rng, T, scale, axes)
    if pos2:
        L = L + pos2 * component(field2, ctx, rng, T, scale2)
    if gene:
        L = L + gene * component("gene", ctx, rng, T)
    L = sigma * unit_avg(L)
    M = mix_of(ctx, mix, temp) if confusion else None
    ti = {tg: k for k, tg in enumerate(ctx["spec"]["targets"])}
    out = np.zeros((S, len(ctx["exp_target"])))
    for e, tg in enumerate(ctx["exp_target"]):
        lam = L[:, ti[tg]]
        if tau:
            lam = lam + tau * rng.normal(size=S)
        p = ctx["n_s"] * np.exp(lam - lam.max())
        n = rng.multinomial(ctx["N_e"][e], p / p.sum())
        if confusion:
            moved = rng.binomial(n, confusion)
            for s in np.where(moved > 0)[0]:
                if M[s].sum() > 0:
                    n[s] -= moved[s]
                    n = n + rng.multinomial(moved[s], M[s])
        out[:, e] = n
    return out


def sim_profiles(ctx, rng, sigma, **kw):
    return profiles(simulate(ctx, rng, sigma, **kw), ctx["exp_target"], ctx["exp_sex"], ctx["spec"])


def sim_similarity(ctx, rng, sigma, **kw):
    W, L = sim_profiles(ctx, rng, sigma, **kw)
    return W @ W.T, L


# --- the null: spatial surrogates ----------------------------------------------------------
#
# Relabelling gene labels (D-099's null) does not hold its size here: when projections follow
# position more finely than the covariates capture, the residual gene and projection
# similarities share spatial structure that relabelling destroys (the partial-Mantel
# problem; `power` shows its false-pass rate). The null is instead built from surrogate
# projection data that follow position only, as strongly as the data under test do, sampled
# with the real experiment totals; the gene similarity is the real one throughout.

FAMILY = (("gp", 100.0), ("gp", 250.0), ("gp", 500.0), ("gp", 1000.0), ("parcel", None), ("linear", None))
STRENGTHS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0)
SIGMAS = (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0)
TAUS = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5)
CURVE_REPS = 100
N_SURR = 2000            # surrogates per family member in the run
N_SURR_SECONDARY = 500
CONF_RATES = (0.0, 0.05, 0.10, 0.20, 0.30)
CONF_REPS = 1000
POWER_OUTER = 1000       # synthetic datasets per model in `power`
POWER_NULL = 2000        # surrogates per family member and grid strength, shared by them
POWER_PERM = 1000        # relabellings for the D-099-style comparison
QUICK = bool(os.environ.get("D104_QUICK"))   # a fast smoke test; `run` refuses it
if QUICK:
    POWER_OUTER, POWER_NULL, CURVE_REPS, N_SURR, N_SURR_SECONDARY, CONF_REPS, N_PERM = 20, 70, 5, 40, 20, 20, 200


def member_kw(mem):
    return {"field": mem[0]} if mem[1] is None else {"field": mem[0], "scale": mem[1]}


def member_name(mem):
    return {"parcel": "per-structure field", "linear": "linear gradient"}.get(mem[0], "smooth field %d um" % (mem[1] or 0))


def interp_up(y, curve, xs):
    """x at which a rising curve reaches y (clipped to the grid's ends)."""
    c = np.maximum.accumulate(np.asarray(curve, float)) + np.arange(len(curve)) * 1e-12
    return float(np.interp(y, c, xs))


def fit_sigma(ctx, rng, size, tau):
    """Lambda's scale at which synthetic profiles are as large as the measured ones."""
    sizes = [np.mean([size_of(sim_profiles(ctx, rng, s, tau=tau)[1]) for _ in range(CURVE_REPS)]) for s in SIGMAS]
    return interp_up(size, sizes, SIGMAS)


def fit_noise(ctx, rng, size, agree):
    """Experiment-to-experiment noise tau and scale sigma, fitted together: for each tau on the
    grid, sigma is set so synthetic profiles are as large as the measured ones; tau is then
    where the synthetic sexes agree as well as the measured ones do (0 if the measured sexes
    agree better than counting noise alone allows; the grid's top if worse than it covers)."""
    rows = []
    for tau in TAUS:
        s = fit_sigma(ctx, rng, size, tau)
        a = np.mean([sex_agreement(simulate(ctx, rng, s, tau=tau), ctx["exp_target"], ctx["exp_sex"],
                                   ctx["spec"]["targets"])[0] for _ in range(CURVE_REPS)])
        rows.append((tau, s, float(a)))
    tau = interp_up(-agree, [-r[2] for r in rows], TAUS)
    return tau, fit_sigma(ctx, rng, size, tau), rows


def r2_curves(ctx, rng, sigma, tau, des):
    """Per family member: how much of the block-centred similarity the covariates explain, on
    average, at each position strength."""
    return {mem: [np.mean([location_r2(sim_similarity(ctx, rng, sigma, tau=tau, pos=b, **member_kw(mem))[0]
                                       [des["I"], des["J"]], des)
                           for _ in range(CURVE_REPS)]) for b in STRENGTHS] for mem in FAMILY}


def match_pos(ctx, rng, sigma, tau, des, r2, **kw):
    """Strength of a 500 um smooth field at which the covariates explain r2, other terms as given."""
    curve = [np.mean([location_r2(sim_similarity(ctx, rng, sigma, tau=tau, pos=b, field="gp", scale=500.0, **kw)[0]
                                  [des["I"], des["J"]], des) for _ in range(CURVE_REPS)]) for b in STRENGTHS]
    return interp_up(r2, curve, STRENGTHS)


def calibrate(c, seed, tau=None):
    """Fit the synthetic data to one variant's measured profiles, read without genes."""
    rng = np.random.default_rng(seed)
    ctx = c["ctx"]
    c["size"] = size_of(c["L"])
    c["agree"] = sex_agreement(c["counts"], ctx["exp_target"], ctx["exp_sex"], c["spec"]["targets"])
    if tau is None:
        c["tau"], c["sigma"], c["noise_rows"] = fit_noise(ctx, rng, c["size"], c["agree"][0])
    else:
        c["tau"], c["sigma"] = tau, fit_sigma(ctx, rng, c["size"], tau)
    c["curves"] = r2_curves(ctx, rng, c["sigma"], c["tau"], c["des"])
    c["r2"] = location_r2(c["Sw"][c["des"]["I"], c["des"]["J"]], c["des"])
    return c


class Surrogate:
    """Decision statistics (location-adjusted correlation; transfer gain) and their null from
    position-only surrogates. For each member of the family the matched strength b is where the
    covariates explain as much of the surrogates' similarity as of the data's (the top of the
    strength grid where the member cannot reach it). The member's p-value is the largest over
    b and every grid strength above it (position at least as strong as the data imply: matching
    by explained share is noisy where a field's curve is flat), and the p-value is the largest
    over the family."""

    def __init__(self, c):
        self.ctx, self.sigma, self.tau, self.des, self.curves = c["ctx"], c["sigma"], c["tau"], c["des"], c["curves"]
        ident = np.arange(c["des"]["n"])[None, :]
        self.part = Partial(c["Sg"], c["des"], ident)
        self.tr = Transfer(c["Sg"], c["des"], ident)

    def stats(self, Sw):
        return float(self.part.corr(Sw[self.des["I"], self.des["J"]])[0]), float(self.tr.gain(Sw)[0])

    def strengths(self, r2):
        out = []
        for mem in FAMILY:
            b = interp_up(r2, self.curves[mem], STRENGTHS)
            out.append((mem, b, bool(r2 > max(self.curves[mem]))))
        return out

    def draw(self, rng, mem, b):
        return sim_profiles(self.ctx, rng, self.sigma, tau=self.tau, pos=b, **member_kw(mem))

    def test(self, Sw, rng, n):
        r, g = self.stats(Sw)
        r2 = location_r2(Sw[self.des["I"], self.des["J"]], self.des)
        per = []
        for mem, b, capped in self.strengths(r2):
            best = None
            for bb in [b] + [s for s in STRENGTHS if s > b]:      # position at least as strong as matched
                null = np.array([self.stats((lambda W: W @ W.T)(self.draw(rng, mem, bb)[0])) for _ in range(n)])
                nr = null[:, 0][np.isfinite(null[:, 0])]
                x = {"member": mem, "strength": b, "at": bb, "capped": capped, "n": len(nr),
                     "p_r": (1 + (nr >= r).sum()) / (1 + len(nr)) if np.isfinite(r) else float("nan"),
                     "p_g": (1 + (null[:, 1] >= g).sum()) / (1 + n),
                     "r95": float(np.percentile(nr, 95)) if len(nr) else float("nan"),
                     "r50": float(np.median(nr)) if len(nr) else float("nan")}
                if best is None:
                    best = x
                else:
                    if not (x["p_r"] <= best["p_r"]):           # keep the largest p (NaN-safe)
                        best.update({k: x[k] for k in ("at", "n", "p_r", "r95", "r50")})
                    best["p_g"] = max(best["p_g"], x["p_g"])
            per.append(best)
        return {"r": r, "g": g, "r2": r2, "per": per,
                "p_r": float("nan") if not np.isfinite(r) else max(x["p_r"] for x in per),
                "p_g": max(x["p_g"] for x in per)}


def confusion_curve(c, pos):
    """Mean adjusted r when there is no link to genes but a share of each subclass's cells is
    mislabelled as a gene-similar subclass: what mislabelling alone would produce. Fixed seed,
    so `power` and `run` compute the identical curve."""
    sur = Surrogate(c)
    rng = np.random.default_rng(SEED + 9)
    out = []
    for rate in CONF_RATES:
        r = [sur.stats(sim_similarity(c["ctx"], rng, c["sigma"], tau=c["tau"], pos=pos, field="gp", scale=500.0,
                                      confusion=rate, mix="gene")[0])[0] for _ in range(CONF_REPS)]
        out.append(float(np.nanmean(r)))
    return out


# --- power and size: synthetic projection data only --------------------------------------

_W = {}


def _init(state):
    os.environ["OMP_NUM_THREADS"] = "1"
    _W.update(state)
    c = state["case"]
    _W["sur"] = Surrogate(c)
    _W["relab"] = Partial(c["Sg"], c["des"], state["P"])
    _W["relab_tr"] = Transfer(c["Sg"], c["des"], state["P"])


def _null(job):
    """Surrogate statistics for one family member at one grid strength (power only)."""
    m, k, n, seed = job
    rng = np.random.default_rng(seed)
    sur = _W["sur"]
    return np.array([sur.stats((lambda W: W @ W.T)(sur.draw(rng, FAMILY[m], STRENGTHS[k])[0])) for _ in range(n)])


def grid_p(nulls, b, x, col):
    """p-value at strength b: computed against the stored null at the grid strengths either
    side of b and interpolated linearly between them."""
    k = int(np.clip(np.searchsorted(STRENGTHS, b) - 1, 0, len(STRENGTHS) - 2))
    f = (b - STRENGTHS[k]) / (STRENGTHS[k + 1] - STRENGTHS[k])
    out = []
    for kk in (k, k + 1):
        v = nulls[kk][:, col]
        v = v[np.isfinite(v)]
        out.append((1 + (v >= x).sum()) / (1 + len(v)))
    return (1 - f) * out[0] + f * out[1]


def _one(job):
    kw, seed = job
    kw = dict(kw)
    c = _W["case"]
    rng = np.random.default_rng(seed)
    Sw, _ = sim_similarity(c["ctx"], rng, kw.pop("sigma", c["sigma"]), **dict({"tau": c["tau"]}, **kw))
    sur, d = _W["sur"], c["des"]
    r, g = sur.stats(Sw)
    r2 = location_r2(Sw[d["I"], d["J"]], d)
    p_r, p_g = [], []
    for m, (mem, b, _) in enumerate(sur.strengths(r2)):
        above = [k for k, s in enumerate(STRENGTHS) if s > b]
        p_r.append(max([grid_p(_W["nulls"][m], b, r, 0)] + [grid_p(_W["nulls"][m], STRENGTHS[k], r, 0) for k in above])
                   if np.isfinite(r) else np.nan)
        p_g.append(max([grid_p(_W["nulls"][m], b, g, 1)] + [grid_p(_W["nulls"][m], STRENGTHS[k], g, 1) for k in above]))
    pr = D.one_sided_p(_W["relab"].corr(Sw[d["I"], d["J"]]))
    pg = D.one_sided_p(_W["relab_tr"].gain(Sw))
    return r, max(p_r), max(p_g), pr, pg


def power_rows(c, models, seed):
    """Each model: POWER_OUTER synthetic datasets through the test. To make this affordable the
    surrogate null is drawn once per family member at each grid strength (POWER_NULL surrogates)
    and a dataset's p-value is interpolated between the two grid strengths around its matched
    strength; the run itself draws fresh surrogates at the exact matched strength."""
    from concurrent.futures import ProcessPoolExecutor
    state = {"case": c, "P": D.relabellings(c["des"]["zone"], POWER_PERM, SEED + 1)}
    workers = max(1, (os.cpu_count() or 2) - 1)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init, initargs=(state,)) as ex:
        jobs, ss = [], np.random.SeedSequence(seed + 1000).generate_state(len(FAMILY) * len(STRENGTHS) * workers)
        per = -(-POWER_NULL // workers)
        for m in range(len(FAMILY)):
            for k in range(len(STRENGTHS)):
                for w in range(workers):
                    jobs.append((m, k, per, int(ss[len(jobs)])))
        res = list(ex.map(_null, jobs))
        nulls = [[np.vstack(res[(m * len(STRENGTHS) + k) * workers:(m * len(STRENGTHS) + k + 1) * workers])
                  for k in range(len(STRENGTHS))] for m in range(len(FAMILY))]
    state["nulls"] = nulls
    seeds = np.random.SeedSequence(seed).generate_state(len(models) * POWER_OUTER)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init, initargs=(state,)) as ex:
        for m, (label, kw) in enumerate(models):
            jobs = [(kw, int(s)) for s in seeds[m * POWER_OUTER:(m + 1) * POWER_OUTER]]
            out = np.array(list(ex.map(_one, jobs, chunksize=8)), float)
            prim = out[:, 1] < ALPHA
            print("  %-46s %+8.3f %8.3f %8.3f %10.3f %9.3f" % (label, np.nanmean(out[:, 0]), prim.mean(),
                                                             (out[:, 1] < ALPHA / 2).mean(), (out[:, 3] < ALPHA).mean(),
                                                             (prim & (out[:, 2] < ALPHA)).mean()))
            sys.stdout.flush()


def header():
    import pandas
    import scipy
    print("python %s, numpy %s, pandas %s, scipy %s; seed %d%s" % (
        sys.version.split()[0], np.__version__, pandas.__version__, scipy.__version__, SEED,
        "; QUICK (smoke test, numbers meaningless)" if QUICK else ""))


def primary(t):
    idx = select(t)
    return calibrate(case(t, idx), SEED)


def describe(c, t):
    des, idx = c["des"], c["idx"]
    print("subclasses compared %d (>= %d labelled neurons over all experiments, >= %.0f%% of MERFISH cells in HY); "
          "neurons %d, of them %d in the %d primary experiments"
          % (len(idx), MIN_NEURONS, 100 * MIN_HY_SHARE, t["counts"][idx].sum(), c["counts"].sum(), c["counts"].shape[1]))
    for z in ZONES:
        print("  %-4s %s" % (z, "; ".join(c["names"][des["zone"] == z])))
    print("pairs compared %d; not compared (share methylation clusters >= %.2f): %d; covariates %d"
          % (len(des["pairs"]), OVERLAP, int(np.triu(des["overlap"], 1).sum()), des["D"].shape[1]))
    a = c["agree"]
    print("measured profiles, read without genes: size %.3f; the sexes agree at r = %.3f (centred log shares) "
          "and r = %.3f (pair similarities; both sexes pooled, Spearman-Brown: %.3f)"
          % (c["size"], a[0], a[1], 2 * a[1] / (1 + a[1])))
    print("fitted: experiment noise tau %.2f%s, sigma %.3f" % (
        c["tau"], " (top of the grid)" if c["tau"] >= TAUS[-1] else "", c["sigma"]))
    print("  tau: sigma, synthetic between-sex r   " + "  ".join("%.2f: %.2f, %.3f" % r for r in c["noise_rows"]))
    print("covariates explain %.1f%% of the block-centred similarity; by position strength (mean over %d "
          "synthetic datasets) and the strength matched to the measured share:" % (100 * c["r2"], CURVE_REPS))
    sur = Surrogate(c)
    for mem, b, capped in sur.strengths(c["r2"]):
        print("  %-22s %s  -> %.2f%s" % (member_name(mem), " ".join("%4.1f" % (100 * v) for v in c["curves"][mem]), b,
                                         " (cannot reach it: top of the grid)" if capped else ""))
    return sur


def describe_and_power():
    header()
    t = load()
    c = primary(t)
    sur = describe(c, t)
    ctx, des = c["ctx"], c["des"]
    rng = np.random.default_rng(SEED + 4)
    matched = sur.strengths(c["r2"])[2][1]                 # FAMILY[2]: the 500 um smooth field
    gm = {"pos": matched, "field": "gp", "scale": 500.0}
    s_tau0 = fit_sigma(ctx, rng, c["size"], 0.0)
    s_tau2 = fit_sigma(ctx, rng, c["size"], 2 * c["tau"])
    gene_b = {a: match_pos(ctx, rng, c["sigma"], c["tau"], des, c["r2"], gene=a) for a in (0.5, 0.75, 1.0, 1.5, 2.0)}
    print("matched 500 um strength %.2f; re-matched beside a gene term: %s; sigma at tau 0: %.3f, at tau %.2f: %.3f"
          % (matched, ", ".join("%.2f -> %.2f" % kv for kv in gene_b.items()), s_tau0, 2 * c["tau"], s_tau2))

    models = [("noise only", {})]
    models += [("bumps %d um" % sc, {"pos": 1.0, "scale": sc}) for sc in (250.0, 500.0, 1000.0)]
    models += [("bumps %d um strong" % sc, {"pos": 3.0, "scale": sc}) for sc in (250.0, 1000.0)]
    models += [("bumps 500 um at the matched strength", {"pos": matched, "scale": 500.0})]
    models += [("smooth field %d um strong" % sc, {"pos": 3.0, "scale": sc, "field": "gp"}) for sc in (100.0, 250.0, 1000.0)]
    models += [("per-structure field strong", {"pos": 3.0, "field": "parcel"}),
               ("outside the family: smooth 150 um strong", {"pos": 3.0, "scale": 150.0, "field": "gp"}),
               ("outside the family: smooth 2000 um", {"pos": 1.0, "scale": 2000.0, "field": "gp"}),
               ("outside the family: bumps 100 um strong", {"pos": 3.0, "scale": 100.0}),
               ("dorsal-ventral gradient 1.0", {"pos": 1.0, "field": "linear", "axes": (1,)}),
               ("medial-lateral gradient 1.0", {"pos": 1.0, "field": "linear", "axes": (2,)}),
               ("3D linear gradient 1.0", {"pos": 1.0, "field": "linear"}),
               ("dorsal-ventral 0.75 + smooth 500 um 0.74", {"pos": 0.75, "field": "linear", "axes": (1,), "pos2": 0.74}),
               ("matched, experiment noise 0", dict(gm, tau=0.0, sigma=s_tau0)),
               ("matched, experiment noise doubled", dict(gm, tau=2 * c["tau"], sigma=s_tau2))]
    models += [("matched + %d%% moved along shared clusters" % (100 * r), dict(gm, confusion=r)) for r in (0.05, 0.10, 0.20)]
    models += [("matched + %d%% moved to gene neighbours" % (100 * r), dict(gm, confusion=r, mix="gene"))
               for r in (0.05, 0.10, 0.20, 0.30)]
    models += [("matched + 20% moved to gene neighbours, broad", dict(gm, confusion=0.20, mix="gene", temp=0.3))]
    models += [("genes %.2f + re-matched position" % a, {"pos": gene_b[a], "field": "gp", "scale": 500.0, "gene": a})
               for a in (0.5, 0.75, 1.0, 1.5, 2.0)]
    print("rates over %d synthetic datasets, one-sided %.2f; null: %d surrogates per family member at each grid "
          "strength, p interpolated between grid strengths; each member takes its largest p at the matched strength or any grid strength above. primary = surrogate null "
          "(largest p over the family of %d); relabelled = D-099's null; transfer = primary AND transfer gain "
          "against the surrogates" % (POWER_OUTER, ALPHA, POWER_NULL, len(FAMILY)))
    print("  %-46s %8s %8s %8s %10s %9s" % ("model", "mean r", "primary", "p<0.025", "relabelled", "transfer"))
    power_rows(c, models, SEED + 3)

    cz = calibrate(case(t, c["idx"], same_zone=True), SEED, tau=c["tau"])
    mz = Surrogate(cz).strengths(cz["r2"])[2][1]
    print("same-zone pairs only (%d pairs; covariates explain %.1f%%; matched strength %.2f):"
          % (len(cz["des"]["pairs"]), 100 * cz["r2"], mz))
    power_rows(cz, [("same zone: noise only", {}),
                    ("same zone: matched position", {"pos": mz, "field": "gp", "scale": 500.0}),
                    ("same zone: genes 1.00 + matched position", {"pos": mz, "field": "gp", "scale": 500.0, "gene": 1.0})],
               SEED + 5)

    # what leaving a zone out does to r when the link is the same everywhere (reference for s8)
    ident = np.arange(des["n"])[None, :]
    drops = {z: design(t, c["idx"], c["n"], drop_zone=z) for z in ZONES if (des["zone"] == z).any()}
    parts = {z: Partial(c["Sg"], d, ident) for z, d in drops.items()}
    rng = np.random.default_rng(SEED + 6)
    full, part = [], {z: [] for z in drops}
    for _ in range(POWER_OUTER):
        Sw, _ = sim_similarity(ctx, rng, c["sigma"], tau=c["tau"], pos=gene_b[1.5], field="gp", scale=500.0, gene=1.5)
        full.append(sur.stats(Sw)[0])
        for z, d in drops.items():
            part[z].append(float(parts[z].corr(Sw[d["I"], d["J"]])[0]))
    print("each zone left out, a uniform link (genes 1.5): mean r %+.3f with all zones; ratio without %s"
          % (np.nanmean(full), ", ".join("%s %.2f" % (z, np.nanmean(v) / np.nanmean(full)) for z, v in part.items())))
    cc = confusion_curve(c, matched)
    print("no link, cells moved to gene neighbours (%d datasets each): mean r %s"
          % (CONF_REPS, ", ".join("%d%%: %+.3f" % (100 * r, v) for r, v in zip(CONF_RATES, cc))))


# --- the run ------------------------------------------------------------------------

def git(*args):
    r = subprocess.run(["git", "-C", ROOT] + list(args), capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("refusing to run: git %s failed: %s" % (" ".join(args), r.stderr.strip()))
    return r.stdout.strip()


def input_hashes():
    import hashlib
    out = {}
    for name, p in list(RAW.items()) + [("tables", os.path.relpath(TABLES, ROOT)),
                                        ("outside_rule", os.path.relpath(OUTSIDE_RULE, ROOT))]:
        h = hashlib.sha256()
        if not os.path.exists(os.path.join(ROOT, p)):
            out[name] = "missing"
            continue
        with open(os.path.join(ROOT, p), "rb") as f:
            for block in iter(lambda: f.read(1 << 22), b""):
                h.update(block)
        out[name] = h.hexdigest()
    return out


PINNED = {
    "st2": "4855044e38ea6aa782a5a6cc9379d7da5e8bfd83c42dacccdd6b2bf42de7b1cf",
    "stats": "b21ca985652fb25f9608f99005139a40757133a76fbe845ae5b175c5c26a447b",
    "terms": "653606461801393cd0a7708a084e3f9c989ffb2ac696ac30629f6ef6180e76c5",
    "genes_10x": "ec5a211c23fab5cce8a06f01bb0a68af087fb979b647ab42d635201f98de2f69",
    "genes_merfish": "410150d9d48cf9e0ca41a6d0ba085471f693d54b75f67ec67b3635a6cb3afd17",
    "merfish": "1e3ae23cc3f8d3d5839cc8a222798ff812fc66f4bd3c0fc8fca0c0279329874e",
    "structure_graph": "b5f0d024d1df09ee18aef15093f55240f4a70b3924ee9689add612dadff80386",
    "tables": "a8abf4bce1af37240ec9fb48b7c0f17570119d8ffade614e9abed9c002d14802",
    "outside_rule": "40b4775925bb026796aa62a149c9c2996c313ab93a184b67f01e0b0efc8457f5",
}


def analyse(t, idx, seed, tau, n_surr=None, **kw):
    """Decision statistics and surrogate p-values for one variant of the comparison."""
    c = calibrate(case(t, idx, **kw), seed, tau=tau)
    res = Surrogate(c).test(c["Sw"], np.random.default_rng(seed + 2), n_surr or N_SURR_SECONDARY)
    res["pairs"], res["units"] = len(c["des"]["pairs"]), len(set(c["des"]["I"]) | set(c["des"]["J"]))
    return res, c


def synthetic_table(t):
    """For `dry`: one synthetic count table (all subclasses, all experiments) drawn from the
    real table's marginals, with a known gene link, built without reading the real profiles."""
    rng = np.random.default_rng(SEED + 99)
    t0 = dict(t)
    tot = t["counts"].sum()
    t0["counts"] = np.outer(t["counts"].sum(axis=1), t["counts"].sum(axis=0)) / tot       # marginals only
    alln = np.arange(len(t["names"]))
    c = case(t0, alln, targets=TARGETS + SINGLE_SEX)
    c["ctx"]["N_e"] = t["counts"].sum(axis=0).astype(int)
    out = dict(t)
    out["counts"] = simulate(c["ctx"], rng, 2.0, tau=0.5, pos=0.74, field="gp", scale=500.0, gene=1.5).astype(int)
    out["counts_pure"] = rng.binomial(out["counts"], 0.75)
    return out


def row(label, fn):
    try:
        print("  %-52s %s" % (label, fn()))
    except Exception as exc:                      # a secondary row never stops the run
        print("  %-52s could not be computed: %s: %s" % (label, type(exc).__name__, exc))
    sys.stdout.flush()


def fmt(v):
    return "r = %+.3f  p = %.4f  (%d subclasses, %d pairs)" % (v["r"], v["p_r"], v["units"], v["pairs"])


def outside_rule(t, c, sur, rng):
    """s11: a ridge rule from genes to projections fitted on every non-hypothalamic source
    (analysis/d104_rehearsal.py freeze), applied unchanged to the compared subclasses."""
    R = np.load(OUTSIDE_RULE, allow_pickle=True)
    gi = R["gene_idx"]
    assert list(t["symbols"][gi]) == list(R["symbols"]), "gene order differs from the frozen rule"
    cols = [list(R["targets"]).index(tg) for tg in c["spec"]["targets"]]
    X = ((t["means"][c["idx"]][:, gi] - R["mean"]) / R["sd"]) @ R["V"]
    G = X - X.mean(axis=0)
    B = R["B"][:, cols]
    pred = double_centre(G @ (B - B.mean(axis=1, keepdims=True)))
    w = np.repeat(np.sqrt(c["n"] / c["n"].mean()), pred.shape[1])
    T1 = np.column_stack([np.ones(len(c["idx"])), c["des"]["std_pos"]])

    def detrend(A):
        return A - T1 @ np.linalg.lstsq(T1, A, rcond=None)[0]

    out = []
    for label, f in (("", lambda A: A), (", linear position trend removed", lambda A: double_centre(detrend(A)))):
        P = unitrows(f(pred))
        obs = wcorr(P.ravel(), unitrows(f(c["L"])).ravel(), w)
        perm = np.array([wcorr(P[rng.permutation(len(P))].ravel(), unitrows(f(c["L"])).ravel(), w) for _ in range(N_PERM)])
        ps = []
        for mem, b, _ in sur.strengths(c["r2"]):
            nul = np.array([wcorr(P.ravel(), unitrows(f(sur.draw(rng, mem, b)[1])).ravel(), w) for _ in range(N_SURR_SECONDARY)])
            ps.append((1 + (nul >= obs).sum()) / (1 + len(nul)))
        out.append("r_unit = %+.3f  permutation p = %.4f  surrogate p = %.4f%s"
                   % (obs, (1 + (perm >= obs).sum()) / (1 + N_PERM), max(ps), label))
    return ("\n" + " " * 55).join(out) + "  (%d components, ridge %g)" % (int(R["k"]), float(R["lam"]))


def run(dry=False):
    header()
    t = load()
    if dry:
        print("DRY RUN: the measured counts are replaced by one synthetic table (genes 1.5, position 0.74, "
              "experiment noise 0.5). NOT A RESULT. Nothing below describes the measured projections.\n")
        t = synthetic_table(t)
    else:
        if QUICK:
            sys.exit("refusing to run: D104_QUICK is set")
        hashes = input_hashes()
        assert set(PINNED) == set(hashes), "PINNED must hold every input"
        differs = sorted(k for k, v in hashes.items() if PINNED.get(k) != v)
        if differs:
            sys.exit("refusing to run: inputs differ from the pinned ones: %s" % ", ".join(differs))
        dirty = git("status", "--porcelain", "--", "analysis", "src")
        if dirty:
            sys.exit("refusing to run: uncommitted changes to the code:\n" + dirty)
        head = git("rev-parse", "HEAD")
        if not head:
            sys.exit("refusing to run: no commit")
        print("commit %s; every input matches its pinned SHA-256\n" % head)

    c = primary(t)
    sur = describe(c, t)
    des, idx, Sg, Sw, tau = c["des"], c["idx"], c["Sg"], c["Sw"], c["tau"]
    res = sur.test(Sw, np.random.default_rng(SEED + 2), N_SURR)
    p1 = res["p_r"]
    passed = bool(p1 < HOLM[1])           # can pass under Holm (at 0.05 only if D-105's p < 0.025)
    verdict = ("PASS (below %.3f: passes under Holm whatever D-105's p)" % HOLM[0] if p1 < HOLM[0] else
               "PASS only if D-105's primary p < %.3f (Holm); otherwise not significant" % HOLM[0] if passed else
               "not significant")
    print("\nDECISION 1 - a cell type's genes predict where its neurons project, position partialled out")
    print("  Holm with D-105 at family-wise %.2f: the smaller one-sided p must be below %.3f, the other then below %.2f"
          % (ALPHA, HOLM[0], HOLM[1]))
    print("  %d subclasses, %d pairs: r = %+.3f   one-sided p = %.4f (largest p over the family; %d surrogates each)  -> %s"
          % (des["n"], len(des["pairs"]), res["r"], p1, N_SURR, verdict))
    for x in res["per"]:
        print("    %-22s matched strength %.2f%s; largest p at strength %.2f: null median %+.3f, 95th %+.3f, p = %.4f (%d surrogates)"
              % (member_name(x["member"]), x["strength"], " (top of the grid)" if x["capped"] else "", x["at"], x["r50"],
                 x["r95"], x["p_r"], x["n"]))
    matched = sur.strengths(c["r2"])[2][1]
    cc = confusion_curve(c, matched)
    print("  no link, cells moved to gene neighbours: mean r %s"
          % ", ".join("%d%%: %+.3f" % (100 * r, v) for r, v in zip(CONF_RATES, cc)))
    if np.isfinite(res["r"]):
        eq = "above 30%" if res["r"] > max(cc) else "%.0f%%" % (100 * interp_up(res["r"], cc, CONF_RATES))
        print("  mislabelling that would give this r with no link: %s%s" % (eq, (
            "; the pass is NOT separable from mislabelling (r at or below the 20% row)" if res["r"] <= cc[3] else
            "; r is above the 20% row") if passed else ""))
    print("DECISION 2 - transfer from the most gene-similar same-zone subclass beats the zone average")
    print("  mean profile-cosine gain %+.4f   one-sided p = %.4f  -> %s"
          % (res["g"], res["p_g"], ("transfer supported" if res["p_g"] < ALPHA else "transfer not shown")
             + ("" if p1 < HOLM[0] else ", if Decision 1 passes under Holm") if passed
             else "not tested (decision 1 did not pass)"))

    print("\nSECONDARY (descriptive; no decision; %d surrogates per family member, largest p)" % N_SURR_SECONDARY)
    y = Sw[des["I"], des["J"]]
    row("s1 D-099's null (genes relabelled within zone)", lambda: "r = %+.3f  p = %.4f" % (
        res["r"], D.one_sided_p(Partial(Sg, des, D.relabellings(des["zone"], N_PERM, SEED + 10)).corr(y))))
    row("s2 zone blocks only (no position or shared-methylation covariates)", lambda: (lambda v: "r = %+.3f  p = %.4f (relabelled)" % (
        v[0], D.one_sided_p(v)))(Partial(Sg, des, D.relabellings(des["zone"], N_PERM, SEED + 11), D=None).corr(y)))
    variants = [
        ("s3 all 10 targets, single-sex ones sex-adjusted", idx, {"targets": TARGETS + SINGLE_SEX, "sex_adjust": True}),
        ("s3u all 10 targets, unadjusted", idx, {"targets": TARGETS + SINGLE_SEX}),
        ("s3f female experiments only (7 targets)", idx, {"sex": "female"}),
        ("s3m male experiments only (7 targets)", idx, {"sex": "male"}),
        ("s4 centroids from the slab (CCF AP %.1f-%.1f mm)" % SLAB_MM, idx, {"centroid": "centroid_slab"}),
        ("s5 MERFISH panel genes as the fingerprint", idx, {"genes": "panel"}),
        ("s6 subclasses with >= 40 labelled neurons", select(t, min_neurons=40), {}),
        ("s7 the three non-hypothalamic subclasses kept", select(t, min_hy=0.0), {}),
        ("s10 same transmitter (Glut/Gaba) also partialled", idx, {"same_transmitter": True}),
        ("s13 same-zone pairs only", idx, {"same_zone": True}),
        ("s14 %s placed in %s" % NEAR_TIE, idx,
         {"zone": np.where(t["names"][idx] == NEAR_TIE[0], NEAR_TIE[1], t["zone"][idx])}),
        ("s15 cleanly clustered neurons only (>= %d per subclass)" % MIN_PURE,
         idx[t["counts_pure"][idx].sum(axis=1) >= MIN_PURE], {"counts": "counts_pure"}),
    ]
    variants += [("s8 without zone %s" % z, idx, {"drop_zone": z}) for z in ZONES if (des["zone"] == z).any()]
    variants += [("s12 without target %s" % tg, idx, {"targets": tuple(x for x in TARGETS if x != tg)}) for tg in TARGETS]
    for k, (label, sel, kw) in enumerate(variants):
        row(label, lambda sel=sel, kw=kw, k=k: fmt(analyse(t, sel, SEED + 20 + k, tau, **kw)[0]))
    print("  s8 is read against the ratios in reports/d104_power.md (what leaving a zone out does to a uniform link)")

    def s9():
        pos = des["pos"]
        allp = [(i, j) for i in range(des["n"]) for j in range(i + 1, des["n"])]
        close = -np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(axis=2))
        return "; ".join("%s r = %+.3f p = %.4f" % ((lab,) + D.mantel(A, B, allp, N_PERM, SEED + 50)) for lab, A, B in (
            ("genes/projections", Sg, Sw), ("closeness/projections", close, Sw), ("closeness/genes", close, Sg)))
    row("s9 all pairs, nothing partialled (Mantel)", s9)
    row("s11 rule learned outside the hypothalamus", lambda: outside_rule(t, c, sur, np.random.default_rng(SEED + 60)))
    ident = np.arange(des["n"])[None, :]

    def s16():
        out = []
        for u in range(des["n"]):
            d = design(t, idx, c["n"], drop_unit=u)
            out.append(float(Partial(Sg, d, ident).corr(Sw[d["I"], d["J"]])[0]))
        lo, hi = int(np.nanargmin(out)), int(np.nanargmax(out))
        return "r from %+.3f (without %s) to %+.3f (without %s); all: %s" % (
            out[lo], c["names"][lo], out[hi], c["names"][hi], " ".join("%+.2f" % v for v in out))
    row("s16 each subclass left out (r only)", s16)


if __name__ == "__main__":
    cmd = sys.argv[1:2]
    if cmd == ["tables"]:
        build_tables()
    elif cmd == ["power"]:
        describe_and_power()
    elif cmd == ["dry"]:
        run(dry=True)
    elif cmd == ["run"]:
        run()
    else:
        sys.exit(__doc__)
