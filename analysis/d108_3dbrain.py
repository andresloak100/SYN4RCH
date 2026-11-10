"""D-108 3D brain: a reusable 3D assembly of the Allen MERFISH atlas in CCF space (folded onto one
hemisphere at z = 5.7 mm, as d105), built from LOCAL atlas data only:

  data/raw/abc/merfish_cells.csv        3.74 M cells: CCF x/y/z, subclass, parcellation category /
                                        division / structure / substructure (fibre-tract labels included)
  data/raw/atlas/structure_graph_1.json the CCF ontology (used to check every tract label is a
                                        descendant of 'fiber tracts', id 1009)
  data/raw/abc/precomputed_stats.h5     Allen 10x cluster means (only for the target-side mean gene
                                        embedding of F4, through D's frozen whole-brain PCA)
  data/interim/d105_models.npz, d108_prep.npz   frozen D embedding, inferred locations (module M1)

EXCLUDED (projection-derived anatomy): data/raw/atlas/connectivity (Allen tracing experiments) and
every Epi-retro-seq row; this file never opens either.

What the brain contains (data/interim/d108_3d_brain.npz + d108_3d_brain_card.json):
  * a 200 um voxel grid: in-brain mask, cells per voxel, fibre-tract fraction, tract system and
    tract structure label per voxel, a flag for voxels filled by interpolation between sections;
  * a voxel route graph (26-neighbour) and, per target region, route-cost fields: tract voxels
    cheap (cost factor 0.25; sensitivity 0.1 and 0.5) and uniform cost (brain-constrained geodesic);
  * every (source region, subclass) unit's REAL cell cloud (not one centroid) for the 27 non-HY
    Epi-retro-seq sources and for the 19 HY test subclasses in HY (atlas cells only, no projection
    data), and the feature families F2-F7 per unit and per (unit, target).
Every feature carries a tag: M = measured (atlas cells and their registered labels), I = inferred
(with an uncertainty), P = projection-derived (none; excluded by construction).

API (reusable for an unseen region):  B = load(); cloud = cloud_cells(B, subclass, region_spec);
  feats = unit_features(B, cloud_xyz, subclass, targets)
"""
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import d105_benchmark as b5  # noqa: E402

ROOT = b5.ROOT
INTERIM = os.path.join(ROOT, "data/interim")
CELLS = os.path.join(INTERIM, "d108_3d_cells.npz")
BRAIN = os.path.join(INTERIM, "d108_3d_brain.npz")
CARD = os.path.join(INTERIM, "d108_3d_brain_card.json")
SEED = 20261110

VOX = 0.2                 # mm
MASK_R = 0.45             # mm: a voxel is in the brain if >= MASK_MIN atlas cells lie within MASK_R of its centre
MASK_MIN = 3
KNN_VOX = 16              # cells used to label a voxel
TRACT_FRAC = 0.5          # tract voxel: >= half of its labelling cells are fibre-tract cells
COST_TRACT = 0.25         # route cost factor in tract voxels (grey = 1); sensitivity values below
COST_SENS = (0.1, 0.5)
CLOUD_MAX = 3000          # cells per unit cloud used for geometry (seeded subsample)
NB_K = 30                 # neighbourhood: 30 nearest atlas cells
NB_PCS = 6
NB_REF = 30000
TRACT_R = 0.3             # mm: a cell "has access" to a tract if a tract voxel centre is within 0.3 mm
OVL_R = 0.2               # mm: overlap radius (source cell within 0.2 mm of a target-region cell)
FLOOR = 0.05              # mm, as d105 A_MIN_MM
LHAT_SAMPLES = 200
MIN_TRACT_VOX = 20        # tract structures with fewer voxels are not used as corridors

TAGS = {
    "F1": "M  measured centroid of the unit's MERFISH cells and of the target region's cells",
    "F2": "M  measured geometry of the unit's cell cloud relative to the target region's cells",
    "F3": "M  measured local cell-type composition (30 nearest atlas cells) around the unit's cells",
    "F4": "M  target composition (atlas cells) x unit gene embedding (atlas 10x); compatibility is a fitted low-rank term",
    "F5": "I  tract voxels are measured labels (cells registered to CCF fibre-tract parcels); reading them as corridors is an "
          "assumption; uncertainty: voxel label interpolated between sections flagged",
    "F6": "I  route cost through a voxel graph with cheap tract voxels: the cost model is assumed (cost factor 0.25, "
          "sensitivity 0.1 / 0.5 stored as the uncertainty), not a known axon route",
    "F7": "I  placement inferred from genes (module M1, out-of-fold), uncertainty: M1 pred_sd propagated by 200 samples",
}


def say(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------------------- cells
def read_cells():
    """MERFISH cells -> compact table (mm, z folded at the midline as d105)."""
    import pandas as pd
    if os.path.exists(CELLS):
        z = np.load(CELLS, allow_pickle=True)
        return {k: z[k] for k in z.files}
    cols = ["subclass", "x_ccf", "y_ccf", "z_ccf", "parcellation_category", "parcellation_division",
            "parcellation_structure", "parcellation_substructure"]
    m = pd.read_csv(b5.path(b5.RAW["merfish"]), usecols=cols)
    out = {"xyz": np.column_stack([m["x_ccf"], m["y_ccf"], np.abs(m["z_ccf"] - b5.MIDLINE_Z_MM)]).astype(np.float32)}
    sub = m["subclass"].map(b5.strip)
    for name, col in (("sub", sub), ("cat", m["parcellation_category"]), ("div", m["parcellation_division"]),
                      ("st", m["parcellation_structure"]), ("sst", m["parcellation_substructure"])):
        codes, uniq = pd.factorize(col.astype(str), sort=True)
        out[name] = codes.astype(np.int16)
        out[name + "_names"] = np.array(list(uniq))
    np.savez(CELLS, **out)
    return out


class Mer(b5.Merfish):
    """b5.Merfish interface (unit_centroid, target_centroid, _structs) over the compact table, so the
    unit clouds and target regions are selected exactly as d105 / d108 Phase 1 select them."""

    def __init__(self, C):
        self.xyz = C["xyz"].astype(float) * 1000.0
        self.div = C["div_names"][C["div"]]
        self.st = C["st_names"][C["st"]]
        self.allst = sorted(set(C["st_names"]))
        sub = C["sub_names"][C["sub"]]
        import pandas as pd
        self.idx = pd.Series(np.arange(len(sub))).groupby(sub).indices

    def unit_cells(self, sub, src):
        """indices of the unit's cloud: d105 unit_centroid's selection, all cells (not the centroid)."""
        ix = self.idx.get(sub, np.array([], int))
        if src == "HY":
            sel, level = self.div[ix] == "HY", "division"
        elif src in b5.SRC_STRUCT:
            sel, level = np.isin(self.st[ix], self._structs(b5.SRC_STRUCT[src])), "parcel"
        else:
            sel, level = np.isin(self.div[ix], b5.SRC_DIV[src]), "parcel"
        if sel.sum() < 20 and src != "HY":
            sel, level = np.isin(self.div[ix], b5.DIV_CCF[b5.DIVISION[src]]), "division"
        if sel.sum() < 20:
            sel, level = np.ones(len(ix), bool), "brain"
        return ix[sel], level

    def target_cells(self, t):
        import d108_prims as P
        return np.where(P.target_mask(self, t))[0]


# ------------------------------------------------------------------------- the voxel brain
def voxelise(C, ontology_ok):
    from scipy.spatial import cKDTree
    t0 = time.time()
    cat = C["cat_names"][C["cat"]]
    inb = cat != "unassigned"                                   # 'unassigned' = outside the annotated brain
    xyz = C["xyz"][inb].astype(float)
    istract = (cat[inb] == "fiber tracts")
    tdiv = C["div"][inb]
    tst = C["st"][inb]
    shape = tuple(int(np.ceil(xyz[:, a].max() / VOX)) + 1 for a in range(3))
    g = np.stack(np.meshgrid(*[np.arange(n) for n in shape], indexing="ij"), -1).reshape(-1, 3)
    centres = (g + 0.5) * VOX
    tree = cKDTree(xyz)
    d, nn = tree.query(centres, k=KNN_VOX, distance_upper_bound=MASK_R, workers=4)
    found = np.isfinite(d)
    nfound = found.sum(1)
    nnc = np.where(found, nn, 0)
    mask = nfound >= MASK_MIN
    tfrac = np.where(nfound > 0, (istract[nnc] & found).sum(1) / np.maximum(nfound, 1), 0.0)
    tract = mask & (tfrac >= TRACT_FRAC)
    # majority tract system (division label) and structure among the tract cells of each tract voxel
    sysl = np.full(len(g), -1)
    stl = np.full(len(g), -1)
    for v in np.where(tract)[0]:
        sel = found[v] & istract[nnc[v]]
        sysl[v] = np.bincount(tdiv[nnc[v][sel]]).argmax()
        stl[v] = np.bincount(tst[nnc[v][sel]]).argmax()
    # cells per voxel (measured) vs voxels labelled only from cells in neighbouring sections (interpolated)
    vi = np.floor(xyz / VOX).astype(int)
    lin = np.ravel_multi_index(vi.T, shape)
    ncell = np.bincount(lin, minlength=len(g))
    interp = mask & (ncell == 0)
    # keep the largest 26-connected component of the mask
    from scipy import ndimage
    lab, nlab = ndimage.label(mask.reshape(shape), structure=np.ones((3, 3, 3)))
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    keep = (lab == sizes.argmax()).ravel()
    say("voxel grid %s (%.1f mm): %d in-brain voxels (%d components; largest %d kept), %d tract voxels, %d interpolated  [%.0fs]"
        % (shape, VOX, mask.sum(), nlab, keep.sum(), (tract & keep).sum(), (interp & keep).sum(), time.time() - t0))
    mask = keep
    tract &= mask
    # tract connectivity
    tl, ntl = ndimage.label(tract.reshape(shape), structure=np.ones((3, 3, 3)))
    ts = np.sort(np.bincount(tl.ravel())[1:])[::-1]
    divn, stn = C["div_names"], C["st_names"]
    sys_counts = {str(divn[a]): int((sysl[tract] == a).sum()) for a in np.unique(sysl[tract])}
    st_counts = {str(stn[a]): int((stl[tract] == a).sum()) for a in np.unique(stl[tract])}
    info = {"shape": shape, "voxel_mm": VOX, "n_mask": int(mask.sum()), "n_tract": int(tract.sum()),
            "n_interpolated": int((interp & mask).sum()), "n_tract_interpolated": int((interp & tract).sum()),
            "n_components_mask_before": int(nlab), "tract_components": int(ntl),
            "tract_component_sizes_top5": ts[:5].tolist(), "tract_system_voxels": sys_counts,
            "tract_structure_voxels": st_counts, "ontology_check": ontology_ok,
            "n_cells_inbrain": int(inb.sum()), "n_tract_cells": int(istract.sum())}
    return {"shape": shape, "mask": mask, "tract": tract, "tfrac": tfrac.astype(np.float32), "ncell": ncell.astype(np.int32),
            "interp": interp & mask, "sys": sysl.astype(np.int16), "stl": stl.astype(np.int16), "info": info, "tree": tree,
            "xyz_inb": xyz, "inb": inb}


def ontology_check(C):
    """every tract division / structure label in the MERFISH table should be a descendant of 'fiber tracts'."""
    g = json.load(open(os.path.join(ROOT, "data/raw/atlas/structure_graph_1.json")))
    root = g["msg"][0]
    desc = set()

    def walk(n, under):
        u = under or n["acronym"] == "fiber tracts"
        if u:
            desc.add(n["acronym"])
        for c in n.get("children", []):
            walk(c, u)
    walk(root, False)
    cat = C["cat_names"][C["cat"]]
    m = cat == "fiber tracts"
    out = {}
    for name in ("div", "st", "sst"):
        labs = sorted(set(C[name + "_names"][np.unique(C[name][m])]))
        labs = [x for x in labs if not x.endswith("unassigned")]
        out[name] = {"n_labels": len(labs), "not_under_fiber_tracts": [x for x in labs if x not in desc]}
    return out


def graph(vb, cost_tract):
    """26-neighbour voxel graph over the mask; edge = length x mean cost factor of its two voxels."""
    from scipy.sparse import coo_matrix
    shape = vb["shape"]
    mask = vb["mask"].reshape(shape)
    node = -np.ones(shape, int)
    node[mask] = np.arange(mask.sum())
    c = np.where(vb["tract"], cost_tract, 1.0).reshape(shape)
    rows, cols, w = [], [], []
    offs = [(a, b, e) for a in (-1, 0, 1) for b in (-1, 0, 1) for e in (-1, 0, 1) if (a, b, e) > (0, 0, 0)]
    for a, b, e in offs:
        sl0 = tuple(slice(max(0, -o), n - max(0, o)) for o, n in zip((a, b, e), shape))
        sl1 = tuple(slice(max(0, o), n - max(0, -o)) for o, n in zip((a, b, e), shape))
        m = mask[sl0] & mask[sl1]
        L = VOX * np.sqrt(a * a + b * b + e * e)
        rows.append(node[sl0][m])
        cols.append(node[sl1][m])
        w.append(L * 0.5 * (c[sl0][m] + c[sl1][m]))
    r = np.concatenate(rows)
    cc = np.concatenate(cols)
    ww = np.concatenate(w)
    n = int(mask.sum())
    return coo_matrix((ww, (r, cc)), shape=(n, n)).tocsr()


def vox_of(vb, xyz):
    """mask-voxel node id for each point (nearest mask voxel when the point's voxel is outside the mask)."""
    from scipy.spatial import cKDTree
    shape = vb["shape"]
    mid = np.where(vb["mask"])[0]
    node = -np.ones(len(vb["mask"]), int)
    node[mid] = np.arange(len(mid))
    vi = np.clip(np.floor(np.asarray(xyz) / VOX).astype(int), 0, np.array(shape) - 1)
    lin = np.ravel_multi_index(vi.T, shape)
    out = node[lin]
    bad = out < 0
    if bad.any():
        if "_mtree" not in vb:
            vb["_mtree"] = cKDTree((np.array(np.unravel_index(mid, shape)).T + 0.5) * VOX)
        _, j = vb["_mtree"].query(np.asarray(xyz)[bad])
        out[bad] = j
    return out


# ------------------------------------------------------------------------- features
def nb_comp(tree, sub_inb, nsub, pts, k=NB_K):
    """sqrt of the subclass composition of the k nearest atlas cells of each point (self excluded)."""
    _, nn = tree.query(pts, k=k + 1, workers=4)
    nn = nn[:, 1:]
    lab = sub_inb[nn]
    M = np.zeros((len(pts), nsub))
    np.add.at(M, (np.repeat(np.arange(len(pts)), k), lab.ravel()), 1.0 / k)
    return np.sqrt(M)


def unit_features(B, cloud, sub, targets, rng=None):
    """all 3D features of ONE unit from its atlas cell cloud (mm, folded) and its subclass name, for the
    given target regions; uses the atlas brain only (no projection data). -> dict of arrays."""
    from scipy.spatial import cKDTree
    rng = rng or np.random.default_rng(SEED)
    cloud = np.asarray(cloud, float)
    cen = cloud.mean(0)
    sc = cloud if len(cloud) <= CLOUD_MAX else cloud[rng.choice(len(cloud), CLOUD_MAX, replace=False)]
    TI = [list(B["targets"]).index(t) for t in targets]
    f = {"centroid": cen, "n_cells": len(cloud)}
    # F2 cloud geometry vs each target region's cells
    med, p05, ovl = [], [], []
    for a in TI:
        d, _ = B["ttree"][a].query(sc, workers=4)
        med.append(np.median(d))
        p05.append(np.quantile(d, 0.05))
        ovl.append(float((d <= OVL_R).mean()))
    f["cl_med"], f["cl_p05"], f["cl_ovl"] = np.array(med), np.array(p05), np.array(ovl)
    X = sc - sc.mean(0)
    ev, V = np.linalg.eigh(np.cov(X.T) + 1e-6 * np.eye(3))
    ev, V = ev[::-1], V[:, ::-1]
    sd = np.sqrt(np.maximum(ev, 1e-8))
    f["shape6"] = np.array([np.log(sd[0]), np.log(sd[1] / sd[0]), np.log(sd[2] / sd[0]), *np.abs(V[:, 0])])
    # F3 neighbourhood
    nc = sc if len(sc) <= 800 else sc[rng.choice(len(sc), 800, replace=False)]
    comp = nb_comp(B["tree"], B["sub_inb"], B["nsub"], nc).mean(0)
    f["nb_comp"] = comp
    f["nb6"] = (comp - B["nb_mu"]) @ B["nb_V"]
    tcs = B["tcomp_sqrt"][TI]
    f["nbsim"] = (tcs @ comp) / (np.linalg.norm(tcs, axis=1) * np.linalg.norm(comp) + 1e-12)
    # F4 compatibility: unit gene embedding (D's frame) . target mean gene embedding (D's frame)
    if sub in B["emb_names"]:
        u = B["emb_D"][list(B["emb_names"]).index(sub)]
        f["gcompat"] = B["temb_D"][TI] @ u
    else:
        f["gcompat"] = np.full(len(TI), np.nan)
    # F5 corridor access
    ac = sc if len(sc) <= 1000 else sc[rng.choice(len(sc), 1000, replace=False)]
    f["tsys6"] = np.array([float((B["systree"][s].query(ac, workers=4)[0] <= TRACT_R).mean()) if B["systree"][s] is not None else 0.0
                           for s in range(len(B["sys_names"]))])
    A = np.array([float((t.query(ac, workers=4)[0] <= TRACT_R).mean()) for t in B["sttree"]])
    f["tract_access"] = A
    f["corr_shared"] = B["tract_B"][TI] @ A
    # F6 route cost: mean over cloud cells of the target's route field at the cell's voxel
    nodes = vox_of(B, sc if len(sc) <= 1000 else ac)
    R = B["fields"][:, TI][:, :, nodes].mean(-1)                # (n_cost, n_targets)
    f["route_costs"] = R
    f["route"] = -np.log(np.maximum(R[B["cost_names"].index("tract%.2f" % COST_TRACT)], FLOOR))
    f["route_sd"] = np.std([-np.log(np.maximum(R[B["cost_names"].index("tract%.2f" % c)], FLOOR)) for c in (COST_TRACT,) + COST_SENS], axis=0)
    geo = R[B["cost_names"].index("uniform")]
    f["route_gain"] = np.log(np.maximum(geo, FLOOR)) + f["route"]          # log(geodesic / tract-weighted cost)
    eu = np.sqrt(((cen[None] - B["tcen"][TI]) ** 2).sum(-1))
    f["euclid"] = eu
    f["geodesic"] = geo
    return f


def lhat_features(B, sub, targets, rng):
    """F7: -log distance from the gene-inferred location (M1, out-of-fold) to each target centroid,
    averaged over LHAT_SAMPLES draws of the inferred location from N(Lhat, Lhat_sd); sd over draws."""
    a = list(B["lhat_names"]).index(sub)
    TI = [list(B["targets"]).index(t) for t in targets]
    L = B["Lhat"][a] + rng.normal(size=(LHAT_SAMPLES, 3)) * B["Lhat_sd"][a]
    d = np.sqrt(((L[:, None] - B["tcen"][TI][None]) ** 2).sum(-1))
    v = -np.log(np.maximum(d, FLOOR))
    return v.mean(0), v.std(0, ddof=1)


def load(full=True):
    """the saved 3D brain; full=True also rebuilds the KD-trees needed by unit_features (reads the cell table)."""
    from scipy.spatial import cKDTree
    z = np.load(BRAIN, allow_pickle=True)
    B = {k: z[k] for k in z.files}
    B["shape"] = tuple(B["shape"])
    B["targets"] = list(B["targets"])
    B["cost_names"] = list(B["cost_names"])
    if full:
        C = read_cells()
        inb = C["cat_names"][C["cat"]] != "unassigned"
        B["tree"] = cKDTree(C["xyz"][inb].astype(float))
        B["sub_inb"] = C["sub"][inb]
        B["nsub"] = len(C["sub_names"])
        M = Mer(C)
        B["ttree"] = []
        rng = np.random.default_rng(SEED + 5)
        for t in B["targets"]:
            ix = M.target_cells(t)
            if len(ix) > 30000:
                ix = rng.choice(ix, 30000, replace=False)
            B["ttree"].append(cKDTree(C["xyz"][ix].astype(float)))
        B["mer"] = M
        B["cells"] = C
        _trees(B)
    return B


def _trees(B):
    from scipy.spatial import cKDTree
    shape = B["shape"]
    tr = np.where(B["tract"])[0]
    cxyz = (np.array(np.unravel_index(tr, shape)).T + 0.5) * VOX
    B["systree"] = [cKDTree(cxyz[B["sys"][tr] == s]) if (B["sys"][tr] == s).sum() > 0 else None
                    for s in B["sys_codes"]]
    B["sttree"] = [cKDTree(cxyz[B["stl"][tr] == s]) for s in B["st_codes"]]


def cloud_cells(B, sub, src):
    """the unit's cloud (mm): d105's selection of the subclass's cells in the source region (src may be 'HY')."""
    ix, level = B["mer"].unit_cells(sub, src)
    return B["cells"]["xyz"][ix].astype(float), level


# ------------------------------------------------------------------------------ build
def build():
    from scipy.sparse.csgraph import dijkstra
    from scipy.spatial import cKDTree
    import d108_prims as P
    t0 = time.time()
    C = read_cells()
    say("cells: %d  [%.0fs]" % (len(C["xyz"]), time.time() - t0))
    onto = ontology_check(C)
    say("ontology check:", onto)
    vb = voxelise(C, onto)
    M = Mer(C)
    tab, S, targets, names = P.world()                         # non-HY sources only (HY rows never read)
    prep = np.load(P.PREP, allow_pickle=True)
    assert list(prep["targets"]) == targets
    # target regions
    rng = np.random.default_rng(SEED + 5)
    shape = vb["shape"]
    mid = np.where(vb["mask"])[0]
    node = -np.ones(len(vb["mask"]), int)
    node[mid] = np.arange(len(mid))
    nsub = len(C["sub_names"])
    sub_all = C["sub"]
    ttree, tcomp, tvox, tcen = [], [], [], []
    for t in targets:
        ix = M.target_cells(t)
        tcen.append(C["xyz"][ix].astype(float).mean(0))
        tcomp.append(np.bincount(sub_all[ix], minlength=nsub) / len(ix))
        nodes = vox_of(vb, C["xyz"][ix].astype(float))
        cnt = np.bincount(nodes, minlength=len(mid))
        tvox.append(np.where(cnt >= 2)[0])
        if len(ix) > 30000:
            ix = rng.choice(ix, 30000, replace=False)
        ttree.append(cKDTree(C["xyz"][ix].astype(float)))
    tcen = np.array(tcen)
    assert np.allclose(tcen, prep["tcen"], atol=1e-3), np.abs(tcen - prep["tcen"]).max()
    tcomp = np.array(tcomp)
    say("target regions: %d (centroids reproduce d108_prep)  [%.0fs]" % (len(targets), time.time() - t0))
    # route fields
    cost_names = ["uniform"] + ["tract%.2f" % c for c in (COST_TRACT,) + COST_SENS]
    fields = np.zeros((len(cost_names), len(targets), len(mid)), np.float32)
    for ci, cn in enumerate(cost_names):
        G = graph(vb, 1.0 if cn == "uniform" else float(cn[5:]))
        for a in range(len(targets)):
            fields[ci, a] = dijkstra(G, directed=False, indices=tvox[a], min_only=True)
    assert np.isfinite(fields).all()
    say("route fields %s  [%.0fs]" % (fields.shape, time.time() - t0))
    # neighbourhood reference PCA (label-free, brain-wide random grey cells)
    inb = vb["inb"]
    sub_inb = C["sub"][inb]
    cat_inb = C["cat_names"][C["cat"]][inb]
    grey = np.where(cat_inb == "grey")[0]
    ref = rng.choice(grey, NB_REF, replace=False)
    R = nb_comp(vb["tree"], sub_inb, nsub, vb["xyz_inb"][ref])
    nb_mu = R.mean(0)
    _, s_, Vt = np.linalg.svd(R - nb_mu, full_matrices=False)
    nb_V = Vt[:NB_PCS].T
    nb_var = float((s_[:NB_PCS] ** 2).sum() / (s_ ** 2).sum())
    say("neighbourhood PCA on %d random grey cells: %d PCs explain %.3f  [%.0fs]" % (NB_REF, NB_PCS, nb_var, time.time() - t0))
    # F4: gene embeddings in D's frame (frozen d105 whole-brain PCA), all Allen subclasses
    z5 = np.load(os.path.join(INTERIM, "d105_models.npz"), allow_pickle=True)
    abc = b5.load_abc()
    keep = z5["wb_keep"]
    Ec = (((abc["X"][:, keep] - z5["wb_mu"]) / z5["wb_sd"]).astype(np.float32) - z5["wb_zc_mean"]) @ z5["wb_loadings"][:, :32]
    subn = abc["subnames"]
    emb_all = np.zeros((len(subn), 32))
    for a in range(len(subn)):
        m = abc["subc"] == a
        emb_all[a] = abc["w"][m] @ Ec[m] / abc["w"][m].sum()
    del abc, Ec
    ni = [subn.index(u) for u in names]
    assert np.allclose(emb_all[ni], z5["wb_emb"][:, :32], atol=1e-3)
    cen_, scl_ = z5["D_centre"], float(z5["D_scale"])
    assert np.allclose((z5["wb_emb"][:, :32] - cen_) / scl_, z5["F_D"], atol=1e-6)
    emb_D = (emb_all - cen_) / scl_
    mer_names = list(C["sub_names"])
    mi = np.array([subn.index(u) if u in subn else -1 for u in mer_names])
    tc_known = tcomp[:, mi >= 0]
    temb = (tc_known / tc_known.sum(1, keepdims=True)) @ emb_all[mi[mi >= 0]]
    temb_D = (temb - cen_) / scl_
    say("F4: target mean embeddings over %d/%d MERFISH subclasses matched to the 10x taxonomy  [%.0fs]" % ((mi >= 0).sum(), len(mi), time.time() - t0))
    # tracts as corridors: systems and structures (>= MIN_TRACT_VOX voxels)
    tr = vb["tract"]
    sys_codes = [int(a) for a in np.unique(vb["sys"][tr]) if not str(C["div_names"][a]).endswith("unassigned")]
    st_codes = [int(a) for a in np.unique(vb["stl"][tr]) if (vb["stl"][tr] == a).sum() >= MIN_TRACT_VOX
                and not str(C["st_names"][a]).endswith("unassigned")]
    B = {"shape": shape, "mask": vb["mask"], "tract": tr, "sys": vb["sys"], "stl": vb["stl"], "sys_codes": np.array(sys_codes),
         "st_codes": np.array(st_codes), "targets": targets, "cost_names": cost_names, "fields": fields, "tcen": tcen,
         "ttree": ttree, "tree": vb["tree"], "sub_inb": sub_inb, "nsub": nsub, "nb_mu": nb_mu, "nb_V": nb_V,
         "tcomp_sqrt": np.sqrt(tcomp), "emb_names": np.array(subn), "emb_D": emb_D, "temb_D": temb_D,
         "sys_names": [str(C["div_names"][a]) for a in sys_codes]}
    _trees(B)
    # target-side tract access (B_s): share of the target region's cells within TRACT_R of tract structure s
    tract_B = np.zeros((len(targets), len(st_codes)))
    for a, t in enumerate(targets):
        ix = M.target_cells(t)
        if len(ix) > 5000:
            ix = rng.choice(ix, 5000, replace=False)
        p = C["xyz"][ix].astype(float)
        tract_B[a] = [float((tt.query(p, workers=4)[0] <= TRACT_R).mean()) for tt in B["sttree"]]
    B["tract_B"] = tract_B
    # units: 27 non-HY sources (atlas cells only) + HY19 in HY (atlas cells only; reusable for the unseen region)
    units = [(k, u) for k in sorted(S) for u in S[k]["units"]] + [("HY", u) for u in b5.HY19]
    B["lhat_names"] = np.array(prep["names"])
    B["Lhat"], B["Lhat_sd"] = prep["Lhat"], prep["Lhat_sd"]
    ucen_prep = {tuple(k.split("|", 1)): v for k, v in zip(prep["unit_keys"], prep["unit_cen"])}
    feats, levels = {}, {}
    clouds_xyz, clouds_off = [], [0]
    for i, (k, u) in enumerate(units):
        ix, level = M.unit_cells(u, k)
        cloud = C["xyz"][ix].astype(float)
        r = np.random.default_rng(SEED + 100 + i)
        f = unit_features(B, cloud, u, targets, r)
        f["distLhat_u"], f["distLhat_u_sd"] = lhat_features(B, u, targets, r)
        if (k, u) in ucen_prep:
            assert np.allclose(f["centroid"], ucen_prep[(k, u)], atol=1e-3), (k, u)
        feats[(k, u)] = f
        levels[(k, u)] = level
        keepc = cloud if len(cloud) <= CLOUD_MAX else cloud[r.choice(len(cloud), CLOUD_MAX, replace=False)]
        clouds_xyz.append(keepc.astype(np.float32))
        clouds_off.append(clouds_off[-1] + len(keepc))
    say("unit clouds and features: %d units (levels %s); non-HY centroids reproduce d108_prep  [%.0fs]"
        % (len(units), {l: list(levels.values()).count(l) for l in set(levels.values())}, time.time() - t0))
    keys = ["%s|%s" % x for x in units]
    out = dict(shape=np.array(shape), voxel_mm=VOX, mask=vb["mask"], tract=tr, tfrac=vb["tfrac"], ncell=vb["ncell"],
               interp=vb["interp"], sys=vb["sys"], stl=vb["stl"], div_names=C["div_names"], st_names=C["st_names"],
               sys_codes=np.array(sys_codes), st_codes=np.array(st_codes), sys_names=np.array(B["sys_names"]),
               st_names_used=np.array([str(C["st_names"][a]) for a in st_codes]),
               targets=np.array(targets), tcen=tcen, tcomp=tcomp, tcomp_sqrt=np.sqrt(tcomp), sub_names=C["sub_names"],
               tract_B=tract_B, cost_names=np.array(cost_names), fields=fields, tvox_counts=np.array([len(v) for v in tvox]),
               nb_mu=nb_mu, nb_V=nb_V, nb_var=nb_var, emb_names=np.array(subn), emb_D=emb_D, temb_D=temb_D,
               lhat_names=B["lhat_names"], Lhat=B["Lhat"], Lhat_sd=B["Lhat_sd"],
               unit_keys=np.array(keys), unit_level=np.array([levels[x] for x in units]),
               cloud_xyz=np.concatenate(clouds_xyz), cloud_off=np.array(clouds_off),
               tags=np.array(json.dumps(TAGS)), provenance_excluded=np.array("P: data/raw/atlas/connectivity and Epi-retro-seq rows never read"))
    names_u = ("centroid", "n_cells", "shape6", "nb6", "tsys6", "tract_access", "nb_comp")
    names_ut = ("cl_med", "cl_p05", "cl_ovl", "nbsim", "gcompat", "corr_shared", "route", "route_sd", "route_gain",
                "euclid", "geodesic", "distLhat_u", "distLhat_u_sd")
    for n in names_u + names_ut:
        out["U_" + n] = np.array([np.asarray(feats[x][n], float) for x in units])
    out["U_route_costs"] = np.array([feats[x]["route_costs"] for x in units])
    np.savez_compressed(BRAIN, **out)
    card = {"id": "d108_3d_brain", "file": os.path.relpath(BRAIN, ROOT), "built": time.strftime("%Y-%m-%d %H:%M:%S"),
            "space": "Allen CCFv3 (mm), folded onto one hemisphere at z = %.1f mm (as d105)" % b5.MIDLINE_Z_MM,
            "sources_read": ["data/raw/abc/merfish_cells.csv", "data/raw/atlas/structure_graph_1.json",
                             "data/raw/abc/precomputed_stats.h5 (F4 only)", "data/interim/d105_models.npz", "data/interim/d108_prep.npz"],
            "excluded_projection_derived": ["data/raw/atlas/connectivity (Allen tracing experiments)", "Epi-retro-seq tables (any row)"],
            "voxels": vb["info"], "route_costs": cost_names, "tract_cost_factor": COST_TRACT, "sensitivity": list(COST_SENS),
            "tract_systems": B["sys_names"], "tract_structures_used": [str(C["st_names"][a]) for a in st_codes],
            "neighbourhood_pca_var": nb_var, "n_units": len(units), "unit_levels": {l: list(levels.values()).count(l) for l in set(levels.values())},
            "tags": TAGS, "meshes_available": "data/raw/atlas/meshes: 20 OBJ meshes (root + 19 hypothalamic structures); no tract meshes, not used",
            "usage": "import d108_3dbrain as B3; B = B3.load(); xyz, level = B3.cloud_cells(B, subclass, region); "
                     "f = B3.unit_features(B, xyz, subclass, targets) (+ B3.lhat_features for F7)",
            "sec": time.time() - t0}
    json.dump(card, open(CARD, "w"), indent=1)
    say("wrote %s and %s  [%.0fs]" % (os.path.relpath(BRAIN, ROOT), os.path.relpath(CARD, ROOT), time.time() - t0))
    figure(out)


def figure(out):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    shape = tuple(out["shape"])
    mask = out["mask"].reshape(shape)
    tract = out["tract"].reshape(shape)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    for a, (axis, ttl) in enumerate(((2, "sagittal projection (x vs y)"), (1, "horizontal projection (x vs z, folded)"), (0, "coronal projection (y vs z, folded)"))):
        ax[a].imshow(mask.sum(axis).T > 0, cmap="Greys", alpha=0.25, origin="upper")
        ax[a].imshow(np.ma.masked_equal(tract.sum(axis).T, 0), cmap="viridis", origin="upper")
        ax[a].set_title("tract voxels (count along the collapsed axis), %s" % ttl, fontsize=8)
        ax[a].set_xticks([]); ax[a].set_yticks([])
    fig.tight_layout()
    p = os.path.join(INTERIM, "d108_3d_brain_views.png")
    fig.savefig(p, dpi=110)
    say("wrote", os.path.relpath(p, ROOT))


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "build":
        build()
    else:
        sys.exit("usage: python3 -I analysis/d108_3dbrain.py build")
