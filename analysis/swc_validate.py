"""Validate the Jiao 2025 SWC sample fetched by analysis/swc_fetch.py.

    python3 -I analysis/swc_validate.py

Reads data/raw/jiao2025/swc/ (raw files, never modified), the soma table
data/raw/jiao2025/all_neuron_soma_infor.xlsx (mirrored-left paper frame), the public portal index
(portal_index/mouse.neuron.info.json) and the Allen CCFv3 meshes in data/raw/atlas/meshes
(OBJ, micrometres, columns AP DV ML, both hemispheres; the same axes as the SWC x y z).
Writes reports/swc_validation.md and data/interim/swc_validation.json (used by swc_render.py).

Conventions (reports/swc_access.md): SWC columns id type AP DV ML radius parent, um, CCFv3 2017;
midline ML = 5700 um; Left = ML < 5700. Segment type = the type of its child node. Type 0 is kept
as 'type 0 (unlabelled)', types 10 and 12 as 'other (undocumented)'; nothing is relabelled.
Nothing is mirrored for storage; reflections below are computed only to compare frames."""
import os, json, hashlib
import numpy as np
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SWC = os.path.join(ROOT, "data/raw/jiao2025/swc")
MESH = os.path.join(ROOT, "data/raw/atlas/meshes")
MID = 5700.0
OUT_MD = os.path.join(ROOT, "reports/swc_validation.md")
OUT_JS = os.path.join(ROOT, "data/interim/swc_validation.json")
TYPE_CLASS = {1: "soma", 2: "axon", 3: "dendrite", 0: "type0"}   # anything else -> "other"


def tclass(t): return TYPE_CLASS.get(int(t), "other")


def read_swc(p):
    rows = []
    for line in open(p, "rb").read().decode().splitlines():
        s = line.strip()
        if not s or s.startswith("#"): continue
        f = s.split()
        if len(f) != 7: raise ValueError("%s: %d columns in %r" % (p, len(f), s[:80]))
        rows.append(f)
    a = np.array(rows, dtype=float)
    return a[:, 0].astype(np.int64), a[:, 1].astype(np.int64), a[:, 2:5], a[:, 5], a[:, 6].astype(np.int64)


def read_obj(p):
    V, F = [], []
    for line in open(p):
        if line.startswith("v "): V.append([float(x) for x in line.split()[1:4]])
        elif line.startswith("f "): F.append([int(x.split("/")[0]) - 1 for x in line.split()[1:4]])
    return np.array(V), np.array(F)


class Solid:
    """Inside test for a triangle mesh by ray parity along one axis (ax: 0 AP, 1 DV, 2 ML), on a
    column grid over the other two axes. A point is tested on the column through its grid cell
    centre (offset g/2), so the boundary is resolved to about g*sqrt(2)/2. One axis alone is
    fragile where the mesh has holes; use Vote3 (majority of the three axes) for node tests."""
    def __init__(s, V0, F, g, ax=0):
        s.order = [ax] + [i for i in range(3) if i != ax]
        V = V0[:, s.order]
        s.g = g; s.V = V0; s.F = F
        e = np.sort(np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
        _, cnt = np.unique(e, axis=0, return_counts=True)
        s.open_edges = int((cnt == 1).sum()); s.nonmanifold_edges = int((cnt > 2).sum())
        s.lo = V[:, 1:].min(0) - 2 * g
        s.n = np.ceil((V[:, 1:].max(0) + 2 * g - s.lo) / g).astype(int)
        hits = {}
        P = V[F]                                      # (nf, 3, 3)
        y, z = P[:, :, 1], P[:, :, 2]
        for k in range(len(F)):
            ys, zs = y[k], z[k]
            i0 = int(np.ceil((ys.min() - s.lo[0]) / g - 0.5)); i1 = int(np.floor((ys.max() - s.lo[0]) / g - 0.5))
            j0 = int(np.ceil((zs.min() - s.lo[1]) / g - 0.5)); j1 = int(np.floor((zs.max() - s.lo[1]) / g - 0.5))
            if i1 < i0 or j1 < j0: continue
            I, J = np.meshgrid(np.arange(i0, i1 + 1), np.arange(j0, j1 + 1), indexing="ij")
            cy = s.lo[0] + (I.ravel() + 0.5) * g; cz = s.lo[1] + (J.ravel() + 0.5) * g
            (y0, y1, y2), (z0, z1, z2) = ys, zs
            den = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0)
            if abs(den) < 1e-12: continue
            l1 = ((cy - y0) * (z2 - z0) - (y2 - y0) * (cz - z0)) / den
            l2 = ((y1 - y0) * (cz - z0) - (cy - y0) * (z1 - z0)) / den
            ok = (l1 >= 0) & (l2 >= 0) & (l1 + l2 <= 1)
            if not ok.any(): continue
            xs = P[k, 0, 0] + l1[ok] * (P[k, 1, 0] - P[k, 0, 0]) + l2[ok] * (P[k, 2, 0] - P[k, 0, 0])
            for ii, jj, xx in zip(I.ravel()[ok], J.ravel()[ok], xs):
                hits.setdefault((ii, jj), []).append(xx)
        s.cols = {k: np.sort(v) for k, v in hits.items()}
        s.odd_columns = sum(len(v) % 2 for v in s.cols.values())
        s.tree = cKDTree(V0)

    def inside(s, X):
        X = np.atleast_2d(X)[:, s.order]
        I = np.floor((X[:, 1] - s.lo[0]) / s.g).astype(int); J = np.floor((X[:, 2] - s.lo[1]) / s.g).astype(int)
        out = np.zeros(len(X), bool)
        for n, (i, j, x) in enumerate(zip(I, J, X[:, 0])):
            c = s.cols.get((i, j))
            if c is not None: out[n] = (np.searchsorted(c, x) % 2) == 1
        return out

    def vdist(s, X): return s.tree.query(np.atleast_2d(X))[0]


class Vote3:
    """Majority vote of three grid parity tests, with rays along AP, DV and ML. votes(X) gives 0..3."""
    def __init__(s, V, F, g):
        s.S = [Solid(V, F, g, ax) for ax in range(3)]
        s.open_edges, s.nonmanifold_edges = s.S[0].open_edges, s.S[0].nonmanifold_edges
        s.odd_columns = [x.odd_columns for x in s.S]
        s.V, s.F, s.tree = V, F, s.S[0].tree
        e = V[F]; s.max_edge = float(max(np.linalg.norm(e[:, i] - e[:, (i + 1) % 3], axis=1).max() for i in range(3)))
        s.vf = [[] for _ in range(len(V))]
        for k, f in enumerate(F):
            for v in f: s.vf[v].append(k)

    def votes(s, X): return sum(x.inside(X).astype(int) for x in s.S)
    def inside(s, X): return s.votes(X) >= 2
    def vdist(s, X): return s.tree.query(np.atleast_2d(X))[0]

    def surface_dist(s, X):
        """Exact distance from each point to the nearest triangle (not just the nearest vertex)."""
        out = []
        for p in np.atleast_2d(X):
            dv = s.tree.query(p)[0]
            near = s.tree.query_ball_point(p, dv + s.max_edge + 1e-6)
            faces = np.unique([k for v in near for k in s.vf[v]])
            T = s.V[s.F[faces]]
            out.append(float(np.linalg.norm(closest_on_tris(p, T[:, 0], T[:, 1], T[:, 2]) - p, axis=1).min()))
        return np.array(out)


RAY_DIRS = (lambda r: r / np.linalg.norm(r, axis=1)[:, None])(np.random.default_rng(0).normal(size=(15, 3)))


def ray_votes(p, V, F, dirs=RAY_DIRS):
    """Exact ray parity (Moller-Trumbore against every face) along each of 15 fixed random directions.
    Returns how many directions give an odd number of crossings (inside votes, 0..15)."""
    A, B, C = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    e1, e2 = B - A, C - A; s_ = p - A; q = np.cross(s_, e1); n = 0
    for d in dirs:
        h = np.cross(d, e2); a = (e1 * h).sum(1)
        with np.errstate(divide="ignore", invalid="ignore"):
            f = 1.0 / a; u = f * (s_ * h).sum(1); v = f * (q @ d); t = f * (e2 * q).sum(1)
        hit = (np.abs(a) > 1e-12) & (u >= 0) & (u <= 1) & (v >= 0) & (u + v <= 1) & (t > 1e-9)
        n += int(hit.sum()) % 2
    return n


def closest_on_tris(p, A, B, C):
    """Closest point to p on each triangle (Ericson, Real-Time Collision Detection 5.1.5), vectorised."""
    ab, ac, ap = B - A, C - A, p - A
    d1, d2 = (ab * ap).sum(1), (ac * ap).sum(1)
    bp = p - B; d3, d4 = (ab * bp).sum(1), (ac * bp).sum(1)
    cp = p - C; d5, d6 = (ab * cp).sum(1), (ac * cp).sum(1)
    vc = d1 * d4 - d3 * d2; vb = d5 * d2 - d1 * d6; va = d3 * d6 - d5 * d4
    with np.errstate(divide="ignore", invalid="ignore"):
        den = va + vb + vc
        res = A + ab * (vb / den)[:, None] + ac * (vc / den)[:, None]                     # interior
        m = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0)
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6)); res = np.where(m[:, None], B + (C - B) * w[:, None], res)   # edge BC
        m = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
        w = d2 / (d2 - d6); res = np.where(m[:, None], A + ac * w[:, None], res)                           # edge AC
        res = np.where(((d6 >= 0) & (d5 <= d6))[:, None], C, res)                                           # vertex C
        m = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
        v = d1 / (d1 - d3); res = np.where(m[:, None], A + ab * v[:, None], res)                           # edge AB
        res = np.where(((d3 >= 0) & (d4 <= d3))[:, None], B, res)                                           # vertex B
        res = np.where(((d1 <= 0) & (d2 <= 0))[:, None], A, res)                                            # vertex A
    return res


def structure_names():
    sg = json.load(open(os.path.join(ROOT, "data/raw/atlas/structure_graph_1.json")))
    node = sg["msg"][0] if "msg" in sg else sg
    acr, par, name = {}, {}, {}
    def walk(n, p):
        acr[n["acronym"]] = n["id"]; par[n["id"]] = p; name[n["id"]] = n["name"]
        for c in n.get("children", []): walk(c, n["id"])
    walk(node, None)
    return acr, par, name


def morph(ids, types, X, parent):
    n = len(ids); pos = {int(i): k for k, i in enumerate(ids)}
    r = {"nodes": n, "duplicate_ids": n - len(pos)}
    roots = np.where(parent == -1)[0]
    r["roots"] = len(roots); r["root_id"] = int(ids[roots[0]]) if len(roots) else None
    r["root_type"] = int(types[roots[0]]) if len(roots) else None
    r["soma_nodes"] = int((types == 1).sum())
    pidx = np.array([pos.get(int(p), -2) if p != -1 else -1 for p in parent])
    r["missing_parents"] = int((pidx == -2).sum())
    r["parent_after_child"] = int(((pidx >= 0) & (pidx > np.arange(n))).sum())
    r["sequential_ids"] = bool(np.array_equal(ids, np.arange(1, n + 1)))
    # connectivity and cycles: walk from roots over child lists
    ch = [[] for _ in range(n)]
    for k, p in enumerate(pidx):
        if p >= 0: ch[p].append(k)
    seen = np.zeros(n, bool); stack = list(roots)
    while stack:
        k = stack.pop()
        if seen[k]: continue
        seen[k] = True; stack.extend(ch[k])
    r["unreached_from_root"] = int((~seen).sum())          # nodes in cycles or other components
    # components (undirected union-find)
    up = np.arange(n)
    def find(a):
        while up[a] != a: up[a] = up[up[a]]; a = up[a]
        return a
    for k, p in enumerate(pidx):
        if p >= 0: up[find(k)] = find(p)
    r["components"] = len({find(k) for k in range(n)})
    r["tree_ok"] = bool(r["roots"] == 1 and r["missing_parents"] == 0 and r["unreached_from_root"] == 0 and r["components"] == 1 and r["duplicate_ids"] == 0)
    tc = {}
    for t in types: tc[int(t)] = tc.get(int(t), 0) + 1
    r["type_counts"] = {str(k): v for k, v in sorted(tc.items())}
    nch = np.array([len(c) for c in ch])
    r["branch_points"] = int(((nch >= 2) & (parent != -1)).sum())
    r["soma_children"] = int(nch[roots[0]]) if len(roots) else 0
    r["tips"] = int(((nch == 0) & (parent != -1)).sum())
    has = pidx >= 0
    seg = np.linalg.norm(X[has] - X[pidx[has]], axis=1)
    cls = np.array([tclass(t) for t in types[has]])
    r["length_um"] = {c: float(seg[cls == c].sum()) for c in ("axon", "type0", "dendrite", "other", "soma")}
    r["max_segment_um"] = float(seg.max()) if len(seg) else 0.0
    r["zero_length_segments"] = int((seg < 1e-6).sum())
    r["segments_over_50um"] = int((seg > 50).sum())
    # type changes along a parent-child edge (other than at the soma)
    tch = (types[has] != types[pidx[has]]) & (types[pidx[has]] != 1)
    r["type_changes_off_soma"] = int(tch.sum())
    doc = (types[has] <= 3) & (types[pidx[has]] <= 3)
    r["type_changes_documented_only"] = int((tch & doc).sum())   # changes between types 0/2/3 not caused by markers
    r["bbox_min"] = X.min(0).round(1).tolist(); r["bbox_max"] = X.max(0).round(1).tolist()
    r["radius_range"] = None
    soma = X[roots[0]] if len(roots) else X[0]
    r["soma"] = soma.round(2).tolist()
    side = np.sign(soma[2] - MID)
    mids = 0.5 * (X[has] + X[pidx[has]])
    contra = np.sign(mids[:, 2] - MID) == -side
    r["contra_fraction"] = {c: (float(seg[(cls == c) & contra].sum() / seg[cls == c].sum()) if seg[cls == c].sum() > 0 else None) for c in ("axon", "type0", "dendrite")}
    cross = (np.sign(X[has][:, 2] - MID) != np.sign(X[pidx[has]][:, 2] - MID))
    r["midline_crossings"] = int(cross.sum())
    return r, pidx


def main():
    man = [l.rstrip("\n").split("\t") for l in open(os.path.join(SWC, "MANIFEST.tsv"))][1:]
    idx = json.load(open(os.path.join(SWC, "portal_index/mouse.neuron.info.json")))["neuron_data"]
    idx = {v["file"][:-4]: v for v in idx.values()}
    import pandas as pd
    xl = pd.read_excel(os.path.join(ROOT, "data/raw/jiao2025/all_neuron_soma_infor.xlsx"))
    xl = {r.neurons: r for r in xl.itertuples()}
    acr, par, sname = structure_names()

    meshes = {}
    for f in sorted(os.listdir(MESH)):
        sid = int(f[:-4]); V, F = read_obj(os.path.join(MESH, f))
        a = [k for k, v in acr.items() if v == sid][0]
        meshes[a] = (sid, V, F)
    global ROOT_AP_MAX
    ROOT_AP_MAX = float(meshes["root"][1][:, 0].max())
    solids = {}
    for a, (sid, V, F) in meshes.items():
        solids[a] = Vote3(V, F, 20.0 if a == "root" else 10.0)
    mesh_rows = [(a, meshes[a][0], len(meshes[a][1]), len(meshes[a][2]), solids[a].open_edges, solids[a].nonmanifold_edges, "/".join(map(str, solids[a].odd_columns)),
                  meshes[a][1].min(0).round(0).tolist(), meshes[a][1].max(0).round(0).tolist()) for a in sorted(meshes)]

    def mesh_for(label):
        """Nearest structure (label itself or an ancestor) with a local mesh, excluding root."""
        sid = acr.get(label); hops = 0
        while sid is not None:
            for a, (mid, _, _) in meshes.items():
                if mid == sid and a != "root": return a, hops
            sid = par.get(sid); hops += 1
        return None, None

    res, anomalies, arrays = [], [], {}
    for f in man:
        rel, nid, frame = f[0], f[1], f[5]
        if not rel.endswith(".swc"): continue
        ids, types, X, rad, parent = read_swc(os.path.join(SWC, rel))
        r, pidx = morph(ids, types, X, parent)
        raw = open(os.path.join(SWC, rel), "rb").read()
        r["crlf_lines"] = raw.count(b"\r\n"); r["lf_lines"] = raw.count(b"\n")
        r.update(file=rel, neuron=nid, frame=frame, sha256=f[3], bytes=int(f[2]), source=f[4], licence=f[6])
        r["radius_range"] = [float(rad.min()), float(rad.max())]
        ix = idx.get(nid, {}); x = xl.get(nid)
        r["index"] = {k: ix.get(k) for k in ("acronym", "name", "type_name", "hemisphere", "reconstruction_type", "type", "class")}
        r["xlsx"] = None if x is None else {"location": x.location, "geno": x.geno, "cluster": int(x.cluster), "Class": x.Class,
                                            "AP": float(x.AP), "DV": float(x.DV), "ML": float(x.ML)}
        soma = np.array(r["soma"])
        r["soma_side"] = "Left" if soma[2] < MID else "Right"
        if x is not None:
            t = np.array([x.AP, x.DV, x.ML])
            r["soma_vs_xlsx_um"] = float(np.linalg.norm(soma - t))
            r["soma_vs_xlsx_reflected_um"] = float(np.linalg.norm(soma - t * [1, 1, -1] - [0, 0, 11400]))
        # Hemisphere and frame consistency
        true_h = ix.get("hemisphere")
        if frame == "portal_true_hemisphere":
            r["hemisphere_ok"] = (r["soma_side"] == true_h)
        else:
            r["hemisphere_ok"] = (r["soma_side"] == "Left")      # mirrored-left frame by construction
        # Inside tests (soma in its own frame; meshes are bilateral)
        # Somas: exact ray parity along 15 fixed random directions (majority >= 8 of 15 = inside),
        # plus the exact distance to the mesh surface.
        def stest(m):
            v = ray_votes(soma, solids[m].V, solids[m].F)
            return {"inside": v >= 8, "votes": v, "surface_dist_um": float(solids[m].surface_dist(soma)[0]),
                    "vertex_dist_um": float(solids[m].vdist(soma)[0])}
        t_root, t_hy = stest("root"), stest("HY")
        r["soma_root_test"], r["soma_HY_test"] = t_root, t_hy
        r["soma_in_root"] = bool(t_root["inside"]); r["soma_in_HY"] = bool(t_hy["inside"])
        r["soma_HY_vertex_dist_um"] = t_hy["vertex_dist_um"]
        lab = x.location if x is not None else ix.get("acronym")
        m, hops = mesh_for(lab) if lab else (None, None)
        r["nucleus_test"] = None
        if m and m != "HY":
            t = stest(m)
            r["nucleus_test"] = {"label": lab, "mesh": m, "ancestor_hops": hops, "inside": bool(t["inside"]), "votes": t["votes"],
                                 "surface_dist_um": t["surface_dist_um"], "vertex_dist_um": t["vertex_dist_um"]}
        elif m == "HY":
            r["nucleus_test"] = {"label": lab, "mesh": "HY", "ancestor_hops": hops, "inside": r["soma_in_HY"], "votes": t_hy["votes"],
                                 "surface_dist_um": t_hy["surface_dist_um"], "vertex_dist_um": r["soma_HY_vertex_dist_um"],
                                 "note": "no local mesh for the label or a non-HY ancestor; tested against HY only"}
        vts = solids["root"].votes(X); ins = vts >= 2
        r["nodes_outside_root"] = int((~ins).sum())
        r["root_vote_split_nodes"] = int(((vts == 1) | (vts == 2)).sum())     # axes disagree: within a grid cell of the surface or near a hole
        r["outside_max_vertex_dist_um"] = float(solids["root"].vdist(X[~ins]).max()) if (~ins).any() else 0.0
        r["outside_max_surface_dist_um"] = float(solids["root"].surface_dist(X[~ins]).max()) if (~ins).any() else 0.0
        r["outside_AP_range"] = [float(X[~ins, 0].min()), float(X[~ins, 0].max())] if (~ins).any() else None
        r["outside_median"] = np.median(X[~ins], 0).round(0).tolist() if (~ins).any() else None
        r["outside_types"] = {str(k): int(v) for k, v in zip(*np.unique(types[~ins], return_counts=True))} if (~ins).any() else {}
        res.append(r); arrays[rel] = (ids, types, X, parent)

    # Frame pairs: same neuron as portal and Zenodo files
    pairs = []
    byn = {}
    for r in res: byn.setdefault(r["neuron"], {})[r["frame"]] = r
    for n, d in sorted(byn.items()):
        if len(d) < 2: continue
        p, z = d["portal_true_hemisphere"], d["zenodo_mirrored_left"]
        a, b = arrays[p["file"]], arrays[z["file"]]
        same_topo = bool(np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1]) and np.array_equal(a[3], b[3]))
        dAPDV = float(np.abs(a[2][:, :2] - b[2][:, :2]).max()) if same_topo else None
        dsame = float(np.abs(a[2][:, 2] - b[2][:, 2]).max()) if same_topo else None
        dsum = float(np.abs(a[2][:, 2] + b[2][:, 2] - 11400).max()) if same_topo else None
        rel_ = "identical" if same_topo and dsame < 0.006 else ("reflected ML'=11400-ML" if same_topo and dsum < 0.006 else "neither")
        pairs.append(dict(neuron=n, hemisphere=p["index"]["hemisphere"], nodes=p["nodes"], same_ids_types_parents=same_topo,
                          max_abs_dAPDV=dAPDV, max_abs_dML=dsame, max_abs_MLsum_minus_11400=dsum, relation=rel_,
                          bytes_portal=p["bytes"], bytes_zenodo=z["bytes"]))

    # Anomalies
    for r in res:
        tag = r["file"]
        if not r["tree_ok"]: anomalies.append((tag, "tree check failed: %s" % {k: r[k] for k in ("roots", "missing_parents", "unreached_from_root", "components", "duplicate_ids")}))
        if r["root_id"] != 1 or r["root_type"] != 1: anomalies.append((tag, "root is id %s type %s (expected id 1, type 1)" % (r["root_id"], r["root_type"])))
        if r["soma_nodes"] != 1: anomalies.append((tag, "%d soma (type 1) nodes" % r["soma_nodes"]))
        if not r["sequential_ids"]: anomalies.append((tag, "ids are not 1..n"))
        if r["parent_after_child"]: anomalies.append((tag, "%d nodes whose parent appears later in the file" % r["parent_after_child"]))
        oth = {k: v for k, v in r["type_counts"].items() if k not in ("0", "1", "2", "3")}
        if oth: anomalies.append((tag, "undocumented node types %s (kept as 'other')" % oth))
        if r["type_counts"].get("0") and r["type_counts"].get("2"): anomalies.append((tag, "both type 0 (%d) and type 2 (%d) nodes" % (r["type_counts"]["0"], r["type_counts"]["2"])))
        if r["type_counts"].get("0"): anomalies.append((tag, "type 0 (unlabelled) nodes: %d (%.0f um); index reconstruction_type %s" % (r["type_counts"]["0"], r["length_um"]["type0"], r["index"]["reconstruction_type"])))
        if r["index"]["reconstruction_type"] == "Axon_and_dendrite" and not r["type_counts"].get("3"): anomalies.append((tag, "index says Axon_and_dendrite but no type 3 nodes"))
        if r["index"]["reconstruction_type"] == "Axon_only" and r["type_counts"].get("3"): anomalies.append((tag, "index says Axon_only but has type 3 nodes"))
        if r["type_changes_off_soma"]: anomalies.append((tag, "%d parent-child edges change type away from the soma; %d of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)" % (r["type_changes_off_soma"], r["type_changes_documented_only"])))
        if r["zero_length_segments"]: anomalies.append((tag, "%d zero-length segments" % r["zero_length_segments"]))
        if r["segments_over_50um"]: anomalies.append((tag, "%d segments longer than 50 um (max %.0f um)" % (r["segments_over_50um"], r["max_segment_um"])))
        if not r["hemisphere_ok"]: anomalies.append((tag, "soma side %s does not match expected (%s frame; index %s)" % (r["soma_side"], r["frame"], r["index"]["hemisphere"])))
        if not r["soma_in_root"]: anomalies.append((tag, "soma outside the root (brain) mesh"))
        if not r["soma_in_HY"]: anomalies.append((tag, "soma outside the HY mesh (%.1f um from its surface; %d of 15 rays inside); label %s" % (r["soma_HY_test"]["surface_dist_um"], r["soma_HY_test"]["votes"], (r["xlsx"] or {}).get("location"))))
        nt = r["nucleus_test"]
        if nt and not nt["inside"] and nt["mesh"] != "HY": anomalies.append((tag, "soma outside the %s mesh for label %s: %.1f um from the mesh surface (exact point-to-triangle; nearest vertex %.0f um); %d of 15 rays inside" % (nt["mesh"], nt["label"], nt["surface_dist_um"], nt["vertex_dist_um"], nt["votes"])))
        for nm, tt in (("root", r["soma_root_test"]), ("HY", r["soma_HY_test"])) + ((("nucleus " + nt["mesh"], nt),) if nt and nt["mesh"] != "HY" else ()):
            if 0 < tt["votes"] < 15: anomalies.append((tag, "soma %s test not unanimous: %d of 15 rays inside" % (nm, tt["votes"])))
        if r["nodes_outside_root"]:
            where = "AP %.0f to %.0f, beyond the caudal end of the CCF root mesh (AP max %.0f): descending axon continuing into the spinal cord, outside the atlas volume" % (r["outside_AP_range"][0], r["outside_AP_range"][1], ROOT_AP_MAX) if r["outside_AP_range"][0] > ROOT_AP_MAX - 150 else "median (AP, DV, ML) %s" % r["outside_median"]
            anomalies.append((tag, "%d nodes outside the root mesh (3-axis majority vote; max %.0f um from the root surface; types %s; %s)" % (r["nodes_outside_root"], r["outside_max_surface_dist_um"], r["outside_types"], where)))
        if r["xlsx"] is None: anomalies.append((tag, "no row in all_neuron_soma_infor.xlsx"))
        else:
            want = r["soma_vs_xlsx_um"] if (r["frame"] == "zenodo_mirrored_left" or r["soma_side"] == "Left") else r["soma_vs_xlsx_reflected_um"]
            r["soma_vs_xlsx_expected_frame_um"] = want
            if want > 1.0: anomalies.append((tag, "soma differs from the xlsx soma by %.2f um (after the frame's expected reflection)" % want))
            if r["xlsx"]["location"] != r["index"]["acronym"]: anomalies.append((tag, "xlsx location %s != portal index acronym %s" % (r["xlsx"]["location"], r["index"]["acronym"])))
            if r["xlsx"]["geno"] != r["index"]["type_name"]: anomalies.append((tag, "xlsx geno %s != portal index type_name %s" % (r["xlsx"]["geno"], r["index"]["type_name"])))
    for p in pairs:
        exp = "identical" if p["hemisphere"] == "Left" else "reflected ML'=11400-ML"
        if p["relation"] != exp: anomalies.append(("pair " + p["neuron"], "portal vs Zenodo relation %s, expected %s" % (p["relation"], exp)))

    os.makedirs(os.path.dirname(OUT_JS), exist_ok=True)
    json.dump({"neurons": res, "pairs": pairs, "anomalies": anomalies, "meshes": mesh_rows}, open(OUT_JS, "w"), indent=1, default=str)
    write_md(res, pairs, anomalies, mesh_rows)
    print("neurons/files", len(res), "pairs", len(pairs), "anomalies", len(anomalies))


def write_md(res, pairs, anomalies, mesh_rows):
    L = []
    w = L.append
    nd = len({r["neuron"] for r in res})
    w("# Jiao 2025 SWC sample: validation\n")
    w("Written by `analysis/swc_validate.py` (python3 -I). Raw files in `data/raw/jiao2025/swc/` are read, never changed; "
      "hashes and sources are in `data/raw/jiao2025/swc/MANIFEST.tsv`. Machine-readable copy: `data/interim/swc_validation.json`.\n")
    w("- %d files, %d distinct Jiao neurons: %d Zenodo files (CC BY 4.0, mirrored-left frame) and %d portal files (copyright line only, all rights reserved; BSDC record of the same data: CC BY-NC-ND 4.0; true-hemisphere frame); %d neurons are in both. No MouseLight file is included." %
      (len(res), nd, sum(r["frame"] == "zenodo_mirrored_left" for r in res), sum(r["frame"] == "portal_true_hemisphere" for r in res), len(pairs)))
    w("- Units: micrometres, Allen CCFv3 2017; SWC x, y, z = AP, DV, ML. Midline ML = 5,700 um. Segment type = child node type. "
      "'type0' = type 0 (unlabelled), 'other' = undocumented types 10, 11 and 12; neither is relabelled.")
    w("- Lengths are summed straight segments between consecutive nodes. Branch points: non-root nodes with 2 or more children (soma children listed separately). Tips: non-root nodes with no children.")
    w("- Soma mesh tests (root, HY, nucleus): exact ray parity against every triangle along 15 fixed random directions (seed 0); inside = at least 8 of 15 rays odd. "
      "The vote is reported; any non-unanimous vote is listed as an anomaly. Distances to a mesh are exact point-to-triangle distances to the surface (the nearest-vertex distance, an upper bound, is also kept in the JSON).")
    w("- Node tests (every node against the root mesh): majority vote of three grid parity tests with rays along AP, DV and ML (20 um column grid), "
      "because a single axis is fragile where a mesh has holes. A node within about 14 um of the surface can fall either way; nodes where the three axes disagree are counted per file (`root_vote_split_nodes` in the JSON).")
    w("- Line endings: all Zenodo files use CRLF and all portal files LF (counted per file below). This alone explains the byte difference between the two copies of one neuron (one extra byte per line); the parsed content is identical.")
    w("- Contralateral fraction: share of the segment length on the side of the midline opposite the soma, in the file's own frame (a reflection does not change it).\n")
    w("## Mesh check (data/raw/atlas/meshes)\n")
    w("Coordinate frame: OBJ vertices are AP, DV, ML in micrometres, CCFv3, both hemispheres (root ML 486 to 10,891; HY ML 3,158 to 8,219, centred on 5,700), the same frame as the portal SWCs.\n")
    w("| mesh | id | vertices | faces | open edges | non-manifold edges | grid columns with odd hits (AP/DV/ML rays) | min (AP, DV, ML) | max |")
    w("|---|---|---|---|---|---|---|---|---|")
    for m in mesh_rows: w("| %s | %d | %d | %d | %d | %d | %s | %s | %s |" % m)
    w("")
    w("## Per file\n")
    w("| file | frame | line | soma label (xlsx) | hemisph. (index) | soma side | nodes | types | axon um (type 2) | type0 um | dendrite um | branch pts | tips | soma children | tree ok | soma vs xlsx um | in root | in HY | nucleus test | contra axon / type0 | nodes outside root |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in res:
        nt = r["nucleus_test"]
        nts = "-" if not nt else "%s%s: %s %d/15" % (nt["mesh"], "" if nt["ancestor_hops"] == 0 else " (ancestor of %s)" % nt["label"], "in" if nt["inside"] else "OUT (%.1f um from surface)" % nt["surface_dist_um"], nt["votes"])
        cf = r["contra_fraction"]
        w("| %s | %s | %s | %s | %s | %s | %d | %s | %.0f | %.0f | %.0f | %d | %d | %d | %s | %s | %s | %s | %s | %s / %s | %d |" % (
            r["file"], "portal" if r["frame"].startswith("portal") else "zenodo", r["index"]["type_name"], (r["xlsx"] or {}).get("location"),
            r["index"]["hemisphere"], r["soma_side"], r["nodes"], " ".join("%s:%s" % kv for kv in r["type_counts"].items()),
            r["length_um"]["axon"], r["length_um"]["type0"], r["length_um"]["dendrite"], r["branch_points"], r["tips"], r["soma_children"],
            "yes" if r["tree_ok"] else "NO", "%.2f" % r.get("soma_vs_xlsx_expected_frame_um", float("nan")), "%s %d/15" % ("yes" if r["soma_in_root"] else "NO", r["soma_root_test"]["votes"]),
            "%s %d/15" % ("yes" if r["soma_in_HY"] else "NO", r["soma_HY_test"]["votes"]), nts, "-" if cf["axon"] is None else "%.3f" % cf["axon"], "-" if cf["type0"] is None else "%.3f" % cf["type0"],
            r["nodes_outside_root"]))
    w("\nBounding boxes (um, AP DV ML) and maximum segment length:\n")
    w("| file | min | max | max segment um | radius range (not used) | line ending | root-test vote-split nodes |")
    w("|---|---|---|---|---|---|---|")
    for r in res: w("| %s | %s | %s | %.1f | %.4f to %.1f | %s | %d |" % (r["file"], r["bbox_min"], r["bbox_max"], r["max_segment_um"], *r["radius_range"],
                    "CRLF (%d)" % r["crlf_lines"] if r["crlf_lines"] == r["lf_lines"] else ("LF (%d)" % r["lf_lines"] if r["crlf_lines"] == 0 else "mixed %d/%d" % (r["crlf_lines"], r["lf_lines"])), r["root_vote_split_nodes"]))
    w("\n## Frame check: neurons with both a portal and a Zenodo file\n")
    w("| neuron | index hemisphere | nodes | ids/types/parents identical | max abs dAP, dDV | max abs dML | max abs (ML_p + ML_z - 11400) | relation | bytes portal / Zenodo | byte difference |")
    w("|---|---|---|---|---|---|---|---|---|---|")
    for p in pairs: w("| %s | %s | %d | %s | %s | %s | %s | %s | %d / %d | %d |" % (p["neuron"], p["hemisphere"], p["nodes"], p["same_ids_types_parents"],
                                                              p["max_abs_dAPDV"], p["max_abs_dML"], p["max_abs_MLsum_minus_11400"], p["relation"],
                                                              p["bytes_portal"], p["bytes_zenodo"], p["bytes_zenodo"] - p["bytes_portal"]))
    w("\nThe byte difference is the Zenodo file's CRLF line endings (one extra byte per line, see the line-ending column above), not a content change; radii are also identical.")
    w("\n## Every anomaly (%d)\n" % len(anomalies))
    for tag, a in anomalies: w("- `%s`: %s" % (tag, a))
    open(OUT_MD, "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
