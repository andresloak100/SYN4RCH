"""Canonical axon cache for the ION hippocampus set. Same treatment for every neuron, train or test.

Per neuron: parse SWC from the zip (no extraction), order topologically from the single soma root, mirror Right
hemisphere neurons to the left (ML -> 11400 - ML, same rule as the hypothalamic adapter), then simplify along the
arbor: keep root, branch points, tips, and any node at >= 40 um path length from the last kept ancestor; finally
subdivide any edge > 50 um. Every kept point lies on the original arbor. Output: xyz_um (AP,DV,ML), parent, types.
"""
from __future__ import annotations

import csv
import json
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ZIP = ROOT / 'data/raw/hippocampus_neuroxiv/ION_Hipp_full_SWC_CCFv3.zip'
INV = ROOT / 'data/interim/hippocampus/inventory.csv'
CACHE = ROOT / 'data/interim/hippocampus/cache'
KEEP_ARC_UM = 40.0
MAX_EDGE_UM = 50.0
ML_MAX = 11400.0


def build(args):
    nid, hemi = args
    out = CACHE / f'{nid}.npz'
    if out.exists():
        return nid, None
    with zipfile.ZipFile(ZIP) as z:
        text = z.read(f'ION_Hipp_full_SWC_CCFv3/{nid}.swc').decode()
    a = np.loadtxt([ln for ln in text.splitlines() if ln.strip() and not ln.startswith('#')], ndmin=2)
    ids = a[:, 0].astype(np.int64)
    xyz = a[:, 2:5].copy()
    pos = {int(i): k for k, i in enumerate(ids)}
    par = np.array([-1 if p == -1 else pos[int(p)] for p in a[:, 6]], np.int64)
    if hemi == 'Right':
        xyz[:, 2] = ML_MAX - xyz[:, 2]
    n = len(par)
    root = int(np.flatnonzero(par == -1)[0])
    kids = [[] for _ in range(n)]
    for c in range(n):
        if par[c] >= 0:
            kids[par[c]].append(c)
    # DFS from root: decide kept nodes, link each kept node to its last kept ancestor
    new_xyz, new_par = [xyz[root]], [-1]
    stack = [(root, 0, 0.0)]          # (node, index of last kept ancestor in new arrays, arc since it)
    while stack:
        node, last, arc = stack.pop()
        for c in kids[node]:
            d = arc + float(np.linalg.norm(xyz[c] - xyz[node]))
            nk = len(kids[c])
            if nk != 1 or d >= KEEP_ARC_UM:
                p0 = new_xyz[last]
                L = float(np.linalg.norm(xyz[c] - p0))
                steps = max(1, int(np.ceil(L / MAX_EDGE_UM - 1e-9)))
                prev = last
                for s in range(1, steps):
                    new_xyz.append(p0 + (xyz[c] - p0) * s / steps)
                    new_par.append(prev)
                    prev = len(new_xyz) - 1
                new_xyz.append(xyz[c])
                new_par.append(prev)
                stack.append((c, len(new_xyz) - 1, 0.0))
            else:
                stack.append((c, last, d))
    X = np.asarray(new_xyz)
    P = np.asarray(new_par, np.int64)
    T = np.full(len(P), 2, np.int64)
    T[0] = 1
    np.savez(out, xyz_um=X, parent=P, types=T)
    child = np.flatnonzero(par >= 0)
    orig_cable = float(np.linalg.norm(xyz[child] - xyz[par[child]], axis=1).sum())
    new_cable = float(np.linalg.norm(X[1:] - X[P[1:]], axis=1).sum())
    return nid, {'orig_nodes': n, 'nodes': len(P), 'orig_cable_um': orig_cable, 'cable_um': new_cable}


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(INV)))
    jobs = [(r['neuron_id'], r['hemisphere']) for r in rows]
    stats = {}
    with ProcessPoolExecutor(6) as ex:
        for k, (nid, s) in enumerate(ex.map(build, jobs, chunksize=8), 1):
            if s:
                stats[nid] = s
            if k % 1000 == 0:
                print(f'  {k}/{len(jobs)}', flush=True)
    p = CACHE.parent / 'cache_stats.json'
    old = json.loads(p.read_text()) if p.exists() else {}
    old.update(stats)
    p.write_text(json.dumps(old))
    r = np.array([v['cable_um'] / v['orig_cable_um'] for v in old.values() if v['orig_cable_um'] > 0])
    nodes = np.array([v['nodes'] for v in old.values()])
    print(len(old), 'cached; cable kept ratio p1/p50/min', np.quantile(r, .01), np.median(r), r.min(),
          '; nodes p50/p99/max', np.median(nodes), np.quantile(nodes, .99), nodes.max())


if __name__ == '__main__':
    main()
