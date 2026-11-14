"""Consensus-blend axon forecaster (hippocampus, development on VAL mice only; the ION test mice are spent).

Idea: the M most typical neighbours (by medoid score) are each translated to the query soma. A point of one
neighbour's axon counts as "agreed" if at least a fraction TAU of the other M-1 neighbours have axon within 200 um
of it. The forecast keeps agreed parts only, built greedily in typicality order:
  - neighbour 1 (the medoid): keep agreed nodes still connected to the soma;
  - each further neighbour: add agreed nodes that are not already covered (no forecast point within COVER um),
    attaching to their own kept parent, else to the nearest forecast point if within GAP um, else dropped.
Edges are subdivided to <= 50 um. Every forecast point lies on some real neighbour's axon.

  python analysis/hipp_consensus.py dev      # grid over M, TAU on VAL; writes consensus_val.json
"""
from __future__ import annotations

import importlib.util
import json
import pickle
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('hdm', ROOT / 'analysis/hipp_donor_model.py')
M_ = importlib.util.module_from_spec(spec)
sys.modules['hdm'] = M_
spec.loader.exec_module(M_)
WORK = M_.WORK
R = 200.0
COVER = 50.0
GAP = 150.0
GRID_M = (4, 8, 16)
GRID_TAU = (0.25, 0.5, 0.75)


def rel(nid):
    xyz, par, _ = M_.tree(nid)
    return xyz - xyz[0], par


def consensus_tree(donors, m, tau):
    trees = [rel(d) for d in donors[:m]]
    kds = [cKDTree(t[0]) for t in trees]
    need = tau * (m - 1)
    out_xyz, out_par = [np.zeros(3)], [-1]
    built = None
    for ci, (X, P) in enumerate(trees):
        support = np.zeros(len(X))
        for cj, kd in enumerate(kds):
            if cj != ci:
                support += kd.query(X, distance_upper_bound=R)[0] <= R
        agreed = support >= need
        agreed[0] = True
        if built is not None:
            covered = built.query(X, distance_upper_bound=COVER)[0] <= COVER
        else:
            covered = np.zeros(len(X), bool)
        built_pts = np.asarray(out_xyz)
        bkd = cKDTree(built_pts)
        mapped = np.full(len(X), -1)
        mapped[0] = 0
        for i in range(1, len(X)):
            if not agreed[i]:
                continue
            p = P[i]
            if covered[i]:
                # already represented: route children through the nearest existing point
                d, j = bkd.query(X[i])
                mapped[i] = j if d <= COVER else -1
                continue
            if mapped[p] >= 0:
                anchor = mapped[p]
            else:
                d, j = bkd.query(X[i])
                if d > GAP:
                    continue
                anchor = j
            a = out_xyz[anchor]
            L = float(np.linalg.norm(X[i] - a))
            steps = max(1, int(np.ceil(L / 50. - 1e-9)))
            prev = anchor
            for s in range(1, steps):
                out_xyz.append(a + (X[i] - a) * s / steps)
                out_par.append(prev)
                prev = len(out_xyz) - 1
            out_xyz.append(X[i])
            out_par.append(prev)
            mapped[i] = len(out_xyz) - 1
        built = cKDTree(np.asarray(out_xyz))
    xyz = np.asarray(out_xyz)
    par = np.asarray(out_par, np.int64)
    return xyz, par, np.r_[1, np.full(len(par) - 1, 2)].astype(np.int64)


def score(tree, truth_id):
    tw, tk = M_.tstats(truth_id)
    T = M_.tree(truth_id)[0] - M_.tree(truth_id)[0][0]
    pw = M_.G.tree_stats(tree)['weights']
    P = tree[0]
    pd = tk.query(P)[0]
    td = cKDTree(P).query(T)[0]
    pr, rc = float(np.average(pd <= R, weights=pw)), float(np.average(td <= R, weights=tw))
    return {'distance_um': (float(np.average(pd, weights=pw)) + float(np.average(td, weights=tw))) / 2,
            'f1_200': 2 * pr * rc / (pr + rc) if pr + rc else 0., 'precision_200': pr, 'recall_200': rc,
            'cable_um': float(pw.sum())}


def dev_job(args):
    q, donors = args
    out = {}
    for m in GRID_M:
        for tau in GRID_TAU:
            out[f'M{m}_tau{tau}'] = score(consensus_tree(donors, m, tau), q)
    return q, out


def ordered(L, q):
    """Candidates in typicality order (lowest medoid score first)."""
    c = L['cva'][q][:16]
    return [c[j] for j in np.argsort(L['med'][q])]


def phase_dev():
    split, meta = M_.load_meta()
    L = pickle.load(open(WORK / 'labels.pkl', 'rb'))
    val = sorted(L['cva'])
    jobs = [(q, ordered(L, q)) for q in val]
    res = {}
    with ProcessPoolExecutor(6) as ex:
        for i, (q, o) in enumerate(ex.map(dev_job, jobs, chunksize=4), 1):
            res[q] = o
            if i % 250 == 0:
                print(f'  {i}/{len(val)}', flush=True)
    A = [meta[q]['animal'] for q in val]
    table = {}
    for key in res[val[0]]:
        table[key] = {k: M_.animal_mean([res[q][key][k] for q in val], A)[0]
                      for k in ('distance_um', 'f1_200', 'precision_200', 'recall_200')}
    med = {k: M_.animal_mean([L['pva'][q][int(np.argmin(L['med'][q]))][k] for q in val], A)[0]
           for k in ('distance_um', 'f1_200', 'precision_200', 'recall_200')}
    table['medoid16_reference'] = med
    (WORK / 'consensus_val.json').write_text(json.dumps({'val_mice': len(set(A)), 'table': table}, indent=1))
    pickle.dump(res, open(WORK / 'consensus_val_perquery.pkl', 'wb'))
    for k, v in sorted(table.items(), key=lambda kv: kv[1]['distance_um']):
        print(f'{k:22s} dist {v["distance_um"]:7.1f}  F1 {v["f1_200"]:.4f}  P {v["precision_200"]:.3f}  R {v["recall_200"]:.3f}')


if __name__ == '__main__':
    {'dev': phase_dev}[sys.argv[1]]()
