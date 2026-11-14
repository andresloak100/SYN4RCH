"""Can biological identity close the hippocampal gap (medoid ~366 um vs oracle ~186 um on VAL)?

DIAGNOSTIC on VAL mice only (train donors). Per neuron: its projection profile = fraction of axon cable in each
Allen region (rolled up to ontology depth <= 8), split ipsi/contra, from the local CCFv3 annotation_25.nrrd.

Rows labelled CEILING use the query's own axon (unavailable at prediction time); they measure how much the
information WOULD help. Rows labelled FORECAST use only the soma, so they are legitimate predictions.

Cell type: no per-neuron transcriptomic type exists for these fMOST neurons. Soma region / layer is already an
input. So "cell type" here can only mean projection class (defined from the axon) or soma subregion.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import pickle
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('hdm', ROOT / 'analysis/hipp_donor_model.py')
M_ = importlib.util.module_from_spec(spec)
sys.modules['hdm'] = M_
spec.loader.exec_module(M_)
WORK = M_.WORK
ANNOT = ROOT / 'data/raw/axon_pretrain_v1/calibration/annotation_25.nrrd'
GRAPH = ROOT / 'data/raw/atlas/structure_graph_1.json'
DEPTH = 8
MID_ML = 5700.0
SEED = 20261114


def read_nrrd(path):
    raw = path.read_bytes()
    head, body = raw.split(b'\n\n', 1)
    h = dict(l.split(': ', 1) for l in head.decode().splitlines() if ': ' in l)
    sizes = [int(s) for s in h['sizes'].split()]
    assert h['type'] == 'unsigned int' and h['encoding'] == 'gzip' and h['endian'] == 'little'
    arr = np.frombuffer(gzip.decompress(body), '<u4')
    return arr.reshape(sizes, order='F')          # NRRD: first axis fastest -> index [AP, DV, ML]


def ontology():
    g = json.loads(GRAPH.read_text())
    roots = g['msg'] if 'msg' in g else g
    up, depth, acr = {}, {}, {}
    def walk(n, path):
        p = path + [n['id']]
        depth[n['id']] = len(p) - 1
        acr[n['id']] = n['acronym']
        up[n['id']] = p[min(len(p) - 1, DEPTH)]
        for c in n.get('children', []):
            walk(c, p)
    for r in roots:
        walk(r, [])
    return up, acr


def profile(nid, ann, up, cols):
    xyz, par, typ = M_.tree(nid)
    w = M_.G.tree_stats((xyz, par, typ))['weights']
    v = np.clip(np.floor(xyz / 25.).astype(int), 0, np.array(ann.shape) - 1)
    ids = ann[v[:, 0], v[:, 1], v[:, 2]]
    side = (xyz[:, 2] > MID_ML).astype(int)        # canonical-left frame: ML > midline = contralateral
    out = Counter()
    for i, s, ww in zip(ids, side, w):
        if i == 0:
            continue
        out[(up.get(int(i), int(i)), s)] += ww
    for k in out:
        cols.setdefault(k, len(cols))
    return out


def hellinger_rows(profs, cols):
    X = np.zeros((len(profs), len(cols)))
    for r, p in enumerate(profs):
        for k, v in p.items():
            if k in cols:
                X[r, cols[k]] = v
    X /= np.maximum(X.sum(1, keepdims=True), 1e-12)
    return np.sqrt(X)


def main():
    from sklearn.cluster import KMeans
    from sklearn.ensemble import HistGradientBoostingClassifier
    split, meta = M_.load_meta()
    L = pickle.load(open(WORK / 'labels.pkl', 'rb'))
    train, val = sorted(L['ctr']), sorted(L['cva'])
    ann = read_nrrd(ANNOT)
    up, acr = ontology()
    # sanity: soma voxel region vs NeuroXiv soma-region label
    agree = []
    for q in val[:400]:
        s = np.clip(np.floor(meta[q]['soma'] / 25.).astype(int), 0, np.array(ann.shape) - 1)
        agree.append(acr.get(int(ann[tuple(s)]), '?') == meta[q]['region'])
    print('annotation orientation check: soma label agreement', float(np.mean(agree)), flush=True)
    assert np.mean(agree) > 0.6, 'annotation axis order looks wrong'

    cols = {}
    need = sorted(set(train) | set(val))
    prof = {n: profile(n, ann, up, cols) for n in need}
    # keep columns carrying >= 0.1% of total train cable
    tot = Counter()
    for n in train:
        tot.update(prof[n])
    s = sum(tot.values())
    keep = {k: i for i, k in enumerate(k for k, v in tot.items() if v / s >= 1e-3)}
    H = {n: row for n, row in zip(need, hellinger_rows([prof[n] for n in need], keep))}
    print('projection columns kept', len(keep), flush=True)

    A = [meta[q]['animal'] for q in val]
    def am(vals):
        return M_.animal_mean(vals, A)[0]
    def pick(q, j):
        return L['pva'][q][j]
    rows = {}
    def add(name, choose, kind):
        d = am([pick(q, choose(q))['distance_um'] for q in val])
        f = am([pick(q, choose(q))['f1_200'] for q in val])
        rows[name] = {'distance_um': d, 'f1_200': f, 'kind': kind}
    medi = {q: int(np.argmin(L['med'][q])) for q in val}
    add('medoid16', lambda q: medi[q], 'FORECAST (baseline)')
    add('ORACLE_best_of_16', lambda q: int(np.argmin([m['distance_um'] for m in L['pva'][q][:16]])), 'CEILING (geometry)')
    add('projection_profile_match', lambda q: int(np.argmin([np.linalg.norm(H[c] - H[q]) for c in L['cva'][q][:16]])),
        'CEILING (knows the query projection profile)')

    Xtr = np.stack([H[n] for n in train])
    feats = lambda n: np.r_[meta[n]['soma'], [float(meta[n]['coarse'] == c) for c in M_.COARSE_LIST]]
    for k in (8, 20, 40):
        km = KMeans(k, n_init=10, random_state=SEED).fit(Xtr)
        lab = {n: int(c) for n, c in zip(train, km.labels_)}
        vlab = {q: int(c) for q, c in zip(val, km.predict(np.stack([H[q] for q in val])))}
        clf = HistGradientBoostingClassifier(max_iter=300, random_state=SEED).fit(
            np.stack([feats(n) for n in train]), [lab[n] for n in train])
        vpred = dict(zip(val, clf.predict(np.stack([feats(q) for q in val]))))
        acc = float(np.mean([vpred[q] == vlab[q] for q in val]))
        def class_medoid(q, cls):
            cands = L['cva'][q][:16]
            same = [j for j, c in enumerate(cands) if lab[c] == cls]
            return min(same, key=lambda j: L['med'][q][j]) if same else medi[q]
        add(f'class{k}_true_then_medoid', lambda q: class_medoid(q, vlab[q]),
            f'CEILING (knows the query projection class, k={k})')
        add(f'class{k}_predicted_then_medoid', lambda q: class_medoid(q, vpred[q]),
            f'FORECAST (class predicted from soma, k={k}, val accuracy {acc:.2f})')
        rows[f'class{k}_predicted_then_medoid']['class_accuracy'] = acc
        rows[f'class{k}_true_then_medoid']['share_with_class_candidate'] = float(
            np.mean([any(lab[c] == vlab[q] for c in L['cva'][q][:16]) for q in val]))
    base, orc = rows['medoid16']['distance_um'], rows['ORACLE_best_of_16']['distance_um']
    for r in rows.values():
        r['gap_closed'] = (base - r['distance_um']) / (base - orc)
    (WORK / 'projection_info_val.json').write_text(json.dumps(rows, indent=1))
    for k, v in rows.items():
        print(f'{k:34s} {v["distance_um"]:7.1f}  F1 {v["f1_200"]:.3f}  gap closed {100*v["gap_closed"]:5.1f}%  [{v["kind"]}]')


if __name__ == '__main__':
    main()
