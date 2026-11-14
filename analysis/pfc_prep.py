"""Structural preparation of ION prefrontal neurons for the frozen test (see data/interim/pfc/PROTOCOL_FROZEN.json).
Axon = type-2 nodes + ancestors; mirror Right; same simplification as hipp_cache.py. Nothing is scored here."""
import csv, json, re, sys, zipfile, hashlib, importlib.util
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
ZIP = ROOT / 'data/raw/pfc_neuroxiv/ION_PFC_full_SWC_CCFv3.zip'
META = ROOT / 'data/raw/axon_pretrain_v1/metadata/neuroxiv_metadata.csv'
CACHE = ROOT / 'data/interim/hippocampus/cache'          # shared tree cache; ids are globally unique
OUT = ROOT / 'data/interim/pfc'
KEEP_ARC_UM, MAX_EDGE_UM = 40.0, 50.0


def build(args):
    nid, member, hemi = args
    text = zipfile.ZipFile(ZIP).read(member).decode()
    rows = [l for l in text.splitlines() if l.strip() and not l.startswith('#')]
    try:
        a = np.loadtxt(rows, ndmin=2)
    except Exception as e:
        return nid, {'ok': False, 'reason': 'parse'}
    if a.shape[1] < 7 or not np.isfinite(a).all():
        return nid, {'ok': False, 'reason': 'fields_or_nonfinite'}
    ids, typ, xyz = a[:, 0].astype(np.int64), a[:, 1].astype(np.int64), a[:, 2:5].copy()
    pos = {int(i): k for k, i in enumerate(ids)}
    try:
        par = np.array([-1 if p == -1 else pos[int(p)] for p in a[:, 6]], np.int64)
    except KeyError:
        return nid, {'ok': False, 'reason': 'missing_parent'}
    roots = np.flatnonzero(par == -1)
    if len(roots) != 1 or typ[roots[0]] != 1:
        return nid, {'ok': False, 'reason': 'root'}
    if np.abs(xyz).max() > 1e6:
        return nid, {'ok': False, 'reason': 'coordinate_magnitude'}
    keep = typ == 2
    if not keep.any():
        return nid, {'ok': False, 'reason': 'no_axon'}
    for s in np.flatnonzero(keep):
        n = par[s]
        while n != -1 and not keep[n]:
            keep[n] = True
            n = par[n]
    keep[roots[0]] = True
    idx = np.flatnonzero(keep)
    remap = np.full(len(par), -1); remap[idx] = np.arange(len(idx))
    sp = np.where(par[idx] >= 0, remap[np.maximum(par[idx], 0)], -1)
    sx = xyz[idx]
    if hemi == 'Right':
        sx[:, 2] = 11400. - sx[:, 2]
    n = len(sp); root = int(np.flatnonzero(sp == -1)[0])
    kids = [[] for _ in range(n)]
    for c in range(n):
        if sp[c] >= 0:
            kids[sp[c]].append(c)
    nx, npar = [sx[root]], [-1]
    stack = [(root, 0, 0.0)]
    while stack:
        node, last, arc = stack.pop()
        for c in kids[node]:
            d = arc + float(np.linalg.norm(sx[c] - sx[node]))
            if len(kids[c]) != 1 or d >= KEEP_ARC_UM:
                p0 = nx[last]; L = float(np.linalg.norm(sx[c] - p0))
                steps = max(1, int(np.ceil(L / MAX_EDGE_UM - 1e-9))); prev = last
                for s in range(1, steps):
                    nx.append(p0 + (sx[c] - p0) * s / steps); npar.append(prev); prev = len(nx) - 1
                nx.append(sx[c]); npar.append(prev)
                stack.append((c, len(nx) - 1, 0.0))
            else:
                stack.append((c, last, d))
    X = np.asarray(nx); Pp = np.asarray(npar, np.int64)
    cable = float(np.linalg.norm(X[1:] - X[Pp[1:]], axis=1).sum())
    if cable <= 0:
        return nid, {'ok': False, 'reason': 'zero_cable'}
    T = np.full(len(Pp), 2, np.int64); T[0] = 1
    np.savez(CACHE / f'{nid}.npz', xyz_um=X, parent=Pp, types=T)
    return nid, {'ok': True, 'soma': X[0].tolist(), 'nodes': len(Pp), 'cable_um': cable}


def main():
    meta = {m['Morphology ID']: m for m in csv.DictReader(open(META))}
    Z = zipfile.ZipFile(ZIP)
    jobs = []
    for i in Z.infolist():
        if i.filename.endswith('.swc'):
            nid = i.filename.split('/')[-1][:-4]
            jobs.append((nid, i.filename, meta[nid]['Hemisphere']))
    res = {}
    with ProcessPoolExecutor(6) as ex:
        for k, (nid, r) in enumerate(ex.map(build, jobs, chunksize=8), 1):
            res[nid] = r
            if k % 1000 == 0:
                print(f'  {k}/{len(jobs)}', flush=True)
    queries, excluded = [], {}
    for nid in sorted(res):
        r = res[nid]
        if not r['ok']:
            excluded[nid] = r['reason']; continue
        reg = meta[nid]['Soma region (based on soma coordinates in ARA brain structures)']
        group = re.sub(r'(1|2/3|4|5|6a|6b)$', '', reg) or reg
        queries.append({'id': nid, 'animal': nid.split('_')[2], 'soma': r['soma'], 'region': reg, 'region_group': group})
    out = {'schema': 'synarch_region_queries_v1', 'dataset': 'pfc',
           'zip_sha256': hashlib.sha256(ZIP.read_bytes()).hexdigest(), 'files': len(jobs), 'excluded': excluded,
           'queries': queries, 'note': 'Structural preparation only; nothing scored.'}
    (OUT / 'queries.json').write_text(json.dumps(out, indent=1))
    nodes = [res[q['id']]['nodes'] for q in queries]
    print(len(queries), 'queries;', len({q['animal'] for q in queries}), 'mice; excluded', Counter(excluded.values()),
          '; nodes p50/p99', int(np.median(nodes)), int(np.quantile(nodes, .99)))
    print(Counter(q['region_group'] for q in queries).most_common(12))


if __name__ == '__main__':
    main()
