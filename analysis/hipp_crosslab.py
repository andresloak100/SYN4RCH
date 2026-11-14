"""Fresh cross-lab test of the hippocampal forecasters on MouseLight + SEU-Allen HPF neurons (NeuroXiv, CC BY-NC).

Queries come from other labs (MouseLight, SEU-Allen); donors stay the ION TRAIN pool (6,274 neurons), exactly as
developed. Nothing was tuned on these neurons.

  prep    select HPF axon neurons, drop any previously exposed MouseLight neuron and any cross-source duplicate
          (soma within 10 um of another query or of any ION neuron), extract axon (type 2 + ancestors), cache with
          the same simplification as hipp_cache.py. Structural information only; nothing is scored.
  freeze  write the protocol with hashes.
  test    single scoring.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import pickle
import sys
import zipfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]


def load(name, rel):
    s = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m
    s.loader.exec_module(m)
    return m


M_ = load('hdm', 'analysis/hipp_donor_model.py')
C_ = load('hc', 'analysis/hipp_consensus.py')
P_ = load('hpi', 'analysis/hipp_projection_info.py')
RAW = ROOT / 'data/raw/hippocampus_crosslab'
META = ROOT / 'data/raw/axon_pretrain_v1/metadata/neuroxiv_metadata.csv'
ML_IDS = ROOT / 'data/raw/axon_pretrain_v1/metadata/mouselight_latest_ids.json'
OUT = ROOT / 'data/interim/hippocampus/crosslab'
PROTOCOL = OUT / 'PROTOCOL_FROZEN.json'
ZIPS = {'MouseLight': 'MouseLight_full_SWC_CCFv3.zip', 'SEU-ALLEN': 'SEU-ALLEN_full_SWC_CCFv3.zip'}
DUP_UM = 10.0


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def exposed_mouselight():
    ids = set()
    for p in (ROOT / 'data/raw/mouselight_external_test').rglob('*'):
        if p.suffix == '.swc':
            ids.add(p.stem)
    for p in (ROOT / 'web').rglob('*.json'):
        try:
            txt = p.read_text()
        except Exception:
            continue
        import re
        ids |= set(re.findall(r'\bA[A-Z]\d{4}\b', txt))
    return ids


def build_one(args):
    nid, source, member, hemi = args
    text = zipfile.ZipFile(RAW / ZIPS[source]).read(member).decode()
    a = np.loadtxt([l for l in text.splitlines() if l.strip() and not l.startswith('#')], ndmin=2)
    ids = a[:, 0].astype(np.int64)
    typ = a[:, 1].astype(np.int64)
    xyz = a[:, 2:5].copy()
    pos = {int(i): k for k, i in enumerate(ids)}
    par = np.array([-1 if p == -1 else pos[int(p)] for p in a[:, 6]], np.int64)
    roots = np.flatnonzero(par == -1)
    if len(roots) != 1 or typ[roots[0]] != 1:
        return nid, {'ok': False, 'reason': f'roots={len(roots)}'}
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
    remap = np.full(len(par), -1)
    remap[idx] = np.arange(len(idx))
    sub_par = np.where(par[idx] >= 0, remap[np.maximum(par[idx], 0)], -1)
    sub_xyz = xyz[idx]
    if hemi == 'Right':
        sub_xyz[:, 2] = 11400. - sub_xyz[:, 2]
    # write a types-1/2 SWC-like array and reuse the exact hipp_cache simplification
    H = load('hcache', 'analysis/hipp_cache.py')
    n = len(sub_par)
    root = int(np.flatnonzero(sub_par == -1)[0])
    kids = [[] for _ in range(n)]
    for c in range(n):
        if sub_par[c] >= 0:
            kids[sub_par[c]].append(c)
    new_xyz, new_par = [sub_xyz[root]], [-1]
    stack = [(root, 0, 0.0)]
    while stack:
        node, last, arc = stack.pop()
        for c in kids[node]:
            d = arc + float(np.linalg.norm(sub_xyz[c] - sub_xyz[node]))
            if len(kids[c]) != 1 or d >= H.KEEP_ARC_UM:
                p0 = new_xyz[last]
                L = float(np.linalg.norm(sub_xyz[c] - p0))
                steps = max(1, int(np.ceil(L / H.MAX_EDGE_UM - 1e-9)))
                prev = last
                for s in range(1, steps):
                    new_xyz.append(p0 + (sub_xyz[c] - p0) * s / steps)
                    new_par.append(prev)
                    prev = len(new_xyz) - 1
                new_xyz.append(sub_xyz[c])
                new_par.append(prev)
                stack.append((c, len(new_xyz) - 1, 0.0))
            else:
                stack.append((c, last, d))
    X = np.asarray(new_xyz)
    Pp = np.asarray(new_par, np.int64)
    T = np.full(len(Pp), 2, np.int64)
    T[0] = 1
    np.savez(M_.CACHE / f'{nid}.npz', xyz_um=X, parent=Pp, types=T)
    cable = float(np.linalg.norm(X[1:] - X[Pp[1:]], axis=1).sum())
    return nid, {'ok': cable > 0, 'reason': '' if cable > 0 else 'zero_cable', 'soma': X[0].tolist(), 'nodes': len(Pp),
                 'cable_um': cable}


def phase_prep():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = {m['Morphology ID']: m for m in csv.DictReader(open(META))}
    hpf = sorted(k for k, m in meta.items() if m['Gray matter division'] == 'HPF'
                 and m['Data source'] in ZIPS and m['with axonal arbor'] == '1')
    members = {}
    for src, z in ZIPS.items():
        for i in zipfile.ZipFile(RAW / z).infolist():
            if i.filename.endswith('.swc'):
                members[i.filename.split('/')[-1][:-4]] = (src, i.filename)
    exposed = exposed_mouselight()
    sample = {n['idString']: n['sample']['id'] for n in json.loads(ML_IDS.read_text())['data']['searchNeurons']['neurons']}
    log = {}
    jobs = []
    for k in hpf:
        src, mem = members[k]
        short = k.split('_')[2] if src == 'MouseLight' else None
        if short and short in exposed:
            log[k] = 'previously_exposed_mouselight'
            continue
        jobs.append((k, src, mem, meta[k]['Hemisphere']))
    built = {}
    with ProcessPoolExecutor(6) as ex:
        for nid, r in ex.map(build_one, jobs):
            built[nid] = r
            if not r['ok']:
                log[nid] = r['reason']
    ok = [n for n in built if built[n]['ok']]
    # duplicates: against ION (all 10,023 somas, canonical frame) and within the query set
    _, ionmeta = M_.load_meta()
    ion_ids = sorted(ionmeta)
    kd = cKDTree(np.stack([ionmeta[i]['soma'] for i in ion_ids]))
    S = {n: np.asarray(built[n]['soma']) for n in ok}
    keep = []
    for n in sorted(ok, key=lambda n: (meta[n]['Data source'] != 'MouseLight', n)):
        if kd.query(S[n])[0] <= DUP_UM:
            log[n] = 'duplicate_of_ION'
            continue
        if any(np.linalg.norm(S[n] - S[m]) <= DUP_UM for m in keep):
            log[n] = 'cross_source_duplicate'
            continue
        keep.append(n)
    queries = []
    for n in sorted(keep):
        m = meta[n]
        src = m['Data source']
        group = (('ML:' + sample.get(n.split('_')[2], 'unknown')) if src == 'MouseLight'
                 else ('SEU:' + m['fMOST Brain ID']))
        queries.append({'id': n, 'source': src, 'group': group, 'soma': S[n].tolist(),
                        'region': m['Soma region (based on soma coordinates in ARA brain structures)']})
    groups = defaultdict(int)
    for q in queries:
        groups[q['group']] += 1
    res = {'schema': 'synarch_hipp_crosslab_queries_v1', 'created_utc': datetime.now(timezone.utc).isoformat(),
           'zips_sha256': {z: sha(RAW / z) for z in ZIPS.values()}, 'hpf_candidates': len(hpf),
           'excluded': log, 'queries': queries, 'groups': dict(groups),
           'note': 'Structural preparation only. No forecast has been compared with any of these neurons.'}
    with (OUT / 'queries.json').open('w') as f:
        json.dump(res, f, indent=1)
    from collections import Counter
    print('HPF candidates', len(hpf), '-> queries', len(queries), '; exclusions', Counter(log.values()))
    print('groups', len(groups), 'by source', Counter(q['source'] for q in queries), 'sizes', sorted(groups.values()))


def qmeta(q):
    return {'animal': q['group'], 'part': 'crosslab', 'soma': np.asarray(q['soma']), 'region': q['region'],
            'coarse': M_.coarse(q['region'])}


def phase_freeze():
    val_c = json.loads((M_.WORK / 'consensus_val.json').read_text())
    val_o = json.loads((M_.WORK / 'one_target_val.json').read_text())
    proto = {
        'schema': 'synarch_hipp_crosslab_protocol_v1', 'frozen_utc': datetime.now(timezone.utc).isoformat(),
        'queries_sha256': sha(OUT / 'queries.json'), 'script_sha256': sha(__file__),
        'consensus_script_sha256': sha(ROOT / 'analysis/hipp_consensus.py'),
        'donor_script_sha256': sha(ROOT / 'analysis/hipp_donor_model.py'),
        'projection_script_sha256': sha(ROOT / 'analysis/hipp_projection_info.py'),
        'labels_sha256': sha(M_.WORK / 'labels.pkl'), 'model_sha256': sha(M_.WORK / 'model.pkl'),
        'donor_pool': 'ION TRAIN neurons (6,274); candidates = 16 soma-nearest (32 for the gated model features)',
        'methods': {
            'consensus': 'M=16, TAU=0.25 (chosen on ION VAL: 358.7 um / F1 0.666)',
            'medoid16': 'most typical of the 16 nearest', 'nearest': 'soma-nearest donor',
            'gated_model': 'frozen ION model (delta from ION val)', 'star_2000': 'soma-only 26-ray star, 2000 um',
            'one_target_medoid': 'CONDITIONAL: query top target region given (from its own axon, i.e. as if known from '
                                 'a retrograde tracer); most typical of the 16 nearest that send >= 5% cable there '
                                 '(ION VAL: 329.3 um)',
            'ORACLE_best_of_16': 'ceiling'},
        'val_reference': {'consensus': val_c['table']['M16_tau0.25'], 'medoid16': val_c['table']['medoid16_reference'],
                          'one_target': val_o['rows']['0.05']},
        'primary': 'Group-equal mean distance_um (group = MouseLight sample or SEU brain). Consensus WINS if its '
                   'distance is below medoid16 with the paired-group bootstrap 95% interval of (consensus - medoid16) '
                   'entirely below 0. Secondary: F1@200; consensus vs nearest and star; one_target_medoid vs medoid16 '
                   '(conditional task, reported separately). Gated model reported for transfer only.',
        'single_attempt': 'Scored once. No retune, refit or re-run after seeing results.',
        'scope': 'Hippocampal neurons; ION donors, MouseLight/SEU-Allen queries; CC BY-NC; local lab use only.',
    }
    with PROTOCOL.open('x') as f:
        json.dump(proto, f, indent=1)
    print('protocol sha256', sha(PROTOCOL))


def test_job(args):
    q, cands16, cands32, one_t_ok = args
    pairs = [M_.metric(d, q) for d in cands16]
    med = M_.med_job((q, cands16))[1]
    order = [cands16[j] for j in np.argsort(med)]
    cons = C_.score(C_.consensus_tree(order, 16, 0.25), q)
    return q, {'pairs': pairs, 'med': med.tolist(), 'consensus': cons}


def phase_test():
    proto = json.loads(PROTOCOL.read_text())
    if (OUT / 'TEST_STARTED.json').exists():
        raise SystemExit('cross-lab test already attempted')
    for k, p in (('queries_sha256', OUT / 'queries.json'), ('script_sha256', Path(__file__)),
                 ('consensus_script_sha256', ROOT / 'analysis/hipp_consensus.py'),
                 ('donor_script_sha256', ROOT / 'analysis/hipp_donor_model.py'),
                 ('model_sha256', M_.WORK / 'model.pkl'), ('labels_sha256', M_.WORK / 'labels.pkl')):
        if sha(p) != proto[k]:
            raise SystemExit(f'{k} changed since freeze')
    (OUT / 'TEST_STARTED.json').write_text(json.dumps({'started_utc': datetime.now(timezone.utc).isoformat(),
                                                       'protocol_sha256': sha(PROTOCOL)}))
    Q = json.loads((OUT / 'queries.json').read_text())['queries']
    _, meta = M_.load_meta()
    for q in Q:
        meta[q['id']] = qmeta(q)
    L = pickle.load(open(M_.WORK / 'labels.pkl', 'rb'))
    model = pickle.load(open(M_.WORK / 'model.pkl', 'rb'))
    delta = json.loads((M_.WORK / 'PROTOCOL_FROZEN.json').read_text())['gate_delta']
    train = sorted(L['ctr'])
    ids = [q['id'] for q in Q]
    c32 = M_.candidates(meta, ids, train, M_.K)
    res = {}
    with ProcessPoolExecutor(6) as ex:
        for q, r in ex.map(test_job, [(q, c32[q][:16], c32[q], None) for q in ids], chunksize=2):
            res[q] = r
    # gated model, stars, one-target conditional
    stars = {q: M_.star_job((q, (2000.,)))[1]['2000'] for q in ids}
    ann = P_.read_nrrd(P_.ANNOT)
    up, acr = P_.ontology()
    cols = {}
    def frac(n):
        p = P_.profile(n, ann, up, cols)
        t = sum(p.values())
        return {k: v / t for k, v in p.items()}
    def home(n):
        sv = np.clip(np.floor(meta[n]['soma'] / 25.).astype(int), 0, np.array(ann.shape) - 1)
        i = int(ann[tuple(sv)])
        return (up.get(i, i), 0)
    per = defaultdict(list)
    detail = {}
    for q in ids:
        r = res[q]
        med = np.asarray(r['med'])
        mi = int(np.argmin(med))
        F = M_.features(meta, L['shapes'], L['central'], q, c32[q], med)
        gi = M_.gated_choice(model.predict(F), med, delta)
        fq = frac(q)
        fq = {k: v for k, v in fq.items() if k != home(q)}
        t = max(fq, key=fq.get) if fq else None
        ok = [j for j, c in enumerate(c32[q][:16]) if t is not None and frac(c).get(t, 0) >= 0.05]
        oi = min(ok, key=lambda j: med[j]) if ok else mi
        rows = {'consensus': r['consensus'], 'medoid16': r['pairs'][mi], 'nearest': r['pairs'][0],
                'gated_model': r['pairs'][gi], 'star_2000': stars[q], 'one_target_medoid': r['pairs'][oi],
                'ORACLE_best_of_16': min(r['pairs'], key=lambda m: m['distance_um'])}
        for k, v in rows.items():
            per[k].append(v)
        detail[q] = {'medoid_donor': c32[q][mi], 'consensus': r['consensus'], 'medoid16': r['pairs'][mi],
                     'one_target': acr.get(t[0], str(t[0])) + ('-contra' if t[1] else '') if t else None}
    G = [meta[q]['animal'] for q in ids]
    summary, by = {}, {}
    for k, rows in per.items():
        d, dby = M_.animal_mean([r['distance_um'] for r in rows], G)
        f, fby = M_.animal_mean([r['f1_200'] for r in rows], G)
        summary[k] = {'distance_um': d, 'f1_200': f}
        by[k] = (dby, fby)
    comps = {}
    for a, b in (('consensus', 'medoid16'), ('consensus', 'nearest'), ('consensus', 'star_2000'),
                 ('one_target_medoid', 'medoid16'), ('gated_model', 'medoid16')):
        dd = {g: by[a][0][g] - by[b][0][g] for g in by[a][0]}
        ff = {g: by[a][1][g] - by[b][1][g] for g in by[a][1]}
        dm, dci = M_.boot_ci(dd)
        fm, fci = M_.boot_ci(ff)
        comps[f'{a}_minus_{b}'] = {'distance': dm, 'distance_ci95': dci, 'f1': fm, 'f1_ci95': fci,
                                   'groups_better_distance': int(sum(v < 0 for v in dd.values())), 'groups': len(dd)}
    win = comps['consensus_minus_medoid16']['distance_ci95'][1] < 0
    out = {'schema': 'synarch_hipp_crosslab_results_v1', 'protocol_sha256': sha(PROTOCOL), 'queries': len(ids),
           'groups': len(set(G)), 'summary': summary, 'comparisons': comps, 'primary_win': win, 'per_query': detail,
           'finished_utc': datetime.now(timezone.utc).isoformat()}
    with (OUT / 'TEST_RESULTS.json').open('x') as f:
        json.dump(out, f, indent=1)
    print('TEST_RESULTS sha256', sha(OUT / 'TEST_RESULTS.json'), '| queries', len(ids), 'groups', len(set(G)))
    for k, v in summary.items():
        print(f'{k:20s} {v["distance_um"]:8.1f}  F1 {v["f1_200"]:.3f}')
    for k, v in comps.items():
        print(f'{k:34s} dist {v["distance"]:+7.1f} {np.round(v["distance_ci95"], 1)}  F1 {v["f1"]:+.3f} '
              f'{np.round(v["f1_ci95"], 3)}  better in {v["groups_better_distance"]}/{v["groups"]}')
    print('PRIMARY WIN:', win)


if __name__ == '__main__':
    {'prep': phase_prep, 'freeze': phase_freeze, 'test': phase_test}[sys.argv[1]]()
