"""Frozen leave-one-mouse-out scoring of the retrieval forecasters on a NEW region (first use: ION prefrontal cortex).

The method is exactly the one developed on hippocampus (analysis/hipp_lomo.py), with no fitted parameters:
  candidates        16 soma-nearest neurons from OTHER mice of the same dataset
  typical           candidate with the lowest mean distance to the other 15 (medoid)
  known targets     target = region (Allen ontology rolled up to depth <= 8, ipsi/contra) holding >= 5% of the
                    query's axon cable, excluding the soma's own region; only identity is known. Pick the candidates
                    sending >= 5% of their cable to the most known targets; break ties by typicality; fall back to
                    typical when none match.
  main_k / random_k the k largest targets / k targets drawn at random (fixed seed) from the query's real targets
Secondary, pre-specified: the same at ontology depth <= 6 (area level; in isocortex depth 8 is layer level).

  selftest   run the identical code on 12 hippocampal validation mice (already-exposed data) to prove it executes
  score      single scoring of the prepared dataset; refuses to run twice or if anything frozen changed
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
s = importlib.util.spec_from_file_location('hpi', ROOT / 'analysis/hipp_projection_info.py')
P = importlib.util.module_from_spec(s)
sys.modules['hpi'] = P
s.loader.exec_module(P)
M_ = P.M_
THR = 0.05
K = 16
KS = (1, 2, 3)
SEED = 20261114
STAR_UM = 2000.0
_ANN = None


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def job(args):
    q, cands = args
    return q, [M_.metric(d, q) for d in cands], M_.med_job((q, cands))[1], M_.star_job((q, (STAR_UM,)))[1][str(int(STAR_UM))]


def prof_job(n):
    """Raw cable per (leaf annotation id, side); rolled up later for each ontology depth."""
    global _ANN
    if _ANN is None:
        _ANN = P.read_nrrd(P.ANNOT)
    xyz, par, typ = M_.tree(n)
    w = M_.G.tree_stats((xyz, par, typ))['weights']
    v = np.clip(np.floor(xyz / 25.).astype(int), 0, np.array(_ANN.shape) - 1)
    ids = _ANN[v[:, 0], v[:, 1], v[:, 2]]
    side = (xyz[:, 2] > P.MID_ML).astype(int)
    out = defaultdict(float)
    for i, sd, ww in zip(ids, side, w):
        if i:
            out[(int(i), int(sd))] += float(ww)
    return n, dict(out)


def rollup_map(depth):
    g = json.loads(P.GRAPH.read_text())
    up = {}
    def walk(n, path):
        p = path + [n['id']]
        up[n['id']] = p[min(len(p) - 1, depth)]
        for c in n.get('children', []):
            walk(c, p)
    for r in (g['msg'] if 'msg' in g else g):
        walk(r, [])
    return up


def run(meta, label):
    ids = sorted(meta)
    A = [meta[q]['animal'] for q in ids]
    cands = M_.candidates(meta, ids, ids, K)
    pairs, med, star = {}, {}, {}
    with ProcessPoolExecutor(6) as ex:
        for i, (q, p, m, st) in enumerate(ex.map(job, [(q, cands[q]) for q in ids], chunksize=4), 1):
            pairs[q], med[q], star[q] = p, m, st
            if i % 1000 == 0:
                print(f'  {label} {i}/{len(ids)}', flush=True)
    with ProcessPoolExecutor(6) as ex:
        raw = dict(ex.map(prof_job, ids, chunksize=32))
    ann = P.read_nrrd(P.ANNOT)
    medi = {q: int(np.argmin(med[q])) for q in ids}
    choices = {'typical_neighbour': medi, 'nearest': {q: 0 for q in ids},
               'ORACLE_best_of_16': {q: int(np.argmin([m['distance_um'] for m in pairs[q]])) for q in ids}}
    ntargets = {}
    for depth in (8, 6):
        up = rollup_map(depth)
        F = {}
        for n in ids:
            acc = defaultdict(float)
            for (i, sd), w in raw[n].items():
                acc[(up.get(i, i), sd)] += w
            t = sum(acc.values())
            F[n] = {k: v / t for k, v in acc.items()} if t > 0 else {}
        def home(n):
            sv = np.clip(np.floor(meta[n]['soma'] / 25.).astype(int), 0, np.array(ann.shape) - 1)
            i = int(ann[tuple(sv)])
            return (up.get(i, i), 0)
        targets = {}
        for q in ids:
            h = home(q)
            t = {k: v for k, v in F[q].items() if k != h and v >= THR}
            targets[q] = sorted(t, key=lambda k: (-t[k], k))
        rng = np.random.default_rng(SEED)
        rand = {q: [targets[q][j] for j in rng.permutation(len(targets[q]))] for q in ids}
        ntargets[depth] = np.quantile([len(targets[q]) for q in ids], [.1, .5, .9]).tolist()
        def pick(q, known):
            if not known:
                return medi[q]
            hits = [sum(F[c].get(t, 0) >= THR for t in known) for c in cands[q]]
            b = max(hits)
            if b == 0:
                return medi[q]
            return min((j for j, h in enumerate(hits) if h == b), key=lambda j: med[q][j])
        tag = '' if depth == 8 else '_area_level'
        for k in KS:
            choices[f'main_targets_{k}{tag}'] = {q: pick(q, targets[q][:k]) for q in ids}
            choices[f'random_targets_{k}{tag}'] = {q: pick(q, rand[q][:k]) for q in ids}
    table, byd, byf = {}, {}, {}
    def add(name, rows):
        d, dby = M_.animal_mean([r['distance_um'] for r in rows], A)
        f, fby = M_.animal_mean([r['f1_200'] for r in rows], A)
        table[name] = {'distance_um': d, 'f1_200': f}
        byd[name], byf[name] = dby, fby
    for name, ch in choices.items():
        add(name, [pairs[q][ch[q]] for q in ids])
    add('star_2000', [star[q] for q in ids])
    for name, v in table.items():
        if name != 'typical_neighbour':
            dm, ci = M_.boot_ci({a: byd[name][a] - byd['typical_neighbour'][a] for a in byd[name]})
            fm, fci = M_.boot_ci({a: byf[name][a] - byf['typical_neighbour'][a] for a in byf[name]})
            v.update(vs_typical_distance=dm, vs_typical_distance_ci95=ci, vs_typical_f1=fm, vs_typical_f1_ci95=fci)
    reg = defaultdict(list)
    for q in ids:
        reg[meta[q]['region_group']].append(q)
    per_region = {r: {'neurons': len(qs), 'mice': len({meta[q]['animal'] for q in qs}),
                      **{nm: float(np.median([pairs[q][choices[nm][q]]['distance_um'] for q in qs]))
                         for nm in ('typical_neighbour', 'main_targets_1', 'main_targets_2', 'ORACLE_best_of_16')}}
                  for r, qs in sorted(reg.items(), key=lambda kv: -len(kv[1]))}
    primary = {
        'one_main_target_beats_typical': table['main_targets_1']['vs_typical_distance_ci95'][1] < 0,
        'typical_beats_nearest': table['nearest']['vs_typical_distance_ci95'][0] > 0,
        'typical_beats_star': table['star_2000']['vs_typical_distance_ci95'][0] > 0}
    primary['all'] = all(primary.values())
    return {'neurons': len(ids), 'mice': len(set(A)), 'table': table, 'primary': primary,
            'real_targets_per_neuron_p10_p50_p90': ntargets, 'per_region_median_um': per_region,
            'per_query': {q: {'typical': pairs[q][medi[q]]['distance_um'],
                              'main1': pairs[q][choices['main_targets_1'][q]]['distance_um'],
                              'oracle': pairs[q][choices['ORACLE_best_of_16'][q]]['distance_um']} for q in ids}}


def show(res):
    print(f"{res['neurons']} neurons, {res['mice']} mice")
    for k, v in res['table'].items():
        ci = (f"  vs typical {v['vs_typical_distance']:+7.1f} [{v['vs_typical_distance_ci95'][0]:+.1f}, "
              f"{v['vs_typical_distance_ci95'][1]:+.1f}]") if 'vs_typical_distance' in v else ''
        print(f"{k:30s} {v['distance_um']:7.1f} um  F1 {v['f1_200']:.3f}{ci}")
    print('targets per neuron p10/p50/p90 (depth 8, 6):', res['real_targets_per_neuron_p10_p50_p90'])
    print('per region (median um): region n mice typical main1 main2 oracle')
    for r, v in list(res['per_region_median_um'].items())[:14]:
        print(f"  {r:8s} {v['neurons']:5d} {v['mice']:4d} {v['typical_neighbour']:7.0f} {v['main_targets_1']:7.0f} "
              f"{v['main_targets_2']:7.0f} {v['ORACLE_best_of_16']:7.0f}")
    print('PRIMARY:', res['primary'])


def phase_selftest():
    split, meta = M_.load_meta()
    mice = split['animals']['val'][:12]
    sub = {n: dict(m, region_group=m['coarse']) for n, m in meta.items() if m['animal'] in mice}
    show(run(sub, 'selftest'))


def phase_score(dataset):
    D = ROOT / 'data/interim' / dataset
    proto = json.loads((D / 'PROTOCOL_FROZEN.json').read_text())
    if (D / 'TEST_STARTED.json').exists():
        raise SystemExit('already attempted')
    for k, p in (('scoring_script_sha256', Path(__file__)), ('projection_script_sha256', ROOT / 'analysis/hipp_projection_info.py'),
                 ('donor_script_sha256', ROOT / 'analysis/hipp_donor_model.py'),
                 ('geometry_sha256', ROOT / 'analysis/brm_geometry_audit.py')):
        if sha(p) != proto[k]:
            raise SystemExit(f'{k} changed since freeze')
    Q = json.loads((D / 'queries.json').read_text())
    (D / 'TEST_STARTED.json').write_text(json.dumps({'started_utc': datetime.now(timezone.utc).isoformat(),
                                                     'protocol_sha256': sha(D / 'PROTOCOL_FROZEN.json'),
                                                     'queries_sha256': sha(D / 'queries.json')}))
    meta = {q['id']: {'animal': q['animal'], 'soma': np.asarray(q['soma']), 'region_group': q['region_group']}
            for q in Q['queries']}
    res = run(meta, dataset)
    res.update(schema='synarch_region_lomo_test_v1', dataset=dataset, protocol_sha256=sha(D / 'PROTOCOL_FROZEN.json'),
               queries_sha256=sha(D / 'queries.json'), finished_utc=datetime.now(timezone.utc).isoformat())
    with (D / 'TEST_RESULTS.json').open('x') as f:
        json.dump(res, f, indent=1)
    print('TEST_RESULTS sha256', sha(D / 'TEST_RESULTS.json'))
    show(res)


if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        phase_selftest()
    else:
        phase_score(sys.argv[2])
