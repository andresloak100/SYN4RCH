"""Leave-one-mouse-out development over ALL 10,023 ION hippocampal neurons (240 mice). DEVELOPMENT ONLY.

The ION test split is already spent, so nothing here is a test. Every neuron is forecast from the 16 soma-nearest
neurons of OTHER mice (whole library as donors). Reports mouse-equal means.

  compute   pair metrics + medoid scores + projection profiles for all neurons (slow; cached to lomo.pkl)
  analyze   baselines, known-target curve, paper-class ceiling, per-region table
"""
from __future__ import annotations

import importlib.util
import json
import pickle
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
s = importlib.util.spec_from_file_location('hpi', ROOT / 'analysis/hipp_projection_info.py')
P = importlib.util.module_from_spec(s)
sys.modules['hpi'] = P
s.loader.exec_module(P)
M_ = P.M_
OUT = M_.H / 'lomo'
NEURON_INFO = ROOT / 'data/raw/hippocampus_neuroxiv/portal_meta/mouse.neuron.info.json'
THR = 0.05
SEED = 20261114
_ANN = None


def job(args):
    q, cands = args
    pairs = [M_.metric(d, q) for d in cands]
    med = M_.med_job((q, cands))[1]
    return q, pairs, med


def prof_job(n):
    global _ANN
    if _ANN is None:
        _ANN = (P.read_nrrd(P.ANNOT), P.ontology()[0])
    ann, up = _ANN
    cols = {}
    p = P.profile(n, ann, up, cols)
    t = sum(p.values())
    return n, {k: v / t for k, v in p.items()}


def phase_compute():
    OUT.mkdir(parents=True, exist_ok=True)
    _, meta = M_.load_meta()
    ids = sorted(meta)
    cands = M_.candidates(meta, ids, ids, 16)
    pairs, med = {}, {}
    with ProcessPoolExecutor(6) as ex:
        for i, (q, p, m) in enumerate(ex.map(job, [(q, cands[q]) for q in ids], chunksize=4), 1):
            pairs[q], med[q] = p, m
            if i % 1000 == 0:
                print(f'  pairs {i}/{len(ids)}', flush=True)
    with ProcessPoolExecutor(6) as ex:
        prof = dict(ex.map(prof_job, ids, chunksize=32))
    pickle.dump({'cands': cands, 'pairs': pairs, 'med': med, 'prof': prof}, open(OUT / 'lomo.pkl', 'wb'))
    print('saved', OUT / 'lomo.pkl')


def phase_analyze():
    _, meta = M_.load_meta()
    L = pickle.load(open(OUT / 'lomo.pkl', 'rb'))
    ids = sorted(L['pairs'])
    A = [meta[q]['animal'] for q in ids]
    ann = P.read_nrrd(P.ANNOT)
    up, acr = P.ontology()
    F = L['prof']
    info = json.loads(NEURON_INFO.read_text())['neuron_data']
    info = info.values() if isinstance(info, dict) else info
    pclass = {'ION_full_' + x['file'][:-4] + '_CCFv3': str(x.get('class')) for x in info}

    def home(n):
        sv = np.clip(np.floor(meta[n]['soma'] / 25.).astype(int), 0, np.array(ann.shape) - 1)
        i = int(ann[tuple(sv)])
        return (up.get(i, i), 0)

    targets = {}
    for q in ids:
        h = home(q)
        t = {k: v for k, v in F[q].items() if k != h and v >= THR}
        targets[q] = sorted(t, key=t.get, reverse=True)
    rng = np.random.default_rng(SEED)
    rand = {q: [targets[q][j] for j in rng.permutation(len(targets[q]))] for q in ids}
    medi = {q: int(np.argmin(L['med'][q])) for q in ids}
    orc = {q: int(np.argmin([m['distance_um'] for m in L['pairs'][q]])) for q in ids}

    def by_targets(q, known):
        if not known:
            return medi[q]
        hits = [sum(F[c].get(t, 0) >= THR for t in known) for c in L['cands'][q]]
        b = max(hits)
        if b == 0:
            return medi[q]
        return min((j for j, h in enumerate(hits) if h == b), key=lambda j: L['med'][q][j])

    def by_class(q, cls):
        ok = [j for j, c in enumerate(L['cands'][q]) if pclass.get(c) == cls]
        return min(ok, key=lambda j: L['med'][q][j]) if ok else medi[q]

    def majority_class(q):
        c = Counter(pclass.get(c) for c in L['cands'][q])
        return c.most_common(1)[0][0]

    choices = {'typical_neighbour': medi, 'nearest': {q: 0 for q in ids}, 'ORACLE_best_of_16': orc}
    for k in (1, 2, 3):
        choices[f'main_targets_{k}'] = {q: by_targets(q, targets[q][:k]) for q in ids}
        choices[f'random_targets_{k}'] = {q: by_targets(q, rand[q][:k]) for q in ids}
    choices['paper_class_known_CEILING'] = {q: by_class(q, pclass.get(q)) for q in ids}
    choices['paper_class_guessed_from_neighbours'] = {q: by_class(q, majority_class(q)) for q in ids}

    table, bym = {}, {}
    for name, ch in choices.items():
        d, dby = M_.animal_mean([L['pairs'][q][ch[q]]['distance_um'] for q in ids], A)
        f, _ = M_.animal_mean([L['pairs'][q][ch[q]]['f1_200'] for q in ids], A)
        table[name] = {'distance_um': d, 'f1_200': f}
        bym[name] = dby
    base, o = table['typical_neighbour']['distance_um'], table['ORACLE_best_of_16']['distance_um']
    for name, v in table.items():
        v['gap_closed'] = (base - v['distance_um']) / (base - o)
        if name not in ('typical_neighbour',):
            dm, ci = M_.boot_ci({a: bym[name][a] - bym['typical_neighbour'][a] for a in bym[name]})
            v['vs_typical'] = dm
            v['vs_typical_ci95'] = ci
    guess_acc = float(np.mean([majority_class(q) == pclass.get(q) for q in ids]))
    # per soma region
    reg = defaultdict(list)
    for q in ids:
        reg[meta[q]['coarse']].append(q)
    per_region = {}
    for r, qs in sorted(reg.items(), key=lambda kv: -len(kv[1])):
        per_region[r] = {'neurons': len(qs), 'mice': len({meta[q]['animal'] for q in qs}),
                         **{name: float(np.median([L['pairs'][q][choices[name][q]]['distance_um'] for q in qs]))
                            for name in ('typical_neighbour', 'main_targets_1', 'main_targets_2', 'ORACLE_best_of_16')}}
    out = {'status': 'DEVELOPMENT_ONLY_leave_one_mouse_out', 'neurons': len(ids), 'mice': len(set(A)), 'table': table,
           'paper_classes': len(set(pclass.values())), 'paper_class_guess_accuracy': guess_acc,
           'per_region_median_um': per_region}
    (OUT / 'lomo_report.json').write_text(json.dumps(out, indent=1))
    print(f'{len(ids)} neurons, {len(set(A))} mice (leave-one-mouse-out, development only)')
    for k, v in table.items():
        ci = f"  vs typical {v['vs_typical']:+6.1f} [{v['vs_typical_ci95'][0]:+.1f}, {v['vs_typical_ci95'][1]:+.1f}]" if 'vs_typical' in v else ''
        print(f"{k:38s} {v['distance_um']:6.1f} um  F1 {v['f1_200']:.3f}  gap {100*v['gap_closed']:5.1f}%{ci}")
    print('paper classes:', out['paper_classes'], '| class guessed from neighbours correct:', round(guess_acc, 3))
    print('per region (median um): region n mice typical main1 main2 oracle')
    for r, v in per_region.items():
        print(f"  {r:7s} {v['neurons']:5d} {v['mice']:4d} {v['typical_neighbour']:7.0f} {v['main_targets_1']:7.0f} {v['main_targets_2']:7.0f} {v['ORACLE_best_of_16']:7.0f}")


if __name__ == '__main__':
    {'compute': phase_compute, 'analyze': phase_analyze}[sys.argv[1]]()
