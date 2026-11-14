"""How much does each additional known target buy? (ION VAL mice only; development, not a test.)

A "known target" = a region (Allen ontology depth <= 8, ipsi/contra) that receives >= 5% of the query's axon cable
outside its home region; only its IDENTITY is known (as from a retrograde tracer), not the amount.

Two ways the targets become known:
  top_k     the k largest targets (optimistic: the lab happened to inject the main targets)
  random_k  k targets drawn at random from the query's real targets (realistic: the injected site is one of them)

Retrieval with exactly that information: among the candidate donors (16 or 32 soma-nearest from other mice), keep
those sending >= 5% cable to the most of the known targets; break ties by typicality (16: medoid score among the 16;
32: train-only donor centrality). Compared with the same-pool medoid (no target information) and the oracle.
"""
import importlib.util
import json
import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
s = importlib.util.spec_from_file_location('hpi', ROOT / 'analysis/hipp_projection_info.py')
P = importlib.util.module_from_spec(s)
sys.modules['hpi'] = P
s.loader.exec_module(P)
M_ = P.M_
THR = 0.05
KS = (1, 2, 3, 5)
SEED = 20261114


def main():
    split, meta = M_.load_meta()
    L = pickle.load(open(M_.WORK / 'labels.pkl', 'rb'))
    val = sorted(L['cva'])
    ann = P.read_nrrd(P.ANNOT)
    up, acr = P.ontology()
    need = sorted(set(val) | {c for q in val for c in L['cva'][q]})
    cols = {}
    F = {}
    for n in need:
        p = P.profile(n, ann, up, cols)
        t = sum(p.values())
        F[n] = {k: v / t for k, v in p.items()}

    def home(n):
        sv = np.clip(np.floor(meta[n]['soma'] / 25.).astype(int), 0, np.array(ann.shape) - 1)
        i = int(ann[tuple(sv)])
        return (up.get(i, i), 0)

    targets = {}
    for q in val:
        t = {k: v for k, v in F[q].items() if k != home(q) and v >= THR}
        targets[q] = sorted(t, key=t.get, reverse=True)
    rng = np.random.default_rng(SEED)
    rand_order = {q: [targets[q][j] for j in rng.permutation(len(targets[q]))] for q in val}
    A = [meta[q]['animal'] for q in val]

    def choose(q, known, pool):
        cands = L['cva'][q][:pool]
        if not known:
            return None
        hits = [sum(F[c].get(t, 0) >= THR for t in known) for c in cands]
        best = max(hits)
        if best == 0:
            return None
        ok = [j for j, h in enumerate(hits) if h == best]
        if pool == 16:
            return min(ok, key=lambda j: L['med'][q][j])
        return min(ok, key=lambda j: L['central'].get(cands[j], np.inf))

    def base(q, pool):
        if pool == 16:
            return int(np.argmin(L['med'][q]))
        cands = L['cva'][q][:pool]
        return min(range(pool), key=lambda j: L['central'].get(cands[j], np.inf))

    def am(js, key):
        return M_.animal_mean([L['pva'][q][j][key] for q, j in zip(val, js)], A)[0]

    out = {'n_val': len(val), 'mice': len(set(A)),
           'real_targets_per_neuron_p10_p50_p90': np.quantile([len(targets[q]) for q in val], [.1, .5, .9]).tolist()}
    for pool in (16, 32):
        b = [base(q, pool) for q in val]
        orc = [int(np.argmin([m['distance_um'] for m in L['pva'][q][:pool]])) for q in val]
        bd, od = am(b, 'distance_um'), am(orc, 'distance_um')
        rows = {'no_target_typical': {'distance_um': bd, 'f1_200': am(b, 'f1_200')},
                'ORACLE': {'distance_um': od, 'f1_200': am(orc, 'f1_200')}}
        for mode, order in (('top', targets), ('random', rand_order)):
            for k in KS:
                js = []
                for q, bj in zip(val, b):
                    j = choose(q, order[q][:k], pool)
                    js.append(bj if j is None else j)
                d = am(js, 'distance_um')
                rows[f'{mode}_{k}'] = {'distance_um': d, 'f1_200': am(js, 'f1_200'),
                                       'gap_closed': (bd - d) / (bd - od),
                                       'neurons_with_k_targets': int(sum(len(order[q]) >= k for q in val))}
        out[f'pool{pool}'] = rows
        print(f'\n--- candidate pool {pool} (no target {bd:.1f} um, oracle {od:.1f} um)')
        for k, v in rows.items():
            extra = f"  gap closed {100 * v['gap_closed']:5.1f}%  (neurons with >= k real targets: {v['neurons_with_k_targets']})" if 'gap_closed' in v else ''
            print(f"{k:18s} {v['distance_um']:7.1f} um  F1 {v['f1_200']:.3f}{extra}")
    print('\nreal targets per neuron p10/p50/p90:', out['real_targets_per_neuron_p10_p50_p90'])
    (M_.WORK / 'multi_target_val.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
