"""Realistic-lab variant of the projection diagnostic (VAL only): suppose ONE target region of the query neuron is
known (e.g. it was retrogradely labelled from that region). Target = the region (ontology depth <= 8, ipsi/contra)
holding the most axon cable outside the soma's own region. Pick the most typical neighbour among candidates that
also send >= 5% of cable to that target; fall back to the medoid if none."""
import importlib.util, json, pickle, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
s = importlib.util.spec_from_file_location('hpi', ROOT / 'analysis/hipp_projection_info.py')
P = importlib.util.module_from_spec(s); sys.modules['hpi'] = P; s.loader.exec_module(P)
M_ = P.M_
split, meta = M_.load_meta()
L = pickle.load(open(M_.WORK / 'labels.pkl', 'rb'))
val = sorted(L['cva'])
ann = P.read_nrrd(P.ANNOT); up, acr = P.ontology()
need = sorted(set(val) | {c for q in val for c in L['cva'][q][:16]})
cols = {}
prof = {n: P.profile(n, ann, up, cols) for n in need}
def frac(n):
    t = sum(prof[n].values()); return {k: v / t for k, v in prof[n].items()}
F = {n: frac(n) for n in need}
def home(n):
    sv = np.clip(np.floor(meta[n]['soma'] / 25.).astype(int), 0, np.array(ann.shape) - 1)
    i = int(ann[tuple(sv)]); return (up.get(i, i), 0)
def top_target(n):
    f = {k: v for k, v in F[n].items() if k != home(n)}
    return max(f, key=f.get) if f else None
A = [meta[q]['animal'] for q in val]
medi = {q: int(np.argmin(L['med'][q])) for q in val}
def choose(q, thr):
    t = top_target(q)
    if t is None:
        return medi[q]
    ok = [j for j, c in enumerate(L['cva'][q][:16]) if F[c].get(t, 0) >= thr]
    return min(ok, key=lambda j: L['med'][q][j]) if ok else medi[q]
base = M_.animal_mean([L['pva'][q][medi[q]]['distance_um'] for q in val], A)[0]
orc = M_.animal_mean([min(m['distance_um'] for m in L['pva'][q][:16]) for q in val], A)[0]
out = {}
for thr in (0.02, 0.05, 0.1, 0.2):
    d = M_.animal_mean([L['pva'][q][choose(q, thr)]['distance_um'] for q in val], A)[0]
    f = M_.animal_mean([L['pva'][q][choose(q, thr)]['f1_200'] for q in val], A)[0]
    out[thr] = {'distance_um': d, 'f1_200': f, 'gap_closed': (base - d) / (base - orc)}
    print(f'one known target, donor needs >= {thr:.0%} cable there: {d:7.1f} um  F1 {f:.3f}  gap closed {100*(base-d)/(base-orc):5.1f}%')
print('neurons with no target outside home region:', sum(top_target(q) is None for q in val), 'of', len(val))
tt = [acr.get(top_target(q)[0], str(top_target(q)[0])) + ('-contra' if top_target(q)[1] else '') for q in val if top_target(q)]
from collections import Counter; print('most common top targets:', Counter(tt).most_common(8))
(M_.WORK / 'one_target_val.json').write_text(json.dumps({'medoid': base, 'oracle16': orc, 'rows': out}, indent=1))
