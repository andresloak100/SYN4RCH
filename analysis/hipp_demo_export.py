"""Export a non-cherry-picked sample of TEST neurons for the local 3D lab demo (CC BY-NC: local only).

Selection: sort the 2,079 test neurons by the model's distance error and take 24 at evenly spaced ranks, from best
to worst, so the demo shows typical and failing cases, not just good ones. For each: truth, model forecast,
nearest-copy forecast and the oracle-best-of-16, all in micrometres relative to the soma (soma-anchored).
Also exports every soma position (all splits) as faint context for the hippocampus shape.
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
import importlib.util as u
spec = u.spec_from_file_location('hdm', ROOT / 'analysis/hipp_donor_model.py')
M = u.module_from_spec(spec)
sys.modules['hdm'] = M
spec.loader.exec_module(M)

DEMO = ROOT / 'data/interim/hippocampus/demo'
N_SHOW = 24


def pack(nid):
    xyz, par, _ = M.tree(nid)
    rel = np.round(xyz - xyz[0]).astype(int)
    return {'xyz': rel.ravel().tolist(), 'parent': par.tolist()}


def main():
    DEMO.mkdir(parents=True, exist_ok=True)
    res = json.loads((M.WORK / 'TEST_RESULTS.json').read_text())
    split, meta = M.load_meta()
    L = __import__('pickle').load(open(M.WORK / 'labels.pkl', 'rb'))
    train = sorted(L['ctr'])
    per = res['per_query']
    order = sorted(per, key=lambda q: per[q]['model']['distance_um'])
    idx = np.linspace(0, len(order) - 1, N_SHOW).round().astype(int)
    picks = [order[i] for i in idx]
    cands = M.candidates(meta, picks, train, M.K)
    neurons = []
    for rank, q in zip(idx, picks):
        c = cands[q][:16]
        oracle = min(range(16), key=lambda j: M.metric(c[j], q)['distance_um'])
        r = per[q]
        neurons.append({
            'id': q, 'mouse': meta[q]['animal'], 'region': meta[q]['region'],
            'soma': np.round(meta[q]['soma']).astype(int).tolist(),
            'percentile': round(100 * rank / (len(order) - 1)),
            'truth': pack(q),
            'model': dict(pack(r['chosen_donor']), donor=r['chosen_donor'], score=r['model']),
            'nearest': dict(pack(c[0]), donor=c[0], score=r['nearest']),
            'oracle': dict(pack(c[oracle]), donor=c[oracle], score=M.metric(c[oracle], q)),
        })
    somas = np.round([meta[n]['soma'] for n in sorted(meta)]).astype(int).ravel().tolist()
    summary = res['summary']
    data = {'neurons': neurons, 'somas': somas, 'summary': summary,
            'note': 'TEST mice only (never used for fitting). 24 neurons at evenly spaced error ranks, best to worst. '
                    'ION hippocampus via NeuroXiv, CC BY-NC 4.0: local lab use only, do not redistribute. '
                    'Cite NeuroXiv (SEU) and Qiu et al. 2024 Science.'}
    (DEMO / 'demo.json').write_text(json.dumps(data, separators=(',', ':')))
    for f in ('three.r128.min.js', 'THREE-LICENSE.txt'):
        shutil.copy(ROOT / 'web/synarch-3d/vendor' / f, DEMO / f)
    print('wrote', DEMO / 'demo.json', (DEMO / 'demo.json').stat().st_size // 1024, 'KB;', len(neurons), 'neurons')




def add_targets():
    """Known-target forecaster for the demo neurons: for every target region that one of the 16 nearest donors
    sends >= 5% of its cable to, the forecast is the most typical donor among those that do. The neuron's real top
    target (from its own axon) is pre-selected. ION test neurons are already spent, so these scores are display only."""
    spec2 = u.spec_from_file_location('hpi', ROOT / 'analysis/hipp_projection_info.py')
    P = u.module_from_spec(spec2)
    sys.modules['hpi'] = P
    spec2.loader.exec_module(P)
    import pickle
    data = json.loads((DEMO / 'demo.json').read_text())
    split, meta = M.load_meta()
    L = pickle.load(open(M.WORK / 'labels.pkl', 'rb'))
    train = sorted(L['ctr'])
    ann = P.read_nrrd(P.ANNOT)
    up, acr = P.ontology()
    g = json.loads(P.GRAPH.read_text())
    names = {}
    def walk(n):
        names[n['id']] = n['name']
        for c in n.get('children', []):
            walk(c)
    for r in (g['msg'] if 'msg' in g else g):
        walk(r)
    cols = {}
    def frac(n):
        p = P.profile(n, ann, up, cols)
        t = sum(p.values())
        return {k: v / t for k, v in p.items()}
    def home(n):
        sv = np.clip(np.floor(meta[n]['soma'] / 25.).astype(int), 0, np.array(ann.shape) - 1)
        i = int(ann[tuple(sv)])
        return (up.get(i, i), 0)
    def label(k):
        return (acr.get(k[0], str(k[0])) + (' (other side)' if k[1] else ''), names.get(k[0], '') + (' – contralateral' if k[1] else ''))
    for n in data['neurons']:
        q = n['id']
        c16 = M.candidates(meta, [q], train, 16)[q]
        med = M.med_job((q, c16))[1]
        fr = {c: frac(c) for c in c16}
        fq = {k: v for k, v in frac(q).items() if k != home(q)}
        true_t = max(fq, key=fq.get) if fq else None
        opts = {}
        for c in c16:
            for k, v in fr[c].items():
                if k != home(q) and v >= 0.05:
                    opts.setdefault(k, []).append(c16.index(c))
        donors, targets = {}, []
        for k, js in opts.items():
            j = min(js, key=lambda j: med[j])
            d = c16[j]
            if d not in donors:
                donors[d] = pack(d)
            a, full = label(k)
            targets.append({'acr': a, 'name': full, 'donor': d, 'n_support': len(js),
                            'truth_fraction': round(fq.get(k, 0.0), 3), 'is_true_top': k == true_t,
                            'score': M.metric(d, q)})
        targets.sort(key=lambda t: (-t['is_true_top'], -t['truth_fraction'], t['acr']))
        n['targets'] = targets
        n['target_donors'] = donors
        n['true_top_target'] = label(true_t)[0] if true_t else None
    xl = json.loads((ROOT / 'data/interim/hippocampus/crosslab/TEST_RESULTS.json').read_text())
    data['crosslab'] = {'summary': xl['summary'], 'queries': xl['queries'], 'groups': xl['groups'],
                        'one_target_vs_medoid': xl['comparisons']['one_target_medoid_minus_medoid16']}
    (DEMO / 'demo.json').write_text(json.dumps(data, separators=(',', ':')))
    print('targets added;', (DEMO / 'demo.json').stat().st_size // 1024, 'KB; options per neuron',
          [len(n['targets']) for n in data['neurons']])


if __name__ == '__main__':
    if sys.argv[1:] == ['targets']:
        add_targets()
    else:
        main()
