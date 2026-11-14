"""Donor-selection axon model for the ION hippocampus set (lab demo; CC BY-NC data, local only).

Task: given only a new neuron's soma position (and the atlas region at that position), forecast its full axon.
Method: take the K=32 soma-nearest TRAIN neurons from other mice, translate each so its soma sits on the query
soma, and pick the one a learned scorer predicts will be closest to the truth. The scorer is gradient boosting on
pair features (soma offset, regions, donor size/shape, donor typicality against the local population). It is
trained on TRAIN queries with leave-mouse-out candidates and tuned on VAL mice. TEST mice are scored once.

Metric: the frozen hypothalamic geometry_metrics (cable-weighted symmetric distance; F1 within 200 um), averaged
within mouse then across mice. Distance is primary: the hypothalamic diagnostic showed F1@200 alone is nearly
satisfiable by a soma-only star.

  python analysis/hipp_donor_model.py labels      # train+val pair labels and medoid/centrality features
  python analysis/hipp_donor_model.py fit         # fit on train, tune on val, save model, report val
  python analysis/hipp_donor_model.py freeze      # write protocol with hashes; required before test
  python analysis/hipp_donor_model.py test        # single scoring of test mice
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import pickle
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
H = ROOT / 'data/interim/hippocampus'
CACHE = H / 'cache'
WORK = H / 'donor_model'
SPLIT = H / 'split.json'
INV = H / 'inventory.csv'
PROTOCOL = WORK / 'PROTOCOL_FROZEN.json'
K = 32
K_MEDOID = 16
K_SELECT = 16          # the scorer chooses among the 16 soma-nearest, with their medoid scores as a feature
K_CENTRAL = 8
STAR_RAYS = 26
STAR_LENGTHS = (500., 1000., 2000., 3000., 4000., 6000.)
BOOT = 10000
SEED = 20261114
COARSE = {'CA1': 'CA1', 'CA2': 'CA2', 'CA3': 'CA3', 'SUB': 'SUB', 'ProS': 'ProS', 'PRE': 'PRESUB', 'POST': 'PRESUB',
          'PAR': 'PRESUB', 'HATA': 'OTHER'}
COARSE_LIST = ['CA1', 'CA2', 'CA3', 'DG', 'SUB', 'ProS', 'PRESUB', 'ENT', 'FIBER', 'OTHER']


def module(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


G = module('hipp_geometry', 'analysis/brm_geometry_audit.py')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def coarse(region):
    if region.startswith('DG'):
        return 'DG'
    if region.startswith('ENT'):
        return 'ENT'
    if region in ('alv', 'ec', 'fi', 'df', 'dhc', 'fp', 'hc', 'ccs', 'cing', 'or', 'scwm', 'amc', 'fiber tracts'):
        return 'FIBER'
    return COARSE.get(region, 'OTHER')


@lru_cache(maxsize=4000)
def tree(nid):
    with np.load(CACHE / f'{nid}.npz') as d:
        return d['xyz_um'], d['parent'], d['types']


@lru_cache(maxsize=4000)
def tstats(nid):
    t = tree(nid)
    s = G.tree_stats(t)
    return s['weights'], cKDTree(t[0] - t[0][0])     # KD tree in soma-relative frame


def metric(pred_id, truth_id):
    """geometry_metrics(pred translated onto truth soma, truth), using soma-relative frames (translation-exact)."""
    pw, pk = tstats(pred_id)
    tw, tk = tstats(truth_id)
    P = tree(pred_id)[0] - tree(pred_id)[0][0]
    T = tree(truth_id)[0] - tree(truth_id)[0][0]
    pd = tk.query(P)[0]
    td = pk.query(T)[0]
    fw, bw = float(np.average(pd, weights=pw)), float(np.average(td, weights=tw))
    pr, rc = float(np.average(pd <= 200., weights=pw)), float(np.average(td <= 200., weights=tw))
    return {'distance_um': (fw + bw) / 2, 'f1_200': 2 * pr * rc / (pr + rc) if pr + rc else 0.,
            'precision_200': pr, 'recall_200': rc}


def pair_job(args):
    q, donors = args
    return q, [metric(d, q) for d in donors]


def med_job(args):
    """Mean distance of each of the first K_SELECT candidates to the other candidates (medoid score)."""
    q, ds = args
    n = len(ds)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = metric(ds[i], ds[j])['distance_um']
    return q, D.sum(1) / (n - 1)


def run_med(cands, label):
    out = {}
    with ProcessPoolExecutor(6) as ex:
        for i, (q, v) in enumerate(ex.map(med_job, [(q, ds[:K_SELECT]) for q, ds in cands.items()], chunksize=4), 1):
            out[q] = v
            if i % 500 == 0:
                print(f'  {label} {i}/{len(cands)}', flush=True)
    return out


def fib_dirs(k):
    i = np.arange(k) + .5
    phi = np.arccos(1 - 2 * i / k)
    th = np.pi * (1 + 5 ** .5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


def star_job(args):
    q, lengths = args
    tw, tk = tstats(q)
    T = tree(q)[0] - tree(q)[0][0]
    out = {}
    for L in lengths:
        steps = int(np.ceil(L / 50.))
        pts = [np.zeros(3)] + [d * L * s / steps for d in fib_dirs(STAR_RAYS) for s in range(1, steps + 1)]
        P = np.asarray(pts)
        par = [-1]
        for r in range(STAR_RAYS):
            for s in range(steps):
                par.append(0 if s == 0 else len(par) - 1)
        pw = G.tree_stats((P, np.asarray(par, np.int64), np.full(len(par), 2, np.int64)))['weights']
        pd = tk.query(P)[0]
        td = cKDTree(P).query(T)[0]
        pr, rc = float(np.average(pd <= 200., weights=pw)), float(np.average(td <= 200., weights=tw))
        out[str(int(L))] = {'distance_um': (np.average(pd, weights=pw) + np.average(td, weights=tw)) / 2,
                            'f1_200': 2 * pr * rc / (pr + rc) if pr + rc else 0.}
    return q, out


# ------------------------------------------------------------------ data
def load_meta():
    split = json.loads(SPLIT.read_text())
    part = {a: k for k, v in split['animals'].items() for a in v}
    rows = {r['neuron_id']: r for r in csv.DictReader(open(INV))}
    meta = {}
    for nid, r in rows.items():
        soma = np.array([float(r['soma_x']), float(r['soma_y']), float(r['soma_z'])])
        if r['hemisphere'] == 'Right':
            soma[2] = 11400. - soma[2]
        meta[nid] = {'animal': r['animal'], 'part': part[r['animal']], 'soma': soma,
                     'region': r['soma_region'], 'coarse': coarse(r['soma_region'])}
    return split, meta


def candidates(meta, queries, pool, k):
    ids = np.array(pool)
    S = np.stack([meta[i]['soma'] for i in ids])
    A = np.array([meta[i]['animal'] for i in ids])
    kd = cKDTree(S)
    out = {}
    for q in queries:
        _, idx = kd.query(meta[q]['soma'], k=k + 300)
        sel = [ids[j] for j in idx if A[j] != meta[q]['animal']][:k]
        out[q] = sel
    return out


def run_pairs(cands, label):
    jobs = [(q, ds) for q, ds in cands.items()]
    res = {}
    with ProcessPoolExecutor(6) as ex:
        for i, (q, ms) in enumerate(ex.map(pair_job, jobs, chunksize=8), 1):
            res[q] = ms
            if i % 500 == 0:
                print(f'  {label} {i}/{len(jobs)}', flush=True)
    return res


def donor_shape(nid):
    t = tree(nid)
    X = t[0] - t[0][0]
    w = G.tree_stats(t)['weights']
    c = np.average(X, axis=0, weights=w)
    spread = np.sqrt(np.average((X - c) ** 2, axis=0, weights=w))
    return {'centroid': c, 'spread': spread, 'cable': float(w.sum()), 'extent': float(np.linalg.norm(X, axis=1).max()),
            'tips': int((np.bincount(t[1][1:], minlength=len(t[1])) == 0).sum())}


def shape_job(nid):
    return nid, donor_shape(nid)


# ------------------------------------------------------------------ features
def features(meta, shapes, central, q, donors, med):
    donors = donors[:K_SELECT]
    med = np.asarray(med)
    medrank = np.argsort(np.argsort(med))
    qs = meta[q]['soma']
    cents = np.stack([shapes[d]['centroid'] for d in donors])
    consensus = cents.mean(0)
    rows = []
    for rank, d in enumerate(donors):
        sh = shapes[d]
        delta = meta[d]['soma'] - qs
        f = [*delta, np.linalg.norm(delta), rank,
             float(meta[d]['region'] == meta[q]['region']), float(meta[d]['coarse'] == meta[q]['coarse']),
             np.log(sh['cable'] + 1), np.log(sh['extent'] + 1), np.log(sh['tips'] + 1), *sh['centroid'], *sh['spread'],
             np.linalg.norm(sh['centroid'] - consensus), central.get(d, np.nan),
             med[rank], med[rank] - med.min(), np.log(med[rank] / med.min()), medrank[rank]]
        f += [float(meta[q]['coarse'] == c) for c in COARSE_LIST]
        f += [float(meta[d]['coarse'] == c) for c in COARSE_LIST]
        rows.append(f)
    return np.asarray(rows, float)


def animal_mean(values, animals):
    by = defaultdict(list)
    for v, a in zip(values, animals):
        by[a].append(v)
    keys = sorted(by)
    return float(np.mean([np.mean(by[a]) for a in keys])), {a: float(np.mean(by[a])) for a in keys}


def boot_ci(diff_by_animal):
    keys = sorted(diff_by_animal)
    v = np.array([diff_by_animal[a] for a in keys])
    rng = np.random.default_rng(SEED)
    idx = rng.integers(0, len(v), (BOOT, len(v)))
    return float(v.mean()), np.quantile(v[idx].mean(1), [.025, .975]).tolist()


# ------------------------------------------------------------------ phases
def phase_labels():
    WORK.mkdir(parents=True, exist_ok=True)
    split, meta = load_meta()
    train = sorted(n for n in meta if meta[n]['part'] == 'train')
    val = sorted(n for n in meta if meta[n]['part'] == 'val')
    print('shapes for train donors', flush=True)
    with ProcessPoolExecutor(6) as ex:
        shapes = dict(ex.map(shape_job, train, chunksize=32))
    # donor centrality: mean distance to its K_CENTRAL nearest train donors from other mice (train-only geometry)
    cc = candidates(meta, train, train, K_CENTRAL)
    cres = run_pairs(cc, 'centrality')
    central = {d: float(np.mean([m['distance_um'] for m in cres[d]])) for d in train}
    ctr = candidates(meta, train, train, K)
    cva = candidates(meta, val, train, K)
    ptr = run_pairs(ctr, 'train pairs')
    pva = run_pairs(cva, 'val pairs')
    # medoid control among K_MEDOID nearest, for val
    med_jobs = {}
    for q in val:
        ds = cva[q][:K_MEDOID]
        for i, a in enumerate(ds):
            med_jobs[a + '|' + q] = (a, [b for b in ds if b != a])
    mres = {}
    keyed = list(med_jobs.items())
    with ProcessPoolExecutor(6) as ex:
        outs = list(ex.map(pair_job, [v for _, v in keyed], chunksize=8))
    for (key, _), (_, ms) in zip(keyed, outs):
        mres[key] = float(np.mean([m['distance_um'] for m in ms]))
    medoid_val = {q: min(cva[q][:K_MEDOID], key=lambda a: mres[a + '|' + q]) for q in val}
    with ProcessPoolExecutor(6) as ex:
        stars = dict(ex.map(star_job, [(q, STAR_LENGTHS) for q in val], chunksize=16))
    pickle.dump({'shapes': shapes, 'central': central, 'ctr': ctr, 'cva': cva, 'ptr': ptr, 'pva': pva,
                 'medoid_val': medoid_val, 'stars_val': stars}, open(WORK / 'labels.pkl', 'wb'))
    print('labels saved', sha(WORK / 'labels.pkl'))


def make_xy(meta, L, cands, pairs, queries):
    X, y, g = [], [], []
    for q in queries:
        F = features(meta, L['shapes'], L['central'], q, cands[q], L['med'][q])
        X.append(F)
        y.append(np.log([m['distance_um'] for m in pairs[q][:K_SELECT]]))
        g += [q] * len(F)
    return np.vstack(X), np.concatenate(y), g


DELTAS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, float('inf'))


def gated_choice(pred, med, delta):
    """Default to the medoid; switch to the scorer's pick only if it predicts >= delta lower log distance."""
    m = int(np.argmin(med))
    j = int(np.argmin(pred))
    return j if pred[m] - pred[j] > delta else m


def evaluate_choice(meta, queries, cands, pairs, choose):
    rows = []
    for q in queries:
        j = choose(q)
        rows.append(pairs[q][j])
    A = [meta[q]['animal'] for q in queries]
    d, dby = animal_mean([r['distance_um'] for r in rows], A)
    f, fby = animal_mean([r['f1_200'] for r in rows], A)
    return {'distance_um': d, 'f1_200': f}, dby, fby


def phase_fit():
    from sklearn.ensemble import HistGradientBoostingRegressor
    split, meta = load_meta()
    L = pickle.load(open(WORK / 'labels.pkl', 'rb'))
    train = sorted(L['ctr'])
    val = sorted(L['cva'])
    X, y, _ = make_xy(meta, L, L['ctr'], L['ptr'], train)
    # per-query centring: the scorer only needs to rank donors within a query
    yc = y.copy()
    o = 0
    for q in train:
        yc[o:o + K_SELECT] -= yc[o:o + K_SELECT].mean()
        o += K_SELECT
    best = None
    grid = [(it, lr, lf) for it in (300, 600) for lr in (.05, .1) for lf in (15, 31)]
    Fval = {q: features(meta, L['shapes'], L['central'], q, L['cva'][q], L['med'][q]) for q in val}
    for it, lr, lf in grid:
        m = HistGradientBoostingRegressor(max_iter=it, learning_rate=lr, max_leaf_nodes=lf, random_state=SEED)
        m.fit(X, yc)
        res, _, _ = evaluate_choice(meta, val, L['cva'], L['pva'], lambda q: int(np.argmin(m.predict(Fval[q]))))
        print(f'  grid it={it} lr={lr} leaves={lf}: val {res}', flush=True)
        if best is None or res['distance_um'] < best[0]['distance_um']:
            best = (res, (it, lr, lf), m)
    res, hp, model = best
    pickle.dump(model, open(WORK / 'model.pkl', 'wb'))
    Pval = {q: model.predict(Fval[q]) for q in val}
    gate = {}
    for dl in DELTAS:
        r, _, _ = evaluate_choice(meta, val, L['cva'], L['pva'], lambda q: gated_choice(Pval[q], L['med'][q], dl))
        sw = float(np.mean([gated_choice(Pval[q], L['med'][q], dl) != int(np.argmin(L['med'][q])) for q in val]))
        gate[str(dl)] = dict(r, switch_rate=sw)
        print(f'  gate delta={dl}: {r} switch_rate={sw:.3f}', flush=True)
    delta = min(DELTAS, key=lambda dl: gate[str(dl)]['distance_um'])
    A = [meta[q]['animal'] for q in val]
    controls = {
        'model': {k: gate[str(delta)][k] for k in ('distance_um', 'f1_200')},
        'scorer_ungated': res,
        'nearest': evaluate_choice(meta, val, L['cva'], L['pva'], lambda q: 0)[0],
        'medoid16': evaluate_choice(meta, val, L['cva'], L['pva'], lambda q: int(np.argmin(L['med'][q])))[0],
        'ORACLE_best_of_16': evaluate_choice(meta, val, L['cva'], L['pva'],
                                             lambda q: int(np.argmin([m['distance_um'] for m in L['pva'][q][:16]])))[0],
        'local_uniform32_expected': {k: animal_mean([np.mean([m[k] for m in L['pva'][q]]) for q in val], A)[0]
                                     for k in ('distance_um', 'f1_200')},
        'ORACLE_best_of_32': evaluate_choice(meta, val, L['cva'], L['pva'],
                                             lambda q: int(np.argmin([m['distance_um'] for m in L['pva'][q]])))[0],
    }
    for Ls in STAR_LENGTHS:
        controls[f'star_{int(Ls)}'] = {k: animal_mean([L['stars_val'][q][str(int(Ls))][k] for q in val], A)[0]
                                       for k in ('distance_um', 'f1_200')}
    best_star = min((k for k in controls if k.startswith('star_')), key=lambda k: controls[k]['distance_um'])
    out = {'hyperparameters': dict(zip(('max_iter', 'learning_rate', 'max_leaf_nodes'), hp)), 'val': controls,
           'gate_delta': delta, 'gate_sweep_val': gate,
           'best_star_on_val': best_star, 'model_sha256': sha(WORK / 'model.pkl'),
           'val_mice': len(set(A)), 'val_neurons': len(val)}
    (WORK / 'val_report.json').write_text(json.dumps(out, indent=1))
    for k, v in controls.items():
        print(f'{k:28s} dist {v["distance_um"]:8.1f}  F1 {v["f1_200"]:.4f}')
    print('best star on val:', best_star)


def phase_freeze():
    val = json.loads((WORK / 'val_report.json').read_text())
    proto = {
        'schema': 'synarch_hipp_donor_protocol_v1', 'frozen_utc': datetime.now(timezone.utc).isoformat(),
        'split_sha256': sha(SPLIT), 'inventory_sha256': sha(INV), 'script_sha256': sha(__file__),
        'cache_script_sha256': sha(ROOT / 'analysis/hipp_cache.py'),
        'geometry_sha256': sha(ROOT / 'analysis/brm_geometry_audit.py'),
        'labels_sha256': sha(WORK / 'labels.pkl'), 'model_sha256': sha(WORK / 'model.pkl'),
        'hyperparameters': val['hyperparameters'], 'star_length_from_val': val['best_star_on_val'],
        'gate_delta': val['gate_delta'],
        'rule': 'model = medoid16 unless the scorer predicts a candidate at least gate_delta lower in log distance.',
        'donor_pool': 'all TRAIN neurons; candidates are the K=32 soma-nearest from other mice',
        'methods': ['model', 'nearest', 'medoid16', 'local_uniform32_expected', 'best_val_star', 'ORACLE_best_of_16',
                    'ORACLE_best_of_32'],
        'development_history': 'v2 ungated scorer still lost to medoid16 on VAL (381.2 vs 366.1). v3 gates: medoid by '
                               'default, switch only on a confident predicted gain; gate tuned on VAL. Earlier: v1 scorer (32 candidates, no medoid feature) lost to medoid16 on VAL (387.3 vs 366.1 um). '
                               'v2 adds the medoid score as a feature and chooses among 16. This revision used VAL, so '
                               'VAL numbers are optimistic; TEST mice were never scored.',
        'primary': 'Animal-equal mean distance_um over the 48 TEST mice. The model WINS if its distance is lower than '
                   'each of nearest, medoid16, local_uniform32_expected and best_val_star, with the paired-mouse '
                   'bootstrap 95% interval of (model - control) entirely below 0. F1@200 is reported as secondary; '
                   'a model F1 below the best copying control is reported, not hidden.',
        'single_attempt': 'Test is scored once. No refit, re-tune or re-run after seeing test numbers.',
        'scope': 'Hippocampal neurons, ION/NeuroXiv, CC BY-NC. Local lab demo only. No public site, no claim beyond '
                 'this dataset and these mice.',
    }
    with PROTOCOL.open('x') as f:
        json.dump(proto, f, indent=1)
    print('protocol sha256', sha(PROTOCOL))


def phase_test():
    proto = json.loads(PROTOCOL.read_text())
    out = WORK / 'TEST_RESULTS.json'
    if out.exists() or (WORK / 'TEST_STARTED.json').exists():
        raise SystemExit('test already attempted')
    for key, path in (('split_sha256', SPLIT), ('script_sha256', Path(__file__)), ('model_sha256', WORK / 'model.pkl'),
                      ('labels_sha256', WORK / 'labels.pkl'), ('cache_script_sha256', ROOT / 'analysis/hipp_cache.py')):
        if sha(path) != proto[key]:
            raise SystemExit(f'{key} changed since freeze')
    (WORK / 'TEST_STARTED.json').write_text(json.dumps({'started_utc': datetime.now(timezone.utc).isoformat(),
                                                        'protocol_sha256': sha(PROTOCOL)}))
    split, meta = load_meta()
    L = pickle.load(open(WORK / 'labels.pkl', 'rb'))
    model = pickle.load(open(WORK / 'model.pkl', 'rb'))
    train = sorted(L['ctr'])
    test = sorted(n for n in meta if meta[n]['part'] == 'test')
    cte = candidates(meta, test, train, K)
    pte = run_pairs(cte, 'test pairs')
    med = run_med(cte, 'test medoid')
    Ls = proto['star_length_from_val'].split('_')[1]
    with ProcessPoolExecutor(6) as ex:
        stars = dict(ex.map(star_job, [(q, (float(Ls),)) for q in test], chunksize=16))
    A = [meta[q]['animal'] for q in test]
    dl = float(proto['gate_delta'])
    choice = {q: gated_choice(model.predict(features(meta, L['shapes'], L['central'], q, cte[q], med[q])), med[q], dl)
              for q in test}
    per = {}
    per['model'] = [pte[q][choice[q]] for q in test]
    per['nearest'] = [pte[q][0] for q in test]
    per['medoid16'] = [pte[q][int(np.argmin(med[q]))] for q in test]
    per['ORACLE_best_of_16'] = [min(pte[q][:16], key=lambda m: m['distance_um']) for q in test]
    per['local_uniform32_expected'] = [{k: float(np.mean([m[k] for m in pte[q]])) for k in ('distance_um', 'f1_200')}
                                       for q in test]
    per['best_val_star'] = [stars[q][Ls] for q in test]
    per['ORACLE_best_of_32'] = [min(pte[q], key=lambda m: m['distance_um']) for q in test]
    summary, by = {}, {}
    for m, rows in per.items():
        d, dby = animal_mean([r['distance_um'] for r in rows], A)
        f, fby = animal_mean([r['f1_200'] for r in rows], A)
        summary[m] = {'distance_um': d, 'f1_200': f}
        by[m] = (dby, fby)
    comps = {}
    for c in ('nearest', 'medoid16', 'local_uniform32_expected', 'best_val_star'):
        dd = {a: by['model'][0][a] - by[c][0][a] for a in by['model'][0]}
        ff = {a: by['model'][1][a] - by[c][1][a] for a in by['model'][1]}
        dm, dci = boot_ci(dd)
        fm, fci = boot_ci(ff)
        comps[c] = {'distance_model_minus_control': dm, 'distance_ci95': dci, 'f1_model_minus_control': fm,
                    'f1_ci95': fci, 'mice_model_better_distance': int(sum(v < 0 for v in dd.values())),
                    'distance_win': dci[1] < 0}
    win = all(v['distance_win'] for v in comps.values())
    res = {'schema': 'synarch_hipp_donor_test_v1', 'protocol_sha256': sha(PROTOCOL), 'test_mice': len(set(A)),
           'test_neurons': len(test), 'summary': summary, 'comparisons': comps, 'primary_win': win,
           'per_query': {q: {'chosen_donor': cte[q][choice[q]], 'model': per['model'][i],
                             'nearest': per['nearest'][i], 'oracle': per['ORACLE_best_of_32'][i]}
                         for i, q in enumerate(test)},
           'finished_utc': datetime.now(timezone.utc).isoformat()}
    with out.open('x') as f:
        json.dump(res, f, indent=1)
    print('TEST_RESULTS sha256', sha(out))
    for k, v in summary.items():
        print(f'{k:28s} dist {v["distance_um"]:8.1f}  F1 {v["f1_200"]:.4f}')
    print(json.dumps(comps, indent=1))
    print('PRIMARY WIN:', win)


def phase_medoid():
    L = pickle.load(open(WORK / 'labels.pkl', 'rb'))
    med = run_med(L['ctr'], 'train medoid')
    med.update(run_med(L['cva'], 'val medoid'))
    L['med'] = med
    pickle.dump(L, open(WORK / 'labels.pkl', 'wb'))
    print('labels+medoid saved', sha(WORK / 'labels.pkl'))


if __name__ == '__main__':
    {'medoid': phase_medoid, 'labels': phase_labels, 'fit': phase_fit, 'freeze': phase_freeze, 'test': phase_test}[sys.argv[1]]()
