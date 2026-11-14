"""Read-only vacuity diagnostic on the closed 262-query exploratory axon comparison.

POST-HOC, ON ALREADY-EXPOSED DATA. Nothing here is a test, a method, or a promotion candidate. It asks one
question: how much of the frozen metric (cable-weighted symmetric distance, F1 at 200 um) can be earned by
forecasts that know nothing except the soma position?

It never writes inside the sealed batch or the parent output, never calls the evaluator's scoring function,
never refits anything. Truth is rebuilt from the stored raw SWC with the frozen adapter, and every rebuilt truth
must match the per-query digest recorded in SCORE_RECEIPT.json, or the run stops.

  /opt/anaconda3/bin/python -I -u analysis/axon_vacuity_diagnostic.py
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/interim/axon_pretrain_hybrid/source_qc_exploratory'
BATCH = OUT / 'batch'
RAW = ROOT / 'data/raw/axon_pretrain_hybrid/external_test'
COHORT = ROOT / 'data/raw/axon_pretrain_v1/metadata/spcd_hy_candidate_cohort.json'
REPORT = ROOT / 'reports/axon_vacuity_diagnostic.json'
EXPECTED_RESULTS_SHA = 'af722dfd'
RADIUS = 200.0
STAR_RAYS = 26
STAR_LENGTHS = (250., 500., 1000., 2000., 4000.)
TRUNCATE_RADII = (250., 500., 1000., 2000., 4000.)
BANDS = (0., 250., 500., 1000., 2000., 4000., np.inf)
ORACLE_K = 16


def module(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


G = module('vac_geometry', 'analysis/brm_geometry_audit.py')
S = module('vac_adapter', 'analysis/axon_pretrain_truth_adapter.py')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def digest(tree):
    h = hashlib.sha256()
    for a in tree:
        v = np.asarray(a)
        h.update(f'{v.dtype}{v.shape}'.encode())
        h.update(np.ascontiguousarray(v).tobytes())
    return h.hexdigest()


def load_truth(row):
    full = S.parse_swc_text((RAW / row['file']).read_text())
    axon = S.extract_axon(full, row['source_reconstruction_type'])
    t = S.for_evaluator(axon, row['source_hemisphere'])
    return (t['xyz_um'], t['parent'], t['types'])


def fib_dirs(k):
    i = np.arange(k) + .5
    phi = np.arccos(1 - 2 * i / k)
    th = np.pi * (1 + 5 ** .5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


def star(soma, length):
    """Soma plus STAR_RAYS straight rays of fixed length, sampled at <= 50 um. Knows only the soma."""
    steps = int(np.ceil(length / 50.))
    xyz, par = [np.asarray(soma, float)], [-1]
    for d in fib_dirs(STAR_RAYS):
        prev = 0
        for s in range(1, steps + 1):
            xyz.append(np.asarray(soma) + d * length * s / steps)
            par.append(prev)
            prev = len(xyz) - 1
    return np.asarray(xyz), np.asarray(par, np.int64), np.full(len(par), 2, np.int64)


def truncate(tree, r):
    xyz, par, typ = tree
    within = np.linalg.norm(xyz - xyz[0], axis=1) <= r
    keep = np.zeros(len(par), bool)
    keep[0] = True
    for i in range(1, len(par)):
        keep[i] = within[i] and keep[par[i]]
    idx = np.flatnonzero(keep)
    remap = np.full(len(par), -1)
    remap[idx] = np.arange(len(idx))
    newpar = np.where(par[idx] >= 0, remap[np.maximum(par[idx], 0)], -1)
    return xyz[idx], newpar.astype(np.int64), typ[idx]


def banded_recall(pred, truth):
    """Truth-cable fraction and recall@200 in radial bands from the soma."""
    from scipy.spatial import cKDTree
    w = G.tree_stats(truth)['weights']
    d = cKDTree(pred[0]).query(truth[0])[0]
    r = np.linalg.norm(truth[0] - truth[0][0], axis=1)
    out = []
    for lo, hi in zip(BANDS[:-1], BANDS[1:]):
        m = (r >= lo) & (r < hi)
        cable = float(w[m].sum())
        out.append({'band_um': [lo, None if np.isinf(hi) else hi], 'truth_weight': cable,
                    'covered_weight': float(w[m & (d <= RADIUS)].sum())})
    return out


def shifted(tree, soma):
    return G.shifted(tree, soma)


def per_query(job):
    qid, row, soma, forecast_paths, expected_digest = job
    truth = load_truth(row)
    if digest(truth) != expected_digest:
        raise RuntimeError(f'{qid}: rebuilt truth differs from SCORE_RECEIPT digest')
    res = {'query_id': qid, 'truth_extent_um': float(np.linalg.norm(truth[0] - truth[0][0], axis=1).max()),
           'truth_cable_um': G.tree_stats(truth)['cable_um'], 'trivial': {}, 'truncated_E1': {}, 'bands': {},
           'recheck': {}}
    res['trivial']['soma_only'] = G.geometry_metrics(
        (np.asarray([soma], float), np.asarray([-1], np.int64), np.asarray([1], np.int64)), truth)
    for L in STAR_LENGTHS:
        res['trivial'][f'star_{int(L)}'] = G.geometry_metrics(star(soma, L), truth)
    for m, p in forecast_paths.items():
        f = G.load_tree(p)
        res['recheck'][m] = G.geometry_metrics(f, truth)
        res['bands'][m] = banded_recall(f, truth)
        if m == 'E1_local_uniform':
            for r in TRUNCATE_RADII:
                res['truncated_E1'][f'r{int(r)}'] = G.geometry_metrics(truncate(f, r), truth)
    return res, truth


def copy_eval(job):
    """Copy other cohort truths onto this soma. Labelled ORACLE where truth-selected."""
    qid, truth_i, soma_i, cands = job
    out = []
    for cid, animal, tree in cands:
        out.append((cid, animal, G.geometry_metrics(shifted(tree, soma_i), truth_i)))
    return qid, out


def animal_mean(vals, animals):
    a = sorted(set(animals))
    return float(np.mean([np.mean([v for v, b in zip(vals, animals) if b == x]) for x in a]))


def main():
    results_sha = sha(BATCH / 'results.json')
    assert results_sha.startswith(EXPECTED_RESULTS_SHA), results_sha
    results = json.loads((BATCH / 'results.json').read_text())
    receipt = json.loads((OUT / 'SCORE_RECEIPT.json').read_text())
    decl = json.loads((BATCH / 'declaration.json').read_text())
    sealed = json.loads((BATCH / 'SEALED.json').read_text())
    rows = {r['neuron_id']: r for r in json.loads(COHORT.read_text())['records']}
    paths = {}
    for r in sealed['records']:
        rec = json.loads((BATCH / r['record_path']).read_text())
        paths.setdefault(r['query_id'], {})[r['method']] = BATCH / rec['tree_path']
    queries = decl['queries']
    animal = {q['id']: q['animal'] for q in queries}
    soma = {q['id']: q['soma_um'] for q in queries}
    jobs = [(q['id'], rows[q['id']], q['soma_um'], paths[q['id']], receipt['truth_cache_digests'][q['id']])
            for q in queries]
    print(f'{len(jobs)} queries; rebuilding truth and scoring trivial forecasts', flush=True)
    per, truths = {}, {}
    with ProcessPoolExecutor(6) as ex:
        for i, (res, tr) in enumerate(ex.map(per_query, jobs, chunksize=4), 1):
            per[res['query_id']] = res
            truths[res['query_id']] = tr
            if i % 50 == 0:
                print(f'  {i}/{len(jobs)}', flush=True)

    # recheck: recomputed metrics must equal the sealed results exactly
    sealed_rows = {(r['query_id'], r['method']): r for r in results['per_query']}
    worst = max(abs(per[q]['recheck'][m][k] - sealed_rows[(q, m)][k])
                for q in per for m in per[q]['recheck'] for k in ('distance_um', 'f1_200'))
    assert worst == 0.0, f'recheck mismatch {worst}'
    print('recheck: all 1572 recomputed distance/F1 values identical to sealed results', flush=True)

    # copying from other cohort truths (exposed data; diagnostics only)
    ids = [q['id'] for q in queries]
    S_xyz = np.asarray([soma[q] for q in ids])
    cjobs = []
    for i, q in enumerate(ids):
        d = np.linalg.norm(S_xyz - S_xyz[i], axis=1)
        d[i] = np.inf
        order = [j for j in np.argsort(d) if np.isfinite(d[j])]
        other = [j for j in order if animal[ids[j]] != animal[q]]
        pick = list(dict.fromkeys(order[:ORACLE_K] + other[:1]))
        cjobs.append((q, truths[q], soma[q], [(ids[j], animal[ids[j]], truths[ids[j]]) for j in pick]))
    print('scoring cohort-truth copies', flush=True)
    copies = {}
    with ProcessPoolExecutor(6) as ex:
        for q, out in ex.map(copy_eval, cjobs, chunksize=2):
            copies[q] = out

    A = [animal[q] for q in ids]
    methods = results['methods']
    table = {}

    def add(name, getter, note):
        dist = [getter(q)['distance_um'] for q in ids]
        f1 = [getter(q)['f1_200'] for q in ids]
        rec = [getter(q)['recall_200'] for q in ids]
        pre = [getter(q)['precision_200'] for q in ids]
        table[name] = {'distance_um': animal_mean(dist, A), 'f1_200': animal_mean(f1, A),
                       'recall_200': animal_mean(rec, A), 'precision_200': animal_mean(pre, A), 'note': note}

    for m in methods:
        add(m, lambda q, m=m: per[q]['recheck'][m], 'sealed forecast, recomputed (identical)')
    add('soma_only', lambda q: per[q]['trivial']['soma_only'], 'single point at the soma')
    for L in STAR_LENGTHS:
        add(f'star_{int(L)}', lambda q, L=L: per[q]['trivial'][f'star_{int(L)}'],
            f'{STAR_RAYS} straight rays of {int(L)} um; soma only; SWEEP, not a method')
    for r in TRUNCATE_RADII:
        add(f'E1_truncated_{int(r)}', lambda q, r=r: per[q]['truncated_E1'][f'r{int(r)}'],
            f'E1 forecast with everything beyond {int(r)} um of the soma removed')
    add('cohort_soma_nearest_any', lambda q: copies[q][0][2],
        'truth of the soma-nearest OTHER cohort neuron (any mouse), translated; exposed-data copy')
    add('cohort_soma_nearest_other_mouse',
        lambda q: next(m for c, a, m in copies[q] if a != animal[q]),
        'truth of the soma-nearest cohort neuron from a DIFFERENT mouse, translated')
    add('ORACLE_best_of_16_by_F1', lambda q: max((m for _, _, m in copies[q][:ORACLE_K]), key=lambda m: m['f1_200']),
        'ORACLE: best of the 16 soma-nearest cohort truths, selected using the query truth. Ceiling, not a method.')
    add('ORACLE_best_of_16_by_distance',
        lambda q: min((m for _, _, m in copies[q][:ORACLE_K]), key=lambda m: m['distance_um']),
        'ORACLE: as above, selected by distance.')

    # radial decomposition of truth cable and coverage, cohort-pooled per band, animal-equal fraction
    bands = {}
    for m in methods:
        rowsb = []
        for b in range(len(BANDS) - 1):
            tw = [per[q]['bands'][m][b]['truth_weight'] for q in ids]
            cw = [per[q]['bands'][m][b]['covered_weight'] for q in ids]
            tot = [sum(x['truth_weight'] for x in per[q]['bands'][m]) for q in ids]
            rowsb.append({'band_um': per[ids[0]]['bands'][m][b]['band_um'],
                          'truth_cable_share': animal_mean([t / s for t, s in zip(tw, tot)], A),
                          'recall_contribution': animal_mean([c / s for c, s in zip(cw, tot)], A),
                          'recall_within_band_pooled': float(sum(cw) / sum(tw)) if sum(tw) else None})
        bands[m] = rowsb

    # query-versus-method variance split on the six sealed methods
    def split(k):
        M = np.asarray([[sealed_rows[(q, m)][k] for m in methods] for q in ids])
        g = M.mean()
        ss_q = M.shape[1] * ((M.mean(1) - g) ** 2).sum()
        ss_m = M.shape[0] * ((M.mean(0) - g) ** 2).sum()
        ss_t = ((M - g) ** 2).sum()
        return {'query_share': float(ss_q / ss_t), 'method_share': float(ss_m / ss_t),
                'interaction_residual_share': float((ss_t - ss_q - ss_m) / ss_t)}

    e1_d = np.asarray([sealed_rows[(q, 'E1_local_uniform')]['distance_um'] for q in ids])
    ext = np.asarray([per[q]['truth_extent_um'] for q in ids])
    report = {
        'schema': 'synarch_axon_vacuity_diagnostic_v1',
        'status': 'POST_HOC_DIAGNOSTIC_ON_EXPOSED_DATA',
        'not_a_test': 'Exploratory arithmetic on the closed 262-query comparison. All 5 mice and the 262 truths '
                      'were already exposed. No result here is a method, a confirmation, or promotable. Sweeps '
                      'and oracles are labelled and must not be read as forecasts.',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'inputs': {'results_sha256': results_sha, 'score_receipt_sha256': sha(OUT / 'SCORE_RECEIPT.json'),
                   'adapter_sha256': sha(ROOT / 'analysis/axon_pretrain_truth_adapter.py'),
                   'geometry_sha256': sha(ROOT / 'analysis/brm_geometry_audit.py'),
                   'script_sha256': sha(__file__)},
        'integrity': {'truth_digests_match_score_receipt': True, 'recheck_max_abs_difference': worst,
                      'queries': len(ids), 'animals': len(set(A))},
        'aggregation': 'animal-equal means over 5 mice, same as the frozen summariser',
        'table': table,
        'radial_bands_truth_to_prediction': bands,
        'variance_split_over_six_sealed_methods': {'distance_um': split('distance_um'), 'f1_200': split('f1_200')},
        'E1_distance_vs_truth_extent_spearman': float(
            np.corrcoef(np.argsort(np.argsort(e1_d)), np.argsort(np.argsort(ext)))[0, 1]),
        'truth_extent_um_quantiles': dict(zip(('p10', 'p50', 'p90'), np.quantile(ext, [.1, .5, .9]).tolist())),
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    with REPORT.open('x') as f:
        json.dump(report, f, indent=1)
    print('report', REPORT, sha(REPORT))
    for k, v in table.items():
        print(f'{k:34s} dist {v["distance_um"]:9.1f}  F1 {v["f1_200"]:.4f}  R {v["recall_200"]:.4f}  P {v["precision_200"]:.4f}')
    print(json.dumps(report['variance_split_over_six_sealed_methods'], indent=1))
    print('E1 dist vs extent spearman', report['E1_distance_vs_truth_extent_spearman'])
    for m in ('E1_local_uniform', 'candidate'):
        print(m, [(b['band_um'], round(b['truth_cable_share'], 3), round(b['recall_contribution'], 4)) for b in bands[m]])


if __name__ == '__main__':
    main()
