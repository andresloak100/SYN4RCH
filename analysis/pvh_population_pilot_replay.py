"""Independent no-refit replay of the exploratory PVH population pilot.

Reads only completed extraction/atlas/forecast outputs. Does not open RNA RDS
objects, select models or refit anything. Collection labels are not absence calls.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data/interim/pvh_population_pilot'
META = ROOT / 'data/interim/paired_projection_feasibility/pvh_li2026'
ALPHAS = {'logistic_frozen_recipe': .1, 'hgbt_fixed': .1, 'molecular_knn31': .1,
          'pool_weighted_knn31': .1, 'gene_copy_1nn': .5, 'train_prior': 0., 'technical_logistic': .1}
TOL = 1e-12
ATLAS_SOURCE_COMMIT = 'c85cc7324ed93b5e412c2fc2c2c6ec97f41eeab7'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def rows(path, delimiter=','):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def near(expected, actual, name):
    a, b = np.asarray(expected), np.asarray(actual)
    require(a.shape == b.shape and np.isfinite(b).all() and
            np.allclose(a, b, atol=TOL, rtol=0), name + ' differs')


def historical_blob(commit, relative_path):
    result = subprocess.run(['git', '-C', str(ROOT), 'show', f'{commit}:{relative_path}'],
                            capture_output=True, check=True)
    return result.stdout


def verify_run(folder, protocol_path, historical_commit=None):
    complete = load(folder / 'RUN_COMPLETE.json')
    for name, digest in complete['files_sha256'].items():
        p = Path(name)
        require(not p.is_absolute() and '..' not in p.parts, 'unsafe output path')
        require(sha(folder / name) == digest, f'output hash differs: {folder.name}/{name}')
    started = load(folder / 'STARTED.json')
    require(started['protocol_snapshot_sha256'] == sha(folder / 'protocol_snapshot.json'), 'protocol snapshot binding differs')
    if historical_commit is None:
        require(sha(folder / 'protocol_snapshot.json') == sha(protocol_path), 'current protocol differs')
        require(started['source_sha256'] == sha(ROOT / 'analysis/pvh_population_pilot.py'), 'primary source differs')
        require(started['extractor_sha256'] == sha(ROOT / 'analysis/pvh_population_extract.R'), 'extractor differs')
    else:
        for field, name in (('source_sha256', 'analysis/pvh_population_pilot.py'),
                            ('extractor_sha256', 'analysis/pvh_population_extract.R'),
                            ('protocol_snapshot_sha256', 'docs/pvh_population_pilot_protocol.json')):
            digest = hashlib.sha256(historical_blob(historical_commit, name)).hexdigest()
            require(started[field] == digest, 'historical atlas source/protocol binding differs: ' + name)
    return complete, started


def source_cohort(protocol):
    raw = ROOT / 'data/raw/paired_projection_feasibility/pvh_li2026'
    for name, pin in protocol['source']['objects'].items():
        require(sha(raw / name) == pin['sha256'], 'source RDS byte hash differs: ' + name)
    for name, digest in protocol['source']['metadata'].items():
        require(sha(META / name) == digest, 'source metadata differs: ' + name)
    for rel, digest in protocol['imports_sha256'].items():
        require(sha(ROOT / rel) == digest, 'imported source differs: ' + rel)
    records = load(META / 'paired_nucleus_metadata_manifest.json')
    wanted = [r for r in records if r['original_study'] is True and r['exact_first_MUX_barcode_join'] is True]
    cohort = [{'object': r['observed_positive_target'], 'cell_id': r['cell_id'], 'pool': r['encoded_pool_id'],
               'sex': r['sex_from_encoded_pool'], 'target': r['observed_positive_target']} for r in wanted]
    cohort.sort(key=lambda r: (r['object'], r['cell_id']))
    require(len(cohort) == 265 and len({(r['object'], r['cell_id']) for r in cohort}) == 265, 'cohort IDs differ')
    require(dict(Counter(r['pool'] for r in cohort)) == protocol['fit_cohort']['expected_pools'], 'pool counts differ')
    require(all(r['individual_mouse_id'] is None for r in wanted), 'unexpected individual animal IDs')
    return records, cohort


def verify_atlas(folder, records):
    actual = rows(folder / 'atlas_nuclei.csv')
    original = {(r['observed_positive_target'], r['cell_id']): r for r in records if r['original_study'] is True}
    excluded = [r for r in records if r['original_study'] is not True]
    require(len(original) == len(actual) == 1225 and len(excluded) == 2498, 'atlas original/reused counts differ')
    meta = {obj: {r['cell_id']: r for r in rows(META / f'{obj}_cell_metadata.csv')} for obj in ('PB', 'SC')}
    seen = set()
    grouped = {'author_inferred_type': Counter(), 'encoded_pool_id': Counter(), 'object_batch': Counter()}
    for r in actual:
        key = r['object'], r['cell_id']
        require(key in original and key not in seen, 'atlas missing/duplicated source nucleus')
        seen.add(key); src = original[key]
        require(r['nucleus_key'] == '|'.join(key), 'atlas nucleus key differs')
        require(r['positive_collection_target'] == key[0] and r['source_study'] == src['source_study'], 'atlas source differs')
        require(r['author_inferred_type'] == meta[key[0]][key[1]]['predicted.celltype'], 'author inferred type differs')
        require(r['encoded_pool_id'] == (src['encoded_pool_id'] or '') and r['object_batch'] == src['object_batch'],
                'atlas pool/batch differs')
        require(all(r[k] == '' for k in ('individual_animal_id', 'measured_position', 'other_target_status')),
                'atlas invented measurement or other-target status')
        for group in grouped:
            grouped[group][r[group], key[0]] += 1
    require(seen == set(original), 'atlas source key set differs')
    count_rows = rows(folder / 'atlas_counts.csv')
    expected = {(group, g, t): n for group, counts in grouped.items() for (g, t), n in counts.items()}
    observed = {(r['grouping'], r['group'], r['positive_collection_target']): int(r['n_nuclei']) for r in count_rows}
    require(expected == observed and len(observed) == len(count_rows), 'atlas grouped counts differ')
    graph = load(folder / 'evidence_graph.json')
    nodes = {n['id']: n for n in graph['nodes']}
    require(len(nodes) == len(graph['nodes']), 'graph node IDs duplicated')
    require(all(n.get('coordinates') is None for n in nodes.values()), 'graph invented coordinates')
    types = [n for n in nodes.values() if n['kind'] == 'author_inferred_type']
    require(len(types) == 12 and all(n['inferred'] is True for n in types), 'inferred types treated as measurements')
    edges = [e for e in graph['edges'] if e['kind'] == 'collected_in']
    target_n = Counter(r['positive_collection_target'] for r in actual)
    want_edges = {(f'type:{g}', f'target:{t}'): n for (g, t), n in grouped['author_inferred_type'].items()}
    require({(e['source'], e['target']): e['n_nuclei'] for e in edges} == want_edges, 'graph evidence counts differ')
    for e in edges:
        t = e['target'].split(':', 1)[1]
        near(e['n_nuclei'] / target_n[t], e['fraction_of_target_collection'], 'sampling fraction')
        require(e['other_target_status'] is None and 'not prevalence' in e['meaning'], 'graph overstated edge meaning')
    require(graph['unknown_pool_nuclei'] == sum(r['encoded_pool_id'] == '' for r in actual), 'unknown pool count differs')
    return {'original_nuclei': 1225, 'reused_nuclei_excluded': 2498, 'author_inferred_types': 12,
            'positive_sampling_counts': dict(target_n), 'individual_animals_or_positions_invented': False}


def excluded_gene(name):
    low = name.casefold(); tokens = re.split(r'[^a-z0-9]+', low)
    aliases = set('gfp egfp yfp eyfp cfp ecfp rfp tdtomato tomato dtomato mcherry cherry mscarlet scarlet tdt '
                  'cre icre flp flpo dre wpre sv40 bgh tva rabies aav raav retroaav cav2 hsv dsred bfp ebfp '
                  'ebfp2 mrfp mrfp1 totalseq hto'.split())
    return bool(set(tokens) & aliases or re.match(r'^ercc[-_]', low) or any(
        word in low for word in ('target', 'label', 'captur', 'class', 'cluster', 'projection', 'barcode',
                                 'reporter', 'vector', 'qc_pass', 'assay')) or low in {'cell_id', 'group_id'})


def extraction_counts(folder, cohort):
    require(rows(folder / 'cohort_full.csv') == cohort, 'extracted cohort differs from pinned source')
    per = {}
    for obj in ('PB', 'SC'):
        base = folder / obj; genes = (base / 'genes.txt').read_text().splitlines(); cs = rows(base / 'cells.tsv', '\t')
        ids = [r['cell_id'] for r in cs]
        require(len(set(genes)) == len(genes) and len(set(ids)) == len(ids), 'export feature/cell IDs duplicated')
        wanted = [r['cell_id'] for r in cohort if r['object'] == obj]
        require(ids == wanted, 'exported cell order/set differs')
        trip = np.loadtxt(base / 'triplets.tsv', skiprows=1, ndmin=2)
        require(trip.shape[1] == 3 and np.isfinite(trip).all() and (trip == np.round(trip)).all(), 'invalid sparse counts')
        gi, ci, val = trip[:, 0].astype(int)-1, trip[:, 1].astype(int)-1, trip[:, 2]
        require((gi >= 0).all() and (gi < len(genes)).all() and (ci >= 0).all() and (ci < len(ids)).all()
                and (val >= 0).all(), 'sparse index/count out of range')
        mat = np.zeros((len(ids), len(genes))); np.add.at(mat, (ci, gi), val)
        rec = dict(line.split('\t', 1) for line in (base / 'receipt.tsv').read_text().splitlines())
        require(rec['assay'] == 'RNA' and rec['layer'] == 'counts' and rec['packages_loaded_by_script'] == 'none',
                'wrong RNA layer or author/package code used')
        per[obj] = genes, ids, mat
    common = sorted(set(per['PB'][0]) & set(per['SC'][0])); excluded = [g for g in common if excluded_gene(g)]
    panel = [g for g in common if g not in set(excluded)]
    matrices = {}
    for obj, (genes, ids, mat) in per.items():
        ix = {g: j for j, g in enumerate(genes)}
        matrices[obj] = {cell: mat[i, [ix[g] for g in panel]] for i, cell in enumerate(ids)}
    counts = np.vstack([matrices[r['object']][r['cell_id']] for r in cohort])
    require((counts.sum(axis=1) > 0).all(), 'zero-library nucleus')
    audit = load(folder / 'panel.json')
    require(audit['common_genes'] == len(panel) and audit['excluded_gene_ids'] == excluded and
            audit['panel_sha256'] == hashlib.sha256('\n'.join(panel).encode()).hexdigest(), 'panel receipt differs')
    return counts, panel, excluded


def scores(forecasts, cohort, folds):
    out = {}
    for method in ALPHAS:
        by_fold = {}
        for fold, spec in folds.items():
            pool_scores = {}
            for pool in spec['test']:
                yy, pp = [], []
                for r in cohort:
                    if r['pool'] == pool:
                        yy.append(r['target'] == 'PB'); pp.append(forecasts[fold, r['object'], r['cell_id'], method]['p'])
                y, p = np.asarray(yy, float), np.asarray(pp); q = np.clip(p, 1e-6, 1-1e-6)
                pool_scores[pool] = {'n': len(y), 'brier': float(np.mean((y-p)**2)),
                                     'log_loss': float(np.mean(-y*np.log(q)-(1-y)*np.log1p(-q)))}
            by_fold[fold] = {'per_pool': pool_scores, **{
                m: float(np.mean([v[m] for v in pool_scores.values()])) for m in ('brier', 'log_loss')}}
        out[method] = {'per_fold': by_fold, **{m+'_macro': float(np.mean([v[m] for v in by_fold.values()]))
                                            for m in ('brier', 'log_loss')}}
    return out


def verify_fit(folder, cohort, counts, panel, protocol, done):
    folds = {k: v for k, v in protocol['folds'].items() if k.startswith('fold_')}
    forecasts = {}
    expected = {}
    for fold, spec in folds.items():
        require(not set(spec['train']) & set(spec['test']), 'pool overlap')
        for r in cohort:
            if r['pool'] in spec['test']:
                for method in ALPHAS:
                    expected[fold, r['object'], r['cell_id'], method] = r
    for r in rows(folder / 'predictions.csv'):
        key = r['fold'], r['object'], r['cell_id'], r['method']
        require(key in expected and key not in forecasts, 'unexpected or duplicate prediction key')
        src = expected[key]
        require((r['pool'], r['sex'], r['collection_target']) == (src['pool'], src['sex'], src['target']), 'forecast identity differs')
        p, raw = float(r['probability']), float(r['raw_probability'])
        require(np.isfinite([p, raw]).all() and 0 <= p <= 1 and 0 <= raw <= 1, 'invalid forecast probability')
        forecasts[key] = {'p': p, 'raw': raw}
    require(set(forecasts) == set(expected) and len(forecasts) == 1855, 'forecasts missing')
    rebuilt = scores(forecasts, cohort, folds); reported = load(folder / 'metrics.json')
    for method, sc in rebuilt.items():
        for metric in ('brier', 'log_loss'):
            near(sc[metric+'_macro'], reported['per_method'][method][metric+'_macro'], 'macro '+metric)
        for fold, fs in sc['per_fold'].items():
            for metric in ('brier', 'log_loss'):
                near(fs[metric], reported['per_method'][method]['per_fold'][fold][metric], 'fold '+metric)
                for pool, ps in fs['per_pool'].items():
                    near(ps[metric], reported['per_method'][method]['per_fold'][fold]['per_pool'][pool][metric], 'pool '+metric)
    import joblib
    sys.path.insert(0, str(ROOT))
    from analysis import pvh_population_pilot as primary
    pools = np.array([r['pool'] for r in cohort]); y = np.array([r['target'] == 'PB' for r in cohort], int)
    lx = np.log1p(counts/counts.sum(axis=1, keepdims=True)*10000)
    technical = np.column_stack([np.log1p(counts.sum(axis=1)), np.log1p((counts > 0).sum(axis=1))])
    worst = 0.; details = {}
    for fold, spec in folds.items():
        tr, te = np.isin(pools, spec['train']), np.isin(pools, spec['test'])
        audit = load(folder / f'fold_{fold}.json'); variance = lx[tr].var(axis=0)
        selected = sorted([j for j, v in enumerate(variance) if v > 0], key=lambda j: (-variance[j], panel[j]))[:500]
        names = [panel[j] for j in selected]; mu = lx[tr][:, selected].mean(axis=0); sd = lx[tr][:, selected].std(axis=0)
        sd[sd == 0] = 1
        require(audit['selected_genes'] == names, 'TRAIN-only selected genes differ')
        near(mu, audit['scaler_mean'], 'training gene mean'); near(sd, audit['scaler_scale'], 'training gene scale')
        train_keys = [f"{cohort[i]['object']}|{cohort[i]['cell_id']}" for i in np.flatnonzero(tr)]
        require(audit['train_cells'] == train_keys and audit['train_pools'] == spec['train'] and
                audit['test_pools'] == spec['test'] and audit['alphas'] == ALPHAS, 'fold fit identity/alphas differ')
        q = float(np.mean([np.mean(y[tr & (pools == pool)]) for pool in spec['train']]))
        near(.5, q, 'balanced collection prior'); near(q, audit['q_pool_balanced'], 'saved TRAIN prior')
        model_name = f'models_{fold}.joblib'; require(model_name in done['files_sha256'], 'model not hash-pinned')
        artifact = joblib.load(folder / model_name); state = artifact['state']
        require(artifact['panel'] == panel and artifact['train_cell_keys'] == train_keys and
                artifact['features']['selected_index'] == selected, 'artifact feature/fit identity differs')
        near(state['q'], q, 'artifact prior'); near(state['ytr'], y[tr], 'artifact TRAIN labels')
        inverse = np.array([1/np.sum(pools[tr] == p) for p in pools[tr]]); inverse /= inverse.mean()
        near(inverse, state['pool_weights'], 'training pool weights')
        x = (lx[:, selected]-mu)/sd
        near(x[tr], state['knn31']._fit_X, 'stored KNN training features'); near(x[tr], state['nn1']._fit_X, 'copy training features')
        near(technical[tr].mean(axis=0), state['tech_scaler'].mean_, 'technical TRAIN mean')
        hs = state['hgbt_adapter']; near(np.median(x[tr], axis=0), hs['median'], 'HGBT TRAIN median')
        near(x[tr].mean(axis=0), hs['scaler_mean'], 'HGBT TRAIN mean'); near(x[tr].std(axis=0), hs['scaler_scale'], 'HGBT TRAIN scale')
        z = (x[te]-np.asarray(hs['scaler_mean']))/np.asarray(hs['scaler_scale'])
        donors = state['knn31'].kneighbors(x[te], return_distance=False); ww = inverse[donors]
        raw = {'logistic_frozen_recipe': state['logistic'].predict_proba(x[te])[:, 1],
               'hgbt_fixed': state['hgbt'].predict_proba(z)[:, 1], 'molecular_knn31': state['knn31'].predict_proba(x[te])[:, 1],
               'pool_weighted_knn31': (ww*y[tr][donors]).sum(axis=1)/ww.sum(axis=1),
               'gene_copy_1nn': y[tr][state['nn1'].kneighbors(x[te], return_distance=False)[:, 0]].astype(float),
               'train_prior': np.full(te.sum(), q),
               'technical_logistic': state['technical'].predict_proba(state['tech_scaler'].transform(technical[te]))[:, 1]}
        api_raw, api_prob = primary.artifact_predict(artifact, panel, counts[te])
        for m, rp in raw.items():
            prob = (1-ALPHAS[m])*rp+ALPHAS[m]*q
            saved = [forecasts[fold, cohort[i]['object'], cohort[i]['cell_id'], m] for i in np.flatnonzero(te)]
            near(rp, [v['raw'] for v in saved], 'saved raw inference '+m); near(prob, [v['p'] for v in saved], 'saved shrunk inference '+m)
            near(rp, api_raw[m], 'features-only API raw '+m); near(prob, api_prob[m], 'features-only API probability '+m)
            worst = max(worst, float(np.max(np.abs(prob-api_prob[m]))))
        details[fold] = {'training_cells': int(tr.sum()), 'evaluation_cells': int(te.sum()), 'selected_genes': len(selected), 'prior': q}
    return {'prediction_rows': 1855, 'per_method': rebuilt, 'fold_checks': details,
            'artifact_inference_max_abs_difference': worst, 'refits': 0}


def replay(base=BASE):
    pp = ROOT / 'docs/pvh_population_pilot_protocol.json'; protocol = load(pp)
    atlas = base / protocol['atlas']['run_id']; extract = base / protocol['extraction']['run_id']
    fit = base / protocol['fit_run']['run_id']
    ac, _ = verify_run(atlas, pp, historical_commit=ATLAS_SOURCE_COMMIT)
    ec, _ = verify_run(extract, pp); fc, started = verify_run(fit, pp)
    recovery = None
    if extract.name == 'pvh_rna_extract_v2':
        failed = base / 'pvh_rna_extract_v1'
        require((failed / 'STARTED.json').is_file() and not (failed / 'RUN_COMPLETE.json').exists(),
                'failed extraction v1 not preserved as incomplete')
        require(not (failed / 'PB').exists() and not (failed / 'SC').exists(), 'failed v1 unexpectedly wrote RNA exports')
        recovery = {'failed_run': failed.name, 'failed_started_sha256': sha(failed / 'STARTED.json'),
                    'completed_run': extract.name, 'failed_run_preserved': True,
                    'meaning': 'Selected-cohort integer check replaces overbroad non-cohort guard before any model fit.'}
    require(started['extraction_run_complete_sha256'] == sha(extract / 'RUN_COMPLETE.json'), 'extraction completion pin differs')
    records, cohort = source_cohort(protocol); atlas_result = verify_atlas(atlas, records)
    counts, panel, excluded = extraction_counts(extract, cohort)
    fit_result = verify_fit(fit, cohort, counts, panel, protocol, fc)
    return {'status': 'PASS', 'classification': 'EXPLORATORY COLLECTION-LABEL POOL TRANSFER; NOT INDIVIDUAL ANIMAL OR CIRCUIT VALIDATION',
            'atlas': atlas_result, 'fit': fit_result, 'n_common_genes': len(panel), 'excluded_gene_ids': excluded,
            'source_sha256': sha(ROOT / 'analysis/pvh_population_pilot.py'), 'protocol_sha256': sha(pp),
            'atlas_historical_source_commit': ATLAS_SOURCE_COMMIT, 'extraction_recovery': recovery,
            'completion_sha256': {p.name: sha(p / 'RUN_COMPLETE.json') for p in (atlas, extract, fit)},
            'replay_source_sha256': sha(__file__), 'no_raw_RDS_expression_read': True, 'refits_performed': 0}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path)
    args = parser.parse_args(); receipt = replay(); output = args.output or BASE / 'pvh_population_pilot_v1/ROOT_REPLAY.json'
    with output.open('x') as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True, allow_nan=False); stream.write('\n')
    print(json.dumps({'status': 'PASS', 'receipt': str(output), 'sha256': sha(output), 'refits': 0}))
