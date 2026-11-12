"""Audit the saved bounded development run without fitting or rerunning any predictor.

Recomputes scores/selection from saved probabilities and checks source identities,
provenance, shrinkage, and the all-training deployment artifact.  The only call
into development code is its features-only deploy_predict using saved coefficients.
"""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / 'data/interim/paired_projection_logloss_development/logloss_dev_barseq2_a1_v1'
PROTOCOL = ROOT / 'docs/paired_projection_logloss_development_protocol.json'
CODE = ROOT / 'analysis/paired_projection_logloss_development.py'
FREEZE = '74bd3111ca7613bfa32a4d90c5e8acf9389caf67'
TOL = 1e-12


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def rows(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def git(*args):
    return subprocess.check_output(['git', '-C', str(ROOT), *args])


def replay():
    # Do not silently replace a completed audit.
    receipt_path = RUN / 'ROOT_REPLAY.json'
    if receipt_path.exists():
        raise RuntimeError('ROOT_REPLAY.json already exists; exclusive audit output required')
    complete = read_json(RUN / 'RUN_COMPLETE.json')
    started = read_json(RUN / 'STARTED.json')
    proto = read_json(PROTOCOL)
    pins = proto['data_status']['pins']
    mp = ROOT / pins['input_manifest']
    manifest = read_json(mp)
    summary = read_json(RUN / 'development_result.json')
    recipe = read_json(RUN / 'frozen_recipe_candidate.json')
    art = read_json(RUN / 'deployment_artifact_a1.json')
    paths = {'manifest': mp, 'code': CODE, 'protocol': PROTOCOL,
             'frozen_runner': ROOT / 'analysis/paired_projection_real_benchmark.py',
             'reference': ROOT / 'analysis/paired_projection_reference.py'}
    for key, path in paths.items():
        assert sha(path) == complete[key + '_sha256'] == started[key + '_sha256'], key
    for name in ('STARTED.json', 'development_result.json', 'frozen_recipe_candidate.json',
                 'deployment_artifact_a1.json', 'predictions.csv'):
        assert sha(RUN / name) == complete[name + '_sha256'], name
    for key, value in started.items():
        assert complete[key] == value, ('started/complete mismatch', key)
    assert complete['run_id'] == RUN.name == 'logloss_dev_barseq2_a1_v1'
    assert started['synthetic'] is False and started['exposure'] == 'EXPOSED_DEVELOPMENT_ONLY'
    assert proto['status'] == 'FROZEN_DEVELOPMENT_AUTHORIZED'
    assert pins['input_manifest_sha256'] == sha(mp)
    assert pins['frozen_runner_sha256'] == sha(paths['frozen_runner'])
    assert pins['reference_sha256'] == sha(paths['reference'])
    for path in (CODE, PROTOCOL):
        rel = str(path.relative_to(ROOT))
        assert hashlib.sha256(git('show', f'{FREEZE}:{rel}')).hexdigest() == sha(path), rel
    commit_time = datetime.fromisoformat(git('show', '-s', '--format=%cI', FREEZE).decode().strip())
    start_time = datetime.fromisoformat(started['started_utc'])
    assert commit_time < start_time < datetime.fromisoformat(complete['finished_utc'])
    ap = mp.parent / manifest['admission']['record_path']
    admission = read_json(ap)
    assert sha(ap) == pins['admission_sha256'] == manifest['admission']['record_sha256']
    assert admission['status'] == 'ADMITTED' and admission['dataset_id'] == manifest['dataset_id']
    tables = {}
    for name, pin in manifest['files'].items():
        path = mp.parent / pin['path']
        assert sha(path) == pin['sha256'] == admission['file_sha256'][name], name
        tables[name] = rows(path)
    assert sha(ROOT / manifest['protocol']['path']) == manifest['protocol']['sha256']

    cells = sorted(tables['cells'], key=lambda r: r['cell_id'])
    ids = [r['cell_id'] for r in cells]
    assert len(ids) == len(set(ids)) == 598
    index = {cid: i for i, cid in enumerate(ids)}
    groups = np.array([r['group_id'] for r in cells])
    animals = sorted(set(groups))
    targets = list(manifest['targets'])
    assert animals == ['XC95', 'XC97'] and len(targets) == 11
    counts = {a: int(np.sum(groups == a)) for a in animals}
    assert counts == {'XC95': 304, 'XC97': 294}
    assert summary['folds'] == [[a] for a in animals]
    mask = {(r['group_id'], r['target']): int(r['assayed']) for r in tables['mask']}
    assert len(mask) == len(tables['mask']) == len(animals) * len(targets)
    assert all(mask[(a, t)] == 1 for a in animals for t in targets)
    truth = {(r['cell_id'], r['target']): int(r['captured']) for r in tables['labels']}
    assert len(truth) == len(tables['labels']) == len(ids) * len(targets)
    assert all(int(r['qc_pass']) == 1 for r in tables['labels'])
    y = np.array([[truth[(cid, t)] for cid in ids] for t in targets])
    assert set(np.unique(y)) == {0, 1}
    cohort_sha = hashlib.sha256('\n'.join(ids).encode()).hexdigest()
    for target in targets:
        assert manifest['frozen_cohort'][target] == {'n': len(ids), 'sha256': cohort_sha}

    # Define the declared grid independently of the development runner.
    logit = ('candidate_logistic', 'location_logistic', 'technical')
    average = ('location_knn', 'molecular_knn', 'matched_knn')
    copy = ('copy_1nn_location', 'copy_1nn_matched')
    families = list(logit + average + copy + ('train_frequency', 'constant_0.5'))
    cs, ks, alphas = (0.01, 0.1, 1.0), (5, 15, 31), (0.0, 0.1, 0.25, 0.5)
    assert proto['grid']['C'] == complete['C_grid'] == list(cs)
    assert proto['grid']['k'] == complete['k_grid'] == list(ks)
    assert proto['grid']['alpha'] == complete['alpha_grid'] == list(alphas)
    assert complete['brier_guard_tol'] == 0.0
    configs, params, baselines = [], {}, {}
    for fam in families:
        if fam in logit:
            settings = [(f'{fam}|C={c}|alpha={a}', {'C': c, 'alpha': a}) for c in cs for a in alphas]
            baselines[fam] = f'{fam}|C=1.0|alpha=0.0'
        elif fam in average:
            settings = [(f'{fam}|k={k}|alpha={a}', {'k': k, 'alpha': a}) for k in ks for a in alphas]
            baselines[fam] = f'{fam}|k=15|alpha=0.0'
        elif fam in copy:
            settings = [(f'{fam}|alpha={a}', {'alpha': a}) for a in alphas]
            baselines[fam] = f'{fam}|alpha=0.0'
        else:
            settings = [(fam, {})]
            baselines[fam] = fam
        for key, par in settings:
            configs.append(key)
            params[key] = dict(family=fam, **par)
    assert len(configs) == 82 and summary['families'] == families
    assert set(summary['all_configurations']) == set(configs)
    ci = {k: i for i, k in enumerate(configs)}
    ti = {t: i for i, t in enumerate(targets)}
    fold = {a: i for i, a in enumerate(animals)}
    probabilities = np.empty((len(configs), len(targets), len(ids)))
    seen = np.zeros(probabilities.shape, dtype=bool)
    nrows = 0
    with (RUN / 'predictions.csv').open(newline='') as f:
        reader = csv.DictReader(f)
        assert reader.fieldnames == ['cell_id', 'group_id', 'endpoint', 'fold', 'config', 'probability', 'captured']
        for row in reader:
            cell, target, key = row['cell_id'], row['endpoint'], row['config']
            i, j, k = ci[key], ti[target], index[cell]
            assert not seen[i, j, k], ('duplicate saved prediction', cell, target, key)
            assert row['group_id'] == groups[k] and int(row['fold']) == fold[groups[k]]
            assert int(row['captured']) == truth[(cell, target)]
            p = float(row['probability'])
            assert np.isfinite(p) and 0 <= p <= 1
            probabilities[i, j, k] = p
            seen[i, j, k] = True
            nrows += 1
    assert seen.all() and nrows == len(configs) * len(targets) * len(ids) == 539396

    errors = {'score': [], 'shrinkage': [], 'artifact_scaler': [], 'artifact_q': [], 'deployment_predict': []}
    # Every shrinkage q is reconstructed from the OTHER animal's source labels.
    for target in targets:
        j = ti[target]
        for animal in animals:
            ev = groups == animal
            tr = ~ev
            q = float(np.mean(y[j, tr]))
            saved_fold = [r for r in summary['fold_log'] if r['endpoint'] == target and r['fold'] == fold[animal]]
            assert len(saved_fold) == 1
            fl = saved_fold[0]
            assert fl['eval_groups'] == [animal] and fl['n_train'] == int(tr.sum()) and fl['n_eval'] == int(ev.sum())
            errors['shrinkage'].append(abs(q - fl['q_fit']))
            for key, par in params.items():
                p = probabilities[ci[key], j, ev]
                if 'alpha' in par:
                    raw_key = key.rsplit('|alpha=', 1)[0] + '|alpha=0.0'
                    want = (1 - par['alpha']) * probabilities[ci[raw_key], j, ev] + par['alpha'] * q
                else:
                    want = np.full(ev.sum(), q if key == 'train_frequency' else 0.5)
                errors['shrinkage'].append(float(np.max(np.abs(p - want))))

    # Independent proper-score calculations, then the declared macro weighting.
    scores = {}
    for key in configs:
        p = probabilities[ci[key]]
        pc = np.clip(p, 1e-6, 1 - 1e-6)
        brier = (p - y) ** 2
        loss = -(y * np.log(pc) + (1 - y) * np.log(1 - pc))
        per_target = {}
        per_animal = {a: {} for a in animals}
        for target in targets:
            j = ti[target]
            pa = {a: {'n': counts[a], 'brier': float(np.mean(brier[j, groups == a])),
                      'log_loss': float(np.mean(loss[j, groups == a]))} for a in animals}
            per_target[target] = {'brier': float(np.mean([v['brier'] for v in pa.values()])),
                                  'log_loss': float(np.mean([v['log_loss'] for v in pa.values()])), 'per_animal': pa}
        for a in animals:
            per_animal[a] = {s: float(np.mean([per_target[t]['per_animal'][a][s] for t in targets]))
                             for s in ('brier', 'log_loss')}
        scores[key] = {'brier_target_equal': float(np.mean([per_target[t]['brier'] for t in targets])),
                       'log_loss_target_equal': float(np.mean([per_target[t]['log_loss'] for t in targets])),
                       'per_animal_macro': per_animal, 'per_target': per_target}

    def same_tree(want, got, label):
        if isinstance(want, dict):
            assert set(want) == set(got), label
            for key, value in want.items():
                same_tree(value, got[key], label + '/' + str(key))
        elif isinstance(want, (int, float)) and not isinstance(want, bool):
            assert isinstance(got, (int, float)) and np.isfinite(got), label
            delta = abs(want - got)
            errors['score'].append(float(delta))
            assert delta <= TOL, (label, want, got, delta)
        else:
            assert want == got, (label, want, got)

    same_tree(scores, summary['all_configurations'], 'scores')

    # Reapply the frozen guard and deterministic tie-rule on reconstructed scores.
    def order(key):
        p = params[key]
        return (scores[key]['log_loss_target_equal'], scores[key]['brier_target_equal'],
                -p.get('alpha', 0.0), p.get('C', 0.0), p.get('k', 0))

    chosen = {}
    for family in families:
        candidates = [key for key in configs if params[key]['family'] == family]
        base = baselines[family]
        valid = [key for key in candidates if scores[key]['brier_target_equal'] <= scores[base]['brier_target_equal']]
        selected = min(valid, key=order)
        unrestricted = min(candidates, key=order)
        chosen[family] = {'selected': selected, 'baseline': base, 'n_configs': len(candidates), 'n_pass_guard': len(valid),
                           'log_loss_target_equal': scores[selected]['log_loss_target_equal'],
                           'brier_target_equal': scores[selected]['brier_target_equal'],
                           'baseline_log_loss_target_equal': scores[base]['log_loss_target_equal'],
                           'baseline_brier_target_equal': scores[base]['brier_target_equal'],
                           'unrestricted_log_loss_best': unrestricted,
                           'guard_rejected_unrestricted_best': unrestricted != selected,
                           'unrestricted_log_loss_target_equal': scores[unrestricted]['log_loss_target_equal'],
                           'unrestricted_brier_target_equal': scores[unrestricted]['brier_target_equal']}
    same_tree(chosen, summary['selected'], 'selection')
    candidate = chosen['candidate_logistic']['selected']
    comparison = summary['comparison']
    assert comparison['candidate_config'] == candidate
    pairs = {}
    for family in families:
        if family == 'candidate_logistic':
            continue
        for label, key in (('tuned', chosen[family]['selected']), ('baseline', chosen[family]['baseline'])):
            pa = {a: {s: scores[key]['per_animal_macro'][a][s] - scores[candidate]['per_animal_macro'][a][s]
                      for s in ('brier', 'log_loss')} for a in animals}
            pairs[family + '|' + label] = {
                'control_config': key,
                'log_loss_target_equal_diff': scores[key]['log_loss_target_equal'] - scores[candidate]['log_loss_target_equal'],
                'brier_target_equal_diff': scores[key]['brier_target_equal'] - scores[candidate]['brier_target_equal'],
                'per_animal_macro_diff': pa,
                'animals_improved_log_loss': sum(v['log_loss'] > 0 for v in pa.values()),
                'animals_improved_brier': sum(v['brier'] > 0 for v in pa.values()),
                'targets_improved_log_loss': sum(scores[key]['per_target'][t]['log_loss'] > scores[candidate]['per_target'][t]['log_loss'] for t in targets),
                'targets_improved_brier': sum(scores[key]['per_target'][t]['brier'] > scores[candidate]['per_target'][t]['brier'] for t in targets),
                'n_targets': len(targets)}
    same_tree(pairs, comparison['pairs'], 'comparisons')

    # Audit the deployment artifact with source arrays; no parameter fitting.
    genes_by_id = {r['cell_id']: r for r in tables['genes']}
    assert len(genes_by_id) == len(tables['genes']) == len(ids) and set(genes_by_id) == set(ids)
    gene_names, pos_names = manifest['endogenous_genes'], manifest['position_columns']
    counts_raw = np.array([[float(genes_by_id[cid][g]) for g in gene_names] for cid in ids])
    positions = np.array([[float(r[p]) for p in pos_names] for r in cells])
    X = np.hstack([np.log1p(counts_raw), positions])
    assert art['gene_features'] == gene_names and art['position_features'] == pos_names
    assert art['features'] == gene_names + pos_names and art['gene_transform'] == 'log1p'
    assert art['training_animals'] == animals and set(art['targets']) == set(targets)
    assert art['setting'] == candidate and art['C'] == recipe['C'] == params[candidate]['C']
    assert art['alpha'] == recipe['alpha'] == params[candidate]['alpha']
    assert recipe['controls_frozen_with_it'] == {f: chosen[f]['selected'] for f in families}
    mu, scale = X.mean(axis=0), X.std(axis=0)
    scale[scale == 0] = 1.0
    independent = {}
    for target in targets:
        t = art['targets'][target]
        assert t['n_fit'] == len(ids) and t['n_fit_per_animal'] == counts
        errors['artifact_scaler'].extend(np.abs(mu - t['scaler_mean']).tolist())
        errors['artifact_scaler'].extend(np.abs(scale - t['scaler_scale']).tolist())
        q = float(np.mean([np.mean(y[ti[target], groups == a]) for a in animals]))
        errors['artifact_q'].append(abs(q - t['q_fit']))
        assert np.asarray(t['coef']).shape == (X.shape[1],) and np.isfinite(t['coef']).all()
        assert np.isfinite(t['intercept'])
        z = ((X - np.asarray(t['scaler_mean'])) / np.asarray(t['scaler_scale'])) @ np.asarray(t['coef']) + t['intercept']
        independent[target] = (1 - art['alpha']) * expit(z) + art['alpha'] * t['q_fit']
    # The only production-code call: inference with saved coefficients, never fit.
    sys.path.insert(0, str(ROOT))
    from analysis import paired_projection_logloss_development as D
    predicted = D.deploy_predict(art, gene_names, pos_names, counts_raw, positions)
    assert set(predicted) == set(independent)
    for target in targets:
        assert predicted[target].shape == (len(ids),) and np.isfinite(predicted[target]).all()
        errors['deployment_predict'].append(float(np.max(np.abs(predicted[target] - independent[target]))))
    maxima = {key: float(max(values, default=0.0)) for key, values in errors.items()}
    assert all(v <= TOL for v in maxima.values()), maxima
    record = {'schema': 'synarch_logloss_development_saved_prediction_replay_v1', 'decision': 'PASS',
              'replay_code_sha256': sha(__file__), 'freeze_commit': FREEZE,
              'run_complete_sha256': sha(RUN / 'RUN_COMPLETE.json'),
              'started_sha256': sha(RUN / 'STARTED.json'), 'manifest_sha256': sha(mp),
              'development_code_sha256': sha(CODE), 'protocol_sha256': sha(PROTOCOL),
              'development_result_sha256': sha(RUN / 'development_result.json'),
              'predictions_sha256': sha(RUN / 'predictions.csv'),
              'recipe_sha256': sha(RUN / 'frozen_recipe_candidate.json'),
              'deployment_artifact_sha256': sha(RUN / 'deployment_artifact_a1.json'),
              'prediction_rows': nrows, 'configurations': len(configs), 'cells': len(ids),
              'animals': counts, 'targets': targets, 'score_tolerance': TOL,
              'maximum_absolute_differences': maxima,
              'exact_source_ids_labels_groups_folds_and_configuration_cohorts': True,
              'all_configuration_scores_and_all_selections_verified': True,
              'all_comparisons_verified': True, 'fit_only_shrinkage_verified': True,
              'deployment_all_training_scalers_frequencies_counts_and_inference_verified': True,
              'selected_settings': {f: chosen[f]['selected'] for f in families},
              'selected_macro_scores': {f: {s: scores[chosen[f]['selected']][s]
                                           for s in ('brier_target_equal', 'log_loss_target_equal')} for f in families},
              'predictor_fits_performed_by_this_check': 0,
              'scope': 'EXPOSED DEVELOPMENT ONLY; selected scores are optimistic. Saved-prediction replay and '
                       'saved-coefficient inference only; no new test, predictor fit, development rerun, or promotion.',
              'checked_utc': datetime.now(timezone.utc).isoformat()}
    with receipt_path.open('x') as f:
        json.dump(record, f, indent=2, sort_keys=True, allow_nan=False)
    print(json.dumps({'decision': record['decision'], 'prediction_rows': nrows, 'configurations': len(configs),
                      'maximum_absolute_differences': maxima, 'selected_candidate': candidate,
                      'receipt': str(receipt_path), 'receipt_sha256': sha(receipt_path)}, sort_keys=True))
    return record


if __name__ == '__main__':
    replay()
