"""Independent, no-refit replay of nonlinear A1 development and fresh simulations.

Does not import the primary runner, its selectors or its metric routines. Reads
saved forecasts and audits only; writes a new receipt on explicit CLI request.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CLIP = (1e-6, 1 - 1e-6)
TOL = 1e-12
FAMILIES = ('hgbt_gene_position', 'hgbt_position')
ALPHAS = (0.0, 0.1, 0.25)
NEW_SEEDS = {'location_only': (876211, 876212, 876213),
             'complementarity': (876221, 876222, 876223),
             'latent_unobserved': (876231, 876232, 876233),
             'batch_shift_missing': (876241, 876242, 876243)}
FROZEN_LR = 'candidate_logistic|C=0.1|alpha=0.1'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def close_tree(expected, actual, path='root'):
    """Verify each independently reconstructed field, allowing metadata extras."""
    if isinstance(expected, dict):
        require(isinstance(actual, dict), path + ': expected mapping')
        for k, v in expected.items():
            require(k in actual, path + ': missing ' + k)
            close_tree(v, actual[k], path + '.' + k)
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(expected) == len(actual), path + ': list differs')
        for i, (e, a) in enumerate(zip(expected, actual)):
            close_tree(e, a, path + f'[{i}]')
    elif isinstance(expected, bool) or expected is None or isinstance(expected, str):
        require(type(expected) is type(actual) and expected == actual, path + ': value differs')
    elif isinstance(expected, (float, int)):
        require(isinstance(actual, (float, int)) and not isinstance(actual, bool), path + ': not numeric')
        require(np.isfinite(actual) and abs(expected - actual) <= TOL, path + ': number differs')
    else:
        require(expected == actual, path + ': value differs')


def verify_completion(run_dir):
    done = load(run_dir / 'RUN_COMPLETE.json')
    for name, digest in done['files_sha256'].items():
        require(Path(name).name == name, 'output hash path must be a basename')
        require(sha(run_dir / name) == digest, f'output hash differs: {name}')
    require('predictions.csv' in done['files_sha256'], 'forecasts not covered by completion receipt')
    if 'protocol_snapshot.json' in done['files_sha256']:
        started = load(run_dir/'STARTED.json')
        require(sha(run_dir/'protocol_snapshot.json') == started['protocol_snapshot_sha256'],
                'protocol snapshot differs from STARTED')
        require(started['protocol_sha256'] == started['protocol_snapshot_sha256'],
                'executing protocol and saved protocol differ')
    return done


def adapter_from_train(x):
    x = np.asarray(x, float)
    med = np.array([np.median(v[np.isfinite(v)]) if np.isfinite(v).any() else 0. for v in x.T])
    filled = np.where(np.isfinite(x), x, med)
    scale = np.std(filled, axis=0)
    scale[scale == 0] = 1.
    return {'median': med.tolist(), 'scaler_mean': np.mean(filled, axis=0).tolist(),
            'scaler_scale': scale.tolist()}


def check_training_audit(audit, x, y, groups):
    """Recompute eligibility, group prior and adapter without fitting any predictor."""
    groups = np.asarray(groups)
    unique = sorted(set(map(str, groups)))
    q = float(np.mean([np.mean(y[groups == group]) for group in unique]))
    expected = {'q_fit': q, 'n_train': len(y), 'n_train_positive': int(np.sum(y)),
                'train_groups': unique, 'adapter': adapter_from_train(x)}
    close_tree(expected, audit)
    return abs(q-audit['q_fit'])


def read_table(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def a1_inputs(protocol):
    mp = ROOT/protocol['pins']['a1_input_manifest']
    require(sha(mp) == protocol['pins']['a1_input_manifest_sha256'], 'A1 input manifest pin differs')
    manifest = load(mp)
    tables = {}
    for name, pin in manifest['files'].items():
        path = mp.parent/pin['path']
        require(sha(path) == pin['sha256'], 'A1 source table pin differs: '+name)
        tables[name] = read_table(path)
    cells = tables['cells']; genes = {r['cell_id']: r for r in tables['genes']}
    ids = [r['cell_id'] for r in cells]
    require(len(ids) == len(set(ids)) == 598, 'A1 input cells differ')
    panel, positions = manifest['endogenous_genes'], manifest['position_columns']
    g = np.array([[float(genes[cell][name]) for name in panel] for cell in ids])
    pos = np.array([[float(row[name]) for name in positions] for row in cells])
    labels = {(r['cell_id'], r['target']): int(r['captured']) for r in tables['labels']}
    group = np.asarray([r['group_id'] for r in cells])
    return manifest, ids, g, pos, group, labels


def verify_a1_training(data, audit, protocol):
    manifest, ids, genes, pos, groups, labels = a1_inputs(protocol)
    views = {FAMILIES[0]: np.column_stack([np.log1p(genes), pos]), FAMILIES[1]: pos}
    require(set(data) == set(labels), 'A1 prediction cell/target keys differ from source labels')
    index = {cell: i for i, cell in enumerate(ids)}
    for key, obs in data.items():
        require(obs['y'] == labels[key] and obs['group'] == groups[index[key[0]]], 'A1 source labels or animals differ')
    worst = 0.
    for key, families in audit['per_endpoint_fold'].items():
        target, fold = key.rsplit('|', 1)
        held = audit['folds'][int(fold)]
        tr = ~np.isin(groups, held)
        y = np.array([labels[(cell, target)] for cell in ids])
        tr &= y >= 0
        require(not set(groups[tr]) & set(held), 'A1 train/test animal overlap')
        for family in FAMILIES:
            worst = max(worst, check_training_audit(families[family], views[family][tr], y[tr], groups[tr]))
    return {'training_prior_max_abs_error': worst, 'source_labels_checked': len(labels),
            'fold_adapters_and_priors_recomputed': True}


def verify_sim_training(run_dir, data, audit, case):
    scenario, seed = case.split('|')
    path = run_dir/f'SYNTHETIC_{scenario}_{seed}.npz'
    with np.load(path, allow_pickle=False) as z:
        d = {k: z[k] for k in z.files}
    require(len(d['cell_id']) == 1000 and len(set(d['cell_id'])) == 1000, 'synthetic input cell IDs differ')
    require(int(d['seed']) == int(seed) and set(d['case_id']) == {case}, 'NPZ case identity differs')
    tr, ev = d['split'] == 'train', d['split'] == 'test'
    require(int(tr.sum()) == 600 and int(ev.sum()) == 400 and int((tr | ev).sum()) == 1000,
            'NPZ split row counts differ')
    tg, eg = sorted(set(d['animal'][tr])), sorted(set(d['animal'][ev]))
    require(len(tg) == 6 and len(eg) == 4 and not set(tg) & set(eg), 'NPZ animal splits differ')
    require(all(np.sum(d['animal'] == a) == 100 for a in tg+eg), 'synthetic cells per animal differ')
    close_tree({'train_animals': tg, 'test_animals': eg, 'test_labels_read_by_fitter': False}, audit['comparators'])
    expected = {(str(d['cell_id'][i]), str(t)): (i, ti) for i in np.flatnonzero(ev)
                for ti, t in enumerate(d['target_names'])}
    require(set(data) == set(expected), 'NPZ and forecast key sets differ')
    for key, obs in data.items():
        i, ti = expected[key]
        require(obs['y'] == d['y'][i, ti] and obs['assayed'] == bool(d['assayed'][i, ti]) and
                obs['group'] == d['animal'][i], 'NPZ label/assay/animal differs from forecast')
    genes = np.log1p(d['genes'])
    views = {FAMILIES[0]: np.column_stack([genes, d['xyz']]), FAMILIES[1]: d['xyz']}
    comp_views = {'candidate_logistic': views[FAMILIES[0]], 'location_logistic': d['xyz'],
                  'molecular': genes, 'technical': d['tech']}
    worst = 0.
    for ti, target in enumerate(d['target_names']):
        fit = tr & (d['y'][:, ti] >= 0)
        y, groups = d['y'][fit, ti], d['animal'][fit]
        q = float(np.mean([np.mean(y[groups == a]) for a in sorted(set(groups))]))
        ca = audit['comparators']['per_target'][str(target)]
        close_tree({'q_fit': q, 'n_train': int(fit.sum())}, ca)
        for name, x in comp_views.items():
            close_tree(adapter_from_train(x[fit]), ca['transforms'][name])
        for family, x in views.items():
            worst = max(worst, check_training_audit(audit['hgbt'][f'{family}|{target}'], x[fit], y, groups))
    return {'labels_checked': len(expected), 'training_prior_max_abs_error': worst,
            'npz_sha256': sha(path), 'train_animals': tg, 'test_animals': eg,
            'all_training_adapters_and_priors_recomputed': True}


def verify_artifacts(run_dir, selection, protocol):
    """Load only our known completed, hash-pinned joblib files; never score training predictions."""
    import joblib
    manifest, ids, genes, pos, groups, labels = a1_inputs(protocol)
    out = {}
    for family in FAMILIES:
        pin = selection['artifacts'][family]
        require(Path(pin['metadata']).name == pin['metadata'] and Path(pin['model_file']).name == pin['model_file'],
                'artifact path must be a basename')
        meta_path, model_path = run_dir/pin['metadata'], run_dir/pin['model_file']
        require(sha(meta_path) == pin['metadata_sha256'] and sha(model_path) == pin['model_sha256'],
                'artifact hash differs from known selection receipt')
        meta = load(meta_path)
        require(meta['model_file'] == model_path.name and meta['model_sha256'] == pin['model_sha256'],
                'artifact metadata model binding differs')
        require(meta['kind'] == 'DEPLOYMENT_REFIT_NOT_OOF' and meta['family'] == family,
                'deployment artifact kind differs')
        require(meta['alpha'] == selection['alphas'][family], 'artifact selected alpha differs')
        gnames = manifest['endogenous_genes'] if family == FAMILIES[0] else []
        require(meta['gene_features'] == gnames and meta['position_features'] == manifest['position_columns'] and
                meta['feature_order'] == gnames+manifest['position_columns'], 'artifact feature order differs')
        x = np.column_stack([np.log1p(genes), pos]) if family == FAMILIES[0] else pos
        require(set(meta['targets']) == {t for _, t in labels}, 'artifact target list differs')
        models = joblib.load(model_path)
        require(set(models) == set(meta['targets']), 'artifact fitted model target list differs')
        for target, target_meta in meta['targets'].items():
            y = np.array([labels[(cell, target)] for cell in ids]); fit = y >= 0
            require(target_meta['train_cell_ids'] == sorted(np.asarray(ids)[fit]) and
                    target_meta['train_groups'] == sorted(set(groups[fit])) and target_meta['n_fit'] == int(fit.sum()),
                    'deployment fit identity differs')
            q = float(np.mean([np.mean(y[fit & (groups == group)]) for group in sorted(set(groups[fit]))]))
            close_tree({'q_fit': q, 'adapter': adapter_from_train(x[fit])}, target_meta)
            require(models[target].n_features_in_ == x.shape[1], 'artifact model feature width differs')
            close_tree(meta['hgbt'], models[target].get_params())
            adapter = target_meta['adapter']
            probe = x[:8].copy()
            z = (probe-np.asarray(adapter['scaler_mean']))/np.asarray(adapter['scaler_scale'])
            p = (1-meta['alpha'])*models[target].predict_proba(z)[:, 1]+meta['alpha']*q
            require(np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all(), 'artifact inference invalid')
        out[family] = {'metadata_sha256': sha(meta_path), 'model_sha256': sha(model_path),
                       'targets': len(models), 'feature_only_inference_finite': True,
                       'fit_identity_adapter_and_prior_verified': True, 'scored_as_oof': False}
    return out


def parse_rows(path, simulation=False):
    """Key every row, reject duplicate/colliding identities and invalid probabilities."""
    cases = {}
    count = 0
    with Path(path).open(newline='') as stream:
        for row in csv.DictReader(stream):
            case = row['case_id'] if simulation else 'A1'
            if simulation:
                scenario, seed = row['scenario'], int(row['seed'])
                require(scenario in NEW_SEEDS and seed in NEW_SEEDS[scenario], 'undeclared synthetic case')
                require(case == f'{scenario}|{seed}', 'case_id disagrees with scenario/seed')
            group = row['animal_id'] if simulation else row['group_id']
            target = row['target'] if simulation else row['endpoint']
            y = int(row['observed_label'] if simulation else row['captured'])
            require(y in (-1, 0, 1), 'invalid captured label')
            assayed = row['assayed'] == 'True' if simulation else True
            if simulation:
                require(row['assayed'] in ('True', 'False'), 'invalid assay flag')
                require((y >= 0) == assayed, 'assay flag and label disagree')
            p = float(row['probability'])
            require(np.isfinite(p) and 0 <= p <= 1, 'invalid probability')
            raw = None if row['raw_probability'] == '' else float(row['raw_probability'])
            require(raw is None or (np.isfinite(raw) and 0 <= raw <= 1), 'invalid raw probability')
            family = row['method'].split('|')[0]
            require((raw is not None) == (family in FAMILIES), 'raw probability availability differs')
            fold = None if simulation else int(row['fold'])
            obs_key = (row['cell_id'], target)
            data = cases.setdefault(case, {})
            observation = data.setdefault(obs_key, {'group': group, 'target': target, 'y': y,
                                                     'assayed': assayed, 'fold': fold, 'methods': {}})
            for attr, value in [('group', group), ('y', y), ('assayed', assayed), ('fold', fold)]:
                require(observation[attr] == value, f'inconsistent observation {case}/{obs_key}/{attr}')
            require(row['method'] not in observation['methods'], 'duplicate forecast key')
            observation['methods'][row['method']] = {'p': p, 'raw': raw}
            count += 1
    require(count > 0, 'empty forecast file')
    for case, data in cases.items():
        method_sets = {tuple(sorted(v['methods'])) for v in data.values()}
        require(len(method_sets) == 1, f'{case}: missing method forecast')
        cell_groups = {}
        for (cell, _), obs in data.items():
            require(cell_groups.setdefault(cell, obs['group']) == obs['group'], 'cell has inconsistent animal')
    return cases, count


def score_case(data):
    names = sorted(next(iter(data.values()))['methods'])
    targets = sorted({o['target'] for o in data.values()})
    animals = sorted({o['group'] for o in data.values()})
    result = {}
    for name in names:
        per_t = {}
        for target in targets:
            per_a = {}
            for animal in animals:
                obs = [o for o in data.values() if o['target'] == target and o['group'] == animal and o['assayed']]
                require(bool(obs), f'animal {animal} lacks assayed target {target}')
                y = np.asarray([o['y'] for o in obs], float)
                p = np.asarray([o['methods'][name]['p'] for o in obs], float)
                pc = np.clip(p, *CLIP)
                per_a[animal] = {'brier': float(np.mean((p-y)**2)),
                                 'log_loss': float(np.mean(-y*np.log(pc)-(1-y)*np.log1p(-pc)))}
            per_t[target] = {metric: float(np.mean([v[metric] for v in per_a.values()]))
                             for metric in ('brier', 'log_loss')}
            per_t[target]['per_animal'] = per_a
        result[name] = {metric+'_target_equal': float(np.mean([v[metric] for v in per_t.values()]))
                        for metric in ('brier', 'log_loss')}
        result[name]['per_target'] = per_t
        result[name]['per_animal_macro'] = {
            a: {metric: float(np.mean([per_t[t]['per_animal'][a][metric] for t in targets]))
                for metric in ('brier', 'log_loss')} for a in animals}
    return result


def select(scores, family):
    names = [f'{family}|alpha={a}' for a in ALPHAS]
    baseline = scores[names[0]]['brier_target_equal']
    order = lambda name: (scores[name]['log_loss_target_equal'], scores[name]['brier_target_equal'],
                           -float(name.split('=')[-1]))
    eligible = [name for name in names if scores[name]['brier_target_equal'] <= baseline]
    best, free = min(eligible, key=order), min(names, key=order)
    return {'selected': best, 'alpha': float(best.split('=')[-1]), 'n_pass_guard': len(eligible),
            'unconstrained_best': free, 'guard_rejected_unconstrained_best': best != free,
            'log_loss_target_equal': scores[best]['log_loss_target_equal'],
            'brier_target_equal': scores[best]['brier_target_equal'],
            'alpha0_log_loss_target_equal': scores[names[0]]['log_loss_target_equal'],
            'alpha0_brier_target_equal': scores[names[0]]['brier_target_equal']}


def gate(scores, candidate, frozen, depth, controls):
    required = [frozen, depth] + list(controls)
    checks = {name: {'lower_log_loss': scores[candidate]['log_loss_target_equal'] < scores[name]['log_loss_target_equal'],
                     'brier_not_increased': scores[candidate]['brier_target_equal'] <= scores[name]['brier_target_equal']}
              for name in required}
    fails = {name: check for name, check in checks.items() if not all(check.values())}
    return {'met': not fails, 'required': required, 'per_comparator': checks, 'failures': fails}


def paired(scores, candidate, controls):
    out = {}
    for control in controls:
        diff = {'log_loss_diff': scores[control]['log_loss_target_equal']-scores[candidate]['log_loss_target_equal'],
                'brier_diff': scores[control]['brier_target_equal']-scores[candidate]['brier_target_equal'],
                'per_animal_macro_diff': {
                    a: {m: scores[control]['per_animal_macro'][a][m]-scores[candidate]['per_animal_macro'][a][m]
                        for m in ('brier', 'log_loss')} for a in scores[candidate]['per_animal_macro']},
                'n_targets': len(scores[control]['per_target'])}
        for metric in ('brier', 'log_loss'):
            diff['targets_lower_'+metric] = sum(scores[control]['per_target'][t][metric] > scores[candidate]['per_target'][t][metric]
                                                for t in scores[control]['per_target'])
        out[control] = diff
    return out


def check_shrinkage(data, audit, simulation=False):
    worst = 0.0
    for (_, target), obs in data.items():
        for name, vals in obs['methods'].items():
            family = name.split('|')[0]
            if family not in FAMILIES:
                continue
            alpha = float(name.split('=')[-1])
            if simulation:
                q = audit['hgbt'][f'{family}|{target}']['q_fit']
            else:
                q = audit['per_endpoint_fold'][f'{target}|{obs["fold"]}'][family]['q_fit']
            require(np.isfinite(q) and 0 <= q <= 1, 'invalid TRAIN prior')
            error = abs(vals['p']-((1-alpha)*vals['raw']+alpha*q))
            worst = max(error, worst)
            require(error <= TOL, 'raw to shrunk forecast differs')
            for other, more in obs['methods'].items():
                if other.split('|')[0] == family:
                    require(abs(more['raw']-vals['raw']) <= TOL, 'raw forecasts differ across alphas')
    return worst


def verify_old_oof(data, path, names):
    expected = {(cell, target, name): obs for (cell, target), obs in data.items() for name in names}
    seen = set()
    worst = 0.0
    with Path(path).open(newline='') as stream:
        for row in csv.DictReader(stream):
            if row['config'] not in names:
                continue
            key = (row['cell_id'], row['endpoint'], row['config'])
            require(key in expected and key not in seen, 'old OOF key unexpected or duplicated')
            seen.add(key)
            obs = expected[key]
            require(obs['group'] == row['group_id'] and obs['fold'] == int(row['fold']) and obs['y'] == int(row['captured']),
                    'old OOF identity, fold or label differs')
            error = abs(obs['methods'][row['config']]['p']-float(row['probability']))
            worst = max(error, worst)
            require(error <= TOL, 'old OOF probability differs')
    require(seen == set(expected), 'old OOF keys missing')
    return {'keys_checked': len(seen), 'max_abs_diff': worst, 'tolerance': TOL}


def sim_name(name):
    return {'candidate_logistic|C=1.0|alpha=0.0': 'candidate_original|C=1.0|alpha=0.0',
            'copy_1nn_location|alpha=0.0': 'copy_1nn_location_raw|alpha=0.0',
            'copy_1nn_matched|alpha=0.0': 'copy_1nn_matched_raw|alpha=0.0'}.get(name, name)


def null_expectation(data):
    """Exact expected proper-score risk for hidden-state null, P(Y|X)=0.5.

    These are prediction-conditional expected risks, not empirical lower bounds.
    Apply the same assay/animal/target weighting as the reported empirical loss.
    """
    out = {}
    for name in next(iter(data.values()))['methods']:
        groups = {}
        for obs in data.values():
            if not obs['assayed']:
                continue
            p = obs['methods'][name]['p']; pc = np.clip(p, *CLIP)
            groups.setdefault((obs['target'], obs['group']), []).append(
                (.25+(p-.5)**2, -.5*np.log(pc)-.5*np.log1p(-pc)))
        means = {}
        for (target, _), values in groups.items():
            means.setdefault(target, []).append(np.mean(values, axis=0))
        br, ll = np.mean([np.mean(v, axis=0) for v in means.values()], axis=0)
        require(br >= .25-TOL and ll >= np.log(2)-TOL, 'null expected proper-score bound violated')
        out[name] = {'expected_brier': float(br), 'expected_log_loss': float(ll)}
    return out


def aggregate_cases(scores, candidate):
    names = sorted(next(iter(scores.values())))
    def aggregate(keys):
        result = {}
        for name in names:
            if name == candidate:
                continue
            diffs = {m: [scores[k][name][m+'_target_equal']-scores[k][candidate][m+'_target_equal'] for k in keys]
                     for m in ('brier', 'log_loss')}
            r = {'n_cases': len(keys)}
            for metric, values in diffs.items():
                r.update({metric+'_diff_mean': float(np.mean(values)), metric+'_diff_min': float(np.min(values)),
                          metric+'_diff_max': float(np.max(values))})
            r['cases_candidate_lower_log_loss'] = sum(v > 0 for v in diffs['log_loss'])
            r['cases_candidate_brier_not_higher'] = sum(v >= 0 for v in diffs['brier'])
            result[name] = r
        return result
    return {'candidate': candidate,
            'per_scenario': {s: aggregate([k for k in scores if k.startswith(s+'|')]) for s in NEW_SEEDS},
            'overall': aggregate(list(scores))}


def descriptive_metrics(scores):
    """Every method's per-case values and per-scenario mean/range; no ranking inference."""
    per_case = {case: {name: {m: sc[m+'_target_equal'] for m in ('brier', 'log_loss')}
                       for name, sc in methods.items()} for case, methods in scores.items()}
    scenarios = {}
    for scenario in NEW_SEEDS:
        keys = sorted(k for k in scores if k.startswith(scenario+'|'))
        scenarios[scenario] = {}
        for name in sorted(next(iter(scores.values()))):
            scenarios[scenario][name] = {}
            for metric in ('brier', 'log_loss'):
                values = [scores[k][name][metric+'_target_equal'] for k in keys]
                scenarios[scenario][name][metric] = {'mean': float(np.mean(values)),
                    'min': float(np.min(values)), 'max': float(np.max(values)),
                    'values': {k: float(scores[k][name][metric+'_target_equal']) for k in keys}}
    return {'per_case': per_case, 'per_scenario_all_methods': scenarios}


def replay_a1(run_dir, protocol):
    done = verify_completion(run_dir)
    cases, n = parse_rows(run_dir/'predictions.csv')
    data = cases['A1']
    require(len(data) == 598*11 and n == 598*11*19, 'A1 forecast row count differs')
    require({o['group'] for o in data.values()} == {'XC95', 'XC97'}, 'A1 animals differ')
    require(len({c for c, _ in data}) == 598 and len({t for _, t in data}) == 11, 'A1 cells/targets differ')
    result = load(run_dir/'result.json'); audit = load(run_dir/'fit_audit.json')
    scores = score_case(data)
    close_tree(scores, result['all_configurations'])
    require(set(scores) == set(protocol['comparators']['a1_keys']) |
            {f'{f}|alpha={a}' for f in FAMILIES for a in ALPHAS}, 'A1 method list differs')
    selections = {f: select(scores, f) for f in FAMILIES}
    close_tree(selections, result['selection'])
    candidate = selections[FAMILIES[0]]['selected']; depth = selections[FAMILIES[1]]['selected']
    require(result['candidate'] == candidate, 'selected candidate differs')
    checks = gate(scores, candidate, FROZEN_LR, depth, protocol['comparators']['strong_copy_averaging'])
    close_tree(checks, result['gate'])
    close_tree(paired(scores, candidate, protocol['comparators']['a1_keys']+[depth]), result['paired_vs_candidate'])
    old = ROOT/protocol['pins']['development_predictions']
    require(sha(old) == protocol['pins']['development_predictions_sha256'], 'old OOF pin differs')
    prior = verify_old_oof(data, old, protocol['comparators']['a1_keys'])
    close_tree(prior, result['baseline_replay'])
    sel = load(run_dir/'selection.json')
    close_tree({'alphas': {f: selections[f]['alpha'] for f in FAMILIES},
                'selected_keys': {f: selections[f]['selected'] for f in FAMILIES}, 'gate_met': checks['met']}, sel)
    shrink = check_shrinkage(data, audit)
    training = verify_a1_training(data, audit, protocol)
    artifacts = verify_artifacts(run_dir, sel, protocol)
    close_tree({'result_sha256': sha(run_dir/'result.json'), 'predictions_sha256': sha(run_dir/'predictions.csv'),
                'run_complete_sha256': sha(run_dir/'RUN_COMPLETE.json'),
                'protocol_snapshot_sha256': sha(run_dir/'protocol_snapshot.json')}, sel)
    return {'phase': 'A1_EXPOSED_DEVELOPMENT_ONLY', 'prediction_rows': n,
            'methods': len(scores), 'scores_selection_gate_replayed': True,
            'baseline_oof_replay': prior, 'raw_shrinkage_max_abs_error': shrink,
            'selected': selections, 'gate': checks, 'scores': scores,
            'training_audit_replay': training, 'deployment_artifact_checks': artifacts,
            'completion_sha256': sha(run_dir/'RUN_COMPLETE.json'), 'refits_performed': 0}


def replay_sim(run_dir, protocol):
    done = verify_completion(run_dir)
    cases, n = parse_rows(run_dir/'predictions.csv', simulation=True)
    expected_cases = {f'{s}|{seed}' for s, seeds in NEW_SEEDS.items() for seed in seeds}
    require(set(cases) == expected_cases, 'synthetic cases missing or unexpected')
    require(n == 12*400*4*15, 'synthetic forecast row count differs')
    result = load(run_dir/'metrics.json')
    require(set(result['cases']) == expected_cases, 'reported synthetic cases differ')
    require(result['settings_selected_on_these_outcomes'] is False, 'test selection disclosure differs')
    record_path = ROOT/protocol['selection_record']['path']
    record = load(record_path)
    started = load(run_dir/'STARTED.json')
    require(sha(record_path) == started['selection_record_sha256'], 'committed selection record pin differs')
    a1_dir = run_dir.parent/'nonlinear_a1_dev_v1'
    require(sha(a1_dir/'selection.json') == record['a1_selection_sha256'] and
            sha(a1_dir/'RUN_COMPLETE.json') == record['a1_run_complete_sha256'], 'A1 selection binding differs')
    selection = load(a1_dir/'selection.json')
    require(selection['alphas'] == record['alphas'], 'frozen selected alpha record differs')
    alphas = record['alphas']
    candidate = f'{FAMILIES[0]}|alpha={float(alphas[FAMILIES[0]])}'
    depth = f'{FAMILIES[1]}|alpha={float(alphas[FAMILIES[1]])}'
    all_scores, all_gates, nulls, training = {}, {}, {}, {}
    worst = 0.0
    names = {sim_name(n) for n in protocol['comparators']['a1_keys']} | {candidate, depth}
    for case in sorted(cases):
        data = cases[case]; scenario, seed = case.split('|')
        require(len(data) == 400*4 and len({c for c, _ in data}) == 400, 'synthetic cells/rows differ')
        require(len({o['group'] for o in data.values()}) == 4 and len({t for _, t in data}) == 4,
                'synthetic animals/targets differ')
        scores = score_case(data)
        require(set(scores) == names, 'synthetic method list differs')
        reported = result['cases'][case]
        for name, sc in scores.items():
            close_tree({m: sc[m] for m in ('brier_target_equal', 'log_loss_target_equal')}, reported['scores'][name])
            close_tree(sc['per_target'], reported['per_target'][name])
        checks = gate(scores, candidate, FROZEN_LR, depth,
                      [sim_name(n) for n in protocol['comparators']['strong_copy_averaging']])
        close_tree(checks, reported['gate'])
        audit = load(run_dir/f'fit_audit_{scenario}_{seed}.json')
        training[case] = verify_sim_training(run_dir, data, audit, case)
        worst = max(worst, check_shrinkage(data, audit, simulation=True))
        all_scores[case], all_gates[case] = scores, checks
        if scenario == 'latent_unobserved':
            nulls[case] = null_expectation(data)
    summary = aggregate_cases(all_scores, candidate)
    close_tree(summary, result['summary'])
    return {'phase': 'FRESH_GENERIC_SYNTHETIC_ONLY', 'prediction_rows': n,
            'cases': len(cases), 'methods': len(names), 'scores_and_gates_replayed': True,
            'summary': summary, 'all_method_metrics': descriptive_metrics(all_scores),
            'raw_shrinkage_max_abs_error': worst,
            'per_case_gate': all_gates, 'hidden_state_expected_risk': nulls,
            'training_data_audit': training, 'frozen_selection_record_sha256': sha(record_path),
            'meaning': 'Three seeds per scenario are bounded sensitivity, not biological validation.',
            'completion_sha256': sha(run_dir/'RUN_COMPLETE.json'), 'refits_performed': 0}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase', choices=('a1', 'sim'))
    p.add_argument('--run-dir', type=Path)
    p.add_argument('--protocol', type=Path, default=ROOT/'docs/paired_projection_nonlinear_protocol.json')
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    rd = args.run_dir or ROOT/'data/interim/paired_projection_nonlinear'/(
        'nonlinear_a1_dev_v1' if args.phase == 'a1' else 'nonlinear_sim_v1')
    protocol = load(args.protocol)
    receipt = replay_a1(rd, protocol) if args.phase == 'a1' else replay_sim(rd, protocol)
    receipt.update({'status': 'PASS', 'replay_source_sha256': sha(__file__),
                    'predictions_sha256': sha(rd/'predictions.csv'), 'protocol_sha256': sha(args.protocol)})
    output = args.output or rd/'ROOT_REPLAY.json'
    with output.open('x') as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True, allow_nan=False); stream.write('\n')
    print(json.dumps({'status': 'PASS', 'receipt': str(output), 'sha256': sha(output),
                      'prediction_rows': receipt['prediction_rows'], 'refits_performed': 0}))


if __name__ == '__main__':
    main()
