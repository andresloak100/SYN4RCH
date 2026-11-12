"""Generic SYNTHETIC animal-held-out sensitivity benchmark; never HY validation.

Only cached, previously exposed TRAIN metadata set coordinate units/scale.
Genes, targets, missingness and their couplings are explicitly invented. The
generator is a spatial mixture plus nonlinear radial fields, latent classes,
Gamma-Poisson counts and Bernoulli assay noise, not the logistic predictor.

Run: python -I analysis/paired_projection_simulation.py run
No downloads, author code, real-data fits, parameter search or test adaptation.
"""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / 'docs/paired_projection_simulation_protocol.json'
OUT = ROOT / 'data/interim/paired_projection_simulation'
SCENARIOS = ('location_only', 'complementarity', 'latent_unobserved', 'batch_shift_missing')
LOGIT = dict(solver='lbfgs', max_iter=2000, tol=1e-6, random_state=20261112)
CLIP = (1e-6, 1 - 1e-6)
LABEL = 'GENERIC_SYNTHETIC_SENSITIVITY_ONLY_NOT_HYPOTHALAMIC_VALIDATION'
PINS = {
    'data/interim/training_donor_selector_development/input/metadata.json':
        'e2397a54fc0655ae4b33708bc65835ea5921c3d45131b6ee6c52a4b08ceae6f8',
    'data/interim/training_donor_selector_development/input/whitelist.json':
        'f8b5aa6f4ea0e2ed36797676fff84b2a73db985f95bdfb8239c63d7084fa77db',
    'data/interim/training_donor_selector_development/regional_profiles/profile_input.json':
        '86d867e6502472b4422e47e507dd35d4ba1ea8042a7ecdc4d7931ce7095a489c',
    'data/interim/paired_projection_logloss_development/logloss_dev_barseq2_a1_v1/frozen_recipe_candidate.json':
        '06e510963264fd3d915410805ffbef07a81e19637ec2b6d9adc460882c8ba02e',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def anchors():
    """Read exactly three cached TRAIN metadata files and the authorized recipe.

    No full Jiao source, regional profile values, arbors, split manifest or old
    prediction/score output is opened, including for hash verification.
    """
    for rel, expected in PINS.items():
        if sha(ROOT / rel) != expected:
            raise ValueError(f'pinned metadata changed: {rel}')
    names = list(PINS)
    rows = json.loads((ROOT / names[0]).read_text())
    whitelist = json.loads((ROOT / names[1]).read_text())
    profile_pin = json.loads((ROOT / names[2]).read_text())
    recipe = json.loads((ROOT / names[3]).read_text())
    if len(rows) != 37 or {r['id'] for r in rows} != set(whitelist):
        raise ValueError('cached metadata are not exactly the exposed TRAIN whitelist')
    if profile_pin['whitelist_sha256'] != PINS[names[1]] or profile_pin['rows'] != 37:
        raise ValueError('TRAIN provenance crosscheck failed')
    xyz = np.asarray([r['soma_um'] for r in rows], dtype=float)
    if xyz.shape != (37, 3) or not np.isfinite(xyz).all():
        raise ValueError('invalid TRAIN coordinate metadata')
    receipt = {
        'classification': LABEL,
        'scope': 'TRAIN HY coordinate units/scale only; remaining simulation is generic',
        'source_rows': 37,
        'source_animal_proxies': len({r['animal'] for r in rows}),
        'coordinate_axes': ['AP_um', 'DV_um', 'ML_um'],
        'coordinate_frame': 'CCFv3 mirrored-left, from existing TRAIN preparation contract',
        'coordinate_centre_um': xyz.mean(axis=0).tolist(),
        'coordinate_scale_um': xyz.std(axis=0).tolist(),
        'coordinate_min_um': xyz.min(axis=0).tolist(),
        'coordinate_max_um': xyz.max(axis=0).tolist(),
        'peptide_tag_categories': len({r['peptide'] for r in rows}),
        'region_vocabulary_size': len(profile_pin['regions']),
        'region_profile_values_read': False,
        'coupling_anchor': 'NONE: no same-cell multigene-to-projection coupling estimated',
        'gene_count_assay_or_missingness_anchor': 'NONE: all numeric values are sensitivity hypotheses',
        'population_warning': 'The 37 selected TRAIN neurons are not a representative HY population sample.',
        'raw_source_pin_recorded_not_opened': profile_pin['source_sha256'],
        'read_files_sha256': dict(PINS),
    }
    return receipt, recipe


def _nonlinear_probability(scenario, z, molecular_state, hidden_state):
    """Independent generator: nonlinear fields and discrete mixture responses."""
    centres = np.array([[-.9, -.4, .3], [.8, -.6, -.3], [.2, .8, .5], [-.3, .1, -.8]])
    radial = np.exp(-.65 * ((z[:, None, :] - centres[None, :, :]) ** 2).sum(axis=2))
    wave = (np.sin(1.4 * z[:, 0]) * np.tanh(z[:, 1]) + 1) / 2
    if scenario == 'latent_unobserved':
        # Each target has exactly two high and two low latent states. Since H
        # is uniform and independent of everything measured, P(Y|X,G)=0.5.
        hidden = np.array([[1, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 1]])
        return .08 + .84 * hidden[hidden_state]
    if scenario == 'location_only':
        # Genes and all latent classes are excluded from the response law.
        return .08 + .58 * radial + .16 * wave[:, None]
    states = np.array([[1, 0, 1, 0], [0, 1, 1, 0], [1, 0, 0, 1], [0, 1, 0, 1]])
    molecular = states[molecular_state]
    return .06 + .32 * radial + .42 * molecular + .12 * radial * molecular + .05 * wave[:, None]


def generate(scenario, seed, anchor, n_train=6, n_test=4, cells_per_animal=100, n_genes=16):
    """Make new simulated animals. No real cell or real gene/target IDs reused."""
    if scenario not in SCENARIOS:
        raise ValueError('unknown scenario')
    if min(n_train, n_test, cells_per_animal) < 1 or n_genes < 4:
        raise ValueError('invalid synthetic design')
    centre = np.asarray(anchor['coordinate_centre_um'])
    scale = np.asarray(anchor['coordinate_scale_um'])
    if centre.shape != (3,) or scale.shape != (3,) or (scale <= 0).any():
        raise ValueError('invalid anchor scale')
    seeds = np.random.SeedSequence(seed).spawn(n_train + n_test + 1)
    parameter_rng = np.random.default_rng(seeds[0])
    # Shared biological gene parameters are fixed independently of model fits.
    gene_base = parameter_rng.uniform(np.log(1.2), np.log(4.0), n_genes)
    spatial_centres = np.array([[-.7, -.4, .2], [.7, -.2, -.4], [.1, .7, .5]])
    pieces = []
    for a, animal_seed in enumerate(seeds[1:]):
        rng = np.random.default_rng(animal_seed)
        test = a >= n_train
        animal_id = f'SYN_{scenario}_{"TEST" if test else "TRAIN"}_{a:02d}'
        n = cells_per_animal
        position_class = rng.choice(3, size=n, p=rng.dirichlet(np.full(3, 12.)))
        z = spatial_centres[position_class] + rng.normal(0, [.55, .5, .45], (n, 3))
        z += rng.normal(0, .12, 3)
        molecular = rng.integers(0, 4, n)
        hidden = rng.integers(0, 4, n)
        p_true = _nonlinear_probability(scenario, z, molecular, hidden)
        y_true = (rng.random((n, 4)) < p_true).astype(np.int8)
        library = rng.lognormal(0, .3, n)
        marker = (np.arange(n_genes)[None, :] % 4 == molecular[:, None])
        log_mean = gene_base[None, :] + 2.3 * marker
        if scenario in ('complementarity', 'batch_shift_missing'):
            log_mean = log_mean + .18 * np.tanh(z[:, np.arange(n_genes) % 3])
        observed_z = z.copy()
        if scenario == 'batch_shift_missing':
            batch = rng.normal(0, .3, n_genes)
            if test:
                batch += np.where(np.arange(n_genes) % 2 == 0, .9, -.7)
            log_mean += batch
            observed_z += rng.normal(0, .40 if test else .12, (n, 3))
            observed_z += rng.normal(.35 if test else 0, .25 if test else .08, 3)
        mean = np.exp(log_mean) * library[:, None]
        # Gamma-Poisson overdispersion, followed by independent detection dropout.
        counts = rng.poisson(rng.gamma(shape=2., scale=mean / 2.)).astype(float)
        counts[rng.random(counts.shape) < .06] = 0.
        xyz = centre + observed_z * scale
        assayed = np.ones((n, 4), dtype=bool)
        y = y_true.copy()
        p_assay = p_true.copy()
        if scenario == 'batch_shift_missing':
            gene_missing = .15 if test else .04
            counts[rng.random(counts.shape) < gene_missing] = np.nan
            # Two declared features are missing for an entire animal, rather
            # than silently changing gene identity or imputing from test data.
            counts[:, rng.choice(n_genes, 2, replace=False)] = np.nan
            xyz[rng.random(xyz.shape) < (.18 if test else .04)] = np.nan
            sensitivity, false_positive = (.72, .08) if test else (.92, .02)
            y = np.where(y_true == 1, rng.random(y.shape) < sensitivity,
                         rng.random(y.shape) < false_positive).astype(np.int8)
            p_assay = false_positive + (sensitivity - false_positive) * p_true
            assayed = rng.random(y.shape) >= (.1 if test else .04)
            y[~assayed] = -1
        tech = np.column_stack([np.log(library), np.isnan(counts).mean(axis=1),
                                np.isnan(xyz).mean(axis=1)])
        pieces.append(dict(
            animal=np.full(n, animal_id), split=np.full(n, 'test' if test else 'train'),
            cell_id=np.asarray([f'{animal_id}_CELL_{i:04d}' for i in range(n)]),
            genes=counts, xyz=xyz, tech=tech, y=y, assayed=assayed,
            generating_probability=p_assay, true_projection=y_true,
            latent_molecular=molecular, latent_hidden=hidden,
        ))
    d = {k: np.concatenate([piece[k] for piece in pieces]) for k in pieces[0]}
    d['gene_names'] = np.asarray([f'SYN_G{i + 1:02d}' for i in range(n_genes)])
    d['target_names'] = np.asarray([f'SYN_T{i + 1}' for i in range(4)])
    d['tech_names'] = np.asarray(['SYN_log_library_factor', 'SYN_gene_missing_fraction',
                                 'SYN_coordinate_missing_fraction'])
    d['classification'] = np.asarray(LABEL)
    return d


def group_weights(groups):
    _, inv, count = np.unique(groups, return_inverse=True, return_counts=True)
    w = 1. / count[inv]
    return w / w.mean()


def fit_view(train, evaluation):
    """Fixed missingness adapter: train median then train StandardScaler.

    This adapter is an explicitly declared simulation extension to a frozen
    recipe that requires complete features. It is identical across methods.
    No missingness indicators are added to the gene/position feature panels.
    """
    from sklearn.preprocessing import StandardScaler
    median = np.array([np.median(col[np.isfinite(col)]) if np.isfinite(col).any() else 0.
                       for col in train.T])
    a = np.where(np.isfinite(train), train, median)
    b = np.where(np.isfinite(evaluation), evaluation, median)
    scaler = StandardScaler().fit(a)
    return scaler.transform(a), scaler.transform(b), {
        'median': median.tolist(), 'scaler_mean': scaler.mean_.tolist(),
        'scaler_scale': scaler.scale_.tolist(),
    }


def method_configs(recipe):
    """All development-frozen families, old candidate and two raw-copy anchors."""
    result = []
    for family, key in recipe['controls_frozen_with_it'].items():
        cfg = {'family': family, 'name': key}
        for field in key.split('|')[1:]:
            k, value = field.split('=')
            cfg[k] = int(value) if k == 'k' else float(value)
        result.append(cfg)
    result.extend([
        {'family': 'candidate_logistic', 'name': 'candidate_original|C=1.0|alpha=0.0', 'C': 1., 'alpha': 0.},
        {'family': 'copy_1nn_location', 'name': 'copy_1nn_location_raw|alpha=0.0', 'alpha': 0.},
        {'family': 'copy_1nn_matched', 'name': 'copy_1nn_matched_raw|alpha=0.0', 'alpha': 0.},
    ])
    if len({c['name'] for c in result}) != len(result):
        raise ValueError('duplicated frozen methods')
    return result


def fit_predict(d, recipe):
    """Fit once on synthetic TRAIN animals, predict TEST features without labels.

    Test labels, generative probabilities and all true/latent states are not
    read here. Training and test IDs must be strictly disjoint.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
    split = d['split']
    train, evaluation = split == 'train', split == 'test'
    train_animals, test_animals = set(d['animal'][train]), set(d['animal'][evaluation])
    if not train_animals or not test_animals or train_animals & test_animals:
        raise ValueError('synthetic animal split is empty or overlapping')
    configs = method_configs(recipe)
    genes = np.log1p(d['genes'])
    views = {'candidate_logistic': np.column_stack([genes, d['xyz']]),
             'location_logistic': d['xyz'], 'molecular': genes, 'technical': d['tech']}
    spaces = {'location_knn': 'location_logistic', 'molecular_knn': 'molecular',
              'matched_knn': 'candidate_logistic', 'copy_1nn_location': 'location_logistic',
              'copy_1nn_matched': 'candidate_logistic'}
    predictions = {c['name']: np.empty((evaluation.sum(), 4)) for c in configs}
    audit = {'train_animals': sorted(train_animals), 'test_animals': sorted(test_animals),
             'test_labels_read_by_fitter': False, 'per_target': {}}
    for ti, target in enumerate(d['target_names']):
        # Index TRAIN first: no test label values are examined even for eligibility.
        eligible = np.zeros(len(split), dtype=bool)
        eligible[train] = d['y'][train, ti] >= 0
        tr = train & eligible
        y_train = d['y'][tr, ti]
        if len(np.unique(y_train)) != 2 or tr.sum() < 31:
            raise ValueError('training target has one class or fewer than 31 cells')
        w = group_weights(d['animal'][tr])
        q = float(np.average(y_train, weights=w))
        z, transformation = {}, {}
        for name, x in views.items():
            a, b, fitted = fit_view(x[tr], x[evaluation])
            z[name] = (a, b)
            transformation[name] = fitted
        target_audit = dict(q_fit=q, n_train=int(tr.sum()), transforms=transformation, models={}, donors={})
        cache = {}
        for cfg in configs:
            family, key = cfg['family'], cfg['name']
            alpha = cfg.get('alpha', 0.)
            if family.endswith('logistic') or family == 'technical':
                a, b = z[family]
                model = LogisticRegression(**LOGIT, C=cfg['C']).fit(a, y_train, sample_weight=w)
                raw = model.predict_proba(b)[:, 1]
                target_audit['models'][key] = {'coef': model.coef_[0].tolist(),
                                               'intercept': float(model.intercept_[0])}
            elif family in spaces:
                a, b = z[spaces[family]]
                cache_key = (family, cfg.get('k', 1))
                if cache_key not in cache:
                    if family.startswith('copy_'):
                        donor = NearestNeighbors(n_neighbors=1).fit(a).kneighbors(b, return_distance=False)[:, 0]
                        raw = y_train[donor].astype(float)
                        donor_groups = d['animal'][tr][donor]
                        if set(donor_groups) & test_animals:
                            raise ValueError('copy donor from test animal')
                        target_audit['donors'][family] = sorted(set(donor_groups))
                    else:
                        raw = KNeighborsClassifier(n_neighbors=cfg['k']).fit(a, y_train).predict_proba(b)[:, 1]
                    cache[cache_key] = raw
                raw = cache[cache_key]
            elif family == 'train_frequency':
                raw = np.full(evaluation.sum(), q)
            elif family == 'constant_0.5':
                raw = np.full(evaluation.sum(), .5)
            else:
                raise ValueError(f'unknown frozen method {family}')
            probability = (1 - alpha) * raw + alpha * q
            if not np.isfinite(probability).all() or ((probability < 0) | (probability > 1)).any():
                raise ValueError('invalid model probability')
            predictions[key][:, ti] = probability
        audit['per_target'][str(target)] = target_audit
    return predictions, audit


def score(d, predictions):
    """Animal-equal within target, then target-equal. Unknown labels excluded."""
    ev = d['split'] == 'test'
    y, animals = d['y'][ev], d['animal'][ev]
    result = {}
    for method, probability in predictions.items():
        per_target = {}
        for ti, target in enumerate(d['target_names']):
            per_animal = {}
            for animal in sorted(set(animals)):
                ok = (animals == animal) & (y[:, ti] >= 0)
                if not ok.any():
                    raise ValueError('test animal lacks an assayed target')
                truth, p = y[ok, ti], probability[ok, ti]
                clipped = np.clip(p, *CLIP)
                per_animal[str(animal)] = {
                    'n_assayed': int(ok.sum()), 'n_positive': int(truth.sum()),
                    'brier': float(np.mean((p - truth) ** 2)),
                    'log_loss': float(np.mean(-truth * np.log(clipped) - (1 - truth) * np.log1p(-clipped))),
                }
            per_target[str(target)] = {
                metric: float(np.mean([v[metric] for v in per_animal.values()]))
                for metric in ('brier', 'log_loss')}
            per_target[str(target)]['per_animal'] = per_animal
        result[method] = {metric + '_target_equal': float(np.mean([v[metric] for v in per_target.values()]))
                          for metric in ('brier', 'log_loss')}
        result[method]['per_target'] = per_target
    return result


def run(output=OUT / 'run_v1', protocol_path=PROTOCOL):
    import scipy
    import sklearn
    protocol = json.loads(Path(protocol_path).read_text())
    if protocol['status'] != 'FROZEN_SYNTHETIC_SENSITIVITY_ONLY':
        raise ValueError('protocol is not frozen for synthetic sensitivity work')
    if protocol['design'] != dict(train_animals=6, test_animals=4, cells_per_animal=100, genes=16, targets=4):
        raise ValueError('bounded design changed; do not tune against synthetic tests')
    if tuple(protocol['scenarios']) != SCENARIOS:
        raise ValueError('frozen scenario list changed')
    if protocol['pins'] != PINS:
        raise ValueError('metadata pins differ from frozen protocol')
    anchor, recipe = anchors()
    configs = method_configs(recipe)
    if [c['name'] for c in configs] != protocol['methods']:
        raise ValueError('frozen comparison settings changed')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / 'STARTED.json', {
        'classification': LABEL, 'started_utc': datetime.now(timezone.utc).isoformat(),
        'protocol_sha256': sha(protocol_path), 'code_sha256': sha(__file__),
        'seeds': protocol['seeds'], 'parameters_selected_on_synthetic_tests': False,
        'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                        'scipy': scipy.__version__, 'scikit_learn': sklearn.__version__},
    })
    write_json(output / 'anchor_receipt.json', anchor)
    write_json(output / 'frozen_recipe_input.json', recipe)
    all_scores, summary_rows = {}, []
    with (output / 'predictions.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['classification', 'scenario', 'cell_id', 'animal_id', 'target',
                         'assayed', 'observed_label', 'method', 'probability'])
        for scenario in SCENARIOS:
            d = generate(scenario, protocol['seeds'][scenario], anchor)
            np.savez_compressed(output / f'SYNTHETIC_{scenario}.npz', **d)
            p, fitted = fit_predict(d, recipe)
            write_json(output / f'fit_audit_{scenario}.json', fitted)
            scores = score(d, p)
            diagnostics = {
                'gene_missing_fraction': float(np.isnan(d['genes']).mean()),
                'coordinate_missing_fraction': float(np.isnan(d['xyz']).mean()),
                'unassayed_label_fraction': float((d['y'] < 0).mean()),
                'train_animals': fitted['train_animals'], 'test_animals': fitted['test_animals'],
                'generative_hypothesis': protocol['scenarios'][scenario],
            }
            if scenario == 'latent_unobserved':
                diagnostics['exact_expected_observable_bayes_risk'] = {'brier': .25, 'log_loss': float(np.log(2))}
                diagnostics['exact_expected_latent_oracle_risk'] = {
                    'brier': .08 * .92, 'log_loss': float(-.08 * np.log(.08) - .92 * np.log(.92))}
                diagnostics['risk_warning'] = 'Analytic expectations, not empirical finite-sample lower bounds.'
            all_scores[scenario] = {'methods': scores, 'diagnostics': diagnostics}
            for method, metrics in scores.items():
                summary_rows.append(dict(classification=LABEL, scenario=scenario, method=method,
                                         brier=metrics['brier_target_equal'], log_loss=metrics['log_loss_target_equal']))
            ev = np.flatnonzero(d['split'] == 'test')
            for i, row in enumerate(ev):
                for ti, target in enumerate(d['target_names']):
                    for method, probability in p.items():
                        writer.writerow([LABEL, scenario, d['cell_id'][row], d['animal'][row], target,
                                         bool(d['assayed'][row, ti]), int(d['y'][row, ti]), method,
                                         float(probability[i, ti])])
    write_json(output / 'metrics.json', {
        'classification': LABEL, 'scenarios': all_scores,
        'scope': 'Software and sensitivity findings only; no biological validation or method selection.',
        'target_names': ['SYN_T1', 'SYN_T2', 'SYN_T3', 'SYN_T4'],
        'parameters_selected_on_synthetic_tests': False,
    })
    with (output / 'summary.csv').open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    files = {p.name: {'sha256': sha(p), 'bytes': p.stat().st_size}
             for p in sorted(output.iterdir()) if p.is_file()}
    done = {'classification': LABEL, 'status': 'SYNTHETIC_RUN_COMPLETE',
            'completed_utc': datetime.now(timezone.utc).isoformat(),
            'code_sha256': sha(__file__), 'protocol_sha256': sha(protocol_path),
            'files': files, 'test_parameter_searches': 0, 'scenarios': len(SCENARIOS),
            'train_animals_per_scenario': 6, 'test_animals_per_scenario': 4,
            'biological_validation_claimed': False, 'real_training_fits_performed': 0}
    write_json(output / 'RUN_COMPLETE.json', done)
    return done


if __name__ == '__main__':
    if sys.argv[1:] != ['run']:
        raise SystemExit(__doc__)
    result = run()
    print(json.dumps({k: result[k] for k in ('classification', 'status', 'scenarios', 'test_parameter_searches')}))
