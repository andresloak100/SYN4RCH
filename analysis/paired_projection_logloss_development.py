"""Bounded log-loss development round: stronger L2 + train-only shrinkage, with equally tuned controls.

Protocol: docs/paired_projection_logloss_development_protocol.json (DRAFT until ROOT freezes it).

* EXPOSED DEVELOPMENT ONLY on the BARseq2 A1 598-cell inputs (XC95, XC97); not a test, no nested claim.
* Imports the frozen runner (analysis/paired_projection_real_benchmark.py) and the reviewed reference
  module for loading, validation and metrics; edits neither.
* Shrinkage p' = (1 - alpha) p + alpha q_fit, q_fit = group-balanced captured rate in the fit animals only.
* 12-point opportunity for every non-copy family:
    logistic families (candidate, location, technical): C {0.01, 0.1, 1} x alpha {0, 0.1, 0.25, 0.5}
    averaging families (location, molecular, matched KNN): k {5, 15, 31} x alpha
    copy families (1NN location, 1NN matched): alpha (alpha=0 is the raw 1NN anchor)
* One global setting per family by leave-one-animal-out OOF target-equal log loss, subject to that family's
  own Brier non-degradation guard (tolerance 0). Every configuration is reported, plus the unrestricted
  log-loss best when the guard rejects it.
* Deployment: selected settings are refit on ALL admitted training animals (both mice) and frozen as a
  recipe plus an A1-specific artifact (coefficients, scaler, q_fit). Future predictions take no labels.
  A different gene panel must refit the recipe on its own training animals; A1 coefficients do not transfer.

Usage:
  /opt/anaconda3/bin/python -I analysis/paired_projection_logloss_development.py synthetic
  /opt/anaconda3/bin/python -I analysis/paired_projection_logloss_development.py run   (refused until frozen)
"""
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis import paired_projection_real_benchmark as B   # noqa: E402  frozen runner (import only)
from analysis import paired_projection_reference as REF      # noqa: E402  reviewed helpers (import only)

PROTOCOL = ROOT / 'docs/paired_projection_logloss_development_protocol.json'
OUT = ROOT / 'data/interim/paired_projection_logloss_development'
DEV_RUN_ID = 'logloss_dev_barseq2_a1_v1'
FROZEN_DEV_STATUS = 'FROZEN_DEVELOPMENT_AUTHORIZED'
C_GRID = (0.01, 0.1, 1.0)
K_GRID = (5, 15, 31)
ALPHA_GRID = (0.0, 0.1, 0.25, 0.5)
BRIER_GUARD_TOL = 0.0
LOGISTIC_FAMILIES = ('candidate_logistic', 'location_logistic', 'technical')
AVERAGING_FAMILIES = ('location_knn', 'molecular_knn', 'matched_knn')
COPY_FAMILIES = ('copy_1nn_location', 'copy_1nn_matched')
REFERENCE_FAMILIES = ('train_frequency', 'constant_0.5')
BASELINE = {'candidate_logistic': dict(C=1.0, alpha=0.0), 'location_logistic': dict(C=1.0, alpha=0.0),
            'technical': dict(C=1.0, alpha=0.0), 'location_knn': dict(k=15, alpha=0.0),
            'molecular_knn': dict(k=15, alpha=0.0), 'matched_knn': dict(k=15, alpha=0.0),
            'copy_1nn_location': dict(alpha=0.0), 'copy_1nn_matched': dict(alpha=0.0),
            'train_frequency': {}, 'constant_0.5': {}}
Refusal = REF.Refusal
sha256_file = REF.sha256_file


def cfg_key(family, C=None, k=None, alpha=None):
    parts = [family]
    if C is not None:
        parts.append(f'C={C}')
    if k is not None:
        parts.append(f'k={k}')
    if alpha is not None:
        parts.append(f'alpha={alpha}')
    return '|'.join(parts)


def parse_key(key):
    fam, *rest = key.split('|')
    out = {'family': fam}
    for r in rest:
        n, v = r.split('=')
        out[n] = int(v) if n == 'k' else float(v)
    return out


def baseline_key(family):
    return cfg_key(family, **BASELINE[family])


def all_config_keys(families):
    keys = []
    for f in families:
        if f in LOGISTIC_FAMILIES:
            keys += [cfg_key(f, C=C, alpha=a) for C in C_GRID for a in ALPHA_GRID]
        elif f in AVERAGING_FAMILIES:
            keys += [cfg_key(f, k=k, alpha=a) for k in K_GRID for a in ALPHA_GRID]
        elif f in COPY_FAMILIES:
            keys += [cfg_key(f, alpha=a) for a in ALPHA_GRID]
        else:
            keys.append(cfg_key(f))
    return keys


def shrink(p, q, alpha):
    return (1.0 - alpha) * p + alpha * q


# ---------------------------------------------------------------- pins (real development run only)

def check_dev_pins(manifest_path, repo=ROOT, protocol_path=PROTOCOL, code_path=None):
    """Committed frozen protocol + committed code + exact input/admission/runner/reference hashes."""
    proto = json.loads(Path(protocol_path).read_text())
    if proto.get('status') != FROZEN_DEV_STATUS:
        raise Refusal(f"development protocol status {proto.get('status')!r} is not {FROZEN_DEV_STATUS}")
    code_path = Path(code_path or __file__)
    for p, what in ((protocol_path, 'development protocol'), (code_path, 'development code')):
        if not B.protocol_committed(p, repo):
            raise Refusal(f'{what} is not committed in git or has uncommitted changes')
    pins = proto['data_status']['pins']
    checks = {'input manifest': (manifest_path, pins['input_manifest_sha256']),
              'frozen runner': (B.__file__, pins['frozen_runner_sha256']),
              'reference module': (REF.__file__, pins['reference_sha256'])}
    for what, (p, want) in checks.items():
        if sha256_file(p) != want:
            raise Refusal(f'{what} SHA-256 differs from the development protocol pin')
    m = json.loads(Path(manifest_path).read_text())
    if sha256_file(Path(manifest_path).parent / m['admission']['record_path']) != pins['admission_sha256']:
        raise Refusal('admission record SHA-256 differs from the development protocol pin')
    return proto


def load(manifest_path):
    """Validate and assemble with the frozen runner's own checks (real pins for real manifests)."""
    manifest_path = Path(manifest_path)
    m = json.loads(manifest_path.read_text())
    B.validate_manifest(m, manifest_path.parent)
    return m, B.assemble(m, manifest_path.parent)


# ---------------------------------------------------------------- fitting one fold

def families_for(d):
    fams = ['candidate_logistic', 'location_logistic'] + (['technical'] if d['tech'].shape[1] else [])
    return fams + list(AVERAGING_FAMILIES) + list(COPY_FAMILIES) + list(REFERENCE_FAMILIES)


def views(d):
    return {'candidate_logistic': np.hstack([d['genes'], d['xyz']]), 'location_logistic': d['xyz'],
            'technical': d['tech'], 'molecular': d['genes']}


def fit_predict_grid(d, tr, ev, y_tr):
    """Every configuration's probabilities for the evaluation rows, fitted on training rows only."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
    from sklearn.preprocessing import StandardScaler
    w = REF.group_balanced_weights(d['group'][tr])
    q = float(np.average(y_tr, weights=w))                       # q_fit: training animals only
    V = views(d)
    sc = {name: StandardScaler().fit(X[tr]) for name, X in V.items() if X.shape[1]}
    Z = {name: (sc[name].transform(V[name][tr]), sc[name].transform(V[name][ev])) for name in sc}
    raw = {}
    for fam in LOGISTIC_FAMILIES:
        if fam not in Z:
            continue
        A, Bm = Z[fam]
        for C in C_GRID:
            raw[cfg_key(fam, C=C)] = LogisticRegression(**dict(B.LOGIT, C=C)).fit(A, y_tr, sample_weight=w) \
                .predict_proba(Bm)[:, 1]
    knn_space = {'location_knn': 'location_logistic', 'molecular_knn': 'molecular', 'matched_knn': 'candidate_logistic'}
    for fam, space in knn_space.items():
        A, Bm = Z[space]
        for k in K_GRID:
            raw[cfg_key(fam, k=k)] = KNeighborsClassifier(n_neighbors=k).fit(A, y_tr).predict_proba(Bm)[:, 1]
    donors = {}
    for fam, space in {'copy_1nn_location': 'location_logistic', 'copy_1nn_matched': 'candidate_logistic'}.items():
        A, Bm = Z[space]
        idx = NearestNeighbors(n_neighbors=1).fit(A).kneighbors(Bm, return_distance=False)[:, 0]
        raw[fam] = y_tr[idx].astype(float)
        donors[fam] = d['group'][tr][idx]
    out = {}
    for base, p in raw.items():
        for a in ALPHA_GRID:
            out[f'{base}|alpha={a}'] = shrink(p, q, a)
    n = int(ev.sum())
    out[cfg_key('train_frequency')] = np.full(n, q)
    out[cfg_key('constant_0.5')] = np.full(n, 0.5)
    return out, donors, q, sc['candidate_logistic']


# ---------------------------------------------------------------- evaluation of every configuration

def evaluate_grid(m, d):
    fams = families_for(d)
    keys = all_config_keys(fams)
    folds = B.folds_for(m, d['group'])
    preds, fold_log = {}, []
    for ep in d['endpoints']:
        y_all = d['y'][ep]
        elig = y_all >= 0
        P = {k: np.full(len(y_all), np.nan) for k in keys}
        fold_of = np.full(len(y_all), -1)
        for fi, eval_groups in enumerate(folds):
            train_groups = set(d['group']) - eval_groups
            REF.assert_disjoint(train_groups, eval_groups)
            tr = elig & np.isin(d['group'], list(train_groups))
            ev = elig & np.isin(d['group'], list(eval_groups))
            if not ev.any():
                continue
            y_tr = y_all[tr]
            if tr.sum() < max(K_GRID) or len(set(y_tr)) < 2:
                raise Refusal(f'{ep} fold {fi}: fewer than {max(K_GRID)} training cells or one class; FAIL')
            out, donors, q, s_full = fit_predict_grid(d, tr, ev, y_tr)
            if not np.allclose(s_full.mean_, views(d)['candidate_logistic'][tr].mean(axis=0)):
                raise Refusal('scaler not fitted on training animals only')
            for fam, dg in donors.items():
                if np.isin(dg, list(eval_groups)).any():
                    raise Refusal(f'{fam}: donor from an evaluation animal')
            if set(out) != set(keys):
                raise Refusal(f'{ep} fold {fi}: configuration set mismatch {sorted(set(out) ^ set(keys))[:4]}')
            for k in keys:
                P[k][ev] = REF.checked_probability(out[k], int(ev.sum()), f'{ep} fold {fi} {k}')
            fold_of[ev] = fi
            fold_log.append({'endpoint': ep, 'fold': fi, 'eval_groups': sorted(eval_groups),
                             'n_train': int(tr.sum()), 'n_eval': int(ev.sum()), 'q_fit': q})
        for k in keys:
            if not np.isfinite(P[k][elig]).all():
                raise Refusal(f'{ep}: {k} missing a probability; FAIL')
        preds[ep] = {'p': P, 'fold_of': fold_of}
    return {'families': fams, 'keys': keys, 'folds': [sorted(f) for f in folds], 'fold_log': fold_log,
            'preds': preds}


def score_configs(d, grid):
    """Per configuration: per target animal-equal scores, target-equal macro, per-animal macro."""
    eps = list(grid['preds'])
    scores = {}
    for k in grid['keys']:
        per_t = {}
        for ep in eps:
            y = d['y'][ep]; e = y >= 0
            pg = B.per_group_scores(grid['preds'][ep]['p'][k][e], y[e], d['group'][e])
            per_t[ep] = {'brier': float(np.mean([v['brier'] for v in pg.values()])),
                         'log_loss': float(np.mean([v['log_loss'] for v in pg.values()])), 'per_animal': pg}
        animals = sorted({a for ep in eps for a in per_t[ep]['per_animal']})
        per_a = {a: {s: float(np.mean([per_t[ep]['per_animal'][a][s] for ep in eps if a in per_t[ep]['per_animal']]))
                     for s in ('brier', 'log_loss')} for a in animals}
        scores[k] = {'log_loss_target_equal': float(np.mean([per_t[ep]['log_loss'] for ep in eps])),
                     'brier_target_equal': float(np.mean([per_t[ep]['brier'] for ep in eps])),
                     'per_animal_macro': per_a, 'per_target': per_t}
    return scores


def _order(scores, key):
    """Min log loss; tie-break lower Brier, larger alpha, smaller C, smaller k (deterministic)."""
    p = parse_key(key)
    return (scores[key]['log_loss_target_equal'], scores[key]['brier_target_equal'], -p.get('alpha', 0.0),
            p.get('C', 0.0), p.get('k', 0))


def select(scores, families):
    """One global setting per family: min target-equal log loss subject to its own Brier guard."""
    chosen = {}
    for f in families:
        cands = [k for k in scores if parse_key(k)['family'] == f]
        base = baseline_key(f)
        limit = scores[base]['brier_target_equal'] + BRIER_GUARD_TOL
        ok = [k for k in cands if scores[k]['brier_target_equal'] <= limit]
        best = min(ok, key=lambda k: _order(scores, k))
        unrestricted = min(cands, key=lambda k: _order(scores, k))
        chosen[f] = {'selected': best, 'baseline': base, 'n_configs': len(cands), 'n_pass_guard': len(ok),
                     'log_loss_target_equal': scores[best]['log_loss_target_equal'],
                     'brier_target_equal': scores[best]['brier_target_equal'],
                     'baseline_log_loss_target_equal': scores[base]['log_loss_target_equal'],
                     'baseline_brier_target_equal': scores[base]['brier_target_equal'],
                     'unrestricted_log_loss_best': unrestricted,
                     'guard_rejected_unrestricted_best': unrestricted != best,
                     'unrestricted_log_loss_target_equal': scores[unrestricted]['log_loss_target_equal'],
                     'unrestricted_brier_target_equal': scores[unrestricted]['brier_target_equal']}
    return chosen


def compare(scores, chosen, against='candidate_logistic'):
    """Selected candidate vs each selected (tuned) control and each untuned baseline; control minus candidate."""
    cand = chosen[against]['selected']
    out = {}
    for f, c in chosen.items():
        if f == against:
            continue
        for label, k in (('tuned', c['selected']), ('baseline', c['baseline'])):
            pa = {a: {s: scores[k]['per_animal_macro'][a][s] - scores[cand]['per_animal_macro'][a][s]
                      for s in ('brier', 'log_loss')} for a in scores[cand]['per_animal_macro']}
            out[f'{f}|{label}'] = {
                'control_config': k,
                'log_loss_target_equal_diff': scores[k]['log_loss_target_equal'] - scores[cand]['log_loss_target_equal'],
                'brier_target_equal_diff': scores[k]['brier_target_equal'] - scores[cand]['brier_target_equal'],
                'per_animal_macro_diff': pa,
                'animals_improved_log_loss': int(sum(v['log_loss'] > 0 for v in pa.values())),
                'animals_improved_brier': int(sum(v['brier'] > 0 for v in pa.values())),
                'targets_improved_log_loss': int(sum(scores[k]['per_target'][t]['log_loss'] >
                                                     scores[cand]['per_target'][t]['log_loss'] for t in scores[k]['per_target'])),
                'targets_improved_brier': int(sum(scores[k]['per_target'][t]['brier'] >
                                                  scores[cand]['per_target'][t]['brier'] for t in scores[k]['per_target'])),
                'n_targets': len(scores[k]['per_target'])}
    return {'candidate_config': cand, 'meaning': 'control minus candidate; positive = candidate better; '
            'leave-one-animal-out OOF development estimates, optimistic (selected on the same predictions)',
            'pairs': out}


# ---------------------------------------------------------------- deployment: refit on all training animals

def deploy_fit(d, setting, gene_names, position_names, endpoints=None):
    """Refit the selected candidate setting on ALL admitted training animals; returns a frozen artifact."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    p = parse_key(setting)
    if p['family'] != 'candidate_logistic':
        raise Refusal('the deployment artifact is for the candidate family only')
    X = views(d)['candidate_logistic']
    art = {'schema': 'synarch_paired_projection_deployment_artifact_v1', 'setting': setting,
           'C': p['C'], 'alpha': p['alpha'], 'features': list(gene_names) + list(position_names),
           'gene_features': list(gene_names), 'position_features': list(position_names),
           'gene_transform': 'log1p', 'training_animals': sorted(set(d['group'])), 'targets': {}}
    for ep in endpoints or d['endpoints']:
        y = d['y'][ep]; e = y >= 0
        w = REF.group_balanced_weights(d['group'][e])
        sc = StandardScaler().fit(X[e])
        lr = LogisticRegression(**dict(B.LOGIT, C=p['C'])).fit(sc.transform(X[e]), y[e], sample_weight=w)
        art['targets'][ep] = {'scaler_mean': sc.mean_.tolist(), 'scaler_scale': sc.scale_.tolist(),
                              'coef': lr.coef_[0].tolist(), 'intercept': float(lr.intercept_[0]),
                              'q_fit': float(np.average(y[e], weights=w)), 'n_fit': int(e.sum()),
                              'n_fit_per_animal': {g: int((d['group'][e] == g).sum()) for g in sorted(set(d['group'][e]))}}
    return art


def deploy_predict(art, gene_names, position_names, gene_counts, positions):
    """Predict for new cells from features only (no labels are accepted). Panel must match exactly."""
    from scipy.special import expit
    if list(gene_names) != art['gene_features'] or list(position_names) != art['position_features']:
        raise Refusal('feature panel differs from the frozen artifact; refit the frozen recipe on the new '
                      'panel training animals instead (no A1 coefficient transfer)')
    arrays = {}
    for name, x, width in (('gene_counts', gene_counts, len(art['gene_features'])),
                           ('positions', positions, len(art['position_features']))):
        if isinstance(x, (str, bytes)):
            raise Refusal(f'{name} must be a numeric 2-D array')
        a = np.asarray(x)
        if a.dtype == object or not (np.issubdtype(a.dtype, np.integer) or np.issubdtype(a.dtype, np.floating)) \
                or a.dtype == bool:
            raise Refusal(f'{name} must be numeric (integer or float), got {a.dtype}')
        if a.ndim != 2 or a.shape[0] < 1 or a.shape[1] != width:
            raise Refusal(f'{name} must have shape (n_cells, {width}), got {a.shape}')
        a = a.astype(float)
        if not np.isfinite(a).all():
            raise Refusal(f'{name} contains non-finite values')
        arrays[name] = a
    if arrays['gene_counts'].shape[0] != arrays['positions'].shape[0]:
        raise Refusal('gene_counts and positions have different numbers of rows')
    if (arrays['gene_counts'] < 0).any():
        raise Refusal('negative gene counts refused before log1p')
    X = np.hstack([np.log1p(arrays['gene_counts']), arrays['positions']])
    n = X.shape[0]
    out = {}
    for ep, t in art['targets'].items():
        z = ((X - np.array(t['scaler_mean'])) / np.array(t['scaler_scale'])) @ np.array(t['coef']) + t['intercept']
        p = shrink(expit(z), t['q_fit'], art['alpha'])
        out[ep] = REF.checked_probability(p, n, f'deploy {ep}')
    return out


def recipe(chosen):
    p = parse_key(chosen['candidate_logistic']['selected'])
    return {'schema': 'synarch_paired_projection_frozen_recipe_v1',
            'status': 'DEVELOPMENT_SELECTED_FREEZE_AFTER_ROOT_REVIEW',
            'family': 'candidate_logistic', 'C': p['C'], 'alpha': p['alpha'],
            'model': 'sklearn LogisticRegression L2 (default), lbfgs, max_iter 2000, tol 1e-6, random_state 20261112',
            'inputs': 'log1p measured endogenous genes + measured position (declared axes); StandardScaler fitted on '
                      'the training animals only',
            'fit_weights': 'group-balanced 1/n_group, mean-normalised; no class_weight',
            'shrinkage': "p' = (1 - alpha) p + alpha q_fit, q_fit = group-balanced captured rate in the training animals only",
            'controls_frozen_with_it': {f: v['selected'] for f, v in chosen.items()},
            'selection_basis': 'leave-one-animal-out OOF on the 2 exposed A1 animals; development-optimistic',
            'deployment': 'refit on all admitted training animals; predictions take features only, never labels',
            'not_transferable': 'A1 23-gene coefficients are not carried over; a new (e.g. hypothalamic) panel refits '
                                'this recipe only on its own admitted training animals before any held-out test',
            'future_rule': 'no new-animal parameter selection; any calibration from fit-side data only'}


def development_round(m, d):
    grid = evaluate_grid(m, d)
    scores = score_configs(d, grid)
    chosen = select(scores, grid['families'])
    art = deploy_fit(d, chosen['candidate_logistic']['selected'], m['endogenous_genes'], m['position_columns'])
    return grid, scores, chosen, compare(scores, chosen), recipe(chosen), art


# ---------------------------------------------------------------- durable execution

def execute(manifest_path, out_root=OUT, run_id=DEV_RUN_ID, repo=ROOT, protocol_path=PROTOCOL, code_path=None):
    manifest_path = Path(manifest_path)
    m = json.loads(manifest_path.read_text())
    if m.get('synthetic') is not True:
        # real development run: fixed identity, canonical output root, the actually executing code
        if run_id != DEV_RUN_ID:
            raise Refusal(f'real development run_id must be {DEV_RUN_ID}')
        if Path(out_root).resolve() != OUT.resolve():
            raise Refusal(f'real development output root must be {OUT}')
        if code_path is not None and Path(code_path).resolve() != Path(__file__).resolve():
            raise Refusal('code_path override refused for real runs: the executing file is what is pinned')
        check_dev_pins(manifest_path, repo, protocol_path, Path(__file__))   # before any data file is parsed
    run_dir = Path(out_root) / run_id
    try:
        run_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        raise Refusal(f'development run_id {run_id} already exists; duplicate runs are refused')
    prov = {'run_id': run_id, 'synthetic': m.get('synthetic'), 'exposure': 'EXPOSED_DEVELOPMENT_ONLY',
            'manifest_sha256': sha256_file(manifest_path), 'code_sha256': sha256_file(code_path or __file__),
            'protocol_sha256': sha256_file(protocol_path), 'frozen_runner_sha256': sha256_file(B.__file__),
            'reference_sha256': sha256_file(REF.__file__), 'C_grid': C_GRID, 'k_grid': K_GRID,
            'alpha_grid': ALPHA_GRID, 'brier_guard_tol': BRIER_GUARD_TOL,
            'started_utc': datetime.now(timezone.utc).isoformat()}
    with (run_dir / 'STARTED.json').open('x') as f:
        json.dump(prov, f, indent=2, sort_keys=True)
    m, d = load(manifest_path)
    grid, scores, chosen, comp, rec, art = development_round(m, d)
    summary = {'exposure': 'EXPOSED DEVELOPMENT ONLY; leave-one-animal-out OOF selection is optimistic; not a test',
               'families': grid['families'], 'folds': grid['folds'], 'fold_log': grid['fold_log'],
               'selected': chosen, 'comparison': comp,
               'all_configurations': {k: {'log_loss_target_equal': v['log_loss_target_equal'],
                                          'brier_target_equal': v['brier_target_equal'],
                                          'per_animal_macro': v['per_animal_macro'],
                                          'per_target': {t: {'brier': x['brier'], 'log_loss': x['log_loss'],
                                                             'per_animal': x['per_animal']}
                                                         for t, x in v['per_target'].items()}}
                                      for k, v in scores.items()}}
    files = {'development_result.json': summary, 'frozen_recipe_candidate.json': rec,
             'deployment_artifact_a1.json': art}
    for name, obj in files.items():
        with (run_dir / name).open('x') as f:
            json.dump(obj, f, indent=2, sort_keys=True, allow_nan=False)
    with (run_dir / 'predictions.csv').open('x', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['cell_id', 'group_id', 'endpoint', 'fold', 'config', 'probability', 'captured'])
        rows = []
        for ep, g in grid['preds'].items():
            y = d['y'][ep]
            for i in np.flatnonzero(y >= 0):
                for k in grid['keys']:
                    rows.append((d['cell_id'][i], d['group'][i], ep, int(g['fold_of'][i]), k,
                                 repr(float(g['p'][k][i])), int(y[i])))
        rows.sort()
        w.writerows(rows)
    done = dict(prov, **{f'{n}_sha256': sha256_file(run_dir / n) for n in
                         list(files) + ['predictions.csv', 'STARTED.json']},
                finished_utc=datetime.now(timezone.utc).isoformat())
    with (run_dir / 'RUN_COMPLETE.json').open('x') as f:
        json.dump(done, f, indent=2, sort_keys=True)
    return run_dir, summary, (m, d, grid, scores, chosen, art)


# ---------------------------------------------------------------- synthetic regression

def synthetic_dataset(folder, sizes=(160, 240), seed=B.SEED, dims=1):
    """Two simulated animals, 1 measured axis, planted gene + position signal (mirrors the A1 shape)."""
    m, t = B.make_synthetic(sizes=sizes, n_genes=8, dims=dims, with_tech=True, seed=seed, unassayed_t2=False)
    return B.write_dataset(m, t, folder)


def synthetic(out_root=OUT / 'synthetic'):
    out_root = Path(out_root)
    mp = synthetic_dataset(out_root / 'data')
    run_dir, summary, (m, d, grid, scores, chosen, art) = execute(mp, out_root=out_root / 'runs', run_id='synthetic_dev_v2')
    ep = d['endpoints'][0]
    y = d['y'][ep]; e = y >= 0
    frozen = np.full(len(y), np.nan)
    for ev_groups in B.folds_for(m, d['group']):
        tr = e & ~np.isin(d['group'], list(ev_groups)); ev = e & np.isin(d['group'], list(ev_groups))
        out, _, _ = B.fit_predict(d, tr, ev, y[tr])
        frozen[ev] = out['candidate_logistic']
    same = bool(np.array_equal(frozen[e], grid['preds'][ep]['p'][cfg_key('candidate_logistic', C=1.0, alpha=0.0)][e]))
    checks = {
        'candidate_C1_alpha0_equals_frozen_runner': same,
        'every_configuration_reported': len(summary['all_configurations']) == len(all_config_keys(grid['families'])),
        'config_count_expected_82': len(grid['keys']) == 3 * 12 + 3 * 12 + 2 * 4 + 2,
        'twelve_point_opportunity': all(chosen[f]['n_configs'] == 12 for f in LOGISTIC_FAMILIES + AVERAGING_FAMILIES),
        'selection_respects_own_guard': all(v['brier_target_equal'] <= v['baseline_brier_target_equal'] + BRIER_GUARD_TOL
                                            for v in chosen.values()),
        'copy_1nn_shrinkage_lowers_log_loss': scores[cfg_key('copy_1nn_location', alpha=0.25)]['log_loss_target_equal'] <
                                              scores[cfg_key('copy_1nn_location', alpha=0.0)]['log_loss_target_equal'],
        'deployment_refit_uses_all_training_animals': all(t['n_fit'] == int((d['y'][ep2] >= 0).sum())
                                                          for ep2, t in art['targets'].items()),
        'two_folds_leave_one_animal_out': len(grid['folds']) == 2,
    }
    receipt = {'schema': 'synarch_logloss_development_synthetic_receipt_v2',
               'meaning': 'SYNTHETIC software regression only; no exposed or new data read',
               'software_status': 'PASS' if all(checks.values()) else 'FAIL', 'checks': checks,
               'selected': chosen, 'run_complete_sha256': sha256_file(run_dir / 'RUN_COMPLETE.json'),
               'code_sha256': sha256_file(__file__), 'protocol_sha256': sha256_file(PROTOCOL),
               'finished_utc': datetime.now(timezone.utc).isoformat()}
    p = out_root / 'synthetic_receipt.json'
    p.write_text(json.dumps(receipt, indent=2, sort_keys=True))
    print('software', receipt['software_status'], checks)
    for f, v in chosen.items():
        print(f'{f:20s} {v["selected"]:36s} LL {v["baseline_log_loss_target_equal"]:.4f} -> {v["log_loss_target_equal"]:.4f}'
              f'  guard_rejected_best={v["guard_rejected_unrestricted_best"]}')
    print('wrote', p, sha256_file(p))
    return receipt


if __name__ == '__main__':
    a = sys.argv[1:]
    if a == ['synthetic']:
        synthetic()
    elif a == ['run']:
        mp = ROOT / json.loads(PROTOCOL.read_text())['data_status']['pins']['input_manifest']
        run_dir, summary, _ = execute(mp)
        print('development run complete', run_dir, sha256_file(run_dir / 'RUN_COMPLETE.json'))
    else:
        sys.exit(__doc__)
