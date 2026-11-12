"""Fixed shallow nonlinear (HGBT) gene-and-position model vs the frozen model and copying controls.

Protocol: docs/paired_projection_nonlinear_protocol.json (one frozen status covers both stages; DRAFT refuses).

* One fixed HistGradientBoostingClassifier structure per target; only alpha {0, 0.1, 0.25} shrinkage toward
  the training-animal q_fit. Families: hgbt_gene_position (log1p genes + position) and hgbt_position.
* Same training-only adapter as the frozen simulation (median fill then StandardScaler), group-balanced
  weights, no class_weight, no structure search.
* Stage 1, A1 (EXPOSED DEVELOPMENT): exact frozen LOAO folds; the 13 comparators are reproduced with the frozen
  development code and must equal the saved OOF predictions (unique exact keys, within 1e-12); the
  all-training deployment artifact is never scored as OOF. One global alpha per family.
* Deployment artifact per selected family: refit on all A1 development rows; own joblib file, hash-verified
  before loading; features-only prediction API.
* Stage 2, fresh simulation: only after docs/paired_projection_nonlinear_selection.json (written after ROOT's
  replay of stage 1) is committed and matches the completed A1 run and the executing source; new seeds;
  TRAIN-only fitting; TEST labels/masks/oracles never reach the fitter.

Usage:
  /opt/anaconda3/bin/python -I analysis/paired_projection_nonlinear.py a1     (refused until frozen)
  /opt/anaconda3/bin/python -I analysis/paired_projection_nonlinear.py sim    (refused until frozen + committed selection)
"""
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis import paired_projection_logloss_development as D   # noqa: E402  frozen dev code (import only)
from analysis import paired_projection_real_benchmark as B        # noqa: E402  frozen runner (import only)
from analysis import paired_projection_reference as REF           # noqa: E402  reviewed helpers (import only)
from analysis import paired_projection_simulation as SIM          # noqa: E402  frozen simulation (import only)

CANONICAL_SOURCE = ROOT / 'analysis/paired_projection_nonlinear.py'
PROTOCOL = ROOT / 'docs/paired_projection_nonlinear_protocol.json'
SELECTION_RECORD = ROOT / 'docs/paired_projection_nonlinear_selection.json'
OUT = ROOT / 'data/interim/paired_projection_nonlinear'
A1_RUN_ID = 'nonlinear_a1_dev_v1'
SIM_RUN_ID = 'nonlinear_sim_v1'
FROZEN_STATUS = 'FROZEN_NONLINEAR_STAGED_AUTHORIZED'
HGBT = dict(loss='log_loss', learning_rate=0.05, max_iter=100, max_leaf_nodes=4, max_depth=2,
            min_samples_leaf=20, l2_regularization=1.0, early_stopping=False, categorical_features=None,
            class_weight=None, random_state=20261112)
ALPHAS = (0.0, 0.1, 0.25)
FAMILIES = ('hgbt_gene_position', 'hgbt_position')
CANDIDATE = 'hgbt_gene_position'
FROZEN_LR = 'candidate_logistic|C=0.1|alpha=0.1'
A1_COMPARATORS = ['candidate_logistic|C=0.1|alpha=0.1', 'constant_0.5', 'copy_1nn_location|alpha=0.5',
                  'copy_1nn_matched|alpha=0.5', 'location_knn|k=31|alpha=0.5', 'location_logistic|C=0.1|alpha=0.0',
                  'matched_knn|k=31|alpha=0.1', 'molecular_knn|k=31|alpha=0.1', 'technical|C=1.0|alpha=0.0',
                  'train_frequency', 'candidate_logistic|C=1.0|alpha=0.0', 'copy_1nn_location|alpha=0.0',
                  'copy_1nn_matched|alpha=0.0']
A1_TO_SIM = dict(zip(A1_COMPARATORS, A1_COMPARATORS[:10] + ['candidate_original|C=1.0|alpha=0.0',
                                                            'copy_1nn_location_raw|alpha=0.0',
                                                            'copy_1nn_matched_raw|alpha=0.0']))
STRONG_COPY_AVG = ['location_knn|k=31|alpha=0.5', 'molecular_knn|k=31|alpha=0.1', 'matched_knn|k=31|alpha=0.1',
                   'copy_1nn_location|alpha=0.5', 'copy_1nn_matched|alpha=0.5',
                   'copy_1nn_location|alpha=0.0', 'copy_1nn_matched|alpha=0.0']
NEW_SEEDS = {'location_only': (876211, 876212, 876213), 'complementarity': (876221, 876222, 876223),
             'latent_unobserved': (876231, 876232, 876233), 'batch_shift_missing': (876241, 876242, 876243)}
REPLAY_TOL = 1e-12
Refusal = REF.Refusal
sha = REF.sha256_file


def key(family, alpha):
    return f'{family}|alpha={alpha}'


def all_hgbt_keys():
    return [key(f, a) for f in FAMILIES for a in ALPHAS]


def write_json(path, obj):
    with Path(path).open('x') as f:
        json.dump(obj, f, indent=2, sort_keys=True, allow_nan=False)
        f.write('\n')


# ---------------------------------------------------------------- the shared fitter (no evaluation labels)

def _fit_model(x_train, y_train, groups_train):
    from sklearn.ensemble import HistGradientBoostingClassifier
    y_train = np.asarray(y_train)
    if len(np.unique(y_train)) != 2 or not set(np.unique(y_train)) <= {0, 1}:
        raise Refusal('training labels must contain both classes 0 and 1')
    a, _, adapter = SIM.fit_view(np.asarray(x_train, float), np.asarray(x_train, float)[:1])
    w = REF.group_balanced_weights(np.asarray(groups_train))
    q = float(np.average(y_train, weights=w))
    model = HistGradientBoostingClassifier(**HGBT).fit(a, y_train, sample_weight=w)
    return model, adapter, q


def _apply_adapter(adapter, x):
    med = np.asarray(adapter['median']); mu = np.asarray(adapter['scaler_mean']); sd = np.asarray(adapter['scaler_scale'])
    return (np.where(np.isfinite(x), x, med) - mu) / sd


def hgbt_fit_predict(x_train, y_train, groups_train, x_eval):
    """Fit one target on TRAINING rows only and predict evaluation features. Returns raw, shrunk and audit.

    Evaluation labels, masks or oracle values are not parameters of this function."""
    model, adapter, q = _fit_model(x_train, y_train, groups_train)
    b = _apply_adapter(adapter, np.asarray(x_eval, float))
    raw = REF.checked_probability(model.predict_proba(b)[:, 1], len(b), 'hgbt raw')
    shrunk = {a_: REF.checked_probability((1 - a_) * raw + a_ * q, len(b), f'hgbt alpha {a_}') for a_ in ALPHAS}
    audit = {'q_fit': q, 'n_train': int(len(y_train)), 'n_train_positive': int(np.sum(y_train)),
             'adapter': adapter, 'n_iter': int(model.n_iter_), 'train_groups': sorted(set(map(str, groups_train)))}
    return raw, shrunk, audit


# ---------------------------------------------------------------- scoring, selection, comparison, gate

def score(pred, y, groups, keys):
    """pred[endpoint][key] -> animal-equal per target, target-equal macro, per-animal macro."""
    out = {}
    eps = list(pred)
    for k in keys:
        per_t = {}
        for ep in eps:
            e = y[ep] >= 0
            pg = B.per_group_scores(pred[ep][k][e], y[ep][e], groups[e])
            per_t[ep] = {'brier': float(np.mean([v['brier'] for v in pg.values()])),
                         'log_loss': float(np.mean([v['log_loss'] for v in pg.values()])), 'per_animal': pg}
        animals = sorted({a for ep in eps for a in per_t[ep]['per_animal']})
        out[k] = {'log_loss_target_equal': float(np.mean([per_t[ep]['log_loss'] for ep in eps])),
                  'brier_target_equal': float(np.mean([per_t[ep]['brier'] for ep in eps])),
                  'per_animal_macro': {a: {s: float(np.mean([per_t[ep]['per_animal'][a][s] for ep in eps
                                                             if a in per_t[ep]['per_animal']]))
                                           for s in ('brier', 'log_loss')} for a in animals},
                  'per_target': per_t}
    return out


def select_alpha(scores, family):
    """Min target-equal log loss, subject to Brier <= own alpha=0 Brier; ties lower Brier then larger alpha."""
    cands = [key(family, a) for a in ALPHAS]
    limit = scores[key(family, 0.0)]['brier_target_equal']
    order = lambda k: (scores[k]['log_loss_target_equal'], scores[k]['brier_target_equal'], -float(k.split('=')[-1]))
    ok = [k for k in cands if scores[k]['brier_target_equal'] <= limit]
    best, free = min(ok, key=order), min(cands, key=order)
    return {'selected': best, 'alpha': float(best.split('=')[-1]), 'n_pass_guard': len(ok),
            'unconstrained_best': free, 'guard_rejected_unconstrained_best': free != best,
            'log_loss_target_equal': scores[best]['log_loss_target_equal'],
            'brier_target_equal': scores[best]['brier_target_equal'],
            'alpha0_log_loss_target_equal': scores[key(family, 0.0)]['log_loss_target_equal'],
            'alpha0_brier_target_equal': scores[key(family, 0.0)]['brier_target_equal']}


def paired(scores, cand, controls):
    """Control minus candidate; positive = candidate better."""
    out = {}
    for c in controls:
        pa = {a: {s: scores[c]['per_animal_macro'][a][s] - scores[cand]['per_animal_macro'][a][s]
                  for s in ('brier', 'log_loss')} for a in scores[cand]['per_animal_macro']}
        out[c] = {'log_loss_diff': scores[c]['log_loss_target_equal'] - scores[cand]['log_loss_target_equal'],
                  'brier_diff': scores[c]['brier_target_equal'] - scores[cand]['brier_target_equal'],
                  'per_animal_macro_diff': pa,
                  'targets_lower_log_loss': int(sum(scores[c]['per_target'][t]['log_loss'] > scores[cand]['per_target'][t]['log_loss']
                                                    for t in scores[c]['per_target'])),
                  'targets_lower_brier': int(sum(scores[c]['per_target'][t]['brier'] > scores[cand]['per_target'][t]['brier']
                                                 for t in scores[c]['per_target'])),
                  'n_targets': len(scores[c]['per_target'])}
    return out


def gate(scores, cand, frozen_lr, depth_selected, strong):
    """Lower log loss AND non-increased Brier vs each required comparator; failures listed separately."""
    req = [frozen_lr, depth_selected] + list(strong)
    rows = {c: {'lower_log_loss': scores[cand]['log_loss_target_equal'] < scores[c]['log_loss_target_equal'],
                'brier_not_increased': scores[cand]['brier_target_equal'] <= scores[c]['brier_target_equal']}
            for c in req}
    failures = {c: r for c, r in rows.items() if not (r['lower_log_loss'] and r['brier_not_increased'])}
    return {'met': not failures, 'required': req, 'per_comparator': rows, 'failures': failures,
            'decision': 'nonlinear candidate meets the development gate' if not failures else
                        'gate not met: frozen operational model retained; nonlinear candidate reported and still tested'}


# ---------------------------------------------------------------- pins shared by both stages

def source_hashes():
    return {'nonlinear': sha(__file__), 'development': sha(D.__file__), 'frozen_runner': sha(B.__file__),
            'reference': sha(REF.__file__), 'simulation': sha(SIM.__file__), 'simulation_protocol': sha(SIM.PROTOCOL)}


def check_recipes(pins):
    """Published recipe (hash) and archived candidate recipe (hash + settings) pinned separately; must agree."""
    pub_path = ROOT / pins['published_recipe']
    if sha(pub_path) != pins['published_recipe_sha256']:
        raise Refusal('published recipe differs from the pin')
    arch_path = SIM.ROOT / [p for p in SIM.PINS if p.endswith('frozen_recipe_candidate.json')][0]
    if sha(arch_path) != pins['archived_candidate_recipe_sha256']:
        raise Refusal('archived candidate recipe differs from the pin')
    pub, arch = json.loads(pub_path.read_text()), json.loads(arch_path.read_text())
    if pub['controls_frozen_with_it'] != arch['controls_frozen_with_it'] or (pub['C'], pub['alpha']) != (arch['C'], arch['alpha']):
        raise Refusal('published and archived recipes disagree on frozen settings')
    if sorted(pub['controls_frozen_with_it'].values()) != sorted(A1_COMPARATORS[:10]):
        raise Refusal('comparator list differs from the frozen recipe controls')
    return pub, arch


def check_common(repo=ROOT, protocol_path=PROTOCOL, source_path=None):
    """Both stages: frozen committed protocol and source, canonical executing source, all pins and constants."""
    proto = json.loads(Path(protocol_path).read_text())
    if proto.get('status') != FROZEN_STATUS:
        raise Refusal(f"nonlinear protocol status {proto.get('status')!r} is not {FROZEN_STATUS}")
    src = Path(source_path or __file__).resolve()
    if src != CANONICAL_SOURCE.resolve():
        raise Refusal('executing source is not the canonical analysis/paired_projection_nonlinear.py')
    for p, what in ((protocol_path, 'nonlinear protocol'), (src, 'nonlinear source')):
        if not B.protocol_committed(p, repo):
            raise Refusal(f'{what} is not committed in git or has uncommitted changes')
    pins = proto['pins']
    for name, p in (('frozen_runner_sha256', B.__file__), ('reference_sha256', REF.__file__),
                    ('development_code_sha256', D.__file__), ('simulation_code_sha256', SIM.__file__),
                    ('simulation_protocol_sha256', SIM.PROTOCOL)):
        if sha(p) != pins[name]:
            raise Refusal(f'{name} differs from the nonlinear protocol pin')
    check_recipes(pins)
    check_protocol_constants(proto)
    return proto


def check_protocol_constants(proto):
    if proto['candidate']['fixed_parameters'] != HGBT or tuple(proto['candidate']['alpha_grid']) != ALPHAS:
        raise Refusal('protocol fixed parameters or alphas differ from the code')
    if proto['a1_phase']['run_id'] != A1_RUN_ID or proto['simulation_phase']['run_id'] != SIM_RUN_ID:
        raise Refusal('protocol run identities differ from the code')
    if {s: tuple(v) for s, v in proto['simulation_phase']['new_seeds'].items()} != NEW_SEEDS:
        raise Refusal('protocol seeds differ from the code')
    if proto['comparators']['a1_keys'] != A1_COMPARATORS:
        raise Refusal('protocol comparators differ from the code')
    return True


def _begin(run_dir, prov, protocol_path):
    try:
        run_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        raise Refusal(f'{run_dir.name} already exists; duplicate runs are refused')
    shutil.copyfile(protocol_path, run_dir / 'protocol_snapshot.json')
    write_json(run_dir / 'STARTED.json', dict(prov, protocol_snapshot_sha256=sha(run_dir / 'protocol_snapshot.json')))


def provenance(run_id, extra=None):
    import scipy
    import sklearn
    p = {'run_id': run_id, 'started_utc': datetime.now(timezone.utc).isoformat(),
         'executing_source': str(Path(__file__).resolve()), 'source_sha256': source_hashes(),
         'hgbt': HGBT, 'alphas': ALPHAS,
         'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                         'scipy': scipy.__version__, 'scikit_learn': sklearn.__version__}}
    p.update(extra or {})
    return p


def finish(run_dir, extra):
    files = {p.name: sha(p) for p in sorted(run_dir.iterdir()) if p.is_file()}
    write_json(run_dir / 'RUN_COMPLETE.json', dict(extra, files_sha256=files,
                                                    finished_utc=datetime.now(timezone.utc).isoformat()))
    return sha(run_dir / 'RUN_COMPLETE.json')


# ---------------------------------------------------------------- stage 1: A1

def check_a1_pins(proto):
    pins = proto['pins']
    for name, p in (('a1_input_manifest_sha256', ROOT / pins['a1_input_manifest']),
                    ('development_predictions_sha256', ROOT / pins['development_predictions'])):
        if sha(p) != pins[name]:
            raise Refusal(f'{name} differs from the nonlinear protocol pin')
    m = json.loads((ROOT / pins['a1_input_manifest']).read_text())
    if sha((ROOT / pins['a1_input_manifest']).parent / m['admission']['record_path']) != pins['a1_admission_sha256']:
        raise Refusal('A1 admission record differs from the pin')


def a1_round(m, d, saved_predictions=None):
    """HGBT grid + reproduced comparators on the frozen LOAO folds."""
    folds = B.folds_for(m, d['group'])
    feats = a1_features(d)
    keys = A1_COMPARATORS + all_hgbt_keys()
    pred, raw, fold_of, audit = {}, {}, {}, {'folds': [sorted(f) for f in folds], 'per_endpoint_fold': {}}
    for ep in d['endpoints']:
        y = d['y'][ep]; e = y >= 0
        pred[ep] = {k: np.full(len(y), np.nan) for k in keys}
        raw[ep] = {f: np.full(len(y), np.nan) for f in FAMILIES}
        fold_of[ep] = np.full(len(y), -1)
        for fi, evg in enumerate(folds):
            tr = e & ~np.isin(d['group'], list(evg)); ev = e & np.isin(d['group'], list(evg))
            REF.assert_disjoint(set(d['group'][tr]), set(evg))
            if not ev.any():
                continue
            grid, _, _, _ = D.fit_predict_grid(d, tr, ev, y[tr])          # frozen comparators, same fold
            for k in A1_COMPARATORS:
                pred[ep][k][ev] = grid[k]
            fa = {}
            for f, X in feats.items():
                r, s, au = hgbt_fit_predict(X[tr], y[tr], d['group'][tr], X[ev])
                raw[ep][f][ev] = r
                for a_ in ALPHAS:
                    pred[ep][key(f, a_)][ev] = s[a_]
                fa[f] = au
            fold_of[ep][ev] = fi
            audit['per_endpoint_fold'][f'{ep}|{fi}'] = fa
        for k in keys:
            if not np.isfinite(pred[ep][k][e]).all():
                raise Refusal(f'{ep}: {k} missing a probability')
    replay = verify_saved(pred, d, fold_of, saved_predictions) if saved_predictions else None
    scores = score(pred, d['y'], d['group'], keys)
    sel = {f: select_alpha(scores, f) for f in FAMILIES}
    cand = sel[CANDIDATE]['selected']
    controls = A1_COMPARATORS + [sel['hgbt_position']['selected']]
    result = {'selection': sel, 'candidate': cand, 'paired_vs_candidate': paired(scores, cand, controls),
              'gate': gate(scores, cand, FROZEN_LR, sel['hgbt_position']['selected'], STRONG_COPY_AVG),
              'baseline_replay': replay}
    return pred, raw, fold_of, audit, scores, result


def a1_features(d):
    return {'hgbt_gene_position': D.views(d)['candidate_logistic'], 'hgbt_position': d['xyz']}


SAVED_COLUMNS = ['cell_id', 'group_id', 'endpoint', 'fold', 'config', 'probability', 'captured']


def verify_saved(pred, d, fold_of, path):
    """Unique exact keys (cell, group, fold, endpoint, config); finite [0, 1]; labels equal; within 1e-12."""
    idx = {c: i for i, c in enumerate(d['cell_id'])}
    want = {(c, str(d['group'][idx[c]]), int(fold_of[ep][idx[c]]), ep, k)
            for ep in pred for k in A1_COMPARATORS for c in d['cell_id'][d['y'][ep] >= 0]}
    seen, worst = set(), 0.0
    with Path(path).open(newline='') as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != SAVED_COLUMNS:
            raise Refusal(f'saved predictions columns {reader.fieldnames} != {SAVED_COLUMNS}')
        for r in reader:
            if r['config'] not in A1_COMPARATORS:
                continue
            t = (r['cell_id'], r['group_id'], int(r['fold']), r['endpoint'], r['config'])
            if t in seen:
                raise Refusal(f'duplicate saved prediction key {t}')
            if t not in want:
                raise Refusal(f'saved prediction for an unexpected key {t}')
            seen.add(t)
            p = float(r['probability'])
            if not np.isfinite(p) or p < 0 or p > 1:
                raise Refusal(f'saved probability outside [0, 1] or non-finite at {t}')
            i = idx[r['cell_id']]
            if int(r['captured']) != d['y'][r['endpoint']][i]:
                raise Refusal(f'saved label differs at {t}')
            worst = max(worst, abs(p - pred[r['endpoint']][r['config']][i]))
    if seen != want or worst > REPLAY_TOL:
        raise Refusal(f'baseline replay failed: {len(seen)}/{len(want)} exact keys, '
                      f'{len(want - seen)} missing, max diff {worst}')
    return {'keys_checked': len(seen), 'max_abs_diff': worst, 'tolerance': REPLAY_TOL}


# ---------------------------------------------------------------- deployment artifact (selected families)

def deploy_fit(d, family, alpha, gene_names, position_names):
    """Refit on ALL A1 development rows (deployment, not out-of-fold) for every target."""
    X = a1_features(d)[family]
    feats = (list(gene_names) if family == 'hgbt_gene_position' else []) + list(position_names)
    models, meta = {}, {}
    for ep in d['endpoints']:
        e = d['y'][ep] >= 0
        model, adapter, q = _fit_model(X[e], d['y'][ep][e], d['group'][e])
        models[ep] = model
        meta[ep] = {'adapter': adapter, 'q_fit': q, 'n_fit': int(e.sum()),
                    'train_cell_ids': sorted(map(str, d['cell_id'][e])),
                    'train_groups': sorted(set(map(str, d['group'][e])))}
    return {'schema': 'synarch_nonlinear_deployment_artifact_v1', 'kind': 'DEPLOYMENT_REFIT_NOT_OOF',
            'family': family, 'alpha': alpha, 'hgbt': HGBT,
            'gene_features': list(gene_names) if family == 'hgbt_gene_position' else [],
            'position_features': list(position_names), 'feature_order': feats, 'gene_transform': 'log1p',
            'targets': meta}, models


def save_artifact(folder, meta, models):
    import joblib
    folder = Path(folder)
    model_path = folder / f"artifact_{meta['family']}.joblib"
    with model_path.open('xb') as f:
        joblib.dump(models, f)
    meta = dict(meta, model_file=model_path.name, model_sha256=sha(model_path))
    write_json(folder / f"artifact_{meta['family']}.json", meta)
    return folder / f"artifact_{meta['family']}.json"


def load_artifact(meta_path, expected_meta_sha256=None):
    """Load our own artifact only after verifying the metadata and model file hashes (no third-party pickles)."""
    import joblib
    meta_path = Path(meta_path)
    if expected_meta_sha256 is not None and sha(meta_path) != expected_meta_sha256:
        raise Refusal('artifact metadata hash differs from the expected value')
    meta = json.loads(meta_path.read_text())
    if meta.get('schema') != 'synarch_nonlinear_deployment_artifact_v1':
        raise Refusal('not a SYN4RCH nonlinear artifact')
    model_path = meta_path.parent / meta['model_file']
    if sha(model_path) != meta['model_sha256']:
        raise Refusal('artifact model file hash differs; refusing to load')
    return meta, joblib.load(model_path)


def _numeric_2d(name, x, width):
    if isinstance(x, (str, bytes)):
        raise Refusal(f'{name} must be a numeric 2-D array')
    a = np.asarray(x)
    if a.dtype == object or a.dtype == bool or not (np.issubdtype(a.dtype, np.integer) or np.issubdtype(a.dtype, np.floating)):
        raise Refusal(f'{name} must be numeric (integer or float), got {a.dtype}')
    if a.ndim != 2 or a.shape[0] < 1 or a.shape[1] != width:
        raise Refusal(f'{name} must have shape (n_cells, {width}), got {a.shape}')
    a = a.astype(float)
    if not np.isfinite(a).all():
        raise Refusal(f'{name} contains non-finite values')
    return a


def artifact_predict(meta, models, gene_names, position_names, gene_counts, positions):
    """Features-only prediction (no labels accepted). Panel and widths must match the artifact exactly."""
    if list(gene_names) != meta['gene_features'] or list(position_names) != meta['position_features']:
        raise Refusal('feature panel differs from the artifact; refit the frozen recipe on the new panel')
    pos = _numeric_2d('positions', positions, len(meta['position_features']))
    if meta['family'] == 'hgbt_gene_position':
        g = _numeric_2d('gene_counts', gene_counts, len(meta['gene_features']))
        if g.shape[0] != pos.shape[0]:
            raise Refusal('gene_counts and positions have different numbers of rows')
        if (g < 0).any():
            raise Refusal('negative gene counts refused before log1p')
        X = np.hstack([np.log1p(g), pos])
    else:
        if gene_counts is not None:
            raise Refusal('the position-only artifact takes no gene counts')
        X = pos
    out = {}
    for ep, t in meta['targets'].items():
        raw = models[ep].predict_proba(_apply_adapter(t['adapter'], X))[:, 1]
        out[ep] = REF.checked_probability((1 - meta['alpha']) * raw + meta['alpha'] * t['q_fit'], X.shape[0], f'artifact {ep}')
    return out


# ---------------------------------------------------------------- stage 1 execution

def execute_a1(manifest_path=None, out_root=OUT, run_id=A1_RUN_ID, saved_predictions=None, repo=ROOT,
               protocol_path=PROTOCOL):
    real = manifest_path is None
    if real:
        if run_id != A1_RUN_ID or Path(out_root).resolve() != OUT.resolve() or Path(protocol_path).resolve() != PROTOCOL.resolve():
            raise Refusal(f'real A1 run must use run_id {A1_RUN_ID}, the canonical output root and protocol')
        proto = check_common(repo)                                     # before any data file is parsed
        check_a1_pins(proto)
        manifest_path = ROOT / proto['pins']['a1_input_manifest']
        saved_predictions = ROOT / proto['pins']['development_predictions']
    else:
        if json.loads(Path(manifest_path).read_text()).get('synthetic') is not True:
            raise Refusal('non-default manifests must be synthetic')
    run_dir = Path(out_root) / run_id
    _begin(run_dir, provenance(run_id, {'phase': 'A1_EXPOSED_DEVELOPMENT', 'manifest_sha256': sha(manifest_path),
                                        'saved_predictions_sha256': sha(saved_predictions) if saved_predictions else None}),
           protocol_path)
    m, d = D.load(manifest_path)
    pred, raw, fold_of, audit, scores, result = a1_round(m, d, saved_predictions)
    write_json(run_dir / 'result.json', dict(result, exposure='EXPOSED DEVELOPMENT ONLY; not a test',
                                             all_configurations=scores))
    write_json(run_dir / 'fit_audit.json', audit)
    with (run_dir / 'predictions.csv').open('x', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['cell_id', 'group_id', 'endpoint', 'fold', 'method', 'raw_probability', 'probability', 'captured'])
        rows = []
        for ep in pred:
            y = d['y'][ep]
            for i in np.flatnonzero(y >= 0):
                for k in pred[ep]:
                    fam = k.split('|')[0]
                    rows.append((d['cell_id'][i], d['group'][i], ep, int(fold_of[ep][i]), k,
                                 repr(float(raw[ep][fam][i])) if fam in FAMILIES else '',
                                 repr(float(pred[ep][k][i])), int(y[i])))
        rows.sort()
        w.writerows(rows)
    artifacts = {}
    for f in FAMILIES:
        meta, models = deploy_fit(d, f, result['selection'][f]['alpha'], m['endogenous_genes'], m['position_columns'])
        mp_ = save_artifact(run_dir, meta, models)
        saved = json.loads(mp_.read_text())
        artifacts[f] = {'metadata': mp_.name, 'metadata_sha256': sha(mp_), 'model_file': saved['model_file'],
                        'model_sha256': saved['model_sha256']}
    done_sha = finish(run_dir, {'run_id': run_id, 'phase': 'A1_EXPOSED_DEVELOPMENT', 'gate_met': result['gate']['met']})
    sel = {'schema': 'synarch_nonlinear_a1_selection_v1', 'run_id': run_id, 'hgbt': HGBT,
           'alphas': {f: result['selection'][f]['alpha'] for f in FAMILIES},
           'selected_keys': {f: result['selection'][f]['selected'] for f in FAMILIES},
           'gate_met': result['gate']['met'], 'source_sha256': source_hashes(),
           'executing_source': str(Path(__file__).resolve()),
           'protocol_snapshot_sha256': sha(run_dir / 'protocol_snapshot.json'),
           'result_sha256': sha(run_dir / 'result.json'), 'predictions_sha256': sha(run_dir / 'predictions.csv'),
           'run_complete_sha256': done_sha, 'artifacts': artifacts,
           'basis': 'A1 leave-one-animal-out OOF on 2 exposed animals (development-optimistic)'}
    write_json(run_dir / 'selection.json', sel)
    return run_dir, result, scores


# ---------------------------------------------------------------- stage 2: fresh simulation

def generate_case(scenario, seed, anchor, authorized=False):
    if any(seed in s for s in NEW_SEEDS.values()) and not authorized:
        raise Refusal('new prespecified simulation seeds are only generated by the authorized simulation run')
    return SIM.generate(scenario, seed, anchor)


def sim_case(d, recipe, alphas):
    """Comparators via the frozen simulation fitter; HGBT on TRAIN only. TEST labels/oracles never passed."""
    comps, comp_audit = SIM.fit_predict(d, recipe)
    train, test = d['split'] == 'train', d['split'] == 'test'
    genes = np.log1p(d['genes'])
    feats = {'hgbt_gene_position': np.column_stack([genes, d['xyz']]), 'hgbt_position': d['xyz']}
    preds = dict(comps)
    raws, audit = {}, {}
    for f, X in feats.items():
        P = np.empty((test.sum(), len(d['target_names']))); R = np.empty_like(P)
        for ti, t in enumerate(d['target_names']):
            ok = np.zeros(len(train), bool); ok[train] = d['y'][train, ti] >= 0     # TRAIN labels only
            r, s, au = hgbt_fit_predict(X[ok], d['y'][ok, ti], d['animal'][ok], X[test])
            R[:, ti], P[:, ti] = r, s[alphas[f]]
            audit[f'{f}|{t}'] = au
        name = key(f, alphas[f])
        preds[name], raws[name] = P, R
    return preds, raws, {'comparators': comp_audit, 'hgbt': audit}


def check_selection_record(record_path=SELECTION_RECORD, a1_dir=OUT / A1_RUN_ID, repo=ROOT):
    """Committed selection record must match the completed A1 run, its files and the executing source."""
    record_path, a1_dir = Path(record_path), Path(a1_dir)
    if not record_path.is_file() or not B.protocol_committed(record_path, repo):
        raise Refusal('selection record is missing, not committed, or modified')
    rec = json.loads(record_path.read_text())
    sel_path = a1_dir / 'selection.json'
    if not sel_path.is_file() or sha(sel_path) != rec.get('a1_selection_sha256'):
        raise Refusal('selection record does not match the A1 selection.json hash')
    sel = json.loads(sel_path.read_text())
    done_path = a1_dir / 'RUN_COMPLETE.json'
    if not done_path.is_file() or sha(done_path) != sel['run_complete_sha256'] or sha(done_path) != rec.get('a1_run_complete_sha256'):
        raise Refusal('A1 RUN_COMPLETE missing or hash mismatch')
    for name, h in json.loads(done_path.read_text())['files_sha256'].items():
        if sha(a1_dir / name) != h:
            raise Refusal(f'A1 file {name} changed after completion')
    if rec.get('alphas') != sel['alphas'] or set(sel['alphas']) != set(FAMILIES) or \
            any(a not in ALPHAS for a in sel['alphas'].values()):
        raise Refusal('selection record alphas differ from the A1 selection or the grid')
    if sel['source_sha256'] != source_hashes():
        raise Refusal('executing source or imported sources differ from those that produced the A1 selection')
    if sel['hgbt'] != HGBT:
        raise Refusal('fixed HGBT structure differs from the A1 selection')
    return {f: float(a) for f, a in sel['alphas'].items()}


def sim_summary(cases):
    names = sorted(next(iter(cases.values()))['scores'])
    cand = [n for n in names if n.startswith(CANDIDATE)][0]
    out = {'candidate': cand, 'per_scenario': {}, 'overall': {}}

    def agg(sel):
        res = {}
        for c in names:
            if c == cand:
                continue
            ll = [cases[k]['scores'][c]['log_loss_target_equal'] - cases[k]['scores'][cand]['log_loss_target_equal'] for k in sel]
            br = [cases[k]['scores'][c]['brier_target_equal'] - cases[k]['scores'][cand]['brier_target_equal'] for k in sel]
            res[c] = {'log_loss_diff_mean': float(np.mean(ll)), 'log_loss_diff_min': float(np.min(ll)),
                      'log_loss_diff_max': float(np.max(ll)), 'brier_diff_mean': float(np.mean(br)),
                      'brier_diff_min': float(np.min(br)), 'brier_diff_max': float(np.max(br)),
                      'cases_candidate_lower_log_loss': int(sum(v > 0 for v in ll)),
                      'cases_candidate_brier_not_higher': int(sum(v >= 0 for v in br)), 'n_cases': len(sel)}
        return res
    for s in SIM.SCENARIOS:
        out['per_scenario'][s] = agg([k for k in cases if k.startswith(s + '|')])
    out['overall'] = agg(list(cases))
    out['meaning'] = 'control minus candidate across cases; positive = candidate better; no biological claim'
    return out


def save_case_npz(path, d, case_id, seed):
    """allow_pickle=False compatible: only numeric/bool/unicode arrays, plus case_id."""
    arrays = {k: np.asarray(v) for k, v in d.items()}
    arrays['case_id'] = np.full(len(d['cell_id']), case_id)
    arrays['seed'] = np.asarray(seed)
    for k, v in arrays.items():
        if v.dtype == object:
            raise Refusal(f'case array {k} has object dtype; refusing a pickle-dependent NPZ')
    with Path(path).open('xb') as f:
        np.savez_compressed(f, **arrays)
    with np.load(path, allow_pickle=False) as z:                     # roundtrip check
        for k, v in arrays.items():
            if not np.array_equal(z[k], v, equal_nan=v.dtype.kind == 'f'):
                raise Refusal(f'NPZ roundtrip differs for {k}')


def execute_sim(out_root=OUT, run_id=SIM_RUN_ID, repo=ROOT, protocol_path=PROTOCOL):
    if run_id != SIM_RUN_ID or Path(out_root).resolve() != OUT.resolve() or Path(protocol_path).resolve() != PROTOCOL.resolve():
        raise Refusal(f'simulation run must use run_id {SIM_RUN_ID}, the canonical output root and protocol')
    check_common(repo)                                                 # before anything is generated
    alphas = check_selection_record(SELECTION_RECORD, OUT / A1_RUN_ID, repo)
    run_dir = Path(out_root) / run_id
    _begin(run_dir, provenance(run_id, {'phase': 'FRESH_SYNTHETIC', 'frozen_alphas': alphas,
                                        'selection_record_sha256': sha(SELECTION_RECORD),
                                        'seeds': {s: list(v) for s, v in NEW_SEEDS.items()}}), protocol_path)
    anchor, recipe = SIM.anchors()
    write_json(run_dir / 'anchor_receipt.json', anchor)
    cases = {}
    with (run_dir / 'predictions.csv').open('x', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['case_id', 'scenario', 'seed', 'cell_id', 'animal_id', 'target', 'assayed', 'observed_label',
                    'method', 'raw_probability', 'probability'])
        for s in SIM.SCENARIOS:
            for seed in NEW_SEEDS[s]:
                case_id = f'{s}|{seed}'
                d = generate_case(s, seed, anchor, authorized=True)
                save_case_npz(run_dir / f'SYNTHETIC_{s}_{seed}.npz', d, case_id, seed)
                preds, raws, audit = sim_case(d, recipe, alphas)
                sc = SIM.score(d, preds)
                gt = gate(sc, key(CANDIDATE, alphas[CANDIDATE]), A1_TO_SIM[FROZEN_LR], key('hgbt_position', alphas['hgbt_position']),
                          [A1_TO_SIM[c] for c in STRONG_COPY_AVG])
                cases[case_id] = {'scenario': s, 'seed': seed, 'gate': gt,
                                  'scores': {k: {kk: vv for kk, vv in v.items() if kk != 'per_target'} for k, v in sc.items()},
                                  'per_target': {k: v['per_target'] for k, v in sc.items()}}
                write_json(run_dir / f'fit_audit_{s}_{seed}.json', audit)
                ev = np.flatnonzero(d['split'] == 'test')
                for i, row in enumerate(ev):
                    for ti, t in enumerate(d['target_names']):
                        for meth, P in preds.items():
                            w.writerow([case_id, s, seed, d['cell_id'][row], d['animal'][row], t,
                                        bool(d['assayed'][row, ti]), int(d['y'][row, ti]), meth,
                                        repr(float(raws[meth][i, ti])) if meth in raws else '', repr(float(P[i, ti]))])
    write_json(run_dir / 'metrics.json', {'classification': SIM.LABEL, 'cases': cases, 'summary': sim_summary(cases),
                                          'settings_selected_on_these_outcomes': False})
    finish(run_dir, {'run_id': run_id, 'phase': 'FRESH_SYNTHETIC', 'cases': len(cases)})
    return run_dir, cases


if __name__ == '__main__':
    a = sys.argv[1:]
    if a == ['a1']:
        rd, res, _ = execute_a1()
        print('A1 development run complete', rd, sha(rd / 'selection.json'), 'gate met:', res['gate']['met'])
    elif a == ['sim']:
        rd, cases = execute_sim()
        print('simulation run complete', rd, sha(rd / 'RUN_COMPLETE.json'), len(cases), 'cases')
    else:
        sys.exit(__doc__)
