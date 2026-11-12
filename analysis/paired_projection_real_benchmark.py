"""Real-neuron paired projection benchmark runner (separate from the reviewed synthetic-only reference).

Protocol: docs/paired_projection_real_protocol.json (DRAFT_FOR_ROOT_REVIEW_NO_REAL_RUN_AUTHORIZED).

Question: does measured endogenous gene expression + measured position predict an observed captured
projection label, on held-out animals, better than location-only and copying controls?

* Imports pure helpers from analysis/paired_projection_reference.py (reviewed); does not edit or bypass it.
* A real (non-synthetic) manifest is refused before any data file is parsed unless it pins an ADMITTED
  ROOT admission record and this committed, unmodified protocol. Exploratory runs may not carry a power
  approval and never PASS; a two-animal transfer is exploratory only.
* Candidate: L2 logistic (C=1, lbfgs, max_iter 2000, tol 1e-6, seed 20261112) on log1p measured
  markers + measured position (1-3 declared axes), fit-only scaling, group-balanced fit weights.
* Controls: constant 0.5, train frequency, location logistic, location KNN15, molecular KNN15,
  matched-input KNN15, 1NN copy on location, 1NN copy on the candidate's exact representation, and a
  technical logistic only when authentic technical variables exist. Donors come from training animals only.

Usage:
  /opt/anaconda3/bin/python -I analysis/paired_projection_real_benchmark.py synthetic
  /opt/anaconda3/bin/python -I analysis/paired_projection_real_benchmark.py run MANIFEST.json
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis import paired_projection_reference as REF  # noqa: E402  (pure helpers only)

OUT = ROOT / 'data/interim/paired_projection_real_benchmark'
PROTOCOL = ROOT / 'docs/paired_projection_real_protocol.json'
SEED = REF.SEED
CLIP = REF.CLIP
K_NN = 15
LOGIT = dict(REF.LOGIT)
Refusal = REF.Refusal
FILE_KEYS = ['cells', 'genes', 'labels', 'mask']
REQUIRED = ['synthetic', 'dataset_id', 'run_class', 'fold_scheme', 'files', 'endpoint_mode', 'group_unit',
            'groups', 'independence_evidence', 'paired_measurement_evidence', 'per_cell_join',
            'endogenous_genes', 'position_columns', 'position_dims', 'position_source', 'position_measured',
            'technical_columns', 'targets', 'label_semantics', 'frozen_cohort']
INFERRED_WORDS = ('tangram', 'infer', 'impute', 'predict', 'mapped', 'registered_from', 'projected', 'estimated',
                  'nearest_atlas', 'transfer')
BASE_METHODS = ['candidate_logistic', 'constant_0.5', 'train_frequency', 'location_logistic', 'location_knn15',
                'molecular_knn15', 'matched_knn15', 'copy_1nn_location', 'copy_1nn_matched']
LOC_COPY_CONTROLS = ['location_logistic', 'location_knn15', 'copy_1nn_location', 'copy_1nn_matched',
                     'matched_knn15', 'molecular_knn15', 'train_frequency']


def sha256_file(p):
    return REF.sha256_file(p)


# ---------------------------------------------------------------- real-run pins (admission + committed protocol)

def protocol_committed(path, repo=ROOT):
    """True only if `path` is tracked in git and identical to HEAD (no working-tree or index change)."""
    rel = str(Path(path).resolve().relative_to(Path(repo).resolve()))
    def git(*a):
        return subprocess.run(['git', '-C', str(repo), *a], capture_output=True, text=True)
    if git('ls-files', '--error-unmatch', rel).returncode != 0:
        return False
    if git('diff', '--quiet', 'HEAD', '--', rel).returncode != 0:
        return False
    return git('diff', '--cached', '--quiet', '--', rel).returncode == 0


CONFIG_KEYS = ['run_id', 'dataset_id', 'run_class', 'fold_scheme', 'endpoint_mode', 'targets', 'group_unit', 'groups',
               'per_cell_join', 'endogenous_genes', 'position_columns', 'position_dims', 'position_source',
               'position_measured', 'technical_columns', 'label_semantics', 'frozen_cohort']
FROZEN_STATUS = 'FROZEN_EXPLORATORY_AUTHORIZED'


def analysis_config(m):
    """Canonical complete analysis configuration (target order kept; groups sorted by ID)."""
    cfg = {k: m.get(k) if k == 'run_id' else m[k] for k in CONFIG_KEYS}
    cfg['groups'] = sorted(({'group_id': g['group_id'], 'kind': g['kind'], 'animals': list(g['animals'])}
                            for g in m['groups']), key=lambda g: g['group_id'])
    cfg['file_sha256'] = {k: m['files'][k]['sha256'] for k in FILE_KEYS}
    return cfg


def analysis_config_sha256(m):
    return hashlib.sha256(json.dumps(analysis_config(m), sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def check_real_pins(m, base, repo=ROOT, protocol_path=PROTOCOL, runner_path=None, reference_path=None):
    """Admission record + exact config + code hashes + frozen committed protocol. Parses no data file."""
    adm = m.get('admission')
    if not isinstance(adm, dict) or not {'record_path', 'record_sha256'} <= set(adm):
        raise Refusal('real run needs admission.record_path and admission.record_sha256 (ROOT admission)')
    rp = (Path(base) / adm['record_path']).resolve()
    if not rp.is_file() or sha256_file(rp) != adm['record_sha256']:
        raise Refusal('admission record missing or SHA-256 mismatch')
    rec = json.loads(rp.read_text())
    if rec.get('status') != 'ADMITTED' or rec.get('dataset_id') != m['dataset_id']:
        raise Refusal('admission record is not ADMITTED for this dataset_id')
    pinned = rec.get('file_sha256', {})
    for k in FILE_KEYS:
        if pinned.get(k) != m['files'][k]['sha256']:
            raise Refusal(f'admission record does not pin files.{k} at the manifest hash')
    try:
        cfg_sha = analysis_config_sha256(m)
    except (KeyError, TypeError) as e:
        raise Refusal(f'malformed analysis configuration: {e}')
    if rec.get('analysis_config_sha256') != cfg_sha:
        raise Refusal('analysis configuration differs from the one bound in the admission record')
    runner_path = Path(runner_path or __file__)
    reference_path = Path(reference_path or REF.__file__)
    if rec.get('runner_sha256') != sha256_file(runner_path):
        raise Refusal('executing runner SHA-256 differs from the admission record')
    if rec.get('reference_sha256') != sha256_file(reference_path):
        raise Refusal('imported reference module SHA-256 differs from the admission record')
    if rec.get('protocol_sha256') != sha256_file(protocol_path):
        raise Refusal('protocol SHA-256 differs from the admission record')
    pr = m.get('protocol')
    if not isinstance(pr, dict) or pr.get('sha256') != sha256_file(protocol_path):
        raise Refusal('manifest must pin the protocol SHA-256 exactly')
    if not protocol_committed(protocol_path, repo):
        raise Refusal('protocol is not committed in git or has uncommitted changes')
    proto = json.loads(Path(protocol_path).read_text())
    if proto.get('status') != FROZEN_STATUS:
        raise Refusal(f"protocol status {proto.get('status')!r} is not {FROZEN_STATUS}; real runs refused")
    if m['run_class'] != 'exploratory':
        raise Refusal('confirmatory path is disabled; only exploratory real runs are possible')
    return rec


def check_run_class(m):
    rc = m['run_class']
    if rc not in ('confirmatory', 'exploratory'):
        raise Refusal("run_class must be 'confirmatory' or 'exploratory'")
    pa = m.get('power_approval')
    if rc == 'exploratory' and pa is not None:
        raise Refusal('exploratory runs must not carry a power approval (no fictitious power)')
    if rc == 'confirmatory' and not m['synthetic']:
        raise Refusal('confirmatory path is disabled; only exploratory real runs are possible')
    if rc == 'confirmatory':
        if not isinstance(pa, dict) or not REF._nonblank_str(pa.get('approved_by', '')) or \
                not REF._nonblank_str(pa.get('reference', '')) or not isinstance(pa.get('min_groups'), int) or \
                pa['min_groups'] < 3:
            raise Refusal('confirmatory run needs power_approval (approved_by, reference, min_groups >= 3)')
        if len(m['groups']) < pa['min_groups']:
            raise Refusal('fewer groups than the approved power requires')
        if len(m['groups']) <= 2:
            raise Refusal('a two-animal transfer is exploratory only')
    fs = m['fold_scheme']
    if not isinstance(fs, dict) or fs.get('kind') not in ('leave_one_group_out', 'grouped_kfold'):
        raise Refusal("fold_scheme.kind must be 'leave_one_group_out' or 'grouped_kfold'")
    if fs['kind'] == 'grouped_kfold' and (not isinstance(fs.get('k'), int) or not 2 <= fs['k'] <= len(m['groups'])):
        raise Refusal('grouped_kfold needs an integer k between 2 and the number of groups')


# ---------------------------------------------------------------- manifest validation (no data parsed)

def validate_manifest(m, base, repo=ROOT, protocol_path=PROTOCOL, runner_path=None, reference_path=None):
    missing = [k for k in REQUIRED if k not in m]
    if missing:
        raise Refusal(f'manifest missing required fields: {missing}')
    if not isinstance(m['synthetic'], bool):
        raise Refusal('synthetic must be a boolean')
    for k in FILE_KEYS:
        if not isinstance(m['files'], dict) or k not in m['files'] or not {'path', 'sha256'} <= set(m['files'][k]):
            raise Refusal(f'files.{k} with path and sha256 required (explicit assay mask included)')
    check_run_class(m)
    if not m['synthetic']:
        check_real_pins(m, base, repo, protocol_path, runner_path, reference_path)   # before data are touched
    if m['endpoint_mode'] not in REF.ENDPOINT_MODES:
        raise Refusal(f'endpoint_mode must be one of {REF.ENDPOINT_MODES}')
    if m['endpoint_mode'] == 'exclusive_two_target' and len(m['targets']) != 2:
        raise Refusal('exclusive_two_target needs exactly two targets')
    if m['group_unit'] not in ('animal', 'pool'):
        raise Refusal("group_unit must be 'animal' or 'pool'")
    if m['per_cell_join'] != 'exact_cell_id':
        raise Refusal('per_cell_join must be exact_cell_id')
    if m['label_semantics'] != 'captured_label':
        raise Refusal("label_semantics must be 'captured_label'")
    for k in ('independence_evidence', 'paired_measurement_evidence'):
        if not REF._nonblank_str(m[k]):
            raise Refusal(f'{k} must be a non-empty evidence string')
    # position: measured in the same cell, 1-3 declared axes
    if m['position_measured'] is not True or m['position_source'] != 'measured_in_same_cell':
        raise Refusal("position must be measured in the same cell (position_source 'measured_in_same_cell')")
    prov = str(m.get('position_provenance', '')).lower()
    if any(w in prov for w in INFERRED_WORDS):
        raise Refusal('position provenance indicates inferred/imputed/mapped coordinates; refused')
    pos = m['position_columns']
    if not isinstance(m['position_dims'], int) or not 1 <= m['position_dims'] <= 3 or \
            not isinstance(pos, list) or len(pos) != m['position_dims']:
        raise Refusal('position_dims must be 1-3 and equal the number of position_columns')
    # groups
    groups = m['groups']
    if not isinstance(groups, list) or len(groups) < 2:
        raise Refusal('at least two groups required')
    seen_ids, seen_animals = set(), set()
    for g in groups:
        if not {'group_id', 'kind', 'animals', 'independence_evidence'} <= set(g):
            raise Refusal(f'group missing fields: {g}')
        if not REF._nonblank_str(g['group_id']):
            raise Refusal('blank or padded group_id')
        if g['kind'] in REF.BAD_UNITS or g['kind'] != m['group_unit']:
            raise Refusal(f"group {g['group_id']}: kind {g['kind']!r} cannot stand for an animal")
        a = g['animals']
        if not isinstance(a, list) or not a or not all(REF._nonblank_str(x) for x in a) or len(set(a)) != len(a):
            raise Refusal(f"group {g['group_id']}: true animal identities required (non-blank, unpadded, unique)")
        if g['kind'] == 'animal' and len(a) != 1:
            raise Refusal(f"group {g['group_id']}: an animal group lists exactly one animal")
        if not REF._nonblank_str(g['independence_evidence']):
            raise Refusal(f"group {g['group_id']}: independence evidence required")
        if g['group_id'] in seen_ids or seen_animals.intersection(a):
            raise Refusal('groups must have unique IDs and disjoint animals')
        seen_ids.add(g['group_id'])
        seen_animals.update(a)
    # names
    tech = m['technical_columns']
    if not isinstance(tech, list):
        raise Refusal('technical_columns must be a list (empty when no authentic technical variables exist)')
    if tech and not REF._nonblank_str(m.get('technical_evidence', '')):
        raise Refusal('technical columns need technical_evidence (authentic variables only)')
    for role, names, can_be_empty in (('endogenous_genes', m['endogenous_genes'], False),
                                      ('position_columns', pos, False), ('technical_columns', tech, True),
                                      ('targets', m['targets'], False)):
        if not isinstance(names, list) or (not names and not can_be_empty):
            raise Refusal(f'{role} must be a non-empty list')
        if not all(REF._nonblank_str(n) for n in names) or len({n.lower() for n in names}) != len(names):
            raise Refusal(f'{role}: names must be unique, non-empty and unpadded')
    predictors = m['endogenous_genes'] + pos + tech
    if len({n.lower() for n in predictors}) != len(predictors):
        raise Refusal('predictor roles must be disjoint')
    for n in predictors:
        why = REF.forbidden_name(n, m['targets'])
        if why:
            raise Refusal(f'forbidden predictor name {n!r}: {why}')
    fc = m['frozen_cohort']
    if not isinstance(fc, dict) or not fc:
        raise Refusal('frozen_cohort (per endpoint sha256 + n of eligible cell IDs) required before fitting')
    for k in FILE_KEYS:
        p = (Path(base) / m['files'][k]['path']).resolve()
        if not p.is_file() or sha256_file(p) != m['files'][k]['sha256']:
            raise Refusal(f'files.{k} missing or SHA-256 mismatch')
    return True


# ---------------------------------------------------------------- data checks (cells/labels/mask first, genes last)

def read_table(m, base, key):
    import pandas as pd
    return pd.read_csv(Path(base) / m['files'][key]['path'], dtype=str, keep_default_na=False)


def cohort_hash(ids):
    return hashlib.sha256('\n'.join(sorted(ids)).encode()).hexdigest()


def check_light(m, cells, labels, mask):
    """Reuses the reviewed light checks with this runner's position width and group rules."""
    shim = dict(m, power_approval={'min_groups': (m.get('power_approval') or {}).get('min_groups', 2)})
    pos = m['position_columns']
    need = ['cell_id', 'group_id'] + pos + m['technical_columns']
    if list(cells.columns) != need:
        raise Refusal(f'cells columns must be exactly {need}')
    # REF.check_light expects exactly 3 position columns; pad a read-only view only for its column check.
    light = _light_any_dims(shim, cells, labels, mask)
    for ep, ye in light['y'].items():
        ids = light['cell_id'][ye >= 0]
        pin = m['frozen_cohort'].get(ep)
        if not isinstance(pin, dict) or pin.get('n') != len(ids) or pin.get('sha256') != cohort_hash(ids):
            raise Refusal(f'endpoint {ep}: eligible cohort differs from the frozen cohort; refused before fitting')
    if set(m['frozen_cohort']) != set(light['y']):
        raise Refusal('frozen_cohort endpoints differ from the derived endpoints')
    return light


def _light_any_dims(m, cells, labels, mask):
    """Same checks as REF.check_light, for 1-3 position columns (REF fixes 3)."""
    npos = len(m['position_columns'])
    rename = {c: f'__pos{i}' for i, c in enumerate(m['position_columns'])}
    c2 = cells.rename(columns=rename)
    extra = [f'__pad{i}' for i in range(3 - npos)]
    for e in extra:
        c2.insert(2 + npos + extra.index(e), e, '0')
    m2 = dict(m, position_columns=[*rename.values(), *extra])
    light = REF.check_light(m2, c2, labels, mask)
    light['xyz'] = light['xyz'][:, :npos]
    return light


def assemble(m, base):
    light = check_light(m, *(read_table(m, base, k) for k in ('cells', 'labels', 'mask')))
    return REF.check_genes(m, read_table(m, base, 'genes'), light)


# ---------------------------------------------------------------- folds, fitting, prediction

def folds_for(m, groups):
    ug = sorted(set(groups))
    if m['fold_scheme']['kind'] == 'leave_one_group_out':
        return [{g} for g in ug]
    order = np.random.default_rng(SEED).permutation(len(ug))
    k = m['fold_scheme']['k']
    return [set(np.array(ug)[order[i::k]]) for i in range(k)]


def methods_for(d):
    return BASE_METHODS + (['technical'] if d['tech'].shape[1] else [])


def fit_predict(d, tr, ev, y_tr):
    from sklearn.linear_model import LogisticRegression
    from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
    from sklearn.preprocessing import StandardScaler
    w = REF.group_balanced_weights(d['group'][tr])
    full = np.hstack([d['genes'], d['xyz']])
    s_full = StandardScaler().fit(full[tr]); A, B = s_full.transform(full[tr]), s_full.transform(full[ev])
    s_pos = StandardScaler().fit(d['xyz'][tr]); P, Q = s_pos.transform(d['xyz'][tr]), s_pos.transform(d['xyz'][ev])
    s_gen = StandardScaler().fit(d['genes'][tr]); G, H = s_gen.transform(d['genes'][tr]), s_gen.transform(d['genes'][ev])
    n = int(ev.sum())
    out = {
        'candidate_logistic': LogisticRegression(**LOGIT).fit(A, y_tr, sample_weight=w).predict_proba(B)[:, 1],
        'constant_0.5': np.full(n, 0.5),
        'train_frequency': np.full(n, float(np.average(y_tr, weights=w))),
        'location_logistic': LogisticRegression(**LOGIT).fit(P, y_tr, sample_weight=w).predict_proba(Q)[:, 1],
        'location_knn15': KNeighborsClassifier(n_neighbors=K_NN).fit(P, y_tr).predict_proba(Q)[:, 1],
        'molecular_knn15': KNeighborsClassifier(n_neighbors=K_NN).fit(G, y_tr).predict_proba(H)[:, 1],
        'matched_knn15': KNeighborsClassifier(n_neighbors=K_NN).fit(A, y_tr).predict_proba(B)[:, 1],
    }
    train_groups = d['group'][tr]
    donors = {}
    for name, (X, Y) in {'copy_1nn_location': (P, Q), 'copy_1nn_matched': (A, B)}.items():
        nn = NearestNeighbors(n_neighbors=1).fit(X)
        idx = nn.kneighbors(Y, return_distance=False)[:, 0]
        out[name] = y_tr[idx].astype(float)
        donors[name] = train_groups[idx]
    if d['tech'].shape[1]:
        s_t = StandardScaler().fit(d['tech'][tr])
        out['technical'] = LogisticRegression(**LOGIT).fit(s_t.transform(d['tech'][tr]), y_tr, sample_weight=w) \
            .predict_proba(s_t.transform(d['tech'][ev]))[:, 1]
    return out, donors, s_full


def evaluate(m, d):
    methods = methods_for(d)
    folds = folds_for(m, d['group'])
    per_ep, fold_log = {}, []
    for ep in d['endpoints']:
        y_all = d['y'][ep]
        elig = y_all >= 0
        preds = {k: np.full(len(y_all), np.nan) for k in methods}
        fold_of = np.full(len(y_all), -1)
        for fi, eval_groups in enumerate(folds):
            train_groups = set(d['group']) - eval_groups
            REF.assert_disjoint(train_groups, eval_groups)
            tr = elig & np.isin(d['group'], list(train_groups))
            ev = elig & np.isin(d['group'], list(eval_groups))
            if not ev.any():
                continue
            y_tr = y_all[tr]
            if tr.sum() < K_NN or len(set(y_tr)) < 2:
                raise Refusal(f'{ep} fold {fi}: fewer than {K_NN} training cells or one class; FAIL')
            out, donors, s_full = fit_predict(d, tr, ev, y_tr)
            if not np.allclose(s_full.mean_, np.hstack([d['genes'], d['xyz']])[tr].mean(axis=0)):
                raise Refusal('scaler not fitted on training animals only')
            for name, dg in donors.items():
                if np.isin(dg, list(eval_groups)).any():
                    raise Refusal(f'{name}: a donor came from an evaluation animal')
            if set(out) != set(methods):
                raise Refusal(f'{ep} fold {fi}: method set mismatch')
            n_eval = int(ev.sum())
            for k in methods:
                preds[k][ev] = REF.checked_probability(out[k], n_eval, f'{ep} fold {fi} {k}')
            fold_of[ev] = fi
            fold_log.append({'endpoint': ep, 'fold': fi, 'eval_groups': sorted(eval_groups),
                             'n_train': int(tr.sum()), 'n_eval': n_eval})
        for k in methods:
            if not np.isfinite(preds[k][elig]).all():
                raise Refusal(f'{ep}: {k} missing a probability for an eligible cell; FAIL')
        per_ep[ep] = summarise(preds, y_all, elig, d['group'], d['cell_id'], methods)
        per_ep[ep]['_fold_of'] = fold_of
    return {'endpoint_mode': d['endpoint_mode'], 'label_status': d['status'], 'eligible_groups': d['eligible_groups'],
            'methods': methods, 'macro': macro_summary(per_ep, methods), 'technical_control': 'present' if 'technical' in methods else
            'not available (no authentic technical variables declared)',
            'folds': [sorted(f) for f in folds], 'fold_log': fold_log, 'per_endpoint': per_ep}


def per_group_scores(p, y, groups):
    pc = np.clip(p, *CLIP)
    out = {}
    for g in sorted(set(groups)):
        s = groups == g
        out[g] = {'n': int(s.sum()), 'brier': float(np.mean((p[s] - y[s]) ** 2)),
                  'log_loss': float(np.mean(-(y[s] * np.log(pc[s]) + (1 - y[s]) * np.log(1 - pc[s]))))}
    return out


def summarise(preds, y_all, elig, groups, cell_ids, methods):
    y, g = y_all[elig], groups[elig]
    res = {'n_eligible_cells': int(elig.sum()), 'eligible_cell_ids_sha256': cohort_hash(cell_ids[elig]),
           'methods': {}, 'paired_vs_candidate': {}}
    pg = {}
    for k in methods:
        p = preds[k][elig]
        pg[k] = per_group_scores(p, y, g)
        br, ll = REF.group_weighted(p, y, g)
        res['methods'][k] = {'brier_animal_equal': br, 'log_loss_animal_equal': ll, 'per_animal': pg[k],
                             'auc_pooled': REF.auc(p, y), 'calibration': REF.calibration(p, y, g)}
    for c in [x for x in LOC_COPY_CONTROLS if x in methods]:
        rows = {a: {'brier_diff': pg[c][a]['brier'] - pg['candidate_logistic'][a]['brier'],
                    'log_loss_diff': pg[c][a]['log_loss'] - pg['candidate_logistic'][a]['log_loss']} for a in pg[c]}
        bd = np.array([r['brier_diff'] for r in rows.values()])
        ld = np.array([r['log_loss_diff'] for r in rows.values()])
        res['paired_vs_candidate'][c] = {
            'meaning': 'control minus candidate per animal; positive = candidate lower (better); descriptive only',
            'per_animal': rows, 'brier_mean': float(bd.mean()), 'brier_min': float(bd.min()),
            'brier_max': float(bd.max()), 'log_loss_mean': float(ld.mean()), 'log_loss_min': float(ld.min()),
            'log_loss_max': float(ld.max()), 'animals_improved_brier': int((bd > 0).sum()),
            'animals_improved_log_loss': int((ld > 0).sum()), 'n_animals': len(rows)}
    res['_preds'] = preds
    return res


def macro_summary(per_ep, methods):
    """Target-equal macro of animal-equal scores, and per-animal macro across the same targets."""
    eps = list(per_ep)
    animals = sorted({a for e in eps for a in per_ep[e]['methods'][methods[0]]['per_animal']})
    coverage = {a: [e for e in eps if a in per_ep[e]['methods'][methods[0]]['per_animal']] for a in animals}
    out = {'targets': eps, 'n_targets': len(eps),
           'coverage': {a: {'targets': c, 'n_targets': len(c)} for a, c in coverage.items()},
           'complete_coverage': all(len(c) == len(eps) for c in coverage.values()),
           'per_animal_rule': 'mean over the targets for which the animal has eligible cells (coverage declared)',
           'methods': {}, 'paired_vs_candidate': {}}
    for k in methods:
        pa = {a: {'brier': float(np.mean([per_ep[e]['methods'][k]['per_animal'][a]['brier'] for e in coverage[a]])),
                  'log_loss': float(np.mean([per_ep[e]['methods'][k]['per_animal'][a]['log_loss'] for e in coverage[a]]))}
              for a in animals}
        out['methods'][k] = {
            'brier_target_equal': float(np.mean([per_ep[e]['methods'][k]['brier_animal_equal'] for e in eps])),
            'log_loss_target_equal': float(np.mean([per_ep[e]['methods'][k]['log_loss_animal_equal'] for e in eps])),
            'per_animal_macro': pa}
    for c in [x for x in LOC_COPY_CONTROLS if x in methods]:
        pa = {a: {'brier_diff': out['methods'][c]['per_animal_macro'][a]['brier'] -
                                out['methods']['candidate_logistic']['per_animal_macro'][a]['brier'],
                  'log_loss_diff': out['methods'][c]['per_animal_macro'][a]['log_loss'] -
                                   out['methods']['candidate_logistic']['per_animal_macro'][a]['log_loss']}
              for a in animals}
        out['paired_vs_candidate'][c] = {
            'meaning': 'control minus candidate; positive = candidate lower (better); descriptive only',
            'brier_target_equal_diff': float(np.mean([per_ep[e]['paired_vs_candidate'][c]['brier_mean'] for e in eps])),
            'log_loss_target_equal_diff': float(np.mean([per_ep[e]['paired_vs_candidate'][c]['log_loss_mean'] for e in eps])),
            'per_animal_macro_diff': pa,
            'animals_improved_brier': int(sum(v['brier_diff'] > 0 for v in pa.values())),
            'animals_improved_log_loss': int(sum(v['log_loss_diff'] > 0 for v in pa.values())),
            'n_animals': len(animals)}
    return out


def gate(m):
    if m['synthetic']:
        return 'SYNTHETIC_SOFTWARE_ONLY'
    if m['run_class'] == 'exploratory':
        return 'EXPLORATORY_NOT_POWERED: descriptive only, never PASS'
    return 'CONFIRMATORY: apply the committed confirmatory_margins to the paired results (decision recorded by ROOT)'


def run_manifest(path, repo=ROOT, protocol_path=PROTOCOL):
    """In-memory run for synthetic regression only; real runs must use execute() (durable, single run_id)."""
    path = Path(path)
    m = json.loads(path.read_text())
    if m.get('synthetic') is not True:
        raise Refusal('real manifests run only through execute() (durable run_id, provenance, keyed predictions)')
    validate_manifest(m, path.parent, repo, protocol_path)
    d = assemble(m, path.parent)
    res = evaluate(m, d)
    res['gate'] = gate(m)
    res['run_class'] = m['run_class']
    return m, d, res


def strip(res):
    r = {k: v for k, v in res.items() if k != 'per_endpoint'}
    r['per_endpoint'] = {e: {k: v for k, v in pt.items() if not k.startswith('_')} for e, pt in res['per_endpoint'].items()}
    return json.loads(json.dumps(r))


def keyed_predictions(d, res):
    """One row per (cell, endpoint, method) for every eligible cell: replayable out-of-fold predictions."""
    rows = []
    for ep, pt in res['per_endpoint'].items():
        y = d['y'][ep]
        for i in np.flatnonzero(y >= 0):
            for k in res['methods']:
                rows.append((d['cell_id'][i], d['group'][i], ep, int(pt['_fold_of'][i]), k,
                             repr(float(pt['_preds'][k][i])), int(y[i])))
    rows.sort()
    return rows


RUN_ID_OK = __import__('re').compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]{2,80}$')


def execute(manifest_path, out_root=OUT, repo=ROOT, protocol_path=PROTOCOL, runner_path=None, reference_path=None):
    """One durable run: exclusive STARTED marker before fitting; immutable outputs; no silent rerun."""
    manifest_path = Path(manifest_path)
    m = json.loads(manifest_path.read_text())
    run_id = m.get('run_id')
    if not isinstance(run_id, str) or not RUN_ID_OK.match(run_id):
        raise Refusal('manifest needs a fixed run_id ([A-Za-z0-9_.-], 3-81 chars)')
    validate_manifest(m, manifest_path.parent, repo, protocol_path, runner_path, reference_path)
    run_dir = Path(out_root) / run_id
    try:
        run_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        raise Refusal(f'run_id {run_id} already exists (started or finished); duplicate runs are refused')
    adm = m.get('admission') or {}
    prov = {'run_id': run_id, 'synthetic': m['synthetic'], 'run_class': m['run_class'],
            'manifest_sha256': sha256_file(manifest_path),
            'analysis_config_sha256': analysis_config_sha256(m),
            'admission_record_sha256': adm.get('record_sha256'),
            'source_sha256': (json.loads((manifest_path.parent / adm['record_path']).read_text()).get('source_sha256')
                              if adm.get('record_path') else None),
            'file_sha256': {k: m['files'][k]['sha256'] for k in FILE_KEYS},
            'runner_sha256': sha256_file(runner_path or __file__),
            'reference_sha256': sha256_file(reference_path or REF.__file__),
            'protocol_sha256': sha256_file(protocol_path),
            'started_utc': datetime.now(timezone.utc).isoformat()}
    with (run_dir / 'STARTED.json').open('x') as f:
        json.dump(prov, f, indent=2, sort_keys=True)
    d = assemble(m, manifest_path.parent)
    res = evaluate(m, d)
    res['gate'], res['run_class'] = gate(m), m['run_class']
    with (run_dir / 'result.json').open('x') as f:
        json.dump(strip(res), f, indent=2, sort_keys=True, allow_nan=False)
    import csv
    with (run_dir / 'predictions.csv').open('x', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['cell_id', 'group_id', 'endpoint', 'fold', 'method', 'probability', 'captured'])
        w.writerows(keyed_predictions(d, res))
    done = dict(prov, result_sha256=sha256_file(run_dir / 'result.json'),
                predictions_sha256=sha256_file(run_dir / 'predictions.csv'),
                started_sha256=sha256_file(run_dir / 'STARTED.json'), gate=res['gate'],
                finished_utc=datetime.now(timezone.utc).isoformat())
    with (run_dir / 'RUN_COMPLETE.json').open('x') as f:
        json.dump(done, f, indent=2, sort_keys=True)
    return run_dir, res


# ---------------------------------------------------------------- synthetic regression

def make_synthetic(sizes=(140, 260, 180, 320, 200, 240), n_genes=10, dims=3, with_tech=True, seed=SEED,
                   unassayed_t2=True):
    """Synthetic only: additive gene + position signal, unequal animals, T2 unassayed in one animal."""
    import pandas as pd
    rng = np.random.default_rng(seed)
    groups = [f'M{g:02d}' for g in range(len(sizes))]
    genes = [f'Gene{j:02d}' for j in range(n_genes)]
    pos = ['x_um', 'y_um', 'z_um'][:dims]
    tech = ['log_total_counts'] if with_tech else []
    unassayed = {(groups[0], 'T2')} if unassayed_t2 else set()
    rows, grows, labs = [], [], []
    for g, n in zip(groups, sizes):
        off = rng.normal(0, 0.1, 3)
        for c in range(n):
            cid = f'{g}_c{c:04d}'
            xyz = rng.uniform(0, 1, 3) + off
            lam = np.exp(rng.normal(1.0, 0.3, n_genes)) * rng.lognormal(0, 0.3)
            lam[0] *= np.exp(1.2 * rng.normal()); lam[1] *= np.exp(1.2 * rng.normal())
            counts = rng.poisson(lam)
            lg = np.log1p(counts)
            rows.append([cid, g, *(xyz[:dims] * 1000), *([np.log(counts.sum() + 1)] if with_tech else [])])
            grows.append([cid, *counts])
            for tg, (wg, wx) in {'T1': ((1.8, -1.2), (2.5, 0.0)), 'T2': ((-1.0, 1.6), (-2.5, 0.0))}.items():
                if (g, tg) in unassayed:
                    continue
                z = wg[0] * (lg[0] - 1.3) + wg[1] * (lg[1] - 1.3) + wx[0] * (xyz[0] - 0.5) * 2
                labs.append([cid, tg, int(rng.uniform() < 0.8 / (1 + np.exp(-(z - 0.5)))), int(rng.uniform() > 0.03)])
    tables = {'cells': pd.DataFrame(rows, columns=['cell_id', 'group_id', *pos, *tech]),
              'genes': pd.DataFrame(grows, columns=['cell_id'] + genes),
              'labels': pd.DataFrame(labs, columns=['cell_id', 'target', 'captured', 'qc_pass']),
              'mask': pd.DataFrame([[g, t, 0 if (g, t) in unassayed else 1] for g in groups for t in ('T1', 'T2')],
                                   columns=['group_id', 'target', 'assayed'])}
    m = {'synthetic': True, 'dataset_id': 'synthetic_real_runner_v1', 'run_class': 'exploratory',
         'fold_scheme': {'kind': 'leave_one_group_out'}, 'endpoint_mode': 'capture_per_target',
         'group_unit': 'animal',
         'groups': [{'group_id': g, 'kind': 'animal', 'animals': [g], 'independence_evidence': 'synthetic animal'}
                    for g in groups],
         'independence_evidence': 'synthetic software test only',
         'paired_measurement_evidence': 'synthetic: same simulated cell', 'per_cell_join': 'exact_cell_id',
         'endogenous_genes': genes, 'position_columns': pos, 'position_dims': dims,
         'position_source': 'measured_in_same_cell', 'position_measured': True,
         'position_provenance': 'synthetic measured coordinates', 'technical_columns': tech,
         'targets': ['T1', 'T2'], 'label_semantics': 'captured_label'}
    if with_tech:
        m['technical_evidence'] = 'synthetic library size'
    m['frozen_cohort'] = freeze_cohort(m, tables)
    return m, tables


def freeze_cohort(m, tables):
    """Compute per-endpoint eligible-ID pins from cells/labels/mask only (no genes, no fitting)."""
    import pandas as pd
    t = {k: pd.read_csv(__import__('io').StringIO(tables[k].to_csv(index=False)), dtype=str, keep_default_na=False)
         for k in ('cells', 'labels', 'mask')}
    shim = dict(m, power_approval={'min_groups': 2})
    light = _light_any_dims(shim, t['cells'], t['labels'], t['mask'])
    return {ep: {'n': int((ye >= 0).sum()), 'sha256': cohort_hash(light['cell_id'][ye >= 0])}
            for ep, ye in light['y'].items()}


def write_dataset(m, tables, folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    m = json.loads(json.dumps(m))
    m['files'] = {}
    for k in FILE_KEYS:
        p = folder / f'{k}.csv'
        tables[k].to_csv(p, index=False)
        m['files'][k] = {'path': p.name, 'sha256': sha256_file(p)}
    mp = folder / 'manifest.json'
    mp.write_text(json.dumps(m, indent=2, sort_keys=True))
    return mp


def synthetic():
    syn = OUT / 'synthetic'
    m, t = make_synthetic()
    _, d, res = run_manifest(write_dataset(m, t, syn / 'clean'))
    m2, t2 = make_synthetic(dims=2, with_tech=False)
    _, d2, res2 = run_manifest(write_dataset(m2, t2, syn / 'two_axes_no_tech'))
    pv = res['per_endpoint']['T1']['paired_vs_candidate']
    checks = {
        'candidate_beats_every_location_copy_control_T1_brier': all(v['brier_mean'] > 0 for v in pv.values()),
        'candidate_beats_every_location_copy_control_T2_brier':
            all(v['brier_mean'] > 0 for v in res['per_endpoint']['T2']['paired_vs_candidate'].values()),
        'leave_one_animal_out_folds': len(res['folds']) == 6 and all(len(f) == 1 for f in res['folds']),
        'technical_absent_when_not_declared': 'technical' not in res2['methods'],
        'two_axis_position_supported': d2['xyz'].shape[1] == 2,
        'gate_synthetic_only': res['gate'] == 'SYNTHETIC_SOFTWARE_ONLY',
        'target_equal_macro_beats_every_location_copy_control':
            all(v['brier_target_equal_diff'] > 0 for v in res['macro']['paired_vs_candidate'].values()),
        'coverage_declared_incomplete_for_M00': res['macro']['coverage']['M00']['n_targets'] == 1,
    }
    receipt = {'schema': 'synarch_paired_projection_real_runner_synthetic_receipt_v2',
               'meaning': 'SYNTHETIC software regression only; no real data loaded or fitted.',
               'software_status': 'PASS' if all(checks.values()) else 'FAIL', 'checks': checks,
               'result_3d_with_tech': strip(res), 'result_2d_no_tech': strip(res2),
               'code_sha256': sha256_file(__file__), 'protocol_sha256': sha256_file(PROTOCOL),
               'reference_code_sha256': sha256_file(REF.__file__),
               'finished_utc': datetime.now(timezone.utc).isoformat()}
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / 'synthetic_receipt.json'
    p.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False))
    print('software', receipt['software_status'], checks)
    for ep in ('T1', 'T2'):
        print(ep, {k: round(v['brier_animal_equal'], 4) for k, v in res['per_endpoint'][ep]['methods'].items()})
    print('wrote', p, sha256_file(p))
    return receipt


if __name__ == '__main__':
    a = sys.argv[1:]
    if a == ['synthetic']:
        synthetic()
    elif len(a) == 2 and a[0] == 'run':
        run_dir, res = execute(a[1])                    # refuses without pins; refuses duplicate run_id
        print('run complete', run_dir, sha256_file(run_dir / 'RUN_COMPLETE.json'), res['gate'])
    else:
        sys.exit(__doc__)
