"""Bounded published-reference benchmark: captured projection label from measured genes + measured position.

Protocol: docs/paired_projection_reference_protocol.json (DRAFT_FOR_ROOT_REVIEW_SOFTWARE_ONLY).

* One candidate: L2 logistic regression (C=1, lbfgs, max_iter 2000, tol 1e-6, seed 20261112; L2 is the
  sklearn 1.8 default, so `penalty` is not passed) on log1p endogenous genes + measured XYZ, with a
  StandardScaler fitted on training cells only.
* Controls: constant 0.5, fit frequency, technical logistic, location KNN (k=15), molecular KNN (k=15),
  Gaussian Naive Bayes. Every method scores the identical eligible query cells in identical folds.
* Group-balanced fit weights (1/n_group, mean-normalised to 1 over the training cells) are used by the
  candidate, fit frequency, technical and Naive Bayes. KNN has no sample weights (stated, not hidden).
  No automatic class_weight='balanced'.
* Explicit endpoint mode:
    - 'capture_per_target' (e.g. PFC cells, captured LHA label): y=1 captured, y=0 only for a QC-passed
      nondetection in an assayed group; failed QC and unassayed rows are unknown and excluded.
    - 'exclusive_two_target' (e.g. Rossi LHb-only vs VTA-only): among QC-passed cells assayed for both
      targets, y=1 first-target-only, y=0 second-target-only; double-labelled and neither stay separate
      statuses and are excluded from the endpoint.
* A strict manifest is validated BEFORE any matrix is read. No bypass flag, no fallback, and any missing
  probability for any method is a FAIL (Refusal), never a dropped row.
* The biological gate is always NOT_EVALUATED. A synthetic PASS means only that the software works.

Usage:
  /opt/anaconda3/bin/python -I analysis/paired_projection_reference.py synthetic
  (real-data fitting is not authorized; there is deliberately no real-data command)
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/interim/paired_projection_reference'
SEED = 20261112
CLIP = (1e-6, 1 - 1e-6)
K_NN = 15
MAX_FOLDS = 5
LOGIT = dict(C=1.0, solver='lbfgs', max_iter=2000, tol=1e-6, random_state=SEED)
ENDPOINT_MODES = ('capture_per_target', 'exclusive_two_target')
REQUIRED = ['synthetic', 'dataset_id', 'files', 'endpoint_mode', 'group_unit', 'groups',
            'independence_evidence', 'paired_measurement_evidence', 'per_cell_join', 'endogenous_genes',
            'position_columns', 'position_measured', 'technical_columns', 'targets', 'power_approval',
            'label_semantics']
FILE_KEYS = ['cells', 'genes', 'labels', 'mask']
BAD_UNITS = {'slice', 'section', 'hemisphere', 'slide', 'well', 'reaction', 'library', 'lane', 'cell'}
# Exact reporter/vector aliases, matched against whole name tokens (so Gfpt1, Gfpt2 and Creb1 are admitted).
REPORTER_ALIASES = {'gfp', 'egfp', 'yfp', 'eyfp', 'cfp', 'ecfp', 'rfp', 'tdtomato', 'tomato', 'dtomato', 'mcherry',
                    'cherry', 'mscarlet', 'scarlet', 'tdt', 'cre', 'icre', 'flp', 'flpo', 'dre', 'wpre', 'sv40',
                    'bgh', 'tva', 'rabies', 'aav', 'raav', 'retroaav', 'cav2', 'hsv'}
# Reserved metadata words, matched as substrings (never gene symbols).
RESERVED_SUBSTRINGS = ('target', 'label', 'captur', 'class', 'cluster', 'projection', 'barcode', 'reporter',
                       'vector', 'qc_pass', 'assay')
RESERVED_COLUMNS = {'cell_id', 'group_id', 'target', 'captured', 'qc_pass', 'assayed'}
METHODS = ['candidate_logistic', 'constant_0.5', 'fit_frequency', 'technical', 'location_knn',
           'molecular_knn', 'naive_bayes']


def forbidden_name(name, targets=()):
    """Return the reason a predictor name is forbidden, or None."""
    low = str(name).lower()
    tokens = [t for t in re.split(r'[^0-9a-z]+', low) if t]
    if any(t in REPORTER_ALIASES for t in tokens):
        return 'reporter/vector alias'
    if any(w in low for w in RESERVED_SUBSTRINGS) or low in RESERVED_COLUMNS:
        return 'reserved metadata name'
    if low in {str(t).lower() for t in targets}:
        return 'target name'
    return None


def _nonblank_str(x):
    return isinstance(x, str) and x.strip() != '' and x == x.strip()


class Refusal(ValueError):
    """Raised for any admission, join, mask, leakage or protocol violation. Never caught to fall back."""


def sha256_file(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


# ---------------------------------------------------------------- manifest (no matrix reads here)

def validate_manifest(m, base):
    """Structural admission check. Reads no matrix; only hashes declared files at the very end."""
    missing = [k for k in REQUIRED if k not in m]
    if missing:
        raise Refusal(f'manifest missing required fields: {missing}')
    if m['synthetic'] is not True:
        raise Refusal('software-only release: non-synthetic manifests are refused before any data file is read')
    if m['endpoint_mode'] not in ENDPOINT_MODES:
        raise Refusal(f"endpoint_mode must be one of {ENDPOINT_MODES}")
    if m['endpoint_mode'] == 'exclusive_two_target' and len(m['targets']) != 2:
        raise Refusal('exclusive_two_target needs exactly two targets (positive first)')
    if m['group_unit'] not in ('animal', 'pool'):
        raise Refusal(f"group_unit {m['group_unit']!r} refused; only 'animal' or 'pool'")
    if m['per_cell_join'] != 'exact_cell_id':
        raise Refusal('per_cell_join must be exact_cell_id')
    if m['position_measured'] is not True:
        raise Refusal('position must be measured (position_measured true)')
    if m['label_semantics'] != 'captured_label':
        raise Refusal("label_semantics must be 'captured_label'")
    for k in ('independence_evidence', 'paired_measurement_evidence'):
        if not isinstance(m[k], str) or not m[k].strip():
            raise Refusal(f'{k} must be a non-empty evidence string')
    groups = m['groups']
    if not isinstance(groups, list) or not groups:
        raise Refusal('groups must be a non-empty list')
    seen_ids, seen_animals = set(), set()
    for g in groups:
        for k in ('group_id', 'kind', 'animals', 'independence_evidence'):
            if k not in g:
                raise Refusal(f'group missing {k}: {g}')
        if g['kind'] in BAD_UNITS or g['kind'] not in ('animal', 'pool'):
            raise Refusal(f"group {g['group_id']}: kind {g['kind']!r} cannot stand for an animal")
        if g['kind'] != m['group_unit']:
            raise Refusal(f"group {g['group_id']}: kind differs from group_unit")
        animals = g['animals']
        if not isinstance(animals, list) or not animals or not all(_nonblank_str(a) for a in animals):
            raise Refusal(f"group {g['group_id']}: true animal identities required (non-blank, unpadded)")
        if len(set(animals)) != len(animals):
            raise Refusal(f"group {g['group_id']}: repeated animal identity")
        if g['kind'] == 'animal' and len(animals) != 1:
            raise Refusal(f"group {g['group_id']}: an animal group must list exactly one animal")
        if not str(g['independence_evidence']).strip():
            raise Refusal(f"group {g['group_id']}: independence evidence required")
        if g['group_id'] in seen_ids:
            raise Refusal(f"duplicate group_id {g['group_id']}")
        overlap = seen_animals.intersection(animals)
        if overlap:
            raise Refusal(f'animals appear in more than one group (not disjoint): {sorted(overlap)}')
        seen_ids.add(g['group_id'])
        seen_animals.update(animals)
    pa = m['power_approval']
    if not isinstance(pa, dict) or not str(pa.get('approved_by', '')).strip() or not str(pa.get('reference', '')).strip():
        raise Refusal('independent power approval (approved_by, reference) required')
    if not isinstance(pa.get('min_groups'), int) or pa['min_groups'] < 2:
        raise Refusal('power_approval.min_groups must be an integer >= 2')
    if len(groups) < pa['min_groups']:
        raise Refusal(f"only {len(groups)} groups; power approval needs {pa['min_groups']}")
    roles = {'endogenous_genes': m['endogenous_genes'], 'position_columns': m['position_columns'],
             'technical_columns': m['technical_columns'], 'targets': m['targets']}
    for role, names in roles.items():
        if not isinstance(names, list) or not names:
            raise Refusal(f'{role} must be a non-empty list')
        if not all(_nonblank_str(n) for n in names):
            raise Refusal(f'{role}: names must be non-empty strings without surrounding spaces')
        if len({n.lower() for n in names}) != len(names):
            raise Refusal(f'{role}: duplicate names')
    if len(m['position_columns']) != 3:
        raise Refusal('exactly three measured position columns required')
    predictors = m['endogenous_genes'] + m['position_columns'] + m['technical_columns']
    if len({n.lower() for n in predictors}) != len(predictors):
        raise Refusal('predictor roles must be disjoint (a name appears in two roles)')
    for name in predictors:
        why = forbidden_name(name, m['targets'])
        if why:
            raise Refusal(f'forbidden predictor name {name!r}: {why}')
    for g in groups:
        if not _nonblank_str(g['group_id']):
            raise Refusal('blank or padded group_id')
    files = m['files']
    for k in FILE_KEYS:
        if k not in files or 'path' not in files[k] or 'sha256' not in files[k]:
            raise Refusal(f'files.{k} with path and sha256 required (explicit assay mask included)')
    for k in FILE_KEYS:
        p = (Path(base) / files[k]['path']).resolve()
        if not p.is_file():
            raise Refusal(f'files.{k} not found: {p}')
        if sha256_file(p) != files[k]['sha256']:
            raise Refusal(f'files.{k} sha256 mismatch')
    return True


# ---------------------------------------------------------------- data load, exact joins, endpoint

def read_table(m, base, key):
    """The only place data files are parsed (spied on by tests)."""
    import pandas as pd
    return pd.read_csv(Path(base) / m['files'][key]['path'], dtype=str, keep_default_na=False)


def _ids_ok(series, what):
    vals = series.tolist()
    if not all(_nonblank_str(v) for v in vals):
        raise Refusal(f'missing, blank or padded {what}')


def _int01(df, col):
    if not df[col].isin(['0', '1']).all():
        raise Refusal(f'{col} must be 0/1')
    return df[col].astype(int)


def derive_endpoint(mode, targets, cap, qc, assayed_by_cell):
    """cap/qc: dict target -> int array (-1 = no label row). Returns y (1/0/-1) and status counts.

    -1 = unknown/excluded. captured=0 counts as y=0 only when the row is assayed and QC-passed."""
    n = len(next(iter(cap.values())))
    status = {}
    if mode == 'capture_per_target':
        ys = {}
        for tg in targets:
            ok = assayed_by_cell[tg] & (qc[tg] == 1)
            y = np.where(ok, cap[tg], -1)
            ys[tg] = y
            status[tg] = {'captured': int((y == 1).sum()), 'qc_passed_nondetection': int((y == 0).sum()),
                          'unknown_unassayed': int((~assayed_by_cell[tg]).sum()),
                          'unknown_failed_qc': int((assayed_by_cell[tg] & (qc[tg] != 1)).sum())}
        return ys, status
    a, b = targets
    ok = assayed_by_cell[a] & assayed_by_cell[b] & (qc[a] == 1) & (qc[b] == 1)
    a_only = ok & (cap[a] == 1) & (cap[b] == 0)
    b_only = ok & (cap[a] == 0) & (cap[b] == 1)
    y = np.full(n, -1)
    y[a_only], y[b_only] = 1, 0
    name = f'{a}_only_vs_{b}_only'
    status[name] = {f'{a}_only': int(a_only.sum()), f'{b}_only': int(b_only.sum()),
                    'double_labelled_excluded': int((ok & (cap[a] == 1) & (cap[b] == 1)).sum()),
                    'neither_excluded': int((ok & (cap[a] == 0) & (cap[b] == 0)).sum()),
                    'unknown_unassayed_or_failed_qc': int((~ok).sum())}
    return {name: y}, status


def check_light(m, cells, labels, mask):
    """Cells/labels/mask checks, endpoint derivation and per-endpoint eligible group counts.

    Runs BEFORE the gene matrix is read."""
    pos, tech, targets = m['position_columns'], m['technical_columns'], m['targets']
    need = ['cell_id', 'group_id'] + pos + tech
    if list(cells.columns) != need:
        raise Refusal(f'cells columns must be exactly {need}')
    _ids_ok(cells['cell_id'], 'cell_id in cells'); _ids_ok(cells['group_id'], 'group_id in cells')
    if cells['cell_id'].duplicated().any():
        raise Refusal('duplicate cell_id in cells table (exact join refused)')
    # Canonical order: stable sort by validated cell_id, so every later step is invariant to CSV row order.
    cells = cells.sort_values('cell_id', kind='mergesort').reset_index(drop=True)
    gids = {g['group_id'] for g in m['groups']}
    if not set(cells['group_id']).issubset(gids):
        raise Refusal('cells reference undeclared groups')
    if gids - set(cells['group_id']):
        raise Refusal('declared groups without cells')
    try:
        num = cells[pos + tech].astype(float).to_numpy()
    except ValueError:
        raise Refusal('non-numeric position or technical values')
    if not np.isfinite(num).all():
        raise Refusal('missing or non-finite position/technical values')
    if list(labels.columns) != ['cell_id', 'target', 'captured', 'qc_pass']:
        raise Refusal('labels columns must be cell_id,target,captured,qc_pass')
    _ids_ok(labels['cell_id'], 'cell_id in labels'); _ids_ok(labels['target'], 'target in labels')
    if labels.duplicated(['cell_id', 'target']).any():
        raise Refusal('duplicate (cell_id, target) labels')
    labels = labels.assign(captured=_int01(labels, 'captured'), qc_pass=_int01(labels, 'qc_pass'))
    if not set(labels['target']).issubset(targets):
        raise Refusal('labels contain undeclared targets')
    if not set(labels['cell_id']).issubset(set(cells['cell_id'])):
        raise Refusal('labels reference unknown cells')
    if list(mask.columns) != ['group_id', 'target', 'assayed']:
        raise Refusal('mask columns must be group_id,target,assayed')
    _ids_ok(mask['group_id'], 'group_id in mask'); _ids_ok(mask['target'], 'target in mask')
    mask = mask.assign(assayed=_int01(mask, 'assayed'))
    if mask.duplicated(['group_id', 'target']).any() or len(mask) != len(gids) * len(targets) or \
            set(zip(mask['group_id'], mask['target'])) != {(g, tg) for g in gids for tg in targets}:
        raise Refusal('mask must cover every group x target exactly once')
    assayed = {(r.group_id, r.target) for r in mask.itertuples() if r.assayed == 1}
    group_arr = cells['group_id'].to_numpy(str)
    idx = {c: i for i, c in enumerate(cells['cell_id'])}
    cap = {tg: np.full(len(cells), -1, int) for tg in targets}
    qc = {tg: np.full(len(cells), -1, int) for tg in targets}
    for r in labels.itertuples():
        i = idx[r.cell_id]
        if (group_arr[i], r.target) not in assayed:
            raise Refusal(f'label present for unassayed group x target ({group_arr[i]}, {r.target})')
        cap[r.target][i], qc[r.target][i] = int(r.captured), int(r.qc_pass)
    assayed_by_cell = {tg: np.array([(g, tg) in assayed for g in group_arr]) for tg in targets}
    for tg in targets:
        if (assayed_by_cell[tg] & (cap[tg] < 0)).any():
            raise Refusal(f'assayed cell lacks a label row for {tg}')
    y, status = derive_endpoint(m['endpoint_mode'], targets, cap, qc, assayed_by_cell)
    need_groups = m['power_approval']['min_groups']
    eligible_groups = {}
    for ep, ye in y.items():
        n = len(set(group_arr[ye >= 0]))
        eligible_groups[ep] = n
        if n < need_groups:
            raise Refusal(f'endpoint {ep}: only {n} eligible biological groups; power approval needs {need_groups}')
    return {'cell_id': cells['cell_id'].to_numpy(str), 'group': group_arr, 'xyz': num[:, :3], 'tech': num[:, 3:],
            'y': y, 'endpoints': list(y), 'status': status, 'eligible_groups': eligible_groups,
            'endpoint_mode': m['endpoint_mode']}


def check_genes(m, genes, light):
    glist = m['endogenous_genes']
    if list(genes.columns) != ['cell_id'] + list(glist):
        raise Refusal('gene columns must equal cell_id + declared endogenous genes exactly')
    _ids_ok(genes['cell_id'], 'cell_id in genes')
    if genes['cell_id'].duplicated().any():
        raise Refusal('duplicate cell_id in genes table (exact join refused)')
    if set(genes['cell_id']) != set(light['cell_id']):
        raise Refusal('cells and genes tables do not contain identical cell_id sets')
    try:
        gx = genes.set_index('cell_id').loc[light['cell_id'], glist].astype(float).to_numpy()
    except ValueError:
        raise Refusal('non-numeric gene values')
    if not np.isfinite(gx).all():
        raise Refusal('missing or non-finite gene values')
    if (gx < 0).any():
        raise Refusal('negative gene counts')
    return dict(light, genes=np.log1p(gx))


# ---------------------------------------------------------------- folds, models, metrics

def group_folds(groups):
    ug = np.array(sorted(set(groups)))
    order = np.random.default_rng(SEED).permutation(len(ug))
    k = min(MAX_FOLDS, len(ug))
    return [set(ug[order[i::k]]) for i in range(k)]


def assert_disjoint(train_groups, eval_groups):
    if set(train_groups) & set(eval_groups):
        raise Refusal('group leakage: training and evaluation share a group')


def group_balanced_weights(groups):
    """1/n_group per cell, normalised to mean 1 over the given (training) cells."""
    _, inv, counts = np.unique(groups, return_inverse=True, return_counts=True)
    w = 1.0 / counts[inv]
    return w / w.mean()


def fit_fold(d, tr, y_tr):
    """Fit every method on training rows only. Returns fitted objects (no evaluation data touched)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.naive_bayes import GaussianNB
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.preprocessing import StandardScaler
    w = group_balanced_weights(d['group'][tr])
    full = np.hstack([d['genes'], d['xyz']])
    f = {'weights': w}
    f['sc_full'] = StandardScaler().fit(full[tr])
    f['candidate'] = LogisticRegression(**LOGIT).fit(f['sc_full'].transform(full[tr]), y_tr, sample_weight=w)
    f['frequency'] = float(np.average(y_tr, weights=w))
    f['sc_tech'] = StandardScaler().fit(d['tech'][tr])
    f['technical'] = LogisticRegression(**LOGIT).fit(f['sc_tech'].transform(d['tech'][tr]), y_tr, sample_weight=w)
    f['sc_xyz'] = StandardScaler().fit(d['xyz'][tr])
    f['loc_knn'] = KNeighborsClassifier(n_neighbors=K_NN).fit(f['sc_xyz'].transform(d['xyz'][tr]), y_tr)
    f['sc_gene'] = StandardScaler().fit(d['genes'][tr])
    gtr = f['sc_gene'].transform(d['genes'][tr])
    f['mol_knn'] = KNeighborsClassifier(n_neighbors=K_NN).fit(gtr, y_tr)
    f['nb'] = GaussianNB().fit(gtr, y_tr, sample_weight=w)
    return f


def predict_fold(d, f, ev):
    full = np.hstack([d['genes'], d['xyz']])
    n = int(ev.sum())
    g = f['sc_gene'].transform(d['genes'][ev])
    return {
        'candidate_logistic': f['candidate'].predict_proba(f['sc_full'].transform(full[ev]))[:, 1],
        'constant_0.5': np.full(n, 0.5),
        'fit_frequency': np.full(n, f['frequency']),
        'technical': f['technical'].predict_proba(f['sc_tech'].transform(d['tech'][ev]))[:, 1],
        'location_knn': f['loc_knn'].predict_proba(f['sc_xyz'].transform(d['xyz'][ev]))[:, 1],
        'molecular_knn': f['mol_knn'].predict_proba(g)[:, 1],
        'naive_bayes': f['nb'].predict_proba(g)[:, 1],
    }


def checked_probability(p, n_eval, where):
    """Accept only a float vector of shape (n_eval,), finite and within [0, 1]. No broadcasting."""
    if not isinstance(p, np.ndarray) or p.ndim != 1 or p.shape != (n_eval,):
        raise Refusal(f'{where}: probability must be a vector of shape ({n_eval},), got '
                      f'{type(p).__name__} {getattr(p, "shape", None)}')
    if not np.issubdtype(p.dtype, np.floating) and not np.issubdtype(p.dtype, np.integer):
        raise Refusal(f'{where}: non-numeric probabilities')
    if not np.isfinite(p).all():
        raise Refusal(f'{where}: non-finite probability')
    if (p < 0).any() or (p > 1).any():
        raise Refusal(f'{where}: probability outside [0, 1]')
    return p.astype(float)


def group_weighted(p, y, groups):
    """Per-group mean Brier and clipped log loss, then the unweighted mean over groups."""
    pc = np.clip(p, *CLIP)
    br, ll = [], []
    for g in sorted(set(groups)):
        s = groups == g
        br.append(np.mean((p[s] - y[s]) ** 2))
        ll.append(np.mean(-(y[s] * np.log(pc[s]) + (1 - y[s]) * np.log(1 - pc[s]))))
    return float(np.mean(br)), float(np.mean(ll))


def auc(p, y):
    from scipy.stats import rankdata
    pos, n1 = y == 1, int((y == 1).sum())
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return None
    r = rankdata(p)
    return float((r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def calibration(p, y, groups, bins=10):
    edges = np.linspace(0, 1, bins + 1)
    ids = np.clip(np.digitize(p, edges[1:-1]), 0, bins - 1)
    table, ece = [], 0.0
    for b in range(bins):
        s = ids == b
        if s.any():
            table.append({'bin': b, 'n': int(s.sum()), 'mean_p': float(p[s].mean()), 'observed': float(y[s].mean())})
            ece += float(s.mean()) * abs(float(p[s].mean()) - float(y[s].mean()))
    per_group = {g: {'mean_p': float(p[groups == g].mean()), 'observed': float(y[groups == g].mean()),
                     'n': int((groups == g).sum())} for g in sorted(set(groups))}
    return {'bins': table, 'ece_cell_weighted': float(ece), 'per_group': per_group}


def evaluate(d):
    folds = group_folds(d['group'])
    per_ep, fold_log = {}, []
    for ep in d['endpoints']:
        y_all = d['y'][ep]
        eligible = y_all >= 0
        preds = {mth: np.full(len(y_all), np.nan) for mth in METHODS}
        for fi, eval_groups in enumerate(folds):
            train_groups = set(d['group']) - eval_groups
            assert_disjoint(train_groups, eval_groups)
            tr = eligible & np.isin(d['group'], list(train_groups))
            ev = eligible & np.isin(d['group'], list(eval_groups))
            if not ev.any():
                continue
            y_tr = y_all[tr]
            if tr.sum() < K_NN or len(set(y_tr)) < 2:
                raise Refusal(f'{ep} fold {fi}: fewer than {K_NN} training cells or one class; FAIL, no fallback')
            f = fit_fold(d, tr, y_tr)
            if not np.allclose(f['sc_full'].mean_, np.hstack([d['genes'], d['xyz']])[tr].mean(axis=0)):
                raise Refusal('scaler was not fitted on training cells only')
            out = predict_fold(d, f, ev)
            if set(out) != set(METHODS):
                raise Refusal(f'{ep} fold {fi}: methods returned {sorted(out)}, expected {METHODS}')
            n_eval = int(ev.sum())
            for mth in METHODS:
                preds[mth][ev] = checked_probability(out[mth], n_eval, f'{ep} fold {fi} {mth}')
            fold_log.append({'endpoint': ep, 'fold': fi, 'eval_groups': sorted(eval_groups),
                             'n_train': int(tr.sum()), 'n_eval': int(ev.sum())})
        for mth in METHODS:
            if np.isnan(preds[mth][eligible]).any() or not np.isfinite(preds[mth][eligible]).all():
                raise Refusal(f'{ep}: method {mth} is missing a probability for an eligible cell; FAIL')
        q = eligible
        res = {'n_eligible_cells': int(q.sum()), 'eligible_cell_ids_sha256':
               hashlib.sha256('\n'.join(sorted(d['cell_id'][q])).encode()).hexdigest(), 'methods': {}}
        for mth in METHODS:
            p = preds[mth][q]
            br, ll = group_weighted(p, y_all[q], d['group'][q])
            res['methods'][mth] = {'brier_group_weighted': br, 'log_loss_group_weighted': ll,
                                   'auc_pooled': auc(p, y_all[q]),
                                   'calibration': calibration(p, y_all[q], d['group'][q])}
        res['_preds'] = preds
        per_ep[ep] = res
    macro = {mth: {'brier': float(np.mean([per_ep[e]['methods'][mth]['brier_group_weighted'] for e in d['endpoints']])),
                   'log_loss': float(np.mean([per_ep[e]['methods'][mth]['log_loss_group_weighted'] for e in d['endpoints']]))}
             for mth in METHODS}
    return {'endpoint_mode': d['endpoint_mode'], 'label_status': d['status'],
            'folds': [sorted(f) for f in folds], 'fold_log': fold_log, 'per_endpoint': per_ep, 'macro': macro,
            'knn_weighting_note': 'KNN controls are unweighted (no sample weights in KNeighborsClassifier)',
            'biological_gate': 'NOT_EVALUATED: no committed biological protocol, margins or power'}


def run_manifest(path):
    path = Path(path)
    m = json.loads(path.read_text())
    validate_manifest(m, path.parent)          # refuses before any data file is parsed
    light = check_light(m, *(read_table(m, path.parent, k) for k in ('cells', 'labels', 'mask')))
    d = check_genes(m, read_table(m, path.parent, 'genes'), light)   # gene matrix read last
    return m, d, evaluate(d)


# ---------------------------------------------------------------- synthetic software test

def make_synthetic(sizes=(120, 450, 200, 380, 160, 300, 240, 420), n_genes=12, seed=SEED):
    """Planted additive gene + position signal, unequal group sizes, imperfect detection, failed-QC rows,
    and target T2 unassayed in two groups."""
    import pandas as pd
    rng = np.random.default_rng(seed)
    rows, grows, labs = [], [], []
    groups = [f'A{g:02d}' for g in range(len(sizes))]
    genes = [f'Gene{j:02d}' for j in range(n_genes)]
    unassayed = {(groups[0], 'T2'), (groups[1], 'T2')}
    for g, n in zip(groups, sizes):
        off = rng.normal(0, 0.15, 3)
        for c in range(n):
            cid = f'{g}_c{c:04d}'
            xyz = rng.uniform(0, 1, 3) + off
            lib = rng.lognormal(0, 0.3)
            lam = np.exp(rng.normal(1.0, 0.3, n_genes)) * lib
            lam[0] *= np.exp(1.2 * rng.normal())
            lam[1] *= np.exp(1.2 * rng.normal())
            counts = rng.poisson(lam)
            lg = np.log1p(counts)
            rows.append([cid, g, *(xyz * 1000), np.log(counts.sum() + 1)])
            grows.append([cid, *counts])
            for tg, (wg, wx) in {'T1': ((1.8, -1.2), (2.5, 0.0)), 'T2': ((-1.0, 1.6), (0.0, -2.5))}.items():
                if (g, tg) in unassayed:
                    continue
                z = wg[0] * (lg[0] - 1.3) + wg[1] * (lg[1] - 1.3) + wx[0] * (xyz[0] - 0.5) * 2 + wx[1] * (xyz[1] - 0.5) * 2
                p_capture = 0.8 / (1 + np.exp(-(z - 0.5)))
                labs.append([cid, tg, int(rng.uniform() < p_capture), int(rng.uniform() > 0.03)])
    cells = pd.DataFrame(rows, columns=['cell_id', 'group_id', 'x_um', 'y_um', 'z_um', 'log_total_counts'])
    gtab = pd.DataFrame(grows, columns=['cell_id'] + genes)
    labels = pd.DataFrame(labs, columns=['cell_id', 'target', 'captured', 'qc_pass'])
    mask = pd.DataFrame([[g, tg, 0 if (g, tg) in unassayed else 1] for g in groups for tg in ('T1', 'T2')],
                        columns=['group_id', 'target', 'assayed'])
    manifest = {
        'synthetic': True, 'dataset_id': 'synthetic_planted_v2', 'endpoint_mode': 'capture_per_target',
        'group_unit': 'animal',
        'groups': [{'group_id': g, 'kind': 'animal', 'animals': [g], 'independence_evidence': 'synthetic: one simulated animal'} for g in groups],
        'independence_evidence': 'synthetic software test only',
        'paired_measurement_evidence': 'synthetic: genes, position and labels generated for the same simulated cell',
        'per_cell_join': 'exact_cell_id', 'endogenous_genes': genes,
        'position_columns': ['x_um', 'y_um', 'z_um'], 'position_measured': True,
        'technical_columns': ['log_total_counts'], 'targets': ['T1', 'T2'],
        'power_approval': {'approved_by': 'synthetic software test only (not a biological power approval)',
                           'reference': 'docs/paired_projection_reference_protocol.json', 'min_groups': 6},
        'label_semantics': 'captured_label'}
    return manifest, {'cells': cells, 'genes': gtab, 'labels': labels, 'mask': mask}


def write_dataset(manifest, tables, folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    m = json.loads(json.dumps(manifest))
    m['files'] = {}
    for k in FILE_KEYS:
        p = folder / f'{k}.csv'
        tables[k].to_csv(p, index=False)
        m['files'][k] = {'path': p.name, 'sha256': sha256_file(p)}
    mp = folder / 'manifest.json'
    mp.write_text(json.dumps(m, indent=2, sort_keys=True))
    return mp


def refusal_cases(manifest, tables, folder):
    """Each planted defect must raise Refusal. Returns {case: message}."""
    import pandas as pd
    folder = Path(folder)
    out = {}

    def expect(name, man, tabs=None):
        mp = write_dataset(man, tabs or tables, folder / name)
        try:
            run_manifest(mp)
        except Refusal as e:
            out[name] = str(e)
            return
        raise AssertionError(f'planted defect not refused: {name}')

    def mod(**kw):
        m = json.loads(json.dumps(manifest))
        m.update(kw)
        return m

    def groups(fn):
        g = json.loads(json.dumps(manifest['groups']))
        fn(g)
        return mod(groups=g)

    m = json.loads(json.dumps(manifest)); m.pop('power_approval'); expect('no_power_approval', m)
    expect('no_endpoint_mode', mod(endpoint_mode='auto'))
    expect('slices_as_animals', mod(group_unit='slice'))
    expect('hemisphere_group', groups(lambda g: g[0].update(kind='hemisphere')))
    expect('missing_animal_identity', groups(lambda g: g[0].update(animals=[])))
    expect('two_slices_of_one_animal_as_two_groups', groups(lambda g: g[1].update(animals=list(g[0]['animals']))))
    pools = [{'group_id': 'P1', 'kind': 'pool', 'animals': ['m1', 'm2'], 'independence_evidence': 'x'},
             {'group_id': 'P2', 'kind': 'pool', 'animals': ['m2', 'm3'], 'independence_evidence': 'x'}]
    expect('overlapping_pool_members', mod(group_unit='pool', groups=pools))
    expect('no_independence_evidence', mod(independence_evidence=''))
    expect('inexact_join_declared', mod(per_cell_join='nearest_cell'))
    expect('incomplete_mask', manifest, dict(tables, mask=tables['mask'].iloc[1:]))
    expect('duplicate_cell_join', manifest, dict(tables, cells=pd.concat([tables['cells'], tables['cells'].iloc[:1]])))
    expect('unmatched_cell_join', manifest, dict(tables, genes=tables['genes'].iloc[:-1]))
    gg = tables['genes'].copy(); gg['tdTomato'] = 1
    expect('reporter_gene', mod(endogenous_genes=manifest['endogenous_genes'] + ['tdTomato']), dict(tables, genes=gg))
    cc = tables['cells'].copy(); cc['cluster_id'] = 1
    expect('class_leakage_column', mod(technical_columns=['log_total_counts', 'cluster_id']), dict(tables, cells=cc))
    lab = tables['labels']
    unass = tables['cells'].loc[tables['cells'].group_id == 'A00', 'cell_id'].iloc[0]
    expect('label_in_unassayed', manifest,
           dict(tables, labels=pd.concat([lab, pd.DataFrame([[unass, 'T2', 0, 1]], columns=lab.columns)])))
    m = mod(); m['power_approval'] = dict(m['power_approval'], min_groups=50)
    expect('underpowered_groups', m)
    expect('non_synthetic_manifest', mod(synthetic=False))
    expect('reporter_as_position', mod(position_columns=['EGFP', 'y_um', 'z_um']),
           dict(tables, cells=tables['cells'].rename(columns={'x_um': 'EGFP'})))
    mk = tables['mask'].copy()
    keep2 = mk['group_id'].isin(['A02', 'A03'])
    mk.loc[(mk.target == 'T2') & ~keep2, 'assayed'] = 0
    lab3 = tables['labels'].merge(tables['cells'][['cell_id', 'group_id']], on='cell_id')
    lab3 = lab3[(lab3.target == 'T1') | lab3.group_id.isin(['A02', 'A03'])].drop(columns='group_id')
    expect('endpoint_assayed_in_two_groups', manifest, dict(tables, mask=mk, labels=lab3))
    lab2 = tables['labels'].copy(); lab2.loc[lab2.target == 'T1', 'captured'] = 0
    expect('one_class_training_fails_no_fallback', mod(targets=['T1']),
           dict(tables, labels=lab2[lab2.target == 'T1'],
                mask=tables['mask'][tables['mask'].target == 'T1']))
    try:
        assert_disjoint({'A01', 'A02'}, {'A02'})
        raise AssertionError('group leakage not refused')
    except Refusal as e:
        out['group_leakage'] = str(e)
    return out


def hand_parity(res, d, ep='T1'):
    """Recompute candidate Brier/log loss with plain loops (unequal groups) and compare."""
    p = res['per_endpoint'][ep]['_preds']['candidate_logistic']
    y = d['y'][ep]
    keep = y >= 0
    sums = {}
    for pi, yi, gi in zip(p[keep], y[keep], d['group'][keep]):
        q = min(max(pi, CLIP[0]), CLIP[1])
        s = sums.setdefault(gi, [0.0, 0.0, 0])
        s[0] += (pi - yi) ** 2
        s[1] += -(yi * np.log(q) + (1 - yi) * np.log(1 - q))
        s[2] += 1
    br = sum(v[0] / v[2] for v in sums.values()) / len(sums)
    ll = sum(v[1] / v[2] for v in sums.values()) / len(sums)
    m = res['per_endpoint'][ep]['methods']['candidate_logistic']
    return {'group_sizes': {g: v[2] for g, v in sorted(sums.items())},
            'brier_code': m['brier_group_weighted'], 'brier_hand': br,
            'log_loss_code': m['log_loss_group_weighted'], 'log_loss_hand': ll,
            'max_abs_diff': max(abs(br - m['brier_group_weighted']), abs(ll - m['log_loss_group_weighted']))}


def preprocessor_invariance(d, ep='T1'):
    """Refit fold 0 after changing hidden evaluation-group features and labels; fitted objects must be identical."""
    folds = group_folds(d['group'])
    y = d['y'][ep]
    elig = y >= 0
    tr = elig & ~np.isin(d['group'], list(folds[0]))
    ev_rows = np.isin(d['group'], list(folds[0]))
    f1 = fit_fold(d, tr, y[tr])
    d2 = dict(d, genes=d['genes'].copy(), xyz=d['xyz'].copy(), tech=d['tech'].copy(), y=dict(d['y']))
    rng = np.random.default_rng(1)
    d2['genes'][ev_rows] = rng.normal(50, 10, d2['genes'][ev_rows].shape)
    d2['xyz'][ev_rows] += 1e4
    d2['tech'][ev_rows] *= -3
    y2 = y.copy(); y2[ev_rows & elig] = 1 - y2[ev_rows & elig]; d2['y'][ep] = y2
    f2 = fit_fold(d2, tr, y2[tr])
    same = all(np.array_equal(getattr(f1[k], a), getattr(f2[k], a))
               for k in ('sc_full', 'sc_tech', 'sc_xyz', 'sc_gene') for a in ('mean_', 'scale_'))
    same &= np.array_equal(f1['candidate'].coef_, f2['candidate'].coef_)
    same &= np.array_equal(f1['technical'].coef_, f2['technical'].coef_)
    same &= f1['frequency'] == f2['frequency']
    same &= np.array_equal(f1['nb'].theta_, f2['nb'].theta_)
    return bool(same)


def synthetic_checks(res, d, refusals, parity, invariant):
    mc = res['macro']
    c = mc['candidate_logistic']
    sizes = parity['group_sizes'].values()
    return {
        'candidate_below_constant': bool(c['brier'] < mc['constant_0.5']['brier'] and c['log_loss'] < mc['constant_0.5']['log_loss']),
        'candidate_below_frequency': bool(c['brier'] < mc['fit_frequency']['brier'] and c['log_loss'] < mc['fit_frequency']['log_loss']),
        'candidate_below_location_knn': bool(c['brier'] < mc['location_knn']['brier']),
        'candidate_below_molecular_knn': bool(c['brier'] < mc['molecular_knn']['brier']),
        'candidate_below_technical': bool(c['brier'] < mc['technical']['brier']),
        'unequal_group_sizes': len(set(sizes)) > 1,
        'metric_parity_1e-12_unequal_groups': bool(parity['max_abs_diff'] < 1e-12),
        'fitted_preprocessors_invariant_to_hidden_eval_values': invariant,
        'unassayed_and_failed_qc_unknown': bool(d['status']['T2']['unknown_unassayed'] == 570 and
                                                d['status']['T1']['unknown_failed_qc'] > 0),
        'identical_eligible_ids_all_methods': True,   # enforced: evaluate() raises on any missing probability
        'all_planted_refusals_raised': len(refusals) == 21,
    }


def strip(res):
    r = {k: v for k, v in res.items() if k != 'per_endpoint'}
    r['per_endpoint'] = {e: {k: v for k, v in pt.items() if k != '_preds'} for e, pt in res['per_endpoint'].items()}
    return json.loads(json.dumps(r))


def synthetic():
    syn = OUT / 'synthetic'
    manifest, tables = make_synthetic()
    mp = write_dataset(manifest, tables, syn / 'clean')
    m, d, res = run_manifest(mp)
    refusals = refusal_cases(manifest, tables, syn / 'refusals')
    parity = hand_parity(res, d)
    invariant = preprocessor_invariance(d)
    checks = synthetic_checks(res, d, refusals, parity, invariant)
    receipt = {
        'schema': 'synarch_paired_projection_reference_synthetic_receipt_v3',
        'meaning': 'SYNTHETIC software test only; says nothing about biology. Biological gate NOT_EVALUATED.',
        'software_status': 'PASS' if all(checks.values()) else 'FAIL',
        'checks': checks, 'parity': parity, 'refusals': refusals, 'result': strip(res),
        'code_sha256': sha256_file(__file__),
        'protocol_sha256': sha256_file(ROOT / 'docs/paired_projection_reference_protocol.json'),
        'manifest_sha256': sha256_file(mp), 'seed': SEED,
        'finished_utc': datetime.now(timezone.utc).isoformat(),
    }
    p = OUT / 'synthetic_receipt.json'
    p.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False))
    print('software', receipt['software_status'], checks)
    print('macro', json.dumps(res['macro']))
    print('wrote', p, sha256_file(p))
    return receipt


if __name__ == '__main__':
    if sys.argv[1:] == ['synthetic']:
        synthetic()
    else:
        sys.exit(__doc__)
