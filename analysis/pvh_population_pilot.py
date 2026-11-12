"""PVH population-to-target pilot on PUBLIC EXISTING data (Li, Resch et al. 2026). EXPLORATORY only.

Protocol: docs/pvh_population_pilot_protocol.json (DRAFT refuses every real run).

Three fixed runs once frozen:
  atlas    metadata-only atlas and machine-readable evidence graph of the 1225 original-study nuclei
           (counts by AUTHOR_INFERRED type x positive collection target; unknowns stay null; no coordinates)
  extract  own base-R reader exports RNA counts for the exact 265 first-MUX nuclei (no author code)
  fit      PB-vs-SC collection-label pilot with two predeclared cross-sex pool folds and fixed methods

The label is membership of a retrograde COLLECTION group, never a projection probability or a negative for any
other target. Pools are not animals. No location comparator exists because no positions were measured.

Usage (refused until the protocol is frozen and committed):
  /opt/anaconda3/bin/python -I analysis/pvh_population_pilot.py atlas
  /opt/anaconda3/bin/python -I analysis/pvh_population_pilot.py extract
  /opt/anaconda3/bin/python -I analysis/pvh_population_pilot.py fit
"""
import csv
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis import paired_projection_nonlinear as N      # noqa: E402  fixed HGBT fitter (import only)
from analysis import paired_projection_reference as REF    # noqa: E402  reviewed helpers (import only)
from analysis import paired_projection_simulation as SIM   # noqa: E402  imported by N (pinned)
from analysis import paired_projection_real_benchmark as B  # noqa: E402  git commit check (pinned)

CANONICAL_SOURCE = ROOT / 'analysis/pvh_population_pilot.py'
EXTRACTOR = ROOT / 'analysis/pvh_population_extract.R'
PROTOCOL = ROOT / 'docs/pvh_population_pilot_protocol.json'
OUT = ROOT / 'data/interim/pvh_population_pilot'
RAW = ROOT / 'data/raw/paired_projection_feasibility/pvh_li2026'
META = ROOT / 'data/interim/paired_projection_feasibility/pvh_li2026'
ATLAS_RUN_ID, EXTRACT_RUN_ID, FIT_RUN_ID = 'pvh_population_atlas_v1', 'pvh_rna_extract_v1', 'pvh_population_pilot_v1'
FROZEN_STATUS = 'FROZEN_PVH_PILOT_AUTHORIZED'
RSCRIPT = '/usr/local/bin/Rscript'
POOLS = {'Mm_PVH_11F_PBprojecting': ('PB', 'female', 144), 'Mm_PVH_12M_PBprojecting': ('PB', 'male', 53),
         'Mm_PVH_13F_SCprojecting': ('SC', 'female', 39), 'Mm_PVH_14M_SCprojecting': ('SC', 'male', 29)}
FOLDS = {'fold_F_to_M': (('Mm_PVH_11F_PBprojecting', 'Mm_PVH_13F_SCprojecting'),
                         ('Mm_PVH_12M_PBprojecting', 'Mm_PVH_14M_SCprojecting')),
         'fold_M_to_F': (('Mm_PVH_12M_PBprojecting', 'Mm_PVH_14M_SCprojecting'),
                         ('Mm_PVH_11F_PBprojecting', 'Mm_PVH_13F_SCprojecting'))}
N_TOP = 500
LOGIT = dict(C=0.1, solver='lbfgs', max_iter=2000, tol=1e-6, random_state=20261112)
ALPHA = {'logistic_frozen_recipe': 0.1, 'hgbt_fixed': 0.1, 'molecular_knn31': 0.1, 'gene_copy_1nn': 0.5,
         'train_prior': 0.0, 'technical_logistic': 0.1, 'pool_weighted_knn31': 0.1}
METHODS = list(ALPHA)
CLIP = (1e-6, 1 - 1e-6)
OBJECT_MD5 = {'PB': '54b226c47a1dc52373d7a0dfa72386d5', 'SC': '40f516689feace69b6d72898edd1b60a'}
Refusal = REF.Refusal
sha = REF.sha256_file


def write_json(path, obj):
    with Path(path).open('x') as f:
        json.dump(obj, f, indent=2, sort_keys=True, allow_nan=False)
        f.write('\n')


# ---------------------------------------------------------------- cohort and atlas (metadata only)

def object_of(record):
    t = record['observed_positive_target']
    if t not in ('PB', 'SC'):
        raise Refusal(f'unexpected positive target {t!r}')
    return t


def build_cohort(records):
    """Exact 265 first-MUX original-study nuclei with consistent pool, sex and target."""
    rows, seen = [], set()
    for r in records:
        if r.get('exact_first_MUX_barcode_join') is not True or r.get('original_study') is not True:
            continue
        pool = r.get('encoded_pool_id')
        if pool not in POOLS:
            raise Refusal(f'joined nucleus with unknown pool {pool!r}')
        target, sex, _ = POOLS[pool]
        if object_of(r) != target or r.get('sex_from_encoded_pool') != sex:
            raise Refusal(f'pool/target/sex inconsistency for {r["cell_id"]}')
        if r.get('individual_mouse_id') is not None:
            raise Refusal('individual mouse IDs are not expected in this cohort')
        key = (target, r['cell_id'])
        if key in seen:
            raise Refusal(f'duplicate cohort key {key}')
        seen.add(key)
        rows.append({'object': target, 'cell_id': r['cell_id'], 'pool': pool, 'sex': sex, 'target': target})
    counts = Counter(r['pool'] for r in rows)
    if dict(counts) != {p: v[2] for p, v in POOLS.items()}:
        raise Refusal(f'cohort pool counts {dict(counts)} differ from the protocol')
    rows.sort(key=lambda r: (r['object'], r['cell_id']))
    return rows


def cohort_sha256(rows):
    import hashlib
    text = '\n'.join(f"{r['object']}\t{r['cell_id']}\t{r['pool']}\t{r['sex']}\t{r['target']}" for r in rows)
    return hashlib.sha256(text.encode()).hexdigest()


def read_metadata_csv(path):
    with Path(path).open(newline='') as f:
        rows = list(csv.DictReader(f))
    ids = [r['cell_id'] for r in rows]
    if len(set(ids)) != len(ids):
        raise Refusal(f'duplicate cell IDs in {path}')
    return {r['cell_id']: r for r in rows}


def build_atlas(records, pb_meta, sc_meta):
    """1225 original-study nuclei: records, counts and an evidence graph. Unknowns stay None."""
    meta = {'PB': pb_meta, 'SC': sc_meta}
    nuclei = []
    for r in records:
        if r.get('original_study') is not True:
            continue
        obj = object_of(r)
        m = meta[obj].get(r['cell_id'])
        if m is None:
            raise Refusal(f'original nucleus {obj}|{r["cell_id"]} has no metadata row')
        pool = r.get('encoded_pool_id')
        nuclei.append({'nucleus_key': f"{obj}|{r['cell_id']}", 'object': obj, 'cell_id': r['cell_id'],
                       'source_study': r['source_study'], 'object_batch': r['object_batch'],
                       'encoded_pool_id': pool, 'sex_from_encoded_pool': r.get('sex_from_encoded_pool'),
                       'positive_collection_target': obj,
                       'author_inferred_type': m.get('predicted.celltype') or None,
                       'individual_animal_id': None, 'measured_position': None,
                       'other_target_status': None})
    keys = [n['nucleus_key'] for n in nuclei]
    if len(set(keys)) != len(keys):
        raise Refusal('duplicate atlas nucleus keys')
    by_target = Counter(n['positive_collection_target'] for n in nuclei)
    type_target = Counter((n['author_inferred_type'], n['positive_collection_target']) for n in nuclei)
    pool_target = Counter((n['encoded_pool_id'], n['positive_collection_target']) for n in nuclei)
    batch_target = Counter((n['object_batch'], n['positive_collection_target']) for n in nuclei)
    counts = [{'grouping': 'author_inferred_type', 'group': g, 'positive_collection_target': t, 'n_nuclei': c}
              for (g, t), c in sorted(type_target.items(), key=lambda kv: (str(kv[0][0]), kv[0][1]))]
    counts += [{'grouping': 'encoded_pool_id', 'group': g, 'positive_collection_target': t, 'n_nuclei': c}
               for (g, t), c in sorted(pool_target.items(), key=lambda kv: (str(kv[0][0]), kv[0][1]))]
    counts += [{'grouping': 'object_batch', 'group': g, 'positive_collection_target': t, 'n_nuclei': c}
               for (g, t), c in sorted(batch_target.items(), key=lambda kv: (str(kv[0][0]), kv[0][1]))]
    graph = evidence_graph(nuclei, type_target, pool_target, by_target)
    return nuclei, counts, graph


def evidence_graph(nuclei, type_target, pool_target, by_target):
    nodes = [{'id': 'region:PVH', 'kind': 'source_region', 'label': 'paraventricular hypothalamus (PVH)',
              'coordinates': None},
             {'id': 'target:PB', 'kind': 'collection_target', 'label': 'parabrachial retrograde collection group',
              'n_nuclei_collected': by_target.get('PB', 0)},
             {'id': 'target:SC', 'kind': 'collection_target', 'label': 'spinal cord retrograde collection group',
              'n_nuclei_collected': by_target.get('SC', 0)}]
    for t in sorted({k[0] for k in type_target}, key=str):
        nodes.append({'id': f'type:{t}', 'kind': 'author_inferred_type', 'label': t, 'inferred': True,
                      'note': 'author predicted.celltype; inferred identity, not a measurement'})
    for p in sorted({k[0] for k in pool_target if k[0] is not None}):
        nodes.append({'id': f'pool:{p}', 'kind': 'encoded_pool', 'label': p, 'sex': POOLS[p][1],
                      'is_individual_animal': False, 'note': 'pool of 4-10 mice; not an animal'})
    edges = []
    for (t, tgt), c in sorted(type_target.items(), key=lambda kv: (str(kv[0][0]), kv[0][1])):
        edges.append({'source': f'type:{t}', 'target': f'target:{tgt}', 'kind': 'collected_in',
                      'n_nuclei': c, 'n_target_collected': by_target[tgt],
                      'fraction_of_target_collection': c / by_target[tgt],
                      'meaning': 'sampling count among collected nuclei; not prevalence or projection probability',
                      'other_target_status': None})
    for (p, tgt), c in sorted(pool_target.items(), key=lambda kv: (str(kv[0][0]), kv[0][1])):
        if p is not None:
            edges.append({'source': f'pool:{p}', 'target': f'target:{tgt}', 'kind': 'pool_collected_for',
                          'n_nuclei': c})
    edges += [{'source': 'region:PVH', 'target': f'target:{t}', 'kind': 'retrograde_collection_from',
               'n_nuclei': by_target[t]} for t in ('PB', 'SC')]
    return {'schema': 'synarch_pvh_evidence_graph_v1', 'nodes': nodes, 'edges': edges,
            'unknown_pool_nuclei': sum(1 for n in nuclei if n['encoded_pool_id'] is None),
            'scope': 'original-study nuclei only (reused SC nuclei excluded); counts are sampling counts',
            'not_present': ['coordinates or positions', 'individual animals', 'other-target negatives',
                            'dual projection', 'axon paths or synapses']}


# ---------------------------------------------------------------- extraction outputs and features

def read_extraction(folder, cohort):
    """Read per-object R exports; return the common panel count matrix aligned to the cohort order."""
    folder = Path(folder)
    per = {}
    for obj in ('PB', 'SC'):
        d = folder / obj
        genes = (d / 'genes.txt').read_text().splitlines()
        if len(set(genes)) != len(genes):
            raise Refusal(f'duplicate gene IDs in {obj}')
        with (d / 'cells.tsv').open() as f:
            cells = [r['cell_id'] for r in csv.DictReader(f, delimiter='\t')]
        trip = np.loadtxt(d / 'triplets.tsv', delimiter='\t', skiprows=1, ndmin=2)
        mat = np.zeros((len(cells), len(genes)))
        if trip.size:
            gi, ci, cv = trip[:, 0].astype(int) - 1, trip[:, 1].astype(int) - 1, trip[:, 2]
            if (cv < 0).any() or (cv != np.round(cv)).any():
                raise Refusal(f'{obj} exported counts are not non-negative whole numbers')
            np.add.at(mat, (ci, gi), cv)
        per[obj] = (genes, cells, mat)
    common = sorted(set(per['PB'][0]) & set(per['SC'][0]))
    excluded = [g for g in common if excluded_gene_reason(g)]
    panel = [g for g in common if g not in set(excluded)]
    rows = []
    for r in cohort:
        genes, cells, mat = per[r['object']]
        hits = [i for i, c in enumerate(cells) if c == r['cell_id']]
        if len(hits) != 1:
            raise Refusal(f"cohort cell {r['object']}|{r['cell_id']} found {len(hits)} times in the export")
        idx = {g: j for j, g in enumerate(genes)}
        rows.append(mat[hits[0], [idx[g] for g in panel]])
    exported = {o: set(per[o][1]) for o in per}
    wanted = {o: {r['cell_id'] for r in cohort if r['object'] == o} for o in per}
    if exported != wanted:
        raise Refusal('export contains cells outside the cohort or misses cohort cells')
    return np.vstack(rows), panel, excluded


def excluded_gene_reason(gene):
    """Fixed endogenous-panel exclusions; preserve native Creb/Creld/Gfpt genes."""
    reason = REF.forbidden_name(gene)
    tokens = re.split(r'[^a-z0-9]+', str(gene).casefold())
    if reason:
        return reason
    if re.match(r'^ercc[-_]', str(gene), flags=re.I) or any(
            t in {'dsred', 'bfp', 'ebfp', 'ebfp2', 'mrfp', 'mrfp1', 'totalseq', 'hto'} for t in tokens):
        return 'reporter/vector/spike feature'
    return None


def checked_counts(counts):
    a = np.asarray(counts)
    if a.dtype == bool or not (np.issubdtype(a.dtype, np.integer) or np.issubdtype(a.dtype, np.floating)):
        raise Refusal('counts must be a numeric 2-D matrix')
    if a.ndim != 2 or min(a.shape) < 1:
        raise Refusal('counts must be a nonempty 2-D matrix')
    if not np.isfinite(a).all() or (a < 0).any() or (a != np.round(a)).any():
        raise Refusal('counts must be finite non-negative whole numbers')
    return a.astype(float, copy=False)


def cp10k_log1p(counts):
    counts = checked_counts(counts)
    tot = counts.sum(axis=1, keepdims=True)
    if (tot <= 0).any():
        raise Refusal('a nucleus has zero counts on the common panel')
    return np.log1p(counts / tot * 1e4)


def technical(counts):
    counts = checked_counts(counts)
    return np.column_stack([np.log1p(counts.sum(axis=1)), np.log1p((counts > 0).sum(axis=1))])


def select_genes(train_logx, panel, n_top=N_TOP):
    """Top min(n_top, non-constant) by TRAINING variance; ties by gene ID (lexical)."""
    var = train_logx.var(axis=0)
    nonconst = [j for j in range(len(panel)) if var[j] > 0]
    order = sorted(nonconst, key=lambda j: (-var[j], panel[j]))
    return order[:min(n_top, len(order))]


def fit_features(train_counts, panel):
    from sklearn.preprocessing import StandardScaler
    lx = cp10k_log1p(train_counts)
    if lx.shape[1] != len(panel) or len(set(panel)) != len(panel):
        raise Refusal('feature panel width or unique gene IDs differ from counts')
    sel = select_genes(lx, panel)
    if not sel:
        raise Refusal('no non-constant training genes')
    scaler = StandardScaler().fit(lx[:, sel])
    return {'selected_index': sel, 'selected_genes': [panel[j] for j in sel], 'n_panel_genes': len(panel),
            'scaler_mean': scaler.mean_.tolist(), 'scaler_scale': scaler.scale_.tolist()}


def apply_features(state, counts):
    counts = checked_counts(counts)
    if counts.shape[1] != state['n_panel_genes']:
        raise Refusal('counts width differs from fitted common panel')
    lx = cp10k_log1p(counts)[:, state['selected_index']]
    return (lx - np.asarray(state['scaler_mean'])) / np.asarray(state['scaler_scale'])


# ---------------------------------------------------------------- fixed methods (features-only predict)

def fit_methods(Xtr, Ttr, ytr, pools_tr):
    """Fit every fixed method on TRAINING rows only. Returns fitted state (no evaluation data touched)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
    from sklearn.preprocessing import StandardScaler
    ytr = np.asarray(ytr)
    if set(np.unique(ytr)) != {0, 1}:
        raise Refusal('training labels must contain both collection groups')
    w = REF.group_balanced_weights(np.asarray(pools_tr))
    q = float(np.average(ytr, weights=w))
    tsc = StandardScaler().fit(Ttr)
    hgbt, hgbt_adapter, hq = N._fit_model(Xtr, ytr, pools_tr)
    return {'q': q, 'hgbt_q': hq,
            'logistic': LogisticRegression(**LOGIT).fit(Xtr, ytr, sample_weight=w),
            'hgbt': hgbt, 'hgbt_adapter': hgbt_adapter,
            'knn31': KNeighborsClassifier(n_neighbors=31).fit(Xtr, ytr),
            'nn1': NearestNeighbors(n_neighbors=1).fit(Xtr), 'ytr': ytr.copy(), 'pool_weights': w.copy(),
            'tech_scaler': tsc,
            'technical': LogisticRegression(**LOGIT).fit(tsc.transform(Ttr), ytr, sample_weight=w)}


def predict_methods(state, Xev, Tev):
    """Features only; no labels are accepted. Returns raw and shrunk probabilities per method."""
    n = len(Xev)
    donors31 = state['knn31'].kneighbors(Xev, return_distance=False)
    donor_weights = state['pool_weights'][donors31]
    raw = {'logistic_frozen_recipe': state['logistic'].predict_proba(Xev)[:, 1],
           'hgbt_fixed': state['hgbt'].predict_proba(N._apply_adapter(state['hgbt_adapter'], Xev))[:, 1],
           'molecular_knn31': state['knn31'].predict_proba(Xev)[:, 1],
           'pool_weighted_knn31': (donor_weights * state['ytr'][donors31]).sum(axis=1) / donor_weights.sum(axis=1),
           'gene_copy_1nn': state['ytr'][state['nn1'].kneighbors(Xev, return_distance=False)[:, 0]].astype(float),
           'train_prior': np.full(n, state['q']),
           'technical_logistic': state['technical'].predict_proba(state['tech_scaler'].transform(Tev))[:, 1]}
    shrunk = {}
    for m in METHODS:
        r = REF.checked_probability(np.asarray(raw[m], float), n, f'{m} raw')
        shrunk[m] = REF.checked_probability((1 - ALPHA[m]) * r + ALPHA[m] * state['q'], n, m)
        raw[m] = r
    return raw, shrunk


def artifact_predict(artifact, gene_names, counts):
    """Reuse one fitted fold from raw measured counts; evaluation labels are not accepted."""
    if list(gene_names) != artifact['panel']:
        raise Refusal('gene identities/order differ from the fitted panel')
    return predict_methods(artifact['state'], apply_features(artifact['features'], counts), technical(counts))


def load_artifact(path, expected_sha256):
    """Load only an own, receipt-pinned local model; never third-party pickles."""
    import joblib
    if sha(path) != expected_sha256:
        raise Refusal('model artifact differs from the expected receipt hash')
    return joblib.load(path)


def run_folds(counts, cohort, panel):
    """Two predeclared cross-sex pool folds. Returns keyed predictions, fold states and metrics."""
    pools = np.array([r['pool'] for r in cohort])
    y = np.array([1 if r['target'] == 'PB' else 0 for r in cohort])
    T = technical(counts)
    preds, states = {}, {}
    for fold, (train_pools, test_pools) in FOLDS.items():
        if set(train_pools) & set(test_pools):
            raise Refusal('pool overlap between training and test')
        tr, te = np.isin(pools, train_pools), np.isin(pools, test_pools)
        fs = fit_features(counts[tr], panel)
        st = fit_methods(apply_features(fs, counts[tr]), T[tr], y[tr], pools[tr])
        raw, shr = predict_methods(st, apply_features(fs, counts[te]), T[te])
        states[fold] = {'features': fs, 'state': st, 'train_index': np.flatnonzero(tr).tolist(),
                        'test_index': np.flatnonzero(te).tolist()}
        preds[fold] = {'index': np.flatnonzero(te), 'raw': raw, 'shrunk': shr}
    return preds, states, metrics(preds, y, pools)


def metrics(preds, y, pools):
    """Cells -> pool -> fold -> macro, for log loss and Brier; AUROC per fold descriptive only."""
    from sklearn.metrics import roc_auc_score
    out = {'per_method': {}}
    for m in METHODS:
        per_fold = {}
        for fold, p in preds.items():
            idx = p['index']; prob = p['shrunk'][m]
            per_pool = {}
            for pool in sorted(set(pools[idx])):
                s = pools[idx] == pool
                yy, pp = y[idx][s], prob[s]
                pc = np.clip(pp, *CLIP)
                per_pool[pool] = {'n': int(s.sum()), 'brier': float(np.mean((pp - yy) ** 2)),
                                  'log_loss': float(np.mean(-(yy * np.log(pc) + (1 - yy) * np.log(1 - pc))))}
            per_fold[fold] = {'per_pool': per_pool,
                              'log_loss': float(np.mean([v['log_loss'] for v in per_pool.values()])),
                              'brier': float(np.mean([v['brier'] for v in per_pool.values()])),
                              'auroc_descriptive': float(roc_auc_score(y[idx], prob)) if len(set(y[idx])) == 2 else None}
        out['per_method'][m] = {'per_fold': per_fold,
                                'log_loss_macro': float(np.mean([v['log_loss'] for v in per_fold.values()])),
                                'brier_macro': float(np.mean([v['brier'] for v in per_fold.values()]))}
    out['location_comparator'] = 'NONE AVAILABLE: no measured positions in these RNA data'
    out['note'] = 'exploratory; no cell-level p-values; no promotion regardless of score'
    return out


# ---------------------------------------------------------------- pins and durable runs

def check_frozen(repo=ROOT, protocol_path=PROTOCOL, source_path=None):
    proto = json.loads(Path(protocol_path).read_text())
    if proto.get('status') != FROZEN_STATUS:
        raise Refusal(f"PVH protocol status {proto.get('status')!r} is not {FROZEN_STATUS}")
    src = Path(source_path or __file__).resolve()
    if src != CANONICAL_SOURCE.resolve():
        raise Refusal('executing source is not the canonical analysis/pvh_population_pilot.py')
    for p, what in ((protocol_path, 'protocol'), (src, 'source'), (EXTRACTOR, 'extractor')):
        if not B.protocol_committed(p, repo):
            raise Refusal(f'{what} is not committed in git or has uncommitted changes')
    for rel, h in proto['imports_sha256'].items():
        if sha(ROOT / rel) != h:
            raise Refusal(f'imported module {rel} differs from the pin')
    for name, info in proto['source']['objects'].items():
        if sha(RAW / name) != info['sha256']:
            raise Refusal(f'{name} differs from the pin')
    for name, h in proto['source']['metadata'].items():
        if sha(META / name) != h:
            raise Refusal(f'{name} differs from the pin')
    if proto['fit_cohort']['expected_pools'] != {p: v[2] for p, v in POOLS.items()}:
        raise Refusal('protocol pools differ from the code')
    if {k: (tuple(v['train']), tuple(v['test'])) for k, v in proto['folds'].items() if k.startswith('fold_')} != FOLDS:
        raise Refusal('protocol folds differ from the code')
    return proto


def _begin(run_id, out_root, extra):
    run_dir = Path(out_root) / run_id
    try:
        run_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        raise Refusal(f'{run_id} already exists; duplicate runs are refused')
    shutil.copyfile(PROTOCOL, run_dir / 'protocol_snapshot.json')
    import sklearn
    write_json(run_dir / 'STARTED.json', dict(extra, run_id=run_id, started_utc=datetime.now(timezone.utc).isoformat(),
                                              source_sha256=sha(__file__), extractor_sha256=sha(EXTRACTOR),
                                              protocol_snapshot_sha256=sha(run_dir / 'protocol_snapshot.json'),
                                              environment={'python': platform.python_version(),
                                                           'numpy': np.__version__, 'scikit_learn': sklearn.__version__}))
    return run_dir


def _finish(run_dir, extra):
    files = {}
    for p in sorted(run_dir.rglob('*')):
        if p.is_file():
            files[str(p.relative_to(run_dir))] = sha(p)
    write_json(run_dir / 'RUN_COMPLETE.json', dict(extra, files_sha256=files,
                                                    finished_utc=datetime.now(timezone.utc).isoformat()))


def _canonical(run_id, expected, out_root):
    if run_id != expected or Path(out_root).resolve() != OUT.resolve():
        raise Refusal(f'real run must use run_id {expected} and the canonical output root')


def load_records():
    return json.loads((META / 'paired_nucleus_metadata_manifest.json').read_text())


def execute_atlas(run_id=ATLAS_RUN_ID, out_root=OUT, repo=ROOT):
    _canonical(run_id, ATLAS_RUN_ID, out_root)
    check_frozen(repo)
    run_dir = _begin(run_id, out_root, {'phase': 'ATLAS_METADATA_ONLY'})
    records = load_records()
    nuclei, counts, graph = build_atlas(records, read_metadata_csv(META / 'PB_cell_metadata.csv'),
                                        read_metadata_csv(META / 'SC_cell_metadata.csv'))
    _write_csv(run_dir / 'atlas_nuclei.csv', nuclei)
    _write_csv(run_dir / 'atlas_counts.csv', counts)
    write_json(run_dir / 'evidence_graph.json', dict(graph, provenance={
        'manifest_sha256': sha(META / 'paired_nucleus_metadata_manifest.json'),
        'PB_metadata_sha256': sha(META / 'PB_cell_metadata.csv'), 'SC_metadata_sha256': sha(META / 'SC_cell_metadata.csv'),
        'run_id': run_id}))
    _finish(run_dir, {'run_id': run_id, 'n_nuclei': len(nuclei)})
    return run_dir


def _write_csv(path, rows):
    with Path(path).open('x', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def execute_extract(run_id=EXTRACT_RUN_ID, out_root=OUT, repo=ROOT):
    _canonical(run_id, EXTRACT_RUN_ID, out_root)
    check_frozen(repo)
    cohort = build_cohort(load_records())
    run_dir = _begin(run_id, out_root, {'phase': 'RNA_COUNTS_EXTRACTION', 'cohort_sha256': cohort_sha256(cohort)})
    with (run_dir / 'cohort.tsv').open('x') as f:
        f.write('object\tcell_id\n')
        for r in cohort:
            f.write(f"{r['object']}\t{r['cell_id']}\n")
    _write_csv(run_dir / 'cohort_full.csv', cohort)
    for obj in ('PB', 'SC'):
        cmd = [RSCRIPT, '--vanilla', str(EXTRACTOR), str(RAW / f'{obj}_projecting.rds'), obj, OBJECT_MD5[obj],
               str(run_dir / 'cohort.tsv'), str(run_dir / obj)]
        res = subprocess.run(cmd, capture_output=True, text=True)
        (run_dir / f'extract_{obj}.log').write_text(res.stdout + res.stderr)
        if res.returncode != 0:
            raise Refusal(f'R extraction refused for {obj}: {res.stderr.strip()[-300:]}')
    counts, panel, excluded = read_extraction(run_dir, cohort)
    write_json(run_dir / 'panel.json', {'common_genes': len(panel), 'excluded_gene_ids': excluded,
                                        'panel_sha256': __import__('hashlib').sha256('\n'.join(panel).encode()).hexdigest()})
    _finish(run_dir, {'run_id': run_id, 'cohort_sha256': cohort_sha256(cohort), 'n_cells': len(cohort)})
    return run_dir


def check_extraction(ext_dir, cohort):
    done = json.loads((Path(ext_dir) / 'RUN_COMPLETE.json').read_text())
    for name, h in done['files_sha256'].items():
        if sha(Path(ext_dir) / name) != h:
            raise Refusal(f'extraction file {name} changed after completion')
    if done['cohort_sha256'] != cohort_sha256(cohort):
        raise Refusal('extraction cohort differs from the current cohort')
    return done


def fit_and_save(run_dir, counts, cohort, panel, excluded):
    import joblib
    preds, states, met = run_folds(counts, cohort, panel)
    y = [1 if r['target'] == 'PB' else 0 for r in cohort]
    with (run_dir / 'predictions.csv').open('x', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['fold', 'object', 'cell_id', 'pool', 'sex', 'collection_target', 'method', 'raw_probability',
                    'probability'])
        for fold, p in preds.items():
            for k, i in enumerate(p['index']):
                r = cohort[i]
                for m in METHODS:
                    w.writerow([fold, r['object'], r['cell_id'], r['pool'], r['sex'], r['target'], m,
                                repr(float(p['raw'][m][k])), repr(float(p['shrunk'][m][k]))])
    for fold, s in states.items():
        st = s['state']
        write_json(run_dir / f'fold_{fold}.json', {
            'train_pools': list(FOLDS[fold][0]), 'test_pools': list(FOLDS[fold][1]),
            'train_cells': [f"{cohort[i]['object']}|{cohort[i]['cell_id']}" for i in s['train_index']],
            'selected_genes': s['features']['selected_genes'], 'scaler_mean': s['features']['scaler_mean'],
            'scaler_scale': s['features']['scaler_scale'], 'q_pool_balanced': st['q'],
            'hgbt_adapter': st['hgbt_adapter'], 'alphas': ALPHA, 'logistic': LOGIT, 'hgbt': N.HGBT})
        with (run_dir / f'models_{fold}.joblib').open('xb') as f:
            joblib.dump({'schema': 'synarch_pvh_fold_inference_v1', 'panel': list(panel),
                         'features': s['features'], 'state': st,
                         'train_cell_keys': [f"{cohort[i]['object']}|{cohort[i]['cell_id']}" for i in s['train_index']]}, f)
        path = run_dir / f'models_{fold}.joblib'
        loaded = load_artifact(path, sha(path))
        rraw, rpred = artifact_predict(loaded, panel, counts[s['test_index']])
        errors = [float(np.max(np.abs(rpred[m] - preds[fold]['shrunk'][m]))) for m in METHODS]
        errors += [float(np.max(np.abs(rraw[m] - preds[fold]['raw'][m]))) for m in METHODS]
        if max(errors) > 1e-12:
            raise Refusal('saved model features-only replay differs from its forecasts')
        write_json(run_dir / f'artifact_replay_{fold}.json', {'max_abs_difference': max(errors),
                                                            'model_sha256': sha(path), 'refits': 0})
    write_json(run_dir / 'metrics.json', dict(met, n_cells=len(cohort), n_panel_genes=len(panel),
                                              excluded_gene_ids=excluded, labels_pb=int(sum(y))))
    return met


def execute_fit(run_id=FIT_RUN_ID, out_root=OUT, repo=ROOT):
    _canonical(run_id, FIT_RUN_ID, out_root)
    check_frozen(repo)
    cohort = build_cohort(load_records())
    ext = OUT / EXTRACT_RUN_ID
    done = check_extraction(ext, cohort)
    run_dir = _begin(run_id, out_root, {'phase': 'EXPLORATORY_FIT', 'cohort_sha256': cohort_sha256(cohort),
                                        'extraction_run_complete_sha256': sha(ext / 'RUN_COMPLETE.json')})
    counts, panel, excluded = read_extraction(ext, cohort)
    met = fit_and_save(run_dir, counts, cohort, panel, excluded)
    _finish(run_dir, {'run_id': run_id, 'extraction': done['run_id']})
    return run_dir, met


if __name__ == '__main__':
    a = sys.argv[1:]
    if a == ['atlas']:
        print('atlas complete', execute_atlas())
    elif a == ['extract']:
        print('extraction complete', execute_extract())
    elif a == ['fit']:
        rd, met = execute_fit()
        print('fit complete', rd, {m: round(v['log_loss_macro'], 4) for m, v in met['per_method'].items()})
    else:
        sys.exit(__doc__)
