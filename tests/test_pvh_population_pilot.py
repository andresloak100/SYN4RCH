import csv
import hashlib
import json
import subprocess

import numpy as np
import pytest

from analysis import paired_projection_reference as REF
from analysis import pvh_population_pilot as P

TOY_POOL_N = {p: v[2] for p, v in P.POOLS.items()}


# ---------------- toy inputs (no real RNA values are read anywhere in these tests)

def toy_records(extra_unjoined=5, reused=4):
    recs = []
    for pool, (target, sex, n) in P.POOLS.items():
        for k in range(n):
            recs.append({'cell_id': f'{pool[7:10]}_{k:04d}', 'encoded_pool_id': pool, 'exact_first_MUX_barcode_join': True,
                         'original_study': True, 'observed_positive_target': target, 'sex_from_encoded_pool': sex,
                         'individual_mouse_id': None, 'object_batch': '11' if target == 'PB' else 'R1_SCProj',
                         'source_study': 'Li et al. 2025/2026', 'other_projection_targets': 'unknown; not verified negatives'})
    for k in range(extra_unjoined):
        recs.append(dict(recs[0], cell_id=f'UNJ_{k}', encoded_pool_id=None, exact_first_MUX_barcode_join=False,
                         sex_from_encoded_pool=None, object_batch='HTO_Opt'))
    for k in range(reused):
        recs.append(dict(recs[-1], cell_id=f'REUSED_{k}', observed_positive_target='SC', original_study=False,
                         source_study='Winter et al. 2023', object_batch='Hyso1'))
    return recs


def toy_meta(records):
    meta = {'PB': {}, 'SC': {}}
    for i, r in enumerate(records):
        meta[r['observed_positive_target']][r['cell_id']] = {'cell_id': r['cell_id'], 'predicted.celltype': f'Seq_S{i % 3}'}
    return meta


def toy_counts(cohort, n_genes=300, seed=0, signal=True):
    rng = np.random.default_rng(seed)
    genes = [f'G{j:04d}' for j in range(n_genes)]
    lam = rng.uniform(0.5, 5, n_genes)
    rows = []
    for r in cohort:
        l = lam.copy()
        if signal and r['target'] == 'PB':
            l[:10] *= 4.0                                   # planted PB-associated genes
        rows.append(rng.poisson(l * rng.lognormal(0, 0.2)))
    return np.asarray(rows, float), genes


def write_toy_extraction(folder, cohort, counts, genes, extra_gene_sc='EGFP', extra_cell=None):
    for obj in ('PB', 'SC'):
        d = folder / obj
        d.mkdir(parents=True)
        g = genes + ([extra_gene_sc] if obj == 'SC' else ['PBonly_gene'])
        idx = [i for i, r in enumerate(cohort) if r['object'] == obj]
        cells = [cohort[i]['cell_id'] for i in idx] + ([extra_cell] if extra_cell and obj == 'PB' else [])
        (d / 'genes.txt').write_text('\n'.join(g) + '\n')
        with (d / 'cells.tsv').open('w') as f:
            f.write('cell_index\tcell_id\tlayer_column\n')
            for k, c in enumerate(cells, 1):
                f.write(f'{k}\t{c}\t{k}\n')
        with (d / 'triplets.tsv').open('w') as f:
            f.write('gene_index\tcell_index\tcount\n')
            for k, i in enumerate(idx, 1):
                for j in np.flatnonzero(counts[i]):
                    f.write(f'{j + 1}\t{k}\t{int(counts[i, j])}\n')
                f.write(f'{len(g)}\t{k}\t3\n')               # the object-specific extra gene


@pytest.fixture(scope='module')
def cohort():
    return P.build_cohort(toy_records())


# ---------------- cohort and atlas

def test_cohort_exact_pools_and_keys(cohort):
    assert len(cohort) == 265 and {r['pool'] for r in cohort} == set(P.POOLS)
    assert all(r['target'] == P.POOLS[r['pool']][0] and r['sex'] == P.POOLS[r['pool']][1] for r in cohort)


@pytest.mark.parametrize('bug', ['count', 'sex', 'target', 'duplicate', 'mouse_id'])
def test_cohort_refusals(bug):
    recs = toy_records()
    if bug == 'count':
        recs = recs[1:]
    elif bug == 'sex':
        recs[0]['sex_from_encoded_pool'] = 'male'
    elif bug == 'target':
        recs[0]['observed_positive_target'] = 'SC'
    elif bug == 'duplicate':
        recs.append(dict(recs[0]))
    elif bug == 'mouse_id':
        recs[0]['individual_mouse_id'] = 'm1'
    with pytest.raises(REF.Refusal):
        P.build_cohort(recs)


def test_atlas_original_only_unknowns_null_no_coordinates():
    recs = toy_records()
    meta = toy_meta(recs)
    nuclei, counts, graph = P.build_atlas(recs, meta['PB'], meta['SC'])
    assert len(nuclei) == 265 + 5 and not any(n['source_study'] != 'Li et al. 2025/2026' for n in nuclei)
    assert all(n['individual_animal_id'] is None and n['measured_position'] is None and n['other_target_status'] is None
               for n in nuclei)
    tt = [c for c in counts if c['grouping'] == 'author_inferred_type']
    assert sum(c['n_nuclei'] for c in tt) == len(nuclei)
    assert graph['unknown_pool_nuclei'] == 5
    for n in graph['nodes']:
        assert n.get('coordinates') is None and 'xyz' not in n and 'position' not in n
    for e in graph['edges']:
        if e['kind'] == 'collected_in':
            assert 'not prevalence' in e['meaning'] and e['other_target_status'] is None
    types = [n for n in graph['nodes'] if n['kind'] == 'author_inferred_type']
    assert types and all(n['inferred'] for n in types)


def test_atlas_refuses_missing_metadata():
    recs = toy_records()
    meta = toy_meta(recs)
    meta['PB'].pop(recs[0]['cell_id'])
    with pytest.raises(REF.Refusal, match='no metadata row'):
        P.build_atlas(recs, meta['PB'], meta['SC'])


# ---------------- extraction outputs, panel and features

def test_read_extraction_common_panel_and_reporter_exclusion(cohort, tmp_path):
    counts, genes = toy_counts(cohort)
    write_toy_extraction(tmp_path, cohort, counts, genes + ['EGFP'], extra_gene_sc='SConly_gene')
    mat, panel, excluded = P.read_extraction(tmp_path, cohort)
    assert excluded == ['EGFP'] and panel == sorted(genes) and mat.shape == (265, len(genes))
    assert np.array_equal(mat, counts)


def test_read_extraction_refuses_extra_cell(cohort, tmp_path):
    counts, genes = toy_counts(cohort)
    write_toy_extraction(tmp_path, cohort, counts, genes, extra_cell='NOT_IN_COHORT')
    with pytest.raises(REF.Refusal, match='outside the cohort'):
        P.read_extraction(tmp_path, cohort)


def test_gene_selection_train_only_and_lexical_ties():
    panel = ['b', 'a', 'c', 'd']
    x = np.array([[0, 0, 1, 5], [2, 2, 1, 5], [0, 0, 1, 5]], float)
    assert P.select_genes(x, panel, n_top=1) == [1]                 # a and b tie: 'a' first
    assert P.select_genes(x, panel, n_top=10) == [1, 0]             # constant genes excluded


def test_feature_state_ignores_evaluation_cells(cohort):
    counts, genes = toy_counts(cohort)
    pools = np.array([r['pool'] for r in cohort])
    tr = np.isin(pools, P.FOLDS['fold_F_to_M'][0])
    s1 = P.fit_features(counts[tr], genes)
    c2 = counts.copy(); c2[~tr] *= 7; c2[~tr, :5] += 100
    s2 = P.fit_features(c2[tr], genes)
    assert s1 == s2 and len(s1['selected_genes']) == min(500, len(genes))


# ---------------- methods, folds and leakage

def test_known_signal_is_recovered(cohort):
    counts, genes = toy_counts(cohort, signal=True)
    _, _, met = P.run_folds(counts, cohort, genes)
    prior = met['per_method']['train_prior']['log_loss_macro']
    assert met['per_method']['logistic_frozen_recipe']['log_loss_macro'] < prior - 0.1
    assert met['per_method']['hgbt_fixed']['log_loss_macro'] < prior
    assert met['location_comparator'].startswith('NONE AVAILABLE')


def test_no_signal_no_material_gain(cohort):
    counts, genes = toy_counts(cohort, signal=False, seed=3)
    _, _, met = P.run_folds(counts, cohort, genes)
    prior = met['per_method']['train_prior']['log_loss_macro']
    assert met['per_method']['logistic_frozen_recipe']['log_loss_macro'] > prior - 0.05


def test_test_label_swap_leaves_predictions_unchanged(cohort):
    counts, genes = toy_counts(cohort)
    p1, s1, _ = P.run_folds(counts, cohort, genes)
    for fold, (train_pools, test_pools) in P.FOLDS.items():
        swapped = [dict(r, target={'PB': 'SC', 'SC': 'PB'}[r['target']]) if r['pool'] in test_pools else r for r in cohort]
        p2, s2, _ = P.run_folds(counts, swapped, genes)
        for m in P.METHODS:
            assert np.array_equal(p1[fold]['shrunk'][m], p2[fold]['shrunk'][m])
            assert np.array_equal(p1[fold]['raw'][m], p2[fold]['raw'][m])
        assert s1[fold]['features'] == s2[fold]['features'] and s1[fold]['state']['q'] == s2[fold]['state']['q']


def test_pool_overlap_refused(cohort, monkeypatch):
    counts, genes = toy_counts(cohort)
    monkeypatch.setattr(P, 'FOLDS', {'bad': (('Mm_PVH_11F_PBprojecting', 'Mm_PVH_13F_SCprojecting'),
                                             ('Mm_PVH_11F_PBprojecting', 'Mm_PVH_14M_SCprojecting'))})
    with pytest.raises(REF.Refusal, match='pool overlap'):
        P.run_folds(counts, cohort, genes)


def test_folds_are_cross_sex_pool_transfers():
    for train, test in P.FOLDS.values():
        assert {P.POOLS[p][1] for p in train} != {P.POOLS[p][1] for p in test}
        assert {P.POOLS[p][0] for p in train} == {P.POOLS[p][0] for p in test} == {'PB', 'SC'}


def test_predict_is_features_only_and_valid(cohort):
    import inspect
    counts, genes = toy_counts(cohort)
    pools = np.array([r['pool'] for r in cohort]); y = np.array([r['target'] == 'PB' for r in cohort], int)
    tr = np.isin(pools, P.FOLDS['fold_F_to_M'][0])
    fs = P.fit_features(counts[tr], genes)
    st = P.fit_methods(P.apply_features(fs, counts[tr]), P.technical(counts[tr]), y[tr], pools[tr])
    raw, shr = P.predict_methods(st, P.apply_features(fs, counts[~tr]), P.technical(counts[~tr]))
    assert set(shr) == set(P.METHODS)
    assert not any(n in ('y', 'labels', 'target') for n in inspect.signature(P.predict_methods).parameters)
    assert np.allclose(shr['gene_copy_1nn'], 0.5 * raw['gene_copy_1nn'] + 0.5 * st['q'])
    assert st['q'] == pytest.approx(0.5)                       # one PB and one SC pool, pool-balanced


def test_metrics_cells_pool_fold_macro_by_hand():
    y = np.array([1, 0, 1, 1, 0, 0])
    pools = np.array(['a', 'a', 'b', 'b', 'c', 'd'])
    prob = np.array([0.8, 0.3, 0.6, 0.9, 0.2, 0.4])
    preds = {'f1': {'index': np.array([0, 1, 2, 3]), 'shrunk': {m: prob[[0, 1, 2, 3]] for m in P.METHODS}},
             'f2': {'index': np.array([4, 5]), 'shrunk': {m: prob[[4, 5]] for m in P.METHODS}}}
    met = P.metrics(preds, y, pools)['per_method']['train_prior']
    a = np.mean([(0.8 - 1) ** 2, 0.3 ** 2]); b = np.mean([(0.6 - 1) ** 2, (0.9 - 1) ** 2])
    c, d = 0.2 ** 2, 0.4 ** 2
    assert met['brier_macro'] == pytest.approx(np.mean([np.mean([a, b]), np.mean([c, d])]), abs=1e-15)


def test_non_feature_metadata_cannot_change_predictions(cohort):
    """Sex labels, cell IDs and author-type-like fields are not features: changing them leaves predictions unchanged."""
    import inspect
    counts, genes = toy_counts(cohort)
    p1, _, _ = P.run_folds(counts, cohort, genes)
    relabelled = [dict(r, sex='unknown', cell_id=f'X{i}', author_type='Seq_S99') for i, r in enumerate(cohort)]
    p2, _, _ = P.run_folds(counts, relabelled, genes)
    for fold in P.FOLDS:
        for m in P.METHODS:
            assert np.array_equal(p1[fold]['shrunk'][m], p2[fold]['shrunk'][m])
    assert list(inspect.signature(P.fit_features).parameters) == ['train_counts', 'panel']
    assert list(inspect.signature(P.technical).parameters) == ['counts']


# ---------------- R extractor (toy RDS mimicking the Assay5 layout; base R only)

R_TOY = r'''
genes <- c("Gm1", "Gm2", "Gm3", "Xist")
cells <- c("c1", "c2", "c3")
x <- c(5, 1, 2, 7, 3); i <- c(0L, 2L, 1L, 3L, 0L); p <- c(0L, 2L, 4L, 5L)
counts <- structure(list(), i = i, p = p, x = x * MULT, Dim = c(4L, 3L))
data <- structure(list(), i = i, p = p, x = log1p(x), Dim = c(4L, 3L))
fm <- matrix(TRUE, 4, 2, dimnames = list(genes, c("counts", "data")))
cm <- matrix(TRUE, 3, 2, dimnames = list(CELLORDER, c("counts", "data")))
rna <- structure(list(), layers = list(counts = counts, data = data), features = fm, cells = cm)
obj <- structure(list(), meta.data = data.frame(a = 1:3, row.names = cells), assays = list(RNA = rna))
saveRDS(obj, OUTPATH)
'''


def make_toy_rds(path, mult='1', cell_order='cells'):
    code = R_TOY.replace('MULT', mult).replace('CELLORDER', cell_order).replace('OUTPATH', f'"{path}"')
    subprocess.run([P.RSCRIPT, '--vanilla', '-e', code], check=True, capture_output=True)
    return hashlib.md5(path.read_bytes()).hexdigest()


def run_extractor(rds, md5, cohort_path, out):
    return subprocess.run([P.RSCRIPT, '--vanilla', str(P.EXTRACTOR), str(rds), 'PB', md5, str(cohort_path), str(out)],
                          capture_output=True, text=True)


@pytest.fixture
def cohort_tsv(tmp_path):
    p = tmp_path / 'cohort.tsv'
    p.write_text('object\tcell_id\nPB\tc3\nPB\tc1\nSC\tother\n')
    return p


def test_extractor_exports_exact_cohort_counts(tmp_path, cohort_tsv):
    rds = tmp_path / 'toy.rds'; md5 = make_toy_rds(rds)
    res = run_extractor(rds, md5, cohort_tsv, tmp_path / 'out')
    assert res.returncode == 0, res.stderr
    genes = (tmp_path / 'out/genes.txt').read_text().split()
    cells = list(csv.DictReader((tmp_path / 'out/cells.tsv').open(), delimiter='\t'))
    trip = [tuple(r.values()) for r in csv.DictReader((tmp_path / 'out/triplets.tsv').open(), delimiter='\t')]
    assert genes == ['Gm1', 'Gm2', 'Gm3', 'Xist'] and [c['cell_id'] for c in cells] == ['c3', 'c1']
    assert sorted(trip) == sorted([('1', '1', '3'), ('1', '2', '5'), ('3', '2', '1')])
    rec = dict(l.split('\t') for l in (tmp_path / 'out/receipt.tsv').read_text().splitlines())
    assert rec['layer'] == 'counts' and rec['packages_loaded_by_script'] == 'none' and rec['exported_cells'] == '2'


def test_extractor_ignores_fractional_values_outside_selected_cohort(tmp_path, cohort_tsv):
    # c2 is not requested: its fractional abundance must not invalidate c1/c3 count exports.
    code = R_TOY.replace('x <- c(5, 1, 2, 7, 3)', 'x <- c(5, 1, 2.25, 7.5, 3)')
    path = tmp_path / 'outside.rds'
    code = code.replace('MULT', '1').replace('CELLORDER', 'cells').replace('OUTPATH', f'"{path}"')
    subprocess.run([P.RSCRIPT, '--vanilla', '-e', code], check=True, capture_output=True)
    result = run_extractor(path, hashlib.md5(path.read_bytes()).hexdigest(), cohort_tsv, tmp_path / 'out')
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader((tmp_path / 'out/triplets.tsv').open(), delimiter='\t'))
    assert all(float(row['count']).is_integer() for row in rows)


@pytest.mark.parametrize('case,msg', [('noninteger', 'whole numbers'), ('cell_order', 'cell order'),
                                      ('md5', 'md5'), ('exists', 'already exists')])
def test_extractor_refusals(tmp_path, cohort_tsv, case, msg):
    rds = tmp_path / 'toy.rds'
    md5 = make_toy_rds(rds, mult='1.5' if case == 'noninteger' else '1',
                       cell_order='rev(cells)' if case == 'cell_order' else 'cells')
    out = tmp_path / 'out'
    if case == 'exists':
        out.mkdir()
    res = run_extractor(rds, '0' * 32 if case == 'md5' else md5, cohort_tsv, out)
    assert res.returncode == 3 and msg in res.stderr
    if case != 'exists':
        assert not out.exists()


def test_extractor_uses_no_packages_or_author_code():
    code = '\n'.join(l for l in P.EXTRACTOR.read_text().splitlines() if not l.lstrip().startswith('#'))
    for bad in ('library(', 'require(', 'source(', 'install.packages', 'Seurat', 'Matrix', 'loadNamespace'):
        assert bad not in code
    assert code.count('::') == code.count('tools::md5sum') == 1


# ---------------- real runs refused while the protocol is a draft

def test_real_runs_refused_while_draft(monkeypatch, tmp_path):
    draft = json.loads(P.PROTOCOL.read_text())
    draft['status'] = 'DRAFT_FOR_ROOT_REVIEW_NO_REAL_RUNS'
    path = tmp_path / 'draft.json'
    path.write_text(json.dumps(draft))
    check = P.check_frozen
    monkeypatch.setattr(P, 'check_frozen', lambda repo=None: check(repo=repo or P.ROOT, protocol_path=path))
    monkeypatch.setattr(subprocess, 'run', lambda *a, **k: (_ for _ in ()).throw(AssertionError('R called')))
    for fn, rid in ((P.execute_atlas, P.ATLAS_RUN_ID), (P.execute_extract, P.EXTRACT_RUN_ID), (P.execute_fit, P.FIT_RUN_ID)):
        before = (P.OUT / rid).exists()
        with pytest.raises(REF.Refusal, match=P.FROZEN_STATUS):
            fn()
        assert (P.OUT / rid).exists() == before


@pytest.mark.parametrize('fn,rid', [('execute_atlas', 'pvh_population_atlas_v2'), ('execute_extract', 'x'),
                                    ('execute_fit', 'pvh_population_pilot_v2')])
def test_real_run_identity_bound(fn, rid, tmp_path):
    with pytest.raises(REF.Refusal, match='canonical output root'):
        getattr(P, fn)(run_id=rid)
    with pytest.raises(REF.Refusal, match='canonical output root'):
        getattr(P, fn)(out_root=tmp_path)


def test_protocol_matches_code_constants():
    proto = json.loads(P.PROTOCOL.read_text())
    assert proto['fit_cohort']['expected_pools'] == TOY_POOL_N
    assert {k: (tuple(v['train']), tuple(v['test'])) for k, v in proto['folds'].items() if k.startswith('fold_')} == P.FOLDS
    for rel, h in proto['imports_sha256'].items():
        assert REF.sha256_file(P.ROOT / rel) == h


def test_reporter_exclusions_keep_native_symbols():
    for g in ['ERCC-001', 'ERCC_002', 'DsRed', 'EBFP2', 'HTO', 'TotalSeq', 'Cre', 'EGFP']:
        assert P.excluded_gene_reason(g)
    for g in ['Creb1', 'Creld1', 'Gfpt1', 'Gfpt2', 'Sim1']:
        assert P.excluded_gene_reason(g) is None


def test_pool_weighted_neighbor_control_by_hand(cohort):
    counts, genes = toy_counts(cohort)
    pools = np.array([r['pool'] for r in cohort]); y = np.array([r['target'] == 'PB' for r in cohort], int)
    tr = np.isin(pools, P.FOLDS['fold_F_to_M'][0])
    fs = P.fit_features(counts[tr], genes)
    st = P.fit_methods(P.apply_features(fs, counts[tr]), P.technical(counts[tr]), y[tr], pools[tr])
    xe = P.apply_features(fs, counts[~tr]); te = P.technical(counts[~tr])
    raw, _ = P.predict_methods(st, xe, te)
    donors = st['knn31'].kneighbors(xe, return_distance=False)
    inverse_sizes = np.array([1 / np.sum(pools[tr] == p) for p in pools[tr]])
    expect = (inverse_sizes[donors] * y[tr][donors]).sum(axis=1) / inverse_sizes[donors].sum(axis=1)
    np.testing.assert_allclose(raw['pool_weighted_knn31'], expect, rtol=0, atol=1e-15)


def test_saved_inference_roundtrip_from_counts(cohort, tmp_path):
    import joblib
    counts, genes = toy_counts(cohort)
    pools = np.array([r['pool'] for r in cohort]); y = np.array([r['target'] == 'PB' for r in cohort], int)
    tr = np.isin(pools, P.FOLDS['fold_F_to_M'][0])
    fs = P.fit_features(counts[tr], genes)
    st = P.fit_methods(P.apply_features(fs, counts[tr]), P.technical(counts[tr]), y[tr], pools[tr])
    path = tmp_path / 'own.joblib'
    joblib.dump({'panel': genes, 'features': fs, 'state': st}, path)
    artifact = P.load_artifact(path, REF.sha256_file(path))
    raw, pred = P.artifact_predict(artifact, genes, counts[~tr])
    direct_raw, direct = P.predict_methods(st, P.apply_features(fs, counts[~tr]), P.technical(counts[~tr]))
    for m in P.METHODS:
        np.testing.assert_array_equal(pred[m], direct[m])
        np.testing.assert_array_equal(raw[m], direct_raw[m])
    with pytest.raises(REF.Refusal, match='receipt hash'):
        P.load_artifact(path, '0' * 64)
    with pytest.raises(REF.Refusal, match='identities'):
        P.artifact_predict(artifact, genes[::-1], counts[~tr])


@pytest.mark.parametrize('bad', [np.array([[True, False]]), np.array([['1', '2']]),
                               np.array([1, 2]), np.array([[np.nan, 1]]), np.array([[np.inf, 1]]),
                               np.array([[-1, 1]]), np.array([[.5, 1]]), np.empty((0, 2))])
def test_invalid_raw_counts_refused(bad):
    with pytest.raises(REF.Refusal):
        P.cp10k_log1p(bad)
