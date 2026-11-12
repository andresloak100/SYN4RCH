import json

import numpy as np
import pytest

from analysis import paired_projection_reference as R


@pytest.fixture(scope='module')
def small(tmp_path_factory):
    m, t = R.make_synthetic(sizes=(60, 140, 90, 120, 70, 110), n_genes=6)
    m['power_approval']['min_groups'] = 4
    folder = tmp_path_factory.mktemp('syn')
    return m, t, folder


def test_clean_synthetic_runs_and_gate_not_evaluated(small):
    m, t, folder = small
    mp = R.write_dataset(m, t, folder / 'clean')
    _, d, res = R.run_manifest(mp)
    assert res['biological_gate'].startswith('NOT_EVALUATED')
    assert set(res['macro']) == set(R.METHODS)
    ids = {e: r['eligible_cell_ids_sha256'] for e, r in res['per_endpoint'].items()}
    assert len(ids) == 2
    assert d['status']['T2']['unknown_unassayed'] == 60 + 140


def test_manifest_refused_before_any_matrix_read(small, monkeypatch):
    m, t, folder = small
    mp = R.write_dataset(dict(m, group_unit='slice'), t, folder / 'slice')
    def boom(*a, **k):
        raise AssertionError('data file parsed before manifest refusal')
    monkeypatch.setattr(R, 'read_table', boom)
    with pytest.raises(R.Refusal):
        R.run_manifest(mp)


@pytest.mark.parametrize('field', ['endpoint_mode', 'power_approval', 'independence_evidence', 'groups', 'files'])
def test_missing_required_field_refused(small, field):
    m, t, folder = small
    mp = R.write_dataset(m, t, folder / 'full')
    full = json.loads(mp.read_text())
    assert R.validate_manifest(full, mp.parent)          # the complete manifest is admitted
    bad = json.loads(json.dumps(full)); bad.pop(field)
    with pytest.raises(R.Refusal):
        R.validate_manifest(bad, mp.parent)


def test_overlapping_pool_members_and_hemispheres_refused(small):
    m, t, folder = small
    mp = R.write_dataset(m, t, folder / 'full2')
    m = json.loads(mp.read_text()); folder = mp.parent
    pools = [{'group_id': 'P1', 'kind': 'pool', 'animals': ['a', 'b'], 'independence_evidence': 'x'},
             {'group_id': 'P2', 'kind': 'pool', 'animals': ['b', 'c'], 'independence_evidence': 'x'}]
    with pytest.raises(R.Refusal, match='not disjoint'):
        R.validate_manifest(dict(m, group_unit='pool', groups=pools), folder)
    g = json.loads(json.dumps(m['groups'])); g[0]['kind'] = 'hemisphere'
    with pytest.raises(R.Refusal, match='cannot stand for an animal'):
        R.validate_manifest(dict(m, groups=g), folder)


def test_group_balanced_weights():
    w = R.group_balanced_weights(np.array(['a', 'a', 'a', 'b']))
    assert np.allclose(w.mean(), 1.0)
    assert np.isclose(w[:3].sum(), w[3])            # each group carries equal total weight


def test_hand_metric_unequal_groups():
    p = np.array([0.9, 0.2, 0.6, 0.0, 1.0])
    y = np.array([1, 0, 0, 0, 1])
    g = np.array(['a', 'a', 'a', 'b', 'b'])
    br, ll = R.group_weighted(p, y, g)
    br_a = (0.01 + 0.04 + 0.36) / 3
    br_b = 0.0
    ll_a = -(np.log(0.9) + np.log(0.8) + np.log(0.4)) / 3
    ll_b = -(np.log(1 - 1e-6) + np.log(1 - 1e-6)) / 2
    assert abs(br - (br_a + br_b) / 2) < 1e-15
    assert abs(ll - (ll_a + ll_b) / 2) < 1e-12


def test_exclusive_two_target_endpoint_statuses():
    cap = {'LHb': np.array([1, 0, 1, 0, 1, -1]), 'VTA': np.array([0, 1, 1, 0, 0, -1])}
    qc = {'LHb': np.array([1, 1, 1, 1, 0, -1]), 'VTA': np.array([1, 1, 1, 1, 1, -1])}
    assayed = {'LHb': np.array([1, 1, 1, 1, 1, 0], bool), 'VTA': np.array([1, 1, 1, 1, 1, 0], bool)}
    y, st = R.derive_endpoint('exclusive_two_target', ['LHb', 'VTA'], cap, qc, assayed)
    assert list(y['LHb_only_vs_VTA_only']) == [1, 0, -1, -1, -1, -1]
    s = st['LHb_only_vs_VTA_only']
    assert s['double_labelled_excluded'] == 1 and s['neither_excluded'] == 1
    assert s['unknown_unassayed_or_failed_qc'] == 2


def test_capture_mode_zero_only_when_assayed_and_qc_passed():
    cap = {'T': np.array([1, 0, 0, -1])}
    qc = {'T': np.array([1, 1, 0, -1])}
    assayed = {'T': np.array([1, 1, 1, 0], bool)}
    y, st = R.derive_endpoint('capture_per_target', ['T'], cap, qc, assayed)
    assert list(y['T']) == [1, 0, -1, -1]
    assert st['T']['unknown_failed_qc'] == 1 and st['T']['unknown_unassayed'] == 1


def test_group_leakage_refused():
    with pytest.raises(R.Refusal):
        R.assert_disjoint({'a', 'b'}, {'b'})


def test_preprocessor_invariance(small):
    m, t, folder = small
    mp = R.write_dataset(m, t, folder / 'inv')
    _, d, _ = R.run_manifest(mp)
    assert R.preprocessor_invariance(d)


def test_missing_probability_is_fail_not_row_drop(small, monkeypatch):
    m, t, folder = small
    mp = R.write_dataset(m, t, folder / 'nan')
    real = R.predict_fold
    def holed(d, f, ev):
        out = real(d, f, ev)
        out['naive_bayes'] = out['naive_bayes'].copy(); out['naive_bayes'][0] = np.nan
        return out
    monkeypatch.setattr(R, 'predict_fold', holed)
    with pytest.raises(R.Refusal, match='non-finite probability'):
        R.run_manifest(mp)


def test_no_real_data_command():
    assert 'real' not in open(R.__file__).read().split("if __name__ == '__main__':")[1]


def _spy(monkeypatch):
    seen = []
    real = R.read_table
    def spy(m, base, key):
        seen.append(key)
        return real(m, base, key)
    monkeypatch.setattr(R, 'read_table', spy)
    return seen


def test_non_synthetic_refused_before_any_file_parsed_or_hashed(small, monkeypatch):
    m, t, folder = small
    mp = R.write_dataset(dict(m, synthetic=False), t, folder / 'real')
    seen = _spy(monkeypatch)
    monkeypatch.setattr(R, 'sha256_file', lambda p: (_ for _ in ()).throw(AssertionError('hashed')))
    with pytest.raises(R.Refusal, match='non-synthetic'):
        R.run_manifest(mp)
    assert seen == []


def test_endpoint_group_count_refused_before_gene_matrix(small, monkeypatch):
    m, t, folder = small
    mk = t['mask'].copy()
    mk.loc[(mk.target == 'T2') & ~mk.group_id.isin(['A02', 'A03']), 'assayed'] = 0
    lab = t['labels'].merge(t['cells'][['cell_id', 'group_id']], on='cell_id')
    lab = lab[(lab.target == 'T1') | lab.group_id.isin(['A02', 'A03'])].drop(columns='group_id')
    mp = R.write_dataset(m, dict(t, mask=mk, labels=lab), folder / 'two_groups')
    seen = _spy(monkeypatch)
    with pytest.raises(R.Refusal, match='only 2 eligible biological groups'):
        R.run_manifest(mp)
    assert 'genes' not in seen and set(seen) == {'cells', 'labels', 'mask'}


def test_mask_content_refused_before_gene_matrix(small, monkeypatch):
    m, t, folder = small
    mk = t['mask'].copy(); mk.loc[0, 'assayed'] = 2
    mp = R.write_dataset(m, dict(t, mask=mk), folder / 'badmask')
    seen = _spy(monkeypatch)
    with pytest.raises(R.Refusal, match='assayed must be 0/1'):
        R.run_manifest(mp)
    assert 'genes' not in seen


def test_genes_read_last_on_clean_run(small, monkeypatch):
    m, t, folder = small
    mp = R.write_dataset(m, t, folder / 'order')
    seen = _spy(monkeypatch)
    R.run_manifest(mp)
    assert seen == ['cells', 'labels', 'mask', 'genes']


@pytest.mark.parametrize('bad', [np.float64(0.3), np.array([0.3]), np.full((5, 1), 0.3),
                                 np.array([0.1, 0.2, 1.2, 0.3, 0.4]), np.array([0.1, -0.1, 0.2, 0.3, 0.4]),
                                 np.array([0.1, np.inf, 0.2, 0.3, 0.4]), [0.1, 0.2, 0.3, 0.4, 0.5],
                                 np.array(['a'] * 5)])
def test_checked_probability_rejects_adversarial(bad):
    with pytest.raises(R.Refusal):
        R.checked_probability(bad, 5, 'test')


def test_checked_probability_accepts_valid():
    assert R.checked_probability(np.array([0.0, 0.5, 1.0]), 3, 'ok').dtype == float


@pytest.mark.parametrize('scalar', [0.5, 1.7])
def test_scalar_or_out_of_range_method_output_fails_run(small, monkeypatch, scalar):
    m, t, folder = small
    mp = R.write_dataset(m, t, folder / f'scalar{scalar}')
    real = R.predict_fold
    def bad(d, f, ev):
        out = real(d, f, ev)
        out['technical'] = np.float64(scalar) if scalar <= 1 else np.full(int(ev.sum()), scalar)
        return out
    monkeypatch.setattr(R, 'predict_fold', bad)
    with pytest.raises(R.Refusal):
        R.run_manifest(mp)


@pytest.mark.parametrize('name', ['Gfpt1', 'Gfpt2', 'Creb1', 'Crem', 'Tomm20'])
def test_endogenous_lookalikes_admitted(name):
    assert R.forbidden_name(name) is None


@pytest.mark.parametrize('name', ['EGFP', 'eGFP', 'GFP', 'tdTomato', 'mCherry', 'Cre', 'H2B-GFP', 'WPRE',
                                  'reporter_umi', 'target_LHA', 'class_id', 'subclass', 'cluster_id',
                                  'captured', 'projection_LHA', 'barcode_count', 'LHA'])
def test_reporter_and_metadata_names_refused(name):
    assert R.forbidden_name(name, targets=['LHA']) is not None


@pytest.mark.parametrize('role,value', [
    ('position_columns', ['EGFP', 'y_um', 'z_um']),
    ('position_columns', ['x_um', 'class_y', 'z_um']),
    ('technical_columns', ['tdTomato_counts']),
    ('endogenous_genes', ['Gene00', 'target_T1']),
    ('position_columns', ['x_um', 'x_um', 'z_um']),
    ('technical_columns', ['Gene00']),                  # overlaps the gene role
    ('targets', ['T1', 'T1']),
    ('targets', ['T1', ' ']),
    ('endogenous_genes', ['Gene00', '']),
])
def test_predictor_roles_checked(small, role, value):
    m, t, folder = small
    mp = R.write_dataset(m, t, folder / 'roles')
    full = json.loads(mp.read_text())
    full[role] = value
    with pytest.raises(R.Refusal):
        R.validate_manifest(full, mp.parent)


def test_blank_ids_refused(small):
    m, t, folder = small
    c = t['cells'].copy(); c.loc[3, 'cell_id'] = ''
    mp = R.write_dataset(m, dict(t, cells=c), folder / 'blank')
    with pytest.raises(R.Refusal, match='blank'):
        R.run_manifest(mp)


@pytest.mark.parametrize('animals', [['A00 '], [' A00'], [''], ['   ']])
def test_blank_or_padded_animal_ids_refused(small, animals):
    m, t, folder = small
    mp = R.write_dataset(m, t, folder / 'animals')
    full = json.loads(mp.read_text())
    full['groups'][1]['animals'] = animals          # 'A00 ' would otherwise look disjoint from 'A00'
    with pytest.raises(R.Refusal, match='animal identities'):
        R.validate_manifest(full, mp.parent)


def test_end_to_end_invariant_to_csv_row_order(small):
    m, t, folder = small
    _, d1, r1 = R.run_manifest(R.write_dataset(m, t, folder / 'ordered'))
    shuffled = {k: v.sample(frac=1, random_state=111).reset_index(drop=True) for k, v in t.items()}
    _, d2, r2 = R.run_manifest(R.write_dataset(m, shuffled, folder / 'shuffled'))
    assert list(d1['cell_id']) == list(d2['cell_id'])
    assert r1['macro'] == r2['macro']                      # exact, every method
    for ep in r1['per_endpoint']:
        a, b = r1['per_endpoint'][ep], r2['per_endpoint'][ep]
        assert a['eligible_cell_ids_sha256'] == b['eligible_cell_ids_sha256']
        assert json.dumps(a['methods'], sort_keys=True) == json.dumps(b['methods'], sort_keys=True)
        for mth in R.METHODS:
            assert np.array_equal(a['_preds'][mth], b['_preds'][mth], equal_nan=True)
    assert r1['folds'] == r2['folds'] and r1['fold_log'] == r2['fold_log']
