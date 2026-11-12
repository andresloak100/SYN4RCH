import json
import shutil
import subprocess

import numpy as np
import pytest

from analysis import paired_projection_logloss_development as D
from analysis import paired_projection_real_benchmark as B
from analysis import paired_projection_reference as REF

REAL_MANIFEST = D.ROOT / 'data/interim/paired_projection_real_benchmark/barseq2_a1_inputs/manifest.json'


@pytest.fixture(scope='module')
def dev(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('dev')
    mp = D.synthetic_dataset(tmp / 'data', sizes=(120, 170))
    run_dir, summary, parts = D.execute(mp, out_root=tmp / 'runs', run_id='t_dev_01')
    return mp, run_dir, summary, parts, tmp


def test_shrink_formula():
    p = np.array([0.0, 0.2, 1.0])
    assert np.allclose(D.shrink(p, 0.3, 0.25), 0.75 * p + 0.25 * 0.3)
    assert np.array_equal(D.shrink(p, 0.3, 0.0), p)


def test_grid_is_fixed_and_complete(dev):
    _, _, summary, (m, d, grid, scores, chosen, art), _ = dev
    assert D.C_GRID == (0.01, 0.1, 1.0) and D.K_GRID == (5, 15, 31) and D.ALPHA_GRID == (0.0, 0.1, 0.25, 0.5)
    assert len(grid['keys']) == 82 and set(summary['all_configurations']) == set(grid['keys'])
    for ep, g in grid['preds'].items():
        e = d['y'][ep] >= 0
        for k in grid['keys']:
            assert np.isfinite(g['p'][k][e]).all() and ((g['p'][k][e] >= 0) & (g['p'][k][e] <= 1)).all()


def test_one_global_setting_per_family_no_target_choice(dev):
    _, _, summary, (m, d, grid, scores, chosen, art), _ = dev
    assert set(chosen) == set(grid['families'])
    for f, v in chosen.items():
        assert D.parse_key(v['selected'])['family'] == f
        if f in D.LOGISTIC_FAMILIES + D.AVERAGING_FAMILIES:
            assert v['n_configs'] == 12
        elif f in D.COPY_FAMILIES:
            assert v['n_configs'] == 4
        assert not any(t in v['selected'] for t in d['endpoints'] if len(t) > 1)


def test_candidate_baseline_equals_frozen_runner(dev):
    _, _, _, (m, d, grid, scores, chosen, art), _ = dev
    for ep in d['endpoints']:
        y = d['y'][ep]; e = y >= 0
        frozen = np.full(len(y), np.nan)
        for evg in B.folds_for(m, d['group']):
            tr = e & ~np.isin(d['group'], list(evg)); ev = e & np.isin(d['group'], list(evg))
            frozen[ev] = B.fit_predict(d, tr, ev, y[tr])[0]['candidate_logistic']
        assert np.array_equal(frozen[e], grid['preds'][ep]['p'][D.cfg_key('candidate_logistic', C=1.0, alpha=0.0)][e])


def test_q_fit_and_fits_invariant_to_hidden_eval_labels(dev):
    _, _, _, (m, d, grid, scores, chosen, art), _ = dev
    ep = d['endpoints'][0]
    y = d['y'][ep]; e = y >= 0
    evg = B.folds_for(m, d['group'])[0]
    tr = e & ~np.isin(d['group'], list(evg)); ev = e & np.isin(d['group'], list(evg))
    a, _, qa, _ = D.fit_predict_grid(d, tr, ev, y[tr])
    d2 = dict(d, genes=d['genes'].copy(), xyz=d['xyz'].copy(), tech=d['tech'].copy())
    d2['genes'][ev] += 7.0; d2['xyz'][ev] *= 3; d2['tech'][ev] -= 2
    b, _, qb, _ = D.fit_predict_grid(d2, tr, ev, y[tr])
    assert qa == qb == float(np.average(y[tr], weights=REF.group_balanced_weights(d['group'][tr])))
    assert np.array_equal(a[D.cfg_key('train_frequency')], b[D.cfg_key('train_frequency')])


def test_selection_rule_guard_and_tie_break():
    def sc(ll, br):
        return {'log_loss_target_equal': ll, 'brier_target_equal': br, 'per_animal_macro': {}, 'per_target': {}}
    f = 'location_knn'
    scores = {D.cfg_key(f, k=15, alpha=0.0): sc(0.60, 0.200),                 # own baseline
              D.cfg_key(f, k=5, alpha=0.1): sc(0.50, 0.201),                  # best LL but degrades Brier
              D.cfg_key(f, k=15, alpha=0.25): sc(0.55, 0.199),
              D.cfg_key(f, k=31, alpha=0.5): sc(0.55, 0.199),                 # LL and Brier tie: larger alpha
              D.cfg_key(f, k=5, alpha=0.5): sc(0.55, 0.199)}                  # full tie: smaller k
    ch = D.select(scores, [f])[f]
    assert ch['selected'] == D.cfg_key(f, k=5, alpha=0.5) and ch['n_pass_guard'] == 4
    assert ch['unrestricted_log_loss_best'] == D.cfg_key(f, k=5, alpha=0.1) and ch['guard_rejected_unrestricted_best']
    scores[D.cfg_key(f, k=31, alpha=0.1)] = sc(0.55, 0.198)                  # LL tie: lower Brier beats larger alpha
    assert D.select(scores, [f])[f]['selected'] == D.cfg_key(f, k=31, alpha=0.1)


def test_comparison_hand_arithmetic(dev):
    _, _, summary, (m, d, grid, scores, chosen, art), _ = dev
    comp = summary['comparison']
    cand = comp['candidate_config']
    pair = comp['pairs']['location_logistic|tuned']
    k = pair['control_config']
    assert abs(pair['log_loss_target_equal_diff'] -
               (scores[k]['log_loss_target_equal'] - scores[cand]['log_loss_target_equal'])) < 1e-15
    a = sorted(scores[cand]['per_animal_macro'])[0]
    hand = np.mean([scores[k]['per_target'][t]['per_animal'][a]['log_loss'] for t in d['endpoints']]) - \
        np.mean([scores[cand]['per_target'][t]['per_animal'][a]['log_loss'] for t in d['endpoints']])
    assert abs(hand - pair['per_animal_macro_diff'][a]['log_loss']) < 1e-12
    assert {p.split('|')[1] for p in comp['pairs']} == {'tuned', 'baseline'}


def test_keyed_predictions_replay(dev):
    import csv
    _, run_dir, summary, (m, d, grid, scores, chosen, art), _ = dev
    k = chosen['candidate_logistic']['selected']
    ep = d['endpoints'][0]
    per = {}
    with (run_dir / 'predictions.csv').open(newline='') as f:
        for r in csv.DictReader(f):
            if r['config'] == k and r['endpoint'] == ep:
                p = min(max(float(r['probability']), REF.CLIP[0]), REF.CLIP[1]); y = int(r['captured'])
                per.setdefault(r['group_id'], []).append(-(y * np.log(p) + (1 - y) * np.log(1 - p)))
    replay = np.mean([np.mean(v) for v in per.values()])
    assert abs(replay - scores[k]['per_target'][ep]['log_loss']) < 1e-12


def test_duplicate_run_refused(dev):
    mp, run_dir, _, _, tmp = dev
    with pytest.raises(REF.Refusal, match='already exists'):
        D.execute(mp, out_root=tmp / 'runs', run_id='t_dev_01')


def test_row_shuffle_invariant_selection(tmp_path):
    m, t = B.make_synthetic(sizes=(110, 150), n_genes=6, dims=1, unassayed_t2=False)
    a = D.execute(B.write_dataset(m, t, tmp_path / 'a'), out_root=tmp_path / 'r', run_id='sa')[1]
    s = {k: v.sample(frac=1, random_state=111).reset_index(drop=True) for k, v in t.items()}
    b = D.execute(B.write_dataset(m, s, tmp_path / 'b'), out_root=tmp_path / 'r', run_id='sb')[1]
    assert json.dumps(a['all_configurations'], sort_keys=True) == json.dumps(b['all_configurations'], sort_keys=True)
    assert a['selected'] == b['selected']


# ---------------- real development run refusals (no exposed data parsed)

def _no_read(monkeypatch):
    seen = []
    monkeypatch.setattr(B, 'read_table', lambda *a, **k: seen.append(a) or (_ for _ in ()).throw(AssertionError('read')))
    return seen


def test_repo_dev_protocol_frozen():
    assert json.loads(D.PROTOCOL.read_text())['status'] == D.FROZEN_DEV_STATUS


@pytest.mark.parametrize('status', ['DRAFT_FOR_ROOT_REVIEW_NO_DEVELOPMENT_RUN', 'FROZEN', ''])
def test_real_run_refused_with_non_frozen_protocol_copy(tmp_path, monkeypatch, status):
    proto = json.loads(D.PROTOCOL.read_text()); proto['status'] = status
    pp = tmp_path / 'p.json'; pp.write_text(json.dumps(proto))
    seen = _no_read(monkeypatch)
    before = (D.OUT / D.DEV_RUN_ID).exists()
    with pytest.raises(REF.Refusal, match='FROZEN_DEVELOPMENT_AUTHORIZED'):
        D.execute(REAL_MANIFEST, protocol_path=pp)
    assert seen == [] and (D.OUT / D.DEV_RUN_ID).exists() == before


@pytest.mark.parametrize('kw,msg', [({'run_id': 'logloss_dev_barseq2_a1_v2'}, 'run_id must be'),
                                    ({'out_root': 'TMP'}, 'output root must be'),
                                    ({'code_path': 'COPY'}, 'code_path override refused')])
def test_real_run_identity_and_code_bound(tmp_path, monkeypatch, kw, msg):
    if kw.get('out_root') == 'TMP':
        kw = {'out_root': tmp_path / 'elsewhere'}
    if kw.get('code_path') == 'COPY':
        c = tmp_path / 'dev_copy.py'; shutil.copy(D.__file__, c); kw = {'code_path': c}
    seen = _no_read(monkeypatch)
    before = (D.OUT / D.DEV_RUN_ID).exists()
    with pytest.raises(REF.Refusal, match=msg):
        D.execute(REAL_MANIFEST, **kw)
    assert seen == [] and (D.OUT / D.DEV_RUN_ID).exists() == before and not (tmp_path / 'elsewhere').exists()


def _git(repo, *a):
    return subprocess.run(['git', '-C', str(repo), *a], capture_output=True, text=True, check=True)


def _pinned_repo(tmp_path, manifest_sha, admission_sha, commit_code=True):
    r = tmp_path / 'repo'
    (r / 'docs').mkdir(parents=True); (r / 'analysis').mkdir()
    proto = json.loads(D.PROTOCOL.read_text())
    proto['status'] = D.FROZEN_DEV_STATUS
    proto['data_status']['pins'].update(input_manifest_sha256=manifest_sha, admission_sha256=admission_sha,
                                        frozen_runner_sha256=REF.sha256_file(B.__file__),
                                        reference_sha256=REF.sha256_file(REF.__file__))
    (r / 'docs/p.json').write_text(json.dumps(proto))
    shutil.copy(D.__file__, r / 'analysis/dev.py')
    _git(r, 'init', '-q'); _git(r, 'config', 'user.email', 't@example.invalid'); _git(r, 'config', 'user.name', 't')
    _git(r, 'add', 'docs/p.json');
    if commit_code:
        _git(r, 'add', 'analysis/dev.py')
    _git(r, 'commit', '-q', '-m', 'p')
    return r


def _fake_manifest(tmp_path):
    d = tmp_path / 'm'; d.mkdir()
    (d / 'adm.json').write_text('{"status": "ADMITTED"}')
    (d / 'manifest.json').write_text(json.dumps({'synthetic': False, 'admission': {'record_path': 'adm.json'}}))
    return d / 'manifest.json', REF.sha256_file(d / 'manifest.json'), REF.sha256_file(d / 'adm.json')


def test_pins_pass_when_everything_matches(tmp_path):
    mp, msha, asha = _fake_manifest(tmp_path)
    r = _pinned_repo(tmp_path, msha, asha)
    assert D.check_dev_pins(mp, r, r / 'docs/p.json', r / 'analysis/dev.py')['status'] == D.FROZEN_DEV_STATUS


@pytest.mark.parametrize('case', ['manifest_pin', 'admission_pin', 'code_uncommitted', 'protocol_modified'])
def test_pins_refuse(tmp_path, case):
    mp, msha, asha = _fake_manifest(tmp_path)
    r = _pinned_repo(tmp_path, '0' * 64 if case == 'manifest_pin' else msha,
                     '0' * 64 if case == 'admission_pin' else asha, commit_code=case != 'code_uncommitted')
    if case == 'protocol_modified':
        (r / 'docs/p.json').write_text((r / 'docs/p.json').read_text() + ' ')
    with pytest.raises(REF.Refusal):
        D.check_dev_pins(mp, r, r / 'docs/p.json', r / 'analysis/dev.py')


def test_frozen_files_untouched():
    assert REF.sha256_file(B.__file__) == '7b23157a704cd6e34eb887a7b9ec59f3e2050e86435f9cd174efb0ec1250baac'
    assert REF.sha256_file(REF.__file__) == '9ae3dd30a7fad55c6b7d232e8faad4ec678bf406c105faa61cade544039723a3'
    assert REF.sha256_file(B.PROTOCOL) == 'd40320c3f47ac5158c6bec3f734ab67d56816d60ba12207f4d5dca352960259b'


# ---------------- deployment artifact: refit on all training animals; features-only predictions

def test_deploy_refit_matches_manual_fit(dev):
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    _, _, _, (m, d, grid, scores, chosen, art), _ = dev
    ep = d['endpoints'][0]
    y = d['y'][ep]; e = y >= 0
    X = np.hstack([d['genes'], d['xyz']])[e]
    sc = StandardScaler().fit(X)
    lr = LogisticRegression(**dict(B.LOGIT, C=art['C'])).fit(sc.transform(X), y[e],
                                                              sample_weight=REF.group_balanced_weights(d['group'][e]))
    t = art['targets'][ep]
    assert np.array_equal(np.array(t['coef']), lr.coef_[0]) and t['n_fit'] == int(e.sum())
    assert set(t['n_fit_per_animal']) == set(d['group'])                          # both training animals


def test_future_predictions_ignore_future_labels(dev):
    import inspect
    _, _, _, (m, d, grid, scores, chosen, art), _ = dev
    m2, t2 = B.make_synthetic(sizes=(70, 80), n_genes=8, dims=1, seed=4242, unassayed_t2=False)
    genes = t2['genes'][m['endogenous_genes']].to_numpy(float)
    pos = t2['cells'][m['position_columns']].to_numpy(float)
    p1 = D.deploy_predict(art, m['endogenous_genes'], m['position_columns'], genes, pos)
    t2['labels']['captured'] = 1 - t2['labels']['captured']                      # flip every future label
    p2 = D.deploy_predict(art, m['endogenous_genes'], m['position_columns'], genes, pos)
    assert all(np.array_equal(p1[k], p2[k]) for k in p1)
    assert not any('label' in n or n in ('y', 'captured') for n in inspect.signature(D.deploy_predict).parameters)


def test_deploy_predict_matches_shrunk_logistic(dev):
    _, _, _, (m, d, grid, scores, chosen, art), _ = dev
    ep = d['endpoints'][0]
    t = art['targets'][ep]
    genes = np.expm1(d['genes'][:5]); pos = d['xyz'][:5]
    z = ((np.hstack([d['genes'][:5], pos]) - t['scaler_mean']) / t['scaler_scale']) @ t['coef'] + t['intercept']
    want = (1 - art['alpha']) / (1 + np.exp(-z)) + art['alpha'] * t['q_fit']
    got = D.deploy_predict(art, m['endogenous_genes'], m['position_columns'], genes, pos)[ep]
    assert np.allclose(got, want, rtol=0, atol=1e-12)


def test_different_panel_refused(dev):
    _, _, _, (m, d, grid, scores, chosen, art), _ = dev
    genes = list(m['endogenous_genes']); genes[0] = 'Hcrt'
    with pytest.raises(REF.Refusal, match='no A1 coefficient transfer'):
        D.deploy_predict(art, genes, m['position_columns'], np.ones((2, len(genes))), np.ones((2, 1)))


def test_artifact_written_and_hashed(dev):
    _, run_dir, _, (m, d, grid, scores, chosen, art), _ = dev
    done = json.loads((run_dir / 'RUN_COMPLETE.json').read_text())
    assert done['deployment_artifact_a1.json_sha256'] == REF.sha256_file(run_dir / 'deployment_artifact_a1.json')
    assert json.loads((run_dir / 'deployment_artifact_a1.json').read_text())['setting'] == chosen['candidate_logistic']['selected']


@pytest.mark.parametrize('case', ['strings', 'object', 'bool', '1d_genes', '3d_pos', 'wrong_gene_width',
                                  'wrong_pos_width', 'row_mismatch', 'nan', 'inf_pos', 'negative_gene', 'empty'])
def test_deploy_predict_strict_inputs(dev, case):
    _, _, _, (m, d, grid, scores, chosen, art), _ = dev
    g, p = np.expm1(d['genes'][:4]), d['xyz'][:4]
    ng = len(m['endogenous_genes'])
    if case == 'strings':
        g = np.full((4, ng), 'a')
    elif case == 'object':
        g = g.astype(object)
    elif case == 'bool':
        g = g > 1
    elif case == '1d_genes':
        g = g[0]
    elif case == '3d_pos':
        p = p[:, :, None]
    elif case == 'wrong_gene_width':
        g = g[:, :-1]
    elif case == 'wrong_pos_width':
        p = np.hstack([p, p])
    elif case == 'row_mismatch':
        p = p[:3]
    elif case == 'nan':
        g = g.copy(); g[0, 0] = np.nan
    elif case == 'inf_pos':
        p = p.copy(); p[1, 0] = np.inf
    elif case == 'negative_gene':
        g = g.copy(); g[2, 3] = -1
    elif case == 'empty':
        g, p = g[:0], p[:0]
    with pytest.raises(REF.Refusal):
        D.deploy_predict(art, m['endogenous_genes'], m['position_columns'], g, p)


def test_deploy_predict_stable_for_extreme_inputs(dev):
    _, _, _, (m, d, grid, scores, chosen, art), _ = dev
    g = np.full((3, len(m['endogenous_genes'])), 1e12); p = np.array([[1e9], [-1e9], [0.0]])
    out = D.deploy_predict(art, m['endogenous_genes'], m['position_columns'], g, p)
    for v in out.values():
        assert np.isfinite(v).all() and ((v >= 0) & (v <= 1)).all()
