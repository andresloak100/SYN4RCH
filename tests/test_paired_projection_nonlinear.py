import csv
import json

import numpy as np
import pytest

from analysis import paired_projection_logloss_development as D
from analysis import paired_projection_nonlinear as N
from analysis import paired_projection_real_benchmark as B
from analysis import paired_projection_reference as REF
from analysis import paired_projection_simulation as SIM


# ---------------- fixed structure, fitter, leakage

def test_fixed_structure_and_grid():
    assert N.HGBT == dict(loss='log_loss', learning_rate=0.05, max_iter=100, max_leaf_nodes=4, max_depth=2,
                          min_samples_leaf=20, l2_regularization=1.0, early_stopping=False,
                          categorical_features=None, class_weight=None, random_state=20261112)
    assert N.ALPHAS == (0.0, 0.1, 0.25) and N.FAMILIES == ('hgbt_gene_position', 'hgbt_position')
    proto = json.loads(N.PROTOCOL.read_text())
    assert proto['candidate']['fixed_parameters'] == N.HGBT and tuple(proto['candidate']['alpha_grid']) == N.ALPHAS
    assert proto['comparators']['a1_keys'] == N.A1_COMPARATORS
    assert {s: tuple(v) for s, v in proto['simulation_phase']['new_seeds'].items()} == N.NEW_SEEDS


def _toy(n=300, seed=0):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 3))
    y = (rng.random(n) < 1 / (1 + np.exp(-(x[:, 0] * x[:, 1] * 2)))).astype(int)
    g = np.repeat(['a', 'b', 'c'], n // 3)
    return x, y, g


def test_fitter_shrinkage_q_and_determinism():
    x, y, g = _toy()
    r1, s1, a1 = N.hgbt_fit_predict(x[:200], y[:200], g[:200], x[200:])
    r2, s2, a2 = N.hgbt_fit_predict(x[:200], y[:200], g[:200], x[200:])
    assert np.array_equal(r1, r2) and a1 == a2
    q = float(np.average(y[:200], weights=REF.group_balanced_weights(g[:200])))
    assert a1['q_fit'] == q
    for a in N.ALPHAS:
        assert np.allclose(s1[a], (1 - a) * r1 + a * q, rtol=0, atol=0)


def test_fitter_preprocessing_invariant_to_evaluation_rows():
    x, y, g = _toy()
    xe = x[200:].copy()
    r1, _, a1 = N.hgbt_fit_predict(x[:200], y[:200], g[:200], xe)
    xe2 = xe.copy(); xe2[:50] = 1e6; xe2[50:60] = np.nan                # change other evaluation rows only
    r2, _, a2 = N.hgbt_fit_predict(x[:200], y[:200], g[:200], xe2)
    assert a1['adapter'] == a2['adapter'] and np.array_equal(r1[60:], r2[60:])


def test_fitter_refuses_one_class():
    x, y, g = _toy()
    with pytest.raises(REF.Refusal):
        N.hgbt_fit_predict(x[:50], np.zeros(50, int), g[:50], x[50:60])


# ---------------- selection, paired, gate (hand cases)

def _s(ll, br, pa=None):
    return {'log_loss_target_equal': ll, 'brier_target_equal': br,
            'per_animal_macro': pa or {'A': {'log_loss': ll, 'brier': br}}, 'per_target': {'T': {'log_loss': ll, 'brier': br}}}


def test_select_alpha_guard_and_ties():
    f = 'hgbt_gene_position'
    sc = {N.key(f, 0.0): _s(0.50, 0.100), N.key(f, 0.1): _s(0.40, 0.101), N.key(f, 0.25): _s(0.45, 0.099)}
    r = N.select_alpha(sc, f)
    assert r['selected'] == N.key(f, 0.25) and r['unconstrained_best'] == N.key(f, 0.1) and r['guard_rejected_unconstrained_best']
    sc = {N.key(f, 0.0): _s(0.45, 0.100), N.key(f, 0.1): _s(0.45, 0.099), N.key(f, 0.25): _s(0.45, 0.099)}
    assert N.select_alpha(sc, f)['selected'] == N.key(f, 0.25)            # LL tie, Brier tie -> larger alpha
    sc[N.key(f, 0.1)] = _s(0.45, 0.098)
    assert N.select_alpha(sc, f)['selected'] == N.key(f, 0.1)             # LL tie -> lower Brier


def test_gate_requires_both_metrics_against_every_comparator():
    sc = {'cand': _s(0.30, 0.09), 'lr': _s(0.31, 0.09), 'depth': _s(0.35, 0.10), 'knn': _s(0.29, 0.08)}
    g = N.gate(sc, 'cand', 'lr', 'depth', ['knn'])
    assert not g['met'] and set(g['failures']) == {'knn'} and 'retained' in g['decision']
    sc['knn'] = _s(0.32, 0.09)
    assert N.gate(sc, 'cand', 'lr', 'depth', ['knn'])['met']
    sc['lr'] = _s(0.31, 0.0899)                                             # LR has lower Brier -> fail
    assert set(N.gate(sc, 'cand', 'lr', 'depth', ['knn'])['failures']) == {'lr'}


def test_paired_arithmetic():
    sc = {'c': _s(0.3, 0.1, {'A': {'log_loss': 0.3, 'brier': 0.1}, 'B': {'log_loss': 0.2, 'brier': 0.05}}),
          'x': _s(0.4, 0.12, {'A': {'log_loss': 0.5, 'brier': 0.13}, 'B': {'log_loss': 0.1, 'brier': 0.04}})}
    p = N.paired(sc, 'c', ['x'])['x']
    assert abs(p['log_loss_diff'] - 0.1) < 1e-15 and abs(p['per_animal_macro_diff']['B']['log_loss'] + 0.1) < 1e-15


# ---------------- A1 phase on a synthetic two-animal dataset

@pytest.fixture(scope='module')
def a1(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('a1')
    mp = D.synthetic_dataset(tmp / 'data', sizes=(120, 160))
    dev_dir, _, _ = D.execute(mp, out_root=tmp / 'dev', run_id='dev_for_nonlinear')
    run_dir, result, scores = N.execute_a1(manifest_path=mp, out_root=tmp / 'nl', run_id='nl_test_a1',
                                           saved_predictions=dev_dir / 'predictions.csv')
    return tmp, mp, dev_dir, run_dir, result, scores


def test_a1_comparators_replay_saved_development_predictions(a1):
    _, _, _, _, result, _ = a1
    r = result['baseline_replay']
    assert r['max_abs_diff'] <= 1e-12 and r['keys_checked'] > 0


def test_a1_comparator_scores_equal_saved_development_scores(a1):
    _, _, dev_dir, _, _, scores = a1
    dev = json.loads((dev_dir / 'development_result.json').read_text())['all_configurations']
    for k in N.A1_COMPARATORS:
        assert abs(scores[k]['log_loss_target_equal'] - dev[k]['log_loss_target_equal']) < 1e-12
        assert abs(scores[k]['brier_target_equal'] - dev[k]['brier_target_equal']) < 1e-12


def test_a1_replay_detects_any_changed_probability(a1, tmp_path):
    tmp, mp, dev_dir, _, _, _ = a1
    rows = list(csv.DictReader((dev_dir / 'predictions.csv').open(newline='')))
    for r in rows:
        if r['config'] == N.FROZEN_LR:
            r['probability'] = repr(float(r['probability']) + 1e-9)
            break
    bad = tmp_path / 'bad.csv'
    with bad.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n'); w.writeheader(); w.writerows(rows)
    with pytest.raises(REF.Refusal, match='baseline replay failed'):
        N.execute_a1(manifest_path=mp, out_root=tmp_path / 'o', run_id='bad_replay', saved_predictions=bad)


def test_a1_never_scores_deployment_artifact(a1, tmp_path, monkeypatch):
    _, mp, dev_dir, _, _, _ = a1
    for name in ('deploy_fit', 'deploy_predict'):
        monkeypatch.setattr(D, name, lambda *a, **k: (_ for _ in ()).throw(AssertionError('artifact used')))
    N.execute_a1(manifest_path=mp, out_root=tmp_path / 'o', run_id='no_artifact',
                 saved_predictions=dev_dir / 'predictions.csv')


def test_a1_outputs_serialized_and_replayable(a1):
    _, _, _, run_dir, result, scores = a1
    done = json.loads((run_dir / 'RUN_COMPLETE.json').read_text())
    for name, h in done['files_sha256'].items():
        assert REF.sha256_file(run_dir / name) == h
    sel = json.loads((run_dir / 'selection.json').read_text())
    assert sel['hgbt'] == N.HGBT and set(sel['alphas']) == set(N.FAMILIES)
    cand = result['candidate']
    per = {}
    with (run_dir / 'predictions.csv').open(newline='') as f:
        for r in csv.DictReader(f):
            if r['method'] == cand:
                assert r['raw_probability'] != ''
                p = min(max(float(r['probability']), REF.CLIP[0]), REF.CLIP[1]); y = int(r['captured'])
                per.setdefault(r['endpoint'], {}).setdefault(r['group_id'], []).append(
                    -(y * np.log(p) + (1 - y) * np.log(1 - p)))
    replay = np.mean([np.mean([np.mean(v) for v in g.values()]) for g in per.values()])
    assert abs(replay - scores[cand]['log_loss_target_equal']) < 1e-12


def test_a1_reports_all_configurations_and_gate(a1):
    _, _, _, run_dir, result, scores = a1
    res = json.loads((run_dir / 'result.json').read_text())
    assert set(res['all_configurations']) == set(N.A1_COMPARATORS + N.all_hgbt_keys())
    assert set(result['gate']['required']) >= {N.FROZEN_LR, *N.STRONG_COPY_AVG}
    assert set(result['paired_vs_candidate']) == set(N.A1_COMPARATORS + [result['selection']['hgbt_position']['selected']])


def test_a1_duplicate_run_refused(a1):
    tmp, mp, dev_dir, _, _, _ = a1
    with pytest.raises(REF.Refusal, match='already exists'):
        N.execute_a1(manifest_path=mp, out_root=tmp / 'nl', run_id='nl_test_a1', saved_predictions=dev_dir / 'predictions.csv')


def _no_read(monkeypatch):
    seen = []
    monkeypatch.setattr(B, 'read_table', lambda *a, **k: seen.append(a) or (_ for _ in ()).throw(AssertionError('read')))
    return seen


def _use_draft_protocol(monkeypatch, tmp_path):
    """Exercise a draft explicitly even after the production protocol is frozen."""
    draft = json.loads(N.PROTOCOL.read_text())
    draft['status'] = 'DRAFT_FOR_ROOT_REVIEW_NO_RUNS'
    path = tmp_path / 'draft_protocol.json'
    path.write_text(json.dumps(draft))
    check = N.check_common
    monkeypatch.setattr(N, 'check_common', lambda repo=None: check(repo=repo or N.ROOT, protocol_path=path))


def test_real_a1_refused_while_draft(monkeypatch, tmp_path):
    _use_draft_protocol(monkeypatch, tmp_path)
    seen = _no_read(monkeypatch)
    before = (N.OUT / N.A1_RUN_ID).exists()
    with pytest.raises(REF.Refusal, match=N.FROZEN_STATUS):
        N.execute_a1()
    assert seen == [] and (N.OUT / N.A1_RUN_ID).exists() == before


@pytest.mark.parametrize('kw', [{'run_id': 'nonlinear_a1_dev_v2'}, {'out_root': 'TMP'}])
def test_real_a1_identity_bound(tmp_path, monkeypatch, kw):
    if kw.get('out_root') == 'TMP':
        kw = {'out_root': tmp_path / 'x'}
    seen = _no_read(monkeypatch)
    with pytest.raises(REF.Refusal, match='canonical output root'):
        N.execute_a1(**kw)
    assert seen == []


def test_non_synthetic_manifest_override_refused(tmp_path):
    m = tmp_path / 'm.json'; m.write_text(json.dumps({'synthetic': False}))
    with pytest.raises(REF.Refusal, match='must be synthetic'):
        N.execute_a1(manifest_path=m, out_root=tmp_path, run_id='x')


# ---------------- simulation phase (old exposed seeds / tiny arbitrary seeds only)

@pytest.fixture(scope='module')
def anchor_recipe():
    return SIM.anchors()


def test_new_seeds_refused_without_authorization(anchor_recipe):
    anchor, _ = anchor_recipe
    for s, seeds in N.NEW_SEEDS.items():
        for seed in seeds:
            with pytest.raises(REF.Refusal, match='only generated by the authorized'):
                N.generate_case(s, seed, anchor)


def test_real_sim_refused_while_draft(monkeypatch, tmp_path):
    _use_draft_protocol(monkeypatch, tmp_path)
    monkeypatch.setattr(SIM, 'generate', lambda *a, **k: (_ for _ in ()).throw(AssertionError('generated')))
    before = (N.OUT / N.SIM_RUN_ID).exists()
    with pytest.raises(REF.Refusal, match=N.FROZEN_STATUS):
        N.execute_sim()
    assert (N.OUT / N.SIM_RUN_ID).exists() == before


def test_sim_identity_bound(tmp_path):
    with pytest.raises(REF.Refusal, match='canonical output root'):
        N.execute_sim(out_root=tmp_path)
    with pytest.raises(REF.Refusal, match='canonical output root'):
        N.execute_sim(run_id='nonlinear_sim_v2')


def test_generator_ids_repeat_across_seeds_so_case_id_needed(anchor_recipe):
    anchor, _ = anchor_recipe
    a = SIM.generate('location_only', 11, anchor, n_train=1, n_test=1, cells_per_animal=5)
    b = SIM.generate('location_only', 12, anchor, n_train=1, n_test=1, cells_per_animal=5)
    assert list(a['cell_id']) == list(b['cell_id']) and not np.array_equal(a['genes'], b['genes'])


@pytest.mark.parametrize('scenario', ['complementarity', 'batch_shift_missing'])
def test_sim_case_ignores_test_labels_masks_and_oracles(anchor_recipe, scenario):
    anchor, recipe = anchor_recipe
    seed = json.loads(SIM.PROTOCOL.read_text())['seeds'][scenario]          # old exposed seed
    d = SIM.generate(scenario, seed, anchor)
    alphas = {'hgbt_gene_position': 0.1, 'hgbt_position': 0.0}
    p1, r1, a1 = N.sim_case(d, recipe, alphas)
    d2 = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in d.items()}
    te = d2['split'] == 'test'
    rng = np.random.default_rng(3)
    d2['y'][te] = rng.integers(-1, 2, d2['y'][te].shape)
    d2['assayed'][te] = ~d2['assayed'][te]
    for k in ('generating_probability', 'true_projection', 'latent_molecular', 'latent_hidden'):
        d2[k] = rng.permutation(d2[k])
    p2, r2, a2 = N.sim_case(d2, recipe, alphas)
    assert set(p1) == set(p2) and all(np.array_equal(p1[k], p2[k]) for k in p1)
    assert all(np.array_equal(r1[k], r2[k]) for k in r1)
    assert json.dumps(a1['hgbt'], sort_keys=True) == json.dumps(a2['hgbt'], sort_keys=True)
    d3 = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in d.items()}
    tr = d3['split'] == 'train'
    d3['y'][tr] = np.where(d3['y'][tr] >= 0, 1 - d3['y'][tr], d3['y'][tr])  # TRAIN labels do matter
    p3, _, _ = N.sim_case(d3, recipe, alphas)
    assert not np.array_equal(p1[N.key('hgbt_gene_position', 0.1)], p3[N.key('hgbt_gene_position', 0.1)])


def test_sim_case_comparators_equal_frozen_simulation(anchor_recipe):
    anchor, recipe = anchor_recipe
    d = SIM.generate('location_only', json.loads(SIM.PROTOCOL.read_text())['seeds']['location_only'], anchor)
    p, _, _ = N.sim_case(d, recipe, {'hgbt_gene_position': 0.0, 'hgbt_position': 0.0})
    frozen, _ = SIM.fit_predict(d, recipe)
    assert all(np.array_equal(p[k], frozen[k]) for k in frozen)
    assert set(frozen) == set(N.A1_TO_SIM.values())


def test_recipes_pinned_separately_and_consistent():
    proto = json.loads(N.PROTOCOL.read_text())
    pub, arch = N.check_recipes(proto['pins'])
    assert REF.sha256_file(N.ROOT / proto['pins']['published_recipe']) != proto['pins']['archived_candidate_recipe_sha256']
    assert pub['controls_frozen_with_it'] == arch['controls_frozen_with_it']


def test_sim_summary_arithmetic():
    def c(cand_ll, x_ll):
        return {'scores': {'hgbt_gene_position|alpha=0.1': {'log_loss_target_equal': cand_ll, 'brier_target_equal': 0.1},
                           'x': {'log_loss_target_equal': x_ll, 'brier_target_equal': 0.12}}}
    cases = {f'{s}|{i}': c(0.3, 0.3 + 0.01 * i) for s in SIM.SCENARIOS for i in (1, 2, 3)}
    s = N.sim_summary(cases)
    assert abs(s['overall']['x']['log_loss_diff_mean'] - 0.02) < 1e-12 and s['overall']['x']['n_cases'] == 12
    assert abs(s['per_scenario']['location_only']['x']['log_loss_diff_min'] - 0.01) < 1e-12


def test_frozen_files_untouched():
    pins = json.loads(N.PROTOCOL.read_text())['pins']
    assert REF.sha256_file(B.__file__) == pins['frozen_runner_sha256']
    assert REF.sha256_file(REF.__file__) == pins['reference_sha256']
    assert REF.sha256_file(D.__file__) == pins['development_code_sha256']
    assert REF.sha256_file(SIM.__file__) == pins['simulation_code_sha256']
    assert REF.sha256_file(SIM.PROTOCOL) == pins['simulation_protocol_sha256']


# ---------------- deployment artifact

def test_artifact_roundtrip_hash_verified_and_all_rows(a1):
    tmp, mp, _, run_dir, result, _ = a1
    m, d = D.load(mp)
    sel = json.loads((run_dir / 'selection.json').read_text())
    for fam in N.FAMILIES:
        info = sel['artifacts'][fam]
        meta, models = N.load_artifact(run_dir / info['metadata'], info['metadata_sha256'])
        assert meta['kind'] == 'DEPLOYMENT_REFIT_NOT_OOF' and meta['alpha'] == result['selection'][fam]['alpha']
        assert meta['hgbt'] == N.HGBT
        for ep, t in meta['targets'].items():
            e = d['y'][ep] >= 0
            assert t['n_fit'] == int(e.sum()) and set(t['train_groups']) == set(d['group'])
            assert t['train_cell_ids'] == sorted(d['cell_id'][e])
        meta2, models2 = N.deploy_fit(d, fam, meta['alpha'], m['endogenous_genes'], m['position_columns'])
        g = np.expm1(d['genes'][:7]) if fam == 'hgbt_gene_position' else None
        gn = m['endogenous_genes'] if fam == 'hgbt_gene_position' else []
        a = N.artifact_predict(meta, models, gn, m['position_columns'], g, d['xyz'][:7])
        b = N.artifact_predict(dict(meta2, model_file='x', model_sha256='x'), models2, gn, m['position_columns'], g, d['xyz'][:7])
        assert all(np.allclose(a[k], b[k], rtol=0, atol=1e-12) for k in a)


def test_artifact_refuses_tampered_files(a1, tmp_path):
    import shutil
    _, _, _, run_dir, _, _ = a1
    sel = json.loads((run_dir / 'selection.json').read_text())
    info = sel['artifacts']['hgbt_gene_position']
    with pytest.raises(REF.Refusal, match='metadata hash'):
        N.load_artifact(run_dir / info['metadata'], '0' * 64)
    for name in (info['metadata'], info['model_file']):
        shutil.copy(run_dir / name, tmp_path / name)
    with (tmp_path / info['model_file']).open('ab') as f:
        f.write(b'tamper')
    with pytest.raises(REF.Refusal, match='model file hash'):
        N.load_artifact(tmp_path / info['metadata'])


@pytest.fixture(scope='module')
def art(a1):
    _, mp, _, run_dir, _, _ = a1
    m, d = D.load(mp)
    info = json.loads((run_dir / 'selection.json').read_text())['artifacts']
    return m, d, {f: N.load_artifact(run_dir / info[f]['metadata'], info[f]['metadata_sha256']) for f in N.FAMILIES}


@pytest.mark.parametrize('case', ['strings', 'bool', 'object', '1d', 'wrong_width', 'rows', 'nan', 'inf',
                                  'negative', 'empty', 'panel', 'genes_for_position'])
def test_artifact_predict_strict_inputs(art, case):
    m, d, arts = art
    fam = 'hgbt_position' if case == 'genes_for_position' else 'hgbt_gene_position'
    meta, models = arts[fam]
    gn = list(m['endogenous_genes']) if fam == 'hgbt_gene_position' else []
    g, p = np.expm1(d['genes'][:4]), d['xyz'][:4].copy()
    if case == 'strings':
        g = np.full(g.shape, 'a')
    elif case == 'bool':
        g = g > 1
    elif case == 'object':
        g = g.astype(object)
    elif case == '1d':
        p = p[:, 0]
    elif case == 'wrong_width':
        g = g[:, :-1]
    elif case == 'rows':
        p = p[:3]
    elif case == 'nan':
        g = g.copy(); g[0, 0] = np.nan
    elif case == 'inf':
        p[1, 0] = np.inf
    elif case == 'negative':
        g = g.copy(); g[2, 1] = -1
    elif case == 'empty':
        g, p = g[:0], p[:0]
    elif case == 'panel':
        gn = gn[::-1]
    with pytest.raises(REF.Refusal):
        N.artifact_predict(meta, models, gn, m['position_columns'], g, p)


def test_artifact_predict_features_only_and_row_independent(art):
    import inspect
    m, d, arts = art
    meta, models = arts['hgbt_gene_position']
    g, p = np.expm1(d['genes'][:6]), d['xyz'][:6].copy()
    a = N.artifact_predict(meta, models, m['endogenous_genes'], m['position_columns'], g, p)
    g2, p2 = g.copy(), p.copy(); g2[3:] *= 50; p2[3:] += 999
    b = N.artifact_predict(meta, models, m['endogenous_genes'], m['position_columns'], g2, p2)
    assert all(np.array_equal(a[k][:3], b[k][:3]) for k in a)
    assert not any('label' in n or n in ('y', 'captured') for n in inspect.signature(N.artifact_predict).parameters)


# ---------------- selection binding and the stage-2 selection record

def test_selection_binds_sources_snapshot_and_outputs(a1):
    _, _, _, run_dir, _, _ = a1
    sel = json.loads((run_dir / 'selection.json').read_text())
    assert sel['source_sha256'] == N.source_hashes()
    assert sel['protocol_snapshot_sha256'] == REF.sha256_file(run_dir / 'protocol_snapshot.json')
    assert REF.sha256_file(run_dir / 'protocol_snapshot.json') == REF.sha256_file(N.PROTOCOL)
    for k, f in (('result_sha256', 'result.json'), ('predictions_sha256', 'predictions.csv'),
                 ('run_complete_sha256', 'RUN_COMPLETE.json')):
        assert sel[k] == REF.sha256_file(run_dir / f)
    started = json.loads((run_dir / 'STARTED.json').read_text())
    assert started['protocol_snapshot_sha256'] == sel['protocol_snapshot_sha256']


def _git(repo, *a):
    import subprocess
    return subprocess.run(['git', '-C', str(repo), *a], capture_output=True, text=True, check=True)


def _record_repo(tmp_path, run_dir, alphas=None, commit=True):
    sel = json.loads((run_dir / 'selection.json').read_text())
    r = tmp_path / 'repo'; (r / 'docs').mkdir(parents=True)
    rec = {'schema': 'synarch_nonlinear_selection_record_v1',
           'a1_selection_sha256': REF.sha256_file(run_dir / 'selection.json'),
           'a1_run_complete_sha256': sel['run_complete_sha256'], 'alphas': alphas or sel['alphas']}
    (r / 'docs/sel.json').write_text(json.dumps(rec))
    _git(r, 'init', '-q'); _git(r, 'config', 'user.email', 't@example.invalid'); _git(r, 'config', 'user.name', 't')
    if commit:
        _git(r, 'add', 'docs/sel.json'); _git(r, 'commit', '-q', '-m', 's')
    else:
        _git(r, 'commit', '-q', '--allow-empty', '-m', 'e')
    return r


def test_selection_record_accepts_matching_committed_record(a1, tmp_path):
    _, _, _, run_dir, result, _ = a1
    r = _record_repo(tmp_path, run_dir)
    alphas = N.check_selection_record(r / 'docs/sel.json', run_dir, r)
    assert alphas == {f: result['selection'][f]['alpha'] for f in N.FAMILIES}


@pytest.mark.parametrize('case', ['uncommitted', 'alphas', 'a1_file_changed', 'source_changed', 'selection_hash'])
def test_selection_record_refusals(a1, tmp_path, monkeypatch, case):
    import shutil
    _, _, _, run_dir, _, _ = a1
    a1_copy = tmp_path / 'a1'; shutil.copytree(run_dir, a1_copy)
    alphas = None
    if case == 'alphas':
        alphas = {'hgbt_gene_position': 0.25 if json.loads((run_dir / 'selection.json').read_text())['alphas']['hgbt_gene_position'] != 0.25 else 0.0,
                  'hgbt_position': 0.0}
    r = _record_repo(tmp_path, a1_copy, alphas=alphas, commit=case != 'uncommitted')
    if case == 'a1_file_changed':
        with (a1_copy / 'predictions.csv').open('a') as f:
            f.write('x\n')
    elif case == 'source_changed':
        real = N.source_hashes()
        monkeypatch.setattr(N, 'source_hashes', lambda: dict(real, simulation='0' * 64))
    elif case == 'selection_hash':
        rec = json.loads((r / 'docs/sel.json').read_text()); rec['a1_selection_sha256'] = '0' * 64
        (r / 'docs/sel.json').write_text(json.dumps(rec)); _git(r, 'commit', '-qam', 'x')
    with pytest.raises(REF.Refusal):
        N.check_selection_record(r / 'docs/sel.json', a1_copy, r)


def test_sim_refused_without_selection_record_even_if_protocol_frozen(monkeypatch, tmp_path):
    monkeypatch.setattr(N, 'check_common', lambda repo=None: {})
    monkeypatch.setattr(SIM, 'generate', lambda *a, **k: (_ for _ in ()).throw(AssertionError('generated')))
    monkeypatch.setattr(N, 'SELECTION_RECORD', tmp_path / 'missing_selection.json')
    before = (N.OUT / N.SIM_RUN_ID).exists()
    with pytest.raises(REF.Refusal, match='selection record'):
        N.execute_sim()
    assert (N.OUT / N.SIM_RUN_ID).exists() == before


# ---------------- common checks

def test_common_refuses_substitute_source(tmp_path):
    proto = json.loads(N.PROTOCOL.read_text()); proto['status'] = N.FROZEN_STATUS
    pp = tmp_path / 'p.json'; pp.write_text(json.dumps(proto))
    with pytest.raises(REF.Refusal, match='canonical'):
        N.check_common(protocol_path=pp, source_path=tmp_path / 'paired_projection_nonlinear.py')


@pytest.mark.parametrize('field', ['params', 'alphas', 'run_id', 'seeds', 'comparators'])
def test_protocol_constants_must_match_code(field):
    proto = json.loads(N.PROTOCOL.read_text())
    if field == 'params':
        proto['candidate']['fixed_parameters']['max_depth'] = 3
    elif field == 'alphas':
        proto['candidate']['alpha_grid'] = [0.0, 0.1]
    elif field == 'run_id':
        proto['simulation_phase']['run_id'] = 'nonlinear_sim_v2'
    elif field == 'seeds':
        proto['simulation_phase']['new_seeds']['location_only'] = [1, 2, 3]
    elif field == 'comparators':
        proto['comparators']['a1_keys'] = proto['comparators']['a1_keys'][:-1]
    with pytest.raises(REF.Refusal):
        N.check_protocol_constants(proto)
    assert N.check_protocol_constants(json.loads(N.PROTOCOL.read_text()))


def test_recipe_pins_include_published_hash(tmp_path):
    pins = dict(json.loads(N.PROTOCOL.read_text())['pins'], published_recipe_sha256='0' * 64)
    with pytest.raises(REF.Refusal, match='published recipe'):
        N.check_recipes(pins)


# ---------------- strict saved-prediction verification

@pytest.mark.parametrize('case', ['duplicate_and_missing', 'extra_key', 'wrong_fold', 'wrong_group', 'prob_range',
                                  'prob_nan', 'label', 'columns'])
def test_verify_saved_strict(a1, tmp_path, case):
    tmp, mp, dev_dir, _, _, _ = a1
    rows = list(csv.DictReader((dev_dir / 'predictions.csv').open(newline='')))
    cols = list(rows[0])
    idx = [i for i, r in enumerate(rows) if r['config'] in N.A1_COMPARATORS]
    i0, i1 = idx[0], idx[1]
    if case == 'duplicate_and_missing':
        rows[i1] = dict(rows[i0])                         # same count: one key duplicated, one missing
    elif case == 'extra_key':
        rows.append(dict(rows[i0], cell_id='NOT_A_CELL'))
    elif case == 'wrong_fold':
        rows[i0]['fold'] = str(1 - int(rows[i0]['fold']))
    elif case == 'wrong_group':
        rows[i0]['group_id'] = 'OTHER'
    elif case == 'prob_range':
        rows[i0]['probability'] = '1.5'
    elif case == 'prob_nan':
        rows[i0]['probability'] = 'nan'
    elif case == 'label':
        rows[i0]['captured'] = str(1 - int(rows[i0]['captured']))
    elif case == 'columns':
        cols = cols[::-1]
    bad = tmp_path / 'bad.csv'
    with bad.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, lineterminator='\n'); w.writeheader(); w.writerows(rows)
    with pytest.raises(REF.Refusal):
        N.execute_a1(manifest_path=mp, out_root=tmp_path / 'o', run_id=f'bad_{case}', saved_predictions=bad)


# ---------------- simulation case files

def test_case_npz_no_pickle_roundtrip(anchor_recipe, tmp_path):
    anchor, _ = anchor_recipe
    d = SIM.generate('batch_shift_missing', json.loads(SIM.PROTOCOL.read_text())['seeds']['batch_shift_missing'], anchor)
    p = tmp_path / 'case.npz'
    N.save_case_npz(p, d, 'batch_shift_missing|764104', 764104)
    with np.load(p, allow_pickle=False) as z:
        assert set(z.files) >= {'case_id', 'seed', 'y', 'assayed', 'xyz', 'genes', 'tech', 'cell_id', 'animal',
                                'split', 'generating_probability', 'true_projection', 'latent_molecular', 'latent_hidden'}
        assert (z['case_id'] == 'batch_shift_missing|764104').all() and np.isnan(z['genes']).any()
        assert np.array_equal(z['y'], d['y'])
    with pytest.raises(FileExistsError):
        N.save_case_npz(p, d, 'x', 1)
