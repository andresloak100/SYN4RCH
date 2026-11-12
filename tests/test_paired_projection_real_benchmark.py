import json
import shutil
import subprocess

import numpy as np
import pytest

from analysis import paired_projection_real_benchmark as B
from analysis import paired_projection_reference as REF

REVIEWED_REFERENCE_SHA = '9ae3dd30a7fad55c6b7d232e8faad4ec678bf406c105faa61cade544039723a3'


@pytest.fixture(scope='module')
def syn():
    return B.make_synthetic(sizes=(60, 110, 80, 130, 90), n_genes=6)


def write(syn, tmp_path, name='d', **over):
    m, t = syn
    m = json.loads(json.dumps(m)); m.update(over)
    return B.write_dataset(m, t, tmp_path / name)


def spy(monkeypatch):
    seen = []
    real = B.read_table
    monkeypatch.setattr(B, 'read_table', lambda m, base, key: (seen.append(key), real(m, base, key))[1])
    return seen


def test_reviewed_reference_unchanged():
    assert REF.sha256_file(REF.__file__) == REVIEWED_REFERENCE_SHA


def test_clean_synthetic_logo(syn, tmp_path):
    _, d, res = B.run_manifest(write(syn, tmp_path))
    assert res['gate'] == 'SYNTHETIC_SOFTWARE_ONLY'
    assert res['folds'] == [[g] for g in sorted(set(d['group']))]
    assert set(res['methods']) >= {'train_frequency', 'location_logistic', 'location_knn15', 'molecular_knn15',
                                   'matched_knn15', 'copy_1nn_location', 'copy_1nn_matched'}
    pv = res['per_endpoint']['T1']['paired_vs_candidate']
    assert set(pv) == set(B.LOC_COPY_CONTROLS)


def test_paired_differences_match_hand(syn, tmp_path):
    _, d, res = B.run_manifest(write(syn, tmp_path))
    ep = res['per_endpoint']['T1']
    y = d['y']['T1']; e = y >= 0
    pc, pl = ep['_preds']['candidate_logistic'][e], ep['_preds']['location_logistic'][e]
    g = d['group'][e]
    diffs = [np.mean((pl[g == a] - y[e][g == a]) ** 2) - np.mean((pc[g == a] - y[e][g == a]) ** 2) for a in sorted(set(g))]
    assert abs(np.mean(diffs) - ep['paired_vs_candidate']['location_logistic']['brier_mean']) < 1e-12


def test_copy_controls_are_hard_copies(syn, tmp_path):
    _, _, res = B.run_manifest(write(syn, tmp_path))
    for k in ('copy_1nn_location', 'copy_1nn_matched'):
        p = res['per_endpoint']['T1']['_preds'][k]
        assert set(np.unique(p[np.isfinite(p)])) <= {0.0, 1.0}


def test_donor_from_eval_animal_refused(syn, tmp_path, monkeypatch):
    real = B.fit_predict
    def leaky(d, tr, ev, y_tr):
        out, donors, s = real(d, tr, ev, y_tr)
        donors['copy_1nn_location'] = d['group'][ev]
        return out, donors, s
    monkeypatch.setattr(B, 'fit_predict', leaky)
    with pytest.raises(REF.Refusal, match='donor came from an evaluation animal'):
        B.run_manifest(write(syn, tmp_path))


@pytest.mark.parametrize('bad', ['scalar', 'range', 'nan'])
def test_invalid_probability_fails(syn, tmp_path, monkeypatch, bad):
    real = B.fit_predict
    def broken(d, tr, ev, y_tr):
        out, donors, s = real(d, tr, ev, y_tr)
        n = int(ev.sum())
        out['matched_knn15'] = {'scalar': np.float64(0.5), 'range': np.full(n, 1.5),
                                'nan': np.full(n, np.nan)}[bad]
        return out, donors, s
    monkeypatch.setattr(B, 'fit_predict', broken)
    with pytest.raises(REF.Refusal):
        B.run_manifest(write(syn, tmp_path))


def test_row_shuffle_invariance(syn, tmp_path):
    m, t = syn
    _, _, a = B.run_manifest(B.write_dataset(m, t, tmp_path / 'a'))
    s = {k: v.sample(frac=1, random_state=111).reset_index(drop=True) for k, v in t.items()}
    _, _, b = B.run_manifest(B.write_dataset(m, s, tmp_path / 'b'))
    assert json.dumps(B.strip(a), sort_keys=True) == json.dumps(B.strip(b), sort_keys=True)


def test_frozen_cohort_mismatch_refused_before_genes(syn, tmp_path, monkeypatch):
    m, t = syn
    fc = json.loads(json.dumps(m['frozen_cohort'])); fc['T1']['n'] += 1
    mp = write(syn, tmp_path, frozen_cohort=fc)
    seen = spy(monkeypatch)
    with pytest.raises(REF.Refusal, match='frozen cohort'):
        B.run_manifest(mp)
    assert 'genes' not in seen


def test_labels_changed_after_freeze_refused(syn, tmp_path):
    m, t = syn
    lab = t['labels'].copy(); lab.loc[lab.index[:5], 'qc_pass'] = 1 - lab.loc[lab.index[:5], 'qc_pass']
    mp = B.write_dataset(m, dict(t, labels=lab), tmp_path / 'chg')
    with pytest.raises(REF.Refusal, match='frozen cohort'):
        B.run_manifest(mp)


def test_one_and_two_axis_positions(tmp_path):
    for dims in (1, 2):
        m, t = B.make_synthetic(sizes=(60, 90, 70, 80), n_genes=5, dims=dims, with_tech=False)
        _, d, res = B.run_manifest(B.write_dataset(m, t, tmp_path / f'd{dims}'))
        assert d['xyz'].shape[1] == dims and 'technical' not in res['methods']


@pytest.mark.parametrize('over,msg', [
    ({'position_source': 'tangram_mapped'}, 'measured in the same cell'),
    ({'position_provenance': 'Tangram inferred XYZ'}, 'inferred'),
    ({'position_provenance': 'imputed from atlas'}, 'inferred'),
    ({'position_dims': 2}, 'position_dims'),
    ({'position_measured': False}, 'measured'),
    ({'technical_evidence': ''}, 'technical_evidence'),
    ({'run_class': 'exploratory', 'power_approval': {'approved_by': 'x', 'reference': 'y', 'min_groups': 3}}, 'no fictitious'),
    ({'run_class': 'confirmatory'}, 'power_approval'),
    ({'fold_scheme': {'kind': 'random_cells'}}, 'fold_scheme'),
    ({'fold_scheme': {'kind': 'grouped_kfold', 'k': 1}}, 'grouped_kfold'),
    ({'group_unit': 'slice'}, "'animal' or 'pool'"),
    ({'frozen_cohort': {}}, 'frozen_cohort'),
])
def test_manifest_refusals_before_any_parse(syn, tmp_path, monkeypatch, over, msg):
    mp = write(syn, tmp_path, **over)
    seen = spy(monkeypatch)
    with pytest.raises(REF.Refusal, match=msg):
        B.run_manifest(mp)
    assert seen == []


@pytest.mark.parametrize('role,val', [('position_columns', ['EGFP', 'y_um', 'z_um']),
                                      ('technical_columns', ['cluster_size']),
                                      ('endogenous_genes', ['Gene00', 'tdTomato'])])
def test_forbidden_names_any_role(syn, tmp_path, role, val):
    mp = write(syn, tmp_path, **{role: val})
    with pytest.raises(REF.Refusal):
        B.run_manifest(mp)


def test_padded_and_overlapping_animals(syn, tmp_path):
    m, _ = syn
    g = json.loads(json.dumps(m['groups'])); g[1]['animals'] = ['M00 ']
    with pytest.raises(REF.Refusal, match='animal identities'):
        B.run_manifest(write(syn, tmp_path, 'pad', groups=g))
    g = json.loads(json.dumps(m['groups'])); g[1]['animals'] = ['M00']
    with pytest.raises(REF.Refusal, match='disjoint'):
        B.run_manifest(write(syn, tmp_path, 'ovl', groups=g))


def test_location_only_signal_gives_no_spurious_candidate_gain(tmp_path):
    m, t = B.make_synthetic(sizes=(150, 200, 170, 220, 190), n_genes=6)
    rng = np.random.default_rng(5)
    gtab = t['genes'].copy()
    for c in gtab.columns[1:]:
        gtab[c] = rng.permutation(gtab[c].to_numpy())          # genes carry no label information
    _, _, res = B.run_manifest(B.write_dataset(m, dict(t, genes=gtab), tmp_path / 'loc'))
    for ep in ('T1', 'T2'):
        diff = res['per_endpoint'][ep]['paired_vs_candidate']['location_logistic']['brier_mean']
        assert diff < 0.01                                       # no material gain over location alone


# ---------------- real-run pins (temporary git repo; synthetic files posing as an admitted dataset)

def _git(repo, *a):
    return subprocess.run(['git', '-C', str(repo), *a], capture_output=True, text=True, check=True)


def _repo(tmp_path, status):
    r = tmp_path / f'repo_{status}'
    (r / 'docs').mkdir(parents=True)
    proto = json.loads(B.PROTOCOL.read_text()); proto['status'] = status
    (r / 'docs/protocol.json').write_text(json.dumps(proto))
    _git(r, 'init', '-q'); _git(r, 'config', 'user.email', 't@example.invalid'); _git(r, 'config', 'user.name', 't')
    _git(r, 'add', 'docs/protocol.json'); _git(r, 'commit', '-q', '-m', 'p')
    return r


@pytest.fixture
def repo(tmp_path):
    return _repo(tmp_path, B.FROZEN_STATUS)


def real_manifest(syn, tmp_path, repo, status='ADMITTED', run_id='run_test_01', **over):
    m, t = syn
    mp = B.write_dataset(dict(json.loads(json.dumps(m)), synthetic=False, run_id=run_id), t, tmp_path / 'real')
    man = json.loads(mp.read_text())
    man.update(over)
    proto = repo / 'docs/protocol.json'
    rec = {'status': status, 'dataset_id': man['dataset_id'],
           'file_sha256': {k: man['files'][k]['sha256'] for k in B.FILE_KEYS},
           'analysis_config_sha256': B.analysis_config_sha256(man),
           'runner_sha256': B.sha256_file(B.__file__), 'reference_sha256': B.sha256_file(REF.__file__),
           'protocol_sha256': B.sha256_file(proto), 'source_sha256': {'pooled.mat': 'f' * 64}}
    (mp.parent / 'admission.json').write_text(json.dumps(rec))
    man['admission'] = {'record_path': 'admission.json', 'record_sha256': B.sha256_file(mp.parent / 'admission.json')}
    man['protocol'] = {'sha256': B.sha256_file(proto)}
    mp.write_text(json.dumps(man))
    return mp


def run(mp, repo, out, **kw):
    return B.execute(mp, out_root=out, repo=repo, protocol_path=repo / 'docs/protocol.json', **kw)


def test_real_execute_writes_durable_keyed_outputs(syn, tmp_path, repo):
    mp = real_manifest(syn, tmp_path, repo)
    run_dir, res = run(mp, repo, tmp_path / 'out')
    assert res['gate'].startswith('EXPLORATORY_NOT_POWERED')
    done = json.loads((run_dir / 'RUN_COMPLETE.json').read_text())
    assert done['result_sha256'] == B.sha256_file(run_dir / 'result.json')
    assert done['predictions_sha256'] == B.sha256_file(run_dir / 'predictions.csv')
    assert done['source_sha256'] == {'pooled.mat': 'f' * 64}
    lines = (run_dir / 'predictions.csv').read_text().splitlines()
    n_elig = sum(pt['n_eligible_cells'] for pt in res['per_endpoint'].values())
    assert len(lines) - 1 == n_elig * len(res['methods'])
    # replay: keyed predictions reproduce the animal-equal Brier for one method/endpoint
    import csv
    rows = [r for r in csv.DictReader(lines) if r['endpoint'] == 'T1' and r['method'] == 'location_logistic']
    per = {}
    for r in rows:
        per.setdefault(r['group_id'], []).append((float(r['probability']) - int(r['captured'])) ** 2)
    replay = np.mean([np.mean(v) for v in per.values()])
    assert abs(replay - res['per_endpoint']['T1']['methods']['location_logistic']['brier_animal_equal']) < 1e-12


def test_duplicate_run_id_refused(syn, tmp_path, repo):
    mp = real_manifest(syn, tmp_path, repo)
    run(mp, repo, tmp_path / 'out')
    with pytest.raises(REF.Refusal, match='already exists'):
        run(mp, repo, tmp_path / 'out')


def test_partial_failure_keeps_started_marker_and_blocks_rerun(syn, tmp_path, repo, monkeypatch):
    mp = real_manifest(syn, tmp_path, repo)
    def boom(m, d):
        raise RuntimeError('simulated crash')
    monkeypatch.setattr(B, 'evaluate', boom)
    with pytest.raises(RuntimeError):
        run(mp, repo, tmp_path / 'out')
    rd = tmp_path / 'out' / 'run_test_01'
    assert (rd / 'STARTED.json').exists() and not (rd / 'RUN_COMPLETE.json').exists()
    monkeypatch.undo()
    with pytest.raises(REF.Refusal, match='already exists'):
        run(mp, repo, tmp_path / 'out')


def test_draft_protocol_refused(syn, tmp_path, monkeypatch):
    draft = _repo(tmp_path, 'DRAFT_FOR_ROOT_REVIEW_NO_REAL_RUN_AUTHORIZED')
    mp = real_manifest(syn, tmp_path, draft)
    seen = spy(monkeypatch)
    with pytest.raises(REF.Refusal, match='not FROZEN_EXPLORATORY_AUTHORIZED'):
        run(mp, draft, tmp_path / 'out')
    assert seen == [] and not (tmp_path / 'out').exists()


def test_repo_protocol_frozen_for_one_exploratory_run():
    p = json.loads(B.PROTOCOL.read_text())
    assert p['status'] == B.FROZEN_STATUS
    assert p['dataset']['run_id'] == 'barseq2_a1_598_logo_v1' and len(p['dataset']['targets']) == 11
    assert p['confirmatory_margins'] is None


@pytest.mark.parametrize('status', ['DRAFT_FOR_ROOT_REVIEW_NO_REAL_RUN_AUTHORIZED', 'DRAFT', '', 'FROZEN'])
def test_copies_with_non_frozen_status_refuse(syn, tmp_path, monkeypatch, status):
    r = _repo(tmp_path, status or 'x')
    if not status:
        pj = json.loads((r / 'docs/protocol.json').read_text()); pj.pop('status')
        (r / 'docs/protocol.json').write_text(json.dumps(pj)); _git(r, 'commit', '-qam', 'nostatus')
    mp = real_manifest(syn, tmp_path, r)
    seen = spy(monkeypatch)
    with pytest.raises(REF.Refusal, match='FROZEN_EXPLORATORY_AUTHORIZED'):
        run(mp, r, tmp_path / 'out')
    assert seen == []


@pytest.mark.parametrize('mutate', ['run_id', 'fold_scheme', 'targets_order', 'group_mapping', 'genes', 'frozen_cohort',
                                    'endpoint_mode', 'position_columns'])
def test_config_changed_after_admission_refused(syn, tmp_path, repo, monkeypatch, mutate):
    mp = real_manifest(syn, tmp_path, repo)
    man = json.loads(mp.read_text())
    if mutate == 'run_id':
        man['run_id'] = 'barseq2_a1_598_logo_v2'
    elif mutate == 'fold_scheme':
        man['fold_scheme'] = {'kind': 'grouped_kfold', 'k': 2}
    elif mutate == 'targets_order':
        man['targets'] = list(reversed(man['targets']))
    elif mutate == 'group_mapping':
        man['groups'][0]['animals'] = ['M99']
    elif mutate == 'genes':
        man['endogenous_genes'] = man['endogenous_genes'][:-1]
    elif mutate == 'frozen_cohort':
        man['frozen_cohort']['T2']['sha256'] = '0' * 64
    elif mutate == 'endpoint_mode':
        man['endpoint_mode'] = 'exclusive_two_target'
    elif mutate == 'position_columns':
        man['position_columns'] = ['y_um', 'x_um', 'z_um']
    mp.write_text(json.dumps(man))
    seen = spy(monkeypatch)
    with pytest.raises(REF.Refusal, match='analysis configuration differs'):
        run(mp, repo, tmp_path / 'out')
    assert seen == []


def test_independence_text_is_not_part_of_config(syn, tmp_path, repo):
    m, _ = syn
    a = B.analysis_config_sha256(dict(m, files={k: {'sha256': 'x'} for k in B.FILE_KEYS}))
    g = json.loads(json.dumps(m['groups'])); g[0]['independence_evidence'] = 'other wording'
    b = B.analysis_config_sha256(dict(m, groups=g, files={k: {'sha256': 'x'} for k in B.FILE_KEYS}))
    assert a == b


@pytest.mark.parametrize('which', ['runner', 'reference'])
def test_mutated_code_refused(syn, tmp_path, repo, which):
    mp = real_manifest(syn, tmp_path, repo)
    src = B.__file__ if which == 'runner' else REF.__file__
    mutated = tmp_path / f'{which}.py'
    mutated.write_text(open(src).read() + '\n# changed\n')
    kw = {'runner_path': mutated} if which == 'runner' else {'reference_path': mutated}
    with pytest.raises(REF.Refusal, match='SHA-256 differs'):
        run(mp, repo, tmp_path / 'out', **kw)


@pytest.mark.parametrize('case', ['no_admission', 'not_admitted', 'admission_hash', 'protocol_hash',
                                  'protocol_uncommitted', 'confirmatory'])
def test_real_run_refused_before_any_parse(syn, tmp_path, repo, monkeypatch, case):
    status = 'PENDING' if case == 'not_admitted' else 'ADMITTED'
    over = {'run_class': 'confirmatory', 'power_approval': {'approved_by': 'ROOT', 'reference': 'r', 'min_groups': 4}} \
        if case == 'confirmatory' else {}
    mp = real_manifest(syn, tmp_path, repo, status=status, **over)
    man = json.loads(mp.read_text())
    if case == 'no_admission':
        man.pop('admission')
    elif case == 'admission_hash':
        man['admission']['record_sha256'] = '0' * 64
    elif case == 'protocol_hash':
        man['protocol']['sha256'] = '0' * 64
    elif case == 'protocol_uncommitted':
        p = repo / 'docs/protocol.json'
        p.write_text(p.read_text() + ' ')
        man['protocol']['sha256'] = B.sha256_file(p)
    mp.write_text(json.dumps(man))
    seen = spy(monkeypatch)
    with pytest.raises(REF.Refusal):
        run(mp, repo, tmp_path / 'out')
    assert seen == []


def test_run_manifest_refuses_real(syn, tmp_path, repo):
    mp = real_manifest(syn, tmp_path, repo)
    with pytest.raises(REF.Refusal, match='only through execute'):
        B.run_manifest(mp, repo=repo, protocol_path=repo / 'docs/protocol.json')


def test_macro_hand_arithmetic(syn, tmp_path):
    _, _, res = B.run_manifest(write(syn, tmp_path))
    eps = list(res['per_endpoint'])
    mac = res['macro']
    for k in res['methods']:
        te = np.mean([res['per_endpoint'][e]['methods'][k]['brier_animal_equal'] for e in eps])
        assert abs(te - mac['methods'][k]['brier_target_equal']) < 1e-15
    # per-animal macro with declared coverage (M00 has no T2 in the fixture)
    assert mac['coverage']['M00']['targets'] == ['T1'] and not mac['complete_coverage']
    a = 'M01'
    hand = np.mean([res['per_endpoint'][e]['methods']['candidate_logistic']['per_animal'][a]['log_loss'] for e in eps])
    assert abs(hand - mac['methods']['candidate_logistic']['per_animal_macro'][a]['log_loss']) < 1e-15
    c = 'copy_1nn_location'
    d_hand = mac['methods'][c]['per_animal_macro'][a]['brier'] - mac['methods']['candidate_logistic']['per_animal_macro'][a]['brier']
    assert abs(d_hand - mac['paired_vs_candidate'][c]['per_animal_macro_diff'][a]['brier_diff']) < 1e-15
    te_diff = np.mean([res['per_endpoint'][e]['paired_vs_candidate'][c]['brier_mean'] for e in eps])
    assert abs(te_diff - mac['paired_vs_candidate'][c]['brier_target_equal_diff']) < 1e-15
    assert set(mac['paired_vs_candidate']) == set(B.LOC_COPY_CONTROLS)


def test_two_animal_confirmatory_refused(tmp_path):
    m, t = B.make_synthetic(sizes=(80, 90), n_genes=5, unassayed_t2=False)
    m = dict(m, run_class='confirmatory', power_approval={'approved_by': 'ROOT', 'reference': 'r', 'min_groups': 3})
    with pytest.raises(REF.Refusal, match='fewer groups|two-animal'):
        B.run_manifest(B.write_dataset(m, t, tmp_path / 'two'))
    m2 = dict(m, run_class='exploratory'); m2.pop('power_approval')
    _, _, res = B.run_manifest(B.write_dataset(m2, t, tmp_path / 'two_expl'))
    assert len(res['folds']) == 2 and res['gate'] == 'SYNTHETIC_SOFTWARE_ONLY'


def test_predictions_csv_quotes_ids_with_commas(syn, tmp_path, repo):
    import csv
    m, t = syn
    odd = {c: c.replace('_c', ',"c') for c in t['cells']['cell_id']}           # commas and quotes in IDs
    t2 = {k: v.copy() for k, v in t.items()}
    for k in ('cells', 'genes', 'labels'):
        t2[k]['cell_id'] = t2[k]['cell_id'].map(odd)
    m2 = dict(m, frozen_cohort=B.freeze_cohort(m, t2))
    run_dir, res = run(real_manifest((m2, t2), tmp_path, repo), repo, tmp_path / 'out')
    with (run_dir / 'predictions.csv').open(newline='') as f:
        rows = list(csv.DictReader(f))
    assert all(len(r) == 7 for r in rows)
    ids = {r['cell_id'] for r in rows}
    assert ids <= set(odd.values()) and any(',' in i and '"' in i for i in ids)
    n_elig = sum(pt['n_eligible_cells'] for pt in res['per_endpoint'].values())
    assert len(rows) == n_elig * len(res['methods'])
