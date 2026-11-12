import csv

import numpy as np
import pytest

from analysis import paired_projection_nonlinear_replay as R


def observation(group, target, y, ps, fold=0):
    return {'group': group, 'target': target, 'y': y, 'assayed': y >= 0, 'fold': fold,
            'methods': {n: {'p': p, 'raw': None} for n, p in ps.items()}}


def test_animal_equal_then_target_equal_is_not_pooled():
    data = {('a0', 'T1'): observation('A', 'T1', 1, {'m': .1}),
            ('b0', 'T1'): observation('B', 'T1', 1, {'m': .9}),
            ('b1', 'T1'): observation('B', 'T1', 1, {'m': .9}),
            ('a0', 'T2'): observation('A', 'T2', 0, {'m': .2}),
            ('b0', 'T2'): observation('B', 'T2', 0, {'m': .2})}
    score = R.score_case(data)['m']
    assert score['brier_target_equal'] == pytest.approx(((.81+.01)/2+.04)/2)
    expected = ((-np.log(.1)-np.log(.9))/2-np.log(.8))/2
    assert score['log_loss_target_equal'] == pytest.approx(expected)


def test_alpha_selection_guard_and_larger_alpha_tiebreak():
    f = R.FAMILIES[0]
    scores = {f+'|alpha=0.0': {'log_loss_target_equal': .5, 'brier_target_equal': .1},
              f+'|alpha=0.1': {'log_loss_target_equal': .3, 'brier_target_equal': .11},
              f+'|alpha=0.25': {'log_loss_target_equal': .4, 'brier_target_equal': .09}}
    selected = R.select(scores, f)
    assert selected['alpha'] == .25
    assert selected['guard_rejected_unconstrained_best'] is True
    scores[f+'|alpha=0.1'] = dict(scores[f+'|alpha=0.25'])
    assert R.select(scores, f)['alpha'] == .25


def test_gate_rejects_one_control_or_brier_tradeoff():
    scores = {'m': {'log_loss_target_equal': .3, 'brier_target_equal': .1},
              'old': {'log_loss_target_equal': .4, 'brier_target_equal': .11},
              'depth': {'log_loss_target_equal': .35, 'brier_target_equal': .10},
              'copy': {'log_loss_target_equal': .4, 'brier_target_equal': .09}}
    g = R.gate(scores, 'm', 'old', 'depth', ['copy'])
    assert g['met'] is False and set(g['failures']) == {'copy'}


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def row(**changes):
    result = dict(cell_id='x', group_id='A', endpoint='T', fold='0', method='copy',
                  raw_probability='', probability='.2', captured='1')
    result.update(changes)
    return result


def test_duplicate_forecast_and_missing_method_rejected(tmp_path):
    p = tmp_path/'p.csv'
    write_csv(p, [row(), row()])
    with pytest.raises(ValueError, match='duplicate'):
        R.parse_rows(p)
    write_csv(p, [row(), row(method='other'), row(cell_id='y')])
    with pytest.raises(ValueError, match='missing method'):
        R.parse_rows(p)


def test_baseline_replay_requires_identity_not_only_probability(tmp_path):
    data = {('x', 'T'): observation('A', 'T', 1, {'copy': .2})}
    p = tmp_path/'old.csv'
    write_csv(p, [dict(cell_id='x', group_id='B', endpoint='T', fold='0', config='copy',
                       probability='.2', captured='1')])
    with pytest.raises(ValueError, match='identity'):
        R.verify_old_oof(data, p, ['copy'])


def test_shrinkage_checks_raw_and_train_prior():
    f = R.FAMILIES[0]
    obs = observation('A', 'T', 1, {})
    obs['methods'] = {f+'|alpha=0.0': {'p': .8, 'raw': .8},
                      f+'|alpha=0.25': {'p': .7, 'raw': .8}}
    data = {('x', 'T'): obs}
    audit = {'per_endpoint_fold': {'T|0': {f: {'q_fit': .4}}}}
    assert R.check_shrinkage(data, audit) < 1e-15
    audit['per_endpoint_fold']['T|0'][f]['q_fit'] = .5
    with pytest.raises(ValueError, match='shrunk'):
        R.check_shrinkage(data, audit)


def test_hidden_null_bound_is_expected_not_empirical():
    data = {('x', 'T'): observation('A', 'T', 1, {'constant': .5, 'lucky': 1.0})}
    empirical = R.score_case(data)
    expected = R.null_expectation(data)
    assert empirical['lucky']['brier_target_equal'] == 0
    assert expected['lucky']['expected_brier'] == .5
    assert expected['constant']['expected_brier'] == .25
    assert expected['constant']['expected_log_loss'] == pytest.approx(np.log(2))


def test_case_seed_identity_rejects_collision(tmp_path):
    p = tmp_path/'p.csv'
    r = dict(case_id='location_only|876211', scenario='location_only', seed='876212',
             cell_id='x', animal_id='A', target='T', assayed='True', observed_label='1',
             method='copy', raw_probability='', probability='.5')
    write_csv(p, [r])
    with pytest.raises(ValueError, match='case_id'):
        R.parse_rows(p, simulation=True)


def test_completion_hash_detects_changed_forecast(tmp_path):
    import json
    p = tmp_path/'predictions.csv'; p.write_text('original')
    (tmp_path/'RUN_COMPLETE.json').write_text(json.dumps({'files_sha256': {'predictions.csv': R.sha(p)}}))
    R.verify_completion(tmp_path)
    p.write_text('changed')
    with pytest.raises(ValueError, match='hash differs'):
        R.verify_completion(tmp_path)


def test_train_adapter_and_prior_are_independently_weighted():
    x = np.array([[np.nan, 1.], [2., 1.], [4., 1.], [6., 1.]])
    y = np.array([1, 1, 1, 0]); groups = np.array(['A', 'A', 'A', 'B'])
    audit = {'q_fit': .5, 'n_train': 4, 'n_train_positive': 3, 'train_groups': ['A', 'B'],
             'adapter': {'median': [4., 1.], 'scaler_mean': [4., 1.],
                         'scaler_scale': [np.sqrt(2), 1.]}}
    assert R.check_training_audit(audit, x, y, groups) == 0
    audit['q_fit'] = .75
    with pytest.raises(ValueError, match='number differs'):
        R.check_training_audit(audit, x, y, groups)


def test_all_method_descriptives_keep_seed_values():
    scores = {f'{scenario}|{seed}': {'m': {'brier_target_equal': .1+i*.01,
                                         'log_loss_target_equal': .3+i*.02}}
              for scenario, seeds in R.NEW_SEEDS.items() for i, seed in enumerate(seeds)}
    result = R.descriptive_metrics(scores)
    first = result['per_scenario_all_methods']['location_only']['m']
    assert first['brier']['mean'] == pytest.approx(.11)
    assert first['log_loss']['min'] == .3 and len(first['log_loss']['values']) == 3
