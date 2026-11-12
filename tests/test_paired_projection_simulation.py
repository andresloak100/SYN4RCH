"""Leakage, null-law, provenance and metric checks; no required model wins."""
import inspect
import json
from pathlib import Path

import numpy as np
import pytest

from analysis import paired_projection_simulation as S


@pytest.fixture(scope='module')
def fixed():
    anchor, recipe = S.anchors()  # cached authorized TRAIN metadata, never real fit
    d = S.generate('batch_shift_missing', 90817, anchor, n_train=3, n_test=2,
                   cells_per_animal=48, n_genes=8)
    p, audit = S.fit_predict(d, recipe)
    return anchor, recipe, d, p, audit


def test_frozen_protocol_and_exact_development_settings():
    protocol = json.loads(S.PROTOCOL.read_text())
    _, recipe = S.anchors()
    assert protocol['status'] == 'FROZEN_SYNTHETIC_SENSITIVITY_ONLY'
    assert protocol['design'] == dict(train_animals=6, test_animals=4, cells_per_animal=100, genes=16, targets=4)
    assert tuple(protocol['scenarios']) == S.SCENARIOS
    assert protocol['seeds'] == dict(zip(S.SCENARIOS, (764101, 764102, 764103, 764104)))
    assert protocol['pins'] == S.PINS
    names = [c['name'] for c in S.method_configs(recipe)]
    assert names == protocol['methods'] and len(names) == 13
    assert set(recipe['controls_frozen_with_it'].values()).issubset(names)
    assert 'location_logistic|C=0.1|alpha=0.0' in names
    assert 'matched_knn|k=31|alpha=0.1' in names
    assert 'location_knn|k=31|alpha=0.5' in names


def test_training_anchor_is_metadata_only_and_not_gene_projection_coupling(fixed):
    anchor = fixed[0]
    assert anchor['source_rows'] == 37 and anchor['source_animal_proxies'] == 32
    assert anchor['region_profile_values_read'] is False
    assert anchor['coupling_anchor'].startswith('NONE')
    assert anchor['classification'] == S.LABEL
    source = inspect.getsource(S.anchors)
    assert 'np.load' not in source and 'all_neurons_projection_length' not in source
    assert set(Path(k).suffix for k in S.PINS) == {'.json'}


def test_generator_is_reproducible_and_groups_disjoint(fixed):
    anchor, _, d, _, _ = fixed
    replay = S.generate('batch_shift_missing', 90817, anchor, n_train=3, n_test=2,
                        cells_per_animal=48, n_genes=8)
    for key in d:
        if d[key].dtype.kind == 'f':
            assert np.array_equal(d[key], replay[key], equal_nan=True)
        else:
            assert np.array_equal(d[key], replay[key])
    train = set(d['animal'][d['split'] == 'train'])
    test = set(d['animal'][d['split'] == 'test'])
    assert len(train) == 3 and len(test) == 2 and not train & test
    assert len(set(d['cell_id'])) == len(d['cell_id'])
    assert all(name.startswith('SYN_') for name in d['cell_id'])


def test_null_response_excludes_molecular_and_hidden_states():
    rng = np.random.default_rng(18)
    z = rng.normal(size=(40, 3))
    a = S._nonlinear_probability('location_only', z, np.zeros(40, int), np.zeros(40, int))
    b = S._nonlinear_probability('location_only', z, np.full(40, 3), np.full(40, 2))
    assert np.array_equal(a, b)
    # Nonlinear law is not a logistic plane: radial centres have a local peak.
    centre = np.array([[-.9, -.4, .3]])
    left, right = centre + [-2., 0., 0.], centre + [2., 0., 0.]
    p = [S._nonlinear_probability('location_only', v, np.zeros(1, int), np.zeros(1, int))[0, 0]
         for v in (left, centre, right)]
    assert p[1] > max(p[0], p[2])


def test_hidden_state_nonidentifiability_is_an_exact_mixture():
    z = np.zeros((4, 3))
    p = S._nonlinear_probability('latent_unobserved', z, np.arange(4), np.arange(4))
    assert np.allclose(p.mean(axis=0), .5)
    assert np.allclose(np.sort(np.unique(p)), [.08, .92])
    changed_genes = S._nonlinear_probability('latent_unobserved', z + 100., np.arange(4)[::-1], np.arange(4))
    assert np.array_equal(p, changed_genes)


def test_missing_and_unassayed_are_explicit_not_zero_labels(fixed):
    d = fixed[2]
    assert np.isnan(d['genes']).any() and np.isnan(d['xyz']).any()
    assert np.array_equal(d['y'] == -1, ~d['assayed'])
    assert set(np.unique(d['y'])) == {-1, 0, 1}
    assert np.nanmin(d['genes']) >= 0


def test_imputer_and_scaler_ignore_evaluation_distribution():
    train = np.array([[1., np.nan, 4.], [3., np.nan, np.nan], [5., np.nan, 8.]])
    a, _, fit1 = S.fit_view(train, np.array([[np.nan, 5., 100.]]))
    b, _, fit2 = S.fit_view(train, np.array([[1e8, np.nan, -1e8]]))
    assert fit1 == fit2 and np.array_equal(a, b)
    assert fit1['median'] == [3., 0., 6.]
    assert np.allclose(fit1['scaler_mean'], [3., 0., 6.])


def test_predictions_and_training_artifacts_ignore_test_labels_and_oracles(fixed):
    _, recipe, d, p1, audit1 = fixed
    changed = {key: value.copy() for key, value in d.items()}
    ev = changed['split'] == 'test'
    changed['y'][ev] = -1
    changed['assayed'][ev] = False
    changed['generating_probability'][:] = .999
    changed['true_projection'][:] = 1
    changed['latent_molecular'][:] = 0
    changed['latent_hidden'][:] = 0
    p2, audit2 = S.fit_predict(changed, recipe)
    assert audit1 == audit2
    assert all(np.array_equal(p1[k], p2[k]) for k in p1)
    assert audit1['test_labels_read_by_fitter'] is False


def test_animal_overlap_is_refused_before_fitting(fixed):
    _, recipe, d, _, _ = fixed
    bad = dict(d, animal=d['animal'].copy())
    bad['animal'][bad['split'] == 'test'] = bad['animal'][bad['split'] == 'train'][0]
    with pytest.raises(ValueError, match='overlapping'):
        S.fit_predict(bad, recipe)


def test_fixed_candidate_equals_manual_train_only_recipe(fixed):
    from sklearn.linear_model import LogisticRegression
    _, _, d, probabilities, audit = fixed
    tr = (d['split'] == 'train') & (d['y'][:, 0] >= 0)
    ev = d['split'] == 'test'
    x = np.column_stack([np.log1p(d['genes']), d['xyz']])
    a, b, _ = S.fit_view(x[tr], x[ev])
    w = S.group_weights(d['animal'][tr])
    model = LogisticRegression(**S.LOGIT, C=.1).fit(a, d['y'][tr, 0], sample_weight=w)
    q = np.average(d['y'][tr, 0], weights=w)
    want = .9 * model.predict_proba(b)[:, 1] + .1 * q
    assert np.array_equal(want, probabilities['candidate_logistic|C=0.1|alpha=0.1'][:, 0])
    assert audit['per_target']['SYN_T1']['q_fit'] == q
    for target in audit['per_target'].values():
        assert all(not (set(groups) & set(audit['test_animals'])) for groups in target['donors'].values())


def test_scoring_is_animal_equal_target_equal_and_excludes_unknown():
    d = dict(split=np.array(['test'] * 4), animal=np.array(['A', 'A', 'A', 'B']),
             y=np.array([[0, 0, 0, -1], [0, 0, 0, -1], [0, 0, 0, 0], [1, 1, 1, 1]]),
             target_names=np.array(['T1', 'T2', 'T3', 'T4']))
    score = S.score(d, {'method': np.full((4, 4), .2)})['method']
    assert score['brier_target_equal'] == pytest.approx((.2 ** 2 + .8 ** 2) / 2)
    assert score['log_loss_target_equal'] == pytest.approx((-np.log(.8) - np.log(.2)) / 2)
    assert score['per_target']['T4']['per_animal']['A']['n_assayed'] == 1


def test_all_probabilities_finite_and_raw_copy_anchors_exist(fixed):
    probabilities = fixed[3]
    for value in probabilities.values():
        assert value.shape == (96, 4) and np.isfinite(value).all()
        assert ((value >= 0) & (value <= 1)).all()
    for key in ['copy_1nn_location_raw|alpha=0.0', 'copy_1nn_matched_raw|alpha=0.0']:
        assert set(np.unique(probabilities[key])).issubset({0., 1.})
