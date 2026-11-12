"""Synthetic-only contract tests. No biological library is opened here."""
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("axon_completion_contract_test", ROOT / "analysis/axon_completion_contract.py")
C = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = C
spec.loader.exec_module(C)


def tree(xyz, parent):
    return SimpleNamespace(xyz=np.asarray(xyz, float), parent=np.asarray(parent, int),
                           secret_full_length=123456, types=np.arange(len(parent)))


def straight():
    return tree([[0, 0, 0], [150, 0, 0], [450, 0, 0], [600, 0, 0]], [-1, 0, 1, 2])


def test_exact_edge_boundary_and_suffix_conservation():
    obs, truth = C.split_arbor(straight())
    assert np.array_equal(obs.xyz, [[0, 0, 0], [150, 0, 0], [300, 0, 0]])
    assert np.array_equal(obs.parent, [-1, 0, 1])
    assert np.array_equal(obs.frontier_indices(), [2])
    assert truth.eligible and truth.forest.cable() == 300
    assert np.isclose(C.sample_forest(truth.forest)[1].sum(), 300)
    assert np.isclose(np.linalg.norm(np.diff(obs.xyz, axis=0), axis=1).sum() + truth.forest.cable(), 600)


def test_branched_mask_interpolates_each_boundary_and_hides_descendants():
    source = tree([[0, 0, 0], [200, 0, 0], [400, 0, 0], [200, 300, 0], [700, 0, 0]], [-1, 0, 1, 1, 2])
    obs, truth = C.split_arbor(source)
    assert len(obs.xyz) == 4
    assert any(np.array_equal(p, [300, 0, 0]) for p in obs.xyz)
    assert any(np.array_equal(p, [200, 100, 0]) for p in obs.xyz)
    assert len(obs.frontier_indices()) == 2
    assert not any(np.array_equal(p, [700, 0, 0]) for p in obs.xyz)
    assert truth.forest.cable() == 600


def test_hidden_geometry_swap_and_zero_length_hidden_branching_are_invisible():
    a = tree([[0, 0, 0], [150, 0, 0], [300, 0, 0], [600, 0, 0]], [-1, 0, 1, 2])
    b = tree([[0, 0, 0], [150, 0, 0], [300, 0, 0], [300, 0, 0], [300, 0, 0],
              [300, 900, 0], [300, -900, 0]], [-1, 0, 1, 2, 2, 3, 4])
    assert C.observe(a).to_dict() == C.observe(b).to_dict()


def test_visible_order_does_not_depend_on_source_ids_or_hidden_row_order():
    a = tree([[0, 0, 0], [200, 0, 0], [500, 0, 0], [0, 200, 0], [0, 500, 0]], [-1, 0, 1, 0, 3])
    b = tree([[0, 0, 0], [0, 200, 0], [0, 500, 0], [200, 0, 0], [500, 0, 0]], [-1, 0, 1, 0, 3])
    assert C.observe(a).sha256 == C.observe(b).sha256


def test_public_observation_has_no_hidden_fields_aliases_or_terminal_flags():
    source = straight(); obs = C.observe(source)
    assert set(obs.to_dict()) == {"schema", "xyz_um", "parent", "path_budget_um"}
    assert set(obs.__slots__) == {"xyz", "parent", "path_budget_um"}
    source.xyz[:] = 999
    assert obs.xyz[0, 0] == 0
    with pytest.raises(ValueError):
        obs.xyz[0, 0] = 1
    with pytest.raises(ValueError):
        obs.xyz.setflags(write=True)
    with pytest.raises(ValueError):
        C.Observation.from_dict({**obs.to_dict(), "full_length": 10000})


def test_perfect_completion_scores_only_hidden_cable():
    obs, truth = C.split_arbor(straight())
    pred = C.Completion(obs, [[450, 0, 0], [600, 0, 0]], np.array([2, 3]))
    m = C.score_completion(pred, truth)
    assert m["distance_um"] == 0 and m["f1_200"] == 1
    assert m["valid_for_comparison"]
    assert m["hidden_cable_um"] == m["extra_cable_um"] == 300
    assert m["covered_hidden_cable_200_um"] == 300
    assert np.array_equal(pred.xyz[:len(obs.xyz)], obs.xyz)
    assert np.array_equal(pred.parent[:len(obs.xyz)], obs.parent)


def test_observation_only_is_failed_not_free_recall_or_dropped():
    obs, truth = C.split_arbor(straight())
    m = C.score_completion(C.observation_only(obs), truth)
    assert m["eligible"] and m["empty_forecast"] and not m["valid_for_comparison"]
    assert np.isinf(m["distance_um"]) and m["f1_200"] == 0 and m["recall_200"] == 0
    encoded = C.metrics_to_json(m)
    assert encoded["distance_um"] is None and encoded["status"] == "failed_empty_forecast"
    json.dumps(encoded, allow_nan=False)


def test_extra_connector_near_observation_is_not_removed_or_unscored():
    obs, truth = C.split_arbor(straight())
    pred = C.Completion(obs, [[600, 0, 0], [0, 100, 0]], np.array([2, 0]))
    m = C.score_completion(pred, truth)
    assert m["extra_cable_um"] == 400
    assert m["precision_200"] < 1 and m["unmatched_extra_cable_200_um"] > 0
    assert m["cable_signed_relative_error"] == pytest.approx(1 / 3)


def test_duplicate_cable_cannot_pass_even_if_raw_nearest_point_f1_is_perfect():
    obs, truth = C.split_arbor(straight())
    pred = C.Completion(obs, [[450, 0, 0], [600, 0, 0], [450, 0, 0], [600, 0, 0]], np.array([2, 3, 2, 5]))
    m = C.score_completion(pred, truth)
    assert m["f1_200"] == 1 and m["duplicate_extra_cable_um"] == 300
    assert m["extra_cable_um"] == 600 and m["cable_signed_relative_error"] == 1
    assert not m["nonduplicate_valid"] and not m["valid_for_comparison"]


def test_generated_copy_of_an_observed_edge_is_extra_not_observed():
    obs, truth = C.split_arbor(straight())
    pred = C.Completion(obs, [[600, 0, 0], [150, 0, 0]], np.array([2, 0]))
    m = C.score_completion(pred, truth)
    assert m["extra_cable_um"] == 450 and m["copied_observed_extra_cable_um"] == 150
    assert not m["nonduplicate_valid"]


def test_forests_have_no_fictitious_connecting_cable():
    f = C.CableForest(np.array([[0, 0, 0], [50, 0, 0], [1000, 1000, 0], [1050, 1000, 0]]),
                      np.array([-1, 0, -1, 2]), np.array([300, 350, 300, 350]))
    x, w = C.sample_forest(f)
    assert f.cable() == w.sum() == 100
    assert not np.any((x[:, 0] > 50) & (x[:, 0] < 1000))


def test_global_path_grid_and_mass_are_segmentation_invariant():
    a = C.CableForest(np.array([[0, 0, 0], [175, 0, 0]]), np.array([-1, 0]), np.array([125, 300]))
    b = C.CableForest(np.array([[0, 0, 0], [25, 0, 0], [80, 0, 0], [175, 0, 0]]),
                      np.array([-1, 0, 1, 2]), np.array([125, 150, 205, 300]))
    ax, aw = C.sample_forest(a); bx, bw = C.sample_forest(b)
    assert np.array_equal(ax, bx) and np.array_equal(aw, bw)
    assert np.array_equal(ax[:, 0], [0, 25, 75, 125, 175])


def test_full_observation_is_predeclared_ineligible():
    obs, truth = C.split_arbor(tree([[0, 0, 0], [200, 0, 0]], [-1, 0]))
    m = C.score_completion(C.observation_only(obs), truth)
    assert not m["eligible"] and m["status"] == "ineligible_no_hidden_cable"
    assert m["distance_um"] is None and m["f1_200"] is None


def test_adapters_preserve_prefix_and_reject_disconnected_or_forged_forecasts():
    obs, truth = C.split_arbor(straight())
    pred = C.Completion(obs, [[600, 0, 0]], np.array([2]))
    assert C.Completion.from_dict(obs, pred.to_dict()).to_dict() == pred.to_dict()
    assert C.Observation.from_dict(obs.to_dict()).sha256 == obs.sha256
    assert C.from_arbor(obs, pred).to_dict() == pred.to_dict()
    forged = C.Observation(obs.xyz + [0, 1, 0], obs.parent)
    with pytest.raises(ValueError):
        C.score_completion(C.observation_only(forged), truth)
    with pytest.raises(ValueError):
        C.from_arbor(obs, tree(pred.xyz + [0, 1, 0], pred.parent))
    with pytest.raises(ValueError):
        C.Completion(obs, [[600, 0, 0]], np.array([-1]))
    with pytest.raises(ValueError):
        C.Completion.from_dict(obs, {**pred.to_dict(), "extra_parent": [1.5]})


def test_known_registration_discrepancy_does_not_move_measured_prefix():
    obs = C.observe(straight()); before = obs.xyz.copy()
    r = C.registration(obs, [3, 4, 0])
    assert r["distance_um"] == 5 and r["table_minus_trace_xyz_um"] == [3, 4, 0]
    assert np.array_equal(before, obs.xyz)
