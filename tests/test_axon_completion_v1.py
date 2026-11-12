"""Focused tests of axon completion v1: synthetic mask, leakage, shared decoder, gate. TRAIN data only."""
import importlib.util
import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(v, "1")

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("axon_completion_v1_under_test", ROOT / "analysis/axon_completion_v1.py")
A = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = A
spec.loader.exec_module(A)
C = A.C()
G = A.G()


def arbor(xyz, parent):
    xyz = np.asarray(xyz, float)
    types = np.full(len(xyz), 2)
    types[0] = 1
    return G.Arbor(xyz, np.asarray(parent, int), types).validate()


def line(start, direction, n, step=50.0):
    d = np.asarray(direction, float) / np.linalg.norm(direction)
    return [np.asarray(start, float) + d * step * (i + 1) for i in range(n)]


def forked(origin=(0, 0, 0), exact_vertex=False):
    """Trunk to ~250 um, fork into two branches that cross 300 um; optionally a vertex exactly at 300 um with two
    hidden children."""
    o = np.asarray(origin, float)
    xyz, parent = [o], [-1]
    for p in line(o, (1, 0, 0), 5):                       # depth 50..250
        xyz.append(p); parent.append(len(xyz) - 2)
    trunk = len(xyz) - 1
    for direction in ((1, 1, 0), (1, -1, 0)):
        prev = trunk
        for p in line(xyz[trunk], direction, 12, step=40.0):
            xyz.append(p); parent.append(prev); prev = len(xyz) - 1
    if exact_vertex:
        base = xyz[0] + np.array([0, 0, 1.0]) * 300.0     # straight branch: vertex exactly at depth 300
        prev = 0
        for p in line(xyz[0], (0, 0, 1), 6):
            xyz.append(p); parent.append(prev); prev = len(xyz) - 1
        assert np.allclose(xyz[prev], base)
        for direction in ((0, 1, 1), (0, -1, 1)):
            q = prev
            for p in line(base, direction, 5):
                xyz.append(p); parent.append(q); q = len(xyz) - 1
    return arbor(xyz, parent)


# --- synthetic mask and donor cut ------------------------------------------------------------------------------------

@pytest.mark.parametrize("exact", [False, True])
def test_donor_cut_matches_the_contract_split_exactly(exact):
    a = forked(exact_vertex=exact)
    obs, truth = C.split_arbor(a, A.BUDGET_UM)
    donor = A.DonorSuffix(a)
    assert math.isclose(donor.total_cable, truth.forest.cable(), rel_tol=1e-12, abs_tol=1e-9)
    assert math.isclose(donor.total_cable + C.CableForest(obs.xyz, obs.parent, A.obs_depth(obs)).cable(), a.cable(),
                        rel_tol=1e-12)
    fidx, F, _ = A.frontiers(obs)
    assert len(F) == donor.k == (3 if exact else 2)
    for r in donor.roots:
        assert np.min(np.linalg.norm(F - r, axis=1)) < 1e-9
    if exact:
        assert sorted(len(c.children()[0]) for c in donor.components) == [1, 1, 2]


def test_self_donor_reconstructs_the_hidden_suffix():
    a = forked(exact_vertex=True)
    obs, truth = C.split_arbor(a, A.BUDGET_UM)
    completion, diag = A.decode(obs, A.DonorSuffix(a))
    m = C.score_completion(completion, truth)
    assert m["status"] == "scored" and m["distance_um"] < 1e-6 and math.isclose(m["f1_200"], 1.0)
    assert math.isclose(m["extra_cable_um"], m["hidden_cable_um"], rel_tol=1e-9)
    assert set(diag["modes"]) <= {"attach_500", "attach_1000", "rigid"}


def test_public_observation_is_blind_to_the_hidden_suffix():
    a = forked()
    b = forked()
    b.xyz[-5:] += np.array([400.0, -300.0, 900.0])        # change only hidden distal geometry
    oa, ob = C.observe(a), C.observe(b)
    assert oa.sha256 == ob.sha256
    assert set(oa.to_dict()) == {"schema", "xyz_um", "parent", "path_budget_um"}     # no IDs, types or flags
    for fa, fb in zip(A.frontiers(oa), A.frontiers(ob)):
        assert np.array_equal(fa, fb)
    pa, pb = A.prefix_samples(oa), A.prefix_samples(ob)
    assert np.array_equal(pa[0], pb[0]) and np.array_equal(pa[1], pb[1])


# --- the shared decoder ---------------------------------------------------------------------------------------------

def shifted_donor(offset, exact=False):
    a = forked(exact_vertex=exact)
    return arbor(a.xyz + np.asarray(offset, float), a.parent)


@pytest.mark.parametrize("offset", [(30, 0, 0), (0, 120, 40), (900, -500, 300)])
def test_decoder_preserves_observation_uses_each_component_once_and_never_inflates(offset):
    q = forked()
    obs, truth = C.split_arbor(q, A.BUDGET_UM)
    donor = A.DonorSuffix(shifted_donor(offset, exact=True))       # 3 components onto 2 frontiers
    completion, diag = A.decode(obs, donor)
    n = len(obs.xyz)
    assert np.array_equal(completion.xyz[:n], obs.xyz) and np.array_equal(completion.parent[:n], obs.parent)
    fidx, _, _ = A.frontiers(obs)
    roots = {int(p) for p in completion.extra_parent if p < n}
    assert roots <= set(fidx.tolist())                              # merged into frontiers: no bridge to the root
    assert sum(diag["modes"].values()) == donor.k                   # every component exactly once
    extra = C.extra_forest(completion).cable()
    assert extra <= donor.total_cable * (1 + 1e-9) + 1e-6           # cable guard
    m = C.score_completion(completion, truth)
    assert m["valid_for_comparison"] and m["duplicate_extra_cable_um"] == 0.0
    if offset[0] == 900:
        assert diag["modes"] == {"rigid": 3} and math.isclose(extra, donor.total_cable, rel_tol=1e-12)


def test_more_frontiers_than_components_leaves_frontiers_unserved():
    q = forked(exact_vertex=True)                                   # 3 frontiers
    obs, _ = C.split_arbor(q, A.BUDGET_UM)
    straight = arbor([[0, 0, 0]] + line((0, 0, 0), (1, 0, 0), 20), [-1] + list(range(20)))
    donor = A.DonorSuffix(straight)
    completion, diag = A.decode(obs, donor)
    assert donor.k == 1 and sum(diag["modes"].values()) == 1
    assert len({int(p) for p in completion.extra_parent if p < len(obs.xyz)}) == 1


def test_compact_attachment_only_when_cable_does_not_grow():
    q = forked()
    obs, _ = C.split_arbor(q, A.BUDGET_UM)
    donor = A.DonorSuffix(shifted_donor((20, 10, 0)))
    for j, comp in enumerate(donor.components):
        out, mode = A.graft(comp, A.frontiers(obs)[1][j])
        if mode.startswith("attach"):
            assert out.cable() <= comp.cable() * (1 + 1e-12) + 1e-9
        else:
            assert math.isclose(out.cable(), comp.cable(), rel_tol=1e-12)


def test_decoder_refuses_truth_and_non_public_inputs():
    a = forked()
    obs, truth = C.split_arbor(a, A.BUDGET_UM)
    with pytest.raises(A.LeakError):
        A.decode(obs, truth)
    with pytest.raises(A.LeakError):
        A.decode(a, A.DonorSuffix(a))


def test_empty_forecast_is_a_failure_that_stays_in_the_mean():
    a = forked()
    obs, truth = C.split_arbor(a, A.BUDGET_UM)
    m = C.metrics_to_json(C.score_completion(C.observation_only(obs), truth))
    assert m["status"] == "failed_empty_forecast" and not m["valid_for_comparison"]

    class L:
        animal = np.array(["a", "a", "b"])
    rows = {0: {"distance_um": 100.0, "f1_200": 0.5, "valid": True}, 1: {"distance_um": None, "f1_200": 0.0, "valid": False},
            2: {"distance_um": 50.0, "f1_200": 0.4, "valid": True}}
    means = A.sample_means(L, rows)
    assert means["a"][0] == math.inf and means["a"][1] == 0.25 and means["b"] == (50.0, 0.4)


def test_invalid_but_finite_forecast_contributes_inf_and_zero():
    """Regression: a duplicate-cable forecast has finite diagnostic D/F but must count as +inf / 0."""
    class L:
        animal = np.array(["a", "a"])
    rows = {0: {"distance_um": 100.0, "f1_200": 0.5, "valid": True},
            1: {"distance_um": 80.0, "f1_200": 0.9, "valid": False, "status": "failed_duplicate_extra_cable"}}
    means = A.sample_means(L, rows)
    assert means["a"] == (math.inf, 0.25)


# --- choice rule and gate -------------------------------------------------------------------------------------------

def test_residual_choice_switches_only_when_predicted_better_on_both_heads():
    cand = np.array([3, 5, 9])
    assert A.residual_choice(cand, 5, [0.2, 0.0, 0.1], [0.01, 0.0, 0.05]) == (3, "switched")
    assert A.residual_choice(cand, 5, [0.05, 0.0, 0.1], [0.01, 0.0, 0.05]) == (9, "switched")
    assert A.residual_choice(cand, 5, [0.2, 0.0, -0.1], [-0.01, 0.0, 0.05]) == (5, "fallback")
    assert A.residual_choice(cand, 5, [np.nan, 0, 0], [0, 0, 0]) == (5, "non_finite_prediction")


def test_root_duplicate_fixture_merges_identical_children_of_the_same_parent():
    """ROOT's fixture: two parallel donor components relocated onto one frontier must not emit the branch twice."""
    import json
    fx = json.loads((ROOT / "data/interim/axon_completion_v1_root/duplicate_decoder_fixture.json").read_text())
    q = arbor(fx["query_xyz"], fx["query_parent"])
    d = arbor(fx["donor_xyz"], fx["donor_parent"])
    obs, truth = C.split_arbor(q, A.BUDGET_UM)
    donor = A.DonorSuffix(d)
    completion, diag = A.decode(obs, donor)
    n = len(obs.xyz)
    assert np.array_equal(completion.xyz[:n], obs.xyz) and np.array_equal(completion.parent[:n], obs.parent)
    m = C.score_completion(completion, truth)
    assert donor.k == 2 and diag["modes"] == {"rigid": 2}
    assert m["duplicate_extra_cable_um"] == 0.0 and m["valid_for_comparison"]
    assert math.isclose(m["extra_cable_um"], 300.0) and math.isclose(diag["merged_duplicate_cable_um"], 300.0)
    assert math.isclose(m["extra_cable_um"] + diag["merged_duplicate_cable_um"], donor.total_cable)


def test_merge_never_joins_different_parents():
    """Two components with identical child offsets on DIFFERENT frontiers stay separate (no rejoining)."""
    q = forked()
    obs, truth = C.split_arbor(q, A.BUDGET_UM)
    completion, diag = A.decode(obs, A.DonorSuffix(shifted_donor((900, -500, 300))))
    assert diag["merged_duplicate_cable_um"] == 0.0
    assert math.isclose(C.extra_forest(completion).cable(), A.DonorSuffix(shifted_donor((900, -500, 300))).total_cable,
                        rel_tol=1e-12)


def _gate_inputs(delta_d=-20.0, delta_f=0.02, ci=True):
    summary = {m: {"distance_um": 600.0, "f1_200": 0.3} for m in A.METHODS}
    versus = {c: {"distance_delta_um": delta_d, "f1_delta": delta_f, "nonfinite": False,
                  "distance_bootstrap95_um": [delta_d - 5, delta_d + 5 if not ci else -1.0],
                  "f1_bootstrap95": [0.001 if ci else -0.001, delta_f + 0.01]} for c in A.CONTROLS + A.DIAGNOSTIC}
    validity = {m: {"invalid": 0, "inflated_share": 0.1, "mean_abs_rel_cable_error": 0.5} for m in A.METHODS}
    return summary, versus, validity


def test_gate_requires_resolved_practical_margins_over_every_control():
    s, v, val = _gate_inputs()
    assert A.verdict(s, v, val, 0.3)["pass"]
    s, v, val = _gate_inputs(delta_d=-5.0)                       # below the 10 um practical margin
    assert not A.verdict(s, v, val, 0.3)["pass"]
    s, v, val = _gate_inputs(ci=False)                           # intervals cross zero
    assert not A.verdict(s, v, val, 0.3)["pass"]
    s, v, val = _gate_inputs()
    v["prefix_abs_nearest"] = dict(v["prefix_abs_nearest"], f1_delta=0.004)   # loses to the strong prefix control
    out = A.verdict(s, v, val, 0.3)
    assert not out["pass"] and not out["margins"]["prefix_abs_nearest"]["beaten"]
    s, v, val = _gate_inputs()
    val[A.SELECTOR] = {"invalid": 1, "inflated_share": 0.1, "mean_abs_rel_cable_error": 0.5}          # one invalid forecast fails the aggregate
    assert not A.verdict(s, v, val, 0.3)["pass"]
    s, v, val = _gate_inputs()
    val[A.SELECTOR] = {"invalid": 0, "inflated_share": 0.2, "mean_abs_rel_cable_error": 0.5}   # inflation guard
    assert not A.verdict(s, v, val, 0.3)["pass"]
    s, v, val = _gate_inputs()
    val[A.SELECTOR] = {"invalid": 0, "inflated_share": 0.1, "mean_abs_rel_cable_error": 0.51}  # more cable error
    out = A.verdict(s, v, val, 0.3)
    assert not out["pass"] and any("relative cable error" in r for r in out["reasons"])
    s, v, val = _gate_inputs()
    val["observation_only"] = {"invalid": 805, "inflated_share": 0.0, "mean_abs_rel_cable_error": 1.0}
    v["observation_only"] = {"distance_delta_um": None, "f1_delta": 0.3, "nonfinite": True,
                             "distance_bootstrap95_um": None, "f1_bootstrap95": None}
    assert A.verdict(s, v, val, 0.3)["pass"]                         # the diagnostic failure control is not required
    s, v, val = _gate_inputs()
    assert not A.verdict(s, v, val, 0.0)["pass"]                     # never switched


def test_one_invalid_strong_control_is_inconclusive_not_a_win():
    """Regression: a single invalid forecast must not waive a strong control whose F1 beats the selector."""
    s, v, val = _gate_inputs()
    val["prefix_abs_nearest"] = {"invalid": 1, "inflated_share": 0.1, "mean_abs_rel_cable_error": 0.5}
    v["prefix_abs_nearest"] = {"distance_delta_um": None, "f1_delta": -0.03, "nonfinite": True,
                               "distance_bootstrap95_um": None, "f1_bootstrap95": None}
    out = A.verdict(s, v, val, 0.3)
    assert not out["pass"]
    assert out["margins"]["prefix_abs_nearest"]["unevaluable"] and not out["margins"]["prefix_abs_nearest"]["beaten"]
    s, v, val = _gate_inputs()                                       # invalid but finite comparison: still fails
    val["E2_local_kernel"] = {"invalid": 1, "inflated_share": 0.1, "mean_abs_rel_cable_error": 0.5}
    assert not A.verdict(s, v, val, 0.3)["pass"]


def test_folds_are_whole_animals():
    animals = [f"a{i:03d}" for i in range(101)]
    plan = A.fold_plan(animals)
    ev = [set(f["eval_animals"]) for f in plan]
    assert set().union(*ev) == set(animals) and sum(map(len, ev)) == 101
    for f in plan:
        inner = [set(i["eval_animals"]) for i in f["inner"]]
        assert set().union(*inner) == set(f["fit_animals"]) and sum(map(len, inner)) == len(f["fit_animals"])


# --- real TRAIN arbors (read-only) -------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def lib():
    return A.V2().Library()


def test_real_donor_cut_parity_and_decode(lib):
    for nid in lib.ids[::400]:
        a = lib.load(nid)
        obs, truth = C.split_arbor(a, A.BUDGET_UM)
        donor = A.DonorSuffix(a)
        assert math.isclose(donor.total_cable, truth.forest.cable(), rel_tol=1e-9)
        fidx, F, _ = A.frontiers(obs)
        assert all(np.min(np.linalg.norm(F - r, axis=1)) < 1e-6 for r in donor.roots)
    q, d = lib.ids[10], next(n for n in lib.ids if lib.animal[lib.index[n]] != lib.animal[10])
    obs, truth = C.split_arbor(lib.load(q), A.BUDGET_UM)
    completion, diag = A.decode(obs, A.DonorSuffix(lib.load(d)))
    m = C.score_completion(completion, truth)
    assert m["valid_for_comparison"] and m["extra_cable_um"] <= diag["donor_hidden_cable_um"] * (1 + 1e-9) + 1e-6


def test_label_door_translates_the_real_frozen_v2_leak_class(lib):
    animals = sorted(set(lib.animal.tolist()))
    door = A.LabelDoor("t", lib, animals[:50], A.Labels(lib), A.V2().Cache(lib))
    inside = int(np.flatnonzero(np.isin(lib.animal, animals[:50]))[0])
    outside = int(np.flatnonzero(~np.isin(lib.animal, animals[:50]))[0])
    assert A.V2().LeakError is not A.LeakError
    with pytest.raises(A.LeakError) as info:
        door.label(outside, inside)
    assert isinstance(info.value.__cause__, A.V2().LeakError)
    with pytest.raises(A.LeakError):
        door.D([inside], [outside])
    try:                                                     # the probe pattern used in the run
        door.label(outside, inside)
    except A.leak_types():
        pass
    else:
        raise AssertionError("probe served")
    raw = A.V2().FitDoor("raw", lib, animals[:50])            # the raw foreign class is caught by leak_types too
    with pytest.raises(A.leak_types()):
        raw.label(outside, inside)


def test_label_job_refuses_a_changed_observation_and_same_animal(lib):
    A._worker_init(False)
    q = lib.ids[20]
    with pytest.raises(RuntimeError):
        A.label_job({"job": "t", "query": q, "obs_sha": "0" * 64, "donors": []})
    obs = C.observe(lib.load(q))
    same = next(n for n in lib.ids if n != q and lib.animal[lib.index[n]] == lib.animal[lib.index[q]])
    with pytest.raises(A.LeakError):
        A.label_job({"job": "t", "query": q, "obs_sha": obs.sha256, "donors": [same]})


def test_metadata_features_have_no_observation_input_and_no_identity():
    names_meta = ["q_soma_ap", "diff_ap", "soma_distance", "donor_log_cable", "w_E1", "is_fallback"]
    A.V2().check_feature_names(names_meta)
    assert not any(n.startswith(A.OBSERVATION_FEATURE_PREFIXES) for n in names_meta)
    with pytest.raises(Exception):
        A.V2().check_feature_names(["query_id"])
