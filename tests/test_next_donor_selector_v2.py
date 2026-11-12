"""Focused tests of the next donor selector v2 (leakage, fallback, numerical parity). TRAIN data only."""
import importlib.util
import math
import os
from pathlib import Path

import numpy as np
import pytest

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(v, "1")

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("nds_v2", ROOT / "analysis/next_donor_selector_v2.py")
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)


@pytest.fixture(scope="module")
def lib():
    return S.Library()


@pytest.fixture(scope="module")
def cache(lib):
    return S.Cache(lib)


@pytest.fixture(scope="module")
def plan(lib):
    return S.fold_plan(lib)


# --- combined rule and fallback -------------------------------------------------------------------------------------

def test_combined_rule_switches_only_when_better_on_both_heads():
    pd_ = np.array([0.5, 0.4, 0.6])
    pf = np.array([0.2, 0.3, 0.1])
    assert S.combined_choice(pd_, pf, 0) == (1, "switched")
    pf2 = np.array([0.35, 0.3, 0.1])                       # top (1) predicted worse F1 than the fallback (0)
    assert S.combined_choice(pd_ - np.array([0, .2, 0]), pf2, 0)[1] in ("top_is_fallback", "not_better_on_both_heads")
    assert S.combined_choice(np.array([0.1, 0.4]), np.array([0.5, 0.3]), 0) == (0, "top_is_fallback")
    assert S.combined_choice(np.array([np.nan, 0.1]), np.array([0.2, 0.9]), 0) == (0, "non_finite_prediction")


def test_fallback_only_or_copying_only_result_is_never_a_win():
    summary = {m: {"distance_um": 600.0, "f1_200": 0.28} for m in (S.SELECTOR,) + S.CONTROLS}
    versus = {c: {"distance_bootstrap95_um": [0, 0], "f1_bootstrap95": [0, 0]} for c in S.CONTROLS}
    v = S.verdict(summary, versus, switch_rate=0.0)
    assert not v["gate_met"] and not v["beats_E1_E2_and_adaptive_fallback_on_both"]
    better = dict(summary)
    better[S.SELECTOR] = {"distance_um": 500.0, "f1_200": 0.4}
    assert not S.verdict(better, versus, switch_rate=0.0)["gate_met"]          # no switch: not a learned win
    assert S.verdict(better, versus, switch_rate=0.1)["gate_met"]


def test_adaptive_fallback_win_with_harmful_switches_fails():
    """Adaptive fallback beats fixed E1/E2; the learned switches make the selector worse than that fallback."""
    summary = {m: {"distance_um": 650.0, "f1_200": 0.26} for m in (S.SELECTOR,) + S.CONTROLS}
    summary["E1_local_uniform"] = {"distance_um": 640.0, "f1_200": 0.27}
    summary["E2_local_kernel"] = {"distance_um": 638.0, "f1_200": 0.271}
    summary["adaptive_fallback"] = {"distance_um": 600.0, "f1_200": 0.30}
    summary[S.SELECTOR] = {"distance_um": 610.0, "f1_200": 0.29}          # beats E1/E2 and every copy, not the fallback
    versus = {c: {"distance_bootstrap95_um": [-1, -0.5], "f1_bootstrap95": [0.001, 0.01]} for c in S.CONTROLS}
    v = S.verdict(summary, versus, switch_rate=0.2)
    assert v["beats"]["E1_local_uniform"] and v["beats"]["E2_local_kernel"]
    assert not v["beats"]["adaptive_fallback"]
    assert not v["beats_E1_E2_and_adaptive_fallback_on_both"] and not v["gate_met"]
    summary[S.SELECTOR] = {"distance_um": 590.0, "f1_200": 0.31}
    assert S.verdict(summary, versus, switch_rate=0.2)["gate_met"]


def test_fallback_emits_exact_e1_arrays(lib):
    q = lib.index[lib.ids[10]]
    mask = lib.animal != lib.animal[q]
    door = S.FitDoor("t", lib, sorted(set(lib.animal.tolist())), cache=S.Cache(lib))
    loc = S.local_choices(lib, door, q, mask)
    d = lib.ids[int(loc["cand"][loc["choice"]["E1_local_uniform"]])]
    a1, e1 = S.attach(lib.load(d), lib.soma[q])
    a2, e2 = S.emit(lib.load(d), lib.soma[q], "attached")
    assert e1 == e2 and S.arbor_sha(a1) == S.arbor_sha(a2)


def test_attachment_preserves_distal_coordinates_and_translation_is_exact(lib):
    d = lib.load(lib.ids[0])
    near = d.xyz[0] + np.array([30.0, 0, 0])
    out, how = S.attach(d, near)
    assert how == "attach_500"
    depth = d.path_lengths()
    distal = d.xyz[depth >= 500.0]
    assert all(any(np.array_equal(p, o) for o in out.xyz) for p in distal[:50])
    far = d.xyz[0] + np.array([800.0, 0, 0])
    out, how = S.attach(d, far)
    assert how == "translation" and np.array_equal(out.xyz, d.shifted(far).xyz)


# --- leakage --------------------------------------------------------------------------------------------------------

def test_fit_door_refuses_outside_and_same_animal(lib, cache):
    animals = sorted(set(lib.animal.tolist()))
    fit = animals[:50]
    store = S.LabelStore(lib)
    door = S.FitDoor("t", lib, fit, labels=store, cache=cache)
    inside = [i for i in range(lib.n) if lib.animal[i] in fit]
    outside = [i for i in range(lib.n) if lib.animal[i] not in fit]
    a, b = inside[0], next(i for i in inside if lib.animal[i] != lib.animal[inside[0]])
    store.add(a, b, "attached", (1.0, 0.5, "translation", "x"))
    assert door.label(a, b)[0] == 1.0
    with pytest.raises(S.LeakError):
        door.label(outside[0], b)
    with pytest.raises(S.LeakError):
        door.label(a, outside[0])
    same = next(i for i in inside if i != a and lib.animal[i] == lib.animal[a])
    with pytest.raises(S.LeakError):
        door.label(a, same)
    with pytest.raises(S.LeakError):
        door.D([a], [outside[0]])
    assert sum(door.refused.values()) == 4


def test_query_animal_never_in_its_own_pool(lib, cache):
    q = 5
    with pytest.raises(S.LeakError):
        S.local_candidates(lib, q, np.ones(lib.n, bool))


def test_folds_are_whole_animals_and_inner_folds_stay_in_outer_fit(plan, lib):
    animals = set(lib.animal.tolist())
    evals = [set(f["eval_animals"]) for f in plan]
    assert set().union(*evals) == animals and sum(len(e) for e in evals) == len(animals)
    for f in plan:
        assert set(f["fit_animals"]) == animals - set(f["eval_animals"])
        ie = [set(i["eval_animals"]) for i in f["inner"]]
        assert set().union(*ie) == set(f["fit_animals"]) and sum(len(e) for e in ie) == len(f["fit_animals"])
        for i in f["inner"]:
            assert set(i["fit_animals"]) | set(i["eval_animals"]) == set(f["fit_animals"])


def test_every_call_draws_candidates_only_from_its_fit_animals(lib, cache, plan):
    fold = plan[0]
    inner = fold["inner"][0]
    for name, pool_animals, queries in (
            ("outer_eval", fold["fit_animals"], S.neurons_of(lib, fold["eval_animals"])[:60]),
            ("inner_eval", inner["fit_animals"], S.neurons_of(lib, inner["eval_animals"])[:60]),
            ("inner_fit", inner["fit_animals"], S.neurons_of(lib, inner["fit_animals"])[:60])):
        locs, _ = S.call_locals(lib, cache, name, pool_animals, queries)
        for q, loc in locs.items():
            assert len(loc["cand"]) == S.K
            assert set(lib.animal[loc["cand"]].tolist()) <= set(pool_animals)
            assert lib.animal[q] not in set(lib.animal[loc["cand"]].tolist())


def test_training_rows_refuse_a_held_out_query_and_features_have_no_identity(lib, cache, plan):
    fold = plan[0]
    inner = fold["inner"][0]
    store = S.LabelStore(lib)
    door = S.FitDoor("t", lib, inner["fit_animals"], labels=store, cache=cache)
    feat = S.Featurizer(lib, door)
    ev = S.neurons_of(lib, inner["eval_animals"])[:3]
    locs, _ = S.call_locals(lib, cache, "x", inner["fit_animals"], ev)
    with pytest.raises(S.LeakError):
        S.training_rows(lib, door, feat, locs)               # eval query labels are never served to a fit
    X = feat.rows(locs[ev[0]], lib.known(ev[0]))              # forecasting features need known inputs only
    assert X.shape[0] == S.K and np.isfinite(X).all()
    S.check_feature_names(feat.names)
    assert not any(n.lower().startswith(("id", "animal", "neuron")) for n in feat.names)
    assert set(lib.known(ev[0])) == set(S.QUERY_KNOWN_FIELDS)
    # a candidate outside the fit animals has no centrality: refused
    other = S.FitDoor("o", lib, fold["eval_animals"], cache=cache)
    q_fit = S.neurons_of(lib, inner["fit_animals"])[0]
    loc_bad = S.local_choices(lib, other, q_fit, np.isin(lib.animal, fold["eval_animals"]) & (lib.animal != lib.animal[q_fit]))
    with pytest.raises(S.LeakError):
        feat.rows(loc_bad, lib.known(q_fit))


def test_check_feature_names_refuses_identity_columns():
    with pytest.raises(S.LeakError):
        S.check_feature_names(["soma_distance", "animal_code"])
    with pytest.raises(S.LeakError):
        S.check_feature_names(["query_id"])


# --- numerical parity -----------------------------------------------------------------------------------------------

def test_parity_with_the_engines_development_forecasts(lib, cache):
    out = S.parity(lib, cache, n_arrays=4, log=lambda m: None)
    assert out["all_equal"], out["mismatches"]
    assert out["E1_donor_equal"] == out["E2_donor_equal"] == out["queries"] == 919


def test_exact_label_is_the_frozen_metric_of_the_emitted_geometry_at_the_table_soma(lib):
    S._worker_init()
    q = lib.ids[100]
    d = next(n for n in lib.ids if lib.animal[lib.index[n]] != lib.animal[lib.index[q]])
    res = S.label_job({"job": "t", "query": q, "pairs": [[d, "attached"], [d, "translated"]]})
    truth = lib.load(q)
    for (did, mode, dist, f1, executed, h) in res["pairs"]:
        arbor, ex = S.emit(lib.load(did), lib.soma[lib.index[q]], mode)
        m = S.template().G.geometry_metrics(arbor, truth)
        assert (dist, f1, executed, h) == (m["distance_um"], m["f1_200"], ex, S.arbor_sha(arbor))


def test_bootstrap_and_objective_are_deterministic():
    per_a = {"a": (100.0, 0.3), "b": (200.0, 0.2), "c": (300.0, 0.4)}
    per_b = {"a": (110.0, 0.25), "b": (210.0, 0.2), "c": (290.0, 0.3)}
    x = S.paired_bootstrap(per_a, per_b, ["a", "b", "c"])
    assert x == S.paired_bootstrap(per_a, per_b, ["a", "b", "c"])
    assert math.isclose(x["distance_delta_um"], -10 / 3)
    assert S.objective([("a", 200.0, 0.5)]) == pytest.approx(math.log(2) + 0.5)
