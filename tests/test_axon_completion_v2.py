"""Focused synthetic checks only; no biological fitting, arbors or score file reads."""
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("axon_completion_v2_test", ROOT / "analysis/axon_completion_v2.py")
B = importlib.util.module_from_spec(spec); sys.modules[spec.name] = B; spec.loader.exec_module(B)
A = B.V1()


def synthetic_context(n=80):
    rng = np.random.default_rng(37)
    lib = SimpleNamespace(n=n, ids=[f"n{i:03d}" for i in range(n)],
                          animal=np.array([f"a{i // 2:03d}" for i in range(n)]),
                          soma=rng.normal(size=(n, 3)) * 100,
                          peptide=np.array(["p0"] * n), nucleus=np.array(["N"] * n), hemisphere=np.array(["L"] * n))
    lib.index = {s: i for i, s in enumerate(lib.ids)}
    lib.known = lambda q: {"soma": lib.soma[q].copy(), "peptide": lib.peptide[q],
                           "nucleus": lib.nucleus[q], "hemisphere": lib.hemisphere[q]}
    D = np.linalg.norm(lib.soma[:, None] - lib.soma[None, :], axis=2)
    cache = SimpleNamespace(D=D)
    prefix_cost = D.astype(np.float32).copy()  # Public prefix geometry is independent of hidden full-D.
    idx = SimpleNamespace(eligible=np.ones(n, bool), root=lib.soma.copy(),
                          F=[s[None] + [300, 0, 0] for s in lib.soma],
                          V=[np.array([[1., 0, 0]]) for _ in range(n)],
                          donor_R=[s[None] + [300, 0, 0] for s in lib.soma],
                          donor_U=[np.array([[1., 0, 0]]) for _ in range(n)],
                          obs_cable=np.full(n, 300.), obs_vertices=np.full(n, 7.), registration=np.zeros(n),
                          donor_hidden_cable=np.linspace(1000, 3000, n), donor_k=np.ones(n),
                          summary={"log_cable": np.linspace(3, 6, n)})
    idx.prefix_costs = lambda: (prefix_cost, prefix_cost)
    return lib, cache, idx


def door(lib, cache, fit=None, labels=None):
    fit = sorted(set(lib.animal[:68].tolist())) if fit is None else fit
    return A.LabelDoor("synthetic", lib, fit, labels, cache)


def test_landmarks_are_fit_only_deterministic_and_one_per_sample():
    lib, cache, _ = synthetic_context(); d = door(lib, cache)
    lm = B.GeometryLandmarks(lib, d)
    assert len(lm.indices) == 32 and len(set(lib.animal[lm.indices])) == 32
    assert all(d.mask[lm.indices]) and lm.ids == B.GeometryLandmarks(lib, door(lib, cache)).ids
    fit = np.flatnonzero(d.mask)
    count = {a: np.sum(lib.animal[fit] == a) for a in set(lib.animal[fit])}
    w = np.array([1 / count[a] for a in lib.animal[fit]]); w /= w.sum()
    risk = cache.D[np.ix_(fit, fit)] @ w
    first = min(range(len(fit)), key=lambda i: (risk[i], lib.ids[fit[i]]))
    assert lm.indices[0] == fit[first]
    for j in range(1, len(lm.indices)):
        used = set(lib.animal[lm.indices[:j]])
        possible = [i for i in fit if lib.animal[i] not in used]
        expected = min(possible, key=lambda i: (-cache.D[i, lm.indices[:j]].min(), lib.ids[i]))
        assert lm.indices[j] == expected


def test_landmark_ties_use_id_and_cap_only_by_sample_count():
    lib, cache, _ = synthetic_context(12); cache.D[:] = 0
    lm = B.GeometryLandmarks(lib, door(lib, cache, sorted(set(lib.animal))))
    assert list(lm.indices) == [0, 2, 4, 6, 8, 10]


def test_heldout_landmark_or_donor_reads_use_interoperable_leak_class():
    lib, cache, _ = synthetic_context(); d = door(lib, cache); lm = B.GeometryLandmarks(lib, d)
    with pytest.raises(A.LeakError):
        lm.rows(np.array([70]))
    with pytest.raises(A.LeakError):
        d.D([0], [70])
    raw = A.V2().FitDoor("raw", lib, sorted(set(lib.animal[:68])), cache=cache)
    with pytest.raises(A.leak_types()):
        raw.D([0], [70])


def test_both_columns_for_the_query_landmark_sample_are_masked():
    lib, cache, idx = synthetic_context(); d = door(lib, cache); f = B.Featurizer(lib, idx, d)
    q = int(f.landmarks.indices[0]); mask = d.mask & (lib.animal != lib.animal[q])
    loc = A.Local(lib, idx, d.door, q, mask); x = f.rows(loc, loc.choice["E1_local_uniform"])
    base, count = len(f.base.names), len(f.landmarks.indices)
    masked = lib.animal[f.landmarks.indices] == lib.animal[q]
    assert np.isnan(x[:, base:base + count][:, masked]).all()
    assert np.isnan(x[:, base + count:][:, masked]).all()
    assert np.isfinite(x[:, base:base + count][:, ~masked]).all()
    assert np.isfinite(x[:, base + count:][:, ~masked]).all()


@pytest.fixture(scope="module")
def toy_model():
    rng = np.random.default_rng(20261113)
    n = 2500; route = (rng.random(n) > .65).astype(int); donor = rng.integers(0, 2, n)
    hazard = (rng.random(n) < .35).astype(int); fallback = (rng.random(n) < .20).astype(int)
    qdir = np.where(route, 1., -1.) + rng.uniform(-.15, .15, n)
    X = np.column_stack([qdir, donor, 1 - donor, hazard, fallback, np.full(n, 300.)])
    match = route == donor
    gd = np.where(match, .2 + .4 * hazard, -.2)
    gf = np.where(match, .08 + .1 * hazard, -.08)
    gc = np.where(hazard, -.5, .2)
    gd[fallback == 1] = 0; gf[fallback == 1] = 0; gc[fallback == 1] = 0
    return B.ResidualHeads().fit(X, gd, gf, gc, np.ones(n))


def test_fixed_hgb_recovers12_new_toy_routes_and_vetoes_attractive_cable_harm(toy_model):
    correct = 0
    for qdir in [-1.13, -1.08, -1.02, -.98, -.92, -.87, .87, .92, .98, 1.02, 1.08, 1.13]:
        # Known prefix direction; donor landmark coordinates; donor cable-size marker; known fallback flag.
        X = np.array([[qdir, 0, 1, 0, 1, 300], [qdir, 0, 1, 0, 0, 300],
                      [qdir, 1, 0, 0, 0, 300], [qdir, 0, 1, 1, 0, 300], [qdir, 1, 0, 1, 0, 300]])
        gd, gf, gc = toy_model.predict(X)
        dangerous = 3 if qdir < 0 else 4
        assert gc[dangerous] < 0 and gd[dangerous] + gf[dangerous] > gd[1 if qdir < 0 else 2] + gf[1 if qdir < 0 else 2]
        selected, why = B.residual_choice(np.arange(5), 0, gd, gf, gc)
        assert why == "switched" and selected == (1 if qdir < 0 else 2)
        correct += 1
    assert correct == 12


def test_hidden_outer_geometry_cache_rows_do_not_change_landmarks_features_or_predictions(toy_model):
    lib, cache, idx = synthetic_context(); d = door(lib, cache)
    loc = A.Local(lib, idx, d.door, 70, d.mask.copy())
    fb = loc.choice["E1_local_uniform"]; f = B.Featurizer(lib, idx, d); x = f.rows(loc, fb)
    ids = f.landmarks.ids
    cache.D[~d.mask, :] = 1e9; cache.D[:, ~d.mask] = 1e9
    for values in idx.summary.values():
        values[~d.mask] = -1e12
    idx.donor_hidden_cable[~d.mask] = 1e12
    for i in np.flatnonzero(~d.mask):
        idx.donor_R[i] = np.array([[9e8, 9e8, 9e8]])
    # Public observation and known metadata stay exactly fixed.
    f2 = B.Featurizer(lib, idx, door(lib, cache)); x2 = f2.rows(loc, fb)
    assert ids == f2.landmarks.ids and np.array_equal(x, x2)
    assert x.shape[1] == len(f.base.names) + 64
    for a, b in zip(toy_model.predict(x[:, :6]), toy_model.predict(x2[:, :6])):
        assert np.array_equal(a, b)


def test_third_head_sign_strict_switches_and_exact_fallback():
    cand = np.array([7, 3, 5])
    assert B.residual_choice(cand, 5, [.4, .4, 0], [.1, .1, 0], [0, 0, 0]) == (3, "switched")
    assert B.residual_choice(cand, 5, [.4, -.2, 0], [.1, .1, 0], [-.01, .2, 0]) == (5, "fallback")
    assert B.residual_choice(cand, 5, [.4, -.2, 0], [.1, .1, 0], [np.nan, .2, 0]) == (5, "non_finite_prediction")


def test_new_random_and_permuted_controls_are_not_old_switch_sets():
    lib, _, _ = synthetic_context(8)
    loc = SimpleNamespace(q=6, cand=np.array([0, 1, 2, 4]))
    old = {"donors": {B.SELECTOR: lib.ids[2], "adaptive_fallback": lib.ids[1], "switch_random": lib.ids[4],
                      "permuted_labels": lib.ids[4], "metadata_only": lib.ids[0]}}
    controls = B.make_controls(lib, loc, 1, 1, "fallback", 0, old)
    assert controls["switch_random"] == lib.ids[1] and controls["permuted_labels"] == lib.ids[0]
    assert controls[B.V1_CONTROL] == old["donors"][B.SELECTOR] and controls["metadata_only"] == lib.ids[0]
    switched = B.make_controls(lib, loc, 1, 0, "switched", 2, old)
    assert switched["switch_random"] != lib.ids[1]


def test_frozen_local_parity_is_checked_but_full_copies_may_be_outside_union():
    lib, cache, idx = synthetic_context(); d = door(lib, cache)
    loc = A.Local(lib, idx, d.door, 70, d.mask.copy())
    old = {"donors": {m: lib.ids[v] for m, v in loc.choice.items()}}
    old["donors"].update({m: lib.ids[int(loc.cand[0])] for m in
                          (B.SELECTOR, "metadata_only", "switch_random", "permuted_labels")})
    old["donors"]["medoid"] = lib.ids[79]
    B.assert_local_parity(lib, loc, old)
    old["donors"]["prefix_abs_nearest"] = lib.ids[79]
    with pytest.raises(RuntimeError, match="local control"):
        B.assert_local_parity(lib, loc, old)


def test_training_targets_weights_and_invalid_rows_match_v1():
    lib, cache, _ = synthetic_context(6)
    labels = A.Labels(lib)
    for q, donor, distance, f1, extra, valid in [(0, 2, 400, .2, 200, True), (0, 4, 200, .4, 100, False),
                                                 (2, 0, 300, .3, 100, True), (2, 4, 200, .5, 150, True)]:
        labels.add(lib.ids[q], {"donor": lib.ids[donor], "valid": valid, "distance_um": distance,
                                "f1_200": f1, "hidden_cable_um": 100., "extra_cable_um": extra})
    d = door(lib, cache, sorted(set(lib.animal)), labels)
    locals_ = {0: SimpleNamespace(cand=np.array([2, 4])), 2: SimpleNamespace(cand=np.array([0, 4]))}
    class F:
        def rows(self, loc, fb, *, candidates):
            return np.asarray(candidates)[:, None].astype(float)
    fb = {0: 2, 2: 0}
    x, gd, gf, gc, w, skip = B.training_rows(lib, d, F(), locals_, fb)
    vx, vd, vf, vw, vs = A.training_rows(lib, d, F(), locals_, fb, labels, metadata=False)
    assert np.array_equal(x, vx) and np.array_equal(gd, vd) and np.array_equal(gf, vf)
    assert np.array_equal(w, vw) and skip == vs == 1
    assert np.allclose(gc, [0, 0, -.5]) and np.allclose(w, [1.5, .75, .75])
    with pytest.raises(ValueError):
        B.cable_error({"valid": True, "hidden_cable_um": 0, "extra_cable_um": 1})


def sealed_fixture(tmp_path):
    A.write_json(tmp_path / "forecasts.json", {})
    A.write_json(tmp_path / "fold_info.json", {})
    A.write_json(tmp_path / "RUN_STARTED.json", {})
    A.write_json(tmp_path / "FORECASTS_SEALED.json", {
        "forecasts_sha256": A.sha_file(tmp_path / "forecasts.json"),
        "fold_info_sha256": A.sha_file(tmp_path / "fold_info.json"),
        "run_started_sha256": A.sha_file(tmp_path / "RUN_STARTED.json")})


def test_eval_labels_cannot_be_loaded_before_seal_and_changed_seal_blocks_resume(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Evaluation labels admitted before seal")
    monkeypatch.setattr(A, "Labels", forbidden)
    with pytest.raises(FileNotFoundError):
        B.load_evaluation_labels(None, tmp_path)
    sealed_fixture(tmp_path); B.evaluation_admitted(tmp_path)
    (tmp_path / "forecasts.json").write_text('{"changed":true}')
    with pytest.raises(RuntimeError, match="Forecast seal"):
        B.evaluation_admitted(tmp_path)
    sealed_fixture(tmp_path); (tmp_path / "fold_info.json").write_text('{"changed":true}')
    with pytest.raises(RuntimeError, match="provenance"):
        B.evaluation_admitted(tmp_path)


def test_extended_gate_requires_frozen_v1_selector_and_unchanged_cable_rule():
    methods = A.METHODS + (B.V1_CONTROL,)
    s = {m: {"distance_um": 600., "f1_200": .3} for m in methods}
    v = {m: {"distance_delta_um": -20., "f1_delta": .02, "nonfinite": False,
             "distance_bootstrap95_um": [-30, -10], "f1_bootstrap95": [.01, .03]} for m in A.CONTROLS + A.DIAGNOSTIC + (B.V1_CONTROL,)}
    val = {m: {"invalid": 0, "inflated_share": .1, "mean_abs_rel_cable_error": .5} for m in methods}
    assert B.verdict(s, v, val, .2)["pass"]
    v[B.V1_CONTROL]["f1_delta"] = .005
    assert not B.verdict(s, v, val, .2)["pass"]
    v[B.V1_CONTROL]["f1_delta"] = .02
    val[B.SELECTOR]["mean_abs_rel_cable_error"] = .5001
    assert not B.verdict(s, v, val, .2)["pass"]
