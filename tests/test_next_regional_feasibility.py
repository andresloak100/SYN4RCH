"""Tests of analysis/next_regional_core.py and analysis/next_regional_feasibility.py (synthetic data only)."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]


def load(name, rel):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


nc = load("next_regional_core", "analysis/next_regional_core.py")
nr = load("next_regional_feasibility_t", "analysis/next_regional_feasibility.py")


def brute_crps(x, y):
    x = np.asarray(x, float)
    return np.abs(x - y).mean() - 0.5 * np.abs(x[:, None] - x[None, :]).mean()


def test_crps_ensemble_matches_brute_force():
    rng = np.random.default_rng(0)
    ens = rng.gamma(1.0, 1.0, size=(7, 3, 11)) * (rng.random((7, 3, 11)) > 0.4)
    y = rng.gamma(1.0, 1.0, size=(7, 3)) * (rng.random((7, 3)) > 0.5)
    want = np.array([[brute_crps(ens[i, j], y[i, j]) for j in range(3)] for i in range(7)])
    np.testing.assert_allclose(nc.crps_ensemble(ens, y), want, atol=1e-12)


def test_quantile_crps_approximates_exact_crps():
    rng = np.random.default_rng(4)
    sample = np.sort(rng.gamma(2.0, 1.0, 5000) * (rng.random(5000) > 0.6))
    q = nc.to_quantiles(sample[None])[0]
    for y in (0.0, 0.5, 2.0, 6.0):
        assert nc.quantile_crps(q, np.array(y)) == pytest.approx(brute_crps(sample[::5], y), abs=0.01)


def test_quantile_crps_of_point_mass_is_absolute_error():
    q = np.full(nc.K, 2.0)
    assert nc.quantile_crps(q, np.array(0.5)) == pytest.approx(1.5)


def test_quantile_at_interpolates_and_is_flat_beyond_the_ends():
    q = np.arange(nc.K, dtype=float)                      # Q(level) = level * K - 0.5 inside
    assert nc.quantile_at(q, 0.5) == pytest.approx(99.0)
    assert nc.quantile_at(q, 0.0) == 0.0 and nc.quantile_at(q, 1.0) == nc.K - 1
    assert nc.quantile_at(q, 0.05) == pytest.approx(0.05 * nc.K - 0.5)


def test_conditional_nonzero_interval_is_calibrated_for_a_mixture():
    """Truths drawn from the forecast itself: marginal and conditional non-zero coverage are both near 0.90,
    while a zero-dominated forecast's marginal interval covers non-zero truths far less (why leg 3b exists)."""
    rng = np.random.default_rng(5)
    p = 0.15
    big = np.where(rng.random(200000) < p, rng.gamma(3.0, 0.8, 200000), 0.0)
    q = nc.to_quantiles(big[None])[0]
    y = np.where(rng.random(20000) < p, rng.gamma(3.0, 0.8, 20000), 0.0)
    Q = np.tile(q, (len(y), 1))[:, None, :]
    P = np.full((len(y), 1), p)
    scores, mats = nc.score_forecast(np.zeros((len(y), 1)), Q, P, y[:, None])
    assert scores["coverage"].mean() == pytest.approx(0.90, abs=0.06)       # marginal: zeros covered for free
    nzcov = np.nanmean(mats["inside_nz"])
    assert nzcov == pytest.approx(0.90, abs=0.02)
    assert np.nanmean(np.where(mats["nonzero"], mats["inside"], np.nan)) < 0.75


def test_conditional_interval_without_nonzero_mass_never_covers_a_nonzero_truth():
    q = np.zeros((1, 1, nc.K))
    scores, mats = nc.score_forecast(np.zeros((1, 1)), q, np.zeros((1, 1)), np.array([[1.0]]))
    assert scores["coverage_nonzero_conditional"][0] == 0.0


def test_animal_mean_counts_every_sample_once():
    overall, ids, values = nc.animal_mean([1, 1, 1, 1, 5], ["a", "a", "a", "a", "b"])
    assert overall == pytest.approx(3.0) and list(values) == [1.0, 5.0]


def test_minimum_detectable_shrinks_with_samples():
    a, b = nc.minimum_detectable(1.0, 25), nc.minimum_detectable(1.0, 40)
    assert b < a and b == pytest.approx(2.874 / np.sqrt(40), rel=0.01)


def synthetic(n_samples=12, per=15, seed=3, regions=4):
    rng = np.random.default_rng(seed)
    group = np.repeat([f"s{i:02d}" for i in range(n_samples)], per)
    label = np.array(["A" if int(g[1:]) % 2 else "B" for g in group])
    soma = rng.normal(0, 100, size=(len(group), 3))
    Y = np.maximum(0.0, 1.0 + 0.01 * soma[:, [0]] + (label == "A")[:, None] * 0.5
                   + rng.normal(0, 0.1, size=(len(group), regions)))
    return soma, label, group, Y


def test_controls_forecast_shapes_and_fallbacks():
    soma, label, group, Y = synthetic()
    c = nc.Controls(soma, label, group, Y)
    for name, setting in (("MEAN", {"weighting": "sample"}), ("IDENTITY", {"weighting": "neuron"}),
                          ("KNN", {"k": 10, "pool": "same_label", "scaling": "z"}),
                          ("RIDGE", {"lambda": 1.0, "degree": 2}), ("KNN-CONF", {"k": 10, "pool": "all", "scaling": "raw"})):
        point, q, p = c.forecast(name, setting, soma[:5], np.array(["Z"] * 5))
        assert point.shape == (5, 4) and q.shape == (5, 4, nc.K) and p.shape == (5, 4)
        assert (np.diff(q, axis=-1) >= -1e-12).all() and (q >= 0).all() and ((p >= 0) & (p <= 1)).all()
    np.testing.assert_allclose(c.identity_point(np.array(["Z"] * 5), "neuron"), c.mean_point(5, "neuron"))


def test_ridge_residuals_are_out_of_fold():
    """In-sample residuals of an over-flexible fit shrink; out-of-fold ones do not."""
    soma, label, group, Y = synthetic()
    c = nc.Controls(soma, label, group, Y)
    oof = c.ridge_residuals(1e-3, 2)
    model = nc._ridge_solve(c._design(soma, label, 2), Y, 1e-3, len(c.labels))
    insample = Y - nc._ridge_apply(model, c._design(soma, label, 2))
    assert np.abs(oof).mean() > np.abs(insample).mean()
    assert set(np.unique(c.subfold)) == {0, 1, 2}
    for s in set(group):                                     # whole samples per residual fold
        assert len(set(c.subfold[group == s])) == 1


def test_same_label_pool_with_too_few_members_uses_those_available():
    soma, label, group, Y = synthetic()
    label = label.astype("<U8")
    label[0] = "rare"
    c = nc.Controls(soma, label, group, Y)
    order, valid = c.neighbours(soma[:1], np.array(["rare"]), "raw", "same_label", 5)
    assert valid[0].sum() == 1 and order[0, 0] == 0
    point, q, p = c.forecast("KNN", {"k": 5, "pool": "same_label", "scaling": "raw"}, soma[:1], np.array(["rare"]))
    np.testing.assert_allclose(point, Y[:1])


def test_choose_settings_uses_whole_sample_folds_and_returns_every_control():
    soma, label, group, Y = synthetic()
    folds = {"0": ["s00", "s01", "s02", "s03"], "1": ["s04", "s05", "s06", "s07"], "2": ["s08", "s09", "s10", "s11"]}
    chosen = nc.choose_settings(soma, label, group, Y, folds)
    assert set(chosen) == set(nc.TUNED) | {"KNN-CONF"}
    assert chosen["KNN-CONF"]["setting"] == chosen["KNN"]["setting"]
    assert chosen["RIDGE@crps"]["grid_size"] == 14 and chosen["KNN@crps"]["criterion"] == "crps"


def test_decide_needs_every_leg():
    animal = {k: {"crps": 0.50, "log_mse": 1.30} for k in set(nc.POINT_CONTROLS) | set(nc.DIST_CONTROLS)}
    animal["KNN"]["log_mse"] = 1.25
    good = {"crps": 0.48, "log_mse": 1.26, "coverage": 0.90, "coverage_nonzero_conditional": 0.89}
    region_ok = {"regions_below_floor": 2, "regions_presence_off": 0}
    assert nc.decide({**animal, "CANDIDATE": good}, region_ok)["outcome"] == "pass"
    assert nc.decide({**animal, "CANDIDATE": {**good, "crps": 0.495}}, region_ok)["outcome"] == "no_gain_over_the_controls"
    assert nc.decide({**animal, "CANDIDATE": {**good, "log_mse": 1.27}}, region_ok)["outcome"] == "better_crps_but_worse_points"
    assert nc.decide({**animal, "CANDIDATE": {**good, "coverage_nonzero_conditional": 0.80}},
                     region_ok)["outcome"] == "better_crps_without_valid_calibration"
    assert nc.decide({**animal, "CANDIDATE": good}, {**region_ok, "regions_below_floor": 6})["outcome"] != "pass"
    assert nc.decide({**animal, "CANDIDATE": good}, {**region_ok, "regions_presence_off": 6})["outcome"] != "pass"


def test_region_checks_count_only_eligible_regions():
    n = 100
    nz = np.zeros((n, 2), bool)
    nz[:40, 0] = True                                         # region 0: 40 non-zero truths in 8 samples
    nz[:10, 1] = True                                         # region 1: too few non-zero truths
    groups = np.repeat([f"g{i}" for i in range(20)], 5)
    inside_nz = np.where(nz, 0.0, np.nan)                     # nothing covered
    mats = {"nonzero": nz, "inside_nz": inside_nz, "inside": np.ones((n, 2), bool), "p_nonzero": nz.mean(axis=0)[None].repeat(n, 0)}
    out = nc.region_checks(mats, groups)
    assert out["eligible_regions"] == 1 and out["regions_below_floor"] == 1 and out["regions_presence_off"] == 0


def test_grids_are_the_declared_ones():
    assert len(nc.knn_grid()) == 20 and len(nc.ridge_grid()) == 14 and len(nc.mean_grid()) == 2
    assert max(nc.CONFIG["ridge"]["lambda_grid"]) == 1000.0
    assert nr.CONFIG["salt"] != "SYNARCH-recon-2026-11-10"


def test_zero_table():
    Y = np.array([[0.0, 1.0], [0.0, 2.0], [3.0, 0.0]])
    rows = nr.zero_table(Y, np.array(["a", "a", "b"]), ["R1", "R2"])
    assert rows[0]["zero_fraction"] == pytest.approx(2 / 3) and rows[0]["nonzero_samples"] == 1
    assert rows[1]["zero_fraction_animal_mean"] == pytest.approx(0.5)


def test_sources_read_only_through_the_sealed_loader():
    for rel in ("analysis/next_regional_feasibility.py", "analysis/next_regional_core.py"):
        text = (ROOT / rel).read_text()
        for forbidden in ("traces/", "all_neurons_projection_length", "brm_real/truth", "open(world.regional_csv",
                          "axon_engine", "paircache"):
            assert forbidden not in text, (rel, forbidden)
    assert "load_training_table" in (ROOT / "analysis/next_regional_feasibility.py").read_text()
