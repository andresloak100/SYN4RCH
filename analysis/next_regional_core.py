"""Shared, frozen-to-be pieces of the DRAFT protocol reports/next_regional_protocol.md.

One scorer for every number: the training feasibility (analysis/next_regional_feasibility.py) and the run-once
commands (analysis/next_regional_run.py) both use the functions below, so a frozen number and a feasibility
number are computed the same way.

  * the forecast format: per neuron and region a point, K = 199 quantiles at levels (i - 0.5) / 199 and
    p_nonzero (the probability that the published length is > 0);
  * the declared scorer: quantile CRPS = 2 x mean pinball loss over the 199 levels; the 90% central interval
    [Q(0.05), Q(0.95)]; the conditional non-zero interval [Q(F0 + 0.05 (1 - F0)), Q(F0 + 0.95 (1 - F0))] with
    F0 = 1 - p_nonzero; Q = linear interpolation of the 199 quantiles in the level, flat beyond the end levels;
  * the controls MEAN, IDENTITY, KNN, RIDGE (out-of-fold residual distribution) and KNN-CONF (KNN point with
    a Mondrian out-of-fold residual distribution), fitted on training rows only;
  * the decision rule.

Pure numpy; nothing here opens a file of the cohort.
"""
from __future__ import annotations

import hashlib
import math

import numpy as np

K = 199
LEVELS = (np.arange(1, K + 1) - 0.5) / K
LEVEL = 0.90

CONFIG = {
    "quantiles": K,
    "quantile_method": "numpy.quantile method='linear' on the method's sample (ensemble, pool or point + residuals)",
    "level": LEVEL,
    "mean_weightings": ["sample", "neuron"],
    "knn": {"k_grid": [5, 10, 20, 40, 80], "pools": ["all", "same_label"], "scalings": ["raw", "z"]},
    "ridge": {"lambda_grid": [1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0], "degrees": [1, 2]},
    "residual_folds": 3,                 # whole-sample out-of-fold residuals inside the fitting set
    "residual_salt": "SYNARCH-next-regional-residuals-2026-11-11",
    "conformal_bins": [0.0, 1.0 / 3.0, 2.0 / 3.0],   # KNN neighbours' non-zero share per region
    "conformal_min_bin": 20,             # fewer residuals in a bin: the region's pooled residuals
}

RULE = {
    "crps_gain": 0.02,                   # leg 1: CRPS at least 2% below every distribution control
    "mse_noninferiority": 0.01,          # leg 2: log MSE at most 1% above the best point control
    "point_gain_label": 0.02,            # descriptive label "also better points": 2% below every point control
    "coverage_band": [0.85, 0.95],       # leg 3a: animal-mean marginal 90% coverage
    "nonzero_coverage_band": [0.85, 0.95],   # leg 3b: animal-mean conditional non-zero coverage
    "region_nonzero_floor": 0.75,        # leg 3c: pooled conditional non-zero coverage per eligible region
    "region_min_nonzero_truths": 30,
    "region_min_nonzero_samples": 8,
    "region_allowance": 5,
    "presence_tolerance": 0.10,          # leg 3d: |mean p_nonzero - observed non-zero rate| per region
    "presence_allowance": 5,
}

POINT_CONTROLS = ("MEAN", "IDENTITY", "KNN", "RIDGE")
DIST_CONTROLS = ("MEAN", "IDENTITY", "KNN@crps", "RIDGE@crps", "KNN-CONF")


# ----------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def animal_mean(per_neuron, groups):
    """Mean over neurons within each brain sample (NaN ignored), then over samples with a value."""
    per_neuron = np.asarray(per_neuron, dtype=float)
    groups = np.asarray(groups)
    ids = np.array(sorted(set(groups.tolist())))
    values = np.array([np.nanmean(per_neuron[groups == g]) if np.isfinite(per_neuron[groups == g]).any() else np.nan
                       for g in ids])
    return float(np.nanmean(values)), ids, values


def crps_ensemble(ens, y):
    """Exact CRPS of an ensemble's empirical distribution (reference only; the declared scorer is quantile_crps)."""
    ens = np.asarray(ens, dtype=float)
    y = np.asarray(y, dtype=float)
    s = np.sort(ens, axis=-1)
    valid = np.isfinite(s)
    c = valid.sum(axis=-1)
    if (c == 0).any():
        raise ValueError("an ensemble has no finite member")
    sz = np.where(valid, s, 0.0)
    term1 = (np.abs(sz - y[..., None]) * valid).sum(axis=-1) / c
    i = np.arange(1, s.shape[-1] + 1)
    w = (2 * i - c[..., None] - 1) * valid
    spread = 2.0 * (w * sz).sum(axis=-1) / (c.astype(float) ** 2)
    return term1 - 0.5 * spread


def to_quantiles(samples):
    """samples (..., n) with NaN for absent members -> (..., K) quantiles, numpy 'linear' method."""
    samples = np.asarray(samples, dtype=float)
    if not np.isnan(samples).any():
        return np.moveaxis(np.quantile(samples, LEVELS, axis=-1), 0, -1)
    flat = samples.reshape(-1, samples.shape[-1])
    out = np.empty((flat.shape[0], K))
    for i, row in enumerate(flat):
        row = row[np.isfinite(row)]
        if not len(row):
            raise ValueError("a predictive sample has no finite member")
        out[i] = np.quantile(row, LEVELS)
    return out.reshape(samples.shape[:-1] + (K,))


def share_nonzero(samples):
    samples = np.asarray(samples, dtype=float)
    finite = np.isfinite(samples)
    return (np.where(finite, samples > 0, False).sum(axis=-1) / finite.sum(axis=-1)).astype(float)


def quantile_at(q, level):
    """Q(level) by linear interpolation of the K quantiles in the level; flat beyond the end levels.
    q (..., K) non-decreasing; level scalar or broadcastable to q.shape[:-1]."""
    q = np.asarray(q, dtype=float)
    k = q.shape[-1]
    level = np.broadcast_to(np.asarray(level, dtype=float), q.shape[:-1])
    pos = np.clip(level * k - 0.5, 0.0, k - 1.0)
    lo = np.floor(pos).astype(int)
    hi = np.minimum(lo + 1, k - 1)
    frac = pos - lo
    a = np.take_along_axis(q, lo[..., None], axis=-1)[..., 0]
    b = np.take_along_axis(q, hi[..., None], axis=-1)[..., 0]
    return a + frac * (b - a)


def quantile_crps(q, y):
    """Declared CRPS: 2 x the mean pinball loss over the K levels (the quantile approximation of CRPS)."""
    q = np.asarray(q, dtype=float)
    d = np.asarray(y, dtype=float)[..., None] - q
    tau = LEVELS if q.shape[-1] == K else (np.arange(1, q.shape[-1] + 1) - 0.5) / q.shape[-1]
    return 2.0 * np.maximum(tau * d, (tau - 1.0) * d).mean(axis=-1)


def interval_score(lo, hi, y, alpha):
    lo, hi, y = (np.asarray(v, dtype=float) for v in (lo, hi, y))
    return (hi - lo) + (2.0 / alpha) * np.maximum(lo - y, 0.0) + (2.0 / alpha) * np.maximum(y - hi, 0.0)


def intervals(q, p_nonzero, level=LEVEL):
    """Marginal central interval and the conditional non-zero interval (NaN where p_nonzero is 0)."""
    a = (1.0 - level) / 2.0
    lo, hi = quantile_at(q, a), quantile_at(q, 1.0 - a)
    p = np.asarray(p_nonzero, dtype=float)
    f0 = 1.0 - p
    nz_lo = quantile_at(q, f0 + a * p)
    nz_hi = quantile_at(q, f0 + (1.0 - a) * p)
    none = p <= 0
    nz_lo, nz_hi = np.where(none, np.nan, nz_lo), np.where(none, np.nan, nz_hi)
    return lo, hi, nz_lo, nz_hi


def score_forecast(point, q, p_nonzero, Y, level=LEVEL):
    """Per-neuron scores (means over the regions) and the neuron x region matrices the per-region checks pool."""
    point, q, p, Y = (np.asarray(v, dtype=float) for v in (point, q, p_nonzero, Y))
    lo, hi, nz_lo, nz_hi = intervals(q, p, level)
    nz = Y > 0
    inside = (Y >= lo - 1e-9) & (Y <= hi + 1e-9)
    inside_nz = (Y >= nz_lo - 1e-9) & (Y <= nz_hi + 1e-9)          # NaN bounds compare False: not covered
    crps = quantile_crps(q, Y)
    with np.errstate(invalid="ignore", divide="ignore"):
        cov_nz = np.where(nz.any(axis=1), (inside_nz & nz).sum(axis=1) / nz.sum(axis=1), np.nan)
        marg_nz = np.where(nz.any(axis=1), (inside & nz).sum(axis=1) / nz.sum(axis=1), np.nan)
    per_neuron = {
        "log_mse": ((point - Y) ** 2).mean(axis=1), "log_mae": np.abs(point - Y).mean(axis=1),
        "crps": crps.mean(axis=1), "coverage": inside.mean(axis=1),
        "coverage_nonzero_conditional": cov_nz, "coverage_nonzero_marginal": marg_nz,
        "width": (hi - lo).mean(axis=1),
        "interval_score": interval_score(lo, hi, Y, 1 - level).mean(axis=1),
        "presence_brier": ((p - nz) ** 2).mean(axis=1),
    }
    matrices = {"inside": inside, "inside_nz": np.where(nz, inside_nz, np.nan), "p_nonzero": p, "nonzero": nz,
                "crps": crps}
    return per_neuron, matrices


def region_checks(matrices, groups, rule=RULE):
    """Legs 3c and 3d on pooled neurons: per region, conditional non-zero coverage (eligible regions only) and
    presence calibration in the large."""
    groups = np.asarray(groups)
    nz = matrices["nonzero"]
    inside_nz = matrices["inside_nz"]
    rows, low, eligible, presence_bad = [], [], 0, 0
    for r in range(nz.shape[1]):
        n_nz = int(nz[:, r].sum())
        n_samples = len(set(groups[nz[:, r]].tolist()))
        ok = n_nz >= rule["region_min_nonzero_truths"] and n_samples >= rule["region_min_nonzero_samples"]
        cov = float(np.nanmean(inside_nz[:, r])) if n_nz else float("nan")
        gap = float(matrices["p_nonzero"][:, r].mean() - nz[:, r].mean())
        eligible += ok
        if ok and cov < rule["region_nonzero_floor"]:
            low.append(r)
        if abs(gap) > rule["presence_tolerance"]:
            presence_bad += 1
        rows.append({"nonzero_truths": n_nz, "nonzero_samples": n_samples, "eligible": bool(ok),
                     "coverage_nonzero_conditional": cov,
                     "coverage_all": float(np.mean(matrices["inside"][:, r])), "presence_gap": gap})
    return {"per_region": rows, "eligible_regions": int(eligible), "regions_below_floor": len(low),
            "regions_below_floor_index": low, "regions_presence_off": int(presence_bad)}


def decide(animal, region, rule=RULE, candidate="CANDIDATE"):
    """The decision rule (section 8). animal: {method: {metric: animal mean}}; region: region_checks of the
    candidate. Returns every leg and the outcome."""
    c = animal[candidate]
    crps_legs = {k: c["crps"] <= (1 - rule["crps_gain"]) * animal[k]["crps"] for k in DIST_CONTROLS}
    best_point = min(animal[k]["log_mse"] for k in POINT_CONTROLS)
    mse_ok = c["log_mse"] <= (1 + rule["mse_noninferiority"]) * best_point
    point_label = all(c["log_mse"] <= (1 - rule["point_gain_label"]) * animal[k]["log_mse"] for k in POINT_CONTROLS)
    b, nb = rule["coverage_band"], rule["nonzero_coverage_band"]
    calib = {
        "3a_marginal_coverage": b[0] <= c["coverage"] <= b[1],
        "3b_nonzero_conditional_coverage": nb[0] <= c["coverage_nonzero_conditional"] <= nb[1],
        "3c_regions_below_floor": region["regions_below_floor"] <= rule["region_allowance"],
        "3d_presence_calibration": region["regions_presence_off"] <= rule["presence_allowance"],
    }
    leg1, leg3 = all(crps_legs.values()), all(calib.values())
    if leg1 and mse_ok and leg3:
        outcome = "pass"
    elif leg1 and mse_ok:
        outcome = "better_crps_without_valid_calibration"
    elif leg1:
        outcome = "better_crps_but_worse_points"
    else:
        outcome = "no_gain_over_the_controls"
    return {"leg1_crps": crps_legs, "leg1": leg1, "leg2_mse_noninferior": bool(mse_ok),
            "best_point_control_log_mse": best_point, "leg3": calib, "leg3_all": leg3,
            "also_better_points_label": bool(point_label), "outcome": outcome}


def stable_folds(samples, n, salt):
    """Whole samples -> n folds by SHA-256(salt:sample) order, dealt in turn."""
    ordered = sorted(set(samples), key=lambda s: hashlib.sha256(f"{salt}:{s}".encode()).hexdigest())
    return {s: i % n for i, s in enumerate(ordered)}


# ----------------------------------------------------------------------------------------------
# the controls
# ----------------------------------------------------------------------------------------------


def _ridge_solve(F, Y, lam, n_labels):
    centre, spread = F.mean(axis=0), F.std(axis=0)
    spread = np.where(spread > 0, spread, 1.0)
    if n_labels:
        spread[-n_labels:] = 1.0
    G = (F - centre) / spread
    ybar = Y.mean(axis=0)
    beta = np.linalg.solve(G.T @ G + lam * len(Y) * np.eye(G.shape[1]), G.T @ (Y - ybar))
    return centre, spread, ybar, beta


def _ridge_apply(model, F):
    centre, spread, ybar, beta = model
    return np.maximum(ybar + ((F - centre) / spread) @ beta, 0.0)


def _residual_distribution(point, resid_sorted):
    """point (m,) one region; resid_sorted (n,) sorted. -> q (m, K), p (m,)."""
    rq = np.quantile(resid_sorted, LEVELS)
    q = np.maximum(point[:, None] + rq[None, :], 0.0)
    above = len(resid_sorted) - np.searchsorted(resid_sorted, -point, side="right")
    return q, above / len(resid_sorted)


class Controls:
    """MEAN, IDENTITY, KNN, RIDGE and KNN-CONF fitted on one set of TRAINING rows.

    A query enters only through its soma position (AP, DV, ML; mirrored frame) and its source label.
      MEAN      point: training mean ('sample' or 'neuron' weighting); distribution: all training values.
      IDENTITY  the same within the query's source label (unseen label: MEAN).
      KNN       k nearest training neurons by soma ('raw' um or training 'z'), pool 'all' or 'same_label';
                point: their mean; distribution: their values.
      RIDGE     ridge on standardised soma terms (degree 1 or 2) + one-hot label, clipped at 0; distribution:
                point + the region's OUT-OF-FOLD residuals (3 whole-sample folds of the fitting set), clipped at 0.
      KNN-CONF  KNN point (given setting); distribution: point + out-of-fold KNN residuals of the same region and
                the same bin of the neighbours' non-zero share ([0, 1/3), [1/3, 2/3), [2/3, 1]), clipped at 0;
                a bin with fewer than 20 residuals uses the region's pooled residuals.
    """

    def __init__(self, soma, label, group, Y, config=CONFIG):
        self.soma, self.label, self.group, self.Y = (np.asarray(soma, float), np.asarray(label).astype(str),
                                                     np.asarray(group).astype(str), np.asarray(Y, float))
        self.config = config
        self.n = len(self.Y)
        self.mu, sd = self.soma.mean(axis=0), self.soma.std(axis=0)
        self.sd = np.where(sd > 0, sd, 1.0)
        samples = sorted(set(self.group.tolist()))
        self.mean_by = {"sample": np.mean([self.Y[self.group == s].mean(axis=0) for s in samples], axis=0),
                        "neuron": self.Y.mean(axis=0)}
        self.labels = sorted(set(self.label.tolist()))
        self.label_mean_by = {"sample": {}, "neuron": {}}
        for value in self.labels:
            rows = self.label == value
            self.label_mean_by["neuron"][value] = self.Y[rows].mean(axis=0)
            self.label_mean_by["sample"][value] = np.mean(
                [self.Y[rows & (self.group == s)].mean(axis=0) for s in samples if (rows & (self.group == s)).any()],
                axis=0)
        fold_of = stable_folds(samples, config["residual_folds"], config["residual_salt"])
        self.subfold = np.array([fold_of[g] for g in self.group])
        self._ridge, self._ridge_resid, self._conf = {}, {}, {}

    # -- MEAN / IDENTITY ---------------------------------------------------------------------------
    def mean_point(self, m, weighting):
        return np.tile(self.mean_by[weighting], (m, 1))

    def identity_point(self, qlabel, weighting):
        return np.vstack([self.label_mean_by[weighting].get(str(l), self.mean_by[weighting]) for l in qlabel])

    def pool_forecast(self, qlabel):
        """(q, p) of MEAN (qlabel None) or IDENTITY."""
        m = 1 if qlabel is None else len(qlabel)
        R = self.Y.shape[1]
        q, p = np.empty((m, R, K)), np.empty((m, R))
        keys = [None] if qlabel is None else sorted(set(np.asarray(qlabel).astype(str).tolist()))
        for key in keys:
            rows = np.arange(self.n) if key is None else np.flatnonzero(self.label == key)
            if not len(rows):
                rows = np.arange(self.n)
            sel = np.arange(m) if key is None else np.flatnonzero(np.asarray(qlabel).astype(str) == key)
            pool = self.Y[rows]
            q[sel] = np.quantile(pool, LEVELS, axis=0).T[None]
            p[sel] = (pool > 0).mean(axis=0)[None]
        return q, p

    # -- KNN ---------------------------------------------------------------------------------------
    def _scaled(self, X, scaling):
        return X if scaling == "raw" else (X - self.mu) / self.sd

    def neighbours(self, qsoma, qlabel, scaling, pool, k_max, rows=None):
        """Nearest training neurons (restricted to `rows` if given). Returns indices into the full set."""
        from scipy.spatial.distance import cdist
        rows = np.arange(self.n) if rows is None else np.asarray(rows)
        D = cdist(self._scaled(np.asarray(qsoma, float), scaling), self._scaled(self.soma[rows], scaling), "sqeuclidean")
        if pool == "same_label":
            same = np.asarray(qlabel).astype(str)[:, None] == self.label[rows][None, :]
            restricted = np.where(same, D, np.inf)
            empty = ~same.any(axis=1)
            restricted[empty] = D[empty]
            D = restricted
        elif pool != "all":
            raise ValueError(pool)
        k_max = min(int(k_max), len(rows))
        order = np.argsort(D, axis=1, kind="stable")[:, :k_max]
        valid = np.isfinite(np.take_along_axis(D, order, axis=1))
        return rows[order], valid

    def knn_ensemble(self, order, valid, k):
        idx, ok = order[:, :k], valid[:, :k]
        ens = self.Y[idx].transpose(0, 2, 1).copy()        # (m, R, k)
        ens[np.broadcast_to(~ok[:, None, :], ens.shape)] = np.nan
        return ens

    def knn_point(self, qsoma, qlabel, setting, rows=None):
        order, valid = self.neighbours(qsoma, qlabel, setting["scaling"], setting["pool"], setting["k"], rows)
        ens = self.knn_ensemble(order, valid, setting["k"])
        return np.nanmean(ens, axis=-1), share_nonzero(ens)

    # -- RIDGE -------------------------------------------------------------------------------------
    def _design(self, soma, labels, degree):
        Z = (np.asarray(soma, float) - self.mu) / self.sd
        parts = [Z]
        if degree == 2:
            p = Z.shape[1]
            parts.append(Z ** 2)
            parts.extend(Z[:, [i]] * Z[:, [j]] for i in range(p) for j in range(i + 1, p))
        parts.append((np.asarray(labels).astype(str)[:, None] == np.array(self.labels)[None, :]).astype(float))
        return np.column_stack(parts)

    def ridge_point(self, qsoma, qlabel, lam, degree):
        key = (float(lam), int(degree))
        if key not in self._ridge:
            self._ridge[key] = _ridge_solve(self._design(self.soma, self.label, degree), self.Y, lam, len(self.labels))
        return _ridge_apply(self._ridge[key], self._design(qsoma, qlabel, degree))

    def ridge_residuals(self, lam, degree):
        """Out-of-fold residuals (n, R), sorted per region, of RIDGE inside this fitting set."""
        key = (float(lam), int(degree))
        if key not in self._ridge_resid:
            F = self._design(self.soma, self.label, degree)
            resid = np.empty_like(self.Y)
            for h in range(self.config["residual_folds"]):
                held = self.subfold == h
                model = _ridge_solve(F[~held], self.Y[~held], lam, len(self.labels))
                resid[held] = self.Y[held] - _ridge_apply(model, F[held])
            self._ridge_resid[key] = np.sort(resid, axis=0)
        return self._ridge_resid[key]

    def ridge_forecast(self, qsoma, qlabel, lam, degree):
        point = self.ridge_point(qsoma, qlabel, lam, degree)
        resid = self.ridge_residuals(lam, degree)
        q, p = np.empty(point.shape + (K,)), np.empty(point.shape)
        for r in range(point.shape[1]):
            q[:, r], p[:, r] = _residual_distribution(point[:, r], resid[:, r])
        return point, q, p

    # -- KNN-CONF ----------------------------------------------------------------------------------
    def _bin(self, share):
        edges = self.config["conformal_bins"]
        return np.searchsorted(np.asarray(edges[1:]), share, side="right")

    def conformal_residuals(self, setting):
        key = tuple(sorted(setting.items()))
        if key not in self._conf:
            resid, bins = np.empty_like(self.Y), np.empty(self.Y.shape, dtype=int)
            for h in range(self.config["residual_folds"]):
                held = np.flatnonzero(self.subfold == h)
                fit = np.flatnonzero(self.subfold != h)
                point, share = self.knn_point(self.soma[held], self.label[held], setting, rows=fit)
                resid[held], bins[held] = self.Y[held] - point, self._bin(share)
            table = {}
            for r in range(self.Y.shape[1]):
                pooled = np.sort(resid[:, r])
                for b in range(len(self.config["conformal_bins"])):
                    vals = np.sort(resid[bins[:, r] == b, r])
                    table[(r, b)] = vals if len(vals) >= self.config["conformal_min_bin"] else pooled
            self._conf[key] = table
        return self._conf[key]

    def knn_conf_forecast(self, qsoma, qlabel, setting):
        point, share = self.knn_point(qsoma, qlabel, setting)
        bins = self._bin(share)
        table = self.conformal_residuals(setting)
        q, p = np.empty(point.shape + (K,)), np.empty(point.shape)
        for r in range(point.shape[1]):
            for b in np.unique(bins[:, r]):
                sel = bins[:, r] == b
                q[sel, r], p[sel, r] = _residual_distribution(point[sel, r], table[(r, int(b))])
        return point, q, p

    # -- one entry point ----------------------------------------------------------------------------
    def forecast(self, name, setting, qsoma, qlabel):
        """(point, q, p_nonzero) of one control with one setting; name without the '@crps' suffix."""
        m = len(qsoma)
        if name == "MEAN":
            q, p = self.pool_forecast(None)
            return self.mean_point(m, setting["weighting"]), np.repeat(q, m, axis=0), np.repeat(p, m, axis=0)
        if name == "IDENTITY":
            q, p = self.pool_forecast(np.asarray(qlabel))
            return self.identity_point(qlabel, setting["weighting"]), q, p
        if name == "KNN":
            order, valid = self.neighbours(qsoma, qlabel, setting["scaling"], setting["pool"], setting["k"])
            ens = self.knn_ensemble(order, valid, setting["k"])
            return np.nanmean(ens, axis=-1), to_quantiles(ens), share_nonzero(ens)
        if name == "RIDGE":
            return self.ridge_forecast(qsoma, qlabel, setting["lambda"], setting["degree"])
        if name == "KNN-CONF":
            return self.knn_conf_forecast(qsoma, qlabel, setting)
        raise ValueError(name)


def mean_grid():
    return [{"weighting": w} for w in CONFIG["mean_weightings"]]


def knn_grid():
    c = CONFIG["knn"]
    return [{"k": k, "pool": p, "scaling": s} for s in c["scalings"] for p in c["pools"] for k in c["k_grid"]]


def ridge_grid():
    c = CONFIG["ridge"]
    return [{"lambda": l, "degree": d} for d in c["degrees"] for l in c["lambda_grid"]]


GRIDS = {"MEAN": mean_grid, "IDENTITY": mean_grid, "KNN": knn_grid, "RIDGE": ridge_grid}
TUNED = ("MEAN", "IDENTITY", "KNN", "RIDGE", "KNN@crps", "RIDGE@crps")


def grid_criteria(ctrl, qsoma, qlabel, qgroup, Y):
    """Animal-mean log MSE of every setting of every control, and animal-mean CRPS of every KNN and RIDGE
    setting, on one validation block."""
    m = len(qsoma)
    out = {}
    for name in ("MEAN", "IDENTITY"):
        pts = [ctrl.mean_point(m, g["weighting"]) if name == "MEAN" else ctrl.identity_point(qlabel, g["weighting"])
               for g in GRIDS[name]()]
        out[name] = [animal_mean(((p - Y) ** 2).mean(axis=1), qgroup)[0] for p in pts]
    kmax = max(CONFIG["knn"]["k_grid"])
    cache = {}
    mse, crps = [], []
    for g in knn_grid():
        key = (g["scaling"], g["pool"])
        if key not in cache:
            order, valid = ctrl.neighbours(qsoma, qlabel, g["scaling"], g["pool"], kmax)
            cache[key] = ctrl.knn_ensemble(order, valid, kmax)
        ens = cache[key][..., :g["k"]]
        mse.append(animal_mean(((np.nanmean(ens, axis=-1) - Y) ** 2).mean(axis=1), qgroup)[0])
        crps.append(animal_mean(quantile_crps(to_quantiles(ens), Y).mean(axis=1), qgroup)[0])
    out["KNN"], out["KNN@crps"] = mse, crps
    mse, crps = [], []
    for g in ridge_grid():
        point, q, _ = ctrl.ridge_forecast(qsoma, qlabel, g["lambda"], g["degree"])
        mse.append(animal_mean(((point - Y) ** 2).mean(axis=1), qgroup)[0])
        crps.append(animal_mean(quantile_crps(q, Y).mean(axis=1), qgroup)[0])
    out["RIDGE"], out["RIDGE@crps"] = mse, crps
    return out


def choose_settings(soma, label, group, Y, folds):
    """Settings of every tuned control by whole-sample CV over `folds` ({fold: [held samples]}), first minimum
    of the fold-averaged animal-mean criterion (ties to the earlier grid entry). KNN-CONF takes KNN's setting."""
    soma, label, group, Y = np.asarray(soma, float), np.asarray(label).astype(str), np.asarray(group).astype(str), Y
    totals = {}
    for held in folds.values():
        v = np.isin(group, held)
        if set(group[v]) & set(group[~v]):
            raise RuntimeError("a brain sample is on both sides of a fold")
        ctrl = Controls(soma[~v], label[~v], group[~v], Y[~v])
        for key, values in grid_criteria(ctrl, soma[v], label[v], group[v], Y[v]).items():
            totals.setdefault(key, np.zeros(len(values)))
            totals[key] += np.asarray(values) / len(folds)
    chosen = {}
    for key in TUNED:
        base = key.split("@")[0]
        j = int(np.argmin(totals[key]))
        chosen[key] = {"setting": GRIDS[base]()[j], "criterion": "crps" if "@" in key else "log_mse",
                       "value": float(totals[key][j]), "grid_index": j, "grid_size": len(totals[key])}
    chosen["KNN-CONF"] = {"setting": dict(chosen["KNN"]["setting"]), "criterion": "KNN's log-MSE setting"}
    return chosen


def method_forecast(ctrl, key, setting, qsoma, qlabel):
    return ctrl.forecast(key.split("@")[0], setting, qsoma, qlabel)


def minimum_detectable(sd, n, alpha=0.05, power=0.80):
    from scipy import stats
    if n < 2:
        return math.inf
    return float((stats.t.ppf(1 - alpha / 2, n - 1) + stats.t.ppf(power, n - 1)) * sd / math.sqrt(n))


def bootstrap_mean_interval(values, resamples, seed):
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(values), size=(resamples, len(values)))
    means = values[idx].mean(axis=1)
    return [float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))]
