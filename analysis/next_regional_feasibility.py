"""Feasibility numbers for the DRAFT protocol reports/next_regional_protocol.md (TRAINING partition only).

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
        python3 -I analysis/next_regional_feasibility.py run [--out data/interim/next_regional]

What it computes, on the 101 training brain samples (4,448 neurons) of the locked split and nothing else:
  * per-region zero fraction and how many regions have enough non-zero training neurons and samples;
  * whole-brain-sample nested cross-validation (5 outer x 3 inner folds, NEW salt) of the declared controls
    MEAN, IDENTITY, KNN, RIDGE (each tuned by inner log MSE; KNN and RIDGE also by inner CRPS) and KNN-CONF,
    all scored by the declared scorer of analysis/next_regional_core.py (199-quantile CRPS, intervals from
    the quantiles, conditional non-zero coverage, per-region checks of the decision rule's calibration leg);
  * paired per-sample differences, detectable effects at 25 / 40 / 65 samples, the spread of near-KNN point
    differences (KNN-ridge blends of the controls' points, the reviewer's exploratory check reproduced), the
    chance that the CRPS leg passes when the true gain is zero, and the chance that a perfectly calibrated
    forecast fails the per-region calibration checks at test size.

What it never does: it reads the regional table, the soma tables and the manifest ONLY through the general
layer's sealed loader (analysis/synarch_recon.py: rows of validation and test brain samples are dropped before
any value on them is parsed). It opens no trace, runs no engine command, runs no candidate (N1, N2, N3), and
writes only under --out. The general layer's files are loaded read-only (no bytecode is written next to them).

Everything here is EXPLORATORY: the 101 training samples were already used by the general layer's run
(reports/synarch_recon_results.md), by the other agent's regional model and by this protocol's first draft and
its review; the numbers below inform the design of a protocol, they are not evidence for any method.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True        # never leave bytecode next to files this script does not own

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/interim/next_regional"


def _load(name, path):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


nc = _load("next_regional_core", ROOT / "analysis/next_regional_core.py")

CONFIG = {
    "schema": "next_regional_feasibility_v2",
    "salt": "SYNARCH-next-regional-2026-11-11",      # NEW salt: not the general layer's fold assignment
    "outer_folds": 5,
    "inner_folds": 3,
    "scale": "log10(1 + published regional axon length)",
    "core": nc.CONFIG,
    "rule": nc.RULE,
    "bootstrap": {"resamples": 10000, "seed": 20261111},
    "power": {"alpha": 0.05, "power": 0.80, "n_samples": [25, 40, 65]},
    "blend_weights": [0.25, 0.5, 0.75],
    "calibration_null": {"draws": 2000, "samples": 40, "seed": 20261112, "nominal": 0.90},
    "enough_nonzero": [{"neurons": 50, "samples": 10}, {"neurons": 100, "samples": 20}, {"neurons": 200, "samples": 30}],
}
METHODS = ("MEAN", "IDENTITY", "KNN", "RIDGE", "KNN@crps", "RIDGE@crps", "KNN-CONF")

# re-exported for the tests
animal_mean, minimum_detectable, bootstrap_mean_interval = nc.animal_mean, nc.minimum_detectable, nc.bootstrap_mean_interval


def zero_table(Y, groups, regions):
    rows = []
    for r, name in enumerate(regions):
        nz = Y[:, r] > 0
        rows.append({"region": name, "zero_fraction": float(1 - nz.mean()),
                     "zero_fraction_animal_mean": nc.animal_mean((~nz).astype(float), groups)[0],
                     "nonzero_neurons": int(nz.sum()), "nonzero_samples": len({g for g, z in zip(groups, nz) if z})})
    return rows


def load_general_core():
    return _load("next_regional_synarch_recon_ro", ROOT / "analysis/synarch_recon.py")


def ccf_check(regions):
    path = ROOT / "data/raw/atlas/structure_graph_1.json"
    if not path.exists():
        return {"checked": False, "reason": "structure graph not on disk"}
    acronyms = set()

    def walk(node):
        acronyms.add(node["acronym"])
        for child in node.get("children", []):
            walk(child)
    walk(json.loads(path.read_text())["msg"][0])
    return {"checked": True, "structure_graph_sha256": nc.sha_file(path), "missing": [r for r in regions if r not in acronyms]}


def nested_cv(soma, label, group, Y, folds):
    """Outer folds: settings from the inner folds (choose_settings), then every method forecast on the held
    samples and scored by the declared scorer. Returns selections, per-neuron scores, region matrices and
    the outer points of KNN and RIDGE (for the blend check)."""
    n, R = Y.shape
    per_neuron = {m: {} for m in METHODS}
    matrices = {m: {} for m in METHODS}
    points = {m: np.full((n, R), np.nan) for m in ("KNN", "RIDGE")}
    selected = {m: {} for m in METHODS}
    for f in range(CONFIG["outer_folds"]):
        held = np.isin(group, folds["outer"][str(f)])
        tr, te = np.flatnonzero(~held), np.flatnonzero(held)
        if set(group[tr]) & set(group[te]):
            raise RuntimeError("a brain sample is on both sides of an outer fold")
        chosen = nc.choose_settings(soma[tr], label[tr], group[tr], Y[tr], folds["inner"][str(f)])
        ctrl = nc.Controls(soma[tr], label[tr], group[tr], Y[tr])
        for m in METHODS:
            selected[m][str(f)] = chosen[m]
            point, q, p = nc.method_forecast(ctrl, m, chosen[m]["setting"], soma[te], label[te])
            if m in points:
                points[m][te] = point
            scores, mats = nc.score_forecast(point, q, p, Y[te])
            for k, v in scores.items():
                per_neuron[m].setdefault(k, np.full(n, np.nan))[te] = v
            for k, v in mats.items():
                matrices[m].setdefault(k, np.full((n, R), np.nan))[te] = v
            del q
    for m in METHODS:
        matrices[m]["nonzero"] = matrices[m]["nonzero"].astype(bool)
        matrices[m]["inside"] = matrices[m]["inside"].astype(bool)
    return selected, per_neuron, matrices, points


def summarise(per_neuron, group):
    ids = np.array(sorted(set(group.tolist())))
    out, per_sample = {}, {}
    for m, metrics in per_neuron.items():
        out[m], per_sample[m] = {}, {}
        for metric, values in metrics.items():
            overall, _, vals = nc.animal_mean(values, group)
            out[m][metric], per_sample[m][metric] = overall, vals
    return out, per_sample, ids


def calibration_null(Y, group, rule, cfg):
    """Chance that a PERFECTLY calibrated forecast fails legs 3c / 3d at test size: draw 40 training samples,
    take their real per-region non-zero counts, draw each eligible region's conditional coverage as
    Binomial(n_nonzero, 0.90) / n_nonzero and each region's observed non-zero rate as Binomial(n, rate) / n
    around the forecast's rate. Neurons are treated as independent (clustering within a sample would widen
    these; stated)."""
    rng = np.random.default_rng(cfg["seed"])
    samples = np.array(sorted(set(group.tolist())))
    nz = Y > 0
    by_sample = {s: np.flatnonzero(group == s) for s in samples}
    below, presence, eligible = [], [], []
    for _ in range(cfg["draws"]):
        pick = rng.choice(samples, size=cfg["samples"], replace=False)
        rows = np.concatenate([by_sample[s] for s in pick])
        cnt = nz[rows].sum(axis=0)
        nsamp = np.array([len({group[i] for i in rows[nz[rows, r]]}) for r in range(Y.shape[1])])
        ok = (cnt >= rule["region_min_nonzero_truths"]) & (nsamp >= rule["region_min_nonzero_samples"])
        cov = rng.binomial(np.maximum(cnt, 1), cfg["nominal"]) / np.maximum(cnt, 1)
        below.append(int((ok & (cov < rule["region_nonzero_floor"])).sum()))
        rate = nz[rows].mean(axis=0)
        gap = rng.binomial(len(rows), rate) / len(rows) - rate
        presence.append(int((np.abs(gap) > rule["presence_tolerance"]).sum()))
        eligible.append(int(ok.sum()))
    below, presence, eligible = map(np.array, (below, presence, eligible))
    return {"draws": cfg["draws"], "samples_per_draw": cfg["samples"],
            "eligible_regions": {"min": int(eligible.min()), "median": float(np.median(eligible)), "max": int(eligible.max())},
            "regions_below_floor": {"mean": float(below.mean()), "max": int(below.max()),
                                    "p_exceeds_allowance": float((below > rule["region_allowance"]).mean())},
            "regions_presence_off": {"mean": float(presence.mean()), "max": int(presence.max()),
                                     "p_exceeds_allowance": float((presence > rule["presence_allowance"]).mean())},
            "assumption": "neurons independent within a region; clustering by sample would widen both"}


def run(out_dir):
    from scipy import stats
    start = time.time()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    core = load_general_core()
    world = core.real_world(ROOT)
    table = core.load_training_table(world)                       # sealed: training rows only
    group = np.asarray(table.sample).astype(str)
    label = np.asarray(table.label).astype(str)
    soma = np.asarray(table.soma, float)
    Y = table.profile()
    regions = list(table.regions)
    rule = CONFIG["rule"]
    samples = sorted(set(group.tolist()))
    label_of = {s: label[group == s][0] for s in samples}
    folds = core.nested_folds(samples, label_of, outer=CONFIG["outer_folds"], inner=CONFIG["inner_folds"],
                              salt=CONFIG["salt"])

    # 1. zeros and support
    ztab = zero_table(Y, group, regions)
    enough = []
    for t in CONFIG["enough_nonzero"]:
        names = [r["region"] for r in ztab if r["nonzero_neurons"] >= t["neurons"] and r["nonzero_samples"] >= t["samples"]]
        enough.append({**t, "regions": len(names)})
    nonzero_per_neuron = (Y > 0).sum(axis=1)

    # 2. nested CV of the controls, declared scorer
    selected, per_neuron, matrices, points = nested_cv(soma, label, group, Y, folds)
    summary, per_sample, ids = summarise(per_neuron, group)
    regions_check = {m: nc.region_checks(matrices[m], group, rule) for m in METHODS}
    band, nband = rule["coverage_band"], rule["nonzero_coverage_band"]
    calibration_legs = {m: {
        "3a_marginal_coverage": band[0] <= summary[m]["coverage"] <= band[1],
        "3b_nonzero_conditional_coverage": nband[0] <= summary[m]["coverage_nonzero_conditional"] <= nband[1],
        "3c_regions_below_floor": regions_check[m]["regions_below_floor"],
        "3c_passes": regions_check[m]["regions_below_floor"] <= rule["region_allowance"],
        "3d_regions_presence_off": regions_check[m]["regions_presence_off"],
        "3d_passes": regions_check[m]["regions_presence_off"] <= rule["presence_allowance"],
        "eligible_regions_pooled_101_samples": regions_check[m]["eligible_regions"],
        "regions_marginal_coverage_below_0.80": int(sum(r["coverage_all"] < 0.80 for r in regions_check[m]["per_region"])),
    } for m in METHODS}

    # 3. paired differences and detectable effects (sample = unit)
    pc, bs = CONFIG["power"], CONFIG["bootstrap"]
    pairs = {"log_mse": [("KNN", "RIDGE"), ("KNN", "IDENTITY"), ("KNN", "MEAN"), ("RIDGE", "IDENTITY"), ("IDENTITY", "MEAN")],
             "crps": [("KNN@crps", "RIDGE@crps"), ("KNN@crps", "IDENTITY"), ("KNN@crps", "MEAN"), ("RIDGE@crps", "IDENTITY"),
                      ("IDENTITY", "MEAN"), ("KNN-CONF", "KNN@crps"), ("KNN-CONF", "RIDGE@crps")]}
    comparisons = {}
    for metric, plist in pairs.items():
        for a, b in plist:
            d = per_sample[a][metric] - per_sample[b][metric]
            comparisons[f"{metric}:{a}-{b}"] = {
                "mean_difference": float(d.mean()), "relative_to_b": float(d.mean() / per_sample[b][metric].mean()),
                "sd_per_sample": float(np.std(d, ddof=1)), "samples_favouring_a": int((d < 0).sum()), "samples": int(len(d)),
                "bootstrap95": nc.bootstrap_mean_interval(d, bs["resamples"], bs["seed"])}
    # near-KNN point spreads: blends of the controls' own outer points (reproduces the reviewer's check)
    near = {}
    knn_s = np.array([((points["KNN"][group == g] - Y[group == g]) ** 2).mean() for g in ids])
    for w in CONFIG["blend_weights"]:
        blend = w * points["KNN"] + (1 - w) * points["RIDGE"]
        s = np.array([((blend[group == g] - Y[group == g]) ** 2).mean() for g in ids])
        d = s - knn_s
        near[f"blend_w{w}"] = {"log_mse": float(s.mean()), "gain_vs_KNN": float(-d.mean() / knn_s.mean()),
                               "sd_diff_vs_KNN": float(np.std(d, ddof=1))}
    dcrps = per_sample["KNN-CONF"]["crps"] - per_sample["KNN@crps"]["crps"]
    near["KNN-CONF_vs_KNN@crps_crps"] = {"gain": float(-dcrps.mean() / per_sample["KNN@crps"]["crps"].mean()),
                                         "sd_diff": float(np.std(dcrps, ddof=1))}
    sds = {"log_mse_between_controls": sorted(v["sd_per_sample"] for k, v in comparisons.items() if k.startswith("log_mse:")),
           "crps_between_controls": sorted(v["sd_per_sample"] for k, v in comparisons.items() if k.startswith("crps:")),
           "log_mse_near_KNN": sorted(v["sd_diff_vs_KNN"] for k, v in near.items() if k.startswith("blend")),
           "crps_near_KNN": [near["KNN-CONF_vs_KNN@crps_crps"]["sd_diff"]]}
    ref = {"log_mse": summary["KNN"]["log_mse"], "crps": summary["KNN@crps"]["crps"]}
    power = {}
    for key, values in sds.items():
        metric = "log_mse" if key.startswith("log_mse") else "crps"
        power[key] = {str(n): [minimum_detectable(sd, n, pc["alpha"], pc["power"]) / ref[metric] for sd in (min(values), max(values))]
                      for n in pc["n_samples"]}
    two_percent = {m: {str(n): 0.02 * ref[m] / minimum_detectable(1.0, n, pc["alpha"], pc["power"]) for n in pc["n_samples"]}
                   for m in ("log_mse", "crps")}
    # chance the CRPS leg passes against ONE control when the true gain is zero (normal approximation)
    null_pass = {}
    for key in ("crps_near_KNN", "crps_between_controls"):
        for sd in sorted(set([min(sds[key]), max(sds[key])])):
            null_pass[f"{key}:sd={sd:.4f}"] = {str(n): float(1 - stats.norm.cdf(rule["crps_gain"] * ref["crps"] * np.sqrt(n) / sd))
                                              for n in pc["n_samples"]}
    for key in ("log_mse_near_KNN",):
        for sd in sorted(set([min(sds[key]), max(sds[key])])):
            null_pass[f"{key}(old 2% point leg):sd={sd:.4f}"] = {
                str(n): float(1 - stats.norm.cdf(0.02 * ref["log_mse"] * np.sqrt(n) / sd)) for n in pc["n_samples"]}

    cal_null = calibration_null(Y, group, rule, CONFIG["calibration_null"])

    label_counts = {l: {"samples": int(len({s for s in samples if label_of[s] == l})), "neurons": int((label == l).sum())}
                    for l in sorted(set(label.tolist()))}
    result = {
        "schema": CONFIG["schema"], "config": CONFIG, "exploratory": True,
        "statement": ("TRAINING partition only (101 brain samples). Exploratory: these samples were used before by "
                      "the general layer's run, the other agent's regional model and this protocol's first draft and "
                      "review. Nothing here is evidence for any method; it sizes the DRAFT protocol."),
        "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 "
                   "python3 -I analysis/next_regional_feasibility.py run",
        "code_sha256": {"feasibility": nc.sha_file(Path(__file__)), "core": nc.sha_file(ROOT / "analysis/next_regional_core.py")},
        "inputs": {"loader": "analysis/synarch_recon.py load_training_table (sealed)",
                   "loader_sha256": nc.sha_file(ROOT / "analysis/synarch_recon.py"), "audit": table.audit},
        "training": {"neurons": int(len(group)), "samples": len(samples), "labels": label_counts,
                     "training_samples_sha256": table.samples_sha256},
        "regions": regions, "ccf_check": ccf_check(regions),
        "zeros": {"overall_zero_fraction": float((Y == 0).mean()), "all_zero_profiles": int((Y.max(axis=1) == 0).sum()),
                  "nonzero_regions_per_neuron": {"median": float(np.median(nonzero_per_neuron)),
                                                 "p10": float(np.quantile(nonzero_per_neuron, 0.1)),
                                                 "p90": float(np.quantile(nonzero_per_neuron, 0.9))},
                  "enough_nonzero": enough, "per_region": ztab},
        "folds": {"salt": CONFIG["salt"], "outer": folds["outer"]},
        "controls_nested_cv": {"animal_mean": summary, "selected": selected},
        "calibration_legs_of_controls": calibration_legs,
        "per_region_checks": {m: regions_check[m]["per_region"] for m in ("KNN", "KNN@crps", "RIDGE@crps", "KNN-CONF", "MEAN")},
        "paired_differences": comparisons, "near_knn": near, "sd_ranges": sds,
        "minimum_detectable_relative": power, "two_percent_detectable_if_sd_below": two_percent,
        "null_pass_rate_crps_leg_one_control": null_pass, "calibration_null": cal_null,
    }
    text = json.dumps(result, indent=1, sort_keys=True, default=lambda v: v.tolist() if hasattr(v, "tolist") else str(v))
    (out_dir / "feasibility.json").write_text(text + "\n")
    with open(out_dir / "per_region.tsv", "w") as handle:
        cols = ["region", "zero_fraction", "zero_fraction_animal_mean", "nonzero_neurons", "nonzero_samples"]
        extra = [(m, k) for m in ("KNN", "KNN-CONF", "RIDGE@crps")
                 for k in ("coverage_all", "coverage_nonzero_conditional", "presence_gap")]
        handle.write("\t".join(cols + [f"{m}:{k}" for m, k in extra]) + "\n")
        for i, row in enumerate(ztab):
            vals = [row[c] if isinstance(row[c], (str, int)) else f"{row[c]:.4f}" for c in cols]
            vals += [f"{regions_check[m]['per_region'][i][k]:.4f}" for m, k in extra]
            handle.write("\t".join(str(v) for v in vals) + "\n")
    manifest = {name: nc.sha_file(out_dir / name) for name in ("feasibility.json", "per_region.tsv")}
    (out_dir / "SHA256SUMS.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    result["seconds"] = round(time.time() - start, 1)
    return result


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        if os.environ.get(var) != "1":
            print(f"warning: {var} is not 1", file=sys.stderr)
    res = run(args.out)
    print(json.dumps({"animal_mean": res["controls_nested_cv"]["animal_mean"],
                      "calibration_legs": res["calibration_legs_of_controls"], "seconds": res["seconds"]},
                     indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
