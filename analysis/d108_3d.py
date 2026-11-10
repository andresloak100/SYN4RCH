"""D-108 3D: does a constrained choice among SPATIAL feature families, measured in the 3D brain
(analysis/d108_3dbrain.py: real cell clouds, tract voxels, voxel route graph), add to model D?

Leakage: projection data only through d108_prims.world() -> d105_benchmark.read_st2_training(), which
drops Source == "HY" right after reading; HY19 removed from every training source. The Allen tracing
experiments (data/raw/atlas/connectivity) are never opened. Every 3D feature comes from atlas data.

Models (all one weighted ridge, D's head, lambda fixed at D's 0.3; nothing tuned on labels except the
family choice of the search):
  D            gD32 x one-hot (model D)
  D+basic      D + F1 (centroid -log distance, unit location x target centroid)      fixed
  D+3D-all     D + F1..F7                                                             fixed
  D+3D-search  D + families chosen by greedy forward selection inside the outer training divisions
               (inner leave-one-division-out mean r_unit; add a family only if the inner score rises by
               >= ADD_MIN; at most MAX_FAM families; D always kept)
Evaluation: nested leave-one-division-out over the 27 non-HY sources (d105 folds), r_unit per held-out
division (d105 form), mean over 7 divisions.
MODES  python3 -I analysis/d108_3d.py {pilot|declare|planted|null0|null2|real|report}
"""
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import d105_benchmark as b5  # noqa: E402
import d108_prims as P  # noqa: E402
import d108_discovery as DS  # noqa: E402

ROOT = b5.ROOT
OUT = os.path.join(ROOT, "data/interim")
BRAIN = os.path.join(OUT, "d108_3d_brain.npz")
DECL = os.path.join(OUT, "d108_3d_declared.json")
SEED = 20261110
LAM = 0.3
ADD_MIN = 0.01
MAX_FAM = 3
WORKERS = int(os.environ.get("D108_WORKERS", "4"))
FLOOR = 0.05
D_TERM = ("K", "gD32", "onehot")

FAMILIES = {
    "F1_basic": [("P", "distLm"), ("K", "Lm", "tc")],
    "F2_cloud": [("P", "cl_med"), ("P", "cl_p05"), ("P", "cl_ovl"), ("K", "shape6", "tc")],
    "F3_neighbourhood": [("K", "nb6", "tc"), ("P", "nbsim")],
    "F4_compat": [("K", "gD8", "tcomp"), ("P", "gcompat")],
    "F5_corridor": [("P", "corr_shared"), ("K", "tsys6", "tc")],
    "F6_route": [("P", "route"), ("P", "route_gain")],
    "F7_inferred_place": [("P", "distLhat_u"), ("K", "Lhat", "tc")],
}
FAM_TAG = {"F1_basic": "M", "F2_cloud": "M", "F3_neighbourhood": "M", "F4_compat": "M (fitted compatibility)",
           "F5_corridor": "I (tract labels measured, corridor reading assumed)", "F6_route": "I (cost model assumed; sd over cost factors)",
           "F7_inferred_place": "I (M1 out-of-fold, bootstrap sd propagated)"}
NEW_PAIRS = ("cl_med", "cl_p05", "cl_ovl", "nbsim", "gcompat", "corr_shared", "route", "route_gain", "distLhat_u")
NEW_NODES = ("shape6", "nb6", "tsys6")
MODELS = ("D", "basic", "all", "search")
PLANTED3 = {"route": ("P", "route"), "nbhd": ("K", "nb6", "tc")}
PLANTED_FAMILY = {"route": "F6_route", "nbhd": "F3_neighbourhood"}
PLANTED_R = 0.40

_G = {}


def say(*a):
    print(*a, flush=True)


class Engine3D(P.Engine):
    PAIRS = P.Engine.PAIRS + NEW_PAIRS

    def __init__(self, S):
        z = np.load(BRAIN, allow_pickle=True)
        self.b_targets = list(z["targets"])
        keys = [tuple(k.split("|", 1)) for k in z["unit_keys"]]
        self.UI = {k: a for a, k in enumerate(keys)}
        self.U = {n: z["U_" + n] for n in NEW_PAIRS + NEW_NODES}
        # unit nodes: centred, one global scale over the 27 sources' units (label-free, as Phase 1's unit_scale)
        rows = [self.UI[(k, u)] for k, s in S.items() for u in s["units"]]
        self.Un = {}
        for n in NEW_NODES:
            M = self.U[n]
            mu = M[rows].mean(0)
            sd = M[rows].std(0)
            sd[sd == 0] = 1
            Mz = (M - mu) / sd
            c = np.sqrt(((Mz[rows] - Mz[rows].mean(0)) ** 2).sum(1).mean())
            self.Un[n] = Mz / c
        super().__init__(tab=None, S=S)
        assert self.b_targets == self.targets

    def _raw(self, node, src, u, ctx=None):
        if node in NEW_NODES:
            return self.Un[node][self.UI[(src, u)]]
        return super()._raw(node, src, u, ctx)

    def _pair_raw(self, kind, src, units, targets):
        if kind not in NEW_PAIRS:
            return super()._pair_raw(kind, src, units, targets)
        tix = [self.TI[t] for t in targets]
        M = np.array([self.U[kind][self.UI[(src, u)]][tix] for u in units], float)
        if kind in ("cl_med", "cl_p05"):
            M = -np.log(np.maximum(M, FLOOR))
        return b5.dc(M)


def fam_terms(fams):
    return [D_TERM] + [t for f in fams for t in FAMILIES[f]]


def ctx():
    if not _G:
        tab, S, targets, names = P.world()
        assert "HY" not in S
        _G["S"], _G["targets"] = S, targets
        _G["eng"] = Engine3D(S)
        _G["data"] = {}
        _G["pilot1"] = json.load(open(os.path.join(OUT, "d108_pilot.json")))
    return _G


# ------------------------------------------------------------------------------ datasets
def d_truth(S):
    eng = ctx()["eng"]
    f = P.Fold(eng, S, S, use_E=False)
    pr = f.predict([D_TERM], f.solve([D_TERM], [0.3])[0], side="tr")
    pr = {k: b5.dc(v) for k, v in pr.items()}
    rms = np.sqrt(np.mean(np.concatenate([v.ravel() for v in pr.values()]) ** 2))
    return {k: v / rms for k, v in pr.items()}


def planted_signal(S, kind):
    eng = ctx()["eng"]
    term = PLANTED3[kind]
    sig = {}
    if term[0] == "K":
        rng = np.random.default_rng(SEED + 31)
        k0 = next(iter(S))
        Bm = rng.normal(size=(eng.unit_X(term[1], k0, S[k0]["units"]).shape[1], eng.enc[term[2]].shape[1]))
        for k, s in S.items():
            sig[k] = b5.dc(eng.unit_X(term[1], k, s["units"]) @ Bm @ eng.enc_Q(term[2], s["targets"]).T)
    else:
        for k, s in S.items():
            sig[k] = eng.pair_D(term[1], k, s["units"], s["targets"])
    rms = np.sqrt(np.mean(np.concatenate([v.ravel() for v in sig.values()]) ** 2))
    return {k: v / rms for k, v in sig.items()}


def calibrate(S, kind, target=PLANTED_R):
    sig = planted_signal(S, kind)
    lo, hi = 0.0, 30.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if DS.true_r(S, sig, mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def dataset(spec):
    G = ctx()
    key = json.dumps(spec, sort_keys=True)
    if key in G["data"]:
        return G["data"][key]
    S = G["S"]
    kind = spec["kind"]
    if kind == "real":
        D = S
    elif kind == "N0":                                  # same seeds as Phase 1 -> the same shuffled datasets
        rng = np.random.default_rng(SEED + 1000 + spec["rep"])
        D = DS.with_counts(S, {k: s["counts"][rng.permutation(len(s["units"]))] for k, s in S.items()})
    elif kind == "N2":                                  # same seeds and beta as Phase 1 -> the same D-as-truth datasets
        rng = np.random.default_rng(SEED + 2000 + spec["rep"])
        D = DS.with_counts(S, DS.sample(S, DS.d2_logits(S, d_truth(S), spec["beta"], rng), rng))
    elif kind == "planted":
        rng = np.random.default_rng(SEED + 5000 + 100 * list(PLANTED3).index(spec["truth"]) + spec["rep"])
        D = DS.with_counts(S, DS.sample(S, DS.planted_logits(S, planted_signal(S, spec["truth"]), spec["beta"]), rng))
    else:
        raise KeyError(kind)
    G["data"] = {key: D}
    return D


# --------------------------------------------------------------------------------- models
def inner_score(folds, terms):
    return float(DS.inner_scores(folds, terms, (LAM,))[0][0])


def search(folds):
    """greedy forward selection of families on top of D (D always kept), fixed lambda."""
    sel = []
    best = inner_score(folds, fam_terms([]))
    base = best
    step1, path = {}, []
    for step in range(MAX_FAM):
        rows = []
        for f in FAMILIES:
            if f in sel:
                continue
            sc = inner_score(folds, fam_terms(sel + [f]))
            rows.append((sc, f))
            if step == 0:
                step1[f] = sc
        rows.sort(key=lambda x: -x[0])
        sc, f = rows[0]
        if sc < best + ADD_MIN:
            break
        sel.append(f)
        best = sc
        path.append((list(sel), sc))
    return sel, best, base, step1, path


def outer_task(args):
    spec, d, full = args
    t0 = time.time()
    G = ctx()
    eng = G["eng"]
    S = dataset(spec)
    _, tr, te = [x for x in P.folds(S) if x[0] == d][0]
    inner = [P.Fold(eng, itr, ite, use_E=False) for _, itr, ite in P.folds(tr)]
    sel, inner_best, inner_D, step1, path = search(inner)
    fo = P.Fold(eng, tr, te, use_E=False)
    terms = {"D": fam_terms([]), "basic": fam_terms(["F1_basic"]), "all": fam_terms(list(FAMILIES)), "search": fam_terms(sel)}
    pred = {m: fo.predict(t, fo.solve(t, [LAM])[0]) for m, t in terms.items()}
    res = {"spec": spec, "division": d, "families": sel, "inner": inner_best, "inner_D": inner_D, "step1": step1,
           "path": path, "inner_all": inner_score(inner, terms["all"]), "inner_basic": inner_score(inner, terms["basic"])}
    for m in MODELS:
        res["r_" + m] = P.r_division(te, pred[m])
    if spec["kind"] == "planted":
        ot = [PLANTED3[spec["truth"]]]
        res["r_oracle"] = P.r_division(te, fo.predict(ot, fo.solve(ot, [LAM])[0]))
    if full:
        res["per_source"] = {k: {m: P.r_division({k: te[k]}, pred[m]) for m in MODELS} for k in te}
        uc = {m: P.unit_cos(te, pred[m]) for m in MODELS}
        res["unit_cos"] = {"%s|%s" % k: {m: uc[m][k] for m in MODELS} for k in uc["D"]}
        # descriptive (not used for any decision): D + each single family on the held-out division
        res["single_family_outer"] = {f: P.r_division(te, fo.predict(fam_terms([f]), fo.solve(fam_terms([f]), [LAM])[0]))
                                      for f in FAMILIES}
        res["ablation"] = {}
        for f in sel:
            rest = fam_terms([x for x in sel if x != f])
            res["ablation"][f] = res["r_search"] - P.r_division(te, fo.predict(rest, fo.solve(rest, [LAM])[0]))
        res["pred_search"] = {k: v.tolist() for k, v in pred["search"].items()}
        res["pred_all"] = {k: v.tolist() for k, v in pred["all"].items()}
    res["sec"] = time.time() - t0
    return res


def run_specs(specs, full=False, label=""):
    from multiprocessing import get_context
    G = ctx()
    divs = sorted({s["division"] for s in G["S"].values()})
    jobs = [(sp, d, full) for sp in specs for d in divs]
    out = {}
    t0 = time.time()
    with get_context("spawn").Pool(WORKERS) as pool:
        for i, r in enumerate(pool.imap_unordered(outer_task, jobs)):
            out.setdefault(json.dumps(r["spec"], sort_keys=True), []).append(r)
            if (i + 1) % 28 == 0 or i + 1 == len(jobs):
                say("%s: %d/%d outer folds done  [%.0fs]" % (label, i + 1, len(jobs), time.time() - t0))
    for k in out:
        out[k].sort(key=lambda r: r["division"])
    return out


def save(name, obj):
    p = os.path.join(OUT, "d108_3d_" + name + ".json")
    json.dump(obj, open(p, "w"))
    say("wrote", os.path.relpath(p, ROOT))


def load(name):
    return json.load(open(os.path.join(OUT, "d108_3d_" + name + ".json")))


def mean7(folds, m):
    return float(np.mean([f["r_" + m] for f in folds]))


# ---------------------------------------------------------------------------------- modes
def pilot():
    """timing on SHUFFLED labels only (N0 replicate 999) + reproduction checks + planted calibration."""
    G = ctx()
    S = G["S"]
    t0 = time.time()
    res = run_specs([{"kind": "N0", "rep": 999}], label="pilot")
    wall = time.time() - t0
    folds = list(res.values())[0]
    out = {"wall_s": wall, "cpu_s": sum(f["sec"] for f in folds), "workers": WORKERS,
           "mean7": {m: mean7(folds, m) for m in MODELS}, "families": [f["families"] for f in folds]}
    say("pilot N0#999: %s; wall %.0fs, fold CPU %.0fs" % ({m: round(v, 3) for m, v in out["mean7"].items()}, wall, out["cpu_s"]))
    # reproduction: D on Phase 1's shuffled replicate 0 must equal Phase 1's r_D fold by fold
    n0 = json.load(open(os.path.join(OUT, "d108_null0.json")))[json.dumps({"kind": "N0", "rep": 0}, sort_keys=True)]
    r0 = list(run_specs([{"kind": "N0", "rep": 0}], label="repro").values())[0]
    diff = max(abs(a["r_D"] - b["r_D"]) for a, b in zip(r0, n0))
    out["repro_N0rep0_maxdiff_rD"] = diff
    say("reproduction of Phase 1 D on N0 replicate 0: max |diff| %.2e" % diff)
    assert diff < 1e-9
    out["planted_beta"] = {}
    for k in PLANTED3:
        bt = calibrate(S, k)
        out["planted_beta"][k] = bt
        say("planted %s: beta %.3f -> true r %.3f" % (k, bt, DS.true_r(S, planted_signal(S, k), bt)))
    out["N2_beta"] = G["pilot1"]["N2_beta"]
    save("pilot", out)


def declare():
    if os.path.exists(DECL):
        sys.exit("already declared (%s)" % DECL)
    pl = load("pilot")
    dec = {"declared_at": time.strftime("%Y-%m-%d %H:%M:%S"), "workers": WORKERS, "lambda": LAM, "add_min": ADD_MIN,
           "max_families": MAX_FAM, "families": {f: [list(t) for t in v] for f, v in FAMILIES.items()}, "family_tags": FAM_TAG,
           "n_N0": 60, "n_N2": 60, "N2_beta": pl["N2_beta"], "n_planted_per_truth": 8,
           "planted": {k: {"term": list(PLANTED3[k]), "family": PLANTED_FAMILY[k], "beta": pl["planted_beta"][k], "true_r": PLANTED_R}
                       for k in PLANTED3},
           "pilot": {"wall_s": pl["wall_s"], "cpu_s": pl["cpu_s"]}}
    json.dump(dec, open(DECL, "w"), indent=1)
    say(json.dumps(dec, indent=1))


def declared():
    return json.load(open(DECL))


def run_planted():
    dec = declared()
    specs = [{"kind": "planted", "truth": k, "beta": v["beta"], "rep": r} for k, v in dec["planted"].items()
             for r in range(dec["n_planted_per_truth"])]
    save("planted", run_specs(specs, label="planted"))


def run_null0():
    save("null0", run_specs([{"kind": "N0", "rep": r} for r in range(declared()["n_N0"])], label="N0"))


def run_null2():
    dec = declared()
    save("null2", run_specs([{"kind": "N2", "beta": dec["N2_beta"], "rep": r} for r in range(dec["n_N2"])], label="N2"))


def run_real():
    declared()
    marker = os.path.join(OUT, "d108_3d_real_opened.json")
    if os.path.exists(marker):
        sys.exit("real-label run already done (%s)" % marker)
    for n in ("planted", "null0", "null2"):
        assert os.path.exists(os.path.join(OUT, "d108_3d_%s.json" % n)), n
    json.dump({"time": time.strftime("%Y-%m-%d %H:%M:%S")}, open(marker, "w"))
    save("real", run_specs([{"kind": "real"}], full=True, label="real"))


# --------------------------------------------------------------------------------- report
def verdict_vs_D(dl, n_better, p_div, p_n2, drop_ok):
    if dl <= 0:
        return "NO GAIN"
    if p_div < 0.05 and p_n2 < 0.05 and n_better >= 5 and drop_ok:
        return "RELIABLE GAIN"
    if n_better >= 5:
        return "ENCOURAGING, NOT RELIABLE"
    return "NOT SHOWN"


def pair_stats(real, a, b, seed=0):
    da = np.array([f["r_" + a] - f["r_" + b] for f in real])
    src = [v[a] - v[b] for f in real for v in f["per_source"].values()]
    un = [v[a] - v[b] for f in real for v in f["unit_cos"].values()]
    drop = [float(np.delete(da, i).mean()) for i in range(len(da))]
    return {"per_division": da.tolist(), "mean_delta": float(da.mean()), "n_better": int((da > 0).sum()),
            "p_division": DS.signflip(da), "drop_one_min": float(min(drop)), "drop_one_all_positive": bool(min(drop) > 0),
            "n_sources": len(src), "mean_delta_source": float(np.mean(src)), "p_source": DS.signflip(src, seed=seed + 1),
            "n_units": len(un), "mean_delta_unit": float(np.mean(un)), "p_unit": DS.signflip(un, seed=seed + 2)}


def report():
    dec = declared()
    real = load("real")[json.dumps({"kind": "real"}, sort_keys=True)]
    n0 = load("null0")
    n2 = load("null2")
    pl = load("planted")
    out = {"declared": dec}
    out["per_division"] = {f["division"]: {m: f["r_" + m] for m in MODELS} for f in real}
    out["mean7"] = {m: mean7(real, m) for m in MODELS}
    out["families_per_fold"] = [{"division": f["division"], "families": f["families"], "inner": f["inner"], "inner_D": f["inner_D"],
                                 "step1": f["step1"], "ablation": f["ablation"], "single_family_outer": f["single_family_outer"]}
                                for f in real]
    cnt = {}
    for f in real:
        for x in f["families"]:
            cnt[x] = cnt.get(x, 0) + 1
    out["family_counts"] = cnt
    # nulls
    def dist(sets, m, ref="D"):
        return [mean7(v, m) - (mean7(v, ref) if ref else 0) for v in sets.values()]
    out["N0"] = {m: {"mean7_mean": float(np.mean(dist(n0, m, None))), "mean7_95": float(np.percentile(dist(n0, m, None), 95)),
                     "p_N0": float((1 + sum(x >= out["mean7"][m] for x in dist(n0, m, None))) / (1 + len(n0)))} for m in MODELS}
    mc = []
    for v in n0.values():
        c = {}
        for f in v:
            for x in f["families"]:
                c[x] = c.get(x, 0) + 1
        mc.append(max(c.values()) if c else 0)
    out["N0"]["maxcount_95"] = float(np.percentile(mc, 95))
    out["N0"]["maxcount_dist"] = np.bincount(mc, minlength=8).tolist()
    out["N0"]["frac_folds_any_family"] = float(np.mean([len(f["families"]) > 0 for v in n0.values() for f in v]))
    out["N0"]["sec_total"] = float(sum(f["sec"] for v in n0.values() for f in v))
    out["N2"] = {}
    for m in ("basic", "all", "search"):
        d2 = dist(n2, m)
        obs = out["mean7"][m] - out["mean7"]["D"]
        out["N2"][m] = {"delta_mean": float(np.mean(d2)), "delta_sd": float(np.std(d2, ddof=1)), "delta_95": float(np.percentile(d2, 95)),
                        "frac_gt0": float(np.mean(np.array(d2) > 0)), "p_N2": float((1 + sum(x >= obs for x in d2)) / (1 + len(d2)))}
    out["N2"]["D_mean7_mean"] = float(np.mean([mean7(v, "D") for v in n2.values()]))
    out["N2"]["frac_folds_any_family"] = float(np.mean([len(f["families"]) > 0 for v in n2.values() for f in v]))
    out["N2"]["sec_total"] = float(sum(f["sec"] for v in n2.values() for f in v))
    # paired tests
    out["tests"] = {}
    for a, b in (("basic", "D"), ("all", "D"), ("search", "D"), ("all", "basic"), ("search", "basic"), ("search", "all")):
        out["tests"]["%s_vs_%s" % (a, b)] = pair_stats(real, a, b)
    # verdicts
    out["verdict"] = {}
    for m in ("all", "search", "basic"):
        t = out["tests"]["%s_vs_D" % m]
        out["verdict"]["%s_vs_D" % m] = verdict_vs_D(t["mean_delta"], t["n_better"], t["p_division"], out["N2"][m]["p_N2"], t["drop_one_all_positive"])
    for m in ("all", "search"):
        t = out["tests"]["%s_vs_basic" % m]
        out["verdict"]["%s_vs_basic" % m] = ("NO GAIN" if t["mean_delta"] <= 0 else
                                              "RELIABLE GAIN" if (t["p_division"] < 0.05 and t["n_better"] >= 5 and t["drop_one_all_positive"]) else
                                              "ENCOURAGING, NOT RELIABLE" if t["n_better"] >= 5 else "NOT SHOWN")
    # the user's prediction
    pred = {}
    mD = out["mean7"]["D"]
    for m in ("all", "search"):
        dl = out["mean7"][m] - mD
        t = out["tests"]["%s_vs_D" % m]
        pred[m] = {"delta": dl, "relative_gain": dl / mD,
                   "P1_similar_or_slightly_better": bool(-0.02 <= dl <= 0.05),
                   "P2_reliable_gain_10_20pct": bool(out["verdict"]["%s_vs_D" % m] == "RELIABLE GAIN" and 0.10 <= dl / mD <= 0.20),
                   "P3_improves_in_5_6_of_7": bool(t["n_better"] in (5, 6)), "n_better": t["n_better"]}
    out["user_prediction"] = pred
    # planted
    prow = {}
    for k, v in pl.items():
        sp = json.loads(k)
        fam = PLANTED_FAMILY[sp["truth"]]
        c = sum(fam in f["families"] for f in v)
        prow.setdefault(sp["truth"], []).append({"rep": sp["rep"], "count": c, "recovered": c >= 4,
                                                 **{m: mean7(v, m) for m in MODELS}, "oracle": mean7(v, "oracle"),
                                                 "other": sorted({x for f in v for x in f["families"] if x != fam})})
    out["planted"] = {t: {"n": len(r), "counts": [x["count"] for x in sorted(r, key=lambda y: y["rep"])],
                          "frac_recovered": float(np.mean([x["recovered"] for x in r])),
                          **{m + "_mean": float(np.mean([x[m] for x in r])) for m in MODELS + ("oracle",)},
                          "p_search_gt_D_reps": DS.signflip([x["search"] - x["D"] for x in r]),
                          "p_all_gt_D_reps": DS.signflip([x["all"] - x["D"] for x in r])} for t, r in prow.items()}
    out["planted_sec"] = float(sum(f["sec"] for v in pl.values() for f in v))
    out["real_sec"] = float(sum(f["sec"] for f in real))
    save("summary", out)
    say(json.dumps({k: v for k, v in out.items() if k not in ("declared", "families_per_fold")}, indent=1))
    for f in out["families_per_fold"]:
        say(f["division"], f["families"], "inner %.3f (D %.3f)" % (f["inner"], f["inner_D"]),
            {k: round(v, 3) for k, v in f["single_family_outer"].items()}, {k: round(v, 3) for k, v in f["ablation"].items()},
            {k: round(v, 3) for k, v in sorted(f["step1"].items(), key=lambda x: -x[1])})


def main():
    modes = {"pilot": pilot, "declare": declare, "planted": run_planted, "null0": run_null0, "null2": run_null2,
             "real": run_real, "report": report}
    if len(sys.argv) != 2 or sys.argv[1] not in modes:
        sys.exit("usage: python3 -I analysis/d108_3d.py {%s}" % "|".join(modes))
    modes[sys.argv[1]]()


if __name__ == "__main__":
    main()
