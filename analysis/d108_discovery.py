"""D-108 Phase 1: can an automatic search discover feature transformations that predict where a
cell type's neurons project better than model D?  (Epi-retro-seq, 27 non-hypothalamic sources.)

No hypothalamic-source row is read anywhere (d108_prims.world -> d105_benchmark.read_st2_training,
which drops Source == "HY" right after reading; HY19 removed from every training source).

The search (per outer fold, using only the outer fold's training divisions):
  typed primitives  unit-side nodes (gene bases: D's whole-brain PCA 32/8, top-2000-gene PCs, B's 10 PCs,
                    NMF programmes, autoencoder latent; nonlinear expansions: random Fourier features,
                    degree-2 polynomials; measured location; inferred location (module M1: genes ->
                    CCF), inferred division composition (M2); inferred methylation L1 composition
                    (G -> E-hat, cross-fitted in the fold, OPTIONAL like every other node); gene x
                    location products) x target-side encodings (one-hot = D's head, target
                    centroid, target gene profile, target cell-type composition), plus pairwise
                    unit x target features (-log distance from measured or inferred location).
                    This is a FIXED typed library of 66 pre-built terms: chains (genes -> inferred
                    location -> distance to target or -> RFF) and cross-modality products are
                    hand-enumerated nodes; the search selects and combines terms additively, it does
                    not compose new chains or products (no gene x methylation or location x methylation
                    term, no product or multiple-kernel blocks).
  search            greedy forward selection of up to MAX_TERMS terms; each candidate pipeline is
                    fitted jointly (one ridge) and scored by mean r_unit over leave-one-division-out
                    folds INSIDE the outer training divisions; lambda chosen from LAMS by the same
                    inner score; a term is added only if the inner score rises by >= ADD_MIN.
  outer score       r_unit (d105 devkit form) of the selected pipeline on the held-out division;
                    mean over the 7 divisions.
MODES  python3 -I analysis/d108_discovery.py MODE
  pilot | planted | null0 | null2 | real | final | report   (see main())
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

ROOT = b5.ROOT
OUT = os.path.join(ROOT, "data/interim")
REPORT = os.path.join(ROOT, "reports/d108_phase1.md")
DECLARED = os.path.join(OUT, "d108_declared.json")
SEED = 20261110
MAX_TERMS = 3
ADD_MIN = 0.01
WORKERS = int(os.environ.get("D108_WORKERS", "5"))
TAU = 0.35                      # experiment-level noise of the D-as-truth simulation (d105 POWER_TAU)
PLANTED_R = 0.40

_G = {}


def say(*a):
    print(*a, flush=True)


def ctx():
    """per-process world and engine (built once)."""
    if not _G:
        tab, S, targets, names = P.world()
        _G["tab"], _G["S"], _G["targets"] = tab, S, targets
        eng = P.Engine(tab, S)
        z = np.load(os.path.join(OUT, "d105_models.npz"), allow_pickle=True)
        eng.F["Kraw"] = z["F_K"]
        eng.F["Braw"] = z["F_B"]
        _G["eng"] = eng
        _G["terms"] = [("K", n, e) for n in eng.UNIT_NODES for e in P.ENCODINGS] + [("P", k) for k in eng.PAIRS]
        _G["data"] = {}
    return _G


D_TERM = ("K", "gD32", "onehot")


# ------------------------------------------------------------------------------ datasets
def with_counts(S, newc):
    out = {}
    for k, s in S.items():
        c = newc[k]
        out[k] = dict(s, counts=c, n=c.sum(axis=1))
    return out


def shares(c):
    return (c + b5.PSEUDO) / (c.sum(axis=0) + b5.PSEUDO * c.shape[0])


def softmax_cols(E):
    E = E - E.max(axis=0)
    q = np.exp(E)
    return q / q.sum(axis=0)


def sample(S, logits, rng):
    """counts ~ Multinomial(real experiment total, softmax over units of the logits)."""
    newc = {}
    for k, s in S.items():
        q = softmax_cols(logits[k])
        N = s["counts"].sum(axis=0)
        c = np.zeros_like(s["counts"])
        for e in range(len(N)):
            c[:, e] = rng.multinomial(int(N[e]), q[:, e])
        newc[k] = c
    return newc


def expected(S, logits):
    return {k: softmax_cols(logits[k]) * s["counts"].sum(axis=0) for k, s in S.items()}


def d_truth(S):
    """D refitted on all 27 sources, its in-sample profile per source (double-centred), RMS 1."""
    eng = ctx()["eng"]
    f = P.Fold(eng, S, S, use_E=False)
    th = f.solve([D_TERM], [0.3])[0]
    pr = f.predict([D_TERM], th, side="tr")
    pr = {k: b5.dc(v) for k, v in pr.items()}
    rms = np.sqrt(np.mean(np.concatenate([v.ravel() for v in pr.values()]) ** 2))
    return {k: v / rms for k, v in pr.items()}


def planted_signal(S, kind):
    """the planted truths (fixed coefficients, seed per kind), double-centred, RMS 1 over all observations."""
    eng = ctx()["eng"]
    rng = np.random.default_rng(SEED + {"GxL_tc": 11, "chain_distLhat": 12, "nmf_tg": 13, "Ehat_tc": 14}[kind])
    term = PLANTED[kind]
    sig = {}
    if term[0] == "K":
        # E-hat (supplementary truth): the conversion fitted on all 27 sources, cross-fitted over subclasses
        cx = P.Fold(eng, S, {}).ctx if term[1] == "Ehat" else None
        X0 = eng.unit_X(term[1], next(iter(S)), S[next(iter(S))]["units"], cx)
        Bm = rng.normal(size=(X0.shape[1], eng.enc[term[2]].shape[1]))
        for k, s in S.items():
            X = eng.unit_X(term[1], k, s["units"], cx)
            Q = eng.enc_Q(term[2], s["targets"])
            sig[k] = b5.dc(X @ Bm @ Q.T)
    else:
        for k, s in S.items():
            sig[k] = eng.pair_D(term[1], k, s["units"], s["targets"])
    rms = np.sqrt(np.mean(np.concatenate([v.ravel() for v in sig.values()]) ** 2))
    return {k: v / rms for k, v in sig.items()}


PLANTED = {"GxL_tc": ("K", "gD4xLm", "tc"), "chain_distLhat": ("P", "distLhat"), "nmf_tg": ("K", "nmf16", "tg"),
           "Ehat_tc": ("K", "Ehat", "tc")}
PLANTED_DECLARED = ("GxL_tc", "chain_distLhat", "nmf_tg")
# supplementary (post-audit) planted truths: an E-hat (methylation) route alone at true r 0.40, and
# "D + route" truths: D's truth as in N2 plus one extra route of modest size (route true r alone)
SUPP_ALONE = (("Ehat_tc", 0.40),)
SUPP_DPLUS = (("chain_distLhat", 0.15), ("chain_distLhat", 0.25), ("Ehat_tc", 0.25))
ORACLE_MARGIN = 0.10


def lodo_mean(S, preds):
    divs = sorted({s["division"] for s in S.values()})
    return float(np.mean([P.r_division({k: S[k] for k in S if S[k]["division"] == d}, preds) for d in divs]))


def planted_logits(S, sig, beta):
    out = {}
    for k, s in S.items():
        base = np.log(shares(s["counts"]))
        tix = [s["targets"].index(t) for t in s["exp_target"]]
        out[k] = base + beta * sig[k][:, tix]
    return out


def d2_logits(S, truth, beta, rng):
    out = {}
    for k, s in S.items():
        a = np.log(shares(s["counts"]).mean(axis=1))
        tix = [s["targets"].index(t) for t in s["exp_target"]]
        out[k] = a[:, None] + beta * truth[k][:, tix] + rng.normal(0, TAU, size=s["counts"].shape)
    return out


def true_r(S, sig, beta):
    ex = with_counts(S, expected(S, planted_logits(S, sig, beta)))
    return lodo_mean(ex, sig)


def calibrate_planted(S, kind, target=PLANTED_R):
    sig = planted_signal(S, kind)
    lo, hi = 0.0, 30.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if true_r(S, sig, mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def d_lodo(S):
    eng = ctx()["eng"]
    preds = {}
    for d, tr, te in P.folds(S):
        f = P.Fold(eng, tr, te, use_E=False)
        preds.update(f.predict([D_TERM], f.solve([D_TERM], [0.3])[0]))
    return lodo_mean(S, preds)


def dataset(spec):
    G = ctx()
    key = json.dumps(spec, sort_keys=True)
    if key in G["data"]:
        return G["data"][key]
    S = G["S"]
    kind = spec["kind"]
    if kind == "real":
        D = S
    elif kind == "N0":
        rng = np.random.default_rng(SEED + 1000 + spec["rep"])
        D = with_counts(S, {k: s["counts"][rng.permutation(len(s["units"]))] for k, s in S.items()})
    elif kind == "N2":
        rng = np.random.default_rng(SEED + 2000 + spec["rep"])
        truth = d_truth(S)
        D = with_counts(S, sample(S, d2_logits(S, truth, spec["beta"], rng), rng))
    elif kind == "planted":
        rng = np.random.default_rng(SEED + 3000 + 100 * list(PLANTED).index(spec["truth"]) + spec["rep"])
        sig = planted_signal(S, spec["truth"])
        D = with_counts(S, sample(S, planted_logits(S, sig, spec["beta"]), rng))
    elif kind == "plantedD":
        # D's truth exactly as in N2 (beta_D, experiment noise TAU) plus one planted route
        rng = np.random.default_rng(SEED + 4000 + 100 * [list(x) for x in SUPP_DPLUS].index([spec["truth"], spec["r"]]) + spec["rep"])
        truth = d_truth(S)
        sig = planted_signal(S, spec["truth"])
        lg = d2_logits(S, truth, spec["betaD"], rng)
        for k, s in S.items():
            tix = [s["targets"].index(t) for t in s["exp_target"]]
            lg[k] = lg[k] + spec["beta"] * sig[k][:, tix]
        D = with_counts(S, sample(S, lg, rng))
    else:
        raise KeyError(kind)
    G["data"] = {key: D}            # keep one dataset per process
    return D


# --------------------------------------------------------------------------------- search
def inner_scores(folds, terms, lams=P.LAMS):
    """mean over the inner held-out divisions of r_unit, for each lambda."""
    R = np.zeros((len(folds), len(lams)))
    for a, f in enumerate(folds):
        for b, th in enumerate(f.solve(terms, lams)):
            R[a, b] = P.r_division(f.test, f.predict(terms, th))
    return R.mean(axis=0), R


def greedy(folds, cands):
    """forward selection; -> (terms, lambda, inner score, path, step-1 table)."""
    sel, best, lam = [], -np.inf, None
    path, step1 = [], {}
    for step in range(MAX_TERMS):
        rows = []
        for t in cands:
            if t in sel:
                continue
            m, _ = inner_scores(folds, sel + [t])
            j = int(np.argmax(m))
            rows.append((float(m[j]), t, P.LAMS[j]))
            if step == 0:
                step1[t] = (float(m[j]), P.LAMS[j])
        rows.sort(key=lambda x: -x[0])
        sc, t, l = rows[0]
        if step > 0 and sc < best + ADD_MIN:
            break
        sel.append(t)
        best, lam = sc, l
        path.append((list(sel), sc, l))
    return sel, lam, best, path, step1


def outer_task(args):
    """one outer fold of one dataset: the nested search, the comparators and ablations."""
    spec, d, full = args
    t0 = time.time()
    G = ctx()
    eng = G["eng"]
    S = dataset(spec)
    fold = [(dd, tr, te) for dd, tr, te in P.folds(S) if dd == d][0]
    _, tr, te = fold
    saved_scales = (eng.unit_scale, eng.pair_scale)
    if spec.get("foldscale"):
        # audit item 11 sensitivity: label-free scale constants from the outer training sources only
        eng.unit_scale, eng.pair_scale = eng.scales_for(tr)
    try:
        return _outer_task(spec, d, full, t0, G, eng, S, tr, te)
    finally:
        eng.unit_scale, eng.pair_scale = saved_scales


def _outer_task(spec, d, full, t0, G, eng, S, tr, te):
    inner = [P.Fold(eng, itr, ite) for _, itr, ite in P.folds(tr)]
    sel, lam, inner_best, path, step1 = greedy(inner, G["terms"])
    fo = P.Fold(eng, tr, te)
    res = {"spec": spec, "division": d, "terms": [list(t) for t in sel], "lambda": lam, "inner": inner_best,
           "path": [[[list(t) for t in p], s, l] for p, s, l in path]}
    th = fo.solve(sel, [lam])[0]
    pred = fo.predict(sel, th)
    res["r_search"] = P.r_division(te, pred)
    pD = fo.predict([D_TERM], fo.solve([D_TERM], [0.3])[0])
    res["r_D"] = P.r_division(te, pD)
    # D with lambda tuned in the inner loop (the search restricted to D's recipe)
    dn = step1[D_TERM]
    pDn = fo.predict([D_TERM], fo.solve([D_TERM], [dn[1]])[0])
    res["r_Dnested"] = P.r_division(te, pDn)
    res["Dnested_lambda"] = dn[1]
    res["step1"] = [[list(t), v[0], v[1]] for t, v in sorted(step1.items(), key=lambda x: -x[1][0])]
    if spec["kind"] in ("planted", "plantedD"):
        # oracle: the true term(s) alone, lambda tuned in the same inner loop, scored on the held-out division
        ot = [tuple(PLANTED[spec["truth"]])] + ([D_TERM] if spec["kind"] == "plantedD" else [])
        mO, _ = inner_scores(inner, ot)
        lo = P.LAMS[int(np.argmax(mO))]
        res["r_oracle"] = P.r_division(te, fo.predict(ot, fo.solve(ot, [lo])[0]))
        res["oracle_lambda"] = lo
    if full:
        res["per_source"] = {k: [P.r_division({k: te[k]}, pred), P.r_division({k: te[k]}, pD)] for k in te}
        res["unit_cos"] = {"%s|%s" % k: [v, P.unit_cos(te, pD)[k]] for k, v in P.unit_cos(te, pred).items()}
        res["pred"] = {k: v.tolist() for k, v in pred.items()}
        # ablation: drop each selected term, same lambda, refit on the outer training data
        res["ablation"] = {}
        for t in sel:
            rest = [x for x in sel if x != t]
            if rest:
                res["ablation"][json.dumps(list(t))] = res["r_search"] - P.r_division(te, fo.predict(rest, fo.solve(rest, [lam])[0]))
        # comparators on the same fold: K, B (frozen d105 recipes), A (location), regularised regression
        # on all gene PCs with lambda tuned in the inner loop (R_nested)
        H = G.setdefault("H", b5.load_frozen()[3])
        for m in ("K", "B", "A"):
            res["r_" + m] = P.r_division(te, H.predict(m, H.fit(m, tr), te))
        KT = ("K", "Kraw", "onehot")
        lamsR = (0.1, 1.0, 10.0, 100.0, 1000.0)
        mR, _ = inner_scores(inner, [KT], lamsR)
        lr = lamsR[int(np.argmax(mR))]
        res["r_Rnested"] = P.r_division(te, fo.predict([KT], fo.solve([KT], [lr])[0]))
        res["Rnested_lambda"] = lr
        # uncertainty of the cross-fitted (inferred) profile: refits on bootstrap-resampled training sources
        rng = np.random.default_rng(SEED + 7 + list(sorted({s["division"] for s in S.values()})).index(d))
        bs = []
        keys = sorted(tr)
        for _ in range(20):
            pick = rng.integers(0, len(keys), len(keys))
            trb = {"%s#%d" % (keys[i], j): tr[keys[i]] for j, i in enumerate(pick)}
            fb = BootFold(eng, trb, te, fo)
            bs.append(fb.predict(sel, fb.solve(sel, [lam])[0]))
        res["pred_sd"] = {k: np.std([b_[k] for b_ in bs], axis=0, ddof=1).tolist() for k in te}
    res["sec"] = time.time() - t0
    return res


class BootFold(P.Fold):
    """a Fold whose training sources are renamed copies (bootstrap); pieces come from the source name."""

    def __init__(self, eng, train, test, parent):
        self.eng = eng
        self.train = train
        self.test = test
        self.ytr = {k: P.src_obs(s) for k, s in train.items()}
        self.ctx = parent.ctx
        self.cache = {}
        self.parent = parent

    def piece(self, term, k, side):
        if side == "tr":
            key = (term, k, side)
            if key not in self.cache:
                base = k.split("#")[0]
                s = self.train[k]
                if term[0] == "K":
                    self.cache[key] = (self.eng.unit_X(term[1], base, s["units"], self.ctx), self.eng.enc_Q(term[2], s["targets"]))
                else:
                    self.cache[key] = self.eng.pair_D(term[1], base, s["units"], s["targets"])
            return self.cache[key]
        return self.parent.piece(term, k, side)


def run_specs(specs, full=False, label=""):
    """all outer folds of all specs over a process pool; -> {spec json: [7 fold results]}."""
    from multiprocessing import get_context
    G = ctx()
    divs = sorted({s["division"] for s in G["S"].values()})
    jobs = [(sp, d, full) for sp in specs for d in divs]
    out = {}
    t0 = time.time()
    with get_context("spawn").Pool(WORKERS) as pool:
        for i, r in enumerate(pool.imap_unordered(outer_task, jobs)):
            out.setdefault(json.dumps(r["spec"], sort_keys=True), []).append(r)
            if (i + 1) % 7 == 0 or i + 1 == len(jobs):
                say("%s: %d/%d outer folds done  [%.0fs]" % (label, i + 1, len(jobs), time.time() - t0))
    for k in out:
        out[k].sort(key=lambda r: r["division"])
    return out


def summary(folds):
    rs = np.array([f["r_search"] for f in folds])
    rd = np.array([f["r_D"] for f in folds])
    return float(rs.mean()), float(rd.mean()), float((rs - rd).mean())


# ---------------------------------------------------------------------------------- modes
def save(name, obj):
    p = os.path.join(OUT, "d108_" + name + ".json")
    json.dump(obj, open(p, "w"))
    say("wrote", os.path.relpath(p, ROOT))


def load(name):
    return json.load(open(os.path.join(OUT, "d108_" + name + ".json")))


def pilot():
    """timing: one full nested search on shuffled labels (N0 replicate 999), all 7 outer folds."""
    G = ctx()
    t0 = time.time()
    res = run_specs([{"kind": "N0", "rep": 999}], full=False, label="pilot")
    wall = time.time() - t0
    folds = list(res.values())[0]
    cpu = sum(f["sec"] for f in folds)
    s, d, dl = summary(folds)
    say("pilot N0: search %+.3f, D %+.3f; wall %.0fs, summed fold CPU %.0fs, workers %d" % (s, d, wall, cpu, WORKERS))
    # calibration of the simulations (planted betas and the D-as-truth beta)
    S = G["S"]
    cal = {"pilot": {"wall_s": wall, "cpu_s": cpu, "search": s, "D": d, "workers": WORKERS,
                     "terms_per_fold": [f["terms"] for f in folds]}}
    rD = d_lodo(S)
    cal["real_D_lodo"] = rD
    truth = d_truth(S)
    grid = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]
    rows = []
    for bt in grid:
        rr = []
        for rep in range(3):
            rng = np.random.default_rng(SEED + 9000 + rep)
            rr.append(d_lodo(with_counts(S, sample(S, d2_logits(S, truth, bt, rng), rng))))
        rows.append((bt, float(np.mean(rr))))
        say("N2 calibration: beta %.2f -> mean D LODO r %+.3f" % rows[-1])
    xs, ys = np.array([r[0] for r in rows]), np.array([r[1] for r in rows])
    beta2 = float(np.interp(rD, ys, xs)) if ys[0] < rD < ys[-1] else float(xs[np.argmin(np.abs(ys - rD))])
    cal["N2_grid"] = rows
    cal["N2_beta"] = beta2
    cal["planted_beta"] = {}
    for kind in PLANTED:
        bt = calibrate_planted(S, kind)
        cal["planted_beta"][kind] = bt
        say("planted %s: beta %.3f gives true r %.3f" % (kind, bt, true_r(S, planted_signal(S, kind), bt)))
    save("pilot", cal)


def declared():
    return json.load(open(DECLARED))


def run_planted():
    cal = load("pilot")
    dec = declared()
    specs = [{"kind": "planted", "truth": k, "beta": cal["planted_beta"][k], "rep": r}
             for k in PLANTED_DECLARED for r in range(dec["n_planted_per_truth"])]
    save("planted", run_specs(specs, label="planted"))


def run_null0():
    dec = declared()
    save("null0", run_specs([{"kind": "N0", "rep": r} for r in range(dec["n_N0"])], label="N0"))


def run_null2():
    cal = load("pilot")
    dec = declared()
    save("null2", run_specs([{"kind": "N2", "beta": cal["N2_beta"], "rep": r} for r in range(dec["n_N2"])], label="N2"))


def run_real():
    marker = os.path.join(OUT, "d108_real_opened.json")
    if os.path.exists(marker) and not os.environ.get("D108_RERUN"):
        sys.exit("the real-label nested search has already been run (%s)" % marker)
    json.dump({"time": time.strftime("%Y-%m-%d %H:%M:%S")}, open(marker, "w"))
    save("real", run_specs([{"kind": "real"}], full=True, label="real"))


SUPP = os.path.join(OUT, "d108_declared_supp.json")


def supp_declare():
    """post-audit supplementary declaration, written BEFORE the supplementary runs: betas of the
    added planted truths, replicate counts and the rules they are judged by."""
    if os.path.exists(SUPP):
        sys.exit("supplementary declaration already written (%s)" % SUPP)
    G = ctx()
    S = G["S"]
    dec = declared()
    t0 = time.time()
    beta = {}
    for kind, r in SUPP_ALONE + SUPP_DPLUS:
        key = "%s@%.2f" % (kind, r)
        if key not in beta:
            beta[key] = calibrate_planted(S, kind, r)
            say("supp %s: route beta %.3f gives route-alone true r %.3f  [%.0fs]" % (
                key, beta[key], true_r(S, planted_signal(S, kind), beta[key]), time.time() - t0))
    out = {"declared_at": time.strftime("%Y-%m-%d %H:%M:%S"), "workers": WORKERS,
           "reason": "post-audit additions (audit items 2 and 11); the main declaration and its rules are unchanged",
           "n_per_truth": 8,
           "alone": [{"truth": k, "true_r": r, "beta": beta["%s@%.2f" % (k, r)]} for k, r in SUPP_ALONE],
           "D_plus_route": [{"truth": k, "route_true_r_alone": r, "beta": beta["%s@%.2f" % (k, r)], "betaD": dec["N2_beta"]}
                            for k, r in SUPP_DPLUS],
           "rerun_declared_planted": "the 3 declared planted truths are rerun unchanged (same seeds) to add the oracle score; their terms must reproduce exactly",
           "rules": {
               "oracle": "oracle = the true term(s) alone (for D + route: D's term + the route term), lambda tuned in the same inner loop, scored on the held-out division",
               "recovery_extended": "planted term selected in >= 4 of 7 outer folds AND search mean7 >= oracle mean7 - %.2f (the expected-range clause)" % ORACLE_MARGIN,
               "D_plus_route_power": "reported per truth: fraction of replicates with the route recovered (extended rule); fraction with search mean7 > D mean7; fraction with exact division sign-flip p(search > D) < 0.05; mean oracle gain over D",
               "foldscale": "one nested real-label search with the label-free scale constants (Lm, rff_Lm, gD4xLm, distLm, distLhat) recomputed from the outer training sources only; reported as a sensitivity, the main result stays the declared run"}}
    json.dump(out, open(SUPP, "w"), indent=1)
    say(json.dumps(out, indent=1))


def run_planted2():
    sup = json.load(open(SUPP))
    cal = load("pilot")
    n = sup["n_per_truth"]
    specs = [{"kind": "planted", "truth": k, "beta": cal["planted_beta"][k], "rep": r} for k in PLANTED_DECLARED for r in range(8)]
    specs += [{"kind": "planted", "truth": a["truth"], "beta": a["beta"], "rep": r} for a in sup["alone"] for r in range(n)]
    specs += [{"kind": "plantedD", "truth": a["truth"], "r": a["route_true_r_alone"], "beta": a["beta"], "betaD": a["betaD"], "rep": r}
              for a in sup["D_plus_route"] for r in range(n)]
    res = run_specs(specs, label="planted2")
    old = load("planted")
    for k, v in old.items():
        new = res[k]
        assert [f["terms"] for f in new] == [f["terms"] for f in v], k
        assert np.allclose([f["r_search"] for f in new], [f["r_search"] for f in v], atol=1e-12), k
    say("declared planted runs reproduced exactly (terms and outer r)")
    save("planted2", res)


def run_foldscale():
    save("real_foldscale", run_specs([{"kind": "real", "foldscale": True}], full=False, label="foldscale"))


def extras():
    """O2 comparator per division; nested out-of-division predictions with an honest sd; module-card
    coverage of the atlas modules."""
    G = ctx()
    S = G["S"]
    H = b5.load_frozen()[3]
    o2 = {}
    for d, tr, te in P.folds(S):
        o2[d] = P.r_division(te, H.predict("O2", H.fit("O2", tr), te))
    say("O2 per division:", {k: round(v, 3) for k, v in o2.items()}, "mean7 %+.3f" % np.mean(list(o2.values())))
    real = load("real")[json.dumps({"kind": "real"}, sort_keys=True)]
    # residuals of the nested predictions on each held-out division (observed double-centred profile - prediction)
    res = {}
    for fo in real:
        for k, M in fo["pred"].items():
            Y = P.src_obs(S[k])[0]
            res.setdefault(fo["division"], []).append((np.asarray(Y) - np.asarray(M)).ravel())
    rv = {d: float(np.mean(np.concatenate(v) ** 2)) for d, v in res.items()}
    inf, cover, cover_boot = [], [], []
    for fo in real:
        d = fo["division"]
        other = float(np.mean([rv[x] for x in rv if x != d]))     # residual variance from the OTHER 6 held-out divisions
        for k, M in fo["pred"].items():
            Y = P.src_obs(S[k])[0]
            sdb = np.asarray(fo["pred_sd"][k])
            sdt = np.sqrt(sdb ** 2 + other)
            cover.append((np.abs(Y - np.asarray(M)) <= 1.96 * sdt).ravel())
            cover_boot.append((np.abs(Y - np.asarray(M)) <= 1.96 * sdb).ravel())
            for a, u in enumerate(S[k]["units"]):
                for b, t in enumerate(S[k]["targets"]):
                    inf.append({"source": k, "subclass": u, "target": t, "value": M[a][b], "sd": float(sdt[a, b]),
                                "sd_refit_only": float(sdb[a, b]), "division_held_out": d,
                                "provenance": "I (inferred: nested search, out-of-division; pipeline selected inside that outer fold)",
                                "predictive": False,
                                "heldout_r_unit_of_this_division": fo["r_search"]})
    cover = float(np.concatenate(cover).mean())
    cover_boot = float(np.concatenate(cover_boot).mean())
    old = os.path.join(P.MODDIR, "S1_inferred_profiles.json")
    if os.path.exists(old):
        os.remove(old)
    newp = os.path.join(P.MODDIR, "nested_search_oof_predictions.json")
    json.dump({"description": "out-of-division predictions of the NESTED search: each value comes from the pipeline selected "
                              "inside the outer fold that held its division out (7 different pipelines), NOT from S1. "
                              "Held-out mean7 r_unit -0.062: NOT PREDICTIVE, do not use in place of D.",
               "sd": "sqrt(variance over 20 bootstrap refits of that fold's pipeline (selection and E-hat fixed) + mean squared "
                     "out-of-division residual of the nested predictions on the OTHER 6 divisions (covers selection and residual error))",
               "coverage_1.96sd": {"sd": cover, "sd_refit_only": cover_boot},
               "residual_var_per_division": rv, "records": inf}, open(newp, "w"))
    say("nested predictions: %d records; 1.96 sd coverage %.3f (refit-only sd %.3f)" % (len(inf), cover, cover_boot))
    # atlas module cards: coverage of pred_sd on the out-of-fold errors
    cp = os.path.join(P.MODDIR, "atlas_module_cards.json")
    cards = json.load(open(cp))
    cov = {}
    for c in cards:
        z = np.load(os.path.join(ROOT, c["file"]), allow_pickle=True)
        err = np.abs(z["Y"] - z["pred_oof"])
        rsd = (z["Y"] - z["pred_oof"]).std(0, ddof=1)[None, :]
        cv = {"coverage_1.96_pred_sd": float((err <= 1.96 * z["pred_sd"]).mean()),
              "coverage_1.96_resid_sd_only": float((err <= 1.96 * rsd).mean()),
              "calibration_factor_for_95pct": float(np.quantile(err / z["pred_sd"], 0.95) / 1.96)}
        c["uncertainty_calibration"] = dict(cv, note="pred_sd is conservative: bootstrap variance plus the out-of-fold residual "
                                                     "variance partly double-count; multiply pred_sd by calibration_factor_for_95pct "
                                                     "for nominal 95% intervals on in-atlas subclasses (factor measured on the same out-of-fold errors)")
        cov[c["id"]] = cv
    json.dump(cards, open(cp, "w"), indent=1)
    say("atlas module coverage:", cov)
    save("extras", {"O2_per_division": o2, "O2_mean7": float(np.mean(list(o2.values()))),
                    "nested_pred_coverage": cover, "nested_pred_coverage_refit_only": cover_boot,
                    "nested_resid_var": rv, "atlas_coverage": cov})


def signflip(d, exact_max=20, nmc=200000, seed=0):
    """one-sided paired sign-flip p for mean(d) > 0 (exact when len(d) <= exact_max)."""
    d = np.asarray(d, float)
    obs = d.mean()
    n = len(d)
    if n <= exact_max:
        import itertools
        S = np.array(list(itertools.product((1, -1), repeat=n)))
        return float(((S @ d) / n >= obs - 1e-12).mean())
    rng = np.random.default_rng(seed)
    S = rng.choice((1.0, -1.0), size=(nmc, n))
    return float((1 + ((S @ d) / n >= obs - 1e-12).sum()) / (1 + nmc))


def report():
    """all numbers for reports/d108_phase1.md -> data/interim/d108_summary.json (+ printed)."""
    import d108_symbolic as SR
    dec = declared()
    out = {"declared": dec}
    real = load("real")[json.dumps({"kind": "real"}, sort_keys=True)]
    divs = [f["division"] for f in real]
    cols = ("r_search", "r_D", "r_Dnested", "r_K", "r_B", "r_A", "r_Rnested")
    tab = {c: [f[c] for f in real] for c in cols}
    out["real_per_division"] = {"divisions": divs, **tab}
    out["real_mean7"] = {c: float(np.mean(v)) for c, v in tab.items()}
    rs = np.array(tab["r_search"])
    delta = float(rs.mean() - np.mean(tab["r_D"]))
    out["delta_search_minus_D"] = delta
    tests = {}
    for c in cols[1:]:
        dd = rs - np.array(tab[c])
        ps = [v[0] - v[1] for f in real for v in f["per_source"].values()] if c == "r_D" else None
        tests[c] = {"mean_delta": float(dd.mean()), "p_division": signflip(dd), "n_div_search_better": int((dd > 0).sum())}
    # source and subclass level (vs D)
    src_d = [v[0] - v[1] for f in real for v in f["per_source"].values()]
    unit_d = [v[0] - v[1] for f in real for v in f["unit_cos"].values()]
    tests["r_D"].update({"n_sources": len(src_d), "mean_delta_source": float(np.mean(src_d)), "p_source": signflip(src_d, seed=1),
                         "n_units": len(unit_d), "mean_delta_unit_cos": float(np.mean(unit_d)), "p_unit": signflip(unit_d, seed=2)})
    out["tests_search_vs"] = tests
    out["real_folds"] = [{"division": f["division"], "terms": f["terms"], "lambda": f["lambda"], "inner": f["inner"],
                          "path": f["path"], "ablation": f["ablation"], "Dnested_lambda": f["Dnested_lambda"],
                          "Rnested_lambda": f["Rnested_lambda"], "D_step1_rank": [json.dumps(x[0]) for x in f["step1"]].index(json.dumps(list(D_TERM))) + 1,
                          "step1_top3": f["step1"][:3], "sec": f["sec"]} for f in real]
    cnt = {}
    abl = {}
    for f in real:
        for t in f["terms"]:
            cnt[json.dumps(t)] = cnt.get(json.dumps(t), 0) + 1
        for t, v in f["ablation"].items():
            abl.setdefault(t, []).append(v)
        if len(f["terms"]) == 1:      # post-audit: a lone term's ablation = r(term) - r(no prediction, 0)
            abl.setdefault(json.dumps(f["terms"][0]), []).append(f["r_search"])
    # N0
    n0 = load("null0")
    n0_mean = [float(np.mean([f["r_search"] for f in v])) for v in n0.values()]
    n0_D = [float(np.mean([f["r_D"] for f in v])) for v in n0.values()]
    n0_maxcnt = []
    for v in n0.values():
        c = {}
        for f in v:
            for t in f["terms"]:
                c[json.dumps(t)] = c.get(json.dumps(t), 0) + 1
        n0_maxcnt.append(max(c.values()))
    thr = float(np.percentile(n0_maxcnt, 95))
    out["N0"] = {"n": len(n0_mean), "search_mean7": [float(np.mean(n0_mean)), float(np.percentile(n0_mean, 95)), float(np.max(n0_mean))],
                 "D_mean7_mean": float(np.mean(n0_D)), "p_N0": float((1 + sum(x >= rs.mean() for x in n0_mean)) / (1 + len(n0_mean))),
                 "maxcount_95pct": thr, "maxcount_dist": np.bincount(n0_maxcnt, minlength=8).tolist(),
                 "sec_total": float(sum(f["sec"] for v in n0.values() for f in v))}
    out["selection"] = {t: {"count": c, "ablation_mean": float(np.mean(abl[t])) if t in abl else None,
                            "ablation_n": len(abl.get(t, [])),
                            "discovered_route": bool(c >= 4 and c > thr and t in abl and np.mean(abl[t]) > 0)}
                        for t, c in sorted(cnt.items(), key=lambda x: -x[1])}
    # N2
    n2 = load("null2")
    n2_delta = [summary(v)[2] for v in n2.values()]
    n2_D = [summary(v)[1] for v in n2.values()]
    n2_S = [summary(v)[0] for v in n2.values()]
    p_n2 = float((1 + sum(x >= delta for x in n2_delta)) / (1 + len(n2_delta)))
    out["N2"] = {"n": len(n2_delta), "delta_mean": float(np.mean(n2_delta)), "delta_95pct": float(np.percentile(n2_delta, 95)),
                 "delta_sd": float(np.std(n2_delta, ddof=1)), "D_mean7_mean": float(np.mean(n2_D)), "search_mean7_mean": float(np.mean(n2_S)),
                 "frac_delta_gt0": float(np.mean(np.array(n2_delta) > 0)), "p_N2": p_n2,
                 "frac_D_term_selected_first": float(np.mean([f["terms"][0] == list(D_TERM) for v in n2.values() for f in v])),
                 "frac_folds_D_term_selected": float(np.mean([list(D_TERM) in f["terms"] for v in n2.values() for f in v])),
                 "sec_total": float(sum(f["sec"] for v in n2.values() for f in v))}
    p_div = tests["r_D"]["p_division"]
    if delta <= 0:
        verdict = "NO"
    elif delta > out["N2"]["delta_95pct"] and p_n2 < 0.05 and p_div < 0.05:
        verdict = "YES"
    else:
        verdict = "NOT SHOWN"
    out["verdict_search_beats_D"] = verdict
    # planted
    pl = load("planted")
    prow = {}
    for k, v in pl.items():
        sp = json.loads(k)
        term = json.dumps(list(PLANTED[sp["truth"]]))
        c = sum(json.dumps(t) == term or t == list(PLANTED[sp["truth"]]) for f in v for t in f["terms"])
        s_, d_, dl = summary(v)
        prow.setdefault(sp["truth"], []).append({"rep": sp["rep"], "count": c, "recovered": c >= 4, "search": s_, "D": d_, "delta": dl})
    out["planted"] = {t: {"n": len(r), "frac_recovered": float(np.mean([x["recovered"] for x in r])),
                          "counts": [x["count"] for x in sorted(r, key=lambda y: y["rep"])],
                          "search_mean": float(np.mean([x["search"] for x in r])), "D_mean": float(np.mean([x["D"] for x in r])),
                          "delta_mean": float(np.mean([x["delta"] for x in r])), "p_delta_reps": signflip([x["delta"] for x in r])}
                      for t, r in prow.items()}
    out["planted_underpowered"] = all(v["frac_recovered"] < 0.5 for v in out["planted"].values())
    # symbolic
    sr = load("sr_real")[json.dumps({"kind": "real"}, sort_keys=True)]
    srn = load("sr_null0")
    fr = np.array([f["r_formula"] for f in sr])
    lf = np.array([f["r_leaf"] for f in sr])
    srn_mean = [float(np.mean([f["r_formula"] for f in v])) for v in srn.values()]
    p_sr = signflip(fr - lf)
    out["SR"] = {"catalogue_size": len(SR.catalogue()[0]), "PEN": SR.PEN,
                 "folds": [{k: f[k] for k in ("division", "formula", "complexity", "inner", "r_formula", "leaf", "r_leaf", "r_ld", "inner_top5", "sec")} for f in sr],
                 "formula_mean7": float(fr.mean()), "leaf_mean7": float(lf.mean()), "ld_mean7": float(np.mean([f["r_ld"] for f in sr])),
                 "p_division_formula_vs_leaf": p_sr, "N0_n": len(srn_mean), "N0_95pct": float(np.percentile(srn_mean, 95)),
                 "N0_mean": float(np.mean(srn_mean)), "p_N0": float((1 + sum(x >= fr.mean() for x in srn_mean)) / (1 + len(srn_mean))),
                 "N0_sec_total": float(sum(f["sec"] for v in srn.values() for f in v))}
    out["SR"]["beats_simpler"] = bool(fr.mean() > lf.mean() and p_sr < 0.05 and fr.mean() > out["SR"]["N0_95pct"])
    out["SR_final"] = load("sr_final")
    out["final"] = load("final")
    out["search_space"] = {"n_terms": len(ctx()["terms"]), "unit_nodes": list(P.Engine.UNIT_NODES), "encodings": list(P.ENCODINGS),
                           "pairs": list(P.Engine.PAIRS), "lambdas": list(P.LAMS), "max_terms": MAX_TERMS,
                           "pipelines_up_to_max_terms": int(sum(__import__("math").comb(len(ctx()["terms"]), k) for k in range(1, MAX_TERMS + 1)))}
    out["audit"] = audit_numbers(real, n0, n2, rs.mean(), delta)
    save("summary", out)
    say(json.dumps({k: v for k, v in out.items() if k not in ("declared", "real_folds", "SR", "SR_final")}, indent=1))


T1_REHEARSAL = {"amygdala-pallidum": 0.227, "hindbrain": 0.071, "hippocampal": 0.343, "isocortex": 0.479,
                "midbrain": 0.296, "olfactory": 0.293, "thalamus": 0.251}   # reports/d105_rehearsal.md (same folds; not rerun)


def audit_numbers(real, n0, n2, real_mean, delta):
    """post-audit numbers: selection optimism, extra comparators, supplementary planted truths, fold-scale sensitivity."""
    a = {}
    Dk = list(D_TERM)

    def d_inner(f):
        return [x[1] for x in f["step1"] if x[0] == Dk][0]
    for name, sets in (("N0", n0), ("N2", n2)):
        fl = [f for v in sets.values() for f in v]
        a["optimism_" + name] = {"selected_inner_mean": float(np.mean([f["inner"] for f in fl])),
                                 "D_term_inner_mean": float(np.mean([d_inner(f) for f in fl])),
                                 "selected_outer_mean": float(np.mean([f["r_search"] for f in fl])),
                                 "frac_folds_D_term_selected": float(np.mean([Dk in f["terms"] for f in fl]))}
    a["optimism_real"] = {"selected_inner_mean": float(np.mean([f["inner"] for f in real])),
                          "D_term_inner": [d_inner(f) for f in real],
                          "selected_outer_mean": float(real_mean)}
    n2d = np.array([summary(v)[2] for v in n2.values()])
    a["N2_delta_min"] = float(n2d.min())
    a["N2_n_delta_below_observed"] = int((n2d < delta).sum())
    rs = np.array([f["r_search"] for f in real])
    divs = [f["division"] for f in real]
    ex = load("extras")
    for m, vals in (("O2", [ex["O2_per_division"][d] for d in divs]), ("T1_rehearsal", [T1_REHEARSAL[d] for d in divs])):
        dd = rs - np.array(vals)
        a["vs_" + m] = {"per_division": vals, "mean7": float(np.mean(vals)), "mean_delta": float(dd.mean()),
                        "p_division": signflip(dd), "n_div_search_better": int((dd > 0).sum())}
    a["extras"] = ex
    # supplementary planted truths
    sup = json.load(open(SUPP))
    pl2 = load("planted2")
    rows = {}
    for k, v in pl2.items():
        sp = json.loads(k)
        name = sp["truth"] + ("" if sp["kind"] == "planted" else " (D + route, route r %.2f)" % sp["r"])
        if sp["kind"] == "planted" and sp["truth"] == "Ehat_tc":
            name = "Ehat_tc (alone, r 0.40)"
        term = list(PLANTED[sp["truth"]])
        c = sum(t == term for f in v for t in f["terms"])
        s_, d_, dl = summary(v)
        orc = float(np.mean([f["r_oracle"] for f in v]))
        dv = np.array([f["r_search"] - f["r_D"] for f in v])
        rows.setdefault(name, []).append({"rep": sp["rep"], "count": c, "rec_declared": c >= 4,
                                          "rec_extended": bool(c >= 4 and s_ >= orc - ORACLE_MARGIN),
                                          "search": s_, "D": d_, "oracle": orc, "delta": dl,
                                          "p_div": signflip(dv), "D_term_count": sum(Dk in f["terms"] for f in v)})
    a["planted2"] = {n: {"n": len(r), "counts": [x["count"] for x in sorted(r, key=lambda y: y["rep"])],
                         "frac_rec_declared": float(np.mean([x["rec_declared"] for x in r])),
                         "frac_rec_extended": float(np.mean([x["rec_extended"] for x in r])),
                         "search_mean": float(np.mean([x["search"] for x in r])), "D_mean": float(np.mean([x["D"] for x in r])),
                         "oracle_mean": float(np.mean([x["oracle"] for x in r])),
                         "oracle_gain_over_D": float(np.mean([x["oracle"] - x["D"] for x in r])),
                         "delta_mean": float(np.mean([x["delta"] for x in r])),
                         "frac_search_gt_D": float(np.mean([x["delta"] > 0 for x in r])),
                         "frac_p_div_lt_05": float(np.mean([x["p_div"] < 0.05 for x in r])),
                         "p_delta_reps": signflip([x["delta"] for x in r]),
                         "mean_D_term_folds": float(np.mean([x["D_term_count"] for x in r]))} for n, r in rows.items()}
    a["supp_declared"] = sup
    fs = load("real_foldscale")[json.dumps({"foldscale": True, "kind": "real"}, sort_keys=True)]
    a["foldscale"] = {"per_division": [f["r_search"] for f in fs], "mean7": float(np.mean([f["r_search"] for f in fs])),
                      "same_terms_as_declared_run": [f["terms"] == g["terms"] for f, g in zip(fs, real)],
                      "max_abs_diff_r": float(np.max(np.abs(np.array([f["r_search"] for f in fs]) - rs))),
                      "sec": float(sum(f["sec"] for f in fs))}
    a["planted2_sec"] = float(sum(f["sec"] for v in pl2.values() for f in v))
    return a


def run_final():
    """the reusable discovered module: the same greedy search, now with all 27 sources and the 7
    leave-one-division-out folds as the selection folds; the selected pipeline is fitted on all
    sources; uncertainty = refits on 200 bootstrap resamples of sources. Its held-out accuracy is the
    NESTED outer score from the real run (d108_real.json), never the score of this selection."""
    G = ctx()
    eng, S = G["eng"], G["S"]
    t0 = time.time()
    sel_folds = [P.Fold(eng, tr, te) for _, tr, te in P.folds(S)]
    sel, lam, best, path, step1 = greedy(sel_folds, G["terms"])
    say("final selection on all 27 sources: %s lambda %s, LODO-7 selection score %.3f (optimistic: selected on it)  [%.0fs]"
        % (sel, lam, best, time.time() - t0))
    fa = P.Fold(eng, S, {})          # no held-out side: E-hat cross-fitted over all training subclasses
    th = fa.solve(sel, [lam])[0]
    rng = np.random.default_rng(SEED + 77)
    keys = sorted(S)
    ths = []
    for _ in range(200):
        pick = rng.integers(0, len(keys), len(keys))
        trb = {"%s#%d" % (keys[i], j): S[keys[i]] for j, i in enumerate(pick)}
        fb = BootFold(eng, trb, S, fa)
        ths.append(fb.solve(sel, [lam])[0])
    ths = np.array(ths)
    real = load("real")[json.dumps({"kind": "real"}, sort_keys=True)]
    rs = [f["r_search"] for f in real]
    os.makedirs(P.MODDIR, exist_ok=True)
    f = os.path.join(P.MODDIR, "S1_discovered_pipeline.npz")
    units = [[k, u] for k in keys for u in S[k]["units"]]
    # self-contained bundle (audit item 5): the feature bases and target encodings the terms use
    bundle = {}
    for t in sel:
        if t[0] != "K" or t[1] not in P.GENE_NODES:
            raise SystemExit("S1 bundle supports gene-node x encoding terms only: %s" % (t,))
        bundle["F_" + t[1]] = eng.F[t[1]]
        bundle["enc_" + t[2]] = eng.enc[t[2]]
    np.savez(f, terms=np.array(json.dumps([list(t) for t in sel])), lam=lam, theta=th, theta_boot=ths,
             theta_boot_sd=ths.std(0, ddof=1), train_units=np.array(["%s|%s" % tuple(x) for x in units]),
             names=np.array(eng.names), targets=np.array(eng.targets), provenance=np.array("I"), **bundle)
    # check: the bundle alone reproduces the engine's in-sample prediction for every source
    pin = fa.predict(sel, th, side="tr")
    err = max(float(np.abs(P.apply_bundle(f, S[k]["units"], S[k]["targets"]) - pin[k]).max()) for k in S)
    assert err < 1e-9, err
    say("S1 bundle reproduces the engine prediction (max abs diff %.1e)" % err)
    card = {"id": "S1_discovered_pipeline", "file": os.path.relpath(f, ROOT),
            "description": "projection-profile model selected by the D-108 search (greedy additive selection from a fixed typed library of 66 terms) on all 27 non-HY sources",
            "terms": [list(t) for t in sel], "lambda": lam, "selection_path": [[[list(t) for t in p], s, l] for p, s, l in path],
            "training_data": "Epi-retro-seq st2, 27 non-HY sources, HY19 removed; %d (source, subclass) units" % len(units),
            "heldout_accuracy": {"nested_outer_r_unit_per_division": dict(zip([x["division"] for x in real], rs)),
                                 "mean7": float(np.mean(rs)),
                                 "note": "nested: the search was rerun inside each outer fold; the selection score on all 27 sources (%.3f) is optimistic" % best},
            "uncertainty": "theta_boot: 200 refits on bootstrap resamples of the 27 sources (selection fixed, so selection variance is NOT included; the nested held-out r of -0.062 is the honest accuracy)",
            "usage": "analysis/d108_prims.apply_bundle(npz, units, targets) for one source: X = F_<node>[units] centred over the "
                     "source's units; Q = enc_<encoding>[targets] centred over the source's targets; prediction = sum over terms of "
                     "X @ theta_term.reshape(dim X, dim Q) @ Q.T (double-centred log-share profile). Units must be Allen subclasses in "
                     "`names` (feature rows come from the atlas: D's whole-brain PCA basis, NMF and autoencoder fitted on all 5,322 clusters).",
            "not_predictive": "nested held-out mean7 -0.062 is below D (+0.228); do not use S1 in place of D",
            "flag": "every output is INFERRED (provenance I)", "sec": time.time() - t0}
    cards_p = os.path.join(P.MODDIR, "discovered_module_cards.json")
    json.dump([card], open(cards_p, "w"), indent=1)
    save("final", {"terms": card["terms"], "lambda": lam, "selection_score": best, "path": card["selection_path"],
                   "step1": [[list(t), v[0], v[1]] for t, v in sorted(step1.items(), key=lambda x: -x[1][0])][:15],
                   "theta_boot_sd_mean": float(ths.std(0, ddof=1).mean()), "theta_abs_mean": float(np.abs(th).mean()),
                   "sec": time.time() - t0})


def main():
    modes = {"pilot": pilot, "planted": run_planted, "null0": run_null0, "null2": run_null2, "real": run_real,
             "final": run_final, "report": report, "supp_declare": supp_declare, "planted2": run_planted2,
             "foldscale": run_foldscale, "extras": extras}
    try:
        import d108_symbolic as SR
        modes.update(SR.MODES)
    except ImportError:
        pass
    if len(sys.argv) != 2 or sys.argv[1] not in modes:
        sys.exit("usage: python3 -I analysis/d108_discovery.py {%s}" % "|".join(modes))
    modes[sys.argv[1]]()


if __name__ == "__main__":
    main()
