"""D-108 Phase 1 primitives: the primitive library, the reusable conversion modules and the fast
ridge engine used by analysis/d108_discovery.py.

LEAKAGE RULE (same as D-105, enforced here)
  Projection data come only from d105_benchmark.read_st2_training(), which drops every
  Supplementary Table 2 row with Source == "HY" right after the sheet is read; sources() then
  removes the 19 hypothalamic test subclasses (HY19) from every training source. Nothing in D-108
  reads a hypothalamic-source row. Unlabelled atlas data (Allen 10x cluster means, MERFISH cells
  and CCF coordinates, hypothalamic atlas cells included) may be used in unsupervised, target-side
  and unit-placement steps, the same way for D and for the search (D's whole-brain PCA was fitted
  on all 5,322 clusters, so every atlas transform here is too).
"""
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import d105_benchmark as b5  # noqa: E402

ROOT = b5.ROOT
INTERIM = os.path.join(ROOT, "data/interim")
PREP = os.path.join(INTERIM, "d108_prep.npz")
MODDIR = os.path.join(INTERIM, "d108_modules")
SEED = 20261110

GENE_NODES = ("gD32", "gD8", "gK32", "gB10", "nmf16", "ae16")
ENCODINGS = ("onehot", "tc", "tg", "tcomp")
LAMS = (0.1, 0.3, 1.0, 3.0, 10.0)
A_MIN_MM = 0.05


def say(*a):
    print(*a, flush=True)


def norm_rows(M):
    M = np.asarray(M, float)
    M = M - M.mean(axis=0)
    c = np.sqrt((M ** 2).sum(axis=1).mean())
    return M / (c if c > 0 else 1.0)


def rff(X, D, seed, bw=None):
    """random Fourier features of an RBF kernel, bandwidth = median pairwise distance."""
    X = np.asarray(X, float)
    if bw is None:
        d = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
        bw = float(np.median(d[np.triu_indices(len(X), 1)]))
    r = np.random.default_rng(seed)
    W = r.normal(size=(X.shape[1], D)) / bw
    bb = r.uniform(0, 2 * np.pi, D)
    return np.sqrt(2.0 / D) * np.cos(X @ W + bb), bw, W, bb


def poly2(X):
    k = X.shape[1]
    cols = [X] + [X[:, [a]] * X[:, [b]] for a in range(k) for b in range(a, k)]
    return np.hstack(cols)


# ------------------------------------------------------------------ the world (no HY rows)
def world():
    tab = b5.read_st2_training()
    if (tab["Source"] == "HY").any():
        sys.exit("leakage guard: HY rows in the training table")
    S = b5.sources(tab)
    assert "HY" not in S and len(S) == 27, sorted(S)
    for s in S.values():
        assert not set(s["units"]) & set(b5.HY19)
    targets = sorted({t for s in S.values() for t in s["targets"]})
    names = sorted({u for s in S.values() for u in s["units"]} | set(b5.HY19))
    return tab, S, targets, names


# ------------------------------------------------------------------------------ KRR modules
def krr_fit(X, Y, seed=0):
    """RBF kernel ridge; alpha and gamma by 5-fold CV inside the training data."""
    from sklearn.kernel_ridge import KernelRidge
    from sklearn.model_selection import GridSearchCV, KFold
    d = X.shape[1]
    grid = {"alpha": [0.01, 0.1, 1.0], "gamma": [0.1 / d, 0.3 / d, 1.0 / d]}
    g = GridSearchCV(KernelRidge(kernel="rbf"), grid, cv=KFold(5, shuffle=True, random_state=seed),
                     scoring="neg_mean_squared_error")
    g.fit(X, Y)
    return g.best_estimator_, g.best_params_


def krr_module(X, Y, names, label, nboot=30, seed=SEED):
    """A reusable nonlinear conversion: out-of-fold predictions (5-fold over subclasses, CV for
    hyper-parameters inside each training fold), held-out accuracy, a full fit, and a bootstrap
    (over subclasses) spread for every subclass's prediction."""
    from sklearn.kernel_ridge import KernelRidge
    from sklearn.model_selection import KFold
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1
    Xs = (X - mu) / sd
    ym = Y.mean(0)
    oof = np.zeros_like(Y)
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(Xs):
        m, _ = krr_fit(Xs[tr], Y[tr] - ym, seed)
        oof[te] = m.predict(Xs[te]) + ym
    full, hp = krr_fit(Xs, Y - ym, seed)
    rng = np.random.default_rng(seed + 7)
    boots = []
    for _ in range(nboot):
        ix = rng.integers(0, len(Xs), len(Xs))
        m = KernelRidge(kernel="rbf", **hp).fit(Xs[ix], Y[ix] - ym)
        boots.append(m.predict(Xs) + ym)
    bsd = np.std(boots, axis=0, ddof=1)
    resid = Y - oof
    r2 = 1 - (resid ** 2).sum(0) / ((Y - Y.mean(0)) ** 2).sum(0)
    base = Y - Y.mean(0)
    acc = {"r2_per_output_oof": [float(x) for x in r2],
           "r2_overall_oof": float(1 - (resid ** 2).sum() / (base ** 2).sum()),
           "rmse_oof": float(np.sqrt((resid ** 2).mean())), "rmse_mean_baseline": float(np.sqrt((base ** 2).mean())),
           "resid_sd_per_output": [float(x) for x in resid.std(0, ddof=1)]}
    pred_sd = np.sqrt(bsd ** 2 + resid.var(0, ddof=1)[None, :])
    return {"label": label, "names": list(names), "X": X, "Y": Y, "x_mu": mu, "x_sd": sd, "y_mean": ym,
            "hp": hp, "oof": oof, "full": full.predict(Xs) + ym, "boot_sd": bsd, "pred_sd": pred_sd,
            "acc": acc, "dual_coef": full.dual_coef_, "Xs": Xs, "model": full, "bootstrap_models": None}


def save_module(mod, mid, desc, inputs, outputs, extra=None):
    os.makedirs(MODDIR, exist_ok=True)
    f = os.path.join(MODDIR, mid + ".npz")
    np.savez(f, names=np.array(mod["names"]), X=mod["X"], Y=mod["Y"], x_mu=mod["x_mu"], x_sd=mod["x_sd"],
             y_mean=mod["y_mean"], train_Xs=mod["Xs"], dual_coef=mod["dual_coef"],
             hp=np.array(json.dumps(mod["hp"])), pred_oof=mod["oof"], pred_full=mod["full"],
             boot_sd=mod["boot_sd"], pred_sd=mod["pred_sd"], provenance=np.array("I"))
    card = {"id": mid, "file": os.path.relpath(f, ROOT), "description": desc, "inputs": inputs, "outputs": outputs,
            "kernel": "rbf kernel ridge", "hyperparameters": mod["hp"], "n_train": len(mod["names"]),
            "heldout_accuracy": mod["acc"],
            "uncertainty": "pred_sd = sqrt(bootstrap-over-subclasses variance (30 refits) + out-of-fold residual variance)",
            "usage": "standardise x with x_mu/x_sd; y = y_mean + K(x, train_Xs) @ dual_coef with K = exp(-gamma |x - x'|^2)",
            "flag": "every output of this module is INFERRED (provenance I)"}
    if extra:
        card.update(extra)
    return card


# --------------------------------------------------------------------------------- prep
def target_mask(mer, t):
    sel = None
    cbn = [x for x in mer.allst if x in ("FN", "IP", "DN", "CBN")]
    if t == "CBN" and cbn and np.isin(mer.st, cbn).sum() > 50:
        sel = np.isin(mer.st, cbn)
    elif t in b5.TGT_DIV:
        sel = np.isin(mer.div, b5.TGT_DIV[t])
    if t not in b5.TGT_DIV:
        sel = np.isin(mer.st, mer._structs(b5.TGT_STRUCT[t]))
    return sel


def prep():
    import torch
    torch.set_num_threads(2)
    t0 = time.time()
    tab, S, targets, names = world()
    z = np.load(os.path.join(INTERIM, "d105_models.npz"), allow_pickle=True)
    assert list(z["names"]) == names and list(z["targets"]) == targets
    say("world: %d sources, %d units, NAMES %d, targets %d  [%.0fs]" % (
        len(S), sum(len(s["units"]) for s in S.values()), len(names), len(targets), time.time() - t0))

    abc = b5.load_abc()
    X = abc["X"]
    keep = z["wb_keep"]
    Z = ((X[:, keep] - z["wb_mu"]) / z["wb_sd"]).astype(np.float32)
    Ec = (Z - z["wb_zc_mean"]) @ z["wb_loadings"]                      # D's whole-brain PCA scores, all clusters
    subn = abc["subnames"]
    nsub = len(subn)

    def sub_avg(M):
        out = np.zeros((nsub, M.shape[1]))
        for a in range(nsub):
            m = abc["subc"] == a
            out[a] = abc["w"][m] @ M[m] / abc["w"][m].sum()
        return out
    emb_all = sub_avg(Ec)
    ni = [subn.index(u) for u in names]
    assert np.allclose(emb_all[ni], z["wb_emb"], atol=1e-3), "D's PCA embedding not reproduced"
    say("whole-brain PCA embedding reproduced for NAMES  [%.0fs]" % (time.time() - t0))

    # NMF programmes (rank 16) on non-negative log expression of D's 4,000 genes, all clusters
    from sklearn.decomposition import NMF
    Xn = X[:, keep].astype(np.float64)
    Xn = Xn / np.maximum(Xn.std(0), 1e-6)
    nmf = NMF(n_components=16, init="nndsvda", max_iter=400, random_state=SEED)
    Wn = nmf.fit_transform(Xn)
    nmf_sub = sub_avg(Wn)
    say("NMF rank 16: reconstruction error %.1f, iterations %d  [%.0fs]" % (nmf.reconstruction_err_, nmf.n_iter_, time.time() - t0))

    # small autoencoder latent (16) on the same standardised genes, all clusters, one seed
    nn = torch.nn
    torch.manual_seed(SEED % 10000)
    g = torch.Generator().manual_seed(SEED % 10000)
    Zt = torch.tensor(Z)
    enc = nn.Sequential(nn.Linear(Z.shape[1], 256), nn.GELU(), nn.Linear(256, 16))
    dec = nn.Sequential(nn.Linear(16, 256), nn.GELU(), nn.Linear(256, Z.shape[1]))
    opt = torch.optim.AdamW(list(enc.parameters()) + list(dec.parameters()), lr=2e-3, weight_decay=1e-4)
    EP = 40
    sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, EP)
    curve = []
    for ep in range(EP):
        perm = torch.randperm(len(Zt), generator=g)
        tot = 0.0
        for i in range(0, len(Zt), 256):
            x = Zt[perm[i:i + 256]]
            m = torch.rand(x.shape, generator=g) < 0.2
            loss = ((dec(enc(x.masked_fill(m, 0.0))) - x) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += float(loss) * len(x)
        sch.step()
        curve.append(tot / len(Zt))
    with torch.no_grad():
        Hae = enc(Zt).numpy().astype(float)
    ae_sub = sub_avg(Hae)
    say("autoencoder 4000-256-16: loss %.3f -> %.3f  [%.0fs]" % (curve[0], curve[-1], time.time() - t0))
    del Z, Zt, Xn, Wn, abc, X

    # gene bases over NAMES (normalised as D: centred over NAMES, one global scale)
    F = {"gD32": z["F_D"], "gD8": norm_rows(z["wb_emb"][:, :8]), "gK32": norm_rows(z["F_K"][:, :32]),
         "gB10": norm_rows(z["F_B"]), "nmf16": norm_rows(np.log1p(nmf_sub[ni])), "ae16": norm_rows(ae_sub[ni])}

    # MERFISH: unit placement, target-side features, subclass location for the modules
    mer = b5.Merfish()
    say("MERFISH read  [%.0fs]" % (time.time() - t0))
    cen, lv = {}, {}
    for src, s in S.items():
        for u in s["units"]:
            c, l = mer.unit_centroid(u, src)
            cen[(src, u)] = c / 1000.0
            lv[l] = lv.get(l, 0) + 1
    cell_sub = np.full(len(mer.xyz), -1)
    for a, u in enumerate(subn):
        ix = mer.idx.get(u)
        if ix is not None:
            cell_sub[ix] = a
    tcen = np.array([mer.target_centroid(t) / 1000.0 for t in targets])
    tcomp = np.zeros((len(targets), nsub))
    for a, t in enumerate(targets):
        m = target_mask(mer, t) & (cell_sub >= 0)
        tcomp[a] = np.bincount(cell_sub[m], minlength=nsub) / m.sum()
    # subclass whole-brain centroid (folded) and division composition, subclasses with >= 50 cells
    divs = [d for d in ("Isocortex", "OLF", "HPF", "CTXsp", "STR", "PAL", "TH", "HY", "MB", "P", "MY", "CB")]
    sub_ok, sub_xyz, sub_div = [], [], []
    for a, u in enumerate(subn):
        ix = mer.idx.get(u)
        if ix is None or len(ix) < 50:
            continue
        sub_ok.append(a)
        sub_xyz.append(mer.xyz[ix].mean(0) / 1000.0)
        dv = mer.div[ix]
        sub_div.append([(dv == d).mean() for d in divs])
    sub_ok = np.array(sub_ok)
    sub_xyz = np.array(sub_xyz)
    sub_div = np.array(sub_div)
    del mer
    say("placements: units by level %s; %d subclasses with >= 50 MERFISH cells  [%.0fs]" % (lv, len(sub_ok), time.time() - t0))

    # target-side encodings
    def enc_norm(M):
        M = M - M.mean(0)
        M = M / np.maximum(M.std(0), 1e-9)
        return M / np.sqrt((M ** 2).sum(1).mean())
    tg_raw = tcomp @ emb_all[:, :32]
    U, s_, _ = np.linalg.svd(tg_raw - tg_raw.mean(0), full_matrices=False)
    tg = enc_norm((U * s_)[:, :4])
    U2, s2, _ = np.linalg.svd(np.sqrt(tcomp) - np.sqrt(tcomp).mean(0), full_matrices=False)
    tco = enc_norm((U2 * s2)[:, :4])
    tc = enc_norm(tcen)

    # modules M1 (genes -> CCF location) and M2 (genes -> division composition), all atlas subclasses
    Xg = emb_all[sub_ok, :32]
    gnames = [subn[a] for a in sub_ok]
    m1 = krr_module(Xg, sub_xyz, gnames, "genes->location")
    m2 = krr_module(Xg, sub_div, gnames, "genes->division composition")
    say("M1 genes->location: OOF R2 %s, RMSE %.2f mm (mean baseline %.2f)  [%.0fs]" % (
        ["%.2f" % x for x in m1["acc"]["r2_per_output_oof"]], m1["acc"]["rmse_oof"], m1["acc"]["rmse_mean_baseline"], time.time() - t0))
    say("M2 genes->division composition: OOF R2 overall %.2f  [%.0fs]" % (m2["acc"]["r2_overall_oof"], time.time() - t0))
    cards = [save_module(m1, "M1_genes_to_location",
                         "Allen subclass whole-brain PCA embedding (32, D's basis) -> MERFISH centroid in CCF (mm, folded)",
                         "wb_emb[:, :32] of an Allen subclass", ["x_ccf_mm", "y_ccf_mm", "z_ccf_folded_mm"],
                         {"training_data": "Allen MERFISH C57BL6J-638850 cells, subclasses with >= 50 cells (unlabelled atlas, all divisions)"}),
             save_module(m2, "M2_genes_to_division_composition",
                         "Allen subclass whole-brain PCA embedding (32) -> share of its MERFISH cells in 12 CCF divisions",
                         "wb_emb[:, :32] of an Allen subclass", divs,
                         {"training_data": "Allen MERFISH C57BL6J-638850 cells, subclasses with >= 50 cells (unlabelled atlas)"})]
    # inferred (out-of-fold) location / composition for NAMES; subclasses without enough cells use the full fit
    pos = {u: a for a, u in enumerate(gnames)}
    missing = [u for u in names if u not in pos]
    say("NAMES not in the module training set (< 50 MERFISH cells; full-fit prediction used): %s" % missing)

    def infer(m, u):
        if u in pos:
            return m["oof"][pos[u]], m["pred_sd"][pos[u]]
        x = (emb_all[subn.index(u), :32] - m["x_mu"]) / m["x_sd"]
        return m["model"].predict(x[None])[0] + m["y_mean"], np.sqrt(m["boot_sd"].mean(0) ** 2 + np.array(m["acc"]["resid_sd_per_output"]) ** 2)
    Lhat = np.array([infer(m1, u)[0] for u in names])
    Lhat_sd = np.array([infer(m1, u)[1] for u in names])
    Lcomp = np.array([infer(m2, u)[0] for u in names])

    keys = sorted(cen)
    np.savez(PREP, names=np.array(names), targets=np.array(targets),
             **{"F_" + k: v for k, v in F.items()},
             Lhat=Lhat, Lhat_sd=Lhat_sd, Lcomp=Lcomp,
             unit_keys=np.array(["%s|%s" % k for k in keys]), unit_cen=np.array([cen[k] for k in keys]),
             tcen=tcen, enc_tc=tc, enc_tg=tg, enc_tcomp=tco, divs=np.array(divs),
             ae_curve=np.array(curve), nmf_err=np.array(nmf.reconstruction_err_))
    json.dump(cards, open(os.path.join(MODDIR, "atlas_module_cards.json"), "w"), indent=1)
    say("wrote %s and %d module cards  [%.0fs]" % (os.path.relpath(PREP, ROOT), len(cards), time.time() - t0))


# ============================================================================ the engine
class Engine:
    """Typed primitives and the shared-coefficient ridge (D's head generalised).

    A pipeline is a list of TERMS; each term maps a (unit, target) pair to features:
      ("K", node, enc): kron(unit-side node features, target-side encoding); unit features are
                        centred within the source, the encoding is centred over the source's targets
                        (encoding "onehot" with node "gD32" is exactly model D's design).
      ("P", kind):      a pairwise unit x target feature (double-centred within the source), one coefficient.
    All terms are fitted jointly by one ridge (penalty lambda on every coefficient) with D's
    sqrt(n) unit weights; terms are pre-scaled to about one unit of Gram trace per observation."""

    UNIT_NODES = ("gD32", "gD8", "gK32", "gB10", "nmf16", "ae16", "rff_gD8", "poly_gD4",
                  "Lm", "rff_Lm", "Lhat", "rff_Lhat", "Lcomp", "Ehat", "gD4xLm", "gD4xLhat")
    PAIRS = ("distLm", "distLhat")

    def __init__(self, tab=None, S=None):
        p = np.load(PREP, allow_pickle=True)
        self.names = list(p["names"])
        self.NI = {u: a for a, u in enumerate(self.names)}
        self.targets = list(p["targets"])
        self.TI = {t: a for a, t in enumerate(self.targets)}
        self.F = {k: p["F_" + k] for k in GENE_NODES}
        self.tcen = p["tcen"]
        self.enc = {"onehot": np.eye(len(self.targets)), "tc": p["enc_tc"], "tg": p["enc_tg"], "tcomp": p["enc_tcomp"]}
        self.cen = {tuple(k.split("|", 1)): v for k, v in zip(p["unit_keys"], p["unit_cen"])}
        self.Lhat_mm = p["Lhat"]
        self.Lhat_sd = p["Lhat_sd"]
        allc = np.array(list(self.cen.values()))
        self.lm_mu, self.lm_sc = allc.mean(0), np.sqrt(((allc - allc.mean(0)) ** 2).sum(1).mean())
        self.F["Lhat"] = norm_rows(self.Lhat_mm)
        self.F["Lcomp"] = norm_rows(np.sqrt(np.clip(p["Lcomp"], 0, None)))
        r1, *_ = rff(self.F["gD8"], 32, SEED + 1)
        self.F["rff_gD8"] = norm_rows(r1)
        self.F["poly_gD4"] = norm_rows(poly2(self.F["gD32"][:, :4] / np.std(self.F["gD32"][:, :4], axis=0)))
        r3, *_ = rff(self.F["Lhat"], 16, SEED + 3)
        self.F["rff_Lhat"] = norm_rows(r3)
        lm = (allc - self.lm_mu) / self.lm_sc
        _, self.rff_lm_bw, self.rff_lm_W, self.rff_lm_b = rff(lm, 16, SEED + 2)
        g4 = self.F["gD32"][:, :4] / np.std(self.F["gD32"][:, :4], axis=0)
        self.g4 = g4
        self.F["gD4xLhat"] = norm_rows(np.einsum("ua,ub->uab", g4, self.F["Lhat"] / self.F["Lhat"].std(0)).reshape(len(g4), -1))
        # global scales of the per-unit (source-dependent) nodes, over all training units
        keys = list(self.cen)
        self.unit_scale = {}
        for node in ("Lm", "rff_Lm", "gD4xLm"):
            M = np.array([self._raw(node, k[0], k[1]) for k in keys])
            M = M - M.mean(0)
            self.unit_scale[node] = np.sqrt((M ** 2).sum(1).mean())
        self.pair_scale = {}
        if S is not None:
            for kind in self.PAIRS:
                v = np.concatenate([self._pair_raw(kind, k, s["units"], s["targets"]).ravel() for k, s in S.items()])
                self.pair_scale[kind] = np.sqrt((v ** 2).mean())
        self.Ecache = {}
        self.ecomp = None
        if tab is not None:
            self._prep_E(tab)

    def scales_for(self, S_sub):
        """the label-free scale constants (unit_scale, pair_scale) recomputed over the sources of
        S_sub only (audit item 11 sensitivity: keep the outer held-out division out of them)."""
        keys = [(k, u) for k, s in S_sub.items() for u in s["units"]]
        us = {}
        for node in ("Lm", "rff_Lm", "gD4xLm"):
            M = np.array([self._raw(node, k, u) for k, u in keys])
            M = M - M.mean(0)
            us[node] = np.sqrt((M ** 2).sum(1).mean())
        ps = {}
        for kind in self.PAIRS:
            v = np.concatenate([self._pair_raw(kind, k, s["units"], s["targets"]).ravel() for k, s in S_sub.items()])
            ps[kind] = np.sqrt((v ** 2).mean())
        return us, ps

    # ---- methylation (L1) composition per (source, subclass), for the cross-fitted G -> E-hat conversion
    def _prep_E(self, tab):
        k = tab[~tab["sub"].isin(b5.HY19)]
        assert (k["Source"] != "HY").all()
        self.L1 = sorted(k["L1"].astype(str).unique())
        LI = {c: a for a, c in enumerate(self.L1)}
        ec = {}
        for (src, sub, ex), g in k.groupby(["Source", "sub", "Experiment"]):
            v = np.bincount([LI[c] for c in g["L1"].astype(str)], minlength=len(self.L1)).astype(float)
            ec.setdefault((src, sub), []).append(v)
        self.ecomp = ec

    def E_features(self, train_srcs, held_subs, train_units, test_units):
        """G -> E-hat (L1 composition) for one fold: compositions from the training sources' neurons
        only (each experiment with >= 3 neurons of the subclass weighted equally), an RBF kernel
        ridge from gD32 fitted on the training subclasses, CROSS-FITTED (5 folds over subclasses)
        for training subclasses, full fit for test subclasses; output: top-8 principal components of
        sqrt(composition) over the training subclasses."""
        key = (tuple(sorted(train_srcs)), tuple(sorted(held_subs)))
        if key in self.Ecache:
            return self.Ecache[key]
        from sklearn.kernel_ridge import KernelRidge
        from sklearn.model_selection import KFold
        comp = {}
        for (src, sub), vs in self.ecomp.items():
            if src not in train_srcs or sub in held_subs or sub not in self.NI:
                continue
            comp.setdefault(sub, []).extend(vs)
        subs = sorted(u for u in comp if u in train_units)
        Y = []
        for u in subs:
            vs = [v / v.sum() for v in comp[u] if v.sum() >= 3]
            if not vs:
                vs = [np.sum(comp[u], axis=0) / np.sum(comp[u])]
            Y.append(np.sqrt(np.mean(vs, axis=0)))
        Y = np.array(Y)
        X = self.F["gD32"][[self.NI[u] for u in subs]]
        hp = {"alpha": 0.1, "gamma": 0.3 / X.shape[1]}
        ym = Y.mean(0)
        oof = np.zeros_like(Y)
        for tr, te in KFold(5, shuffle=True, random_state=SEED).split(X):
            oof[te] = KernelRidge(kernel="rbf", **hp).fit(X[tr], Y[tr] - ym).predict(X[te]) + ym
        full = KernelRidge(kernel="rbf", **hp).fit(X, Y - ym)
        U, s, Vt = np.linalg.svd(Y - ym, full_matrices=False)
        V8 = Vt[:8].T
        out = {u: (oof[a] - ym) @ V8 for a, u in enumerate(subs)}
        others = sorted(set(test_units) | (set(train_units) - set(subs)))
        if others:
            P = full.predict(self.F["gD32"][[self.NI[u] for u in others]]) + ym
            for a, u in enumerate(others):
                out[u] = (P[a] - ym) @ V8
        M = np.array([out[u] for u in subs])
        sc = np.sqrt(((M - M.mean(0)) ** 2).sum(1).mean())
        out = {u: v / sc for u, v in out.items()}
        self.Ecache[key] = out
        return out

    # ---- raw unit features
    def _raw(self, node, src, u, ctx=None):
        if node in self.F:
            return self.F[node][self.NI[u]]
        if node == "Lm":
            return (self.cen[(src, u)] - self.lm_mu) / self.lm_sc
        if node == "rff_Lm":
            lm = (self.cen[(src, u)] - self.lm_mu) / self.lm_sc
            return np.sqrt(2.0 / 16) * np.cos(lm @ self.rff_lm_W + self.rff_lm_b)
        if node == "gD4xLm":
            lm = (self.cen[(src, u)] - self.lm_mu) / self.lm_sc
            return np.outer(self.g4[self.NI[u]], lm).ravel()
        if node == "Ehat":
            return ctx["E"][u]
        raise KeyError(node)

    def unit_X(self, node, src, units, ctx=None):
        X = np.array([self._raw(node, src, u, ctx) for u in units])
        if node in self.unit_scale:
            X = X / self.unit_scale[node]
        return X - X.mean(0)

    def _pair_raw(self, kind, src, units, targets):
        if kind == "distLm":
            P = np.array([self.cen[(src, u)] for u in units])
        else:
            P = self.Lhat_mm[[self.NI[u] for u in units]]
        T = self.tcen[[self.TI[t] for t in targets]]
        d = np.sqrt(((P[:, None] - T[None]) ** 2).sum(-1))
        return b5.dc(-np.log(np.maximum(d, A_MIN_MM)))

    def pair_D(self, kind, src, units, targets):
        return self._pair_raw(kind, src, units, targets) / self.pair_scale[kind]

    def enc_Q(self, enc, targets):
        Q = self.enc[enc][[self.TI[t] for t in targets]]
        return Q - Q.mean(0)


def folds(S):
    """d105 folds (leave one division out; held-out subclasses removed from the training sources),
    with the full held-out source records (counts included) as the test side."""
    for d, train, test in b5.folds(S):
        yield d, train, {k: S[k] for k in test}


def src_obs(s):
    """observed double-centred profile and D's weights (w^2 = n / mean n) for one source."""
    Y = b5.dc(b5.src_profile(s))
    w2 = s["n"] / s["n"].mean()
    return Y, w2


class Fold:
    """Cached per-source pieces of one train/test split for the engine."""

    def __init__(self, eng, train, test, test_obs=None, use_E=True):
        self.eng = eng
        self.train = train
        self.test = test
        self.ytr = {k: src_obs(s) for k, s in train.items()}
        self.obs = test_obs
        held = {u for s in test.values() for u in s["units"]}
        tu = {u for s in train.values() for u in s["units"]}
        teu = {u for s in test.values() for u in s["units"]}
        self.ctx = {}
        if use_E and eng.ecomp is not None:
            self.ctx["E"] = eng.E_features(set(train), held, tu, teu)
        self.cache = {}

    def piece(self, term, k, side):
        key = (term, k, side)
        if key not in self.cache:
            s = (self.train if side == "tr" else self.test)[k]
            if term[0] == "K":
                self.cache[key] = (self.eng.unit_X(term[1], k, s["units"], self.ctx), self.eng.enc_Q(term[2], s["targets"]))
            else:
                self.cache[key] = self.eng.pair_D(term[1], k, s["units"], s["targets"])
        return self.cache[key]

    def dim(self, term):
        if term[0] == "P":
            return 1
        k = next(iter(self.train))
        X, Q = self.piece(term, k, "tr")
        return X.shape[1] * Q.shape[1]

    def block(self, ta, tb):
        key = ("G", ta, tb)
        if key in self.cache:
            return self.cache[key]
        if ta[0] == "P" and tb[0] == "K":
            out = self.block(tb, ta).T
            self.cache[key] = out
            return out
        acc = None
        if ta[0] == "K" and tb[0] == "K":
            As, Cs = [], []
            for k in self.train:
                Xa, Qa = self.piece(ta, k, "tr")
                Xb, Qb = self.piece(tb, k, "tr")
                w2 = self.ytr[k][1]
                As.append((Xa.T * w2) @ Xb)
                Cs.append(Qa.T @ Qb)
            pa, pb = As[0].shape
            qa, qb = Cs[0].shape
            M = np.stack([a.ravel() for a in As]).T @ np.stack([c.ravel() for c in Cs])
            acc = M.reshape(pa, pb, qa, qb).transpose(0, 2, 1, 3).reshape(pa * qa, pb * qb)
        elif ta[0] == "K":
            acc = 0
            for k in self.train:
                Xa, Qa = self.piece(ta, k, "tr")
                D = self.piece(tb, k, "tr")
                w2 = self.ytr[k][1]
                acc = acc + ((Xa.T * w2) @ D @ Qa).ravel()
            acc = acc[:, None]
        else:
            acc = 0.0
            for k in self.train:
                w2 = self.ytr[k][1]
                acc = acc + float((w2[:, None] * self.piece(ta, k, "tr") * self.piece(tb, k, "tr")).sum())
            acc = np.array([[acc]])
        self.cache[key] = acc
        return acc

    def rhs(self, t):
        key = ("r", t)
        if key not in self.cache:
            acc = 0
            for k in self.train:
                Y, w2 = self.ytr[k]
                if t[0] == "K":
                    X, Q = self.piece(t, k, "tr")
                    acc = acc + ((X.T * w2) @ Y @ Q).ravel()
                else:
                    acc = acc + np.array([float((w2[:, None] * self.piece(t, k, "tr") * Y).sum())])
            self.cache[key] = acc
        return self.cache[key]

    def solve(self, terms, lams):
        from scipy.linalg import cho_factor, cho_solve
        G = np.block([[self.block(a, b) for b in terms] for a in terms])
        r = np.concatenate([self.rhs(t) for t in terms])
        out = []
        for lam in lams:
            c = cho_factor(G + lam * np.eye(len(G)), lower=False, check_finite=False)
            out.append(cho_solve(c, r, check_finite=False))
        return out

    def predict(self, terms, theta, side="te"):
        S = self.test if side == "te" else self.train
        preds = {}
        for k in S:
            o = 0
            P = 0
            for t in terms:
                if t[0] == "K":
                    X, Q = self.piece(t, k, side)
                    d = X.shape[1] * Q.shape[1]
                    P = P + X @ theta[o:o + d].reshape(X.shape[1], Q.shape[1]) @ Q.T
                    o += d
                else:
                    P = P + theta[o] * self.piece(t, k, side)
                    o += 1
            preds[k] = P
        return preds


def r_division(S_test, preds):
    """devkit r_unit for one held-out division (sources pooled, sqrt(n / mean n) unit weights)."""
    srcs = sorted(S_test)
    nn = np.concatenate([S_test[k]["n"] for k in srcs]).astype(float)
    P, O, W = [], [], []
    for k in srcs:
        o = b5.dc(b5.src_profile(S_test[k]))
        P.append(b5.unitrows(b5.dc(np.asarray(preds[k], float))).ravel())
        O.append(b5.unitrows(o).ravel())
        W.append((np.sqrt(S_test[k]["n"] / nn.mean())[:, None] * np.ones((1, o.shape[1]))).ravel())
    return b5.wcorr(np.concatenate(P), np.concatenate(O), np.concatenate(W))


def unit_cos(S_test, preds):
    """per unit: cosine of the predicted and measured double-centred rows -> {(src, unit): cos}."""
    out = {}
    for k in S_test:
        o = b5.unitrows(b5.dc(b5.src_profile(S_test[k])))
        p = b5.unitrows(b5.dc(np.asarray(preds[k], float)))
        for a, u in enumerate(S_test[k]["units"]):
            out[(k, u)] = float(o[a] @ p[a])
    return out


def apply_bundle(path, units, targets):
    """apply a saved discovered pipeline (S1_discovered_pipeline.npz) to ONE source: its subclasses
    `units` (Allen subclass names present in the bundle's `names`) and its `targets` (names in the
    bundle's `targets`). Reproduces Fold.predict for gene-node x target-encoding terms: unit features
    are the bundle's rows for these units centred over the source's units, the target encoding is
    centred over the source's targets, prediction = sum_terms X theta_term Q^T (a double-centred
    log-share profile, INFERRED)."""
    z = np.load(path, allow_pickle=True)
    terms = json.loads(str(z["terms"]))
    names = list(z["names"])
    tg = list(z["targets"])
    th = z["theta"]
    pred, o = 0, 0
    for t in terms:
        if t[0] != "K":
            raise NotImplementedError("bundle holds gene-node x encoding terms only")
        X = z["F_" + t[1]][[names.index(u) for u in units]]
        X = X - X.mean(0)
        Q = z["enc_" + t[2]][[tg.index(x) for x in targets]]
        Q = Q - Q.mean(0)
        d = X.shape[1] * Q.shape[1]
        pred = pred + X @ th[o:o + d].reshape(X.shape[1], Q.shape[1]) @ Q.T
        o += d
    assert o == len(th)
    return pred


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "prep":
        prep()
    else:
        sys.exit("usage: python3 -I analysis/d108_prims.py prep")
