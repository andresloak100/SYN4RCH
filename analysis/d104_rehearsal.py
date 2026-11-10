"""D-104 rehearsal and the frozen outside rule: every Epi-retro-seq source EXCEPT the hypothalamus.

    python3 -I analysis/d104_rehearsal.py transfer [nperm]   # ridge 'learn elsewhere' rule, leave one division out
    python3 -I analysis/d104_rehearsal.py similarity [nperm] # the similarity test (design A) on each non-HY division
    python3 -I analysis/d104_rehearsal.py freeze             # fit the rule on all non-HY sources -> d104_outside_rule.npz
    python3 -I analysis/d104_rehearsal.py report             # all three, written to reports/d104_rehearsal.md

No hypothalamic-source neuron is read: Supplementary Table 2 is filtered to Source != "HY" as
soon as it is loaded, and the 19 hypothalamic test subclasses (HY19) are removed from every
other source. Nothing here relates genes to the projections the D-104 test is scored on. Gene
means include the 19 test subclasses (genes only), so the rule's gene features are defined on
them. d104_tables.npz is opened only for its gene symbols.

`freeze` writes the rule that d104_epiretro.outside_rule() applies (s11). This module imports
d104_epiretro for its loaders and statistics; d104_epiretro does not import this module.
"""
import os
import re
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, os.path.join(ROOT, "src"))
import d104_epiretro as E  # noqa: E402

GENES_CACHE = os.path.join(ROOT, "data/interim/d104_rehearsal_genes.npz")
REPORT = os.path.join(ROOT, "reports/d104_rehearsal.md")

HY19 = ["PVH-SO-PVa Otp Glut", "AVPV-MEPO-SFO Tbr1 Glut", "DMH Hmx2 Gaba", "LHA-AHN-PVH Otp Trh Glut",
        "MPN-MPO-PVpo Hmx2 Glut", "PVpo-VMPO-MPN Hmx2 Gaba", "AHN Onecut3 Gaba", "AHN-RCH-LHA Otp Fezf1 Glut",
        "AHN-SBPV-PVHd Pdrm12 Gaba", "BST-MPN Six3 Nrgn Gaba", "PH-ant-LHA Otp Bsx Glut",
        "PVHd-SBPV Six3 Prox1 Gaba", "VMH Nr5a1 Glut", "DMH-LHA Gsx1 Gaba", "DMH-LHA Vgll2 Glut",
        "LHA Barhl2 Glut", "PVHd-DMH Lhx6 Gaba", "STN-PSTN Pitx2 Glut", "ZI Pax6 Gaba"]
HY_TARGETS = ("PFC", "MOB", "STR", "TH", "SC", "P", "MY", "AMY", "PAL", "VTA")
DIVISION = {"MOp": "isocortex", "ACA": "isocortex", "SSp": "isocortex", "AUD": "isocortex", "PTLp": "isocortex",
            "AI": "isocortex", "VIS": "isocortex", "RSP": "isocortex",
            "MOB": "olfactory", "AON": "olfactory", "PIRa": "olfactory", "PIRp": "olfactory",
            "ENT": "hippocampal", "CAa": "hippocampal", "CAp": "hippocampal", "DGa": "hippocampal", "DGp": "hippocampal",
            "AMY": "amygdala-pallidum", "PAL": "amygdala-pallidum", "STR": "amygdala-pallidum",
            "THl": "thalamus", "THm": "thalamus", "THp": "thalamus",
            "VTA": "midbrain", "MRN": "midbrain", "SC": "midbrain", "PAG": "midbrain", "IC": "midbrain",
            "P": "hindbrain", "MY": "hindbrain", "CBX": "cerebellum"}
SUBCORTICAL = ("amygdala-pallidum", "thalamus", "midbrain", "hindbrain")
MIN_NEURONS = 20
PSEUDO = 0.5
GRID = [(k, lam) for k in (10, 20, 40) for lam in (100.0, 1000.0, 10000.0)]
FROZEN = {"k": 10, "lam": 1000.0, "wpow": 0.5}


def strip(s):
    return re.sub(r"^\d+\s+", "", str(s))


# --- data: non-hypothalamic sources only ----------------------------------------------------------

_C = {}


def table():
    """Supplementary Table 2, kept neurons from non-hypothalamic sources with a subclass."""
    if "table" not in _C:
        import pandas as pd
        df = pd.read_excel(E.path("st2"), sheet_name="TableS2", header=1)
        k = df[(df["PassTargetFilter"] == 1) & (df["Source"] != "HY") & df["Subclass"].notna()].copy()
        del df
        assert (k["Source"] == "HY").sum() == 0
        k["sub"] = k["Subclass"].map(strip)
        _C["table"] = k
    return _C["table"]


def sources(min_neurons=MIN_NEURONS, drop=HY19):
    """Per source: units (subclasses with >= min_neurons labelled neurons there, test subclasses
    removed), experiments, their targets and sexes, counts[unit, experiment]."""
    key = ("sources", min_neurons, tuple(drop))
    if key in _C:
        return _C[key]
    k = table()
    k = k[~k["sub"].isin(drop)]
    out = {}
    for s, g in k.groupby("Source"):
        vc = g["sub"].value_counts()
        units = sorted(u for u, n in vc.items() if n >= min_neurons)
        exps = sorted(g["Experiment"].unique())
        et = np.array([g.loc[g["Experiment"] == e, "Target"].iloc[0] for e in exps])
        ex = np.array([g.loc[g["Experiment"] == e, "Sex"].iloc[0] for e in exps])
        ct = np.array([[int(((g["sub"] == u) & (g["Experiment"] == e)).sum()) for e in exps] for u in units])
        targets = sorted(set(et))
        if len(units) >= 2 and len(targets) >= 2:
            out[s] = {"units": units, "exps": exps, "exp_target": et, "exp_sex": ex, "counts": ct,
                      "targets": targets, "n": ct.sum(axis=1), "division": DIVISION[s]}
    _C[key] = out
    return out


def names():
    S = sources()
    return sorted({u for s in S.values() for u in s["units"]} | set(HY19))


def genes(nm):
    """log2(CPM+1) subclass means (Allen 10x cluster means, cell-weighted), cached."""
    if os.path.exists(GENES_CACHE):
        z = np.load(GENES_CACHE, allow_pickle=True)
        if list(z["names"]) == list(nm):
            return z["means"], z["symbols"]
    means, sym = E.subclass_means(list(nm))
    np.savez(GENES_CACHE, names=np.array(list(nm)), means=means, symbols=sym)
    return means, sym


def dc(A):
    return A - A.mean(axis=1, keepdims=True) - A.mean(axis=0, keepdims=True) + A.mean()


def log_shares(S):
    """Shares among the source's units, PSEUDO added, averaged over a target's experiments, logged."""
    c = S["counts"]
    q = (c + PSEUDO) / (c.sum(axis=0) + PSEUDO * c.shape[0])
    return np.log(np.column_stack([q[:, S["exp_target"] == t].mean(axis=1) for t in S["targets"]]))


def feature_basis(means, sym, n_genes=2000, k=None):
    """The most variable genes across the named subclasses (no mt-/Rpl/Rps), z-scored; with k,
    their top-k principal component scores. Returns scores and what reproduces them."""
    var = means.var(axis=0)
    var[np.array([str(g).startswith(("mt-", "Rpl", "Rps")) for g in sym])] = -1
    keep = np.argsort(var)[::-1][:n_genes]
    F = means[:, keep]
    mu = F.mean(axis=0)
    sd = F.std(axis=0)
    sd[sd == 0] = 1
    Z = (F - mu) / sd
    V = None
    if k:
        U, s, Vt = np.linalg.svd(Z, full_matrices=False)
        Z = U[:, :k] * s[:k]
        V = Vt[:k].T
    return Z, {"gene_idx": keep, "mean": mu, "sd": sd, "V": V}


def features(nm, means, sym, n_genes=2000, k=None):
    Z, _ = feature_basis(means, sym, n_genes, k)
    return dict(zip(nm, Z))


# --- the leave-one-division-out protocol ------------------------------------------------------------

def unitrows(A):
    n = np.linalg.norm(A, axis=1, keepdims=True)
    return A / np.where(n > 0, n, 1.0)


def subset(s, keep):
    return dict(s, units=[s["units"][a] for a in keep], counts=s["counts"][keep], n=s["n"][keep])


def folds(S):
    """Leave one division out. Held out: every source of the division. Training: every other
    source, with any subclass that appears among the held-out units removed."""
    for d in sorted({s["division"] for s in S.values()}):
        test_src = [k for k in sorted(S) if S[k]["division"] == d]
        held = {u for k in test_src for u in S[k]["units"]}
        train = {}
        for k in sorted(S):
            if S[k]["division"] == d:
                continue
            keep = [a for a, u in enumerate(S[k]["units"]) if u not in held]
            if len(keep) >= 2:
                train[k] = subset(S[k], keep)
        test = {k: {"units": list(S[k]["units"]), "targets": list(S[k]["targets"]), "division": d} for k in test_src}
        yield d, train, test


def _wcorr(a, b, w):
    a = a * w
    b = b * w
    den = np.sqrt((a @ a) * (b @ b))
    return float(a @ b / den) if den > 0 else 0.0


def lodo(fit_predict, nperm=2000, seed=0, S=None):
    """Per held-out division, pooled over its sources: r_unit (weighted correlation of unit-row
    normalised double-centred prediction and measured profile, weights sqrt(neurons)), r_entry
    (the same without row normalisation), p (predicted rows permuted among units within source).
    Selection metric: the mean of r_unit over the divisions."""
    S = S or sources()
    rng = np.random.default_rng(seed)
    out = {}
    for d, train, test in folds(S):
        pred = fit_predict(train, test)
        P, O, Pe, Oe, W = [], [], [], [], []
        nn = np.concatenate([S[k]["n"] for k in test]).astype(float)
        for k in test:
            p = np.asarray(pred[k], float)
            assert p.shape == (len(S[k]["units"]), len(S[k]["targets"])), (k, p.shape)
            o = dc(log_shares(S[k]))
            w = np.sqrt(S[k]["n"] / nn.mean())[:, None] * np.ones((1, o.shape[1]))
            P.append(unitrows(dc(p))); O.append(unitrows(o)); Pe.append(dc(p)); Oe.append(o); W.append(w)

        def cat(L):
            return np.concatenate([x.ravel() for x in L])
        r = _wcorr(cat(P), cat(O), cat(W))
        re_ = _wcorr(cat(Pe), cat(Oe), cat(W))
        null = np.empty(nperm)
        for b in range(nperm):
            null[b] = _wcorr(cat([x[rng.permutation(len(x))] for x in P]), cat(O), cat(W))
        out[d] = {"r_unit": r, "r_entry": re_, "p": float((1 + (null >= r).sum()) / (1 + nperm)),
                  "units": int(len(nn)), "sources": list(test)}
    divs = sorted(out)
    rs = np.array([out[d]["r_unit"] for d in divs])
    sub = np.array([out[d]["r_unit"] for d in SUBCORTICAL])
    out["_summary"] = {"mean_r_unit": float(rs.mean()), "min_r_unit": float(rs.min()),
                       "positive_divisions": int((rs > 0).sum()), "n_divisions": len(divs),
                       "subcortical_mean_r_unit": float(sub.mean())}
    return out


# --- the ridge rule ---------------------------------------------------------------------------------

_G = {}


def ridge_inputs():
    if not _G:
        S = sources()
        _G["S"] = S
        _G["TARGETS"] = sorted({t for s in S.values() for t in s["targets"]})
        _G["TI"] = {t: a for a, t in enumerate(_G["TARGETS"])}
        _G["NAMES"] = names()
        _G["MEANS"], _G["SYM"] = genes(_G["NAMES"])
    return _G


def make(k=20, lam=1000.0, wpow=0.5, keep_b=None):
    """Ridge from gene PCs to the double-centred log-share profile, one coefficient vector per
    target shared by every source. keep_b: a dict that receives the fitted coefficients."""
    g = ridge_inputs()
    TARGETS, TI = g["TARGETS"], g["TI"]
    F = features(g["NAMES"], g["MEANS"], g["SYM"], k=k)

    def tmap(targets):
        M = np.zeros((len(targets), len(TARGETS)))
        for a, t in enumerate(targets):
            M[a, TI[t]] = 1.0
        return M - M.mean(axis=0, keepdims=True)

    def fit_predict(train, test):
        A, y = [], []
        for s in train.values():
            Y = dc(log_shares(s))
            G = np.array([F[u] for u in s["units"]]); G = G - G.mean(axis=0)
            w = np.repeat((s["n"] / s["n"].mean()) ** wpow, Y.shape[1])
            A.append(np.einsum("ij,tk->itjk", G, tmap(s["targets"])).reshape(Y.size, -1) * w[:, None])
            y.append(Y.ravel() * w)
        A = np.vstack(A); y = np.concatenate(y)
        b = np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ y)
        if keep_b is not None:
            keep_b["b"] = b
        out = {}
        for k_, s in test.items():
            G = np.array([F[u] for u in s["units"]]); G = G - G.mean(axis=0)
            At = np.einsum("ij,tk->itjk", G, tmap(s["targets"])).reshape(len(s["units"]) * len(s["targets"]), -1)
            out[k_] = (At @ b).reshape(len(s["units"]), len(s["targets"]))
        return out
    return fit_predict


def transfer(nperm=2000, verbose=True):
    rows = []
    for k, lam in GRID:
        t0 = time.time()
        o = lodo(make(k, lam), nperm=nperm)
        rows.append((k, lam, o))
        if verbose:
            divs = sorted(d for d in o if d[0] != "_")
            sm = o["_summary"]
            print("PCs %2d lam %6.0f: " % (k, lam)
                  + "  ".join("%s %+.2f(p%.3f)" % (d[:5], o[d]["r_unit"], o[d]["p"]) for d in divs)
                  + " | mean %+.3f, subcortical %+.3f, positive %d/%d  [%.0fs]"
                  % (sm["mean_r_unit"], sm["subcortical_mean_r_unit"], sm["positive_divisions"], sm["n_divisions"],
                     time.time() - t0), flush=True)
    return rows


# --- the similarity rehearsal (design A, location partialled out) -----------------------------------

DIV_CCF = {"isocortex": ["Isocortex"], "olfactory": ["OLF"], "hippocampal": ["HPF"],
           "amygdala-pallidum": ["CTXsp", "STR", "PAL"], "thalamus": ["TH"], "midbrain": ["MB"],
           "hindbrain": ["P", "MY"]}
# a source's position parcel is a reading of its name; sub-dissections share one parcel
SRC_STRUCT = {"MOp": ["MOp"], "ACA": ["ACAd", "ACAv"], "SSp": "SSp", "AUD": "AUD", "PTLp": ["VISa", "VISrl"],
              "AI": ["AId", "AIv", "AIp"], "VIS": ["VISp", "VISl", "VISal", "VISam", "VISpm", "VISpl", "VISli", "VISpor"],
              "RSP": ["RSPd", "RSPv", "RSPagl"], "MOB": ["MOB"], "AON": ["AON"], "PIRa": ["PIR"], "PIRp": ["PIR"],
              "ENT": ["ENTl", "ENTm"], "CAa": ["CA1", "CA2", "CA3"], "CAp": ["CA1", "CA2", "CA3"], "DGa": ["DG"],
              "DGp": ["DG"], "AMY": ["BLA", "BMA", "LA", "PA", "CEA", "MEA", "AAA", "IA", "BA", "COAa", "COAp"],
              "SC": ["SCm", "SCs"], "PAG": ["PAG"], "IC": ["IC"], "MRN": ["MRN"],
              "VTA": ["VTA", "SNc", "SNr", "IF", "RL", "CLI", "IPN", "RR", "PN", "RN"]}
SRC_DIV = {"PAL": ["PAL"], "STR": ["STR"], "THl": ["TH"], "THm": ["TH"], "THp": ["TH"], "P": ["P"], "MY": ["MY"]}

_M = {}


def merfish():
    if not _M:
        import pandas as pd
        m = E.merfish_cells()
        _M["IDX"] = m.groupby("sub").indices
        _M["XYZ"] = np.column_stack([m["x_ccf"], m["y_ccf"], np.abs(m["z_ccf"] - E.MIDLINE_Z_MM)]) * 1000.0
        _M["DIVCOL"] = m["parcellation_division"].values.astype(str)
        _M["STRCOL"] = m["parcellation_structure"].values.astype(str)
        _M["ALLSTR"] = sorted(set(_M["STRCOL"]))
        g = ridge_inputs()
        _M["ROW"] = {n: a for a, n in enumerate(g["NAMES"])}
        _M["PANEL"] = np.array(sorted(set(pd.read_csv(E.path("genes_merfish"))["gene_symbol"]) & set(g["SYM"])))
    return _M


def structs(src):
    v = SRC_STRUCT[src]
    return [x for x in merfish()["ALLSTR"] if x.startswith(v)] if isinstance(v, str) else v


def unit_cells(sub, src, rng, limit=1000):
    M = merfish()
    ix = M["IDX"].get(sub, np.array([], int))
    sel = np.isin(M["STRCOL"][ix], structs(src)) if src in SRC_STRUCT else np.isin(M["DIVCOL"][ix], SRC_DIV[src])
    level = "parcel"
    if sel.sum() < 20:
        sel, level = np.isin(M["DIVCOL"][ix], DIV_CCF[DIVISION[src]]), "division"
    if sel.sum() < 20:
        sel, level = np.ones(len(ix), bool), "brain"
    xyz = M["XYZ"][ix[sel]]
    if len(xyz) > limit:
        xyz = xyz[rng.choice(len(xyz), limit, replace=False)]
    return xyz, level


def div_share(sub, division):
    M = merfish()
    ix = M["IDX"].get(sub, np.array([], int))
    return float(np.isin(M["DIVCOL"][ix], DIV_CCF[division]).mean()) if len(ix) else 0.0


def overlap(cells):
    n = len(cells)
    Kr = np.zeros((len(E.KERNEL_UM), n, n))
    for a in range(n):
        for b in range(a, n):
            d2 = ((cells[a][:, None, :] - cells[b][None, :, :]) ** 2).sum(axis=2)
            for k, sc in enumerate(E.KERNEL_UM):
                Kr[k, a, b] = Kr[k, b, a] = np.exp(-d2 / (4 * sc ** 2)).mean()
    for k in range(len(E.KERNEL_UM)):
        dg = np.sqrt(np.diag(Kr[k]))
        Kr[k] = Kr[k] / np.outer(dg, dg)
    return Kr


def build(division, min_share=0.5, both_sex=True, seed=5, grad=True):
    """The division's design: units (>= 20 neurons in a source, >= min_share of MERFISH cells in
    the division), pairs within source, 14 covariates (8 + six coordinate products), block =
    source, weights 1/(1/n_i + 1/n_j), measured profile similarity."""
    S = sources()
    TAB = table()
    rng = np.random.default_rng(seed)
    units, I, J, block, Drows, w, y, notes = [], [], [], [], [], [], [], []
    for src in [k for k in sorted(S) if S[k]["division"] == division]:
        s = S[src]
        targets = [t for t in s["targets"] if len(set(s["exp_sex"][s["exp_target"] == t])) >= 2] if both_sex else list(s["targets"])
        ecols = np.isin(s["exp_target"], targets)
        keep = [a for a, u in enumerate(s["units"]) if div_share(u, division) >= min_share
                and s["counts"][a, ecols].sum() >= 10]
        if len(keep) < 3 or len(targets) < 3:
            notes.append("%s skipped (%d units, %d both-sex targets)" % (src, len(keep), len(targets)))
            continue
        counts = s["counts"][keep][:, ecols]
        et = s["exp_target"][ecols]
        nm = [s["units"][a] for a in keep]
        n = counts.sum(axis=1).astype(float)
        P, _ = E.profiles(counts, et, s["exp_sex"][ecols], {"targets": tuple(targets)})
        Sp = P @ P.T
        got = [unit_cells(u, src, rng) for u in nm]
        cells = [g[0] for g in got]
        cen = np.array([c.mean(axis=0) for c in cells])
        kemb = overlap(cells)
        g = TAB[TAB["Source"] == src]
        l4 = sorted(g["L4"].astype(str).unique())
        H = np.array([[float(((g["sub"] == u) & (g["L4"].astype(str) == c)).sum()) for c in l4] for u in nm])
        H /= H.sum(axis=1, keepdims=True)
        share = np.minimum(H[:, None, :], H[None, :, :]).sum(axis=2)
        off = len(units)
        cs = cen - cen.mean(axis=0)
        sd = cs.std(axis=0)
        cs = cs / np.where(sd > 0, sd, 1.0)
        for i in range(len(nm)):
            for j in range(i + 1, len(nm)):
                d = float(np.sqrt(((cen[i] - cen[j]) ** 2).sum()))
                Drows.append([-d] + [np.exp(-d ** 2 / (2 * sc ** 2)) for sc in E.KERNEL_UM]
                             + [kemb[k][i, j] for k in range(len(E.KERNEL_UM))] + [share[i, j]]
                             + ([cs[i, a] * cs[j, b] + cs[i, b] * cs[j, a] for a in range(3) for b in range(a, 3)] if grad else []))
                I.append(off + i); J.append(off + j); block.append(src)
                w.append(1.0 / (1.0 / n[i] + 1.0 / n[j])); y.append(Sp[i, j])
        units += [(src, u) for u in nm]
        notes.append("%s %d units, %d targets (%s), position from %s" % (
            src, len(nm), len(targets), ",".join(targets),
            "/".join("%s:%d" % (lv, sum(1 for x in got if x[1] == lv)) for lv in ("parcel", "division", "brain") if any(x[1] == lv for x in got))))
    if not units:
        return None
    w = np.array(w)
    return {"units": units, "I": np.array(I), "J": np.array(J), "block": np.array(block), "D": np.array(Drows),
            "w": w / w.mean(), "y": np.array(y), "notes": notes}


def gene_sim(des, kernel="raw"):
    g, M = ridge_inputs(), merfish()
    subs = sorted({u for _, u in des["units"]})
    rows = np.array([M["ROW"][u] for u in subs])
    t = {"means": g["MEANS"][rows], "symbols": g["SYM"], "panel": M["PANEL"]}
    if kernel == "raw":
        F = E.fingerprints(t, np.arange(len(subs)))
    elif kernel == "panel":
        F = E.fingerprints(t, np.arange(len(subs)), genes="panel")
    else:
        raise ValueError(kernel)
    Sg = F @ F.T
    pos = np.array([subs.index(u) for _, u in des["units"]])
    return Sg[np.ix_(pos, pos)]


def perms(des, n, rng):
    N = len(des["units"])
    src = np.array([s for s, _ in des["units"]])
    P = np.tile(np.arange(N), (n + 1, 1))
    for b in range(1, n + 1):
        for s in sorted(set(src)):          # sorted: set order follows the per-run string hash
            ix = np.where(src == s)[0]
            P[b, ix] = ix[rng.permutation(len(ix))]
    return P


def sim_test(des, kernel="raw", nperm=2000, seed=11):
    Sg = gene_sim(des, kernel)
    P = perms(des, nperm, np.random.default_rng(seed))
    out = {}
    for name, D in (("raw", None), ("adj", "use")):
        r = E.Partial(Sg, des, P, D=D).corr(des["y"])
        nul = r[1:][np.isfinite(r[1:])]
        out[name] = (float(r[0]), float((r[0] - nul.mean()) / nul.std()), float((1 + (nul >= r[0]).sum()) / (1 + len(nul))))
    out["r2_proj"] = E.location_r2(des["y"], des)
    out["r2_gene"] = E.location_r2(Sg[des["I"], des["J"]], des)
    return out


def similarity(nperm=2000, verbose=True):
    designs = {d: build(d) for d in sorted(DIV_CCF)}
    rows = []
    for d, des in designs.items():
        if verbose:
            print(d, "|", "; ".join(des["notes"]) if des else "nothing to compare", flush=True)
    if verbose:
        print("  %-18s %5s %5s %4s  %-22s %-26s %s" % ("division", "units", "pairs", "cov", "no position: r z p",
                                                     "position partialled: r z p", "position explains (proj / genes)"))
    for d, des in designs.items():
        if des is None:
            continue
        o = sim_test(des, nperm=nperm)
        rows.append((d, des, o))
        if verbose:
            print("  %-18s %5d %5d %4d  %+.2f %+5.2f %.3f      %+.2f %+5.2f %.3f      %4.0f%% / %2.0f%%" % (
                d, len(des["units"]), len(des["I"]), des["D"].shape[1], *o["raw"], *o["adj"],
                100 * o["r2_proj"], 100 * o["r2_gene"]), flush=True)
    Z = [(d, o["raw"][1], o["adj"][1]) for d, _, o in rows]
    sub = [z for z in Z if z[0] in SUBCORTICAL]
    summ = (np.mean([z[1] for z in Z]), np.mean([z[2] for z in Z]), np.mean([z[1] for z in sub]), np.mean([z[2] for z in sub]))
    if verbose:
        print("  mean z: all %+.2f raw, %+.2f adjusted | four subcortical %+.2f raw, %+.2f adjusted" % summ)
    return designs, rows, summ


# --- freeze the outside rule ------------------------------------------------------------------------

def freeze(verbose=True):
    g = ridge_inputs()
    S, TARGETS = g["S"], g["TARGETS"]
    k, lam, wpow = FROZEN["k"], FROZEN["lam"], FROZEN["wpow"]
    missing = [t for t in HY_TARGETS if t not in TARGETS]
    assert not missing, "hypothalamic targets absent from the outside rule: %s" % missing
    for s in S.values():
        assert not set(s["units"]) & set(HY19)
    tsym = np.load(E.TABLES, allow_pickle=True)["symbols"]      # gene symbols only
    assert list(tsym) == list(g["SYM"]), "subclass_means symbols differ from d104_tables.npz"

    Z, basis = feature_basis(g["MEANS"], g["SYM"], k=k)
    kb = {}
    fp = make(k, lam, wpow, keep_b=kb)
    check = sorted(S)[:: max(1, len(S) // 4)][:4]
    via_fit = fp(S, {c: {"units": S[c]["units"], "targets": S[c]["targets"]} for c in check})
    B = kb["b"].reshape(k, len(TARGETS))
    R = {"gene_idx": basis["gene_idx"].astype(np.int64), "symbols": g["SYM"][basis["gene_idx"]],
         "mean": basis["mean"], "sd": basis["sd"], "V": basis["V"], "B": B,
         "targets": np.array(TARGETS), "k": k, "lam": lam}

    # self-checks: the npz reproduces the PC scores and, through outside_rule's formula, the fit
    assert np.allclose(((g["MEANS"][:, R["gene_idx"]] - R["mean"]) / R["sd"]) @ R["V"], Z, atol=1e-8)
    row = {n: a for a, n in enumerate(g["NAMES"])}
    worst = 0.0
    for c in check:
        idx = [row[u] for u in S[c]["units"]]
        cols = [list(R["targets"]).index(t) for t in S[c]["targets"]]
        X = ((g["MEANS"][idx][:, R["gene_idx"]] - R["mean"]) / R["sd"]) @ R["V"]
        G = X - X.mean(axis=0)
        Bc = R["B"][:, cols]
        pred = G @ (Bc - Bc.mean(axis=1, keepdims=True))
        diff = float(np.abs(dc(pred) - dc(via_fit[c])).max())
        worst = max(worst, diff)
        assert diff < 1e-8, (c, diff)

    np.savez(E.OUTSIDE_RULE, **R)
    summ = {"k": k, "lam": lam, "wpow": wpow, "genes": len(R["gene_idx"]), "targets": TARGETS,
            "sources": len(S), "units": int(sum(len(s["units"]) for s in S.values())),
            "subclasses": len({u for s in S.values() for u in s["units"]}), "names": len(g["NAMES"]),
            "checked": check, "max_diff": worst}
    if verbose:
        print("froze %s: %d PCs of %d genes, ridge %g, wpow %g; %d sources, %d units (%d subclasses); "
              "%d targets (%s); self-check on %s max |diff| %.1e"
              % (os.path.relpath(E.OUTSIDE_RULE, ROOT), k, summ["genes"], lam, wpow, summ["sources"], summ["units"],
                 summ["subclasses"], len(TARGETS), ",".join(TARGETS), ",".join(check), worst))
    return summ


# --- report -----------------------------------------------------------------------------------------

def report(nperm=2000):
    t0 = time.time()
    tr = transfer(nperm)
    t1 = time.time()
    designs, sim, summ = similarity(nperm)
    t2 = time.time()
    fz = freeze()
    S = sources()
    divs = sorted({s["division"] for s in S.values()})
    L = ["# D-104 rehearsal outside the hypothalamus", "",
         "Produced by `python3 -I analysis/d104_rehearsal.py report`. No hypothalamic-source projection "
         "was read: Supplementary Table 2 is filtered to non-hypothalamic sources on loading, and the 19 "
         "hypothalamic test subclasses are removed from every source (their gene means are used, genes only).", "",
         "Data: %d non-hypothalamic sources in %d divisions, %d units (subclass in a source with >= %d labelled neurons)."
         % (len(S), len(divs), sum(len(s["units"]) for s in S.values()), MIN_NEURONS), "",
         "## 1. Transfer: the ridge rule, leave one division out", "",
         "Ridge from gene principal components (2,000 most variable genes) to the double-centred log-share "
         "profile, sample weights (n / mean n)^0.5. Each division is held out in turn, with its subclasses "
         "removed from training. r_unit: weighted correlation of unit-row-normalised prediction and measured "
         "profile; p from %d within-source permutations. Selection metric: mean r_unit over divisions (%.0f s)."
         % (nperm, t1 - t0), "",
         "| PCs | lambda | " + " | ".join(divs) + " | mean | subcortical | positive |",
         "|---|---|" + "---|" * len(divs) + "---|---|---|"]
    for k, lam, o in tr:
        sm = o["_summary"]
        L.append("| %d | %g | " % (k, lam) + " | ".join("%+.3f (p %.4f)" % (o[d]["r_unit"], o[d]["p"]) for d in divs)
                 + " | %+.3f | %+.3f | %d/%d |" % (sm["mean_r_unit"], sm["subcortical_mean_r_unit"],
                                                  sm["positive_divisions"], sm["n_divisions"]))
    best = max(tr, key=lambda x: x[2]["_summary"]["mean_r_unit"])
    L += ["", "Best mean r_unit: PCs %d, lambda %g (%+.3f)." % (best[0], best[1], best[2]["_summary"]["mean_r_unit"]), "",
          "## 2. Similarity test (design A) with location partialled out", "",
          "Pairs of units within a source; gene similarity from d104_epiretro.fingerprints over the division's "
          "subclasses; profile similarity from d104_epiretro.profiles on targets labelled in both sexes. Units need "
          ">= 50%% of their MERFISH cells in the division. %d covariates (centroid distance, three Gaussian "
          "proximities, three intermingling scales, shared methylation clusters, six coordinate products), "
          "block = source, weights 1/(1/n_i + 1/n_j). Null: gene labels permuted within source, %d permutations (%.0f s)."
          % (sim[0][1]["D"].shape[1] if sim else 0, nperm, t2 - t1), "",
          "| division | units | pairs | r raw | z raw | p raw | r adjusted | z adjusted | p adjusted | position explains proj / genes |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for d, des, o in sim:
        L.append("| %s | %d | %d | %+.3f | %+.2f | %.4f | %+.3f | %+.2f | %.4f | %.0f%% / %.0f%% |" % (
            d, len(des["units"]), len(des["I"]), *o["raw"], *o["adj"], 100 * o["r2_proj"], 100 * o["r2_gene"]))
    L += ["", "Mean z over all divisions: %+.2f raw, %+.2f adjusted; four subcortical: %+.2f raw, %+.2f adjusted." % summ,
          "", "Sources per division:", ""]
    for d, des in designs.items():
        L.append("- %s: %s" % (d, "; ".join(des["notes"]) if des else "nothing to compare"))
    L += ["", "## 3. The frozen outside rule", "",
          "`data/interim/d104_outside_rule.npz`, loaded by `d104_epiretro.outside_rule()`: %d PCs of %d genes, "
          "ridge %g, weight power %g, fitted on all %d non-hypothalamic sources (%d units, %d subclasses). "
          "%d targets: %s. Self-check on %s: the stored rule reproduces the fitted predictions (max difference %.1e)."
          % (fz["k"], fz["genes"], fz["lam"], fz["wpow"], fz["sources"], fz["units"], fz["subclasses"],
             len(fz["targets"]), ", ".join(fz["targets"]), ", ".join(fz["checked"]), fz["max_diff"]), "",
          "Total run time %.0f s." % (time.time() - t0), ""]
    with open(REPORT, "w") as f:
        f.write("\n".join(L))
    print("wrote %s" % os.path.relpath(REPORT, ROOT))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    t0 = time.time()
    if cmd == "transfer":
        transfer(n)
    elif cmd == "similarity":
        similarity(n)
    elif cmd == "freeze":
        freeze()
    elif cmd == "report":
        report(n)
    else:
        sys.exit(__doc__)
    print("[%.0f s]" % (time.time() - t0))
