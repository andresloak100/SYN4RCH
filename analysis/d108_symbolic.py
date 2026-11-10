"""D-108 Phase 1, symbolic layer: an exhaustive small symbolic-regression search (written here; no
PySR/gplearn) over a dozen learned and measured variables, scored under the same nested
leave-one-division-out scheme as the discovery search.

Variables for a (unit u in source s, target t) pair, all standardised:
  z1 z2 z3  learned gene axes: the top-3 left singular vectors of D's coefficient matrix (gD32 x
            one-hot ridge, lambda 0.3) fitted on the CURRENT training sources only; z = gD32 @ U[:, :3]
  ux uy uz  measured location of the unit (MERFISH centroid of its subclass in its source)
  tx ty tz  target centroid; tg1 = first principal component of the target's gene profile
  ld        -log distance(unit location, target centroid)   (model A's feature)
  ldh       -log distance(inferred location from genes, module M1, target centroid)
Expressions (depth <= 2): leaves v; op(f(a), g(b)) with f, g in {id, sq, tanh} and op in {+, -, *, /};
op1(op2(a, b), c) with op1 in {*, /}, op2 in {+, -, *, /}. Division is protected: x / (y + 0.1 sign y).
Model: profile[u, t] = c * dc_s(expr)  (one coefficient; r_unit depends only on the sign of c,
which is fitted on the training sources). Selection: mean inner r_unit - PEN * complexity.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import d105_benchmark as b5  # noqa: E402
import d108_prims as P  # noqa: E402

PEN = 0.005
VARS = ("z1", "z2", "z3", "ux", "uy", "uz", "tx", "ty", "tz", "tg1", "ld", "ldh")
UN = ("id", "sq", "tanh")
OPS = ("+", "-", "*", "/")


def _un(f, x):
    return x if f == "id" else (x * x if f == "sq" else np.tanh(x))


def _op(o, a, b):
    if o == "+":
        return a + b
    if o == "-":
        return a - b
    if o == "*":
        return a * b
    return a / (b + 0.1 * np.where(b >= 0, 1.0, -1.0))


def catalogue():
    """names and complexities of every expression, in evaluation order."""
    names, comp = [], []
    for v in VARS:
        names.append(v)
        comp.append(1)
    U = [(f, v) for f in UN for v in VARS]
    un = lambda f, v: v if f == "id" else "%s(%s)" % (f, v)
    for o in OPS:
        for f1, a in U:
            for f2, b in U:
                names.append("(%s %s %s)" % (un(f1, a), o, un(f2, b)))
                comp.append(3 + (f1 != "id") + (f2 != "id"))
    for o1 in ("*", "/"):
        for o2 in OPS:
            for a in VARS:
                for b in VARS:
                    for c in VARS:
                        names.append("((%s %s %s) %s %s)" % (a, o2, b, o1, c))
                        comp.append(5)
    return names, np.array(comp)


def evaluate_all(V):
    """V: (nvar, nobs) -> generator of (nobs,) chunks (n_chunk, nobs) in catalogue order."""
    yield V
    Uv = np.concatenate([_un(f, V) for f in UN])            # (36, nobs)
    for o in OPS:
        yield _op(o, Uv[:, None, :], Uv[None, :, :]).reshape(-1, V.shape[1])
    for o1 in ("*", "/"):
        for o2 in OPS:
            W = _op(o2, V[:, None, :], V[None, :, :])        # (12, 12, nobs)
            yield _op(o1, W[:, :, None, :], V[None, None, :, :]).reshape(-1, V.shape[1])


class SRFold:
    """training sources -> learned axes, training signs; test division -> signed r_unit per expression."""

    def __init__(self, eng, train, test, scales):
        self.eng = eng
        f = P.Fold(eng, train, train, use_E=False)
        Bm = f.solve([("K", "gD32", "onehot")], [0.3])[0].reshape(32, len(eng.targets))
        U, s, _ = np.linalg.svd(Bm, full_matrices=False)
        zw = U[:, :3]
        tu = sorted({u for x in train.values() for u in x["units"]})
        Z = eng.F["gD32"][[eng.NI[u] for u in tu]] @ zw
        self.zw = zw / Z.std(0)
        self.sc = scales
        self.train, self.test = train, test

    def vars_of(self, k, s):
        eng = self.eng
        n, T = len(s["units"]), len(s["targets"])
        ui = [eng.NI[u] for u in s["units"]]
        z = eng.F["gD32"][ui] @ self.zw
        lm = np.array([(eng.cen[(k, u)] - eng.lm_mu) / eng.lm_sc for u in s["units"]]) / self.sc["lm"]
        ti = [eng.TI[t] for t in s["targets"]]
        tc = eng.enc["tc"][ti] / self.sc["tc"]
        tg = eng.enc["tg"][ti, 0] / self.sc["tg"]
        L = np.array([eng.cen[(k, u)] for u in s["units"]])
        Lh = eng.Lhat_mm[ui]
        Tc = eng.tcen[ti]
        ld = (-np.log(np.maximum(np.sqrt(((L[:, None] - Tc[None]) ** 2).sum(-1)), 0.05)) - self.sc["ld"][0]) / self.sc["ld"][1]
        ldh = (-np.log(np.maximum(np.sqrt(((Lh[:, None] - Tc[None]) ** 2).sum(-1)), 0.05)) - self.sc["ldh"][0]) / self.sc["ldh"][1]
        uu = lambda x: np.repeat(x[:, None], T, 1)
        tt = lambda x: np.repeat(x[None, :], n, 0)
        cols = [uu(z[:, 0]), uu(z[:, 1]), uu(z[:, 2]), uu(lm[:, 0]), uu(lm[:, 1]), uu(lm[:, 2]),
                tt(tc[:, 0]), tt(tc[:, 1]), tt(tc[:, 2]), tt(tg), ld, ldh]
        return np.stack([c.ravel() for c in cols])

    def score(self):
        """-> (training sign, held-out r_unit) per expression (signed r = sign * r)."""
        blocks, V = [], []
        o = 0
        for side, SS in (("tr", self.train), ("te", self.test)):
            for k in sorted(SS):
                s = SS[k]
                v = self.vars_of(k, s)
                blocks.append((side, k, o, len(s["units"]), len(s["targets"])))
                V.append(v)
                o += v.shape[1]
        V = np.concatenate(V, axis=1)
        obs = {}
        nn = np.concatenate([self.test[k]["n"] for k in sorted(self.test)]).astype(float)
        for side, k, *_ in blocks:
            s = (self.train if side == "tr" else self.test)[k]
            Y = b5.dc(b5.src_profile(s))
            if side == "tr":
                obs[(side, k)] = (Y, s["n"] / s["n"].mean())
            else:
                obs[(side, k)] = (b5.unitrows(Y), np.sqrt(s["n"] / nn.mean()))
        sgn, rr = [], []
        for E in evaluate_all(V):
            E = np.nan_to_num(E, nan=0.0, posinf=0.0, neginf=0.0)
            num = np.zeros(len(E))
            pp = np.zeros(len(E))
            po = np.zeros(len(E))
            oo = 0.0
            for side, k, off, n, T in blocks:
                X = E[:, off:off + n * T].reshape(len(E), n, T)
                X = X - X.mean(2, keepdims=True) - X.mean(1, keepdims=True) + X.mean((1, 2), keepdims=True)
                X[np.abs(X) < 1e-9] = 0.0           # a unit-only or target-only expression double-centres to 0
                Y, w = obs[(side, k)]
                if side == "tr":
                    num += np.einsum("bij,ij,i->b", X, Y, w)
                else:
                    nr = np.linalg.norm(X, axis=2, keepdims=True)
                    X = X / np.where(nr > 1e-9, nr, np.inf)
                    Xw = X * w[None, :, None]
                    Ow = Y * w[:, None]
                    pp += (Xw * Xw).sum((1, 2))
                    po += (Xw * Ow[None]).sum((1, 2))
                    oo += float((Ow * Ow).sum())
            den = np.sqrt(pp * oo)
            r = np.where(den > 0, po / np.where(den > 0, den, 1), 0.0)
            s_ = np.sign(num)
            sgn.append(s_)
            rr.append(s_ * r)
        return np.concatenate(sgn), np.concatenate(rr)


def scales(eng, S):
    allc = np.array([eng.cen[(k, u)] for k, s in S.items() for u in s["units"]])
    lm = (allc - eng.lm_mu) / eng.lm_sc
    ld, ldh = [], []
    for k, s in S.items():
        L = np.array([eng.cen[(k, u)] for u in s["units"]])
        Lh = eng.Lhat_mm[[eng.NI[u] for u in s["units"]]]
        Tc = eng.tcen[[eng.TI[t] for t in s["targets"]]]
        ld.append(-np.log(np.maximum(np.sqrt(((L[:, None] - Tc[None]) ** 2).sum(-1)), 0.05)).ravel())
        ldh.append(-np.log(np.maximum(np.sqrt(((Lh[:, None] - Tc[None]) ** 2).sum(-1)), 0.05)).ravel())
    ld, ldh = np.concatenate(ld), np.concatenate(ldh)
    return {"lm": lm.std(0), "tc": eng.enc["tc"].std(0), "tg": eng.enc["tg"][:, 0].std(),
            "ld": (ld.mean(), ld.std()), "ldh": (ldh.mean(), ldh.std())}


def sr_outer(args):
    """one outer fold: inner selection over the catalogue, outer score of the chosen formula and of
    the best leaf (the simplest alternative), chosen the same way."""
    import d108_discovery as X
    spec, d = args
    t0 = time.time()
    G = X.ctx()
    eng = G["eng"]
    S = X.dataset(spec)
    sc = G.setdefault("sr_scales", scales(eng, G["S"]))
    names, comp = G.setdefault("sr_cat", catalogue())
    _, tr, te = [x for x in P.folds(S) if x[0] == d][0]
    inner = np.mean([SRFold(eng, itr, ite, sc).score()[1] for _, itr, ite in P.folds(tr)], axis=0)
    pen = inner - PEN * comp
    j = int(np.argmax(pen))
    leaf = int(np.argmax(np.where(comp == 1, inner, -np.inf)))
    sg, outer = SRFold(eng, tr, te, sc).score()
    a = VARS.index("ld")
    return {"spec": spec, "division": d, "formula": names[j], "complexity": int(comp[j]), "inner": float(inner[j]),
            "r_formula": float(outer[j]), "sign": float(sg[j]), "leaf": names[leaf], "r_leaf": float(outer[leaf]),
            "r_ld": float(outer[a]), "inner_top5": [[names[i], float(inner[i]), int(comp[i])] for i in np.argsort(-pen)[:5]],
            "sec": time.time() - t0}


def run_sr(specs, label):
    from multiprocessing import get_context
    import d108_discovery as X
    G = X.ctx()
    divs = sorted({s["division"] for s in G["S"].values()})
    jobs = [(sp, d) for sp in specs for d in divs]
    out = {}
    t0 = time.time()
    with get_context("spawn").Pool(X.WORKERS) as pool:
        for i, r in enumerate(pool.imap_unordered(sr_outer, jobs)):
            out.setdefault(json.dumps(r["spec"], sort_keys=True), []).append(r)
            if (i + 1) % 14 == 0 or i + 1 == len(jobs):
                X.say("%s: %d/%d  [%.0fs]" % (label, i + 1, len(jobs), time.time() - t0))
    for k in out:
        out[k].sort(key=lambda r: r["division"])
    return out


def sr_final():
    """selection on all 27 sources (7-division LODO), coefficient by weighted least squares with a
    bootstrap over sources; returns the module card."""
    import d108_discovery as X
    G = X.ctx()
    eng, S = G["eng"], G["S"]
    sc = scales(eng, S)
    names, comp = catalogue()
    inner = np.mean([SRFold(eng, tr, te, sc).score()[1] for _, tr, te in P.folds(S)], axis=0)
    j = int(np.argmax(inner - PEN * comp))
    f = SRFold(eng, S, S, sc)
    # the formula's values on every training unit (double-centred), coefficient by weighted LS
    vals = {}
    for k in sorted(S):
        V = f.vars_of(k, S[k])
        E = None
        for ch_start, E_ in _iter_index(V, j):
            E = E_
        vals[k] = b5.dc(np.nan_to_num(E).reshape(len(S[k]["units"]), len(S[k]["targets"])))

    def coef(keys):
        num = den = 0.0
        for k in keys:
            Y, w2 = P.src_obs(S[k])
            num += float((w2[:, None] * vals[k] * Y).sum())
            den += float((w2[:, None] * vals[k] ** 2).sum())
        return num / den
    c = coef(sorted(S))
    rng = np.random.default_rng(P.SEED + 55)
    ks = sorted(S)
    boots = [coef([ks[i] for i in rng.integers(0, len(ks), len(ks))]) for _ in range(200)]
    zw = f.zw
    return {"formula": names[j], "complexity": int(comp[j]), "inner_lodo_all7": float(inner[j]), "coef": c,
            "coef_boot_sd": float(np.std(boots, ddof=1)), "coef_boot_95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
            "z_axes_gD32": zw.tolist(), "scales": {k: (np.asarray(v).tolist()) for k, v in sc.items()},
            "top10": [[names[i], float(inner[i]), int(comp[i])] for i in np.argsort(-(inner - PEN * comp))[:10]]}


def _iter_index(V, j):
    """yield the single expression j's values (walks the catalogue chunks)."""
    o = 0
    for E in evaluate_all(V):
        if o <= j < o + len(E):
            yield o, E[j - o]
            return
        o += len(E)


def mode_sr_real():
    import d108_discovery as X
    X.save("sr_real", run_sr([{"kind": "real"}], "SR real"))
    card = sr_final()
    X.save("sr_final", card)


def mode_sr_null0():
    import d108_discovery as X
    dec = X.declared()
    X.save("sr_null0", run_sr([{"kind": "N0", "rep": r} for r in range(dec["n_SR_N0"])], "SR N0"))


def mode_sr_pilot():
    import d108_discovery as X
    t0 = time.time()
    r = run_sr([{"kind": "N0", "rep": 999}], "SR pilot")
    X.save("sr_pilot", r)
    X.say("SR pilot wall %.0fs, catalogue %d expressions" % (time.time() - t0, len(catalogue()[0])))


MODES = {"sr_pilot": mode_sr_pilot, "sr_real": mode_sr_real, "sr_null0": mode_sr_null0}
