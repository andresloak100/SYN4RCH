"""D-101, ABANDONED DESIGN -- kept for the record and for its gene programmes.

The held-out-prediction test below was set aside before it touched the measured wiring.
On random wiring both of its nulls passed far more than 5% of datasets: relabelling the
gene kernel (every programme passed, p = 0.0033 at the floor, in a dry run) and the
Freedman-Lane residual relabelling (17-35% of 40 random datasets). The gain in held-out
R-squared from adding a gene kernel depends on how that kernel lines up with location,
and neither null keeps that alignment while removing a gene-wiring link. The test actually
run is analysis/d101_programmes.py, which reuses D-099's validated null. This file still
provides programme_genes() and representations().

Original question: does a wiring-relevant gene programme predict an untraced population's
wiring better than its location alone?

    python3 analysis/d101_conversions.py programs   # build the gene programmes from GO
    python3 analysis/d101_conversions.py power      # size and power on synthetic wiring
    python3 analysis/d101_conversions.py run        # the pre-registered test, run once

The user's "conversion" idea, made testable: rather than comparing populations on all
2,000 most variable genes, convert each population's expression into a programme chosen
from biology before any result -- axon guidance and cell-cell adhesion, or transcription
factors -- and ask whether that programme improves the prediction of held-out wiring
beyond what location gives.

Populations, wiring fingerprints, the overlap relation and injection positions are D-099's
(analysis/d099_gene_wiring.py), unchanged. Gene programmes come from the Gene Ontology
(data/raw/go: go-basic.obo and mgi.gaf.gz): a term and all its is_a / part_of
descendants, annotations without a NOT qualifier and with evidence other than ND.
"""
import collections
import gzip
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
import d099_gene_wiring as D  # noqa: E402

OBO = os.path.join(ROOT, "data/raw/go/go-basic.obo")
GAF = os.path.join(ROOT, "data/raw/go/mgi.gaf.gz")
PROGRAMS = collections.OrderedDict([
    ("all_genes", None),                                   # D-091's 2,000 most variable genes
    ("guidance_adhesion", ["GO:0007411", "GO:0098609"]),   # axon guidance; cell-cell adhesion (GO's
                                                           # replacement for obsolete GO:0098742)
    ("transcription_factors", ["GO:0003700"]),             # DNA-binding transcription factor activity
])
LAMBDA = 1.0            # ridge, on kernels with unit diagonal; fixed, not tuned on outcomes
N_PERM = 5000
SEED = 20261111
FAMILY_ALPHA = 0.05     # Holm over the three programmes


# --- gene programmes --------------------------------------------------------------------

def go_descendants(roots):
    children = collections.defaultdict(set)
    cur = None
    for line in open(OBO, encoding="utf-8"):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = {}
        elif line.startswith("[") and line.endswith("]"):
            cur = None
        elif cur is not None:
            if line.startswith("id: "):
                cur["id"] = line[4:]
            elif line.startswith("is_a: "):
                children[line[6:].split(" ! ")[0]].add(cur["id"])
            elif line.startswith("relationship: part_of "):
                children[line.split()[2]].add(cur["id"])
    out = set(roots)
    stack = list(roots)
    while stack:
        for c in children.get(stack.pop(), ()):
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


def programme_genes():
    terms = {k: go_descendants(v) for k, v in PROGRAMS.items() if v}
    genes = {k: set() for k in terms}
    with gzip.open(GAF, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if "NOT" in f[3] or f[6] == "ND":
                continue
            for k, ts in terms.items():
                if f[4] in ts:
                    genes[k].add(f[2])
    return genes


def representations(placed, keys):
    """Per programme, each population's mean over its placed cell types of the programme's
    genes, z-scored across HypoMap's 465 cell types, unit-normalised."""
    import d095_pair_test as T95
    z = np.load(D.GENES_NPZ, allow_pickle=True)
    means = z["means"]
    names = [str(c) for c in z["clusters"]]
    symbols = [str(g) for g in z["genes"]]
    gi = {g: i for i, g in enumerate(symbols)}
    idx = {n: i for i, n in enumerate(names)}
    out, sizes = collections.OrderedDict(), {}
    gp = programme_genes()
    for k in PROGRAMS:
        if k == "all_genes":
            F, _ = T95.atlas_matrix(D.GENES_NPZ)
        else:
            cols = sorted(gi[g] for g in gp[k] if g in gi)
            F = means[:, cols]
            sd = F.std(axis=0)
            F = (F[:, sd > 0] - F[:, sd > 0].mean(axis=0)) / sd[sd > 0]
        sizes[k] = F.shape[1]
        out[k] = np.array([D.unit(F[[idx[l] for l in placed[p]["leaf_clusters"]]].mean(axis=0)) for p in keys])
    return out, sizes


# --- the prediction and its null ------------------------------------------------------

MIN_TRAIN = 10          # a held-out population is graded only if this many others may train it


def training_sets(overlap):
    """Held-out population p is predicted from every population that is neither p nor
    overlaps p (shares cells), so no near-duplicate leaks into its own prediction. A
    population with fewer than MIN_TRAIN such others is not graded (None)."""
    n = overlap.shape[0]
    out = []
    for p in range(n):
        tr = np.array([j for j in range(n) if j != p and not overlap[p, j]], dtype=int)
        out.append(tr if len(tr) >= MIN_TRAIN else None)
    return out


def location_kernel(pos):
    d = np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(axis=2))
    h = np.median(d[np.triu_indices(len(pos), 1)])      # bandwidth from positions alone
    return np.exp(-d ** 2 / (2 * h ** 2)), h


def solve_weights(K, train, perms=None, lam=None):
    """For each held-out p: c_p (over its training set) such that the prediction is
    sum_j c_pj w_j -- kernel ridge on mean-centred wiring, mean added back. With `perms`
    (P, n) the gene part of K is relabelled: the caller passes K as a function of a perm."""
    lam = LAMBDA if lam is None else lam
    out = []
    for p, tr in enumerate(train):
        if tr is None:
            out.append(None)
            continue
        m = len(tr)
        if perms is None:
            Ks = K[np.ix_(tr, tr)][None]
            kp = K[p, tr][None]
        else:
            Ks = np.stack([K(pi)[np.ix_(tr, tr)] for pi in perms])
            kp = np.stack([K(pi)[p, tr] for pi in perms])
        v = np.linalg.solve(Ks + lam * np.eye(m)[None], kp[:, :, None])[:, :, 0]
        c = v + ((1.0 - v.sum(axis=1)) / m)[:, None]
        out.append(c)
    return out


def r2(C, train, Gw):
    """Pooled held-out R^2 for every row of C (rows = relabellings): 1 - sum_p ||w_hat - w||^2
    / sum_p ||w - mean_train(w)||^2, from the wiring Gram matrix Gw alone."""
    sse = 0.0
    sst = 0.0
    for p, (c, tr) in enumerate(zip(C, train)):
        if tr is None:
            continue
        Gtt = Gw[np.ix_(tr, tr)]
        sse = sse + np.einsum("km,mn,kn->k", c, Gtt, c) - 2 * c @ Gw[tr, p] + Gw[p, p]
        u = np.full(len(tr), 1.0 / len(tr))
        sst += u @ Gtt @ u - 2 * u @ Gw[tr, p] + Gw[p, p]
    return 1.0 - sse / sst


def relabel(region, n_perm, seed):
    return D.relabellings(region, n_perm, seed)


class FreedmanLane:
    """The null for "genes add to location" (Freedman and Lane 1983). The reduced model --
    location only, fitted in-sample -- gives fitted wiring F and residuals E. Under each
    within-region relabelling pi of the residuals, synthetic wiring F + E[pi] is scored by the
    same held-out models. The gene kernels are never permuted, so their correlation with
    location is kept: permuting the genes themselves (D-101's first draft) broke it and
    passed every programme on random wiring."""

    def __init__(self, W, KL, P, lam=None):
        lam = LAMBDA if lam is None else lam
        mean = W.mean(axis=0)
        alpha = np.linalg.solve(KL + lam * np.eye(len(KL)), W - mean)
        F = mean + KL @ alpha
        E = W - F
        self.FF, self.FE, self.EE, self.P = F @ F.T, F @ E.T, E @ E.T, P


def r2_null(c, train, fl, chunk=1000):
    """Held-out R^2 of the single-row weights `c` on F + E[pi] for every row pi of fl.P
    (row 0, the identity, is the observed wiring)."""
    K = fl.P.shape[0]
    sse, sst = np.zeros(K), np.zeros(K)
    for p, (cp, tr) in enumerate(zip(c, train)):
        if tr is None:
            continue
        cp = cp[0]
        m = len(tr)
        u = np.full(m, 1.0 / m)
        FFtt = fl.FF[np.ix_(tr, tr)]
        for a in range(0, K, chunk):
            Pc = fl.P[a:a + chunk]
            Pt, pp = Pc[:, tr], Pc[:, p]
            A = fl.FE[tr][:, Pt].transpose(1, 0, 2)                    # F_i . E_pi(j)
            Gtt = FFtt[None] + A + A.transpose(0, 2, 1) + fl.EE[Pt[:, :, None], Pt[:, None, :]]
            gtp = fl.FF[tr, p][None] + fl.FE[tr][:, pp].T + fl.FE[p][Pt] + fl.EE[Pt, pp[:, None]]
            gpp = fl.FF[p, p] + 2 * fl.FE[p, pp] + fl.EE[pp, pp]
            sse[a:a + chunk] += np.einsum("m,kmn,n->k", cp, Gtt, cp) - 2 * gtp @ cp + gpp
            sst[a:a + chunk] += np.einsum("m,kmn,n->k", u, Gtt, u) - 2 * gtp @ u + gpp
    return 1.0 - sse / sst


def gain_test(W, KL, KG, train, P):
    """Observed gain in held-out R^2 from adding the gene kernel, its Freedman-Lane null,
    and the one-sided p."""
    fl = FreedmanLane(W, KL, P)
    cL = solve_weights(KL, train)
    cLG = solve_weights(KL + KG, train)
    d = r2_null(cLG, train, fl) - r2_null(cL, train, fl)
    return d, (1 + (d[1:] >= d[0]).sum()) / (1 + len(d) - 1)


def holm(ps):
    order = np.argsort(ps)
    adj = np.empty(len(ps))
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (len(ps) - rank) * ps[i])
        adj[i] = min(1.0, running)
    return adj


# --- data ----------------------------------------------------------------------------

def load():
    placed = json.load(open(D.POPS))
    keys, G, W, n_struct = D.features(placed)
    region, overlap, cand, eligible, _ = D.structure(placed, keys)
    _, pos = D.location(placed, keys, cand)
    reps, sizes = representations(placed, keys)
    return placed, keys, W, region, overlap, pos, reps, sizes


# --- power: synthetic wiring --------------------------------------------------------------

POWER_SIMS = 100
POWER_PERM = 200
LATENT = 50


def synth(rep, region, pos, b, s_pos, rng, KL_h):
    """w = regional + noise + s_pos * smooth field over position + b * programme through a
    random projection."""
    n = len(region)
    regs = sorted(set(region))
    U = D.unitrows(rng.normal(size=(len(regs), LATENT)))
    W = U[[regs.index(r) for r in region]] + D.unitrows(rng.normal(size=(n, LATENT)))
    if s_pos:
        centres = pos[rng.integers(n, size=40)] + rng.normal(0, KL_h / 2, size=(40, 3))
        field = np.exp(-((pos[:, None, :] - centres[None]) ** 2).sum(axis=2) / (2 * (KL_h / 2) ** 2)) @ rng.normal(size=(40, LATENT))
        W = W + s_pos * D.unitrows(field)
    if b:
        lam, U2 = np.linalg.eigh(rep @ rep.T)
        C = U2 * np.sqrt(np.clip(lam, 0, None))
        W = W + b * D.unitrows(C @ rng.normal(size=(n, LATENT)))
    return D.unitrows(W)


def power():
    placed, keys, W, region, overlap, pos, reps, sizes = load()
    train = training_sets(overlap)
    KL, h = location_kernel(pos)
    graded = [t for t in train if t is not None]
    print("populations %d; graded %d (training sets of at least %d); median training set %d"
          % (len(keys), len(graded), MIN_TRAIN, int(np.median([len(t) for t in graded]))))
    print("location kernel bandwidth (median pairwise distance) %.0f um" % h)
    for k, s in sizes.items():
        print("programme %-22s %4d genes" % (k, s))
    Gw = W @ W.T
    cL = solve_weights(KL, train)
    real_loc = r2(cL, train, Gw)[0]
    print("real wiring, location only (no genes): held-out R^2 %.3f" % real_loc)
    P = relabel(region, POWER_PERM, SEED + 1)
    rng = np.random.default_rng(SEED + 2)
    first = next(iter(reps.values()))
    match = min((0.5, 1.0, 1.5, 2.0, 3.0, 4.0), key=lambda s: abs(np.mean(
        [r2(cL, train, (lambda Wx: Wx @ Wx.T)(synth(first, region, pos, 0.0, s, rng, h)))[0]
         for _ in range(40)]) - real_loc))
    print("synthetic position strength matched to the real location-only R^2: %.1f" % match)
    print("rates over %d synthetic datasets, Freedman-Lane null with %d relabellings; a pass is "
          "p < alpha/3, the strictest Holm step" % (POWER_SIMS, POWER_PERM))
    for k, rep in reps.items():
        KG = rep @ rep.T
        print("programme %s:" % k)
        settings = [("iid", 0.0, 0.0), ("regional", 0.0, 0.0), ("position", 0.0, match),
                    ("position", 0.0, 2 * match), ("genes", 0.5, match), ("genes", 1.0, match),
                    ("genes", 1.5, match)]
        for label, b, s_pos in settings:
            hits, gains = 0, []
            for _ in range(POWER_SIMS):
                if label == "iid":
                    Ws = D.unitrows(rng.normal(size=W.shape))
                else:
                    Ws = synth(rep, region, pos, b, s_pos, rng, h)
                d, p = gain_test(Ws, KL, KG, train, P)
                gains.append(d[0])
                hits += p < FAMILY_ALPHA / 3
            print("  %-9s programme %.1f, position %.1f: mean R^2 gain %+.4f   passes %.3f"
                  % (label, b, s_pos, float(np.mean(gains)), hits / POWER_SIMS))


# --- the run ------------------------------------------------------------------------

def run():
    placed, keys, W, region, overlap, pos, reps, sizes = load()
    # pinning and clean-tree checks reuse D-099's, plus the GO inputs
    hashes = D.input_hashes(placed)
    import hashlib
    for name, f in (("go_obo", OBO), ("go_gaf", GAF)):
        hashes[name] = "%s (1 file)" % hashlib.sha256(open(f, "rb").read()).hexdigest()
    differs = sorted(k for k, v in hashes.items() if PINNED.get(k) != v)
    if differs:
        sys.exit("refusing to run: inputs differ from the pinned ones: %s" % ", ".join(differs))
    dirty = D.git("status", "--porcelain", "--", "analysis", "src")
    if dirty:
        sys.exit("refusing to run: uncommitted changes to the code:\n" + dirty)
    print("commit %s; every input matches its pinned SHA-256\n" % D.git("rev-parse", "HEAD"))

    train = training_sets(overlap)
    KL, h = location_kernel(pos)
    Gw = W @ W.T
    cL = solve_weights(KL, train)
    base = r2(cL, train, Gw)[0]
    P = relabel(region, N_PERM, SEED)
    print("held-out R^2, location only: %.4f" % base)
    rows = []
    for k, rep in reps.items():
        d, p = gain_test(W, KL, rep @ rep.T, train, P)
        rows.append((k, d[0], np.median(d[1:]), np.percentile(d[1:], 95), p))
    adj = holm(np.array([r[4] for r in rows]))
    print("\nDECISION - a programme improves held-out wiring prediction beyond location "
          "(Holm-adjusted one-sided p < %.2f)" % FAMILY_ALPHA)
    for (k, d0, med, q95, p), pa in zip(rows, adj):
        print("  %-22s R^2 gain %+.4f   null median %+.4f, 95th %+.4f   p = %.4f   Holm %.4f  -> %s"
              % (k, d0, med, q95, p, pa, "PASS" if pa < FAMILY_ALPHA else "not significant"))
    print("\nSECONDARY (no decision)")
    for k, rep in reps.items():
        cG = solve_weights(rep @ rep.T, train)
        print("  %-22s held-out R^2, genes alone %.4f" % (k, r2(cG, train, Gw)[0]))
    regmean = []
    for p, tr in enumerate(train):
        if tr is None:
            regmean.append(None)
            continue
        same = [j for j in tr if region[j] == region[p]]
        regmean.append(np.full(len(tr), 0.0))
        if same:
            regmean[-1][[list(tr).index(j) for j in same]] = 1.0 / len(same)
        else:
            regmean[-1][:] = 1.0 / len(tr)
    print("  region mean of the other populations      held-out R^2 %.4f"
          % r2([None if c is None else c[None] for c in regmean], train, Gw)[0])
    for lam in (0.1, 10.0):
        b = r2(solve_weights(KL, train, lam=lam), train, Gw)[0]
        gains = ["%s %+.4f" % (k, r2(solve_weights(KL + rep @ rep.T, train, lam=lam), train, Gw)[0] - b)
                 for k, rep in reps.items()]
        print("  ridge %.1f: location only %.4f; gains %s" % (lam, b, ", ".join(gains)))


PINNED = {
    "genes_npz": "5b103fecaacc8d814966db2e2f638e120c7de2d907f29e4755ee44af5155a496 (1 file)",
    "cluster_tables": "8e4a611df5b37abbe820ff59d03f0f561690c8202abdc8d4b6bfd3fba64b6522 (1 file)",
    "populations": "bd8ecebc70b43052d57e37b643b639de6024367a3ef622b0c6279685143636ae (1 file)",
    "hy_injections": "031aa47061ec084737e5ef0128e3027ccf1d9f1dba86864be308352f55ae9fc0 (1 file)",
    "structure_graph": "b5f0d024d1df09ee18aef15093f55240f4a70b3924ee9689add612dadff80386 (1 file)",
    "projections": "a6ef25acf06a4aa11ea7d3895d3f3692c4fda60bf8cb621b40a91252b4bb6e75 (117 files)",
    "go_obo": "b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc (1 file)",
    "go_gaf": "9b7613f786f9f805531a8f6768f310220d28c82095f938e1212c86f39596ce30 (1 file)",
}


if __name__ == "__main__":
    cmd = sys.argv[1:2]
    if cmd == ["programs"]:
        gp = programme_genes()
        placed = json.load(open(D.POPS))
        _, sizes = representations(placed, sorted(placed))
        for k in PROGRAMS:
            print("%-22s %s genes in HypoMap with variance (GO annotated: %s)"
                  % (k, sizes[k], len(gp[k]) if k in gp else "n/a"))
    elif cmd == ["power"]:
        power()
    elif cmd == ["run"]:
        run()
    else:
        sys.exit(__doc__)
