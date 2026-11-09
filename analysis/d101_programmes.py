"""D-101: does converting a population's genes into a wiring-related programme reveal a
link to wiring that the all-gene fingerprint (D-099) did not?

    python3 analysis/d101_programmes.py power   # size and power on synthetic wiring
    python3 analysis/d101_programmes.py run     # the pre-registered test, run once

D-099's validated test, unchanged except for the gene similarity: the same 73 Cre-defined
populations, wiring fingerprints, overlap relation, 86 eligible within-region pairs,
location covariates and exact within-region relabelling null. Gene similarity is computed
from one of two programmes chosen from the Gene Ontology before any result
(analysis/d101_conversions.py builds them): axon guidance + cell-cell adhesion, and DNA-binding
transcription factors.

The family is three tests: D-099's all-gene test, already run (one-sided p = 0.2057), and
the two programmes. Holm over the three.

analysis/d101_conversions.py holds the held-out-prediction design this replaced. It was
abandoned before it touched the real wiring: on random wiring its two candidate nulls
(gene relabelling; Freedman-Lane residual relabelling) passed 17-100% of datasets.
"""
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
import d099_gene_wiring as D  # noqa: E402
import d101_conversions as C  # noqa: E402

PROGRAMMES = ("guidance_adhesion", "transcription_factors")
D099_P = 0.2057                     # D-099's decision-1 p, all genes, already observed
ALPHA = 0.05
N_PERM = 10000
SEED = 20261112
POWER_SIMS = 400
POWER_PERM = 1000


def holm(ps):
    order = np.argsort(ps)
    adj, running = np.empty(len(ps)), 0.0
    for rank, i in enumerate(order):
        running = max(running, (len(ps) - rank) * ps[i])
        adj[i] = min(1.0, running)
    return adj


def setup():
    placed = json.load(open(D.POPS))
    keys, G_all, W, _ = D.features(placed)
    region, overlap, cand, eligible, _ = D.structure(placed, keys)
    Dloc, pos = D.location(placed, keys, cand)
    reps, sizes = C.representations(placed, keys)
    return placed, keys, W, region, overlap, cand, eligible, Dloc, reps, sizes


def power():
    placed, keys, W, region, overlap, cand, eligible, Dloc, reps, sizes = setup()
    I = np.array([p[0] for p in cand])
    J = np.array([p[1] for p in cand])
    P = D.relabellings(region, POWER_PERM, SEED + 1)
    rng = np.random.default_rng(SEED + 2)
    print("populations %d; eligible within-region pairs %d" % (len(keys), len(eligible)))
    for k in PROGRAMMES:
        rep = reps[k]
        Sg = rep @ rep.T
        prim = D.Partial(Sg, overlap, cand, region, Dloc, P)
        tr = D.Transfer(Sg, overlap, region, P)
        ctx = D.context(placed, keys, rep, region)
        print("programme %s (%d genes); a pass needs p < %.4f (the first Holm step of three)"
              % (k, sizes[k], ALPHA / 3))
        models = [("iid", None), ("gene", 0.0), ("gene", 0.5), ("gene", 0.75), ("gene", 1.0),
                  ("position", 0.6), ("position", 1.0), ("cells", 1.0)]
        for model in models:
            hp = ht = 0
            rs = []
            for _ in range(POWER_SIMS):
                if model[0] == "iid":
                    Wx = D.unitrows(rng.normal(size=W.shape))
                else:
                    Wx = D.simulate(model, ctx, rng)
                Sw = Wx @ Wx.T
                r = prim.corr(Sw[I, J])
                rs.append(r[0])
                pp = D.one_sided_p(r) < ALPHA / 3
                hp += pp
                ht += pp and D.one_sided_p(tr.gain(Sw)) < ALPHA
            print("  %-9s %4s   mean adjusted r %+.3f   decision 1 %.3f   transfer %.3f"
                  % (model[0], "" if model[1] is None else "%.2f" % model[1], float(np.nanmean(rs)),
                     hp / POWER_SIMS, ht / POWER_SIMS))


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


def run():
    placed = json.load(open(D.POPS))
    hashes = D.input_hashes(placed)
    import hashlib
    for name, f in (("go_obo", C.OBO), ("go_gaf", C.GAF)):
        hashes[name] = "%s (1 file)" % hashlib.sha256(open(f, "rb").read()).hexdigest()
    differs = sorted(k for k, v in hashes.items() if PINNED.get(k) != v)
    if differs:
        sys.exit("refusing to run: inputs differ from the pinned ones: %s" % ", ".join(differs))
    dirty = D.git("status", "--porcelain", "--", "analysis", "src")
    if dirty:
        sys.exit("refusing to run: uncommitted changes to the code:\n" + dirty)
    print("commit %s; every input matches its pinned SHA-256\n" % D.git("rev-parse", "HEAD"))

    placed, keys, W, region, overlap, cand, eligible, Dloc, reps, sizes = setup()
    Sw = W @ W.T
    I = np.array([p[0] for p in cand])
    J = np.array([p[1] for p in cand])
    P = D.relabellings(region, N_PERM, SEED)
    rows = []
    for k in PROGRAMMES:
        Sg = reps[k] @ reps[k].T
        r = D.Partial(Sg, overlap, cand, region, Dloc, P).corr(Sw[I, J])
        g = D.Transfer(Sg, overlap, region, P).gain(Sw)
        rows.append((k, r, D.one_sided_p(r), g, D.one_sided_p(g)))
    adj = holm(np.array([D099_P] + [row[2] for row in rows]))
    print("FAMILY: D-099 all genes (p = %.4f, already run) and the two programmes; Holm-adjusted "
          "one-sided p < %.2f passes" % (D099_P, ALPHA))
    print("  all_genes              (D-099)                                        Holm %.4f" % adj[0])
    for (k, r, p, g, pg), pa in zip(rows, adj[1:]):
        ok = pa < ALPHA
        print("  %-22s r = %+.3f  null 95th %+.3f  p = %.4f  Holm %.4f  -> %s"
              % (k, r[0], np.nanpercentile(r[1:], 95), p, pa, "PASS" if ok else "not significant"))
        print("  %-22s transfer gain %+.4f  null 95th %+.4f  p = %.4f  -> %s"
              % ("", g[0], np.nanpercentile(g[1:], 95), pg,
                 "transfer supported" if ok and pg < ALPHA else
                 "transfer not shown" if ok else "not tested"))
    print("\nSECONDARY (no decision)")
    for k in PROGRAMMES:
        Sg = reps[k] @ reps[k].T
        r = D.Partial(Sg, overlap, cand, region, None, P).corr(Sw[I, J])
        print("  %-22s within region, location NOT partialled  r = %+.3f  p = %.4f"
              % (k, r[0], D.one_sided_p(r)))
        Ie = np.array([p[0] for p in eligible])
        Je = np.array([p[1] for p in eligible])
        print("  %-22s correlation with the all-gene similarity over the eligible pairs %+.3f"
              % ("", np.corrcoef(Sg[Ie, Je], (lambda G: G @ G.T)(D.features(placed)[1])[Ie, Je])[0, 1]))


if __name__ == "__main__":
    cmd = sys.argv[1:2]
    if cmd == ["power"]:
        power()
    elif cmd == ["run"]:
        run()
    else:
        sys.exit(__doc__)
