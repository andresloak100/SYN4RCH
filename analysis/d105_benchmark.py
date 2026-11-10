"""D-105: a pre-registered transfer benchmark on Epi-retro-seq (Zhou et al. 2023).

Models are trained on the projections of NON-hypothalamic source regions and then predict the
projection profiles of 19 hypothalamic Allen subclasses (HY19) over 7 targets. The predictions
are frozen first; the measured hypothalamic profiles are opened once, in `run`, and scored.

LEAKAGE RULE (absolute)
  Only `run` reads hypothalamus-SOURCE projection data (rows of Supplementary Table 2 with
  Source == "HY"). The training loader drops those rows immediately after reading the sheet,
  before any other processing, and removes the 19 test subclasses from every other source.
  `dry` never reads the hypothalamic rows: it builds one synthetic hypothalamic table with the
  real layout (19 subclasses, 7 targets x one female and one male experiment) from the
  generator in `synthetic_hy_table` and runs every line of the evaluation on it. `power` uses
  synthetic profiles with the stated marginals (constants below), never the real counts.
  Gene expression (Allen 10x cluster means) and MERFISH positions of any subclass, the 19
  included, may be used in every mode. The gene inputs (z-scores, PCA) are transductive over
  the genes of all subclasses including the 19, exactly as in the rehearsal.

MODES   python3 -I analysis/d105_benchmark.py MODE
  freeze    build every frozen model from non-hypothalamic data only; write
            data/interim/d105_models.npz, data/interim/d105_t1_seed{20..25}.pt and
            data/interim/d105_predictions.npz (every model's prediction for the 19 units x 7
            targets, spatial strata, permutation index matrix); print SHA-256 of each file.
  rehearse  the frozen configurations under leave-one-division-out on non-hypothalamic
            sources; writes reports/d105_rehearsal.md.
  power     size and power of the primary decision on synthetic hypothalamic data;
            writes reports/d105_power.md.
  dry       the whole of `run` on one synthetic hypothalamic table. NOT A RESULT.
  run       the pre-registered evaluation, executed once. Refuses if D105_QUICK is set, if
            PINNED is empty or any input SHA-256 differs from it, or if
            `git status --porcelain -- analysis src data/interim/d105_*` is not empty, or if this
            file is not the committed analysis/d105_benchmark.py (blob equal to HEAD's). It writes
            data/interim/d105_run_opened.json before the hypothalamic rows are read and refuses
            if that marker exists, unless D105_RERUN=1 (the repeat is announced in the output).
D105_QUICK=1 shrinks every long step for smoke tests and writes d105_quick_* files instead.
"""
import hashlib
import json
import os
import re
import subprocess
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np  # noqa: E402

# ------------------------------------------------------------------------------- constants
SEED = 20261110
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = {
    "st2": "data/raw/epiretro/st2_cells.xlsx",
    "stats": "data/raw/abc/precomputed_stats.h5",
    "terms": "data/raw/abc/term_counts.csv",
    "genes_10x": "data/raw/abc/genes_10x.csv",
    "merfish": "data/raw/abc/merfish_cells.csv",
}
QUICK = bool(os.environ.get("D105_QUICK"))
TAG = "d105_quick_" if QUICK else "d105_"
INTERIM = os.path.join(ROOT, "data/interim")
MODELS = os.path.join(INTERIM, TAG + "models.npz")
PREDS = os.path.join(INTERIM, TAG + "predictions.npz")
T1_PT = os.path.join(INTERIM, TAG + "t1_seed%d.pt")
REHEARSAL_MD = os.path.join(ROOT, "reports/d105_rehearsal.md")
POWER_MD = os.path.join(ROOT, "reports/d105_power.md")
RUN_MARKER = os.path.join(INTERIM, "d105_run_opened.json")   # written by `run` before the HY rows are read
SCRIPT_REL = "analysis/d105_benchmark.py"

HY19 = ["PVH-SO-PVa Otp Glut", "AVPV-MEPO-SFO Tbr1 Glut", "DMH Hmx2 Gaba", "LHA-AHN-PVH Otp Trh Glut",
        "MPN-MPO-PVpo Hmx2 Glut", "PVpo-VMPO-MPN Hmx2 Gaba", "AHN Onecut3 Gaba", "AHN-RCH-LHA Otp Fezf1 Glut",
        "AHN-SBPV-PVHd Pdrm12 Gaba", "BST-MPN Six3 Nrgn Gaba", "PH-ant-LHA Otp Bsx Glut",
        "PVHd-SBPV Six3 Prox1 Gaba", "VMH Nr5a1 Glut", "DMH-LHA Gsx1 Gaba", "DMH-LHA Vgll2 Glut",
        "LHA Barhl2 Glut", "PVHd-DMH Lhx6 Gaba", "STN-PSTN Pitx2 Glut", "ZI Pax6 Gaba"]
PRIMARY = ("PFC", "MOB", "STR", "TH", "SC", "P", "MY")      # hypothalamic targets, both sexes
SEXES = ("female", "male")
DIVISION = {"MOp": "isocortex", "ACA": "isocortex", "SSp": "isocortex", "AUD": "isocortex", "PTLp": "isocortex",
            "AI": "isocortex", "VIS": "isocortex", "RSP": "isocortex",
            "MOB": "olfactory", "AON": "olfactory", "PIRa": "olfactory", "PIRp": "olfactory",
            "ENT": "hippocampal", "CAa": "hippocampal", "CAp": "hippocampal", "DGa": "hippocampal", "DGp": "hippocampal",
            "AMY": "amygdala-pallidum", "PAL": "amygdala-pallidum", "STR": "amygdala-pallidum",
            "THl": "thalamus", "THm": "thalamus", "THp": "thalamus",
            "VTA": "midbrain", "MRN": "midbrain", "SC": "midbrain", "PAG": "midbrain", "IC": "midbrain",
            "P": "hindbrain", "MY": "hindbrain", "CBX": "cerebellum"}
SUB4 = ("amygdala-pallidum", "thalamus", "midbrain", "hindbrain")
MIN_NEURONS = 20                # labelled neurons of a subclass in a training source
PSEUDO = 0.5                    # added to every unit count in an experiment
MIDLINE_Z_MM = 5.7              # MERFISH cells folded onto one hemisphere

# gene features (devkit.features; pretrained/pretrain.py)
N_GENES = 2000                  # most variable genes across NAMES (models K, B)
WB_GENES = 4000                 # most variable genes across all Allen clusters (models D, O2, T1)
WB_KEEP = 64                    # PCA loadings stored
WB_PCS = 32                     # PCs used by D and O2
B_PCS = 10
# heads (all: one coefficient vector per target shared by all sources, sqrt(n) unit weights)
LAM = {"D": 0.3, "O2": 0.3, "K": 1.0, "B": 1000.0, "A": 0.1, "T1": 1.0}
KT_SELF, KT_LEN_MM = 0.5, 2.0   # O2: Kt = 0.5 I + 0.5 exp(-d / 2 mm)
A_MIN_MM = 0.05                 # A: distances floored at 50 um
# T1 (r2_finetune/DECLARATION.txt)
T1_SEEDS = (20, 21) if QUICK else (20, 21, 22, 23, 24, 25)
DAE_DIM, DAE_HIDDEN, DAE_EPOCHS, DAE_BATCH, DAE_LR, DAE_WD = 32, 256, (2 if QUICK else 40), 512, 2e-3, 1e-4
FT_STEPS, FT_MU, FT_RANK, FT_LR = (20 if QUICK else 500), 10.0, 4, 1e-3
TORCH_THREADS = 2               # as pretrain.py / pre.py / common.py
# decision
N_PERM = 10000
N_BOOT = 2000                   # interval: r_D +/- t(0.975, 18) * sqrt(19/18) * bootstrap SD (units resampled)
STRATUM_MIN, STRATUM_MAX = 3, 4 # size-constrained Ward strata of the 19 HY centroids (3-4 units each)
HOLM = (0.025, 0.05)            # D-104 p >= 0.025 or unknown -> 0.025; D-104 p < 0.025 -> 0.05
CEIL_LOW, CEIL_HIGH, CI_REF, CI_NEG = 0.6, 0.8, 0.38, 0.25
SECONDARY = ("O2", "K", "B", "A", "T1")
MAXT = ("D", "O2", "K", "B")
MODEL_ORDER = ("D", "O2", "K", "B", "A", "T1")
# rehearsal targets (round 2: mean over 7 divisions, subcortical-4 mean)
EXPECTED = {"D": (0.228, 0.161), "O2": (0.247, 0.177), "K": (0.176, 0.087), "B": (0.161, None),
            "T1": (0.280, 0.211)}
REHEARSE_NPERM = 50 if QUICK else 2000
# power (synthetic hypothalamus; no hypothalamic count is read)
POWER_TOTAL = 958               # labelled neurons of the 19 units over 14 experiments (spec-supplied, as p4_power.py)
POWER_N_PER_UNIT = POWER_TOTAL / 19.0   # expected labelled neurons per subclass (relative sizes only)
POWER_LOGN_SD = 0.6
POWER_EXP_TOTAL = (25, 150)     # neurons per experiment, uniform integers, then rescaled to sum to POWER_TOTAL
POWER_SD_TRUE, POWER_TAU = 0.85, 0.35   # as protocol/p4_power.py
POWER_RHOS = (0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4)
POWER_NSIM = 20 if QUICK else 400
POWER_NSIM_NULL = 20 if QUICK else 1000
POWER_NBOOT = 200 if QUICK else 1000
POWER_LS_MM = (0.25, 0.5, 1.0)  # spatially smooth truths (Gaussian process over HY centroids)

# Filled after `freeze` from its printout; an empty dict makes `run` refuse.
PINNED = {
    'st2': '4855044e38ea6aa782a5a6cc9379d7da5e8bfd83c42dacccdd6b2bf42de7b1cf',
    'stats': 'b21ca985652fb25f9608f99005139a40757133a76fbe845ae5b175c5c26a447b',
    'terms': '653606461801393cd0a7708a084e3f9c989ffb2ac696ac30629f6ef6180e76c5',
    'genes_10x': 'ec5a211c23fab5cce8a06f01bb0a68af087fb979b647ab42d635201f98de2f69',
    'merfish': '1e3ae23cc3f8d3d5839cc8a222798ff812fc66f4bd3c0fc8fca0c0279329874e',
    'models': '413cfeb2a0628a0da385aa6362b9cf5a07a8c01c0e2f8a21b040bfc6dfe09c46',
    'predictions': '77d4ddd71cccbb34736996e41f12741e02ae06b86c37a76864a3ff4ea4b53a51',
    't1_seed20': '51671d38612ccbc084231e93218983db0d9b4108f5e686b37eda6943d49d80e8',
    't1_seed21': 'c735c5f239d0266d701af27bb7b144f37aab99a69bf3d794957734f8e26999ca',
    't1_seed22': 'f9533fc28297709a4b42c270d8d61e8254f00004a194a6e67369d16d865f67fc',
    't1_seed23': '2fff8cfb6328db9fc50cd31e25131dd6777636ceec9a8e27dd6daf6f0932a11e',
    't1_seed24': 'a9975f45b6d8dfa29f3ed96cae67ed8685fe6d28e6172fa88896878840e768f9',
    't1_seed25': '31245b5b8fdee05cb18bb4bf6e22ccce89618607842864bef87818b2ef9d5270',
}


def strip(name):
    return re.sub(r"^\d+\s+", "", str(name))


def path(rel):
    return os.path.join(ROOT, rel)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def say(*a):
    print(*a, flush=True)


# ------------------------------------------------------------------- projection tables
def read_st2_training():
    """Supplementary Table 2 WITHOUT the hypothalamic source. The Source == "HY" rows are
    dropped right after the sheet is read, before anything else touches the table."""
    import pandas as pd
    df = pd.read_excel(path(RAW["st2"]), sheet_name="TableS2", header=1)
    df = df[df["Source"] != "HY"].copy()
    if (df["Source"] == "HY").sum() != 0:
        sys.exit("leakage guard: hypothalamic rows survived the training filter")
    k = df[(df["PassTargetFilter"] == 1) & df["Subclass"].notna()].copy()
    k["sub"] = k["Subclass"].map(strip)
    return k


def read_st2_hypothalamus():
    """THE ONLY READER OF HYPOTHALAMIC-SOURCE ROWS. Called by `run` only."""
    import pandas as pd
    if MODE != "run":
        sys.exit("leakage guard: hypothalamic rows may be read in `run` only")
    df = pd.read_excel(path(RAW["st2"]), sheet_name="TableS2", header=1)
    return df[df["Source"] == "HY"].copy()


def sources(k):
    """devkit.sources(): per non-HY source, units = subclasses with >= MIN_NEURONS labelled
    neurons (HY19 removed first), experiments, their targets and sexes, counts[unit, exp]."""
    import pandas as pd
    if (k["Source"] == "HY").sum() != 0:
        sys.exit("leakage guard: hypothalamic rows reached sources()")
    k = k[~k["sub"].isin(HY19)]
    out = {}
    for s, g in k.groupby("Source"):
        vc = g["sub"].value_counts()
        units = sorted(u for u, n in vc.items() if n >= MIN_NEURONS)
        exps = sorted(g["Experiment"].unique())
        first = g.groupby("Experiment").first()
        et = np.array([first.loc[e, "Target"] for e in exps])
        ex = np.array([first.loc[e, "Sex"] for e in exps])
        ct = pd.crosstab(g["sub"], g["Experiment"]).reindex(index=units, columns=exps, fill_value=0).values.astype(int)
        targets = sorted(set(et))
        if len(units) >= 2 and len(targets) >= 2:
            out[s] = {"units": units, "exps": exps, "exp_target": et, "exp_sex": ex, "counts": ct,
                      "targets": targets, "n": ct.sum(axis=1), "division": DIVISION[s]}
    if "HY" in out:
        sys.exit("leakage guard: HY among the training sources")
    return out


def dc(A):
    return A - A.mean(axis=1, keepdims=True) - A.mean(axis=0, keepdims=True) + A.mean()


def log_shares(c, exp_target, targets):
    """shares among the units of each experiment (PSEUDO added), averaged over a target's
    experiments, logged: units x targets."""
    q = (c + PSEUDO) / (c.sum(axis=0) + PSEUDO * c.shape[0])
    return np.log(np.column_stack([q[:, np.asarray(exp_target) == t].mean(axis=1) for t in targets]))


def src_profile(s):
    return log_shares(s["counts"], s["exp_target"], s["targets"])


def unitrows(A):
    n = np.linalg.norm(A, axis=-1, keepdims=True)
    return A / np.where(n > 0, n, 1.0)


def wcorr(a, b, w):
    """devkit._wcorr: uncentred correlation of a*w and b*w."""
    a = a * w
    b = b * w
    den = np.sqrt((a @ a) * (b @ b))
    return float(a @ b / den) if den > 0 else 0.0


# --------------------------------------------------------------------- gene features
def load_abc():
    """Allen whole-brain 10x cluster means (5,322 clusters), cluster -> subclass map, weights."""
    import h5py
    import pandas as pd
    f = h5py.File(path(RAW["stats"]), "r")
    tree = json.loads(f["taxonomy_tree"][()].decode())
    c2r = json.loads(f["cluster_to_row"][()].decode())
    genes = json.loads(f["col_names"][()].decode())
    sym = pd.read_csv(path(RAW["genes_10x"])).set_index("gene_identifier")["gene_symbol"]
    sym = np.array([str(sym.get(g, g)) for g in genes])
    cells = pd.read_csv(path(RAW["terms"])).set_index("label")["number_of_cells"]
    X = f["sum"][:, :]                                  # cluster means (n_cells is 1 throughout)
    n = X.shape[0]
    subkeys = list(tree["CCN20230722_SUBC"])
    subnames = [strip(tree["name_mapper"]["CCN20230722_SUBC"][k]["name"]) for k in subkeys]
    sup2sub = {s: a for a, k in enumerate(subkeys) for s in tree["CCN20230722_SUBC"][k]}
    subc = np.full(n, -1)
    w = np.zeros(n)
    for sp in tree["CCN20230722_SUPT"]:
        for c in tree["CCN20230722_SUPT"][sp]:
            r = c2r[c]
            subc[r] = sup2sub[sp]
            w[r] = cells[c]
    assert (subc >= 0).all() and len(set(subnames)) == len(subnames)
    return {"X": X, "sym": sym, "genes": np.array(genes), "subc": subc, "w": w, "subnames": subnames}


def subclass_means(abc, names):
    """mean log2(CPM+1) per subclass: cluster means weighted by cells (d104 subclass_means)."""
    out = np.zeros((len(names), abc["X"].shape[1]))
    for a, s in enumerate(names):
        m = abc["subc"] == abc["subnames"].index(s)
        out[a] = abc["w"][m] @ abc["X"][m] / abc["w"][m].sum()
    return out


def gene_zscores(means, sym):
    """devkit.features: z-scored N_GENES most variable genes across the named subclasses."""
    var = means.var(axis=0)
    var[np.array([str(g).startswith(("mt-", "Rpl", "Rps")) for g in sym])] = -1
    keep = np.argsort(var)[::-1][:N_GENES]
    F = means[:, keep]
    mu = F.mean(axis=0)
    sd = F.std(axis=0)
    sd[sd == 0] = 1
    return (F - mu) / sd, keep, mu, sd


def norm_rows(M):
    """evalkit.norm: centre over NAMES, one global scale (mean squared row norm 1)."""
    m = M.mean(axis=0)
    M = M - m
    c = np.sqrt((M ** 2).sum(axis=1).mean())
    return M / c, m, c


def wb_pca(abc, names):
    """pretrained/pretrain.py (1): whole-brain PCA of the WB_GENES most variable genes across all
    Allen clusters; subclass embedding = cell-weighted mean of its clusters' scores."""
    import torch
    X, sym = abc["X"], abc["sym"]
    var = X.var(axis=0)
    var[np.array([g.startswith(("mt-", "Rpl", "Rps")) for g in sym])] = -1
    keep = np.sort(np.argsort(var)[::-1][:WB_GENES])
    Z = X[:, keep]
    mu = Z.mean(0)
    sd = Z.std(0)
    sd[sd == 0] = 1
    Z = ((Z - mu) / sd).astype(np.float32)
    zc_mean = Z.mean(0)
    Zc = torch.tensor(Z - zc_mean)
    ev, V = torch.linalg.eigh(Zc.T @ Zc)
    s2 = np.clip(ev.numpy()[::-1], 0, None)
    load = torch.flip(V, [1])[:, :WB_KEEP]
    Ec = (Zc @ load).numpy().astype(float)
    emb = np.zeros((len(names), WB_KEEP))
    for a, u in enumerate(names):
        m = abc["subc"] == abc["subnames"].index(u)
        emb[a] = abc["w"][m] @ Ec[m] / abc["w"][m].sum()
    return {"Z": Z, "keep": keep, "mu": mu, "sd": sd, "zc_mean": zc_mean, "loadings": load.numpy(),
            "emb": emb, "var_explained": float(s2[:WB_PCS].sum() / s2.sum())}


# --------------------------------------------------------------------------- positions
TGT_STRUCT = {"PFC": ["PL", "ILA", "ACAd", "ACAv", "ORBl", "ORBm", "ORBvl", "MOs"], "MOB": ["MOB"],
              "SC": ["SCm", "SCs"], "VTA": ["VTA"], "MOp": ["MOp"], "ACA": ["ACAd", "ACAv"],
              "SSp": "SSp", "AUD": "AUD", "VIS": "VIS", "RSP": "RSP", "PTLp": ["VISa", "VISrl"],
              "AI": ["AId", "AIv", "AIp"], "PIR": ["PIR"], "ENT": ["ENTl", "ENTm"], "PAG": ["PAG"],
              "IC": ["IC"], "MRN": ["MRN"],
              "AMY": ["BLA", "BMA", "LA", "PA", "CEA", "MEA", "AAA", "IA", "BA", "COAa", "COAp"],
              "AUDp": ["AUDp"], "VISp": ["VISp"]}
TGT_DIV = {"STR": ["STR"], "TH": ["TH"], "P": ["P"], "MY": ["MY"], "PAL": ["PAL"], "HY": ["HY"],
           "HPF": ["HPF"], "CB": ["CB"], "CBX": ["CB"], "CBN": ["CB"], "MB": ["MB"]}
DIV_CCF = {"isocortex": ["Isocortex"], "olfactory": ["OLF"], "hippocampal": ["HPF"],
           "amygdala-pallidum": ["CTXsp", "STR", "PAL"], "thalamus": ["TH"], "midbrain": ["MB"],
           "hindbrain": ["P", "MY"]}
SRC_STRUCT = {"MOp": ["MOp"], "ACA": ["ACAd", "ACAv"], "SSp": "SSp", "AUD": "AUD", "PTLp": ["VISa", "VISrl"],
              "AI": ["AId", "AIv", "AIp"], "VIS": ["VISp", "VISl", "VISal", "VISam", "VISpm", "VISpl", "VISli", "VISpor"],
              "RSP": ["RSPd", "RSPv", "RSPagl"], "MOB": ["MOB"], "AON": ["AON"], "PIRa": ["PIR"], "PIRp": ["PIR"],
              "ENT": ["ENTl", "ENTm"], "CAa": ["CA1", "CA2", "CA3"], "CAp": ["CA1", "CA2", "CA3"], "DGa": ["DG"],
              "DGp": ["DG"], "AMY": ["BLA", "BMA", "LA", "PA", "CEA", "MEA", "AAA", "IA", "BA", "COAa", "COAp"],
              "SC": ["SCm", "SCs"], "PAG": ["PAG"], "IC": ["IC"], "MRN": ["MRN"],
              "VTA": ["VTA", "SNc", "SNr", "IF", "RL", "CLI", "IPN", "RR", "PN", "RN"]}
SRC_DIV = {"PAL": ["PAL"], "STR": ["STR"], "THl": ["TH"], "THm": ["TH"], "THp": ["TH"], "P": ["P"], "MY": ["MY"]}


class Merfish:
    """Allen MERFISH cells in CCF (um), folded at the midline."""

    def __init__(self):
        import pandas as pd
        cols = ["subclass", "x_ccf", "y_ccf", "z_ccf", "parcellation_division", "parcellation_structure"]
        m = pd.read_csv(path(RAW["merfish"]), usecols=cols)
        sub = m["subclass"].map(strip).values
        self.xyz = np.column_stack([m["x_ccf"], m["y_ccf"], np.abs(m["z_ccf"] - MIDLINE_Z_MM)]) * 1000.0
        self.div = m["parcellation_division"].values.astype(str)
        self.st = m["parcellation_structure"].values.astype(str)
        self.allst = sorted(set(self.st))
        self.idx = pd.Series(np.arange(len(sub))).groupby(sub).indices
        del m

    def _structs(self, v):
        return [x for x in self.allst if x.startswith(v)] if isinstance(v, str) else v

    def unit_centroid(self, sub, src):
        """rehearse_pos.unit_cells: cells in the source's parcel, else its division, else all;
        here ALL such cells (the rehearsal subsampled 400)."""
        ix = self.idx.get(sub, np.array([], int))
        if src in SRC_STRUCT:
            sel = np.isin(self.st[ix], self._structs(SRC_STRUCT[src]))
        else:
            sel = np.isin(self.div[ix], SRC_DIV[src])
        level = "parcel"
        if sel.sum() < 20:
            sel, level = np.isin(self.div[ix], DIV_CCF[DIVISION[src]]), "division"
        if sel.sum() < 20:
            sel, level = np.ones(len(ix), bool), "brain"
        return self.xyz[ix[sel]].mean(axis=0), level

    def target_centroid(self, t):
        """location/pos.py target_xyz."""
        sel = None
        cbn = [x for x in self.allst if x in ("FN", "IP", "DN", "CBN")]
        if t == "CBN" and cbn and np.isin(self.st, cbn).sum() > 50:
            sel = np.isin(self.st, cbn)
        elif t in TGT_DIV:
            sel = np.isin(self.div, TGT_DIV[t])
        if t not in TGT_DIV:
            sel = np.isin(self.st, self._structs(TGT_STRUCT[t]))
        assert sel.sum() > 50, (t, int(sel.sum()))
        return self.xyz[sel].mean(axis=0)

    def hy_cells(self, sub):
        ix = self.idx.get(sub, np.array([], int))
        return self.xyz[ix[self.div[ix] == "HY"]]


# ------------------------------------------------------------------------------ heads
class Heads:
    """Shared-target ridge heads (design/model_ridge.py, evalkit.ridge, r2_objective.ls) and the
    proximity model A (location/common.py, loc = prox)."""

    def __init__(self, targets, feats, cen, tcen, RT):
        self.T = list(targets)
        self.TI = {t: a for a, t in enumerate(self.T)}
        self.F = feats                       # model -> {subclass: vector}
        self.cen = cen                       # (source, unit) -> centroid in mm ("HY", unit) included
        self.tcen = tcen                     # target -> centroid in mm
        self.RT = RT                         # Kt^(1/2) over self.T

    def tmap(self, targets):
        M = np.zeros((len(targets), len(self.T)))
        for a, t in enumerate(targets):
            M[a, self.TI[t]] = 1.0
        return M - M.mean(axis=0, keepdims=True)

    def G(self, model, units):
        G = np.array([self.F[model][u] for u in units])
        return G - G.mean(axis=0)

    def fit_ls(self, model, train, lam, R=None):
        NT = len(self.T)
        R_ = np.eye(NT) if R is None else R
        A, y = [], []
        for s in train.values():
            Y = dc(src_profile(s))
            w = np.repeat((s["n"] / s["n"].mean()) ** 0.5, Y.shape[1])
            A.append(np.einsum("ij,tk->itjk", self.G(model, s["units"]), self.tmap(s["targets"]) @ R_)
                     .reshape(Y.size, -1) * w[:, None])
            y.append(Y.ravel() * w)
        A = np.vstack(A)
        y = np.concatenate(y)
        b = np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ y)
        return b.reshape(-1, NT) @ R_

    def predict_ls(self, model, B, test):
        return {k: self.G(model, s["units"]) @ B @ self.tmap(s["targets"]).T for k, s in test.items()}

    def prox(self, src, units, targets):
        X = np.array([self.cen[(src, u)] for u in units])
        T = np.array([self.tcen[t] for t in targets])
        d = np.sqrt(((X[:, None, :] - T[None, :, :]) ** 2).sum(-1))
        return dc(-np.log(np.maximum(d, A_MIN_MM)))

    def fit_A(self, train):
        a, y = [], []
        for src, s in train.items():
            Y = dc(src_profile(s))
            w = np.repeat((s["n"] / s["n"].mean()) ** 0.5, Y.shape[1])
            a.append(self.prox(src, s["units"], s["targets"]).ravel() / np.sqrt(LAM["A"]) * w)
            y.append(Y.ravel() * w)
        a = np.concatenate(a)
        y = np.concatenate(y)
        return float((a @ y) / (a @ a + 1.0) / np.sqrt(LAM["A"]))

    def predict_A(self, c, test):
        return {k: c * self.prox(k, s["units"], s["targets"]) for k, s in test.items()}

    def fit(self, model, train):
        """-> frozen parameters of a deterministic model."""
        if model == "A":
            return self.fit_A(train)
        if model == "O2":
            return self.fit_ls("D", train, LAM["O2"], R=self.RT)
        return self.fit_ls(model, train, LAM[model])

    def predict(self, model, par, test):
        if model == "A":
            return self.predict_A(par, test)
        return self.predict_ls("D" if model == "O2" else model, par, test)


def target_kernel_root(tcen, targets):
    TC = np.array([tcen[t] for t in targets])
    DT = np.sqrt(((TC[:, None] - TC[None]) ** 2).sum(-1))
    KT = KT_SELF * np.eye(len(targets)) + (1 - KT_SELF) * np.exp(-DT / KT_LEN_MM)
    e, V = np.linalg.eigh(KT)
    return KT, (V * np.sqrt(e)) @ V.T


# --------------------------------------------------------------------------------- T1
def torch_setup():
    import torch
    torch.set_num_threads(TORCH_THREADS)
    return torch


def pretrain_dae(Z, seed):
    """r2_finetune/pre.py train_dae(32, seed): masked autoencoder on Allen cluster profiles only."""
    torch = torch_setup()
    nn = torch.nn
    Zt = torch.tensor(Z)
    n, ng = Zt.shape

    def mask(x, g):
        p = 0.15 + 0.15 * torch.rand(x.shape[0], 1, generator=g)
        m = torch.rand(x.shape, generator=g) < p
        return x.masked_fill(m, 0.0), m

    torch.manual_seed(seed)
    g = torch.Generator().manual_seed(seed)
    enc = nn.Sequential(nn.Linear(ng, DAE_HIDDEN), nn.GELU(), nn.Linear(DAE_HIDDEN, DAE_DIM))
    dec = nn.Sequential(nn.Linear(DAE_DIM, DAE_HIDDEN), nn.GELU(), nn.Linear(DAE_HIDDEN, ng))
    opt = torch.optim.AdamW(list(enc.parameters()) + list(dec.parameters()), lr=DAE_LR, weight_decay=DAE_WD)
    sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, DAE_EPOCHS)
    curve = []
    for ep in range(DAE_EPOCHS):
        perm = torch.randperm(n, generator=g)
        tot = 0.0
        for i in range(0, n, DAE_BATCH):
            x = Zt[perm[i:i + DAE_BATCH]]
            xm, m = mask(x, g)
            rec = dec(enc(xm))
            loss = ((rec - x) ** 2)[m].mean() + 0.1 * ((rec - x) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += float(loss.detach()) * len(x)
        sch.step()
        curve.append(tot / n)
    enc.eval()
    return enc, np.array(curve)


def make_enc():
    torch = torch_setup()
    nn = torch.nn
    return nn.Sequential(nn.Linear(WB_GENES, DAE_HIDDEN), nn.GELU(), nn.Linear(DAE_HIDDEN, DAE_DIM))


def cluster_agg(abc_subc, abc_w, subnames, names):
    rows = [np.where(abc_subc == subnames.index(u))[0] for u in names]
    used = np.concatenate(rows)
    Agg = np.zeros((len(names), len(used)), dtype=np.float32)
    o = 0
    for a, r in enumerate(rows):
        Agg[a, o:o + len(r)] = abc_w[r] / abc_w[r].sum()
        o += len(r)
    return used, Agg


def t1_base(enc, Zu, Agg):
    """frozen embedding E0 and normalised hidden activations Hn over NAMES (run_ft.py)."""
    torch = torch_setup()
    with torch.no_grad():
        Hc = torch.nn.functional.gelu(Zu @ enc[0].weight.T + enc[0].bias)
        Hraw = (Agg @ Hc).double()
        Eraw = (Agg @ enc[2](Hc)).double()

    def normc(M):
        m = M.mean(0, keepdim=True)
        return m, torch.sqrt(((M - m) ** 2).sum(1).mean())
    m0, c0 = normc(Eraw)
    mh, ch = normc(Hraw)
    return (Eraw - m0) / c0, (Hraw - mh) / ch


class T1:
    """LoRA rank-4 fine-tune of the last encoder layer with the shared ridge head
    (r2_finetune common.py / run_ft.py, variant T1)."""

    def __init__(self, heads, names):
        self.h = heads
        self.NI = {u: a for a, u in enumerate(names)}

    def prep(self, train):
        torch = torch_setup()
        T64 = lambda a: torch.tensor(np.asarray(a, dtype=np.float64))
        out = []
        for k, s in sorted(train.items()):
            Y = dc(src_profile(s))
            w = (s["n"] / s["n"].mean()) ** 0.5
            out.append((torch.tensor([self.NI[u] for u in s["units"]]), T64(self.h.tmap(s["targets"])), T64(Y),
                        T64(w)[:, None]))
        return out

    @staticmethod
    def sse(E, Bh, P):
        tot = 0.0
        for idx, Tc, Y, w in P:
            G = E[idx]
            G = G - G.mean(0, keepdim=True)
            tot = tot + ((w * (G @ Bh @ Tc.T - Y)) ** 2).sum()
        return tot

    def ridge_B(self, E, lam, P):
        torch = torch_setup()
        A, y = [], []
        for idx, Tc, Y, w in P:
            G = E[idx]
            G = G - G.mean(0, keepdim=True)
            A.append((torch.einsum("ij,tk->itjk", G, Tc) * w[:, :, None, None]).reshape(len(idx) * Tc.shape[0], -1))
            y.append((Y * w).reshape(-1))
        A = torch.cat(A)
        y = torch.cat(y)
        b = torch.linalg.solve(A.T @ A + lam * torch.eye(A.shape[1], dtype=A.dtype), A.T @ y)
        return b.reshape(E.shape[1], len(self.h.T))

    def fit(self, E0, Hn, train, seed):
        torch = torch_setup()
        P = self.prep(train)
        lam = LAM["T1"]
        torch.manual_seed(1000 * seed + 7)
        U = torch.zeros(DAE_DIM, FT_RANK, dtype=torch.float64, requires_grad=True)
        Vv = (torch.randn(DAE_HIDDEN, FT_RANK, dtype=torch.float64) / 16.0).requires_grad_(True)
        embed = lambda: E0 + Hn @ Vv @ U.T
        Bh = self.ridge_B(E0, lam, P).clone().requires_grad_(True)
        opt = torch.optim.Adam([{"params": [Bh, U, Vv], "lr": FT_LR}])
        sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, FT_STEPS)
        sse0 = None
        for it in range(FT_STEPS):
            E = embed()
            data = self.sse(E, Bh, P)
            loss = data + lam * (Bh ** 2).sum() + FT_MU * ((Vv @ U.T) ** 2).sum()
            if it == 0:
                sse0 = float(data.detach())
            opt.zero_grad()
            loss.backward()
            opt.step()
            sch.step()
        with torch.no_grad():
            E = embed()
            diag = {"sse0": sse0, "sse1": float(self.sse(E, Bh, P)), "rel_dE": float((E - E0).norm() / E0.norm())}
        return E.detach(), Bh.detach(), U.detach(), Vv.detach(), diag

    def predict(self, E, Bh, test):
        torch = torch_setup()
        out = {}
        for k, s in test.items():
            G = E[torch.tensor([self.NI[u] for u in s["units"]])]
            G = G - G.mean(0, keepdim=True)
            out[k] = (G @ Bh @ torch.tensor(self.h.tmap(s["targets"])).T).numpy()
        return out


# --------------------------------------------------------------------- spatial strata
def strata(X):
    """Size-constrained Ward clustering of the unit centroids into strata of STRATUM_MIN to
    STRATUM_MAX units (the spec's 3-4; the p5 rule round(19/3.5) Ward clusters gave [3, 8, 2, 6]).
    1. Agglomerate as Ward does (merge the pair with the smallest increase in within-cluster
       sum of squares, |A||B|/(|A|+|B|) * |cA - cB|^2), but only pairs whose merged size is at
       most STRATUM_MAX; stop when no such pair is left.
    2. Every unit of a stratum still smaller than STRATUM_MIN moves, in unit order, to the
       nearest (by centroid) stratum that has room. Ties break on the lowest index."""
    cl = [[a] for a in range(len(X))]
    while True:
        best = None
        for i in range(len(cl)):
            for j in range(i + 1, len(cl)):
                if len(cl[i]) + len(cl[j]) > STRATUM_MAX:
                    continue
                ci, cj = X[cl[i]].mean(0), X[cl[j]].mean(0)
                cost = len(cl[i]) * len(cl[j]) / (len(cl[i]) + len(cl[j])) * ((ci - cj) ** 2).sum()
                if best is None or cost < best[0] - 1e-15:
                    best = (cost, i, j)
        if best is None:
            break
        _, i, j = best
        cl[i] = sorted(cl[i] + cl[j])
        del cl[j]
    small = [c for c in cl if len(c) < STRATUM_MIN]
    cl = [c for c in cl if len(c) >= STRATUM_MIN]
    for a in sorted(u for c in small for u in c):
        room = [i for i, c in enumerate(cl) if len(c) < STRATUM_MAX]
        if not room:
            sys.exit("strata: no stratum has room for unit %d" % a)
        d = [((X[cl[i]].mean(0) - X[a]) ** 2).sum() for i in room]
        cl[room[int(np.argmin(d))]].append(a)
    cl = sorted(sorted(c) for c in cl)
    lab = np.empty(len(X), int)
    for s, c in enumerate(cl):
        lab[c] = s
    sizes = np.bincount(lab)
    if sizes.min() < STRATUM_MIN or sizes.max() > STRATUM_MAX:
        sys.exit("strata: sizes %s outside %d-%d" % (sizes.tolist(), STRATUM_MIN, STRATUM_MAX))
    return lab


def perms_within(lab, nperm, rng):
    out = np.tile(np.arange(len(lab)), (nperm, 1))
    for b in range(nperm):
        for c in np.unique(lab):
            ix = np.where(lab == c)[0]
            out[b, ix] = ix[rng.permutation(len(ix))]
    return out


# ----------------------------------------------------------------------------- freeze
def build_world():
    """Everything the models need, from non-hypothalamic projections, genes and positions."""
    t0 = time.time()
    tab = read_st2_training()
    S = sources(tab)
    targets = sorted({t for s in S.values() for t in s["targets"]})
    assert set(PRIMARY) <= set(targets), set(PRIMARY) - set(targets)
    names = sorted({u for s in S.values() for u in s["units"]} | set(HY19))
    say("training sources %d (%s), units %d, neurons %d, targets %d, NAMES %d  [%.0fs]" % (
        len(S), ",".join(sorted(S)), sum(len(s["units"]) for s in S.values()), sum(int(s["n"].sum()) for s in S.values()),
        len(targets), len(names), time.time() - t0))
    return S, targets, names


def freeze():
    torch = torch_setup()
    t0 = time.time()
    S, targets, names = build_world()
    abc = load_abc()
    say("Allen clusters %d x genes %d  [%.0fs]" % (abc["X"].shape + (time.time() - t0,)))
    means = subclass_means(abc, names)
    Zg, gkeep, gmu, gsd = gene_zscores(means, abc["sym"])
    U, s, _ = np.linalg.svd(Zg, full_matrices=False)
    FK, k_m, k_c = norm_rows(U * s)
    FB = U[:, :B_PCS] * s[:B_PCS]
    wb = wb_pca(abc, names)
    FD, d_m, d_c = norm_rows(wb["emb"][:, :WB_PCS])
    say("gene features: K %d PCs of %d genes; B %d PCs; D %d of %d whole-brain PCs (var explained %.3f)  [%.0fs]" % (
        FK.shape[1], N_GENES, B_PCS, WB_PCS, WB_KEEP, wb["var_explained"], time.time() - t0))

    mer = Merfish()
    cen, levels = {}, {}
    for src, sd_ in S.items():
        for u in sd_["units"]:
            c, lv = mer.unit_centroid(u, src)
            cen[(src, u)] = c / 1000.0
            levels[lv] = levels.get(lv, 0) + 1
    hy_n = []
    for u in HY19:
        x = mer.hy_cells(u)
        assert len(x) > 0, u
        cen[("HY", u)] = x.mean(axis=0) / 1000.0
        hy_n.append(len(x))
    tcen = {t: mer.target_centroid(t) / 1000.0 for t in targets}
    del mer
    say("centroids: training units by level %s; HY19 MERFISH cells in HY %s  [%.0fs]" % (levels, hy_n, time.time() - t0))
    KT, RT = target_kernel_root(tcen, targets)

    feats = {"D": dict(zip(names, FD)), "K": dict(zip(names, FK)), "B": dict(zip(names, FB))}
    H = Heads(targets, feats, cen, tcen, RT)
    test = {"HY": {"units": list(HY19), "targets": list(PRIMARY)}}
    par, pred = {}, {}
    for m in ("D", "O2", "K", "B", "A"):
        par[m] = H.fit(m, S)
        pred[m] = H.predict(m, par[m], test)["HY"]
    say("A coefficient on the double-centred -log distance: %+.4f" % par["A"])

    # T1: pretrain the masked autoencoder per seed (genes only), then LoRA fine-tune on all 27 sources
    used, Agg = cluster_agg(abc["subc"], abc["w"], abc["subnames"], names)
    Zu = torch.tensor(wb["Z"][used])
    Aggt = torch.tensor(Agg)
    t1 = T1(H, names)
    t1_pred, t1_info, t1_B, t1_U, t1_V, t1_files = [], {}, [], [], [], []
    for seed in T1_SEEDS:
        enc, curve = pretrain_dae(wb["Z"], seed)
        fpt = T1_PT % seed
        torch.save(enc.state_dict(), fpt)
        t1_files.append(fpt)
        E0, Hn = t1_base(enc, Zu, Aggt)
        E, Bh, Uu, Vv, diag = t1.fit(E0, Hn, S, seed)
        t1_pred.append(t1.predict(E, Bh, test)["HY"])
        t1_B.append(Bh.numpy()); t1_U.append(Uu.numpy()); t1_V.append(Vv.numpy())
        t1_info[seed] = {"pretrain_loss_first": float(curve[0]), "pretrain_loss_last": float(curve[-1]),
                         "pretrain_curve": [round(float(x), 5) for x in curve], **diag}
        say("T1 seed %d: pretrain loss %.4f -> %.4f; fine-tune SSE %.1f -> %.1f, relative embedding change %.3f  [%.0fs]" % (
            seed, curve[0], curve[-1], diag["sse0"], diag["sse1"], diag["rel_dE"], time.time() - t0))
    pred["T1"] = np.mean(t1_pred, axis=0)

    hyX = np.array([cen[("HY", u)] for u in HY19])
    lab = strata(hyX)
    perms = perms_within(lab, N_PERM, np.random.default_rng(SEED))
    say("spatial strata (size-constrained Ward, %d): sizes %s; %s" % (lab.max() + 1, np.bincount(lab).tolist(),
        "; ".join(", ".join(HY19[a] for a in np.where(lab == s)[0]) for s in range(lab.max() + 1))))

    config = {"seed": SEED, "lambda": LAM, "n_genes": N_GENES, "wb_genes": WB_GENES, "wb_pcs": WB_PCS, "b_pcs": B_PCS,
              "kt": [KT_SELF, KT_LEN_MM], "a_min_mm": A_MIN_MM, "min_neurons": MIN_NEURONS, "pseudo": PSEUDO,
              "t1": {"seeds": list(T1_SEEDS), "dae": [DAE_DIM, DAE_HIDDEN, DAE_EPOCHS, DAE_BATCH, DAE_LR, DAE_WD],
                     "ft": [FT_STEPS, FT_MU, FT_RANK, FT_LR], "torch_threads": TORCH_THREADS, "info": t1_info},
              "quick": QUICK, "training_sources": sorted(S), "n_perm": N_PERM, "stratum_sizes": [STRATUM_MIN, STRATUM_MAX]}
    np.savez(MODELS, names=np.array(names), targets=np.array(targets), primary=np.array(PRIMARY), hy19=np.array(HY19),
             config=np.array(json.dumps(config)),
             gene_keep=gkeep, gene_symbols=abc["sym"][gkeep], gene_mu=gmu, gene_sd=gsd,
             F_K=FK, K_centre=k_m, K_scale=k_c, F_B=FB,
             wb_keep=wb["keep"], wb_symbols=abc["sym"][wb["keep"]], wb_mu=wb["mu"], wb_sd=wb["sd"],
             wb_zc_mean=wb["zc_mean"], wb_loadings=wb["loadings"], wb_emb=wb["emb"], F_D=FD, D_centre=d_m, D_scale=d_c,
             B_D=par["D"], B_O2=par["O2"], B_K=par["K"], B_B=par["B"], A_coef=np.array(par["A"]),
             KT=KT, RT=RT, target_centroids_mm=np.array([tcen[t] for t in targets]),
             unit_centroid_keys=np.array(["%s|%s" % k for k in cen]), unit_centroids_mm=np.array(list(cen.values())),
             T1_B=np.array(t1_B), T1_U=np.array(t1_U), T1_V=np.array(t1_V), T1_seeds=np.array(T1_SEEDS),
             T1_cluster_rows=used)
    np.savez(PREDS, hy19=np.array(HY19), primary=np.array(PRIMARY),
             **{"pred_" + m: pred[m] for m in MODEL_ORDER}, pred_T1_seeds=np.array(t1_pred),
             hy_centroids_mm=hyX, hy_merfish_cells=np.array(hy_n), strata=lab, perms=perms.astype(np.int16),
             perm_seed=np.array(SEED))
    say("\nfrozen predictions, unit-row cosine between models over the 19 units (mean):")
    for a in MODEL_ORDER:
        say("  %-3s " % a + " ".join("%+.2f" % (unitrows(dc(pred[a])) * unitrows(dc(pred[b]))).sum(1).mean() for b in MODEL_ORDER))
    say("\nSHA-256 of written files:")
    for p in [MODELS, PREDS] + t1_files:
        say("  %s  %s" % (sha256(p), os.path.relpath(p, ROOT)))
    say("\nPINNED for `run` (inputs and frozen files):")
    hs = input_hashes()
    say("PINNED = {\n" + "".join("    %r: %r,\n" % kv for kv in hs.items()) + "}")
    say("freeze done  [%.0fs]" % (time.time() - t0))


# ---------------------------------------------------------------------------- rehearse
def folds(S):
    """devkit.folds: leave one division out; training sources lose every held-out subclass."""
    for d in sorted({s["division"] for s in S.values()}):
        test_src = [k for k in sorted(S) if S[k]["division"] == d]
        held = {u for k in test_src for u in S[k]["units"]}
        train = {}
        for k in sorted(S):
            if S[k]["division"] == d:
                continue
            keep = [a for a, u in enumerate(S[k]["units"]) if u not in held]
            if len(keep) >= 2:
                train[k] = dict(S[k], units=[S[k]["units"][a] for a in keep], counts=S[k]["counts"][keep],
                                n=S[k]["n"][keep])
        test = {k: {"units": list(S[k]["units"]), "targets": list(S[k]["targets"]), "division": d} for k in test_src}
        yield d, train, test


def lodo_score(S, pred, nperm, seed=0):
    """devkit.lodo scoring of stored predictions: r_unit per held-out division, permutation p."""
    rng = np.random.default_rng(seed)
    out = {}
    divs = sorted({s["division"] for s in S.values()})
    for d in divs:
        srcs = [k for k in sorted(S) if S[k]["division"] == d]
        nn = np.concatenate([S[k]["n"] for k in srcs]).astype(float)
        P, O, W = [], [], []
        for k in srcs:
            o = dc(src_profile(S[k]))
            P.append(unitrows(dc(np.asarray(pred[k], float))))
            O.append(unitrows(o))
            W.append(np.sqrt(S[k]["n"] / nn.mean())[:, None] * np.ones((1, o.shape[1])))
        cat = lambda L: np.concatenate([x.ravel() for x in L])
        r = wcorr(cat(P), cat(O), cat(W))
        null = np.empty(nperm)
        for b in range(nperm):
            null[b] = wcorr(cat([x[rng.permutation(len(x))] for x in P]), cat(O), cat(W))
        out[d] = (r, float((1 + (null >= r).sum()) / (1 + nperm)))
    rs = np.array([out[d][0] for d in divs])
    return out, float(rs.mean()), float(np.mean([out[d][0] for d in SUB4]))


def load_frozen():
    z = np.load(MODELS, allow_pickle=True)
    names = list(z["names"])
    targets = list(z["targets"])
    cen = {tuple(k.split("|", 1)): v for k, v in zip(z["unit_centroid_keys"], z["unit_centroids_mm"])}
    tcen = dict(zip(targets, z["target_centroids_mm"]))
    feats = {"D": dict(zip(names, z["F_D"])), "K": dict(zip(names, z["F_K"])), "B": dict(zip(names, z["F_B"]))}
    return z, names, targets, Heads(targets, feats, cen, tcen, z["RT"])


def rehearse():
    torch = torch_setup()
    t0 = time.time()
    S, targets, names = build_world()
    z, names_f, targets_f, H = load_frozen()
    assert names_f == names and targets_f == targets, "training world differs from the frozen one"
    preds = {m: {} for m in ("D", "O2", "K", "B", "A")}
    for d, train, test in folds(S):
        for m in preds:
            preds[m].update(H.predict(m, H.fit(m, train), test))
    say("deterministic models refitted per fold  [%.0fs]" % (time.time() - t0))
    # T1: frozen encoders from freeze, fine-tuned inside each fold
    abc = load_abc()
    used, Agg = cluster_agg(abc["subc"], abc["w"], abc["subnames"], names)
    assert np.array_equal(used, z["T1_cluster_rows"])
    keep = z["wb_keep"]
    Zu = torch.tensor(((abc["X"][used][:, keep] - z["wb_mu"]) / z["wb_sd"]).astype(np.float32))
    del abc
    t1 = T1(H, names)
    seeds = [int(s) for s in z["T1_seeds"]]
    t1_seed_preds, diags = [], []
    for seed in seeds:
        enc = make_enc()
        enc.load_state_dict(torch.load(T1_PT % seed))
        enc.eval()
        E0, Hn = t1_base(enc, Zu, torch.tensor(Agg))
        ps = {}
        for d, train, test in folds(S):
            E, Bh, _, _, dg = t1.fit(E0, Hn, train, seed)
            ps.update(t1.predict(E, Bh, test))
            diags.append(dg)
        t1_seed_preds.append(ps)
        say("T1 seed %d fine-tuned in 7 folds  [%.0fs]" % (seed, time.time() - t0))
    preds["T1"] = {k: np.mean([p[k] for p in t1_seed_preds], axis=0) for k in S}

    divs = sorted({s["division"] for s in S.values()})
    rows, summ = [], {}
    for m in MODEL_ORDER:
        o, m7, s4 = lodo_score(S, preds[m], REHEARSE_NPERM)
        summ[m] = (m7, s4, o)
    single = [lodo_score(S, p, 0)[1:] for p in t1_seed_preds]
    sm = np.array([x[0] for x in single])
    lines = ["# D-105 rehearsal: the frozen configurations, leave one division out",
             "",
             "Generated by `python3 -I analysis/d105_benchmark.py rehearse`%s. No hypothalamic-source neuron is read: "
             "the 27 non-hypothalamic sources are scored, each division held out in turn, with every held-out "
             "subclass removed from the training sources (devkit.folds). The models are refitted per fold with "
             "the frozen configurations and the frozen gene features of `data/interim/%smodels.npz`; T1 uses "
             "the frozen encoders and is fine-tuned inside each fold. Statistic: devkit r_unit; p = predicted "
             "rows permuted among a source's units (%d permutations, seed 0)." % (" (QUICK)" if QUICK else "", TAG, REHEARSE_NPERM),
             "",
             "| model | " + " | ".join(d for d in divs) + " | mean7 | sub4 | round 2 | check |",
             "|---|" + "---|" * len(divs) + "---|---|---|---|"]
    for m in MODEL_ORDER:
        m7, s4, o = summ[m]
        ex = EXPECTED.get(m)
        tol = 0.005 if m != "T1" else 2 * float(sm.std(ddof=1) if len(sm) > 1 else 0.04) / np.sqrt(max(len(sm), 1))
        if ex is None:
            chk, exs = "nominal (expected about 0)", "-"
        else:
            ok = abs(m7 - ex[0]) <= tol and (ex[1] is None or abs(s4 - ex[1]) <= tol)
            chk = ("reproduced" if ok else "DIFFERS") + " (tolerance %.3f)" % tol
            exs = "%+.3f / %s" % (ex[0], "-" if ex[1] is None else "%+.3f" % ex[1])
        lines.append("| %s | " % m + " | ".join("%+.3f (p %.3f)" % o[d] for d in divs) + " | %+.3f | %+.3f | %s | %s |" % (m7, s4, exs, chk))
        say("%-3s " % m + " ".join("%s %+.3f(p%.3f)" % (d[:5], o[d][0], o[d][1]) for d in divs) + " | mean7 %+.3f sub4 %+.3f | round 2 %s | %s" % (m7, s4, exs, chk))
    lines += ["",
              "T1 single seeds (%s): mean7 %s; mean %+.3f, sd %.3f. Round 2: single seeds mean +0.257, sd 0.041; ensemble +0.280 / +0.211. "
              "The T1 tolerance is two standard errors of the seed mean (seed sd / sqrt(seeds) x 2): the encoders are "
              "retrained from scratch, so the ensemble differs from round 2 by seed noise."
              % (",".join(map(str, seeds)), " ".join("%+.3f" % x for x in sm), sm.mean(), sm.std(ddof=1) if len(sm) > 1 else 0.0),
              "T1 fine-tune diagnostics (mean over folds x seeds): SSE %.1f -> %.1f, relative embedding change %.3f."
              % (np.mean([d_["sse0"] for d_ in diags]), np.mean([d_["sse1"] for d_ in diags]), np.mean([d_["rel_dE"] for d_ in diags])),
              "",
              "Deviation: model A's training-unit centroids use all of a subclass's MERFISH cells in the source's parcel "
              "(the round-1 location code subsampled 400 cells with a seeded generator); A has no round-2 value to reproduce.",
              ""]
    say("T1 single seeds: " + " ".join("%+.3f" % x for x in sm) + "  mean %+.3f sd %.3f" % (sm.mean(), sm.std(ddof=1) if len(sm) > 1 else 0))
    if QUICK:
        say("QUICK: report not written")
    else:
        open(REHEARSAL_MD, "w").write("\n".join(lines))
        say("wrote", os.path.relpath(REHEARSAL_MD, ROOT))
    say("rehearse done  [%.0fs]" % (time.time() - t0))


# --------------------------------------------------------------------------- evaluation
def hy_layout(hy):
    """Counts of the 19 units in the 14 primary experiments (shares among the 19)."""
    import pandas as pd
    hy = hy.copy()
    for col in ("Source", "Target", "Sex", "Experiment"):     # whitespace and type only; labels must still match
        hy[col] = hy[col].map(lambda v: str(v).strip() if (isinstance(v, str) or col == "Experiment") and pd.notna(v) else v)
    k = hy[(hy["Source"] == "HY") & (hy["PassTargetFilter"] == 1) & hy["Subclass"].notna()].copy()
    k["sub"] = k["Subclass"].map(lambda v: strip(v).strip())
    k = k[k["Target"].isin(PRIMARY)]
    exps = sorted(k["Experiment"].unique())
    first = k.groupby("Experiment").first()
    et = np.array([first.loc[e, "Target"] for e in exps])
    ex = np.array([first.loc[e, "Sex"] for e in exps])
    problems = []
    if len(exps) != 2 * len(PRIMARY):
        problems.append("%d primary experiments, expected %d" % (len(exps), 2 * len(PRIMARY)))
    for t in PRIMARY:
        if sorted(ex[et == t]) != sorted(SEXES):
            problems.append("target %s: sexes %s" % (t, list(ex[et == t])))
    if problems:
        sys.exit("refusing to score: layout differs from the pre-registered one: " + "; ".join(problems))
    k = k[k["sub"].isin(HY19)]
    c = pd.crosstab(k["sub"], k["Experiment"]).reindex(index=HY19, columns=exps, fill_value=0).values.astype(int)
    return {"counts": c, "exp_target": et, "exp_sex": ex, "exps": exps}


def measured(L):
    O = dc(log_shares(L["counts"], L["exp_target"], PRIMARY))
    n = L["counts"].sum(axis=1).astype(float)
    sexp = {}
    for sx in SEXES:
        m = L["exp_sex"] == sx
        sexp[sx] = dc(log_shares(L["counts"][:, m], L["exp_target"][m], PRIMARY))
    return O, n, sexp


def unit_weights(n):
    return np.sqrt(n / n.mean())


def r_of(P, Orow, w):
    """devkit r_unit for one stack: P, Orow unit-row normalised (units x targets), w per unit."""
    return wcorr(P.ravel(), Orow.ravel(), np.repeat(w, P.shape[1]))


def r_perm(P, Orow, w, perms):
    a = P[perms] * w[None, :, None]
    b = Orow * w[:, None]
    den = np.sqrt((a * a).sum((1, 2)) * (b * b).sum())
    return np.where(den > 0, (a * b[None]).sum((1, 2)) / np.where(den > 0, den, 1), 0.0)


def semi_partial_perm(PX, PA, Orow, w, perms=None):
    """plib.partial_r: X's prediction with A's regressed out (weighted, one slope), vs measured."""
    X = PX[None] if perms is None else PX[perms]
    W = w[None, :, None]
    xa = (X * W * PA[None] * W).sum((1, 2))
    aa = ((PA * w[:, None]) ** 2).sum()
    res = X - PA[None] * (xa / aa if aa > 0 else 0.0)[:, None, None]
    a = res * W
    b = Orow * w[:, None]
    den = np.sqrt((a * a).sum((1, 2)) * (b * b).sum())
    out = np.where(den > 0, (a * b[None]).sum((1, 2)) / np.where(den > 0, den, 1), 0.0)
    return out if perms is not None else float(out[0])


def pval(null, obs):
    return float((1 + (null >= obs - 1e-12).sum()) / (1 + len(null)))


def bootstrap_ci(P, Orow, n, nboot, rng):
    """95% interval of r_D from a bootstrap over units: r +/- t(0.975, m-1) * sqrt(m/(m-1)) * SD of
    the bootstrap draws. The plain percentile interval covered only about 0.87-0.92 in simulation
    (19 units); this one about 0.90-0.94 (reports/d105_power.md). Draws whose resampled units all
    have n = 0 are dropped (NaN); -> (lo, hi, number of NaN draws)."""
    from scipy.stats import t as tdist
    m = len(n)
    idx = rng.integers(0, m, size=(nboot, m))
    nb = n[idx]
    with np.errstate(invalid="ignore", divide="ignore"):
        w = np.sqrt(nb / nb.mean(axis=1, keepdims=True))
        a = P[idx] * w[:, :, None]
        b = Orow[idx] * w[:, :, None]
        den = np.sqrt((a * a).sum((1, 2)) * (b * b).sum((1, 2)))
        r = np.where(den > 0, (a * b).sum((1, 2)) / np.where(den > 0, den, 1), np.nan)
    bad = int((~np.isfinite(r)).sum())
    r0 = r_of(P, Orow, unit_weights(n)) if n.sum() > 0 else float("nan")
    if nboot - bad < 2 or not np.isfinite(r0):
        return float("nan"), float("nan"), bad
    half = tdist.ppf(0.975, m - 1) * np.sqrt(m / (m - 1.0)) * np.nanstd(r, ddof=1)
    return r0 - half, r0 + half, bad


def ceiling(sexp, w):
    """between-sex reliability in the r_unit form (same unit weights as the test), Spearman-Brown
    corrected; also the plain Pearson correlation of the two double-centred profiles (descriptive).
    It estimates the reliability of the measured profile; the largest attainable r_D is about its
    square root."""
    r = r_of(unitrows(sexp["female"]), unitrows(sexp["male"]), w)
    a, b = sexp["female"].ravel(), sexp["male"].ravel()
    plain = float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else float("nan")
    return r, (2 * r / (1 + r) if r > -1 else float("nan")), plain


def decide(rD, pD, rA, ceil_sb, ci):
    """-> (primary verdict at each Holm step, reading if not pass). A non-finite ceiling or
    interval reads INSUFFICIENT EVIDENCE."""
    beats = rD > max(rA, 0.0)
    step = {a: bool(pD < a and beats) for a in HOLM}
    if not (np.isfinite(ceil_sb) and np.isfinite(ci[0]) and np.isfinite(ci[1])):
        reading = "INSUFFICIENT EVIDENCE"
    elif ceil_sb < CEIL_LOW or (ci[0] <= CI_REF <= ci[1]):
        reading = "INSUFFICIENT EVIDENCE"
    elif ceil_sb >= CEIL_HIGH and ci[1] < CI_NEG:
        reading = "NEGATIVE"
    else:
        reading = "NOT SHOWN"
    return step, reading, beats


def evaluate(hy, label):
    """Every statistic of the pre-registered evaluation; `hy` holds hypothalamic-source rows
    (real in `run`, synthetic in `dry`). The layout and the primary decision are printed first;
    each secondary item is computed and printed on its own, so a failure there cannot hide D1."""
    pr = np.load(PREDS, allow_pickle=True)
    if list(pr["hy19"]) != HY19 or tuple(pr["primary"]) != PRIMARY:
        sys.exit("refusing to score: frozen predictions do not match HY19 / PRIMARY")
    perms = pr["perms"].astype(int)
    lab = pr["strata"]
    pred = {m: pr["pred_" + m] for m in MODEL_ORDER}
    L = hy_layout(hy)
    O, n, sexp = measured(L)
    w = unit_weights(n)
    Orow = unitrows(O)
    P = {m: unitrows(dc(pred[m])) for m in MODEL_ORDER}

    # ---- layout and primary decision
    r = {m: r_of(P[m], Orow, w) for m in ("D", "A")}
    p = {"D": pval(r_perm(P["D"], Orow, w, perms), r["D"])}
    r_fm, ceil_sb, r_fm_plain = ceiling(sexp, w)
    ci = bootstrap_ci(P["D"], Orow, n, N_BOOT, np.random.default_rng(SEED + 1))
    step, reading, beats = decide(r["D"], p["D"], r["A"], ceil_sb, ci)
    if step[0.025]:
        verdict = "PASS"
    elif step[0.05]:
        verdict = "PASS only if D-104's Decision 1 p < 0.025; otherwise %s" % reading
    else:
        verdict = reading
    out = ["# D-105 result: transfer of gene-to-projection models to the hypothalamus" if label == "run" else
           "# DRY RUN - NOT A RESULT (synthetic hypothalamic table)", ""]
    out += ["## Layout", "",
            "19 units (HY19), %d experiments: %s." % (len(L["exps"]), "; ".join(
                "%s %s %s (%d of the 19)" % (e, t, sx, c) for e, t, sx, c in zip(L["exps"], L["exp_target"], L["exp_sex"], L["counts"].sum(0)))),
            "Neurons of the 19 units over the 14 experiments: %d (per unit %s)." % (n.sum(), " ".join("%d" % x for x in n)),
            "Spatial strata (fixed in freeze): sizes %s." % np.bincount(lab).tolist()]
    if (n == 0).any():
        out.append("WARNING: units with no labelled neuron (weight 0): %s." % ", ".join(u for u, x in zip(HY19, n) if x == 0))
    out += ["", "## Primary decision (D1)", "",
            "r_unit(D) = %+.4f; within-strata permutation p = %.4f (%d permutations, seed %d, one-sided; ties within 1e-12 "
            "count as at least as large)." % (r["D"], p["D"], len(perms), SEED),
            "r_unit(A) = %+.4f; r_D > max(r_A, 0): %s." % (r["A"], "yes" if beats else "no"),
            "Holm with D-104 Decision 1 (family-wise 0.05): pass at the Holm step means p < 0.025 if D-104's p is larger "
            "or not yet known; p < 0.05 if D-104's p < 0.025.",
            "  threshold 0.025 (D-104 p >= 0.025 or unknown): %s" % ("PASS" if step[0.025] else "not pass"),
            "  threshold 0.05  (D-104 p <  0.025):            %s" % ("PASS" if step[0.05] else "not pass"),
            "Noise ceiling (between-sex reliability, no genes; female-only vs male-only profile in the r_unit form with the "
            "test's unit weights): r %+.4f, Spearman-Brown %+.4f. This is a reliability; the largest attainable r_D is "
            "about its square root (%.2f). Plain Pearson of the two double-centred profiles (descriptive): %+.4f." % (
                r_fm, ceil_sb, np.sqrt(ceil_sb) if ceil_sb > 0 else float("nan"), r_fm_plain),
            "95%% interval of r_D (units resampled, %d draws, seed %d; r_D +/- t(0.975, 18) x sqrt(19/18) x bootstrap SD): "
            "[%+.4f, %+.4f]%s." % (N_BOOT, SEED + 1, ci[0], ci[1], "" if not ci[2] else "; %d draws undefined, dropped" % ci[2]),
            "Readings if not PASS: INSUFFICIENT EVIDENCE if the ceiling < %.1f, the interval contains %.2f, or either is "
            "undefined; NEGATIVE if the ceiling >= %.1f and the interval's upper bound < %.2f; otherwise NOT SHOWN." % (
                CEIL_LOW, CI_REF, CEIL_HIGH, CI_NEG),
            "", "READING: %s" % verdict, ""]
    say("\n".join(out))

    # ---- secondary items, each on its own
    res = {"r": r, "p": p, "ceiling": ceil_sb, "ci": ci, "reading": verdict}

    def secondary(name, fn):
        try:
            say("\n".join(fn()))
        except Exception as e:      # a secondary failure must not stop the report
            say("SECONDARY %s FAILED: %r\n" % (name, e))

    null = {}

    def table():
        from scipy.stats import binom  # noqa: F401  (import checked here, used below)
        for m in MODEL_ORDER:
            r[m] = r_of(P[m], Orow, w)
            null[m] = r_perm(P[m], Orow, w, perms)
            p[m] = pval(null[m], r[m])
        maxnull = np.max([null[m] for m in MAXT], axis=0)
        pfam = {m: pval(maxnull, r[m]) for m in MAXT}
        top = {m: (np.argmax(dc(pred[m]), axis=1) == np.argmax(O, axis=1)) for m in MODEL_ORDER}
        res["top"] = top
        lines = ["## Secondary (no decision)", "",
                 "| model | r_unit | p (within strata) | max-T family p | top-target hits |",
                 "|---|---|---|---|---|"]
        for m in MODEL_ORDER:
            lines.append("| %s | %+.4f | %.4f%s | %s | %d/19 = %.2f |" % (
                m, r[m], p[m], " (a)" if m == "A" else "", "%.4f" % pfam[m] if m in pfam else "-", top[m].sum(), top[m].mean()))
        best = max(r[m] for m in MAXT)
        lines += ["",
                  "(a) A's permutation p is not calibrated when the true profile is a smooth function of position (A is itself "
                  "smooth in position; its size under smooth truths is in reports/d105_power.md); only r_A enters the decision.",
                  "max-T over {D, O2, K, B}: family p of the best model %.4f (max observed %+.4f). Max-T can exceed its "
                  "nominal size under smooth position-only truths (driven by K and B; sizes in reports/d105_power.md)." % (
                      pval(maxnull, best), best),
                  "Top target = the target with the largest value of the double-centred profile, for prediction and measurement.", ""]
        return lines
    secondary("model table and max-T", table)

    def semi():
        sp = semi_partial_perm(P["D"], P["A"], Orow, w)
        sp_null = semi_partial_perm(P["D"], P["A"], Orow, w, perms)
        return ["Semi-partial r of D after removing A's prediction: %+.4f, p = %.4f (same permutations, D's rows only)." % (
            sp, pval(sp_null, sp))]
    secondary("semi-partial", semi)

    def toptarget():
        from scipy.stats import binom
        top = np.argmax(dc(pred["D"]), axis=1) == np.argmax(O, axis=1)
        k = int(top.sum())
        raw = np.argmax(log_shares(L["counts"], L["exp_target"], PRIMARY), axis=1)
        return ["Top-target accuracy of D: %d/19 = %.2f vs chance 1/7 = 0.14; one-sided binomial p %.4f." % (
                    k, k / 19.0, float(binom.sf(k - 1, 19, 1 / 7.0))),
                "Descriptive: units whose measured top target is the same on raw log shares as on the double-centred "
                "profile: %d/19." % int((raw == np.argmax(O, axis=1)).sum())]
    secondary("top-target", toptarget)

    def t1seeds():
        return ["T1 single seeds r_unit: %s." % " ".join("%+.3f" % r_of(unitrows(dc(x)), Orow, w) for x in pr["pred_T1_seeds"]), ""]
    secondary("T1 seeds", t1seeds)

    def cosines():
        lines = ["## Per-unit cosines (prediction vs measured, double-centred rows)", "",
                 "| unit | n | stratum | measured top | " + " | ".join(MODEL_ORDER) + " |",
                 "|---|---|---|---|" + "---|" * len(MODEL_ORDER)]
        for a, u in enumerate(HY19):
            lines.append("| %s | %d | %d | %s | " % (u, n[a], lab[a], PRIMARY[int(np.argmax(O[a]))]) +
                         " | ".join("%+.2f" % (P[m][a] @ Orow[a]) for m in MODEL_ORDER) + " |")
        return lines + [""]
    secondary("per-unit cosines", cosines)
    if label != "run":
        say("DRY RUN - NOT A RESULT")
    return res


# --------------------------------------------------------------------------------- dry
def synthetic_hy_table(rng):
    """One synthetic hypothalamic table with the real layout: 19 units, 7 targets x (female,
    male) experiments. Unit sizes 958/19 x lognormal(0, 0.6); experiment totals uniform 25-150;
    true profile = independent N(0, 0.85^2) double-centred (no relation to any prediction);
    experiment noise 0.35. Adds non-HY19 neurons and filtered-out rows to exercise the filters."""
    import pandas as pd
    size = POWER_N_PER_UNIT * rng.lognormal(0.0, POWER_LOGN_SD, len(HY19))
    T = dc(rng.standard_normal((len(HY19), len(PRIMARY))))
    T *= POWER_SD_TRUE / T.std()
    rows = []
    for t_i, t in enumerate(PRIMARY):
        for sx in SEXES:
            e = "SYN-%s-%s" % (t, sx[0].upper())
            tot = int(rng.integers(POWER_EXP_TOTAL[0], POWER_EXP_TOTAL[1] + 1))
            pp = size * np.exp(T[:, t_i] + POWER_TAU * rng.standard_normal(len(HY19)))
            c = rng.multinomial(tot, pp / pp.sum())
            for u_i, u in enumerate(HY19):
                rows += [(e, t, sx, "%03d %s" % (u_i, u), 1)] * int(c[u_i])
            rows += [(e, t, sx, "999 Other hypothalamic Glut", 1)] * 3 + [(e, t, sx, "%03d %s" % (0, HY19[0]), 0)] * 2
    rows += [("SYN-VTA-M", "VTA", "male", "%03d %s" % (1, HY19[1]), 1)] * 5      # single-sex target: ignored
    return pd.DataFrame(rows, columns=["Experiment", "Target", "Sex", "Subclass", "PassTargetFilter"]).assign(Source="HY")


def dry():
    say("DRY RUN - NOT A RESULT")
    guards(enforce=False)
    hy = synthetic_hy_table(np.random.default_rng(SEED + 7))
    evaluate(hy, "dry")
    say("DRY RUN - NOT A RESULT")


# ---------------------------------------------------------------------------- guards
def git(*args, enforce=True):
    r = subprocess.run(["git", "-C", ROOT] + list(args), capture_output=True, text=True)
    if r.returncode != 0:
        msg = "git %s failed: %s" % (" ".join(args), r.stderr.strip())
        if enforce:
            sys.exit("refusing to run: " + msg)
        say("  (note) " + msg)
        return None
    return r.stdout.strip()


def frozen_files():
    out = [("models", MODELS), ("predictions", PREDS)]
    seeds = T1_SEEDS
    if os.path.exists(MODELS):
        seeds = [int(s) for s in np.load(MODELS, allow_pickle=True)["T1_seeds"]]
    return out + [("t1_seed%d" % s, T1_PT % s) for s in seeds]


def input_hashes():
    out = {}
    for k, p in [(k, path(v)) for k, v in RAW.items()] + frozen_files():
        out[k] = sha256(p) if os.path.exists(p) else "missing"
    return out


def guards(enforce):
    problems = []
    if QUICK:
        problems.append("D105_QUICK is set")
    hashes = input_hashes()
    if not PINNED:
        problems.append("PINNED is empty")
    else:
        differ = sorted(k for k in set(hashes) | set(PINNED) if PINNED.get(k) != hashes.get(k))
        if differ:
            problems.append("SHA-256 differs from PINNED for: " + ", ".join(differ))
        missing = sorted(k for k in set(hashes) | set(PINNED) if "missing" in (PINNED.get(k), hashes.get(k)))
        if missing:
            problems.append("missing files or 'missing' pins: " + ", ".join(missing))
    # the script itself: must be the tracked analysis/d105_benchmark.py, identical to HEAD's blob
    me = os.path.realpath(os.path.abspath(__file__))
    if me != os.path.realpath(path(SCRIPT_REL)):
        problems.append("this script is %s, not %s" % (me, SCRIPT_REL))
    blob_head = git("rev-parse", "HEAD:" + SCRIPT_REL, enforce=False)
    blob_work = git("hash-object", me, enforce=False)
    if not blob_head or blob_head != blob_work:
        problems.append("script differs from HEAD:%s (or is not committed): HEAD blob %s, working blob %s" % (
            SCRIPT_REL, blob_head, blob_work))
    dirty = git("status", "--porcelain", "--", "analysis", "src", "data/interim/d105_*", enforce=enforce)
    if dirty:
        problems.append("uncommitted changes:\n" + dirty)
    commit = git("rev-parse", "HEAD", enforce=enforce)
    say("commit %s" % commit)
    say("  script %s  sha256 %s  git blob %s (HEAD blob %s)" % (SCRIPT_REL, sha256(me), blob_work, blob_head))
    for k, v in hashes.items():
        say("  sha256 %s  %s" % (v, k))
    if problems:
        if enforce:
            sys.exit("refusing to run: " + "; ".join(problems))
        say("  (dry) `run` would refuse now: " + "; ".join(problems))


def run():
    guards(enforce=True)
    # executed once: a marker is written before the hypothalamic rows are opened
    rerun = bool(os.environ.get("D105_RERUN"))
    if os.path.exists(RUN_MARKER):
        prev = open(RUN_MARKER).read()
        if not rerun:
            sys.exit("refusing to run: the hypothalamic rows were already opened (%s):\n%s\n"
                     "Set D105_RERUN=1 to repeat; the repeat is announced in the output." % (os.path.relpath(RUN_MARKER, ROOT), prev))
        say("REPEAT RUN (D105_RERUN=1). Earlier opening recorded: %s" % prev)
    rec = {"opened_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "commit": git("rev-parse", "HEAD"),
           "hashes": input_hashes(), "repeat": rerun}
    hist = []
    if os.path.exists(RUN_MARKER):
        old = json.load(open(RUN_MARKER))
        hist = old.get("history", []) + [{k: v for k, v in old.items() if k != "history"}]
    rec["history"] = hist
    with open(RUN_MARKER, "w") as f:
        json.dump(rec, f, indent=1)
    say("hypothalamic rows opened %s (recorded in %s)" % (rec["opened_utc"], os.path.relpath(RUN_MARKER, ROOT)))
    hy = read_st2_hypothalamus()
    evaluate(hy, "run")


# ----------------------------------------------------------------------------- power
def power():
    """Size and power of D1 on synthetic hypothalamic data (protocol/p4_power.py, p5_nullsize.py).
    No hypothalamic count is read: unit sizes and experiment totals come from the constants."""
    t0 = time.time()
    pr = np.load(PREDS, allow_pickle=True)
    MD = pr["pred_D"]
    PM = {m: unitrows(dc(pr["pred_" + m])) for m in MODEL_ORDER}
    PD, PA = PM["D"], PM["A"]
    perms = pr["perms"].astype(int)
    lab = pr["strata"]
    X = pr["hy_centroids_mm"]
    rng = np.random.default_rng(SEED + 3)
    nu, nt = len(HY19), len(PRIMARY)
    zM = dc(MD)
    zM = zM / zM.std()
    tg = np.repeat(np.arange(nt), 2)
    sx = np.tile(np.array(SEXES), nt)
    full = np.array([rng.permutation(nu) for _ in range(len(perms))])   # comparator: unrestricted permutation
    d2 = ((X[:, None] - X[None]) ** 2).sum(-1)
    lss = sorted(set(POWER_LS_MM) | {2.0})
    chol = {ls: np.linalg.cholesky(np.exp(-d2 / (2 * ls ** 2)) + 1e-6 * np.eye(nu)) for ls in lss}

    def one(rho, noise):
        size = POWER_N_PER_UNIT * rng.lognormal(0.0, POWER_LOGN_SD, nu)
        et = rng.integers(POWER_EXP_TOTAL[0], POWER_EXP_TOTAL[1] + 1, 2 * nt)
        et = np.maximum(1, np.round(et * POWER_TOTAL / et.sum())).astype(int)
        E = rng.standard_normal((nu, nt)) if noise is None else chol[noise] @ rng.standard_normal((nu, nt))
        E = dc(E)
        E = E / E.std()
        T = rho * zM + np.sqrt(1 - rho ** 2) * E
        T = dc(T)
        T *= POWER_SD_TRUE / T.std()
        c = np.zeros((nu, 2 * nt), int)
        for e in range(2 * nt):
            pp = size * np.exp(T[:, tg[e]] + POWER_TAU * rng.standard_normal(nu))
            c[:, e] = rng.multinomial(et[e], pp / pp.sum())
        L = {"counts": c, "exp_target": np.array(PRIMARY)[tg], "exp_sex": sx}
        O, n, sexp = measured(L)
        w = unit_weights(n)
        Orow = unitrows(O)
        true_r = r_of(PD, unitrows(T), w)
        rD = r_of(PD, Orow, w)
        rA = r_of(PA, Orow, w)
        nulls = {m: r_perm(PM[m], Orow, w, perms) for m in MAXT + ("A",)}
        robs = {m: r_of(PM[m], Orow, w) for m in MAXT + ("A",)}
        pD = pval(nulls["D"], rD)
        pA = pval(nulls["A"], robs["A"])
        pT = pval(np.max([nulls[m] for m in MAXT], axis=0), max(robs[m] for m in MAXT))
        pF = pval(r_perm(PD, Orow, w, full), rD)
        _, csb, _ = ceiling(sexp, w)
        ci = bootstrap_ci(PD, Orow, n, POWER_NBOOT, rng)
        step, reading, _ = decide(rD, pD, rA, csb, ci)
        return dict(true=true_r, rD=rD, pD=pD, pF=pF, pA=pA, pT=pT, ceil=csb, ci=ci[:2], pass1=step[0.025],
                    pass2=step[0.05], reading=reading, n=n.sum())

    scen = [("independent noise", None)] + [("spatial noise %.2f mm" % ls, ls) for ls in lss]
    lines = ["# D-105 power: size and power of the primary decision on synthetic hypothalamic data", "",
             "Generated by `python3 -I analysis/d105_benchmark.py power`%s. No hypothalamic count is read. "
             "The prediction is model D's frozen prediction for the 19 units; model A's frozen prediction enters the "
             "PASS condition r_D > max(r_A, 0). True profile T = rho * D's double-centred prediction + sqrt(1 - rho^2) * "
             "noise, double-centred, sd %.2f; rho is a mixing weight, and the realised r_unit of D against the noiseless T "
             "('true r') is about 0.96 x rho. The noise is independent N(0, 1) entries or a Gaussian process over the 19 "
             "MERFISH centroids (length %s mm, independent per target). Counts: 7 targets x (female, male); experiment "
             "totals uniform %d-%d, then rescaled so the 14 experiments hold %d neurons of the 19 units (the total given "
             "in the spec, as in protocol/p4_power.py); unit sizes %.1f x lognormal(0, %.1f), shares proportional to size x "
             "exp(T + %.2f experiment noise). Every dataset is scored exactly as `run` does: r_unit, within-strata "
             "permutation p with the frozen %d permutations (strata sizes %s), the noise ceiling, the interval (%d bootstrap "
             "draws here, %d in `run`) and the reading. 'p full' uses unrestricted row permutations (comparator only). "
             "'CI covers mean r_D': share of datasets whose 95%% interval contains the mean of r_D over that row's datasets "
             "(a calibration check at a fixed truth)." % (
                 " (QUICK)" if QUICK else "", POWER_SD_TRUE, " / ".join("%g" % x for x in lss), POWER_EXP_TOTAL[0],
                 POWER_EXP_TOTAL[1], POWER_TOTAL, POWER_N_PER_UNIT, POWER_LOGN_SD, POWER_TAU, len(perms),
                 np.bincount(lab).tolist(), POWER_NBOOT, N_BOOT), ""]
    hdr = ("| noise | rho (mixing) | sims | true r | measured r_D | ceiling (SB) | pass p<0.025 | pass p<0.05 | "
           "p<0.025 (full perm) | max-T p<0.05 | A p<0.05 | CI covers mean r_D | INSUFFICIENT | NEGATIVE | NOT SHOWN |")
    lines += [hdr, "|" + "---|" * 15]
    say(hdr)
    for name, noise in scen:
        rhos = POWER_RHOS if noise in (None, 0.5) else (0.0,)
        for rho in rhos:
            nsim = POWER_NSIM_NULL if rho == 0.0 else POWER_NSIM
            R = [one(rho, noise) for _ in range(nsim)]
            f = lambda k: np.array([x[k] for x in R])
            rd = f("reading")
            mu = f("rD").mean()
            ci = f("ci")
            row = ("| %s | %.2f | %d | %+.3f | %+.3f | %.2f | %.3f | %.3f | %.3f | %.3f | %.3f | %.2f | %.2f | %.2f | %.2f |" % (
                name, rho, nsim, f("true").mean(), mu, np.nanmedian(f("ceil")), f("pass1").mean(), f("pass2").mean(),
                (f("pF") < 0.025).mean(), (f("pT") < 0.05).mean(), (f("pA") < 0.05).mean(),
                ((ci[:, 0] <= mu) & (mu <= ci[:, 1])).mean(), (rd == "INSUFFICIENT EVIDENCE").mean(),
                (rd == "NEGATIVE").mean(), (rd == "NOT SHOWN").mean()))
            lines.append(row)
            say(row + "  [%.0fs]" % (time.time() - t0))
    se = np.sqrt(0.025 * 0.975 / POWER_NSIM_NULL)
    lines += ["",
              "Size: the pass rates at rho = 0 should be at or below 0.025 (Monte Carlo SE %.3f with %d datasets); "
              "the pass rule also needs r_D > max(r_A, 0), which can only lower them. With spatially smooth truths, "
              "compare the within-strata column with the unrestricted-permutation column. The max-T and A columns at "
              "rho = 0 are the sizes of those secondary tests at nominal 0.05 (A's p is not calibrated under smooth "
              "truths, since A is itself a smooth function of position; only r_A enters the decision)." % (se, POWER_NSIM_NULL),
              "Interval: r_D +/- t(0.975, 18) x sqrt(19/18) x bootstrap SD. Resampling units treats them as exchangeable, "
              "so the interval still undercovers somewhat when the truth is spatially smooth.",
              "Readings are given for datasets that do not pass; a dataset that passes still has a reading computed, so the "
              "three reading columns describe the interval and ceiling rules on every dataset.",
              "The ceiling is the between-sex reliability (Spearman-Brown) of the measured profile; the attainable r_D is "
              "about its square root.",
              "Labelled neurons of the 19 units per synthetic dataset: %d (spec-supplied total)." % POWER_TOTAL,
              ""]
    if QUICK:
        say("QUICK: report not written")
    else:
        open(POWER_MD, "w").write("\n".join(lines))
        say("wrote", os.path.relpath(POWER_MD, ROOT))
    say("power done  [%.0fs]" % (time.time() - t0))


MODE = None

if __name__ == "__main__":
    modes = {"freeze": freeze, "rehearse": rehearse, "power": power, "dry": dry, "run": run}
    if len(sys.argv) != 2 or sys.argv[1] not in modes:
        sys.exit("usage: python3 -I analysis/d105_benchmark.py {%s}" % "|".join(modes))
    MODE = sys.argv[1]
    if MODE == "run" and QUICK:
        sys.exit("refusing to run: D105_QUICK is set")
    modes[MODE]()
