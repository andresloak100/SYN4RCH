"""Axon completion v1: prefix-conditioned donor choice + one shared frontier-graft decoder. TRAIN-only DEVELOPMENT.

Protocol: reports/axon_completion_v1_protocol.md. Output: data/interim/axon_completion_v1/.

Task (synthetic masking of TRAIN arbors, not whole-neuron or synaptic-circuit prediction): each query reveals ONLY
the first 300 um of root-path distance of its canonical TRAIN axon arbor, cut exactly at the budget by ROOT's contract
(analysis/axon_completion_contract.py, read-only, pinned). A forecast is an append-only Completion of that public
Observation; ROOT's score_completion scores the NEW cable only against the hidden suffix (50 um sampling, 200 um
tolerance). MissingTruth is created only inside a scoring job and is never an argument of any predictor.

Every method (candidate and controls) chooses ONE donor arbor from TRAIN fit animals and passes it through the SAME
decoder `decode` with the SAME cable guard:
  * the donor is cut at the same 300 um budget; each donor suffix component (one per donor boundary vertex) is used
    at most once (here: exactly once) and is grafted onto one visible query frontier;
  * frontiers and components are matched by visible position and direction (Hungarian on the first min(m, k), the
    remaining components to their cheapest frontier);
  * the component root is MERGED into the observed frontier vertex: no connector edge, no bridge to the forest root;
  * the frozen compact attachment (500 um, then 1,000 um) is used only while the component's cable does not exceed its
    original donor cable; otherwise rigid relocation (moves the component, preserves its length exactly);
  * guard: total new cable <= the donor's total hidden cable; the observation is preserved exactly (Completion).

TRAINING DATA ONLY: the 4,368 included TRAIN neurons of 101 brain samples. No validation or test row is read.

Subcommands (all output under data/interim/axon_completion_v1/ unless --out is given):
  index   public observations (prefix samples, frontiers) and donor-side suffix summaries; eligibility (no scores)
  plan    folds, calls, candidate sets and the label-pair union (no truth opened)
  run     the experiment, continuously and resumably; refuses unless protocol FINAL, code committed, contract pinned
  smoke   a small end-to-end pass into a scratch folder (never a result)
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT_DEFAULT = ROOT / "data/interim/axon_completion_v1"
PROTOCOL = ROOT / "reports/axon_completion_v1_protocol.md"
V2_PATH = ROOT / "analysis/next_donor_selector_v2.py"
CONTRACT_PATH = ROOT / "analysis/axon_completion_contract.py"

V2_SHA = "90567075fdf4c482074a2325e22c9d5a25dfdbd6f286d3e3a2c8232192403418"     # frozen selector v2 (read-only reuse)
CONTRACT_SHA = "13e6cd60b2021e0d673407ec7aaf5906a3652b41310a11a41053148f9c41c8bd"   # ROOT-finalized contract (read-only)

BUDGET_UM = 300.0
STEP_UM = 50.0
TOLERANCE_UM = 200.0
K_SOMA = 16            # the engine's E1/E2 local set (table soma)
K_PREFIX = 8           # + 8 nearest by absolute prefix cost + 8 nearest by soma-centred prefix cost
DIRECTION_UM = 100.0   # matching cost = |frontier - component root| + 100 um * (1 - cos(direction))
LADDER_UM = (500.0, 1000.0)
INFLATED_RATIO = 2.0   # a forecast is cable-inflated when new cable > 2 x hidden cable

SALT = "axon_completion_v1|20261113"
SEED = 20261113
BOOTSTRAP_SEED, BOOTSTRAP_DRAWS = 20261110, 10000
OUTER_FOLDS, INNER_FOLDS = 5, 4
HGB = dict(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=50)
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")

# predeclared gate (protocol section 7)
WIN_DISTANCE_UM = 10.0
WIN_F1 = 0.010
INFLATION_SLACK = 0.02

SELECTOR = "completion_selector"
FALLBACK_SET = ("E1_local_uniform", "E2_local_kernel", "local_prefix_nearest", "prefix_abs_nearest",
                "prefix_centred_nearest", "prefix_peptide_nearest")
METADATA_FALLBACKS = ("E1_local_uniform", "E2_local_kernel")
COPIES = ("nearest", "medoid", "peptide_nearest", "peptide_medoid")
CONTROLS = ("E1_local_uniform", "E2_local_kernel", "adaptive_fallback", "local_prefix_nearest", "prefix_abs_nearest",
            "prefix_centred_nearest", "prefix_peptide_nearest", "metadata_only", "nearest", "medoid",
            "peptide_nearest", "peptide_medoid", "switch_random", "permuted_labels")
DIAGNOSTIC = ("observation_only",)          # a failure control: must be reported as failed (empty forecasts)
METHODS = (SELECTOR,) + CONTROLS + DIAGNOSTIC


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value, *, exclusive=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x" if exclusive else "w") as f:
        f.write(json.dumps(value, sort_keys=True, indent=1, allow_nan=False) + "\n")
    return sha_file(path)


def read_json(path):
    return json.loads(Path(path).read_text())


def stable_seed(*parts):
    return int(hashlib.sha256("|".join(str(p) for p in (SEED,) + parts).encode()).hexdigest()[:8], 16)


class LeakError(RuntimeError):
    """A read outside the fit animals of a fit call, a held-out row, or truth passed to a predictor."""


def leak_types():
    """Both leak classes: this module's and the frozen v2 FitDoor's (a different class object)."""
    return (LeakError, V2().LeakError)


class CableGuardError(RuntimeError):
    """The shared decoder would emit more new cable than the donor's hidden cable."""


# ------------------------------------------------------------------------------------------------------------------
# pinned read-only modules
# ------------------------------------------------------------------------------------------------------------------

_MODS = {}


def _load(name, path, pin, *, allow_unpinned=False):
    if name not in _MODS:
        actual = sha_file(path)
        if pin is None and not allow_unpinned:
            raise RuntimeError(f"{path.name} is not pinned yet")
        if pin is not None and actual != pin:
            raise RuntimeError(f"{path.name} changed (sha {actual[:12]}, pinned {pin[:12]})")
        spec = importlib.util.spec_from_file_location(f"axon_completion_v1_{name}", path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        _MODS[name] = module
    return _MODS[name]


_DEV = {"allow_unpinned_contract": False}


def V2():
    return _load("v2", V2_PATH, V2_SHA)


def C():
    return _load("contract", CONTRACT_PATH, CONTRACT_SHA, allow_unpinned=_DEV["allow_unpinned_contract"])


def T():
    return V2().template()


def G():
    return T().G


# ------------------------------------------------------------------------------------------------------------------
# public observation features and donor suffixes
# ------------------------------------------------------------------------------------------------------------------

def unit(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return np.where(n > 0, v / np.where(n > 0, n, 1.0), 0.0)


def obs_depth(obs):
    xyz, parent = obs.xyz, obs.parent
    out = np.zeros(len(xyz))
    for i in range(1, len(xyz)):
        out[i] = out[parent[i]] + np.linalg.norm(xyz[i] - xyz[parent[i]])
    return out


def frontiers(obs):
    """Visible frontier positions and incoming directions, from the public observation only."""
    if not isinstance(obs, C().Observation):
        raise LeakError("frontiers need a public Observation")
    idx = np.asarray(obs.frontier_indices(), dtype=int)
    idx = idx[idx > 0]
    pos = obs.xyz[idx]
    direction = unit(obs.xyz[idx] - obs.xyz[obs.parent[idx]])
    return idx, pos, direction


def prefix_samples(obs):
    """50 um arc-length samples of the visible prefix (contract sampler), absolute coordinates, with masses."""
    Cm = C()
    forest = Cm.CableForest(obs.xyz, obs.parent, obs_depth(obs))
    pts, w = Cm.sample_forest(forest, STEP_UM)
    if len(pts) == 0:
        pts, w = obs.xyz[:1].copy(), np.ones(1)
    return pts, w


class DonorSuffix:
    """The suffix forest of a DONOR arbor beyond the same 300 um budget (whole donor geometry is allowed for fit
    donors). One component per donor boundary vertex: root = the exact boundary point; the visible direction into it."""

    def __init__(self, arbor, budget=BUDGET_UM):
        xyz, parent = np.asarray(arbor.xyz, float), np.asarray(arbor.parent)
        depth = arbor.path_lengths()
        n = len(xyz)
        comp = np.full(n, -1)                  # component of each hidden node
        local = np.full(n, -1)                 # its index inside the component
        roots, dirs, comps = [], [], []        # comps: lists of (xyz rows, parent rows)
        boundary_of_vertex = {}
        for i in range(1, n):
            p = int(parent[i])
            if depth[i] <= budget and comp[p] < 0:
                continue                       # visible
            if comp[p] >= 0:                   # inside an existing component
                c = comp[p]
                comp[i], local[i] = c, len(comps[c][0])
                comps[c][0].append(xyz[i]); comps[c][1].append(local[p])
                continue
            # p is visible, i is hidden: the edge crosses the budget, or p sits exactly on it
            if depth[p] >= budget:              # boundary vertex p itself (exact-budget vertex): one shared root
                if p not in boundary_of_vertex:
                    boundary_of_vertex[p] = len(comps)
                    roots.append(xyz[p].copy())
                    dirs.append(unit(xyz[p] - xyz[parent[p]]) if p > 0 else np.zeros(3))
                    comps.append(([xyz[p].copy()], [-1]))
                c = boundary_of_vertex[p]
            else:                               # exact interpolation on the crossing edge
                fraction = (budget - depth[p]) / (depth[i] - depth[p])
                cut = xyz[p] + fraction * (xyz[i] - xyz[p])
                c = len(comps)
                roots.append(cut)
                dirs.append(unit(xyz[i] - xyz[p]))
                comps.append(([cut], [-1]))
            comp[i], local[i] = c, len(comps[c][0])
            comps[c][0].append(xyz[i]); comps[c][1].append(0)
        Gm = G()
        self.components = []
        for cx, cp in comps:
            types = np.full(len(cx), 2, dtype=int)
            types[0] = 1
            self.components.append(Gm.Arbor(np.asarray(cx, float), np.asarray(cp, dtype=int), types).validate())
        self.roots = np.asarray(roots, float).reshape(-1, 3)
        self.dirs = np.asarray(dirs, float).reshape(-1, 3)
        self.cable = np.array([a.cable() for a in self.components], dtype=float)
        self.total_cable = float(self.cable.sum())

    @property
    def k(self):
        return len(self.components)


def assignment(F, V, R, U):
    """Frontier/component matching on visible position and direction. Returns (frontier of each component, cost)."""
    from scipy.optimize import linear_sum_assignment
    m, k = len(F), len(R)
    if m == 0 or k == 0:
        return np.full(k, -1), np.zeros((m, k))
    dist = np.linalg.norm(F[:, None, :] - R[None, :, :], axis=2)
    cos = np.clip(V @ U.T, -1.0, 1.0)
    cost = dist + DIRECTION_UM * (1.0 - cos)
    rows, cols = linear_sum_assignment(cost)
    target = np.full(k, -1)
    target[cols] = rows
    for j in np.flatnonzero(target < 0):
        target[j] = int(np.lexsort((np.arange(m), cost[:, j]))[0])
    return target, cost


def plan_features(F, V, R, U):
    """Decoder-plan diagnostics (no grafting, no truth): used as prefix-aware features."""
    target, cost = assignment(F, V, R, U)
    m, k = len(F), len(R)
    if k == 0 or m == 0:
        return {"plan_cost_mean": 1e4, "plan_cost_max": 1e4, "plan_shift_mean": 1e4, "plan_shift_max": 1e4,
                "plan_cos_mean": -1.0, "plan_unserved": float(m), "plan_share_500": 0.0, "plan_share_1000": 0.0,
                "plan_max_per_frontier": 0.0}
    c = cost[target, np.arange(k)]
    shift = np.linalg.norm(F[target] - R, axis=1)
    cos = np.einsum("ij,ij->i", V[target], U)
    served = np.bincount(target, minlength=m)
    return {"plan_cost_mean": float(c.mean()), "plan_cost_max": float(c.max()), "plan_shift_mean": float(shift.mean()),
            "plan_shift_max": float(shift.max()), "plan_cos_mean": float(cos.mean()),
            "plan_unserved": float((served == 0).sum()), "plan_share_500": float((shift <= LADDER_UM[0] / 3).mean()),
            "plan_share_1000": float((shift <= LADDER_UM[1] / 3).mean()), "plan_max_per_frontier": float(served.max())}


def graft(component, target_xyz):
    """Root of `component` moved onto the frontier: compact attachment (500, then 1,000 um) while its cable does not
    exceed the original component cable, else rigid relocation (exact length)."""
    Tm = T()
    original = component.cable()
    for radius in LADDER_UM:
        try:
            out, _ = Tm.compact_attachment(component, target_xyz, radius)
        except Tm.UnsupportedAttachment:
            continue
        if out.cable() <= original * (1 + 1e-12) + 1e-9:
            return out, f"attach_{int(radius)}"
    return component.shifted(target_xyz), "rigid"


def decode(obs, donor):
    """THE shared decoder. Inputs: a public Observation and a DonorSuffix (fit donor). Returns (Completion, diag)."""
    Cm = C()
    if isinstance(donor, Cm.MissingTruth) or not isinstance(donor, DonorSuffix):
        raise LeakError("the decoder takes a donor suffix, never MissingTruth")
    if not isinstance(obs, Cm.Observation):
        raise LeakError("the decoder takes a public Observation")
    fidx, F, V = frontiers(obs)
    modes = Counter()
    if len(F) == 0 or donor.k == 0:
        return Cm.observation_only(obs), {"modes": {}, "components": donor.k, "frontiers": len(F), "extra_cable_um": 0.0,
                                          "donor_hidden_cable_um": donor.total_cable, "no_graft": True}
    target, _ = assignment(F, V, donor.roots, donor.dirs)
    n = len(obs.xyz)
    ex, ep = [], []
    child_at = {}                                 # (new parent index, exact child xyz) -> generated node index
    merged_cable = 0.0
    for j, comp in enumerate(donor.components):
        f = int(target[j])
        out, mode = graft(comp, F[f])
        modes[mode] += 1
        mapping = np.empty(len(out.xyz), dtype=int)
        mapping[0] = fidx[f]                      # the component root is MERGED into the observed frontier vertex
        for t in range(1, len(out.xyz)):
            parent_new = int(mapping[out.parent[t]])
            key = (parent_new, out.xyz[t].astype("<f8").tobytes())
            if key in child_at:                   # identical generated child of the SAME parent: shared trunk
                mapping[t] = child_at[key]
                merged_cable += float(np.linalg.norm(out.xyz[t] - (obs.xyz[parent_new] if parent_new < n else ex[parent_new - n])))
                continue
            mapping[t] = n + len(ex)
            child_at[key] = mapping[t]
            ex.append(out.xyz[t])
            ep.append(parent_new)
    completion = Cm.Completion(obs, np.asarray(ex, float).reshape(-1, 3), np.asarray(ep, dtype=int))
    extra = Cm.extra_forest(completion).cable()
    if extra > donor.total_cable * (1 + 1e-9) + 1e-6:
        raise CableGuardError(f"new cable {extra:.3f} > donor hidden cable {donor.total_cable:.3f}")
    return completion, {"modes": dict(modes), "components": donor.k, "frontiers": len(F), "extra_cable_um": extra,
                        "donor_hidden_cable_um": donor.total_cable, "merged_duplicate_cable_um": merged_cable,
                        "no_graft": False}


def completion_sha(completion):
    h = hashlib.sha256()
    h.update(completion.observation.sha256.encode())
    for a, kind in ((completion.extra_xyz, "<f8"), (completion.extra_parent, "<i8")):
        data = np.ascontiguousarray(np.asarray(a).astype(kind, copy=False))
        h.update(str(data.shape).encode())
        h.update(data.tobytes())
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------------------------
# the index: public observation summaries (query side) and donor suffix summaries (donor side), per TRAIN neuron
# ------------------------------------------------------------------------------------------------------------------

_W = {}


def _worker_init(allow_unpinned=False):
    for v in THREAD_VARS:
        os.environ.setdefault(v, "1")
    _DEV["allow_unpinned_contract"] = allow_unpinned
    _W["lib"] = V2().Library()
    C(); T()


def index_job(nid):
    lib = _W["lib"]
    arbor = lib.load(nid)
    obs = C().observe(arbor, BUDGET_UM)
    fidx, F, V = frontiers(obs)
    pts, w = prefix_samples(obs)
    donor = DonorSuffix(arbor)
    s = V2().summarise(arbor)
    return {"id": nid, "obs_sha": obs.sha256, "root": obs.xyz[0].tolist(), "F": F.tolist(), "V": V.tolist(),
            "pts": pts.tolist(), "w": w.tolist(), "obs_cable": float(C().CableForest(obs.xyz, obs.parent, obs_depth(obs)).cable()),
            "obs_vertices": int(len(obs.xyz)),
            "donor_R": donor.roots.tolist(), "donor_U": donor.dirs.tolist(), "donor_cable": donor.cable.tolist(),
            "donor_hidden_cable": donor.total_cable, "summary": s,
            "registration_um": float(np.linalg.norm(lib.soma[lib.index[nid]] - obs.xyz[0]))}


class Index:
    """Per-neuron public observation summaries + donor-side suffix summaries. `query_*` fields derive from the public
    observation only; `donor_*` fields are read only for donors of a fit pool."""

    def __init__(self, lib, rows):
        if [r["id"] for r in rows] != lib.ids:
            raise RuntimeError("index order differs from the library")
        self.lib = lib
        self.n = lib.n
        self.obs_sha = [r["obs_sha"] for r in rows]
        self.root = np.array([r["root"] for r in rows])
        self.F = [np.asarray(r["F"], float).reshape(-1, 3) for r in rows]
        self.V = [np.asarray(r["V"], float).reshape(-1, 3) for r in rows]
        self.pts = [np.asarray(r["pts"], float).reshape(-1, 3) for r in rows]
        self.w = [np.asarray(r["w"], float) for r in rows]
        self.obs_cable = np.array([r["obs_cable"] for r in rows])
        self.obs_vertices = np.array([r["obs_vertices"] for r in rows], float)
        self.donor_R = [np.asarray(r["donor_R"], float).reshape(-1, 3) for r in rows]
        self.donor_U = [np.asarray(r["donor_U"], float).reshape(-1, 3) for r in rows]
        self.donor_k = np.array([len(r["donor_cable"]) for r in rows], float)
        self.donor_hidden_cable = np.array([r["donor_hidden_cable"] for r in rows])
        keys = list(rows[0]["summary"])
        self.summary = {k: np.array([r["summary"][k] for r in rows], float) for k in keys}
        self.registration = np.array([r["registration_um"] for r in rows])
        # fixed eligibility rule (contract MissingTruth.eligible: hidden cable > 1e-9), decided before any score
        self.eligible = self.donor_hidden_cable > 1e-9
        self._abs = self._centred = None

    def prefix_costs(self):
        """Symmetric arc-length-weighted chamfer between visible prefixes, absolute and soma(root)-centred.
        Held in memory only (never serialised)."""
        if self._abs is None:
            counts = np.array([len(p) for p in self.pts])
            starts = np.r_[0, np.cumsum(counts)[:-1]]
            owner = np.repeat(np.arange(self.n), counts)
            P = np.vstack(self.pts)
            W = np.concatenate(self.w)
            Pc = P - self.root[owner]
            wsum = np.bincount(owner, W, minlength=self.n)
            out = []
            for allp in (P, Pc):
                M = np.empty((self.n, self.n), dtype=np.float32)
                for q in range(self.n):
                    qp = self.pts[q] - (self.root[q] if allp is Pc else 0.0)
                    qw = self.w[q] / self.w[q].sum()
                    centre = qp[0]                     # centred at the query's first sample: exact where it matters
                    a_, b_ = qp - centre, allp - centre
                    D2 = (a_ ** 2).sum(axis=1)[:, None] + (b_ ** 2).sum(axis=1)[None, :] - 2.0 * (a_ @ b_.T)
                    D = np.sqrt(np.maximum(D2, 0.0))
                    to_d = np.minimum.reduceat(D, starts, axis=1)          # (m, n): query point -> nearest donor point
                    a = qw @ to_d
                    b = np.bincount(owner, W * D.min(axis=0), minlength=self.n) / wsum
                    M[q] = 0.5 * (a + b)
                out.append(M)
            self._abs, self._centred = out
        return self._abs, self._centred


def build_index(lib, pool, log, out):
    path = out / "index.jsonl"
    done = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                done[r["id"]] = r
    todo = [n for n in lib.ids if n not in done]
    log(f"index: {len(done)} done, {len(todo)} to build")
    with open(path, "a") as f:
        for i, r in enumerate(pool.imap_unordered(index_job, todo, chunksize=8)):
            f.write(json.dumps(r, separators=(",", ":")) + "\n")
            done[r["id"]] = r
            if (i + 1) % 500 == 0:
                f.flush()
                log(f"index: {i + 1}/{len(todo)}")
    return Index(lib, [done[n] for n in lib.ids])


# ------------------------------------------------------------------------------------------------------------------
# folds, pools, candidate sets and control rules
# ------------------------------------------------------------------------------------------------------------------

def folds(animals, salt, k):
    order = sorted(set(animals), key=lambda a: (hashlib.sha256(f"{SALT}|{salt}|{a}".encode()).hexdigest(), a))
    return [sorted(order[j::k]) for j in range(k)]


def fold_plan(animals):
    animals = sorted(set(animals))
    plan = []
    for k, ev in enumerate(folds(animals, "outer", OUTER_FOLDS)):
        fit = sorted(set(animals) - set(ev))
        inner = folds(fit, f"inner|{k}", INNER_FOLDS)
        plan.append({"outer": k, "eval_animals": ev, "fit_animals": fit,
                     "inner": [{"inner": j, "eval_animals": ie, "fit_animals": sorted(set(fit) - set(ie))}
                               for j, ie in enumerate(inner)]})
    return plan


def neurons_of(lib, idx, animals):
    return [int(i) for i in np.flatnonzero(np.isin(lib.animal, list(animals)) & idx.eligible)]


def _argmin(values, members):
    """Smallest value; ties by library position (= ID order)."""
    order = np.lexsort((members, values))
    return int(members[order[0]])


class Local:
    """One query in one pool: the candidate union and every control's donor (all by known inputs/observation only)."""

    def __init__(self, lib, idx, door, q, pool_mask):
        V2m = V2()
        if pool_mask[q] or (pool_mask & (lib.animal == lib.animal[q])).any():
            raise LeakError("the query's animal is in its own pool")
        pool_mask = pool_mask & idx.eligible
        loc = V2m.local_choices(lib, door, q, pool_mask)          # soma16, E1/E2 weights, risks, choices (via door)
        self.q = q
        self.soma16 = loc["cand"]
        self.soma_dist16 = loc["dist"]
        self.weights, self.risk, self.h = loc["weights"], loc["risk"], loc["h"]
        A, Cc = idx.prefix_costs()
        pool = np.flatnonzero(pool_mask)
        a, c = A[q, pool].astype(float), Cc[q, pool].astype(float)
        abs8 = pool[np.lexsort((pool, a))[:K_PREFIX]]
        cen8 = pool[np.lexsort((pool, c))[:K_PREFIX]]
        pep = lib.peptide[q]
        same = pool[lib.peptide[pool] == pep] if pep else pool[:0]
        pep_choice = _argmin(A[q, same].astype(float), same) if len(same) else _argmin(a, pool)
        self.choice = {
            "E1_local_uniform": int(self.soma16[loc["choice"]["E1_local_uniform"]]),
            "E2_local_kernel": int(self.soma16[loc["choice"]["E2_local_kernel"]]),
            "local_prefix_nearest": _argmin(A[q, self.soma16].astype(float), self.soma16),
            "prefix_abs_nearest": _argmin(a, pool),
            "prefix_centred_nearest": _argmin(c, pool),
            "prefix_peptide_nearest": pep_choice,
        }
        self.cand = np.array(sorted(set(self.soma16.tolist()) | set(abs8.tolist()) | set(cen8.tolist())
                                    | set(self.choice.values())), dtype=int)
        self.in_soma16 = np.isin(self.cand, self.soma16)
        self.in_abs8 = np.isin(self.cand, abs8)
        self.in_cen8 = np.isin(self.cand, cen8)
        self.abs_cost = A[q, self.cand].astype(float)
        self.cen_cost = Cc[q, self.cand].astype(float)
        self.abs_rank = np.searchsorted(np.sort(a), self.abs_cost, side="left").astype(float)
        self.cen_rank = np.searchsorted(np.sort(c), self.cen_cost, side="left").astype(float)


def pool_mask(lib, animals, q):
    return np.isin(lib.animal, list(animals)) & (lib.animal != lib.animal[q])


# ------------------------------------------------------------------------------------------------------------------
# features
# ------------------------------------------------------------------------------------------------------------------

class Featurizer:
    """Fitted on the fit animals of one call. prefix=False gives the metadata-only control (no observation input)."""

    def __init__(self, lib, idx, door, *, prefix):
        self.lib, self.idx, self.prefix = lib, idx, prefix
        fit = np.flatnonzero(door.mask)
        self.peptides = sorted(set(lib.peptide[fit].tolist()) - {""})
        self.nuclei = sorted(k for k, c in Counter(lib.nucleus[fit].tolist()).items() if k and c >= 20)
        self.hemis = sorted(set(lib.hemisphere[fit].tolist()) - {""})
        D = door.D(fit, fit)
        count = Counter(lib.animal[fit].tolist())
        w = np.array([1.0 / count[a] for a in lib.animal[fit]])
        self.centrality = np.full(lib.n, np.nan)
        self.centrality[fit] = D @ (w / w.sum())
        self.names = None

    def query_obs(self, q):
        """Observation-derived query features (public prefix only)."""
        idx = self.idx
        F, V = idx.F[q], idx.V[q]
        rel = F - idx.root[q] if len(F) else np.zeros((1, 3))
        return {"q_frontiers": float(len(F)), "q_obs_cable": idx.obs_cable[q], "q_obs_vertices": idx.obs_vertices[q],
                **{f"q_front_off_{a}": float(rel[:, i].mean()) for i, a in enumerate("xyz")},
                **{f"q_front_dir_{a}": float(V[:, i].mean()) if len(V) else 0.0 for i, a in enumerate("xyz")},
                "q_front_spread": float(np.linalg.norm(rel - rel.mean(axis=0), axis=1).mean())}

    def rows(self, local, fallback_donor, *, candidates=None):
        lib, idx = self.lib, self.idx
        q = local.q
        cand = local.cand if candidates is None else candidates
        sel = np.searchsorted(local.cand, cand)
        if not np.array_equal(local.cand[sel], cand):
            raise RuntimeError("candidates outside the local union")
        if not np.isfinite(self.centrality[cand]).all() or not np.isfinite(self.centrality[fallback_donor]):
            raise LeakError("a candidate donor outside the fit animals")
        known = lib.known(q)
        qs = np.asarray(known["soma"], float)
        n = len(cand)
        cols = {}
        for a, name in enumerate(("ap", "dv", "ml")):
            cols[f"q_soma_{name}"] = np.full(n, qs[a])
            cols[f"diff_{name}"] = lib.soma[cand, a] - qs[a]
        sd = np.linalg.norm(lib.soma[cand] - qs, axis=1)
        cols["soma_distance"] = sd
        cols["soma_distance_vs_fallback"] = sd - np.linalg.norm(lib.soma[fallback_donor] - qs)
        for p in self.peptides:
            cols[f"q_pep_{p}"] = np.full(n, float(known["peptide"] == p))
        for nuc in self.nuclei:
            cols[f"q_nuc_{nuc}"] = np.full(n, float(known["nucleus"] == nuc))
        for hm in self.hemis:
            cols[f"q_hemi_{hm}"] = np.full(n, float(known["hemisphere"] == hm))
        cols["same_peptide"] = (lib.peptide[cand] == known["peptide"]).astype(float)
        cols["same_nucleus"] = (lib.nucleus[cand] == known["nucleus"]).astype(float)
        cols["same_hemisphere"] = (lib.hemisphere[cand] == known["hemisphere"]).astype(float)
        for key, values in idx.summary.items():
            cols[f"donor_{key}"] = values[cand]
        cols["donor_suffix_log_cable"] = np.log1p(idx.donor_hidden_cable[cand])
        cols["donor_suffix_log_cable_vs_fallback"] = cols["donor_suffix_log_cable"] - math.log1p(idx.donor_hidden_cable[fallback_donor])
        cols["donor_components"] = idx.donor_k[cand]
        cols["donor_centrality"] = self.centrality[cand]
        cols["donor_centrality_vs_fallback"] = self.centrality[cand] - self.centrality[fallback_donor]
        pos16 = {int(d): i for i, d in enumerate(local.soma16)}
        e1 = local.choice["E1_local_uniform"]
        e2 = local.choice["E2_local_kernel"]
        for name, key in (("w_E1", "E1_local_uniform"), ("w_E2", "E2_local_kernel")):
            cols[name] = np.array([local.weights[key][pos16[int(d)]] if int(d) in pos16 else 0.0 for d in cand])
        for name, key in (("risk_E1_rel", "E1_local_uniform"), ("risk_E2_rel", "E2_local_kernel")):
            r = local.risk[key]
            cols[name] = np.array([r[pos16[int(d)]] - r.min() if int(d) in pos16 else -1.0 for d in cand])
        cols["in_soma16"] = local.in_soma16[sel].astype(float)
        cols["is_E1_choice"] = (cand == e1).astype(float)
        cols["is_E2_choice"] = (cand == e2).astype(float)
        cols["is_fallback"] = (cand == fallback_donor).astype(float)
        if self.prefix:
            for k, v in self.query_obs(q).items():
                cols[k] = np.full(n, v)
            cols["registration_offset"] = np.full(n, idx.registration[q])
            cols["in_abs8"] = local.in_abs8[sel].astype(float)
            cols["in_centred8"] = local.in_cen8[sel].astype(float)
            cols["prefix_abs_cost"] = local.abs_cost[sel]
            cols["prefix_centred_cost"] = local.cen_cost[sel]
            cols["prefix_abs_rank"] = local.abs_rank[sel]
            cols["prefix_centred_rank"] = local.cen_rank[sel]
            A, Cc = idx.prefix_costs()
            cols["prefix_abs_cost_vs_fallback"] = local.abs_cost[sel] - float(A[q, fallback_donor])
            cols["prefix_centred_cost_vs_fallback"] = local.cen_cost[sel] - float(Cc[q, fallback_donor])
            plans = [plan_features(idx.F[q], idx.V[q], idx.donor_R[d], idx.donor_U[d]) for d in cand]
            fb_plan = plan_features(idx.F[q], idx.V[q], idx.donor_R[fallback_donor], idx.donor_U[fallback_donor])
            for key in plans[0]:
                cols[key] = np.array([p[key] for p in plans])
            cols["plan_cost_mean_vs_fallback"] = cols["plan_cost_mean"] - fb_plan["plan_cost_mean"]
            cols["components_minus_frontiers"] = idx.donor_k[cand] - len(idx.F[q])
        names = list(cols)
        if self.names is None:
            V2().check_feature_names(names)
            self.names = names
        elif names != self.names:
            raise RuntimeError("feature columns changed within a fit call")
        X = np.column_stack([cols[k] for k in names]).astype(float)
        if not np.isfinite(X).all():
            raise RuntimeError("non-finite feature")
        return X


OBSERVATION_FEATURE_PREFIXES = ("q_frontiers", "q_obs", "q_front", "registration", "in_abs8", "in_centred8", "prefix_",
                                "plan_", "components_minus")


# ------------------------------------------------------------------------------------------------------------------
# labels: exact contract scores of the shared decoder's completion (truth opened only inside the scoring job)
# ------------------------------------------------------------------------------------------------------------------

def label_job(job):
    """One query: split its TRAIN arbor (MissingTruth stays inside this function), decode each donor from the PUBLIC
    observation only, score with the contract."""
    lib = _W["lib"]
    Cm = C()
    qid = job["query"]
    q = lib.index[qid]
    obs, truth = Cm.split_arbor(lib.load(qid), BUDGET_UM)
    if obs.sha256 != job["obs_sha"]:
        raise RuntimeError(f"{qid}: observation differs from the index")
    out = []
    for did in job["donors"]:
        if did is None:
            completion, diag = Cm.observation_only(obs), {"modes": {}, "no_graft": True}
        else:
            if lib.animal[lib.index[did]] == lib.animal[q]:
                raise LeakError("same-animal pair in a label job")
            try:
                completion, diag = decode(obs, DonorSuffix(lib.load(did)))
            except CableGuardError as error:
                out.append({"donor": did, "status": "failed_cable_guard", "valid": False, "error": str(error)})
                continue
        m = Cm.metrics_to_json(Cm.score_completion(completion, truth, step_um=STEP_UM, tolerance_um=TOLERANCE_UM))
        out.append({"donor": did, "status": m["status"], "valid": bool(m["valid_for_comparison"]),
                    "eligible": bool(m["eligible"]), "distance_um": m["distance_um"], "f1_200": m["f1_200"],
                    "precision_200": m["precision_200"], "recall_200": m["recall_200"],
                    "extra_cable_um": m["extra_cable_um"], "hidden_cable_um": m["hidden_cable_um"],
                    "empty": bool(m["empty_forecast"]), "duplicate_um": m["duplicate_extra_cable_um"],
                    "modes": diag.get("modes", {}), "merged_duplicate_cable_um": diag.get("merged_duplicate_cable_um", 0.0),
                    "cable_relative_error": m.get("cable_relative_error"),
                    "completion_sha256": completion_sha(completion)})
    return {"job": job["job"], "query": qid, "pairs": out}


def decode_job(job):
    """Forecast only (no truth): decode each chosen donor from the public observation; return completion hashes."""
    lib = _W["lib"]
    Cm = C()
    qid = job["query"]
    obs = Cm.Observation.from_dict(job["observation"])
    if obs.sha256 != job["obs_sha"]:
        raise RuntimeError("observation hash")
    out = {}
    for did in job["donors"]:
        if did is None:
            out["none"] = {"completion_sha256": completion_sha(Cm.observation_only(obs)), "modes": {}, "extra_cable_um": 0.0}
            continue
        try:
            completion, diag = decode(obs, DonorSuffix(lib.load(did)))
        except CableGuardError as error:
            out[did] = {"completion_sha256": None, "error": str(error)}
            continue
        out[did] = {"completion_sha256": completion_sha(completion), "modes": diag["modes"],
                    "extra_cable_um": diag["extra_cable_um"]}
    return {"job": job["job"], "query": qid, "decoded": out}


def observation_job(nid):
    lib = _W["lib"]
    return nid, C().observe(lib.load(nid), BUDGET_UM).to_dict()


class Labels:
    def __init__(self, lib):
        self.lib, self.data = lib, {}

    def add(self, qid, r):
        key = (self.lib.index[qid], None if r["donor"] is None else self.lib.index[r["donor"]])
        if key in self.data and (self.data[key].get("distance_um"), self.data[key].get("f1_200")) != (r.get("distance_um"), r.get("f1_200")):
            raise RuntimeError(f"label {key} recomputed differently")
        self.data[key] = r

    def get(self, q, d):
        return self.data[(q, d)]

    def has(self, q, d):
        return (q, d) in self.data

    def load(self, path):
        done = set()
        if Path(path).exists():
            for line in Path(path).read_text().splitlines():
                if line.strip():
                    r = json.loads(line)
                    for p in r["pairs"]:
                        self.add(r["query"], p)
                    done.add(r["job"])
        return done


def run_jobs(fn, jobs, path, pool, log, what, store=None):
    done = set()
    if Path(path).exists():
        for line in Path(path).read_text().splitlines():
            if line.strip():
                done.add(json.loads(line)["job"])
    todo = [j for j in jobs if j["job"] not in done]
    total = sum(len(j["donors"]) for j in jobs)
    log(f"{what}: {len(jobs)} jobs, {total} pairs; {len(jobs) - len(todo)} jobs already done; {len(todo)} to run")
    start, last, n_pairs = time.time(), time.time(), 0
    with open(path, "a") as f:
        for i, result in enumerate(pool.imap_unordered(fn, todo, chunksize=2)):
            f.write(json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n")
            f.flush()
            n_pairs += len(result.get("pairs", result.get("decoded", [])))
            now = time.time()
            if now - last > 60 or i + 1 == len(todo):
                rate = n_pairs / max(now - start, 1e-9)
                left = total - (total - sum(len(j["donors"]) for j in todo)) - n_pairs
                log(f"{what}: {i + 1}/{len(todo)} jobs, {rate:.1f} pairs/s, ETA {left / max(rate, 1e-9) / 60:.1f} min")
                last = now
    if store is not None:
        store.load(path)


# ------------------------------------------------------------------------------------------------------------------
# the learner: one fixed HGB configuration, two heads, fitted to residual gains over the inner-chosen fallback
# ------------------------------------------------------------------------------------------------------------------

def make_hgb():
    from sklearn.ensemble import HistGradientBoostingRegressor
    return HistGradientBoostingRegressor(random_state=SEED, **HGB)


def joint_loss(distance, f1):
    if distance is None or f1 is None or not math.isfinite(distance):
        return math.inf
    return math.log1p(distance / 200.0) + (1.0 - f1)


def logd(distance):
    return math.log1p(distance / 200.0)


class ResidualHeads:
    """gain_d = log1p(d_fallback/200) - log1p(d/200); gain_f = F1 - F1_fallback (positive = better than fallback)."""

    def fit(self, X, gd, gf, w):
        self.d, self.f = make_hgb(), make_hgb()
        self.d.fit(X, gd, sample_weight=w)
        self.f.fit(X, gf, sample_weight=w)
        return self

    def predict(self, X):
        return self.d.predict(X), self.f.predict(X)


def residual_choice(cand, fallback, gd, gf):
    """Switch to the candidate with the largest predicted gd + gf among those predicted better than the fallback on
    BOTH heads (gd > 0 and gf > 0); ties by ID position; otherwise the fallback."""
    gd, gf = np.asarray(gd, float), np.asarray(gf, float)
    if not (np.isfinite(gd).all() and np.isfinite(gf).all()):
        return int(fallback), "non_finite_prediction"
    ok = (gd > 0) & (gf > 0) & (cand != fallback)
    if not ok.any():
        return int(fallback), "fallback"
    score = np.where(ok, gd + gf, -np.inf)
    top = int(np.lexsort((np.arange(len(cand)), -score))[0])
    return int(cand[top]), "switched"


def training_rows(lib, door, feat, locals_, fallback_of, labels, *, metadata, permute_seed=None):
    X, gd, gf, w, skipped = [], [], [], [], 0
    per_animal = Counter(lib.animal[q] for q in locals_)
    for q in sorted(locals_):
        loc = locals_[q]
        fb = fallback_of[q]
        cand = loc.soma16 if metadata else loc.cand
        cand = np.sort(cand)
        lf = door.label(q, fb)
        if not lf["valid"]:
            skipped += len(cand)
            continue
        rows = feat.rows(loc, fb, candidates=cand)
        keep = []
        for i, d in enumerate(cand):
            ld = door.label(q, int(d))
            if not ld["valid"]:
                skipped += 1
                continue
            keep.append(i)
            gd.append(logd(lf["distance_um"]) - logd(ld["distance_um"]))
            gf.append(ld["f1_200"] - lf["f1_200"])
        X.append(rows[keep])
        w.extend([1.0 / (per_animal[lib.animal[q]] * max(len(keep), 1))] * len(keep))
    X = np.vstack(X)
    gd, gf, w = np.array(gd), np.array(gf), np.array(w)
    w = w / w.mean()
    if permute_seed is not None:
        perm = np.random.default_rng(permute_seed).permutation(len(gd))
        gd, gf = gd[perm], gf[perm]
    return X, gd, gf, w, skipped


class LabelDoor:
    """FitDoor for labels (reuses the frozen v2 FitDoor: both endpoints in the call's fit animals, no same animal)."""

    def __init__(self, name, lib, fit_animals, labels, cache):
        self.door = V2().FitDoor(name, lib, fit_animals, labels=labels, cache=cache)
        self.mask = self.door.mask

    def label(self, q, d):
        try:
            return self.door.label(q, d)
        except V2().LeakError as error:
            raise LeakError(str(error)) from error

    def D(self, rows, cols):
        try:
            return self.door.D(rows, cols)
        except V2().LeakError as error:
            raise LeakError(str(error)) from error

    @property
    def refused(self):
        return self.door.refused


# ------------------------------------------------------------------------------------------------------------------
# aggregation, intervals, gate
# ------------------------------------------------------------------------------------------------------------------

def sample_means(lib, rows):
    """rows: q -> label dict. Neuron values within a brain sample, then equal samples. Failures stay in the
    denominator: a missing/non-finite distance counts as +inf and a missing F1 as 0 (never dropped)."""
    by = defaultdict(list)
    for q, r in rows.items():
        if not r.get("valid", False):              # invalid FIRST: finite diagnostics never enter the aggregate
            by[lib.animal[q]].append((math.inf, 0.0))
            continue
        dist, f1 = r.get("distance_um"), r.get("f1_200")
        by[lib.animal[q]].append((math.inf if dist is None else float(dist), 0.0 if f1 is None else float(f1)))
    return {a: (float(np.mean([v[0] for v in vals])), float(np.mean([v[1] for v in vals]))) for a, vals in by.items()}


def finite_json(value):
    """Strict JSON: non-finite floats become null with an explicit marker (never silently)."""
    if isinstance(value, dict):
        return {k: finite_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite_json(v) for v in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        return {"nonfinite": str(float(value))}
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def paired_bootstrap(per_a, per_b, animals):
    a = np.array([per_a[x] for x in animals], float)
    b = np.array([per_b[x] for x in animals], float)
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        return {"distance_delta_um": None, "f1_delta": float((a[:, 1] - b[:, 1]).mean()), "distance_bootstrap95_um": None,
                "f1_bootstrap95": None, "samples_both_better": None, "nonfinite": True}
    diff = a - b
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    idx = rng.integers(0, len(animals), size=(BOOTSTRAP_DRAWS, len(animals)))
    boot = diff[idx].mean(axis=1)
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    return {"distance_delta_um": float(diff[:, 0].mean()), "f1_delta": float(diff[:, 1].mean()),
            "distance_bootstrap95_um": [float(lo[0]), float(hi[0])], "f1_bootstrap95": [float(lo[1]), float(hi[1])],
            "samples_both_better": int(((diff[:, 0] < 0) & (diff[:, 1] > 0)).sum()), "nonfinite": False}


def verdict(summary, versus, validity, switch_rate):
    """The predeclared gate (protocol section 7). Every criterion must hold; otherwise FAIL with the reasons."""
    reasons = []
    sel_v = validity[SELECTOR]
    if sel_v["invalid"] > 0:
        reasons.append(f"selector has {sel_v['invalid']} invalid/empty/duplicate/guard-failed forecasts")
    if switch_rate <= 0:
        reasons.append("selector never switched from its fallback (not a learned result)")
    margin = {}
    for c in CONTROLS:
        v = versus[c]
        # an unevaluable required control (any invalid forecast or a non-finite comparison) is never "beaten":
        # the joint result is FAIL (inconclusive); its failures stay in its denominators and are reported
        if validity[c]["invalid"] > 0 or v.get("nonfinite") or v["distance_delta_um"] is None:
            margin[c] = {"beaten": False, "unevaluable": True, "invalid": validity[c]["invalid"],
                         "f1_delta": v.get("f1_delta")}
            reasons.append(f"required control {c} is unevaluable ({validity[c]['invalid']} invalid forecasts): "
                           f"inconclusive, not a win")
        else:
            ok = (v["distance_delta_um"] <= -WIN_DISTANCE_UM and v["f1_delta"] >= WIN_F1
                  and v["distance_bootstrap95_um"][1] < 0 and v["f1_bootstrap95"][0] > 0)
            margin[c] = {"beaten": bool(ok), "unevaluable": False, "distance_delta_um": v["distance_delta_um"],
                         "f1_delta": v["f1_delta"]}
            if not ok:
                reasons.append(f"not a resolved practical win over {c} (needs <= -{WIN_DISTANCE_UM} um and >= +{WIN_F1} F1, "
                               f"both 95% intervals excluding zero)")
        if sel_v["inflated_share"] > validity[c]["inflated_share"] + INFLATION_SLACK:
            reasons.append(f"cable guard: inflated share {sel_v['inflated_share']:.3f} exceeds {c} "
                           f"({validity[c]['inflated_share']:.3f}) + {INFLATION_SLACK}")
    fb_err = validity["adaptive_fallback"]["mean_abs_rel_cable_error"]
    if not sel_v["mean_abs_rel_cable_error"] <= fb_err:
        reasons.append(f"cable guard: sample-equal mean absolute relative cable error {sel_v['mean_abs_rel_cable_error']:.4f} "
                       f"exceeds the adaptive fallback's {fb_err:.4f}")
    passed = not reasons
    return {"pass": passed, "sentence": "PASS: development completion gate met" if passed else "FAIL: development completion gate not met",
            "reasons": reasons, "margins": margin, "switch_rate": switch_rate,
            "scope": "synthetic-masking DEVELOPMENT completion on TRAIN only (101 samples already used in design); "
                     "not whole-neuron or synaptic-circuit prediction; no viewer promotion"}


# ------------------------------------------------------------------------------------------------------------------
# the experiment
# ------------------------------------------------------------------------------------------------------------------

class Log:
    def __init__(self, path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def __call__(self, message):
        line = f"[{utc()}] {message}"
        print(line, flush=True)
        with self.path.open("a") as f:
            f.write(line + "\n")


def code_pins():
    files = ["analysis/axon_completion_v1.py", "reports/axon_completion_v1_protocol.md", "tests/test_axon_completion_v1.py",
             "analysis/axon_completion_contract.py", "tests/test_axon_completion_contract.py", "analysis/axon_completion_preflight.py",
             "analysis/next_donor_selector_v2.py", "analysis/brm_template.py",
             "analysis/d110_geometry.py", "analysis/brm_geometry.py"]
    return {f: sha_file(ROOT / f) for f in files if (ROOT / f).exists()}


def git_clean(paths):
    import subprocess
    status = subprocess.run(["git", "status", "--porcelain", "--", *paths], cwd=ROOT, capture_output=True, text=True).stdout
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    return status.strip() == "", head


def experiment(out, workers, log, *, animals=None, dev=False):
    import multiprocessing as mp
    lib = V2().Library()
    cache = V2().Cache(lib)
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers, initializer=_worker_init, initargs=(dev,)) as pool:
        idx = build_index(lib, pool, log, out)
        all_animals = sorted(set(lib.animal.tolist()))
        if animals is not None:
            all_animals = all_animals[:animals]
        eligible_all = int(idx.eligible.sum())
        cohort = [q for q in range(lib.n) if lib.animal[q] in set(all_animals)]
        excluded = [lib.ids[q] for q in cohort if not idx.eligible[q]]
        write_json(out / "cohort.json", {"neurons": len(cohort), "eligible": len(cohort) - len(excluded),
                                         "excluded_fully_observed": excluded, "eligible_library": eligible_all,
                                         "rule": "exclude when contract hidden cable <= 1e-9 um (fully observed); fixed before scoring",
                                         "registration_offset_um_percentiles_50_90_99_100":
                                             np.percentile(idx.registration[cohort], [50, 90, 99, 100]).tolist(),
                                         "frontiers_histogram": np.bincount([len(idx.F[q]) for q in cohort]).tolist()})
        log(f"cohort: {len(cohort)} neurons, {len(excluded)} excluded as fully observed")
        t0 = time.time()
        idx.prefix_costs()
        log(f"prefix costs in memory ({time.time() - t0:.0f} s)")
        plan = fold_plan(all_animals)
        write_json(out / "folds.json", plan)

        # ---- build every call's local sets (known inputs + public observations only) ----
        calls = {}
        for fold in plan:
            k = fold["outer"]
            door = V2().FitDoor(f"outer{k}", lib, fold["fit_animals"], cache=cache)
            for q in neurons_of(lib, idx, fold["fit_animals"]):
                calls[("fit", k, q)] = Local(lib, idx, door, q, pool_mask(lib, fold["fit_animals"], q))
            for q in neurons_of(lib, idx, fold["eval_animals"]):
                calls[("eval", k, q)] = Local(lib, idx, door, q, pool_mask(lib, fold["fit_animals"], q))
            for inner in fold["inner"]:
                idoor = V2().FitDoor(f"outer{k}_inner{inner['inner']}", lib, inner["fit_animals"], cache=cache)
                for q in neurons_of(lib, idx, inner["eval_animals"]):
                    calls[("inner", k, q)] = Local(lib, idx, idoor, q, pool_mask(lib, inner["fit_animals"], q))
        log(f"local sets: {len(calls)} (query, call) entries")

        # ---- fit label pairs: outer-fit candidate unions + inner-eval fallback-set donors ----
        pairs = defaultdict(set)
        for (kind, k, q), loc in calls.items():
            if kind == "fit":
                pairs[q].update(int(d) for d in loc.cand)
            elif kind == "inner":
                pairs[q].update(loc.choice.values())
        n_pairs = sum(len(v) for v in pairs.values())
        write_json(out / "plan_counts.json", {"fit_label_pairs": n_pairs, "queries": len(pairs),
                                              "mean_candidates": float(np.mean([len(l.cand) for l in calls.values()]))})
        log(f"fit label pairs: {n_pairs}")
        jobs = [{"job": f"fit|{lib.ids[q]}", "query": lib.ids[q], "obs_sha": idx.obs_sha[q],
                 "donors": [lib.ids[d] for d in sorted(ds)]} for q, ds in sorted(pairs.items())]
        labels = Labels(lib)
        run_jobs(label_job, jobs, out / "labels_fit.jsonl", pool, log, "fit labels", store=labels)

        # ---- per outer fold: inner fallback choice, residual fits, outer-eval forecasts ----
        forecasts, fold_info = {}, {}
        for fold in plan:
            k = fold["outer"]
            door = LabelDoor(f"outer{k}_labels", lib, fold["fit_animals"], labels, cache)
            inner_locs = {q: l for (kind, kk, q), l in calls.items() if kind == "inner" and kk == k}
            losses = {}
            for m in FALLBACK_SET:
                rows = {q: door.label(q, l.choice[m]) for q, l in inner_locs.items()}
                per = defaultdict(list)
                for q, r in rows.items():
                    per[lib.animal[q]].append(joint_loss(r["distance_um"], r["f1_200"]) if r["valid"] else math.inf)
                losses[m] = float(np.mean([np.mean(v) for v in per.values()]))
            fallback = min(FALLBACK_SET, key=lambda m: (losses[m], FALLBACK_SET.index(m)))
            meta_fallback = min(METADATA_FALLBACKS, key=lambda m: (losses[m], METADATA_FALLBACKS.index(m)))
            try:                                                   # the door must refuse an outer-eval label read
                ev0 = neurons_of(lib, idx, fold["eval_animals"])[0]
                door.label(ev0, calls[("fit", k, neurons_of(lib, idx, fold["fit_animals"])[0])].cand[0])
            except leak_types():
                pass
            else:
                raise RuntimeError("probe read of an outer-eval label was served")
            fit_locs = {q: l for (kind, kk, q), l in calls.items() if kind == "fit" and kk == k}
            feat = Featurizer(lib, idx, door, prefix=True)
            mfeat = Featurizer(lib, idx, door, prefix=False)
            fb_of = {q: l.choice[fallback] for q, l in fit_locs.items()}
            mfb_of = {q: l.choice[meta_fallback] for q, l in fit_locs.items()}
            X, gd, gf, w, skipped = training_rows(lib, door, feat, fit_locs, fb_of, labels, metadata=False)
            model = ResidualHeads().fit(X, gd, gf, w)
            Xm, gdm, gfm, wm, skm = training_rows(lib, door, mfeat, fit_locs, mfb_of, labels, metadata=True)
            mmodel = ResidualHeads().fit(Xm, gdm, gfm, wm)
            perm = np.random.default_rng(stable_seed("permuted", k)).permutation(len(gd))   # same rows, permuted targets
            pmodel = ResidualHeads().fit(X, gd[perm], gf[perm], w)
            fold_info[k] = {"fallback": fallback, "metadata_fallback": meta_fallback, "inner_losses": losses,
                            "train_rows": int(len(gd)), "skipped_invalid_rows": int(skipped),
                            "metadata_train_rows": int(len(gdm)), "metadata_skipped": int(skm),
                            "door_refused": dict(door.refused), "features": feat.names, "metadata_features": mfeat.names}
            log(f"outer {k}: fallback {fallback}, metadata fallback {meta_fallback}, rows {len(gd)}")
            copy_door = V2().FitDoor(f"outer{k}_copies", lib, fold["fit_animals"], cache=cache)
            copies = V2().Copies(lib, copy_door, np.isin(lib.animal, fold["fit_animals"]) & idx.eligible)
            for q in neurons_of(lib, idx, fold["eval_animals"]):
                loc = calls[("eval", k, q)]
                fb = loc.choice[fallback]
                gdp_, gfp_ = model.predict(feat.rows(loc, fb))
                choice, why = residual_choice(loc.cand, fb, gdp_, gfp_)
                mfb = loc.choice[meta_fallback]
                mc = np.sort(loc.soma16)
                gdm_, gfm_ = mmodel.predict(mfeat.rows(loc, mfb, candidates=mc))
                mchoice, mwhy = residual_choice(mc, mfb, gdm_, gfm_)
                gdq, gfq = pmodel.predict(feat.rows(loc, fb))
                pchoice, _ = residual_choice(loc.cand, fb, gdq, gfq)
                if why == "switched":
                    others = [int(d) for d in loc.cand if int(d) != fb]
                    rchoice = others[np.random.default_rng(stable_seed("switch_random", lib.ids[q])).integers(len(others))]
                else:
                    rchoice = fb
                cp = copies.choose(q)
                donors = {SELECTOR: choice, "adaptive_fallback": fb, "metadata_only": mchoice,
                          "permuted_labels": pchoice, "switch_random": rchoice, "observation_only": None,
                          **{m: loc.choice[m] for m in FALLBACK_SET}, **{m: cp[m] for m in COPIES}}
                forecasts[lib.ids[q]] = {"outer": k, "why": why, "metadata_why": mwhy, "obs_sha": idx.obs_sha[q],
                                         "donors": {m: (None if d is None else lib.ids[d]) for m, d in donors.items()}}
        write_json(out / "fold_info.json", finite_json(fold_info))

        # ---- decode the forecasts (no truth), hash, seal ----
        observations = dict(pool.imap_unordered(observation_job, list(forecasts), chunksize=8))
        djobs = [{"job": f"decode|{qid}", "query": qid, "obs_sha": f["obs_sha"], "observation": observations[qid],
                  "donors": sorted(set(f["donors"].values()), key=lambda d: (d is not None, d or ""))}
                 for qid, f in sorted(forecasts.items())]
        run_jobs(decode_job, djobs, out / "decoded.jsonl", pool, log, "decode forecasts")
        decoded = {}
        for line in (out / "decoded.jsonl").read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                decoded[r["query"]] = r["decoded"]
        for qid, f in forecasts.items():
            f["completion_sha256"] = {m: decoded[qid]["none" if d is None else d]["completion_sha256"] for m, d in f["donors"].items()}
        fsha = write_json(out / "forecasts.json", forecasts)
        seal = out / "FORECASTS_SEALED.json"
        if seal.exists():
            if read_json(seal)["forecasts_sha256"] != fsha:
                raise RuntimeError("sealed forecasts differ on resume")
        else:
            write_json(seal, {"utc": utc(), "forecasts_sha256": fsha, "queries": len(forecasts)}, exclusive=True)
        log(f"forecasts sealed ({len(forecasts)} queries, sha {fsha[:12]})")

        # ---- score after the seal ----
        sjobs = [{"job": f"score|{qid}", "query": qid, "obs_sha": f["obs_sha"],
                  "donors": sorted(set(f["donors"].values()), key=lambda d: (d is not None, d or ""))}
                 for qid, f in sorted(forecasts.items())]
        scores = Labels(lib)
        run_jobs(label_job, sjobs, out / "labels_score.jsonl", pool, log, "score forecasts", store=scores)
    return summarise(lib, idx, forecasts, scores, fold_info, out, log)


def summarise(lib, idx, forecasts, scores, fold_info, out, log):
    per_method, validity, executed = {}, {}, {}
    for m in METHODS:
        rows = {}
        bad = Counter()
        invalid_diag = []
        cable_err = defaultdict(list)
        merged = 0.0
        modes = Counter()
        ratios = []
        for qid, f in forecasts.items():
            q = lib.index[qid]
            d = f["donors"][m]
            r = scores.get(q, None if d is None else lib.index[d])
            if r.get("completion_sha256") is not None and r["completion_sha256"] != f["completion_sha256"][m]:
                raise RuntimeError(f"{qid} {m}: scored completion differs from the sealed forecast")
            rows[q] = r
            if not r["valid"]:
                bad[r["status"]] += 1
                invalid_diag.append({"query": qid, "status": r["status"], "distance_um": r.get("distance_um"),
                                     "f1_200": r.get("f1_200"), "duplicate_um": r.get("duplicate_um")})
            for k2, v in r.get("modes", {}).items():
                modes[k2] += v
            merged += float(r.get("merged_duplicate_cable_um") or 0.0)
            if r.get("hidden_cable_um"):
                ratios.append((r.get("extra_cable_um") or 0.0) / r["hidden_cable_um"])
                # every forecast counts (an empty one has error 1; a guard failure, with no cable recorded, inf)
                if r.get("extra_cable_um") is None:
                    cable_err[lib.animal[q]].append(math.inf)
                else:
                    cable_err[lib.animal[q]].append(abs(r["extra_cable_um"] - r["hidden_cable_um"]) / r["hidden_cable_um"])
            else:
                cable_err[lib.animal[q]].append(math.inf)
        per_method[m] = rows
        ratios = np.array(ratios)
        validity[m] = {"invalid": int(sum(bad.values())), "by_status": dict(bad), "invalid_diagnostics": invalid_diag,
                       "inflated_share": float((ratios > INFLATED_RATIO).mean()) if len(ratios) else 0.0,
                       "median_cable_ratio": float(np.median(ratios)) if len(ratios) else 0.0,
                       "mean_log2_cable_ratio": float(np.mean(np.log2(np.maximum(ratios, 1e-12)))) if len(ratios) else 0.0,
                       "mean_abs_rel_cable_error": float(np.mean([np.mean(v) for v in cable_err.values()])),
                       "merged_duplicate_cable_um": merged}
        executed[m] = dict(modes)
    means = {m: sample_means(lib, rows) for m, rows in per_method.items()}
    animals = sorted(means[SELECTOR])
    summary = {m: {"distance_um": float(np.mean([means[m][a][0] for a in animals])),
                   "f1_200": float(np.mean([means[m][a][1] for a in animals]))} for m in METHODS}
    versus = {c: paired_bootstrap(means[SELECTOR], means[c], animals) for c in CONTROLS + DIAGNOSTIC}
    switch_rate = float(np.mean([f["why"] == "switched" for f in forecasts.values()]))
    v = verdict(summary, versus, validity, switch_rate)
    result = {"utc": utc(), "neurons": len(forecasts), "samples": len(animals), "summary": summary, "versus": versus,
              "validity": validity, "executed_grafts": executed, "verdict": v, "fold_info": fold_info,
              "why_counts": dict(Counter(f["why"] for f in forecasts.values())),
              "metadata_why_counts": dict(Counter(f["metadata_why"] for f in forecasts.values()))}
    write_json(out / "results.json", finite_json(result))
    log(f"{v['sentence']}")
    return result


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "smoke"):
        p = sub.add_parser(name)
        p.add_argument("--workers", type=int, default=4)
        p.add_argument("--out", type=Path, default=None)
        p.add_argument("--animals", type=int, default=None)
    args = ap.parse_args(argv)
    for v in THREAD_VARS:
        if os.environ.get(v) != "1":
            raise SystemExit(f"set {v}=1")
    if args.cmd == "smoke":
        if args.out is None or OUT_DEFAULT in [args.out, *args.out.parents]:
            raise SystemExit("smoke writes only to a scratch folder outside data/interim/axon_completion_v1")
        _DEV["allow_unpinned_contract"] = True
        log = Log(args.out / "run.log")
        experiment(args.out, args.workers, log, animals=args.animals or 15, dev=True)
        return
    out = OUT_DEFAULT
    if args.out is not None or args.animals is not None:
        raise SystemExit("run uses the fixed output folder and the whole TRAIN cohort")
    if CONTRACT_SHA is None:
        raise SystemExit("the contract is not pinned: wait for ROOT's contract check")
    if "STATUS: FINAL" not in PROTOCOL.read_text():
        raise SystemExit("protocol is not FINAL")
    owned = ["analysis/axon_completion_v1.py", "tests/test_axon_completion_v1.py", "reports/axon_completion_v1_protocol.md",
             "analysis/axon_completion_contract.py", "tests/test_axon_completion_contract.py", "analysis/axon_completion_preflight.py"]
    clean, head = git_clean(owned)
    if not clean:
        raise SystemExit("owned files must be committed and clean before the run")
    out.mkdir(parents=True, exist_ok=True)
    if (out / "RUN_COMPLETE.json").exists():
        raise SystemExit("RUN_COMPLETE exists: never repeated")
    log = Log(out / "run.log")
    pins = code_pins()
    started = out / "RUN_STARTED.json"
    if started.exists():
        if read_json(started)["code_sha256"] != pins:
            raise SystemExit("code/protocol/contract changed since RUN_STARTED; resume refused")
        log("resuming")
    else:
        write_json(started, {"utc": utc(), "git_head": head, "code_sha256": pins, "workers": args.workers,
                             "salt": SALT, "seed": SEED}, exclusive=True)
    try:
        experiment(out, args.workers, log)
    except leak_types() as error:
        write_json(out / "RUN_FAILED.json", {"utc": utc(), "reason": str(error)}, exclusive=True)
        raise
    write_json(out / "RUN_COMPLETE.json", {"utc": utc(), "results_sha256": sha_file(out / "results.json")}, exclusive=True)


if __name__ == "__main__":
    main()
