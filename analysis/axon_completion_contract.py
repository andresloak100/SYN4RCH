"""Observation and missing-cable scoring for a NEW partial-axon task.

Only ``Observation`` goes to a predictor. ``split_arbor`` and ``MissingTruth``
belong to acquisition/evaluation. The original full-Arbor metric is untouched.
Coordinates stay in the input frame; no soma translation is performed.

The 50-um sampler uses a soma-path grid, with component/fork/tip anchors and
trapezoidal ARC-LENGTH mass. Distance is the symmetric weighted nearest-sample
distance and F1 uses a 200-um tolerance, as in the old metric conventions. The
new metric scores suffix forests, not a fictitious tree connecting their roots.
It is not numerically interchangeable with the old full-Arbor metric.
"""
from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
import hashlib
import json
import math

import numpy as np
from scipy.spatial import cKDTree

BUDGET_UM = 300.0
STEP_UM = 50.0
TOLERANCE_UM = 200.0
_EPS = 1e-9


def _array(value, dtype):
    source = np.array(value, dtype=dtype, copy=True)
    # Immutable bytes backing also refuses setflags(write=True), not just writes.
    return np.frombuffer(source.tobytes(), dtype=source.dtype).reshape(source.shape)


def _tree_arrays(value):
    xyz = np.asarray(value.xyz, dtype=float)
    parent = np.asarray(value.parent)
    if xyz.ndim != 2 or xyz.shape[1:] != (3,) or len(xyz) == 0 or not np.isfinite(xyz).all():
        raise ValueError("A finite, nonempty XYZ tree is required")
    if parent.shape != (len(xyz),) or not np.issubdtype(parent.dtype, np.integer):
        raise ValueError("One integer parent per vertex is required")
    if parent[0] != -1 or np.any(parent[1:] < 0) or np.any(parent[1:] >= np.arange(1, len(xyz))):
        raise ValueError("One root and topological parent order are required")
    return xyz, parent


def _depth(xyz, parent):
    out = np.zeros(len(xyz))
    lengths = np.linalg.norm(xyz[1:] - xyz[parent[1:]], axis=1)
    for i in range(1, len(xyz)):
        out[i] = out[parent[i]] + lengths[i - 1]
    return out


def _children(parent):
    children = [[] for _ in parent]
    for i, p in enumerate(parent):
        if p >= 0:
            children[int(p)].append(i)
    return children


@dataclass(frozen=True, slots=True)
class Observation:
    """Measured prefix only. Leaves are not labelled as true biological tips."""

    xyz: np.ndarray
    parent: np.ndarray
    path_budget_um: float = BUDGET_UM

    def __post_init__(self):
        xyz, parent = _tree_arrays(self)
        budget = float(self.path_budget_um)
        if not np.isfinite(budget) or budget <= 0:
            raise ValueError("Positive fixed path budget required")
        if _depth(xyz, parent).max() > budget + _EPS:
            raise ValueError("Observation extends beyond its acquisition budget")
        object.__setattr__(self, "xyz", _array(xyz, float))
        object.__setattr__(self, "parent", _array(parent, np.int64))
        object.__setattr__(self, "path_budget_um", budget)

    def to_dict(self):
        return {"schema": "axon_observation_v1", "xyz_um": self.xyz.tolist(),
                "parent": self.parent.tolist(), "path_budget_um": self.path_budget_um}

    @classmethod
    def from_dict(cls, data):
        if set(data) != {"schema", "xyz_um", "parent", "path_budget_um"} or data["schema"] != "axon_observation_v1":
            raise ValueError("Unexpected observation schema or hidden fields")
        return cls(np.asarray(data["xyz_um"]), np.asarray(data["parent"]), data["path_budget_um"])

    @property
    def sha256(self):
        raw = json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
        return hashlib.sha256(raw).hexdigest()

    def frontier_indices(self):
        """Acquisition-boundary vertices, derived only from visible geometry."""
        return np.flatnonzero(np.abs(_depth(self.xyz, self.parent) - self.path_budget_um) <= _EPS)


@dataclass(frozen=True, slots=True)
class Completion:
    """Append-only forecast. New parents index [observed vertices, new vertices]."""

    observation: Observation
    extra_xyz: np.ndarray
    extra_parent: np.ndarray

    def __post_init__(self):
        if not isinstance(self.observation, Observation):
            raise TypeError("A public Observation is required")
        xyz = np.asarray(self.extra_xyz, float)
        parent = np.asarray(self.extra_parent)
        if xyz.size == 0:
            xyz = xyz.reshape(0, 3)
        if xyz.ndim != 2 or xyz.shape[1:] != (3,) or not np.isfinite(xyz).all():
            raise ValueError("New XYZ must be finite (n,3)")
        if parent.shape != (len(xyz),) or not np.issubdtype(parent.dtype, np.integer):
            raise ValueError("New parents must be integer indices")
        n = len(self.observation.xyz)
        if np.any(parent < 0) or np.any(parent >= n + np.arange(len(xyz))):
            raise ValueError("Every new node must attach to an earlier observed/generated node")
        object.__setattr__(self, "extra_xyz", _array(xyz, float))
        object.__setattr__(self, "extra_parent", _array(parent, np.int64))

    @property
    def xyz(self):
        return np.vstack([self.observation.xyz, self.extra_xyz])

    @property
    def parent(self):
        return np.r_[self.observation.parent, self.extra_parent]

    def to_dict(self):
        return {"schema": "axon_completion_v1", "observation_sha256": self.observation.sha256,
                "extra_xyz_um": self.extra_xyz.tolist(), "extra_parent": self.extra_parent.tolist()}

    @classmethod
    def from_dict(cls, observation, data):
        if set(data) != {"schema", "observation_sha256", "extra_xyz_um", "extra_parent"} or data["schema"] != "axon_completion_v1":
            raise ValueError("Unexpected completion schema")
        if data["observation_sha256"] != observation.sha256:
            raise ValueError("Completion belongs to another observation")
        parent = np.asarray(data["extra_parent"]) if data["extra_parent"] else np.empty(0, dtype=int)
        return cls(observation, np.asarray(data["extra_xyz_um"]).reshape(-1, 3), parent)


def observation_only(observation):
    return Completion(observation, np.empty((0, 3)), np.empty(0, dtype=int))


def from_arbor(observation, arbor):
    """Adapter only: the complete observed block must be retained byte-for-value."""
    xyz, parent = _tree_arrays(arbor)
    n = len(observation.xyz)
    if len(xyz) < n or not np.array_equal(xyz[:n], observation.xyz) or not np.array_equal(parent[:n], observation.parent):
        raise ValueError("Observed coordinates/topology were changed or omitted")
    return Completion(observation, xyz[n:], parent[n:])


@dataclass(frozen=True, slots=True)
class CableForest:
    """Evaluator/sampler geometry, with no artificial edges between components."""

    xyz: np.ndarray
    parent: np.ndarray
    path_um: np.ndarray

    def __post_init__(self):
        xyz, parent, depth = np.asarray(self.xyz, float), np.asarray(self.parent), np.asarray(self.path_um, float)
        if xyz.shape != (len(parent), 3) or depth.shape != (len(parent),) or not np.issubdtype(parent.dtype, np.integer):
            raise ValueError("Invalid cable-forest shapes")
        if not np.isfinite(xyz).all() or not np.isfinite(depth).all() or np.any(depth < 0):
            raise ValueError("Invalid cable-forest geometry")
        linked = parent >= 0
        if np.any(parent < -1) or np.any(parent[linked] >= np.flatnonzero(linked)):
            raise ValueError("Forest parents must be roots or earlier vertices")
        lengths = np.linalg.norm(xyz[linked] - xyz[parent[linked]], axis=1)
        if not np.allclose(depth[linked] - depth[parent[linked]], lengths, atol=1e-8, rtol=1e-10):
            raise ValueError("Forest path coordinate disagrees with its cable")
        object.__setattr__(self, "xyz", _array(xyz, float))
        object.__setattr__(self, "parent", _array(parent, np.int64))
        object.__setattr__(self, "path_um", _array(depth, float))

    def cable(self):
        linked = self.parent >= 0
        return float(np.linalg.norm(self.xyz[linked] - self.xyz[self.parent[linked]], axis=1).sum())


def extra_forest(completion):
    """All newly generated edges, including connectors near the known prefix."""
    xyz, parent = completion.xyz, completion.parent
    depth = _depth(xyz, parent)
    n = len(completion.observation.xyz)
    extra_parent = parent[n:]
    anchors = np.unique(extra_parent[extra_parent < n])
    count = len(anchors)
    observed_map = np.full(n, -1, dtype=int)
    observed_map[anchors] = np.arange(count)
    new_parent = np.empty(len(extra_parent), dtype=int)
    attached = extra_parent < n
    new_parent[attached] = observed_map[extra_parent[attached]]
    new_parent[~attached] = count + extra_parent[~attached] - n
    return CableForest(np.vstack([xyz[anchors], xyz[n:]]), np.r_[np.full(count, -1, dtype=int), new_parent],
                       np.r_[depth[anchors], depth[n:]])


@dataclass(frozen=True, slots=True)
class MissingTruth:
    """Evaluator-only object. Never a prediction argument or feature source."""

    observation_sha256: str
    forest: CableForest

    @property
    def eligible(self):
        return self.forest.cable() > _EPS


def _canonical_visible(xyz, parent):
    """Visible-only order; hidden/source row indices cannot affect public IDs."""
    xyz, parent = np.asarray(xyz, float), np.asarray(parent, int)
    children = _children(parent)
    signatures = [None] * len(xyz)
    for i in range(len(xyz) - 1, -1, -1):
        raw = xyz[i].astype("<f8").tobytes() + b"".join(sorted(signatures[c] for c in children[i]))
        signatures[i] = hashlib.sha256(raw).digest()
    order, queue = [], deque([0])
    while queue:
        i = queue.popleft(); order.append(i)
        queue.extend(sorted(children[i], key=lambda c: signatures[c]))
    remap = {old: new for new, old in enumerate(order)}
    return xyz[order], np.asarray([-1] + [remap[int(parent[i])] for i in order[1:]]), remap


def split_arbor(arbor, budget_um=BUDGET_UM):
    """Evaluator/acquisition only: exact edge split at a fixed root-path budget.

    No outgoing edge, even a zero-length marker edge, is exposed beyond a
    boundary vertex. Thus hidden branch multiplicity cannot leak through zero
    length children at the cut. Zero-length hidden edges carry zero score mass.
    """
    xyz, parent = _tree_arrays(arbor)
    budget = float(budget_um)
    if not np.isfinite(budget) or budget <= 0:
        raise ValueError("Positive fixed path budget required")
    depth = _depth(xyz, parent)
    observed_xyz, observed_parent = [xyz[0]], [-1]
    original_map, boundary_map = {0: 0}, {}
    for i in range(1, len(xyz)):
        p = int(parent[i])
        if p not in original_map or depth[p] >= budget:
            continue
        if depth[i] <= budget:
            original_map[i] = len(observed_xyz)
            observed_xyz.append(xyz[i]); observed_parent.append(original_map[p])
        else:
            fraction = (budget - depth[p]) / (depth[i] - depth[p])
            boundary_map[i] = len(observed_xyz)
            observed_xyz.append(xyz[p] + fraction * (xyz[i] - xyz[p]))
            observed_parent.append(original_map[p])
    ox, op, remap = _canonical_visible(observed_xyz, observed_parent)
    observation = Observation(ox, op, budget)
    original_map = {k: remap[v] for k, v in original_map.items()}
    boundary_map = {k: remap[v] for k, v in boundary_map.items()}
    extra_xyz, extra_parent = [], []
    n = len(observation.xyz)
    for i in range(1, len(xyz)):
        if i in original_map:
            continue
        p = boundary_map[i] if i in boundary_map else original_map[int(parent[i])]
        original_map[i] = n + len(extra_xyz)
        extra_xyz.append(xyz[i]); extra_parent.append(p)
    completion = Completion(observation, np.asarray(extra_xyz, float).reshape(-1, 3), np.asarray(extra_parent, int))
    return observation, MissingTruth(observation.sha256, extra_forest(completion))


def observe(arbor, budget_um=BUDGET_UM):
    """Acquisition convenience: return the public prefix and discard truth."""
    return split_arbor(arbor, budget_um)[0]


def registration(observation, table_soma):
    """Known registration discrepancy; never translates measured observations."""
    table_soma = np.asarray(table_soma, float)
    if table_soma.shape != (3,) or not np.isfinite(table_soma).all():
        raise ValueError("Finite table soma XYZ required")
    delta = table_soma - observation.xyz[0]
    return {"trace_root_xyz_um": observation.xyz[0].tolist(), "table_soma_um": table_soma.tolist(),
            "table_minus_trace_xyz_um": delta.tolist(), "distance_um": float(np.linalg.norm(delta))}


def sample_forest(forest, step_um=STEP_UM):
    """Root-path grid, branch/tip anchors, exact arc-length quadrature mass."""
    step = float(step_um)
    if not np.isfinite(step) or step <= 0:
        raise ValueError("Positive sampling step required")
    children = _children(forest.parent)
    linked = forest.parent >= 0
    edge_lengths = np.zeros(len(forest.parent))
    edge_lengths[linked] = np.linalg.norm(forest.xyz[linked] - forest.xyz[forest.parent[linked]], axis=1)
    critical = {i for i, p in enumerate(forest.parent) if p == -1 or len(children[i]) != 1}
    points, masses = [], []
    for start in sorted(critical):
        for child in children[start]:
            run = [start, child]
            while run[-1] not in critical:
                run.append(children[run[-1]][0])
            curve = forest.xyz[run]
            origin = float(forest.path_um[start])
            first = (math.floor(origin / step) + 1) * step - origin
            if len(run) == 2:
                total = float(edge_lengths[child])
                if total <= 0:
                    continue
                if first >= total:
                    points.extend((curve[0], curve[1])); masses.extend((total / 2, total / 2))
                    continue
                if first > 0 and first + step >= total:
                    mid = curve[0] + (first / total) * (curve[1] - curve[0])
                    points.extend((curve[0], mid, curve[1]))
                    masses.extend((first / 2, total / 2, (total - first) / 2))
                    continue
            cumulative = np.r_[0., np.cumsum(np.linalg.norm(np.diff(curve, axis=0), axis=1))]
            total = float(cumulative[-1])
            if total <= 0:
                continue
            grid = np.arange(first, total, step)
            locations = np.unique(np.r_[0., grid[(grid > 0) & (grid < total)], total])
            sampled = np.column_stack([np.interp(locations, cumulative, curve[:, k]) for k in range(3)])
            weight = np.zeros(len(locations)); gaps = np.diff(locations)
            weight[:-1] += gaps / 2; weight[1:] += gaps / 2
            points.extend(sampled); masses.extend(weight)
    return np.asarray(points, float).reshape(-1, 3), np.asarray(masses, float)


def _exact_duplicate_lengths(completion):
    """Endpoint-exact edge duplicates, including reversed and observed copies.

    This deliberately does not claim to detect every partially overlapping or
    differently segmented curve. Cable inflation must remain a separate gate.
    """
    xyz, parent = completion.xyz, completion.parent
    n = len(completion.observation.xyz)
    lengths = np.linalg.norm(xyz[1:] - xyz[parent[1:]], axis=1)
    seen, observed, duplicate, copied_observed = set(), set(), 0., 0.
    for i in range(1, len(xyz)):
        a, b = tuple(xyz[parent[i]].tolist()), tuple(xyz[i].tolist())
        length = float(lengths[i - 1])
        if length <= 0:
            continue
        key = tuple(sorted((a, b)))
        if i < n:
            observed.add(key)
        elif key in seen:
            duplicate += length
            if key in observed:
                copied_observed += length
        seen.add(key)
    return duplicate, copied_observed


def score_completion(prediction, truth, *, step_um=STEP_UM, tolerance_um=TOLERANCE_UM):
    """Score NEW cable only. Failed/empty forecasts stay in the denominator."""
    if not isinstance(prediction, Completion) or not isinstance(truth, MissingTruth):
        raise TypeError("Completion and evaluator-only MissingTruth required")
    if prediction.observation.sha256 != truth.observation_sha256:
        raise ValueError("Prediction changed or substituted the measured prefix")
    tolerance = float(tolerance_um)
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("Positive scoring tolerance required")
    suffix = extra_forest(prediction)
    true_length, predicted_length = truth.forest.cable(), suffix.cable()
    duplicate, copied_observed = _exact_duplicate_lengths(prediction)
    out = {"eligible": truth.eligible, "hidden_cable_um": true_length, "extra_cable_um": predicted_length,
           "empty_forecast": predicted_length <= _EPS, "duplicate_extra_cable_um": duplicate,
           "copied_observed_extra_cable_um": copied_observed, "nonduplicate_valid": duplicate <= _EPS,
           "step_um": float(step_um), "tolerance_um": tolerance}
    if not truth.eligible:
        out.update(status="ineligible_no_hidden_cable", valid_for_comparison=False, distance_um=None,
                   precision_200=None, recall_200=None, f1_200=None)
        return out
    out.update(cable_signed_relative_error=(predicted_length - true_length) / true_length,
               cable_relative_error=abs(predicted_length - true_length) / true_length)
    if predicted_length <= _EPS:
        out.update(status="failed_empty_forecast", valid_for_comparison=False, distance_um=math.inf,
                   prediction_to_truth_um=math.inf, truth_to_prediction_um=math.inf,
                   precision_200=0., recall_200=0., f1_200=0., matched_extra_cable_200_um=0.,
                   unmatched_extra_cable_200_um=0., covered_hidden_cable_200_um=0.,
                   uncovered_hidden_cable_200_um=true_length)
        return out
    px, pw = sample_forest(suffix, step_um)
    tx, tw = sample_forest(truth.forest, step_um)
    pd = cKDTree(tx).query(px)[0]; td = cKDTree(px).query(tx)[0]
    forward, backward = float(np.average(pd, weights=pw)), float(np.average(td, weights=tw))
    precision = float(np.average(pd <= tolerance, weights=pw))
    recall = float(np.average(td <= tolerance, weights=tw))
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.
    out.update(status="scored" if out["nonduplicate_valid"] else "failed_duplicate_extra_cable",
               valid_for_comparison=out["nonduplicate_valid"], distance_um=(forward + backward) / 2,
               prediction_to_truth_um=forward, truth_to_prediction_um=backward,
               precision_200=precision, recall_200=recall, f1_200=f1,
               matched_extra_cable_200_um=precision * predicted_length,
               unmatched_extra_cable_200_um=(1 - precision) * predicted_length,
               covered_hidden_cable_200_um=recall * true_length,
               uncovered_hidden_cable_200_um=(1 - recall) * true_length)
    return out


def metrics_to_json(metrics):
    """Strict-JSON adapter: null never conceals an empty/failed forecast."""
    out = dict(metrics)
    nonfinite = [k for k, v in out.items() if isinstance(v, (float, np.floating)) and not np.isfinite(v)]
    for k in nonfinite:
        out[k] = None
    if nonfinite:
        out["nonfinite_metric_fields"] = nonfinite
        out["valid_for_comparison"] = False
    json.dumps(out, allow_nan=False)
    return out
