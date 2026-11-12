"""D-110: exploratory animal-held-out prediction of explicit axon geometry.

The fitted object accepts only training arbors and soma XYZ. Query geometry is
loaded by score_predictions(), after every prediction is written. See the
predeclared protocol in reports/d110_protocol.md. No network or git writes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter, deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data/interim/d110_geometry"
SCRIPT = Path(__file__).resolve()
PROTOCOL = ROOT / "reports/d110_protocol.md"
VALIDATION = ROOT / "data/interim/swc_validation.json"
RAW = ROOT / "data/raw/jiao2025/swc"
CONFIG = {
    "step_um": 50.0, "neighbors": 30, "reservoir_per_animal": 512,
    "position_scale_um": 1000.0, "path_scale_um": 3000.0,
    "samples": 5, "max_nodes": 5000, "seed": 20261110,
    "f1_tolerance_um": 200.0, "bootstrap_replicates": 10000,
    "frame": "mirrored-left", "coordinate_system": "CCFv3_AP_DV_ML_um",
    "model": "animal_balanced_empirical_branching_growth_v1",
}


def canonical(obj):
    return json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n"


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(canonical(obj))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stable_seed(*parts):
    raw = "/".join(map(str, (CONFIG["seed"],) + parts)).encode()
    return int.from_bytes(hashlib.sha256(raw).digest()[:4], "little")


@dataclass
class Arbor:
    xyz: np.ndarray
    parent: np.ndarray
    types: np.ndarray

    def children(self):
        out = [[] for _ in self.parent]
        for i, p in enumerate(self.parent):
            if p >= 0:
                out[int(p)].append(i)
        return out

    def validate(self):
        n = len(self.xyz)
        assert self.xyz.shape == (n, 3) and self.parent.shape == (n,)
        assert self.types.shape == (n,) and n >= 1 and np.isfinite(self.xyz).all()
        assert self.parent[0] == -1 and np.count_nonzero(self.parent == -1) == 1
        for i, p in enumerate(self.parent[1:], 1):
            assert 0 <= p < i, "topological parent ordering is required"
        return self

    def cable(self):
        if len(self.xyz) == 1:
            return 0.0
        return float(np.linalg.norm(self.xyz[1:] - self.xyz[self.parent[1:]], axis=1).sum())

    def branches(self):
        # The root is a soma, not an axonal bifurcation.
        return sum(len(c) >= 2 for c in self.children()[1:])

    def path_lengths(self):
        out = np.zeros(len(self.xyz))
        for i, p in enumerate(self.parent[1:], 1):
            out[i] = out[p] + np.linalg.norm(self.xyz[i] - self.xyz[p])
        return out

    def weights(self):
        # Trapezoidal arc-length mass at sample nodes.
        out = np.zeros(len(self.xyz))
        if len(self.xyz) == 1:
            out[0] = 1.0
            return out
        lengths = np.linalg.norm(self.xyz[1:] - self.xyz[self.parent[1:]], axis=1)
        out[1:] += lengths / 2
        np.add.at(out, self.parent[1:], lengths / 2)
        if out.sum() <= 0:
            out[0] = 1.0
        return out

    def shifted(self, soma):
        return Arbor(self.xyz - self.xyz[0] + np.asarray(soma), self.parent.copy(), self.types.copy())


def save_arbor(path, arbor):
    arbor.validate()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, xyz_um=arbor.xyz, parent=arbor.parent, types=arbor.types)


def load_arbor(path):
    with np.load(path, allow_pickle=False) as d:
        return Arbor(d["xyz_um"], d["parent"], d["types"]).validate()


def extract_explicit_axon(path, reflected=False):
    """Retain type-2 nodes and their ancestor connector paths to the soma.

    Undocumented marker ancestors remain connector geometry. A type-0/3 path
    into an explicit axon is refused, not silently reclassified.
    """
    rows = np.loadtxt(path, comments="#", ndmin=2)
    assert rows.shape[1] == 7 and np.isfinite(rows).all()
    ids = rows[:, 0].astype(int)
    parent_ids = rows[:, 6].astype(int)
    types = rows[:, 1].astype(int)
    assert len(np.unique(ids)) == len(ids)
    pos = {int(v): i for i, v in enumerate(ids)}
    roots = np.flatnonzero(parent_ids == -1)
    assert len(roots) == 1 and types[roots[0]] == 1
    root = int(roots[0])
    keep = {root}
    explicit = np.flatnonzero(types == 2)
    assert len(explicit) > 0, "primary data require explicit type-2 axons"
    for i in explicit:
        seen = set()
        while int(i) not in keep:
            i = int(i)
            assert i not in seen, "cycle"
            seen.add(i)
            assert types[i] in (1, 2, 10, 11, 12), "ambiguous ancestor in explicit axon"
            keep.add(i)
            p = int(parent_ids[i])
            assert p in pos, "missing axon ancestor"
            i = pos[p]
    # A deterministic topological traversal, independent of raw row ordering.
    children = {i: [] for i in keep}
    for i in keep - {root}:
        children[pos[int(parent_ids[i])]].append(i)
    order, queue = [], deque([root])
    while queue:
        i = queue.popleft()
        order.append(i)
        queue.extend(sorted(children[i], key=lambda k: ids[k]))
    assert len(order) == len(keep)
    remap = {old: new for new, old in enumerate(order)}
    xyz = rows[order, 2:5].copy()
    if reflected:
        xyz[:, 2] = 11400.0 - xyz[:, 2]
    parent = np.array([-1] + [remap[pos[int(parent_ids[i])]] for i in order[1:]], dtype=int)
    return Arbor(xyz, parent, types[order]).validate()


def resample(arbor, step=None):
    """Resample branch-to-branch polylines, preserving every fork and tip."""
    step = float(step or CONFIG["step_um"])
    children = arbor.children()
    critical = {0} | {i for i, c in enumerate(children) if len(c) != 1}
    xyz, parent, types = [arbor.xyz[0].copy()], [-1], [1]
    anchors = {0: 0}
    queue = deque([0])
    while queue:
        start = queue.popleft()
        for child in children[start]:
            run = [start, child]
            while run[-1] not in critical:
                run.append(children[run[-1]][0])
            points = arbor.xyz[run]
            d = np.linalg.norm(np.diff(points, axis=0), axis=1)
            cum = np.r_[0.0, np.cumsum(d)]
            total = float(cum[-1])
            dist = np.r_[np.arange(step, total, step), total] if total > 1e-10 else np.array([0.0])
            previous = anchors[start]
            for t in dist:
                p = np.array([np.interp(t, cum, points[:, k]) for k in range(3)])
                xyz.append(p)
                parent.append(previous)
                types.append(2)
                previous = len(xyz) - 1
            anchors[run[-1]] = previous
            if children[run[-1]]:
                queue.append(run[-1])
    return Arbor(np.asarray(xyz), np.asarray(parent, dtype=int), np.asarray(types, dtype=int)).validate()


def unit(v):
    norm = np.linalg.norm(v)
    return np.asarray(v, dtype=float) / norm if norm > 1e-12 else np.zeros(3)


def rotate_match(vectors, old, new):
    """Rigidly transport a training turn/fork to the current incoming tangent."""
    vectors = np.asarray(vectors, dtype=float).reshape(-1, 3)
    old, new = unit(old), unit(new)
    if np.linalg.norm(old) == 0 or np.linalg.norm(new) == 0:
        return vectors.copy()
    cross = np.cross(old, new)
    sine, cosine = np.linalg.norm(cross), float(np.clip(old @ new, -1, 1))
    if sine < 1e-10:
        if cosine > 0:
            return vectors.copy()
        basis = np.eye(3)[np.argmin(np.abs(old))]
        axis = unit(np.cross(old, basis))
        return -vectors + 2 * (vectors @ axis)[:, None] * axis
    axis = cross / sine
    return (vectors * cosine + np.cross(axis, vectors) * sine
            + (vectors @ axis)[:, None] * axis * (1 - cosine))


def state(soma, point, incoming, path_um):
    return np.r_[np.asarray(soma) / CONFIG["position_scale_um"],
                 np.asarray(point) / CONFIG["position_scale_um"], incoming,
                 path_um / CONFIG["path_scale_um"]]


class GrowthModel:
    """A conditional empirical transition distribution with an animal-balanced prior.

    The fitted reservoir contains exactly 512 transitions per training animal.
    Within an animal, a neuron then one of its events are sampled uniformly.
    Prediction takes a 3-vector soma only: no endpoint/length/topology arguments.
    """
    def __init__(self, training, seed):
        assert training
        self.training_ids = tuple(t["neuron"] for t in training)
        self.training_animals = tuple(sorted({t["animal"] for t in training}))
        rng = np.random.default_rng(seed)
        events, root_events, stats, by_animal = [], [], [], {}
        self.transition_counts = {}
        for t in training:
            a = t["arbor"].validate()
            children, depth = a.children(), a.path_lengths()
            soma = a.xyz[0]
            local = []
            for i in range(1, len(a.xyz)):
                inc = unit(a.xyz[i] - a.xyz[a.parent[i]])
                outgoing = np.array([unit(a.xyz[j] - a.xyz[i]) for j in children[i]]).reshape(-1, 3)
                local.append((state(soma, a.xyz[i], inc, depth[i]), inc, outgoing))
            assert local
            by_animal.setdefault(t["animal"], []).append(local)
            self.transition_counts[t["neuron"]] = len(local)
            root_events.append((soma.copy(), np.array([unit(a.xyz[j] - soma) for j in children[0]]).reshape(-1, 3), t["animal"]))
            stats.append((t["animal"], a.cable(), float(depth.max())))
        for animal in self.training_animals:
            neurons = by_animal[animal]
            for _ in range(CONFIG["reservoir_per_animal"]):
                ev = neurons[int(rng.integers(len(neurons)))]
                events.append(ev[int(rng.integers(len(ev)))])
        self.states = np.array([e[0] for e in events])
        self.events = events
        self.tree = cKDTree(self.states)
        self.roots = root_events
        self.root_tree = cKDTree(np.array([e[0] for e in root_events]) / CONFIG["position_scale_um"])
        count = Counter(e[2] for e in root_events)
        self.root_prior = np.array([1 / count[e[2]] for e in root_events])
        # Each animal contributes one within-animal median to each budget.
        self.cable_budget = float(np.median([np.median([s[1] for s in stats if s[0] == g]) for g in self.training_animals]))
        self.path_budget = float(np.median([np.median([s[2] for s in stats if s[0] == g]) for g in self.training_animals]))
        self.cable_budget = max(CONFIG["step_um"], self.cable_budget)
        self.path_budget = max(CONFIG["step_um"], self.path_budget)

    @staticmethod
    def _choose(tree, x, rng, prior=None):
        d, idx = tree.query(x, k=min(CONFIG["neighbors"], tree.n))
        d, idx = np.atleast_1d(d), np.atleast_1d(idx).astype(int)
        # Stable conditional Gaussian kernel; the common minimum cancels.
        logw = -0.5 * (d * d - d.min() ** 2)
        w = np.exp(np.maximum(logw, -700))
        if prior is not None:
            w *= prior[idx]
        w /= w.sum()
        return int(rng.choice(idx, p=w))

    def predict(self, soma, seed):
        soma = np.asarray(soma, dtype=float)
        assert soma.shape == (3,) and np.isfinite(soma).all()
        rng = np.random.default_rng(seed)
        root_index = self._choose(self.root_tree, soma / CONFIG["position_scale_um"], rng, self.root_prior)
        xyz, parent = [soma.copy()], [-1]
        active = deque()
        cable = 0.0
        flags = {"node_cap_reached": False, "cable_budget_reached": False, "path_budget_reached": False}

        def extend(p, direction, old_depth):
            nonlocal cable
            remaining = min(CONFIG["step_um"], self.cable_budget - cable, self.path_budget - old_depth)
            if len(xyz) >= CONFIG["max_nodes"]:
                flags["node_cap_reached"] = True
                return
            if remaining <= 1e-8:
                flags["cable_budget_reached"] |= self.cable_budget - cable <= 1e-8
                flags["path_budget_reached"] |= self.path_budget - old_depth <= 1e-8
                return
            direction = unit(direction)
            if np.linalg.norm(direction) == 0:
                return
            xyz.append(xyz[p] + remaining * direction)
            parent.append(p)
            cable += remaining
            active.append((len(xyz) - 1, direction, old_depth + remaining))

        for direction in self.roots[root_index][1]:
            extend(0, direction, 0.0)
        while active:
            i, incoming, depth = active.popleft()
            if cable >= self.cable_budget - 1e-8:
                flags["cable_budget_reached"] = True
                break
            if len(xyz) >= CONFIG["max_nodes"]:
                flags["node_cap_reached"] = True
                break
            if depth >= self.path_budget - 1e-8:
                flags["path_budget_reached"] = True
                continue
            event_index = self._choose(self.tree, state(soma, xyz[i], incoming, depth), rng)
            _, train_incoming, outgoing = self.events[event_index]
            directions = rotate_match(outgoing, train_incoming, incoming)
            for direction in directions:
                extend(i, direction, depth)
        arbor = Arbor(np.asarray(xyz), np.asarray(parent, dtype=int), np.r_[1, np.full(len(xyz) - 1, 2, dtype=int)]).validate()
        return arbor, {**flags, "truncated": bool(any(flags.values())), "budget_um": self.cable_budget,
                       "path_budget_um": self.path_budget, "training_template": None}


def geometry_metrics(prediction, truth):
    pw, tw = prediction.weights(), truth.weights()
    pd = cKDTree(truth.xyz).query(prediction.xyz)[0]
    td = cKDTree(prediction.xyz).query(truth.xyz)[0]
    forward = float(np.average(pd, weights=pw))
    backward = float(np.average(td, weights=tw))
    tol = CONFIG["f1_tolerance_um"]
    precision = float(np.average(pd <= tol, weights=pw))
    recall = float(np.average(td <= tol, weights=tw))
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    true_cable = truth.cable()
    return {
        "distance_um": (forward + backward) / 2,
        "prediction_to_truth_um": forward, "truth_to_prediction_um": backward,
        "precision_200": precision, "recall_200": recall, "f1_200": f1,
        "cable_relative_error": abs(prediction.cable() - true_cable) / max(true_cable, 1e-12),
        "cable_signed_relative_error": (prediction.cable() - true_cable) / max(true_cable, 1e-12),
        "branch_abs_error": float(abs(prediction.branches() - truth.branches())),
    }


def mean_metrics(rows):
    return {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}


def training_medoid(training):
    centred = [t["arbor"].shifted(np.zeros(3)) for t in training]
    count = Counter(t["animal"] for t in training)
    weights = np.array([1 / count[t["animal"]] for t in training], dtype=float)
    weights /= weights.sum()
    losses = []
    for p in centred:
        losses.append(sum(w * geometry_metrics(p, q)["distance_um"] for w, q in zip(weights, centred)))
    return int(np.argmin(losses))


def predict_query(training, soma, query_id, fold):
    """Only training records and soma are accepted. All outputs are independent of truth."""
    assert not any(t["neuron"] == query_id or t["animal"] == fold for t in training)
    model = GrowthModel(training, stable_seed("fit", fold))
    predictions = []
    for sample in range(CONFIG["samples"]):
        seed = stable_seed("growth", query_id, sample)
        a, info = model.predict(soma, seed)
        predictions.append(("growth", sample, seed, a, info))
    nearest = int(np.argmin([np.linalg.norm(t["arbor"].xyz[0] - soma) for t in training]))
    medoid = training_medoid(training)
    for name, choice in (("nearest", nearest), ("medoid", medoid)):
        a = training[choice]["arbor"].shifted(soma)
        predictions.append((name, 0, None, a, {"truncated": False, "node_cap_reached": False,
                            "budget_um": None, "training_template": training[choice]["neuron"]}))
    a = Arbor(np.asarray(soma, dtype=float).reshape(1, 3), np.array([-1]), np.array([1]))
    predictions.append(("soma", 0, None, a, {"truncated": False, "node_cap_reached": False,
                        "budget_um": None, "training_template": None}))
    return predictions


def current_declaration():
    data = json.loads(VALIDATION.read_text())
    choices = {}
    for row in data["neurons"]:
        if int(row["type_counts"].get("2", 0)) <= 0:
            continue
        neuron = row["neuron"]
        if neuron not in choices or row["frame"] == "portal_true_hemisphere":
            choices[neuron] = row
    queries = []
    for neuron, d in sorted(choices.items()):
        path = RAW / d["file"]
        assert sha(path) == d["sha256"], "raw input hash differs from validated record"
        reflected = d["frame"] == "portal_true_hemisphere" and d["index"]["hemisphere"] == "Right"
        soma = np.array(d["soma"], dtype=float)
        if reflected:
            soma[2] = 11400.0 - soma[2]
        queries.append({"neuron": neuron, "animal": neuron.split("_")[0], "tag": d["index"]["type_name"],
                        "hemisphere": d["index"]["hemisphere"], "soma_um": soma.tolist(),
                        "file": str(path.relative_to(ROOT)), "sha256": sha(path),
                        "source_frame": d["frame"], "reflect_to_left": reflected})
    animals = sorted({q["animal"] for q in queries})
    assert len(queries) == 16 and len(animals) == 14, "predeclared primary cohort changed"
    folds = [{"animal": a, "train": [q["neuron"] for q in queries if q["animal"] != a],
              "test": [q["neuron"] for q in queries if q["animal"] == a]} for a in animals]
    return {"status": "declared_unscored", "config": CONFIG, "queries": queries, "folds": folds,
            "code_sha256": sha(SCRIPT), "protocol_sha256": sha(PROTOCOL),
            "validation_sha256": sha(VALIDATION), "display_case_ids": [q["neuron"] for q in queries[:3]]}


def prepare():
    d = current_declaration()
    path = BASE / "declaration.json"
    if path.exists() and json.loads(path.read_text()) != d:
        raise RuntimeError("existing declaration differs; preserve it and explicitly remove it before redeclaring")
    write_json(path, d)
    print(f"Declared {len(d['queries'])} neurons / {len(d['folds'])} animals, no geometry scored: {path}")
    return d


def synthetic_training():
    out = []
    for animal in range(4):
        soma = np.array([5000.0 + 100 * animal, 5500.0, 4500.0])
        # Two independent neurons from animal 0 tests animal balancing.
        for n in range(2 if animal == 0 else 1):
            xyz = np.array([[0, 0, 0], [50, 0, 0], [100, 0, 0], [150, 25, 0],
                            [150, -25, 0], [200, 50, 0], [200, -50, 0]], dtype=float) + soma
            a = resample(Arbor(xyz, np.array([-1, 0, 1, 2, 2, 3, 4]), np.array([1, 2, 2, 2, 2, 2, 2])))
            out.append({"neuron": f"A{animal}_{n}", "animal": f"A{animal}", "arbor": a})
    return out


def synthetic_checks():
    training = synthetic_training()
    model = GrowthModel(training, 123)
    query = np.array([5450.0, 5500.0, 4500.0])
    a, info = model.predict(query, 987)
    b, _ = model.predict(query.copy(), 987)
    assert np.array_equal(a.xyz, b.xyz) and np.array_equal(a.parent, b.parent)
    assert np.array_equal(a.xyz[0], query) and len(a.xyz) <= CONFIG["max_nodes"]
    assert a.cable() <= model.cable_budget + 1e-6
    assert a.path_lengths().max() <= model.path_budget + 1e-6
    assert geometry_metrics(a, a)["distance_um"] == 0
    assert geometry_metrics(a, a)["f1_200"] == 1
    # A far displacement must be detected, not erased by centering in the scorer.
    shifted = Arbor(a.xyz + [3000, 0, 0], a.parent, a.types)
    assert geometry_metrics(a, shifted)["distance_um"] > 2000
    assert geometry_metrics(a, shifted)["f1_200"] == 0
    # Branch/tip-preserving resampling retains topology and the exact endpoints.
    source = training[0]["arbor"]
    sampled = resample(source, 20)
    assert source.branches() == sampled.branches()
    for i, c in enumerate(source.children()):
        if len(c) != 1:
            assert cKDTree(sampled.xyz).query(source.xyz[i])[0] < 1e-8
    # Transport preserves lengths, pairwise fork angles and incoming direction.
    v = rotate_match(np.eye(3), [1, 0, 0], [0, 1, 0])
    assert np.allclose(v @ v.T, np.eye(3)) and np.allclose(v[0], [0, 1, 0])
    # Query and training animal overlap must fail before fitting.
    try:
        predict_query(training, query, "NEW", "A0")
    except AssertionError:
        pass
    else:
        raise AssertionError("animal leakage guard failed")
    # Prediction has no access path to a query truth file or truth argument.
    import inspect
    assert tuple(inspect.signature(GrowthModel.predict).parameters) == ("self", "soma", "seed")
    return {"passed": True, "synthetic_nodes": len(a.xyz), "synthetic_branches": a.branches(),
            "invariants": ["repeatable_seed", "soma_only_query", "animal_disjoint_guard", "fold_only_budgets",
                           "topology_preserving_resampling", "absolute_coordinate_scoring", "distance_identity",
                           "fork_rotation", "finite_connected_tree"]}


def check():
    result = {**synthetic_checks(), "code_sha256": sha(SCRIPT), "protocol_sha256": sha(PROTOCOL)}
    write_json(BASE / "check.json", result)
    print(canonical(result))


def freeze(review_path):
    declaration = prepare()
    review_path = Path(review_path).resolve()
    review = review_path.read_text()
    for token in ("D110_REVIEW: PASS", declaration["code_sha256"], declaration["protocol_sha256"]):
        if token not in review:
            raise RuntimeError(f"independent review lacks required token: {token}")
    check_path = BASE / "check.json"
    checked = json.loads(check_path.read_text())
    assert checked["passed"] and checked["code_sha256"] == sha(SCRIPT) and checked["protocol_sha256"] == sha(PROTOCOL)
    frozen = {"status": "frozen", "declaration_sha256": sha(BASE / "declaration.json"),
              "review_path": str(review_path), "review_sha256": sha(review_path),
              "check_sha256": sha(check_path), "code_sha256": sha(SCRIPT), "protocol_sha256": sha(PROTOCOL)}
    write_json(BASE / "frozen.json", frozen)
    print(f"Frozen without scoring: {BASE / 'frozen.json'}")


def verify_frozen(frozen_path):
    frozen_path = Path(frozen_path).resolve()
    f = json.loads(frozen_path.read_text())
    assert f["status"] == "frozen"
    assert f["code_sha256"] == sha(SCRIPT) and f["protocol_sha256"] == sha(PROTOCOL)
    assert f["review_sha256"] == sha(f["review_path"])
    assert f["check_sha256"] == sha(BASE / "check.json")
    assert f["declaration_sha256"] == sha(BASE / "declaration.json")
    declared = json.loads((BASE / "declaration.json").read_text())
    assert declared == current_declaration(), "inputs, splits or configuration changed"
    return declared


def verify_prediction_bundle(manifest, declared, base=None):
    """Reject incomplete or altered outputs before the scorer reads any truth.

    The optional base exists for synthetic tamper tests; production always uses
    BASE. Every prediction is hashed here, before any query's truth is opened.
    """
    base = Path(base or BASE).resolve()
    manifest_path = base / "prediction_manifest.json"
    marker = (base / "PREDICTIONS_COMPLETE").read_text().strip()
    assert marker == sha(manifest_path), "prediction completion marker/hash mismatch"
    assert manifest_path.read_text() == canonical(manifest), "passed and saved prediction manifests differ"
    assert (base / "declaration.json").read_text() == canonical(declared), "declaration differs from saved canonical version"
    assert manifest["declaration_sha256"] == sha(base / "declaration.json"), "prediction declaration hash mismatch"
    assert manifest["frozen_manifest_sha256"] == sha(base / "frozen.json"), "frozen manifest hash mismatch"
    assert manifest["code_sha256"] == declared["code_sha256"]
    assert manifest["protocol_sha256"] == declared["protocol_sha256"]
    assert manifest["config"] == declared["config"] == CONFIG
    expected = {q["neuron"]: q for q in declared["queries"]}
    assert len(expected) == len(declared["queries"]) == 16, "declared cohort must have16unique queries"
    actual_ids = [q["neuron"] for q in manifest["queries"]]
    assert len(actual_ids) == len(set(actual_ids)) == 16 and set(actual_ids) == set(expected), "query IDs are missing, duplicated or unexpected"
    combinations = {("growth", i) for i in range(CONFIG["samples"])} | {("nearest", 0), ("medoid", 0), ("soma", 0)}
    seen_paths = set()
    for query in manifest["queries"]:
        q = expected[query["neuron"]]
        for key in ("animal", "tag", "hemisphere", "soma_um"):
            assert query[key] == q[key], f"query metadata differs: {key}"
        assert query["fold"] == q["animal"]
        fold = next(f for f in declared["folds"] if f["animal"] == q["animal"])
        train_animals = sorted({expected[n]["animal"] for n in fold["train"]})
        assert query["train_animals"] == train_animals and q["animal"] not in train_animals
        keys = [(p["model"], p["sample"]) for p in query["predictions"]]
        assert len(keys) == len(combinations) and len(set(keys)) == len(keys) and set(keys) == combinations, "model/sample set is incomplete, duplicated or unexpected"
        for p in query["predictions"]:
            assert isinstance(p["model"], str) and type(p["sample"]) is int
            relative = f"predictions/{q['neuron']}/{p['model']}_{p['sample']}.npz"
            assert p["path"] == relative and relative not in seen_paths, "prediction path aliases another output"
            seen_paths.add(relative)
            path = (base / relative).resolve()
            assert path.is_relative_to(base) and path.is_file()
            assert sha(path) == p["sha256"], f"prediction bytes changed: {relative}"


def verify_truth_hash(record, root=None):
    """Run immediately before opening each withheld raw geometry for scoring."""
    path = Path(root or ROOT) / record["file"]
    assert sha(path) == record["sha256"], f"withheld raw input changed: {record['neuron']}"
    return path


def score_predictions(manifest):
    """Read truth only after validating the complete, frozen prediction bundle."""
    rows = []
    declared = json.loads((BASE / "declaration.json").read_text())
    verify_prediction_bundle(manifest, declared)
    records = {q["neuron"]: q for q in declared["queries"]}
    for query in manifest["queries"]:
        q = records[query["neuron"]]
        truth_path_raw = verify_truth_hash(q)
        truth = resample(extract_explicit_axon(truth_path_raw, q["reflect_to_left"]))
        truth_path = f"truth/{q['neuron']}.npz"
        save_arbor(BASE / truth_path, truth)
        models = {}
        for record in query["predictions"]:
            assert sha(BASE / record["path"]) == record["sha256"], "prediction bytes changed during scoring"
            metrics = geometry_metrics(load_arbor(BASE / record["path"]), truth)
            models.setdefault(record["model"], {"samples": []})["samples"].append(metrics)
        for values in models.values():
            values["mean"] = mean_metrics(values["samples"])
        rows.append({"neuron": q["neuron"], "animal": q["animal"], "tag": q["tag"],
                     "hemisphere": q["hemisphere"], "soma_um": q["soma_um"],
                     "truth_path": truth_path, "truth_sha256": sha(BASE / truth_path), "models": models})
    names = ("growth", "nearest", "medoid", "soma")
    animal_metrics = []
    for animal in sorted({r["animal"] for r in rows}):
        r = [q for q in rows if q["animal"] == animal]
        animal_metrics.append({"animal": animal, "neurons": len(r),
                               "models": {m: mean_metrics([q["models"][m]["mean"] for q in r]) for m in names}})
    summary = {m: mean_metrics([r["models"][m] for r in animal_metrics]) for m in names}
    delta = np.array([r["models"]["growth"]["distance_um"] - r["models"]["nearest"]["distance_um"] for r in animal_metrics])
    rng = np.random.default_rng(stable_seed("bootstrap"))
    boot = np.array([rng.choice(delta, size=len(delta), replace=True).mean() for _ in range(CONFIG["bootstrap_replicates"])])
    result = {"config": CONFIG, "queries": rows, "animal_metrics": animal_metrics, "summary": summary,
              "comparison": {"growth_minus_nearest_distance_um": float(delta.mean()),
                             "animal_differences": [{"animal": r["animal"], "distance_um": float(d)} for r, d in zip(animal_metrics, delta)],
                             "bootstrap95_um": np.quantile(boot, [0.025, 0.975]).tolist(),
                             "animals_growth_better": int((delta < 0).sum()), "animals": len(delta)},
              "limitations": ["Exploratory convenience sample:16neurons/14animals; not an independent confirmatory cohort.",
                              "Only explicit type2 axons; unlabelled type0 neurons excluded.",
                              "Soma-conditioned growth; no transcriptome or peptide input, no synaptic prediction.",
                              "Nearest-sample metrics use50um resampling; geometry below that scale is not tested.",
                              "Bootstrap is descriptive; overlapping CV training sets make its interval approximate.",
                              "No anatomical obstacle constraints; generated paths may be biologically implausible.",
                              "Five generated samples are all scored; no best-of-sample selection."]}
    write_json(BASE / "metrics.json", result)
    return result


def run(frozen_path):
    declared = verify_frozen(frozen_path)
    if (BASE / "RUN_STARTED").exists():
        raise RuntimeError("D-110 already started; preserve this run and use a new experiment for revisions")
    (BASE / "RUN_STARTED").write_text(sha(frozen_path) + "\n")
    records = {q["neuron"]: q for q in declared["queries"]}
    queries = []
    # Per-fold training loader is explicitly passed only train IDs. Query raw
    # geometry can be training data in another CV fold, never in its own fit.
    for fold in declared["folds"]:
        assert not set(fold["train"]) & set(fold["test"])
        training = []
        for neuron in fold["train"]:
            q = records[neuron]
            assert q["animal"] != fold["animal"]
            a = resample(extract_explicit_axon(ROOT / q["file"], q["reflect_to_left"]))
            training.append({"neuron": neuron, "animal": q["animal"], "arbor": a})
        for neuron in fold["test"]:
            q = records[neuron]
            soma = np.asarray(q["soma_um"], dtype=float)
            predicted = predict_query(training, soma, neuron, fold["animal"])
            out = []
            for model, sample, seed, arbor, info in predicted:
                path = f"predictions/{neuron}/{model}_{sample}.npz"
                save_arbor(BASE / path, arbor)
                out.append({"model": model, "sample": sample, "seed": seed, "path": path,
                            "sha256": sha(BASE / path), "nodes": len(arbor.xyz),
                            "cable_um": arbor.cable(), "branches": arbor.branches(), **info})
            queries.append({"neuron": neuron, "animal": q["animal"], "tag": q["tag"],
                            "hemisphere": q["hemisphere"], "soma_um": q["soma_um"],
                            "fold": fold["animal"], "train_animals": sorted({t["animal"] for t in training}),
                            "predictions": out})
            print(f"Predictions saved for {neuron} (held-out animal {fold['animal']})", flush=True)
    import platform
    import scipy
    manifest = {"config": CONFIG, "code_sha256": sha(SCRIPT), "protocol_sha256": sha(PROTOCOL),
                "declaration_sha256": sha(BASE / "declaration.json"),
                "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                            "scipy": scipy.__version__, "executable": sys.executable},
                "frozen_manifest_sha256": sha(frozen_path), "queries": queries}
    write_json(BASE / "prediction_manifest.json", manifest)
    (BASE / "PREDICTIONS_COMPLETE").write_text(sha(BASE / "prediction_manifest.json") + "\n")
    result = score_predictions(manifest)
    print(canonical({"summary": result["summary"], "comparison": result["comparison"]}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("prepare")
    sub.add_parser("check")
    p = sub.add_parser("freeze")
    p.add_argument("--review", required=True)
    p = sub.add_parser("run")
    p.add_argument("--frozen", required=True)
    args = parser.parse_args()
    if args.command == "prepare": prepare()
    elif args.command == "check": check()
    elif args.command == "freeze": freeze(args.review)
    elif args.command == "run": run(args.frozen)


if __name__ == "__main__":
    main()
