"""SYNARCH general reconstruction layer: modules and simple baselines (library, no biological I/O).

  * JointState / JointMixtureModule: the UNCHANGED JointModel of analysis/brm_joint.py wrapped as a
    module that answers any conditional among its fields. A fitted model is read out once into plain
    arrays; batch inference then reproduces JointModel.infer() exactly (tested), without its per-query cost.
    A JointState is read-only, and its fingerprint binds everything that is saved with it: parameters and the
    field each belongs to, field order and types, transforms, standardisation constants, representation,
    fit settings, SUPPORT RANGES and support margin, interval settings, training units and groups.
  * TaskBaselines: the simple comparators that are always reported (MEAN, IDENTITY, COPY and its
    same-label variant) and the two declared references (k-nearest-neighbour average, ridge). MEAN and
    IDENTITY have one declared setting: whether every brain sample or every neuron counts once.
  * ReferenceAtlasModule: population-level reference alignment; its answers are flagged 'reference'.
  * AxonArborAdapter: the interface through which axon reconstruction becomes one module. It declares
    'axon_arbor' with status 'blocked: independent validation pending' (nothing passed to it unblocks it),
    rejects a measured regional profile as input, and admits an inferred input, whatever its variable name,
    only after the provenance check of the core (hidden_axon_problem). It does not import the axon engine.

Nothing here opens a data file. Rows arrive as arrays from analysis/synarch_recon.py loaders.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
import types
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import cdist
from scipy.special import logsumexp, ndtr, ndtri


def _sibling(name):
    """Load a sibling file by path under a fixed name (python3 -I does not put this folder on sys.path)."""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


core = _sibling("synarch_recon")

LOG_2PI = math.log(2 * math.pi)
SUPPORT_MARGIN = 0.05      # a continuous input may lie this fraction of the training range outside it
QUANTILE_ITERATIONS = 60   # bisection steps of every mixture quantile (medians and interval end points)
STATE_FORMAT = "synarch_recon_state_v2"   # format of a saved JointState; bump it with any change of the reader
COMPARATORS = ("MEAN", "IDENTITY", "COPY", "COPY_SAME_LABEL", "KNN", "RIDGE")
ALWAYS_REPORTED = ("MEAN", "IDENTITY", "COPY", "COPY_SAME_LABEL")
EXTRA_REFERENCES = ("KNN", "RIDGE")
TUNED = ("MEAN", "IDENTITY", "KNN", "RIDGE")      # comparators with a setting chosen on inner folds
WEIGHTINGS = ("sample", "neuron")                 # MEAN / IDENTITY: every brain sample, or every neuron, counts once
DEFAULT_WEIGHTING = "sample"

# ----------------------------------------------------------------------------------------------
# fields and rows
# ----------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class FieldLayout:
    """Registry variables -> scalar model fields. A vector contributes one field per coordinate."""
    namespace: str
    variables: tuple
    fields: tuple          # (field, variable, coordinate or None, 'continuous' | 'categorical', transform)

    @classmethod
    def from_registry(cls, registry, namespace, variables, coordinates=None):
        coordinates = coordinates or {}
        fields = []
        for name in variables:
            spec = registry.variable(namespace, name)
            if spec.predictor_ban:
                raise core.ForbiddenInputError(f"{name} can never be a model field: {spec.predictor_ban}")
            if spec.kind == "categorical":
                fields.append((name, name, None, "categorical", "none"))
            elif spec.kind == "continuous":
                fields.append((name, name, None, "continuous", spec.transform))
            elif spec.kind == "vector":
                chosen = tuple(coordinates.get(name, spec.components))
                if not chosen or not set(chosen) <= set(spec.components):
                    raise core.RegistryError(f"{name}: coordinates must be registered coordinates")
                fields.extend((f"{name}:{c}", name, c, "continuous", spec.transform) for c in chosen)
            else:
                raise core.RegistryError(f"{name} of kind {spec.kind} cannot be a model field")
        return cls(namespace, tuple(variables), tuple(fields))

    def field_types(self):
        return {f: kind for f, _, _, kind, _ in self.fields}

    def fields_of(self, variable, coordinates=None):
        chosen = [(f, c) for f, v, c, _, _ in self.fields if v == variable]
        if coordinates is None:
            return [f for f, _ in chosen]
        by = {c: f for f, c in chosen}
        return [by[c] for c in coordinates]

    def coordinates_of(self, variable):
        return tuple(c for _, v, c, _, _ in self.fields if v == variable and c is not None)

    def info(self, field_name):
        for f, v, c, kind, transform in self.fields:
            if f == field_name:
                return v, c, kind, transform
        raise KeyError(field_name)


@dataclass
class Rows:
    """Units as arrays on the MODELLING scale: continuous fields are floats (NaN = missing),
    categorical fields are object arrays (None = missing)."""
    unit: np.ndarray
    group: np.ndarray
    cont: dict
    cat: dict

    def __len__(self):
        return len(self.unit)

    def take(self, index):
        index = np.asarray(index)
        return Rows(self.unit[index], self.group[index], {f: v[index] for f, v in self.cont.items()},
                    {f: v[index] for f, v in self.cat.items()})

    def matrix(self, fields):
        if not fields:
            return np.zeros((len(self), 0))
        return np.column_stack([self.cont[f] for f in fields])

    def frozen_copy(self):
        """A private copy whose arrays cannot be written to."""
        def frozen(array):
            out = np.array(array)
            out.setflags(write=False)
            return out

        return Rows(frozen(self.unit), frozen(self.group), {f: frozen(v) for f, v in self.cont.items()},
                    {f: frozen(v) for f, v in self.cat.items()})


def jiao_layout(registry, regions):
    return FieldLayout.from_registry(registry, "jiao", ("soma_position", "source_label", "regional_axon_length"),
                                     {"regional_axon_length": tuple(regions)})


def rows_from_table(table, layout):
    """The training table as model rows (lengths on the declared log10(1 + L) scale)."""
    cont, cat = {}, {}
    profile = None if table.length is None else table.profile()
    for name, variable, coordinate, _, _ in layout.fields:
        if variable == "soma_position":
            cont[name] = table.soma[:, core.SOMA_AXES.index(coordinate)].astype(float)
        elif variable == "source_label":
            cat[name] = table.label.astype(object)
        elif variable == "regional_axon_length":
            if profile is None:
                raise core.ReconError("the table was loaded without outcomes")
            cont[name] = profile[:, table.regions.index(coordinate)].astype(float)
        else:
            raise core.ReconError(f"no table column for {variable}")
    return Rows(table.neuron.copy(), table.sample.copy(), cont, cat)


REPRESENTATIONS = ("log", "hurdle")
ON, OFF = "on", "off"


def encoding_for(layout, representation):
    """How each continuous field enters the joint model.

    'log'    : one Gaussian field per coordinate, on the modelling scale (53 fields for the Jiao layout).
    'hurdle' : a non-negative, zero-inflated coordinate (modelling transform log10(1 + x)) becomes two
               fields of the SAME unchanged model: a categorical presence field ('on' when the value is
               above zero) and a Gaussian level field that is observed only when the value is above zero
               and missing otherwise. Zero is then an exact outcome with its own probability, instead of
               a value a Gaussian has to stretch over.
    """
    if representation not in REPRESENTATIONS:
        raise ValueError(f"unknown representation {representation}")
    return {name: ("hurdle" if representation == "hurdle" and transform == "log10_1p" else "gaussian")
            for name, _, _, kind, transform in layout.fields if kind == "continuous"}


def model_fields(name, how):
    return (name + "|on", name + "|level") if how == "hurdle" else (name,)


def model_field_types(layout, encoding):
    types = {}
    for name, _, _, kind, _ in layout.fields:
        if kind == "categorical":
            types[name] = "categorical"
        elif encoding[name] == "hurdle":
            types[name + "|on"], types[name + "|level"] = "categorical", "continuous"
        else:
            types[name] = "continuous"
    return types


def joint_rows(rows, layout, encoding=None):
    """Mappings for JointModel.fit: one per unit, observed fields only."""
    encoding = encoding or encoding_for(layout, "log")
    out = []
    for i in range(len(rows)):
        row = {}
        for name, _, _, kind, _ in layout.fields:
            if kind == "continuous":
                value = rows.cont.get(name)
                if value is None or not np.isfinite(value[i]):
                    continue
                if encoding[name] == "hurdle":
                    present = value[i] > 0
                    row[name + "|on"] = ON if present else OFF
                    if present:
                        row[name + "|level"] = float(value[i])
                else:
                    row[name] = float(value[i])
            else:
                value = rows.cat.get(name)
                if value is not None and value[i] is not None:
                    row[name] = str(value[i])
        out.append(row)
    return out


# ----------------------------------------------------------------------------------------------
# the joint model as a fitted state with exact batch conditionals
# ----------------------------------------------------------------------------------------------


def _mixture_quantiles(weights, means, scales, probabilities, chunk=256, iterations=QUANTILE_ITERATIONS):
    """Quantiles of univariate Gaussian mixtures. weights (n, K); means, scales (K, T) -> (n, T, Q).

    Same bracket and edge rules as JointModel._quantile (smallest and largest component quantile),
    solved by bisection for all rows and targets at once.
    """
    n, targets = weights.shape[0], means.shape[1]
    out = np.empty((n, targets, len(probabilities)))
    for q, p in enumerate(probabilities):
        individual = means + float(ndtri(p)) * scales
        lower, upper = individual.min(axis=0), individual.max(axis=0)
        for start in range(0, n, chunk):
            w = weights[start:start + chunk]

            def cdf(x, w=w):
                return np.einsum("ck,ckt->ct", w, ndtr((x[:, None, :] - means[None]) / scales[None]))

            lo = np.tile(lower, (len(w), 1))
            hi = np.tile(upper, (len(w), 1))
            at_lower = cdf(lo) - p >= 0
            at_upper = (cdf(hi) - p <= 0) & ~at_lower
            a, b = lo.copy(), hi.copy()
            for _ in range(iterations):
                mid = 0.5 * (a + b)
                below = cdf(mid) < p
                a = np.where(below, mid, a)
                b = np.where(below, b, mid)
            solved = 0.5 * (a + b)
            out[start:start + chunk, :, q] = np.where(at_lower, lo, np.where(at_upper, hi, solved))
    return out


def _hurdle_quantiles(weights, present, means, scales, probabilities, chunk=256, iterations=QUANTILE_ITERATIONS):
    """Quantiles of zero-inflated mixtures: per component a point mass at zero (1 - present) and a
    Gaussian (present). weights (n, K); present, means, scales (K, T) -> (n, T, Q).

    The quantile is the smallest x whose mixture CDF reaches the probability; where the mass at zero
    already reaches it, the quantile is exactly zero.
    """
    n, targets = weights.shape[0], means.shape[1]
    out = np.empty((n, targets, len(probabilities)))
    for q, p in enumerate(probabilities):
        individual = means + float(ndtri(p)) * scales
        lower = np.minimum(individual.min(axis=0), 0.0) - 1e-6
        upper = np.maximum(individual.max(axis=0), 0.0)
        for start in range(0, n, chunk):
            w = weights[start:start + chunk]

            def cdf(x, w=w):
                part = (1.0 - present)[None] * (x[:, None, :] >= 0) + present[None] * ndtr(
                    (x[:, None, :] - means[None]) / scales[None])
                return np.einsum("ck,ckt->ct", w, part)

            a, b = np.tile(lower, (len(w), 1)), np.tile(upper, (len(w), 1))
            for _ in range(iterations):
                mid = 0.5 * (a + b)
                below = cdf(mid) < p
                a = np.where(below, mid, a)
                b = np.where(below, b, mid)
            out[start:start + chunk, :, q] = np.where(np.abs(b) < 1e-9, 0.0, b)
    return out


_ARRAY_KINDS = ("weights", "means", "variances", "tables")
_BLOCK_KEYS = ("name", "fields", "n_rows", "n_components", "converged", "iterations", "objective", "log_likelihood")
_SETTING_KEYS = ("max_components", "seed", "n_init", "max_iter", "tol", "variance_floor", "dirichlet_alpha",
                 "representation")
_HEAD_KEYS = ("format", "schema", "core_sha256", "reader", "layout", "encoding", "field_meta", "settings", "units",
              "groups", "ranges", "index", "blocks")


def _plain(value):
    """A JSON-plain deep copy: mappings -> dict, sequences -> list, numpy scalars -> Python numbers."""
    if isinstance(value, Mapping):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def _frozen(value):
    """A read-only deep copy of a JSON-plain value: mappings -> mapping proxies, lists -> tuples."""
    if isinstance(value, Mapping):
        return types.MappingProxyType({k: _frozen(v) for k, v in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_frozen(v) for v in value)
    return value


def _read_only(array):
    out = np.array(array, dtype=float)
    out.setflags(write=False)
    return out


def _finite_float(value):
    return type(value) is float and math.isfinite(value)


class JointState:
    """A fitted JointModel read out into arrays, with the training provenance it must carry.

    Callers speak in LOGICAL fields (one per coordinate of a registry variable). `encoding` says how
    each continuous logical field was presented to the model: 'gaussian' (itself) or 'hurdle'
    (presence + level). All values go in and come out on the modelling scale of the logical field.

    A state is READ-ONLY once built (its arrays, mappings and attributes cannot be changed in place), and it
    is checked against the form a fitted state must have. Its fingerprint, state_hash(), is computed from
    exactly what save() writes and is recomputed on every call (it is never cached).

    WHAT THE FINGERPRINT BINDS (everything in a saved state that can change an answer, a refusal or a notice):
       1. the state format and the layer's schema name;
       2. the SHA-256 of the unchanged joint-model source that fitted it;
       3. the reader's settings: the SUPPORT MARGIN (how far outside the training range a continuous input
          may lie before it is refused or marginalised) and the bisection steps of every interval quantile;
       4. the layout: namespace, variables, and every field (name, variable, coordinate, kind, modelling
          transform) in order;
       5. the encoding of every continuous field (representation: 'gaussian' or presence / level);
       6. per model field: type, supported or not, block, the standardisation constants (training mean and
          scale), the categories and their order;
       7. the fit settings: number of components requested, seed, restarts, iteration limit, tolerance,
          variance floor, Dirichlet prior, representation;
       8. per block: name, fields in order, rows, number of components, convergence record;
       9. every parameter array (weights, means, variances, tables): values, shape, and WHICH block, kind and
          field it belongs to;
      10. the SUPPORT RANGES of the continuous fields (they decide refusals and marginalisation);
      11. the training units and the training groups, as lists (they decide the leakage refusal, the
          seen-group notice, the coverage flag and the query-group exclusion of a hidden-axon forecast).
    NOT in a state, and therefore not in its fingerprint: the policy and the interval level (run-time settings,
    returned with every result), the task definitions and the validation evidence (bound by the pre-run
    declaration and by run_complete.json), the comparator bank (it has its own hash, TaskBaselines).
    """

    def __init__(self, *, layout, encoding, field_meta, blocks, settings, units, groups, ranges, support_margin=None,
                 quantile_iterations=None, model=None):
        self.layout = layout
        self.encoding = _frozen({str(k): v for k, v in dict(encoding).items()})
        self.field_meta = _frozen(_plain(field_meta))
        self.settings = _frozen(_plain(settings))
        frozen = []
        for b in blocks:
            if set(b) != set(_BLOCK_KEYS) | set(_ARRAY_KINDS):
                raise core.ReconError("not a valid joint state: a block holds other entries than a fitted block has")
            entry = {k: _frozen(_plain(b[k])) for k in _BLOCK_KEYS}
            entry["weights"] = _read_only(b["weights"])
            for kind in ("means", "variances", "tables"):
                entry[kind] = types.MappingProxyType({str(n): _read_only(a) for n, a in dict(b[kind]).items()})
            frozen.append(types.MappingProxyType(entry))
        self.blocks = tuple(frozen)
        self.units = frozenset(str(u) for u in units)
        self.groups = tuple(sorted({str(g) for g in groups}))
        # continuous logical field -> (min, max) of the training values
        self.ranges = types.MappingProxyType({str(k): (tuple(v) if isinstance(v, (list, tuple)) else v)
                                              for k, v in dict(ranges).items()})
        self.support_margin = SUPPORT_MARGIN if support_margin is None else support_margin
        self.quantile_iterations = QUANTILE_ITERATIONS if quantile_iterations is None else quantile_iterations
        self.model = model                      # the JointModel object, when just fitted (never saved)
        self._validate()
        self._model_block = {f: i for i, b in enumerate(self.blocks) for f in b["fields"]}
        self._kind = {}
        for name, _, _, kind, _ in layout.fields:
            self._kind[name] = "categorical" if kind == "categorical" else self.encoding[name]
        self.units_sha256 = core.ids_sha(sorted(self.units))
        self.groups_sha256 = core.ids_sha(self.groups)
        self.units_folded = frozenset(core.fold_id(u) for u in self.units)
        self.groups_folded = frozenset(core.fold_id(g) for g in self.groups)
        self._sealed = True

    def __setattr__(self, name, value):
        if self.__dict__.get("_sealed"):
            raise AttributeError("a JointState is read-only once built; fit or load another state instead")
        object.__setattr__(self, name, value)

    def _validate(self):
        """The form every fitted state has. A loaded file that does not have it is not a state."""
        def need(condition, why):
            if not condition:
                raise core.ReconError(f"not a valid joint state: {why}")

        layout = self.layout
        need(isinstance(layout, FieldLayout) and isinstance(layout.namespace, str) and layout.namespace
             and isinstance(layout.variables, tuple) and isinstance(layout.fields, tuple), "layout")
        need(all(isinstance(v, str) and v for v in layout.variables)
             and len(set(layout.variables)) == len(layout.variables), "layout variables")
        for entry in layout.fields:
            need(isinstance(entry, tuple) and len(entry) == 5, "layout field entry")
            name, variable, coordinate, kind, transform = entry
            need(isinstance(name, str) and name and variable in layout.variables
                 and (coordinate is None or isinstance(coordinate, str))
                 and kind in ("continuous", "categorical") and transform in core.TRANSFORMS, f"layout field {name!r}")
        names = [entry[0] for entry in layout.fields]
        need(names and len(set(names)) == len(names), "field names")
        continuous = {entry[0] for entry in layout.fields if entry[3] == "continuous"}
        need(set(self.settings) == set(_SETTING_KEYS) and self.settings["representation"] in REPRESENTATIONS, "settings")
        need(dict(self.encoding) == encoding_for(layout, self.settings["representation"]),
             "the encoding is not the one of the stated representation")
        expected = model_field_types(layout, self.encoding)
        need(set(self.field_meta) == set(expected), "field metadata of other fields than the model has")
        supported = set()
        for f, kind in expected.items():
            meta = self.field_meta[f]
            need(isinstance(meta, Mapping) and meta.get("type") == kind and isinstance(meta.get("supported"), bool),
                 f"metadata of {f}")
            if not meta["supported"]:
                need(set(meta) == {"type", "supported", "block"} and meta["block"] is None, f"metadata of {f}")
                continue
            supported.add(f)
            need(isinstance(meta["block"], str), f"block of {f}")
            if kind == "continuous":
                need(set(meta) == {"type", "supported", "block", "mean", "scale"} and _finite_float(meta["mean"])
                     and _finite_float(meta["scale"]) and meta["scale"] > 0, f"standardisation of {f}")
            else:
                categories = meta.get("categories")
                need(set(meta) == {"type", "supported", "block", "categories"} and isinstance(categories, tuple)
                     and categories and all(isinstance(c, str) for c in categories)
                     and len(set(categories)) == len(categories), f"categories of {f}")
        block_names = [b["name"] for b in self.blocks]
        need(block_names and all(isinstance(n, str) for n in block_names)
             and len(set(block_names)) == len(block_names), "block names")
        in_blocks = []
        for b in self.blocks:
            k, fields = b["n_components"], b["fields"]
            need(type(k) is int and k >= 1 and type(b["n_rows"]) is int and b["n_rows"] >= 1
                 and type(b["iterations"]) is int and type(b["converged"]) is bool
                 and _finite_float(b["objective"]) and _finite_float(b["log_likelihood"]), f"record of {b['name']}")
            need(isinstance(fields, tuple) and len(set(fields)) == len(fields)
                 and all(f in supported and self.field_meta[f]["block"] == b["name"] for f in fields),
                 f"fields of {b['name']}")
            cont = {f for f in fields if expected[f] == "continuous"}
            cat = set(fields) - cont
            w = b["weights"]
            need(w.shape == (k,) and np.isfinite(w).all() and (w > 0).all() and abs(float(w.sum()) - 1.0) < 1e-6,
                 f"weights of {b['name']}")
            need(set(b["means"]) == cont and set(b["variances"]) == cont and set(b["tables"]) == cat,
                 f"the parameters of {b['name']} are not those of its fields")
            for f in cont:
                need(b["means"][f].shape == (k,) and np.isfinite(b["means"][f]).all()
                     and b["variances"][f].shape == (k,) and np.isfinite(b["variances"][f]).all()
                     and (b["variances"][f] > 0).all(), f"parameters of {f}")
            for f in cat:
                table = b["tables"][f]
                need(table.shape == (k, len(self.field_meta[f]["categories"])) and np.isfinite(table).all()
                     and (table > 0).all(), f"table of {f}")
            in_blocks.extend(fields)
        need(sorted(in_blocks) == sorted(supported), "every supported field belongs to exactly one block")
        observed = {name for name in continuous
                    if self.field_meta[name + "|on" if self.encoding[name] == "hurdle" else name]["supported"]}
        need(set(self.ranges) == observed,
             "the support ranges are not exactly those of the continuous fields that were observed in training")
        for name, r in self.ranges.items():
            need(isinstance(r, tuple) and len(r) == 2 and all(_finite_float(v) for v in r) and r[0] <= r[1],
                 f"the support range of {name} is not two finite numbers in order")
        need(self.units and all(core.identifier_problem(u) is None for u in self.units), "training units")
        need(self.groups and all(core.identifier_problem(g) is None for g in self.groups), "training groups")
        need(_finite_float(self.support_margin) and self.support_margin >= 0, "support margin")
        need(type(self.quantile_iterations) is int and 1 <= self.quantile_iterations <= 200, "quantile iterations")

    # -- construction --------------------------------------------------------------------------
    @classmethod
    def from_model(cls, model, layout, encoding, rows, settings):
        meta = {}
        for name, m in model.field_metadata_.items():
            entry = {"type": m["type"], "supported": bool(m["supported"]), "block": m["block"]}
            if m["supported"] and m["type"] == "continuous":
                entry.update(mean=float(m["training_mean"]), scale=float(m["training_scale"]))
            elif m["supported"]:
                if any(t != "str" for t in m["category_types"]):
                    raise NotImplementedError("the batch reader handles string categories only")
                entry["categories"] = [str(c) for c in m["categories"]]
            meta[name] = entry
        blocks = []
        for b in model.blocks_:
            blocks.append({"name": b.name, "fields": list(b.fields), "weights": np.array(b.params.weights, float),
                           "means": {f: np.array(v, float) for f, v in b.params.means.items()},
                           "variances": {f: np.array(v, float) for f, v in b.params.variances.items()},
                           "tables": {f: np.array(v, float) for f, v in b.params.tables.items()},
                           "n_rows": int(b.n_rows), "n_components": int(b.n_components),
                           "converged": bool(b.converged), "iterations": int(b.iterations),
                           "objective": float(b.objective), "log_likelihood": float(b.log_likelihood)})
        ranges = {}
        for name, _, _, kind, _ in layout.fields:
            if kind == "continuous" and name in rows.cont:
                values = rows.cont[name][np.isfinite(rows.cont[name])]
                if len(values):
                    ranges[name] = (float(values.min()), float(values.max()))
        return cls(layout=layout, encoding=encoding, field_meta=meta, blocks=blocks, settings=dict(settings),
                   units=rows.unit, groups=rows.group, ranges=ranges, model=model)

    def diagnostics(self):
        return [{"name": b["name"], "fields": len(b["fields"]), "n_rows": b["n_rows"],
                 "n_components": b["n_components"], "converged": b["converged"], "iterations": b["iterations"],
                 "objective": b["objective"], "log_likelihood": b["log_likelihood"]} for b in self.blocks]

    # -- fingerprint, save, load ---------------------------------------------------------------
    def _saved(self):
        """Everything that is saved, in its saved form: (head without the hash, arrays in index order).
        The fingerprint covers exactly this, so no saved item can change without changing it."""
        index, arrays = [], []
        for bi, b in enumerate(self.blocks):
            index.append([bi, "weights", None])
            arrays.append(b["weights"])
            for kind in ("means", "variances", "tables"):
                for name in sorted(b[kind]):
                    index.append([bi, kind, name])
                    arrays.append(b[kind][name])
        head = {"format": STATE_FORMAT, "schema": core.SCHEMA, "core_sha256": core.BRM_JOINT_SHA256,
                "reader": {"support_margin": self.support_margin, "quantile_iterations": self.quantile_iterations},
                "layout": {"namespace": self.layout.namespace, "variables": list(self.layout.variables),
                           "fields": [list(f) for f in self.layout.fields]},
                "encoding": self.encoding, "field_meta": self.field_meta, "settings": self.settings,
                "units": sorted(self.units), "groups": list(self.groups),       # the lists themselves, not a digest
                "ranges": {k: list(v) for k, v in sorted(self.ranges.items())}, "index": index,
                "blocks": [{k: b[k] for k in _BLOCK_KEYS} for b in self.blocks]}
        return head, arrays

    @staticmethod
    def _fingerprint(head, arrays):
        return core.sha_bytes((core.canonical(head) + core.array_sha(*arrays)).encode())

    def state_hash(self):
        """SHA-256 of everything that is saved (see the class docstring for the list). Never cached."""
        return self._fingerprint(*self._saved())

    def save(self, path):
        head, arrays = self._saved()
        head["state_hash"] = self._fingerprint(head, arrays)
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(path, head=np.array(core.canonical(head)), **{f"a{i}": a for i, a in enumerate(arrays)})
        return path

    @classmethod
    def load(cls, path, expected_hash):
        """Read a saved state; refuse anything that is not exactly the EXPECTED state as save() writes it.

        `expected_hash` is required: the fingerprint that was recorded elsewhere when the state was made (the
        fit record, selection.json, serving.json). The hash stored inside the file only shows that the file
        agrees with itself, and anybody who edits the file can rewrite it; so a state is never loaded on the
        file's own word. A file that holds another state, however consistent, is refused."""
        def bad(why):
            raise core.FrozenSourceChanged(f"{path}: {why}")

        def no_constant(token):
            raise ValueError(f"non-finite number {token}")

        try:
            with np.load(path, allow_pickle=False) as data:
                if "head" not in data.files:
                    bad("no head")
                text = str(data["head"])
                stored = {k: data[k] for k in data.files if k != "head"}
            head = json.loads(text, parse_constant=no_constant)
        except core.FrozenSourceChanged:
            raise
        except Exception as error:      # noqa: BLE001  (an unreadable file is a changed file)
            bad(f"unreadable state file ({type(error).__name__}: {error})")
        if not isinstance(head, dict) or core.canonical(head) != text:
            bad("the stored head is not in canonical form")
        recorded = head.pop("state_hash", None)
        if not isinstance(recorded, str) or set(head) != set(_HEAD_KEYS):
            bad("the stored head does not hold exactly the entries of a saved state")
        if not isinstance(expected_hash, str) or recorded != expected_hash:
            bad(f"this file holds state {recorded}, not the expected state {expected_hash}")
        index = head["index"]
        if not isinstance(index, list) or set(stored) != {f"a{i}" for i in range(len(index))}:
            bad("the stored arrays are not the ones the index lists")
        arrays = [stored[f"a{i}"] for i in range(len(index))]
        if any(a.dtype != np.dtype("<f8") for a in arrays):
            bad("a stored array is not 64-bit floating point")
        if cls._fingerprint(head, arrays) != recorded:          # 1. the file as written
            bad("stored state does not reproduce its recorded hash")
        try:
            blocks = [dict(b, weights=None, means={}, variances={}, tables={}) for b in head["blocks"]]
            seen = set()
            for i, entry in enumerate(index):
                bi, kind, name = entry
                if type(bi) is not int or not 0 <= bi < len(blocks) or kind not in _ARRAY_KINDS \
                        or (kind == "weights") != (name is None) or (bi, kind, name) in seen:
                    bad("malformed array index")
                seen.add((bi, kind, name))
                if kind == "weights":
                    blocks[bi]["weights"] = arrays[i]
                else:
                    blocks[bi][kind][name] = arrays[i]
            if any(b["weights"] is None for b in blocks):
                bad("a block has no weights")
            layout = FieldLayout(head["layout"]["namespace"], tuple(head["layout"]["variables"]),
                                 tuple(tuple(f) for f in head["layout"]["fields"]))
            state = cls(layout=layout, encoding=head["encoding"], field_meta=head["field_meta"], blocks=blocks,
                        settings=head["settings"], units=head["units"], groups=head["groups"], ranges=head["ranges"],
                        support_margin=head["reader"]["support_margin"],
                        quantile_iterations=head["reader"]["quantile_iterations"])
        except core.FrozenSourceChanged:
            raise
        except (core.ReconError, AttributeError, KeyError, IndexError, TypeError, ValueError) as error:
            bad(f"not a state as save() writes it ({type(error).__name__}: {error})")
        # 2. the rebuilt state, under the running format, schema and joint-model hash, in normal form
        again, _ = state._saved()
        if core.canonical(again) != core.canonical(head) or state.state_hash() != recorded:
            bad("stored state does not reproduce its recorded hash under the running code")
        return state

    # -- logical fields ------------------------------------------------------------------------
    def _lead(self, name):
        """The model field that says whether a logical field was ever observed, and in which block."""
        return name + "|on" if self._kind[name] == "hurdle" else name

    def kind_of(self, name):
        """'categorical', 'gaussian' or 'hurdle' for a logical field of this state; None for another name."""
        return self._kind.get(name)

    def supported(self, name):
        return name in self._kind and bool(self.field_meta.get(self._lead(name), {}).get("supported"))

    def block_of(self, name):
        return self._model_block.get(self._lead(name)) if name in self._kind else None

    def categories(self, name):
        return list(self.field_meta[name]["categories"])

    def _presence(self, block, name):
        """Per component: probability that a hurdle field is above zero, and its level mean / sd."""
        k = len(block["weights"])
        on_meta, level_meta = self.field_meta[name + "|on"], self.field_meta[name + "|level"]
        present = np.zeros(k)
        if ON in on_meta["categories"]:
            present = block["tables"][name + "|on"][:, on_meta["categories"].index(ON)]
        if level_meta["supported"]:
            mu = level_meta["mean"] + level_meta["scale"] * block["means"][name + "|level"]
            sd = level_meta["scale"] * np.sqrt(block["variances"][name + "|level"])
        else:                                   # never above zero in training: only the mass at zero remains
            present, mu, sd = np.zeros(k), np.zeros(k), np.ones(k)
        return present, mu, sd

    # -- inference -----------------------------------------------------------------------------
    def posterior(self, block_index, rows, evidence):
        """Component probabilities of every row given its observed evidence fields in one block."""
        block = self.blocks[block_index]
        n = len(rows)
        joint = np.tile(np.log(block["weights"]), (n, 1))
        used, unknown = [], {}
        for name, _, _, kind, _ in self.layout.fields:
            if name not in evidence or self.block_of(name) != block_index or not self.supported(name):
                continue
            how = self._kind[name]
            if how == "categorical":
                if name not in rows.cat:
                    continue
                position = {c: i for i, c in enumerate(self.field_meta[name]["categories"])}
                index = np.array([-2 if v is None else position.get(str(v), -1) for v in rows.cat[name]])
                unknown[name] = index == -1
                mask = index >= 0
                if mask.any():
                    joint[mask] += np.log(block["tables"][name][:, index[mask]]).T
            else:
                if name not in rows.cont:
                    continue
                x = np.asarray(rows.cont[name], float)
                finite = np.isfinite(x)
                if not finite.any():
                    continue
                if how == "gaussian":
                    meta = self.field_meta[name]
                    z = (x[finite] - meta["mean"]) / meta["scale"]
                    mu, var = block["means"][name], block["variances"][name]
                    joint[finite] += -0.5 * (LOG_2PI + np.log(var) + (z[:, None] - mu) ** 2 / var)
                else:
                    on_meta, level_meta = self.field_meta[name + "|on"], self.field_meta[name + "|level"]
                    table = block["tables"][name + "|on"]
                    above, zero = finite & (x > 0), finite & (x <= 0)
                    if OFF in on_meta["categories"]:
                        joint[zero] += np.log(table[:, on_meta["categories"].index(OFF)])
                    if ON in on_meta["categories"]:
                        joint[above] += np.log(table[:, on_meta["categories"].index(ON)])
                        if level_meta["supported"] and above.any():
                            z = (x[above] - level_meta["mean"]) / level_meta["scale"]
                            mu, var = block["means"][name + "|level"], block["variances"][name + "|level"]
                            joint[above] += -0.5 * (LOG_2PI + np.log(var) + (z[:, None] - mu) ** 2 / var)
                    else:                           # a value above zero where training never saw one
                        unknown[name] = above
            used.append(name)
        shifted = joint - joint.max(axis=1, keepdims=True)
        return np.exp(shifted - logsumexp(shifted, axis=1, keepdims=True)), used, unknown

    def conditional(self, rows, evidence, *, cont_targets=(), cat_targets=(), quantiles=()):
        """Exact fitted conditionals for every row: means, mixture quantiles and class probabilities.

        `evidence` lists the logical fields that may condition; a field that is missing for a row, an
        unseen category, or a field of another co-observation block never conditions (as in
        JointModel.infer). Values are returned on the modelling scale of the fields.
        """
        evidence = set(evidence)
        targets = list(cont_targets) + list(cat_targets)
        if evidence & set(targets):
            raise ValueError("a target field cannot also be evidence")
        for name in targets:
            if not self.supported(name):
                raise core.ReconError(f"{name} has no training observations in this state")
        n = len(rows)
        mean = np.full((n, len(cont_targets)), np.nan)
        variance = np.full((n, len(cont_targets)), np.nan)
        quant = np.full((n, len(cont_targets), len(quantiles)), np.nan)
        proba, unknown, conditioned = {}, {}, {}
        for bi, block in enumerate(self.blocks):
            mine = [j for j, f in enumerate(cont_targets) if self.block_of(f) == bi]
            mine_cat = [f for f in cat_targets if self.block_of(f) == bi]
            if not mine and not mine_cat:
                continue
            weights, used, unk = self.posterior(bi, rows, evidence)
            conditioned[block["name"]] = used
            unknown.update(unk)
            plain = [j for j in mine if self._kind[cont_targets[j]] == "gaussian"]
            hurdle = [j for j in mine if self._kind[cont_targets[j]] == "hurdle"]
            if plain:
                names = [cont_targets[j] for j in plain]
                mu = np.column_stack([block["means"][f] for f in names])           # (K, T), standardised
                var = np.column_stack([block["variances"][f] for f in names])
                offset = np.array([self.field_meta[f]["mean"] for f in names])
                scale = np.array([self.field_meta[f]["scale"] for f in names])
                m = weights @ mu
                mean[:, plain] = offset + scale * m
                variance[:, plain] = scale ** 2 * np.maximum(weights @ (var + mu ** 2) - m ** 2, 0.0)
                if len(quantiles):
                    q = _mixture_quantiles(weights, mu, np.sqrt(var), list(quantiles),
                                           iterations=self.quantile_iterations)
                    quant[:, plain, :] = offset[None, :, None] + scale[None, :, None] * q
            if hurdle:
                parts = [self._presence(block, cont_targets[j]) for j in hurdle]
                present = np.column_stack([p[0] for p in parts])                    # (K, T)
                mu = np.column_stack([p[1] for p in parts])                         # modelling scale
                sd = np.column_stack([p[2] for p in parts])
                m = weights @ (present * mu)
                mean[:, hurdle] = m
                variance[:, hurdle] = np.maximum(weights @ (present * (sd ** 2 + mu ** 2)) - m ** 2, 0.0)
                if len(quantiles):
                    quant[:, hurdle, :] = _hurdle_quantiles(weights, present, mu, sd, list(quantiles),
                                                            iterations=self.quantile_iterations)
            for name in mine_cat:
                proba[name] = (weights @ block["tables"][name], self.categories(name))
        return {"mean": mean, "variance": variance, "quantiles": quant, "proba": proba, "unknown": unknown,
                "conditioned": conditioned}

    def outside_range(self, name, value, margin=None):
        """True when a continuous input lies outside the training range by more than the support margin
        OF THIS STATE (a fraction of the range)."""
        if name not in self.ranges:
            return False
        low, high = self.ranges[name]
        slack = (self.support_margin if margin is None else margin) * (high - low)
        return bool(value < low - slack or value > high + slack)


def fit_joint(rows, layout, *, components, seed, budget, representation="log"):
    """Fit the unchanged JointModel on the given rows and read it out. Returns (state, diagnostics)."""
    encoding = encoding_for(layout, representation)
    constructor = {"max_components": int(components), "seed": int(seed), "n_init": int(budget["n_init"]),
                   "max_iter": int(budget["max_iter"]), "tol": float(budget["tol"]),
                   "variance_floor": float(budget.get("variance_floor", 1e-3)),
                   "dirichlet_alpha": float(budget.get("dirichlet_alpha", 1.1))}
    model = core.JOINT.JointModel(model_field_types(layout, encoding), **constructor)
    model.fit(joint_rows(rows, layout, encoding))
    state = JointState.from_model(model, layout, encoding, rows, dict(constructor, representation=representation))
    return state, state.diagnostics()


# ----------------------------------------------------------------------------------------------
# tasks
# ----------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class TaskSpec:
    task: str
    title: str
    kind: str                       # 'profile' | 'position' | 'label'
    target_variable: str
    target_coordinates: tuple       # () for a categorical target
    target_fields: tuple
    evidence_cont: tuple            # continuous evidence fields, in the column order used for distances
    evidence_cat: tuple             # categorical evidence fields
    requires: tuple                 # ((variable, coordinates or None), ...)
    primary: str
    higher_is_better: bool
    copy_scaling: str               # 'raw' | 'z'
    knn_scalings: tuple
    ridge_degrees: tuple
    declared_measured_axon_input: bool
    clip_min: float | None          # physical support of the target on the modelling scale
    scale: str
    description: str

    def evidence_fields(self):
        return tuple(self.evidence_cont) + tuple(self.evidence_cat)

    def conditional(self):
        return core.Conditional(self.task + "_" + self.title, self.target_variable,
                                self.target_coordinates or None, self.requires, self.task,
                                self.declared_measured_axon_input, self.description)


PRIMARY = {"profile": ("mse", False), "position": ("euclid_mean", False), "label": ("accuracy", True)}


def make_task(layout, registry, *, task, title, target, target_coordinates=None, evidence, copy_scaling="z",
              knn_scalings=("z",), ridge_degrees=(1,), declared_measured_axon_input=False, description=""):
    """A task = one conditional: `target` (variable, optional coordinates) given `evidence` variables."""
    spec = registry.variable(layout.namespace, target)
    if spec.kind == "categorical":
        kind, coordinates, clip, scale = "label", (), None, "nominal"
    else:
        coordinates = tuple(target_coordinates) if target_coordinates is not None else layout.coordinates_of(target)
        kind = "profile" if spec.transform == "log10_1p" else "position"
        clip = 0.0 if spec.transform == "log10_1p" else None
        scale = "log10(1 + value)" if spec.transform == "log10_1p" else spec.units
    cont, cat, requires = [], [], []
    for variable, coords in evidence:
        e = registry.variable(layout.namespace, variable)
        if e.kind == "categorical":
            cat.extend(layout.fields_of(variable))
            requires.append((variable, None))
        else:
            chosen = tuple(coords) if coords is not None else layout.coordinates_of(variable)
            cont.extend(layout.fields_of(variable, chosen if e.kind == "vector" else None))
            whole = e.kind != "vector" or set(chosen) == set(layout.coordinates_of(variable))
            requires.append((variable, None if whole else chosen))
    primary, higher = PRIMARY[kind]
    return TaskSpec(task, title, kind, target, coordinates,
                    tuple(layout.fields_of(target, coordinates if spec.kind == "vector" else None)), tuple(cont),
                    tuple(cat), tuple(requires), primary, higher, copy_scaling, tuple(knn_scalings),
                    tuple(ridge_degrees), declared_measured_axon_input, clip, scale, description)


def build_tasks(layout, registry, t4_observed):
    """The four declared tasks of the masked-variable benchmark."""
    regions = layout.coordinates_of("regional_axon_length")
    observed = tuple(r for r in regions if r in set(t4_observed))
    hidden = tuple(r for r in regions if r not in set(t4_observed))
    if not observed or not hidden:
        raise core.ReconError("T4 needs a non-empty observed block and a non-empty hidden block")
    return {
        "T1": make_task(layout, registry, task="T1", title="profile_from_soma_label", target="regional_axon_length",
                        evidence=(("soma_position", None), ("source_label", None)), copy_scaling="raw",
                        knn_scalings=("raw", "z"), ridge_degrees=(1, 2),
                        description="regional profile (all regions) from soma position and source label; the "
                                    "neuron's axon is entirely hidden"),
        "T2": make_task(layout, registry, task="T2", title="soma_from_profile_label", target="soma_position",
                        evidence=(("regional_axon_length", None), ("source_label", None)), copy_scaling="raw",
                        knn_scalings=("raw", "z"), declared_measured_axon_input=True,
                        description="soma position from the MEASURED regional profile and the source label "
                                    "(a declared profile-input task; not an axon forecast)"),
        "T3": make_task(layout, registry, task="T3", title="label_from_soma_profile", target="source_label",
                        evidence=(("soma_position", None), ("regional_axon_length", None)), copy_scaling="z",
                        knn_scalings=("z",), declared_measured_axon_input=True,
                        description="source label from soma position and the MEASURED regional profile "
                                    "(a declared profile-input task)"),
        "T4": make_task(layout, registry, task="T4", title="outside_block_from_inside_block_soma_label",
                        target="regional_axon_length", target_coordinates=hidden,
                        evidence=(("regional_axon_length", observed), ("soma_position", None), ("source_label", None)),
                        copy_scaling="z", knn_scalings=("z",), declared_measured_axon_input=True,
                        description="the declared block of regions outside the hypothalamus from the MEASURED "
                                    "hypothalamic regions, soma position and source label (completion of a partly "
                                    "measured profile; not a hidden-axon forecast)"),
    }


def engine_predict(state, task, rows, *, level=None, median=False):
    """The joint engine's answer to a task for every row.

    Continuous targets: 'point' is the conditional mean (clipped to the target's physical support),
    'median' the conditional median, 'lower'/'upper' the central interval at `level`.
    Categorical targets: 'proba' over 'classes'.
    """
    evidence = task.evidence_fields()
    if task.kind == "label":
        out = state.conditional(rows, evidence, cat_targets=task.target_fields)
        proba, classes = out["proba"][task.target_fields[0]]
        return {"proba": proba, "classes": classes, "unknown": out["unknown"], "conditioned": out["conditioned"]}
    probabilities = ([0.5] if median else []) + ([(1 - level) / 2, (1 + level) / 2] if level is not None else [])
    out = state.conditional(rows, evidence, cont_targets=task.target_fields, quantiles=probabilities)

    def support(x):
        return x if task.clip_min is None else np.maximum(x, task.clip_min)

    result = {"point": support(out["mean"]), "unknown": out["unknown"], "conditioned": out["conditioned"]}
    position = 0
    if median:
        result["median"] = support(out["quantiles"][:, :, 0])
        position = 1
    if level is not None:
        result["lower"] = support(out["quantiles"][:, :, position])
        result["upper"] = support(out["quantiles"][:, :, position + 1])
    return result


# ----------------------------------------------------------------------------------------------
# simple baselines and references
# ----------------------------------------------------------------------------------------------


def sample_weighted_mean(values, group):
    """Mean of per-brain-sample means: every brain sample has equal weight."""
    values, group = np.asarray(values, float), np.asarray(group)
    return np.mean([values[group == g].mean(axis=0) for g in sorted(set(group.tolist()))], axis=0)


def _softmax(scores):
    scores = scores - scores.max(axis=1, keepdims=True)
    e = np.exp(scores)
    return e / e.sum(axis=1, keepdims=True)


class TaskBaselines:
    """Every simple comparator of one task, fitted on ONE training set.

    MEAN            training mean (label target: label frequencies). One declared setting, chosen by the caller
                    on inner folds like k and lambda: weighting 'sample' (every brain sample counts once) or
                    'neuron' (every neuron counts once; this is what a one-component mixture predicts).
    IDENTITY        the source label's training mean with the same weighting choice; global fallback (MEAN with
                    that weighting) for an unseen label.
    COPY            the value of the nearest training unit of ANOTHER brain sample in the observed space.
    COPY_SAME_LABEL the same, among training units with the query's source label (fallback: COPY).
    KNN             average over the k nearest training units of other brain samples (k, pool, scaling chosen
                    by the caller on inner folds).
    RIDGE           ridge regression on the standardised observed variables and the one-hot label
                    (label target: multinomial logistic regression with a ridge penalty).
    Distances are Euclidean, on the modelling scale ('raw') or on training-standardised coordinates ('z').
    Ties go to the first training unit in the stored order.
    """

    def __init__(self, task, train, *, vote_alpha=1.0, logistic_max_iter=500):
        self.task, self.n = task, len(train)
        self.group = np.asarray(train.group)
        self.unit = np.asarray(train.unit)
        self.X = train.matrix(task.evidence_cont)
        self.label = (np.array([str(v) for v in train.cat[task.evidence_cat[0]]]) if task.evidence_cat else None)
        self.vote_alpha, self.logistic_max_iter = float(vote_alpha), int(logistic_max_iter)
        self.mu = self.X.mean(axis=0) if self.X.shape[1] else np.zeros(0)
        sd = self.X.std(axis=0) if self.X.shape[1] else np.zeros(0)
        self.sd = np.where(sd > 0, sd, 1.0)
        samples = sorted(set(self.group.tolist()))
        if task.kind == "label":
            truth = np.asarray(train.cat[task.target_fields[0]], dtype=object)
            self.classes = sorted({str(t) for t in truth})
            position = {c: i for i, c in enumerate(self.classes)}
            self.Y = np.zeros((self.n, len(self.classes)))
            self.Y[np.arange(self.n), [position[str(t)] for t in truth]] = 1.0
        else:
            self.classes = None
            self.Y = train.matrix(task.target_fields)
        self.groups = tuple(str(g) for g in samples)
        self.groups_sha256 = core.ids_sha(self.groups)
        # 'sample': every brain sample counts once; 'neuron': every training neuron counts once
        self.mean_by = {"sample": np.mean([self.Y[self.group == s].mean(axis=0) for s in samples], axis=0),
                        "neuron": self.Y.mean(axis=0)}
        self.mean = self.mean_by["sample"]      # the prior of the label votes (COPY, KNN) and the KNN fallback
        self.label_mean_by, self.label_rows = {"sample": {}, "neuron": {}}, {}
        if self.label is not None:
            for value in sorted(set(self.label.tolist())):
                rows_of = self.label == value
                self.label_rows[value] = np.flatnonzero(rows_of)
                self.label_mean_by["sample"][value] = np.mean(
                    [self.Y[rows_of & (self.group == s)].mean(axis=0)
                     for s in samples if (rows_of & (self.group == s)).any()], axis=0)
                self.label_mean_by["neuron"][value] = self.Y[rows_of].mean(axis=0)
        self.label_mean = self.label_mean_by["sample"]
        self._ridge_cache, self._logistic_cache = {}, {}

    @staticmethod
    def weighting(params):
        value = (params or {}).get("weighting", DEFAULT_WEIGHTING)
        if value not in WEIGHTINGS:
            raise ValueError(f"unknown weighting {value}")
        return value

    def state_hash(self):
        """SHA-256 of the bank: the task it answers, its settings and every training row it holds."""
        t = self.task
        return core.digest({"schema": core.SCHEMA, "task": t.task, "kind": t.kind, "target_fields": t.target_fields,
                            "evidence_cont": t.evidence_cont, "evidence_cat": t.evidence_cat,
                            "copy_scaling": t.copy_scaling, "clip_min": t.clip_min, "vote_alpha": self.vote_alpha,
                            "logistic_max_iter": self.logistic_max_iter,
                            "rows": core.array_sha(self.unit, self.group, self.X, self.Y,
                                                   self.label if self.label is not None else np.zeros(0))})

    def answer_hash(self, name, params=None):
        """SHA-256 of ONE comparator's answer rule: the bank, the comparator and every setting it used.
        Two answers of different comparators, or of one comparator with different settings, never share it."""
        if name not in COMPARATORS:
            raise ValueError(f"unknown comparator {name}")
        if name in ("MEAN", "IDENTITY"):
            params = {"weighting": self.weighting(params)}
        elif name in ("COPY", "COPY_SAME_LABEL"):
            params = {"scaling": self.task.copy_scaling, "k": 1}
        return core.digest({"bank": self.state_hash(), "comparator": name, "params": params})

    def applicable(self, name):
        if name not in COMPARATORS:
            return False
        if name in ("IDENTITY", "COPY_SAME_LABEL"):
            return self.label is not None
        if name in ("COPY", "COPY_SAME_LABEL", "KNN", "RIDGE"):
            return self.X.shape[1] > 0 or (name == "RIDGE" and self.label is not None)
        return True

    # -- distances and neighbours ----------------------------------------------------------------
    def _scaled(self, X, scaling):
        if scaling == "raw":
            return X
        if scaling == "z":
            return (X - self.mu) / self.sd
        raise ValueError(f"unknown scaling {scaling}")

    def _distances(self, query, scaling, pool):
        D = cdist(self._scaled(query.matrix(self.task.evidence_cont), scaling), self._scaled(self.X, scaling),
                  "sqeuclidean")
        D[np.asarray(query.group)[:, None] == self.group[None, :]] = np.inf      # another brain sample only
        if pool == "same_label":
            qlabel = np.array([str(v) for v in query.cat[self.task.evidence_cat[0]]])
            restricted = np.where(qlabel[:, None] == self.label[None, :], D, np.inf)
            empty = ~np.isfinite(restricted).any(axis=1)
            restricted[empty] = D[empty]                                            # fallback: every label
            D = restricted
        elif pool != "all":
            raise ValueError(f"unknown pool {pool}")
        return D

    def _neighbours(self, query, scaling, pool, k_max):
        D = self._distances(query, scaling, pool)
        k_max = min(int(k_max), self.n)
        order = np.argsort(D, axis=1, kind="stable")[:, :k_max]
        valid = np.isfinite(np.take_along_axis(D, order, axis=1))
        return order, valid

    def _knn(self, order, valid, k):
        """Average of the k nearest (label target: one-hot votes smoothed by the sample-weighted prior)."""
        k = min(int(k), order.shape[1])
        idx, ok = order[:, :k], valid[:, :k]
        values = self.Y[idx] * ok[:, :, None]
        count = ok.sum(axis=1)
        if self.task.kind == "label":
            return (values.sum(axis=1) + self.vote_alpha * self.mean) / (count[:, None] + self.vote_alpha)
        out = values.sum(axis=1) / np.maximum(count, 1)[:, None]
        out[count == 0] = self.mean
        return out

    # -- ridge and logistic ------------------------------------------------------------------------
    def _design(self, X, labels, degree):
        Z = (X - self.mu) / self.sd
        parts = [Z]
        if degree == 2:
            p = Z.shape[1]
            parts.append(Z ** 2)
            parts.extend(Z[:, [i]] * Z[:, [j]] for i in range(p) for j in range(i + 1, p))
        elif degree != 1:
            raise ValueError("ridge degree is 1 or 2")
        if self.label is not None:
            values = np.array(sorted(self.label_mean))
            parts.append((np.array([str(v) for v in labels])[:, None] == values[None, :]).astype(float))
        return np.column_stack(parts)

    def _ridge(self, degree):
        if degree not in self._ridge_cache:
            F = self._design(self.X, self.label, degree)
            centre = F.mean(axis=0)
            spread = F.std(axis=0)
            spread = np.where(spread > 0, spread, 1.0)
            if self.label is not None:                      # one-hot columns stay 0/1
                spread[-len(self.label_mean):] = 1.0
            G = (F - centre) / spread
            target_mean = self.Y.mean(axis=0)
            evals, evecs = np.linalg.eigh(G.T @ G)
            self._ridge_cache[degree] = (centre, spread, target_mean, np.maximum(evals, 0.0), evecs,
                                         evecs.T @ (G.T @ (self.Y - target_mean)))
        return self._ridge_cache[degree]

    def _ridge_predict(self, query, lam, degree):
        centre, spread, target_mean, evals, evecs, projected = self._ridge(degree)
        coefficients = evecs @ (projected / (evals + lam * self.n)[:, None])
        labels = query.cat[self.task.evidence_cat[0]] if self.label is not None else None
        G = (self._design(query.matrix(self.task.evidence_cont), labels, degree) - centre) / spread
        return target_mean + G @ coefficients

    def _logistic(self, lam):
        if lam not in self._logistic_cache:
            Z = (self.X - self.mu) / self.sd
            n, p = Z.shape
            c = self.Y.shape[1]

            def objective(theta):
                W, b = theta[:p * c].reshape(p, c), theta[p * c:]
                P = _softmax(Z @ W + b)
                loss = -np.log(np.maximum((P * self.Y).sum(axis=1), 1e-300)).mean() + lam * (W ** 2).sum()
                G = (P - self.Y) / n
                return loss, np.concatenate([(Z.T @ G + 2 * lam * W).ravel(), G.sum(axis=0)])

            fit = minimize(objective, np.zeros(p * c + c), jac=True, method="L-BFGS-B",
                           options={"maxiter": self.logistic_max_iter})
            self._logistic_cache[lam] = (fit.x[:p * c].reshape(p, c), fit.x[p * c:], bool(fit.success), int(fit.nit))
        return self._logistic_cache[lam]

    def logistic_report(self):
        """Did the optimiser of each fitted logistic reference stop by itself, and after how many iterations."""
        return {str(lam): {"converged": ok, "iterations": nit} for lam, (_, _, ok, nit) in self._logistic_cache.items()}

    # -- public --------------------------------------------------------------------------------------
    def predict_many(self, name, query, grid=(None,)):
        """Predictions for every parameter setting in `grid` (shared distance work is done once)."""
        if not self.applicable(name):
            raise core.ReconError(f"{name} does not apply to task {self.task.task}")
        m = len(query)
        out = []
        if name == "MEAN":
            return [np.tile(self.mean_by[self.weighting(params)], (m, 1)) for params in grid]
        if name == "IDENTITY":
            labels = query.cat[self.task.evidence_cat[0]]
            return [np.vstack([self.label_mean_by[self.weighting(params)].get(str(v), self.mean_by[self.weighting(params)])
                               for v in labels]) for params in grid]
        if name in ("COPY", "COPY_SAME_LABEL"):
            order, valid = self._neighbours(query, self.task.copy_scaling,
                                            "same_label" if name == "COPY_SAME_LABEL" else "all", 1)
            prediction = self._knn(order, valid, 1)
            return [prediction for _ in grid]
        if name == "KNN":
            cache = {}
            for params in grid:
                key = (params["scaling"], params["pool"])
                if key not in cache:
                    cache[key] = self._neighbours(query, params["scaling"], params["pool"],
                                                  max(p["k"] for p in grid if (p["scaling"], p["pool"]) == key))
                out.append(self._knn(*cache[key], params["k"]))
            return out
        for params in grid:                                       # RIDGE
            if self.task.kind == "label":
                W, b, _, _ = self._logistic(float(params["lambda"]))
                Z = (query.matrix(self.task.evidence_cont) - self.mu) / self.sd
                out.append(_softmax(Z @ W + b))
            else:
                out.append(self._ridge_predict(query, float(params["lambda"]), int(params.get("degree", 1))))
        return out

    def predict(self, name, query, params=None):
        return self.predict_many(name, query, (params,))[0]

    def interval(self, name, query, level):
        """Empirical central interval of the training values (MEAN: all rows; IDENTITY: the label's rows)."""
        if self.task.kind == "label" or name not in ("MEAN", "IDENTITY"):
            return None
        probabilities = [(1 - level) / 2, (1 + level) / 2]
        overall = np.quantile(self.Y, probabilities, axis=0)
        lower = np.tile(overall[0], (len(query), 1))
        upper = np.tile(overall[1], (len(query), 1))
        if name == "IDENTITY" and self.label is not None:
            labels = query.cat[self.task.evidence_cat[0]]
            cache = {}
            for i, value in enumerate(labels):
                value = str(value)
                if value in self.label_rows and len(self.label_rows[value]) >= 2:
                    if value not in cache:
                        cache[value] = np.quantile(self.Y[self.label_rows[value]], probabilities, axis=0)
                    lower[i], upper[i] = cache[value]
        return lower, upper


# ----------------------------------------------------------------------------------------------
# the joint mixture as a registered module
# ----------------------------------------------------------------------------------------------

def _prediction_set(classes, probabilities, level):
    order = np.argsort(-np.asarray(probabilities), kind="stable")
    chosen, total = [], 0.0
    for i in order:
        chosen.append(classes[i])
        total += float(probabilities[i])
        if total >= level:
            break
    return tuple(chosen)


class JointMixtureModule(core.Module):
    """The unchanged JointModel as a module: any conditional among its fields.

    Declared conditionals (tasks) carry their benchmark evidence. Any other conditional among the
    fields can be computed, but it is 'not_validated' and the declared policy decides what is returned.

    `train` are the rows the states were fitted on: the comparators that give baseline-level answers use
    exactly these rows, so the leakage guard and the training-group provenance of a state describe the
    comparator bank too. A module whose rows are not those of its states is refused at construction.
    """
    kind = "estimator"

    def __init__(self, registry, layout, *, states, train, tasks=None, evidence=None, name="joint_mixture",
                 default_state=None, vote_alpha=1.0, logistic_max_iter=500):
        self.registry, self.layout, self.name = registry, layout, name
        self.namespace = layout.namespace
        self.grain = registry.namespace(layout.namespace).grain
        self.tasks = dict(tasks or {})
        self.states = dict(states)
        self.default_state = default_state or next(iter(self.states))
        self.train = train.frozen_copy()            # the comparators' rows: a private copy that cannot be written to
        self._evidence = dict(evidence or {})
        self._conditionals = {t.conditional().name: t for t in self.tasks.values()}
        self.vote_alpha, self.logistic_max_iter = vote_alpha, logistic_max_iter
        self._baselines = {}
        train_units = frozenset(str(u) for u in train.unit)
        train_groups = tuple(sorted({str(g) for g in train.group}))
        self._registered, self._state_groups, self._bank_groups = {}, {}, {}
        for key, state in self.states.items():
            if state.layout != layout:
                raise core.ReconError(f"state {key} was fitted on another field layout than this module's")
            if state.units != train_units or state.groups != train_groups:
                raise core.ReconError(f"the comparator rows are not the rows state {key} was fitted on: the units or "
                                      "brain samples differ, so a baseline-level answer could not be guarded")
            self._registered[key] = state.state_hash()
            self._state_groups[self._registered[key]] = frozenset(state.groups)

    # -- registration --------------------------------------------------------------------------
    def produces(self):
        return tuple(self.layout.variables)

    def consumes(self):
        return tuple(self.layout.variables)

    def conditionals(self):
        return tuple(t.conditional() for t in self.tasks.values())

    def evidence_for(self, conditional_name):
        return self._evidence.get(conditional_name, core.ValidationEvidence())

    def set_evidence(self, conditional_name, evidence):
        if conditional_name not in self._conditionals:
            raise core.RegistryError(f"{conditional_name} is not a declared conditional of {self.name}")
        self._evidence[conditional_name] = evidence

    def state_groups(self, state_hash):
        """Training groups of one of this module's states, or of a comparator bank that gave a baseline-level
        answer in this process (by the hash that answer reported)."""
        found = self._state_groups.get(state_hash)
        return found if found is not None else self._bank_groups.get(state_hash)

    def _state_key(self, task):
        return task.task if task.task in self.states else self.default_state

    def state_for(self, task):
        return self.states[self._state_key(task)]

    # -- planning ------------------------------------------------------------------------------
    def _match(self, target, coordinates, supplied):
        for task in self.tasks.values():
            if task.target_variable != target.variable or set(task.target_coordinates) != set(coordinates):
                continue
            if set(v for v, _ in task.requires) != set(supplied):
                continue
            if all(set(supplied[v]) == set(c if c is not None else self.layout.coordinates_of(v))
                   for v, c in task.requires):
                return task
        return None

    def plan(self, request, target, evidence, policy):
        layout, registry = self.layout, self.registry
        spec = registry.variable(self.namespace, target.variable)
        usable, ignored, notices, values = {}, [], [], {}
        for variable, record in evidence.items():
            if variable not in layout.variables or record.status == "reference":
                ignored.append((variable, "not a field of this module"))
                continue
            if record.variable != variable:
                return core.Refusal("invalid_request", f"a record of {record.variable} was supplied as {variable}")
            v_spec = registry.variable(self.namespace, variable)
            # every record states its scale and frame; its value reaches the modelling scale by ONE conversion
            if v_spec.kind == "categorical":
                problem = core.statement_problem(v_spec, record.scale, record.transform, record.frame)
                if problem:
                    return core.Refusal("invalid_request", problem)
                values[variable] = str(record.value)
            else:
                try:
                    values[variable] = core.model_values(v_spec, record)
                except core.ScaleError as error:
                    return core.Refusal("invalid_request", str(error))
                if v_spec.kind == "vector":
                    known = set(layout.coordinates_of(variable))
                    if not record.coordinates or not set(record.coordinates) <= known \
                            or len(set(record.coordinates)) != len(record.coordinates) \
                            or len(values[variable]) != len(record.coordinates):
                        return core.Refusal("invalid_request", f"{variable}: coordinates outside this module's fields")
                elif len(values[variable]) != 1:
                    return core.Refusal("invalid_request", f"{variable}: one value expected")
            usable[variable] = record
        # An inferred input that was computed only from measurements this request also supplies adds no
        # information about the unit; conditioning on it would count the same evidence twice. It never conditions.
        complete = {v for v, r in usable.items() if r.status == "measured"
                    and (registry.variable(self.namespace, v).kind != "vector"
                         or set(r.coordinates) >= set(layout.coordinates_of(v)))}
        for variable in [v for v, r in usable.items() if r.status == "inferred"]:
            roots = core.lineage_roots(usable[variable])
            if all(status == "measured" and root in complete for root, status in roots):
                ignored.append((variable, "inferred only from measurements that this request supplies: it adds no "
                                          "information and does not condition the answer"))
                notices.append(core.Notice("redundant_inferred_input",
                                           f"the inferred {variable} was computed only from measurements supplied in "
                                           "this request; it is not used", variable, "ignored"))
                del usable[variable]
        statuses = {v: r.status for v, r in usable.items()}
        have = set(usable[target.variable].coordinates) if target.variable in usable else set()
        if spec.kind == "vector":
            coordinates = tuple(target.coordinates) if target.coordinates is not None else \
                tuple(c for c in layout.coordinates_of(target.variable) if c not in have)
            if not coordinates or not set(coordinates) <= set(layout.coordinates_of(target.variable)):
                return core.Refusal("invalid_request", f"{target.variable}: nothing left to infer, or coordinates "
                                    "outside this module's fields")
        else:
            coordinates = ()
        supplied = {v: (tuple(r.coordinates) if r.coordinates else ()) for v, r in usable.items()}
        task = self._match(target, coordinates, supplied)
        if task is None:
            adhoc = [(v, supplied[v] or None) for v in layout.variables if v in supplied]
            task_like = make_task(layout, registry, task="adhoc", title="undeclared", target=target.variable,
                                  target_coordinates=coordinates or None, evidence=adhoc)
        else:
            task_like = task
            # a declared task answers in ITS OWN coordinate order: names and values stay aligned whatever
            # the order in which the request listed the coordinates
            coordinates = tuple(task.target_coordinates)
        key = self._state_key(task_like)
        state = self.states[key]
        if state.state_hash() != self._registered[key]:
            raise core.FrozenSourceChanged(f"state {key} of module {self.name} is no longer the state that was "
                                           "registered: it was changed in memory")
        cont, cat = {}, {}
        for variable, record in usable.items():
            v_spec = registry.variable(self.namespace, variable)
            if v_spec.kind == "categorical":
                cat[layout.fields_of(variable)[0]] = values[variable]
                continue
            names = layout.fields_of(variable, record.coordinates) if v_spec.kind == "vector" else \
                layout.fields_of(variable)
            if record.status == "inferred" and any(state.kind_of(n) == "hurdle" for n in names):
                return core.Refusal("invalid_request",
                                    f"an inferred {variable} cannot be evidence under the presence / level encoding "
                                    "of this state: it is a mean over 'absent' and 'present', not an observation of "
                                    "either")
            for name, value in zip(names, values[variable]):
                cont[name] = float(value)
        # support: unseen categories and continuous inputs outside the training range
        for name in list(cat):
            if not state.supported(name) or cat[name] not in state.categories(name):
                variable = layout.info(name)[0]
                if policy.on_unseen_category == "refuse":
                    return core.Refusal("out_of_support", f"{variable} = '{cat[name]}' was never seen in training")
                notices.append(core.Notice("out_of_support", f"{variable} = '{cat[name]}' was never seen in "
                                           "training; the answer marginalises over it", variable, "marginalised"))
                ignored.append((variable, "a value never seen in training; marginalised"))
                del cat[name]
        outside = defaultdict(list)
        for name, value in cont.items():
            if state.outside_range(name, value):
                outside[layout.info(name)[0]].append(name)
        for variable, names in outside.items():
            detail = (f"{variable}: {len(names)} coordinate(s) outside the training range "
                      f"(margin {state.support_margin:.0%})")
            if policy.on_out_of_range == "refuse":
                return core.Refusal("out_of_support", detail)
            notices.append(core.Notice("out_of_support", detail + "; the variable is dropped", variable, "marginalised"))
            ignored.append((variable, "outside the training range; marginalised"))
            for name in layout.fields_of(variable):
                cont.pop(name, None)
        target_fields = task_like.target_fields
        for name in target_fields:
            if not state.supported(name):
                return core.Refusal("unidentifiable_cross_modal",
                                    f"{target.variable}: no training unit of this state measured it")
        evidence_fields = set(cont) | set(cat)
        blocks = {state.block_of(f) for f in target_fields}
        linked_fields = {f for f in evidence_fields if state.block_of(f) in blocks}
        if evidence_fields and not linked_fields:
            return core.Refusal("unidentifiable_cross_modal",
                                f"no training unit measured {target.variable} together with any supplied variable "
                                "(disconnected co-observation blocks): the conditional is not identified",
                                tuple(sorted(str(b) for b in blocks | {state.block_of(f) for f in evidence_fields})))
        # what conditions the answer, exactly: a variable is 'used' only through fields linked to the target;
        # everything supplied that does not condition is listed as ignored, with the reason, and never as used
        used = []
        for variable in layout.variables:
            if variable not in usable:
                continue
            mine = [f for f in layout.fields_of(variable) if f in evidence_fields]
            linked = [f for f in mine if f in linked_fields]
            if linked:
                used.append((variable, statuses[variable]))
            if mine and len(linked) < len(mine):
                unlinked = [layout.info(f)[1] for f in mine if f not in linked_fields]
                label = variable if not linked else f"{variable}[{', '.join(str(c) for c in unlinked)}]"
                ignored.append((label, "measured only with other units than the target; cannot condition"))
                notices.append(core.Notice("evidence_not_linked",
                                           f"{label} was never measured together with {target.variable} in this "
                                           "state, so it does not condition the answer", variable, "ignored"))
        implied = registry.implied_group(self.namespace, request.unit_id)
        in_state = bool(core.id_in(request.unit_id, state.units, state.units_folded))
        # the unit's brain sample as the TRAINING rows spell it (if it is among them): comparators that copy or
        # average donors exclude the query's own brain sample by this name, whatever the spelling of the request
        donor_group = None
        for candidate in (request.group.group_id, implied):
            if core.id_in(candidate, state.groups, state.groups_folded):
                donor_group = next(g for g in state.groups if core.fold_id(g) == core.fold_id(candidate))
                break
        return core.Plan(self, target, coordinates, None if task is None else task.conditional(), tuple(used),
                         tuple(dict.fromkeys(ignored)), notices, in_state, donor_group is not None,
                         {"task": task_like, "state": state, "state_key": key,
                          "cont": {f: v for f, v in cont.items() if f in linked_fields},
                          "cat": {f: v for f, v in cat.items() if f in linked_fields},
                          "group": donor_group if donor_group is not None else request.group.group_id,
                          "unit": request.unit_id, "declared": task is not None})

    def provenance(self, plan):
        state = plan.payload["state"]
        return {"state_hash": state.state_hash(), "training_groups_hash": state.groups_sha256,
                "training_groups": state.groups}

    # -- answering -----------------------------------------------------------------------------
    def _query_rows(self, plan):
        payload = plan.payload
        task = payload["task"]
        cont = {f: np.array([payload["cont"].get(f, np.nan)]) for f in task.evidence_cont}
        cat = {f: np.array([payload["cat"].get(f)], dtype=object) for f in task.evidence_cat}
        return Rows(np.array([payload["unit"]]), np.array([str(payload["group"])]), cont, cat)

    def _estimate(self, task, variable, coordinates, point, median=None, lower=None, upper=None, level=None,
                  proba=None, classes=None, note=""):
        spec = self.registry.variable(self.namespace, variable)
        if task.kind == "label":
            best = classes[int(np.argmax(proba))]
            return core.Estimate(variable, "categorical", (), core.scale_statement(spec, "identity"), best,
                                 probabilities={c: float(p) for c, p in zip(classes, proba)},
                                 prediction_set=_prediction_set(classes, proba, level), note=note,
                                 transform="identity", frame=spec.frame)
        interval = None if lower is None else {"level": level, "lower": tuple(float(v) for v in lower),
                                               "upper": tuple(float(v) for v in upper)}
        original = None
        if spec.transform == "log10_1p":
            original = {"units": spec.units, "note": "back-transformed quantiles are exact; a back-transformed "
                        "mean would NOT be an expected length, so none is given"}
            if median is not None:
                original["median"] = tuple(float(v) for v in core.inverse_transform("log10_1p", median))
            if lower is not None:
                original["lower"] = tuple(float(v) for v in core.inverse_transform("log10_1p", lower))
                original["upper"] = tuple(float(v) for v in core.inverse_transform("log10_1p", upper))
        # point, median and interval are on the MODELLING scale; the estimate says so by tag, text and frame
        tag = spec.transform if spec.transform in core.MODEL_SCALE_TEXT else "identity"
        return core.Estimate(variable, spec.kind, tuple(coordinates), core.scale_statement(spec, tag),
                             tuple(float(v) for v in point),
                             None if median is None else tuple(float(v) for v in median), interval,
                             original_units=original, note=note, transform=tag, frame=spec.frame)

    def estimate(self, plan, level):
        task, state = plan.payload["task"], plan.payload["state"]
        rows = self._query_rows(plan)
        out = engine_predict(state, task, rows, level=level, median=True)
        note = "conditional mean, median and central interval of the fitted mixture; fitted distribution only"
        if task.kind == "label":
            return self._estimate(task, plan.target.variable, (), None, level=level, proba=out["proba"][0],
                                  classes=out["classes"], note=note)
        return self._estimate(task, plan.target.variable, plan.coordinates, out["point"][0], out["median"][0],
                              out["lower"][0], out["upper"][0], level, note=note)

    def baseline(self, plan, level, evidence):
        """The baseline-level answer: the benchmark's best comparator for a declared conditional,
        otherwise the label mean (global fallback). The hash it reports binds the bank's rows, the comparator
        that answered and every setting of it."""
        task = plan.payload["task"]
        key = (task.task, task.target_fields, task.evidence_fields())
        if key not in self._baselines:
            self._baselines[key] = TaskBaselines(task, self.train, vote_alpha=self.vote_alpha,
                                                 logistic_max_iter=self.logistic_max_iter)
        bank, state = self._baselines[key], plan.payload["state"]
        if bank.groups != state.groups or frozenset(str(u) for u in bank.unit) != state.units:
            raise core.FrozenSourceChanged("the comparator rows are no longer the rows the state was fitted on")
        name, params = evidence.best_comparator, evidence.best_comparator_params
        if not plan.payload["declared"] or name is None or not bank.applicable(name):
            name, params = ("IDENTITY" if bank.applicable("IDENTITY") else "MEAN"), None
        rows = self._query_rows(plan)
        if name not in ("MEAN", "IDENTITY") and not np.isfinite(rows.matrix(task.evidence_cont)).all():
            name, params = ("IDENTITY" if bank.applicable("IDENTITY") else "MEAN"), None    # an input was dropped
        # an unseen or missing label falls back inside the comparators (global mean, any-label donors)
        prediction = bank.predict(name, rows, params)[0]
        if task.clip_min is not None:
            prediction = np.maximum(prediction, task.clip_min)
        note = f"baseline-level answer ({name}" + (f", {params}" if params else "") + "); not the learned engine"
        coverage = bool(((evidence.coverage or {}).get("baselines", {}).get(name, {}) or {}).get("validated", False))
        if task.kind == "label":
            estimate = self._estimate(task, plan.target.variable, (), None, level=level, proba=prediction,
                                      classes=bank.classes, note=note)
        else:
            interval = bank.interval(name, rows, level)
            lower, upper = (interval[0][0], interval[1][0]) if interval is not None else (None, None)
            if lower is not None and task.clip_min is not None:
                lower, upper = np.maximum(lower, task.clip_min), np.maximum(upper, task.clip_min)
            estimate = self._estimate(task, plan.target.variable, plan.coordinates, prediction, None, lower, upper,
                                      level if lower is not None else None, note=note)
        answer_hash = bank.answer_hash(name, params)
        self._bank_groups[answer_hash] = frozenset(bank.groups)
        return core.BaselineAnswer(estimate, f"baseline:{name}", answer_hash, coverage, bank.groups_sha256, bank.groups)


# ----------------------------------------------------------------------------------------------
# reference atlas alignment: population level only
# ----------------------------------------------------------------------------------------------


class ReferenceAtlasModule(core.Module):
    """A declared population-level alignment to a reference atlas.

    It answers with the reference value of the unit's GROUP (for example the atlas mean of the unit's
    assigned subclass), status 'reference'. It never returns a measured or unit-level inferred value,
    and it exists only where the data contract declares a population-level reference edge.
    """
    kind = "reference"

    def __init__(self, registry, store, *, namespace, from_namespace, key_variable, variable, name=None):
        rule = registry.reference_rule(from_namespace, namespace)
        if rule is None:
            raise core.PairingError(f"the data contract declares no population-level reference alignment from "
                                    f"{from_namespace} to {namespace}")
        registry.variable(from_namespace, variable)
        registry.variable(namespace, key_variable)
        self.registry, self.store, self.rule = registry, store, rule
        self.namespace, self.grain = namespace, registry.namespace(namespace).grain
        self.from_namespace, self.key_variable, self.variable = from_namespace, key_variable, variable
        self.name = name or f"reference:{from_namespace}->{namespace}:{variable}"

    def produces(self):
        return (self.variable,)

    def consumes(self):
        return (self.key_variable,)

    def reference_answer(self, request, target, evidence):
        if target.variable != self.variable:
            return None
        key = evidence.get(self.key_variable)
        if key is None:
            return core.PropertyResult(target.variable, (), "refused", module=self.name, refusal=core.Refusal(
                "unidentifiable_cross_modal",
                f"{self.variable} exists only as a population-level reference of {self.from_namespace}; the unit's "
                f"{self.key_variable} is needed to attach its group's reference value"))
        view = self.store.unit(self.namespace, "cell_type_population", f"{self.key_variable}={key.value}",
                               statuses=("reference",))
        record = view.get(self.variable)
        if record is None:
            return core.PropertyResult(target.variable, (), "refused", module=self.name, refusal=core.Refusal(
                "unidentifiable_cross_modal",
                f"no reference value of {self.variable} is attached to the group {self.key_variable}={key.value}"))
        estimate = core.Estimate(self.variable, "vector" if record.coordinates else "continuous",
                                 tuple(record.coordinates), record.scale, record.value,
                                 note="population-level reference of the unit's group; other cells and animals",
                                 transform=record.transform, frame=record.frame)
        notice = core.Notice("population_reference",
                             f"reference value of the group {self.key_variable}={key.value} from {self.from_namespace}; "
                             "not a measurement of this unit and not a unit-level inference", self.variable)
        return core.PropertyResult(target.variable, tuple(record.coordinates), "reference", estimate=estimate,
                                   module=self.name, validation="reference_only", reference=record.reference,
                                   evidence_used=((self.key_variable, key.status),), notices=(notice,))


# ----------------------------------------------------------------------------------------------
# axon reconstruction as one module: the adapter
# ----------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class AxonQuery:
    """What an axon backend may see of a query neuron. Nothing else is passed.

    Every value comes with the scale statement, value tag and frame of its record."""
    unit_id: str
    group_id: str | None
    inputs: dict            # variable -> {'value', 'scale', 'transform', 'frame'}, permitted measured inputs only
    inferred: dict          # variable -> the same plus 'state_hash', for inferred properties that passed the checks


class AxonArborAdapter(core.Module):
    """Integration interface for axon reconstruction.

    * declares 'axon_arbor' with status 'blocked: independent validation pending'. NOTHING passed to this
      class unblocks it: no backend, and no evidence object (this layer holds exploratory evidence only).
      Unblocking is a reviewed change of this file, made together with an independently validated backend;
    * permitted measured query inputs (reports/axon_dataset_plan.md, section 8): soma position, source label,
      soma label, hemisphere, reconstruction type (core.HIDDEN_AXON_INPUTS);
    * inferred inputs, WHATEVER their variable name: only values that the reconstruction layer itself issued
      for this unit, whose whole ancestry (transitively) consists of permitted measured inputs and of states
      that are registered and do not contain the query's brain sample (core.hidden_axon_problem);
    * never inputs: any measured axon-derived value of the unit (its regional profile, in whole or in part,
      or its trace), and anything that descends from one. A measured regional profile is rejected with
      ForbiddenInputError;
    * this file does not import the axon engine.

    Limit that no check here can remove: records are not signed. A value that somebody re-enters by hand as
    a 'measured' record is taken as a measurement. The guarantee covers measured records as the loaders make
    them and inferred records exactly as reconstruct() issued them.
    """
    kind = "adapter"
    STATUS = "blocked: independent validation pending"
    PERMITTED = core.HIDDEN_AXON_INPUTS

    def __init__(self, registry, *, namespace="jiao", backend=None, validation=None, name="axon_arbor_adapter",
                 ledger=None):
        self.registry, self.namespace, self.name = registry, namespace, name
        self.grain = registry.namespace(namespace).grain
        self.backend = backend                    # kept for a later, validated integration; never called while blocked
        self.offered_validation = validation      # recorded only: it does not change the status
        self._ledger, self._layer = ledger, None
        self._states = {}
        self._conditional = core.Conditional("axon_arbor_from_permitted_inputs", "axon_arbor", None,
                                             (("soma_position", None), ("source_label", None)), None, False,
                                             "3D axon of a neuron whose axon is entirely hidden")

    def produces(self):
        return ("axon_arbor",)

    def consumes(self):
        return tuple(v for v in self.PERMITTED if self.registry.has_variable(self.namespace, v))

    def conditionals(self):
        return (self._conditional,)

    def evidence_for(self, conditional_name):
        """Always blocked. A backend or an evidence object passed to the constructor changes nothing."""
        return core.ValidationEvidence(self.STATUS, note="no independently validated axon engine is attached; "
                                       "nothing passed to the adapter can unblock it")

    def bind(self, reconstructor):
        """Registered in a reconstruction layer: inferred inputs are verified against THAT layer's ledger and
        against the states of its registered modules."""
        if self._layer is not None and self._layer is not reconstructor:
            raise core.RegistryError(f"{self.name} is already registered in another reconstruction layer")
        self._layer = reconstructor
        if self._ledger is None:
            self._ledger = reconstructor.ledger

    def register_state(self, state, training_groups=None, evidence_variables=None):
        """Declare a regional state whose inferred outputs may be passed on.

        Pass the state itself (anything with state_hash() and groups, such as a JointState): its hash and its
        training groups are then read from the same object. The older form register_state(hash, groups,
        inputs) is still accepted. Either way a declaration is believed only as far as it agrees with what the
        layer recorded when it issued the inferred value: check_inputs refuses a value whose producing state
        had other training groups than the registered ones. An empty group list, and a second registration
        of the same hash with other groups or inputs, are refused."""
        if hasattr(state, "state_hash") and hasattr(state, "groups"):
            state_hash, groups = state.state_hash(), frozenset(str(g) for g in state.groups)
            if training_groups is not None and frozenset(str(g) for g in training_groups) != groups:
                raise core.RegistryError("the declared training groups are not the groups of the state")
        else:
            if training_groups is None:
                raise core.RegistryError("register a state object, or a state hash with its training groups")
            state_hash, groups = str(state), frozenset(str(g) for g in training_groups)
        if not groups:
            raise core.RegistryError("a state without training groups cannot be registered")
        entry = (groups, None if evidence_variables is None else tuple(evidence_variables))
        if state_hash in self._states and self._states[state_hash] != entry:
            raise core.RegistryError("this state hash is already registered with other training groups or inputs")
        self._states[state_hash] = entry
        return state_hash

    def _resolve(self, state_hash):
        if state_hash in self._states:
            return self._states[state_hash][0]
        return None if self._layer is None else self._layer.resolve_state(state_hash)

    def _inputs_of(self, state_hash):
        return self._states[state_hash][1] if state_hash in self._states else None

    def check_inputs(self, group_id, evidence, unit_id=None):
        """Raise ForbiddenInputError unless every supplied record is a permitted input of a hidden-axon
        forecast of ONE unit (core.hidden_axon_problem: the same rule the general layer applies).

        `evidence` maps variable names to records. They must all be records of one unit (the pairing rule),
        of `unit_id` when it is given, and each must state the registered scale and frame of its variable."""
        keys = {record.key() for record in evidence.values()}
        if len(keys) > 1:
            raise core.ForbiddenInputError("the supplied records belong to different units: "
                                           + ", ".join(sorted(str(k) for k in keys)))
        key = next(iter(keys)) if keys else None
        if unit_id is not None:
            wanted = (self.namespace, self.grain, str(unit_id))
            if key is not None and key != wanted:
                raise core.ForbiddenInputError(f"the supplied records belong to {key}, not to the query unit {wanted}")
            key = wanted
        if key is not None and key[:2] != (self.namespace, self.grain):
            raise core.ForbiddenInputError(f"the supplied records are not {self.namespace} {self.grain} records")
        problem = None if key is None else core.identifier_problem(key[2], "the query's unit ID")
        if problem is None:
            problem = core.hidden_axon_problem(self.registry, self.namespace, key, group_id, evidence,
                                               ledger=self._ledger, resolve=self._resolve,
                                               state_inputs=self._inputs_of)
        if problem is None:
            for variable, record in evidence.items():
                problem = core.statement_problem(self.registry.variable(self.namespace, variable), record.scale,
                                                 record.transform, record.frame)
                if problem:
                    break
        if problem:
            raise core.ForbiddenInputError(problem)
        return True

    def plan(self, request, target, evidence, policy):
        self.check_inputs(request.group.group_id, evidence, request.unit_id)
        missing = [v for v, _ in self._conditional.requires if v not in evidence]
        if missing:
            return core.Refusal("invalid_request", f"an axon forecast needs {', '.join(missing)}")
        return core.Plan(self, target, (), self._conditional,
                         tuple((v, r.status) for v, r in evidence.items()), (), [], False, False,
                         {"unit": request.unit_id, "group": request.group.group_id, "evidence": evidence})

    def state_groups(self, state_hash):
        """The training groups of a state that was registered here explicitly (register_state)."""
        entry = self._states.get(state_hash)
        return None if entry is None else entry[0]

    def provenance(self, plan):
        return {"state_hash": None, "training_groups_hash": None, "training_groups": None}

    def query_for(self, plan):
        """What a backend would be given: values with their scale statement, tag and frame."""
        def stated(record, **extra):
            return dict(value=record.value, scale=record.scale, transform=record.transform, frame=record.frame, **extra)

        evidence = plan.payload["evidence"]
        return AxonQuery(plan.payload["unit"], plan.payload["group"],
                         {v: stated(r) for v, r in evidence.items() if r.status == "measured"},
                         {v: stated(r, state_hash=r.uncertainty.state_hash) for v, r in evidence.items()
                          if r.status == "inferred"})

    def estimate(self, plan, level):
        """Never reached through reconstruct() (the gate refuses a blocked conditional first), and refused here
        too: while the adapter is blocked no backend is called."""
        raise core.ReconError(self.STATUS)
