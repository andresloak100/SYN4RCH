"""SYNARCH general reconstruction layer: core (library, no biological I/O on import).

Given the measurements that exist for a unit, estimate requested properties that were not measured,
say how uncertain the estimate is, and REFUSE when the data cannot identify the answer.

What this file holds
  * a variable registry (kind, units, scale, frame, modelling transform) seeded from the local data
    contract (data/interim/brm_data_contract.json): namespaces, grains and the joins it allows;
  * provenance-aware observation records and the PAIRING RULE: two records are paired only when
    namespace, grain and unit ID are identical; a cross-dataset link exists only as a declared
    population-level REFERENCE alignment, flagged 'reference', never 'measured';
  * sealed loaders for the Jiao et al. 2025 TRAINING partition (rows of validation and test brain
    samples are dropped unparsed, immediately after they are read);
  * whole-brain-sample splits by a salted hash, and metrics aggregated per brain sample;
  * the reconstruction-request API: reconstruct(request) returns, per requested property, an
    INFERRED value with provenance or a structured REFUSAL;
  * the provenance ledger and the ONE rule for inputs of a hidden-axon forecast (hidden_axon_problem):
    every inferred input, whatever its variable name, is checked against what this layer itself issued.

What the provenance checks can and cannot see: records are plain values and are not signed. An inferred
value is accepted in a hidden-axon forecast only if this layer issued exactly that record (ledger) and
every state in its ancestry is registered and excludes the query's group. A value that somebody re-enters
by hand as a 'measured' record cannot be told apart from a measurement; no check here can.

On import it verifies the SHA-256 of analysis/brm_joint.py (the unchanged joint mixture model)
against the constant below and against the pin in data/interim/brm_template_method/frozen.json,
then executes exactly the verified bytes. It reads no biological table on import.

Concrete modules and baselines: analysis/synarch_recon_modules.py.
Benchmark commands: analysis/synarch_recon_benchmark.py. Protocol: reports/synarch_recon_protocol.md.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
import types
import unicodedata
from collections import Counter, defaultdict
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "synarch_recon_v1"

# ----------------------------------------------------------------------------------------------
# hashing and JSON
# ----------------------------------------------------------------------------------------------


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def jsonable(value):
    if isinstance(value, Mapping):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        items = sorted(value) if isinstance(value, (set, frozenset)) else value
        return [jsonable(v) for v in items]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, np.floating):
        value = float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "__dataclass_fields__"):
        return {k: jsonable(getattr(value, k)) for k in value.__dataclass_fields__}
    return value


def canonical(value):
    """Canonical JSON text: sorted keys, compact, trailing newline (non-finite numbers become null)."""
    return json.dumps(jsonable(value), sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"


def digest(value):
    return sha_bytes(canonical(value).encode())


def ids_sha(ids):
    """SHA-256 of a list of identifiers in the given order.

    The encoding is injective (a JSON list of strings): two different lists never share a digest.
    Joining the identifiers by a separator would not be, because an identifier may contain it."""
    return sha_bytes(json.dumps([str(i) for i in ids], ensure_ascii=True, separators=(",", ":")).encode())


def array_sha(*arrays):
    """SHA-256 of a sequence of arrays. Every array is framed by its kind, shape and byte length, so
    values, shapes, order and the number of arrays are all bound (an injective encoding)."""
    h = hashlib.sha256()
    h.update(len(arrays).to_bytes(8, "big"))
    for a in arrays:
        a = np.asarray(a)
        if a.dtype.kind in "US" or a.dtype == object:
            kind = "text"
            body = json.dumps([None if v is None else str(v) for v in a.ravel().tolist()], ensure_ascii=True,
                              separators=(",", ":")).encode()
        else:
            kind = "f8"
            body = np.ascontiguousarray(a, dtype="<f8").tobytes()
        frame = json.dumps([kind, list(a.shape), len(body)], separators=(",", ":")).encode()
        h.update(len(frame).to_bytes(8, "big"))
        h.update(frame)
        h.update(body)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value, *, exclusive=False, pretty=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = (json.dumps(jsonable(value), sort_keys=True, indent=1, allow_nan=False) + "\n") if pretty else canonical(value)
    with path.open("x" if exclusive else "w") as handle:
        handle.write(text)
    return path


def stable_seed(*parts):
    """A 32-bit seed that depends only on the given parts (no global random state)."""
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:4], "big")


# ----------------------------------------------------------------------------------------------
# errors
# ----------------------------------------------------------------------------------------------


class ReconError(RuntimeError):
    pass


class FrozenSourceChanged(ReconError):
    """A pinned source file differs from its recorded SHA-256."""


class SealedRowError(ReconError):
    """A row or record of a brain sample outside the training partition reached a loader or a fit."""


class PairingError(ReconError):
    """Records of different units, grains or namespaces were about to be joined as paired."""


class RegistryError(ReconError):
    pass


class ForbiddenInputError(ReconError):
    """An input that must never be a predictor for the requested property was supplied."""


class ScaleError(ReconError, ValueError):
    """A value does not state the registered scale and frame of its variable, or cannot be converted."""


# ----------------------------------------------------------------------------------------------
# the unchanged joint model, verified at import
# ----------------------------------------------------------------------------------------------

BRM_JOINT_REL = "analysis/brm_joint.py"
BRM_JOINT_SHA256 = "ae7f5ec432bcbd2b0eebb65595549bc12158b1b7c047f5ace167c7840c0c937e"
BRM_JOINT_PIN_FILE = "data/interim/brm_template_method/frozen.json"


JOINT_INFO = {}


def load_joint_core(path=None, *, expected=BRM_JOINT_SHA256, pin_file=None, module_name="synarch_recon_brm_joint"):
    """Verify analysis/brm_joint.py against its pinned SHA-256 and execute exactly the verified bytes.

    The file is read once; the hash and the executed code come from the same bytes, so the file
    cannot change between the check and the import. Nothing in the module is modified; the check's
    record is kept in JOINT_INFO.
    """
    path = Path(path) if path is not None else ROOT / BRM_JOINT_REL
    source = path.read_bytes()
    actual = sha_bytes(source)
    if actual != expected:
        raise FrozenSourceChanged(f"{path} has SHA-256 {actual}, the pinned joint model is {expected}")
    pin_path = ROOT / BRM_JOINT_PIN_FILE if pin_file is None else Path(pin_file)
    pin_status = "pin file not found; constant checked only"
    if pin_path.exists():
        pinned = json.loads(pin_path.read_text()).get("source_sha256", {}).get(BRM_JOINT_REL)
        if pinned != expected:
            raise FrozenSourceChanged(f"{pin_path} pins {BRM_JOINT_REL} at {pinned}, this layer expects {expected}")
        pin_status = "matches " + str(BRM_JOINT_PIN_FILE)
    module = types.ModuleType(module_name)
    module.__file__ = str(path)
    sys.modules[module_name] = module
    exec(compile(source, str(path), "exec"), module.__dict__)
    JOINT_INFO[module_name] = {"path": str(path), "sha256": actual, "pin": pin_status}
    return module


JOINT = load_joint_core()

# ----------------------------------------------------------------------------------------------
# vocabulary
# ----------------------------------------------------------------------------------------------

KINDS = ("continuous", "categorical", "vector", "structured")
STATUSES = ("measured", "reference", "inferred", "missing")
GRAINS = ("neuron", "nucleus", "cell", "cell_type_population", "brain_sample", "region", "study_record")
POPULATION_GRAINS = ("cell_type_population", "brain_sample", "region")
GROUP_KINDS = ("brain_sample", "donor", "experiment", "study", "none")
TRANSFORMS = ("identity", "log10_1p", "none")
UNCERTAINTY_KINDS = ("none_stated", "sd", "interval", "inferred")
JOIN_KINDS = ("paired_same_unit", "reference_population", "feature_axis_only", "hypothesis_only", "forbidden")
REFUSAL_CODES = ("unsupported_variable", "unidentifiable_cross_modal", "not_validated", "no_better_than_baseline",
                 "out_of_support", "forbidden_input", "invalid_request", "leakage_guard", "invalid_estimate")
VALIDATION_STATES = ("not_validated", "no_better_than_baseline", "exploratory_supported")
BLOCKED_PREFIX = "blocked"


def forward_transform(name, values):
    """Measurement scale -> modelling scale."""
    x = np.asarray(values, dtype=float)
    if name == "identity":
        return x
    if name == "log10_1p":
        if (x < 0).any():
            raise ValueError("log10(1 + x) needs non-negative values")
        return np.log10(1.0 + x)
    raise ValueError(f"no numerical transform named '{name}'")


def inverse_transform(name, values):
    """Modelling scale -> measurement scale. Exact for quantiles; NOT an expected value for means."""
    y = np.asarray(values, dtype=float)
    if name == "identity":
        return y
    if name == "log10_1p":
        return np.power(10.0, y) - 1.0
    raise ValueError(f"no numerical transform named '{name}'")


# ----------------------------------------------------------------------------------------------
# one statement of scale and frame per value, and one conversion
# ----------------------------------------------------------------------------------------------

MODEL_SCALE_TEXT = {"log10_1p": "log10(1 + value)"}     # value tag -> the text that states it
VALUE_TAGS = ("identity",) + tuple(MODEL_SCALE_TEXT)


def scale_statement(spec, transform="identity"):
    """THE one text by which a value of a registered variable states its scale.

    transform 'identity': the value is in the registered measurement scale and units.
    The variable's own modelling transform (for example 'log10_1p'): the value is already on that scale.
    Any other tag, and any other text, is refused wherever a record or an estimate is checked."""
    if transform == "identity":
        return f"{spec.scale}; {spec.units}"
    if transform in MODEL_SCALE_TEXT and transform == spec.transform:
        return MODEL_SCALE_TEXT[transform]
    raise ScaleError(f"{spec.name}: no registered scale goes with the value tag {transform!r}")


def statement_problem(spec, scale, transform, frame):
    """None when the scale text, the value tag and the frame are the registered ones; else what is wrong.

    A missing statement is a mismatch: an empty or unknown scale, a frame that is not stated for a variable
    that has one, units other than the registered ones, a tag that contradicts the text."""
    try:
        expected = scale_statement(spec, transform)
    except ScaleError as error:
        return str(error)
    if scale != expected:
        return (f"{spec.name}: the stated scale {scale!r} is not the registered statement {expected!r} "
                f"for the value tag {transform!r}")
    if frame != spec.frame:
        return (f"{spec.name}: the stated frame is not the registered frame"
                + (" (no frame is stated)" if frame is None else ""))
    return None


def model_values(spec, record):
    """The value of a record on the MODELLING scale of its variable, converted exactly once.

    The record must state the registered scale and frame (statement_problem). A value tagged 'identity' is in
    measurement units and is transformed once; a value tagged with the variable's modelling transform is
    already there and is returned unchanged. Anything else raises ScaleError."""
    problem = statement_problem(spec, record.scale, record.transform, record.frame)
    if problem:
        raise ScaleError(problem)
    wanted = spec.transform if spec.transform in MODEL_SCALE_TEXT else "identity"
    try:
        values = np.atleast_1d(np.asarray(record.value, dtype=float))
    except (TypeError, ValueError):
        raise ScaleError(f"{spec.name}: the value is not numeric") from None
    if not np.isfinite(values).all():
        raise ScaleError(f"{spec.name}: non-finite value")
    if record.transform == "identity" and wanted != "identity":
        if wanted == "log10_1p" and (values < 0).any():
            raise ScaleError(f"{spec.name}: a negative value cannot be put on the scale log10(1 + value)")
        values = forward_transform(wanted, values)
    if wanted == "log10_1p" and (values < 0).any():
        raise ScaleError(f"{spec.name}: a value below zero on the scale log10(1 + value)")
    return values


# ----------------------------------------------------------------------------------------------
# identifiers
# ----------------------------------------------------------------------------------------------


def identifier_problem(value, what="identifier"):
    """None for a usable identifier: a non-empty string without surrounding blanks or control characters."""
    if not isinstance(value, str) or not value:
        return f"{what} is empty or not a string"
    if value != value.strip() or any(ord(ch) < 32 or ord(ch) == 127 for ch in value):
        return f"{what} {value!r} has leading or trailing blanks or control characters"
    return None


def fold_id(value):
    """Comparison key of an identifier: compatibility-normalised, stripped, case-folded."""
    return unicodedata.normalize("NFKC", str(value)).strip().casefold()


def id_in(value, ids, folded=None):
    """'exact' when the identifier is one of `ids`, 'folded' when it only matches after fold_id, else ''.

    Membership of a unit or group in a training state is decided on the folded form: a query that spells a
    training sample with another case or with blanks is treated as that sample, never as a new one."""
    if value is None:
        return ""
    if value in ids:
        return "exact"
    folded = folded if folded is not None else {fold_id(i) for i in ids}
    return "folded" if fold_id(value) in folded else ""


# ----------------------------------------------------------------------------------------------
# registry: namespaces, variables, joins
# ----------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class NamespaceSpec:
    name: str
    grain: str
    unit_key: str
    group_kind: str
    group_evidence: str
    species: str
    readiness: str = "ready"
    note: str = ""
    group_rule: str = ""          # 'unit_prefix': the unit ID itself names its group (text before the first '_')


@dataclass(frozen=True)
class VariableSpec:
    name: str
    namespace: str
    kind: str                     # continuous | categorical | vector | structured
    units: str
    scale: str
    frame: str | None = None      # coordinate frame, if any
    transform: str = "identity"   # transformation used for modelling
    components: tuple = ()        # coordinate names of a vector
    n_components: int | None = None   # when the coordinates are too many to list (genes)
    family: str = ""              # 'axon' for everything derived from the axon
    derived_from: tuple = ()      # variables this one is a function of
    predictor_ban: str | None = None  # reason, when the variable may never be a predictor
    measured_input_needs_declared_task: bool = False
    note: str = ""


@dataclass(frozen=True)
class JoinRule:
    left: tuple
    right: tuple
    kind: str          # one of JOIN_KINDS
    level: str
    key: str | None
    contract_from: str
    contract_to: str
    contract_allowed: object

    def touches(self, a, b):
        return (a in self.left and b in self.right) or (b in self.left and a in self.right)


# Contract dataset key -> (grain, unit key, group kind, group evidence, species, readiness).
_NAMESPACE_SEED = {
    "jiao": ("neuron", "neuron ID <brain sample>_<n>; row key = published full ID", "brain_sample",
             "ID prefix = 'fMOST Mouse ID' of Nature Supplementary Table 1 (166 rows, 166 unique IDs, "
             "reports/d109_jiao_conventions.md section 6). The paper reports 166 mice and the data hold 166 sample "
             "IDs; no statement was found that no mouse gave two samples, so a brain sample is not asserted to be "
             "a distinct animal. Every sample carries exactly one source label.", "Mus musculus (male)", "ready"),
    "epiretro": ("nucleus", "Cell Name", "experiment",
                 "'Experiment' label only; there is no explicit animal-ID column, so an experiment is not an animal",
                 "Mus musculus", "ready"),
    "ABC_RNA": ("cell_type_population", "CCN20230722 cluster ID", "none",
                "reference cluster means of other cells and mice; no donor attached to a mean", "Mus musculus", "ready"),
    "MERFISH_metadata": ("cell", "cell_label within donor", "donor", "single donor C57BL6J-638850 (donor_label column)",
                         "Mus musculus", "ready"),
    "HypoMap": ("cell", "obs/_index in its own namespace", "donor",
                "Sample_ID / donor_id annotation columns; not verified by this layer", "Mus musculus", "ready"),
    "Allen_population_connectivity": ("brain_sample", "Allen experiment ID", "experiment",
                                      "one tracer injection experiment; a population, not a neuron", "Mus musculus",
                                      "ready"),
    "curated_literature": ("study_record", "record_id within study_id", "study", "study_id of the curated record",
                           "as recorded per study", "ready"),
    "Rossi_2021": ("cell", "cell barcode (unreadable locally)", "none",
                   "local HDF5 matrices cannot be opened; no unit can be verified", "Mus musculus",
                   "blocked_invalid_local_HDF5_signature"),
}

# Contract linkage node -> namespaces it speaks about.
_CONTRACT_NODES = {
    "Jiao_soma_tag": ("jiao",), "Jiao_regional_lengths": ("jiao",), "Jiao_soma_tag_regional_lengths": ("jiao",),
    "Jiao_local_trace": ("jiao",), "Jiao_tag": ("jiao",), "Jiao_axon_geometry": ("jiao",), "Jiao_neuron_ID": ("jiao",),
    "Jiao_terminal_row_number": ("jiao_terminal_rows",),
    "Epi_global_methylation_and_clusters": ("epiretro",), "Epi_retro_target": ("epiretro",), "Epi_nucleus": ("epiretro",),
    "ABC_RNA_means": ("ABC_RNA",), "ABC_RNA_mean": ("ABC_RNA",), "ABC_RNA": ("ABC_RNA",),
    "MERFISH_position": ("MERFISH_metadata",), "HypoMap_RNA": ("HypoMap",),
    "HypoMap_or_ABC_RNA": ("HypoMap", "ABC_RNA"),
    "any_local_individual_transcriptome": ("HypoMap", "Rossi_2021", "MERFISH_metadata"),
}

_REFERENCE_LEVELS = ("inferred_taxonomic_reference", "taxonomic_reference")


def _join_kind(edge):
    allowed, level = edge.get("allowed"), edge.get("level")
    if allowed is True and isinstance(level, str) and level.startswith("same_"):
        return "paired_same_unit"
    if allowed is False:
        return "forbidden"
    if isinstance(allowed, str) and level in _REFERENCE_LEVELS:
        return "reference_population"
    if isinstance(allowed, str) and level == "shared_feature_namespace":
        return "feature_axis_only"
    if isinstance(allowed, str) and level == "coarse_marker_inference":
        return "hypothesis_only"
    raise RegistryError(f"the data contract holds a linkage edge this layer cannot classify: {edge}")


JIAO_REGIONS = tuple(("ARH PMv TU TMv MM AHN VMH DMH PVH MPO MPN PH LHA LPO ZI ME MRN PAG SCm SNr VTA DR MY IO NTS PARN "
                      "IRN PB PRNr PCG PA BMA BLA MS BST SI MEA CEA LS ACB LH PVT MD PT AM SPFm HPF Isocortex OLF").split())
# The first 16 published columns are the 16 descendants of HY among the 49 (checked against
# data/raw/atlas/structure_graph_1.json, metadata only); the other 33 lie outside the hypothalamus.
JIAO_HYPOTHALAMIC = JIAO_REGIONS[:16]
SOMA_AXES = ("AP", "DV", "ML")
JIAO_FRAME = ("Allen CCFv3 (2017), 10 um grid, micrometres, as published: +AP posterior, +DV ventral, mirrored to "
              "ML <= 5,700 (axes, unit and mirroring are INFERRED for this table; reports/d109_jiao_conventions.md)")
# Namespaces whose unit ID encodes its group (the Jiao neuron ID is <brain sample>_<n>; the loaders check it).
_GROUP_RULES = {"jiao": "unit_prefix"}
# The ONLY measured inputs of a hidden-axon forecast (reports/axon_dataset_plan.md, section 8, rule 4), and the
# only measured ancestors an inferred input of such a forecast may have. Everything else is refused.
HIDDEN_AXON_INPUTS = ("soma_position", "source_label", "soma_label_fine", "hemisphere", "reconstruction_type")
INFERRED_SOURCE_PREFIX = "reconstruct:"          # source path of every record made from an inferred value


class Registry:
    """Namespaces, variables and join rules. Nothing here holds data."""

    def __init__(self):
        self.namespaces = {}
        self.variables = {}
        self.joins = []
        self.contract = None

    # -- construction --------------------------------------------------------------------------
    def add_namespace(self, spec):
        if spec.grain not in GRAINS or spec.group_kind not in GROUP_KINDS:
            raise RegistryError(f"namespace {spec.name}: unknown grain or group kind")
        if spec.name in self.namespaces:
            raise RegistryError(f"namespace {spec.name} is already registered")
        self.namespaces[spec.name] = spec
        return spec

    def add_variable(self, spec):
        if spec.namespace not in self.namespaces:
            raise RegistryError(f"variable {spec.name}: namespace {spec.namespace} is not registered")
        if spec.kind not in KINDS or spec.transform not in TRANSFORMS:
            raise RegistryError(f"variable {spec.name}: unknown kind or transform")
        if spec.kind == "vector" and not spec.components and not spec.n_components:
            raise RegistryError(f"vector variable {spec.name} needs its coordinates")
        if spec.kind != "vector" and spec.components:
            raise RegistryError(f"variable {spec.name}: only vectors have coordinates")
        if len(set(spec.components)) != len(spec.components):
            raise RegistryError(f"variable {spec.name}: duplicate coordinate names")
        key = (spec.namespace, spec.name)
        if key in self.variables:
            raise RegistryError(f"variable {spec.name} is already registered in {spec.namespace}")
        self.variables[key] = spec
        return spec

    def add_join(self, rule):
        if rule.kind not in JOIN_KINDS:
            raise RegistryError(f"unknown join kind {rule.kind}")
        for name in rule.left + rule.right:
            if name not in self.namespaces:
                raise RegistryError(f"join rule mentions an unregistered namespace {name}")
        self.joins.append(rule)
        return rule

    # -- queries -------------------------------------------------------------------------------
    def namespace(self, name):
        if name not in self.namespaces:
            raise RegistryError(f"namespace {name} is not registered")
        return self.namespaces[name]

    def has_variable(self, namespace, name):
        return (namespace, name) in self.variables

    def variable(self, namespace, name):
        if (namespace, name) not in self.variables:
            raise RegistryError(f"variable {name} is not registered in namespace {namespace}")
        return self.variables[(namespace, name)]

    def namespaces_measuring(self, name):
        return sorted(ns for ns, var in self.variables if var == name)

    def implied_group(self, namespace, unit_id):
        """The group that a unit ID itself names, where the namespace's IDs encode it; None otherwise."""
        if namespace in self.namespaces and self.namespaces[namespace].group_rule == "unit_prefix":
            return str(unit_id).split("_", 1)[0]
        return None

    def joins_between(self, a, b):
        return [rule for rule in self.joins if rule.touches(a, b)]

    def reference_rule(self, from_namespace, to_namespace):
        for rule in self.joins_between(from_namespace, to_namespace):
            if rule.kind == "reference_population":
                return rule
        return None

    def describe(self):
        return {"contract": self.contract,
                "namespaces": {k: jsonable(v) for k, v in sorted(self.namespaces.items())},
                "variables": [jsonable(v) for _, v in sorted(self.variables.items())],
                "joins": [jsonable(j) for j in self.joins]}

    # -- seeding from the data contract --------------------------------------------------------
    @classmethod
    def from_contract(cls, path=None):
        """Namespaces, grains and allowed joins exactly as the local data contract states them."""
        path = Path(path) if path is not None else ROOT / "data/interim/brm_data_contract.json"
        contract = json.loads(path.read_text())
        registry = cls()
        registry.contract = {"path": str(path), "sha256": sha_file(path), "status": contract.get("status")}
        datasets = contract["datasets"]
        if set(datasets) != set(_NAMESPACE_SEED):
            raise RegistryError("the data contract's datasets differ from the namespaces declared in this layer: "
                                f"{sorted(set(datasets) ^ set(_NAMESPACE_SEED))}")
        for name, (grain, key, group_kind, evidence, species, readiness) in _NAMESPACE_SEED.items():
            registry.add_namespace(NamespaceSpec(name, grain, key, group_kind, evidence, species, readiness,
                                                 note=str(datasets[name].get("unit", "")),
                                                 group_rule=_GROUP_RULES.get(name, "")))
        # The terminal-count table has no unit key at all; it gets its own namespace so nothing joins it.
        registry.add_namespace(NamespaceSpec("jiao_terminal_rows", "neuron", "numerical row index only (unverified order)",
                                             "none", "no neuron key in the file", "Mus musculus (male)", "no_unit_key"))
        for edge in contract["linkage_edges"]:
            if edge["from"] not in _CONTRACT_NODES or edge["to"] not in _CONTRACT_NODES:
                raise RegistryError(f"unknown linkage node in the data contract: {edge['from']} -> {edge['to']}")
            registry.add_join(JoinRule(_CONTRACT_NODES[edge["from"]], _CONTRACT_NODES[edge["to"]], _join_kind(edge),
                                       str(edge.get("level")), edge.get("key"), edge["from"], edge["to"],
                                       edge.get("allowed")))
        jiao = datasets["jiao"]
        if tuple(jiao["regional_columns"]) != JIAO_REGIONS:
            raise RegistryError("the contract's 49 regional columns differ from the declared region list")
        _seed_variables(registry, jiao)
        return registry


def _seed_variables(registry, jiao):
    add = registry.add_variable
    add(VariableSpec("soma_position", "jiao", "vector", "um", "linear position", JIAO_FRAME, "identity", SOMA_AXES,
                     note="published soma table; ML is the mirrored coordinate"))
    add(VariableSpec("source_label", "jiao", "categorical", "label", "nominal", None, "none",
                     note="the brain sample's Cre line, or the orexin promoter label 'Wt(orexin)' (portal: 'Orexin'). "
                          "A genetic label, not measured gene expression; one label per brain sample, so a label "
                          "effect is confounded with sample and batch."))
    add(VariableSpec("soma_label_fine", "jiao", "categorical", "label", "nominal (59 CCF acronyms)", None, "none",
                     derived_from=("soma_position",),
                     note="authors' fine soma label (all_neuron_soma_infor.xlsx, column location); how it was "
                          "assigned is unresolved; a function of the soma position, so never evidence for it"))
    add(VariableSpec("regional_axon_length", "jiao", "vector", "um (indicated by the source, not resolved)",
                     "length >= 0 per region; both hemispheres summed (indicated, not stated)",
                     "49 CCF acronyms (25 leaves, 24 parents; MY is the ancestor of IO, NTS, PARN, IRN)", "log10_1p",
                     JIAO_REGIONS, family="axon", measured_input_needs_declared_task=True,
                     note="published per-region axon length; the columns are not a partition and do not sum to the "
                          "total axon length"))
    add(VariableSpec("axon_arbor", "jiao", "structured", "um", "rooted 3D tree",
                     "Allen CCFv3 (2017), mirrored-left frame, complete measured tree", "none", family="axon",
                     note="3D trace; produced only through the axon adapter"))
    add(VariableSpec("hemisphere", "jiao", "categorical", "label", "nominal (Left / Right)", None, "none",
                     note="true soma hemisphere, portal index only; not loaded by this layer"))
    add(VariableSpec("reconstruction_type", "jiao", "categorical", "label", "nominal", None, "none",
                     note="Axon_and_dendrite or Axon_only, portal index only; not loaded by this layer"))
    banned = {"axon-derived_cluster": "projectome_cluster", "projectome_type": "projectome_type",
              "Class": "projectome_class"}
    if set(jiao.get("excluded_predictors", [])) != set(banned):
        raise RegistryError("the contract's excluded predictors differ from the ones declared in this layer")
    for contract_name, name in banned.items():
        add(VariableSpec(name, "jiao", "categorical", "label", "nominal", None, "none", family="axon",
                         predictor_ban="derived from the axons by the authors; never a predictor",
                         note=f"contract excluded predictor '{contract_name}'; never loaded"))
    add(VariableSpec("retrograde_target", "epiretro", "categorical", "label", "nominal", None, "none",
                     note="one sampled retrograde target per nucleus; other targets are unobserved, not negatives"))
    add(VariableSpec("methylation_global", "epiretro", "vector", "fraction", "fraction in [0, 1]", None, "identity",
                     ("mCGFrac", "mCHFrac"), note="global fractions; no per-locus matrix exists locally"))
    add(VariableSpec("assigned_subclass", "epiretro", "categorical", "label", "nominal", None, "none",
                     note="assigned through methylome-RNA co-clustering: an inferred label, not measured RNA"))
    add(VariableSpec("gene_expression", "ABC_RNA", "vector", "log2(CPM+1)", "cluster mean", None, "identity",
                     n_components=32285, note="reference means of other cells and mice"))
    add(VariableSpec("gene_expression", "HypoMap", "vector", "stored scale not established", "per cell", None,
                     "identity", n_components=51216, note="own cell namespace; no position, projection or methylome"))
    add(VariableSpec("gene_expression", "Rossi_2021", "vector", "unreadable", "per cell", None, "identity",
                     n_components=1, note="local matrices cannot be opened"))
    add(VariableSpec("ccf_position", "MERFISH_metadata", "vector", "mm", "linear position", "Allen CCF (mm)", "identity",
                     ("x", "y", "z")))
    add(VariableSpec("regional_tracer_signal", "Allen_population_connectivity", "vector", "as published",
                     "population tracer measurement", None, "identity", n_components=1))
    add(VariableSpec("readout", "curated_literature", "structured", "study specific", "study specific", None, "none"))


# ----------------------------------------------------------------------------------------------
# observation records and the pairing rule
# ----------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Group:
    """The unit of independence a record belongs to, as verified in the source."""
    group_id: str | None
    kind: str
    evidence: str

    def __post_init__(self):
        if self.kind not in GROUP_KINDS:
            raise ValueError(f"unknown group kind {self.kind}")


NO_GROUP = Group(None, "none", "no group identifier verified in the source")


@dataclass(frozen=True)
class SourceRef:
    path: str
    sha256: str | None
    row_key: str | None


@dataclass(frozen=True)
class Uncertainty:
    kind: str = "none_stated"
    sd: object = None
    lower: object = None
    upper: object = None
    level: float | None = None
    producing_module: str | None = None
    state_hash: str | None = None
    training_groups_hash: str | None = None
    coverage: str | None = None           # 'validated' | 'unvalidated', inferred values only
    note: str = ""

    def __post_init__(self):
        if self.kind not in UNCERTAINTY_KINDS:
            raise ValueError(f"unknown uncertainty kind {self.kind}")
        if self.kind == "inferred":
            if not (self.producing_module and self.state_hash and self.training_groups_hash):
                raise ValueError("an inferred value names its producing module, state hash and training-group hash")
            if self.coverage not in ("validated", "unvalidated"):
                raise ValueError("an inferred value says whether its interval coverage was validated")
        elif self.producing_module or self.state_hash or self.training_groups_hash or self.coverage:
            raise ValueError("only inferred values carry a producing module, state hash or coverage flag")


@dataclass(frozen=True)
class ReferenceLink:
    """A declared population-level alignment to another dataset (never a unit-level pairing)."""
    from_namespace: str
    key_name: str
    key_value: str
    level: str
    contract_allowed: str


@dataclass(frozen=True)
class LineageStep:
    """One inference in the ancestry of an inferred value: what was produced, by which state, from what,
    and for which unit. A lineage is what a record SAYS about itself; what the layer itself recorded when it
    issued the value is in its ledger (ProvenanceLedger), and that is what the hidden-axon rule reads."""
    variable: str
    module: str
    state_hash: str
    training_groups_hash: str
    evidence: tuple                 # ((variable, status of that evidence), ...)
    unit: tuple = ()                # (namespace, grain, unit ID) the inference was made for
    group_id: str | None = None     # that unit's group


@dataclass(frozen=True)
class Observation:
    species: str
    namespace: str
    grain: str
    unit_id: str
    group: Group
    variable: str
    status: str
    value: object
    scale: str
    coordinates: tuple = ()
    frame: str | None = None
    uncertainty: Uncertainty = Uncertainty()
    source: SourceRef = SourceRef("", None, None)
    reference: ReferenceLink | None = None
    derived_from: tuple = ()        # for inferred values: the direct evidence ((variable, status), ...)
    note: str = ""
    transform: str | None = None    # REQUIRED for every record that holds a value (there is no default):
                                    # 'identity': the value is in the registered measurement scale and units;
                                    # 'log10_1p': the value is already on the modelling scale log10(1 + value).
                                    # `scale` must be the registered text for that tag (scale_statement).
    lineage: tuple = ()             # for inferred values: this inference first, then those of its inferred inputs

    def __post_init__(self):
        if self.status not in STATUSES:
            raise ValueError(f"unknown status {self.status}")
        if self.grain not in GRAINS:
            raise ValueError(f"unknown grain {self.grain}")
        if not isinstance(self.unit_id, str) or not self.unit_id:
            raise ValueError("a record needs a non-empty unit ID")
        if self.status == "missing" and self.value is not None:
            raise ValueError("a missing record holds no value (a missing value is not a zero)")
        if self.status != "missing" and self.value is None:
            raise ValueError("a record without a value has status 'missing'")
        if self.status == "inferred" and self.uncertainty.kind != "inferred":
            raise ValueError("an inferred record carries inferred-value uncertainty (module, state, coverage)")
        if self.status != "inferred" and self.uncertainty.kind == "inferred":
            raise ValueError("only inferred records carry inferred-value uncertainty")
        if self.status == "reference":
            if self.reference is None:
                raise ValueError("a reference record names its declared alignment")
            if self.grain not in POPULATION_GRAINS:
                raise ValueError("a reference alignment is population-level: it attaches to a group, not to a unit")
        elif self.reference is not None:
            raise ValueError("only reference records carry a reference alignment; it is never 'measured'")
        if self.coordinates and self.status != "missing":
            if not isinstance(self.value, tuple) or len(self.value) != len(self.coordinates):
                raise ValueError("a vector value is a tuple aligned with its coordinates")
        if self.status != "missing":
            # one statement of scale: the tag is required, and the text may not contradict it
            if self.transform not in VALUE_TAGS:
                raise ValueError("a record states the scale of its value: transform='identity' (registered "
                                 "measurement scale and units) or the variable's modelling transform; got "
                                 f"{self.transform!r}")
            on_model_scale = self.transform != "identity"
            if not isinstance(self.scale, str) or (self.scale in MODEL_SCALE_TEXT.values()) != on_model_scale \
                    or (on_model_scale and self.scale != MODEL_SCALE_TEXT[self.transform]):
                raise ValueError(f"the scale text {self.scale!r} contradicts the value tag {self.transform!r}")
        if self.status == "inferred":
            first, u = (self.lineage[0] if self.lineage else None), self.uncertainty
            if not isinstance(first, LineageStep) or first.variable != self.variable \
                    or first.state_hash != u.state_hash or first.module != u.producing_module \
                    or first.training_groups_hash != u.training_groups_hash \
                    or tuple(first.evidence) != tuple(self.derived_from):
                raise ValueError("an inferred record carries its lineage: its own inference first, with the module, "
                                 "state, training groups and evidence that its uncertainty and derived_from state")
            for position, step in enumerate(self.lineage):
                if not isinstance(step, LineageStep):
                    raise ValueError("a lineage is made of LineageStep entries")
                if tuple(step.unit) != self.key() or step.group_id != self.group.group_id:
                    raise ValueError("every step of a lineage belongs to the record's own unit and group: a value "
                                     "inferred for one unit never becomes a record of another")
                later = {s.variable for s in self.lineage[position + 1:] if isinstance(s, LineageStep)}
                for ancestor, ancestor_status in step.evidence:
                    if ancestor_status == "inferred" and ancestor not in later:
                        raise ValueError(f"the lineage is incomplete: the inference that produced {ancestor} is missing")
        else:
            if self.lineage:
                raise ValueError("only inferred records carry a lineage")
            if self.derived_from:
                raise ValueError("only inferred records are derived from other records")
            if isinstance(self.source, SourceRef) and str(self.source.path).startswith(INFERRED_SOURCE_PREFIX):
                raise ValueError(f"a record whose source is a reconstruction is 'inferred', not '{self.status}'")

    def key(self):
        return (self.namespace, self.grain, self.unit_id)


def paired(a, b):
    """THE PAIRING RULE: identical namespace, grain and unit ID, and nothing else."""
    return a.key() == b.key()


def lineage_roots(record):
    """{(variable, status)} of every NON-inferred item in the ancestry of an inferred record: what it was
    ultimately computed from, as its lineage states it."""
    return {(variable, status) for step in record.lineage for variable, status in step.evidence
            if status != "inferred"}


def record_digest(record):
    """SHA-256 of everything a later provenance check relies on in a record: unit, group, variable, value,
    scale statement, producing module and state, direct evidence and the whole lineage."""
    u = record.uncertainty
    return digest({"key": list(record.key()), "group": [record.group.group_id, record.group.kind],
                   "variable": record.variable, "status": record.status, "coordinates": list(record.coordinates),
                   "value": record.value, "scale": record.scale, "transform": record.transform, "frame": record.frame,
                   "module": u.producing_module, "state_hash": u.state_hash,
                   "training_groups_hash": u.training_groups_hash, "coverage": u.coverage,
                   "interval": [u.lower, u.upper, u.level], "derived_from": record.derived_from,
                   "lineage": record.lineage})


class ObservationStore:
    """Records by unit. There is no operation that joins records of different units."""

    def __init__(self, registry):
        self.registry = registry
        self._units = defaultdict(dict)

    def __len__(self):
        return sum(len(v) for v in self._units.values())

    def add(self, record):
        namespace = self.registry.namespace(record.namespace)
        if record.status == "reference":
            rule = self.registry.reference_rule(record.reference.from_namespace, record.namespace)
            if rule is None:
                raise PairingError(f"no population-level reference alignment from {record.reference.from_namespace} "
                                   f"to {record.namespace} is declared in the data contract")
            spec = self.registry.variable(record.reference.from_namespace, record.variable)
        else:
            if record.grain != namespace.grain:
                raise PairingError(f"{record.namespace} records are at grain {namespace.grain}, not {record.grain}")
            spec = self.registry.variable(record.namespace, record.variable)
            if spec.predictor_ban and record.status != "missing":
                raise ForbiddenInputError(f"{record.variable} is never loaded or used: {spec.predictor_ban}")
            if record.group.kind not in (namespace.group_kind, "none"):
                raise PairingError(f"{record.namespace} units are grouped by {namespace.group_kind}, not by "
                                   f"{record.group.kind}")
            implied = self.registry.implied_group(record.namespace, record.unit_id)
            if implied is not None and record.group.group_id is not None and record.group.group_id != implied:
                raise PairingError(f"by its ID, unit {record.unit_id} belongs to group {implied}, not to "
                                   f"{record.group.group_id}")
        problem = identifier_problem(record.unit_id, "unit ID") or (
            None if record.group.group_id is None else identifier_problem(record.group.group_id, "group ID"))
        if problem:
            raise ValueError(f"{record.variable}: {problem}")
        if record.coordinates and spec.components and not set(record.coordinates) <= set(spec.components):
            raise ValueError(f"{record.variable}: coordinates outside the registered ones")
        if record.status != "missing":
            problem = statement_problem(spec, record.scale, record.transform, record.frame)
            if problem:
                raise ScaleError(problem)
        slot = self._units[record.key()]
        if (record.variable, record.status) in slot:
            raise ValueError(f"duplicate record for {record.key()} {record.variable} [{record.status}]")
        for other in slot.values():
            if (other.group.group_id, other.group.kind) != (record.group.group_id, record.group.kind):
                raise PairingError(f"unit {record.key()} is already filed under group {other.group.group_id}; one "
                                   f"unit belongs to one group, not also to {record.group.group_id}")
        slot[(record.variable, record.status)] = record
        return record

    def units(self, namespace, grain):
        return sorted(u for (ns, g, u) in self._units if ns == namespace and g == grain)

    def unit(self, namespace, grain, unit_id, *, statuses=("measured", "reference", "inferred")):
        """The paired view of ONE unit: variable -> record (first available status in `statuses`)."""
        slot = self._units.get((namespace, grain, unit_id), {})
        view = {}
        for status in statuses:
            for (variable, s), record in slot.items():
                if s == status and variable not in view:
                    view[variable] = record
        return view

    def join(self, records):
        """Return variable -> record for records of one unit; refuse anything else."""
        records = list(records)
        if not records:
            return {}
        first = records[0]
        view = {}
        for record in records:
            if not paired(first, record):
                raise PairingError(f"{record.key()} and {first.key()} are different units: records of different "
                                   "cells, grains or datasets are never joined as paired")
            if (record.group.group_id, record.group.kind) != (first.group.group_id, first.group.kind):
                raise PairingError(f"records of {first.key()} name two groups ({first.group.group_id}, "
                                   f"{record.group.group_id}): one unit belongs to one group")
            if record.variable in view:
                raise PairingError(f"two records of {record.variable} for one unit")
            view[record.variable] = record
        return view

    def attach_reference(self, *, from_namespace, to_namespace, key_name, key_value, variable, value, coordinates=(),
                         scale=None, species, source, to_grain="cell_type_population", uncertainty=Uncertainty(),
                         note=""):
        """The ONLY cross-dataset link: a population-level reference value attached to a group.

        The value is on the registered scale of the reference variable; `scale`, when given, must be exactly
        that statement (scale_statement of the variable in `from_namespace`)."""
        rule = self.registry.reference_rule(from_namespace, to_namespace)
        if rule is None:
            kinds = sorted({r.kind for r in self.registry.joins_between(from_namespace, to_namespace)}) or ["none"]
            raise PairingError(f"the data contract declares no population-level reference alignment from "
                               f"{from_namespace} to {to_namespace} (rules found: {', '.join(kinds)})")
        spec = self.registry.variable(from_namespace, variable)
        statement = scale_statement(spec, "identity")
        if scale is not None and scale != statement:
            raise ScaleError(f"{variable}: the stated scale {scale!r} is not the registered statement {statement!r}")
        record = Observation(species, to_namespace, to_grain, f"{key_name}={key_value}",
                             Group(str(key_value), "none", f"group of {to_namespace} units sharing {key_name}"),
                             variable, "reference", value, statement, tuple(coordinates), spec.frame, uncertainty,
                             source, ReferenceLink(from_namespace, key_name, str(key_value), rule.level,
                                                   str(rule.contract_allowed)), note=note, transform="identity")
        return self.add(record)


# ----------------------------------------------------------------------------------------------
# worlds, the sealed training partition and the Jiao loaders
# ----------------------------------------------------------------------------------------------

REAL_PINS = {
    "manifest": "649f829a592cd3d13bcda5cffdbd049e759568d38bacb144f6a61220610a1ca4",
    "soma_csv": "cc01d048bf5a7c9c9f46869fb750f2d67a98ca934903d2e2e71c1a67b1d97cb9",
    "soma_xlsx": "7b7ae75b94688619dbb367d5f9e71d8b157615624932375a00e80c796f21f162",
    "regional_csv": "2072bc20550607c195be0a42cee9fcbe0e0ab5aaebfaee0a1016b2a43d35e63c",
}
SOMA_CSV_COLUMNS = ("file", "neurons", "geno", "location", "AP", "DV", "ML")     # never 'type' or 'name'
SOMA_XLSX_COLUMNS = ("neurons", "location")                                       # never 'cluster' or 'Class'
MIDLINE_ML = 5700.0


@dataclass(frozen=True)
class World:
    """Where the tables of one cohort live. The real cohort, or a synthetic one of the same format."""
    name: str
    synthetic: bool
    manifest: Path
    soma_csv: Path
    soma_xlsx: Path
    regional_csv: Path
    out: Path
    pins: tuple = ()
    expected: tuple = ()          # (("train_samples", n), ("train_neurons", n))
    regions: tuple = JIAO_REGIONS
    t4_observed: tuple = JIAO_HYPOTHALAMIC
    label_alias: tuple = (("Wt(orexin)", "Orexin"),)
    species: str = "Mus musculus (male)"

    def pin(self, name):
        return dict(self.pins).get(name)


def real_world(root=ROOT):
    root = Path(root)
    j = root / "data/raw/jiao2025"
    return World("jiao2025_training_partition", False, root / "reports/axon_dataset_manifest.tsv",
                 j / "All_soma_infor_HY_neurons.csv", j / "all_neuron_soma_infor.xlsx",
                 j / "all_neurons_projection_length_fig1.csv", root / "data/interim/synarch_recon",
                 tuple(sorted(REAL_PINS.items())), (("train_neurons", 4448), ("train_samples", 101)))


def synthetic_world(directory):
    directory = Path(directory)
    meta = read_json(directory / "world.json")
    if meta.get("synthetic") is not True:
        raise ReconError("world.json does not describe a synthetic world")
    return World(meta["name"], True, directory / "manifest.tsv", directory / "soma.csv", directory / "soma.xlsx",
                 directory / "regional.csv", directory / "out", (), tuple(sorted(meta["expected"].items())),
                 tuple(meta["regions"]), tuple(meta["t4_observed"]), tuple(map(tuple, meta.get("label_alias", []))),
                 meta.get("species", "synthetic"))


class TrainingPartition:
    """The training brain samples of the locked split, and nothing about the others.

    Manifest rows of validation and test samples are dropped the moment they are read; not even
    their IDs are kept. Everything downstream is a whitelist: a sample that is not listed here is
    refused, whatever it is.
    """

    def __init__(self, world):
        if not Path(world.manifest).exists():
            raise ReconError(f"{world.manifest} is missing; regenerate it with "
                             "python3 -I analysis/axon_dataset_audit.py")
        self.manifest_sha256 = sha_file(world.manifest)
        if world.pin("manifest") and self.manifest_sha256 != world.pin("manifest"):
            raise FrozenSourceChanged("the cohort manifest differs from the SHA-256 pinned by the locked plan")
        neuron_sample, neuron_label, neuron_soma_label, sample_label = {}, {}, {}, {}
        with open(world.manifest, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["part"] != "train":
                    continue                      # sealed: nothing of this row is kept
                neuron, sample = row["neuron"], row["animal"]
                if neuron in neuron_sample or neuron.split("_", 1)[0] != sample:
                    raise ReconError(f"manifest row {neuron} is duplicated or disagrees with its sample ID")
                if sample_label.setdefault(sample, row["source_label"]) != row["source_label"]:
                    raise ReconError(f"brain sample {sample} carries more than one source label")
                neuron_sample[neuron], neuron_label[neuron] = sample, row["source_label"]
                neuron_soma_label[neuron] = row["soma_label"]
        self.neuron_sample, self.neuron_label, self.neuron_soma_label = neuron_sample, neuron_label, neuron_soma_label
        self.sample_label = dict(sorted(sample_label.items()))
        self.train_samples = frozenset(sample_label)
        self.train_neurons = frozenset(neuron_sample)
        expected = dict(world.expected)
        if expected and (len(self.train_samples) != expected["train_samples"]
                         or len(self.train_neurons) != expected["train_neurons"]):
            raise ReconError(f"training partition has {len(self.train_samples)} samples / {len(self.train_neurons)} "
                             f"neurons; expected {expected['train_samples']} / {expected['train_neurons']}")
        self.samples_sha256 = ids_sha(sorted(self.train_samples))

    def is_training(self, sample):
        return sample in self.train_samples

    def require_training(self, samples, context):
        """Refuse any brain sample that is not a training sample (the canary guard)."""
        outside = sorted({str(s) for s in samples} - self.train_samples)
        if outside:
            raise SealedRowError(f"{context}: brain sample(s) {', '.join(outside[:5])} are not in the training "
                                 "partition; validation and test rows are sealed")


def _check_file_pin(world, name, path):
    actual = sha_file(path)
    if world.pin(name) and actual != world.pin(name):
        raise FrozenSourceChanged(f"{path} differs from its pinned SHA-256")
    return actual


def read_soma_training(world, partition):
    """Allow-listed columns of the soma table, training samples only. 'type' and 'name' are never touched."""
    sha = _check_file_pin(world, "soma_csv", world.soma_csv)
    rows, dropped = {}, 0
    with open(world.soma_csv, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        index = {name: header.index(name) for name in SOMA_CSV_COLUMNS}
        for row in reader:
            neuron = row[index["neurons"]]
            if not partition.is_training(neuron.split("_", 1)[0]):
                dropped += 1                      # sealed: no other cell of this row is read
                continue
            if neuron in rows:
                raise ReconError(f"duplicate soma row {neuron}")
            rows[neuron] = {"full_id": row[index["file"]], "geno": row[index["geno"]],
                            "location": row[index["location"]],
                            "soma": (float(row[index["AP"]]), float(row[index["DV"]]), float(row[index["ML"]]))}
    return rows, {"sha256": sha, "training_rows": len(rows), "sealed_rows_dropped_unparsed": dropped}


def read_fine_labels_training(world, partition):
    """Fine soma label (xlsx column 'location'), training samples only; 'cluster' and 'Class' are never read."""
    import openpyxl

    sha = _check_file_pin(world, "soma_xlsx", world.soma_xlsx)
    workbook = openpyxl.load_workbook(world.soma_xlsx, read_only=True, data_only=True)
    try:
        sheet = workbook.worksheets[0]
        header = [c for c in next(sheet.iter_rows(min_row=1, max_row=1, values_only=True))]
        cols = {name: header.index(name) + 1 for name in SOMA_XLSX_COLUMNS}
        low, high = min(cols.values()), max(cols.values())
        between = [h for h in header[low - 1:high] if h not in SOMA_XLSX_COLUMNS]
        if between:
            raise ReconError(f"columns {between} lie between the allowed xlsx columns; refusing to read across them")
        labels, dropped = {}, 0
        for cells in sheet.iter_rows(min_row=2, min_col=low, max_col=high, values_only=True):
            neuron = cells[cols["neurons"] - low]
            if neuron is None:
                continue
            neuron = str(neuron)
            if not partition.is_training(neuron.split("_", 1)[0]):
                dropped += 1
                continue
            if neuron in labels:
                raise ReconError(f"duplicate xlsx row {neuron}")
            labels[neuron] = str(cells[cols["location"] - low])
    finally:
        workbook.close()
    return labels, {"sha256": sha, "training_rows": len(labels), "sealed_rows_dropped_unparsed": dropped}


def read_regional_training(world, partition):
    """Published regional lengths of TRAINING samples only.

    The file is read line by line; a line whose brain sample is not a training sample is dropped
    before any of its values is parsed, so no validation or test value ever exists in memory.
    """
    sha = _check_file_pin(world, "regional_csv", world.regional_csv)
    lengths, dropped = {}, 0
    with open(world.regional_csv, newline="") as handle:
        header = handle.readline().rstrip("\r\n").split(",")
        if tuple(header[1:]) != tuple(world.regions):
            raise ReconError("the regional table does not hold the declared region columns in the declared order")
        for line in handle:
            cut = line.find(",")
            full_id = line[:cut]
            if not partition.is_training(full_id.split("_", 1)[0]):
                dropped += 1                      # sealed: never parsed, never stored
                continue
            values = np.array(line[cut + 1:].rstrip("\r\n").split(","), dtype=float)
            if values.shape != (len(world.regions),) or not np.isfinite(values).all() or (values < 0).any():
                raise ReconError(f"invalid regional row {full_id}")
            if full_id in lengths:
                raise ReconError(f"duplicate regional row {full_id}")
            lengths[full_id] = values
    return lengths, {"sha256": sha, "training_rows": len(lengths), "sealed_rows_dropped_unparsed": dropped}


@dataclass(frozen=True)
class TrainingTable:
    """Dense arrays of the training partition, ordered by neuron ID. Built only by assemble_training_table."""
    world: str
    neuron: np.ndarray
    full_id: np.ndarray
    sample: np.ndarray
    label: np.ndarray
    fine_label: np.ndarray
    csv_location: np.ndarray
    soma: np.ndarray
    length: np.ndarray | None
    regions: tuple
    samples_sha256: str
    audit: dict = field(default_factory=dict)

    def __len__(self):
        return len(self.neuron)

    def profile(self):
        """Regional profile on the declared modelling scale, log10(1 + published length)."""
        if self.length is None:
            raise ReconError("this table was loaded without outcomes")
        return forward_transform("log10_1p", self.length)

    def subset(self, index):
        index = np.asarray(index)
        return TrainingTable(self.world, self.neuron[index], self.full_id[index], self.sample[index],
                             self.label[index], self.fine_label[index], self.csv_location[index], self.soma[index],
                             None if self.length is None else self.length[index], self.regions, self.samples_sha256,
                             self.audit)

    def sha256(self):
        arrays = [self.neuron, self.sample, self.label, self.soma]
        if self.length is not None:
            arrays.append(self.length)
        return array_sha(*arrays)


def assemble_training_table(world, partition, soma_rows, fine_labels, lengths, audit=None):
    """Join the three tables by neuron ID (never by position) and refuse any non-training sample."""
    partition.require_training({n.split("_", 1)[0] for n in soma_rows}, "soma table")
    partition.require_training({n.split("_", 1)[0] for n in fine_labels}, "fine-label table")
    if lengths is not None:
        partition.require_training({i.split("_", 1)[0] for i in lengths}, "regional table")
    if set(soma_rows) != partition.train_neurons or set(fine_labels) != partition.train_neurons:
        raise ReconError("the soma tables do not hold every training neuron exactly once")
    alias = dict(world.label_alias)
    order = sorted(partition.train_neurons)
    full_ids, by_full = [], {}
    for neuron in order:
        row = soma_rows[neuron]
        sample = partition.neuron_sample[neuron]
        if not row["full_id"].startswith(neuron + "_"):
            raise ReconError(f"published full ID and neuron ID disagree for {neuron}")
        if alias.get(row["geno"], row["geno"]) != partition.neuron_label[neuron]:
            raise ReconError(f"table label and manifest label disagree for {neuron}")
        if fine_labels[neuron] != partition.neuron_soma_label[neuron]:
            raise ReconError(f"fine soma label and manifest soma label disagree for {neuron}")
        if not all(math.isfinite(v) for v in row["soma"]) or row["soma"][2] > MIDLINE_ML + 1e-6:
            raise ReconError(f"soma of {neuron} is not finite or not in the mirrored frame (ML <= 5,700)")
        if row["full_id"] in by_full:
            raise ReconError(f"published full ID {row['full_id']} is not unique")
        if row["full_id"] != f"{neuron}_{row['geno']}_{fine_labels[neuron]}":
            raise ReconError(f"published full ID of {neuron} is not <neuron>_<source label>_<fine soma label>")
        if sample != neuron.split("_", 1)[0]:
            raise ReconError(f"brain sample of {neuron} is not its ID prefix")
        by_full[row["full_id"]] = neuron
        full_ids.append(row["full_id"])
    length = None
    if lengths is not None:
        if set(lengths) != set(by_full):
            raise ReconError("the regional table does not hold every training neuron exactly once")
        length = np.vstack([lengths[i] for i in full_ids])
    sample = np.array([partition.neuron_sample[n] for n in order])
    labels = np.array([soma_rows[n]["geno"] for n in order])
    per_sample = defaultdict(set)
    for s, l in zip(sample, labels):
        per_sample[s].add(l)
    if any(len(v) != 1 for v in per_sample.values()):
        raise ReconError("a brain sample carries more than one source label in the soma table")
    table = TrainingTable(world.name, np.array(order), np.array(full_ids), sample, labels,
                          np.array([fine_labels[n] for n in order]),
                          np.array([soma_rows[n]["location"] for n in order]),
                          np.array([soma_rows[n]["soma"] for n in order], dtype=float), length, tuple(world.regions),
                          partition.samples_sha256, dict(audit or {}))
    partition.require_training(table.sample, "assembled training table")
    return table


def load_training_table(world, *, with_outcomes=True, partition=None):
    """The training partition as dense arrays. Sealed rows are dropped unparsed by every reader."""
    partition = partition or TrainingPartition(world)
    soma_rows, soma_audit = read_soma_training(world, partition)
    fine, fine_audit = read_fine_labels_training(world, partition)
    lengths, regional_audit = (None, None)
    if with_outcomes:
        lengths, regional_audit = read_regional_training(world, partition)
    audit = {"manifest_sha256": partition.manifest_sha256, "training_samples": len(partition.train_samples),
             "training_neurons": len(partition.train_neurons), "training_samples_sha256": partition.samples_sha256,
             "soma_csv": soma_audit, "soma_xlsx": fine_audit, "regional_csv": regional_audit}
    return assemble_training_table(world, partition, soma_rows, fine, lengths, audit)


def jiao_group(sample, registry):
    return Group(str(sample), "brain_sample", registry.namespace("jiao").group_evidence)


def jiao_records(world, table, registry, store=None):
    """The training neurons as provenance-aware records (four measured records per neuron)."""
    store = store if store is not None else ObservationStore(registry)
    namespace = registry.namespace("jiao")
    spec = {name: registry.variable("jiao", name) for name in
            ("soma_position", "source_label", "soma_label_fine", "regional_axon_length")}
    files = {"csv": str(world.soma_csv), "xlsx": str(world.soma_xlsx), "regional": str(world.regional_csv)}
    shas = {"csv": table.audit.get("soma_csv", {}).get("sha256"), "xlsx": table.audit.get("soma_xlsx", {}).get("sha256"),
            "regional": (table.audit.get("regional_csv") or {}).get("sha256")}
    for i, neuron in enumerate(table.neuron):
        neuron, full_id = str(neuron), str(table.full_id[i])
        group = jiao_group(table.sample[i], registry)

        def record(variable, value, coordinates, which, row_key, note=""):
            v = spec[variable]
            return Observation(namespace.species, "jiao", "neuron", neuron, group, variable, "measured", value,
                               scale_statement(v, "identity"), coordinates, v.frame, Uncertainty("none_stated"),
                               SourceRef(files[which], shas[which], row_key), note=note, transform="identity")

        store.add(record("soma_position", tuple(float(v) for v in table.soma[i]), SOMA_AXES, "csv", full_id))
        store.add(record("source_label", str(table.label[i]), (), "csv", full_id,
                         "sample-level genetic label carried by every neuron of the brain sample"))
        store.add(record("soma_label_fine", str(table.fine_label[i]), (), "xlsx", neuron))
        if table.length is not None:
            store.add(record("regional_axon_length", tuple(float(v) for v in table.length[i]), tuple(table.regions),
                             "regional", full_id))
    return store


# ----------------------------------------------------------------------------------------------
# whole-brain-sample splits by a salted hash
# ----------------------------------------------------------------------------------------------


def hash_order(salt, sample):
    return hashlib.sha256(f"{salt}:{sample}".encode()).hexdigest()


def assign_folds(samples, label_of, n_folds, salt):
    """Whole samples -> folds. Inside each source label the samples are ordered by SHA-256(salt:sample);
    the labels are taken in sorted order and the samples are dealt in turn, so fold sizes differ by at most
    one sample and every label is spread as evenly as its sample count allows. Nothing but the salt, the
    sample IDs and the labels enters."""
    samples = sorted(set(samples))
    if len(samples) < n_folds:
        raise ReconError("fewer brain samples than folds")
    ordered = []
    for label in sorted({label_of[s] for s in samples}):
        ordered.extend(sorted((s for s in samples if label_of[s] == label), key=lambda s: hash_order(salt, s)))
    return {sample: position % n_folds for position, sample in enumerate(ordered)}


def nested_folds(samples, label_of, *, outer, inner, salt):
    """Outer folds over all samples; inner folds over the outer-training samples only."""
    samples = sorted(set(samples))
    outer_of = assign_folds(samples, label_of, outer, salt)
    result = {"salt": salt, "outer": {}, "inner": {}}
    for fold in range(outer):
        held = sorted(s for s in samples if outer_of[s] == fold)
        training = sorted(s for s in samples if outer_of[s] != fold)
        inner_of = assign_folds(training, label_of, inner, f"{salt}:inner:{fold}")
        result["outer"][str(fold)] = held
        result["inner"][str(fold)] = {str(i): sorted(s for s in training if inner_of[s] == i) for i in range(inner)}
    check_nested_folds(result, samples)
    return result


def check_nested_folds(splits, samples):
    samples = sorted(set(samples))
    seen = [s for held in splits["outer"].values() for s in held]
    if sorted(seen) != samples:
        raise ReconError("outer folds are not a partition of the brain samples (a sample is missing or in two folds)")
    for fold, held in splits["outer"].items():
        training = set(samples) - set(held)
        inner = [s for part in splits["inner"][fold].values() for s in part]
        if sorted(inner) != sorted(training):
            raise ReconError(f"inner folds of outer fold {fold} are not a partition of its training samples")
        if set(inner) & set(held):
            raise ReconError(f"an inner fold of outer fold {fold} contains a held-out sample")
    return True


# ----------------------------------------------------------------------------------------------
# metrics: per neuron, then per brain sample, then over brain samples
# ----------------------------------------------------------------------------------------------


def profile_scores(prediction, truth):
    """Per-neuron MSE, MAE and Pearson correlation across the target coordinates.

    Correlation is undefined (NaN) for a constant truth profile; a constant prediction scores 0.
    """
    prediction, truth = np.asarray(prediction, float), np.asarray(truth, float)
    error = prediction - truth
    pc = prediction - prediction.mean(axis=1, keepdims=True)
    tc = truth - truth.mean(axis=1, keepdims=True)
    ps, ts = np.sqrt((pc**2).sum(axis=1)), np.sqrt((tc**2).sum(axis=1))
    corr = np.zeros(len(truth))
    ok = (ps > 0) & (ts > 0)
    corr[ok] = (pc[ok] * tc[ok]).sum(axis=1) / (ps[ok] * ts[ok])
    corr[ts == 0] = np.nan
    return {"mse": (error**2).mean(axis=1), "mae": np.abs(error).mean(axis=1), "corr": corr}


def position_scores(prediction, truth, axes=SOMA_AXES):
    prediction, truth = np.asarray(prediction, float), np.asarray(truth, float)
    error = prediction - truth
    out = {"euclid": np.sqrt((error**2).sum(axis=1))}
    for j, axis in enumerate(axes):
        out["abs_" + axis] = np.abs(error[:, j])
    return out


def label_scores(probability, classes, truth, floor=1e-12):
    """Per-neuron correctness and log-loss. NaN where the true label is not among `classes`."""
    probability = np.asarray(probability, float)
    index = {c: i for i, c in enumerate(classes)}
    position = np.array([index.get(t, -1) for t in truth])
    known = position >= 0
    correct = np.full(len(truth), np.nan)
    logloss = np.full(len(truth), np.nan)
    if known.any():
        best = probability[known].argmax(axis=1)
        correct[known] = (best == position[known]).astype(float)
        p = probability[known, position[known]]
        logloss[known] = -np.log(np.maximum(p, floor))
    return {"accuracy": correct, "logloss": logloss}


def by_sample(values, sample, *, how="mean"):
    """{brain sample: mean (or median) of its neurons' values}, NaN entries ignored; all-NaN samples dropped."""
    values, sample = np.asarray(values, float), np.asarray(sample)
    out = {}
    for s in sorted(set(sample.tolist())):
        v = values[sample == s]
        v = v[np.isfinite(v)]
        if len(v):
            out[s] = float(np.median(v) if how == "median" else v.mean())
    return out


def over_samples(per_sample):
    """Equal weight per brain sample."""
    return float(np.mean(list(per_sample.values()))) if per_sample else float("nan")


def paired_bootstrap(differences, seed, *, samples=10000, level=0.95):
    """Percentile interval of the mean paired difference, resampling brain samples. Descriptive only."""
    d = np.asarray(list(differences), float)
    if len(d) == 0:
        return {"n": 0, "mean": float("nan"), "lower": float("nan"), "upper": float("nan"), "level": level}
    rng = np.random.default_rng(seed)
    draws = d[rng.integers(0, len(d), size=(samples, len(d)))].mean(axis=1)
    alpha = (1 - level) / 2
    return {"n": int(len(d)), "mean": float(d.mean()), "lower": float(np.quantile(draws, alpha)),
            "upper": float(np.quantile(draws, 1 - alpha)), "level": level, "resamples": samples}


# ----------------------------------------------------------------------------------------------
# reconstruction-request API
# ----------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Notice:
    code: str
    detail: str
    variable: str | None = None
    handling: str | None = None


@dataclass(frozen=True)
class Refusal:
    code: str
    reason: str
    detail: tuple = ()

    def __post_init__(self):
        if self.code not in REFUSAL_CODES:
            raise ValueError(f"unknown refusal code {self.code}")


@dataclass(frozen=True)
class Estimate:
    variable: str
    kind: str
    coordinates: tuple
    scale: str                       # scale of point, median and interval: the registered statement for `transform`
    point: object
    median: object = None
    interval: dict | None = None     # {'level', 'lower', 'upper'} on `scale`
    probabilities: dict | None = None
    prediction_set: tuple | None = None
    original_units: dict | None = None   # back-transformed quantiles (exact for quantiles, not for means)
    note: str = ""
    transform: str | None = None         # REQUIRED (no default): 'identity', or the variable's modelling transform
                                         # when point / median / interval are on the modelling scale
    frame: str | None = None             # the registered frame of the variable (None only when it has none)


@dataclass(frozen=True)
class ValidationEvidence:
    """What a benchmark established for one declared conditional of a module."""
    status: str = "not_validated"
    task: str | None = None
    benchmark: dict | None = None        # pointer: metrics file, its SHA-256, the pre-run declaration
    primary_metric: str | None = None
    coverage_validated: bool = False
    coverage: dict | None = None
    best_comparator: str | None = None
    best_comparator_params: dict | None = None
    confirmatory: bool = False
    note: str = ""

    def __post_init__(self):
        if self.status not in VALIDATION_STATES and not self.status.startswith(BLOCKED_PREFIX):
            raise ValueError(f"unknown validation status {self.status}")
        if self.confirmatory:
            raise ValueError("no confirmatory evidence exists in this layer; exploratory benchmarks only")
        if self.status == "exploratory_supported" and not self.benchmark:
            raise ValueError("a supported status points to its benchmark evidence")


@dataclass(frozen=True)
class Conditional:
    """One declared conditional of a module: target given evidence."""
    name: str
    target: str
    target_coordinates: tuple | None            # None: the whole variable
    requires: tuple                             # ((variable, coordinates or None), ...)
    task: str | None = None
    declared_measured_axon_input: bool = False  # True only for tasks declared as taking a measured profile
    description: str = ""


@dataclass(frozen=True)
class Target:
    variable: str
    coordinates: tuple | None = None

    @property
    def label(self):
        return self.variable if self.coordinates is None else f"{self.variable}[{len(self.coordinates)}]"


@dataclass(frozen=True)
class Policy:
    """The declared policy for answers that the evidence does not support. Run-time settings: they are part
    of every result (ReconstructionResult.policy), not of any saved state."""
    on_no_better_than_baseline: str = "baseline"    # 'baseline' (flagged baseline-level answer) | 'refuse'
    on_not_validated: str = "refuse"                # 'refuse' | 'baseline' | 'engine_flagged'
    on_unseen_category: str = "marginalise"         # 'marginalise' (drop it and say so) | 'refuse'
    on_out_of_range: str = "refuse"                 # 'refuse' | 'marginalise'
    interval_level: float = 0.90

    def __post_init__(self):
        if self.on_no_better_than_baseline not in ("baseline", "refuse"):
            raise ValueError("on_no_better_than_baseline is 'baseline' or 'refuse'")
        if self.on_not_validated not in ("refuse", "baseline", "engine_flagged"):
            raise ValueError("on_not_validated is 'refuse', 'baseline' or 'engine_flagged'")
        if self.on_unseen_category not in ("marginalise", "refuse") or self.on_out_of_range not in ("refuse", "marginalise"):
            raise ValueError("support policies are 'marginalise' or 'refuse'")
        if not 0 < self.interval_level < 1:
            raise ValueError("interval_level lies strictly between 0 and 1")


@dataclass(frozen=True)
class ReconstructionRequest:
    namespace: str
    grain: str
    unit_id: str
    group: Group
    observed: tuple
    requested: tuple
    species: str = ""
    interval_level: float | None = None

    @classmethod
    def from_values(cls, registry, *, namespace, unit_id, observed, requested, group_id=None, interval_level=None):
        """Wrap plain values as measured records OF THIS UNIT, in the registered measurement scale, units and
        frame of each variable. Vector values are a sequence over all coordinates or a mapping
        coordinate -> value (a partially measured vector)."""
        spec_ns = registry.namespace(namespace)
        group = Group(None if group_id is None else str(group_id), spec_ns.group_kind if group_id is not None else "none",
                      spec_ns.group_evidence if group_id is not None else "no group supplied with the request")
        records = []
        for variable, value in observed.items():
            if isinstance(value, Observation):
                records.append(value)
                continue
            if not registry.has_variable(namespace, variable):
                raise RegistryError(f"{variable} is not a variable of namespace {namespace}; pass a full record "
                                    "if it comes from another dataset")
            spec = registry.variable(namespace, variable)
            coordinates = ()
            if spec.kind == "vector":
                if isinstance(value, Mapping):
                    coordinates = tuple(c for c in spec.components if c in value)
                    if len(coordinates) != len(value):
                        raise ValueError(f"{variable}: unknown coordinate names")
                    value = tuple(float(value[c]) for c in coordinates)
                else:
                    coordinates, value = tuple(spec.components), tuple(float(v) for v in value)
            elif spec.kind == "continuous":
                value = float(value)
            records.append(Observation(spec_ns.species, namespace, spec_ns.grain, str(unit_id), group, variable,
                                       "measured", value, scale_statement(spec, "identity"), coordinates, spec.frame,
                                       source=SourceRef("request", None, str(unit_id)), transform="identity"))
        targets = tuple(t if isinstance(t, Target) else
                        (Target(t[0], tuple(t[1])) if isinstance(t, (tuple, list)) else Target(str(t)))
                        for t in requested)
        return cls(namespace, spec_ns.grain, str(unit_id), group, tuple(records), targets, spec_ns.species,
                   interval_level)


@dataclass
class Plan:
    """How a module would answer one target: which conditional, which evidence, which caveats."""
    module: object
    target: Target
    coordinates: tuple
    conditional: Conditional | None
    evidence_used: tuple = ()        # ((variable, status), ...): exactly what conditions the answer
    evidence_ignored: tuple = ()     # ((variable or variable[coordinates], reason), ...): supplied, not used
    notices: list = field(default_factory=list)
    unit_in_training: bool = False
    group_in_training: bool = False
    payload: object = None


@dataclass(frozen=True)
class BaselineAnswer:
    """A baseline-level answer of a module, with the provenance of the comparator bank that made it."""
    estimate: Estimate
    module: str                      # 'baseline:<comparator>'
    state_hash: str                  # binds the bank's rows and every setting of the comparator that answered
    coverage_validated: bool = False
    training_groups_hash: str | None = None
    training_groups: tuple = ()


@dataclass(frozen=True)
class PropertyResult:
    variable: str
    coordinates: tuple
    status: str                          # 'inferred' | 'refused' | 'reference'
    estimate: Estimate | None = None
    refusal: Refusal | None = None
    module: str | None = None
    state_hash: str | None = None
    training_groups_hash: str | None = None
    coverage_validated: bool = False
    validation: str | None = None
    conditional: str | None = None
    benchmark: dict | None = None
    evidence_used: tuple = ()
    evidence_ignored: tuple = ()
    notices: tuple = ()
    baseline_level: bool = False
    reference: ReferenceLink | None = None
    lineage: tuple = ()
    unit: tuple | None = None            # (namespace, grain, unit ID) this result was computed for
    group: Group | None = None           # that unit's group

    def codes(self):
        """Every structured code attached to this result (the refusal's and the notices')."""
        codes = {n.code for n in self.notices}
        if self.refusal is not None:
            codes.add(self.refusal.code)
        return codes

    def to_observation(self, request):
        """The inferred value as a record of THE UNIT IT WAS INFERRED FOR: status 'inferred', never 'measured',
        with the scale statement and frame of the estimate. A request for another unit or group is refused."""
        if self.status != "inferred" or self.estimate is None:
            raise ReconError("only an inferred value becomes an inferred record")
        key = (request.namespace, request.grain, request.unit_id)
        if self.unit is None or tuple(self.unit) != key:
            raise PairingError(f"this value was inferred for {self.unit}; it cannot become a record of {key}")
        group = self.group if self.group is not None else request.group
        if request.group.group_id is not None and request.group.group_id != group.group_id:
            raise PairingError(f"this value was inferred for a unit of group {group.group_id}; the request names "
                               f"group {request.group.group_id}")
        e = self.estimate
        point = tuple(e.point) if e.coordinates else e.point
        interval = e.interval or {}
        uncertainty = Uncertainty("inferred", None, interval.get("lower"), interval.get("upper"), interval.get("level"),
                                  self.module, self.state_hash, self.training_groups_hash,
                                  "validated" if self.coverage_validated else "unvalidated",
                                  "; ".join(sorted(self.codes())))
        return Observation(request.species or "unspecified", request.namespace, request.grain, request.unit_id,
                           group, self.variable, "inferred", point, e.scale, tuple(e.coordinates), e.frame,
                           uncertainty, SourceRef(INFERRED_SOURCE_PREFIX + str(self.module), self.state_hash,
                                                  request.unit_id),
                           derived_from=tuple(self.evidence_used), note=self.validation or "",
                           transform=e.transform, lineage=tuple(self.lineage))


@dataclass(frozen=True)
class ReconstructionResult:
    request: ReconstructionRequest
    properties: dict
    policy: Policy

    def __getitem__(self, variable):
        return self.properties[variable]

    def to_json(self):
        return {"unit": {"namespace": self.request.namespace, "grain": self.request.grain,
                         "unit_id": self.request.unit_id, "group": jsonable(self.request.group)},
                "policy": jsonable(self.policy),
                "properties": {k: jsonable(v) for k, v in self.properties.items()}}


class Module:
    """Base class of everything that can answer a reconstruction request.

    A module REGISTERS what it consumes and produces, its grain and namespace, and for every declared
    conditional its validation status with a pointer to the benchmark evidence.
    """
    name = ""
    namespace = ""
    grain = ""
    kind = "estimator"

    def produces(self):
        return ()

    def consumes(self):
        return ()

    def conditionals(self):
        return ()

    def evidence_for(self, conditional_name):
        return ValidationEvidence()

    def bind(self, reconstructor):
        """Called once when the module is registered in a Reconstructor."""

    def state_groups(self, state_hash):
        """The training groups of a state (or comparator bank) this module answers with, by its hash;
        None when the module holds no state with that hash."""
        return None

    def plan(self, request, target, evidence, policy):
        """Return a Plan, or a Refusal when this module cannot answer."""
        raise NotImplementedError

    def estimate(self, plan, level):
        raise NotImplementedError

    def baseline(self, plan, level, evidence):
        """A baseline-level answer for the plan: a BaselineAnswer, or None when there is none."""
        return None

    def provenance(self, plan):
        return {"state_hash": None, "training_groups_hash": None, "training_groups": None}

    def spec(self):
        return {"name": self.name, "kind": self.kind, "namespace": self.namespace, "grain": self.grain,
                "consumes": list(self.consumes()), "produces": list(self.produces()),
                "conditionals": [{"name": c.name, "target": c.target,
                                  "target_coordinates": None if c.target_coordinates is None else len(c.target_coordinates),
                                  "requires": [v for v, _ in c.requires], "task": c.task,
                                  "declared_measured_axon_input": c.declared_measured_axon_input,
                                  "validation": jsonable(self.evidence_for(c.name))} for c in self.conditionals()]}


# ----------------------------------------------------------------------------------------------
# provenance: what this layer itself issued, and the one rule for hidden-axon inputs
# ----------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class IssuedFact:
    """What the layer recorded when it issued one inferred value (not what the record says about itself)."""
    digest: str                      # record_digest of the issued record
    unit: tuple                      # (namespace, grain, unit ID)
    group_id: str | None
    variable: str
    module: str
    state_hash: str
    state_groups: frozenset | None   # training groups of the producing state, read from the module at issue time
    evidence: tuple                  # ((variable, status, digest of that inferred input or None), ...)
    conditional: str | None = None
    training_groups_hash: str | None = None      # the hash the module reported for those groups


class ProvenanceLedger:
    """The inferred values a reconstruction layer has issued, keyed by the digest of the issued record.

    reconstruct() writes one entry per inferred result, and nothing else writes entries (`_issue` is not part of
    the public interface). A later check finds an entry only for a record that is exactly the one that was
    issued: any change of its unit, group, value, scale statement, producing state or lineage gives another
    digest. Several Reconstructors may share one ledger (pass it to their constructors). The ledger lives in
    memory: a record that was saved and reloaded in another process is not in it."""

    def __init__(self):
        self._facts = {}

    def __len__(self):
        return len(self._facts)

    def _issue(self, record, *, state_groups, evidence, conditional=None):
        if record.status != "inferred":
            raise ReconError("only inferred records are issued")
        u = record.uncertainty
        fact = IssuedFact(record_digest(record), record.key(), record.group.group_id, record.variable,
                          str(u.producing_module), str(u.state_hash),
                          None if state_groups is None else frozenset(str(g) for g in state_groups),
                          tuple((str(v), str(s), d) for v, s, d in evidence), conditional,
                          None if u.training_groups_hash is None else str(u.training_groups_hash))
        self._facts[fact.digest] = fact
        return fact

    def get(self, record_sha):
        return self._facts.get(record_sha)

    def fact(self, record):
        """The entry of exactly this record, or None."""
        try:
            return self._facts.get(record_digest(record))
        except (AttributeError, TypeError, ValueError):
            return None

    def ancestry(self, record):
        """Every issued entry in the ancestry of an inferred record, its own first; None when the record or
        one of its inferred ancestors was not issued here (or was changed since)."""
        first = self.fact(record)
        if first is None:
            return None
        out, stack, seen = [], [first], set()
        while stack:
            fact = stack.pop()
            if fact.digest in seen:
                continue
            seen.add(fact.digest)
            out.append(fact)
            for _, status, child in fact.evidence:
                if status == "inferred":
                    child_fact = self._facts.get(child)
                    if child_fact is None:
                        return None
                    stack.append(child_fact)
        return out


def _inferred_input_problem(registry, namespace, unit_key, group_id, record, ledger, resolve, state_inputs,
                            axon_hidden=True):
    """Why an INFERRED record may not feed an answer about a unit whose axon is hidden; None when it may.

    The ancestry is read from the ledger (what the layer recorded when it issued each value), transitively.
    For the record and every inferred ancestor: it was issued here, unchanged, for this unit and group; its
    producing state is registered now with the training groups recorded at issue; those groups do not contain
    the query's group; and, when the answer itself is axon-derived (`axon_hidden`), every measured ancestor is
    one of HIDDEN_AXON_INPUTS and not axon-derived."""
    name = f"inferred {record.variable}"
    if group_id is None:
        return f"{name}: the query's brain sample is unknown, so no producing state can be shown to exclude it"
    if ledger is None:
        return f"{name}: no ledger of issued values is attached, so its provenance cannot be verified"
    first = ledger.fact(record)
    if first is None:
        return (f"{name}: it was not issued by this reconstruction layer, or it was changed after it was issued; "
                "its provenance cannot be verified")
    stack, seen = [first], set()
    while stack:
        fact = stack.pop()
        if fact.digest in seen:
            continue
        seen.add(fact.digest)
        if unit_key is not None and tuple(fact.unit) != tuple(unit_key):
            return f"{name}: its ancestry holds a value inferred for another unit ({fact.unit[2]})"
        if fact.group_id != group_id:
            return (f"{name}: {fact.variable} was inferred for a unit of group {fact.group_id}, the query names "
                    f"group {group_id}")
        registered = resolve(fact.state_hash) if resolve is not None else None
        if registered is None:
            return (f"{name}: the state that produced {fact.variable} is not registered, so it cannot be shown to "
                    "exclude the query's brain sample")
        registered = frozenset(str(g) for g in registered)
        if fact.state_groups is None or registered != fact.state_groups:
            return (f"{name}: the registered training groups of the state that produced {fact.variable} are not "
                    "the ones recorded when the value was issued")
        if fact.training_groups_hash != ids_sha(sorted(registered)):
            return (f"{name}: the training-group hash reported with {fact.variable} is not the hash of the "
                    "registered training groups of its producing state")
        if id_in(group_id, registered):
            return f"{name}: the state that produced {fact.variable} contains the query's brain sample"
        allowed = state_inputs(fact.state_hash) if state_inputs is not None else None
        for ancestor, status, child in fact.evidence:
            if not registry.has_variable(namespace, ancestor):
                return f"{name} descends from {ancestor}, which is not a registered variable of {namespace}"
            ancestor_spec = registry.variable(namespace, ancestor)
            if ancestor_spec.predictor_ban:
                return f"{name} descends from {ancestor}, which may never be a predictor"
            if allowed is not None and ancestor not in allowed:
                return f"{name}: {ancestor} is not a declared input of the registered state"
            if status == "measured":
                if axon_hidden and (ancestor_spec.family == "axon" or ancestor not in HIDDEN_AXON_INPUTS):
                    return (f"{name} descends from a measured {ancestor}, which is not permitted as an input of a "
                            "hidden-axon forecast")
            elif status == "inferred":
                child_fact = ledger.get(child)
                if child_fact is None:
                    return (f"{name}: the ancestry of its inferred input {ancestor} is missing: that value was not "
                            "issued by this layer, or was changed since")
                stack.append(child_fact)
            else:
                return f"{name} descends from a {status} {ancestor}"
    return None


def hidden_axon_problem(registry, namespace, unit_key, group_id, evidence, *, ledger, resolve, state_inputs=None,
                        axon_hidden=True):
    """Why `evidence` may not feed an answer about ONE unit whose axon is hidden; None when every record is
    admissible. ONE rule for every caller.

    `axon_hidden=True`: a HIDDEN-AXON FORECAST. The axon adapter, and every conditional of the general layer
    whose target is axon-derived and that is not declared as taking a measured profile (T1 and every undeclared
    conditional with such a target):
      * a measured record: a registered variable of the namespace, one of HIDDEN_AXON_INPUTS, not axon-derived,
        a record of this unit and group, supplied under its own name, and carrying no inferred provenance;
      * an inferred record, WHATEVER its variable name (a soma position or a source label can be inferred
        from the hidden profile as well as anything else): _inferred_input_problem, transitively;
      * any other status is refused.
    `axon_hidden=False`: every other conditional that is not a declared profile-input task. The same provenance
    and state rule for every inferred record (issued here, for this unit, by registered states that do not hold
    the unit's group, group known); the whitelist of measured inputs is not applied, because the answer is not
    axon-derived.

    `evidence` maps the name a record was supplied under to the record; `resolve(state_hash)` gives the training
    groups of a registered state or None; `state_inputs(state_hash)` optionally gives its declared inputs."""
    if group_id is not None:
        problem = identifier_problem(group_id, "the query's group ID")
        if problem:
            return problem
        implied = None if unit_key is None else registry.implied_group(namespace, unit_key[2])
        if implied is not None and implied != group_id:
            return f"by its ID, unit {unit_key[2]} belongs to group {implied}, not to {group_id}"
    for name, record in evidence.items():
        variable = record.variable
        if variable != name:
            return f"a record of {variable} was supplied under the name {name}"
        if record.status == "reference" and not axon_hidden:
            continue                                    # a population-level reference of the unit's own group
        if not registry.has_variable(namespace, variable):
            return f"{variable} is not a variable of {namespace}"
        spec = registry.variable(namespace, variable)
        if spec.predictor_ban:
            return f"{variable} may never be a predictor: {spec.predictor_ban}"
        if unit_key is not None and record.key() != tuple(unit_key):
            return f"the {variable} record belongs to {record.key()}, not to the query unit {tuple(unit_key)}"
        if group_id is not None and record.group.group_id is not None and record.group.group_id != group_id:
            return f"the {variable} record is filed under group {record.group.group_id}, the query names {group_id}"
        if record.status == "measured":
            if record.derived_from or record.lineage or record.uncertainty.kind == "inferred" \
                    or str(record.source.path).startswith(INFERRED_SOURCE_PREFIX):
                return f"the {variable} record says 'measured' but carries the provenance of an inferred value"
            if axon_hidden and spec.family == "axon":
                return (f"a measured {variable} of the query neuron is never an input to a hidden-axon forecast: "
                        "the neuron's own axon, dendrite and regional profile are for scoring only")
            if axon_hidden and variable not in HIDDEN_AXON_INPUTS:
                return f"{variable} is not a permitted input of a hidden-axon forecast"
        elif record.status == "inferred":
            # EVERY inferred input, including permitted names such as soma_position and source_label
            problem = _inferred_input_problem(registry, namespace, unit_key, group_id, record, ledger, resolve,
                                              state_inputs, axon_hidden)
            if problem:
                return problem
        elif axon_hidden:
            return f"a {record.status} {variable} is not a permitted input of a hidden-axon forecast"
    return None


class Reconstructor:
    """The general layer: checks the request against the registry and the pairing rule, finds the
    module, applies the validation gate and the declared policy, and returns values with provenance
    or structured refusals. Every inferred value it returns is entered in its ledger."""

    def __init__(self, registry, policy=None, ledger=None):
        self.registry = registry
        self.policy = policy or Policy()
        self.ledger = ledger if ledger is not None else ProvenanceLedger()
        self.modules = []

    def register(self, module):
        self.registry.namespace(module.namespace)
        if any(m.name == module.name for m in self.modules):
            raise RegistryError(f"module {module.name} is already registered")
        if module.kind != "reference":
            for variable in module.produces():
                self.registry.variable(module.namespace, variable)
        module.bind(self)
        self.modules.append(module)
        return module

    def module(self, name):
        for module in self.modules:
            if module.name == name:
                return module
        raise RegistryError(f"no module named {name} is registered")

    def describe(self):
        return [m.spec() for m in self.modules]

    def resolve_state(self, state_hash):
        """Training groups of a state held by a registered module, by its hash; None when no module holds it."""
        for module in self.modules:
            groups = module.state_groups(state_hash)
            if groups is not None:
                return frozenset(str(g) for g in groups)
        return None

    # -- request-level checks ------------------------------------------------------------------
    def _check_request(self, request):
        """(refusal or None, the request with its effective group).

        Identifiers must be usable; every record must be the unit's own (pairing rule) or the reference record
        of the unit's own group; every value states the registered scale and frame of its variable; one unit
        belongs to one group, and where the unit ID names its group (Jiao) the stated group must be that one."""
        registry = self.registry

        def refuse(code, reason, detail=()):
            return Refusal(code, reason, tuple(detail)), request

        if request.namespace not in registry.namespaces:
            return refuse("invalid_request", f"namespace {request.namespace} is not registered")
        namespace = registry.namespace(request.namespace)
        if namespace.grain != request.grain:
            return refuse("invalid_request", f"{request.namespace} units are at grain {namespace.grain}, not "
                          f"{request.grain}")
        problem = identifier_problem(request.unit_id, "the unit ID")
        if problem:
            return refuse("invalid_request", problem)
        if request.interval_level is not None and not 0 < request.interval_level < 1:
            return refuse("invalid_request", "interval_level lies strictly between 0 and 1")
        key = (request.namespace, request.grain, request.unit_id)
        seen, groups, own, references = set(), set(), {}, []
        for record in request.observed:
            mine = record.key() == key
            group_reference = (record.status == "reference" and record.namespace == request.namespace)
            if not mine and not group_reference:
                return refuse("unidentifiable_cross_modal",
                              f"observed record {record.variable} belongs to {record.key()}, not to the requested "
                              f"unit {key}: records of different cells, grains or datasets are never joined as paired",
                              (f"{r.contract_from} -> {r.contract_to}: {r.kind}"
                               for r in registry.joins_between(record.namespace, request.namespace)))
            if record.variable in seen:
                return refuse("invalid_request", f"{record.variable} is supplied twice")
            seen.add(record.variable)
            if record.status == "reference":
                references.append(record)
                continue
            if not registry.has_variable(request.namespace, record.variable):
                return refuse("invalid_request", f"{record.variable} is not a variable of {request.namespace}")
            spec = registry.variable(request.namespace, record.variable)
            if spec.predictor_ban and record.status != "missing":
                return refuse("forbidden_input", f"{record.variable} may never be a predictor: {spec.predictor_ban}")
            if record.group.group_id is not None:
                groups.add(record.group.group_id)
            if record.status == "missing":
                continue
            try:
                replace(record)                 # the record's own invariants, checked again on the object as it is now
            except (TypeError, ValueError) as error:
                return refuse("invalid_request", f"{record.variable}: not a valid record ({error})")
            problem = statement_problem(spec, record.scale, record.transform, record.frame)
            if problem:
                return refuse("invalid_request", problem)
            own[record.variable] = record
        stated = request.group.group_id
        named = groups | ({stated} if stated is not None else set())
        if len(named) > 1:
            return refuse("invalid_request", f"the request and its records name the groups {sorted(map(str, named))}: "
                          "one unit belongs to one group")
        effective = next(iter(named)) if named else None
        if effective is not None:
            problem = identifier_problem(effective, "the group ID")
            if problem:
                return refuse("invalid_request", problem)
            implied = registry.implied_group(request.namespace, request.unit_id)
            if implied is not None and implied != effective:
                return refuse("invalid_request", f"by its ID, unit {request.unit_id} belongs to group {implied}, not "
                              f"to {effective}")
        for record in references:               # a reference value attaches only to the unit's OWN group
            link = record.reference
            key_record = own.get(link.key_name)
            if key_record is None or str(key_record.value) != link.key_value \
                    or record.unit_id != f"{link.key_name}={link.key_value}":
                return refuse("unidentifiable_cross_modal",
                              f"the reference record of {record.variable} belongs to the group "
                              f"{link.key_name}={link.key_value}; the unit's own {link.key_name} is "
                              f"{'not supplied' if key_record is None else repr(key_record.value)}: a reference "
                              "value is attached only to the unit's own group")
            if registry.reference_rule(link.from_namespace, request.namespace) is None \
                    or not registry.has_variable(link.from_namespace, record.variable):
                return refuse("unidentifiable_cross_modal",
                              f"no population-level reference alignment of {record.variable} from "
                              f"{link.from_namespace} to {request.namespace} is declared in the data contract")
            problem = statement_problem(registry.variable(link.from_namespace, record.variable), record.scale,
                                        record.transform, record.frame)
            if problem:
                return refuse("invalid_request", problem)
        if stated is None and effective is not None:
            request = replace(request, group=Group(effective, namespace.group_kind, namespace.group_evidence))
        return None, request

    def reconstruct(self, request):
        level = request.interval_level or self.policy.interval_level
        blanket, effective = self._check_request(request)
        evidence = {r.variable: r for r in effective.observed if r.status != "missing"}
        asked = Counter(t.variable for t in request.requested)
        properties = {}
        for target in request.requested:
            if blanket is not None:
                properties[target.variable] = PropertyResult(target.variable, target.coordinates or (), "refused",
                                                             refusal=blanket)
            elif asked[target.variable] > 1:
                properties[target.variable] = self._refused(
                    target, "invalid_request", f"{target.variable} is requested more than once; ask for a variable "
                    "once and list every wanted coordinate there")
            else:
                properties[target.variable] = self._one(effective, target, evidence, level)
        return ReconstructionResult(request, properties, self.policy)

    # -- one requested property ----------------------------------------------------------------
    def _refused(self, target, code, reason, detail=(), **extra):
        return PropertyResult(target.variable, target.coordinates or (), "refused", refusal=Refusal(code, reason, tuple(detail)),
                              **extra)

    def _cross_modal(self, request, target, evidence):
        holders = self.registry.namespaces_measuring(target.variable)
        for module in self.modules:
            if module.kind == "reference" and module.namespace == request.namespace:
                answer = module.reference_answer(request, target, evidence)
                if answer is not None:
                    return answer
        detail = []
        for holder in holders:
            rules = self.registry.joins_between(holder, request.namespace)
            if not rules:
                detail.append(f"{holder}: no linkage to {request.namespace} in the data contract")
            for rule in rules:
                detail.append(f"{rule.contract_from} -> {rule.contract_to}: {rule.kind} ({rule.contract_allowed})")
        return self._refused(target, "unidentifiable_cross_modal",
                             f"{target.variable} is measured only in {', '.join(holders)}, on other units than "
                             f"{request.namespace} {request.grain}s; no verified pairing links them at this grain",
                             detail)

    def _one(self, request, target, evidence, level):
        registry = self.registry
        if not registry.namespaces_measuring(target.variable):
            return self._refused(target, "unsupported_variable", f"{target.variable} is not a registered variable")
        if not registry.has_variable(request.namespace, target.variable):
            return self._cross_modal(request, target, evidence)
        spec = registry.variable(request.namespace, target.variable)
        if target.coordinates is not None and (len(set(target.coordinates)) != len(target.coordinates)
                                               or not set(target.coordinates) <= set(spec.components)):
            return self._refused(target, "invalid_request", f"{target.variable}: unknown or repeated coordinates requested")
        supplied = evidence.get(target.variable)
        if supplied is not None:
            # A vector may be partly observed: then the request is for its other coordinates.
            have = set(supplied.coordinates)
            overlap = (spec.kind != "vector"
                       or (target.coordinates is not None and have & set(target.coordinates))
                       or (target.coordinates is None and have >= set(spec.components)))
            if overlap:
                return self._refused(target, "invalid_request",
                                     f"{target.variable} (or part of what is requested) was supplied as observed")
        for variable, record in evidence.items():
            if record.status == "reference" or not registry.has_variable(request.namespace, variable):
                continue
            if target.variable in registry.variable(request.namespace, variable).derived_from:
                return self._refused(target, "forbidden_input",
                                     f"{variable} is derived from {target.variable} and cannot be evidence for it")
        candidates = [m for m in self.modules if m.kind != "reference" and m.namespace == request.namespace
                      and m.grain == request.grain and target.variable in m.produces()]
        if not candidates:
            return self._refused(target, "unsupported_variable",
                                 f"no registered module produces {target.variable} for {request.namespace} {request.grain}s")
        first_refusal = None
        for module in candidates:
            try:
                plan = module.plan(request, target, evidence, self.policy)
            except ForbiddenInputError as error:
                plan = Refusal("forbidden_input", str(error))
            if isinstance(plan, Refusal):
                first_refusal = first_refusal or (module, plan)
                continue
            return self._gate(request, target, plan, level, evidence)
        module, refusal = first_refusal
        return PropertyResult(target.variable, target.coordinates or (), "refused", refusal=refusal, module=module.name)

    def _inferred_input_notices(self, variable, record, group_id):
        """What the ledger says about an inferred input of a DECLARED profile-input task (there the measured
        profile is an input by declaration, so the input is flagged, not refused)."""
        facts = self.ledger.ancestry(record)
        if facts is None:
            return [Notice("inferred_input_unverified",
                           f"{variable} was not issued by this layer, or it (or one of its ancestors) was changed "
                           "since: its provenance is only what the record says about itself", variable)]
        notices = []
        for fact in facts:
            if fact.state_groups is None:
                notices.append(Notice("inferred_input_unverified",
                                      f"the training groups of the state that produced {fact.variable} are not known "
                                      "to this layer", variable))
            elif id_in(group_id, fact.state_groups):
                notices.append(Notice("group_in_training_state",
                                      f"the state that produced the inferred {fact.variable} holds other units of "
                                      "this unit's group", variable))
        return notices

    def _gate(self, request, target, plan, level, observed):
        registry, module = self.registry, plan.module
        key = (request.namespace, request.grain, request.unit_id)
        group_id = request.group.group_id
        spec = registry.variable(request.namespace, target.variable)
        evidence = (module.evidence_for(plan.conditional.name) if plan.conditional is not None
                    else ValidationEvidence("not_validated", note="this conditional is not declared or benchmarked"))
        provenance = module.provenance(plan)
        common = dict(module=module.name, state_hash=provenance.get("state_hash"),
                      training_groups_hash=provenance.get("training_groups_hash"), validation=evidence.status,
                      conditional=None if plan.conditional is None else plan.conditional.name,
                      benchmark=evidence.benchmark, evidence_used=tuple(plan.evidence_used),
                      evidence_ignored=tuple(plan.evidence_ignored), unit=key, group=request.group)
        notices = list(plan.notices)
        for variable, status in plan.evidence_used:
            if status in ("inferred", "reference"):
                notices.append(Notice(f"conditioned_on_{status}_input",
                                      f"{variable} is {status}, not measured; its uncertainty is not propagated and "
                                      "the interval coverage is not validated for such inputs", variable))
        if plan.group_in_training:
            notices.append(Notice("group_in_training_state",
                                  "other units of this unit's group are in the module state; the benchmark evidence "
                                  "refers to groups outside the state only"))
        elif group_id is None and provenance.get("training_groups_hash"):
            notices.append(Notice("group_unknown",
                                  "no group was given for this unit, so it cannot be shown to lie outside the groups "
                                  "the module state was fitted on"))
        benchmarked = ((evidence.coverage or {}).get("engine") or {}).get("level")
        if benchmarked is not None and abs(float(benchmarked) - level) > 1e-12:
            notices.append(Notice("interval_level_not_benchmarked",
                                  f"coverage was checked at level {benchmarked}, not at the requested {level}"))
        declared = plan.conditional is not None and plan.conditional.declared_measured_axon_input
        used_inferred = [variable for variable, status in plan.evidence_used if status == "inferred"]

        def refuse(code, reason):
            return self._refused(target, code, reason, notices=tuple(notices), **common)

        def needs_declaration(variable):
            return (registry.has_variable(request.namespace, variable)
                    and registry.variable(request.namespace, variable).measured_input_needs_declared_task)

        # 1. a measured axon-derived input is an input only in tasks declared as taking it
        for variable, status in plan.evidence_used:
            if needs_declaration(variable) and status == "measured" and not declared:
                return refuse("forbidden_input",
                              f"a measured {variable} is an input only in tasks declared as taking it; it is "
                              "never an input to a hidden-axon forecast or to an undeclared conditional")
        # 2. an inferred input brings its whole ancestry: a permitted name does not launder what it came from
        for variable in used_inferred:
            for step in observed[variable].lineage:
                for ancestor, ancestor_status in step.evidence:
                    if ancestor == target.variable:
                        return refuse("forbidden_input",
                                      f"the inferred {variable} was itself derived from {target.variable}: "
                                      "the requested property cannot be among the ancestors of its evidence")
                    if needs_declaration(ancestor) and ancestor_status == "measured" and not declared:
                        return refuse("forbidden_input",
                                      f"the inferred {variable} descends from a measured {ancestor}, which is "
                                      "an input only in tasks declared as taking it")
        # 3. The unit's axon is hidden in every request that is not one of the declared profile-input tasks.
        #    ONE rule then covers every record of the request, used by the module or not: an inferred input,
        #    whatever its name, must be a value this layer issued for this unit, by registered states that do
        #    not hold the unit's group; and when the answer itself is axon-derived (a hidden-axon FORECAST), the
        #    measured inputs and the measured ancestors are held to the whitelist.
        if not declared:
            problem = hidden_axon_problem(registry, request.namespace, key, group_id, observed, ledger=self.ledger,
                                          resolve=self.resolve_state, axon_hidden=spec.family == "axon")
            if problem:
                return refuse("forbidden_input", problem)
        else:
            for variable in used_inferred:
                for notice in self._inferred_input_notices(variable, observed[variable], group_id):
                    if notice not in notices:
                        notices.append(notice)
        if plan.unit_in_training:
            return refuse("leakage_guard", "this unit is part of the module state; use a state that does not contain it")
        if evidence.status.startswith(BLOCKED_PREFIX):
            return refuse("not_validated", evidence.status)
        degraded = any(n.code in ("out_of_support", "conditioned_on_inferred_input", "conditioned_on_reference_input",
                                  "group_in_training_state", "group_unknown", "interval_level_not_benchmarked",
                                  "evidence_not_linked", "inferred_input_unverified")
                       for n in notices)

        def inferred(estimate, *, coverage, baseline_level=False, training_groups=None, **override):
            fields = dict(common)
            fields.update(override)
            # every answer states the registered scale and frame of its variable, in the planned coordinates
            problem = statement_problem(spec, estimate.scale, estimate.transform, estimate.frame)
            if problem is None and tuple(estimate.coordinates) != tuple(plan.coordinates):
                problem = f"{target.variable}: the answer's coordinates are not the planned ones"
            if problem:
                return refuse("invalid_estimate", f"the answer of {fields['module']} was withheld: {problem}")
            # this inference first, then the whole lineage of every inferred input: nothing is erased
            lineage = [LineageStep(target.variable, str(fields["module"]), str(fields["state_hash"]),
                                   str(fields["training_groups_hash"]), tuple(plan.evidence_used), key, group_id)]
            for variable in used_inferred:
                lineage.extend(observed[variable].lineage)
            result = PropertyResult(target.variable, tuple(plan.coordinates), "inferred", estimate=estimate,
                                    coverage_validated=bool(coverage) and not degraded, notices=tuple(notices),
                                    baseline_level=baseline_level, lineage=tuple(lineage), **fields)
            try:
                record = result.to_observation(request)
            except (TypeError, ValueError) as error:
                return refuse("invalid_estimate", f"the answer of {fields['module']} was withheld: it does not make "
                              f"a valid record ({error})")
            self.ledger._issue(record, state_groups=training_groups, conditional=fields["conditional"],
                              evidence=tuple((variable, status,
                                              record_digest(observed[variable]) if status == "inferred" else None)
                                             for variable, status in plan.evidence_used))
            return result

        def baseline_answer(code, reason):
            answer = module.baseline(plan, level, evidence)
            if answer is None:
                return refuse(code, reason + "; no baseline-level answer is available")
            notices.append(Notice(code, reason, handling=f"baseline-level answer from {answer.module}"))
            return inferred(answer.estimate, coverage=answer.coverage_validated, baseline_level=True,
                            training_groups=answer.training_groups, module=answer.module,
                            state_hash=answer.state_hash, training_groups_hash=answer.training_groups_hash)

        engine_groups = provenance.get("training_groups")
        if evidence.status == "exploratory_supported":
            return inferred(module.estimate(plan, level), coverage=evidence.coverage_validated,
                            training_groups=engine_groups)
        if evidence.status == "no_better_than_baseline":
            reason = ("the benchmark of this conditional did not beat the simple baselines"
                      + (f" (best comparator: {evidence.best_comparator})" if evidence.best_comparator else ""))
            if self.policy.on_no_better_than_baseline == "refuse":
                return refuse("no_better_than_baseline", reason)
            return baseline_answer("no_better_than_baseline", reason)
        reason = evidence.note or "no benchmark evidence is registered for this conditional"
        if self.policy.on_not_validated == "refuse":
            return refuse("not_validated", reason)
        if self.policy.on_not_validated == "baseline":
            return baseline_answer("not_validated", reason)
        notices.append(Notice("not_validated", reason, handling="engine answer returned by policy, flagged"))
        return inferred(module.estimate(plan, level), coverage=False, training_groups=engine_groups)
