"""SYN4RCH axon prediction engine (library).

Task: predict the 3D axon of a hidden hypothalamic neuron (Jiao et al. 2025) from what is known
about it, and decide honestly whether that beats copying a measured neuron. The binding rules are
in reports/axon_dataset_plan.md and its amendments 1 and 2 (cohort, sealing, inclusion, no clipping; a
neuron the source does not serve; no test trace on disk before the freeze; this engine downloads nothing), in
reports/prediction_independent_evaluation_contract.md (the independent referee) and in
reports/axon_engine_protocol.md (candidates, controls, gate, verdict). Commands live in
analysis/axon_engine_cli.py.

What this module does and does not do on import: it verifies the pinned SHA-256 of the other
agent's frozen prototype (analysis/brm_template.py and its dependencies) against both the freeze
record data/interim/brm_template_method/frozen.json and the constants below, then imports those
modules read-only. It reads no trace and no regional row on import. The independent referee
(analysis/prediction_independent_evaluator.py) is loaded read-only on first use and never edited.

Evidence boundary (plan sections 7 and 8):
  * a forecast takes only the query's table soma (mirrored-left frame), source label, soma label,
    hemisphere and reconstruction type, plus states and donors of training animals;
  * regional rows of validation and test animals are skipped unparsed by the only loader;
  * held-out traces are opened first by `admit`, after the complete forecast manifest is committed
    and bound to the current method; admission applies the mechanical inclusion rules and reads no
    forecast; scoring follows (validation: this engine's scorer; test: the referee's score_once).

Donor pool of every method: the whole training library (for a training query in development, minus
its own animal). Pair costs come from a pinned all-pairs cache of the frozen metric; no method uses
a sampled panel.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict, deque
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "axon_engine_v1"

# ----------------------------------------------------------------------------------------------
# hashing and JSON
# ----------------------------------------------------------------------------------------------


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    """Same canonical JSON as the frozen code (sorted keys, compact, newline)."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"


def digest(value):
    return sha_bytes(canonical(value).encode())


def ids_sha(ids):
    return sha_bytes("\n".join(ids).encode())


def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value, *, exclusive=False, pretty=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    value = jsonable(value)
    text = json.dumps(value, sort_keys=True, indent=1, allow_nan=False) + "\n" if pretty else canonical(value)
    with path.open("x" if exclusive else "w") as handle:
        handle.write(text)
    return path


def utc_now():
    return datetime.now(timezone.utc).isoformat()


# ----------------------------------------------------------------------------------------------
# errors
# ----------------------------------------------------------------------------------------------


class EngineError(RuntimeError):
    pass


class FrozenSourceChanged(EngineError):
    """A pinned file of the frozen prototype differs from its recorded SHA-256."""


class SealError(EngineError):
    """Sealed data (a held-out trace or regional row) was requested outside its stage."""


class LeakError(EngineError):
    """A query animal (or a held-out animal) is in a regional state or a donor set."""


class StageBlocked(EngineError):
    """The stage cannot run on the data that are on disk; the message states the blocker."""


class OnceOnly(EngineError):
    """A once-only step was requested a second time."""


class FreezeMismatch(EngineError):
    """Code, data or parameters differ from the freeze record."""


class Exclusion(Exception):
    """A mechanical inclusion rule of plan section 9 is violated."""

    def __init__(self, rule, detail):
        super().__init__(f"rule {rule}: {detail}")
        self.rule, self.detail = rule, detail


# ----------------------------------------------------------------------------------------------
# the frozen prototype: pinned, imported read-only
# ----------------------------------------------------------------------------------------------

FREEZE_RECORD = "data/interim/brm_template_method/frozen.json"
FREEZE_RECORD_SHA = "aa71c5a618686810343110de059a716df0aa70797466f69c715b10f944684049"
PINNED_SOURCES = {
    "analysis/brm_template.py": "0a3da130c0e0049e8c9cc5df94523442f542f44a6815e82d3ad947a960a8d49a",
    "analysis/brm_geometry.py": "8de317bc4f261e8061b0bccad1dbba537b2a8daf6f6834d94d303645d6611187",
    "analysis/brm_joint.py": "ae7f5ec432bcbd2b0eebb65595549bc12158b1b7c047f5ace167c7840c0c937e",
    "analysis/brm_real.py": "5225ea5833b32fee57e735600eda87a9784cb2b1798cd326ee4203e1b6997648",
    "analysis/d110_geometry.py": "3b011643b1b224256fdce1ad8a31014e01e644c6a794bd1d49af5e6543090a5e",
}


def verify_frozen_sources(root=ROOT):
    """Every file pinned by the other agent's freeze record must be unchanged; the five modules
    this engine imports must in addition match the constants written in this file."""
    root = Path(root)
    record_path = root / FREEZE_RECORD
    record_sha = sha_file(record_path)
    if record_sha != FREEZE_RECORD_SHA:
        raise FrozenSourceChanged(f"{FREEZE_RECORD} changed: {record_sha}")
    record = read_json(record_path)
    pins = record["source_sha256"]
    for relative, expected in PINNED_SOURCES.items():
        if pins.get(relative) != expected:
            raise FrozenSourceChanged(f"freeze record and engine disagree on the pin of {relative}")
    if record.get("algorithm_source_sha256") != PINNED_SOURCES["analysis/brm_template.py"]:
        raise FrozenSourceChanged("freeze record names another algorithm source")
    changed = [relative for relative, expected in sorted(pins.items()) if sha_file(root / relative) != expected]
    if changed:
        raise FrozenSourceChanged("pinned frozen files changed: " + ", ".join(changed))
    return {"freeze_record": FREEZE_RECORD, "freeze_record_sha256": record_sha, "source_sha256": dict(sorted(pins.items()))}


FROZEN_PINS = verify_frozen_sources()


def _load_frozen_template():
    name = "axon_engine_frozen_brm_template"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / "analysis/brm_template.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


T = _load_frozen_template()      # analysis/brm_template.py   (five policies, validation gate)
Geo = T.Geo                      # analysis/brm_geometry.py   (Query, RegionalState, state_posterior, assert_query_safe)
G = T.G                          # analysis/d110_geometry.py  (Arbor, resample, geometry_metrics, training_medoid, save_arbor)
R = Geo.R                        # analysis/brm_real.py       (fit_joint, fit_outer_fold, snapshot_model)
REGIONS = tuple(R.REGIONS)

# ----------------------------------------------------------------------------------------------
# declared configuration (pinned by `freeze`)
# ----------------------------------------------------------------------------------------------

PLAN_MANIFEST_SHA = "649f829a592cd3d13bcda5cffdbd049e759568d38bacb144f6a61220610a1ca4"
PLAN_ANIMALS_SHA = "11a0680bfcb53b204e7224a1f982174bbf11b25fa6aedb3e9476ed465f69df9f"
SOMA_TABLE_SHA = "cc01d048bf5a7c9c9f46869fb750f2d67a98ca934903d2e2e71c1a67b1d97cb9"
REGIONAL_TABLE_SHA = "2072bc20550607c195be0a42cee9fcbe0e0ab5aaebfaee0a1016b2a43d35e63c"

MIDLINE_ML = 5700.0
MIRROR_ML = 11400.0
BOX = np.array([13200.0, 8000.0, 11400.0])
STEP = 50.0
ALLOWED_TYPES = (0, 1, 2, 3, 4, 10, 11, 12)
DENDRITE_TYPES = (3, 4)
MARKER_TYPES = (10, 11, 12)
STAGES = ("develop", "validate", "test")
PART_OF_STAGE = {"develop": "train", "validate": "validation", "test": "test"}


# Copying controls. The four official ones are the referee contract's (exact names); the engine obtains them
# from the referee's own copying_baselines and keeps an independent replica that must agree with it.
OFFICIAL_COPIES = ("nearest", "medoid", "peptide_nearest", "peptide_medoid")
FROZEN_CONTROLS = ("C1_medoid_copy", "C4_uniform_attach_500", "C5_uniform_attach_1000")
CONTROLS = OFFICIAL_COPIES + FROZEN_CONTROLS
COPY_BASELINES = OFFICIAL_COPIES + ("C1_medoid_copy",)      # one measured neuron, translated (referee role copy_baselines)
ATTACHMENT_CONTROLS = ("C4_uniform_attach_500", "C5_uniform_attach_1000")     # referee role matched_controls, always
FROZEN_CANDIDATES = ("F1_conditional_attach_500", "F2_conditional_attach_1000")
ENGINE_CANDIDATES = ("E1_local_uniform", "E2_local_kernel", "E3_local_kernel_label",
                     "E4_local_kernel_regional", "E5_local_kernel_label_regional", "E6_local_consensus")
CANDIDATES = FROZEN_CANDIDATES + ENGINE_CANDIDATES
METHODS = CONTROLS + CANDIDATES
FROZEN_POLICY_OF = {"C1_medoid_copy": "baseline", "C4_uniform_attach_500": "uniform_anchor_500",
                    "F1_conditional_attach_500": "conditional_anchor_500",
                    "C5_uniform_attach_1000": "uniform_anchor_1000",
                    "F2_conditional_attach_1000": "conditional_anchor_1000"}
METHOD_OF_FROZEN_POLICY = {policy: name for name, policy in FROZEN_POLICY_OF.items()}
# Conditional candidates (label- or regional-conditional weights) and the controls they must also beat.
MATCHED_CONTROLS = {
    "F1_conditional_attach_500": ("C4_uniform_attach_500",),
    "F2_conditional_attach_1000": ("C5_uniform_attach_1000",),
    "E3_local_kernel_label": ("E1_local_uniform", "E2_local_kernel"),
    "E4_local_kernel_regional": ("E1_local_uniform", "E2_local_kernel"),
    "E5_local_kernel_label_regional": ("E1_local_uniform", "E2_local_kernel"),
}
# Methods whose forecast arrays depend on the numbers of a fitted regional state. The frozen baseline and
# uniform policies consult the state for membership checks only: their arrays do not depend on it (tested).
USES_FITTED_STATE = ("F1_conditional_attach_500", "F2_conditional_attach_1000",
                     "E4_local_kernel_regional", "E5_local_kernel_label_regional")
FAMILIES = {"all": CANDIDATES, "no_fitted_state": tuple(c for c in CANDIDATES if c not in USES_FITTED_STATE)}
USER_LABEL = {"C1_medoid_copy": "C1", "medoid": "C1 by the referee's rule", "nearest": "C2", "peptide_nearest": "C3",
              "C4_uniform_attach_500": "C4", "C5_uniform_attach_1000": "C5", "peptide_medoid": "same-label representative"}
METHOD_TEXT = {
    "nearest": "official copy (referee contract): the training neuron with the nearest table soma, translated to the query soma (the user's C2)",
    "medoid": "official copy (referee contract): the training neuron minimising the animal-balanced mean soma-centred geometric distance to the whole training pool, translated (the user's C1 by the referee's arithmetic)",
    "peptide_nearest": "official copy (referee contract): nearest among training neurons with the query's source label; falls back to 'nearest' when the label is missing or has no donor (the user's C3)",
    "peptide_medoid": "official copy (referee contract): medoid among training neurons with the query's source label; falls back to 'medoid'",
    "C1_medoid_copy": "frozen prototype's baseline: animal-balanced geometric medoid of the whole training pool by the frozen function, translated (the user's C1; equals 'medoid' except at floating-point ties)",
    "C4_uniform_attach_500": "frozen uniform-attachment control, 500 um, whole training pool",
    "C5_uniform_attach_1000": "frozen uniform-attachment control, 1,000 um, whole training pool",
    "F1_conditional_attach_500": "frozen prototype: conditional medoid of the whole training pool, 500 um local attachment",
    "F2_conditional_attach_1000": "frozen prototype: conditional medoid of the whole training pool, 1,000 um local attachment",
    "E1_local_uniform": "local medoid of the K nearest donors, uniform animal-balanced weights (matched uniform control of E3 to E5)",
    "E2_local_kernel": "local medoid, weights from an adaptive kernel on soma distance",
    "E3_local_kernel_label": "local medoid, kernel times a source-label factor",
    "E4_local_kernel_regional": "local medoid, kernel times the frozen conditional (regional) weights",
    "E5_local_kernel_label_regional": "local medoid, kernel times label factor times frozen conditional weights",
    "E6_local_consensus": "the E2 template with the cable that fewer than a quarter of its (kernel-weighted) neighbours share within 200 um pruned away: not a copy of any single neuron",
}


def method_role(model):
    return "official_copy" if model in OFFICIAL_COPIES else ("frozen_control" if model in FROZEN_CONTROLS else "candidate")


CONFIG = {
    "schema": SCHEMA,
    "frame": "mirrored-left (ML' = 11400 - ML for right-hemisphere neurons, whole neuron); CCFv3 AP DV ML in um",
    "step_um": STEP,
    "known_query_soma": "paper soma table, mirrored-left frame; a query trace is never opened for prediction",
    "inclusion": {
        "rules": "plan section 9, rules 1 to 6, and rule 7 as replaced by amendment 1 (no clipping)",
        "soma_tolerance_um": 5.0, "min_axon_cable_um": 1000.0,
        "allowed_types": list(ALLOWED_TYPES), "dendrite_types": list(DENDRITE_TYPES),
        "marker_types": list(MARKER_TYPES),
        "axon_rule": "all cable that is not dendrite and not the soma; a marker takes the type of its nearest non-marker ancestor",
        "terminal_marker_twigs": "kept as axon cable (plan rule; the frozen extractor drops them)",
        "non_axon_connectors": ("ancestors of axon nodes that are not axon (a marker inside a path; a dendritic stem that an axon "
                                "leaves) are kept as connector geometry of the scored tree and their length is logged; rule 6 "
                                "counts axon cable only"),
        "atlas_box_um": BOX.tolist(),
        "clipping": ("none (plan amendment 1): truth, donors and forecasts are complete trees in the same full coordinate "
                     "domain; the length of cable outside the atlas box is logged per neuron and reported per part, and is "
                     "never removed, moved or reconnected"),
        "copy_preference": "true-hemisphere copy first, then mirrored-left copy",
        "unavailable_trace_rule": ("rule 8 (plan amendment 2): a neuron the download listing records as unavailable_at_source is "
                                   "excluded, logged and counted per part and per source label; decided by the source alone"),
        "test_traces": "no test trace may be listed on disk before the freeze (plan amendment 2); the real freeze refuses otherwise",
    },
    "metrics": {"distance": "distance_um (symmetric arc-weighted nearest-sample distance, frozen geometry_metrics)",
                "f1": "f1_200", "aggregation": "per neuron, mean within animal, mean over animals"},
    "donor_pools": {
        "library": "every training neuron whose trace passes the inclusion rules; for a training query (development) minus its own animal",
        "nearest, peptide_nearest": "the whole library (peptide_nearest: its neurons with the query's source label)",
        "medoid, peptide_medoid": "the whole library, all pairs (peptide_medoid: its neurons with the query's source label)",
        "C1, C4, C5, F1, F2": "the whole library, handed to the unmodified frozen predict_policy; pair costs read from the pinned all-pairs cache",
        "E1 to E6": "the K nearest neurons of the whole library by table-soma distance",
        "sampled panels": "none",
    },
    "pair_cache": {"what": "soma-centred distance_um and 1 - f1_200 between every two library arbors, frozen geometry_metrics, lower neuron ID first",
                   "storage": "dense float64 matrices in the order of the sorted library IDs, with the ID list, the arbor hashes and the metric code hashes",
                   "block_pairs": 2000, "spot_check_pairs": 24},
    "official_copies": {
        "names": list(OFFICIAL_COPIES),
        "source": ("the referee's own copying_baselines(training trees, known query metadata); its arrays are the official four "
                   "controls. The engine's replica (its own code on the pinned cache, same arithmetic as the contract) must give "
                   "the same donors and arrays"),
        "heldout_inline_api_pair_limit": 250000,
        "above_the_limit": "forecasts carry the replica; `copies-check` must then show exact equality with copying_baselines before any truth is opened",
        "develop_api_pair_budget": 80000,
        "develop": "replica for every query; copying_baselines is called for the first held-out animals (ID order) within the pair budget",
    },
    "out_of_fold": {"folds": 5, "salt": "SYN4RCH-axon-engine-oof-2026-11-10",
                    "rule": "training animals dealt to folds in SHA-256 order inside each source label"},
    "regional_state": "frozen brm_real.fit_outer_fold (K grid and inner folds of the frozen code), training animals only",
    "engine": {"k_neighbours": 16, "kernel_h0_um": 100.0, "kernel_rank": 8, "label_other_factor": 0.25,
               "attachment_ladder_um": [500.0, 1000.0], "attachment_last_resort": "translation of the same template",
               "consensus_radius_um": 200.0, "consensus_min_support": 0.25, "consensus_min_cable_um": 1000.0,
               "neighbour_distance": "Euclidean between table somas; ties by neuron ID"},
    "controls": list(CONTROLS), "official_copy_controls": list(OFFICIAL_COPIES), "frozen_controls": list(FROZEN_CONTROLS),
    "frozen_candidates": list(FROZEN_CANDIDATES), "engine_candidates": list(ENGINE_CANDIDATES),
    "matched_controls": {k: list(v) for k, v in MATCHED_CONTROLS.items()},
    "uses_fitted_state": list(USES_FITTED_STATE), "families": {k: list(v) for k, v in FAMILIES.items()},
    "gate": ("a candidate of the declared family is eligible when its animal-mean distance is strictly lower AND its "
             "animal-mean F1 strictly higher than EVERY control (the four official copies, the frozen baseline and both "
             "frozen uniform-attachment controls) on the whole admitted validation set; a conditional candidate must also "
             "beat its matched controls on both; a frozen conditional policy is admitted only when the frozen gate itself "
             "selected it; among eligible candidates lowest distance, then highest F1, then the fixed order; none eligible "
             "means 'no method beyond copying'"),
    "referee_roles": {"candidate": "the frozen selection (or, when the gate selected nothing, the predeclared descriptive reference)",
                      "copy_baselines": list(COPY_BASELINES),
                      "matched_controls": "C4 and C5 always, plus the matched controls of a conditional candidate",
                      "success": "beats every copy baseline AND every matched control on both animal-mean measures"},
    "validation_scoring": ("this engine's scorer on the admitted validation queries, all methods; every score is computed by the "
                           "frozen metric and by the referee's metric module and must agree bit for bit; scored once"),
    "test_scoring": "the referee's freeze_declaration, copying_baselines, seal_predictions and score_once; scored once",
    "subgroups": {"always": "per animal", "by": ["source_label", "soma_label"], "sparse_below_animals": 5,
                  "interpretation": "every subgroup row is descriptive; rows with fewer than 5 animals are marked sparse",
                  "unassigned_soma_labels": "HY and fibre-tract labels (lower-case initial) are shown as unassigned",
                  "counts": "every row shows neurons and animals (brain sample IDs)"},
    "unit_of_independence": "animal = ID prefix = the portal's brain sample ID (166 IDs; the paper reports 166 mice)",
}
CONFIG_SHA = digest(CONFIG)
SUBGROUP_SPEC_SHA = digest(CONFIG["subgroups"])

ENGINE_FILES = ("analysis/axon_engine.py", "analysis/axon_engine_cli.py")
REFEREE_FILES = {"evaluator": "analysis/prediction_independent_evaluator.py", "geometry_math": "analysis/brm_geometry_audit.py"}
REFEREE_CONTRACT = "reports/prediction_independent_evaluation_contract.md"
REFEREE_FINDING = "reports/prediction_referee_interop_finding.md"


def engine_code_sha(root=ROOT):
    return {relative: sha_file(Path(root) / relative) for relative in ENGINE_FILES if (Path(root) / relative).exists()}


def code_hashes(root=ROOT):
    return {"engine": engine_code_sha(root), "frozen": dict(FROZEN_PINS["source_sha256"]),
            "freeze_record_sha256": FROZEN_PINS["freeze_record_sha256"]}


def referee_hashes(root=ROOT):
    """Recorded with every manifest; the referee is not part of the method and is never edited here."""
    return {name: sha_file(Path(root) / relative) for name, relative in REFEREE_FILES.items()}


_REFEREE = {}


def referee():
    """The independent referee, imported read-only from its own file on first use."""
    module = _REFEREE.get("module")
    if module is None:
        name = "axon_engine_independent_referee"
        spec = importlib.util.spec_from_file_location(name, ROOT / REFEREE_FILES["evaluator"])
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        if module.source_hashes() != referee_hashes():
            raise EngineError("the referee's own source hashes differ from the files on disk")
        if tuple(module.COPIES) != OFFICIAL_COPIES:
            raise EngineError(f"the referee's mandatory copies are {module.COPIES}, the engine declares {OFFICIAL_COPIES}")
        _REFEREE["module"] = module
        _REFEREE["sha256"] = referee_hashes()
    elif _REFEREE["sha256"] != referee_hashes():
        raise EngineError("the referee's files changed while this process was running")
    return module


def arbor_sha(arbor):
    """Content hash of the three arrays (reproducible; an NPZ file hash is not, because of zip timestamps)."""
    h = hashlib.sha256()
    for array, kind in ((arbor.xyz, "<f8"), (arbor.parent, "<i8"), (arbor.types, "<i8")):
        data = np.ascontiguousarray(np.asarray(array).astype(kind, copy=False))
        h.update(str(data.shape).encode())
        h.update(data.tobytes())
    return h.hexdigest()


def same_arbor(a, b):
    return (np.array_equal(a.xyz, b.xyz) and np.array_equal(a.parent, b.parent) and np.array_equal(a.types, b.types))


def as_tree(arbor):
    """The referee's tree format: (xyz_um, parent, types)."""
    return (arbor.xyz, arbor.parent, arbor.types)


def same_tree(arbor, tree):
    return (np.array_equal(arbor.xyz, tree[0]) and np.array_equal(arbor.parent, tree[1]) and np.array_equal(arbor.types, tree[2]))



# ----------------------------------------------------------------------------------------------
# world: where the data live
# ----------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class World:
    name: str
    kind: str                    # "real" or "synthetic"
    manifest: str                # one row per neuron: part, labels (plan manifest)
    local_manifest: str          # the same plus hemisphere and reconstruction type
    soma_table: str              # file, neurons, geno, location, AP, DV, ML (mirrored-left frame)
    regional_table: str          # axon length in 49 regions per neuron; held-out rows are sealed
    trace_manifests: tuple       # ((manifest path, base directory, layout), ...)
    out: str                     # engine outputs
    pins: tuple = ()             # (("manifest", sha256), ...)
    expected: tuple = ()         # (("train", animals, neurons), ...)
    protocol: str = ""
    development_snapshot: bool = False   # a fixed partial training library (local37): `develop` may run on it

    def to_json(self):
        return {"name": self.name, "kind": self.kind, "manifest": self.manifest, "local_manifest": self.local_manifest,
                "soma_table": self.soma_table, "regional_table": self.regional_table,
                "trace_manifests": [list(x) for x in self.trace_manifests], "out": self.out,
                "pins": [list(x) for x in self.pins], "expected": [list(x) for x in self.expected],
                "protocol": self.protocol, "development_snapshot": self.development_snapshot}

    @staticmethod
    def from_json(value):
        return World(value["name"], value["kind"], value["manifest"], value["local_manifest"], value["soma_table"],
                     value["regional_table"], tuple(tuple(x) for x in value["trace_manifests"]), value["out"],
                     tuple(tuple(x) for x in value["pins"]), tuple(tuple(x) for x in value["expected"]),
                     value.get("protocol", ""), bool(value.get("development_snapshot", False)))


def real_world(root=ROOT, out=None, legacy_only=False):
    """The locked cohort. With legacy_only the trace store is the fixed set of 37 training traces under
    data/raw/jiao2025/swc (what was on disk when the engine was built): the development snapshot `local37`.
    It ignores data/raw/jiao2025/traces altogether, so a fetch in progress cannot change what it sees."""
    root = Path(root)
    j = root / "data/raw/jiao2025"
    legacy = (str(j / "swc/MANIFEST.tsv"), str(j / "swc"), "legacy_training")
    return World(
        name="jiao2025_local37" if legacy_only else "jiao2025", kind="real",
        manifest=str(root / "reports/axon_dataset_manifest.tsv"),
        local_manifest=str(root / "data/interim/axon_dataset/local_manifest.tsv"),
        soma_table=str(j / "All_soma_infor_HY_neurons.csv"),
        regional_table=str(j / "all_neurons_projection_length_fig1.csv"),
        trace_manifests=(legacy,) if legacy_only else (legacy, (str(j / "traces/MANIFEST.tsv"), str(j / "traces"), "by_part")),
        out=str(out or root / ("data/interim/axon_engine/real_local37" if legacy_only else "data/interim/axon_engine/real")),
        pins=(("manifest", PLAN_MANIFEST_SHA), ("soma_table", SOMA_TABLE_SHA), ("regional_table", REGIONAL_TABLE_SHA)),
        expected=(("train", 101, 4448), ("validation", 25, 1185), ("test", 40, 1547)),
        protocol=str(root / "reports/axon_engine_protocol.md"), development_snapshot=bool(legacy_only))


def synthetic_world(directory, out=None):
    """A generated world keeps the real layout: manifest, soma table, regional table, traces by part."""
    directory = Path(directory)
    data = directory / "data"
    described = read_json(data / "world.json")
    return World(
        name=described["name"], kind="synthetic",
        manifest=str(data / "manifest.tsv"), local_manifest=str(data / "manifest.tsv"),
        soma_table=str(data / "soma.csv"), regional_table=str(data / "regional.csv"),
        trace_manifests=((str(data / "traces/MANIFEST.tsv"), str(data / "traces"), "by_part"),),
        out=str(out or directory / "engine"), protocol=str(ROOT / "reports/axon_engine_protocol.md"))


# ----------------------------------------------------------------------------------------------
# cohort (metadata and permitted inputs only)
# ----------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class NeuronRow:
    neuron: str
    animal: str
    full_id: str
    source_label: str     # manifest label (one per animal)
    tag: str              # the soma table's own spelling, as the frozen regional model was built on it
    soma_label: str       # fine soma label of the manifest
    hemisphere: str
    recon_type: str
    part: str
    klass: str
    soma: tuple           # table soma, mirrored-left frame, um (AP, DV, ML)


SOMA_COLUMNS = ("file", "neurons", "geno", "location", "AP", "DV", "ML")   # never 'type'


def soma_label_kind(label):
    """'nucleus', or why the label assigns no nucleus."""
    if label is None or label == "":
        return "unknown"
    if label == "HY":
        return "unassigned_HY"
    if label[:1].islower():
        return "unassigned_fibre_tract"
    return "nucleus"


class Cohort:
    def __init__(self, world):
        self.world = world
        self.sha = {"manifest": sha_file(world.manifest), "local_manifest": sha_file(world.local_manifest),
                    "soma_table": sha_file(world.soma_table)}
        pins = dict(world.pins)
        for key in ("manifest", "soma_table"):
            if key in pins and self.sha[key] != pins[key]:
                raise EngineError(f"{key} differs from the pinned SHA-256 of the locked plan")
        with open(world.manifest, newline="") as handle:
            plan_rows = list(csv.DictReader(handle, delimiter="\t"))
        with open(world.local_manifest, newline="") as handle:
            local = {r["neuron"]: r for r in csv.DictReader(handle, delimiter="\t")}
        soma = {}
        with open(world.soma_table, newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader)
            index = {name: header.index(name) for name in SOMA_COLUMNS}
            for row in reader:
                neuron = row[index["neurons"]]
                if neuron in soma:
                    raise EngineError(f"duplicate soma row {neuron}")
                soma[neuron] = (row[index["file"]], row[index["geno"]],
                                (float(row[index["AP"]]), float(row[index["DV"]]), float(row[index["ML"]])))
        if len(local) != len(plan_rows) or set(soma) != set(local):
            raise EngineError("manifest, local manifest and soma table do not list the same neurons")
        rows, label_of, part_of, never_heldout = {}, {}, {}, set()
        for r in plan_rows:
            neuron, animal = r["neuron"], r["animal"]
            if any(r.get(flag) == "1" for flag in ("animal_inspected", "animal_development", "animal_in_li_2024")):
                never_heldout.add(animal)      # plan section 4, rule 1: these animals are training only
            l = local[neuron]
            for key in ("animal", "source_label", "soma_label", "part"):
                if l[key] != r[key]:
                    raise EngineError(f"local manifest differs from the plan manifest for {neuron}: {key}")
            full_id, tag, xyz = soma[neuron]
            if animal != neuron.split("_", 1)[0] or not full_id.startswith(neuron + "_"):
                raise EngineError(f"identity fields disagree for {neuron}")
            if r["part"] not in ("train", "validation", "test"):
                raise EngineError(f"unknown part for {neuron}")
            if ("Orexin" if tag == "Wt(orexin)" else tag) != r["source_label"]:
                raise EngineError(f"table tag and manifest label disagree for {neuron}")
            if label_of.setdefault(animal, r["source_label"]) != r["source_label"] or part_of.setdefault(animal, r["part"]) != r["part"]:
                raise EngineError(f"animal {animal} has more than one source label or part")
            if not all(math.isfinite(v) for v in xyz) or xyz[2] > MIDLINE_ML + 1e-6:
                raise EngineError(f"table soma of {neuron} is not in the mirrored-left frame")
            rows[neuron] = NeuronRow(neuron, animal, full_id, r["source_label"], tag, r["soma_label"], l["hemisphere"],
                                     l["reconstruction_type"], r["part"], r.get("class", ""), xyz)
        self.rows = dict(sorted(rows.items()))
        self.order = list(self.rows)
        self.animal_label, self.animal_part = label_of, part_of
        self.part_animals = {p: sorted(a for a, q in part_of.items() if q == p) for p in ("train", "validation", "test")}
        self.part_neurons = {p: [n for n in self.order if self.rows[n].part == p] for p in ("train", "validation", "test")}
        self.heldout_animals = frozenset(self.part_animals["validation"]) | frozenset(self.part_animals["test"])
        self.never_heldout = sorted(never_heldout)
        if never_heldout & self.heldout_animals:
            raise EngineError("an animal flagged as inspected, development or Li-shared is in validation or test")
        self.animal_neurons = defaultdict(list)
        for n in self.order:
            self.animal_neurons[self.rows[n].animal].append(n)
        for part, animals, neurons in world.expected:
            if len(self.part_animals[part]) != animals or len(self.part_neurons[part]) != neurons:
                raise EngineError(f"part {part}: {len(self.part_animals[part])} animals / {len(self.part_neurons[part])} "
                                  f"neurons, the locked plan has {animals} / {neurons}")
        self.cohort_sha = digest([[r.neuron, r.animal, r.full_id, r.source_label, r.tag, r.soma_label, r.hemisphere,
                                   r.recon_type, r.part, list(r.soma)] for r in self.rows.values()])

    def summary(self):
        return {"neurons": len(self.order), "animals": len(self.animal_part),
                "parts": {p: {"animals": len(self.part_animals[p]), "neurons": len(self.part_neurons[p])}
                          for p in ("train", "validation", "test")},
                "sha256": dict(self.sha), "cohort_sha256": self.cohort_sha}


# ----------------------------------------------------------------------------------------------
# sealed regional table: training rows only, held-out lines are skipped unparsed
# ----------------------------------------------------------------------------------------------


def load_training_profiles(world, cohort):
    """log10(1 + published length) for the 49 regions of TRAINING animals' neurons only.

    The file is read line by line. A line whose animal is not a training animal of the locked split
    is dropped before any of its values is parsed, so no validation or test regional value ever
    exists in memory, and nothing can be fitted on it. Returns (profiles by full ID, audit record).
    """
    pins = dict(world.pins)
    table_sha = sha_file(world.regional_table)
    if "regional_table" in pins and table_sha != pins["regional_table"]:
        raise EngineError("regional table differs from its pinned SHA-256")
    train_animals = frozenset(cohort.part_animals["train"])
    wanted = {cohort.rows[n].full_id for n in cohort.part_neurons["train"]}
    profiles, skipped_sealed = {}, 0
    with open(world.regional_table, newline="") as handle:
        header = handle.readline().rstrip("\r\n").split(",")
        if tuple(header[1:]) != REGIONS:
            raise EngineError("regional table does not have the 49 declared region columns")
        for line in handle:
            cut = line.find(",")
            full_id = line[:cut]
            if full_id.split("_", 1)[0] not in train_animals:
                skipped_sealed += 1          # sealed: never parsed, never stored
                continue
            values = np.array(line[cut + 1:].rstrip("\r\n").split(","), dtype=float)
            if values.shape != (len(REGIONS),) or not np.isfinite(values).all() or (values < 0).any():
                raise EngineError(f"invalid regional row {full_id}")
            if full_id in profiles or full_id not in wanted:
                raise EngineError(f"unexpected or duplicate training regional row {full_id}")
            profiles[full_id] = np.log10(1.0 + values)
    if set(profiles) != wanted:
        raise EngineError("regional table does not hold every training neuron exactly once")
    audit = {"regional_table_sha256": table_sha, "training_rows_loaded": len(profiles),
             "sealed_rows_skipped_unparsed": skipped_sealed,
             "training_animals_sha256": ids_sha(sorted(train_animals))}
    return profiles, audit


# ----------------------------------------------------------------------------------------------
# traces on disk and the stage guard
# ----------------------------------------------------------------------------------------------


def frame_is_true_hemisphere(frame):
    if frame.endswith("true_hemisphere"):
        return True
    if frame.endswith("mirrored_left"):
        return False
    raise EngineError(f"unknown trace frame '{frame}'")


class TraceStore:
    """Lists the traces on disk from their manifests (names and hashes only) and opens one only when
    the stage and the purpose allow it. Listing never opens a trace."""

    def __init__(self, world, cohort, stage):
        if stage not in STAGES:
            raise EngineError(f"unknown stage {stage}")
        self.world, self.cohort, self.stage = world, cohort, stage
        self.copies = defaultdict(list)
        self.manifest_sha = {}
        self._truth_unlocked = False
        for manifest, base, layout in world.trace_manifests:
            if not os.path.exists(manifest):
                continue
            self.manifest_sha[manifest] = sha_file(manifest)
            with open(manifest, newline="") as handle:
                for r in csv.DictReader(handle, delimiter="\t"):
                    unavailable = r.get("frame") == UNAVAILABLE_FRAME
                    if not unavailable and not r["file"].endswith(".swc"):
                        continue
                    neuron = r["neuron"]
                    row = cohort.rows.get(neuron)
                    if row is None:
                        raise EngineError(f"trace listing row '{r['file'] or neuron}' is not a neuron of the cohort")
                    if layout == "legacy_training" and row.part != "train":
                        raise SealError(f"{r['file']}: a held-out trace is listed in the training-only directory")
                    if layout == "by_part" and r["file"] and Path(r["file"]).parts[0] != row.part:
                        raise SealError(f"{r['file']}: stored under another part than the locked split says ({row.part})")
                    if unavailable:
                        # rule 8 (plan amendment 2): the source serves no trace; listed with its reason, never silently absent
                        self.copies[neuron].append({"path": None, "file": r["file"] or f"{row.part}/{row.animal}/{neuron}.swc",
                                                    "sha256": UNAVAILABLE_FRAME, "frame": r["frame"], "true_hemisphere": False,
                                                    "unavailable": True,
                                                    "reason": " | ".join(v for v in (r.get("how"), r.get("source")) if v)})
                        continue
                    self.copies[neuron].append({"path": str(Path(base) / r["file"]), "file": r["file"],
                                                "sha256": r["sha256"], "frame": r["frame"],
                                                "true_hemisphere": frame_is_true_hemisphere(r["frame"])})

    def available(self, part):
        """(neurons with a listed trace or listed as unavailable under rule 8, neurons of the part)."""
        neurons = self.cohort.part_neurons[part]
        return sum(1 for n in neurons if self.copies.get(n)), len(neurons)

    def choose(self, neuron):
        copies = self.copies.get(neuron)
        if not copies:
            return None
        real = [c for c in copies if not c.get("unavailable")]
        return sorted(real, key=lambda c: (not c["true_hemisphere"], c["file"]))[0] if real else copies[0]

    def unlock_truth(self, commitment_sha):
        """Called by `admit` only, after the committed forecast bundle has been verified."""
        if self.stage == "develop":
            raise SealError("the develop stage never opens held-out truth")
        if not commitment_sha:
            raise SealError("truth can be opened only against a prediction commitment")
        self._truth_unlocked = True

    def open(self, neuron, purpose):
        """Return the chosen copy of a trace, or refuse. purpose is 'donor' or 'truth'."""
        row = self.cohort.rows[neuron]
        if purpose not in ("donor", "truth"):
            raise EngineError("purpose must be 'donor' or 'truth'")
        if row.part == "train":
            if purpose == "truth" and self.stage != "develop":
                raise SealError("training truth is scored in the develop stage only")
        else:
            wanted = PART_OF_STAGE[self.stage]
            if purpose == "donor":
                raise SealError(f"{neuron}: a {row.part} animal can never be a donor")
            if row.part != wanted:
                raise SealError(f"{neuron}: stage '{self.stage}' may not open a {row.part} trace")
            if not self._truth_unlocked:
                raise SealError(f"{neuron}: held-out truth is opened only by `admit`, after the forecasts are committed")
        copy = self.choose(neuron)
        if copy is None:
            raise StageBlocked(f"no trace of {neuron} is on disk")
        return copy

    def require(self, stage, *, need_truth):
        """Refuse a sealed stage while its traces are not all on disk, and say what is missing."""
        have_train, total_train = self.available("train")
        if stage == "develop":
            if have_train == 0:
                raise StageBlocked("no training trace is on disk")
            if have_train < total_train and self.world.kind == "real" and not self.world.development_snapshot:
                raise StageBlocked(
                    f"stage 'develop' on the whole cohort needs the final training listing: {have_train} of {total_train} "
                    "training neurons are listed (a file row or an unavailable_at_source row each); the download may still be "
                    "running. The traces have not been fetched completely. The fixed development snapshot of the 37 traces "
                    "under data/raw/jiao2025/swc is `--world local37`")
            return
        part = PART_OF_STAGE[stage]
        have, total = self.available(part)
        problems = []
        if have_train < total_train:
            problems.append(f"{have_train} of {total_train} training traces are on disk (the dense donor library needs all)")
        if need_truth and have < total:
            problems.append(f"{have} of {total} {part} traces are on disk (scoring needs all)")
        if problems:
            where = "; ".join(m for m, _, _ in self.world.trace_manifests)
            raise StageBlocked(
                f"stage '{stage}' cannot run on {self.world.name}: " + "; ".join(problems) + ". "
                "Blocker: the traces have not been fetched completely. Plan amendment 2 (2026-11-11): the user chose the "
                "portal's per-neuron address; stage A (training + validation) is downloaded by the other agent's "
                "analysis/axon_trace_download.py; stage B (test) only after the freeze. This engine downloads nothing. "
                "Every neuron of a needed part must be listed in data/raw/jiao2025/traces/MANIFEST.tsv (columns file, neuron, "
                "bytes, sha256, source, frame) with a file <part>/<animal>/<neuron>.swc (frame ending in true_hemisphere or "
                f"mirrored_left) or as {UNAVAILABLE_FRAME}. Trace manifests looked for: {where}")


# ----------------------------------------------------------------------------------------------
# ingestion: plan section 9, frame rule, axon rule, canonical 50 um arbor
# ----------------------------------------------------------------------------------------------

RULE_TEXT = {
    "1": "tree: node 1 is the soma (type 1, parent -1), every parent exists, one connected tree, no cycle",
    "2": "soma within 5 um of the table soma after the frame rule",
    "3": "soma hemisphere equals the index",
    "4": "node types within {0,1,2,3,4,10,11,12}; the axon is coded as type 2 or type 0, not both",
    "5": "Axon_only files have no dendrite node; Axon_and_dendrite files have at least one",
    "6": "axon cable of at least 1,000 um",
    "7": "no clipping (amendment 1): the complete measured tree is used; cable outside the atlas box is logged, never removed",
    "8": ("unavailable at the source (plan amendment 2): the download listing records the neuron as unavailable_at_source; "
          "excluded, counted per part and per source label, never silently absent"),
}
UNAVAILABLE_FRAME = "unavailable_at_source"


def unavailable_log(neuron, row, copy):
    """The log entry of a neuron the download listing records as unavailable at the source (rule 8, plan amendment 2)."""
    return {"neuron": neuron, "animal": row.animal, "part": row.part, "file": copy["file"], "source_file": copy["file"],
            "manifest_sha256": copy["sha256"], "status": "excluded", "rules_violated": ["8"],
            "details": {"8": "listed as unavailable at the source" + (f": {copy['reason']}" if copy.get("reason") else "")}}


def rule8_report(entries, cohort, part):
    """Plan amendment 2, reporting duty: the neurons the source does not serve, per part and per source label, and
    whether a whole brain sample or more than 5% of a held-out part is lost (then the report says so in its first lines)."""
    lost = sorted(n for n, e in entries.items() if "8" in e.get("rules_violated", []))
    declared, lost_set = cohort.part_neurons[part], set(n for n, e in entries.items() if "8" in e.get("rules_violated", []))
    whole = sorted(a for a in cohort.part_animals[part] if all(n in lost_set for n in cohort.animal_neurons[a]))
    share = len(lost) / len(declared) if declared else 0.0
    return {"part": part, "neurons": len(lost), "of": len(declared), "share": share,
            "by_source_label": dict(sorted(Counter(cohort.rows[n].source_label for n in lost).items())),
            "brain_samples_lost_entirely": whole, "neuron_ids": lost,
            "scored_cohort_smaller_than_locked_notice": bool(part != "train" and (whole or share > 0.05))}


def read_swc(path):
    try:
        rows = np.loadtxt(str(path), comments="#", ndmin=2)
    except Exception as error:            # unreadable numeric table
        raise Exclusion("1", f"not a numeric SWC table ({type(error).__name__})")
    if rows.ndim != 2 or rows.shape[0] < 1 or rows.shape[1] != 7 or not np.isfinite(rows).all():
        raise Exclusion("1", "not a finite 7-column SWC table")
    return rows


def tree_structure(rows):
    """Rule 1. Returns ids, types, parent index, and the breadth-first order with children by ascending ID
    (the frozen extractor's own traversal, so that node order is the same)."""
    n = len(rows)
    ids, types, parent_ids = rows[:, 0].astype(np.int64), rows[:, 1].astype(np.int64), rows[:, 6].astype(np.int64)
    if not (np.array_equal(ids, rows[:, 0]) and np.array_equal(types, rows[:, 1]) and np.array_equal(parent_ids, rows[:, 6])):
        raise Exclusion("1", "node ID, type or parent is not an integer")
    id_list = ids.tolist()
    position = {v: i for i, v in enumerate(id_list)}
    if len(position) != n:
        raise Exclusion("1", "duplicate node IDs")
    roots = np.flatnonzero(parent_ids == -1)
    if len(roots) != 1:
        raise Exclusion("1", f"{len(roots)} root nodes")
    root = int(roots[0])
    if id_list[root] != 1 or types[root] != 1:
        raise Exclusion("1", "node 1 is not the type-1 root")
    parent = np.empty(n, dtype=np.int64)
    children = [[] for _ in range(n)]
    for i, p in enumerate(parent_ids.tolist()):
        if i == root:
            parent[i] = -1
            continue
        j = position.get(p)
        if j is None:
            raise Exclusion("1", f"node {id_list[i]} names a parent that does not exist")
        parent[i] = j
        children[j].append(i)
    order, queue = [], deque([root])
    while queue:
        i = queue.popleft()
        order.append(i)
        kids = children[i]
        if len(kids) > 1:
            kids.sort(key=id_list.__getitem__)
        queue.extend(kids)
    if len(order) != n:
        raise Exclusion("1", "not one connected tree without a cycle")
    return {"ids": ids, "types": types, "parent": parent, "order": order, "root": root}


def outside_box_cable(xyz, parent):
    """Length of cable outside the atlas box, exactly, edge by edge (segment against box). It is only
    logged: plan amendment 1 forbids removing, moving or reconnecting any measured cable."""
    inside_node = ((xyz >= 0.0) & (xyz <= BOX)).all(axis=1)
    info = {"nodes_outside": int((~inside_node).sum()), "outside_um": 0.0, "soma_outside_box": bool(not inside_node[0])}
    if len(xyz) < 2 or inside_node.all():
        return info
    a, b = xyz[parent[1:]], xyz[1:]
    d = b - a
    length = np.linalg.norm(d, axis=1)
    t0, t1 = np.zeros(len(d)), np.ones(len(d))
    with np.errstate(divide="ignore", invalid="ignore"):
        for k in range(3):
            still = d[:, k] == 0.0
            lo, hi = (0.0 - a[:, k]) / d[:, k], (BOX[k] - a[:, k]) / d[:, k]
            enter, leave = np.minimum(lo, hi), np.maximum(lo, hi)
            miss = still & ((a[:, k] < 0.0) | (a[:, k] > BOX[k]))
            enter = np.where(still, np.where(miss, np.inf, -np.inf), enter)
            leave = np.where(still, np.where(miss, -np.inf, np.inf), leave)
            t0, t1 = np.maximum(t0, enter), np.minimum(t1, leave)
    inside_fraction = np.clip(t1 - t0, 0.0, 1.0)
    info["outside_um"] = float((length * (1.0 - inside_fraction)).sum())
    return info


def extract_axon(rows, structure, *, reflect, mode="plan"):
    """The axon as a soma-rooted tree in the mirrored-left frame.

    mode 'plan'   the benchmark's rule (plan section 2): every node that is not dendrite (3, 4) and not the
                  soma is axon; a marker (10, 11, 12) takes the type of its nearest non-marker ancestor.
    mode 'frozen' exactly the node set of the frozen extractor d110_geometry.extract_explicit_axon: nodes of
                  raw type 2 and their ancestors (used only to prove that the two paths agree).
    Ancestors that are not axon are kept as connector geometry, as the frozen extractor does for markers.
    """
    types, parent, order, root = structure["types"], structure["parent"], structure["order"], structure["root"]
    n = len(types)
    effective = types.copy()
    for i in order[1:]:
        if types[i] in MARKER_TYPES:
            effective[i] = effective[parent[i]]
    present = set(np.unique(types).tolist())
    if not present <= set(ALLOWED_TYPES):
        raise Exclusion("4", f"node types outside the allowed set: {sorted(present - set(ALLOWED_TYPES))}")
    has2, has0 = 2 in present, 0 in present
    if has2 and has0:
        raise Exclusion("4", "the axon is coded both as type 2 and as type 0")
    code = 2 if has2 else (0 if has0 else None)
    explicit = (types == code) if code is not None else np.zeros(n, dtype=bool)
    if mode == "frozen":
        if code != 2:
            raise EngineError("the frozen extractor accepts explicit type-2 axons only")
        axon = explicit.copy()
    else:
        axon = np.isin(effective, (0, 2))
    axon[root] = False

    def closure(flags):
        keep = flags.copy()
        keep[root] = True
        for i in reversed(order):
            if keep[i] and i != root:
                keep[parent[i]] = True
        return keep

    keep = closure(axon)
    explicit_keep = closure(explicit) if mode == "plan" else keep     # the frozen extractor's node set
    kept_order = [i for i in order if keep[i]]
    if mode == "frozen" and not all(types[i] in (1, 2, 10, 11, 12) for i in kept_order):
        raise EngineError("ambiguous ancestor in explicit axon (the frozen extractor refuses this file)")
    remap = np.full(n, -1, dtype=np.int64)
    remap[kept_order] = np.arange(len(kept_order))
    xyz = rows[kept_order, 2:5].copy()
    if reflect:
        xyz[:, 2] = MIRROR_ML - xyz[:, 2]
    new_parent = np.array([-1] + [int(remap[parent[i]]) for i in kept_order[1:]], dtype=int)
    seg = np.zeros(len(kept_order))
    if len(kept_order) > 1:
        seg[1:] = np.linalg.norm(xyz[1:] - xyz[new_parent[1:]], axis=1)
    kept = np.asarray(kept_order)
    is_axon = axon[kept]
    connector = ~is_axon
    connector[0] = False
    twig = is_axon & ~explicit_keep[kept]
    stats = {"axon_coding": code, "axon_nodes": int(is_axon.sum()), "axon_cable_um": float(seg[is_axon].sum()),
             "terminal_marker_twig_nodes": int(twig.sum()), "terminal_marker_twig_um": float(seg[twig].sum()),
             "connector_nodes": int(connector.sum()), "connector_um": float(seg[connector].sum()),
             "dendrite_connector_nodes": int(np.isin(effective[kept][connector], DENDRITE_TYPES).sum())}
    return xyz, new_parent, types[kept_order], stats


def outside_summary(entries):
    """Per part: how much measured cable lies outside the atlas box (logged, never removed)."""
    values = [e["outside_atlas_box"]["outside_um"] for e in entries if e.get("status", "included") == "included"]
    return {"neurons": len(values), "neurons_with_cable_outside": int(sum(v > 0 for v in values)),
            "total_um": float(sum(values)), "max_um": float(max(values, default=0.0))}


def connector_summary(entries):
    """Per part: non-axon connector geometry kept in the scored trees (markers inside a path; dendritic stems)."""
    kept = [e for e in entries if e.get("status", "included") == "included"]
    stems = [e for e in kept if e.get("dendrite_connector_nodes", 0) > 0]
    return {"neurons": len(kept), "neurons_with_connector_nodes": int(sum(e.get("connector_nodes", 0) > 0 for e in kept)),
            "connector_um_total": float(sum(e.get("connector_um", 0.0) for e in kept)),
            "neurons_whose_axon_leaves_a_dendrite": len(stems),
            "dendritic_stem_nodes_total": int(sum(e["dendrite_connector_nodes"] for e in stems))}


def canonical_arbor(xyz, parent, types):
    """The frozen code's own canonical form: branch-to-branch resampling at 50 um."""
    return G.resample(G.Arbor(xyz, parent, types).validate(), STEP)


def frozen_replay_arbor(path, reflect):
    """This engine's extraction restricted to the frozen extractor's node set, no clipping."""
    rows = read_swc(path)
    xyz, parent, types, _ = extract_axon(rows, tree_structure(rows), reflect=reflect, mode="frozen")
    return canonical_arbor(xyz, parent, types)


def ingest_trace(path, *, frame, hemisphere, recon_type, table_soma, expected_sha=None):
    """Read one trace, apply plan section 9, return the log and (when included) the canonical arbor.

    Every violated rule is listed. Nothing here depends on any model or on any forecast.
    """
    path = Path(path)
    log = {"file": path.name, "file_sha256": sha_file(path), "frame": frame, "status": "excluded", "rules_violated": [],
           "details": {}}
    if expected_sha is not None and log["file_sha256"] != expected_sha:
        raise EngineError(f"{path}: bytes differ from the trace manifest")
    try:
        rows = read_swc(path)
        structure = tree_structure(rows)
    except Exclusion as error:
        log["rules_violated"].append(error.rule)
        log["details"][error.rule] = error.detail
        return log, None
    types = structure["types"]
    log["nodes"] = int(len(rows))
    log["type_counts"] = {str(k): int(v) for k, v in sorted(Counter(types.tolist()).items())}
    true_frame = frame_is_true_hemisphere(frame)
    reflect = bool(true_frame and hemisphere == "Right")
    soma_raw = rows[structure["root"], 2:5].copy()
    file_hemisphere = "Right" if soma_raw[2] > MIDLINE_ML else "Left"
    soma = soma_raw.copy()
    if reflect:
        soma[2] = MIRROR_ML - soma[2]
    log.update({"mirrored_to_left": reflect, "file_hemisphere": file_hemisphere, "index_hemisphere": hemisphere,
                "trace_soma_um": soma.tolist(), "table_soma_um": [float(v) for v in table_soma]})
    violations = {}
    if hemisphere not in ("Left", "Right"):
        violations["3"] = f"index hemisphere '{hemisphere}' is neither Left nor Right"
    elif true_frame and file_hemisphere != hemisphere:
        violations["3"] = f"file soma is on the {file_hemisphere.lower()}, the index says {hemisphere}"
    elif not true_frame and file_hemisphere != "Left":
        violations["3"] = "a mirrored-left file has its soma on the right"
    delta = float(np.linalg.norm(soma - np.asarray(table_soma, dtype=float)))
    log["soma_table_distance_um"] = delta
    if delta > CONFIG["inclusion"]["soma_tolerance_um"]:
        violations["2"] = f"trace soma is {delta:.2f} um from the table soma"
    dendrite_nodes = int(np.isin(types, DENDRITE_TYPES).sum())
    log["dendrite_nodes"] = dendrite_nodes
    if recon_type == "Axon_only":
        if dendrite_nodes:
            violations["5"] = "Axon_only file has dendrite nodes"
    elif recon_type == "Axon_and_dendrite":
        if not dendrite_nodes:
            violations["5"] = "Axon_and_dendrite file has no dendrite node (its dendrite could be counted as axon)"
    else:
        violations["5"] = f"unknown reconstruction type '{recon_type}'"
    try:
        xyz, parent, kept_types, stats = extract_axon(rows, structure, reflect=reflect, mode="plan")
    except Exclusion as error:
        violations[error.rule] = error.detail
        log["rules_violated"] = sorted(violations)
        log["details"] = violations
        return log, None
    log.update(stats)
    if stats["axon_cable_um"] < CONFIG["inclusion"]["min_axon_cable_um"]:
        violations["6"] = f"axon cable {stats['axon_cable_um']:.1f} um"
    if violations:
        log["rules_violated"] = sorted(violations)
        log["details"] = violations
        return log, None
    arbor = canonical_arbor(xyz, parent, kept_types)                 # the complete measured tree: nothing is clipped
    log.update({"status": "included", "outside_atlas_box": outside_box_cable(xyz, parent),
                "canonical_nodes": int(len(arbor.xyz)),
                "canonical_cable_um": arbor.cable(), "canonical_branches": int(arbor.branches()),
                "arbor_sha256": arbor_sha(arbor)})
    return log, arbor


# ----------------------------------------------------------------------------------------------
# all-pairs cache of the frozen metric on soma-centred training arbors, and the memo that lets the
# unmodified frozen prototype read it
# ----------------------------------------------------------------------------------------------

_WORKER = {}


def pair_metric_identity():
    return {"definition": ("frozen d110_geometry.geometry_metrics(a.shifted(0), b.shifted(0)) with the lower neuron ID first: "
                           "distance_um and 1 - f1_200, exactly as frozen brm_template.pair_costs fills its matrices"),
            "analysis/d110_geometry.py": PINNED_SOURCES["analysis/d110_geometry.py"],
            "analysis/brm_template.py": PINNED_SOURCES["analysis/brm_template.py"]}


def pair_cost(a_centred, b_centred):
    """One cache entry: the frozen metric itself, nothing else."""
    metrics = G.geometry_metrics(a_centred, b_centred)
    return metrics["distance_um"], 1.0 - metrics["f1_200"]


def _pair_block_init(library_dir, ids):
    _WORKER["pair_library"], _WORKER["pair_ids"], _WORKER["pair_centred"] = library_dir, ids, {}


def _pair_centred(i):
    centred = _WORKER["pair_centred"]
    arbor = centred.get(i)
    if arbor is None:
        arbor = G.load_arbor(Path(_WORKER["pair_library"]) / f"{_WORKER['pair_ids'][i]}.npz").shifted(np.zeros(3))
        centred[i] = arbor
    return arbor


def _pair_block_task(block):
    out = []
    for i, js in block:
        a = _pair_centred(i)
        d, f = np.empty(len(js)), np.empty(len(js))
        for k, j in enumerate(np.asarray(js).tolist()):
            d[k], f[k] = pair_cost(a, _pair_centred(j))
        out.append((int(i), np.asarray(js, dtype=np.int64), d, f))
    return out


class PairCache:
    """TRAIN-ONLY all-pairs pair costs of one training library, on disk, pinned and resumable.

    D[i, j] = distance_um and F[i, j] = 1 - f1_200 of the frozen metric between the soma-centred canonical
    arbors of library neurons i and j (sorted neuron IDs; the lower ID is the first argument, as in the
    frozen pair_costs). A cache belongs to one library fingerprint (ordered IDs and arbor content hashes) and
    one metric code hash. K marks the pairs already computed, so a build can stop and resume; entries of an
    older cache whose two arbors are unchanged are imported instead of recomputed."""

    SCHEMA = "axon_engine_pair_cache_v1"

    def __init__(self, root, library_dir, ids, arbor_hashes, *, siblings=()):
        self.root, self.library_dir = Path(root), Path(library_dir)
        self.ids = list(ids)
        if self.ids != sorted(self.ids) or len(set(self.ids)) != len(self.ids):
            raise EngineError("pair cache: library IDs must be distinct and sorted")
        self.hashes = [arbor_hashes[n] for n in self.ids]
        self.n = len(self.ids)
        self.pairs = self.n * (self.n - 1) // 2
        self.fingerprint = digest([[n, h] for n, h in zip(self.ids, self.hashes)])
        self.dir = self.root / self.fingerprint[:16]
        self.index = {n: i for i, n in enumerate(self.ids)}
        self.siblings = [Path(s) for s in siblings]
        self._D = self._F = None
        self.record = None

    # -- files -----------------------------------------------------------------------------
    def index_record(self):
        return {"schema": self.SCHEMA, "n": self.n, "pairs": self.pairs, "ids": self.ids, "arbor_sha256": self.hashes,
                "ids_sha256": ids_sha(self.ids), "library_fingerprint": self.fingerprint, "metric": pair_metric_identity(),
                "layout": "D.npy, F.npy: float64 n x n, symmetric, zero diagonal; K.npy: uint8 n x n, 1 where the pair is computed"}

    def _memmap(self, name, mode):
        return np.load(self.dir / name, mmap_mode=mode)

    def exists(self):
        return (self.dir / "index.json").exists()

    def is_complete(self):
        return (self.dir / "COMPLETE.json").exists()

    def _create(self, log=None):
        self.dir.mkdir(parents=True, exist_ok=True)
        if self.exists():
            if read_json(self.dir / "index.json") != jsonable(self.index_record()):
                raise EngineError(f"pair cache {self.dir}: its index does not describe this library and metric")
            return
        for name, dtype in (("D.npy", np.float64), ("F.npy", np.float64), ("K.npy", np.uint8)):
            array = np.lib.format.open_memmap(self.dir / name, mode="w+", dtype=dtype, shape=(self.n, self.n))
            array.flush()
            del array
        imported = self._import_known(log)
        write_json(self.dir / "index.json", self.index_record(), pretty=True)
        write_json(self.dir / "import.json", {"imported_pairs": imported}, pretty=True)

    def _import_known(self, log=None):
        """Entries of other caches (same metric code) whose two arbors are bit for bit the same."""
        total = 0
        candidates = []
        for root in [self.root] + self.siblings:
            if root.exists():
                candidates += [p for p in sorted(root.iterdir()) if p.is_dir() and p != self.dir and (p / "index.json").exists()]
        if not candidates:
            return 0
        D, F, K = self._memmap("D.npy", "r+"), self._memmap("F.npy", "r+"), self._memmap("K.npy", "r+")
        for other in candidates:
            index = read_json(other / "index.json")
            if index.get("schema") != self.SCHEMA or index.get("metric") != pair_metric_identity():
                continue
            theirs = {n: (i, h) for i, (n, h) in enumerate(zip(index["ids"], index["arbor_sha256"]))}
            shared = [n for n, h in zip(self.ids, self.hashes) if n in theirs and theirs[n][1] == h]
            if len(shared) < 2:
                continue
            mine = np.array([self.index[n] for n in shared])
            old = np.array([theirs[n][0] for n in shared])
            oK = np.load(other / "K.npy", mmap_mode="r")[np.ix_(old, old)].astype(bool)
            new = oK & ~K[np.ix_(mine, mine)].astype(bool)
            if not new.any():
                continue
            for name, target in (("D.npy", D), ("F.npy", F)):
                block = target[np.ix_(mine, mine)]
                block[new] = np.load(other / name, mmap_mode="r")[np.ix_(old, old)][new]
                target[np.ix_(mine, mine)] = block
            D.flush()
            F.flush()
            block = K[np.ix_(mine, mine)]
            block[new] = 1
            K[np.ix_(mine, mine)] = block
            K.flush()
            total += int(np.triu(new, 1).sum())
        if total and log:
            log(f"pair cache: {total} pairs imported from caches of unchanged arbors")
        return total

    # -- build (resumable) -----------------------------------------------------------------
    def progress(self):
        if not self.exists():
            return {"exists": False, "complete": False, "pairs": self.pairs, "pairs_known": 0}
        K = self._memmap("K.npy", "r")
        known = int(np.triu(np.asarray(K), 1).sum()) if self.n else 0
        return {"exists": True, "complete": self.is_complete(), "pairs": self.pairs, "pairs_known": known,
                "directory": str(self.dir)}

    def build(self, *, workers=1, pool_factory=None, minutes=None, log=None):
        """Compute the missing pairs. Stops cleanly when the time budget is used; call again to resume."""
        started = time.time()
        if self.is_complete():
            return {**self.progress(), "computed_now": 0, "seconds": 0.0}
        self._create(log)
        D, F, K = self._memmap("D.npy", "r+"), self._memmap("F.npy", "r+"), self._memmap("K.npy", "r+")
        size = CONFIG["pair_cache"]["block_pairs"]
        todo = self.pairs - (int(np.triu(np.asarray(K), 1).sum()) if self.n else 0)

        def blocks():
            # Row i contributes its not yet computed pairs (i, j > i). A result only ever marks entries of rows
            # this generator has already passed (upper triangle) or of the lower triangle, so nothing is skipped.
            block, count = [], 0
            for i in range(self.n - 1):
                js = (np.flatnonzero(np.asarray(K[i, i + 1:]) == 0) + i + 1).astype(np.int32)
                while len(js):
                    take, js = js[:size - count], js[size - count:]
                    block.append((i, take))
                    count += len(take)
                    if count >= size:
                        yield block
                        block, count = [], 0
            if block:
                yield block

        if log and todo:
            log(f"pair cache: {todo} of {self.pairs} soma-centred pairs to compute")
        deadline = None if minutes is None else started + 60.0 * float(minutes)
        computed, pending = 0, []

        def absorb(result):
            nonlocal computed
            for i, js, d, f in result:
                D[i, js] = d
                D[js, i] = d
                F[i, js] = f
                F[js, i] = f
                pending.append((i, js))
                computed += len(js)

        def commit():
            if not pending:
                return
            D.flush()
            F.flush()                      # values first, the 'computed' marks only afterwards
            for i, js in pending:
                K[i, js] = 1
                K[js, i] = 1
            K.flush()
            pending.clear()

        if todo:
            if workers > 1 and todo > 2 * size and pool_factory is not None:
                with pool_factory(workers, _pair_block_init, (str(self.library_dir), self.ids)) as pool:
                    for done, result in enumerate(pool.imap_unordered(_pair_block_task, blocks()), 1):
                        absorb(result)
                        if done % 8 == 0:
                            commit()
                        if deadline is not None and time.time() > deadline and computed < todo:
                            pool.terminate()
                            break
            else:
                _pair_block_init(str(self.library_dir), self.ids)
                for done, b in enumerate(blocks(), 1):
                    absorb(_pair_block_task(b))
                    if done % 8 == 0:
                        commit()
                    if deadline is not None and time.time() > deadline and computed < todo:
                        break
            commit()
        seconds = time.time() - started
        status = {"computed_now": computed, "seconds": round(seconds, 2), "workers": workers,
                  "ms_per_pair_wall": (1000.0 * seconds / computed) if computed else None}
        if computed == todo:
            self._finalise(status)
        return {**self.progress(), **status}

    def _finalise(self, status):
        K = np.asarray(self._memmap("K.npy", "r"))
        if self.n > 1 and int(np.triu(K, 1).sum()) != self.pairs:
            raise EngineError("pair cache: not every pair is computed")
        D, F = np.asarray(self._memmap("D.npy", "r")), np.asarray(self._memmap("F.npy", "r"))
        if not (np.isfinite(D).all() and np.isfinite(F).all() and (D >= 0).all() and np.array_equal(D, D.T)
                and np.array_equal(F, F.T) and not D.diagonal().any() and not F.diagonal().any()):
            raise EngineError("pair cache: matrices are not finite, symmetric and zero on the diagonal")
        record = {"schema": self.SCHEMA, "library_fingerprint": self.fingerprint, "ids_sha256": ids_sha(self.ids),
                  "n": self.n, "pairs": self.pairs, "metric": pair_metric_identity(),
                  "index_sha256": sha_file(self.dir / "index.json"), "D_sha256": sha_file(self.dir / "D.npy"),
                  "F_sha256": sha_file(self.dir / "F.npy")}
        record["cache_sha256"] = digest(record)
        record["bytes"] = {name: (self.dir / name).stat().st_size for name in ("D.npy", "F.npy", "K.npy")}
        record["last_build"] = status
        record["completed_utc"] = utc_now()
        write_json(self.dir / "COMPLETE.json", record, pretty=True)

    # -- use -------------------------------------------------------------------------------
    def load(self, *, verify_files=True):
        if not self.is_complete():
            p = self.progress()
            raise EngineError(f"the all-pairs cache of this training library is not complete ({p['pairs_known']} of {p['pairs']} "
                              "pairs); run `cache` (it resumes)")
        if read_json(self.dir / "index.json") != jsonable(self.index_record()):
            raise EngineError("pair cache: its index does not describe this library and metric")
        record = read_json(self.dir / "COMPLETE.json")
        pinned = {k: record[k] for k in ("schema", "library_fingerprint", "ids_sha256", "n", "pairs", "metric", "index_sha256",
                                         "D_sha256", "F_sha256")}
        if (record["cache_sha256"] != digest(pinned) or record["library_fingerprint"] != self.fingerprint
                or record["metric"] != pair_metric_identity() or record["index_sha256"] != sha_file(self.dir / "index.json")):
            raise EngineError("pair cache: its completion record does not describe this library, metric and index")
        if verify_files and (sha_file(self.dir / "D.npy") != record["D_sha256"] or sha_file(self.dir / "F.npy") != record["F_sha256"]):
            raise EngineError("pair cache: a matrix file changed after the cache was completed")
        self._D, self._F = self._memmap("D.npy", "r"), self._memmap("F.npy", "r")
        self.record = record
        return self

    @property
    def sha256(self):
        return self.record["cache_sha256"]

    def sub(self, ids):
        """(D, F) for the given library IDs, in that order, read-only."""
        index = [self.index[n] for n in ids]
        if index == list(range(self.n)):
            D, F = np.asarray(self._D), np.asarray(self._F)
        else:
            grid = np.ix_(index, index)
            D, F = np.asarray(self._D[grid]), np.asarray(self._F[grid])
            D.flags.writeable = False
            F.flags.writeable = False
        return D, F

    def sample_pairs(self, count):
        """A fixed sample of pairs, chosen by hash of the library fingerprint."""
        if self.pairs <= count:
            return [(i, j) for i in range(self.n) for j in range(i + 1, self.n)]
        chosen, k = [], 0
        while len(chosen) < count:
            h = hashlib.sha256(f"{self.fingerprint}:{k}".encode()).digest()
            i, j = int.from_bytes(h[:8], "big") % self.n, int.from_bytes(h[8:16], "big") % self.n
            k += 1
            if i != j and (min(i, j), max(i, j)) not in chosen:
                chosen.append((min(i, j), max(i, j)))
        return chosen

    def prove(self, count=None, *, with_referee=True):
        """Recompute sampled pairs with the frozen pair_costs and with the referee's own loader and metric
        (the code its copying_baselines uses) and compare with the stored numbers bit for bit."""
        count = CONFIG["pair_cache"]["spot_check_pairs"] if count is None else count
        pairs = self.sample_pairs(count)
        Ref = referee() if with_referee else None
        frozen_ok = referee_ok = True
        seconds = {"frozen": 0.0, "referee": 0.0}
        for i, j in pairs:
            paths = [self.library_dir / f"{self.ids[k]}.npz" for k in (i, j)]
            donors = [{"arbor": G.load_arbor(p)} for p in paths]
            t = time.time()
            D, F = T.pair_costs(donors)                                   # the frozen function itself
            seconds["frozen"] += time.time() - t
            frozen_ok &= bool(D[0, 1] == self._D[i, j] == self._D[j, i] and F[0, 1] == self._F[i, j] == self._F[j, i])
            if Ref is not None:
                trees = [Ref.G.shifted(Ref.G.load_tree(p), np.zeros(3)) for p in paths]
                t = time.time()
                m = Ref.G.geometry_metrics(trees[0], trees[1])
                seconds["referee"] += time.time() - t
                referee_ok &= bool(m["distance_um"] == self._D[i, j] and 1.0 - m["f1_200"] == self._F[i, j])
        out = {"pairs_checked": len(pairs), "of": self.pairs, "frozen_pair_costs_bit_identical": bool(frozen_ok),
               "ms_per_pair_frozen": 1000.0 * seconds["frozen"] / max(1, len(pairs))}
        if Ref is not None:
            out.update({"referee_metric_bit_identical": bool(referee_ok),
                        "ms_per_pair_referee": 1000.0 * seconds["referee"] / max(1, len(pairs))})
        return out


class FrozenMemo:
    """Hands the stored pair costs to the UNMODIFIED frozen functions.

    Inside `frozen_memo()` the frozen brm_template.pair_costs and d110_geometry.training_medoid are
    replaced by look-ups into the pinned cache for donors whose arbors are the registered library arbors;
    any other donor list is passed to the original functions. The medoid uses the frozen function's own
    arithmetic (animal weights, left-to-right sum) on the stored distances. `prove_memo` shows that trees and
    records of the frozen predict_policy are bit-identical with and without the memo."""

    def __init__(self, cache_getter):
        self._cache_getter = cache_getter
        self._arbor, self._pool, self._medoid = {}, (None, None, None), {}

    @property
    def cache(self):
        return self._cache_getter()

    def register(self, neuron, arbor):
        known = self._arbor.get(neuron)
        if known is not None and known is not arbor:
            raise EngineError(f"another arbor is already registered for {neuron}")
        for array in (arbor.xyz, arbor.parent, arbor.types):
            array.flags.writeable = False
        self._arbor[neuron] = arbor

    def _registered(self, records):
        return all(self._arbor.get(r["neuron"]) is r["arbor"] for r in records)

    def matrices(self, ids):
        """(D, F) of one donor pool; the most recent pool is kept."""
        key = ids_sha(ids)
        if self._pool[0] != key:
            D, F = self.cache.sub(ids)
            self._pool = (key, D, F)
        return self._pool[1], self._pool[2]

    # -- drop-in replacements with the frozen functions' signatures --------------------------
    def _pair_costs(self, donors):
        if not donors:
            raise ValueError("No geometry donors.")
        if not self._registered(donors):
            return self._original_pair_costs(donors)
        return self.matrices([d["neuron"] for d in donors])

    def _training_medoid(self, training):
        if not training or not self._registered(training):
            return self._original_training_medoid(training)
        ids = [t["neuron"] for t in training]
        key = ids_sha(ids)
        index = self._medoid.get(key)
        if index is None:
            D, _ = self.matrices(ids)
            count = Counter(t["animal"] for t in training)
            weights = np.array([1 / count[t["animal"]] for t in training], dtype=float)
            weights /= weights.sum()
            losses = []
            for i in range(len(ids)):                # same left-to-right sum as the frozen function
                losses.append(sum(w * d for w, d in zip(weights, D[i].tolist())))
            index = int(np.argmin(losses))
            self._medoid[key] = index
        return index

    @contextmanager
    def frozen_memo(self):
        self._original_pair_costs, self._original_training_medoid = T.pair_costs, G.training_medoid
        T.pair_costs, G.training_medoid = self._pair_costs, self._training_medoid
        try:
            yield self
        finally:
            T.pair_costs, G.training_medoid = self._original_pair_costs, self._original_training_medoid


# ----------------------------------------------------------------------------------------------
# regional states (frozen joint model, training animals only)
# ----------------------------------------------------------------------------------------------


def oof_fold_of(animals, label_of, config=CONFIG):
    """Deterministic folds of TRAINING animals, dealt in hash order inside each source label."""
    spec = config["out_of_fold"]
    folds, offset = {}, 0
    for label in sorted({label_of[a] for a in animals}):
        group = sorted((a for a in animals if label_of[a] == label),
                       key=lambda a: sha_bytes((spec["salt"] + ":" + a).encode()))
        for a in group:
            folds[a] = offset % spec["folds"]
            offset += 1
    return folds


def fit_regional_state(state_id, cohort, profiles, train_animals, served_neurons):
    """Fit the frozen joint regional model on the given TRAINING animals and wrap it as the frozen
    RegionalState with the membership proof the frozen code demands.

    Refuses when any fitting animal is not a training animal of the locked split. The fit itself is
    the unmodified brm_real.fit_outer_fold (K chosen by its own inner whole-animal folds inside the
    fitting animals). Its `query` argument only fills an interface slot: a fixed placeholder is
    passed, so the state cannot depend on which neurons it later serves.
    """
    train_animals = sorted(set(train_animals))
    leak = [a for a in train_animals if a in cohort.heldout_animals or cohort.animal_part.get(a) != "train"]
    if leak:
        raise LeakError(f"regional state '{state_id}' would be fitted on animals outside the training part: {leak[:5]}")
    served = [cohort.rows[n] for n in served_neurons]
    overlap = {r.animal for r in served} & set(train_animals)
    if overlap:
        raise LeakError(f"regional state '{state_id}' would serve its own fitting animals: {sorted(overlap)[:5]}")
    inside = set(train_animals)
    rows = [cohort.rows[n] for n in cohort.order if cohort.rows[n].animal in inside]
    ids = np.array([r.full_id for r in rows])
    train = R.Inputs(ids, np.array([r.animal for r in rows]), np.array([r.soma for r in rows], dtype=float),
                     np.array([r.tag for r in rows]))
    Y = np.array([profiles[r.full_id] for r in rows])
    slot = R.Inputs(np.array(["QUERYSLOT_000_slot_NA"]), np.array(["QUERYSLOT"]),
                    train.X.mean(axis=0, keepdims=True), np.array(["slot"]))
    inner = R.inner_folds(train, np.arange(len(ids)), state_id)
    _, record = R.fit_outer_fold(train, Y, slot, inner, fold_name=state_id)
    model = record["outer_model_state"]
    proof = {"name": state_id, "model_state_sha256": T._digest(model), "train_ids": ids.tolist(),
             "train_animals": train_animals, "test_ids": [r.full_id for r in served],
             "test_animals": sorted({r.animal for r in served})}
    fit = {"selected_k": record["selected_k"],
           "selection": [{"k": s["k"], "mean_animal_mse": s["mean_animal_mse"]} for s in record["selection"]],
           "converged": [bool(b["converged"]) for b in record["outer_fit_diagnostics"]["blocks"]],
           "iterations": [int(b["iterations"]) for b in record["outer_fit_diagnostics"]["blocks"]]}
    return Geo.RegionalState(model, proof), fit


def _state_task(args):
    world_json, state_id, train_animals, served = args
    engine = Engine(World.from_json(world_json), "develop", workers=1, log=None)
    return engine.states.build(state_id, train_animals, served)


class StateBook:
    """Fits, stores and serves the regional states; keeps the state manifest."""

    def __init__(self, engine):
        self.engine = engine
        self.dir = Path(engine.world.out) / "states"
        self._states, self._sets = {}, {}
        cohort = engine.cohort
        self.fold = oof_fold_of(cohort.part_animals["train"], cohort.animal_label)
        self.plan = {"all_training": {"role": "all training animals; serves validation and test queries",
                                      "train_animals": list(cohort.part_animals["train"]),
                                      "served": cohort.part_neurons["validation"] + cohort.part_neurons["test"]}}
        for k in range(CONFIG["out_of_fold"]["folds"]):
            held = sorted(a for a, f in self.fold.items() if f == k)
            if not held:
                continue
            self.plan[f"oof_{k}"] = {
                "role": f"out-of-fold state {k}: training animals outside fold {k}; serves the training animals of fold {k}",
                "train_animals": sorted(a for a, f in self.fold.items() if f != k),
                "served": [n for a in held for n in cohort.animal_neurons[a]]}

    def state_id_for(self, neuron):
        row = self.engine.cohort.rows[neuron]
        return "all_training" if row.part != "train" else f"oof_{self.fold[row.animal]}"

    def _key(self, state_id):
        plan = self.plan[state_id]
        return digest({"state": state_id, "train_animals_sha256": ids_sha(plan["train_animals"]),
                       "served_sha256": ids_sha(plan["served"]), "cohort": self.engine.cohort.cohort_sha,
                       "regional_table": sha_file(self.engine.world.regional_table),
                       "frozen": FROZEN_PINS["source_sha256"], "wrapper": engine_code_sha().get("analysis/axon_engine.py")})

    def build(self, state_id, train_animals=None, served=None):
        """Fit (or reuse) one state; returns its manifest entry."""
        plan = self.plan[state_id]
        train_animals = plan["train_animals"] if train_animals is None else train_animals
        served = plan["served"] if served is None else served
        meta_path = self.dir / f"{state_id}.meta.json"
        if meta_path.exists():
            meta = read_json(meta_path)
            if meta.get("cache_key") == self._key(state_id) and (self.dir / f"{state_id}.json").exists():
                return meta["entry"]
        started = time.time()
        state, fit = fit_regional_state(state_id, self.engine.cohort, self.engine.profiles, train_animals, served)
        self.dir.mkdir(parents=True, exist_ok=True)
        write_json(self.dir / f"{state_id}.json", state.model)
        write_json(self.dir / f"{state_id}.proof.json", state.proof)
        if sha_file(self.dir / f"{state_id}.json") != state.proof["model_state_sha256"]:
            raise EngineError("saved regional state does not hash to its recorded identity")
        heldout = self.engine.cohort.heldout_animals
        entry = {"state_id": state_id, "role": plan["role"], "model_state_sha256": state.proof["model_state_sha256"],
                 "animals_sha256": ids_sha(sorted(state.proof["train_animals"])),
                 "animals": len(state.proof["train_animals"]), "neurons": len(state.proof["train_ids"]),
                 "training_ids_sha256": state.model["training_ids_sha256"],
                 "served_animals_sha256": ids_sha(state.proof["test_animals"]),
                 "served_neurons": len(state.proof["test_ids"]),
                 "validation_or_test_animals_in_state": sorted(set(state.proof["train_animals"]) & heldout),
                 "fit": fit, "seconds": round(time.time() - started, 1),
                 "code_sha256": {"wrapper": engine_code_sha().get("analysis/axon_engine.py"),
                                 "analysis/brm_real.py": PINNED_SOURCES["analysis/brm_real.py"],
                                 "analysis/brm_joint.py": PINNED_SOURCES["analysis/brm_joint.py"],
                                 "analysis/brm_geometry.py": PINNED_SOURCES["analysis/brm_geometry.py"]}}
        if entry["validation_or_test_animals_in_state"]:
            raise LeakError(f"state {state_id} contains held-out animals")
        write_json(meta_path, {"cache_key": self._key(state_id), "entry": entry}, pretty=True)
        return entry

    def build_many(self, state_ids, workers=1):
        entries = {}
        todo = list(state_ids)
        if workers > 1 and len(todo) > 1 and self.engine.can_parallel():
            with self.engine.pool(min(workers, len(todo))) as pool:
                jobs = [(self.engine.world.to_json(), s, None, None) for s in todo]
                for s, entry in zip(todo, pool.map(_state_task, jobs)):
                    entries[s] = entry
        else:
            for s in todo:
                entries[s] = self.build(s)
                if self.engine.log:
                    self.engine.log(f"  state {s}: {entries[s]['animals']} animals, K={entries[s]['fit']['selected_k']}")
        manifest = self.manifest()
        return {s: manifest["states"][s] for s in todo if s in manifest["states"]}

    def manifest(self):
        states = {}
        for state_id in self.plan:
            meta_path = self.dir / f"{state_id}.meta.json"
            if meta_path.exists():
                meta = read_json(meta_path)
                if meta.get("cache_key") == self._key(state_id):
                    states[state_id] = meta["entry"]
        value = {"schema": SCHEMA, "world": self.engine.world.name, "states": states,
                 "fold_of_training_animal": dict(sorted(self.fold.items()))}
        write_json(Path(self.engine.world.out) / "state_manifest.json", value, pretty=True)
        return value

    def get(self, state_id):
        state = self._states.get(state_id)
        if state is None:
            meta_path = self.dir / f"{state_id}.meta.json"
            if not meta_path.exists() or read_json(meta_path).get("cache_key") != self._key(state_id):
                raise EngineError(f"regional state '{state_id}' has not been fitted for these inputs; run `states`")
            model, proof = read_json(self.dir / f"{state_id}.json"), read_json(self.dir / f"{state_id}.proof.json")
            entry, plan, cohort = read_json(meta_path)["entry"], self.plan[state_id], self.engine.cohort
            bad = (set(proof["train_animals"]) | {i.split("_", 1)[0] for i in proof["train_ids"]}) & cohort.heldout_animals
            if bad:
                raise LeakError(f"regional state '{state_id}' contains validation or test animals: {sorted(bad)[:5]}")
            # parameters, proof, state manifest and the planned fitting and served sets must all agree
            if not (T._digest(model) == sha_file(self.dir / f"{state_id}.json") == proof["model_state_sha256"]
                    == entry["model_state_sha256"]):
                raise EngineError(f"regional state '{state_id}': parameters, proof and state manifest do not carry one hash")
            fitted = set(plan["train_animals"])
            planned_ids = [cohort.rows[n].full_id for n in cohort.order if cohort.rows[n].animal in fitted]
            served_ids = [cohort.rows[n].full_id for n in plan["served"]]
            if (proof.get("name") != state_id or proof["train_animals"] != sorted(fitted) or proof["train_ids"] != planned_ids
                    or model["training_ids_sha256"] != ids_sha(planned_ids)):
                raise LeakError(f"regional state '{state_id}': its proof does not list exactly the planned fitting animals and neurons")
            if proof["test_ids"] != served_ids or proof["test_animals"] != sorted({cohort.rows[n].animal for n in plan["served"]}):
                raise EngineError(f"regional state '{state_id}': its proof does not list exactly the planned served neurons")
            state = Geo.RegionalState(model, proof)
            self._states[state_id] = state
            self._sets[state_id] = (frozenset(proof["train_animals"]), frozenset(i.split("_", 1)[0] for i in proof["train_ids"]))
        return state

    def animal_sets(self, state_id):
        self.get(state_id)
        return self._sets[state_id]

    def file_hashes(self, state_id):
        """The state's identity and the bytes of its two files, as they are on disk now (pinned by manifests and the freeze)."""
        paths = {"state": self.dir / f"{state_id}.json", "proof": self.dir / f"{state_id}.proof.json"}
        if not all(p.exists() for p in paths.values()):
            raise EngineError(f"regional state '{state_id}' has not been fitted; run `states`")
        state = self.get(state_id)                      # parameters, proof, manifest entry and plan agree, or this refuses
        return {"model_state_sha256": state.proof["model_state_sha256"], "state_file_sha256": sha_file(paths["state"]),
                "proof_file_sha256": sha_file(paths["proof"]),
                "animals_sha256": ids_sha(sorted(state.proof["train_animals"])),
                "training_ids_sha256": state.model["training_ids_sha256"]}

    def refit_sha(self, state_id):
        """Fit the state again from the training rows and return its identity: a stored state whose hashes were
        all rewritten consistently cannot pass this (the fit is deterministic)."""
        plan = self.plan[state_id]
        state, _ = fit_regional_state(state_id, self.engine.cohort, self.engine.profiles, plan["train_animals"], plan["served"])
        return state.proof["model_state_sha256"]


# ----------------------------------------------------------------------------------------------
# scoring, aggregation, gate, verdict, subgroups
# ----------------------------------------------------------------------------------------------

SCORE_KEYS = ("distance_um", "f1_200")
KEPT_METRICS = ("distance_um", "f1_200", "precision_200", "recall_200", "prediction_to_truth_um",
                "truth_to_prediction_um", "cable_signed_relative_error")
CLAIM_LIMIT = "A numerical win on these animals is not a confirmatory population claim."


def score_arbor(prediction, truth):
    """Exactly the frozen geometry_metrics."""
    metrics = G.geometry_metrics(prediction, truth)
    return {key: float(metrics[key]) for key in KEPT_METRICS}


def score_arbor_checked(prediction, truth):
    """The frozen metric, and the referee's metric module on the same arrays: they must agree bit for bit."""
    mine = score_arbor(prediction, truth)
    theirs = referee().G.geometry_metrics(as_tree(prediction), as_tree(truth))
    if any(mine[key] != float(theirs[key]) for key in KEPT_METRICS):
        raise EngineError("the frozen metric and the referee's metric disagree on a forecast: nothing is scored")
    return mine


def aggregate(rows, methods):
    """rows: [{'animal', 'metrics': {method: {...}}}]. Per neuron -> mean within animal -> mean over animals."""
    animals = sorted({r["animal"] for r in rows})
    keys = sorted({k for r in rows for m in methods for k in r["metrics"][m]}) if rows else []
    per_animal = []
    for a in animals:
        mine = [r for r in rows if r["animal"] == a]
        per_animal.append({"animal": a, "neurons": len(mine),
                           "methods": {m: {k: float(np.mean([r["metrics"][m][k] for r in mine])) for k in keys}
                                       for m in methods}})
    summary = {m: {k: float(np.mean([a["methods"][m][k] for a in per_animal])) for k in keys}
               for m in methods} if per_animal else {}
    return {"animals": len(animals), "neurons": len(rows), "per_animal": per_animal, "summary": summary}


def better(summary, left, right):
    return (summary[left]["distance_um"] < summary[right]["distance_um"]
            and summary[left]["f1_200"] > summary[right]["f1_200"])


def strongest_controls(summary, controls=CONTROLS):
    controls = list(controls)
    return {"distance_um": min(controls, key=lambda k: (summary[k]["distance_um"], controls.index(k))),
            "f1_200": max(controls, key=lambda k: (summary[k]["f1_200"], -controls.index(k)))}


def engine_gate(summary, frozen_selected_policy, family="all"):
    """The declared validation rule for the candidate family (protocol section 7)."""
    if family not in FAMILIES:
        raise EngineError(f"unknown candidate family '{family}'")
    members = FAMILIES[family]
    table = {}
    for c in CANDIDATES:
        versus_controls = {k: better(summary, c, k) for k in CONTROLS}
        versus_matched = {k: better(summary, c, k) for k in MATCHED_CONTROLS.get(c, ())}
        admitted, note = True, None
        if c not in members:
            admitted, note = False, f"outside the declared candidate family '{family}'"
        elif c in FROZEN_CANDIDATES and FROZEN_POLICY_OF[c] != frozen_selected_policy:
            admitted, note = False, "the frozen prototype's own validation gate did not select this policy"
        table[c] = {"distance_um": summary[c]["distance_um"], "f1_200": summary[c]["f1_200"],
                    "conditional": c in MATCHED_CONTROLS, "uses_fitted_state": c in USES_FITTED_STATE,
                    "beats_control": versus_controls, "beats_every_control": all(versus_controls.values()),
                    "beats_every_official_copy": all(versus_controls[k] for k in OFFICIAL_COPIES),
                    "beats_matched_control": versus_matched, "beats_every_matched_control": all(versus_matched.values()),
                    "admitted": admitted, "note": note}
        table[c]["eligible"] = bool(admitted and table[c]["beats_every_control"] and table[c]["beats_every_matched_control"])
    rank = lambda c: (summary[c]["distance_um"], -summary[c]["f1_200"], CANDIDATES.index(c))
    eligible = [c for c in CANDIDATES if table[c]["eligible"]]
    chosen = min(eligible, key=rank) if eligible else None
    best_distance = min(METHODS, key=lambda m: (summary[m]["distance_um"], METHODS.index(m)))
    best_f1 = max(METHODS, key=lambda m: (summary[m]["f1_200"], -METHODS.index(m)))
    gate = {"schema": "axon_engine_gate_v2", "rule": CONFIG["gate"], "family": family, "family_members": list(members),
            "controls": list(CONTROLS), "candidates": table, "eligible": eligible, "selected": chosen,
            "verdict": "candidate_selected" if chosen else "no method beyond copying",
            # fills the referee's single candidate slot when nothing is selected; it carries no claim
            "descriptive_reference": None if chosen else min(members, key=rank),
            "strongest_control": strongest_controls(summary), "best_method_on_distance": best_distance,
            "best_method_on_f1": best_f1,
            "best_method_is_a_copying_control": {"distance_um": best_distance in CONTROLS, "f1_200": best_f1 in CONTROLS},
            "frozen_gate_selected_policy": frozen_selected_policy,
            "control_summary": {k: {s: summary[k][s] for s in SCORE_KEYS} for k in CONTROLS},
            "unseen_test_gain_guaranteed": False}
    gate["gate_sha256"] = digest(gate)
    return gate


def paired_bootstrap(deltas, seed_text, samples=10000):
    deltas = np.asarray(deltas, dtype=float)
    if len(deltas) < 2:
        return None
    rng = np.random.default_rng(int.from_bytes(hashlib.sha256(seed_text.encode()).digest()[:4], "big"))
    draws = deltas[rng.integers(0, len(deltas), size=(samples, len(deltas)))].mean(axis=1)
    return [float(v) for v in np.quantile(draws, [0.025, 0.975])]


def compare(per_animal, left, right, seed_text):
    d = np.array([a["methods"][left]["distance_um"] - a["methods"][right]["distance_um"] for a in per_animal])
    f = np.array([a["methods"][left]["f1_200"] - a["methods"][right]["f1_200"] for a in per_animal])
    return {"distance_delta_um": float(d.mean()), "f1_delta": float(f.mean()),
            "animals_better_on_distance": int((d < 0).sum()), "animals_better_on_f1": int((f > 0).sum()),
            "animals_better_on_both": int(((d < 0) & (f > 0)).sum()), "animals": int(len(d)),
            "distance_bootstrap95_um": paired_bootstrap(d, seed_text + "|d"),
            "f1_bootstrap95": paired_bootstrap(f, seed_text + "|f"), "bootstrap": "paired over animals, descriptive only"}


def _interval(values, digits):
    return "no interval" if not values else f"paired 95% bootstrap interval {values[0]:+.{digits}f} to {values[1]:+.{digits}f}"


def margin_text(selected, strongest, versus):
    """Paired differences against the strongest control on each measure, with their intervals."""
    d, f = versus[strongest["distance_um"]], versus[strongest["f1_200"]]
    return (f"Against the strongest control on distance ({strongest['distance_um']}): {d['distance_delta_um']:+.1f} um "
            f"({_interval(d.get('distance_bootstrap95_um'), 1)}). Against the strongest control on F1 "
            f"({strongest['f1_200']}): {f['f1_delta']:+.4f} ({_interval(f.get('f1_bootstrap95'), 4)}). " + CLAIM_LIMIT)


def verdict_for(result, selected, *, context):
    """The success rule on one set scored by this engine (validation; development for information)."""
    summary, per_animal = result["summary"], result["per_animal"]
    strongest = strongest_controls(summary)
    out = {"context": context, "selected_candidate": selected, "animals": result["animals"], "neurons": result["neurons"],
           "unit": "animals are brain sample IDs (ID prefix); see the protocol for the caveat", "scored_by": "engine",
           "strongest_control": strongest,
           "control_summary": {k: {s: summary[k][s] for s in SCORE_KEYS} for k in CONTROLS}}
    if selected is None:
        out.update({"success": False, "verdict": "no method beyond copying",
                    "statement": "No candidate passed the validation gate, so no method beyond copying is claimed. "
                                 "Controls and candidates are reported descriptively."})
        return out
    versus = {k: {"beats_on_both": better(summary, selected, k), **compare(per_animal, selected, k, f"{context}|{selected}|{k}")}
              for k in CONTROLS}
    matched = {k: {"beats_on_both": better(summary, selected, k), **compare(per_animal, selected, k, f"{context}|{selected}|{k}")}
               for k in MATCHED_CONTROLS.get(selected, ()) if k not in CONTROLS}
    beats_controls = all(v["beats_on_both"] for v in versus.values())
    beats_matched = all(v["beats_on_both"] for v in matched.values())
    out.update(_verdict_text(selected, beats_controls, beats_matched, [k for k, v in {**versus, **matched}.items() if not v["beats_on_both"]],
                             bool(matched), margin_text(selected, strongest, versus)))
    out.update({"selected_summary": {s: summary[selected][s] for s in SCORE_KEYS},
                "versus_control": versus, "versus_matched_control": matched,
                "beats_every_copying_control": beats_controls, "beats_every_matched_control": beats_matched,
                "conditional_candidate": selected in MATCHED_CONTROLS})
    return out


def _verdict_text(selected, beats_controls, beats_matched, failed, has_matched, margins):
    success = bool(beats_controls and beats_matched)
    if success:
        verdict = "beats copying on both measures"
        statement = (f"{selected} has lower animal-mean distance AND higher animal-mean F1 than every predeclared copying "
                     "control" + (" and than its matched controls" if has_matched else "") + " on these animals. " + margins)
    elif beats_controls:
        verdict = "rule not met: better than every copying control, not better than its matched control"
        statement = (f"{selected} is better on both measures than every copying control, but NOT than its matched control(s) "
                     + ", ".join(failed) + ". By the predeclared rule this is not a success: nothing is shown for the "
                     "conditional weights; whether an unconditional local template beats copying was not the frozen question. "
                     + margins)
    else:
        verdict = "does not beat copying"
        statement = f"{selected} FAILS the success rule: it is not better on both measures than: " + ", ".join(failed) + ". " + margins
    return {"success": success, "verdict": verdict, "statement": statement}


def verdict_from_referee(result, roles, selected, *, context="test"):
    """The success rule read from the independent referee's own result (its animal means, its booleans, its
    paired differences). `selected` is the frozen selection; None means the validation gate selected nothing
    and the candidate slot holds the descriptive reference, for which nothing can be claimed."""
    overall = result["overall"]
    summary, candidate = overall["summary"], roles["candidate"]
    copies, matched = list(roles["copy_baselines"]), list(roles["matched_controls"])
    controls = [m for m in CONTROLS if m in copies + matched]
    strongest = strongest_controls(summary, controls)
    out = {"context": context, "selected_candidate": selected, "referee_candidate_slot": candidate,
           "animals": overall["independent_animals"], "neurons": overall["neurons"], "scored_by": "independent referee (score_once)",
           "unit": "animals are brain sample IDs (ID prefix); see the protocol for the caveat",
           "referee_roles": roles, "strongest_control": strongest,
           "control_summary": {k: {s: summary[k][s] for s in SCORE_KEYS} for k in controls},
           "candidate_summary": {s: summary[candidate][s] for s in SCORE_KEYS},
           "referee_beats_all_copying": bool(overall["beats_all_copying"]),
           "referee_beats_all_matched": bool(overall["beats_all_matched"]),
           "versus": {m: {k: v for k, v in overall["against"][m].items() if k != "animal_differences"} for m in copies + matched}}
    beat = lambda m: better(summary, candidate, m)
    if [m for m in copies if not beat(m)] and overall["beats_all_copying"]:
        raise EngineError("the referee's verdict and its own animal means disagree")
    if selected is None:
        out.update({"success": False, "verdict": "no method beyond copying",
                    "statement": ("The validation gate selected no candidate, so no method beyond copying is claimed, whatever "
                                  f"this table shows. {candidate} only fills the referee's candidate slot as the predeclared "
                                  "descriptive reference.")})
        return out
    if candidate != selected:
        raise EngineError("the referee's candidate is not the frozen selection")
    beats_controls = bool(overall["beats_all_copying"]) and all(beat(m) for m in ATTACHMENT_CONTROLS)
    own_matched = [m for m in matched if m not in ATTACHMENT_CONTROLS]
    beats_matched = all(beat(m) for m in own_matched)
    failed = [m for m in copies + matched if not beat(m)]
    out.update(_verdict_text(selected, beats_controls, beats_matched, failed, bool(own_matched),
                             margin_text(selected, strongest, overall["against"])))
    if out["success"] != bool(overall["beats_all_copying"] and overall["beats_all_matched"]):
        raise EngineError("the engine's reading of the rule differs from the referee's booleans")
    out.update({"beats_every_copying_control": beats_controls, "beats_every_matched_control": beats_matched,
                "conditional_candidate": selected in MATCHED_CONTROLS})
    return out


def subgroup_report(rows, methods, selected, *, controls=None, sparse_below=5):
    """Per animal always; by source label; by soma nucleus. Every row: neurons AND animals. Frozen spec.
    Every subgroup row is descriptive; rows with fewer than `sparse_below` animals are marked sparse."""
    controls = [m for m in (CONTROLS if controls is None else controls) if m in methods]
    whole = aggregate(rows, methods)
    label_of = {r["animal"]: r["source_label"] for r in rows}
    per_animal = [{"animal": a["animal"], "source_label": label_of[a["animal"]], "neurons": a["neurons"],
                   "methods": {m: {k: a["methods"][m][k] for k in SCORE_KEYS} for m in methods}}
                  for a in whole["per_animal"]]

    def group(field, value_of):
        out = []
        for value in sorted({value_of(r) for r in rows}):
            mine = [r for r in rows if value_of(r) == value]
            agg = aggregate(mine, methods)
            sparse = agg["animals"] < sparse_below
            entry = {field: value, "neurons": agg["neurons"], "animals": agg["animals"],
                     "descriptive_only": True, "sparse_descriptive_only": sparse,
                     "interpretation": "sparse_descriptive_only" if sparse else "descriptive_only",
                     "note": (f"fewer than {sparse_below} animals: descriptive only, no population claim" if sparse
                              else "descriptive, no population claim"),
                     "methods": {m: {k: agg["summary"][m][k] for k in SCORE_KEYS} for m in methods}}
            if field == "soma_label":
                kind = soma_label_kind(value if value != "unknown" else "")
                entry["soma_label_kind"] = kind
                entry["nucleus_assigned"] = kind == "nucleus"
            if selected is not None and selected in methods:
                wins = {k: better(agg["summary"], selected, k) for k in controls}
                entry["descriptive_comparison"] = {"selected_better_than_every_control_on_both": all(wins.values()),
                                                   "controls_not_beaten": [k for k, v in wins.items() if not v]}
            else:
                entry["descriptive_comparison"] = None
            out.append(entry)
        return out

    return {"spec_sha256": SUBGROUP_SPEC_SHA, "selected_candidate": selected, "sparse_below_animals": sparse_below,
            "interpretation": "every subgroup row is descriptive and carries no population claim",
            "unit": "animals are brain sample IDs (ID prefix)", "per_animal": per_animal,
            "by_source_label": group("source_label", lambda r: r["source_label"] or "unknown"),
            "by_soma_label": group("soma_label", lambda r: r["soma_label"] or "unknown")}


# ----------------------------------------------------------------------------------------------
# the referee: preflight on a three-donor fixture (its real functions, a temporary directory)
# ----------------------------------------------------------------------------------------------


def referee_preflight():
    """Can the referee on disk seal a bundle that honestly declares an upstream fit on TRAINING animals, and
    does it refuse one whose forecasts name a held-out animal? Its own four functions are run on a tiny
    synthetic fixture; nothing of the cohort is involved. The answer is cached per referee source hash."""
    import tempfile
    key = digest(referee_hashes())
    if _REFEREE.get("preflight_key") == key:
        return _REFEREE["preflight"]
    Ref = referee()

    def line(soma, axis):
        xyz = np.zeros((3, 3))
        xyz[:, axis] = [0.0, 50.0, 100.0]
        return xyz + np.array(soma), np.array([-1, 0, 1]), np.array([1, 2, 2])

    def attempt(declared, manifest_level, row_level):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp) / "bundle"
            try:
                train, donors = [], []
                for name, animal, soma, tag, axis in (("d1", "A", [0.0, 0.0, 0.0], "P", 0), ("d2", "A", [0.0, 100.0, 0.0], "P", 1),
                                                      ("d3", "B", [200.0, 0.0, 0.0], "Q", 1)):
                    tree, relative = line(soma, axis), f"train/{name}.npz"
                    Ref.save_tree(base / relative, tree)
                    row = {"id": name, "animal": animal, "soma_um": soma, "peptide": tag, "arbor_path": relative,
                           "arbor_sha256": Ref.sha(base / relative)}
                    train.append(row)
                    donors.append({**row, "arbor": tree})
                queries = [{"id": "q1", "animal": "T1", "soma_um": [10.0, 0.0, 0.0], "peptide": "P", "nucleus": "N"},
                           {"id": "q2", "animal": "T2", "soma_um": [190.0, 0.0, 0.0], "peptide": None, "nucleus": None}]
                method = base / "method/fixture.txt"
                method.parent.mkdir(parents=True, exist_ok=True)
                method.write_text("referee preflight fixture\n")
                models = ["candidate", *Ref.COPIES, "matched"]
                declaration = {
                    "geometry_contract": {"units": "um", "axes": ["AP", "DV", "ML"], "frame": "fixture", "admission": "fixture",
                                          "preprocessing": "fixture", "clipping": "none"},
                    "method_artifacts": {"method/fixture.txt": Ref.sha(method)},
                    "splits": {"train": train, "validation": [{"id": "v1", "animal": "V", "soma_um": [0.0, 0.0, 0.0]}], "test": queries},
                    "excluded_ids": [], "upstream_fit_animals": list(declared),
                    "models": {"candidate": "candidate", "copy_baselines": list(Ref.COPIES), "matched_controls": ["matched"]},
                    "config": dict(Ref.CONFIG)}
                Ref.freeze_declaration(base, declaration)
                copies = Ref.copying_baselines(donors, queries)
                rows = []
                for q in queries:
                    for model in models:
                        source = model if model in Ref.COPIES else "medoid"
                        relative = f"predictions/{model}/{q['id']}.npz"
                        Ref.save_tree(base / relative, copies["predictions"][source][q["id"]])
                        rows.append({"query_id": q["id"], "model": model, "path": relative, "sha256": Ref.sha(base / relative),
                                     "donor_ids": [copies["donor_proof"][source][q["id"]]],
                                     "upstream_fit_animals": list(row_level) if model == "candidate" else []})
                Ref.seal_predictions(base, {"declaration_sha256": Ref.sha(base / "declaration.json"), "predictions": rows,
                                            "upstream_fit_animals": list(manifest_level)})
                return {"sealed": True, "error": None}
            except Exception as error:             # the exact refusal is the finding
                return {"sealed": False, "error": f"{type(error).__name__}: {error}"}

    plain = attempt([], [], [])
    honest = attempt(["A", "B"], ["A", "B"], ["A", "B"])
    leak = attempt([], ["T1"], ["T1"])
    out = {"referee_sha256": referee_hashes(),
           "no_upstream_fit_sealable": plain["sealed"], "no_upstream_fit_error": plain["error"],
           "fitted_state_sealable": honest["sealed"], "fitted_state_error": honest["error"],
           "heldout_upstream_refused": not leak["sealed"], "heldout_upstream_error": leak["error"],
           "what": ("three attempts on a fixture: no upstream fit; an upstream fit on the two TRAINING animals declared at all three "
                    "levels; a TEST animal named in the manifest and candidate rows"),
           "written_up_in": REFEREE_FINDING}
    _REFEREE["preflight_key"], _REFEREE["preflight"] = key, out
    return out


# ----------------------------------------------------------------------------------------------
# the engine
# ----------------------------------------------------------------------------------------------


def consensus_prune(template, others, weights):
    """Keep the template's cable that enough neighbours share (declared radius and support), plus the
    path from each kept node back to the soma. Vertices are never moved; when too little would remain
    the template is returned whole. Uses donors' arbors only."""
    from scipy.spatial import cKDTree
    spec = CONFIG["engine"]
    weights = np.asarray(weights, dtype=float)
    n = len(template.xyz)
    if not others or weights.sum() <= 0 or n < 2:
        return template, {"consensus_kept_cable_fraction": 1.0, "fallback_reason": "consensus_without_neighbours"}
    support = np.zeros(n)
    for arbor, w in zip(others, weights):
        support += w * (cKDTree(arbor.xyz).query(template.xyz)[0] <= spec["consensus_radius_um"])
    keep = support / weights.sum() >= spec["consensus_min_support"]
    keep[0] = True
    parent = template.parent
    for i in range(n - 1, 0, -1):
        if keep[i]:
            keep[parent[i]] = True
    lengths = np.r_[0.0, np.linalg.norm(template.xyz[1:] - template.xyz[parent[1:]], axis=1)]
    total, kept = float(lengths.sum()), float(lengths[keep].sum())
    info = {"consensus_kept_cable_fraction": kept / total if total > 0 else 1.0, "consensus_kept_cable_um": kept}
    if kept < spec["consensus_min_cable_um"]:
        return template, {**info, "fallback_reason": "consensus_cable_below_minimum"}
    remap = np.cumsum(keep) - 1
    index = np.flatnonzero(keep)
    pruned = G.Arbor(template.xyz[index].copy(), np.r_[-1, remap[parent[index[1:]]]].astype(int), template.types[index].copy())
    return pruned.validate(), info


def _predict_init(world_json, stage, options):
    _WORKER["engine"] = Engine(World.from_json(world_json), stage, workers=1, log=None, **options)
    _WORKER["engine"].prepare_prediction(verify_files=False)


def _predict_task(neurons):
    engine = _WORKER["engine"]
    return [engine.forecast_and_save(n) for n in neurons]


def _ingest_task(args):
    path, frame, hemisphere, recon_type, table_soma, expected_sha, out_path = args[:7]
    alternates = args[7] if len(args) > 7 else ()
    log, arbor = ingest_trace(path, frame=frame, hemisphere=hemisphere, recon_type=recon_type, table_soma=table_soma,
                              expected_sha=expected_sha)
    if arbor is not None:
        G.save_arbor(out_path, arbor)
        log["canonical_file_sha256"] = sha_file(out_path)
    if alternates:            # a neuron with more than one copy on disk: the copies are checked against each other
        others = []
        for alt_path, alt_frame, alt_sha in alternates:
            alt_log, _ = ingest_trace(alt_path, frame=alt_frame, hemisphere=hemisphere, recon_type=recon_type,
                                      table_soma=table_soma, expected_sha=alt_sha)
            others.append({"file": Path(alt_path).name, "frame": alt_frame, "status": alt_log["status"],
                           "rules_violated": alt_log["rules_violated"], "arbor_sha256": alt_log.get("arbor_sha256")})
        log["copy_check"] = {"copies": 1 + len(others), "others": others,
                             "all_copies_give_the_same_canonical_arbor": all(
                                 o["status"] == log["status"] and o["arbor_sha256"] == log.get("arbor_sha256") for o in others)}
    return log


def alternate_copies(store, neuron, chosen):
    return tuple((c["path"], c["frame"], c["sha256"]) for c in store.copies.get(neuron, [])
                 if c is not chosen and not c.get("unavailable"))


def _score_task(args):
    truth_path, items = args
    truth = G.load_arbor(truth_path)
    return [(model, score_arbor_checked(G.load_arbor(path), truth)) for model, path in items]


class Engine:
    """One world, one stage. Lazy parts: library (training arbors), all-pairs cache, profiles (training
    regional rows), states (regional states)."""

    def __init__(self, world, stage, *, workers=1, log=print, develop_per_animal=None, family="all"):
        if family not in FAMILIES:
            raise EngineError(f"unknown candidate family '{family}' (declared: {sorted(FAMILIES)})")
        self.world, self.stage, self.workers, self.log = world, stage, max(1, int(workers)), log
        self.develop_per_animal = develop_per_animal      # development only: a fixed subsample of queries per animal
        self.family = family                              # declared when validation forecasts are made; pinned by the freeze
        self.out = Path(world.out)
        self.cohort = Cohort(world)
        self.store = TraceStore(world, self.cohort, stage)
        self.states = StateBook(self)
        self.memo = FrozenMemo(lambda: self.cache)
        self._profiles = self._profile_audit = None
        self._library = self._cache = None
        self._donors, self._scopes = {}, {}
        self.stage_dir = self.out / stage

    # -- small helpers ---------------------------------------------------------------------
    def say(self, text):
        if self.log:
            self.log(text)

    @staticmethod
    def can_parallel():
        module = sys.modules.get("axon_engine")
        return module is not None and getattr(module, "__file__", None) == __file__

    def pool(self, workers, initializer=None, initargs=()):
        import multiprocessing
        return multiprocessing.get_context("spawn").Pool(workers, initializer=initializer, initargs=initargs)

    def options(self):
        return {"develop_per_animal": self.develop_per_animal, "family": self.family}

    @property
    def profiles(self):
        if self._profiles is None:
            self._profiles, self._profile_audit = load_training_profiles(self.world, self.cohort)
        return self._profiles

    @property
    def profile_audit(self):
        _ = self.profiles
        return self._profile_audit

    def label(self):
        if self.world.kind == "synthetic":
            return f"SYNTHETIC world '{self.world.name}', stage {self.stage}: not biological evidence"
        if self.stage == "develop" and self.world.name.endswith("_local37"):
            return ("DEVELOPMENT on the 37 training traces of data/raw/jiao2025/swc (leave one animal out): "
                    "not a test, no success claim")
        return {"develop": "DEVELOPMENT on training animals (leave one animal out): not a test, no success claim",
                "validate": "VALIDATION (selection of the candidate; not the test)",
                "test": "INDEPENDENT TEST (run once)"}[self.stage]

    def protocol_state(self):
        protocol = self.world.protocol
        text = Path(protocol).read_text() if protocol and os.path.exists(protocol) else ""
        return {"path": protocol, "sha256": sha_bytes(text.encode()) if text else None, "draft": (not text) or "STATUS: DRAFT" in text}

    # -- library ---------------------------------------------------------------------------
    def training_listing(self):
        """[neuron, file, sha256, other copies] for every training neuron listed on disk: the chosen copy and, when
        there are more, the others (they are checked against it at ingestion). Opens no trace."""
        out = []
        for neuron in self.cohort.part_neurons["train"]:
            copy = self.store.choose(neuron)
            if copy is not None:
                others = sorted([c["file"], c["sha256"]] for c in self.store.copies[neuron] if c is not copy)
                out.append([neuron, copy["file"], copy["sha256"], others])
        return out

    def build_library(self):
        """Ingest every training trace on disk into the canonical library (cached by file hash)."""
        self.store.require(self.stage, need_truth=False)
        directory = self.out / "library"
        directory.mkdir(parents=True, exist_ok=True)
        log_path = self.out / "ingestion_log.json"
        code = engine_code_sha().get("analysis/axon_engine.py")
        previous = read_json(log_path) if log_path.exists() else {}
        entries = previous.get("neurons", {}) if previous.get("engine_sha256") == code else {}
        jobs, kept = [], {}
        for neuron in self.cohort.part_neurons["train"]:
            copy = self.store.choose(neuron)
            if copy is None:
                continue
            row = self.cohort.rows[neuron]
            old = entries.get(neuron)
            if copy.get("unavailable"):
                kept[neuron] = unavailable_log(neuron, row, copy)
                continue
            self.store.open(neuron, "donor")
            target = directory / f"{neuron}.npz"
            others = sorted([c["file"], c["sha256"]] for c in self.store.copies[neuron] if c is not copy)
            if old and old.get("manifest_sha256") == copy["sha256"] and old.get("source_file") == copy["file"] and (
                    old["status"] == "excluded" or target.exists()) and old.get("other_copies", []) == others:
                kept[neuron] = old
                continue
            jobs.append((neuron, copy, (copy["path"], copy["frame"], row.hemisphere, row.recon_type, list(row.soma),
                                        copy["sha256"], str(target), alternate_copies(self.store, neuron, copy))))
        if jobs:
            self.say(f"ingesting {len(jobs)} training traces ({len(kept)} cached)")
            if self.workers > 1 and len(jobs) > 8 and self.can_parallel():
                with self.pool(self.workers) as pool:
                    logs = pool.map(_ingest_task, [j[2] for j in jobs], chunksize=8)
            else:
                logs = [_ingest_task(j[2]) for j in jobs]
            for (neuron, copy, _), log in zip(jobs, logs):
                row = self.cohort.rows[neuron]
                log.update({"neuron": neuron, "animal": row.animal, "part": row.part, "source_file": copy["file"],
                            "manifest_sha256": copy["sha256"], "copies_on_disk": len(self.store.copies[neuron]),
                            "other_copies": sorted([c["file"], c["sha256"]] for c in self.store.copies[neuron] if c is not copy)})
                kept[neuron] = log
        kept = dict(sorted(kept.items()))
        record = {"schema": SCHEMA, "world": self.world.name, "engine_sha256": code, "rules": RULE_TEXT,
                  "listing_sha256": digest(self.training_listing()), "neurons": kept}
        write_json(log_path, record, pretty=True)
        self._library = self._cache = None
        self._donors, self._scopes = {}, {}
        self.memo = FrozenMemo(lambda: self.cache)
        return self.library_summary()

    def library_problem(self):
        """Why the stored library is not the library of the traces on disk now (None when it is)."""
        log_path = self.out / "ingestion_log.json"
        if not log_path.exists():
            return "the training library has not been built"
        record = read_json(log_path)
        if record.get("engine_sha256") != engine_code_sha().get("analysis/axon_engine.py"):
            return "the training library was built by another version of the engine"
        listing = self.training_listing()
        if record.get("listing_sha256") != digest(listing) or sorted(record["neurons"]) != [row[0] for row in listing]:
            return (f"the training library holds {len(record['neurons'])} neurons but {len(listing)} training traces are listed "
                    "on disk, or their files changed")
        for neuron, file, sha, _ in listing:
            entry = record["neurons"][neuron]
            if entry.get("source_file") != file or entry.get("manifest_sha256") != sha:
                return f"the library entry of {neuron} was made from another file"
            if entry["status"] == "included" and not (self.out / "library" / f"{neuron}.npz").exists():
                return f"the canonical arbor of {neuron} is missing"
        return None

    def _load_library(self):
        problem = self.library_problem()
        if problem:
            raise EngineError(problem + "; run `ingest`")
        record = read_json(self.out / "ingestion_log.json")
        included = {n: e for n, e in record["neurons"].items() if e["status"] == "included"}
        for neuron in included:
            if self.cohort.rows[neuron].part != "train":
                raise LeakError(f"{neuron}: a held-out neuron is in the donor library")
        ids = sorted(included)
        rows = [self.cohort.rows[n] for n in ids]
        return {"ids": ids, "entries": included, "all": record["neurons"],
                "soma": np.array([r.soma for r in rows], dtype=float).reshape(-1, 3),
                "animal": np.array([r.animal for r in rows]), "label": np.array([r.source_label for r in rows]),
                "index": {n: i for i, n in enumerate(ids)},
                "listing_sha256": record["listing_sha256"],
                "fingerprint": digest([[n, included[n]["arbor_sha256"]] for n in ids]),
                "training_ids_sha256": ids_sha(ids)}

    @property
    def library(self):
        if self._library is None:
            self._library = self._load_library()
        return self._library

    def ensure_library(self):
        """The library must be the one of the traces on disk; a stale or missing one is rebuilt (cached entries
        make this cheap), never used."""
        problem = self.library_problem()
        if problem:
            self.say(f"{problem}: ingesting")
            self.build_library()
        return self.library

    def library_summary(self):
        lib = self.library
        excluded = {n: e for n, e in lib["all"].items() if e["status"] != "included"}
        by_rule = Counter(rule for e in excluded.values() for rule in e["rules_violated"])
        coding = Counter(str(e.get("axon_coding")) for e in lib["entries"].values())
        return {"training_traces_listed": len(lib["all"]), "included": len(lib["ids"]), "excluded": len(excluded),
                "excluded_by_rule": dict(sorted(by_rule.items())),
                "exclusions": {n: {"rules": e["rules_violated"], "details": e["details"]} for n, e in sorted(excluded.items())},
                "animals": len(set(lib["animal"].tolist())), "axon_coding": dict(sorted(coding.items())),
                "cable_outside_atlas_box": outside_summary(lib["entries"].values()),
                "connector_geometry": connector_summary(lib["entries"].values()),
                "unavailable_at_source": rule8_report(lib["all"], self.cohort, "train"),
                "second_copies": {"neurons_with_more_than_one_copy": sum(1 for e in lib["all"].values() if e.get("copy_check")),
                                  "copies_that_disagree": sorted(n for n, e in lib["all"].items() if e.get("copy_check")
                                                                 and not e["copy_check"]["all_copies_give_the_same_canonical_arbor"])},
                "listing_sha256": lib["listing_sha256"],
                "fingerprint": lib["fingerprint"], "training_ids_sha256": lib["training_ids_sha256"]}

    # -- all-pairs cache ---------------------------------------------------------------------
    def pair_cache(self, siblings=()):
        lib = self.library
        return PairCache(self.out / "paircache", self.out / "library", lib["ids"],
                         {n: lib["entries"][n]["arbor_sha256"] for n in lib["ids"]}, siblings=siblings)

    def build_cache(self, *, minutes=None, siblings=()):
        self.ensure_library()
        cache = self.pair_cache(siblings)
        status = cache.build(workers=self.workers, pool_factory=self.pool if self.can_parallel() else None,
                             minutes=minutes, log=self.log)
        self._cache = None
        self.memo = FrozenMemo(lambda: self.cache)
        for neuron, record in self._donors.items():
            self.memo.register(neuron, record["arbor"])
        return status

    def open_cache(self, *, verify_files=True):
        if self._cache is None:
            self._cache = self.pair_cache().load(verify_files=verify_files)
        return self._cache

    @property
    def cache(self):
        return self.open_cache()

    # -- donors and donor pools ----------------------------------------------------------------
    def donor(self, neuron):
        """A donor record in the frozen prototype's format (training animals only)."""
        record = self._donors.get(neuron)
        if record is None:
            row = self.cohort.rows[neuron]
            if row.part != "train":
                raise LeakError(f"{neuron}: a {row.part} neuron cannot be a donor")
            entry = self.library["entries"][neuron]
            arbor = G.load_arbor(self.out / "library" / f"{neuron}.npz")
            if arbor_sha(arbor) != entry["arbor_sha256"]:
                raise EngineError(f"library arbor of {neuron} changed on disk")
            self.memo.register(neuron, arbor)
            record = {"neuron": neuron, "animal": row.animal, "full_id": row.full_id, "tag": row.tag,
                      "source_label": row.source_label, "soma_label": row.soma_label,
                      "table_soma": np.array(row.soma, dtype=float), "swc_soma": arbor.xyz[0].copy(),
                      "profile": self.profiles[row.full_id], "arbor": arbor}
            self._donors[neuron] = record
        return record

    def scope(self, exclude_animal):
        """The donor pool of a query: the whole training library, minus the query's own animal when the
        query is itself a training animal (development). No panel, no sampling."""
        scope = self._scopes.get(exclude_animal)
        if scope is None:
            lib = self.library
            keep = np.ones(len(lib["ids"]), dtype=bool) if exclude_animal is None else lib["animal"] != exclude_animal
            ids = [n for n, k in zip(lib["ids"], keep) if k]
            if not ids:
                raise EngineError("no donor is available for this query")
            animals = lib["animal"][keep].tolist()
            scope = {"ids": ids, "soma": lib["soma"][keep], "label": lib["label"][keep], "animal": animals,
                     "animal_set": frozenset(animals), "ids_sha256": ids_sha(ids),
                     "animals_sha256": ids_sha(sorted(set(animals))), "exclude_animal": exclude_animal,
                     "pairs": len(ids) * (len(ids) - 1) // 2, "medoids": {}, "donors": None}
            self._scopes[exclude_animal] = scope
        return scope

    def pool_donors(self, scope):
        if scope["donors"] is None:
            scope["donors"] = [self.donor(n) for n in scope["ids"]]
        return scope["donors"]

    # -- the four official copies: the engine's replica of the contract's rules ---------------------
    def replica_medoid(self, scope, label):
        """Index (in the pool) of the animal-balanced geometric medoid of the pool, or of its neurons with the
        given source label: argmin of (distances @ animal weights), the contract's arithmetic, on the pinned cache."""
        index = scope["medoids"].get(label)
        if index is None:
            D, _ = self.memo.matrices(scope["ids"])
            members = list(range(len(scope["ids"]))) if label is None else np.flatnonzero(scope["label"] == label).tolist()
            animals = [scope["animal"][i] for i in members]
            counts = Counter(animals)
            weights = np.array([1 / counts[a] for a in animals])
            weights /= weights.sum()
            distances = D[np.ix_(members, members)]
            index = members[int(np.argmin(distances @ weights))]
            scope["medoids"][label] = index
        return index

    def replica_copies(self, row, scope):
        soma = np.array(row.soma, dtype=float)
        d2 = ((scope["soma"] - soma) ** 2).sum(axis=1)
        nearest = int(np.argmin(d2))
        medoid = self.replica_medoid(scope, None)
        same = np.flatnonzero(scope["label"] == row.source_label) if row.source_label else np.array([], dtype=int)
        if len(same):
            label_nearest, label_medoid, fallback = int(same[int(np.argmin(d2[same]))]), self.replica_medoid(scope, row.source_label), False
        else:
            label_nearest, label_medoid, fallback = nearest, medoid, True
        return {"nearest": nearest, "medoid": medoid, "peptide_nearest": label_nearest, "peptide_medoid": label_medoid,
                "fallback": fallback, "d2": d2, "same_label_donors": int(len(same))}

    def known_row(self, neuron):
        """What the referee may know about a neuron: identity, table soma, source label, soma label."""
        row = self.cohort.rows[neuron]
        return {"id": neuron, "animal": row.animal, "soma_um": [float(v) for v in row.soma],
                "peptide": row.source_label or None, "nucleus": row.soma_label or None}

    def referee_copies(self, scope, neurons):
        """The referee's own copying_baselines on this donor pool and these queries (known metadata only)."""
        Ref = referee()
        donors = [{**self.known_row(n), "arbor": as_tree(self.donor(n)["arbor"])} for n in scope["ids"]]
        return Ref.copying_baselines(donors, [self.known_row(n) for n in neurons])

    # -- isolation -------------------------------------------------------------------------
    def assert_isolated(self, row, state_id, scope):
        """Refuse when the query's animal is in the regional state or among the donors, or when a
        held-out animal is in either. Runs for every forecast, whatever the method."""
        animals, id_animals = self.states.animal_sets(state_id)
        if row.animal in animals or row.animal in id_animals:
            raise LeakError(f"query {row.neuron}: its animal {row.animal} is in regional state '{state_id}'")
        held = self.cohort.heldout_animals
        if animals & held or id_animals & held:
            raise LeakError(f"regional state '{state_id}' contains validation or test animals")
        if row.animal in scope["animal_set"]:
            raise LeakError(f"query {row.neuron}: its animal {row.animal} is in the donor set")
        if scope["animal_set"] & held:
            raise LeakError("the donor set contains validation or test animals")

    # -- candidates ------------------------------------------------------------------------
    def attach(self, arbor, soma):
        """Frozen local attachment at the first supported radius of the ladder; else translation."""
        tried = []
        for radius in CONFIG["engine"]["attachment_ladder_um"]:
            try:
                out, details = T.compact_attachment(arbor, soma, radius)
            except T.UnsupportedAttachment as unsupported:
                tried.append({"radius_um": radius, "reason": unsupported.reason})
                continue
            return out, {"executed": f"local_attachment_{int(radius)}", "radius_um": radius,
                         "soma_shift_um": details["soma_shift_um"], "edge_ratio_min": details.get("edge_ratio_min"),
                         "edge_ratio_max": details.get("edge_ratio_max"), "unsupported": tried,
                         "distal_vertices_kept": details.get("distal_vertex_count")}
        shift = float(np.linalg.norm(np.asarray(soma) - arbor.xyz[0]))
        return arbor.shifted(soma), {"executed": "translation", "radius_um": None, "soma_shift_um": shift, "unsupported": tried}

    def forecast(self, neuron):
        """Every control and candidate for one query, from permitted inputs only.

        Inputs taken from the query: table soma, source label, table tag, identity. Nothing else.
        Returns a list of dicts (model, arbor, provenance) in the order of METHODS, and the prediction-log entry.
        """
        row = self.cohort.rows[neuron]
        soma = np.array(row.soma, dtype=float)
        scope = self.scope(row.animal if row.part == "train" else None)
        state_id = self.states.state_id_for(neuron)
        regional = self.states.get(state_id)
        self.assert_isolated(row, state_id, scope)
        query = T.Query(row.neuron, row.animal, row.full_id, soma.copy(), soma.copy(), row.tag)
        pool = self.pool_donors(scope)
        out = {}

        def put(model, arbor, **info):
            out[model] = {"model": model, "arbor": arbor, **info}

        # the four official copies (replica of the contract's rules; `predict` sets them against the referee's API)
        replica = self.replica_copies(row, scope)
        d2 = replica["d2"]
        for model in OFFICIAL_COPIES:
            i = replica[model]
            label_rule = model.startswith("peptide_")
            fallback = bool(label_rule and replica["fallback"])
            put(model, pool[i]["arbor"].shifted(soma), family="official_copy", requested_policy=model,
                executed_policy=model[len("peptide_"):] if fallback else model, fallback=fallback,
                fallback_reason="no_training_donor_with_the_query_source_label" if fallback else None,
                selected_donor=scope["ids"][i], donor_ids=[scope["ids"][i]],
                details={"soma_distance_um": float(np.sqrt(d2[i])), "same_label_donors": replica["same_label_donors"],
                         "pool_donors": len(scope["ids"])})

        # frozen prototype, unmodified predict_policy on the whole pool (C1, C4, C5 controls; F1, F2 candidates)
        with self.memo.frozen_memo():
            for policy in T.POLICIES:
                arbor, info = T.predict_policy(query, pool, regional, policy)
                attachment = info.get("attachment") or {}
                put(METHOD_OF_FROZEN_POLICY[policy], arbor, family="frozen", requested_policy=policy,
                    executed_policy=info["executed_policy"], fallback=bool(info["fallback"]) and policy != T.BASELINE,
                    fallback_reason=info["fallback_reason"], selected_donor=info["selected_donor"],
                    donor_ids=[info["selected_donor"]],
                    details={"baseline_donor": info["baseline_donor"], "attempted_donor": info.get("attempted_donor"),
                             "radius_um": info.get("radius_um"), "soma_shift_um": attachment.get("soma_shift_um"),
                             "edge_ratio_min": attachment.get("edge_ratio_min"), "edge_ratio_max": attachment.get("edge_ratio_max"),
                             "pool_donors": len(pool)})
                del info

        # engine family: weighted geometric medoid among the K nearest donors of the whole library
        spec = CONFIG["engine"]
        n = len(scope["ids"])
        k = min(spec["k_neighbours"], n)
        nearest_k = np.lexsort((np.arange(n), d2))[:k]
        index = np.sort(nearest_k)                                   # ascending neuron ID: the tie rule
        donors = [pool[i] for i in index]
        ids = [d["neuron"] for d in donors]
        D, F = self.cache.sub(ids)
        distance = np.sqrt(d2[index])
        animals = [d["animal"] for d in donors]
        count = Counter(animals)
        prior = np.array([1.0 / count[a] for a in animals]) / len(count)
        h = max(spec["kernel_h0_um"], float(np.sort(distance)[min(spec["kernel_rank"], k) - 1]))
        kernel = np.exp(-distance ** 2 / (2.0 * h * h))
        same = np.array([d["source_label"] == row.source_label for d in donors])
        label = np.where(same, 1.0, spec["label_other_factor"])
        Geo.assert_query_safe(query, regional)
        posterior = Geo.state_posterior(regional.model, soma[None], [row.tag])[0]
        responsibilities = Geo.state_posterior(regional.model, np.array([d["table_soma"] for d in donors]),
                                               [d["tag"] for d in donors],
                                               measured_profiles=np.array([d["profile"] for d in donors]))
        regional_weights = T.conditional_weights(animals, responsibilities, posterior)["weights"]
        raw = {"E1_local_uniform": prior, "E2_local_kernel": prior * kernel, "E3_local_kernel_label": prior * kernel * label,
               "E4_local_kernel_regional": regional_weights * kernel,
               "E5_local_kernel_label_regional": regional_weights * kernel * label}
        weights = {name: w / w.sum() for name, w in raw.items()}
        choice = {name: T.weighted_medoid(D, w)[0] for name, w in weights.items()}
        w2 = weights["E2_local_kernel"]
        choice["E6_local_consensus"] = choice["E2_local_kernel"]
        attached = {}

        def attach(i):
            if i not in attached:
                attached[i] = self.attach(donors[i]["arbor"], soma)
            return attached[i]

        for name in ENGINE_CANDIDATES:
            i = choice[name]
            arbor, how = attach(i)
            extra = {}
            if name == "E6_local_consensus":
                others = [j for j in range(k) if j != i]
                arbor, extra = consensus_prune(arbor, [attach(j)[0] for j in others], w2[others])
            unsupported = "; ".join(f"{int(t['radius_um'])} um: {t['reason']}" for t in how["unsupported"])
            reasons = [r for r in (unsupported, extra.get("fallback_reason")) if r]
            put(name, arbor, family="engine", requested_policy=name + "+local_attachment_500",
                executed_policy=name + "+" + how["executed"] + ("+unpruned" if extra.get("fallback_reason") else ""),
                fallback=bool(reasons), fallback_reason="; ".join(reasons) or None,
                selected_donor=ids[i], donor_ids=list(ids),           # every neighbour enters the choice (E6: also the pruning)
                details={"k": k, "kernel_bandwidth_um": h, "neighbour_radius_um": float(distance.max()),
                         "nearest_neighbour_um": float(distance.min()), "same_label_neighbours": int(same.sum()),
                         "label_factor_active": bool(same.any() and not same.all()),
                         "selected_weight": float(weights.get(name, w2)[i]), "soma_shift_um": how["soma_shift_um"],
                         "attachment_radius_um": how["radius_um"], "edge_ratio_min": how.get("edge_ratio_min"),
                         "edge_ratio_max": how.get("edge_ratio_max"), **{key: v for key, v in extra.items() if key != "fallback_reason"}})
        if set(out) != set(METHODS):
            raise EngineError("a forecast is missing for a declared method")
        log = {"query": neuron, "animal": row.animal, "state_id": state_id,
               "state_sha256": regional.proof["model_state_sha256"],
               "state_animals_sha256": ids_sha(sorted(regional.proof["train_animals"])),
               "donor_animals_sha256": scope["animals_sha256"], "donor_ids_sha256": scope["ids_sha256"],
               "query_animal_in_state": False, "query_animal_in_donors": False}
        return [out[m] for m in METHODS], log

    # -- prediction run --------------------------------------------------------------------
    def queries(self):
        if self.stage == "develop":
            ids = list(self.library["ids"])           # training neurons with an included trace (truth is needed)
            if self.develop_per_animal:               # the first neurons of each animal in a fixed hash order
                by_animal = defaultdict(list)
                for n in ids:
                    by_animal[self.cohort.rows[n].animal].append(n)
                order = lambda n: sha_bytes(("SYN4RCH-axon-engine-develop:" + n).encode())
                ids = sorted(n for v in by_animal.values() for n in sorted(v, key=order)[:self.develop_per_animal])
            return ids
        return list(self.cohort.part_neurons[PART_OF_STAGE[self.stage]])   # every declared query, before any QC

    def needed_states(self):
        return sorted({self.states.state_id_for(n) for n in self.queries()})

    def prepare_prediction(self, *, verify_files=True):
        _ = self.library, self.profiles
        self.open_cache(verify_files=verify_files)
        for state_id in self.needed_states():
            self.states.get(state_id)

    def record_of(self, neuron, f, log):
        row = self.cohort.rows[neuron]
        arbor = f["arbor"]
        relative = f"predictions/{f['model']}/{neuron}.npz"
        path = self.stage_dir / relative
        G.save_arbor(path, arbor)
        uses_state = f["model"] in USES_FITTED_STATE
        return {
            "query_id": neuron, "neuron": neuron, "animal": row.animal, "full_id": row.full_id, "part": row.part,
            "known_soma_um": list(row.soma), "source_label": row.source_label, "table_tag": row.tag,
            "soma_label": row.soma_label, "hemisphere": row.hemisphere, "reconstruction_type": row.recon_type,
            "model": f["model"], "role": method_role(f["model"]), "family": f["family"],
            "requested_policy": f["requested_policy"], "executed_policy": f["executed_policy"],
            "fallback": bool(f["fallback"]), "fallback_reason": f["fallback_reason"],
            "selected_donor": f["selected_donor"], "donor_ids": list(f["donor_ids"]),
            "donor_pool_sha256": log["donor_ids_sha256"],
            "regional_state_id": log["state_id"] if uses_state else None,
            "regional_state_sha256": log["state_sha256"] if uses_state else None,
            "upstream_fit_animals_sha256": log["state_animals_sha256"] if uses_state else None,
            "details": jsonable(f["details"]),
            "path": relative, "sha256": sha_file(path), "arrays_sha256": arbor_sha(arbor), "nodes": int(len(arbor.xyz))}

    def forecast_and_save(self, neuron):
        row = self.cohort.rows[neuron]
        forecasts, log = self.forecast(neuron)
        soma = np.array(row.soma, dtype=float)
        records = []
        for f in forecasts:
            arbor = f["arbor"].validate()
            if not np.array_equal(arbor.xyz[0], soma):
                raise EngineError(f"{neuron}/{f['model']}: the forecast is not rooted at the known soma")
            if len(arbor.xyz) > 1 and np.linalg.norm(arbor.xyz[1:] - arbor.xyz[arbor.parent[1:]], axis=1).max() > STEP + 1e-6:
                raise EngineError(f"{neuron}/{f['model']}: an edge is longer than the canonical step")
            records.append(self.record_of(neuron, f, log))
        return {"query": neuron, "records": records, "log": log}

    # -- pins: everything a forecast depends on ------------------------------------------------------
    def method_pins(self):
        """Hashes of everything the forecasts of this stage depend on, recomputed from what is on disk now.
        A committed manifest carries them; they are compared again before any truth is opened and at the freeze."""
        lib = self.library
        listing = self.training_listing()
        return {"code_sha256": code_hashes(), "config_sha256": CONFIG_SHA, "subgroup_spec_sha256": SUBGROUP_SPEC_SHA,
                "protocol_sha256": self.protocol_state()["sha256"] if self.world.kind == "real" and self.stage != "develop" else None,
                "cohort": {**self.cohort.sha, "cohort_sha256": self.cohort.cohort_sha},
                "regional_table_sha256": sha_file(self.world.regional_table),
                "training_traces": {"listed": len(listing), "listing_sha256": digest(listing),
                                    "included": len(lib["ids"]), "excluded": len(lib["all"]) - len(lib["ids"])},
                "training_library_fingerprint": lib["fingerprint"], "training_ids_sha256": lib["training_ids_sha256"],
                "pair_cache_sha256": self.cache.sha256,
                "regional_states": {s: self.states.file_hashes(s) for s in self.needed_states()}}

    def copy_groups(self, queries):
        """[(donor pool, its queries)]: one pool for a held-out stage, one per held-out animal in development."""
        if self.stage != "develop":
            return [(self.scope(None), list(queries))]
        by_animal = defaultdict(list)
        for n in queries:
            by_animal[self.cohort.rows[n].animal].append(n)
        return [(self.scope(a), by_animal[a]) for a in sorted(by_animal)]

    def apply_referee_copies(self, groups, records):
        """Call the referee's copying_baselines for the given pools and make ITS arrays the official four
        controls. `records` maps (query, model) to the forecast record written by the replica; a record whose
        donor or arrays differ is overwritten with the referee's forecast and counted."""
        started, compared, differing, fallback_differing, pairs = time.time(), 0, [], 0, 0
        for scope, neurons in groups:
            result = self.referee_copies(scope, neurons)
            pairs += scope["pairs"]
            for n in neurons:
                for model in OFFICIAL_COPIES:
                    record = records[(n, model)]
                    tree, donor = result["predictions"][model][n], result["donor_proof"][model][n]
                    path = self.stage_dir / record["path"]
                    agrees = record["selected_donor"] == donor and same_tree(G.load_arbor(path), tree)
                    compared += 1
                    if not agrees:
                        differing.append([n, model, record["selected_donor"], donor])
                        arbor = G.Arbor(np.array(tree[0]), np.array(tree[1]), np.array(tree[2])).validate()
                        G.save_arbor(path, arbor)
                        record.update({"selected_donor": donor, "donor_ids": [donor], "sha256": sha_file(path),
                                       "arrays_sha256": arbor_sha(arbor), "nodes": int(len(arbor.xyz))})
                    record["official_source"] = "referee copying_baselines"
                    record["replica_agreed"] = bool(agrees)
                if bool(result["same_peptide_fallback"][n]) != bool(records[(n, "peptide_nearest")]["fallback"]):
                    fallback_differing += 1
        return {"source": "referee copying_baselines", "referee_sha256": referee_hashes(), "pools": len(groups),
                "queries": sum(len(q) for _, q in groups), "pairs_recomputed_by_the_referee": pairs, "forecasts_compared": compared,
                "replica_differences": len(differing), "replica_difference_examples": differing[:10],
                "fallback_flag_differences": fallback_differing, "seconds": round(time.time() - started, 1)}

    def predict(self, *, minutes=None):
        """Forecasts for every query x every control and candidate, the manifest, and its commitment.
        Opens no held-out trace and no held-out regional row. Resumable until the commitment is written:
        with a time budget it stops cleanly and the same command continues."""
        stage = self.stage
        self.store.require(stage, need_truth=False)
        freeze = self.verify_freeze() if stage == "test" else None
        if freeze is not None:
            self.family = freeze["family"]
        commitment = self.stage_dir / "PREDICTION_COMMITMENT.json"
        if stage != "develop" and commitment.exists():
            raise OnceOnly(f"{stage}: forecasts are already committed ({commitment}); the stage is not repeated")
        if stage != "develop" and self.world.kind == "real" and self.protocol_state()["draft"]:
            raise EngineError("the biological protocol is still marked 'STATUS: DRAFT'; it is finalised (and its hash pinned) "
                              "before any held-out forecast is committed")
        self.ensure_library()
        queries = self.queries()
        missing = [s for s in self.needed_states() if s not in self.states.manifest()["states"]]
        if missing:
            raise EngineError(f"regional states not fitted: {missing}; run `states`")
        pins = self.method_pins()
        proof = self.cache.prove()                  # sampled pairs again by the frozen pair_costs and by the referee's metric
        if not (proof["frozen_pair_costs_bit_identical"] and proof["referee_metric_bit_identical"]):
            raise EngineError("the all-pairs cache does not hold the frozen pair costs and the referee's metric on its sampled pairs")
        partial = self.stage_dir / "partial"
        if self.stage_dir.exists():
            same = (partial / "PINS.json").exists() and read_json(partial / "PINS.json") == jsonable({"pins": pins, "options": self.options()})
            if not same or commitment.exists():           # another method, or a finished development pass: start again
                if (self.stage_dir / "TRUTH_OPENED.json").exists():
                    raise OnceOnly(f"{stage}: truth was opened; forecasts are never redone")
                shutil.rmtree(self.stage_dir)
        partial.mkdir(parents=True, exist_ok=True)
        if not (partial / "PINS.json").exists():
            write_json(partial / "PINS.json", {"pins": pins, "options": self.options()})
        done = {}
        for path in sorted(partial.glob("q_*.json")):
            item = read_json(path)
            if all((self.stage_dir / r["path"]).exists() and sha_file(self.stage_dir / r["path"]) == r["sha256"] for r in item["records"]):
                done[item["query"]] = item
        todo = [n for n in queries if n not in done]
        started = time.time()
        deadline = None if minutes is None else started + 60.0 * float(minutes)
        self.say(f"{self.label()}: forecasting {len(todo)} of {len(queries)} queries x {len(METHODS)} methods "
                 f"(donor pool: the whole training library, {len(self.library['ids'])} neurons)")

        def absorb(items):
            for item in items:
                write_json(partial / f"q_{item['query']}.json", item)
                done[item["query"]] = item

        if todo:
            self.prepare_prediction()
            if self.workers > 1 and len(todo) > 2 * self.workers and self.can_parallel():
                size = max(1, min(8, len(todo) // (self.workers * 4)))
                chunks = [todo[i:i + size] for i in range(0, len(todo), size)]
                with self.pool(self.workers, _predict_init, (self.world.to_json(), stage, self.options())) as pool:
                    for count, items in enumerate(pool.imap(_predict_task, chunks), 1):
                        absorb(items)
                        if count % 20 == 0:
                            self.say(f"  {len(done)} of {len(queries)} queries forecast ({time.time() - started:.0f} s)")
                        if deadline is not None and time.time() > deadline and len(done) < len(queries):
                            pool.terminate()
                            break
            else:
                for neuron in todo:
                    absorb([self.forecast_and_save(neuron)])
                    if deadline is not None and time.time() > deadline and len(done) < len(queries):
                        break
        if len(done) < len(queries):
            self.say(f"  stopped at the time budget: {len(done)} of {len(queries)} queries forecast; run the same command to continue")
            return {"status": "incomplete", "stage": stage, "queries_forecast": len(done), "queries": len(queries),
                    "seconds": round(time.time() - started, 1)}
        predictions, logs = [], []
        for neuron in queries:
            logs.append(done[neuron]["log"])
            predictions.extend(done[neuron]["records"])
        records = {(p["query_id"], p["model"]): p for p in predictions}
        for p in predictions:
            if p["model"] in OFFICIAL_COPIES:
                p["official_source"], p["replica_agreed"] = "engine replica of the contract's rule (cache)", None
        # the four official copies through the referee's own API
        groups = self.copy_groups(queries)
        spec = CONFIG["official_copies"]
        if stage == "develop":
            chosen, used = [], 0
            for scope, neurons in groups:
                if used + scope["pairs"] > spec["develop_api_pair_budget"] and chosen:
                    break
                chosen.append((scope, neurons))
                used += scope["pairs"]
            official = self.apply_referee_copies(chosen, records)
            official.update({"api_verified": len(chosen) == len(groups), "pools_in_stage": len(groups),
                             "animals_with_referee_copies": [s["exclude_animal"] for s, _ in chosen],
                             "rule": spec["develop"]})
        elif groups[0][0]["pairs"] <= spec["heldout_inline_api_pair_limit"]:
            official = self.apply_referee_copies(groups, records)
            official["api_verified"] = True
        else:
            official = {"source": "engine replica (pinned cache, the contract's arithmetic)", "api_verified": False,
                        "pairs": groups[0][0]["pairs"], "rule": spec["above_the_limit"]}
        lib = self.library
        states = self.states.manifest()["states"]
        needed = self.needed_states()
        donor_lists = {scope["ids_sha256"]: scope["ids"] for scope, _ in groups}
        upstream = {ids_sha(sorted(self.states.get(s).proof["train_animals"])): sorted(self.states.get(s).proof["train_animals"])
                    for s in needed}
        manifest = {
            "schema": "axon_engine_prediction_manifest_v2", "label": self.label(), "world": self.world.name,
            "world_kind": self.world.kind, "stage": stage, "part": PART_OF_STAGE[stage],
            "develop_queries_per_animal": self.develop_per_animal if stage == "develop" else None,
            "family": self.family, "family_members": list(FAMILIES[self.family]),
            "pins": pins, "code_sha256": code_hashes(), "referee_sha256": referee_hashes(),
            "config": CONFIG, "config_sha256": CONFIG_SHA,
            "cohort_manifest_sha256": self.cohort.sha["manifest"], "cohort_sha256": self.cohort.cohort_sha,
            "soma_table_sha256": self.cohort.sha["soma_table"],
            "known_query_soma": CONFIG["known_query_soma"],
            "training_ids_sha256": lib["training_ids_sha256"], "training_ids": lib["ids"],
            "training_library_fingerprint": lib["fingerprint"], "training_traces": pins["training_traces"],
            "pair_cache": {"sha256": self.cache.sha256, "pairs": self.cache.pairs, "metric": pair_metric_identity(), "proof": proof},
            "donor_pool": CONFIG["donor_pools"], "official_copies": official,
            "regional_rows": self.profile_audit,
            "regional_states": [{"state_id": s, "model_state_sha256": states[s]["model_state_sha256"],
                                 "animals_sha256": states[s]["animals_sha256"], "animals": states[s]["animals"],
                                 "training_ids_sha256": states[s]["training_ids_sha256"],
                                 "files_sha256": self.states.file_hashes(s)} for s in needed],
            "upstream_fit_animals": sorted({a for s in needed for a in self.states.get(s).proof["train_animals"]}),
            "upstream_lists": upstream,
            "methods": [{"name": m, "role": method_role(m), "what": METHOD_TEXT[m], "user_label": USER_LABEL.get(m),
                         "uses_fitted_state": m in USES_FITTED_STATE,
                         "matched_controls": list(MATCHED_CONTROLS.get(m, ()))} for m in METHODS],
            "frozen_selection": None if freeze is None else {"selected": freeze["selected"], "freeze_sha256": freeze["_sha256"],
                                                             "descriptive_reference": freeze["descriptive_reference"]},
            "queries": queries, "donor_lists": donor_lists, "prediction_log": logs, "predictions": predictions}
        self.audit_manifest(manifest)
        path = write_json(self.stage_dir / "prediction_manifest.json", manifest)
        commit = {"schema": "axon_engine_prediction_commitment_v2", "stage": stage, "world": self.world.name,
                  "manifest_sha256": sha_file(path), "predictions": len(predictions), "queries": len(queries),
                  "methods": list(METHODS), "official_copies_api_verified": official["api_verified"],
                  "committed_utc": utc_now(),
                  "statement": ("development on training animals, leave one animal out: every forecast was written before it "
                                "was scored and used neither donors nor regional rows of its own animal (checked for each "
                                "forecast); the same traces serve as donors for the queries of other animals")
                  if stage == "develop" else
                  "written after every forecast was saved and before any held-out trace of this stage was opened"}
        write_json(commitment, commit, exclusive=stage != "develop", pretty=True)
        shutil.rmtree(partial)
        self.say(f"  {len(predictions)} forecasts committed: {commit['manifest_sha256']} ({time.time() - started:.0f} s)")
        return {"status": "committed", **commit}

    def audit_manifest(self, manifest):
        """Completeness and the automatic check of the prediction log (plan section 8, rule 6)."""
        queries = manifest["queries"]
        keys = [(p["query_id"], p["model"]) for p in manifest["predictions"]]
        if len(keys) != len(set(keys)) or set(keys) != {(q, m) for q in queries for m in METHODS}:
            raise EngineError("the manifest is not complete: every query x every control and candidate is required")
        if len(queries) != len(set(queries)) or set(queries) != set(self.queries()):
            raise EngineError("the manifest does not list exactly the declared queries")
        held = self.cohort.heldout_animals
        train_ids = set(manifest["training_ids"])
        if any(self.cohort.rows[n].part != "train" for n in train_ids):
            raise LeakError("a held-out neuron is in the training ID list")
        if set(manifest["upstream_fit_animals"]) & held or any(set(v) & held for v in manifest["upstream_lists"].values()):
            raise LeakError("a regional state used here was fitted on validation or test animals")
        for log in manifest["prediction_log"]:
            animals, id_animals = self.states.animal_sets(log["state_id"])
            if log["animal"] in animals or log["animal"] in id_animals or animals & held:
                raise LeakError(f"prediction log: query {log['query']} is not isolated from state {log['state_id']}")
            if ids_sha(sorted(animals)) != log["state_animals_sha256"]:
                raise EngineError("prediction log: state animal hash differs")
        pools, pool_animals = {}, {}
        for key, ids in manifest["donor_lists"].items():
            if ids_sha(ids) != key or not set(ids) <= train_ids:
                raise LeakError("a declared donor pool is not a list of training-library neurons")
            pools[key] = set(ids)
            pool_animals[key] = {self.cohort.rows[n].animal for n in ids}
        for p in manifest["predictions"]:
            pool = pools[p["donor_pool_sha256"]]
            if p["selected_donor"] not in p["donor_ids"] or not set(p["donor_ids"]) <= pool:
                raise LeakError(f"{p['query_id']}/{p['model']}: a donor is not in the declared donor pool")
            if p["animal"] in pool_animals[p["donor_pool_sha256"]] or p["animal"] != self.cohort.rows[p["query_id"]].animal:
                raise LeakError(f"{p['query_id']}/{p['model']}: the query's animal supplied a donor")
            if p["model"] in COPY_BASELINES and len(p["donor_ids"]) != 1:
                raise EngineError(f"{p['query_id']}/{p['model']}: a copy has exactly one donor")
            uses = p["model"] in USES_FITTED_STATE
            if uses != (p["regional_state_id"] is not None) or (uses and p["upstream_fit_animals_sha256"] not in manifest["upstream_lists"]):
                raise EngineError(f"{p['query_id']}/{p['model']}: its fitted-state provenance is not recorded")
        return True

    # -- the committed bundle, bound to the current method ---------------------------------------
    def verify_commitment(self):
        """The committed forecasts are complete and unchanged AND were made by the method that is on disk now:
        code, configuration, protocol, cohort tables, training traces and library, all-pairs cache and regional
        state bytes. Runs before any truth is opened, before scoring and before the freeze."""
        commitment_path = self.stage_dir / "PREDICTION_COMMITMENT.json"
        manifest_path = self.stage_dir / "prediction_manifest.json"
        if not commitment_path.exists() or not manifest_path.exists():
            raise EngineError(f"{self.stage}: there is no committed forecast bundle; run `predict` first")
        commit = read_json(commitment_path)
        if sha_file(manifest_path) != commit["manifest_sha256"]:
            raise EngineError("the prediction manifest differs from its commitment")
        manifest = read_json(manifest_path)
        self.family = manifest["family"]
        problem = self.library_problem()
        if problem:
            raise EngineError(f"the method changed since these forecasts were committed (training library): {problem}")
        try:
            now = jsonable(self.method_pins())
        except LeakError:
            raise
        except EngineError as error:
            raise EngineError(f"the method changed since these forecasts were committed ({error})")
        changed = sorted(key for key in set(now) | set(manifest["pins"]) if canonical(now.get(key)) != canonical(manifest["pins"].get(key)))
        if changed or manifest["config_sha256"] != CONFIG_SHA:
            raise EngineError("the method changed since these forecasts were committed (" + ", ".join(changed or ["config"]) +
                              "): forecasts are scored and frozen only for the method that made them")
        self.audit_manifest(manifest)
        files = {str(p.relative_to(self.stage_dir)) for p in (self.stage_dir / "predictions").rglob("*.npz")}
        if files != {p["path"] for p in manifest["predictions"]}:
            raise EngineError("forecast files on disk differ from the committed manifest")
        for p in manifest["predictions"]:
            if sha_file(self.stage_dir / p["path"]) != p["sha256"]:
                raise EngineError(f"forecast {p['path']} changed after the commitment")
        return manifest, commit

    def copies_check(self):
        """For a pool too large for the inline call: run the referee's copying_baselines on the whole stage and
        record whether the committed four copies are exactly its donors and arrays. Opens no truth."""
        manifest, commit = self.verify_commitment()
        path = self.stage_dir / "COPIES_VERIFIED.json"
        if manifest["official_copies"]["api_verified"]:
            return {"status": "already verified when the forecasts were made", **manifest["official_copies"]}
        if (self.stage_dir / "TRUTH_OPENED.json").exists():
            raise OnceOnly("truth was opened; the copies are checked before that")
        started, differing = time.time(), []
        scope = self.scope(None)
        result = self.referee_copies(scope, manifest["queries"])
        records = {(p["query_id"], p["model"]): p for p in manifest["predictions"]}
        for n in manifest["queries"]:
            for model in OFFICIAL_COPIES:
                record = records[(n, model)]
                tree = result["predictions"][model][n]
                if record["selected_donor"] != result["donor_proof"][model][n] or not same_tree(G.load_arbor(self.stage_dir / record["path"]), tree):
                    differing.append([n, model, record["selected_donor"], result["donor_proof"][model][n]])
        record = {"stage": self.stage, "manifest_sha256": commit["manifest_sha256"], "referee_sha256": referee_hashes(),
                  "queries": len(manifest["queries"]), "forecasts_compared": 4 * len(manifest["queries"]),
                  "pairs_recomputed_by_the_referee": scope["pairs"], "differences": len(differing),
                  "difference_examples": differing[:20], "all_equal": not differing,
                  "seconds": round(time.time() - started, 1), "checked_utc": utc_now()}
        write_json(path, record, pretty=True)
        return record

    def require_copies_verified(self, manifest, commit):
        if manifest["official_copies"]["api_verified"]:
            return True
        path = self.stage_dir / "COPIES_VERIFIED.json"
        record = read_json(path) if path.exists() else {}
        if record.get("manifest_sha256") != commit["manifest_sha256"] or not record.get("all_equal"):
            raise EngineError("the four official copies of this bundle have not been shown equal to the referee's own "
                              "copying_baselines; run `copies-check` (no truth is opened before it agrees)")
        return True

    # -- admission: the mechanical inclusion rules on held-out traces; reads no forecast -----------------
    def heldout_listing(self):
        part = PART_OF_STAGE[self.stage]
        out = []
        for neuron in self.cohort.part_neurons[part]:
            copy = self.store.choose(neuron)
            out.append([neuron, None, None] if copy is None else [neuron, copy["file"], copy["sha256"]])
        return out

    def admit(self):
        """Open the held-out traces of this stage for the first time and apply plan section 9 (rules 1 to 6).
        Needs the committed, method-bound forecasts but reads none of them; writes the exclusion log and the
        canonical truth with hashes. The post-admission query list is what is scored."""
        stage = self.stage
        if stage == "develop":
            raise EngineError("development scores the training library; there is nothing to admit")
        self.store.require(stage, need_truth=True)
        freeze = self.verify_freeze() if stage == "test" else None
        manifest, commit = self.verify_commitment()
        self.require_copies_verified(manifest, commit)
        if freeze is not None:
            self.require_referee_can_seal(self.referee_roles(freeze))
        done = self.stage_dir / "ADMISSION_COMPLETE.json"
        if done.exists():
            raise OnceOnly(f"{stage}: the held-out traces were already admitted ({done})")
        if (self.stage_dir / "SCORE_STARTED.json").exists() or (self.stage_dir / "referee").exists():
            raise OnceOnly(f"{stage}: a later step has started; admission is not repeated")
        listing = self.heldout_listing()
        absent = [n for n in manifest["queries"]
                  if not (self.store.choose(n) or {}).get("unavailable") and not os.path.exists((self.store.choose(n) or {}).get("path") or "")]
        if absent:
            raise StageBlocked(f"stage '{stage}': {len(absent)} held-out traces are listed in the trace manifest but their files are "
                               f"not on disk (first: {absent[0]}); the traces have not been fetched completely. Nothing was opened")
        marker = self.stage_dir / "TRUTH_OPENED.json"
        record = {"stage": stage, "world": self.world.name, "manifest_sha256": commit["manifest_sha256"],
                  "truth_listing_sha256": digest(listing), "opened_utc": utc_now(),
                  "freeze_sha256": None if freeze is None else freeze["_sha256"],
                  "purpose": "mechanical admission (plan section 9, rules 1 to 6); no forecast is read and nothing is scored"}
        resumed = False
        if marker.exists():
            # Only an admission that died before completing may be redone, on the same commitment and the same
            # trace files: it reads no forecast, scores nothing, and has one outcome.
            old = read_json(marker)
            if old.get("manifest_sha256") != record["manifest_sha256"] or old.get("truth_listing_sha256") != record["truth_listing_sha256"]:
                raise OnceOnly(f"{stage}: truth was opened against another commitment or other trace files; this bundle is not admitted")
            resumed = True
        else:
            write_json(marker, record, exclusive=True, pretty=True)
        truth_dir = self.stage_dir / "truth"
        if truth_dir.exists():
            shutil.rmtree(truth_dir)
        truth_dir.mkdir(parents=True)
        self.store.unlock_truth(commit["manifest_sha256"])
        jobs, truth_log = [], {}
        for neuron in manifest["queries"]:
            row = self.cohort.rows[neuron]
            copy = self.store.open(neuron, "truth")
            if copy.get("unavailable"):
                truth_log[neuron] = unavailable_log(neuron, row, copy)
                continue
            jobs.append((neuron, copy, (copy["path"], copy["frame"], row.hemisphere, row.recon_type, list(row.soma),
                                        copy["sha256"], str(truth_dir / f"{neuron}.npz"), alternate_copies(self.store, neuron, copy))))
        if self.workers > 1 and len(jobs) > 8 and self.can_parallel():
            with self.pool(self.workers) as pool:
                logs = pool.map(_ingest_task, [j[2] for j in jobs], chunksize=8)
        else:
            logs = [_ingest_task(j[2]) for j in jobs]
        for (neuron, copy, _), log in zip(jobs, logs):
            row = self.cohort.rows[neuron]
            log.update({"neuron": neuron, "animal": row.animal, "part": row.part, "source_file": copy["file"],
                        "manifest_sha256": copy["sha256"]})
            truth_log[neuron] = log
        truth_log = dict(sorted(truth_log.items()))
        log_path = write_json(self.stage_dir / "truth_ingestion_log.json",
                              {"part": PART_OF_STAGE[stage], "rules": RULE_TEXT, "neurons": truth_log}, pretty=True)
        admitted = [n for n in manifest["queries"] if truth_log[n]["status"] == "included"]
        excluded = {n: {"rules": l["rules_violated"], "details": l["details"]} for n, l in truth_log.items() if l["status"] != "included"}
        if not admitted:
            raise EngineError(f"{stage}: no held-out neuron passed the inclusion rules")
        summary = {"stage": stage, "world": self.world.name, "manifest_sha256": commit["manifest_sha256"],
                   "truth_listing_sha256": record["truth_listing_sha256"],
                   "declared": len(manifest["queries"]), "admitted": len(admitted), "admitted_ids": admitted,
                   "admitted_ids_sha256": ids_sha(admitted),
                   "admitted_animals": len({self.cohort.rows[n].animal for n in admitted}),
                   "declared_animals": len({self.cohort.rows[n].animal for n in manifest["queries"]}),
                   "excluded": excluded, "excluded_by_rule": dict(sorted(Counter(r for e in excluded.values() for r in e["rules"]).items())),
                   "truth": {n: {"sha256": truth_log[n]["canonical_file_sha256"], "arrays_sha256": truth_log[n]["arbor_sha256"]}
                             for n in admitted},
                   "cable_outside_atlas_box": outside_summary(truth_log[n] for n in admitted),
                   "connector_geometry": connector_summary(truth_log[n] for n in admitted),
                   "unavailable_at_source": rule8_report(truth_log, self.cohort, PART_OF_STAGE[stage]),
                   "truth_ingestion_log_sha256": sha_file(log_path), "resumed_after_interruption": resumed,
                   "forecasts_read": 0, "completed_utc": utc_now()}
        write_json(done, summary, exclusive=True, pretty=True)
        self.say(f"{self.label()}: admitted {len(admitted)} of {len(manifest['queries'])} declared queries "
                 f"({len(excluded)} excluded by the inclusion rules: {summary['excluded_by_rule'] or 'none'})")
        return summary

    def require_admission(self, commit):
        path = self.stage_dir / "ADMISSION_COMPLETE.json"
        if not path.exists():
            raise EngineError(f"{self.stage}: the held-out traces have not been admitted; run `admit` first")
        admission = read_json(path)
        if admission["manifest_sha256"] != commit["manifest_sha256"]:
            raise EngineError("the admission belongs to another commitment")
        if admission["truth_listing_sha256"] != digest(self.heldout_listing()):
            raise EngineError("the held-out trace files changed after the admission")
        on_disk = {p.stem for p in (self.stage_dir / "truth").glob("*.npz")}
        if on_disk != set(admission["admitted_ids"]):
            raise EngineError("canonical truth files differ from the admission record")
        for neuron, entry in admission["truth"].items():
            if sha_file(self.stage_dir / "truth" / f"{neuron}.npz") != entry["sha256"]:
                raise EngineError(f"canonical truth of {neuron} changed after the admission")
        admission["_sha256"] = sha_file(path)
        return admission

    # -- scoring ---------------------------------------------------------------------------
    def score(self):
        """Development: the training library as truth. Validation: this engine's scorer, once, on the admitted
        queries, then both gates. Test: the independent referee's score_once on the sealed bundle, once."""
        if self.stage == "test":
            return self.score_test()
        stage = self.stage
        self.store.require(stage, need_truth=stage != "develop")
        manifest, commit = self.verify_commitment()
        if stage == "develop":
            scored_ids = list(manifest["queries"])
            truth_of = lambda n: self.out / "library" / f"{n}.npz"
            truth_entries = [self.library["entries"][n] for n in scored_ids]
            excluded, admission = {}, None
        else:
            admission = self.require_admission(commit)
            started = self.stage_dir / "SCORE_STARTED.json"
            if started.exists() or (self.stage_dir / "metrics.json").exists() or (self.stage_dir / "SCORE_COMPLETE.json").exists():
                raise OnceOnly(f"{stage}: scoring was already started ({started}); a held-out set is scored once and a started "
                               "scoring is never repeated or resumed")
            # everything that could still refuse is tried before the single attempt is spent: the frozen gate must
            # accept its proof (training pool, admitted validation neurons, test animals, regional state)
            gate_proof = self.gate_proof([{"neuron": n} for n in admission["admitted_ids"]])
            try:
                T.select_policy([{"neuron": n, "animal": self.cohort.rows[n].animal, "full_id": self.cohort.rows[n].full_id,
                                  "metrics": {policy: {"distance_um": 1.0, "f1_200": 0.5} for policy in T.POLICIES}}
                                 for n in admission["admitted_ids"]], gate_proof)
            except ValueError as error:
                raise EngineError(f"{stage}: the frozen gate does not accept its selection proof ({error}); nothing was scored")
            write_json(started, {"stage": stage, "manifest_sha256": commit["manifest_sha256"], "admission_sha256": admission["_sha256"],
                                 "started_utc": utc_now()}, exclusive=True, pretty=True)
            scored_ids = list(admission["admitted_ids"])
            truth_of = lambda n: self.stage_dir / "truth" / f"{n}.npz"
            log = read_json(self.stage_dir / "truth_ingestion_log.json")["neurons"]
            truth_entries = [log[n] for n in scored_ids]
            excluded = admission["excluded"]
        by_query = defaultdict(dict)
        for p in manifest["predictions"]:
            by_query[p["query_id"]][p["model"]] = p
        jobs = [(n, (str(truth_of(n)), [(m, str(self.stage_dir / by_query[n][m]["path"])) for m in METHODS])) for n in scored_ids]
        self.say(f"{self.label()}: scoring {len(jobs)} of {len(manifest['queries'])} queries "
                 f"({len(manifest['queries']) - len(jobs)} excluded by the inclusion rules)")
        if self.workers > 1 and len(jobs) > 4 * self.workers and self.can_parallel():
            with self.pool(self.workers) as pool:
                scored = pool.map(_score_task, [j[1] for j in jobs], chunksize=4)
        else:
            scored = [_score_task(j[1]) for j in jobs]
        rows = []
        for (neuron, _), items in zip(jobs, scored):
            row = self.cohort.rows[neuron]
            rows.append({"neuron": neuron, "animal": row.animal, "source_label": row.source_label,
                         "soma_label": row.soma_label, "class": row.klass, "metrics": dict(items)})
        result = aggregate(rows, METHODS)
        output = {
            "schema": "axon_engine_metrics_v2", "label": self.label(), "world": self.world.name, "stage": stage,
            "scored_by": "engine (frozen metric; every score also computed by the referee's metric module, bit-identical)",
            "manifest_sha256": commit["manifest_sha256"], "code_sha256": code_hashes(), "referee_sha256": referee_hashes(),
            "config_sha256": CONFIG_SHA, "family": manifest["family"],
            "official_copies": manifest["official_copies"],
            "declared_queries": len(manifest["queries"]), "scored_queries": len(rows),
            "declared_animals": len({self.cohort.rows[n].animal for n in manifest["queries"]}),
            "scored_animals": result["animals"],
            "excluded_by_rule": dict(sorted(Counter(r for e in excluded.values() for r in e["rules"]).items())),
            "excluded": excluded, "rules": RULE_TEXT,
            "cable_outside_atlas_box": outside_summary(truth_entries),
            "unavailable_at_source": (rule8_report(self.library["all"], self.cohort, "train") if admission is None
                                      else admission["unavailable_at_source"]),
            "unit": CONFIG["unit_of_independence"], "summary": result["summary"], "per_animal": result["per_animal"],
            "fallbacks": fallback_counts(manifest["predictions"], set(scored_ids)),
            "queries": rows}
        if stage == "develop":
            output["development_only"] = True
            output["rule_shown_for_information"] = engine_gate(result["summary"], None, manifest["family"])
            output["paired_differences"] = development_differences(result)
            dev = [r for r in rows if r["class"] == "development"]
            if dev:
                sub = aggregate(dev, METHODS)
                output["development_16_subset"] = {"neurons": sub["neurons"], "animals": sub["animals"], "summary": sub["summary"]}
            output["subgroups"] = subgroup_report(rows, METHODS, None)
        else:
            frozen_gate = self.run_frozen_gate(rows, gate_proof)
            gate = engine_gate(result["summary"], frozen_gate["selected_policy"], manifest["family"])
            write_json(self.stage_dir / "frozen_gate.json", frozen_gate)
            write_json(self.stage_dir / "engine_gate.json", gate, pretty=True)
            output["admission_sha256"] = admission["_sha256"]
            output["frozen_gate"] = {"selected_policy": frozen_gate["selected_policy"], "gate_sha256": frozen_gate["gate_sha256"],
                                     "eligible_nonbaseline_policies": frozen_gate["eligible_nonbaseline_policies"],
                                     "conditional_beats_matching_uniform": frozen_gate["conditional_beats_matching_uniform"]}
            output["engine_gate"] = gate
            output["verdict"] = verdict_for(result, gate["selected"], context="validation")
            output["subgroups"] = subgroup_report(rows, METHODS, gate["selected"])
        path = write_json(self.stage_dir / "metrics.json", output, exclusive=stage != "develop")
        (self.stage_dir / "tables.md").write_text(render_tables(output))
        if stage != "develop":
            write_json(self.stage_dir / "SCORE_COMPLETE.json", {"metrics_sha256": sha_file(path), "completed_utc": utc_now(),
                                                                "manifest_sha256": commit["manifest_sha256"]}, exclusive=True, pretty=True)
        return output

    # -- the frozen prototype's own gate -------------------------------------------------------
    def gate_proof(self, validation_rows):
        pool = self.pool_donors(self.scope(None))
        state = self.states.get("all_training")
        sha = state.proof["model_state_sha256"]
        identity = lambda n: {"neuron": n, "animal": self.cohort.rows[n].animal, "full_id": self.cohort.rows[n].full_id}
        validation = [r["neuron"] for r in validation_rows]
        return {"training_animals": sorted({d["animal"] for d in pool}),
                "validation_animals": sorted({self.cohort.rows[n].animal for n in validation}),
                "test_animals": list(self.cohort.part_animals["test"]),
                "training_ids": T.donor_provenance(pool),
                "validation_ids": [identity(n) for n in validation],
                "test_ids": [identity(n) for n in self.cohort.part_neurons["test"]],
                "regional_states": [{"model": state.model, "proof": state.proof}],
                "state_by_validation_neuron": {n: sha for n in validation},
                "state_by_test_animal": {a: sha for a in self.cohort.part_animals["test"]}}

    def run_frozen_gate(self, validation_rows, proof=None):
        """The unmodified brm_template.select_policy on the five frozen policies."""
        records = [{"neuron": r["neuron"], "animal": r["animal"], "full_id": self.cohort.rows[r["neuron"]].full_id,
                    "metrics": {policy: {k: r["metrics"][METHOD_OF_FROZEN_POLICY[policy]][k] for k in SCORE_KEYS}
                                for policy in T.POLICIES}} for r in validation_rows]
        return T.select_policy(records, self.gate_proof(validation_rows) if proof is None else proof)

    def verify_frozen_gate(self, gate):
        """The checks of the frozen predict_gated, done once for the whole test run."""
        if gate.get("config") != T.CONFIG or gate.get("source_provenance") != T.source_provenance():
            raise FreezeMismatch("frozen gate: configuration or source provenance differs")
        saved = dict(gate)
        recorded = saved.pop("gate_sha256", None)
        if recorded != T._digest(saved) or Geo.canonical(T.select_policy(gate["validation_records"], gate["proof"])) != Geo.canonical(gate):
            raise FreezeMismatch("frozen gate: decision, metrics or proof changed or is not reproducible")
        proof = gate["proof"]
        pool = self.pool_donors(self.scope(None))
        if Geo.canonical(T.donor_provenance(pool)) != Geo.canonical(sorted(proof["training_ids"], key=lambda r: r["neuron"])):
            raise FreezeMismatch("frozen gate: the training pool differs from the declared training set")
        state = self.states.get("all_training")
        declared = [s for s in proof["regional_states"] if s["proof"]["model_state_sha256"] == state.proof["model_state_sha256"]]
        if len(declared) != 1 or Geo.canonical({"model": state.model, "proof": state.proof}) != Geo.canonical(declared[0]):
            raise FreezeMismatch("frozen gate: the regional state differs from the declared one")
        return True

    # -- freeze ----------------------------------------------------------------------------
    def git_state(self):
        """Engine, protocol and tests must be committed and unmodified for a real freeze and test. The pinned
        dependencies the commit does not hold (frozen prototype, referee, ignored records) are listed with hashes."""
        pinned = list(ENGINE_FILES) + ["reports/axon_engine_protocol.md", "tests/test_axon_engine.py"]
        dependencies = sorted(FROZEN_PINS["source_sha256"]) + [FREEZE_RECORD] + sorted(REFEREE_FILES.values()) + [REFEREE_CONTRACT]
        try:
            run = lambda *args: subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, timeout=30).stdout
            head = run("rev-parse", "HEAD").strip()
            dirty = run("status", "--porcelain", "--", *pinned).strip()
            tracked = set(run("ls-files", "--", *dependencies).split())
            modified = {line[3:] for line in run("status", "--porcelain", "--", *dependencies).splitlines()}
            return {"head": head or None, "engine_files_committed_and_clean": bool(head) and dirty == "",
                    "pinned_dependencies_not_reproduced_by_the_commit": {
                        d: sha_file(ROOT / d) for d in dependencies if d not in tracked or d in modified}}
        except Exception as error:
            return {"head": None, "engine_files_committed_and_clean": False, "error": type(error).__name__}

    def freeze_pins(self):
        """What the freeze pins: the method pins of the validated forecasts (with the all-training state) and the test queries."""
        pins = self.method_pins()
        return {**pins, "family": self.family,
                "test_queries_sha256": ids_sha(self.cohort.part_neurons["test"]),
                "test_queries": len(self.cohort.part_neurons["test"])}

    def freeze(self):
        """Pin code, data, parameters and the validation decision. Written once, before any test step."""
        if self.stage != "validate":
            raise EngineError("the freeze is made from the validation stage")
        validate = self.stage_dir
        path = self.out / "freeze" / "frozen.json"
        if path.exists():
            raise OnceOnly(f"the freeze record exists ({path}); it is not rewritten")
        if (self.out / "test").exists():
            raise OnceOnly("a test directory exists; the freeze must precede every test step")
        for name in ("PREDICTION_COMMITMENT.json", "TRUTH_OPENED.json", "ADMISSION_COMPLETE.json", "SCORE_STARTED.json",
                     "SCORE_COMPLETE.json", "metrics.json", "frozen_gate.json", "engine_gate.json"):
            if not (validate / name).exists():
                raise EngineError(f"freeze needs a completed validation ({name} is missing)")
        test_listed = self.store.available("test")[0]
        if self.world.kind == "real" and test_listed:
            raise EngineError(f"plan amendment 2: no test trace may be on disk before the freeze, and {test_listed} test neurons are "
                              "listed in the trace manifests; the method is not frozen in this state")
        git = self.git_state()
        if self.world.kind == "real":
            if self.protocol_state()["draft"]:
                raise EngineError("the biological protocol is still marked 'STATUS: DRAFT'; finalise it before the freeze")
            if not git["engine_files_committed_and_clean"]:
                raise EngineError("the freeze needs a committed tree: engine, protocol and tests committed and unmodified")
        manifest, commit = self.verify_commitment()        # the validated forecasts were made by the method on disk now
        metrics = read_json(validate / "metrics.json")
        metrics_sha = sha_file(validate / "metrics.json")
        if read_json(validate / "SCORE_COMPLETE.json")["metrics_sha256"] != metrics_sha or metrics["manifest_sha256"] != commit["manifest_sha256"]:
            raise EngineError("validation metrics changed after scoring or belong to another commitment")
        if metrics["code_sha256"] != jsonable(code_hashes()):
            raise EngineError("the method changed since the validation forecasts were committed (scoring code)")
        gate, frozen_gate = read_json(validate / "engine_gate.json"), read_json(validate / "frozen_gate.json")
        saved = dict(gate)
        if saved.pop("gate_sha256") != digest(saved) or gate != metrics["engine_gate"] or gate["family"] != manifest["family"]:
            raise EngineError("the engine gate record changed")
        # the all-training state is fitted again from the training rows: a stored state is not trusted on its hashes alone
        refit = self.states.refit_sha("all_training")
        stored = self.states.file_hashes("all_training")["model_state_sha256"]
        if refit != stored:
            raise EngineError("the stored all-training regional state is not what a fresh fit on the training animals gives")
        selected = gate["selected"]
        reference = gate["descriptive_reference"]
        roles = self.referee_roles({"selected": selected, "descriptive_reference": reference})
        preflight = referee_preflight()
        bundle_uses_state = sorted(m for m in self.bundle_models(roles) if m in USES_FITTED_STATE)
        record = {"schema": "axon_engine_freeze_v2", "world": self.world.name, "world_kind": self.world.kind,
                  "frozen_utc": utc_now(), "pins": self.freeze_pins(), "git": git,
                  "family": manifest["family"], "selected": selected, "descriptive_reference": reference,
                  "validation_method_equals_frozen_method": True,
                  "all_training_state_refit_sha256": refit,
                  "validation": {"manifest_sha256": commit["manifest_sha256"], "metrics_sha256": metrics_sha,
                                 "admission_sha256": sha_file(validate / "ADMISSION_COMPLETE.json"),
                                 "admitted": metrics["scored_queries"], "admitted_animals": metrics["scored_animals"],
                                 "frozen_gate_sha256": frozen_gate["gate_sha256"],
                                 "frozen_gate_file_sha256": sha_file(validate / "frozen_gate.json"),
                                 "frozen_gate_selected_policy": frozen_gate["selected_policy"],
                                 "engine_gate_sha256": gate["gate_sha256"], "verdict": gate["verdict"]},
                  "referee": {"sha256_at_freeze": referee_hashes(), "roles_at_test": roles,
                              "bundle_methods_with_a_fitted_state": bundle_uses_state, "preflight": preflight,
                              "can_seal_the_test_bundle": bool(preflight["no_upstream_fit_sealable"] and
                                                               (not bundle_uses_state or preflight["fitted_state_sealable"]))},
                  "test_rule": ("the selected candidate must be better on both animal-mean measures than every copy baseline and "
                                "every matched control of the referee declaration, on the admitted test animals; scored once by "
                                "the referee. With no selected candidate the verdict is 'no method beyond copying' whatever the "
                                "test shows"),
                  "test_access": {"test_stage_directory_exists": False, "test_traces_listed_on_disk_at_freeze": test_listed,
                                  "strongest_seal_no_test_trace_on_disk": test_listed == 0,
                                  "note": ("plan amendment 2: no test trace may be on disk before the freeze; on the real cohort the "
                                           "freeze refuses when one is listed (a synthetic world holds its generated test traces "
                                           "from the start, which this record states). The engine opens a test trace only in "
                                           "`admit --stage test`, after the test forecasts are committed")}}
        write_json(path, record, exclusive=True, pretty=True)
        return {**record, "_sha256": sha_file(path)}

    def verify_freeze(self):
        path = self.out / "freeze" / "frozen.json"
        if not path.exists():
            raise FreezeMismatch("there is no freeze record; run `freeze` after validation and before any test step")
        record = read_json(path)
        self.family = record["family"]
        problem = self.library_problem()
        if problem:
            raise FreezeMismatch(f"since the freeze these changed: training library ({problem})")
        try:
            now = jsonable(self.freeze_pins())      # validation and test queries share the all-training state and the pool
        except EngineError as error:
            raise FreezeMismatch(f"since the freeze these changed: {error}")
        changed = sorted(key for key in set(now) | set(record["pins"]) if canonical(now.get(key)) != canonical(record["pins"].get(key)))
        if changed:
            raise FreezeMismatch("since the freeze these changed: " + ", ".join(changed))
        if self.world.kind == "real":
            git = self.git_state()
            if not git["engine_files_committed_and_clean"]:
                raise FreezeMismatch("the test needs a committed tree: engine, protocol and tests must be committed and unmodified")
        validate = self.out / "validate"
        if sha_file(validate / "frozen_gate.json") != record["validation"]["frozen_gate_file_sha256"]:
            raise FreezeMismatch("the frozen gate file changed since the freeze")
        gate = read_json(validate / "engine_gate.json")
        if gate["gate_sha256"] != record["validation"]["engine_gate_sha256"] or gate["selected"] != record["selected"]:
            raise FreezeMismatch("the engine gate changed since the freeze")
        if sha_file(validate / "ADMISSION_COMPLETE.json") != record["validation"]["admission_sha256"]:
            raise FreezeMismatch("the validation admission record changed since the freeze")
        self.verify_frozen_gate(read_json(validate / "frozen_gate.json"))
        return {**record, "_sha256": sha_file(path)}

    # -- the test, through the independent referee ---------------------------------------------------
    @staticmethod
    def referee_roles(freeze):
        candidate = freeze["selected"] or freeze["descriptive_reference"]
        matched = list(ATTACHMENT_CONTROLS) + [m for m in MATCHED_CONTROLS.get(candidate, ()) if m not in ATTACHMENT_CONTROLS]
        return {"candidate": candidate, "copy_baselines": list(COPY_BASELINES), "matched_controls": matched}

    @staticmethod
    def bundle_models(roles):
        return [roles["candidate"], *roles["copy_baselines"], *roles["matched_controls"]]

    def require_referee_can_seal(self, roles):
        """Before any test truth is opened: the referee on disk must be able to seal this bundle honestly."""
        preflight = referee_preflight()
        if not preflight["no_upstream_fit_sealable"]:
            raise StageBlocked(f"the referee on disk cannot seal even a bundle without an upstream fit: {preflight['no_upstream_fit_error']}")
        fitted = sorted(m for m in self.bundle_models(roles) if m in USES_FITTED_STATE)
        if fitted and not preflight["fitted_state_sealable"]:
            raise StageBlocked(
                f"the test bundle holds methods with a fitted regional state ({', '.join(fitted)}) whose forecasts must declare the "
                "training animals behind that state, and the referee on disk refuses such a declaration after it has frozen it: "
                f"{preflight['fitted_state_error']}. Its held-out set depends on the key order of 'splits' ({REFEREE_FINDING}). "
                "The engine does not declare an empty list for a fitted state and reorders nothing. No test truth was opened. "
                "Next step: the referee's owner fixes that check; then run this command again")
        return preflight

    def build_declaration(self, base, manifest, admission, roles, freeze):
        """Everything the referee freezes before forecasts: cohort after admission, exclusions, known metadata,
        training trees with hashes, roles, geometry contract, method artifacts, evaluation settings."""
        Ref = referee()
        cohort, lib = self.cohort, self.library
        train = []
        for n in lib["ids"]:
            relative = f"training/{n}.npz"
            target = base / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(self.out / "library" / f"{n}.npz", target)
            train.append({**self.known_row(n), "arbor_path": relative, "arbor_sha256": Ref.sha(target)})
        validation_admission = read_json(self.out / "validate" / "ADMISSION_COMPLETE.json")
        validation = [self.known_row(n) for n in validation_admission["admitted_ids"]]
        test = [self.known_row(n) for n in admission["admitted_ids"]]
        included = set(lib["ids"]) | set(validation_admission["admitted_ids"]) | set(admission["admitted_ids"])
        artifacts = {}

        def bind(source, relative):
            target = base / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
            artifacts[relative] = Ref.sha(target)

        for relative in list(ENGINE_FILES) + ["tests/test_axon_engine.py"] + sorted(FROZEN_PINS["source_sha256"]):
            bind(ROOT / relative, f"method/source/{relative}")
        if self.world.protocol and os.path.exists(self.world.protocol):
            bind(self.world.protocol, "method/protocol/axon_engine_protocol.md")
        bind(self.out / "freeze" / "frozen.json", "method/freeze/frozen.json")
        for name in ("engine_gate.json", "frozen_gate.json", "ADMISSION_COMPLETE.json", "PREDICTION_COMMITMENT.json"):
            bind(self.out / "validate" / name, f"method/validation/{name}")
        for name in ("all_training.json", "all_training.proof.json"):
            bind(self.states.dir / name, f"method/regional_state/{name}")
        bind(self.out / "state_manifest.json", "method/regional_state/state_manifest.json")
        bind(self.cache.dir / "COMPLETE.json", "method/pair_cache/COMPLETE.json")
        bind(self.cache.dir / "index.json", "method/pair_cache/index.json")
        for name in ("PREDICTION_COMMITMENT.json", "prediction_manifest.json", "ADMISSION_COMPLETE.json", "truth_ingestion_log.json"):
            bind(self.stage_dir / name, f"method/test/{name}")
        fitted = sorted(m for m in self.bundle_models(roles) if m in USES_FITTED_STATE)
        state_animals = sorted(self.states.get("all_training").proof["train_animals"])
        return {
            "geometry_contract": {
                "units": "um", "axes": ["AP", "DV", "ML"],
                "frame": "Allen CCFv3; mirrored-left frame (ML' = 11400 - ML for right-hemisphere neurons, whole neuron); midline ML 5700",
                "admission": ("plan section 9 rules 1 to 6 (reports/axon_dataset_plan.md), applied mechanically by the engine to every "
                              "trace before this declaration; test queries were admitted after the forecasts were committed "
                              "(method/test/PREDICTION_COMMITMENT.json) by a step that reads no forecast; exclusions are in excluded_ids"),
                "preprocessing": ("axon = all cable that is not dendrite (types 3, 4) and not the soma; markers (10, 11, 12) inherit "
                                  "the type of their nearest non-marker ancestor; non-axon ancestors are kept as connector geometry; "
                                  "branch-to-branch resampling at 50 um by the frozen d110_geometry.resample; forks and tips kept"),
                "clipping": "none (plan amendment 1): truth, copies and candidate forecasts are complete measured trees in one coordinate domain",
                "known_soma": "paper soma table, mirrored-left frame, for training rows and queries; a forecast is rooted exactly at it",
                "truth_root": "the trace's own soma node, within 5 um of the table soma (inclusion rule 2)"},
            "method_artifacts": artifacts,
            "splits": {"train": train, "validation": validation, "test": test},
            "excluded_ids": [n for n in cohort.order if n not in included],
            "test_excluded_ids": [], "test_excluded_animals": list(cohort.never_heldout),
            "upstream_fit_animals": state_animals if fitted else [],
            "models": {"candidate": roles["candidate"], "copy_baselines": list(roles["copy_baselines"]),
                       "matched_controls": list(roles["matched_controls"])},
            "config": dict(Ref.CONFIG)}

    def seal(self):
        """Build the referee's bundle for the test and obtain its seal: freeze_declaration, copying_baselines,
        the committed forecast files copied byte for byte, seal_predictions. Opens no truth."""
        if self.stage != "test":
            raise EngineError("only the test is sealed and scored by the referee; validation is scored by this engine (protocol section 7)")
        freeze = self.verify_freeze()
        roles = self.referee_roles(freeze)
        self.require_referee_can_seal(roles)
        Ref = referee()
        started = time.time()
        base, declaration, bundle_manifest, info = self.prepare_bundle(freeze, roles)
        seal = Ref.seal_predictions(base, bundle_manifest)
        record = {"stage": "test", "world": self.world.name, **info,
                  "referee_manifest_sha256": seal["manifest_sha256"],
                  "referee_seal_sha256": sha_file(base / "PREDICTION_SEAL.json"), "referee_sha256": seal["source_sha256"],
                  "referee_functions_called": ["freeze_declaration", "copying_baselines" if info["copying_baselines_called_for_the_bundle"] else
                                               "copying_baselines (inside seal_predictions only; the committed copies were verified by `copies-check`)",
                                               "seal_predictions"],
                  "seconds": round(time.time() - started, 1), "sealed_utc": utc_now()}
        write_json(self.stage_dir / "REFEREE_SEAL.json", record, pretty=True)
        self.say(f"{self.label()}: referee seal {record['referee_seal_sha256'][:16]} on {record['forecasts']} forecasts "
                 f"({record['queries']} admitted queries x {len(record['models'])} declared methods)")
        return record

    def prepare_bundle(self, freeze, roles):
        """The referee's bundle up to the seal: its declaration frozen by freeze_declaration, the four copies
        once more from its copying_baselines (small pools), every committed forecast file of the admitted
        queries copied byte for byte, and the manifest the seal will receive."""
        manifest, commit = self.verify_commitment()
        admission = self.require_admission(commit)
        Ref = referee()
        base = self.stage_dir / "referee"
        if (base / "PREDICTION_SEAL.json").exists() or (base / "SCORE_STARTED.json").exists():
            raise OnceOnly(f"test: the referee's bundle is already sealed ({base}); it is not rebuilt")
        attempts_path = self.stage_dir / "SEAL_ATTEMPTS.json"
        attempts = read_json(attempts_path) if attempts_path.exists() else []
        if base.exists():                 # an unsealed, unscored bundle is rebuilt from the same committed inputs
            shutil.rmtree(base)
        attempts.append({"started_utc": utc_now(), "manifest_sha256": commit["manifest_sha256"], "admission_sha256": admission["_sha256"]})
        write_json(attempts_path, attempts, pretty=True)
        base.mkdir(parents=True)
        declaration = self.build_declaration(base, manifest, admission, roles, freeze)
        frozen = Ref.freeze_declaration(base, declaration)
        queries, models = admission["admitted_ids"], self.bundle_models(roles)
        records = {(p["query_id"], p["model"]): p for p in manifest["predictions"]}
        fitted_animals = sorted(self.states.get("all_training").proof["train_animals"])
        # the four official copies once more through the referee's own function (small pools), on the admitted queries
        scope = self.scope(None)
        api_in_seal = scope["pairs"] <= CONFIG["official_copies"]["heldout_inline_api_pair_limit"]
        if api_in_seal:
            donors = [{**row, "arbor": Ref.G.load_tree(base / row["arbor_path"])} for row in declaration["splits"]["train"]]
            copies = Ref.copying_baselines(donors, declaration["splits"]["test"])
        rows = []
        for q in queries:
            for model in models:
                record = records[(q, model)]
                relative = f"predictions/{model}/{q}.npz"
                target = base / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(self.stage_dir / record["path"], target)
                if Ref.sha(target) != record["sha256"]:
                    raise EngineError(f"{relative}: the copied forecast is not the committed one")
                if api_in_seal and model in OFFICIAL_COPIES:
                    tree = Ref.G.load_tree(target)
                    if record["donor_ids"] != [copies["donor_proof"][model][q]] or not all(
                            np.array_equal(a, b) for a, b in zip(tree, copies["predictions"][model][q])):
                        raise EngineError(f"{q}/{model}: the committed copy is not the referee's own copy; nothing is sealed")
                rows.append({"query_id": q, "model": model, "path": relative, "sha256": record["sha256"],
                             "donor_ids": list(record["donor_ids"]),
                             "upstream_fit_animals": fitted_animals if model in USES_FITTED_STATE else []})
        bundle_manifest = {"declaration_sha256": Ref.sha(base / "declaration.json"), "predictions": rows,
                           "upstream_fit_animals": list(declaration["upstream_fit_animals"])}
        info = {"engine_commitment_manifest_sha256": commit["manifest_sha256"], "admission_sha256": admission["_sha256"],
                "freeze_sha256": freeze["_sha256"], "roles": roles, "models": models, "queries": len(queries),
                "forecasts": len(rows), "declaration_sha256": frozen["declaration_sha256"],
                "excluded_ids": len(declaration["excluded_ids"]), "copying_baselines_called_for_the_bundle": bool(api_in_seal)}
        return base, declaration, bundle_manifest, info

    def score_test(self):
        """The referee's score_once on the sealed bundle, then the verdict and the subgroup tables read from
        ITS result. Run once: a started referee scoring is never repeated."""
        freeze = self.verify_freeze()
        manifest, commit = self.verify_commitment()
        admission = self.require_admission(commit)
        base = self.stage_dir / "referee"
        seal_path = self.stage_dir / "REFEREE_SEAL.json"
        if not (base / "PREDICTION_SEAL.json").exists() or not seal_path.exists():
            raise EngineError("test: the referee's bundle is not sealed; run `seal` first")
        if (base / "SCORE_STARTED.json").exists() or (self.stage_dir / "SCORE_COMPLETE.json").exists():
            raise OnceOnly("test: the referee's scoring was already started; the test is scored once and never repeated")
        Ref = referee()
        sealed = read_json(seal_path)
        roles, models = sealed["roles"], sealed["models"]
        if (sealed["engine_commitment_manifest_sha256"] != commit["manifest_sha256"] or sealed["admission_sha256"] != admission["_sha256"]
                or sealed["freeze_sha256"] != freeze["_sha256"] or sealed["referee_seal_sha256"] != sha_file(base / "PREDICTION_SEAL.json")
                or roles != self.referee_roles(freeze)):
            raise EngineError("test: the sealed bundle does not belong to this commitment, admission and freeze")
        # dry checks that would otherwise spend the referee's single attempt: read-only calls of its own validators
        declaration, _ = Ref.frozen_declaration(base)
        Ref.validate_bundle(base, read_json(base / "prediction_manifest.json"), declaration)
        if [q["id"] for q in declaration["splits"]["test"]] != admission["admitted_ids"]:
            raise EngineError("test: the declared test queries are not the admitted ones")
        truth_dir = self.stage_dir / "truth"
        opened = []

        def loader(query_id):
            opened.append(query_id)
            return truth_dir / f"{query_id}.npz"

        self.say(f"{self.label()}: the referee scores {len(admission['admitted_ids'])} admitted queries x {len(models)} methods")
        result = Ref.score_once(base, loader)
        rows = []
        for q in result["queries"]:
            row = self.cohort.rows[q["id"]]
            if (q["truth_sha256"] != Ref.sha(base / q["truth_path"])
                    or arbor_sha(G.load_arbor(base / q["truth_path"])) != admission["truth"][q["id"]]["arrays_sha256"]):
                raise EngineError("the truth the referee stored is not the admitted truth")
            rows.append({"neuron": q["id"], "animal": row.animal, "source_label": row.source_label, "soma_label": row.soma_label,
                         "class": row.klass, "metrics": {m: {k: float(q["models"][m][k]) for k in KEPT_METRICS} for m in models}})
        mine = aggregate(rows, models)
        for m in models:
            if any(mine["summary"][m][k] != result["overall"]["summary"][m][k] for k in SCORE_KEYS):
                raise EngineError("the engine's animal means differ from the referee's")
        selected = freeze["selected"]
        verdict = verdict_from_referee(result, roles, selected)
        # the committed forecasts of the other candidates, scored for the record only, after the referee has finished
        others = [m for m in CANDIDATES if m not in models]
        by_query = defaultdict(dict)
        for p in manifest["predictions"]:
            by_query[p["query_id"]][p["model"]] = p
        other_rows = []
        for q in result["queries"]:
            truth = G.load_arbor(base / q["truth_path"])
            other_rows.append({"animal": q["animal"], "metrics": {
                m: score_arbor_checked(G.load_arbor(self.stage_dir / by_query[q["id"]][m]["path"]), truth) for m in others}})
        other = aggregate(other_rows, others) if others else {"summary": {}}
        output = {
            "schema": "axon_engine_metrics_v2", "label": self.label(), "world": self.world.name, "stage": "test",
            "scored_by": "independent referee (analysis/prediction_independent_evaluator.py score_once)",
            "referee": {"sha256": result["source_sha256"], "sha256_at_freeze": freeze["referee"]["sha256_at_freeze"],
                        "changed_since_the_freeze": result["source_sha256"] != freeze["referee"]["sha256_at_freeze"],
                        "declaration_sha256": result["declaration_sha256"],
                        "manifest_sha256": result["manifest_sha256"], "seal_sha256": result["seal_sha256"],
                        "metrics_file_sha256": sha_file(base / "metrics.json"), "status": result["status"], "roles": roles,
                        "beats_all_copying": bool(result["overall"]["beats_all_copying"]),
                        "beats_all_matched": bool(result["overall"]["beats_all_matched"]), "limits": result["limits"]},
            "manifest_sha256": commit["manifest_sha256"], "code_sha256": code_hashes(), "config_sha256": CONFIG_SHA,
            "family": freeze["family"], "frozen_selection": {"selected": selected, "descriptive_reference": freeze["descriptive_reference"],
                                                             "freeze_sha256": freeze["_sha256"]},
            "official_copies": manifest["official_copies"],
            "declared_queries": len(manifest["queries"]), "scored_queries": len(rows),
            "declared_animals": admission["declared_animals"], "scored_animals": mine["animals"],
            "excluded_by_rule": admission["excluded_by_rule"], "excluded": admission["excluded"], "rules": RULE_TEXT,
            "cable_outside_atlas_box": admission["cable_outside_atlas_box"],
            "unavailable_at_source": admission["unavailable_at_source"],
            "unit": CONFIG["unit_of_independence"], "methods_scored_by_the_referee": models,
            "summary": mine["summary"], "per_animal": mine["per_animal"],
            "fallbacks": fallback_counts([p for p in manifest["predictions"] if p["model"] in models], set(admission["admitted_ids"])),
            "verdict": verdict,
            "subgroups": subgroup_report(rows, models, selected, controls=[m for m in models if m in CONTROLS]),
            "referee_subgroups": {"by_peptide": {k: {x: v[x] for x in ("neurons", "independent_animals", "interpretation")}
                                                 for k, v in result["by_peptide"].items()},
                                  "by_nucleus": {k: {x: v[x] for x in ("neurons", "independent_animals", "interpretation")}
                                                 for k, v in result["by_nucleus"].items()}},
            "other_candidates_descriptive": {
                "note": ("committed forecasts of the candidates that were not selected, scored by this engine after the referee "
                         "had finished, against the truth the referee stored; outside the referee's bundle, outside the rule, "
                         "no claim and no selection"),
                "summary": {m: {k: other["summary"][m][k] for k in SCORE_KEYS} for m in others}},
            "queries": rows}
        path = write_json(self.stage_dir / "metrics.json", output, exclusive=True)
        (self.stage_dir / "tables.md").write_text(render_tables(output))
        write_json(self.stage_dir / "SCORE_COMPLETE.json", {"metrics_sha256": sha_file(path), "completed_utc": utc_now(),
                                                            "manifest_sha256": commit["manifest_sha256"],
                                                            "referee_metrics_sha256": sha_file(base / "metrics.json"),
                                                            "truth_files_opened_by_the_referee": len(opened)}, exclusive=True, pretty=True)
        return output

    # -- status ----------------------------------------------------------------------------
    def status(self):
        out = {"world": self.world.name, "kind": self.world.kind, "cohort": self.cohort.summary(),
               "traces_on_disk": {p: dict(zip(("have", "of"), self.store.available(p))) for p in ("train", "validation", "test")},
               "stages": {}, "protocol_draft": self.protocol_state()["draft"]}
        for stage in STAGES:
            try:
                self.store.require(stage, need_truth=stage != "develop")
                out["stages"][stage] = "data present"
            except StageBlocked as blocked:
                out["stages"][stage] = "BLOCKED: " + str(blocked)
        problem = self.library_problem()
        out["training_library"] = problem or "current"
        if not problem:
            out["pair_cache"] = self.pair_cache().progress()
        out["referee"] = {"sha256": referee_hashes(), "preflight": referee_preflight()}
        return out


def fallback_counts(predictions, queries):
    """How often each method executed something else than it was asked for, on the scored queries."""
    out = {}
    for p in predictions:
        if p["query_id"] in queries:
            entry = out.setdefault(p["model"], {"forecasts": 0, "fallbacks": 0, "reasons": Counter()})
            entry["forecasts"] += 1
            if p["fallback"]:
                entry["fallbacks"] += 1
                entry["reasons"][p["fallback_reason"]] += 1
    return {m: {"forecasts": e["forecasts"], "fallbacks": e["fallbacks"], "reasons": dict(sorted(e["reasons"].items()))}
            for m, e in out.items()}


def development_differences(result):
    """Paired animal differences of every candidate against the strongest control on each measure, with
    bootstrap intervals. Development information only."""
    summary, per_animal = result["summary"], result["per_animal"]
    strongest = strongest_controls(summary)
    out = {"strongest_control": strongest, "candidates": {}}
    for c in CANDIDATES:
        out["candidates"][c] = {
            "versus_strongest_on_distance": {"control": strongest["distance_um"],
                                             **compare(per_animal, c, strongest["distance_um"], f"develop|{c}|d")},
            "versus_strongest_on_f1": {"control": strongest["f1_200"], **compare(per_animal, c, strongest["f1_200"], f"develop|{c}|f")}}
    return out


# ----------------------------------------------------------------------------------------------
# tables
# ----------------------------------------------------------------------------------------------


def _short(model):
    return model.split("_")[0] if model[:1] in "CEF" and model[1:2].isdigit() else model


def render_tables(output):
    """Markdown tables of one scored stage (written next to metrics.json)."""
    lines = [f"# {output['label']}", ""]
    if output.get("development_only"):
        lines += ["**Development data, not a test. No success claim is drawn from this table.**", ""]
    lost = output.get("unavailable_at_source") or {"neurons": 0, "scored_cohort_smaller_than_locked_notice": False}
    if lost["scored_cohort_smaller_than_locked_notice"]:      # plan amendment 2: said in the first lines
        lines += [f"**The scored cohort is smaller than the locked one: {lost['neurons']} of {lost['of']} neurons of this part are "
                  f"unavailable at the source (rule 8); brain samples lost entirely: "
                  f"{', '.join(lost['brain_samples_lost_entirely']) or 'none'}.**", ""]
    lines += [f"Scored by: {output.get('scored_by', 'engine')}.", "",
              f"Scored: {output['scored_queries']} of {output['declared_queries']} declared queries, "
              f"{output['scored_animals']} of {output['declared_animals']} animals (brain sample IDs). "
              f"Excluded by the mechanical inclusion rules: {sum(len(e['rules']) > 0 for e in output['excluded'].values())} "
              f"(by rule: {output['excluded_by_rule'] or 'none'})."
              + (f" Unavailable at the source (rule 8): {lost['neurons']} (by source label: {lost['by_source_label']})."
                 if lost["neurons"] else ""), ""]
    summary = output["summary"]
    shown = [m for m in METHODS if m in summary]
    lines += ["| method | role | animal-mean distance, um | animal-mean F1 at 200 um | fallbacks |", "|---|---|---:|---:|---|"]
    fallbacks = output.get("fallbacks", {})
    for m in shown:
        fb = fallbacks.get(m)
        note = "" if not fb else f"{fb['fallbacks']} of {fb['forecasts']}"
        role = method_role(m).replace("_", " ") + (f" ({USER_LABEL[m]})" if m in USER_LABEL else "")
        lines.append(f"| {m} | {role} | {summary[m]['distance_um']:.1f} | {summary[m]['f1_200']:.4f} | {note} |")
    lines.append("")
    other = output.get("other_candidates_descriptive")
    if other and other["summary"]:
        lines += [f"Candidates that were not selected ({other['note']}):", "",
                  "| method | animal-mean distance, um | animal-mean F1 |", "|---|---:|---:|"]
        lines += [f"| {m} | {v['distance_um']:.1f} | {v['f1_200']:.4f} |" for m, v in other["summary"].items()]
        lines.append("")
    if "development_16_subset" in output:
        sub = output["development_16_subset"]
        lines += [f"The {sub['neurons']} earlier development neurons only ({sub['animals']} animals), same forecasts:", "",
                  "| method | animal-mean distance, um | animal-mean F1 |", "|---|---:|---:|"]
        lines += [f"| {m} | {sub['summary'][m]['distance_um']:.1f} | {sub['summary'][m]['f1_200']:.4f} |" for m in shown]
        lines.append("")
    paired = output.get("paired_differences")
    if paired:
        s = paired["strongest_control"]
        lines += [f"Paired animal differences, candidate minus the strongest control on each measure (distance: {s['distance_um']}; "
                  f"F1: {s['f1_200']}), with paired 95% bootstrap intervals over animals. Descriptive.", "",
                  "| candidate | distance difference, um | interval | animals better | F1 difference | interval | animals better |",
                  "|---|---:|---|---:|---:|---|---:|"]
        for c, v in paired["candidates"].items():
            d, f = v["versus_strongest_on_distance"], v["versus_strongest_on_f1"]
            di, fi = d["distance_bootstrap95_um"] or [float("nan")] * 2, f["f1_bootstrap95"] or [float("nan")] * 2
            lines.append(f"| {c} | {d['distance_delta_um']:+.1f} | {di[0]:+.1f} to {di[1]:+.1f} | {d['animals_better_on_distance']} of {d['animals']} "
                         f"| {f['f1_delta']:+.4f} | {fi[0]:+.4f} to {fi[1]:+.4f} | {f['animals_better_on_f1']} of {f['animals']} |")
        lines.append("")
    gate = output.get("engine_gate") or output.get("rule_shown_for_information")
    if gate:
        title = "Rule shown for information (development: no selection is made)" if output.get("development_only") else "Validation gate"
        lines += [f"## {title}", "", f"Family: {gate['family']}. Strongest control on distance: {gate['strongest_control']['distance_um']}; "
                  f"on F1: {gate['strongest_control']['f1_200']}. Selected: {gate['selected'] or 'none'} ({gate['verdict']}).", "",
                  "| candidate | beats every control | controls not beaten | beats matched controls | admitted | eligible |", "|---|---|---|---|---|---|"]
        for c, v in gate["candidates"].items():
            lost = ", ".join(_short(k) for k, w in v["beats_control"].items() if not w) or "-"
            lines.append(f"| {c} | {'yes' if v['beats_every_control'] else 'no'} | {lost} | "
                         f"{'yes' if v['beats_every_matched_control'] else 'no'} | {'yes' if v['admitted'] else 'no: ' + str(v['note'])} | "
                         f"{'yes' if v['eligible'] else 'no'} |")
        lines.append("")
    verdict = output.get("verdict")
    if verdict:
        lines += [f"**Verdict ({verdict['context']}): {verdict['verdict']}.** {verdict['statement']}", ""]
    groups = output.get("subgroups")
    if groups:
        selected = groups["selected_candidate"]
        columns = [m for m in ([selected] if selected else []) + [m for m in CONTROLS if m in summary] if m]
        extra = [m for m in shown if m not in columns and m in output.get("methods_scored_by_the_referee", [])]
        columns += extra
        head = "| " + " | ".join(f"{_short(m)} dist | {_short(m)} F1" for m in columns) + " |"
        rule = "|" + "---:|" * (2 * len(columns))
        cells = lambda methods: " | ".join(f"{methods[m]['distance_um']:.0f} | {methods[m]['f1_200']:.3f}" for m in columns)
        lines += ["## Per animal (brain sample ID)", "", "| animal | source label | neurons " + head, "|---|---|---:" + rule]
        for a in groups["per_animal"]:
            lines.append(f"| {a['animal']} | {a['source_label']} | {a['neurons']} | {cells(a['methods'])} |")
        for title, key, field in (("By source label", "by_source_label", "source_label"), ("By soma label", "by_soma_label", "soma_label")):
            lines += ["", f"## {title}", "", f"| {field} | neurons | animals | note " + head, "|---|---:|---:|---" + rule]
            for g in groups[key]:
                note = "sparse: descriptive only" if g["sparse_descriptive_only"] else "descriptive"
                if key == "by_soma_label" and not g["nucleus_assigned"]:
                    note += "; " + g["soma_label_kind"].replace("_", " ")
                lines.append(f"| {g[field]} | {g['neurons']} | {g['animals']} | {note} | {cells(g['methods'])} |")
        lines += ["", f"Every subgroup row is descriptive and carries no population claim; rows with fewer than "
                      f"{groups['sparse_below_animals']} animals are marked sparse."]
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------------------------
# checks on the training traces that are on disk (run by `selfcheck`)
# ----------------------------------------------------------------------------------------------


def check_frozen_equivalence(engine):
    """On every training trace on disk that codes its axon as type 2: this engine's extraction,
    restricted to the frozen extractor's node set, gives the frozen canonical arbor bit for bit; and
    what the plan's axon rule (terminal marker nodes are axon cable) changes. Nothing is clipped."""
    rows = []
    for neuron in engine.cohort.part_neurons["train"]:
        for copy in engine.store.copies.get(neuron, []):
            if copy.get("unavailable"):
                continue
            engine.store.open(neuron, "donor")
            row = engine.cohort.rows[neuron]
            reflect = bool(copy["true_hemisphere"] and row.hemisphere == "Right")
            raw = read_swc(copy["path"])
            if 2 not in set(raw[:, 1].astype(int).tolist()):
                continue
            frozen = G.resample(G.extract_explicit_axon(copy["path"], reflect), STEP)
            replay = frozen_replay_arbor(copy["path"], reflect)
            log, plan = ingest_trace(copy["path"], frame=copy["frame"], hemisphere=row.hemisphere,
                                     recon_type=row.recon_type, table_soma=row.soma, expected_sha=copy["sha256"])
            entry = {"neuron": neuron, "file": copy["file"], "bit_identical_to_frozen": bool(same_arbor(frozen, replay)),
                     "frozen_nodes": int(len(frozen.xyz)), "frozen_cable_um": frozen.cable(), "plan_status": log["status"]}
            if plan is not None:
                m = G.geometry_metrics(plan, frozen)
                entry.update({"plan_nodes": int(len(plan.xyz)), "plan_cable_um": plan.cable(),
                              "terminal_marker_twigs": log["terminal_marker_twig_nodes"],
                              "terminal_marker_twig_um": log["terminal_marker_twig_um"],
                              "outside_atlas_box_um": log["outside_atlas_box"]["outside_um"],
                              "plan_equals_frozen": bool(same_arbor(plan, frozen)),
                              "plan_vs_frozen_distance_um": m["distance_um"], "plan_vs_frozen_f1_200": m["f1_200"]})
            rows.append(entry)
    return {"files": len(rows), "all_bit_identical": bool(rows) and all(r["bit_identical_to_frozen"] for r in rows),
            "plan_rule_changes_arbor_in": sum(1 for r in rows if not r.get("plan_equals_frozen", True)),
            "max_terminal_marker_twig_um": max((r.get("terminal_marker_twig_um", 0.0) for r in rows), default=0.0),
            "max_outside_atlas_box_um": max((r.get("outside_atlas_box_um", 0.0) for r in rows), default=0.0),
            "max_plan_vs_frozen_distance_um": max((r.get("plan_vs_frozen_distance_um", 0.0) for r in rows), default=0.0),
            "rows": rows}


def check_two_copies(engine):
    """Neurons with a true-hemisphere copy and a mirrored-left copy: after the frame rule the two files
    hold the same numbers, and both give the same canonical arbor."""
    rows = []
    for neuron, copies in sorted(engine.store.copies.items()):
        true_copies = [c for c in copies if c.get("true_hemisphere") and not c.get("unavailable")]
        mirrored = [c for c in copies if not c.get("true_hemisphere") and not c.get("unavailable")]
        if not true_copies or not mirrored:
            continue
        row = engine.cohort.rows[neuron]
        engine.store.open(neuron, "donor")
        a, b = read_swc(true_copies[0]["path"]), read_swc(mirrored[0]["path"])
        if row.hemisphere == "Right":
            a = a.copy()
            a[:, 4] = MIRROR_ML - a[:, 4]
        arbors = []
        for copy in (true_copies[0], mirrored[0]):
            _, arbor = ingest_trace(copy["path"], frame=copy["frame"], hemisphere=row.hemisphere,
                                    recon_type=row.recon_type, table_soma=row.soma, expected_sha=copy["sha256"])
            arbors.append(arbor)
        rows.append({"neuron": neuron, "hemisphere": row.hemisphere, "nodes": int(len(a)),
                     "parsed_tables_equal_after_frame_rule": bool(a.shape == b.shape and np.array_equal(a, b)),
                     "canonical_arbors_equal": bool(arbors[0] is not None and arbors[1] is not None and same_arbor(*arbors))})
    return {"neurons": len(rows), "all_equal": bool(rows) and all(r["parsed_tables_equal_after_frame_rule"] and r["canonical_arbors_equal"] for r in rows),
            "rows": rows}


def replay_earlier_medoid(root=ROOT, engine=None):
    """Recompute the earlier representative-copy figure on the 16 development neurons with (a) the
    frozen extractor and (b) this engine's extraction in frozen mode: same donors (the other
    development neurons), same anchoring (trace soma). The earlier run reported 922.49 um / 0.1978.
    With an engine, also (c) the same 16-neuron experiment under the plan's axon rule (terminal markers)
    with the forecast anchored at the table soma, to show what those two choices change. History only:
    these figures are never a threshold."""
    root = Path(root)
    declared = read_json(root / "data/interim/brm_geometry/declaration.json")["queries"]
    saved = read_json(root / "data/interim/brm_geometry/metrics.json")["summary"]["medoid"]
    out = {"earlier": {"distance_um": saved["distance_um"], "f1_200": saved["f1_200"]}}
    variants = [("frozen_extractor", lambda q: G.resample(G.extract_explicit_axon(root / q["file"], q["reflect_to_left"]), STEP),
                 lambda q: np.array(q["swc_soma_um"], dtype=float)),
                ("engine_extraction_frozen_mode", lambda q: frozen_replay_arbor(root / q["file"], q["reflect_to_left"]),
                 lambda q: np.array(q["swc_soma_um"], dtype=float))]
    if engine is not None:
        variants.append(("plan_rules_and_table_soma", lambda q: G.load_arbor(engine.out / "library" / f"{q['neuron']}.npz"),
                         lambda q: np.array(engine.cohort.rows[q["neuron"]].soma, dtype=float)))
    for name, loader, anchor in variants:
        records = sorted(({"neuron": q["neuron"], "animal": q["animal"], "arbor": loader(q), "soma": anchor(q)}
                          for q in declared), key=lambda r: r["neuron"])
        rows = []
        for q in records:
            donors = [d for d in records if d["animal"] != q["animal"]]
            prediction = donors[G.training_medoid(donors)]["arbor"].shifted(q["soma"])
            rows.append({"animal": q["animal"], "metrics": {"medoid": score_arbor(prediction, q["arbor"])}})
        summary = aggregate(rows, ["medoid"])["summary"]["medoid"]
        out[name] = {"distance_um": summary["distance_um"], "f1_200": summary["f1_200"],
                     "equals_earlier_exactly": bool(summary["distance_um"] == saved["distance_um"] and summary["f1_200"] == saved["f1_200"])}
    return out


def prove_memo(engine, queries=3):
    """The pinned cache and the memoised medoid give the frozen functions' numbers bit for bit on the WHOLE
    donor pool of a query, and the unmodified predict_policy returns identical trees and records with and
    without the memo. (Without the memo the frozen functions recompute every pair: use small pools.)"""
    lib = engine.library
    chosen = lib["ids"][:: max(1, len(lib["ids"]) // queries)][:queries]
    rows = []
    for neuron in chosen:
        row = engine.cohort.rows[neuron]
        scope = engine.scope(row.animal if row.part == "train" else None)
        pool = engine.pool_donors(scope)
        regional = engine.states.get(engine.states.state_id_for(neuron))
        soma = np.array(row.soma, dtype=float)
        query = T.Query(row.neuron, row.animal, row.full_id, soma.copy(), soma.copy(), row.tag)
        D0, F0 = T.pair_costs(pool)
        medoid0 = G.training_medoid(pool)
        with engine.memo.frozen_memo():
            D1, F1 = T.pair_costs(pool)
            medoid1 = G.training_medoid(pool)
        entry = {"query": neuron, "donors": len(pool), "pair_distance_bit_identical": bool(np.array_equal(D0, D1)),
                 "pair_f1_loss_bit_identical": bool(np.array_equal(F0, F1)), "medoid_identical": medoid0 == medoid1,
                 "policies": {}}
        for policy in T.POLICIES:
            plain, info0 = T.predict_policy(query, pool, regional, policy)
            with engine.memo.frozen_memo():
                memo, info1 = T.predict_policy(query, pool, regional, policy)
            entry["policies"][policy] = bool(same_arbor(plain, memo) and Geo.canonical(jsonable(info0)) == Geo.canonical(jsonable(info1)))
        rows.append(entry)
    ok = all(r["pair_distance_bit_identical"] and r["pair_f1_loss_bit_identical"] and r["medoid_identical"]
             and all(r["policies"].values()) for r in rows)
    return {"queries": len(rows), "all_bit_identical": bool(rows) and ok, "rows": rows}


def check_official_copies(engine, max_pairs=None):
    """Leave one animal out on the training library: the engine's replica of the four official copies against
    the referee's own copying_baselines, donors and arrays, for every query of every pool within the pair budget."""
    budget = CONFIG["official_copies"]["develop_api_pair_budget"] if max_pairs is None else max_pairs
    lib = engine.library
    animals = sorted(set(lib["animal"].tolist()))
    used, pools, compared, differing = 0, 0, 0, []
    started = time.time()
    for animal in animals:
        scope = engine.scope(animal)
        if pools and used + scope["pairs"] > budget:
            break
        neurons = [n for n in lib["ids"] if engine.cohort.rows[n].animal == animal]
        result = engine.referee_copies(scope, neurons)
        pool = engine.pool_donors(scope)
        for n in neurons:
            row = engine.cohort.rows[n]
            replica = engine.replica_copies(row, scope)
            for model in OFFICIAL_COPIES:
                donor = scope["ids"][replica[model]]
                arbor = pool[replica[model]]["arbor"].shifted(np.array(row.soma, dtype=float))
                compared += 1
                if donor != result["donor_proof"][model][n] or not same_tree(arbor, result["predictions"][model][n]):
                    differing.append([n, model, donor, result["donor_proof"][model][n]])
            if bool(result["same_peptide_fallback"][n]) != bool(replica["fallback"]):
                differing.append([n, "same_label_fallback_flag", replica["fallback"], result["same_peptide_fallback"][n]])
        used += scope["pairs"]
        pools += 1
    return {"pools_checked": pools, "pools": len(animals), "forecasts_compared": compared, "differences": len(differing),
            "difference_examples": differing[:10], "all_equal": compared > 0 and not differing,
            "pairs_recomputed_by_the_referee": used, "seconds": round(time.time() - started, 1)}


def cost_probe(engine, sizes=(500, 1000, 2000, 4448), pair_sample=120):
    """Measured cost on the training arbors that are on disk, and what it means for the whole training library.

    (a) time per pair of the frozen metric (the cache's cost) and of the referee's metric (the cost inside its
    copying_baselines and its seal); (b) one call of the UNMODIFIED frozen predict_policy on donor pools of the
    given sizes, with memoised pair costs. For (b) the pool is the local arbors repeated under made-up IDs and the
    pair-cost matrix holds placeholder numbers: only the per-call cost of the frozen function is timed, no
    forecast is made or kept."""
    import resource
    Ref = referee()
    lib = engine.library
    arbors = [engine.donor(n)["arbor"] for n in lib["ids"]]
    centred = [a.shifted(np.zeros(3)) for a in arbors]
    trees = [Ref.G.shifted(as_tree(a), np.zeros(3)) for a in arbors]
    pairs = [(i, j) for i in range(len(arbors)) for j in range(i + 1, len(arbors))]
    step = max(1, len(pairs) // pair_sample)
    sample = pairs[::step][:pair_sample]
    t = time.time()
    frozen = [pair_cost(centred[i], centred[j]) for i, j in sample]
    frozen_ms = 1000.0 * (time.time() - t) / len(sample)
    t = time.time()
    theirs = [Ref.G.geometry_metrics(trees[i], trees[j]) for i, j in sample]
    referee_ms = 1000.0 * (time.time() - t) / len(sample)
    identical = all(a[0] == b["distance_um"] and a[1] == 1.0 - b["f1_200"] for a, b in zip(frozen, theirs))
    nodes = [len(a.xyz) for a in arbors]
    full = 4448
    full_pairs = full * (full - 1) // 2
    out = {"what": "measured on the local training arbors; projections assume the whole library has arbors of the same size",
           "local_arbors": len(arbors), "canonical_nodes_median": float(np.median(nodes)), "canonical_nodes_mean": float(np.mean(nodes)),
           "canonical_nodes_max": int(max(nodes)), "pairs_timed": len(sample),
           "ms_per_pair_frozen_metric": frozen_ms, "ms_per_pair_referee_metric": referee_ms,
           "frozen_and_referee_bit_identical_on_the_timed_pairs": bool(identical),
           "projection_4448_training_neurons": {
               "pairs": full_pairs, "cache_cpu_hours": frozen_ms * full_pairs / 3.6e6,
               "cache_wall_hours_4_workers": frozen_ms * full_pairs / 3.6e6 / 4,
               "cache_matrices_bytes": 2 * full * full * 8 + full * full,
               "referee_copying_baselines_hours_single_process_per_call": referee_ms * full_pairs / 3.6e6,
               "referee_calls": "one inside seal_predictions at the test, one per `copies-check` (validation and test)"}}
    row = engine.cohort.rows[lib["ids"][0]]
    regional = engine.states.get(engine.states.state_id_for(row.neuron))
    soma = np.array(row.soma, dtype=float)
    query = T.Query(row.neuron, row.animal, row.full_id, soma.copy(), soma.copy(), row.tag)
    donors_real = [engine.donor(n) for n in lib["ids"] if engine.cohort.rows[n].animal != row.animal]
    rng = np.random.default_rng(0)
    calls = []
    for n in sizes:
        pool = []
        for i in range(n):
            d = donors_real[i % len(donors_real)]
            animal = f"9{i // 44:05d}"                         # about 44 donors per made-up animal, as in the cohort
            neuron = f"{animal}_{i % 44:03d}"
            pool.append(dict(d, neuron=neuron, animal=animal, full_id=neuron + "_x_y"))
        pool.sort(key=lambda d: d["neuron"])
        values = rng.uniform(300.0, 3000.0, size=(n, n))
        D = np.triu(values, 1)
        D = D + D.T
        F = np.full((n, n), 0.5)
        np.fill_diagonal(F, 0.0)
        memo = FrozenMemo(lambda: None)
        for d in pool:
            memo._arbor[d["neuron"]] = d["arbor"]
        memo._pool = (ids_sha([d["neuron"] for d in pool]), D, F)
        entry = {"donors": n}
        with memo.frozen_memo():
            t = time.time()
            T.predict_policy(query, pool, regional, "baseline")
            entry["first_call_with_the_medoid_s"] = time.time() - t
            for policy in T.POLICIES:
                t = time.time()
                _, info = T.predict_policy(query, pool, regional, policy)
                entry[policy + "_s"] = time.time() - t
                del info
        entry["five_policies_s"] = sum(entry[p + "_s"] for p in T.POLICIES)
        entry["peak_rss_mb_so_far"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1e6 if sys.platform == "darwin" else 1e3)
        calls.append(entry)
        del pool, D, F, memo, values
    out["frozen_predict_policy_per_query"] = calls
    last = calls[-1]
    if last["donors"] == full:
        per_query = last["five_policies_s"]
        out["projection_4448_training_neurons"].update({
            "frozen_five_policies_seconds_per_query": per_query,
            "frozen_validation_1185_queries_cpu_hours": per_query * 1185 / 3600,
            "frozen_test_1547_queries_cpu_hours": per_query * 1547 / 3600,
            "note": "`predict` is resumable (`--minutes`), so no single run has to exceed a quarter of an hour"})
    return out



# ----------------------------------------------------------------------------------------------
# synthetic worlds (mechanism checks only: synthetic success is not biological success)
# ----------------------------------------------------------------------------------------------

SYNTHETIC_KINDS = {
    "informative": "source label and soma position carry route information: the dominant route of a neuron is f(label, soma zone), drawn with probability 0.75, any other route otherwise",
    "none": "neither label nor soma position carries route information (route drawn uniformly at random per neuron)",
    "dense": "the soma neighbourhood carries route information, the source label does not: dominant route f(soma zone) with probability 0.75; labels assigned at random",
    "informative_deterministic": "route = f(label, soma zone) with probability 1: the regime in which copying the nearest same-label neuron is already nearly optimal",
}
SYNTHETIC_SEEDS = {"informative": 2026111001, "none": 2026111002, "dense": 2026111003, "informative_deterministic": 2026111004}
SYNTHETIC_DOMINANCE = 0.75
SYNTHETIC_LABELS = ("SynA", "SynB", "SynC", "SynD")
SYNTHETIC_SIZES = {
    "standard": {"train": (8, 7, 5, 4), "validation": (3, 2, 2, 1), "test": (5, 4, 2, 1), "neurons": (8, 14)},
    "tiny": {"train": (3, 3, 2, 2), "validation": (1, 1, 1, 1), "test": (1, 1, 1, 1), "neurons": (4, 5)},
}
SYNTHETIC_PLANTED = {"train": ("soma_off", "short_axon", "both_codes"),
                     "validation": ("broken_tree", "wrong_hemisphere", "no_dendrite"),
                     "test": ("soma_off", "both_codes", "no_dendrite", "short_axon", "wrong_hemisphere")}
_PRIMARY = np.array([[3600.0, 3800.0, 4600.0], [12950.0, 5600.0, 4900.0], [6400.0, 2600.0, 4300.0], [6800.0, 6200.0, 7900.0]])
_SECONDARY = np.array([[5200.0, 5000.0, 3600.0], [8200.0, 4800.0, 3800.0], [5600.0, 7000.0, 4200.0], [7800.0, 7100.0, 5200.0],
                       [4800.0, 6200.0, 5000.0], [8800.0, 6600.0, 4400.0], [6200.0, 4200.0, 5400.0], [7400.0, 3800.0, 3400.0]])


def _unit(v):
    return v / np.linalg.norm(v)


def _polyline(a, b, rng, *, wobble=True, step=20.0):
    length = float(np.linalg.norm(b - a))
    n = max(2, int(math.ceil(length / step)) + 1)
    t = np.linspace(0.0, 1.0, n)
    points = a + t[:, None] * (b - a)
    if wobble and length > 300.0:
        u = _unit(np.cross((b - a) / length, rng.normal(size=3)))
        amplitude, waves, phase = rng.uniform(30.0, 80.0), length / rng.uniform(600.0, 1200.0), rng.uniform(0.0, 2 * math.pi)
        points = points + (amplitude * np.sin(math.pi * t) * np.sin(2 * math.pi * waves * t + phase))[:, None] * u
    return points


def _synthetic_route(kind, label_index, soma, rng):
    """Route family 0..7 (primary target = family // 2, secondary target = family)."""
    zone_ml, zone_ap = int(soma[2] < 5000.0), int(soma[0] >= 6500.0)
    if kind == "none":
        family = int(rng.integers(8))
    else:
        dominant = 2 * label_index + zone_ml if kind.startswith("informative") else 2 * (2 * zone_ap + zone_ml)
        if kind == "informative_deterministic" or rng.random() < SYNTHETIC_DOMINANCE:
            family = dominant
        else:
            family = [f for f in range(8) if f != dominant][int(rng.integers(7))]
    return family // 2, family


def _synthetic_neuron(kind, label_index, soma, rng, *, axon_code, dendrites=True, planted=None):
    """One synthetic neuron in the mirrored-left frame: lists of xyz, parent index and node type."""
    xyz, parent, types = [soma.copy()], [-1], [1]

    def path(points, start, kind_code):
        previous, indices = start, []
        for p in points[1:]:
            xyz.append(p)
            parent.append(previous)
            types.append(kind_code)
            previous = len(xyz) - 1
            indices.append(previous)
        return indices

    primary, secondary = _synthetic_route(kind, label_index, soma, rng)
    junction = soma + 600.0 * _unit(_PRIMARY[primary] - soma) + rng.normal(0.0, 40.0, 3)
    trunk = path(_polyline(soma, junction, rng, wobble=False), 0, axon_code)
    if planted == "short_axon":
        del xyz[21:], parent[21:], types[21:]                 # about 400 um of axon
    else:
        main_end = _PRIMARY[primary] + rng.normal(0.0, 120.0, 3)
        main = path(_polyline(junction, main_end, rng), trunk[-1], axon_code)
        for _ in range(3):
            path(_polyline(main_end, main_end + rng.uniform(250.0, 500.0) * _unit(rng.normal(size=3)), rng, wobble=False), main[-1], axon_code)
        second_end = _SECONDARY[secondary] + rng.normal(0.0, 120.0, 3)
        second = path(_polyline(junction, second_end, rng), trunk[-1], axon_code)
        for _ in range(2):
            path(_polyline(second_end, second_end + rng.uniform(250.0, 500.0) * _unit(rng.normal(size=3)), rng, wobble=False), second[-1], axon_code)
        if rng.random() < 0.5:                                 # an idiosyncratic collateral
            start = main[int(len(main) * rng.uniform(0.3, 0.8))]
            path(_polyline(xyz[start], xyz[start] + rng.uniform(1000.0, 2000.0) * _unit(rng.normal(size=3)), rng), start, axon_code)
        if rng.random() < 0.3:                                 # a marker inside the path (connector)
            types[main[len(main) // 2]] = 10
        if rng.random() < 0.3:                                 # a marker at a tip (terminal twig)
            xyz.append(xyz[second[-1]] + np.array([3.0, 2.0, 1.0]))
            parent.append(second[-1])
            types.append(12)
    if planted == "both_codes":
        other = 0 if axon_code == 2 else 2
        for i in range(len(types) // 2, len(types)):
            if types[i] == axon_code:
                types[i] = other
    if dendrites and planted != "no_dendrite":
        for _ in range(3):
            path(_polyline(soma, soma + 150.0 * _unit(rng.normal(size=3)), rng, wobble=False), 0, 3)
    return np.array(xyz), np.array(parent), np.array(types), {"primary": primary, "secondary": secondary}


def make_synthetic_world(directory, kind, *, size="standard", seed=None):
    """Write a synthetic world in the real layout. The generator is declared here and uses nothing of
    the engine's candidate code. Returns the path of the world's data directory."""
    if kind not in SYNTHETIC_KINDS:
        raise EngineError(f"unknown synthetic world '{kind}'")
    seed = SYNTHETIC_SEEDS[kind] if seed is None else int(seed)
    spec = SYNTHETIC_SIZES[size]
    directory = Path(directory)
    data = directory / "data"
    if data.exists():
        raise EngineError(f"{data} exists; a synthetic world is not overwritten")
    (data / "traces").mkdir(parents=True)
    rng = np.random.default_rng(seed)
    centroids = np.random.default_rng(20261110).uniform([2500.0, 1500.0, 2500.0], [13200.0, 7800.0, 8500.0], size=(len(REGIONS), 3))
    animals = []
    for part in ("train", "validation", "test"):
        for label_index, count in enumerate(spec[part]):
            animals.extend((part, label_index) for _ in range(count))
    numbers = 700000 + rng.permutation(len(animals)) * 7 + 11
    manifest, soma_rows, regional_rows, trace_rows, truth = [], [], [], [], {}
    planted_counter = Counter()
    for animal_index, ((part, label_index), number) in enumerate(zip(animals, numbers)):
        animal, label = str(int(number)), SYNTHETIC_LABELS[label_index]
        axon_code = 2 if animal_index % 2 == 0 else 0
        mirrored_file = animal_index % 5 == 0
        for k in range(int(rng.integers(spec["neurons"][0], spec["neurons"][1] + 1))):
            neuron = f"{animal}_{k + 1:03d}"
            soma = rng.uniform([6200.0, 5700.0, 4700.0], [6800.0, 6300.0, 5300.0])
            planted = None
            if k == 1 and planted_counter[part] < len(SYNTHETIC_PLANTED[part]):
                planted = SYNTHETIC_PLANTED[part][planted_counter[part]]
                planted_counter[part] += 1
            hemisphere = "Right" if (rng.random() < 0.5 or planted == "wrong_hemisphere") else "Left"
            axon_only = rng.random() < 0.1 and planted is None
            xyz, parent, types, route = _synthetic_neuron(kind, label_index, soma, rng, axon_code=axon_code,
                                                         dendrites=not axon_only, planted=planted)
            xyz = np.round(xyz, 2)
            # regional lengths from the generated axon (mirrored frame): cable per nearest centroid
            lengths = np.zeros(len(REGIONS))
            axon_nodes = np.flatnonzero(np.isin(types, (0, 2, 10, 12)))
            seg = np.linalg.norm(xyz[axon_nodes] - xyz[parent[axon_nodes]], axis=1)
            nearest = np.argmin(((xyz[axon_nodes][:, None, :] - centroids[None]) ** 2).sum(axis=2), axis=1)
            np.add.at(lengths, nearest, seg)
            zone = 2 * int(soma[0] >= 6500.0) + int(soma[2] < 5000.0)
            draw = rng.random()
            soma_label = "HY" if draw < 0.15 else ("fx" if draw < 0.19 else ("SNa", "SNb", "SNc", "SNd")[zone])
            offset = _unit(rng.normal(size=3)) * (rng.uniform(0.5, 1.9) if rng.random() < 0.15 else 0.002)
            table_soma = xyz[0] + offset + (np.array([30.0, 0.0, 0.0]) if planted == "soma_off" else 0.0)
            table_soma[2] = min(table_soma[2], MIDLINE_ML)
            full_id = f"{neuron}_{label}_{soma_label}"
            written = xyz.copy()
            if hemisphere == "Right" and not mirrored_file:
                written[:, 2] = MIRROR_ML - written[:, 2]
            ids = np.arange(1, len(xyz) + 1)
            parent_ids = np.where(parent < 0, -1, parent + 1)
            if planted == "broken_tree":
                parent_ids[len(parent_ids) // 2] = 99999999
            relative = f"{part}/{animal}/{neuron}.swc"
            path = data / "traces" / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("w") as handle:
                handle.write("# synthetic trace: id type AP DV ML radius parent\n")
                for i in range(len(xyz)):
                    handle.write(f"{ids[i]} {types[i]} {written[i, 0]:.2f} {written[i, 1]:.2f} {written[i, 2]:.2f} 1.0 {parent_ids[i]}\n")
            manifest.append([neuron, animal, label, soma_label, "synthetic", "0", "0", "0", part, f"synthetic {kind} world",
                             "Left" if planted == "wrong_hemisphere" else hemisphere,
                             "Axon_only" if axon_only else "Axon_and_dendrite", soma_label])
            soma_rows.append([full_id, neuron, label, soma_label] + [f"{v:.6f}" for v in table_soma])
            regional_rows.append([full_id] + [f"{v:.2f}" for v in lengths])
            trace_rows.append([relative, neuron, str(path.stat().st_size), sha_file(path), "synthetic generator",
                               "synthetic_mirrored_left" if mirrored_file else "synthetic_true_hemisphere"])
            truth[neuron] = {"part": part, "planted": planted, **route}
    order = np.argsort([m[0] for m in manifest])
    head = ["neuron", "animal", "source_label", "soma_label", "class", "animal_inspected", "animal_development",
            "animal_in_li_2024", "part", "part_reason", "hemisphere", "reconstruction_type", "index_soma_acronym"]
    with (data / "manifest.tsv").open("w") as handle:
        handle.write("\t".join(head) + "\n")
        for i in order:
            handle.write("\t".join(manifest[i]) + "\n")
    with (data / "soma.csv").open("w") as handle:
        handle.write("file,neurons,geno,location,AP,DV,ML\n")
        for i in order:
            handle.write(",".join(soma_rows[i]) + "\n")
    with (data / "regional.csv").open("w") as handle:
        handle.write("," + ",".join(REGIONS) + "\n")
        for i in order:
            handle.write(",".join(regional_rows[i]) + "\n")
    with (data / "traces/MANIFEST.tsv").open("w") as handle:
        handle.write("file\tneuron\tbytes\tsha256\tsource\tframe\n")
        for i in order:
            handle.write("\t".join(trace_rows[i]) + "\n")
    write_json(data / "world.json", {"name": f"synthetic_{kind}", "kind": kind, "what": SYNTHETIC_KINDS[kind], "seed": seed,
                                     "size": size, "sizes": spec, "labels": list(SYNTHETIC_LABELS),
                                     "planted_inclusion_failures": SYNTHETIC_PLANTED,
                                     "note": "mechanism check only; synthetic success is not biological success"}, pretty=True)
    write_json(directory / "generator_truth.json", truth, pretty=True)     # never read by the engine
    return data


def scramble_outcomes(directory, *, parts=("validation", "test"), animals=None, seed=7):
    """Canary helper: replace the traces and regional rows of the given parts (or animals) by other
    values, leaving every permitted input (manifest, soma table) untouched. Trees stay valid."""
    data = Path(directory) / "data"
    rng = np.random.default_rng(seed)
    with (data / "manifest.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    chosen = {r["neuron"] for r in rows if (r["animal"] in animals if animals is not None else r["part"] in parts)}
    with (data / "traces/MANIFEST.tsv").open(newline="") as handle:
        traces = list(csv.DictReader(handle, delimiter="\t"))
    for t in traces:
        if t["neuron"] not in chosen:
            continue
        path = data / "traces" / t["file"]
        lines = [l for l in path.read_text().splitlines() if l and not l.startswith("#")]
        out = ["# scrambled outcome"]
        amplitude, phase = rng.uniform(300.0, 900.0, 3), rng.uniform(0.0, 6.28, 3)
        for i, line in enumerate(lines):
            f = line.split()
            if i > 0 and f[1] != "3":
                for k in range(3):
                    f[2 + k] = f"{float(f[2 + k]) + amplitude[k] * math.sin(phase[k] + i / 37.0) * min(1.0, i / 40.0):.2f}"
            out.append(" ".join(f))
        path.write_text("\n".join(out) + "\n")
        t["bytes"], t["sha256"] = str(path.stat().st_size), sha_file(path)
    with (data / "traces/MANIFEST.tsv").open("w") as handle:
        handle.write("file\tneuron\tbytes\tsha256\tsource\tframe\n")
        for t in traces:
            handle.write("\t".join(t[k] for k in ("file", "neuron", "bytes", "sha256", "source", "frame")) + "\n")
    lines = (data / "regional.csv").read_text().splitlines()
    for i in range(1, len(lines)):
        full_id = lines[i].split(",", 1)[0]
        if "_".join(full_id.split("_", 2)[:2]) in chosen:
            lines[i] = full_id + "," + ",".join(f"{v:.2f}" for v in rng.uniform(0.0, 5000.0, len(REGIONS)))
    (data / "regional.csv").write_text("\n".join(lines) + "\n")
    return sorted(chosen)


def hidden_outcome_canary(directory, stage, *, workers=1, log=None):
    """Plan section 7 canary. Copy the world, replace every held-out trace and regional row by other
    values, rebuild library, pair cache and states from scratch there, forecast again, and compare every
    forecast with the committed one array for array. (The all-pairs cache is a function of the training arbors
    only: entries of unchanged arbors are imported from the original cache, the rest recomputed.)"""
    directory = Path(directory)
    copy = directory.parent / (directory.name + "__canary")
    if copy.exists():
        shutil.rmtree(copy)
    shutil.copytree(directory / "data", copy / "data")
    changed = scramble_outcomes(copy)
    original = Path(synthetic_world(directory).out)
    manifest = read_json(original / stage / "prediction_manifest.json")
    engine = Engine(synthetic_world(copy), stage, workers=workers, log=log, family=manifest["family"])
    engine.build_library()
    engine.build_cache(siblings=[original / "paircache"])
    engine.states.build_many(engine.needed_states())
    engine.prepare_prediction()
    committed = {(p["query_id"], p["model"]): p for p in manifest["predictions"]}
    differing, skipped = [], 0
    for neuron in engine.queries():
        forecasts, _ = engine.forecast(neuron)
        for f in forecasts:
            record = committed[(neuron, f["model"])]
            if record.get("replica_agreed") is False:       # the committed copy is the referee's; the replica differed at a tie
                skipped += 1
            elif arbor_sha(f["arbor"]) != record["arrays_sha256"]:
                differing.append((neuron, f["model"]))
    shutil.rmtree(copy)
    return {"stage": stage, "held_out_neurons_scrambled": len(changed), "forecasts_compared": len(committed) - skipped,
            "official_copies_skipped_because_the_replica_differed": skipped,
            "forecasts_differing": len(differing), "bit_identical": not differing, "examples": differing[:5]}


def run_synthetic(kind, directory, *, size="standard", workers=1, log=print, develop=True, canary=True, family="all"):
    """The whole pipeline on one synthetic world: library, all-pairs cache, states, development pass,
    validation (forecasts, admission, scoring, both gates), freeze, test forecasts, canaries, admission, the
    referee's seal and single scoring, and the refusal to repeat any once-only step."""
    directory = Path(directory)
    started = time.time()
    make_synthetic_world(directory, kind, size=size)
    world = synthetic_world(directory)
    say = log or (lambda text: None)
    brief = lambda summary: {m: {k: summary[m][k] for k in SCORE_KEYS} for m in summary}
    out = {"world": world.name, "what": SYNTHETIC_KINDS[kind], "size": size, "family": family,
           "note": "mechanism check only; synthetic success is not biological success"}
    dev = Engine(world, "develop", workers=workers, log=log)
    out["library"] = dev.build_library()
    built = dev.build_cache()
    out["pair_cache"] = {"pairs": built["pairs"], "seconds": built["seconds"], "sha256": dev.cache.sha256, "proof": dev.cache.prove()}
    states = dev.states.build_many(list(dev.states.plan) if develop else ["all_training"], workers=workers)
    out["states"] = {s: {"animals": e["animals"], "selected_k": e["fit"]["selected_k"],
                         "held_out_animals_in_state": e["validation_or_test_animals_in_state"]} for s, e in states.items()}
    if develop:
        commit = dev.predict()
        d = dev.score()
        out["develop"] = {"queries": d["scored_queries"], "animals": d["scored_animals"], "summary": brief(d["summary"]),
                          "official_copies": {k: d["official_copies"][k] for k in ("api_verified", "pools", "pools_in_stage",
                                                                                   "forecasts_compared", "replica_differences")}}
    val = Engine(world, "validate", workers=workers, log=log, family=family)
    val.predict()
    out["validation_canary"] = hidden_outcome_canary(directory, "validate", log=None) if canary else None
    admitted = val.admit()
    v = val.score()
    out["validate"] = {"declared": v["declared_queries"], "scored": v["scored_queries"], "animals": v["scored_animals"],
                       "excluded_by_rule": v["excluded_by_rule"], "summary": brief(v["summary"]),
                       "official_copies": {k: v["official_copies"].get(k) for k in ("source", "api_verified", "forecasts_compared", "replica_differences")},
                       "frozen_gate_selected_policy": v["frozen_gate"]["selected_policy"], "family": v["engine_gate"]["family"],
                       "eligible": v["engine_gate"]["eligible"], "selected": v["engine_gate"]["selected"],
                       "descriptive_reference": v["engine_gate"]["descriptive_reference"],
                       "verdict": v["engine_gate"]["verdict"], "strongest_control": v["engine_gate"]["strongest_control"],
                       "admission_read_forecasts": admitted["forecasts_read"]}
    frozen = val.freeze()
    out["freeze"] = {"selected": frozen["selected"], "descriptive_reference": frozen["descriptive_reference"],
                     "family": frozen["family"], "sha256": frozen["_sha256"], "referee_roles_at_test": frozen["referee"]["roles_at_test"],
                     "referee_can_seal_the_test_bundle": frozen["referee"]["can_seal_the_test_bundle"],
                     "bundle_methods_with_a_fitted_state": frozen["referee"]["bundle_methods_with_a_fitted_state"]}
    test = Engine(world, "test", workers=workers, log=log)
    commit = test.predict()
    out["test_commitment"] = commit["manifest_sha256"]
    out["test_canary"] = hidden_outcome_canary(directory, "test", log=None) if canary else None
    try:
        admission = test.admit()
        sealed = test.seal()
        t = test.score()
        out["test"] = {"declared": t["declared_queries"], "admitted": admission["admitted"], "scored": t["scored_queries"],
                       "animals": t["scored_animals"], "admission_read_forecasts": admission["forecasts_read"],
                       "excluded_by_rule": t["excluded_by_rule"], "scored_by": t["scored_by"],
                       "referee_functions_called": sealed["referee_functions_called"] + ["score_once"],
                       "referee_seal_sha256": sealed["referee_seal_sha256"], "referee_roles": sealed["roles"],
                       "excluded_ids_declared_before_the_seal": sealed["excluded_ids"],
                       "summary": brief(t["summary"]), "other_candidates_descriptive": t["other_candidates_descriptive"]["summary"],
                       "selected": t["verdict"]["selected_candidate"], "success": t["verdict"]["success"],
                       "verdict": t["verdict"]["verdict"], "statement": t["verdict"]["statement"],
                       "strongest_control": t["verdict"]["strongest_control"],
                       "referee_beats_all_copying": t["referee"]["beats_all_copying"],
                       "referee_beats_all_matched": t["referee"]["beats_all_matched"],
                       "source_labels": [{k: g[k] for k in ("source_label", "neurons", "animals", "sparse_descriptive_only")}
                                         for g in t["subgroups"]["by_source_label"]],
                       "soma_labels": [{k: g[k] for k in ("soma_label", "neurons", "animals", "sparse_descriptive_only", "soma_label_kind")}
                                       for g in t["subgroups"]["by_soma_label"]]}
    except StageBlocked as blocked:
        out["test"] = {"blocked": str(blocked), "test_truth_opened": (test.stage_dir / "TRUTH_OPENED.json").exists(),
                       "verdict": "not scored: blocked before any test truth was opened", "selected": frozen["selected"]}
    refusals = {}
    steps = [("second_test_predict", lambda: Engine(world, "test", workers=1, log=None).predict()),
             ("second_test_admit", lambda: Engine(world, "test", workers=1, log=None).admit()),
             ("second_test_seal", lambda: Engine(world, "test", workers=1, log=None).seal()),
             ("second_test_score", lambda: Engine(world, "test", workers=1, log=None).score()),
             ("second_validation_predict", lambda: Engine(world, "validate", workers=1, log=None).predict()),
             ("second_validation_admit", lambda: Engine(world, "validate", workers=1, log=None).admit()),
             ("second_validation_score", lambda: Engine(world, "validate", workers=1, log=None).score()),
             ("second_freeze", lambda: Engine(world, "validate", workers=1, log=None).freeze())]
    for name, action in steps:
        try:
            action()
            refusals[name] = "NOT REFUSED"
        except OnceOnly as refused:
            refusals[name] = "refused: " + str(refused)[:90]
        except StageBlocked as blocked:
            refusals[name] = "blocked: " + str(blocked)[:90]
        except EngineError as missing:               # a step whose prerequisite never ran because the stage is blocked
            refusals[name] = "not possible: " + str(missing)[:90]
    out["once_only"] = refusals
    out["seconds"] = round(time.time() - started, 1)
    write_json(directory / "summary.json", out, pretty=True)
    say(f"synthetic '{kind}' (family {family}): validation selected {out['validate']['selected']} ({out['validate']['verdict']}); "
        f"test: {out['test']['verdict']} ({out['seconds']} s)")
    return out



# ----------------------------------------------------------------------------------------------
# development diagnostics that need no trace (training regional rows; soma tables)
# ----------------------------------------------------------------------------------------------


def regional_proxy(engine, k_grid=(8, 16, 32)):
    """DEVELOPMENT diagnostic at regional resolution, TRAINING animals only, leave one animal out.

    Not geometry and not a test. For each training neuron its 49-region log profile is 'predicted' by the
    profile analogue of the four official copying controls (global medoid, nearest soma, same-label nearest,
    same-label medoid) and of the engine's local medoids (K nearest somas of other animals; uniform, kernel, kernel
    times label factor), with Euclidean distance between profiles in place of the geometric distance.
    It asks, on the dense table that is on disk, what soma position and source label carry.
    """
    from scipy.spatial.distance import cdist
    cohort, spec = engine.cohort, CONFIG["engine"]
    rows = [cohort.rows[n] for n in cohort.part_neurons["train"]]
    if any(r.part != "train" for r in rows):
        raise LeakError("regional proxy: a held-out neuron was requested")
    Y = np.array([engine.profiles[r.full_id] for r in rows])
    X = np.array([r.soma for r in rows], dtype=float)
    label, animal, tags = np.array([r.source_label for r in rows]), np.array([r.animal for r in rows]), [r.tag for r in rows]
    index = {r.neuron: i for i, r in enumerate(rows)}
    P = cdist(Y, Y)
    methods = ["global_medoid (medoid, C1 analogue)", "nearest_soma (nearest analogue)", "label_nearest (peptide_nearest analogue)",
               "label_medoid (peptide_medoid analogue)"]
    for k in k_grid:
        methods += [f"local_uniform_K{k} (E1 analogue)", f"local_kernel_K{k} (E2 analogue)", f"local_kernel_label_K{k} (E3 analogue)"]
    methods += ["joint_model_out_of_fold (frozen regional state; a mean, not a copy)", "local_kernel_mean_K16 (a mean, not a copy)"]
    pred = {m: np.zeros_like(Y) for m in methods}
    for state_id, plan in engine.states.plan.items():
        if not state_id.startswith("oof_"):
            continue
        state = engine.states.get(state_id)
        served = np.array([index[n] for n in plan["served"]])
        posterior = Geo.state_posterior(state.model, X[served], [tags[i] for i in served])
        pred[methods[-2]][served] = Geo.predicted_profile(state.model, posterior)
    neighbour_stats = []
    for a in sorted(set(animal.tolist())):
        q, d = np.flatnonzero(animal == a), np.flatnonzero(animal != a)
        count = Counter(animal[d].tolist())
        w = np.zeros(len(rows))
        w[d] = [1.0 / count[x] for x in animal[d]]
        w /= w.sum()
        risk = P @ w
        medoid = d[int(np.argmin(risk[d]))]
        soma_d = cdist(X[q], X[d])
        same = d[label[d] == label[q[0]]]
        if len(same):
            ws = np.zeros(len(rows))
            cs = Counter(animal[same].tolist())
            ws[same] = [1.0 / cs[x] for x in animal[same]]
            label_medoid = same[int(np.argmin((P @ (ws / ws.sum()))[same]))]
        else:
            label_medoid = medoid
        same_mask = label[d] == label[q[0]]
        for row, i in enumerate(q):
            dist = soma_d[row]
            nearest = int(np.argmin(dist))
            pred[methods[0]][i], pred[methods[1]][i] = Y[medoid], Y[d[nearest]]
            pred[methods[2]][i] = Y[d[np.flatnonzero(same_mask)[int(np.argmin(dist[same_mask]))]]] if same_mask.any() else Y[d[nearest]]
            pred[methods[3]][i] = Y[label_medoid]
            order = np.lexsort((np.arange(len(d)), dist))
            for k in k_grid:
                nb = np.sort(order[:k])
                g, dk = d[nb], dist[nb]
                ck = Counter(animal[g].tolist())
                prior = np.array([1.0 / ck[x] for x in animal[g]]) / len(ck)
                h = max(spec["kernel_h0_um"], float(np.sort(dk)[min(spec["kernel_rank"], k) - 1]))
                kernel = np.exp(-dk ** 2 / (2.0 * h * h))
                lab = np.where(label[g] == label[i], 1.0, spec["label_other_factor"])
                sub = P[np.ix_(g, g)]
                for name, wk in ((f"local_uniform_K{k} (E1 analogue)", prior), (f"local_kernel_K{k} (E2 analogue)", prior * kernel),
                                 (f"local_kernel_label_K{k} (E3 analogue)", prior * kernel * lab)):
                    pred[name][i] = Y[g[int(np.argmin(sub @ (wk / wk.sum())))]]
                if k == 16:
                    wk = prior * kernel
                    pred[methods[-1]][i] = (wk / wk.sum()) @ Y[g]
                    neighbour_stats.append((float(dk.min()), float(dk.max()), h, int((label[g] == label[i]).sum())))
    animals = sorted(set(animal.tolist()))
    table = {}
    for m in methods:
        mse = ((pred[m] - Y) ** 2).mean(axis=1)
        r = np.array([np.corrcoef(pred[m][i], Y[i])[0, 1] if np.ptp(pred[m][i]) > 0 and np.ptp(Y[i]) > 0 else 0.0 for i in range(len(Y))])
        table[m] = {"animal_mean_mse_log10": float(np.mean([mse[animal == a].mean() for a in animals])),
                    "animal_mean_profile_r": float(np.mean([r[animal == a].mean() for a in animals]))}
    stats = np.array(neighbour_stats)
    return {"what": "DEVELOPMENT diagnostic, regional resolution, training animals only, leave one animal out; not geometry, not a test",
            "neurons": len(rows), "animals": len(animals), "regional_rows": engine.profile_audit,
            "settings": {"k_grid": list(k_grid), **{k: spec[k] for k in ("kernel_h0_um", "kernel_rank", "label_other_factor")}},
            "neighbourhood_K16": {"nearest_um_median": float(np.median(stats[:, 0])), "radius_um_median": float(np.median(stats[:, 1])),
                                  "bandwidth_um_median": float(np.median(stats[:, 2])),
                                  "same_label_neighbours_median": float(np.median(stats[:, 3])),
                                  "same_label_neighbours_mean": float(np.mean(stats[:, 3]))},
            "methods": table}



def coverage_diagnostic(engine):
    """Inputs only (soma tables, labels) of TRAINING and VALIDATION animals: how the dense training library
    covers a query, and how often the frozen prototype's local attachment can apply at all. Training queries
    are taken leave one animal out. No test animal's soma or label is used here (requirement D: test animals
    stay out of design choices), no trace and no regional row is read. Assumes every training trace passes the
    inclusion rules."""
    from scipy.spatial.distance import cdist
    cohort, spec = engine.cohort, CONFIG["engine"]
    train = [cohort.rows[n] for n in cohort.part_neurons["train"]]
    T_xyz = np.array([r.soma for r in train], dtype=float)
    T_label = np.array([r.source_label for r in train])
    T_animal = np.array([r.animal for r in train])
    k = spec["k_neighbours"]
    pct = lambda v: [float(np.percentile(v, p)) for p in (10, 50, 90)]
    out = {"what": ("inputs only: table somas and labels of training animals (leave one animal out) and validation animals; "
                    "no test animal; assumes every training trace is included"),
           "training_neurons": len(train), "training_animals": len(set(T_animal.tolist()))}

    def block(rows, D):
        ordered = np.sort(D, axis=1)
        nb = np.argsort(D, axis=1)[:, :k]
        same = np.array([(T_label[nb[i]] == rows[i].source_label).sum() for i in range(len(rows))])
        reach500, reach1000 = D <= 500.0 / 3, D <= 1000.0 / 3
        return {"queries": len(rows),
                "nearest_training_soma_um_p10_p50_p90": pct(ordered[:, 0]),
                "kernel_rank_distance_um_p10_p50_p90": pct(ordered[:, spec["kernel_rank"] - 1]),
                "neighbourhood_radius_um_p10_p50_p90": pct(ordered[:, k - 1]),
                "same_label_among_k_nearest_p10_p50_p90": pct(same), "queries_with_no_same_label_neighbour": int((same == 0).sum()),
                "share_of_queries_whose_nearest_donor_is_within_167um": float((ordered[:, 0] <= 500.0 / 3).mean()),
                "share_of_donor_query_pairs_within_167um": float(reach500.mean()),
                "share_of_donor_query_pairs_within_333um": float(reach1000.mean()),
                # one template for every query (C1, C4, C5): the best any single training neuron could reach
                "best_single_training_donor_share_of_queries_within_167um": float(reach500.mean(axis=0).max()),
                "best_single_training_donor_share_of_queries_within_333um": float(reach1000.mean(axis=0).max())}

    D = cdist(T_xyz, T_xyz)
    D[T_animal[:, None] == T_animal[None, :]] = np.inf           # leave one animal out
    out["training_leave_one_animal_out"] = block(train, D)
    rows = [cohort.rows[n] for n in cohort.part_neurons["validation"]]
    if rows:
        out["validation"] = block(rows, cdist(np.array([r.soma for r in rows], dtype=float), T_xyz))
    return out
