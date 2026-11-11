"""SYNARCH axon bridge: the general reconstruction layer as a supplier of PERMITTED QUERY INPUTS for the axon
engine (analysis/axon_engine.py). Library; the engine imports it lazily, it never imports the engine.

What it does, for one query neuron of any part of the cohort:

  1. it makes the MEASURED records of the permitted inputs only: table soma position, source label, and (where
     the layer's registry has them) fine soma label, hemisphere and reconstruction type, each with the registered
     scale text, value tag and frame (core.scale_statement, checked by core.statement_problem). The neuron's
     measured regional profile is never made into a record here, not even for a training neuron in development;
  2. it asks the layer for the neuron's regional axon-length profile from soma position and source label (the
     layer's declared task T1), through Reconstructor.reconstruct under the layer's DEFAULT policy. The
     reconstructor comes ONLY from synarch_recon_benchmark.build_reconstructor over a finished, verified run:
     fold=None (the all-training serving states) for a neuron outside the training partition, fold=k for a
     training neuron of the layer's outer fold k, so that the producing state never saw the neuron's brain sample;
  3. it turns the answer into the layer's own inferred record (PropertyResult.to_observation) and passes the
     measured records and that record through the axon adapter OF THE SAME reconstructor (AxonArborAdapter.plan,
     which calls check_inputs, then query_for). What comes back is an AxonQuery. If the layer refuses, or the
     adapter does not admit the value, there is no inferred profile and the reason is returned. Nothing here works
     around a refusal;
  4. it returns a provenance record with every answer: the layer's codes and notices, whether the value is the
     engine's estimate or a flagged baseline-level value (and which comparator), the fingerprint and training
     groups of the producing state, the SHA-256 of run_complete.json and of the code files.

The adapter's status stays 'blocked: independent validation pending'. Nothing here unblocks it, and nothing here
calls a backend: the bridge only reads what the adapter would hand to one.

Worlds. The engine's real world maps to the layer's real world (core.real_world). On the real cohort this file
NEVER freezes, runs, scores or serves the layer's benchmark: it only reads a finished run, and when there is none
it says 'general_layer_run_missing'. An engine SYNTHETIC world maps to a layer world made from that world's own
tables, TRAINING ROWS ONLY (convert_synthetic copies no held-out row of any table), on which the layer's own
benchmark machinery is run with the layer's own small declared grid (ensure_synthetic_run).

It holds no fitted state of its own, opens no trace, and reads no regional row of a validation or test brain
sample: the only regional rows it touches are training rows (the converter copies their text; the layer's sealed
loader parses them).

    python3 -I analysis/synarch_axon_bridge.py status [--world DIR]      # read-only
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from collections import namedtuple
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


def _sibling(name):
    """Load a sibling file by path under a fixed name (python3 -I does not put this folder on sys.path)."""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


bench = _sibling("synarch_recon_benchmark")       # the layer's commands; it loads the core and the modules
core, mods = bench.core, bench.mods

ROOT = core.ROOT
SCHEMA = "synarch_axon_bridge_v1"
BRIDGE_FILE = "analysis/synarch_axon_bridge.py"
LAYER_CODE_FILES = ("analysis/synarch_recon.py", "analysis/synarch_recon_modules.py", "analysis/synarch_recon_benchmark.py")
NAMESPACE = "jiao"
TARGET = "regional_axon_length"                    # the inferred property the axon engine may receive
TASK = "T1"                                        # the layer's declared task: profile from soma position + source label
T1_INPUTS = ("soma_position", "source_label")
MODULE = "joint_mixture"
ADAPTER = "axon_arbor_adapter"
BLOCKED = "blocked: independent validation pending"
RUN_MISSING = "general_layer_run_missing"          # the reason code when no finished run exists for a world
LAYER_FOLDER = "general_layer"                     # the layer world of an engine synthetic world, under the engine's output
AUTO = object()

Stated = namedtuple("Stated", "value scale transform frame coordinates")      # what core.model_values reads of a value


class BridgeRefused(core.ReconError):
    """The bridge cannot be built or used as asked: the finished run does not verify, the two worlds do not
    correspond, or a forbidden record was offered. Never a silent fallback."""


class BridgeLeak(BridgeRefused):
    """A validation or test brain sample, or the query's own brain sample, is in a state that would produce an
    inferred input."""


# ----------------------------------------------------------------------------------------------
# code identity
# ----------------------------------------------------------------------------------------------


def code_hashes(root=ROOT):
    """SHA-256 of this file and of the three code files of the general layer, as they are on disk now."""
    return {relative: core.sha_file(Path(root) / relative) for relative in (BRIDGE_FILE,) + LAYER_CODE_FILES}


def layer_code_hashes(root=ROOT):
    return {relative: core.sha_file(Path(root) / relative) for relative in LAYER_CODE_FILES}


# ----------------------------------------------------------------------------------------------
# world mapping
# ----------------------------------------------------------------------------------------------


def layer_directory(engine_world):
    """Where the layer world of an engine SYNTHETIC world lives: inside the engine's own output folder, so a
    rebuilt or copied engine world carries (or lacks) its own layer run."""
    return Path(engine_world.out) / LAYER_FOLDER


def layer_world_for(engine_world):
    """The layer world that belongs to an engine world, or None when a synthetic world has not been converted.

    Real cohort: core.real_world(), after checking that both worlds read the same three tables. Synthetic: the
    converted world under <engine output>/general_layer."""
    if engine_world.kind == "real":
        world = core.real_world()
        pairs = ((engine_world.manifest, world.manifest), (engine_world.soma_table, world.soma_csv),
                 (engine_world.regional_table, world.regional_csv))
        for mine, theirs in pairs:
            if Path(mine).resolve() != Path(theirs).resolve():
                raise BridgeRefused(f"the engine's real world reads {mine}, the layer's real world reads {theirs}: "
                                    "the two do not describe the same cohort")
        return world
    if engine_world.kind != "synthetic":
        raise BridgeRefused(f"unknown kind of engine world '{engine_world.kind}'")
    directory = layer_directory(engine_world)
    if not (directory / "world.json").exists():
        return None
    return core.synthetic_world(directory)


# ----------------------------------------------------------------------------------------------
# the converter: an engine synthetic world's TRAINING rows in the layer's file formats
# ----------------------------------------------------------------------------------------------

MANIFEST_HEAD = ("neuron", "animal", "source_label", "soma_label", "class", "animal_inspected", "animal_development",
                 "animal_in_li_2024", "part", "part_reason")


def synthetic_config():
    """The grid declared for every engine synthetic world: the layer's OWN small configuration, unchanged."""
    return json.loads(json.dumps(bench.TINY_CONFIG))


def training_content(engine_world):
    """The training partition of an engine synthetic world as text, read without parsing any held-out row.

    Returns the manifest rows (dicts, training only), the soma table header and lines, the regional header and
    lines, and a SHA-256 of all of it. A line of the soma or regional table is kept only when its brain sample (the
    ID prefix) is a training sample; every other line is dropped as soon as its ID is read, unparsed."""
    with open(engine_world.manifest, newline="") as handle:
        manifest = [{k: row[k] for k in MANIFEST_HEAD} for row in csv.DictReader(handle, delimiter="\t") if row["part"] == "train"]
    samples = {row["animal"] for row in manifest}
    neurons = {row["neuron"] for row in manifest}
    if any(row["neuron"].split("_", 1)[0] != row["animal"] for row in manifest) or len(neurons) != len(manifest):
        raise BridgeRefused("the engine world's manifest has a duplicated neuron or an ID that disagrees with its sample")
    with open(engine_world.soma_table, newline="") as handle:
        soma_head = handle.readline().rstrip("\r\n")
        columns = soma_head.split(",")
        if tuple(columns) != ("file", "neurons", "geno", "location", "AP", "DV", "ML"):
            raise BridgeRefused("the engine world's soma table does not have the declared columns")
        soma = [line.rstrip("\r\n") for line in handle if line.split(",", 2)[1] in neurons]
    with open(engine_world.regional_table, newline="") as handle:
        regional_head = handle.readline().rstrip("\r\n")
        regional = []
        for line in handle:
            if line[:max(line.find(","), 0)].split("_", 1)[0] in samples:      # sealed otherwise: never parsed, never kept
                regional.append(line.rstrip("\r\n"))
    if len(soma) != len(manifest) or len(regional) != len(manifest):
        raise BridgeRefused("the engine world's soma or regional table does not hold every training neuron exactly once")
    sha = core.digest({"manifest": manifest, "soma": [soma_head] + soma, "regional": [regional_head] + regional})
    return {"manifest": manifest, "soma_head": soma_head, "soma": soma, "regional_head": regional_head,
            "regional": regional, "samples": sorted(samples), "neurons": sorted(neurons), "sha256": sha}


def convert_synthetic(engine_world, *, config=None):
    """Write the layer world of an engine SYNTHETIC world: its training rows, and nothing of a held-out row.

    Files (the layer's formats): manifest.tsv, soma.csv, soma.xlsx (columns neurons, location), regional.csv,
    world.json (regions, the declared configuration), conversion.json (what was copied, with hashes). The folder
    is created once and never overwritten."""
    import openpyxl

    if engine_world.kind != "synthetic":
        raise BridgeRefused("only an engine SYNTHETIC world is converted; the real cohort is read through core.real_world()")
    directory = layer_directory(engine_world)
    if directory.exists():
        raise BridgeRefused(f"{directory} exists; a layer world is converted once and never overwritten")
    content = training_content(engine_world)
    regions = tuple(content["regional_head"].split(",")[1:])
    if len(set(regions)) != len(regions) or not set(regions) <= set(core.JIAO_REGIONS):
        raise BridgeRefused("the engine world's regional columns are not registered region names")
    observed = tuple(r for r in regions if r in core.JIAO_HYPOTHALAMIC)
    directory.mkdir(parents=True)
    with open(directory / "manifest.tsv", "w", newline="") as handle:
        handle.write("\t".join(MANIFEST_HEAD) + "\n")
        for row in content["manifest"]:
            handle.write("\t".join(row[k] for k in MANIFEST_HEAD) + "\n")
    (directory / "soma.csv").write_text("\n".join([content["soma_head"]] + content["soma"]) + "\n")
    (directory / "regional.csv").write_text("\n".join([content["regional_head"]] + content["regional"]) + "\n")
    workbook = openpyxl.Workbook(write_only=True)
    sheet = workbook.create_sheet("Sheet1")
    sheet.append(["neurons", "location"])
    for row in content["manifest"]:
        sheet.append([row["neuron"], row["soma_label"]])
    workbook.save(directory / "soma.xlsx")
    meta = {"name": f"axon_bridge_{engine_world.name}", "synthetic": True, "regions": list(regions),
            # the one spelling difference the engine's cohort knows: soma table 'Wt(orexin)', manifest 'Orexin'
            "t4_observed": list(observed), "label_alias": [["Wt(orexin)", "Orexin"]], "species": "synthetic",
            "expected": {"train_samples": len(content["samples"]), "train_neurons": len(content["neurons"])},
            "config": config if config is not None else synthetic_config(),
            "made_by": f"{BRIDGE_FILE} from the TRAINING rows of the engine synthetic world '{engine_world.name}'"}
    core.write_json(directory / "world.json", meta, pretty=True)
    files = ("manifest.tsv", "soma.csv", "soma.xlsx", "regional.csv", "world.json")
    record = {"schema": SCHEMA, "engine_world": engine_world.name, "training_content_sha256": content["sha256"],
              "training_samples": len(content["samples"]), "training_neurons": len(content["neurons"]),
              "held_out_rows_copied": 0, "what": "training rows only, of every table; held-out rows are not copied",
              "config_is_the_layers_own_small_grid": config is None,
              "files_sha256": {name: core.sha_file(directory / name) for name in files}}
    core.write_json(directory / "conversion.json", record, pretty=True)
    return core.synthetic_world(directory)


def verify_conversion(engine_world):
    """The layer world is still the conversion of THIS engine world's training rows, and its files are the ones
    that were written. Anything else is refused: a state fitted on other rows would not be this world's state."""
    directory = layer_directory(engine_world)
    path = directory / "conversion.json"
    if not path.exists():
        raise BridgeRefused(f"{directory} was not made by the bridge's converter (conversion.json is missing)")
    record = core.read_json(path)
    for name, sha in record["files_sha256"].items():
        if not (directory / name).exists() or core.sha_file(directory / name) != sha:
            raise BridgeRefused(f"{directory / name} is not the file the converter wrote")
    content = training_content(engine_world)
    if content["sha256"] != record["training_content_sha256"]:
        raise BridgeRefused("the layer world is not the conversion of this engine world's training rows: the training "
                            "tables changed after the conversion (the layer's states were fitted on other rows)")
    listed = (directory / "manifest.tsv").read_text().splitlines()[1:]
    if sorted(line.split("\t")[0] for line in listed) != content["neurons"]:
        raise BridgeRefused("the converted manifest does not list exactly the engine world's training neurons")
    return record


def ensure_synthetic_run(engine_world, *, log=None):
    """On an engine SYNTHETIC world: convert (once), freeze, and run the layer's benchmark once (its run fits the
    all-training serving states and binds them in its completion record), with the layer's own commands.
    Returns run_status(). Refused on the real cohort, always."""
    if engine_world.kind != "synthetic":
        raise BridgeRefused("the axon bridge never freezes, runs, scores or serves the general layer's benchmark on the "
                            "real cohort; it only reads a finished run (the run is the layer owner's single run)")
    quiet = log or (lambda *parts: None)
    directory = layer_directory(engine_world)
    if not directory.exists():
        convert_synthetic(engine_world)
    verify_conversion(engine_world)
    world = core.synthetic_world(directory)
    out = Path(world.out)
    if not (out / "pre_run_declaration.json").exists():
        bench.freeze(world, "synthetic mechanism check through the axon bridge: no review applies", log=lambda *a: None)
    if not (out / "run_complete.json").exists():
        if (out / "run_started.json").exists():
            raise BridgeRefused(f"{out}: a run was started and did not complete; the layer's benchmark runs once per world")
        # the layer's single run: folds, selection, the all-training states, scoring, the completion record that
        # binds all of it (post=False leaves out its demonstration and its report, which are not needed here)
        bench.run(world, workers=1, directory=str(directory), log=quiet, post=False)
    return run_status(engine_world)


# ----------------------------------------------------------------------------------------------
# the finished run: what is pinned, and whether it verifies
# ----------------------------------------------------------------------------------------------


def _fold_key(fold):
    return "serving" if fold is None else f"outer{int(fold)}"


def _serving_record(out):
    """The record of the all-training states that the completion record of a finished run BINDS (it is held inside
    run_complete.json; serving.json is only a copy of it), or None when the run binds none."""
    key = getattr(bench, "SERVING_KEY", None)
    if key is None:
        raise BridgeRefused("this bridge was written for a general layer whose completion record binds its all-training "
                            "states; the layer on disk has no such record")
    return core.read_json(out / "run_complete.json").get(key)


def _recorded_state(out, fold, selection=None):
    """The T1 state of a fold as the finished run RECORDED it: fingerprint, training groups, file. For the
    all-training states the record is the one inside run_complete.json, never a loose file."""
    if fold is None:
        serving = _serving_record(out)
        if serving is None:
            raise BridgeRefused("the completion record of the finished run binds no all-training states: held-out "
                                "queries cannot be answered by it")
        chosen = serving["tasks"][TASK]
        if chosen.get("failed"):
            raise BridgeRefused(f"the all-training fit for {TASK} failed in the run: {chosen['failed']}")
    else:
        selection = selection if selection is not None else core.read_json(out / "selection.json")
        chosen = selection["engine_used"][TASK][str(int(fold))]
    name = bench.state_file("serving" if fold is None else f"outer{int(fold)}", chosen["representation"], chosen["components"])
    return {"state_hash": chosen["state_hash"], "training_groups_sha256": chosen["training_groups_sha256"],
            "representation": chosen["representation"], "components": chosen["components"], "file": f"states/{name}"}


def run_pins(engine_world, folds=(None,)):
    """Everything of the general layer that an E7 forecast of this world depends on, as it is on disk NOW.

    {'run': None, 'reason': 'general_layer_run_missing'} when no finished run exists. Otherwise: the SHA-256 of
    run_complete.json (after the layer's own check that every file it binds still has its recorded hash and that
    serving.json, if present, is the exact copy of the bound record), of the declaration, the selection and the
    evidence; the T1 status and best comparator; and for every needed fold the recorded state fingerprint with
    the SHA-256 of its file. `folds` is a list of outer folds (None = the all-training serving state), or 'all':
    every outer fold, and the serving state when the run bound one. Loads no outcome and fits nothing."""
    world = layer_world_for(engine_world)
    out = None if world is None else Path(world.out)
    if out is None or not (out / "run_complete.json").exists():
        return {"run": None, "reason": RUN_MISSING}
    try:
        record = bench.verify_run(world)
        evidence = core.read_json(out / "evidence.json")["tasks"][TASK]
        selection = core.read_json(out / "selection.json")
        if isinstance(folds, str):
            if folds != "all":
                raise BridgeRefused("folds is a list of outer folds or 'all'")
            serving = _serving_record(out)
            usable = serving is not None and not serving["tasks"][TASK].get("failed")
            folds = sorted(int(k) for k in selection["engine_used"][TASK]) + ([None] if usable else [])
        states = {}
        for fold in folds:
            entry = _recorded_state(out, fold, selection)
            file = out / entry["file"]
            if not file.exists():
                raise BridgeRefused(f"{file}: the state file of {_fold_key(fold)} is missing")
            states[_fold_key(fold)] = {**entry, "file_sha256": core.sha_file(file)}
    except BridgeRefused:
        raise
    except core.ReconError as error:
        raise BridgeRefused(f"the finished run of the general layer does not verify: {error}") from error
    pins = {"world": world.name, "task": TASK, "run_complete_sha256": core.sha_file(out / "run_complete.json"),
            "pre_run_declaration_sha256": record["sha256"]["pre_run_declaration.json"],
            "selection_sha256": record["sha256"]["selection.json"], "evidence_sha256": record["sha256"]["evidence.json"],
            "status": evidence["status"], "best_comparator": evidence["best_comparator"],
            "best_comparator_params": evidence["best_comparator_params"], "states": states}
    if engine_world.kind == "synthetic":
        pins["conversion_sha256"] = core.sha_file(layer_directory(engine_world) / "conversion.json")
    return {"run": pins, "reason": None}


def run_status(engine_world):
    """A readable summary for reports and commands: is there a finished run, does it verify, what does T1 return."""
    info = {"engine_world": engine_world.name, "kind": engine_world.kind, "code_sha256": code_hashes()}
    try:
        world = layer_world_for(engine_world)
        info["layer_world"] = None if world is None else world.name
        info["layer_out"] = None if world is None else str(world.out)
        pins = run_pins(engine_world, folds=())
        info["finished_run"] = pins["run"] is not None
        if pins["run"] is None:
            info["reason"] = pins["reason"]
            return info
        if engine_world.kind == "synthetic":
            verify_conversion(engine_world)
        bench.verify_freeze(world)
        info.update({"verified": True, **{k: pins["run"][k] for k in ("run_complete_sha256", "status", "best_comparator",
                                                                      "best_comparator_params")}})
    except core.ReconError as error:
        info.update({"finished_run": info.get("finished_run", True), "verified": False, "problem": str(error)})
    return info


def run_artifacts(engine_world, folds=(None,)):
    """[(name, path)] of the files of the finished run that a sealed bundle should carry: the completion record
    (it holds the record of the all-training states), the declaration, the selection, the evidence and the state
    files of the needed folds."""
    world = layer_world_for(engine_world)
    out = None if world is None else Path(world.out)
    if out is None or not (out / "run_complete.json").exists():
        return []
    names = ["run_complete.json", "pre_run_declaration.json", "selection.json", "evidence.json"]
    names += [_recorded_state(out, fold)["file"] for fold in folds]
    return [(name, out / name) for name in names]


# ----------------------------------------------------------------------------------------------
# answers
# ----------------------------------------------------------------------------------------------


@dataclass
class Answer:
    """What the bridge returns for one query neuron.

    query       the AxonQuery the adapter stated (permitted measured inputs and the admitted inferred profile), or
                None when there is no inferred profile;
    profile     None, or {'names', 'value', 'scale', 'transform', 'frame', 'state_hash'}: the inferred regional
                profile exactly as the adapter stated it, with the region NAME of every value;
    reason      None, or the reason code (the layer's refusal code, or 'general_layer_run_missing');
    provenance  small JSON-plain record (see provenance_fields); upstream_groups: the training brain samples of the
                state that was consulted (empty when no state was consulted)."""
    query: object
    profile: dict | None
    reason: str | None
    provenance: dict
    upstream_groups: tuple = ()
    record: object = field(default=None, repr=False)        # the layer's own inferred record (for checks; not saved)


def _record_problem(spec, record):
    """The layer's ONE check of a record against its registered variable: its statement of scale, tag and frame,
    and (where the layer on disk has the check) its own invariants and its value against the frame it states."""
    check = getattr(core, "record_problem", None)
    if check is not None:
        return check(spec, record)
    return core.statement_problem(spec, record.scale, record.transform, record.frame)


def _short(text, limit=240):
    text = " ".join(str(text).split())
    return text if len(text) <= limit else text[:limit - 3] + "..."


class Session:
    """The bridge of ONE engine world. Built once per engine process; it holds reconstructors, never a fitted state
    of its own.

    `heldout_groups` are the validation and test brain sample IDs of the engine's cohort (permitted metadata): a
    state that holds one of them is refused, whatever its files say about themselves.
    `expected_run_complete_sha256`, when given, is the anchor: the SHA-256 of the completion record that the caller
    registered earlier. The layer's own file checks cannot tell a finished run from one whose every file was
    rewritten consistently; a completion record with another SHA-256 than the registered one is refused here."""

    def __init__(self, engine_world, *, heldout_groups=(), expected_run_complete_sha256=None):
        self.engine_world = engine_world
        self.expected_run = expected_run_complete_sha256
        self.heldout = frozenset(str(g) for g in heldout_groups)
        self.heldout_folded = frozenset(core.fold_id(g) for g in self.heldout)
        self.ledger = core.ProvenanceLedger()
        self.code_sha256 = code_hashes()
        self.world = layer_world_for(engine_world)
        self._ctx, self._layers, self._source_sha = None, {}, {}
        self.out = None if self.world is None else Path(self.world.out)
        self.ready = self.out is not None and (self.out / "run_complete.json").exists()
        self.reason = None if self.ready else RUN_MISSING
        self.run_complete_sha256 = None
        self.fold_of_group, self.training_groups = {}, frozenset()
        if self.ready:
            self._open()

    # -- verification ----------------------------------------------------------------------------
    def _open(self):
        """Verify the finished run before anything is asked of it: the conversion (synthetic), the declaration
        against the code and inputs on disk, every file the completion record binds. Loads no outcome."""
        try:
            if self.engine_world.kind == "synthetic":
                verify_conversion(self.engine_world)
            declaration = bench.verify_freeze(self.world)
            bench.verify_run(self.world)
        except BridgeRefused:
            raise
        except core.ReconError as error:
            raise BridgeRefused(f"the finished run of the general layer does not verify: {error}") from error
        if TASK not in declaration["config"]["tasks"]:
            raise BridgeRefused(f"the finished run did not benchmark {TASK}: there is no declared profile forecast to ask for")
        self.run_complete_sha256 = core.sha_file(self.out / "run_complete.json")
        if self.expected_run is not None and self.run_complete_sha256 != self.expected_run:
            raise BridgeRefused("the completion record of the finished run is not the registered one (run_complete.json has "
                                f"SHA-256 {self.run_complete_sha256}, registered: {self.expected_run}): a forged or "
                                "replaced run is refused")
        for fold, held in declaration["splits"]["outer"].items():
            for group in held:
                self.fold_of_group[str(group)] = int(fold)
        self.training_groups = frozenset(self.fold_of_group)
        leak = sorted(g for g in self.training_groups if core.fold_id(g) in self.heldout_folded)
        if leak:
            raise BridgeLeak(f"the general layer's training partition holds validation or test brain samples of the "
                             f"engine's cohort: {leak[:5]}")
        self.selection = core.read_json(self.out / "selection.json")

    def layer(self, fold):
        """(reconstructor, joint module, T1 state, adapter) of a fold, from build_reconstructor ONLY."""
        if not self.ready:
            raise BridgeRefused(f"no finished run of the general layer exists for {self.engine_world.name}")
        if fold not in self._layers:
            if fold is not None and int(fold) not in set(self.fold_of_group.values()):
                raise BridgeRefused(f"the finished run has no outer fold {fold}")
            try:
                recorded = _recorded_state(self.out, fold, self.selection)
                reconstructor, ctx = bench.build_reconstructor(self.world, fold=fold, ctx=self._ctx, ledger=self.ledger)
            except BridgeRefused:
                raise
            except core.ReconError as error:
                raise BridgeRefused(f"the general layer refuses to serve the finished run: {error}") from error
            self._ctx = ctx
            if reconstructor.policy != core.Policy():
                raise BridgeRefused("the reconstructor does not run under the layer's default policy")
            module, adapter = reconstructor.module(MODULE), reconstructor.module(ADAPTER)
            state = module.states[TASK]
            if state.state_hash() != recorded["state_hash"] or state.groups_sha256 != recorded["training_groups_sha256"]:
                raise BridgeRefused(f"the {TASK} state of {_fold_key(fold)} is not the state the finished run recorded")
            groups = {core.fold_id(g) for g in state.groups}
            if groups & self.heldout_folded:
                raise BridgeLeak(f"the {TASK} state of {_fold_key(fold)} holds validation or test brain samples: "
                                 f"{sorted(g for g in state.groups if core.fold_id(g) in self.heldout_folded)[:5]}")
            if not set(state.groups) <= self.training_groups:
                raise BridgeLeak(f"the {TASK} state of {_fold_key(fold)} holds brain samples outside the layer's training partition")
            if fold is not None and any(self.fold_of_group[g] == int(fold) for g in state.groups):
                raise BridgeLeak(f"the state of outer fold {fold} holds a brain sample of that fold")
            if adapter.STATUS != BLOCKED or adapter.evidence_for(adapter.conditionals()[0].name).status != BLOCKED:
                raise BridgeRefused("the axon adapter is not in its blocked state; this bridge was written for a blocked adapter")
            self._layers[fold] = (reconstructor, module, state, adapter)
        return self._layers[fold]

    @property
    def registry(self):
        if self._ctx is not None:
            return self._ctx.registry
        if not hasattr(self, "_registry"):
            self._registry = core.Registry.from_contract(ROOT / bench.CONTRACT)
        return self._registry

    # -- which state may answer for a neuron ------------------------------------------------------
    def fold_for(self, row):
        """None (the all-training serving states) for a neuron outside the training partition; the layer's outer
        fold that holds the neuron's brain sample for a training neuron. The engine's part and the layer's
        partition must agree about the brain sample, or nothing is served."""
        animal = str(row.animal)
        in_layer = animal in self.training_groups
        if (row.part == "train") != in_layer:
            raise BridgeLeak(f"{row.neuron}: the engine's cohort has brain sample {animal} in part '{row.part}', the "
                             f"general layer's partition {'holds' if in_layer else 'does not hold'} it as a training sample")
        return self.fold_of_group[animal] if in_layer else None

    def assert_isolated(self, row, state, fold):
        """The producing state never saw the query's brain sample, and holds no held-out brain sample."""
        animal = core.fold_id(row.animal)
        if animal in {core.fold_id(g) for g in state.groups} or core.id_in(row.neuron, state.units, state.units_folded):
            raise BridgeLeak(f"query {row.neuron}: its brain sample {row.animal} is in the {TASK} state of {_fold_key(fold)}")
        if {core.fold_id(g) for g in state.groups} & self.heldout_folded:
            raise BridgeLeak(f"the {TASK} state of {_fold_key(fold)} holds validation or test brain samples")

    # -- records ---------------------------------------------------------------------------------
    def _sha(self, path):
        path = str(path)
        if path not in self._source_sha:
            self._source_sha[path] = core.sha_file(path)
        return self._source_sha[path]

    def group_of(self, row):
        return core.jiao_group(row.animal, self.registry)

    def measured_records(self, row, *, group=None):
        """The measured records of the PERMITTED inputs of one neuron, and nothing else: {variable: record}."""
        registry = self.registry
        namespace = registry.namespace(NAMESPACE)
        group = self.group_of(row) if group is None else group
        world = self.engine_world
        values = {"soma_position": (tuple(float(v) for v in row.soma), core.SOMA_AXES, world.soma_table),
                  "source_label": (str(row.tag), (), world.soma_table),          # the soma table's own spelling
                  "soma_label_fine": (row.soma_label, (), world.manifest),
                  "hemisphere": (row.hemisphere, (), world.local_manifest),
                  "reconstruction_type": (row.recon_type, (), world.local_manifest)}
        records = {}
        for variable in core.HIDDEN_AXON_INPUTS:
            value, coordinates, source = values[variable]
            if not registry.has_variable(NAMESPACE, variable) or value is None or value == "":
                continue
            spec = registry.variable(NAMESPACE, variable)
            record = core.Observation(self.world.species, NAMESPACE, namespace.grain, str(row.neuron), group, variable,
                                      "measured", value, core.scale_statement(spec, "identity"), tuple(coordinates),
                                      spec.frame, core.Uncertainty("none_stated"),
                                      core.SourceRef(str(source), self._sha(source), str(row.full_id)),
                                      transform="identity")
            problem = core.statement_problem(spec, record.scale, record.transform, record.frame)
            if problem:
                raise BridgeRefused(f"a permitted input does not state its registered scale and frame: {problem}")
            records[variable] = record
        return records

    def _request(self, row, records, target, group=None):
        namespace = self.registry.namespace(NAMESPACE)
        group = self.group_of(row) if group is None else group
        return core.ReconstructionRequest(NAMESPACE, namespace.grain, str(row.neuron), group, tuple(records),
                                          (core.Target(target),), self.world.species)

    # -- the two steps ---------------------------------------------------------------------------
    def offered_problem(self, row, records, *, group=None, fold=None):
        """Why a set of records may not be offered for a hidden-axon forecast of this neuron: (code, reason), or
        None when it may. The layer's ONE rule (core.hidden_axon_problem: code 'forbidden_input') and its
        statement check of scale, tag and frame (code 'invalid_request'), applied to every record."""
        reconstructor = self.layer(fold)[0]
        group = self.group_of(row) if group is None else group
        key = (NAMESPACE, self.registry.namespace(NAMESPACE).grain, str(row.neuron))
        problem = core.hidden_axon_problem(self.registry, NAMESPACE, key, group.group_id, dict(records),
                                           ledger=self.ledger, resolve=reconstructor.resolve_state, axon_hidden=True)
        if problem:
            return "forbidden_input", problem
        for variable, record in records.items():
            problem = _record_problem(self.registry.variable(NAMESPACE, variable), record)
            if problem:
                return "invalid_request", problem
        return None

    def ask(self, row, fold, *, records=None, group=None):
        """The layer's answer to T1 for one neuron by the states of `fold`: (request, PropertyResult), or
        (None, Refusal) when what was offered is not exactly the permitted T1 evidence.

        `records` defaults to the neuron's measured soma position and source label. Anything else offered here is
        held to the hidden-axon rule first: a measured regional profile of the neuron, whole or partial, is
        refused before the layer is asked (with it the layer would answer another, declared profile-input task)."""
        reconstructor = self.layer(fold)[0]
        if records is None:
            measured = self.measured_records(row, group=group)
            records = {v: measured[v] for v in T1_INPUTS if v in measured}
        problem = self.offered_problem(row, records, group=group, fold=fold)
        if problem:
            return None, core.Refusal(*problem)
        if set(records) != set(T1_INPUTS):
            return None, core.Refusal("invalid_request", f"the profile forecast takes exactly {' and '.join(T1_INPUTS)}; "
                                      f"offered: {', '.join(sorted(records)) or 'nothing'}")
        request = self._request(row, records.values(), TARGET, group)
        return request, reconstructor.reconstruct(request)[TARGET]

    def admit(self, row, fold, measured, inferred, *, group=None):
        """Pass the measured records and ONE inferred record through the axon adapter of the fold's reconstructor.
        Returns (AxonQuery, None), or (None, why) when the adapter does not admit them."""
        reconstructor, _, _, adapter = self.layer(fold)
        evidence = dict(measured)
        if inferred is not None:
            if inferred.variable in evidence:
                return None, f"{inferred.variable} is offered twice (measured and inferred)"
            evidence[inferred.variable] = inferred
        request = self._request(row, evidence.values(), "axon_arbor", group)
        target = core.Target("axon_arbor")
        try:
            plan = adapter.plan(request, target, evidence, reconstructor.policy)
        except core.ForbiddenInputError as error:
            return None, str(error)
        if isinstance(plan, core.Refusal):
            return None, plan.reason
        if adapter.STATUS != BLOCKED:
            raise BridgeRefused("the axon adapter changed its status while the bridge was running")
        return adapter.query_for(plan), None

    # -- the whole chain for one neuron ---------------------------------------------------------------
    def base_provenance(self):
        return {"schema": SCHEMA, "used": False, "reason_code": self.reason, "refused_by": None if self.ready else "bridge",
                "reason": None if self.ready else f"no finished run of the general layer exists for {self.engine_world.name}",
                "task": TASK, "conditional": None, "fold": None, "level": None, "comparator": None, "comparator_params": None,
                "validation": None, "codes": [], "notices": [], "coverage_validated": False, "module": None,
                "value_state_hash": None, "state_sha256": None, "training_groups_sha256": None, "training_groups": 0,
                "record_sha256": None, "adapter_status": BLOCKED, "run_complete_sha256": self.run_complete_sha256,
                "code_sha256": core.digest(self.code_sha256)}

    def infer(self, row, *, fold=AUTO, group=None):
        """The inferred regional profile of one query neuron as a permitted query input, or the reason there is none.

        `row` is the engine's cohort row (identity, part, table soma, labels, hemisphere, reconstruction type): the
        bridge takes nothing else of the neuron. `fold` and `group` exist for checks; the engine passes neither."""
        provenance = self.base_provenance()
        if not self.ready:
            return Answer(None, None, RUN_MISSING, provenance)
        chosen = self.fold_for(row) if fold is AUTO else fold
        reconstructor, module, state, adapter = self.layer(chosen)
        if group is None:
            self.assert_isolated(row, state, chosen)
        groups = tuple(state.groups)
        provenance.update({"fold": _fold_key(chosen), "state_sha256": state.state_hash(),
                           "training_groups_sha256": state.groups_sha256, "training_groups": len(groups)})
        measured = self.measured_records(row, group=group)
        request, result = self.ask(row, chosen, group=group)

        def refused(by, code, why):
            provenance.update({"reason_code": code, "refused_by": by, "reason": _short(why)})
            return Answer(None, None, code, provenance, groups)

        if request is None:
            return refused("bridge", result.code, result.reason)
        provenance.update({"conditional": result.conditional, "validation": result.validation,
                           "codes": sorted(result.codes()), "module": result.module,
                           "notices": [[n.code, n.handling] for n in result.notices],
                           "coverage_validated": bool(result.coverage_validated), "value_state_hash": result.state_hash})
        if result.status != "inferred":
            return refused("layer", result.refusal.code, result.refusal.reason)
        declared = module.tasks[TASK].conditional().name
        used = dict(result.evidence_used)
        if result.conditional != declared or not set(used) <= set(T1_INPUTS) or set(used.values()) - {"measured"}:
            return refused("bridge", "invalid_request", f"the layer did not answer its declared task {TASK} from measured "
                                                        f"{' and '.join(T1_INPUTS)} (it answered {result.conditional})")
        record = result.to_observation(request)                 # the layer's own inferred record
        level = "baseline" if result.baseline_level else "engine"
        comparator = result.module.split(":", 1)[1] if result.baseline_level else None
        evidence = module.evidence_for(declared)
        provenance.update({"level": level, "comparator": comparator, "record_sha256": core.record_digest(record),
                           "comparator_params": evidence.best_comparator_params
                           if comparator is not None and comparator == evidence.best_comparator else None})
        if level == "engine" and result.state_hash != state.state_hash():
            return refused("bridge", "invalid_estimate", "the engine-level value does not name the registered state")
        query, why = self.admit(row, chosen, measured, record, group=group)
        if query is None:
            return refused("adapter", "forbidden_input", why)
        stated = query.inferred.get(TARGET)
        if stated is None or set(query.inferred) != {TARGET} or tuple(stated["value"]) != tuple(record.value) \
                or stated["state_hash"] != result.state_hash:
            return refused("bridge", "invalid_estimate", "the adapter did not state exactly the issued profile")
        if set(query.inputs) != set(measured) or len(record.coordinates) != len(stated["value"]):
            return refused("bridge", "invalid_estimate", "the adapter's query does not hold exactly the permitted inputs")
        profile = {"names": tuple(record.coordinates), **stated}
        provenance.update({"used": True, "reason_code": None, "refused_by": None, "reason": None,
                           "adapter_status": adapter.STATUS})
        return Answer(query, profile, None, provenance, groups, record)


def open_session(engine_world, *, heldout_groups=(), expected_run_complete_sha256=None):
    return Session(engine_world, heldout_groups=heldout_groups, expected_run_complete_sha256=expected_run_complete_sha256)


# ----------------------------------------------------------------------------------------------
# the similarity between an inferred profile and measured donor profiles (what E7 multiplies its weights by)
# ----------------------------------------------------------------------------------------------


def _target_spec(registry=None):
    registry = registry if registry is not None else _default_registry()
    return registry.variable(NAMESPACE, TARGET)


_REGISTRY = {}


def _default_registry():
    if "registry" not in _REGISTRY:
        _REGISTRY["registry"] = core.Registry.from_contract(ROOT / bench.CONTRACT)
    return _REGISTRY["registry"]


def model_values_by_name(stated, names, *, registry=None):
    """{region name: value on the layer's MODEL scale} of one stated profile.

    The values are converted by the layer's single conversion function (core.model_values), exactly once and
    according to the value's own stated tag: a value tagged 'identity' is in published units and is transformed;
    a value tagged 'log10_1p' is already on the model scale and is returned as it is. A missing, unknown or
    contradictory scale statement, other units or another frame raise core.ScaleError. Region names must be
    distinct: values are matched by NAME, never by position."""
    names = tuple(str(n) for n in names)
    if len(set(names)) != len(names):
        raise BridgeRefused("a region name is duplicated: values are matched by name and the match would be ambiguous")
    try:
        count = len(stated.value)
    except TypeError:
        raise BridgeRefused("a profile is one value per region name") from None
    if count != len(names):
        raise BridgeRefused(f"{count} values for {len(names)} region names")
    values = core.model_values(_target_spec(registry), stated._replace(coordinates=names))
    return dict(zip(names, (float(v) for v in values)))


def profile_similarity(inferred, donor_values, donor_names, *, donor_tag="log10_1p", registry=None):
    """The factor s_j by which E7 multiplies the weight of every neighbour.

    inferred      {'names', 'value', 'scale', 'transform', 'frame'}: the query's inferred profile as the adapter
                  stated it (Answer.profile);
    donor_values  one row per neighbour: its MEASURED regional profile, in the order of `donor_names`;
    donor_tag     the value tag of the donor rows: 'log10_1p' (the engine's sealed loader stores log10(1 + published
                  length)) or 'identity' (published lengths).

    Both sides are brought to the model scale by model_values_by_name, each exactly once. Regions are matched by
    name; a missing, extra or duplicated name is an error. The distance is computed in the sorted order of the
    region names, so the result does not depend on the order in which either side lists its regions.
    d_j = Euclidean distance over the regions; tau = median of d_j; s_j = exp(-d_j^2 / (2 tau^2)), and every
    s_j is 1 when tau is 0."""
    spec = _target_spec(registry)
    p = model_values_by_name(Stated(inferred["value"], inferred["scale"], inferred["transform"], inferred["frame"], ()),
                             inferred["names"], registry=registry)
    donor_names = tuple(str(n) for n in donor_names)
    if len(set(donor_names)) != len(donor_names):
        raise BridgeRefused("a donor region name is duplicated: values are matched by name")
    if set(p) != set(donor_names):
        missing, extra = sorted(set(donor_names) - set(p)), sorted(set(p) - set(donor_names))
        raise BridgeRefused(f"the inferred profile and the donor profiles do not hold the same regions by name "
                            f"(missing in the inferred profile: {missing[:5]}; not in the donor profiles: {extra[:5]})")
    order = sorted(donor_names)
    scale = core.scale_statement(spec, donor_tag)
    rows = []
    for values in donor_values:
        q = model_values_by_name(Stated(tuple(float(v) for v in values), scale, donor_tag, spec.frame, ()), donor_names,
                                 registry=registry)
        rows.append([q[name] for name in order])
    if not rows:
        raise BridgeRefused("no donor profile was given")
    Q = np.array(rows, dtype=float)
    target = np.array([p[name] for name in order], dtype=float)
    d = np.sqrt(((Q - target[None, :]) ** 2).sum(axis=1))
    tau = float(np.median(d))
    s = np.ones(len(d)) if tau == 0.0 else np.exp(-d ** 2 / (2.0 * tau * tau))
    return {"s": s, "d": d, "tau": tau, "regions": len(order)}


# ----------------------------------------------------------------------------------------------
# command line (read-only)
# ----------------------------------------------------------------------------------------------


def main(argv=None):
    parser = argparse.ArgumentParser(description="Status of the general layer's finished run as the axon bridge sees it (read-only).")
    parser.add_argument("command", choices=["status"])
    parser.add_argument("--world", default="real", help="'real' (default) or the directory of an engine synthetic world")
    args = parser.parse_args(argv)
    World = namedtuple("World", "name kind manifest soma_table regional_table out")
    if args.world == "real":
        real = core.real_world()
        world = World("jiao2025", "real", str(real.manifest), str(real.soma_csv), str(real.regional_csv),
                      str(ROOT / "data/interim/axon_engine/real"))
    else:
        data = Path(args.world) / "data"
        world = World(core.read_json(data / "world.json")["name"], "synthetic", str(data / "manifest.tsv"),
                      str(data / "soma.csv"), str(data / "regional.csv"), str(Path(args.world) / "engine"))
    print(json.dumps(core.jsonable(run_status(world)), indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
