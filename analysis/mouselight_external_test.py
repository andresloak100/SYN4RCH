"""Prospective, soma-only MouseLight transfer of frozen Jiao TRAIN geometry recipes.

No network, fitting, old-cohort writes, or external geometry I/O on import, audit,
freeze or predict. Production lifecycle (fixed paths/run ID): audit; ROOT resolves
the draft admission; freeze; predict; external acquisition; score (once).

This is a new experiment. E4 failed its previous full gate; it is not a validated
winner. E6 is an additional fixed candidate. Missing molecular labels stay missing.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
RUN_ID = "mouselight_external_geometry_v1"
OUT = Path("data/interim/mouselight_external_test") / RUN_ID
RAW = Path("data/raw/mouselight_external_test") / RUN_ID
PROTOCOL = Path("docs/mouselight_external_protocol.json")
ACQUISITION_CODE = "analysis/mouselight_external_acquire.py"
ACQUISITION_TESTS = "tests/test_mouselight_external_acquire.py"
METADATA = Path("data/interim/paired_projection_search_round6/image_routes/mouselight_hy_source_manifest.json")
LIBRARY = Path("data/interim/axon_engine/real")
METHODS = ("nearest_copy", "global_medoid_copy", "E1_local_uniform", "E2_local_kernel",
           "E4_local_kernel_regional", "E6_local_consensus")
CANDIDATES = METHODS[4:]
CONTROLS = METHODS[:4]
SOURCE_PINS = {
    "analysis/axon_engine.py": "d9eac6f220691b33adc729ecc24a39d3290418039b7717bb052900741da971a8",
    "analysis/brm_template.py": "0a3da130c0e0049e8c9cc5df94523442f542f44a6815e82d3ad947a960a8d49a",
    "analysis/brm_geometry.py": "8de317bc4f261e8061b0bccad1dbba537b2a8daf6f6834d94d303645d6611187",
    "analysis/brm_joint.py": "ae7f5ec432bcbd2b0eebb65595549bc12158b1b7c047f5ace167c7840c0c937e",
    "analysis/brm_real.py": "5225ea5833b32fee57e735600eda87a9784cb2b1798cd326ee4203e1b6997648",
    "analysis/d110_geometry.py": "3b011643b1b224256fdce1ad8a31014e01e644c6a794bd1d49af5e6543090a5e",
    "data/interim/axon_engine/real/freeze/frozen.json": "a05ac68b58dd9bde92f535c4e097e573b5c35a20696d9cdccc3710c1e50b0202",
    "data/interim/axon_engine/real/ingestion_log.json": "2826ba68b02b946957826b7229712272cd9b97406e6d60a5642e8699e5245c0f",
    "data/interim/axon_engine/real/states/all_training.json": "b690e97ea65a48c9737229dbffe0501f97e85c9bdc387a0f692cec09dda98962",
    "data/interim/axon_engine/real/states/all_training.proof.json": "40ef843b647317224e31348e17314c808d85bf0c9280e535789c6e81b0367e13",
}


class Refusal(ValueError):
    pass


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write_new(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as f:
        f.write(canonical(value))


def now():
    return datetime.now(timezone.utc).isoformat()


def local_path(root, relative):
    """Only canonical, contained repository paths; no symlink escape or ../ alias."""
    p = Path(relative)
    if p.is_absolute() or ".." in p.parts or p.as_posix() != str(relative):
        raise Refusal("a canonical repository-relative path is required")
    result = Path(root) / p
    if result.resolve() != Path(root).resolve() / p:
        raise Refusal("symlink or noncanonical path refused")
    return result


def frozen_module():
    """Read-only import; none of the Engine's write commands is called."""
    name = "mouselight_external_frozen_engine"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, ROOT / "analysis/axon_engine.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


@dataclass(frozen=True)
class Query:
    neuron: str
    animal: str
    soma: tuple                         # mirrored-left AP, DV, ML micrometres
    hemisphere: str


def normalize_soma(ccfv3_xyz):
    if set(ccfv3_xyz) != {"x", "y", "z"}:
        raise Refusal("only the measured CCFv3 x/y/z soma record is accepted")
    xyz = np.array([ccfv3_xyz[k] for k in ("x", "y", "z")], dtype=float)
    if not np.isfinite(xyz).all() or (xyz < 0).any() or (xyz > [11400, 8000, 13200]).any():
        raise Refusal("CCFv3 soma is outside the declared atlas coordinate bounds")
    x, y, z = xyz
    right = x > 5700.0
    return (float(z), float(y), float(11400.0 - x if right else x)), "Right" if right else "Left"


def load_queries(root, protocol):
    """Parse permitted metadata fields only; no trace/path/image fields become features."""
    path = local_path(root, protocol["metadata"]["path"])
    if sha(path) != protocol["metadata"]["sha256"]:
        raise Refusal("metadata source pin differs")
    data = read(path)
    if data["count"] != 97 or data["sample_count"] != 22 or len(data["records"]) != 97:
        raise Refusal("source cohort differs from the metadata-only 97-neuron/22-sample snapshot")
    excluded_animals = set(protocol["admission"]["excluded_animals"])
    excluded_ids = set(protocol["admission"]["excluded_neurons"])
    out, ids, sample_numbers = [], set(), {}
    for row in data["records"]:
        neuron, sample = row["neuron"], row["sample_id"]
        if not isinstance(neuron, str) or not neuron.startswith("AA") or not neuron[2:].isdigit() or neuron in ids:
            raise Refusal("invalid or duplicate MouseLight neuron identifier")
        ids.add(neuron)
        if not isinstance(sample, str) or not sample:
            raise Refusal("explicit source sample ID required")
        number = row["sample_number"]
        if sample_numbers.setdefault(sample, number) != number:
            raise Refusal("sample UUID/number identity conflict")
        soma, hemisphere = normalize_soma(row["ccfv3_soma_um"])
        if sample not in excluded_animals and neuron not in excluded_ids:
            out.append(Query(neuron, "mouselight:" + sample, soma, hemisphere))
    if len(sample_numbers) != 22 or not out:
        raise Refusal("sample count differs or all queries were excluded")
    unknown = excluded_ids - ids
    if unknown or excluded_animals - set(sample_numbers):
        raise Refusal("declared exclusions do not identify members of the source cohort")
    return sorted(out, key=lambda q: q.neuron)


def query_record(query):
    return {"neuron": query.neuron, "animal": query.animal, "soma_ap_dv_ml_um": list(query.soma),
            "hemisphere": query.hemisphere, "source_label": None}


class FrozenLibrary:
    """Read-only adapter to admitted TRAIN arbors/cache/state. No old query proof is edited."""
    def __init__(self, root=ROOT):
        if Path(root).resolve() != ROOT.resolve():
            raise Refusal("production TRAIN assets must be the actual repository")
        verify_sources(root)                 # verify executable bytes before importing them
        self.E = frozen_module()
        self.engine = self.E.Engine(self.E.real_world(), "develop", workers=1, log=None)
        self.lib = self.engine.library       # refuses a missing/stale ingestion; never ensure/rebuild
        if len(self.lib["ids"]) != 4368 or len(set(self.lib["animal"])) != 101:
            raise Refusal("frozen TRAIN donor count changed")
        self.cache = self.engine.open_cache()  # verifies matrix bytes, no writes
        self.ids = tuple(self.lib["ids"])
        self.somas, self.animals = self.lib["soma"], tuple(self.lib["animal"].tolist())
        self.model = read(Path(root) / LIBRARY / "states/all_training.json")
        proof = read(Path(root) / LIBRARY / "states/all_training.proof.json")
        if set(proof["train_animals"]) != set(self.animals) or len(proof["train_ids"]) != 4448:
            raise Refusal("regional state and TRAIN donor animals do not agree")
        self.training_animals = frozenset(proof["train_animals"])
        self.engine.profiles                 # sealed TRAIN-only loader verifies the regional source
        self.global_medoid = None

    def donor(self, index):
        return self.engine.donor(self.ids[index])

    def distances(self, indices):
        return self.cache.sub([self.ids[i] for i in indices])[0]

    def medoid_index(self):
        if self.global_medoid is None:
            # Same ordering and matrix-vector arithmetic as the official copying control.
            weights = animal_prior(self.animals)
            self.global_medoid = int(np.argmin(self.distances(range(len(self.ids))) @ weights))
        return self.global_medoid

    def pins(self):
        self.cache.load(verify_files=True)
        actual = {}
        for n in self.ids:
            expected = self.lib["entries"][n]["canonical_file_sha256"]
            actual[n] = sha(ROOT / LIBRARY / "library" / f"{n}.npz")
            if actual[n] != expected:
                raise Refusal("canonical TRAIN arbor file pin differs: " + n)
        return {"donor_ids_sha256": self.E.ids_sha(self.ids),
                "training_library_fingerprint": self.lib["fingerprint"],
                "training_animals_sha256": self.E.ids_sha(sorted(self.training_animals)),
                "pair_cache_sha256": self.cache.sha256,
                "arbor_file_sha256": actual}


def animal_prior(animals):
    count = Counter(animals)
    return np.array([1.0 / count[a] for a in animals]) / len(count)


def attach(arbor, soma, E):
    tried = []
    for radius in (500.0, 1000.0):
        try:
            result, details = E.T.compact_attachment(arbor, soma, radius)
            return result, {"executed": "local_attachment", "radius_um": radius, "unsupported": tried,
                            "soma_shift_um": details["soma_shift_um"]}
        except E.T.UnsupportedAttachment as error:
            tried.append({"radius_um": radius, "reason": error.reason})
    return arbor.shifted(soma), {"executed": "translation", "radius_um": None, "unsupported": tried}


def predict_one(query, library):
    """The only query argument has no hidden geometry/profile/gene/label slot."""
    E = library.E
    if query.animal in library.training_animals or query.animal.removeprefix("mouselight:") in library.training_animals:
        raise Refusal("query animal is in TRAIN")
    soma = np.asarray(query.soma, float)
    d2 = ((library.somas - soma) ** 2).sum(axis=1)
    n, k = len(library.ids), min(16, len(library.ids))
    nearest_k = np.lexsort((np.arange(n), d2))[:k]
    indices = np.sort(nearest_k)       # donor IDs sorted globally; ties by ID
    donors = [library.donor(i) for i in indices]
    animals = [d["animal"] for d in donors]
    prior = animal_prior(animals)
    distance = np.sqrt(d2[indices])
    h = max(100.0, float(np.sort(distance)[min(8, k) - 1]))
    kernel = np.exp(-distance ** 2 / (2.0 * h * h))
    # Missing source tag is marginalized by the frozen posterior, not supplied as WT/Cre/nucleus.
    posterior = E.Geo.state_posterior(library.model, soma[None], [None])[0]
    responsibility = E.Geo.state_posterior(
        library.model, np.array([d["table_soma"] for d in donors]), [d["tag"] for d in donors],
        measured_profiles=np.array([d["profile"] for d in donors]))
    regional = E.T.conditional_weights(animals, responsibility, posterior)["weights"]
    raw = {"E1_local_uniform": prior, "E2_local_kernel": prior * kernel,
           "E4_local_kernel_regional": regional * kernel}
    weights = {name: w / w.sum() for name, w in raw.items()}
    D = library.distances(indices)
    choice = {name: E.T.weighted_medoid(D, w)[0] for name, w in weights.items()}
    attached = {}

    def local(i):
        if i not in attached:
            attached[i] = attach(donors[i]["arbor"], soma, E)
        return attached[i]

    result = {}
    for name, index in (("nearest_copy", int(np.argmin(d2))), ("global_medoid_copy", library.medoid_index())):
        result[name] = {"arbor": library.donor(index)["arbor"].shifted(soma),
                        "selected_donor": library.ids[index], "donor_ids": [library.ids[index]], "attachment": "translation"}
    for name in weights:
        i = choice[name]
        arbor, how = local(i)
        result[name] = {"arbor": arbor, "selected_donor": library.ids[indices[i]],
                        "donor_ids": [library.ids[j] for j in indices], "weights": weights[name].tolist(),
                        "attachment": how, "kernel_bandwidth_um": h}
    i = choice["E2_local_kernel"]
    others = [j for j in range(k) if j != i]
    arbor, details = E.consensus_prune(local(i)[0], [local(j)[0] for j in others], weights["E2_local_kernel"][others])
    result["E6_local_consensus"] = {"arbor": arbor, "selected_donor": library.ids[indices[i]],
                                    "donor_ids": [library.ids[j] for j in indices], "consensus": details,
                                    "weights": weights["E2_local_kernel"].tolist(), "attachment": local(i)[1]}
    return {name: result[name] for name in METHODS}


def protocol_of(root):
    p = read(local_path(root, str(PROTOCOL)))
    if p["schema"] != "synarch_mouselight_external_protocol_v1" or p["run_id"] != RUN_ID:
        raise Refusal("wrong experiment schema/run ID")
    if p["output_dir"] != str(OUT) or p["raw_dir"] != str(RAW) or tuple(p["methods"]) != METHODS:
        raise Refusal("canonical paths or fixed method family changed")
    return p


def verify_sources(root):
    for name, expected in SOURCE_PINS.items():
        if sha(local_path(root, name)) != expected:
            raise Refusal("frozen source/state pin differs: " + name)


def calibration_of(root, protocol, queries):
    spec = protocol["admission"]["calibration"]
    if not spec.get("receipt_path") or not spec.get("receipt_sha256"):
        raise Refusal("source-format calibration is unresolved")
    path = local_path(root, spec["receipt_path"])
    if sha(path) != spec["receipt_sha256"]:
        raise Refusal("calibration receipt pin differs")
    receipt = read(path)
    expected = {"schema": "mouselight_source_calibration_v1", "validated": True, "non_hy": True,
                "frame": "CCFv3", "axes": ["ML", "DV", "AP"], "units": "um",
                "trace_domain": "separate_axon_tracing", "root": "one_type1_root",
                "nonroot_semantics": "all_nodes_in_separate_axon_tracing_are_axon"}
    if any(receipt.get(k) != v for k, v in expected.items()):
        raise Refusal("calibration has not validated the supported native CCFv3 axon export")
    if not receipt.get("selection_rule") or receipt["selection_rule"] != spec.get("selection_rule"):
        raise Refusal("calibration selection was not declared prospectively")
    if not receipt.get("animal") or "mouselight:" + receipt["animal"] in {q.animal for q in queries}:
        raise Refusal("calibration animal overlaps evaluation or is unspecified")
    if receipt.get("neuron") in {q.neuron for q in queries}:
        raise Refusal("calibration neuron overlaps evaluation")
    if not receipt.get("allowed_node_types") or set(receipt["allowed_node_types"]) - {0, 1, 2, 5, 6}:
        raise Refusal("calibration axon node coding is unresolved")
    if not receipt.get("evidence"):
        raise Refusal("calibration evidence files are missing")
    for item in receipt["evidence"]:
        if sha(local_path(root, item["path"])) != item["sha256"]:
            raise Refusal("calibration evidence pin differs")
    return receipt


def exposure_of(root, protocol, queries):
    spec = protocol["admission"]["exposure_audit"]
    if not spec.get("path") or not spec.get("sha256"):
        raise Refusal("actual geometry-exposure audit is unresolved")
    path = local_path(root, spec["path"])
    if sha(path) != spec["sha256"]:
        raise Refusal("exposure-audit receipt pin differs")
    receipt = read(path)
    if receipt.get("schema") != "mouselight_exposure_audit_v1" or not receipt.get("passed"):
        raise Refusal("exposure audit did not pass")
    if not receipt.get("test_geometry_unopened") or receipt.get("overlapping_query_ids") != [] or receipt.get("overlapping_sample_ids") != []:
        raise Refusal("actual geometry exposure overlaps the external test or is unresolved")
    if receipt.get("query_ids_sha256") != digest([q.neuron for q in queries]):
        raise Refusal("exposure receipt does not bind the exact evaluation cohort")
    if not receipt.get("evidence"):
        raise Refusal("exposure evidence is missing")
    for item in receipt["evidence"]:
        if sha(local_path(root, item["path"])) != item["sha256"]:
            raise Refusal("exposure evidence pin differs")
    return receipt


def acquisition_recipe_files(root, protocol):
    """Freeze the exact acquisition/conversion program before any test request."""
    recipe = protocol.get("acquisition", {})
    if not recipe.get("root_reviewed") or recipe.get("root_admission") != "READY":
        raise Refusal("acquisition recipe has not passed ROOT review")
    if recipe.get("code") != ACQUISITION_CODE or recipe.get("tests") != ACQUISITION_TESTS:
        raise Refusal("acquisition source/test paths differ from the declared recipe")
    files = [ACQUISITION_CODE, ACQUISITION_TESTS]
    for role in ("source_manifest", "constants"):
        item = recipe.get(role, {})
        if not item.get("path") or not item.get("sha256"):
            raise Refusal("acquisition " + role + " pin is unresolved")
        if sha(local_path(root, item["path"])) != item["sha256"]:
            raise Refusal("acquisition " + role + " pin differs")
        files.append(item["path"])
    for relative in files:
        if not local_path(root, relative).is_file():
            raise Refusal("missing acquisition recipe file: " + relative)
    return files


def runtime_pins(root, protocol, library):
    verify_sources(root)
    # This module, its test/protocol, every imported frozen module and every TRAIN arbor are bound.
    files = dict(SOURCE_PINS)
    for name in ("analysis/mouselight_external_test.py", "tests/test_mouselight_external_test.py", str(PROTOCOL),
                 protocol["metadata"]["path"], *acquisition_recipe_files(root, protocol)):
        files[name] = sha(local_path(root, name))
    calibration = protocol["admission"]["calibration"]
    if calibration.get("receipt_path"):
        files[calibration["receipt_path"]] = sha(local_path(root, calibration["receipt_path"]))
    exposure = protocol["admission"]["exposure_audit"]
    if exposure.get("path"):
        files[exposure["path"]] = sha(local_path(root, exposure["path"]))
    data_paths = (library.engine.world.manifest, library.engine.world.soma_table,
                  library.engine.world.regional_table, library.engine.world.local_manifest)
    data_paths += tuple(item[0] for item in library.engine.world.trace_manifests)
    for name in data_paths:
        relative = str(Path(name).relative_to(Path(root)))
        files[relative] = sha(local_path(root, relative))
    # The old engine verifies its own frozen dependencies too; bind that complete list.
    for name, expected in library.E.FROZEN_PINS["source_sha256"].items():
        if sha(local_path(root, name)) != expected:
            raise Refusal("inherited frozen dependency changed: " + name)
        files[name] = expected
    return {"files": files, "library": library.pins(), "numerical_versions": {
        "numpy": np.__version__, "scipy": scipy.__version__, "python": sys.version, "threads": {n: os.environ.get(n) for n in
        ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}}}


def audit(root=ROOT):
    """Metadata-only readiness report; does not construct the donor library or write anything."""
    p = protocol_of(root)
    verify_sources(root)
    queries = load_queries(root, p)
    return {"run_id": RUN_ID, "status": p["status"], "queries": len(queries),
            "animals": len({q.animal for q in queries}), "blockers": p["admission"]["blockers"],
            "metadata_sha256": sha(local_path(root, p["metadata"]["path"])),
            "query_inputs_sha256": digest([query_record(q) for q in queries]), "geometry_opened": False}


def reviewed_commit(root):
    """Freeze only committed adapter/protocol/tests, as required by ROOT's review step."""
    command = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, check=True)
    head = command.stdout.decode().strip()
    for relative in ("analysis/mouselight_external_test.py", "tests/test_mouselight_external_test.py", str(PROTOCOL),
                     ACQUISITION_CODE, ACQUISITION_TESTS):
        result = subprocess.run(["git", "show", f"{head}:{relative}"], cwd=root, capture_output=True)
        if result.returncode or hashlib.sha256(result.stdout).hexdigest() != sha(local_path(root, relative)):
            raise Refusal("reviewed adapter/protocol/tests must be committed before freeze: " + relative)
    return head


def freeze(root=ROOT):
    p = protocol_of(root)
    if p["status"] != "READY_FOR_ROOT_FREEZE" or p["admission"]["blockers"] or not p["admission"].get("root_reviewed"):
        raise Refusal("draft admission is unresolved; ROOT must finalize the protocol before freeze")
    queries = load_queries(root, p)
    calibration_of(root, p, queries)
    exposure_of(root, p, queries)
    acquisition_recipe_files(root, p)
    head = reviewed_commit(root)
    raw = local_path(root, str(RAW))
    if raw.exists() and any(raw.rglob("*")):
        raise Refusal("external test files exist before the freeze")
    library = FrozenLibrary(root)
    pins = runtime_pins(root, p, library)
    out = local_path(root, str(OUT))
    write_new(out / "FREEZE_STARTED.json", {"run_id": RUN_ID, "utc": now()})
    record = {"run_id": RUN_ID, "utc": now(), "pins": pins, "reviewed_commit": head,
              "queries": [query_record(q) for q in queries], "queries_sha256": digest([query_record(q) for q in queries]),
              "external_truth_absent": True, "source_label": None,
              "old_experiment": "closed; no old query proofs or artifacts changed"}
    write_new(out / "FROZEN.json", record)
    return {"run_id": RUN_ID, "frozen_sha256": sha(out / "FROZEN.json"), "queries": len(queries)}


def verify_freeze(root, library):
    p = protocol_of(root)
    out = local_path(root, str(OUT))
    frozen = read(out / "FROZEN.json")
    queries = load_queries(root, p)
    calibration_of(root, p, queries)
    exposure_of(root, p, queries)
    if frozen["run_id"] != RUN_ID or frozen["pins"] != runtime_pins(root, p, library):
        raise Refusal("method/protocol/input/library changed after freeze")
    if frozen["queries"] != [query_record(q) for q in queries] or frozen["queries_sha256"] != digest(frozen["queries"]):
        raise Refusal("query cohort changed after freeze")
    return frozen, queries


def predict(root=ROOT):
    library = FrozenLibrary(root)
    frozen, queries = verify_freeze(root, library)
    out, raw = local_path(root, str(OUT)), local_path(root, str(RAW))
    if raw.exists() and any(raw.rglob("*")):
        raise Refusal("external truth must be absent while forecasting")
    write_new(out / "PREDICT_STARTED.json", {"run_id": RUN_ID, "utc": now(), "freeze_sha256": sha(out / "FROZEN.json")})
    rows = []
    for query in queries:
        for method, result in predict_one(query, library).items():
            relative = f"predictions/{query.neuron}/{method}.npz"
            target = out / relative
            if target.exists():
                raise Refusal("prediction output already exists")
            library.E.G.save_arbor(target, result["arbor"])
            rows.append({"query": query.neuron, "animal": query.animal, "method": method,
                         "path": relative, "sha256": sha(target), "arrays_sha256": library.E.arbor_sha(result["arbor"]),
                         **{k: v for k, v in result.items() if k != "arbor"}})
    manifest = {"run_id": RUN_ID, "freeze_sha256": sha(out / "FROZEN.json"), "queries_sha256": frozen["queries_sha256"],
                "methods": list(METHODS), "predictions": rows}
    write_new(out / "prediction_manifest.json", manifest)
    verify_freeze(root, library)               # detects concurrent source changes before commitment
    if raw.exists() and any(raw.rglob("*")):
        raise Refusal("external truth arrived during forecasting; prediction commitment refused")
    commit = {"run_id": RUN_ID, "utc": now(), "manifest_sha256": sha(out / "prediction_manifest.json"),
              "freeze_sha256": sha(out / "FROZEN.json"), "prediction_count": len(rows), "truth_absent": True}
    write_new(out / "PREDICTIONS_COMMITTED.json", commit)
    return commit


def verify_predictions(root, library, queries):
    out = local_path(root, str(OUT))
    commit, manifest = read(out / "PREDICTIONS_COMMITTED.json"), read(out / "prediction_manifest.json")
    if commit["run_id"] != RUN_ID or commit["manifest_sha256"] != sha(out / "prediction_manifest.json"):
        raise Refusal("prediction commitment differs")
    freeze_sha = sha(out / "FROZEN.json")
    if commit["freeze_sha256"] != freeze_sha or manifest["freeze_sha256"] != freeze_sha:
        raise Refusal("commitment is not bound to this freeze")
    expected = {(q.neuron, m) for q in queries for m in METHODS}
    animal_of = {q.neuron: q.animal for q in queries}
    if tuple(manifest["methods"]) != METHODS or manifest["queries_sha256"] != read(out / "FROZEN.json")["queries_sha256"]:
        raise Refusal("method family or query binding differs")
    seen = set()
    for row in manifest["predictions"]:
        key = (row["query"], row["method"])
        relative = f"predictions/{key[0]}/{key[1]}.npz"
        if key not in expected or key in seen or row["path"] != relative or row["animal"] != animal_of.get(key[0]):
            raise Refusal("duplicate, missing, unknown or noncanonical prediction")
        seen.add(key)
        path = local_path(out, relative)
        if sha(path) != row["sha256"] or library.E.arbor_sha(library.E.G.load_arbor(path)) != row["arrays_sha256"]:
            raise Refusal("prediction file changed after commitment")
    if seen != expected or commit["prediction_count"] != len(expected):
        raise Refusal("not every declared prediction is sealed")
    return commit, manifest


def acquisition_of(root, commit, queries):
    """Read receipts/hashes only. Geometry is opened later, after SCORE_STARTED is exclusive."""
    out = local_path(root, str(OUT))
    path = out / "ACQUISITION.json"
    receipt = read(path)
    if receipt["run_id"] != RUN_ID or receipt["prediction_commitment_sha256"] != sha(out / "PREDICTIONS_COMMITTED.json"):
        raise Refusal("acquisition was not authorized by this prediction commitment")
    if datetime.fromisoformat(receipt["started_utc"]) <= datetime.fromisoformat(commit["utc"]):
        raise Refusal("external acquisition began before prediction commitment")
    expected = {q.neuron for q in queries}
    if set(receipt["files"]) != expected:
        raise Refusal("acquisition receipt must cover every declared query")
    for neuron, item in receipt["files"].items():
        if item["status"] == "unavailable_at_source":
            if not item.get("reason") or not item.get("source_url"):
                raise Refusal("source-unavailable exclusion needs an explicit source reason")
            continue
        if item["status"] != "downloaded" or item["path"] != str(RAW / f"{neuron}.swc"):
            raise Refusal("unrecognized or noncanonical acquired file")
        if item.get("frame") != "CCFv3" or item.get("trace_domain") != "separate_axon_tracing":
            raise Refusal("sample-space/combined trace cannot masquerade as native CCFv3 axon truth")
        if not item.get("source_url") or not item.get("license"):
            raise Refusal("source route and per-neuron reuse terms must be recorded")
        file = local_path(root, item["path"])
        if sha(file) != item["sha256"] or file.stat().st_size != item["bytes"]:
            raise Refusal("acquired file pin differs")
    return receipt, sha(path)


def ingest_truth(path, query, calibration, E):
    """Only supported AFTER sealed prediction: a calibrated separate CCFv3 axon SWC.

    Whole-tree axon geometry is retained, including connectors and all out-of-box
    cable. Source 0/2/5/6 codes are interpreted only within the separately exported
    axon domain confirmed by calibration; no dendrite tracing is read.
    """
    try:
        rows = E.read_swc(path)
        structure = E.tree_structure(rows)
        if set(structure["types"].tolist()) - set(calibration["allowed_node_types"]):
            raise Refusal("unvalidated node type in separate axon tracing")
        rows = rows.copy()
        rows[:, 2:5] = rows[:, [4, 3, 2]]      # native ML DV AP -> AP DV ML
        if query.hemisphere == "Right":
            rows[:, 4] = 11400.0 - rows[:, 4]
        root = structure["root"]
        delta = float(np.linalg.norm(rows[root, 2:5] - np.asarray(query.soma)))
        if delta > 5.0:
            raise Refusal(f"metadata/trace soma mismatch: {delta:.6f} um > 5 um")
        rows[:, 1] = 2
        rows[root, 1] = 1
        xyz, parent, types, stats = E.extract_axon(rows, E.tree_structure(rows), reflect=False, mode="plan")
        if stats["axon_cable_um"] < 1000.0:
            raise Refusal("complete axon cable is below 1000 um")
        arbor = E.canonical_arbor(xyz, parent, types)
        return arbor, {"status": "included", "axon_cable_um": stats["axon_cable_um"],
                       "soma_error_um": delta, "outside_atlas_box": E.outside_box_cable(xyz, parent),
                       "arrays_sha256": E.arbor_sha(arbor), "clipping": "none"}
    except (E.Exclusion, Refusal) as error:
        return None, {"status": "excluded", "reason": str(error)}


def equal_animal_summary(rows):
    by = {}
    for row in rows:
        by.setdefault((row["animal"], row["method"]), []).append(row)
    animal = [{"animal": a, "method": m, "neurons": len(v),
               **{s: float(np.mean([r[s] for r in v])) for s in ("distance_um", "f1_200")}}
              for (a, m), v in sorted(by.items())]
    summary = {}
    for method in METHODS:
        subset = [r for r in animal if r["method"] == method]
        if subset:
            summary[method] = {"animals": len(subset),
                               **{s: float(np.mean([r[s] for r in subset])) for s in ("distance_um", "f1_200")}}
    return animal, summary


def cohort_accounting(queries, admitted_ids):
    admitted_ids = set(admitted_ids)
    declared = {q.neuron for q in queries}
    if not admitted_ids <= declared:
        raise Refusal("an admitted neuron is outside the declaration")
    declared_animals = sorted({q.animal for q in queries})
    admitted_animals = sorted({q.animal for q in queries if q.neuron in admitted_ids})
    lost_animals = sorted(set(declared_animals) - set(admitted_animals))
    return {"declared_neurons": len(queries), "admitted_neurons": len(admitted_ids),
            "excluded_neurons": len(declared - admitted_ids),
            "declared_animals": len(declared_animals), "admitted_animals": len(admitted_animals),
            "lost_animals": len(lost_animals), "declared_animal_ids": declared_animals,
            "admitted_animal_ids": admitted_animals, "lost_animal_ids": lost_animals,
            "cohort_smaller_than_declared": admitted_ids != declared}


def score(root=ROOT):
    library = FrozenLibrary(root)
    frozen, queries = verify_freeze(root, library)
    commit, manifest = verify_predictions(root, library, queries)
    receipt, receipt_sha = acquisition_of(root, commit, queries)
    p = protocol_of(root)
    calibration = calibration_of(root, p, queries)
    out = local_path(root, str(OUT))
    write_new(out / "SCORE_STARTED.json", {"run_id": RUN_ID, "utc": now(), "acquisition_sha256": receipt_sha})
    admission, truth = {}, {}
    # Admission is truth-only; no prediction geometry/metrics is consulted for exclusions.
    for query in queries:
        item = receipt["files"][query.neuron]
        if item["status"] == "unavailable_at_source":
            admission[query.neuron] = {"status": "excluded", "reason": item["reason"], "source_unavailable": True}
        else:
            arbor, admission[query.neuron] = ingest_truth(local_path(root, item["path"]), query, calibration, library.E)
            if arbor is not None:
                truth[query.neuron] = arbor
    write_new(out / "ADMISSION.json", {"run_id": RUN_ID, "queries": admission, "prediction_metrics_read": False})
    if not truth:
        raise Refusal("no admissible truth; score attempt spent and no accuracy computed")
    rows = []
    for row in manifest["predictions"]:
        if row["query"] in truth:
            metric = library.E.G.geometry_metrics(library.E.G.load_arbor(out / row["path"]), truth[row["query"]])
            rows.append({"query": row["query"], "animal": row["animal"], "method": row["method"],
                         "distance_um": metric["distance_um"], "f1_200": metric["f1_200"]})
    animals, summary = equal_animal_summary(rows)
    verdict = {candidate: all(summary[candidate]["distance_um"] < summary[c]["distance_um"] and
                              summary[candidate]["f1_200"] > summary[c]["f1_200"] for c in CONTROLS)
               for candidate in CANDIDATES}
    result = {"run_id": RUN_ID, "summary": summary, "per_animal": animals, "per_neuron": rows,
              "beats_every_fixed_control_on_both_means": verdict,
              "interpretation": "descriptive external transfer; no population-level confirmation or viewer promotion",
              **cohort_accounting(queries, truth), "acquisition_sha256": receipt_sha}
    # Fail if a sealed input, calibration, TRAIN asset or source changed while scoring.
    verify_freeze(root, library)
    verify_predictions(root, library, queries)
    if sha(out / "ACQUISITION.json") != receipt_sha:
        raise Refusal("acquisition receipt changed during scoring")
    for item in receipt["files"].values():
        if item["status"] == "downloaded" and sha(local_path(root, item["path"])) != item["sha256"]:
            raise Refusal("truth file changed during scoring")
    write_new(out / "metrics.json", result)
    write_new(out / "SCORE_COMPLETE.json", {"run_id": RUN_ID, "utc": now(), "metrics_sha256": sha(out / "metrics.json"),
                                           "freeze_sha256": sha(out / "FROZEN.json"), "scored_once": True})
    return {k: v for k, v in result.items() if k not in ("per_animal", "per_neuron")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("audit", "freeze", "predict", "score"))
    args = parser.parse_args()
    try:
        print(canonical(globals()[args.command]()))
    except (Refusal, FileExistsError) as error:
        print(canonical({"run_id": RUN_ID, "refused": str(error)}))
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
