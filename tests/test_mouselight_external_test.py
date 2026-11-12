"""Synthetic geometry and metadata-only checks; no MouseLight trace is read/downloaded.

Numerical parity uses the frozen engine on artificial donor trees and its already
available TRAIN-only regional state. No model fitting, external forecast or score.
"""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("mouselight_external_test_adapter", ROOT / "analysis/mouselight_external_test.py")
M = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = M
spec.loader.exec_module(M)


@pytest.fixture(scope="module")
def E():
    return M.frozen_module()


@pytest.fixture(scope="module")
def toy(E):
    """Artificial trees only; real frozen parameter state is permitted TRAIN evidence."""
    model = M.read(ROOT / M.LIBRARY / "states/all_training.json")
    proof = M.read(ROOT / M.LIBRARY / "states/all_training.proof.json")
    tags = model["field_metadata"][E.R.TAG]["categories"]
    rng = np.random.default_rng(81)
    records = []
    for i in range(18):
        animal = f"toy{i // 2:02d}"
        neuron = f"{animal}_{i:03d}"
        soma = np.array([6600, 5800, 4800], float) + rng.normal(0, 230, 3)
        # Simple rooted branch with two distal paths; resampled by frozen code.
        xyz = np.array([soma, soma + [600, 10 * i, 0], soma + [1800, 30 * i, 90],
                        soma + [700, 650 + 9 * i, 120], soma + [800, 1300 + 11 * i, 200]])
        arbor = E.G.resample(E.G.Arbor(xyz, np.array([-1, 0, 1, 1, 3]), np.array([1, 2, 2, 2, 2])), 50)
        records.append({"neuron": neuron, "animal": animal, "full_id": neuron + "_tag_HY",
                        "tag": tags[i % len(tags)], "source_label": tags[i % len(tags)], "soma_label": "HY",
                        "table_soma": soma, "swc_soma": soma.copy(), "profile": rng.uniform(0.2, 4.0, 49), "arbor": arbor})
    records.sort(key=lambda d: d["neuron"])
    D, F = E.T.pair_costs(records)

    class Library:
        def __init__(self):
            self.E, self.model = E, model
            self.records = records
            self.ids = tuple(d["neuron"] for d in records)
            self.somas = np.array([d["table_soma"] for d in records])
            self.animals = tuple(d["animal"] for d in records)
            self.training_animals = frozenset(self.animals)

        def donor(self, index):
            return self.records[index]

        def distances(self, indices):
            return D[np.ix_(list(indices), list(indices))]

        def medoid_index(self):
            return int(np.argmin(D @ M.animal_prior(self.animals)))

    library = Library()
    return library, proof, D, F


def test_coordinate_order_and_whole_hemisphere_mirror():
    assert M.normalize_soma({"x": 4200, "y": 5900, "z": 6700}) == ((6700., 5900., 4200.), "Left")
    assert M.normalize_soma({"x": 7200, "y": 5900, "z": 6700}) == ((6700., 5900., 4200.), "Right")
    with pytest.raises(M.Refusal, match="bounds"):
        M.normalize_soma({"x": 12000, "y": 5900, "z": 6700})
    with pytest.raises(M.Refusal, match="only"):
        M.normalize_soma({"x": 4200, "y": 5900, "z": 6700, "axon": 1})


def test_query_interface_has_no_hidden_outcome_or_label_slot():
    assert set(M.Query.__dataclass_fields__) == {"neuron", "animal", "soma", "hemisphere"}
    with pytest.raises(TypeError):
        M.Query("AA0001", "mouselight:sample", (6600., 5900., 4500.), "Left", profile=np.ones(49))


def test_metadata_only_audit_is_read_only_and_reports_live_status():
    result = M.audit()
    assert (result["queries"], result["animals"]) == (97, 22)
    assert result["geometry_opened"] is False
    protocol = M.read(ROOT / M.PROTOCOL)
    assert result["blockers"] == protocol["admission"]["blockers"]
    assert result["status"] == protocol["status"]


def test_draft_refuses_before_library_forecast_or_geometry(monkeypatch):
    protocol = M.read(ROOT / M.PROTOCOL)
    protocol["status"] = "DRAFT_ADMISSION_BLOCKED"
    monkeypatch.setattr(M, "protocol_of", lambda root: protocol)
    def forbidden(*args, **kwargs):
        raise AssertionError("must not construct a library/forecast")
    monkeypatch.setattr(M, "FrozenLibrary", forbidden)
    with pytest.raises(M.Refusal, match="draft admission"):
        M.freeze()


def test_metadata_derived_paths_and_hidden_descriptors_do_not_enter_query(tmp_path):
    protocol = M.read(ROOT / M.PROTOCOL)
    data = M.read(ROOT / M.METADATA)
    target = tmp_path / M.METADATA
    target.parent.mkdir(parents=True)
    target.write_text(M.canonical(data))
    protocol["metadata"]["sha256"] = M.sha(target)
    first = M.load_queries(tmp_path, protocol)
    for row in data["records"]:
        row["image_metadata_matches"] = [{"consensus": "hidden_other_tree.swc", "axon_length": 9999999}]
        row["hidden_profile"] = [99] * 49
        row["region"] = "arbitrary_unused_label"
    target.write_text(M.canonical(data))
    with pytest.raises(M.Refusal, match="source pin"):
        M.load_queries(tmp_path, protocol)
    # Re-authorized metadata on a synthetic copy still cannot become predictor inputs.
    protocol["metadata"]["sha256"] = M.sha(target)
    assert first == M.load_queries(tmp_path, protocol)


def test_numeric_parity_with_unchanged_engine_forecast_on_artificial_donors(toy, E):
    library, original_proof, D, F = toy
    q = M.Query("external_001", "external", (6650., 5875., 4750.), "Left")
    portable = M.predict_one(q, library)
    # A throwaway in-memory synthetic served-ID proof only; no old proof file is changed.
    proof = deepcopy(original_proof)
    full_id = "external_001_unknown_HY"
    proof["test_ids"], proof["test_animals"] = [full_id], [q.animal]
    regional = E.Geo.RegionalState(library.model, proof)
    index = {n: i for i, n in enumerate(library.ids)}

    class Cache:
        def sub(self, ids):
            selected = [index[n] for n in ids]
            return D[np.ix_(selected, selected)], F[np.ix_(selected, selected)]

    engine = E.Engine.__new__(E.Engine)
    engine._cache = Cache()
    engine.memo = E.FrozenMemo(lambda: engine._cache)
    for d in library.records:
        engine.memo.register(d["neuron"], d["arbor"])
    row = E.NeuronRow(q.neuron, q.animal, full_id, None, None, "HY", "Left", "Axon_and_dendrite", "test", "HY", q.soma)
    engine.cohort = SimpleNamespace(rows={q.neuron: row}, heldout_animals={q.animal})
    engine.states = SimpleNamespace(state_id_for=lambda n: "fixture", get=lambda n: regional,
                                   animal_sets=lambda s: (set(proof["train_animals"]), set(proof["train_animals"])))
    scope = {"ids": list(library.ids), "soma": library.somas, "animal": list(library.animals),
             "label": np.array([d["source_label"] for d in library.records]), "medoids": {},
             "animal_set": set(library.animals), "animals_sha256": E.ids_sha(sorted(set(library.animals))),
             "ids_sha256": E.ids_sha(library.ids)}
    engine.scope = lambda excluded: scope
    engine.pool_donors = lambda s: library.records
    provenance = {"fold": None, "state_sha256": None, "used": False, "reason_code": "fixture_no_label"}
    engine.layer_answer = lambda r: SimpleNamespace(profile=None, reason="fixture_no_label", upstream_groups=(), provenance=provenance)
    reference, _ = engine.forecast(q.neuron)
    reference = {r["model"]: r for r in reference}
    rename = {"nearest_copy": "nearest", "global_medoid_copy": "medoid"}
    for name in M.METHODS:
        expected = reference[rename.get(name, name)]
        got = portable[name]
        assert got["selected_donor"] == expected["selected_donor"]
        assert E.same_arbor(got["arbor"], expected["arbor"]), name


def test_donor_selection_cannot_see_hidden_truth(toy, tmp_path):
    library, _, _, _ = toy
    q = M.Query("AA9999", "mouselight:never-fitted", (6650., 5875., 4750.), "Left")
    hidden = tmp_path / "AA9999.swc"
    hidden.write_text("malformed outcome canary; not a source of query inputs")
    a = M.predict_one(q, library)
    hidden.write_text("a completely different hidden outcome and regional profile")
    b = M.predict_one(q, library)
    for name in M.METHODS:
        assert a[name]["selected_donor"] == b[name]["selected_donor"]
        assert library.E.same_arbor(a[name]["arbor"], b[name]["arbor"])


def test_missing_source_label_is_none_not_an_invented_category(toy, monkeypatch):
    library, _, _, _ = toy
    original = library.E.Geo.state_posterior
    calls = []

    def record(model, coordinates, tags, **kwargs):
        calls.append((tags, set(kwargs)))
        return original(model, coordinates, tags, **kwargs)
    monkeypatch.setattr(library.E.Geo, "state_posterior", record)
    M.predict_one(M.Query("AA9999", "mouselight:unseen", (6650., 5875., 4750.), "Left"), library)
    assert calls[0] == ([None], set())
    assert "measured_profiles" in calls[1][1]       # measured profiles belong to TRAIN donors only


def test_query_training_animal_overlap_refused(toy):
    library, _, _, _ = toy
    with pytest.raises(M.Refusal, match="in TRAIN"):
        M.predict_one(M.Query("AA9999", "mouselight:" + library.animals[0], (6650., 5875., 4750.), "Left"), library)


def test_whole_animal_aggregation_not_pooled_neurons():
    rows = [{"animal": "a", "method": m, "distance_um": 1., "f1_200": .8}
            for m in M.METHODS for _ in range(10)]
    rows += [{"animal": "b", "method": m, "distance_um": 9., "f1_200": .2} for m in M.METHODS]
    animals, summary = M.equal_animal_summary(rows)
    assert len(animals) == 2 * len(M.METHODS)
    assert summary["nearest_copy"]["distance_um"] == 5.
    assert summary["nearest_copy"]["f1_200"] == .5


def test_cohort_report_identifies_a_completely_lost_animal_and_partial_loss():
    queries = [M.Query("AA0001", "mouselight:a", (6000., 5000., 4000.), "Left"),
               M.Query("AA0002", "mouselight:a", (6000., 5000., 4000.), "Left"),
               M.Query("AA0003", "mouselight:b", (6000., 5000., 4000.), "Left")]
    result = M.cohort_accounting(queries, {"AA0001"})
    assert (result["declared_neurons"], result["admitted_neurons"], result["excluded_neurons"]) == (3, 1, 2)
    assert (result["declared_animals"], result["admitted_animals"], result["lost_animals"]) == (2, 1, 1)
    assert result["declared_animal_ids"] == ["mouselight:a", "mouselight:b"]
    assert result["admitted_animal_ids"] == ["mouselight:a"]
    assert result["lost_animal_ids"] == ["mouselight:b"]
    assert result["cohort_smaller_than_declared"]


def test_exclusive_markers_and_contained_canonical_paths(tmp_path):
    path = tmp_path / "STARTED.json"
    M.write_new(path, {"run": 1})
    with pytest.raises(FileExistsError):
        M.write_new(path, {"run": 2})
    assert M.read(path) == {"run": 1}
    for bad in ("../elsewhere", "/tmp/a", "a/../b", "a//b"):
        with pytest.raises(M.Refusal):
            M.local_path(tmp_path, bad)
    (tmp_path / "linked").symlink_to(tmp_path.parent)
    with pytest.raises(M.Refusal, match="symlink"):
        M.local_path(tmp_path, "linked/escape")


def test_source_byte_changes_are_refused_before_import(tmp_path, monkeypatch):
    file = tmp_path / "analysis/frozen.py"
    file.parent.mkdir()
    file.write_text("original source")
    monkeypatch.setattr(M, "SOURCE_PINS", {"analysis/frozen.py": M.sha(file)})
    M.verify_sources(tmp_path)
    file.write_text("modified source")
    with pytest.raises(M.Refusal, match="source/state pin"):
        M.verify_sources(tmp_path)


def calibration_receipt(tmp_path, protocol, queries):
    evidence = tmp_path / "calibration/native_source_schema.json"
    evidence.parent.mkdir(parents=True)
    evidence.write_text("synthetic calibration evidence only")
    receipt = {"schema": "mouselight_source_calibration_v1", "validated": True, "non_hy": True,
               "frame": "CCFv3", "axes": ["ML", "DV", "AP"], "units": "um", "trace_domain": "separate_axon_tracing",
               "root": "one_type1_root", "nonroot_semantics": "all_nodes_in_separate_axon_tracing_are_axon",
               "animal": "separate-calibration-sample", "neuron": "AA9998", "selection_rule": "synthetic selection",
               "allowed_node_types": [0, 1, 2, 5, 6],
               "evidence": [{"path": "calibration/native_source_schema.json", "sha256": M.sha(evidence)}]}
    target = tmp_path / "calibration/receipt.json"
    target.write_text(M.canonical(receipt))
    protocol["admission"]["calibration"] = {"receipt_path": "calibration/receipt.json", "receipt_sha256": M.sha(target),
                                               "selection_rule": "synthetic selection"}
    return receipt, target


def test_calibration_requires_separate_animal_and_pinned_evidence(tmp_path):
    p = M.read(ROOT / M.PROTOCOL)
    queries = [M.Query("AA9999", "mouselight:test-sample", (6000., 5000., 4000.), "Left")]
    receipt, target = calibration_receipt(tmp_path, p, queries)
    assert M.calibration_of(tmp_path, p, queries)["validated"]
    receipt["animal"] = "test-sample"
    target.write_text(M.canonical(receipt))
    p["admission"]["calibration"]["receipt_sha256"] = M.sha(target)
    with pytest.raises(M.Refusal, match="animal overlaps"):
        M.calibration_of(tmp_path, p, queries)


def test_calibrated_synthetic_truth_keeps_outside_cable_and_marker_connector(tmp_path, E):
    q = M.Query("AA9999", "mouselight:test", (12800., 5000., 4000.), "Left")
    # Native x/y/z=ML/DV/AP; 0 path and 5 branch markers are within the verified separate axon domain.
    path = tmp_path / "synthetic.swc"
    path.write_text("1 1 4000 5000 12800 1 -1\n2 0 4000 5000 13000 1 1\n3 5 4000 5000 13500 1 2\n4 6 4000 5000 14200 1 3\n")
    calibration = {"allowed_node_types": [0, 1, 2, 5, 6]}
    arbor, log = M.ingest_truth(path, q, calibration, E)
    assert log["status"] == "included" and log["clipping"] == "none"
    assert arbor.xyz[:, 0].max() == 14200.
    assert log["axon_cable_um"] == 1400.
    assert log["outside_atlas_box"]["outside_um"] > 0
    assert arbor.cable() == 1400.
    # Soma mismatch and combined-domain dendrite code are mechanical exclusions.
    wrong = M.Query(q.neuron, q.animal, (12820., 5000., 4000.), "Left")
    assert M.ingest_truth(path, wrong, calibration, E)[1]["status"] == "excluded"
    path.write_text(path.read_text().replace("2 0 ", "2 3 "))
    assert M.ingest_truth(path, q, calibration, E)[1]["status"] == "excluded"


def test_acquisition_requires_commitment_and_post_commit_time(tmp_path):
    out = tmp_path / M.OUT
    out.mkdir(parents=True)
    commit = {"utc": "2026-11-12T12:00:00+00:00"}
    (out / "PREDICTIONS_COMMITTED.json").write_text(M.canonical(commit))
    q = M.Query("AA9999", "mouselight:test", (6000., 5000., 4000.), "Left")
    receipt = {"run_id": M.RUN_ID, "prediction_commitment_sha256": M.sha(out / "PREDICTIONS_COMMITTED.json"),
               "started_utc": "2026-11-12T11:59:59+00:00", "files": {q.neuron: {
                   "status": "unavailable_at_source", "reason": "source receipt", "source_url": "https://example.test"}}}
    (out / "ACQUISITION.json").write_text(M.canonical(receipt))
    with pytest.raises(M.Refusal, match="before prediction"):
        M.acquisition_of(tmp_path, commit, [q])
    receipt["started_utc"] = "2026-11-12T12:00:01+00:00"
    (out / "ACQUISITION.json").write_text(M.canonical(receipt))
    assert M.acquisition_of(tmp_path, commit, [q])[0] == receipt


def test_sealed_prediction_bytes_and_query_binding_cannot_change(tmp_path, toy):
    library, _, _, _ = toy
    q = M.Query("AA9999", "mouselight:unseen", (6650., 5875., 4750.), "Left")
    out = tmp_path / M.OUT
    out.mkdir(parents=True)
    query_sha = M.digest([M.query_record(q)])
    (out / "FROZEN.json").write_text(M.canonical({"queries_sha256": query_sha}))
    rows = []
    for method, result in M.predict_one(q, library).items():
        relative = f"predictions/{q.neuron}/{method}.npz"
        library.E.G.save_arbor(out / relative, result["arbor"])
        rows.append({"query": q.neuron, "animal": q.animal, "method": method, "path": relative,
                     "sha256": M.sha(out / relative), "arrays_sha256": library.E.arbor_sha(result["arbor"])})
    manifest = {"run_id": M.RUN_ID, "freeze_sha256": M.sha(out / "FROZEN.json"), "queries_sha256": query_sha,
                "methods": list(M.METHODS), "predictions": rows}
    (out / "prediction_manifest.json").write_text(M.canonical(manifest))
    commit = {"run_id": M.RUN_ID, "manifest_sha256": M.sha(out / "prediction_manifest.json"),
              "freeze_sha256": M.sha(out / "FROZEN.json"), "prediction_count": len(rows)}
    (out / "PREDICTIONS_COMMITTED.json").write_text(M.canonical(commit))
    M.verify_predictions(tmp_path, library, [q])
    target = out / rows[0]["path"]
    target.write_bytes(target.read_bytes() + b"changed")
    with pytest.raises(M.Refusal, match="changed after commitment"):
        M.verify_predictions(tmp_path, library, [q])


@pytest.fixture
def prospective_fixture(tmp_path, toy, monkeypatch):
    """Complete new lifecycle on artificial assets in a temporary root, never real outputs."""
    library, _, D, _ = toy
    p = M.read(ROOT / M.PROTOCOL)
    p["status"] = "READY_FOR_ROOT_FREEZE"
    p["admission"]["root_reviewed"] = True
    p["admission"]["blockers"] = []
    metadata = M.read(ROOT / M.METADATA)
    keep = {metadata["records"][0]["neuron"], metadata["records"][1]["neuron"]}
    p["admission"]["excluded_neurons"] = sorted(r["neuron"] for r in metadata["records"] if r["neuron"] not in keep)
    target = tmp_path / M.METADATA
    target.parent.mkdir(parents=True)
    target.write_bytes((ROOT / M.METADATA).read_bytes())
    queries = M.load_queries(tmp_path, p)
    calibration_receipt(tmp_path, p, queries)
    exposure = {"schema": "mouselight_exposure_audit_v1", "passed": True, "test_geometry_unopened": True,
                "overlapping_query_ids": [], "overlapping_sample_ids": [],
                "query_ids_sha256": M.digest([q.neuron for q in queries]),
                "evidence": [{"path": "calibration/native_source_schema.json",
                              "sha256": M.sha(tmp_path / "calibration/native_source_schema.json")}]}
    exposure_path = tmp_path / "calibration/exposure.json"
    exposure_path.write_text(M.canonical(exposure))
    p["admission"]["exposure_audit"] = {"path": "calibration/exposure.json", "sha256": M.sha(exposure_path)}
    for name in ("analysis/mouselight_external_test.py", "tests/test_mouselight_external_test.py"):
        file = tmp_path / name
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_bytes((ROOT / name).read_bytes())
    # The acquisition code is pinned but never executed by this synthetic adapter fixture.
    for name in (M.ACQUISITION_CODE, M.ACQUISITION_TESTS):
        file = tmp_path / name
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_text("# artificial acquisition source fixture; no network\n")
    p["acquisition"]["root_reviewed"] = True
    p["acquisition"]["root_admission"] = "READY"
    for role in ("source_manifest", "constants"):
        relative = "fixture/" + role + ".json"
        file = tmp_path / relative
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_text(M.canonical({"synthetic_role": role}))
        p["acquisition"][role] = {"path": relative, "sha256": M.sha(file)}
    (tmp_path / M.PROTOCOL).parent.mkdir(parents=True)
    (tmp_path / M.PROTOCOL).write_text(M.canonical(p))
    state = tmp_path / "fixture/train_state.json"
    state.parent.mkdir(exist_ok=True)
    state.write_text("synthetic TRAIN fixture")
    monkeypatch.setattr(M, "SOURCE_PINS", {"fixture/train_state.json": M.sha(state)})
    proxy_E = SimpleNamespace(**library.E.__dict__)
    proxy_E.FROZEN_PINS = {"source_sha256": {}}
    world = SimpleNamespace(manifest=str(state), local_manifest=str(state), soma_table=str(state),
                            regional_table=str(state), trace_manifests=())
    portable = SimpleNamespace(E=proxy_E, model=library.model, ids=library.ids, somas=library.somas,
                               animals=library.animals, training_animals=library.training_animals,
                               donor=library.donor, distances=library.distances, medoid_index=library.medoid_index,
                               engine=SimpleNamespace(world=world),
                               pins=lambda: {"fixture_library": M.digest(D.tolist())})
    monkeypatch.setattr(M, "FrozenLibrary", lambda root: portable)
    monkeypatch.setattr(M, "reviewed_commit", lambda root: "synthetic-reviewed-commit")
    return tmp_path, p, queries, portable


def test_full_synthetic_lifecycle_is_exclusive_and_preserves_every_method(prospective_fixture):
    root, p, queries, library = prospective_fixture
    M.freeze(root)
    with pytest.raises(FileExistsError):
        M.freeze(root)
    commit = M.predict(root)
    assert commit["prediction_count"] == len(queries) * len(M.METHODS)
    with pytest.raises(FileExistsError):
        M.predict(root)
    out = root / M.OUT
    receipt = {"run_id": M.RUN_ID, "prediction_commitment_sha256": M.sha(out / "PREDICTIONS_COMMITTED.json"),
               "started_utc": "2099-01-01T00:00:00+00:00", "files": {}}
    for q in queries:
        ap, dv, ml = q.soma
        if q.hemisphere == "Right":
            ml = 11400.0 - ml
        file = root / M.RAW / f"{q.neuron}.swc"
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_text(f"1 1 {ml} {dv} {ap} 1 -1\n2 2 {ml} {dv} {ap + 600} 1 1\n3 6 {ml} {dv} {ap + 1400} 1 2\n")
        receipt["files"][q.neuron] = {"status": "downloaded", "path": str(file.relative_to(root)),
                                       "sha256": M.sha(file), "bytes": file.stat().st_size,
                                       "source_url": "https://example.test/synthetic.swc", "license": "synthetic fixture",
                                       "frame": "CCFv3", "trace_domain": "separate_axon_tracing"}
    (out / "ACQUISITION.json").write_text(M.canonical(receipt))
    result = M.score(root)
    assert result["admitted_neurons"] == len(queries)
    assert set(result["summary"]) == set(M.METHODS)
    assert M.read(out / "SCORE_COMPLETE.json")["scored_once"] is True
    with pytest.raises(FileExistsError):
        M.score(root)


@pytest.mark.parametrize("changed", ["protocol", "library", "source", "metadata", "acquirer", "acquirer_tests",
                                     "acquisition_manifest", "constants"])
def test_freeze_binds_protocol_library_source_and_metadata(prospective_fixture, changed):
    root, p, queries, library = prospective_fixture
    M.freeze(root)
    if changed == "protocol":
        p["scores"]["inference"] = "modified after freeze"
        (root / M.PROTOCOL).write_text(M.canonical(p))
    elif changed == "library":
        library.pins = lambda: {"fixture_library": "different TRAIN donors"}
    elif changed == "source":
        (root / "fixture/train_state.json").write_text("different state")
    elif changed == "metadata":
        (root / M.METADATA).write_text("different metadata")
    elif changed == "acquirer":
        (root / M.ACQUISITION_CODE).write_text("different acquisition program")
    elif changed == "acquirer_tests":
        (root / M.ACQUISITION_TESTS).write_text("different acquisition tests")
    else:
        role = "source_manifest" if changed == "acquisition_manifest" else "constants"
        (root / p["acquisition"][role]["path"]).write_text("different acquisition input")
    with pytest.raises(M.Refusal):
        M.verify_freeze(root, library)


def test_truth_presence_before_prediction_is_refused(prospective_fixture):
    root, p, queries, library = prospective_fixture
    M.freeze(root)
    file = root / M.RAW / "unopened_canary.swc"
    file.parent.mkdir(parents=True)
    file.write_text("never opened")
    with pytest.raises(M.Refusal, match="truth must be absent"):
        M.predict(root)


def test_truth_arriving_during_prediction_prevents_commitment(prospective_fixture, monkeypatch):
    root, p, queries, library = prospective_fixture
    M.freeze(root)
    original = M.predict_one

    def inject_unopened_canary(query, donors):
        result = original(query, donors)
        file = root / M.RAW / "unopened_canary.swc"
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_text("arrived during synthetic forecasting; never parsed")
        return result

    monkeypatch.setattr(M, "predict_one", inject_unopened_canary)
    with pytest.raises(M.Refusal, match="arrived during forecasting"):
        M.predict(root)
    assert (root / M.OUT / "PREDICT_STARTED.json").exists()
    assert (root / M.OUT / "prediction_manifest.json").exists()
    assert not (root / M.OUT / "PREDICTIONS_COMMITTED.json").exists()
