"""Synthetic acquisition checks only; no live requests or real axon files."""
from copy import deepcopy
from datetime import datetime, timedelta, timezone
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, ROOT / file)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


A = load("mouselight_acquirer_tests", "analysis/mouselight_external_acquire.py")
M = load("mouselight_acquirer_adapter_contract", "analysis/mouselight_external_test.py")


def save(root, name, value):
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(A.canonical(value))
    return path


@pytest.fixture
def sealed(tmp_path):
    records = []
    for i in range(97):
        records.append({"neuron": f"AA{i+1:04d}", "sample_uuid": f"sample{i%22:02d}",
                        "axon_tracing_uuid": f"axon{i:04d}", "tracing_structure": {"name": "axon", "value": 1},
                        "trace_domain": "separate_axon_tracing", "doi": f"10.synthetic/{i}",
                        "license": {"name": "synthetic"}, "license_metadata_sha256": "a" * 64,
                        "ccfv3_soma_um": {"x": 1000.+i, "y": 2000., "z": 3000.}})
    constants = {"data": {"structureIdentifiers": [{"id": "soma", "name": "soma", "value": 1},
                 {"id": "path", "name": "path", "value": 0}, {"id": "branch", "name": "branch point", "value": 5},
                 {"id": "end", "name": "end point", "value": 6}, {"id": "dendrite", "name": "dendrite", "value": 3}]}}
    save(tmp_path, A.SOURCE, {"schema": "mouselight_acquisition_source_v1", "endpoint": A.ENDPOINT, "records": records})
    save(tmp_path, A.CONSTANTS, constants)
    save(tmp_path, "soma_metadata.json", {"synthetic": True})
    for name in (A.CODE, A.TESTS, A.ADAPTER, A.ADAPTER_TESTS):
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((ROOT / name).read_bytes())
    p = {"run_id": A.RUN_ID, "status": "READY_FOR_ROOT_FREEZE", "admission": {"root_reviewed": True, "blockers": []},
         "output_dir": str(A.OUT), "raw_dir": str(A.RAW), "methods": list(A.METHODS),
         "metadata": {"path": "soma_metadata.json"}, "acquisition": {"root_admission": "READY", "root_reviewed": True,
         "code": A.CODE, "tests": A.TESTS, "source_manifest": {"path": A.SOURCE, "sha256": A.sha(tmp_path/A.SOURCE)},
         "constants": {"path": A.CONSTANTS, "sha256": A.sha(tmp_path/A.CONSTANTS)}}}
    save(tmp_path, A.PROTOCOL, p)
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "add", A.CODE, A.TESTS, A.ADAPTER, A.ADAPTER_TESTS, str(A.PROTOCOL)], cwd=tmp_path, check=True)
    subprocess.run(["git", "-c", "user.name=Synthetic", "-c", "user.email=synthetic@example.invalid", "commit", "-qm", "synthetic sealed source"], cwd=tmp_path, check=True)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=tmp_path).decode().strip()
    # Use the actual adapter's record constructor: no parallel mock schema.
    queries = [M.query_record(M.Query(r["neuron"], "mouselight:"+r["sample_uuid"],
               (r["ccfv3_soma_um"]["z"], r["ccfv3_soma_um"]["y"], r["ccfv3_soma_um"]["x"]), "Left")) for r in records]
    names = (A.CODE, A.TESTS, A.ADAPTER, A.ADAPTER_TESTS, str(A.PROTOCOL), A.SOURCE, A.CONSTANTS, "soma_metadata.json")
    frozen = {"run_id": A.RUN_ID, "queries": queries, "queries_sha256": A.digest(queries), "external_truth_absent": True,
              "reviewed_commit": revision, "pins": {"files": {n: A.sha(tmp_path/n) for n in names}}}
    save(tmp_path, A.OUT/"FROZEN.json", frozen)
    rows = []
    for q in queries:
        for method in A.METHODS:
            relative = f"predictions/{q['neuron']}/{method}.npz"
            file = tmp_path / A.OUT / relative
            file.parent.mkdir(parents=True, exist_ok=True)
            file.write_bytes(b"synthetic sealed forecast")
            rows.append({"query": q["neuron"], "animal": q["animal"], "method": method, "path": relative, "sha256": A.sha(file)})
    manifest = {"run_id": A.RUN_ID, "freeze_sha256": A.sha(tmp_path/A.OUT/"FROZEN.json"), "queries_sha256": frozen["queries_sha256"], "methods": list(A.METHODS), "predictions": rows}
    save(tmp_path, A.OUT/"prediction_manifest.json", manifest)
    save(tmp_path, A.OUT/"PREDICTIONS_COMMITTED.json", {"run_id": A.RUN_ID, "utc": (datetime.now(timezone.utc)-timedelta(minutes=1)).isoformat(),
         "freeze_sha256": manifest["freeze_sha256"], "manifest_sha256": A.sha(tmp_path/A.OUT/"prediction_manifest.json"), "prediction_count": len(rows), "truth_absent": True})
    return tmp_path, records, constants


def payload(record):
    xyz = record["ccfv3_soma_um"]
    root = {**xyz, "radius": 1., "sampleNumber": 1, "parentNumber": -1, "structureIdentifierId": "soma"}
    tip = {**root, "x": xyz["x"]+1500, "sampleNumber": 2, "parentNumber": 1, "structureIdentifierId": "end"}
    return {"tracings": [{"id": record["axon_tracing_uuid"], "nodes": [tip, root]}], "timing": {}}


class Response(io.BytesIO):
    def __init__(self, body, headers=None, status=200, url=A.ENDPOINT):
        super().__init__(body)
        self.headers = headers if headers is not None else {"Content-Length": str(len(body))}
        self.status, self.url = status, url

    def geturl(self):
        return self.url


def fetcher(records, seen, fail_at=None, empty_at=None):
    by_uuid = {r["axon_tracing_uuid"]: r for r in records}
    def fetch(request):
        body = json.loads(request.data)
        assert request.full_url == A.ENDPOINT and request.method == "POST"
        assert body["ccfVersion"] == "CCFV30" and len(body["ids"]) == 1
        row = by_uuid[body["ids"][0]]
        seen.append(row["neuron"])
        if fail_at == row["neuron"]:
            raise TimeoutError("synthetic transport timeout")
        obj = {"tracings": []} if empty_at == row["neuron"] else payload(row)
        return Response(A.canonical(obj).encode())
    return fetch


def test_check_is_read_only_and_exact_adapter_schema(sealed, monkeypatch):
    root, _, _ = sealed
    before = sorted(str(p.relative_to(root)) for p in root.rglob("*"))
    monkeypatch.setattr(A, "open_response", lambda *a: pytest.fail("network from check"))
    result = A.gate(root)
    assert len(result["records"]) == 97
    assert before == sorted(str(p.relative_to(root)) for p in root.rglob("*"))
    assert not (root/A.RAW).exists()


@pytest.mark.parametrize("mutation", ["draft", "not_ready", "code", "protocol", "prediction", "manifest", "freeze", "raw", "missing_prediction", "executing_copy"])
def test_gate_refuses_before_network_or_raw_creation(sealed, monkeypatch, mutation):
    root, _, _ = sealed
    if mutation in {"draft", "not_ready", "protocol"}:
        p = A.read(root/A.PROTOCOL)
        if mutation == "draft": p["status"] = "DRAFT"
        elif mutation == "not_ready": p["acquisition"]["root_admission"] = "DRAFT"
        else: p["new_field"] = True
        save(root, A.PROTOCOL, p)
    elif mutation == "code": (root/A.CODE).write_text("changed")
    elif mutation == "executing_copy": monkeypatch.setattr(A, "__file__", root/"alternate.py"); (root/"alternate.py").write_text("different code")
    elif mutation == "raw": (root/A.RAW).mkdir(parents=True)
    else:
        names = {"prediction": A.OUT/"predictions/AA0001/nearest_copy.npz", "manifest": A.OUT/"prediction_manifest.json", "freeze": A.OUT/"FROZEN.json", "missing_prediction": A.OUT/"predictions/AA0001/nearest_copy.npz"}
        file = root/names[mutation]
        if mutation == "missing_prediction": file.unlink()
        else: file.write_bytes(file.read_bytes()+b" ")
    monkeypatch.setattr(A, "open_response", lambda *a: pytest.fail("network before gate"))
    with pytest.raises((A.Refusal, FileNotFoundError)):
        A.acquire(root)
    assert not (root/A.OUT/"ACQUISITION_STARTED.json").exists()
    if mutation != "raw": assert not (root/A.RAW).exists()


def test_committed_guard_rejects_uncommitted_files(sealed):
    root, _, _ = sealed
    frozen = A.read(root/A.OUT/"FROZEN.json")
    (root/A.TESTS).write_text("edited after commit")
    with pytest.raises(A.Refusal, match="uncommitted or changed"):
        A.committed(root, frozen["reviewed_commit"], [A.TESTS])


def test_complete_sequential_acquisition_on_fresh_parent_matches_adapter_receipt(sealed):
    root, records, _ = sealed
    assert not (root/A.RAW).parent.exists()
    seen = []
    result = A.acquire(root, fetch=fetcher(records, seen, empty_at="AA0005"))
    assert seen == [r["neuron"] for r in records]
    assert result["statuses"] == {"downloaded": 96, "unavailable_at_source": 1}
    commit = A.read(root/A.OUT/"PREDICTIONS_COMMITTED.json")
    queries = [M.Query(q["neuron"], q["animal"], tuple(q["soma_ap_dv_ml_um"]), q["hemisphere"]) for q in A.read(root/A.OUT/"FROZEN.json")["queries"]]
    receipt, _ = M.acquisition_of(root, commit, queries)
    assert len(receipt["files"]) == 97
    with pytest.raises(A.Refusal, match="already complete"):
        A.acquire(root, fetch=lambda *a: pytest.fail("repeated acquisition"))


def test_timeout_stops_without_exclusion_then_hash_only_resume(sealed, monkeypatch):
    root, records, _ = sealed
    first = []
    with pytest.raises(A.Refusal, match="no neuron excluded"):
        A.acquire(root, fetch=fetcher(records, first, fail_at="AA0003"))
    assert first == ["AA0001", "AA0002", "AA0003"]
    assert not (root/A.OUT/"ACQUISITION.json").exists()
    assert A.budget_used(root/A.OUT) >= A.RESERVATION
    original = A.convert
    def convert(body, record, constants):
        assert record["neuron"] not in {"AA0001", "AA0002"}
        return original(body, record, constants)
    monkeypatch.setattr(A, "convert", convert)
    resumed = []
    result = A.acquire(root, resume=True, fetch=fetcher(records, resumed))
    assert resumed == [r["neuron"] for r in records[2:]]
    assert result["statuses"] == {"downloaded": 97}


@pytest.mark.parametrize("field", ["raw", "swc", "status", "provenance"])
def test_resume_refuses_modified_checkpoint_or_bytes(sealed, field):
    root, records, _ = sealed
    with pytest.raises(A.Refusal): A.acquire(root, fetch=fetcher(records, [], fail_at="AA0002"))
    checkpoint = root/A.OUT/"acquisition_checkpoints/AA0001.json"
    item = A.read(checkpoint)
    if field in {"raw", "swc"}:
        name = item["raw_json" if field == "raw" else "swc"]["path"]
        (root/name).write_bytes(b"changed")
    else:
        item["status" if field == "status" else "sample_uuid"] = "unavailable_at_source" if field == "status" else "other"
        checkpoint.write_text(A.canonical(item))
    with pytest.raises(A.Refusal): A.acquire(root, resume=True, fetch=lambda *a: pytest.fail("network with tampered checkpoint"))


def test_converter_preserves_values_and_generic_axon_roles(sealed):
    _, records, constants = sealed
    obj = payload(records[0])
    result = A.convert(A.canonical(obj).encode(), records[0], constants)
    rows = [s.split() for s in result.decode().splitlines() if not s.startswith("#")]
    assert rows[0] == ["1", "1", "1000", "2000", "3000", "1", "-1"]
    assert rows[1] == ["2", "6", "2500", "2000", "3000", "1", "1"]
    assert A.convert(b'{"tracings":[]}', records[0], constants) is None


@pytest.mark.parametrize("mutation", ["mixed", "wrong_uuid", "dendrite", "cycle", "second_root", "nonfinite", "wrong_soma", "bool_id"])
def test_converter_refuses_unreviewed_source_semantics(sealed, mutation):
    _, records, constants = sealed
    obj = payload(records[0]); nodes = obj["tracings"][0]["nodes"]
    if mutation == "mixed": obj["tracings"].append(deepcopy(obj["tracings"][0]))
    elif mutation == "wrong_uuid": obj["tracings"][0]["id"] = "different"
    elif mutation == "dendrite": nodes[0]["structureIdentifierId"] = "dendrite"
    elif mutation == "cycle": nodes[0]["parentNumber"] = 2
    elif mutation == "second_root": nodes[0]["parentNumber"] = -1
    elif mutation == "nonfinite": nodes[0]["x"] = float("nan")
    elif mutation == "wrong_soma": nodes[1]["x"] += 10
    elif mutation == "bool_id": nodes[0]["sampleNumber"] = True
    body = json.dumps(obj).encode()
    with pytest.raises(A.Refusal): A.convert(body, records[0], constants)


@pytest.mark.parametrize("kind", ["oversize", "false_length", "encoded", "redirect", "http404"])
def test_transport_faults_stop_without_source_exclusion(sealed, monkeypatch, kind):
    root, records, _ = sealed
    def fetch(request):
        body = A.canonical(payload(records[0])).encode()
        if kind == "oversize": return Response(body, {"Content-Length": str(A.MAX_RESPONSE_BYTES+1)})
        if kind == "false_length": return Response(body, {"Content-Length": str(len(body)+2)})
        if kind == "encoded": return Response(body, {"Content-Encoding": "gzip"})
        if kind == "redirect": return Response(body, url="https://other.invalid/tracings")
        return Response(body, status=404)
    with pytest.raises(A.Refusal, match="no neuron excluded"): A.acquire(root, fetch=fetch)
    assert not (root/A.OUT/"ACQUISITION.json").exists()
    assert not list((root/A.OUT/"acquisition_checkpoints").glob("*.json"))


def test_reservation_budget_includes_unfinished_requests(sealed, monkeypatch):
    root, records, _ = sealed
    monkeypatch.setattr(A, "MAX_TOTAL_BYTES", A.RESERVATION)
    with pytest.raises(A.Refusal): A.acquire(root, fetch=fetcher(records, [], fail_at="AA0001"))
    attempts = root/A.OUT/"acquisition_attempts"
    end = next(attempts.glob("*.end.json"))
    end.unlink()                            # simulate a process death before terminal receipt
    assert A.budget_used(root/A.OUT) == A.RESERVATION
    with pytest.raises(A.Refusal, match="budget"):
        A.acquire(root, resume=True, fetch=lambda *a: pytest.fail("request exceeding durable budget"))


def test_canonical_paths_and_exclusive_running_lock(sealed):
    root, records, _ = sealed
    with pytest.raises(A.Refusal): A.local(root, "../elsewhere")
    (root/"alias").symlink_to(root/A.OUT, target_is_directory=True)
    with pytest.raises(A.Refusal): A.local(root, "alias/file")
    (root/A.OUT/"ACQUISITION_RUNNING.lock").write_text("another worker")
    with pytest.raises(FileExistsError): A.acquire(root, fetch=lambda *a: pytest.fail("concurrent request"))
    assert not (root/A.RAW).exists()
