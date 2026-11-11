"""Synthetic files only: metadata gates precede every TEST raw byte read."""
import copy
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys

import pytest

P = Path(__file__).resolve().parents[1] / "analysis/axon_trace_download_stage_b_audit.py"
spec = importlib.util.spec_from_file_location("b_audit", P)
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + "\n")


def tsv(rows, fields):
    f = io.StringIO(newline="")
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return f.getvalue().encode()


@pytest.fixture
def bundle(tmp_path, monkeypatch):
    root = tmp_path
    base = root / "data/raw/jiao2025/traces"
    b = base / "stage_b"
    b.mkdir(parents=True)
    plan = [{"neuron": f"{100001+i:06d}_001", "animal": str(100001+i), "part": p}
            for i, p in enumerate(("train", "train", "validation", "test", "test"))]
    raw_plan = tsv(plan, ("neuron", "animal", "part"))
    plan_path = root / "reports/axon_dataset_manifest.tsv"
    plan_path.parent.mkdir(parents=True)
    plan_path.write_bytes(raw_plan)
    local = root / "data/interim/axon_dataset/local_manifest.tsv"
    local.parent.mkdir(parents=True)
    local.write_bytes(raw_plan)
    plan_sha, local_sha = A.sha_file(plan_path), A.sha_file(local)
    monkeypatch.setattr(A, "PLAN_SHA", plan_sha)
    monkeypatch.setattr(A, "COUNTS", {"train": 2, "validation": 1, "test": 2})
    monkeypatch.setattr(A, "ANIMALS", {"train": 2, "validation": 1, "test": 2})
    source_pins = {}
    for rel in ("analysis/axon_trace_download.py", "analysis/axon_trace_download_stage_b.py"):
        p = root / rel
        p.parent.mkdir(exist_ok=True)
        p.write_text("# synthetic source\n")
        source_pins[rel] = A.sha_file(p)
    a_rows = []
    a_receipts = {}
    for r in plan[:3]:
        n = r["neuron"]
        row = dict(file=f"{r['part']}/{r['animal']}/{n}.swc", neuron=n, bytes=3,
                   sha256=hashlib.sha256(b"raw").hexdigest(), source="legacy", frame="portal_true_hemisphere", licence="private", how="copy")
        a_rows.append(row)
        path = base / "receipts" / (n + ".json")
        put(path, row)
        a_receipts[n] = A.sha_file(path)
    prefix = tsv(a_rows, A.FIELDS)
    (b / "stage_a_manifest.tsv").write_bytes(prefix)
    put(base / "summary.json", {"stage": "A", "status": "complete"})
    (base / "acquisition.jsonl").write_text("{}\n")
    a_metadata = {p: A.sha_file(base / p) for p in ("summary.json", "acquisition.jsonl")}
    wanted = [r["neuron"] for r in plan[3:]]
    test_sha = hashlib.sha256("\n".join(wanted).encode()).hexdigest()
    engine = root / "data/interim/axon_engine/real"
    paths = {"freeze": "freeze/frozen.json", "commitment": "test/PREDICTION_COMMITMENT.json",
             "prediction_manifest": "test/prediction_manifest.json", "copies_verified": "test/COPIES_VERIFIED.json"}
    protocol = root / "reports/axon_engine_protocol.md"
    protocol.write_text("STATUS: FINAL\n")
    pins = {"training_traces": {"listing_sha256": "train"}, "cohort": {"manifest": plan_sha, "local_manifest": local_sha},
            "code_sha256": {"engine": {}, "frozen": {}, "general_layer": {}}, "protocol_sha256": A.sha_file(protocol)}
    freeze = {"schema": "axon_engine_freeze_v2", "world_kind": "real", "validation_method_equals_frozen_method": True,
              "git": {"engine_files_committed_and_clean": True},
              "test_access": {"strongest_seal_no_test_trace_on_disk": True, "test_traces_listed_on_disk_at_freeze": 0},
              "pins": dict(pins, family="all", test_queries=2, test_queries_sha256=test_sha), "family": "all",
              "selected": "E2_local_kernel", "descriptive_reference": "E2_local_kernel"}
    put(engine / paths["freeze"], freeze)
    predictions = []
    for n in wanted:
        for m in A.METHODS:
            rel = f"predictions/{n}/{m}.npz"
            p = engine / "test" / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes((n + m).encode())
            predictions.append({"query_id": n, "model": m, "path": rel, "sha256": A.sha_file(p)})
    manifest = {"schema": "axon_engine_prediction_manifest_v2", "stage": "test", "part": "test", "world_kind": "real",
                "queries": wanted, "family": "all", "pins": pins, "predictions": predictions,
                "frozen_selection": {"freeze_sha256": A.sha_file(engine / paths["freeze"]),
                                     "selected": freeze["selected"], "descriptive_reference": freeze["descriptive_reference"]}}
    put(engine / paths["prediction_manifest"], manifest)
    manifest_sha = A.sha_file(engine / paths["prediction_manifest"])
    put(engine / paths["commitment"], {"schema": "axon_engine_prediction_commitment_v2", "stage": "test", "manifest_sha256": manifest_sha,
                                       "queries": 2, "predictions": 2 * len(A.METHODS), "methods": list(A.METHODS)})
    put(engine / paths["copies_verified"], {"stage": "test", "manifest_sha256": manifest_sha, "all_equal": True,
                                            "differences": 0, "queries": 2, "forecasts_compared": 8})
    auth = {"schema": "synarch_stage_b_authorization_v1", "proof_paths": {k: str((engine / p).relative_to(root)) for k, p in paths.items()},
            "proof_sha256": {k: A.sha_file(engine / p) for k, p in paths.items()}, "plan_sha256": plan_sha, "local_manifest_sha256": local_sha,
            "stage_a_manifest_sha256": hashlib.sha256(prefix).hexdigest(), "stage_a_manifest_bytes": len(prefix),
            "stage_a_receipts_sha256": a_receipts, "stage_a_metadata_sha256": a_metadata, "source_sha256": source_pins}
    put(b / "authorization.json", auth)
    auth_sha = A.sha_file(b / "authorization.json")
    receipts, events = {}, [{"run_id": "r1", "event": "run_start"}]
    for i, r in enumerate(plan[3:]):
        n = r["neuron"]
        available = i == 0
        rel = f"test/{r['animal']}/{n}.swc" if available else ""
        data = b"# synthetic SWC bytes only\n"
        if available:
            p = base / rel
            p.parent.mkdir(parents=True)
            p.write_bytes(data)
        receipt = dict(schema="synarch_stage_b_download_v1", stage="B", plan_sha256=plan_sha, local_manifest_sha256=local_sha,
                       authorization_sha256=auth_sha, proof_sha256=auth["proof_sha256"], neuron=n, file=rel,
                       bytes=len(data) if available else 0, sha256=hashlib.sha256(data).hexdigest() if available else A.UNAVAILABLE,
                       source=A.PORTAL.format(n), url=A.PORTAL.format(n), final_url=A.PORTAL.format(n),
                       frame="portal_true_hemisphere" if available else A.UNAVAILABLE, status=200 if available else 404,
                       availability="available" if available else A.UNAVAILABLE, licence="private", how="public source 200" if available else "HTTP 404 source unavailable")
        put(b / "receipts" / (n + ".json"), receipt)
        receipts[n] = receipt
        events.extend([{"run_id": "r1", "event": "request_start", "neuron": n, "url": A.PORTAL.format(n), "authorization_sha256": auth_sha, "monotonic_start": 10 + i * 2, "attempt": 1},
                       {"run_id": "r1", "event": "checkpoint", "neuron": n, "bytes": receipt["bytes"], "sha256": receipt["sha256"]}])
    body = tsv([{k: r[k] for k in A.FIELDS} for r in receipts.values()], A.FIELDS).split(b"\n", 1)[1]
    (base / "MANIFEST.tsv").write_bytes(prefix + body)
    run = {"run_id": "r1", "stage": "B", "authorization_sha256": auth_sha, "completed": 2, "unfinished_ids": [],
           "status": "complete", "failures": {}, "acquired": 1, "unavailable_at_source": 1, "utc": "2026-11-11T20:00:00+00:00",
           "unavailable_ids": [wanted[1]], "unavailable_by_animal": {"100004": 0, "100005": 1},
           "unavailable_by_source_label": {"unknown": 1}, "expected_test_neurons": 2,
           "expected_test_animals": 2, "completed_by_part": {"test": 2}}
    put(b / "runs/r1.json", run)
    put(b / "summary.json", run)
    events.append({"run_id": "r1", "event": "run_complete", "completed": 2, "status": "complete"})
    (b / "acquisition.jsonl").write_text("".join(json.dumps(e) + "\n" for e in events))
    return root, base, b, wanted


def run(bundle):
    return A.audit(bundle[0], authorized_test_bytes=True)


def forbid_test_raw(monkeypatch, base):
    original = A.sha_file
    def guarded(path):
        assert not Path(path).is_relative_to(base / "test"), "TEST bytes read before metadata failure"
        return original(path)
    monkeypatch.setattr(A, "sha_file", guarded)


def test_complete_available_and_explicit_unavailable(bundle):
    result = run(bundle)
    assert result["status"] == "PASS", result["issues"]
    assert result["checks"]["fully_rehashed_test_files"] == 1
    assert result["checks"]["unavailable"] == 1


def test_explicit_authorization_required(bundle):
    with pytest.raises(ValueError, match="authorization"):
        A.audit(bundle[0])


@pytest.mark.parametrize("mutation", ["forecast", "freeze", "copy_proof", "a_prefix", "a_receipt", "a_source", "test_url", "request_spacing", "missing_receipt", "extra_test_file", "unavailable_error"])
def test_metadata_corruption_refuses_before_test_byte_reads(bundle, monkeypatch, mutation):
    root, base, b, wanted = bundle
    if mutation == "forecast":
        next((root / "data/interim/axon_engine/real/test/predictions").rglob("*.npz")).write_bytes(b"changed")
    elif mutation == "freeze":
        p = root / "data/interim/axon_engine/real/freeze/frozen.json"
        value = A.read_json(p); value["selected"] = None; put(p, value)
    elif mutation == "copy_proof":
        p = root / "data/interim/axon_engine/real/test/COPIES_VERIFIED.json"
        value = A.read_json(p); value["all_equal"] = False; put(p, value)
    elif mutation == "a_prefix":
        p = base / "MANIFEST.tsv"; p.write_bytes(p.read_bytes().replace(b"legacy", b"edited", 1))
    elif mutation == "a_receipt":
        next((base / "receipts").glob("*.json")).write_text("{}")
    elif mutation == "a_source":
        (root / "analysis/axon_trace_download.py").write_text("changed")
    elif mutation == "test_url":
        with (b / "acquisition.jsonl").open("a") as f:
            f.write(json.dumps({"run_id": "r1", "event": "request_queued", "neuron": "100001_001", "url": A.PORTAL.format("100001_001")}) + "\n")
    elif mutation == "request_spacing":
        p = b / "acquisition.jsonl"
        es = [json.loads(line) for line in p.read_text().splitlines()]
        es[3]["monotonic_start"] = 10.1
        p.write_text("".join(json.dumps(e) + "\n" for e in es))
    elif mutation == "missing_receipt":
        (b / "receipts" / (wanted[0] + ".json")).unlink()
    elif mutation == "extra_test_file":
        (base / "test/extra.swc").write_text("not expected")
    else:
        p = b / "receipts" / (wanted[1] + ".json")
        value = A.read_json(p); value["status"] = 500; put(p, value)
    forbid_test_raw(monkeypatch, base)
    result = run(bundle)
    assert result["status"] == "FAIL"
    assert result["checks"]["fully_rehashed_test_files"] == 0


def test_actual_raw_corruption_detected(bundle):
    next((bundle[1] / "test").rglob("*.swc")).write_bytes(b"bad")
    result = run(bundle)
    assert result["status"] == "FAIL"
    assert any("raw byte/hash mismatch" in s for s in result["issues"])


def test_path_escape_rejected():
    with pytest.raises(ValueError, match="path escape"):
        A.within(Path("/private/tmp/a"), "../outside")


@pytest.mark.parametrize("target", ["a_receipt", "unavailable_b_receipt", "run", "a_snapshot", "extra_receipt", "extra_file"])
def test_late_metadata_or_file_set_change_is_not_accepted(bundle, monkeypatch, target):
    _, base, b, wanted = bundle
    original = A.sha_file
    changed = False
    def changing_sha(path):
        nonlocal changed
        if Path(path).is_relative_to(base / "test") and not changed:
            changed = True
            if target == "a_receipt":
                next((base / "receipts").glob("*.json")).write_text("{}")
            elif target == "unavailable_b_receipt":
                (b / "receipts" / (wanted[1] + ".json")).write_text("{}")
            elif target == "run":
                (b / "runs/r1.json").write_text("{}")
            elif target == "a_snapshot":
                (b / "stage_a_manifest.tsv").write_text("changed")
            elif target == "extra_receipt":
                put(b / "receipts/999999_001.json", {})
            else:
                (base / "test/extra.swc").write_text("extra")
        return original(path)
    monkeypatch.setattr(A, "sha_file", changing_sha)
    result = run(bundle)
    assert changed
    assert result["status"] == "FAIL"
    assert any("changed during audit" in s for s in result["issues"])


def test_actual_downloader_synthetic_resume_bundle_matches_auditor_schema(tmp_path, monkeypatch):
    """Fake transport, real NEW downloader output, independent no-import audit math."""
    path = Path(__file__).with_name("test_axon_trace_download_stage_b.py")
    spec = importlib.util.spec_from_file_location("synthetic_b_schema_helper", path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    rows, plan_sha, freeze_sha, commit_sha = helper.fixture(tmp_path)
    project = Path(__file__).resolve().parents[1]
    for relative in ("analysis/axon_trace_download.py", "analysis/axon_trace_download_stage_b.py"):
        (tmp_path / relative).write_bytes((project / relative).read_bytes())
    clock = helper.H.Clock()
    def downloader():
        transport = helper.H.Transport(clock)
        return helper.B.Downloader(tmp_path, expected_plan_sha=plan_sha, freeze_sha=freeze_sha,
                                   commitment_sha=commit_sha, download=True, opener=transport,
                                   free_space=lambda: 20_000_000_000, clock=clock.now,
                                   sleep=clock.sleep, checker=lambda *_: None)
    assert downloader().run(max_new=1)["status"] == "partial"
    assert downloader().run()["status"] == "complete"
    monkeypatch.setattr(A, "PLAN_SHA", plan_sha)
    monkeypatch.setattr(A, "COUNTS", {"train": 1, "validation": 1, "test": 2})
    monkeypatch.setattr(A, "ANIMALS", {"train": 1, "validation": 1, "test": 2})
    manifest = A.read_json(tmp_path / "data/interim/axon_engine/real/test/prediction_manifest.json")
    monkeypatch.setattr(A, "METHODS", tuple(m["name"] for m in manifest["methods"]))
    result = A.audit(tmp_path, authorized_test_bytes=True)
    assert result["status"] == "PASS", result["issues"]
    assert result["checks"]["runs"] == 2
    assert result["checks"]["fully_rehashed_test_files"] == 2
