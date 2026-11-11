"""Synthetic proof bundles and fake portal only. No real TEST request or read."""
from copy import deepcopy
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import ssl
import sys
import urllib.error

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("stage_b_downloader_tested", ROOT / "analysis/axon_trace_download_stage_b.py")
B = importlib.util.module_from_spec(SPEC); sys.modules[SPEC.name] = B; SPEC.loader.exec_module(B)
HELPER_SPEC = importlib.util.spec_from_file_location("stage_a_synthetic_helpers", ROOT / "tests/test_axon_trace_download.py")
H = importlib.util.module_from_spec(HELPER_SPEC); HELPER_SPEC.loader.exec_module(H)
BODY = H.BODY


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True); path.write_text(B.A.canonical(value))


def fixture(root):
    rows = [{"neuron": n, "animal": n.split("_")[0], "part": part, "source_label": label}
            for n, part, label in (("100001_001", "train", "A"), ("100002_001", "validation", "A"),
                                   ("100003_001", "test", "A"), ("100004_001", "test", "B"))]
    plan = root / "reports/axon_dataset_manifest.tsv"; local = root / "data/interim/axon_dataset/local_manifest.tsv"
    H.tsv(plan, rows); H.tsv(local, rows); plan_sha = B.sha(plan)
    base = root / "data/raw/jiao2025/traces"; a_rows = []
    for r in rows[:2]:
        relative = f"{r['part']}/{r['animal']}/{r['neuron']}.swc"; path = base / relative
        path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(BODY)
        receipt = dict(file=relative, neuron=r["neuron"], bytes=len(BODY), sha256=B.sha(path), source="synthetic_A",
                       frame="portal_true_hemisphere", licence="synthetic", how="already acquired Stage A")
        a_rows.append(receipt); save(base / "receipts" / (r["neuron"] + ".json"), receipt)
    H.tsv(base / "MANIFEST.tsv", a_rows, B.FIELDS)
    a_summary = dict(stage="A", run_id="runA", status="complete", completed=2)
    save(base / "summary.json", a_summary); save(base / "runs/runA.json", a_summary)
    (base / "acquisition.jsonl").write_text(B.A.canonical(dict(run_id="runA", event="run_complete")))
    for relative in ("analysis/axon_engine.py", "reports/axon_engine_protocol.md", "tests/test_axon_engine.py", "tests/test_synarch_axon_bridge.py"):
        p = root / relative; p.parent.mkdir(parents=True, exist_ok=True); p.write_text("fixed synthetic immutable method\n")
    wanted = sorted(r["neuron"] for r in rows if r["part"] == "test")
    pins = dict(cohort={"manifest": plan_sha, "local_manifest": B.sha(local)},
                test_queries=len(wanted), test_queries_sha256=B.ids_sha(wanted), family="all",
                code_sha256={"engine": {"analysis/axon_engine.py": B.sha(root / "analysis/axon_engine.py")}, "frozen": {}, "general_layer": {}},
                protocol_sha256=B.sha(root / "reports/axon_engine_protocol.md"))
    freeze = dict(schema="axon_engine_freeze_v2", world_kind="real", validation_method_equals_frozen_method=True,
                  git={"engine_files_committed_and_clean": True}, pins=pins, family="all", selected="candidate", descriptive_reference="candidate",
                  test_access={"test_traces_listed_on_disk_at_freeze": 0, "strongest_seal_no_test_trace_on_disk": True})
    save(root / B.FROZEN, freeze); freeze_sha = B.sha(root / B.FROZEN)
    test_dir = (root / B.COMMITMENT).parent; methods = ["candidate", "nearest", "medoid", "peptide_nearest", "peptide_medoid"]
    predictions = []
    for n in wanted:
        for model in methods:
            rel = f"predictions/{model}/{n}.npz"; path = test_dir / rel; path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"synthetic committed PREDICTION only; never parsed")
            predictions.append(dict(query_id=n, model=model, path=rel, sha256=B.sha(path)))
    manifest = dict(schema="axon_engine_prediction_manifest_v2", stage="test", part="test", world_kind="real", queries=wanted,
                    pins={k: v for k, v in pins.items() if k not in ("family", "test_queries", "test_queries_sha256")},
                    code_sha256=pins["code_sha256"], family="all", methods=[{"name": m} for m in methods], predictions=predictions,
                    frozen_selection={"selected": "candidate", "freeze_sha256": freeze_sha, "descriptive_reference": "candidate"})
    save(test_dir / "prediction_manifest.json", manifest); msha = B.sha(test_dir / "prediction_manifest.json")
    commit = dict(schema="axon_engine_prediction_commitment_v2", stage="test", manifest_sha256=msha, queries=len(wanted), predictions=len(predictions), methods=methods)
    save(root / B.COMMITMENT, commit)
    save(test_dir / "COPIES_VERIFIED.json", dict(stage="test", manifest_sha256=msha, all_equal=True, differences=0, queries=len(wanted), forecasts_compared=4 * len(wanted)))
    return rows, plan_sha, freeze_sha, B.sha(root / B.COMMITMENT)


@pytest.fixture
def setup(tmp_path):
    rows, plan_sha, freeze_sha, commit_sha = fixture(tmp_path); clock = H.Clock()
    def make(script=None, *, download=True, checker=lambda *_a: None):
        transport = H.Transport(clock, script)
        model = B.Downloader(tmp_path, expected_plan_sha=plan_sha, freeze_sha=freeze_sha, commitment_sha=commit_sha,
                             download=download, opener=transport, free_space=lambda: 20_000_000_000,
                             clock=clock.now, sleep=clock.sleep, checker=checker)
        return model, transport
    return tmp_path, rows, plan_sha, freeze_sha, commit_sha, make


@pytest.mark.parametrize("missing", ["go", "freeze", "commitment", "copies"])
def test_no_transport_or_individual_TEST_URL_before_complete_proof(setup, missing, monkeypatch):
    root, _, _, _, _, make = setup
    if missing != "go":
        path = {"freeze": root / B.FROZEN, "commitment": root / B.COMMITMENT,
                "copies": (root / B.COMMITMENT).parent / "COPIES_VERIFIED.json"}[missing]; path.unlink()
    monkeypatch.setattr(B.A.urllib.request, "Request", lambda *_a, **_k: pytest.fail("Constructed TEST request before proof"))
    with pytest.raises((B.A.AcquisitionError, FileNotFoundError)): make(download=missing != "go")
    assert not (root / "data/raw/jiao2025/traces/test").exists()


def test_exact_TEST_membership_and_no_A_URLs(setup):
    root, rows, _, _, _, make = setup; d, transport = make()
    assert set(d.rows) == {r["neuron"] for r in rows if r["part"] == "test"}
    with pytest.raises(B.A.AcquisitionError): d.url(rows[0]["neuron"])
    with pytest.raises(KeyError): d.target(rows[0]["neuron"])
    assert not transport.requests


def test_fake_full_download_preserves_A_bytes_metadata_receipts_and_rate(setup):
    root, rows, _, _, _, make = setup; base = root / "data/raw/jiao2025/traces"
    a_paths = [base / "MANIFEST.tsv", base / "summary.json", base / "acquisition.jsonl", base / "runs/runA.json",
               *list((base / "receipts").glob("*.json")), *list((base / "train").rglob("*.swc")), *list((base / "validation").rglob("*.swc"))]
    before = {str(p): p.read_bytes() for p in a_paths}; d, transport = make(); result = d.run()
    assert result["status"] == "complete" and result["completed"] == result["acquired"] == 2
    assert {n for n, _, _ in transport.requests} == set(d.rows)
    assert all(b[1] - a[1] >= 2 for a, b in zip(transport.requests, transport.requests[1:]))
    for p in a_paths:
        if p.name == "MANIFEST.tsv": assert p.read_bytes().startswith(before[str(p)])
        else: assert p.read_bytes() == before[str(p)]
    assert len(B.A.read_tsv(base / "MANIFEST.tsv")) == 4
    for n in d.rows:
        receipt = B.read(d.state / "receipts" / (n + ".json"))
        assert d.target(n).read_bytes() == BODY and receipt["availability"] == "available"
        assert receipt["authorization_sha256"] == B.sha(d.state / "authorization.json")
    assert B.read(d.state / "summary.json") == B.read(d.state / "runs" / (d.run_id + ".json"))


def test_partial_then_verified_resume_repairs_B_suffix_only(setup):
    _, _, _, _, _, make = setup; d, _ = make(); first = d.run(max_new=1)
    assert first["status"] == "partial"
    (d.base / "MANIFEST.tsv").write_bytes((d.state / "stage_a_manifest.tsv").read_bytes())
    resumed, transport = make(); result = resumed.run()
    assert result["status"] == "complete" and len(transport.requests) == 1
    assert len(B.A.read_tsv(d.base / "MANIFEST.tsv")) == 4


@pytest.mark.parametrize("status", [404, 410])
def test_source_unavailable_has_explicit_receipt_and_manifest_row(setup, status):
    _, _, _, _, _, make = setup
    n = "100003_001"; error = urllib.error.HTTPError("synthetic", status, "no trace", {}, None)
    d, _ = make({n: [error]}); result = d.run()
    assert result["status"] == "complete" and result["completed"] == 2
    assert result["unavailable_at_source"] == 1 and result["unavailable_by_animal"]["100003"] == 1
    r = B.read(d.state / "receipts" / (n + ".json"))
    assert r["frame"] == r["sha256"] == B.UNAVAILABLE and r["bytes"] == 0 and not r["file"] and not d.target(n).exists()
    assert next(r for r in B.A.read_tsv(d.base / "MANIFEST.tsv") if r["neuron"] == n)["frame"] == B.UNAVAILABLE


def test_changed_freeze_source_or_prediction_rejected_before_transport(setup):
    root, _, _, _, _, make = setup; d, transport = make()
    (root / "analysis/axon_engine.py").write_text("changed method\n")
    with pytest.raises(B.A.StopAcquisition): d.url(d.order[0])
    assert not transport.requests
    with pytest.raises(B.A.AcquisitionError): make()


def test_corrupt_B_raw_file_not_silently_refetched(setup):
    _, _, _, _, _, make = setup; d, _ = make(); d.run(max_new=1)
    d.target(d.order[0]).write_bytes(BODY + b"changed")
    resumed, transport = make(); result = resumed.run(max_new=0)
    assert result["status"] == "failed" and d.order[0] in result["failures"] and not transport.requests


@pytest.mark.parametrize("problem", ["401", "ssl", "wrapped_ssl", "html_header", "html_body", "redirect", "long_retry"])
def test_fatal_access_responses_stop_workers_without_alternate_route(setup, problem):
    _, _, _, _, _, make = setup
    if problem == "401": action = urllib.error.HTTPError("synthetic", 401, "auth", {}, None)
    elif problem == "ssl": action = ssl.SSLError("certificate")
    elif problem == "wrapped_ssl": action = urllib.error.URLError(ssl.SSLError("certificate"))
    elif problem == "html_header": action = lambda u: H.Response(u, headers={"Content-Type": "text/html"})
    elif problem == "html_body": action = lambda u: H.Response(u, b"<html>login</html>")
    elif problem == "redirect": action = lambda u: H.Response("https://synthetic.invalid/login")
    else: action = urllib.error.HTTPError("synthetic", 429, "slow", {"Retry-After": "3600"}, None)
    d, transport = make({n: [action] for n in ("100003_001", "100004_001")}); result = d.run()
    assert result["status"] == "failed" and d.stop.is_set() and result["unavailable_at_source"] == 0
    # Body-only failures are discovered after headers release the rate gate;
    # another already-started transfer may precede discovery. New starts after
    # discovery must still be impossible. Header/access/TLS failures stop inside
    # the gate and therefore admit exactly one start.
    assert len(transport.requests) in ((1, 2) if problem == "html_body" else (1,))
    starts = len(transport.requests); d.job(d.order[0]); assert len(transport.requests) == starts


def test_retry_is_bounded_and_shared_cooldown_then_success(setup):
    _, _, _, _, _, make = setup
    action = urllib.error.HTTPError("synthetic", 429, "wait", {"Retry-After": "7"}, None)
    d, transport = make({"100003_001": [action]}); result = d.run()
    assert result["status"] == "complete" and len(transport.requests) == 3
    first = next(i for i, r in enumerate(transport.requests) if r[0] == "100003_001")
    at = transport.requests[first][1]
    assert all(t >= at + 7 for _, t, _ in transport.requests[first + 1:])


def test_dirty_committed_engine_guard_is_mandatory(setup):
    _, _, _, _, _, make = setup
    def dirty(*args): raise B.A.AcquisitionError("uncommitted/dirty method")
    with pytest.raises(B.A.AcquisitionError, match="dirty"): make(checker=dirty)


def test_plan_is_metadata_only_without_freeze_or_test_URLs(setup, monkeypatch):
    root, _, plan_sha, _, _, _ = setup
    original_metadata = B.metadata
    monkeypatch.setattr(B, "metadata", lambda r: original_metadata(r, plan_sha))
    monkeypatch.setattr(B.A.urllib.request, "Request", lambda *_a, **_k: pytest.fail("plan made a request"))
    (root / B.FROZEN).unlink(); (root / B.COMMITMENT).unlink()
    result = B.plan(root)
    assert result["requests"] == 0 and result["test_neurons"] == 2 and result["test_animals"] == 2
    assert not (root / "data/raw/jiao2025/traces/stage_b").exists()


@pytest.mark.parametrize("problem", ["sixteen_wrong_names", "wrong_family", "duplicate"])
def test_production_exact_method_names_and_family_before_url_or_transport(setup, monkeypatch, problem):
    root, rows, plan_sha, freeze_sha, _, _ = setup
    manifest_path = (root / B.COMMITMENT).parent / "prediction_manifest.json"
    manifest = B.read(manifest_path); commit = B.read(root / B.COMMITMENT)
    if problem == "sixteen_wrong_names": names = [f"wrong{i}" for i in range(16)]
    elif problem == "duplicate": names = [B.METHODS[0]] * 16
    else: names = list(B.METHODS); manifest["family"] = "no_fitted_state"
    manifest["methods"] = [{"name": name} for name in names]; commit["methods"] = names
    if problem == "wrong_family":
        freeze = B.read(root / B.FROZEN); freeze["family"] = manifest["family"]
        save(root / B.FROZEN, freeze); freeze_sha = B.sha(root / B.FROZEN)
        manifest["frozen_selection"]["freeze_sha256"] = freeze_sha
    save(manifest_path, manifest); commit["manifest_sha256"] = B.sha(manifest_path); save(root / B.COMMITMENT, commit)
    # Force the production-method branch on synthetic metadata only. Its error
    # must precede the deliberately incomplete32-forecast matrix/source checks.
    monkeypatch.setattr(B, "PLAN_SHA", plan_sha)
    expected = "Duplicate" if problem == "duplicate" else "approved16methods"
    with pytest.raises(B.A.AcquisitionError, match=expected):
        B.verify_proofs(root, {r["neuron"]: r for r in rows}, B.sha(root / "data/interim/axon_dataset/local_manifest.tsv"),
                        B.FROZEN, freeze_sha, B.COMMITMENT, B.sha(root / B.COMMITMENT), checker=lambda *_a: None)
    B.method_contract(list(B.METHODS), "all", production=True)
