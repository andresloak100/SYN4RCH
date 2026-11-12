"""Acquire sealed MouseLight test truth, never fit or score it.

Import/check are read-only. Only the explicit acquire/resume commands can request
the exact frozen separate-axon UUIDs. An interrupted request consumes its full
durable reservation; completed checkpoints are resumed by hashes, not geometry.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import closing
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import urllib.error
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
RUN_ID = "mouselight_external_geometry_v1"
OUT = Path("data/interim/mouselight_external_test") / RUN_ID
RAW = Path("data/raw/mouselight_external_test") / RUN_ID
PROTOCOL = Path("docs/mouselight_external_protocol.json")
CODE = "analysis/mouselight_external_acquire.py"
TESTS = "tests/test_mouselight_external_acquire.py"
ADAPTER = "analysis/mouselight_external_test.py"
ADAPTER_TESTS = "tests/test_mouselight_external_test.py"
SOURCE = "data/interim/paired_projection_search_round6/image_routes/calibration/acquisition_source_manifest.json"
CONSTANTS = "data/interim/paired_projection_search_round6/image_routes/calibration/official_constants_metadata.json"
ENDPOINT = "https://ml-neuronbrowser.janelia.org/tracings"
METHODS = ("nearest_copy", "global_medoid_copy", "E1_local_uniform", "E2_local_kernel",
           "E4_local_kernel_regional", "E6_local_consensus")
MAX_TOTAL_BYTES = 4 * 1024 ** 3
MAX_RESPONSE_BYTES = 256 * 1024 ** 2
# One extra byte detects an over-limit chunked response, and is reserved too.
RESERVATION = MAX_RESPONSE_BYTES + 1


class Refusal(ValueError):
    pass


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def now():
    return datetime.now(timezone.utc).isoformat()


def local(root, relative):
    p = Path(relative)
    if p.is_absolute() or ".." in p.parts or p.as_posix() != str(relative):
        raise Refusal("noncanonical repository-relative path")
    result = Path(root) / p
    if result.resolve() != Path(root).resolve() / p:
        raise Refusal("symlink or path alias")
    return result


def write_new(path, value):
    with Path(path).open("x") as stream:
        stream.write(canonical(value))
        stream.flush()
        os.fsync(stream.fileno())


def committed(root, revision, paths):
    if not isinstance(revision, str) or len(revision) != 40 or any(c not in "0123456789abcdef" for c in revision):
        raise Refusal("invalid reviewed commit")
    for name in paths:
        result = subprocess.run(["git", "show", f"{revision}:{name}"], cwd=root, capture_output=True)
        if result.returncode or hashlib.sha256(result.stdout).hexdigest() != sha(local(root, name)):
            raise Refusal("reviewed acquisition source/protocol is uncommitted or changed: " + name)


def source_records(root, protocol):
    spec = protocol.get("acquisition", {})
    if spec.get("root_admission") != "READY" or spec.get("root_reviewed") is not True:
        raise Refusal("explicit ROOT acquisition admission is not READY")
    for key, path in (("source_manifest", SOURCE), ("constants", CONSTANTS)):
        item = spec.get(key, {})
        if item.get("path") != path or item.get("sha256") != sha(local(root, path)):
            raise Refusal("acquisition metadata/constants pin differs")
    if spec.get("code") != CODE or spec.get("tests") != TESTS:
        raise Refusal("acquisition implementation identity differs")
    source = read(local(root, SOURCE))
    if source.get("schema") != "mouselight_acquisition_source_v1" or source.get("endpoint") != ENDPOINT:
        raise Refusal("unreviewed acquisition route")
    records = source.get("records", [])
    if len(records) != 97 or len({r["neuron"] for r in records}) != 97:
        raise Refusal("source UUID cohort must contain exactly97 distinct neurons")
    if len({r["sample_uuid"] for r in records}) != 22 or len({r["axon_tracing_uuid"] for r in records}) != 97:
        raise Refusal("source animal/tracing cohort differs")
    if [r["neuron"] for r in records] != sorted(r["neuron"] for r in records):
        raise Refusal("source request order is not canonical")
    for row in records:
        if not row["neuron"].startswith("AA") or not row["neuron"][2:].isdigit():
            raise Refusal("invalid source neuron identity")
        if row.get("trace_domain") != "separate_axon_tracing" or row.get("tracing_structure") != {"name": "axon", "value": 1}:
            raise Refusal("source UUID is not an explicitly identified separate axon")
        if not row.get("license") or not row.get("doi") or not row.get("license_metadata_sha256"):
            raise Refusal("per-neuron provenance/licence is unresolved")
    return records


def gate(root=ROOT, resume=False):
    """Verify only metadata, sealed hashes and commits; do not parse truth or NPZs."""
    p = read(local(root, str(PROTOCOL)))
    if p.get("run_id") != RUN_ID or p.get("status") != "READY_FOR_ROOT_FREEZE" or p.get("admission", {}).get("root_reviewed") is not True or p["admission"].get("blockers"):
        raise Refusal("external protocol admission is unresolved")
    if p.get("output_dir") != str(OUT) or p.get("raw_dir") != str(RAW) or tuple(p.get("methods", [])) != METHODS:
        raise Refusal("fixed experiment identity differs")
    records = source_records(root, p)
    out, raw = local(root, str(OUT)), local(root, str(RAW))
    if (out / "SCORE_STARTED.json").exists() or (out / "ACQUISITION.json").exists():
        raise Refusal("acquisition already complete or scoring started")
    frozen = read(out / "FROZEN.json")
    commit = read(out / "PREDICTIONS_COMMITTED.json")
    manifest = read(out / "prediction_manifest.json")
    if frozen.get("run_id") != RUN_ID or commit.get("run_id") != RUN_ID or manifest.get("run_id") != RUN_ID:
        raise Refusal("sealed run identity differs")
    if commit.get("truth_absent") is not True or frozen.get("external_truth_absent") is not True:
        raise Refusal("external truth was not absent when predictions were committed")
    if commit.get("freeze_sha256") != sha(out / "FROZEN.json") or manifest.get("freeze_sha256") != commit["freeze_sha256"] or commit.get("manifest_sha256") != sha(out / "prediction_manifest.json"):
        raise Refusal("prediction/freeze commitment differs")
    files = frozen["pins"]["files"]
    required = (CODE, TESTS, ADAPTER, ADAPTER_TESTS, str(PROTOCOL), SOURCE, CONSTANTS, p["metadata"]["path"])
    if any(name not in files for name in required):
        raise Refusal("frozen source pins do not include every acquisition dependency")
    for name, expected in files.items():
        if sha(local(root, name)) != expected:
            raise Refusal("frozen source/input pin differs: " + name)
    if sha(Path(__file__)) != files[CODE]:
        raise Refusal("the executing acquisition code differs from the frozen code")
    committed(root, frozen.get("reviewed_commit"), (CODE, TESTS, ADAPTER, ADAPTER_TESTS, str(PROTOCOL)))
    queries = frozen["queries"]
    if frozen.get("queries_sha256") != digest(queries) or manifest.get("queries_sha256") != frozen["queries_sha256"]:
        raise Refusal("frozen query commitment differs")
    expected_queries = {r["neuron"]: r for r in records}
    if len(queries) != len(records) or {q["neuron"] for q in queries} != set(expected_queries):
        raise Refusal("source UUID manifest differs from frozen cohort")
    for q in queries:
        row = expected_queries[q["neuron"]]
        xyz = row["ccfv3_soma_um"]
        expected_soma = [xyz["z"], xyz["y"], 11400.0 - xyz["x"] if xyz["x"] > 5700 else xyz["x"]]
        if q["animal"] != "mouselight:" + row["sample_uuid"] or list(q["soma_ap_dv_ml_um"]) != expected_soma or q["hemisphere"] != ("Right" if xyz["x"] > 5700 else "Left"):
            raise Refusal("source sample/soma differs from frozen query")
    if tuple(manifest.get("methods", [])) != METHODS:
        raise Refusal("prediction method family differs")
    expected_pairs = {(q["neuron"], m) for q in queries for m in METHODS}
    seen = set()
    animals = {q["neuron"]: q["animal"] for q in queries}
    for row in manifest["predictions"]:
        pair = (row["query"], row["method"])
        path = f"predictions/{pair[0]}/{pair[1]}.npz"
        if pair not in expected_pairs or pair in seen or row["path"] != path or row["animal"] != animals.get(pair[0]) or sha(local(out, path)) != row["sha256"]:
            raise Refusal("missing, duplicated or changed sealed prediction")
        seen.add(pair)
    if seen != expected_pairs or commit.get("prediction_count") != len(seen):
        raise Refusal("not every frozen forecast was committed")
    binding = {"run_id": RUN_ID, "prediction_commitment_sha256": sha(out / "PREDICTIONS_COMMITTED.json"),
               "freeze_sha256": sha(out / "FROZEN.json"), "source_manifest_sha256": sha(local(root, SOURCE)),
               "constants_sha256": sha(local(root, CONSTANTS))}
    started = out / "ACQUISITION_STARTED.json"
    if resume:
        old = read(started)
        if any(old.get(k) != v for k, v in binding.items()) or not raw.is_dir():
            raise Refusal("resume belongs to another commitment or source")
        if datetime.fromisoformat(old["started_utc"]) <= datetime.fromisoformat(commit["utc"]):
            raise Refusal("acquisition start predates prediction commitment")
    elif raw.exists() or started.exists():
        raise Refusal("canonical RAW/start marker already exists; only explicit hash-checked resume is allowed")
    return {"protocol": p, "records": records, "binding": binding, "commit": commit, "out": out, "raw": raw}


def convert(raw_bytes, record, constants):
    """The calibrated conversion: values, IDs and role codes preserved exactly."""
    obj = json.loads(raw_bytes)
    tracings = obj.get("tracings")
    if tracings == [] and set(obj) <= {"tracings", "timing"}:
        return None
    if not isinstance(tracings, list) or len(tracings) != 1 or tracings[0].get("id") != record["axon_tracing_uuid"]:
        raise Refusal("response does not contain exactly the requested separate axon UUID")
    nodes = tracings[0].get("nodes")
    if not isinstance(nodes, list) or not nodes:
        raise Refusal("source axon has no node array")
    mapping = {v["id"]: v["value"] for v in constants["data"]["structureIdentifiers"]}
    by_id, children, roots = {}, {}, []
    for row in nodes:
        if not isinstance(row, dict) or type(row.get("sampleNumber")) is not int or type(row.get("parentNumber")) is not int:
            raise Refusal("invalid source node/parent identity")
        i = row["sampleNumber"]
        if i <= 0 or i in by_id or mapping.get(row.get("structureIdentifierId")) not in {0, 1, 2, 5, 6}:
            raise Refusal("duplicate node or uncalibrated source node role")
        if any(type(row.get(k)) not in (int, float) or not math.isfinite(row[k]) for k in ("x", "y", "z", "radius")):
            raise Refusal("nonfinite/non-numeric source geometry")
        by_id[i], children[i] = row, []
        if row["parentNumber"] == -1:
            roots.append(i)
    if roots != [1] and set(roots) != {1}:
        raise Refusal("source tree does not have one node1 root")
    if len(roots) != 1 or mapping[by_id[1]["structureIdentifierId"]] != 1:
        raise Refusal("source root is not a single type1 soma")
    soma = record["ccfv3_soma_um"]
    if math.dist([by_id[1][k] for k in ("x", "y", "z")], [soma[k] for k in ("x", "y", "z")]) > 5:
        raise Refusal("source root differs from pinned atlas soma by more than5um")
    for i, row in by_id.items():
        if row["parentNumber"] != -1:
            if row["parentNumber"] not in by_id or mapping[row["structureIdentifierId"]] == 1:
                raise Refusal("missing parent or second soma")
            children[row["parentNumber"]].append(i)
    seen, stack = set(), [1]
    while stack:
        i = stack.pop()
        if i in seen:
            raise Refusal("cycle in source tree")
        seen.add(i)
        stack.extend(children[i])
    if len(seen) != len(nodes):
        raise Refusal("disconnected source tree")
    lines = ["# Derived from the official separate axon JSON; unchanged source ML,DV,AP micrometres.",
             "# sampleNumber structureIdentifier.value x y z radius parentNumber"]
    for i in sorted(by_id):
        row = by_id[i]
        lines.append(" ".join([str(i), str(mapping[row["structureIdentifierId"]]),
                               *(format(row[k], ".17g") for k in ("x", "y", "z", "radius")), str(row["parentNumber"])]))
    return ("\n".join(lines) + "\n").encode()


def checkpoint_files(root, state):
    """Resume completed transactions using hashes only, never by parsing axons."""
    folder = state["out"] / "acquisition_checkpoints"
    result = {}
    expected = {r["neuron"]: r for r in state["records"]}
    for path in sorted(folder.glob("*.json")):
        item = read(path)
        neuron = path.stem
        seal = path.with_suffix(".sha256")
        if not seal.is_file() or seal.read_text().strip() != sha(path):
            raise Refusal("checkpoint metadata/status seal differs")
        if neuron not in expected or item.get("neuron") != neuron or item.get("source_manifest_sha256") != state["binding"]["source_manifest_sha256"] or item.get("axon_tracing_uuid") != expected[neuron]["axon_tracing_uuid"]:
            raise Refusal("checkpoint source identity differs")
        record = expected[neuron]
        for key in ("sample_uuid", "doi", "license", "license_metadata_sha256"):
            if item.get(key) != record[key]:
                raise Refusal("checkpoint sample/licence provenance differs")
        if item.get("frame") != "CCFv3" or item.get("trace_domain") != "separate_axon_tracing" or item.get("source_url") != ENDPOINT:
            raise Refusal("checkpoint source frame/domain differs")
        if item.get("status") not in {"downloaded", "unavailable_at_source"}:
            raise Refusal("unrecognized checkpoint status")
        attempt = item.get("attempt_id")
        if not isinstance(attempt, str) or len(attempt) != 6 or not attempt.isdigit():
            raise Refusal("checkpoint transport identity differs")
        endpath = state["out"] / "acquisition_attempts" / f"{attempt}.end.json"
        startpath = state["out"] / "acquisition_attempts" / f"{attempt}.start.json"
        end, start = read(endpath), read(startpath)
        if sha(endpath) != item.get("transport_receipt_sha256") or end.get("start_sha256") != sha(startpath) or start.get("neuron") != neuron or start.get("axon_tracing_uuid") != record["axon_tracing_uuid"] or end.get("status") != "transport_complete" or end.get("raw_json") != item.get("raw_json"):
            raise Refusal("checkpoint transport receipt differs")
        if item["raw_json"]["path"] != str(RAW / "raw_json" / f"{neuron}.{attempt}.json"):
            raise Refusal("noncanonical raw response path")
        if item["status"] == "downloaded":
            pin = item.get("swc", {})
            if pin.get("path") != str(RAW / f"{neuron}.swc") or any(item.get(k) != pin.get(k) for k in ("path", "sha256", "bytes")):
                raise Refusal("canonical SWC pin differs")
        elif not item.get("reason") or (state["raw"] / f"{neuron}.swc").exists():
            raise Refusal("source-unavailable checkpoint conflicts with acquired SWC")
        for kind in ("raw_json",) + (("swc",) if item["status"] == "downloaded" else ()):
            pin = item[kind]
            file = local(root, pin["path"])
            if sha(file) != pin["sha256"] or file.stat().st_size != pin["bytes"]:
                raise Refusal("checkpoint bytes changed")
        result[neuron] = item
    for path in state["raw"].glob("*.swc"):
        if path.stem not in result:
            raise Refusal("uncheckpointed canonical SWC requires ROOT review; never overwrite")
    return result


def budget_used(out):
    used = 0
    for path in sorted((out / "acquisition_attempts").glob("*.start.json")):
        start = read(path)
        if start.get("reservation_bytes") != RESERVATION:
            raise Refusal("attempt budget reservation differs")
        complete = path.with_name(path.name.replace(".start.json", ".end.json"))
        if complete.exists():
            end = read(complete)
            if end.get("start_sha256") != sha(path) or type(end.get("charged_bytes")) is not int or not 0 <= end["charged_bytes"] <= RESERVATION:
                raise Refusal("attempt budget receipt differs")
            used += end["charged_bytes"]
        else:
            used += RESERVATION                 # conservatively charge interrupted requests
    return used


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def open_response(request):
    return urllib.request.build_opener(NoRedirect).open(request, timeout=60)


def request_one(root, state, record, fetch):
    out, raw = state["out"], state["raw"]
    if budget_used(out) + RESERVATION > MAX_TOTAL_BYTES:
        raise Refusal("cumulative transport budget cannot reserve another response; stop without excluding a neuron")
    folder = out / "acquisition_attempts"
    sequence = len(list(folder.glob("*.start.json"))) + 1
    start = folder / f"{sequence:06d}.start.json"
    payload = {"ids": [record["axon_tracing_uuid"]], "ccfVersion": "CCFV30"}
    write_new(start, {**state["binding"], "neuron": record["neuron"], "axon_tracing_uuid": record["axon_tracing_uuid"],
                      "utc": now(), "reservation_bytes": RESERVATION, "endpoint": ENDPOINT, "request": payload})
    relative = RAW / "raw_json" / f"{record['neuron']}.{sequence:06d}.json"
    target = local(root, str(relative))
    received = 0
    headers, status = {}, None
    try:
        request = urllib.request.Request(ENDPOINT, data=canonical(payload).encode(),
                                        headers={"Content-Type": "application/json", "Accept-Encoding": "identity"}, method="POST")
        with closing(fetch(request)) as response, target.open("xb") as stream:
            status, headers = response.status, dict(response.headers)
            if response.geturl() != ENDPOINT or status != 200:
                raise Refusal("redirected or non200 source response")
            encoding = response.headers.get("Content-Encoding", "identity").lower()
            if encoding != "identity":
                raise Refusal("encoded transport is not the reviewed uncompressed source route")
            declared = response.headers.get("Content-Length")
            if declared is not None and (int(declared) < 0 or int(declared) > MAX_RESPONSE_BYTES):
                raise Refusal("source Content-Length exceeds response budget")
            while True:
                block = response.read(min(1 << 20, RESERVATION - received))
                if not block:
                    break
                received += len(block)
                stream.write(block)
                if received > MAX_RESPONSE_BYTES:
                    raise Refusal("source response exceeds256MiB; stop without excluding a neuron")
            if declared is not None and received != int(declared):
                raise Refusal("source Content-Length does not match received bytes")
            stream.flush()
            os.fsync(stream.fileno())
        end = {"start_sha256": sha(start), "charged_bytes": received, "status": "transport_complete",
               "http_status": status, "headers": headers, "raw_json": {"path": str(relative), "bytes": received, "sha256": sha(target)}}
        write_new(start.with_name(start.name.replace(".start.json", ".end.json")), end)
        return target.read_bytes(), {**end, "attempt_id": f"{sequence:06d}", "receipt_sha256": sha(start.with_name(start.name.replace(".start.json", ".end.json")))}
    except Exception as error:
        endpath = start.with_name(start.name.replace(".start.json", ".end.json"))
        if not endpath.exists():
            # A failed open/stream is charged its full allowance, including unknown transfer.
            write_new(endpath, {"start_sha256": sha(start), "charged_bytes": RESERVATION, "status": "transport_error",
                                "http_status": getattr(error, "code", status), "headers": headers, "received_bytes": received,
                                "error": str(error), "classification": "unresolved_transport_not_source_unavailability"})
        raise Refusal("source request failed; no neuron excluded: " + str(error)) from error


def acquire(root=ROOT, resume=False, fetch=None):
    state = gate(root, resume=resume)
    out, raw = state["out"], state["raw"]
    lock = out / "ACQUISITION_RUNNING.lock"
    with lock.open("x") as stream:
        stream.write(canonical({"run_id": RUN_ID, "pid": os.getpid(), "utc": now()}))
    try:
        if not resume:
            started = {**state["binding"], "started_utc": now(), "max_total_bytes": MAX_TOTAL_BYTES,
                       "max_response_bytes": MAX_RESPONSE_BYTES, "reservation_bytes": RESERVATION}
            if datetime.fromisoformat(started["started_utc"]) <= datetime.fromisoformat(state["commit"]["utc"]):
                raise Refusal("acquisition clock must follow prediction commitment")
            write_new(out / "ACQUISITION_STARTED.json", started)
            raw.parent.mkdir(parents=True, exist_ok=True)  # only after every gate and exclusive start
            raw.mkdir()                     # no existing RAW accepted
            (raw / "raw_json").mkdir()
            (out / "acquisition_attempts").mkdir()
            (out / "acquisition_checkpoints").mkdir()
        completed = checkpoint_files(root, state)
        constants = read(local(root, CONSTANTS))
        for record in state["records"]:
            if record["neuron"] in completed:
                continue                    # only hashes above, no cached answer parsing
            body, transport = request_one(root, state, record, fetch or open_response)
            swc = convert(body, record, constants)
            item = {"neuron": record["neuron"], "source_manifest_sha256": state["binding"]["source_manifest_sha256"],
                    "axon_tracing_uuid": record["axon_tracing_uuid"], "sample_uuid": record["sample_uuid"],
                    "source_url": ENDPOINT, "doi": record["doi"], "license": record["license"],
                    "license_metadata_sha256": record["license_metadata_sha256"], "frame": "CCFv3", "trace_domain": "separate_axon_tracing",
                    "raw_json": transport["raw_json"], "attempt_id": transport["attempt_id"], "transport_receipt_sha256": transport["receipt_sha256"]}
            if swc is None:
                item.update(status="unavailable_at_source", reason="HTTP200 official response contains an empty tracings array for the exact pinned axonUUID; raw response retained.")
            else:
                file = raw / f"{record['neuron']}.swc"
                with file.open("xb") as stream:
                    stream.write(swc)
                    stream.flush()
                    os.fsync(stream.fileno())
                item.update(status="downloaded", path=str(RAW / file.name), sha256=sha(file), bytes=len(swc),
                            swc={"path": str(RAW / file.name), "sha256": sha(file), "bytes": len(swc)})
            checkpoint = out / "acquisition_checkpoints" / f"{record['neuron']}.json"
            write_new(checkpoint, item)
            with checkpoint.with_suffix(".sha256").open("x") as stream:
                stream.write(sha(checkpoint) + "\n")
                stream.flush()
                os.fsync(stream.fileno())
        completed = checkpoint_files(root, state)
        if set(completed) != {r["neuron"] for r in state["records"]}:
            raise Refusal("acquisition is incomplete")
        gate(root, resume=True)              # detect changed seals/source before completion
        receipt = {**state["binding"], "started_utc": read(out / "ACQUISITION_STARTED.json")["started_utc"],
                   "completed_utc": now(), "transport_charged_bytes": budget_used(out), "files": completed,
                   "geometry_scores_computed": False, "conversion": "preserved source sampleNumber,parentNumber,XYZ,radius and official role map; no clipping"}
        write_new(out / "ACQUISITION.json", receipt)
        return {"run_id": RUN_ID, "queries": len(completed), "acquisition_sha256": sha(out / "ACQUISITION.json"),
                "statuses": dict(Counter(v["status"] for v in completed.values())), "charged_bytes": receipt["transport_charged_bytes"]}
    finally:
        lock.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("check", "acquire", "resume"))
    args = parser.parse_args()
    try:
        if args.command == "check":
            result = gate()
            print(canonical({"run_id": RUN_ID, "ready": True, "queries": len(result["records"]), "network_requested": False}))
        else:
            print(canonical(acquire(resume=args.command == "resume")))
    except (Refusal, FileExistsError, FileNotFoundError, KeyError, json.JSONDecodeError) as error:
        print(canonical({"run_id": RUN_ID, "refused": str(error)}))
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
