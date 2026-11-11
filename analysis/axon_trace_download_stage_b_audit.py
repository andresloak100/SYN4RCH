"""Independent Stage B byte-acquisition audit; no downloader, engine or geometry imports.

TEST metadata and bytes are examined only when the caller explicitly authorizes this
audit after the freeze and global test-forecast commitment. Recorded logs are not an
OS access trace. Stage A metadata/receipts are rechecked; its raw files are not replayed.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import datetime as dt
import hashlib
import io
import json
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
PLAN_SHA = "649f829a592cd3d13bcda5cffdbd049e759568d38bacb144f6a61220610a1ca4"
COUNTS = {"train": 4448, "validation": 1185, "test": 1547}
ANIMALS = {"train": 101, "validation": 25, "test": 40}
PORTAL = "https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/hy/hy_neuron_download/{}.swc"
FIELDS = ("file", "neuron", "bytes", "sha256", "source", "frame", "licence", "how")
METHODS = ("nearest", "medoid", "peptide_nearest", "peptide_medoid", "C1_medoid_copy",
           "C4_uniform_attach_500", "C5_uniform_attach_1000", "F1_conditional_attach_500",
           "F2_conditional_attach_1000", "E1_local_uniform", "E2_local_kernel",
           "E3_local_kernel_label", "E4_local_kernel_regional", "E5_local_kernel_label_regional",
           "E6_local_consensus", "E7_local_kernel_inferred_profile")
UNAVAILABLE = "unavailable_at_source"


def sha_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def rows(data):
    return list(csv.DictReader(io.StringIO(data.decode()), delimiter="\t"))


def within(base, relative):
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError("relative artifact path required")
    p = Path(base) / relative
    if p.is_symlink() or not p.resolve().is_relative_to(Path(base).resolve()):
        raise ValueError("artifact path escape or symlink")
    return p


def audit(root=ROOT, *, authorized_test_bytes=False):
    if not authorized_test_bytes:
        raise ValueError("explicit post-freeze TEST byte-audit authorization required")
    root = Path(root).resolve()
    base = root / "data/raw/jiao2025/traces"
    b = base / "stage_b"
    issues, inputs, files, checks = [], {}, [], {}

    def require(ok, message):
        if not ok:
            issues.append(message)

    def pin(p):
        inputs[str(p.relative_to(root))] = sha_file(p)
        return inputs[str(p.relative_to(root))]

    # Fail before any raw TEST byte read if metadata, proof or file-set binding fails.
    try:
        plan_path = root / "reports/axon_dataset_manifest.tsv"
        require(pin(plan_path) == PLAN_SHA, "locked plan SHA")
        all_rows = rows(plan_path.read_bytes())
        plan = {r["neuron"]: r for r in all_rows}
        require(len(plan) == len(all_rows) == sum(COUNTS.values()), "exact unique plan IDs")
        require(Counter(r["part"] for r in all_rows) == COUNTS, "locked split counts")
        groups = {s: {r["animal"] for r in all_rows if r["part"] == s} for s in COUNTS}
        require({s: len(a) for s, a in groups.items()} == ANIMALS, "locked animal counts")
        require(not any(groups[a] & groups[c] for a, c in (("train", "validation"), ("train", "test"), ("validation", "test"))), "animal split overlap")
        require(all(re.fullmatch(r"\d{6}_\d{3}", n) and r["animal"] == n.split("_")[0] for n, r in plan.items()), "plan ID/animal binding")
        wanted = {n: r for n, r in plan.items() if r["part"] == "test"}
        stage_a = set(plan) - set(wanted)
        local_path = root / "data/interim/axon_dataset/local_manifest.tsv"
        local_sha = pin(local_path)
        local_rows = rows(local_path.read_bytes())
        local = {r["neuron"]: r for r in local_rows}
        require(len(local) == len(local_rows) and set(local) == set(plan), "local exact plan join")
        require(all(all(local.get(n, {}).get(k) == r[k] for k in ("animal", "part")) for n, r in plan.items()), "local split/animal join")

        auth_path = b / "authorization.json"
        auth_sha = pin(auth_path)
        auth = read_json(auth_path)
        require(auth.get("schema") == "synarch_stage_b_authorization_v1" and auth.get("plan_sha256") == PLAN_SHA
                and auth.get("local_manifest_sha256") == local_sha, "authorization declaration binding")
        proof_names = {"freeze", "commitment", "prediction_manifest", "copies_verified"}
        require(set(auth["proof_paths"]) == set(auth["proof_sha256"]) == proof_names, "exact authorization proofs")
        expected_proofs = {"freeze": "freeze/frozen.json", "commitment": "test/PREDICTION_COMMITMENT.json",
                           "prediction_manifest": "test/prediction_manifest.json", "copies_verified": "test/COPIES_VERIFIED.json"}
        require(auth["proof_paths"] == {k: "data/interim/axon_engine/real/" + v for k, v in expected_proofs.items()}, "designated real engine proof paths")
        proofs = {}
        for name in proof_names:
            p = within(root, auth["proof_paths"][name])
            require(pin(p) == auth["proof_sha256"][name], f"changed {name} proof")
            proofs[name] = read_json(p)
        freeze, commit, manifest, copies = [proofs[k] for k in ("freeze", "commitment", "prediction_manifest", "copies_verified")]
        require(freeze.get("schema") == "axon_engine_freeze_v2" and freeze.get("world_kind") == "real"
                and freeze.get("validation_method_equals_frozen_method") is True
                and freeze.get("git", {}).get("engine_files_committed_and_clean") is True, "real engine freeze proof")
        access = freeze["test_access"]
        require(access.get("strongest_seal_no_test_trace_on_disk") is True and access.get("test_traces_listed_on_disk_at_freeze") == 0, "TEST existed at freeze")
        ids_sha = hashlib.sha256("\n".join(sorted(wanted)).encode()).hexdigest()
        require(freeze["pins"]["test_queries"] == len(wanted) and freeze["pins"]["test_queries_sha256"] == ids_sha, "freeze exact TEST set")
        require(commit.get("schema") == "axon_engine_prediction_commitment_v2" and commit.get("stage") == "test"
                and commit.get("manifest_sha256") == auth["proof_sha256"]["prediction_manifest"], "test commitment hash/stage")
        require(manifest.get("schema") == "axon_engine_prediction_manifest_v2" and manifest.get("stage") == "test"
                and manifest.get("part") == "test" and manifest.get("world_kind") == "real", "test forecast manifest schema")
        require(len(manifest["queries"]) == len(wanted) and set(manifest["queries"]) == set(wanted), "forecast exact TEST set")
        require(manifest["family"] == freeze["family"] == "all", "fixed all-candidate family")
        frozen_selection = manifest["frozen_selection"]
        require(frozen_selection["freeze_sha256"] == auth["proof_sha256"]["freeze"]
                and all(frozen_selection[k] == freeze[k] for k in ("selected", "descriptive_reference")), "frozen selection binding")
        frozen_pins = {k: v for k, v in freeze["pins"].items() if k not in ("family", "test_queries", "test_queries_sha256")}
        require(manifest["pins"] == frozen_pins, "frozen method changed")
        require(frozen_pins["cohort"]["manifest"] == PLAN_SHA and frozen_pins["cohort"]["local_manifest"] == local_sha, "frozen cohort input pins")
        method_sources = {}
        for section in ("engine", "frozen", "general_layer"):
            method_sources.update(frozen_pins["code_sha256"][section])
        method_sources["reports/axon_engine_protocol.md"] = frozen_pins["protocol_sha256"]
        for relative, sha in method_sources.items():
            require(pin(within(root, relative)) == sha, "frozen engine source/protocol changed")
        predictions = manifest["predictions"]
        keys = [(p["query_id"], p["model"]) for p in predictions]
        require(len(keys) == len(set(keys)) == len(wanted) * len(METHODS)
                and set(keys) == {(n, m) for n in wanted for m in METHODS}, "complete global TEST forecast matrix")
        require(commit.get("queries") == len(wanted) and commit.get("predictions") == len(keys)
                and set(commit.get("methods", [])) == set(METHODS), "commit forecast counts")
        require(copies.get("stage") == "test" and copies.get("manifest_sha256") == commit["manifest_sha256"]
                and copies.get("all_equal") is True and copies.get("differences") == 0
                and copies.get("queries") == len(wanted) and copies.get("forecasts_compared") == 4 * len(wanted), "four copies not independently verified")
        prediction_base = within(root, auth["proof_paths"]["prediction_manifest"]).parent
        forecast_paths = {str(p.relative_to(prediction_base)) for p in (prediction_base / "predictions").rglob("*.npz")}
        require(forecast_paths == {p["path"] for p in predictions}, "extra/missing forecast files")
        for p in predictions:
            require(sha_file(within(prediction_base, p["path"])) == p["sha256"], "changed forecast bytes")

        pin(b / "stage_a_manifest.tsv")
        prefix = (b / "stage_a_manifest.tsv").read_bytes()
        require(len(prefix) == auth["stage_a_manifest_bytes"] and hashlib.sha256(prefix).hexdigest() == auth["stage_a_manifest_sha256"], "A snapshot bytes/hash")
        a_rows = rows(prefix)
        require(len(a_rows) == len(stage_a) and {r["neuron"] for r in a_rows} == stage_a, "exact A snapshot IDs")
        pin(base / "MANIFEST.tsv")
        acquired = (base / "MANIFEST.tsv").read_bytes()
        require(acquired.startswith(prefix), "A manifest prefix changed")
        manifest_rows = rows(acquired)
        manifest_rows_by_id = {r["neuron"]: r for r in manifest_rows}
        require(len(manifest_rows_by_id) == len(manifest_rows) == len(plan) and set(manifest_rows_by_id) == set(plan), "common manifest exact A+B set")
        require(manifest_rows[:len(a_rows)] == a_rows, "A manifest rows changed")
        require(set(auth["stage_a_receipts_sha256"]) == stage_a, "A receipt snapshot exact IDs")
        require({p.stem for p in (base / "receipts").glob("*.json")} == stage_a, "A receipt file set changed")
        for n, sha in auth["stage_a_receipts_sha256"].items():
            require(pin(base / "receipts" / (n + ".json")) == sha, "A receipt changed")
        require(bool(auth["stage_a_metadata_sha256"]), "missing A metadata snapshot")
        expected_a_metadata = {"summary.json", "acquisition.jsonl"} | {str(p.relative_to(base)) for p in (base / "runs").glob("*.json")}
        require(set(auth["stage_a_metadata_sha256"]) == expected_a_metadata, "exact A metadata snapshot")
        for rel, sha in auth["stage_a_metadata_sha256"].items():
            require(pin(within(base, rel)) == sha, "A metadata changed")
        require(set(auth["source_sha256"]) == {"analysis/axon_trace_download.py", "analysis/axon_trace_download_stage_b.py"}, "exact acquisition source pins")
        for rel, sha in auth["source_sha256"].items():
            require(pin(within(root, rel)) == sha, "acquisition source changed")

        receipt_paths = {p.stem: p for p in (b / "receipts").glob("*.json")}
        require(set(receipt_paths) == set(wanted), "exact B receipt IDs")
        for p in receipt_paths.values():
            pin(p)
        receipts = {n: read_json(p) for n, p in receipt_paths.items()}
        expected_files = set()
        for n, r in receipts.items():
            require(n in wanted, "non-TEST receipt")
            if n not in wanted:
                continue
            require(r.get("schema") == "synarch_stage_b_download_v1" and r.get("stage") == "B"
                    and r.get("plan_sha256") == PLAN_SHA and r.get("local_manifest_sha256") == local_sha
                    and r.get("authorization_sha256") == auth_sha and r.get("proof_sha256") == auth["proof_sha256"], "B receipt proof binding")
            require(r.get("neuron") == n and all(str(r.get(k)) == manifest_rows_by_id[n].get(k) for k in FIELDS), "B receipt/manifest binding")
            require(r.get("source") == r.get("url") == PORTAL.format(n), "B source URL")
            if r.get("availability") == "available":
                rel = f"test/{wanted[n]['animal']}/{n}.swc"
                expected_files.add(rel)
                require(r.get("status") == 200 and r.get("final_url") == PORTAL.format(n)
                        and r.get("frame") == "portal_true_hemisphere" and r.get("file") == rel
                        and type(r.get("bytes")) is int and 0 < r["bytes"] <= 128 * 1024 * 1024
                        and bool(re.fullmatch(r"[0-9a-f]{64}", r.get("sha256", ""))), "available B receipt schema")
            else:
                require(r.get("availability") == UNAVAILABLE and r.get("status") in (404, 410)
                        and r.get("file") == "" and r.get("bytes") == 0 and r.get("sha256") == r.get("frame") == UNAVAILABLE
                        and bool(r.get("how")), "unavailable B receipt must be explicit 404/410")
        actual_files = {str(p.relative_to(base)) for p in (base / "test").rglob("*.swc")}
        require(actual_files == expected_files and not (base / "test").is_symlink()
                and not any(p.is_symlink() for p in (base / "test").rglob("*")), "unexpected/missing/symlink TEST files")
        require(not list((base / "test").rglob("*.part")), "unfinished B partial file")

        run_paths = {p.stem: p for p in (b / "runs").glob("*.json")}
        for p in run_paths.values():
            pin(p)
        runs = {rid: read_json(p) for rid, p in run_paths.items()}
        for rid, run in runs.items():
            require(run.get("run_id") == rid and run.get("stage") == "B" and run.get("authorization_sha256") == auth_sha, "B run identity/proof")
            require(run.get("completed", -1) + len(run.get("unfinished_ids", [])) == len(wanted)
                    and len(set(run.get("unfinished_ids", []))) == len(run.get("unfinished_ids", []))
                    and set(run.get("unfinished_ids", [])) <= set(wanted), "B run counts/membership")
        summary_path = b / "summary.json"
        pin(summary_path)
        summary = read_json(summary_path)
        require(summary.get("run_id") in runs and summary == runs.get(summary.get("run_id")), "B summary differs from durable run")
        require(summary.get("status") == "complete" and summary.get("completed") == len(wanted)
                and summary.get("unfinished_ids") == [] and summary.get("failures") == {}, "B final summary incomplete")
        require(bool(runs) and summary["utc"] == max(r["utc"] for r in runs.values()), "B summary not latest run")
        available = sum(r.get("availability") == "available" for r in receipts.values())
        require(summary.get("acquired") == available and summary.get("unavailable_at_source") == len(wanted) - available, "B availability summary")
        unavailable_ids = sorted(n for n, r in receipts.items() if r.get("availability") == UNAVAILABLE)
        require(len(summary.get("unavailable_ids", [])) == len(unavailable_ids) and set(summary.get("unavailable_ids", [])) == set(unavailable_ids), "B explicit unavailable ID ledger")
        by_animal = {animal: sum(wanted[n]["animal"] == animal for n in unavailable_ids) for animal in sorted(groups["test"])}
        labels = sorted({r.get("source_label", "unknown") for r in wanted.values()})
        by_label = {label: sum(wanted[n].get("source_label", "unknown") == label for n in unavailable_ids) for label in labels}
        require(summary.get("unavailable_by_animal") == by_animal and summary.get("unavailable_by_source_label") == by_label, "B unavailable per-animal/label counts")
        require(summary.get("expected_test_neurons") == len(wanted) and summary.get("expected_test_animals") == len(groups["test"])
                and summary.get("completed_by_part") == {"test": len(wanted)}, "B expected/completed cohort counts")
        events, request_ids, checkpoints, starts, endings = Counter(), Counter(), Counter(), Counter(), Counter()
        last_request = {}
        log_path = b / "acquisition.jsonl"
        pin(log_path)
        for line in log_path.read_text().splitlines():
            e = json.loads(line)
            events[e["event"]] += 1
            require(e.get("run_id") in runs, "B log unknown run")
            n = e.get("neuron")
            require(n is None or n in wanted, "B log non-TEST ID")
            for key in ("url", "source", "final_url"):
                value = e.get(key)
                if isinstance(value, str) and "/hy_neuron_download/" in value:
                    match = re.fullmatch(re.escape(PORTAL).replace(re.escape("{}"), r"(\d{6}_\d{3})"), value)
                    require(bool(match) and match[1] in wanted, "B log non-TEST URL")
            if e["event"] == "request_start":
                require(n in wanted and e.get("url") == PORTAL.format(n) and e.get("authorization_sha256") == auth_sha, "B request URL/authorization")
                start = float(e["monotonic_start"])
                require(math.isfinite(start) and (e["run_id"] not in last_request or start - last_request[e["run_id"]] >= 2 - 1e-6), "B recorded request spacing")
                last_request[e["run_id"]] = start
                require(type(e.get("attempt")) is int and 1 <= e["attempt"] <= 3, "B retry attempt bound")
                request_ids[n] += 1
            if e["event"] == "checkpoint":
                checkpoints[n] += 1
                require(n in receipts and e.get("sha256") == receipts[n]["sha256"] and e.get("bytes") == receipts[n]["bytes"], "B checkpoint mismatch")
            if e["event"] == "run_start":
                starts[e["run_id"]] += 1
            if e["event"] == "run_complete":
                endings[e["run_id"]] += 1
                require(e.get("completed") == runs[e["run_id"]]["completed"] and e.get("status") == runs[e["run_id"]]["status"], "B log ending/summary mismatch")
        require(set(request_ids) == set(wanted), "B requests must cover exact TEST set")
        require(set(checkpoints) == set(wanted) and all(v == 1 for v in checkpoints.values()), "B durable checkpoints exact-once")
        require(set(starts) == set(endings) == set(runs) and all(starts[k] == endings[k] == 1 for k in runs), "B run log lifecycle")
        checks.update(declared_test=len(wanted), declared_test_animals=len(groups["test"]), available=available,
                      unavailable=len(wanted) - available, preserved_stage_a_rows=len(a_rows),
                      unavailable_ids=unavailable_ids, unavailable_by_animal=by_animal, unavailable_by_source_label=by_label,
                      forecast_files_verified=len(predictions), recorded_requests=sum(request_ids.values()),
                      runs=len(runs), events=dict(events), metadata_passed=not issues)
        if not issues:
            def verify(n):
                r = receipts[n]
                p = within(base, r["file"])
                size, actual = p.stat().st_size, sha_file(p)
                return {"neuron": n, "file": r["file"], "bytes": size, "sha256": actual,
                        "receipt_sha256": sha_file(receipt_paths[n])}, size == r["bytes"] and actual == r["sha256"]
            with ThreadPoolExecutor(max_workers=2) as pool:
                for row, ok in pool.map(verify, sorted(n for n, r in receipts.items() if r["availability"] == "available")):
                    files.append(row)
                    require(ok, "TEST raw byte/hash mismatch " + row["neuron"])
            require(sha_file(auth_path) == auth_sha and sha_file(log_path) == inputs[str(log_path.relative_to(root))], "acquisition changed during audit")
            for name, rel in auth["proof_paths"].items():
                require(sha_file(within(root, rel)) == auth["proof_sha256"][name], "proof changed during audit")
            for relative, sha in inputs.items():
                require(sha_file(within(root, relative)) == sha, "metadata changed during audit")
            require({p.stem for p in (base / "receipts").glob("*.json")} == stage_a
                    and {p.stem for p in (b / "receipts").glob("*.json")} == set(wanted)
                    and {p.stem for p in (b / "runs").glob("*.json")} == set(runs), "receipt/run set changed during audit")
            require({str(p.relative_to(base)) for p in (base / "test").rglob("*.swc")} == expected_files
                    and not (base / "test").is_symlink() and not any(p.is_symlink() for p in (base / "test").rglob("*")), "TEST file set changed during audit")
    except (OSError, ValueError, KeyError, TypeError) as e:
        issues.append("incomplete/invalid acquisition metadata: " + type(e).__name__)
    checks.update(fully_rehashed_test_files=len(files), total_test_bytes=sum(r["bytes"] for r in files))
    return {"schema": "synarch_stage_b_acquisition_audit_v1", "status": "PASS" if not issues else "FAIL",
            "utc": dt.datetime.now(dt.timezone.utc).isoformat(), "audit_source_sha256": sha_file(__file__),
            "inputs_sha256": inputs, "checks": checks, "issues": issues, "files": files,
            "limitations": ["Local byte hashes, not independent upstream checksums or geometry admission.",
                            "No geometry parsing, model fits or scoring; Stage A raw bytes are not rehashed.",
                            "Recorded acquisition/proof checks are not OS isolation or protection against rewriting every trusted pin."]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=ROOT)
    p.add_argument("--authorized-test-byte-audit", action="store_true")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    result = audit(a.root, authorized_test_bytes=a.authorized_test_byte_audit)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("x") as f:
        json.dump(result, f, sort_keys=True, indent=2, allow_nan=False)
        f.write("\n")
    print(json.dumps({"status": result["status"], "checks": result["checks"], "issues": result["issues"]}))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
