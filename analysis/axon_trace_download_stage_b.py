"""Post-freeze TEST acquisition only; no geometry parsing or requests on import.

--plan is metadata-only. Actual acquisition requires --download and externally
pinned freeze/commitment hashes. The old Stage A files/code are never rewritten;
only its common manifest gains an exact-prefix-preserving TEST suffix.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import csv
import fcntl
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import ssl
import subprocess
import threading
import time
import urllib.error
import urllib.request
import uuid

ROOT = Path(__file__).resolve().parents[1]
A_PATH = ROOT / "analysis/axon_trace_download.py"
A_SHA = "2e7315543a2f0ed469346348da184e3914ac2d36f8308e9afe90a4bb9f0afd58"
if hashlib.sha256(A_PATH.read_bytes()).hexdigest() != A_SHA: raise RuntimeError("Frozen Stage A utilities changed")
SPEC = importlib.util.spec_from_file_location("stage_b_acquisition_utilities", A_PATH)
A = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(A)
PLAN_SHA = A.PLAN_SHA
SCHEMA = "synarch_stage_b_download_v1"
AUTH_SCHEMA = "synarch_stage_b_authorization_v1"
FIELDS = A.MANIFEST_FIELDS
FROZEN = "data/interim/axon_engine/real/freeze/frozen.json"
COMMITMENT = "data/interim/axon_engine/real/test/PREDICTION_COMMITMENT.json"
UNAVAILABLE = "unavailable_at_source"
METHODS = ("nearest", "medoid", "peptide_nearest", "peptide_medoid", "C1_medoid_copy",
           "C4_uniform_attach_500", "C5_uniform_attach_1000", "F1_conditional_attach_500",
           "F2_conditional_attach_1000", "E1_local_uniform", "E2_local_kernel",
           "E3_local_kernel_label", "E4_local_kernel_regional", "E5_local_kernel_label_regional",
           "E6_local_consensus", "E7_local_kernel_inferred_profile")


def under(base, relative):
    if not isinstance(relative, str) or Path(relative).is_absolute(): raise A.AcquisitionError("Explicit relative artifact path required")
    result = (Path(base) / relative).resolve()
    if not result.is_relative_to(Path(base).resolve()): raise A.AcquisitionError("Artifact escaped its declared root")
    return result


def read(path): return json.loads(Path(path).read_text())
def sha(path): return A.sha_file(path)
def ids_sha(ids): return hashlib.sha256("\n".join(ids).encode()).hexdigest()


def metadata(root, expected_plan_sha=PLAN_SHA):
    root = Path(root); plan_path = root / "reports/axon_dataset_manifest.tsv"; local_path = root / "data/interim/axon_dataset/local_manifest.tsv"
    if sha(plan_path) != expected_plan_sha: raise A.AcquisitionError("Locked split SHA differs")
    rows = A.Downloader._unique(A.read_tsv(plan_path), "locked plan")
    local = A.Downloader._unique(A.read_tsv(local_path), "local manifest")
    if set(rows) != set(local) or any(any(rows[n][k] != local[n][k] for k in ("part", "animal")) for n in rows):
        raise A.AcquisitionError("Plan/local exact identity and split join differs")
    parts = {p: {r["animal"] for r in rows.values() if r["part"] == p} for p in ("train", "validation", "test")}
    if any(parts[a] & parts[b] for a, b in (("train", "validation"), ("train", "test"), ("validation", "test"))): raise A.AcquisitionError("Animal split overlap")
    if expected_plan_sha == PLAN_SHA and (len(rows) != 7180 or len(parts["test"]) != 40 or sum(r["part"] == "test" for r in rows.values()) != 1547):
        raise A.AcquisitionError("Production TEST cohort is not1547neurons/40animals")
    return rows, sha(local_path)


def plan(root=ROOT):
    rows, local_sha = metadata(root); wanted = [r for r in rows.values() if r["part"] == "test"]
    return {"stage": "B", "test_neurons": len(wanted), "test_animals": len({r["animal"] for r in wanted}),
            "plan_sha256": PLAN_SHA, "local_manifest_sha256": local_sha, "requests": 0,
            "requires": "ROOTGO, frozen clean method, committed TEST forecasts and verified copying controls"}


def git_clean(root, paths):
    if subprocess.run(["git", "--no-optional-locks", "-C", str(root), "ls-files", "--error-unmatch", "--", *paths], capture_output=True, timeout=30).returncode:
        raise A.AcquisitionError("Engine/protocol/tests must be committed before TEST acquisition")
    status = subprocess.run(["git", "--no-optional-locks", "-C", str(root), "status", "--porcelain", "--", *paths], capture_output=True, text=True, timeout=30)
    if status.returncode or status.stdout.strip(): raise A.AcquisitionError("Frozen engine/protocol/tests are modified")


def method_contract(methods, family, *, production):
    if len(methods) != len(set(methods)):
        raise A.AcquisitionError("Duplicate forecast methods")
    if production and (set(methods) != set(METHODS) or family != "all"):
        raise A.AcquisitionError("Production Stage B requires the exact approved16methods and all family")


def verify_proofs(root, rows, local_sha, freeze, freeze_sha, commitment, commitment_sha, *, checker=git_clean):
    """Read metadata/source/predicted bytes only; never a measured TEST arbor."""
    root = Path(root); freeze_path, commit_path = under(root, freeze), under(root, commitment)
    if sha(freeze_path) != freeze_sha or sha(commit_path) != commitment_sha: raise A.AcquisitionError("Explicit freeze/commitment hash mismatch")
    f, c = read(freeze_path), read(commit_path); test_dir = commit_path.parent
    manifest_path, copies_path = test_dir / "prediction_manifest.json", test_dir / "COPIES_VERIFIED.json"
    m, copies = read(manifest_path), read(copies_path); wanted = sorted(n for n, r in rows.items() if r["part"] == "test")
    if (f.get("schema") != "axon_engine_freeze_v2" or f.get("world_kind") != "real" or not f.get("validation_method_equals_frozen_method")
        or not f.get("git", {}).get("engine_files_committed_and_clean") or f.get("test_access", {}).get("test_traces_listed_on_disk_at_freeze") != 0
        or f.get("test_access", {}).get("strongest_seal_no_test_trace_on_disk") is not True): raise A.AcquisitionError("Missing real pre-TEST freeze seal")
    pins = f["pins"]; cohort = pins["cohort"]
    if cohort["manifest"] != sha(root / "reports/axon_dataset_manifest.tsv") or cohort["local_manifest"] != local_sha or pins["test_queries"] != len(wanted) or pins["test_queries_sha256"] != ids_sha(wanted):
        raise A.AcquisitionError("Freeze does not bind exact locked TEST IDs")
    if (c.get("schema") != "axon_engine_prediction_commitment_v2" or c.get("stage") != "test" or sha(manifest_path) != c["manifest_sha256"]
        or m.get("schema") != "axon_engine_prediction_manifest_v2" or m.get("stage") != "test" or m.get("part") != "test" or m.get("world_kind") != "real"
        or m.get("queries") != wanted or c.get("queries") != len(wanted)): raise A.AcquisitionError("TEST forecast commitment incomplete or changed")
    expected_pins = {k: v for k, v in pins.items() if k not in ("family", "test_queries", "test_queries_sha256")}
    if m["pins"] != expected_pins or m["code_sha256"] != pins["code_sha256"] or m["family"] != f["family"] or m["frozen_selection"] != {
        "selected": f["selected"], "freeze_sha256": freeze_sha, "descriptive_reference": f["descriptive_reference"]}:
        raise A.AcquisitionError("Forecasts belong to another frozen method")
    methods = [r["name"] for r in m["methods"]]
    method_contract(methods, m["family"], production=cohort["manifest"] == PLAN_SHA)
    if methods != c["methods"]: raise A.AcquisitionError("Declared forecast method matrix differs")
    seen, predicted = set(), []
    for p in m["predictions"]:
        key = p["query_id"], p["model"]
        if key in seen or key[0] not in wanted or key[1] not in methods: raise A.AcquisitionError("Extra/duplicate forecast")
        relative = p["path"]
        if not relative.startswith("predictions/") or not relative.endswith(".npz"): raise A.AcquisitionError("Only committed predictedNPZ byte checks are permitted")
        path = under(test_dir, relative)
        if sha(path) != p["sha256"]: raise A.AcquisitionError("Committed prediction bytes changed")
        predicted.append(path); seen.add(key)
    if seen != {(n, model) for n in wanted for model in methods} or c["predictions"] != len(seen): raise A.AcquisitionError("Missing TEST forecasts")
    if copies.get("stage") != "test" or copies.get("manifest_sha256") != c["manifest_sha256"] or copies.get("all_equal") is not True or copies.get("differences", 0) != 0 or copies.get("queries") != len(wanted) or copies.get("forecasts_compared") != 4 * len(wanted):
        raise A.AcquisitionError("Copy-control equality must be verified before TEST download")
    source_pins = {}
    for section in ("engine", "frozen", "general_layer"):
        source_pins.update(pins["code_sha256"].get(section, {}))
    source_pins["reports/axon_engine_protocol.md"] = pins["protocol_sha256"]
    for relative, expected in source_pins.items():
        if sha(under(root, relative)) != expected: raise A.AcquisitionError("Pinned method source/protocol changed")
    tracked = sorted(set(pins["code_sha256"]["engine"]) | {"reports/axon_engine_protocol.md", "tests/test_axon_engine.py", "tests/test_synarch_axon_bridge.py"})
    checker(root, tracked)
    proof_paths = {"freeze": freeze, "commitment": commitment, "prediction_manifest": str(manifest_path.relative_to(root)), "copies_verified": str(copies_path.relative_to(root))}
    proof_sha = {k: sha(under(root, p)) for k, p in proof_paths.items()}
    watched = [under(root, p) for p in [*proof_paths.values(), *source_pins, *tracked, "reports/axon_dataset_manifest.tsv", "data/interim/axon_dataset/local_manifest.tsv"]] + predicted
    return proof_paths, proof_sha, {str(p): (p.stat().st_size, p.stat().st_mtime_ns, p.stat().st_ctime_ns, p.stat().st_ino) for p in watched}


class Downloader:
    def __init__(self, root=ROOT, *, freeze_sha, commitment_sha, download=False, freeze=FROZEN, commitment=COMMITMENT,
                 expected_plan_sha=PLAN_SHA, opener=None, free_space=None, clock=time.monotonic, sleep=None, checker=git_clean):
        if not download: raise A.AcquisitionError("Explicit ROOTGO/--download required before TEST acquisition")
        self.root = Path(root).resolve(); self.base = self.root / "data/raw/jiao2025/traces"; self.state = self.base / "stage_b"
        self.all_rows, self.local_sha = metadata(self.root, expected_plan_sha); self.plan_sha = expected_plan_sha
        self.proof_paths, self.proof_sha, self.watched = verify_proofs(self.root, self.all_rows, self.local_sha, freeze, freeze_sha, commitment, commitment_sha, checker=checker)
        for path in (Path(__file__).resolve(), A_PATH):
            self.watched[str(path)] = (path.stat().st_size, path.stat().st_mtime_ns, path.stat().st_ctime_ns, path.stat().st_ino)
        self.rows = {n: r for n, r in self.all_rows.items() if r["part"] == "test"}; self.order = sorted(self.rows)
        self.stop, self.lock = threading.Event(), threading.RLock(); self.clock = clock
        self.limiter = A.StartLimiter(self.stop, clock=clock, sleep=sleep)
        self.opener = opener or urllib.request.build_opener(A.NoRedirect()).open
        self.free_space = free_space or (lambda: shutil.disk_usage(self.root).free)
        self.receipts, self.failures = {}, {}; self.run_id = A.utc().replace(":", "").replace("-", "") + "_" + uuid.uuid4().hex[:8]
        self.authorization_sha = None

    def barrier(self):
        for name, expected in self.watched.items():
            p = Path(name); current = (p.stat().st_size, p.stat().st_mtime_ns, p.stat().st_ctime_ns, p.stat().st_ino)
            if current != expected: raise A.StopAcquisition("Frozen authorization/source/prediction metadata changed; resume must reverify full hashes")
        if self.authorization_sha is not None and sha(self.state / "authorization.json") != self.authorization_sha:
            raise A.StopAcquisition("Stage B authorization changed")

    def url(self, neuron):
        self.barrier()
        if neuron not in self.rows: raise A.AcquisitionError("Only exact TEST IDs can reach URL construction")
        return A.SOURCE.format(neuron)

    def target(self, neuron):
        r = self.rows[neuron]
        return self.base / "test" / r["animal"] / (neuron + ".swc")

    def log(self, event, **fields):
        with self.lock:
            self.state.mkdir(parents=True, exist_ok=True)
            with (self.state / "acquisition.jsonl").open("a") as f:
                f.write(A.canonical(dict(run_id=self.run_id, utc=A.utc(), event=event, **fields))); f.flush(); os.fsync(f.fileno())

    def initialize(self):
        snapshot = self.state / "stage_a_manifest.tsv"; common = self.base / "MANIFEST.tsv"
        a_bytes = snapshot.read_bytes() if snapshot.exists() else common.read_bytes()
        if not a_bytes.endswith(b"\n") or not common.read_bytes().startswith(a_bytes): raise A.AcquisitionError("Stage A manifest prefix changed")
        a_rows = list(csv.DictReader(io.StringIO(a_bytes.decode()), delimiter="\t")); wanted = {n for n, r in self.all_rows.items() if r["part"] != "test"}
        if len(a_rows) != len(wanted) or {r["neuron"] for r in a_rows} != wanted: raise A.AcquisitionError("Stage A prefix is not the exact TRAIN/validation cohort")
        for r in a_rows:
            part = self.all_rows[r["neuron"]]["part"]
            if not r["file"].startswith(part + "/"): raise A.AcquisitionError("Stage A path/split differs")
        metadata_paths = [self.base / "summary.json", self.base / "acquisition.jsonl", *sorted((self.base / "runs").glob("*.json"))]
        receipt_sha = {n: sha(self.base / "receipts" / (n + ".json")) for n in sorted(wanted)}
        source_pins = {"analysis/axon_trace_download_stage_b.py": sha(Path(__file__)), "analysis/axon_trace_download.py": A_SHA}
        auth = {"schema": AUTH_SCHEMA, "proof_paths": self.proof_paths, "proof_sha256": self.proof_sha,
                "plan_sha256": self.plan_sha, "local_manifest_sha256": self.local_sha,
                "stage_a_manifest_sha256": hashlib.sha256(a_bytes).hexdigest(), "stage_a_manifest_bytes": len(a_bytes),
                "stage_a_receipts_sha256": receipt_sha,
                "stage_a_metadata_sha256": {str(p.relative_to(self.base)): sha(p) for p in metadata_paths}, "source_sha256": source_pins}
        auth_path = self.state / "authorization.json"
        if auth_path.exists():
            old = read(auth_path)
            if {k: v for k, v in old.items() if k != "created_utc"} != auth: raise A.AcquisitionError("Existing B authorization or unchanged Stage A records differ")
        else:
            if any((self.base / "test").rglob("*.swc")): raise A.AcquisitionError("Unreceipted TEST files predate B authorization")
            self.state.mkdir(parents=True, exist_ok=True); A.atomic_bytes(snapshot, a_bytes)
            A.atomic_bytes(auth_path, A.canonical({**auth, "created_utc": A.utc()}).encode())
        self.a_bytes = a_bytes; self.authorization_sha = sha(auth_path)

    def manifest(self):
        if not (self.base / "MANIFEST.tsv").read_bytes().startswith(self.a_bytes): raise A.AcquisitionError("Stage A prefix changed during B")
        text = io.StringIO(newline=""); writer = csv.DictWriter(text, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        for n in self.order:
            if n in self.receipts: writer.writerow({k: self.receipts[n][k] for k in FIELDS})
        A.atomic_bytes(self.base / "MANIFEST.tsv", self.a_bytes + text.getvalue().encode())

    def resume(self, n):
        path = self.state / "receipts" / (n + ".json")
        if not path.exists(): return False
        r = read(path)
        if r.get("schema") != SCHEMA or r.get("neuron") != n or r.get("authorization_sha256") != self.authorization_sha or r.get("proof_sha256") != self.proof_sha or r.get("plan_sha256") != self.plan_sha or r.get("local_manifest_sha256") != self.local_sha:
            raise A.AcquisitionError("B receipt declaration differs")
        url = self.url(n)
        if r.get("url") != url or r.get("source") != url or r.get("final_url") != url: raise A.AcquisitionError("Receipt source URL differs")
        target = self.target(n)
        if r.get("availability") == UNAVAILABLE:
            if r.get("status") not in (404, 410) or r.get("file") != "" or r.get("frame") != UNAVAILABLE or r.get("sha256") != UNAVAILABLE or r.get("bytes") != 0 or target.exists(): raise A.AcquisitionError("Invalid unavailable receipt")
        elif r.get("availability") == "available":
            if r.get("status") != 200 or not 0 < r.get("bytes", 0) <= A.MAX_FILE or target.is_symlink() or r.get("file") != str(target.relative_to(self.base)) or r.get("frame") != "portal_true_hemisphere" or target.stat().st_size != r.get("bytes") or sha(target) != r.get("sha256"): raise A.AcquisitionError("B raw byte/hash mismatch; no silent redownload")
        else: raise A.AcquisitionError("Unknown receipt availability")
        self.receipts[n] = r; self.log("resume_verified", neuron=n, sha256=r["sha256"]); return True

    def checkpoint(self, n, *, count, checksum, status, headers, unavailable=False):
        self.barrier(); url = self.url(n); reason = f"HTTP{status}: source serves no trace" if unavailable else None
        row = dict(schema=SCHEMA, stage="B", neuron=n, availability=UNAVAILABLE if unavailable else "available",
                   file="" if unavailable else str(self.target(n).relative_to(self.base)), bytes=count, sha256=checksum,
                   frame=UNAVAILABLE if unavailable else "portal_true_hemisphere", source=url, url=url, final_url=url,
                   licence="no portal licence stated; private research acquisition; do not redistribute",
                   how=reason or "Stage B public portal GET; raw bytes unchanged", reason=reason, status=status,
                   fetched_utc=A.utc(), etag=headers.get("ETag"), last_modified=headers.get("Last-Modified"),
                   plan_sha256=self.plan_sha, local_manifest_sha256=self.local_sha,
                   authorization_sha256=self.authorization_sha, proof_sha256=self.proof_sha)
        with self.lock:
            A.atomic_bytes(self.state / "receipts" / (n + ".json"), A.canonical(row).encode()); self.receipts[n] = row; self.manifest()
            self.log("checkpoint", neuron=n, availability=row["availability"], bytes=count, sha256=checksum)

    def stream(self, n, response):
        target = self.target(n); target.parent.mkdir(parents=True, exist_ok=True)
        if target.is_symlink(): raise A.AcquisitionError("Symlink TEST target refused")
        partial = target.with_suffix(".swc.part"); partial.unlink(missing_ok=True)
        expected = response.headers.get("Content-Length"); expected = int(expected) if expected is not None else None
        if expected is not None and not 0 < expected <= A.MAX_FILE: raise A.StopAcquisition("Invalid/oversized Content-Length")
        count, h = 0, hashlib.sha256()
        with partial.open("xb") as f:
            while True:
                if self.stop.is_set(): raise A.StopAcquisition("Interrupted B transfer; partial retained")
                block = response.read(A.CHUNK)
                if not block: break
                if count == 0 and (block.lstrip().lower().startswith((b"<", b"{", b"[")) or b"please login" in block[:4096].lower()):
                    self.stop.set(); raise A.StopAcquisition("HTML/login/JSON response refused")
                count += len(block)
                if count > A.MAX_FILE or self.free_space() < 1_000_000_000: raise A.StopAcquisition("File-size/free-space limit")
                f.write(block); h.update(block)
            f.flush(); os.fsync(f.fileno())
        if count == 0 or (expected is not None and count != expected): raise A.AcquisitionError("Empty/incomplete B file")
        partial.replace(target); A.sync_dir(target.parent); return count, h.hexdigest()

    def fetch(self, n):
        url = self.url(n)
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "synarch-stage-b/1.0", "Accept": "text/plain"})
                def start():
                    self.barrier(); self.log("request_start", neuron=n, url=url, attempt=attempt + 1, monotonic_start=self.clock(), authorization_sha256=self.authorization_sha)
                    try:
                        response = self.opener(request, timeout=120)
                        if response.status != 200:
                            status, headers = response.status, response.headers; response.close()
                            raise urllib.error.HTTPError(url, status, "unexpected status", headers, None)
                        if response.geturl() != url or any(k in response.headers.get("Content-Type", "").lower() for k in ("html", "json")):
                            response.close(); self.stop.set(); raise A.StopAcquisition("Redirect/HTML/login refused")
                        return response
                    except urllib.error.HTTPError as e:
                        if e.code == 429 or 500 <= e.code < 600:
                            try: self.limiter.defer(A.Downloader.retry_delay(e.headers.get("Retry-After") if e.headers else None, attempt))
                            except A.StopAcquisition: self.stop.set(); raise
                        if e.code in (401, 403) or 300 <= e.code < 400: self.stop.set()
                        raise
                    except (ssl.SSLError, ssl.CertificateError): self.stop.set(); raise
                    except urllib.error.URLError as e:
                        if isinstance(e.reason, (ssl.SSLError, ssl.CertificateError)): self.stop.set()
                        raise
                    except A.StopAcquisition: self.stop.set(); raise
                with self.limiter.call(start) as response:
                    if response.status != 200: raise urllib.error.HTTPError(url, response.status, "unexpected status", response.headers, None)
                    if response.geturl() != url or any(k in response.headers.get("Content-Type", "").lower() for k in ("html", "json")): raise A.StopAcquisition("Redirect/HTML/login refused")
                    count, checksum = self.stream(n, response); self.checkpoint(n, count=count, checksum=checksum, status=200, headers=response.headers)
                return
            except urllib.error.HTTPError as e:
                self.log("http_error", neuron=n, status=e.code, attempt=attempt + 1)
                if e.code in (404, 410): self.checkpoint(n, count=0, checksum=UNAVAILABLE, status=e.code, headers=e.headers or {}, unavailable=True); return
                if e.code in (401, 403) or 300 <= e.code < 400: raise A.StopAcquisition("Access/redirect response; no alternate route") from e
                if not (e.code == 429 or 500 <= e.code < 600) or attempt == 2: raise A.AcquisitionError(f"HTTP{e.code}; remains unfinished") from e
            except (ssl.SSLError, ssl.CertificateError) as e: raise A.StopAcquisition("TLS verification failure; remains enabled") from e
            except urllib.error.URLError as e:
                if isinstance(e.reason, (ssl.SSLError, ssl.CertificateError)): raise A.StopAcquisition("TLS verification failure") from e
                raise A.AcquisitionError("Network failure; remains unfinished") from e

    def job(self, n):
        try:
            if self.stop.is_set(): raise A.StopAcquisition("Stopped before request")
            self.fetch(n)
        except Exception as error:
            if isinstance(error, A.StopAcquisition): self.stop.set()
            with self.lock: self.failures[n] = {"type": type(error).__name__, "error": str(error)}
            self.log("file_failed", neuron=n, **self.failures[n])

    def run(self, max_new=None):
        if max_new is not None and (isinstance(max_new, bool) or not isinstance(max_new, int) or max_new < 0): raise ValueError("Nonnegative max-new required")
        self.base.mkdir(parents=True, exist_ok=True)
        with (self.base / ".download.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB); self.barrier(); self.initialize()
            self.log("run_start", stage="B", max_new=max_new, authorization_sha256=self.authorization_sha)
            try:
                if self.free_space() <= A.MIN_FREE: raise A.StopAcquisition("More than12GB free space required")
                pending = []
                for n in self.order:
                    try:
                        if not self.resume(n): pending.append(n)
                    except Exception as error:
                        self.failures[n] = {"type": type(error).__name__, "error": str(error)}
                        self.log("resume_failed", neuron=n, **self.failures[n])
                self.manifest(); selected = pending if max_new is None else pending[:max_new]
                with ThreadPoolExecutor(max_workers=2) as pool:
                    for future in as_completed([pool.submit(self.job, n) for n in selected]): future.result()
                self.barrier()
            except Exception as error:
                self.stop.set(); self.failures["__run__"] = {"type": type(error).__name__, "error": str(error)}
                self.log("run_stopped", **self.failures["__run__"])
            unfinished = [n for n in self.order if n not in self.receipts]
            unavailable = sorted(n for n, r in self.receipts.items() if r["availability"] == UNAVAILABLE)
            result = dict(schema=SCHEMA, stage="B", run_id=self.run_id, utc=A.utc(), plan_sha256=self.plan_sha,
                local_manifest_sha256=self.local_sha, authorization_sha256=self.authorization_sha,
                completed=len(self.receipts), acquired=len(self.receipts) - len(unavailable), unavailable_at_source=len(unavailable),
                unavailable_ids=unavailable, unavailable_by_source_label={label: sum(self.rows[n].get("source_label", "unknown") == label for n in unavailable) for label in sorted({r.get("source_label", "unknown") for r in self.rows.values()})},
                unavailable_by_animal={animal: sum(self.rows[n]["animal"] == animal for n in unavailable) for animal in sorted({r["animal"] for r in self.rows.values()})},
                expected_test_neurons=len(self.rows), expected_test_animals=len({r["animal"] for r in self.rows.values()}),
                completed_by_part={"test": len(self.receipts)}, unfinished_ids=unfinished, failures=self.failures,
                status="failed" if self.failures else "partial" if unfinished else "complete")
            A.atomic_bytes(self.state / "runs" / (self.run_id + ".json"), A.canonical(result).encode())
            A.atomic_bytes(self.state / "summary.json", A.canonical(result).encode())
            self.log("run_complete", status=result["status"], completed=result["completed"], unavailable=len(unavailable), failed=len(self.failures))
            return result


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--root", type=Path, default=ROOT)
    action = parser.add_mutually_exclusive_group(required=True); action.add_argument("--plan", action="store_true"); action.add_argument("--download", action="store_true")
    parser.add_argument("--freeze", default=FROZEN); parser.add_argument("--commitment", default=COMMITMENT)
    parser.add_argument("--freeze-sha"); parser.add_argument("--commitment-sha"); parser.add_argument("--max-new", type=int)
    args = parser.parse_args()
    if args.plan: result = plan(args.root)
    else:
        if not args.freeze_sha or not args.commitment_sha: parser.error("--freeze-sha and --commitment-sha are required for --download")
        result = Downloader(args.root, freeze=args.freeze, commitment=args.commitment, freeze_sha=args.freeze_sha,
                            commitment_sha=args.commitment_sha, download=True).run(args.max_new)
    print(A.canonical({k: v for k, v in result.items() if k not in ("unfinished_ids", "failures", "unavailable_ids")}), end="")
    return 2 if result.get("status") == "failed" else 0


if __name__ == "__main__": raise SystemExit(main())
