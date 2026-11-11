"""Verify the trace download listing before the engine consumes it. Read-only; hashes bytes, never parses a trace.

    python3 -I analysis/axon_trace_verify.py --stage A [--allow-partial] [--out FILE]
    python3 -I analysis/axon_trace_verify.py --stage B [--allow-partial] [--out FILE]

The download belongs to the other agent (analysis/axon_trace_download.py, with the listing, receipts and
logs under data/raw/jiao2025/traces/). This script changes nothing there. It checks the listing against
the locked cohort (reports/axon_dataset_manifest.tsv, pinned by the plan) and against the bytes on disk:

  * listing columns and frames; every path is <part>/<brain sample>/<neuron>.swc under the neuron's locked part;
  * stage A: every training and validation neuron is listed (as a file or as unavailable at the source) and
    NO test neuron, test path or test address appears anywhere (listing, receipts, log, directory tree);
    stage B: the same, plus every test neuron listed;
  * every listed file exists with the listed size and SHA-256, and does not start like a web page;
  * every file row has a receipt that agrees with it and with the cohort manifest's SHA-256;
  * the downloader's own summary says the run is complete, without failures or unfinished neurons.

Exit code 0 only if everything holds (with --allow-partial, incompleteness alone is reported, not failed).
Stdlib only. No trace content is read beyond its first non-blank byte.
"""
import argparse
import hashlib
import json
import os
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "reports/axon_dataset_plan.md"
COHORT = ROOT / "reports/axon_dataset_manifest.tsv"
TRACES = ROOT / "data/raw/jiao2025/traces"
COLUMNS = ["file", "neuron", "bytes", "sha256", "source", "frame"]
FILE_FRAMES = ("portal_true_hemisphere", "zenodo_mirrored_left")
UNAVAILABLE = "unavailable_at_source"
PARTS = {"A": ("train", "validation"), "B": ("train", "validation", "test")}


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def first_byte(path):
    with open(path, "rb") as f:
        head = f.read(256).lstrip()
    return head[:1]


def verify(stage, traces=TRACES, cohort=COHORT, plan=PLAN, allow_partial=False):
    traces, problems, notes = Path(traces), [], []
    cohort_sha = sha(cohort)
    if cohort_sha not in Path(plan).read_text():
        problems.append("the cohort manifest is not the one pinned by the locked plan")
    rows = [line.rstrip("\n").split("\t") for line in open(cohort)]
    head = rows[0]
    locked = {r[head.index("neuron")]: (r[head.index("animal")], r[head.index("part")]) for r in rows[1:]}
    allowed, sealed = set(PARTS[stage]), {"train", "validation", "test"} - set(PARTS[stage])
    sealed_animals = {a for a, p in locked.values() if p in sealed}

    listing = traces / "MANIFEST.tsv"
    table = [line.rstrip("\n").split("\t") for line in open(listing)]
    columns = table[0]
    if any(c not in columns for c in COLUMNS):
        problems.append("listing lacks columns: " + ", ".join(c for c in COLUMNS if c not in columns))
        return {"ok": False, "problems": problems}
    listed, frames, checked = {}, Counter(), 0
    for r in table[1:]:
        row = dict(zip(columns, r))
        neuron, frame, file = row["neuron"], row["frame"], row["file"]
        if neuron not in locked:
            problems.append("%s: not a neuron of the cohort" % neuron)
            continue
        animal, part = locked[neuron]
        frames[part, frame] += 1
        if part not in allowed:
            problems.append("%s: a %s neuron is listed in stage %s" % (neuron, part, stage))
            continue
        listed.setdefault(neuron, []).append(row)
        if frame == UNAVAILABLE:
            continue
        if frame not in FILE_FRAMES:
            problems.append("%s: unknown frame %s" % (neuron, frame))
        if file != "%s/%s/%s.swc" % (part, animal, neuron):
            problems.append("%s: stored as %s, expected %s/%s/%s.swc" % (neuron, file, part, animal, neuron))
            continue
        path = traces / file
        if not path.is_file():
            problems.append("%s: listed file is missing" % file)
            continue
        if path.stat().st_size != int(row["bytes"]) or sha(path) != row["sha256"]:
            problems.append("%s: size or SHA-256 differs from the listing" % file)
        if first_byte(path) in (b"<", b"{", b""):
            problems.append("%s: does not start like a trace (web page, JSON or empty)" % file)
        receipt = traces / "receipts" / (neuron + ".json")
        if not receipt.is_file():
            problems.append("%s: no receipt" % neuron)
        else:
            got = json.loads(receipt.read_text())
            for key in ("file", "sha256", "frame", "source"):
                if str(got.get(key)) != row[key]:
                    problems.append("%s: receipt disagrees with the listing on %s" % (neuron, key))
            if int(got.get("bytes", -1)) != int(row["bytes"]):
                problems.append("%s: receipt disagrees with the listing on bytes" % neuron)
            if got.get("plan_sha256") != cohort_sha:
                problems.append("%s: receipt was made for another cohort manifest" % neuron)
            if got.get("origin") == "portal" and got.get("status") != 200:
                problems.append("%s: portal receipt without HTTP 200" % neuron)
        checked += 1

    need = {n for n, (_, p) in locked.items() if p in allowed}
    missing = sorted(need - set(listed))
    per_part = {p: {"neurons": sum(1 for n in need if locked[n][1] == p),
                    "with_file": sum(1 for n in listed if locked[n][1] == p and any(r["frame"] != UNAVAILABLE for r in listed[n])),
                    "unavailable_at_source": sum(1 for n in listed if locked[n][1] == p
                                                 and all(r["frame"] == UNAVAILABLE for r in listed[n]))}
                for p in PARTS[stage]}
    if missing:
        (notes if allow_partial else problems).append("%d neurons of stage %s are not listed yet (first: %s)"
                                                      % (len(missing), stage, ", ".join(missing[:3])))

    # nothing of a sealed part anywhere: directory tree, receipts, log
    for part in sealed:
        if (traces / part).exists():
            problems.append("a directory for the sealed part '%s' exists" % part)
    stray = []
    for base, _, files in os.walk(traces):
        for name in files:
            stem = name.split(".")[0]
            if stem in locked and locked[stem][1] in sealed:
                stray.append(os.path.join(base, name))
            elif name.endswith(".swc") and (stem not in locked or not any(r["file"].endswith("/" + name) for r in listed.get(stem, []))):
                stray.append(os.path.join(base, name))
    if stray:
        problems.append("%d files of sealed neurons or unlisted traces on disk (first: %s)" % (len(stray), stray[0]))
    log = traces / "acquisition.jsonl"
    requested_sealed = 0
    if log.exists():
        for line in open(log):
            try:
                event = json.loads(line)
            except ValueError:
                continue
            neuron = str(event.get("neuron", ""))
            if (neuron in locked and locked[neuron][1] in sealed) or neuron.split("_")[0] in sealed_animals:
                requested_sealed += 1
    if requested_sealed:
        problems.append("%d log events name a neuron of a sealed part" % requested_sealed)

    summary = {}
    try:
        summary = json.loads((traces / "summary.json").read_text())
    except (OSError, ValueError):
        problems.append("the downloader's summary.json is missing or unreadable")
    complete = (summary.get("status") == "complete" and not summary.get("failures") and not summary.get("unfinished_ids"))
    if summary and summary.get("plan_sha256") != cohort_sha:
        problems.append("the downloader's summary was made for another cohort manifest")
    if not complete:
        (notes if allow_partial else problems).append("the downloader's summary does not say complete (status %s, %d failures, "
                                                      "%d unfinished)" % (summary.get("status"), len(summary.get("failures") or {}),
                                                                          len(summary.get("unfinished_ids") or [])))
    downloader = ROOT / "analysis/axon_trace_download.py"
    return {"schema": "synarch_trace_verification_v1", "stage": stage, "ok": not problems, "problems": problems, "notes": notes,
            "listing_sha256": sha(listing), "listing_rows": len(table) - 1, "files_hashed": checked, "per_part": per_part,
            "rows_by_part_and_frame": {"%s/%s" % k: v for k, v in sorted(frames.items())},
            "cohort_manifest_sha256": cohort_sha, "downloader_sha256": sha(downloader) if downloader.exists() else None,
            "summary_status": summary.get("status"), "summary_run": summary.get("run_id"), "allow_partial": allow_partial,
            "integrity_scope": "local receipts and hashes show the files are the bytes that were fetched; the source publishes "
                               "no per-file checksum, so this is not a check against an upstream reference"}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("A", "B"), required=True)
    ap.add_argument("--allow-partial", action="store_true")
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    report = verify(args.stage, allow_partial=args.allow_partial)
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("stage %s: %s | listing rows %d, files hashed %d | %s" % (
        args.stage, "OK" if report["ok"] else "PROBLEMS", report.get("listing_rows", 0), report.get("files_hashed", 0),
        json.dumps(report.get("per_part"))))
    for line in report["problems"][:20]:
        print("PROBLEM:", line)
    for line in report.get("notes", []):
        print("note:", line)
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
