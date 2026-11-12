"""Structural TRAIN-only check of the frozen completion mask; no model scoring."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CONTRACT_SHA = "13e6cd60b2021e0d673407ec7aaf5906a3652b41310a11a41053148f9c41c8bd"
REAL = ROOT / "data/interim/axon_engine/real"
OUT = ROOT / "data/interim/axon_completion_v1_root"
_CONTRACT = None


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def contract():
    global _CONTRACT
    if _CONTRACT is None:
        path = ROOT / "analysis/axon_completion_contract.py"
        if sha(path) != CONTRACT_SHA:
            raise RuntimeError("Completion contract changed before structural preflight")
        spec = importlib.util.spec_from_file_location("completion_root_preflight_contract", path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        _CONTRACT = module
    return _CONTRACT


def check_one(item):
    nid, info = item
    path = REAL / "library" / f"{nid}.npz"
    if sha(path) != info["canonical_file_sha256"]:
        raise RuntimeError(f"Changed canonical TRAIN arbor: {nid}")
    with np.load(path, allow_pickle=False) as z:
        arbor = SimpleNamespace(xyz=z["xyz_um"], parent=z["parent"])
    c = contract()
    obs, hidden = c.split_arbor(arbor)
    original = float(np.linalg.norm(arbor.xyz[1:] - arbor.xyz[arbor.parent[1:]], axis=1).sum())
    observed = float(np.linalg.norm(obs.xyz[1:] - obs.xyz[obs.parent[1:]], axis=1).sum())
    missing = hidden.forest.cable()
    conservation_error = abs(observed + missing - original)
    if conservation_error > 1e-7 + original * 1e-12:
        raise RuntimeError(f"Mask loses or adds cable: {nid}")
    points, mass = c.sample_forest(hidden.forest)
    mass_error = abs(float(mass.sum()) - missing)
    if mass_error > 1e-7 + missing * 1e-12:
        raise RuntimeError(f"Suffix quadrature loses or adds cable: {nid}")
    frontiers = len(obs.frontier_indices())
    if hidden.eligible and not frontiers:
        raise RuntimeError(f"Missing cable has no observed boundary: {nid}")
    return {"neuron": nid, "brain_sample": info["animal"], "eligible": bool(hidden.eligible),
            "observed_nodes": len(obs.xyz), "frontiers": frontiers,
            "original_cable_um": original, "observed_cable_um": observed,
            "hidden_cable_um": missing, "observed_cable_fraction": observed / original,
            "observation_sha256": obs.sha256, "sampled_hidden_points": len(points),
            "conservation_error_um": conservation_error, "quadrature_mass_error_um": mass_error,
            "registration_offset_um": float(np.linalg.norm(np.asarray(info["table_soma_um"]) - obs.xyz[0]))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=2, choices=(1, 2))
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    output = OUT / "structural_preflight.json"
    if output.exists():
        raise RuntimeError("Completed structural receipt exists; preserve it")
    contract()
    ingestion = REAL / "ingestion_log.json"
    data = json.loads(ingestion.read_text())["neurons"]
    items = sorted((n, d) for n, d in data.items() if d["part"] == "train" and d["status"] == "included")
    if len(items) != 4368 or len({d["animal"] for _, d in items}) != 101:
        raise RuntimeError("Unexpected original TRAIN cohort")
    started = time.monotonic()
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for i, row in enumerate(executor.map(check_one, items, chunksize=8), 1):
            rows.append(row)
            if i % 500 == 0:
                print(f"Structural mask check: {i}/{len(items)}", flush=True)
    rows.sort(key=lambda r: r["neuron"])
    cohort = OUT / "structural_cohort.csv"
    with cohort.open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    eligible = [r for r in rows if r["eligible"]]
    receipt = {"schema": "axon_completion_root_structural_preflight_v1",
               "utc": datetime.now(timezone.utc).isoformat(), "passed": True,
               "scope": "Original TRAIN only. Structural masking, eligibility and cable-mass checks; no forecasts or model scores.",
               "contract_sha256": CONTRACT_SHA, "script_sha256": sha(__file__),
               "ingestion_sha256": sha(ingestion), "cohort_csv_sha256": sha(cohort),
               "neurons": len(rows), "brain_samples": 101, "eligible_neurons": len(eligible),
               "eligible_brain_samples": len({r["brain_sample"] for r in eligible}),
               "excluded_fully_observed": len(rows) - len(eligible),
               "excluded_by_sample": dict(Counter(r["brain_sample"] for r in rows if not r["eligible"])),
               "frontiers_histogram": dict(Counter(r["frontiers"] for r in rows)),
               "observed_cable_fraction_percentiles_0_50_90_100": np.percentile([r["observed_cable_fraction"] for r in eligible], [0, 50, 90, 100]).tolist(),
               "max_conservation_error_um": max(r["conservation_error_um"] for r in rows),
               "max_quadrature_mass_error_um": max(r["quadrature_mass_error_um"] for r in rows),
               "seconds": time.monotonic() - started}
    output.write_text(json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False) + "\n")
    print(json.dumps(receipt, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
