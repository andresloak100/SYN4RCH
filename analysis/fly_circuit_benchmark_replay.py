"""Independently replay saved fly scores; never fit a model or change forecasts.

Uses raw table filtering and scikit-learn metrics, independently of the
benchmark's graph loader and metric implementation. One-brain development only.
"""
from pathlib import Path
import hashlib
import json

import numpy as np
import pyarrow.parquet as pq
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data/interim/fly_circuit_benchmark"


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def replay():
    frozen = read(BASE / "FROZEN.json")
    protocol = read(ROOT / "docs/fly_circuit_benchmark_protocol.json")
    commitment = read(BASE / "PREDICTIONS_COMMITTED.json")
    result = read(ROOT / "reports/fly_circuit_benchmark.json")
    complete = read(BASE / "SCORE_COMPLETE.json")
    assert sha(ROOT / "reports/fly_circuit_benchmark.json") == complete["report_sha256"]
    assert sha(BASE / "FROZEN.json") == commitment["freeze_sha256"] == result["freeze_sha256"]
    assert sha(BASE / "PREDICTIONS_COMMITTED.json") == result["prediction_commitment_sha256"]
    for path, expected in frozen["pins"].items():
        assert sha(ROOT / path) == expected, path

    all_ids = np.asarray(frozen["node_ids"], dtype=np.int64)
    sources = all_ids[frozen["query_indices"]]
    targets = all_ids[frozen["train_indices"]]
    assert (np.diff(sources) > 0).all() and (np.diff(targets) > 0).all()
    assert not np.intersect1d(sources, targets).size
    truth = np.zeros((len(sources), len(targets)), dtype=bool)
    admitted_rows = 0
    columns = ["Presynaptic_ID", "Postsynaptic_ID", "Connectivity"]
    for batch in pq.ParquetFile(ROOT / protocol["edges_path"]).iter_batches(
            batch_size=131072, columns=columns):
        pre, post = (batch.column(i).to_numpy() for i in (0, 1))
        a, b = np.searchsorted(sources, pre), np.searchsorted(targets, post)
        aa, bb = np.minimum(a, len(sources) - 1), np.minimum(b, len(targets) - 1)
        take = ((a < len(sources)) & (b < len(targets))
                & (sources[aa] == pre) & (targets[bb] == post))
        counts = batch.column(2).to_numpy()[take]
        assert np.isfinite(counts).all() and (counts > 0).all()
        assert (counts == np.floor(counts)).all()
        truth[a[take], b[take]] = True
        admitted_rows += int(take.sum())

    scores, differences = {}, {}
    for name, pin in commitment["predictions"].items():
        path = ROOT / pin["path"]
        assert sha(path) == pin["sha256"], name
        p = np.load(path, allow_pickle=False)
        assert p.shape == truth.shape and np.isfinite(p).all()
        assert ((p >= 0) & (p <= 1)).all()
        y, pp = truth.ravel(), p.ravel()
        scores[name] = {
            "average_precision": float(average_precision_score(y, pp)),
            "log_loss": float(log_loss(y, np.clip(pp, 1e-8, 1 - 1e-8), labels=[0, 1])),
            "brier": float(brier_score_loss(y, pp)),
        }
        differences[name] = {}
        for key, value in scores[name].items():
            difference = abs(value - result["scores"][name][key])
            assert difference < 1e-12, (name, key, difference)
            differences[name][key] = difference
        assert result["scores"][name]["pairs"] == truth.size
        assert result["scores"][name]["positives"] == int(truth.sum())

    candidate = scores["shared_input_logistic"]
    gates = {
        name: (candidate["average_precision"] > values["average_precision"]
               and candidate["log_loss"] < values["log_loss"]
               and candidate["brier"] < values["brier"])
        for name, values in scores.items() if name != "shared_input_logistic"
    }
    assert gates == result["candidate_beats_control_on_AP_logloss_Brier"]
    assert all(gates.values()) == result["beats_every_control"]
    report = {
        "schema": "synarch_fly_circuit_independent_replay_v1", "status": "PASS",
        "source_report_sha256": sha(ROOT / "reports/fly_circuit_benchmark.json"),
        "source_count_table_sha256": sha(ROOT / protocol["edges_path"]),
        "query_sources": len(sources), "targets": len(targets), "pairs": int(truth.size),
        "positive_pairs": int(truth.sum()), "raw_rows_admitted": admitted_rows,
        "independent_metrics": scores, "absolute_metric_differences": differences,
        "candidate_beats_every_control": all(gates.values()),
        "method": "Independent raw Arrow-table filtering and scikit-learn AP/log-loss/Brier; no model fit or forecast regeneration.",
        "scope": "Single-fly masked graph-completion development. AP is not percentage accuracy. No mouse improvement, neural-control efficacy or longevity benefit established.",
    }
    output = ROOT / "reports/fly_circuit_benchmark_replay.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return {"status": report["status"], "pairs": report["pairs"],
            "positive_pairs": report["positive_pairs"], "all_control_gate": all(gates.values()),
            "report_sha256": sha(output)}


if __name__ == "__main__":
    print(json.dumps(replay(), sort_keys=True))
