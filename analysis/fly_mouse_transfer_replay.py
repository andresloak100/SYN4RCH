"""Independent CSV/metric replay of the saved cortical transfer development run.

No model fitting, forecast changes or new validation cohort is performed.
"""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data/interim/fly_mouse_transfer"


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
    commitment = read(BASE / "PREDICTIONS_COMMITTED.json")
    protocol = read(ROOT / "docs/fly_mouse_transfer_protocol.json")
    result = read(ROOT / "reports/fly_mouse_transfer.json")
    completed = read(BASE / "SCORE_COMPLETE.json")
    assert sha(ROOT / "reports/fly_mouse_transfer.json") == completed["report_sha256"]
    assert sha(BASE / "FROZEN.json") == commitment["freeze_sha256"] == result["freeze_sha256"]
    assert sha(BASE / "PREDICTIONS_COMMITTED.json") == result["prediction_commitment_sha256"]
    for path, expected in frozen["pins"].items():
        assert sha(ROOT / path) == expected, path
    ids = frozen["node_ids"]
    source_ids = {str(ids[i]): k for k, i in enumerate(frozen["query_indices"])}
    target_ids = {str(ids[i]): k for k, i in enumerate(frozen["train_indices"])}
    assert not set(source_ids) & set(target_ids)
    y = np.zeros((len(source_ids), len(target_ids)), dtype=bool)
    participants, self_rows, rows = set(), 0, 0
    with (ROOT / protocol["synapse_path"]).open() as stream:
        for row in csv.DictReader(stream):
            pre, post = row["pre_root_id"], row["post_root_id"]
            assert pre in ids and post in ids
            participants.update((pre, post))
            rows += 1
            self_rows += pre == post
            if pre in source_ids and post in target_ids:
                y[source_ids[pre], target_ids[post]] = True
    assert len(ids) == 362 and len(set(ids) - participants) == 28
    assert rows == 1961 and self_rows == 2
    scores, differences = {}, {}
    for name, pin in commitment["predictions"].items():
        file = ROOT / pin["path"]
        assert sha(file) == pin["sha256"], name
        p = np.load(file, allow_pickle=False)
        assert p.shape == y.shape and np.isfinite(p).all()
        assert ((p >= 0) & (p <= 1)).all()
        yy, pp = y.ravel(), p.ravel()
        values = {
            "average_precision": float(average_precision_score(yy, pp)),
            "log_loss": float(log_loss(yy, np.clip(pp, 1e-8, 1 - 1e-8), labels=[0, 1])),
            "brier": float(brier_score_loss(yy, pp)),
        }
        scores[name] = values
        differences[name] = {}
        for key, value in values.items():
            difference = abs(value - result["scores"][name][key])
            assert difference < 1e-12, (name, key, difference)
            differences[name][key] = difference
        assert result["scores"][name]["pairs"] == y.size
        assert result["scores"][name]["positives"] == int(y.sum())
    controls = set(scores) - {"frozen_fly_residual", "fly_informed_prior"}
    verdicts = {}
    for transfer in ("frozen_fly_residual", "fly_informed_prior"):
        a = scores[transfer]
        verdicts[transfer] = {
            control: (a["average_precision"] > scores[control]["average_precision"]
                      and a["log_loss"] < scores[control]["log_loss"]
                      and a["brier"] < scores[control]["brier"])
            for control in sorted(controls)
        }
    assert verdicts == result["transfer_beats_each_mouse_control_on_AP_LL_Brier"]
    report = {
        "schema": "synarch_fly_mouse_transfer_independent_replay_v1", "status": "PASS",
        "result_sha256": sha(ROOT / "reports/fly_mouse_transfer.json"),
        "nodes": len(ids), "isolates_retained": len(set(ids) - participants),
        "source_synapse_rows": rows, "excluded_autapses": self_rows,
        "queries": len(source_ids), "targets": len(target_ids), "pairs": int(y.size),
        "positive_pairs": int(y.sum()), "independent_metrics": scores,
        "absolute_metric_differences": differences,
        "transfer_beats_each_mouse_control_on_all_metrics": verdicts,
        "no_fit_or_forecast_regeneration": True,
        "scope": "Nonblind development on one visual-cortex mouse. Frozen fly residual performs worse; the fly-prior gain is tiny and descriptive. No hypothalamic, independent-animal or longevity result.",
    }
    output = ROOT / "reports/fly_mouse_transfer_replay.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return {"status": report["status"], "pairs": int(y.size), "positive_pairs": int(y.sum()),
            "report_sha256": sha(output)}


if __name__ == "__main__":
    print(json.dumps(replay(), sort_keys=True))
