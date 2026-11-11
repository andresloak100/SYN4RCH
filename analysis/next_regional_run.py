"""Run-once commands of the DRAFT protocol reports/next_regional_protocol.md: query loader, predict, score.

NOT TO BE RUN ON HELD-OUT DATA BEFORE THE FREEZE. Until the protocol is frozen the only real-data command
that may run is `check-loader` (training rows only). Every other command is exercised on synthetic worlds by
tests/test_next_regional_run.py.

    python3 -I analysis/next_regional_run.py check-loader [--out data/interim/next_regional]
    python3 -I analysis/next_regional_run.py select --out DIR                     (freeze step; training only)
    python3 -I analysis/next_regional_run.py predict --part test --selection DIR/selection.json --out DIR
    python3 -I analysis/next_regional_run.py score --part test --forecasts DIR --commitment SHA \
        --selected METHOD|none --after-engine-score FILE --out DIR

The barrier between forecasts and outcomes is mechanical:
  * `load_query_inputs` reads the manifest and the allow-listed soma-table columns (file, neurons, geno,
    location, AP, DV, ML) of ONE part; rows of every other part are dropped unparsed. It never opens the
    regional table.
  * `predict` reads training outcomes only through the general layer's sealed loader (validation and test
    lines of the regional table are dropped before any value is parsed) and the query inputs above. The
    hidden-outcome canary (tests) shows the forecasts are bit-identical when the query samples' regional
    rows are perturbed, deleted, or unparseable.
  * `score` refuses to run unless the forecast manifest hashes to the recorded commitment, every forecast
    file and every code file matches the manifest, the current evaluation's test score file exists (real
    cohort), and no scoring was started before. Only then does it read the regional rows of the scored part,
    and only of the committed query neurons.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import io
import json
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/interim/next_regional"
QUERY_COLUMNS = ("file", "neurons", "geno", "location", "AP", "DV", "ML")   # never 'type' or 'name'
CODE_FILES = ("analysis/next_regional_core.py", "analysis/next_regional_run.py")
EXPOSURE_LABEL = ("held-out samples of the locked split; their regional outcomes were modelled in aggregate before the "
                  "split was locked (brm_real) and summarised over all samples in the published paper; not pristine")
VIEWER_EVIDENCE_LEVEL = "held_out_previously_modelled_in_aggregate"   # not in the current viewer: export refuses
CONTROL_METHODS = ("MEAN", "IDENTITY", "KNN", "RIDGE", "KNN@crps", "RIDGE@crps", "KNN-CONF")
CANDIDATES = {}        # name -> callable(ctrl_table, query) -> (point, q, p); N1 and N2 are registered before the freeze


class RunRefused(RuntimeError):
    pass


def _load(name, path):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


nc = _load("next_regional_core", ROOT / "analysis/next_regional_core.py")


def general_core():
    return _load("next_regional_synarch_recon_ro", ROOT / "analysis/synarch_recon.py")


# ----------------------------------------------------------------------------------------------
# query inputs (never outcomes)
# ----------------------------------------------------------------------------------------------


def load_query_inputs(world, part):
    """Permitted inputs of the neurons of one part: manifest (neuron, sample, source label, soma label) and the
    allow-listed soma-table columns. Rows of other parts are dropped before any of their cells is parsed."""
    core = general_core()
    audit = {"part": part}
    audit["manifest_sha256"] = nc.sha_file(world.manifest)
    if world.pin("manifest") and audit["manifest_sha256"] != world.pin("manifest"):
        raise RunRefused("the cohort manifest differs from its pinned SHA-256")
    rows = {}
    with open(world.manifest, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["part"] != part:
                continue
            rows[row["neuron"]] = {"sample": row["animal"], "label": row["source_label"], "soma_label": row["soma_label"]}
    audit["soma_csv_sha256"] = nc.sha_file(world.soma_csv)
    if world.pin("soma_csv") and audit["soma_csv_sha256"] != world.pin("soma_csv"):
        raise RunRefused("the soma table differs from its pinned SHA-256")
    alias = dict(world.label_alias)
    found, dropped = {}, 0
    with open(world.soma_csv, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        index = {name: header.index(name) for name in QUERY_COLUMNS}
        for cells in reader:
            neuron = cells[index["neurons"]]
            if neuron not in rows:
                dropped += 1
                continue
            if neuron in found:
                raise RunRefused(f"duplicate soma row {neuron}")
            geno = cells[index["geno"]]
            if alias.get(geno, geno) != rows[neuron]["label"]:
                raise RunRefused(f"table label and manifest label disagree for {neuron}")
            soma = tuple(float(cells[index[a]]) for a in ("AP", "DV", "ML"))
            if not all(np.isfinite(soma)) or soma[2] > 5700.0 + 1e-6:
                raise RunRefused(f"soma of {neuron} is not finite or not in the mirrored frame")
            found[neuron] = {"full_id": cells[index["file"]], "geno": geno, "soma": soma}
    missing = sorted(set(rows) - set(found))
    if missing:
        raise RunRefused(f"{len(missing)} query neurons have no soma row: {missing[:5]}")
    order = sorted(rows)
    query = {"neuron": np.array(order), "full_id": np.array([found[n]["full_id"] for n in order]),
             "sample": np.array([rows[n]["sample"] for n in order]),
             "label": np.array([found[n]["geno"] for n in order]),      # the table's spelling, as in training
             "soma_label": np.array([rows[n]["soma_label"] for n in order]),
             "soma": np.array([found[n]["soma"] for n in order], dtype=float)}
    audit.update({"neurons": len(order), "samples": len(set(query["sample"].tolist())), "other_rows_dropped_unparsed": dropped,
                  "query_ids_sha256": core.ids_sha(order)})
    return query, audit


def check_loader(world, out_dir):
    """Training rows only: the query loader, run on part 'train', must agree with the sealed loader."""
    core = general_core()
    table = core.load_training_table(world, with_outcomes=False)
    query, audit = load_query_inputs(world, "train")
    same = (np.array_equal(query["neuron"], np.asarray(table.neuron)) and np.array_equal(query["soma"], table.soma)
            and np.array_equal(query["label"], np.asarray(table.label).astype(str))
            and np.array_equal(query["sample"], np.asarray(table.sample).astype(str)))
    result = {"agrees_with_sealed_loader": bool(same), "audit": audit, "sealed_audit": table.audit,
              "code_sha256": {f: nc.sha_file(ROOT / f) for f in CODE_FILES}}
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    (Path(out_dir) / "loader_check.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n")
    if not same:
        raise RunRefused("the query loader disagrees with the sealed loader on training rows")
    return result


# ----------------------------------------------------------------------------------------------
# select (freeze step) and predict
# ----------------------------------------------------------------------------------------------


def training_arrays(world):
    core = general_core()
    table = core.load_training_table(world)
    return (table, np.asarray(table.soma, float), np.asarray(table.label).astype(str),
            np.asarray(table.sample).astype(str), table.profile())


def select(world, out_dir, salt="SYNARCH-next-regional-2026-11-11", folds_n=5):
    """Control settings for the test forecasts: whole-sample CV over the 5 outer folds of all training samples."""
    core = general_core()
    table, soma, label, group, Y = training_arrays(world)
    samples = sorted(set(group.tolist()))
    label_of = {s: label[group == s][0] for s in samples}
    folds = core.nested_folds(samples, label_of, outer=folds_n, inner=3, salt=salt)["outer"]
    chosen = nc.choose_settings(soma, label, group, Y, folds)
    record = {"schema": "next_regional_selection_v1", "salt": salt, "folds": folds, "controls": chosen,
              "candidates": "registered before the freeze (N1, N2); none now",
              "training_samples_sha256": table.samples_sha256,
              "code_sha256": {f: nc.sha_file(ROOT / f) for f in CODE_FILES}}
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    path = Path(out_dir) / "selection.json"
    with open(path, "x") as handle:
        handle.write(json.dumps(record, indent=1, sort_keys=True) + "\n")
    return record


def predict(world, part, selection, out_dir, *, allow_train_part=False):
    """Sealed forecasts of every method for the neurons of `part`. Writes forecasts/<method>.npz and
    manifest.json (no timestamps); returns the manifest SHA-256 (the commitment)."""
    if part == "train" and not allow_train_part:
        raise RunRefused("predict is for held-out parts")
    out_dir = Path(out_dir)
    if (out_dir / "manifest.json").exists():
        raise RunRefused("forecasts already exist here; a commitment is never rewritten")
    table, soma, label, group, Y = training_arrays(world)
    query, qaudit = load_query_inputs(world, part)
    overlap = set(query["sample"].tolist()) & set(group.tolist())
    if overlap and not allow_train_part:
        raise RunRefused(f"query samples overlap the training samples: {sorted(overlap)[:5]}")
    ctrl = nc.Controls(soma, label, group, Y)
    (out_dir / "forecasts").mkdir(parents=True, exist_ok=True)
    files = {}
    methods = {m: selection["controls"][m]["setting"] for m in CONTROL_METHODS}
    for name in sorted(methods) + sorted(CANDIDATES):
        if name in CANDIDATES:
            point, q, p = CANDIDATES[name](ctrl, query)
        else:
            point, q, p = nc.method_forecast(ctrl, name, methods[name], query["soma"], query["label"])
        buf = io.BytesIO()
        np.savez_compressed(buf, point=np.asarray(point, np.float64), q=np.asarray(q, np.float32),
                            p_nonzero=np.asarray(p, np.float64))
        path = out_dir / "forecasts" / f"{name}.npz"
        path.write_bytes(buf.getvalue())
        files[name] = nc.sha_file(path)
    manifest = {"schema": "next_regional_forecasts_v1", "part": part, "regions": list(table.regions),
                "query_neurons": query["neuron"].tolist(), "query_full_ids": query["full_id"].tolist(),
                "query_samples": query["sample"].tolist(), "inputs_used": ["soma AP/DV/ML", "source label"],
                "query_audit": qaudit, "training_audit": table.audit, "training_samples_sha256": table.samples_sha256,
                "settings": methods, "candidates": sorted(CANDIDATES), "forecast_sha256": files,
                "quantile_levels": "(i - 0.5) / 199, i = 1..199",
                "code_sha256": {f: nc.sha_file(ROOT / f) for f in CODE_FILES},
                "loader_sha256": nc.sha_file(ROOT / "analysis/synarch_recon.py")}
    text = json.dumps(manifest, indent=1, sort_keys=True) + "\n"
    with open(out_dir / "manifest.json", "x") as handle:
        handle.write(text)
    return nc.sha_file(out_dir / "manifest.json")


# ----------------------------------------------------------------------------------------------
# score (after the commitment and after the current evaluation's test score)
# ----------------------------------------------------------------------------------------------


def read_regional_rows(world, full_ids, regions):
    """Regional rows of exactly these published full IDs; every other line is dropped unparsed."""
    wanted = set(full_ids)
    found = {}
    with open(world.regional_csv, newline="") as handle:
        header = handle.readline().rstrip("\r\n").split(",")
        if tuple(header[1:]) != tuple(regions):
            raise RunRefused("the regional table does not hold the committed regions in order")
        for line in handle:
            cut = line.find(",")
            if line[:cut] not in wanted:
                continue
            values = np.array(line[cut + 1:].rstrip("\r\n").split(","), dtype=float)
            if values.shape != (len(regions),) or not np.isfinite(values).all() or (values < 0).any():
                raise RunRefused(f"invalid regional row {line[:cut]}")
            found[line[:cut]] = values
    return found


def verify(forecast_dir, commitment):
    forecast_dir = Path(forecast_dir)
    if nc.sha_file(forecast_dir / "manifest.json") != commitment:
        raise RunRefused("the forecast manifest does not hash to the recorded commitment")
    manifest = json.loads((forecast_dir / "manifest.json").read_text())
    for name, sha in manifest["forecast_sha256"].items():
        if nc.sha_file(forecast_dir / "forecasts" / f"{name}.npz") != sha:
            raise RunRefused(f"forecast file of {name} changed after the commitment")
    for f, sha in manifest["code_sha256"].items():
        if nc.sha_file(ROOT / f) != sha:
            raise RunRefused(f"{f} changed after the commitment")
    return manifest


def score(world, part, forecast_dir, commitment, selected, out_dir, after_engine_score=None, seed=20261111,
          resamples=10000):
    out_dir = Path(out_dir)
    manifest = verify(forecast_dir, commitment)
    if manifest["part"] != part:
        raise RunRefused("the committed forecasts are for another part")
    engine = None
    if not world.synthetic:
        if not after_engine_score or not Path(after_engine_score).is_file():
            raise RunRefused("the current evaluation's test score file must exist before this scoring")
        engine = {"path": str(after_engine_score), "sha256": nc.sha_file(after_engine_score)}
    if selected not in (None, "none") and selected not in manifest["candidates"]:
        raise RunRefused("the selected candidate has no committed forecast")
    out_dir.mkdir(parents=True, exist_ok=True)
    started = {"commitment": commitment, "selected": selected, "after_engine_score": engine}
    with open(out_dir / "scoring_started.json", "x") as handle:     # a started scoring is never repeated
        handle.write(json.dumps(started, indent=1, sort_keys=True) + "\n")
    if selected in (None, "none"):
        result = {"outcome": "not_run", "reason": "no candidate qualified on training; no regional row was opened"}
        (out_dir / "scores.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n")
        return result
    regions = manifest["regions"]
    rows = read_regional_rows(world, manifest["query_full_ids"], regions)
    have = np.array([f in rows for f in manifest["query_full_ids"]])
    missing = [f for f, h in zip(manifest["query_full_ids"], have) if not h]
    Y = np.log10(1.0 + np.vstack([rows[f] for f, h in zip(manifest["query_full_ids"], have) if h]))
    group = np.array(manifest["query_samples"])[have]
    animal, per_sample, regional, filled = {}, {}, {}, 0
    fc = {}
    for name in manifest["forecast_sha256"]:
        with np.load(Path(forecast_dir) / "forecasts" / f"{name}.npz") as z:
            fc[name] = (z["point"][have], z["q"][have].astype(np.float64), z["p_nonzero"][have])
    for name, (point, q, p) in fc.items():
        if name == selected:
            bad = ~(np.isfinite(point).all(axis=1) & np.isfinite(q).all(axis=(1, 2)) & np.isfinite(p).all(axis=1))
            if bad.any():                                    # a missing forecast counts as MEAN's (cannot help)
                filled = int(bad.sum())
                mp, mq, mpp = fc["MEAN"]
                point, q, p = (np.where(bad[:, None], mp, point), np.where(bad[:, None, None], mq, q),
                               np.where(bad[:, None], mpp, p))
        scores, mats = nc.score_forecast(point, q, p, Y)
        animal[name], per_sample[name] = {}, {}
        for k, v in scores.items():
            animal[name][k], _, per_sample[name][k] = nc.animal_mean(v, group)
        regional[name] = nc.region_checks(mats, group)
    animal_c = dict(animal)
    animal_c["CANDIDATE"] = animal[selected]
    decision = nc.decide(animal_c, regional[selected])
    paired = {}
    for ctrl in sorted(set(nc.POINT_CONTROLS) | set(nc.DIST_CONTROLS)):
        for metric in ("log_mse", "crps"):
            d = per_sample[selected][metric] - per_sample[ctrl][metric]
            paired[f"{metric}:{selected}-{ctrl}"] = {"mean": float(np.mean(d)), "samples_favouring_candidate": int((d < 0).sum()),
                                                     "bootstrap95": nc.bootstrap_mean_interval(d, resamples, seed)}
    result = {"schema": "next_regional_scores_v1", "part": part, "commitment": commitment, "selected": selected,
              "exposure_label": EXPOSURE_LABEL, "missing_regional_rows": missing, "candidate_forecasts_filled_by_MEAN": filled,
              "animal_mean": animal, "decision": decision, "paired_descriptive": paired,
              "per_region": {selected: regional[selected]["per_region"]}, "after_engine_score": engine}
    (out_dir / "scores.json").write_text(json.dumps(result, indent=1, sort_keys=True, default=float) + "\n")
    if decision["outcome"] == "pass":
        point, q, p = fc[selected]
        lo, hi, _, _ = nc.intervals(q, p)
        payload = viewer_input(manifest, have, point, lo, hi, decision, commitment, nc.sha_file(out_dir / "scores.json"))
        (out_dir / "viewer_regional_input.json").write_text(json.dumps(payload, sort_keys=True) + "\n")
    return result


def viewer_input(manifest, have, point, lo, hi, decision, protocol_sha, results_sha):
    """The input contract of analysis/prediction_viewer_export.py export_regional, with an evidence level the
    current viewer does not know: the export refuses until the viewer's owner adds that level."""
    neurons = np.array(manifest["query_neurons"])[have]
    samples = np.array(manifest["query_samples"])[have]
    return {"task": "next_regional_profile", "scale": "log10(1 + published regional axon length)",
            "regions": manifest["regions"],
            "status": {"beats_all_baselines": decision["outcome"] == "pass", "evidence_level": VIEWER_EVIDENCE_LEVEL,
                       "benchmark": EXPOSURE_LABEL, "protocol_sha256": protocol_sha, "results_sha256": results_sha,
                       "model_sha256": manifest["code_sha256"]["analysis/next_regional_core.py"]},
            "units": [{"unit_id": str(n), "group_id": str(s), "predicted": [float(v) for v in point[i]],
                       "interval_low": [float(v) for v in lo[i]], "interval_high": [float(v) for v in hi[i]],
                       "coverage_validated": bool(decision["leg3_all"])} for i, (n, s) in enumerate(zip(neurons, samples))]}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check-loader")
    c.add_argument("--out", default=str(OUT))
    s = sub.add_parser("select")
    s.add_argument("--out", required=True)
    p = sub.add_parser("predict")
    p.add_argument("--part", required=True, choices=["validation", "test"])
    p.add_argument("--selection", required=True)
    p.add_argument("--out", required=True)
    r = sub.add_parser("score")
    r.add_argument("--part", required=True, choices=["validation", "test"])
    r.add_argument("--forecasts", required=True)
    r.add_argument("--commitment", required=True)
    r.add_argument("--selected", required=True)
    r.add_argument("--after-engine-score", required=True)
    r.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    world = general_core().real_world(ROOT)
    if args.cmd == "check-loader":
        print(json.dumps(check_loader(world, args.out)["audit"], indent=1, sort_keys=True))
    elif args.cmd == "select":
        select(world, args.out)
    elif args.cmd == "predict":
        print(predict(world, args.part, json.loads(Path(args.selection).read_text()), args.out))
    elif args.cmd == "score":
        print(json.dumps(score(world, args.part, args.forecasts, args.commitment, args.selected, args.out,
                               args.after_engine_score)["decision"], indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
