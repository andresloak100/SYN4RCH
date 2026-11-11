"""Tests of analysis/next_regional_run.py on synthetic worlds: query loader, hidden-outcome canary, scorer.

No real cohort file is read here."""
import importlib.util
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]


def load(name, rel):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


run = load("next_regional_run_t", "analysis/next_regional_run.py")
nc = run.nc
REGIONS = ("ARH", "VMH", "LHA", "ZI", "PAG", "MY")
POISON = "SEALED_CANARY"


def make_world(directory, seed=0):
    """Manifest, soma.csv, soma.xlsx, regional.csv in the real formats. Validation rows carry poison in every
    soma value (a loader of the test part must never parse them)."""
    import openpyxl
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    manifest, soma_rows, xlsx_rows, regional_rows = [], [], [], []
    counts = {"train_samples": 0, "train_neurons": 0}
    number = 400000
    plan = [("Pomc", "train", 4), ("Sst", "train", 4), ("Pomc", "test", 2), ("Sst", "test", 1), ("Sst", "validation", 1)]
    for label, part, n_samples in plan:
        for _ in range(n_samples):
            number += 1
            sample = str(number)
            for i in range(1, 16 if part == "train" else 8):
                neuron = f"{sample}_{i:03d}"
                fine = "ARH" if label == "Pomc" else "VMH"
                soma = np.array([6000.0, 6000.0, 4800.0]) + rng.normal(0, 150, 3) + (300.0 if label == "Sst" else 0.0)
                on = rng.random(len(REGIONS)) < (0.6 if label == "Pomc" else 0.3)
                length = np.where(on, 10 ** rng.uniform(1.5, 3.5, len(REGIONS)), 0.0)
                full = f"{neuron}_{label}_{fine}"
                manifest.append([neuron, sample, label, fine, "synthetic", "0", "0", "0", part, "synthetic"])
                vals = [f"{v:.6f}" for v in soma] if part != "validation" else [POISON] * 3
                soma_rows.append([full, "FORBIDDEN_TYPE", neuron, "FORBIDDEN_NAME", label, fine] + vals)
                xlsx_rows.append([neuron, fine])
                regional_rows.append([full] + ([f"{v:.4f}" for v in length] if part != "validation" else [POISON] * len(REGIONS)))
            if part == "train":
                counts["train_samples"] += 1
                counts["train_neurons"] += 15
    with open(directory / "manifest.tsv", "w") as h:
        h.write("\t".join(["neuron", "animal", "source_label", "soma_label", "class", "animal_inspected",
                           "animal_development", "animal_in_li_2024", "part", "part_reason"]) + "\n")
        for row in sorted(manifest):
            h.write("\t".join(row) + "\n")
    with open(directory / "soma.csv", "w") as h:
        h.write("file,type,neurons,name,geno,location,AP,DV,ML\n")
        for row in soma_rows:
            h.write(",".join(row) + "\n")
    with open(directory / "regional.csv", "w") as h:
        h.write("," + ",".join(REGIONS) + "\n")
        for row in regional_rows:
            h.write(",".join(row) + "\n")
    wb = openpyxl.Workbook(write_only=True)
    sheet = wb.create_sheet("Sheet1")
    sheet.append(["neurons", "location"])
    for row in xlsx_rows:
        sheet.append(row)
    wb.save(directory / "soma.xlsx")
    meta = {"name": f"next_regional_synthetic_{seed}", "synthetic": True, "regions": list(REGIONS),
            "t4_observed": list(REGIONS[:2]), "expected": counts, "label_alias": []}
    (directory / "world.json").write_text(json.dumps(meta))
    return run.general_core().synthetic_world(directory)


def rewrite_test_rows(world, how):
    """Change the regional rows of the TEST samples only: 'perturb', 'delete' or 'poison'."""
    test_neurons = {l.split("\t")[0] for l in Path(world.manifest).read_text().splitlines()[1:] if l.split("\t")[8] == "test"}
    lines = Path(world.regional_csv).read_text().splitlines()
    out = [lines[0]]
    for line in lines[1:]:
        full = line.split(",", 1)[0]
        if "_".join(full.split("_")[:2]) in test_neurons:
            if how == "delete":
                continue
            values = line.split(",")[1:]
            values = [f"{float(v) * 3 + 7:.4f}" for v in values] if how == "perturb" else [POISON] * len(values)
            line = ",".join([full] + values)
        out.append(line)
    Path(world.regional_csv).write_text("\n".join(out) + "\n")


@pytest.fixture(scope="module")
def base(tmp_path_factory):
    d = tmp_path_factory.mktemp("world")
    world = make_world(d / "w")
    selection = run.select(world, d / "sel")
    return d, world, selection


def test_query_loader_reads_one_part_and_never_the_regional_table(base, tmp_path):
    d, world, _ = base
    copy = tmp_path / "w"
    shutil.copytree(Path(world.manifest).parent, copy)
    (copy / "regional.csv").unlink()                      # no regional table at all
    w2 = run.general_core().synthetic_world(copy)
    query, audit = run.load_query_inputs(w2, "test")      # validation rows carry poison: never parsed
    assert audit["neurons"] == 3 * 7 and audit["samples"] == 3
    assert set(query["label"]) == {"Pomc", "Sst"} and query["soma"].shape == (21, 3)
    with pytest.raises(ValueError):
        run.load_query_inputs(w2, "validation")            # parsing poison fails: that part was never parsed above


def test_check_loader_agrees_with_the_sealed_loader(base, tmp_path):
    _, world, _ = base
    assert run.check_loader(world, tmp_path)["agrees_with_sealed_loader"]


@pytest.mark.parametrize("how", ["perturb", "delete", "poison"])
def test_hidden_outcome_canary(base, tmp_path, how):
    """Forecasts are bit-identical whatever the test samples' regional rows hold, or if they are absent."""
    d, world, selection = base
    ref_dir = d / "ref"
    if not (ref_dir / "manifest.json").exists():
        run.predict(world, "test", selection, ref_dir)
    ref = json.loads((ref_dir / "manifest.json").read_text())["forecast_sha256"]
    copy = tmp_path / "w"
    shutil.copytree(Path(world.manifest).parent, copy)
    w2 = run.general_core().synthetic_world(copy)
    rewrite_test_rows(w2, how)
    run.predict(w2, "test", selection, tmp_path / "f")
    got = json.loads((tmp_path / "f" / "manifest.json").read_text())["forecast_sha256"]
    assert got == ref and set(got) == set(run.CONTROL_METHODS)


def test_predict_refuses_training_part_and_rewrites(base, tmp_path):
    d, world, selection = base
    with pytest.raises(run.RunRefused):
        run.predict(world, "train", selection, tmp_path / "x")
    run.predict(world, "test", selection, tmp_path / "y")
    with pytest.raises(run.RunRefused):
        run.predict(world, "test", selection, tmp_path / "y")


def test_score_checks_commitment_runs_once_and_decides(base, tmp_path, monkeypatch):
    d, world, selection = base
    monkeypatch.setitem(run.CANDIDATES, "CAND",
                        lambda ctrl, q: ctrl.forecast("KNN", {"k": 10, "pool": "all", "scaling": "raw"}, q["soma"], q["label"]))
    commitment = run.predict(world, "test", selection, tmp_path / "f")
    with pytest.raises(run.RunRefused):
        run.score(world, "test", tmp_path / "f", "0" * 64, "CAND", tmp_path / "s0")
    with pytest.raises(run.RunRefused):
        run.score(world, "test", tmp_path / "f", commitment, "UNKNOWN", tmp_path / "s0")
    res = run.score(world, "test", tmp_path / "f", commitment, "CAND", tmp_path / "s1")
    assert res["decision"]["outcome"] in ("pass", "better_crps_without_valid_calibration", "better_crps_but_worse_points",
                                          "no_gain_over_the_controls")
    assert res["missing_regional_rows"] == [] and set(res["animal_mean"]) == set(run.CONTROL_METHODS) | {"CAND"}
    with pytest.raises(FileExistsError):
        run.score(world, "test", tmp_path / "f", commitment, "CAND", tmp_path / "s1")     # never repeated
    # tampering after the commitment is refused
    path = tmp_path / "f" / "forecasts" / "CAND.npz"
    path.write_bytes(path.read_bytes() + b"x")
    with pytest.raises(run.RunRefused):
        run.score(world, "test", tmp_path / "f", commitment, "CAND", tmp_path / "s2")


def test_score_not_run_never_opens_regional_rows(base, tmp_path, monkeypatch):
    d, world, selection = base
    commitment = run.predict(world, "test", selection, tmp_path / "f")

    def boom(*a, **k):
        raise AssertionError("regional rows opened")
    monkeypatch.setattr(run, "read_regional_rows", boom)
    assert run.score(world, "test", tmp_path / "f", commitment, "none", tmp_path / "s")["outcome"] == "not_run"


def test_viewer_export_refuses_the_new_evidence_level(base, tmp_path):
    view = load("prediction_viewer_export_t", "analysis/prediction_viewer_export.py")
    payload = {"task": "t", "scale": "s", "regions": list(REGIONS),
               "status": {"beats_all_baselines": True, "evidence_level": run.VIEWER_EVIDENCE_LEVEL, "benchmark": run.EXPOSURE_LABEL,
                          "protocol_sha256": "a", "results_sha256": "b", "model_sha256": "c"}, "units": []}
    (tmp_path / "p.json").write_text(json.dumps(payload))
    with pytest.raises(view.ExportRefused):
        view.export_regional(tmp_path / "p.json", tmp_path / "out.json")
    assert not (tmp_path / "out.json").exists()


def test_run_source_never_names_traces_or_engine_outputs():
    text = (ROOT / "analysis/next_regional_run.py").read_text()
    for forbidden in ("traces/", "brm_real/truth", "paircache", "axon_engine.py"):
        assert forbidden not in text
