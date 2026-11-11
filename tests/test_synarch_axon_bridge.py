"""Checks of the bridge between the general reconstruction layer and the axon engine
(analysis/synarch_axon_bridge.py) and of the engine candidate that uses it (E7_local_kernel_inferred_profile).

Everything statistical here runs on SYNTHETIC worlds made in temporary directories. Nothing here is biological
evidence. The few checks that name the real cohort only show a refusal, or the recorded absence of a finished run:
they open no trace and no regional row, and never touch data/raw/jiao2025/traces.

The three guarantees of the general layer are checked here IN THE CONTEXT OF THE AXON WORKFLOW:
  K1  what may feed a forecast of a hidden axon: only permitted measured inputs, and inferred values that the layer
      itself issued for this neuron from states that never saw its brain sample;
  K2  every value states its scale, tag and frame and is converted exactly once;
  K3  a saved state is what the finished run recorded, and tampering is refused.

    python3 -I -m pytest tests/test_synarch_axon_bridge.py -q
"""
import csv
import dataclasses
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

# every process that calls the referee runs with exactly these four thread variables (late item L6), set before numpy
# is imported, because numpy's thread pools are fixed when it loads (round-2 finding: the engine checks the pools)
os.environ.update({v: "1" for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")})

import numpy as np  # noqa: E402
import pytest  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import axon_engine as E  # noqa: E402

B = E.bridge()                      # the bridge, loaded the way the engine loads it
core, mods, bench = B.core, B.mods, B.bench
TARGET, E7, E2 = B.TARGET, "E7_local_kernel_inferred_profile", "E2_local_kernel"
BLOCKED = "blocked: independent validation pending"
REAL = (ROOT / "reports/axon_dataset_manifest.tsv").exists() and (ROOT / "data/raw/jiao2025/swc/MANIFEST.tsv").exists()
real_only = pytest.mark.skipif(not REAL, reason="the real cohort tables are not on disk")


# ------------------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------------------


def load(path):
    return json.loads(Path(path).read_text())


def clone(source, target):
    shutil.copytree(source, target)
    return target


def engine_of(directory, stage, **kw):
    return E.Engine(E.synthetic_world(directory), stage, workers=1, log=None, **kw)


def session_of(directory):
    world = E.synthetic_world(directory)
    cohort = E.Cohort(world)
    return B.open_session(world, heldout_groups=cohort.heldout_animals), cohort


def layer_out(directory):
    return Path(directory) / "engine" / B.LAYER_FOLDER / "out"


def forecasts_of(engine, neuron):
    forecasts, log = engine.forecast(neuron)
    return {f["model"]: f for f in forecasts}, log


def array_hashes(engine, neurons):
    out = {}
    for neuron in neurons:
        for model, f in forecasts_of(engine, neuron)[0].items():
            out[(neuron, model)] = E.arbor_sha(f["arbor"])
    return out


def profile_record(session, row, names, values, *, group=None):
    """A MEASURED regional profile of a neuron as a record (what must never reach a forecast of its hidden axon)."""
    registry = session.registry
    spec, namespace = registry.variable("jiao", TARGET), registry.namespace("jiao")
    return core.Observation(session.world.species, "jiao", namespace.grain, str(row.neuron),
                            session.group_of(row) if group is None else group, TARGET, "measured",
                            tuple(float(v) for v in values), core.scale_statement(spec, "identity"), tuple(names), spec.frame,
                            source=core.SourceRef("test", None, str(row.neuron)), transform="identity")


def regional_rows(directory):
    """{full ID: published lengths} of the TRAINING rows of a synthetic world's regional table, and the region names."""
    world = E.synthetic_world(directory)
    cohort = E.Cohort(world)
    train = {cohort.rows[n].full_id for n in cohort.part_neurons["train"]}
    lines = Path(world.regional_table).read_text().splitlines()
    names = lines[0].split(",")[1:]
    rows = {}
    for line in lines[1:]:
        full_id, rest = line.split(",", 1)
        if full_id in train:
            rows[full_id] = np.array(rest.split(","), dtype=float)
    return rows, names


def bound_serving(directory):
    """The record of the all-training states that the completion record of the layer's finished run binds."""
    return load(layer_out(directory) / "run_complete.json")[bench.SERVING_KEY]


def read_state(path):
    with np.load(path, allow_pickle=False) as data:
        return json.loads(str(data["head"])), {k: data[k] for k in data.files if k != "head"}


def write_state(path, head, arrays, *, rehash=False):
    """Write a state file back. rehash=True is what an editor who knows the format would do: the hash stored in the
    same file is recomputed, so the file agrees with itself again."""
    if rehash:
        body = {k: v for k, v in head.items() if k != "state_hash"}
        head = dict(body, state_hash=mods.JointState._fingerprint(body, [arrays[f"a{i}"] for i in range(len(body["index"]))]))
    np.savez(path, head=np.array(core.canonical(head)), **arrays)


def scramble_regional_rows(directory, animal, seed=3):
    """Other numbers in the regional rows of ONE animal; nothing else is touched."""
    path = Path(directory) / "data/regional.csv"
    rng = np.random.default_rng(seed)
    lines = path.read_text().splitlines()
    changed = 0
    for i in range(1, len(lines)):
        full_id = lines[i].split(",", 1)[0]
        if full_id.split("_", 1)[0] == animal:
            lines[i] = full_id + "," + ",".join(f"{v:.2f}" for v in rng.uniform(0.0, 5000.0, len(E.REGIONS)))
            changed += 1
    path.write_text("\n".join(lines) + "\n")
    return changed


# ------------------------------------------------------------------------------------------
# fixtures
# ------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def plain(tmp_path_factory):
    """A tiny synthetic world with the engine's library, all-pairs cache and regional states, and NO run of the
    general layer."""
    directory = tmp_path_factory.mktemp("plain") / "informative"
    E.make_synthetic_world(directory, "informative", size="tiny")
    engine = engine_of(directory, "develop")
    engine.build_library()
    engine.build_cache()
    engine.states.build_many(list(engine.states.plan))
    return directory


@pytest.fixture(scope="module")
def layered(plain, tmp_path_factory):
    """The same world with the general layer's own benchmark run on its training rows (through the bridge)."""
    directory = clone(plain, tmp_path_factory.mktemp("layered") / "informative")
    summary = engine_of(directory, "develop").prepare_general_layer()
    assert summary["finished_run"] and summary["verified"] is True
    return directory


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    """The general layer's OWN small planted world (soma, label and profile share latent classes), written in the
    engine's table layout, with the layer's benchmark run on it through the bridge. In that world the layer's
    benchmark supports its engine for the profile forecast (reports/synarch_recon_protocol.md, section 17, first
    row), so the bridge receives an ENGINE-level value here. The world has no traces: the bridge needs none.
    Held-out neurons are made up (permitted metadata only; their regional rows hold text that cannot be parsed)."""
    base = tmp_path_factory.mktemp("planted")
    source = bench.make_synthetic_world(base / "layer_source", "planted", size="tiny", seed=0)
    table = core.load_training_table(source)                         # the layer's sealed loader: training rows only
    manifest_label = {"Wt(orexin)": "Orexin"}
    data = base / "world" / "data"
    data.mkdir(parents=True)
    head = ["neuron", "animal", "source_label", "soma_label", "class", "animal_inspected", "animal_development",
            "animal_in_li_2024", "part", "part_reason", "hemisphere", "reconstruction_type", "index_soma_acronym"]
    manifest, soma, regional = [], [], []
    for i, neuron in enumerate(table.neuron):
        neuron, sample, label, fine = str(neuron), str(table.sample[i]), str(table.label[i]), str(table.fine_label[i])
        manifest.append([neuron, sample, manifest_label.get(label, label), fine, "synthetic", "0", "0", "0", "train", "planted",
                         "Right" if i % 2 else "Left", "Axon_and_dendrite", fine])
        soma.append([str(table.full_id[i]), neuron, label, fine] + [f"{v:.6f}" for v in table.soma[i]])
        regional.append([str(table.full_id[i])] + [f"{v:.4f}" for v in table.length[i]])
    rng = np.random.default_rng(5)
    labels = sorted(set(table.label.tolist()))
    for part, first in (("validation", 900001), ("test", 900101)):
        for a in range(3):
            sample, label = str(first + a), labels[a % len(labels)]
            for k in range(1, 4):
                neuron, fine = f"{sample}_{k:03d}", "ARH"
                xyz = table.soma[int(rng.integers(len(table)))] + rng.normal(0.0, 30.0, 3)
                xyz[2] = min(xyz[2], 5690.0)
                full = f"{neuron}_{label}_{fine}"
                manifest.append([neuron, sample, manifest_label.get(label, label), fine, "synthetic", "0", "0", "0", part,
                                 "made up", "Left", "Axon_only", fine])
                soma.append([full, neuron, label, fine] + [f"{v:.6f}" for v in xyz])
                regional.append([full] + ["SEALED_ROW_NEVER_PARSED"] * len(table.regions))
    (data / "manifest.tsv").write_text("\n".join("\t".join(r) for r in [head] + sorted(manifest)) + "\n")
    (data / "soma.csv").write_text("\n".join(",".join(r) for r in [["file", "neurons", "geno", "location", "AP", "DV", "ML"]] + soma) + "\n")
    (data / "regional.csv").write_text("\n".join(",".join(r) for r in [[""] + list(table.regions)] + regional) + "\n")
    E.write_json(data / "world.json", {"name": "synthetic_layer_planted", "kind": "planted"}, pretty=True)
    world = E.synthetic_world(base / "world")
    status = B.ensure_synthetic_run(world)
    assert status["verified"] is True and status["status"] == "exploratory_supported", status
    return base / "world"


@pytest.fixture(scope="module")
def e7_world(tmp_path_factory):
    """The whole pipeline, once, on the synthetic world in which the validation gate selects E7: 27 soma zones of
    200 um set the dominant route, more zones than the frozen regional model has components. Standard size. The
    generator seed of this world was chosen so that E7 is selected (analysis/axon_engine.py, SYNTHETIC_SEEDS)."""
    directory = tmp_path_factory.mktemp("e7_world") / "fine_zones"
    summary = E.run_synthetic("fine_zones", directory, size="standard", workers=1, log=None, family="all", develop=False,
                              canary=False)
    return directory, summary


@pytest.fixture(scope="module")
def passing_world(tmp_path_factory):
    """A synthetic world whose test PASSES (the soma zone sets the route; family without a fitted state; the kernel
    candidate E2 is selected and beats every control), for the viewer export gate. Standard size."""
    directory = tmp_path_factory.mktemp("passing_world") / "dense"
    summary = E.run_synthetic("dense", directory, size="standard", workers=1, log=None, family="no_fitted_state", develop=False,
                              canary=False)
    return directory, summary


# ------------------------------------------------------------------------------------------
# the converter and the world mapping
# ------------------------------------------------------------------------------------------


def test_the_converter_copies_training_rows_only_and_no_held_out_outcome(plain, tmp_path):
    directory = clone(plain, tmp_path / "a")
    world = E.synthetic_world(directory)
    cohort = E.Cohort(world)
    held = [n for part in ("validation", "test") for n in cohort.part_neurons[part]]
    train = cohort.part_neurons["train"]
    layer_world = B.convert_synthetic(world)
    folder = B.layer_directory(world)
    record = load(folder / "conversion.json")
    assert record["held_out_rows_copied"] == 0 and record["training_neurons"] == len(train)
    text = {name: (folder / name).read_text() for name in ("manifest.tsv", "soma.csv", "regional.csv")}
    assert not any(n in body for n in held for body in text.values())           # no held-out neuron, in any table
    assert not any(cohort.rows[n].animal in body for n in held for body in text.values())
    source = {line.split(",", 1)[0]: line for line in Path(world.regional_table).read_text().splitlines()[1:]}
    copied = text["regional.csv"].splitlines()
    assert copied[1:] == [source[cohort.rows[n].full_id] for n in sorted(train, key=lambda n: list(source).index(cohort.rows[n].full_id))]
    assert len(copied) - 1 == len(train) and copied[0] == "," + ",".join(E.REGIONS)
    import openpyxl
    sheet = openpyxl.load_workbook(folder / "soma.xlsx", read_only=True).worksheets[0]
    cells = list(sheet.iter_rows(values_only=True))
    assert cells[0] == ("neurons", "location") and {c[0] for c in cells[1:]} == set(train)
    # the layer's own sealed loader reads exactly the training partition from it
    table = core.load_training_table(layer_world)
    assert sorted(table.neuron.tolist()) == sorted(train) and table.length.shape == (len(train), len(E.REGIONS))
    assert load(folder / "world.json")["config"] == bench.TINY_CONFIG       # the layer's own small declared grid, unchanged
    # held-out outcomes replaced, sealed rows made unparseable, held-out traces deleted: the conversion is the same, byte for byte
    second = clone(plain, tmp_path / "b")
    E.scramble_outcomes(second)
    lines = (second / "data/regional.csv").read_text().splitlines()
    held_full = {cohort.rows[n].full_id for n in held}
    (second / "data/regional.csv").write_text("\n".join(
        l if l.split(",", 1)[0] not in held_full else l.split(",", 1)[0] + ",THIS,ROW,IS,SEALED" for l in lines) + "\n")
    for part in ("validation", "test"):
        shutil.rmtree(second / "data/traces" / part)
    B.convert_synthetic(E.synthetic_world(second))
    other = B.layer_directory(E.synthetic_world(second))
    assert all((other / name).read_bytes() == (folder / name).read_bytes() for name in text)
    assert load(other / "conversion.json")["training_content_sha256"] == record["training_content_sha256"]
    with pytest.raises(B.BridgeRefused, match="never overwritten"):
        B.convert_synthetic(world)
    # a training row that changes after the conversion makes the layer world stale: refused, not used
    assert B.verify_conversion(world)["training_neurons"] == len(train)
    scramble_regional_rows(directory, cohort.rows[train[0]].animal)
    with pytest.raises(B.BridgeRefused, match="not the conversion of this engine world's training rows"):
        B.verify_conversion(world)


def test_on_the_real_cohort_the_bridge_never_starts_the_layers_benchmark(tmp_path):
    """The engine's real world maps to the layer's real world; the bridge reads a finished run or says there is
    none. It refuses to convert, freeze, run or serve there. Nothing of the cohort is opened by this test."""
    real = E.real_world(out=tmp_path, legacy_only=True)
    with pytest.raises(B.BridgeRefused, match="never freezes, runs, scores or serves"):
        B.ensure_synthetic_run(real)
    with pytest.raises(B.BridgeRefused, match="only an engine SYNTHETIC world is converted"):
        B.convert_synthetic(real)
    layer = B.layer_world_for(real)
    assert layer == core.real_world() and not layer.synthetic and Path(layer.out) == ROOT / "data/interim/synarch_recon"
    wrong = dataclasses.replace(real, regional_table=str(tmp_path / "another.csv"))
    with pytest.raises(B.BridgeRefused, match="do not describe the same cohort"):
        B.layer_world_for(wrong)
    assert not (tmp_path / B.LAYER_FOLDER).exists()


@real_only
def test_real_cohort_as_it_is_now_the_finished_run_is_not_named_in_the_protocol_so_held_out_forecasts_are_refused(tmp_path):
    """Review finding 8: the state of the real cohort on disk now. The general layer's owner has finished its single run
    (data/interim/synarch_recon/run_complete.json), and the protocol does not name it on its anchor line: a held-out
    forecast is refused for both candidate families, and the engine cannot start, replace or serve the layer's
    benchmark. Only the 37 local training traces are listed (legacy_only); every engine output goes to tmp_path; no
    trace, no regional row and nothing under data/raw/jiao2025/traces is opened."""
    real = E.real_world(out=tmp_path, legacy_only=True)
    run = ROOT / "data/interim/synarch_recon/run_complete.json"
    if not run.exists():
        pytest.skip("no finished run of the general layer on disk")
    sha = E.sha_file(run)
    for family in ("all", "no_fitted_state"):
        engine = E.Engine(real, "validate", workers=1, log=None, family=family)
        with pytest.raises(E.EngineError, match="never freezes, runs, scores or serves"):
            engine.prepare_general_layer()
        with pytest.raises(E.EngineError, match=f"{sha}.*is not named on the protocol's anchor line"):
            engine.require_general_layer()
        with pytest.raises(E.StageBlocked, match="have not been fetched"):
            engine.predict()
    status = B.run_status(real)
    assert status["finished_run"] is True and status["run_complete_sha256"] == sha
    assert E.Engine(real, "validate", workers=1, log=None).protocol_anchor() is None      # '(not written yet)'
    assert not (tmp_path / "validate").exists() and not (tmp_path / "states").exists()  # nothing written, nothing registered


# ------------------------------------------------------------------------------------------
# K1 in context: what may feed a forecast of the hidden axon
# ------------------------------------------------------------------------------------------


def test_K1a_the_querys_measured_regional_profile_is_never_evidence(planted):
    session, cohort = session_of(planted)
    regions, inside = session.world.regions, session.world.t4_observed
    for part in ("validation", "train"):                                 # also not for a training neuron in development
        row = cohort.rows[cohort.part_neurons[part][0]]
        fold = session.fold_for(row)
        measured = session.measured_records(row)
        assert set(measured) == set(core.HIDDEN_AXON_INPUTS) and all(r.status == "measured" for r in measured.values())
        t1 = {v: measured[v] for v in B.T1_INPUTS}
        whole = profile_record(session, row, regions, [120.0] * len(regions))
        partial = profile_record(session, row, inside, [120.0] * len(inside))
        for offered in (whole, partial):
            request, refusal = session.ask(row, fold, records={**t1, TARGET: offered})
            assert request is None and refusal.code == "forbidden_input"
            assert "never an input to a hidden-axon forecast" in refusal.reason
            query, why = session.admit(row, fold, {**measured, TARGET: offered}, None)
            assert query is None and "never an input to a hidden-axon forecast" in why
        # why the bridge refuses BEFORE it asks: with a partly measured profile the layer would answer another,
        # declared task (completion of a measured profile), which is not a forecast of a hidden axon
        reconstructor = session.layer(fold)[0]
        direct = reconstructor.reconstruct(session._request(row, [*t1.values(), partial], TARGET))[TARGET]
        assert direct.status == "inferred" and direct.conditional.startswith("T4_")
        answer = session.infer(row)
        assert answer.profile is not None and TARGET not in answer.query.inputs
        assert dict(answer.record.derived_from) == {"soma_position": "measured", "source_label": "measured"}


def test_K1b_a_training_query_is_answered_by_the_fold_state_that_never_saw_its_brain_sample(planted):
    session, cohort = session_of(planted)
    row = cohort.rows[cohort.part_neurons["train"][0]]
    k = session.fold_for(row)
    j = next(f for f in sorted(set(session.fold_of_group.values())) if f != k)
    state_k, state_j, serving = session.layer(k)[2], session.layer(j)[2], session.layer(None)[2]
    assert row.animal not in state_k.groups and row.animal in state_j.groups and row.animal in serving.groups
    answer = session.infer(row)                                          # the bridge picks fold k by itself
    assert answer.profile is not None and answer.provenance["fold"] == f"outer{k}"
    assert answer.provenance["state_sha256"] == state_k.state_hash() and row.animal not in answer.upstream_groups
    for wrong in (j, None):                                              # another fold's state, or the all-training state
        with pytest.raises(B.BridgeLeak, match="its brain sample .* is in the T1 state"):
            session.infer(row, fold=wrong)
        request, result = session.ask(row, wrong)                        # the layer itself: the neuron is in that state
        assert result.status == "refused" and result.refusal.code == "leakage_guard"
    # a NEW neuron of the same brain sample is answered by state j with a warning; as an input of the axon forecast
    # that value is refused for the state that produced it
    sibling = dataclasses.replace(row, neuron=f"{row.animal}_99999", full_id=f"{row.animal}_99999_x_y")
    request, result = session.ask(sibling, j)
    assert result.status == "inferred" and "group_in_training_state" in result.codes()
    record = result.to_observation(request)
    query, why = session.admit(sibling, j, session.measured_records(sibling), record)
    assert query is None and "contains the query's brain sample" in why
    request, result = session.ask(sibling, k)                            # the same neuron by the state that excludes the sample
    query, why = session.admit(sibling, k, session.measured_records(sibling), result.to_observation(request))
    assert query is not None and why is None
    # a neuron outside the training partition is answered by the all-training states, and only by them
    held = cohort.rows[cohort.part_neurons["test"][0]]
    assert session.fold_for(held) is None and session.infer(held).provenance["fold"] == "serving"
    moved = dataclasses.replace(held, part="train")                      # the two partitions must agree about a brain sample
    with pytest.raises(B.BridgeLeak, match="does not hold it as a training sample"):
        session.fold_for(moved)


def test_K1c_a_value_that_descends_from_the_measured_profile_is_refused_directly_and_through_a_chain(planted):
    session, cohort = session_of(planted)
    row = cohort.rows[cohort.part_neurons["validation"][0]]
    reconstructor = session.layer(None)[0]
    measured = session.measured_records(row)
    regions, inside = session.world.regions, session.world.t4_observed
    without = lambda name: {k: v for k, v in measured.items() if k != name}               # noqa: E731
    # directly: the soma position inferred from the MEASURED profile (the layer's declared task T2)
    whole = profile_record(session, row, regions, [150.0] * len(regions))
    request = session._request(row, [whole, measured["source_label"]], "soma_position")
    t2 = reconstructor.reconstruct(request)["soma_position"]
    assert t2.status == "inferred" and not t2.baseline_level
    soma_record = t2.to_observation(request)
    query, why = session.admit(row, None, without("soma_position"), soma_record)
    assert query is None and "descends from a measured regional_axon_length" in why
    # directly: the rest of the profile inferred from a measured part of it (T4)
    partial = profile_record(session, row, inside, [150.0] * len(inside))
    request = session._request(row, [measured["soma_position"], measured["source_label"], partial], TARGET)
    t4 = reconstructor.reconstruct(request)[TARGET]
    assert t4.status == "inferred" and t4.conditional.startswith("T4_")
    query, why = session.admit(row, None, measured, t4.to_observation(request))
    assert query is None and "descends from a measured regional_axon_length" in why
    # the finished run's own reconstructor issues no second step from such a value: whatever is asked next from the
    # tainted soma is refused (the label is among its ancestors; a profile from it would be a hidden-axon forecast)
    for target in ("source_label", TARGET):
        others = [soma_record] if target == "source_label" else [soma_record, measured["source_label"]]
        refused = reconstructor.reconstruct(session._request(row, others, target))[target]
        assert refused.status == "refused" and refused.refusal.code == "forbidden_input", target
    # through a chain: measured profile -> inferred soma -> inferred source label. The label's own step names only an
    # inferred soma. To have such a chain ISSUED at all, a second reconstruction layer is built over the SAME states
    # of the finished run with two test-only declared steps (soma from the profile; label from a soma), sharing the
    # ledger of issued values, as the general layer's own chain test does.
    _, module, state, _ = session.layer(None)
    registry, layout = session.registry, module.layout
    steps = dict(module.tasks)
    steps["T2b"] = mods.make_task(layout, registry, task="T2b", title="soma_from_profile", target="soma_position",
                                  evidence=((TARGET, None),), declared_measured_axon_input=True)
    steps["T3b"] = mods.make_task(layout, registry, task="T3b", title="label_from_soma_declared", target="source_label",
                                  evidence=(("soma_position", None),), declared_measured_axon_input=True)
    pointer = {"metrics": "test-only evidence", "metrics_sha256": "0" * 64}
    evidence = {t.conditional().name: core.ValidationEvidence(status="exploratory_supported", task=n, benchmark=pointer,
                                                              primary_metric=t.primary) for n, t in steps.items()}
    chained = core.Reconstructor(registry, None, session.ledger)
    chained.register(mods.JointMixtureModule(registry, layout, states={n: module.states.get(n, state) for n in steps},
                                             train=module.train, tasks=steps, evidence=evidence))
    chained.register(mods.AxonArborAdapter(registry))
    request = session._request(row, [whole], "soma_position")
    hop1 = chained.reconstruct(request)["soma_position"]
    assert hop1.status == "inferred", hop1.refusal
    from_profile = hop1.to_observation(request)
    request = session._request(row, [from_profile], "source_label")
    hop2 = chained.reconstruct(request)["source_label"]
    assert hop2.status == "inferred", hop2.refusal
    label_record = hop2.to_observation(request)
    assert dict(label_record.derived_from) == {"soma_position": "inferred"}                 # no profile in its own step
    assert [(s.variable, dict(s.evidence)) for s in label_record.lineage] == [
        ("source_label", {"soma_position": "inferred"}), ("soma_position", {TARGET: "measured"})]
    query, why = session.admit(row, None, without("source_label"), label_record)
    assert query is None and "descends from a measured regional_axon_length" in why
    # and the profile forecast itself is never asked with such an input
    for records in ({"soma_position": soma_record, "source_label": measured["source_label"]},
                    {"soma_position": measured["soma_position"], "source_label": label_record}):
        request, refusal = session.ask(row, None, records=records)
        assert request is None and refusal.code == "forbidden_input"
        assert "descends from a measured regional_axon_length" in refusal.reason


def test_K1d_a_hand_made_inferred_record_and_a_value_issued_for_another_neuron_are_refused(planted):
    session, cohort = session_of(planted)
    first, second = (cohort.rows[n] for n in cohort.part_neurons["validation"][:2])
    measured = session.measured_records(first)
    answer = session.infer(first)
    issued = answer.record
    state = session.layer(None)[2]
    spec = session.registry.variable("jiao", TARGET)
    key = ("jiao", "neuron", first.neuron)
    evidence = (("soma_position", "measured"), ("source_label", "measured"))
    uncertainty = core.Uncertainty("inferred", None, None, None, None, B.MODULE, state.state_hash(), state.groups_sha256, "unvalidated")
    handmade = core.Observation(session.world.species, "jiao", "neuron", first.neuron, session.group_of(first), TARGET, "inferred",
                                tuple(1.0 for _ in issued.coordinates), core.scale_statement(spec, "log10_1p"),
                                tuple(issued.coordinates), spec.frame, uncertainty,
                                core.SourceRef("reconstruct:" + B.MODULE, state.state_hash(), first.neuron),
                                derived_from=evidence, transform="log10_1p",
                                lineage=(core.LineageStep(TARGET, B.MODULE, state.state_hash(), state.groups_sha256, evidence, key,
                                                          first.animal),))
    assert handmade.lineage == issued.lineage and handmade.derived_from == issued.derived_from      # it looks exactly right
    query, why = session.admit(first, None, measured, handmade)
    assert query is None and "not issued by this reconstruction layer" in why
    changed = dataclasses.replace(issued, value=tuple(v + 0.25 for v in issued.value))             # issued, then edited
    query, why = session.admit(first, None, measured, changed)
    assert query is None and "changed after it was issued" in why
    query, why = session.admit(first, None, measured, issued)                                      # exactly as issued: admitted
    assert query is not None and why is None
    # a value the layer issued for ANOTHER neuron
    request_2, result_2 = session.ask(second, None)
    assert result_2.status == "inferred"
    with pytest.raises(core.PairingError, match="it cannot become a record"):
        result_2.to_observation(session._request(first, measured.values(), TARGET))
    theirs = result_2.to_observation(request_2)
    query, why = session.admit(first, None, measured, theirs)
    assert query is None and "different units" in why
    query, why = session.admit(first, None, {}, theirs)
    assert query is None and "not to the query unit" in why


def test_K1e_a_legitimate_chain_passes_and_the_adapter_stays_blocked(planted, layered):
    """Engine-level value (the planted world) and flagged baseline-level value (the tiny engine world): both are the
    layer's own issued records, both pass the adapter, and the adapter is still blocked afterwards."""
    for directory, level in ((planted, "engine"), (layered, "baseline")):
        session, cohort = session_of(directory)
        for part in ("validation", "test", "train"):
            row = cohort.rows[cohort.part_neurons[part][0]]
            fold = session.fold_for(row)
            reconstructor, module, state, adapter = session.layer(fold)
            answer = session.infer(row)
            p = answer.provenance
            assert answer.reason is None and p["used"] and p["level"] == level and p["conditional"] == "T1_profile_from_soma_label"
            assert set(answer.query.inputs) == set(session.measured_records(row)) == set(core.HIDDEN_AXON_INPUTS)
            assert set(answer.query.inferred) == {TARGET} and answer.query.unit_id == row.neuron and answer.query.group_id == row.animal
            stated = answer.query.inferred[TARGET]
            assert (stated["scale"], stated["transform"]) == ("log10(1 + value)", "log10_1p")
            assert stated["frame"] == session.registry.variable("jiao", TARGET).frame
            assert answer.profile["names"] == tuple(session.world.regions) and len(answer.profile["value"]) == len(session.world.regions)
            assert answer.query.inputs["soma_position"]["value"] == tuple(float(v) for v in row.soma)
            assert answer.query.inputs["source_label"]["value"] == row.tag and answer.query.inputs["hemisphere"]["value"] == row.hemisphere
            assert p["state_sha256"] == state.state_hash() and row.animal not in answer.upstream_groups
            assert tuple(answer.upstream_groups) == tuple(state.groups) and p["run_complete_sha256"] == session.run_complete_sha256
            if level == "engine":
                assert p["value_state_hash"] == p["state_sha256"] and p["comparator"] is None and p["module"] == B.MODULE
            else:
                assert p["comparator"] in mods.COMPARATORS and p["module"] == "baseline:" + p["comparator"]
                assert "no_better_than_baseline" in p["codes"] and p["value_state_hash"] != p["state_sha256"]
            # the adapter: exactly as blocked as before; through the request API the inputs pass and the arbor is withheld
            assert p["adapter_status"] == BLOCKED == adapter.STATUS == adapter.evidence_for(adapter.conditionals()[0].name).status
            arbor = reconstructor.reconstruct(session._request(row, [*session.measured_records(row).values(), answer.record],
                                                               "axon_arbor"))["axon_arbor"]
            assert arbor.status == "refused" and arbor.refusal.code == "not_validated" and arbor.refusal.reason == BLOCKED
            assert arbor.estimate is None
            with pytest.raises(core.ReconError, match="blocked"):
                adapter.estimate(None, 0.9)
        assert reconstructor.policy == core.Policy()                      # the layer's default policy, nothing else
    source = (ROOT / B.BRIDGE_FILE).read_text()
    assert "import axon_engine" not in source and "backend=" not in source and ".backend" not in source     # nothing is handed to a backend


# ------------------------------------------------------------------------------------------
# K2 in context: scale, tag and frame; one conversion
# ------------------------------------------------------------------------------------------


def independent_weights(directory, engine, neuron, donor_ids, *, convert_p_again=False, convert_q_again=False, e2=False):
    """E7's weights recomputed from the tables, with every value brought to log10(1 + length) exactly once: the
    inferred profile arrives on that scale already; the donors' published lengths are transformed here, once."""
    cohort = engine.cohort
    row = cohort.rows[neuron]
    answer = engine.layer_session().infer(row)
    assert answer.profile["transform"] == "log10_1p"
    p = dict(zip(answer.profile["names"], answer.profile["value"]))
    rows, names = regional_rows(directory)
    spec = E.CONFIG["engine"]
    library = [n for n in engine.library["ids"] if cohort.rows[n].animal != row.animal or row.part != "train"]
    distance = {n: float(np.linalg.norm(np.array(cohort.rows[n].soma) - np.array(row.soma))) for n in library}
    neighbours = sorted(sorted(library, key=lambda n: (distance[n], n))[:spec["k_neighbours"]])
    assert neighbours == list(donor_ids)
    d = np.array([distance[n] for n in neighbours])
    animals = [cohort.rows[n].animal for n in neighbours]
    prior = np.array([1.0 / animals.count(a) for a in animals]) / len(set(animals))
    h = max(spec["kernel_h0_um"], float(np.sort(d)[spec["kernel_rank"] - 1]))
    kernel = np.exp(-d ** 2 / (2.0 * h * h))
    target = np.array([p[name] for name in names])
    if convert_p_again:
        target = np.log10(1.0 + target)
    Q = np.array([np.log10(1.0 + rows[cohort.rows[n].full_id]) for n in neighbours])
    if convert_q_again:
        Q = np.log10(1.0 + Q)
    gap = np.sqrt(((Q - target[None, :]) ** 2).sum(axis=1))
    tau = float(np.median(gap))
    s = np.ones(len(gap)) if tau == 0 else np.exp(-gap ** 2 / (2.0 * tau * tau))
    raw = prior * kernel * (np.ones(len(gap)) if e2 else s)
    return raw / raw.sum(), s


def test_K2a_the_weights_come_from_values_converted_exactly_once(layered):
    engine = engine_of(layered, "validate")
    engine.prepare_prediction()
    checked = 0
    for neuron in engine.cohort.part_neurons["validation"][:6]:
        forecasts, _ = forecasts_of(engine, neuron)
        details = forecasts[E7]["details"]
        if not details["inferred_profile_used"]:
            continue
        mine, similarity = independent_weights(layered, engine, neuron, forecasts[E7]["donor_ids"])
        assert np.allclose(details["weights"], mine, rtol=1e-12, atol=0.0)
        assert np.allclose(details["similarity"], similarity, rtol=1e-12, atol=0.0)
        # review finding 7: E2's rows carry no weights, so E2's weights are recomputed independently (s_j = 1) and the
        # similarity must really move E7's weights away from them
        e2_weights, _ = independent_weights(layered, engine, neuron, forecasts[E7]["donor_ids"], e2=True)
        assert "weights" not in forecasts[E2]["details"] and float(np.abs(similarity - 1.0).max()) > 1e-3
        assert float(np.abs(np.asarray(details["weights"]) - e2_weights).max()) > 1e-6
        # a deliberate second conversion, of the inferred profile or of the donor profiles, gives other weights
        for again in ({"convert_p_again": True}, {"convert_q_again": True}):
            twice, _ = independent_weights(layered, engine, neuron, forecasts[E7]["donor_ids"], **again)
            assert float(np.abs(twice - mine).max()) > 1e-4, again
        checked += 1
    assert checked >= 3
    # the same through the bridge's own function: donors given as published lengths (tag 'identity') are transformed
    # once; donors given on the model scale (tag 'log10_1p', what the engine's sealed loader stores) not at all
    neuron = engine.cohort.part_neurons["validation"][0]
    forecasts, _ = forecasts_of(engine, neuron)
    answer = engine.layer_session().infer(engine.cohort.rows[neuron])
    rows, names = regional_rows(layered)
    published = [rows[engine.cohort.rows[n].full_id] for n in forecasts[E7]["donor_ids"]]
    stored = [engine.profiles[engine.cohort.rows[n].full_id] for n in forecasts[E7]["donor_ids"]]
    once = B.profile_similarity(answer.profile, published, names, donor_tag="identity")
    kept = B.profile_similarity(answer.profile, stored, E.REGIONS, donor_tag="log10_1p")
    assert np.array_equal(once["s"], kept["s"]) and once["tau"] == kept["tau"]
    assert np.allclose(kept["s"], forecasts[E7]["details"]["similarity"], rtol=0, atol=0)
    twice = B.profile_similarity(answer.profile, stored, E.REGIONS, donor_tag="identity")      # a wrong tag: converted again
    assert float(np.abs(twice["s"] - kept["s"]).max()) > 1e-4
    with pytest.raises(core.ScaleError, match="a length is not negative|negative value|below zero"):
        B.profile_similarity(answer.profile, [-np.ones(len(E.REGIONS))], E.REGIONS)


def test_K2b_a_wrong_frame_other_units_and_a_contradictory_scale_statement_are_refused(planted, layered):
    session, cohort = session_of(planted)
    row = cohort.rows[cohort.part_neurons["validation"][0]]
    reconstructor = session.layer(None)[0]
    measured = session.measured_records(row)
    soma, label = measured["soma_position"], measured["source_label"]
    smuggled = dataclasses.replace(soma)
    object.__setattr__(smuggled, "transform", "log10_1p")              # the record class itself refuses to build this one
    with pytest.raises(ValueError, match="contradicts the value tag"):
        dataclasses.replace(soma, transform="log10_1p")
    right_side = tuple(soma.value[:2]) + (11400.0 - soma.value[2],)       # true right-hemisphere coordinates, not mirrored
    cases = ((dataclasses.replace(soma, frame=None), "no frame is stated"),
             (dataclasses.replace(soma, frame="another atlas, left-right flipped"), "not the registered frame"),
             (dataclasses.replace(soma, scale="linear position; mm"), "is not the registered statement"),
             (dataclasses.replace(soma, value=right_side), "lies outside the registered frame"),
             (smuggled, "contradicts the value tag"))
    for bad, text in cases:
        request, refusal = session.ask(row, None, records={"soma_position": bad, "source_label": label})
        assert request is None and refusal.code in ("invalid_request", "forbidden_input") and text in refusal.reason, text
        query, why = session.admit(row, None, {**measured, "soma_position": bad}, None)
        assert query is None and text in why
        direct = reconstructor.reconstruct(session._request(row, [bad, label], TARGET))[TARGET]       # the layer itself says the same
        assert direct.status == "refused" and direct.refusal.code == "invalid_request" and text in direct.refusal.reason
    # the inferred PROFILE, at the point where E7 would use it: other units, a tag that contradicts the text, no frame
    answer = session.infer(row)
    profile, names = answer.profile, answer.profile["names"]
    donors = [np.full(len(names), 2.0), np.full(len(names), 1.0)]
    assert B.profile_similarity(profile, donors, names)["regions"] == len(names)
    spec = session.registry.variable("jiao", TARGET)
    for change, text in (({"scale": "length >= 0 per region; mm"}, "is not the registered statement"),
                         ({"transform": "identity"}, "is not the registered statement"),
                         ({"scale": core.scale_statement(spec, "identity")}, "is not the registered statement"),
                         ({"frame": None}, "no frame is stated"),
                         ({"transform": "ln"}, "no registered scale goes with the value tag")):
        with pytest.raises(core.ScaleError, match=text):
            B.profile_similarity({**profile, **change}, donors, names)
    with pytest.raises(core.ScaleError):                                 # donors that claim another scale than they are in
        B.profile_similarity(profile, donors, names, donor_tag="ln")
    # in the engine a profile that does not state the registered scale stops the forecast; it is not silently dropped
    engine = engine_of(layered, "validate")
    engine.prepare_prediction()
    live = engine.layer_session()
    real = live.infer

    def in_other_units(row, **kw):
        answer = real(row, **kw)
        answer.profile = {**answer.profile, "scale": "length >= 0 per region; mm"}
        return answer

    live.infer = in_other_units
    with pytest.raises(E.EngineError, match="is not the registered statement"):
        engine.forecast(engine.cohort.part_neurons["validation"][0])


def test_K2c_regions_are_matched_by_name_never_by_position(layered):
    engine = engine_of(layered, "validate")
    engine.prepare_prediction()
    neuron = engine.cohort.part_neurons["validation"][0]
    reference, _ = forecasts_of(engine, neuron)
    answer = engine.layer_session().infer(engine.cohort.rows[neuron])
    names, values = list(answer.profile["names"]), list(answer.profile["value"])
    donors = np.array([engine.profiles[engine.cohort.rows[n].full_id] for n in reference[E7]["donor_ids"]])
    base = B.profile_similarity(answer.profile, donors, E.REGIONS)
    order = np.random.default_rng(1).permutation(len(names))
    assert not np.array_equal(order, np.arange(len(names)))
    permuted = {**answer.profile, "names": tuple(names[i] for i in order), "value": tuple(values[i] for i in order)}
    one = B.profile_similarity(permuted, donors, E.REGIONS)                                 # the inferred side reordered
    two = B.profile_similarity(answer.profile, donors[:, order], [E.REGIONS[i] for i in order])     # the donor side reordered
    assert np.array_equal(one["s"], base["s"]) and np.array_equal(two["s"], base["s"]) and one["tau"] == two["tau"] == base["tau"]
    by_position = B.profile_similarity({**answer.profile, "value": tuple(values[i] for i in order)}, donors, E.REGIONS)
    assert not np.allclose(by_position["s"], base["s"])                  # values moved without their names: another answer
    with pytest.raises(B.BridgeRefused, match="do not hold the same regions by name"):          # a region missing on one side
        B.profile_similarity({**answer.profile, "names": tuple(names[:-1]), "value": tuple(values[:-1])}, donors, E.REGIONS)
    with pytest.raises(core.ScaleError, match="not coordinates of the registered variable"):   # a name that is no region at all
        B.profile_similarity({**answer.profile, "names": tuple(names[:-1]) + ("NotARegion",)}, donors, E.REGIONS)
    with pytest.raises(B.BridgeRefused, match="duplicated"):
        B.profile_similarity({**answer.profile, "names": tuple(names[:-1]) + (names[0],)}, donors, E.REGIONS)
    with pytest.raises(B.BridgeRefused, match="values for"):
        B.profile_similarity({**answer.profile, "names": tuple(names[:-1])}, donors, E.REGIONS)
    with pytest.raises(B.BridgeRefused, match="duplicated"):
        B.profile_similarity(answer.profile, donors, list(E.REGIONS[:-1]) + [E.REGIONS[0]])
    # the whole forecast with the layer's answer handed over in another column order: the same weights and arrays
    live = engine.layer_session()
    real = live.infer

    def reordered(row, **kw):
        answer = real(row, **kw)
        names, values = answer.profile["names"], answer.profile["value"]
        answer.profile = {**answer.profile, "names": tuple(names[i] for i in order), "value": tuple(values[i] for i in order)}
        return answer

    live.infer = reordered
    again, _ = forecasts_of(engine, neuron)
    assert again[E7]["details"]["weights"] == reference[E7]["details"]["weights"]
    assert E.arbor_sha(again[E7]["arbor"]) == E.arbor_sha(reference[E7]["arbor"])


# ------------------------------------------------------------------------------------------
# K3 in context: the states and the files of the finished run
# ------------------------------------------------------------------------------------------


def state_tampers(heldout_animal):
    """One answer-affecting item of a saved state changed at a time."""
    AP = "soma_position:AP"

    def weights(h, a):
        w = a["a0"].copy()
        w[0] += 1e-9
        a["a0"] = w

    return (("a support range", lambda h, a: h["ranges"].__setitem__(AP, [h["ranges"][AP][0] - 5000.0, h["ranges"][AP][1] + 5000.0])),
            ("the support margin (a reader setting)", lambda h, a: h["reader"].__setitem__("support_margin", 0.5)),
            ("a fit setting", lambda h, a: h["settings"].__setitem__("seed", h["settings"]["seed"] + 1)),
            ("a parameter", weights),
            ("the training groups (one removed)", lambda h, a: h["groups"].pop(0)),
            ("the training groups (a held-out brain sample added)", lambda h, a: h["groups"].append(heldout_animal)))


def test_K3a_a_changed_state_file_makes_the_bridge_refuse_to_build(layered, tmp_path):
    directory = clone(layered, tmp_path / "w")
    world = E.synthetic_world(directory)
    cohort = E.Cohort(world)
    out = layer_out(directory)
    pins = B.run_pins(world, folds="all")["run"]["states"]
    held_row = cohort.rows[cohort.part_neurons["validation"][0]]
    train_row = cohort.rows[cohort.part_neurons["train"][0]]
    fold = B.open_session(world, heldout_groups=cohort.heldout_animals).fold_for(train_row)
    for key, row in (("serving", held_row), (f"outer{fold}", train_row)):
        path = out / pins[key]["file"]
        original = path.read_bytes()
        assert B.open_session(world, heldout_groups=cohort.heldout_animals).infer(row).profile is not None
        for what, change in state_tampers(cohort.part_animals["test"][0]):
            for rehash, text in ((False, "does not reproduce its recorded hash"), (True, "not the expected state")):
                head, arrays = read_state(path)
                change(head, arrays)
                write_state(path, head, arrays, rehash=rehash)
                session = B.open_session(world, heldout_groups=cohort.heldout_animals)
                with pytest.raises(B.BridgeRefused, match=text) as refused:
                    session.infer(row)
                assert "refuses to serve the finished run" in str(refused.value), (key, what)
                with pytest.raises(E.EngineError, match="general reconstruction layer"):       # the engine: no forecast at all
                    engine_of(directory, "validate" if key == "serving" else "develop").forecast(row.neuron)
                path.write_bytes(original)
        assert B.open_session(world, heldout_groups=cohort.heldout_animals).infer(row).profile is not None


def test_K3b_an_edited_evidence_or_selection_file_of_the_finished_run_is_refused(layered, tmp_path):
    directory = clone(layered, tmp_path / "w")
    world = E.synthetic_world(directory)
    out = layer_out(directory)
    assert B.open_session(world).ready
    for name, edit in (("evidence.json", lambda d: d["tasks"]["T1"].__setitem__("status", "exploratory_supported")),
                       ("selection.json", lambda d: d["serving"]["T1"].__setitem__("components", 1)),
                       ("serving.json", lambda d: d["tasks"]["T1"].__setitem__("state_hash", "0" * 64)),
                       ("metrics.json", lambda d: d.__setitem__("exploratory", False)),
                       ("pre_run_declaration.json", lambda d: d["config"].__setitem__("min_relative_gain", 0.0))):
        original = (out / name).read_bytes()
        record = load(out / name)
        edit(record)
        (out / name).write_text(core.canonical(record))
        with pytest.raises(B.BridgeRefused, match="does not verify"):
            B.open_session(world)
        with pytest.raises(B.BridgeRefused, match="does not verify"):
            B.run_pins(world, folds="all")
        with pytest.raises(E.EngineError, match="general reconstruction layer"):
            engine_of(directory, "validate").layer_session()
        (out / name).write_bytes(original)
    assert B.open_session(world).ready
    # a converted table edited after the conversion: the layer's declaration no longer matches its inputs
    regional = B.layer_directory(world) / "regional.csv"
    regional.write_text(regional.read_text().replace(".00,", ".50,", 1))
    with pytest.raises(B.BridgeRefused, match="not the file the converter wrote"):
        B.open_session(world)


def test_K3c_every_e7_forecast_row_names_the_state_the_finished_run_recorded(layered, tmp_path):
    directory = clone(layered, tmp_path / "w")
    out = layer_out(directory)
    serving, selection = bound_serving(directory), load(out / "selection.json")
    assert load(out / "serving.json") == {**serving, "run_complete_sha256": E.sha_file(out / "run_complete.json")}      # a copy, nothing more
    session, cohort = session_of(directory)
    run_sha = E.sha_file(out / "run_complete.json")
    for stage in ("develop", "validate"):
        engine_of(directory, stage).predict()
        manifest = load(directory / "engine" / stage / "prediction_manifest.json")
        rows = [p for p in manifest["predictions"] if p["model"] == E7]
        assert len(rows) == len(manifest["queries"]) and all("general_layer" not in p for p in manifest["predictions"] if p["model"] != E7)
        for p in rows:
            layer = p["general_layer"]
            row = cohort.rows[p["query_id"]]
            if stage == "validate":
                assert layer["fold"] == "serving" and layer["state_sha256"] == serving["tasks"]["T1"]["state_hash"]
            else:
                fold = session.fold_of_group[row.animal]
                assert layer["fold"] == f"outer{fold}" and layer["state_sha256"] == selection["engine_used"]["T1"][str(fold)]["state_hash"]
            assert layer["run_complete_sha256"] == run_sha == manifest["pins"]["general_layer"]["run"]["run_complete_sha256"]
            assert manifest["pins"]["general_layer"]["run"]["states"][layer["fold"]]["state_hash"] == layer["state_sha256"]
            upstream = manifest["upstream_lists"][p["upstream_fit_animals_sha256"]]      # also on fallback rows
            assert row.animal not in upstream and not set(upstream) & cohort.heldout_animals
            assert upstream == sorted(session.layer(None if stage == "validate" else session.fold_of_group[row.animal])[2].groups)
            assert layer["adapter_status"] == BLOCKED and p["regional_state_id"] is None
            assert len(json.dumps(layer)) < 2000                         # small fields in the row; lists and file hashes by hash
            assert p["donor_ids"] == sorted(p["donor_ids"]) and len(p["donor_ids"]) == E.CONFIG["engine"]["k_neighbours"]
            if layer["used"]:
                assert (layer["level"] == "engine") == (layer["value_state_hash"] == layer["state_sha256"])
            else:
                assert p["fallback"] and f"no_inferred_profile:{layer['reason_code']}" in p["fallback_reason"]
        use = manifest["general_layer"]["use_by_E7"]
        assert use["forecasts"] == len(rows) == use["engine_level_value"] + use["baseline_level_value"] + use["fallback_to_E2_weights"]
        assert manifest["general_layer"]["registration"]["run_complete_sha256"] == run_sha
        assert set(manifest["general_layer"]["upstream_lists"]) <= set(manifest["upstream_lists"])
        assert set(manifest["pins"]["code_sha256"]["general_layer"]) == set(E.LAYER_CODE_FILES)
        assert E.BRIDGE_FILE in manifest["pins"]["code_sha256"]["engine"]
    state_manifest = load(directory / "engine/state_manifest.json")
    assert state_manifest["general_layer"]["run_complete_sha256"] == run_sha
    assert state_manifest["general_layer"]["state_fingerprints"]["serving"] == serving["tasks"]["T1"]["state_hash"]


# ------------------------------------------------------------------------------------------
# hidden outcomes
# ------------------------------------------------------------------------------------------


def rebuilt(directory):
    """Everything derived from a world's tables and traces, made from scratch: library, pair cache, regional states,
    the general layer's run."""
    engine = engine_of(directory, "develop")
    engine.build_library()
    engine.build_cache()
    engine.states.build_many(list(engine.states.plan))
    engine.prepare_general_layer()
    return directory


def held_out_hashes(directory):
    out = {}
    for stage, part in (("validate", "validation"), ("test", "test")):
        engine = engine_of(directory, stage)
        engine.prepare_prediction()
        out.update(array_hashes(engine, engine.cohort.part_neurons[part]))
    return out


def test_hidden_outcomes_scrambled_or_deleted_change_no_forecast_of_any_method(layered, tmp_path):
    reference = held_out_hashes(layered)
    assert {m for _, m in reference} == set(E.METHODS) and E7 in E.METHODS
    session, _ = session_of(layered)
    assert session.ready                                                 # E7 really used the layer here
    copy = tmp_path / "scrambled"
    shutil.copytree(layered / "data", copy / "data")
    cohort = E.Cohort(E.synthetic_world(layered))
    held = {n for part in ("validation", "test") for n in cohort.part_neurons[part]}
    assert set(E.scramble_outcomes(copy)) == held                        # every held-out trace and regional row: other values
    truth = load(layered / "generator_truth.json")
    order = sorted(held)
    shuffled = dict(zip(order, [truth[n] for n in order[1:] + order[:1]]))
    E.write_json(copy / "generator_truth.json", {**truth, **shuffled}, pretty=True)      # the generator's truth, shuffled too
    assert (copy / "data/regional.csv").read_text() != (layered / "data/regional.csv").read_text()
    rebuilt(copy)
    again = held_out_hashes(copy)
    assert again == reference and len(again) == len(held) * len(E.METHODS)
    # the held-out traces do not exist at all, and the held-out regional rows cannot be parsed: still identical
    gone = tmp_path / "deleted"
    shutil.copytree(layered / "data", gone / "data")
    for part in ("validation", "test"):
        for trace in (gone / "data/traces" / part).rglob("*.swc"):
            trace.unlink()
    held_full = {cohort.rows[n].full_id for n in held}
    lines = (gone / "data/regional.csv").read_text().splitlines()
    (gone / "data/regional.csv").write_text("\n".join(
        l if l.split(",", 1)[0] not in held_full else l.split(",", 1)[0] + ",SEALED,NEVER,PARSED" for l in lines) + "\n")
    rebuilt(gone)
    assert held_out_hashes(gone) == reference
    assert load(B.layer_directory(E.synthetic_world(gone)) / "conversion.json")["training_content_sha256"] == \
        load(B.layer_directory(E.synthetic_world(layered)) / "conversion.json")["training_content_sha256"]


def test_one_training_animals_regional_rows_scrambled_nothing_crashes_and_only_dependent_forecasts_change(layered, tmp_path):
    """Development, leave one animal out. The regional rows of ONE training animal are replaced and everything is
    rebuilt. Methods without a fitted state do not move. A forecast may change where its fitted state or its
    neighbourhood used that animal. The animal's own queries: their regional state and their general-layer state
    never saw the animal and no donor is its own, so their forecasts stay, with one stated exception for E7: which
    kind of value the layer issues (its engine's or a comparator's) is a verdict of the layer's benchmark over ALL
    training samples, the query's own included; where that verdict changes, E7 may change."""
    before = engine_of(layered, "develop")
    before.prepare_prediction()
    neurons = before.queries()
    reference = {n: forecasts_of(before, n)[0] for n in neurons}
    animal = before.cohort.rows[neurons[0]].animal
    copy = tmp_path / "w"
    shutil.copytree(layered / "data", copy / "data")
    assert scramble_regional_rows(copy, animal) > 0
    rebuilt(copy)                                                        # nothing crashes: states refitted, the layer run again
    after = engine_of(copy, "develop")
    after.prepare_prediction()
    again = {n: forecasts_of(after, n)[0] for n in neurons}
    stateless = [m for m in E.METHODS if m not in E.USES_FITTED_STATE]
    same = lambda n, m: E.arbor_sha(again[n][m]["arbor"]) == E.arbor_sha(reference[n][m]["arbor"])     # noqa: E731
    assert all(same(n, m) for n in neurons for m in stateless)
    own = [n for n in neurons if before.cohort.rows[n].animal == animal]
    assert own and all(same(n, m) for n in own for m in E.USES_REGIONAL_STATE)       # its own regional state never saw it
    evidence = lambda d: {k: load(layer_out(d) / "evidence.json")["tasks"]["T1"][k]                    # noqa: E731
                          for k in ("status", "best_comparator", "best_comparator_params")}
    for n in own:
        a, b = reference[n][E7]["general_layer"], again[n][E7]["general_layer"]
        assert a["fold"] == b["fold"] and a["state_sha256"] == b["state_sha256"]     # the state that never saw the animal: the same state
    if evidence(layered) == evidence(copy):
        assert all(same(n, E7) for n in own)
        assert all(reference[n][E7]["details"]["weights"] == again[n][E7]["details"]["weights"] for n in own)
    # queries of other animals with a donor of that animal among their neighbours: the measured donor profile changed,
    # so E7's similarity factor moved (whatever the template choice did)
    touched = [n for n in neurons if n not in own and reference[n][E7]["details"]["inferred_profile_used"]
               and again[n][E7]["details"]["inferred_profile_used"]
               and any(before.cohort.rows[d].animal == animal for d in reference[n][E7]["donor_ids"])]
    assert touched and all(reference[n][E7]["details"]["similarity"] != again[n][E7]["details"]["similarity"] for n in touched)


# ------------------------------------------------------------------------------------------
# planted leaks
# ------------------------------------------------------------------------------------------


def test_planted_leak_a_held_out_animal_inside_a_serving_state_is_refused(layered, tmp_path):
    directory = clone(layered, tmp_path / "w")
    session, cohort = session_of(directory)
    row = cohort.rows[cohort.part_neurons["validation"][0]]
    intruder = cohort.part_animals["test"][0]
    reconstructor, module, state, adapter = session.layer(None)
    assert session.infer(row).profile is not None
    # (a) planted in memory, behind every file check: the bridge's own test of the state refuses it
    twin = mods.JointState(layout=state.layout, encoding=state.encoding, field_meta=state.field_meta, blocks=state.blocks,
                           settings=state.settings, units=state.units, groups=tuple(state.groups) + (intruder,),
                           ranges=state.ranges, support_margin=state.support_margin, quantile_iterations=state.quantile_iterations)
    session._layers[None] = (reconstructor, module, twin, adapter)
    with pytest.raises(B.BridgeLeak, match="holds validation or test brain samples"):
        session.infer(row)
    # (b) planted on disk CONSISTENTLY: the state file with its hash rewritten, the completion record rewritten to
    #     expect that state, and serving.json made the exact copy again. Every file check of the layer passes.
    engine_of(directory, "validate").layer_pins()                        # the engine has registered the genuine run
    out = layer_out(directory)
    world = E.synthetic_world(directory)
    path = out / B.run_pins(world, folds=[None])["run"]["states"]["serving"]["file"]
    head, arrays = read_state(path)
    head["groups"] = sorted(head["groups"] + [intruder])
    write_state(path, head, arrays, rehash=True)
    forged, _ = read_state(path)
    complete = load(out / "run_complete.json")
    genuine = complete[bench.SERVING_KEY]["tasks"]["T1"]["state_hash"]
    for entry in complete[bench.SERVING_KEY]["tasks"].values():          # every task that is served by this state file
        if entry["state_hash"] == genuine:
            entry["state_hash"], entry["training_groups_sha256"] = forged["state_hash"], core.ids_sha(forged["groups"])
    (out / "run_complete.json").write_text(core.canonical(complete))
    (out / "serving.json").write_text(core.canonical(bench.serving_file_content(B.layer_world_for(world))))
    assert bench.verify_run(B.layer_world_for(world))                    # the layer's file checks: passed
    with pytest.raises(B.BridgeRefused) as refused:                      # the layer's module: not the rows the state was fitted on
        B.open_session(world, heldout_groups=cohort.heldout_animals).infer(row)
    assert "comparator rows are not the rows" in str(refused.value) or "cannot be shown to exclude" in str(refused.value) \
        or "holds validation or test" in str(refused.value)
    with pytest.raises(E.EngineError, match="is not the one registered in this world's state manifest"):
        engine_of(directory, "validate").layer_pins()                    # the engine: not the completion record it registered
    with pytest.raises(E.EngineError):
        engine_of(directory, "validate").predict()
    assert not (directory / "engine/validate/PREDICTION_COMMITMENT.json").exists()
    # (c) a held-out brain sample in the layer's training partition is caught when the session opens
    with pytest.raises(B.BridgeLeak, match="training partition holds validation or test brain samples"):
        B.open_session(E.synthetic_world(layered), heldout_groups=set(cohort.heldout_animals) | {cohort.part_animals["train"][0]})


def test_planted_leak_the_query_animal_inside_its_own_state_is_refused(layered, tmp_path):
    directory = clone(layered, tmp_path / "w")
    engine = engine_of(directory, "develop")
    engine.prepare_prediction()
    neuron = engine.queries()[0]
    row = engine.cohort.rows[neuron]
    session = engine.layer_session()
    k = session.fold_for(row)
    assert len(engine.forecast(neuron)[0]) == len(E.METHODS)             # clean: forecasts are produced
    wrong = next(f for f in sorted(set(session.fold_of_group.values())) if f != k)
    session.fold_for = lambda row: wrong                                 # the bridge is made to pick a state that saw the animal
    with pytest.raises(E.LeakError, match="its brain sample .* is in the T1 state"):
        engine.forecast(neuron)
    session.fold_for = lambda row: None                                  # ... or the all-training state
    with pytest.raises(E.LeakError, match="is in the T1 state"):
        engine.forecast(neuron)
    # the engine's own second check, should a bridge ever hand back such a state without complaining
    fresh = engine_of(directory, "develop")
    fresh.prepare_prediction()
    live = fresh.layer_session()
    real = live.infer

    def lying(row, **kw):
        answer = real(row, **kw)
        answer.upstream_groups = tuple(answer.upstream_groups) + (row.animal,)
        return answer

    live.infer = lying
    with pytest.raises(E.LeakError, match="holds its animal or a held-out animal"):
        fresh.forecast(neuron)
    # and the manifest audit: an E7 row whose upstream list holds the query's animal
    engine_of(directory, "develop").predict()
    auditor = engine_of(directory, "develop")
    manifest = load(directory / "engine/develop/prediction_manifest.json")
    assert auditor.audit_manifest(manifest)
    victim = next(p for p in manifest["predictions"] if p["model"] == E7 and p["upstream_fit_animals_sha256"])
    planted = sorted(manifest["upstream_lists"][victim["upstream_fit_animals_sha256"]] + [victim["animal"]])
    manifest["upstream_lists"][E.ids_sha(planted)] = planted
    victim["upstream_fit_animals_sha256"] = E.ids_sha(planted)
    with pytest.raises(E.LeakError, match="holds the query's animal or a held-out animal"):
        auditor.audit_manifest(manifest)
    manifest = load(directory / "engine/develop/prediction_manifest.json")
    victim = next(p for p in manifest["predictions"] if p["model"] == E7)
    del victim["general_layer"]
    with pytest.raises(E.EngineError, match="fitted-state provenance is not recorded"):
        auditor.audit_manifest(manifest)


def test_planted_leak_a_held_out_regional_row_handed_to_the_bridge_is_refused(layered):
    session, cohort = session_of(layered)
    row = cohort.rows[cohort.part_neurons["test"][0]]
    line = next(l for l in (layered / "data/regional.csv").read_text().splitlines() if l.startswith(row.full_id + ","))
    values = [float(v) for v in line.split(",")[1:]]                    # the held-out row itself (a synthetic world)
    measured = session.measured_records(row)
    stolen = profile_record(session, row, E.REGIONS, values)
    request, refusal = session.ask(row, None, records={"soma_position": measured["soma_position"],
                                                       "source_label": measured["source_label"], TARGET: stolen})
    assert request is None and refusal.code == "forbidden_input" and "never an input to a hidden-axon forecast" in refusal.reason
    query, why = session.admit(row, None, {**measured, TARGET: stolen}, None)
    assert query is None and "never an input to a hidden-axon forecast" in why
    # relabelled as 'inferred' by hand: it was never issued by the layer
    spec = session.registry.variable("jiao", TARGET)
    state = session.layer(None)[2]
    evidence = (("soma_position", "measured"), ("source_label", "measured"))
    disguised = core.Observation(
        session.world.species, "jiao", "neuron", row.neuron, session.group_of(row), TARGET, "inferred",
        tuple(float(v) for v in np.log10(1.0 + np.array(values))), core.scale_statement(spec, "log10_1p"), tuple(E.REGIONS), spec.frame,
        core.Uncertainty("inferred", None, None, None, None, B.MODULE, state.state_hash(), state.groups_sha256, "unvalidated"),
        core.SourceRef("reconstruct:" + B.MODULE, state.state_hash(), row.neuron), derived_from=evidence, transform="log10_1p",
        lineage=(core.LineageStep(TARGET, B.MODULE, state.state_hash(), state.groups_sha256, evidence,
                                  ("jiao", "neuron", row.neuron), row.animal),))
    query, why = session.admit(row, None, measured, disguised)
    assert query is None and "not issued by this reconstruction layer" in why
    # the bridge's own interface takes a cohort row and nothing else of the neuron
    assert set(session.measured_records(row)) == set(core.HIDDEN_AXON_INPUTS)
    assert "profile" not in E.NeuronRow.__dataclass_fields__


def test_planted_leak_a_forged_run_complete_is_refused(layered, tmp_path):
    """The reviewer's forgery in this context: the evidence is edited so that the layer's engine looks supported, and
    run_complete.json is rewritten to bind the edited file. The layer's own file check then passes; what refuses is
    the binding of the serving states to the completion record, the engine's registration of the run, and (after
    forecasts were committed) the engine's pins."""
    directory = clone(layered, tmp_path / "w")
    world = E.synthetic_world(directory)
    cohort = E.Cohort(world)
    out = layer_out(directory)
    engine_of(directory, "validate").predict()                           # forecasts committed against the genuine run
    genuine = E.sha_file(out / "run_complete.json")
    evidence = load(out / "evidence.json")
    evidence["tasks"]["T1"]["status"] = "exploratory_supported"
    (out / "evidence.json").write_text(core.canonical(evidence))
    record = load(out / "run_complete.json")
    record["sha256"]["evidence.json"] = record["evidence_sha256"] = core.sha_file(out / "evidence.json")
    (out / "run_complete.json").write_text(core.canonical(record))
    row = cohort.rows[cohort.part_neurons["validation"][0]]
    with pytest.raises(B.BridgeRefused, match="serving.json is not the record"):       # the copy still names the genuine record
        B.open_session(world, heldout_groups=cohort.heldout_animals)
    (out / "serving.json").write_text(core.canonical(bench.serving_file_content(B.layer_world_for(world))))
    assert bench.verify_run(B.layer_world_for(world))                    # now every file check of the layer passes
    assert E.sha_file(out / "run_complete.json") != genuine
    # without an anchor the forgery is believed: the layer now answers with its engine, as the edited evidence says.
    # This is the limit the layer states itself; the anchor is the SHA-256 of the completion record
    believed = B.open_session(world, heldout_groups=cohort.heldout_animals).infer(row)
    assert believed.provenance["level"] == "engine"
    with pytest.raises(B.BridgeRefused, match="is not the registered one"):            # the bridge, given the registered SHA-256
        B.open_session(world, heldout_groups=cohort.heldout_animals, expected_run_complete_sha256=genuine)
    with pytest.raises(E.EngineError, match="is not the one registered in this world's state manifest"):
        engine_of(directory, "validate").layer_pins()                    # the engine: it registered the genuine run
    with pytest.raises(E.EngineError, match="is not the one registered"):
        engine_of(directory, "validate").forecast(row.neuron)
    with pytest.raises(E.EngineError, match="method changed since these forecasts were committed"):
        engine_of(directory, "validate").admit()
    assert not (directory / "engine/validate/TRUTH_OPENED.json").exists()
    # the run removed altogether after it was registered
    shutil.rmtree(B.layer_directory(world))
    with pytest.raises(E.EngineError, match="no finished run is on disk now"):
        engine_of(directory, "validate").layer_pins()


def forge_run(directory):
    """The provenance attacker's forgery: evidence edited so that the layer's engine looks supported, run_complete.json and
    serving.json rewritten so that every file check of the layer passes."""
    world = E.synthetic_world(directory)
    out = layer_out(directory)
    evidence = load(out / "evidence.json")
    evidence["tasks"]["T1"]["status"] = "exploratory_supported"
    (out / "evidence.json").write_text(core.canonical(evidence))
    record = load(out / "run_complete.json")
    record["sha256"]["evidence.json"] = record["evidence_sha256"] = core.sha_file(out / "evidence.json")
    (out / "run_complete.json").write_text(core.canonical(record))
    (out / "serving.json").write_text(core.canonical(bench.serving_file_content(B.layer_world_for(world))))
    assert bench.verify_run(B.layer_world_for(world))
    return E.sha_file(out / "run_complete.json")


def test_provenance_a_forged_run_with_the_registration_deleted_is_never_registered_at_a_held_out_stage(layered, tmp_path):
    """The provenance attacker's finding: the run was forged consistently and the engine's registration file deleted. A
    held-out stage never registers a run: it refuses, and nothing is forecast or committed. Registration happens only
    in the development stage, in `layer` and when the synthetic run is built."""
    directory = clone(layered, tmp_path / "w")
    genuine = load(directory / "engine/states/general_layer.json")["run_complete_sha256"]
    forged = forge_run(directory)
    assert forged != genuine
    (directory / "engine/states/general_layer.json").unlink()
    # round 2: the world's anchor (written by the build step) refuses the forged run first; the rule that a held-out stage
    # never registers a run is shown on its own just below, with the GENUINE run and no registration
    for stage in ("validate", "test"):
        with pytest.raises(E.EngineError, match="no registered finished run|is not the run anchored for this world"):
            engine_of(directory, stage).layer_pins()
    with pytest.raises(E.EngineError, match="no registered finished run|is not the run anchored for this world"):
        engine_of(directory, "validate").predict()
    honest = clone(layered, tmp_path / "honest")
    (honest / "engine/states/general_layer.json").unlink()
    for stage in ("validate", "test"):
        with pytest.raises(E.EngineError, match="no registered finished run"):
            engine_of(honest, stage).layer_pins()
    assert not (honest / "engine/states/general_layer.json").exists()
    assert not (directory / "engine/validate/PREDICTION_COMMITMENT.json").exists()
    assert not (directory / "engine/states/general_layer.json").exists()            # nothing was registered by the refusal
    # round 2: development no longer registers it either (it did before: 'trust on first use'); the run that the build
    # step anchored next to the world's description is the only one any stage uses
    with pytest.raises(E.EngineError, match="is not the run anchored for this world"):
        engine_of(directory, "develop").layer_pins()
    assert not (directory / "engine/states/general_layer.json").exists()


# ------------------------------------------------------------------------------------------
# fallback: no inferred profile means E2's weights, bit for bit, and the row says why
# ------------------------------------------------------------------------------------------


def assert_is_e2(forecasts, reason):
    e7, e2 = forecasts[E7], forecasts[E2]
    assert E.same_arbor(e7["arbor"], e2["arbor"]) and E.arbor_sha(e7["arbor"]) == E.arbor_sha(e2["arbor"])
    assert e7["selected_donor"] == e2["selected_donor"] and e7["donor_ids"] == e2["donor_ids"]
    assert e7["fallback"] is True and f"no_inferred_profile:{reason}" in e7["fallback_reason"]
    assert e7["executed_policy"].endswith("+E2_weights") and e7["details"]["inferred_profile_used"] is False
    assert e7["general_layer"]["used"] is False and e7["general_layer"]["reason_code"] == reason
    assert e7["general_layer"]["adapter_status"] == BLOCKED


def test_fallback_when_no_finished_run_exists_e7_is_e2_and_says_general_layer_run_missing(plain, tmp_path):
    directory = clone(plain, tmp_path / "w")
    assert B.run_pins(E.synthetic_world(directory)) == {"run": None, "reason": "general_layer_run_missing"}
    for stage, part in (("develop", "train"), ("validate", "validation")):
        engine = engine_of(directory, stage)
        engine.prepare_prediction()
        for neuron in engine.queries()[:4]:
            forecasts, log = forecasts_of(engine, neuron)
            assert_is_e2(forecasts, "general_layer_run_missing")
            assert forecasts[E7]["general_layer"]["refused_by"] == "bridge" and forecasts[E7]["general_layer"]["state_sha256"] is None
            assert forecasts[E7]["general_layer"]["level"] is None
            assert log["general_layer"] == {"fold": None, "state_sha256": None, "state_animals_sha256": None, "used": False,
                                            "reason_code": "general_layer_run_missing", "query_animal_in_state": False}
    engine_of(directory, "validate").predict()
    manifest = load(directory / "engine/validate/prediction_manifest.json")
    rows = {(p["query_id"], p["model"]): p for p in manifest["predictions"]}
    assert manifest["pins"]["general_layer"] == {"run": None, "reason": "general_layer_run_missing"}
    assert manifest["general_layer"]["use_by_E7"]["fallback_by_reason"] == {"general_layer_run_missing": len(manifest["queries"])}
    assert manifest["general_layer"]["registration"] is None and manifest["general_layer"]["upstream_lists"] == []
    for q in manifest["queries"]:
        assert rows[(q, E7)]["arrays_sha256"] == rows[(q, E2)]["arrays_sha256"]
        assert rows[(q, E7)]["upstream_fit_animals_sha256"] is None                   # no fitted state touched this forecast
    # an E7 that is E2 can never pass the gate: it must beat E2 strictly on both measures
    summary = {m: {"distance_um": 900.0, "f1_200": 0.2} for m in E.METHODS}
    summary[E2] = summary[E7] = {"distance_um": 600.0, "f1_200": 0.5}
    gate = E.engine_gate(summary, "baseline", "all")
    assert gate["selected"] == E2 and not gate["candidates"][E7]["eligible"]
    assert gate["candidates"][E7]["beats_matched_control"] == {"E1_local_uniform": True, E2: False}


def test_fallback_when_the_layer_refuses_e7_is_e2_and_the_row_carries_the_layers_reason(layered):
    engine = engine_of(layered, "validate")
    engine.prepare_prediction()
    neuron = engine.cohort.part_neurons["validation"][0]
    row = engine.cohort.rows[neuron]
    clean, _ = forecasts_of(engine, neuron)
    assert clean[E7]["details"]["inferred_profile_used"] and clean[E7]["general_layer"]["used"]
    state = engine.layer_session().layer(None)[2]
    # (a) out of support: a soma far outside the range the state was fitted on. The layer refuses; nothing is extrapolated
    far = (row.soma[0], row.soma[1], row.soma[2] - 3000.0)
    engine.cohort.rows[neuron] = dataclasses.replace(row, soma=far)
    forecasts, log = forecasts_of(engine, neuron)
    assert_is_e2(forecasts, "out_of_support")
    assert forecasts[E7]["general_layer"]["refused_by"] == "layer" and "outside the training range" in forecasts[E7]["general_layer"]["reason"]
    assert forecasts[E7]["general_layer"]["level"] is None               # the layer issued nothing
    assert forecasts[E7]["general_layer"]["state_sha256"] == state.state_hash()              # the state was consulted: it is named
    assert log["general_layer"]["state_animals_sha256"] == E.ids_sha(sorted(state.groups))   # ... and so are its training animals
    engine.cohort.rows[neuron] = row
    # (b) the brain sample of the neuron is unknown to the layer: no state can be shown to exclude it
    session = engine.layer_session()
    real = session.infer
    session.infer = lambda row, **kw: real(row, group=core.NO_GROUP)
    forecasts, _ = forecasts_of(engine, neuron)
    assert_is_e2(forecasts, "forbidden_input")
    assert forecasts[E7]["general_layer"]["refused_by"] == "adapter" and "brain sample is unknown" in forecasts[E7]["general_layer"]["reason"]
    assert "group_unknown" in forecasts[E7]["general_layer"]["codes"]    # the layer answered with that warning; the adapter did not admit it
    session.infer = real
    # (c) a source label the state never saw is NOT a refusal under the layer's default policy: the value is issued from
    #     the soma alone and flagged; E7 uses what the layer issues and records the flag
    engine.cohort.rows[neuron] = dataclasses.replace(row, tag="LabelNeverSeen")
    forecasts, _ = forecasts_of(engine, neuron)
    layer = forecasts[E7]["general_layer"]
    assert layer["used"] and "out_of_support" in layer["codes"] and not layer["coverage_validated"]
    engine.cohort.rows[neuron] = row
    again, _ = forecasts_of(engine, neuron)
    assert E.arbor_sha(again[E7]["arbor"]) == E.arbor_sha(clean[E7]["arbor"]) and again[E7]["details"]["weights"] == clean[E7]["details"]["weights"]


# ------------------------------------------------------------------------------------------
# determinism
# ------------------------------------------------------------------------------------------


def test_two_processes_and_one_or_three_workers_give_identical_forecasts(layered, tmp_path):
    results = {}
    for name, workers in (("one", 1), ("three", 3)):
        directory = clone(layered, tmp_path / name)
        done = subprocess.run([sys.executable, "-I", str(ROOT / "analysis/axon_engine_cli.py"), "--world", str(directory),
                               "--workers", str(workers), "predict", "--stage", "validate"], capture_output=True, text=True, timeout=900)
        assert done.returncode == 0, done.stdout + done.stderr
        manifest = load(directory / "engine/validate/prediction_manifest.json")
        results[name] = {(p["query_id"], p["model"]): (p["arrays_sha256"], p["selected_donor"], json.dumps(p.get("general_layer"), sort_keys=True),
                                                       json.dumps(p["details"], sort_keys=True)) for p in manifest["predictions"]}
    assert results["one"] == results["three"] and len(results["one"]) == len({q for q, _ in results["one"]}) * len(E.METHODS)
    engine = engine_of(layered, "validate")                              # and this process gives the same arrays
    engine.prepare_prediction()
    mine = array_hashes(engine, engine.cohort.part_neurons["validation"])
    official = set(E.OFFICIAL_COPIES)
    assert all(results["one"][key][0] == sha for key, sha in mine.items() if key[1] not in official)
    assert any(json.loads(v[2])["used"] for (q, m), v in results["one"].items() if m == E7)


# ------------------------------------------------------------------------------------------
# binding: what the committed forecasts pin, and what makes admission refuse
# ------------------------------------------------------------------------------------------


def test_after_predict_a_changed_bridge_layer_file_run_or_state_makes_admit_refuse_before_any_truth(layered, tmp_path, monkeypatch):
    directory = clone(layered, tmp_path / "w")
    engine_of(directory, "validate").predict()
    base = directory / "engine/validate"
    out = layer_out(directory)
    manifest = load(base / "prediction_manifest.json")
    pinned = manifest["pins"]
    assert pinned["code_sha256"]["engine"][E.BRIDGE_FILE] == E.sha_file(ROOT / E.BRIDGE_FILE)
    assert pinned["code_sha256"]["general_layer"] == {f: E.sha_file(ROOT / f) for f in E.LAYER_CODE_FILES}
    run = pinned["general_layer"]["run"]
    assert run["run_complete_sha256"] == E.sha_file(out / "run_complete.json") and set(run["states"]) >= {"serving", "outer0"}
    assert run["states"]["serving"]["file_sha256"] == E.sha_file(out / run["states"]["serving"]["file"])
    assert "general_layer" in E.CONFIG and manifest["config"]["general_layer"]["bridge"] == E.BRIDGE_FILE
    refuses = lambda text: pytest.raises(E.EngineError, match=text)      # noqa: E731
    # (a) the bridge file changed after the forecasts were committed
    real_engine = E.engine_code_sha
    monkeypatch.setattr(E, "engine_code_sha", lambda root=E.ROOT: {**real_engine(root), E.BRIDGE_FILE: "0" * 64})
    with refuses(r"method changed since these forecasts were committed \(code_sha256"):
        engine_of(directory, "validate").admit()
    monkeypatch.undo()
    # (b) a code file of the general layer changed
    for relative in E.LAYER_CODE_FILES:
        real_layer = E.general_layer_code_sha
        monkeypatch.setattr(E, "general_layer_code_sha", lambda root=E.ROOT, r=relative, f=real_layer: {**f(root), r: "1" * 64})
        with refuses(r"method changed since these forecasts were committed \(code_sha256"):
            engine_of(directory, "validate").admit()
        monkeypatch.undo()
    # (c) a file of the finished run changed
    original = (out / "metrics.json").read_bytes()
    (out / "metrics.json").write_bytes(original + b" ")
    with refuses("method changed since these forecasts were committed"):
        engine_of(directory, "validate").admit()
    (out / "metrics.json").write_bytes(original)
    # (d) a state file changed: the all-training state that answered, and a fold state that this stage did not use
    for key in ("serving", "outer0"):
        path = out / run["states"][key]["file"]
        original = path.read_bytes()
        head, arrays = read_state(path)
        head["ranges"]["soma_position:AP"][1] += 1000.0
        write_state(path, head, arrays, rehash=True)
        with refuses("method changed since these forecasts were committed"):
            engine_of(directory, "validate").admit()
        path.write_bytes(original)
    # (e) the copy of the all-training record (serving.json) edited
    original = (out / "serving.json").read_bytes()
    copy = load(out / "serving.json")
    copy["tasks"]["T1"]["components"] += 1
    (out / "serving.json").write_text(core.canonical(copy))
    with refuses("method changed since these forecasts were committed"):
        engine_of(directory, "validate").admit()
    (out / "serving.json").write_bytes(original)
    assert not (base / "TRUTH_OPENED.json").exists() and not (base / "truth").exists()
    # everything restored: the same committed forecasts are admitted
    assert engine_of(directory, "validate").admit()["admitted"] > 0
    # and the freeze pins the same things
    engine_of(directory, "validate").score()
    record = engine_of(directory, "validate").freeze()
    assert record["pins"]["general_layer"]["run"]["run_complete_sha256"] == run["run_complete_sha256"]
    assert record["general_layer"]["finished_run_verified_at_freeze"] is True
    path = out / run["states"]["serving"]["file"]
    head, arrays = read_state(path)
    head["ranges"]["soma_position:AP"][1] += 1000.0
    write_state(path, head, arrays, rehash=True)
    with pytest.raises(E.FreezeMismatch, match="general_layer"):
        engine_of(directory, "test").verify_freeze()


def test_real_cohort_rule_a_used_run_must_be_named_in_the_protocol(layered, tmp_path):
    """Held-out stages of the REAL cohort: without a finished run the family 'all' is refused; with one, its
    completion record's SHA-256 must be written in the protocol. Shown on a synthetic world treated as the real
    cohort (no file of the real cohort is involved)."""
    directory = clone(layered, tmp_path / "w")
    protocol = tmp_path / "protocol.md"
    protocol.write_text("# a protocol that does not name the run\n")
    world = dataclasses.replace(E.synthetic_world(directory), protocol=str(protocol))
    as_real = dataclasses.replace(world, kind="real")

    class Probe(E.Engine):                                               # the bridge keeps reading the synthetic layer world
        def layer_run_on_disk(self):
            return E.Engine.layer_run_on_disk(E.Engine(world, self.stage, workers=1, log=None, family=self.family))

        def layer_pins(self, *, fresh=False):
            return E.Engine.layer_pins(E.Engine(world, self.stage, workers=1, log=None, family=self.family), fresh=fresh)

        def layer_session(self):
            return E.Engine.layer_session(E.Engine(world, self.stage, workers=1, log=None, family=self.family))

    probe = Probe(as_real, "validate", workers=1, log=None, family="all")
    sha = E.sha_file(layer_out(directory) / "run_complete.json")
    with pytest.raises(E.EngineError, match="is not named on the protocol's anchor line"):
        probe.require_general_layer()
    # review finding 10 and the attackers' suspicion: the full SHA-256 anywhere else in the text is not an acceptance
    protocol.write_text(f"# protocol\n\nThe run {sha} is NOT accepted.\n\n"
                        "**general layer run_complete.json SHA-256: (not written yet)**\n")
    with pytest.raises(E.EngineError, match="is not named on the protocol's anchor line"):
        probe.require_general_layer()
    protocol.write_text(f"# protocol\n\nThe run {sha} is NOT accepted (no anchor line at all).\n")
    with pytest.raises(E.EngineError, match="is not named on the protocol's anchor line"):
        probe.require_general_layer()
    protocol.write_text(f"**general layer run_complete.json SHA-256: {sha}**\n**general layer run_complete.json SHA-256: {'0' * 64}**\n")
    with pytest.raises(E.EngineError, match="more than one anchor line"):
        probe.require_general_layer()
    # round 2: only the exact whole line counts (the earlier spelling with a capital and back quotes is now refused)
    protocol.write_text(f"# protocol\n\nGeneral layer run_complete.json SHA-256: `{sha}`\n")
    with pytest.raises(E.EngineError, match="stands outside an anchor line"):
        probe.require_general_layer()
    protocol.write_text(f"# protocol\n\n**general layer run_complete.json SHA-256: {sha}**\n")
    probe.require_general_layer()                                        # named: accepted, and the serving state is verified
    E.Engine(world, "validate", workers=1, log=None).require_general_layer()               # a synthetic world: no such rule
    # no finished run at all
    bare = clone(layered, tmp_path / "bare")
    shutil.rmtree(B.layer_directory(E.synthetic_world(bare)))
    (bare / "engine/states/general_layer.json").unlink()
    bare_world = E.synthetic_world(bare)

    class Bare(E.Engine):
        def layer_run_on_disk(self):
            return E.Engine.layer_run_on_disk(E.Engine(bare_world, self.stage, workers=1, log=None, family=self.family))

        def layer_pins(self, *, fresh=False):
            return E.Engine.layer_pins(E.Engine(bare_world, self.stage, workers=1, log=None, family=self.family), fresh=fresh)

    with pytest.raises(E.EngineError, match="general_layer_run_missing"):
        Bare(dataclasses.replace(bare_world, kind="real"), "validate", workers=1, log=None, family="all").require_general_layer()
    Bare(dataclasses.replace(bare_world, kind="real"), "validate", workers=1, log=None, family="no_fitted_state").require_general_layer()


# ------------------------------------------------------------------------------------------
# the candidate itself, the family, the gate
# ------------------------------------------------------------------------------------------


def test_e7_is_the_ninth_candidate_conditional_with_a_fitted_state_and_last_in_the_gate_order():
    assert E.CANDIDATES[-1] == E7 and len(E.CANDIDATES) == 9 and len(E.METHODS) == 16
    assert E.FAMILIES["all"] == E.CANDIDATES and len(E.FAMILIES["all"]) == 9
    assert E.FAMILIES["no_fitted_state"] == ("E1_local_uniform", E2, "E3_local_kernel_label", "E6_local_consensus")
    assert E.MATCHED_CONTROLS[E7] == ("E1_local_uniform", E2) and E7 in E.USES_FITTED_STATE
    assert E.USES_LAYER_STATE == (E7,) and E7 not in E.USES_REGIONAL_STATE
    roles = E.Engine.referee_roles({"selected": E7, "descriptive_reference": None})
    assert roles["candidate"] == E7 and roles["matched_controls"] == ["C4_uniform_attach_500", "C5_uniform_attach_1000", "E1_local_uniform", E2]
    base = {m: {"distance_um": 900.0, "f1_200": 0.20} for m in E.METHODS}
    win = {**base, E7: {"distance_um": 600.0, "f1_200": 0.40}}
    assert E.engine_gate(win, "baseline", "all")["selected"] == E7
    assert E.engine_gate(win, "baseline", "no_fitted_state")["selected"] is None            # outside the family without a fitted state
    # it must beat E1 and E2 as well as every control
    for rival in ("E1_local_uniform", E2):
        gate = E.engine_gate({**win, rival: {"distance_um": 590.0, "f1_200": 0.41}}, "baseline", "all")
        assert not gate["candidates"][E7]["eligible"] and gate["selected"] == rival
    # exact ties with an earlier candidate go to the earlier one: E7 is last in the fixed order
    tie = {**base, E7: {"distance_um": 600.0, "f1_200": 0.40}, "E4_local_kernel_regional": {"distance_um": 600.0, "f1_200": 0.40}}
    assert E.engine_gate(tie, "baseline", "all")["selected"] == "E4_local_kernel_regional"
    assert "E7" in E.CONFIG["general_layer"]["candidate"] and E.CONFIG["engine_candidates"][-1] == E7
    assert "by inheritance" in E.CONFIG["inclusion"]["non_axon_connectors"]


# ------------------------------------------------------------------------------------------
# a world in which E7 is selected: the test goes through the independent referee, and the viewer gate
# ------------------------------------------------------------------------------------------


def test_e7_selected_bundle_is_frozen_sealed_and_scored_by_the_referee_with_honest_upstream_lists(e7_world, tmp_path):
    """The world in which the validation gate selects E7. The test goes through the referee's four functions with the
    training brain samples of the general layer's state declared at all three levels and E1 and E2 as matched
    controls. Whether E7 then meets the rule is read from the referee's result; it is not assumed."""
    directory, summary = e7_world
    out = Path(E.synthetic_world(directory).out)
    engine = engine_of(directory, "test")
    Ref = E.referee()
    assert summary["validate"]["selected"] == E7 and E7 in summary["validate"]["eligible"]
    freeze = load(out / "freeze/frozen.json")
    assert freeze["selected"] == E7 and freeze["family"] == "all" and freeze["referee"]["bundle_methods_with_a_fitted_state"] == [E7]
    assert freeze["general_layer"]["finished_run_verified_at_freeze"] and freeze["referee"]["can_seal_the_test_bundle"]
    assert summary["test"]["referee_functions_called"] == ["freeze_declaration", "copying_baselines", "seal_predictions", "score_once"]
    base = out / "test/referee"
    for name in ("declaration.json", "COHORT_FROZEN.json", "prediction_manifest.json", "PREDICTION_SEAL.json", "SCORE_STARTED.json",
                 "metrics.json", "SCORE_COMPLETE.json"):
        assert (base / name).exists(), name
    declaration, frozen = Ref.frozen_declaration(base)
    bundle, committed = load(base / "prediction_manifest.json"), load(out / "test/prediction_manifest.json")
    roles = declaration["models"]
    assert roles == {"candidate": E7, "copy_baselines": list(E.COPY_BASELINES),
                     "matched_controls": ["C4_uniform_attach_500", "C5_uniform_attach_1000", "E1_local_uniform", E2]}
    # the upstream list at all three levels: the training animals of the general layer's all-training state
    session, cohort = session_of(directory)
    state = session.layer(None)[2]
    animals = sorted(state.groups)
    assert animals == sorted(cohort.part_animals["train"]) and not set(animals) & cohort.heldout_animals
    assert declaration["upstream_fit_animals"] == bundle["upstream_fit_animals"] == animals
    by_key = {(p["query_id"], p["model"]): p for p in committed["predictions"]}
    admitted = load(out / "test/ADMISSION_COMPLETE.json")["admitted_ids"]
    models = [roles["candidate"], *roles["copy_baselines"], *roles["matched_controls"]]
    assert {(r["query_id"], r["model"]) for r in bundle["predictions"]} == {(q, m) for q in admitted for m in models}
    for r in bundle["predictions"]:
        mine = by_key[(r["query_id"], r["model"])]
        assert r["sha256"] == mine["sha256"] == Ref.sha(base / r["path"])
        assert r["hemisphere"] == cohort.rows[r["query_id"]].hemisphere in ("Left", "Right")
        if r["model"] == E7:
            assert set(r) == {"query_id", "model", "path", "sha256", "donor_ids", "upstream_fit_animals", "hemisphere", "general_layer"}
            assert r["upstream_fit_animals"] == animals and r["general_layer"] == mine["general_layer"]
            assert r["general_layer"]["state_sha256"] == state.state_hash() and r["general_layer"]["fold"] == "serving"
            assert len(r["donor_ids"]) == E.CONFIG["engine"]["k_neighbours"]
        else:
            assert set(r) == {"query_id", "model", "path", "sha256", "donor_ids", "upstream_fit_animals", "hemisphere"}
            assert r["upstream_fit_animals"] == []
    # what the method was, bound into the bundle: the bridge, the layer's code, the finished run and its T1 state
    artifacts = declaration["method_artifacts"]
    needed = [f"method/source/{E.BRIDGE_FILE}", *(f"method/source/{f}" for f in E.LAYER_CODE_FILES),
              "method/general_layer/run_complete.json", "method/general_layer/pre_run_declaration.json",
              "method/general_layer/evidence.json", "method/general_layer/selection.json",
              "method/general_layer/engine_registration.json", "method/regional_state/state_manifest.json"]
    for name in needed:
        assert artifacts[name] == Ref.sha(base / name), name
    assert any(name.startswith("method/general_layer/states/serving_") for name in artifacts)       # the all-training T1 state file
    bound = load(base / "method/general_layer/run_complete.json")[bench.SERVING_KEY]["tasks"]["T1"]  # its record, inside the completion record
    assert bound["state_hash"] == state.state_hash()
    assert artifacts["method/general_layer/run_complete.json"] == freeze["pins"]["general_layer"]["run"]["run_complete_sha256"]
    # scored once by the referee; the engine's table is the referee's
    result, mine = load(base / "metrics.json"), load(out / "test/metrics.json")
    assert result["status"] == "scored_once" and mine["scored_by"].startswith("independent referee")
    for m in models:
        for key in E.SCORE_KEYS:
            assert mine["summary"][m][key] == result["overall"]["summary"][m][key]
    verdict = mine["verdict"]
    assert verdict["selected_candidate"] == E7 and verdict["referee_candidate_slot"] == E7
    assert verdict["success"] == bool(result["overall"]["beats_all_copying"] and result["overall"]["beats_all_matched"])
    assert set(verdict["versus"]) == set(models) - {E7} and {"E1_local_uniform", E2} <= set(verdict["versus"])
    assert all(text.startswith("refused") for text in summary["once_only"].values())
    use = summary["general_layer_use_by_E7"]["test"]
    assert use["forecasts"] == len(admitted) == use["engine_level_value"] + use["baseline_level_value"] + use["fallback_to_E2_weights"]
    # a held-out animal in an E7 row's upstream list: the referee refuses the bundle
    forged = json.loads(json.dumps(bundle))
    next(r for r in forged["predictions"] if r["model"] == E7)["upstream_fit_animals"] = animals + [cohort.part_animals["test"][0]]
    with pytest.raises(ValueError, match="Upstream fit saw validation/test animals"):
        Ref.validate_bundle(base, forged, declaration)
    # the viewer export gate (read-only use), with a plan file that pins this synthetic manifest: a bundle is exported
    # only when the candidate beat every copying control AND its matched controls on both measures
    import prediction_viewer_export as V
    plan = tmp_path / "plan.md"
    cohort_manifest = directory / "data/manifest.tsv"
    plan.write_text(f"# scratch plan for a synthetic world\n\nmanifest SHA-256 `{V.sha(cohort_manifest)}`\n")
    target = tmp_path / "export.json"
    if verdict["success"]:
        export = V.export_paths(base, target, plan=plan, manifest=cohort_manifest)
        assert export["roles"]["candidate"] == E7 and export["benchmark_status"]["beats_all_matched"] is True
        assert all(n["hemisphere"] == cohort.rows[n["id"]].hemisphere in ("Left", "Right") for n in export["neurons"])
        layer = export["neurons"][0]["layers"]["prediction"]
        assert layer["model"] == E7 and layer["provenance"]["general_layer"]["adapter_status"] == BLOCKED
        assert layer["provenance_fields_omitted_as_large"] == [] and layer["upstream_fit_animals"]["count"] == len(animals)
    else:
        # what this world shows today: E7 is better than every copying control on the test animals and NOT better than
        # its matched control E2. The rule is not met, nothing is claimed, and nothing is exported.
        assert verdict["verdict"].startswith("rule not met") or verdict["verdict"] == "does not beat copying"
        with pytest.raises(V.ExportRefused, match="success condition not met"):
            V.export_paths(base, target, plan=plan, manifest=cohort_manifest)
        assert not target.exists()


def test_the_viewer_export_gate_accepts_a_passing_bundle_and_reads_the_hemisphere_from_its_rows(passing_world, tmp_path):
    directory, summary = passing_world
    assert summary["test"]["success"] is True and summary["test"]["selected"] == E2
    out = Path(E.synthetic_world(directory).out)
    base = out / "test/referee"
    cohort = E.Cohort(E.synthetic_world(directory))
    bundle = load(base / "prediction_manifest.json")
    assert all(r["hemisphere"] == cohort.rows[r["query_id"]].hemisphere in ("Left", "Right") for r in bundle["predictions"])
    assert all("general_layer" not in r and r["upstream_fit_animals"] == [] for r in bundle["predictions"])     # no fitted state in this bundle
    import prediction_viewer_export as V
    plan = tmp_path / "plan.md"
    cohort_manifest = directory / "data/manifest.tsv"
    plan.write_text(f"# scratch plan for a synthetic world\n\nmanifest SHA-256 `{V.sha(cohort_manifest)}`\n")
    export = V.export_paths(base, tmp_path / "export.json", plan=plan, manifest=cohort_manifest)
    assert export["roles"]["candidate"] == E2 and export["benchmark_status"]["beats_all_copying"] is True
    assert export["benchmark_status"]["beats_all_matched"] is True and (tmp_path / "export.json").exists()
    admitted = load(out / "test/ADMISSION_COMPLETE.json")["admitted_ids"]
    assert [n["id"] for n in export["neurons"]] == sorted(admitted)
    hemispheres = {n["id"]: n["hemisphere"] for n in export["neurons"]}
    assert "unknown" not in hemispheres.values() and set(hemispheres.values()) == {"Left", "Right"}
    assert all(hemispheres[n] == cohort.rows[n].hemisphere for n in admitted)
    # the three layers of a neuron and the comparison the viewer shows
    first = export["neurons"][0]
    assert set(first["layers"]) == {"measured", "prediction", "copy_baseline"} and first["layers"]["prediction"]["model"] == E2
    assert first["layers"]["prediction"]["provenance"]["hemisphere"] == first["hemisphere"]
    # E7 and the general layer in this world: forecast for every query, outside the declared family, never selectable
    assert summary["freeze"]["family"] == "no_fitted_state" and E7 not in E.FAMILIES["no_fitted_state"]
    use = summary["general_layer_use_by_E7"]["test"]
    assert use["forecasts"] == len(admitted)


def test_finding9_the_layer_command_says_that_it_registers_the_run_on_first_sight():
    """Review finding 9: `layer` writes the world's registration of the finished run (states/general_layer.json) when it
    sees a run for the first time; its help text and the protocol say so instead of calling it read-only."""
    text = (ROOT / "analysis/axon_engine_cli.py").read_text()
    assert "registers the finished run on first sight" in text
    assert "registers the finished run on first sight" in (ROOT / "reports/axon_engine_protocol.md").read_text()


# ------------------------------------------------------------------------------------------
# second attack round (round 2): the protocol's anchor line, and the synthetic anchor of the finished run
# ------------------------------------------------------------------------------------------


def test_round2_only_a_whole_visible_anchor_line_in_the_protocol_counts(layered, tmp_path):
    """Round-2 finding: the anchor phrase was matched anywhere, in any letter case, so a hidden HTML comment, a code
    block or a sentence that REJECTS the run was read as the acceptance. Now only one whole line in exactly the stated
    spelling counts, and the phrase anywhere else makes the engine refuse the protocol."""
    world = E.synthetic_world(layered)
    sha = E.sha_file(layer_out(layered) / "run_complete.json")
    protocol = tmp_path / "protocol.md"

    def anchor_of(text):
        protocol.write_text(text)
        return E.Engine(dataclasses.replace(world, protocol=str(protocol)), "develop", workers=1, log=None).protocol_anchor()

    outside = "stands outside an anchor line"
    for text in (f"# p\n\nThe lead does NOT accept this run: general layer run_complete.json SHA-256: {sha}\n",
                 f"# p\n\nThe line **general layer run_complete.json SHA-256: {sha}** is NOT accepted.\n",
                 f"# p\n\n<!-- general layer run_complete.json SHA-256: {sha} -->\n",
                 f"# p\n\n<!--\n**general layer run_complete.json SHA-256: {sha}**\n-->\n",
                 f"# p\n\n```\n**general layer run_complete.json SHA-256: {sha}**\n```\n",
                 f"# p\n\n~~~\ngeneral layer run_complete.json sha-256: {sha}\n~~~\n",
                 f"# p\n\n**General layer run_complete.json SHA-256: {sha}**\n",
                 f"# p\n\n**general layer run_complete.json SHA-256: {sha}** (accepted)\n",
                 f"# p\n\n**general layer run_complete.json SHA-256: {sha.upper()}**\n"):
        with pytest.raises(E.EngineError, match=outside):
            anchor_of(text)
    with pytest.raises(E.EngineError, match="more than one anchor line"):
        anchor_of(f"**general layer run_complete.json SHA-256: {sha}**\n**general layer run_complete.json SHA-256: (not written yet)**\n")
    assert anchor_of("# p\n\n**general layer run_complete.json SHA-256: (not written yet)**\n") is None
    assert anchor_of(f"# p\n\nThe run {sha} is mentioned, not anchored.\n") is None
    assert anchor_of(f"# p\n\n<!-- a comment -->\n\n  **general layer run_complete.json SHA-256: {sha}**  \n") == sha
    # the accepted protocol as it is now: the placeholder, read without error
    assert E.Engine(E.real_world(legacy_only=True), "develop", workers=1, log=None).protocol_anchor() is None


def test_round2_an_edited_or_remade_registration_cannot_bring_a_forged_run_in(layered, tmp_path):
    """Round-2 finding: on a synthetic world the registration (engine/states/general_layer.json) was the only anchor,
    so editing its SHA-256 to a forged run, or deleting it and running development, let a forged run's ENGINE-level
    values into committed validation forecasts. The step that builds the run now writes an anchor next to the world's
    own description (data/general_layer_anchor.json); every stage refuses a run that the anchor does not name."""
    ANCHOR = "general_layer_anchor.json"
    genuine = E.sha_file(layer_out(layered) / "run_complete.json")
    # (a) the registration's SHA-256 edited to the forged run
    a = clone(layered, tmp_path / "a")
    forged = forge_run(a)
    registration = load(a / "engine/states/general_layer.json")
    registration["run_complete_sha256"] = forged
    E.write_json(a / "engine/states/general_layer.json", registration, pretty=True)
    for stage in ("develop", "validate", "test"):
        with pytest.raises(E.EngineError, match="is not the run anchored for this world"):
            engine_of(a, stage).layer_pins()
    with pytest.raises(E.EngineError, match="is not the run anchored for this world"):
        engine_of(a, "validate").predict()
    assert not (a / "engine/validate/PREDICTION_COMMITMENT.json").exists()
    # (b) the registration deleted, then development run again: nothing is registered, every stage refuses
    b = clone(layered, tmp_path / "b")
    forge_run(b)
    (b / "engine/states/general_layer.json").unlink()
    with pytest.raises(E.EngineError, match="is not the run anchored for this world"):
        engine_of(b, "develop").layer_pins()
    with pytest.raises(E.EngineError, match="is not the run anchored for this world"):
        engine_of(b, "develop").predict()
    assert not (b / "engine/states/general_layer.json").exists()
    with pytest.raises(E.EngineError, match="no registered finished run|is not the run anchored"):
        engine_of(b, "validate").predict()
    assert not (b / "engine/validate/PREDICTION_COMMITMENT.json").exists()
    # (c) the anchor deleted as well: a run that the build step did not anchor is never used, in development either
    c = clone(layered, tmp_path / "c")
    forge_run(c)
    (c / "engine/states/general_layer.json").unlink()
    (c / "data" / ANCHOR).unlink()
    with pytest.raises(E.EngineError, match="is not anchored for this world"):
        engine_of(c, "develop").layer_pins()
    # (d) ... and the build step does not anchor a run it finds on disk instead of building it
    with pytest.raises(E.EngineError, match="is not anchored for this world"):
        engine_of(c, "develop").prepare_general_layer()
    assert not (c / "data" / ANCHOR).exists() and not (c / "engine/states/general_layer.json").exists()
    # (e) the genuine run: the anchor is pinned with the forecasts, and an anchor changed afterwards makes admit refuse
    assert E.LAYER_ANCHOR_FILE == ANCHOR and load(layered / "data" / ANCHOR)["runs"]["engine"]["run_complete_sha256"] == genuine
    e = clone(layered, tmp_path / "e")
    engine_of(e, "validate").predict()
    pins = load(e / "engine/validate/prediction_manifest.json")["pins"]["general_layer"]
    assert pins["anchor"]["run_complete_sha256"] == pins["run"]["run_complete_sha256"] == genuine
    anchor = load(e / "data" / ANCHOR)
    anchor["runs"]["engine"]["run_complete_sha256"] = "0" * 64
    E.write_json(e / "data" / ANCHOR, anchor, pretty=True)
    with pytest.raises(E.EngineError, match="is not the run anchored for this world"):
        engine_of(e, "validate").admit()
    assert not (e / "engine/validate/TRUTH_OPENED.json").exists()
