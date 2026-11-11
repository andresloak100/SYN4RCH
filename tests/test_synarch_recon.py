"""SYNARCH reconstruction layer: synthetic checks, plus checks of the real training partition.

Everything statistical here runs on SYNTHETIC data. The few tests that touch the real cohort read the
training brain samples only, and compute nothing that relates an input to a regional value: they check
identifiers, counts, statuses and file handling. No test can start the real benchmark.
"""
import csv
import dataclasses
import importlib.util
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / "analysis" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


core = _load("synarch_recon")
mods = _load("synarch_recon_modules")
bench = _load("synarch_recon_benchmark")
REGISTRY = core.Registry.from_contract()

REGIONS = ("ARH", "VMH", "LHA", "PAG", "MY", "BST")
INSIDE = ("ARH", "VMH", "LHA")
LABELS = ("A", "B", "C")
CENTRES = np.array([[5800., 5600., 4700.], [6600., 6300., 5100.], [7300., 5500., 4500.], [6900., 6600., 5400.]])
LEVELS = np.array([[3.5, 1.2, 2.5, 1.0, 1.5, 3.0], [1.2, 3.4, 1.4, 2.8, 3.2, 1.1],
                   [2.4, 2.2, 3.6, 3.3, 1.1, 2.0], [3.1, 3.0, 1.2, 1.4, 2.6, 3.4]])
TABLE = np.array([[.8, .1, .1], [.1, .8, .1], [.1, .1, .8], [.45, .45, .1]])
PRESENT = np.array([[.95, .05, .9, .05, .1, .9], [.05, .95, .1, .9, .9, .05],
                    [.9, .9, .95, .9, .05, .5], [.9, .9, .05, .1, .9, .95]])
SD_SOMA, SD_LEVEL = 150.0, 0.3
POINTER = {"metrics": "synthetic test evidence", "metrics_sha256": "0" * 64}


# ----------------------------------------------------------------------------------------------
# fixtures
# ----------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def registry():
    return REGISTRY


@pytest.fixture(scope="module")
def layout(registry):
    return mods.jiao_layout(registry, REGIONS)


@pytest.fixture(scope="module")
def tasks(registry, layout):
    return mods.build_tasks(layout, registry, INSIDE)


def planted_rows(n, seed, prefix="u", groups=40):
    """An exact latent-class world: class -> soma (Gaussian), label (table), log-lengths (Gaussian).
    A MIRRORED cohort, like the published table: a soma beyond the midline is mirrored (ML' = 11,400 - ML), so
    every ML is <= 5,700 (about 2% of the neurons of the most medial class are affected)."""
    rng = np.random.default_rng(seed)
    cls = rng.integers(0, 4, n)
    soma = CENTRES[cls] + rng.normal(0, SD_SOMA, (n, 3))
    soma[:, 2] = core.MIDLINE_ML - np.abs(core.MIDLINE_ML - soma[:, 2])
    y = LEVELS[cls] + rng.normal(0, SD_LEVEL, (n, 6))
    label = np.array([LABELS[rng.choice(3, p=TABLE[c])] for c in cls], dtype=object)
    cont = {f"soma_position:{a}": soma[:, j] for j, a in enumerate(core.SOMA_AXES)}
    cont.update({f"regional_axon_length:{r}": y[:, j] for j, r in enumerate(REGIONS)})
    group = np.array([f"{prefix}g{i % groups:03d}" for i in range(n)])
    unit = np.array([f"{g}_{i:05d}" for i, g in enumerate(group)])      # the Jiao convention: <brain sample>_<n>
    rows = mods.Rows(unit, group, cont, {"source_label": label})
    return rows, soma, y, label


def zero_inflated_rows(n, seed, prefix="z", groups=40):
    """The same latent classes, but a region is reached only with a class-specific probability: exact zeros."""
    rows, soma, y, label = planted_rows(n, seed, prefix, groups)
    rng = np.random.default_rng(seed + 10_000)
    cls = nearest_class_of_rows(n, seed)
    reached = rng.random((n, 6)) < PRESENT[cls]
    y = np.where(reached, np.maximum(y, 0.05), 0.0)
    for j, r in enumerate(REGIONS):
        rows.cont[f"regional_axon_length:{r}"] = y[:, j]
    return rows, soma, y, label


def nearest_class_of_rows(n, seed):
    return np.random.default_rng(seed).integers(0, 4, n)       # the first draw of planted_rows: the true classes


def class_posterior(log_terms):
    log_terms = log_terms - log_terms.max(axis=1, keepdims=True)
    p = np.exp(log_terms)
    return p / p.sum(axis=1, keepdims=True)


def true_profile_given_soma_label(soma, label):
    index = np.array([LABELS.index(l) for l in label])
    log_terms = np.log(TABLE[:, index]).T - 0.5 * (((soma[:, None, :] - CENTRES[None]) / SD_SOMA) ** 2).sum(axis=2)
    return class_posterior(log_terms) @ LEVELS


def nearest_class(soma):
    return np.argmin(((soma[:, None, :] - CENTRES[None]) ** 2).sum(axis=2), axis=1)


def true_soma_given_profile_label(y, label):
    index = np.array([LABELS.index(l) for l in label])
    log_terms = np.log(TABLE[:, index]).T - 0.5 * (((y[:, None, :] - LEVELS[None]) / SD_LEVEL) ** 2).sum(axis=2)
    return class_posterior(log_terms) @ CENTRES


@pytest.fixture(scope="module")
def planted(layout):
    rows, soma, y, label = planted_rows(2400, 11)
    state, diagnostics = mods.fit_joint(rows, layout, components=4, seed=5,
                                        budget={"n_init": 5, "max_iter": 200, "tol": 1e-7})
    return {"rows": rows, "state": state, "diagnostics": diagnostics}


def evidence_for(tasks, status, *, coverage=True, best="IDENTITY"):
    out = {}
    for name, task in tasks.items():
        out[task.conditional().name] = core.ValidationEvidence(
            status=status, task=name, benchmark=POINTER, primary_metric=task.primary, coverage_validated=coverage,
            best_comparator=best if task.kind != "label" else "MEAN")
    return out


def reconstructor(registry, layout, tasks, planted, status="exploratory_supported", policy=None, coverage=True,
                  ledger=None):
    module = mods.JointMixtureModule(registry, layout, states={n: planted["state"] for n in tasks},
                                     train=planted["rows"], tasks=tasks,
                                     evidence=evidence_for(tasks, status, coverage=coverage))
    r = core.Reconstructor(registry, policy, ledger)
    r.register(module)
    r.register(mods.AxonArborAdapter(registry))
    return r, module


def request(registry, observed, requested, unit="fresh_0001", group="fresh"):
    """A request about a unit of the brain sample 'fresh' (Jiao IDs are <brain sample>_<n>)."""
    return core.ReconstructionRequest.from_values(registry, namespace="jiao", unit_id=unit, group_id=group,
                                                  observed=observed, requested=requested)


def measured(variable, value, coordinates=(), unit="fresh_0001", group="fresh", **changes):
    """A measured Jiao record in the registered scale and frame; `changes` override single fields."""
    spec, ns = REGISTRY.variable("jiao", variable), REGISTRY.namespace("jiao")
    fields = dict(species=ns.species, namespace="jiao", grain="neuron", unit_id=unit,
                  group=core.Group(group, "brain_sample", ns.group_evidence), variable=variable, status="measured",
                  value=value, scale=core.scale_statement(spec), coordinates=tuple(coordinates), frame=spec.frame,
                  transform="identity")
    fields.update(changes)
    return core.Observation(**fields)


def query_rows(soma_row, label, unit="q", group="fresh"):
    return mods.Rows(np.array([unit]), np.array([group]),
                     {f"soma_position:{a}": np.array([soma_row[j]]) for j, a in enumerate(core.SOMA_AXES)},
                     {"source_label": np.array([label], dtype=object)})


def lengths(y_row, names=REGIONS):
    """Log-scale values -> published-scale lengths, as a request supplies them."""
    return {r: float(10 ** v - 1) for r, v in zip(names, y_row)}


@pytest.fixture(scope="module")
def null_world(tmp_path_factory):
    return bench.run_synthetic(tmp_path_factory.mktemp("null") / "world", "null", size="tiny", seed=0,
                               log=lambda *a: None)


@pytest.fixture(scope="module")
def planted_world(tmp_path_factory):
    return bench.run_synthetic(tmp_path_factory.mktemp("planted") / "world", "planted", size="tiny", seed=0,
                               log=lambda *a: None)


@pytest.fixture()
def bare_world(tmp_path):
    return bench.make_synthetic_world(tmp_path / "world", "planted", size="tiny", seed=3)


# ----------------------------------------------------------------------------------------------
# the unchanged joint model is verified at import
# ----------------------------------------------------------------------------------------------


def test_joint_core_hash_is_verified_at_import(tmp_path):
    info = core.JOINT_INFO["synarch_recon_brm_joint"]
    assert info["sha256"] == core.BRM_JOINT_SHA256 == core.sha_file(ROOT / "analysis/brm_joint.py")
    assert hasattr(core.JOINT, "JointModel")
    tampered = tmp_path / "brm_joint.py"
    tampered.write_bytes((ROOT / "analysis/brm_joint.py").read_bytes() + b"\n# changed\n")
    with pytest.raises(core.FrozenSourceChanged):
        core.load_joint_core(tampered, module_name="tampered_joint")
    wrong_pin = tmp_path / "frozen.json"
    wrong_pin.write_text(json.dumps({"source_sha256": {"analysis/brm_joint.py": "0" * 64}}))
    with pytest.raises(core.FrozenSourceChanged):
        core.load_joint_core(pin_file=wrong_pin, module_name="wrong_pin_joint")


# ----------------------------------------------------------------------------------------------
# registry seeded from the data contract
# ----------------------------------------------------------------------------------------------


def test_registry_holds_the_contract_namespaces_grains_and_joins(registry):
    grains = {name: spec.grain for name, spec in registry.namespaces.items()}
    assert grains["jiao"] == "neuron" and grains["epiretro"] == "nucleus"
    assert grains["ABC_RNA"] == "cell_type_population" and grains["curated_literature"] == "study_record"
    assert registry.namespace("jiao").group_kind == "brain_sample"
    assert "not asserted to be a distinct animal" in registry.namespace("jiao").group_evidence
    assert registry.namespace("epiretro").group_kind == "experiment"
    kinds = lambda a, b: {r.kind for r in registry.joins_between(a, b)}
    assert kinds("jiao", "jiao") == {"paired_same_unit"}
    assert kinds("epiretro", "ABC_RNA") == {"reference_population"}
    assert kinds("jiao", "HypoMap") == {"hypothesis_only", "forbidden"}
    assert registry.reference_rule("ABC_RNA", "epiretro") is not None
    assert registry.reference_rule("ABC_RNA", "jiao") is None and registry.reference_rule("HypoMap", "jiao") is None
    assert registry.contract["sha256"] == core.sha_file(ROOT / "data/interim/brm_data_contract.json")


def test_variable_registry_declares_kind_units_frame_and_transform(registry):
    profile = registry.variable("jiao", "regional_axon_length")
    assert (profile.kind, profile.transform, len(profile.components)) == ("vector", "log10_1p", 49)
    assert profile.family == "axon" and profile.measured_input_needs_declared_task
    assert "not resolved" in profile.units
    soma = registry.variable("jiao", "soma_position")
    assert soma.components == ("AP", "DV", "ML") and "CCFv3" in soma.frame and soma.transform == "identity"
    assert registry.variable("jiao", "soma_label_fine").derived_from == ("soma_position",)
    for banned in ("projectome_cluster", "projectome_type", "projectome_class"):
        assert registry.variable("jiao", banned).predictor_ban
        with pytest.raises(core.ForbiddenInputError):
            mods.FieldLayout.from_registry(registry, "jiao", ("soma_position", banned))
    assert registry.namespaces_measuring("gene_expression") == ["ABC_RNA", "HypoMap", "Rossi_2021"]
    np.testing.assert_allclose(core.inverse_transform("log10_1p", core.forward_transform("log10_1p", [0, 9, 99])),
                               [0, 9, 99])


# ----------------------------------------------------------------------------------------------
# pairing rule
# ----------------------------------------------------------------------------------------------


def _record(namespace, grain, unit, variable, value, coordinates=(), group=None):
    """A measured record that states the registered scale and frame of its variable."""
    spec, ns = REGISTRY.variable(namespace, variable), REGISTRY.namespace(namespace)
    if group is None:
        group = core.Group(REGISTRY.implied_group(namespace, unit) or "g1", ns.group_kind, "test")
    return core.Observation("Mus musculus", namespace, grain, unit, group, variable, "measured", value,
                            core.scale_statement(spec), coordinates, spec.frame, transform="identity")


def test_records_pair_only_with_identical_namespace_grain_and_unit(registry):
    a = _record("jiao", "neuron", "s1_001", "soma_position", (1.0, 2.0, 3.0), core.SOMA_AXES)
    same = _record("jiao", "neuron", "s1_001", "source_label", "Pomc")
    other_unit = _record("jiao", "neuron", "s1_002", "source_label", "Pomc")
    other_namespace = _record("epiretro", "nucleus", "s1_001", "retrograde_target", "VTA")
    assert core.paired(a, same)
    assert not core.paired(a, other_unit) and not core.paired(a, other_namespace)
    store = core.ObservationStore(registry)
    for record in (a, same, other_unit, other_namespace):
        store.add(record)
    assert set(store.join([a, same])) == {"soma_position", "source_label"}
    with pytest.raises(core.PairingError):
        store.join([a, other_unit])
    with pytest.raises(core.PairingError):
        store.join([a, other_namespace])          # same ID string, different dataset: never the same unit
    assert set(store.unit("jiao", "neuron", "s1_001")) == {"soma_position", "source_label"}
    assert set(store.unit("epiretro", "nucleus", "s1_001")) == {"retrograde_target"}
    with pytest.raises(core.PairingError):          # a jiao record cannot sit at another grain
        store.add(_record("jiao", "nucleus", "s1_003", "source_label", "Pomc"))
    # one unit belongs to one group: a second group for a stored unit is refused, by the store and by join()
    elsewhere = core.Group("s2", "brain_sample", "test")
    with pytest.raises(core.PairingError, match="by its ID"):          # the Jiao neuron ID names its brain sample
        store.add(_record("jiao", "neuron", "s1_001", "hemisphere", "Left", group=elsewhere))
    with pytest.raises(core.PairingError, match="one unit belongs to one group"):
        store.add(_record("epiretro", "nucleus", "s1_001", "assigned_subclass", "x",
                          group=core.Group("another_experiment", "experiment", "test")))
    with pytest.raises(core.PairingError, match="grouped by"):         # an experiment is not a brain sample
        store.add(_record("epiretro", "nucleus", "n9", "assigned_subclass", "x", group=elsewhere))
    with pytest.raises(core.PairingError, match="one unit belongs to one group"):
        store.join([a, dataclasses.replace(same, group=elsewhere)])
    with pytest.raises(ValueError, match="blanks"):                    # identifiers are compared as given: no blanks
        store.add(_record("epiretro", "nucleus", "n10 ", "assigned_subclass", "x"))


def test_reference_alignment_is_population_level_and_flagged_reference(registry):
    store = core.ObservationStore(registry)
    source = core.SourceRef("synthetic atlas", None, "subclass 001")
    record = store.attach_reference(from_namespace="ABC_RNA", to_namespace="epiretro", key_name="assigned_subclass",
                                    key_value="001 CLA-EPd-CTX Car3 Glut", variable="gene_expression",
                                    value=(1.5, 0.2), coordinates=("GeneA", "GeneB"),
                                    species="Mus musculus", source=source)
    assert record.status == "reference" and record.grain == "cell_type_population"
    assert record.scale == "cluster mean; log2(CPM+1)" and record.transform == "identity"    # the registered statement
    with pytest.raises(core.ScaleError):            # a reference value in another stated scale is refused
        store.attach_reference(from_namespace="ABC_RNA", to_namespace="epiretro", key_name="assigned_subclass",
                               key_value="002 other", variable="gene_expression", value=(1.0,), coordinates=("GeneA",),
                               scale="raw counts", species="Mus musculus", source=source)
    assert record.reference.from_namespace == "ABC_RNA" and "group-level" in record.reference.contract_allowed
    # never 'measured', never attached to an individual unit
    with pytest.raises(ValueError):
        dataclasses.replace(record, status="measured")
    with pytest.raises(ValueError):
        dataclasses.replace(record, grain="nucleus")
    # no reference alignment exists towards Jiao neurons: genes of another dataset cannot be attached
    for origin in ("ABC_RNA", "HypoMap"):
        with pytest.raises(core.PairingError):
            store.attach_reference(from_namespace=origin, to_namespace="jiao", key_name="source_label",
                                   key_value="Pomc", variable="gene_expression", value=(1.0,), coordinates=("GeneA",),
                                   species="Mus musculus", source=source)
    # the nuclei of the group see it only as a reference of their group
    module = mods.ReferenceAtlasModule(registry, store, namespace="epiretro", from_namespace="ABC_RNA",
                                       key_variable="assigned_subclass", variable="gene_expression")
    r = core.Reconstructor(registry)
    r.register(module)
    ask = core.ReconstructionRequest.from_values(
        registry, namespace="epiretro", unit_id="nucleus_17", group_id="exp_3",
        observed={"assigned_subclass": "001 CLA-EPd-CTX Car3 Glut"}, requested=["gene_expression"])
    answer = r.reconstruct(ask)["gene_expression"]
    assert answer.status == "reference" and "population_reference" in answer.codes()
    assert answer.estimate.point == (1.5, 0.2) and answer.reference.from_namespace == "ABC_RNA"
    with pytest.raises(core.ReconError):
        answer.to_observation(ask)                # a reference value never becomes an inferred unit record
    missing_group = core.ReconstructionRequest.from_values(
        registry, namespace="epiretro", unit_id="nucleus_18", group_id="exp_3",
        observed={"assigned_subclass": "999 unknown"}, requested=["gene_expression"])
    assert r.reconstruct(missing_group)["gene_expression"].refusal.code == "unidentifiable_cross_modal"
    # a reference value attaches only to the unit's OWN group: the record of subclass 001 handed in as evidence
    # for a nucleus that is assigned to another subclass, or to none, is refused
    for own in ({"assigned_subclass": "002 another subclass"}, {}):
        foreign = core.ReconstructionRequest.from_values(
            registry, namespace="epiretro", unit_id="nucleus_19", group_id="exp_3",
            observed=dict(own, gene_expression=record), requested=["retrograde_target"])
        refused = r.reconstruct(foreign)["retrograde_target"]
        assert refused.refusal.code == "unidentifiable_cross_modal" and "own group" in refused.refusal.reason
    same_group = core.ReconstructionRequest.from_values(
        registry, namespace="epiretro", unit_id="nucleus_20", group_id="exp_3",
        observed={"assigned_subclass": "001 CLA-EPd-CTX Car3 Glut", "gene_expression": record},
        requested=["retrograde_target"])
    assert r.reconstruct(same_group)["retrograde_target"].refusal.code == "unsupported_variable"   # accepted as evidence
    with pytest.raises(core.PairingError):
        mods.ReferenceAtlasModule(registry, store, namespace="jiao", from_namespace="ABC_RNA",
                                  key_variable="source_label", variable="gene_expression")


def test_outcome_derived_columns_cannot_be_stored_or_used(registry):
    store = core.ObservationStore(registry)
    with pytest.raises(core.ForbiddenInputError):
        store.add(_record("jiao", "neuron", "s1_001", "projectome_cluster", "c7"))


# ----------------------------------------------------------------------------------------------
# status and uncertainty
# ----------------------------------------------------------------------------------------------


def test_record_status_and_uncertainty_invariants():
    group = core.Group("s1", "brain_sample", "test")
    with pytest.raises(ValueError):                 # a missing value is not a zero
        core.Observation("m", "jiao", "neuron", "u", group, "soma_position", "missing", (0.0, 0.0, 0.0), "um",
                         core.SOMA_AXES)
    core.Observation("m", "jiao", "neuron", "u", group, "soma_position", "missing", None, "um")
    with pytest.raises(ValueError):                 # inferred needs module, state, training groups, coverage flag
        core.Observation("m", "jiao", "neuron", "u", group, "source_label", "inferred", "A", "nominal",
                         transform="identity")
    with pytest.raises(ValueError, match="states the scale"):          # a value without its scale tag is no record
        core.Observation("m", "jiao", "neuron", "u", group, "source_label", "measured", "A", "nominal; label")
    with pytest.raises(ValueError, match="derived from"):              # only inferred records have evidence
        measured("source_label", "A", derived_from=(("soma_position", "measured"),))
    with pytest.raises(ValueError, match="'inferred', not 'measured'"):
        measured("source_label", "A", source=core.SourceRef("reconstruct:joint_mixture", "s", "u"))
    with pytest.raises(ValueError):
        core.Uncertainty("inferred", producing_module="m", state_hash="s", training_groups_hash="g")
    inferred = core.Uncertainty("inferred", producing_module="m", state_hash="s", training_groups_hash="g",
                                coverage="unvalidated")
    with pytest.raises(ValueError):                 # a measured record cannot carry a model's uncertainty
        core.Observation("m", "jiao", "neuron", "u", group, "source_label", "measured", "A", "nominal",
                         uncertainty=inferred, transform="identity")
    with pytest.raises(ValueError):
        core.Uncertainty("sd", sd=1.0, producing_module="m")
    with pytest.raises(ValueError):
        core.ValidationEvidence(status="exploratory_supported", benchmark=POINTER, confirmatory=True)


def test_inferred_value_carries_module_state_training_hash_and_coverage(registry, layout, tasks, planted):
    r, module = reconstructor(registry, layout, tasks, planted)
    fresh, soma, y, label = planted_rows(5, 99, prefix="f")
    ask = request(registry, {"soma_position": soma[0], "source_label": label[0]}, ["regional_axon_length"])
    answer = r.reconstruct(ask)["regional_axon_length"]
    state = planted["state"]
    assert answer.status == "inferred" and answer.refusal is None and not answer.baseline_level
    assert answer.module == "joint_mixture" and answer.state_hash == state.state_hash()
    assert answer.training_groups_hash == state.groups_sha256 == core.ids_sha(state.groups)
    assert answer.coverage_validated and answer.validation == "exploratory_supported"
    assert answer.benchmark == POINTER and answer.conditional == "T1_profile_from_soma_label"
    assert dict(answer.evidence_used) == {"soma_position": "measured", "source_label": "measured"}
    e = answer.estimate
    assert e.coordinates == REGIONS and e.scale == "log10(1 + value)" and e.interval["level"] == 0.90
    assert all(lo <= p <= hi for lo, p, hi in zip(e.interval["lower"], e.point, e.interval["upper"]))
    np.testing.assert_allclose(e.original_units["median"], 10 ** np.array(e.median) - 1)
    record = answer.to_observation(ask)
    assert record.status == "inferred" and record.uncertainty.kind == "inferred"
    assert record.uncertainty.producing_module == "joint_mixture" and record.uncertainty.coverage == "validated"
    assert record.uncertainty.state_hash == state.state_hash()
    assert record.uncertainty.training_groups_hash == state.groups_sha256
    assert record.derived_from == answer.evidence_used and record.unit_id == "fresh_0001"
    assert record.group.group_id == "fresh" and answer.unit == ("jiao", "neuron", "fresh_0001")
    # scale statement, tag and frame travel with the value; the lineage names the unit; the layer keeps an entry
    profile_spec = registry.variable("jiao", "regional_axon_length")
    assert (record.scale, record.transform, record.frame) == ("log10(1 + value)", "log10_1p", profile_spec.frame)
    assert (e.scale, e.transform, e.frame) == (record.scale, record.transform, record.frame)
    assert record.lineage[0].unit == ("jiao", "neuron", "fresh_0001") and record.lineage[0].group_id == "fresh"
    fact = r.ledger.fact(record)
    assert fact.state_hash == state.state_hash() and fact.state_groups == frozenset(state.groups)
    assert fact.evidence == (("soma_position", "measured", None), ("source_label", "measured", None))
    with pytest.raises(ValueError):                 # and it can never be relabelled as measured
        dataclasses.replace(record, status="measured")
    with pytest.raises(core.PairingError):          # nor become a record of another unit, or of another group
        answer.to_observation(request(registry, {}, [], unit="fresh_0002"))
    with pytest.raises(core.PairingError):
        answer.to_observation(request(registry, {}, [], unit="fresh_0001", group="elsewhere"))
    # the single request equals the batch conditional of the same state
    batch = mods.engine_predict(state, tasks["T1"], fresh.take([0]), level=0.90, median=True)
    np.testing.assert_allclose(e.point, batch["point"][0], rtol=1e-10)
    np.testing.assert_allclose(e.interval["upper"], batch["upper"][0], rtol=1e-10)


def test_inferred_input_is_flagged_and_coverage_is_not_claimed(registry, layout, tasks, planted):
    lenient, _ = reconstructor(registry, layout, tasks, planted, policy=core.Policy(on_not_validated="engine_flagged"))
    _, soma, y, label = planted_rows(3, 98, prefix="f")
    first = lenient.reconstruct(request(registry, {"source_label": label[0]}, ["soma_position"]))
    soma_record = first["soma_position"].to_observation(first.request)          # inferred from the measured label
    assert soma_record.status == "inferred" and soma_record.uncertainty.coverage == "unvalidated"
    chained = lenient.reconstruct(request(registry, {"soma_position": soma_record},
                                          ["regional_axon_length"]))["regional_axon_length"]
    assert chained.status == "inferred" and "conditioned_on_inferred_input" in chained.codes()
    assert not chained.coverage_validated and chained.conditional is None
    assert dict(chained.evidence_used) == {"soma_position": "inferred"}
    record = chained.to_observation(first.request)
    assert record.uncertainty.coverage == "unvalidated"
    assert [(step.variable, dict(step.evidence)) for step in record.lineage] == [
        ("regional_axon_length", {"soma_position": "inferred"}), ("soma_position", {"source_label": "measured"})]
    assert [f.variable for f in lenient.ledger.ancestry(record)] == ["regional_axon_length", "soma_position"]


def test_an_inferred_input_that_repeats_the_request_never_conditions(registry, layout, tasks, planted, zero_inflated):
    """A value inferred only from measurements that the request supplies anyway adds nothing. Conditioning on it
    would count the same evidence twice (and treat a point estimate as a measurement), so it is ignored, said
    so, and the answer is exactly the answer without it."""
    lenient, _ = reconstructor(registry, layout, tasks, planted, policy=core.Policy(on_not_validated="engine_flagged"))
    _, soma, y, label = planted_rows(3, 84, prefix="f")
    first = lenient.reconstruct(request(registry, {"soma_position": soma[0]}, ["regional_axon_length"]))
    profile = first["regional_axon_length"].to_observation(first.request)        # from the measured soma only
    alone = lenient.reconstruct(request(registry, {"soma_position": soma[0]}, ["source_label"]))["source_label"]
    both = lenient.reconstruct(request(registry, {"soma_position": soma[0], "regional_axon_length": profile},
                                       ["source_label"]))["source_label"]
    assert both.estimate.probabilities == alone.estimate.probabilities          # identical, not merely close
    assert dict(both.evidence_used) == {"soma_position": "measured"} and both.conditional is None
    assert [v for v, _ in both.evidence_ignored] == ["regional_axon_length"]
    notice = next(n for n in both.notices if n.code == "redundant_inferred_input")
    assert notice.variable == "regional_axon_length" and notice.handling == "ignored"
    assert "conditioned_on_inferred_input" not in both.codes()
    # the same holds for a label inferred from the soma: with the measured soma it does not turn the request
    # into the benchmarked conditional 'profile from soma + MEASURED label'
    label_hop = lenient.reconstruct(request(registry, {"soma_position": soma[0]}, ["source_label"]))
    label_record = label_hop["source_label"].to_observation(label_hop.request)
    strict, _ = reconstructor(registry, layout, tasks, planted, ledger=lenient.ledger)
    answer = strict.reconstruct(request(registry, {"soma_position": soma[0], "source_label": label_record},
                                        ["regional_axon_length"]))["regional_axon_length"]
    assert answer.status == "refused" and answer.refusal.code == "not_validated" and answer.conditional is None
    assert [v for v, _ in answer.evidence_ignored] == ["source_label"]
    # with only PART of the soma measured the inferred profile does carry more, so it is used (and flagged)
    partial = lenient.reconstruct(request(registry, {"soma_position": {"AP": float(soma[0][0])},
                                                     "regional_axon_length": profile}, ["source_label"]))["source_label"]
    assert dict(partial.evidence_used) == {"soma_position": "measured", "regional_axon_length": "inferred"}
    assert "conditioned_on_inferred_input" in partial.codes()
    # under the presence / level encoding an inferred profile is a mean over 'absent' and 'present': it cannot
    # be read as an observation of either, so it is refused as evidence
    hurdle = mods.JointMixtureModule(registry, layout, states={n: zero_inflated["hurdle"] for n in tasks},
                                     train=zero_inflated["rows"], tasks=tasks,
                                     evidence=evidence_for(tasks, "exploratory_supported"))
    zr = core.Reconstructor(registry, core.Policy(on_not_validated="engine_flagged"))
    zr.register(hurdle)
    made = zr.reconstruct(request(registry, {"source_label": label[0]}, ["regional_axon_length"]))
    mean_profile = made["regional_axon_length"].to_observation(made.request)
    refused = zr.reconstruct(request(registry, {"regional_axon_length": mean_profile}, ["soma_position"]))["soma_position"]
    assert refused.status == "refused" and refused.refusal.code == "invalid_request"
    assert "presence / level" in refused.refusal.reason


def test_an_inferred_input_cannot_launder_the_profile_or_the_target(registry, layout, tasks, planted):
    """T1 is the hidden-axon forecast. A soma or a label that was inferred FROM the measured profile carries the
    profile in its ancestry and is refused there; and no evidence may descend from the requested property."""
    r, _ = reconstructor(registry, layout, tasks, planted)
    _, soma, y, label = planted_rows(3, 86, prefix="f")
    profile = lengths(y[0])
    from_profile = r.reconstruct(request(registry, {"regional_axon_length": profile, "source_label": label[0]},
                                         ["soma_position"]))
    soma_record = from_profile["soma_position"].to_observation(from_profile.request)
    assert from_profile["soma_position"].coverage_validated and soma_record.status == "inferred"
    laundered = r.reconstruct(request(registry, {"soma_position": soma_record, "source_label": label[0]},
                                      ["regional_axon_length"]))["regional_axon_length"]
    assert laundered.status == "refused" and laundered.refusal.code == "forbidden_input"
    assert "derived from regional_axon_length" in laundered.refusal.reason and laundered.estimate is None
    label_from_both = r.reconstruct(request(registry, {"soma_position": soma[0], "regional_axon_length": profile},
                                            ["source_label"]))
    label_record = label_from_both["source_label"].to_observation(label_from_both.request)
    for target in ("regional_axon_length",):          # the label knows the profile: not an input of the forecast
        answer = r.reconstruct(request(registry, {"soma_position": soma[0], "source_label": label_record},
                                       [target]))[target]
        assert answer.status == "refused" and answer.refusal.code == "forbidden_input"
    # circular even inside a declared profile task: a label inferred from the soma cannot help to infer the soma
    circular = r.reconstruct(request(registry, {"regional_axon_length": profile, "source_label": label_record},
                                     ["soma_position"]))["soma_position"]
    assert circular.status == "refused" and "derived from soma_position" in circular.refusal.reason
    # one unit belongs to one group
    other_group = dataclasses.replace(from_profile.request.observed[1], group=core.Group("elsewhere", "brain_sample", "x"))
    mismatch = r.reconstruct(request(registry, {"regional_axon_length": profile, "source_label": other_group},
                                     ["soma_position"]))["soma_position"]
    assert mismatch.status == "refused" and mismatch.refusal.code == "invalid_request"


# ----------------------------------------------------------------------------------------------
# structured refusals
# ----------------------------------------------------------------------------------------------


def test_unsupported_variable(registry, layout, tasks, planted):
    r, _ = reconstructor(registry, layout, tasks, planted)
    _, soma, _, label = planted_rows(1, 97, prefix="f")
    observed = {"soma_position": soma[0], "source_label": label[0]}
    result = r.reconstruct(request(registry, observed, ["firing_rate", "hemisphere", "regional_axon_length"]))
    assert result["firing_rate"].refusal.code == "unsupported_variable"           # not a registered variable
    assert result["hemisphere"].refusal.code == "unsupported_variable"            # registered, no module produces it
    assert result["hemisphere"].estimate is None
    assert result["regional_axon_length"].status == "inferred"                     # refusals are per property


def test_gene_expression_of_a_single_jiao_neuron_is_unidentifiable(registry, layout, tasks, planted):
    r, _ = reconstructor(registry, layout, tasks, planted)
    _, soma, _, label = planted_rows(1, 96, prefix="f")
    answer = r.reconstruct(request(registry, {"soma_position": soma[0], "source_label": label[0]},
                                   ["gene_expression"]))["gene_expression"]
    assert answer.status == "refused" and answer.refusal.code == "unidentifiable_cross_modal"
    detail = " ".join(answer.refusal.detail)
    assert "hypothesis_only" in detail and "forbidden" in detail and answer.estimate is None


def test_evidence_from_another_cell_or_dataset_is_refused(registry, layout, tasks, planted):
    r, _ = reconstructor(registry, layout, tasks, planted)
    _, soma, _, label = planted_rows(1, 95, prefix="f")
    neighbour = _record("jiao", "neuron", "some_other_neuron", "soma_position", tuple(soma[0]), core.SOMA_AXES)
    foreign = _record("epiretro", "nucleus", "fresh0001", "retrograde_target", "VTA")
    for intruder in (neighbour, foreign):
        result = r.reconstruct(request(registry, {"x": intruder, "source_label": label[0]}, ["regional_axon_length"]))
        answer = result["regional_axon_length"]
        assert answer.refusal.code == "unidentifiable_cross_modal" and "never joined as paired" in answer.refusal.reason


def test_disconnected_measurement_blocks_are_unidentifiable(registry, layout, tasks):
    """Units that measured soma + label and OTHER units that measured the profile: no conditional links them."""
    rows, soma, y, label = planted_rows(400, 21)
    half = np.arange(400) < 200
    for name in rows.cont:
        hide = half if name.startswith("regional") else ~half
        rows.cont[name] = np.where(hide, np.nan, rows.cont[name])
    rows.cat["source_label"] = np.array([l if h else None for l, h in zip(label, half)], dtype=object)
    state, _ = mods.fit_joint(rows, layout, components=2, seed=1, budget={"n_init": 2, "max_iter": 50, "tol": 1e-6})
    assert len(state.blocks) == 2
    module = mods.JointMixtureModule(registry, layout, states={n: state for n in tasks}, train=rows, tasks=tasks,
                                     evidence=evidence_for(tasks, "exploratory_supported"))
    r = core.Reconstructor(registry)
    r.register(module)
    answer = r.reconstruct(request(registry, {"soma_position": soma[0], "source_label": label[0]},
                                   ["regional_axon_length"]))["regional_axon_length"]
    assert answer.status == "refused" and answer.refusal.code == "unidentifiable_cross_modal"
    assert "disconnected" in answer.refusal.reason


def test_evidence_that_cannot_condition_is_ignored_and_said_so(registry, layout, tasks):
    """The label was only ever measured on OTHER units than soma and profile: it cannot condition the answer.
    It must appear as ignored, never as used, and the coverage claim is dropped."""
    rows, soma, y, label = planted_rows(600, 22)
    label_only = np.arange(600) >= 400
    for name in rows.cont:
        rows.cont[name] = np.where(label_only, np.nan, rows.cont[name])
    rows.cat["source_label"] = np.array([l if only else None for l, only in zip(label, label_only)], dtype=object)
    state, _ = mods.fit_joint(rows, layout, components=3, seed=1, budget={"n_init": 2, "max_iter": 50, "tol": 1e-6})
    assert len(state.blocks) == 2
    module = mods.JointMixtureModule(registry, layout, states={n: state for n in tasks}, train=rows, tasks=tasks,
                                     evidence=evidence_for(tasks, "exploratory_supported"))
    r = core.Reconstructor(registry)
    r.register(module)
    answer = r.reconstruct(request(registry, {"soma_position": soma[0], "source_label": label[0]},
                                   ["regional_axon_length"]))["regional_axon_length"]
    assert answer.status == "inferred" and "evidence_not_linked" in answer.codes() and not answer.coverage_validated
    assert dict(answer.evidence_used) == {"soma_position": "measured"}
    assert [v for v, _ in answer.evidence_ignored] == ["source_label"]
    assert not set(dict(answer.evidence_used)) & {v for v, _ in answer.evidence_ignored}
    query = mods.Rows(np.array(["q"]), np.array(["fresh_sample"]),
                      {f"soma_position:{a}": np.array([soma[0][j]]) for j, a in enumerate(core.SOMA_AXES)},
                      {"source_label": np.array([None], dtype=object)})
    np.testing.assert_allclose(answer.estimate.point, mods.engine_predict(state, tasks["T1"], query)["point"][0],
                               rtol=1e-12)


def test_undeclared_conditional_follows_the_declared_policy(registry, layout, tasks, planted):
    _, soma, _, label = planted_rows(1, 94, prefix="f")
    ask = request(registry, {"soma_position": soma[0]}, ["regional_axon_length"])       # no label: not a declared task
    refuse, _ = reconstructor(registry, layout, tasks, planted)
    answer = refuse.reconstruct(ask)["regional_axon_length"]
    assert answer.status == "refused" and answer.refusal.code == "not_validated" and answer.conditional is None
    baseline, _ = reconstructor(registry, layout, tasks, planted, policy=core.Policy(on_not_validated="baseline"))
    answer = baseline.reconstruct(ask)["regional_axon_length"]
    assert answer.status == "inferred" and answer.baseline_level and "not_validated" in answer.codes()
    assert answer.module == "baseline:MEAN" and not answer.coverage_validated
    engine, _ = reconstructor(registry, layout, tasks, planted, policy=core.Policy(on_not_validated="engine_flagged"))
    answer = engine.reconstruct(ask)["regional_axon_length"]
    assert answer.status == "inferred" and "not_validated" in answer.codes() and not answer.coverage_validated
    assert answer.module == "joint_mixture" and not answer.baseline_level


def test_no_better_than_baseline_returns_flagged_baseline_or_refuses(registry, layout, tasks, planted):
    _, soma, _, label = planted_rows(1, 93, prefix="f")
    ask = request(registry, {"soma_position": soma[0], "source_label": label[0]}, ["regional_axon_length"])
    default, _ = reconstructor(registry, layout, tasks, planted, status="no_better_than_baseline")
    answer = default.reconstruct(ask)["regional_axon_length"]
    assert answer.status == "inferred" and answer.baseline_level and answer.module == "baseline:IDENTITY"
    assert "no_better_than_baseline" in answer.codes() and answer.validation == "no_better_than_baseline"
    bank = mods.TaskBaselines(tasks["T1"], planted["rows"])
    query = mods.Rows(np.array(["q"]), np.array(["fresh_sample"]),
                      {f"soma_position:{a}": np.array([soma[0][j]]) for j, a in enumerate(core.SOMA_AXES)},
                      {"source_label": np.array([label[0]], dtype=object)})
    np.testing.assert_allclose(answer.estimate.point, bank.predict("IDENTITY", query)[0])
    strict, _ = reconstructor(registry, layout, tasks, planted, status="no_better_than_baseline",
                              policy=core.Policy(on_no_better_than_baseline="refuse"))
    answer = strict.reconstruct(ask)["regional_axon_length"]
    assert answer.status == "refused" and answer.refusal.code == "no_better_than_baseline" and answer.estimate is None


def test_unseen_source_label_is_marginalised_and_said_so(registry, layout, tasks, planted):
    _, soma, _, _ = planted_rows(1, 92, prefix="f")
    r, _ = reconstructor(registry, layout, tasks, planted)
    answer = r.reconstruct(request(registry, {"soma_position": soma[0], "source_label": "NeverSeen"},
                                   ["regional_axon_length"]))["regional_axon_length"]
    assert answer.status == "inferred" and "out_of_support" in answer.codes() and not answer.coverage_validated
    notice = next(n for n in answer.notices if n.code == "out_of_support")
    assert notice.variable == "source_label" and notice.handling == "marginalised"
    assert dict(answer.evidence_used) == {"soma_position": "measured"}
    assert [v for v, _ in answer.evidence_ignored] == ["source_label"]       # supplied, not used: listed, with why
    soma_only, _ = reconstructor(registry, layout, tasks, planted, policy=core.Policy(on_not_validated="engine_flagged"))
    reference = soma_only.reconstruct(request(registry, {"soma_position": soma[0]}, ["regional_axon_length"]))
    np.testing.assert_allclose(answer.estimate.point, reference["regional_axon_length"].estimate.point, rtol=1e-12)
    strict, _ = reconstructor(registry, layout, tasks, planted, policy=core.Policy(on_unseen_category="refuse"))
    refused = strict.reconstruct(request(registry, {"soma_position": soma[0], "source_label": "NeverSeen"},
                                         ["regional_axon_length"]))["regional_axon_length"]
    assert refused.status == "refused" and refused.refusal.code == "out_of_support"


def test_input_far_outside_the_training_range_is_out_of_support(registry, layout, tasks, planted):
    r, _ = reconstructor(registry, layout, tasks, planted)
    ask = request(registry, {"soma_position": (100.0, 100.0, 100.0), "source_label": "A"}, ["regional_axon_length"])
    assert r.reconstruct(ask)["regional_axon_length"].refusal.code == "out_of_support"
    lenient, _ = reconstructor(registry, layout, tasks, planted, policy=core.Policy(on_out_of_range="marginalise"))
    answer = lenient.reconstruct(ask)["regional_axon_length"]
    assert answer.status == "inferred" and "out_of_support" in answer.codes() and not answer.coverage_validated
    assert dict(answer.evidence_used) == {"source_label": "measured"}
    assert [v for v, _ in answer.evidence_ignored] == ["soma_position"]
    assert not set(dict(answer.evidence_used)) & {v for v, _ in answer.evidence_ignored}


def test_a_variable_derived_from_the_target_cannot_be_its_evidence(registry, layout, tasks, planted):
    r, _ = reconstructor(registry, layout, tasks, planted)
    _, _, y, label = planted_rows(1, 91, prefix="f")
    answer = r.reconstruct(request(registry, {"regional_axon_length": lengths(y[0]), "source_label": label[0],
                                              "soma_label_fine": "ARH"}, ["soma_position"]))["soma_position"]
    assert answer.status == "refused" and answer.refusal.code == "forbidden_input"
    banned = r.reconstruct(request(registry, {"source_label": label[0], "projectome_cluster": "c3"},
                                   ["soma_position"]))["soma_position"]
    assert banned.refusal.code == "forbidden_input"


def test_units_and_groups_of_the_state_are_guarded(registry, layout, tasks, planted):
    r, _ = reconstructor(registry, layout, tasks, planted)
    rows = planted["rows"]
    soma = [float(rows.cont[f"soma_position:{a}"][0]) for a in core.SOMA_AXES]
    label = rows.cat["source_label"][0]
    observed = {"soma_position": soma, "source_label": label}
    unit, group = str(rows.unit[0]), str(rows.group[0])
    ask = lambda u, g: r.reconstruct(request(registry, observed, ["regional_axon_length"], unit=u,
                                             group=g))["regional_axon_length"]
    inside = ask(unit, group)
    assert inside.status == "refused" and inside.refusal.code == "leakage_guard"
    sibling = ask(f"{group}_99999", group)
    assert sibling.status == "inferred" and "group_in_training_state" in sibling.codes()
    assert not sibling.coverage_validated
    fresh = ask("fresh_0001", "fresh")
    assert fresh.status == "inferred" and fresh.coverage_validated and not fresh.codes()
    # an identifier is never made 'new' by its spelling: blanks are refused, another case still matches the
    # training unit or sample, and where the unit ID names its group (Jiao) the stated group must be that one
    for bad_unit, bad_group in ((unit + " ", group), (f"{group}_99999", group + " "), (f"{group}_99999", ""),
                                (f"{group}_99999", "another_sample"), (" " + unit, group)):
        refused = ask(bad_unit, bad_group)
        assert refused.status == "refused" and refused.refusal.code == "invalid_request", (bad_unit, bad_group)
    shouted = ask(unit.upper(), group.upper())
    assert shouted.status == "refused" and shouted.refusal.code == "leakage_guard"
    shouted_sibling = ask(f"{group.upper()}_99999", group.upper())
    assert "group_in_training_state" in shouted_sibling.codes() and not shouted_sibling.coverage_validated
    # no group stated: unknown, flagged; but a unit whose ID names a training sample is still recognised
    unknown = r.reconstruct(core.ReconstructionRequest.from_values(
        registry, namespace="jiao", unit_id="fresh_0001", observed=observed,
        requested=["regional_axon_length"]))["regional_axon_length"]
    assert "group_unknown" in unknown.codes() and not unknown.coverage_validated
    unstated = r.reconstruct(core.ReconstructionRequest.from_values(
        registry, namespace="jiao", unit_id=f"{group}_99999", observed=observed,
        requested=["regional_axon_length"]))["regional_axon_length"]
    assert "group_in_training_state" in unstated.codes() and not unstated.coverage_validated
    # a request whose own records name a group, although the request itself names none, is taken at their word
    record_says = core.ReconstructionRequest(
        "jiao", "neuron", f"{group}_99999", core.NO_GROUP,
        (measured("soma_position", tuple(soma), core.SOMA_AXES, unit=f"{group}_99999", group=group),
         measured("source_label", label, unit=f"{group}_99999", group=group)),
        (core.Target("regional_axon_length"),))
    assert "group_in_training_state" in r.reconstruct(record_says)["regional_axon_length"].codes()
    # asking twice for one variable, or for a coordinate twice, is not answered
    twice = r.reconstruct(request(registry, observed, [("regional_axon_length", REGIONS[:2]),
                                                       ("regional_axon_length", REGIONS[2:])]))
    assert twice["regional_axon_length"].refusal.code == "invalid_request"
    repeated = r.reconstruct(request(registry, observed, [("regional_axon_length", ("ARH", "ARH"))]))
    assert repeated["regional_axon_length"].refusal.code == "invalid_request"


def test_copy_donors_never_come_from_the_query_brain_sample_whatever_its_spelling(registry, layout, tasks, planted):
    """A baseline-level COPY answer takes the nearest training neuron of ANOTHER brain sample. A new neuron of a
    training sample sits exactly on a training neuron of that sample here; that neuron must not be its donor, also
    when the request spells the sample in another case or does not state it at all."""
    rows = planted["rows"]
    module = mods.JointMixtureModule(registry, layout, states={n: planted["state"] for n in tasks}, train=rows,
                                     tasks=tasks, evidence=evidence_for(tasks, "no_better_than_baseline", best="COPY"))
    r = core.Reconstructor(registry)
    r.register(module)
    group = str(rows.group[0])
    soma = [float(rows.cont[f"soma_position:{a}"][0]) for a in core.SOMA_AXES]
    label = str(rows.cat["source_label"][0])
    own_profile = np.array([rows.cont[f"regional_axon_length:{c}"][0] for c in REGIONS])     # the neuron it sits on
    bank = mods.TaskBaselines(tasks["T1"], rows)
    expected = bank.predict("COPY", query_rows(soma, label, group=group))[0]                 # another sample's neuron
    assert not np.allclose(expected, own_profile)
    np.testing.assert_allclose(bank.predict("COPY", query_rows(soma, label, group="fresh"))[0], own_profile)

    def copied(unit, stated):
        ask = core.ReconstructionRequest.from_values(registry, namespace="jiao", unit_id=unit, group_id=stated,
                                                     observed={"soma_position": soma, "source_label": label},
                                                     requested=["regional_axon_length"])
        answer = r.reconstruct(ask)["regional_axon_length"]
        assert answer.module == "baseline:COPY" and answer.baseline_level, answer
        return answer

    for unit, stated in ((f"{group}_99999", group), (f"{group.upper()}_99999", group.upper()), (f"{group}_99999", None)):
        answer = copied(unit, stated)
        np.testing.assert_allclose(answer.estimate.point, expected)
        assert "group_in_training_state" in answer.codes() and not answer.coverage_validated
    np.testing.assert_allclose(copied("fresh_0001", "fresh").estimate.point, own_profile)     # a truly new sample may


def test_requested_coordinates_in_any_order_get_the_same_values(registry, layout, tasks, planted):
    """A declared task answers in its own coordinate order; the names and the values stay together whatever the
    order of the request, for the engine and for a baseline-level answer, and in the record made from it."""
    import itertools

    _, soma, y, label = planted_rows(1, 71, prefix="f")
    hidden = REGIONS[3:]
    t4 = {"regional_axon_length": lengths(y[0][:3], INSIDE), "soma_position": soma[0], "source_label": label[0]}
    t1 = {"soma_position": soma[0], "source_label": label[0]}
    for status in ("exploratory_supported", "no_better_than_baseline"):
        r, _ = reconstructor(registry, layout, tasks, planted, status=status)
        for observed, coordinates in ((t4, hidden), (t1, REGIONS)):
            reference, orders = None, list(itertools.permutations(coordinates))
            for order in orders[::max(1, len(orders) // 8)]:
                ask = request(registry, observed, [("regional_axon_length", order)])
                answer = r.reconstruct(ask)["regional_axon_length"]
                assert answer.status == "inferred" and answer.baseline_level == (status != "exploratory_supported")
                assert set(answer.estimate.coordinates) == set(coordinates) == set(answer.coordinates)
                values = dict(zip(answer.estimate.coordinates, answer.estimate.point))
                reference = reference or values
                assert values == reference, (status, order)
                record = answer.to_observation(ask)
                assert dict(zip(record.coordinates, record.value)) == reference
    # an undeclared conditional (a subset of the regions) keeps the order of the request, names and values together
    lenient, _ = reconstructor(registry, layout, tasks, planted, policy=core.Policy(on_not_validated="engine_flagged"))
    forward = lenient.reconstruct(request(registry, t1, [("regional_axon_length", ("ARH", "PAG"))]))["regional_axon_length"]
    backward = lenient.reconstruct(request(registry, t1, [("regional_axon_length", ("PAG", "ARH"))]))["regional_axon_length"]
    assert forward.estimate.coordinates == ("ARH", "PAG") and backward.estimate.coordinates == ("PAG", "ARH")
    assert dict(zip(forward.estimate.coordinates, forward.estimate.point)) == \
        dict(zip(backward.estimate.coordinates, backward.estimate.point))


def test_comparator_rows_must_be_the_rows_of_the_state(registry, layout, tasks, planted):
    """A baseline-level answer comes from the comparator bank. Its leakage guard and its provenance are those of
    the state only if the bank holds exactly the state's rows, so anything else is refused at construction; and
    the answer reports the bank's own training groups and a hash of the comparator that answered."""
    other_rows, _, _, _ = planted_rows(300, 5, prefix="x", groups=3)
    with pytest.raises(core.ReconError, match="not the rows state"):
        mods.JointMixtureModule(registry, layout, states={n: planted["state"] for n in tasks}, train=other_rows,
                                tasks=tasks)
    fewer = planted["rows"].take(np.arange(len(planted["rows"]) - 1))
    with pytest.raises(core.ReconError, match="not the rows state"):
        mods.JointMixtureModule(registry, layout, states={n: planted["state"] for n in tasks}, train=fewer, tasks=tasks)
    other_layout = mods.jiao_layout(registry, REGIONS[:5])
    with pytest.raises(core.ReconError, match="another field layout"):
        mods.JointMixtureModule(registry, other_layout, states={n: planted["state"] for n in tasks},
                                train=planted["rows"], tasks=mods.build_tasks(other_layout, registry, INSIDE))
    r, module = reconstructor(registry, layout, tasks, planted, status="no_better_than_baseline")
    rows = planted["rows"]
    observed = {"soma_position": [float(rows.cont[f"soma_position:{a}"][0]) for a in core.SOMA_AXES],
                "source_label": rows.cat["source_label"][0]}
    own = r.reconstruct(request(registry, observed, ["regional_axon_length"], unit=str(rows.unit[0]),
                                group=str(rows.group[0])))["regional_axon_length"]
    assert own.status == "refused" and own.refusal.code == "leakage_guard"      # a row of the bank is not answered
    _, soma, _, label = planted_rows(1, 72, prefix="f")
    answer = r.reconstruct(request(registry, {"soma_position": soma[0], "source_label": label[0]},
                                   ["regional_axon_length"]))["regional_axon_length"]
    bank = mods.TaskBaselines(tasks["T1"], rows)
    assert answer.baseline_level and answer.module == "baseline:IDENTITY"
    assert answer.training_groups_hash == bank.groups_sha256 == planted["state"].groups_sha256
    assert answer.state_hash == bank.answer_hash("IDENTITY") != planted["state"].state_hash()
    assert module.state_groups(answer.state_hash) == frozenset(bank.groups)
    # the module keeps its own copy of the rows, and that copy cannot be written to ...
    assert module.train is not rows and not module.train.unit.flags.writeable
    assert not module.train.cont["soma_position:AP"].flags.writeable
    with pytest.raises(ValueError):
        module.train.cont["soma_position:AP"][0] = 0.0
    # ... and rows that are swapped in behind the module are noticed when a comparator bank is built from them
    module.train = other_rows.frozen_copy()
    _, _, y, _ = planted_rows(1, 72, prefix="f")
    with pytest.raises(core.FrozenSourceChanged, match="no longer the rows"):
        r.reconstruct(request(registry, {"regional_axon_length": lengths(y[0]), "source_label": label[0]},
                              ["soma_position"]))
    assert not hasattr(core.ProvenanceLedger, "issue")              # only reconstruct() writes to the ledger


def test_measured_profile_is_an_input_only_in_declared_tasks(registry, layout, tasks, planted):
    engine, _ = reconstructor(registry, layout, tasks, planted, policy=core.Policy(on_not_validated="engine_flagged"))
    _, soma, y, label = planted_rows(1, 90, prefix="f")
    # declared: soma from the full profile + label (T2); label from soma + profile (T3); outside from inside (T4)
    t2 = engine.reconstruct(request(registry, {"regional_axon_length": lengths(y[0]), "source_label": label[0]},
                                    ["soma_position"]))["soma_position"]
    t3 = engine.reconstruct(request(registry, {"regional_axon_length": lengths(y[0]), "soma_position": soma[0]},
                                    ["source_label"]))["source_label"]
    t4 = engine.reconstruct(request(registry, {"regional_axon_length": lengths(y[0][:3], INSIDE),
                                               "soma_position": soma[0], "source_label": label[0]},
                                    ["regional_axon_length"]))["regional_axon_length"]
    assert (t2.conditional, t3.conditional) == ("T2_soma_from_profile_label", "T3_label_from_soma_profile")
    assert t4.conditional == "T4_outside_block_from_inside_block_soma_label" and t4.estimate.coordinates == REGIONS[3:]
    assert t2.status == t3.status == t4.status == "inferred"
    assert set(t3.estimate.probabilities) == set(LABELS) and t3.estimate.point in t3.estimate.prediction_set
    # undeclared conditionals that would take a measured profile are refused whatever the policy
    partial = engine.reconstruct(request(registry, {"regional_axon_length": lengths(y[0][:5], REGIONS[:5]),
                                                    "soma_position": soma[0], "source_label": label[0]},
                                         ["regional_axon_length"]))["regional_axon_length"]
    profile_only = engine.reconstruct(request(registry, {"regional_axon_length": lengths(y[0])},
                                              ["source_label"]))["source_label"]
    for answer in (partial, profile_only):
        assert answer.status == "refused" and answer.refusal.code == "forbidden_input"
    # asking for something that was supplied is not a reconstruction
    again = engine.reconstruct(request(registry, {"regional_axon_length": lengths(y[0]), "source_label": label[0]},
                                       ["regional_axon_length"]))["regional_axon_length"]
    assert again.refusal.code == "invalid_request"


def test_modules_register_what_they_consume_and_produce_with_evidence(registry, layout, tasks, planted):
    r, module = reconstructor(registry, layout, tasks, planted)
    listing = {m["name"]: m for m in r.describe()}
    joint = listing["joint_mixture"]
    assert (joint["namespace"], joint["grain"]) == ("jiao", "neuron")
    assert set(joint["produces"]) == set(joint["consumes"]) == {"soma_position", "source_label", "regional_axon_length"}
    by_task = {c["task"]: c for c in joint["conditionals"]}
    assert set(by_task) == {"T1", "T2", "T3", "T4"}
    assert by_task["T1"]["validation"]["status"] == "exploratory_supported"
    assert by_task["T1"]["validation"]["benchmark"] == POINTER and not by_task["T1"]["validation"]["confirmatory"]
    assert not by_task["T1"]["declared_measured_axon_input"] and by_task["T2"]["declared_measured_axon_input"]
    assert listing["axon_arbor_adapter"]["conditionals"][0]["validation"]["status"].startswith("blocked")
    with pytest.raises(core.RegistryError):
        r.register(module)
    with pytest.raises(core.RegistryError):
        module.set_evidence("not_a_conditional", core.ValidationEvidence())


# ----------------------------------------------------------------------------------------------
# sealed rows
# ----------------------------------------------------------------------------------------------


def _manifest_rows(world, part):
    with open(world.manifest, newline="") as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t") if row["part"] == part]


def test_sealed_rows_are_dropped_unparsed_and_forbidden_columns_never_read(bare_world):
    world = bare_world
    partition = core.TrainingPartition(world)
    table = core.load_training_table(world, partition=partition)          # sealed rows hold unparseable values
    held = _manifest_rows(world, "validation") + _manifest_rows(world, "test")
    assert held and not {r["animal"] for r in held} & partition.train_samples
    assert not {r["neuron"] for r in held} & set(table.neuron.tolist())
    assert set(table.sample.tolist()) == set(partition.train_samples)
    assert len(table) == dict(world.expected)["train_neurons"] and table.length.shape == (len(table), len(world.regions))
    for audit in (table.audit["soma_csv"], table.audit["soma_xlsx"], table.audit["regional_csv"]):
        assert audit["training_rows"] == len(table) and audit["sealed_rows_dropped_unparsed"] == len(held)
    everything = " ".join(map(str, table.fine_label.tolist() + table.csv_location.tolist() + table.label.tolist()
                              + table.full_id.tolist()))
    assert "FORBIDDEN" not in everything and bench.SEALED not in everything
    assert np.isfinite(table.soma).all() and (table.soma[:, 2] <= core.MIDLINE_ML).all()
    assert "Wt(orexin)" in table.label                                     # table spelling kept, alias checked
    inputs_only = core.load_training_table(world, with_outcomes=False, partition=partition)
    assert inputs_only.length is None and inputs_only.audit["regional_csv"] is None
    with pytest.raises(core.ReconError):
        inputs_only.profile()


def test_a_canary_with_a_validation_sample_id_is_rejected(bare_world):
    world = bare_world
    partition = core.TrainingPartition(world)
    canary = _manifest_rows(world, "validation")[0]
    sample, neuron = canary["animal"], canary["neuron"]
    assert not partition.is_training(sample)
    with pytest.raises(core.SealedRowError):
        partition.require_training([sample], "canary")
    soma_rows, _ = core.read_soma_training(world, partition)
    fine, _ = core.read_fine_labels_training(world, partition)
    lengths_, _ = core.read_regional_training(world, partition)
    assert neuron not in soma_rows and neuron not in fine and not any(k.startswith(neuron) for k in lengths_)
    full = f"{neuron}_Pomc_ARH"
    planted_soma = dict(soma_rows)
    planted_soma[neuron] = {"full_id": full, "geno": "Pomc", "location": "ARH", "soma": (6000.0, 6000.0, 5000.0)}
    with pytest.raises(core.SealedRowError):
        core.assemble_training_table(world, partition, planted_soma, fine, lengths_)
    planted_lengths = dict(lengths_)
    planted_lengths[full] = np.zeros(len(world.regions))
    with pytest.raises(core.SealedRowError):
        core.assemble_training_table(world, partition, soma_rows, fine, planted_lengths)
    planted_fine = dict(fine)
    planted_fine[neuron] = "ARH"
    with pytest.raises(core.SealedRowError):
        core.assemble_training_table(world, partition, soma_rows, planted_fine, lengths_)
    ctx = bench.build_context(world)
    ctx.rows.group[0] = sample                       # a sealed sample smuggled into the rows of a fit
    with pytest.raises(core.SealedRowError):
        bench._fit(ctx, np.arange(len(ctx.rows)), "log", 1, ("canary",))
    outer = int(ctx.outer_of[0] + 1) % ctx.config["outer_folds"]       # folds whose training rows hold row 0
    inner = int(ctx.inner_of[outer][0] + 1) % ctx.config["inner_folds"]
    with pytest.raises(core.SealedRowError):
        bench.job_baseline_inner(ctx, outer, inner)
    with pytest.raises(core.SealedRowError):
        bench.job_baseline_outer(ctx, outer, {})


def test_changing_sealed_rows_changes_nothing(bare_world, tmp_path):
    world = bare_world
    before = core.load_training_table(world)
    copy = tmp_path / "scrambled"
    shutil.copytree(Path(world.manifest).parent, copy)
    rng = np.random.default_rng(0)
    for name in ("regional.csv", "soma.csv"):
        lines = (copy / name).read_text().splitlines()
        out = [lines[0]]
        for line in lines[1:]:
            while bench.SEALED in line:
                line = line.replace(bench.SEALED, f"{rng.uniform(0, 9000):.3f}", 1)
            out.append(line)
        (copy / name).write_text("\n".join(out) + "\n")
    after = core.load_training_table(core.synthetic_world(copy))
    assert core.sha_file(copy / "regional.csv") != core.sha_file(world.regional_csv)
    assert after.sha256() == before.sha256()
    np.testing.assert_array_equal(after.length, before.length)
    np.testing.assert_array_equal(after.soma, before.soma)


def test_loaders_refuse_bad_tables(bare_world, tmp_path):
    import openpyxl

    world = bare_world
    partition = core.TrainingPartition(world)
    # a forbidden column lying between the two allowed xlsx columns: refuse to read across it
    bad = tmp_path / "bad.xlsx"
    workbook = openpyxl.Workbook()
    workbook.active.append([None, "neurons", "cluster", "location"])
    workbook.active.append(["x", "300001_001", "FORBIDDEN", "ARH"])
    workbook.save(bad)
    with pytest.raises(core.ReconError):
        core.read_fine_labels_training(dataclasses.replace(world, soma_xlsx=bad), partition)
    # pinned hashes and expected counts
    with pytest.raises(core.FrozenSourceChanged):
        core.TrainingPartition(dataclasses.replace(world, pins=(("manifest", "0" * 64),)))
    with pytest.raises(core.FrozenSourceChanged):
        core.read_regional_training(dataclasses.replace(world, pins=(("regional_csv", "0" * 64),)), partition)
    with pytest.raises(core.ReconError):
        core.TrainingPartition(dataclasses.replace(world, expected=(("train_neurons", 1), ("train_samples", 1))))
    # a negative length, or a wrong region header, is not silently accepted
    text = Path(world.regional_csv).read_text().splitlines()
    first_training = next(i for i, line in enumerate(text[1:], 1) if bench.SEALED not in line)
    cells = text[first_training].split(",")
    cells[1] = "-5.0"
    negative = tmp_path / "negative.csv"
    negative.write_text("\n".join(text[:first_training] + [",".join(cells)] + text[first_training + 1:]) + "\n")
    with pytest.raises(core.ReconError):
        core.read_regional_training(dataclasses.replace(world, regional_csv=negative), partition)
    with pytest.raises(core.ReconError):
        core.read_regional_training(dataclasses.replace(world, regions=tuple(reversed(world.regions))), partition)


def test_training_neurons_become_measured_paired_records(bare_world, registry):
    world = bare_world
    table = core.load_training_table(world)
    store = core.jiao_records(world, table, registry)
    assert len(store) == 4 * len(table)
    unit = str(table.neuron[5])
    view = store.unit("jiao", "neuron", unit)
    assert set(view) == {"soma_position", "source_label", "soma_label_fine", "regional_axon_length"}
    for record in view.values():
        assert record.status == "measured" and record.uncertainty.kind == "none_stated"
        assert record.group.kind == "brain_sample" and record.group.group_id == str(table.sample[5])
        assert record.source.sha256 and record.source.row_key in (unit, str(table.full_id[5]))
    assert view["regional_axon_length"].coordinates == tuple(world.regions)
    assert view["regional_axon_length"].value == tuple(float(v) for v in table.length[5])
    assert view["regional_axon_length"].source.path == str(world.regional_csv)
    assert store.join(view.values()).keys() == view.keys()
    other = store.unit("jiao", "neuron", str(table.neuron[6]))
    with pytest.raises(core.PairingError):
        store.join([view["soma_position"], other["regional_axon_length"]])


# ----------------------------------------------------------------------------------------------
# whole-sample splits
# ----------------------------------------------------------------------------------------------


def test_folds_are_whole_samples_nested_and_fixed_by_the_salt():
    samples = [f"s{i:03d}" for i in range(53)]
    label_of = {s: "ABCDE"[i % 5] if i < 50 else "RARE" for i, s in enumerate(samples)}
    splits = core.nested_folds(samples, label_of, outer=5, inner=3, salt="salt-1")
    outer = splits["outer"]
    assert sorted(s for held in outer.values() for s in held) == samples          # each sample in exactly one fold
    assert max(map(len, outer.values())) - min(map(len, outer.values())) <= 1
    for fold, held in outer.items():
        training = set(samples) - set(held)
        inner = splits["inner"][fold]
        assert sorted(s for part in inner.values() for s in part) == sorted(training)   # inner folds: outer-training only
        assert not set(held) & {s for part in inner.values() for s in part}
        for label in "ABCDE":                                                      # labels spread evenly
            assert sum(label_of[s] == label for s in held) == 2
    assert core.nested_folds(list(reversed(samples)), label_of, outer=5, inner=3, salt="salt-1") == splits
    assert core.nested_folds(samples, label_of, outer=5, inner=3, salt="salt-2")["outer"] != outer
    broken = json.loads(json.dumps(splits))
    broken["outer"]["0"].append(broken["outer"]["1"][0])                           # a sample in two folds
    with pytest.raises(core.ReconError):
        core.check_nested_folds(broken, samples)
    leaking = json.loads(json.dumps(splits))
    leaking["inner"]["0"]["0"].append(leaking["outer"]["0"][0])                    # a held-out sample in an inner fold
    with pytest.raises(core.ReconError):
        core.check_nested_folds(leaking, samples)
    with pytest.raises(core.ReconError):
        core.assign_folds(samples[:3], label_of, 5, "salt-1")


def test_benchmark_folds_never_split_a_brain_sample(bare_world):
    ctx = bench.build_context(bare_world)
    group = ctx.rows.group
    for outer in range(ctx.config["outer_folds"]):
        held = set(group[ctx.outer_of == outer])
        assert held and not held & set(group[ctx.outer_of != outer])
        inner = ctx.inner_of[outer]
        assert (inner[ctx.outer_of == outer] == -1).all() and (inner[ctx.outer_of != outer] >= 0).all()
        for i in range(ctx.config["inner_folds"]):
            assert not set(group[inner == i]) & set(group[(inner != i) & (inner >= 0)])
            assert not set(group[inner == i]) & held
    one_label = {s: set(ctx.table.label[ctx.table.sample == s]) for s in set(ctx.table.sample)}
    assert all(len(v) == 1 for v in one_label.values())


# ----------------------------------------------------------------------------------------------
# baselines on a toy set
# ----------------------------------------------------------------------------------------------


def _toy(registry):
    """Six training units in three brain samples, one observed coordinate, one target coordinate."""
    toy_layout = mods.FieldLayout.from_registry(registry, "jiao", ("soma_position", "source_label", "regional_axon_length"),
                                                {"soma_position": ("AP",), "regional_axon_length": ("ARH",)})
    task = mods.make_task(toy_layout, registry, task="toy", title="toy", target="regional_axon_length",
                          evidence=(("soma_position", None), ("source_label", None)), copy_scaling="raw",
                          knn_scalings=("raw",))

    def rows(units, groups, x, labels, y):
        return mods.Rows(np.array(units), np.array(groups), {"soma_position:AP": np.array(x, float),
                                                             "regional_axon_length:ARH": np.array(y, float)},
                         {"source_label": np.array(labels, dtype=object)})

    train = rows(["a1", "a2", "a3", "b1", "c1", "c2"], ["s1", "s1", "s1", "s2", "s3", "s3"],
                 [0.0, 1.0, 2.0, 10.0, 20.0, 21.0], ["A", "A", "A", "A", "B", "B"], [1.0, 2.0, 3.0, 10.0, 40.0, 60.0])
    return toy_layout, task, rows, train


def test_simple_baselines_on_a_toy_set(registry):
    _, task, rows, train = _toy(registry)
    bank = mods.TaskBaselines(task, train)
    query = rows(["q1", "q2", "q3"], ["new", "new", "s3"], [1.4, 18.0, 20.4], ["A", "B", "Z"], [0.0, 0.0, 0.0])
    # MEAN, weighting 'sample' (the default): every brain sample counts once -> mean(2, 10, 50)
    np.testing.assert_allclose(bank.predict("MEAN", query)[:, 0], [62 / 3] * 3)
    np.testing.assert_allclose(bank.predict("MEAN", query, {"weighting": "sample"})[:, 0], [62 / 3] * 3)
    # MEAN, weighting 'neuron': every training neuron counts once -> the row mean 116 / 6
    np.testing.assert_allclose(bank.predict("MEAN", query, {"weighting": "neuron"})[:, 0], [116 / 6] * 3)
    # IDENTITY: label mean over that label's samples (A: mean(2, 10) = 6; B: 50); unseen label -> global fallback
    np.testing.assert_allclose(bank.predict("IDENTITY", query)[:, 0], [6.0, 50.0, 62 / 3])
    # ... or over that label's neurons (A: 16 / 4 = 4; B: 50); unseen label -> the neuron-weighted global mean
    np.testing.assert_allclose(bank.predict("IDENTITY", query, {"weighting": "neuron"})[:, 0], [4.0, 50.0, 116 / 6])
    both = bank.predict_many("MEAN", query, bench.mean_grid(bench.CONFIG))
    np.testing.assert_allclose([p[0, 0] for p in both], [62 / 3, 116 / 6])
    with pytest.raises(ValueError):
        bank.predict("MEAN", query, {"weighting": "median"})
    # COPY: nearest training unit of ANOTHER brain sample (q3 belongs to s3, so c1/c2 are not donors)
    np.testing.assert_allclose(bank.predict("COPY", query)[:, 0], [2.0, 40.0, 10.0])
    # same-label COPY: q2 (label B) copies c1; q3's label is unseen -> falls back to plain COPY
    np.testing.assert_allclose(bank.predict("COPY_SAME_LABEL", query)[:, 0], [2.0, 40.0, 10.0])
    # ties go to the first training unit in the stored order (0.5 is equally far from a1 and a2)
    tie = rows(["t"], ["new"], [0.5], ["A"], [0.0])
    np.testing.assert_allclose(bank.predict("COPY", tie)[:, 0], [1.0])
    # k nearest: q1 -> a2, a3 (k=2); same-label pool for q2 holds only c1, c2
    grid = ({"k": 2, "pool": "all", "scaling": "raw"}, {"k": 5, "pool": "same_label", "scaling": "raw"})
    two, pooled = bank.predict_many("KNN", query, grid)
    np.testing.assert_allclose(two[:, 0], [2.5, 50.0, 6.5])
    np.testing.assert_allclose(pooled[:, 0], [(1 + 2 + 3 + 10) / 4, 50.0, (1 + 2 + 3 + 10) / 4])
    # empirical intervals of the baselines
    lower, upper = bank.interval("IDENTITY", query, 0.5)
    assert lower[0, 0] == np.quantile([1, 2, 3, 10], .25) and upper[1, 0] == np.quantile([40, 60], .75)
    assert bank.interval("COPY", query, 0.5) is None
    assert all(bank.applicable(name) for name in mods.COMPARATORS)
    # the hash reported with a baseline-level answer: the bank, the comparator and every setting it used
    hashes = {bank.answer_hash("MEAN"), bank.answer_hash("MEAN", {"weighting": "neuron"}), bank.answer_hash("IDENTITY"),
              bank.answer_hash("COPY"), bank.answer_hash("COPY_SAME_LABEL"), bank.answer_hash("KNN", grid[0]),
              bank.answer_hash("KNN", grid[1]), bank.answer_hash("RIDGE", {"lambda": 0.1, "degree": 1}),
              bank.answer_hash("RIDGE", {"lambda": 1.0, "degree": 1})}
    assert len(hashes) == 9 and bank.answer_hash("MEAN") == bank.answer_hash("MEAN", {"weighting": "sample"})
    assert mods.TaskBaselines(task, train, vote_alpha=7.0).state_hash() != bank.state_hash()
    relabelled = rows(["a1", "a2", "a3", "b1", "c1", "c2"], ["s1", "s1", "s1", "s2", "s3", "s3"],
                      [0.0, 1.0, 2.0, 10.0, 20.0, 21.0], ["A", "A", "A", "B", "B", "B"], [1.0, 2.0, 3.0, 10.0, 40.0, 60.0])
    assert mods.TaskBaselines(task, relabelled).state_hash() != bank.state_hash()      # the labels are rows too


def test_ridge_reference_matches_the_closed_form(registry, layout, tasks):
    rows, _, _, _ = planted_rows(300, 31)
    task = tasks["T2"]
    bank = mods.TaskBaselines(task, rows.take(np.arange(200)))
    query = rows.take(np.arange(200, 300))
    lam = 0.1
    prediction = bank.predict("RIDGE", query, {"lambda": lam, "degree": 1})

    def design(r):
        z = (r.matrix(task.evidence_cont) - bank.mu) / bank.sd
        onehot = (np.array([str(v) for v in r.cat["source_label"]])[:, None] == np.array(LABELS)[None, :]).astype(float)
        return np.column_stack([z, onehot])

    F, Y = design(rows.take(np.arange(200))), rows.take(np.arange(200)).matrix(task.target_fields)
    centre = F.mean(axis=0)
    spread = np.concatenate([F[:, :6].std(axis=0), np.ones(3)])
    G = (F - centre) / spread
    beta = np.linalg.solve(G.T @ G + lam * 200 * np.eye(G.shape[1]), G.T @ (Y - Y.mean(axis=0)))
    np.testing.assert_allclose(prediction, Y.mean(axis=0) + ((design(query) - centre) / spread) @ beta, rtol=1e-8)
    quadratic = mods.TaskBaselines(tasks["T1"], rows.take(np.arange(200))).predict(
        "RIDGE", query, {"lambda": 1e-3, "degree": 2})
    assert quadratic.shape == (100, 6) and np.isfinite(quadratic).all()


def test_label_baselines_give_probabilities(registry, layout, tasks):
    rows, _, _, label = planted_rows(600, 32)
    task = tasks["T3"]
    train, query = rows.take(np.arange(450)), rows.take(np.arange(450, 600))
    bank = mods.TaskBaselines(task, train)
    assert bank.classes == list(LABELS) and not bank.applicable("IDENTITY") and not bank.applicable("COPY_SAME_LABEL")
    prior = bank.predict("MEAN", query)
    np.testing.assert_allclose(prior[0], mods.sample_weighted_mean(bank.Y, train.group))
    frequencies = bank.predict("MEAN", query, {"weighting": "neuron"})
    np.testing.assert_allclose(frequencies[0], [np.mean(train.cat["source_label"] == c) for c in LABELS])
    truth = query.cat["source_label"]
    accuracy = {}
    for name, params in (("MEAN", None), ("COPY", None), ("KNN", {"k": 15, "pool": "all", "scaling": "z"}),
                         ("RIDGE", {"lambda": 1e-3})):
        p = bank.predict(name, query, params)
        np.testing.assert_allclose(p.sum(axis=1), 1.0)
        assert (p > 0).all()
        accuracy[name] = np.nanmean(core.label_scores(p, bank.classes, truth)["accuracy"])
    assert accuracy["KNN"] > accuracy["MEAN"] + 0.2 and accuracy["RIDGE"] > accuracy["MEAN"] + 0.2
    copy = bank.predict("COPY", query)                 # one vote smoothed by the prior: (1 + prior) / 2
    np.testing.assert_allclose(np.sort(copy, axis=1)[:, -1], (1 + prior[np.arange(150), copy.argmax(axis=1)]) / 2)


# ----------------------------------------------------------------------------------------------
# metrics
# ----------------------------------------------------------------------------------------------


def test_metrics_are_aggregated_per_brain_sample_then_over_samples():
    truth = np.array([[0., 0, 0], [1, 2, 3], [1, 2, 3], [1, 2, 4]])
    prediction = np.array([[1., 1, 1], [1, 2, 3], [2, 2, 2], [3, 2, 1]])
    sample = np.array(["s1", "s2", "s2", "s2"])
    scores = core.profile_scores(prediction, truth)
    np.testing.assert_allclose(scores["mse"], [1, 0, 2 / 3, 13 / 3])
    assert np.isnan(scores["corr"][0])                 # constant truth: undefined, excluded from correlation only
    assert scores["corr"][1] == pytest.approx(1.0) and scores["corr"][2] == 0.0    # constant prediction scores 0
    per = core.by_sample(scores["mse"], sample)
    assert per == {"s1": pytest.approx(1.0), "s2": pytest.approx(5 / 3)}
    assert core.over_samples(per) == pytest.approx((1 + 5 / 3) / 2)                 # not the pooled mean 1.5
    assert scores["mse"].mean() == pytest.approx(1.5)
    assert set(core.by_sample(scores["corr"], sample)) == {"s2"}                    # s1 has no defined correlation
    position = core.position_scores(np.array([[3., 4, 0], [0, 0, 0]]), np.zeros((2, 3)))
    np.testing.assert_allclose(position["euclid"], [5, 0])
    np.testing.assert_allclose(position["abs_AP"], [3, 0])
    assert core.by_sample([1.0, 2.0, 30.0], np.array(["a", "a", "a"]), how="median") == {"a": 2.0}
    labels = core.label_scores(np.array([[.7, .3], [.2, .8], [.5, .5]]), ["A", "B"], ["A", "A", "C"])
    np.testing.assert_allclose(labels["accuracy"][:2], [1, 0])
    np.testing.assert_allclose(labels["logloss"][:2], [-np.log(.7), -np.log(.2)])
    assert np.isnan(labels["accuracy"][2]) and np.isnan(labels["logloss"][2])       # label never seen: unscoreable
    boot = core.paired_bootstrap([-1.0, -2.0, -1.5, -0.5, -1.2, -0.8], seed=3, samples=2000)
    again = core.paired_bootstrap([-1.0, -2.0, -1.5, -0.5, -1.2, -0.8], seed=3, samples=2000)
    assert boot == again and boot["lower"] < boot["mean"] < boot["upper"] < 0 and boot["n"] == 6
    cmp_ = bench.compare("mse", {"a": 1.0, "b": 1.0, "c": 1.0}, {"a": 2.0, "b": 3.0, "c": 1.5, "d": 9.0}, 1,
                         bench.TINY_CONFIG)
    assert cmp_["samples"] == 3 and cmp_["engine_better_strict"] and cmp_["samples_favouring_engine"] == 3
    assert cmp_["engine_minus_comparator"]["mean"] == pytest.approx(-7 / 6) and cmp_["advantage_mean"] == pytest.approx(7 / 6)
    assert cmp_["relative_gain"] == pytest.approx((7 / 6) / (6.5 / 3)) and cmp_["engine_better_clear"]
    higher = bench.compare("accuracy", {"a": .9, "b": .8}, {"a": .5, "b": .6}, 1, bench.TINY_CONFIG)
    assert higher["engine_better_strict"] and higher["advantage_interval"][0] > 0 and higher["engine_better_interval"]
    # a small but perfectly consistent advantage: the interval excludes zero, the declared minimum gain does not
    samples = [f"s{i}" for i in range(40)]
    comparator = {s: 1.0 + 0.01 * i for i, s in enumerate(samples)}
    slightly = bench.compare("mse", {s: v * 0.995 for s, v in comparator.items()}, comparator, 1, bench.TINY_CONFIG)
    assert slightly["engine_better_strict"] and slightly["engine_better_interval"]
    assert slightly["relative_gain"] == pytest.approx(0.005) and not slightly["engine_better_clear"]
    assert bench.CONFIG["min_relative_gain"] == 0.02 == slightly["min_relative_gain"]
    enough = bench.compare("mse", {s: v * 0.97 for s, v in comparator.items()}, comparator, 1, bench.TINY_CONFIG)
    assert enough["engine_better_clear"] and enough["relative_gain"] == pytest.approx(0.03)
    # two predictors that agree to rounding error are tied, whatever the sign of the rounding
    noise = bench.compare("mse", {s: v * (1 - 1e-15) for s, v in comparator.items()}, comparator, 1, bench.TINY_CONFIG)
    assert not noise["engine_better_strict"] and not noise["engine_better_interval"] and noise["advantage_mean"] == 0
    assert bench.compare("mse", {}, comparator, 1, bench.TINY_CONFIG)["samples"] == 0


# ----------------------------------------------------------------------------------------------
# the joint model as a module
# ----------------------------------------------------------------------------------------------


def test_batch_conditionals_equal_the_unchanged_model(layout, tasks, planted):
    """Means, interval endpoints and probabilities equal JointModel.infer(), with missing and unseen inputs."""
    state = planted["state"]
    rows, _, _, _ = planted_rows(40, 55, prefix="q")
    rng = np.random.default_rng(4)
    rows.cat["source_label"][::7] = "NeverSeen"
    rows.cat["source_label"][3::9] = None
    for name in list(rows.cont)[:5]:
        rows.cont[name] = np.where(rng.random(40) < 0.2, np.nan, rows.cont[name])
    for task in tasks.values():
        out = mods.engine_predict(state, task, rows, level=0.9, median=True)
        for i in range(len(rows)):
            observed = {f: float(rows.cont[f][i]) for f in task.evidence_cont if np.isfinite(rows.cont[f][i])}
            observed.update({f: rows.cat[f][i] for f in task.evidence_cat if rows.cat[f][i] is not None})
            reference = state.model.infer(observed, interval_level=0.9)["predictions"]
            if task.kind == "label":
                expected = {d["value"]: d["probability"] for d in reference["source_label"]["probabilities"]}
                np.testing.assert_allclose(out["proba"][i], [expected[c] for c in out["classes"]], atol=1e-12)
                value = out["classes"][0]
                assert np.log(out["proba"][i, 0]) == pytest.approx(
                    state.model.conditional_log_prob("source_label", value, observed), abs=1e-9)
                continue
            for j, name in enumerate(task.target_fields):
                mean, lower, upper = (reference[name]["mean"], reference[name]["interval"]["lower"],
                                      reference[name]["interval"]["upper"])
                if task.clip_min is not None:
                    mean, lower, upper = max(mean, 0.0), max(lower, 0.0), max(upper, 0.0)
                scale = max(1.0, abs(mean))
                assert abs(out["point"][i, j] - mean) <= 1e-9 * scale
                assert abs(out["lower"][i, j] - lower) <= 1e-7 * scale
                assert abs(out["upper"][i, j] - upper) <= 1e-7 * scale
                assert out["lower"][i, j] <= out["median"][i, j] <= out["upper"][i, j]
    unknown = mods.engine_predict(state, tasks["T1"], rows)["unknown"]["source_label"]
    assert unknown.sum() == sum(v == "NeverSeen" for v in rows.cat["source_label"]) > 0


def _encoded(rows, i, fields):
    """One row's evidence as the hurdle model sees it: presence, and the level where the value is above zero."""
    observed = {}
    for name in fields:
        if name in rows.cat:
            observed[name] = rows.cat[name][i]
        elif name.startswith("regional"):
            value = float(rows.cont[name][i])
            observed[name + "|on"] = mods.ON if value > 0 else mods.OFF
            if value > 0:
                observed[name + "|level"] = value
        else:
            observed[name] = float(rows.cont[name][i])
    return observed


@pytest.fixture(scope="module")
def zero_inflated(layout):
    rows, soma, y, label = zero_inflated_rows(3000, 12)
    budget = {"n_init": 5, "max_iter": 200, "tol": 1e-7}
    hurdle, _ = mods.fit_joint(rows, layout, components=4, seed=5, budget=budget, representation="hurdle")
    log, _ = mods.fit_joint(rows, layout, components=4, seed=5, budget=budget, representation="log")
    return {"rows": rows, "hurdle": hurdle, "log": log}


def test_presence_level_representation_is_the_same_unchanged_model(layout, tasks, zero_inflated):
    """'hurdle' only changes how a zero-inflated coordinate is presented: presence (categorical) and level
    (Gaussian, missing at zero). Posterior weights equal JointModel.infer(); means and quantiles of the
    zero-inflated mixture equal what the model's own sampler produces."""
    state = zero_inflated["hurdle"]
    assert state.settings["representation"] == "hurdle" and state.model.max_components == 4
    assert set(state.encoding.values()) == {"gaussian", "hurdle"}
    assert len(mods.model_field_types(layout, state.encoding)) == 3 + 1 + 2 * 6
    assert len(mods.model_field_types(layout, mods.encoding_for(layout, "log"))) == 3 + 1 + 6
    rows, _, y, _ = zero_inflated_rows(30, 61, prefix="q")
    assert (y == 0).mean() > 0.3
    for task in (tasks["T1"], tasks["T2"], tasks["T4"]):
        weights, used, _ = state.posterior(0, rows, set(task.evidence_fields()))
        assert set(used) == set(task.evidence_fields())
        for i in range(len(rows)):
            reference = state.model.infer(_encoded(rows, i, task.evidence_fields()))["blocks"]["block_0"]
            np.testing.assert_allclose(weights[i], reference["posterior_weights"], atol=1e-10)
    # label probabilities from zero-inflated evidence
    t3 = mods.engine_predict(state, tasks["T3"], rows)
    for i in range(len(rows)):
        reference = state.model.infer(_encoded(rows, i, tasks["T3"].evidence_fields()))["predictions"]["source_label"]
        expected = {d["value"]: d["probability"] for d in reference["probabilities"]}
        np.testing.assert_allclose(t3["proba"][i], [expected[c] for c in t3["classes"]], atol=1e-10)
    # the zero-inflated predictive distribution against the model's own sampler
    out = mods.engine_predict(state, tasks["T1"], rows.take([0, 1]), level=0.9, median=True)
    for i in (0, 1):
        draws = state.model.sample_missing(_encoded(rows, i, tasks["T1"].evidence_fields()), n_samples=20000,
                                           seed=7 + i)["samples"]
        for j, name in enumerate(tasks["T1"].target_fields):
            sampled = np.array([d[name + "|level"] if d[name + "|on"] == mods.ON else 0.0 for d in draws])
            assert abs(out["point"][i, j] - np.maximum(sampled, 0).mean()) < 0.03 + 4 * sampled.std() / np.sqrt(len(sampled))
            zero_mass = (sampled <= 0).mean()
            for key, probability in (("lower", 0.05), ("median", 0.5), ("upper", 0.95)):
                if zero_mass > probability + 0.02:
                    assert out[key][i, j] == 0.0                      # the mass at zero reaches the quantile
                elif zero_mass < probability - 0.02:
                    assert out[key][i, j] > 0 and abs((sampled <= out[key][i, j]).mean() - probability) < 0.02
    assert (out["lower"] <= out["median"]).all() and (out["median"] <= out["upper"]).all() and (out["lower"] >= 0).all()


def test_presence_level_engine_recovers_a_zero_inflated_world(registry, layout, tasks, zero_inflated, tmp_path):
    hurdle, log = zero_inflated["hurdle"], zero_inflated["log"]
    fresh, soma, y, label = zero_inflated_rows(1500, 78, prefix="f")
    index = np.array([LABELS.index(l) for l in label])
    log_terms = np.log(TABLE[:, index]).T - 0.5 * (((soma[:, None, :] - CENTRES[None]) / SD_SOMA) ** 2).sum(axis=2)
    truth = class_posterior(log_terms) @ (PRESENT * LEVELS)                 # E[log-length | soma, label]
    predicted = mods.engine_predict(hurdle, tasks["T1"], fresh, level=0.9, median=True)
    plain = mods.engine_predict(log, tasks["T1"], fresh, level=0.9, median=True)
    error = lambda p: np.sqrt(((p - truth) ** 2).mean())
    # with as many components as classes both representations find the conditional MEAN ...
    assert error(predicted["point"]) < 0.1 * truth.std() and error(predicted["point"]) < 1.1 * error(plain["point"])
    # ... but only the presence/level form states the zeros: an exact-zero median where zero is the likelier
    # outcome, a lower absolute error for it, and intervals that hold for zeros and non-zeros alike
    inside = (y >= predicted["lower"]) & (y <= predicted["upper"])
    assert 0.90 <= inside.mean() < 0.97 and inside[y > 0].mean() > 0.85 and inside[y == 0].mean() > 0.95
    assert abs((predicted["median"] == 0).mean() - (y == 0).mean()) < 0.05 and (plain["median"] == 0).mean() == 0
    absolute = lambda p: np.abs(p - y).mean()
    assert absolute(predicted["median"]) < 0.9 * absolute(predicted["point"])
    assert absolute(predicted["median"]) < 0.9 * absolute(plain["median"])
    # zeros in a measured profile are evidence too (soma from profile), also through the request API
    t2 = mods.engine_predict(hurdle, tasks["T2"], fresh)
    assert np.sqrt(((t2["point"] - CENTRES[nearest_class(soma)]) ** 2).sum(axis=1)).mean() < 40.0
    module = mods.JointMixtureModule(registry, layout, states={n: hurdle for n in tasks}, train=zero_inflated["rows"],
                                     tasks=tasks, evidence=evidence_for(tasks, "exploratory_supported"))
    r = core.Reconstructor(registry)
    r.register(module)
    answer = r.reconstruct(request(registry, {"regional_axon_length": lengths(y[0]), "source_label": label[0]},
                                   ["soma_position"]))["soma_position"]
    np.testing.assert_allclose(answer.estimate.point, t2["point"][0], rtol=1e-9)
    profile = r.reconstruct(request(registry, {"soma_position": soma[0], "source_label": label[0]},
                                    ["regional_axon_length"]))["regional_axon_length"]
    np.testing.assert_allclose(profile.estimate.point, predicted["point"][0], rtol=1e-9)
    np.testing.assert_allclose(profile.estimate.original_units["upper"], 10 ** predicted["upper"][0] - 1, rtol=1e-9)
    # state identity
    loaded = mods.JointState.load(hurdle.save(tmp_path / "hurdle.npz"), hurdle.state_hash())
    assert loaded.state_hash() == hurdle.state_hash() != log.state_hash() and loaded.encoding == hurdle.encoding
    np.testing.assert_array_equal(mods.engine_predict(loaded, tasks["T1"], fresh.take(np.arange(20)))["point"],
                                  predicted["point"][:20])


def test_joint_module_recovers_a_planted_conditional(registry, layout, tasks, planted):
    """In a world that IS a latent-class mixture, the fitted conditionals approach the true ones."""
    state = planted["state"]
    assert all(b["converged"] for b in planted["diagnostics"])
    fresh, soma, y, label = planted_rows(1500, 77, prefix="f")
    t1 = mods.engine_predict(state, tasks["T1"], fresh, level=0.9, median=True)
    truth1 = true_profile_given_soma_label(soma, label)
    error1 = np.sqrt(((t1["point"] - truth1) ** 2).mean())
    assert error1 < 0.15 * truth1.std() and np.corrcoef(t1["point"].ravel(), truth1.ravel())[0, 1] > 0.99
    t2 = mods.engine_predict(state, tasks["T2"], fresh, level=0.9)
    truth2 = true_soma_given_profile_label(y, label)
    assert np.sqrt(((t2["point"] - truth2) ** 2).sum(axis=1)).mean() < 25.0        # um; class centres are ~1 mm apart
    t3 = mods.engine_predict(state, tasks["T3"], fresh)
    accuracy = np.mean(np.array(t3["classes"])[t3["proba"].argmax(axis=1)] == label)
    best_possible = np.mean(np.array(LABELS)[TABLE[nearest_class(soma)].argmax(axis=1)] == label)
    assert accuracy > best_possible - 0.03          # the label is only 80% (or 45%) determined by the class
    # the engine beats the simple rules here, and its intervals cover as stated
    bank = mods.TaskBaselines(tasks["T1"], planted["rows"])
    mse = lambda p: ((p - y) ** 2).mean()
    assert mse(t1["point"]) < 0.8 * min(mse(bank.predict(name, fresh)) for name in ("MEAN", "IDENTITY", "COPY"))
    coverage1 = ((y >= t1["lower"]) & (y <= t1["upper"])).mean()
    coverage2 = ((soma >= t2["lower"]) & (soma <= t2["upper"])).mean()
    assert abs(coverage1 - 0.9) < 0.03 and abs(coverage2 - 0.9) < 0.03
    # and the same answer comes out of the request API
    r, _ = reconstructor(registry, layout, tasks, planted)
    answer = r.reconstruct(request(registry, {"regional_axon_length": lengths(y[0]), "source_label": label[0]},
                                   ["soma_position"]))["soma_position"]
    np.testing.assert_allclose(answer.estimate.point, t2["point"][0], rtol=1e-9)
    assert answer.estimate.scale == "linear position; um" and answer.estimate.transform == "identity"
    assert answer.estimate.frame == core.JIAO_FRAME and answer.estimate.original_units is None


def test_no_signal_world_returns_the_baseline_flag(null_world):
    """Where soma, label and profile are independent, the benchmark says so and the API falls back."""
    world, metrics = null_world
    evidence = core.read_json(Path(world.out) / "evidence.json")
    assert evidence["exploratory"] and not evidence["confirmatory"]
    for name in ("T1", "T2", "T3", "T4"):
        assert evidence["tasks"][name]["status"] == "no_better_than_baseline"
        assert metrics["tasks"][name]["verdict"]["all_comparators_clear"] != "every"
        assert evidence["tasks"][name]["min_relative_gain"] == 0.02
    t1 = metrics["tasks"]["T1"]["methods"]
    assert abs(t1["ENGINE"]["mse"]["value"] - t1["MEAN"]["mse"]["value"]) < 0.05 * t1["MEAN"]["mse"]["value"]
    assert t1["COPY"]["mse"]["value"] > 1.5 * t1["MEAN"]["mse"]["value"]            # copying noise doubles the error
    r, ctx = bench.build_reconstructor(world, fold=0)
    i = int(np.flatnonzero(ctx.outer_of == 0)[0])
    make = lambda: core.ReconstructionRequest.from_values(
        ctx.registry, namespace="jiao", unit_id=str(ctx.table.neuron[i]), group_id=str(ctx.table.sample[i]),
        observed={"soma_position": ctx.table.soma[i], "source_label": str(ctx.table.label[i])},
        requested=["regional_axon_length"])
    answer = r.reconstruct(make())["regional_axon_length"]
    assert answer.status == "inferred" and answer.baseline_level and "no_better_than_baseline" in answer.codes()
    assert answer.module.startswith("baseline:") and answer.validation == "no_better_than_baseline"
    assert answer.benchmark["metrics_sha256"] == core.sha_file(Path(world.out) / "metrics.json")
    assert "group_in_training_state" not in answer.codes()                         # fold-0 states never saw this sample
    strict, _ = bench.build_reconstructor(world, fold=0, policy=core.Policy(on_no_better_than_baseline="refuse"), ctx=ctx)
    refused = strict.reconstruct(make())["regional_axon_length"]
    assert refused.status == "refused" and refused.refusal.code == "no_better_than_baseline"
    serving, _ = bench.build_reconstructor(world, ctx=ctx)                          # all-training state: unit is inside
    assert serving.reconstruct(make())["regional_axon_length"].refusal.code == "leakage_guard"


def test_benchmark_pipeline_runs_once_and_is_bound_to_its_freeze(planted_world, tmp_path, monkeypatch):
    world, metrics = planted_world
    out = Path(world.out)
    for name in ("pre_run_declaration.json", "run_started.json", "predictions_complete.json", "run_complete.json",
                 "metrics.json", "evidence.json", "selection.json", "fits.json", "serving.json", "api_demo.json",
                 "results.md"):
        assert (out / name).exists(), name
    assert all(bench.status(world)[name] for name in ("run_complete.json", "serving.json", "api_demo.json"))
    assert bench.status(world)["freeze_matches_current_files"] is True
    assert metrics["code_files_changed_since_freeze"] == []
    declaration = core.read_json(out / "pre_run_declaration.json")
    assert declaration["exploratory"] and not declaration["confirmatory_claim"]
    assert set(declaration["tasks"]) == {"T1", "T2", "T3", "T4"} and declaration["splits_sha256"]
    assert declaration["code_sha256"]["analysis/brm_joint.py"] == core.BRM_JOINT_SHA256
    ctx = bench.build_context(world)
    selection = core.read_json(out / "selection.json")
    for name, task in ctx.tasks.items():
        with np.load(out / "predictions" / f"{name}.npz") as data:
            assert list(data["neuron"]) == list(ctx.table.neuron) and (data["fold"] == ctx.outer_of).all()
            assert np.isfinite(data["ENGINE"]).all()                               # every neuron predicted once
            if task.kind != "label":
                assert (data["ENGINE_lower"] <= data["ENGINE_upper"]).all()
            if task.kind == "profile":
                assert all((data[m] >= 0).all() for m in ("ENGINE", "ENGINE_MEDIAN", "RIDGE", "ENGINE_lower"))
        for fold in range(ctx.config["outer_folds"]):                              # by states that never saw its sample
            used = selection["engine_used"][name][str(fold)]
            state = mods.JointState.load(out / "states" / bench.state_file(f"outer{fold}", used["representation"],
                                                                           used["components"]), used["state_hash"])
            assert state.settings["representation"] == used["representation"]
            assert state.state_hash() == used["state_hash"]
            assert not set(state.groups) & set(ctx.rows.group[ctx.outer_of == fold])
            assert not state.units & set(ctx.rows.unit[ctx.outer_of == fold])
            ctx.partition.require_training(state.groups, "state")
    t1 = metrics["tasks"]["T1"]
    assert t1["methods"]["ENGINE"]["mse"]["value"] < t1["methods"]["MEAN"]["mse"]["value"]    # planted signal is found
    assert t1["comparisons"]["mse"]["MEAN"]["engine_better_interval"]
    assert t1["samples"] == dict(world.expected)["train_samples"] and t1["neurons"] == len(ctx.table)
    assert set(t1["by_source_label"]) == set(ctx.classes)
    crh = t1["by_source_label"]["Crh"]
    assert crh["samples"] == 1 and crh["descriptive_only"] and "comparisons" not in crh
    assert "comparisons" in t1["by_source_label"]["Pomc"]
    assert metrics["tasks"]["T3"]["unseen_label"]["samples"] == 1                   # the single Crh sample
    assert t1["rows_with_label_unseen_in_training"] == crh["neurons"]
    assert set(t1["coverage"]) == {"ENGINE", "MEAN", "IDENTITY"} and 0 <= t1["coverage"]["ENGINE"]["empirical"] <= 1
    # runs once
    with pytest.raises(core.ReconError, match="already been started"):
        bench.run(world)
    # the finished run binds its files; a copy of the folder verifies as well as the original
    complete = core.read_json(out / "run_complete.json")
    assert set(complete["sha256"]) == set(bench.BOUND_FILES) and complete["completed_by"] == "run"
    assert complete["sha256"]["metrics.json"] == core.sha_file(out / "metrics.json")
    # the completion record HOLDS the record of the all-training states; serving.json is a copy of it
    bound = complete[bench.SERVING_KEY]
    assert set(bound["tasks"]) == {"T1", "T2", "T3", "T4"} and all(e["state_hash"] for e in bound["tasks"].values())
    assert core.read_json(out / "serving.json") == dict(bound, run_complete_sha256=core.sha_file(out / "run_complete.json"))
    for name, entry in bound["tasks"].items():                 # each fitted on every training sample
        file = out / "states" / bench.state_file("serving", entry["representation"], entry["components"])
        state = mods.JointState.load(file, entry["state_hash"])
        assert set(state.groups) == set(ctx.rows.group) and state.groups_sha256 == entry["training_groups_sha256"]
        assert [entry["representation"], entry["components"]] == [selection["serving"][name][k]
                                                                  for k in ("representation", "components")]
    copy = tmp_path / "copy"
    shutil.copytree(Path(world.manifest).parent, copy)
    twin = core.synthetic_world(copy)
    bench.verify_run(twin)
    # scoring is reproducible from the saved predictions; a re-score never overwrites the run's own metrics
    rescored = bench.score(twin, log=lambda *a: None)
    assert core.canonical(rescored["tasks"]["T1"]["methods"]) == core.canonical(metrics["tasks"]["T1"]["methods"])
    assert rescored["code_files_changed_since_freeze"] == [] and rescored["rescore"]["number"] == 1
    assert (copy / "out" / "rescores" / "001" / "metrics.json").exists()
    assert core.sha_file(copy / "out" / "metrics.json") == core.sha_file(out / "metrics.json")
    assert core.sha_file(copy / "out" / "evidence.json") == core.sha_file(out / "evidence.json")
    bench.verify_run(twin)
    with pytest.raises(core.ReconError, match="already complete"):
        bench.finish(twin, how="again")
    # scoring or reporting with other code is possible; it is said in the re-score and in the FIRST lines of the report
    real_hashes = bench.code_hashes
    monkeypatch.setattr(bench, "code_hashes", lambda: dict(real_hashes(), **{"analysis/synarch_recon_benchmark.py": "x"}))
    again = bench.score(twin, log=lambda *a: None)
    assert again["code_files_changed_since_freeze"] == ["analysis/synarch_recon_benchmark.py"]
    assert again["rescore"]["number"] == 2 and (copy / "out" / "rescores" / "002" / "metrics.json").exists()
    report = bench.render_report(twin).splitlines()
    assert "differs from the frozen code" in report[2] and "synarch_recon_benchmark.py" in report[2]
    monkeypatch.undo()
    assert "differs from the frozen code" not in bench.render_report(world)
    text = bench.render_report(world)
    assert "Multiplicity: 4 tasks" in text and "minimum gain is 2%" in text and "K = 1 was chosen in" in text
    assert "The verdict of a task is that last comparison" in text
    assert f"SHA-256 `{core.sha_file(out / 'run_complete.json')}`" in text        # the report names its anchor
    with np.load(copy / "out" / "predictions" / "T1.npz") as data:
        arrays = {k: data[k] for k in data.files}
    arrays["ENGINE"] = arrays["ENGINE"] * 0.5
    np.savez_compressed(copy / "out" / "predictions" / "T1.npz", **arrays)
    with pytest.raises(core.FrozenSourceChanged):
        bench.score(twin, log=lambda *a: None)
    # the freeze binds configuration, inputs and code
    bench.verify_freeze(world)
    meta = core.read_json(copy / "world.json")
    meta["config"]["engine"]["k_grid"] = [1, 2]
    core.write_json(copy / "world.json", meta)
    with pytest.raises(core.FrozenSourceChanged, match="config"):
        bench.verify_freeze(twin)
    with pytest.raises(FileExistsError):
        bench.freeze(world, "again", log=lambda *a: None)
    with pytest.raises(core.ReconError):
        bench.freeze(bench.make_synthetic_world(tmp_path / "unreviewed", "null", size="tiny", seed=5), "  ")


def test_a_broken_start_can_be_set_aside_only_before_outcomes_were_used(bare_world, planted_world, tmp_path):
    world = bare_world
    out = Path(world.out)
    with pytest.raises(core.ReconError):
        bench.abort(world, "nothing was frozen or started")
    bench.freeze(world, "synthetic", log=lambda *a: None)
    unfrozen = bench.abort(world, "a late review comment: freeze withdrawn before any run", log=lambda *a: None)
    record = core.read_json(unfrozen / "abort.json")
    assert not record["run_had_started"] and record["outcome_files_present"] == []
    assert record["pre_run_declaration_sha256"] == core.sha_file(unfrozen / "pre_run_declaration.json")
    # a freeze that was withdrawn before any run binds nothing: the design may still change
    assert bench.freeze(world, "synthetic, second freeze", log=lambda *a: None)["supersedes"] is None
    core.write_json(out / "run_started.json", {"started_utc": "simulated crash"}, exclusive=True)
    with pytest.raises(core.ReconError, match="already been started"):
        bench.run(world, log=lambda *a: None)
    with pytest.raises(core.ReconError):
        bench.abort(world, " ")
    folder = bench.abort(world, "simulated crash before any prediction", log=lambda *a: None)
    assert (folder / "pre_run_declaration.json").exists() and (folder / "run_started.json").exists()
    assert not (out / "pre_run_declaration.json").exists() and not (out / "run_started.json").exists()
    assert core.read_json(folder / "abort.json")["run_had_started"]
    with pytest.raises(core.ReconError, match="not frozen"):
        bench.run(world, log=lambda *a: None)
    # after a STARTED run was set aside the next freeze must be the same benchmark: other tasks, another grid
    # or another salt (other folds) are refused; only code and protocol text may differ
    meta_path = Path(world.manifest).parent / "world.json"
    original = meta_path.read_text()
    for change in (lambda c: c.__setitem__("tasks", ["T1", "T2"]),
                   lambda c: c["engine"].__setitem__("k_grid", [1, 2]),
                   lambda c: c.__setitem__("salt", "another-salt"),
                   lambda c: c.__setitem__("min_relative_gain", 0.0)):
        meta = json.loads(original)
        change(meta["config"])
        core.write_json(meta_path, meta, pretty=True)
        with pytest.raises(core.ReconError, match="must keep its design"):
            bench.freeze(core.synthetic_world(meta_path.parent), "another design", log=lambda *a: None)
        assert not (out / "pre_run_declaration.json").exists()
    meta_path.write_text(original)
    declaration = bench.freeze(world, "synthetic, after the abort", log=lambda *a: None)
    assert [a["reason"] for a in declaration["aborted_starts_before_this_freeze"]] == [
        "a late review comment: freeze withdrawn before any run", "simulated crash before any prediction"]
    previous = declaration["supersedes"]
    assert previous["design_differences"] == [] and previous["changed_files"] == []
    assert previous["pre_run_declaration_sha256"] == core.sha_file(folder / "pre_run_declaration.json")
    metrics = bench.run(world, log=lambda *a: None)
    assert set(metrics["tasks"]) == {"T1", "T2", "T3", "T4"} and "set aside" in (out / "results.md").read_text()
    for finished in (world, planted_world[0]):              # results exist: they stand
        with pytest.raises(core.ReconError, match="cannot be set aside"):
            bench.abort(finished, "I did not like the result")
    # A start that got as far as writing predictions or inner-fold scores has USED outcomes: it cannot be set
    # aside, even though no metric file exists yet; it is finished from its saved predictions instead.
    second = bench.make_synthetic_world(tmp_path / "second", "planted", size="tiny", seed=4)
    out2 = Path(second.out)
    bench.freeze(second, "synthetic", log=lambda *a: None)
    core.write_json(out2 / "run_started.json", {"started_utc": "simulated"}, exclusive=True)
    ctx = bench.build_context(second)
    result = bench.execute(ctx, log=lambda *a: None)
    manifest = bench.save_predictions(ctx, result)           # the 'crash' happens right after this line
    assert (out2 / "fits.json").exists() and not (out2 / "predictions_complete.json").exists()
    with pytest.raises(core.ReconError, match="cannot be set aside"):
        bench.abort(second, "simulated crash after the predictions were written")
    assert (out2 / "pre_run_declaration.json").exists() and (out2 / "fits.json").exists()
    with pytest.raises(core.ReconError, match="not saved completely"):
        bench.finish(second, how="test")
    core.write_json(out2 / "predictions_complete.json", {"utc": "simulated", "predictions": manifest})
    finished = bench.finish(second, how="score command, after an interrupted run", log=lambda *a: None)
    assert set(finished["tasks"]) == {"T1", "T2", "T3", "T4"} and (out2 / "run_complete.json").exists()
    assert core.read_json(out2 / "run_complete.json")["completed_by"].startswith("score command")
    bench.verify_run(second)


def test_a_failed_all_training_fit_does_not_touch_the_benchmark(bare_world, monkeypatch):
    world = bare_world
    total = dict(world.expected)["train_neurons"]
    real_fit = mods.fit_joint

    def no_serving(rows, layout, **kwargs):
        if len(rows) == total:                       # only the all-training fit uses every training neuron
            raise FloatingPointError("planted failure of the all-training fit")
        return real_fit(rows, layout, **kwargs)

    monkeypatch.setattr(mods, "fit_joint", no_serving)
    bench.freeze(world, "synthetic", log=lambda *a: None)
    metrics = bench.run(world, log=lambda *a: None)
    out = Path(world.out)
    assert (out / "run_complete.json").exists() and set(metrics["tasks"]) == {"T1", "T2", "T3", "T4"}
    serving = core.read_json(out / "serving.json")
    assert all("planted failure" in entry["failed"] for entry in serving["tasks"].values())
    assert core.read_json(out / "run_complete.json")[bench.SERVING_KEY]["tasks"] == serving["tasks"]
    with pytest.raises(core.ReconError, match="all-training fit"):
        bench.build_reconstructor(world)
    held_out, _ = bench.build_reconstructor(world, fold=0)          # out-of-fold states are unaffected
    assert len(held_out.describe()) == 2
    # The record of the all-training states is part of the completion record. By hand they can be RESTORED, not
    # changed: a fit that failed in the run stays failed for that run, and no file of the run is rewritten.
    monkeypatch.undo()
    before = {name: core.sha_file(out / name) for name in ("run_complete.json", "serving.json", "metrics.json")}
    kept = bench.fit_serving(world, log=lambda *a: None)
    assert all("planted failure" in entry["failed"] for entry in kept["tasks"].values())
    assert before == {name: core.sha_file(out / name) for name in before} and not list((out / "states").glob("serving_*"))
    with pytest.raises(core.ReconError, match="all-training fit"):
        bench.build_reconstructor(world)
    bench.verify_run(world)
    with pytest.raises(core.ReconError):
        bench.fit_serving(bench.make_synthetic_world(Path(world.manifest).parent.parent / "other", "null", seed=9))
    # any OTHER failure of the step is recorded too, for every task it left without a state: the step never
    # stops the benchmark
    ctx = bench.build_context(world)
    chosen = {"T1": {"representation": "log", "components": 1}, "T2": {"representation": "log", "components": 2}}

    def broken_pool(*args, **kwargs):
        raise RuntimeError("the worker pool died")

    monkeypatch.setattr(bench, "run_jobs", broken_pool)
    record = bench.all_training_record(ctx, chosen, log=lambda *a: None)
    assert record["fits"] == [] and set(record["tasks"]) == {"T1", "T2"}
    for entry in record["tasks"].values():
        assert entry["state_hash"] is None and "did not finish" in entry["failed"] and "pool died" in entry["failed"]


def test_failed_fit_is_recorded_and_the_next_setting_is_used(bare_world, monkeypatch):
    ctx = bench.build_context(bare_world)
    real_fit = mods.fit_joint

    def flaky(rows, layout, *, components, seed, budget, representation="log"):
        if components == 4:
            raise FloatingPointError("planted failure")
        return real_fit(rows, layout, components=components, seed=seed, budget=budget, representation=representation)

    monkeypatch.setattr(mods, "fit_joint", flaky)
    inner = bench.job_engine_inner(ctx, 0, 0, "hurdle", 4)
    assert "failed" in inner and "planted failure" in inner["failed"]
    records = [bench.job_engine_inner(ctx, o, i, r, k) for o in range(5) for i in range(3)
               for r, k in bench.engine_grid(ctx.config)]
    selection = bench.select_engine(ctx, records)
    rankings = [selection[t][o]["ranking"] for t in selection for o in selection[t]]
    assert all(len(ranking) == 4 and all(k != 4 for _, k in ranking) for ranking in rankings)
    outer = bench.job_engine_outer(ctx, 0, "log", 4, ["T1"])
    assert "failed" in outer and outer["tasks"] == ["T1"]


def test_the_real_cohort_cannot_be_run_or_opened_before_its_freeze(tmp_path):
    """Pointed at an empty output folder, the real cohort refuses: no declaration, no outcome is read."""
    world = dataclasses.replace(core.real_world(), out=tmp_path / "empty")
    assert not world.synthetic
    with pytest.raises(core.ReconError, match="not frozen"):
        bench.run(world)
    with pytest.raises(core.ReconError, match="not frozen"):
        bench.build_context(world)
    with pytest.raises(core.ReconError):
        bench.report(world)
    assert not (tmp_path / "empty" / "run_started.json").exists()
    with pytest.raises(core.ReconError):
        bench.freeze(world, "")


# ----------------------------------------------------------------------------------------------
# the axon adapter
# ----------------------------------------------------------------------------------------------


def test_axon_adapter_is_blocked_and_rejects_a_measured_profile(registry, layout, tasks, planted):
    r, _ = reconstructor(registry, layout, tasks, planted)
    _, soma, y, label = planted_rows(1, 89, prefix="f")
    permitted = {"soma_position": soma[0], "source_label": label[0], "soma_label_fine": "ARH"}
    blocked = r.reconstruct(request(registry, permitted, ["axon_arbor"]))["axon_arbor"]
    assert blocked.status == "refused" and blocked.refusal.code == "not_validated"
    assert blocked.validation == "blocked: independent validation pending" == blocked.refusal.reason
    assert blocked.module == "axon_arbor_adapter" and blocked.estimate is None
    for profile in (lengths(y[0]), lengths(y[0][:3], INSIDE)):                      # whole or partial: both rejected
        observed = dict(permitted, regional_axon_length=profile)
        answer = r.reconstruct(request(registry, observed, ["axon_arbor"]))["axon_arbor"]
        assert answer.status == "refused" and answer.refusal.code == "forbidden_input"
        assert "never an input to a hidden-axon forecast" in answer.refusal.reason
    adapter = r.module("axon_arbor_adapter")
    assert adapter.STATUS == "blocked: independent validation pending"
    assert adapter.PERMITTED == core.HIDDEN_AXON_INPUTS == ("soma_position", "source_label", "soma_label_fine",
                                                           "hemisphere", "reconstruction_type")
    measured_profile = request(registry, {"regional_axon_length": lengths(y[0])}, []).observed[0]
    with pytest.raises(core.ForbiddenInputError, match="never an input to a hidden-axon forecast"):
        adapter.check_inputs("fresh", {"regional_axon_length": measured_profile})
    with pytest.raises(core.ForbiddenInputError, match="never be a predictor"):
        adapter.check_inputs("fresh", {"hemisphere": measured("hemisphere", "Left"),
                                       "projectome_type": measured("projectome_type", "t1")})
    assert adapter.check_inputs("fresh", {"hemisphere": measured("hemisphere", "Left"),
                                          "reconstruction_type": measured("reconstruction_type", "Axon_only")})
    listing = {m["name"]: m for m in r.describe()}["axon_arbor_adapter"]
    assert listing["produces"] == ["axon_arbor"] and listing["conditionals"][0]["validation"]["status"] == adapter.STATUS
    # NOTHING unblocks it: not a backend, and not an evidence object (this layer holds exploratory evidence only)
    calls = []
    backend = lambda query: calls.append(query) or "a tree"                         # noqa: E731
    for extra in (dict(backend=backend),
                  dict(backend=backend, validation=core.ValidationEvidence("exploratory_supported", benchmark=POINTER))):
        layer = core.Reconstructor(registry)
        attached = layer.register(mods.AxonArborAdapter(registry, **extra))
        still = layer.reconstruct(request(registry, permitted, ["axon_arbor"]))["axon_arbor"]
        assert still.status == "refused" and still.refusal.code == "not_validated" and still.estimate is None
        assert still.validation == adapter.STATUS == attached.evidence_for("axon_arbor_from_permitted_inputs").status
        with pytest.raises(core.ReconError, match="blocked"):
            attached.estimate(None, 0.9)
    assert calls == []                                                              # no backend was ever called
    with pytest.raises(ValueError):                 # and no evidence object of this layer can claim confirmation
        core.ValidationEvidence("exploratory_supported", benchmark=POINTER, confirmatory=True)
    source = (ROOT / "analysis/synarch_recon_modules.py").read_text() + (ROOT / "analysis/synarch_recon.py").read_text()
    assert "import axon_engine" not in source and "axon_engine.py" not in source


# ==============================================================================================
# Regression: review blockers K1 K2 K3
# ==============================================================================================
#
# One focused test per item of the three acceptance criteria. Every test below FAILS on a copy of the code in
# which only that criterion's fix is switched off, and PASSES on the code as it is; both runs, and the
# switching-off edits, are saved under data/interim/synarch_recon/fix_evidence/round2/ (for the code as it is
# now; fix_evidence/ itself holds the same proof for the code of the first round).
#
#   K1  provenance and transitive ancestry are validated for EVERY inferred input of EVERY request, including
#       inferred soma positions and source labels; hidden measured-axon ancestry and producing states that
#       contain the query's brain sample are refused.            K1a-K1f; second round: K1g, K1h, K1i
#   K2  every record states its scale and frame, and its value lies in the frame it states; a value reaches
#       the modelling scale by exactly one conversion; a mismatch is refused.
#                                                                  K2a-K2e; second round: K2f, K2g
#   K3  the fingerprint of a saved state binds the support ranges and every other answer-affecting item, the
#       expected fingerprint of every served state is behind the completion record; tampering is rejected.
#                                                                  K3a-K3d and the run files; second round: K3e
#
# Everything here is synthetic.


class _Chain:
    """What the K1 tests share (see the `chain` fixture)."""


@pytest.fixture(scope="module")
def chain(registry, layout, tasks, planted):
    """Two states for one query neuron. State A was fitted WITH the neuron's brain sample; state B without it.
    Four reconstruction layers share one ledger of issued values:
        held_out    state B, default policy, the declared tasks: the layer that would make the forecast;
        lenient_B   state B, undeclared conditionals answered and flagged;
        with_A      state A, undeclared conditionals answered and flagged;
        chained     state B, plus two test-only conditionals DECLARED as profile-input tasks (soma from the
                    measured profile; label from a soma), so that an honest two-step chain exists;
        baseline_B  state B, every task 'no better than baseline': baseline-level answers.
    `adapter` is the axon adapter registered in held_out; state A is declared to it (and to the adapter of
    lenient_B), so that it can be resolved there."""
    c = _Chain()
    rows = planted["rows"]
    c.group = str(rows.group[0])                                   # the query's brain sample: in A, not in B
    c.unit = f"{c.group}_99999"                                    # a NEW neuron of that sample
    rows_B = rows.take(np.flatnonzero(rows.group != c.group))
    c.state_A = planted["state"]
    c.state_B, _ = mods.fit_joint(rows_B, layout, components=4, seed=5,
                                  budget={"n_init": 2, "max_iter": 120, "tol": 1e-6})
    assert c.group in c.state_A.groups and c.group not in c.state_B.groups
    declared = dict(tasks)
    declared["T2b"] = mods.make_task(layout, registry, task="T2b", title="soma_from_profile", target="soma_position",
                                     evidence=(("regional_axon_length", None),), declared_measured_axon_input=True)
    declared["T3b"] = mods.make_task(layout, registry, task="T3b", title="label_from_soma_declared",
                                     target="source_label", evidence=(("soma_position", None),),
                                     declared_measured_axon_input=True)
    c.ledger = core.ProvenanceLedger()
    lenient = core.Policy(on_not_validated="engine_flagged")

    def layer(state, train, policy=None, task_set=tasks, status="exploratory_supported"):
        module = mods.JointMixtureModule(registry, layout, states={n: state for n in task_set}, train=train,
                                         tasks=task_set, evidence=evidence_for(task_set, status))
        r = core.Reconstructor(registry, policy, c.ledger)
        r.register(module)
        r.register(mods.AxonArborAdapter(registry))
        return r

    c.held_out, c.lenient_B = layer(c.state_B, rows_B), layer(c.state_B, rows_B, lenient)
    c.with_A = layer(c.state_A, rows, lenient)
    c.chained = layer(c.state_B, rows_B, task_set=declared)
    c.baseline_B = layer(c.state_B, rows_B, status="no_better_than_baseline")
    c.adapter = c.held_out.module("axon_arbor_adapter")
    c.adapter.register_state(c.state_A)
    c.lenient_B.module("axon_arbor_adapter").register_state(c.state_A)       # so that state A is known there too
    _, soma, y, label = planted_rows(1, 83, prefix="f")
    c.soma, c.label, c.profile = tuple(float(v) for v in soma[0]), str(label[0]), lengths(y[0])
    c.ask = lambda observed, requested, unit=None, group=None: request(                      # noqa: E731
        registry, observed, requested, unit=unit or c.unit, group=c.group if group is None else group)
    c.m_soma, c.m_label = c.ask({"soma_position": c.soma, "source_label": c.label}, []).observed

    def handmade(variable, value, evidence, *, state_hash, groups_hash, coordinates=(), module="joint_mixture"):
        """An inferred record of the query neuron that nobody issued: built by hand, with a one-step lineage."""
        spec = registry.variable("jiao", variable)
        key = ("jiao", "neuron", c.unit)
        uncertainty = core.Uncertainty("inferred", None, None, None, None, module, state_hash, groups_hash, "unvalidated")
        return core.Observation("synthetic", "jiao", "neuron", c.unit, c.m_soma.group, variable, "inferred", value,
                                core.scale_statement(spec), tuple(coordinates), spec.frame, uncertainty,
                                core.SourceRef("reconstruct:" + module, state_hash, c.unit), derived_from=tuple(evidence),
                                transform="identity",
                                lineage=(core.LineageStep(variable, module, state_hash, groups_hash, tuple(evidence),
                                                          key, c.group),))

    c.handmade = handmade
    return c


def test_K1a_inferred_source_label_descending_from_the_measured_profile_is_refused(registry, chain):
    """[K1 a] The other agent's counterexample
    (data/interim/synarch_recon_codex_review/lineage_counterexample.json): an INFERRED source_label whose direct
    evidence is the neuron's MEASURED regional profile, produced by a state that is not registered. The reviewed
    adapter accepted it, because 'source_label' is a permitted name and the checks ran for axon-family names
    only. It must be refused: as built by hand, with its state registered, with a state that holds the query's
    brain sample, and as the layer itself issues it (task T3: label from soma + measured profile)."""
    c = chain
    from_profile = (("regional_axon_length", "measured"),)
    forged = c.handmade("source_label", "A", from_profile, state_hash="unregistered-state", groups_hash="x")
    assert forged.status == "inferred" and forged.derived_from == from_profile
    with pytest.raises(core.ForbiddenInputError):                              # 1. exactly the counterexample
        c.adapter.check_inputs(c.group, {"source_label": forged})
    registered = mods.AxonArborAdapter(registry, ledger=c.ledger)              # 2. its state registered, the profile
    registered.register_state("unregistered-state", ["another_sample"],        #    even declared as an input of it
                              ("soma_position", "source_label", "regional_axon_length"))
    with pytest.raises(core.ForbiddenInputError):
        registered.check_inputs(c.group, {"source_label": forged})
    holding = mods.AxonArborAdapter(registry, ledger=c.ledger)                 # 3. a state that holds the query's sample
    holding.register_state("unregistered-state", [c.group, "another_sample"])
    with pytest.raises(core.ForbiddenInputError):
        holding.check_inputs(c.group, {"source_label": forged})
    # 4. the honest case: the layer itself infers the label from soma + MEASURED profile (a valid T3 answer of a
    #    registered state that excludes the brain sample). As an input of the axon forecast it is refused for
    #    what it descends from.
    t3 = c.held_out.reconstruct(c.ask({"soma_position": c.soma, "regional_axon_length": c.profile}, ["source_label"]))
    assert t3["source_label"].status == "inferred" and t3["source_label"].coverage_validated
    label_record = t3["source_label"].to_observation(t3.request)
    assert dict(label_record.derived_from)["regional_axon_length"] == "measured"
    with pytest.raises(core.ForbiddenInputError, match="descends from a measured regional_axon_length"):
        c.adapter.check_inputs(c.group, {"source_label": label_record})
    with pytest.raises(core.ForbiddenInputError, match="descends from a measured regional_axon_length"):
        c.adapter.check_inputs(c.group, {"soma_position": c.m_soma, "source_label": label_record})
    refused = c.held_out.reconstruct(c.ask({"soma_position": c.soma, "source_label": label_record},
                                           ["axon_arbor"]))["axon_arbor"]
    assert refused.status == "refused" and refused.refusal.code == "forbidden_input"
    assert "descends from a measured regional_axon_length" in refused.refusal.reason


def test_K1b_inferred_soma_position_descending_from_the_measured_profile_is_refused(registry, chain):
    """[K1 b] The same through the other permitted name: a soma position inferred from the measured profile
    (task T2) is an inferred value of 'soma_position', and it knows the hidden axon."""
    c = chain
    t2 = c.held_out.reconstruct(c.ask({"regional_axon_length": c.profile, "source_label": c.label}, ["soma_position"]))
    assert t2["soma_position"].status == "inferred" and t2["soma_position"].coverage_validated
    soma_record = t2["soma_position"].to_observation(t2.request)
    with pytest.raises(core.ForbiddenInputError, match="descends from a measured regional_axon_length"):
        c.adapter.check_inputs(c.group, {"soma_position": soma_record})
    with pytest.raises(core.ForbiddenInputError, match="descends from a measured regional_axon_length"):
        c.adapter.check_inputs(c.group, {"soma_position": soma_record, "source_label": c.m_label})
    refused = c.held_out.reconstruct(c.ask({"soma_position": soma_record, "source_label": c.label},
                                           ["axon_arbor"]))["axon_arbor"]
    assert refused.status == "refused" and refused.refusal.code == "forbidden_input"
    forged = c.handmade("soma_position", c.soma, (("regional_axon_length", "measured"),), coordinates=core.SOMA_AXES,
                        state_hash=c.state_B.state_hash(), groups_hash=c.state_B.groups_sha256)
    with pytest.raises(core.ForbiddenInputError):
        c.adapter.check_inputs(c.group, {"soma_position": forged, "source_label": c.m_label})
    # the general layer's own hidden-axon forecast (T1) refuses the same input
    t1 = c.held_out.reconstruct(c.ask({"soma_position": soma_record, "source_label": c.label},
                                      ["regional_axon_length"]))["regional_axon_length"]
    assert t1.status == "refused" and t1.refusal.code == "forbidden_input" and t1.estimate is None


def test_K1c_transitive_chain_profile_to_soma_to_label_is_refused(registry, chain):
    """[K1 c] measured profile -> inferred soma -> inferred source label -> adapter. The label's own evidence is
    only an inferred soma: nothing in its first step names the profile. The ancestry is followed transitively,
    in what the layer recorded when it issued each value, and it cannot be cut off."""
    c = chain
    hop1 = c.chained.reconstruct(c.ask({"regional_axon_length": c.profile}, ["soma_position"]))
    soma_record = hop1["soma_position"].to_observation(hop1.request)
    hop2 = c.chained.reconstruct(c.ask({"soma_position": soma_record}, ["source_label"]))
    assert hop1["soma_position"].status == hop2["source_label"].status == "inferred"
    label_record = hop2["source_label"].to_observation(hop2.request)
    assert [(s.variable, dict(s.evidence)) for s in label_record.lineage] == [
        ("source_label", {"soma_position": "inferred"}), ("soma_position", {"regional_axon_length": "measured"})]
    assert dict(label_record.derived_from) == {"soma_position": "inferred"}       # no profile in its own step
    with pytest.raises(core.ForbiddenInputError, match="descends from a measured regional_axon_length"):
        c.adapter.check_inputs(c.group, {"source_label": label_record})
    with pytest.raises(core.ForbiddenInputError, match="descends from a measured regional_axon_length"):
        c.adapter.check_inputs(c.group, {"soma_position": c.m_soma, "source_label": label_record})
    for target in ("axon_arbor", "regional_axon_length"):                          # the adapter, and T1
        refused = c.held_out.reconstruct(c.ask({"soma_position": c.soma, "source_label": label_record}, [target]))[target]
        assert refused.status == "refused" and refused.refusal.code == "forbidden_input", target
    # the chain cannot be shortened: a record refuses a truncated lineage, and a lineage cut behind the
    # constructor is refused as well (the record's invariants are checked again on the object as it is; and it
    # is no longer the record that was issued)
    with pytest.raises(ValueError, match="incomplete"):
        dataclasses.replace(label_record, lineage=label_record.lineage[:1])
    cut = dataclasses.replace(label_record)
    object.__setattr__(cut, "lineage", label_record.lineage[:1])
    assert c.ledger.fact(cut) is None and c.ledger.fact(label_record) is not None
    with pytest.raises(core.ForbiddenInputError, match="lineage is incomplete|changed after it was issued"):
        c.adapter.check_inputs(c.group, {"source_label": cut})
    emptied = dataclasses.replace(label_record)
    object.__setattr__(emptied, "lineage", ())
    with pytest.raises(core.ForbiddenInputError):
        c.adapter.check_inputs(c.group, {"source_label": emptied})


def test_K1d_inferred_input_from_a_state_that_holds_the_query_brain_sample_is_refused(registry, chain):
    """[K1 d] A label inferred from the neuron's own measured soma is a permitted kind of input - unless the
    state that produced it was fitted with the neuron's brain sample. The training groups are those of the
    registered state as the layer recorded them, not what the record or a registration claims, and the sample
    is not hidden by its spelling."""
    c = chain
    hop = c.with_A.reconstruct(c.ask({"soma_position": c.soma}, ["source_label"]))             # by state A
    assert hop["source_label"].status == "inferred" and "group_in_training_state" in hop["source_label"].codes()
    label_A = hop["source_label"].to_observation(hop.request)
    assert label_A.uncertainty.state_hash == c.state_A.state_hash() and c.group in c.state_A.groups
    with pytest.raises(core.ForbiddenInputError, match="contains the query's brain sample"):
        c.adapter.check_inputs(c.group, {"source_label": label_A})
    with pytest.raises(core.ForbiddenInputError, match="contains the query's brain sample"):
        c.adapter.check_inputs(c.group, {"soma_position": c.m_soma, "source_label": label_A})
    for target in ("axon_arbor", "regional_axon_length"):                                      # the adapter, and T1
        refused = c.held_out.reconstruct(c.ask({"soma_position": c.soma, "source_label": label_A}, [target]))[target]
        assert refused.status == "refused" and refused.refusal.code == "forbidden_input", target
        assert "contains the query's brain sample" in refused.refusal.reason
    # the same inference by state B, fitted without the sample, is admissible: the state is the whole difference
    hop_B = c.lenient_B.reconstruct(c.ask({"soma_position": c.soma}, ["source_label"]))
    label_B = hop_B["source_label"].to_observation(hop_B.request)
    assert "group_in_training_state" not in hop_B["source_label"].codes()
    assert c.adapter.check_inputs(c.group, {"soma_position": c.m_soma, "source_label": label_B})
    # spelling the sample differently does not make it another sample
    for spelled in (c.group + " ", " " + c.group, c.group + "\n", "", c.group.upper()):
        with pytest.raises(core.ForbiddenInputError):
            c.adapter.check_inputs(spelled, {"source_label": label_A})
    # The state rule is not only the adapter's: it holds for every answer that is not one of the declared
    # profile-input tasks, also one that is not axon-derived (here: the label from an inferred profile).
    by_A = c.with_A.reconstruct(c.ask({"soma_position": c.soma}, ["regional_axon_length"]))
    by_B = c.lenient_B.reconstruct(c.ask({"soma_position": c.soma}, ["regional_axon_length"]))
    profile_A = by_A["regional_axon_length"].to_observation(by_A.request)
    profile_B = by_B["regional_axon_length"].to_observation(by_B.request)
    refused = c.lenient_B.reconstruct(c.ask({"regional_axon_length": profile_A}, ["source_label"]))["source_label"]
    assert refused.status == "refused" and refused.refusal.code == "forbidden_input"
    assert "contains the query's brain sample" in refused.refusal.reason
    admitted = c.lenient_B.reconstruct(c.ask({"regional_axon_length": profile_B}, ["source_label"]))["source_label"]
    assert admitted.status == "inferred" and "conditioned_on_inferred_input" in admitted.codes()
    # The rule has no exception for a request that has the shape of a DECLARED profile-input task: what is
    # supplied there as the 'profile' is an inferred value like any other, held to the same state rule. The
    # same chain by state B is answered, flagged, and without the coverage flag.
    from_label_A = c.with_A.reconstruct(c.ask({"source_label": c.label}, ["regional_axon_length"]))
    profile_by_A = from_label_A["regional_axon_length"].to_observation(from_label_A.request)
    for reason in ("is not registered", "contains the query's brain sample"):
        declared = c.chained.reconstruct(c.ask({"regional_axon_length": profile_by_A}, ["soma_position"]))["soma_position"]
        assert declared.status == "refused" and declared.refusal.code == "forbidden_input"
        assert declared.conditional == "T2b_soma_from_profile" and declared.estimate is None
        assert reason in declared.refusal.reason
        c.chained.module("axon_arbor_adapter").register_state(c.state_A)      # known to that layer from now on
    from_label_B = c.lenient_B.reconstruct(c.ask({"source_label": c.label}, ["regional_axon_length"]))
    flagged = c.chained.reconstruct(c.ask(
        {"regional_axon_length": from_label_B["regional_axon_length"].to_observation(from_label_B.request)},
        ["soma_position"]))["soma_position"]
    assert flagged.status == "inferred" and flagged.conditional == "T2b_soma_from_profile"
    assert "conditioned_on_inferred_input" in flagged.codes() and not flagged.coverage_validated
    # a registration that leaves the sample out is not believed; nor can a registration be silently replaced
    without = [g for g in c.state_A.groups if g != c.group]
    liar = mods.AxonArborAdapter(registry, ledger=c.ledger)
    liar.register_state(c.state_A.state_hash(), without, ("soma_position", "source_label"))
    with pytest.raises(core.ForbiddenInputError, match="not the ones recorded"):
        liar.check_inputs(c.group, {"source_label": label_A})
    with pytest.raises(core.RegistryError):
        c.adapter.register_state(c.state_A.state_hash(), without)
    with pytest.raises(core.RegistryError):
        c.adapter.register_state(c.state_A, without)
    with pytest.raises(core.RegistryError):
        mods.AxonArborAdapter(registry).register_state("some-state", [])


def test_K1e_missing_unregistered_or_unverifiable_provenance_and_unknown_group_are_refused(registry, chain):
    """[K1 e] An inferred input is admitted only when its provenance can be verified: it is exactly a value this
    layer issued, for this neuron, by a registered state, and the query's brain sample is known."""
    c = chain
    hop = c.lenient_B.reconstruct(c.ask({"soma_position": c.soma}, ["source_label"]))
    clean = hop["source_label"].to_observation(hop.request)
    assert c.adapter.check_inputs(c.group, {"source_label": clean})                 # admissible exactly as issued
    with pytest.raises(core.ForbiddenInputError, match="brain sample is unknown"):  # unknown query group
        c.adapter.check_inputs(None, {"source_label": clean})
    with pytest.raises(core.ForbiddenInputError, match="not registered"):           # unregistered producing state
        mods.AxonArborAdapter(registry, ledger=c.ledger).check_inputs(c.group, {"source_label": clean})
    with pytest.raises(core.ForbiddenInputError, match="no ledger"):                # nothing to verify against
        mods.AxonArborAdapter(registry).check_inputs(c.group, {"source_label": clean})
    elsewhere = core.Reconstructor(registry)                                        # a layer that never issued it
    stranger = elsewhere.register(mods.AxonArborAdapter(registry))
    stranger.register_state(c.state_B)
    with pytest.raises(core.ForbiddenInputError, match="not issued by this reconstruction layer"):
        stranger.check_inputs(c.group, {"source_label": clean})
    # never issued: built by hand with a lineage that looks clean and names the registered state B
    other_value = next(v for v in LABELS if v != clean.value)
    looks_clean = c.handmade("source_label", other_value, (("soma_position", "measured"),),
                             state_hash=c.state_B.state_hash(), groups_hash=c.state_B.groups_sha256)
    assert looks_clean.lineage == clean.lineage and looks_clean.derived_from == clean.derived_from
    with pytest.raises(core.ForbiddenInputError, match="not issued by this reconstruction layer"):
        c.adapter.check_inputs(c.group, {"source_label": looks_clean})
    for target in ("axon_arbor", "regional_axon_length"):                           # the adapter, and T1
        refused = c.held_out.reconstruct(c.ask({"soma_position": c.soma, "source_label": looks_clean}, [target]))[target]
        assert refused.status == "refused" and refused.refusal.code == "forbidden_input", target
    # ... and an answer that is not axon-derived (the label from a hand-made 'inferred' profile) is refused as well;
    # so is a request that has the shape of a declared profile-input task: no task takes an inferred input on its
    # own word
    made_up = c.handmade("regional_axon_length", tuple(c.profile[r] for r in REGIONS), (("soma_position", "measured"),),
                         coordinates=REGIONS, state_hash=c.state_B.state_hash(), groups_hash=c.state_B.groups_sha256)
    unverified = c.lenient_B.reconstruct(c.ask({"regional_axon_length": made_up}, ["source_label"]))["source_label"]
    assert unverified.status == "refused" and unverified.refusal.code == "forbidden_input"
    assert "not issued by this reconstruction layer" in unverified.refusal.reason
    from_label = c.handmade("regional_axon_length", tuple(c.profile[r] for r in REGIONS), (("source_label", "measured"),),
                            coordinates=REGIONS, state_hash=c.state_B.state_hash(), groups_hash=c.state_B.groups_sha256)
    declared = c.chained.reconstruct(c.ask({"regional_axon_length": from_label}, ["soma_position"]))["soma_position"]
    assert declared.status == "refused" and declared.refusal.code == "forbidden_input"
    assert declared.conditional == "T2b_soma_from_profile"
    assert "not issued by this reconstruction layer" in declared.refusal.reason
    # issued, then changed: another value, or a lineage rewritten behind the constructor
    with pytest.raises(core.ForbiddenInputError, match="changed after it was issued"):
        c.adapter.check_inputs(c.group, {"source_label": dataclasses.replace(clean, value=other_value)})
    rewritten = dataclasses.replace(clean)
    object.__setattr__(rewritten, "lineage", (dataclasses.replace(clean.lineage[0], state_hash="another-state"),))
    with pytest.raises(core.ForbiddenInputError):
        c.adapter.check_inputs(c.group, {"source_label": rewritten})
    # a value inferred for ANOTHER neuron cannot be re-issued for this one, and is refused if it is forced
    neighbour = c.lenient_B.reconstruct(c.ask({"soma_position": c.soma}, ["source_label"], unit=f"{c.group}_88888"))
    with pytest.raises(core.PairingError):
        neighbour["source_label"].to_observation(hop.request)
    theirs = neighbour["source_label"].to_observation(neighbour.request)
    with pytest.raises(core.ForbiddenInputError):
        c.adapter.check_inputs(c.group, {"source_label": theirs}, unit_id=c.unit)
    with pytest.raises(core.ForbiddenInputError, match="different units"):
        c.adapter.check_inputs(c.group, {"soma_position": c.m_soma, "source_label": theirs})
    # a record that says 'measured' cannot carry inferred provenance, and a name does not relabel a record
    with pytest.raises(ValueError):
        dataclasses.replace(c.m_soma, derived_from=(("regional_axon_length", "measured"),))
    smuggled = dataclasses.replace(c.m_soma)
    object.__setattr__(smuggled, "derived_from", (("regional_axon_length", "measured"),))
    with pytest.raises(core.ForbiddenInputError, match="says 'measured' but carries"):
        c.adapter.check_inputs(c.group, {"soma_position": smuggled})
    profile_record = c.ask({"regional_axon_length": c.profile}, []).observed[0]
    with pytest.raises(core.ForbiddenInputError, match="supplied under the name soma_position"):
        c.adapter.check_inputs(c.group, {"soma_position": profile_record})
    # the neuron's brain sample is the one its ID names: another one cannot be claimed
    claimed = c.with_A.reconstruct(c.ask({"soma_position": c.soma}, ["source_label"], group="another_sample"))
    assert claimed["source_label"].status == "refused" and claimed["source_label"].refusal.code == "invalid_request"
    with pytest.raises(core.ForbiddenInputError, match="by its ID"):
        c.adapter.check_inputs("another_sample", {"source_label": clean})


def test_K1f_a_legitimate_chain_passes_the_check_and_the_adapter_stays_blocked(registry, chain):
    """[K1 f] The check admits what the contract permits: the neuron's measured permitted inputs, and values
    inferred from them (in one step or in two, by the engine or by a comparator bank) by states that do not hold
    its brain sample. Passing the check changes nothing about the adapter's status: it stays blocked. The twin
    of the same chain that touches the measured profile is refused."""
    c = chain
    assert c.adapter.check_inputs(c.group, {"soma_position": c.m_soma, "source_label": c.m_label})
    hop1 = c.lenient_B.reconstruct(c.ask({"soma_position": c.soma}, ["source_label"]))
    label_record = hop1["source_label"].to_observation(hop1.request)               # label <- measured soma
    assert c.adapter.check_inputs(c.group, {"soma_position": c.m_soma, "source_label": label_record})
    hop2 = c.lenient_B.reconstruct(c.ask({"source_label": label_record}, ["regional_axon_length"]))
    assert hop2["regional_axon_length"].status == "inferred"                        # profile <- that inferred label
    profile_record = hop2["regional_axon_length"].to_observation(hop2.request)
    assert [s.variable for s in profile_record.lineage] == ["regional_axon_length", "source_label"]
    assert lineage_roots_of(profile_record) == {("soma_position", "measured")}
    assert c.adapter.check_inputs(c.group, {"soma_position": c.m_soma, "source_label": label_record,
                                            "regional_axon_length": profile_record})
    base = c.baseline_B.reconstruct(c.ask({"soma_position": c.soma, "source_label": c.label}, ["regional_axon_length"]))
    assert base["regional_axon_length"].baseline_level                              # a comparator bank's value
    bank_record = base["regional_axon_length"].to_observation(base.request)
    assert c.baseline_B.module("axon_arbor_adapter").check_inputs(c.group, {"regional_axon_length": bank_record,
                                                                           "soma_position": c.m_soma})
    # through the request API the inputs are accepted and the answer is withheld: the adapter is blocked
    for observed in ({"soma_position": c.soma, "source_label": c.label},
                     {"soma_position": c.soma, "source_label": label_record},
                     {"soma_position": c.soma, "source_label": label_record, "regional_axon_length": profile_record}):
        answer = c.lenient_B.reconstruct(c.ask(observed, ["axon_arbor"]))["axon_arbor"]
        assert answer.status == "refused" and answer.refusal.code == "not_validated" and answer.estimate is None
        assert answer.refusal.reason == "blocked: independent validation pending" == c.adapter.STATUS
    # the twin that touches the hidden axon: the label inferred from soma + MEASURED profile; same state, same
    # neuron, same names - refused as a forbidden input, not merely withheld
    t3 = c.held_out.reconstruct(c.ask({"soma_position": c.soma, "regional_axon_length": c.profile}, ["source_label"]))
    dirty = t3["source_label"].to_observation(t3.request)
    twin = c.lenient_B.reconstruct(c.ask({"soma_position": c.soma, "source_label": dirty}, ["axon_arbor"]))["axon_arbor"]
    assert twin.status == "refused" and twin.refusal.code == "forbidden_input"
    with pytest.raises(core.ForbiddenInputError):
        c.adapter.check_inputs(c.group, {"soma_position": c.m_soma, "source_label": dirty})


# ---------------------------------------------------------------------------------------------- K1, second round
#
# Three ways in which a forbidden input still got through the first fix (found by an independent attack on it).
# Each test below fails on a copy of the code in which only its repair is switched off
# (data/interim/synarch_recon/fix_evidence/round2/).


def test_K1g_a_request_shaped_like_a_declared_profile_task_gets_no_exception_for_inferred_inputs(registry, chain):
    """[K1 g] A request with the SHAPE of task T4 (regions outside the block from the block inside + soma + label)
    that holds no measured profile at all: its three inputs are INFERRED. It is a forecast of hidden axon
    regions. The first fix skipped the whole provenance rule for every conditional declared as taking a measured
    profile, without looking whether a measured profile was there: so this forecast was answered from inputs
    inferred by a state that holds the neuron's brain sample, and from hand-built records that nobody issued.
    Now the rule applies to every inferred input of every request; only the MEASURED profile block of a declared
    task is an input by declaration."""
    c = chain
    outside = REGIONS[3:]
    t4 = "T4_outside_block_from_inside_block_soma_label"

    def issued(layer, observed, target):
        hop = layer.reconstruct(c.ask(observed, [target]))
        name = target if isinstance(target, str) else target[0]
        assert hop[name].status == "inferred", hop[name]
        return hop[name].to_observation(hop.request)

    def three(layer):          # inside block <- measured soma; soma <- measured label; label <- measured soma
        return {"regional_axon_length": issued(layer, {"soma_position": c.soma}, ("regional_axon_length", INSIDE)),
                "soma_position": issued(layer, {"source_label": c.label}, "soma_position"),
                "source_label": issued(layer, {"soma_position": c.soma}, "source_label")}

    by_A, by_B = three(c.with_A), three(c.lenient_B)
    assert {r.status for r in by_A.values()} == {"inferred"}
    assert all(r.uncertainty.state_hash == c.state_A.state_hash() for r in by_A.values()) and c.group in c.state_A.groups
    # 1. the reported case: three inputs inferred by a state that holds the query's brain sample
    answer = c.held_out.reconstruct(c.ask(by_A, [("regional_axon_length", outside)]))["regional_axon_length"]
    assert answer.conditional == t4                               # the request did match the declared task's shape
    assert answer.status == "refused" and answer.refusal.code == "forbidden_input" and answer.estimate is None
    assert "contains the query's brain sample" in answer.refusal.reason
    # 2. three hand-built 'inferred' records that the layer never issued; the soma value was computed from the
    #    neuron's WHOLE measured profile (task T2), and the records say nothing of it
    t2 = c.held_out.reconstruct(c.ask({"regional_axon_length": c.profile, "source_label": c.label}, ["soma_position"]))
    state, groups = c.state_B.state_hash(), c.state_B.groups_sha256
    forged = {"soma_position": c.handmade("soma_position", tuple(t2["soma_position"].estimate.point),
                                          (("source_label", "measured"),), coordinates=core.SOMA_AXES,
                                          state_hash=state, groups_hash=groups),
              "source_label": c.handmade("source_label", "B", (("soma_position", "measured"),), state_hash=state,
                                         groups_hash=groups),
              "regional_axon_length": c.handmade("regional_axon_length", tuple(c.profile[r] for r in INSIDE),
                                                 (("soma_position", "measured"),), coordinates=INSIDE,
                                                 state_hash=state, groups_hash=groups)}
    answer = c.held_out.reconstruct(c.ask(forged, [("regional_axon_length", outside)]))["regional_axon_length"]
    assert answer.conditional == t4 and answer.status == "refused" and answer.refusal.code == "forbidden_input"
    assert "not issued by this reconstruction layer" in answer.refusal.reason
    # 3. the same holds when the profile block IS measured (the true task T4) and soma and label are inferred:
    #    a declared task takes its measured block by declaration, and nothing else on its word
    block = {r: c.profile[r] for r in INSIDE}
    for inputs, reason in ((by_A, "contains the query's brain sample"), (forged, "not issued by this reconstruction layer")):
        observed = {"regional_axon_length": block, "soma_position": inputs["soma_position"],
                    "source_label": inputs["source_label"]}
        answer = c.held_out.reconstruct(c.ask(observed, [("regional_axon_length", outside)]))["regional_axon_length"]
        assert answer.conditional == t4 and answer.status == "refused" and answer.refusal.code == "forbidden_input"
        assert reason in answer.refusal.reason
    # 4. what stays permitted. The true task with measured inputs: answered, coverage flag kept, no notice.
    measured_t4 = c.held_out.reconstruct(c.ask({"regional_axon_length": block, "soma_position": c.soma,
                                                "source_label": c.label},
                                               [("regional_axon_length", outside)]))["regional_axon_length"]
    assert measured_t4.status == "inferred" and measured_t4.conditional == t4
    assert measured_t4.coverage_validated and not measured_t4.codes()
    #    The same shapes with inputs inferred by state B, which never saw the brain sample: answered, each inferred
    #    input flagged, and no coverage claimed.
    for observed in (by_B, {"regional_axon_length": block, "soma_position": by_B["soma_position"],
                            "source_label": by_B["source_label"]}):
        legit = c.held_out.reconstruct(c.ask(observed, [("regional_axon_length", outside)]))["regional_axon_length"]
        assert legit.status == "inferred" and legit.conditional == t4 and not legit.coverage_validated
        assert legit.codes() == {"conditioned_on_inferred_input"}
        assert {n.variable for n in legit.notices} == {v for v, r in observed.items() if isinstance(r, core.Observation)}
    #    And an inferred value under the name of the profile is never taken for the measured block: the exception
    #    is for a measured record only.
    assert core.hidden_axon_problem(registry, "jiao", ("jiao", "neuron", c.unit), c.group,
                                    {"regional_axon_length": by_A["regional_axon_length"]}, ledger=c.ledger,
                                    resolve=c.held_out.resolve_state, axon_hidden=True,
                                    declared_measured=("regional_axon_length",))
    # 5. an issued value whose statement of scale is rewritten afterwards is no longer the record that was issued.
    #    In a request shaped like a declared task it used to be flagged only, and was then converted a second time.
    hop = c.lenient_B.reconstruct(c.ask({"source_label": c.label}, ["regional_axon_length"]))
    as_issued = hop["regional_axon_length"].to_observation(hop.request)
    spec = registry.variable("jiao", "regional_axon_length")
    relabelled = dataclasses.replace(as_issued, scale=core.scale_statement(spec, "identity"), transform="identity")
    taken = c.chained.reconstruct(c.ask({"regional_axon_length": as_issued}, ["soma_position"]))["soma_position"]
    assert taken.status == "inferred" and taken.conditional == "T2b_soma_from_profile"
    rewritten = c.chained.reconstruct(c.ask({"regional_axon_length": relabelled}, ["soma_position"]))["soma_position"]
    assert rewritten.status == "refused" and rewritten.refusal.code == "forbidden_input"
    assert "changed after it was issued" in rewritten.refusal.reason


def test_K1h_a_record_that_says_measured_but_names_a_producing_state_is_refused(registry, chain):
    """[K1 h] A record whose status says 'measured' while its uncertainty names the module, the state and the
    training groups that produced its value (its uncertainty KIND still says 'none stated'). The value here is
    the source label that was inferred WITH the neuron's measured profile (task T3). The first fix looked at the
    uncertainty's kind only, and re-validated the record but not its uncertainty: the adapter accepted the
    record, and T1 answered from it with the coverage flag and no notice."""
    c = chain
    t3 = c.held_out.reconstruct(c.ask({"soma_position": c.soma, "regional_axon_length": c.profile}, ["source_label"]))
    dirty = t3["source_label"].to_observation(t3.request)

    def relabelled(**fields):
        """A 'measured' record of the dirty value whose uncertainty carries `fields` behind the constructor."""
        record = c.ask({"source_label": dirty.value}, []).observed[0]
        uncertainty = core.Uncertainty("none_stated")
        for name, value in fields.items():
            object.__setattr__(uncertainty, name, value)
        object.__setattr__(record, "uncertainty", uncertainty)
        return record

    u = dirty.uncertainty
    forged = relabelled(producing_module=u.producing_module, state_hash=u.state_hash,
                        training_groups_hash=u.training_groups_hash, coverage="unvalidated")
    assert forged.status == "measured" and forged.value == dirty.value and forged.uncertainty.kind == "none_stated"
    assert forged.uncertainty.state_hash == c.state_B.state_hash() and not forged.lineage and not forged.derived_from
    with pytest.raises(core.ForbiddenInputError, match="says 'measured' but carries"):   # as issued it was refused too
        c.adapter.check_inputs(c.group, {"soma_position": c.m_soma, "source_label": forged})
    with pytest.raises(core.ForbiddenInputError, match="descends from a measured regional_axon_length"):
        c.adapter.check_inputs(c.group, {"soma_position": c.m_soma, "source_label": dirty})
    for target in ("axon_arbor", "regional_axon_length"):                               # the adapter, and T1
        answer = c.held_out.reconstruct(c.ask({"soma_position": c.soma, "source_label": forged}, [target]))[target]
        assert answer.status == "refused" and answer.refusal.code == "invalid_request" and answer.estimate is None
        assert "only inferred values carry" in answer.refusal.reason
    with pytest.raises(ValueError, match="only inferred values carry"):
        core.ObservationStore(registry).add(forged)
    assert "only inferred values carry" in core.invariant_problem(forged)
    # each of the four provenance fields alone is enough
    for name, value in (("producing_module", "joint_mixture"), ("state_hash", u.state_hash),
                        ("training_groups_hash", u.training_groups_hash), ("coverage", "validated")):
        one = relabelled(**{name: value})
        with pytest.raises(core.ForbiddenInputError, match="says 'measured' but carries"):
            c.adapter.check_inputs(c.group, {"source_label": one})
        refused = c.held_out.reconstruct(c.ask({"soma_position": c.soma, "source_label": one},
                                               ["regional_axon_length"]))["regional_axon_length"]
        assert refused.status == "refused" and refused.refusal.code == "invalid_request", name
    # the same family: a source that is a reconstruction, however the prefix is spelled; a source whose hash is a
    # fitted state of this layer; a part of the record that is not what it should be
    for path in ("Reconstruct:joint_mixture", " RECONSTRUCT:joint_mixture"):
        with pytest.raises(ValueError, match="'inferred', not 'measured'"):
            measured("source_label", dirty.value, unit=c.unit, group=c.group, source=core.SourceRef(path, None, c.unit))
    from_state = measured("source_label", dirty.value, unit=c.unit, group=c.group,
                          source=core.SourceRef("tables/labels.csv", c.state_B.state_hash(), c.unit))
    with pytest.raises(core.ForbiddenInputError, match="says 'measured' but carries"):
        c.adapter.check_inputs(c.group, {"source_label": from_state})
    refused = c.held_out.reconstruct(c.ask({"soma_position": c.soma, "source_label": from_state},
                                           ["regional_axon_length"]))["regional_axon_length"]
    assert refused.status == "refused" and refused.refusal.code == "forbidden_input"
    hijacked = c.ask({"source_label": c.label}, []).observed[0]
    object.__setattr__(hijacked, "key", lambda: ("jiao", "neuron", c.unit))        # a method replaced on the object
    assert "not record fields" in core.invariant_problem(hijacked)
    with pytest.raises(core.ForbiddenInputError):
        c.adapter.check_inputs(c.group, {"source_label": hijacked})
    not_a_group = c.ask({"source_label": c.label}, []).observed[0]
    object.__setattr__(not_a_group, "group", (c.group, "brain_sample"))
    assert core.invariant_problem(not_a_group)
    assert c.held_out.reconstruct(core.ReconstructionRequest(
        "jiao", "neuron", c.unit, c.m_soma.group, (c.m_soma, not_a_group),
        (core.Target("regional_axon_length"),)))["regional_axon_length"].refusal.code == "invalid_request"
    # the plain measured label, from a table, is what it says it is
    plain = measured("source_label", c.label, unit=c.unit, group=c.group,
                     source=core.SourceRef("tables/labels.csv", "ab" * 32, c.unit))
    assert core.invariant_problem(plain) is None
    assert c.adapter.check_inputs(c.group, {"soma_position": c.m_soma, "source_label": plain})
    answered = c.held_out.reconstruct(c.ask({"soma_position": c.soma, "source_label": plain},
                                            ["regional_axon_length"]))["regional_axon_length"]
    assert answered.status == "inferred" and answered.coverage_validated


def test_K1i_a_state_whose_training_units_contradict_its_group_list_is_refused(registry, layout, tasks, planted):
    """[K1 i] Whether a state holds the query's brain sample is read from the state's group list. A Jiao unit ID
    names its brain sample (<brain sample>_<n>); the first fix enforced that for the query and never for a
    state's own units. A state fitted on rows whose units are '<sample>_<n>' but whose group column says 'zz'
    lists no such sample: a label it inferred for a new neuron of that sample passed as coming from a state that
    'excludes' the sample, and its direct answer carried no seen-group notice. Such rows and such a state are
    now refused where a state enters the layer: a module, and the adapter's register."""
    rows = planted["rows"]
    sample = str(rows.group[0])
    relabelled = mods.Rows(rows.unit.copy(), np.where(rows.group == sample, "zz", rows.group), dict(rows.cont),
                           dict(rows.cat))
    hidden = sum(str(u).startswith(sample + "_") for u in relabelled.unit)
    state_M, _ = mods.fit_joint(relabelled, layout, components=2, seed=1, budget={"n_init": 1, "max_iter": 30, "tol": 1e-5})
    assert hidden > 0 and sample not in state_M.groups and "zz" in state_M.groups       # what the list hides
    assert sum(u.startswith(sample + "_") for u in state_M.units) == hidden
    with pytest.raises(core.ReconError, match="by its ID it belongs to"):                # as a module's state
        mods.JointMixtureModule(registry, layout, states={n: state_M for n in tasks}, train=relabelled, tasks=tasks)
    with pytest.raises(core.RegistryError, match="cannot be shown to exclude"):          # in the adapter's register
        mods.AxonArborAdapter(registry).register_state(state_M)
    assert "not in both" in mods.group_naming_problem(registry, "jiao", sorted(state_M.units), state_M.groups)
    # two neurons of two brain samples exchange their group: the SET of groups is unchanged, the rows are wrong,
    # and a comparator that excludes 'the query's own brain sample' would copy from it
    i, j = 0, int(np.flatnonzero(rows.group != rows.group[0])[0])
    exchanged_groups = rows.group.copy()
    exchanged_groups[[i, j]] = rows.group[[j, i]]
    exchanged = mods.Rows(rows.unit.copy(), exchanged_groups, dict(rows.cont), dict(rows.cat))
    assert set(exchanged.group) == set(rows.group)
    assert mods.group_naming_problem(registry, "jiao", exchanged.unit, exchanged.group, aligned=True)
    with pytest.raises(core.ReconError, match="by its ID it belongs to"):
        mods.JointMixtureModule(registry, layout, states={n: planted["state"] for n in tasks}, train=exchanged,
                                tasks=tasks)
    # consistent rows and their state are accepted as before; a namespace whose IDs name no group has no rule
    assert mods.group_naming_problem(registry, "jiao", rows.unit, rows.group, aligned=True) is None
    assert mods.group_naming_problem(registry, "jiao", sorted(planted["state"].units), planted["state"].groups) is None
    mods.JointMixtureModule(registry, layout, states={n: planted["state"] for n in tasks}, train=rows, tasks=tasks)
    assert mods.AxonArborAdapter(registry).register_state(planted["state"]) == planted["state"].state_hash()
    assert mods.group_naming_problem(registry, "epiretro", ["n1", "n2"], ["experiment_7"]) is None


def test_identifiers_with_invisible_characters_and_malformed_adapter_records_are_refused(registry, chain):
    """Hardening next to K1. An identifier with a zero-width space, a byte-order mark, a soft hyphen, a joiner
    or a C1 control character anywhere in it is refused: with one of them a training sample looked like a new
    one. And the adapter refuses a record whose coordinates are not coordinates of its variable."""
    c = chain
    for invisible in ("​", "﻿", "­", "‍", "\u0080", "\x00", "\x7f"):
        spelled = c.group[:2] + invisible + c.group[2:]
        assert core.identifier_problem(spelled) and core.identifier_problem(spelled + "_99999", "the unit ID")
        answer = c.with_A.reconstruct(c.ask({"soma_position": c.soma, "source_label": c.label},
                                            ["regional_axon_length"], unit=f"{spelled}_99999",
                                            group=spelled))["regional_axon_length"]
        assert answer.status == "refused" and answer.refusal.code == "invalid_request", repr(invisible)
        with pytest.raises(core.ForbiddenInputError):
            c.adapter.check_inputs(spelled, {"soma_position": c.m_soma, "source_label": c.m_label})
    assert core.identifier_problem(c.group) is None and core.identifier_problem("001 CLA-EPd-CTX Car3 Glut") is None
    six = measured("soma_position", tuple(c.profile[r] for r in REGIONS), REGIONS, unit=c.unit, group=c.group)
    with pytest.raises(core.ForbiddenInputError, match="not coordinates of the registered variable"):
        c.adapter.check_inputs(c.group, {"soma_position": six, "source_label": c.m_label})
    lettered = measured("hemisphere", ("L",) * 6, REGIONS, unit=c.unit, group=c.group)
    with pytest.raises(core.ForbiddenInputError, match="only a vector has coordinates"):
        c.adapter.check_inputs(c.group, {"hemisphere": lettered})
    bare = measured("soma_position", c.soma, (), unit=c.unit, group=c.group)
    with pytest.raises(core.ForbiddenInputError, match="not coordinates of the registered variable"):
        c.adapter.check_inputs(c.group, {"soma_position": bare})

    class Shouting(str):                    # a label whose str() is not its characters
        def __str__(self):
            return "B"

    assert core.plain_label(Shouting("A")) == "A" and core.plain_label(7) == "7"


def lineage_roots_of(record):
    return core.lineage_roots(record)


# ---------------------------------------------------------------------------------------------- K2


class FixedProfileModule(core.Module):
    """A second registered module: 'regional profile from the source label', always the same value on the
    modelling scale log10(1 + length). `stated=False` makes it state its scale by text only, without the value
    tag and the frame, as a careless module would."""
    namespace, grain, kind = "jiao", "neuron", "estimator"
    STATE, GROUPS = "f" * 64, ("some_other_sample",)

    def __init__(self, registry, coordinates, value, *, stated=True, name="fixed_profile"):
        self.registry, self.coordinates, self.value, self.stated, self.name = registry, tuple(coordinates), value, stated, name
        self.conditional = core.Conditional("fixed_profile_from_label", "regional_axon_length", None,
                                            (("source_label", None),))

    def produces(self):
        return ("regional_axon_length",)

    def conditionals(self):
        return (self.conditional,)

    def evidence_for(self, conditional_name):
        return core.ValidationEvidence("exploratory_supported", benchmark=POINTER)

    def state_groups(self, state_hash):
        return frozenset(self.GROUPS) if state_hash == self.STATE else None

    def plan(self, request, target, evidence, policy):
        if set(evidence) != {"source_label"}:
            return core.Refusal("invalid_request", "this module answers from the source label only")
        return core.Plan(self, target, self.coordinates, self.conditional, (("source_label", "measured"),))

    def provenance(self, plan):
        return {"state_hash": self.STATE, "training_groups_hash": core.ids_sha(self.GROUPS),
                "training_groups": self.GROUPS}

    def estimate(self, plan, level):
        spec = self.registry.variable("jiao", "regional_axon_length")
        point = (self.value,) * len(self.coordinates)
        if not self.stated:
            return core.Estimate("regional_axon_length", "vector", self.coordinates, "log10(1 + value)", point)
        return core.Estimate("regional_axon_length", "vector", self.coordinates,
                             core.scale_statement(spec, "log10_1p"), point, transform="log10_1p", frame=spec.frame)


def _conditioning(module, layout, ask, variable, coordinates=None, policy=None):
    """The values the joint module would condition on for `variable`, on the modelling scale."""
    plan = module.plan(ask, ask.requested[0], {o.variable: o for o in ask.observed}, policy or core.Policy())
    assert isinstance(plan, core.Plan), plan
    return [plan.payload["cont"][f] for f in layout.fields_of(variable, coordinates)]


def _model_means(state, observations, fields, clip=False):
    """The UNCHANGED JointModel.infer() at model-space observations: conditional means of `fields`."""
    predictions = state.model.infer(observations, interval_level=0.9)["predictions"]
    return [max(predictions[f]["mean"], 0.0) if clip else predictions[f]["mean"] for f in fields]


@pytest.fixture(scope="module")
def scale_api(registry, layout, tasks, planted):
    """The joint module behind a second module that returns a fixed profile, undeclared conditionals answered."""
    policy = core.Policy(on_not_validated="engine_flagged")
    module = mods.JointMixtureModule(registry, layout, states={n: planted["state"] for n in tasks},
                                     train=planted["rows"], tasks=tasks,
                                     evidence=evidence_for(tasks, "exploratory_supported"))
    api = core.Reconstructor(registry, policy)
    api.register(FixedProfileModule(registry, REGIONS, 2.0))
    api.register(module)
    return api, module, policy


def test_K2a_round_trip_estimate_to_record_to_evidence_keeps_the_model_space_value(registry, layout, tasks, planted,
                                                                                  scale_api):
    """[K2 a] A measured length of 99 um is 2.0 on the modelling scale log10(1 + length). An ESTIMATE of 2.0, made
    into a record and used as evidence, must condition the next answer at 2.0: the same number on the same
    scale. The reviewed code transformed the record again and conditioned at log10(1 + 2.0) = 0.477121."""
    api, module, policy = scale_api
    spec = registry.variable("jiao", "regional_axon_length")
    fields = layout.fields_of("regional_axon_length")
    first = api.reconstruct(request(registry, {"source_label": "A"}, ["regional_axon_length"]))
    result = first["regional_axon_length"]
    assert result.status == "inferred" and result.module == "fixed_profile" and result.estimate.point == (2.0,) * 6
    record = result.to_observation(first.request)
    assert record.value == (2.0,) * 6 and record.status == "inferred"
    assert (record.scale, record.transform, record.frame) == ("log10(1 + value)", "log10_1p", spec.frame)
    chained = request(registry, {"regional_axon_length": record}, ["soma_position"])
    assert _conditioning(module, layout, chained, "regional_axon_length", policy=policy) == [2.0] * 6    # exactly
    as_measured = request(registry, {"regional_axon_length": {r: 99.0 for r in REGIONS}}, ["soma_position"])
    np.testing.assert_allclose(_conditioning(module, layout, as_measured, "regional_axon_length", policy=policy),
                               2.0, rtol=0, atol=1e-12)
    # the chained answer is the unchanged model at 2.0, and visibly not the model at 0.477121
    answer = api.reconstruct(chained)["soma_position"]
    soma_fields = layout.fields_of("soma_position")
    at_two = _model_means(planted["state"], {f: 2.0 for f in fields}, soma_fields)
    at_twice = _model_means(planted["state"], {f: float(np.log10(3.0)) for f in fields}, soma_fields)
    assert answer.status == "inferred"
    np.testing.assert_allclose(answer.estimate.point, at_two, rtol=1e-9)
    assert np.abs(np.array(at_two) - np.array(at_twice)).max() > 1.0               # micrometres apart
    # the same round trip for the joint module's own estimate: profile from the label, then soma from that profile
    own = api.module("joint_mixture")
    alone = core.Reconstructor(registry, policy)
    alone.register(mods.JointMixtureModule(registry, layout, states=own.states, train=own.train, tasks=tasks))
    made = alone.reconstruct(request(registry, {"source_label": "B"}, ["regional_axon_length"]))
    estimate = made["regional_axon_length"].estimate
    back = made["regional_axon_length"].to_observation(made.request)
    again = request(registry, {"regional_axon_length": back}, ["soma_position"])
    assert _conditioning(module, layout, again, "regional_axon_length", policy=policy) == list(estimate.point)
    np.testing.assert_allclose(alone.reconstruct(again)["soma_position"].estimate.point,
                               _model_means(planted["state"], dict(zip(fields, estimate.point)), soma_fields), rtol=1e-9)
    # a module that states its scale by text only (no value tag, no frame) is stopped at the source
    careless = core.Reconstructor(registry, policy)
    careless.register(FixedProfileModule(registry, REGIONS, 2.0, stated=False))
    withheld = careless.reconstruct(request(registry, {"source_label": "A"}, ["regional_axon_length"]))["regional_axon_length"]
    assert withheld.status == "refused" and withheld.refusal.code == "invalid_estimate" and withheld.estimate is None


def test_K2b_a_measurement_scale_record_is_transformed_once_and_a_model_scale_record_not_at_all(registry, layout, tasks,
                                                                                               planted, scale_api):
    """[K2 b] One conversion, decided by what the record states: a value in length units is transformed once;
    a value that states the modelling scale is used as it is; a variable modelled on its own scale is never
    transformed."""
    api, module, policy = scale_api
    spec, soma_spec = registry.variable("jiao", "regional_axon_length"), registry.variable("jiao", "soma_position")
    lengths_um, on_model_scale = (99.0, 999.0, 99.0, 999.0, 99.0, 999.0), (2.0, 3.0, 2.0, 3.0, 2.0, 3.0)
    in_length_units = measured("regional_axon_length", lengths_um, REGIONS)
    already_there = measured("regional_axon_length", on_model_scale, REGIONS,
                             scale=core.scale_statement(spec, "log10_1p"), transform="log10_1p")
    np.testing.assert_allclose(core.model_values(spec, in_length_units), on_model_scale, rtol=0, atol=1e-12)
    assert core.model_values(spec, already_there).tolist() == list(on_model_scale)             # untouched
    answers = []
    for record in (in_length_units, already_there):
        ask = request(registry, {"regional_axon_length": record, "source_label": "A"}, ["soma_position"])
        np.testing.assert_allclose(_conditioning(module, layout, ask, "regional_axon_length", policy=policy),
                                   on_model_scale, rtol=0, atol=1e-12)
        answer = api.reconstruct(ask)["soma_position"]
        assert answer.status == "inferred" and answer.conditional == "T2_soma_from_profile_label"
        answers.append(answer.estimate.point)
    np.testing.assert_allclose(answers[0], answers[1], rtol=1e-9)
    np.testing.assert_allclose(answers[0], _model_means(planted["state"], dict(
        zip(layout.fields_of("regional_axon_length"), on_model_scale), source_label="A"),
        layout.fields_of("soma_position")), rtol=1e-9)
    somewhere = measured("soma_position", (6000.0, 6000.0, 5000.0), core.SOMA_AXES)
    assert core.model_values(soma_spec, somewhere).tolist() == [6000.0, 6000.0, 5000.0]
    ask = request(registry, {"soma_position": somewhere, "source_label": "A"}, ["regional_axon_length"])
    assert _conditioning(module, layout, ask, "soma_position", policy=policy) == [6000.0, 6000.0, 5000.0]


def test_K2c_a_record_with_a_missing_or_mismatched_scale_or_frame_is_refused(registry, layout, tasks, planted, scale_api):
    """[K2 c] Every record states the registered scale and frame of its variable. A statement that is missing,
    unknown, in other units, in another frame, or that contradicts the value tag is refused: by the request API,
    by the planner, by the store and by the axon adapter. Bad numbers are refusals too, not exceptions."""
    api, module, policy = scale_api
    soma_value = (6000.0, 6000.0, 5000.0)
    good_profile = measured("regional_axon_length", (99.0,) * 6, REGIONS)
    good_soma = measured("soma_position", soma_value, core.SOMA_AXES)
    label_record = measured("source_label", "A")
    bad_profiles = {
        "other units (0.099 mm is 99 um)": dict(value=(0.099,) * 6, scale="length >= 0 per region; mm"),
        "scale not stated": dict(scale=""),
        "unknown scale": dict(scale="z-score of the cohort"),
        "another frame": dict(frame="true hemispheres, 12 regions"),
        "frame not stated": dict(frame=None),
        "negative length": dict(value=(-5.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
        "non-finite length": dict(value=(float("nan"), 1.0, 1.0, 1.0, 1.0, 1.0)),
        "below zero on the modelling scale": dict(value=(-0.5,) * 6, scale="log10(1 + value)", transform="log10_1p"),
    }
    bad_somas = {
        "frame not stated": dict(frame=None),
        "another frame": dict(frame="true hemisphere, not mirrored"),
        "other units": dict(scale="linear position; mm"),
        "scale not stated": dict(scale=""),
        "a tag the variable does not have": dict(scale="log10(1 + value)", transform="log10_1p"),
    }
    plain_statement = ("scale not stated", "unknown scale", "other units", "other units (0.099 mm is 99 um)",
                       "another frame", "frame not stated", "a tag the variable does not have")
    cases = [("regional_axon_length", why, measured("regional_axon_length", **dict(dict(value=(99.0,) * 6,
                                                                                    coordinates=REGIONS), **change)),
              {"source_label": "A"}, "soma_position") for why, change in bad_profiles.items()]
    cases += [("soma_position", why, measured("soma_position", **dict(dict(value=soma_value, coordinates=core.SOMA_AXES),
                                                                     **change)),
               {"source_label": "A"}, "regional_axon_length") for why, change in bad_somas.items()]
    cases.append(("source_label", "a tag the variable does not have",
                  measured("source_label", "A", scale="log10(1 + value)", transform="log10_1p"),
                  {"soma_position": soma_value}, "regional_axon_length"))
    store = core.ObservationStore(registry)
    for variable, why, record, rest, target in cases:
        ask = request(registry, dict(rest, **{variable: record}), [target])
        refused = api.reconstruct(ask)[target]
        assert refused.status == "refused" and refused.refusal.code == "invalid_request", (variable, why, refused)
        plan = module.plan(ask, ask.requested[0], {o.variable: o for o in ask.observed}, policy)
        assert isinstance(plan, core.Refusal) and plan.code == "invalid_request", (variable, why)
        if why in plain_statement:                               # a statement problem: the store refuses the record too
            with pytest.raises(core.ScaleError):
                store.add(record)
            assert core.statement_problem(registry.variable("jiao", variable), record.scale, record.transform,
                                          record.frame)
    # the good twins of the same records are accepted everywhere
    assert api.reconstruct(request(registry, {"regional_axon_length": good_profile, "source_label": "A"},
                                   ["soma_position"]))["soma_position"].status == "inferred"
    assert api.reconstruct(request(registry, {"soma_position": good_soma, "source_label": label_record},
                                   ["regional_axon_length"]))["regional_axon_length"].status == "inferred"
    for record in (good_profile, good_soma, label_record):
        store.add(record)
    # a record cannot be built with a tag that contradicts its text, or without a tag; forced, it is refused
    for contradiction in (dict(transform="identity", scale="log10(1 + value)"),
                          dict(transform="log10_1p"), dict(transform=None), dict(transform="standardised")):
        with pytest.raises(ValueError):
            measured("regional_axon_length", (2.0,) * 6, REGIONS, **contradiction)
    forced = dataclasses.replace(good_profile)
    object.__setattr__(forced, "scale", "log10(1 + value)")       # says modelling scale, tagged measurement units
    refused = api.reconstruct(request(registry, {"regional_axon_length": forced, "source_label": "A"},
                                      ["soma_position"]))["soma_position"]
    assert refused.status == "refused" and refused.refusal.code == "invalid_request"
    # the axon adapter applies the same statement check to its inputs
    adapter = core.Reconstructor(registry).register(mods.AxonArborAdapter(registry))
    assert adapter.check_inputs("fresh", {"soma_position": good_soma, "source_label": label_record})
    for why in ("frame not stated", "another frame", "other units"):
        bad = measured("soma_position", **dict(dict(value=soma_value, coordinates=core.SOMA_AXES), **bad_somas[why]))
        with pytest.raises(core.ForbiddenInputError):
            adapter.check_inputs("fresh", {"soma_position": bad, "source_label": label_record})


def test_K2d_a_two_hop_chain_keeps_the_value_at_every_hop(registry, layout, tasks, planted):
    """[K2 d] label -> regions inside the hypothalamic block -> regions outside it -> soma. At every hop the
    planner conditions on exactly the numbers of the previous estimate, and every answer is the unchanged model
    at those numbers."""
    policy = core.Policy(on_not_validated="engine_flagged")
    r, module = reconstructor(registry, layout, tasks, planted, policy=policy)
    outside = REGIONS[3:]
    state = planted["state"]
    inside_fields, outside_fields = (layout.fields_of("regional_axon_length", c) for c in (INSIDE, outside))
    hop1 = r.reconstruct(request(registry, {"source_label": "B"}, [("regional_axon_length", INSIDE)]))
    a = hop1["regional_axon_length"]
    record_a = a.to_observation(hop1.request)
    assert a.status == "inferred" and record_a.coordinates == INSIDE and record_a.transform == "log10_1p"
    ask2 = request(registry, {"regional_axon_length": record_a}, [("regional_axon_length", outside)])
    assert _conditioning(module, layout, ask2, "regional_axon_length", INSIDE, policy) == list(a.estimate.point)
    hop2 = r.reconstruct(ask2)
    b = hop2["regional_axon_length"]
    assert b.status == "inferred" and b.estimate.coordinates == outside
    np.testing.assert_allclose(b.estimate.point, _model_means(state, dict(zip(inside_fields, a.estimate.point)),
                                                              outside_fields, clip=True), rtol=1e-9)
    record_b = b.to_observation(hop2.request)
    assert record_b.value == b.estimate.point and record_b.transform == "log10_1p"
    ask3 = request(registry, {"regional_axon_length": record_b}, ["soma_position"])
    assert _conditioning(module, layout, ask3, "regional_axon_length", outside, policy) == list(b.estimate.point)
    hop3 = r.reconstruct(ask3)["soma_position"]
    np.testing.assert_allclose(hop3.estimate.point, _model_means(state, dict(zip(outside_fields, b.estimate.point)),
                                                                 layout.fields_of("soma_position")), rtol=1e-9)
    assert hop3.estimate.transform == "identity" and hop3.estimate.scale == "linear position; um"
    record_c = hop3.to_observation(request(registry, {}, []))
    assert [s.variable for s in record_c.lineage] == ["soma_position", "regional_axon_length", "regional_axon_length"]
    # had a value been transformed again on the way, the last answer would be somewhere else
    twice = _model_means(state, dict(zip(outside_fields, np.log10(1.0 + np.array(b.estimate.point)))),
                         layout.fields_of("soma_position"))
    assert np.abs(np.array(twice) - np.array(hop3.estimate.point)).max() > 1.0


def test_K2e_the_inferred_input_warning_is_kept_and_nothing_is_called_validated(registry, layout, tasks, planted):
    """[K2 e] Converting a chained estimate correctly does not make it a measurement. The answer still says that
    it was conditioned on an inferred input whose uncertainty is not propagated, and it does not claim validated
    coverage, even for a conditional whose benchmark did validate it for measured inputs."""
    declared = dict(tasks)
    declared["T2b"] = mods.make_task(layout, registry, task="T2b", title="soma_from_profile", target="soma_position",
                                     evidence=(("regional_axon_length", None),), declared_measured_axon_input=True)
    module = mods.JointMixtureModule(registry, layout, states={n: planted["state"] for n in declared},
                                     train=planted["rows"], tasks=declared,
                                     evidence=evidence_for(declared, "exploratory_supported", coverage=True))
    r = core.Reconstructor(registry, core.Policy(on_not_validated="engine_flagged"))
    r.register(module)
    _, soma, y, label = planted_rows(1, 81, prefix="f")
    from_measured = r.reconstruct(request(registry, {"regional_axon_length": lengths(y[0])}, ["soma_position"]))["soma_position"]
    assert from_measured.conditional == "T2b_soma_from_profile" and from_measured.coverage_validated
    assert not from_measured.codes()
    made = r.reconstruct(request(registry, {"source_label": label[0]}, ["regional_axon_length"]))
    record = made["regional_axon_length"].to_observation(made.request)
    assert record.uncertainty.coverage == "unvalidated"
    ask = request(registry, {"regional_axon_length": record}, ["soma_position"])
    chained = r.reconstruct(ask)["soma_position"]
    assert chained.status == "inferred" and chained.conditional == "T2b_soma_from_profile"
    assert chained.validation == "exploratory_supported" and not chained.coverage_validated
    notice = next(n for n in chained.notices if n.code == "conditioned_on_inferred_input")
    assert notice.variable == "regional_axon_length" and "uncertainty is not propagated" in notice.detail
    assert dict(chained.evidence_used) == {"regional_axon_length": "inferred"}
    assert chained.to_observation(ask).uncertainty.coverage == "unvalidated"
    # and the value it was conditioned on is the estimate's own: the answer is the unchanged model there
    fields = layout.fields_of("regional_axon_length")
    np.testing.assert_allclose(chained.estimate.point, _model_means(planted["state"], dict(zip(fields, record.value)),
                                                                    layout.fields_of("soma_position")), rtol=1e-9)
    assert _conditioning(module, layout, ask, "regional_axon_length") == list(made["regional_axon_length"].estimate.point)


# ---------------------------------------------------------------------------------------------- K2, second round
#
# The conversion itself held under attack. What did not: the FRAME of a record was validated as text only.


class FixedSomaModule(core.Module):
    """A registered module that answers 'soma position from the source label' with one fixed point."""
    namespace, grain, kind, name = "jiao", "neuron", "estimator", "fixed_soma"
    STATE, GROUPS = "e" * 64, ("some_other_sample",)

    def __init__(self, registry, point):
        self.registry, self.point = registry, tuple(point)
        self.conditional = core.Conditional("fixed_soma_from_label", "soma_position", None, (("source_label", None),))

    def produces(self):
        return ("soma_position",)

    def conditionals(self):
        return (self.conditional,)

    def evidence_for(self, conditional_name):
        return core.ValidationEvidence("exploratory_supported", benchmark=POINTER)

    def state_groups(self, state_hash):
        return frozenset(self.GROUPS) if state_hash == self.STATE else None

    def plan(self, request, target, evidence, policy):
        return core.Plan(self, target, core.SOMA_AXES, self.conditional, (("source_label", "measured"),))

    def provenance(self, plan):
        return {"state_hash": self.STATE, "training_groups_hash": core.ids_sha(self.GROUPS), "training_groups": self.GROUPS}

    def estimate(self, plan, level):
        spec = self.registry.variable("jiao", "soma_position")
        return core.Estimate("soma_position", "vector", core.SOMA_AXES, core.scale_statement(spec), self.point,
                             transform="identity", frame=spec.frame)


def mirrored_cohort(registry):
    """A small MIRRORED cohort (every ML <= 5,700) with a class of neurons at the midline, and its state."""
    regions = ("PVH", "DMH", "ZI")
    layout = mods.jiao_layout(registry, regions)
    rng = np.random.default_rng(7)
    n = 900
    cls = rng.integers(0, 3, n)                              # classes at the midline, 140 um and 1,500 um lateral
    raw = np.array([5690.0, 5560.0, 4200.0])[cls] + rng.normal(0, 1.0, n) * np.array([40.0, 40.0, 160.0])[cls]
    cont = {"soma_position:AP": rng.normal(6400, 120, n), "soma_position:DV": rng.normal(6100, 120, n),
            "soma_position:ML": core.MIDLINE_ML - np.abs(core.MIDLINE_ML - raw)}
    logs = np.array([[3.6, 1.0, 2.0], [1.0, 3.6, 2.0], [2.0, 2.0, 3.6]])[cls] + rng.normal(0, 0.2, (n, 3))
    cont.update({f"regional_axon_length:{r}": np.maximum(logs[:, j], 0.02) for j, r in enumerate(regions)})
    group = np.array([f"s{i % 30:02d}" for i in range(n)])
    rows = mods.Rows(np.array([f"{g}_{i:04d}" for i, g in enumerate(group)]), group, cont,
                     {"source_label": np.array(["Cre"] * n, dtype=object)})
    state, _ = mods.fit_joint(rows, layout, components=3, seed=1, budget={"n_init": 3, "max_iter": 200, "tol": 1e-7})
    return regions, layout, rows, state


def test_K2f_a_value_outside_the_frame_its_record_states_is_refused(registry, bare_world):
    """[K2 f] A record states its frame in words, and the words were all that was checked. A soma in TRUE
    right-hemisphere coordinates (ML > 5,700) carried the registered statement 'mirrored to ML <= 5,700' and
    was answered at the unmirrored number: with coverage claimed and no notice, as long as it lay inside the
    state's support margin. The value is now held to the frame it states, by the one rule the loader applies to
    a training row: by the request API, the planner, the record store and the axon adapter; and nothing the
    layer issues under a frame's statement lies outside that frame."""
    regions, lay, rows, state = mirrored_cohort(registry)
    tasks = mods.build_tasks(lay, registry, regions[:1])
    evidence = {t.conditional().name: core.ValidationEvidence("exploratory_supported", task=k, benchmark=POINTER,
                                                              coverage_validated=True) for k, t in tasks.items()}
    api = core.Reconstructor(registry)
    module = api.register(mods.JointMixtureModule(registry, lay, states={k: state for k in tasks}, train=rows,
                                                  tasks=tasks, evidence=evidence))
    adapter = api.register(mods.AxonArborAdapter(registry))
    spec = registry.variable("jiao", "soma_position")

    def t1(ml):
        ask = core.ReconstructionRequest.from_values(
            registry, namespace="jiao", unit_id="q1_0001", group_id="q1", requested=["regional_axon_length"],
            observed={"soma_position": (6400.0, 6100.0, ml), "source_label": "Cre"})
        return ask, api.reconstruct(ask)["regional_axon_length"]

    low, high = state.ranges["soma_position:ML"]
    beyond = core.MIDLINE_ML + 75.0                         # 75 um into the right hemisphere, NOT mirrored
    assert high <= core.MIDLINE_ML and state.support_margin * (high - low) > 75.0
    assert not state.outside_range("soma_position:ML", beyond)          # the support check alone lets it through
    ask, refused = t1(beyond)
    soma, label = ask.observed
    assert "mirrored to ML <= 5,700" in soma.frame and soma.value[2] == beyond
    assert core.statement_problem(spec, soma.scale, soma.transform, soma.frame) is None      # its words are right
    assert refused.status == "refused" and refused.refusal.code == "invalid_request" and refused.estimate is None
    assert "ML = 5775 lies outside the registered frame" in refused.refusal.reason           # its number is not
    assert "ML = 5775 lies outside the registered frame" in core.value_problem(spec, soma)
    assert "ML' = 11,400 - ML" in refused.refusal.reason                # and the reason says what to do
    plan = module.plan(ask, ask.requested[0], {o.variable: o for o in ask.observed}, api.policy)
    assert isinstance(plan, core.Refusal) and plan.code == "invalid_request"
    with pytest.raises(core.ScaleError, match="outside the registered frame"):
        core.model_values(spec, soma)
    with pytest.raises(core.ScaleError, match="outside the registered frame"):
        core.ObservationStore(registry).add(soma)
    with pytest.raises(core.ForbiddenInputError, match="outside the registered frame"):
        adapter.check_inputs("q1", {"soma_position": soma, "source_label": label}, "q1_0001")
    far = dataclasses.replace(soma, value=(6400.0, 6100.0, 9000.0))     # no support range exists in these two
    with pytest.raises(core.ScaleError):
        core.ObservationStore(registry).add(far)
    with pytest.raises(core.ForbiddenInputError):
        adapter.check_inputs("q1", {"soma_position": far, "source_label": label}, "q1_0001")
    # at the midline, and mirrored (ML' = 11,400 - ML), the same neuron is answered as before
    answers = {}
    for ml in (core.MIDLINE_ML, 2 * core.MIDLINE_ML - beyond):
        ask_ok, answers[ml] = t1(ml)
        assert answers[ml].status == "inferred" and answers[ml].coverage_validated and not answers[ml].codes()
        assert adapter.check_inputs("q1", dict(zip(("soma_position", "source_label"), ask_ok.observed)), "q1_0001")
        core.ObservationStore(registry).add(ask_ok.observed[0])
    # what was returned before the fix: the model at the unmirrored number, far from the mirrored neuron's answer
    query = mods.Rows(np.array(["q"]), np.array(["q1"]),
                      {f"soma_position:{a}": np.array([v]) for a, v in zip(core.SOMA_AXES, (6400.0, 6100.0, beyond))},
                      {"source_label": np.array(["Cre"], dtype=object)})
    at_the_wrong_number = mods.engine_predict(state, tasks["T1"], query)["point"][0]
    assert np.abs(at_the_wrong_number - np.array(answers[2 * core.MIDLINE_ML - beyond].estimate.point)).max() > 0.5
    # the loader gives a training row the same verdict, by the same function
    assert "outside the registered frame" in core.bounds_problem("soma_position", core.SOMA_BOUNDS,
                                                                 (6400.0, 6100.0, beyond), core.SOMA_AXES)
    partition = core.TrainingPartition(bare_world)
    soma_rows, _ = core.read_soma_training(bare_world, partition)
    fine, _ = core.read_fine_labels_training(bare_world, partition)
    core.assemble_training_table(bare_world, partition, soma_rows, fine, None)
    neuron = sorted(soma_rows)[0]
    moved = dict(soma_rows)
    moved[neuron] = dict(soma_rows[neuron], soma=soma_rows[neuron]["soma"][:2] + (beyond,))
    with pytest.raises(core.ReconError, match="mirrored frame.*ML = 5775 lies outside the registered frame"):
        core.assemble_training_table(bare_world, partition, moved, fine, None)
    # the other statements of the frame and of the scales are numbers too: the atlas volume, a length, a fraction
    for value in ((6400.0, 6100.0, -10.0), (14000.0, 6100.0, 5000.0), (6400.0, 8200.0, 5000.0)):
        outside = dataclasses.replace(soma, value=value)
        assert "outside the registered frame" in core.value_problem(spec, outside)
        with pytest.raises(core.ForbiddenInputError):
            adapter.check_inputs("q1", {"soma_position": outside, "source_label": label}, "q1_0001")
    profile_spec = registry.variable("jiao", "regional_axon_length")
    negative = measured("regional_axon_length", (10.0, -1.0, 10.0), regions)
    assert "DMH = -1 lies outside the registered frame" in core.value_problem(profile_spec, negative)
    fraction = core.Observation("m", "epiretro", "nucleus", "n1", core.Group("e1", "experiment", "test"),
                                "methylation_global", "measured", (0.8, 1.5), "fraction in [0, 1]; fraction",
                                ("mCGFrac", "mCHFrac"), transform="identity")
    with pytest.raises(core.ScaleError, match="mCHFrac = 1.5 lies outside the registered scale"):
        core.ObservationStore(registry).add(fraction)
    core.ObservationStore(registry).add(dataclasses.replace(fraction, value=(0.8, 0.02)))
    # a value that is not one plain number per coordinate is a refusal, not an exception out of reconstruct()
    for bad in ((("1.0", "2.0", "3.0")), ((1.0, 2.0), (3.0, 4.0), (5.0, 6.0)), (True, 1.0, 2.0), (1.0, 2.0, float("inf"))):
        shaped = dataclasses.replace(soma, value=bad)
        answer = api.reconstruct(dataclasses.replace(ask, observed=(shaped, label)))["regional_axon_length"]
        assert answer.status == "refused" and answer.refusal.code == "invalid_request", bad
    with pytest.raises(ValueError, match="not a plain number"):
        core.ReconstructionRequest.from_values(registry, namespace="jiao", unit_id="q1_0001", group_id="q1",
                                               observed={"soma_position": ("6400", "6100", "5000")}, requested=[])
    # ESTIMATES. The interval of a soma estimate is cut at the midline (the fitted mixture has a tail beyond it;
    # no soma of the mirrored frame lies there), the point is not moved, and the answer is what it was otherwise.
    medial = int(np.flatnonzero(rows.cont["soma_position:ML"] > 5690.0)[0])
    lengths_um = {r: float(10 ** rows.cont[f"regional_axon_length:{r}"][medial] - 1) for r in regions}
    t2 = api.reconstruct(core.ReconstructionRequest.from_values(
        registry, namespace="jiao", unit_id="q1_0002", group_id="q1", requested=["soma_position"],
        observed={"regional_axon_length": lengths_um, "source_label": "Cre"}))["soma_position"]
    raw = mods.engine_predict(state, tasks["T2"], rows.take([medial]), level=0.9, median=True)
    assert raw["upper"][0][2] > core.MIDLINE_ML                         # the uncut interval does leave the frame
    assert t2.status == "inferred" and t2.estimate.interval["upper"][2] == core.MIDLINE_ML
    assert "cut at the bounds of the registered frame" in t2.estimate.note and t2.coverage_validated
    np.testing.assert_allclose(t2.estimate.point, raw["point"][0], rtol=1e-9)
    np.testing.assert_allclose(t2.estimate.interval["lower"], raw["lower"][0], rtol=1e-9)
    np.testing.assert_allclose(t2.estimate.interval["upper"][:2], raw["upper"][0][:2], rtol=1e-9)
    assert max(t2.estimate.point[2], t2.estimate.median[2]) <= core.MIDLINE_ML
    # ... and a module whose POINT lies beyond the midline has its answer withheld; one inside is issued
    for point, code in (((6400.0, 6100.0, beyond), "invalid_estimate"), ((6400.0, 6100.0, 5625.0), None)):
        layer = core.Reconstructor(registry)
        layer.register(FixedSomaModule(registry, point))
        answer = layer.reconstruct(request(registry, {"source_label": "A"}, ["soma_position"]))["soma_position"]
        if code is None:
            assert answer.status == "inferred" and answer.estimate.point == point
            assert core.value_problem(spec, answer.to_observation(request(registry, {}, []))) is None
        else:
            assert answer.status == "refused" and answer.refusal.code == code and answer.estimate is None
            assert "point: soma_position: ML = 5775 lies outside the registered frame" in answer.refusal.reason
    assert core.estimate_problem(spec, dataclasses.replace(t2.estimate, interval=dict(
        t2.estimate.interval, upper=(9000.0, 7000.0, beyond)))).startswith("upper interval end")


def test_K2g_a_registry_that_disagrees_with_the_layout_of_a_state_is_refused(registry, layout, tasks, planted):
    """[K2 g] The ONE conversion is decided by the module's registry; the state was fitted on the scale its own
    layout records. Nothing checked that the two agree: with a registry that models the profile on its own scale
    and a state fitted on log10(1 + length), a measured length conditioned unconverted and the state's log
    numbers were issued labelled as micrometres. A module now refuses such a registry at construction."""
    state = planted["state"]
    assert {f[4] for f in state.layout.fields if f[1] == "regional_axon_length"} == {"log10_1p"}
    key = ("jiao", "regional_axon_length")

    def changed(**fields):
        other = core.Registry.from_contract()
        other.variables[key] = dataclasses.replace(other.variables[key], **fields)
        return other

    other = changed(transform="identity")
    assert other.variable(*key).transform == "identity"
    problem = mods.layout_problem(other, layout)
    assert "regional_axon_length:ARH" in problem and "'log10_1p' in the layout" in problem and "'identity' in the registry" in problem
    with pytest.raises(core.ReconError, match="disagrees with the field layout of its states"):
        mods.JointMixtureModule(other, layout, states={n: state for n in tasks}, train=planted["rows"], tasks=tasks)
    # what it would do: the same record, 300 um per region, conditions at 300 instead of log10(301)
    record = measured("regional_axon_length", (300.0,) * 6, REGIONS)
    assert core.model_values(other.variable(*key), record).tolist() == [300.0] * 6
    np.testing.assert_allclose(core.model_values(registry.variable(*key), record), np.log10(301.0))
    # a registry in which a field's coordinates or its kind differ is refused as well
    for other in (changed(components=REGIONS[:5] + ("ZZZ",)), changed(kind="categorical", components=())):
        assert mods.layout_problem(other, layout)
        with pytest.raises(core.ReconError, match="disagrees with the field layout of its states"):
            mods.JointMixtureModule(other, layout, states={n: state for n in tasks}, train=planted["rows"], tasks=tasks)
    # the registry the layout was made from is accepted, also as a fresh object
    assert mods.layout_problem(registry, layout) is None and mods.layout_problem(core.Registry.from_contract(), layout) is None
    mods.JointMixtureModule(core.Registry.from_contract(), layout, states={n: state for n in tasks},
                            train=planted["rows"], tasks=tasks)


# ---------------------------------------------------------------------------------------------- K3


def _read_saved(path):
    with np.load(path, allow_pickle=False) as data:
        return json.loads(str(data["head"])), {k: data[k] for k in data.files if k != "head"}


def _write_saved(path, head, arrays, raw=None):
    np.savez(path, head=np.array(raw if raw is not None else core.canonical(head)), **arrays)
    return path


def _rehash(head, arrays):
    """What an editor who knows the algorithm would do: rewrite the hash stored in the same file."""
    head = {k: v for k, v in head.items() if k != "state_hash"}
    head["state_hash"] = mods.JointState._fingerprint(head, [arrays[f"a{i}"] for i in range(len(head["index"]))])
    return head


def _rejected(path, expected_hash):
    """Why JointState.load refuses the file when the caller expects the state `expected_hash`; None if it loads."""
    try:
        mods.JointState.load(path, expected_hash)
    except core.FrozenSourceChanged as error:
        return str(error)
    return None


def _slot(head, kind, name=None, block=0):
    for i, (bi, kd, nm) in enumerate(head["index"]):
        if bi == block and kd == kind and (name is None or nm == name):
            return f"a{i}", i
    raise KeyError((kind, name))


def _relabel(kind):
    def mutate(h, a):                       # every array of this kind is attached to the NEXT sorted field name
        slots = [i for i, (_, kd, _) in enumerate(h["index"]) if kd == kind]
        names = [h["index"][i][2] for i in slots]
        for i, name in zip(slots, ["!pad"] + names[:-1]):
            h["index"][i][2] = name
    return mutate


def _state_tampers(head):
    """ONE saved item changed at a time; the recorded hash is left as it was. (tag, mutate(head, arrays))."""
    AP, DV = "soma_position:AP", "soma_position:DV"
    first_region = next(k for k in head["encoding"] if k.startswith("regional"))
    group, unit = head["groups"][3], head["units"][5]
    tampers = []

    def add(tag, function):
        tampers.append((tag, function))

    def on_ranges(tag, function):
        add("ranges: " + tag, lambda h, a: function(h["ranges"]))

    on_ranges("one range widened", lambda r: r.__setitem__(AP, [r[AP][0] - 5000.0, r[AP][1] + 5000.0]))
    on_ranges("one bound moved by one ulp", lambda r: r.__setitem__(AP, [r[AP][0], float(np.nextafter(r[AP][1], np.inf))]))
    on_ranges("one range deleted", lambda r: r.__delitem__(AP))
    on_ranges("a range added for a label", lambda r: r.__setitem__("source_label", [0.0, 1.0]))
    on_ranges("a third number in a range", lambda r: r.__setitem__(AP, r[AP] + [1e9]))
    on_ranges("bounds as text", lambda r: r.__setitem__(AP, [repr(r[AP][0]), repr(r[AP][1])]))
    on_ranges("bounds swapped", lambda r: r.__setitem__(AP, r[AP][::-1]))
    on_ranges("a bound removed", lambda r: r.__setitem__(AP, [r[AP][0], None]))
    add("ranges: all removed", lambda h, a: h.__setitem__("ranges", {}))
    add("reader: support margin 0.05 -> 0.5", lambda h, a: h["reader"].__setitem__("support_margin", 0.5))
    add("reader: quantile steps 60 -> 5", lambda h, a: h["reader"].__setitem__("quantile_iterations", 5))
    add("reader: an entry removed", lambda h, a: h["reader"].__delitem__("support_margin"))
    add("format name", lambda h, a: h.__setitem__("format", "synarch_recon_state_v1"))
    add("schema name", lambda h, a: h.__setitem__("schema", "another_schema"))
    add("joint-model source hash", lambda h, a: h.__setitem__("core_sha256", "0" * 64))
    add("recorded hash replaced", lambda h, a: h.__setitem__("state_hash", "0" * 64))
    add("recorded hash removed", lambda h, a: h.__delitem__("state_hash"))
    add("an unknown entry added", lambda h, a: h.__setitem__("note", "an unknown top-level entry"))
    add("an extra array", lambda h, a: a.__setitem__("a999", np.arange(3.0)))
    add("an array removed", lambda h, a: a.__delitem__(_slot(h, "means", AP)[0]))
    # parameters: values, shapes, type, and which field each array belongs to
    add("weights: one entry + 1e-9", lambda h, a: a.__setitem__(_slot(h, "weights")[0], a[_slot(h, "weights")[0]] + np.eye(1, len(a[_slot(h, "weights")[0]]))[0] * 1e-9))
    add("means: one entry + 1e-9", lambda h, a: a.__setitem__(_slot(h, "means", AP)[0], a[_slot(h, "means", AP)[0]] + 1e-9))
    add("variances: x (1 + 1e-9)", lambda h, a: a.__setitem__(_slot(h, "variances", AP)[0], a[_slot(h, "variances", AP)[0]] * (1 + 1e-9)))
    add("tables: one table + 1e-9", lambda h, a: a.__setitem__(_slot(h, "tables")[0], a[_slot(h, "tables")[0]] + 1e-9))
    add("means: reshaped (K,) -> (K, 1)", lambda h, a: a.__setitem__(_slot(h, "means", AP)[0], a[_slot(h, "means", AP)[0]].reshape(-1, 1)))
    add("means: stored as 32-bit", lambda h, a: a.__setitem__(_slot(h, "means", AP)[0], a[_slot(h, "means", AP)[0]].astype(np.float32)))
    add("weights: components reversed", lambda h, a: a.__setitem__(_slot(h, "weights")[0], a[_slot(h, "weights")[0]][::-1].copy()))

    def swap_names(h, a):
        i, j = _slot(h, "means", AP)[1], _slot(h, "means", DV)[1]
        h["index"][i][2], h["index"][j][2] = h["index"][j][2], h["index"][i][2]

    def swap_kinds(h, a):
        i, j = _slot(h, "means", AP)[1], _slot(h, "variances", AP)[1]
        h["index"][i][1], h["index"][j][1] = "variances", "means"

    add("index: two means arrays exchange their fields", swap_names)
    add("index: means and variances of one field exchanged", swap_kinds)
    for kind in ("means", "variances", "tables"):
        add(f"index: every {kind} array attached to the next field", _relabel(kind))
    # layout: field order, kinds, transforms, names, namespace, variables
    fields = lambda h: h["layout"]["fields"]                                        # noqa: E731
    add("layout: two fields exchanged", lambda h, a: fields(h).__setitem__(slice(0, 2), fields(h)[1::-1]))
    add("layout: a kind changed", lambda h, a: fields(h)[0].__setitem__(3, "categorical"))
    add("layout: a transform changed", lambda h, a: next(f for f in fields(h) if f[4] == "log10_1p").__setitem__(4, "identity"))
    add("layout: a field's variable changed", lambda h, a: fields(h)[0].__setitem__(1, "regional_axon_length"))
    add("layout: a field renamed", lambda h, a: fields(h)[0].__setitem__(0, fields(h)[0][0] + "x"))
    add("layout: a coordinate changed", lambda h, a: fields(h)[0].__setitem__(2, "DV"))
    add("layout: namespace", lambda h, a: h["layout"].__setitem__("namespace", "epiretro"))
    add("layout: variables dropped", lambda h, a: h["layout"].__setitem__("variables", ["soma_position", "regional_axon_length"]))
    add("layout: variables reordered", lambda h, a: h["layout"].__setitem__("variables", h["layout"]["variables"][::-1]))
    # standardisation constants, categories, support
    add("standardisation: training mean + 1", lambda h, a: h["field_meta"][AP].__setitem__("mean", h["field_meta"][AP]["mean"] + 1.0))
    add("standardisation: training scale x 1.01", lambda h, a: h["field_meta"][AP].__setitem__("scale", h["field_meta"][AP]["scale"] * 1.01))
    add("categories: order", lambda h, a: h["field_meta"]["source_label"].__setitem__("categories", h["field_meta"]["source_label"]["categories"][::-1]))
    add("field: supported -> not supported", lambda h, a: h["field_meta"][AP].__setitem__("supported", False))
    add("field: block", lambda h, a: h["field_meta"][AP].__setitem__("block", "another_block"))
    add("field: type", lambda h, a: h["field_meta"][AP].__setitem__("type", "categorical"))
    # fit settings, representation, component count, block records
    for key, value in (("max_components", 9), ("seed", 123456), ("n_init", 9), ("max_iter", 9), ("tol", 0.5),
                       ("variance_floor", 0.5), ("dirichlet_alpha", 2.5)):
        add(f"settings: {key}", lambda h, a, key=key, value=value: h["settings"].__setitem__(key, value))
    add("settings: representation", lambda h, a: h["settings"].__setitem__(
        "representation", "hurdle" if h["settings"]["representation"] == "log" else "log"))
    add("encoding: one field", lambda h, a: h["encoding"].__setitem__(
        first_region, "hurdle" if h["encoding"][first_region] == "gaussian" else "gaussian"))
    add("encoding: an extra field", lambda h, a: h["encoding"].__setitem__("extra_field", "gaussian"))
    add("block: number of components", lambda h, a: h["blocks"][0].__setitem__("n_components", h["blocks"][0]["n_components"] + 1))
    add("block: number of rows", lambda h, a: h["blocks"][0].__setitem__("n_rows", h["blocks"][0]["n_rows"] + 1))
    add("block: name", lambda h, a: h["blocks"][0].__setitem__("name", "renamed"))
    add("block: field order", lambda h, a: h["blocks"][0].__setitem__("fields", h["blocks"][0]["fields"][::-1]))
    add("block: convergence flag", lambda h, a: h["blocks"][0].__setitem__("converged", not h["blocks"][0]["converged"]))
    add("block: iterations", lambda h, a: h["blocks"][0].__setitem__("iterations", h["blocks"][0]["iterations"] + 1))
    add("block: objective", lambda h, a: h["blocks"][0].__setitem__("objective", h["blocks"][0]["objective"] + 1.0))
    add("block: log-likelihood", lambda h, a: h["blocks"][0].__setitem__("log_likelihood", h["blocks"][0]["log_likelihood"] + 1.0))
    # training groups and units: as lists, in an encoding that cannot be made to collide
    position = lambda items, value: max(items.index(value), 1)                      # noqa: E731
    add("groups: one removed", lambda h, a: h["groups"].remove(group))
    add("groups: one added", lambda h, a: h["groups"].append("zz_added"))
    add("groups: one renamed", lambda h, a: h["groups"].__setitem__(h["groups"].index(group), group + "x"))
    add("groups: two IDs joined by a newline", lambda h, a: h["groups"].__setitem__(
        slice(position(h["groups"], group) - 1, position(h["groups"], group) + 1),
        ["\n".join(h["groups"][position(h["groups"], group) - 1:position(h["groups"], group) + 1])]))
    add("groups: the whole list joined", lambda h, a: h.__setitem__("groups", ["\n".join(h["groups"])]))
    add("groups: duplicates appended", lambda h, a: h["groups"].extend(h["groups"][:3]))
    add("units: one removed", lambda h, a: h["units"].remove(unit))
    add("units: one added", lambda h, a: h["units"].append("zz_added"))
    add("units: two IDs joined by a newline", lambda h, a: h["units"].__setitem__(
        slice(position(h["units"], unit) - 1, position(h["units"], unit) + 1),
        ["\n".join(h["units"][position(h["units"], unit) - 1:position(h["units"], unit) + 1])]))
    add("units: the whole list joined", lambda h, a: h.__setitem__("units", ["\n".join(h["units"])]))
    return tampers


def _twin(state, **changes):
    """The same state with single constructor arguments replaced."""
    arguments = dict(layout=state.layout, encoding=state.encoding, field_meta=state.field_meta, blocks=state.blocks,
                     settings=state.settings, units=state.units, groups=state.groups, ranges=state.ranges,
                     support_margin=state.support_margin, quantile_iterations=state.quantile_iterations)
    arguments.update(changes)
    return mods.JointState(**arguments)


def test_K3a_the_fingerprint_changes_when_support_ranges_or_any_setting_changes(layout, planted):
    """[K3 a] The support ranges decide refusals and marginalisation; the margin, the interval settings and the
    fit settings change answers too. Each of them is part of the state's fingerprint."""
    state = planted["state"]
    original = state.state_hash()
    assert _twin(state).state_hash() == original                                     # the twin itself is faithful
    AP = "soma_position:AP"
    low, high = state.ranges[AP]
    changed = {
        "a support range, widened": dict(ranges=dict(state.ranges, **{AP: (low - 5000.0, high + 5000.0)})),
        "a support range, by one ulp": dict(ranges=dict(state.ranges, **{AP: (low, float(np.nextafter(high, np.inf)))})),
        "the support margin": dict(support_margin=0.5),
        "the quantile steps": dict(quantile_iterations=30),
        "the training units": dict(units=sorted(state.units)[1:]),
        "the training groups": dict(groups=state.groups[1:]),
    }
    for key, value in (("max_components", 9), ("seed", 77), ("n_init", 9), ("max_iter", 9), ("tol", 0.5),
                       ("variance_floor", 0.5), ("dirichlet_alpha", 2.5)):
        changed[f"the setting {key}"] = dict(settings=dict(state.settings, **{key: value}))
    hashes = {what: _twin(state, **change).state_hash() for what, change in changed.items()}
    same = sorted(what for what, value in hashes.items() if value == original)
    assert not same, f"the fingerprint does not change with: {same}"
    assert len(set(hashes.values())) == len(hashes)
    # and the ranges do decide: the same input is refused by the state and answered by its widened twin
    wide = _twin(state, ranges=dict(state.ranges, **{AP: (low - 5000.0, high + 5000.0)}))
    assert state.outside_range(AP, high + 3000.0) and not wide.outside_range(AP, high + 3000.0)
    assert _twin(state, support_margin=0.5).outside_range(AP, high + 0.4 * (high - low)) is False
    assert state.outside_range(AP, high + 0.4 * (high - low)) is True
    # a state without the range of a field it was fitted on, or with a malformed one, is not a state
    for broken in ({k: v for k, v in state.ranges.items() if k != AP}, dict(state.ranges, **{AP: (low, high, 1e9)}),
                   dict(state.ranges, **{AP: (high, low)}), dict(state.ranges, source_label=(0.0, 1.0))):
        with pytest.raises(core.ReconError, match="not a valid joint state"):
            _twin(state, ranges=broken)


def test_K3b_a_saved_state_whose_ranges_were_edited_is_rejected_on_load(layout, planted, tmp_path):
    """[K3 b] The reviewed defect: ranges are saved and loaded but were not in the fingerprint, so an edited
    file verified. Now every edit of the saved ranges is rejected; and a file whose stored hash was rewritten to
    fit the edit is rejected when the caller says which state it expects."""
    state = planted["state"]
    original = state.state_hash()
    path = state.save(tmp_path / "state.npz")
    assert _rejected(path, original) is None
    with pytest.raises(TypeError):                           # a state is never loaded on the file's own word
        mods.JointState.load(path)
    assert "not the expected state" in _rejected(path, None)
    head, arrays = _read_saved(path)
    AP = "soma_position:AP"
    edits = {
        "widened": lambda r: r.__setitem__(AP, [r[AP][0] - 5000.0, r[AP][1] + 5000.0]),
        "one ulp": lambda r: r.__setitem__(AP, [r[AP][0], float(np.nextafter(r[AP][1], np.inf))]),
        "one deleted": lambda r: r.__delitem__(AP),
        "one added": lambda r: r.__setitem__("source_label", [0.0, 1.0]),
        "a third number": lambda r: r.__setitem__(AP, r[AP] + [1e9]),
        "bounds as text": lambda r: r.__setitem__(AP, [repr(r[AP][0]), repr(r[AP][1])]),
        "all removed": lambda r: r.clear(),
    }
    accepted = []
    for what, edit in edits.items():
        tampered = json.loads(json.dumps(head))
        edit(tampered["ranges"])
        if _rejected(_write_saved(tmp_path / "edited.npz", tampered, arrays), original) is None:
            accepted.append(what)
    assert not accepted, f"a saved state with edited ranges was accepted: {accepted}"
    # the editor also rewrites the hash that is stored in the same file: it is no longer the expected state
    widened = json.loads(json.dumps(head))
    edits["widened"](widened["ranges"])
    widened = _rehash(widened, arrays)
    rewritten = _write_saved(tmp_path / "rewritten.npz", widened, arrays)
    assert widened["state_hash"] != original
    assert "not the expected state" in _rejected(rewritten, original)
    other = mods.JointState.load(rewritten, widened["state_hash"])      # asked for by ITS hash it is a state: another one
    assert other.state_hash() == widened["state_hash"] != original and other.ranges[AP] != state.ranges[AP]
    # a malformed or missing range does not become a state even when the caller asks for the rewritten hash
    for what in ("a third number", "bounds as text", "one deleted", "one added", "all removed"):
        malformed = json.loads(json.dumps(head))
        edits[what](malformed["ranges"])
        malformed = _rehash(malformed, arrays)
        assert "not a state as save() writes it" in _rejected(
            _write_saved(tmp_path / "malformed.npz", malformed, arrays), malformed["state_hash"]), what


def test_K3c_every_answer_affecting_item_of_a_saved_state_is_bound(layout, tasks, planted, zero_inflated, tmp_path):
    """[K3 c] The same for EVERY item of a saved state that can change an answer, a refusal or a notice. One item
    is changed at a time and the recorded hash is left as it was; every changed file must be rejected. The list
    (it is also in the JointState docstring and in the protocol):

       1. state format and schema name                      7. fit settings: components requested, seed, restarts,
       2. SHA-256 of the unchanged joint-model source           iteration limit, tolerance, variance floor,
       3. reader settings: SUPPORT MARGIN, quantile steps       Dirichlet prior, representation
       4. layout: namespace, variables, every field         8. per block: name, fields in order, rows, number of
          (name, variable, coordinate, kind, transform)         components, convergence record
          in order                                           9. every parameter array: values, shape, type, and
       5. encoding of every continuous field                    which block / kind / field it belongs to
       6. per model field: type, supported, block,          10. SUPPORT RANGES
          standardisation mean and scale, categories        11. training units and training groups, as lists
          and their order
    Also: nothing can be changed in memory behind the fingerprint (a state is read-only and its hash is never
    cached), and a module refuses a state that is no longer the one that was registered."""
    for representation, state in (("log", planted["state"]), ("hurdle", zero_inflated["hurdle"])):
        original = state.state_hash()
        path = state.save(tmp_path / f"{representation}.npz")
        assert _rejected(path, original) is None
        head, arrays = _read_saved(path)
        accepted = []
        for tag, mutate in _state_tampers(head):
            h, a = json.loads(json.dumps(head)), {k: v.copy() for k, v in arrays.items()}
            mutate(h, a)
            if _rejected(_write_saved(tmp_path / "tampered.npz", h, a), original) is None:
                accepted.append(tag)
        assert not accepted, f"{representation}: tampered files were accepted: {accepted}"
        # JSON writes a missing value and a not-a-number alike; the stored text itself must be the canonical one
        h = json.loads(json.dumps(head))
        h["layout"]["fields"] = [[f[0], f[1], float("nan") if f[2] is None else f[2], f[3], f[4]]
                                 for f in h["layout"]["fields"]]
        raw = json.dumps(h, sort_keys=True, separators=(",", ":"), allow_nan=True) + "\n"
        assert _rejected(_write_saved(tmp_path / "nan.npz", h, arrays, raw), original) is not None
        assert _rejected(_write_saved(tmp_path / "spaced.npz", head, arrays, json.dumps(head, sort_keys=True)),
                         original) is not None
        # edits that are not states at all stay rejected when the stored hash is rewritten to fit them AND the
        # caller asks for that rewritten hash
        for tag, mutate in _state_tampers(head):
            if not tag.startswith(("index: every", "groups: two", "groups: the whole", "units: two", "units: the whole",
                                   "ranges: a third", "ranges: bounds as text", "ranges: one range deleted",
                                   "ranges: a range added", "ranges: all removed", "ranges: a bound removed",
                                   "block: number of components", "layout: a kind", "field: type",
                                   "encoding: one field", "settings: representation", "reader: an entry removed")):
                continue
            h, a = json.loads(json.dumps(head)), {k: v.copy() for k, v in arrays.items()}
            mutate(h, a)
            h = _rehash(h, a)
            assert _rejected(_write_saved(tmp_path / "rehashed.npz", h, a), h["state_hash"]) is not None, (representation, tag)
    # in memory: read-only, and the fingerprint is recomputed every time it is asked for
    state = mods.JointState.load(tmp_path / "log.npz", planted["state"].state_hash())
    before = state.state_hash()
    AP = "soma_position:AP"
    for attempt in (lambda: state.ranges.__setitem__(AP, (-1e9, 1e9)),
                    lambda: state.blocks[0]["weights"].__setitem__(slice(None), 1.0),
                    lambda: state.blocks[0]["means"][AP].__setitem__(0, 9.0),
                    lambda: state.field_meta[AP].__setitem__("mean", 0.0),
                    lambda: state.settings.__setitem__("seed", 1),
                    lambda: setattr(state, "ranges", {}),
                    lambda: setattr(state, "support_margin", 0.5),
                    lambda: state.blocks.__setitem__(0, None)):
        with pytest.raises((TypeError, ValueError, AttributeError)):
            attempt()
    assert state.state_hash() == before
    object.__setattr__(state, "ranges", dict(state.ranges, **{AP: (-1e9, 1e9)}))       # forced past the guard
    assert state.state_hash() != before                                               # the fingerprint follows
    registry = REGISTRY
    module = mods.JointMixtureModule(registry, layout, states={n: planted["state"] for n in tasks},
                                     train=planted["rows"], tasks=tasks,
                                     evidence=evidence_for(tasks, "exploratory_supported"))
    r = core.Reconstructor(registry)
    r.register(module)
    _, soma, _, label = planted_rows(1, 73, prefix="f")
    ask = request(registry, {"soma_position": soma[0], "source_label": label[0]}, ["regional_axon_length"])
    assert r.reconstruct(ask)["regional_axon_length"].status == "inferred"
    swapped = dict(module.states)
    module.states["T1"] = _twin(planted["state"], ranges=dict(planted["state"].ranges, **{AP: (-1e9, 1e9)}))
    with pytest.raises(core.FrozenSourceChanged, match="no longer the state that was registered"):
        r.reconstruct(ask)
    module.states.update(swapped)
    assert r.reconstruct(ask)["regional_axon_length"].status == "inferred"


def test_K3d_an_untampered_state_round_trips_with_identical_answers(registry, layout, tasks, planted, zero_inflated,
                                                                   tmp_path):
    """[K3 d] Saving and loading changes nothing: the same fingerprint, the same saved bytes, the same numbers,
    the same refusals - while the same file with one widened support range is not accepted."""
    rows, soma, _, label = planted_rows(25, 56, prefix="q")
    for representation, state, train in (("log", planted["state"], planted["rows"]),
                                         ("hurdle", zero_inflated["hurdle"], zero_inflated["rows"])):
        path = state.save(tmp_path / f"{representation}.npz")
        loaded = mods.JointState.load(path, expected_hash=state.state_hash())
        assert loaded.state_hash() == state.state_hash() and loaded.groups == state.groups and loaded.units == state.units
        assert dict(loaded.ranges) == dict(state.ranges) and loaded.support_margin == state.support_margin
        assert loaded.groups_sha256 == state.groups_sha256 == core.ids_sha(state.groups)
        again = loaded.save(tmp_path / f"{representation}_again.npz")
        assert _read_saved(again)[0] == _read_saved(path)[0]
        assert all(np.array_equal(a, _read_saved(path)[1][k]) for k, a in _read_saved(again)[1].items())
        for task in tasks.values():
            level = None if task.kind == "label" else 0.9
            a = mods.engine_predict(state, task, rows, level=level, median=level is not None)
            b = mods.engine_predict(loaded, task, rows, level=level, median=level is not None)
            for key in ("proba",) if task.kind == "label" else ("point", "median", "lower", "upper"):
                np.testing.assert_array_equal(a[key], b[key])
        answers = []
        for used in (state, loaded):
            module = mods.JointMixtureModule(registry, layout, states={n: used for n in tasks}, train=train, tasks=tasks,
                                             evidence=evidence_for(tasks, "exploratory_supported"))
            r = core.Reconstructor(registry)
            r.register(module)
            far = (state.ranges["soma_position:AP"][1] + 3000.0, float(soma[0][1]), float(soma[0][2]))
            results = [r.reconstruct(request(registry, observed, ["regional_axon_length"], unit=unit, group=group))
                       ["regional_axon_length"] for observed, unit, group in (
                ({"soma_position": soma[0], "source_label": label[0]}, "fresh_0001", "fresh"),
                ({"soma_position": far, "source_label": label[0]}, "fresh_0001", "fresh"),
                ({"soma_position": soma[0], "source_label": "NeverSeen"}, "fresh_0001", "fresh"),
                ({"soma_position": soma[0], "source_label": label[0]}, str(train.unit[0]), str(train.group[0])),
                ({"soma_position": soma[0], "source_label": label[0]}, f"{train.group[0]}_99999", str(train.group[0])))]
            answers.append([(a.status, sorted(a.codes()), a.coverage_validated, a.state_hash,
                             None if a.estimate is None else (a.estimate.point, a.estimate.interval["lower"],
                                                              a.estimate.interval["upper"])) for a in results])
        assert answers[0] == answers[1]
        assert [a[0] for a in answers[0]] == ["inferred", "refused", "inferred", "refused", "inferred"]
        head, arrays = _read_saved(path)                      # the control: one widened range, and it is not accepted
        name = "soma_position:AP"
        head["ranges"][name] = [head["ranges"][name][0] - 5000.0, head["ranges"][name][1] + 5000.0]
        assert _rejected(_write_saved(tmp_path / "wider.npz", head, arrays), state.state_hash()) is not None, representation


def test_run_files_around_the_states_are_bound_to_the_finished_run(planted_world, tmp_path, monkeypatch):
    """[K3, the files around a state] A state file is served only if it is the state whose fingerprint was
    recorded when the run fitted it; the evidence that decides between the engine's answer and a baseline-level
    answer, the selection and the metrics only if they are the files the completion record binds; the
    all-training states only for the run they were fitted for; and nothing when the code has changed."""
    world, _ = planted_world
    source = Path(world.manifest).parent

    def twin(name):
        shutil.copytree(source, tmp_path / name)
        return core.synthetic_world(tmp_path / name), tmp_path / name / "out"

    clean, _ = twin("clean")
    for fold in (None, 0):
        r, ctx = bench.build_reconstructor(clean, fold=fold)
        assert [m.name for m in r.modules] == ["joint_mixture", "axon_arbor_adapter"]
    # 1. a state file replaced by another state that agrees with itself (fitted on EVERY sample: it holds fold 0)
    replaced, out = twin("replaced")
    used = core.read_json(out / "selection.json")["engine_used"]["T1"]["0"]
    everything, _ = mods.fit_joint(ctx.rows, ctx.layout, components=used["components"], seed=1,
                                   budget=ctx.config["engine"], representation=used["representation"])
    file = out / "states" / bench.state_file("outer0", used["representation"], used["components"])
    everything.save(file)
    assert mods.JointState.load(file, everything.state_hash()).state_hash() != used["state_hash"]
    with pytest.raises(core.FrozenSourceChanged, match="not the expected state"):
        bench.build_reconstructor(replaced, fold=0)
    # 2. the evidence edited by hand: every task 'supported', coverage 'validated'
    edited, out = twin("evidence")
    evidence = core.read_json(out / "evidence.json")
    for entry in evidence["tasks"].values():
        entry["status"] = "exploratory_supported"
        entry["coverage"]["engine"]["validated"] = True
    core.write_json(out / "evidence.json", evidence)
    with pytest.raises(core.FrozenSourceChanged, match="evidence.json"):
        bench.build_reconstructor(edited, fold=0)
    with pytest.raises(core.FrozenSourceChanged, match="evidence.json"):
        bench.evidence_objects(edited)
    # 3. the selection edited so that it points at the replaced state
    pointed, out = twin("selection")
    selection = core.read_json(out / "selection.json")
    selection["engine_used"]["T1"]["0"]["state_hash"] = everything.state_hash()
    core.write_json(out / "selection.json", selection)
    with pytest.raises(core.FrozenSourceChanged, match="selection.json"):
        bench.build_reconstructor(pointed, fold=0)
    # 4. the metrics, or the declaration, changed after the run
    for name in ("metrics.json", "pre_run_declaration.json"):
        changed, out = twin("changed_" + name)
        record = core.read_json(out / name)
        record["note"] = "edited after the run"
        core.write_json(out / name, record, pretty=name.startswith("pre_run"))
        with pytest.raises(core.ReconError):
            bench.build_reconstructor(changed, fold=0)
        with pytest.raises(core.FrozenSourceChanged, match=name):          # nor is a report rendered from them
            bench.render_report(changed)
    # 5. serving.json is a copy of the record that the completion record holds: a copy that says it belongs to
    #    another finished run is an edited file of this run, and then nothing is served, for any fold
    foreign, out = twin("serving")
    serving = core.read_json(out / "serving.json")
    serving["run_complete_sha256"] = "0" * 64
    core.write_json(out / "serving.json", serving)
    for fold in (None, 0):
        with pytest.raises(core.FrozenSourceChanged, match="serving.json"):
            bench.build_reconstructor(foreign, fold=fold)
    (out / "serving.json").unlink()                                   # without the copy the bound record still serves
    for fold in (None, 0):
        bench.build_reconstructor(foreign, fold=fold)
    # 6. the code is no longer the frozen code
    real_hashes = bench.code_hashes
    monkeypatch.setattr(bench, "code_hashes", lambda: dict(real_hashes(), **{"analysis/synarch_recon_modules.py": "x"}))
    with pytest.raises(core.FrozenSourceChanged, match="code_sha256"):
        bench.build_reconstructor(clean, fold=0)


# ---------------------------------------------------------------------------------------------- K3, second round
#
# The fingerprint of a state held under attack. What did not: for the all-training states the EXPECTED
# fingerprint came from serving.json, a file that nothing bound.


def _edit_saved_range(file):
    """What an editor who knows the algorithm does to a saved state: widen ONE support range and rewrite the
    hash stored in the same file. Returns (the state's fingerprint before, after)."""
    head, arrays = _read_saved(file)
    before = head["state_hash"]
    head["ranges"]["soma_position:AP"][1] += 1e6
    head = _rehash(head, arrays)
    _write_saved(file, head, arrays)
    return before, head["state_hash"]


def test_K3e_the_all_training_states_are_bound_by_the_completion_record(planted_world, tmp_path, monkeypatch):
    """[K3 e] The all-training states answer for units outside the training partition. Their expected
    fingerprints were read from serving.json, which was written after run_complete.json and bound by nothing: a
    state with one widened support range was served once the copy of its fingerprint in serving.json was
    replaced too, while the documented anchor (the SHA-256 of run_complete.json) and `status` still verified.
    A request that had been refused as out of support was then answered.

    Now the run fits these states before it scores, and its completion record HOLDS their record. The API takes
    every expected fingerprint from behind that one anchor; serving.json is a copy that must be exactly the copy;
    by hand the states can be restored, never changed."""
    world, _ = planted_world
    source = Path(world.manifest).parent
    quiet = lambda *a: None                                                          # noqa: E731

    def twin(name):
        shutil.copytree(source, tmp_path / name)
        return core.synthetic_world(tmp_path / name), tmp_path / name / "out"

    def bound(out):
        return core.read_json(out / "run_complete.json")[bench.SERVING_KEY]

    def state_path(out, entry):
        return out / "states" / bench.state_file("serving", entry["representation"], entry["components"])

    def answers(w):
        """T1 for a new neuron: in the middle of the training somas, and 20% of the AP range above them."""
        r, ctx = bench.build_reconstructor(w)
        soma = ctx.table.soma
        low, high = soma[:, 0].min(), soma[:, 0].max()
        found = []
        for ap in (np.median(soma[:, 0]), high + 0.2 * (high - low)):
            ask = core.ReconstructionRequest.from_values(
                ctx.registry, namespace="jiao", unit_id="999999_001", group_id="999999", requested=["regional_axon_length"],
                observed={"soma_position": (float(ap), float(np.median(soma[:, 1])), float(np.median(soma[:, 2]))),
                          "source_label": str(ctx.table.label[0])})
            answer = r.reconstruct(ask)["regional_axon_length"]
            found.append(answer.status if answer.refusal is None else answer.refusal.code)
        return found

    clean, out = twin("clean")
    anchor = core.sha_file(out / "run_complete.json")
    assert answers(clean) == ["inferred", "out_of_support"]
    assert bench.status(clean)["run_files_match_the_completion_record"] is True
    assert bench.status(clean)["run_complete_sha256"] == anchor          # `status` shows the anchor
    assert core.read_json(out / "serving.json") == dict(bound(out), run_complete_sha256=anchor)

    # (a) the reported break: one saved item of the all-training state edited, the hash stored in that file
    #     rewritten, and the copy of the fingerprint in serving.json replaced. Nothing the anchor covers changes.
    edited, out = twin("edited")
    entry = bound(out)["tasks"]["T1"]
    old, new = _edit_saved_range(state_path(out, entry))
    assert old == entry["state_hash"] != new
    with pytest.raises(core.FrozenSourceChanged, match="not the expected state"):
        bench.build_reconstructor(edited)                               # the state file alone: refused before, too
    text = (out / "serving.json").read_text()
    assert old in text
    (out / "serving.json").write_text(text.replace(old, new))
    assert core.sha_file(out / "run_complete.json") == anchor
    for fold in (None, 0):
        with pytest.raises(core.FrozenSourceChanged, match="serving.json is not the record"):
            bench.build_reconstructor(edited, fold=fold)
    with pytest.raises(core.FrozenSourceChanged, match="serving.json is not the record"):
        bench.verify_run(edited)
    with pytest.raises(core.FrozenSourceChanged):
        bench.fit_serving(edited, log=quiet)
    report = bench.status(edited)["run_files_match_the_completion_record"]
    assert report is not True and "serving.json" in report             # `status` says it as well
    (out / "serving.json").unlink()                                     # and without the copy: the expected state is
    with pytest.raises(core.FrozenSourceChanged, match="not the expected state"):    # the bound one, not the file's
        bench.build_reconstructor(edited)

    # (b) no state edited by hand: serving.json points task T1 at ANOTHER genuine all-training state, of a
    #     setting that the selection of this run did not choose (fitted by the layer's own job)
    swapped, out = twin("swapped")
    ctx = bench.build_context(swapped)
    entry = bound(out)["tasks"]["T1"]
    other = next((r, k) for r, k in bench.engine_grid(ctx.config) if (r, k) != (entry["representation"], entry["components"]))
    fit = bench.job_serving(ctx, other[0], other[1], out / "states")
    assert fit["state_hash"] and fit["state_hash"] != entry["state_hash"]
    serving = core.read_json(out / "serving.json")
    serving["tasks"]["T1"].update(representation=other[0], components=other[1], state_hash=fit["state_hash"],
                                  training_groups_sha256=fit["training_groups_sha256"])
    core.write_json(out / "serving.json", serving)
    assert core.sha_file(out / "run_complete.json") == anchor
    with pytest.raises(core.FrozenSourceChanged, match="serving.json is not the record"):
        bench.build_reconstructor(swapped)

    # (c) WHERE THE GUARANTEE ENDS, shown rather than claimed: an editor who also rewrites the completion record
    #     makes every file agree again. The edited state is then served, under ANOTHER SHA-256 of
    #     run_complete.json: the anchor that `run` printed and the report states is what detects it.
    rewritten, out = twin("rewritten")
    old, new = _edit_saved_range(state_path(out, bound(out)["tasks"]["T1"]))
    complete = core.read_json(out / "run_complete.json")
    for task_entry in complete[bench.SERVING_KEY]["tasks"].values():
        if task_entry["state_hash"] == old:
            task_entry["state_hash"] = new
    core.write_json(out / "run_complete.json", complete)
    core.write_json(out / "serving.json", bench.serving_file_content(rewritten))
    assert answers(rewritten) == ["inferred", "inferred"]              # the widened range answers what was refused
    assert core.sha_file(out / "run_complete.json") != anchor
    #     ... and even then the all-training state of a task must be the setting that the bound selection chose
    pointed, out = twin("pointed")
    fit = bench.job_serving(bench.build_context(pointed), other[0], other[1], out / "states")
    complete = core.read_json(out / "run_complete.json")
    complete[bench.SERVING_KEY]["tasks"]["T1"].update(representation=other[0], components=other[1],
                                                      state_hash=fit["state_hash"],
                                                      training_groups_sha256=fit["training_groups_sha256"])
    core.write_json(out / "run_complete.json", complete)
    core.write_json(out / "serving.json", bench.serving_file_content(pointed))
    bench.verify_run(pointed)
    with pytest.raises(core.FrozenSourceChanged, match="not the engine setting that the selection"):
        bench.build_reconstructor(pointed)

    # (d) by hand the all-training states are RESTORED from the completion record, never changed
    lost, out = twin("lost")
    files = sorted({state_path(out, e) for e in bound(out)["tasks"].values()})
    for file in files:
        file.unlink()
    (out / "serving.json").unlink()
    with pytest.raises(core.ReconError, match="is missing"):
        bench.build_reconstructor(lost)
    bench.build_reconstructor(lost, fold=0)                             # the out-of-fold states are not affected
    real_fit = mods.fit_joint

    def reseeded(rows, layout, *, components, seed, budget, representation="log"):
        return real_fit(rows, layout, components=components, seed=seed + 1, budget=budget, representation=representation)

    monkeypatch.setattr(mods, "fit_joint", reseeded)                    # a re-fit that is not the recorded state
    with pytest.raises(core.FrozenSourceChanged, match="not the state that the finished run recorded"):
        bench.fit_serving(lost, log=quiet)
    assert not any(f.exists() for f in files) and not (out / "serving.json").exists()
    assert not (out / "states" / "restore_in_progress").exists()
    monkeypatch.setattr(mods, "fit_joint", real_fit)
    restored = bench.fit_serving(lost, log=quiet)
    assert restored == dict(bound(out), run_complete_sha256=anchor) and all(f.exists() for f in files)
    assert answers(lost) == ["inferred", "out_of_support"] and core.sha_file(out / "run_complete.json") == anchor
    with pytest.raises(core.ReconError, match="only after the benchmark has completed"):
        bench.fit_serving(bench.make_synthetic_world(tmp_path / "never_run", "null", size="tiny", seed=9), log=quiet)

    # (e) the order of the run: the all-training states exist before the first metric is computed, and the
    #     completion record and its copy are written after it, also without the steps that follow a run
    fresh = bench.make_synthetic_world(tmp_path / "order", "planted", size="tiny", seed=6)
    bench.freeze(fresh, "synthetic", log=quiet)
    seen, real_score = {}, bench.score

    def watched(w, log=quiet):
        folder = Path(w.out)
        seen["states"] = sorted(p.name for p in (folder / "states").glob("serving_*.npz"))
        seen["files"] = {name: (folder / name).exists() for name in ("metrics.json", "run_complete.json", "serving.json")}
        return real_score(w, log=log)

    monkeypatch.setattr(bench, "score", watched)
    bench.run(fresh, log=quiet, post=False)
    folder = Path(fresh.out)
    assert seen["files"] == {"metrics.json": False, "run_complete.json": False, "serving.json": False}
    assert seen["states"] == sorted({state_path(folder, e).name for e in bound(folder)["tasks"].values()})
    assert (folder / "serving.json").exists() and not (folder / "api_demo.json").exists()
    bench.verify_run(fresh)
    bench.build_reconstructor(fresh)


def test_K3f_state_files_and_derived_lookups(registry, layout, tasks, planted, tmp_path):
    """Smaller defects next to K3. A stored member that is not an array is a changed file, not a crash; a state
    is written to exactly the path it is given; the lookups that are derived from a state cannot be edited in
    place, and a state in which one of them was replaced has no fingerprint and answers no request."""
    import zipfile

    state = planted["state"]
    path = state.save(tmp_path / "state.npz")
    raw = tmp_path / "raw_member.npz"
    with zipfile.ZipFile(path) as source, zipfile.ZipFile(raw, "w") as target:
        for member in source.namelist():
            if member == "a0.npy":
                target.writestr("a0", b"not an array")
            else:
                target.writestr(member, source.read(member))
    assert "64-bit floating point" in _rejected(raw, state.state_hash())
    odd = state.save(tmp_path / "state_without_suffix")
    assert odd == tmp_path / "state_without_suffix" and odd.exists()
    assert not (tmp_path / "state_without_suffix.npz").exists()
    assert mods.JointState.load(odd, state.state_hash()).state_hash() == state.state_hash()
    loaded = mods.JointState.load(path, state.state_hash())
    with pytest.raises(TypeError):
        loaded._model_block["soma_position:AP"] = 1
    with pytest.raises(TypeError):
        loaded._kind["soma_position:AP"] = "hurdle"
    assert loaded.state_hash() == state.state_hash()
    for name, value in (("units_folded", frozenset()), ("groups_folded", frozenset()), ("groups_sha256", "0" * 64),
                        ("_kind", dict(loaded._kind)), ("_model_block", {}),
                        ("units", frozenset({"x_1"})), ("groups", ("zz",))):     # the lists they are derived from
        forced = mods.JointState.load(path, state.state_hash())
        object.__setattr__(forced, name, value)                         # past the read-only guard
        with pytest.raises(core.FrozenSourceChanged, match="derived lookup"):
            forced.state_hash()
    # a training unit spelled in another case is recognised through the folded identifiers; with them emptied
    # behind the guard the state no longer answers at all
    used = mods.JointState.load(path, state.state_hash())
    module = mods.JointMixtureModule(registry, layout, states={n: used for n in tasks}, train=planted["rows"],
                                     tasks=tasks, evidence=evidence_for(tasks, "exploratory_supported"))
    r = core.Reconstructor(registry)
    r.register(module)
    rows = planted["rows"]
    observed = {"soma_position": [float(rows.cont[f"soma_position:{a}"][0]) for a in core.SOMA_AXES],
                "source_label": rows.cat["source_label"][0]}
    shouted = request(registry, observed, ["regional_axon_length"], unit=str(rows.unit[0]).upper(),
                      group=str(rows.group[0]).upper())
    assert r.reconstruct(shouted)["regional_axon_length"].refusal.code == "leakage_guard"
    object.__setattr__(used, "units_folded", frozenset())
    with pytest.raises(core.FrozenSourceChanged, match="derived lookup"):
        r.reconstruct(shouted)


# ----------------------------------------------------------------------------------------------
# the comparison is fair to a model that learned nothing, and the status needs a real gain
# ----------------------------------------------------------------------------------------------


WIDE = core.JIAO_REGIONS[:30]


def _no_signal_rows(seed, groups=40):
    """Soma, label and profile mutually independent; very uneven numbers of neurons per brain sample (half of
    the samples hold one to three neurons, the other half 40 to 160); 30 regions."""
    rng = np.random.default_rng(seed)
    sizes = np.concatenate([rng.integers(1, 4, groups // 2), rng.integers(40, 160, groups - groups // 2)])
    group = np.repeat([f"n{seed}g{i:03d}" for i in range(groups)], sizes)
    n = len(group)
    soma = np.array([6500.0, 6000.0, 4800.0]) + rng.normal(0, [600, 400, 350], (n, 3))
    present = rng.random((n, len(WIDE))) < rng.uniform(0.2, 0.8, len(WIDE))
    y = np.where(present, np.maximum(rng.uniform(2.0, 4.0, len(WIDE)) + rng.normal(0, 0.5, (n, len(WIDE))), 0.05), 0.0)
    cont = {f"soma_position:{a}": soma[:, j] for j, a in enumerate(core.SOMA_AXES)}
    cont.update({f"regional_axon_length:{r}": y[:, j] for j, r in enumerate(WIDE)})
    labels = np.array([LABELS[i % 3] for i in range(groups)], dtype=object)[np.repeat(np.arange(groups), sizes)]
    return mods.Rows(np.array([f"{g}_{i:05d}" for i, g in enumerate(group)]), group, cont, {"source_label": labels}), y


def test_a_one_component_engine_is_the_neuron_weighted_mean_and_does_not_beat_the_tuned_mean(registry):
    """With one component the engine has learned nothing from its inputs: its prediction IS the neuron-weighted
    training mean. Against a MEAN baseline that may only weight brain samples equally, such an engine is 'better,
    with an interval excluding zero' in most no-signal worlds: with uneven sample sizes the sample-weighted mean
    is simply the noisier estimate of the same thing. So the weighting of MEAN (and IDENTITY) is a declared
    setting chosen on inner folds, like k and lambda; against that MEAN the one-component engine wins nothing."""
    layout = mods.jiao_layout(registry, WIDE)
    task = mods.build_tasks(layout, registry, WIDE[:10])["T1"]
    grid, config = bench.mean_grid(bench.CONFIG), bench.TINY_CONFIG
    assert grid == [{"weighting": "sample"}, {"weighting": "neuron"}]
    beats_sample_weighted = beats_tuned = chose_neuron = 0
    worlds = 6
    for seed in range(worlds):
        rows, y = _no_signal_rows(seed)
        number = np.array([int(g[-3:]) for g in rows.group])
        fold = number % 5
        engine, by_sample, tuned = (np.zeros_like(y) for _ in range(3))
        for k in range(5):
            outer_train, test = np.flatnonzero(fold != k), rows.take(np.flatnonzero(fold == k))
            # the selection the benchmark makes: mean inner-fold primary metric, ties to the earlier setting
            inner, scored = (number[outer_train] // 5) % 3, []
            cells = []
            for i in range(3):
                bank = mods.TaskBaselines(task, rows.take(outer_train[inner != i]))
                held = rows.take(outer_train[inner == i])
                cells.append([bench.primary_score(task, p, held.matrix(task.target_fields), held.group)
                              for p in bank.predict_many("MEAN", held, grid)])
            for j, params in enumerate(grid):
                scored.append((params, (float(np.mean([c[j][0] for c in cells])), 0.0)))
            chosen = bench.pick(task, scored)[0]
            chose_neuron += chosen == {"weighting": "neuron"}
            bank = mods.TaskBaselines(task, rows.take(outer_train))
            by_sample[fold == k] = bank.predict("MEAN", test, {"weighting": "sample"})
            tuned[fold == k] = bank.predict("MEAN", test, chosen)
            engine[fold == k] = bank.predict("MEAN", test, {"weighting": "neuron"})
            if seed == 0:                                    # the fitted one-component engine IS that mean
                state, _ = mods.fit_joint(rows.take(outer_train), layout, components=1, seed=1,
                                          budget={"n_init": 1, "max_iter": 50, "tol": 1e-7})
                np.testing.assert_allclose(mods.engine_predict(state, task, test)["point"], engine[fold == k],
                                           rtol=0, atol=1e-9)
        assert np.abs(by_sample - engine).max() > 1e-3                               # the other mean is another
        per = lambda p: core.by_sample(((p - y) ** 2).mean(axis=1), rows.group)      # noqa: E731
        beats_sample_weighted += bench.compare("mse", per(engine), per(by_sample), seed, config)["engine_better_interval"]
        against_tuned = bench.compare("mse", per(engine), per(tuned), seed, config)
        beats_tuned += against_tuned["engine_better_interval"]
        assert not against_tuned["engine_better_clear"]
    assert beats_sample_weighted >= worlds - 2      # the defect of a MEAN that can only weight samples equally
    assert beats_tuned == 0 and chose_neuron == 5 * worlds
    assert mods.TUNED == ("MEAN", "IDENTITY", "KNN", "RIDGE") and mods.WEIGHTINGS == ("sample", "neuron")
    assert [bench.tuned_grid(bench.CONFIG, task, m, True) for m in ("MEAN", "IDENTITY")] == [grid, grid]


def test_the_weighting_is_selected_on_inner_folds_and_reported(null_world, planted_world):
    for world, metrics in (null_world, planted_world):
        out = Path(world.out)
        selection = core.read_json(out / "selection.json")
        declaration = core.read_json(out / "pre_run_declaration.json")
        assert declaration["config"]["mean"] == {"weightings": ["sample", "neuron"]}
        assert declaration["config"]["min_relative_gain"] == 0.02 and declaration["protocol_sha256"] is None
        assert declaration["state_reader"] == {"state_format": mods.STATE_FORMAT, "support_margin": 0.05,
                                               "quantile_iterations": 60}
        for name in ("T1", "T2", "T3", "T4"):
            for fold in map(str, range(5)):
                chosen = selection["comparators"][name][fold]
                assert chosen["MEAN"]["selected"] in bench.mean_grid(bench.CONFIG) and len(chosen["MEAN"]["grid"]) == 2
                assert ("IDENTITY" in chosen) == (name != "T3")
            assert set(selection["serving"][name]["comparators"]) == set(selection["comparators"][name]["0"])
        evidence = core.read_json(out / "evidence.json")
        for name, entry in evidence["tasks"].items():
            best = entry["best_comparator"]
            if best in mods.TUNED:
                assert entry["best_comparator_params"] == selection["serving"][name]["comparators"][best]
    # in the no-signal world the engine does not clearly beat MEAN on any task
    world, metrics = null_world
    for name, result in metrics["tasks"].items():
        assert not result["comparisons"][result["primary_metric"]]["MEAN"]["engine_better_clear"], name


def test_the_status_needs_every_interval_and_the_minimum_gain(registry, layout, tasks, tmp_path):
    """exploratory_supported = better than EVERY comparator, every interval excluding zero, AND a relative gain of
    at least the declared minimum over each. A consistent advantage of half a percent is not enough."""
    task, rng = tasks["T2"], np.random.default_rng(0)
    n = 400
    sample = np.array([f"s{i % 40:02d}" for i in range(n)])
    truth = rng.normal(6000, 500, (n, 3))
    error = rng.normal(0, 100, (n, 3))

    def data(factor):
        arrays = {"truth": truth, "sample": sample, "label": np.array(["A"] * n), "classes": np.array(["A"]),
                  "ENGINE": truth + error, "ENGINE_MEDIAN": truth + error, "ENGINE_lower": truth - 300,
                  "ENGINE_upper": truth + 300, "unknown_label": np.zeros(n, bool), "outside_range": np.zeros(n, bool)}
        for comparator in mods.COMPARATORS:
            arrays[comparator] = truth + error * factor                    # every comparator: the same, slightly worse
        for comparator in ("MEAN", "IDENTITY"):
            arrays[comparator + "_lower"], arrays[comparator + "_upper"] = truth - 400, truth + 400
        return arrays

    config = json.loads(json.dumps(bench.TINY_CONFIG))
    world = dataclasses.replace(core.real_world(), synthetic=True, out=tmp_path)
    core.write_json(tmp_path / "selection.json", {"serving": {"T2": {"comparators": {}}}})
    outcomes = {}
    for factor in (1.005, 1.05):
        result = bench.score_task("T2", task, data(factor), config, "salt")
        verdict = result["verdict"]
        assert verdict["all_comparators"] == verdict["all_comparators_interval"] == "every"        # in both cases
        assert result["best_comparator"]["relative_gain_of_engine"] == pytest.approx(1 - 1 / factor, rel=1e-6)
        for name in ("metrics.json", "evidence.json"):
            (tmp_path / name).unlink(missing_ok=True)
        core.write_json(tmp_path / "metrics.json", {"tasks": {"T2": result}})
        evidence = bench.write_evidence(world, {"tasks": {"T2": result}}, config)["tasks"]["T2"]
        outcomes[factor] = (verdict["all_comparators_clear"], verdict["always_reported_baselines_clear"], evidence["status"])
        assert evidence["min_relative_gain"] == 0.02 and "min_relative_gain" in evidence["status_rule"]
        assert evidence["relative_gain_over_best_comparator"] == result["best_comparator"]["relative_gain_of_engine"]
    assert outcomes[1.005] == ("none", "none", "no_better_than_baseline")          # 0.5%: every interval, no status
    assert outcomes[1.05] == ("every", "every", "exploratory_supported")           # 4.8%: supported
    config["min_relative_gain"] = 0.0                                              # without the minimum, 0.5% would pass
    assert bench.score_task("T2", task, data(1.005), config, "salt")["verdict"]["all_comparators_clear"] == "every"


def test_partly_linked_evidence_is_reported_exactly(registry, layout, tasks):
    """One soma coordinate was measured together with the profile; the other two only on other units. The answer
    conditions on the first, says so, lists the other two as ignored under their own name, and equals the
    conditional that never saw them."""
    rows, soma, y, label = planted_rows(600, 23)
    elsewhere = np.arange(600) >= 400                            # these units measured DV and ML, and nothing else
    for name in rows.cont:
        only_there = name in ("soma_position:DV", "soma_position:ML")
        rows.cont[name] = np.where(elsewhere if not only_there else ~elsewhere, np.nan, rows.cont[name])
    rows.cat["source_label"] = np.array([None if e else l for l, e in zip(label, elsewhere)], dtype=object)
    state, _ = mods.fit_joint(rows, layout, components=3, seed=1, budget={"n_init": 2, "max_iter": 50, "tol": 1e-6})
    assert len(state.blocks) == 2
    module = mods.JointMixtureModule(registry, layout, states={n: state for n in tasks}, train=rows, tasks=tasks,
                                     evidence=evidence_for(tasks, "exploratory_supported"))
    r = core.Reconstructor(registry)
    r.register(module)
    answer = r.reconstruct(request(registry, {"soma_position": soma[0], "source_label": label[0]},
                                   ["regional_axon_length"]))["regional_axon_length"]
    assert answer.status == "inferred" and "evidence_not_linked" in answer.codes() and not answer.coverage_validated
    assert dict(answer.evidence_used) == {"soma_position": "measured", "source_label": "measured"}
    assert [v for v, _ in answer.evidence_ignored] == ["soma_position[DV, ML]"]
    assert not set(dict(answer.evidence_used)) & {v for v, _ in answer.evidence_ignored}
    query = query_rows(soma[0], label[0])
    np.testing.assert_allclose(answer.estimate.point, mods.engine_predict(state, tasks["T1"], query)["point"][0], rtol=1e-12)
    for hidden in ("soma_position:DV", "soma_position:ML"):                # the answer does not depend on them
        query.cont[hidden] = np.array([np.nan])
    np.testing.assert_allclose(answer.estimate.point, mods.engine_predict(state, tasks["T1"], query)["point"][0], rtol=1e-12)


# ----------------------------------------------------------------------------------------------
# the real training partition: identifiers, counts and file handling only
# ----------------------------------------------------------------------------------------------

REAL = core.real_world()
needs_real = pytest.mark.skipif(not all(Path(p).exists() for p in (REAL.manifest, REAL.soma_csv, REAL.soma_xlsx)),
                                reason="the Jiao tables or the cohort manifest are not on disk")


@needs_real
def test_real_training_partition_inputs_and_a_real_canary():
    partition = core.TrainingPartition(REAL)
    assert (len(partition.train_samples), len(partition.train_neurons)) == (101, 4448)
    table = core.load_training_table(REAL, with_outcomes=False, partition=partition)      # inputs only
    assert table.length is None and len(table) == 4448 and len(set(table.sample.tolist())) == 101
    assert set(table.sample.tolist()) == set(partition.train_samples)
    assert len(set(table.label.tolist())) == 16 and "Wt(orexin)" in table.label
    assert (table.soma[:, 2] <= core.MIDLINE_ML).all() and np.isfinite(table.soma).all()
    assert all(len(set(table.label[table.sample == s].tolist())) == 1 for s in partition.train_samples)
    assert table.audit["soma_csv"]["sealed_rows_dropped_unparsed"] == 7180 - 4448
    with open(REAL.manifest, newline="") as handle:                                       # identifiers only
        sealed = next(row["animal"] for row in csv.DictReader(handle, delimiter="\t") if row["part"] == "validation")
    assert not partition.is_training(sealed)
    with pytest.raises(core.SealedRowError):
        partition.require_training([sealed], "canary with a real validation sample ID")
    with pytest.raises(core.SealedRowError):
        core.assemble_training_table(REAL, partition, {f"{sealed}_001": {"full_id": f"{sealed}_001_x_y", "geno": "x",
                                                                         "location": "y", "soma": (1.0, 2.0, 3.0)}},
                                     {}, None)
    splits = bench.make_splits(table, bench.CONFIG)
    assert sorted(s for held in splits["outer"].values() for s in held) == sorted(partition.train_samples)
    assert sorted(len(held) for held in splits["outer"].values()) == [20, 20, 20, 20, 21]


@needs_real
def test_real_training_neurons_load_as_measured_records():
    """Loads the training rows of the regional table to build records; nothing is fitted or cross-tabulated."""
    if not Path(REAL.regional_csv).exists():
        pytest.skip("regional table not on disk")
    registry = core.Registry.from_contract()
    table = core.load_training_table(REAL)
    store = core.jiao_records(REAL, table, registry)
    assert len(store) == 4 * 4448 and len(store.units("jiao", "neuron")) == 4448
    assert table.length.shape == (4448, 49) and np.isfinite(table.length).all() and (table.length >= 0).all()
    assert table.audit["regional_csv"]["sealed_rows_dropped_unparsed"] == 7180 - 4448
    view = store.unit("jiao", "neuron", str(table.neuron[0]))
    assert {r.status for r in view.values()} == {"measured"} and len(view) == 4
    assert view["regional_axon_length"].source.sha256 == core.REAL_PINS["regional_csv"]
    assert view["soma_position"].group.kind == "brain_sample"
    assert view["soma_position"].group.group_id == str(table.neuron[0]).split("_")[0]
