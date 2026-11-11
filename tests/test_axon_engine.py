"""Checks of the axon prediction engine (analysis/axon_engine.py).

Synthetic worlds are generated in temporary directories. A few checks read the 37 TRAINING traces and
the training rows of the regional table that are on disk (stage 'develop'); they never touch a
validation or test trace or regional row, and they are skipped when the raw data are absent.
The independent referee (analysis/prediction_independent_evaluator.py) is imported read-only by the
engine and really called; it is never edited or patched here (one test wraps its functions to count calls).

    python3 -I -m pytest tests/test_axon_engine.py -q
"""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Every process that calls the referee runs with exactly these four thread variables (late item L6): the referee's
# audited distance cache records them, and the engine refuses the referee without them. They are set BEFORE numpy is
# imported, because numpy's thread pools are fixed when it loads (round-2 finding: the engine checks the pools).
THREADS = {v: "1" for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
os.environ.update(THREADS)

import numpy as np  # noqa: E402
import pytest  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import axon_engine as E  # noqa: E402

T, G, Geo = E.T, E.G, E.Geo
REAL = (ROOT / "data/raw/jiao2025/swc/MANIFEST.tsv").exists() and (ROOT / "reports/axon_dataset_manifest.tsv").exists()
real_only = pytest.mark.skipif(not REAL, reason="raw training data are not on disk")


# ------------------------------------------------------------------------------------------
# fixtures
# ------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def pipeline(tmp_path_factory):
    """The whole pipeline, once, on a tiny synthetic world (canaries included), with the candidate family that
    has no fitted regional state: the test goes through the referee's four functions."""
    directory = tmp_path_factory.mktemp("pipeline") / "informative"
    summary = E.run_synthetic("informative", directory, size="tiny", workers=1, log=None, family="no_fitted_state")
    return directory, summary


@pytest.fixture(scope="module")
def pipeline_all(tmp_path_factory):
    """The same tiny world with the whole candidate family (fitted-state candidates included); no development pass, no canary."""
    directory = tmp_path_factory.mktemp("pipeline_all") / "informative"
    summary = E.run_synthetic("informative", directory, size="tiny", workers=1, log=None, family="all", develop=False, canary=False)
    return directory, summary


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """A second tiny world with its library, all-pairs cache and regional states built; nothing forecast yet."""
    directory = tmp_path_factory.mktemp("built") / "dense"
    E.make_synthetic_world(directory, "dense", size="tiny")
    engine = E.Engine(E.synthetic_world(directory), "develop", workers=1, log=None)
    engine.build_library()
    engine.build_cache()
    engine.states.build_many(list(engine.states.plan))
    return directory


@pytest.fixture(scope="module")
def frozen_world(built, tmp_path_factory):
    """The `built` world taken through validation and the freeze, with its test forecasts committed; no test truth opened."""
    directory = clone(built, tmp_path_factory.mktemp("frozen") / "w")
    validated(directory)
    engine_of(directory, "validate").freeze()
    engine_of(directory, "test").predict()
    return directory


def clone(source, target):
    shutil.copytree(source, target)
    return target


def engine_of(directory, stage, **kw):
    return E.Engine(E.synthetic_world(directory), stage, workers=1, log=None, **kw)


def validated(directory, family="all"):
    engine_of(directory, "validate", family=family).predict()
    engine_of(directory, "validate").admit()
    return engine_of(directory, "validate").score()


def arrays(manifest):
    return {(p["query_id"], p["model"]): p["arrays_sha256"] for p in manifest["predictions"]}


def replace_traces(directory, part):
    """Other trees in every trace file of one part (and their hashes in the trace manifest); nothing else is touched."""
    path = directory / "data/traces/MANIFEST.tsv"
    lines = path.read_text().splitlines()
    changed = 0
    for i, line in enumerate(lines[1:], 1):
        fields = line.split("\t")
        if fields[0].startswith(part + "/"):
            trace = directory / "data/traces" / fields[0]
            rows = trace.read_text().splitlines()
            for k in range(len(rows) - 1, 0, -1):                       # move the last axon node by 40 um
                f = rows[k].split()
                if not rows[k].startswith("#") and f[1] in ("0", "2"):
                    f[2] = f"{float(f[2]) + 40.0:.2f}"
                    rows[k] = " ".join(f)
                    break
            trace.write_text("\n".join(rows) + "\n")
            fields[2], fields[3] = str(trace.stat().st_size), E.sha_file(trace)
            lines[i] = "\t".join(fields)
            changed += 1
    path.write_text("\n".join(lines) + "\n")
    return changed


def load(path):
    return json.loads(Path(path).read_text())


def write_swc(path, nodes):
    """nodes: (id, type, AP, DV, ML, parent)."""
    path.write_text("# test\n" + "".join(f"{i} {t} {a:.2f} {d:.2f} {m:.2f} 1.0 {p}\n" for i, t, a, d, m, p in nodes))
    return path


def small_tree(code=2, *, soma=(6500.0, 6000.0, 5000.0), main=70, with_dendrite=True, markers=True):
    """Soma; a main axon path along +AP with a marker inside it; a side branch; a marker at the tip;
    a dendrite with its own marker."""
    a, d, m = soma
    nodes = [(1, 1, a, d, m, -1)]
    for k in range(1, main + 1):
        kind = 10 if (markers and k == main // 2) else code
        nodes.append((k + 1, kind, a + 20.0 * k, d, m, k))
    last_main = main + 1
    start = 21
    for k in range(1, 11):
        nodes.append((last_main + k, code, a + 20.0 * 20, d + 20.0 * k, m, start if k == 1 else last_main + k - 1))
    n = last_main + 10
    if markers:
        n += 1
        nodes.append((n, 12, a + 20.0 * main + 3.0, d, m, last_main))          # terminal marker on the axon
    if with_dendrite:
        first = n + 1
        for k in range(3):
            nodes.append((first + k, 3, a, d - 20.0 * (k + 1), m, 1 if k == 0 else first + k - 1))
        n = first + 2
        if markers:
            n += 1
            nodes.append((n, 11, a, d - 85.0, m, n - 1))                        # marker inside the dendrite
    return nodes


def ingest(path, **kw):
    options = {"frame": "test_true_hemisphere", "hemisphere": "Left", "recon_type": "Axon_and_dendrite",
               "table_soma": (6500.0, 6000.0, 5000.0)}
    options.update(kw)
    return E.ingest_trace(path, **options)


# ------------------------------------------------------------------------------------------
# the frozen prototype is pinned and imported unchanged
# ------------------------------------------------------------------------------------------


def test_frozen_pins_are_verified_and_a_changed_file_is_detected(tmp_path):
    pins = E.verify_frozen_sources()
    assert set(E.PINNED_SOURCES) <= set(pins["source_sha256"])
    for relative in list(pins["source_sha256"]) + [E.FREEZE_RECORD]:
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / relative, target)
    assert E.verify_frozen_sources(tmp_path)["source_sha256"] == pins["source_sha256"]
    with (tmp_path / "analysis/brm_geometry.py").open("a") as handle:
        handle.write("\n# one changed byte\n")
    with pytest.raises(E.FrozenSourceChanged, match="brm_geometry"):
        E.verify_frozen_sources(tmp_path)


def test_engine_uses_the_frozen_modules_five_policies_and_the_contracts_four_copies():
    assert Path(T.__file__).name == "brm_template.py" and Path(G.__file__).name == "d110_geometry.py"
    assert tuple(T.POLICIES) == ("baseline", "uniform_anchor_500", "conditional_anchor_500", "uniform_anchor_1000", "conditional_anchor_1000")
    assert set(E.FROZEN_POLICY_OF.values()) == set(T.POLICIES)
    assert len(E.ENGINE_CANDIDATES) <= 7 and len(E.METHODS) == len(set(E.METHODS)) == 16      # E1 to E7; 7 controls + 9 candidates
    assert E.OFFICIAL_COPIES == ("nearest", "medoid", "peptide_nearest", "peptide_medoid") == tuple(E.referee().COPIES)
    assert set(E.CONTROLS) == set(E.OFFICIAL_COPIES) | {"C1_medoid_copy", "C4_uniform_attach_500", "C5_uniform_attach_1000"}
    assert "sampled panels" in E.CONFIG["donor_pools"] and E.CONFIG["donor_pools"]["sampled panels"] == "none"
    assert not hasattr(E, "select_panel")
    # the referee is the file in the repository, loaded read-only, and its hashes are recorded
    assert Path(E.referee().__file__) == ROOT / "analysis/prediction_independent_evaluator.py"
    assert E.referee().source_hashes() == E.referee_hashes()


# ------------------------------------------------------------------------------------------
# ingestion: raw adapter, axon rule, markers, frame, resampling, inclusion rules, no clipping
# ------------------------------------------------------------------------------------------


def test_axon_rule_gives_the_same_arbor_for_type2_and_type0_coding_and_matches_the_frozen_extractor(tmp_path):
    two, zero = write_swc(tmp_path / "two.swc", small_tree(2)), write_swc(tmp_path / "zero.swc", small_tree(0))
    log2, a2 = ingest(two)
    log0, a0 = ingest(zero)
    assert log2["status"] == log0["status"] == "included"
    assert (log2["axon_coding"], log0["axon_coding"]) == (2, 0)
    assert E.same_arbor(a2, a0)
    # the frozen extractor (type 2 only) and this engine in frozen mode: bit-identical
    frozen = G.resample(G.extract_explicit_axon(two, False), 50.0)
    assert E.same_arbor(frozen, E.frozen_replay_arbor(two, False))
    # plan rule: the marker inside the path is axon by inheritance (no connector node), the terminal marker is axon cable (3 um more)
    assert log2["terminal_marker_twig_nodes"] == 1 and log2["terminal_marker_twig_um"] == pytest.approx(3.0)
    assert a2.cable() == pytest.approx(frozen.cable() + 3.0)
    # nothing of the dendrite or of its marker is in the axon
    assert a2.xyz[:, 1].min() >= 6000.0 - 1e-9 and log2["connector_nodes"] == 0


def test_a_marker_inherits_the_type_of_its_nearest_non_marker_ancestor(tmp_path):
    """Markers 10, 11, 12: below an axon node they are axon, below a dendrite node dendrite, also through a chain of markers."""
    nodes = small_tree(0, markers=False)
    n = max(i for i, *_ in nodes)
    axon_tip, dendrite_tip = 71, n                                     # end of the main axon path; end of the dendrite
    nodes += [(n + 1, 10, 6500.0 + 1410.0, 6000.0, 5000.0, axon_tip), (n + 2, 11, 6500.0 + 1420.0, 6000.0, 5000.0, n + 1),
              (n + 3, 12, 6500.0, 5920.0, 5000.0, dendrite_tip), (n + 4, 10, 6500.0, 5900.0, 5000.0, n + 3)]
    log, arbor = ingest(write_swc(tmp_path / "m.swc", nodes))
    plain_log, plain = ingest(write_swc(tmp_path / "p.swc", small_tree(0, markers=False)))
    assert log["status"] == "included" and log["terminal_marker_twig_nodes"] == 2                # the two below the axon tip
    assert log["axon_cable_um"] == pytest.approx(plain_log["axon_cable_um"] + 20.0)
    assert arbor.xyz[:, 0].max() == pytest.approx(6500.0 + 1420.0)
    assert arbor.xyz[:, 1].min() >= 6000.0 - 1e-9                      # the markers below the dendrite stayed dendrite
    assert log["dendrite_nodes"] == plain_log["dendrite_nodes"] == 3   # raw type-3 nodes; markers are counted by inheritance


def test_an_axon_that_leaves_a_dendrite_keeps_the_dendritic_stem_as_logged_connector_geometry(tmp_path):
    """Fixed now, before any validation trace exists: the stem is part of the scored tree, it is not axon cable
    (rule 6 counts axon only), and it is logged. Applied alike to donors and truth."""
    a, d, m = 6500.0, 6000.0, 5000.0
    nodes = [(1, 1, a, d, m, -1), (2, 3, a, d - 30.0, m, 1), (3, 3, a, d - 60.0, m, 2)]         # a dendritic stem of 60 um
    for k in range(1, 61):
        nodes.append((3 + k, 0, a + 20.0 * k, d - 60.0, m, 3 if k == 1 else 2 + k))             # 1,200 um of axon from the stem
    nodes.append((64, 3, a, d - 90.0, m, 3))                                                     # the dendrite continues
    log, arbor = ingest(write_swc(tmp_path / "x.swc", nodes))
    assert log["status"] == "included" and log["axon_cable_um"] == pytest.approx(1200.0)
    assert log["dendrite_connector_nodes"] == 2 and log["connector_um"] == pytest.approx(60.0)
    # stem kept, the rest of the dendrite not; the 50 um samples cut the right-angle corner by one chord (40, 10) instead of 50 um
    assert arbor.cable() == pytest.approx(1260.0 - (50.0 - np.hypot(40.0, 10.0))) and arbor.xyz[:, 1].min() == pytest.approx(d - 60.0)
    short = nodes[:3] + nodes[3:3 + 49] + [(64, 3, a, d - 90.0, m, 3)]                           # 980 um of axon + 60 um of stem
    log, arbor = ingest(write_swc(tmp_path / "s.swc", short))
    assert arbor is None and log["rules_violated"] == ["6"]                                      # the stem does not count as axon
    summary = E.connector_summary([{"status": "included", "connector_nodes": 2, "connector_um": 60.0, "dendrite_connector_nodes": 2}])
    assert summary["neurons_whose_axon_leaves_a_dendrite"] == 1 and summary["connector_um_total"] == 60.0


def test_atlas_frame_constants_and_the_mirror_changes_ml_only(tmp_path):
    assert (E.MIDLINE_ML, E.MIRROR_ML, E.STEP) == (5700.0, 11400.0, 50.0) and E.BOX.tolist() == [13200.0, 8000.0, 11400.0]
    assert E.CONFIG["inclusion"]["soma_tolerance_um"] == 5.0 and E.CONFIG["inclusion"]["min_axon_cable_um"] == 1000.0
    left = small_tree(0, soma=(6500.0, 6000.0, 5100.0))
    right = [(i, t, a, d, 11400.0 - m, p) for i, t, a, d, m, p in left]
    table = (6500.0, 6000.0, 5100.0)
    log_l, arbor_l = ingest(write_swc(tmp_path / "l.swc", left), table_soma=table)
    log_r, arbor_r = ingest(write_swc(tmp_path / "r.swc", right), hemisphere="Right", table_soma=table)
    assert not log_l["mirrored_to_left"] and log_r["mirrored_to_left"] and log_r["file_hemisphere"] == "Right"
    assert np.array_equal(arbor_l.xyz[:, :2], arbor_r.xyz[:, :2])                     # AP and DV untouched
    assert np.allclose(arbor_l.xyz[:, 2], arbor_r.xyz[:, 2]) and arbor_r.xyz[0, 2] == pytest.approx(5100.0)
    # a soma exactly on the midline is left; just beyond it is right
    assert ingest(write_swc(tmp_path / "m.swc", small_tree(0, soma=(6500.0, 6000.0, 5700.0))), table_soma=(6500.0, 6000.0, 5700.0))[0]["file_hemisphere"] == "Left"
    beyond = ingest(write_swc(tmp_path / "b.swc", small_tree(0, soma=(6500.0, 6000.0, 5700.5))), table_soma=(6500.0, 6000.0, 5699.5), hemisphere="Right")[0]
    assert beyond["file_hemisphere"] == "Right" and beyond["trace_soma_um"][2] == pytest.approx(5699.5)


def test_resampling_keeps_forks_and_tips_and_no_edge_exceeds_50_um(tmp_path):
    """The canonical arbor on the engine's own path: frozen branch-to-branch resampling at 50 um."""
    nodes = small_tree(0, markers=False, with_dendrite=False)
    log, arbor = ingest(write_swc(tmp_path / "x.swc", nodes), recon_type="Axon_only")
    edges = np.linalg.norm(arbor.xyz[1:] - arbor.xyz[arbor.parent[1:]], axis=1)
    assert edges.max() <= 50.0 + 1e-9 and arbor.parent[0] == -1 and (arbor.parent[1:] < np.arange(1, len(arbor.xyz))).all()
    produced = {tuple(np.round(v, 6)) for v in arbor.xyz.tolist()}
    fork, tip_main, tip_side = (6500.0 + 400.0, 6000.0, 5000.0), (6500.0 + 1400.0, 6000.0, 5000.0), (6500.0 + 400.0, 6200.0, 5000.0)
    assert {fork, tip_main, tip_side} <= produced and tuple(arbor.xyz[0]) == (6500.0, 6000.0, 5000.0)
    assert arbor.branches() == 1 and int((np.bincount(arbor.parent[1:], minlength=len(arbor.xyz)) == 0).sum()) == 2
    assert arbor.cable() == pytest.approx(1400.0 + 200.0) == pytest.approx(log["axon_cable_um"])   # straight runs: no cable lost
    assert log["canonical_nodes"] == len(arbor.xyz) and log["arbor_sha256"] == E.arbor_sha(arbor)
    # the same arbor, shifted, seen by the referee's own checks: a valid canonical tree
    E.referee().tree_check(E.as_tree(arbor), [6500.0, 6000.0, 5000.0], 5.0)


def test_true_hemisphere_copy_mirrored_equals_the_mirrored_left_copy(tmp_path):
    left = small_tree(0, soma=(6500.0, 6000.0, 5000.0))
    right = [(i, t, a, d, 11400.0 - m, p) for i, t, a, d, m, p in left]
    _, mirrored = ingest(write_swc(tmp_path / "m.swc", left), frame="x_mirrored_left", hemisphere="Right")
    log, true_frame = ingest(write_swc(tmp_path / "t.swc", right), frame="x_true_hemisphere", hemisphere="Right")
    assert log["mirrored_to_left"] and log["file_hemisphere"] == "Right" and E.same_arbor(mirrored, true_frame)
    assert true_frame.xyz[0, 2] == 5000.0


@pytest.mark.parametrize("case,rule", [
    ("missing_parent", "1"), ("two_roots", "1"), ("cycle", "1"), ("root_not_soma", "1"), ("not_a_table", "1"), ("soma_far", "2"),
    ("hemisphere", "3"), ("bad_type", "4"), ("both_codes", "4"), ("no_dendrite", "5"), ("axon_only_with_dendrite", "5"),
    ("short_axon", "6")])
def test_every_inclusion_rule_excludes_and_is_logged(tmp_path, case, rule):
    nodes, kw = small_tree(2), {}
    if case == "missing_parent":
        nodes[30] = nodes[30][:5] + (987654,)
    elif case == "two_roots":
        nodes[30] = nodes[30][:5] + (-1,)
    elif case == "cycle":
        nodes[30] = nodes[30][:5] + (nodes[32][0],)
    elif case == "root_not_soma":
        nodes[0] = (1, 2) + nodes[0][2:]
    elif case == "soma_far":
        kw["table_soma"] = (6500.0, 6000.0, 5006.0)
    elif case == "hemisphere":
        kw["hemisphere"] = "Right"                      # the file's soma is on the left
    elif case == "bad_type":
        nodes[5] = (nodes[5][0], 7) + nodes[5][2:]
    elif case == "both_codes":
        nodes[5] = (nodes[5][0], 0) + nodes[5][2:]
    elif case == "no_dendrite":
        nodes = small_tree(2, with_dendrite=False)
    elif case == "axon_only_with_dendrite":
        kw["recon_type"] = "Axon_only"
    elif case == "short_axon":
        nodes = small_tree(2, main=30)
    path = write_swc(tmp_path / "x.swc", nodes)
    if case == "not_a_table":
        path.write_text("this is not an SWC file\n1 2 three\n")
    log, arbor = ingest(path, **kw)
    assert arbor is None and log["status"] == "excluded" and rule in log["rules_violated"] and rule in log["details"]


def test_axon_only_without_dendrite_is_included(tmp_path):
    log, arbor = ingest(write_swc(tmp_path / "x.swc", small_tree(0, with_dendrite=False)), recon_type="Axon_only")
    assert log["status"] == "included" and arbor is not None and log["dendrite_nodes"] == 0


def test_no_clipping_a_path_that_leaves_the_atlas_and_returns_keeps_exactly_its_measured_cable(tmp_path):
    """Plan amendment 1, first counterexample of the review: exit and re-entry at different points."""
    nodes = [(1, 1, 12800.0, 4000.0, 5000.0, -1)]
    path = [(12800.0 + 50.0 * k, 4000.0) for k in range(1, 11)] + [(13300.0, 4050.0), (13300.0, 4100.0)] + \
           [(13300.0 - 50.0 * k, 4100.0) for k in range(1, 30)]
    for k, (ap, dv) in enumerate(path, 2):
        nodes.append((k, 2, ap, dv, 5000.0, k - 1))
    log, arbor = ingest(write_swc(tmp_path / "x.swc", nodes), recon_type="Axon_only", table_soma=(12800.0, 4000.0, 5000.0))
    measured = 500.0 + 100.0 + 1450.0                                  # out along +AP, across, back along -AP
    assert log["status"] == "included" and arbor.cable() == pytest.approx(measured) == pytest.approx(log["axon_cable_um"])
    assert arbor.xyz[:, 0].max() == pytest.approx(13300.0)             # the cable outside the box is still there
    outside = log["outside_atlas_box"]
    assert outside["outside_um"] == pytest.approx(100.0 + 100.0 + 100.0) and outside["nodes_outside"] == 5
    assert {(13300.0, 4050.0), (13300.0, 4100.0)} <= {(round(x, 6), round(y, 6)) for x, y, _ in arbor.xyz.tolist()}


def test_no_clipping_an_edge_between_two_outside_points_that_crosses_the_atlas_is_kept_whole(tmp_path):
    """Plan amendment 1, second counterexample: both ends outside, 13,200 um of the edge inside."""
    nodes = [(1, 1, 100.0, 4000.0, 5000.0, -1), (2, 2, 50.0, 4000.0, 5000.0, 1), (3, 2, -50.0, 4000.0, 5000.0, 2),
             (4, 2, 13250.0, 4000.0, 5000.0, 3)]
    log, arbor = ingest(write_swc(tmp_path / "x.swc", nodes), recon_type="Axon_only", table_soma=(100.0, 4000.0, 5000.0))
    assert log["status"] == "included" and arbor.cable() == pytest.approx(50.0 + 100.0 + 13300.0)
    assert log["outside_atlas_box"]["outside_um"] == pytest.approx(50.0 + 50.0 + 50.0)      # exact, edge by edge
    assert arbor.xyz[:, 0].min() == pytest.approx(-50.0) and arbor.xyz[:, 0].max() == pytest.approx(13250.0)
    assert not hasattr(E, "clip_to_box") and "amendment 1" in E.RULE_TEXT["7"]
    inside_only = E.outside_box_cable(np.array([[1.0, 1.0, 1.0], [50.0, 1.0, 1.0]]), np.array([-1, 0]))
    assert inside_only == {"nodes_outside": 0, "outside_um": 0.0, "soma_outside_box": False}


@pytest.mark.parametrize("case", ["exit_then_reentry", "outside_endpoints_cross_box"])
def test_the_reviews_two_clipping_counterexamples_return_exactly_the_measured_cable(tmp_path, case):
    """The coordinates of data/interim/axon_engine_codex_review/clipping_counterexamples.json (the draft's clipping
    returned 1,400 um and 100 um). Amendment 1: ingestion returns the measured tree, whole."""
    stored = {"exit_then_reentry": ([[13000.0, 1000.0, 1000.0], [14000.0, 1000.0, 1000.0], [14000.0, 2000.0, 1000.0], [13000.0, 2000.0, 1000.0]], 3000.0, 2600.0),
              "outside_endpoints_cross_box": ([[100.0, 1000.0, 1000.0], [-100.0, 1000.0, 1000.0], [13300.0, 1000.0, 1000.0]], 13600.0, 300.0)}
    xyz, measured, outside = stored[case]
    record = ROOT / "data/interim/axon_engine_codex_review/clipping_counterexamples.json"
    if record.exists():                                                 # the review's own file, read only
        saved = next(c for c in load(record)["cases"] if c["case"] == case)
        assert saved["input_xyz_um"] == xyz and saved["returned_cable_um"] != measured
    nodes = [(i + 1, 1 if i == 0 else 2, a, d, m, -1 if i == 0 else i) for i, (a, d, m) in enumerate(xyz)]
    log, arbor = ingest(write_swc(tmp_path / "x.swc", nodes), recon_type="Axon_only", table_soma=tuple(xyz[0]))
    assert log["status"] == "included" and arbor.cable() == pytest.approx(measured) and log["axon_cable_um"] == pytest.approx(measured)
    assert log["outside_atlas_box"]["outside_um"] == pytest.approx(outside)
    produced = {tuple(v) for v in arbor.xyz.tolist()}
    assert all(tuple(v) in produced for v in xyz)                       # every measured vertex is still there, unmoved


def test_forecasts_and_truth_share_the_full_coordinate_domain(pipeline):
    """A translated or attached template is never clipped either: same domain for truth, controls and candidates."""
    directory, summary = pipeline
    engine = engine_of(directory, "test")
    lib = engine.library
    outside = [n for n in lib["ids"] if lib["entries"][n]["outside_atlas_box"]["outside_um"] > 0]
    assert outside and summary["library"]["cable_outside_atlas_box"]["neurons_with_cable_outside"] == len(outside)
    manifest = load(Path(engine.world.out) / "test" / "prediction_manifest.json")
    copies = [p for p in manifest["predictions"] if p["selected_donor"] in outside and p["model"] in E.COPY_BASELINES]
    assert copies
    for p in copies[:5]:
        donor = engine.donor(p["selected_donor"])["arbor"]
        assert (donor.xyz > E.BOX).any() or (donor.xyz < 0).any()
        forecast = G.load_arbor(Path(engine.world.out) / "test" / p["path"])
        assert len(forecast.xyz) == len(donor.xyz) and forecast.cable() == pytest.approx(donor.cable())
    metrics = load(Path(engine.world.out) / "test" / "metrics.json")
    assert set(metrics["cable_outside_atlas_box"]) == {"neurons", "neurons_with_cable_outside", "total_um", "max_um"}
    declaration = load(Path(engine.world.out) / "test" / "referee" / "declaration.json")
    assert declaration["geometry_contract"]["clipping"].startswith("none (plan amendment 1)")


@real_only
def test_real_type2_traces_give_the_frozen_canonical_arbor_bit_for_bit(tmp_path):
    engine = E.Engine(E.real_world(out=tmp_path, legacy_only=True), "develop", workers=1, log=None)
    result = E.check_frozen_equivalence(engine)
    assert result["files"] >= 16 and result["all_bit_identical"]
    # what the plan's rules add on these files is small and is measured, not hidden
    assert result["max_terminal_marker_twig_um"] < 50.0 and result["max_plan_vs_frozen_distance_um"] < 5.0


@real_only
def test_real_portal_copy_mirrored_equals_zenodo_copy(tmp_path):
    engine = E.Engine(E.real_world(out=tmp_path, legacy_only=True), "develop", workers=1, log=None)
    result = E.check_two_copies(engine)
    assert result["neurons"] == 4 and result["all_equal"]


@real_only
def test_real_replay_reproduces_the_earlier_medoid_figure_exactly():
    result = E.replay_earlier_medoid()
    assert result["earlier"]["distance_um"] == pytest.approx(922.4851270347687)
    assert result["frozen_extractor"]["equals_earlier_exactly"] and result["engine_extraction_frozen_mode"]["equals_earlier_exactly"]


# ------------------------------------------------------------------------------------------
# metrics
# ------------------------------------------------------------------------------------------


def test_scoring_is_the_frozen_metric_and_the_referees_metric_averaged_within_animal_then_over_animals():
    a = G.resample(G.Arbor(np.array([[0.0, 0, 0], [600.0, 0, 0], [600.0, 400.0, 0]]), np.array([-1, 0, 1]), np.array([1, 2, 2])), 50.0)
    b = a.shifted(np.array([0.0, 150.0, 0.0]))
    frozen = G.geometry_metrics(a, b)
    assert E.score_arbor(a, b) == {k: frozen[k] for k in E.KEPT_METRICS} == E.score_arbor_checked(a, b)
    theirs = E.referee().G.geometry_metrics(E.as_tree(a), E.as_tree(b))
    assert all(frozen[k] == theirs[k] for k in E.KEPT_METRICS)
    row = lambda animal, value: {"animal": animal, "metrics": {"m": {k: value for k in E.KEPT_METRICS}}}
    result = E.aggregate([row("A", 10.0), row("A", 30.0), row("B", 100.0)], ["m"])
    assert result["animals"] == 2 and result["neurons"] == 3
    assert result["summary"]["m"]["distance_um"] == pytest.approx(60.0)          # not the neuron mean 46.7


def test_a_disagreement_between_the_two_metric_implementations_stops_the_scoring(monkeypatch):
    a = G.resample(G.Arbor(np.array([[0.0, 0, 0], [600.0, 0, 0]]), np.array([-1, 0]), np.array([1, 2])), 50.0)
    Ref = E.referee()
    real = Ref.G.geometry_metrics
    monkeypatch.setattr(Ref.G, "geometry_metrics", lambda p, t: {**real(p, t), "distance_um": real(p, t)["distance_um"] + 1e-9})
    with pytest.raises(E.EngineError, match="disagree"):
        E.score_arbor_checked(a, a.shifted(np.array([0.0, 10.0, 0.0])))


# ------------------------------------------------------------------------------------------
# sealing
# ------------------------------------------------------------------------------------------


def garble_sealed_rows(directory):
    world = E.synthetic_world(directory)
    cohort = E.Cohort(world)
    held = {cohort.rows[n].full_id for part in ("validation", "test") for n in cohort.part_neurons[part]}
    lines = Path(world.regional_table).read_text().splitlines()
    for i in range(1, len(lines)):
        full_id = lines[i].split(",", 1)[0]
        if full_id in held:
            lines[i] = full_id + ",THIS,ROW,IS,SEALED,AND,MUST,NEVER,BE,PARSED"
    Path(world.regional_table).write_text("\n".join(lines) + "\n")
    return len(held)


def test_sealed_regional_rows_are_never_parsed(built, tmp_path):
    directory = clone(built, tmp_path / "w")
    world = E.synthetic_world(directory)
    cohort = E.Cohort(world)
    before, _ = E.load_training_profiles(world, cohort)
    held = garble_sealed_rows(directory)
    after, audit = E.load_training_profiles(world, cohort)
    assert audit["sealed_rows_skipped_unparsed"] == held and audit["training_rows_loaded"] == len(cohort.part_neurons["train"])
    assert set(after) == {cohort.rows[n].full_id for n in cohort.part_neurons["train"]}
    assert all(np.array_equal(before[k], after[k]) for k in before)
    # a damaged TRAINING row is not skipped: it is an error
    lines = Path(world.regional_table).read_text().splitlines()
    target = next(i for i in range(1, len(lines)) if lines[i].split(",", 1)[0] in after)
    lines[target] = lines[target].split(",", 1)[0] + ",not,numbers"
    Path(world.regional_table).write_text("\n".join(lines) + "\n")
    with pytest.raises(Exception):
        E.load_training_profiles(world, cohort)


@real_only
def test_real_regional_loader_keeps_training_animals_only():
    world = E.real_world()
    cohort = E.Cohort(world)
    profiles, audit = E.load_training_profiles(world, cohort)
    assert audit["training_rows_loaded"] == 4448 and audit["sealed_rows_skipped_unparsed"] == 2732
    assert {full_id.split("_", 1)[0] for full_id in profiles} == set(cohort.part_animals["train"])
    assert not {full_id.split("_", 1)[0] for full_id in profiles} & cohort.heldout_animals
    assert len(cohort.never_heldout) == 33 and not set(cohort.never_heldout) & cohort.heldout_animals    # 32 inspected + the Li-shared one


def test_stage_guard_on_traces(built):
    world = E.synthetic_world(built)
    cohort = E.Cohort(world)
    train, validation, test = (cohort.part_neurons[p][0] for p in ("train", "validation", "test"))
    develop = E.TraceStore(world, cohort, "develop")
    assert develop.open(train, "donor")["path"].endswith(".swc") and develop.open(train, "truth")
    for neuron in (validation, test):
        for purpose in ("donor", "truth"):
            with pytest.raises(E.SealError):
                develop.open(neuron, purpose)
    with pytest.raises(E.SealError):
        develop.unlock_truth("abc")
    validate = E.TraceStore(world, cohort, "validate")
    with pytest.raises(E.SealError, match="after the forecasts are committed"):
        validate.open(validation, "truth")
    validate.unlock_truth("commitment")
    assert validate.open(validation, "truth")
    with pytest.raises(E.SealError):
        validate.open(test, "truth")                    # the validation stage never opens a test trace
    with pytest.raises(E.SealError):
        validate.open(validation, "donor")              # a held-out neuron is never a donor
    with pytest.raises(E.SealError):
        validate.open(train, "truth")


def test_a_held_out_trace_listed_in_the_training_directory_is_refused(built, tmp_path):
    world = E.synthetic_world(built)
    cohort = E.Cohort(world)
    neuron = cohort.part_neurons["validation"][0]
    legacy = tmp_path / "MANIFEST.tsv"
    legacy.write_text("file\tneuron\tbytes\tsha256\tsource\tframe\n" + f"x/{neuron}.swc\t{neuron}\t1\tabc\tx\tportal_true_hemisphere\n")
    planted = E.World(world.name, world.kind, world.manifest, world.local_manifest, world.soma_table, world.regional_table,
                      world.trace_manifests + ((str(legacy), str(tmp_path), "legacy_training"),), world.out)
    with pytest.raises(E.SealError, match="training-only directory"):
        E.TraceStore(planted, cohort, "develop")


@real_only
def test_real_validation_and_test_refuse_and_state_the_blocker(tmp_path):
    for stage in ("validate", "test"):
        engine = E.Engine(E.real_world(out=tmp_path / stage, legacy_only=True), stage, workers=1, log=None)
        for action in (engine.predict, engine.admit, engine.build_library):
            with pytest.raises(E.StageBlocked, match="have not been fetched"):
                action()
        assert not (tmp_path / stage / stage).exists()                    # nothing was written
    status = E.Engine(E.real_world(out=tmp_path / "s", legacy_only=True), "develop", workers=1, log=None).status()
    assert status["traces_on_disk"]["validation"]["have"] == 0 and status["traces_on_disk"]["test"]["have"] == 0
    assert status["traces_on_disk"]["train"]["have"] == 37 and status["stages"]["develop"] == "data present"


# ------------------------------------------------------------------------------------------
# regional states
# ------------------------------------------------------------------------------------------


def test_regional_states_hold_training_animals_only_and_satisfy_the_frozen_contract(built):
    engine = engine_of(built, "validate")
    cohort = engine.cohort
    manifest = engine.states.manifest()["states"]
    assert set(manifest) == set(engine.states.plan) and "all_training" in manifest
    for state_id, entry in manifest.items():
        state = engine.states.get(state_id)
        assert entry["validation_or_test_animals_in_state"] == []
        assert not set(state.proof["train_animals"]) & cohort.heldout_animals
        assert entry["animals_sha256"] == E.ids_sha(sorted(state.proof["train_animals"]))
        assert T._digest(state.model) == state.proof["model_state_sha256"] == entry["model_state_sha256"]
        hashes = engine.states.file_hashes(state_id)
        assert hashes["state_file_sha256"] == E.sha_file(engine.states.dir / f"{state_id}.json") == entry["model_state_sha256"]
    # the frozen safety check accepts a validation query against the all-training state ...
    row = cohort.rows[cohort.part_neurons["validation"][0]]
    soma = np.array(row.soma)
    Geo.assert_query_safe(T.Query(row.neuron, row.animal, row.full_id, soma, soma, row.tag), engine.states.get("all_training"))
    # ... and a training query only against the out-of-fold state that excludes its animal
    train_row = cohort.rows[cohort.part_neurons["train"][0]]
    own = engine.states.get(engine.states.state_id_for(train_row.neuron))
    assert train_row.animal not in own.proof["train_animals"]
    q = T.Query(train_row.neuron, train_row.animal, train_row.full_id, np.array(train_row.soma), np.array(train_row.soma), train_row.tag)
    Geo.assert_query_safe(q, own)
    other = next(s for s in manifest if s.startswith("oof_") and s != engine.states.state_id_for(train_row.neuron))
    with pytest.raises(ValueError):
        Geo.assert_query_safe(q, engine.states.get(other))


def test_fitting_a_state_on_a_held_out_animal_is_refused_and_states_ignore_whom_they_serve(built):
    engine = engine_of(built, "develop")
    cohort = engine.cohort
    train = cohort.part_animals["train"]
    for part in ("validation", "test"):
        with pytest.raises(E.LeakError, match="outside the training part"):
            E.fit_regional_state("planted", cohort, engine.profiles, train + [cohort.part_animals[part][0]],
                                 cohort.part_neurons["test"][:2])
    with pytest.raises(E.LeakError, match="its own fitting animals"):
        E.fit_regional_state("planted", cohort, engine.profiles, train, cohort.part_neurons["train"][:1])
    a, _ = E.fit_regional_state("same", cohort, engine.profiles, train, cohort.part_neurons["validation"])
    b, _ = E.fit_regional_state("same", cohort, engine.profiles, train, cohort.part_neurons["test"][:3])
    assert a.proof["model_state_sha256"] == b.proof["model_state_sha256"] and a.proof["test_ids"] != b.proof["test_ids"]
    # a fresh fit reproduces the stored state exactly (the freeze relies on this)
    assert engine.states.refit_sha("all_training") == engine.states.file_hashes("all_training")["model_state_sha256"]


def test_the_frozen_baseline_and_uniform_policies_do_not_depend_on_the_numbers_of_the_regional_state(built):
    """C1, C4 and C5 consult the state for membership only: with another valid state their arrays are the same, so
    their forecasts declare no upstream fit. The conditional policies do depend on it."""
    engine = engine_of(built, "validate")
    engine.prepare_prediction()
    neuron = engine.cohort.part_neurons["validation"][0]
    row = engine.cohort.rows[neuron]
    soma = np.array(row.soma, dtype=float)
    query = T.Query(row.neuron, row.animal, row.full_id, soma.copy(), soma.copy(), row.tag)
    pool = engine.pool_donors(engine.scope(None))
    state = engine.states.get("all_training")
    model = json.loads(json.dumps(state.model))
    block = model["blocks"][0]
    block["weights"] = list(reversed(block["weights"]))                     # another mixture: a different numeric state
    for field in block["means"]:
        block["means"][field] = [v + 0.37 for v in block["means"][field]]
    other = Geo.RegionalState(model, {**state.proof, "model_state_sha256": T._digest(model)})
    changed = False
    with engine.memo.frozen_memo():
        for policy in T.POLICIES:
            one, _ = T.predict_policy(query, pool, state, policy)
            two, info = T.predict_policy(query, pool, other, policy)
            if policy.startswith("conditional_"):
                changed |= not E.same_arbor(one, two) or info["weights"] != T.predict_policy(query, pool, state, policy)[1]["weights"]
            else:
                assert E.same_arbor(one, two), policy
    assert changed                                                           # the conditional weights did see the other state
    assert set(E.USES_REGIONAL_STATE) == {"F1_conditional_attach_500", "F2_conditional_attach_1000",
                                          "E4_local_kernel_regional", "E5_local_kernel_label_regional"}
    # E7's fitted state is the general reconstruction layer's, not this regional state (tests/test_synarch_axon_bridge.py)
    assert set(E.USES_FITTED_STATE) == set(E.USES_REGIONAL_STATE) | {"E7_local_kernel_inferred_profile"}


# ------------------------------------------------------------------------------------------
# planted-leak canaries
# ------------------------------------------------------------------------------------------


def test_planted_leak_a_held_out_animal_in_the_state_makes_the_engine_refuse(built, tmp_path):
    directory = clone(built, tmp_path / "w")
    engine = engine_of(directory, "validate")
    query = engine.cohort.rows[engine.cohort.part_neurons["validation"][0]]
    proof_path = Path(engine.world.out) / "states" / "all_training.proof.json"
    proof = load(proof_path)
    proof["train_animals"].append(query.animal)             # the query's animal is put into the state
    proof_path.write_text(json.dumps(proof))
    fresh = engine_of(directory, "validate")
    with pytest.raises(E.LeakError, match="validation or test animals"):
        fresh.states.get("all_training")
    with pytest.raises(E.LeakError):
        fresh.predict()
    assert not (Path(fresh.world.out) / "validate" / "PREDICTION_COMMITMENT.json").exists()


def test_planted_leak_a_training_query_in_its_own_state_or_donor_set_is_refused(built, tmp_path):
    directory = clone(built, tmp_path / "w")
    engine = engine_of(directory, "develop")
    neuron = engine.library["ids"][0]
    row = engine.cohort.rows[neuron]
    assert len(engine.forecast(neuron)[0]) == len(E.METHODS)          # clean: forecasts are produced
    # (a) the query's own animal among the donors
    leaky = engine_of(directory, "develop")
    whole = leaky.scope(None)
    leaky.scope = lambda exclude_animal: whole
    with pytest.raises(E.LeakError, match="is in the donor set"):
        leaky.forecast(neuron)
    # (b) the query's own animal in its out-of-fold state, planted in memory: the per-forecast check refuses
    state_id = engine.states.state_id_for(neuron)
    in_memory = engine_of(directory, "develop")
    animals, id_animals = in_memory.states.animal_sets(state_id)
    in_memory.states._sets[state_id] = (animals | {row.animal}, id_animals)
    with pytest.raises(E.LeakError, match="is in regional state"):
        in_memory.forecast(neuron)
    # (c) the same leak planted on disk: the loader refuses a proof that is not the planned fitting set
    proof_path = Path(engine.world.out) / "states" / f"{state_id}.proof.json"
    proof = load(proof_path)
    proof["train_animals"] = sorted(proof["train_animals"] + [row.animal])
    proof_path.write_text(json.dumps(proof))
    with pytest.raises(E.LeakError, match="planned fitting animals"):
        engine_of(directory, "develop").forecast(neuron)


def test_a_state_whose_bytes_proof_or_manifest_disagree_is_refused(built, tmp_path):
    directory = clone(built, tmp_path / "w")
    states = directory / "engine/states"
    assert engine_of(directory, "validate").states.get("all_training")
    model = load(states / "all_training.json")
    model["blocks"][0]["weights"] = list(reversed(model["blocks"][0]["weights"]))          # another model, same proof
    (states / "all_training.json").write_text(E.canonical(model))
    with pytest.raises(E.EngineError, match="do not carry one hash"):
        engine_of(directory, "validate").states.get("all_training")
    second = clone(built, tmp_path / "w2")
    proof_path = second / "engine/states/all_training.proof.json"
    proof = load(proof_path)
    proof["test_ids"] = proof["test_ids"][:-1]                                             # serves another set than planned
    proof_path.write_text(json.dumps(proof))
    with pytest.raises(E.EngineError, match="planned served neurons"):
        engine_of(second, "validate").states.get("all_training")


def test_a_held_out_neuron_can_never_be_a_donor(built):
    engine = engine_of(built, "validate")
    for part in ("validation", "test"):
        with pytest.raises(E.LeakError):
            engine.donor(engine.cohort.part_neurons[part][0])
    assert all(engine.cohort.rows[n].part == "train" for n in engine.library["ids"])


def test_the_manifest_audit_catches_a_planted_donor_from_the_query_animal(pipeline):
    directory, _ = pipeline
    engine = engine_of(directory, "develop")
    path = Path(engine.world.out) / "develop" / "prediction_manifest.json"
    manifest = load(path)
    assert engine.audit_manifest(manifest)
    row = manifest["predictions"][0]
    own = next(n for n in manifest["training_ids"] if n.split("_")[0] == row["animal"])
    pool = sorted(manifest["donor_lists"][row["donor_pool_sha256"]] + [own])
    for p in manifest["predictions"]:
        if p["donor_pool_sha256"] == row["donor_pool_sha256"]:
            p["donor_pool_sha256"] = E.ids_sha(pool)
    manifest["donor_lists"][E.ids_sha(pool)] = pool
    with pytest.raises(E.LeakError, match="the query's animal supplied a donor"):
        engine.audit_manifest(manifest)
    manifest = load(path)
    manifest["predictions"].pop()
    with pytest.raises(E.EngineError, match="not complete"):
        engine.audit_manifest(manifest)
    manifest = load(path)
    victim = next(p for p in manifest["predictions"] if p["model"] == "E4_local_kernel_regional")
    victim["regional_state_id"] = None                                    # a fitted-state forecast without its provenance
    with pytest.raises(E.EngineError, match="fitted-state provenance"):
        engine.audit_manifest(manifest)


# ------------------------------------------------------------------------------------------
# hidden-outcome invariance
# ------------------------------------------------------------------------------------------


def test_hidden_outcome_canary_every_forecast_is_bit_identical(pipeline):
    """Plan section 7: held-out traces and regional rows replaced by other values, everything rebuilt
    from scratch, every forecast compared with the committed one."""
    _, summary = pipeline
    for key, part in (("validation_canary", "validate"), ("test_canary", "test")):
        canary = summary[key]
        assert canary["bit_identical"] and canary["forecasts_differing"] == 0
        assert canary["forecasts_compared"] == summary[part]["declared"] * len(E.METHODS)
        assert canary["held_out_neurons_scrambled"] > 0 and canary["official_copies_skipped_because_the_replica_differed"] == 0


def test_forecasts_need_neither_held_out_traces_nor_held_out_regional_rows(built, tmp_path):
    first, second = clone(built, tmp_path / "a"), clone(built, tmp_path / "b")
    engine_of(first, "validate").predict()
    reference = load(first / "engine/validate/prediction_manifest.json")
    removed = 0
    for part in ("validation", "test"):
        for trace in (second / "data/traces" / part).rglob("*.swc"):
            trace.unlink()                                   # the held-out traces do not exist at all
            removed += 1
    garble_sealed_rows(second)
    engine = engine_of(second, "validate")
    state_before = load(first / "engine/state_manifest.json")["states"]["all_training"]["model_state_sha256"]
    engine.states.build_many(["all_training"])               # the table's bytes changed: the state is refitted
    assert engine.states.manifest()["states"]["all_training"]["model_state_sha256"] == state_before
    engine.predict()
    again = load(second / "engine/validate/prediction_manifest.json")
    assert removed > 0 and arrays(again) == arrays(reference) and len(arrays(again)) == len(reference["queries"]) * len(E.METHODS)
    with pytest.raises(E.StageBlocked, match="not on disk"):   # admission does need them, says so, and opens nothing
        engine_of(second, "validate").admit()
    assert not (second / "engine/validate/TRUTH_OPENED.json").exists()


def test_development_forecast_ignores_the_query_animals_own_traces_and_rows(built, tmp_path):
    directory = clone(built, tmp_path / "w")
    engine = engine_of(directory, "develop")
    engine.predict()
    before = load(directory / "engine/develop/prediction_manifest.json")
    animal = engine.cohort.rows[engine.library["ids"][0]].animal
    E.scramble_outcomes(directory, animals={animal})
    rebuilt = engine_of(directory, "develop")
    with pytest.raises(E.EngineError, match="run `ingest`"):
        rebuilt.library                                       # a library that is not the one of the traces on disk is never used
    rebuilt.build_library()
    with pytest.raises(E.EngineError, match="run `cache`"):
        rebuilt.cache                                         # nor a pair cache of another library
    status = rebuilt.build_cache()
    assert status["complete"] and 0 < status["computed_now"] < status["pairs"]          # unchanged pairs were imported
    rebuilt.states.build_many(list(rebuilt.states.plan))
    rebuilt.predict()
    after = load(directory / "engine/develop/prediction_manifest.json")
    mine = [k for k in arrays(before) if k[0].split("_")[0] == animal]
    assert mine and all(k in arrays(after) for k in mine)
    assert all(arrays(after)[k] == arrays(before)[k] for k in mine)          # its own forecasts: bit-identical
    others = [k for k in arrays(before) if k[0].split("_")[0] != animal and k in arrays(after)]
    assert any(arrays(after)[k] != arrays(before)[k] for k in others)        # it is a donor for the others


def test_forecasts_are_deterministic(built):
    neuron = E.Cohort(E.synthetic_world(built)).part_neurons["test"][0]
    one = {f["model"]: E.arbor_sha(f["arbor"]) for f in engine_of(built, "test").forecast(neuron)[0]}
    two = {f["model"]: E.arbor_sha(f["arbor"]) for f in engine_of(built, "test").forecast(neuron)[0]}
    assert one == two and set(one) == set(E.METHODS)


# ------------------------------------------------------------------------------------------
# attachment and consensus
# ------------------------------------------------------------------------------------------


def long_arbor(soma=(6000.0, 5000.0, 5000.0)):
    points = np.array([[0.0, 0, 0], [900.0, 0, 0], [1800.0, 300.0, 0], [900.0, 700.0, 0]]) + np.array(soma)
    return G.resample(G.Arbor(points, np.array([-1, 0, 1, 1]), np.array([1, 2, 2, 2])), 50.0)


@pytest.mark.parametrize("shift,expected,radius", [(100.0, "local_attachment_500", 500.0), (250.0, "local_attachment_1000", 1000.0),
                                                    (400.0, "translation", None)])
def test_distal_vertices_are_not_moved_by_the_attachment(built, shift, expected, radius):
    engine = engine_of(built, "develop")
    arbor = long_arbor()
    soma = arbor.xyz[0] + np.array([0.0, shift, 0.0])
    out, how = engine.attach(arbor, soma)
    assert how["executed"] == expected and np.array_equal(out.xyz[0], soma)
    if radius is None:
        assert np.allclose(out.xyz - arbor.xyz, soma - arbor.xyz[0]) and [t["radius_um"] for t in how["unsupported"]] == [500.0, 1000.0]
        return
    produced = {tuple(v) for v in out.xyz.tolist()}
    depth = arbor.path_lengths()
    distal = arbor.xyz[depth >= radius]
    assert len(distal) > 5 and all(tuple(v) in produced for v in distal.tolist())       # exactly the measured coordinates
    assert out.branches() == arbor.branches()
    assert np.linalg.norm(out.xyz[1:] - out.xyz[out.parent[1:]], axis=1).max() <= 50.0 + 1e-6


def test_attached_forecasts_keep_the_measured_distal_coordinates_of_their_donor(pipeline):
    directory, _ = pipeline
    engine = engine_of(directory, "test")
    manifest = load(Path(engine.world.out) / "test" / "prediction_manifest.json")
    checked = 0
    for p in manifest["predictions"]:
        radius = p["details"].get("attachment_radius_um") if p["family"] == "engine" else (
            p["details"].get("radius_um") if p["family"] == "frozen" and not p["fallback"] and p["executed_policy"] != "baseline" else None)
        if not radius or p["model"] == "E6_local_consensus":
            continue
        donor = engine.donor(p["selected_donor"])["arbor"]
        forecast = G.load_arbor(Path(engine.world.out) / "test" / p["path"])
        produced = {tuple(v) for v in forecast.xyz.tolist()}
        assert all(tuple(v) in produced for v in donor.xyz[donor.path_lengths() >= radius].tolist())
        checked += 1
    assert checked > 10


def test_consensus_pruning_keeps_shared_cable_only_and_never_moves_a_vertex():
    template = long_arbor()                                             # trunk, then two branches
    trunk_and_first = G.resample(G.Arbor(template.xyz[[0]].repeat(3, axis=0) + np.array([[0.0, 0, 0], [900.0, 0, 0], [1800.0, 300.0, 0]]),
                                         np.array([-1, 0, 1]), np.array([1, 2, 2])), 50.0)
    pruned, info = E.consensus_prune(template, [trunk_and_first, trunk_and_first, trunk_and_first], np.ones(3))
    kept = {tuple(v) for v in pruned.xyz.tolist()}
    assert kept <= {tuple(v) for v in template.xyz.tolist()} and len(pruned.xyz) < len(template.xyz)
    assert 0.5 < info["consensus_kept_cable_fraction"] < 1.0 and "fallback_reason" not in info
    assert pruned.xyz[:, 1].max() < 5000.0 + 300.0 + 200.0 + 1e-6       # the unshared branch (up to +700 in DV) is gone
    pruned.validate()
    # nobody shares anything: too little would remain, the template is returned whole
    far = template.shifted(template.xyz[0] + np.array([0.0, 0.0, 4000.0]))
    whole, info = E.consensus_prune(template, [far], np.ones(1))
    assert whole is template and info["fallback_reason"] == "consensus_cable_below_minimum"


# ------------------------------------------------------------------------------------------
# the all-pairs cache, and the frozen prototype on the whole training pool
# ------------------------------------------------------------------------------------------


def test_the_pair_cache_holds_the_frozen_pair_costs_and_the_referees_metric_bit_for_bit(built):
    engine = engine_of(built, "develop")
    cache, lib = engine.cache, engine.library
    assert cache.n == len(lib["ids"]) and cache.ids == lib["ids"] and cache.fingerprint == lib["fingerprint"]
    assert cache.pairs == cache.n * (cache.n - 1) // 2
    pool = engine.pool_donors(engine.scope(None))                       # every library neuron, sorted by ID
    D0, F0 = T.pair_costs(pool)                                         # the unmodified frozen function, every pair
    D, F = cache.sub(lib["ids"])
    assert np.array_equal(D0, D) and np.array_equal(F0, F) and not D.flags.writeable
    proof = cache.prove(count=cache.pairs)                              # every pair again, also by the referee's loader and metric
    assert proof["pairs_checked"] == cache.pairs and proof["frozen_pair_costs_bit_identical"] and proof["referee_metric_bit_identical"]
    index = load(cache.dir / "index.json")
    assert index["ids"] == lib["ids"] and index["metric"]["analysis/d110_geometry.py"] == E.PINNED_SOURCES["analysis/d110_geometry.py"]
    assert index["arbor_sha256"] == [lib["entries"][n]["arbor_sha256"] for n in lib["ids"]]
    sub_ids = lib["ids"][3:20:2]
    Ds, Fs = cache.sub(sub_ids)
    assert np.array_equal(Ds, T.pair_costs([engine.donor(n) for n in sub_ids])[0])


def test_the_pair_cache_build_stops_at_its_budget_resumes_and_gives_the_same_pinned_cache(built, tmp_path, monkeypatch):
    engine = engine_of(built, "develop")
    reference = engine.cache.sha256
    monkeypatch.setitem(E.CONFIG["pair_cache"], "block_pairs", 60)
    lib = engine.library
    hashes = {n: lib["entries"][n]["arbor_sha256"] for n in lib["ids"]}
    cache = E.PairCache(tmp_path / "paircache", Path(engine.world.out) / "library", lib["ids"], hashes)
    first = cache.build(minutes=0.0)                                     # the budget is used up after the first blocks
    assert not first["complete"] and 0 < first["pairs_known"] < first["pairs"]
    with pytest.raises(E.EngineError, match="not complete"):
        cache.load()
    second = cache.build(minutes=0.0)
    assert second["pairs_known"] > first["pairs_known"]
    last = cache.build()
    assert last["complete"] and last["computed_now"] == last["pairs"] - second["pairs_known"]
    assert cache.load().sha256 == reference                              # same matrices, same pinned identity
    assert cache.build()["computed_now"] == 0


def test_a_changed_pair_cache_or_another_library_is_refused(built, tmp_path):
    directory = clone(built, tmp_path / "w")
    engine = engine_of(directory, "develop")
    path = engine.cache.dir / "D.npy"
    data = bytearray(path.read_bytes())
    data[-8:] = np.float64(123.456).tobytes()                            # one stored distance changed after completion
    path.write_bytes(bytes(data))
    with pytest.raises(E.EngineError, match="changed after the cache was completed"):
        engine_of(directory, "develop").cache
    with pytest.raises(E.EngineError, match="changed after the cache was completed"):
        engine_of(directory, "validate").predict()
    assert not (directory / "engine/validate/PREDICTION_COMMITMENT.json").exists()
    with pytest.raises(E.EngineError, match="distinct and sorted"):
        E.PairCache(tmp_path / "x", tmp_path, ["b", "a"], {"a": "1", "b": "2"})


def test_the_unmodified_frozen_policies_on_the_whole_pool_are_bit_identical_with_and_without_the_cache(built):
    engine = engine_of(built, "develop")
    original_pair, original_medoid = T.pair_costs, G.training_medoid
    proof = E.prove_memo(engine, queries=3)
    assert proof["queries"] == 3 and proof["all_bit_identical"]
    assert all(all(r["policies"].values()) and len(r["policies"]) == 5 for r in proof["rows"])
    assert all(r["donors"] >= len(engine.library["ids"]) - 6 for r in proof["rows"])       # the whole library but the query's animal
    assert T.pair_costs is original_pair and G.training_medoid is original_medoid          # the patch is undone
    # donors the memo does not know are handed to the frozen functions themselves
    strangers = [{"neuron": f"X{i}_001", "animal": f"X{i}", "arbor": long_arbor((6000.0 + 300 * i, 5000.0, 5000.0))} for i in range(3)]
    with engine.memo.frozen_memo():
        D, F = T.pair_costs(strangers)
        index = G.training_medoid(strangers)
    D0, F0 = T.pair_costs(strangers)
    assert np.array_equal(D, D0) and np.array_equal(F, F0) and index == G.training_medoid(strangers)


def test_forecasting_stops_at_its_budget_and_resumes_to_the_same_committed_forecasts(built, tmp_path, monkeypatch):
    """The frozen prototype on the whole pool is slow at full size, so `predict` is resumable until its commitment."""
    reference, directory = clone(built, tmp_path / "ref"), clone(built, tmp_path / "w")
    engine_of(reference, "validate").predict()
    first = engine_of(directory, "validate").predict(minutes=0.0)         # the budget is used up after the first query
    base = directory / "engine/validate"
    assert first["status"] == "incomplete" and 0 < first["queries_forecast"] < first["queries"]
    assert not (base / "PREDICTION_COMMITMENT.json").exists() and not (base / "prediction_manifest.json").exists()
    with pytest.raises(E.EngineError, match="no committed forecast bundle"):
        engine_of(directory, "validate").admit()                          # nothing can be opened against a partial run
    second = engine_of(directory, "validate").predict(minutes=0.0)
    assert second["queries_forecast"] == first["queries_forecast"] + 1    # earlier queries were kept, not redone
    done = engine_of(directory, "validate").predict()
    assert done["status"] == "committed" and not (base / "partial").exists()
    assert arrays(load(base / "prediction_manifest.json")) == arrays(load(reference / "engine/validate/prediction_manifest.json"))
    # a partial run made by another method is discarded, never completed
    third = clone(built, tmp_path / "w3")
    engine_of(third, "validate").predict(minutes=0.0)
    kept = sorted(p.name for p in (third / "engine/validate/partial").glob("q_*.json"))
    monkeypatch.setattr(E, "CONFIG_SHA", "d" * 64)
    engine_of(third, "validate").predict(minutes=0.0)
    assert sorted(p.name for p in (third / "engine/validate/partial").glob("q_*.json")) == kept      # started again from the first query
    assert load(third / "engine/validate/partial/PINS.json")["pins"]["config_sha256"] == "d" * 64


def test_every_method_draws_on_the_whole_training_library(pipeline):
    directory, _ = pipeline
    engine = engine_of(directory, "test")
    manifest = load(Path(engine.world.out) / "test" / "prediction_manifest.json")
    assert list(manifest["donor_lists"].values()) == [manifest["training_ids"]] == [engine.library["ids"]]
    assert manifest["pair_cache"]["sha256"] == engine.cache.sha256 and manifest["donor_pool"]["sampled panels"] == "none"
    proof = manifest["pair_cache"]["proof"]                               # sampled pairs re-proved when the forecasts were made
    assert proof["frozen_pair_costs_bit_identical"] and proof["referee_metric_bit_identical"] and proof["pairs_checked"] == 24
    for p in manifest["predictions"]:
        assert p["details"].get("pool_donors", len(manifest["training_ids"])) == len(manifest["training_ids"])
        if p["family"] == "engine":
            assert len(p["donor_ids"]) == E.CONFIG["engine"]["k_neighbours"] and p["selected_donor"] in p["donor_ids"]
        else:
            assert p["donor_ids"] == [p["selected_donor"]]
    develop = load(Path(engine.world.out) / "develop" / "prediction_manifest.json")
    for p in develop["predictions"][:200]:
        pool = develop["donor_lists"][p["donor_pool_sha256"]]
        assert pool == [n for n in develop["training_ids"] if n.split("_")[0] != p["animal"]]      # all other animals


# ------------------------------------------------------------------------------------------
# the four official copies: the referee's own copying_baselines
# ------------------------------------------------------------------------------------------


def test_the_engines_replica_of_the_four_copies_equals_the_referees_copying_baselines(built):
    engine = engine_of(built, "develop")
    result = E.check_official_copies(engine, max_pairs=10 ** 9)
    assert result["pools_checked"] == result["pools"] == len(set(engine.library["animal"].tolist()))
    assert result["all_equal"] and result["differences"] == 0 and result["forecasts_compared"] == 4 * len(engine.library["ids"])
    # held-out queries: one pool, the whole library
    held = engine_of(built, "test")
    held.prepare_prediction()
    scope, neurons = held.scope(None), held.cohort.part_neurons["test"]
    theirs = held.referee_copies(scope, neurons)
    for n in neurons:
        row = held.cohort.rows[n]
        mine = held.replica_copies(row, scope)
        for model in E.OFFICIAL_COPIES:
            assert scope["ids"][mine[model]] == theirs["donor_proof"][model][n]
            arbor = held.donor(scope["ids"][mine[model]])["arbor"].shifted(np.array(row.soma))
            assert E.same_tree(arbor, theirs["predictions"][model][n]) and np.array_equal(arbor.xyz[0], np.array(row.soma))
        assert theirs["same_peptide_fallback"][n] == mine["fallback"]


def test_at_a_floating_point_tie_the_official_medoid_follows_the_referee_and_the_frozen_baseline_keeps_its_own_sum():
    """The review's counterexample (data/interim/axon_engine_codex_review/interop_probe.json): six straight donors whose
    third and fourth are exactly tied as medoid. The frozen function sums left to right, the referee uses a matrix product,
    and rounding makes them pick different donors. The engine keeps both: its memo reproduces the frozen choice (C1), its
    replica of the official copy reproduces the referee's choice (`medoid`). Neither is changed to match the other."""
    from types import SimpleNamespace
    lengths = [31.001035906967545, 33.887605838294576, 36.49199583972045, 38.601467154664064, 38.819810240162646, 39.93453987841596]
    record = ROOT / "data/interim/axon_engine_codex_review/interop_probe.json"
    if record.exists():
        assert load(record)["medoid_rounding_counterexample"]["lengths_um"] == lengths
    training = [{"neuron": f"{i:02d}_1", "animal": f"{i:02d}",
                 "arbor": G.Arbor(np.array([[0.0, 0, 0], [value, 0, 0]]), np.array([-1, 0]), np.array([1, 2]))} for i, value in enumerate(lengths)]
    ids, animals = [t["neuron"] for t in training], [t["animal"] for t in training]
    D, F = T.pair_costs(training)
    assert np.array_equal(D, np.abs(np.subtract.outer(lengths, lengths)) / 2)        # two straight lines: half the length difference
    frozen = G.training_medoid(training)                                             # the frozen function itself
    referee = E.referee().G.geometry_medoid(D, list(range(6)), animals)              # the referee's function itself
    memo = E.FrozenMemo(lambda: None)
    for t in training:
        memo._arbor[t["neuron"]] = t["arbor"]
    memo._pool = (E.ids_sha(ids), D, F)
    engine = object.__new__(E.Engine)
    engine.memo = SimpleNamespace(matrices=lambda _ids: (D, F))
    scope = {"ids": ids, "label": np.array(["a"] * 6), "animal": animals, "medoids": {}}
    assert memo._training_medoid(training) == frozen and E.Engine.replica_medoid(engine, scope, None) == referee
    assert {frozen, referee} <= {2, 3}                                               # the two tied donors; on this machine 2 and 3


@real_only
def test_real_local_traces_the_replica_equals_the_referees_copying_baselines(tmp_path):
    """Leave one animal out on the 37 local training traces (the first pools within a pair budget; `selfcheck` runs all 32)."""
    engine = E.Engine(E.real_world(out=tmp_path, legacy_only=True), "develop", workers=1, log=None)
    engine.build_library()
    engine.build_cache()
    proof = engine.cache.prove(count=40)
    assert proof["frozen_pair_costs_bit_identical"] and proof["referee_metric_bit_identical"]
    result = E.check_official_copies(engine, max_pairs=4000)
    assert result["pools_checked"] >= 5 and result["all_equal"] and result["differences"] == 0


def test_the_referees_arrays_are_the_official_copies_even_where_the_replica_would_differ(built, tmp_path, monkeypatch):
    directory = clone(built, tmp_path / "w")
    real = E.Engine.replica_copies

    def wrong(self, row, scope):                               # a replica that picks another medoid
        out = real(self, row, scope)
        out["medoid"] = (out["medoid"] + 1) % len(scope["ids"])
        return out

    monkeypatch.setattr(E.Engine, "replica_copies", wrong)
    engine = engine_of(directory, "validate")
    engine.predict()
    monkeypatch.undo()
    manifest = load(directory / "engine/validate/prediction_manifest.json")
    official = manifest["official_copies"]
    assert official["source"] == "referee copying_baselines" and official["api_verified"]
    assert official["replica_differences"] == len(manifest["queries"]) and official["forecasts_compared"] == 4 * len(manifest["queries"])
    check = engine_of(directory, "validate")
    check.prepare_prediction()
    theirs = check.referee_copies(check.scope(None), manifest["queries"])
    for p in manifest["predictions"]:
        if p["model"] in E.OFFICIAL_COPIES:
            assert p["official_source"] == "referee copying_baselines" and p["replica_agreed"] == (p["model"] != "medoid")
            assert p["donor_ids"] == [theirs["donor_proof"][p["model"]][p["query_id"]]]
            saved = G.load_arbor(directory / "engine/validate" / p["path"])
            assert E.same_tree(saved, theirs["predictions"][p["model"]][p["query_id"]]) and E.arbor_sha(saved) == p["arrays_sha256"]


def test_a_pool_too_large_for_the_inline_call_needs_copies_check_before_any_truth(built, tmp_path, monkeypatch):
    directory = clone(built, tmp_path / "w")
    monkeypatch.setitem(E.CONFIG["official_copies"], "heldout_inline_api_pair_limit", 10)      # pretend the pool is large
    monkeypatch.setattr(E, "CONFIG_SHA", E.digest(E.CONFIG))
    commit = engine_of(directory, "validate").predict()
    assert commit["official_copies_api_verified"] is False
    with pytest.raises(E.EngineError, match="copies-check"):
        engine_of(directory, "validate").admit()
    assert not (directory / "engine/validate/TRUTH_OPENED.json").exists()
    record = engine_of(directory, "validate").copies_check()
    assert record["all_equal"] and record["differences"] == 0 and record["forecasts_compared"] == 4 * record["queries"]
    assert engine_of(directory, "validate").admit()["admitted"] > 0


# ------------------------------------------------------------------------------------------
# controls, manifest and output contract
# ------------------------------------------------------------------------------------------


def test_prediction_manifest_is_complete_and_arbors_are_standard(pipeline):
    directory, summary = pipeline
    engine = engine_of(directory, "test")
    base = Path(engine.world.out) / "test"
    manifest = load(base / "prediction_manifest.json")
    commitment = load(base / "PREDICTION_COMMITMENT.json")
    assert commitment["manifest_sha256"] == E.sha_file(base / "prediction_manifest.json") == summary["test_commitment"]
    queries = engine.cohort.part_neurons["test"]
    assert manifest["queries"] == queries                                 # every declared query, QC failures included
    keys = {(p["query_id"], p["model"]) for p in manifest["predictions"]}
    assert keys == {(q, m) for q in queries for m in E.METHODS} and len(manifest["predictions"]) == len(keys)
    assert manifest["training_ids_sha256"] == E.ids_sha(manifest["training_ids"])
    assert not set(manifest["upstream_fit_animals"]) & engine.cohort.heldout_animals
    needed = {"query_id", "animal", "full_id", "part", "known_soma_um", "source_label", "soma_label", "model", "role",
              "requested_policy", "executed_policy", "fallback", "fallback_reason", "selected_donor", "donor_ids",
              "donor_pool_sha256", "regional_state_id", "regional_state_sha256", "upstream_fit_animals_sha256",
              "path", "sha256", "arrays_sha256", "nodes"}
    states = {s["state_id"]: s["model_state_sha256"] for s in manifest["regional_states"]}
    referee = E.referee()
    for p in manifest["predictions"]:
        assert needed <= set(p) and p["part"] == "test" and p["path"] == f"predictions/{p['model']}/{p['query_id']}.npz"
        row = engine.cohort.rows[p["query_id"]]
        assert p["known_soma_um"] == list(row.soma)                       # the table soma, never a trace soma
        with np.load(base / p["path"], allow_pickle=False) as data:
            assert set(data.files) == {"xyz_um", "parent", "types"}
        arbor = G.load_arbor(base / p["path"])
        assert np.array_equal(arbor.xyz[0], np.array(row.soma)) and arbor.parent[0] == -1 and len(arbor.xyz) == p["nodes"]
        assert E.sha_file(base / p["path"]) == p["sha256"] and E.arbor_sha(arbor) == p["arrays_sha256"]
        referee.tree_check(referee.G.load_tree(base / p["path"]), list(row.soma), 0.0)     # the referee's own tree checks, zero tolerance
        pool = manifest["donor_lists"][p["donor_pool_sha256"]]
        assert E.ids_sha(pool) == p["donor_pool_sha256"] and p["selected_donor"] in p["donor_ids"] and set(p["donor_ids"]) <= set(pool)
        assert all(engine.cohort.rows[d].part == "train" for d in pool)
        uses = p["model"] in E.USES_REGIONAL_STATE                        # the engine's own regional state (E7 names another one)
        assert (p["regional_state_id"] is not None) == uses
        if uses:
            assert states[p["regional_state_id"]] == p["regional_state_sha256"]
            assert manifest["upstream_lists"][p["upstream_fit_animals_sha256"]] == manifest["upstream_fit_animals"]
        if p["fallback"]:
            assert p["fallback_reason"]
    assert manifest["code_sha256"]["frozen"] == E.FROZEN_PINS["source_sha256"] and manifest["referee_sha256"] == E.referee_hashes()
    assert manifest["cohort_manifest_sha256"] == engine.cohort.sha["manifest"]
    assert manifest["pins"]["pair_cache_sha256"] == engine.cache.sha256
    assert manifest["pins"]["training_traces"]["listed"] == summary["library"]["training_traces_listed"]
    assert {log["query"] for log in manifest["prediction_log"]} == set(queries)
    assert all(not log["query_animal_in_state"] and not log["query_animal_in_donors"] for log in manifest["prediction_log"])


def test_controls_are_what_they_are_declared_to_be(pipeline):
    directory, _ = pipeline
    engine = engine_of(directory, "test")
    engine.prepare_prediction()
    lib = engine.library
    pool = engine.pool_donors(engine.scope(None))
    medoid = pool[G.training_medoid(pool)]                                # the frozen function itself, every pair recomputed
    referee = E.referee()
    D = T.pair_costs(pool)[0]
    for neuron in engine.cohort.part_neurons["test"][:4]:
        row = engine.cohort.rows[neuron]
        soma = np.array(row.soma)
        forecasts = {f["model"]: f for f in engine.forecast(neuron)[0]}
        d = np.linalg.norm(lib["soma"] - soma, axis=1)
        assert forecasts["nearest"]["selected_donor"] == lib["ids"][int(np.argmin(d))]
        assert E.same_arbor(forecasts["nearest"]["arbor"], engine.donor(lib["ids"][int(np.argmin(d))])["arbor"].shifted(soma))
        same = np.flatnonzero(lib["label"] == row.source_label)
        label_nearest, label_medoid = forecasts["peptide_nearest"], forecasts["peptide_medoid"]
        if len(same):
            assert label_nearest["selected_donor"] == lib["ids"][int(same[np.argmin(d[same])])] and not label_nearest["fallback"]
            pick = same[referee.G.geometry_medoid(D, same.tolist(), [lib["animal"][i] for i in same])]
            assert label_medoid["selected_donor"] == lib["ids"][pick] and not label_medoid["fallback"]
        else:
            assert label_nearest["fallback"] and label_nearest["selected_donor"] == forecasts["nearest"]["selected_donor"]
        whole = referee.G.geometry_medoid(D, list(range(len(pool))), lib["animal"].tolist())
        assert forecasts["medoid"]["selected_donor"] == lib["ids"][whole]
        c1 = forecasts["C1_medoid_copy"]
        assert c1["selected_donor"] == medoid["neuron"] and E.same_arbor(c1["arbor"], medoid["arbor"].shifted(soma))
        for name in ("C4_uniform_attach_500", "C5_uniform_attach_1000"):
            assert forecasts[name]["requested_policy"] == E.FROZEN_POLICY_OF[name]
            if forecasts[name]["fallback"]:                               # unsupported attachment: exact baseline
                assert E.same_arbor(forecasts[name]["arbor"], c1["arbor"])


def test_same_label_copies_fall_back_and_say_so_when_the_label_has_no_donor(built):
    engine = engine_of(built, "test")
    engine.prepare_prediction()
    neuron = engine.cohort.part_neurons["test"][0]
    row = engine.cohort.rows[neuron]
    for label in ("LabelWithoutDonor", ""):
        engine.cohort.rows[neuron] = E.NeuronRow(row.neuron, row.animal, row.full_id, label, row.tag, row.soma_label,
                                                 row.hemisphere, row.recon_type, row.part, row.klass, row.soma)
        forecasts = {f["model"]: f for f in engine.forecast(neuron)[0]}
        for model, plain in (("peptide_nearest", "nearest"), ("peptide_medoid", "medoid")):
            assert forecasts[model]["fallback"] and "no_training_donor" in forecasts[model]["fallback_reason"]
            assert forecasts[model]["selected_donor"] == forecasts[plain]["selected_donor"] and forecasts[model]["executed_policy"] == plain
        theirs = engine.referee_copies(engine.scope(None), [neuron])                 # the referee falls back the same way
        assert theirs["same_peptide_fallback"][neuron] is True
        assert theirs["donor_proof"]["peptide_medoid"][neuron] == forecasts["medoid"]["selected_donor"]


# ------------------------------------------------------------------------------------------
# gates and verdict
# ------------------------------------------------------------------------------------------


def summary_of(**override):
    base = {m: {"distance_um": 900.0, "f1_200": 0.20} for m in E.METHODS}
    for name, (distance, f1) in override.items():
        base[name] = {"distance_um": distance, "f1_200": f1}
    return base


def test_engine_gate_needs_both_measures_against_every_control():
    gate = E.engine_gate(summary_of(E2_local_kernel=(700.0, 0.30)), "baseline")
    assert gate["selected"] == "E2_local_kernel" and gate["verdict"] == "candidate_selected" and gate["descriptive_reference"] is None
    # better distance, F1 only equal to one control: not eligible
    gate = E.engine_gate(summary_of(E2_local_kernel=(700.0, 0.30), nearest=(1200.0, 0.30)), "baseline")
    assert gate["selected"] is None and gate["verdict"] == "no method beyond copying"
    assert not gate["candidates"]["E2_local_kernel"]["beats_control"]["nearest"]
    assert gate["strongest_control"] == {"distance_um": "medoid", "f1_200": "nearest"}
    assert gate["descriptive_reference"] == "E2_local_kernel"             # fills the referee's candidate slot; no claim
    # each of the seven controls counts: the four official copies, the frozen baseline, both uniform attachments
    for control in E.CONTROLS:
        gate = E.engine_gate(summary_of(E2_local_kernel=(700.0, 0.30), **{control: (650.0, 0.10)}), "baseline")
        assert gate["selected"] is None and gate["candidates"]["E2_local_kernel"]["beats_control"][control] is False
    # the best method on validation is itself a copying control: said so
    gate = E.engine_gate(summary_of(peptide_medoid=(500.0, 0.50)), "baseline")
    assert gate["verdict"] == "no method beyond copying" and gate["best_method_is_a_copying_control"] == {"distance_um": True, "f1_200": True}


def test_conditional_candidates_must_also_beat_their_matched_controls():
    # the label candidate beats every copy but not the uniform local control: the uniform one is chosen
    gate = E.engine_gate(summary_of(E3_local_kernel_label=(600.0, 0.35), E1_local_uniform=(650.0, 0.36)), "baseline")
    assert not gate["candidates"]["E3_local_kernel_label"]["eligible"] and gate["selected"] == "E1_local_uniform"
    gate = E.engine_gate(summary_of(E3_local_kernel_label=(600.0, 0.37), E1_local_uniform=(650.0, 0.36), E2_local_kernel=(640.0, 0.33)), "baseline")
    assert gate["selected"] == "E3_local_kernel_label" and gate["eligible"] == ["E1_local_uniform", "E2_local_kernel", "E3_local_kernel_label"]
    # a frozen conditional policy is admitted only when the frozen gate itself chose it, and must beat its uniform twin
    best = summary_of(F1_conditional_attach_500=(500.0, 0.50))
    assert E.engine_gate(best, "baseline")["selected"] is None
    assert E.engine_gate(best, "conditional_anchor_500")["selected"] == "F1_conditional_attach_500"
    twin = summary_of(F1_conditional_attach_500=(500.0, 0.50), C4_uniform_attach_500=(450.0, 0.20))
    assert E.engine_gate(twin, "conditional_anchor_500")["selected"] is None
    # ties: lowest distance, then highest F1, then the fixed order
    tie = E.engine_gate(summary_of(E1_local_uniform=(700.0, 0.30), E2_local_kernel=(700.0, 0.31)), "baseline")
    assert tie["selected"] == "E2_local_kernel"


def test_the_declared_family_limits_what_the_gate_may_select():
    summary = summary_of(E4_local_kernel_regional=(500.0, 0.50), E2_local_kernel=(700.0, 0.30))
    assert E.engine_gate(summary, "baseline", "all")["selected"] == "E4_local_kernel_regional"
    gate = E.engine_gate(summary, "baseline", "no_fitted_state")
    assert gate["selected"] == "E2_local_kernel" and gate["family_members"] == list(E.FAMILIES["no_fitted_state"])
    assert not gate["candidates"]["E4_local_kernel_regional"]["admitted"] and "family" in gate["candidates"]["E4_local_kernel_regional"]["note"]
    assert not set(E.FAMILIES["no_fitted_state"]) & set(E.USES_FITTED_STATE) and set(E.FAMILIES["all"]) == set(E.CANDIDATES)
    with pytest.raises(E.EngineError):
        E.engine_gate(summary, "baseline", "something_else")


def scored(per_animal_values):
    rows = []
    for animal, values in per_animal_values.items():
        metrics = {m: {k: 0.0 for k in E.KEPT_METRICS} for m in E.METHODS}
        for m in E.METHODS:
            distance, f1 = values.get(m, (900.0, 0.20))
            metrics[m]["distance_um"], metrics[m]["f1_200"] = distance, f1
        rows.append({"neuron": animal + "_001", "animal": animal, "source_label": "L", "soma_label": "N", "metrics": metrics})
    return E.aggregate(rows, E.METHODS)


def test_verdict_states_success_failure_the_unmet_matched_control_and_the_paired_margins():
    win = E.verdict_for(scored({a: {"E2_local_kernel": (700.0, 0.30)} for a in "ABCDEF"}), "E2_local_kernel", context="validation")
    assert win["success"] and win["beats_every_copying_control"] and win["verdict"] == "beats copying on both measures"
    assert win["versus_control"]["medoid"]["animals_better_on_both"] == 6 and set(win["versus_control"]) == set(E.CONTROLS)
    assert "paired 95% bootstrap interval" in win["statement"] and E.CLAIM_LIMIT in win["statement"]
    assert "-200.0 um" in win["statement"] and "+0.1000" in win["statement"]
    lose = E.verdict_for(scored({a: {"E2_local_kernel": (700.0, 0.19)} for a in "ABCDEF"}), "E2_local_kernel", context="validation")
    assert not lose["success"] and lose["verdict"] == "does not beat copying" and "FAILS" in lose["statement"]
    half = E.verdict_for(scored({a: {"E4_local_kernel_regional": (700.0, 0.30), "E1_local_uniform": (690.0, 0.31)} for a in "ABCDEF"}),
                         "E4_local_kernel_regional", context="validation")
    assert not half["success"] and half["beats_every_copying_control"] and not half["beats_every_matched_control"]
    assert half["verdict"].startswith("rule not met") and "E1_local_uniform" in half["statement"]
    none = E.verdict_for(scored({a: {} for a in "ABC"}), None, context="validation")
    assert not none["success"] and none["verdict"] == "no method beyond copying"


def referee_result(candidate, roles, values):
    """A result shaped like the referee's score_once output, for the verdict reader."""
    Ref = E.referee()
    models = [roles["candidate"], *roles["copy_baselines"], *roles["matched_controls"]]
    records = [{"id": f"{a}_1", "animal": a, "peptide": "P", "nucleus": "N",
                "models": {m: {"distance_um": values.get(m, (900.0, 0.20))[0], "f1_200": values.get(m, (900.0, 0.20))[1]} for m in models}}
               for a in "ABCDEF"]
    return {"overall": Ref.summarize(records, models, roles), "queries": records}


def test_the_test_verdict_is_read_from_the_referees_own_result():
    roles = E.Engine.referee_roles({"selected": "E3_local_kernel_label", "descriptive_reference": None})
    assert roles == {"candidate": "E3_local_kernel_label", "copy_baselines": list(E.COPY_BASELINES),
                     "matched_controls": ["C4_uniform_attach_500", "C5_uniform_attach_1000", "E1_local_uniform", "E2_local_kernel"]}
    assert E.Engine.referee_roles({"selected": "F1_conditional_attach_500", "descriptive_reference": None})["matched_controls"] == list(E.ATTACHMENT_CONTROLS)
    assert E.Engine.referee_roles({"selected": "E2_local_kernel", "descriptive_reference": None})["matched_controls"] == list(E.ATTACHMENT_CONTROLS)
    win = E.verdict_from_referee(referee_result("E3", roles, {"E3_local_kernel_label": (700.0, 0.30)}), roles, "E3_local_kernel_label")
    assert win["success"] and win["referee_beats_all_copying"] and win["referee_beats_all_matched"] and win["scored_by"].startswith("independent referee")
    assert E.CLAIM_LIMIT in win["statement"] and win["animals"] == 6
    matched = E.verdict_from_referee(referee_result("E3", roles, {"E3_local_kernel_label": (700.0, 0.30), "E1_local_uniform": (650.0, 0.35)}),
                                     roles, "E3_local_kernel_label")
    assert not matched["success"] and matched["verdict"].startswith("rule not met") and "E1_local_uniform" in matched["statement"]
    attach = E.verdict_from_referee(referee_result("E3", roles, {"E3_local_kernel_label": (700.0, 0.30), "C5_uniform_attach_1000": (650.0, 0.35)}),
                                    roles, "E3_local_kernel_label")
    assert not attach["success"] and attach["verdict"] == "does not beat copying" and "C5_uniform_attach_1000" in attach["statement"]
    copy = E.verdict_from_referee(referee_result("E3", roles, {"E3_local_kernel_label": (700.0, 0.30), "peptide_medoid": (650.0, 0.35)}),
                                  roles, "E3_local_kernel_label")
    assert not copy["success"] and not copy["referee_beats_all_copying"] and "peptide_medoid" in copy["statement"]
    # nothing was selected on validation: no claim, whatever the test numbers of the reference are
    nothing = E.verdict_from_referee(referee_result("E3", roles, {"E3_local_kernel_label": (100.0, 0.90)}), roles, None)
    assert not nothing["success"] and nothing["verdict"] == "no method beyond copying" and "descriptive reference" in nothing["statement"]
    with pytest.raises(E.EngineError, match="not the frozen selection"):
        E.verdict_from_referee(referee_result("E3", roles, {}), roles, "E2_local_kernel")


def test_the_frozen_gate_runs_unmodified_on_the_whole_pool_and_its_decision_reproduces_predict_gated(frozen_world):
    engine = engine_of(frozen_world, "test")
    engine.prepare_prediction()
    gate = load(Path(engine.world.out) / "validate" / "frozen_gate.json")
    assert gate["schema"] == "brm_template_gate_v1" and gate["selected_policy"] in T.POLICIES
    assert len(gate["proof"]["training_ids"]) == len(engine.library["ids"])                # every training donor is bound
    assert engine.verify_frozen_gate(gate)
    method = E.METHOD_OF_FROZEN_POLICY[gate["selected_policy"]]
    pool, regional = engine.pool_donors(engine.scope(None)), engine.states.get("all_training")
    for neuron in engine.cohort.part_neurons["test"][:2]:
        row = engine.cohort.rows[neuron]
        soma = np.array(row.soma)
        expected, info = T.predict_gated(T.Query(row.neuron, row.animal, row.full_id, soma, soma, row.tag), pool, regional, gate)
        saved = G.load_arbor(Path(engine.world.out) / "test" / "predictions" / method / f"{neuron}.npz")
        assert E.same_arbor(expected, saved) and info["gate_sha256"] == gate["gate_sha256"]
    changed = dict(gate, selected_policy="baseline" if gate["selected_policy"] != "baseline" else "uniform_anchor_500")
    with pytest.raises(E.FreezeMismatch):
        engine.verify_frozen_gate(changed)


# ------------------------------------------------------------------------------------------
# the test goes through the independent referee
# ------------------------------------------------------------------------------------------


def test_the_test_was_declared_sealed_and_scored_by_the_referee(pipeline):
    directory, summary = pipeline
    engine = engine_of(directory, "test")
    stage = Path(engine.world.out) / "test"
    base = stage / "referee"
    Ref = E.referee()
    for name in ("declaration.json", "COHORT_FROZEN.json", "prediction_manifest.json", "PREDICTION_SEAL.json", "SCORE_STARTED.json",
                 "metrics.json", "SCORE_COMPLETE.json"):
        assert (base / name).exists(), name
    assert summary["test"]["referee_functions_called"] == ["freeze_declaration", "copying_baselines", "seal_predictions", "score_once"]
    assert summary["test"]["scored_by"].startswith("independent referee")
    declaration, frozen = Ref.frozen_declaration(base)                    # its own re-check of the frozen declaration
    seal, result, mine = load(base / "PREDICTION_SEAL.json"), load(base / "metrics.json"), load(stage / "metrics.json")
    assert seal["manifest_sha256"] == Ref.sha(base / "prediction_manifest.json") and result["status"] == "scored_once"
    assert seal["source_sha256"] == Ref.source_hashes() == E.referee_hashes() == mine["referee"]["sha256"]
    # declaration: cohort after admission, exclusions, known metadata only, training trees with hashes, roles, contract, artifacts
    admission = load(stage / "ADMISSION_COMPLETE.json")
    assert [q["id"] for q in declaration["splits"]["test"]] == admission["admitted_ids"]
    assert all(set(q) == {"id", "animal", "soma_um", "peptide", "nucleus"} for q in declaration["splits"]["test"] + declaration["splits"]["validation"])
    assert set(admission["excluded"]) <= set(declaration["excluded_ids"]) and len(admission["excluded"]) == 4
    cohort = engine.cohort
    listed = {r["id"] for split in declaration["splits"].values() for r in split} | set(declaration["excluded_ids"])
    assert listed == set(cohort.order)                                    # every cohort neuron is in a split or excluded
    assert [r["id"] for r in declaration["splits"]["train"]] == engine.library["ids"]
    for r in declaration["splits"]["train"][:5]:
        assert Ref.sha(base / r["arbor_path"]) == r["arbor_sha256"] == E.sha_file(Path(engine.world.out) / "library" / f"{r['id']}.npz")
        assert r["soma_um"] == list(cohort.rows[r["id"]].soma) and r["peptide"] == cohort.rows[r["id"]].source_label
    assert declaration["config"] == Ref.CONFIG and declaration["upstream_fit_animals"] == []
    assert declaration["geometry_contract"]["units"] == "um" and declaration["geometry_contract"]["axes"] == ["AP", "DV", "ML"]
    artifacts = declaration["method_artifacts"]
    for needed in ("method/source/analysis/axon_engine.py", "method/source/analysis/brm_template.py", "method/protocol/axon_engine_protocol.md",
                   "method/freeze/frozen.json", "method/regional_state/all_training.json", "method/pair_cache/COMPLETE.json",
                   "method/test/PREDICTION_COMMITMENT.json", "method/test/ADMISSION_COMPLETE.json"):
        assert artifacts[needed] == Ref.sha(base / needed)
    roles = declaration["models"]
    assert roles == summary["test"]["referee_roles"] and set(E.OFFICIAL_COPIES) <= set(roles["copy_baselines"])
    assert roles["matched_controls"][:2] == list(E.ATTACHMENT_CONTROLS)
    # bundle: one forecast per admitted query per declared method, byte for byte the committed files
    models = [roles["candidate"], *roles["copy_baselines"], *roles["matched_controls"]]
    bundle, committed = load(base / "prediction_manifest.json"), load(stage / "prediction_manifest.json")
    by_key = {(p["query_id"], p["model"]): p for p in committed["predictions"]}
    assert {(r["query_id"], r["model"]) for r in bundle["predictions"]} == {(q, m) for q in admission["admitted_ids"] for m in models}
    for r in bundle["predictions"]:
        assert set(r) == {"query_id", "model", "path", "sha256", "donor_ids", "upstream_fit_animals", "hemisphere"}
        assert r["hemisphere"] == cohort.rows[r["query_id"]].hemisphere in ("Left", "Right")      # permitted metadata, kept for the viewer
        assert r["sha256"] == by_key[(r["query_id"], r["model"])]["sha256"] == Ref.sha(base / r["path"])
        assert r["path"] == f"predictions/{r['model']}/{r['query_id']}.npz" and r["upstream_fit_animals"] == []
        assert len(r["donor_ids"]) == 1 or r["model"] not in roles["copy_baselines"]
    assert bundle["declaration_sha256"] == frozen["declaration_sha256"] == mine["referee"]["declaration_sha256"]
    # the four copies in the bundle are what the referee recreates by itself
    donors = [{**r, "arbor": Ref.G.load_tree(base / r["arbor_path"])} for r in declaration["splits"]["train"]]
    copies = Ref.copying_baselines(donors, declaration["splits"]["test"])
    for r in bundle["predictions"]:
        if r["model"] in E.OFFICIAL_COPIES:
            assert r["donor_ids"] == [copies["donor_proof"][r["model"]][r["query_id"]]]
            assert all(np.array_equal(a, b) for a, b in zip(Ref.G.load_tree(base / r["path"]), copies["predictions"][r["model"]][r["query_id"]]))
    # scoring: the engine's table is the referee's, bit for bit; truth is stored and hashed by the referee
    assert mine["scored_by"].startswith("independent referee") and mine["methods_scored_by_the_referee"] == models
    for m in models:
        for key in E.SCORE_KEYS:
            assert mine["summary"][m][key] == result["overall"]["summary"][m][key]
    assert mine["scored_queries"] == result["overall"]["neurons"] == len(admission["admitted_ids"])
    assert mine["scored_animals"] == result["overall"]["independent_animals"]
    for q in result["queries"]:                                           # the truth the referee stored is the admitted truth
        assert Ref.sha(base / q["truth_path"]) == q["truth_sha256"]
        assert E.arbor_sha(G.load_arbor(base / q["truth_path"])) == admission["truth"][q["id"]]["arrays_sha256"]
    verdict = mine["verdict"]
    assert verdict["referee_beats_all_copying"] == result["overall"]["beats_all_copying"]
    assert verdict["success"] == bool(mine["frozen_selection"]["selected"] and result["overall"]["beats_all_copying"] and result["overall"]["beats_all_matched"])
    assert set(mine["other_candidates_descriptive"]["summary"]) == set(E.CANDIDATES) - set(models)
    assert "no claim" in mine["other_candidates_descriptive"]["note"]


def test_the_four_referee_functions_are_really_called_in_the_contracts_order(frozen_world, tmp_path, monkeypatch):
    directory = clone(frozen_world, tmp_path / "w")
    freeze = load(directory / "engine/freeze/frozen.json")
    if not freeze["referee"]["can_seal_the_test_bundle"]:
        pytest.skip("the frozen selection of this world has a fitted state, which the referee on disk cannot seal (tested below)")
    Ref = E.referee()
    calls = []
    for name in ("freeze_declaration", "copying_baselines", "seal_predictions", "score_once"):
        original = getattr(Ref, name)
        monkeypatch.setattr(Ref, name, lambda *a, _o=original, _n=name, **k: (calls.append(_n), _o(*a, **k))[1])
    engine_of(directory, "test").admit()
    assert calls == []                                                    # admission involves no referee function and no forecast
    engine_of(directory, "test").seal()
    assert calls == ["freeze_declaration", "copying_baselines", "seal_predictions", "copying_baselines"]   # the last: inside its seal
    output = engine_of(directory, "test").score()
    assert calls[-1] == "score_once" and calls.count("score_once") == 1 and output["scored_by"].startswith("independent referee")
    with pytest.raises(E.OnceOnly):
        engine_of(directory, "test").score()
    assert calls.count("score_once") == 1


def test_admission_reads_no_forecast_and_fixes_the_exclusions_before_the_referee_declaration(frozen_world, tmp_path):
    directory = clone(frozen_world, tmp_path / "w")
    stage = directory / "engine/test"
    engine = engine_of(directory, "test")
    freeze = load(directory / "engine/freeze/frozen.json")
    if not freeze["referee"]["can_seal_the_test_bundle"]:
        pytest.skip("the frozen selection of this world has a fitted state, which the referee on disk cannot seal")
    summary = engine.admit()
    assert summary["forecasts_read"] == 0 and summary["declared"] - summary["admitted"] == len(summary["excluded"]) == 4
    assert summary["excluded_by_rule"] == {"2": 1, "4": 1, "5": 1, "6": 1}
    # the admission is a function of the trace files alone: the same rules applied here, one by one, give the same log
    log = load(stage / "truth_ingestion_log.json")["neurons"]
    store = E.TraceStore(engine.world, engine.cohort, "test")
    store.unlock_truth("check")
    for neuron in engine.cohort.part_neurons["test"]:
        row, copy = engine.cohort.rows[neuron], store.open(neuron, "truth")
        entry, arbor = E.ingest_trace(copy["path"], frame=copy["frame"], hemisphere=row.hemisphere, recon_type=row.recon_type, table_soma=row.soma)
        assert entry["status"] == log[neuron]["status"] and entry["rules_violated"] == log[neuron]["rules_violated"]
        if arbor is not None:
            assert E.arbor_sha(arbor) == summary["truth"][neuron]["arrays_sha256"] and (stage / "truth" / f"{neuron}.npz").exists()
        else:
            assert neuron in summary["excluded"] and not (stage / "truth" / f"{neuron}.npz").exists()
    marker = load(stage / "TRUTH_OPENED.json")
    assert marker["manifest_sha256"] == load(stage / "PREDICTION_COMMITMENT.json")["manifest_sha256"] and "no forecast is read" in marker["purpose"]
    sealed = engine_of(directory, "test").seal()
    declaration = load(stage / "referee/declaration.json")
    assert set(summary["excluded"]) <= set(declaration["excluded_ids"]) and sealed["queries"] == summary["admitted"]
    assert not any((stage / "referee/predictions").rglob(f"{next(iter(summary['excluded']))}.npz"))
    result = engine_of(directory, "test").score()                         # the single referee attempt is not wasted on a missing truth
    assert result["scored_queries"] == summary["admitted"] and (stage / "referee/SCORE_COMPLETE.json").exists()


def test_referee_preflight_and_the_block_before_any_test_truth(frozen_world, tmp_path):
    """The referee on disk is probed with its own functions. While it refuses an honest training-animal upstream list
    (reports/prediction_referee_interop_finding.md), a bundle with a fitted-state method is blocked BEFORE truth is opened."""
    preflight = E.referee_preflight()
    assert preflight["no_upstream_fit_sealable"] is True and preflight["referee_sha256"] == E.referee_hashes()
    engine = engine_of(clone(frozen_world, tmp_path / "w"), "test")
    stateless = E.Engine.referee_roles({"selected": "E2_local_kernel", "descriptive_reference": None})
    fitted = E.Engine.referee_roles({"selected": "E4_local_kernel_regional", "descriptive_reference": None})
    assert engine.require_referee_can_seal(stateless)
    if preflight["fitted_state_sealable"]:
        assert engine.require_referee_can_seal(fitted) and preflight["heldout_upstream_refused"]
        return
    assert "Upstream fit saw validation/test animals" in preflight["fitted_state_error"]
    with pytest.raises(E.StageBlocked, match="key order of 'splits'"):
        engine.require_referee_can_seal(fitted)
    freeze = load(Path(engine.world.out) / "freeze/frozen.json")
    assert freeze["referee"]["preflight"]["fitted_state_sealable"] is False
    if freeze["referee"]["bundle_methods_with_a_fitted_state"]:
        assert freeze["referee"]["can_seal_the_test_bundle"] is False
        for action in (engine.admit, engine.seal):
            with pytest.raises(E.StageBlocked, match="No test truth was opened"):
                action()
        assert not (Path(engine.world.out) / "test" / "TRUTH_OPENED.json").exists()


def test_a_fitted_state_bundle_declares_its_training_animals_at_all_three_levels_and_meets_the_referees_checks(frozen_world, tmp_path, monkeypatch):
    """A candidate with a fitted regional state, forced here whatever the gate chose. The bundle the engine builds lists
    the training animals behind the state in the declaration, in the manifest and in the candidate's rows, never an empty
    list, and passes the referee's own validators. Its seal accepts it with the referee as corrected on 2026-11-11; with the
    earlier referee (held-out set dependent on the key order of 'splits') the seal refused it after its JSON round trip."""
    preflight = E.referee_preflight()
    directory = clone(frozen_world, tmp_path / "w")
    forced = E.Engine.referee_roles({"selected": "E4_local_kernel_regional", "descriptive_reference": None})
    monkeypatch.setattr(E.Engine, "referee_roles", staticmethod(lambda freeze: forced))
    if not preflight["fitted_state_sealable"]:
        monkeypatch.setattr(E.Engine, "require_referee_can_seal", lambda self, roles: preflight)   # look past the block, in this test only
    engine_of(directory, "test").admit()
    engine = engine_of(directory, "test")
    base, declaration, bundle, info = engine.prepare_bundle(engine.verify_freeze(), forced)
    Ref = E.referee()
    state_animals = sorted(engine.states.get("all_training").proof["train_animals"])
    assert declaration["upstream_fit_animals"] == bundle["upstream_fit_animals"] == state_animals
    assert set(state_animals) == set(engine.cohort.part_animals["train"]) and not set(state_animals) & engine.cohort.heldout_animals
    for r in bundle["predictions"]:
        assert r["upstream_fit_animals"] == (state_animals if r["model"] in E.USES_FITTED_STATE else [])
    assert {r["model"] for r in bundle["predictions"] if r["upstream_fit_animals"]} == {"E4_local_kernel_regional"}
    trees = Ref.validate_bundle(base, bundle, declaration, fresh_baselines=True)       # the referee's validators, declaration as declared
    assert len(trees) == info["forecasts"]
    if preflight["fitted_state_sealable"]:
        seal = Ref.seal_predictions(base, bundle)
        assert seal["manifest_sha256"] == Ref.sha(base / "prediction_manifest.json") and (base / "PREDICTION_SEAL.json").exists()
        # a held-out animal in any upstream list is refused by the referee, and by the engine's own audit before that
        bad = clone(frozen_world, tmp_path / "bad")
        engine_of(bad, "test").admit()
        other = engine_of(bad, "test")
        base2, declaration2, bundle2, _ = other.prepare_bundle(other.verify_freeze(), forced)
        bundle2["predictions"][0]["upstream_fit_animals"] = [other.cohort.part_animals["test"][0]]
        with pytest.raises(ValueError, match="Upstream fit saw validation/test animals"):
            Ref.seal_predictions(base2, bundle2)
    else:
        with pytest.raises(ValueError, match="Upstream fit saw validation/test animals"):
            Ref.seal_predictions(base, bundle)                                         # its seal, after its JSON round trip
        assert not (base / "PREDICTION_SEAL.json").exists()
    assert not (base / "SCORE_STARTED.json").exists()


def test_with_the_whole_family_a_fitted_state_method_is_scored_by_the_referee_or_blocked_before_truth(pipeline_all):
    """The whole candidate family on a tiny world, end to end. Whatever fills the referee's candidate slot here has a
    fitted regional state: the bundle goes through the referee with honest upstream lists, or, if the referee on disk
    cannot seal it, the test is blocked before any test truth is opened."""
    directory, summary = pipeline_all
    out = Path(E.synthetic_world(directory).out)
    freeze = load(out / "freeze/frozen.json")
    fitted = freeze["referee"]["bundle_methods_with_a_fitted_state"]
    assert freeze["family"] == "all" and fitted, "this fixture is meant to put a fitted-state method into the bundle"
    if not freeze["referee"]["can_seal_the_test_bundle"]:
        assert "blocked" in summary["test"] and summary["test"]["test_truth_opened"] is False
        assert not (out / "test/TRUTH_OPENED.json").exists() and not (out / "test/referee").exists()
        return
    assert summary["test"]["referee_functions_called"] == ["freeze_declaration", "copying_baselines", "seal_predictions", "score_once"]
    declaration, bundle = load(out / "test/referee/declaration.json"), load(out / "test/referee/prediction_manifest.json")
    engine = engine_of(directory, "test")
    state_animals = sorted(engine.states.get("all_training").proof["train_animals"])
    assert declaration["upstream_fit_animals"] == bundle["upstream_fit_animals"] == state_animals == sorted(engine.cohort.part_animals["train"])
    for r in bundle["predictions"]:
        assert r["upstream_fit_animals"] == (state_animals if r["model"] in fitted else [])
    metrics = load(out / "test/metrics.json")
    assert metrics["scored_by"].startswith("independent referee") and metrics["referee"]["status"] == "scored_once"
    if freeze["selected"] is None:
        assert metrics["verdict"]["verdict"] == "no method beyond copying" and metrics["verdict"]["success"] is False
        assert metrics["verdict"]["referee_candidate_slot"] == freeze["descriptive_reference"]
    assert all(text.startswith("refused") for text in summary["once_only"].values())


# ------------------------------------------------------------------------------------------
# once only, commitment, binding of validation to the frozen method, freeze
# ------------------------------------------------------------------------------------------


def test_the_pipeline_ran_end_to_end_and_refuses_every_repeat(pipeline):
    directory, summary = pipeline
    assert len(summary["once_only"]) == 8 and all(text.startswith("refused") for text in summary["once_only"].values())
    out = Path(E.synthetic_world(directory).out)
    for stage in ("validate", "test"):
        for name in ("PREDICTION_COMMITMENT.json", "TRUTH_OPENED.json", "ADMISSION_COMPLETE.json", "SCORE_COMPLETE.json", "metrics.json", "tables.md"):
            assert (out / stage / name).exists()
        opened = load(out / stage / "TRUTH_OPENED.json")
        assert opened["manifest_sha256"] == load(out / stage / "PREDICTION_COMMITMENT.json")["manifest_sha256"]
    assert (out / "validate/SCORE_STARTED.json").exists() and (out / "test/referee/SCORE_STARTED.json").exists()
    freeze = load(out / "freeze/frozen.json")
    assert freeze["selected"] == summary["validate"]["selected"] and freeze["family"] == "no_fitted_state"
    assert freeze["test_access"]["test_stage_directory_exists"] is False
    assert freeze["test_access"]["test_traces_listed_on_disk_at_freeze"] == summary["test"]["declared"]      # said plainly: they were on disk
    assert freeze["test_access"]["strongest_seal_no_test_trace_on_disk"] is False
    assert freeze["pins"]["config_sha256"] == E.CONFIG_SHA and freeze["pins"]["code_sha256"]["frozen"] == E.FROZEN_PINS["source_sha256"]
    assert freeze["all_training_state_refit_sha256"] == freeze["pins"]["regional_states"]["all_training"]["model_state_sha256"]
    assert freeze["referee"]["can_seal_the_test_bundle"] is True and freeze["referee"]["sha256_at_freeze"] == E.referee_hashes()
    # every planted inclusion failure of the generator is logged with its rule, per part
    assert summary["library"]["excluded_by_rule"] == {"2": 1, "4": 1, "6": 1}
    assert summary["validate"]["excluded_by_rule"] == {"1": 1, "2": 1, "3": 1, "5": 1}
    assert summary["test"]["excluded_by_rule"] == {"2": 1, "4": 1, "5": 1, "6": 1}        # four animals, four planted failures
    assert summary["test"]["declared"] - summary["test"]["scored"] == 4 and summary["validate"]["declared"] - summary["validate"]["scored"] == 3
    assert summary["test"]["admission_read_forecasts"] == 0 and summary["validate"]["admission_read_forecasts"] == 0
    assert summary["pair_cache"]["proof"]["frozen_pair_costs_bit_identical"] and summary["pair_cache"]["proof"]["referee_metric_bit_identical"]
    assert summary["validate"]["official_copies"] == {"source": "referee copying_baselines", "api_verified": True,
                                                      "forecasts_compared": 4 * summary["validate"]["declared"], "replica_differences": 0}


def test_a_forecast_changed_after_the_commitment_is_refused_before_any_truth_is_opened(built, tmp_path):
    directory = clone(built, tmp_path / "w")
    engine_of(directory, "validate").predict()
    base = directory / "engine/validate"
    manifest = load(base / "prediction_manifest.json")
    victim, donor = manifest["predictions"][0], manifest["predictions"][-1]
    shutil.copy(base / donor["path"], base / victim["path"])
    with pytest.raises(E.EngineError, match="changed after the commitment"):
        engine_of(directory, "validate").admit()
    assert not (base / "TRUTH_OPENED.json").exists() and not (base / "metrics.json").exists()
    with pytest.raises(E.OnceOnly):                              # and the forecasts cannot simply be redone
        engine_of(directory, "validate").predict()


def test_validation_is_bound_to_the_method_before_truth_is_opened(built, tmp_path, monkeypatch):
    """R10: code, configuration, sources, training traces and library, pair cache and state bytes of the committed
    forecasts are compared with the method on disk BEFORE the opened marker is written."""
    directory = clone(built, tmp_path / "w")
    engine_of(directory, "validate").predict()
    base = directory / "engine/validate"
    # (a) the engine or frozen code changed
    real_hashes = E.code_hashes()
    monkeypatch.setattr(E, "code_hashes", lambda root=E.ROOT: {**real_hashes, "engine": {"analysis/axon_engine.py": "0" * 64}})
    with pytest.raises(E.EngineError, match=r"method changed since these forecasts were committed \(code_sha256"):
        engine_of(directory, "validate").admit()
    monkeypatch.undo()
    # (b) a setting changed
    monkeypatch.setattr(E, "CONFIG_SHA", "f" * 64)
    with pytest.raises(E.EngineError, match="config"):
        engine_of(directory, "validate").admit()
    monkeypatch.undo()
    # (c) the bytes of the regional state changed (same parameters, another file)
    state_path = directory / "engine/states/all_training.proof.json"
    original = state_path.read_text()
    state_path.write_text(json.dumps(json.loads(original), indent=3))
    with pytest.raises(E.EngineError, match="regional_states"):
        engine_of(directory, "validate").admit()
    state_path.write_text(original)
    assert not (base / "TRUTH_OPENED.json").exists()
    # (d) the training library changed after the forecasts were committed (the reviewer's probe)
    animal = E.Cohort(E.synthetic_world(directory)).part_animals["train"][0]
    E.scramble_outcomes(directory, animals={animal})
    with pytest.raises(E.EngineError, match="training library"):
        engine_of(directory, "validate").admit()
    with pytest.raises(E.EngineError, match="training library"):
        engine_of(directory, "validate").score()
    assert not (base / "TRUTH_OPENED.json").exists() and not (base / "SCORE_STARTED.json").exists() and not (base / "truth").exists()


def test_a_stale_or_partial_library_is_never_used_for_held_out_forecasts(built, tmp_path):
    """The reviewer's probe: an ingestion log that holds only part of the training traces listed on disk."""
    directory = clone(built, tmp_path / "w")
    log_path = directory / "engine/ingestion_log.json"
    record = load(log_path)
    whole = len(record["neurons"])
    for neuron in sorted(record["neurons"])[::3]:
        del record["neurons"][neuron]
    log_path.write_text(json.dumps(record))
    stale = engine_of(directory, "validate")
    assert "training traces are listed" in stale.library_problem()
    with pytest.raises(E.EngineError, match="run `ingest`"):
        stale.library
    commit = engine_of(directory, "validate").predict()          # predict ingests again instead of forecasting from the stale library
    manifest = load(directory / "engine/validate/prediction_manifest.json")
    traces = manifest["training_traces"]
    assert commit["status"] == "committed" and traces["listed"] == whole and traces["included"] + traces["excluded"] == whole
    assert len(manifest["training_ids"]) == traces["included"] and engine_of(directory, "validate").library_problem() is None
    reference = clone(built, tmp_path / "ref")
    engine_of(reference, "validate").predict()
    assert arrays(manifest) == arrays(load(reference / "engine/validate/prediction_manifest.json"))


def test_a_started_validation_scoring_is_never_repeated_or_resumed(built, tmp_path):
    directory = clone(built, tmp_path / "w")
    output = validated(directory)
    base = directory / "engine/validate"
    assert output["scored_by"].startswith("engine") and (base / "SCORE_STARTED.json").exists()
    with pytest.raises(E.OnceOnly):
        engine_of(directory, "validate").score()
    # the reviewer's probe: the completion marker is removed and the truth traces replaced; nothing runs again
    (base / "SCORE_COMPLETE.json").unlink()
    before = (base / "metrics.json").read_bytes()
    with pytest.raises(E.OnceOnly, match="never repeated or resumed"):
        engine_of(directory, "validate").score()
    assert replace_traces(directory, "validation") > 0
    for action in (engine_of(directory, "validate").score, engine_of(directory, "validate").admit, engine_of(directory, "validate").predict):
        with pytest.raises(E.EngineError):
            action()
    assert (base / "metrics.json").read_bytes() == before
    second = clone(built, tmp_path / "w2")
    engine_of(second, "validate").predict()
    engine_of(second, "validate").admit()
    (second / "engine/validate/SCORE_STARTED.json").write_text("{}")       # a scoring that died right after it started
    with pytest.raises(E.OnceOnly, match="never repeated or resumed"):
        engine_of(second, "validate").score()
    assert not (second / "engine/validate/metrics.json").exists()


def test_an_interrupted_admission_is_completed_only_on_the_same_commitment_and_trace_files(built, tmp_path):
    directory = clone(built, tmp_path / "w")
    commit = engine_of(directory, "validate").predict()
    base = directory / "engine/validate"
    engine = engine_of(directory, "validate")
    listing = E.digest(engine.heldout_listing())
    (base / "TRUTH_OPENED.json").write_text(json.dumps({"stage": "validate", "manifest_sha256": "f" * 64, "truth_listing_sha256": listing}))
    with pytest.raises(E.OnceOnly, match="another commitment or other trace files"):
        engine_of(directory, "validate").admit()
    (base / "TRUTH_OPENED.json").write_text(json.dumps({"stage": "validate", "manifest_sha256": commit["manifest_sha256"], "truth_listing_sha256": "0" * 64}))
    with pytest.raises(E.OnceOnly, match="another commitment or other trace files"):
        engine_of(directory, "validate").admit()
    (base / "TRUTH_OPENED.json").write_text(json.dumps({"stage": "validate", "manifest_sha256": commit["manifest_sha256"], "truth_listing_sha256": listing}))
    (base / "truth").mkdir()
    (base / "truth/leftover.npz").write_bytes(b"half written")
    summary = engine_of(directory, "validate").admit()            # died before completing: the same mechanical step, same inputs
    assert summary["resumed_after_interruption"] and not (base / "truth/leftover.npz").exists()
    with pytest.raises(E.OnceOnly):
        engine_of(directory, "validate").admit()
    # truth files changed after the admission: scoring refuses and starts nothing
    assert replace_traces(directory, "validation") > 0
    with pytest.raises(E.EngineError, match="trace files changed after the admission"):
        engine_of(directory, "validate").score()
    assert not (base / "SCORE_STARTED.json").exists()


def test_the_freeze_refuses_a_method_other_than_the_one_that_was_validated(built, tmp_path, monkeypatch):
    directory = clone(built, tmp_path / "w")
    validated(directory)
    changed = json.loads(json.dumps(E.CONFIG))
    changed["engine"]["k_neighbours"] = 8                       # a setting edited after validation
    monkeypatch.setattr(E, "CONFIG", changed)
    monkeypatch.setattr(E, "CONFIG_SHA", E.digest(changed))
    with pytest.raises(E.EngineError, match="method changed since these forecasts were committed"):
        engine_of(directory, "validate").freeze()
    assert not (directory / "engine/freeze/frozen.json").exists()
    monkeypatch.undo()
    record = engine_of(directory, "validate").freeze()          # the validated method itself can be frozen
    assert record["validation_method_equals_frozen_method"]
    assert set(record["pins"]["regional_states"]["all_training"]) >= {"model_state_sha256", "state_file_sha256", "proof_file_sha256"}
    assert record["pins"]["pair_cache_sha256"] and record["pins"]["training_traces"]["listing_sha256"]
    assert record["git"]["pinned_dependencies_not_reproduced_by_the_commit"] is not None
    with pytest.raises(E.OnceOnly):
        engine_of(directory, "validate").freeze()


def test_a_consistently_forged_regional_state_does_not_survive_the_refit_at_the_freeze(built, tmp_path):
    """The reviewer's forgery: another model under the all-training name with every stored hash rewritten to agree."""
    directory = clone(built, tmp_path / "w")
    states = directory / "engine/states"
    genuine, other = load(states / "all_training.json"), load(states / "oof_0.json")
    other["training_ids_sha256"] = genuine["training_ids_sha256"]
    (states / "all_training.json").write_text(E.canonical(other))
    sha = E.sha_file(states / "all_training.json")
    assert sha == T._digest(other) != T._digest(genuine)
    proof = load(states / "all_training.proof.json")
    proof["model_state_sha256"] = sha
    (states / "all_training.proof.json").write_text(E.canonical(proof))
    meta = load(states / "all_training.meta.json")
    meta["entry"]["model_state_sha256"] = sha
    (states / "all_training.meta.json").write_text(json.dumps(meta))
    assert engine_of(directory, "validate").states.get("all_training")          # every stored hash agrees: the loader accepts it
    validated(directory)
    with pytest.raises(E.EngineError, match="not what a fresh fit"):
        engine_of(directory, "validate").freeze()
    assert not (directory / "engine/freeze/frozen.json").exists()


def test_test_stage_needs_the_freeze_and_the_freeze_needs_a_finished_validation(built, tmp_path):
    directory = clone(built, tmp_path / "w")
    for action in (engine_of(directory, "test").predict, engine_of(directory, "test").admit, engine_of(directory, "test").seal,
                   engine_of(directory, "test").score):
        with pytest.raises(E.FreezeMismatch, match="no freeze record"):
            action()
    with pytest.raises(E.EngineError, match="completed validation"):
        engine_of(directory, "validate").freeze()
    with pytest.raises(E.EngineError, match="no committed forecast bundle"):
        engine_of(directory, "validate").admit()
    with pytest.raises(E.EngineError, match="only the test is sealed"):
        engine_of(directory, "validate").seal()
    engine_of(directory, "validate").predict()
    with pytest.raises(E.EngineError, match="have not been admitted"):
        engine_of(directory, "validate").score()
    with pytest.raises(E.EngineError, match="completed validation"):
        engine_of(directory, "validate").freeze()


def test_the_freeze_detects_changed_parameters_library_or_gate(frozen_world, tmp_path, monkeypatch):
    copy = clone(frozen_world, tmp_path / "w")
    assert engine_of(copy, "test").verify_freeze()["selected"] == load(copy / "engine/freeze/frozen.json")["selected"]
    gate_path = copy / "engine/validate/engine_gate.json"
    gate = load(gate_path)
    original = gate_path.read_text()
    gate["selected"] = "E1_local_uniform" if gate["selected"] != "E1_local_uniform" else "E2_local_kernel"
    gate_path.write_text(json.dumps(gate))
    with pytest.raises(E.FreezeMismatch, match="engine gate changed"):
        engine_of(copy, "test").verify_freeze()
    gate_path.write_text(original)
    monkeypatch.setattr(E, "CONFIG_SHA", "e" * 64)
    with pytest.raises(E.FreezeMismatch, match="config"):
        engine_of(copy, "test").verify_freeze()
    monkeypatch.undo()
    animal = E.Cohort(E.synthetic_world(copy)).part_animals["train"][0]
    E.scramble_outcomes(copy, animals={animal})
    with pytest.raises(E.FreezeMismatch, match="training library"):
        engine_of(copy, "test").verify_freeze()


@real_only
def test_real_freeze_is_impossible_today_and_held_out_forecasts_need_a_final_protocol(tmp_path):
    engine = E.Engine(E.real_world(out=tmp_path, legacy_only=True), "validate", workers=1, log=None)
    with pytest.raises(E.EngineError, match="completed validation"):
        engine.freeze()
    assert engine.protocol_state()["draft"] is True               # reports/axon_engine_protocol.md is a draft until validation data exist


# ------------------------------------------------------------------------------------------
# plan amendment 2: a neuron the source does not serve (rule 8); no test trace before the freeze;
# a listing that is still being written; second copies
# ------------------------------------------------------------------------------------------


def mark_unavailable(directory, neurons, *, keep_path=True):
    """Rewrite listing rows as the download listing records a neuron the source does not serve: frame
    'unavailable_at_source', no file on disk (with or without a path in the file column)."""
    path = directory / "data/traces/MANIFEST.tsv"
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines[1:], 1):
        fields = line.split("\t")
        if fields[1] in neurons:
            (directory / "data/traces" / fields[0]).unlink()
            fields[0], fields[2], fields[3], fields[4], fields[5] = (fields[0] if keep_path else ""), "0", "-", "HTTP 404 at the source", E.UNAVAILABLE_FRAME
            lines[i] = "\t".join(fields)
    path.write_text("\n".join(lines) + "\n")


@pytest.mark.parametrize("keep_path", [True, False])
def test_a_neuron_the_source_does_not_serve_is_listed_excluded_under_rule_8_and_counted(built, tmp_path, keep_path):
    directory = clone(built, tmp_path / "w")
    cohort = E.Cohort(E.synthetic_world(directory))
    donor, query = cohort.part_neurons["train"][5], cohort.part_neurons["validation"][0]
    mark_unavailable(directory, {donor, query}, keep_path=keep_path)
    engine = engine_of(directory, "validate")
    assert engine.store.available("train")[0] == len(cohort.part_neurons["train"])      # listed, not silently absent
    summary = engine.build_library()
    assert summary["exclusions"][donor]["rules"] == ["8"] and summary["excluded_by_rule"]["8"] == 1
    assert "HTTP 404 at the source" in summary["exclusions"][donor]["details"]["8"]
    lost = summary["unavailable_at_source"]
    assert lost["neurons"] == 1 and lost["by_source_label"] == {cohort.rows[donor].source_label: 1} and lost["neuron_ids"] == [donor]
    assert lost["scored_cohort_smaller_than_locked_notice"] is False                    # a training neuron is simply not a donor
    assert donor not in engine.library["ids"]
    engine.build_cache()
    engine.states.build_many(["all_training"])
    engine.predict()
    admitted = engine_of(directory, "validate").admit()
    assert admitted["excluded"][query]["rules"] == ["8"] and admitted["excluded_by_rule"]["8"] == 1 and query not in admitted["admitted_ids"]
    report = admitted["unavailable_at_source"]
    assert report["neurons"] == 1 and report["by_source_label"] == {cohort.rows[query].source_label: 1}
    assert report["share"] == pytest.approx(1 / len(cohort.part_neurons["validation"]))
    assert report["scored_cohort_smaller_than_locked_notice"] == (report["share"] > 0.05) and report["brain_samples_lost_entirely"] == []
    output = engine_of(directory, "validate").score()
    assert output["unavailable_at_source"] == report
    tables = (directory / "engine/validate/tables.md").read_text()
    assert "Unavailable at the source (rule 8): 1" in tables
    assert tables.splitlines()[2].startswith("**The scored cohort is smaller than the locked one") == report["scored_cohort_smaller_than_locked_notice"]
    assert "plan amendment 2" in E.RULE_TEXT["8"] and E.UNAVAILABLE_FRAME == "unavailable_at_source"


def test_rule_8_reporting_duty_a_whole_held_out_brain_sample_or_more_than_five_percent(built):
    cohort = E.Cohort(E.synthetic_world(built))
    part = cohort.part_neurons["test"]
    lost = lambda neurons: {n: {"rules_violated": ["8"] if n in neurons else []} for n in part}
    animal = cohort.part_animals["test"][0]
    whole = E.rule8_report(lost(set(cohort.animal_neurons[animal])), cohort, "test")
    assert whole["brain_samples_lost_entirely"] == [animal] and whole["scored_cohort_smaller_than_locked_notice"]
    none = E.rule8_report(lost(set()), cohort, "test")
    assert none["neurons"] == 0 and not none["scored_cohort_smaller_than_locked_notice"] and none["of"] == len(part)
    rows = {f"t{i}": None for i in range(100)}
    fake = type("C", (), {"part_neurons": {"validation": list(rows)}, "part_animals": {"validation": ["a", "b"]},
                          "animal_neurons": {"a": list(rows)[:50], "b": list(rows)[50:]},
                          "rows": {n: type("R", (), {"source_label": "L"})() for n in rows}})()
    five = E.rule8_report({n: {"rules_violated": ["8"] if i < 5 else []} for i, n in enumerate(rows)}, fake, "validation")
    six = E.rule8_report({n: {"rules_violated": ["8"] if i < 6 else []} for i, n in enumerate(rows)}, fake, "validation")
    assert not five["scored_cohort_smaller_than_locked_notice"] and six["scored_cohort_smaller_than_locked_notice"]      # more than 5%
    assert six["by_source_label"] == {"L": 6}


def test_on_the_real_cohort_the_freeze_refuses_while_a_test_trace_is_listed_on_disk(built, tmp_path):
    """Plan amendment 2: no test file may be on disk before the freeze. A synthetic world carries its generated test
    traces from the start (its freeze records that); the same directory treated as the real cohort is refused."""
    import dataclasses
    directory = clone(built, tmp_path / "w")
    validated(directory)
    as_real = dataclasses.replace(E.synthetic_world(directory), kind="real")
    with pytest.raises(E.EngineError, match="no test trace may be on disk before the freeze"):
        E.Engine(as_real, "validate", workers=1, log=None).freeze()
    assert not (directory / "engine/freeze/frozen.json").exists()
    record = engine_of(directory, "validate").freeze()
    assert record["test_access"]["test_traces_listed_on_disk_at_freeze"] > 0 and not record["test_access"]["strongest_seal_no_test_trace_on_disk"]


def test_development_on_the_whole_real_cohort_needs_the_final_training_listing(built, tmp_path):
    """A download in progress lists only part of the training traces: `develop` on the real cohort waits for the final
    listing; the fixed development snapshot (local37) is declared as partial and may run."""
    import dataclasses
    directory = clone(built, tmp_path / "w")
    path = directory / "data/traces/MANIFEST.tsv"
    lines = path.read_text().splitlines()
    missing = [l for l in lines[1:] if l.startswith("train/")][:3]                          # three training rows not listed yet
    path.write_text("\n".join(l for l in lines if l not in missing) + "\n")
    world = E.synthetic_world(directory)
    cohort = E.Cohort(world)
    assert E.TraceStore(world, cohort, "develop").available("train")[0] == len(cohort.part_neurons["train"]) - 3
    as_real = dataclasses.replace(world, kind="real")
    with pytest.raises(E.StageBlocked, match="final training listing"):
        E.TraceStore(as_real, cohort, "develop").require("develop", need_truth=False)
    snapshot = dataclasses.replace(as_real, development_snapshot=True)
    E.TraceStore(snapshot, cohort, "develop").require("develop", need_truth=False)          # declared snapshot: allowed
    E.TraceStore(world, cohort, "develop").require("develop", need_truth=False)             # synthetic worlds are unaffected
    for stage in ("validate", "test"):
        with pytest.raises(E.StageBlocked, match="have not been fetched completely"):
            E.TraceStore(as_real, cohort, stage).require(stage, need_truth=False)
    assert E.real_world(legacy_only=True).development_snapshot and not E.real_world().development_snapshot
    assert len(E.real_world(legacy_only=True).trace_manifests) == 1 and len(E.real_world().trace_manifests) == 2


def test_a_second_copy_of_a_training_trace_is_checked_against_the_chosen_copy(built, tmp_path):
    directory = clone(built, tmp_path / "w")
    path = directory / "data/traces/MANIFEST.tsv"
    lines = path.read_text().splitlines()
    first = next(l for l in lines[1:] if l.startswith("train/")).split("\t")
    neuron = first[1]
    copy_file = first[0].replace(".swc", "_zcopy.swc")
    shutil.copy(directory / "data/traces" / first[0], directory / "data/traces" / copy_file)
    path.write_text("\n".join(lines + ["\t".join([copy_file] + first[1:])]) + "\n")
    engine = engine_of(directory, "develop")
    assert "training traces are listed" in engine.library_problem()                 # a new copy: the library is ingested again
    summary = engine.build_library()
    entry = engine.library["all"][neuron]
    assert entry["source_file"] == first[0] and entry["copies_on_disk"] == 2
    assert entry["copy_check"]["copies"] == 2 and entry["copy_check"]["all_copies_give_the_same_canonical_arbor"]
    assert summary["second_copies"] == {"neurons_with_more_than_one_copy": 1, "copies_that_disagree": []}
    assert summary["fingerprint"] == engine_of(built, "develop").library["fingerprint"]      # the chosen arbors are unchanged
    # the second copy now holds another tree: the disagreement is logged and reported, the chosen copy is still used
    copy_path = directory / "data/traces" / copy_file
    rows = copy_path.read_text().splitlines()
    k = max(i for i, row in enumerate(rows) if not row.startswith("#") and row.split()[1] in ("0", "2"))      # the last axon node
    fields = rows[k].split()
    fields[2] = f"{float(fields[2]) + 77.0:.2f}"
    rows[k] = " ".join(fields)
    copy_path.write_text("\n".join(rows) + "\n")
    lines = path.read_text().splitlines()
    last = lines[-1].split("\t")
    last[2], last[3] = str(copy_path.stat().st_size), E.sha_file(copy_path)
    path.write_text("\n".join(lines[:-1] + ["\t".join(last)]) + "\n")
    summary = engine_of(directory, "develop").build_library()
    assert summary["second_copies"]["copies_that_disagree"] == [neuron] and neuron in engine_of(directory, "develop").library["ids"]


# ------------------------------------------------------------------------------------------
# subgroup reporting (frozen before test access) and wording
# ------------------------------------------------------------------------------------------


def test_subgroup_report_shows_neurons_and_animals_and_marks_every_group_descriptive():
    rows = []
    plan = [("Big", "PVH", 6), ("Small", "HY", 4), ("Small2", "fx", 1), ("", "", 1)]
    for label, soma_label, animals in plan:
        for a in range(animals):
            for k in range(2):
                metrics = {m: {key: 0.0 for key in E.KEPT_METRICS} for m in E.METHODS}
                for m in E.METHODS:
                    metrics[m]["distance_um"], metrics[m]["f1_200"] = (500.0, 0.5) if m == "E2_local_kernel" else (900.0, 0.2)
                rows.append({"neuron": f"{label or 'U'}{a}_{k}", "animal": f"{label or 'U'}{a}", "source_label": label,
                             "soma_label": soma_label, "metrics": metrics})
    report = E.subgroup_report(rows, E.METHODS, "E2_local_kernel")
    assert report["spec_sha256"] == E.SUBGROUP_SPEC_SHA and len(report["per_animal"]) == 12
    assert all(set(a["methods"]) == set(E.METHODS) and a["neurons"] == 2 for a in report["per_animal"])
    by_label = {g["source_label"]: g for g in report["by_source_label"]}
    big, small = by_label["Big"], by_label["Small"]
    assert (big["neurons"], big["animals"], big["sparse_descriptive_only"], big["interpretation"]) == (12, 6, False, "descriptive_only")
    assert (small["neurons"], small["animals"], small["sparse_descriptive_only"], small["interpretation"]) == (8, 4, True, "sparse_descriptive_only")
    assert all(g["descriptive_only"] and "no population claim" in g["note"] for g in report["by_source_label"] + report["by_soma_label"])
    assert not any("population_statement" in g for g in report["by_source_label"])          # no group carries a population statement
    assert big["descriptive_comparison"] == {"selected_better_than_every_control_on_both": True, "controls_not_beaten": []}
    assert "unknown" in by_label and by_label["unknown"]["sparse_descriptive_only"]
    by_soma = {g["soma_label"]: g for g in report["by_soma_label"]}
    assert by_soma["PVH"]["nucleus_assigned"] and by_soma["PVH"]["soma_label_kind"] == "nucleus"
    assert not by_soma["HY"]["nucleus_assigned"] and by_soma["HY"]["soma_label_kind"] == "unassigned_HY"
    assert not by_soma["fx"]["nucleus_assigned"] and by_soma["fx"]["soma_label_kind"] == "unassigned_fibre_tract"
    assert by_soma["unknown"]["soma_label_kind"] == "unknown"
    text = E.render_tables({"label": "x", "scored_queries": 24, "declared_queries": 24, "scored_animals": 12, "declared_animals": 12,
                            "excluded": {}, "excluded_by_rule": {}, "summary": E.aggregate(rows, E.METHODS)["summary"], "subgroups": report})
    assert "sparse: descriptive only" in text and "unassigned HY" in text and "unassigned fibre tract" in text
    assert "carries no population claim" in text


def test_pipeline_metrics_carry_the_frozen_subgroup_report_and_the_development_table_its_paired_differences(pipeline):
    directory, _ = pipeline
    out = Path(E.synthetic_world(directory).out)
    metrics = load(out / "test" / "metrics.json")
    groups = metrics["subgroups"]
    assert groups["spec_sha256"] == E.SUBGROUP_SPEC_SHA
    assert sum(a["neurons"] for a in groups["per_animal"]) == metrics["scored_queries"]
    assert len(groups["per_animal"]) == metrics["scored_animals"]
    models = metrics["methods_scored_by_the_referee"]
    assert all(set(a["methods"]) == set(models) for a in groups["per_animal"])           # the candidate slot and every control, per animal
    for key in ("by_source_label", "by_soma_label"):
        assert sum(g["neurons"] for g in groups[key]) == metrics["scored_queries"]
        assert all({"neurons", "animals", "sparse_descriptive_only", "interpretation"} <= set(g) for g in groups[key])
        assert all(g["sparse_descriptive_only"] == (g["animals"] < 5) for g in groups[key])
    assert metrics["verdict"]["context"] == "test" and "brain sample IDs" in metrics["verdict"]["unit"]
    assert set(metrics["referee_subgroups"]) == {"by_peptide", "by_nucleus"}
    develop = load(out / "develop" / "metrics.json")
    assert develop["development_only"] and set(develop["paired_differences"]["candidates"]) == set(E.CANDIDATES)
    one = develop["paired_differences"]["candidates"]["E2_local_kernel"]["versus_strongest_on_distance"]
    assert {"control", "distance_delta_um", "distance_bootstrap95_um", "animals_better_on_distance", "animals"} <= set(one)
    tables = (out / "develop" / "tables.md").read_text()
    assert "Development data, not a test" in tables and "Paired animal differences" in tables and "no selection is made" in tables
    assert develop["official_copies"]["api_verified"] and develop["official_copies"]["replica_differences"] == 0
    validation = load(out / "validate" / "metrics.json")
    assert validation["scored_by"].startswith("engine") and "bit-identical" in validation["scored_by"]
    assert set(validation["summary"]) == set(E.METHODS) and "Validation gate" in (out / "validate" / "tables.md").read_text()


# ------------------------------------------------------------------------------------------
# diagnostics
# ------------------------------------------------------------------------------------------


def test_the_coverage_diagnostic_uses_training_and_validation_inputs_only(built):
    engine = engine_of(built, "develop")
    before = E.coverage_diagnostic(engine)
    assert set(before) >= {"training_leave_one_animal_out", "validation"} and "test" not in before
    poisoned = engine_of(built, "develop")
    for neuron in poisoned.cohort.part_neurons["test"]:              # any use of a test soma or label would change the result
        row = poisoned.cohort.rows[neuron]
        poisoned.cohort.rows[neuron] = E.NeuronRow(row.neuron, row.animal, row.full_id, "Poison", row.tag, row.soma_label,
                                                   row.hemisphere, row.recon_type, row.part, row.klass, (1.0e9, -1.0e9, 3.0e9))
    assert E.coverage_diagnostic(poisoned) == before
    assert before["training_leave_one_animal_out"]["queries"] == len(engine.cohort.part_neurons["train"])


def test_cost_probe_times_the_frozen_and_referee_metrics_and_the_unmodified_policies(built):
    engine = engine_of(built, "develop")
    engine.prepare_prediction()
    result = E.cost_probe(engine, sizes=(40, 80), pair_sample=20)
    assert result["frozen_and_referee_bit_identical_on_the_timed_pairs"] and result["ms_per_pair_frozen_metric"] > 0
    assert result["projection_4448_training_neurons"]["pairs"] == 9890128
    assert [c["donors"] for c in result["frozen_predict_policy_per_query"]] == [40, 80]
    assert all(c["five_policies_s"] > 0 for c in result["frozen_predict_policy_per_query"])


def test_commands_run_on_a_synthetic_world(built, tmp_path, capsys):
    sys.path.insert(0, str(ROOT / "analysis"))
    import axon_engine_cli as cli
    directory = clone(built, tmp_path / "w")
    assert cli.main(["--world", str(directory), "--workers", "1", "status"]) == 0
    assert cli.main(["--world", str(directory), "--workers", "1", "cache"]) == 0
    assert cli.main(["--world", str(directory), "--workers", "1", "predict", "--stage", "validate", "--family", "no_fitted_state"]) == 0
    assert cli.main(["--world", str(directory), "--workers", "1", "predict", "--stage", "validate"]) == 4          # once only
    assert cli.main(["--world", str(directory), "--workers", "1", "score", "--stage", "validate"]) == 0          # admits, then scores
    assert cli.main(["--world", str(directory), "--workers", "1", "score", "--stage", "validate"]) == 4
    assert cli.main(["--world", str(directory), "--workers", "1", "freeze"]) == 0
    assert cli.main(["--world", str(directory), "--workers", "1", "predict", "--stage", "test"]) == 0
    assert cli.main(["--world", str(directory), "--workers", "1", "score", "--stage", "test"]) == 0              # admit, seal, referee
    assert cli.main(["--world", str(directory), "--workers", "1", "seal"]) == 4
    text = capsys.readouterr().out
    assert "independent referee" in text and "REFUSED (OnceOnly)" in text
    assert load(directory / "engine/freeze/frozen.json")["family"] == "no_fitted_state"
    assert (directory / "engine/test/referee/SCORE_COMPLETE.json").exists()


# ------------------------------------------------------------------------------------------
# fix round of 2026-11-11 (late items L2, L3, L6; review finding 1). Each test below was written before its fix and
# failed then; the failing output is kept in data/interim/synarch_axon_bridge/regression_before_fix.txt
# ------------------------------------------------------------------------------------------


def test_L2_the_engine_hashes_the_same_three_referee_sources_as_the_referee_and_pins_them(built, tmp_path):
    assert set(E.REFEREE_FILES) == {"evaluator", "geometry_math", "training_distance_cache"}
    assert E.REFEREE_FILES["training_distance_cache"] == "analysis/prediction_distance_cache.py"
    assert E.referee().source_hashes() == E.referee_hashes()          # like with like: the guard compares three entries
    directory = clone(built, tmp_path / "w")
    engine_of(directory, "validate").predict()
    manifest = load(directory / "engine/validate/prediction_manifest.json")
    assert manifest["referee_sha256"] == E.referee_hashes() and len(manifest["referee_sha256"]) == 3
    commit = load(directory / "engine/validate/PREDICTION_COMMITMENT.json")
    assert commit["referee_sha256"] == E.referee_hashes()


def test_L3_development_honours_the_pair_budget_even_for_its_first_copy_group(built, tmp_path, monkeypatch):
    """L3: with a budget smaller than the first animal's pool, no copy group is sent to the referee in development, and
    the manifest says so. Only a redundant development check changes: every forecast array is the same."""
    reference = clone(built, tmp_path / "ref")
    engine_of(reference, "develop", develop_per_animal=1).predict()
    directory = clone(built, tmp_path / "w")
    monkeypatch.setitem(E.CONFIG["official_copies"], "develop_api_pair_budget", 1)
    monkeypatch.setattr(E, "CONFIG_SHA", E.digest(E.CONFIG))
    Ref = E.referee()
    calls = []
    real = Ref.copying_baselines
    monkeypatch.setattr(Ref, "copying_baselines", lambda *a, **k: calls.append(1) or real(*a, **k))
    engine_of(directory, "develop", develop_per_animal=1).predict()
    official = load(directory / "engine/develop/prediction_manifest.json")["official_copies"]
    assert calls == [], f"{len(calls)} copy groups were sent to the referee although every group is over the budget"
    assert official["api_verified"] is False and official["animals_with_referee_copies"] == []
    assert official["api_verification"] == "no development API verification (pair budget)"
    assert arrays(load(directory / "engine/develop/prediction_manifest.json")) == arrays(load(reference / "engine/develop/prediction_manifest.json"))


def test_L6_the_engine_refuses_the_referee_unless_the_four_thread_variables_are_one(monkeypatch):
    E.referee()
    for name in THREADS:
        monkeypatch.setenv(name, "4")
        with pytest.raises(E.EngineError, match=name):
            E.referee()
        monkeypatch.delenv(name)
        with pytest.raises(E.EngineError, match=name):
            E.referee()
        monkeypatch.setenv(name, "1")
    E.referee()


def test_L6_the_command_line_sets_the_four_thread_variables_before_numpy_is_imported(built, tmp_path):
    directory = clone(built, tmp_path / "w")
    env = {k: v for k, v in os.environ.items() if k not in THREADS}
    probe = ("import sys, runpy; sys.argv = ['axon_engine_cli.py', '--help']; import os\n"
             "try:\n    runpy.run_path('analysis/axon_engine_cli.py', run_name='__main__')\nexcept SystemExit:\n    pass\n"
             "print('THREADS', [os.environ.get(v) for v in %r]); print('NUMPY_FIRST', 'numpy' in sys.modules)" % (list(THREADS),))
    out = subprocess.run([sys.executable, "-I", "-c", probe], cwd=ROOT, env=env, capture_output=True, text=True, timeout=300)
    assert "THREADS ['1', '1', '1', '1']" in out.stdout, out.stdout + out.stderr
    text = (ROOT / "analysis/axon_engine_cli.py").read_text()
    assert text.index("VECLIB_MAXIMUM_THREADS") < text.index("import axon_engine as E")
    done = subprocess.run([sys.executable, "-I", "analysis/axon_engine_cli.py", "--world", str(directory), "--workers", "1", "status"],
                          cwd=ROOT, env=env, capture_output=True, text=True, timeout=600)
    assert done.returncode == 0, done.stdout[-2000:] + done.stderr[-2000:]


def audit_bundle(bundle, minutes=None):
    """The lead's audit command, exactly as the referee's contract states it (one worker here)."""
    command = [sys.executable, "-I", "analysis/prediction_distance_cache.py", "--base", str(bundle), "--declaration",
               "draft_declaration.json", "--workers", "1"] + ([] if minutes is None else ["--minutes", str(minutes)])
    out = subprocess.run(command, cwd=ROOT, env={**os.environ, **THREADS}, capture_output=True, text=True, timeout=1200)
    assert out.returncode == 0, out.stdout[-2000:] + out.stderr[-3000:]
    return json.loads(out.stdout.strip().splitlines()[-1])


def test_L6_prepare_audit_and_reuse_the_referees_training_distance_cache_end_to_end(built, tmp_path, monkeypatch):
    """L6 on a synthetic world: the engine prepares a training-only bundle and a draft declaration from its completed
    native cache; the referee's helper audits every pair; the engine's referee calls on the EXACT declared training pool
    then use the audited cache, with bit-identical copies, and every other pool never does. A corrupted matrix or another
    donor pool is refused; the receipt is bound in the commitment, the freeze and the referee's declaration."""
    directory = clone(built, tmp_path / "w")
    dev = engine_of(directory, "develop")
    assert dev.referee_cache_state()["status"] == "absent"
    prepared = dev.prepare_referee_cache()
    bundle = Path(prepared["bundle"])
    draft = load(bundle / "draft_declaration.json")
    descriptor = draft["training_distance_cache"]
    assert descriptor == E.REFEREE_CACHE_DESCRIPTOR and descriptor["schema"] == "referee_training_distance_cache_v1"
    paths = [descriptor[k] for k in ("index_path", "complete_path", "distance_path", "f1_path", "mask_path", "audit_path")]
    assert set(paths) <= set(draft["method_artifacts"])
    assert draft["method_artifacts"][descriptor["audit_path"]] == E.REFEREE_CACHE_RECEIPT_PLACEHOLDER
    assert all(E.sha_file(bundle / p) == draft["method_artifacts"][p] for p in paths[:-1])
    train_ids = [r["id"] for r in draft["splits"]["train"]]
    assert train_ids == dev.library["ids"] and set(draft["splits"]) == {"train", "validation", "test"}
    held = set(dev.cohort.part_neurons["validation"]) | set(dev.cohort.part_neurons["test"])
    assert not held & set(train_ids) and not any((bundle / "training" / f"{n}.npz").exists() for n in held)
    assert all(set(r) <= {"id", "animal", "soma_um", "peptide", "nucleus"} for p in ("validation", "test") for r in draft["splits"][p])
    assert dev.referee_cache_state()["status"] == "prepared_audit_pending"
    assert dev.prepare_referee_cache()["status"] == "prepared_audit_pending"          # idempotent; nothing rebuilt
    # the audit, by the referee's own helper and command line; then the receipt unlocks the loader
    result = audit_bundle(bundle)
    assert result["status"] == "complete" and result["pairs_verified"] == dev.cache.pairs
    val = engine_of(directory, "validate")
    state = val.referee_cache_state()
    assert state["status"] == "complete" and state["receipt_sha256"] == E.sha_file(bundle / descriptor["audit_path"])
    # bit-identical copies with and without the cache, on the exact declared training pool
    Ref = E.referee()
    scope = val.scope(None)
    queries = val.cohort.part_neurons["validation"]
    donors, cached = val.referee_training_donors(scope)
    assert cached is not None and cached["cache_base"] == bundle
    plain = Ref.copying_baselines(donors, [val.known_row(n) for n in queries])
    reused = Ref.copying_baselines(donors, [val.known_row(n) for n in queries], **cached)
    assert plain["donor_proof"] == reused["donor_proof"]
    assert all(np.array_equal(a, b) for m in plain["predictions"] for q in plain["predictions"][m]
               for a, b in zip(plain["predictions"][m][q], reused["predictions"][m][q]))
    # a smaller or another donor pool is refused by the referee's loader; a development subpool never gets the cache
    with pytest.raises(ValueError, match="EXACT sorted frozen TRAIN pool"):
        Ref.copying_baselines(donors[1:], [val.known_row(n) for n in queries], **cached)
    animal = val.cohort.rows[dev.library["ids"][0]].animal
    assert engine_of(directory, "develop").referee_training_donors(engine_of(directory, "develop").scope(animal))[1] is None
    assert engine_of(directory, "develop").referee_training_donors(engine_of(directory, "develop").scope(None))[1] is None
    # a corrupted matrix is refused (by the engine before the referee, and by the referee's loader itself)
    matrix = bundle / descriptor["distance_path"]
    original = matrix.read_bytes()
    D = np.load(matrix)
    D[0, 1] = D[1, 0] = D[0, 1] + 1e-9
    np.save(matrix, D)
    with pytest.raises(E.EngineError, match="audited training distance cache"):
        engine_of(directory, "validate").referee_cache_state()
    with pytest.raises(ValueError, match="changed after completion"):
        Ref.copying_baselines(donors, [val.known_row(n) for n in queries], **cached)
    matrix.write_bytes(original)
    # through the engine's own steps: a pool above the inline limit is checked by copies-check WITH the audited cache
    monkeypatch.setitem(E.CONFIG["official_copies"], "heldout_inline_api_pair_limit", 10)
    monkeypatch.setattr(E, "CONFIG_SHA", E.digest(E.CONFIG))
    calls = []
    real = Ref.copying_baselines
    monkeypatch.setattr(Ref, "copying_baselines", lambda *a, **k: calls.append(k.get("training_distance_cache")) or real(*a, **k))
    commit = engine_of(directory, "validate").predict()
    assert commit["official_copies_api_verified"] is False
    assert commit["referee_training_cache_receipt_sha256"] == state["receipt_sha256"]
    pins = load(directory / "engine/validate/prediction_manifest.json")["pins"]["referee_training_cache"]
    assert pins["status"] == "complete" and pins["receipt_sha256"] == state["receipt_sha256"]
    record = engine_of(directory, "validate").copies_check()
    assert record["all_equal"] and record["referee_training_cache_receipt_sha256"] == state["receipt_sha256"]
    assert calls[-1] == descriptor
    # the receipt is bound: changed after the commitment, admit refuses before any truth
    receipt = bundle / descriptor["audit_path"]
    kept = receipt.read_bytes()
    receipt.write_bytes(kept + b" ")
    with pytest.raises(E.EngineError, match="referee_training_cache|audited training distance cache"):
        engine_of(directory, "validate").admit()
    assert not (directory / "engine/validate/TRUTH_OPENED.json").exists()
    receipt.write_bytes(kept)
    engine_of(directory, "validate").admit()
    engine_of(directory, "validate").score()
    frozen = engine_of(directory, "validate").freeze()
    assert frozen["pins"]["referee_training_cache"]["receipt_sha256"] == state["receipt_sha256"]
    assert frozen["referee_training_cache_receipt_sha256"] == state["receipt_sha256"]
    engine_of(directory, "test").predict()
    engine_of(directory, "test").copies_check()
    engine_of(directory, "test").admit()
    calls.clear()
    seal = engine_of(directory, "test").seal()
    assert seal["referee_training_cache_receipt_sha256"] == state["receipt_sha256"]
    assert calls and all(c == descriptor for c in calls)               # inside seal_predictions: the declared cache
    declaration = load(directory / "engine/test/referee/declaration.json")
    assert declaration["training_distance_cache"] == descriptor
    assert all(p in declaration["method_artifacts"] for p in paths)
    assert declaration["method_artifacts"][descriptor["audit_path"]] == state["receipt_sha256"]
    assert engine_of(directory, "test").score()["scored_by"].startswith("independent referee")


def test_L6_on_the_real_cohort_a_validation_commitment_is_refused_while_the_audit_receipt_is_missing(built, tmp_path):
    import dataclasses
    directory = clone(built, tmp_path / "w")
    protocol = tmp_path / "protocol.md"
    protocol.write_text("# a final protocol (no draft line)\n")
    as_real = dataclasses.replace(E.synthetic_world(directory), kind="real", protocol=str(protocol))
    for stage in ("validate", "test"):
        engine = E.Engine(as_real, stage, workers=1, log=None, family="no_fitted_state")
        with pytest.raises(E.EngineError, match="audit receipt"):
            engine.require_referee_training_cache()
    with pytest.raises(E.EngineError, match="audit receipt"):
        E.Engine(as_real, "validate", workers=1, log=None, family="no_fitted_state").predict()
    assert not (directory / "engine/validate/PREDICTION_COMMITMENT.json").exists()
    E.Engine(E.synthetic_world(directory), "validate", workers=1, log=None).require_referee_training_cache()   # synthetic: optional


def test_finding1_the_library_is_keyed_to_its_ingestion_code_not_to_the_whole_engine_file(built, tmp_path):
    """Review finding 1: a library made by an engine file whose INGESTION code is the same is accepted; a library whose
    ingestion code differs is refused."""
    directory = clone(built, tmp_path / "w")
    log_path = directory / "engine/ingestion_log.json"
    record = load(log_path)
    record["engine_sha256"] = "0" * 64                     # another version of the engine file, same ingestion code
    log_path.write_text(json.dumps(record))
    assert engine_of(directory, "validate").library_problem() is None
    record["ingestion_code_sha256"] = "1" * 64              # other ingestion code: refused
    log_path.write_text(json.dumps(record))
    assert "ingestion code" in engine_of(directory, "validate").library_problem()
    del record["ingestion_code_sha256"]                     # an old record without the key, from an unknown engine file
    log_path.write_text(json.dumps(record))
    assert "ingestion code" in engine_of(directory, "validate").library_problem()


def test_finding1_the_committed_engine_has_the_same_ingestion_code_and_its_real_library_is_accepted():
    """The engine that built the full-cohort training library (git f46a87e, file SHA-256 e575ddc3...) has exactly the
    current ingestion code, function by function; its ingestion log is therefore accepted without a re-ingest. The real
    ingestion log is only READ here (its record); no trace and no library arbor is opened."""
    text = subprocess.run(["git", "show", "f46a87e:analysis/axon_engine.py"], cwd=ROOT, capture_output=True, check=True).stdout
    assert E.sha_bytes(text) == "e575ddc32d5d5e71b30a88744a012e509ee68ce6953f6d5f37596cd676df7d23"
    current = E.ingestion_code_sha((ROOT / "analysis/axon_engine.py").read_text())
    assert E.ingestion_code_sha(text.decode()) == current == E.INGESTION_CODE_SHA
    assert E.LEGACY_INGESTION_ENGINES == {E.sha_bytes(text): current}
    real_log = ROOT / "data/interim/axon_engine/real/ingestion_log.json"
    if real_log.exists():
        record = load(real_log)
        # the log may have been rewritten by a later engine (re-ingested at 12:37 on 2026-11-11 by d790bef, e4201e5c...):
        # what binds it is the ingestion code, which must be exactly the current one
        assert record["ingestion_code_sha256"] == current and record["rules"] == E.jsonable(E.RULE_TEXT)
        if record["engine_sha256"] not in E.LEGACY_INGESTION_ENGINES:
            later = subprocess.run(["git", "show", "d790bef:analysis/axon_engine.py"], cwd=ROOT, capture_output=True, check=True).stdout
            assert record["engine_sha256"] == E.sha_bytes(later) and E.ingestion_code_sha(later.decode()) == current
        assert E.ingestion_record_problem(record) is None


def test_finding1_the_fingerprinted_ingestion_code_is_closed_under_what_it_calls():
    import ast
    import builtins
    text = (ROOT / "analysis/axon_engine.py").read_text()
    tree = ast.parse(text)
    top = {n.name: n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
    covered = set(E.INGESTION_FUNCTIONS) | set(E.INGESTION_CONSTANTS)
    imported = {a.asname or a.name.split(".")[0] for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom)) for a in n.names}
    allowed = covered | imported | set(dir(builtins)) | {"G", "T", "CONFIG", "EngineError"}
    for name in E.INGESTION_FUNCTIONS:
        used = {n.id for n in ast.walk(top[name]) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
        local = {a.arg for n in ast.walk(top[name]) if isinstance(n, ast.arguments) for a in n.args + n.kwonlyargs}
        local |= {n.id for n in ast.walk(top[name]) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store)}
        local |= {n.name for n in ast.walk(top[name]) if isinstance(n, ast.FunctionDef)}
        local |= {n.name for n in ast.walk(top[name]) if isinstance(n, ast.ExceptHandler) and n.name}
        assert used - local <= allowed, (name, sorted(used - local - allowed))
    uses_config = [n for n in E.INGESTION_FUNCTIONS if "CONFIG[" in ast.get_source_segment(text, top[n])]
    for name in uses_config:
        keys = {s for s in ast.get_source_segment(text, top[name]).split('CONFIG["inclusion"]["')[1:]}
        assert {k.split('"')[0] for k in keys} <= set(E.INGESTION_CONFIG_KEYS), name


# ------------------------------------------------------------------------------------------
# second attack round (round 2): the referee cache receipt, the thread guard
# ------------------------------------------------------------------------------------------


def distance_cache_helper():
    """The referee's own cache helper, imported read-only (as an attacker would use it to write a receipt)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("round2_pdc", ROOT / "analysis/prediction_distance_cache.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def forge_receipt(directory, bundle):
    """The attackers' route: change D.npy (one donor made close to everybody), make COMPLETE.json agree with it, and write a
    'complete' audit receipt with the helper's own binding and digest, WITHOUT running the audit."""
    C, desc = distance_cache_helper(), E.REFEREE_CACHE_DESCRIPTOR
    path = bundle / desc["distance_path"]
    D = np.load(path)
    k = D.shape[0] // 2
    D[k, :] *= 0.01
    D[:, k] *= 0.01
    D[k, k] = 0.0
    np.save(path, np.ascontiguousarray(D))
    complete = load(bundle / desc["complete_path"])
    complete["D_sha256"] = C.G.sha(path)
    complete["cache_sha256"] = C.digest({key: complete[key] for key in ("schema", "library_fingerprint", "ids_sha256", "n", "pairs",
                                                                        "metric", "index_sha256", "D_sha256", "F_sha256")})
    (bundle / desc["complete_path"]).write_text(C.canonical(complete))
    val = engine_of(directory, "validate")
    donors = [{**val.known_row(n), "arbor": E.as_tree(val.donor(n)["arbor"]), "arbor_path": f"training/{n}.npz",
               "arbor_sha256": E.sha_file(bundle / f"training/{n}.npz")} for n in val.scope(None)["ids"]]
    _, binding, _ = C._binding(bundle, desc, donors)
    n, total = binding["n"], binding["pairs"]
    body = {"schema": C.AUDIT_SCHEMA, "status": "complete", "scope": "ALL_upper_triangle_D_bitwise;F_hash_and_structure_only",
            "binding": binding, "pairs_verified": total, "pairs_expected": total,
            "row_pairs_verified": [n - i - 1 for i in range(n)], "bitwise_mismatches": 0}
    (bundle / desc["audit_path"]).write_text(C.canonical({**body, "receipt_sha256": C.digest(body)}))
    return C.G.sha(path)


def test_round2_a_forged_audit_receipt_for_a_changed_matrix_is_refused_and_nothing_is_committed(built, tmp_path):
    """Round-2 finding (major): the receipt is a self-digest; anyone who can write the bundle can make a 'complete' one
    without running the audit. The engine now also requires the five cache files to be byte for byte its own completed
    all-pairs cache and the hashes its draft pinned at preparation. The forged bundle is refused before any referee call,
    and no forecast is committed. Editing the draft's pins to the forged matrix as well does not help."""
    directory = clone(built, tmp_path / "w")
    bundle = Path(engine_of(directory, "develop").prepare_referee_cache()["bundle"])
    pinned = load(bundle / "draft_declaration.json")["method_artifacts"][E.REFEREE_CACHE_DESCRIPTOR["distance_path"]]
    forged = forge_receipt(directory, bundle)
    assert forged != pinned
    with pytest.raises(E.EngineError, match="not byte for byte the engine's own completed all-pairs cache"):
        engine_of(directory, "validate").referee_cache_state()
    with pytest.raises(E.EngineError, match="not byte for byte the engine's own completed all-pairs cache"):
        engine_of(directory, "validate").predict()
    assert not (directory / "engine/validate/PREDICTION_COMMITMENT.json").exists()
    draft = load(bundle / "draft_declaration.json")
    draft["method_artifacts"][E.REFEREE_CACHE_DESCRIPTOR["distance_path"]] = forged
    E.write_json(bundle / "draft_declaration.json", draft, pretty=True)
    with pytest.raises(E.EngineError, match=r"\['distance_path'\] are not\s+byte for byte|not byte for byte"):
        engine_of(directory, "validate").predict()
    assert not (directory / "engine/validate/PREDICTION_COMMITMENT.json").exists()


def test_round2_a_referee_copy_through_the_audited_cache_that_differs_from_the_replica_is_refused(built, tmp_path, monkeypatch):
    """Round-2 finding (major), second half: through the audited cache the referee and the engine's replica read the same
    bitwise matrix, so a difference means tampering. It used to be overwritten with the referee's copy and only counted;
    now it is an error and nothing is committed. Without the tampering the same steps commit with 0 differences."""
    directory = clone(built, tmp_path / "w")
    bundle = Path(engine_of(directory, "develop").prepare_referee_cache()["bundle"])
    assert audit_bundle(bundle)["status"] == "complete"
    honest = clone(directory, tmp_path / "honest")
    commit = engine_of(honest, "validate").predict()
    official = load(honest / "engine/validate/prediction_manifest.json")["official_copies"]
    assert commit["referee_training_cache_receipt_sha256"] and official["api_verified"] and official["replica_differences"] == 0
    Ref = E.referee()
    real = Ref.copying_baselines
    probe = engine_of(directory, "validate")

    def tampered(donors, queries, **kw):
        result = real(donors, queries, **kw)
        if kw.get("training_distance_cache"):                  # the referee's answer through the cache, changed for one query
            q = queries[0]["id"]
            other = next(d["id"] for d in donors if d["id"] != result["donor_proof"]["medoid"][q])
            result["donor_proof"]["medoid"][q] = other
            result["predictions"]["medoid"][q] = E.as_tree(probe.donor(other)["arbor"])
        return result

    monkeypatch.setattr(Ref, "copying_baselines", tampered)
    with pytest.raises(E.EngineError, match="differs from the engine's own copy"):
        engine_of(directory, "validate").predict()
    assert not (directory / "engine/validate/PREDICTION_COMMITMENT.json").exists()


def run_probe(code, env):
    out = subprocess.run([sys.executable, "-I", "-c", code], cwd=ROOT, env=env, capture_output=True, text=True, timeout=600)
    return out.stdout + out.stderr


def test_round2_the_thread_guard_checks_the_values_at_import_and_the_real_thread_pools(tmp_path):
    """Round-2 finding: importing the engine loads analysis/brm_real.py, which writes '1' into the four variables after
    numpy is loaded, so os.environ alone proved nothing and a process with 8 BLAS threads passed. The guard now checks
    the values when the engine started and the thread pools themselves."""
    bare = {k: v for k, v in os.environ.items() if k not in THREADS}
    head = f"import sys, os; sys.path.insert(0, {str(ROOT / 'analysis')!r})\n"
    tail = ("print('ENV', [os.environ.get(v) for v in E.THREAD_VARIABLES])\n"
            "try:\n    E.referee(); print('ACCEPTED')\nexcept E.EngineError as e:\n    print('REFUSED', e)\n")
    numpy_first = run_probe(head + "import numpy\nimport axon_engine as E\n" + tail, bare)
    assert "ENV ['1', '1', '1', '1']" in numpy_first, numpy_first       # the variables look right afterwards ...
    assert "REFUSED" in numpy_first and "ACCEPTED" not in numpy_first, numpy_first      # ... and the referee is still refused
    engine_first = run_probe(head + "import axon_engine as E\n" + tail, bare)
    assert "REFUSED" in engine_first and "ACCEPTED" not in engine_first, engine_first
    late = run_probe(head + "import numpy\nos.environ.update(" + repr(THREADS) + ")\nimport axon_engine as E\n" + tail, bare)
    assert "REFUSED" in late and "thread pools" in late, late           # set after numpy: the pools still run many threads
    proper = run_probe(head + "os.environ.update(" + repr(THREADS) + ")\nimport numpy\nimport axon_engine as E\n" + tail, bare)
    assert "ACCEPTED" in proper, proper
    assert all(n == 1 for _, n in (E.thread_pools() or [])), E.thread_pools()          # this test process itself


# ------------------------------------------------------------------------------------------
# the engine's compact adapter of the frozen predict_policy (omits only its two discarded diagnostics)
# ------------------------------------------------------------------------------------------


def _route(sign, soma=(0., 0., 0.), bend=100.):
    points = np.array([[0., 0., 0.], [500. * sign, 0., 0.], [1200. * sign, bend, 0.], [500. * sign, 350., 0.]]) + np.asarray(soma)
    return G.resample(G.Arbor(points, np.array([-1, 0, 1, 1]), np.array([1, 2, 2, 2])), 50.)


def _policy_fixture(arbors):
    """The frozen module's own synthetic fixture shape (tests/test_brm_template.py), built on the engine's copy of it."""
    tags = ["p", "q", "p", "q", "p", "q"]
    donors = []
    for i, arbor in enumerate(arbors):
        neuron, animal, tag = f"T{i}_001", f"T{i}", tags[i % len(tags)]
        donors.append({"neuron": neuron, "animal": animal, "full_id": f"{neuron}_{tag}_HY", "tag": tag,
                       "table_soma": np.zeros(3), "swc_soma": np.zeros(3), "profile": np.full(49, 4. * (tag == "q")), "arbor": arbor})
    queries = [T.Query("E0_001", "E0", "E0_001_p_HY", np.zeros(3), np.zeros(3), "p"),
               T.Query("E1_001", "E1", "E1_001_q_HY", np.zeros(3), np.zeros(3), "q")]
    continuous = [*Geo.R.POSITION, *Geo.R.target_fields()]
    train_ids = [d["full_id"] for d in donors]
    import hashlib
    model = {"core_sha256": Geo.JOINT_SHA, "target_fields": list(Geo.R.target_fields()),
             "field_types": {**{f: "continuous" for f in continuous}, Geo.R.TAG: "categorical"},
             "training_ids_sha256": hashlib.sha256("\n".join(train_ids).encode()).hexdigest(),
             "field_metadata": {f: {"training_mean": 0., "training_scale": 1.} for f in continuous},
             "blocks": [{"name": "block_0", "fields": [*continuous, Geo.R.TAG], "weights": [.5, .5],
                         "means": {f: [0., 0.] if f in Geo.R.POSITION else [0., 4.] for f in continuous},
                         "variances": {f: [1., 1.] if f in Geo.R.POSITION else [.01, .01] for f in continuous},
                         "tables": {Geo.R.TAG: [[.999, .001], [.001, .999]]}}]}
    model["field_metadata"][Geo.R.TAG] = {"categories": ["p", "q"]}
    proof = {"name": "synthetic_train_only", "train_ids": train_ids, "train_animals": [d["animal"] for d in donors],
             "test_ids": [q.full_id for q in queries], "test_animals": [q.animal for q in queries],
             "model_state_sha256": T._digest(model)}
    return donors, Geo.RegionalState(model, proof), queries


def _all_parity(query, donors, regional):
    out = [E.compact_parity(query, donors, regional, policy) for policy in T.POLICIES]
    for r in out:
        assert r["arbor_bitwise_equal"] and r["omitted_keys_exactly_the_two"] and r["records_equal_without_the_two"], r
    return out


def test_compact_adapter_equals_the_frozen_policy_for_all_five_policies_with_valid_attachment():
    donors, regional, queries = _policy_fixture([_route(1), _route(-1), _route(1, bend=180.), _route(-1, bend=40.)])
    for query in queries:
        out = _all_parity(query, donors, regional)
        assert {r["executed_policy"] for r in out} == set(T.POLICIES)              # every attachment was supported
    # reversed input order: the adapter sorts exactly as the frozen function does
    out = _all_parity(queries[1], list(reversed(donors)), regional)
    assert len(out) == 5
    _, info = E.compact_predict_policy(queries[1], donors, regional, "conditional_anchor_500")
    _, frozen = T.predict_policy(queries[1], donors, regional, "conditional_anchor_500")
    assert set(frozen) - set(info) == set(E.COMPACT_OMITTED_INFO_KEYS) and not set(info) - set(frozen)
    assert info["f1_risk"] == frozen["f1_risk"] and info["distance_risk"] == frozen["distance_risk"]


def test_compact_adapter_falls_back_to_the_exact_baseline_like_the_frozen_policy():
    donors, regional, queries = _policy_fixture([_route(1), _route(-1), _route(1, bend=180.)])
    far = T.Query("E1_001", "E1", "E1_001_q_HY", np.array([1000., 0., 0.]), np.array([1000., 0., 0.]), "q")
    out = _all_parity(far, donors, regional)
    assert {r["fallback_reason"] for r in out[1:]} == {"soma_shift_exceeds_radius_over_three"}
    assert {r["executed_policy"] for r in out} == {T.BASELINE}


def test_compact_adapter_keeps_the_frozen_tie_rule():
    same = _route(1)
    donors, regional, queries = _policy_fixture([same, G.Arbor(same.xyz.copy(), same.parent.copy(), same.types.copy()),
                                                 _route(-1), G.Arbor(_route(-1).xyz.copy(), _route(-1).parent.copy(), _route(-1).types.copy())])
    for query in queries:
        for r in _all_parity(query, donors, regional):
            assert r["executed_policy"] in T.POLICIES
    for policy in T.POLICIES:
        _, a = E.compact_predict_policy(queries[0], donors, regional, policy)
        _, b = T.predict_policy(queries[0], donors, regional, policy)
        assert a["selected_donor"] == b["selected_donor"] and a.get("attempted_donor") == b.get("attempted_donor")


def test_compact_adapter_refuses_every_leak_and_bad_input_exactly_like_the_frozen_policy():
    donors, regional, queries = _policy_fixture([_route(1), _route(-1), _route(1, bend=180.)])
    own = dict(donors[0], neuron="E0_002", animal="E0", full_id="E0_002_p_HY")
    cases = [(queries[0], donors + [own], "conditional_anchor_500"),               # the query's animal among the donors
             (queries[0], donors + [dict(donors[0])], "uniform_anchor_500"),       # a duplicated donor
             (queries[0], [], T.BASELINE),                                          # no donors
             (queries[0], donors, "oracle_best_geometry")]                          # outside the fixed family
    for query, pool, policy in cases:
        with pytest.raises(Exception) as frozen:
            T.predict_policy(query, pool, regional, policy)
        with pytest.raises(Exception) as compact:
            E.compact_predict_policy(query, pool, regional, policy)
        assert type(frozen.value) is type(compact.value) and str(frozen.value) == str(compact.value)
    with pytest.raises(TypeError):
        E.compact_predict_policy(queries[0], donors, regional, "baseline", query_profile=np.zeros(49))


def test_compact_adapter_inside_the_frozen_memo_gives_the_same_sixteen_records_and_log(built):
    """On the synthetic world's training queries: the five policies through the pinned cache, and the engine's whole
    forecast (sixteen records and the prediction log) through the adapter and through the unmodified function."""
    engine = engine_of(built, "develop")
    neurons = list(engine.library["ids"])[::max(1, len(engine.library["ids"]) // 4)][:4]
    proof = E.prove_compact_parity(engine, neurons)
    assert proof["queries"] == len(neurons) and proof["all_identical"], proof
    assert all(r["methods"] == len(E.METHODS) for r in proof["rows"])
    assert E.compact_predict_policy is not T.predict_policy                          # the proof's routing is undone
    with pytest.raises(E.EngineError):
        E.prove_compact_parity(engine_of(built, "validate"), neurons)
