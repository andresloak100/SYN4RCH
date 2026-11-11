"""Synthetic tests of the protocol-2 feasibility helpers; no biological file is read."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("next_donor_selector_feasibility_tested",
                                              ROOT / "analysis/next_donor_selector_feasibility.py")
M = importlib.util.module_from_spec(SPEC); sys.modules[SPEC.name] = M; SPEC.loader.exec_module(M)


def animals(n=23):
    return [f"A{i:03d}" for i in range(n)]


def test_folds_are_whole_animal_and_disjoint():
    folds = M.nested_folds(animals())
    seen = []
    for f in folds:
        assert not set(f["fit_animals"]) & set(f["evaluation_animals"])
        seen += f["evaluation_animals"]
        inner_eval = []
        for i in f["inner"]:
            assert set(i["fit_animals"]) | set(i["evaluation_animals"]) == set(f["fit_animals"])
            assert not set(i["fit_animals"]) & set(i["evaluation_animals"])
            assert not set(i["fit_animals"]) & set(f["evaluation_animals"])
            inner_eval += i["evaluation_animals"]
        assert sorted(inner_eval) == sorted(f["fit_animals"])
    assert sorted(seen) == animals()


def test_folds_are_deterministic_and_salted():
    assert M.nested_folds(animals()) == M.nested_folds(list(reversed(animals())))
    assert M.assign_folds(animals(), 5, "outer") != M.assign_folds(animals(), 5, "other")


def test_pair_guard_refuses_evaluation_endpoints():
    animal_of = {"q1": "A", "d1": "B", "e1": "C", "s1": "A"}
    assert M.check_fit_pairs([("q1", "d1")], animal_of, ["A", "B"], ["C"]) == 1
    with pytest.raises(ValueError):
        M.check_fit_pairs([("q1", "e1")], animal_of, ["A", "B"], ["C"])          # donor in evaluation animal
    with pytest.raises(ValueError):
        M.check_fit_pairs([("e1", "d1")], animal_of, ["A", "B"], ["C"])          # query in evaluation animal
    with pytest.raises(ValueError):
        M.check_fit_pairs([("q1", "s1")], animal_of, ["A", "B"], ["C"])          # same animal
    with pytest.raises(ValueError):
        M.check_fit_pairs([("q1", "d1")], animal_of, ["A", "B", "C"], ["C"])     # overlap


def test_candidate_pairs_exclude_own_animal_and_respect_pool():
    rng = np.random.default_rng(0)
    ids = [f"n{i:02d}" for i in range(30)]
    animal_of = {n: f"A{i % 6}" for i, n in enumerate(ids)}
    soma = {n: rng.normal(size=3) * 100 for n in ids}
    pool = ids[:20]
    pairs = M.candidate_pairs(ids[20:], pool, animal_of, soma, 4)
    assert len(pairs) == 10 * 4
    for q, d in pairs:
        assert d in pool and animal_of[q] != animal_of[d]
    # the chosen donors are the 4 nearest other-animal pool donors
    q = ids[25]
    want = sorted((np.linalg.norm(soma[d] - soma[q]), d) for d in pool if animal_of[d] != animal_of[q])[:4]
    assert [d for qq, d in pairs if qq == q] == [d for _, d in want]


def test_fallback_rule_needs_joint_dominance():
    pred = {"s": (500., .30), "m": (600., .25)}
    assert M.fallback_choice("s", "m", pred)[0] == "s"
    assert M.fallback_choice("s", "m", {"s": (500., .20), "m": (600., .25)})[0] == "m"   # better distance only
    assert M.fallback_choice("s", "m", {"s": (700., .40), "m": (600., .25)})[0] == "m"   # better F1 only
    assert M.fallback_choice("s", "m", {"s": (np.nan, .40), "m": (600., .25)})[0] == "m"
    assert M.fallback_choice("m", "m", pred) == ("m", "selector_equals_fallback")


def test_mdd_scales_with_root_n():
    assert M.mdd(100., 100) == pytest.approx(M.mdd(100., 25) / 2)
    assert M.mdd(100., 100) == pytest.approx((1.959964 + 0.841621) * 10., rel=1e-5)


def test_metric_is_translation_invariant_and_offset_bounded():
    """The fact the cache-meaning check rests on, on a synthetic arbor: translating BOTH trees leaves the metric
    unchanged up to rounding; translating only the prediction by delta moves the distance by at most |delta|."""
    _, G, R = M.frozen_modules()
    rng = np.random.default_rng(1)

    def tree(n, seed):
        r = np.random.default_rng(seed)
        xyz = [np.zeros(3)]; parent = [-1]
        for i in range(1, n):
            p = int(r.integers(0, i)); xyz.append(xyz[p] + r.normal(size=3) * 20.); parent.append(p)
        return (np.array(xyz) + r.normal(size=3) * 1000., np.array(parent), np.r_[1, np.full(n - 1, 2)])
    a, b = tree(60, 2), tree(80, 3)
    centred = R.geometry_metrics(R.shifted(a, np.zeros(3)), R.shifted(b, np.zeros(3)))
    moved = R.geometry_metrics(R.shifted(a, b[0][0]), b)
    assert moved["distance_um"] == pytest.approx(centred["distance_um"], abs=1e-9)
    assert moved["f1_200"] == pytest.approx(centred["f1_200"], abs=1e-12)
    delta = rng.normal(size=3); delta *= 4.8 / np.linalg.norm(delta)
    off = R.geometry_metrics(R.shifted(a, b[0][0] + delta), b)
    assert abs(off["distance_um"] - centred["distance_um"]) <= 4.8 + 1e-9
    # the frozen d110 metric and the referee metric agree on the centred pair
    g = G.geometry_metrics(G.Arbor(*a).shifted(np.zeros(3)), G.Arbor(*b).shifted(np.zeros(3)))
    assert g["distance_um"] == centred["distance_um"] and g["f1_200"] == centred["f1_200"]


# ---------------------------------------------------------------------------------------------
# review round 1: final CV, fit-only cache door, fallback in the candidate set, switch controls,
# verdict sentence, known-input loader
# ---------------------------------------------------------------------------------------------

def test_final_folds_and_fit_configurations():
    an = animals(23)
    ff = M.final_folds(an)
    assert sorted(a for f in ff for a in f["evaluation_animals"]) == an
    for f in ff:
        assert not set(f["fit_animals"]) & set(f["evaluation_animals"])
    configs = M.fit_configurations(an)
    assert len(configs) == 5 + 20 + 5 + 1
    for name, fit, forbidden, evaluation in configs:
        assert set(evaluation) <= set(forbidden) and not set(fit) & set(forbidden)
        if "_inner" in name:
            outer = int(name[5])
            assert set(M.nested_folds(an)[outer]["evaluation_animals"]) <= set(forbidden)
    assert configs[-1][0] == "final_all" and configs[-1][1] == an and configs[-1][2] == []


def synthetic_cache(n=24, n_animals=6, seed=4):
    rng = np.random.default_rng(seed)
    ids = [f"n{i:02d}" for i in range(n)]
    animal_of = {x: f"A{i % n_animals}" for i, x in enumerate(ids)}
    X = rng.normal(size=(n, 3)) * 100
    D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1)); F = D / D.max()
    return ids, animal_of, D, F


def test_fit_cache_refuses_every_non_fit_neuron_and_counts():
    ids, animal_of, D, F = synthetic_cache()
    fc = M.FitCache(D, F, ids, animal_of, ["A0", "A1", "A2"])
    fit = [x for x in ids if animal_of[x] in {"A0", "A1", "A2"}]
    other = [x for x in ids if x not in fit]
    assert fc.block(fit, fit).shape == (len(fit), len(fit)) and fc.reads == len(fit) ** 2
    for bad in ([other[0]], ):
        with pytest.raises(ValueError):
            fc.block(bad, fit)
        with pytest.raises(ValueError):
            fc.block(fit, bad, "F")
        with pytest.raises(ValueError):
            fc.pairs([(fit[0], bad[0])])
    assert fc.refused == 3


def test_fallback_is_the_referee_medoid_without_the_query_animal():
    _, _, R = M.frozen_modules()
    ids, animal_of, D, F = synthetic_cache()
    fit_animals = ["A0", "A1", "A2", "A3"]
    pool = [x for x in ids if animal_of[x] in fit_animals]
    fc = M.FitCache(D, F, ids, animal_of, fit_animals)
    fb = M.fallback_donors(fc, pool, pool, R.geometry_medoid)
    for q in pool:
        donors = sorted(x for x in pool if animal_of[x] != animal_of[q])
        ix = [ids.index(x) for x in donors]
        want = donors[R.geometry_medoid(D[np.ix_(ix, ix)], list(range(len(donors))), [animal_of[x] for x in donors])]
        assert fb[q] == want and animal_of[fb[q]] != animal_of[q]


def test_candidate_sets_always_hold_the_fallback_once():
    pairs = [("q1", "a"), ("q1", "b"), ("q2", "c"), ("q2", "f")]
    sets = M.candidate_sets(pairs, {"q1": "f", "q2": "f"})
    assert sets == {"q1": ["a", "b", "f"], "q2": ["c", "f"]}


def test_matched_switch_controls_switch_the_same_queries():
    final = {"q1": "s1", "q2": "f", "q3": "s3"}
    fallback = {"q1": "f", "q2": "f", "q3": "f"}
    nearest = {"q1": "f", "q2": "a", "q3": "b"}
    cands = {"q1": ["a", "s1", "f"], "q2": ["a", "f"], "q3": ["b", "s3", "f"]}
    sn, sr, cov = M.matched_switch_controls(final, fallback, nearest, cands)
    assert cov == {"queries": 3, "switched": 2, "switched_fraction": 2 / 3}
    assert sn == {"q1": "a", "q2": "f", "q3": "b"}           # q1: nearest is the fallback, so the next candidate
    assert sr["q2"] == "f" and sr["q1"] in {"a", "s1"} and sr["q3"] in {"b", "s3"}
    assert M.matched_switch_controls(final, fallback, nearest, cands)[1] == sr     # seeded


def test_verdict_sentence():
    ok = {"medoid": (-50., .02), "nearest": (-200., .05)}
    tight = {"medoid": ((-90., -10.), (.01, .03)), "nearest": ((-300., -100.), (.02, .08))}
    wide = {"medoid": ((-150., 40.), (-.01, .05)), "nearest": ((-300., -100.), (.02, .08))}
    assert M.verdict_sentence(ok, tight) == "gate met and resolved against every control"
    assert M.verdict_sentence(ok, wide) == "gate met, not resolved"
    assert M.verdict_sentence({**ok, "kernel_all": (5., .03)}, {**tight, "kernel_all": ((0, 9), (.01, .05))}) == "gate not met"


def write_world(tmp_path):
    man = tmp_path / "manifest.tsv"; loc = tmp_path / "local.tsv"; soma = tmp_path / "soma.csv"
    regional = tmp_path / "all_neurons_projection_length_fig1.csv"
    regional.write_text("neuron,REGION\nA1_1,999\n")
    man.write_text("neuron\tanimal\tsource_label\tsoma_label\tpart\n"
                   "A1_1\tA1\tOrexin\tLHA\ttrain\nB2_1\tB2\tAgrp\tARH\ttest\n")
    loc.write_text("neuron\tanimal\tsource_label\tsoma_label\tpart\themisphere\treconstruction_type\n"
                   "A1_1\tA1\tOrexin\tLHA\ttrain\tLeft\tAxon_and_dendrite\nB2_1\tB2\tAgrp\tARH\ttest\tRight\tAxon_only\n")
    soma.write_text("file,type,neurons,name,geno,location,AP,DV,ML\n"
                    "A1_1_x,SECRET_TYPE,A1_1,nm,Wt(orexin),LHA,1.0,2.0,3.0\n"
                    "B2_1_x,SECRET_TYPE,B2_1,nm,Agrp,ARH,4.0,5.0,6.0\n")
    return man, loc, soma, regional


def test_known_inputs_never_open_the_regional_table_and_return_known_fields(tmp_path, monkeypatch):
    import builtins, io
    man, loc, soma, regional = write_world(tmp_path)
    opened = []
    real_open = builtins.open
    def spy(file, *a, **k):
        opened.append(str(file)); return real_open(file, *a, **k)
    monkeypatch.setattr(builtins, "open", spy)
    monkeypatch.setattr(io, "open", spy)
    train = M.known_query_inputs(("train",), man, loc, soma, pins=False)
    test = M.known_query_inputs(("test",), man, loc, soma, pins=False)
    assert str(regional) not in opened and not any("projection_length" in p for p in opened)
    assert list(train) == ["A1_1"] and list(test) == ["B2_1"]
    for row in (train["A1_1"], test["B2_1"]):
        assert tuple(row) == M.KNOWN_FIELDS
        assert "SECRET_TYPE" not in repr(row)
        assert tuple(M.referee_row(row)) == M.REFEREE_KNOWN_FIELDS
    assert train["A1_1"]["peptide"] == "Orexin" and test["B2_1"]["reconstruction_type"] == "Axon_only"
    np.testing.assert_array_equal(test["B2_1"]["soma_um"], [4., 5., 6.])
    with pytest.raises(ValueError):
        M.known_query_inputs(("train",), man, loc, regional, pins=False)      # a regional file is refused by name


def test_known_inputs_refuse_unmirrored_soma_and_label_mismatch(tmp_path):
    man, loc, soma, _ = write_world(tmp_path)
    soma.write_text(soma.read_text().replace("4.0,5.0,6.0", "4.0,5.0,6000.0"))
    with pytest.raises(ValueError):
        M.known_query_inputs(("test",), man, loc, soma, pins=False)
    assert M.known_query_inputs(("train",), man, loc, soma, pins=False)       # the bad row is skipped unparsed
    loc.write_text(loc.read_text().replace("A1_1\tA1\tOrexin", "A1_1\tA1\tAgrp"))
    with pytest.raises(ValueError):
        M.known_query_inputs(("train",), man, loc, soma, pins=False)
