"""Typed graph safeguards, source fidelity and reuse of the existing corpus; no network."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("circuit_evidence_layer", ROOT / "analysis/circuit_evidence_layer.py")
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)
MOUSE, FLY = "mus_musculus", "drosophila_melanogaster"


def node(key, kind, species=MOUSE, **kw):
    return {"id": key, "kind": kind, "label": key, "species": species, **kw}


def edge(key, source, target, kind, species=MOUSE, **kw):
    return {"id": key, "source": source, "target": target, "kind": kind, "species": species,
            "study_id": "mouse_study" if species == MOUSE else "fly_study", "measurement_level": "population",
            "source_locator": "Results, figure 1", "conditions": {"sex": None, "n_animals": None}, **kw}


def source():
    return {"schema": C.SOURCES_SCHEMA,
            "studies": [{"study_id": "mouse_study", "species": [MOUSE], "title": "Mouse source", "url": "https://example.test/mouse"},
                        {"study_id": "fly_study", "species": [FLY], "title": "Fly source", "url": "https://example.test/fly"}],
            "nodes": [node("mouse_population", "population", region_id="dmh", marker_genes=["Ppp1r17"]),
                      node("mouse_region", "region", region_id="dmh"), node("intervention", "intervention"),
                      node("outcome", "outcome", outcome_domain="physiology"), node("gene", "gene"),
                      node("fly_population", "population", FLY, marker_genes=["shared_symbol"])], "edges": []}


ATLAS = {"dmh": {"region_id": "dmh", "status": "resolved", "allen_id": 830}}


def test_intervention_outcomes_and_controls_do_not_become_anatomy():
    data = source()
    data["edges"] = [edge("target", "intervention", "mouse_population", "intervention_target"),
                     edge("effect", "intervention", "outcome", "intervention_outcome", direction="increase",
                          perturbation_class="molecular", perturbation_polarity="increase", controls=["vehicle"])]
    graph = C.build_graph(data, atlas=ATLAS)
    assert all(e["layer"] == "functional" and not e["is_anatomical_edge"] for e in graph["edges"])
    assert graph["coverage"]["edges_by_kind"] == {"intervention_outcome": 1, "intervention_target": 1}
    effect = next(e for e in graph["edges"] if e["id"] == "effect")
    assert effect["perturbation_class"] == "molecular" and effect["perturbation_polarity"] == "increase"
    assert effect["controls"] == ["vehicle"]


def test_regions_can_be_placed_but_no_individual_position_is_invented():
    graph = C.build_graph(source(), atlas=ATLAS)
    pop = next(n for n in graph["nodes"] if n["id"] == "mouse_population")
    assert pop["atlas_placement"] == {"level": "region_only", "region_id": "dmh", "allen_id": 830, "individual_cell_position": False}
    assert "coordinates" not in pop and graph["cell_coordinates_inferred"] == 0
    fly = next(n for n in graph["nodes"] if n["id"] == "fly_population")
    assert fly["atlas_placement"] is None


def test_marker_label_is_not_a_gene_assay_or_join():
    graph = C.build_graph(source(), atlas=ATLAS)
    assert graph["edges"] == [] and graph["automatic_marker_or_proximity_joins"] == 0
    pop = next(n for n in graph["nodes"] if n["id"] == "mouse_population")
    assert pop["gene_measurement_scope"] == "identity_annotation_only"
    data = source()
    data["edges"] = [edge("fake_gene", "mouse_population", "gene", "measured_gene_expression")]
    with pytest.raises(C.ValidationError, match="same-cell"):
        C.build_graph(data)


def test_actual_same_cell_gene_assay_stays_separate_from_marker_identity():
    data = source()
    data["nodes"].append(node("cell", "individual_cell"))
    data["edges"] = [edge("measurement", "cell", "gene", "measured_gene_expression", measurement_level="individual_cell",
                          paired_measurement={"cell_id": "cell", "endogenous": True, "assay_id": "assay1", "measurement_id": "measured1"})]
    graph = C.build_graph(data)
    assert graph["edges"][0]["layer"] == "molecular" and not graph["edges"][0]["is_anatomical_edge"]


def test_fly_edge_cannot_be_imported_as_mouse_anatomy():
    data = source()
    data["edges"] = [edge("wrong", "fly_population", "mouse_region", "measured_axon_projection")]
    with pytest.raises(C.ValidationError, match="cross-species"):
        C.build_graph(data)


@pytest.mark.parametrize("basis", sorted(C.FORBIDDEN_BASES))
def test_gene_or_proximity_population_join_is_rejected(basis):
    data = source()
    data["nodes"].append(node("second_population", "population", marker_genes=["Ppp1r17"]))
    data["edges"] = [edge("join", "mouse_population", "second_population", "population_mapping", mapping_basis=basis)]
    with pytest.raises(C.ValidationError, match="shared genes/proximity"):
        C.build_graph(data)


def test_population_mapping_is_explicit_and_never_cell_equivalence():
    data = source()
    data["nodes"].append(node("second_population", "population"))
    data["edges"] = [edge("join", "mouse_population", "second_population", "population_mapping",
                          mapping_basis="explicit_population_correspondence")]
    graph = C.build_graph(data)
    assert graph["edges"][0]["layer"] == "mapping" and not graph["edges"][0]["identifies_synapse"]
    data["nodes"][-1]["kind"] = "individual_cell"
    with pytest.raises(C.ValidationError, match="population mapping"):
        C.build_graph(data)


def test_projection_is_not_an_identified_synapse_or_cell_pair():
    data = source()
    data["edges"] = [edge("projection", "mouse_population", "mouse_region", "measured_axon_projection", evidence_class="anatomical")]
    graph = C.build_graph(data)
    assert graph["edges"][0]["is_anatomical_edge"] and not graph["edges"][0]["identifies_synapse"]
    data["edges"][0]["kind"] = "identified_synapse"
    with pytest.raises(C.ValidationError, match="synapse needs"):
        C.build_graph(data)


def test_identified_synapse_requires_measured_same_species_cell_pair():
    data = source()
    data["nodes"] += [node("pre", "individual_cell"), node("post", "individual_cell")]
    data["edges"] = [edge("synapse", "pre", "post", "identified_synapse", measurement_level="synapse",
                          evidence_class="anatomical", synapse_identification={"method": "electron_microscopy",
                          "presynaptic_cell_id": "pre", "postsynaptic_cell_id": "post", "n_synapses": 3})]
    graph = C.build_graph(data)
    assert graph["edges"][0]["identifies_synapse"]
    data["edges"][0]["synapse_identification"]["method"] = "shared_marker"
    with pytest.raises(C.ValidationError, match="synapse needs"):
        C.build_graph(data)


def test_no_transmitter_name_to_synaptic_sign_transfer():
    data = source()
    data["nodes"].append(node("fly_target", "population", FLY, neurotransmitter="glutamate"))
    data["edges"] = [edge("fly_model", "fly_population", "fly_target", "inferred_connection", FLY,
                          synaptic_sign={"species": FLY, "basis": "source_specific_model", "study_id": "fly_study",
                                         "source_locator": "Fly model assumptions", "value": "inhibitory"})]
    graph = C.build_graph(data)
    assert graph["edges"][0]["synaptic_sign"]["value"] == "inhibitory"
    data["edges"] = [edge("mouse", "mouse_population", "mouse_region", "inferred_connection",
                          synaptic_sign=graph["edges"][0]["synaptic_sign"])]
    with pytest.raises(C.ValidationError, match="sign requires"):
        C.build_graph(data)
    data["edges"][0]["synaptic_sign"]["species"] = MOUSE
    with pytest.raises(C.ValidationError, match="sign requires"):
        C.build_graph(data)                 # relabeling the rule cannot change its fly-only source study
    data["edges"] = [edge("mouse", "mouse_population", "mouse_region", "inferred_connection")]
    assert "synaptic_sign" not in C.build_graph(data)["edges"][0]


def test_missing_sex_n_and_unreviewed_status_are_visible_not_filled():
    data = source()
    data["edges"] = [edge("e", "intervention", "outcome", "intervention_outcome")]
    graph = C.build_graph(data)
    assert graph["edges"][0]["conditions"] == {"sex": None, "n_animals": None}
    assert graph["edges"][0]["review_status"] == "unreviewed"
    assert graph["coverage"]["missing_condition_fields"] == {"sex": 1, "n_animals": 1}


def test_numeric_confidence_and_invalid_pooled_animal_n_are_rejected():
    data = source()
    data["edges"] = [edge("e", "intervention", "outcome", "intervention_outcome", confidence=.95)]
    with pytest.raises(C.ValidationError, match="numeric confidence"):
        C.build_graph(data)
    data["edges"][0].pop("confidence")
    data["edges"][0]["conditions"]["n_animals"] = 2.5
    with pytest.raises(C.ValidationError, match="animal N"):
        C.build_graph(data)


def test_pvh_sampling_preserves_inferred_type_and_pooled_units():
    pvh = {"schema": "synarch_pvh_evidence_graph_v1",
           "unknown_pool_nuclei": 6, "scope": "Original-study nuclei only; sampling counts",
           "not_present": ["individual animals", "other-target negatives"], "provenance": {"run_id": "source_run"},
           "nodes": [{"id": "region:PVH", "kind": "source_region", "label": "PVH"},
                     {"id": "type:GeneA", "kind": "author_inferred_type", "label": "GeneA"},
                     {"id": "pool:1", "kind": "encoded_pool", "label": "pool1", "sex": "female"},
                     {"id": "target:SC", "kind": "collection_target", "label": "SC collection"}],
           "edges": [{"kind": "collected_in", "source": "type:GeneA", "target": "target:SC", "n_nuclei": 4,
                      "fraction_of_target_collection": .4},
                     {"kind": "pool_collected_for", "source": "pool:1", "target": "target:SC", "n_nuclei": 10}]}
    graph = C.build_graph(source(), pvh_graph=pvh)
    assert all(e["kind"] == "retrograde_collection" and e["layer"] == "sampling" for e in graph["edges"])
    assert all(e["conditions"]["n_animals"] is None for e in graph["edges"])
    pop = next(n for n in graph["nodes"] if n["id"] == "pvh_sampling:type:GeneA")
    assert pop["identity_status"] == "author_inferred_type" and pop["marker_genes"] == []
    assert not any(e["kind"] == "population_mapping" for e in graph["edges"])
    assert graph["pvh_sampling_metadata"] == {key: pvh[key] for key in
                                              ("schema", "scope", "unknown_pool_nuclei", "not_present", "provenance")}


def test_existing_corpus_is_reused_without_upgrading_or_defaulting_species():
    records = [json.loads(x) for x in (ROOT / "data/curated/records.jsonl").read_text().splitlines() if x.strip()]
    populations = C.read(ROOT / "data/vocab/populations.json")["entries"]
    atlas = {r["region_id"]: r for r in C.read(ROOT / "viewer/data/allen_mapping.json")["regions"]}
    empty = {"schema": C.SOURCES_SCHEMA, "nodes": [], "edges": [], "studies": []}
    graph = C.build_graph(empty, records, populations, atlas)
    assert graph["legacy_curated_records_imported"] == len(records)
    assert len(graph["edges"]) == 2 * len(records)
    assert not any(e["is_anatomical_edge"] for e in graph["edges"])
    assert all(e["review_status"] == "unreviewed" for e in graph["edges"])
    unknown = [e for e in graph["edges"] if e["species"] is None]
    assert len(unknown) == 2 * sum(r.get("species") in {None, "other"} for r in records)
    assert any(e["source_species_value"] == "other" for e in unknown)
    assert len(graph["studies"]) == len({r["study_id"] for r in records})
    nodes = {n["id"]: n for n in graph["nodes"]}
    unresolved_records = {r["record_id"] for r in records if r.get("population_id") not in populations}
    for e in graph["edges"]:
        if e["kind"] == "intervention_target" and e["record_id"] in unresolved_records:
            assert nodes[e["target"]]["kind"] == "region"


def test_new_curated_nodes_cannot_silently_reuse_pvh_or_legacy_identity():
    data = source()
    data["nodes"][0]["id"] = "pvh_sampling:type:Ppp1r17"
    with pytest.raises(C.ValidationError, match="explicit population mapping"):
        C.build_graph(data)


def test_shared_outcome_label_is_canonical_and_each_raw_timepoint_stays_on_its_edge():
    records = []
    for rid, sid, raw in (("r1", "study1", "Food intake increased at 1 hour"),
                          ("r2", "study2", "Food intake decreased at 24 hours")):
        records.append({"record_id": rid, "study_id": sid, "species": MOUSE, "region_id": "dmh",
                        "population_id": None, "population_raw": "DMH neurons", "outcome_id": "food_intake",
                        "outcome_raw": raw, "source_locator": "Figure 1", "direction": "increase",
                        "method": "chemogenetic", "valence": "activation"})
    graph = C.build_graph({"schema": C.SOURCES_SCHEMA, "studies": [], "nodes": [], "edges": []}, records,
                          outcomes={"food_intake": {"label": "Food intake / feeding"}})
    shared = [n for n in graph["nodes"] if n["kind"] == "outcome"]
    assert len(shared) == 1 and shared[0]["label"] == "Food intake / feeding"
    observed = {e["record_id"]: e["source_outcome_raw"] for e in graph["edges"] if e["kind"] == "intervention_outcome"}
    assert observed == {"r1": records[0]["outcome_raw"], "r2": records[1]["outcome_raw"]}
    fallback = C.build_graph({"schema": C.SOURCES_SCHEMA, "studies": [], "nodes": [], "edges": []}, records)
    assert next(n for n in fallback["nodes"] if n["kind"] == "outcome")["label"] == "food intake"


def test_distinct_molecular_targets_methods_and_reported_statistics_keep_source_context():
    records = []
    for rid, target, p_value, significant in (("gene_a", "GeneA", .01, True), ("gene_b", "GeneB", .2, False)):
        records.append({"record_id": rid, "study_id": "study1", "species": MOUSE, "population_id": "dmh_population",
                        "population_raw": "DMH population", "region_id": "dmh", "outcome_id": "food_intake",
                        "outcome_raw": "Reported food intake", "source_locator": "Results " + rid,
                        "method": "conditional_knockout", "actuator": None, "perturbation_class": "molecular",
                        "perturbation_polarity": "decrease", "target_raw": target, "target_id": target.lower(),
                        "projection_target_raw": "PVH" if rid == "gene_a" else None,
                        "projection_target_id": "pvh" if rid == "gene_a" else None,
                        "p_value": p_value, "statistical_test": "source_reported_test", "significant": significant})
    graph = C.build_graph({"schema": C.SOURCES_SCHEMA, "studies": [], "nodes": [], "edges": []}, records,
                          populations={"dmh_population": {"label": "DMH population", "nucleus": "dmh"}})
    targets = {e["record_id"]: e for e in graph["edges"] if e["kind"] == "intervention_target"}
    assert targets["gene_a"]["target"] == targets["gene_b"]["target"]
    assert targets["gene_a"]["source"] != targets["gene_b"]["source"]
    for record in records:
        result = targets[record["record_id"]]
        assert result["method"] == "conditional_knockout" and result["actuator"] is None
        assert result["molecular_target"] == {"raw": record["target_raw"], "id": record["target_id"]}
        assert result["reported_statistics"] == {key: record[key] for key in ("p_value", "statistical_test", "significant")}
        assert result["projection_target_raw"] == record["projection_target_raw"]
        assert result["projection_target_id"] == record["projection_target_id"]


def test_new_curated_statistics_align_without_inventing_p_or_significance():
    data = source()
    data["edges"] = [edge("e", "intervention", "outcome", "intervention_outcome", method="chemogenetic",
                          statistical_result={"test": "log_rank", "p_value": .0406})]
    result = C.build_graph(data)["edges"][0]
    assert result["method"] == "chemogenetic" and result["actuator"] is None
    assert result["reported_statistics"] == {"p_value": .0406, "statistical_test": "log_rank", "significant": None}
    assert result["statistical_result"] == data["edges"][0]["statistical_result"]
    data["edges"][0]["statistical_result"] = {"test": "source_test", "reported_p_threshold": .05}
    result = C.build_graph(data)["edges"][0]
    assert result["reported_statistics"]["p_value"] is None
    assert result["reported_statistics"]["significant"] is None
    assert result["reported_statistics"]["reported_p_threshold"] == .05


def test_reported_cohort_dispersion_and_readout_context_are_preserved_exactly():
    record = {"record_id": "context1", "study_id": "study1", "species": MOUSE,
              "population_id": None, "population_raw": "DMH neurons", "region_id": "dmh",
              "outcome_id": "food_intake", "outcome_raw": "Source-specific food intake readout",
              "source_locator": "Figure 2", "method": "conditional_knockout", "valence": "unclear",
              "background_genotype": "reported background genotype", "genetic_background": "reported genetic background",
              "strain": "reported strain", "challenge_detail": "reported metabolic challenge context",
              "dispersion_type": "sem", "dispersion_value": 1.25, "assay": "reported readout assay",
              "assay_units": "reported units", "tissue_raw": "reported tissue", "tissue_id": None}
    graph = C.build_graph({"schema": C.SOURCES_SCHEMA, "studies": [], "nodes": [], "edges": []}, [record])
    for e in graph["edges"]:
        for field in ("background_genotype", "genetic_background", "strain", "challenge_detail"):
            assert e["conditions"][field] == record[field]
    outcome = next(e for e in graph["edges"] if e["kind"] == "intervention_outcome")
    assert outcome["effect"]["dispersion_type"] == "sem" and outcome["effect"]["dispersion_value"] == 1.25
    assert outcome["readout_context"] == {k: record[k] for k in ("assay", "assay_units", "tissue_raw", "tissue_id")}
    assert outcome["conditions"]["sex"] is None and outcome["conditions"]["n_animals"] is None


def test_local_cli_build_validate_and_source_hashes(tmp_path):
    curated, records, atlas, populations = [tmp_path / x for x in ("curated.json", "records.jsonl", "atlas.json", "populations.json")]
    curated.write_text(json.dumps(source()))
    records.write_text("")
    atlas.write_text(json.dumps({"regions": list(ATLAS.values())}))
    populations.write_text(json.dumps({"entries": {}}))
    graph = C.build(curated, records, None, atlas, populations, tmp_path / "output")
    assert (tmp_path / "output/coverage.json").exists()
    assert graph["sources"][0]["sha256"] == C.sha(curated)
    result = subprocess.run([sys.executable, "-I", str(ROOT / "analysis/circuit_evidence_layer.py"), "validate",
                             str(tmp_path / "output/graph.json")], capture_output=True, text=True)
    assert result.returncode == 0 and json.loads(result.stdout)["valid"]
