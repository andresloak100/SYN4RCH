"""Build/validate a local, typed anatomical and functional evidence layer.

Reuses SYN4RCH's Tier-A conditions, provenance and perturbation axes. Literature
observations remain observations, including their unreviewed status; controls do
not become anatomical edges. No network, model fitting, cell-position inference,
species transfer, construct details or experimental parameters are implemented.

CLI: python3 -I analysis/circuit_evidence_layer.py build --curated-sources PATH
     python3 -I analysis/circuit_evidence_layer.py validate PATH_TO_GRAPH
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from syn4rch.atlas.signaling import axis_class_of, axis_of  # noqa: E402
from syn4rch.schema.enums import Direction, PerturbationClass, PerturbationPolarity, Sex  # noqa: E402

SCHEMA = "synarch_circuit_evidence_graph_v1"
SOURCES_SCHEMA = "synarch_circuit_sources_v1"
OUT = ROOT / "data/processed/circuit_evidence_v1"
NODE_KINDS = {"individual_cell", "population", "region", "gene", "intervention", "outcome", "organ",
              "collection_group", "animal_pool"}
EDGE_LAYERS = {
    "measured_axon_projection": "anatomical", "identified_synapse": "anatomical",
    "inferred_connection": "inference", "retrograde_label": "sampling", "retrograde_collection": "sampling",
    "intervention_target": "functional", "intervention_outcome": "functional", "functional_association": "functional",
    "marker_identity": "identity", "measured_gene_expression": "molecular", "population_mapping": "mapping",
}
LEVELS = {"individual_cell", "synapse", "population", "region", "animal", "pooled_animals", "nucleus",
          "study", "population_annotation", "inferred", "organ"}
BIOLOGICAL = {"individual_cell", "population", "region", "organ"}
FORBIDDEN_BASES = {"shared_gene", "shared_gene_name", "shared_marker", "marker_only", "proximity", "spatial_proximity"}
LIMITATIONS = [
    "Evidence integration, not aging reversal, clinical efficacy or a complete connectome.",
    "Intervention outcomes and their controls are not axon projections or identified synapses.",
    "Anatomical projections, retrograde collection, identified synapses and inferred connections remain different edge types.",
    "Identity-marker names are not actual paired endogenous gene measurements in an individual source cell.",
    "Species and study populations remain explicit; shared genes and spatial proximity cannot identify cells or import fly edges into mice.",
    "Transmitter labels do not automatically determine synaptic sign; sign evidence is species/source specific.",
    "Atlas placement is regional only. Example Jiao traces are not identified as new study populations.",
    "Null means unknown/not curated. Not_reported is preserved when a source was checked. Missing edges are evidence gaps, not biological negatives.",
]


class ValidationError(ValueError):
    pass


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def enum_values(cls):
    return {x.value for x in cls}


def identifier(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValidationError(label + " must be a nonempty identifier")
    return value


def reject_numeric_confidence(value):
    if isinstance(value, float) and not math.isfinite(value):
        raise ValidationError("non-finite numeric evidence is not valid JSON data")
    if isinstance(value, dict):
        for key, item in value.items():
            if key.lower() in {"confidence", "confidence_score", "confidence_probability", "edge_probability"} and isinstance(item, (int, float)):
                raise ValidationError("numeric confidence/probability must not be invented")
            reject_numeric_confidence(item)
    elif isinstance(value, list):
        for item in value:
            reject_numeric_confidence(item)


def validate_graph(graph):
    if graph.get("schema") != SCHEMA:
        raise ValidationError("unsupported circuit graph schema")
    reject_numeric_confidence(graph)
    nodes, studies = {}, {}
    for node in graph.get("nodes", []):
        key = identifier(node.get("id"), "node id")
        if key in nodes or node.get("kind") not in NODE_KINDS:
            raise ValidationError("duplicate node or unknown node kind: " + key)
        identifier(node.get("label"), "node label")
        species = node.get("species")
        if species is not None and (not isinstance(species, str) or species in {"other", "mouse", "fly"}):
            raise ValidationError("species must be explicit or null")
        placement = node.get("atlas_placement")
        if placement and (species != "mus_musculus" or placement.get("level") != "region_only"
                          or placement.get("individual_cell_position") is not False
                          or placement.get("region_id") != node.get("region_id")):
            raise ValidationError("mouse atlas placement cannot become another species or an individual cell location")
        if node.get("marker_genes") and node.get("gene_measurement_scope") != "identity_annotation_only":
            raise ValidationError("marker names must be labeled as identity annotations")
        nodes[key] = node
    for study in graph.get("studies", []):
        key = identifier(study.get("study_id"), "study id")
        if key in studies or not isinstance(study.get("species", []), list):
            raise ValidationError("duplicate study or non-list species")
        studies[key] = study
    seen = set()
    for edge in graph.get("edges", []):
        key = identifier(edge.get("id"), "edge id")
        if key in seen:
            raise ValidationError("duplicate edge: " + key)
        seen.add(key)
        kind = edge.get("kind")
        if kind not in EDGE_LAYERS or edge.get("layer") != EDGE_LAYERS[kind]:
            raise ValidationError("unknown edge type or incorrect evidence layer")
        if edge.get("is_anatomical_edge") != (kind in {"measured_axon_projection", "identified_synapse"}):
            raise ValidationError("functional/control evidence cannot be relabeled as anatomical")
        if edge.get("source") not in nodes or edge.get("target") not in nodes:
            raise ValidationError("edge endpoint is missing")
        source, target = nodes[edge["source"]], nodes[edge["target"]]
        species = edge.get("species")
        if source.get("species") != species or target.get("species") != species:
            raise ValidationError("cross-species edge import/identity transfer is forbidden")
        study = studies.get(edge.get("study_id"))
        if study is None or (species is not None and species not in study.get("species", [])):
            raise ValidationError("edge must retain a matching source-study species")
        if edge.get("measurement_level") not in LEVELS:
            raise ValidationError("measurement level must be explicit")
        if not edge.get("source_locator"):
            raise ValidationError("every statement needs a source locator")
        conditions = edge.get("conditions")
        if not isinstance(conditions, dict) or conditions.get("sex") not in enum_values(Sex) | {None}:
            raise ValidationError("source sex/unknown must be preserved")
        n = conditions.get("n_animals")
        if n is not None and (type(n) is not int or n < 1):
            raise ValidationError("known animal N must be a positive integer; pooled nuclei are not animals")
        if edge.get("direction") not in enum_values(Direction) | {None}:
            raise ValidationError("invalid outcome direction")
        if edge.get("perturbation_class") not in enum_values(PerturbationClass) | {None}:
            raise ValidationError("invalid perturbation class")
        if edge.get("perturbation_polarity") not in enum_values(PerturbationPolarity) | {None}:
            raise ValidationError("invalid perturbation polarity")
        for field in ("join_basis", "mapping_basis", "identity_basis"):
            if edge.get(field) in FORBIDDEN_BASES:
                raise ValidationError("shared genes/proximity cannot justify a biological join")
        sk, tk = source["kind"], target["kind"]
        if kind in {"measured_axon_projection", "identified_synapse"} and edge.get("evidence_class") not in {None, "anatomical"}:
            raise ValidationError("causal/functional evidence is not measured anatomy")
        if kind == "measured_axon_projection" and (sk not in BIOLOGICAL or tk not in {"region", "organ"}):
            raise ValidationError("projection evidence cannot identify an individual postsynaptic cell")
        if kind == "identified_synapse":
            proof = edge.get("synapse_identification") or {}
            if (sk != "individual_cell" or tk != "individual_cell" or edge["measurement_level"] != "synapse"
                    or proof.get("method") not in {"electron_microscopy", "correlated_light_electron_microscopy"}
                    or proof.get("presynaptic_cell_id") != edge["source"] or proof.get("postsynaptic_cell_id") != edge["target"]
                    or type(proof.get("n_synapses")) is not int or proof["n_synapses"] < 1):
                raise ValidationError("identified synapse needs explicit measured EM cell-pair/count evidence")
        if kind == "intervention_target" and (sk != "intervention" or tk not in BIOLOGICAL | {"gene"}):
            raise ValidationError("intervention target is not a connectome edge")
        if kind == "intervention_outcome" and (sk != "intervention" or tk != "outcome"):
            raise ValidationError("outcome evidence must connect an intervention to an outcome")
        if kind == "inferred_connection" and (sk not in BIOLOGICAL or tk not in BIOLOGICAL):
            raise ValidationError("inferred connection endpoints must be biological entities")
        if kind in {"retrograde_label", "retrograde_collection"} and (sk not in BIOLOGICAL | {"animal_pool"} or tk not in {"region", "organ", "collection_group"}):
            raise ValidationError("retrograde label/collection must retain its sampling domain")
        if kind == "marker_identity" and not ((sk in {"population", "individual_cell"} and tk == "gene")
                or (sk == "gene" and tk in {"population", "individual_cell"})):
            raise ValidationError("marker identity is an annotation, not a projection or gene assay")
        if kind == "measured_gene_expression":
            proof = edge.get("paired_measurement") or {}
            if (sk != "individual_cell" or tk != "gene" or proof.get("cell_id") != edge["source"]
                    or proof.get("endogenous") is not True or not proof.get("assay_id") or not proof.get("measurement_id")):
                raise ValidationError("actual paired gene measurement requires a measured same-cell endogenous assay")
        if kind == "population_mapping" and (sk != "population" or tk != "population"
                or edge.get("mapping_basis") != "explicit_population_correspondence"):
            raise ValidationError("population mapping requires explicit supported population correspondence; never cell equivalence")
        sign = edge.get("synaptic_sign")
        if sign is not None:
            sign_study = studies.get(sign.get("study_id")) if isinstance(sign, dict) else None
            if (kind not in {"identified_synapse", "inferred_connection", "functional_association"}
                    or sk not in BIOLOGICAL or tk not in BIOLOGICAL
                    or not isinstance(sign, dict) or sign.get("species") != species
                    or sign.get("basis") not in {"measured", "source_specific_model"}
                    or sign_study is None or species not in sign_study.get("species", [])
                    or sign.get("value") not in {"excitatory", "inhibitory", "mixed", "unknown"}
                    or not sign.get("source_locator")):
                raise ValidationError("synaptic sign requires source/species-specific evidence, not a transmitter-name lookup")
    return {"valid": True, "nodes": len(nodes), "edges": len(seen), "studies": len(studies)}


def normalized_node(node, atlas):
    node = deepcopy(node)
    node.setdefault("region_id", None)
    node.setdefault("marker_genes", [])
    node.setdefault("identity_status", "not_reported")
    node.setdefault("study_ids", [])
    node["gene_measurement_scope"] = "identity_annotation_only" if node["marker_genes"] else None
    node["atlas_placement"] = None
    region = atlas.get(node["region_id"])
    if node.get("species") == "mus_musculus" and region and region.get("status") == "resolved":
        node["atlas_placement"] = {"level": "region_only", "region_id": node["region_id"],
                                   "allen_id": region["allen_id"], "individual_cell_position": False}
    return node


def normalized_edge(edge):
    edge = deepcopy(edge)
    edge["layer"] = EDGE_LAYERS.get(edge.get("kind"))
    conditions = {"sex": None, "n_animals": None}
    conditions.update(edge.get("conditions") or {})
    edge["conditions"] = conditions
    for name in ("direction", "effect", "perturbation_class", "perturbation_polarity", "support_quote"):
        edge.setdefault(name, None)
    edge.setdefault("method", None)
    edge.setdefault("actuator", None)
    edge.setdefault("molecular_target", {"raw": edge.get("target_raw"), "id": edge.get("target_id")})
    if "reported_statistics" not in edge:
        # Field aliases only. A reported threshold is never substituted for an exact p,
        # and significance is never inferred from a p value.
        statistical_result = edge.get("statistical_result") or {}
        effect = edge.get("effect") if isinstance(edge.get("effect"), dict) else {}

        def reported(names):
            for values in (edge, statistical_result, effect):
                for name in names:
                    if values.get(name) is not None:
                        return values[name]
            return None

        edge["reported_statistics"] = {"p_value": reported(("p_value",)),
                                       "statistical_test": reported(("statistical_test", "test")),
                                       "significant": reported(("significant",))}
        threshold = reported(("reported_p_threshold",))
        if threshold is not None:
            edge["reported_statistics"]["reported_p_threshold"] = threshold
    edge.setdefault("controls", [])
    edge.setdefault("uncertainty", [])
    edge.setdefault("review_status", "unreviewed")
    edge.setdefault("extraction_method", None)
    edge["is_anatomical_edge"] = edge.get("kind") in {"measured_axon_projection", "identified_synapse"}
    edge["identifies_synapse"] = edge.get("kind") == "identified_synapse"
    return edge


def add_node(nodes, node, atlas):
    node = normalized_node(node, atlas)
    if node["id"] in nodes:
        old = nodes[node["id"]]
        for key in ("kind", "species", "region_id"):
            if old.get(key) != node.get(key):
                raise ValidationError("node identity conflict: " + node["id"])
        old["study_ids"] = sorted(set(old["study_ids"]) | set(node["study_ids"]))
    else:
        nodes[node["id"]] = node


def import_tier_a(records, nodes, edges, studies, populations, atlas, outcomes):
    """Reuse the existing observation grain; do not turn a projection-specific intervention into anatomy."""
    for record in records:
        rid, sid = record["record_id"], record["study_id"]
        species = record.get("species")
        species = None if species == "other" else species
        study = studies.setdefault(sid, {"study_id": sid, "title": record.get("title"), "doi": record.get("doi"),
                                        "url": ("https://doi.org/" + record["doi"]) if record.get("doi") else None,
                                        "species": [], "source_kind": "existing_tier_a_records"})
        if species and species not in study["species"]:
            study["species"].append(species)
        namespace = species or "unknown_species"
        pid = record.get("population_id")
        population = populations.get(pid, {})
        resolved = pid in populations
        kind = "population" if resolved else "region"
        region_id = population.get("nucleus") if resolved else record.get("region_id")
        target = f"tier_a:{namespace}:{kind}:{pid if resolved else region_id or rid}"
        label = (population.get("label") or record.get("population_raw")) if resolved else (
            atlas.get(region_id, {}).get("label") or record.get("region_raw") or "Unresolved intervention region")
        markers = (population.get("marker_genes", []) if species == "mus_musculus" else record.get("marker_genes", [])) if resolved else []
        add_node(nodes, {"id": target, "kind": kind, "species": species,
                        "label": label or "Unresolved intervention region", "region_id": region_id,
                        "marker_genes": markers,
                        "identity_status": "vocabulary_population" if resolved else "region_level_source_record",
                        "legacy_population_id": pid, "study_ids": [sid]}, atlas)
        intervention, outcome = "tier_a:intervention:" + rid, f"tier_a:{namespace}:outcome:{record.get('outcome_id') or rid}"
        cls, polarity = axis_class_of(record), axis_of(record)
        # Legacy valence names are translated only via the repository's already stored paired axes.
        if polarity not in enum_values(PerturbationPolarity):
            polarity = {"activation": "increase", "inhibition": "decrease", "ablation": "remove"}.get(polarity, "unclear")
        add_node(nodes, {"id": intervention, "kind": "intervention", "species": species,
                        "label": f"{cls} perturbation ({polarity})", "study_ids": [sid], "record_id": rid,
                        "method": record.get("method"), "actuator": record.get("actuator"),
                        "molecular_target": {"raw": record.get("target_raw"), "id": record.get("target_id")}}, atlas)
        oid = record.get("outcome_id")
        outcome_label = (outcomes.get(oid, {}).get("label") or oid.replace("_", " ")) if oid else (
            record.get("outcome_raw") or "Unresolved outcome")
        add_node(nodes, {"id": outcome, "kind": "outcome", "species": species,
                        "label": outcome_label,
                        "outcome_id": record.get("outcome_id"), "outcome_domain": "unclassified", "study_ids": [sid]}, atlas)
        common = {"study_id": sid, "species": species, "measurement_level": "animal",
                  "source_locator": record.get("source_locator"), "record_id": rid,
                  "source_collection": "data/curated/records.jsonl", "evidence_class": record.get("evidence_class"),
                  "source_species_value": record.get("species"), "source_region_id": record.get("region_id"),
                  "source_population_raw": record.get("population_raw"),
                  "conditions": {key: record.get(key) for key in ("sex", "n_animals", "age_value", "age_unit", "age_group", "diet", "challenge",
                                                                  "background_genotype", "genetic_background", "strain", "challenge_detail")},
                  "conditions_evidence": record.get("condition_evidence", []), "controls": record.get("control_types", []),
                  "perturbation_class": cls, "perturbation_polarity": polarity,
                  "method": record.get("method"), "actuator": record.get("actuator"),
                  "molecular_target": {"raw": record.get("target_raw"), "id": record.get("target_id")},
                  "reported_statistics": {key: record.get(key) for key in ("p_value", "statistical_test", "significant")},
                  "review_status": record.get("review_status"), "extraction_method": record.get("extraction_method"),
                  "uncertainty": ["Existing curated observation; source review status is retained, not upgraded."]}
        edges.append(normalized_edge({**common, "id": "tier_a:target:" + rid, "source": intervention, "target": target,
                                      "kind": "intervention_target", "projection_specific": record.get("projection_specific"),
                                      "projection_target_id": record.get("projection_target_id"),
                                      "projection_target_raw": record.get("projection_target_raw")}))
        effect = {key: record.get(key) for key in ("effect_size_value", "effect_size_type", "effect_size_units", "ci_low", "ci_high",
                                                   "n_treatment", "n_control", "dispersion_type", "dispersion_value")}
        edges.append(normalized_edge({**common, "id": "tier_a:outcome:" + rid, "source": intervention, "target": outcome,
                                      "kind": "intervention_outcome", "direction": record.get("direction"), "effect": effect,
                                      "readout_context": {key: record.get(key) for key in ("assay", "assay_units", "tissue_raw", "tissue_id")},
                                      "source_outcome_raw": record.get("outcome_raw"),
                                      "source_quote_available_in_original_record": bool(record.get("support_quote"))}))


def import_pvh(graph, nodes, edges, studies, atlas):
    if graph.get("schema") != "synarch_pvh_evidence_graph_v1":
        raise ValidationError("unsupported PVH graph; no silent biological reinterpretation")
    sid = "doi:10.1016/j.celrep.2025.116904"
    studies.setdefault(sid, {"study_id": sid, "title": "A spatial and projection-based transcriptomic atlas of paraventricular hypothalamic cell types",
                            "doi": "10.1016/j.celrep.2025.116904", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC13007250/",
                            "species": ["mus_musculus"], "source_kind": "existing_pvh_population_graph"})
    kinds = {"source_region": "region", "collection_target": "collection_group",
             "author_inferred_type": "population", "encoded_pool": "animal_pool"}
    for original in graph["nodes"]:
        kind = kinds.get(original["kind"])
        if kind is None:
            raise ValidationError("unknown PVH node kind")
        node = {"id": "pvh_sampling:" + original["id"], "kind": kind, "label": original["label"],
                "species": "mus_musculus", "region_id": "pvh" if kind in {"region", "population"} else None,
                "study_ids": [sid], "identity_status": "author_inferred_type" if original["kind"] == "author_inferred_type" else "source_defined",
                "marker_genes": [], "original_node": original,
                "is_individual_animal": False if kind == "animal_pool" else None}
        add_node(nodes, node, atlas)
    for i, original in enumerate(graph["edges"]):
        if original["kind"] not in {"collected_in", "pool_collected_for", "retrograde_collection_from"}:
            raise ValidationError("unknown PVH edge kind")
        edges.append(normalized_edge({"id": "pvh_sampling:edge:" + str(i), "kind": "retrograde_collection",
                                      "source": "pvh_sampling:" + original["source"], "target": "pvh_sampling:" + original["target"],
                                      "species": "mus_musculus", "study_id": sid, "measurement_level": "pooled_animals",
                                      "source_locator": "PVH population atlas evidence_graph.json edge " + str(i),
                                      "conditions": {"sex": nodes["pvh_sampling:" + original["source"]].get("original_node", {}).get("sex"), "n_animals": None},
                                      "n_nuclei": original["n_nuclei"], "sampling_fraction": original.get("fraction_of_target_collection"),
                                      "meaning": "Observed collected-nucleus count, not prevalence, axon geometry, synapses or other-target negatives.",
                                      "original_edge_kind": original["kind"], "review_status": "source_metadata",
                                      "uncertainty": ["Author type is inferred; pools are not individual mice; other-target status is unknown."]}))


def coverage_of(graph):
    edges, nodes = graph["edges"], graph["nodes"]
    missing = {field: sum(e.get("conditions", {}).get(field) is None for e in edges) for field in ("sex", "n_animals")}
    gaps = []
    for node in nodes:
        if node["kind"] != "population":
            continue
        outgoing = [e for e in edges if e["source"] == node["id"]]
        incoming = [e for e in edges if e["target"] == node["id"]]
        absent = []
        if not any(e["kind"] == "measured_axon_projection" for e in outgoing):
            absent.append("no_measured_axon_projection_in_layer")
        if not any(e["kind"] == "identified_synapse" for e in outgoing + incoming):
            absent.append("no_identified_synapse_in_layer")
        absent.append("population_identity_is_not_a_same_individual_cell_multimodal_join")
        if node.get("atlas_placement") is None:
            absent.append("no_mouse_regional_atlas_placement")
        gaps.append({"node_id": node["id"], "species": node.get("species"), "gaps": absent})
    not_reported = {field: sum(e.get("conditions", {}).get(field) == "not_reported" for e in edges) for field in ("sex", "n_animals")}
    return {"schema": "synarch_circuit_evidence_coverage_v1", "nodes": len(nodes), "edges": len(edges),
            "studies": len(graph["studies"]), "nodes_by_kind": dict(Counter(n["kind"] for n in nodes)),
            "edges_by_kind": dict(Counter(e["kind"] for e in edges)), "edges_by_layer": dict(Counter(e["layer"] for e in edges)),
            "edges_by_species": dict(Counter(e.get("species") or "unknown" for e in edges)),
            "edges_by_measurement_level": dict(Counter(e["measurement_level"] for e in edges)),
            "edges_by_review_status": dict(Counter(e.get("review_status") or "unknown" for e in edges)),
            "missing_condition_fields": missing, "checked_not_reported_condition_fields": not_reported, "population_gaps": gaps,
            "counts_meaning": "Counts of source observations/entities, not probabilities, confidence or independent animal replication.",
            "no_measured_edge_is_not_a_negative": True}


def build_graph(curated, records=(), populations=None, atlas=None, pvh_graph=None, sources=None, outcomes=None):
    if curated.get("schema") != SOURCES_SCHEMA:
        raise ValidationError("unsupported local curated-source schema")
    nodes, edges, studies = {}, [], {}
    atlas, populations = atlas or {}, populations or {}
    for study in curated.get("studies", []):
        if study["study_id"] in studies:
            raise ValidationError("duplicate curated study")
        studies[study["study_id"]] = deepcopy(study)
    for node in curated.get("nodes", []):
        if node["id"].startswith(("tier_a:", "pvh_sampling:")):
            raise ValidationError("curated nodes cannot reuse imported IDs; use an explicit population mapping")
        if node["id"] in nodes:
            raise ValidationError("duplicate curated node")
        add_node(nodes, node, atlas)
    for edge in curated.get("edges", []):
        edges.append(normalized_edge(edge))
    import_tier_a(records, nodes, edges, studies, populations, atlas, outcomes or {})
    if pvh_graph is not None:
        import_pvh(pvh_graph, nodes, edges, studies, atlas)
    graph = {"schema": SCHEMA, "schema_version": "1.0.0", "nodes": sorted(nodes.values(), key=lambda n: n["id"]),
             "edges": sorted(edges, key=lambda e: e["id"]), "studies": sorted(studies.values(), key=lambda s: s["study_id"]),
             "edge_layers": EDGE_LAYERS, "sources": sources or [], "limitations": LIMITATIONS,
             "legacy_curated_records_imported": len(records), "automatic_marker_or_proximity_joins": 0,
             "cell_coordinates_inferred": 0, "claim": "Public-source computational research evidence layer",
             "curated_source_scope": curated.get("scope"), "curated_source_limitations": curated.get("limitations", []),
             "curated_source_missing_data": curated.get("missing_data", {}),
             "pvh_sampling_metadata": ({key: deepcopy(pvh_graph.get(key)) for key in
                                         ("schema", "scope", "unknown_pool_nuclei", "not_present", "provenance")}
                                        if pvh_graph is not None else None)}
    validate_graph(graph)
    graph["coverage"] = coverage_of(graph)
    return graph


def build(curated_path, records_path=ROOT / "data/curated/records.jsonl", pvh_path=ROOT / "data/processed/pvh_population_pilot_v1/evidence_graph.json",
          atlas_path=ROOT / "viewer/data/allen_mapping.json", populations_path=ROOT / "data/vocab/populations.json", out_dir=OUT,
          outcomes_path=ROOT / "data/vocab/outcomes.json"):
    paths = [Path(curated_path), Path(records_path), Path(atlas_path), Path(populations_path), Path(outcomes_path)]
    if pvh_path is not None:
        paths.append(Path(pvh_path))
    curated = read(curated_path)
    # Bind acquired source text/receipts when supplied; never fetch a missing source.
    for study in curated.get("studies", []):
        if study.get("full_text_path"):
            path = ROOT / study["full_text_path"]
            if not path.resolve().is_relative_to(ROOT.resolve()) or sha(path) != study.get("full_text_sha256"):
                raise ValidationError("local full-text provenance hash differs: " + study["study_id"])
            paths.append(path)
        for relative in ([study["retrieval_receipt"]] if study.get("retrieval_receipt") else []) + study.get("additional_source_paths", []):
            path = ROOT / relative
            if not path.resolve().is_relative_to(ROOT.resolve()):
                raise ValidationError("source provenance must stay within the local research repository")
            paths.append(path)
    records = [json.loads(line) for line in Path(records_path).read_text().splitlines() if line.strip()]
    atlas = {r["region_id"]: r for r in read(atlas_path)["regions"]}
    graph = build_graph(curated, records, read(populations_path)["entries"], atlas,
                        read(pvh_path) if pvh_path is not None else None,
                        [{"path": str(p), "sha256": sha(p)} for p in paths], read(outcomes_path)["entries"])
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, payload in (("graph.json", graph), ("coverage.json", graph["coverage"])):
        temporary = out_dir / (name + ".tmp")
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
        temporary.replace(out_dir / name)
    return graph


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    create = sub.add_parser("build")
    create.add_argument("--curated-sources", type=Path, required=True)
    create.add_argument("--records", type=Path, default=ROOT / "data/curated/records.jsonl")
    create.add_argument("--pvh-graph", type=Path, default=ROOT / "data/processed/pvh_population_pilot_v1/evidence_graph.json")
    create.add_argument("--without-pvh", action="store_true")
    create.add_argument("--out-dir", type=Path, default=OUT)
    check = sub.add_parser("validate")
    check.add_argument("graph", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            result = validate_graph(read(args.graph))
        else:
            graph = build(args.curated_sources, args.records, None if args.without_pvh else args.pvh_graph, out_dir=args.out_dir)
            result = {"graph": str(args.out_dir / "graph.json"), "coverage": str(args.out_dir / "coverage.json"),
                      **validate_graph(graph)}
        print(json.dumps(result, sort_keys=True))
        return 0
    except (ValidationError, OSError, KeyError, json.JSONDecodeError) as error:
        print(json.dumps({"valid": False, "error": str(error)}), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
