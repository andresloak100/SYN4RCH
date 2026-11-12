"""Summarize the completed evidence-map development work without fitting models.

The full local neuron scenes and raw source files are deliberately not bundled.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads((ROOT / path).read_text())


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def main():
    coverage = read("data/processed/circuit_evidence_v1/coverage.json")
    fly = read("reports/fly_circuit_benchmark.json")
    mouse = read("reports/fly_mouse_transfer.json")
    references = [
        "analysis/circuit_evidence_layer.py",
        "analysis/circuit_evidence_viewer.py",
        "src/syn4rch/atlas/circuit_evidence_template.html",
        "data/interim/circuit_evidence_v1/curated_sources.json",
        "data/interim/circuit_evidence_v1/atlas_context_manifest.json",
        "data/processed/circuit_evidence_v1/graph.json",
        "data/processed/circuit_evidence_v1/coverage.json",
        "reports/circuit_evidence_validation.json",
        "reports/fly_circuit_benchmark.json",
        "reports/fly_circuit_benchmark_replay.json",
        "reports/fly_mouse_transfer.json",
        "reports/fly_mouse_transfer_replay.json",
    ]
    replays = [read(p) for p in references if p.endswith("_replay.json")]
    assert all(v["status"] == "PASS" for v in replays)
    assert read("reports/circuit_evidence_validation.json")["status"] == "PASS"
    report = {
        "schema": "synarch_circuit_evidence_development_summary_v1",
        "status": "COMPLETED_DEVELOPMENT_ONLY",
        "map": {
            "local_viewer": "reports/circuit_evidence_viewer.html",
            "local_scene": "reports/circuit_evidence_neurons.html",
            "measured_example_neurons": 37,
            "nodes": coverage["nodes"],
            "evidence_links": coverage["edges"],
            "studies": coverage["studies"],
            "links_by_kind": coverage["edges_by_kind"],
            "existing_records_reused": 722,
            "existing_studies_reused": 109,
            "identified_synapses": 0,
            "individual_cells_matched_across_datasets": 0,
            "curation_status": "All literature records remain unreviewed; source fidelity was independently audited.",
            "placement": "Mouse region meshes only; no study population is assigned to an example traced neuron. Rat, unknown-species and fly anatomy are not placed in the mouse scene.",
            "counts_meaning": "Evidence links include paired intervention-target and intervention-outcome records. They are not independent observations or synapses.",
        },
        "fly_graph_completion": {
            "specimens": 1,
            "nodes": fly["nodes"],
            "hidden_source_profiles": fly["query_sources"],
            "scores": fly["scores"],
            "beats_all_fixed_controls": fly["beats_every_control"],
            "scope": "Within-specimen development on FlyWire v630; no independent-animal validation. AP is average precision, not percent accuracy.",
            "permitted_inputs": "Known incoming connections; whole hidden outgoing profiles removed before feature construction. No hidden axon geometry, genes or target labels used as query inputs.",
        },
        "fly_mouse_transfer": {
            "mouse_specimens": 1,
            "region": "visual cortex L2/3, not hypothalamus",
            "nodes": mouse["nodes"],
            "hidden_source_profiles": mouse["query_sources"],
            "scores": mouse["scores"],
            "finding": "Direct frozen fly parameter transfer performs worse than mouse-only training. A fly-informed slope prior has only a tiny descriptive gain over mouse-only shared_input_logistic: AP about +0.52% relative; log loss about -0.017%; Brier about -0.018%.",
            "scope": "Nonblind development on one mouse, retaining metadata-defined isolated neurons. No hypothalamic or independent-animal transfer claim.",
            "biological_transfer": "No fly cell identities, anatomical edges or neurotransmitter signs are copied into mouse biology.",
        },
        "source_studies": [
            {
                "doi": "10.1002/advs.78177",
                "scope": "Published mouse LH-to-VTA mechanogenetic circuit evidence; no aging or lifespan endpoint. Functional evidence is not relabeled as identified synapses.",
            },
            {
                "doi": "10.1016/j.cmet.2023.12.011",
                "scope": "Published chemogenetic DMH Ppp1r17 lifespan result in males, n=13 per group, log-rank p=0.0406. Chemogenetic rather than GV/ultrasound; lifespan animals were separate from the physiological-assay cohorts. No gas-vesicle longevity result.",
            },
        ],
        "verification": {
            "evidence_layer_unit_checks": 25,
            "viewer_unit_checks": 41,
            "fly_synthetic_checks": 15,
            "mouse_transfer_synthetic_checks": 17,
            "independent_source_and_metric_replays": "PASS",
            "viewer_browser_checks": [
                "3D scene visible by default at 595x719",
                "DMH selected and highlighted on open",
                "VTA display-context binding highlights the verified Allen region",
                "Original measured scene controls remain accessible",
                "Paper DOI links, p value, N scope and cohort-specific age displayed",
            ],
        },
        "next_scientific_requirements": [
            "Review individual source assertions and retain study-specific controls, animal identities and uncertainty.",
            "Evaluate any functional predictor on held-out intervention outcomes and independent studies or animals, with simple controls.",
            "Establish same-system connection and functional response measurements before interpreting an anatomical path as a predictable intervention outcome.",
        ],
        "limits": [
            "No complete hypothalamic circuit reconstruction or new validated 3D axon predictor is established here.",
            "These data do not demonstrate gas-vesicle-mediated lifespan extension, aging reversal or clinical efficacy.",
            "The full generated 3D viewer remains local because its measured trace licences restrict redistribution.",
            "This summary does not distribute raw papers, raw synapse data or licensed neuron geometry.",
        ],
        "artifact_hashes": {path: sha(path) for path in references},
    }
    out = ROOT / "reports/circuit_evidence_development_summary.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"path": str(out.relative_to(ROOT)), "status": report["status"], "sha256": sha(out.relative_to(ROOT))}))


if __name__ == "__main__":
    main()
