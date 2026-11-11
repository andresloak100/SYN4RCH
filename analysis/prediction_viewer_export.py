"""Export interface between a PASSED biological prediction test and the local 3D viewer.

    python3 -I analysis/prediction_viewer_export.py paths --bundle DIR --out FILE [--n-initial 12]
    python3 -I analysis/prediction_viewer_export.py regional --predictions FILE --out FILE

It draws nothing. It turns an already scored evaluation into one JSON file that the existing local
viewer builder (analysis/swc_render.py) can load as an overlay, and it REFUSES unless the success
condition was actually met.

`paths`: DIR is a bundle scored once by the independent referee
(analysis/prediction_independent_evaluator.py: declaration.json, prediction_manifest.json,
PREDICTION_SEAL.json, metrics.json, SCORE_COMPLETE.json, predictions/<model>/<id>.npz,
truth/<id>.npz). It is exported only if
  * every hash still agrees (declaration, forecasts, seal, metrics, saved truth);
  * the declared parts are the parts of the locked cohort plan, and every TEST neuron of the plan
    is either scored or listed in the declaration's exclusions (no silent dropping);
  * the candidate beat EVERY declared copying control on both animal-mean measures, and a matched
    control was declared and beaten too.
Each neuron then gets three layers: measured reconstruction, model prediction, and the strongest
copying control. All scored neurons are listed, so failures stay inspectable.

`regional`: for the case where only a regional prediction task passes. Predicted regional values
become a per-region heat map labelled as a regional prediction. They are never turned into paths.

Nothing here opens a trace of a held-out animal before it was scored: `paths` reads only what the
referee saved after scoring. Stdlib only.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "reports/axon_dataset_plan.md"
COHORT = ROOT / "reports/axon_dataset_manifest.tsv"
STRUCTURES = ROOT / "data/raw/atlas/structure_graph_1.json"
EXAMPLE_SALT = "SYN4RCH-viewer-examples-v1"
SCHEMA = "synarch_viewer_export_v1"
LARGE = 2000
ROW_CORE = ("query_id", "model", "path", "sha256", "donor_ids", "upstream_fit_animals")
EVIDENCE = {
    "independent_test": "independent test",
    "exploratory_training_partition": "exploratory: cross-validation inside the training part, not an independent test",
}
NOT_SYNAPSES = "A path is not a synapse: no connection between cells is shown or implied."
LOCAL_ONLY = ("Local viewing only. Measured traces of held-out animals come from sources that do not permit "
              "redistribution (portal: all rights reserved; data centre: CC BY-NC-ND 4.0).")
FRAME_RULE = ("The three layers of a neuron share the declared frame and must receive the same display transform. In the "
              "mirrored-left frame (midline ML 5,700 um) a neuron whose soma lies in the Right hemisphere is displayed at "
              "ML' = 11,400 - ML.")


class ExportRefused(Exception):
    """The success condition or an integrity check was not met; nothing is exported."""


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def cohort(plan=PLAN, manifest=COHORT):
    """neuron -> (animal, part), from the cohort manifest whose SHA-256 the committed plan pins."""
    plan, manifest = Path(plan), Path(manifest)
    if not manifest.exists():
        raise ExportRefused("cohort manifest missing: regenerate it with analysis/axon_dataset_audit.py (%s)" % manifest)
    if sha(manifest) not in plan.read_text():
        raise ExportRefused("the cohort manifest is not the one pinned by the locked plan")
    rows = [line.rstrip("\n").split("\t") for line in open(manifest)]
    head = rows[0]
    n, a, p = head.index("neuron"), head.index("animal"), head.index("part")
    return {r[n]: (r[a], r[p]) for r in rows[1:]}


def check_agreed_test(decl, locked):
    """The bundle must be the locked test: declared parts are the plan's parts, nothing dropped silently."""
    for part, rows in decl["splits"].items():
        for row in rows:
            got = locked.get(row["id"])
            if got is None:
                raise ExportRefused("%s is not a neuron of the locked cohort" % row["id"])
            if got != (row["animal"], part):
                raise ExportRefused("%s is declared as %s of brain sample %s, but the locked plan has it in %s of %s"
                                    % (row["id"], part, row["animal"], got[1], got[0]))
    scored = {row["id"] for row in decl["splits"]["test"]}
    excluded = set(decl.get("excluded_ids", []))
    unaccounted = sorted(n for n, (_, part) in locked.items() if part == "test" and n not in scored and n not in excluded)
    if unaccounted:
        raise ExportRefused("not the whole locked test: %d test neurons are neither scored nor in the declared exclusions "
                            "(first: %s)" % (len(unaccounted), ", ".join(unaccounted[:5])))
    return {"plan_test_neurons": sum(part == "test" for _, part in locked.values()),
            "plan_test_brain_samples": len({a for a, part in locked.values() if part == "test"}),
            "scored_neurons": len(scored), "scored_brain_samples": len({row["animal"] for row in decl["splits"]["test"]}),
            "declared_exclusions_in_test": sorted(n for n in excluded if locked.get(n, ("", ""))[1] == "test")}


def order_key(salt, value):
    return hashlib.sha256((salt + ":" + value).encode()).hexdigest()


def initial_examples(queries, n):
    """Fixed rule that uses identifiers only, never a result: brain samples in hash order, and in each
    sample the neuron first in hash order; one neuron per sample, n samples."""
    by_animal = {}
    for q in queries:
        by_animal.setdefault(q["animal"], []).append(q["id"])
    animals = sorted(by_animal, key=lambda a: order_key(EXAMPLE_SALT, a))[:n]
    return [min(by_animal[a], key=lambda i: order_key(EXAMPLE_SALT, i)) for a in animals]


def strongest_copy(summary, copies):
    """Fixed rule: the copying control with the lowest animal-mean distance; ties by higher F1, then name."""
    return min(copies, key=lambda m: (summary[m]["distance_um"], -summary[m]["f1_200"], m))


def verdict(distance_delta, f1_delta):
    if distance_delta < 0 and f1_delta > 0:
        return "prediction better on both measures"
    if distance_delta >= 0 and f1_delta <= 0:
        return "prediction not better on either measure"
    return "mixed: better on one measure only"


def load_scored_bundle(base):
    base = Path(base)
    need = ("declaration.json", "prediction_manifest.json", "PREDICTION_SEAL.json", "metrics.json", "SCORE_COMPLETE.json")
    missing = [n for n in need if not (base / n).exists()]
    if missing:
        raise ExportRefused("not a completed referee evaluation: missing " + ", ".join(missing))
    if (base / "SCORE_FAILED.json").exists():
        raise ExportRefused("the referee recorded a failed scoring (SCORE_FAILED.json)")
    done, seal, metrics = read(base / "SCORE_COMPLETE.json"), read(base / "PREDICTION_SEAL.json"), read(base / "metrics.json")
    if done.get("metrics_sha256") != sha(base / "metrics.json") or done.get("seal_sha256") != sha(base / "PREDICTION_SEAL.json"):
        raise ExportRefused("metrics or seal changed after scoring")
    if seal.get("manifest_sha256") != sha(base / "prediction_manifest.json") or seal.get("declaration_sha256") != sha(base / "declaration.json"):
        raise ExportRefused("manifest or declaration changed after sealing")
    if metrics.get("status") != "scored_once" or metrics.get("seal_sha256") != done.get("seal_sha256"):
        raise ExportRefused("metrics.json does not record one completed scoring of this seal")
    return base, read(base / "declaration.json"), read(base / "prediction_manifest.json"), metrics


def export_paths(bundle, out, n_initial=12, plan=PLAN, manifest=COHORT):
    base, decl, forecasts, metrics = load_scored_bundle(bundle)
    roles, overall = decl["models"], metrics["overall"]
    candidate, copies, matched = roles["candidate"], list(roles["copy_baselines"]), list(roles.get("matched_controls", []))
    coverage = check_agreed_test(decl, cohort(plan, manifest))
    queries = metrics["queries"]
    if {q["id"] for q in queries} != {row["id"] for row in decl["splits"]["test"]}:
        raise ExportRefused("the scored neurons are not the declared test neurons")
    if not overall.get("beats_all_copying"):
        raise ExportRefused("success condition not met: the candidate did not beat every declared copying control on both "
                            "animal-mean distance and F1. Nothing is exported.")
    if not overall.get("matched_controls_available") or not overall.get("beats_all_matched"):
        raise ExportRefused("success condition not met: the agreed test also requires a declared matched control, beaten on "
                            "both measures. Nothing is exported.")
    summary = overall["summary"]
    baseline = strongest_copy(summary, copies)
    best_f1 = max(copies, key=lambda m: (summary[m]["f1_200"], -summary[m]["distance_um"], m))
    rows = {(r["query_id"], r["model"]): r for r in forecasts["predictions"]}
    known = {q["id"]: q for q in decl["splits"]["test"]}
    first = set(initial_examples(queries, n_initial))
    status = {candidate: "model prediction (not measured)"}
    status.update({m: "copying control: a measured TRAINING neuron moved to this soma (not this neuron)" for m in copies})
    status.update({m: "matched control (not measured)" for m in matched})

    def layer(qid, model):
        r = rows[qid, model]
        if sha(base / r["path"]) != r["sha256"]:
            raise ExportRefused("forecast file changed: " + r["path"])
        fit = r.get("upstream_fit_animals")
        extra, omitted = {}, []
        for k in sorted(set(r) - set(ROW_CORE)):
            if len(json.dumps(r[k])) > LARGE:
                omitted.append(k)
            else:
                extra[k] = r[k]
        return {"model": model, "status": status[model], "path": r["path"], "sha256": r["sha256"], "donor_ids": r["donor_ids"],
                "upstream_fit_animals": None if fit is None else
                {"count": len(fit), "sha256": hashlib.sha256("\n".join(sorted(fit)).encode()).hexdigest()},
                "provenance": extra, "provenance_fields_omitted_as_large": omitted}

    neurons = []
    for q in sorted(queries, key=lambda q: q["id"]):
        if sha(base / q["truth_path"]) != q["truth_sha256"]:
            raise ExportRefused("saved truth changed: " + q["truth_path"])
        mine, theirs = q["models"][candidate], q["models"][baseline]
        dd, df = mine["distance_um"] - theirs["distance_um"], mine["f1_200"] - theirs["f1_200"]
        neurons.append({
            "id": q["id"], "held_out_brain_sample": q["animal"], "source_label": q.get("peptide"), "soma_label": q.get("nucleus"),
            "known_soma_um": known[q["id"]]["soma_um"], "hemisphere": rows[q["id"], candidate].get("hemisphere", "unknown"),
            "initial_example": q["id"] in first,
            "layers": {
                "measured": {"status": "measured reconstruction (held-out truth, opened only for scoring)",
                             "path": q["truth_path"], "sha256": q["truth_sha256"]},
                "prediction": layer(q["id"], candidate),
                "copy_baseline": layer(q["id"], baseline),
            },
            "other_controls": {m: layer(q["id"], m) for m in [*copies, *matched] if m != baseline},
            "metrics": q["models"],
            "comparison": {"against": baseline, "distance_delta_um": dd, "f1_delta": df, "verdict": verdict(dd, df)},
            "uncertainty": "uncalibrated: no interval with validated coverage accompanies this forecast",
        })
    export = {
        "schema": SCHEMA, "kind": "paths", "bundle": str(base.resolve()),
        "benchmark_status": {
            "statement": "Independent test of the locked cohort, scored once by the independent referee. Candidate '%s' beat all "
                         "%d declared copying controls and all %d matched controls on animal-mean distance and F1."
                         % (candidate, len(copies), len(matched)),
            "evidence_level": "independent_test", "beats_all_copying": True, "beats_all_matched": True,
            "neurons": overall["neurons"], "brain_samples": overall["independent_animals"], "coverage": coverage,
            "animal_mean": {m: summary[m] for m in [candidate, *copies, *matched]}, "against": overall["against"],
            "interpretation": "descriptive: one test, no claim beyond the declared rule; intervals are bootstrap over brain samples",
            "unit_caveat": "A brain sample ID is treated as one animal: 166 IDs and 166 reported mice, no explicit one-to-one statement.",
            "referee_limits": metrics.get("limits"),
        },
        "roles": {"candidate": candidate, "copy_baselines": copies, "matched_controls": matched, "copy_layer": baseline,
                  "copy_layer_rule": "copying control with the lowest animal-mean distance (ties: higher F1, then name)",
                  "best_f1_copy": best_f1},
        "frame": {"geometry_contract": decl.get("geometry_contract"), "rule": FRAME_RULE},
        "hashes": {"declaration_sha256": sha(base / "declaration.json"), "manifest_sha256": sha(base / "prediction_manifest.json"),
                   "seal_sha256": sha(base / "PREDICTION_SEAL.json"), "metrics_sha256": sha(base / "metrics.json"),
                   "referee_source_sha256": metrics.get("source_sha256"), "method_artifacts": decl.get("method_artifacts"),
                   "cohort_manifest_sha256": sha(manifest)},
        "example_rule": {"salt": EXAMPLE_SALT, "n": n_initial, "all_neurons_listed": True,
                         "rule": "brain samples in SHA-256 order of salt:sample; in each, the neuron first in SHA-256 order of "
                                 "salt:neuron. Identifiers only: no result enters the choice."},
        "labels": {"measured": "Measured reconstruction", "prediction": "Model prediction (not measured)",
                   "copy_baseline": "Copying baseline (a measured training neuron moved to this soma)"},
        "notices": [NOT_SYNAPSES, LOCAL_ONLY, "Predictions and copies are never drawn in the style of measured anatomy."],
        "neurons": neurons,
    }
    _write(out, export)
    return export


def export_regional(predictions, out, structures=STRUCTURES):
    """Input contract, written by a regional task that passed its frozen benchmark:
    {"task": str, "scale": str, "regions": [CCF acronyms],
     "status": {"beats_all_baselines": bool, "evidence_level": "independent_test" | "exploratory_training_partition",
                "benchmark": str, "protocol_sha256": str, "results_sha256": str, "model_sha256": str},
     "units": [{"unit_id", "group_id", "predicted": [...], "interval_low": [...]?, "interval_high": [...]?,
                "coverage_validated": bool?}]}"""
    p = read(predictions)
    status = p.get("status", {})
    if status.get("beats_all_baselines") is not True:
        raise ExportRefused("success condition not met: the regional task did not beat every declared baseline. Nothing is exported.")
    if status.get("evidence_level") not in EVIDENCE:
        raise ExportRefused("the regional result does not say whether it is an independent test or exploratory")
    for k in ("benchmark", "protocol_sha256", "results_sha256", "model_sha256"):
        if not status.get(k):
            raise ExportRefused("the regional result does not carry its benchmark status and hashes (%s)" % k)
    ids = {}

    def walk(node):
        ids[node["acronym"]] = node["id"]
        for child in node.get("children", []):
            walk(child)
    walk(read(structures)["msg"][0])
    unknown = [r for r in p["regions"] if r not in ids]
    if unknown:
        raise ExportRefused("regions that are not CCF structures: " + ", ".join(unknown))
    units = []
    for u in p["units"]:
        for k in ("predicted", "interval_low", "interval_high"):
            if u.get(k) is not None and len(u[k]) != len(p["regions"]):
                raise ExportRefused("a unit's %s values do not match the region list" % k)
        validated = u.get("coverage_validated") is True and u.get("interval_low") is not None and u.get("interval_high") is not None
        units.append({"unit_id": u["unit_id"], "group_id": u["group_id"], "predicted": u["predicted"],
                      "interval_low": u.get("interval_low"), "interval_high": u.get("interval_high"),
                      "uncertainty": "interval with validated coverage" if validated else "uncalibrated: coverage not validated"})
    export = {"schema": SCHEMA, "kind": "regional", "task": p["task"], "scale": p["scale"], "benchmark_status": status,
              "regions": [{"acronym": r, "ccf_structure_id": ids[r]} for r in p["regions"]],
              "labels": {"heatmap": "Regional prediction (%s). Model output per atlas region: not a measured value, not a path."
                         % EVIDENCE[status["evidence_level"]]},
              "notices": ["Regional values are shown on atlas regions only. They are not converted into axon paths.", NOT_SYNAPSES],
              "source_sha256": sha(predictions), "units": units}
    _write(out, export)
    return export


def _write(out, export):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(export, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"
    out.write_text(text)
    Path(str(out) + ".sha256").write_text(hashlib.sha256(text.encode()).hexdigest() + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("paths")
    a.add_argument("--bundle", required=True)
    a.add_argument("--out", required=True)
    a.add_argument("--n-initial", type=int, default=12)
    b = sub.add_parser("regional")
    b.add_argument("--predictions", required=True)
    b.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    try:
        if args.cmd == "paths":
            e = export_paths(args.bundle, args.out, args.n_initial)
            print("exported %d neurons (%d initial examples) to %s"
                  % (len(e["neurons"]), sum(n["initial_example"] for n in e["neurons"]), args.out))
        else:
            e = export_regional(args.predictions, args.out)
            print("exported regional predictions for %d units to %s" % (len(e["units"]), args.out))
    except ExportRefused as refusal:
        print("REFUSED: %s" % refusal)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
