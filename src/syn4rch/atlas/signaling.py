"""Derive the measured-signalling graph the viewer plays back.

What this graph is
------------------
One node on the left is a hypothalamic cell population. One node on the right
is a systemic physiological readout. An edge exists only because *N curated
records* say a perturbation of that population was followed by a measured
change in that readout. The edge is therefore a **count of published
observations**, not a measured flow of signal and not a prediction.

What this graph is **not**
--------------------------
It is not a projection map. ``data/curated/records.jsonl`` carries a
``projection_target_id`` on a handful of records only (the count is recomputed
here as ``n_records_with_projection_target``), so no nucleus-to-nucleus wiring
can be drawn from this dataset and none is. Nothing here is integrated,
simulated or extrapolated: every number in the output is a count of rows.

Derivation rules (each one is a decision, not a convenience)
------------------------------------------------------------
``axis`` -- the signed description of what the perturbation did to the
population. Taken from the record's ``perturbation_polarity`` when that field
exists and is non-null; otherwise it falls back to ``valence``. This module
was written against schema v0.1.0, where no record carried a polarity and
every record took the ``valence`` path; v0.2.0 added the field and the
derivation picked it up with no code change, which is also why the fallback
stays -- a record written before the migration still keys identically.

``unclear`` is carried through untouched and is never recoded to a sign. A
polarity is only ever read on its own class's axis, never translated to
another: a ``conditional_knockout`` carrying ``decrease`` means less of the
target gene product, and this module never restates that as less firing,
because the records do not say what the knockout did to the firing.

:func:`syn4rch.schema.rules.perturbation_polarity_from` is deliberately *not*
used as a second fallback. It can infer a sign from the method alone (a
knockout is a loss of function on the molecular axis), and that inference is
correct on its own axis -- but it is an inference, and this view animates what
the records state. Only a polarity stored on the record is read here.

**Polarity trap.** From schema v0.2.0 ``perturbation_polarity`` takes the
values ``increase`` / ``decrease`` / ``remove`` / ``bidirectional`` /
``unclear``, which are the same two words as the outcome ``direction``
vocabulary while meaning something completely different: ``axis = "increase"``
is more firing (or more gene product) in the *brain*, while
``direction = "increase"`` is a *bodily readout* going up. The two must never
be shown side by side bare, so every edge also carries ``axis_phrase``, which
spells the axis out in its class's own terms, and that is what the page shows.
See docs/DECISIONS.md (D-021, D-022).

``axis_class`` -- taken from the record's ``perturbation_class`` when that
field exists and is non-null, otherwise mapped from ``method`` by
:data:`METHOD_CLASS`.

Edge identity -- ``(population_id, outcome_id, axis, axis_class)``. Axis and
axis class are part of the key, not an average over it: chemogenetically
activating AgRP neurons and knocking a gene out of them are different
perturbations and must not be merged into one arrow.

``consensus`` -- the modal ``direction`` of the edge's records, or ``None``
when two or more directions tie for the most records. It is the direction most
records report, not a statistical majority; read it together with
``agreement`` (= ``max(directions) / sum(directions)``).

``contested`` -- **one rule**: ``agreement < 0.7`` and ``n_records >= 3``
(:data:`CONTESTED_RULE`). Fewer than three records is too few to call a
disagreement; at three or more, an edge whose modal direction holds less than
70% of its records is drawn split rather than averaged.

The obvious alternative -- "two or more directions, each reported by two or
more distinct studies" -- was measured against this corpus and rejected as the
contested rule. It flags only a couple of further edges, and one of them is
``arc_agrp -> food_intake`` under "more firing": 21 records say the mice ate
more, 3 say no change, over 9 studies, and it is the most replicated claim in
the dataset.
Marking that contested would misrepresent agreement as conflict, which is the
same failure as hiding conflict, pointed the other way. The signal is kept
instead as the separate, separately named ``minority_replicated`` flag, which
is reported in the provenance readout and is not drawn as conflict.

``n_studies`` / ``n_labs`` -- counted with the project's own identity
functions, :func:`syn4rch.qc.keys.study_key` and
:func:`syn4rch.qc.keys.lab_key`. A record whose lab identity is unknown
(``lab_key`` returns ``None``) is not counted as a lab and is not lumped in
with other unknowns, so ``n_labs`` can be smaller than ``n_studies``.

Python 3.10 compatible, standard library only (docs/DECISIONS.md D-001).
"""

from __future__ import annotations

import datetime as _dt
import json
import os
from typing import Any, Dict, List, Optional, Tuple

from ..qc.keys import lab_key, population_key, study_key

__all__ = [
    "CONTESTED_RULE",
    "MINORITY_REPLICATED_RULE",
    "METHOD_CLASS",
    "OUTCOME_SYSTEM",
    "SYSTEM_ORDER",
    "RECORD_IDS_CAP",
    "axis_of",
    "axis_class_of",
    "axis_phrase",
    "build_graph",
    "write_graph",
]

#: Readable statement of the contested rule, copied into the output so the
#: page can show the definition next to the flag.
CONTESTED_RULE = "contested = agreement < 0.7 and n_records >= 3"

#: The secondary flag. Not the contested rule, and never drawn as conflict:
#: it says a minority direction was reported by more than one study, which is
#: worth reading next to a high agreement score.
MINORITY_REPLICATED_RULE = (
    "minority_replicated = two or more directions, each reported by two or "
    "more distinct studies"
)

#: ``method`` -> ``axis_class``. Used only when the record carries no
#: ``perturbation_class`` of its own. This mirrors
#: :func:`syn4rch.schema.rules.perturbation_class_from_method`; it is kept
#: local so this module reads one file and imports nothing that is still
#: moving, and ``tests/test_atlas_signaling.py`` fails if the two ever
#: disagree.
METHOD_CLASS = {
    "optogenetic": "circuit_activity",
    "chemogenetic": "circuit_activity",
    "electrical_stimulation": "circuit_activity",
    "conditional_knockout": "molecular",
    "conditional_overexpression": "molecular",
    "lesion": "structural",
    "genetic_ablation": "structural",
    "local_pharmacology": "pharmacological",
}
AXIS_CLASS_FALLBACK = "other"

#: Presentational grouping ONLY. This is a layout aid for the viewer -- it
#: decides which outcomes sit near each other on the page and nothing else. It
#: is **not** a vocabulary claim and it is deliberately not written back into
#: ``data/vocab/outcomes.json``, which owns the outcome vocabulary.
OUTCOME_SYSTEM = {
    "food_intake": "feeding",
    "body_weight": "body_composition",
    "fat_mass": "body_composition",
    "lean_mass": "body_composition",
    "energy_expenditure": "energy_expenditure",
    "rer": "energy_expenditure",
    "glucose_tolerance": "glucose_metabolism",
    "insulin_sensitivity": "glucose_metabolism",
    "blood_glucose": "glucose_metabolism",
    "plasma_insulin": "glucose_metabolism",
    "plasma_glucagon": "glucose_metabolism",
    "hepatic_glucose_production": "glucose_metabolism",
    "hepatic_lipid": "glucose_metabolism",
    "wat_lipolysis": "energy_expenditure",
    "bat_thermogenesis": "thermoregulation",
    "core_temperature": "thermoregulation",
    "sympathetic_nerve_activity": "autonomic",
    "blood_pressure": "autonomic",
    "heart_rate": "autonomic",
    "corticosterone": "endocrine",
    "thyroid_hormone": "endocrine",
    "osmolality_avp": "endocrine",
    "gonadotropins": "reproduction",
    "estrous_cyclicity": "reproduction",
    "locomotor_activity": "behaviour",
    # Added when the circuit curation pass brought the first records with these
    # readouts (D-051). The map is asserted total against the real dataset, so a
    # new outcome fails the suite rather than quietly vanishing from the atlas --
    # which is the guard working, not a bug in the pass.
    "sleep": "behaviour",
    "gut_motility": "autonomic",
    "bone_density": "body_composition",
    # The second circuit batch brought the first inflammatory readouts: PVH CRH
    # inhibition lowering IL-6 and TNF in an LPS model. `inflammatory_markers_other`
    # is already placed; these two are the named cytokines.
    "il6": "endocrine",
    "tnf": "endocrine",
    "inflammatory_markers_other": "endocrine",
    "immune_cell_composition": "endocrine",
}

#: Order the systems are laid out in. Also presentational.
SYSTEM_ORDER = [
    "feeding",
    "body_composition",
    "energy_expenditure",
    "glucose_metabolism",
    "thermoregulation",
    "autonomic",
    "endocrine",
    "reproduction",
    "behaviour",
]

SYSTEM_LABEL = {
    "feeding": "Feeding",
    "body_composition": "Body composition",
    "energy_expenditure": "Energy expenditure",
    "glucose_metabolism": "Glucose and lipid metabolism",
    "thermoregulation": "Thermoregulation",
    "autonomic": "Autonomic outflow",
    "endocrine": "Endocrine axes",
    "reproduction": "Reproduction",
    "behaviour": "Behaviour",
}

#: At most this many ``record_id`` values are listed per edge. The cap is
#: reported in the output as ``record_ids_cap``.
RECORD_IDS_CAP = 25

#: Exclusion reasons. Every input record lands in exactly one of these or in
#: ``n_records_used``.
EX_NO_POP = "no resolved population_id (region-level record, not cell-type-specific)"
EX_NO_OUTCOME = "no outcome_id"
EX_NO_DIRECTION = "no direction reported"
EX_POP_NOT_IN_VOCAB = "population_id is not in the population vocabulary"


class _RecordView:
    """Attribute view over a decoded record, for the ``qc.keys`` functions.

    ``qc.keys`` reads records by attribute. Wrapping the decoded JSON instead
    of building a ``TierARecord`` keeps this module independent of the
    dataclass's enum coercion -- a half-reported row that the container would
    reject still has a perfectly good study and lab identity, and identity is
    all that is wanted here.
    """

    __slots__ = ("study_id", "population_id", "population_raw",
                 "last_author", "last_author_affiliation")

    def __init__(self, rec: Dict[str, Any]) -> None:
        self.study_id = rec.get("study_id")
        self.population_id = rec.get("population_id")
        self.population_raw = rec.get("population_raw")
        self.last_author = rec.get("last_author")
        self.last_author_affiliation = rec.get("last_author_affiliation")


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _read_jsonl(path: str) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def _vocab(path: str) -> Dict[str, Dict[str, Any]]:
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as fh:
        return json.load(fh).get("entries") or {}


def axis_of(rec: Dict[str, Any]) -> Optional[str]:
    """``perturbation_polarity`` when the record carries one, else ``valence``."""
    pol = rec.get("perturbation_polarity")
    if pol:
        return pol
    return rec.get("valence")


def axis_class_of(rec: Dict[str, Any]) -> str:
    """``perturbation_class`` when the record carries one, else from ``method``."""
    cls = rec.get("perturbation_class")
    if cls:
        return cls
    return METHOD_CLASS.get(rec.get("method"), AXIS_CLASS_FALLBACK)


def system_of(outcome_id: str) -> Optional[str]:
    return OUTCOME_SYSTEM.get(outcome_id)


#: ``(axis_class, axis)`` -> plain words for the axis. The schema v0.2.0
#: polarity values ``increase`` and ``decrease`` are the same two words as the
#: outcome ``direction`` values and mean something else entirely, so no edge is
#: ever labelled with the bare value. An unknown combination falls through to
#: the raw value with the word "perturbation" attached, which is wrong-free: it
#: says less rather than guessing.
AXIS_PHRASE = {
    ("circuit_activity", "increase"): "more firing",
    ("circuit_activity", "decrease"): "less firing",
    ("circuit_activity", "activation"): "activation (more firing)",
    ("circuit_activity", "inhibition"): "inhibition (less firing)",
    ("molecular", "increase"): "more of the target gene product",
    ("molecular", "decrease"): "less of the target gene product",
    ("pharmacological", "activation"): "drug-driven activation",
    ("pharmacological", "inhibition"): "drug-driven inhibition",
    ("pharmacological", "increase"): "drug-driven increase",
    ("pharmacological", "decrease"): "drug-driven decrease",
}
AXIS_PHRASE_ANY = {
    "remove": "population or gene product removed",
    "ablation": "population ablated",
    "bidirectional": "run both ways; the record covers both",
    "unclear": "sign of the manipulation not established",
    "activation": "activation",
    "inhibition": "inhibition",
}


def axis_phrase(axis: Optional[str], axis_class: Optional[str]) -> str:
    """The axis in its own class's words -- never the bare polarity value."""
    hit = AXIS_PHRASE.get((axis_class or "", axis or ""))
    if hit:
        return hit
    hit = AXIS_PHRASE_ANY.get(axis or "")
    if hit:
        return hit
    return "perturbation recorded as %r" % (axis,)


def _consensus(directions: Dict[str, int]) -> Tuple[Optional[str], float]:
    """(modal direction or None on a tie, agreement = max/sum)."""
    total = sum(directions.values())
    if not total:
        return None, 0.0
    top = max(directions.values())
    winners = sorted(k for k, v in directions.items() if v == top)
    agreement = top / total
    return (winners[0] if len(winners) == 1 else None), agreement


def _contested(directions: Dict[str, int], agreement: float) -> bool:
    """:data:`CONTESTED_RULE`, and nothing else."""
    return sum(directions.values()) >= 3 and agreement < 0.7


def _minority_replicated(studies_by_direction: Dict[str, set]) -> bool:
    """:data:`MINORITY_REPLICATED_RULE`. Reported, never drawn as conflict."""
    return sum(1 for d in studies_by_direction.values() if len(d) >= 2) >= 2


def build_graph(records_path: str, populations_path: str, outcomes_path: str) -> Dict[str, Any]:
    """Build the graph. Every count comes from ``records_path``."""
    records = _read_jsonl(records_path)
    pop_vocab = _vocab(populations_path)
    out_vocab = _vocab(outcomes_path)

    # Guard the presentational mapping before anything is emitted: a system
    # that is not in the layout order would silently drop its outcomes.
    bad_systems = sorted({s for s in OUTCOME_SYSTEM.values() if s not in SYSTEM_ORDER})
    assert not bad_systems, "system(s) not in SYSTEM_ORDER: %s" % ", ".join(bad_systems)

    excluded: Dict[str, int] = {}

    def drop(reason: str) -> None:
        excluded[reason] = excluded.get(reason, 0) + 1

    pop_acc: Dict[str, Dict[str, Any]] = {}
    out_acc: Dict[str, Dict[str, Any]] = {}
    edge_acc: Dict[Tuple[str, str, str, str], Dict[str, Any]] = {}
    used = 0
    n_projection = 0
    seen_outcomes: set = set()
    dir_seen: Dict[str, int] = {}

    for rec in records:
        if rec.get("projection_target_id"):
            n_projection += 1
        pid = rec.get("population_id")
        oid = rec.get("outcome_id")
        direction = rec.get("direction")
        view = _RecordView(rec)
        sk = study_key(view)

        # One reason per record, first match wins, so the buckets partition the
        # input exactly. A record with neither a population nor an outcome is
        # counted as region-level, which is the more informative of the two.
        reason: Optional[str] = None
        if not pid:
            reason = EX_NO_POP
        elif pop_vocab and pid not in pop_vocab:
            reason = EX_POP_NOT_IN_VOCAB
        elif not oid:
            reason = EX_NO_OUTCOME
        elif not direction:
            reason = EX_NO_DIRECTION

        if oid:
            # Readouts are listed even when every record naming them was
            # excluded, so the page can say "4 records, all region-level"
            # rather than quietly omitting the readout.
            seen_outcomes.add(oid)
            o = out_acc.setdefault(oid, {"records": 0, "studies": set(), "excluded": 0})
            if reason:
                o["excluded"] += 1
            else:
                o["records"] += 1
                o["studies"].add(sk)

        if reason:
            drop(reason)
            continue

        lk = lab_key(view)
        axis = axis_of(rec) or "unclear"
        cls = axis_class_of(rec)
        used += 1

        p = pop_acc.setdefault(pid, {"records": 0, "studies": set()})
        p["records"] += 1
        p["studies"].add(sk)

        key = (pid, oid, axis, cls)
        e = edge_acc.setdefault(key, {
            "directions": {}, "studies": set(), "labs": set(),
            "studies_by_direction": {}, "record_ids": [], "n_records": 0,
        })
        e["n_records"] += 1
        e["directions"][direction] = e["directions"].get(direction, 0) + 1
        dir_seen[direction] = dir_seen.get(direction, 0) + 1
        e["studies"].add(sk)
        if lk is not None:
            e["labs"].add(lk)
        e["studies_by_direction"].setdefault(direction, set()).add(sk)
        rid = rec.get("record_id")
        if rid and len(e["record_ids"]) < RECORD_IDS_CAP:
            e["record_ids"].append(rid)

    # Every outcome that reached the graph must have a system.
    missing = sorted(o for o in seen_outcomes if system_of(o) is None)
    assert not missing, "outcome(s) with no system in OUTCOME_SYSTEM: %s" % ", ".join(missing)

    # Populations: every vocabulary population is listed so the viewer can say
    # "no curated record yet" from this file rather than from nothing. The
    # counts are still only ever counts of rows in records_path.
    populations: List[Dict[str, Any]] = []
    for pid in sorted(set(pop_vocab) | set(pop_acc)):
        v = pop_vocab.get(pid) or {}
        acc = pop_acc.get(pid) or {"records": 0, "studies": set()}
        populations.append({
            "id": pid,
            "label": v.get("label") or pid,
            "region_id": v.get("nucleus"),
            "neurotransmitter": v.get("neurotransmitter") or [],
            "marker_genes": v.get("marker_genes") or [],
            "n_records": acc["records"],
            "n_studies": len(acc["studies"]),
            "in_vocabulary": pid in pop_vocab,
        })
    populations.sort(key=lambda r: (-r["n_records"], r["id"]))

    outcomes: List[Dict[str, Any]] = []
    for oid in sorted(out_acc):
        v = out_vocab.get(oid) or {}
        outcomes.append({
            "id": oid,
            "label": v.get("label") or oid,
            "tier": v.get("tier"),
            "modality": v.get("modality"),
            "system": system_of(oid),
            "n_records": out_acc[oid]["records"],
            "n_studies": len(out_acc[oid]["studies"]),
            "n_records_excluded": out_acc[oid]["excluded"],
            "in_vocabulary": oid in out_vocab,
        })
    outcomes.sort(key=lambda r: (SYSTEM_ORDER.index(r["system"]), -r["n_records"], r["id"]))

    edges: List[Dict[str, Any]] = []
    for (pid, oid, axis, cls), e in edge_acc.items():
        consensus, agreement = _consensus(e["directions"])
        edges.append({
            "population_id": pid,
            "outcome_id": oid,
            "axis": axis,
            "axis_class": cls,
            "axis_phrase": axis_phrase(axis, cls),
            "n_records": e["n_records"],
            "n_studies": len(e["studies"]),
            "n_labs": len(e["labs"]),
            "directions": dict(sorted(e["directions"].items())),
            "consensus": consensus,
            "agreement": round(agreement, 3),
            "contested": _contested(e["directions"], agreement),
            "minority_replicated": _minority_replicated(e["studies_by_direction"]),
            "record_ids": e["record_ids"],
            "record_ids_truncated": e["n_records"] > len(e["record_ids"]),
        })
    edges.sort(key=lambda r: (-r["n_records"], r["population_id"], r["outcome_id"], r["axis"]))

    n_in = len(records)
    assert used + sum(excluded.values()) == n_in, (
        "records are not partitioned: %d used + %d excluded != %d in"
        % (used, sum(excluded.values()), n_in))
    for e in edges:
        assert e["n_records"] >= 1 and sum(e["directions"].values()) == e["n_records"]

    n_unclear = sum(1 for e in edges if e["axis"] == "unclear")
    n_contested = sum(1 for e in edges if e["contested"])
    n_minority = sum(1 for e in edges if e["minority_replicated"] and not e["contested"])
    n_pop_with_records = sum(1 for p in populations if p["n_records"])
    n_excluded_no_pop = excluded.get(EX_NO_POP, 0)

    honesty = [
        "Nothing on this page is predicted. Every arrow, pulse and number is a "
        "count of curated observations already published in %d papers."
        % len({study_key(_RecordView(r)) for r in records}),
        "An arrow is a count of published observations, not a measurement of "
        "signal flow. It says that perturbing this population was followed by a "
        "measured change in this readout in N records; it does not say that a "
        "signal travels this way, and the animation's speed and timing carry no "
        "physiological units.",
        "No neural projection is drawn. Only %d of %d curated records name a "
        "projection target, so this dataset cannot support a wiring diagram and "
        "none is shown: the right-hand nodes are bodily readouts, not brain "
        "regions." % (n_projection, n_in),
        "%d of %d records are region-level (no cell-type-specific population) "
        "and are excluded from this view entirely. They are real evidence about "
        "a nucleus; they just cannot be attached to a population."
        % (n_excluded_no_pop, n_in),
        "%d of the %d edges are unpolarised: the records do not establish "
        "whether the manipulation raised or lowered what it acted on, so they "
        "are drawn as dashed channels with hollow pulses and must not be read "
        "as switching the neurons on or off. The readout itself was still "
        "measured, and its direction is still shown." % (n_unclear, len(edges)),
        "The axis and the direction can read as the same two words and are not "
        "the same statement. The axis is what the experiment did inside the "
        "brain -- more firing, less of a gene product, the population removed. "
        "The direction is what the body did. Every axis on this page is "
        "therefore written out in words rather than left as its bare value.",
        "%d of the %d edges are contested, meaning the records disagree about "
        "the direction. They are drawn split rather than averaged, and the "
        "split is the real record count on each side. Rule: %s. A further %d "
        "edges agree well overall but have a minority direction that more than "
        "one study reported; that is flagged in the readout, not drawn as "
        "conflict." % (n_contested, len(edges), CONTESTED_RULE, n_minority),
        "%d of the %d vocabulary populations carry no curated record at all. "
        "That is a statement about how far curation has got, not about what the "
        "literature contains." % (len(populations) - n_pop_with_records, len(populations)),
    ]

    return {
        "generated_at": _now(),
        "generated_from": "data/curated/records.jsonl",
        "mode_label": "Measured signalling - evidence playback",
        "mode_note": ("Recorded observations animated. Nothing is simulated: no "
                      "equation is integrated and no unmeasured response is "
                      "inferred."),
        "n_records_in": n_in,
        "n_records_used": used,
        "n_records_with_projection_target": n_projection,
        "excluded": dict(sorted(excluded.items())),
        "axis_source": ("perturbation_polarity when the record carries it, "
                        "otherwise valence"),
        "axis_class_source": ("perturbation_class when the record carries it, "
                              "otherwise mapped from method"),
        "edge_identity": "population_id x outcome_id x axis x axis_class",
        "axis_trap_note": ("The axis describes the manipulation in the brain. "
                           "The direction describes the bodily readout. From "
                           "schema v0.2.0 both vocabularies use the words "
                           "'increase' and 'decrease', so the axis is always "
                           "shown as a phrase (axis_phrase) and never as the "
                           "bare value."),
        "directions_seen": dict(sorted(dir_seen.items())),
        "consensus_rule": ("the modal direction; null when two or more "
                           "directions tie for the most records"),
        "agreement_rule": "max(directions) / sum(directions)",
        "contested_rule": CONTESTED_RULE,
        "minority_replicated_rule": MINORITY_REPLICATED_RULE,
        "record_ids_cap": RECORD_IDS_CAP,
        "system_order": list(SYSTEM_ORDER),
        "system_labels": dict(SYSTEM_LABEL),
        "system_note": ("Systems group the readouts for layout only. They are "
                        "not part of the outcome vocabulary and nothing is "
                        "written back to data/vocab/outcomes.json."),
        "populations": populations,
        "outcomes": outcomes,
        "edges": edges,
        "n_edges": len(edges),
        "n_edges_contested": n_contested,
        "n_edges_minority_replicated": n_minority,
        "n_edges_unclear_axis": n_unclear,
        "n_populations_with_records": n_pop_with_records,
        "_honesty": honesty,
    }


def write_graph(graph: Dict[str, Any], path: str) -> str:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(graph, fh, indent=1)
    return path


def main(argv=None) -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(os.path.dirname(os.path.dirname(here)))
    g = build_graph(
        os.path.join(root, "data", "curated", "records.jsonl"),
        os.path.join(root, "data", "vocab", "populations.json"),
        os.path.join(root, "data", "vocab", "outcomes.json"),
    )
    out = write_graph(g, os.path.join(root, "viewer", "data", "signaling.json"))
    print("%d records in, %d used, %d excluded" % (g["n_records_in"], g["n_records_used"],
                                                   sum(g["excluded"].values())))
    for reason, n in g["excluded"].items():
        print("  excluded %4d  %s" % (n, reason))
    print("%d edges (%d contested, %d unclear axis) over %d populations and %d outcomes"
          % (g["n_edges"], g["n_edges_contested"], g["n_edges_unclear_axis"],
             g["n_populations_with_records"], len(g["outcomes"])))
    print("wrote %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
