"""D-097: where the literature is, for the shared readout panel.

    python3 analysis/panel_harvest.py

D-095's data target needs the same readouts measured on the same circuits. This script
asks Europe PMC, for each of the 12 circuits with a connectivity fingerprint and each
candidate readout, one question: which open-access research papers name the circuit and
the readout in their title or abstract and use a circuit-activation method anywhere?

Nothing is read and nothing is scored. The output is where papers exist, which is what
D-097 chooses the panel from, before any of them is read:

    reports/panel_harvest.md            every query, its hit count, what it returned
    data/interim/panel_candidates.json  the papers, with the cells that found them

Uses the project's own Europe PMC client, cache and pacing (syn4rch.harvest).
"""
import collections
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from syn4rch.harvest import europepmc, queries as Q  # noqa: E402
from syn4rch.harvest.cache import Cache, RawLog  # noqa: E402
from syn4rch.schema.io import read_jsonl  # noqa: E402

# The 12 circuits D-095 places by connectivity and that hold findings.
CIRCUITS = ("arc_agrp", "arc_kisspeptin", "arc_pomc", "lha_mch", "lha_orexin", "lha_vgat",
            "poa_galanin", "pvh_crh", "pvh_oxt", "scn_vip", "tmn_histamine", "vmh_sf1")

# One candidate readout per physiological system (D-097 picks five of these).
READOUTS = ("food_intake", "corticosterone", "wakefulness", "locomotor_activity",
            "body_weight", "blood_glucose", "energy_expenditure", "core_temperature",
            "heart_rate")

# Search terms, written out so the query can be read and rerun. Taken from
# data/vocab/outcomes.json, minus acronyms and words too broad to search on.
READOUT_TERMS = {
    "food_intake": ["food intake", "feeding", "food consumption", "caloric intake",
                    "hyperphagia", "hypophagia"],
    "corticosterone": ["corticosterone", "HPA axis", "glucocorticoid", "adrenocorticotropic hormone"],
    "wakefulness": ["wakefulness", "arousal", "sleep", "NREM", "REM sleep"],
    "locomotor_activity": ["locomotor activity", "locomotion", "physical activity",
                           "spontaneous activity", "wheel running", "home cage activity"],
    "body_weight": ["body weight", "body mass", "weight gain", "weight loss"],
    "blood_glucose": ["blood glucose", "glycemia", "glycaemia", "plasma glucose",
                      "hyperglycemia", "hypoglycemia"],
    "energy_expenditure": ["energy expenditure", "oxygen consumption", "metabolic rate",
                           "indirect calorimetry"],
    "core_temperature": ["body temperature", "core temperature", "hypothermia", "hyperthermia",
                         "thermogenesis"],
    "heart_rate": ["heart rate", "bradycardia", "tachycardia", "blood pressure", "cardiovascular"],
}

# Circuit-activity methods (data/vocab/methods.json), short acronyms written out where
# they are specific enough to search on (ChR2).
ACTIVATION_TERMS = ["chemogenetic", "chemogenetically", "DREADD", "DREADDs", "hM3Dq",
                    "Gq-DREADD", "clozapine-N-oxide", "optogenetic", "optogenetically",
                    "channelrhodopsin", "ChR2", "photostimulation", "Chrimson",
                    "optogenetic activation", "chemogenetic activation"]

FILTERS = 'OPEN_ACCESS:y AND NOT PUB_TYPE:"review"'
PAGE = 100


def circuit_terms(pid):
    # 'NPY neurons' names a different population outside the arcuate
    return [t for t in Q.population_terms(pid) if t != "NPY neurons"]


def field(block):
    return "(TITLE:%s OR ABSTRACT:%s)" % (block, block)


def build_query(pid, oid):
    P = Q.epmc_block(circuit_terms(pid))
    O = Q.epmc_block(READOUT_TERMS[oid])
    A = Q.epmc_block(ACTIVATION_TERMS)
    return "%s AND %s AND %s AND %s" % (field(P), A, field(O), FILTERS)


def main():
    cache, raw = Cache(), RawLog()
    curated = {(r.pmcid or "").upper() for r in read_jsonl(os.path.join(ROOT, "data/curated/records.jsonl"))
               if r.pmcid}
    papers = {}
    log = []
    for pid in CIRCUITS:
        for oid in READOUTS:
            qid = "PANEL-%s-%s" % (pid, oid)
            term = build_query(pid, oid)
            try:
                recs, hits = europepmc.search(term, page_size=PAGE, cache=cache, raw_log=raw)
                err = None
            except Exception as exc:                   # recorded, never dropped
                recs, hits, err = [], 0, str(exc)[:200]
            found = []
            for r in recs:
                pmc = (r.get("pmcid") or "").upper()
                if not pmc:
                    continue
                found.append(pmc)
                p = papers.setdefault(pmc, {"pmcid": pmc, "doi": r.get("doi"), "pmid": r.get("pmid"),
                                            "title": r.get("title"), "year": r.get("year") or r.get("pub_year"),
                                            "journal": r.get("journal"), "cells": [], "curated": pmc in curated})
                p["cells"].append([pid, oid])
            log.append({"qid": qid, "circuit": pid, "readout": oid, "term": term, "hits": hits,
                        "returned": len(recs), "with_pmcid": len(found),
                        "new": sum(1 for x in found if x not in curated), "error": err})
            sys.stderr.write("  %-40s hits %4d  new %3d%s\n" % (qid, hits, log[-1]["new"], "  ERROR" if err else ""))
            sys.stderr.flush()

    os.makedirs(os.path.join(ROOT, "data/interim"), exist_ok=True)
    with open(os.path.join(ROOT, "data/interim/panel_candidates.json"), "w", encoding="utf-8") as fh:
        json.dump(sorted(papers.values(), key=lambda p: p["pmcid"]), fh, indent=1, ensure_ascii=False)

    grid = collections.defaultdict(dict)
    for e in log:
        grid[e["readout"]][e["circuit"]] = e
    with open(os.path.join(ROOT, "reports/panel_harvest.md"), "w", encoding="utf-8") as fh:
        fh.write("# Panel harvest (D-097)\n\n")
        fh.write("Produced by `python3 analysis/panel_harvest.py`. One Europe PMC query per circuit and "
                 "readout: the circuit and the readout in the title or abstract, a circuit-activation "
                 "method anywhere, open access, not a review. `new` counts returned papers with a PMCID "
                 "that no curated record cites. Nothing here has been read.\n\n")
        fh.write("## New open-access papers per circuit and readout\n\n")
        fh.write("| readout | " + " | ".join(CIRCUITS) + " |\n")
        fh.write("|---|" + "---|" * len(CIRCUITS) + "\n")
        for oid in READOUTS:
            fh.write("| %s | %s |\n" % (oid, " | ".join(str(grid[oid][c]["new"]) for c in CIRCUITS)))
        fh.write("\n%d distinct papers returned, %d not cited by any curated record.\n\n"
                 % (len(papers), sum(1 for p in papers.values() if not p["curated"])))
        fh.write("## Every query\n\n")
        for e in log:
            fh.write("### %s\n\nhits %d, returned %d, with PMCID %d, new %d%s\n\n```\n%s\n```\n\n"
                     % (e["qid"], e["hits"], e["returned"], e["with_pmcid"], e["new"],
                        (", ERROR: " + e["error"]) if e["error"] else "", e["term"]))
    print("queries %d, errors %d, papers %d, new %d"
          % (len(log), sum(1 for e in log if e["error"]), len(papers),
             sum(1 for p in papers.values() if not p["curated"])))


if __name__ == "__main__":
    main()
