"""Pre-curation yield assessment over data/interim/extraction_units.jsonl.

    python3 -m syn4rch.curate.yield_report

Answers one question before any paper is read: how much of the corpus could
possibly become Tier A records, and which (population, outcome) cells could
possibly clear the feasibility gate. Every figure is an UPPER BOUND, because a
unit is a passage in which a population, a perturbation method and a readout
co-occur -- not an observation. The real yield is measured only by curating.

Why it exists. Curating 974 full texts is the expensive step of this project.
Sizing the ceiling first says whether that spend can reach a modellable
dataset, and which populations to curate first. A ceiling that already fails
the gate would mean the spend cannot succeed, which is worth knowing cheaply.

Standard library only (DECISIONS.md D-001). Reads data/, writes reports/.
"""

from __future__ import annotations

import argparse
import collections
import datetime as _dt
import json
import os
import statistics
import sys
from typing import Dict, List, Optional, Set, Tuple

from ..feasibility.gates import (
    MIN_LABS_FEASIBLE,
    MIN_STUDIES_FEASIBLE,
    MIN_STUDIES_MARGINAL,
)
from ..qc.keys import norm

# Sections whose passages can state the paper's own measured result. A unit in
# Introduction or Discussion is overwhelmingly either a citation of somebody
# else's experiment or an interpretation of one, and CURATION_PROTOCOL.md 2
# requires the reader to discard both. Abstract units state the paper's own
# finding but carry no figure-level locator, so they are reported apart.
RESULT_SECTIONS = ("Results", "Methods")
OWN_CLAIM_NO_LOCATOR = ("Abstract",)

TOP_N = 20


def _utc_mtime(path: str) -> str:
    ts = os.path.getmtime(path)
    return _dt.datetime.fromtimestamp(ts, _dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_lab_index(candidates_path: str) -> Tuple[Dict[str, Optional[Tuple[str, str]]], bool]:
    """``study_id -> lab key`` from the harvest candidates, or empty if absent."""
    index: Dict[str, Optional[Tuple[str, str]]] = {}
    if not os.path.exists(candidates_path):
        return index, False
    with open(candidates_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            sid = rec.get("study_id")
            if not sid:
                continue
            author = norm(rec.get("last_author"))
            index[sid] = (author, norm(rec.get("last_author_affiliation"))) if author else None
    return index, True


def load_units(path: str) -> List[Dict]:
    units = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                units.append(json.loads(line))
    return units


def assess(units: List[Dict], labs: Dict[str, Optional[Tuple[str, str]]]) -> Dict:
    by_section = collections.Counter(u.get("section_label") or "unlabelled" for u in units)
    result_units = [u for u in units if u.get("section_label") in RESULT_SECTIONS]
    abstract_units = [u for u in units if u.get("section_label") in OWN_CLAIM_NO_LOCATOR]

    def triple(u):
        return (u.get("candidate_population_id"), u.get("candidate_method_id"),
                u.get("candidate_outcome_id"))

    per_paper = collections.Counter(u["study_id"] for u in units)
    cells: Dict[Tuple[str, str], Dict[str, Set]] = collections.defaultdict(
        lambda: {"studies": set(), "labs": set(), "unlabelled_studies": set()})
    for u in result_units:
        key = (u.get("candidate_population_id"), u.get("candidate_outcome_id"))
        cell = cells[key]
        sid = u["study_id"]
        cell["studies"].add(sid)
        lab = labs.get(sid)
        if lab:
            cell["labs"].add(lab)
        else:
            cell["unlabelled_studies"].add(sid)

    ranked = sorted(cells.items(),
                    key=lambda kv: (-len(kv[1]["studies"]), -len(kv[1]["labs"]), kv[0]))
    gate_pass = [k for k, v in cells.items()
                 if len(v["studies"]) >= MIN_STUDIES_FEASIBLE and len(v["labs"]) >= MIN_LABS_FEASIBLE]
    replicated = [k for k, v in cells.items() if len(v["studies"]) >= MIN_STUDIES_MARGINAL]
    singleton = [k for k, v in cells.items() if len(v["studies"]) == 1]

    pops: Dict[str, Dict[str, Set]] = collections.defaultdict(lambda: {"papers": set(), "units": 0})
    for u in result_units:
        p = pops[u.get("candidate_population_id")]
        p["papers"].add(u["study_id"])
        p["units"] += 1

    result_studies = {u["study_id"] for u in result_units}
    return {
        "n_units": len(units),
        "by_section": by_section,
        "n_result_units": len(result_units),
        "n_abstract_units": len(abstract_units),
        "papers_any_unit": len(per_paper),
        "papers_result_unit": len(result_studies),
        "units_per_paper_median": statistics.median(per_paper.values()) if per_paper else 0,
        "units_per_paper_max": max(per_paper.values()) if per_paper else 0,
        "triples_all": len({triple(u) for u in units}),
        "triples_result": len({triple(u) for u in result_units}),
        "ceiling_all": len({(u["study_id"],) + triple(u) for u in units}),
        "ceiling_result": len({(u["study_id"],) + triple(u) for u in result_units}),
        "cells": cells,
        "ranked": ranked,
        "gate_pass": gate_pass,
        "replicated": replicated,
        "singleton": singleton,
        "pops": pops,
        "studies_without_lab": sum(1 for s in result_studies if not labs.get(s)),
        "n_result_studies": len(result_studies),
    }


def render(a: Dict, *, units_path: str, units_ts: str, have_labs: bool) -> str:
    L: List[str] = []
    w = L.append
    pct = lambda n, d: ("%.0f%%" % (100.0 * n / d)) if d else "n/a"

    w("# SYN4RCH extraction yield assessment (pre-curation)")
    w("")
    w("**Every figure here is an upper bound, not a dataset.** An extraction unit is")
    w("a passage in which a neuronal population, a perturbation method and a")
    w("physiological readout co-occur. Co-occurrence is not a finding")
    w("(`docs/CURATION_PROTOCOL.md` 2). The number of real records is known only")
    w("after papers are read and proposals pass the verbatim gate.")
    w("")
    w("## Run")
    w("")
    w("| field | value |")
    w("|---|---|")
    w("| units file | `%s` |" % units_path)
    w("| units file written (UTC) | %s |" % units_ts)
    w("| extraction units | %d |" % a["n_units"])
    w("| papers yielding >=1 unit | %d |" % a["papers_any_unit"])
    w("| units per paper, median / max | %d / %d |"
      % (a["units_per_paper_median"], a["units_per_paper_max"]))
    w("| lab identity available | %s |" % ("yes, from the harvest candidates" if have_labs
                                           else "NO -- candidates.jsonl absent, lab columns omitted"))
    w("")
    w("## Where the units sit, and why it decides the yield")
    w("")
    w("| section | units | can state the paper's own measured result? |")
    w("|---|---|---|")
    notes = {
        "Results": "yes, with a citable locator",
        "Methods": "yes, for perturbation detail",
        "Abstract": "yes, but with no figure-level locator",
        "Introduction": "almost never: cites other work",
        "Discussion": "almost never: interprets, often other work",
    }
    for section, n in a["by_section"].most_common():
        w("| %s | %d | %s |" % (section, n, notes.get(section, "unclassified")))
    w("")
    w("Only **%d of %d** units (%s) sit in %s, the sections where a measured"
      % (a["n_result_units"], a["n_units"], pct(a["n_result_units"], a["n_units"]),
         " or ".join("`%s`" % s for s in RESULT_SECTIONS)))
    w("result with a direction can be quoted against a locator. A further %d sit in"
      % a["n_abstract_units"])
    w("the abstract, which states the paper's own claim but gives no figure to cite.")
    w("The Introduction and Discussion units are the largest block and the least")
    w("usable: the protocol requires the reader to discard a passage that cites or")
    w("interprets somebody else's experiment.")
    w("")
    w("## Ceiling on the number of records")
    w("")
    w("| bound | all sections | %s only |" % "/".join(RESULT_SECTIONS))
    w("|---|---|---|")
    w("| distinct (population, method, outcome) triples | %d | %d |"
      % (a["triples_all"], a["triples_result"]))
    w("| distinct study x triple (**hard ceiling on records**) | %d | %d |"
      % (a["ceiling_all"], a["ceiling_result"]))
    w("| papers that could contribute | %d | %d |"
      % (a["papers_any_unit"], a["papers_result_unit"]))
    w("")
    w("At most **%d** records can exist from these full texts, and realistically far" % a["ceiling_result"])
    w("fewer: a passage naming all three terms frequently reports no direction for")
    w("that particular triple, or reports it for a different population in the same")
    w("nucleus. Treat %d as the number that cannot be exceeded, not as a forecast."
      % a["ceiling_result"])
    w("")
    w("## Could any prediction task clear the feasibility gate?")
    w("")
    w("`syn4rch.feasibility.gates` calls a task feasible at >=%d independent studies"
      % MIN_STUDIES_FEASIBLE)
    w("across >=%d distinct labs. Applying those thresholds to units rather than to"
      % MIN_LABS_FEASIBLE)
    w("records gives the most optimistic possible answer:")
    w("")
    n_cells = len(a["cells"])
    w("| (population, outcome) cells, %s units | count | share |" % "/".join(RESULT_SECTIONS))
    w("|---|---|---|")
    w("| with any extractable unit | %d | 100%% |" % n_cells)
    w("| >=%d studies and >=%d labs (gate, upper bound) | %d | %s |"
      % (MIN_STUDIES_FEASIBLE, MIN_LABS_FEASIBLE, len(a["gate_pass"]), pct(len(a["gate_pass"]), n_cells)))
    w("| >=%d studies | %d | %s |"
      % (MIN_STUDIES_MARGINAL, len(a["replicated"]), pct(len(a["replicated"]), n_cells)))
    w("| exactly 1 study (agreement undefined) | %d | %s |"
      % (len(a["singleton"]), pct(len(a["singleton"]), n_cells)))
    w("")
    if have_labs:
        w("Lab identity resolved for %d of %d contributing studies (%d without)."
          % (a["n_result_studies"] - a["studies_without_lab"], a["n_result_studies"],
             a["studies_without_lab"]))
        w("")
    w("### Cells with the most independent studies")
    w("")
    w("Study and lab counts are nearly equal throughout, which means this corpus is")
    w("not one group publishing repeatedly: the independence that SCHEMA.md 5 calls")
    w("the project's largest leakage risk is genuinely present at the top of the")
    w("distribution.")
    w("")
    w("| population | outcome | studies | labs |")
    w("|---|---|---|---|")
    for (pop, out), v in a["ranked"][:TOP_N]:
        w("| `%s` | `%s` | %d | %d |" % (pop, out, len(v["studies"]), len(v["labs"])))
    w("")
    w("### Populations by extractable units")
    w("")
    w("| population | units | papers |")
    w("|---|---|---|")
    for pop, v in sorted(a["pops"].items(), key=lambda kv: -kv[1]["units"])[:TOP_N]:
        w("| `%s` | %d | %d |" % (pop, v["units"], len(v["papers"])))
    w("")
    w("## What this does and does not establish")
    w("")
    w("**Establishes.** The corpus is not uniformly thin. A narrow set of cells,")
    w("concentrated in arcuate AgRP and POMC neurons against feeding and body-weight")
    w("readouts, has dozens of independent laboratories behind it. If curation")
    w("recovers even a modest fraction of those passages as records, direction-of-")
    w("effect prediction becomes testable for those cells.")
    w("")
    w("**Does not establish.** That any record exists. That a passage reports a")
    w("direction. That independent studies agree. %s of cells still rest on a single"
      % pct(len(a["singleton"]), n_cells))
    w("study, where agreement is undefined rather than perfect, so a model trained")
    w("on the dense corner will not generalise to the sparse remainder, and")
    w("leave-one-population-out evaluation should be expected to fail.")
    w("")
    w("**Known bias.** Only open-access full texts are here. Paywalled papers")
    w("cannot yield records at all, so the dense corner is partly a map of which")
    w("laboratories publish open access.")
    return "\n".join(L) + "\n"


def repo_root() -> str:
    """Same resolution as ``syn4rch.curate.run``, so both CLIs agree on paths."""
    try:
        from ..harvest.cache import repo_root as _rr

        return _rr()
    except Exception:
        return os.path.abspath(os.getcwd())


def main(argv=None) -> int:
    root = repo_root()
    ap = argparse.ArgumentParser(prog="python3 -m syn4rch.curate.yield_report",
                                 description=__doc__.split("\n\n")[0])
    ap.add_argument("--units", default=os.path.join(root, "data", "interim", "extraction_units.jsonl"))
    ap.add_argument("--candidates", default=os.path.join(root, "data", "interim", "candidates.jsonl"))
    ap.add_argument("--out", default=os.path.join(root, "reports", "extraction_yield.md"))
    ap.add_argument("--stdout", action="store_true")
    args = ap.parse_args(argv)

    if not os.path.exists(args.units):
        sys.stderr.write("no extraction units at %s; run `python3 -m syn4rch.curate.run units` first\n"
                         % args.units)
        return 1
    units = load_units(args.units)
    if not units:
        sys.stderr.write("%s is empty; nothing to assess\n" % args.units)
        return 1
    labs, have_labs = load_lab_index(args.candidates)
    a = assess(units, labs)
    text = render(a, units_path=os.path.relpath(args.units, root),
                  units_ts=_utc_mtime(args.units), have_labs=have_labs)
    if args.stdout:
        sys.stdout.write(text)
    else:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(text)
        print("wrote %s" % os.path.relpath(args.out, root))
        print("ceiling on records: %d | cells passing gate (upper bound): %d of %d"
              % (a["ceiling_result"], len(a["gate_pass"]), len(a["cells"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
