"""D-097: choose the panel's five readouts and the papers to read, by a fixed rule.

    python3 analysis/panel_select.py

Reads the panel harvest (analysis/panel_harvest.py; its queries are re-issued from the
local cache, which gives Europe PMC's relevance order without a network call) and the
current corpus. Uses only where findings and papers exist, never what any of them says:
a circuit's existing finding counts as filled if it has a defined direction (a tie does
not), which is the only use of a direction here.

The rule, fixed in D-097 before any paper is read:
  1. A circuit is covered on a readout if its activation finding is filled now, or the
     harvest found at least one open-access paper for it that no curated record cites.
  2. The panel is the five readouts covering the most circuits; ties go to the readout
     with more new papers summed over circuits, then alphabetically.
  3. The reading set is, for every panel readout and every circuit not yet filled on it,
     the first CAP new papers in Europe PMC's relevance order. A paper found for several
     cells is read once.
Writes reports/panel_selection.md and data/interim/panel_reading_set.json.
"""
import collections
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
import d095_pair_test as T  # noqa: E402
import panel_harvest as H  # noqa: E402
from syn4rch.harvest import europepmc  # noqa: E402
from syn4rch.harvest.cache import Cache, RawLog  # noqa: E402

PANEL_SIZE = 5
CAP = 3


def main():
    rows, _, _, _, _, _ = T.load_records()
    prim = T.primary_rows(rows)
    filled = collections.defaultdict(set)              # readout -> circuits with a defined finding
    by = collections.defaultdict(list)
    for x in prim:
        if x.cell[1] == "circuit_activity" and x.cell[2] == "increase":
            by[(x.p, x.cell[0])].append(x.direction)
    for (p, o), ds in by.items():
        if T.modal(ds):
            filled[o].add(p)

    cache, raw = Cache(), RawLog()
    # "curated" was set by the harvest: a paper some curated record cites
    papers = {p["pmcid"]: p for p in json.load(open(os.path.join(ROOT, "data/interim/panel_candidates.json")))}
    new_ranked = {}                                    # (circuit, readout) -> new pmcids, relevance order
    for pid in H.CIRCUITS:
        for oid in H.READOUTS:
            recs, _ = europepmc.search(H.build_query(pid, oid), page_size=H.PAGE, cache=cache, raw_log=raw)
            order = [(r.get("pmcid") or "").upper() for r in recs if r.get("pmcid")]
            new_ranked[(pid, oid)] = [x for x in order if x in papers and not papers[x]["curated"]]

    cover = {}
    for oid in H.READOUTS:
        c = {p for p in H.CIRCUITS if p in filled[oid] or new_ranked[(p, oid)]}
        n_new = sum(len(new_ranked[(p, oid)]) for p in H.CIRCUITS)
        cover[oid] = (len(c), n_new, sorted(c))
    ranking = sorted(H.READOUTS, key=lambda o: (-cover[o][0], -cover[o][1], o))
    panel = ranking[:PANEL_SIZE]

    reading, why = [], collections.defaultdict(list)
    for oid in panel:
        for pid in H.CIRCUITS:
            if pid in filled[oid]:
                continue
            for pmc in new_ranked[(pid, oid)][:CAP]:
                why[pmc].append("%s/%s" % (pid, oid))
                if pmc not in reading:
                    reading.append(pmc)

    out = [{"pmcid": pmc, "doi": papers[pmc]["doi"], "title": papers[pmc]["title"],
            "year": papers[pmc]["year"], "for_cells": why[pmc]} for pmc in reading]
    with open(os.path.join(ROOT, "data/interim/panel_reading_set.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)

    with open(os.path.join(ROOT, "reports/panel_selection.md"), "w", encoding="utf-8") as fh:
        fh.write("# Panel selection (D-097)\n\nProduced by `python3 analysis/panel_select.py` before any "
                 "selected paper was read. Rule in the script's docstring and in D-097.\n\n")
        fh.write("## Readouts, by circuits covered\n\n| readout | filled now | covered | new papers |"
                 " circuits covered |\n|---|---|---|---|---|\n")
        for oid in ranking:
            fh.write("| %s%s | %d | %d | %d | %s |\n"
                     % (oid, " **(panel)**" if oid in panel else "", len(filled[oid] & set(H.CIRCUITS)),
                        cover[oid][0], cover[oid][1], ", ".join(cover[oid][2])))
        fh.write("\n## The panel: filled now (F), new papers to read (number), nothing (-)\n\n")
        fh.write("| circuit | " + " | ".join(panel) + " |\n|---|" + "---|" * len(panel) + "\n")
        for pid in H.CIRCUITS:
            cells = []
            for oid in panel:
                if pid in filled[oid]:
                    cells.append("F")
                else:
                    k = min(CAP, len(new_ranked[(pid, oid)]))
                    cells.append(str(k) if k else "-")
            fh.write("| %s | %s |\n" % (pid, " | ".join(cells)))
        n_cells = len(H.CIRCUITS) * len(panel)
        n_filled = sum(1 for oid in panel for p in H.CIRCUITS if p in filled[oid])
        n_try = sum(1 for oid in panel for p in H.CIRCUITS if p not in filled[oid] and new_ranked[(p, oid)])
        fh.write("\n%d cells: %d filled now, %d with papers to read, %d with none.\n"
                 % (n_cells, n_filled, n_try, n_cells - n_filled - n_try))
        fh.write("\n## Reading set: %d papers (at most %d per empty cell)\n\n" % (len(reading), CAP))
        for o in out:
            fh.write("* %s (%s) -- %s -- for %s\n" % (o["pmcid"], o["year"], (o["title"] or "")[:110],
                                                    ", ".join(o["for_cells"])))
    print("panel:", panel)
    print("cells %d, filled %d, to read %d, none %d; reading set %d papers"
          % (n_cells, n_filled, n_try, n_cells - n_filled - n_try, len(reading)))


if __name__ == "__main__":
    main()
