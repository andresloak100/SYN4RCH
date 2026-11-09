"""D-097: retrieve the reading set's full text, then stage each paper's units for a reader.

    python3 analysis/panel_stage.py fetch            # full text into data/raw/fulltext/
    python3 -m syn4rch.curate.run units               # stages 4-5, rebuilds extraction_units.jsonl
    python3 analysis/panel_stage.py stage OUT_DIR    # one JSON per paper, as D-072/D-087 staged them

Staging follows the earlier batches: a paper's extraction units in Results, Abstract and
unlabelled ("Other") passages, highest unit score first, at most 40, written with the
paper's identifiers. The staged file names no panel cell: a reader sees passages, never
the targeting (D-097).
"""
import collections
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from syn4rch.harvest import europepmc  # noqa: E402
from syn4rch.harvest.cache import Cache, RawLog  # noqa: E402

READING_SET = os.path.join(ROOT, "data/interim/panel_reading_set.json")
UNITS = os.path.join(ROOT, "data/interim/extraction_units.jsonl")
RESULTY = {"Results", "Other", "Abstract"}
CAP = 40


def fetch():
    cache, raw = Cache(), RawLog()
    got = missing = 0
    for p in json.load(open(READING_SET)):
        path = europepmc.fetch_fulltext(p["pmcid"], cache=cache, raw_log=raw)
        if path:
            got += 1
        else:
            missing += 1
            print("  no full text: %s" % p["pmcid"])
    print("full text: %d on disk, %d unavailable" % (got, missing))


def stage(out_dir):
    reading = {p["pmcid"]: p for p in json.load(open(READING_SET))}
    units = collections.defaultdict(list)
    with open(UNITS, encoding="utf-8") as fh:
        for line in fh:
            u = json.loads(line)
            if (u.get("pmcid") or "").upper() in reading:
                units[u["pmcid"].upper()].append(u)
    os.makedirs(out_dir, exist_ok=True)
    manifest = []
    for pmcid, meta in reading.items():
        rows = [u for u in units.get(pmcid, []) if u.get("section_label") in RESULTY]
        rows.sort(key=lambda u: (-(u.get("score") or 0), u.get("locator") or ""))
        given = rows[:CAP]
        manifest.append({"pmcid": pmcid, "units_available": len(units.get(pmcid, [])),
                         "units_resulty": len(rows), "units_given": len(given)})
        if not given:
            continue
        payload = {
            "pmcid": pmcid,
            "study_id": given[0].get("study_id"),
            "doi": given[0].get("doi") or meta.get("doi"),
            "title": meta.get("title"),
            "year": meta.get("year"),
            "units": [{k: u.get(k) for k in ("unit_id", "locator", "section_label", "passage_text",
                                              "candidate_population_id", "candidate_method_id",
                                              "candidate_outcome_id", "candidate_region_ids")}
                      for u in given],
        }
        with open(os.path.join(out_dir, "%s.json" % pmcid), "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=1, ensure_ascii=False)
    with open(os.path.join(out_dir, "_manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=1)
    staged = [m for m in manifest if m["units_given"]]
    print("reading set %d: staged %d papers with %d units; %d with no readable unit"
          % (len(manifest), len(staged), sum(m["units_given"] for m in staged), len(manifest) - len(staged)))
    for m in manifest:
        if not m["units_given"]:
            print("  no unit in Results/Abstract/Other: %s (units anywhere: %d)"
                  % (m["pmcid"], m["units_available"]))


if __name__ == "__main__":
    if sys.argv[1:2] == ["fetch"]:
        fetch()
    elif sys.argv[1:2] == ["stage"] and len(sys.argv) == 3:
        stage(sys.argv[2])
    else:
        sys.exit(__doc__)
