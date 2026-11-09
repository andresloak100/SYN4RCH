"""Allen Mouse Brain Connectivity Atlas inputs for D-093 and D-095.

    python3 analysis/fetch_connectivity.py

Three steps, each cached under data/raw/atlas/connectivity/ (not version controlled) and
skipped when its output exists:

1. hy_injections.json -- every experiment whose primary injection is in the
   hypothalamus (structure 1097).
2. circuit_experiments.json -- each circuit's cell-type tracing experiments. Matched on
   the GENE the Cre line drives, not the construct's name (Hcrt-IRES-Cre and Hcrt-Cre
   are one gene; Sf1 and Nr5a1 are one gene), and on the circuit's nucleus. The circuit's
   Cre lines and nucleus come from data/vocab/populations.json: this is the step that
   says where a paper's circuit sits, the analogue of D-065 for HypoMap.
3. proj_<id>.json -- projection density in every structure, one experiment at a time
   with a pause: the Allen API is a shared public service.
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

D = "data/raw/atlas/connectivity"
API = "https://api.brain-map.org/api/v2/data/query.json?criteria="
SAFE = ":,[]$=()'"
HEADERS = {"User-Agent": "syn4rch-atlas/0.1"}

SYN = {"sf1": "nr5a1", "vgat": "slc32a1", "vglut2": "slc17a6"}
NUC = {"arc": {"ARH"}, "lha": {"LHA"}, "pvh": {"PVH", "PVHd", "PVp", "PVa", "PVi"},
       "vmh": {"VMH"}, "dmh": {"DMH"}, "tmn": {"TMv", "TMd", "TM"}, "scn": {"SCH"},
       "poa": {"MPO", "MPN", "LPO", "AVPV"}, "mpo": {"MPO", "MPN"}, "zi": {"ZI"},
       "avpv": {"AVPV", "MPO"}}


def query(criteria, timeout=180):
    req = urllib.request.Request(API + urllib.parse.quote(criteria, safe=SAFE), headers=HEADERS)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def gene(line):
    g = re.split(r"[-_ ]", (line or ""))[0].lower()
    return SYN.get(g, g)


def main():
    os.makedirs(D, exist_ok=True)

    dest = os.path.join(D, "hy_injections.json")
    if not os.path.exists(dest):
        d = query("service::mouse_connectivity_injection_structure"
                  "[injection_structures$eq1097][primary_structure_only$eqtrue]", timeout=120)
        json.dump(d, open(dest, "w"))
    rows = json.load(open(dest))["msg"]

    dest = os.path.join(D, "circuit_experiments.json")
    if not os.path.exists(dest):
        vocab = json.load(open("data/vocab/populations.json"))["entries"]
        match = {}
        for pid, e in vocab.items():
            genes = {gene(c) for c in (e.get("cre_lines") or [])}
            nuc = NUC.get((e.get("nucleus") or "").lower(), set())
            exps = [r["id"] for r in rows
                    if gene(r.get("transgenic-line")) in genes and r.get("structure-abbrev") in nuc]
            if exps:
                match[pid] = exps
        json.dump(match, open(dest, "w"), indent=1)
    match = json.load(open(dest))

    ids = sorted({e for v in match.values() for e in v})
    for i, eid in enumerate(ids, 1):
        dest = os.path.join(D, "proj_%d.json" % eid)
        if os.path.exists(dest):
            continue
        d = query("model::ProjectionStructureUnionize,rma::criteria,[section_data_set_id$eq%d],"
                  "[is_injection$eqfalse],[hemisphere_id$eq3],rma::options[num_rows$eqall]" % eid)
        if not d.get("success"):
            print("FAILED %d: %s" % (eid, str(d.get("msg"))[:120]), file=sys.stderr)
            continue
        json.dump(d["msg"], open(dest + ".part", "w"))
        os.replace(dest + ".part", dest)
        print("  %d/%d  experiment %d: %d structures" % (i, len(ids), eid, len(d["msg"])),
              file=sys.stderr, flush=True)
        time.sleep(1.0)
    print("%d circuits, %d experiments" % (len(match), len(ids)), file=sys.stderr)


if __name__ == "__main__":
    main()
