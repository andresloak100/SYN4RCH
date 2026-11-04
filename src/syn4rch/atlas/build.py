"""Build the SYN4RCH atlas viewer.

    PYTHONPATH=src python3 -m syn4rch.atlas.build [--offline]

1. Download (or reuse the cached) Allen structure graph 1.
2. Resolve every region vocabulary entry to an Allen id (exact matches only).
3. Download Allen CCF 2017 meshes for the whole brain (997) and every
   resolved region; simplify by vertex clustering.
4. Derive the measured-signalling graph from data/curated/records.jsonl.
5. Write viewer/data/{allen_mapping,scene,geometry,signaling}.json and inline
   the scene, the geometry and the signalling graph into viewer/index.html.

If the mesh host cannot be reached the build stops. No geometry is ever
substituted.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
import urllib.error
import urllib.request
from typing import Dict, List

from . import evidence, mesh, ontology, signaling

ONTOLOGY_URL = "https://api.brain-map.org/api/v2/structure_graph_download/1.json"
MESH_URL = ("https://download.alleninstitute.org/informatics-archive/current-release/"
            "mouse_ccf/annotation/ccf_2017/structure_meshes/{id}.obj")
BRAIN_ID = 997

# Vertex budgets after simplification (keeps the inlined page small).
BUDGET = {"outline": 16000, "context": 2600, "nucleus": 1400}
MIN_CELL_UM = {"outline": 40.0, "context": 30.0, "nucleus": 12.0}


def repo_root() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.dirname(os.path.dirname(os.path.dirname(here)))


class HostUnreachable(RuntimeError):
    pass


def fetch(url: str, dest: str, offline: bool = False, timeout: float = 60.0) -> str:
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    if offline:
        raise HostUnreachable("offline and not cached: %s" % url)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "syn4rch-atlas/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            if resp.status != 200:
                raise HostUnreachable("%s returned HTTP %s" % (url, resp.status))
            data = resp.read()
    except (urllib.error.URLError, OSError) as exc:
        raise HostUnreachable("%s: %s" % (url, exc))
    tmp = dest + ".part"
    with open(tmp, "wb") as fh:
        fh.write(data)
    os.replace(tmp, dest)
    return dest


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _dump(obj) -> str:
    return json.dumps(obj, separators=(",", ":"), ensure_ascii=False)


def inline_json(obj) -> str:
    """JSON safe to place inside a <script type=application/json> element."""
    return _dump(obj).replace("<", "\\u003c")


def build(offline: bool = False, root: str = None) -> Dict:
    root = root or repo_root()
    raw = os.path.join(root, "data", "raw", "atlas")
    out_dir = os.path.join(root, "viewer", "data")
    os.makedirs(out_dir, exist_ok=True)

    graph_path = fetch(ONTOLOGY_URL, os.path.join(raw, "structure_graph_1.json"), offline)
    with open(graph_path, encoding="utf-8") as fh:
        graph = json.load(fh)
    nodes = ontology.flatten(graph)

    with open(os.path.join(root, "data", "vocab", "regions.json"), encoding="utf-8") as fh:
        regions = json.load(fh)["entries"]
    with open(os.path.join(root, "data", "vocab", "populations.json"), encoding="utf-8") as fh:
        populations = json.load(fh)["entries"]

    mapping = ontology.resolve_regions(regions, graph)
    by_region = {m["region_id"]: m for m in mapping}

    # Population nuclei are region ids; check each one is a known region.
    pop_rows: List[Dict] = []
    for pid, p in populations.items():
        m = by_region.get(p["nucleus"])
        pop_rows.append({
            "id": pid, "label": p["label"], "nucleus": p["nucleus"],
            "allen_id": m["allen_id"] if m else None,
            "nucleus_status": m["status"] if m else "not_in_region_vocab",
            "marker_genes": p.get("marker_genes") or [],
            "neurotransmitter": p.get("neurotransmitter") or [],
            "cre_lines": p.get("cre_lines") or [],
            "curation_confidence": p.get("curation_confidence"),
            "notes": p.get("notes"),
        })

    # Test the mesh host on the whole-brain mesh first; stop if unreachable.
    mesh_dir = os.path.join(raw, "meshes")
    targets = [(BRAIN_ID, "outline", None)]
    for m in mapping:
        if m["status"] == "resolved":
            role = "context" if m["region_id"] == "hy" else "nucleus"
            targets.append((m["allen_id"], role, m["region_id"]))

    geometry: Dict[str, Dict] = {}
    structures: List[Dict] = []
    for sid, role, rid in targets:
        path = fetch(MESH_URL.format(id=sid), os.path.join(mesh_dir, "%d.obj" % sid), offline)
        with open(path, encoding="utf-8") as fh:
            v, t = mesh.parse_obj(fh)
        v2, t2, cell = mesh.simplify_to(v, t, BUDGET[role], min_cell=MIN_CELL_UM[role])
        enc = mesh.encode(v2, t2)
        geometry[str(sid)] = {"p": enc.pop("p"), "i": enc.pop("i")}
        node = ontology.find_by_id(nodes, sid)
        structures.append(dict(
            enc, allen_id=sid, allen_acronym=node["acronym"], allen_name=node["name"],
            allen_color=node["color"], role=role, region_id=rid,
            source_url=MESH_URL.format(id=sid),
            source_vertex_count=len(v), source_triangle_count=len(t),
            cluster_cell_um=round(cell, 1),
        ))
        print("  mesh %-10s %-6s %6d -> %5d verts  cell %.0f um"
              % (sid, node["acronym"], len(v), len(v2), cell))

    lit = evidence.literature(
        os.path.join(root, "data", "interim", "candidates.jsonl"),
        os.path.join(root, "reports", "harvest_summary.md"))
    regions_lit = evidence.region_literature(
        os.path.join(root, "data", "interim", "candidates.jsonl"),
        [m["region_id"] for m in mapping if m["status"] == "resolved"])
    causal = evidence.causal_records(os.path.join(root, "data", "curated", "records.jsonl"))
    kg_dir = os.path.join(root, "src", "syn4rch", "kg")
    kg_files = [f for f in os.listdir(kg_dir) if f.endswith(".py")] if os.path.isdir(kg_dir) else []

    scene = {
        "built_at": _now(),
        "ccf": {"version": "CCFv3 2017", "units": "micrometres",
                "axes": {"x": "anterior->posterior", "y": "dorsal->ventral", "z": "left->right"},
                "ontology_url": ONTOLOGY_URL, "mesh_url_template": MESH_URL},
        "placement": "nucleus-level: populations are attached to the Allen volume of their nucleus; no single-cell positions exist in this dataset",
        "structures": structures,
        "regions": mapping,
        "populations": pop_rows,
        "evidence": {
            "literature": lit,
            "regions": regions_lit,
            "causal": causal,
            "kg": {"available": False, "path": "src/syn4rch/kg/", "module_files": kg_files,
                   "status": "knowledge graph not built; no brain-to-body links are drawn"},
        },
    }

    with open(os.path.join(out_dir, "allen_mapping.json"), "w", encoding="utf-8") as fh:
        json.dump({"built_at": scene["built_at"], "ontology_url": ONTOLOGY_URL,
                   "rules": ontology.__doc__.strip(), "regions": mapping,
                   "unresolved": [m["region_id"] for m in mapping if m["status"] != "resolved"]},
                  fh, indent=1)
    with open(os.path.join(out_dir, "scene.json"), "w", encoding="utf-8") as fh:
        json.dump(scene, fh, indent=1)
    with open(os.path.join(out_dir, "geometry.json"), "w", encoding="utf-8") as fh:
        fh.write(_dump(geometry))

    # Measured signalling: curated population -> systemic readout counts. Every
    # number in it is recomputed from records.jsonl on each build.
    sig = signaling.build_graph(
        os.path.join(root, "data", "curated", "records.jsonl"),
        os.path.join(root, "data", "vocab", "populations.json"),
        os.path.join(root, "data", "vocab", "outcomes.json"))
    signaling.write_graph(sig, os.path.join(out_dir, "signaling.json"))

    tpl_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "viewer_template.html")
    with open(tpl_path, encoding="utf-8") as fh:
        tpl = fh.read()
    html = (tpl.replace("__SCENE_JSON__", inline_json(scene))
               .replace("__GEOMETRY_JSON__", inline_json(geometry))
               .replace("__SIGNALING_JSON__", inline_json(sig)))
    page = os.path.join(root, "viewer", "index.html")
    with open(page, "w", encoding="utf-8") as fh:
        fh.write(html)
    return {"page": page, "page_bytes": os.path.getsize(page), "mapping": mapping,
            "literature": lit, "signaling": sig}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python3 -m syn4rch.atlas.build", description=__doc__.split("\n\n")[0])
    ap.add_argument("--offline", action="store_true", help="use cached downloads only")
    args = ap.parse_args(argv)
    try:
        res = build(offline=args.offline)
    except HostUnreachable as exc:
        print("STOP: Allen host unreachable, nothing substituted: %s" % exc, file=sys.stderr)
        return 2
    for m in res["mapping"]:
        print("  %-7s %-10s %s" % (m["region_id"], m["status"],
              "%s %s (%s)" % (m["allen_id"], m["allen_acronym"], m["allen_name"]) if m["allen_id"] else m["reason"]))
    print("literature source: %s (%s)" % (res["literature"]["source"], res["literature"]["run_timestamp"]))
    sig = res["signaling"]
    print("signalling: %d of %d records used, %d edges (%d contested, %d unclear axis), "
          "%d populations with records, %d readouts"
          % (sig["n_records_used"], sig["n_records_in"], sig["n_edges"],
             sig["n_edges_contested"], sig["n_edges_unclear_axis"],
             sig["n_populations_with_records"], len(sig["outcomes"])))
    print("wrote %s (%.2f MB)" % (res["page"], res["page_bytes"] / 1e6))
    return 0


if __name__ == "__main__":
    sys.exit(main())
