"""Resolve SYN4RCH region vocabulary entries to Allen structure ids.

Matching rules (exact, never fuzzy):

* Search scope is the Allen Hypothalamus subtree (acronym ``HY``) only, so a
  synonym such as ``ARC`` cannot accidentally hit an unrelated structure
  elsewhere in the brain.
* Tier 1: the entry's ``abbreviation`` equals an Allen acronym (case
  sensitive), or its ``label`` equals an Allen name (normalised).
* Tier 2: any ``synonym`` equals an Allen acronym (case sensitive) or name
  (normalised).
* One distinct tier-1 id wins; other ids hit only by synonyms are reported as
  ``alternates`` so the reviewer can see the ambiguity. With no tier-1 hit,
  exactly one distinct tier-2 id resolves. Anything else is unresolved.
"""

from __future__ import annotations

import re
import unicodedata
from typing import Dict, Iterable, List, Optional

SCOPE_ACRONYM = "HY"


def norm(text: str) -> str:
    t = unicodedata.normalize("NFKD", text)
    t = "".join(c for c in t if not unicodedata.combining(c)).lower()
    t = re.sub(r"[-_/,]+", " ", t)
    return re.sub(r"\s+", " ", t).strip()


def flatten(graph: Dict) -> List[Dict]:
    """Flatten the Allen structure_graph_download JSON into a list of nodes.

    Each node gets ``path`` (list of ancestor acronyms, root first).
    """
    roots = graph["msg"] if "msg" in graph else [graph]
    out: List[Dict] = []

    def walk(node: Dict, path: List[str]) -> None:
        out.append({
            "id": node["id"],
            "acronym": node["acronym"],
            "name": node["name"],
            "parent_id": node.get("parent_structure_id"),
            "color": node.get("color_hex_triplet"),
            "path": list(path),
        })
        for child in node.get("children") or []:
            walk(child, path + [node["acronym"]])

    for r in roots:
        walk(r, [])
    return out


def subtree(nodes: Iterable[Dict], acronym: str) -> List[Dict]:
    nodes = list(nodes)
    return [n for n in nodes if n["acronym"] == acronym or acronym in n["path"]]


def find_by_id(nodes: Iterable[Dict], sid: int) -> Optional[Dict]:
    for n in nodes:
        if n["id"] == sid:
            return n
    return None


def resolve_entry(entry: Dict, scope: List[Dict]) -> Dict:
    by_acr: Dict[str, Dict] = {n["acronym"]: n for n in scope}
    by_name: Dict[str, List[Dict]] = {}
    for n in scope:
        by_name.setdefault(norm(n["name"]), []).append(n)

    tier1: Dict[int, List[str]] = {}
    tier2: Dict[int, List[str]] = {}

    def hit(bucket: Dict[int, List[str]], node: Dict, why: str) -> None:
        bucket.setdefault(node["id"], []).append(why)

    abbr = entry.get("abbreviation")
    if abbr and abbr in by_acr:
        hit(tier1, by_acr[abbr], "abbreviation=acronym:" + abbr)
    label = entry.get("label")
    if label:
        for n in by_name.get(norm(label), []):
            hit(tier1, n, "label=name")
    for syn in entry.get("synonyms") or []:
        if syn in by_acr:
            hit(tier2, by_acr[syn], "synonym=acronym:" + syn)
        for n in by_name.get(norm(syn), []):
            hit(tier2, n, "synonym=name:" + syn)

    result = {
        "region_id": entry["id"],
        "label": label,
        "status": "unresolved",
        "allen_id": None,
        "allen_acronym": None,
        "allen_name": None,
        "matched_by": [],
        "alternates": [],
        "reason": None,
    }
    if len(tier1) == 1:
        chosen = next(iter(tier1))
    elif len(tier1) > 1:
        result["status"] = "ambiguous"
        result["reason"] = "abbreviation and label match different Allen structures"
        result["alternates"] = _describe(tier1, scope)
        return result
    elif len(tier2) == 1:
        chosen = next(iter(tier2))
    elif len(tier2) > 1:
        result["status"] = "ambiguous"
        result["reason"] = "synonyms match more than one Allen structure and no abbreviation/label match decides"
        result["alternates"] = _describe(tier2, scope)
        return result
    else:
        result["reason"] = "no exact acronym or name match in the Allen %s subtree" % SCOPE_ACRONYM
        return result

    node = find_by_id(scope, chosen)
    result.update({
        "status": "resolved",
        "allen_id": node["id"],
        "allen_acronym": node["acronym"],
        "allen_name": node["name"],
        "matched_by": tier1.get(chosen, []) + tier2.get(chosen, []),
        "alternates": _describe({k: v for k, v in tier2.items() if k != chosen}, scope),
    })
    return result


def _describe(bucket: Dict[int, List[str]], scope: List[Dict]) -> List[Dict]:
    out = []
    for sid, why in sorted(bucket.items()):
        n = find_by_id(scope, sid)
        out.append({"allen_id": sid, "allen_acronym": n["acronym"],
                    "allen_name": n["name"], "matched_by": why})
    return out


def resolve_regions(regions: Dict[str, Dict], graph: Dict) -> List[Dict]:
    """Resolve every region vocabulary entry. ``regions`` is the ``entries`` map."""
    nodes = flatten(graph)
    scope = subtree(nodes, SCOPE_ACRONYM)
    if not scope:
        raise ValueError("Allen graph has no %s structure" % SCOPE_ACRONYM)
    return [resolve_entry(e, scope) for _, e in regions.items()]
