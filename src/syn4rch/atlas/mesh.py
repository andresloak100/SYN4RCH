"""Wavefront OBJ parsing and vertex-clustering simplification (stdlib only)."""

from __future__ import annotations

import math
from typing import Dict, Iterable, List, Sequence, Tuple

Vec = Tuple[float, float, float]
Tri = Tuple[int, int, int]


def parse_obj(lines: Iterable[str]) -> Tuple[List[Vec], List[Tri]]:
    """Return vertices and triangles (0-based). Polygons are fan-triangulated."""
    verts: List[Vec] = []
    tris: List[Tri] = []
    for line in lines:
        if line.startswith("v "):
            p = line.split()
            verts.append((float(p[1]), float(p[2]), float(p[3])))
        elif line.startswith("f "):
            idx = []
            for tok in line.split()[1:]:
                i = int(tok.split("/")[0])
                idx.append(i - 1 if i > 0 else len(verts) + i)
            for k in range(1, len(idx) - 1):
                tris.append((idx[0], idx[k], idx[k + 1]))
    return verts, tris


def bbox(verts: Sequence[Vec]) -> Tuple[Vec, Vec]:
    xs, ys, zs = zip(*verts)
    return (min(xs), min(ys), min(zs)), (max(xs), max(ys), max(zs))


def cluster(verts: Sequence[Vec], tris: Sequence[Tri], cell: float) -> Tuple[List[Vec], List[Tri]]:
    """Vertex clustering on a uniform grid of edge ``cell``.

    Every vertex is replaced by the mean of the vertices sharing its grid
    cell; triangles that collapse or duplicate are dropped. Geometry is only
    ever moved by less than one cell diagonal, never invented.
    """
    if cell <= 0:
        raise ValueError("cell must be positive")
    (x0, y0, z0), _ = bbox(verts)
    key_of: Dict[Tuple[int, int, int], int] = {}
    sums: List[List[float]] = []
    remap: List[int] = []
    for (x, y, z) in verts:
        k = (int((x - x0) // cell), int((y - y0) // cell), int((z - z0) // cell))
        j = key_of.get(k)
        if j is None:
            j = len(sums)
            key_of[k] = j
            sums.append([0.0, 0.0, 0.0, 0.0])
        s = sums[j]
        s[0] += x; s[1] += y; s[2] += z; s[3] += 1
        remap.append(j)
    new_verts = [(s[0] / s[3], s[1] / s[3], s[2] / s[3]) for s in sums]
    seen = set()
    new_tris: List[Tri] = []
    for a, b, c in tris:
        a, b, c = remap[a], remap[b], remap[c]
        if a == b or b == c or a == c:
            continue
        key = tuple(sorted((a, b, c)))
        if key in seen:
            continue
        seen.add(key)
        new_tris.append((a, b, c))
    # drop vertices no longer referenced
    used = sorted({i for t in new_tris for i in t})
    pos = {old: new for new, old in enumerate(used)}
    return [new_verts[i] for i in used], [(pos[a], pos[b], pos[c]) for a, b, c in new_tris]


def simplify_to(verts: Sequence[Vec], tris: Sequence[Tri], max_verts: int,
                min_cell: float = 10.0, max_cell: float = 2000.0) -> Tuple[List[Vec], List[Tri], float]:
    """Smallest cell (>= ``min_cell``) whose clustered mesh has <= ``max_verts``.

    Returns (verts, tris, cell). If the mesh is already small enough at
    ``min_cell`` that result is returned.
    """
    lo, hi = min_cell, max_cell
    v, t = cluster(verts, tris, lo)
    if len(v) <= max_verts:
        return v, t, lo
    best = None
    for _ in range(18):
        mid = math.sqrt(lo * hi)
        v, t = cluster(verts, tris, mid)
        if len(v) <= max_verts:
            best = (v, t, mid)
            hi = mid
        else:
            lo = mid
        if hi / lo < 1.03:
            break
    if best is None:
        v, t = cluster(verts, tris, max_cell)
        best = (v, t, max_cell)
    return best


def encode(verts: Sequence[Vec], tris: Sequence[Tri]) -> Dict:
    """Compact JSON-ready geometry: integer CCF micrometres + flat indices."""
    lo, hi = bbox(verts)
    n = len(verts)
    cx = sum(v[0] for v in verts) / n
    cy = sum(v[1] for v in verts) / n
    cz = sum(v[2] for v in verts) / n
    return {
        "units": "um",
        "frame": "CCFv3 (x=AP, y=DV, z=ML)",
        "vertex_count": n,
        "triangle_count": len(tris),
        "bbox_um": [[round(c) for c in lo], [round(c) for c in hi]],
        "vertex_mean_um": [round(cx), round(cy), round(cz)],
        "p": [round(c) for v in verts for c in v],
        "i": [i for t in tris for i in t],
    }
