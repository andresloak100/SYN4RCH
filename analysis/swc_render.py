"""Render the Jiao 2025 SWC sample: interactive three.js viewer + two static PNG views.

    python3 -I analysis/swc_render.py            # full private viewer (portal + Zenodo geometry): local only
    python3 -I analysis/swc_render.py --share    # also writes reports/swc_viewer_zenodo_only.html (18 CC BY 4.0 neurons only)

Inputs: data/raw/jiao2025/swc/ (raw SWCs, MANIFEST.tsv), data/interim/swc_validation.json
(written by analysis/swc_validate.py), data/raw/atlas/meshes/*.obj (CCFv3 um, AP DV ML).
Outputs: reports/swc_viewer.html (self-contained except three.js r128 from cdnjs),
data/interim/swc_view_dorsal.png, data/interim/swc_view_sagittal.png,
data/interim/swc_render.json (decimation statistics).
With --share, also reports/swc_viewer_zenodo_only.html: only the 18 Zenodo CC BY 4.0 Jiao files (each
drawn from its Zenodo file, also for the 4 neurons that have a portal file), the Allen meshes, and
attribution with a statement of the changes made for display. No portal geometry is embedded.

Scene orientation: data (AP, DV, ML) in CCFv3 um are right-handed (AP to posterior, DV to ventral,
ML to the animal's right; the portal index puts Left somas at ML < 5700). The viewer maps them to
scene x = -ML (toward the animal's LEFT), y = -DV (dorsal), z = -AP (anterior). That frame is also
right-handed (determinant +1), so the brain is not shown mirror-imaged.

Measured traces only; nothing is synthesised. One trace per neuron: the portal file when the
neuron has one (true hemisphere), otherwise the Zenodo file (mirrored-left). The viewer shows a
chosen frame by a display-time reflection ML' = 11400 - ML of whole neurons only; the raw files
are never changed.

Display decimation (viewer only; PNGs use every node): within each unbranched run between kept
nodes, Ramer-Douglas-Peucker with tolerance EPS um. Always kept: the soma, every branch point,
every tip, and every node whose type differs from its parent's or any child's (so type styling
is exact). The maximum distance of any dropped node from the drawn polyline is measured and
reported."""
import os, json, base64, colorsys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SWC = os.path.join(ROOT, "data/raw/jiao2025/swc")
MESH = os.path.join(ROOT, "data/raw/atlas/meshes")
OUT = os.path.join(ROOT, "reports/swc_viewer.html")
OUT_SHARE = os.path.join(ROOT, "reports/swc_viewer_zenodo_only.html")
# Subresource integrity for cdnjs three.js r128 three.min.js, as published by the cdnjs API
# (https://api.cdnjs.com/libraries/three.js/r128?fields=sri, read 2026-11-10).
SRI = "sha512-dLxUelApnYxpLt6K2iomGngnHO83iUvZytA3YjDUCjT0HDOHKXnVYdf3hU4JjM8uEhxf9nD1/ey98U3t2vZ0qQ=="
EPS = 1.0
MID = 5700.0
CLS = {3: 0, 2: 1, 0: 2}            # 0 dendrite (type 3), 1 axon (type 2), 2 type 0 (unlabelled); 3 other; 4 soma
CLS_NAMES = ["dendrite (type 3)", "axon (type 2)", "type 0 (unlabelled)", "other (types 10/11/12, undocumented)", "soma (type 1)"]
NUCLEI = ["PVH", "ARH", "VMH", "DMH", "LHA", "ZI", "AHN", "PH", "MPO", "LPO", "SCH", "SO", "ME", "PMv", "TM", "AVPV", "MEPO", "VMPO"]
LIC_Z = "CC BY 4.0, Jiao et al. 2025, Zenodo 14588094"
LIC_P = ("Portal (mouse.digital-brain.cn): copyright line only, 'Copyright 2023 Digital Brain. All Rights Reserved.'; "
         "BSDC record of the same data (doi:10.12412/BSDC.1667278800.20003): CC BY-NC-ND 4.0. Local private viewing with attribution only; "
         "this page embeds simplified and reflected traces, so it must not be shared")


def read_swc(p):
    a = np.array([l.split() for l in open(p, "rb").read().decode().splitlines() if l.strip() and not l.startswith("#")], float)
    return a[:, 0].astype(np.int64), a[:, 1].astype(np.int64), a[:, 2:5], a[:, 6].astype(np.int64)


def read_obj(p):
    V, F = [], []
    for line in open(p):
        if line.startswith("v "): V.append([float(x) for x in line.split()[1:4]])
        elif line.startswith("f "): F.append([int(x.split("/")[0]) - 1 for x in line.split()[1:4]])
    return np.array(V, np.float32), np.array(F, np.int32)


def seg_dist(P, A, B):
    d = B - A; L = float(d @ d)
    if L == 0: return np.linalg.norm(P - A, axis=1)
    t = np.clip((P - A) @ d / L, 0, 1)
    return np.linalg.norm(P - (A + t[:, None] * d), axis=1)


def rdp(P, eps):
    keep = np.zeros(len(P), bool); keep[0] = keep[-1] = True
    st = [(0, len(P) - 1)]
    while st:
        i, j = st.pop()
        if j <= i + 1: continue
        d = seg_dist(P[i + 1:j], P[i], P[j]); k = int(d.argmax())
        if d[k] > eps: keep[i + 1 + k] = True; st += [(i, i + 1 + k), (i + 1 + k, j)]
    return keep


def decimate(ids, types, X, parent):
    n = len(ids); pos = {int(i): k for k, i in enumerate(ids)}
    pidx = np.array([pos[int(p)] if p != -1 else -1 for p in parent])
    ch = [[] for _ in range(n)]
    for k, p in enumerate(pidx):
        if p >= 0: ch[p].append(k)
    must = np.zeros(n, bool)
    for k in range(n):
        if pidx[k] < 0 or len(ch[k]) != 1: must[k] = True
        if pidx[k] >= 0 and types[pidx[k]] != types[k]: must[k] = True; must[pidx[k]] = True
    keep = must.copy(); maxd = 0.0
    for k in np.where(must)[0]:
        for c in ch[k]:
            run = [k, c]
            while not must[run[-1]]: run.append(ch[run[-1]][0])
            kp = rdp(X[run], EPS)
            if not kp.all():
                P = X[run]; ki = np.where(kp)[0]
                for a, b in zip(ki[:-1], ki[1:]):
                    if b > a + 1: maxd = max(maxd, float(seg_dist(P[a + 1:b], P[a], P[b]).max()))
            keep[np.array(run)[kp]] = True
    kidx = np.where(keep)[0]; newi = -np.ones(n, int); newi[kidx] = np.arange(len(kidx))
    kpar = []
    for k in kidx:
        p = pidx[k]
        while p >= 0 and not keep[p]: p = pidx[p]
        kpar.append(newi[p] if p >= 0 else -1)
    cls = np.array([4 if t == 1 else CLS.get(int(t), 3) for t in types[kidx]], np.uint8)
    return X[kidx].astype(np.float32), np.array(kpar, np.int32), cls, maxd, pidx


def edge_types(types, pidx):
    """Type used to style the edge into each node: the node's own type, except that an edge into an
    undocumented marker node (types 10/11/12) takes the type of the nearest non-marker ancestor."""
    et = types.copy()
    for k in np.where((types != 0) & (types != 1) & (types != 2) & (types != 3))[0]:
        q = pidx[k]
        while q >= 0 and types[q] not in (0, 1, 2, 3): q = pidx[q]
        et[k] = types[q] if q >= 0 else types[k]
    return et


def b64(a): return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()


def colours(n):
    out = []
    for i in range(n):
        h = (i * 0.618033988749895) % 1.0
        out.append(h)
    return out


def main(share=False):
    val = json.load(open(os.path.join(ROOT, "data/interim/swc_validation.json")))
    byn = {}
    for r in val["neurons"]: byn.setdefault(r["neuron"], {})[r["frame"]] = r
    pairs = {p["neuron"]: p for p in val["pairs"]}
    if share:   # CC BY 4.0 Zenodo files only; never embed portal geometry in a shareable page
        byn = {n: {"zenodo_mirrored_left": d["zenodo_mirrored_left"]} for n, d in byn.items() if "zenodo_mirrored_left" in d}
    order = sorted(byn, key=lambda n: (list(byn[n].values())[0]["index"]["type_name"], n))
    hues = colours(len(order))
    neurons, stats, full = [], [], {}
    for i, n in enumerate(order):
        d = byn[n]
        r = d.get("portal_true_hemisphere") or d["zenodo_mirrored_left"]
        ids, types, X, parent = read_swc(os.path.join(SWC, r["file"]))
        P, par, cls, maxd, pidx = decimate(ids, types, X, parent)
        full[n] = (X, pidx, types, r)
        ism = ~np.isin(types, (0, 1, 2, 3)) & (pidx >= 0)
        mkl = np.linalg.norm(X[ism] - X[pidx[ism]], axis=1)
        stats.append({"neuron": n, "file": r["file"], "nodes": len(ids), "kept": len(P), "max_displacement_um": maxd,
                      "marker_edges": int(ism.sum()), "marker_edge_um": float(mkl.sum()), "_mk": mkl})
        files = [{"file": x["file"], "frame": x["frame"], "sha": x["sha256"], "bytes": x["bytes"], "source": x["source"],
                  "licence": LIC_Z if x["frame"].startswith("zenodo") else LIC_P} for x in d.values()]
        if share: files[0]["licence"] = LIC_Z + "; changed for display: simplified (RDP 1.0 um) and, for Right neurons in the true-hemisphere frame, reflected ML' = 11400 - ML"
        xl = r["xlsx"] or {}
        neurons.append({"id": n, "line": r["index"]["type_name"], "label": xl.get("location"), "index_acronym": None if share else r["index"]["acronym"],
                        "hemisphere": r["index"]["hemisphere"], "recon": None if share else r["index"]["reconstruction_type"], "frame": r["frame"],
                        "drawn_file": r["file"], "files": files, "hue": hues[i], "nodes": len(ids), "kept": len(P),
                        "maxd": round(maxd, 3), "types": r["type_counts"], "len": {k: round(v) for k, v in r["length_um"].items()},
                        "bp": r["branch_points"], "tips": r["tips"], "soma": r["soma"], "pair": None if share else pairs.get(n, {}).get("relation"),
                        "outside": r["nodes_outside_root"], "P": b64(P), "par": b64(par), "cls": b64(cls)})
    meshes = []
    for f in sorted(os.listdir(MESH)):
        V, F = read_obj(os.path.join(MESH, f))
        sid = int(f[:-4])
        acr = {997: "root", 1097: "HY"}.get(sid)
        if acr is None:
            acr = [m[0] for m in val["meshes"] if m[1] == sid][0]
        meshes.append({"acr": acr, "id": sid, "V": b64(V), "F": b64(F), "nv": len(V), "nf": len(F)})
    tot_nodes = sum(s["nodes"] for s in stats); tot_kept = sum(s["kept"] for s in stats)
    mk = np.concatenate([s.pop("_mk") for s in stats]) if stats else np.zeros(0)
    mk_n, mk_len, mk_max = len(mk), float(mk.sum()), float(mk.max()) if len(mk) else 0.0
    maxd_all = max(s["max_displacement_um"] for s in stats)
    simp = ("Display only: within each unbranched run, Ramer-Douglas-Peucker with tolerance %.1f um. Kept: soma, all branch points, all tips, "
            "all nodes where the type changes. %d of %d nodes drawn (%.1f%%); largest distance of a dropped node from the drawn line %.3f um. "
            "Meshes are drawn in full (no mesh decimation). Lines are 1 pixel wide; SWC radii are not used (0.004 to 56 um, not diameters). "
            "Somas are drawn as a fixed 28 um-radius marker sphere, not a measured soma size. "
            "Undocumented marker nodes (types 10/11/12, %d nodes) are drawn as points, and the %d edges into them (%.1f um in all, longest %.1f um) are drawn "
            "in the style of the nearest non-marker ancestor's type, so arbors have no gaps. "
            "Neurons marked with the reflection sign are mirrored across ML = 5,700 um for the chosen frame (ML' = 11400 - ML, whole neuron). "
            "Raw SWC files are unchanged." % (EPS, tot_kept, tot_nodes, 100.0 * tot_kept / tot_nodes, maxd_all, mk_n, mk_n, mk_len, mk_max))
    if not share:
        json.dump({"eps_um": EPS, "nodes": tot_nodes, "kept": tot_kept, "max_displacement_um": maxd_all,
                   "marker_edges": mk_n, "marker_edge_um": mk_len, "marker_edge_max_um": mk_max, "per_neuron": stats},
                  open(os.path.join(ROOT, "data/interim/swc_render.json"), "w"), indent=1)
    data = {"neurons": neurons, "meshes": meshes, "simp": simp, "cls_names": CLS_NAMES, "nuclei": NUCLEI, "share": share,
            "n_files": len(neurons) if share else len(val["neurons"]), "n_pairs": 0 if share else len(pairs)}
    html = TEMPLATE.replace("__SRI__", SRI).replace("__DATA__", json.dumps(data, separators=(",", ":")))
    out = OUT_SHARE if share else OUT
    open(out, "w").write(html)
    print("viewer", out, os.path.getsize(out), "bytes;", simp)
    if not share: pngs(full, order, hues)


def pngs(full, order, hues):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    root, _ = read_obj(os.path.join(MESH, "997.obj"))
    hy, _ = read_obj(os.path.join(MESH, "1097.obj"))
    views = [("dorsal", 2, 0, "ML (um, CCFv3)", "AP (um, CCFv3)", "Dorsal view (ML x AP), true-hemisphere frame"),
             ("sagittal", 0, 1, "AP (um, CCFv3)", "DV (um, CCFv3)", "Side view (AP x DV), true-hemisphere frame")]
    for name, a, b, xl, yl, title in views:
        fig, ax = plt.subplots(figsize=(10, 9 if name == "dorsal" else 6.2), dpi=150)
        ax.scatter(root[:, a], root[:, b], s=0.05, c="#b8b8b8", alpha=0.25, linewidths=0, rasterized=True)
        ax.scatter(hy[:, a], hy[:, b], s=0.2, c="#7a7a7a", alpha=0.35, linewidths=0, rasterized=True)
        if name == "dorsal": ax.axvline(MID, color="#555555", lw=0.8, ls="--")
        for i, n in enumerate(order):
            X, pidx, types, r = full[n]
            Y = X.copy()
            if r["frame"].startswith("zenodo") and r["index"]["hemisphere"] == "Right":
                Y[:, 2] = 11400 - Y[:, 2]       # display reflection back to the true hemisphere
            h = pidx >= 0
            segs = np.stack([Y[h][:, [a, b]], Y[pidx[h]][:, [a, b]]], 1)
            col = colorsys.hls_to_rgb(hues[i], 0.45, 0.85)
            isd = edge_types(types, pidx)[h] == 3
            ax.add_collection(LineCollection(segs[~isd], colors=[col], linewidths=0.25, alpha=0.75))
            ax.add_collection(LineCollection(segs[isd], colors=[col], linewidths=0.9, alpha=1.0))
            ax.plot(Y[0, a], Y[0, b], "o", ms=3, color=col, mec="black", mew=0.4)
        allx = np.concatenate([root[:, [a, b]]] + [full[n][0][:, [a, b]] for n in order])   # never clip traces
        ax.set_xlim(allx[:, 0].min() - 50, allx[:, 0].max() + 50); ax.set_ylim(allx[:, 1].min() - 50, allx[:, 1].max() + 50)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xlabel(xl); ax.set_ylabel(yl)
        ax.set_title("%s\n%d Jiao et al. 2025 neurons, measured traces only (every node drawn). Thick = dendrite (type 3); thin = type 2 or type 0.\n"
                     "Grey dots: CCFv3 root and HY mesh vertices. Zenodo-only Right neurons reflected back (ML' = 11400 - ML) for display." % (title, len(order)), fontsize=8)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        fig.tight_layout()
        p = os.path.join(ROOT, "data/interim/swc_view_%s.png" % name)
        fig.savefig(p); plt.close(fig); print("png", p, os.path.getsize(p))


TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Hypothalamic Neuron Traces</title>
<style>
:root{--bg:#f6f6f4;--panel:#ffffff;--ink:#1d1d1b;--muted:#5d5d58;--line:#d9d9d4;--accent:#1f5fa8;--warn:#8a4b00}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#121314;--panel:#1c1d1f;--ink:#ececea;--muted:#a3a39e;--line:#34353a;--accent:#7fb2ff;--warn:#f0b35a}}
:root[data-theme="dark"]{--bg:#121314;--panel:#1c1d1f;--ink:#ececea;--muted:#a3a39e;--line:#34353a;--accent:#7fb2ff;--warn:#f0b35a}
*{box-sizing:border-box}html,body{margin:0;height:100%}
body{background:var(--bg);color:var(--ink);font:13px/1.4 system-ui,-apple-system,Segoe UI,sans-serif;display:flex;overflow:hidden}
#side{width:340px;min-width:340px;height:100%;overflow:auto;background:var(--panel);border-right:1px solid var(--line);padding:12px 14px}
#view{flex:1;position:relative;min-width:0}#c{display:block;width:100%;height:100%}
h1{font-size:15px;margin:0 0 4px}h2{font-size:12px;margin:14px 0 4px;text-transform:uppercase;letter-spacing:.04em;color:var(--muted)}
.note{border:1px solid var(--warn);color:var(--warn);padding:6px 8px;border-radius:4px;font-weight:600;margin:6px 0}
.small{font-size:11.5px;color:var(--muted)}label{display:block;cursor:pointer}
.grp{margin:6px 0 2px;font-weight:600;display:flex;justify-content:space-between;align-items:center}
.grp button,.btn{font:inherit;font-size:11px;background:none;border:1px solid var(--line);color:var(--ink);border-radius:3px;padding:1px 6px;cursor:pointer}
.sw{display:inline-block;width:10px;height:10px;border-radius:2px;margin:0 5px 0 2px;vertical-align:-1px}
.nm{cursor:pointer;text-decoration:underline dotted}.nm.sel{color:var(--accent);font-weight:600}
#prov{position:absolute;right:10px;top:10px;width:330px;max-width:calc(100% - 20px);max-height:calc(100% - 20px);overflow:auto;background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:10px 12px;display:none}
#prov table{border-collapse:collapse;width:100%}#prov td{padding:2px 4px;vertical-align:top;border-bottom:1px solid var(--line);word-break:break-all}#prov td:first-child{color:var(--muted);white-space:nowrap;word-break:normal}
#ori{position:absolute;left:8px;top:8px;width:120px;height:120px;pointer-events:none}
#hud{position:absolute;left:10px;bottom:8px;font-size:11px;color:var(--muted);pointer-events:none}
@media (max-width:760px){body{flex-direction:column;overflow:auto}#side{width:100%;min-width:0;height:auto;border-right:0;border-bottom:1px solid var(--line);padding:12px 16px}#view{height:80vh;flex:none}}
</style></head><body>
<div id="side">
<h1>Hypothalamic neuron traces</h1>
<div class="small">Jiao et al. 2025, Nat Neurosci, doi:10.1038/s41593-025-01919-0. Allen CCFv3 (2017), micrometres.</div>
<div class="note">Measured traces only; no synthesised arbors.</div>
<div class="small" id="counts"></div>
<h2>Frame</h2>
<label><input type="radio" name="fr" value="true" checked> True hemisphere (portal frame; measured bilateral)</label>
<label><input type="radio" name="fr" value="mirror"> Mirrored-left (Zenodo / paper-table frame)</label>
<div class="small">Neurons drawn reflected (ML' = 11400 &minus; ML) in the chosen frame are marked &#8646; in the list. For the 4 neurons with both files the reflection was verified exact on every node.</div>
<h2>Show</h2>
<div id="cls"></div>
<label><input type="checkbox" id="mRoot" checked> brain outline (root mesh)</label>
<label><input type="checkbox" id="mHY" checked> hypothalamus (HY mesh)</label>
<label><input type="checkbox" id="mNuc"> 18 hypothalamic nuclei meshes</label>
<label><input type="checkbox" id="mMid" checked> midline plane (ML = 5,700 um)</label>
<label><input type="checkbox" id="soma" checked> soma markers</label>
<div style="margin-top:6px"><button class="btn" data-v="dorsal">dorsal</button> <button class="btn" data-v="lateral">left lateral</button> <button class="btn" data-v="front">frontal (anterior)</button> <button class="btn" id="theme">light/dark</button></div>
<h2>Neurons (by peptide line)</h2>
<div><button class="btn" id="all">all</button> <button class="btn" id="none">none</button></div>
<div id="list"></div>
<h2>Display simplifications</h2>
<div class="small" id="simp"></div>
<div class="small" style="margin-top:6px">Styling: dendrite (type 3) solid and opaque; axon (type 2) translucent; type 0 (unlabelled) translucent and dashed; other (types 10/11/12, undocumented, single nodes inside arbors) as points, with the edge into each drawn in its parent's style. Type 0 is not relabelled as axon. Soma spheres are fixed-size markers. The marker at top left shows A/P, L/R and D/V. Colours identify neurons; click a name for its provenance.</div>
<h2>Attribution and licences</h2>
<div class="small" id="lic"></div>
</div>
<div id="view"><canvas id="c"></canvas><canvas id="ori" width="120" height="120"></canvas><div id="prov"></div><div id="hud">drag: rotate &middot; right-drag or shift-drag: pan &middot; wheel or pinch: zoom</div></div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js" integrity="__SRI__" crossorigin="anonymous" referrerpolicy="no-referrer"></script>
<script>
const D=__DATA__;
function dec(s,T){const b=atob(s),u=new Uint8Array(b.length);for(let i=0;i<b.length;i++)u[i]=b.charCodeAt(i);return new T(u.buffer)}
const root=document.documentElement;
function isDark(){const t=root.getAttribute('data-theme');return t?t==='dark':matchMedia('(prefers-color-scheme: dark)').matches}
const canvas=document.getElementById('c');
const R=new THREE.WebGLRenderer({canvas,antialias:true});R.setPixelRatio(Math.min(2,devicePixelRatio));
const S=new THREE.Scene();
const cam=new THREE.PerspectiveCamera(35,1,10,200000);
// data (AP,DV,ML) are right-handed: AP->posterior, DV->ventral, ML->animal's right. Inside 'world' the
// coordinates are (ML,DV,AP); the world scale (-1,-1,-1) gives scene x = animal's LEFT, y = dorsal,
// z = anterior. The column swap (det -1) times this scale (det -1) has det +1: no mirror image.
const world=new THREE.Group();world.scale.set(-1,-1,-1);S.add(world);
function toW(ap,dv,ml){return [ml,dv,ap]}
const MIDX=5700;
// meshes
const meshObjs={};
for(const m of D.meshes){const V=dec(m.V,Float32Array),F=dec(m.F,Int32Array);const P=new Float32Array(V.length);
 for(let i=0;i<V.length;i+=3){P[i]=V[i+2];P[i+1]=V[i+1];P[i+2]=V[i]}
 const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.BufferAttribute(P,3));g.setIndex(new THREE.BufferAttribute(new Uint32Array(F.buffer),1));g.computeVertexNormals();
 const mat=new THREE.MeshLambertMaterial({transparent:true,opacity:m.acr==='root'?0.035:(m.acr==='HY'?0.07:0.22),depthWrite:false,side:THREE.DoubleSide});
 const o=new THREE.Mesh(g,mat);o.renderOrder=m.acr==='root'?-2:-1;world.add(o);meshObjs[m.acr]=o}
const amb=new THREE.AmbientLight(0xffffff,0.75);S.add(amb);const dl=new THREE.DirectionalLight(0xffffff,0.5);dl.position.set(1,2,1);S.add(dl);
// midline plane
const pg=new THREE.PlaneGeometry(13200,8000);const pm=new THREE.MeshBasicMaterial({transparent:true,opacity:0.08,side:THREE.DoubleSide,depthWrite:false});
const plane=new THREE.Mesh(pg,pm);plane.rotation.y=Math.PI/2;plane.position.set(MIDX,4000,6600);world.add(plane);
// neurons
const N=[];
for(const n of D.neurons){const P=dec(n.P,Float32Array),par=dec(n.par,Int32Array),cls=dec(n.cls,Uint8Array);
 const grp=new THREE.Group();const segs=[[],[],[],[]];const pts=[];
 for(let i=0;i<par.length;i++){const p=par[i];if(cls[i]===3)pts.push(P[3*i+2],P[3*i+1],P[3*i]);if(p<0)continue;let c=cls[i]===4?0:cls[i];
  // edge into an undocumented marker node: style of the nearest non-marker ancestor (no gap in the arbor)
  if(c===3){let q=p;while(q>=0&&cls[q]===3)q=par[q];c=q<0||cls[q]===4?2:cls[q]}
  segs[c].push(P[3*p+2],P[3*p+1],P[3*p],P[3*i+2],P[3*i+1],P[3*i])}
 const parts=[];
 for(let c=0;c<3;c++){if(!segs[c].length){parts.push(null);continue}const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.Float32BufferAttribute(segs[c],3));
  const mat=c===2?new THREE.LineDashedMaterial({dashSize:12,gapSize:6,transparent:true,opacity:0.6}):new THREE.LineBasicMaterial({transparent:c!==0,opacity:c===0?1:0.55});
  const l=new THREE.LineSegments(g,mat);if(c===2)l.computeLineDistances();grp.add(l);parts.push(l)}
 let op=null;if(pts.length){const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.Float32BufferAttribute(pts,3));op=new THREE.Points(g,new THREE.PointsMaterial({size:3,sizeAttenuation:false}));grp.add(op)}
 parts.push(op);
 const sm=new THREE.Mesh(new THREE.SphereGeometry(28,12,8),new THREE.MeshLambertMaterial());sm.position.set(n.soma[2],n.soma[1],n.soma[0]);grp.add(sm);
 world.add(grp);N.push({n,grp,parts,sm,on:true})}
function colOf(n){const c=new THREE.Color();c.setHSL(n.hue,0.85,isDark()?0.62:0.42);return c}
function frame(){return document.querySelector('input[name=fr]:checked').value}
function reflected(n){const f=frame();const right=n.hemisphere==='Right';
 if(n.frame==='portal_true_hemisphere')return f==='mirror'&&right;return f==='true'&&right}
function apply(){const dark=isDark();R.setClearColor(dark?0x121314:0xf6f6f4);
 for(const o of Object.values(meshObjs))o.material.color.set(dark?0xbfc4cc:0x505660);
 meshObjs.root.material.opacity=dark?0.045:0.035;
 pm.color.set(dark?0x7fb2ff:0x1f5fa8);
 const show=[...document.querySelectorAll('#cls input')].map(e=>e.checked);
 for(const x of N){const c=colOf(x.n);x.grp.visible=x.on;
  x.parts.forEach((p,i)=>{if(!p)return;p.material.color.copy(c);p.visible=show[i]});x.sm.material.color.copy(c);x.sm.visible=document.getElementById('soma').checked;
  const r=reflected(x.n);x.grp.scale.x=r?-1:1;x.grp.position.x=r?11400:0;
  const sw=document.getElementById('sw_'+x.n.id);if(sw){sw.style.background='#'+c.getHexString();document.getElementById('rf_'+x.n.id).textContent=r?' ⇆':''}}
 meshObjs.root.visible=document.getElementById('mRoot').checked;meshObjs.HY.visible=document.getElementById('mHY').checked;
 for(const a of D.nuclei)if(meshObjs[a])meshObjs[a].visible=document.getElementById('mNuc').checked;
 plane.visible=document.getElementById('mMid').checked;draw()}
// UI
const counts=D.neurons.reduce((a,n)=>(a[n.hemisphere]=(a[n.hemisphere]||0)+1,a),{});
document.getElementById('counts').textContent=`${D.neurons.length} neurons (${counts.Left||0} Left, ${counts.Right||0} Right soma) from ${D.n_files} raw files; ${new Set(D.neurons.map(n=>n.line)).size} peptide lines.`;
document.getElementById('simp').textContent=D.simp;
document.getElementById('lic').textContent=D.share
 ?'Traces: the 18 Jiao et al. 2025 SWC files of Zenodo record 14588094 (Supplementary Fig. 4), CC BY 4.0, doi:10.1038/s41593-025-01919-0. Changes made for this page: line simplification (Ramer-Douglas-Peucker, 1.0 um tolerance) and, in the true-hemisphere frame, reflection of Right-hemisphere neurons across ML = 5,700 um (ML\' = 11400 - ML); hemisphere per neuron from the public portal index. No portal (all rights reserved) or BSDC (CC BY-NC-ND 4.0) geometry is included. Meshes: Allen Mouse Brain Common Coordinate Framework v3 (2017), Allen Institute for Brain Science, used with citation under the Allen Institute terms of use. No MouseLight neurons are shown.'
 :'Zenodo 14588094 files: CC BY 4.0, Jiao et al. 2025. Portal files (mouse.digital-brain.cn, CEBSIT / Digital Brain): copyright line only, "Copyright 2023 Digital Brain. All Rights Reserved."; the BSDC record of the same data (doi:10.12412/BSDC.1667278800.20003) shows CC BY-NC-ND 4.0. This page embeds simplified and reflected portal traces: keep it local for private viewing with attribution; do not share or publish it (use the Zenodo-only build for sharing). Meshes: Allen Mouse Brain Atlas CCFv3 (2017), Allen Institute for Brain Science. No MouseLight neurons are shown.';
const clsBox=document.getElementById('cls');
D.cls_names.slice(0,4).forEach((t,i)=>{const l=document.createElement('label');l.innerHTML=`<input type="checkbox" checked> ${t}`;clsBox.appendChild(l)});
const lines={};for(const x of N)(lines[x.n.line]=lines[x.n.line]||[]).push(x);
const list=document.getElementById('list');
for(const [ln,xs] of Object.entries(lines).sort()){const g=document.createElement('div');
 g.innerHTML=`<div class="grp"><label><input type="checkbox" checked class="gl"> ${ln} (${xs.length})</label></div>`;
 const gl=g.querySelector('.gl');
 for(const x of xs){const l=document.createElement('label');const n=x.n;
  l.innerHTML=`<input type="checkbox" checked><span class="sw" id="sw_${n.id}"></span><span class="nm" id="nm_${n.id}">${n.id}</span> <span class="small">${n.label||''} ${n.hemisphere[0]} &middot; ${n.frame.startsWith('portal')?'portal':'zenodo'}<span id="rf_${n.id}"></span></span>`;
  const cb=l.querySelector('input');cb.onchange=()=>{x.on=cb.checked;apply()};x.cb=cb;
  l.querySelector('.nm').onclick=e=>{e.preventDefault();showProv(n)};g.appendChild(l)}
 gl.onchange=()=>{for(const x of xs){x.on=gl.checked;x.cb.checked=gl.checked}apply()};list.appendChild(g)}
document.getElementById('all').onclick=()=>{N.forEach(x=>{x.on=true;x.cb.checked=true});document.querySelectorAll('.gl').forEach(e=>e.checked=true);apply()};
document.getElementById('none').onclick=()=>{N.forEach(x=>{x.on=false;x.cb.checked=false});document.querySelectorAll('.gl').forEach(e=>e.checked=false);apply()};
document.querySelectorAll('#side input').forEach(e=>{if(!e.onchange)e.addEventListener('change',apply)});
const prov=document.getElementById('prov');
function esc(s){return String(s).replace(/[&<>"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]))}
function showProv(n){document.querySelectorAll('.nm').forEach(e=>e.classList.toggle('sel',e.id==='nm_'+n.id));
 const fr=n.files.map(f=>`<tr><td>file</td><td>${esc(f.file)}<br><span class="small">${esc(f.frame)} &middot; ${f.bytes.toLocaleString()} bytes<br>sha256 ${f.sha.slice(0,16)}&hellip;<br>${esc(f.source)}<br>licence: ${esc(f.licence)}</span></td></tr>`).join('');
 const L=n.len;prov.innerHTML=`<div style="display:flex;justify-content:space-between"><b>${n.id}</b><button class="btn" onclick="prov.style.display='none'">close</button></div>
 <table><tr><td>peptide line</td><td>${esc(n.line)}</td></tr><tr><td>soma label (xlsx)</td><td>${esc(n.label)} <span class="small">(portal index: ${esc(n.index_acronym)})</span></td></tr>
 <tr><td>hemisphere</td><td>${n.hemisphere} (portal index)</td></tr><tr><td>reconstruction</td><td>${esc(n.recon)}</td></tr>
 <tr><td>drawn from</td><td>${esc(n.drawn_file)}${reflected(n)?' <b>(drawn reflected for this frame)</b>':''}</td></tr>${fr}
 <tr><td>soma (AP,DV,ML)</td><td>${n.soma.map(v=>v.toFixed(1)).join(', ')} um (file frame)</td></tr>
 <tr><td>node types</td><td>${Object.entries(n.types).map(([k,v])=>k+': '+v).join(', ')}</td></tr>
 <tr><td>length (um)</td><td>type 2 ${L.axon.toLocaleString()}, type 0 ${L.type0.toLocaleString()}, type 3 ${L.dendrite.toLocaleString()}</td></tr>
 <tr><td>branch pts / tips</td><td>${n.bp} / ${n.tips}</td></tr><tr><td>drawn nodes</td><td>${n.kept.toLocaleString()} of ${n.nodes.toLocaleString()} (max displacement ${n.maxd} um)</td></tr>
 <tr><td>nodes outside root mesh</td><td>${n.outside}</td></tr>${n.pair?`<tr><td>portal vs Zenodo</td><td>${esc(n.pair)} (all nodes)</td></tr>`:''}</table>`;prov.style.display='block'}
// camera control (inline orbit)
S.updateMatrixWorld(true);const target=new THREE.Vector3();let rad=24000,th=0.6,ph=1.1;const BR=new THREE.Box3().setFromObject(meshObjs.HY);BR.getCenter(target);
function place(){cam.position.set(target.x+rad*Math.sin(ph)*Math.sin(th),target.y+rad*Math.cos(ph),target.z+rad*Math.sin(ph)*Math.cos(th));cam.up.set(0,1,0);cam.lookAt(target)}
let drag=null;canvas.addEventListener('contextmenu',e=>e.preventDefault());
canvas.addEventListener('pointerdown',e=>{canvas.setPointerCapture(e.pointerId);drag={x:e.clientX,y:e.clientY,pan:e.button===2||e.shiftKey}});
canvas.addEventListener('pointerup',()=>drag=null);
canvas.addEventListener('pointermove',e=>{if(!drag||pinch)return;const dx=e.clientX-drag.x,dy=e.clientY-drag.y;drag.x=e.clientX;drag.y=e.clientY;
 if(drag.pan){const s=rad*0.0012;const r=new THREE.Vector3().setFromMatrixColumn(cam.matrix,0),u=new THREE.Vector3().setFromMatrixColumn(cam.matrix,1);target.addScaledVector(r,-dx*s).addScaledVector(u,dy*s)}
 else{th-=dx*0.006;ph=Math.min(Math.PI-0.01,Math.max(0.01,ph-dy*0.006))}place();draw()});
canvas.addEventListener('wheel',e=>{e.preventDefault();rad*=Math.exp(e.deltaY*0.001);rad=Math.min(80000,Math.max(300,rad));place();draw()},{passive:false});
let pinch=null;canvas.addEventListener('touchstart',e=>{if(e.touches.length===2){const[a,b]=e.touches;pinch=Math.hypot(a.clientX-b.clientX,a.clientY-b.clientY)}},{passive:true});
canvas.addEventListener('touchmove',e=>{if(e.touches.length===2&&pinch){const[a,b]=e.touches;const d=Math.hypot(a.clientX-b.clientX,a.clientY-b.clientY);rad*=pinch/d;pinch=d;place();draw()}},{passive:true});
canvas.addEventListener('touchend',()=>{pinch=null});
document.querySelectorAll('[data-v]').forEach(b=>b.onclick=()=>{const v=b.dataset.v;// dorsal: from above, anterior at the top, animal's left on the left; left lateral: from the
 // animal's left (+x), anterior to the left; frontal: from in front (+z), animal's left on the right.
 if(v==='dorsal'){th=Math.PI;ph=0.02}else if(v==='lateral'){th=Math.PI/2;ph=Math.PI/2}else{th=0;ph=Math.PI/2}place();draw()});
document.getElementById('theme').onclick=()=>{root.setAttribute('data-theme',isDark()?'light':'dark');apply()};
matchMedia('(prefers-color-scheme: dark)').addEventListener('change',apply);
function resize(){const v=document.getElementById('view');R.setSize(v.clientWidth,v.clientHeight,false);cam.aspect=v.clientWidth/Math.max(1,v.clientHeight);cam.updateProjectionMatrix();draw()}
// orientation marker: anatomical directions in scene coordinates, projected with the camera rotation
const ORI=[['A',0,0,1],['P',0,0,-1],['L',1,0,0],['R',-1,0,0],['D',0,1,0],['V',0,-1,0]];
const oc=document.getElementById('ori').getContext('2d');
function drawOri(){const dark=isDark(),ink=dark?'#ececea':'#1d1d1b',mut=dark?'#a3a39e':'#5d5d58';oc.setTransform(1,0,0,1,0,0);oc.clearRect(0,0,120,120);
 const q=cam.quaternion.clone().invert();const items=ORI.map(([t,x,y,z])=>{const v=new THREE.Vector3(x,y,z).applyQuaternion(q);return [t,v]}).sort((a,b)=>a[1].z-b[1].z);
 oc.font='bold 12px system-ui,sans-serif';oc.textAlign='center';oc.textBaseline='middle';
 for(const [t,v] of items){const X=60+v.x*38,Y=60-v.y*38,front=v.z>=-0.05;oc.strokeStyle=front?ink:mut;oc.globalAlpha=front?1:0.45;oc.lineWidth=1.5;
  oc.beginPath();oc.moveTo(60,60);oc.lineTo(X,Y);oc.stroke();oc.fillStyle=front?ink:mut;oc.fillText(t,60+v.x*50,60-v.y*50)}oc.globalAlpha=1}
let pend=false;function draw(){if(pend)return;pend=true;requestAnimationFrame(()=>{pend=false;R.render(S,cam);drawOri()})}
addEventListener('resize',resize);place();apply();resize();
</script></body></html>
"""

if __name__ == "__main__":
    import sys
    main()
    if "--share" in sys.argv[1:]: main(share=True)
