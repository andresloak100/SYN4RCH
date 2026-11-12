"""Local 3D circuit-evidence viewer: source-backed evidence graph beside the existing measured neuron scene.

* Reads data/processed/circuit_evidence_v1/graph.json (+ coverage.json), schema synarch_circuit_evidence_graph_v1,
  written by the ROOT-owned evidence layer. This module only checks and displays it; it adds no biology.
* The 3D scene is a DERIVED copy of reports/swc_viewer.html (never modified) with a small appended message API
  that highlights existing Allen region meshes. Neuron geometry, hemispheres, licence and provenance text are
  byte-identical to the original; traces are never labelled with evidence populations.
* Region placement only: populations link to Allen regions present as meshes; regions without a mesh (e.g. VTA)
  are shown as "not available in this 3D scene". Fly data never appear as mouse anatomy.
* No confidence scores or calibrated probabilities are invented; missing fields are shown as "not reported".
* LOCAL ONLY: the scene embeds traces whose licence restricts sharing (see the original viewer).

Usage:
  /opt/anaconda3/bin/python -I analysis/circuit_evidence_viewer.py            (real graph, when it exists)
  /opt/anaconda3/bin/python -I analysis/circuit_evidence_viewer.py --fixture  (clearly labelled synthetic fixture)
"""
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE_VIEWER = ROOT / 'reports/swc_viewer.html'
SOURCE_VIEWER_SHA256 = 'dafd46ddae4d9739cd17d43cc2f4aaf6bce4c04aa1d4fcfd31fc59d1250aa7ca'
TEMPLATE = ROOT / 'src/syn4rch/atlas/circuit_evidence_template.html'
GRAPH = ROOT / 'data/processed/circuit_evidence_v1/graph.json'
COVERAGE = ROOT / 'data/processed/circuit_evidence_v1/coverage.json'
OUT = ROOT / 'reports/circuit_evidence_viewer.html'
NEURON_OUT = ROOT / 'reports/circuit_evidence_neurons.html'
SCHEMA = 'synarch_circuit_evidence_graph_v1'
MOUSE = {'mus musculus', 'mouse'}                 # compared after canonicalisation (underscores -> spaces)
# Same keys as circuit_evidence_layer.reject_numeric_confidence, plus any '*probability' key. Only NUMERIC values are
# refused there; string curation enums (e.g. curation_confidence='high') and other source fields are allowed.
NUMERIC_CONFIDENCE_KEYS = {'confidence', 'confidence_score', 'confidence_probability', 'edge_probability', 'probability'}
CONTEXT_MANIFEST = ROOT / 'data/interim/circuit_evidence_v1/atlas_context_manifest.json'
CONTEXT_PINS = {'VTA': (749, '9dd98e09e05b155b11606af662e1455eaf4e68a9959b714ba6d8abdcc37b0ba6'),
                'PAG': (795, 'e819c2b07bcde880ec68962802365bce9ba437e43d9538347af643b3527ba8a3')}
CCF_EXTENT_UM = (13200.0, 8000.0, 11400.0)          # AP, DV, ML
CONTEXT_MANIFEST_SHA256 = '1f536d41b84196957ba43e623374878894c1f46d7e1da7a8783ad4b1dc97356c'   # ROOT, with 2 bindings
ONTOLOGY = ROOT / 'data/raw/atlas/structure_graph_1.json'
ONTOLOGY_SHA256 = 'b5f0d024d1df09ee18aef15093f55240f4a70b3924ee9689add612dadff80386'
BINDING_LABEL = 'atlas context annotation'


class Refusal(ValueError):
    pass


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------- the existing measured scene

def read_scene(html):
    """Parse the embedded `const D=...;` JSON of the measured viewer."""
    a = html.index('const D=') + len('const D=')
    b = html.index('\n', a)
    return json.loads(html[a:b].rstrip().rstrip(';'))


def scene_regions(scene, context=None):
    """Allen regions present as meshes in the derived scene: original meshes plus pinned context meshes."""
    out = {m['acr']: int(m['id']) for m in scene['meshes'] if m['acr'] != 'root'}
    for c in context or []:
        out[c['acr']] = int(c['id'])
    return out


def read_obj(path):
    """Vertices (AP, DV, ML micrometres, Allen CCFv3 2017) and 0-based triangle indices of an Allen OBJ mesh."""
    import numpy as np
    V, F = [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith('v '):
            V.append([float(x) for x in line.split()[1:4]])
        elif line.startswith('f '):
            idx = [int(t.split('/')[0]) - 1 for t in line.split()[1:]]
            if len(idx) != 3:
                raise Refusal('only triangular faces are supported')
            F.append(idx)
    V, F = np.asarray(V, np.float32), np.asarray(F, np.int32)
    if V.ndim != 2 or V.shape[1] != 3 or F.size == 0 or F.min() < 0 or F.max() >= len(V):
        raise Refusal('malformed OBJ mesh')
    return V, F


def _pinned_manifest(manifest, expected_sha):
    if expected_sha is not None and sha256(manifest) != expected_sha:
        raise Refusal('context manifest differs from its pinned SHA-256')
    return json.loads(Path(manifest).read_text())


def load_context_meshes(manifest=CONTEXT_MANIFEST, expected_sha=CONTEXT_MANIFEST_SHA256):
    """The ROOT-pinned extra context regions (VTA, PAG): hash-checked, frame-checked, regional only."""
    import base64
    import numpy as np
    m = _pinned_manifest(manifest, expected_sha)
    out = []
    for item in m['items']:
        acr = item['acronym']
        if acr not in CONTEXT_PINS or CONTEXT_PINS[acr] != (int(item['allen_id']), item['sha256']):
            raise Refusal(f'context mesh {acr} is not one of the pinned regions')
        path = ROOT / item['path']
        if sha256(path) != item['sha256']:
            raise Refusal(f'context mesh {acr} differs from its pinned SHA-256')
        V, F = read_obj(path)
        if (V < 0).any() or (V > np.asarray(CCF_EXTENT_UM, np.float32)).any():
            raise Refusal(f'context mesh {acr} is outside the CCFv3 extent; frame not as expected')
        out.append({'acr': acr, 'id': int(item['allen_id']), 'nv': len(V), 'nf': len(F),
                    'V': base64.b64encode(V.tobytes()).decode(), 'F': base64.b64encode(F.astype(np.int32).tobytes()).decode(),
                    'sha256': item['sha256'], 'label': 'extra context region (not hypothalamic); regional mesh, not cell locations'})
    return out


BRIDGE = r"""
<script>
/* SYN4RCH evidence-viewer bridge (appended to a derived copy; the measured scene above is unchanged).
   It only shows or highlights EXISTING Allen region meshes. It never labels or recolours traces. */
(function(){
 /* Extra context region meshes (pinned Allen CCFv3 2017 OBJ; same AP,DV,ML frame as D.meshes). Added to the
    derived scene only; D itself is untouched. */
 const CTX=/*__CONTEXT_MESHES__*/[];
 for(const m of CTX){const V=dec(m.V,Float32Array),F=dec(m.F,Int32Array);const P=new Float32Array(V.length);
  for(let i=0;i<V.length;i+=3){P[i]=V[i+2];P[i+1]=V[i+1];P[i+2]=V[i]}
  const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.BufferAttribute(P,3));g.setIndex(new THREE.BufferAttribute(new Uint32Array(F.buffer),1));g.computeVertexNormals();
  const o=new THREE.Mesh(g,new THREE.MeshLambertMaterial({transparent:true,opacity:0.16,depthWrite:false,side:THREE.DoubleSide}));
  o.renderOrder=-1;world.add(o);meshObjs[m.acr]=o}
 if(CTX.length){const l=document.createElement('label');l.innerHTML='<input type="checkbox" id="mCtx" checked> extra context regions: '+CTX.map(m=>m.acr).join(', ')+' (not hypothalamic; regional meshes only)';
  const anchor=document.getElementById('mNuc').parentNode;anchor.parentNode.insertBefore(l,anchor.nextSibling);
  const cb=l.querySelector('input');const vis=()=>{for(const m of CTX)meshObjs[m.acr].visible=cb.checked;draw()};cb.addEventListener('change',vis)}
 const REG={};for(const m of D.meshes){if(m.acr!=='root')REG[m.acr]=m.id}for(const m of CTX)REG[m.acr]=m.id;
 let hl=null;const OP={};for(const k in meshObjs)OP[k]=meshObjs[k].material.opacity;
 function restore(){for(const k in OP)meshObjs[k].material.opacity=OP[k]}
 const note=document.createElement('div');note.className='note';
 note.textContent='Opened from the circuit-evidence viewer: region highlights mark atlas regions only. Traces are not identified with any evidence population.';
 document.getElementById('side').insertBefore(note,document.getElementById('side').children[2]);
 function paint(){if(!hl)return;const o=meshObjs[hl];if(!o)return;o.visible=true;o.material.opacity=0.45;
  o.material.color.set(isDark()?0xffb347:0xc0561a);draw()}
 function focus(acr){const o=meshObjs[acr];const b=new THREE.Box3().setFromObject(o);b.getCenter(target);
  const s=b.getSize(new THREE.Vector3()).length();rad=Math.min(80000,Math.max(1500,s*2.6));place();draw()}
 function reply(msg){try{parent.postMessage(msg,'*')}catch(e){}}
 addEventListener('message',ev=>{const m=ev.data||{};
  if(m.type==='synarch-highlight-region'){
   const acr=m.acronym&&REG[m.acronym]!==undefined?m.acronym:Object.keys(REG).find(k=>REG[k]===Number(m.allen_id));
   if(!acr){hl=null;restore();apply();reply({type:'synarch-region-status',request:m.request,available:false});return}
   hl=null;restore();apply();hl=acr;paint();if(m.focus!==false)focus(acr);
   reply({type:'synarch-region-status',request:m.request,available:true,acronym:acr,allen_id:REG[acr]})}
  else if(m.type==='synarch-clear-region'){hl=null;restore();apply();reply({type:'synarch-region-status',request:m.request,cleared:true})}
 });
 document.querySelectorAll('#side input').forEach(e=>e.addEventListener('change',()=>setTimeout(paint,0)));
 document.getElementById('theme').addEventListener('click',()=>setTimeout(paint,0));
 /* Compact controls: at narrow widths the original page stacks its long control panel above the canvas. Here the
    canvas fills the frame and every original control stays reachable from one toggle (nothing removed). */
 const st=document.createElement('style');st.textContent=
  'body.syn-compact{flex-direction:row;overflow:hidden;height:100%}'
  +'body.syn-compact #view{height:100%;flex:1}'
  +'body.syn-compact #side{display:none}'
  +'body.syn-compact.syn-open #side{display:block;position:fixed;z-index:5;left:0;top:0;bottom:0;width:min(340px,92vw);min-width:0;height:auto;overflow:auto;border-right:1px solid var(--line);box-shadow:0 0 12px rgba(0,0,0,.25)}'
  +'#synCtl{position:fixed;z-index:6;right:8px;bottom:8px;font:inherit;font-size:12px;padding:4px 10px;background:var(--panel);color:var(--ink);border:1px solid var(--line);border-radius:4px;cursor:pointer}';
 document.head.appendChild(st);
 const ctl=document.createElement('button');ctl.id='synCtl';ctl.type='button';document.body.appendChild(ctl);
 function setCompact(c,open){document.body.classList.toggle('syn-compact',c);document.body.classList.toggle('syn-open',c&&open);
  ctl.textContent=!c?'Hide scene controls':open?'Close scene controls':'Scene controls';ctl.setAttribute('aria-expanded',String(!c||open));resize()}
 let compact=matchMedia('(max-width:760px)').matches,open=false;setCompact(compact,open);
 ctl.onclick=()=>{if(!compact){compact=true;open=false}else if(!open){open=true}else{open=false;if(!matchMedia('(max-width:760px)').matches)compact=false}setCompact(compact,open)};
 reply({type:'synarch-scene-ready',regions:REG});
})();
</script></body></html>"""


def _ontology_acronyms(path=ONTOLOGY, expected_sha=ONTOLOGY_SHA256):
    if sha256(path) != expected_sha:
        raise Refusal('Allen ontology differs from its pinned SHA-256')
    out, stack = {}, list(json.loads(Path(path).read_text())['msg'])
    while stack:
        n = stack.pop()
        out[int(n['id'])] = n['acronym']
        stack.extend(n.get('children') or [])
    return out


def load_region_bindings(manifest=CONTEXT_MANIFEST, expected_sha=CONTEXT_MANIFEST_SHA256, ontology=ONTOLOGY):
    """ROOT's explicit display bindings (evidence region node -> pinned context mesh). Display context only: they
    never change the graph, never add an atlas placement and never imply cell coordinates."""
    m = _pinned_manifest(manifest, expected_sha)
    if m.get('ontology_sha256') != ONTOLOGY_SHA256:
        raise Refusal('context manifest names a different ontology')
    acr = _ontology_acronyms(ontology)
    out = {}
    for b in m.get('evidence_node_region_bindings', []):
        pair = (b['acronym'], int(b['allen_id']))
        if b['acronym'] not in CONTEXT_PINS or CONTEXT_PINS[b['acronym']][0] != pair[1]:
            raise Refusal(f"binding {b['node_id']} is not to a pinned context mesh")
        if acr.get(pair[1]) != b['acronym']:
            raise Refusal(f"binding {b['node_id']}: Allen id {pair[1]} is not {b['acronym']} in the pinned ontology")
        if b['node_id'] in out:
            raise Refusal(f"duplicate binding for {b['node_id']}")
        out[b['node_id']] = {'acronym': b['acronym'], 'allen_id': pair[1], 'basis': b.get('basis'),
                             'source_region_label': b.get('source_region_label')}
    return out


def build_neuron_view(src=SOURCE_VIEWER, dst=NEURON_OUT, expected_sha=SOURCE_VIEWER_SHA256, context=None):
    """Write the derived scene: the original bytes with the bridge (and pinned context meshes) appended."""
    if sha256(src) != expected_sha:
        raise Refusal('source neuron viewer differs from its pinned SHA-256; refusing to derive from it')
    html = Path(src).read_text()
    tail = '</script></body></html>'
    if not html.rstrip().endswith(tail):
        raise Refusal('unexpected end of the source viewer')
    cut = html.rstrip()[:-len('</body></html>')]
    ctx = json.dumps(context or []).replace('</', '<\\/')
    bridge = BRIDGE.lstrip('\n').replace('/*__CONTEXT_MESHES__*/[]', ctx, 1)
    derived = cut + bridge.replace('</script></body></html>', '</script></body></html>\n', 1)
    if read_scene(derived) != read_scene(html):
        raise Refusal('derived scene data differ from the source')
    Path(dst).write_text(derived)
    return {'source_sha256': expected_sha, 'derived_sha256': sha256(dst),
            'regions': scene_regions(read_scene(html), context)}


# ---------------------------------------------------------------- the evidence graph

def _walk_items(obj, path=''):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield path + k, k, v
            yield from _walk_items(v, path + k + '.')
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _walk_items(v, f'{path}{i}.')


def canonical_species(species):
    return None if species is None else str(species).strip().lower().replace('_', ' ')


def is_mouse(species):
    return canonical_species(species) in MOUSE


def is_fly(species):
    return 'drosophila' in (canonical_species(species) or '')


def validate_graph(g):
    if g.get('schema') != SCHEMA:
        raise Refusal(f'graph schema must be {SCHEMA}')
    for key in ('nodes', 'edges'):
        if not isinstance(g.get(key), list):
            raise Refusal(f'graph.{key} must be a list')
    for path, k, v in _walk_items(g):
        kl = k.lower()
        numeric = isinstance(v, (int, float)) and not isinstance(v, bool)
        if numeric and (kl in NUMERIC_CONFIDENCE_KEYS or kl.endswith('_probability')):
            raise Refusal(f'numeric confidence/probability must not be invented: {path}')
    ids = [n.get('id') for n in g['nodes']]
    if len(set(ids)) != len(ids) or any(not i for i in ids):
        raise Refusal('node ids must be present and unique')
    for n in g['nodes']:
        if 'species' in n and n['species'] is not None and not str(n['species']).strip():
            raise Refusal(f"node {n['id']}: species must be explicit or null (unknown)")
        ap = n.get('atlas_placement')
        if ap:
            if ap.get('level') != 'region_only' or ap.get('individual_cell_position') not in (False, None):
                raise Refusal(f"node {n['id']}: atlas placement must be region_only without individual cell position")
            if not is_mouse(n.get('species')) and (ap.get('allen_id') or ap.get('region_id')):
                raise Refusal(f"node {n['id']}: non-mouse data cannot be placed in the mouse atlas")
    eids = [e.get('id') for e in g['edges']]
    if len(set(eids)) != len(eids) or any(not i for i in eids):
        raise Refusal('edge ids must be present and unique')
    known = set(ids)
    for e in g['edges']:
        if e.get('source') not in known or e.get('target') not in known:
            raise Refusal(f"edge {e['id']} references an unknown node")
        if not e.get('layer'):
            raise Refusal(f"edge {e['id']} needs an evidence layer")
        if not isinstance(e.get('uncertainty', []), list):
            raise Refusal(f"edge {e['id']}: uncertainty must be a list of stated limitations")
    return True


def attach_scene_regions(g, regions, bindings=None):
    """Mark, per mouse population, whether its atlas region exists as a mesh in the 3D scene.

    `bindings` (ROOT's explicit display bindings) apply only to mouse REGION nodes with no curated atlas placement;
    they are shown as an atlas context annotation, never as a placement or cell coordinates."""
    by_id = {v: k for k, v in regions.items()}
    out = json.loads(json.dumps(g))
    bindings = bindings or {}
    nodes = {n['id']: n for n in out['nodes']}
    for nid, b in bindings.items():
        n = nodes.get(nid)
        if n is None or n.get('kind') != 'region' or not is_mouse(n.get('species')) or n.get('atlas_placement'):
            raise Refusal(f'display binding {nid} must name a mouse region node without a curated atlas placement')
        if regions.get(b['acronym']) != b['allen_id']:
            raise Refusal(f'display binding {nid}: {b["acronym"]} mesh is not in the scene')
    for n in out['nodes']:
        if n['id'] in bindings:
            b = bindings[n['id']]
            n['_scene'] = {'available': True, 'acronym': b['acronym'], 'allen_id': b['allen_id'],
                           'display_binding': BINDING_LABEL, 'basis': b['basis']}
            continue
        ap = n.get('atlas_placement') or {}
        acr = ap.get('region_id') if ap.get('region_id') in regions else by_id.get(ap.get('allen_id'))
        if is_mouse(n.get('species')):
            n['_scene'] = {'available': bool(acr), 'acronym': acr, 'allen_id': regions.get(acr)} if ap \
                else {'available': False, 'no_region': True}
        else:
            n['_scene'] = {'available': False, 'not_mouse_anatomy': True, 'species': n.get('species')}
    return out


def render(g, coverage, template=TEMPLATE, out=OUT, neuron_src='circuit_evidence_neurons.html', fixture=False,
           context=None, bindings=None):
    validate_graph(g)
    regions = scene_regions(read_scene(Path(SOURCE_VIEWER).read_text()), context)
    payload = {'graph': attach_scene_regions(g, regions, bindings), 'coverage': coverage or g.get('coverage') or {},
               'scene_regions': regions, 'fixture': bool(fixture), 'neuron_viewer': neuron_src}
    data = json.dumps(payload, sort_keys=True).replace('</', '<\\/')
    html = Path(template).read_text()
    if html.count('/*__EVIDENCE_DATA__*/null') != 1:
        raise Refusal('template data placeholder missing or repeated')
    Path(out).write_text(html.replace('/*__EVIDENCE_DATA__*/null', data))
    return out


def read_embedded(html):
    m = re.search(r'const EVIDENCE=(.*?);\n', html, re.S)
    return json.loads(m.group(1).replace('<\\/', '</'))


# ---------------------------------------------------------------- synthetic fixture (no biology)

def fixture_graph():
    """Clearly labelled synthetic fixture in the ACTUAL graph schema (analysis/circuit_evidence_layer.py).

    It exercises species (mouse, fly, rat, unknown), every display path and region availability. No node or edge
    is a biological claim; labels start with FIXTURE and the source locator is not a study."""
    def pop(i, species, region=None, allen=None, kind='population'):
        n = {'id': f'fixture_{i}', 'kind': kind, 'label': f'FIXTURE {kind} {i}', 'species': species,
             'region_id': region, 'identity_status': 'fixture placeholder'}
        if kind == 'population':
            n.update(marker_genes=[f'FixtureGene{i}'], gene_measurement_scope='identity_annotation_only')
        if region:
            n['atlas_placement'] = {'level': 'region_only', 'region_id': region, 'allen_id': allen,
                                    'individual_cell_position': False}
        return n
    M, FLY, RAT = 'mus_musculus', 'drosophila_melanogaster', 'rattus_norvegicus'
    nodes = [pop('dmh', M, 'DMH', 830), pop('vta', M, 'VTA', 749, kind='region'), pop('lhb', M, 'LHb', 186),
             pop('fly', FLY), pop('rat', RAT), pop('unknown', None),
             {'id': 'fixture_intervention', 'kind': 'intervention', 'label': 'FIXTURE intervention', 'species': M},
             {'id': 'fixture_outcome', 'kind': 'outcome', 'label': 'FIXTURE outcome measure', 'species': M}]

    def edge(i, src, dst, kind, layer, species, level, **kw):
        e = {'id': f'fx_{i}', 'source': src, 'target': dst, 'kind': kind, 'layer': layer, 'species': species,
             'study_id': 'fixture_study', 'measurement_level': level,
             'source_locator': 'https://example.org/fixture-not-a-study', 'review_status': 'fixture',
             'is_anatomical_edge': kind in {'measured_axon_projection', 'identified_synapse'},
             'conditions': {'sex': None, 'n_animals': None}, 'direction': None, 'effect': None,
             'perturbation_class': None, 'perturbation_polarity': None,
             'uncertainty': ['synthetic fixture: no biological meaning']}
        e.update(kw)
        return e
    edges = [
        edge(1, 'fixture_dmh', 'fixture_vta', 'measured_axon_projection', 'anatomical', M, 'population',
             conditions={'sex': 'female', 'n_animals': 3}),
        edge(2, 'fixture_dmh', 'fixture_lhb', 'inferred_connection', 'inference', M, 'inferred'),
        edge(3, 'fixture_intervention', 'fixture_dmh', 'intervention_target', 'functional', M, 'animal',
             perturbation_class='circuit_activity', perturbation_polarity='increase'),
        edge(4, 'fixture_intervention', 'fixture_outcome', 'intervention_outcome', 'functional', M, 'animal',
             direction='increase', perturbation_class='circuit_activity', perturbation_polarity='increase',
             conditions={'sex': 'not_reported', 'n_animals': 8},
             effect={'effect_size_value': None, 'effect_size_type': None, 'effect_size_units': None,
                     'ci_low': None, 'ci_high': None, 'n_treatment': None, 'n_control': None}),
        edge(5, 'fixture_fly', 'fixture_fly', 'functional_association', 'functional', FLY, 'population'),
        edge(6, 'fixture_rat', 'fixture_rat', 'functional_association', 'functional', RAT, 'animal'),
        edge(7, 'fixture_unknown', 'fixture_unknown', 'functional_association', 'functional', None, 'animal'),
    ]
    return {'schema': SCHEMA, 'nodes': nodes, 'edges': edges,
            'studies': [{'study_id': 'fixture_study', 'title': 'FIXTURE (not a study)', 'species': [M, FLY, RAT],
                         'url': None}],
            'coverage': {'note': 'synthetic fixture'}, 'limitations': ['Synthetic fixture for interface testing only.']}


def main(argv):
    context = load_context_meshes() if CONTEXT_MANIFEST.exists() else []
    scene = build_neuron_view(context=context)
    if argv == ['--fixture']:
        g, cov, fx = fixture_graph(), {'note': 'synthetic fixture'}, True
    elif not argv:
        if not GRAPH.exists():
            raise SystemExit(f'{GRAPH} does not exist yet; use --fixture to preview the interface')
        g = json.loads(GRAPH.read_text())
        cov = json.loads(COVERAGE.read_text()) if COVERAGE.exists() else {}
        fx = False
    else:
        sys.exit(__doc__)
    bindings = load_region_bindings() if (context and not fx) else None
    out = render(g, cov, fixture=fx, context=context, bindings=bindings)
    print('wrote', out, sha256(out), '| scene', scene['derived_sha256'][:16], '| fixture' if fx else '')


if __name__ == '__main__':
    main(sys.argv[1:])
