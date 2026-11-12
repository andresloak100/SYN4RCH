import copy
import json

import pytest

from analysis import circuit_evidence_viewer as V


@pytest.fixture
def g():
    return V.fixture_graph()


# ---------------- graph checks (species, layers, claims)

def test_fixture_is_valid_and_clearly_synthetic(g):
    assert V.validate_graph(g)
    assert all(n['label'].startswith('FIXTURE') for n in g['nodes'])
    assert all(e['review_status'] == 'fixture' for e in g['edges'])


@pytest.mark.parametrize('mutate', [
    lambda g: g.update(schema='other'),
    lambda g: g['edges'][0].update(target='missing_node'),
    lambda g: g['nodes'].append(dict(g['nodes'][0])),                                    # duplicate id
    lambda g: g['edges'][0].pop('layer'),
    lambda g: g['nodes'][5].update(species='  '),                                         # blank, not null
    lambda g: g['nodes'][0].pop('species'),
    lambda g: g['edges'][0].update(uncertainty='low'),
    lambda g: g['edges'][0].update(confidence=0.9),                                       # invented confidence
    lambda g: g['edges'][3]['conditions'].update(calibrated_probability=0.7),
    lambda g: g['nodes'][0].update(edge_probability=1),
    lambda g: g['nodes'][0]['atlas_placement'].update(individual_cell_position=True),
    lambda g: g['nodes'][0]['atlas_placement'].update(level='cell'),
    lambda g: g['nodes'][4].update(atlas_placement={'level': 'region_only', 'region_id': 'DMH', 'allen_id': 830,
                                                    'individual_cell_position': False}),   # rat placed in mouse
    lambda g: g['nodes'][3].update(atlas_placement={'level': 'region_only', 'region_id': 'DMH', 'allen_id': 830,
                                                    'individual_cell_position': False}),   # fly placed in mouse
])
def test_graph_refusals(g, mutate):
    mutate(g)
    with pytest.raises(V.Refusal):
        V.validate_graph(g)


def test_fixture_uses_the_actual_schema(g):
    from analysis import circuit_evidence_layer as L
    assert L.validate_graph(g) is not False                                   # ROOT's validator accepts the fixture
    assert {e['layer'] for e in g['edges']} <= set(L.EDGE_LAYERS.values())
    assert all(L.EDGE_LAYERS[e['kind']] == e['layer'] for e in g['edges'])


@pytest.mark.parametrize('sp,mouse,fly', [
    ('mus_musculus', True, False), ('Mus musculus', True, False), ('mouse', True, False),
    ('rattus_norvegicus', False, False), (None, False, False), ('drosophila_melanogaster', False, True)])
def test_species_canonical_forms(sp, mouse, fly):
    assert V.is_mouse(sp) is mouse and V.is_fly(sp) is fly


def test_string_curation_confidence_is_allowed_numeric_is_not(g):
    g['edges'][0]['curation_confidence'] = 'high'                            # source curation enum: allowed
    g['edges'][0]['conditions']['score_units'] = 'points'                     # ordinary source field: allowed
    assert V.validate_graph(g)
    g['edges'][0]['confidence'] = 0.8
    with pytest.raises(V.Refusal):
        V.validate_graph(g)


def test_scene_availability_by_region_and_species(g):
    regions = V.scene_regions(V.read_scene(V.SOURCE_VIEWER.read_text()))
    assert regions['DMH'] == 830 and regions['LHA'] == 194 and 'VTA' not in regions
    out = {n['id']: n['_scene'] for n in V.attach_scene_regions(g, regions)['nodes']}
    assert out['fixture_dmh'] == {'available': True, 'acronym': 'DMH', 'allen_id': 830}   # mus_musculus is mouse
    assert out['fixture_vta']['available'] is False and out['fixture_vta']['allen_id'] is None
    assert out['fixture_lhb']['available'] is False                            # no mesh, nothing invented
    for nid, sp in (('fixture_fly', 'drosophila_melanogaster'), ('fixture_rat', 'rattus_norvegicus'),
                    ('fixture_unknown', None)):
        assert out[nid] == {'available': False, 'not_mouse_anatomy': True, 'species': sp}


def test_context_meshes_pinned_and_highlightable(g):
    ctx = V.load_context_meshes()
    assert [(c['acr'], c['id']) for c in ctx] == [('VTA', 749), ('PAG', 795)]
    assert all('not hypothalamic' in c['label'] for c in ctx)
    regions = V.scene_regions(V.read_scene(V.SOURCE_VIEWER.read_text()), ctx)
    assert regions['VTA'] == 749 and regions['PAG'] == 795 and regions['DMH'] == 830
    out = {n['id']: n['_scene'] for n in V.attach_scene_regions(g, regions)['nodes']}
    assert out['fixture_vta'] == {'available': True, 'acronym': 'VTA', 'allen_id': 749}
    assert out['fixture_lhb']['available'] is False


def test_context_mesh_with_wrong_pin_is_refused(tmp_path):
    m = json.loads(V.CONTEXT_MANIFEST.read_text())
    m['items'][0]['allen_id'] = 830
    p = tmp_path / 'm.json'
    p.write_text(json.dumps(m))
    with pytest.raises(V.Refusal):
        V.load_context_meshes(p)                                                 # manifest hash differs
    with pytest.raises(V.Refusal):
        V.load_context_meshes(p, expected_sha=None)                              # and the per-item pin differs


def test_region_bindings_are_pinned_and_display_only():
    b = V.load_region_bindings()
    assert {k: (v['acronym'], v['allen_id']) for k, v in b.items()} == \
        {'hou2026:vta': ('VTA', 749), 'tokizane2024:pag': ('PAG', 795)}


@pytest.mark.parametrize('change', [
    lambda m: m['evidence_node_region_bindings'][0].update(allen_id=830),        # not a pinned context mesh
    lambda m: m['evidence_node_region_bindings'][0].update(acronym='PAG'),       # id/acronym disagree
    lambda m: m.update(ontology_sha256='0' * 64),
])
def test_region_binding_mismatch_is_refused(tmp_path, change):
    m = json.loads(V.CONTEXT_MANIFEST.read_text())
    change(m)
    p = tmp_path / 'm.json'
    p.write_text(json.dumps(m))
    with pytest.raises(V.Refusal):
        V.load_region_bindings(p, expected_sha=None)


def test_bindings_only_apply_to_mouse_region_nodes_without_placement(g):
    regions = V.scene_regions(V.read_scene(V.SOURCE_VIEWER.read_text()), V.load_context_meshes())
    g['nodes'][1].pop('atlas_placement')                                          # fixture_vta: region, no placement
    bind = {'fixture_vta': {'acronym': 'VTA', 'allen_id': 749, 'basis': 'test'}}
    before = copy.deepcopy(g)
    out = {n['id']: n for n in V.attach_scene_regions(g, regions, bind)['nodes']}
    assert g == before                                                           # graph facts unchanged
    sc = out['fixture_vta']['_scene']
    assert sc['available'] and sc['acronym'] == 'VTA' and sc['display_binding'] == 'atlas context annotation'
    assert 'atlas_placement' not in out['fixture_vta']
    for bad in ('fixture_dmh', 'fixture_rat', 'fixture_lhb'):                    # placed / not mouse / population
        with pytest.raises(V.Refusal):
            V.attach_scene_regions(g, regions, {bad: bind['fixture_vta']})


def test_region_lookup_by_allen_id_only(g):
    regions = {'DMH': 830}
    g['nodes'][0]['atlas_placement']['region_id'] = None
    out = V.attach_scene_regions(g, regions)['nodes'][0]['_scene']
    assert out['available'] and out['acronym'] == 'DMH'


# ---------------- rendering

def test_render_embeds_graph_exactly_and_safely(g, tmp_path):
    g['nodes'][0]['label'] = 'FIXTURE </script><b>x</b>'                       # must not break the page
    out = V.render(g, {'note': 'c'}, out=tmp_path / 'v.html', fixture=True)
    html = out.read_text()
    data = V.read_embedded(html)
    stripped = copy.deepcopy(data['graph'])
    for n in stripped['nodes']:
        n.pop('_scene')
    assert stripped == g and data['fixture'] is True and data['coverage'] == {'note': 'c'}
    assert html.count('</script>') == 1                                        # only the template's own close tag
    assert html.index('const EVIDENCE=') < html.index('</script>')


def test_render_refuses_invalid_graph(g, tmp_path):
    g['edges'][0]['probability'] = 0.5
    with pytest.raises(V.Refusal):
        V.render(g, {}, out=tmp_path / 'v.html')
    assert not (tmp_path / 'v.html').exists()


def test_template_states_scope_and_has_one_placeholder():
    t = V.TEMPLATE.read_text()
    assert t.count('/*__EVIDENCE_DATA__*/null') == 1
    assert 'filter(Boolean).sort((a,b)=>isMouse' not in t                     # unknown species is not dropped
    for phrase in ("replace(/_/g,' ')", 'unknown species', "includes('drosophila')", 'extra context regions',
                   'No confidence scores are computed', 'Placement is by atlas region only',
                   'never matched to individual traces', 'not mouse anatomy', 'not available in this 3D scene',
                   'not reported', 'must not be shared'):
        assert phrase in t, phrase


# ---------------- measured neuron geometry unchanged

def test_derived_scene_keeps_geometry_and_original_untouched(tmp_path):
    before = V.sha256(V.SOURCE_VIEWER)
    info = V.build_neuron_view(dst=tmp_path / 'n.html')
    assert V.sha256(V.SOURCE_VIEWER) == before == V.SOURCE_VIEWER_SHA256
    src, der = V.SOURCE_VIEWER.read_text(), (tmp_path / 'n.html').read_text()
    assert V.read_scene(der) == V.read_scene(src)                             # neurons, meshes, licence text
    assert der.startswith(src.rstrip()[:-len('</body></html>')])              # original bytes kept verbatim
    bridge = der[len(src.rstrip()) - len('</body></html>'):]
    assert 'synarch-highlight-region' in bridge and 'never labels or recolours traces' in bridge
    for forbidden in ('x.grp', 'N.forEach', '.soma', 'D.neurons'):
        assert forbidden not in bridge                                         # bridge never touches traces
    assert info['regions']['DMH'] == 830


def test_derived_scene_with_context_keeps_original_bytes(tmp_path):
    ctx = V.load_context_meshes()
    info = V.build_neuron_view(dst=tmp_path / 'n.html', context=ctx)
    src, der = V.SOURCE_VIEWER.read_text(), (tmp_path / 'n.html').read_text()
    assert V.sha256(V.SOURCE_VIEWER) == V.SOURCE_VIEWER_SHA256
    assert V.read_scene(der) == V.read_scene(src)                             # neurons, meshes, frame unchanged
    assert der.startswith(src.rstrip()[:-len('</body></html>')])
    bridge = der[len(src.rstrip()) - len('</body></html>'):]
    assert 'extra context regions: ' in bridge and 'not hypothalamic' in bridge and '"acr": "VTA"' in bridge
    for forbidden in ('x.grp', 'N.forEach', '.soma', 'D.neurons', 'Line(', 'LineSegments'):
        assert forbidden not in bridge                                         # no traces touched, no lines drawn
    assert info['regions']['VTA'] == 749 and info['regions']['PAG'] == 795
    assert 'syn-compact' in bridge and "getElementById('side')" in bridge         # controls toggled, not removed
    assert 'remove()' not in bridge and 'innerHTML=\'\'' not in bridge


def test_rendered_page_keeps_unknown_species(g, tmp_path):
    data = V.read_embedded(V.render(g, {}, out=tmp_path / 'v.html').read_text())
    assert sum(e['species'] is None for e in data['graph']['edges']) == 1


def test_derived_scene_refuses_changed_source(tmp_path):
    p = tmp_path / 'swc_viewer.html'
    p.write_text(V.SOURCE_VIEWER.read_text() + ' ')
    with pytest.raises(V.Refusal):
        V.build_neuron_view(src=p, dst=tmp_path / 'n.html')


# ---------------- real graph (skipped until ROOT's build exists)

def test_template_lists_regions_and_labels_inferred_types():
    t = V.TEMPLATE.read_text()
    for phrase in ("Populations and regions", "n.kind==='population'||n.kind==='region'", 'author_inferred_type',
                   'author-inferred cell type, not a measured identity', 'source_outcome_raw',
                   'outcome as written in the source', 'not synapses and not independent observations',
                   'Identified synapses', 'Matched individual cells', 'ppp1r17', 'a gap, not a negative result',
                   'where in the paper', 'reported_statistics', 'not reported (unknown)', 'n_scope', 'what n counts',
                   'cohort age (this item only)', 'display_binding', 'not cell coordinates', 'max-width:900px',
                   'd-filters', 'id="facts"'):
        assert phrase in t, phrase


@pytest.mark.skipif(not V.GRAPH.exists(), reason='real graph not built')
def test_real_graph_renders_with_every_record_and_raw_outcome(tmp_path):
    g = json.loads(V.GRAPH.read_text())
    cov = json.loads(V.COVERAGE.read_text())
    ctx = V.load_context_meshes()
    bind = V.load_region_bindings()
    data = V.read_embedded(V.render(g, cov, out=tmp_path / 'v.html', context=ctx, bindings=bind).read_text())
    assert data['fixture'] is False and len(data['graph']['edges']) == len(g['edges'])
    assert data['graph']['edges'] == g['edges']                                  # every link kept, unchanged
    recs = {e['record_id'] for e in g['edges'] if e.get('record_id') and e.get('source_collection')}
    assert recs == {e['record_id'] for e in data['graph']['edges'] if e.get('record_id') and e.get('source_collection')}
    nodes = {n['id']: n for n in data['graph']['nodes']}
    # region nodes that receive interventions are present (listed under "Populations and regions")
    hit = {e['target'] for e in g['edges'] if e['kind'] == 'intervention_target' and nodes[e['target']]['kind'] == 'region'}
    assert hit and all(nodes[i]['kind'] == 'region' for i in hit)
    for n in nodes.values():                                                    # no non-mouse placement
        if not V.is_mouse(n['species']):
            assert n['_scene']['available'] is False
        if n['_scene'].get('available') and not n['_scene'].get('display_binding'):
            assert n['atlas_placement']['level'] == 'region_only'
        if n['_scene'].get('display_binding'):
            assert n['kind'] == 'region' and n['id'] in bind
    assert cov['edges_by_kind'].get('identified_synapse', 0) == 0
    for nid, acr in (('hou2026:vta', 'VTA'), ('tokizane2024:pag', 'PAG')):
        assert nodes[nid]['_scene']['acronym'] == acr and nodes[nid]['atlas_placement'] is None
    stripped = [{k: v for k, v in n.items() if k != '_scene'} for n in data['graph']['nodes']]
    assert stripped == g['nodes']                                                # biological facts unchanged
    studies = {s['study_id']: s for s in g['studies']}
    assert all(str(studies[e['study_id']].get('url', '')).startswith('https://') for e in g['edges'])
