"""Fixed, review-bound SYNTHETIC ONLY conditional-template benchmark.

No biological data/state is read. Import performs no fitting or data I/O.
Test truth is generated only after every world's test forecast is hash-saved.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import importlib.util
import json
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/interim/brm_template_synthetic'
PROTOCOL='reports/brm_template_protocol.md'
WORLDS=('informative','no_information','misleading')
SPLITS={'train':16,'validation':12,'test':12}
SEED=20261110
METHODS=('baseline','uniform_anchor_500','conditional_anchor_500','uniform_anchor_1000','conditional_anchor_1000','selected')
SPEC={'worlds':list(WORLDS),'rows_per_split':SPLITS,'neurons_per_animal':1,'regional_K':2,
      'seed':SEED,'truth_attachment_radius_um':750.,'origin_variation_bound_um':20.,
      'origins':'informative paired_identical_A_B_positions; no_information_and_misleading_all_zero',
      'base_nodes_relative_um':[[0.,0.,0.],[0.,-400.,0.],[1600.,-800.,0.],[2400.,-1600.,800.],[2400.,-1600.,-800.]],
      'base_parent':[-1,0,1,2,2],'family_B':'reflect_AP_about_zero','base_spacing_max_um':50.,
      'profile_A':[4.]*24+[.25]*25,'profile_B':[.25]*24+[4.]*25,'uninformative_profile':[1.]*49,
      'misleading':'training family=label; validation/test family=opposite label',
      'all_query_profiles_and_families_withheld':True,
      'informative_pass':'selected conditional; correct family12/12; lower distance AND higher F1 than baseline AND matched uniform',
      'no_information_pass':'baseline selected; selected arrays bitwise equal baseline for12/12',
      'misleading_pass':'baseline selected; selected arrays bitwise equal baseline for12/12',
      'exact_donor_claim':'known correct FAMILY among multiple equivalent donors; unique-donor unit fixture is separate',
      'test_truth_generation':'only after global216_forecast_hash_barrier'}


def library():
    name='brm_template_benchmark_library'
    if name not in sys.modules:
        spec=importlib.util.spec_from_file_location(name,ROOT/'analysis/brm_template.py')
        module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return sys.modules[name]


def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n'
def digest(value):return hashlib.sha256(canonical(value).encode()).hexdigest()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def seed(*parts):return int.from_bytes(hashlib.sha256('|'.join(map(str,(SEED,)+parts)).encode()).digest()[:4],'big')


def save(path,value,*,preserve=False):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);body=canonical(value)
    if path.exists():
        if path.read_text()==body:return
        if not preserve:raise RuntimeError(f'Preserve previous artifact: {path}')
        archive=path.parent/'archives'/(path.name+'.'+sha(path));archive.parent.mkdir(exist_ok=True)
        if not archive.exists():archive.write_bytes(path.read_bytes())
    path.write_text(body)


def save_tree(path,arbor):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('xb') as handle:np.savez_compressed(handle,xyz_um=arbor.xyz,parent=arbor.parent,types=arbor.types)


def known_rows(world,split):
    """Fixture index/world are never copied into the predictor Query."""
    count=SPLITS[split];rng=np.random.default_rng(seed(world,split,'origins'))
    origins=rng.uniform(-20.,20.,(count//2,3)) if world=='informative' else np.zeros((count//2,3))
    prefix={'informative':'INF','no_information':'NOI','misleading':'MIS'}[world]+{'train':'T','validation':'V','test':'E'}[split]
    rows=[]
    for i in range(count):
        animal=prefix+f'{i:02d}';neuron=animal+'_001';tag='U' if world=='no_information' else ('A' if i%2==0 else 'B')
        row={'world':world,'split':split,'index':i,'neuron':neuron,'animal':animal,
             'full_id':neuron+'_'+tag+'_SYN','soma_um':origins[i//2].tolist(),'tag':tag}
        if split=='train':row['profile']=list(SPEC['uninformative_profile'] if world=='no_information' else SPEC['profile_A' if tag=='A' else 'profile_B'])
        rows.append(row)
    return rows


def observed_dataset():return {w:{s:known_rows(w,s) for s in SPLITS} for w in WORLDS}


def subdivide(points,parent,types):
    """Independent segment subdivision; every original vertex survives."""
    xyz=[points[0].copy()];parents=[-1];outtypes=[int(types[0])];remap={0:0}
    for i in range(1,len(points)):
        start,end=points[parent[i]],points[i];n=max(1,int(np.ceil(np.linalg.norm(end-start)/50.-1e-12)))
        previous=remap[int(parent[i])]
        for k in range(1,n):
            xyz.append(start+(end-start)*(k/n));parents.append(previous);outtypes.append(2);previous=len(xyz)-1
        xyz.append(end.copy());parents.append(previous);outtypes.append(int(types[i]));remap[i]=len(xyz)-1
    return np.asarray(xyz),np.asarray(parents,dtype=int),np.asarray(outtypes,dtype=int)


def true_family(row):
    family=row['index']%2
    return 1-family if row['world']=='misleading' and row['split']!='train' else family


def make_truth(row):
    """Analytic generator only, independent of the candidate attachment code."""
    T=library();points=np.array(SPEC['base_nodes_relative_um']);family=true_family(row)
    if family:points[:,0]*=-1
    parent=np.array(SPEC['base_parent']);types=np.array([1,2,2,2,2])
    xyz,parent,types=subdivide(points,parent,types)
    depth=np.zeros(len(xyz))
    for i in range(1,len(xyz)):depth[i]=depth[parent[i]]+np.linalg.norm(xyz[i]-xyz[parent[i]])
    t=np.minimum(depth/750.,1.);u=1.-3*t*t+2*t*t*t;u[depth>=750.]=0.
    moved=xyz+u[:,None]*np.array(row['soma_um']);moved[0]=row['soma_um']
    assert np.array_equal(moved[depth>=750.],xyz[depth>=750.])
    xyz,parent,types=subdivide(moved,parent,types)
    return T.G.Arbor(xyz,parent,types).validate(),family


def query(row):
    T=library();soma=np.asarray(row['soma_um'])
    return T.Query(row['neuron'],row['animal'],row['full_id'],soma,soma,row['tag'])


def identity(row):return {key:row[key] for key in ('neuron','animal','full_id')}


def sources():
    T=library()
    paths=['analysis/brm_template.py','tests/test_brm_template.py','analysis/brm_template_benchmark.py',
           'tests/test_brm_template_benchmark.py',PROTOCOL,'analysis/brm_geometry.py','analysis/d110_geometry.py','analysis/brm_real.py','analysis/brm_joint.py']
    return {path:sha(ROOT/path) for path in paths}


def current_declaration():
    T=library();data=observed_dataset()
    return {'status':'declared_synthetic_only_unrun','specification':SPEC,'source_sha256':sources(),
        'library_provenance':T.source_provenance(),'runtime':T.Geo.R.runtime(),'synthetic_data_sha256':digest(data),
        'split_ids':{w:{s:[identity(r) for r in rows] for s,rows in splits.items()} for w,splits in data.items()},
        'seeds':{w:{**{s:seed(w,s,'origins') for s in SPLITS},'regional_fit':seed(w,'regional_fit')} for w in WORLDS}}


def unopened():
    if (BASE/'RUN_STARTED.json').exists():raise RuntimeError('Fixed benchmark already started; no automatic rerun')


def prepare():
    unopened()
    if (BASE/'frozen.json').exists():raise RuntimeError('Frozen metadata cannot change')
    d=current_declaration();save(BASE/'declaration.json',d,preserve=True);return d


def analytic_check():
    data=observed_dataset();animals=[]
    for world,splits in data.items():
        for split,rows in splits.items():
            assert len(rows)==SPLITS[split] and len({r['animal'] for r in rows})==len(rows)
            animals.extend(r['animal'] for r in rows)
            for i in range(0,len(rows),2):assert rows[i]['soma_um']==rows[i+1]['soma_um']
            assert sum(true_family(r)==0 for r in rows)==len(rows)//2
            if split!='train':assert all('profile' not in r and 'family' not in r and 'arbor' not in r for r in rows)
    assert len(animals)==len(set(animals))==120
    for world in WORLDS:
        for row in data[world]['train'][:2]:
            arbor,_=make_truth(row)
            assert np.linalg.norm(arbor.xyz[1:]-arbor.xyz[arbor.parent[1:]],axis=1).max()<=50.+1e-6
    return {'passed':True,'analytic_generator_only':True,'no_model_fit_or_policy_forecast':True}


def check():
    unopened()
    if (BASE/'frozen.json').exists():raise RuntimeError('Frozen metadata cannot change')
    d=read(BASE/'declaration.json');assert d==current_declaration()
    result={**analytic_check(),'declaration_sha256':sha(BASE/'declaration.json'),'source_sha256':d['source_sha256']}
    save(BASE/'check.json',result,preserve=True);return result


def freeze(review):
    unopened();d=read(BASE/'declaration.json');assert d==current_declaration()
    c=read(BASE/'check.json');assert c['passed'] and c['declaration_sha256']==sha(BASE/'declaration.json') and c['source_sha256']==d['source_sha256']
    review=Path(review);review=review if review.is_absolute() else ROOT/review
    text=review.read_text();assert 'BRM_TEMPLATE_REVIEW: PASS' in text.splitlines()
    assert all(h in text for h in d['source_sha256'].values())
    frozen={'status':'reviewed_frozen_synthetic_only','declaration_sha256':sha(BASE/'declaration.json'),
        'check_sha256':sha(BASE/'check.json'),'review_path':str(review.relative_to(ROOT)),
        'review_sha256':sha(review),'source_sha256':d['source_sha256']}
    if (BASE/'frozen.json').exists():raise RuntimeError('Preserve existing freeze')
    save(BASE/'frozen.json',frozen);return frozen


def verify_frozen():
    d=read(BASE/'declaration.json');f=read(BASE/'frozen.json');assert d==current_declaration()
    assert f['declaration_sha256']==sha(BASE/'declaration.json') and f['check_sha256']==sha(BASE/'check.json')
    assert f['source_sha256']==d['source_sha256'] and f['review_sha256']==sha(ROOT/f['review_path'])
    c=read(BASE/'check.json');assert c['passed'] and c['declaration_sha256']==f['declaration_sha256']
    return d,f


def write_sources(world,split,rows):
    result=[]
    for row in rows:
        arbor,family=make_truth(row);relative=f'worlds/{world}/{split}_truth/{row["neuron"]}.npz'
        save_tree(BASE/relative,arbor)
        result.append({**identity(row),'family':family,'path':relative,'sha256':sha(BASE/relative)})
    relative=f'worlds/{world}/{split}_truth.json';save(BASE/relative,result)
    return result,relative


def load_donors(rows,index):
    T=library();lookup={r['neuron']:r for r in index};donors=[]
    for row in rows:
        entry=lookup[row['neuron']];assert sha(BASE/entry['path'])==entry['sha256']
        donors.append({**identity(row),'table_soma':np.array(row['soma_um']),'swc_soma':np.array(row['soma_um']),
            'tag':row['tag'],'profile':np.array(row['profile']),'arbor':T.G.load_arbor(BASE/entry['path'])})
    return sorted(donors,key=lambda r:r['neuron'])


def fit_regional(world,train,validation,test):
    T=library();R=T.Geo.R
    inputs=R.Inputs(np.array([r['full_id'] for r in train]),np.array([r['animal'] for r in train]),
        np.array([r['soma_um'] for r in train]),np.array([r['tag'] for r in train]))
    model=R.fit_joint(inputs,np.array([r['profile'] for r in train]),2,seed(world,'regional_fit'))
    state=R.snapshot_model(model,inputs.ids);relative=f'worlds/{world}/regional_state.json'
    T.Geo.write_json(BASE/relative,state,exclusive=True)
    proof={'name':world+'_train_only','train_ids':inputs.ids.tolist(),'train_animals':sorted(inputs.animal.tolist()),
        'test_ids':[r['full_id'] for r in validation+test],'test_animals':sorted(r['animal'] for r in validation+test),
        'model_state_sha256':T._digest(state)}
    assert proof['model_state_sha256']==sha(BASE/relative)
    diag=f'worlds/{world}/fit_diagnostics.json';save(BASE/diag,model.diagnostics())
    return T.RegionalState(state,proof),relative,diag


def policy_forecasts(world,split,row,donors,regional,gate=None):
    T=library();results=[]
    for method in METHODS if split=='test' else T.POLICIES:
        arbor,info=T.predict_gated(query(row),donors,regional,gate) if method=='selected' else T.predict_policy(query(row),donors,regional,method)
        relative=f'worlds/{world}/{split}_predictions/{row["neuron"]}/{method}.npz';save_tree(BASE/relative,arbor)
        results.append({'method':method,'path':relative,'sha256':sha(BASE/relative),'info':info})
    return {**identity(row),'soma_um':row['soma_um'],'tag':row['tag'],'predictions':results}


def same_tree(a,b):return np.array_equal(a.xyz,b.xyz) and np.array_equal(a.parent,b.parent) and np.array_equal(a.types,b.types)


def validate_forecasts(manifest,d):
    assert (BASE/'PREDICTIONS_COMPLETE').read_text().strip()==sha(BASE/'prediction_manifest.json')
    assert read(BASE/'prediction_manifest.json')==manifest and manifest['declaration_sha256']==sha(BASE/'declaration.json')
    assert manifest['frozen_sha256']==sha(BASE/'frozen.json') and manifest['source_sha256']==d['source_sha256']
    assert sha(BASE/'synthetic_data.json')==d['synthetic_data_sha256']==manifest['synthetic_data_sha256']
    for relative,h in manifest['artifacts'].items():assert sha(BASE/relative)==h
    known=observed_dataset();seen=set()
    assert set(manifest['worlds'])==set(WORLDS)
    for world,record in manifest['worlds'].items():
        assert [identity(r) for r in record['test_queries']]==[identity(r) for r in known[world]['test']]
        known_by_id={r['neuron']:r for r in known[world]['test']}
        for row in record['test_queries']:
            expected=known_by_id[row['neuron']]
            assert row['soma_um']==expected['soma_um'] and row['tag']==expected['tag']
            assert {r['method'] for r in row['predictions']}==set(METHODS) and len(row['predictions'])==6
            for p in row['predictions']:
                assert p['path']==f'worlds/{world}/test_predictions/{row["neuron"]}/{p["method"]}.npz'
                assert p['path'] not in seen and sha(BASE/p['path'])==p['sha256'];seen.add(p['path'])
                a=library().G.load_arbor(BASE/p['path']);a.validate()
                assert np.allclose(a.xyz[0],row['soma_um'],atol=1e-8,rtol=0)
    assert len(seen)==216


def score_test(manifest,d):
    validate_forecasts(manifest,d)  # GLOBAL barrier precedes all test-truth generation/reads.
    T=library();data=read(BASE/'synthetic_data.json');world_results={}
    for world in WORLDS:
        truth_index,truth_path=write_sources(world,'test',data[world]['test'])
        truths={r['neuron']:r for r in truth_index};train_index=read(BASE/f'worlds/{world}/train_truth.json')
        families={r['neuron']:r['family'] for r in train_index};gate=read(BASE/f'worlds/{world}/gate.json')
        qresults=[]
        for row in manifest['worlds'][world]['test_queries']:
            entry=truths[row['neuron']];assert sha(BASE/entry['path'])==entry['sha256']
            truth=T.G.load_arbor(BASE/entry['path']);models={};trees={};correct={}
            for p in row['predictions']:
                assert sha(BASE/p['path'])==p['sha256'];arbor=T.G.load_arbor(BASE/p['path']);trees[p['method']]=arbor
                models[p['method']]=T.G.geometry_metrics(arbor,truth)
                correct[p['method']]=families[p['info']['selected_donor']]==entry['family']
            qresults.append({**identity(row),'true_family':entry['family'],'truth_path':entry['path'],'truth_sha256':entry['sha256'],
                'metrics':models,'correct_family':correct,'selected_bitwise_baseline':same_tree(trees['selected'],trees['baseline'])})
        summaries={m:T.G.mean_metrics([r['metrics'][m] for r in qresults]) for m in METHODS}
        selected=gate['selected_policy'];baseline=summaries['baseline'];actual=summaries['selected']
        beats=lambda x:actual['distance_um']<x['distance_um'] and actual['f1_200']>x['f1_200']
        matched=selected.replace('conditional_','uniform_') if selected.startswith('conditional_') else None
        recovered=sum(r['correct_family']['selected'] for r in qresults);identical=all(r['selected_bitwise_baseline'] for r in qresults)
        passed=(selected.startswith('conditional_') and recovered==12 and beats(baseline) and beats(summaries[matched])) if world=='informative' else selected=='baseline' and identical
        fit=read(BASE/manifest['worlds'][world]['fit_diagnostics_path'])
        fit_status=[{key:block[key] for key in ('name','n_components','converged','iterations','log_likelihood','regularized_objective')} for block in fit['blocks']]
        world_results[world]={'queries':qresults,'summary':summaries,
            'aggregation':'equal_neuron_equals_equal_animal_one_neuron_per_animal','animals':12,
            'selected_policy':selected,'known_family_recovered':recovered,
            'known_family_total':12,'selected_bitwise_baseline_all':identical,'beats_baseline_both':beats(baseline),
            'beats_matching_uniform_both':beats(summaries[matched]) if matched else None,'declared_pass':bool(passed),
            'regional_fit_status':fit_status,'all_regional_blocks_converged':all(block['converged'] for block in fit['blocks']),
            'test_truth_index_path':truth_path,'test_truth_index_sha256':sha(BASE/truth_path)}
    return {'status':'synthetic_only_complete','source_sha256':d['source_sha256'],
        'prediction_manifest_sha256':sha(BASE/'prediction_manifest.json'),'worlds':world_results,
        'all_declared_worlds_pass':all(r['declared_pass'] for r in world_results.values()),
        'limits':['Synthetic fixture only; no biological performance claim.','Known-family recovery is not unique donor identity among equivalent family templates.',
                  'All gate decisions frozen from validation before test-truth generation.','No test-driven retuning or biological run is performed.']}


def run():
    unopened();d,frozen=verify_frozen();start=time.monotonic()
    with (BASE/'RUN_STARTED.json').open('x') as handle:
        handle.write(canonical({'frozen_sha256':sha(BASE/'frozen.json'),'utc':datetime.now(timezone.utc).isoformat()}))
    try:
        T=library();data=observed_dataset();save(BASE/'synthetic_data.json',data)
        manifest={'declaration_sha256':sha(BASE/'declaration.json'),'frozen_sha256':sha(BASE/'frozen.json'),
            'source_sha256':d['source_sha256'],'synthetic_data_sha256':sha(BASE/'synthetic_data.json'),'artifacts':{},'worlds':{}}
        for world in WORLDS:
            rows=data[world];train_index,tp=write_sources(world,'train',rows['train']);val_index,vp=write_sources(world,'validation',rows['validation'])
            donors=load_donors(rows['train'],train_index)
            regional,sp,dp=fit_regional(world,rows['train'],rows['validation'],rows['test'])
            proof={'training_animals':sorted(r['animal'] for r in rows['train']),
                'validation_animals':sorted(r['animal'] for r in rows['validation']),'test_animals':sorted(r['animal'] for r in rows['test']),
                'training_ids':T.donor_provenance(donors),'validation_ids':[identity(r) for r in rows['validation']],
                'test_ids':[identity(r) for r in rows['test']],'regional_states':[{'model':regional.model,'proof':regional.proof}],
                'state_by_validation_neuron':{r['neuron']:regional.proof['model_state_sha256'] for r in rows['validation']},
                'state_by_test_animal':{r['animal']:regional.proof['model_state_sha256'] for r in rows['test']}}
            proofpath=f'worlds/{world}/proof.json';save(BASE/proofpath,proof)
            val_lookup={r['neuron']:r for r in val_index};validation=[];selection=[]
            for row in rows['validation']:
                forecasts=policy_forecasts(world,'validation',row,donors,regional)
                truth=T.G.load_arbor(BASE/val_lookup[row['neuron']]['path']);fullmetrics={}
                for p in forecasts['predictions']:fullmetrics[p['method']]=T.G.geometry_metrics(T.G.load_arbor(BASE/p['path']),truth)
                forecasts.update({'truth_path':val_lookup[row['neuron']]['path'],'truth_sha256':val_lookup[row['neuron']]['sha256'],'metrics':fullmetrics});validation.append(forecasts)
                selection.append({**identity(row),'metrics':{m:{k:fullmetrics[m][k] for k in ('distance_um','f1_200')} for m in T.POLICIES}})
            valpath=f'worlds/{world}/validation_manifest.json';save(BASE/valpath,validation)
            gate=T.select_policy(selection,proof);gatepath=f'worlds/{world}/gate.json';save(BASE/gatepath,gate)
            tests=[policy_forecasts(world,'test',row,donors,regional,gate) for row in rows['test']]
            manifest['worlds'][world]={'train_truth_index':tp,'validation_truth_index':vp,'regional_state_path':sp,
                'fit_diagnostics_path':dp,'proof_path':proofpath,'validation_manifest_path':valpath,'gate_path':gatepath,
                'gate_sha256':gate['gate_sha256'],'test_queries':tests}
            paths=[tp,vp,sp,dp,proofpath,valpath,gatepath,*[r['path'] for r in train_index+val_index],
                   *[p['path'] for r in validation for p in r['predictions']]]
            manifest['artifacts'].update({path:sha(BASE/path) for path in paths})
            print(f'{world}: validation gate and all72 test forecasts saved',flush=True)
        save(BASE/'prediction_manifest.json',manifest)
        with (BASE/'PREDICTIONS_COMPLETE').open('x') as handle:handle.write(sha(BASE/'prediction_manifest.json')+'\n')
        results=score_test(manifest,d);save(BASE/'metrics.json',results)
        save(BASE/'RUN_COMPLETE.json',{'elapsed_seconds':time.monotonic()-start,'prediction_manifest_sha256':sha(BASE/'prediction_manifest.json'),'metrics_sha256':sha(BASE/'metrics.json')})
        return {'status':results['status'],'all_declared_worlds_pass':results['all_declared_worlds_pass'],
            'worlds':{w:{k:r[k] for k in ('selected_policy','known_family_recovered','declared_pass')} for w,r in results['worlds'].items()}}
    except Exception as error:
        save(BASE/'RUN_FAILED.json',{'type':type(error).__name__,'error':str(error),'elapsed_seconds':time.monotonic()-start});raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=('prepare','check','freeze','run'))
    p.add_argument('--review',default='reports/brm_template_review.md');a=p.parse_args()
    value=prepare() if a.stage=='prepare' else check() if a.stage=='check' else freeze(a.review) if a.stage=='freeze' else run()
    print(canonical({'stage':a.stage,'status':value.get('status'),'passed':value.get('passed'),
        **({'all_declared_worlds_pass':value['all_declared_worlds_pass'],'worlds':value['worlds']} if a.stage=='run' else {})}))
