"""Boundary and tiny integration tests, never the fixed synthetic benchmark."""
from copy import deepcopy
import importlib.util
from pathlib import Path
import sys

import numpy as np

PATH=Path(__file__).resolve().parents[1]/'analysis/brm_template_benchmark.py'
SPEC=importlib.util.spec_from_file_location('brm_template_benchmark_tested',PATH)
B=importlib.util.module_from_spec(SPEC);sys.modules[SPEC.name]=B;SPEC.loader.exec_module(B)


def test_paired_positions_and_disjoint_animals_match_declared_information_controls():
    data=B.observed_dataset();all_animals=[]
    for world,splits in data.items():
        for split,rows in splits.items():
            assert len(rows)==B.SPLITS[split]
            all_animals.extend(r['animal'] for r in rows)
            for i in range(0,len(rows),2):assert rows[i]['soma_um']==rows[i+1]['soma_um']
            assert sum(B.true_family(r)==0 for r in rows)==len(rows)//2
            if world=='no_information':
                assert all(r['soma_um']==[0.,0.,0.] and r['tag']=='U' for r in rows)
                if split=='train':assert len({tuple(r['profile']) for r in rows})==1
    assert len(all_animals)==len(set(all_animals))==120


def test_query_cannot_carry_hidden_profiles_families_or_geometry():
    row=B.known_rows('informative','test')[0];first=B.query(row)
    changed=deepcopy(row);changed.update(profile=[1e30]*49,family=999,arbor=object(),truth={'target':-1e30})
    second=B.query(changed)
    assert set(first.__dict__)=={'neuron','animal','full_id','table_soma','swc_soma','tag'}
    for field in first.__dict__:
        assert np.array_equal(getattr(first,field),getattr(second,field))


def test_tiny_synthetic_fit_state_hash_and_allowed_query_prediction(tmp_path,monkeypatch):
    """Separate four-row fixture: no declared world generation or test scoring."""
    monkeypatch.setattr(B,'BASE',tmp_path)
    T=B.library();train=[]
    for i in range(4):
        animal=f'SMKT{i}';tag='A' if i%2==0 else 'B';neuron=animal+'_001'
        train.append({'neuron':neuron,'animal':animal,'full_id':neuron+'_'+tag+'_SYN',
            'soma_um':[0.,0.,0.],'tag':tag,'profile':[float(i%2)*3.]*49})
    def held(animal):return {'neuron':animal+'_001','animal':animal,'full_id':animal+'_001_A_SYN','soma_um':[0.,0.,0.],'tag':'A'}
    val,test=[held('SMKV0')],[held('SMKE0')]
    regional,path,_=B.fit_regional('tiny_fixture',train,val,test)
    assert regional.proof['model_state_sha256']==T._digest(regional.model)==B.sha(tmp_path/path)
    donors=[]
    for row in train:
        arbor=T.G.Arbor(np.array([[0.,0.,0.],[50.,0.,0.]]),np.array([-1,0]),np.array([1,2]))
        donors.append({**B.identity(row),'table_soma':np.zeros(3),'swc_soma':np.zeros(3),
            'tag':row['tag'],'profile':np.array(row['profile']),'arbor':arbor})
    prediction,info=T.predict_policy(B.query(test[0]),donors,regional,'baseline')
    assert np.array_equal(prediction.xyz,donors[0]['arbor'].xyz)
    assert info['executed_policy']=='baseline'
