"""Synthetic tables/graphs/parameters only; no real mouse outcome fitting."""
from copy import deepcopy
import csv
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m


T=module("fly_mouse_transfer_synthetic_tests","analysis/fly_mouse_transfer.py")
F=module("fly_mouse_transfer_frozen_reference_tests","analysis/fly_circuit_benchmark.py")
STATE={"coefficients":[-.7,.1,.2,.3],"training_means":[.01,.02,.8],"training_scales":[.1,.05,.6],"converged":True}


def toy():
    rng=np.random.default_rng(71);a=rng.binomial(1,.22,(20,20))*rng.integers(1,4,(20,20));np.fill_diagonal(a,0)
    a[-3:,:]=0;a[:,-3:]=0                 # metadata-defined assessed isolates
    return a.astype(float)


def table(path,rows,fields):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w") as stream:
        w=csv.DictWriter(stream,fieldnames=fields);w.writeheader();w.writerows(rows)


@pytest.fixture
def world(tmp_path):
    base=tmp_path/T.BASE;src=base/"source";src.mkdir(parents=True)
    ids=[str(1000+i) for i in range(20)];a=toy()
    metadata=[{"pt_root_id":v,"cell_type":"e"} for v in ids]+[{"pt_root_id":"999","cell_type":"i"}]
    synapses=[]
    for i,j in zip(*np.nonzero(a)):
        for _ in range(int(a[i,j])):synapses.append({"id":str(len(synapses)+1),"pre_root_id":ids[i],"post_root_id":ids[j]})
    m=src/"metadata.csv";s=src/"synapses.csv";table(m,metadata,["pt_root_id","cell_type"]);table(s,synapses,["id","pre_root_id","post_root_id"])
    for path in (T.CODE,T.TESTS,T.REFERENCE):
        f=tmp_path/path;f.parent.mkdir(parents=True,exist_ok=True);f.write_bytes((ROOT/path).read_bytes())
    fly_commit=Path("data/interim/fly_circuit_benchmark/PREDICTIONS_COMMITTED.json")
    f=tmp_path/fly_commit;f.parent.mkdir(parents=True,exist_ok=True);f.write_text(T.canonical({"info":{"fits":{"shared_input_logistic":STATE}}}))
    p={"status":"FROZEN_DEVELOPMENT_AUTHORIZED","root_scientific_admission":True,"species":"mus_musculus","assessed_count":20,"synapse_rows":len(synapses),
       "metadata_path":str(m.relative_to(tmp_path)),"synapse_path":str(s.relative_to(tmp_path)),"fly_prediction_commitment":str(fly_commit),
       "inputs":[{"path":str(f.relative_to(tmp_path)),"sha256":T.sha(f)} for f in (m,s)],
       "dependencies":[{"path":v,"sha256":T.sha(tmp_path/v)} for v in (T.REFERENCE,str(fly_commit))],"seed":"synthetic-transfer","query_fraction":.2,"l2":.001}
    f=tmp_path/T.PROTOCOL;f.parent.mkdir(parents=True,exist_ok=True);f.write_text(T.canonical(p))
    return tmp_path,a,ids,m,s


def test_assayed_cohort_is_metadata_defined_and_includes_isolates(world):
    root,a,ids,m,s=world;r=T.validate_tables(m,s,20,int(a.sum()))
    assert r["cohort_ids"]==ids and r["cohort_count"]==20
    assert r["isolates_retained"]==3


@pytest.mark.parametrize("fault",["unknown_endpoint","duplicate_synapse","duplicate_metadata","bad_root_id"])
def test_ambiguous_or_unassessed_tables_refused(world,fault):
    root,a,ids,m,s=world
    if fault in {"unknown_endpoint","duplicate_synapse"}:
        rows=list(csv.DictReader(s.open()))
        if fault=="unknown_endpoint":rows[0]["post_root_id"]="999"
        else:rows[0]["id"]=rows[1]["id"]
        table(s,rows,["id","pre_root_id","post_root_id"])
    else:
        rows=list(csv.DictReader(m.open()))
        if fault=="duplicate_metadata":rows[0]["pt_root_id"]=rows[1]["pt_root_id"]
        else:rows[0]["pt_root_id"]="1000.0"
        table(m,rows,["pt_root_id","cell_type"])
    with pytest.raises(T.Refusal):T.validate_tables(m,s,20,int(a.sum()))


def test_csv_count_orientation_and_query_row_mask(world):
    root,a,ids,m,s=world;q=[2,5,13]
    observed=T.source_matrix(s,ids,q)
    assert np.array_equal(observed,F.mask_outgoing(a,q),equal_nan=True)
    truth=T.source_matrix(s,ids,q,truth=True)
    assert np.array_equal(truth[q],a[q])
    assert not truth[[i for i in range(20) if i not in q]].any()


def test_parameter_conversion_preserves_raw_residual_function():
    means=np.array([.08,.15,1.6]);scales=np.array([.2,.3,.8]);converted=T.converted_prior(STATE,means,scales)
    x=np.random.default_rng(18).normal(size=(40,3));fm,fs,theta=T.validate_parameters(STATE)
    one=theta[0]+((x-fm)/fs)@theta[1:];two=converted[0]+((x-means)/scales)@converted[1:]
    assert np.allclose(one,two,rtol=1e-13,atol=1e-13)


def test_frozen_parameters_and_scalers_are_not_changed_by_mouse_features():
    state=deepcopy(STATE);x=np.random.default_rng(2).normal(size=(3,4,3));offset=np.zeros((3,4))
    first=T.frozen_predict(x,offset,state);second=T.frozen_predict(x,offset+1,state)
    assert state==STATE
    assert (second>first).all()             # only the declared observed mouse target offset changes


def test_zero_centered_objective_exactly_matches_frozen_mouse_only_model():
    a=toy();q=[2,5,13,18];obs=F.mask_outgoing(a,q)
    x,y,offset,qx,qoffset=T.transfer_features(obs,q,F)
    predicted,_=F.fit_predict(obs,q,.001)
    direct,info=T.fit_centered(x,y,offset,qx,qoffset,None,.001)
    assert np.allclose(direct,predicted["shared_input_logistic"],rtol=1e-12,atol=1e-12)
    assert info["intercept_penalized"] is False and info["l2"]==.001


def test_arbitrary_hidden_outgoing_profiles_leave_all_eight_forecasts_unchanged():
    a=toy();q=[2,5,13,18];changed=a.copy();changed[q]=np.random.default_rng(8).integers(0,8,(4,20))
    state=deepcopy(STATE);one,info1=T.fit_all(F.mask_outgoing(a,q),q,state,F);two,info2=T.fit_all(F.mask_outgoing(changed,q),q,state,F)
    assert info1==info2 and state==STATE and set(one)==set(T.METHODS)
    for method in T.METHODS:assert np.array_equal(one[method],two[method])


def test_training_context_still_excludes_each_supervised_sources_own_row():
    a=toy();q=[2,5,13,18];obs=F.mask_outgoing(a,q);t=np.array([i for i in range(20) if i not in q]);b=(obs[t]>0).astype(float)
    reverse,common,degree,offset=F.training_context(b,obs[np.ix_(t,t)],t)
    for i in range(len(t)):
        explicit=b.copy();explicit[i]=0;d=explicit.sum(0);num=explicit[:,t[i]]@explicit[:,t];den=np.sqrt(d[t[i]]*d[t]);expected=np.divide(num,den,out=np.zeros_like(num),where=den>0)
        assert np.allclose(common[i],expected)
        prior=np.clip((explicit[:,t].sum(0)+1)/len(t),1e-8,1-1e-8)
        assert np.allclose(offset[i],np.log(prior/(1-prior)))


@pytest.mark.parametrize("bad",[0.,float("nan"),float("inf")])
def test_invalid_fly_scales_refused(bad):
    state=deepcopy(STATE);state["training_scales"][0]=bad
    with pytest.raises(T.Refusal):T.converted_prior(state,np.zeros(3),np.ones(3))


def test_complete_synthetic_lifecycle_same_pairs_and_no_causal_hypothalamus_claim(world):
    root,a,ids,m,s=world
    result=T.freeze(root);assert result["nodes"]==20 and result["isolates_retained"]==3
    T.predict(root);report=T.score(root)
    assert report["status"]=="DEVELOPMENT_ONLY" and report["mouse_specimens"]==1
    assert set(report["scores"])==set(T.METHODS)
    assert len({v["pairs"] for v in report["scores"].values()})==1
    assert set(report["transfer_beats_mouse_only_on_all_three"])==set(T.TRANSFER_METHODS)
    assert report["fit_and_transfer_info"]["frozen_fly_residual"]["fitted_on_mouse"] is False
    with pytest.raises(FileExistsError):T.score(root)


def test_frozen_fly_dependency_mutation_blocks_before_any_mouse_fit(world):
    root,a,ids,m,s=world;T.freeze(root)
    file=root/"data/interim/fly_circuit_benchmark/PREDICTIONS_COMMITTED.json";file.write_text(file.read_text()+" ")
    with pytest.raises(T.Refusal,match="changed"):T.predict(root)
    assert not (root/T.BASE/"PREDICT_STARTED.json").exists()


def test_forecast_tamper_blocks_before_mouse_hidden_outcome_score(world,monkeypatch):
    root,a,ids,m,s=world;T.freeze(root);T.predict(root)
    f=root/T.BASE/"frozen_fly_residual.npy";f.write_bytes(f.read_bytes()+b"tamper")
    monkeypatch.setattr(T,"source_matrix",lambda *a,**k:pytest.fail("hiddenoutcome beforepredictioncheck"))
    with pytest.raises(T.Refusal,match="prediction pin"):T.score(root)
    assert not (root/T.BASE/"SCORE_STARTED.json").exists()
