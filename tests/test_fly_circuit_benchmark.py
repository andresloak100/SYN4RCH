"""Synthetic directed graphs only; never read the real fly edge outcomes."""
import csv
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("fly_circuit_development_tests",ROOT/"analysis/fly_circuit_benchmark.py")
F=importlib.util.module_from_spec(spec);sys.modules[spec.name]=F;spec.loader.exec_module(F)


def graph(n=20):
    rng=np.random.default_rng(77)
    a=rng.binomial(1,.25,(n,n))*rng.integers(1,8,(n,n))
    a=np.maximum(a,a.T*(rng.random((n,n))<.45))
    np.fill_diagonal(a,0)
    return a.astype(float)


def parquet(path,a,ids):
    i,j=np.nonzero(a)
    pq.write_table(pa.table({"Presynaptic_ID":np.array(ids,dtype=np.int64)[i],"Postsynaptic_ID":np.array(ids,dtype=np.int64)[j],
                            "Connectivity":a[i,j].astype(np.int64),"Excitatory x Connectivity":-a[i,j].astype(np.int64)}),path)


@pytest.fixture
def world(tmp_path):
    base=tmp_path/F.BASE;base.mkdir(parents=True)
    ids=list(range(1000,1020));a=graph()
    metadata=base/"metadata.csv"
    with metadata.open("w") as stream:
        w=csv.writer(stream);w.writerow(["","Completed"]);w.writerows((i,"True") for i in ids)
    edges=base/"edges.parquet";parquet(edges,a,ids)
    for name in (F.CODE,F.TESTS):
        file=tmp_path/name;file.parent.mkdir(parents=True,exist_ok=True);file.write_bytes((ROOT/name).read_bytes())
    inputs=[{"path":str(p.relative_to(tmp_path)),"sha256":F.sha(p)} for p in (metadata,edges)]
    p={"status":"FROZEN_DEVELOPMENT_AUTHORIZED","species":"drosophila_melanogaster","inputs":inputs,
       "metadata_path":inputs[0]["path"],"edges_path":inputs[1]["path"],"node_count":20,"query_fraction":.2,"seed":"synthetic","l2":.001,"sources":[{"synthetic":True}]}
    target=tmp_path/F.PROTOCOL;target.parent.mkdir(parents=True);target.write_text(F.canonical(p))
    return tmp_path,a,ids


def test_id_selection_is_metadata_only_reproducible_and_order_independent():
    ids=list(range(1000,1100))
    first=F.choose_nodes(ids,20,.2,"fixed")
    assert first==F.choose_nodes(ids[::-1],20,.2,"fixed")
    assert len(first[0])==20 and len(first[1])==4
    with pytest.raises(F.Refusal):F.choose_nodes([1,1,2],3)


def test_entire_query_outgoing_swap_leaves_every_forecast_bitwise_unchanged():
    a=graph();q=[2,7,11,17]
    b=a.copy();b[q]=np.random.default_rng(88).integers(0,40,(len(q),len(a)))
    one,info1=F.fit_predict(F.mask_outgoing(a,q),q)
    two,info2=F.fit_predict(F.mask_outgoing(b,q),q)
    assert set(one)==set(F.METHODS)
    assert info1==info2
    for name in F.METHODS:assert np.array_equal(one[name],two[name])


def test_unmasked_query_outgoing_is_refused():
    with pytest.raises(F.Refusal,match="wholly missing"):
        F.fit_predict(graph(),[2,7])


def test_each_training_source_context_matches_explicit_source_row_deletion():
    a=graph();q=[2,7,11];obs=F.mask_outgoing(a,q);t=np.array([i for i in range(len(a)) if i not in q]);b=(obs[t]>0).astype(float)
    reverse,common,degree,offset=F.training_context(b,obs[np.ix_(t,t)],t)
    for i in range(len(t)):
        explicit=b.copy();explicit[i,:]=0
        d=explicit.sum(0);num=explicit[:,t[i]]@explicit[:,t]
        den=np.sqrt(d[t[i]]*d[t]);co=np.divide(num,den,out=np.zeros_like(num),where=den>0)
        prior=np.clip((explicit[:,t].sum(0)+1)/len(t),1e-8,1-1e-8)
        assert np.allclose(common[i],co)
        assert np.allclose(degree[i],np.log1p(d[t[i]]))
        assert np.allclose(offset[i],np.log(prior/(1-prior)))
        raw=obs[np.ix_(t,t)].copy();raw[i,:]=0
        assert np.array_equal(reverse[i],np.log1p(raw[:,i]))


def test_a_training_sources_own_outgoing_row_cannot_change_its_context():
    a=graph();q=[2,7,11];t=np.array([i for i in range(len(a)) if i not in q]);i=4
    obs=F.mask_outgoing(a,q);b=(obs[t]>0).astype(float)
    before=F.training_context(b,obs[np.ix_(t,t)],t)
    obs[t[i],:]=np.random.default_rng(8).integers(0,4,len(a));obs[t[i],t[i]]=0
    after=F.training_context((obs[t]>0).astype(float),obs[np.ix_(t,t)],t)
    for x,y in zip(before,after):assert np.array_equal(x[i],y[i])


def test_block_and_copy_controls_do_not_turn_excluded_autapses_into_negatives():
    a=graph();q=[2,7,11];obs=F.mask_outgoing(a,q);t=np.array([i for i in range(len(a)) if i not in q])
    predicted,info=F.fit_predict(obs,q);b=(obs[t]>0).astype(float);y=b[:,t];degree=b.sum(0)
    prior=(y.sum(0)+1)/(len(t)+1);cuts=np.unique(np.quantile(np.log1p(degree[t]),[.25,.5,.75]))
    groups=np.searchsorted(cuts,np.log1p(degree),side="right")
    for i,source in enumerate(q):
        members=groups[t]==groups[source]
        for j in range(len(t)):
            eligible=members.sum()-int(members[j])
            expected=(y[members,j].sum()+10*prior[j])/(eligible+10)
            assert predicted["incoming_degree_block"][i,j]==pytest.approx(expected)
        donor_position=int(np.flatnonzero(t==info["copy_donor_indices"][i])[0])
        assert predicted["nearest_incoming_copy"][i,donor_position]==prior[donor_position]


def test_streaming_orientation_count_semantics_and_global_query_mask(world):
    root,a,ids=world;q=[2,7,11]
    path=root/F.BASE/"edges.parquet"
    observed=F.stream_matrix(path,[str(v) for v in ids],q)
    expected=F.mask_outgoing(a,q)
    assert np.array_equal(observed,expected,equal_nan=True)
    truth=F.stream_matrix(path,[str(v) for v in ids],q,truth=True)
    assert np.array_equal(truth[q],a[q])
    assert not truth[[i for i in range(len(a)) if i not in q]].any()
    # Negative signed-model weights were deliberately supplied: actual countpositives survive.
    assert (observed[np.isfinite(observed)]>=0).all()


def test_average_precision_ties_and_probability_losses_have_known_values():
    m=F.metrics([1,0,1],[.5,.5,.5])
    assert m["average_precision"]==pytest.approx(2/3)
    assert m["log_loss"]==pytest.approx(np.log(2))
    assert m["brier"]==pytest.approx(.25)
    assert F.metrics([1,0,1],[.5,.5,.1])["average_precision"]==pytest.approx(7/12)
    assert F.metrics([0,0],[.2,.2])["average_precision"] is None


@pytest.mark.parametrize("probabilities",[[float("nan"),.5],[-.1,.5],[1.1,.5]])
def test_invalid_probabilities_refused(probabilities):
    with pytest.raises(F.Refusal):F.metrics([0,1],probabilities)


def test_freeze_reads_metadata_schema_hashes_but_never_edge_values(world,monkeypatch):
    root,a,ids=world
    monkeypatch.setattr(F,"stream_matrix",lambda *a,**k:pytest.fail("edge outcome access duringfreeze"))
    result=F.freeze(root)
    assert result["nodes"]==20 and result["query_sources"]==4
    assert F.read(root/F.BASE/"FROZEN.json")["edge_outcomes_opened_by_freeze"] is False


def test_complete_synthetic_lifecycle_and_species_transfer_is_not_claimed(world):
    root,a,ids=world
    F.freeze(root);F.predict(root);result=F.score(root)
    assert result["status"]=="DEVELOPMENT_ONLY" and result["specimens"]==1
    assert result["transfer"]["mouse_edges_added"]==0
    assert result["transfer"]["mouse_accuracy_gain_demonstrated"] is False
    assert result["transfer"]["fly_signs_transferred"] is False
    assert set(result["scores"])==set(F.METHODS)
    assert len({v["pairs"] for v in result["scores"].values()})==1
    with pytest.raises(FileExistsError):F.score(root)


def test_frozen_source_change_refuses_before_predict(world):
    root,a,ids=world;F.freeze(root)
    file=root/F.BASE/"metadata.csv";file.write_text(file.read_text()+"\n")
    with pytest.raises(F.Refusal,match="frozen input"):F.predict(root)
    assert not (root/F.BASE/"PREDICT_STARTED.json").exists()


def test_prediction_tamper_refuses_before_opening_heldout_edges(world,monkeypatch):
    root,a,ids=world;F.freeze(root);F.predict(root)
    file=root/F.BASE/"global_frequency.npy";file.write_bytes(file.read_bytes()+b"tamper")
    monkeypatch.setattr(F,"stream_matrix",lambda *a,**k:pytest.fail("truth beforepredictionverification"))
    with pytest.raises(F.Refusal,match="prediction bytes"):F.score(root)
    assert not (root/F.BASE/"SCORE_STARTED.json").exists()
