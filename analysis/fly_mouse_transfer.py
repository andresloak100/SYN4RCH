"""Frozen-fly-parameter transfer to one archival mouse cortex graph, DEVELOPMENT.

Only unsigned chemical-synapse counts and masked incoming graph context are used.
No fly cell, edge, neurotransmitter sign or molecular identity joins a mouse graph.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

ROOT=Path(__file__).resolve().parents[1]
BASE=Path("data/interim/fly_mouse_transfer")
PROTOCOL=Path("docs/fly_mouse_transfer_protocol.json")
CODE="analysis/fly_mouse_transfer.py"
TESTS="tests/test_fly_mouse_transfer.py"
REPORT=Path("reports/fly_mouse_transfer.json")
REFERENCE="analysis/fly_circuit_benchmark.py"
MOUSE_METHODS=("global_frequency","target_frequency","incoming_degree_block","nearest_incoming_copy","reciprocity_logistic","shared_input_logistic")
TRANSFER_METHODS=("frozen_fly_residual","fly_informed_prior")
METHODS=MOUSE_METHODS+TRANSFER_METHODS


class Refusal(ValueError):pass


def canonical(value):return json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False)+"\n"


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1<<20),b""):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())
def now():return datetime.now(timezone.utc).isoformat()


def write_new(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("x") as stream:stream.write(canonical(value))


def reference(root):
    p=read(root/PROTOCOL)
    expected=next(v["sha256"] for v in p["dependencies"] if v["path"]==REFERENCE)
    if sha(root/REFERENCE)!=expected:raise Refusal("frozen fly reference code changed")
    spec=importlib.util.spec_from_file_location("frozen_fly_transfer_math",root/REFERENCE)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def validate_tables(metadata_path,synapse_path,expected_count,expected_synapses):
    """Define the assessed cohort from metadata before examining graph participants."""
    with Path(metadata_path).open() as stream:metadata=list(csv.DictReader(stream))
    if not metadata or not {"cell_type","pt_root_id"}<=set(metadata[0]):raise Refusal("metadata schema differs")
    selected=[row["pt_root_id"] for row in metadata if row["cell_type"]=="e"]
    if len(selected)!=expected_count or len(set(selected))!=len(selected):raise Refusal("assessed metadata-e cohort is incomplete/duplicated")
    if any(not v.isdigit() or int(v)<=0 or str(int(v))!=v for v in selected):raise Refusal("invalid neuron identifier")
    selected_set=set(selected)
    for node in selected:
        if sum(row["pt_root_id"]==node for row in metadata)!=1:raise Refusal("ambiguous selected metadata ID")
    with Path(synapse_path).open() as stream:rows=list(csv.DictReader(stream))
    if len(rows)!=expected_synapses or len({r["id"] for r in rows})!=len(rows):raise Refusal("synapse rows duplicated or snapshot differs")
    if any(r["pre_root_id"] not in selected_set or r["post_root_id"] not in selected_set for r in rows):
        raise Refusal("a synapse endpoint is outside the documented assessed cohort; no silent filtering")
    participants={r[key] for r in rows for key in ("pre_root_id","post_root_id")}
    return {"cohort_ids":sorted(selected,key=int),"cohort_count":len(selected),"isolates_retained":len(selected_set-participants),
            "synapse_rows":len(rows),"self_synapse_rows":sum(r["pre_root_id"]==r["post_root_id"] for r in rows)}


def source_matrix(path,ids,queries,truth=False):
    """Drop all held-out source rows before constructing an observed graph."""
    lookup={str(v):i for i,v in enumerate(ids)};q=set(queries);a=np.zeros((len(ids),len(ids)),float)
    with Path(path).open() as stream:
        for row in csv.DictReader(stream):
            i,j=lookup[row["pre_root_id"]],lookup[row["post_root_id"]]
            if i==j:continue                       # autapses excluded prospectively
            if (i in q)!=truth:continue             # never admit query-outgoing counts into observed context
            a[i,j]+=1                              # one unique identified-synapse row, not cleft/spine volume
    if not truth:a[list(q),:]=np.nan
    return a


def transfer_features(observed,queries,F):
    n=len(observed);q=np.array(queries,int);t=np.array([i for i in range(n) if i not in set(queries)])
    if not np.isnan(observed[q]).all() or not np.isfinite(observed[t]).all():raise Refusal("whole query source rows must be absent")
    b=(observed[t]>0).astype(float);y=b[:,t];degree=b.sum(0)
    reverse,common,d,train_offset=F.training_context(b,observed[np.ix_(t,t)],t)
    train_features=np.stack([reverse,common,d],-1);nonself=~np.eye(len(t),dtype=bool)
    shared=b.T@b;norm=np.sqrt(degree[:,None]*degree[None,:]);cos=np.divide(shared,norm,out=np.zeros_like(shared),where=norm>0)
    query_features=np.stack([np.log1p(observed[np.ix_(t,q)].T),cos[np.ix_(q,t)],
                            np.broadcast_to(np.log1p(degree[q])[:,None],(len(q),len(t)))],-1)
    prior=(y.sum(0)+1)/(len(t)+1)
    return (train_features[nonself],y[nonself],train_offset[nonself],query_features,
            np.broadcast_to(np.log(prior/(1-prior)),(len(q),len(t))).copy())


def validate_parameters(state):
    means=np.asarray(state["training_means"],float);scales=np.asarray(state["training_scales"],float);theta=np.asarray(state["coefficients"],float)
    if means.shape!=(3,) or scales.shape!=(3,) or theta.shape!=(4,) or not np.isfinite(np.r_[means,scales,theta]).all() or (scales<=0).any():
        raise Refusal("frozen fly parameter schema is invalid")
    return means,scales,theta


def converted_prior(fly_state,mouse_means,mouse_scales):
    """Exactly preserve the raw-feature residual function across standardizers."""
    fm,fs,theta=validate_parameters(fly_state)
    mm,ms=np.asarray(mouse_means,float),np.asarray(mouse_scales,float)
    if mm.shape!=(3,) or ms.shape!=(3,) or not np.isfinite(np.r_[mm,ms]).all() or (ms<=0).any():raise Refusal("mouse scaler invalid")
    return np.r_[theta[0]+np.sum(theta[1:]*(mm-fm)/fs),theta[1:]*ms/fs]


def frozen_predict(query_features,offset,fly_state):
    means,scales,theta=validate_parameters(fly_state)
    return expit(offset+theta[0]+((query_features-means)/scales)@theta[1:])


def fit_centered(features,y,offset,query_features,query_offset,center_state,l2):
    """Same mouse objective/scaling as the reference, with a nonzero slope prior.

    Intercept is unpenalized in both models. This changes the convex objective,
    rather than merely initializing the same objective with fly weights.
    """
    means=features.mean(0);scales=np.maximum(features.std(0),1e-8)
    center=np.zeros(4) if center_state is None else converted_prior(center_state,means,scales)
    x=np.column_stack([np.ones(len(features)),(features-means)/scales])
    def loss(theta):
        z=offset+x@theta;delta=theta[1:]-center[1:]
        return float(np.mean(np.logaddexp(0,z)-y*z)+l2*np.dot(delta,delta)/2),x.T@(expit(z)-y)/len(y)+np.r_[0,l2*delta]
    fit=minimize(loss,np.zeros(4),jac=True,method="L-BFGS-B",options={"maxiter":150,"ftol":1e-11,"gtol":1e-7})
    if not fit.success or not np.isfinite(fit.x).all():raise Refusal("mouse prior fit did not converge: "+str(fit.message))
    prediction=expit(query_offset+fit.x[0]+((query_features-means)/scales)@fit.x[1:])
    return prediction,{"coefficients":fit.x.tolist(),"training_means":means.tolist(),"training_scales":scales.tolist(),
                       "converted_fly_center":center.tolist(),"intercept_penalized":False,"l2":l2,"converged":True}


def fit_all(observed,queries,fly_state,F,l2=.001):
    predictions,info=F.fit_predict(observed,queries,l2)
    x,y,offset,qx,qoffset=transfer_features(observed,queries,F)
    predictions["frozen_fly_residual"]=frozen_predict(qx,qoffset,fly_state)
    predictions["fly_informed_prior"],prior_info=fit_centered(x,y,offset,qx,qoffset,fly_state,l2)
    info["fly_informed_prior"]=prior_info
    info["frozen_fly_residual"]={"state":fly_state,"fitted_on_mouse":False,"mouse_calibration":"TRAIN target-frequency logit offset only"}
    if set(predictions)!=set(METHODS) or any(not np.isfinite(v).all() or (v<0).any() or (v>1).any() for v in predictions.values()):raise Refusal("invalid complete transfer forecasts")
    return predictions,info


def freeze(root=ROOT):
    p=read(root/PROTOCOL)
    if p["status"]!="FROZEN_DEVELOPMENT_AUTHORIZED" or not p["root_scientific_admission"] or p["species"]!="mus_musculus":raise Refusal("mouse proxy development not admitted")
    pins={CODE:sha(root/CODE),TESTS:sha(root/TESTS),str(PROTOCOL):sha(root/PROTOCOL)}
    for item in p["inputs"]+p["dependencies"]:
        if sha(root/item["path"])!=item["sha256"]:raise Refusal("input/dependency pin differs")
        pins[item["path"]]=item["sha256"]
    admitted=validate_tables(root/p["metadata_path"],root/p["synapse_path"],p["assessed_count"],p["synapse_rows"])
    F=reference(root);ids,queries=F.choose_nodes(admitted["cohort_ids"],p["assessed_count"],p["query_fraction"],p["seed"])
    fly=read(root/p["fly_prediction_commitment"])["info"]["fits"]["shared_input_logistic"]
    validate_parameters(fly)
    write_new(root/BASE/"FROZEN.json",{"utc":now(),"pins":pins,"node_ids":ids,"query_indices":queries,
              "train_indices":[i for i in range(len(ids)) if i not in set(queries)],"fly_state":fly,
              "admission":admitted,"scope":"nonblind masked single-cortex-specimen DEVELOPMENT"})
    return {"nodes":len(ids),"queries":len(queries),"isolates_retained":admitted["isolates_retained"],"frozen_sha256":sha(root/BASE/"FROZEN.json")}


def verify(root):
    f=read(root/BASE/"FROZEN.json")
    for path,expected in f["pins"].items():
        if sha(root/path)!=expected:raise Refusal("frozen source/dependency changed: "+path)
    if sha(Path(__file__))!=f["pins"][CODE]:raise Refusal("executing code differs")
    return f,read(root/PROTOCOL)


def predict(root=ROOT):
    f,p=verify(root);base=root/BASE
    write_new(base/"PREDICT_STARTED.json",{"utc":now(),"freeze_sha256":sha(base/"FROZEN.json")})
    observed=source_matrix(root/p["synapse_path"],f["node_ids"],f["query_indices"])
    predictions,info=fit_all(observed,f["query_indices"],f["fly_state"],reference(root),p["l2"])
    files={}
    for method in METHODS:
        path=base/f"{method}.npy"
        with path.open("xb") as stream:np.save(stream,predictions[method],allow_pickle=False)
        files[method]={"path":str(BASE/path.name),"sha256":sha(path)}
    verify(root)
    write_new(base/"PREDICTIONS_COMMITTED.json",{"utc":now(),"freeze_sha256":sha(base/"FROZEN.json"),"predictions":files,"info":info,"mouse_heldout_profiles_scored":False})
    return {"methods":list(METHODS),"commitment_sha256":sha(base/"PREDICTIONS_COMMITTED.json")}


def score(root=ROOT):
    f,p=verify(root);base=root/BASE;commit_path=base/"PREDICTIONS_COMMITTED.json";commit_sha=sha(commit_path);c=read(commit_path)
    if c["freeze_sha256"]!=sha(base/"FROZEN.json") or set(c["predictions"])!=set(METHODS):raise Refusal("prediction commitment differs")
    for name,pin in c["predictions"].items():
        if pin["path"]!=str(BASE/f"{name}.npy") or sha(root/pin["path"])!=pin["sha256"]:raise Refusal("prediction pin differs")
    write_new(base/"SCORE_STARTED.json",{"utc":now(),"prediction_commitment_sha256":commit_sha})
    truth=source_matrix(root/p["synapse_path"],f["node_ids"],f["query_indices"],truth=True)
    y=truth[np.ix_(f["query_indices"],f["train_indices"])]>0;F=reference(root)
    scores={name:F.metrics(y,np.load(root/pin["path"],allow_pickle=False)) for name,pin in c["predictions"].items()}
    wins={name:{control:(scores[name]["average_precision"] is not None and scores[name]["average_precision"]>scores[control]["average_precision"] and scores[name]["log_loss"]<scores[control]["log_loss"] and scores[name]["brier"]<scores[control]["brier"]) for control in MOUSE_METHODS} for name in TRANSFER_METHODS}
    verify(root)
    if sha(commit_path)!=commit_sha:raise Refusal("prediction commitment changed duringscore")
    for pin in c["predictions"].values():
        if sha(root/pin["path"])!=pin["sha256"]:raise Refusal("forecast changed duringscore")
    report={"schema":"synarch_fly_mouse_transfer_development_v1","status":"DEVELOPMENT_ONLY","species":"mus_musculus","namespace":"microns_pinky_L23_v185",
            "region":"visualcortexL2/3, NOT hypothalamus","mouse_specimens":1,"fly_specimens":1,"nodes":len(f["node_ids"]),"query_sources":len(f["query_indices"]),"target_nodes":len(f["train_indices"]),
            "admission":f["admission"],"scores":scores,"transfer_beats_each_mouse_control_on_AP_LL_Brier":wins,
            "transfer_beats_mouse_only_on_all_three":{m:wins[m]["shared_input_logistic"] for m in TRANSFER_METHODS},"fit_and_transfer_info":c["info"],
            "limits":["Nonblind masked single-brain cortex development, not new-animal or hypothalamic validation.","Known incoming context predicts recorded outgoing edges; no missing-axon-path or gene prediction tested.","Zero counts mean no identifiedsynapse inside the bounded archived volume, not wholebrain biologicalabsence.","No flyedges, cellidentities or neurotransmitter signs added to mousegraph; no causal intervention/agingclaim.","Archivalv185CSV1961rows differs from2022paper1960; bothselfsynapserows excluded prospectively."],
            "freeze_sha256":sha(base/"FROZEN.json"),"prediction_commitment_sha256":commit_sha,"protocol":str(PROTOCOL),"utc":now()}
    write_new(root/REPORT,report);write_new(base/"SCORE_COMPLETE.json",{"utc":now(),"report_sha256":sha(root/REPORT)})
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("command",choices=("freeze","predict","score"));args=parser.parse_args()
    print(canonical(globals()[args.command]()))


if __name__=="__main__":main()
