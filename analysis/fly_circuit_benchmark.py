"""Bounded FlyWire graph-completion DEVELOPMENT; never anatomical mouse transfer.

Known incoming connections predict withheld whole outgoing profiles. Freeze is
metadata/hash-only; predictions are committed before held-out edge scoring.
One female specimen supplies both partitions: this is not animal validation.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[1]
BASE = Path("data/interim/fly_circuit_benchmark")
PROTOCOL = Path("docs/fly_circuit_benchmark_protocol.json")
CODE = "analysis/fly_circuit_benchmark.py"
TESTS = "tests/test_fly_circuit_benchmark.py"
REPORT = Path("reports/fly_circuit_benchmark.json")
METHODS = ("global_frequency", "target_frequency", "incoming_degree_block", "nearest_incoming_copy",
           "reciprocity_logistic", "shared_input_logistic")


class Refusal(ValueError):
    pass


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for b in iter(lambda: stream.read(1<<20), b""):
            h.update(b)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write_new(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        stream.write(canonical(value))


def now():
    return datetime.now(timezone.utc).isoformat()


def choose_nodes(ids, count=2048, fraction=.2, seed="synarch-fly-dev-v1"):
    ids = [str(v) for v in ids]
    if len(ids) != len(set(ids)) or len(ids) < count:
        raise Refusal("metadata IDs are duplicated or fewer than declared")
    rank = lambda tag, v: hashlib.sha256(f"{seed}/{tag}/{v}".encode()).hexdigest()
    selected = sorted(sorted(ids, key=lambda v: rank("universe", v))[:count], key=int)
    holdout = set(sorted(selected, key=lambda v: rank("outgoing_holdout", v))[:max(1, round(count*fraction))])
    return selected, [i for i, v in enumerate(selected) if v in holdout]


def freeze(root=ROOT):
    import pyarrow.parquet as pq
    p = read(root/PROTOCOL)
    if p["status"] != "FROZEN_DEVELOPMENT_AUTHORIZED" or p["species"] != "drosophila_melanogaster":
        raise Refusal("development protocol is not authorized for the fly species")
    for item in p["inputs"]:
        if sha(root/item["path"]) != item["sha256"]:
            raise Refusal("source hash differs")
    with (root/p["metadata_path"]).open() as stream:
        rows = list(csv.DictReader(stream))
    if set(rows[0]) != {"", "Completed"} or any(r["Completed"] != "True" for r in rows):
        raise Refusal("unexpected source neuron-completion metadata")
    ids, queries = choose_nodes([r[""] for r in rows], p["node_count"], p["query_fraction"], p["seed"])
    schema = pq.ParquetFile(root/p["edges_path"]).schema_arrow.names
    if not {"Presynaptic_ID", "Postsynaptic_ID", "Connectivity"} <= set(schema):
        raise Refusal("measured count columns missing; model-signed weights are not substitutes")
    pins = {str(PROTOCOL): sha(root/PROTOCOL), CODE: sha(root/CODE), TESTS: sha(root/TESTS),
            **{v["path"]: v["sha256"] for v in p["inputs"]}}
    record = {"schema": "fly_circuit_development_freeze_v1", "utc": now(), "pins": pins,
              "node_ids": ids, "query_indices": queries, "train_indices": [i for i in range(len(ids)) if i not in set(queries)],
              "edge_outcomes_opened_by_freeze": False, "specimens": 1, "species": p["species"]}
    write_new(root/BASE/"FROZEN.json", record)
    return {"frozen_sha256": sha(root/BASE/"FROZEN.json"), "nodes": len(ids), "query_sources": len(queries)}


def verify(root):
    f = read(root/BASE/"FROZEN.json")
    for path, expected in f["pins"].items():
        if sha(root/path) != expected:
            raise Refusal("frozen input/code/protocol differs: "+path)
    if sha(Path(__file__)) != f["pins"][CODE]:
        raise Refusal("executing code differs from frozen code")
    return f, read(root/PROTOCOL)


def stream_matrix(path, ids, query_indices, truth=False):
    """Drop every query-source record before constructing observed features.

    Only count cells of the explicitly selected universe. An absent pair means
    zero identified synapses in this release, not biological absence.
    """
    import pyarrow.parquet as pq
    keys = np.array([int(v) for v in ids], dtype=np.int64)
    if not np.all(keys[:-1] < keys[1:]):
        raise Refusal("source IDs must be unique numeric-order int64s")
    n = len(ids)
    query = np.zeros(n, bool); query[query_indices] = True
    mat = np.zeros((n, n), float)
    for batch in pq.ParquetFile(path).iter_batches(batch_size=131072, columns=["Presynaptic_ID", "Postsynaptic_ID", "Connectivity"]):
        pre = batch.column(0).to_numpy(); post = batch.column(1).to_numpy()
        i = np.searchsorted(keys, pre); j = np.searchsorted(keys, post)
        valid = (i<n) & (j<n)
        clipped_i, clipped_j = np.minimum(i,n-1), np.minimum(j,n-1)
        valid &= (keys[clipped_i]==pre) & (keys[clipped_j]==post)
        # Apply source holdout before counts/graph features are admitted.
        valid &= query[clipped_i] if truth else ~query[clipped_i]
        weights = batch.column(2).to_numpy()[valid]
        if not np.isfinite(weights).all() or (weights<=0).any() or not np.equal(weights,np.floor(weights)).all():
            raise Refusal("invalid identified-synapse counts")
        np.add.at(mat, (i[valid],j[valid]), weights)
    if not truth:
        mat[query] = np.nan
        train = np.flatnonzero(~query); mat[train,train] = 0  # autapses excluded prospectively
    return mat


def mask_outgoing(matrix, queries):
    result = np.array(matrix, dtype=float, copy=True)
    result[queries,:] = np.nan
    train = np.array([i for i in range(len(result)) if i not in set(queries)])
    result[train,train] = 0
    return result


def logistic(features, labels, offset, l2):
    flat = features.reshape(-1,features.shape[-1]); y = labels.ravel(); off = offset.ravel()
    means, scales = flat.mean(0), flat.std(0)
    scales = np.maximum(scales,1e-8)
    x = np.column_stack([np.ones(len(flat)), (flat-means)/scales])
    def loss(theta):
        z = off+x@theta
        return (float(np.mean(np.logaddexp(0,z)-y*z)+l2*np.dot(theta[1:],theta[1:])/2),
                x.T@(expit(z)-y)/len(y)+np.r_[0,l2*theta[1:]])
    fit = minimize(loss, np.zeros(x.shape[1]), jac=True, method="L-BFGS-B", options={"maxiter":150,"ftol":1e-11,"gtol":1e-7})
    if not fit.success or not np.isfinite(fit.x).all():
        raise Refusal("development fit did not converge: "+str(fit.message))
    def predict(new, new_offset):
        return expit(new_offset+np.column_stack([np.ones(new.size//new.shape[-1]),((new-means)/scales).reshape(-1,new.shape[-1])])@fit.x).reshape(new.shape[:-1])
    return predict, {"coefficients":fit.x.tolist(),"training_means":means.tolist(),"training_scales":scales.tolist(),"converged":True}


def training_context(b, counts, train):
    """For every supervised source i, its whole outgoing row is absent in features.

    Diagonal entries are zero. Removing row i therefore leaves the common-input
    numerator and source incoming degree unchanged, while subtracting b[i,j]
    from each target incoming degree and from its frequency-estimation count.
    This is exactly the explicit source-row mask, without constructing n graphs.
    """
    y=b[:,train]
    degree=b.sum(0)[train]
    shared=b[:,train].T@b[:,train]
    target_degree=degree[None,:]-y
    norm=np.sqrt(degree[:,None]*target_degree)
    common=np.divide(shared,norm,out=np.zeros_like(shared),where=norm>0)
    prior=np.clip((y.sum(0)[None,:]-y+1)/len(train),1e-8,1-1e-8)
    return np.log1p(counts.T),common,np.broadcast_to(np.log1p(degree)[:,None],y.shape),np.log(prior/(1-prior))


def fit_predict(observed, query_indices, l2=.001):
    n = len(observed); q = np.array(query_indices,int); mask = np.ones(n,bool); mask[q]=False; t=np.flatnonzero(mask)
    if not len(q) or len(t)<4 or not np.isnan(observed[q]).all() or not np.isfinite(observed[t]).all() or (observed[t]<0).any():
        raise Refusal("only a graph with wholly missing query outgoing rows is admitted")
    b = (observed[t]>0).astype(float)             # known incoming sources are TRAIN sources only
    y = b[:,t]
    prior = (y.sum(0)+1)/(len(t)+1)             # denominator: nTRAIN-1 nonself pairs plus Beta(1,1)
    offset = np.log(prior/(1-prior))
    degree = b.sum(0)
    shared = b.T@b
    norm = np.sqrt(degree[:,None]*degree[None,:])
    cosine = np.divide(shared,norm,out=np.zeros_like(shared),where=norm>0)
    similarity = cosine[np.ix_(q,t)]
    donors = np.argmax(similarity,axis=1)       # ties by frozen numeric neuron-ID order
    predicted = {"global_frequency":np.full((len(q),len(t)),(y.sum()+1)/(len(t)*(len(t)-1)+2)),
                 "target_frequency":np.broadcast_to(prior,(len(q),len(t))).copy(),
                 "nearest_incoming_copy":.5*y[donors]+.5*prior}
    # A donor's self-target was never eligible; use the TRAIN prior for that slot.
    predicted["nearest_incoming_copy"][np.arange(len(q)),donors]=prior[donors]
    cuts = np.unique(np.quantile(np.log1p(degree[t]),[.25,.5,.75]))
    groups = np.searchsorted(cuts,np.log1p(degree),side="right")
    block = []
    for source in q:
        members = groups[t]==groups[source]
        block.append((y[members].sum(0)+10*prior)/(members.sum()-members.astype(int)+10))
    predicted["incoming_degree_block"] = np.array(block)
    # SOURCE->TARGET reverse connectivity is observed, never the query's hidden row.
    reverse,common,train_degree,train_offset=training_context(b,observed[np.ix_(t,t)],t)
    query_reverse = np.log1p(observed[np.ix_(t,q)].T)
    query_common=cosine[np.ix_(q,t)]
    query_degree = np.broadcast_to(np.log1p(degree[q])[:,None],(len(q),len(t)))
    exclude_self = ~np.eye(len(t),dtype=bool)
    fit_info={}
    for name, ftrain, fquery in (
        ("reciprocity_logistic",np.stack([reverse,train_degree],-1),np.stack([query_reverse,query_degree],-1)),
        ("shared_input_logistic",np.stack([reverse,common,train_degree],-1),np.stack([query_reverse,query_common,query_degree],-1))):
        fn, info=logistic(ftrain[exclude_self],y[exclude_self],train_offset[exclude_self],l2)
        predicted[name]=fn(fquery,np.broadcast_to(offset,(len(q),len(t))).ravel())
        fit_info[name]=info
    if set(predicted)!=set(METHODS) or any(v.shape!=(len(q),len(t)) or not np.isfinite(v).all() or (v<0).any() or (v>1).any() for v in predicted.values()):
        raise Refusal("invalid complete prediction family")
    return predicted,{"fits":fit_info,"copy_donor_indices":t[donors].tolist(),"features":"known reverse count, shared incoming-context cosine, known incoming degree; all query outgoing rows absent"}


def predict(root=ROOT):
    f,p=verify(root); base=root/BASE
    write_new(base/"PREDICT_STARTED.json",{"utc":now(),"freeze_sha256":sha(base/"FROZEN.json")})
    observed=stream_matrix(root/p["edges_path"],f["node_ids"],f["query_indices"])
    predicted,info=fit_predict(observed,f["query_indices"],p["l2"])
    paths={}
    for name in METHODS:
        file=base/f"{name}.npy"
        with file.open("xb") as stream: np.save(stream,predicted[name],allow_pickle=False)
        paths[name]={"path":str(BASE/file.name),"sha256":sha(file)}
    verify(root)
    write_new(base/"PREDICTIONS_COMMITTED.json",{"utc":now(),"freeze_sha256":sha(base/"FROZEN.json"),"predictions":paths,"info":info,"heldout_outcomes_scored":False})
    return {"methods":list(METHODS),"query_sources":len(f["query_indices"]),"commitment_sha256":sha(base/"PREDICTIONS_COMMITTED.json")}


def metrics(labels, probabilities):
    y=np.asarray(labels,dtype=float).ravel(); p=np.asarray(probabilities,dtype=float).ravel()
    if y.shape!=p.shape or not np.isin(y,[0,1]).all() or not np.isfinite(p).all() or (p<0).any() or (p>1).any():
        raise Refusal("metrics require matched binary labels and finite probabilities")
    clipped=np.clip(p,1e-8,1-1e-8)
    ap=None
    if y.sum()>0:
        order=np.argsort(-p,kind="stable"); sy=y[order]; sp=p[order]; tp=np.cumsum(sy)
        ends=np.flatnonzero(np.r_[sp[:-1]!=sp[1:],True]); hits=tp[ends]
        ap=float(np.sum(np.diff(np.r_[0,hits])*(hits/(ends+1)))/y.sum())
    return {"average_precision":ap,"log_loss":float(np.mean(-y*np.log(clipped)-(1-y)*np.log1p(-clipped))),
            "brier":float(np.mean((p-y)**2)),"positives":int(y.sum()),"pairs":len(y)}


def score(root=ROOT):
    f,p=verify(root); base=root/BASE; c=read(base/"PREDICTIONS_COMMITTED.json")
    if c["freeze_sha256"]!=sha(base/"FROZEN.json") or set(c["predictions"])!=set(METHODS): raise Refusal("prediction commitment differs")
    for pin in c["predictions"].values():
        if sha(root/pin["path"])!=pin["sha256"]:raise Refusal("prediction bytes changed")
    write_new(base/"SCORE_STARTED.json",{"utc":now(),"prediction_commitment_sha256":sha(base/"PREDICTIONS_COMMITTED.json")})
    truth=stream_matrix(root/p["edges_path"],f["node_ids"],f["query_indices"],truth=True)
    y=truth[np.ix_(f["query_indices"],f["train_indices"])]>0
    result={name:metrics(y,np.load(root/c["predictions"][name]["path"],allow_pickle=False)) for name in METHODS}
    candidate=result["shared_input_logistic"]; controls=METHODS[:-1]
    wins={name: candidate["average_precision"] is not None and candidate["average_precision"]>result[name]["average_precision"] and candidate["log_loss"]<result[name]["log_loss"] and candidate["brier"]<result[name]["brier"] for name in controls}
    report={"schema":"synarch_fly_circuit_benchmark_v1","species":"drosophila_melanogaster","namespace":"flywire_fafb_v630","specimens":1,
            "status":"DEVELOPMENT_ONLY","nodes":len(f["node_ids"]),"query_sources":len(f["query_indices"]),"target_neurons":len(f["train_indices"]),
            "scores":result,"candidate_beats_control_on_AP_logloss_Brier":wins,"beats_every_control":all(wins.values()),
            "allowed_information":"Observed incoming connectivity from nonquery sources; no genes, no molecular identity, no hidden outgoing edges.",
            "absence_semantics":"Zero identified chemical synapses in the published source table among the selected completed-neuron set; not proof of biological absence or absence of hormonal signaling.",
            "validation":"One specimen; metadata-hash source holdout within one graph; no new-animal, sex, species or causal-intervention validation.",
            "transfer":{"mouse_edges_added":0,"mouse_accuracy_gain_demonstrated":False,"fly_signs_transferred":False,"admitted":"Algorithm/evaluation reference only; any mouse transfer needs paired mouse ground truth."},
            "freeze_sha256":sha(base/"FROZEN.json"),"prediction_commitment_sha256":sha(base/"PREDICTIONS_COMMITTED.json"),"protocol":str(PROTOCOL),"source":p["sources"],"utc":now()}
    verify(root)
    write_new(root/REPORT,report)
    write_new(base/"SCORE_COMPLETE.json",{"utc":now(),"report_sha256":sha(root/REPORT)})
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("command",choices=("freeze","predict","score"));args=parser.parse_args()
    print(canonical(globals()[args.command]()))


if __name__=="__main__":main()
