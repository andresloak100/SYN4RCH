"""SYNARCH masked-variable benchmark on the Jiao et al. 2025 TRAINING partition: commands.

    python3 -I analysis/synarch_recon_benchmark.py inputs                 # input audit, training samples only
    python3 -I analysis/synarch_recon_benchmark.py splits                 # the whole-sample folds (inputs only)
    python3 -I analysis/synarch_recon_benchmark.py synthetic --size tiny  # full pipeline on synthetic worlds
    python3 -I analysis/synarch_recon_benchmark.py timing                 # one full fit at the real size, synthetic
    python3 -I analysis/synarch_recon_benchmark.py freeze --reviewed "..."   # pre-run declaration (after review)
    python3 -I analysis/synarch_recon_benchmark.py run --workers 6        # ONCE, only after the freeze
    python3 -I analysis/synarch_recon_benchmark.py score | demo | report | status
    python3 -I analysis/synarch_recon_benchmark.py serving                # RESTORE the all-training states of a
                                                                          # finished run (the run itself fits them
                                                                          # and its completion record binds them)
    python3 -I analysis/synarch_recon_benchmark.py abort --reason "..."   # set aside a freeze that was never run, or
                                                                          # a start that broke before any outcome-
                                                                          # derived file was written
    python3 -I analysis/synarch_recon_benchmark.py nullrate --size tiny --worlds 100 --workers 6
                                                                          # no-signal synthetic worlds: how often
                                                                          # does the status rule say 'supported'?

Binding rules: reports/synarch_recon_protocol.md. Order of work: everything is developed on synthetic
data and on inputs; `run` refuses the real cohort until `freeze` has written the pre-run declaration,
re-checks every pinned hash, and writes an exclusive start marker so the benchmark runs once.
No command here opens a trace, downloads anything, or reads a row of a validation or test brain sample.
"""
from __future__ import annotations

import os

for _name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import argparse
import importlib.util
import json
import platform
import shutil
import sys
import time
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy


def _sibling(name):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


core = _sibling("synarch_recon")
mods = _sibling("synarch_recon_modules")

ROOT = core.ROOT
CODE_FILES = ("analysis/synarch_recon.py", "analysis/synarch_recon_modules.py", "analysis/synarch_recon_benchmark.py",
              "tests/test_synarch_recon.py", "analysis/brm_joint.py")
PROTOCOL = "reports/synarch_recon_protocol.md"
RESULTS = "reports/synarch_recon_results.md"
CONTRACT = "data/interim/brm_data_contract.json"

EXPOSURE = ("Exploratory. The 49-region table has been on disk since 2026-11-10 and another agent's regional model "
            "(reports/brm_real_results.md) was fitted and scored on all 166 brain samples, training samples included. "
            "The writer of this benchmark read that report before writing the protocol. Nothing here is confirmatory.")

# The configuration that the pre-run declaration freezes for the real cohort.
CONFIG = {
    "schema": core.SCHEMA,
    "salt": "SYNARCH-recon-2026-11-10",
    "outer_folds": 5,
    "inner_folds": 3,
    "interval_level": 0.90,
    "tasks": ["T1", "T2", "T3", "T4"],
    "engine": {"representations": ["log", "hurdle"], "k_grid": [1, 2, 4, 8, 12, 16], "n_init": 5, "max_iter": 300,
               "tol": 1e-7, "variance_floor": 1e-3, "dirichlet_alpha": 1.1},
    "mean": {"weightings": ["sample", "neuron"]},      # MEAN and IDENTITY: who counts once, chosen on inner folds
    "knn": {"k_grid": [5, 10, 20, 40, 80], "pools": ["all", "same_label"], "vote_alpha": 1.0},
    "ridge": {"lambda_grid": [1e-3, 1e-2, 1e-1, 1.0, 10.0]},
    "logistic": {"lambda_grid": [1e-4, 1e-3, 1e-2, 1e-1], "max_iter": 500},
    "bootstrap": {"resamples": 10000, "level": 0.95},
    "min_relative_gain": 0.02,       # the engine must be at least this much better than a comparator, relative
                                     # to the comparator's primary metric, to count as CLEARLY better
    "tie_tolerance": 1e-12,          # paired differences below this (relative to the metric) are exact ties
    "coverage_tolerance": 0.05,
    "min_samples_for_interval": 5,
}
TINY_CONFIG = dict(CONFIG, engine={"representations": ["log", "hurdle"], "k_grid": [1, 2, 4], "n_init": 2,
                                   "max_iter": 60, "tol": 1e-6, "variance_floor": 1e-3, "dirichlet_alpha": 1.1},
                   knn={"k_grid": [3, 9], "pools": ["all", "same_label"], "vote_alpha": 1.0},
                   ridge={"lambda_grid": [1e-2, 1.0]}, logistic={"lambda_grid": [1e-3, 1e-1], "max_iter": 100},
                   bootstrap={"resamples": 500, "level": 0.95})

ENGINE, ENGINE_MEDIAN = "ENGINE", "ENGINE_MEDIAN"
METRICS = {"profile": ("mse", "mae", "corr"),
           "position": ("euclid_mean", "euclid_median", "abs_AP", "abs_DV", "abs_ML"),
           "label": ("accuracy", "logloss")}
HIGHER_IS_BETTER = {"corr", "accuracy"}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def say(*parts):
    print(*parts, flush=True)


# ----------------------------------------------------------------------------------------------
# worlds and configuration
# ----------------------------------------------------------------------------------------------


def get_world(directory=None):
    return core.real_world() if directory is None else core.synthetic_world(directory)


def config_for(world):
    if not world.synthetic:
        return json.loads(json.dumps(CONFIG))
    meta = core.read_json(Path(world.manifest).parent / "world.json")
    return json.loads(json.dumps(meta.get("config") or CONFIG))


def code_hashes():
    return {rel: (core.sha_file(ROOT / rel) if (ROOT / rel).exists() else None) for rel in CODE_FILES}


def runtime():
    import openpyxl
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
            "openpyxl": openpyxl.__version__}


def task_declaration(tasks):
    return {name: {"title": t.title, "kind": t.kind, "target": t.target_variable,
                   "target_coordinates": list(t.target_coordinates),
                   "evidence_continuous_fields": list(t.evidence_cont), "evidence_categorical_fields": list(t.evidence_cat),
                   "primary_metric": t.primary, "higher_is_better": t.higher_is_better,
                   "metrics": list(METRICS[t.kind]), "copy_scaling": t.copy_scaling, "knn_scalings": list(t.knn_scalings),
                   "ridge_degrees": list(t.ridge_degrees), "declared_measured_axon_input": t.declared_measured_axon_input,
                   "target_scale": t.scale, "clip_min": t.clip_min, "description": t.description}
            for name, t in tasks.items()}


def engine_grid(config):
    """Engine settings in the declared order: representation, then number of components."""
    return [[r, int(k)] for r in config["engine"]["representations"] for k in config["engine"]["k_grid"]]


def mean_grid(config):
    return [{"weighting": w} for w in config["mean"]["weightings"]]


def knn_grid(config, task, has_label):
    pools = [p for p in config["knn"]["pools"] if p == "all" or has_label]
    return [{"k": int(k), "pool": pool, "scaling": scaling} for scaling in task.knn_scalings for pool in pools
            for k in config["knn"]["k_grid"]]


def ridge_grid(config, task):
    if task.kind == "label":
        return [{"lambda": float(l)} for l in config["logistic"]["lambda_grid"]]
    return [{"lambda": float(l), "degree": int(d)} for d in task.ridge_degrees for l in config["ridge"]["lambda_grid"]]


def tuned_grid(config, task, method, has_label):
    """The declared grid of one tuned comparator, in the declared order (ties go to the earlier setting)."""
    if method in ("MEAN", "IDENTITY"):
        return mean_grid(config)
    if method == "KNN":
        return knn_grid(config, task, has_label)
    if method == "RIDGE":
        return ridge_grid(config, task)
    raise core.ReconError(f"{method} has no tuned setting")


# ----------------------------------------------------------------------------------------------
# context: the training partition as model rows, with its folds
# ----------------------------------------------------------------------------------------------


@dataclass
class Context:
    world: object
    config: dict
    partition: object
    table: object
    registry: object
    layout: object
    rows: object
    tasks: dict
    splits: dict
    outer_of: np.ndarray       # outer fold of every row
    inner_of: dict             # outer fold -> inner fold of every row (-1 for the rows held out by that outer fold)
    classes: list              # every source label of the training partition, sorted


def fold_arrays(table, splits):
    outer_of = np.full(len(table), -1)
    for fold, held in splits["outer"].items():
        outer_of[np.isin(table.sample, held)] = int(fold)
    inner_of = {}
    for fold, parts in splits["inner"].items():
        a = np.full(len(table), -1)
        for inner, members in parts.items():
            a[np.isin(table.sample, members)] = int(inner)
        if (a[outer_of == int(fold)] != -1).any() or (a[outer_of != int(fold)] == -1).any():
            raise core.ReconError("inner folds do not partition the outer-training rows")
        inner_of[int(fold)] = a
    if (outer_of < 0).any():
        raise core.ReconError("a training row has no outer fold")
    return outer_of, inner_of


def make_splits(table, config):
    label_of = {s: l for s, l in zip(table.sample.tolist(), table.label.tolist())}
    return core.nested_folds(sorted(set(table.sample.tolist())), label_of, outer=config["outer_folds"],
                             inner=config["inner_folds"], salt=config["salt"])


def build_context(world, config=None):
    """Load the training partition WITH outcomes. For the real cohort this is refused until the
    pre-run declaration is verified and the run has been opened (the order of work is binding)."""
    config = config or config_for(world)
    if not world.synthetic:
        verify_freeze(world)
        if not (Path(world.out) / "run_started.json").exists():
            raise core.ReconError("the real outcomes are opened only by `run`, after its start marker is written")
    partition = core.TrainingPartition(world)
    table = core.load_training_table(world, with_outcomes=True, partition=partition)
    registry = core.Registry.from_contract(ROOT / CONTRACT)
    layout = mods.jiao_layout(registry, world.regions)
    rows = mods.rows_from_table(table, layout)
    tasks = {name: t for name, t in mods.build_tasks(layout, registry, world.t4_observed).items()
             if name in config["tasks"]}
    splits = make_splits(table, config)
    outer_of, inner_of = fold_arrays(table, splits)
    return Context(world, config, partition, table, registry, layout, rows, tasks, splits, outer_of, inner_of,
                   sorted(set(table.label.tolist())))


# ----------------------------------------------------------------------------------------------
# scoring helpers
# ----------------------------------------------------------------------------------------------


def target_matrix(task, rows):
    if task.kind == "label":
        return np.asarray(rows.cat[task.target_fields[0]], dtype=object)
    return rows.matrix(task.target_fields)


def align_classes(probability, classes, all_classes):
    out = np.zeros((len(probability), len(all_classes)))
    position = {c: i for i, c in enumerate(all_classes)}
    for j, c in enumerate(classes):
        out[:, position[c]] = probability[:, j]
    return out


def neuron_scores(task, prediction, truth, classes=None, seen=None):
    """Per-neuron scores of one method. Label task: NaN where the true label was absent from training."""
    if task.kind == "profile":
        return core.profile_scores(prediction, truth)
    if task.kind == "position":
        return core.position_scores(prediction, truth, task.target_coordinates)
    scores = core.label_scores(prediction, classes, truth)
    if seen is not None:
        for key in scores:
            scores[key] = np.where(seen, scores[key], np.nan)
    return scores


def sample_metrics(task, scores, sample):
    """{metric: {brain sample: value}} from per-neuron scores."""
    if task.kind == "position":
        out = {"euclid_mean": core.by_sample(scores["euclid"], sample),
               "euclid_median": core.by_sample(scores["euclid"], sample, how="median")}
        for axis in task.target_coordinates:
            out["abs_" + axis] = core.by_sample(scores["abs_" + axis], sample)
        return out
    return {metric: core.by_sample(values, sample) for metric, values in scores.items()}


def primary_score(task, prediction, truth, sample, classes=None, seen=None):
    """(primary metric, tie-break) of a prediction, aggregated per brain sample and then over samples."""
    per = sample_metrics(task, neuron_scores(task, prediction, truth, classes, seen), sample)
    value = core.over_samples(per[task.primary])
    tie = core.over_samples(per["logloss"]) if task.kind == "label" else 0.0
    return float(value), float(tie)


def better(task, a, b):
    """True when score a = (primary, tie-break) is strictly better than b."""
    if a[0] != b[0] and not (np.isnan(a[0]) or np.isnan(b[0])):
        return a[0] > b[0] if task.higher_is_better else a[0] < b[0]
    if np.isnan(b[0]) and not np.isnan(a[0]):
        return True
    return a[1] < b[1]


def pick(task, scored):
    """First best entry of [(params, (primary, tie)), ...] in the declared grid order."""
    best = None
    for params, score in scored:
        if best is None or better(task, score, best[1]):
            best = (params, score)
    return best


# ----------------------------------------------------------------------------------------------
# jobs (each is a pure function of the context and its arguments)
# ----------------------------------------------------------------------------------------------

_CTX = None


def _init_worker(directory):
    global _CTX
    _CTX = build_context(get_world(directory))


def _seen_mask(truth, classes):
    known = set(classes)
    return np.array([t in known for t in truth])


def _fit(ctx, index, representation, components, seed_parts):
    ctx.partition.require_training(ctx.rows.group[index], "joint-model fit")
    budget = ctx.config["engine"]
    seed = core.stable_seed(ctx.config["salt"], "engine", representation, *seed_parts)
    return mods.fit_joint(ctx.rows.take(index), ctx.layout, components=components, seed=seed, budget=budget,
                          representation=representation), seed


def job_engine_inner(ctx, outer, inner, representation, components):
    """Fit on outer-training minus one inner fold; score that inner fold on every task (selection only)."""
    start = time.time()
    record = {"kind": "engine_inner", "outer": outer, "inner": inner, "representation": representation,
              "components": components}
    train = np.flatnonzero((ctx.outer_of != outer) & (ctx.inner_of[outer] != inner))
    valid = np.flatnonzero(ctx.inner_of[outer] == inner)
    try:
        (state, diagnostics), seed = _fit(ctx, train, representation, components, ("inner", outer, inner, components))
    except (FloatingPointError, ValueError, np.linalg.LinAlgError) as error:
        record.update(failed=repr(error), seconds=time.time() - start)
        return record
    query = ctx.rows.take(valid)
    scores = {}
    for name, task in ctx.tasks.items():
        out = mods.engine_predict(state, task, query)
        truth = target_matrix(task, query)
        if task.kind == "label":
            scores[name] = primary_score(task, out["proba"], truth, query.group, out["classes"],
                                         _seen_mask(truth, out["classes"]))
        else:
            scores[name] = primary_score(task, out["point"], truth, query.group)
    record.update(seed=seed, n_train=len(train), n_valid=len(valid), diagnostics=diagnostics, scores=scores,
                  seconds=time.time() - start)
    return record


def job_baseline_inner(ctx, outer, inner):
    """Score every grid setting of every tuned comparator (the weighting of MEAN and IDENTITY, and the two
    tuned references) on one inner fold (selection only)."""
    start = time.time()
    train = np.flatnonzero((ctx.outer_of != outer) & (ctx.inner_of[outer] != inner))
    valid = np.flatnonzero(ctx.inner_of[outer] == inner)
    ctx.partition.require_training(ctx.rows.group[train], "baseline fit")
    query = ctx.rows.take(valid)
    scores = {}
    for name, task in ctx.tasks.items():
        bank = mods.TaskBaselines(task, ctx.rows.take(train), vote_alpha=ctx.config["knn"]["vote_alpha"],
                                  logistic_max_iter=ctx.config["logistic"]["max_iter"])
        truth = target_matrix(task, query)
        seen = _seen_mask(truth, bank.classes) if task.kind == "label" else None
        scores[name] = {}
        for method in mods.TUNED:
            if not bank.applicable(method):
                continue
            grid = tuned_grid(ctx.config, task, method, bank.label is not None)
            predictions = bank.predict_many(method, query, grid)
            scores[name][method] = [primary_score(task, clip(task, p), truth, query.group, bank.classes, seen)
                                    for p in predictions]
    return {"kind": "baseline_inner", "outer": outer, "inner": inner, "scores": scores, "seconds": time.time() - start}


def clip(task, prediction):
    """The physical support rule, applied to every method alike."""
    return prediction if task.clip_min is None else np.maximum(prediction, task.clip_min)


def state_file(prefix, representation, components):
    return f"{prefix}_{representation}_K{components}.npz"


def job_engine_outer(ctx, outer, representation, components, task_names, save_to=None):
    """Fit on the outer-training samples; predict the held-out samples for the tasks that chose this setting."""
    start = time.time()
    record = {"kind": "engine_outer", "outer": outer, "representation": representation, "components": components,
              "tasks": list(task_names)}
    train = np.flatnonzero(ctx.outer_of != outer)
    test = np.flatnonzero(ctx.outer_of == outer)
    if set(ctx.rows.group[train]) & set(ctx.rows.group[test]):
        raise core.ReconError("a brain sample is on both sides of an outer fold")
    try:
        (state, diagnostics), seed = _fit(ctx, train, representation, components, ("outer", outer, components))
    except (FloatingPointError, ValueError, np.linalg.LinAlgError) as error:
        record.update(failed=repr(error), seconds=time.time() - start)
        return record
    query = ctx.rows.take(test)
    if set(state.groups) & set(query.group.tolist()) or state.units & set(query.unit.tolist()):
        raise core.ReconError("the state of an outer fold contains a held-out brain sample or neuron")
    level = ctx.config["interval_level"]
    predictions = {}
    for name in task_names:
        task = ctx.tasks[name]
        if task.kind == "label":
            out = mods.engine_predict(state, task, query)
            truth = target_matrix(task, query)
            predictions[name] = {"proba": align_classes(out["proba"], out["classes"], ctx.classes),
                                 "seen": _seen_mask(truth, out["classes"])}
        else:
            out = mods.engine_predict(state, task, query, level=level, median=True)
            unknown = np.zeros(len(query), bool)
            for field_name in task.evidence_cat:
                unknown |= out["unknown"].get(field_name, np.zeros(len(query), bool))
            predictions[name] = {"point": out["point"], "median": out["median"], "lower": out["lower"],
                                 "upper": out["upper"], "unknown_label": unknown}
        outside = np.zeros(len(query), bool)             # inputs the API would call out of support
        for field_name in task.evidence_cont:
            low, high = state.ranges[field_name]
            slack = state.support_margin * (high - low)
            outside |= (query.cont[field_name] < low - slack) | (query.cont[field_name] > high + slack)
        predictions[name]["outside_range"] = outside
    if save_to is not None:
        state.save(Path(save_to) / state_file(f"outer{outer}", representation, components))
    record.update(seed=seed, n_train=len(train), n_test=len(test), index=test, diagnostics=diagnostics,
                  state_hash=state.state_hash(), training_groups_sha256=state.groups_sha256, predictions=predictions,
                  state_excludes_every_held_out_sample_and_neuron=True, seconds=time.time() - start)
    return record


def job_baseline_outer(ctx, outer, selected):
    """Every comparator on the held-out samples of one outer fold, with the settings chosen on inner folds."""
    start = time.time()
    train = np.flatnonzero(ctx.outer_of != outer)
    test = np.flatnonzero(ctx.outer_of == outer)
    ctx.partition.require_training(ctx.rows.group[train], "baseline fit")
    query = ctx.rows.take(test)
    level = ctx.config["interval_level"]
    predictions, optimiser = {}, {}
    for name, task in ctx.tasks.items():
        bank = mods.TaskBaselines(task, ctx.rows.take(train), vote_alpha=ctx.config["knn"]["vote_alpha"],
                                  logistic_max_iter=ctx.config["logistic"]["max_iter"])
        entry = {}
        for method in mods.COMPARATORS:
            if not bank.applicable(method):
                continue
            prediction = clip(task, bank.predict(method, query, selected.get(name, {}).get(method)))
            entry[method] = align_classes(prediction, bank.classes, ctx.classes) if task.kind == "label" else prediction
            interval = bank.interval(method, query, level)
            if interval is not None:
                entry[method + "_lower"], entry[method + "_upper"] = clip(task, interval[0]), clip(task, interval[1])
        if task.kind == "label":
            entry["seen"] = _seen_mask(target_matrix(task, query), bank.classes)
            optimiser[name] = bank.logistic_report()
        predictions[name] = entry
    return {"kind": "baseline_outer", "outer": outer, "index": test, "predictions": predictions,
            "logistic": optimiser, "seconds": time.time() - start}


def job_serving(ctx, representation, components, save_to):
    """The state that answers requests about units outside the training partition: all training samples."""
    start = time.time()
    index = np.arange(len(ctx.rows))
    try:
        (state, diagnostics), seed = _fit(ctx, index, representation, components, ("serving", components))
    except (FloatingPointError, ValueError, np.linalg.LinAlgError) as error:
        return {"kind": "engine_serving", "representation": representation, "components": components,
                "failed": repr(error), "seconds": time.time() - start}
    state.save(Path(save_to) / state_file("serving", representation, components))
    return {"kind": "engine_serving", "representation": representation, "components": components, "seed": seed,
            "diagnostics": diagnostics,
            "state_hash": state.state_hash(), "training_groups_sha256": state.groups_sha256,
            "seconds": time.time() - start}


def _dispatch(job):
    name, args = job[0], job[1:]
    return globals()[name](_CTX, *args)


def run_jobs(ctx, jobs, workers, directory):
    """Run jobs in order; with workers > 1, in separate processes that each load the same context."""
    global _CTX
    if workers <= 1 or len(jobs) <= 1:
        _CTX = ctx
        return [_dispatch(job) for job in jobs]
    if __name__ != "__main__":
        raise core.ReconError("parallel workers are available from the command line only")
    import multiprocessing
    with multiprocessing.get_context("spawn").Pool(workers, initializer=_init_worker, initargs=(directory,)) as pool:
        return pool.map(_dispatch, jobs, chunksize=1)


# ----------------------------------------------------------------------------------------------
# the benchmark
# ----------------------------------------------------------------------------------------------


def select_engine(ctx, inner_records):
    """Per task and outer fold: the engine setting (representation, K) with the best mean inner-fold
    primary metric. Ties go to the earlier setting in the declared order (log before hurdle, smaller K)."""
    table = {}
    for record in inner_records:
        table[(record["outer"], record["inner"], record["representation"], record["components"])] = record
    k_grid = ctx.config["engine"]["k_grid"]
    selection = {}
    for name, task in ctx.tasks.items():
        selection[name] = {}
        for outer in range(ctx.config["outer_folds"]):
            scored, ranking = [], []
            for setting in engine_grid(ctx.config):
                cells = [table[(outer, inner, *setting)] for inner in range(ctx.config["inner_folds"])]
                if any("failed" in c for c in cells):
                    continue
                primary = float(np.mean([c["scores"][name][0] for c in cells]))
                tie = float(np.mean([c["scores"][name][1] for c in cells]))
                scored.append((setting, (primary, tie)))
            if not scored:
                raise core.ReconError(f"no engine setting could be fitted on every inner fold of outer fold {outer}")
            remaining = list(scored)
            while remaining:
                best = pick(task, remaining)
                ranking.append(best[0])
                remaining = [entry for entry in remaining if entry[0] != best[0]]
            selection[name][outer] = {"selected": ranking[0], "ranking": ranking,
                                      "inner_scores": {f"{r}:K{k}": s[0] for (r, k), s in scored},
                                      "at_grid_edge": ranking[0][1] in (min(k_grid), max(k_grid))}
    return selection


def select_baselines(ctx, inner_records):
    """Per task, outer fold and tuned comparator: the grid setting with the best mean inner-fold primary metric."""
    by = {(r["outer"], r["inner"]): r for r in inner_records}
    selection = {}
    for name, task in ctx.tasks.items():
        selection[name] = {}
        for outer in range(ctx.config["outer_folds"]):
            selection[name][outer] = {}
            cells = [by[(outer, inner)]["scores"][name] for inner in range(ctx.config["inner_folds"])]
            for method in mods.TUNED:
                if method not in cells[0]:
                    continue
                grid = tuned_grid(ctx.config, task, method, bool(task.evidence_cat))
                scored = [(params, (float(np.mean([c[method][g][0] for c in cells])),
                                    float(np.mean([c[method][g][1] for c in cells])))) for g, params in enumerate(grid)]
                best = pick(task, scored)
                selection[name][outer][method] = {"selected": best[0], "inner_score": best[1][0],
                                                  "grid": [[p, s[0]] for p, s in scored]}
    return selection


def _mode(values, order):
    """Most frequent value; ties go to the earliest in `order`."""
    counts = Counter(json.dumps(v, sort_keys=True) for v in values)
    top = max(counts.values())
    for candidate in order:
        if counts.get(json.dumps(candidate, sort_keys=True), 0) == top:
            return candidate
    raise core.ReconError("no mode found")


def execute(ctx, *, workers=1, directory=None, log=say):
    """Nested whole-sample cross-validation of the engine and every comparator. Returns what `run` saves."""
    out = Path(ctx.world.out)
    states = out / "states"
    config, started = ctx.config, time.time()
    folds, inner = range(config["outer_folds"]), range(config["inner_folds"])
    jobs = [("job_engine_inner", o, i, r, k) for o in folds for i in inner for r, k in engine_grid(config)]
    jobs += [("job_baseline_inner", o, i) for o in folds for i in inner]
    log(f"inner stage: {len(jobs)} jobs")
    records = run_jobs(ctx, jobs, workers, directory)
    engine_inner = [r for r in records if r["kind"] == "engine_inner"]
    baseline_inner = [r for r in records if r["kind"] == "baseline_inner"]
    engine_selection = select_engine(ctx, engine_inner)
    baseline_selection = select_baselines(ctx, baseline_inner)

    n = len(ctx.rows)
    store = {}
    for name, task in ctx.tasks.items():
        width = len(ctx.classes) if task.kind == "label" else len(task.target_fields)
        store[name] = {"fold": ctx.outer_of.copy()}
        for key in ((ENGINE,) if task.kind == "label" else (ENGINE, ENGINE_MEDIAN, "ENGINE_lower", "ENGINE_upper")):
            store[name][key] = np.full((n, width), np.nan)
        store[name]["seen" if task.kind == "label" else "unknown_label"] = np.zeros(n, bool)
        store[name]["outside_range"] = np.zeros(n, bool)

    outer_fits, deviations, used = [], [], {name: {} for name in ctx.tasks}
    pending = {o: {name: list(engine_selection[name][o]["ranking"]) for name in ctx.tasks} for o in folds}
    while any(pending[o] for o in folds):
        jobs = []
        for o in folds:
            wanted = {}
            for name, ranking in pending[o].items():
                wanted.setdefault(tuple(ranking[0]), []).append(name)
            jobs += [("job_engine_outer", o, r, k, names, states) for (r, k), names in sorted(wanted.items())]
        log(f"outer engine stage: {len(jobs)} fits")
        for record in run_jobs(ctx, jobs, workers, directory):
            o = record["outer"]
            if "failed" in record:
                for name in record["tasks"]:
                    pending[o][name].pop(0)
                    if not pending[o][name]:
                        raise core.ReconError(f"no component count could be fitted on outer fold {o} for {name}")
                    deviations.append({"outer": o, "task": name,
                                       "failed_setting": [record["representation"], record["components"]],
                                       "error": record["failed"], "next": pending[o][name][0]})
                outer_fits.append({k: v for k, v in record.items()})
                continue
            index = record.pop("index")
            for name, prediction in record.pop("predictions").items():
                task = ctx.tasks[name]
                if task.kind == "label":
                    store[name][ENGINE][index] = prediction["proba"]
                    store[name]["seen"][index] = prediction["seen"]
                else:
                    store[name][ENGINE][index] = prediction["point"]
                    store[name][ENGINE_MEDIAN][index] = prediction["median"]
                    store[name]["ENGINE_lower"][index] = prediction["lower"]
                    store[name]["ENGINE_upper"][index] = prediction["upper"]
                    store[name]["unknown_label"][index] = prediction["unknown_label"]
                store[name]["outside_range"][index] = prediction["outside_range"]
                used[name][o] = {"representation": record["representation"], "components": record["components"],
                                 "state_hash": record["state_hash"],
                                 "training_groups_sha256": record["training_groups_sha256"]}
                del pending[o][name]
            outer_fits.append(record)

    selected = {o: {name: {m: baseline_selection[name][o][m]["selected"] for m in baseline_selection[name][o]}
                    for name in ctx.tasks} for o in folds}
    log("outer comparator stage")
    baseline_outer = run_jobs(ctx, [("job_baseline_outer", o, selected[o]) for o in folds], workers, directory)
    for record in baseline_outer:
        index = record.pop("index")
        for name, entry in record.pop("predictions").items():
            for key, value in entry.items():
                if key == "seen":
                    if (store[name]["seen"][index] != value).any():
                        raise core.ReconError("engine and comparators disagree on which labels were in training")
                    continue
                if key not in store[name]:
                    store[name][key] = np.full((n, value.shape[1]), np.nan)
                store[name][key][index] = value

    serving = {}
    for name, task in ctx.tasks.items():
        r, k = _mode([[used[name][o]["representation"], used[name][o]["components"]] for o in folds],
                     engine_grid(config))
        entry = {"representation": r, "components": k, "comparators": {}}
        for method in baseline_selection[name][0]:
            grid = tuned_grid(config, task, method, bool(task.evidence_cat))
            entry["comparators"][method] = _mode([selected[o][name][method] for o in folds], grid)
        serving[name] = entry

    for name, task in ctx.tasks.items():
        for key, value in store[name].items():
            if value.dtype.kind == "f" and np.isnan(value).any():
                raise core.ReconError(f"{name}/{key}: a training neuron has no out-of-fold prediction")
        store[name]["truth"] = target_matrix(task, ctx.rows)
        if task.kind == "label":
            store[name]["truth"] = np.array([str(t) for t in store[name]["truth"]])
    fits = {"engine_inner": engine_inner, "engine_outer": outer_fits,
            "baseline_inner_seconds": [r["seconds"] for r in baseline_inner],
            "baseline_outer_seconds": [r["seconds"] for r in baseline_outer],
            "logistic_outer": [{"outer": r["outer"], **r["logistic"]} for r in baseline_outer]}
    selection = {"engine": engine_selection, "comparators": baseline_selection, "engine_used": used,
                 "deviations": deviations, "serving": serving}
    return {"store": store, "fits": fits, "selection": selection, "seconds": time.time() - started}


SERVING_KEY = "all_training_states"      # where run_complete.json holds the record of the all-training states


def all_training_record(ctx, chosen, *, workers=1, directory=None, log=say, save_to=None):
    """One all-training state per chosen engine setting: fitted on every training sample, never scored. It is
    what the API uses for units outside the training partition.

    Returns the record that the COMPLETION RECORD binds: per task the setting, the fingerprint of its state and
    of its training groups, or why the fit failed. The step reads no new data and cannot change a metric, and it
    never stops the benchmark: a numerical failure of a fit is recorded for that fit, and any other failure of
    the step is recorded for every task it left without a state."""
    folder = Path(save_to) if save_to is not None else Path(ctx.world.out) / "states"
    settings = sorted({(s["representation"], s["components"]) for s in chosen.values()})
    log(f"all-training states: {len(settings)} fits")
    try:
        fits, stopped = run_jobs(ctx, [("job_serving", r, k, folder) for r, k in settings], workers, directory), None
    except Exception as error:      # noqa: BLE001  (whatever happened here, the benchmark itself still completes)
        fits, stopped = [], f"the all-training step did not finish: {error!r}"
    record = {"fits": fits, "tasks": {}}
    for name, setting in chosen.items():
        fit = next((f for f in fits if (f["representation"], f["components"]) ==
                    (setting["representation"], setting["components"])), {"failed": stopped})
        record["tasks"][name] = {"representation": setting["representation"], "components": setting["components"],
                                 "state_hash": fit.get("state_hash"), "failed": fit.get("failed"),
                                 "training_groups_sha256": fit.get("training_groups_sha256")}
    return record


def serving_file_content(world):
    """serving.json as the completion record determines it: the bound record of the all-training states, plus the
    SHA-256 of the completion record it is a copy from. The file holds nothing of its own; the reconstruction API
    never reads an expected fingerprint from it."""
    out = Path(world.out)
    complete = core.read_json(out / "run_complete.json")
    if SERVING_KEY not in complete:
        raise core.ReconError("the completion record of this run does not bind its all-training states (it was "
                              "written by an earlier version of this layer); they cannot be served")
    return dict(complete[SERVING_KEY], run_complete_sha256=core.sha_file(out / "run_complete.json"))


def fit_serving(world, *, workers=1, directory=None, log=say, ctx=None):
    """RESTORE the all-training states of a FINISHED run (the `serving` command).

    The run itself fits these states, before it writes its completion record, and that record holds their
    fingerprints. This command fits them again (it reads nothing new) and keeps a state file only if the state is
    exactly the one the run recorded; it writes serving.json again when that copy is missing. It never changes
    what the run recorded: a fit that failed in the run stays failed for that run, and a re-fit that gives
    another state is refused and nothing is written."""
    out = Path(world.out)
    if not (out / "run_complete.json").exists():
        raise core.ReconError("the all-training states are fitted by the run itself; they can be restored only "
                              "after the benchmark has completed")
    complete = verify_run(world)
    if SERVING_KEY not in complete:
        raise core.ReconError("the completion record of this run does not bind its all-training states (it was "
                              "written by an earlier version of this layer); they cannot be restored")
    bound = complete[SERVING_KEY]
    wanted = {name: entry for name, entry in bound["tasks"].items() if not entry["failed"]}
    if wanted:
        ctx = ctx or build_context(world)
        scratch = out / "states" / "restore_in_progress"
        if scratch.exists():
            shutil.rmtree(scratch)
        try:
            again = all_training_record(ctx, wanted, workers=workers, directory=directory, log=log, save_to=scratch)
            for name, entry in wanted.items():
                fitted = again["tasks"][name]
                if fitted["state_hash"] != entry["state_hash"] \
                        or fitted["training_groups_sha256"] != entry["training_groups_sha256"]:
                    raise core.FrozenSourceChanged(
                        f"the re-fitted all-training state of {name} is not the state that the finished run "
                        f"recorded ({fitted['failed'] or fitted['state_hash']}); nothing was written")
            for file in sorted(scratch.glob("*.npz")):
                os.replace(file, out / "states" / file.name)
        finally:
            shutil.rmtree(scratch, ignore_errors=True)
    if not (out / "serving.json").exists():
        core.write_json(out / "serving.json", serving_file_content(world), exclusive=True)
    for name, entry in bound["tasks"].items():
        if entry["failed"]:
            log(f"{name}: the all-training fit failed in the run ({entry['failed']}); it stays failed for this run")
    return core.read_json(out / "serving.json")


def save_predictions(ctx, result):
    """Write every out-of-fold prediction, the fit records and the selection. Paths in the manifest are
    relative to the output folder, so a copied folder verifies as well as the original."""
    out = Path(ctx.world.out)
    (out / "predictions").mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name, entry in result["store"].items():
        relative = f"predictions/{name}.npz"
        np.savez_compressed(out / relative, neuron=ctx.table.neuron, sample=ctx.table.sample, label=ctx.table.label,
                            classes=np.array(ctx.classes), **entry)
        manifest[name] = {"path": relative, "sha256": core.sha_file(out / relative), "arrays": sorted(entry)}
    core.write_json(out / "fits.json", result["fits"])
    core.write_json(out / "selection.json", result["selection"])
    core.write_json(out / "prediction_manifest.json", manifest)
    return manifest


# ----------------------------------------------------------------------------------------------
# metrics from saved predictions
# ----------------------------------------------------------------------------------------------


def _summary(per_sample):
    return {metric: {"value": core.over_samples(values), "samples": len(values)} for metric, values in per_sample.items()}


def _advantage(metric, engine, comparator):
    """Positive when the engine is better on this metric."""
    return engine - comparator if metric in HIGHER_IS_BETTER else comparator - engine


def compare(metric, engine, comparator, seed, config):
    """Paired brain-sample differences (engine minus comparator) with a descriptive bootstrap interval.

    Three ways of being 'better', all reported:
      strict    the engine's mean over brain samples is better, by any amount;
      interval  the descriptive interval of the paired difference excludes zero in the engine's favour;
      clear     the interval rule AND a relative gain of at least config['min_relative_gain'] over the
                comparator's own mean (the declared minimum effect).
    Paired differences smaller than config['tie_tolerance'] (relative to the metric) are exact ties: two
    predictors that agree to rounding error are not better than one another."""
    common = sorted(set(engine) & set(comparator))
    e = np.array([engine[s] for s in common], dtype=float)
    c = np.array([comparator[s] for s in common], dtype=float)
    tolerance = float(config.get("tie_tolerance", 0.0)) * np.maximum(1.0, np.maximum(np.abs(e), np.abs(c)))
    difference = np.where(np.abs(e - c) <= tolerance, 0.0, e - c)
    boot = core.paired_bootstrap(difference, seed, samples=config["bootstrap"]["resamples"],
                                 level=config["bootstrap"]["level"])
    sign = 1.0 if metric in HIGHER_IS_BETTER else -1.0
    advantage = sign * difference
    minimum = float(config.get("min_relative_gain", 0.0))
    if not common:
        return {"samples": 0, "engine_minus_comparator": boot, "advantage_mean": float("nan"),
                "advantage_interval": [float("nan"), float("nan")], "samples_favouring_engine": 0,
                "relative_gain": float("nan"), "min_relative_gain": minimum, "engine_better_strict": False,
                "engine_better_interval": False, "engine_better_clear": False, "label": "descriptive"}
    low, high = sorted((sign * boot["lower"], sign * boot["upper"]))
    reference = abs(float(c.mean()))
    relative = float(advantage.mean() / reference) if reference > 0 else 0.0
    return {"samples": len(common), "engine_minus_comparator": boot, "advantage_mean": float(advantage.mean()),
            "advantage_interval": [low, high], "samples_favouring_engine": int((advantage > 0).sum()),
            "relative_gain": relative, "min_relative_gain": minimum,
            "engine_better_strict": bool(advantage.mean() > 0), "engine_better_interval": bool(low > 0),
            "engine_better_clear": bool(low > 0 and relative >= minimum), "label": "descriptive"}


def _verdict(flags):
    flags = list(flags)
    if not flags:
        return "not applicable"
    return "every" if all(flags) else ("some" if any(flags) else "none")


def score_task(name, task, data, config, salt):
    truth, sample, label = data["truth"], data["sample"], data["label"]
    classes = [str(c) for c in data["classes"]]
    seen = data["seen"] if task.kind == "label" else None
    methods = [ENGINE] + ([ENGINE_MEDIAN] if task.kind != "label" else []) + \
        [m for m in mods.COMPARATORS if m in data]
    per = {m: sample_metrics(task, neuron_scores(task, data[m], truth, classes, seen), sample) for m in methods}
    result = {"task": name, "title": task.title, "kind": task.kind, "primary_metric": task.primary,
              "higher_is_better": task.higher_is_better, "neurons": int(len(sample)),
              "samples": int(len(set(sample.tolist()))),
              "methods": {m: _summary(per[m]) for m in methods},
              "per_sample": {m: per[m] for m in methods}}
    if task.kind == "profile":
        constant = int((np.asarray(truth).std(axis=1) == 0).sum())
        result["constant_truth_profiles_excluded_from_correlation_only"] = constant
    if task.kind == "label":
        unseen_samples = sorted(set(sample[~seen].tolist()))
        result["unseen_label"] = {"neurons": int((~seen).sum()), "samples": len(unseen_samples),
                                  "rule": "the true label had no training sample in that fold; unscoreable for "
                                          "every method and excluded from accuracy and log-loss for all of them"}
    else:
        result["rows_with_label_unseen_in_training"] = int(data["unknown_label"].sum())
    result["rows_with_an_input_outside_the_training_range"] = int(data["outside_range"].sum())
    comparators = [m for m in mods.COMPARATORS if m in data]
    result["comparisons"] = {
        metric: {c: compare(metric, per[ENGINE][metric], per[c][metric],
                            core.stable_seed(salt, "bootstrap", name, metric, c), config) for c in comparators}
        for metric in METRICS[task.kind]}
    if task.kind == "profile":
        result["comparisons"]["mae_median_estimate"] = {
            c: compare("mae", per[ENGINE_MEDIAN]["mae"], per[c]["mae"],
                       core.stable_seed(salt, "bootstrap", name, "mae_median", c), config) for c in comparators}
    primary = result["comparisons"][task.primary]
    always = [c for c in mods.ALWAYS_REPORTED if c in primary]
    result["verdict"] = {
        "primary_metric": task.primary,
        "always_reported_baselines": _verdict(primary[c]["engine_better_strict"] for c in always),
        "all_comparators": _verdict(primary[c]["engine_better_strict"] for c in comparators),
        "all_comparators_interval": _verdict(primary[c]["engine_better_interval"] for c in comparators),
        "always_reported_baselines_clear": _verdict(primary[c]["engine_better_clear"] for c in always),
        "all_comparators_clear": _verdict(primary[c]["engine_better_clear"] for c in comparators),
        "min_relative_gain": float(config.get("min_relative_gain", 0.0)),
        "by_metric": {metric: _verdict(v[c]["engine_better_strict"] for c in comparators)
                      for metric, v in result["comparisons"].items()},
        "rule": "strict: the engine's mean over brain samples is better, by any amount; interval: the descriptive "
                "95% bootstrap interval of the paired difference also excludes zero in the engine's favour; clear: "
                "the interval rule and a relative gain of at least the declared minimum over that comparator"}
    best = None
    for c in comparators:
        value = result["methods"][c][task.primary]["value"]
        if best is None or (value > best[1] if task.higher_is_better else value < best[1]):
            best = (c, value)
    result["best_comparator"] = {"name": best[0], "value": best[1],
                                 "relative_gain_of_engine": primary[best[0]]["relative_gain"]}

    # interval coverage: per neuron fraction of target coordinates covered, per sample, over samples
    level = config["interval_level"]
    coverage = {}
    if task.kind != "label":
        for m in [ENGINE] + [c for c in ("MEAN", "IDENTITY") if c + "_lower" in data]:
            lower, upper = data[m + "_lower"], data[m + "_upper"]
            inside = (truth >= lower) & (truth <= upper)
            entry = {"level": level, "empirical": core.over_samples(core.by_sample(inside.mean(axis=1), sample)),
                     "mean_width": core.over_samples(core.by_sample((upper - lower).mean(axis=1), sample)),
                     "per_coordinate": {c: core.over_samples(core.by_sample(inside[:, j].astype(float), sample))
                                        for j, c in enumerate(task.target_coordinates)}}
            if task.kind == "profile":
                nonzero = truth > 0
                share = np.where(nonzero.any(axis=1), (inside & nonzero).sum(axis=1) / np.maximum(nonzero.sum(axis=1), 1),
                                 np.nan)
                entry["share_of_nonzero_truths_inside"] = core.over_samples(core.by_sample(share, sample))
            coverage[m] = entry
        tolerance = config["coverage_tolerance"]
        for m, entry in coverage.items():
            if task.kind == "profile":
                # the truth has a mass at exactly zero, so a central interval cannot hit its level exactly:
                # it must reach the level (within the tolerance) and be no wider than the unconditional interval
                reference = coverage["MEAN"]["mean_width"] if "MEAN" in coverage else float("inf")
                entry["rule"] = "coverage >= level - tolerance and mean width <= width of the MEAN interval"
                entry["validated"] = bool(entry["empirical"] >= level - tolerance and entry["mean_width"] <= reference)
            else:
                entry["rule"] = "|coverage - level| <= tolerance"
                entry["validated"] = bool(abs(entry["empirical"] - level) <= tolerance)
    else:
        probability = data[ENGINE]
        order = np.argsort(-probability, axis=1, kind="stable")
        cumulative = np.cumsum(np.take_along_axis(probability, order, axis=1), axis=1)
        size = (cumulative < level).sum(axis=1) + 1
        position = {c: i for i, c in enumerate(classes)}
        truth_index = np.array([position[t] for t in truth])
        rank = np.argmax(order == truth_index[:, None], axis=1)
        inside = np.where(seen, (rank < size).astype(float), np.nan)
        entry = {"level": level, "empirical": core.over_samples(core.by_sample(inside, sample)),
                 "mean_set_size": core.over_samples(core.by_sample(np.where(seen, size, np.nan), sample)),
                 "kind": "smallest label set holding at least the level"}
        entry["rule"] = "coverage >= level - tolerance (a label set cannot hit a level exactly)"
        entry["validated"] = bool(entry["empirical"] >= level - config["coverage_tolerance"])
        coverage[ENGINE] = entry
    result["coverage"] = coverage

    # by source label
    by_label = {}
    for value in sorted(set(label.tolist())):
        members = sorted(set(sample[label == value].tolist()))
        entry = {"samples": len(members), "neurons": int((label == value).sum()),
                 "descriptive_only": len(members) < config["min_samples_for_interval"], "methods": {}}
        for m in methods:
            entry["methods"][m] = {metric: core.over_samples({s: v for s, v in per[m][metric].items() if s in members})
                                   for metric in per[m]}
            entry["scored_samples"] = len([s for s in members if s in per[m][task.primary]])
        entry["engine_better_than"] = {
            c: bool(_advantage(task.primary, entry["methods"][ENGINE][task.primary],
                               entry["methods"][c][task.primary]) > 0) for c in comparators
            if not np.isnan(entry["methods"][ENGINE][task.primary])}
        if not entry["descriptive_only"]:
            entry["comparisons"] = {
                c: compare(task.primary, {s: v for s, v in per[ENGINE][task.primary].items() if s in members},
                           {s: v for s, v in per[c][task.primary].items() if s in members},
                           core.stable_seed(salt, "bootstrap", name, "label", value, c), config) for c in comparators}
        by_label[value] = entry
    result["by_source_label"] = by_label
    return result


def score(world, log=say):
    """Metrics from the saved predictions only (no table is read).

    The FIRST scoring of a run writes metrics.json and evidence.json, and they are never overwritten. Any
    later call is a re-score: its metrics go to rescores/<n>/metrics.json with a note of the code that made
    them, and nothing that the reconstruction API or the run's report reads is changed."""
    out = Path(world.out)
    declaration = core.read_json(out / "pre_run_declaration.json")
    config, frozen_code = declaration["config"], declaration["code_sha256"]      # the frozen settings, always
    registry = core.Registry.from_contract(ROOT / CONTRACT)
    layout = mods.jiao_layout(registry, world.regions)
    tasks = mods.build_tasks(layout, registry, world.t4_observed)
    manifest = core.read_json(out / "prediction_manifest.json")
    current_code = code_hashes()
    metrics = {"schema": core.SCHEMA, "world": world.name, "synthetic": world.synthetic, "exploratory": True,
               "exposure": EXPOSURE, "unit_of_independence": "brain sample (ID prefix); not asserted to be an animal",
               "code_files_changed_since_freeze": sorted(k for k in frozen_code if frozen_code[k] != current_code.get(k)),
               "tasks": {}}
    for name in config["tasks"]:
        entry = manifest[name]
        if core.sha_file(out / entry["path"]) != entry["sha256"]:
            raise core.FrozenSourceChanged(f"saved predictions of {name} changed after the run")
        with np.load(out / entry["path"], allow_pickle=False) as data:
            arrays = {k: data[k] for k in data.files}
        metrics["tasks"][name] = score_task(name, tasks[name], arrays, config, config["salt"])
        v = metrics["tasks"][name]["verdict"]
        log(f"{name}: engine beats every/some/none of the always-reported baselines: {v['always_reported_baselines']}; "
            f"of all comparators: {v['all_comparators']} (interval rule: {v['all_comparators_interval']}; with the "
            f"minimum gain: {v['all_comparators_clear']})")
    if not (out / "metrics.json").exists():
        core.write_json(out / "metrics.json", metrics, exclusive=True)
        write_evidence(world, metrics, config)
        return metrics
    number = 1 + len([d for d in (out / "rescores").glob("[0-9][0-9][0-9]")]) if (out / "rescores").exists() else 1
    metrics["rescore"] = {"number": number, "utc": utc_now(),
                          "note": "a repeat of the scoring from the saved predictions; the run's own metrics.json "
                                  "and evidence.json are unchanged"}
    core.write_json(out / "rescores" / f"{number:03d}" / "metrics.json", metrics, exclusive=True)
    log(f"re-score {number} written to {out / 'rescores' / f'{number:03d}'}; the run's own metrics are unchanged")
    return metrics


def write_evidence(world, metrics, config):
    """What the benchmark established, in the form the reconstruction API consumes."""
    out = Path(world.out)
    selection = core.read_json(out / "selection.json")
    declaration = out / "pre_run_declaration.json"
    pointer = {"metrics": "metrics.json", "metrics_sha256": core.sha_file(out / "metrics.json"),
               "pre_run_declaration_sha256": core.sha_file(declaration) if declaration.exists() else None,
               "kind": "exploratory nested whole-brain-sample cross-validation inside the training partition"}
    evidence = {"schema": core.SCHEMA, "world": world.name, "exploratory": True, "confirmatory": False,
                "exposure": EXPOSURE, "benchmark": pointer, "tasks": {}}
    for name, result in metrics["tasks"].items():
        # the ONLY decision of the benchmark: better than every comparator, each interval excluding zero,
        # and by at least the declared minimum relative gain
        supported = result["verdict"]["all_comparators_clear"] == "every"
        best = result["best_comparator"]["name"]
        evidence["tasks"][name] = {
            "conditional": f"{name}_{result['title']}",
            "status": "exploratory_supported" if supported else "no_better_than_baseline",
            "status_rule": "better than every comparator on the primary metric, every descriptive 95% interval "
                           "excluding zero, and a relative gain of at least min_relative_gain over each",
            "min_relative_gain": float(config.get("min_relative_gain", 0.0)),
            "relative_gain_over_best_comparator": result["best_comparator"]["relative_gain_of_engine"],
            "primary_metric": result["primary_metric"],
            "engine": result["methods"][ENGINE][result["primary_metric"]]["value"],
            "comparators": {c: result["methods"][c][result["primary_metric"]]["value"]
                            for c in mods.COMPARATORS if c in result["methods"]},
            "verdict": result["verdict"],
            "best_comparator": best,
            "best_comparator_params": selection["serving"][name]["comparators"].get(best),
            "coverage": {"engine": result["coverage"].get(ENGINE),
                         "baselines": {c: v for c, v in result["coverage"].items() if c != ENGINE}},
            "serving": selection["serving"][name],
        }
    core.write_json(out / "evidence.json", evidence, exclusive=True)
    return evidence


def verify_run(world):
    """A FINISHED run: its completion record exists and every file it binds still has the recorded SHA-256
    (the declaration, the prediction manifest, the fit records, the selection, the metrics and the evidence).

    The record of the all-training states is INSIDE the completion record (SERVING_KEY). serving.json is a copy
    of it; when that file exists it must be exactly the copy, so an edited serving.json makes the run fail this
    check instead of passing next to an unchanged completion record."""
    out = Path(world.out)
    if not (out / "run_complete.json").exists():
        raise core.ReconError("no finished benchmark exists for this world (run_complete.json is missing)")
    record = core.read_json(out / "run_complete.json")
    for name, sha in record["sha256"].items():
        if not (out / name).exists() or core.sha_file(out / name) != sha:
            raise core.FrozenSourceChanged(f"{name} is not the file that the finished run recorded")
    if (out / "serving.json").exists():
        if SERVING_KEY not in record:
            raise core.FrozenSourceChanged("serving.json exists, but the completion record of this run binds no "
                                           "all-training states")
        try:
            copy = core.canonical(core.read_json(out / "serving.json"))
        except ValueError:
            copy = None
        if copy != core.canonical(serving_file_content(world)):
            raise core.FrozenSourceChanged("serving.json is not the record of the all-training states that the "
                                           "finished run bound in its completion record")
    return record


def evidence_objects(world):
    """evidence.json -> {conditional name: ValidationEvidence}. The file decides whether the API returns the
    engine's answer, a baseline-level answer or a refusal, so it is read only when it is the file the finished
    run recorded."""
    out = Path(world.out)
    verify_run(world)
    evidence = core.read_json(out / "evidence.json")
    if evidence["benchmark"]["metrics_sha256"] != core.sha_file(out / "metrics.json"):
        raise core.FrozenSourceChanged("evidence.json does not point to the metrics of this run")
    result = {}
    for name, e in evidence["tasks"].items():
        engine_coverage = e["coverage"].get("engine") or {}
        result[e["conditional"]] = core.ValidationEvidence(
            status=e["status"], task=name, benchmark=evidence["benchmark"], primary_metric=e["primary_metric"],
            coverage_validated=bool(engine_coverage.get("validated", False)), coverage=e["coverage"],
            best_comparator=e["best_comparator"], best_comparator_params=e["best_comparator_params"],
            note="exploratory benchmark inside the training partition; not confirmatory")
    return result


# ----------------------------------------------------------------------------------------------
# freeze, run once, status
# ----------------------------------------------------------------------------------------------


def build_declaration(world, reviewed, frozen_utc):
    """Everything the single run is bound to. Built from inputs only: no outcome value is read."""
    config = config_for(world)
    partition = core.TrainingPartition(world)
    table = core.load_training_table(world, with_outcomes=False, partition=partition)
    registry = core.Registry.from_contract(ROOT / CONTRACT)
    layout = mods.jiao_layout(registry, world.regions)
    tasks = {n: t for n, t in mods.build_tasks(layout, registry, world.t4_observed).items() if n in config["tasks"]}
    splits = make_splits(table, config)
    protocol = ROOT / PROTOCOL
    declaration = {
        "schema": core.SCHEMA, "status": "frozen_before_the_single_run", "world": world.name,
        "synthetic": world.synthetic, "frozen_utc": frozen_utc, "reviewed": reviewed,
        "exploratory": True, "confirmatory_claim": False, "exposure": EXPOSURE,
        "unit_of_independence": "brain sample (ID prefix). The paper reports 166 mice and the data hold 166 sample "
                                "IDs; no statement was found that no mouse gave two samples. Every sample carries "
                                "one source label, so a label effect is confounded with sample and batch.",
        "code_sha256": code_hashes(),
        # the protocol binds the real cohort; a synthetic world is a check of the code and does not pin it
        "protocol_sha256": None if world.synthetic else (core.sha_file(protocol) if protocol.exists() else None),
        "state_reader": {"state_format": mods.STATE_FORMAT, "support_margin": mods.SUPPORT_MARGIN,
                         "quantile_iterations": mods.QUANTILE_ITERATIONS},
        "input_sha256": {"manifest": partition.manifest_sha256, "soma_csv": core.sha_file(world.soma_csv),
                         "soma_xlsx": core.sha_file(world.soma_xlsx), "regional_csv": core.sha_file(world.regional_csv),
                         "data_contract": core.sha_file(ROOT / CONTRACT)},
        "training": {"samples": len(partition.train_samples), "neurons": len(partition.train_neurons),
                     "samples_sha256": partition.samples_sha256, "neurons_sha256": core.ids_sha(table.neuron.tolist()),
                     "inputs_sha256": core.array_sha(table.neuron, table.sample, table.label, table.soma)},
        "config": config, "tasks": task_declaration(tasks),
        "comparators": {"always_reported": list(mods.ALWAYS_REPORTED), "extra_references": list(mods.EXTRA_REFERENCES)},
        "regions": list(world.regions), "t4_observed_block": list(world.t4_observed),
        "splits": splits, "splits_sha256": core.digest(splits),
        "runtime": runtime(),
        "aborted_starts_before_this_freeze": aborted_starts(world),
    }
    declaration["supersedes"] = superseded(world, declaration)
    return declaration


# What a benchmark IS: none of this may change between a started run that was set aside and the freeze after it.
DESIGN_KEYS = ("world", "synthetic", "config", "tasks", "comparators", "regions", "t4_observed_block", "splits",
               "splits_sha256", "input_sha256", "training", "state_reader")
# Files that exist only once outcomes have been used (predictions, scores, selection on scores).
OUTCOME_FILES = ("predictions", "prediction_manifest.json", "predictions_complete.json", "fits.json", "selection.json",
                 "metrics.json", "evidence.json", "run_complete.json", "serving.json", "rescores")
# Files that a finished run binds by SHA-256 in run_complete.json. The record of the all-training states is not
# among them because run_complete.json HOLDS it (SERVING_KEY); serving.json is a copy of that record, written
# after the completion record, and verify_run compares the copy with it whenever the file exists.
BOUND_FILES = ("pre_run_declaration.json", "prediction_manifest.json", "predictions_complete.json", "fits.json",
               "selection.json", "metrics.json", "evidence.json")


def aborted_starts(world):
    """Earlier freezes and starts that were set aside before any outcome-derived file existed."""
    folder = Path(world.out) / "aborted"
    if not folder.exists():
        return []
    return [dict(core.read_json(p), folder=p.parent.name) for p in sorted(folder.glob("*/abort.json"))]


def superseded(world, declaration):
    """The last STARTED run that was set aside, and what differs from its declaration. None when no started run
    was set aside (a freeze that was withdrawn before any run binds nothing)."""
    started = [a for a in aborted_starts(world) if a.get("run_had_started")]
    if not started:
        return None
    last = started[-1]
    old = core.read_json(Path(world.out) / "aborted" / last["folder"] / "pre_run_declaration.json")
    design = sorted(k for k in DESIGN_KEYS if core.canonical(old.get(k)) != core.canonical(declaration.get(k)))
    old_code, new_code = old.get("code_sha256", {}), declaration["code_sha256"]
    files = sorted(k for k in set(old_code) | set(new_code) if old_code.get(k) != new_code.get(k))
    if old.get("protocol_sha256") != declaration.get("protocol_sha256"):
        files.append(PROTOCOL)
    return {"pre_run_declaration_sha256": last.get("pre_run_declaration_sha256"), "aborted_utc": last["aborted_utc"],
            "reason": last["reason"], "design_differences": design, "changed_files": files}


def abort(world, reason, log=say):
    """Set aside a freeze that was never run, or a start that broke for a TECHNICAL reason before any
    outcome-derived file was written.

    Refused as soon as a file that holds predictions or scores exists (OUTCOME_FILES: the predictions and their
    manifest, fits.json with its inner-fold scores, selection.json, the metrics, the evidence, the completion
    record): from then on outcomes have been used and the run stands; `score` finishes an interrupted run from
    its saved predictions. The declaration, the start marker and any fitted but unscored states are moved to
    aborted/<time>/ and listed in the next pre-run declaration. Nothing is deleted. After a STARTED run was set
    aside, the next freeze must keep its whole design (see `freeze`).
    """
    out = Path(world.out)
    if not reason or not reason.strip():
        raise core.ReconError("abort needs --reason")
    if not (out / "run_started.json").exists() and not (out / "pre_run_declaration.json").exists():
        raise core.ReconError("there is no freeze and no started run to set aside")
    present = [name for name in OUTCOME_FILES if (out / name).exists()]
    if present:
        raise core.ReconError("outcome-derived files exist (" + ", ".join(present) + "): this run cannot be set "
                              "aside; its results stand")
    stamp = utc_now().replace(":", "").replace("+0000", "Z")
    target = out / "aborted" / stamp
    target.mkdir(parents=True)
    declaration = out / "pre_run_declaration.json"
    declaration_sha = core.sha_file(declaration) if declaration.exists() else None
    moved = []
    for name in ("pre_run_declaration.json", "run_started.json", "states"):
        if (out / name).exists():
            (out / name).rename(target / name)
            moved.append(name)
    core.write_json(target / "abort.json", {"aborted_utc": utc_now(), "reason": reason.strip(), "moved": moved,
                                            "run_had_started": "run_started.json" in moved,
                                            "outcome_files_present": present,
                                            "fitted_but_unscored_states_existed": "states" in moved,
                                            "pre_run_declaration_sha256": declaration_sha}, pretty=True)
    log(f"start set aside: {target}")
    return target


def freeze(world, reviewed, log=say):
    """Write the pre-run declaration (exclusive create).

    After a STARTED run was set aside, a new freeze must be the SAME benchmark: tasks, comparators, grids,
    metrics, salt, folds, inputs and reader settings (DESIGN_KEYS) are compared with the declaration that was set
    aside and any difference is refused. Only code and protocol text may differ (the repair); the new
    declaration records the old one's SHA-256 and the files that changed."""
    out = Path(world.out)
    if not reviewed or not reviewed.strip():
        raise core.ReconError("freeze needs --reviewed: who reviewed the draft protocol and code, and where")
    if not world.synthetic and (not (ROOT / PROTOCOL).exists() or None in code_hashes().values()):
        raise core.ReconError("the protocol or a code file is missing; the real cohort cannot be frozen")
    declaration = build_declaration(world, reviewed.strip(), utc_now())
    previous = declaration["supersedes"]
    if previous and previous["design_differences"]:
        raise core.ReconError("a started run was set aside; the next freeze must keep its design, but these parts "
                              "differ: " + ", ".join(previous["design_differences"]) + ". A different design is a "
                              "new benchmark, reported next to the first, not a restart of it.")
    path = out / "pre_run_declaration.json"
    out.mkdir(parents=True, exist_ok=True)
    core.write_json(path, declaration, exclusive=True, pretty=True)
    log(f"frozen: {path}  sha256 {core.sha_file(path)}")
    log(f"  {declaration['training']['samples']} brain samples, {declaration['training']['neurons']} neurons; tasks "
        f"{', '.join(declaration['tasks'])}; engine settings {len(engine_grid(declaration['config']))}; "
        f"folds {declaration['splits_sha256'][:12]}")
    return declaration


def verify_freeze(world):
    """The declaration exists and still describes this code, protocol, inputs, configuration and split."""
    path = Path(world.out) / "pre_run_declaration.json"
    if not path.exists():
        raise core.ReconError(f"not frozen: {path} does not exist. The benchmark runs only after the freeze.")
    frozen = core.read_json(path)
    current = build_declaration(world, frozen["reviewed"], frozen["frozen_utc"])
    differing = sorted(k for k in set(frozen) | set(current)
                       if core.canonical(frozen.get(k)) != core.canonical(current.get(k)))
    if differing:
        raise core.FrozenSourceChanged("the pre-run declaration no longer matches: " + ", ".join(differing))
    return frozen


def finish(world, *, how, seconds=None, log=say, ctx=None, workers=1, directory=None):
    """Finish a run whose predictions are saved: fit the all-training states, score, and write the completion
    record (once).

    Order: (1) the all-training states are fitted from the saved selection, before this call computes any
    metric (they read no new data, are never scored and cannot change a metric); (2) the metrics and the
    evidence are computed from the saved predictions; (3) run_complete.json is written. It binds by SHA-256 the
    declaration, the prediction manifest, the fit records, the selection, the metrics and the evidence, and it
    HOLDS the record of the all-training states: their settings, fingerprints and training groups. So every
    expected fingerprint that the reconstruction API uses is behind the one anchor, the SHA-256 of
    run_complete.json; (4) serving.json is written as a copy of that record."""
    out = Path(world.out)
    if (out / "run_complete.json").exists():
        raise core.ReconError("this run is already complete")
    if not (out / "predictions_complete.json").exists():
        raise core.ReconError("the predictions of this run were not saved completely; there is nothing to finish")
    if (out / "serving.json").exists():
        raise core.ReconError("serving.json exists although the run is not complete: it is written only as a copy "
                              "of the completion record. Move it away and finish again; nothing was changed")
    selection = core.read_json(out / "selection.json")
    serving = all_training_record(ctx or build_context(world), selection["serving"], workers=workers,
                                  directory=directory, log=log)
    if (out / "metrics.json").exists():             # interrupted between the first scoring and this record
        metrics = core.read_json(out / "metrics.json")
        if not (out / "evidence.json").exists():
            write_evidence(world, metrics, core.read_json(out / "pre_run_declaration.json")["config"])
    else:
        metrics = score(world, log=log)
    hashes = {name: core.sha_file(out / name) for name in BOUND_FILES}
    core.write_json(out / "run_complete.json",
                    {"completed_utc": utc_now(), "completed_by": how, "seconds": seconds, "sha256": hashes,
                     "pre_run_declaration_sha256": hashes["pre_run_declaration.json"],
                     "metrics_sha256": hashes["metrics.json"], "evidence_sha256": hashes["evidence.json"],
                     "deviations": selection["deviations"], SERVING_KEY: serving}, exclusive=True)
    # a copy of what the completion record holds, for readers of the file; never a source of an expected fingerprint
    core.write_json(out / "serving.json", serving_file_content(world), exclusive=True)
    return metrics


def run(world, *, workers=1, directory=None, log=say, post=True):
    """The single run: verify the freeze, write the exclusive start marker, fit, predict, save, fit the
    all-training states, score, and write the completion record that binds all of it.
    `post=False` skips the two steps that follow a finished benchmark (demo, report)."""
    out = Path(world.out)
    frozen = verify_freeze(world)
    declaration_sha = core.sha_file(out / "pre_run_declaration.json")
    try:
        core.write_json(out / "run_started.json", {"started_utc": utc_now(), "pre_run_declaration_sha256": declaration_sha,
                                                    "workers": workers, "platform": platform.platform()},
                        exclusive=True)
    except FileExistsError:
        raise core.ReconError("this benchmark has already been started; it runs once and is not re-run") from None
    ctx = build_context(world, frozen["config"])
    if core.digest(ctx.splits) != frozen["splits_sha256"]:
        raise core.FrozenSourceChanged("the folds differ from the frozen ones")
    result = execute(ctx, workers=workers, directory=directory, log=log)
    manifest = save_predictions(ctx, result)
    core.write_json(out / "predictions_complete.json", {"utc": utc_now(), "predictions": manifest,
                                                         "table_sha256": ctx.table.sha256(),
                                                         "training_samples_sha256": ctx.partition.samples_sha256})
    metrics = finish(world, how="run", seconds=result["seconds"], log=log, ctx=ctx, workers=workers,
                     directory=directory)
    log(f"run complete in {result['seconds']:.0f} s; run_complete.json sha256 "
        f"{core.sha_file(out / 'run_complete.json')}")
    # The benchmark is complete and stands. What follows reads no new data and can be repeated by hand.
    for name, step in (("demo", lambda: demo(world, log=log)), ("report", lambda: report(world, log=log))):
        if not post:
            break
        try:
            step()
        except Exception as error:
            log(f"{name} did not finish ({error!r}); repeat it with the `{name}` command")
    return metrics


def status(world):
    out = Path(world.out)
    info = {"world": world.name, "out": str(out)}
    for name in ("pre_run_declaration.json", "run_started.json", "predictions_complete.json", "run_complete.json",
                 "metrics.json", "evidence.json", "serving.json", "api_demo.json"):
        info[name] = (out / name).exists()
    if info["pre_run_declaration.json"]:
        try:
            verify_freeze(world)
            info["freeze_matches_current_files"] = True
        except core.ReconError as error:
            info["freeze_matches_current_files"] = str(error)
    if info["run_complete.json"]:
        # the anchor: compare it with the SHA-256 that `run` printed and that the report states in its first lines
        info["run_complete_sha256"] = core.sha_file(out / "run_complete.json")
        try:
            verify_run(world)
            info["run_files_match_the_completion_record"] = True
        except core.ReconError as error:
            info["run_files_match_the_completion_record"] = str(error)
    return info


# ----------------------------------------------------------------------------------------------
# input audit and folds (inputs only)
# ----------------------------------------------------------------------------------------------


def input_audit(world, log=say):
    """Training samples only. Soma and labels are inputs; of the regional table only the header and
    per-region MARGINAL summaries are computed. Nothing relates an input to a regional value."""
    partition = core.TrainingPartition(world)
    table = core.load_training_table(world, with_outcomes=True, partition=partition)
    config = config_for(world)
    splits = make_splits(table, config)
    label_of = dict(zip(table.sample.tolist(), table.label.tolist()))
    labels = {}
    for value in sorted(set(table.label.tolist())):
        members = table.label == value
        labels[value] = {"samples": len(set(table.sample[members].tolist())), "neurons": int(members.sum())}
    per_sample = Counter(table.sample.tolist())
    profile = table.profile()
    audit = {
        "schema": core.SCHEMA, "world": world.name, "created_utc": utc_now(),
        "scope": "training brain samples only; inputs, and per-region marginal summaries of the regional table; "
                 "no input was related to any regional value and no model was fitted",
        "sources": table.audit,
        "training": {"samples": len(partition.train_samples), "neurons": len(table),
                     "neurons_per_sample": {"min": min(per_sample.values()), "median": float(np.median(list(per_sample.values()))),
                                            "max": max(per_sample.values())},
                     "one_source_label_per_sample": True},
        "source_labels": labels,
        "soma_um": {axis: {"min": float(table.soma[:, j].min()), "max": float(table.soma[:, j].max()),
                           "mean": float(table.soma[:, j].mean()), "sd": float(table.soma[:, j].std())}
                    for j, axis in enumerate(core.SOMA_AXES)},
        "exact_duplicate_soma_rows": int(len(table) - len({tuple(r) for r in table.soma.tolist()})),
        "fine_soma_labels": len(set(table.fine_label.tolist())),
        "regional_marginals": {
            "regions": len(table.regions), "all_zero_profiles": int((table.length.sum(axis=1) == 0).sum()),
            "zero_fraction_overall": float((table.length == 0).mean()),
            "per_region": {r: {"zero_fraction": float((table.length[:, j] == 0).mean()),
                               "max_published": float(table.length[:, j].max()),
                               "mean_log10_1p": float(profile[:, j].mean()), "sd_log10_1p": float(profile[:, j].std())}
                           for j, r in enumerate(table.regions)}},
        "folds": {"salt": config["salt"],
                  "outer": {fold: {"samples": len(held), "neurons": int(np.isin(table.sample, held).sum()),
                                   "labels": dict(sorted(Counter(label_of[s] for s in held).items()))}
                            for fold, held in splits["outer"].items()},
                  "splits_sha256": core.digest(splits)},
    }
    path = Path(world.out) / "input_audit.json"
    core.write_json(path, audit, pretty=True)
    log(f"input audit written: {path}")
    log(f"training partition: {audit['training']['samples']} brain samples, {audit['training']['neurons']} neurons")
    for fold, entry in audit["folds"]["outer"].items():
        log(f"  outer fold {fold}: {entry['samples']} samples, {entry['neurons']} neurons")
    return audit


# ----------------------------------------------------------------------------------------------
# the reconstruction API on a finished run
# ----------------------------------------------------------------------------------------------


def build_reconstructor(world, *, fold=None, policy=None, ctx=None, ledger=None):
    """A Reconstructor over a FINISHED run. fold=None: the all-training states (for units outside the
    training partition). fold=k: the states that never saw outer fold k (for that fold's own units).

    Nothing is served unless it is what the finished run recorded: the declaration still describes the code
    and the inputs, every file bound by run_complete.json has its recorded SHA-256, and every state file holds
    exactly the state whose fingerprint and training groups were recorded when it was fitted.

    Where the EXPECTED fingerprint of a state comes from: for fold=k, selection.json, which run_complete.json
    binds by SHA-256; for fold=None, the record of the all-training states INSIDE run_complete.json. Never
    serving.json: that file is a copy, and verify_run refuses the run when the copy differs. The all-training
    state of a task must also be the engine setting that the bound selection chose for it."""
    out = Path(world.out)
    if not (out / "run_complete.json").exists():
        raise core.ReconError("no benchmark evidence exists yet for this world")
    frozen = verify_freeze(world)
    complete = verify_run(world)
    ctx = ctx or build_context(world)
    selection = core.read_json(out / "selection.json")
    serving = None
    if fold is None:
        if SERVING_KEY not in complete:
            raise core.ReconError("the completion record of this run does not bind its all-training states (it "
                                  "was written by an earlier version of this layer); they cannot be served")
        serving = complete[SERVING_KEY]                 # the bound record itself
    states, cache = {}, {}
    for name in ctx.tasks:
        if fold is None:
            chosen = serving["tasks"][name]
            if chosen["failed"]:
                raise core.ReconError(f"the all-training fit for {name} failed: {chosen['failed']}")
            selected = selection["serving"][name]
            if (chosen["representation"], chosen["components"]) != (selected["representation"], selected["components"]):
                raise core.FrozenSourceChanged(f"the all-training state recorded for {name} is not the engine setting "
                                               "that the selection of this run chose for it")
        else:
            chosen = selection["engine_used"][name][str(fold)]
        file = out / "states" / state_file("serving" if fold is None else f"outer{fold}", chosen["representation"],
                                           chosen["components"])
        if fold is None and not file.exists():
            raise core.ReconError(f"{file} is missing; the `serving` command restores the all-training states of a "
                                  "finished run from its completion record")
        if file not in cache:
            cache[file] = mods.JointState.load(file, expected_hash=chosen["state_hash"])
        if cache[file].groups_sha256 != chosen["training_groups_sha256"]:
            raise core.FrozenSourceChanged(f"{file}: the training groups are not the recorded ones")
        reader = {"state_format": mods.STATE_FORMAT, "support_margin": cache[file].support_margin,
                  "quantile_iterations": cache[file].quantile_iterations}
        if reader != frozen["state_reader"]:
            raise core.FrozenSourceChanged(f"{file}: the state was not made with the frozen reader settings")
        states[name] = cache[file]
    train = ctx.rows if fold is None else ctx.rows.take(np.flatnonzero(ctx.outer_of != fold))
    module = mods.JointMixtureModule(ctx.registry, ctx.layout, states=states, train=train, tasks=ctx.tasks,
                                     evidence=evidence_objects(world),
                                     default_state="T1" if "T1" in states else next(iter(states)),
                                     vote_alpha=ctx.config["knn"]["vote_alpha"],
                                     logistic_max_iter=ctx.config["logistic"]["max_iter"])
    reconstructor = core.Reconstructor(ctx.registry, policy, ledger)
    reconstructor.register(module)
    reconstructor.register(mods.AxonArborAdapter(ctx.registry))
    return reconstructor, ctx


def demo(world, log=say):
    """Requests through the API: one held-out-style answer per task and the declared refusals."""
    reconstructor, ctx = build_reconstructor(world, fold=0)
    index = int(np.flatnonzero(ctx.outer_of == 0)[0])
    unit, sample = str(ctx.table.neuron[index]), str(ctx.table.sample[index])
    soma = tuple(float(v) for v in ctx.table.soma[index])
    label = str(ctx.table.label[index])
    length = dict(zip(ctx.table.regions, (float(v) for v in ctx.table.length[index])))
    make = lambda observed, requested: core.ReconstructionRequest.from_values(
        ctx.registry, namespace="jiao", unit_id=unit, group_id=sample, observed=observed, requested=requested)
    requests = {
        "T1 profile from soma + label": make({"soma_position": soma, "source_label": label}, ["regional_axon_length"]),
        "T2 soma from profile + label": make({"regional_axon_length": length, "source_label": label}, ["soma_position"]),
        "T3 label from soma + profile": make({"soma_position": soma, "regional_axon_length": length}, ["source_label"]),
        "T4 outside block from inside block": make(
            {"regional_axon_length": {r: length[r] for r in ctx.world.t4_observed}, "soma_position": soma,
             "source_label": label}, ["regional_axon_length"]),
        "gene expression of this neuron": make({"soma_position": soma, "source_label": label}, ["gene_expression"]),
        "3D axon (adapter)": make({"soma_position": soma, "source_label": label}, ["axon_arbor"]),
        "3D axon from a measured profile": make({"soma_position": soma, "source_label": label,
                                                 "regional_axon_length": length}, ["axon_arbor"]),
        "unregistered variable": make({"soma_position": soma}, ["firing_rate"]),
        "unseen source label": make({"soma_position": soma, "source_label": "NotALabel"}, ["regional_axon_length"]),
    }
    rows = []

    def show(title, request):
        answers = reconstructor.reconstruct(request)
        for variable, result in answers.properties.items():
            rows.append({"request": title, "variable": variable, "status": result.status,
                         "codes": sorted(result.codes()), "module": result.module, "validation": result.validation,
                         "coverage_validated": result.coverage_validated, "baseline_level": result.baseline_level,
                         "reason": None if result.refusal is None else result.refusal.reason})
            log(f"{title}: {variable} -> {result.status} {sorted(result.codes())} [{result.module}]")
        return answers

    answers = {title: show(title, request) for title, request in requests.items()}
    # an inferred input brings its ancestry: a soma position inferred from the MEASURED profile is not an input
    # of a forecast of the hidden axon, whatever its name; one inferred from the source label alone is checked
    # and then only withheld because the adapter is blocked
    from_profile = answers["T2 soma from profile + label"]["soma_position"]
    if from_profile.status == "inferred":
        tainted = from_profile.to_observation(requests["T2 soma from profile + label"])
        show("3D axon from a soma inferred from the measured profile",
             make({"soma_position": tainted, "source_label": label}, ["axon_arbor"]))
        show("T1 profile from a soma inferred from the measured profile",
             make({"soma_position": tainted, "source_label": label}, ["regional_axon_length"]))
    core.write_json(Path(world.out) / "api_demo.json", {"unit": unit, "brain_sample": sample,
                                                        "state": "outer fold 0 states (the unit's sample is not in them)",
                                                        "requests": rows}, pretty=True)
    return rows


# ----------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------

METHOD_NAMES = {ENGINE: "Learned engine (joint mixture)",
                ENGINE_MEDIAN: "Learned engine, conditional median", "MEAN": "MEAN (training mean)",
                "IDENTITY": "IDENTITY (source label's training mean)", "COPY": "COPY (nearest training neuron)",
                "COPY_SAME_LABEL": "COPY, same source label", "KNN": "k-nearest-neighbour average",
                "RIDGE": "ridge regression"}
METRIC_NAMES = {"mse": "log MSE", "mae": "log MAE", "corr": "profile r", "euclid_mean": "mean error (um)",
                "mae_median_estimate": "log MAE of the engine's conditional median",
                "euclid_median": "median error (um)", "abs_AP": "AP (um)", "abs_DV": "DV (um)", "abs_ML": "ML (um)",
                "accuracy": "accuracy", "logloss": "log-loss"}
METHOD_ORDER = (ENGINE, ENGINE_MEDIAN) + mods.COMPARATORS
TASK_NAMES = {"T1": "T1: regional profile from soma + source label",
              "T2": "T2: soma position from the measured regional profile + source label",
              "T3": "T3: source label from soma + measured regional profile",
              "T4": "T4: regions outside the hypothalamus from the measured hypothalamic regions + soma + source label"}


def _fmt(value, digits=3):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "n/a"
    return f"{value:.{digits}f}" if abs(value) < 1000 else f"{value:,.0f}"


def render_report(world):
    """The results report, from the files of the run. For a FINISHED run they are read only if they are the
    files that its completion record binds: the report states that record's SHA-256 in its first lines, and no
    result may stand under that anchor that the anchor does not cover. (The table of API requests comes from
    api_demo.json, which the `demo` command writes after the run; the report says that it is an illustration.)"""
    out = Path(world.out)
    if (out / "run_complete.json").exists():
        verify_run(world)
    metrics = core.read_json(out / "metrics.json")
    selection = core.read_json(out / "selection.json")
    fits = core.read_json(out / "fits.json")
    declaration = core.read_json(out / "pre_run_declaration.json")
    lines = ["# SYNARCH reconstruction layer: masked-variable benchmark results", ""]
    # said first: the code that scored or reports differs from the code that was frozen
    now = code_hashes()
    changed = sorted(set(metrics.get("code_files_changed_since_freeze") or [])
                     | {k for k, v in declaration["code_sha256"].items() if now.get(k) != v})
    if changed:
        lines += ["**The code that scored or reported these results differs from the frozen code in: "
                  + ", ".join(changed) + ".** The predictions are the frozen run's; only scoring or reporting was "
                  "repeated.", ""]
    lines += [("**SYNTHETIC WORLD, not biological data.** " if world.synthetic else "")
              + "Generated by `analysis/synarch_recon_benchmark.py report` from the saved predictions.", "",
              "**Not confirmatory.** " + EXPOSURE, ""]
    if (out / "run_complete.json").exists():
        lines += ["Completion record of this run: `run_complete.json`, SHA-256 `"
                  + core.sha_file(out / "run_complete.json") + "`. It binds, by their SHA-256, the pre-run "
                  "declaration, the prediction manifest, the fit records, the selection, the metrics and the "
                  "evidence that the reconstruction API reads, and it holds the fingerprints of the all-training "
                  "states that answer for units outside the training partition.", ""]
    lines += [
              "Unit of independence: the brain sample (the ID prefix). The paper reports 166 mice and the data hold "
              "166 sample IDs; no statement was found that no mouse gave two samples. Every sample carries one "
              "source label, so a label effect is confounded with sample and batch. All intervals are descriptive.",
              "", "## Verdict per task", "",
              "| Task | Primary metric | Engine | Best comparator | Gain over it | Beats always-reported baselines "
              "(any amount) | Beats all comparators (any amount) | ... with interval | ... with interval and the "
              "minimum gain | Status for the API |", "|---|---|---:|---|---:|---|---|---|---|---|"]
    evidence = core.read_json(out / "evidence.json")
    minimum = float(declaration["config"].get("min_relative_gain", 0.0))
    for name, r in sorted(metrics["tasks"].items()):
        v, p = r["verdict"], r["primary_metric"]
        lines.append(f"| {name} | {METRIC_NAMES[p]} | {_fmt(r['methods'][ENGINE][p]['value'])} | "
                     f"{r['best_comparator']['name']} {_fmt(r['best_comparator']['value'])} | "
                     f"{100 * r['best_comparator']['relative_gain_of_engine']:+.2f}% | "
                     f"{v['always_reported_baselines']} | {v['all_comparators']} | {v['all_comparators_interval']} | "
                     f"{v['all_comparators_clear']} | {evidence['tasks'][name]['status']} |")
    lines += ["", "'every / some / none' counts the comparators that the engine beats. 'Any amount' compares means "
              "over brain samples and is satisfied by a difference of any size, also one far below what could "
              "matter. The interval rule also asks the descriptive 95% bootstrap interval of the paired difference "
              f"to exclude zero. The minimum gain is {100 * minimum:g}% of the comparator's own value (declared "
              "before the run). A task is `exploratory_supported` only if the engine beats EVERY comparator with "
              "the interval rule and the minimum gain; otherwise the API returns the baseline-level answer, "
              "flagged, or refuses. **The verdict of a task is that last comparison**: the engine beat every, some "
              "or none of the comparators clearly (interval and minimum gain). The two 'any amount' columns are "
              "shown because they are the literal comparison of means; a model that learned nothing satisfies "
              "them about half of the time.", "",
              f"Multiplicity: {len(metrics['tasks'])} tasks, one primary metric each. The status rule above is the "
              "only decision taken from this benchmark. Every other interval in this report (other metrics, "
              "single comparators, source labels) is an unadjusted description and decides nothing.", ""]
    for name, r in sorted(metrics["tasks"].items()):
        task_metrics = METRICS[r["kind"]]
        lines += [f"## {TASK_NAMES.get(name, name)}", "",
                  f"{r['neurons']} neurons from {r['samples']} brain samples, each predicted once by states that "
                  "never saw its brain sample.", "",
                  "| Method | " + " | ".join(METRIC_NAMES[m] for m in task_metrics) + " |",
                  "|---|" + "---:|" * len(task_metrics)]
        for method in (m for m in METHOD_ORDER if m in r["methods"]):
            summary = r["methods"][method]
            lines.append(f"| {METHOD_NAMES[method]} | " + " | ".join(_fmt(summary[m]["value"]) for m in task_metrics) + " |")
        lines += ["", f"Paired brain-sample differences on {METRIC_NAMES[r['primary_metric']]} (engine minus comparator; "
                  "descriptive 95% bootstrap interval over brain samples):", "",
                  "| Comparator | Difference | Interval | Samples favouring the engine | Engine better |",
                  "|---|---:|---|---:|---|"]
        for c in (m for m in mods.COMPARATORS if m in r["comparisons"][r["primary_metric"]]):
            cmp_ = r["comparisons"][r["primary_metric"]][c]
            b = cmp_["engine_minus_comparator"]
            lines.append(f"| {METHOD_NAMES[c]} | {_fmt(b['mean'], 4)} | [{_fmt(b['lower'], 4)}, {_fmt(b['upper'], 4)}] | "
                         f"{cmp_['samples_favouring_engine']} / {cmp_['samples']} | "
                         f"{'yes' if cmp_['engine_better_strict'] else 'NO'}"
                         f"{' (interval excludes 0)' if cmp_['engine_better_interval'] else ''} |")
        shown = list(METRICS[r["kind"]]) + [m for m in r["verdict"]["by_metric"] if m not in METRICS[r["kind"]]]
        by_metric = ", ".join(f"{METRIC_NAMES.get(m, m)}: {r['verdict']['by_metric'][m]}" for m in shown)
        lines += ["", f"Engine beats every / some / none of the comparators, by metric: {by_metric}.", ""]
        cov = r["coverage"].get(ENGINE)
        if cov:
            extra = (f", mean width {_fmt(cov['mean_width'])}" if "mean_width" in cov
                     else f", mean set size {_fmt(cov['mean_set_size'], 2)}")
            what = "prediction-set" if r["kind"] == "label" else "interval"
            lines.append(f"Engine {int(cov['level'] * 100)}% {what} coverage: {_fmt(cov['empirical'])}{extra}; "
                         f"coverage {'validated' if cov['validated'] else 'NOT validated'} (rule: {cov['rule']}; "
                         "tolerance 0.05)."
                         + (f" Share of the non-zero truths inside their interval: "
                            f"{_fmt(cov['share_of_nonzero_truths_inside'])} (a central interval is exactly [0, 0] "
                            "where zero is very likely, so this is a description, not a calibration target)."
                            if "share_of_nonzero_truths_inside" in cov else ""))
            for c in ("MEAN", "IDENTITY"):
                if c in r["coverage"]:
                    lines.append(f"{c} empirical interval: coverage {_fmt(r['coverage'][c]['empirical'])}, mean width "
                                 f"{_fmt(r['coverage'][c]['mean_width'])}.")
            lines.append("")
        if "unseen_label" in r:
            lines += [f"Unscoreable for every method (true label absent from that fold's training samples): "
                      f"{r['unseen_label']['neurons']} neurons in {r['unseen_label']['samples']} brain samples.", ""]
        elif r.get("rows_with_label_unseen_in_training"):
            lines += [f"{r['rows_with_label_unseen_in_training']} neurons had a source label with no training sample "
                      "in their fold; the engine marginalised the label (out of support) and the label-based "
                      "baselines fell back to their global rule.", ""]
        chosen = Counter(f"{u['representation']} K={u['components']}" for u in selection["engine_used"][name].values())
        one = sum(int(u["components"]) == 1 for u in selection["engine_used"][name].values())
        edge = sum(bool(selection["engine"][name][o]["at_grid_edge"]) for o in selection["engine"][name])
        tuned = "; ".join(f"{m}: {json.dumps(v, sort_keys=True)}"
                          for m, v in sorted(selection["serving"][name]["comparators"].items()))
        if r.get("rows_with_an_input_outside_the_training_range"):
            lines += [f"{r['rows_with_an_input_outside_the_training_range']} held-out neurons had an input outside "
                      "the range of their fold's training samples (the API calls that out of support); they are "
                      "scored like all others.", ""]
        lines += [f"Engine setting chosen on inner folds (count of outer folds): {dict(sorted(chosen.items()))}; "
                  f"component count at a grid edge in {edge} folds. K = 1 was chosen in {one} of "
                  f"{len(selection['engine_used'][name])} folds: with one component the engine has learned nothing "
                  "from its inputs and predicts the neuron-weighted training mean (label target: the neuron-weighted "
                  "label frequencies), up to the small prior smoothing of the model. Most frequent settings of the "
                  f"tuned comparators (MEAN and IDENTITY: who counts once): {tuned}.", "",
                  f"By source label ({METRIC_NAMES[r['primary_metric']]}; groups under five brain samples are "
                  "descriptive only):", "",
                  "| Source label | Samples | Neurons | Engine | " +
                  " | ".join(c for c in mods.COMPARATORS if c in r["methods"]) + " | Engine best | Note |",
                  "|---|---:|---:|---:|" + "---:|" * len([c for c in mods.COMPARATORS if c in r["methods"]]) + "---|---|"]
        for value, entry in r["by_source_label"].items():
            cells = [_fmt(entry["methods"][c][r["primary_metric"]]) for c in mods.COMPARATORS if c in r["methods"]]
            flags = entry["engine_better_than"]
            lines.append(f"| {value} | {entry['samples']} | {entry['neurons']} | "
                         f"{_fmt(entry['methods'][ENGINE][r['primary_metric']])} | " + " | ".join(cells) + " | "
                         f"{'yes' if flags and all(flags.values()) else ('n/a' if not flags else 'no')} | "
                         f"{'descriptive only' if entry['descriptive_only'] else ''} |")
        lines.append("")
    complete = core.read_json(out / "run_complete.json") if (out / "run_complete.json").exists() else {}
    serving_fits = complete.get(SERVING_KEY, {}).get("fits", [])
    every_fit = fits["engine_inner"] + fits["engine_outer"] + serving_fits
    not_converged = sum(1 for f in every_fit if "failed" not in f and any(not b["converged"] for b in f["diagnostics"]))
    failed = sum(1 for f in every_fit if "failed" in f)
    set_aside = aborted_starts(world)
    previous = declaration.get("supersedes")
    lines += ["## Fits", "",
              f"{len(fits['engine_inner'])} inner fits, {len(fits['engine_outer'])} outer fits, "
              f"{len(serving_fits)} all-training fits (fitted before scoring; never scored); {not_converged} stopped at the "
              f"iteration limit, {failed} failed. Deviations from the selected setting: {len(selection['deviations'])}. "
              f"Freezes or starts set aside before any outcome-derived file existed: {len(set_aside)}"
              + ("".join(f" ({a['aborted_utc']}: {a['reason']}; run had started: {a['run_had_started']})"
                         for a in set_aside)) + "."
              + (f" This run follows a started run that was set aside (declaration {previous['pre_run_declaration_sha256']}); "
                 f"the design is unchanged and these files differ: {', '.join(previous['changed_files']) or 'none'}."
                 if previous else ""), ""]
    demo_path = out / "api_demo.json"
    if demo_path.exists():
        demo_rows = core.read_json(demo_path)["requests"]
        lines += ["## The reconstruction API on one unit (states that never saw its brain sample)", "",
                  "An illustration, written after the run by the `demo` command (`api_demo.json`); it is not a "
                  "result of the benchmark and the completion record does not bind it.", "",
                  "| Request | Property | Outcome | Codes | Module |", "|---|---|---|---|---|"]
        for row in demo_rows:
            lines.append(f"| {row['request']} | {row['variable']} | {row['status']} | {', '.join(row['codes']) or '-'} | "
                         f"{row['module'] or '-'} |")
        lines.append("")
    lines += ["## Limits", "",
              "- Exploratory: the regional table was on disk and had been modelled before this benchmark was "
              "written. No result here is confirmatory, and no status is 'validated'.",
              "- Regional lengths are 49 published summaries on the scale log10(1 + length); they are not a "
              "partition of the brain, their unit is indicated but not resolved, and MY contains IO, NTS, PARN, IRN.",
              "- Intervals describe the fitted distribution, not parameter uncertainty.",
              "- A back-transformed log mean is not an expected length; only quantiles back-transform exactly.",
              "- Nothing here is a 3D path, a synapse, gene expression, or a causal effect of a label.",
              "- No validation or test brain sample was read; no trace was opened.", ""]
    return "\n".join(lines)


def report(world, log=say):
    if not (Path(world.out) / "run_complete.json").exists():
        raise core.ReconError("the benchmark has not completed; no results report is written")
    text = render_report(world)
    path = Path(world.out) / "results.md" if world.synthetic else ROOT / RESULTS
    path.write_text(text)
    log(f"report written: {path}")
    return path


# ----------------------------------------------------------------------------------------------
# synthetic worlds (same file formats as the real cohort; sealed rows are unparseable canaries)
# ----------------------------------------------------------------------------------------------

REAL_SHAPE = (("Adcyap1", 12, 709), ("Sst", 13, 813), ("Penk", 10, 304), ("Nts", 9, 542), ("Pdyn", 11, 445),
              ("Tac2", 10, 371), ("Oxt", 8, 290), ("Pomc", 6, 238), ("Tac1", 6, 386), ("Wt(orexin)", 4, 69),
              ("Agrp", 3, 87), ("Crh", 1, 16), ("Pmch", 3, 107), ("Avp", 2, 13), ("Trh", 1, 42), ("Vip", 2, 16))
TINY_SHAPE = (("Pomc", 10, 150), ("Sst", 10, 150), ("Nts", 8, 110), ("Wt(orexin)", 3, 30), ("Crh", 1, 12))
TINY_REGIONS = ("ARH", "VMH", "LHA", "ZI", "PAG", "MY", "BST", "PVT")
TINY_OBSERVED = ("ARH", "VMH", "LHA", "ZI")
SEALED = "SEALED_CANARY"


def _split_counts(total, parts, rng):
    """`parts` positive integers summing to `total`, uneven on purpose."""
    weights = rng.gamma(1.2, 1.0, parts) + 0.05
    counts = np.maximum(1, np.floor(weights / weights.sum() * (total - parts)).astype(int) + 1)
    while counts.sum() < total:
        counts[rng.integers(parts)] += 1
    while counts.sum() > total:
        candidates = np.flatnonzero(counts > 1)
        counts[rng.choice(candidates)] -= 1
    return counts.tolist()


def make_synthetic_world(directory, kind, *, size="tiny", seed=0, config=None):
    """Write a synthetic cohort in the real file formats. kind: 'planted' (soma, label and profile share
    latent classes) or 'null' (profile, soma and label are mutually independent).

    Validation and test samples are present in every file with UNPARSEABLE values, so a loader that
    parsed a sealed row would fail, and the forbidden columns hold poison strings.
    """
    import openpyxl

    if kind not in ("planted", "null"):
        raise ValueError("kind is 'planted' or 'null'")
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(core.stable_seed("synthetic-world", kind, size, seed))
    shape, regions, observed = (REAL_SHAPE, core.JIAO_REGIONS, core.JIAO_HYPOTHALAMIC) if size == "real" else \
        (TINY_SHAPE, TINY_REGIONS, TINY_OBSERVED)
    n_classes = 12 if size == "real" else 4
    labels = [s[0] for s in shape]
    alias = {"Wt(orexin)": "Orexin"}
    centres = np.column_stack([rng.uniform(5200, 7900, n_classes), rng.uniform(5000, 6900, n_classes),
                               rng.uniform(4100, 5500, n_classes)])
    nuclei = ["ARH", "VMH", "LHA", "PVH", "DMH", "ZI", "MPO", "TU", "PH", "AHN", "MM", "HY"][:n_classes]
    present = np.where(rng.random((n_classes, len(regions))) < 0.3, 0.9, 0.04)
    level = rng.uniform(2.0, 4.2, (n_classes, len(regions)))
    class_of_label = rng.dirichlet(np.full(n_classes, 0.25), len(labels))
    base_present = rng.uniform(0.05, 0.6, len(regions))
    base_level = rng.uniform(2.0, 4.0, len(regions))
    if size == "real":
        present[:, -1], base_present[-1] = 0.0, 0.0            # one region that is zero for everyone

    manifest, soma_rows, xlsx_rows, regional_rows = [], [], [], []
    sample_number = 300000
    counts = {"train_samples": 0, "train_neurons": 0}

    def add_sample(label, n_neurons, part):
        nonlocal sample_number
        sample_number += 1
        sample = str(sample_number)
        shift, gain = rng.normal(0, 40, 3), rng.normal(0, 0.08)
        for i in range(1, n_neurons + 1):
            neuron = f"{sample}_{i:03d}"
            if kind == "planted":
                c = int(rng.choice(n_classes, p=class_of_label[labels.index(label)]))
                soma = centres[c] + shift + rng.normal(0, 110, 3)
                on = rng.random(len(regions)) < present[c]
                y = np.where(on, np.maximum(level[c] + gain + rng.normal(0, 0.35, len(regions)), 0.05), 0.0)
            else:
                c = int(rng.integers(n_classes))
                soma = np.array([6500.0, 6000.0, 4800.0]) + rng.normal(0, [600, 400, 350])
                on = rng.random(len(regions)) < base_present
                y = np.where(on, np.maximum(base_level + rng.normal(0, 0.5, len(regions)), 0.05), 0.0)
            soma[2] = min(soma[2], 5695.0)
            fine = nuclei[c]
            full = f"{neuron}_{label}_{fine}"
            manifest.append([neuron, sample, alias.get(label, label), fine, "synthetic", "0", "0", "0", part, "synthetic"])
            if part == "train":
                soma_rows.append([full, "FORBIDDEN_TYPE", neuron, "FORBIDDEN_NAME", label, fine] + [f"{v:.6f}" for v in soma])
                xlsx_rows.append([full, "FORBIDDEN_CLUSTER", neuron, fine, label] + [float(v) for v in soma] + ["FORBIDDEN_CLASS"])
                regional_rows.append([full] + [f"{v:.4f}" for v in np.where(y > 0, 10 ** y - 1, 0.0)])
            else:
                soma_rows.append([full, "FORBIDDEN_TYPE", neuron, "FORBIDDEN_NAME", label, SEALED, SEALED, SEALED, SEALED])
                xlsx_rows.append([full, "FORBIDDEN_CLUSTER", neuron, SEALED, label, SEALED, SEALED, SEALED, "FORBIDDEN_CLASS"])
                regional_rows.append([full] + [SEALED] * len(regions))
        if part == "train":
            counts["train_samples"] += 1
            counts["train_neurons"] += n_neurons

    for label, n_samples, n_neurons in shape:
        for per_sample in _split_counts(n_neurons, n_samples, rng):
            add_sample(label, per_sample, "train")
        for part, n_held in (("validation", 1), ("test", 2 if size == "real" else 1)):
            for _ in range(n_held):
                add_sample(label, int(rng.integers(3, 12)), part)
    if size == "real":                                         # edge cases of the real table
        for row in regional_rows[:13]:
            if row[1] != SEALED:
                row[1:] = ["0.0000"] * len(regions)            # all-zero profiles
        first = next(r for r in soma_rows if r[5] != SEALED)
        second = [r for r in soma_rows if r[5] != SEALED][1]
        second[6:9] = first[6:9]                               # an exact duplicate soma
        for target in xlsx_rows:
            if target[2] == second[2]:
                target[5:8] = [float(v) for v in first[6:9]]
    order = rng.permutation(len(manifest))
    with open(directory / "manifest.tsv", "w", newline="") as handle:
        handle.write("\t".join(["neuron", "animal", "source_label", "soma_label", "class", "animal_inspected",
                                "animal_development", "animal_in_li_2024", "part", "part_reason"]) + "\n")
        for i in sorted(range(len(manifest)), key=lambda i: manifest[i][0]):
            handle.write("\t".join(manifest[i]) + "\n")
    with open(directory / "soma.csv", "w", newline="") as handle:
        handle.write("file,type,neurons,name,geno,location,AP,DV,ML\n")
        for i in order:
            handle.write(",".join(soma_rows[i]) + "\n")
    with open(directory / "regional.csv", "w", newline="") as handle:
        handle.write("," + ",".join(regions) + "\n")
        for i in rng.permutation(len(manifest)):
            handle.write(",".join(regional_rows[i]) + "\n")
    workbook = openpyxl.Workbook(write_only=True)
    sheet = workbook.create_sheet("Sheet1")
    sheet.append([None, "cluster", "neurons", "location", "geno", "AP", "DV", "ML", "Class"])
    for i in rng.permutation(len(manifest)):
        sheet.append(xlsx_rows[i])
    workbook.save(directory / "soma.xlsx")
    meta = {"name": f"synthetic_{kind}_{size}_{seed}", "synthetic": True, "kind": kind, "size": size, "seed": seed,
            "regions": list(regions), "t4_observed": list(observed), "expected": counts,
            "label_alias": [list(a) for a in alias.items()],
            "config": config or (json.loads(json.dumps(CONFIG)) if size == "real" else json.loads(json.dumps(TINY_CONFIG)))}
    core.write_json(directory / "world.json", meta, pretty=True)
    return core.synthetic_world(directory)


def run_synthetic(directory, kind, *, size="tiny", seed=0, workers=1, log=say, config=None):
    """Freeze, run once, score and report on one synthetic world. Returns (world, metrics)."""
    world = make_synthetic_world(directory, kind, size=size, seed=seed, config=config)
    freeze(world, "synthetic check: no review applies", log=lambda *a: None)
    metrics = run(world, workers=workers, directory=str(directory), log=log)
    return world, metrics


def timing(directory, log=say):
    """One full fit of the joint model on synthetic data of the real size, per component count."""
    world = make_synthetic_world(directory, "planted", size="real", seed=1)
    ctx = build_context(world)
    index = np.flatnonzero(ctx.outer_of != 0)
    rows = ctx.rows.take(index)
    out = {"world": world.name, "rows": int(len(rows)), "fields": len(ctx.layout.fields),
           "budget": ctx.config["engine"], "fits": []}
    for representation, components in engine_grid(ctx.config):
        start = time.time()
        state, diagnostics = mods.fit_joint(rows, ctx.layout, components=components, seed=1, budget=ctx.config["engine"],
                                            representation=representation)
        seconds = time.time() - start
        start = time.time()
        query = ctx.rows.take(np.flatnonzero(ctx.outer_of == 0))
        for task in ctx.tasks.values():
            mods.engine_predict(state, task, query, level=None if task.kind == "label" else 0.9,
                                median=task.kind != "label")
        model_fields = len(mods.model_field_types(ctx.layout, state.encoding))
        out["fits"].append({"representation": representation, "components": components, "model_fields": model_fields,
                            "fit_seconds": seconds, "predict_all_tasks_seconds": time.time() - start,
                            "iterations_best_start": diagnostics[0]["iterations"], "converged": diagnostics[0]["converged"]})
        log(f"{representation} K={components}: fit {seconds:.1f} s ({len(rows)} rows x {model_fields} model fields), "
            f"predict {out['fits'][-1]['predict_all_tasks_seconds']:.1f} s")
    core.write_json(Path(directory) / "timing.json", out, pretty=True)
    return out


# ----------------------------------------------------------------------------------------------
# how often does the status rule say 'supported' when there is nothing to learn? (synthetic only)
# ----------------------------------------------------------------------------------------------


NULLRATE_REAL_MAX_K = 4


def _nullrate_world(job):
    """One no-signal synthetic world through the whole frozen pipeline; the world is deleted afterwards."""
    directory, size, seed, config = job
    directory = Path(directory)
    if directory.exists():
        shutil.rmtree(directory)
    quiet = lambda *a: None                                        # noqa: E731
    start = time.time()
    world = make_synthetic_world(directory, "null", size=size, seed=seed, config=config)
    freeze(world, "no-signal simulation: no review applies", log=quiet)
    metrics = run(world, workers=1, directory=str(directory), log=quiet, post=False)
    evidence = core.read_json(Path(world.out) / "evidence.json")
    selection = core.read_json(Path(world.out) / "selection.json")
    tasks = {}
    for name, result in metrics["tasks"].items():
        v = result["verdict"]
        tasks[name] = {"status": evidence["tasks"][name]["status"],
                       "relative_gain_over_best_comparator": result["best_comparator"]["relative_gain_of_engine"],
                       "best_comparator": result["best_comparator"]["name"],
                       "all_comparators_interval": v["all_comparators_interval"],
                       "always_reported_baselines_any_amount": v["always_reported_baselines"],
                       "always_reported_baselines_clear": v["always_reported_baselines_clear"],
                       "engine_used": sorted({f"{u['representation']} K={u['components']}"
                                              for u in selection["engine_used"][name].values()}),
                       "mean_weighting": sorted({c["MEAN"]["selected"]["weighting"]
                                                 for c in selection["comparators"][name].values()})}
    shutil.rmtree(directory)
    return {"seed": seed, "seconds": time.time() - start, "tasks": tasks}


def nullrate(base, *, size="tiny", worlds=100, workers=1, seed0=1000, log=say):
    """How often does the status rule call a task 'exploratory_supported' when there is NOTHING to learn?

    Synthetic no-signal worlds only (soma, source label and regional profile mutually independent), each run
    through the whole frozen pipeline: freeze, nested folds, selection, scoring, status rule.
      size 'tiny': the small cohort (32 training samples) with the full small grid.
      size 'real': the real shape (101 training samples, 4,448 neurons, 49 regions, the real label counts), every
                   declared setting, except that the engine's component grid stops at K = 4 (K in {1, 2, 4}, both
                   representations). The three largest values (8, 12, 16) take most of the time of a run, about
                   ten minutes per world with them; where there is nothing to learn the inner folds choose the
                   smallest K almost always. The six comparators use their full declared grids.
    The summary is written to <base>/nullrate_<size>.json. Nothing here touches the real cohort."""
    base = Path(base)
    config = json.loads(json.dumps(CONFIG if size == "real" else TINY_CONFIG))
    if size == "real":
        config["engine"]["k_grid"] = [k for k in config["engine"]["k_grid"] if k <= NULLRATE_REAL_MAX_K]
    jobs = [(str(base / "nullrate_tmp" / f"{size}_{seed0 + i}"), size, seed0 + i, config) for i in range(worlds)]
    if workers <= 1 or len(jobs) <= 1:
        results = []
        for job in jobs:
            results.append(_nullrate_world(job))
            log(f"world {len(results)}/{len(jobs)}: " + ", ".join(
                f"{t} {r['status']}" for t, r in results[-1]["tasks"].items()))
    else:
        if __name__ != "__main__":
            raise core.ReconError("parallel workers are available from the command line only")
        import multiprocessing
        with multiprocessing.get_context("spawn").Pool(workers) as pool:
            results = pool.map(_nullrate_world, jobs, chunksize=1)
    names = list(results[0]["tasks"])
    supported = lambda r, t: r["tasks"][t]["status"] == "exploratory_supported"          # noqa: E731
    interval_only = lambda r, t: r["tasks"][t]["all_comparators_interval"] == "every"    # noqa: E731
    summary = {"worlds": len(results),
               "per_task": {t: {"supported_by_the_status_rule": sum(supported(r, t) for r in results),
                                "would_be_supported_without_the_minimum_gain": sum(interval_only(r, t) for r in results),
                                "beats_every_always_reported_baseline_by_any_amount": sum(
                                    r["tasks"][t]["always_reported_baselines_any_amount"] == "every" for r in results),
                                "beats_every_always_reported_baseline_clearly": sum(
                                    r["tasks"][t]["always_reported_baselines_clear"] == "every" for r in results),
                                "largest_relative_gain_over_the_best_comparator": max(
                                    r["tasks"][t]["relative_gain_over_best_comparator"] for r in results),
                                "relative_gains_where_the_interval_rule_alone_passed": sorted(
                                    r["tasks"][t]["relative_gain_over_best_comparator"] for r in results
                                    if interval_only(r, t))} for t in names},
               "worlds_with_at_least_one_task_supported": sum(any(supported(r, t) for t in names) for r in results),
               "worlds_with_at_least_one_task_supported_without_the_minimum_gain": sum(
                   any(interval_only(r, t) for t in names) for r in results)}
    record = {"schema": core.SCHEMA, "synthetic": True, "kind": "no-signal worlds: soma, label and profile independent",
              "size": size, "seed0": seed0, "created_utc": utc_now(), "code_sha256": code_hashes(),
              "config": config,
              "engine_component_grid": config["engine"]["k_grid"],
              "engine_component_grid_is_the_declared_one": config["engine"]["k_grid"] == (
                  CONFIG if size == "real" else TINY_CONFIG)["engine"]["k_grid"],
              "status_rule": "better than every comparator on the primary metric, every descriptive 95% interval "
                             "excluding zero, and a relative gain of at least min_relative_gain over each",
              "summary": summary, "per_world": results}
    if (base / "nullrate_tmp").exists() and not any((base / "nullrate_tmp").iterdir()):
        (base / "nullrate_tmp").rmdir()
    core.write_json(base / f"nullrate_{size}.json", record, pretty=True)
    log(f"{len(results)} no-signal worlds ({size}): at least one task 'supported' in "
        f"{summary['worlds_with_at_least_one_task_supported']} worlds "
        f"({summary['worlds_with_at_least_one_task_supported_without_the_minimum_gain']} without the minimum gain)")
    for t in names:
        log(f"  {t}: supported {summary['per_task'][t]['supported_by_the_status_rule']}, without the minimum gain "
            f"{summary['per_task'][t]['would_be_supported_without_the_minimum_gain']}")
    return record


# ----------------------------------------------------------------------------------------------
# command line
# ----------------------------------------------------------------------------------------------


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=["inputs", "splits", "synthetic", "timing", "freeze", "run", "score",
                                            "serving", "report", "status", "demo", "abort", "nullrate"])
    parser.add_argument("--world", default=None, help="directory of a synthetic world (default: the real cohort)")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--reviewed", default="", help="freeze: who reviewed the draft, and where")
    parser.add_argument("--reason", default="", help="abort: the technical reason")
    parser.add_argument("--size", default="tiny", choices=["tiny", "real"])
    parser.add_argument("--kind", default="both", choices=["planted", "null", "both"])
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", default=None, help="synthetic / timing / nullrate: directory for the synthetic worlds")
    parser.add_argument("--worlds", type=int, default=100, help="nullrate: number of no-signal synthetic worlds")
    args = parser.parse_args(argv)
    world = get_world(args.world) if args.command not in ("synthetic", "timing", "nullrate") else None
    if args.command == "inputs":
        input_audit(world)
    elif args.command == "splits":
        partition = core.TrainingPartition(world)
        table = core.load_training_table(world, with_outcomes=False, partition=partition)
        splits = make_splits(table, config_for(world))
        say(json.dumps({"splits_sha256": core.digest(splits),
                        "outer": {f: len(s) for f, s in splits["outer"].items()}}, indent=1))
    elif args.command in ("synthetic", "timing", "nullrate"):
        base = Path(args.out) if args.out else ROOT / "data/interim/synarch_recon/synthetic"
        if args.command == "timing":
            timing(base / "timing_world")
            return 0
        if args.command == "nullrate":
            nullrate(base, size=args.size, worlds=args.worlds, workers=args.workers, seed0=1000 + args.seed)
            return 0
        summary = {}
        for kind in (["planted", "null"] if args.kind == "both" else [args.kind]):
            directory = base / f"{kind}_{args.size}_{args.seed}"
            if directory.exists():
                raise core.ReconError(f"{directory} exists; a synthetic world is run once per directory")
            start = time.time()
            _, metrics = run_synthetic(directory, kind, size=args.size, seed=args.seed, workers=args.workers)
            summary[kind] = {"seconds": time.time() - start,
                             "verdicts": {t: r["verdict"] for t, r in metrics["tasks"].items()}}
        core.write_json(base / f"summary_{args.size}_{args.seed}.json", summary, pretty=True)
    elif args.command == "freeze":
        freeze(world, args.reviewed)
    elif args.command == "run":
        run(world, workers=args.workers, directory=args.world)
    elif args.command == "score":
        out = Path(world.out)
        if (out / "predictions_complete.json").exists() and not (out / "run_complete.json").exists():
            finish(world, how="score command, after an interrupted run")      # the first scoring of this run
        else:
            score(world)                                                      # a re-score: goes to rescores/<n>/
    elif args.command == "serving":
        fit_serving(world, workers=args.workers, directory=args.world)
    elif args.command == "report":
        report(world)
    elif args.command == "status":
        say(json.dumps(status(world), indent=1))
    elif args.command == "demo":
        demo(world)
    elif args.command == "abort":
        abort(world, args.reason)
    return 0


if __name__ == "__main__":
    sys.exit(main())
