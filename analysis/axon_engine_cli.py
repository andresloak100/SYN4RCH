"""Commands of the SYN4RCH axon prediction engine.

    python3 -I analysis/axon_engine_cli.py status
    python3 -I analysis/axon_engine_cli.py ingest [--stage develop|validate|test]
    python3 -I analysis/axon_engine_cli.py cache [--minutes M]        (all-pairs cache of the training library; resumes)
    python3 -I analysis/axon_engine_cli.py states [--which all|oof|all_training]
    python3 -I analysis/axon_engine_cli.py selfcheck
    python3 -I analysis/axon_engine_cli.py cost                       (measured cost; projection to 4,448 training neurons)
    python3 -I analysis/axon_engine_cli.py develop [--per-animal N] [--minutes M]
    python3 -I analysis/axon_engine_cli.py coverage                   (inputs only; training and validation animals)
    python3 -I analysis/axon_engine_cli.py proxy                      (training regional rows only: a development diagnostic)
    python3 -I analysis/axon_engine_cli.py predict --stage validate|test [--family all|no_fitted_state] [--minutes M]
    python3 -I analysis/axon_engine_cli.py copies-check --stage validate|test
    python3 -I analysis/axon_engine_cli.py admit --stage validate|test
    python3 -I analysis/axon_engine_cli.py seal                       (test only: the referee's declaration and seal)
    python3 -I analysis/axon_engine_cli.py score --stage validate|test
    python3 -I analysis/axon_engine_cli.py freeze
    python3 -I analysis/axon_engine_cli.py synthetic --kind informative|none|dense|informative_deterministic|all [--family ...]

`--world DIR` runs the same commands on a synthetic world directory instead of the real cohort.
`--world local37` runs them on the real cohort with the fixed set of 37 training traces under
data/raw/jiao2025/swc only (the development snapshot; outputs in data/interim/axon_engine/real_local37).

`predict` writes every forecast and the commitment and opens no held-out trace. `admit` opens the held-out
traces for the mechanical inclusion rules only and reads no forecast. `score` scores once: validation by this
engine, the test by the independent referee (`seal`, then its score_once). `score` runs `admit` (and, for the
test, `seal`) itself when they have not been run, so `predict` and `score` are enough; the single steps exist
so that the commitment and the seal can be checked in between. Nothing here downloads anything.
"""
import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import axon_engine as E  # noqa: E402


def show(value):
    print(json.dumps(E.jsonable(value), indent=1, sort_keys=True))


def world_of(args):
    if args.world == "real":
        return E.real_world()
    if args.world == "local37":
        return E.real_world(legacy_only=True)
    return E.synthetic_world(args.world)


def summary_table(output):
    lost = output.get("unavailable_at_source") or {}
    notice = ([f"THE SCORED COHORT IS SMALLER THAN THE LOCKED ONE: {lost['neurons']} of {lost['of']} neurons unavailable at the source "
               f"(rule 8); brain samples lost entirely: {lost['brain_samples_lost_entirely'] or 'none'}"]
              if lost.get("scored_cohort_smaller_than_locked_notice") else [])
    lines = notice + [f"{output['label']}", f"scored by: {output.get('scored_by', 'engine')}",
             f"scored {output['scored_queries']} of {output['declared_queries']} queries, {output['scored_animals']} animals (brain sample IDs)",
             f"{'method':34s} {'role':16s} {'distance um':>12s} {'F1':>8s}"]
    for m in E.METHODS:
        s = output["summary"].get(m)
        if s:
            lines.append(f"{m:34s} {E.method_role(m):16s} {s['distance_um']:12.1f} {s['f1_200']:8.4f}")
    other = output.get("other_candidates_descriptive", {}).get("summary", {})
    if other:
        lines.append("not selected (scored by the engine after the referee, descriptive, no claim):")
        lines += [f"{m:34s} {'candidate':16s} {s['distance_um']:12.1f} {s['f1_200']:8.4f}" for m, s in other.items()]
    return "\n".join(lines)


def selfcheck(world, workers, memo_queries):
    engine = E.Engine(world, "develop", workers=workers)
    engine.ensure_library()
    cache = engine.build_cache()
    missing = [s for s in engine.states.plan if s not in engine.states.manifest()["states"]]
    if missing:
        engine.states.build_many(missing, workers=workers)
    result = {"frozen_sources": E.FROZEN_PINS, "referee_sha256": E.referee_hashes(), "cohort": engine.cohort.summary(),
              "library": engine.library_summary(),
              "frozen_equivalence_on_type2_files": E.check_frozen_equivalence(engine),
              "two_copies": E.check_two_copies(engine),
              "pair_cache": {"pairs": cache["pairs"], "sha256": engine.cache.sha256, "bytes": engine.cache.record["bytes"],
                             "proof_against_frozen_pair_costs_and_referee_metric":
                                 engine.cache.prove(count=min(engine.cache.pairs, 2000))},
              "memo_proof": E.prove_memo(engine, queries=memo_queries),
              "official_copies_replica_against_referee_copying_baselines": E.check_official_copies(engine),
              "referee_preflight": E.referee_preflight()}
    if world.kind == "real":
        result["replay_of_earlier_medoid"] = E.replay_earlier_medoid(engine=engine)
    E.write_json(Path(world.out) / "checks" / "selfcheck.json", result, pretty=True)
    brief = dict(result)
    for key in ("frozen_equivalence_on_type2_files", "two_copies", "memo_proof"):
        brief[key] = {k: v for k, v in result[key].items() if k != "rows"}
    brief["library"] = {k: v for k, v in result["library"].items() if k != "exclusions"}
    return brief


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--world", default="real",
                        help="'real' (default: the locked cohort with every trace listed on disk), 'local37' (the same cohort with "
                             "only the 37 training traces of data/raw/jiao2025/swc: the development snapshot), or a synthetic world directory")
    parser.add_argument("--workers", type=int, default=4, help="worker processes (at most 4)")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("status")
    p = sub.add_parser("ingest")
    p.add_argument("--stage", choices=E.STAGES, default="develop")
    p = sub.add_parser("cache")
    p.add_argument("--minutes", type=float, default=None, help="stop cleanly after this many minutes; run again to resume")
    p = sub.add_parser("states")
    p.add_argument("--which", choices=("all", "oof", "all_training"), default="all")
    p = sub.add_parser("selfcheck")
    p.add_argument("--memo-queries", type=int, default=3)
    p = sub.add_parser("cost")
    p.add_argument("--sizes", type=int, nargs="+", default=[500, 1000, 2000, 4448])
    p = sub.add_parser("develop")
    p.add_argument("--per-animal", type=int, default=None,
                   help="forecast and score only the first N neurons of each animal (fixed hash order); default all")
    p.add_argument("--minutes", type=float, default=None, help="time budget of the forecasting step; run again to resume")
    sub.add_parser("coverage")
    sub.add_parser("proxy")
    p = sub.add_parser("predict")
    p.add_argument("--stage", choices=("validate", "test"), required=True)
    p.add_argument("--family", choices=tuple(E.FAMILIES), default="all",
                   help="candidate family admitted by the validation gate; declared when validation forecasts are made, pinned by the freeze")
    p.add_argument("--minutes", type=float, default=None, help="time budget; run the same command again to resume")
    for name in ("copies-check", "admit", "score"):
        p = sub.add_parser(name)
        p.add_argument("--stage", choices=("validate", "test"), required=True)
    sub.add_parser("seal")
    sub.add_parser("freeze")
    p = sub.add_parser("synthetic")
    p.add_argument("--kind", choices=tuple(E.SYNTHETIC_KINDS) + ("all",), default="all")
    p.add_argument("--size", choices=tuple(E.SYNTHETIC_SIZES), default="standard")
    p.add_argument("--family", choices=tuple(E.FAMILIES), default="all")
    p.add_argument("--dir", default=None, help="default: data/interim/axon_engine/synthetic[_<family>]")
    args = parser.parse_args(argv)
    workers = max(1, min(4, args.workers))
    started = time.time()
    try:
        if args.command == "synthetic":
            kinds = list(E.SYNTHETIC_KINDS) if args.kind == "all" else [args.kind]
            base = Path(args.dir) if args.dir else E.ROOT / "data/interim/axon_engine" / (
                "synthetic" if args.family == "all" else f"synthetic_{args.family}")
            for kind in kinds:
                out = E.run_synthetic(kind, base / kind, size=args.size, workers=workers, family=args.family)
                show({k: out[k] for k in ("world", "what", "family", "library", "pair_cache", "validate", "freeze",
                                          "validation_canary", "test_canary", "test", "once_only", "seconds")})
            return 0
        world = world_of(args)
        if args.command == "status":
            show(E.Engine(world, "develop", workers=1).status())
        elif args.command == "ingest":
            show(E.Engine(world, args.stage, workers=workers).build_library())
        elif args.command == "cache":
            engine = E.Engine(world, "develop", workers=workers)
            status = engine.build_cache(minutes=args.minutes)
            if status["complete"]:
                status["sha256"] = engine.cache.sha256
                status["proof"] = engine.cache.prove()
            show(status)
        elif args.command == "states":
            engine = E.Engine(world, "develop", workers=workers)
            names = [s for s in engine.states.plan if args.which == "all" or (args.which == "oof") == s.startswith("oof")]
            entries = engine.states.build_many(names, workers=workers)
            show({s: {k: e[k] for k in ("role", "animals", "neurons", "model_state_sha256", "animals_sha256",
                                        "validation_or_test_animals_in_state", "fit", "seconds")} for s, e in entries.items()})
            show({"regional_rows": engine.profile_audit})
        elif args.command == "selfcheck":
            show(selfcheck(world, workers, args.memo_queries))
        elif args.command == "cost":
            engine = E.Engine(world, "develop", workers=1)
            result = E.cost_probe(engine, sizes=tuple(args.sizes))
            E.write_json(Path(world.out) / "checks" / "cost.json", result, pretty=True)
            show(result)
        elif args.command == "develop":
            engine = E.Engine(world, "develop", workers=workers, develop_per_animal=args.per_animal)
            engine.ensure_library()
            progress = engine.pair_cache().progress()
            if not progress["complete"]:
                missing = progress["pairs"] - progress["pairs_known"]
                if missing > 200000:                 # hours of work: a step of its own, resumable
                    raise E.EngineError(f"the all-pairs cache of the training library lacks {missing} of {progress['pairs']} pairs; "
                                        "run `cache` (with `--minutes` it stops cleanly and resumes)")
                engine.build_cache()
            commit = engine.predict(minutes=args.minutes)
            if commit["status"] != "committed":
                show(commit)
                return 0
            output = engine.score()
            print(summary_table(output))
            sub16 = output.get("development_16_subset")
            if sub16:
                print(f"\nthe {sub16['neurons']} earlier development neurons only ({sub16['animals']} animals):")
                for m in E.METHODS:
                    print(f"{m:34s} {sub16['summary'][m]['distance_um']:12.1f} {sub16['summary'][m]['f1_200']:8.4f}")
            print(f"\nofficial copies: {json.dumps({k: v for k, v in output['official_copies'].items() if k != 'replica_difference_examples'})}")
            print(f"tables: {engine.stage_dir / 'tables.md'}")
        elif args.command == "coverage":
            result = E.coverage_diagnostic(E.Engine(world, "develop", workers=1))
            E.write_json(Path(world.out) / "checks" / "coverage.json", result, pretty=True)
            show(result)
        elif args.command == "proxy":
            result = E.regional_proxy(E.Engine(world, "develop", workers=1))
            E.write_json(Path(world.out) / "checks" / "regional_proxy.json", result, pretty=True)
            print(result["what"])
            show({k: v for k, v in result.items() if k not in ("methods", "what")})
            for name, row in result["methods"].items():
                print(f"{name:72s} MSE {row['animal_mean_mse_log10']:.4f}   profile r {row['animal_mean_profile_r']:.4f}")
        elif args.command == "predict":
            show(E.Engine(world, args.stage, workers=workers, family=args.family).predict(minutes=args.minutes))
        elif args.command == "copies-check":
            show(E.Engine(world, args.stage, workers=workers).copies_check())
        elif args.command == "admit":
            summary = E.Engine(world, args.stage, workers=workers).admit()
            show({k: v for k, v in summary.items() if k not in ("truth", "admitted_ids")})
        elif args.command == "seal":
            show(E.Engine(world, "test", workers=workers).seal())
        elif args.command == "score":
            engine = E.Engine(world, args.stage, workers=workers)
            if not (engine.stage_dir / "ADMISSION_COMPLETE.json").exists():
                engine.admit()
            if args.stage == "test" and not (engine.stage_dir / "referee" / "PREDICTION_SEAL.json").exists():
                engine.seal()
            output = engine.score()
            print(summary_table(output))
            show(output["verdict"])
            print(f"tables: {engine.stage_dir / 'tables.md'}")
        elif args.command == "freeze":
            record = E.Engine(world, "validate", workers=workers).freeze()
            show({"selected": record["selected"], "descriptive_reference": record["descriptive_reference"],
                  "family": record["family"], "validation": record["validation"], "referee": record["referee"],
                  "test_access": record["test_access"], "sha256": record["_sha256"]})
    except E.StageBlocked as blocked:
        print("BLOCKED: " + str(blocked))
        return 3
    except E.EngineError as refused:          # once-only, seal, leak, freeze mismatch, missing prerequisite
        print(f"REFUSED ({type(refused).__name__}): {refused}")
        return 4
    print(f"[{time.time() - started:.1f} s]", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
