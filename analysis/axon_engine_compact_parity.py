"""Compact-adapter parity on real TRAINING development queries (three animals); writes the JSON given.

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    python3 -I analysis/axon_engine_compact_parity.py . data/interim/axon_engine/real/checks/compact_parity_train.json
"""
import json, sys, time
sys.path.insert(0, sys.argv[1] + "/analysis")
import axon_engine as E
engine = E.Engine(E.real_world(), "develop", workers=1, develop_per_animal=10)
engine.ensure_library()
qs = engine.queries()
animals, chosen = set(), []
for n in qs:                                   # training development queries only, three different animals spread over the list
    a = engine.cohort.rows[n].animal
    if a not in animals and len(chosen) < 3 and (len(chosen) == 0 or qs.index(n) > qs.index(chosen[-1]) + len(qs) // 4):
        animals.add(a); chosen.append(n)
t = time.time()
proof = E.prove_compact_parity(engine, chosen)
proof["seconds"] = time.time() - t
proof["stage"] = engine.stage
proof["engine_sha256"] = E.engine_code_sha().get("analysis/axon_engine.py")
E.write_json(sys.argv[2], proof, pretty=True)
print(json.dumps({k: v for k, v in proof.items() if k != "rows"}, indent=1))
for r in proof["rows"]:
    print(r["query"], r["pool_donors"], r["methods"], all(r["sixteen_records_identical"].values()), r["prediction_log_identical"],
          [(p["policy"], p["arbor_bitwise_equal"], p["records_equal_without_the_two"], p["executed_policy"]) for p in r["policies"]])
