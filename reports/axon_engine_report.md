# Axon prediction engine: what was built, what was checked, what is not known

2026-11-11, early morning. Report of the revision that followed the internal review, the other agent's read-only
review and the independent referee's contract. Code: `analysis/axon_engine.py`, `analysis/axon_engine_cli.py`.
Tests: `tests/test_axon_engine.py`. Rules: `reports/axon_engine_protocol.md` (STATUS: DRAFT). Outputs:
`data/interim/axon_engine/` (not version controlled). I committed nothing and downloaded nothing.

## In eight lines

1. **No independent biological test was run, and no biological accuracy is claimed.** The engine has opened no
   validation or test trace, and none was listed on disk when this was written (01:50).
2. The engine no longer uses a sampled donor panel. Every method draws on the **whole training library**, through a
   pinned all-pairs cache of the frozen metric. The frozen prototype runs unmodified on that whole pool.
3. The test now goes **through the independent referee**: its four functions are really called, end to end, on
   synthetic worlds. The four mandatory copying controls are the referee's own arrays.
4. Held-out traces are opened first by a new **admission** step that applies the mechanical inclusion rules and
   reads no forecast. Scoring a held-out set is once only, with no resume. Validation is bound to the method on
   disk **before** any truth is opened.
5. Development pass on the 37 local training traces (development, not a test): **no candidate beats the
   strongest copying control on both measures.** The numbers are those of the first build, because with 37 traces
   the panel already was the whole library.
6. Three things happened in the repository while I worked; none of them is mine. The referee's owner corrected the
   referee at 00:45 (the defect that blocked every method with a fitted regional state). Plan amendment 2 was
   committed (the user chose the access route; a rule for neurons the source does not serve; no test trace before
   the freeze): I implemented it. A download of training traces by the other agent started at 00:55 and is still
   running: I used none of it; every real-data number here comes from the fixed 37-trace snapshot.
7. Open: the protocol is a DRAFT, nothing is committed, the dense library is not complete, and the referee's own
   recomputation of the copies takes about 8 hours per call at full size and cannot be resumed or given a cache.
8. Smallest next action: when the other agent's listing of training traces is final, run `ingest`, `cache`,
   `states` and `develop` on it (training animals only). That is the first real answer about the dense library,
   and it spends neither validation nor test.

## 1. Why there is no biological result, and the exact blocker

At 01:37 on 2026-11-11 the listings held 787 of 4,448 training traces (the 37 early ones and 750 from the running
download), 0 of 1,185 validation traces and 0 of 1,547 test traces. The engine's answer on the whole cohort:

```
develop : BLOCKED: stage 'develop' on the whole cohort needs the final training listing: 787 of 4448 training
          neurons are listed (a file row or an unavailable_at_source row each); the download may still be running.
validate: BLOCKED: stage 'validate' cannot run on jiao2025: 787 of 4448 training traces are on disk (the dense donor
          library needs all); 0 of 1185 validation traces are on disk (scoring needs all). Blocker: the traces have
          not been fetched completely. Plan amendment 2 (2026-11-11): the user chose the portal's per-neuron address;
          stage A (training + validation) is downloaded by the other agent's analysis/axon_trace_download.py; stage B
          (test) only after the freeze. This engine downloads nothing.
test    : BLOCKED (the same, and 0 of 1547 test traces).
freeze  : REFUSED (EngineError): freeze needs a completed validation (PREDICTION_COMMITMENT.json is missing)
```

What stands between today and a test, in order:

1. The download of stage A (training and validation) must finish with a final listing. It is the other agent's.
2. The protocol must lose its DRAFT line (the user's acceptance) and engine, protocol and tests must be committed.
   The engine refuses to commit a held-out forecast while the protocol is a draft, and the real freeze needs a
   committed tree. I committed nothing.
3. Validation, then the freeze, then the test traces (stage B, fetched only after the freeze), then the test.

## 2. What happened in the repository during this work

I state these because they changed what I had to build, and because none came to me as a message.

| when | what | what I did |
|---|---|---|
| 00:45 | `analysis/prediction_independent_evaluator.py` and its tests changed (SHA-256 `eee1993e...` to `534d7e88...`; 24 tests to 51). The owner made the held-out set independent of the key order of `splits`. Until then the referee refused an honest training-animal list for a fitted regional state after it had frozen the declaration, and accepted a test-animal list (`reports/prediction_referee_interop_finding.md`, now marked resolved) | The engine noticed the change in mid-run and stopped ("the referee's files changed while this process was running"). I re-read the referee: only that check changed. I re-ran everything with the new file. The engine probes the referee with its own functions before any test truth is opened, so it would have blocked the affected bundles with the old file and passes them with the new one |
| 00:49, commit d43b654 | `reports/axon_dataset_plan_amendment_2.md`: the user said "download the missing neuron traces so we can run the real test against copying"; route = the portal's per-neuron address; the other agent owns the download; rule 8 (a neuron the source does not serve is excluded, logged, counted per part and per source label); no test trace on disk before the freeze. With it `reports/axon_engine_ingestion_fields.md` (what the listing must contain) | I treated the committed amendment as part of the binding plan, as the first build did with amendment 1. Rule 8 is in force in the engine with its reporting duty; the real freeze refuses while a test neuron is listed; `develop` on the whole cohort waits for the final training listing |
| 00:55, still running | `data/raw/jiao2025/traces/` appeared: a stage A download by `analysis/axon_trace_download.py` (2 workers; training files only so far, about 17 a minute: 787 of 4,448 training neurons listed at 01:37, 1,197 at 02:01) | I did not start, stop or touch it, and I did not ingest any of its files into a result. To keep my checks reproducible while the listing changes, the engine has a fixed development snapshot, `--world local37`: exactly the 37 traces under `data/raw/jiao2025/swc`. Before I noticed the download, one run of my unit tests read training traces through the whole-cohort listing while 297 were listed (three tests; training animals, which the plan allows in development); their outputs went to a temporary directory and are in no result. Those tests now use the snapshot |

The task text I was given says "no bulk fetch is approved". The amendment records a later decision by the user.
I did not act on either beyond making the engine able to read the listing when it is final.

## 3. What was built in this revision

| piece | what it does now |
|---|---|
| donor pool | every method uses the whole training library (for a training query in development: minus its own animal). The sampled panel is gone from code, protocol and tests |
| all-pairs cache | the frozen `geometry_metrics` on soma-centred arbors for every pair of library neurons (distance and 1 - F1), lower ID first, exactly as the frozen `pair_costs`; dense matrices with the ordered ID list, arbor hashes, metric code hashes and a SHA-256; built in parallel, resumable (`cache --minutes`), entries of unchanged arbors imported from older caches; sampled pairs are re-proved against the frozen `pair_costs` and against the referee's loader and metric every time forecasts are made |
| frozen prototype | its unmodified `predict_policy` receives the whole library; its `pair_costs` and `training_medoid` read the cache through a memo; trees and records are bit-identical with and without it |
| four official copies | `nearest`, `medoid`, `peptide_nearest`, `peptide_medoid` are the arrays returned by the referee's own `copying_baselines`. The engine's own computation (the contract's arithmetic on the cache) must agree; a difference is counted and the referee's forecast kept |
| other controls | the frozen baseline (`C1_medoid_copy`) and the two frozen uniform-attachment controls (C4, C5), on the whole library |
| candidates | F1, F2 (frozen conditional) and E1 to E6, unchanged in definition; `donor_ids` of E1 to E6 now list all 16 neighbours |
| candidate family | `all` (default) or `no_fitted_state` (E1, E2, E3, E6: no regional state, none of its exposure); declared when validation forecasts are made, pinned by the freeze |
| `predict` | all 15 methods for every declared query, manifest, commitment; opens no held-out trace; resumable until the commitment (`--minutes`) |
| `copies-check` | for a pool too large for the inline call: the referee's `copying_baselines` on the whole stage, exact equality required before any truth is opened |
| `admit` | new. Checks that the committed forecasts were made by the method on disk now (code, configuration, protocol, cohort tables, training traces and library, cache, regional state bytes), writes the opened marker, applies inclusion rules 1 to 6 and rule 8 to the held-out traces, writes the exclusion log and the canonical truth. Reads no forecast |
| validation scoring | by this engine, once, no resume; every score is also computed by the referee's metric module and must agree bit for bit; frozen gate, then engine gate over the declared family |
| `freeze` | pins everything a forecast depends on; fits the all-training regional state again and requires the stored hash; refuses on the real cohort while a test neuron is listed; records the referee's hashes and which pinned files the commit does not hold |
| test | `seal`: the referee's `freeze_declaration`, `copying_baselines`, the committed forecast files copied byte for byte, `seal_predictions`. `score`: the referee's `score_once`. The verdict and the tables are read from its result |
| reports | per animal, by source label, by soma label; neurons and animals in every row; every subgroup row is descriptive, sparse below 5 animals; the verdict sentence carries the paired differences against the strongest control on each measure with their intervals |

Commands (all `python3 -I analysis/axon_engine_cli.py [--world real|local37|DIR] ...`): `status`, `ingest`,
`cache`, `states`, `selfcheck`, `cost`, `develop`, `coverage`, `proxy`, `predict --stage`, `copies-check --stage`,
`admit --stage`, `seal`, `score --stage`, `freeze`, `synthetic`.

**Why validation is not scored by the referee.** The referee takes one frozen candidate and scores the split it
calls test. Validation has to score all eight candidates to choose one. So validation is scored by the engine
with the frozen metric (bit-identical to the referee's on every score), with the referee's own four copies as
controls; the test, where the claim is made, goes through the referee. The protocol says so in section 7.

**When the gate selects nothing.** The verdict is "no method beyond copying" and stays that. The gate also
records a descriptive reference (the family member with the lowest validation distance) only to fill the
referee's single candidate slot if the test is run to report the controls. Nothing can be claimed for it.

## 4. Checks run and their outcomes

Final code: engine `e575ddc3...`, commands `cd7dc111...`, tests `ff80b942...` (full hashes in section 12);
referee `534d7e88...`.

| check | outcome |
|---|---|
| unit tests, `python3 -I -m pytest tests/test_axon_engine.py -q` | **101 passed** (148 s on a quiet machine; 266 s in the last run, while other processes were busy) |
| the referee's own tests, unchanged | 51 passed |
| frozen files pinned; a changed byte detected | passed; the 9 pinned files and the freeze record match |
| ingestion of the 37 snapshot traces | 37 included, 0 excluded; 16 code the axon as type 2, 21 as type 0; table soma and trace soma differ by at most 1.9 um (35 of 37 below 0.01 um); 5 axons have cable outside the atlas box (464 to 1,932 um, 4,914 um in all), kept; no connector node; the 4 neurons with two copies agree |
| my extraction against the frozen extractor, 17 type-2 file copies | bit-identical canonical arbors on the frozen extractor's node set; the plan's rule adds terminal marker nodes, at most 20.5 um per neuron |
| replay of the earlier figure on the 16 development neurons | 922.4851 um / 0.19782, exactly the earlier value, by the frozen extractor and by my extraction (history only, never a threshold) |
| all-pairs cache of the 37 arbors (666 pairs) | every pair bit-identical to the frozen `pair_costs` and to the referee's metric; SHA-256 `5f11b380...`, the same after a rebuild |
| cache build stops at its time budget, resumes, gives the same pinned cache; a changed matrix file is refused | passed (synthetic) |
| unmodified frozen policies on the whole pool, with and without the cache (3 queries x 5 policies, 36 donors each) | bit-identical distances, F1 losses, medoid, trees and records |
| the four copies: my computation against the referee's `copying_baselines` | 148 of 148 forecasts equal (donor and arrays) on the 37 traces, all 32 leave-one-animal-out pools; equal on every synthetic pool checked |
| floating-point tie between the two medoid arithmetics (the review's counterexample) | reproduced: the frozen sum picks donor 2, the referee's matrix product donor 3; my memo follows the frozen function, my copy follows the referee |
| the referee's four functions really called, in the contract's order | passed on synthetic worlds (a test wraps them to count the calls); the engine's animal means equal the referee's bit for bit |
| a bundle with a fitted regional state | honest training-animal lists at all three levels; sealed and scored by the referee (with the file of 00:45); a test animal in any list is refused |
| admission reads no forecast; exclusions are in the referee declaration before the seal; the referee's single attempt is not spent on a missing truth | passed |
| once only | second `predict`, `admit`, `seal`, `score`, `freeze`: refused. A validation scoring whose completion marker was removed and whose truth was replaced: refused (the reviewer's probe). An interrupted admission: completed only on the same commitment and trace files |
| validation bound to the method before truth | code, configuration, state bytes or training library changed after the forecasts: `admit` refuses and writes no marker (the reviewer's probe) |
| stale library | an ingestion log holding part of the listed traces is never used: `predict` ingests again (the reviewer's probe) |
| forged regional state with every stored hash rewritten | accepted by the loader, refused by the refit at the freeze |
| hidden-outcome invariance | held-out traces and regional rows replaced, everything rebuilt: all forecasts of the synthetic worlds bit-identical; forecasts also unchanged with the held-out traces deleted |
| planted leaks | a held-out animal in a state, a query's animal in its state or among its donors: refused each time |
| rule 8, test traces before the freeze, a listing still being written, a second copy that disagrees | passed (synthetic), section 2 |
| real `predict`, `admit`, `score` for validation and test; real `freeze` | BLOCKED (exit 3) and REFUSED (exit 4), section 1 |

## 5. Development pass on the 37 snapshot traces

**Development data, not a test. No success claim.** Leave one animal out: donors are the traces of the other 31
animals (35 or 36 neurons: the whole library); the regional state is the out-of-fold state that excludes the
query's animal. 37 neurons, 32 animals (brain sample IDs). Forecasts anchored at the table soma. The four official
copies are the referee's arrays for every animal (my own computation agreed on all 148).

| method | role | animal-mean distance, um | animal-mean F1 | the 16 development neurons only (14 animals) |
|---|---|---:|---:|---|
| `nearest` | official copy (C2) | 945.0 | **0.1866** | 932.8 / 0.1806 |
| `medoid` | official copy (C1 by the referee's rule) | **824.6** | 0.1773 | 835.5 / 0.1872 |
| `peptide_nearest` | official copy (C3) | 1,053.1 | 0.1389 | 1,147.1 / 0.1385 |
| `peptide_medoid` | official copy | 1,154.2 | 0.1193 | 1,255.2 / 0.1232 |
| `C1_medoid_copy` | frozen baseline (C1) | 824.6 | 0.1773 | 835.5 / 0.1872 |
| `C4_uniform_attach_500` | frozen control (C4) | 824.6 | 0.1773 | 835.5 / 0.1872 |
| `C5_uniform_attach_1000` | frozen control (C5) | 825.7 | 0.1756 | 838.1 / 0.1833 |
| F1 frozen conditional 500 | candidate | 824.6 | 0.1773 | 835.5 / 0.1872 |
| F2 frozen conditional 1,000 | candidate | 825.7 | 0.1756 | 838.1 / 0.1833 |
| E1 local medoid, uniform | candidate | 818.3 | 0.1856 | 845.7 / 0.1976 |
| E2 local medoid, kernel | candidate | 800.8 | 0.1805 | 813.3 / 0.1964 |
| E3 kernel x label | candidate | 824.0 | 0.1806 | 838.8 / 0.1992 |
| E4 kernel x regional | candidate | 812.5 | 0.1798 | 852.7 / 0.1906 |
| E5 kernel x label x regional | candidate | 828.5 | 0.1778 | 840.5 / 0.1913 |
| E6 consensus-pruned | candidate | 921.9 | 0.1846 | 962.8 / 0.2097 |

Paired animal differences against the strongest control on each measure (distance: `medoid`; F1: `nearest`),
with 95% bootstrap intervals over the 32 animals:

| candidate | distance difference, um | interval | F1 difference | interval |
|---|---:|---|---:|---|
| E1 | -6.2 | -47.6 to +39.4 | -0.0010 | -0.0349 to +0.0340 |
| E2 | -23.7 | -80.6 to +30.8 | -0.0061 | -0.0386 to +0.0281 |
| E3 | -0.5 | -63.9 to +62.8 | -0.0060 | -0.0415 to +0.0297 |
| E4 | -12.1 | -69.8 to +43.7 | -0.0068 | -0.0426 to +0.0304 |
| E5 | +3.9 | -60.3 to +66.9 | -0.0088 | -0.0441 to +0.0260 |
| E6 | +97.3 | +14.1 to +185.1 | -0.0020 | -0.0369 to +0.0337 |
| frozen conditional 500 (candidate F1) | +0.0 | +0.0 to +0.0 | -0.0093 | -0.0479 to +0.0292 |
| frozen conditional 1,000 (candidate F2) | +1.2 | +0.0 to +3.5 | -0.0110 | -0.0485 to +0.0267 |

Read plainly:

* No candidate beats the strongest control on both measures. E1 to E4 have a lower mean distance than every
  control, and none reaches the F1 of `nearest`. By the rule this is "no method beyond copying".
* Every difference between a candidate and the strongest control lies inside its interval around zero, except
  that E6 is worse on distance. E2 is closer than `medoid` in 10 of 32 animals. "E2 lowers the distance by 24 um"
  is within noise.
* The frozen prototype is its baseline here. Its attachment is refused when the template's soma is farther from
  the query than a third of the radius: 37 of 37 forecasts fell back at 500 um and 36 of 37 at 1,000 um. The
  conditional weights picked another template for 10 of 37 queries; none of the ten could be attached.
* The frozen baseline and the referee's `medoid` chose the same neuron for all 37 queries.
* This table says little about the dense library. Here the nearest donor is a median 327 um away (134 to 880) and
  the 16 "nearest" donors reach out to 1,080 um, almost half the library.
* Earlier figures for comparison only: representative copy 922.5 um / 0.198 on the 16 development neurons with a
  15-donor library. The last column forecasts the same 16 neurons from a 36-donor library, so it is another experiment.

## 6. Synthetic worlds (mechanism checks; not biological evidence)

Each world: 24 training, 8 validation, 12 test animals, 8 to 14 neurons each, traces written as SWC files in the
real layout with both axon codings, both frames, markers, dendrites, cable outside the atlas box and planted
violations of the inclusion rules. The whole pipeline ran on each: library, all-pairs cache (about 35,000 pairs),
six states, development, validation (forecasts, canary, admission, scoring, both gates), freeze, test forecasts,
canary, admission, the referee's seal and `score_once`, refusal of every repeat. About 2.5 minutes per run. Two of
the four worlds were also run with the family without a fitted state. In the development pass of a synthetic
world the referee's `copying_baselines` is called for the first 2 of the 24 leave-one-animal-out pools (the
declared pair budget); in validation and test for the whole stage.

| world | family | validation gate | test (once, scored by the referee) |
|---|---|---|---|
| label and soma zone set the dominant route (probability 0.75) | all | only E5 eligible, selected | E5 587.6 um / 0.527; strongest controls 620.6 um (`peptide_medoid`) and 0.425 (`peptide_nearest`); matched controls beaten: **success**. Margins: -33.0 um (interval -83.6 to +23.8) and +0.102 F1 (+0.051 to +0.155) |
| the same world | no_fitted_state | nothing eligible (E3 has 693.7 um against 563.3 um for `peptide_medoid`): "no method beyond copying" | run with E3 in the candidate slot as descriptive reference: 859.9 um / 0.390, not better than the copies; nothing claimed |
| no information (route at random) | all | nothing eligible | reference: the frozen conditional 500 candidate, 1,228.3 um / 0.124, equal to its uniform twin C4; nothing claimed |
| soma zone sets the dominant route, label is noise | all | E1, E2, E4, E5, E6 eligible; E4 selected (lowest distance); E3 (label) not eligible | E4 732.9 um / 0.491 against 1,039.3 um / 0.362 for the strongest control (`peptide_nearest`), but E1 has 731.9 um: **rule not met**, nothing claimed for the regional weights |
| the same world | no_fitted_state | E1, E2, E6 eligible; E2 selected | E2 742.6 um / 0.476 against 1,039.3 um / 0.362: **success**. Margins: -296.7 um (-422.4 to -158.3) and +0.113 F1 (+0.060 to +0.162) |
| label and soma zone set the route with certainty | all | nothing eligible: `peptide_nearest` 209.0 um / 0.695 is already nearly optimal | reference E5 210.3 um / 0.698 happens to be better than every control on the test animals (`peptide_nearest` 220.1 um / 0.659): **nothing is claimed**, because validation selected nothing |

What this shows: the pipeline separates the cases; the matched-control requirement bites (fourth row); a world in
which copying is already optimal is reported as such (last row), and a reference that happens to win on the test
after losing on validation earns nothing; with nothing selected the test still runs through the referee and
claims nothing. The declared family changes what may be selected, not what is forecast or reported. In all six
runs every hidden-outcome canary was bit-identical (20,205 forecasts compared), every repeat of a once-only step
was refused, the cache proof held, the referee sealed bundles with a fitted-state method (rows 1, 3, 4 and 6: the
training animals behind the state declared at all three levels), and my computation of the four copies equalled
the referee's on every validation and test query. The margins are smaller than in the first build, as the
internal review predicted, because the representative copies are now medoids of the whole library: in the first
row the distance margin over the strongest control is 33 um, with an interval that includes zero. Passing here
proves the machinery, nothing about neurons.

## 7. Cost at full size (measured here, projected to 4,448 training neurons)

Measured on the 37 snapshot arbors (median 1,070 canonical nodes, mean 1,330), one process, nothing else running.

| step | measured | projection for the whole library |
|---|---|---|
| frozen metric, one pair | 3.4 to 4.0 ms | all-pairs cache, 9,890,128 pairs: 9.5 to 11 processor hours; about 2.5 hours with 4 workers; resumable in blocks of 2,000 pairs |
| cache on disk | 24 kB for 37 | two float64 matrices of 158 MB and one byte matrix of 20 MB: 336 MB, read by memory map |
| referee's metric, one pair | 2.7 to 3.1 ms | **one call of the referee's `copying_baselines`: 7.5 to 8.6 hours in one process**, not resumable; it is called once per `copies-check` (validation, test) and once inside the seal: three calls |
| unmodified frozen `predict_policy`, pool of 4,448 donors of real size, memoised | 7.0 s for the five policies of one query (0.6 s for the baseline, 1.5 to 1.7 s for each attached policy); 2.8 s once for the medoid; peak memory 2.6 GB | validation, 1,185 queries: 2.3 processor hours; test, 1,547 queries: 3.0; `predict` is resumable, so no single run has to exceed a quarter of an hour |
| engine candidates, scoring | milliseconds per query | minutes |

The pool for the frozen-policy timing is the local arbors repeated under made-up IDs with placeholder pair costs:
it times the frozen function's own work per call (its checks of every donor, and the two 4,448 x 4,448 matrices
it writes into its record as lists). No forecast was made from it.

**No method needs a restricted donor pool.** What remains untested is the real thing: no cache of the real
library exists because the library is not on disk, so every projection above assumes the other 4,411 arbors are
the size of the 37.

**The referee's part is its own.** Its contract says full-cohort runtime "has not been demonstrated" and that a
"training-only reusable distance cache or an independently reviewed exact implementation is needed". The engine's
cache is exactly that matrix, proven bit-identical to the referee's metric on every pair checked. Request to the
referee's owner: a reviewed way to hand it a pinned distance matrix. I did not patch the referee.

## 8. Where the scientific bottleneck is, in my judgement

**a. The answer depends on a library that is still arriving.** Everything measured in 3D so far comes from 37
traces. From the soma tables of training and validation animals only: with all training traces a query's nearest
donor from another animal is a median 60 um away (training, leave one animal out; 59 um for validation), and the
16 nearest lie within a median 173 um. Whether neighbours 60 um apart have similar axons at a 200 um tolerance is
not known. `develop` on the complete training library answers that on training animals.

**b. The one dense signal on disk is regional, and it favours position, not the label.** On the 4,448 training
rows of the regional table, leave one animal out, with profile distance in place of geometry (a proxy, not a test):

| "forecast" of the 49-region log profile | squared error | profile correlation |
|---|---:|---:|
| global medoid copy | 2.12 | 0.31 |
| nearest-soma copy | 2.40 | 0.40 |
| same-label nearest copy | 2.39 | 0.39 |
| same-label medoid copy | 2.21 | 0.31 |
| local medoid of the 16 nearest (uniform / kernel) | 1.96 / 1.98 | 0.46 / 0.46 |
| local medoid with the label factor | 2.03 | 0.45 |
| local medoid, K = 8 / K = 32 (uniform) | 2.07 / 1.90 | 0.45 / 0.47 |
| mean of the 16 nearest (not a copy) | 1.31 | 0.58 |
| frozen joint model, out of fold (not a copy) | 1.44 | 0.49 |

At regional resolution the local medoid beats all four kinds of copy on both columns, which is the pattern the
engine needs in 3D. The label adds nothing on top of position: the 16 nearest somas hold a median 3 same-label
donors, and 26% of training queries (24% of validation queries) have none.

**c. The ceiling of this family is copying itself.** A mean of neighbours beats the best single copy by a wide
margin at regional level (1.31 against 1.96). Every candidate outputs one measured arbor, and arbors cannot be
averaged. E6 is a first, crude step away from that; on the sparse snapshot it trades recall for precision (0.14
/ 0.35 against 0.18 / 0.22 for E2) and loses on distance.

**d. The regional state adds little that the neighbours do not carry, for the frozen model as it is.** In the
proxy the frozen joint model (1.44) is a worse regional predictor than the plain mean of the 16 nearest neurons
(1.31), and on the snapshot E4 is not better than E2. This holds for the frozen model with at most 8 components
only: all six fits select 8, the largest value of the frozen grid, and the inner-fold error is still falling
there (1.547, 1.535, 1.499, 1.444 for 1, 2, 4, 8 components in the all-training fit). A larger grid would be a
new declared candidate, fitted on training animals only; I did not touch the frozen grid.

**e. The frozen prototype can show little on this cohort.** It picks one template from the whole pool regardless
of where the query is, then refuses to attach it when the soma shift exceeds a third of the radius. From soma
tables: the best single training neuron lies within 167 um of 1.6% of validation queries and within 333 um of
7.3%. So the uniform controls C4 and C5 must equal the baseline for at least 98% and 93% of validation queries.
The two frozen conditional candidates change template with the query, but through eight regional components,
not through distance.

**f. The F1 measure is low for everything.** No method reaches 0.19 on the snapshot: at 200 um even the best copy overlaps a
fifth of the cable. Whether the dense library lifts that is the open question of (a).

## 9. The review findings, one by one

Internal review (15 problems):

| | finding | status |
|---|---|---|
| 1 | controls and frozen prototype on a sampled panel; no all-pairs cache | **fixed.** Section 3; no panel remains |
| 2 | the engine never calls the referee; its bundle is refused | **fixed** for the test (four functions, section 4). Validation stays with the engine's scorer, for the reason in section 3; this is a deliberate limit, not an omission |
| 3 | referee defect: held-out set depends on key order | the owner corrected it at 00:45. The engine never declared an empty list for a fitted state and reordered nothing; it probes the referee before any test truth is opened |
| 4 | inclusion rules reached held-out truth only at scoring; the referee drops nothing | **fixed**: `admit` |
| 5 | a finished scoring could be repeated | **fixed**: a started scoring is never repeated or resumed. One narrower rule kept on purpose: an admission that died (no forecast read, nothing scored) may be completed on the same commitment and the same trace files |
| 6 | validation bound to the method only at the freeze | **fixed**: all pins are compared before the opened marker |
| 7 | `predict` on a stale library | **fixed**: a library that is not the one of the listed traces is never used |
| 8 | full-pool cost not measured; nothing resumable | **fixed**: section 7; `cache` and `predict` resume |
| 9 | the coverage diagnostic used test animals' somas and labels | **fixed and disclosed**: it now uses training and validation animals only; the old record is kept as `real_local37/checks/coverage_with_test_inputs_2026-11-10.json`; the exposure statements say so; no setting was changed |
| 10 | statements stronger or older than the evidence | this report replaces them; the referee's soma tolerance is 5 um (as the plan's rule 2), not 2 |
| 11 | one unavailable trace blocks a stage; no rule | the user's amendment 2 created rule 8; **implemented** with its reporting duty |
| 12 | constant freeze flag; commit check too narrow; forged state accepted | **fixed**: section 3 (`freeze`) |
| 13 | subgroup "population statement"; verdict without uncertainty | **fixed**: every subgroup row is descriptive; the verdict carries the paired intervals and the sentence that a numerical win is not a confirmatory population claim |
| 14 | no list of what is frozen with its tests; two tests missing; dendritic-stem rule open; `donor_ids` | **fixed**: protocol section 4 has the list; the stem rule is fixed now. One part **not taken as written**: for the frozen attached policies (C4, C5 and the two frozen conditional candidates) `donor_ids` holds the selected template only. Their template is chosen among the whole pool, 4,448 IDs per row would be noise, and the referee requires exactly one donor for a copy; the whole pool is recorded once by hash |
| 15 | "the regional state adds little" drawn from a model at its grid ceiling | **fixed**: section 8d |

The other agent's review (`reports/axon_engine_codex_review.md`): clipping removed, the two counterexamples are
regression tests on its own coordinates (closed); full copying controls connected to the referee, on the whole
pool, with the tie arithmetic settled by taking the referee's arrays (closed); state provenance and the binding
of validation to the frozen method (closed). Its request for "one consolidated review" is for whoever reviews
this revision; I did not treat its sentences as orders.

## 10. What is open, and what only the user can decide

* **The protocol is a DRAFT.** Accepting it includes two choices: the candidate family (`all`, or `no_fitted_state`
  if the exposure of the regional model class is to be kept out of the claim), and the one permitted reset of K on
  training animals.
* **Nothing is committed.** Engine, protocol and tests must be committed before the real freeze. The frozen
  prototype and the referee are untracked files too; the freeze lists them with their hashes.
* **The dense-library behaviour of every candidate is untested on real traces.** K, the kernel, the label factor
  and the pruning settings were fixed by reasoning, not by data.
* **The referee at full size**: about 8 hours per call, three calls, one process, no resume (section 7).
* **Listing details I cannot see yet**: how the downloader will write a neuron the source does not serve. The
  engine accepts a row with `frame` = `unavailable_at_source` with or without a path. The download's completion
  marker is the downloader's; the engine only requires that every neuron of a needed part is listed.
* **Exposure.** Validation and test animals are new geometry animals, not a pristine resource: the 49-region
  table of all 7,180 neurons has been on disk since 2026-11-10; the other agent's regional model was fitted and
  scored by 5-fold whole-animal cross-validation over all 166 animals; summaries were printed and one raw profile
  was displayed (the last two taken from the earlier audit, not re-verified by me). Soma positions and labels of
  test animals, which are inputs and not outcomes, were summarised in the plan and in the first coverage record.
  The engine loads training rows only; F1, F2, E4 and E5 inherit the regional model class.
* **Unit.** "Animal" means brain sample ID. The paper reports 166 mice and there are 166 IDs; I found no statement
  that no mouse gave two samples.
* Permitted inputs not used by any method: soma label, reconstruction type, hemisphere beyond the frame rule.

## 11. When the listing is final

```
python3 -I analysis/axon_engine_cli.py status
python3 -I analysis/axon_engine_cli.py ingest --stage validate     # all listed training traces
python3 -I analysis/axon_engine_cli.py cache --minutes 14          # repeat until complete (about 2.5 h in all)
python3 -I analysis/axon_engine_cli.py states
python3 -I analysis/axon_engine_cli.py develop --per-animal 10 --minutes 14   # training only: the first dense-library answer
#   accept the protocol (remove 'STATUS: DRAFT'); commit engine, protocol, tests
python3 -I analysis/axon_engine_cli.py predict --stage validate --minutes 14  # repeat until committed
python3 -I analysis/axon_engine_cli.py copies-check --stage validate          # the referee's copies, about 8 h, one process
python3 -I analysis/axon_engine_cli.py admit --stage validate
python3 -I analysis/axon_engine_cli.py score --stage validate      # once; both gates
python3 -I analysis/axon_engine_cli.py freeze                      # refuses while a test neuron is listed on disk
python3 -I analysis/axon_engine_cli.py predict --stage test --minutes 14      # needs no test trace
python3 -I analysis/axon_engine_cli.py copies-check --stage test
#   stage B: the test traces are fetched now, not before
python3 -I analysis/axon_engine_cli.py admit --stage test
python3 -I analysis/axon_engine_cli.py seal                        # the referee's declaration and seal (it recreates the copies)
python3 -I analysis/axon_engine_cli.py score --stage test          # the referee's score_once, once
```

The `develop` line needs only the training traces and spends neither validation nor test.

## 12. Hashes of this revision

| file | SHA-256 |
|---|---|
| `analysis/axon_engine.py` | `e575ddc32d5d5e71b30a88744a012e509ee68ce6953f6d5f37596cd676df7d23` |
| `analysis/axon_engine_cli.py` | `cd7dc1117fb52c176a5908c4559cfda139ffc7efe895b65322ade1e3e353c40f` |
| `tests/test_axon_engine.py` | `ff80b94241ac2d40a5816a4f886cb1900b4dc0b1e33bf2913424f3ffbe8a23c1` |
| `reports/axon_engine_protocol.md` | `71fb18ef0af2edf7a00f0fa2b95d951d5763f4bf6b3868e3f9d04b583c43dad0` |
| `analysis/prediction_independent_evaluator.py` (the referee, not mine, as used) | `534d7e88833542b597b7f3abc529983f48eb22a7d46ab2681be836aca341dc0d` |
| `analysis/brm_geometry_audit.py` (the referee's metric module, not mine) | `9c3abfdb56fc8e1de53358957b9140b3a19be994a8d87e9961b5b878bbdfc998` |

The synthetic bundles hold as a method artifact the protocol text as it was when they were sealed (SHA-256
`5fe3923402c6e943e47abed42f9d1ae9e6984b96a5467038b50a5230d6feee12`); two sentences were made more precise afterwards (connector
geometry; rule 8 in the admission step). No rule changed.

Records: `data/interim/axon_engine/real_local37/` (`ingestion_log.json`, `state_manifest.json`,
`paircache/45d643c6cb50bc69/`, `develop/` with 555 forecasts, `prediction_manifest.json`, `metrics.json` and
`tables.md`, `checks/selfcheck.json`, `checks/cost.json`, `checks/coverage.json`, `checks/regional_proxy.json`,
`checks/coverage_with_test_inputs_2026-11-10.json`); `data/interim/axon_engine/synthetic/<world>/summary.json`
and `data/interim/axon_engine/synthetic_no_fitted_state/<world>/summary.json`, each with its `engine/` directory
(the referee's bundle is `engine/test/referee/`).
