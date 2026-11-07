# SYN4RCH generalisation probes (D-056)

Three questions asked of the result `docs/DECISIONS.md` D-055 established, each with a threshold fixed in D-056 **before** this module was written. The comparison is always `perturbation_axis_cell_majority` against `global_majority` on the slice `circuit_activity__primary_cell_seen`, the same pair and slice D-055 reports, so nothing here is a different measurement of a different thing.

This does **not** reopen D-048. That pre-registration named leave-one-study-out and its primary was met. A worse number under a harder split is a statement about the harder split.

## Run

| field | value |
|---|---|
| generated (UTC) | 2026-11-07T02:55:52Z |
| input | `data/curated/records.jsonl` |
| records read | 709 |
| labellable examples | 708 |
| excluded (counted) | 1 |
| cell | `(population_key, perturbation_axis, outcome_id)` |
| folds requested | 5 |
| bootstrap resamples | 1000 |
| permutations | 200 |
| seed | 20240101 |

Group counts, which are the subject of probe 1: `lab` 88 groups (largest 54 rows), `study_id` 102 groups (largest 34 rows).

## Probe 1 - is the signal cross-laboratory?

`lab` was `(last_author, affiliation)` until D-056. The affiliation is free text, the same department is spelled differently in every paper, and the corpus therefore held 75 studies and 75 "labs" -- so every leave-one-lab-out number this project published was a leave-one-study-out number under another name, on a corpus where four laboratories had published eleven papers between them. Both rows are printed because a key change that moves a result must be shown moving it (D-050).

| grouping | delta | clustered 95% CI | accuracy | rows / groups | McNemar p | verdict |
|---|---|---|---|---|---|---|
| `study_id` | +0.151 | [+0.029, +0.282] | 0.500 vs 0.349 | 192 / 45 | 0.016 | **established** |
| `lab` | +0.167 | [+0.037, +0.295] | 0.500 vs 0.333 | 192 / 40 | 0.007 | **established** |

## Probe 2 - replicated, or merely seen?

How many distinct `lab` groups had reached the test row's own cell in training. `depth_0` is the complement of the headline slice: there the candidate and the floor both fall back to the global mode, so a non-zero delta would mean the instrument is not what this project says it is.

| stratum | delta | clustered 95% CI | accuracy | rows / groups | McNemar p | verdict |
|---|---|---|---|---|---|---|
| `depth_0` | +0.000 | [+0.000, +0.000] | 0.317 vs 0.317 | 287 / 50 | n/a | **no_difference** |
| `depth_1` | +0.118 | [-0.078, +0.306] | 0.452 vs 0.333 | 93 / 25 | 0.177 | **no_difference** |
| `depth_2_plus` | +0.212 | [+0.039, +0.398] | 0.545 vs 0.333 | 99 / 28 | 0.024 | **established** |

Rows where the baseline's own backoff level and the training set disagree about whether the cell was populated: **0**. The two are computed independently and must agree; this line is the check.

## Probe 3 - is the `established` rule calibrated here?

labels permuted within study_id, so every study keeps its own label distribution and the cell structure is unchanged; what is destroyed is the association between a cell and a direction across studies. This is a check on the decision rule, not a p-value for the observed effect.

| field | value |
|---|---|
| permutations | 200 |
| observed delta | +0.151 |
| observed verdict | **established** |
| null deltas reaching the observed one | 0 of 200 |
| empirical p (add-one) | 0.005 |
| null delta range | -0.125 to +0.151 |
| null delta median / mean / 95th | +0.047 / +0.045 / +0.115 |
| null reached `established` | 6 (0.0300) |
| threshold | < 0.05 |
| calibrated | yes |

Verdicts over the null: `established` 6, `no_difference` 180, `suggestive` 14.

**The null is not centred on zero.** Its median delta is +0.047. Permuting within study preserves each study's own label multiset, so a cell that one study dominates keeps that study's labels and cell lookup still recovers study-level label structure. The observed +0.151 has to be read against that origin, not against 0, and the pre-registered threshold this probe tests is about how often the *verdict rule* fires -- not about where the null sits.

## What this report does not claim

* The permutation null is a check on the **decision rule**, under a null that keeps the study clustering. It is not a p-value for the observed effect, and a calibrated rule does not make an effect real.
* Depth is measured in groups that *reached a cell*, not in independent replications of an experiment. Two papers from one laboratory are one group here (D-056) and two papers sharing a mouse cohort are still two, because the schema records no cohort identity (D-013).
* Nothing here measures generalisation to an unseen population or to the molecular half. Both are already answered, negatively and exactly, in `reports/baseline_benchmark.json`: the candidate's delta against the floor is 0.000 under leave-one-population-out, and the molecular subset has no rows whose cell was populated in training.

