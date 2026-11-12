# SYNTHETIC sensitivity benchmark — software results only

**This is a generic simulation, not hypothalamic validation.** It does not establish real projection accuracy, full axon paths, synapses or causal circuits. All gene identities, target identities, couplings and noise parameters are invented. The only empirical anchor is coordinate centre/scale from 37 previously exposed TRAIN neurons (32 animal proxies); that selected sample is not representative of hypothalamic populations.

Each of four fixed scenarios has six simulated training animals and four disjoint simulated test animals, 100 cells per animal, 16 synthetic genes and four synthetic targets. Seeds are fixed at 764101, 764102, 764103 and 764104, respectively. One seed per scenario supplies sensitivity examples, not precise Monte Carlo estimates.

The generator uses spatial mixtures, nonlinear radial and sinusoidal fields, latent classes, Gamma-Poisson gene counts and Bernoulli assay noise. It is independent of the logistic predictor. The location-only scenario excludes molecular state from the outcome law. In the hidden-state scenario, the exact conditional probability given measured features is 0.5 for each target.

## Fixed comparisons

All ten families/settings in the authorized frozen development recipe are retained. The updated gene-plus-position recipe uses C=0.1 and shrinkage alpha=0.1; the original C=1, alpha=0 candidate and two raw 1NN copies are additional descriptive anchors. Location LR remains C=0.1, alpha=0; location averaging uses k=31, alpha=0.5; joint and molecular averaging use k=31, alpha=0.1. The copy controls retain alpha=0.5. No settings or seeds were chosen from these synthetic test results. No A1 expression, labels, coordinates or coefficients entered the simulation.

A fixed missing-data adapter applies training medians after gene log1p and then training StandardScaler values, shared by all feature methods. This is an explicit extension to the frozen complete-input recipe. No test values enter imputation/scaling and no test labels, latent states or oracle probabilities enter fitting. Unknown target labels remain -1 and are excluded from scoring.

## Synthetic test scores

Lower is better. Entries are **Brier / log loss**, averaged equally over simulated animals within each target, then equally over targets. Log loss uses fixed clipping at 1e-6 and 1−1e-6. All 13 methods and individual animal/target scores are retained in the output files.

| Synthetic scenario | Updated genes + position | Original genes + position | Frozen location LR | Frozen location average | Frozen joint average | Constant 0.5 |
|---|---:|---:|---:|---:|---:|---:|
| Location only; genes add no information | 0.21749 / 0.62552 | 0.21940 / 0.63137 | 0.21190 / 0.61218 | 0.21431 / 0.61878 | 0.22179 / 0.63450 | 0.25000 / 0.69315 |
| Complementary molecular and spatial information | 0.19677 / 0.57904 | 0.19801 / 0.58297 | 0.24140 / 0.67580 | 0.24166 / 0.67633 | 0.19152 / 0.56928 | 0.25000 / 0.69315 |
| Unobserved, non-identifiable state | 0.25922 / 0.71245 | 0.26297 / 0.72111 | 0.25057 / 0.69436 | 0.25177 / 0.69679 | 0.25763 / 0.70924 | 0.25000 / 0.69315 |
| Animal batch shift, noise and missing features | 0.22519 / 0.64170 | 0.22985 / 0.65609 | 0.23434 / 0.66185 | 0.23240 / 0.65750 | 0.21597 / 0.62139 | 0.25000 / 0.69315 |

The updated recipe has lower log loss than the original candidate in these four fixed examples. It does not win every comparison. Location LR performs better when genes contain no incremental information. In the unobserved-state scenario, the updated recipe performs worse than constant 0.5. Joint 31-neighbor averaging performs better than the updated recipe in the complementarity and shifted/missing scenarios. These outcomes illustrate dependence on assumptions and available information; they are not estimates of performance on real hypothalamic neurons.

## Verification and reproducibility

- Twelve focused checks passed, covering fixed settings/seeds, reproducibility, null response laws, hidden-state non-identifiability, training-only imputation/scaling, test-label/state invariance, disjoint animals, donor eligibility and metric arithmetic. No test requires a model to win.
- Exactly one declared benchmark run completed. There were zero synthetic-test parameter searches, no scenario/seed deletion and no calibration after results.
- All 14 run-output hashes matched; 83,200 keyed prediction rows were unique; all 832 animal/target/method Brier and log-loss pairs were independently replayed from the saved predictions within 1e-12.
- The complete benchmark was not refit for verification. No real-data model fit, download, contact or account was used.

Files: `run_v1/summary.csv` (all method totals), `run_v1/metrics.json` (individual animals and targets), `run_v1/predictions.csv` (keyed predictions), `run_v1/SYNTHETIC_*.npz` (synthetic inputs/labels/latent audit states), `run_v1/fit_audit_*.json` (training-only transforms/models/donor groups), `run_v1/anchor_receipt.json`, `run_v1/STARTED.json` and `run_v1/RUN_COMPLETE.json`.

No biological success threshold or next model choice is inferred from this run.

## ROOT verification and project state

ROOT independently ran the 12 focused simulation checks and replayed the 83,200 saved predictions, checking exact synthetic identities, assay masks and labels against the four input archives. All 832 per-animal/target/method scores matched exactly; all 14 output hashes and source/protocol hashes matched. No predictor was refitted during ROOT verification.

ROOT verification receipt: `run_v1/ROOT_VERIFICATION.json`, SHA-256 `6bc2ae9ab0402335b32a6bc656795c715e4b5448cf7e80b4f4405177d6ca6cfc`. Full-run receipt SHA-256: `14efdad73b4872e7a5e09eaea755955c6890de84e0039d5570e80501153fad08`. Metrics SHA-256: `d271384c7e64c6c1faea82dc75de2577ddd82a84a1a6ace0369f6e246bbce07c`. Predictions SHA-256: `f66c9c9aa87260a2889a2c6c28260c31d1830a8e013c4d0b3f9ee0351808125a`.

The real A1 model and generic recipe are frozen in [models/paired_projection](../models/paired_projection/) after the [development improvement](paired_projection_logloss_development_results.md). New real hypothalamic accuracy is still untested. The user excludes lab collaboration; existing public data and clearly labeled simulations are the current route. These scripts and reports belong in andresloak100/SYN4RCH.
