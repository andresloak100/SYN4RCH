# D-108 Phase 1: automatic discovery of transformations (Epi-retro-seq, non-hypothalamic sources)

## Declared before the main runs (2026-11-10 12:27:36)

Timing pilot (one full nested search on shuffled labels): 129 s of CPU per nested search (7 outer folds); 58 s wall with 5 workers.

Replicates fixed now: shuffled-label null N0 60 searches; D-as-truth null N2 60 searches (beta 0.259); planted controls 8 per truth x 3 truths (betas GxL_tc 0.805, chain_distLhat 0.502, nmf_tg 0.586); symbolic layer on 60 N0 replicates. 4 worker processes.

Decision rules:

* **search_beats_D:** YES iff delta = mean7(search) - mean7(D) > 0 AND delta above the 95th percentile of delta under N2 (p_N2 = (1 + #N2 delta >= observed)/(1 + n_N2) < 0.05) AND exact one-sided division sign-flip p < 0.05; NO iff delta <= 0; NOT SHOWN otherwise
* **search_beats_zero:** p_N0 = (1 + #N0 search mean7 >= observed)/(1 + n_N0)
* **discovered_route:** a term selected in >= 4 of 7 outer folds AND that count above the 95th percentile of the largest per-term count under N0 AND mean outer ablation delta > 0
* **planted_recovery:** planted term selected in >= 4 of 7 outer folds; if every planted truth is recovered in < 50% of replicates the search is called underpowered
* **symbolic:** a formula beats simpler alternatives iff mean7(formula) > mean7(best leaf, chosen the same nested way) with exact division sign-flip p < 0.05 AND mean7(formula) above the 95th percentile of the SR outer mean7 under N0

Disclosure: During the timing check (before this declaration) one real-label outer fold was run by mistake: isocortex, search r -0.049 (Lhat x target centroid + E-hat x target centroid) vs D +0.310. Nothing in the search was changed after it.

---

## Results (main runs 2026-11-10 13:15 to 14:09; post-audit additions 14:21 to 14:40; nothing in the search changed after the declaration)

**In short.** The automatic search did not beat model D. On held-out divisions it scored r = -0.062 against D's +0.228. That is worse than every comparator, and no better than the same search run on shuffled labels. The search overfits its training divisions. On shuffled labels the pipelines it selects score +0.24 inside the training divisions, against +0.02 for D's term. That optimism from selection alone is about as large as D's entire real signal. Planted controls show what the engine can and cannot find. It reliably finds a simple one-coefficient chain (distance from the inferred location) when the effect is large. It misses a planted methylation (E-hat) route. When D is true and a realistic extra route is added (worth about +0.04 in r), the search gains only +0.01, which cannot be detected. The symbolic-regression layer found no formula that beats a single variable. **search beats D: NO.**

### What was run

* Data: Epi-retro-seq Supplementary Table 2, 27 non-hypothalamic sources, 187 (source, subclass) units. Rows with Source == "HY" are dropped when the sheet is read (`d105_benchmark.read_st2_training`). The 19 HY test subclasses are removed from every training source (asserted in `d108_prims.world`). Unlabelled atlas data (Allen 10x clusters, MERFISH cells, hypothalamic atlas cells included) are used only in unsupervised, target-side and placement steps, the same way for D and for the search. Methylation enters only as one optional node (E-hat).
* Metric: r_unit exactly as in d105 (per held-out division, sources pooled, sqrt(n / mean n) unit weights), averaged over the 7 divisions (mean7). Outer folds leave one division out, and held-out subclasses are removed from the training sources (d105 folds).
* Primitive library: a **fixed typed library of 66 pre-built terms**, 16 unit-side nodes x 4 target encodings plus 2 pairwise unit x target features.
  * Gene bases: D's whole-brain PCA 32 (gD32) and its first 8 (gD8); K's top-2000-gene PCs (gK32); B's 10 PCs (gB10); NMF rank 16 (nmf16); a masked-autoencoder latent, 4000-256-16 (ae16; loss 0.771 to 0.313 over 40 epochs). All were fitted on all 5,322 atlas clusters, as D's PCA was.
  * Nonlinear expansions: random Fourier features of gD8 (rff_gD8) and a degree-2 polynomial of gD4 (poly_gD4).
  * Location: measured MERFISH centroid (Lm, rff_Lm); location inferred from genes (Lhat, rff_Lhat, module M1); inferred division composition (Lcomp, module M2); gene x location products (gD4xLm, gD4xLhat).
  * Methylation (optional): E-hat, a mapping from genes to L1 methylation composition (RBF kernel ridge), cross-fitted over subclasses inside each fold on the training sources only.
  * Target encodings: one-hot (D's head), target centroid (tc), target gene profile (tg) and target cell-type composition (tcomp). Pairwise: -log distance from the measured or inferred location to the target centroid (distLm, distLhat).
  * **What the library does not do (critic change 1 only partly met).** Chains (genes to M1 to inferred location to distance to target, or to RFF) and the cross-modality products are hand-enumerated nodes. The search selects pre-built terms and adds them together. It does not build new chains or products. There is no gene x methylation or location x methylation term, and no product or multiple-kernel blocks. D's exact recipe (gD32 x one-hot, lambda 0.3) is one term. Refitted alone, it reproduces D's d105 rehearsal score (+0.228).
* Search: greedy forward selection of up to 3 terms, fitted jointly by one weighted ridge with lambda in {0.1, 0.3, 1, 3, 10}. Each candidate is scored by mean inner leave-one-division-out r_unit, inside the outer training divisions only. A term is added only if the inner score rises by >= 0.01. The space holds 47,971 pipelines of 1 to 3 terms (x 5 lambdas). Greedy search scores at most 195 term sets (975 term-set x lambda candidates) per outer fold, each over 6 inner folds.
* Comparators on identical folds:
  * Fitted here: frozen D; D_nested (D's recipe with lambda tuned in the inner loop); K, B and A (frozen d105 recipes); O2 (D with the target-kernel penalty, frozen d105 recipe, added after the audit); R_nested (ridge on all of K's gene PCs, lambda tuned in the inner loop).
  * Cited, not rerun: T1 (fine-tuned autoencoder, 6 seeds), per division from `reports/d105_rehearsal.md` on the same folds.
* Tests: paired one-sided sign-flip (exact over the 7 divisions; Monte Carlo with 200,000 draws over the 27 sources and over the 187 units).
  * N0: the whole nested search on 60 unit-shuffled datasets.
  * N2: the whole nested search on 60 datasets simulated with D as the truth (beta 0.259, tau 0.35, calibrated so that D's LODO r matches the real +0.228). The N2 truth is D's in-sample fit on all 27 sources (`d_truth`), not a cross-fitted profile. This is acceptable for a generator, and is disclosed here.
  * Planted controls: 8 replicates for each declared truth, plus the post-audit truths described below.

### Main result (real labels, nested; r_unit per held-out division)

| division | search | D | D_nested | O2 | K | B | A | T1 (d105 rehearsal) |
|---|---|---|---|---|---|---|---|---|
| amygdala-pallidum | +0.260 | +0.215 | +0.215 | +0.236 | +0.162 | +0.335 | +0.352 | +0.227 |
| hindbrain | -0.259 | -0.009 | -0.009 | +0.034 | -0.030 | -0.128 | +0.213 | +0.071 |
| hippocampal | -0.124 | +0.376 | +0.376 | +0.388 | +0.298 | +0.223 | +0.081 | +0.343 |
| isocortex | -0.049 | +0.310 | +0.335 | +0.345 | +0.249 | +0.360 | +0.003 | +0.479 |
| midbrain | -0.014 | +0.111 | +0.111 | +0.116 | +0.162 | +0.183 | +0.059 | +0.296 |
| olfactory | +0.101 | +0.268 | +0.268 | +0.289 | +0.340 | +0.225 | -0.198 | +0.293 |
| thalamus | -0.349 | +0.326 | +0.219 | +0.321 | +0.055 | -0.073 | -0.058 | +0.251 |
| **mean7** | **-0.062** | **+0.228** | +0.216 | +0.247 | +0.176 | +0.161 | +0.065 | +0.280 |

R_nested (regularised regression on all gene PCs, lambda tuned in the inner loop) chose lambda 1.0 in every fold and gave +0.176, identical to K in every fold. It has the same features as K, and its grid contains K's lambda, so it is **not an independent regularised-regression baseline**. T1's rehearsal score is known to be optimistic from selection: on the hypothalamic test in D-105 it scored +0.08.

Search minus comparator, with a one-sided sign-flip p over the 7 divisions and the number of divisions where the search is better:

| vs | mean delta | p | search better |
|---|---|---|---|
| D | -0.290 | 0.992 | 1 of 7 |
| D_nested | -0.278 | 0.992 | 1 of 7 |
| O2 | -0.309 | 0.992 | 1 of 7 |
| T1 (rehearsal) | -0.342 | 0.992 | 1 of 7 |
| K | -0.238 | 0.992 | 1 of 7 |
| B | -0.223 | 1.000 | 0 of 7 |
| A | -0.127 | 0.898 | 1 of 7 |
| R_nested | -0.238 | 0.992 | 1 of 7 |

Against D at finer levels: per source, mean delta -0.340 (27 sources, p = 1.0); per unit (cosine), mean delta -0.185 (187 units, p = 1.0).

The isocortex fold reproduces the disclosed accidental pilot fold exactly (-0.049 vs +0.310), as expected from a deterministic search.

### Nulls and the optimism of selection

* **N0, shuffled labels (60 full searches).** The search's mean7 averaged -0.009 (95th percentile +0.118, maximum +0.259). The real search scored -0.062, so p_N0 = 0.754. On real labels the search does no better than on shuffled labels.
* **N2, D as the truth (60 full searches).** D's mean7 averaged +0.219 and the search's +0.147. Delta (search - D) averaged -0.072 (sd 0.077, 95th percentile +0.053; delta > 0 in 18% of replicates). This is the cost of searching when the simple model is right. The search picked D's own term in only 6.2% of folds (first in 3.8%). The observed delta of -0.290 is in the far lower tail of the N2 distribution: 1 of 60 N2 replicates is lower (minimum -0.318). p_N2 = 0.984.
* **Measured selection optimism.** This is the mean inner (training-division) score of the selected pipeline, against the inner score of D's term, in the same folds:

| data | selected pipeline, inner | D's term, inner | selected pipeline, outer | folds selecting D's term |
|---|---|---|---|---|
| N0 shuffled (420 folds) | +0.240 | +0.021 | -0.009 | 0.5% |
| N2 D is the truth (420 folds) | +0.365 | +0.198 | +0.147 | 6.2% |
| real (7 folds) | +0.410 | +0.123 to +0.261 | -0.062 | 0 of 7 |

  On shuffled labels the selected pipelines score +0.24 inside the training divisions with no signal at all. This is a winner's curse over about 200 term sets x 5 lambdas, scored on 6 inner folds. It is about as large as D's whole real inner signal (+0.123 to +0.261). The real search's inner score (+0.41) exceeds D's by about what N2 shows when D is the truth (+0.365 vs +0.198), so the inner gap is not evidence that a different form fits the real data better.

### Verdict (declared rule)

delta = -0.290 <= 0, so **search_beats_D = NO**. It also fails the other two conditions: p_N2 = 0.984 and division sign-flip p = 0.992.

### Discovered pipelines per fold and stability

| held-out division | selected terms (unit node x target encoding) | lambda | inner score | outer r | rank of D's term at step 1 (of 66) |
|---|---|---|---|---|---|
| amygdala-pallidum | ae16 x tcomp + gD8 x tc | 10 | 0.355 | +0.260 | 12 |
| hindbrain | ae16 x tcomp + poly_gD4 x onehot | 10 | 0.387 | -0.259 | 10 |
| hippocampal | ae16 x tc + Lhat x tc + Lm x tcomp | 10 | 0.438 | -0.124 | 38 |
| isocortex | Lhat x tc + Ehat x tc | 1 | 0.359 | -0.049 | 21 |
| midbrain | rff_gD8 x tcomp + Ehat x tcomp + Lhat x tcomp | 0.1 | 0.504 | -0.014 | 32 |
| olfactory | rff_gD8 x tc + poly_gD4 x tc + rff_gD8 x tg | 10 | 0.388 | +0.101 | 10 |
| thalamus | nmf16 x tcomp + rff_Lhat x tcomp + gK32 x tcomp | 3 | 0.440 | -0.349 | 30 |

* Inner scores (0.35 to 0.50) are far above the outer scores (mean -0.06). The winners are low-dimensional target encodings (tc, tcomp), and the N0 optimism above shows that winners of this size arise without any signal.
* Selection frequency: 16 distinct terms across the 7 folds. The most frequent terms (ae16 x tcomp, Lhat x tc) were each selected in 2 of 7 folds. Under N0 the 95th percentile of the largest per-term count is 4.05. **No route counts as discovered** (rule: selected in >= 4 of 7 folds, count above 4.05, and positive mean ablation).
* Ablation is the outer score lost when a term is dropped, refitted at the same lambda:
  * ae16 x tcomp: mean +0.021 (2 folds).
  * Lhat x tc: mean -0.322 (2 folds: -0.068 in hippocampal and -0.577 in isocortex, so in isocortex dropping it would have raised the outer score by 0.58).
  * Every other term appears in a single fold, with values between -0.183 and +0.084.
  * Post-audit change: a term selected alone now has ablation = its own outer r (against predicting nothing). Before this, a lone term could never count as a discovered route. No real fold selected a single term, so nothing changes here.

Ablation = r(full pipeline) - r(pipeline without the term). A positive value means the term helped on the held-out division; a negative value means it hurt.

### Planted positive controls

Every planted truth sits on top of the real counts. "Oracle" means the true term(s) alone, with lambda tuned in the same inner loop and scored on the held-out division (for D + route truths, the oracle is D's term plus the route). The extended recovery rule (post-audit, declared at 14:21:31 in `d108_declared_supp.json` before these runs) requires the planted term in >= 4 of 7 folds AND search mean7 >= oracle mean7 - 0.10. The 3 declared truths were rerun with the same seeds to add the oracle score, and their selected terms and outer r reproduced exactly.

| planted truth | recovered: declared rule | recovered: extended rule | per-replicate fold counts | search mean7 | D mean7 | oracle mean7 | search - D (p over reps) | reps search > D | reps with division p < 0.05 |
|---|---|---|---|---|---|---|---|---|---|
| gene x measured location, target centroid (gD4xLm x tc), r 0.40 | 6 of 8 | 6 of 8 | 4,4,3,2,4,5,5,5 | +0.200 | +0.191 | +0.116 | +0.009 (0.336) | 4 of 8 | 0 of 8 |
| chain: genes to inferred location to distance (distLhat), r 0.40 | 8 of 8 | 6 of 8 | 7 in every replicate | +0.320 | +0.206 | +0.387 | +0.114 (0.004) | 8 of 8 | 2 of 8 |
| NMF programmes x target gene profile (nmf16 x tg), r 0.40 | 2 of 8 | 2 of 8 | 2,0,2,4,1,3,2,4 | +0.281 | +0.350 | +0.318 | -0.069 (0.992) | 1 of 8 | 0 of 8 |
| *post-audit:* E-hat (methylation) x target centroid, r 0.40 | 0 of 8 | 0 of 8 | 1,2,1,2,3,1,3,0 | +0.298 | +0.291 | +0.360 | +0.007 (0.309) | 3 of 8 | 0 of 8 |
| *post-audit:* D + distLhat route (route r 0.15) | 5 of 8 | 5 of 8 | 5,4,0,3,0,5,7,4 | +0.256 | +0.243 | +0.287 | +0.014 (0.289) | 4 of 8 | 0 of 8 |
| *post-audit:* D + distLhat route (route r 0.25) | 8 of 8 | 7 of 8 | 7 in every replicate | +0.396 | +0.203 | +0.435 | +0.193 (0.004) | 8 of 8 | 5 of 8 |
| *post-audit:* D + E-hat x tc route (route r 0.25) | 0 of 8 | 0 of 8 | 2,3,3,1,2,2,2,1 | +0.353 | +0.346 | +0.387 | +0.008 (0.367) | 4 of 8 | 0 of 8 |

The "true r" of a planted truth is the LODO r of the planted signal alone on expected counts. For D + route truths it is the size of the route alone. D's part is N2's truth (beta 0.259, tau 0.35).

What the controls show, and where they stop:

* By the declared rule the search is not underpowered, because two declared truths were recovered in >= 50% of replicates. That conclusion rests on single-coefficient or low-dimensional terms. The 16 x 4 nmf16 x tg truth is recovered 2 of 8 times and loses to D (-0.069). GxL_tc is recovered 6 of 8 times but does not beat D: its oracle (+0.116) is itself below D, because the real signal underneath is D-like.
* **The engine cannot find a methylation route.** E-hat x tc was recovered in 0 of 16 replicates, whether alone at r 0.40 or added to D at r 0.25. Inside a fold, E-hat is re-estimated from fewer sources, so the planted route is only seen through a noisy copy.
* **Power for a realistic gain over D is low.** With D true plus a distance route worth about +0.044 in oracle r, the search recovered the route in 5 of 8 replicates but gained only +0.014 over D. No replicate reached p < 0.05. Only a large added route (oracle gain +0.23, nearly twice D's whole real signal) was beaten reliably, with p < 0.05 in 5 of 8 replicates.
* Only one effect size was planted per declared truth, and 8 replicates per truth gives wide uncertainty on every fraction above.

### Symbolic-regression layer (written here; no external package)

The catalogue holds 19,020 expressions of depth <= 2 over 12 variables:

* 3 learned gene axes, z1-z3, fitted on the current training sources only;
* measured unit location (ux, uy, uz) and target centroid (tx, ty, tz);
* the target gene axis tg1;
* -log distance to the target from the measured location (ld) and from the inferred location (ldh).

Unary ops are id, sq and tanh; binary ops are + - * / (protected). Expressions are selected by nested inner r_unit minus 0.005 x complexity.

| held-out division | formula chosen inside training divisions | complexity | inner | outer r | best leaf (same rule) | leaf outer r |
|---|---|---|---|---|---|---|
| amygdala-pallidum | `((z1 - uz) * tg1)` | 5 | +0.282 | +0.001 | `ldh` | -0.163 |
| hindbrain | `((z1 + z3) * tg1)` | 5 | +0.366 | -0.054 | `ldh` | +0.212 |
| hippocampal | `(tg1 * sq(z3))` | 4 | +0.322 | -0.018 | `ldh` | +0.036 |
| isocortex | `((ux * ldh) * tg1)` | 5 | +0.244 | -0.077 | `ldh` | +0.080 |
| midbrain | `((ty / uz) * z3)` | 5 | +0.280 | +0.013 | `ldh` | -0.100 |
| olfactory | `((z1 - z3) * tg1)` | 5 | +0.343 | +0.443 | `ldh` | +0.213 |
| thalamus | `((z1 - z3) * tg1)` | 5 | +0.298 | -0.126 | `ld` | -0.058 |

* The chosen formulas averaged mean7 +0.026, against +0.031 for the best leaf (division sign-flip p = 0.523). Under N0 (60 runs) the symbolic layer's outer mean7 averaged +0.014 (95th percentile +0.142), so the real p_N0 = 0.410. **No formula beats the simpler alternatives**: the declared rule fails on all three conditions.
* The recurring shape is (gene axis combination) x tg1, that is, this cell type's gene axis times the target's gene axis. Chosen on all 27 sources it is `(z1 - uz) * tg1`, coefficient -0.680 (bootstrap 95% interval -0.757 to -0.583), with a LODO selection score of 0.259 (optimistic, because it was selected on that score). It does not hold up on held-out divisions.

### Saved reusable modules (data/interim/d108_modules/, every output flagged INFERRED, provenance I)

| id | conversion | training data | held-out accuracy | uncertainty |
|---|---|---|---|---|
| M1_genes_to_location | 32-PC gene embedding to MERFISH CCF centroid (RBF kernel ridge) | 325 Allen subclasses with >= 50 MERFISH cells | out-of-fold R² 0.79 / 0.76 / 0.68 (x, y, z); RMSE 0.78 mm vs 1.65 mm for the mean | pred_sd = sqrt(bootstrap variance from 30 refits + OOF residual variance). It is conservative: 1.96 x pred_sd covers 99.2% of out-of-fold errors. The card gives a calibration factor of 0.80 for nominal 95% |
| M2_genes_to_division_composition | gene embedding to share of cells in 12 CCF divisions | same 325 subclasses | out-of-fold R² 0.67 overall (0.06 to 0.83 per division) | same form; coverage 97.3%, calibration factor 0.82 |
| S1_discovered_pipeline | the search rerun on all 27 sources: ae16 x tc + nmf16 x tg + gB10 x tg, lambda 10 | 27 non-HY sources, 187 units | **nested held-out mean7 -0.062** (worse than D); the selection score of 0.373 is optimistic | 200 bootstrap refits over sources, with selection fixed (mean coefficient sd 0.116 vs mean absolute coefficient 0.216). Marked **not predictive** |
| nested_search_oof_predictions.json (replaces the mislabelled S1_inferred_profiles.json) | out-of-division predictions of the **nested search** for every (source, subclass, target). Each value comes from the pipeline selected in the outer fold that held its division out (7 different pipelines), not from S1 | as above | that division's held-out r (mean7 -0.062); every record has `predictive: false` | sd = sqrt(variance over 20 bootstrap refits + mean squared out-of-division residual on the other 6 divisions). This covers selection and residual error, and 1.96 sd covers 95.2% of held-out values (the refit-only sd covered 27.5%). It still ignores E-hat re-estimation |
| SR formula (d108_sr_final.json) | `(z1 - uz) * tg1` | 27 sources | nested mean7 +0.026 (not better than one variable) | coefficient bootstrap sd 0.045 |

S1 can now be used on its own. The npz bundles the feature rows it uses (ae16, nmf16 and gB10 for all subclass names), the tc and tg target encodings, and the names and targets. `d108_prims.apply_bundle(npz, units, targets)` applies it, using within-source centring of the unit features and the target encoding. On all 27 sources it reproduces the engine's prediction exactly (maximum absolute difference 0). The card in `discovered_module_cards.json` documents its use. M1 and M2 are useful conversions in their own right, because genes predict location well. S1 should not be used for prediction in place of D.

### Sensitivity: scale constants inside the fold (audit item 11)

A few label-free scale constants (for Lm, rff_Lm, gD4xLm, distLm and distLhat) had been computed over all 27 sources, the outer held-out division included. The atlas and placement rule allows this. As a check, I reran the whole real-label nested search with these constants recomputed from the outer training sources only. It selected the same terms in all 7 folds. The outer r changed by at most 0.0023, and mean7 went from -0.0619 to -0.0616. The symbolic layer's variable scales (`d108_symbolic.scales`) are still computed over all 27 sources. They are label-free and were not rerun.

### Timings (4 workers, shared machine)

* **Main chain (13:15 to 14:09):**
  * Real nested search with comparators, ablations and bootstrap: 38 s wall.
  * Final module: 12 s.
  * Symbolic layer, real labels: 20 s.
  * N0, 60 searches: 16.8 min (3,904 s of summed fold time, about 65 s per nested search).
  * N2, 60 searches: 22.4 min (5,316 s).
  * Planted controls, 24 searches: 6.4 min.
  * Symbolic layer on N0, 60 runs: 7.0 min (1,653 s).
* **Post-audit chain (14:21 to 14:40):**
  * Planted controls, 56 searches (24 declared reruns + 32 new): 17.0 min (3,996 s of summed fold time).
  * Fold-scale sensitivity, 1 search: 26 s.
  * Final module rebuild: 13 s.
  * Extras (O2 comparator, prediction file, card calibration): 5 s.

### Disclosures

* After the declaration, two modes were added to `d108_discovery.py` without touching the search: `final` (refit on all 27 sources and save the module) and `report` (statistics).
  * `final` failed on its first run (13:15:57 to 13:16:09). Its all-sources fold used every subclass as held out, which left E-hat with no training subclasses. It was fixed to use an empty held-out side and rerun (13:17).
  * The fix was saved at 13:17:14, after N0 had started (13:16:29). N0 ran on the pre-edit file (its worker processes had already imported it). N2, the planted controls and SR-N0 ran on the edited file.
  * An independent audit reran real-label folds and N0/N2 replicates on the edited file and reproduced every per-fold result exactly, so the search was unchanged.
* **Post-audit code changes (14:17 to 14:40):**
  * New modes: `supp_declare`, `planted2`, `foldscale`, `extras`.
  * An oracle score for planted specs.
  * A self-contained S1 bundle (`final`) and `d108_prims.apply_bundle` / `Engine.scales_for`.
  * The single-term ablation convention and the audit numbers in `report`.
  * Edits to the report functions were made while the `planted2` workers were running. They had already imported the file, and those edits do not touch the search.
  * The 3 declared planted truths reproduced exactly (selected terms and outer r), and so did the S1 selection.
* The code default is 5 workers; every run used `D108_WORKERS=4` as declared.
* The pilot timing (129 s per nested search) was measured at 5 workers. The main runs averaged about 65 s (N0) and 89 s (N2) of summed fold time per nested search.
* **E-hat:** it is cross-fitted over subclasses, not left out by experiment (critic change 3 not adopted). This limits plate-fingerprint leakage but does not remove it. The bootstrap refits reuse the fold's E-hat, so the refit sd ignores E-hat uncertainty.
* **N2 truth:** it is D's in-sample fit on all 27 sources, not a cross-fitted profile.

### Response to the audit

| # | item | action |
|---|---|---|
| 1 | "What this means" overreached | Fixed. Measured selection optimism added (N0 inner +0.240 vs D's term +0.021), and the "closer to D's form" and "not a broken engine" claims removed |
| 2 | planted controls thin | Fixed in part: an E-hat truth, 3 D + route truths at realistic sizes, an oracle score and the expected-range clause (declared before running). Most truths still have one effect size, with 8 replicates each. Stated above |
| 3 | not a true typed-graph search | Wording fixed (fixed typed library). Composition operators were **not added**. Adding them would change the declared search after the result was seen. It would also add candidates to a search whose failure is already selection variance. This is left to a newly declared Phase 2 |
| 4 | S1_inferred_profiles mislabelled, sd too small | Fixed: renamed, `predictive: false` per record, and an sd that covers selection and residual error (coverage 95.2%) |
| 5 | S1 not reusable alone | Fixed: bundle plus `apply_bundle`, verified exact |
| 6 | M1/M2 pred_sd conservative | Fixed: coverage and calibration factors written into the cards (npz unchanged) |
| 7 | N2 truth in-sample | Disclosed |
| 8 | N2-tail wording | Fixed (1 of 60 N2 replicates lower, minimum -0.318) |
| 9 | edit timing | Disclosed with times |
| 10 | narrow comparators | O2 fitted on the same folds and T1 cited from the rehearsal; R_nested = K stated |
| 11 | global scale constants | Sensitivity run: same terms, at most 0.0023 change. The declared run stays the main result |
| 12 | E-hat not leave-experiment-out; bootstrap ignores E-hat | Not changed, because this would alter the declared search. Disclosed |
| 13 | lone term could not be a route | Fixed in `report`; no effect here |

### What this means

On these data, a search over gene bases, learned latents, location, inferred conversions and source x target interactions overfits the 6 training divisions. With about 200 term sets x 5 lambdas scored on 6 inner folds, the best candidate scores +0.24 inside the training divisions even on shuffled labels. That winner's curse is as large as D's whole real signal, and it explains the inner-outer gap on the real data. What the search prefers therefore says nothing about the form of the real signal. The failure is the engine's selection variance, not evidence about biology.

The planted controls bound what the engine can do. It reliably finds a large, simple distance chain. It does not find a methylation route. And it would not have shown a realistic improvement over D (+0.04 in r) even if one existed.

D remains the best tested model (O2 and T1's rehearsal score are slightly higher, but T1's edge did not survive its own held-out test in D-105). Next steps should reduce selection variance before widening the search:

* fewer, lower-variance candidates;
* a larger ADD_MIN;
* or a rule that admits a term only if it beats D inside the inner loop by more than the measured N0 optimism.

Composition operators and methylation interactions should wait until the engine can recover a realistic planted gain over D.
