# Axon completion v1: DEVELOPMENT result (TRAIN, synthetic 300 um masking)

**Verdict (predeclared gate): FAIL: development completion gate not met.**

**Scope:**
* This is synthetic-masking development completion on the 101 TRAIN brain samples, which have already informed
  earlier designs. It is not whole-neuron or synaptic-circuit prediction.
* Validation and test samples were not opened.
* There is no viewer or website promotion.

**Run details:**
* Protocol: `reports/axon_completion_v1_protocol.md` (FINAL, SHA-256 `698a8204…`).
* Code: `analysis/axon_completion_v1.py` (`27cc4250…`). Tests: `tests/test_axon_completion_v1.py` (`26745423…`).
* ROOT's contract `13e6cd60…` and preflight `0bc0693a…`.
* All were committed at b9ecc2b before `RUN_STARTED`. 37 tests passed.
* The run: 2026-11-12 03:57–04:58 UTC, one continuous run, 4 workers, all four thread variables set to 1,
  `/opt/anaconda3/bin/python -I`.

## Cohort and cost

* **Cohort:** 4,368 neurons in 101 samples, all eligible. The fixed fully-observed rule excluded 0, as in ROOT's
  structural preflight.
* **Table-minus-trace soma offset:** median 0.002 um, 99th percentile 2.4 um, maximum 4.8 um.
* **Work done:**
  * 185,293 fit label pairs, at 75 pairs/s (46 min);
  * 37,423 scored forecast pairs after the seal, at 78 pairs/s;
  * prefix costs: 206 s, in memory only;
  * mean candidate union: 27.2 donors.

## Result

Hidden-cable scores come from the contract `score_completion` (new cable only, 50 um sampling, 200 um tolerance).
Neurons are averaged within each sample, then the samples equally.

| method | hidden distance um | hidden F1 | invalid | inflated (>2x) | mean abs rel cable error |
|---|---:|---:|---:|---:|---:|
| **completion_selector** | 620.4 | 0.2699 | 0 | 0.285 | 2.297 |
| E1_local_uniform | 619.7 | 0.2599 | 0 | 0.293 | 2.153 |
| E2_local_kernel | 620.2 | 0.2605 | 0 | 0.292 | 2.154 |
| adaptive_fallback (E1/E2, inner-chosen) | 622.1 | 0.2593 | 0 | 0.293 | 2.197 |
| metadata_only (learned, no prefix input) | 617.0 | 0.2613 | 0 | 0.291 | 2.155 |
| local_prefix_nearest (matched same pool) | 762.0 | 0.2240 | 0 | 0.230 | 2.023 |
| prefix_abs_nearest (full pool) | 769.1 | 0.2217 | 0 | 0.232 | 2.027 |
| prefix_centred_nearest (full pool) | 1103.2 | 0.1100 | 0 | 0.241 | 2.159 |
| prefix_peptide_nearest | 799.8 | 0.2083 | 0 | 0.224 | 1.800 |
| nearest / medoid | 799.8 / 687.1 | 0.2088 / 0.2213 | 0 / 0 | 0.229 / 0.321 | 1.939 / 2.735 |
| peptide_nearest / peptide_medoid | 823.4 / 749.5 | 0.1971 / 0.1896 | 0 / 0 | 0.217 / 0.429 | 1.685 / 2.891 |
| switch_random | 688.5 | 0.2453 | 0 | 0.289 | 2.269 |
| permuted_labels | 622.1 | 0.2593 | 0 | 0.293 | 2.197 |
| observation_only (diagnostic, not required) | +inf | 0.0000 | 4,368 empty | 0.000 | 1.000 |

Selector minus control. Brackets are paired sample-bootstrap 95% intervals (10,000 draws). The gate needs distance
≤ −10 um and F1 ≥ +0.010, with both intervals excluding zero.

| control | distance um | F1 | samples better on both |
|---|---|---|---:|
| E1_local_uniform | +0.6 (−12.0 to +18.9) | +0.0101 (+0.0040 to +0.0177) | 46 / 101 |
| E2_local_kernel | +0.2 (−12.1 to +18.1) | +0.0094 (+0.0034 to +0.0172) | 44 |
| adaptive_fallback | −1.8 (−14.3 to +16.4) | +0.0107 (+0.0049 to +0.0182) | 47 |
| metadata_only | +3.3 (−10.0 to +23.0) | +0.0086 (+0.0016 to +0.0177) | 45 |
| local_prefix_nearest | −141.6 (−167.5 to −111.6) | +0.0459 (+0.0385 to +0.0545) | 92 |
| prefix_abs_nearest | −148.8 (−172.3 to −123.6) | +0.0483 (+0.0403 to +0.0571) | 90 |
| prefix_centred_nearest | −482.9 (−516.1 to −449.8) | +0.1599 (+0.1475 to +0.1720) | 100 |
| prefix_peptide_nearest | −179.5 (−206.3 to −154.1) | +0.0617 (+0.0509 to +0.0734) | 95 |
| nearest / medoid | −179.5 / −66.8 | +0.0611 / +0.0486 | 98 / 70 |
| peptide_nearest / peptide_medoid | −203.0 / −129.1 | +0.0729 / +0.0804 | 91 / 83 |
| switch_random | −68.1 (−85.8 to −51.5) | +0.0247 (+0.0178 to +0.0324) | 84 |
| permuted_labels | −1.8 (−14.3 to +16.4) | +0.0107 (+0.0049 to +0.0182) | 47 |

**Why it fails:**

1. **No distance gain over the local medoids.** Against E1, E2, the adaptive fallback, the metadata-only learner and
   the permuted control, the distance difference is −1.8 to +3.3 um, and every interval crosses zero.
2. **The F1 gain is resolved but small.** It is +0.009 to +0.011, with intervals above zero. It only reaches the
   +0.010 margin against E1 and the adaptive fallback.
3. **The cable guard fails, on both parts:**
   * the selector's inflated share is 0.285, which is more than 0.02 above the shares of the copy and prefix controls
     (0.217–0.241);
   * its mean absolute relative cable error is 2.297, higher than the adaptive fallback's 2.197.
   * So part of its F1 gain may come from choosing donors with more cable.

**Plain reading:**
* The prefix-aware COPYING controls are much worse than the local soma medoids (about 760–1,100 um vs 620 um).
  Matching the first 300 um does not, on its own, pick a donor whose distal route is right.
* The learned selector adds about +0.01 hidden F1 and no distance gain over E1/E2.
* Its prefix input adds little beyond the metadata-only learner: +0.009 F1, distance +3.3 um (a tie).

## Details

* **Fallback chosen by inner groups:** E1, E1, E2, E1, E2 across folds 0–4. The inner joint losses ranked E1 ≈ E2
  below all four prefix-aware rules in every fold. The metadata-only fallback matched in every fold.
* **Switches:**
  * selector: 1,166 of 4,368 (0.267);
  * metadata_only: 979;
  * permuted_labels: 0. It never predicted a gain on both heads, so its forecasts equal the adaptive fallback exactly.
* **Fit-only door:** each outer fold's probe read of an outer-eval label was refused (1 per fold), with no other
  refusal. No invalid training labels were skipped.
* **Validity:** 0 invalid forecasts in every required method. Only the diagnostic observation-only forecasts are
  empty failures.
* **Shared-trunk merge:** the merge removed 0–856 um of exactly duplicated cable per method over all 4,368 forecasts
  (selector 138 um). Without it, those few forecasts would have been duplicate-cable failures.
* **Executed grafts (components):**

  | method | attach 500 | attach 1,000 | rigid |
  |---|---:|---:|---:|
  | selector | 2,698 | 4,527 | 13,683 |
  | E1 | 2,104 | 4,402 | 13,850 |
  | prefix_abs_nearest | 3,798 | 2,634 | 9,632 |

  Most components are relocated rigidly, so their distal coordinates move.
* **Cable:** for the local medoids and the selector, the median new/hidden cable ratio is about 1.27. The copy and
  prefix-nearest controls are about 1.0.

## Limits

* **Development only:** the TRAIN samples informed the design, so this is not evidence of generalisation.
* **The masking is synthetic:** the 300 um cut is a fixed rule, not a real partial reconstruction.
* **Isolation is code-guarded per fit:** all TRAIN outcomes are on disk and loadable by the workers. They were not
  physically absent.
* **The F1 tolerance is lenient:** at 200 um, a forecast can gain F1 by emitting more cable, which is why the cable
  guards exist.

## Files (`data/interim/axon_completion_v1/`, git-ignored)

| file | SHA-256 |
|---|---|
| RUN_STARTED.json | `743654b4a54f8adcf3425f6a43268a81ae47956c6b127539a7cfb831abe7fc35` |
| cohort.json | `5b17eff65fd4b5aa54b976fe17926da204d890e959388a80662e8b06d761918f` |
| folds.json | `0c0413f50690899cc8f6f72ce996887ee6ff79f038caade86755a13f66c47e8d` |
| plan_counts.json | `75989210e1c23ba969c34db21ea2e7a72c55ac91dc407caa6141379c78d7a7e7` |
| fold_info.json | `12d5423e74a1aca607a9a9c4eeae0786fb54846207571ff027803ede3e604916` |
| forecasts.json (sealed) | `118c94184997f202d1ad4350a16ff0680d3e4fa44e27b26743c2a4fef5f6cfa0` |
| FORECASTS_SEALED.json | `db83d7a2b819eec5d84ceca68172df8975b0b813fdc2653bd806e3e603b64b99` |
| results.json | `bc76687b0bdd99063dca2aa82a089bf2522cec21afdd268c8fed5bd5773ba428` |
| RUN_COMPLETE.json | `bcf7e65184033cee175aad1dc7880cc8e2eb0dd0489e8ace8c0483ed24e67865` |

The same folder also holds `index.jsonl`, `labels_fit.jsonl` (185,293 pairs), `labels_score.jsonl` (37,423 pairs),
`decoded.jsonl` and `run.log`. The experiment is complete and is not repeated.
