# Axon completion v2: DEVELOPMENT result (TRAIN, synthetic 300 um masking)

**Verdict (predeclared gate): FAIL: development completion gate not met.** ROOT checked it independently and
confirmed it.

**Scope:**
* This is one error-informed development follow-up to v1, on the same 101 TRAIN brain samples, which are already
  spent. It is not whole-neuron or synaptic-circuit prediction.
* There is no independent confirmation.
* Validation, test and external data were not opened.
* There is no viewer or website promotion.

**Run details:**
* Protocol: `reports/axon_completion_v2_protocol.md` (FINAL/APPROVED, SHA-256 `48d3a777…`).
* Code: `analysis/axon_completion_v2.py` (`f75712f5…`). Tests: `tests/test_axon_completion_v2.py` (`496af39f…`).
* All were committed at 4a1a7ce before `RUN_STARTED`.
* The run: 2026-11-12 05:26–05:36 UTC, one continuous run, 4 workers, all four thread variables set to 1.

## What changed from v1

The v2 candidate added two things at once:
* fit-only donor geometry landmarks: 32 per outer fold, via the verified cache D. Landmark columns that share the
  query's sample are masked;
* a third residual head for the cable-error gain. A switch also needs a predicted cable-error gain of 0 or more.

The mask, decoder, folds, pools, labels and the inner-chosen fallbacks are v1's, unchanged. The two additions are
**bundled**, so this result does not attribute any effect to either one alone.

**Labels and scoring:**
* All fitting labels were reused from v1 (cache plan: 0 missing pairs, and 4,368 of 4,368 local/control parity).
* After the forecast seal, only 32 evaluation pairs were missing from both v1 caches, and only those were scored.
* Every score matched its sealed completion hash.

## Result

4,368 neurons, 101 samples. Neurons are averaged within each sample, then the samples equally.

| method | hidden distance um | hidden F1 | invalid | share > 2x cable | mean abs rel cable error |
|---|---:|---:|---:|---:|---:|
| **completion_selector (v2)** | 623.3 | 0.2636 | 0 | 0.262 | 2.013 |
| completion_selector_v1 (frozen) | 620.4 | 0.2699 | 0 | 0.285 | 2.297 |
| E1_local_uniform | 619.7 | 0.2599 | 0 | 0.293 | 2.153 |
| E2_local_kernel | 620.2 | 0.2605 | 0 | 0.292 | 2.154 |
| adaptive_fallback | 622.1 | 0.2593 | 0 | 0.293 | 2.197 |
| metadata_only | 617.0 | 0.2613 | 0 | 0.291 | 2.155 |
| local_prefix_nearest | 762.0 | 0.2240 | 0 | 0.230 | 2.023 |
| prefix_abs_nearest | 769.1 | 0.2217 | 0 | 0.232 | 2.027 |
| prefix_centred_nearest | 1103.2 | 0.1100 | 0 | 0.241 | 2.159 |
| prefix_peptide_nearest | 799.8 | 0.2083 | 0 | 0.224 | 1.800 |
| nearest / medoid | 799.8 / 687.1 | 0.2088 / 0.2213 | 0 / 0 | 0.229 / 0.321 | 1.939 / 2.735 |
| peptide_nearest / peptide_medoid | 823.4 / 749.5 | 0.1971 / 0.1896 | 0 / 0 | 0.217 / 0.429 | 1.685 / 2.891 |
| switch_random (v2 switches) | 669.7 | 0.2480 | 0 | 0.287 | 2.266 |
| permuted_labels (v2) | 622.1 | 0.2593 | 0 | 0.293 | 2.197 |
| observation_only (diagnostic, not required) | +inf | 0.0000 | 4,368 empty | 0.000 | 1.000 |

v2 selector minus control. Brackets are paired sample-bootstrap 95% intervals (10,000 draws). The gate needs distance
≤ −10 um and F1 ≥ +0.010, with both intervals excluding zero.

| control | distance um | F1 | samples better on both |
|---|---|---|---:|
| completion_selector_v1 | +3.0 (−3.4 to +9.4) | −0.0064 (−0.0140 to −0.0010) | 29 / 101 |
| E1_local_uniform | +3.6 (−6.6 to +18.3) | +0.0037 (−0.0013 to +0.0083) | 45 |
| E2_local_kernel | +3.1 (−6.8 to +17.4) | +0.0031 (−0.0020 to +0.0079) | 40 |
| adaptive_fallback | +1.2 (−8.7 to +15.6) | +0.0043 (−0.0007 to +0.0089) | 47 |
| metadata_only | +6.3 (−4.5 to +22.2) | +0.0023 (−0.0037 to +0.0080) | 37 |
| local_prefix_nearest | −138.7 (−164.2 to −109.9) | +0.0396 (+0.0327 to +0.0463) | 88 |
| prefix_abs_nearest | −145.8 (−169.3 to −120.6) | +0.0419 (+0.0346 to +0.0489) | 88 |
| prefix_centred_nearest | −479.9 (−513.6 to −446.2) | +0.1536 (+0.1409 to +0.1656) | 99 |
| prefix_peptide_nearest | −176.5 (−202.9 to −151.6) | +0.0553 (+0.0460 to +0.0652) | 92 |
| nearest / medoid | −176.5 / −63.8 | +0.0547 / +0.0423 | 97 / 67 |
| peptide_nearest / peptide_medoid | −200.1 / −126.2 | +0.0665 / +0.0740 | 92 / 82 |
| switch_random | −46.4 (−61.5 to −33.7) | +0.0155 (+0.0103 to +0.0216) | 75 |
| permuted_labels | +1.2 (−8.7 to +15.6) | +0.0043 (−0.0007 to +0.0089) | 47 |

**Why it fails:**

1. **No resolved win over the local medoids or the metadata-only learner.** Distance differences are +1.2 to
   +6.3 um, and every distance and F1 interval crosses zero.
2. **Worse than the frozen v1 selector on F1:** −0.0064, with the interval below zero. Distance is a tie.
3. **The cable guard, part one (share above 2x) fails.** The v2 share is 0.262, which is more than 0.02 above the
   copy and prefix controls (0.217–0.241).
4. **The cable guard, part two (mean absolute relative cable error) passes.** v2's 2.013 is below the adaptive
   fallback's 2.197.

**Plain reading:**
* The v2 bundle gave a lower mean absolute relative cable error and fewer forecasts above 2x cable than the v1
  selector and the local medoids. This does not mean it emits less cable in every case.
* It switched less often: 718 of 4,368 (0.164), against v1's 1,166.
* It lost v1's small F1 gain and gained no distance.
* The development evidence does not support adopting it.

## Details

* **Fallback chosen by inner groups:** E1, E1, E2, E1, E2 for folds 0–4, the same as v1 (not reselected).
* **Fit-only door:** each fold's probe read of an outer-eval label was refused (1 per fold). No invalid training
  labels were skipped.
* **Validity:** 0 invalid forecasts in every required method.
* **Shared-trunk merge:** it removed 131 um of exactly duplicated cable in total over the v2 selector's forecasts.
* **Permuted control:** the v2 permuted-labels control never switched, so its forecasts equal the adaptive fallback.

## ROOT verification

ROOT independently checked the saved outputs with `analysis/axon_completion_v2_result_check.py` (SHA-256
`cad5a7a5fa6368d44b551f975dbbebea69fec67f2186fa778279b11c50592fed`). ROOT corrected one output-path assumption in
it: the folds are read from the pinned v1 `folds.json`, which v2 reuses.
* It checked:
  * 4,368 queries in 101 samples;
  * 74,256 score-hash matches and 37,375 selected pairs;
  * 3,650 exact fallbacks.
* The means, bootstrap intervals and gate=false agree, with a maximum error of 2.3e-13.
* The old test is unchanged.

Receipt: `data/interim/axon_completion_v2_root/saved_result_check.json`, SHA-256
`fc3c31b19d76d3bba1726519661a641c1c37f25b1c311571573857b71a8b22ed`.

## Limits

* **Development only:** spent TRAIN samples, a synthetic 300 um mask, and isolation that is code-guarded per fit
  rather than physical absence.
* **Bundled change:** the landmarks and the cable head were not tested separately, so no component attribution is
  made.
* **Closed:** the experiment is closed and will not be repeated or tuned.

## Files (`data/interim/axon_completion_v2/`, git-ignored)

| file | SHA-256 |
|---|---|
| plan.json | `e067fb60d97086f1314a022692e3b18f936cd7932e2464bbedef53406a2b5ba2` |
| RUN_STARTED.json | `89ab32dc98fa680337d0ed6e4c5f1d456e28169b838f2bc38d3852ad3cc3b5e8` |
| fold_info.json | `4546c0276f544062bfcf29b6733763f3a49e46b74b15e6a4c230df4216204748` |
| forecasts.json (sealed) | `60d8a6b4521ad5e3a97f018e9256a34941f886187c5d448af0b8c0f0ca7b40fb` |
| FORECASTS_SEALED.json | `362a9ec5cbbbf429fc3b8d55b9ef87ceddf83eac5830169a98a4711262c750ee` |
| cache_coverage.json | `d9853f25748e18df6d578fd026188f0c30379d15bc4a45048f10d049f27588c0` |
| results.json | `7e0d8312c54e03f3cae83d7704e3951c8f49a06c1d24ee9cdf1cdd6fb9b664ba` |
| RUN_COMPLETE.json | `ccf15bd40d78c04268060bba014bd661caa5772a6506ce9c73f24c2f20652f6d` |

The same folder also holds `outer_0..4.json` with their completion receipts, `decoded_new.jsonl`,
`labels_missing.jsonl` (32 pairs) and `run.log`.
