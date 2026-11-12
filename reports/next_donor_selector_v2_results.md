# Next donor selector v2: DEVELOPMENT result (TRAIN nested CV)

**This is a development result, not an unseen test.** The 101 TRAIN brain samples have already informed the design.
The 25 validation and 40 test samples were not opened or scored. Any claim needs fresh independent data.

* **Protocol:** `reports/next_donor_selector_v2_protocol.md` (FINAL, SHA-256 `fbd05b06…`).
* **Code:** `analysis/next_donor_selector_v2.py` (`90567075…`). **Tests:** `tests/test_next_donor_selector_v2.py` (`5dfe55c2…`, 14 passed).
* All three were committed at 865ec91 before `RUN_STARTED`.
* **Run:** 2026-11-12 01:02–01:30 UTC, one continuous run, 4 workers, all four thread variables set to 1, `/opt/anaconda3/bin/python -I`.

## Result

**Verdict sentence (predeclared rule): "development gate met, not resolved".**

* **Point gate:** the selector has lower distance and higher F1 than every control, and its switch rate (0.435) is
  above 0.
* **E1, E2 and adaptive fallback:** the selector beats all three on both measures, but by very small margins.
  * Distance: the gains are tiny (−3.3, −0.1 and −1.8 µm) and every distance interval crosses zero.
    Against E2 the distance is effectively a tie: 605.9 vs 606.0 µm.
  * F1: the gain is small (+0.005 to +0.007). Its interval excludes zero against E1 and the adaptive fallback,
    but not against E2.
* **Plain reading:** the learned choice adds about +0.006 F1 and no meaningful distance gain over the local medoids.
  The large gains are against the copies and the random or nearest switches, not against E1/E2.

4,368 neurons, 101 brain samples. Neurons are averaged within each sample, then samples are averaged equally.

| method | distance µm | F1 at 200 µm |
|---|---:|---:|
| selector_v2 | 605.9 | 0.2759 |
| E1_local_uniform | 609.2 | 0.2688 |
| E2_local_kernel | 606.0 | 0.2705 |
| adaptive_fallback (inner-selected E1/E2, unchanged) | 607.7 | 0.2697 |
| nearest / medoid / peptide_nearest / peptide_medoid | 783.4 / 705.5 / 777.6 / 773.0 | 0.2216 / 0.2014 / 0.2130 / 0.1699 |
| local_nearest / local_label_nearest | 783.4 / 760.1 | 0.2222 / 0.2262 |
| switch_nearest / switch_random | 689.3 / 688.7 | 0.2518 / 0.2479 |
| permuted_labels | 758.4 | 0.2266 |
| oracle best of 16 (truth-informed ceiling, descriptive) | 397.0 | 0.3838 |

Selector minus control. Brackets are paired brain-sample bootstrap 95% intervals (10,000 draws). Negative distance
and positive F1 favour the selector.

| control | distance µm | F1 | samples better on both |
|---|---|---|---:|
| E1_local_uniform | −3.3 (−10.2 to +3.7) | +0.0070 (+0.0020 to +0.0115) | 45 / 101 |
| E2_local_kernel | −0.1 (−7.8 to +7.9) | +0.0054 (−0.0008 to +0.0108) | 45 / 101 |
| adaptive_fallback | −1.8 (−9.4 to +6.2) | +0.0062 (+0.0005 to +0.0111) | 46 / 101 |
| nearest | −177.5 (−200.2 to −155.2) | +0.0543 (+0.0461 to +0.0620) | 92 |
| medoid | −99.6 (−130.8 to −70.0) | +0.0745 (+0.0610 to +0.0875) | 83 |
| peptide_nearest | −171.7 (−201.6 to −142.7) | +0.0629 (+0.0512 to +0.0747) | 87 |
| peptide_medoid | −167.1 (−204.7 to −132.5) | +0.1060 (+0.0919 to +0.1204) | 90 |
| local_nearest | −177.5 (−200.2 to −155.3) | +0.0536 (+0.0453 to +0.0617) | 92 |
| local_label_nearest | −154.2 (−179.9 to −130.3) | +0.0497 (+0.0397 to +0.0597) | 88 |
| switch_nearest | −83.4 (−100.3 to −68.5) | +0.0241 (+0.0188 to +0.0290) | 88 |
| switch_random | −82.8 (−99.2 to −67.2) | +0.0280 (+0.0226 to +0.0339) | 83 |
| permuted_labels | −152.5 (−177.8 to −128.9) | +0.0493 (+0.0406 to +0.0579) | 93 |

## Switching, inner choices, ablation

* **Combined-rule branches:**
  * switched: 1,901 (0.435);
  * top candidate was the fallback: 1,836;
  * top candidate not predicted better on both heads: 631.
  * The permuted-label control switched on 0.632 of queries.
* **Switched subset** (1,901 neurons, neuron means): selector 576.3 µm / 0.2936. On the same neurons:
  * the adaptive fallback scores 593.9 µm / 0.2691;
  * E1 592.0 / 0.2699;
  * E2 591.7 / 0.2712;
  * switch_nearest 765.0 / 0.2292;
  * switch_random 779.6 / 0.2185.
* **Inner choices by outer fold:**
  * fallback: E2, E1, E1, E2, E1;
  * learner: HGB in all five folds.
* **Fit-only door:** 0 unexpected refusals. Each outer fold's probe read of an outer-eval label was refused.
* **Fixed ablation, attachment effect:** attached minus translated is +0.5 to +0.6 µm in distance and +0.004 to
  +0.005 F1, for selector, E1 and E2.
* **Fixed ablation, donor-choice effect** (both translated):
  * vs E1: −3.1 µm (−10.0 to +3.6), +0.0082 F1 (+0.0041 to +0.0122);
  * vs E2: +0.0 µm (−7.3 to +7.5), +0.0056 F1 (+0.0002 to +0.0105).
  * So the small F1 gain comes from donor choice, not from the attachment.
* **Per-label coverage:** every source label had switches, at rates 0.34 to 0.68. Small labels (fewer than 5 samples)
  are descriptive only. The table is in `tables.md`.

## Executed geometry (counts of 4,368)

| method | attach 500 µm | attach 1,000 µm | translation |
|---|---:|---:|---:|
| selector_v2 | 2,975 | 1,333 | 60 |
| E1_local_uniform | 2,884 | 1,416 | 68 |
| E2_local_kernel | 3,083 | 1,221 | 64 |
| adaptive_fallback | 2,969 | 1,333 | 66 |
| local_nearest | 4,258 | 103 | 7 |
| switch_nearest | 3,532 | 795 | 41 |
| switch_random | 2,928 | 1,376 | 64 |
| permuted_labels | 2,798 | 1,499 | 71 |
| the four copies | 0 | 0 | 4,368 each |

**Wording limitation (disclosed, protocol not edited during the run).** The protocol says distal coordinates are
"preserved exactly". That holds only when the compact attachment succeeds (500 or 1,000 µm). The translation
fallback moves the whole arbor rigidly to the table soma. That happened for 60 selector forecasts, 60–71 for the
local methods, and every copy.

## Caveats

* Development only: TRAIN informed the design, so this is not evidence of generalisation.
* The distance improvement over E1/E2 is not established: all three intervals cross zero, and against E2 it is a tie.
* The F1 improvement is small and not resolved against E2.
* The oracle gap (397 µm) shows room remains, but the permitted inputs recover little of it.
* Labels are exact recomputations for the emitted geometry. Cached D was used only for the medoids and centrality.
  Parity with the engine's committed development forecasts: E1/E2 donors 919/919, copy donors 3,676/3,676, and 240
  of 240 arrays checked, with bit-identical metrics.

## Files (`data/interim/next_donor_selector_v2/`, not committed: data/interim is git-ignored)

| file | SHA-256 |
|---|---|
| RUN_STARTED.json | `54c11868f10d1ef9c0d78183837d139c5884819212d20034272c42ed81e96046` |
| parity.json | `f458f0090c20e27330af9570cb98cab09f03e97b5417aae6b69013c0e2289761` |
| folds.json | `f97a99339e0ee952217d959b721b51f53bedd968388de30dcf609d9dc7e5ca10` |
| forecasts.json (sealed) | `50be784179be2c26e0b61db77ba7699bf94ae5eb81f3f8c6489ed5e34631ed9a` |
| FORECASTS_SEALED.json | `6ba7eeb8902118f1f40ffb7657b7eb12ecb62450431ac2b1e03df2e1fac25f97` |
| results.json | `e38752e0d567d9719ab8c1530efd7325cc3a5ee3c17e4d2f54d8f1f2422e4e9f` |
| tables.md | `6ff9a137e32025524028bf81a39c21d821331a7dcdb823667a2e584e6f3117d8` |
| RUN_COMPLETE.json | `c20f7113057042b784d805634d7ac7203583b7709f2107b83fb80ddb654bbcde` |

Also in that folder: `labels_fit.jsonl` (186,757 exact fit labels), `labels_score.jsonl` (91,941 outer scoring pairs),
`run.log` (progress), and `donor_summaries.npz`.
