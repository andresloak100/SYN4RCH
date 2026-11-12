# Nonlinear gene-and-position model: results (development and synthetic only)

## Headline

* **The fixed shallow nonlinear model (HGBT, gene + depth) beats every copying, averaging and depth-only control
  on the exposed A1 data.**
* **It loses to the frozen logistic model on both log loss (0.3251 vs 0.3190) and Brier (0.0960 vs 0.0952).**
* **The A1 development gate was therefore not met. The frozen logistic model stays the operational model; the
  nonlinear model is archived as experimental and not promoted.**

On 12 fresh synthetic cases, the nonlinear model met the gate in 2. Both were batch-shift/missing cases, and
even that scenario passed only 2 of 3. Its failures are informative:
* **location only:** a depth-only nonlinear model wins;
* **complementarity:** matched-input KNN wins;
* **hidden latent signal:** the constant 0.5 (the null) wins.

**Scope:**
* The real data are 598 cortical cells from 2 mice (XC95, XC97), both already seen during development.
* The outputs are probabilities for 11 captured projection targets.
* No hypothalamic, 3D-geometry, synapse or new-animal claim is made, and there is no cell-level significance
  testing.

## What was run

**Model:**
* One scikit-learn 1.8 `HistGradientBoostingClassifier` per target, with one fixed structure.
* 100 trees, learning rate 0.05, max depth 2, max 4 leaves, min 20 samples per leaf, L2 1.0.
* No early stopping, no class weights, seed 20261112.
* Same training-only preprocessing as the earlier work (median fill, then scaling), group-balanced weights, and
  shrinkage `p' = (1 - α) p + α q_fit` toward the training animals' capture rate.

**Families:**
* gene + depth: log1p of 23 measured genes plus measured laminar depth;
* depth only, which shows what the genes add beyond nonlinear position.

**Selection:**
* α from {0, 0.1, 0.25} only, one global α per family.
* Chosen by target-equal log loss under that family's own α=0 Brier guard.
* No structure search and no retuning.

**Selected:** gene + depth **α = 0.1**; depth only **α = 0.25**. They were frozen before any simulation outcome
was generated.

**Comparators:** the 13 frozen methods from the earlier work:
* the 10 selected recipe controls;
* the original logistic model;
* two raw 1NN copy anchors.

On A1 they were recomputed fold by fold, and all 85,514 keyed out-of-fold probabilities matched the saved
development predictions exactly (max difference 0). The all-animal deployment artifact was never scored as
out-of-fold.

**Folds:** leave one animal out (two folds), with exact frozen cell IDs, targets and masks.

## Stage 1: A1 real development (both mice exposed)

### Every method

Ordered by log loss. Each is the target-equal mean of the animal-equal scores, with each mouse's own macro.

| method | log loss | Brier | XC95 log loss / Brier | XC97 log loss / Brier |
|---|---:|---:|---|---|
| frozen logistic, operational (C=0.1, α=0.1) | 0.3190 | 0.0952 | 0.3461 / 0.1016 | 0.2920 / 0.0888 |
| **gene + depth HGBT (selected α=0.1)** | 0.3251 | 0.0960 | 0.3536 / 0.1013 | 0.2966 / 0.0907 |
| matched KNN (k=31, α=0.1) | 0.3494 | 0.1035 | 0.3895 / 0.1119 | 0.3092 / 0.0950 |
| depth logistic | 0.3509 | 0.1080 | 0.3773 / 0.1150 | 0.3244 / 0.1010 |
| depth KNN (k=31, α=0.5) | 0.3535 | 0.1091 | 0.3785 / 0.1156 | 0.3285 / 0.1027 |
| **depth-only HGBT (selected α=0.25)** | 0.3539 | 0.1092 | 0.3810 / 0.1162 | 0.3268 / 0.1021 |
| original logistic (C=1, α=0) | 0.3561 | 0.1003 | 0.4079 / 0.1085 | 0.3043 / 0.0920 |
| molecular KNN (k=31, α=0.1) | 0.3573 | 0.1062 | 0.3989 / 0.1149 | 0.3156 / 0.0976 |
| training frequency | 0.3733 | 0.1172 | 0.3988 / 0.1239 | 0.3478 / 0.1105 |
| technical logistic | 0.3739 | 0.1162 | 0.4012 / 0.1228 | 0.3467 / 0.1097 |
| copy-1NN matched (α=0.5) | 0.3975 | 0.1242 | 0.4297 / 0.1304 | 0.3652 / 0.1179 |
| copy-1NN depth (α=0.5) | 0.4163 | 0.1317 | 0.4506 / 0.1390 | 0.3819 / 0.1244 |
| constant 0.5 | 0.6931 | 0.2500 | 0.6931 / 0.2500 | 0.6931 / 0.2500 |
| copy-1NN matched, raw | 2.5682 | 0.1859 | 2.5904 / 0.1875 | 2.5461 / 0.1843 |
| copy-1NN depth, raw | 2.7530 | 0.1993 | 2.7805 / 0.2013 | 2.7255 / 0.1973 |

### HGBT grid (raw = α 0, plus the shrinkage options)

| configuration | log loss | Brier | XC95 log loss | XC97 log loss |
|---|---:|---:|---:|---:|
| `hgbt_gene_position|alpha=0.0` (raw, α=0) | 0.3276 | 0.0961 | 0.3589 | 0.2963 |
| `hgbt_gene_position|alpha=0.1` | 0.3251 | 0.0960 | 0.3536 | 0.2966 |
| `hgbt_gene_position|alpha=0.25` | 0.3259 | 0.0968 | 0.3524 | 0.2995 |
| `hgbt_position|alpha=0.0` (raw, α=0) | 0.3601 | 0.1109 | 0.3901 | 0.3301 |
| `hgbt_position|alpha=0.1` | 0.3564 | 0.1099 | 0.3849 | 0.3279 |
| `hgbt_position|alpha=0.25` | 0.3539 | 0.1092 | 0.3810 | 0.3268 |

### Selected gene + depth HGBT against every comparator

Control minus candidate: positive means the nonlinear model is better.

| comparator | log-loss diff | Brier diff | targets: candidate lower log loss / Brier |
|---|---:|---:|---|
| frozen logistic, operational (C=0.1, α=0.1) | -0.0060 | -0.0008 | 3 / 3 of 11 |
| original logistic (C=1, α=0) | +0.0310 | +0.0043 | 11 / 10 of 11 |
| constant 0.5 | +0.3681 | +0.1540 | 11 / 11 of 11 |
| copy-1NN depth, raw | +2.4279 | +0.1032 | 11 / 11 of 11 |
| copy-1NN depth (α=0.5) | +0.0912 | +0.0357 | 10 / 11 of 11 |
| copy-1NN matched, raw | +2.2432 | +0.0899 | 11 / 11 of 11 |
| copy-1NN matched (α=0.5) | +0.0724 | +0.0281 | 9 / 11 of 11 |
| **depth-only HGBT (selected α=0.25)** | +0.0288 | +0.0131 | 7 / 7 of 11 |
| depth KNN (k=31, α=0.5) | +0.0284 | +0.0131 | 6 / 6 of 11 |
| depth logistic | +0.0258 | +0.0120 | 6 / 5 of 11 |
| matched KNN (k=31, α=0.1) | +0.0243 | +0.0074 | 9 / 5 of 11 |
| molecular KNN (k=31, α=0.1) | +0.0322 | +0.0102 | 9 / 5 of 11 |
| technical logistic | +0.0489 | +0.0202 | 7 / 5 of 11 |
| training frequency | +0.0482 | +0.0211 | 7 / 5 of 11 |

**Gate** (lower log loss and no higher Brier than the frozen logistic model, the depth-only HGBT and every strong
copy/averaging control): **NOT MET.**
* The only failure is the frozen logistic model, which is better on both metrics.
* By mouse, the nonlinear model has slightly better Brier in XC95 (0.1013 vs 0.1016) but worse log loss in both
  mice.
* Every depth-only, averaging and copying requirement is met.

## Stage 2: 12 fresh synthetic cases (generic sensitivity, not biology)

**Setup:**
* The frozen generator, unchanged: 4 scenarios × 3 new prespecified seeds.
* Each case has 6 training and 4 test simulated animals, 100 cells per animal, 16 genes and 4 targets.
* Every method was refit on the training animals only. Test labels, masks and oracle values never reached any
  fitter.
* The settings were the ones frozen from stage 1; none were chosen from these outcomes.

### Each case

| case (scenario · seed) | candidate log loss | Brier | constant 0.5 log loss | best log loss (method) | gate | failed against |
|---|---:|---:|---:|---|---|---|
| batch_shift_missing · 876241 | 0.6355 | 0.2217 | 0.6931 | 0.6355 (**gene + depth HGBT (selected α=0.1)**) | **met** | — |
| batch_shift_missing · 876242 | 0.6464 | 0.2272 | 0.6931 | 0.6464 (**gene + depth HGBT (selected α=0.1)**) | **met** | — |
| batch_shift_missing · 876243 | 0.6181 | 0.2144 | 0.6931 | 0.6121 (molecular KNN (k=31, α=0.1)) | not met | matched KNN, molecular KNN |
| complementarity · 876221 | 0.5692 | 0.1926 | 0.6931 | 0.5630 (matched KNN (k=31, α=0.1)) | not met | matched KNN, molecular KNN |
| complementarity · 876222 | 0.5811 | 0.1973 | 0.6931 | 0.5683 (matched KNN (k=31, α=0.1)) | not met | frozen logistic, matched KNN, molecular KNN |
| complementarity · 876223 | 0.5768 | 0.1960 | 0.6931 | 0.5649 (matched KNN (k=31, α=0.1)) | not met | frozen logistic, matched KNN, molecular KNN |
| latent_unobserved · 876231 | 0.7071 | 0.2567 | 0.6931 | 0.6931 (constant 0.5) | not met | depth HGBT, depth KNN |
| latent_unobserved · 876232 | 0.7038 | 0.2552 | 0.6931 | 0.6931 (constant 0.5) | not met | frozen logistic, depth HGBT, depth KNN |
| latent_unobserved · 876233 | 0.7129 | 0.2595 | 0.6931 | 0.6923 (training frequency) | not met | frozen logistic, depth HGBT, depth KNN, matched KNN, molecular KNN |
| location_only · 876211 | 0.6112 | 0.2110 | 0.6931 | 0.6013 (**depth-only HGBT (selected α=0.25)**) | not met | depth HGBT, depth KNN |
| location_only · 876212 | 0.6141 | 0.2115 | 0.6931 | 0.6086 (**depth-only HGBT (selected α=0.25)**) | not met | depth HGBT, depth KNN |
| location_only · 876213 | 0.6130 | 0.2116 | 0.6931 | 0.6099 (**depth-only HGBT (selected α=0.25)**) | not met | depth HGBT |

### By scenario

Control minus candidate: positive means the nonlinear model is better.

**batch_shift_missing**

| comparator | mean log-loss diff [range] | mean Brier diff [range] | cases candidate lower log loss | cases Brier not higher |
|---|---|---|---|---|
| frozen logistic, operational (C=0.1, α=0.1) | +0.0138 [+0.0066, +0.0239] | +0.0052 [+0.0023, +0.0090] | 3/3 | 3/3 |
| **depth-only HGBT (selected α=0.25)** | +0.0225 [+0.0172, +0.0294] | +0.0103 [+0.0079, +0.0130] | 3/3 | 3/3 |
| depth logistic | +0.0219 [+0.0199, +0.0254] | +0.0103 [+0.0093, +0.0114] | 3/3 | 3/3 |
| depth KNN (k=31, α=0.5) | +0.0168 [+0.0118, +0.0254] | +0.0078 [+0.0056, +0.0112] | 3/3 | 3/3 |
| molecular KNN (k=31, α=0.1) | +0.0005 [-0.0060, +0.0053] | +0.0004 [-0.0030, +0.0023] | 2/3 | 2/3 |
| matched KNN (k=31, α=0.1) | +0.0002 [-0.0034, +0.0023] | +0.0005 [-0.0014, +0.0019] | 2/3 | 2/3 |
| copy-1NN depth (α=0.5) | +0.1534 [+0.1234, +0.1781] | +0.0681 [+0.0548, +0.0795] | 3/3 | 3/3 |
| copy-1NN matched (α=0.5) | +0.1021 [+0.0929, +0.1076] | +0.0453 [+0.0415, +0.0480] | 3/3 | 3/3 |
| training frequency | +0.0205 [+0.0141, +0.0282] | +0.0095 [+0.0067, +0.0124] | 3/3 | 3/3 |
| constant 0.5 | +0.0598 [+0.0468, +0.0750] | +0.0289 [+0.0228, +0.0356] | 3/3 | 3/3 |

**complementarity**

| comparator | mean log-loss diff [range] | mean Brier diff [range] | cases candidate lower log loss | cases Brier not higher |
|---|---|---|---|---|
| frozen logistic, operational (C=0.1, α=0.1) | -0.0023 [-0.0075, +0.0035] | -0.0008 [-0.0030, +0.0017] | 1/3 | 1/3 |
| **depth-only HGBT (selected α=0.25)** | +0.0897 [+0.0847, +0.0926] | +0.0410 [+0.0383, +0.0429] | 3/3 | 3/3 |
| depth logistic | +0.0907 [+0.0863, +0.0982] | +0.0415 [+0.0390, +0.0446] | 3/3 | 3/3 |
| depth KNN (k=31, α=0.5) | +0.0901 [+0.0844, +0.0976] | +0.0412 [+0.0383, +0.0444] | 3/3 | 3/3 |
| molecular KNN (k=31, α=0.1) | +0.0001 [-0.0028, +0.0019] | -0.0006 [-0.0014, -0.0001] | 2/3 | 0/3 |
| matched KNN (k=31, α=0.1) | -0.0103 [-0.0128, -0.0061] | -0.0047 [-0.0055, -0.0033] | 0/3 | 0/3 |
| copy-1NN depth (α=0.5) | +0.2210 [+0.2183, +0.2234] | +0.0980 [+0.0965, +0.0998] | 3/3 | 3/3 |
| copy-1NN matched (α=0.5) | +0.1048 [+0.0989, +0.1104] | +0.0459 [+0.0439, +0.0481] | 3/3 | 3/3 |
| training frequency | +0.1022 [+0.0982, +0.1074] | +0.0471 [+0.0458, +0.0491] | 3/3 | 3/3 |
| constant 0.5 | +0.1175 [+0.1120, +0.1240] | +0.0547 [+0.0527, +0.0574] | 3/3 | 3/3 |

**latent_unobserved**

| comparator | mean log-loss diff [range] | mean Brier diff [range] | cases candidate lower log loss | cases Brier not higher |
|---|---|---|---|---|
| frozen logistic, operational (C=0.1, α=0.1) | +0.0004 [-0.0018, +0.0037] | +0.0003 [-0.0009, +0.0018] | 1/3 | 1/3 |
| **depth-only HGBT (selected α=0.25)** | -0.0078 [-0.0115, -0.0022] | -0.0038 [-0.0055, -0.0011] | 0/3 | 0/3 |
| depth logistic | -0.0093 [-0.0131, -0.0034] | -0.0044 [-0.0062, -0.0016] | 0/3 | 0/3 |
| depth KNN (k=31, α=0.5) | -0.0111 [-0.0175, -0.0038] | -0.0053 [-0.0084, -0.0018] | 0/3 | 0/3 |
| molecular KNN (k=31, α=0.1) | -0.0023 [-0.0139, +0.0065] | -0.0012 [-0.0067, +0.0031] | 2/3 | 2/3 |
| matched KNN (k=31, α=0.1) | +0.0003 [-0.0098, +0.0077] | +0.0001 [-0.0048, +0.0038] | 2/3 | 2/3 |
| copy-1NN depth (α=0.5) | +0.1414 [+0.1137, +0.1553] | +0.0608 [+0.0483, +0.0672] | 3/3 | 3/3 |
| copy-1NN matched (α=0.5) | +0.1280 [+0.1105, +0.1511] | +0.0547 [+0.0468, +0.0651] | 3/3 | 3/3 |
| training frequency | -0.0139 [-0.0206, -0.0084] | -0.0067 [-0.0099, -0.0041] | 0/3 | 0/3 |
| constant 0.5 | -0.0148 [-0.0198, -0.0106] | -0.0071 [-0.0095, -0.0052] | 0/3 | 0/3 |

**location_only**

| comparator | mean log-loss diff [range] | mean Brier diff [range] | cases candidate lower log loss | cases Brier not higher |
|---|---|---|---|---|
| frozen logistic, operational (C=0.1, α=0.1) | +0.0153 [+0.0112, +0.0198] | +0.0073 [+0.0053, +0.0093] | 3/3 | 3/3 |
| **depth-only HGBT (selected α=0.25)** | -0.0062 [-0.0100, -0.0031] | -0.0025 [-0.0045, -0.0012] | 0/3 | 0/3 |
| depth logistic | +0.0047 [+0.0011, +0.0098] | +0.0028 [+0.0012, +0.0052] | 3/3 | 3/3 |
| depth KNN (k=31, α=0.5) | -0.0005 [-0.0053, +0.0061] | -0.0001 [-0.0026, +0.0028] | 1/3 | 1/3 |
| molecular KNN (k=31, α=0.1) | +0.0456 [+0.0359, +0.0519] | +0.0206 [+0.0160, +0.0237] | 3/3 | 3/3 |
| matched KNN (k=31, α=0.1) | +0.0194 [+0.0118, +0.0271] | +0.0088 [+0.0054, +0.0120] | 3/3 | 3/3 |
| copy-1NN depth (α=0.5) | +0.1025 [+0.0979, +0.1096] | +0.0441 [+0.0417, +0.0472] | 3/3 | 3/3 |
| copy-1NN matched (α=0.5) | +0.1314 [+0.1294, +0.1351] | +0.0567 [+0.0551, +0.0587] | 3/3 | 3/3 |
| training frequency | +0.0337 [+0.0258, +0.0404] | +0.0156 [+0.0115, +0.0186] | 3/3 | 3/3 |
| constant 0.5 | +0.0804 [+0.0790, +0.0819] | +0.0386 [+0.0384, +0.0390] | 3/3 | 3/3 |

**overall (12 cases)**

| comparator | mean log-loss diff [range] | mean Brier diff [range] | cases candidate lower log loss | cases Brier not higher |
|---|---|---|---|---|
| frozen logistic, operational (C=0.1, α=0.1) | +0.0068 [-0.0075, +0.0239] | +0.0030 [-0.0030, +0.0093] | 8/12 | 8/12 |
| **depth-only HGBT (selected α=0.25)** | +0.0246 [-0.0115, +0.0926] | +0.0113 [-0.0055, +0.0429] | 6/12 | 6/12 |
| depth logistic | +0.0270 [-0.0131, +0.0982] | +0.0125 [-0.0062, +0.0446] | 9/12 | 9/12 |
| depth KNN (k=31, α=0.5) | +0.0238 [-0.0175, +0.0976] | +0.0109 [-0.0084, +0.0444] | 7/12 | 7/12 |
| molecular KNN (k=31, α=0.1) | +0.0110 [-0.0139, +0.0519] | +0.0048 [-0.0067, +0.0237] | 9/12 | 7/12 |
| matched KNN (k=31, α=0.1) | +0.0024 [-0.0128, +0.0271] | +0.0012 [-0.0055, +0.0120] | 7/12 | 7/12 |
| copy-1NN depth (α=0.5) | +0.1546 [+0.0979, +0.2234] | +0.0678 [+0.0417, +0.0998] | 12/12 | 12/12 |
| copy-1NN matched (α=0.5) | +0.1166 [+0.0929, +0.1511] | +0.0507 [+0.0415, +0.0651] | 12/12 | 12/12 |
| training frequency | +0.0357 [-0.0206, +0.1074] | +0.0164 [-0.0099, +0.0491] | 9/12 | 9/12 |
| constant 0.5 | +0.0607 [-0.0198, +0.1240] | +0.0288 [-0.0095, +0.0574] | 9/12 | 9/12 |

### Reading the failures

* **location_only:** projection depends on position alone, and a depth-only model wins.
  * Depth-only HGBT beats the gene + depth model in all 3 seeds, by −0.0062 log loss on average.
  * Extra gene inputs only add noise here.
* **complementarity:** matched-input KNN wins all 3 seeds (−0.0103 log loss on average).
  * The frozen logistic model also wins 2 of 3.
  * The fixed shallow trees capture this gene × position structure less well than local averaging does.
* **latent_unobserved:** the signal comes from a hidden state, so nothing observed predicts it.
  * The expected best observable risk is that of the constant 0.5: log loss log 2 ≈ 0.6931, Brier 0.25.
  * The constant 0.5 or training frequency gives the lowest log loss in every seed.
  * The nonlinear model overfits slightly (0.704–0.713), and depth-only models also beat it.
* **batch_shift_missing:** the gate was met in 2 of 3 seeds.
  * It beats the frozen logistic model and all depth and copy controls in 3 of 3.
  * It roughly ties molecular and matched KNN, which win seed 876243. The result is not consistent.

**Overall:**
* The nonlinear model always beats the copy controls (12/12).
* It beats the frozen logistic model on log loss in 8 of 12 cases, with a mean difference of +0.0068 and a range
  of −0.0075 to +0.0239.
* It meets the full gate in only 2 of 12.

## Decision

* **Kept:** the frozen logistic recipe (`models/paired_projection/logloss_recipe_v1.json`, A1 weights
  `models/paired_projection/a1_logloss_v1.json`) remains the operational model.
* **Archived, not promoted:** the nonlinear artifacts are in `models/paired_projection/nonlinear_development_v1/`,
  marked experimental.
* **Not done:** no further model or parameter search was run.

## Provenance

* **Code and protocol:**
  * Source commit `a76448adb27eb884066b5769c1e9ffada49949ba`, covering the code, tests, protocol and ROOT's
    replay tool.
  * Code `analysis/paired_projection_nonlinear.py`, SHA-256
    `20988e42b544c68df833c7e2967b21e90832287a2982a607db4e90be30db8adc`.
  * Protocol `docs/paired_projection_nonlinear_protocol.json`, SHA-256
    `f7943016ada22a3958e7800db2ca3ecefd8fae8979a3fdd3a713afc0218a28a6`; each run kept an identical snapshot.
* **Selection:** commit `edca36cecc2e8ee0566141435c58ea3eac98ea6d`, record
  `docs/paired_projection_nonlinear_selection.json`, SHA-256
  `ddf6f69814259a57b1a5fc7af6fe2392f39b69eba7146b90f86c6b14e25414d4`.
* **A1 run `nonlinear_a1_dev_v1`:**
  * `selection.json` `d7a8ea2e1409fc1ca09afb688421bfcbaadfeb89004ea2701be175c04f55daa2`
  * `RUN_COMPLETE.json` `ad05b855891f5516385da25bae8fc36c813dc749ceeb9920711f45676784a86c`
  * `result.json` `e7cf0f19c021c6556544199d478349b6bd655938fab3e74fa5295f388ee0daa7`
  * `predictions.csv` `86b173d18e24f78f95cc72ed64f50a5163dc61ffb080d40fa893e1606d4ef816`
* **Simulation run `nonlinear_sim_v1`:**
  * `RUN_COMPLETE.json` `faf4400068474484b12a9c328c47eeb49b7a010b57d4938c638b8abeaec9c177`
  * `metrics.json` `8b0042f00f4677c357849ed68bc62e2688b7b21d8c304b1cea2d9badd730137a`
  * `predictions.csv` `495581c4db21c59759adcde98fe48fa685a529370c6db34f2319905b2b7e416c`
  * The 12 case NPZ files are hashed in `RUN_COMPLETE.json`.
  * The raw CSV and NPZ files stay local and uncommitted.
* **Independent replays (ROOT), with no refits:**
  * A1 `reports/paired_projection_nonlinear_a1_replay.json`
    (`98f1a1840d1f115328085b1f4204831ddad13dc3921d8611871e419a5df54b76`);
  * simulation `reports/paired_projection_nonlinear_simulation_replay.json`
    (`267398ebbdd144a8d39277bdb0168153aedcf104eb892da1b99d371eb6d0a736`).
* **Full metrics:** `reports/paired_projection_nonlinear_metrics.json` contains the complete A1 `result.json`
  and the simulation `metrics.json` in clearly separated `real_development_A1` and `synthetic_simulation`
  sections.

## Sources

* **Data:** BARseq2 A1 pooled data, Sun et al., Nature Neuroscience 2021,
  [DOI 10.1038/s41593-021-00842-4](https://doi.org/10.1038/s41593-021-00842-4). Public deposit:
  [Mendeley Data jnx89bmv4s/2](https://data.mendeley.com/datasets/jnx89bmv4s/2). The extraction and admission
  are recorded in `docs/paired_projection_barseq2_admission.json`.
* **Method:**
  [scikit-learn 1.8 `HistGradientBoostingClassifier`](https://scikit-learn.org/1.8/modules/generated/sklearn.ensemble.HistGradientBoostingClassifier.html)
  (official documentation).

## Limits

* Both A1 mice were exposed during development, so the A1 numbers are development estimates, not a test.
* The simulations are generic sensitivity scenarios with invented gene and target couplings. They are not
  hypothalamic data.
* The probabilities are for captured labels (raw count above 2), not true projection absence.
* There is no claim about full axon geometry, synapses or circuits.
