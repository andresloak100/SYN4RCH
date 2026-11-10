# D-105 `run` output

Produced once by `python3 -I analysis/d105_benchmark.py run` at commit 8e31fdc (hypothalamic rows opened 2026-11-10T14:44:18Z). Verbatim below.

```
commit 8e31fdceaa25673af973dba3073b65cee6329e32
  script analysis/d105_benchmark.py  sha256 464df9b6cbd13a290a6c7638cefa97d4cd8122656ff81e6d937af816bd3cd7da  git blob ed419efbabe36f4622847a3520c098b8f3213e7f (HEAD blob ed419efbabe36f4622847a3520c098b8f3213e7f)
  sha256 4855044e38ea6aa782a5a6cc9379d7da5e8bfd83c42dacccdd6b2bf42de7b1cf  st2
  sha256 b21ca985652fb25f9608f99005139a40757133a76fbe845ae5b175c5c26a447b  stats
  sha256 653606461801393cd0a7708a084e3f9c989ffb2ac696ac30629f6ef6180e76c5  terms
  sha256 ec5a211c23fab5cce8a06f01bb0a68af087fb979b647ab42d635201f98de2f69  genes_10x
  sha256 1e3ae23cc3f8d3d5839cc8a222798ff812fc66f4bd3c0fc8fca0c0279329874e  merfish
  sha256 413cfeb2a0628a0da385aa6362b9cf5a07a8c01c0e2f8a21b040bfc6dfe09c46  models
  sha256 77d4ddd71cccbb34736996e41f12741e02ae06b86c37a76864a3ff4ea4b53a51  predictions
  sha256 51671d38612ccbc084231e93218983db0d9b4108f5e686b37eda6943d49d80e8  t1_seed20
  sha256 c735c5f239d0266d701af27bb7b144f37aab99a69bf3d794957734f8e26999ca  t1_seed21
  sha256 f9533fc28297709a4b42c270d8d61e8254f00004a194a6e67369d16d865f67fc  t1_seed22
  sha256 2fff8cfb6328db9fc50cd31e25131dd6777636ceec9a8e27dd6daf6f0932a11e  t1_seed23
  sha256 a9975f45b6d8dfa29f3ed96cae67ed8685fe6d28e6172fa88896878840e768f9  t1_seed24
  sha256 31245b5b8fdee05cb18bb4bf6e22ccce89618607842864bef87818b2ef9d5270  t1_seed25
hypothalamic rows opened 2026-11-10T14:44:18Z (recorded in data/interim/d105_run_opened.json)
# D-105 result: transfer of gene-to-projection models to the hypothalamus

## Layout

19 units (HY19), 14 experiments: Jf6H7J8K9L MOB female (27 of the 19); Jm6H7J8K9L MOB male (45 of the 19); Pf6H7J8K9L STR female (45 of the 19); Pm6H7J8K9L STR male (43 of the 19); Rf6H7J8K9L TH female (78 of the 19); Rm6H7J8K9L TH male (97 of the 19); Tf6H7J8K9L SC female (152 of the 19); Tm6H7J8K9L SC male (94 of the 19); Vf6H7J8K9L P female (92 of the 19); Vm6H7J8K9L P male (43 of the 19); Wf6H7J8K9L MY female (41 of the 19); Wm6H7J8K9L MY male (81 of the 19); Zf6H7J8K9L PFC female (37 of the 19); Zm6H7J8K9L PFC male (83 of the 19).
Neurons of the 19 units over the 14 experiments: 958 (per unit 39 34 26 28 15 87 22 21 29 26 55 17 119 17 164 46 25 92 96).
Spatial strata (fixed in freeze): sizes [4, 3, 4, 4, 4].

## Primary decision (D1)

r_unit(D) = +0.2041; within-strata permutation p = 0.0630 (10000 permutations, seed 20261110, one-sided; ties within 1e-12 count as at least as large).
r_unit(A) = -0.2656; r_D > max(r_A, 0): yes.
Holm with D-104 Decision 1 (family-wise 0.05): pass at the Holm step means p < 0.025 if D-104's p is larger or not yet known; p < 0.05 if D-104's p < 0.025.
  threshold 0.025 (D-104 p >= 0.025 or unknown): not pass
  threshold 0.05  (D-104 p <  0.025):            not pass
Noise ceiling (between-sex reliability, no genes; female-only vs male-only profile in the r_unit form with the test's unit weights): r +0.5033, Spearman-Brown +0.6696. This is a reliability; the largest attainable r_D is about its square root (0.82). Plain Pearson of the two double-centred profiles (descriptive): +0.5093.
95% interval of r_D (units resampled, 2000 draws, seed 20261111; r_D +/- t(0.975, 18) x sqrt(19/18) x bootstrap SD): [-0.0636, +0.4718].
Readings if not PASS: INSUFFICIENT EVIDENCE if the ceiling < 0.6, the interval contains 0.38, or either is undefined; NEGATIVE if the ceiling >= 0.8 and the interval's upper bound < 0.25; otherwise NOT SHOWN.

READING: INSUFFICIENT EVIDENCE

## Secondary (no decision)

| model | r_unit | p (within strata) | max-T family p | top-target hits |
|---|---|---|---|---|
| D | +0.2041 | 0.0630 | 0.1865 | 5/19 = 0.26 |
| O2 | +0.2108 | 0.0528 | 0.1597 | 7/19 = 0.37 |
| K | +0.1356 | 0.2279 | 0.5050 | 4/19 = 0.21 |
| B | +0.0902 | 0.6540 | 0.7018 | 3/19 = 0.16 |
| A | -0.2656 | 0.4773 (a) | - | 0/19 = 0.00 |
| T1 | +0.0766 | 0.2276 | - | 4/19 = 0.21 |

(a) A's permutation p is not calibrated when the true profile is a smooth function of position (A is itself smooth in position; its size under smooth truths is in reports/d105_power.md); only r_A enters the decision.
max-T over {D, O2, K, B}: family p of the best model 0.1597 (max observed +0.2108). Max-T can exceed its nominal size under smooth position-only truths (driven by K and B; sizes in reports/d105_power.md).
Top target = the target with the largest value of the double-centred profile, for prediction and measurement.

Semi-partial r of D after removing A's prediction: +0.2335, p = 0.0239 (same permutations, D's rows only).
Top-target accuracy of D: 5/19 = 0.26 vs chance 1/7 = 0.14; one-sided binomial p 0.1236.
Descriptive: units whose measured top target is the same on raw log shares as on the double-centred profile: 15/19.
T1 single seeds r_unit: +0.080 +0.181 +0.163 +0.055 +0.055 -0.001.

## Per-unit cosines (prediction vs measured, double-centred rows)

| unit | n | stratum | measured top | D | O2 | K | B | A | T1 |
|---|---|---|---|---|---|---|---|---|---|
| PVH-SO-PVa Otp Glut | 39 | 0 | MY | +0.90 | +0.88 | +0.49 | +0.19 | +0.42 | +0.48 |
| AVPV-MEPO-SFO Tbr1 Glut | 34 | 1 | MY | +0.64 | +0.71 | +0.43 | -0.25 | +0.24 | +0.35 |
| DMH Hmx2 Gaba | 26 | 2 | PFC | +0.25 | +0.21 | -0.09 | +0.03 | -0.21 | +0.09 |
| LHA-AHN-PVH Otp Trh Glut | 28 | 3 | MOB | +0.70 | +0.70 | +0.61 | +0.81 | -0.82 | +0.58 |
| MPN-MPO-PVpo Hmx2 Glut | 15 | 1 | MOB | +0.75 | +0.79 | +0.83 | +0.77 | -0.97 | +0.80 |
| PVpo-VMPO-MPN Hmx2 Gaba | 87 | 1 | MOB | -0.71 | -0.67 | -0.62 | -0.69 | -0.66 | -0.54 |
| AHN Onecut3 Gaba | 22 | 3 | SC | +0.13 | +0.07 | -0.11 | -0.33 | -0.35 | -0.14 |
| AHN-RCH-LHA Otp Fezf1 Glut | 21 | 3 | MOB | +0.47 | +0.44 | +0.22 | +0.32 | -0.27 | +0.37 |
| AHN-SBPV-PVHd Pdrm12 Gaba | 29 | 3 | PFC | -0.69 | -0.69 | -0.62 | -0.37 | -0.33 | -0.68 |
| BST-MPN Six3 Nrgn Gaba | 26 | 0 | MY | -0.68 | -0.67 | -0.57 | -0.33 | +0.34 | -0.52 |
| PH-ant-LHA Otp Bsx Glut | 55 | 4 | P | -0.51 | -0.56 | +0.07 | -0.05 | -0.88 | -0.42 |
| PVHd-SBPV Six3 Prox1 Gaba | 17 | 0 | P | +0.32 | +0.27 | +0.47 | +0.55 | -0.32 | -0.38 |
| VMH Nr5a1 Glut | 119 | 2 | SC | +0.55 | +0.56 | +0.51 | +0.47 | +0.28 | +0.31 |
| DMH-LHA Gsx1 Gaba | 17 | 2 | PFC | +0.59 | +0.74 | -0.30 | -0.07 | +0.70 | +0.88 |
| DMH-LHA Vgll2 Glut | 164 | 2 | MY | +0.20 | +0.19 | +0.26 | +0.20 | -0.04 | -0.01 |
| LHA Barhl2 Glut | 46 | 0 | TH | +0.13 | +0.11 | +0.38 | +0.02 | -0.34 | -0.07 |
| PVHd-DMH Lhx6 Gaba | 25 | 4 | SC | +0.31 | +0.29 | +0.57 | +0.73 | -0.70 | -0.21 |
| STN-PSTN Pitx2 Glut | 92 | 4 | STR | +0.42 | +0.48 | -0.26 | +0.08 | -0.59 | +0.56 |
| ZI Pax6 Gaba | 96 | 4 | SC | +0.53 | +0.52 | +0.36 | +0.12 | -0.64 | +0.25 |

exit 0
```
