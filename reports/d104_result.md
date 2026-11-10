# D-104 `run` output

The pre-registered test, run once from commit 8e31fdc (D-104 code unchanged since 7700132; entry 628658c), by `PYTHONPATH=src python3 -I analysis/d104_epiretro.py run`. Run back to back with D-105 (whose primary p was 0.0630), as co-registered.

```
python 3.12.4, numpy 1.26.4, pandas 2.2.2, scipy 1.13.1; seed 20261113
commit 8e31fdceaa25673af973dba3073b65cee6329e32; every input matches its pinned SHA-256

subclasses compared 19 (>= 20 labelled neurons over all experiments, >= 50% of MERFISH cells in HY); neurons 1112, of them 958 in the 14 primary experiments
  PVZ  PVH-SO-PVa Otp Glut
  PVR  AVPV-MEPO-SFO Tbr1 Glut; DMH Hmx2 Gaba; LHA-AHN-PVH Otp Trh Glut; MPN-MPO-PVpo Hmx2 Glut; PVpo-VMPO-MPN Hmx2 Gaba
  MEZ  AHN Onecut3 Gaba; AHN-RCH-LHA Otp Fezf1 Glut; AHN-SBPV-PVHd Pdrm12 Gaba; BST-MPN Six3 Nrgn Gaba; PH-ant-LHA Otp Bsx Glut; PVHd-SBPV Six3 Prox1 Gaba; VMH Nr5a1 Glut
  LZ   DMH-LHA Gsx1 Gaba; DMH-LHA Vgll2 Glut; LHA Barhl2 Glut; PVHd-DMH Lhx6 Gaba; STN-PSTN Pitx2 Glut; ZI Pax6 Gaba
pairs compared 171; not compared (share methylation clusters >= 0.25): 0; covariates 14
measured profiles, read without genes: size 2.026; the sexes agree at r = 0.509 (centred log shares) and r = 0.267 (pair similarities; both sexes pooled, Spearman-Brown: 0.421)
fitted: experiment noise tau 0.45, sigma 2.271
  tau: sigma, synthetic between-sex r   0.00: 2.34, 0.621  0.20: 2.36, 0.597  0.40: 2.26, 0.531  0.60: 2.15, 0.440  0.80: 1.97, 0.345  1.00: 1.64, 0.216  1.25: 0.93, 0.050  1.50: 0.25, -0.017
covariates explain 21.0% of the block-centred similarity; by position strength (mean over 100 synthetic datasets) and the strength matched to the measured share:
  smooth field 100 um    10.4 10.9 12.2 12.1 14.1 15.3 16.4 16.7 16.9 17.3  -> 6.00 (cannot reach it: top of the grid)
  smooth field 250 um    11.0 10.7 13.5 13.2 19.0 20.8 26.9 27.4 28.5 32.6  -> 1.52
  smooth field 500 um    10.6 10.9 14.5 18.9 25.1 32.3 37.1 42.4 41.3 45.9  -> 0.84
  smooth field 1000 um    9.8 11.5 15.8 21.9 26.3 39.1 42.8 50.0 51.9 57.1  -> 0.72
  per-structure field    10.0 10.7 11.6 12.3 14.1 16.0 16.9 18.0 19.0 19.6  -> 6.00 (cannot reach it: top of the grid)
  linear gradient        10.0 11.9 15.4 21.9 25.8 35.1 38.3 41.8 44.9 46.4  -> 0.72

DECISION 1 - a cell type's genes predict where its neurons project, position partialled out
  Holm with D-105 at family-wise 0.05: the smaller one-sided p must be below 0.025, the other then below 0.05
  19 subclasses, 171 pairs: r = +0.223   one-sided p = 0.0810 (largest p over the family; 2000 surrogates each)  -> not significant
    smooth field 100 um    matched strength 6.00 (top of the grid); largest p at strength 6.00: null median +0.014, 95th +0.227, p = 0.0535 (2000 surrogates)
    smooth field 250 um    matched strength 1.52; largest p at strength 4.00: null median -0.006, 95th +0.197, p = 0.0320 (2000 surrogates)
    smooth field 500 um    matched strength 0.84; largest p at strength 1.50: null median -0.009, 95th +0.174, p = 0.0225 (2000 surrogates)
    smooth field 1000 um   matched strength 0.72; largest p at strength 0.75: null median -0.007, 95th +0.190, p = 0.0310 (2000 surrogates)
    per-structure field    matched strength 6.00 (top of the grid); largest p at strength 6.00: null median +0.036, 95th +0.255, p = 0.0810 (2000 surrogates)
    linear gradient        matched strength 0.72; largest p at strength 1.00: null median -0.010, 95th +0.224, p = 0.0510 (2000 surrogates)
  no link, cells moved to gene neighbours: mean r 0%: +0.000, 5%: +0.025, 10%: +0.057, 20%: +0.119, 30%: +0.175
  mislabelling that would give this r with no link: above 30%
DECISION 2 - transfer from the most gene-similar same-zone subclass beats the zone average
  mean profile-cosine gain -0.0292   one-sided p = 0.9280  -> not tested (decision 1 did not pass)

SECONDARY (descriptive; no decision; 500 surrogates per family member, largest p)
  s1 D-099's null (genes relabelled within zone)       r = +0.223  p = 0.0687
  s2 zone blocks only (no position or shared-methylation covariates) r = +0.214  p = 0.0870 (relabelled)
  s3 all 10 targets, single-sex ones sex-adjusted      r = +0.155  p = 0.1717  (19 subclasses, 171 pairs)
  s3u all 10 targets, unadjusted                       r = +0.175  p = 0.1417  (19 subclasses, 171 pairs)
  s3f female experiments only (7 targets)              r = +0.068  p = 0.3453  (19 subclasses, 171 pairs)
  s3m male experiments only (7 targets)                r = +0.332  p = 0.0100  (19 subclasses, 171 pairs)
  s4 centroids from the slab (CCF AP 5.0-7.4 mm)       r = +0.226  p = 0.0739  (19 subclasses, 171 pairs)
  s5 MERFISH panel genes as the fingerprint            r = +0.240  p = 0.0279  (19 subclasses, 171 pairs)
  s6 subclasses with >= 40 labelled neurons            r = +0.500  p = 0.0898  (9 subclasses, 36 pairs)
  s7 the three non-hypothalamic subclasses kept        r = +0.174  p = 0.1098  (22 subclasses, 231 pairs)
  s10 same transmitter (Glut/Gaba) also partialled     r = +0.212  p = 0.0798  (19 subclasses, 171 pairs)
  s13 same-zone pairs only                             r = +0.193  p = 0.3673  (18 subclasses, 46 pairs)
  s14 PH-ant-LHA Otp Bsx Glut placed in LZ             r = +0.249  p = 0.0539  (19 subclasses, 171 pairs)
  s15 cleanly clustered neurons only (>= 10 per subclass) r = +0.171  p = 0.1776  (17 subclasses, 136 pairs)
  s8 without zone PVZ                                  r = +0.218  p = 0.1018  (18 subclasses, 153 pairs)
  s8 without zone PVR                                  r = +0.198  p = 0.2136  (14 subclasses, 91 pairs)
  s8 without zone MEZ                                  r = +0.405  p = 0.0599  (12 subclasses, 66 pairs)
  s8 without zone LZ                                   r = +0.022  p = 0.5509  (13 subclasses, 78 pairs)
  s12 without target PFC                               r = +0.204  p = 0.1118  (19 subclasses, 171 pairs)
  s12 without target MOB                               r = +0.232  p = 0.0798  (19 subclasses, 171 pairs)
  s12 without target STR                               r = +0.250  p = 0.0499  (19 subclasses, 171 pairs)
  s12 without target TH                                r = +0.072  p = 0.3273  (19 subclasses, 171 pairs)
  s12 without target SC                                r = +0.279  p = 0.0279  (19 subclasses, 171 pairs)
  s12 without target P                                 r = +0.225  p = 0.0918  (19 subclasses, 171 pairs)
  s12 without target MY                                r = +0.211  p = 0.0938  (19 subclasses, 171 pairs)
  s8 is read against the ratios in reports/d104_power.md (what leaving a zone out does to a uniform link)
  s9 all pairs, nothing partialled (Mantel)            genes/projections r = +0.227 p = 0.0116; closeness/projections r = +0.259 p = 0.0017; closeness/genes r = +0.205 p = 0.0112
  s11 rule learned outside the hypothalamus            r_unit = +0.095  permutation p = 0.2007  surrogate p = 0.2096
                                                       r_unit = +0.090  permutation p = 0.1842  surrogate p = 0.2156, linear position trend removed  (10 components, ridge 1000)
  s16 each subclass left out (r only)                  r from +0.051 (without DMH-LHA Vgll2 Glut) to +0.281 (without AHN-RCH-LHA Otp Fezf1 Glut); all: +0.26 +0.28 +0.22 +0.20 +0.23 +0.23 +0.24 +0.05 +0.21 +0.24 +0.22 +0.26 +0.22 +0.17 +0.24 +0.25 +0.27 +0.26 +0.13
```
