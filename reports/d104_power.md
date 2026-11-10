# D-104 `power` output

Size and power of the D-104 test on synthetic projection data, by `PYTHONPATH=src python3 analysis/d104_epiretro.py power`. The real gene fingerprints, MERFISH positions and experiment totals are used; projections are synthetic. The measured projection profiles were read without genes, to set the synthetic profiles' size, experiment noise (from how far the two sexes agree) and position strength. "primary" is the pass rate at one-sided 0.05; "p<0.025" the rate at the first Holm threshold shared with D-105.

```
python 3.12.4, numpy 1.26.4, pandas 2.2.2, scipy 1.13.1; seed 20261113
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
matched 500 um strength 0.84; re-matched beside a gene term: 0.50 -> 0.93, 0.75 -> 1.00, 1.00 -> 1.21, 1.50 -> 1.26, 2.00 -> 1.66; sigma at tau 0: 2.370, at tau 0.90: 1.801
rates over 1000 synthetic datasets, one-sided 0.05; null: 2000 surrogates per family member at each grid strength, p interpolated between grid strengths; each member takes its largest p at the matched strength or any grid strength above. primary = surrogate null (largest p over the family of 6); relabelled = D-099's null; transfer = primary AND transfer gain against the surrogates
  model                                            mean r  primary  p<0.025 relabelled  transfer
  noise only                                       -0.001    0.007    0.004      0.042     0.000
  bumps 250 um                                     +0.003    0.008    0.005      0.049     0.001
  bumps 500 um                                     -0.006    0.009    0.002      0.038     0.000
  bumps 1000 um                                    -0.010    0.014    0.003      0.056     0.000
  bumps 250 um strong                              +0.009    0.012    0.003      0.059     0.000
  bumps 1000 um strong                             -0.029    0.003    0.001      0.028     0.000
  bumps 500 um at the matched strength             +0.001    0.014    0.005      0.062     0.001
  smooth field 100 um strong                       +0.020    0.026    0.014      0.069     0.003
  smooth field 250 um strong                       +0.008    0.025    0.011      0.069     0.004
  smooth field 1000 um strong                      -0.031    0.005    0.002      0.031     0.000
  per-structure field strong                       +0.048    0.047    0.020      0.146     0.008
  outside the family: smooth 150 um strong         +0.011    0.019    0.011      0.066     0.002
  outside the family: smooth 2000 um               -0.013    0.019    0.008      0.055     0.002
  outside the family: bumps 100 um strong          +0.021    0.023    0.012      0.080     0.001
  dorsal-ventral gradient 1.0                      +0.008    0.024    0.007      0.102     0.003
  medial-lateral gradient 1.0                      -0.002    0.037    0.016      0.088     0.000
  3D linear gradient 1.0                           -0.010    0.019    0.005      0.074     0.001
  dorsal-ventral 0.75 + smooth 500 um 0.74         -0.002    0.017    0.006      0.055     0.002
  matched, experiment noise 0                      -0.002    0.009    0.004      0.048     0.000
  matched, experiment noise doubled                +0.002    0.009    0.003      0.053     0.001
  matched + 5% moved along shared clusters         +0.009    0.021    0.009      0.064     0.001
  matched + 10% moved along shared clusters        +0.009    0.012    0.007      0.069     0.001
  matched + 20% moved along shared clusters        +0.017    0.022    0.009      0.074     0.002
  matched + 5% moved to gene neighbours            +0.030    0.027    0.013      0.093     0.004
  matched + 10% moved to gene neighbours           +0.061    0.063    0.024      0.167     0.005
  matched + 20% moved to gene neighbours           +0.115    0.137    0.079      0.289     0.026
  matched + 30% moved to gene neighbours           +0.165    0.255    0.148      0.435     0.085
  matched + 20% moved to gene neighbours, broad    +0.070    0.065    0.029      0.165     0.006
  genes 0.50 + re-matched position                 +0.036    0.034    0.012      0.121     0.001
  genes 0.75 + re-matched position                 +0.075    0.079    0.037      0.192     0.004
  genes 1.00 + re-matched position                 +0.095    0.114    0.061      0.232     0.015
  genes 1.50 + re-matched position                 +0.147    0.202    0.117      0.398     0.024
  genes 2.00 + re-matched position                 +0.167    0.261    0.154      0.429     0.044
same-zone pairs only (46 pairs; covariates explain 49.7%; matched strength 1.10):
  same zone: noise only                            -0.000    0.012    0.007      0.053     0.002
  same zone: matched position                      +0.013    0.012    0.008      0.055     0.001
  same zone: genes 1.00 + matched position         +0.120    0.055    0.031      0.143     0.010
each zone left out, a uniform link (genes 1.5): mean r +0.154 with all zones; ratio without PVZ 1.00, PVR 1.18, MEZ 1.06, LZ 0.61
no link, cells moved to gene neighbours (1000 datasets each): mean r 0%: +0.000, 5%: +0.025, 10%: +0.057, 20%: +0.119, 30%: +0.175
```
