# D-104 rehearsal outside the hypothalamus

Produced by `python3 -I analysis/d104_rehearsal.py report`. No hypothalamic-source projection was read: Supplementary Table 2 is filtered to non-hypothalamic sources on loading, and the 19 hypothalamic test subclasses are removed from every source (their gene means are used, genes only).

Data: 27 non-hypothalamic sources in 7 divisions, 187 units (subclass in a source with >= 20 labelled neurons).

## 1. Transfer: the ridge rule, leave one division out

Ridge from gene principal components (2,000 most variable genes) to the double-centred log-share profile, sample weights (n / mean n)^0.5. Each division is held out in turn, with its subclasses removed from training. r_unit: weighted correlation of unit-row-normalised prediction and measured profile; p from 2000 within-source permutations. Selection metric: mean r_unit over divisions (25 s).

| PCs | lambda | amygdala-pallidum | hindbrain | hippocampal | isocortex | midbrain | olfactory | thalamus | mean | subcortical | positive |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 100 | +0.325 (p 0.0645) | +0.009 (p 0.4838) | +0.199 (p 0.0750) | +0.263 (p 0.0010) | +0.118 (p 0.2409) | +0.112 (p 0.2364) | -0.098 (p 0.5142) | +0.133 | +0.089 | 6/7 |
| 10 | 1000 | +0.335 (p 0.0870) | -0.128 (p 0.8131) | +0.223 (p 0.0630) | +0.360 (p 0.0005) | +0.183 (p 0.1269) | +0.225 (p 0.1314) | -0.073 (p 0.4483) | +0.161 | +0.080 | 5/7 |
| 10 | 10000 | +0.272 (p 0.2374) | -0.309 (p 0.9815) | +0.222 (p 0.1404) | +0.176 (p 0.0240) | +0.184 (p 0.1359) | +0.290 (p 0.0895) | -0.092 (p 0.5167) | +0.106 | +0.014 | 5/7 |
| 20 | 100 | +0.234 (p 0.0865) | -0.096 (p 0.7591) | +0.281 (p 0.0555) | +0.054 (p 0.2369) | +0.002 (p 0.5622) | +0.442 (p 0.0085) | -0.100 (p 0.5992) | +0.117 | +0.010 | 5/7 |
| 20 | 1000 | +0.108 (p 0.2844) | -0.150 (p 0.8586) | +0.327 (p 0.0400) | +0.207 (p 0.0020) | +0.084 (p 0.3093) | +0.319 (p 0.0605) | +0.152 (p 0.1734) | +0.149 | +0.048 | 6/7 |
| 20 | 10000 | -0.009 (p 0.4598) | -0.108 (p 0.7206) | +0.263 (p 0.1079) | +0.120 (p 0.0335) | +0.038 (p 0.4293) | +0.310 (p 0.0760) | +0.022 (p 0.5532) | +0.091 | -0.015 | 5/7 |
| 40 | 100 | +0.304 (p 0.0425) | -0.062 (p 0.6572) | +0.248 (p 0.0840) | +0.073 (p 0.0875) | -0.061 (p 0.6832) | +0.384 (p 0.0345) | -0.222 (p 0.8631) | +0.095 | -0.010 | 4/7 |
| 40 | 1000 | +0.172 (p 0.2179) | -0.022 (p 0.5552) | +0.208 (p 0.1289) | +0.199 (p 0.0015) | +0.153 (p 0.1139) | +0.354 (p 0.0515) | +0.000 (p 0.5002) | +0.152 | +0.076 | 6/7 |
| 40 | 10000 | +0.070 (p 0.3888) | -0.018 (p 0.5607) | +0.180 (p 0.2104) | +0.175 (p 0.0065) | +0.154 (p 0.1539) | +0.355 (p 0.0465) | +0.069 (p 0.4318) | +0.141 | +0.069 | 6/7 |

Best mean r_unit: PCs 10, lambda 1000 (+0.161).

## 2. Similarity test (design A) with location partialled out

Pairs of units within a source; gene similarity from d104_epiretro.fingerprints over the division's subclasses; profile similarity from d104_epiretro.profiles on targets labelled in both sexes. Units need >= 50% of their MERFISH cells in the division. 14 covariates (centroid distance, three Gaussian proximities, three intermingling scales, shared methylation clusters, six coordinate products), block = source, weights 1/(1/n_i + 1/n_j). Null: gene labels permuted within source, 2000 permutations (110 s).

| division | units | pairs | r raw | z raw | p raw | r adjusted | z adjusted | p adjusted | position explains proj / genes |
|---|---|---|---|---|---|---|---|---|---|
| amygdala-pallidum | 16 | 57 | +0.266 | +1.52 | 0.0830 | -0.113 | -0.65 | 0.7506 | 66% / 71% |
| hindbrain | 23 | 123 | +0.452 | +3.83 | 0.0010 | +0.363 | +3.02 | 0.0060 | 24% / 46% |
| hippocampal | 10 | 24 | -0.233 | -0.90 | 0.7901 | -0.495 | -1.19 | 0.8466 | 41% / 76% |
| isocortex | 58 | 185 | +0.620 | +5.94 | 0.0005 | +0.436 | +4.33 | 0.0005 | 31% / 45% |
| midbrain | 21 | 66 | +0.170 | +1.26 | 0.1139 | +0.132 | +0.92 | 0.1759 | 27% / 31% |
| thalamus | 12 | 19 | +0.458 | +1.86 | 0.0360 | -0.482 | -0.66 | 0.6542 | 96% / 97% |

Mean z over all divisions: +2.25 raw, +0.97 adjusted; four subcortical: +2.12 raw, +0.66 adjusted.

Sources per division:

- amygdala-pallidum: AMY 9 units, 7 targets (ENT,HPF,MOB,PFC,STR,TH,VTA), position from parcel:9; PAL 7 units, 5 targets (P,SC,STR,TH,VTA), position from parcel:7
- hindbrain: MY 13 units, 5 targets (CBN,CBX,HY,SC,TH), position from parcel:11/division:2; P 10 units, 5 targets (CBN,CBX,HY,MY,TH), position from parcel:10
- hippocampal: CAp 3 units, 4 targets (ENT,HY,MOB,STR), position from parcel:3; DGp skipped (2 units, 2 both-sex targets); ENT 7 units, 8 targets (ACA,HPF,PFC,PIR,PTLp,RSP,STR,TH), position from parcel:5/division:2
- isocortex: ACA 7 units, 11 targets (AI,AUDp,MOp,P,PTLp,RSP,SC,STR,TH,VISp,VTA), position from parcel:7; AI 8 units, 8 targets (ACA,AUDp,PFC,RSP,SC,SSp,STR,VISp), position from parcel:8; AUD 7 units, 9 targets (ACA,MOp,P,PTLp,RSP,SC,STR,TH,VISp), position from parcel:7; MOp 6 units, 9 targets (ACA,AI,AUDp,P,SC,SSp,STR,TH,VTA), position from parcel:6; PTLp 6 units, 8 targets (ACA,AUDp,MOp,P,RSP,SC,TH,VISp), position from parcel:6; RSP 8 units, 8 targets (ACA,AUDp,MOp,P,PTLp,SC,TH,VISp), position from parcel:8; SSp 7 units, 9 targets (ACA,AUDp,MOp,MY,P,PTLp,RSP,SC,TH), position from parcel:7; VIS 9 units, 8 targets (ACA,AUDp,P,PTLp,RSP,SC,STR,TH), position from parcel:9
- midbrain: IC skipped (0 units, 0 both-sex targets); MRN 6 units, 3 targets (MOB,P,PAL), position from parcel:6; PAG skipped (7 units, 1 both-sex targets); SC 6 units, 3 targets (MY,P,TH), position from parcel:5/division:1; VTA 9 units, 11 targets (ACA,HY,MOB,P,PAL,PTLp,RSP,SC,SSp,STR,TH), position from parcel:9
- olfactory: nothing to compare
- thalamus: THl 3 units, 7 targets (ACA,AUDp,MOp,P,PTLp,RSP,SC), position from parcel:3; THm 5 units, 6 targets (AUDp,MOp,PFC,RSP,SSp,VISp), position from parcel:5; THp 4 units, 7 targets (ACA,AUDp,P,PTLp,SSp,VISp,VTA), position from parcel:4

## 3. The frozen outside rule

`data/interim/d104_outside_rule.npz`, loaded by `d104_epiretro.outside_rule()`: 10 PCs of 2000 genes, ridge 1000, weight power 0.5, fitted on all 27 non-hypothalamic sources (187 units, 94 subclasses). 24 targets: ACA, AI, AMY, AUDp, CBN, CBX, ENT, HPF, HY, MOB, MOp, MY, P, PAL, PFC, PIR, PTLp, RSP, SC, SSp, STR, TH, VISp, VTA. Self-check on ACA, DGp, MY, PTLp: the stored rule reproduces the fitted predictions (max difference 6.0e-15).

Total run time 135 s.
