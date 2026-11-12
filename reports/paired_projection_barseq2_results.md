# BARseq2 gene-and-depth transfer: real-neuron result

2026-11-12. **Positive result against copying, with a trade-off against the strongest location control. The full two-score comparison is not a win.**

Claude executed the frozen model once on 598 paired auditory-cortex neurons from two separate mice, XC95 and XC97. Training and evaluation were separated by whole animal in both directions. Each cell has 23 measured marker genes, measured cortical depth and 11 captured-target labels. All ten methods evaluated exactly the same eligible cells and targets.

## Main result

Both metrics are lower-is-better probability losses. Values average cells within animal, animals equally within target, then the 11 targets equally. They are not percentages of correct circuits.

| Method | Brier | Log loss |
|---|---:|---:|
| Gene + depth logistic | 0.100289 | 0.356117 |
| Depth-only logistic | 0.108138 | 0.351664 |
| Gene + depth, 15-neighbour averaging | 0.104577 | 0.453361 |
| Depth, 15-neighbour averaging | 0.114605 | 0.505398 |
| Genes, 15-neighbour averaging | 0.106942 | 0.471079 |
| Copy one gene + depth nearest neuron | 0.185896 | 2.568248 |
| Copy one depth nearest neuron | 0.199267 | 2.752982 |
| Technical count-proxy logistic | 0.116217 | 0.373950 |
| Training capture frequency | 0.117176 | 0.373284 |
| Constant 0.5 | 0.250000 | 0.693147 |

The candidate has the lowest Brier error of all controls. Relative to depth-only logistic it improves Brier by **7.26%**, but log loss is **1.27% worse**. Relative to the strongest copying/averaging comparator, matched-input 15-neighbour averaging, it improves Brier by **4.10%** and log loss by **21.45%**. These are descriptive differences on this dataset, without significance or population-level claims.

One-neuron copying produces hard probabilities of 0 or 1. Wrong confident calls are heavily penalized by log loss (clipped at 1e-6); beating those controls alone would be a weak bar. The matched-input averaging control is the more informative comparison. Lower Brier or log loss alone does not establish calibrated probabilities.

## Each held-out animal

| Evaluation mouse | Gene + depth Brier / log loss | Depth-only logistic | Matched-input averaging |
|---|---:|---:|---:|
| XC95 | 0.108548 / 0.407941 | 0.115300 / 0.379526 | 0.113324 / 0.551801 |
| XC97 | 0.092029 / 0.304294 | 0.100976 / 0.323802 | 0.095830 / 0.354921 |

Brier improves against every location/copying control in both transfer directions. Both scores improve against matched-input averaging in each mouse. Against depth-only logistic, log loss improves only when XC97 is evaluated; it worsens when XC95 is evaluated.

## Destination-specific behavior

Against depth-only logistic, four targets improve on both scores: VisIp, AudC, Thal and Tect. The other seven worsen on both scores. The aggregate improvement therefore does not apply uniformly across destinations. The largest log-loss penalties arise for Rstr, SS and Cstr; the strongest benefits occur for AudC, Thal and Tect. No post-result exclusions, thresholds, model tuning or recalibration were used.

| Target | Candidate minus depth Brier | Candidate minus depth log loss |
|---|---:|---:|
| OFC | +0.001697 | +0.020077 |
| M | +0.000120 | +0.002933 |
| Rstr | +0.017988 | +0.161527 |
| SS | +0.013599 | +0.088783 |
| Cstr | +0.019755 | +0.078940 |
| Amyg | +0.000841 | +0.014632 |
| VisIp | -0.008133 | -0.015945 |
| VisC | +0.000362 | +0.005999 |
| AudC | -0.059304 | -0.101419 |
| Thal | -0.041730 | -0.100451 |
| Tect | -0.031537 | -0.106090 |

## Verification and provenance

- Freeze commit: [493c192](https://github.com/andresloak100/SYN4RCH/commit/493c1923b623de367a51720eec777bd3b20954a0), committed before the biological run.
- Run: `barseq2_a1_598_logo_v1`, 2026-11-12 11:12:15 EST to 11:12:16 EST. A durable started marker and exclusive outputs prevent a silent repeat.
- Software verification: 121 focused tests pass; ROOT independently reproduced the simulated software checks.
- ROOT replay: 65,780 keyed saved predictions, exact source labels/query IDs/animal folds, and 13,156 nearest-training-donor copying checks. Maximum score difference: 3.77e-15. No predictor was refitted by the replay.
- `result.json`: SHA-256 `9f90168fb62b2acad6e9b403fea4af10f4ac42d83f68dd59ababdf91b6cc6f23`.
- `predictions.csv`: SHA-256 `e4abf8d17bb07ba27ed2d2c67064bf471edb8d204eb013430fdef19b617cadc5`.
- `RUN_COMPLETE.json`: SHA-256 `ff304e74a0df65b444e27afb245950ac7373931147ab2743e6b055de8bc1aac3`.
- `ROOT_REPLAY.json`: SHA-256 `fc7700441e23b976fa6d120aeb0b34ac35c127f5896f56dfa366612b25613874`.
- Local outputs: `data/interim/paired_projection_real_benchmark/barseq2_a1_598_logo_v1/`. Raw data and full predictions remain local; the source/admission and protocol are committed.

## Interpretation and limits

The result supports further work on measured molecular information for destination prediction, and demonstrates better probability predictions than anatomical or matched-feature copying in these two mice. It does not establish a joint improvement over every location control, a new biological discovery or transfer to a population of animals. No confirmatory p-value is reported.

The cohort is selected for marker QC and detectable projection somewhere, with several destinations represented by very few captured cells. The measured anatomical control is one laminar axis, so the comparison cannot establish an advantage beyond full XYZ anatomy. Captured labels and nondetections are assay outcomes, not proof of anatomical absence, synapses or complete axon geometry. These are auditory-cortex neurons, not hypothalamic neurons.

The fixed run is closed. Any later method must be labeled development on these exposed animals and evaluated independently before a broader claim or viewer promotion.

Protocol: [paired_projection_real_protocol.json](../docs/paired_projection_real_protocol.json). Admission: [paired_projection_barseq2_admission.json](../docs/paired_projection_barseq2_admission.json). Source: [Sun et al.](https://doi.org/10.1038/s41593-021-00842-4), [CC BY 4.0 Mendeley version 2](https://data.mendeley.com/datasets/jnx89bmv4s/2).
