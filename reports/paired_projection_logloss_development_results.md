# Log-loss development improvement and frozen model

2026-11-12. **The bounded development round improved both losses and the model is frozen. No new real-animal or hypothalamic validation has run.**

The selected gene-plus-depth setting is **C=0.1, shrinkage alpha=0.1**. On the exposed A1 development animals, log loss decreases from **0.356117 to 0.319043** (10.41% lower); Brier decreases from **0.100289 to 0.095234** (5.04% lower). These are model-selection/development numbers, with optimism from choosing the settings on these same predictions.

## What changed

The model uses stronger L2 regularization and moves predictions 10% toward the target capture frequency measured in the fitting animals. No class-frequency reweighting, target-specific hyperparameters, feature selection or new architecture was added. This improves the probability losses on this dataset; it is not proof of calibrated probabilities or a new biological rule.

The fixed grid contained 82 configurations in total. Candidate, location-only and technical logistic families each received the same 3-by-4 C/shrinkage grid. Three neighbour-averaging families each received a 3-by-4 neighbour/shrinkage grid; two copying families received four shrinkage choices, retaining raw copying anchors. One global setting per family minimized log loss while preserving its own baseline Brier. Every configuration and the unconstrained log-loss minimum are retained in the local results. No grid was expanded after seeing the scores.

## Selected controls

Both metrics are lower-is-better, with animals equal within target and targets equal. Every method used identical cells, targets and animal folds. Control settings below were selected on the same exposed development data.

| Family | Selected setting | Brier | Log loss |
|---|---|---:|---:|
| Genes + depth logistic | C=0.1|alpha=0.1 | 0.095234 | 0.319043 |
| Depth-only logistic | C=0.1|alpha=0.0 | 0.108020 | 0.350867 |
| Depth nearest-neighbour averaging | k=31|alpha=0.5 | 0.109125 | 0.353500 |
| Matched gene + depth averaging | k=31|alpha=0.1 | 0.103470 | 0.349359 |
| Molecular averaging | k=31|alpha=0.1 | 0.106211 | 0.357252 |
| Copy one depth-nearest neuron | alpha=0.5 | 0.131697 | 0.416290 |
| Copy one matched-input neuron | alpha=0.5 | 0.124156 | 0.397459 |
| Technical proxies | C=1.0|alpha=0.0 | 0.116217 | 0.373950 |
| Training capture frequency | reference | 0.117176 | 0.373284 |
| Constant 0.5 | reference | 0.250000 | 0.693147 |

The selected candidate beats every selected control on both aggregate losses. This remains a development result; model and control selection used the same two animals. Some unconstrained control minima were rejected by the shared Brier guard and remain explicitly reported.

| Held-out development animal | New gene + depth Brier / log loss | Selected depth logistic | Selected matched averaging |
|---|---:|---:|---:|
| XC95 | 0.101622 / 0.346122 | 0.114994 / 0.377299 | 0.111945 / 0.389506 |
| XC97 | 0.088846 / 0.291964 | 0.101045 / 0.324436 | 0.094995 / 0.309212 |

Against selected depth-only logistic, the candidate improves Brier on 7/11 targets, log loss on 7/11, and both on 7/11. It does not improve every destination.

## Audit and freeze

- Development grid/code were committed before fitting: [74bd311](https://github.com/andresloak100/SYN4RCH/commit/74bd3111ca7613bfa32a4d90c5e8acf9389caf67).
- One bounded real-data development run began at 2026-11-12 12:17:50 EST. The older one-shot run remains unchanged.
- ROOT independently ran 168 focused tests and reproduced the synthetic software checks.
- Independent saved-prediction replay verified **539,396 rows and 82 configurations**, identities, labels, folds, losses, selections, comparisons, fitting-only shrinkage and the deployment artifact. Maximum score difference was zero. This audit performed no predictor fits or development reruns.
- Fitted A1 model: [a1_logloss_v1.json](../models/paired_projection/a1_logloss_v1.json). It contains actual coefficients, scalers and training capture frequencies, refitted on both training mice.
- Frozen generic recipe and all selected control settings: [logloss_recipe_v1.json](../models/paired_projection/logloss_recipe_v1.json).
- Freeze record and hashes: [paired_projection_logloss_freeze.json](../docs/paired_projection_logloss_freeze.json).
- Original deployment artifact SHA-256: `7cdb560801787680d6458c5a798386b3dfbeb32f9162e7f6ad2dc2caeb39c2f1`.
- Development result SHA-256: `2df329387bec6ea6a1a1aa4c0184bfed236d8c139e657c310151defae4822647`.
- Independent replay SHA-256: `5cfa13b9cfe66d00175924b438ce8e1686c887f8d15e9070830b6fec8c78aa25`.

The pure prediction interface accepts measured feature matrices and no target labels. It refuses incompatible gene panels, negative counts, missing/nonfinite values and invalid dimensions. An incompatible hypothalamic assay requires fitting the frozen recipe on its own separately admitted training animals; A1 coefficients cannot be presented as a hypothalamic model.

## Existing data and simulation route

The user excludes lab collaboration. The strategy is therefore existing public/released data and clearly labeled simulation, rather than generating new wet-lab measurements. The checked sources do not yet admit a new hypothalamic animal test:

- **Coester-2026:** the documentation describes 98 endogenous genes plus two viral probes in spatial cells. Public deposition is deferred; the checked repository has no paired tables or animal/target replication metadata. Its public Shiny browser can display viral probes, but the inspected interface exports plots rather than exact paired count tables with animal identities. [Approved data request, already posted once](https://github.com/perslab/Coester-2026/issues/1).
- **Li/Resch PVH:** there are genuinely traced MERFISH animals and a measured mCherry probe. Slide3_male traces spinal cord and Slide6_male traces PB. One actual spatial animal per target and pooled sequencing do not establish the required repeated-animal test. They must not be described as entirely untraced or assigned synthetic animal identities.
- **Rossi LHA:** [the earlier data request](https://github.com/stuberlab/Rossi-et-al-2021/issues/1) is posted; per-cell HCR tables remain unavailable in the checked sources. Winter/Luo animal-pool and position metadata remain unresolved.
- **BARseq2 M1/XC116 fallback:** it is a separate paired mouse, but detailed dissection/hemisphere definitions leave no verified unchanged A1 endpoint. Similar names are not sufficient to score a transfer test.

Source receipts are in `data/interim/paired_projection_search_round5/`. Project code, models and results stay in **andresloak100/SYN4RCH**. The external issues are data requests only; no further author contact is authorized.

A separate, explicitly synthetic multi-scenario benchmark is being built with disjoint simulated training/test animals. It will compare the frozen recipe and controls under known molecular, spatial, hidden-variable and measurement-error hypotheses. Only supported training-side anatomical statistics can anchor its physical units. Unknown molecular/projection couplings remain assumptions. Synthetic scores will not be recorded as real hypothalamic accuracy, de novo circuit discovery, synapses or complete 3D axon validation.

References for evaluation boundaries: [model-selection bias](https://www.jmlr.org/papers/v11/cawley10a.html), [probability calibration and proper-score limitations](https://scikit-learn.org/1.8/modules/calibration.html).
