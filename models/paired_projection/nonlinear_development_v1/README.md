# Nonlinear development artifacts (EXPERIMENTAL, NOT PROMOTED)

These are exact copies of the two deployment artifacts written by run `nonlinear_a1_dev_v1`
(source commit `a76448a`). They are **not** the operational model. The operational model is still the frozen
logistic recipe in `models/paired_projection/logloss_recipe_v1.json` and its A1 weights in
`models/paired_projection/a1_logloss_v1.json`.

**Why not promoted:** on the exposed A1 development data, the nonlinear gene + depth model lost to the frozen
logistic model on both log loss (0.3251 vs 0.3190) and Brier (0.0960 vs 0.0952). See
`reports/paired_projection_nonlinear_results.md`.

**What they are:**
* One fixed scikit-learn 1.8 `HistGradientBoostingClassifier` per target, for the 11 captured targets:
  * 100 trees, depth 2, 4 leaves, min 20 samples per leaf, L2 1.0, learning rate 0.05, seed 20261112;
  * no early stopping and no class weights.
* Refit on all 598 A1 development cells from mice XC95 and XC97. These are `DEPLOYMENT_REFIT_NOT_OOF`
  artifacts, not out-of-fold estimates.
* Shrinkage toward the training capture rate: alpha 0.1 for the gene + depth model, alpha 0.25 for depth only.

| file | SHA-256 |
|---|---|
| `artifact_hgbt_gene_position.json` | `246788f5610e5e44e125b44d8a80c580c5b3dcf7f9dbd3b54c2ea3caaee24037` |
| `artifact_hgbt_gene_position.joblib` | `42eb7782b391a1273aec36849255f7c1f74f83bf710af1b3489917654bc9f3c5` |
| `artifact_hgbt_position.json` | `6b089e87bfdc4874c9d4c84c59407824baa72cd28da166cda0345e6e09ea5d05` |
| `artifact_hgbt_position.joblib` | `163aa53be00f6397894f5da6c3841eeac59300c08947f65ef6e1e13abfbed2b4` |

**Contents:**
* Each `.json` file holds:
  * the feature order (23 measured genes, log1p, then `depth_norm1200_um`; or depth only);
  * the training-only median and scaler, `q_fit`, the training cell IDs and groups;
  * the fixed parameters and alpha;
  * the SHA-256 of its sibling `.joblib` file.
* Keep the `.json`/`.joblib` pairs together with their names unchanged: the loader finds the model through
  `model_file`.

**Loading and prediction:**

```python
from analysis.paired_projection_nonlinear import load_artifact, artifact_predict
meta, models = load_artifact('models/paired_projection/nonlinear_development_v1/artifact_hgbt_gene_position.json')
probs = artifact_predict(meta, models, meta['gene_features'], meta['position_features'], gene_counts, depth)
```

* `load_artifact` checks the model file's SHA-256 before `joblib.load`. Only load these SYN4RCH files, never
  third-party pickles.
* Prediction takes features only:
  * numeric 2-D arrays with exact widths and matching rows;
  * finite values, and gene counts of at least 0.
* Probabilities are for the captured BARseq2 label: raw count above 2 counts as captured, and a lower count is
  nondetection, not absence.

**Limits:** A1 cortical assay only. Not hypothalamic, not validated on new animals, and no claim about 3D
geometry or synapses. A different gene panel would need a refit on its own admitted training animals.
