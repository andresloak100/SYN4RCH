# PVH collection-label pilot artifacts (research only, not deployed)

These are exact copies of the two fold artifacts written by run `pvh_population_pilot_v1`. They are research
artifacts for an exploratory PB-versus-SC collection-label pilot. They are not deployed models and are not
promoted. No new refit was made for this archive.

**What they are:**
* One artifact per predeclared cross-sex pool fold:
  * `fold_F_to_M` was trained on the two female pools;
  * `fold_M_to_F` was trained on the two male pools.
* Each `.joblib` holds the fitted state of all seven methods, the 30,477-gene common panel, the training-only
  500-gene selection and scaler, the pool-balanced prior, and the training cell keys.
* The matching `fold_*.json` lists the same fold settings in readable form.
* Each artifact predicts membership of the PB retrograde collection group, conditional on a nucleus being in
  this PB-or-SC collected cohort. That is not a projection probability. Other targets are unknown, and a pool
  is 4 to 10 pooled mice, not one animal.

| file | SHA-256 |
|---|---|
| `models_fold_F_to_M.joblib` | `66211535bcae9c784ded379484a0e248681f55da3b7699fa353628998cea153d` |
| `fold_fold_F_to_M.json` | `87f1e6043946fb8fb676f6305a31faaca28c5c5a625662d28428825f55869ca4` |
| `models_fold_M_to_F.joblib` | `aa7bceab8d1089d1c4c2343389e7c7047b38c944a2fc204a8201cf2cec6999c9` |
| `fold_fold_M_to_F.json` | `d6150fa4c602be0279539729562f07f13555e3bd3734e9aad7bc9739ebfa5943` |

**Loading:** only these SYN4RCH files, after a hash check; never third-party pickles.

```python
from analysis.pvh_population_pilot import load_artifact, artifact_predict
art = load_artifact('models/paired_projection/pvh_population_pilot_v1/models_fold_F_to_M.joblib',
                    '66211535bcae9c784ded379484a0e248681f55da3b7699fa353628998cea153d')
raw, shrunk = artifact_predict(art, gene_names, raw_counts)   # features only, exact panel order
```

**Results:** see `reports/pvh_population_pilot_results.md`. No method beat every control on both log loss and
Brier.

**Exposure:** all 265 nuclei and their labels are now exposed, so these artifacts may not be scored on them as
a fresh test. See `docs/pvh_population_pilot_exposure.json`.
