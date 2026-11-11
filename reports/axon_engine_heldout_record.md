# Axon engine: held-out stage record (hashes committed as each step happens)

Protocol `reports/axon_engine_protocol.md` FINAL (commit aac2e87, SHA-256 `13e9b83c6f31412732da89d34ef0a39b2a5e59fa7b7af6b27cf06df6d6c97943`).
Engine `analysis/axon_engine.py` SHA-256 `d9eac6f220691b33adc729ecc24a39d3290418039b7717bb052900741da971a8`.
Referee audit receipt `2eb4d2e4...` (record: `reports/axon_engine_referee_audit_record.md`).

## Validation

1. Forecasts committed 2026-11-11 16:34 EST, before any validation trace was opened: 1,185 queries x 16 methods = 18,960.
   * `data/interim/axon_engine/real/validate/PREDICTION_COMMITMENT.json` SHA-256 `b410470141c8830cae3ced9a3b2e2f3927945aa7b765da4f1b8675ae0f483c0e`
   * `data/interim/axon_engine/real/validate/prediction_manifest.json` SHA-256 (the commitment) `8e2f7c15770757532c78e1fa193d27acb3384e7676bca3408f646c238e921d82`
2. Copy check 16:35 EST: the referee's `copying_baselines` through the audited cache; 4,740 official-copy forecasts
   compared, 0 differences. `COPIES_VERIFIED.json` SHA-256 `98b4a64848a4ff141a2cd7ef3a1b169ca95abe1f411f71b4077bb6dd45c570bf`.
3. Admission 16:35 EST (validation traces opened only now; 0 forecasts read): excluded by rule 2: 60, rule 4: 3, rule 6: 1;
   0 unavailable at source. `ADMISSION_COMPLETE.json` `557ec889df327dbdab75fc29b10f9add9122f18ecfa98cff8439fed73c734186`,
   `TRUTH_OPENED.json` `2ef6941e9ddbb09486095d89c39035ec341953a1e599629775f0b8335fb1d418`,
   `truth_ingestion_log.json` `68ead63594e458d7a745c74ae630a6c4f4b5055f772b45ca4c909273182e9e37`.
4. Scored ONCE 16:37 EST: 1,121 of 1,185 queries, 24 of 25 animals. Selected: **E4_local_kernel_regional**
   (candidate_selected): 631.9 um / F1 0.2728 vs strongest copying control on distance C5 711.3 um (-79.4 um, paired
   95% bootstrap -108.2 to -50.2) and on F1 nearest 0.2243 (+0.0485, +0.0308 to +0.0672); beats its matched controls.
   Eligible also: E1, E2, E6. Validation is selection, not the test.
   Caveats (stated before the test): 64 of 1,185 declared validation neurons (5.4%) were excluded by the mechanical rules,
   and one whole validation brain sample has no scored neuron (24 of 25 scored). E4 beats its matched control E1 by only
   0.009 um in distance and 0.00155 in F1, and both matched-control intervals cross zero: most of the gain over copying
   comes from local donor selection; an added benefit of the regional state is not established. E7 was not eligible
   (its distance did not beat its matched controls). Root's independent check of the arithmetic and the selection rule
   passed: `data/interim/axon_engine/real/checks/validation_arithmetic_gate_check.json` SHA-256
   `31b536a0290c4098aed33ab3641aff22e764edcb3ba88c6a4069b2d554d270a5`.
   * 8e1b4c2d7a762fc979879ee7938b14f5dd3a40b387d08e4f056c8b5a1bfb04fc  SCORE_COMPLETE.json
   * 42516e20b33d55b02d803332c1ca71a45df6f435fab8dd66aa5f88fbf2d1f908  metrics.json
   * beef20acbd06182b487d22758bb6509527556cde96207622fbc38ec2c807c790  engine_gate.json
   * 6935aabeb8ae6ccd262c1d84fe0cf8747723b69af4afb20a43c8ad6cf9afeba2  frozen_gate.json
   * 6e8279cd8ee1910a92ab7fd8800ae8d86758f6b6fe9986ed977dcd21c57d4956  tables.md

## Freeze

Frozen 2026-11-11 16:38 EST, with 0 test traces on disk and no test stage directory. Selected: E4_local_kernel_regional
(the frozen prototype's own gate chose conditional_anchor_1000 among its five policies; recorded, not the selected method).
* `data/interim/axon_engine/real/freeze/frozen.json` SHA-256 `a05ac68b58dd9bde92f535c4e097e573b5c35a20696d9cdccc3710c1e50b0202`

## Test

Forecasting for all 1,547 declared test queries x 16 methods started 16:38 EST, before any test trace is fetched.
