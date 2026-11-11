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
