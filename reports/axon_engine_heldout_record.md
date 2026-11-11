# Axon engine: held-out stage record (hashes committed as each step happens)

Protocol `reports/axon_engine_protocol.md` FINAL (commit aac2e87, SHA-256 `13e9b83c6f31412732da89d34ef0a39b2a5e59fa7b7af6b27cf06df6d6c97943`).
Engine `analysis/axon_engine.py` SHA-256 `d9eac6f220691b33adc729ecc24a39d3290418039b7717bb052900741da971a8`.
Referee audit receipt `2eb4d2e4...` (record: `reports/axon_engine_referee_audit_record.md`).

## Validation

1. Forecasts committed 2026-11-11 16:34 EST, before any validation trace was opened: 1,185 queries x 16 methods = 18,960.
   * `data/interim/axon_engine/real/validate/PREDICTION_COMMITMENT.json` SHA-256 `b410470141c8830cae3ced9a3b2e2f3927945aa7b765da4f1b8675ae0f483c0e`
   * `data/interim/axon_engine/real/validate/prediction_manifest.json` SHA-256 (the commitment) `8e2f7c15770757532c78e1fa193d27acb3384e7676bca3408f646c238e921d82`
