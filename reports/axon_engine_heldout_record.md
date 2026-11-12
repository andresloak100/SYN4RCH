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
* Test forecasts committed 2026-11-11 17:11 EST, with no test trace on disk: 1,547 queries x 16 methods = 24,752.
  * `data/interim/axon_engine/real/test/PREDICTION_COMMITMENT.json` SHA-256 `e913f14ead130b2c127156911c4004e9cae76d4ec33f00413d20b926a416859d`
  * `data/interim/axon_engine/real/test/prediction_manifest.json` SHA-256 (the commitment) `3268960a44bf86ac48ffb763efc818a55eff4a6baa81b3e2d2d418ac58e82f29`
* Copy check 17:13 EST: 6,188 official-copy forecasts compared with the referee's copying_baselines, 0 differences.
  * `data/interim/axon_engine/real/test/COPIES_VERIFIED.json` SHA-256 `cb85c97c3044145daa1086a4902ddaab2947beac389f60e496cef53567ab85d7`
* Stage B (Root): all 1,547 reserved test files from 40 brain samples downloaded, 0 unavailable, 2,233,127,294 bytes.
  Root's independent byte audit PASSED (every test file and all 24,752 forecasts rehashed; Stage A metadata preserved;
  acquisition lifecycle and frozen proofs verified; no Root geometry parsing or scoring).
  * `data/interim/axon_trace_download_audit/stage_b_pre_acquisition_pin.json` SHA-256 `9909bd8a6edf401ae5c071ecd2e1b1314524c45fe58587dc976670989caaaee4`
  * `data/interim/axon_trace_download_audit/stage_b_audit.json` SHA-256 `a035e99ccc9826e745984423db1752d3de45a830b067b832604b5e7faf0a6486`
  * `data/interim/axon_trace_download_audit/stage_b_post_acquisition_pin.json` (Root acceptance) SHA-256 `5c9f4b3b34945f548a1d4000fe7bd77786dcbd1b51685b548fb5e446a89ac63c`
  * Root's metadata duplicate check over all 7,180 raw receipts: no identical bytes across splits or across brain
    samples. Two same-sample identical pairs exist (TRAIN 211275_053/054; TEST 202571_013/095). Limitation recorded;
    the fixed admission rules were not changed and both test neurons were admitted and scored.
    SHA-256 `5d61eb7ed77f8ccd24ca1a249944b21800fa23e18533fadc140b034a4d7e642b`
  * Before admission, freeze, test commitment, engine (`d9eac6f2…`) and protocol (`13e9b83c…`) hashes re-checked: unchanged.
* Admission 2026-11-11 19:04 EST (fixed mechanical rules, 0 forecasts read): 1,475 of 1,547 neurons admitted, 39 of 40
  brain samples. Excluded 72: rule 2 (soma mismatch): 9, rule 4 (axon coded as two types): 5, rule 5 (Axon_only file has
  dendrite nodes): 51, rule 6 (axon cable too short): 8 (one neuron broke rules 5 and 6). Brain sample **201228** is lost
  entirely: all 51 of its neurons broke rule 5. 0 unavailable at source.
  * `ADMISSION_COMPLETE.json` `73399c50222782f2407adce8e15d0a6179c157abe2766811f3ef06c2af691182`
  * `TRUTH_OPENED.json` `2bee455632266ca1062688da1e081795e853dee286e0e86c4e198523f3ae3307`
  * `truth_ingestion_log.json` `216695221f9288f938f03f05209f48081453d7f8be370fd2fd806cdd77b16902`
* Seal 19:06 EST (one attempt): `REFEREE_SEAL.json` `980ca6d09241651df38ae1eb557e2879f11bec2632940f8678a7d29d86f3e151`,
  `referee/PREDICTION_SEAL.json` `9d7b1f08f893445524f1b80bf922eba45cc0627a419edd04dc4de2c20d67488a`,
  `SEAL_ATTEMPTS.json` `88c46be45c555990be819dbe5846805c64159ce4fc51b6d7202cfa73f2e39091`.
* **Scored ONCE 19:11 EST by the independent referee: RULE NOT MET.** Selected E4_local_kernel_regional: 649.9 um,
  F1 0.2708 (animal means, 39 animals, 1,475 neurons).
  E4 minus control (negative distance and positive F1 favour E4; paired 95% animal bootstrap):

  | control | role | distance um | F1 |
  |---|---|---|---|
  | nearest | copy | -190.4 (-215.8 to -165.3) | +0.0586 (+0.0478 to +0.0699) |
  | medoid | copy | -88.2 (-124.9 to -51.7) | +0.0590 (+0.0459 to +0.0718) |
  | peptide_nearest | copy | -145.3 (-188.7 to -103.0) | +0.0488 (+0.0344 to +0.0632) |
  | peptide_medoid | copy | -146.9 (-195.3 to -105.7) | +0.0879 (+0.0683 to +0.1071) |
  | C1_medoid_copy | copy | -88.2 (-124.5 to -51.8) | +0.0590 (+0.0458 to +0.0721) |
  | C4_uniform_attach_500 | matched | -88.2 (-124.6 to -51.0) | +0.0590 (+0.0460 to +0.0718) |
  | C5_uniform_attach_1000 | matched | -88.0 (-124.9 to -50.3) | +0.0584 (+0.0456 to +0.0712) |
  | E1_local_uniform | matched | **+12.3 (-3.7 to +30.0)** | **-0.0016 (-0.0067 to +0.0036)** |
  | E2_local_kernel | matched | **+15.5 (-0.9 to +34.1)** | **-0.00002 (-0.0043 to +0.0043)** |

  Wrapper verdict: `success = false`, "rule not met: better than every copying control, not better than its matched
  control"; frozen selection E4 (freeze a05ac68b) equals the referee's candidate slot. Raw referee: `beats_all_copying =
  true`, `beats_all_matched = false`, status `scored_once`, `changed_since_the_freeze = false`.
  Plain reading: local donor selection (E1/E2) beats copying by about 90-190 um and about +0.05 to +0.09 F1 on unseen
  animals; the conditional/regional part that E4 adds shows no benefit (numerically slightly worse, intervals cross zero).
  No viewer promotion follows from this test. Unselected candidates (descriptive only, no claim): E6 644.0 / 0.2828,
  E7 635.1 / 0.2767, E3 654.5 / 0.2704, E5 666.2 / 0.2698, F1 738.1 / 0.2118, F2 737.9 / 0.2123.
  * `SCORE_COMPLETE.json` `6219aa6a31a6909f4b0d12eae010ede58ac3aed0598ab9e61cd111fe78a1eef1`
  * `metrics.json` (wrapper) `ae18960e168636b9af9b394b9bc7b20098d2a297de594f3cd8e65d5eff2e351b`
  * `referee/metrics.json` `6f87f6c0e95a5246a5d4dc0c6559ac63f697e76fa958ed04bd6d05fff214b153`
  * `referee/declaration.json` `722712fe4fc686625170ba030bc5edfd28e981f5219dec98cde996f8cdb51ae2`
  * `referee/COHORT_FROZEN.json` `9f50f5722c638bb250b66802ecb2bdb3eae81a9e347ae8ee15db3101354ef4fc`
  * `referee/SCORE_STARTED.json` `b066323ee03f3903c497b5340a53bc11e6a4fd4c0a1f670d1b537c4279ca8c76`
  * `referee/SCORE_COMPLETE.json` `4fb245d36c16c6e5849f7a8a82d64d5db8cac11435e2eb7fc1464a2db97fba39`
  * `tables.md` `efe54e56826a961755f0c6fe6593cbc33f4888958fe3ac81f51a52fddf721e05`
  * Logs: `data/interim/axon_engine/logs/test_{admit,seal,score}.log`.
* Next (Root): independent check of the full wrapper gate.
