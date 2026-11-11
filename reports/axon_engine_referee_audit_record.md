# Referee audit of the training distance cache: record (pre-audit pins)

2026-11-11, 12:40 EST. Written before the referee's full audit was started. Training data only.

* Native cache: `data/interim/axon_engine/real/paircache/85f0f247bacfd142/`, all 9,537,528 pairs of the 4,368
  included training neurons (80 excluded by the mechanical rules). Library fingerprint
  `85f0f247bacfd142acc7fe7b9e73ad4b958c4c8d065f951c9eb8f8abd6030bc4`; pair-cache SHA-256
  `b8c089ee46f1915cc71696e33e89d0d11242f72f4e3946c5c72f169a5b930bd0`. Verified by the final engine (commit d790bef):
  24 sampled pairs bit-identical to the frozen pair_costs and to the referee's metric.
* Fresh audit folder, prepared by `axon_engine_cli.py referee-cache --prepare` into a folder that did not exist before
  (no receipt, progress file or verified mask): `data/interim/axon_engine/real/referee_cache/bundle/`

| file | SHA-256 |
|---|---|
| draft_declaration.json | 00d34d771262d78b2b9b8e3f1c3af7398d1a7908f5a1587013c184a5d12afd5e |
| cache/D.npy | fd197604091af2a712d9a9a31de1ae81484c7f860c48612ae7d1766b8363f165 |
| cache/F.npy | 925e543df708f832d1aa59991e9cabcfbb46da7df8c2f5b7b42b41e7fef54145 |
| cache/K.npy | fe3c80c2049196cc2a8291a86fab5db680d403578c85b812e09e0bdf90037ab7 |
| cache/COMPLETE.json | 8a714591e96ac635230870cb849f0847b094b12ca2d65759d2708bc2e6ce269b |
| cache/index.json | 292dc49a5626849f551de18da0f98730253b5c86226752cdd6e9c480468a3a50 |
| training/ (4,368 arbors): SHA-256 of the sorted `shasum -a 256` listing | 8fe411d281d6daa60b8d0267e717b7afc448061257a690a7ee25bd8a91c72279 |

The full recompute starts only after the referee's owner pins this folder. The completed receipt's SHA-256 is
added below before any validation forecast is committed.

## Audit receipt

Pre-audit pin by the referee's owner: `data/interim/axon_trace_download_audit/referee_cache_pre_audit_pin.json`,
SHA-256 `94bb1746886a0dee9d959e1ba7fa644355882ce445565b321ee8ee6f0a5fd2f5` (all six files, all 4,368 arbors and their arrays,
sorted IDs, disjoint splits, no prior receipt, progress or mask). Helper `analysis/prediction_distance_cache.py` SHA-256
`88461504a10fd92b8d093f82cf440bb82fa668f427c017644469b9c3bc2893c4`.

Full audit started 2026-11-11 12:45 EST on this folder, 4 workers, the four thread variables at 1, `/opt/anaconda3/bin/python`
(3.12.4) for every chunk; log `data/interim/axon_engine/logs/referee_audit.log`. Scope: the helper replays every D entry
bitwise; F is checked by hash and structure only, so F values are not independently verified gold.

Receipt: complete, 2026-11-11 15:39 EST (chunk 6). `data/interim/axon_engine/real/referee_cache/bundle/cache/full_audit.json`

| item | value |
|---|---|
| receipt file SHA-256 | `2eb4d2e463723eae88a6e2321c870850f9146d9e1e20d0c79fba8d0f72589498` |
| receipt content SHA-256 (`receipt_sha256`) | `8a60d8f5447f5b69d457f36b80b302d1b9c960ef9bd139bd11ec3d01c02216ac` |
| schema / status | `referee_full_distance_audit_v1` / `complete` |
| pairs verified | 9,537,528 of 9,537,528 (all upper-triangle D entries, bitwise) |
| bitwise mismatches | 0 |
| scope | `ALL_upper_triangle_D_bitwise;F_hash_and_structure_only` (F not independently recomputed) |
| bound files | the six pre-audit pins above, unchanged |

Independent post-audit verification by the referee's owner (actual loader):
`data/interim/axon_trace_download_audit/referee_cache_independent_post_verification.json`, SHA-256
`01f48bc5f71f0cdea43a1705538a8ed130266cb7fdd964275831411a7b0888ed` (verifier `analysis/referee_cache_post_audit_verify.py`,
SHA-256 `b70b9bfe906c2455bccfbcb5d92878041fca8fd0d447a7c9cce7decbf066c9c7`). All six artifacts, five sources and the 4,368
training NPZ bytes unchanged; the frozen loader accepted the cache; no test trace on disk.

Post-audit pin by the referee's owner: `data/interim/axon_trace_download_audit/referee_cache_post_audit_pin.json`, SHA-256
`f1b322653efb8e4e0fa2b3553b7810e73a11dfe5d1b26b24d370d7b6562c0ebb` (acceptance PASS, 2026-11-11).

Committed before any validation forecast. Engine at this point: `analysis/axon_engine.py` SHA-256
`d9eac6f220691b33adc729ecc24a39d3290418039b7717bb052900741da971a8` (commit 5bbd4dd).
