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

(not yet: audit not started)
