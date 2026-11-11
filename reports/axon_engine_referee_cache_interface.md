# For the referee's owner: the engine's training-only distance cache, and when the referee may change

2026-11-11. Read-only information; the referee's files are not edited from this side.

## The cache (being built now from the 4,448 training traces)

* Folder: `data/interim/axon_engine/real/paircache/<first 16 hex of the library fingerprint>/`.
* `index.json`: schema `axon_engine_pair_cache_v1`; `ids` (training neuron IDs, sorted); `arbor_sha256` (content hash
  of each canonical 50 um arbor, same order); `ids_sha256`; `library_fingerprint` = SHA-256 of the canonical JSON of
  `[[id, arbor_sha256], ...]`; `metric` = the definition plus the SHA-256 of `analysis/d110_geometry.py` and
  `analysis/brm_template.py`.
* `D.npy`, `F.npy`: float64, n x n, symmetric, zero diagonal. `D[i, j]` = `distance_um` and `F[i, j]` = `1 - f1_200`
  of the frozen `geometry_metrics(a.shifted(0), b.shifted(0))`, soma-centred, lower neuron ID as first argument
  (exactly what the frozen `pair_costs` computes). `K.npy`: uint8, 1 where a pair is computed (resumable build).
* `COMPLETE.json`: the SHA-256 of `index.json`, `D.npy`, `F.npy`, and `cache_sha256` = SHA-256 of the canonical JSON of
  {schema, library_fingerprint, ids_sha256, n, pairs, metric, index_sha256, D_sha256, F_sha256}. A cache without it
  is incomplete and must not be used.
* Training arbors only: no validation or test neuron is in it. A fixed sample of pairs is re-proved against the frozen
  `pair_costs` and against the referee's own loader and metric every time forecasts are made.

What would remove the repeated 8-hour recomputation: `copying_baselines` (and the copy step inside
`seal_predictions`) accepting a pinned matrix like this one, checking `cache_sha256`, the arbor hashes of the
donors it was given and a random sample of pairs against its own metric, and refusing on any mismatch.

## When the referee source may change

* Open now, until the engine commits its VALIDATION forecasts (`predict --stage validate`, then `copies-check`).
  From that point the referee's SHA-256 is part of the commitment, the freeze and the seal; a change afterwards makes
  the engine refuse to admit, freeze or seal, and would cost the validation run.
* Expected order here: training ingestion and the cache (about 3 to 4 hours from 08:15 EST), development on the
  training animals, then validation. I will say before validation starts. After the test is scored, it may change again.
