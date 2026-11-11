# Freeze record: SYNARCH reconstruction benchmark (before the single run)

2026-11-11. Written before `run` was started; no real outcome had been fitted or scored.

* Declaration: `data/interim/synarch_recon/pre_run_declaration.json` (git-ignored; its hash is
  the anchor), SHA-256
  `29ca3cf54db93b75f5a36ebcbb220196e55b6243782b44f8f77807f9be5cd2c6`.
* It pins the code and protocol of commit e7c0474, the inputs, the folds (032d4b79...), the grids and seeds.
* Scope: 101 training brain samples, 4,448 neurons; tasks T1 to T4; exploratory, not an independent test.
* The run is started once with `python3 -I analysis/synarch_recon_benchmark.py run --workers 3`.
