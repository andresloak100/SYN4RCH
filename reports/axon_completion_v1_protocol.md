# Axon completion v1: protocol (TRAIN-only synthetic-masking development experiment)

STATUS: FINAL, 2026-11-12 (UTC), after ROOT's concrete contract/protocol check and before the full TRAIN run.
Code: `analysis/axon_completion_v1.py`. Tests: `tests/test_axon_completion_v1.py`.
Output: `data/interim/axon_completion_v1/`.
Command: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /opt/anaconda3/bin/python -I analysis/axon_completion_v1.py run --workers 4`.

The E4 evaluation and selector v2 are closed and unchanged. This experiment reads their frozen code read-only (the v2
library, cache and E1/E2/copy rules, pinned by SHA-256) and writes nothing of theirs.

## 1. Question and status

Each query reveals only the first 300 um of root-path distance of its TRAIN axon arbor. Can a prefix-conditioned,
learned donor choice, with one shared frontier-graft decoder, complete the hidden cable better than every control
(E1/E2, metadata-only choice, copying, and the strongest prefix-aware copying) on BOTH hidden distance and hidden F1,
by a practical margin and with valid, non-inflated cable?

This is **synthetic-masking development completion**. It is not whole-neuron or synaptic-circuit prediction.
The 101 TRAIN brain samples have already informed earlier designs. Validation and test samples are never opened.
A pass would need fresh independent data before any claim. A fail is reported as a fail. There is no viewer or
website promotion based on development scores.

## 2. Data, observation and truth

* **Cohort:** the 4,368 included TRAIN neurons of 101 brain samples (v2 `Library`; file hashes checked).
* **Observation (ROOT's contract, `analysis/axon_completion_contract.py`, read-only, finalized by ROOT and pinned in
  the code, SHA-256 `13e6cd60b2021e0d673407ec7aaf5906a3652b41310a11a41053148f9c41c8bd`):** `split_arbor` cuts
  each arbor exactly at 300 um of root path. Crossing edges are split by exact interpolation, and no edge beyond a
  boundary vertex is exposed.
  * The public `Observation` carries coordinates and parents only, in a visible-only canonical order. It has no
    source IDs, types, terminal flags, hidden topology or whole-arbor summaries.
  * The observation keeps the trace-root coordinates. The table-minus-trace soma offset is reported (`cohort.json`)
    and used only as a known scalar.
* **Truth:** `MissingTruth` is never an argument of a predictor, feature or decoder; the decoder refuses it.
  * Scoring jobs (`label_job`) create it next to `score_completion`.
  * The acquisition/index stage may construct and discard it: `observe` calls `split_arbor`, and the donor-side cut
    reads the whole TRAIN arbor. Only public observation summaries are kept as query-side fields.
* **Isolation is code-guarded per fit, not physical absence.** Every TRAIN arbor (and so every outcome) is on disk
  and can be loaded by a worker. Each fit reads labels, cache D and centrality only through the fit-only door, and
  each outer fold makes one refused probe read. Outer-eval forecasts are sealed before they are scored.
* **Fixed exclusion rule:** a query is excluded as fully observed when its contract hidden cable is at most 1e-9 um.
  It is decided in the index stage, before any score, and the count is reported. The same rule makes a donor
  ineligible (it has nothing to graft).
  * ROOT's structural preflight found 0 such exclusions in all 4,368 neurons / 101 samples, with cable and
    quadrature conserved (no model scores). Receipt: `data/interim/axon_completion_v1_root/structural_preflight.json`.
* **Forbidden inputs:** the query's hidden suffix, its regions, endpoints, terminal flags, whole-arbor summaries and
  any outcome.
* **Query inputs:** the public observation, plus the known metadata of v2 (table soma, peptide, nucleus,
  hemisphere).
* **Donor side** (fit animals only): the whole donor geometry is allowed. This covers its suffix components, the v2
  morphology summaries, the suffix cable, and centrality from the cached D (fit-only door).

## 3. The shared decoder and cable guard (identical for every method)

1. **Cut:** the donor is cut at the same 300 um budget. Each donor boundary vertex gives one suffix component, rooted
   at the exact boundary point, with the visible direction into it. Parity with the contract's split is tested.
2. **Matching:** query frontiers (`Observation.frontier_indices`, visible geometry only) and donor components are
   matched by cost = |frontier − component root| + 100 um × (1 − cos of the two incoming directions).
   * The Hungarian assignment covers min(m, k) pairs. Any remaining components go to their cheapest frontier.
   * Every component is used exactly once: no cable multiplication. Frontiers without a component stay unserved.
3. **Graft:** the component root is MERGED into the observed frontier vertex. There is no connector edge and no
   bridge to the forest root.
   * **Shared-trunk merge:** while each component is traversed once, a generated child whose parent is the same new
     node and whose coordinates are exactly identical to an already generated child is merged into it. This prevents
     two parallel components relocated onto one frontier from emitting the same branch twice (ROOT's fixture
     `data/interim/axon_completion_v1_root/duplicate_decoder_fixture.json` is a regression test).
   * Different parents are never merged and observed geometry is never altered.
   * The removed duplicate cable is recorded per forecast and per method.
   * The frozen `compact_attachment` is tried at 500 um, then 1,000 um. It is kept only if the component's cable does
     not exceed its original cable.
   * Otherwise the component is relocated rigidly: it moves, and its length is exact.
4. **Guard and validity:**
   * the observation is preserved exactly (append-only `Completion`);
   * the total new cable must not exceed the donor's hidden cable, else the forecast is a `failed_cable_guard`;
   * empty and duplicate forecasts are flagged by the contract;
   * all new cable counts, including cable near the observation.
5. **Observation-only** (no new cable) is a diagnostic failure control.

## 4. Candidate set and control rules (per query, per pool)

A pool is always the call's fit animals minus the query's own animal, restricted to eligible donors.

* **Prefix costs:** a symmetric, arc-length-weighted chamfer between the 50 um samples of the two visible prefixes.
  * Absolute version: in atlas coordinates.
  * Soma-centred version: each prefix taken relative to its own root.
  * Held in memory only, never serialised.
* **Candidate union:** the 16 table-soma-nearest donors (exactly the engine's E1/E2 set), plus the 8 nearest by
  absolute prefix cost, plus the 8 nearest by centred prefix cost, plus the same-peptide prefix nearest.

| method | donor rule |
|---|---|
| `E1_local_uniform`, `E2_local_kernel` | exact v2/engine local medoids on the 16 soma-nearest (cached D via door) |
| `local_prefix_nearest` | matched same-pool prefix selection: lowest absolute prefix cost among the 16 soma-nearest |
| `prefix_abs_nearest` | full pool, lowest absolute-coordinate prefix cost |
| `prefix_centred_nearest` | full pool, lowest soma-centred prefix cost |
| `prefix_peptide_nearest` | same peptide, lowest absolute prefix cost (else the full pool) |
| `adaptive_fallback` | the exact fallback chosen for that outer fold from the six rows above, by inner groups only |
| `nearest`, `medoid`, `peptide_nearest`, `peptide_medoid` | v2/referee full-pool copying rules |
| `metadata_only` | the same learner and rule, 16 soma-nearest only, no observation-derived feature; its fallback is E1 or E2 (inner groups) |
| `switch_random` | on the selector's switched queries, a seeded random other candidate; elsewhere the fallback |
| `permuted_labels` | the selector pipeline with outer-fit residual targets permuted (negative control) |
| `observation_only` | no new cable (diagnostic failure control) |

Every method's donor goes through the decoder of section 3.

## 5. The candidate: prefix-conditioned residual donor choice

* **Folds:** 5 outer folds of whole brain samples, each with 4 inner folds (salt `axon_completion_v1|20261113`).
* **Fallback:** for each outer fold, the fallback is chosen from the six fallback rules by the inner-eval
  animal-mean joint loss `log1p(d/200) + (1 − F1)`.
  * Each inner-eval query draws donors from its inner-fit pool only.
  * Ties go in table order. Failures count as infinite loss.
* **Learner:** one fixed HGB configuration (v2 settings: max_iter 200, learning_rate 0.05, max_leaf_nodes 31,
  min_samples_leaf 50, random_state 20261113). There is no architecture or hyperparameter search.
  * It is fitted as two heads to the residual gains over the fallback donor:
    `gd = log1p(d_fb/200) − log1p(d/200)` and `gf = F1 − F1_fb`.
  * Rows are every candidate of every outer-fit query (pool = outer fit minus the own animal), with weights
    equalising samples.
  * Invalid training labels are excluded from fitting and counted.
* **Features:**
  * known metadata and pair terms as v2;
  * donor summaries, suffix cable and centrality;
  * the query's observation summaries (frontier count, prefix cable, frontier offsets, directions and spread);
  * prefix costs and ranks;
  * decoder-plan diagnostics (matching cost, shift, direction agreement, unserved frontiers, attachment-feasible
    share);
  * differences against the fallback.
  * A name check refuses identity-like columns.
* **Choice:** switch to the candidate with the largest predicted gd + gf among those predicted better than the
  fallback on BOTH heads. Otherwise emit the fallback.
* **Fit-only door:** labels and cache D are read through the frozen v2 `FitDoor`. Each outer fold makes one probe read
  of an outer-eval label, which must be refused.

## 6. Forecasts before scores

* For every outer-eval query, every method's donor is chosen from known inputs and the public observation only.
* Each is decoded from the public observation alone (`Observation.from_dict`, no truth in the process). The
  completion hashes are written to `forecasts.json` and sealed (`FORECASTS_SEALED.json`).
* Only then are outer-eval queries scored, in a separate `labels_score.jsonl`. The scored completion hash must equal
  the sealed hash.

## 7. Measures and the predeclared gate

* **Measures:** contract `score_completion` on new cable only (50 um sampling, 200 um tolerance): hidden distance
  and hidden F1.
* **Aggregation:** neurons within a brain sample, then samples equally.
* **Failures are never dropped:** a forecast that is not `valid_for_comparison` counts as distance +inf and F1 0 in
  its sample. This covers empty, duplicate-cable, non-finite and guard-failed forecasts, and it applies even when the
  contract reports finite diagnostic values for that forecast. Those diagnostic values are kept separately per method
  (`invalid_diagnostics`) and never enter an aggregate.
* **Intervals:** paired sample bootstrap 95% (10,000 draws, seed 20261110).

**PASS** requires ALL of the following. Otherwise FAIL, with the reasons listed.

1. **Validity:** the selector has 0 invalid forecasts (empty, duplicate, non-finite or guard-failed).
2. **Practical, resolved win over EVERY required control in section 4.** `observation_only` is a diagnostic
   failure control: it is reported but is NOT a required control.
   * distance(selector) − distance(control) ≤ −10 um, with the bootstrap upper bound < 0;
   * F1(selector) − F1(control) ≥ +0.010, with the bootstrap lower bound > 0.
   * A required control that has any invalid forecast, or a non-finite comparison, is UNEVALUABLE. The joint result
     is then FAIL (inconclusive), never a win: a single invalid case cannot waive a strong control. That control's
     failures stay in its denominators and are reported.
3. **Cable guard**, both parts:
   * the selector's share of cable-inflated forecasts (new cable > 2 × hidden cable) does not exceed any required
     control's share by more than 0.02;
   * the selector's sample-equal mean absolute relative cable error |new − hidden| / hidden is no higher than that
     of `adaptive_fallback`.
   * Every forecast counts: an empty forecast has error 1, and one with no recorded cable counts as +inf.
   * This closes the chance of buying F1 through moderately excessive cable.
4. **Switching:** the switch rate is above 0.

**Also reported:**
* the per-method means and validity counts, and the invalid forecasts' diagnostic values (kept separate);
* the cable ratio (median, mean log2), the mean absolute relative cable error, the merged duplicate cable and the
  executed graft modes;
* the fallback chosen per fold and the inner losses;
* the switch counts, the metadata-only result, and the door refusals.

## 8. Run mechanics

* `run` refuses unless: this file is FINAL; the owned files are committed and clean; and the contract SHA is pinned
  in the code after ROOT's check.
* **Before the run:** synthetic-mask, leakage and decoder tests pass, and a scratch smoke run (15 samples, never
  recorded as a result) completes.
* **The run:** one continuous run with 4 workers, all four BLAS thread variables set to 1, and resumable JSONL
  checkpoints (index, fit labels, decoded forecasts, score labels).
* **Pre-run smoke** (15 samples, scratch folder only, never a result; no setting changed because of its scores):
  * The pre-fix decoder emitted a few duplicate-cable forecasts in some controls. After the shared-trunk merge, every
    method had 0 invalid forecasts; only the diagnostic observation-only forecasts are failures, by design.
  * Measured throughput was about 80–110 scored pairs/s with 4 workers, plus about 200 s of in-memory prefix costs.
  * The one pure refactor made after the smoke: the permuted-label control reuses the selector's training rows,
    with identical results.
* **Leak errors:** leak errors from the frozen v2 `FitDoor` (a different exception class) are translated by the
  label door, and are caught together with this module's own class at the probe and in the run handler. Either
  class writes `RUN_FAILED`.
* **Prefix costs:** they are computed in memory by a centred matrix product (agreement with the direct formula
  within 4e-7 um on checked rows) and are never written.
* **Resuming:** a resume is refused if code, protocol or contract changed after `RUN_STARTED`. A leak writes
  `RUN_FAILED`, and after `RUN_COMPLETE` the run is never repeated.
* **Kept small:** no full pairwise matrix and no full-geometry audit is written. Only summaries, choices, hashes and
  per-pair metrics are kept.
