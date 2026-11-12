# Next donor selector v2: protocol (TRAIN-only development experiment)

STATUS: FINAL, 2026-11-12, committed before the TRAIN run. Code: `analysis/next_donor_selector_v2.py`.
Tests: `tests/test_next_donor_selector_v2.py`. Output: `data/interim/next_donor_selector_v2/`.
Command: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /opt/anaconda3/bin/python -I analysis/next_donor_selector_v2.py run --workers 4`.

The completed E4 evaluation (`reports/axon_engine_heldout_record.md`) is closed and unchanged. This experiment reads
its files read-only. The old draft (`reports/next_donor_selector_protocol.md`, `analysis/next_donor_selector_feasibility.py`,
`data/interim/next_donor_selector/`) is preserved and not used.

## 1. Question and status

Can a learned choice among the same 16 local donors that E1 and E2 use beat exact E1, exact E2 and the exact
inner-selected (adaptive) E1/E2 fallback on BOTH animal-mean distance and animal-mean F1 at 200 um? The emitted
geometry is always produced by the same frozen attachment operator.

This is a **development** result. The 101 TRAIN brain samples have already informed the design (engine
development, earlier selectors, the general layer). The 25 validation and 40 test samples are used: this
experiment never opens or scores them, and never calls anything here an unseen test. A development pass would need
fresh independent data before any claim. A fail is reported as a fail.

## 2. Data and inputs

* **Cohort:** the 4,368 included TRAIN neurons of 101 brain samples (`ingestion_log.json`, part train, status
  included). The canonical arbors are `data/interim/axon_engine/real/library/<id>.npz`; each file hash is checked.
* **Query inputs (the only ones):** table soma (AP, DV, ML, mirrored-left frame, as the engine), source label
  (peptide), soma label (nucleus) and hemisphere.
  * `reconstruction_type` is excluded because it describes the trace, not the cell.
  * IDs and brain samples are used for grouping only, never as features. A name check refuses identity-like columns.
* **Never prediction inputs:** the query's measured geometry, endpoints, regional profile or any outcome.
* **Donor side** (TRAIN only):
  * morphology summaries of each donor arbor: log cable, branches, terminals, RMS extent per axis, cable-weighted
    centroid offset, maximum radius and node count;
  * donor centrality: the animal-balanced mean cached D to the fit pool, from fit animals only;
  * the donor's trace soma, used only to compute the attachment's soma shift.
* **Pair and local-structure features:** soma difference and distance, rank, soma shift and ladder feasibility,
  same peptide/nucleus/hemisphere flags, E1 and E2 weights, relative E1 and E2 medoid risks, and whether the donor
  is the E1 or E2 choice.
* **Cache:** the engine's complete TRAIN cache supplies D only (soma-centred, verified). It is used for the E1/E2
  local medoids, the copy medoids and centrality, which is how the engine and the referee use it. The cached F is
  never read.

## 3. Candidate set, emitted geometry and labels

* **Candidate set:** the 16 donors of the call's pool nearest by table soma. Ties go by ID, and the set is sorted by ID
  (exactly the engine's E1–E7 set). The E1 and E2 choices lie in this set, so the selected fallback donor is always
  a candidate.
* **Emitted geometry** of any local method for donor d: the engine's ladder, i.e. frozen `compact_attachment` at
  500 um, then 1,000 um, else translation of the same arbor, with the root at the query's TABLE soma. Distal
  coordinates are preserved exactly. The selector, E1, E2, the adaptive fallback and every local control use this
  same operator.
* **Labels:** exact `d110 geometry_metrics(emitted, measured query arbor)`: distance_um and f1_200, recomputed for
  every (query, candidate) pair that any fit step reads. That union is 186,757 pairs over all 50 calls. No cached
  label and no full-pair audit is used.

## 4. Nested whole-sample cross-validation

* **Outer folds:** 5 outer folds of whole brain samples, ordered by salted hash (`next_donor_selector_v2|20261112`).
* **Inner folds:** 4 inner folds within each outer fit set.
* **Pools:** every call draws donors only from its fit samples, minus the query's own sample.
* **Fit-only door:** labels, cache D and centrality are read only through `FitDoor`. Both endpoints must be in the
  call's fit samples, and a same-sample pair is refused. Each outer fold makes one probe read of an outer-eval label,
  which must be refused. Any other refusal is a `RUN_FAILED`.
* **Inner choices** (fixed small search):
  1. The fallback, E1 or E2 (no parameters), is chosen by the inner-eval animal-mean joint loss
     `log1p(d/200) + (1 - F1)`. Ties go to E1.
  2. The learner, from {ridge alpha 100, ridge alpha 1, one HGB}, is chosen by the same objective of the combined
     rule. Ties go in that order.
  * HGB settings: max_iter 200, learning_rate 0.05, max_leaf_nodes 31, min_samples_leaf 50, random_state 20261112.
  * Heads: `log1p(d/200)` and F1.
  * Pair weights equalise samples.
* **Combined rule:** take the top candidate by predicted log-distance minus predicted F1 (ties by ID). Keep it only
  if it is predicted better than the fallback donor on BOTH heads. Otherwise emit the fallback, whose exact arrays
  are checked to be identical to the adaptive-fallback control. Non-finite predictions fall back.
* **Outer step:** refit on all outer-fit samples. Forecast the outer-eval samples from known inputs only. Write the
  donor choices of every method and the hashes of the emitted arrays, and seal them (`FORECASTS_SEALED.json`). Only
  then are outer-eval traces scored. The scored geometry hash must equal the sealed hash.

## 5. Controls (same queries, same outer-fit pool)

| control | rule |
|---|---|
| `E1_local_uniform`, `E2_local_kernel` | exact engine arithmetic (parity: 919/919 donors, arrays and metrics bit-identical with the engine's committed development forecasts) |
| `adaptive_fallback` | the unchanged inner-selected E1 or E2 of that outer fold (added after a focused review: a pass by adaptive fallback alone, with harmful learned switches, must fail) |
| `nearest`, `medoid`, `peptide_nearest`, `peptide_medoid` | the referee's copying rules on the full outer-fit pool, translated (parity: 3,676/3,676 donors) |
| `local_nearest`, `local_label_nearest` | matched same-pool simple selection (nearest soma; nearest same-peptide, else nearest), attached |
| `switch_nearest`, `switch_random` | on exactly the switched queries: the soma-nearest, or a seeded random candidate other than the fallback; elsewhere the fallback |
| `permuted_labels` | same pipeline and chosen learner, outer-fit labels permuted (negative control) |

Descriptive only: the best of the 16 by truth (an oracle ceiling), and the fixed ablation listed in section 6.

## 6. Measures, gate and reporting

* **Aggregation:** neurons within a brain sample, then equal samples.
* **Intervals:** paired sample-bootstrap 95% (10,000 draws, seed 20261110).
* **Development gate:** the selector has lower distance AND higher F1 than EVERY control in section 5, and its switch
  rate is above 0.
  * The main question is E1, E2 and the adaptive fallback.
  * A copying-only or fallback-only result is not a learned-selector win.
* **Reporting sentence:**
  * "development gate met and resolved against every control" when every interval favours the selector on both
    measures;
  * "met, not resolved" when the point gate holds but some interval does not;
  * otherwise "development gate not met".
* **Also reported:** the switch rate and combined-rule branches, the inner choices by fold, per-label coverage and
  per-label means, the switched-subset means, and the door read counts.
* **Fixed ablation:** selector, E1 and E2 are also scored translated.
  * Attachment effect: attached minus translated, per method.
  * Donor-choice effect: selector minus E1 and selector minus E2, both translated.

## 7. Run mechanics

* **Before the run:** this file is FINAL; code, protocol and tests are committed and clean; parity with the engine's
  development forecasts must hold, otherwise nothing runs.
* **The run:** one continuous run with 4 workers on one long-lived worker pool. Every thread variable is 1.
* **Resuming:** labels are written to resumable JSONL checkpoints (`labels_fit.jsonl`, `labels_score.jsonl`), and
  progress, rate and ETA go to `run.log`. A resume is refused if the code or protocol changed after `RUN_STARTED`.
* **Once only:** after `RUN_COMPLETE` the experiment is never repeated. A leak writes `RUN_FAILED`, and a new attempt
  needs a new, disclosed version.
* **No full-matrix serialization:** only the summaries, choices and hashes are written.
* **Pre-run check:** a 15-sample end-to-end smoke run (scratch folder only, not recorded as a result) checked the
  pipeline. No setting was changed because of it.
