# Protocol 2 (DRAFT for review): a donor selector trained on all 101 training brain samples

Status: DRAFT, revision 2 (after review round 1), 2026-11-11. Not frozen, not run. Nothing here changes the
current evaluation (`reports/axon_engine_protocol.md`, E1 to E7), its code, its cache, its referee or its schedule.
New files only: `analysis/next_donor_selector*`, `tests/test_next_donor_selector*`, `data/interim/next_donor_selector/`.
Feasibility numbers: `data/interim/next_donor_selector/feasibility_all.json`, written by
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -I analysis/next_donor_selector_feasibility.py all`
(training data only; 143 s). Tests: `python3 -I -m pytest -q -p no:cacheprovider tests/test_next_donor_selector_feasibility.py`.

## 1. Question

From known query inputs only, can a selector trained on all 101 training brain samples pick one training donor
whose measured axon, translated to the query's table soma, beats every copying control on BOTH the animal-mean
distance and the animal-mean F1 at 200 um, with copying as the fallback?

The output is always one measured training arbor, translated. No new geometry, no attachment, no deformation.
Any claim is about the combined rule (selector plus fallback), never about the selector alone.

## 2. What is already known, and why this is exploratory until held-out animals answer

* Codex's development selectors (`analysis/training_donor_selector*.py`, runs `run_v1` (ridge) and `run_hgb_v2`
  under `data/interim/training_donor_selector_development/`, `results.json`, 37 previously exposed training
  neurons, 32 animals), selector minus control:

  | control | ridge (run_v1): distance um / F1 | HGB (run_hgb_v2): distance um / F1 |
  |---|---|---|
  | nearest | -117.5 / +0.010 (both better) | -103.8 / +0.028 (both better) |
  | peptide_nearest | -206.9 / +0.048 (both better) | -193.3 / +0.066 (both better) |
  | peptide_medoid | -273.2 / +0.064 (both better) | -259.5 / +0.082 (both better) |
  | medoid | **+7.0** / +0.019 (lost on distance) | **+20.7** / +0.037 (lost on distance) |
  | kernel16 | -2.7 / +0.019 (both better) | **+11.0** / +0.037 (lost on distance) |
  | kernel_all | **+9.7** / +0.013 (lost on distance) | **+23.4** / +0.031 (lost on distance) |

  `beats_all_four_copying` is false in both runs.
* The general layer's frozen run (`reports/synarch_recon_results.md`) lost to KNN/ridge on all four tasks.
* The 548 um best-donor figure (local 16 neighbours, 37 neurons) is truth-informed: the best donor was picked by
  looking at the answer. It is a ceiling diagnostic, not evidence that permitted inputs can find that donor. No
  rule, setting or candidate set in this protocol is chosen from any oracle quantity. The protocol's own oracle
  ceiling is reported only after outer scoring, for description.
* The 101 training animals' geometry has been used for design by Codex (37 neurons), by the engine's development
  step and by the general layer. Every TRAIN number this protocol produces is therefore internal and exploratory,
  even though it is honestly nested. Only sealed forecasts on held-out animals can confirm (section 12).

## 3. Unit, cohort, pool

* Unit of independence: the animal (brain sample ID). Aggregation: neurons within animal, then equal animals.
* Training library: the 4,368 included training neurons of 101 animals (`data/interim/axon_engine/real/ingestion_log.json`,
  80 of 4,448 excluded by the plan's rules). Neurons per animal: 1 to 135, median 40.
* Donor pool of every method for a query: all training neurons of the current fit animals (TRAIN CV) or of all
  101 animals (held-out), never the query's animal. No sampled panel.

## 4. Permitted inputs, and exactly where they come from

**Query inputs (training or held-out), one loader only:** `known_query_inputs(parts)` in the protocol code
(implemented and tested now in `analysis/next_donor_selector_feasibility.py`). It reads exactly three files:

| field | source file | column |
|---|---|---|
| `id`, `animal`, `part` | `reports/axon_dataset_manifest.tsv` (pinned, engine `PLAN_MANIFEST_SHA`) | `neuron`, `animal`, `part` |
| `peptide`, `nucleus` | same | `source_label`, `soma_label` |
| `hemisphere`, `reconstruction_type` | `data/interim/axon_dataset/local_manifest.tsv` (from the portal index; animal, labels and part cross-checked against the plan manifest) | `hemisphere`, `reconstruction_type` |
| `soma_um` (table soma, mirrored-left frame) | `data/raw/jiao2025/All_soma_infor_HY_neurons.csv` (pinned, engine `SOMA_TABLE_SHA`) | `AP`, `DV`, `ML` only, plus `file`, `neurons`, `geno` for the identity and label checks; never `type`, `name` or `location` |

These are the sources and checks the engine's `Cohort` uses (`analysis/axon_engine.py`, `World.local_manifest`,
`World.soma_table`, `SOMA_COLUMNS`, the label and mirrored-frame checks); the loader re-implements them so that
protocol 2 does not import the engine. Rows whose neuron is not in the requested parts are skipped before any of
their fields is parsed. **The regional table is never opened** (no argument names it; a file named like it is
refused; a test spies on every `open` and asserts it). Output fields: `id, animal, soma_um, peptide, nucleus,
hemisphere, reconstruction_type`, nothing else. Rows handed to the referee carry its `known()` fields only
(`id, animal, soma_um, peptide, nucleus`). Never the query's measured arbor, endpoints, regional profile or any
outcome. Checked on the training rows (`known_inputs`): 4,448 rows; table soma equal to the ingestion log's
table soma for 4,368/4,368 included neurons; pins equal to the engine's constants; no held-out row parsed.

* Donors: their known fields and their measured TRAIN arbors (summaries: cable, branch count, axis RMS extent;
  centrality in the fit pool from the centred cache, fit animals only, section 7).
* Training labels: the metric between a fit-animal query and a fit-animal donor (section 5).
* Held-out forecasting: every held-out neuron declared in the plan manifest for that stage is forecast before any
  of that stage's truth is admitted; only the neurons the engine admits are scored, and the number forecast but not
  admitted is reported.

## 5. The native cache: what it means, and which labels may come from it (verified)

Definition (`analysis/axon_engine.py`, `PairCache`, `pair_cost`, read-only): `D[i,j]` = `distance_um`,
`F[i,j]` = `1 - f1_200` of the frozen `geometry_metrics(a.shifted(0), b.shifted(0))`, both arbors moved so their
own root (the TRACE soma) is at the origin, lower ID first.

The referee's label for a translated donor is `geometry_metrics(donor.shifted(table_soma_q), query_measured)`:
donor root on the query's TABLE soma, query arbor in absolute coordinates, rooted at its TRACE soma. The metric
is a symmetric mean of nearest-point distances, so moving both trees by one vector changes nothing, and the two
differ only by the query's table-to-trace offset `delta_q`; the distance can move by at most `|delta_q|`.

Checked by the command above (output section `cache_meaning`):

| comparison | where | result |
|---|---|---|
| cache vs its own formula (d110 metric, centred) | local37, 1,322 directed cross-animal pairs; real cache, 70 pairs | bit-identical: 1,322/1,322; 70/70 |
| cache vs the referee's metric module on centred trees | local37, 1,322 | bit-identical 1,322/1,322 |
| cache vs donor translated to the query's TRACE soma, absolute | local37, 1,322 | equal within 9.1e-13 um, 2.2e-16 F1 (rounding) |
| cache vs referee label (donor at query TABLE soma) | local37, 1,322 | distance max 1.33 um, mean 0.021 um; F1 changed in 33 pairs, max 0.0019 |
| same | real library, 40 random pairs | distance max 0.23 um; F1 changed in 1 pair, max 0.00016 |
| same | real library, 30 stress pairs (6 largest-offset queries, 5 nearest donors each) | distance max 1.69 um; F1 changed in 29/30, max 0.0057 |
| bound `|label - D| <= |delta_q|` | all 1,392 pairs above | 0 violations |
| per-query donor ranking, cache D vs label distance | local37, 37 queries | same best donor 37/37; minimum Spearman 1.0; regret 0 um |
| frozen attachment ladder (500, then 1,000 um, else translation) vs cache | local37, 60 sampled pairs | 59 fell back to translation (diff <= 0.34 um); the 1 attached pair differs by 22.9 um and 0.118 F1 |

The independent reviewer's own check (seed 777, 12 random and 8 large-offset pairs) found the same: D and F equal
the referee metric on centred trees exactly; the label differs by up to 2.53 um and 0.0006 F1 on large-offset
queries; 0 bound violations. Training offsets `|table - trace soma|`, 4,368 neurons: median 0.002 um, 95th
percentile 0.64 um, max 4.84 um; 162 above 1 um. The real cache (`data/interim/axon_engine/real/paircache/85f0f247bacfd142`,
`cache_sha256 b8c089ee...0bd0`) is read only after its COMPLETE.json exists; its D, F and index hashes and every
arbor hash against the ingestion log are verified.

Consequences, fixed now:
* **Fit labels may come from the cache**, as a declared soma-centred target (error against the referee label
  bounded by the query offset, at most 4.84 um in distance; at most 0.006 F1 observed). Every fit label read is
  checked for finiteness; the bound is asserted on every pair that is also recomputed (open decision 2: exact
  recompute instead).
* **Every SCORED number is recomputed exactly** as the referee scores it: donor translated (or, for the attached
  controls, attached) to the query's table soma, against the query's measured arbor, by the referee's metric
  module. This covers inner-CV, final-CV, outer and held-out scores, for the selector and for every control.
* **The medoid, peptide-medoid, C1, kernel, local-medoid and E1/E2 template choices read the cache** (through the
  fit-only door of section 7), because the referee's own medoid is defined on soma-centred distances.
* **Attachment outputs (C4, C5, E1, E2) can never come from the cache**; their arrays are made by the frozen
  `compact_attachment` and scored by recomputation.

## 6. Method (fixed now)

* **Fallback donor of a query (fixed now):** the referee `medoid` (`brm_geometry_audit.geometry_medoid`
  arithmetic: animal-balanced weights, argmin of `D @ w`) of the pool WITHOUT the query's animal. For an
  evaluation or held-out query this is the pool medoid; for a fit query it is the leave-own-animal-out medoid, so
  the referee rule "the query's animal never enters the donors" holds for fit labels too. Measured on TRAIN: 2 or
  3 distinct fallback donors among the fit queries of a fit call, 1 for its evaluation queries (`folds.fit_calls`).
  Reason for `medoid`: on the exposed development runs it was the strongest copy on distance (exposed-data
  reasoning, open decision 1).
* **Candidate set of a query:** the M nearest donors of other fit animals by table soma (ties by ID) **plus the
  query's fallback donor**. M is chosen by CV from {16, 64}. Without the addition the fallback would almost never be
  a trained pair: for M = 64 the fallback lies outside the 64 nearest for 3,297 to 3,541 of the 3,361 to 3,556 fit
  queries of an outer fold (`folds.fit_calls[*].M64_fallback_added_to_candidates`).
* **Fit pairs:** every (fit query, candidate) pair, the (fit query, fallback donor) pair included; both endpoints in
  the fit animals of that fit call (section 7).
* Learner, chosen by CV from three fixed configurations: two-head ridge (alpha 1 or 100, as Codex's
  `training_donor_selector.py`), and one histogram gradient boosting configuration (max_iter 200, learning_rate
  0.05, max_leaf_nodes 31, min_samples_leaf 50, random_state 20261111). Heads: `log1p(D/200)` and `F`.
  Six configurations in all (3 learners x 2 values of M).
* Features: Codex's pair features (query and donor soma, differences, absolute differences, soma distance, source
  and soma labels, same-label flags, known flags) plus query hemisphere and reconstruction type, and donor TRAIN
  summaries (cable, branches, axis RMS extent, centrality = animal-weighted mean cached D to the fit pool, plus a
  flag "is the fallback donor"). All vocabularies, scalers and centrality are fitted on the fit animals of that fit
  call only.
* Pair weights: 1 / count(query animal, donor animal), normalised to mean 1.
* Ranking within the candidate set: score = predicted log distance + predicted F loss; ties by donor ID.
* **Copying fallback rule (fixed now):** the selector's top donor is kept only if its predicted distance is lower
  AND its predicted F1 higher than the predicted values of the fallback donor for that query (both are trained
  pairs); otherwise the fallback donor is copied. Non-finite predictions, or a fit that cannot be made, fall back.
  Every row records which branch was taken (`fallback_choice`, tested).

## 7. Whole-animal cross-validation on TRAIN, and the final fit

* **Outer/inner (assessment):** 5 outer folds of whole animals (salted hash `next_donor_selector|20261111`,
  salt `outer`), 4 inner folds within each outer fit set (salt `inner|k`). Outer: 80 or 81 fit animals (3,361 to
  3,556 neurons), 20 or 21 evaluation animals (812 to 1,007 neurons). Inner: 60 or 61 fit animals, 20 or 21
  evaluation animals. Inner CV picks one configuration per outer fold by the animal-mean joint objective
  `log1p(distance/200) + (1 - F1)` of the combined rule (ties: fewer parameters, then ridge, then smaller M); it is
  refitted on the outer fit animals and forecasts the outer evaluation animals.
* **Final CV (the configuration of the final fit), fixed now:** 5 whole-animal folds over all 101 training animals
  (salt `final`): 80 or 81 fit animals (3,338 to 3,717 neurons), 20 or 21 evaluation animals (651 to 1,030
  neurons). Same six configurations, same objective, same tie rules. Its choice, and only its choice, is the
  configuration fitted on all 101 animals for the held-out forecasts. The outer folds' choices are reported, not
  used. If the final-CV choice differs from the most frequent outer choice, that is reported.
* **Both endpoints of every training pair are fit-only animals of that fit call**, never an outer or inner
  evaluation animal, never the same animal twice. Proven for all 31 fit calls (5 outer, 20 inner, 5 final-CV,
  final-all), label pairs and fallback pairs, by `check_fit_pairs` (0 violations), tested.
* **Fit-only cache door:** every cache read of a fit call (labels, centrality, medoids, fallback donors, kernels,
  local medoids, E1/E2 templates, peptide medoids) goes through one accessor, `FitCache`, built with that call's
  fit animals; it raises on any neuron outside them and counts entries served and refused. Shown for all 31 fit
  calls (a probe read of an evaluation neuron was refused in all 30 calls that have evaluation animals), tested.
  The run logs the read and refusal counts per fit call; any refusal during a run is a `RUN_FAILED`.
* All outer forecasts (selector and every control, for all 4,368 queries) are written and hash-sealed before the
  first outer score is computed; the query's measured arbor is loaded only by the scoring step.

## 8. Controls (same queries, same pool); the TRAIN gate and the held-out gate use the SAME set

| control | rule | TRAIN CV: made by | held-out: made by |
|---|---|---|---|
| `nearest`, `medoid`, `peptide_nearest`, `peptide_medoid` | referee arithmetic, full eligible pool, translated | protocol code on the fit pool (cache through `FitCache`) | the referee itself, inside its seal |
| C1 `C1_medoid_copy` | frozen `training_medoid` of the pool, translated | protocol code, same arithmetic via the cache | engine's committed forecast (read-only, hash-checked) |
| C4 `C4_uniform_attach_500`, C5 `C5_uniform_attach_1000` | the uniform branch of the frozen `predict_policy`: C1's donor attached by `compact_attachment` at 500 (C4) or 1,000 um (C5) to the query's table soma, else C1's translation | protocol code: exactly that branch (the posterior it computes does not change a uniform output); a pre-run test on a small TRAIN pool must show bit-identical arrays with the unmodified `predict_policy`, otherwise the unmodified `predict_policy` is used | engine's committed forecasts |
| E1, E2 | engine arithmetic: weighted medoid of the 16 nearest (cache D), weights `a_j/A` (E1) or `(a_j/A) k_j` (E2), attached by the ladder 500, then 1,000 um, else translated | protocol code | engine's committed forecasts |
| `kernel16`, `kernel_all` | Codex's fixed-kernel medoids, translated | protocol code | protocol code, sealed with protocol 2's forecasts |
| `local_medoid_M` (matched to the chosen M) | referee-arithmetic animal-balanced medoid of the M nearest donors of the chosen configuration (cache D), translated | protocol code | protocol code, sealed |
| `switch_nearest`, `switch_random` (matched to the switch) | on exactly the queries where the combined rule left the fallback: the soma-nearest candidate other than the fallback, or a seeded random candidate other than the fallback (`stable_seed(20261111, "switch_random", query)`); elsewhere the fallback | protocol code | protocol code, sealed |
| `permuted_labels` (negative control) | the same pipeline and chosen configuration, labels permuted among the fit pairs of each fit call (seed 20261111) | protocol code | protocol code, sealed |

* The held-out pool of the engine's controls is the whole training library, the same 4,368 neurons as protocol 2's
  final fit. Their committed forecasts are read only, hash-checked against the engine's commitment.
* Cost on TRAIN (measured, section 10): attachment calls are at most 48 ms each when the donor is near enough to
  attach (30 attached calls on real nearest donors), so the attached controls cost at most about 0.35 processor
  hours; C4 and C5 are therefore computed on TRAIN and are in the TRAIN gate.
* **Coverage is reported**: the fraction of queries where the combined rule left the fallback, overall and per
  animal (`matched_switch_controls`, tested), for TRAIN outer and held-out.

## 9. Measures, gate, decision rule and the reporting sentence

* Measures: `distance_um` and `f1_200` of the referee's metric; neuron then animal means.
* **Gate (both measures, every control of section 8, TRAIN and held-out alike):** the combined rule must have
  lower animal-mean distance AND higher animal-mean F1 than every control. No alpha, no correction: the rule is an
  intersection over all controls and both measures, so adding a control can only make it harder.
* **Reporting sentence (fixed now; `verdict_sentence`, tested):** "gate met and resolved against every control"
  when every paired animal-bootstrap 95% interval (selector minus control) excludes 0 in the favourable direction
  on both measures; "gate met, not resolved" when the point gate holds but some interval spans 0; "gate not met"
  otherwise. Success is decided by the point gate; the sentence always accompanies it. The null pass rate of the
  point gate is not known; the `permuted_labels` outcome is reported as one empirical null draw.
* TRAIN outcome: pass, or "no selector beyond copying on TRAIN". A TRAIN failure ends protocol 2: no held-out
  forecast is scored for a claim (open decision 4: seal them for the record anyway).
* Held-out outcome: success only if the gate holds on the held-out animals; otherwise "does not beat copying".
* Intervals: paired animal bootstrap 95% (10,000 draws, seed 20261110 as the referee), descriptive only.
* Subgroups (source label, soma label): descriptive; sparse below 5 animals.

## 10. Feasibility (TRAIN only; command in the header)

| item | number | source (`feasibility_all.json`) |
|---|---|---|
| training neurons / animals | 4,368 / 101 | `folds.counts` |
| directed cross-animal pairs (all) | 18,781,574 | `folds.all_directed_cross_animal_pairs` |
| unordered pairs in the cache | 9,537,528 (complete) | `folds.all_unordered_pairs_cache`, cache COMPLETE.json |
| full cross-animal fit pairs, outer fit | 11.1 to 12.4 million per fold | `folds.fit_calls[outer*].full_cross_animal_fit_pairs` |
| fit calls | 31 (5 outer, 20 inner, 5 final-CV, 1 final) | `folds.fit_calls` |
| label pairs any step can read (union, all 31 calls, fallback pairs included) | M=16: 202,670; M=64: 691,325 | `folds.union_of_candidate_label_pairs` |
| fallback pairs over all calls | 117,936 | `folds.fallback_pairs_all_calls` |
| frozen metric (d110, centred), one pair, real arbors | 7.0 ms mean (6.0 to 7.0 across three runs; load moves it about 10%) | `cache_meaning.real_library_sample` |
| referee metric (scoring), one pair | 5.4 ms mean | same, `ms_per_pair_referee_label_mean` |
| arbor load | 0.8 ms | same, `ms_per_arbor_load_mean` |
| attachment ladder, one call | 13.7 ms mean, 47.7 ms max when it attaches (30 of 50 calls) | same, `attachment_ladder` |
| exact recompute of the M=64 label union (open decision 2) | 1.35 processor h, 0.67 h on 2 workers | `timing.M64` |
| exact recompute of all directed pairs | 36.7 processor h (not needed) | `timing.full_directed_recompute_all_cross_animal_pairs` |
| all scoring (inner, final-CV, outer, 17 arrays per outer query) plus 6 attachment calls per query, upper bound | 205,296 pairs + 26,208 attachments: 0.66 processor h, 0.33 h on 2 workers | `timing.scoring` |
| model fitting (about 160 ridge/HGB fits), bootstrap | not measured | |

Minimal detectable paired difference (80% power, two-sided 0.05), from the SD of per-animal differences in Codex's
two development runs (32 animals, about 1 neuron each; planning only, `mdd`):

| animals | distance, median SD 416 um | F1, median SD 0.090 |
|---|---|---|
| 101 (TRAIN CV) | 116 um | 0.025 |
| 65 (validation + test) | 145 um | 0.031 |
| 40 (test) | 184 um | 0.040 |
| 25 (validation) | 233 um | 0.050 |

These SDs come from animals with about one neuron each, so they likely overstate the SD of animal means over about
40 neurons; the TRAIN CV gives the right SD before any held-out seal. Codex's selectors differed from the copies by
3 to 273 um: differences under about 100 um will not be resolved on 40 animals, hence the reporting sentence.

## 11. Run-once mechanics

1. Review; finalise this text; the lead commits protocol, code and tests. Pin: code hashes, the cache
   `cache_sha256`, the ingestion log, manifest, local manifest and soma table hashes, folds (outer, inner, final),
   configurations, seeds.
2. **`PROTOCOL_PREPARED`** (hash of everything pinned, and the commit) **is committed before the engine's
   `admit --stage validate`.** At the time of writing the engine has no `validate` directory
   (`data/interim/axon_engine/real/` holds `ingestion_log.json`, `library`, `paircache`, `referee_cache`; listing
   only). Only the TRAIN run and the held-out seal come later. If `PROTOCOL_PREPARED` is committed after
   `admit --stage validate`, the protocol declares that its design may have been informed by validation results;
   this affects only the validation animals, never the test animals.
3. `RUN_STARTED` -> outer forecasts sealed (`FORECASTS_COMPLETE`) -> outer scoring -> final CV -> `RUN_COMPLETE`.
   Never repeated; a failure writes `RUN_FAILED` and the reason; a new attempt needs a new, disclosed protocol
   version.
4. If TRAIN passes: fit the final-CV configuration on all 101 animals; forecast every declared held-out query of the
   stage from `known_query_inputs` only; write and hash the forecasts and protocol 2's own controls; the seal record
   also states the number of test neurons listed in the trace manifests by the engine's read-only listing (the
   check its freeze makes; it must be 0 for a test seal); the lead commits the hash (the seal time is the commit
   time).
5. Score through the independent referee (`analysis/prediction_independent_evaluator.py`, its own new bundle, the
   combined rule as candidate) only at the time section 12 allows, from the engine's admitted canonical truths
   (hash-checked). Protocol 2 never opens a trace itself.

## 12. Exposure: which animals can still be a genuinely unseen test, and when

Facts: protocol 2's forecasts need no held-out trace and no held-out regional row; only the known inputs of
section 4 and the training library. The current evaluation opens the 25 validation traces at
`admit --stage validate` and the 40 test traces after its freeze (stage B fetch, then `admit --stage test`). Test
soma positions and source labels were summarised earlier (engine protocol section 2): an input exposure this
protocol inherits; no outcome.

Options:
* **(a) Sealed before opening.** Protocol 2 is prepared before validation admission (section 11.2), and its
  held-out forecasts are sealed before the animals' traces are opened (for the test: before the stage B fetch).
  Order: `PROTOCOL_PREPARED` committed -> engine validation admitted and scored (protocol 2's TRAIN run at any time
  after `PROTOCOL_PREPARED`) -> engine freeze ->
  protocol 2 test forecasts sealed (commit) -> stage B fetch -> engine test admitted, sealed, scored once ->
  protocol 2 scored once.
  Why this cannot contaminate the current evaluation: its methods are frozen at its commitment and freeze; it never
  reads protocol 2's files; protocol 2 changes none of its files and opens no trace; protocol 2 reads test truth
  only after the referee's `SCORE_COMPLETE` of the current test, from the engine's admitted truths. A protocol 2
  commit cannot break the engine's freeze or test: the engine's `git_state` (axon_engine.py, lines 3907 to 3921,
  used at 4036 to 4038) checks only that its pinned files are committed and clean and never compares HEAD;
  protocol 2 touches no pinned path (engine files, engine protocol and tests, bridge test, frozen pins, freeze
  record, referee files and contract, general-layer code), and `PROTOCOL_PREPARED` asserts that its file list is
  disjoint from that list. Protocol 2 never runs more than 2 worker processes, including while the engine's
  predict or copies-check runs.
* **(b) Consumed re-use.** Any animal whose trace was opened before protocol 2's seal is a re-used outcome:
  reported as such, never "pristine", never for a claim.
* **(c) New animals from another dataset.** Real public resources of complete single mouse neurons registered to
  CCFv3 exist (MouseLight, Winnubst et al. 2019; SEU-ALLEN, Peng et al. 2021), but as far as I know they hold few
  or no hypothalamic neurons of these source labels, and their tracing pipelines differ. Not checked here (no
  network). Not a practical test set for this question now.

**Recommendation: (a) on the 40 test animals**, with `PROTOCOL_PREPARED` committed before validation admission.
The 25 validation animals are (a) only if the held-out forecasts are also sealed before `admit --stage validate`
(which needs the TRAIN run and its pass first, and would ask the engine lead to wait); otherwise they
are (b). Trade-off: 40 animals resolve only differences above roughly 180 um / 0.04 F1 (section 10), but they are
genuinely unseen for protocol 2 and the window is natural (the test traces are not fetched until after the
engine's freeze), so the current evaluation needs no delay.

**If the engine never admits test truths** (it stops before the test, or its test fails before admission):
protocol 2 is not scored under this version. Its sealed forecasts stay sealed. A later scoring needs a new,
disclosed protocol version naming the admission path (still the independent referee); the test animals count as
genuinely unseen for it only if no test trace was opened before protocol 2's seal, otherwise (b).

## 13. What a result would and would not mean

* Positive (held-out gate met): permitted inputs plus a TRAIN-fitted selector, with the medoid fallback, can pick
  translated donors that beat every copying, frozen, fixed-kernel and matched control on both measures, on these
  animals, at this tolerance, with the stated coverage. It does not mean the selected arbor is the right neuron,
  that soma position explains axon routes, or anything synaptic or causal; it does not generalise beyond this
  dataset's labelling and registration. "Not resolved" means the point gate held but the animals cannot separate
  it from chance against at least one control.
* Negative: this selector, with these features, cannot beat copying here. It does not show that no selector can;
  the TRAIN oracle ceiling (descriptive) says how much room was there.
* Any TRAIN-only result: internal and exploratory (section 2).

## 14. Decisions taken in revision 2, and what remains open

Taken (review round 1): query-input sources (section 4); TRAIN gate = held-out gate, C4/C5 included (section 8);
final CV fixed (section 7); fallback donor added to every candidate set and every fit query's labels (section 6);
coverage and matched switch controls (section 8); `local_medoid_M` matched control (section 8); fit-only cache door
(section 7); `PROTOCOL_PREPARED` before validation admission (section 11); `git_state` argument and worker cap
(section 12); no-test-truth rule (section 12); reporting sentence (section 9).

Open (lead or user):
1. Fallback donor: `medoid` (proposed) or another copy fixed now.
2. Fit labels: cache (proposed, bounded error) or exact recompute of the M=64 union (0.67 h on 2 workers).
3. Whether to add M = whole pool (about 11 million rows per outer fit: memory, not compute, is the limit).
4. If TRAIN fails: seal held-out forecasts for the record anyway, or stop.
5. Validation animals: whether the engine lead will wait for a protocol 2 seal before `admit --stage validate`
   (then validation is (a); otherwise (b)). Either way `PROTOCOL_PREPARED` precedes validation admission.
