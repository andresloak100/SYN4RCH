# Predicting a hidden 3D axon, against copying: the biological protocol

STATUS: DRAFT

This line stays until the validation traces exist and the user has accepted the text. On the real cohort the
engine refuses to commit any held-out forecast and to freeze while it is here; once it is removed the file's
SHA-256 is pinned by the validation forecasts and by the freeze. No independent biological test has been run.
Written 2026-11-10, revised 2026-11-11 after an internal review, the other agent's read-only review
(`reports/axon_engine_codex_review.md`) and the independent referee's contract. Implementation:
`analysis/axon_engine.py`, `analysis/axon_engine_cli.py`; tests `tests/test_axon_engine.py`.

Binding documents, none of them changed here: the locked cohort plan `reports/axon_dataset_plan.md` (commit
c96b167) with its amendment 1 `reports/axon_dataset_plan_amendment_1.md` (commit c659eac: no clipping) and its
amendment 2 `reports/axon_dataset_plan_amendment_2.md` (commit d43b654: the access route is the portal's per-neuron
address and the download belongs to the other agent, this engine downloads nothing; rule 8 for a neuron the
source does not serve; no test trace on disk before the freeze); the listing contract
`reports/axon_engine_ingestion_fields.md`; the referee's contract `reports/prediction_independent_evaluation_contract.md` with its code
`analysis/prediction_independent_evaluator.py`; the other agent's frozen prototype `analysis/brm_template.py`
with its freeze record.

## 1. Question

For a hypothalamic neuron of Jiao et al. 2025 whose axon is hidden, does a forecast made from what is known about
it lie closer to the measured axon than a copy of a measured neuron does?

A method **succeeds** only if, on the same admitted animals, its animal-mean geometric distance is strictly lower
**and** its animal-mean F1 is strictly higher than **every** copying control of section 6. A method that uses
conditional weights (source-label or regional) must in addition beat its matched control with uniform weights on
both measures. If the best method on validation is itself a copying control, the result is "no method beyond
copying". Success on synthetic data is not biological success. A failure is reported as a failure. A numerical win
on these animals is not a confirmatory population claim.

## 2. Cohort, unit, exposure

* 7,180 neurons, 166 animals, one source label per animal. Split by whole animal, fixed by a salted hash:
  training 101 animals / 4,448 neurons, validation 25 / 1,185, test 40 / 1,547. Manifest SHA-256
  `649f829a592cd3d13bcda5cffdbd049e759568d38bacb144f6a61220610a1ca4`.
* **Unit of independence.** An "animal" is the ID prefix, the portal's brain sample ID. The paper reports 166 mice
  and the data hold 166 sample IDs; no statement was found that no mouse contributed two samples. Every count
  below is therefore "animals (brain sample IDs)". A source-label effect cannot be separated from the animal
  and its batch (one label per animal; no batch record).
* **On disk when the engine was built and checked:** 37 traces under `data/raw/jiao2025/swc`, all of training
  animals (16 development neurons, 21 others), and no validation or test trace. These 37 are the fixed
  development snapshot `local37` (`--world local37`). The 16 development neurons stay development data.
  (A fetch of training traces into `data/raw/jiao2025/traces` was started by another process on 2026-11-11; the
  snapshot ignores that directory, and the report says what was there.)
* **Exposure, stated exactly.** The validation and test animals are new *geometry* animals: none of their
  traces has been opened by anyone in this project. They are not a pristine multimodal resource. The paper's
  49-region table for all 7,180 neurons, these animals included, has been on disk since 2026-11-10 and earlier
  work in this repository used it: the other agent's regional model was fitted and scored by 5-fold whole-animal
  cross-validation over all 166 animals, so regional outcomes of animals now in validation and test informed
  that model's design (its form, its grid of component counts) and its reported accuracy; outcome-only summaries
  of the table were printed (for example the count of all-zero profiles); and one raw profile was displayed.
  This engine loads the regional rows of training animals only (section 3) and refits that model class on them,
  but the model class itself is inherited with that exposure. It touches the candidates F1, F2, E4 and E5. The
  controls and the candidates E1, E2, E3 and E6 use no regional information.
* **Inputs of test animals that were looked at.** No outcome of a test animal was used. Their soma positions and
  source labels (permitted inputs, not outcomes) were summarised twice before the engine's settings were written
  down: in the plan's section 4 (distance from held-out somas to the nearest training soma) and by the first
  version of the engine's coverage diagnostic (radius of the 16 nearest training somas, same-label neighbours,
  reach of one donor), whose record is kept as
  `data/interim/axon_engine/real_local37/checks/coverage_with_test_inputs_2026-11-10.json`. The diagnostic now uses
  training animals (leave one animal out) and validation animals only, and no setting was changed because of
  either record.
* Apart from that, no test animal entered any choice made here: candidates, settings, preprocessing,
  normalisation statistics, donor selection, regional states. No registration is fitted. Eligibility follows the
  plan's mechanical rules only, never a forecast's performance.

## 3. Evidence boundary

**Permitted query inputs:** table soma position, source label, hemisphere, soma label, reconstruction type, and
what models fitted on training animals infer from them. The engine uses the table soma, the source label (and
the table's spelling of it for the regional model) and the identity. **Never inputs:** the query's own axon,
dendrite, or regional profile.

* **Known query soma.** The paper's soma table, mirrored-left frame, for every query. A query's trace is never
  opened for prediction, not even for its soma node. On the 37 inspected traces the table soma and the trace
  soma differ by at most 1.9 um (35 of 37 by less than 0.01 um). Every forecast is rooted exactly at the table soma.
* **Regional rows.** One loader. It reads the table line by line and drops every line of a validation or test
  animal before any value on it is parsed. On the real table: 4,448 rows loaded, 2,732 skipped unparsed.
* **Regional states.** The frozen joint model (`brm_real.fit_outer_fold`, unmodified, with its own choice of the
  number of components by inner whole-animal folds) fitted on training animals only: one state on all 101
  training animals for validation and test queries, and five out-of-fold states by animal (label-stratified,
  hash order) for training queries used in development. A fit on any non-training animal is refused. A state
  manifest records each state's ID, the hash of its sorted animal list and the code hashes. A state is used
  only when its stored parameters, its proof, the state manifest and the planned fitting and served sets agree;
  the bytes of the state and proof files are pinned by every forecast manifest and by the freeze, and the freeze
  fits the all-training state again and requires the same hash.
* **Donors.** Canonical arbors of training animals; for a training query, of other training animals.
* **Automatic checks, every forecast, every method:** the query's animal is neither in its regional state nor
  among its donors; no validation or test animal is in any state or donor set. A prediction log (query, state ID,
  hash of the donor animal set) is written into the manifest and re-checked before the commitment and again
  before any truth is opened.
* **Stage guard.** `develop` opens training traces only. `predict` (validation or test) opens no held-out trace
  and runs with the held-out traces absent from the machine. Held-out traces are opened first by `admit`, after
  the committed forecasts have been verified and bound to the method on disk (section 7).

## 4. Ingestion and the truth that is scored

Plan section 9 with amendments 1 and 2, applied mechanically to every trace, donors and truth alike; every
violated rule is logged; no rule consults a model, a forecast or an error.

1. Tree: node 1 is the soma (type 1, parent -1); every parent exists; one connected tree; no cycle.
2. Soma within 5 um of the table soma after the frame rule.
3. Soma hemisphere equals the index (for a mirrored-left file: soma on the left).
4. Node types within {0, 1, 2, 3, 4, 10, 11, 12}; axon coded as type 2 or type 0, not both (recorded).
5. `Axon_only` files have no dendrite node; `Axon_and_dendrite` files have at least one.
6. Axon cable of at least 1,000 um (the plan's axon rule, on the complete tree).
7. **No clipping** (amendment 1). The complete measured tree is ingested, used as a donor and scored. **The
   scoring domain is the complete measured tree**: truth, copying controls and candidate forecasts live in the
   same full coordinate domain, including cable outside the atlas box (AP 0 to 13,200, DV 0 to 8,000, ML 0 to
   11,400). The length of cable outside the box is logged per neuron and reported per part; it is never removed,
   moved or reconnected. The engine has no clipping function.
8. **Unavailable at the source** (amendment 2). A neuron the download listing records as `unavailable_at_source`
   is excluded, logged under rule 8 with the listing's reason and counted per part and per source label. It is
   decided by the source alone. A training neuron under rule 8 is not a donor; a held-out neuron under rule 8 is in
   the referee declaration's exclusions before anything is sealed. If a whole held-out brain sample is lost, or more
   than 5% of a held-out part, the tables and the report say in their first lines that the scored cohort is
   smaller than the locked one.

* **Frame.** CCFv3 micrometres, axes AP, DV, ML; midline ML 5,700. A true-hemisphere file of a right-hemisphere
  neuron is mirrored, ML' = 11,400 - ML, whole neuron; AP and DV are untouched. When a neuron has a
  true-hemisphere copy and a mirrored-left copy, the true-hemisphere copy is used.
* **Axon rule (the plan's).** Every node that is not dendrite (types 3, 4) and not the soma is axon; a marker
  node (10, 11, 12) takes the type of its nearest non-marker ancestor.
* **Connector geometry (fixed now, before any validation trace exists).** An ancestor of an axon node that is
  not axon itself (the dendritic stem of an axon that leaves a dendrite; a marker between the soma and the first
  axon node) is kept as connector geometry of the rooted tree that is scored; its length is logged per neuron and
  reported per part; rule 6 counts axon cable only. A marker inside an axon path is axon by inheritance. No
  connector occurs in the 37 traces.
* **Canonical arbor.** The frozen code's own `resample`: branch-to-branch polylines sampled every 50 um; every
  fork and tip is kept at its measured coordinates; between them the curve is a sampled one, so a sharp bend is
  cut by a chord (the canonical cable can be slightly shorter than the measured cable).
* **The one difference from the earlier runs.** Terminal marker nodes below an axon node are axon cable under
  the plan's rule; the frozen extractor drops them. It is at most 20.5 um per neuron on the 37 traces. Apart
  from these nodes the canonical arbors of the 16 development neurons are those of the earlier runs, bit for bit.
* **Cable outside the atlas box** stays where it was measured. In the 37 traces, 5 axons run past AP 13,200.
* A neuron that violates a rule is out for every method: not a donor, not scored. Forecasts are written for
  every declared query before any held-out trace is opened; exclusions are counted per part and per rule.
* **Listing.** The engine reads traces only through the listing files (`data/raw/jiao2025/swc/MANIFEST.tsv` for
  the 37 early traces, `data/raw/jiao2025/traces/MANIFEST.tsv` for the download; fields in
  `reports/axon_engine_ingestion_fields.md`). It never scans directories, never opens a file that is not listed,
  re-hashes every file on open and refuses a mismatch. Every neuron of a needed part must be listed, with a file or
  as unavailable; a listing that is still being written cannot be used for the whole cohort (`develop` on it is
  refused until every training neuron is listed), and the training listing and the held-out listing are pinned by
  hash at the forecasts and at the admission. When a neuron has more than one copy on disk the copies are ingested
  and compared; a disagreement is logged and reported, and the preferred copy is used.

**Frozen before any forecast, each with its test** (`tests/test_axon_engine.py`):

| item of the referee's contract | where | test |
|---|---|---|
| raw SWC adapter (7 numeric columns, integer IDs, one soma root, parents, connectivity) | `read_swc`, `tree_structure` | `test_every_inclusion_rule_excludes_and_is_logged` (cases missing_parent, two_roots, cycle, root_not_soma, not_a_table) |
| type-0 / type-2 interpretation | `extract_axon` | `test_axon_rule_gives_the_same_arbor_for_type2_and_type0_coding_and_matches_the_frozen_extractor`; real: `test_real_type2_traces_give_the_frozen_canonical_arbor_bit_for_bit` |
| marker inheritance | `extract_axon` | `test_a_marker_inherits_the_type_of_its_nearest_non_marker_ancestor` |
| connector geometry, dendritic stem | `extract_axon` | `test_an_axon_that_leaves_a_dendrite_keeps_the_dendritic_stem_as_logged_connector_geometry` |
| hemisphere transform | `ingest_trace` | `test_true_hemisphere_copy_mirrored_equals_the_mirrored_left_copy`; real: `test_real_portal_copy_mirrored_equals_zenodo_copy` |
| atlas frame (midline 5,700; mirror 11,400; ML only; box) | constants | `test_atlas_frame_constants_and_the_mirror_changes_ml_only` |
| resampling (edges at most 50 um, forks and tips kept) | `canonical_arbor` | `test_resampling_keeps_forks_and_tips_and_no_edge_exceeds_50_um` |
| clipping: none | `ingest_trace` | the two `test_no_clipping_...` tests, `test_the_reviews_two_clipping_counterexamples_return_exactly_the_measured_cable`, `test_forecasts_and_truth_share_the_full_coordinate_domain` |
| mechanical QC rules 1 to 6 | `ingest_trace` | `test_every_inclusion_rule_excludes_and_is_logged`, `test_axon_only_without_dendrite_is_included` |
| rule 8, unavailable at the source, and its reporting duty | `TraceStore`, `rule8_report` | `test_a_neuron_the_source_does_not_serve_is_listed_excluded_under_rule_8_and_counted`, `test_rule_8_reporting_duty_a_whole_held_out_brain_sample_or_more_than_five_percent` |
| listing: final before use; second copies compared | `TraceStore.require`, `build_library` | `test_development_on_the_whole_real_cohort_needs_the_final_training_listing`, `test_a_second_copy_of_a_training_trace_is_checked_against_the_chosen_copy` |
| no test trace on disk before the freeze | `freeze` | `test_on_the_real_cohort_the_freeze_refuses_while_a_test_trace_is_listed_on_disk` |
| QC without model errors | `admit` | `test_admission_reads_no_forecast_and_fixes_the_exclusions_before_the_referee_declaration` |
| hidden-outcome invariance | canary | `test_hidden_outcome_canary_every_forecast_is_bit_identical`, `test_forecasts_need_neither_held_out_traces_nor_held_out_regional_rows`, `test_development_forecast_ignores_the_query_animals_own_traces_and_rows` |
| planted leaks | checks | the three `test_planted_leak_...` tests, `test_the_manifest_audit_catches_a_planted_donor_from_the_query_animal` |

## 5. Measures

Exactly the frozen `geometry_metrics`: `distance_um`, the symmetric arc-weighted nearest-sample distance between
forecast and truth, and `f1_200`, the F1 of cable within 200 um. Scored per neuron, averaged within the animal,
then over animals. Lower distance and higher F1 are better. The referee's metric module gives the same numbers
bit for bit (checked on every validation score and on every pair of the local library).

## 6. Methods: seven copying controls, two frozen candidates, six engine candidates

All are deterministic, take only the permitted inputs, and are computed by the engine for every declared query
of every stage. Fixed settings; the SHA-256 of the configuration block is pinned by every manifest and by the freeze.

### Donor pool of every method (no sampled panel anywhere)

The **training library**: every training neuron whose trace passes section 4 (for a training query in
development: minus the query's own animal).

| methods | pool |
|---|---|
| `nearest`, `peptide_nearest` | the whole library (the second: its neurons with the query's source label) |
| `medoid`, `peptide_medoid` | the whole library, all pairs (the second: its neurons with the query's source label) |
| C1, C4, C5, F1, F2 (frozen prototype) | the whole library, handed to the unmodified `predict_policy` |
| E1 to E6 | the 16 nearest neurons of the whole library by table-soma distance |

**All-pairs cache.** The soma-centred `distance_um` and `1 - f1_200` between every two library arbors, by the
frozen `geometry_metrics` with the lower neuron ID first (exactly what the frozen `pair_costs` computes), stored
as dense matrices with the ordered ID list, the arbor content hashes, the metric code hashes and its own SHA-256.
It holds training arbors only, is built in resumable blocks, and is pinned by every manifest and by the freeze.
A fixed sample of pairs is recomputed by the frozen `pair_costs` and by the referee's own loader and metric and
must agree bit for bit. The frozen prototype reads its pair costs and its medoid from this cache through a memo of
its own two functions; on small pools trees and records are bit-identical with and without it.

### Copying controls

| | control | how it is made |
|---|---|---|
| `nearest` (the user's C2) | the training neuron with the nearest table soma, translated to the query soma | referee contract |
| `medoid` (C1 by the referee's arithmetic) | the training neuron minimising the animal-balanced mean soma-centred geometric distance, translated | referee contract |
| `peptide_nearest` (C3) | `nearest` among training neurons with the query's source label; falls back to `nearest` when the label is missing or has no donor (recorded) | referee contract |
| `peptide_medoid` | `medoid` among training neurons with the query's source label; falls back to `medoid` | referee contract |
| `C1_medoid_copy` (C1) | the frozen prototype's `baseline`: its own medoid function on the whole library, translated | unmodified `predict_policy` |
| `C4_uniform_attach_500` (C4) | frozen `uniform_anchor_500` | unmodified `predict_policy` |
| `C5_uniform_attach_1000` (C5) | frozen `uniform_anchor_1000` | unmodified `predict_policy` |

* **The four official copies are the referee's.** For a held-out stage the engine calls the referee's own
  `copying_baselines(training trees, known query metadata)` and writes ITS arrays as the four controls. The engine
  also computes them by itself (its own code on the pinned cache, the contract's arithmetic); the two must give the
  same donors and arrays, and a difference is counted and the referee's forecast kept. When the pool is too large
  for an inline call (more than 250,000 pairs: the real library has 9,890,128), the forecasts carry the engine's
  computation and `copies-check` must then show exact equality with `copying_baselines` on the whole stage
  **before any truth is opened**. At the test the referee recreates the four copies once more inside its seal.
  In development (leave one animal out, one pool per animal) the referee's function is called for the first
  animals within a budget of 80,000 pairs; today that is every animal.
* `C1_medoid_copy` and `medoid` are the same rule in two arithmetics (the frozen function sums left to right, the
  referee uses a matrix product). They give the same neuron except at floating-point ties. Both are controls; the
  frozen prototype is not changed to force agreement.
* C1, C4 and C5 consult the regional state for membership checks only: their arrays do not depend on its numbers
  (tested), so they declare no upstream fit.

### Candidates

**Frozen candidates.** F1 = `conditional_anchor_500`, F2 = `conditional_anchor_1000`, through the unmodified
`brm_template.predict_policy` with the whole library and the query's regional state.

**Engine candidates.** Neighbourhood N(q): the K = 16 training donors with the nearest table somas (other
animals; ties by neuron ID). d_j: soma distance. a_j = 1 / (donors of j's animal in N); A = animals in N.
Kernel k_j = exp(-d_j^2 / 2h^2), h = max(100 um, distance to the 8th nearest donor). Label factor l_j = 1 for
the query's source label, 0.25 otherwise. r_j = the frozen conditional weights (`conditional_weights`) of the
neighbourhood donors under the query's regional state. D = the cached soma-centred pair distance.

| | weights w_j proportional to | template | conditional | fitted state | must also beat |
|---|---|---|---|---|---|
| E1 | a_j / A | weighted medoid: argmin_i sum_j w_j D(i, j) | no | no | |
| E2 | (a_j / A) k_j | weighted medoid | no | no | |
| E3 | (a_j / A) k_j l_j | weighted medoid | label | no | E1 and E2 |
| E4 | r_j k_j | weighted medoid | regional | yes | E1 and E2 |
| E5 | r_j k_j l_j | weighted medoid | label and regional | yes | E1 and E2 |
| E6 | as E2 | the E2 template, pruned: a vertex stays if at least 0.25 of the other neighbours' weight has cable within 200 um of it, or if it leads to such a vertex | no | no | |

The template is attached to the query soma by the frozen `compact_attachment` at 500 um; if unsupported, at
1,000 um; if unsupported, the template is translated. For E6 the neighbours are attached the same way before
support is counted, no vertex is moved, and the unpruned template is returned when less than 1,000 um of cable
would remain. Every fallback is recorded with its reason. E6 is the only method whose output is not a whole
measured arbor. F1 is matched to C4, F2 to C5.

**Candidate family (declared when the validation forecasts are made, pinned by the freeze).** `all`: the eight
candidates. `no_fitted_state`: E1, E2, E3 and E6 only, the candidates that use no regional state and so carry none
of the exposure of section 2. All eight are always forecast, scored and reported; the gate may select only a
member of the declared family. The default is `all`; choosing the other is the user's decision before validation.

**One permitted change before the protocol is finalised, on training animals only.** K may be reset once to 8 or
32 by the development pass on the complete training library (`develop`, leave one animal out), choosing the K with
the lowest animal-mean distance of E2 among those whose F1 is not below that of K = 16. The choice is written here
before this file loses its DRAFT line. Nothing else is tuned, and nothing at all after any held-out forecast is committed.

## 7. Validation, freeze, test

Order of the steps of a held-out stage, each run once:

1. `predict`: forecasts for every declared query by all fifteen methods, then the manifest and its SHA-256 (the
   **commitment**). Opens no held-out trace. Resumable until the commitment exists (the frozen prototype on the
   whole library is slow); a partial run made by another method is discarded.
2. `copies-check`: only when the pool was too large for the inline call (section 6).
3. `admit`: verifies that the committed forecasts are complete, unchanged and **were made by the method on disk
   now** (engine and frozen code, configuration, this protocol, cohort tables, training traces and library,
   all-pairs cache, regional state and proof bytes); only then writes the opened marker and opens the held-out
   traces, applies rules 1 to 6 and rule 8 and writes the exclusion log and the canonical truth with hashes. It
   reads no forecast. An admission that died may be completed on the same commitment and the same trace files; nothing else
   is ever resumed.
4. Scoring, on the admitted queries only.

**Validation** (25 animals) **is scored by this engine**, because the referee takes one frozen candidate and scores
the split it calls test, while validation must score the whole family to select. The engine's scorer is the frozen
metric; every validation score is also computed by the referee's metric module and must agree bit for bit; the
four official copies are the referee's (section 6). A started validation scoring is never repeated or resumed.

* **Frozen gate.** The unmodified `brm_template.select_policy` on the five frozen policies.
* **Engine gate.** A candidate of the declared family is eligible when it beats every one of the seven controls
  on both measures on the whole admitted validation set, and, if conditional, its matched controls as well (F1:
  C4; F2: C5; E3, E4, E5: E1 and E2). A frozen candidate is admitted only if the frozen gate itself selected that
  policy. Among eligible candidates: lowest distance, then highest F1, then the order F1, F2, E1 ... E6.
* **Nothing eligible:** the verdict is **"no method beyond copying"**, and that is the result. The gate then also
  records a *descriptive reference*, the family member with the lowest validation distance (then highest F1, then
  the fixed order), only so that the referee's single candidate slot can be filled if the test is run to report
  the controls; no claim can come from it, whatever the test shows.

**Freeze.** Written once, after validation and before any test step. On the real cohort it refuses while any test
neuron is listed on disk (amendment 2: the test traces are fetched only after the freeze). It verifies again that
the validated forecasts belong to the method on disk, fits the all-training regional state again and requires the stored hash,
and pins: engine and frozen code hashes, this protocol, the cohort tables, the configuration, the training trace
listing and library fingerprint, the all-pairs cache, the regional state (hash and file bytes), both gates, the
family, the selected candidate (or the descriptive reference) and the validation admission. It records the
referee's hashes and whether the referee on disk can seal the planned test bundle, how many test traces were
already on disk, and which pinned dependencies the commit does not hold. The real freeze also needs engine,
protocol and tests committed and unmodified.

**Test** (40 animals), **through the independent referee**:

1. `predict --stage test`: as above; it needs no test trace, so the test traces can be fetched after it.
2. `admit --stage test`: as above. Before it, the engine checks with the referee's own functions on a fixture that
   the referee on disk can seal the planned bundle; if not, it stops before any test truth is opened.
3. `seal`: builds the referee's bundle. `freeze_declaration` receives: the cohort after admission (training
   library, admitted validation queries, admitted test queries; known metadata only: id, animal, table soma,
   source label as `peptide`, soma label as `nucleus`), every other cohort neuron in `excluded_ids`, the animals the
   plan keeps out of validation and test in `test_excluded_animals`, the training trees with path and SHA-256, the
   roles, the geometry contract, the method artifacts (engine, frozen sources, tests, this protocol, the freeze,
   both gates, state and proof files, cache record, the engine's commitment and admission records) and the
   referee's evaluation settings (seed 20261110, 10,000 bootstrap samples, soma tolerance 5 um). Then
   `copying_baselines`, then every committed forecast file of the admitted queries copied byte for byte, then
   `seal_predictions`, which recreates the four copies by itself.
4. `score --stage test`: the referee's `score_once`, once. The verdict and the tables are read from its result.

Roles in the referee's declaration: `candidate` = the frozen selection (or the descriptive reference);
`copy_baselines` = `nearest`, `medoid`, `peptide_nearest`, `peptide_medoid`, `C1_medoid_copy`;
`matched_controls` = C4 and C5 always, plus the matched controls of a conditional candidate. Forecast rows carry
`donor_ids` (the one copied neuron; the 16 neighbours for E1 to E6; the selected template for the frozen attached
policies, whose whole pool is recorded once by hash) and `upstream_fit_animals` (the training animals of the
regional state for F1, F2, E4, E5; an empty list otherwise). The same list is in the declaration and the manifest
when the bundle holds such a method. The candidates that were not selected stay outside the referee's bundle;
after the referee has finished the engine scores their committed forecasts for the record only (no claim, no selection).

**Verdict.** Success: the selected candidate is better on both animal-mean measures than every copy baseline and
every matched control of the declaration, on the admitted test animals. Better than every copying control but not
than its own matched control: "rule not met", nothing is claimed for the conditional weights. Otherwise "does not
beat copying". Nothing selected on validation: "no method beyond copying". The sentence states the paired animal
differences against the strongest control on each measure with their 95% bootstrap intervals; intervals are
descriptive and do not enter the rule. The earlier figures on the 16 development neurons (922 um, F1 0.198) are
history, never a threshold.

## 8. Subgroups (fixed now)

* Per animal, always: every animal with its number of neurons and both measures for the candidate and each control.
* By source label, and by soma label (the manifest's fine label). `HY` and fibre-tract labels (lower-case
  initial) are shown as their own rows, marked unassigned; an empty label is shown as `unknown`.
* Every row shows neurons and animals (brain sample IDs). **Every subgroup row is descriptive and carries no
  population claim**; a group with fewer than 5 animals is in addition marked sparse. From the plan: test labels
  with 5 or more animals are Adcyap1 (6), Sst (5), Penk (5).

## 9. What each forecast leaves behind

For every query and every method: an arbor NPZ written by the frozen `save_arbor` (`xyz_um` N x 3 float,
mirrored-left frame, um; `parent` N int, root first with -1, parents before children; `types` N int), rooted
exactly at the known soma, edges at most 50 um, at `predictions/<method>/<query>.npz`. One manifest per stage
(JSON) with, per forecast: query ID, animal, full ID, part; the known soma used; source label and soma label;
method name and role; requested and executed policy, fallback and reason; the selected donor and the donor ID
list; the SHA-256 of the donor pool (lists stored once); regional state ID and hash and the hash of its animal
list; NPZ path, NPZ SHA-256, a content hash of the three arrays, node count. Once per manifest: the method pins
(section 7), the training ID list and its SHA-256, the library fingerprint, the cohort manifest hash, engine,
frozen and referee code hashes, the configuration, the regional states, the animals behind every upstream fit,
where the four official copies came from, the prediction log. The manifest lists every declared query times every
method; its SHA-256 is the commitment. The referee's bundle for the test is in `test/referee/`.

## 10. Order of work and commands

```
python3 -I analysis/axon_engine_cli.py status
python3 -I analysis/axon_engine_cli.py ingest --stage validate      # training library; needs all training traces
python3 -I analysis/axon_engine_cli.py cache --minutes 14           # all-pairs cache; repeat until complete
python3 -I analysis/axon_engine_cli.py states
python3 -I analysis/axon_engine_cli.py develop --per-animal 10 --minutes 14   # training only; the one permitted choice of K
#   finalise this protocol (remove the DRAFT line); commit engine + protocol + tests
python3 -I analysis/axon_engine_cli.py predict --stage validate --minutes 14  # repeat until committed; needs no validation trace
python3 -I analysis/axon_engine_cli.py copies-check --stage validate          # the referee's copying_baselines on the whole pool
python3 -I analysis/axon_engine_cli.py admit --stage validate       # opens validation traces: inclusion rules only
python3 -I analysis/axon_engine_cli.py score --stage validate       # once; both gates
python3 -I analysis/axon_engine_cli.py freeze
python3 -I analysis/axon_engine_cli.py predict --stage test --minutes 14      # before the test traces are fetched
python3 -I analysis/axon_engine_cli.py copies-check --stage test
#   fetch the test traces (plan stage B)
python3 -I analysis/axon_engine_cli.py admit --stage test
python3 -I analysis/axon_engine_cli.py seal                         # the referee's declaration and seal
python3 -I analysis/axon_engine_cli.py score --stage test           # the referee's score_once, once
```

`score` runs `admit` (and, for the test, `seal`) itself when they have not been run; the single steps exist so
that the commitment and the seal can be checked in between. `--world local37` runs any command on the fixed
development snapshot (outputs in `data/interim/axon_engine/real_local37`); without it the commands use every
trace listed under `data/raw/jiao2025/swc` and `data/raw/jiao2025/traces` (outputs in `data/interim/axon_engine/real`).

Traces are expected under `data/raw/jiao2025/traces/{train,validation,test}/<animal>/<neuron>.swc`, listed in
`data/raw/jiao2025/traces/MANIFEST.tsv` (columns file, neuron, bytes, sha256, source, frame; frame ending in
`true_hemisphere` or `mirrored_left`, or `unavailable_at_source` with no file). A trace stored under another part
than the locked split says is refused. The download is the other agent's (`analysis/axon_trace_download.py`,
stage A = training and validation; stage B = test, after the freeze); this engine starts none.

**Cost at full size** (measured on the 37 local arbors and projected; `reports/axon_engine_report.md` has the
figures): the all-pairs cache and the frozen prototype's forecasts take hours of processor time and are resumable;
the referee's `copying_baselines` recomputes all 9,890,128 pairs in one process each time it is called (once per
`copies-check`, once inside its seal) and cannot be resumed or fed a cache. That part is the referee's own; its
contract says full-cohort runtime has not been demonstrated.

## 11. What a result could and could not say

* The claim, if any, is about sampled-curve geometry at 50 um on these animals: not synapses, not circuits,
  not gene expression, not causes.
* Every method here outputs a measured arbor, or a pruned one. None grows new branches. A success would show
  that the choice of which measured arbor to use (or what to keep of it) carries information beyond copying.
* 40 test animals (brain sample IDs); subgroup tables are descriptive.
* Labels, animals and batches are confounded.
* Regional-conditional candidates inherit the exposure of section 2.
