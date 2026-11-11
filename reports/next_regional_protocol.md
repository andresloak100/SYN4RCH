# Next protocol 1: regional axon profile with calibrated uncertainty, against strong simple methods

STATUS: DRAFT, revision 2 (after the first independent review). Nothing here is frozen. No forecast exists.
No validation or test row, trace or result was read to write it.

Written 2026-11-11. Prefix `next_regional`. Files: this protocol; `analysis/next_regional_core.py` (forecast
format, the one declared scorer, the controls, the decision rule); `analysis/next_regional_feasibility.py`
(training-only feasibility); `analysis/next_regional_run.py` (query loader, select, predict, score);
`tests/test_next_regional_feasibility.py`, `tests/test_next_regional_run.py`; outputs in
`data/interim/next_regional/`. It changes nothing in the current evaluation (`reports/axon_engine_protocol.md`,
E1 to E7), in the general layer (`analysis/synarch_*`), in the referee or in any other file.

## 0. In plain words

* **Question.** For a neuron never seen, from its known inputs only, can a method give a **better predictive
  distribution** of how much axon it sends to each of the 49 published regions than strong simple methods,
  **without worse point predictions**, and with **honest uncertainty**?
* **The method to beat is KNN** (the 40 nearest training neurons by soma position). On the training animals it
  beat the general layer's learned model, ridge and the means (exploratory, section 13).
* **Decision.** One candidate, chosen on training animals only, must have a CRPS at least 2% below every
  distribution control, a log MSE no more than 1% above the best point control, and calibration that holds
  for the zero part and the non-zero part separately, overall and region by region.
* **Who can test it.** No Jiao held-out animal is pristine for regional outcomes (section 12). Recommendation:
  seal forecasts for the 40 test animals before their traces are fetched, score only after the current
  evaluation has finished, label the result "held-out, previously modelled in aggregate; not pristine".
* **Likely outcome, stated now.** On training, simple residual-based distributions (the design of candidate N1)
  were worse than KNN's own neighbour distribution and badly calibrated for presence (section 13). The most
  likely result of this protocol is "not run" or "no gain". That is an acceptable, informative outcome.

## 1. Question

For a hypothalamic neuron of Jiao et al. 2025 whose brain sample was never used to fit anything, predict from its
permitted inputs its 49-region axon-length profile on the scale `log10(1 + published length)` as a predictive
distribution per region (point, 199 quantiles, probability of non-zero length). Is the selected candidate's
distribution better (CRPS) than every declared control's, its point no worse than the best control's, and its
uncertainty calibrated?

Declared secondary: presence or absence per region, scored by the Brier score; and the descriptive label "also
better points" (log MSE at least 2% below every point control).

## 2. Unit, cohort

* **Unit of independence:** the brain sample (ID prefix; "animal" below). Every summary is a mean over regions
  per neuron, then over neurons within a sample, then over samples (each sample counts once).
* **Fitting:** the 101 training samples (4,448 neurons) of the locked split, through the general layer's sealed
  loader (`synarch_recon.load_training_table`): validation and test lines are dropped unparsed (2,732 lines).
* **Scoring set (recommended, section 12):** the 40 test samples (1,547 neurons; `reports/axon_dataset_plan.md`).
  Validation samples are not used.
* The published table holds a row for every cohort neuron, so plan rule 8 (trace not served) does not remove a
  regional row. Any missing row is listed by ID at scoring and the neuron is dropped for every method alike.

## 3. Inputs (decided: soma position and source label only)

* **Used by every method, candidates and controls alike:** table soma position (AP, DV, ML; mirrored-left frame)
  and source label.
* **Not used (decided in this revision, finding 10):** soma label, hemisphere, reconstruction type. Hemisphere is
  already folded into the mirrored frame; the reconstruction type is the soma table's `type` column, which the
  general layer's sealed loader never reads (its allow-list is `file, neurons, geno, location, AP, DV, ML`), and
  this protocol's query loader keeps exactly that allow-list. If the lead wants either input, it must be added
  before the freeze to candidates and controls alike, and the loaders' allow-lists change with it.
* **Training information:** regional rows, soma rows and labels of training samples only. No trace, no geometry.
* **Never inputs:** the query's regional profile, trace, endpoints or geometry, any value inferred from them,
  any validation or test outcome, and any file derived from them (in particular `data/interim/brm_real/truth.npz`,
  `data/interim/brm_real/predictions/`, the current evaluation's admitted truth and scores).
* **Mechanics (finding 5, now written):** `next_regional_run.load_query_inputs(world, part)` reads the manifest
  and the allow-listed soma columns of one part; rows of every other part are dropped before any cell is parsed;
  it never opens the regional table. `predict` reads training outcomes only through the sealed loader.
  Hidden-outcome canary (`tests/test_next_regional_run.py::test_hidden_outcome_canary`): on a synthetic world in
  the real formats, the forecast files are **bit-identical** when the test samples' regional rows are perturbed,
  deleted or replaced by unparseable poison. The query loader also works with the regional file absent. On real
  data, `check-loader` (training rows only) confirms the query loader returns exactly the sealed loader's
  neurons, samples, labels and soma positions (`data/interim/next_regional/loader_check.json`: 4,448 neurons,
  101 samples, 2,732 other rows dropped unparsed, `agrees_with_sealed_loader: true`).

## 4. Outputs (the sealed forecast format) and the declared scorer

Per query neuron, method and region (49): `point` (log scale), `q` = 199 quantiles at levels (i - 0.5) / 199,
`p_nonzero`. Stored as `forecasts/<method>.npz` (point and p_nonzero float64, q float32) with `manifest.json`
(query IDs, samples, inputs used, settings, training-sample hash, loader audit, every forecast file's SHA-256 and
every code file's SHA-256; no timestamp). The manifest's SHA-256 is the **commitment**.

**One scorer for everything (finding 7).** `next_regional_core.score_forecast` is used by the feasibility run and
by the final scorer, for every method alike:
* `Q(level)`: linear interpolation of the 199 quantiles in the level, flat beyond the end levels.
* **CRPS** = 2 x mean pinball loss over the 199 levels (quantile CRPS). For empirical samples (KNN neighbours, a
  label's pool, point + residuals) the 199 quantiles are `numpy.quantile(..., method="linear")`. On training,
  KNN's quantile CRPS is 0.4345 against 0.435 for its exact ensemble CRPS in revision 1.
* **Marginal 90% interval** `[Q(0.05), Q(0.95)]`; coverage = truth inside.
* **Conditional non-zero 90% interval** `[Q(F0 + 0.05 p), Q(F0 + 0.95 p)]` with p = p_nonzero, F0 = 1 - p: the
  central 90% of the forecast's non-zero part. If p = 0 the interval is empty and a non-zero truth is not covered.
  A calibrated forecast covers about 90% of non-zero truths with it, whatever the zero share
  (`test_conditional_nonzero_interval_is_calibrated_for_a_mixture`), while its marginal interval can cover far
  fewer of them. This is why the review's proposed per-region floor on *marginal* coverage of non-zero truths
  was replaced by one on *conditional* coverage (section 8, finding 1).
* Also: log MSE, log MAE, width, interval score (alpha 0.1), presence Brier with p_nonzero, presence gap per
  region (mean p_nonzero minus observed non-zero share).

**Viewer (finding 6, decided).** On a pass the scorer writes `viewer_regional_input.json` in the input contract
of `analysis/prediction_viewer_export.py export_regional`, with `evidence_level` =
`held_out_previously_modelled_in_aggregate` and the exposure label of section 12 in `benchmark`. The current
viewer does not know that level and refuses the export (tested:
`test_viewer_export_refuses_the_new_evidence_level`). **A pass is not shown in the viewer until the viewer's owner
adds that level with a label of its own** (for example "held-out, previously modelled in aggregate; not
pristine"). It is never exported as `independent_test`, whose heat-map label would say "independent test". When
shown, it is a labelled heat map on CCF regions (all 49 names checked against `structure_graph_1.json`, none
missing), never a path.

## 5. Methods

### Controls (fitted per fold on training samples only; settings by whole-sample CV inside the fitting set)

| control | point | predictive distribution | settings |
|---|---|---|---|
| MEAN | training mean | all training neurons' values | weighting: sample or neuron |
| IDENTITY | source label's training mean (unseen label: MEAN) | the label's training neurons | weighting |
| KNN | mean of the k nearest training neurons by soma | their values | k in {5, 10, 20, 40, 80}; pool all or same label; raw um or z |
| RIDGE | ridge on standardised soma terms + one-hot label, clipped at 0 | point + the region's **out-of-fold** residuals (3 whole-sample folds of the fitting set), clipped at 0 | lambda in {1e-3, 1e-2, 0.1, 1, 10, 100, 1000}; degree 1 or 2 |
| KNN-CONF | KNN point (KNN's log-MSE setting) | point + out-of-fold KNN residuals of the same region and the same bin of the neighbours' non-zero share ([0, 1/3), [1/3, 2/3), [2/3, 1]; under 20 residuals: the region's pooled residuals), clipped at 0 | none of its own |

* **Point controls:** MEAN, IDENTITY, KNN, RIDGE, each tuned by inner log MSE.
* **Distribution controls:** MEAN, IDENTITY, KNN@crps and RIDGE@crps (tuned by inner CRPS), and KNN-CONF.
* **Finding 3 (RIDGE):** its distribution now uses out-of-fold residuals, the same construction as N1, and the
  lambda grid reaches 1000. CRPS tuning now picks lambda 1 or 10 (interior of the grid) instead of the old edge.
* **Finding 4 (attribution):** KNN-CONF is "conformalised KNN" (N1 at w = 1) and is a distribution control the
  candidate must beat, so a pass cannot be just KNN with conformal intervals.

### Candidates (fixed now; none has been run)

* **N1, KNN-ridge blend with Mondrian out-of-fold residuals.** Point = w KNN + (1 - w) RIDGE (the controls'
  log-MSE settings of the same fold), w in {0.25, 0.5, 0.75} by inner log MSE. Distribution: the KNN-CONF
  construction applied to the blend's out-of-fold residuals. **Expectation stated before any run:** the same
  construction applied to KNN (KNN-CONF) was 10.7% *worse* than KNN@crps in CRPS on training and overstated
  presence in all 49 regions (section 13). N1 will very probably fail; it stays declared unless the lead decides
  otherwise (open decision 2).
* **N2, hurdle KNN.** p_nonzero = (non-zero count among the k_p nearest + 1 x training presence rate) / (k_p + 1),
  k_p in {10, 20, 40, 80}; level = values of the k_l nearest training neurons with non-zero length in that region,
  k_l in {5, 10, 20, 40}; distribution = 0 with probability 1 - p, else the level values; point = p x mean level.
  (k_p, k_l) by inner CRPS. Zero-aware by construction.
* **N3 is dropped** (revision 1 kept it as "descriptive only, not eligible"; its own frozen benchmark ranks it
  below KNN and the layer's served T1 answer is `baseline:KNN`, a control). Its exploratory numbers stay in
  section 5's history only.
* **Code still to be written before the freeze:** the N1 and N2 forecast functions (registered in
  `next_regional_run.CANDIDATES`) and their nested selection. They pass through the same `predict`, so the
  canary covers them; the canary and the tests are re-run at the freeze.

**Exploratory numbers already known (training, not evidence).** General layer's run (`reports/synarch_recon_results.md`,
T1, its own folds): KNN log MSE 1.259, RIDGE 1.344, IDENTITY 1.536, MEAN 1.547, its learned joint mixture 1.371
(8.9% worse than KNN). This protocol's feasibility: section 13.

## 6. Selection (training only, once, before any test forecast)

Whole-sample nested CV on the 101 training samples: 5 outer x 3 inner folds by `synarch_recon.nested_folds` with
salt `SYNARCH-next-regional-2026-11-11` (`feasibility.json`, `folds`). Settings of every method from the inner
folds (`next_regional_core.choose_settings`); every method scored on the outer folds by the declared scorer.

A candidate is **eligible** if, on the outer folds pooled over the 101 samples, it meets the decision rule of
section 8 against every control. Among eligible: lowest animal-mean CRPS, then lowest log MSE, then N1 before N2.
**If nothing is eligible, the test is not run and no test regional row is ever opened** (`score --selected none`
writes "not_run" without reading the regional table; tested). The selection record is written before any test
forecast. The test forecasts of every method use settings chosen by 5-fold whole-sample CV on all 101 training
samples (`next_regional_run.py select`).

## 7. Measures

* **Primary distribution:** animal-mean quantile CRPS.
* **Primary point (non-inferiority):** animal-mean log MSE.
* **Calibration:** marginal 90% coverage; conditional non-zero 90% coverage (animal mean, and pooled per region);
  presence gap per region.
* Described: log MAE, profile correlation, width, interval score, presence Brier, per region, per source label
  and per soma label tables (labels with fewer than 5 test samples marked sparse; only Adcyap1, Sst, Penk have 5
  or more).

## 8. Decision rule (applied once on the test samples, to the selected candidate only)

**Pass** only if all three legs hold:

1. **CRPS:** at least 2% below each distribution control (MEAN, IDENTITY, KNN@crps, RIDGE@crps, KNN-CONF).
2. **Points not worse:** log MSE at most 1% above the best point control (MEAN, IDENTITY, KNN, RIDGE; log-MSE tuned).
3. **Calibration:**
   * 3a. animal-mean marginal 90% coverage in [0.85, 0.95];
   * 3b. animal-mean conditional non-zero 90% coverage in [0.85, 0.95];
   * 3c. among regions with at least 30 non-zero test truths in at least 8 test samples, at most 5 with pooled
     conditional non-zero coverage below 0.75;
   * 3d. at most 5 of 49 regions with |mean p_nonzero - observed non-zero share| > 0.10 (pooled).

Outcomes, stated as such: **pass**; **better CRPS without valid calibration** (1, 2 hold; 3 fails); **better CRPS
but worse points** (1 holds, 2 fails); **no gain over the controls** (1 fails); **not run**. The descriptive label
**"also better points"** is added when log MSE is at least 2% below every point control; it decides nothing.

Why this rule (finding 2, decided in this revision). Revision 1 required a 2% log-MSE gain over KNN. The
reviewer's training check, reproduced here (`feasibility.json`, `near_knn`), shows the blend points at 0.75 /
0.5 / 0.25 are +0.56% / -0.46% / -3.1% against KNN: a 2% point gain is out of reach for N1, so the old rule tested
almost nothing. The question the user asked is about calibrated uncertainty, so the weight moves to the proper
score with points held non-inferior. **This change was made after seeing training numbers (exploratory) and
before any candidate was run.**

Why leg 3 has teeth (finding 1). In revision 1 the per-region guard used all-truth coverage, which zeros pass for
free. Now: on training, KNN fails 3b (0.825) and 3c (8 regions below 0.75); RIDGE@crps and KNN-CONF fail 3c and
3d (presence overstated in all 49 regions); MEAN fails 3a (0.952); only IDENTITY passes all four (section 13).
A perfectly calibrated forecast fails 3c or 3d at test size with probability about 0 (2,000 draws of 40 training
samples with binomial coverage; neurons treated as independent, which clustering by sample would widen). **The
thresholds (0.85 to 0.95, 0.75, 30 / 8, 0.10, 5 regions) were set after seeing the controls' training
calibration; recorded as such.**

The 2% minimum on CRPS is the general layer's (largest no-signal gain over the best comparator in its synthetic
worlds was 1.0%). Secondary presence: lower Brier than every control by at least 2%; reported only.

## 9. Intervals, multiplicity, chance of a false pass

* Paired per-sample differences (candidate minus each control), mean, samples favouring the candidate and a 95%
  bootstrap interval over samples (10,000 resamples, seed 20261111). **Descriptive; they do not enter the rule.**
* One selected candidate, one primary decision that is the conjunction of legs 1 to 3.
* **Null pass rate (finding 8).** The rule compares point estimates with a threshold and does not require an
  interval to exclude zero. If the candidate's true CRPS equals a control's, the chance of an observed gain of 2%
  or more against that one control at 40 samples is about 1e-6 when the per-sample SD of the difference is 0.012
  (a near-KNN candidate; KNN-CONF minus KNN@crps) and about 0.23 when it is 0.074 (as different from KNN as
  IDENTITY is). Because the gain must hold against every distribution control, including KNN@crps, the rate is
  at most the smallest of these, which for any candidate close to KNN is negligible
  (`feasibility.json`, `null_pass_rate_crps_leg_one_control`; normal approximation).

## 10. Run-once mechanics

1. **Freeze** (after review): write and review N1 and N2 code; re-run the tests (canary included) and
   `check-loader`; run the nested selection (section 6) and `next_regional_run.py select`; pin this protocol, the
   code (`next_regional_*`), the sealed loader's SHA-256, the manifest and table hashes, the selection records.
2. **Predict:** `python3 -I analysis/next_regional_run.py predict --part test --selection .../selection.json --out
   data/interim/next_regional/test_forecasts`. It prints the manifest SHA-256 = the **commitment**, recorded by the
   user or the lead (a commit or a pinned note) **before plan stage B (test trace fetch)**. No regional row of a
   test sample is read by it (canary). **If the commitment is not recorded before stage B starts, the result
   automatically becomes option (b), "re-use of consumed outcomes" (section 12), and is reported only so.**
3. Wait. The current evaluation runs unchanged (its predict, stage B fetch, admit, seal, score).
4. **Score**, only after the current evaluation's `score --stage test` has completed:
   `python3 -I analysis/next_regional_run.py score --part test --forecasts ... --commitment <sha> --selected
   <N1|N2|none> --after-engine-score <the engine's test score file> --out data/interim/next_regional/test_scores`.
   It refuses unless the manifest hashes to the commitment, every forecast and code file is unchanged, the engine
   score file exists (its SHA-256 is recorded), and no scoring was started before (`scoring_started.json` is
   created exclusively). Then it reads the regional rows of the committed test neurons only. Nothing is tuned
   afterwards.
5. **Failure handling.** Before step 2 a bug may be fixed with a new freeze record that says why. After the
   commitment: a crash of the scorer may be fixed only if forecasts and rule are untouched (logged; commitment
   re-verified; the code pin then names the fixed scorer, so a fix needs a written record). A missing candidate
   forecast counts as MEAN's. A label absent from training uses the declared fallback (test: Vip has no test
   sample; Crh and Trh have one training sample each).

## 11. What a result would and would not mean

* **Pass:** from soma position and source label, the candidate gives 49-region predictive distributions that are
  better (CRPS) than KNN, ridge, the means and conformalised KNN, with points no worse and calibration holding for
  zeros and non-zeros, on 40 brain samples held out of every fit. It would not show paths, synapses, causes or a
  label effect (label, sample and batch are confounded), and, given section 12, it is not a pristine test.
* **Fail:** at this sample size the candidate is not shown to beat the simple methods; the system keeps saying
  "no better than baseline" for regional profiles. It would not show that no better method exists.
* **Not run:** training CV already said no; the test samples' regional rows stay unread by this protocol.

## 12. Exposure and which animals can still test this

**Facts.** (i) The paper's 49-region table for all 7,180 neurons has been on disk since 2026-11-10. (ii) The other
agent's regional model (`reports/brm_real_results.md`) was fitted and scored by 5-fold whole-animal CV on **all 166
samples**; its files (`data/interim/brm_real/truth.npz`, `predictions/`) contain held-out outcomes. (iii) The
general layer's author, the engine's author and the writer of this draft read that report's aggregate numbers
(166 samples together; no per-sample or per-neuron held-out value). (iv) **The published paper's own label-level
projection summaries pool all 166 samples** (finding 9), another aggregate exposure of held-out outcomes.
(v) The current evaluation opens validation traces at its validation `admit`, and test traces after plan stage B;
once a trace is opened, its regional profile can be derived from it.
**So no Jiao validation or test sample is pristine for the regional question, whatever is done next.**

| option | what | trade-off |
|---|---|---|
| (a) seal before opening | freeze, commit test forecasts before stage B, score after the current evaluation's test score | 40 samples; contamination limited to aggregate exposure (ii to iv); no change to the current evaluation; must be ready before stage B, otherwise it becomes (b) |
| (a') also validation | the same for validation, sealed before the current evaluation's validation `admit` | 65 samples, but the window is short (hours) and a rushed review is worse than lower power |
| (b) re-use | score any samples later, labelled "re-use of consumed outcomes" | no deadline; weak evidence, never "pristine" |
| (c) new animals | Li et al. 2024 mouse 18696 (115 Oxt neurons, **101 of them hypothalamic**; same laboratory and pipeline; CC BY 4.0 archive; described there as a lactating female, unlike the Jiao cohort); MouseLight (Janelia, 97 hypothalamic neurons from 22 mice, another laboratory, no genetic label); Ding et al. 2025 (266 hypothalamic neurons from 5 mice, same institute, no label, LHA/ZI/PSTN only); counts from `reports/axon_dataset_other_sources.md` | the only genuinely unseen animals; none has the published 49-region table, so profiles must be computed from traces by a re-implementation of the authors' length measure (validated first on training traces); only one animal has a source label; fetching needs the user's decision |

**Recommendation: (a) on the 40 test samples**, labelled everywhere "held-out samples of the locked split; their
regional outcomes were modelled in aggregate before the split was locked (brm_real) and summarised over all
samples in the published paper; not pristine" (`next_regional_run.EXPOSURE_LABEL`). Later, (c) MouseLight + Ding
as a label-free transfer check under its own protocol.

**Why (a) does not contaminate the current evaluation.** Its methods are frozen and its forecasts committed before
stage B; it never reads held-out regional rows before the engine's test score. It runs no engine command, writes
only under `data/interim/next_regional/`, reads before scoring only training rows and the test samples' permitted
soma inputs (the same inputs the engine reads at `predict --stage test`), and opens test regional rows only after
the engine's test score is complete, so no step of the current evaluation can depend on them. One residual: if
this protocol is frozen after the engine's validation result is known, that result (geometry on other samples)
was seen before the freeze; the freeze record says whether it was.

## 13. Feasibility (TRAINING only; exploratory; commands and outputs)

Commands:
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -I
analysis/next_regional_feasibility.py run` (about 7 minutes) -> `data/interim/next_regional/feasibility.json`,
`per_region.tsv`, `SHA256SUMS.json`; and `python3 -I analysis/next_regional_run.py check-loader` ->
`loader_check.json`. Loader audit: 4,448 training rows, 2,732 sealed rows dropped unparsed, regional table SHA-256
`2072bc20...e63c`. These samples were already used by the general layer, by brm_real, and by revision 1 of this
draft and its review.

* **Zeros.** 72.8% of all region values are exactly zero; per region from 18% (LHA) to 96% (NTS, Isocortex);
  median 79%; 22 regions more than 80% zero, 9 more than 90%. A neuron has a median of 13 non-zero regions
  (10th to 90th percentile 5 to 22); 10 training profiles are all zero.
* **Support.** Every region has at least 161 non-zero training neurons in at least 51 samples (fewest: NTS). At
  50 / 10, 100 / 20 and 200 / 30 (neurons / samples): 49, 49 and 47 regions. In 2,000 draws of 40 training
  samples all 49 regions met the 3c eligibility (30 truths in 8 samples) every time.
* **Controls, nested CV, declared scorer (animal mean):**

| control | log MSE | CRPS | marginal cov. | cond. non-zero cov. | regions < 0.75 (3c) | presence off (3d) | presence Brier |
|---|---:|---:|---:|---:|---:|---:|---:|
| KNN (k 40, all, raw, every fold; = KNN@crps) | 1.255 | 0.4345 | 0.944 | 0.825 | 8 | 0 | 0.136 |
| KNN-CONF | 1.255 | 0.481 | 0.900 | 0.888 | 12 | 49 | 0.319 |
| RIDGE (degree 2, lambda 0.01) | 1.346 | 0.525 | 0.898 | 0.883 | 13 | 49 | 0.337 |
| RIDGE@crps (degree 2, lambda 1 or 10) | 1.467 | 0.516 | 0.912 | 0.867 | 14 | 49 | 0.328 |
| IDENTITY (neuron weighting) | 1.532 | 0.513 | 0.940 | 0.853 | 1 | 0 | 0.162 |
| MEAN | 1.544 | 0.517 | 0.952 | 0.894 | 0 | 0 | 0.163 |

  KNN's lowest per-region conditional non-zero coverage: Isocortex 0.44, PA 0.61, HPF 0.61. Residual-based
  distributions (KNN-CONF, RIDGE) overstate presence by 0.11 to 0.52 per region, because adding residuals to a
  positive point leaves little mass at zero. KNN beats RIDGE on log MSE in 86 / 101 samples (difference -0.091,
  95% bootstrap [-0.118, -0.062]) and KNN@crps beats RIDGE@crps on CRPS in 97 / 101. KNN-CONF is worse than
  KNN@crps on CRPS in all 101 samples (+0.047, [0.045, 0.049]).
* **Near-KNN spreads (the reviewer's check, reproduced; finding 2).** Blend points w KNN + (1 - w) RIDGE against
  KNN, log MSE: w 0.75 +0.56% (SD of the per-sample difference 0.040), w 0.5 -0.46% (0.077), w 0.25 -3.1% (0.112).
  CRPS, KNN-CONF minus KNN@crps: SD 0.012.
* **What the test can detect** (paired t, two-sided 5%, power 80%; relative to KNN). Smallest detectable
  difference at 25 / 40 / 65 samples:

| per-sample SD from | log MSE | CRPS |
|---|---|---|
| near-KNN differences | 1.9-5.2% / 1.5-4.0% / 1.1-3.1% | 1.6% / 1.2% / 0.9% |
| differences between controls | 6.8-11.6% / 5.3-9.1% / 4.1-7.0% | 1.6-9.9% / 1.2-7.7% / 0.9-6.0% |

  A 2% CRPS gain is detectable at 40 samples when the SD is at most 0.019, which holds for near-KNN candidates
  (0.012). So the CRPS leg can see a 2% gain for candidates close to KNN. The old point leg could also have seen
  a 2% gain at near-KNN spreads (0.040 is below 0.055), but the training blends show no such point gain to see.

## 14. Changes from revision 1 (review findings)

| finding | change |
|---|---|
| 1 coverage leg | per-region guard now on **conditional** non-zero coverage in eligible regions (30 truths / 8 samples), plus 3b and presence calibration 3d; the review's literal fix (marginal coverage of non-zero truths >= 0.80 per region) was not adopted because a perfectly calibrated forecast fails it wherever more than 80% of values are zero (22 regions); thresholds recorded as set after seeing training numbers |
| 2 N1 vs 2% point leg | rule redeclared: CRPS 2% + points non-inferior (1%) + calibration; "also better points" label descriptive; expectation for N1 stated; decided without running N1 or N2 |
| 3 RIDGE | out-of-fold residuals; lambda grid to 1000; CRPS tuning now interior (1 or 10) |
| 4 conformal KNN | KNN-CONF is a distribution control |
| 5 loader, scorer | written (`next_regional_run.py`), canary and refusal tests, real `check-loader` on training |
| 6 viewer | new evidence level; current export refuses; no export as `independent_test` |
| 7 CRPS definition | one declared quantile scorer for feasibility and test; numbers recomputed with it |
| 8 null pass rate | stated in section 9 |
| 9 exposure | paper's pooled summaries added; 18696: 101 hypothalamic of 115; fallback to (b) if stage B is missed |
| 10 inputs | decided: soma position and source label only, for all methods |

## 15. Open decisions (lead or user)

1. Option (a) on 40 test samples, or (a') with validation (only if the review ends before the current
   evaluation's validation `admit`).
2. N1 as declared is expected to fail (section 5). Keep it (and accept a likely "not run"), drop it, or replace it
   before the freeze by a zero-aware variant (which would be a new design choice made after seeing training
   control numbers, to be recorded as such).
3. The redeclared rule (section 8: CRPS 2%, points non-inferior 1%, calibration) versus revision 1's rule
   (2% on both). Confirm before the freeze.
4. The calibration thresholds of leg 3 (set after seeing training numbers): keep or change before the freeze.
5. The viewer's owner: add an evidence level for "held-out, previously modelled in aggregate".
6. Hemisphere and reconstruction type: decided "not used"; overrule only before the freeze.
7. Who records the commitment before stage B, and who does the second independent review.
