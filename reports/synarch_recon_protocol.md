# SYNARCH reconstruction layer: masked-variable benchmark protocol

**Status: DRAFT. Written 2026-11-10 by the builder; revised 2026-11-11 by the fixer after a full review, the other
agent's review and three attacks on the first fixes (section 20). NOT frozen. Nothing has been fitted or scored on
real outcomes.** This file becomes binding when
`python3 -I analysis/synarch_recon_benchmark.py freeze --reviewed "..."` writes
`data/interim/synarch_recon/pre_run_declaration.json`. That declaration holds the SHA-256 of this file, of the
code, of the tests and of every input, the folds, the grids and the seeds. `run` then checks all of them again
and runs once.

Code: `analysis/synarch_recon.py` (core), `analysis/synarch_recon_modules.py` (modules and baselines),
`analysis/synarch_recon_benchmark.py` (commands), `tests/test_synarch_recon.py`. Results will be written to
`reports/synarch_recon_results.md` by the `report` command, only after the run. No biological result exists yet.

## 1. What is being tested

SYNARCH is meant to estimate properties of a unit that were not measured, say how uncertain the estimate is, and
refuse when the data cannot identify the answer. This benchmark asks one bounded question about the first learned
engine of that layer: **on four masked-variable tasks, does the unchanged joint mixture model
(`analysis/brm_joint.py`) beat simple rules on neurons of brain samples it never saw?**

The answer decides what the reconstruction API returns for each task: the engine's value, or a baseline-level
value that is flagged as such, or a refusal (section 12).

It is not a test of 3D axon paths, synapses, gene expression or causes. The 49 regional lengths are published
summaries of an axon, not the axon.

## 2. Exposure: this is exploratory

* The 49-region table has been on disk since 2026-11-10. Another agent's regional model was fitted and scored on
  all 166 brain samples, training samples included (`reports/brm_real_results.md`).
* The builder read that report before writing this protocol, as instructed. Two choices here were informed by it:
  the component grid goes up to the model's limit of 16 (8 was chosen at the top of that run's grid in all folds),
  and the engine's conditional median is reported next to its mean (copying won on absolute error there).
* The presence/level representation of section 7 was added after a synthetic probe and after seeing one
  marginal fact of the training rows: 72.8% of the regional values are exactly zero.
* The changes of section 20 were made after reviews of the code on synthetic data only. No reviewer, attacker or
  fixer fitted or scored anything on real outcomes.
* No confirmatory claim will be made from this benchmark. The best status a task can reach is
  `exploratory_supported`.

## 3. Data and sealing

**Data: the TRAINING partition of Jiao et al. 2025 only: 101 brain samples, 4,448 neurons**, defined by
`reports/axon_dataset_manifest.tsv` (column `part == 'train'`; locked by `reports/axon_dataset_plan.md`).

| file | what is read | SHA-256 |
|---|---|---|
| `reports/axon_dataset_manifest.tsv` | rows with `part == 'train'`: neuron, sample, source label, soma label | `649f829a592cd3d13bcda5cffdbd049e759568d38bacb144f6a61220610a1ca4` |
| `data/raw/jiao2025/All_soma_infor_HY_neurons.csv` | columns `file, neurons, geno, location, AP, DV, ML` only | `cc01d048bf5a7c9c9f46869fb750f2d67a98ca934903d2e2e71c1a67b1d97cb9` |
| `data/raw/jiao2025/all_neuron_soma_infor.xlsx` | columns `neurons, location` only | `7b7ae75b94688619dbb367d5f9e71d8b157615624932375a00e80c796f21f162` |
| `data/raw/jiao2025/all_neurons_projection_length_fig1.csv` | ID and 49 region columns | `2072bc20550607c195be0a42cee9fcbe0e0ab5aaebfaee0a1016b2a43d35e63c` |
| `data/interim/brm_data_contract.json` | namespaces, grains, allowed joins | `46d2d1b0c9986ea110612c6549d163664e0551e04d8e5650f1d6cee88e27b510` |
| `analysis/brm_joint.py` | the unchanged joint model, hash checked at import | `ae7f5ec432bcbd2b0eebb65595549bc12158b1b7c047f5ace167c7840c0c937e` |

**Sealing.** Rows of validation and test brain samples are sealed. Every reader keeps a row only if its sample ID
(the ID prefix) is a training sample, decided before any other cell of the row is parsed; the regional table is
read line by line and a sealed line is dropped unparsed. The manifest reader keeps nothing of a non-training row,
not even its ID, so everything downstream is a whitelist. The three tables are joined by neuron ID, never by
position, and must each hold exactly the 4,448 training neurons. Every fit re-checks that all its rows belong to
training samples. In the synthetic worlds the sealed rows hold unparseable text, so a reader that parsed one would
fail; changing them leaves the loaded table identical (tested).

**Never read:** the `type` and `name` columns of the CSV, the `cluster` and `Class` columns of the xlsx (derived
from the axons; registered as banned predictors), any trace (`.swc`), any earlier fitted state
(`data/interim/brm_real`, `brm_joint`, `brm_geometry`, `d1*`), the Epi-retro-seq evaluation set. The xlsx reader
asks the library for columns C to D only and refuses if a foreign column lies between the two allowed ones; the
library streams the sheet, so other cells pass through its parser but are never returned to this code.

**What was read from the real cohort before the freeze** (builder 2026-11-10; fixer 2026-11-11, the same reads
again through the test file and the `inputs` command):

* the manifest, and the allowed input columns of the two soma tables, training rows only;
* the header of the regional table;
* the training rows of the regional table, for two purposes only: per-region marginal summaries
  (`data/interim/synarch_recon/input_audit.json`: zero fraction, maximum, mean and SD of log10(1 + length) per
  region, number of all-zero profiles) and to check that the rows load as records (counts, statuses, finite and
  non-negative values);
* the local CCF ontology (`data/raw/atlas/structure_graph_1.json`), to check the T4 block (section 6).

No input was related to a regional value, no task model was fitted, and no fold was scored on real data.

## 4. Unit of independence

The unit is the **brain sample** (the ID prefix). The paper reports 166 mice and the data hold 166 sample IDs; no
statement was found that no mouse gave two samples, so this protocol says "brain samples" and does not equate a
sample with an animal. Every sample carries exactly one source label, so a label effect is confounded with sample
and batch. All splits are by whole sample and all metrics are averaged per sample first.

Training partition (inputs only): 1 to 135 neurons per sample (median 40); 16 source labels:
Adcyap1 12 samples / 709 neurons, Agrp 3 / 87, Avp 2 / 13, Crh 1 / 16, Nts 9 / 542, Oxt 8 / 290, Pdyn 11 / 445,
Penk 10 / 304, Pmch 3 / 107, Pomc 6 / 238, Sst 13 / 813, Tac1 6 / 386, Tac2 10 / 371, Trh 1 / 42, Vip 2 / 16,
Wt(orexin) 4 / 69. `Wt(orexin)` is the orexin promoter label, not a Cre line; a source label is a genetic label,
not measured gene expression.

## 5. Variables and scales

| variable | kind | scale used for modelling and scoring |
|---|---|---|
| `soma_position` | vector AP, DV, ML | micrometres as published (CCFv3 2017 frame, mirrored to ML <= 5,700; axes, unit and mirroring are inferred for this table, `reports/d109_jiao_conventions.md`) |
| `source_label` | categorical, 16 values | the table's own spelling |
| `regional_axon_length` | vector, 49 regions | **log10(1 + published length)**; the length unit is indicated (micrometres), not resolved |
| `soma_label_fine` | categorical | loaded and recorded; **not an input of any task** (it is a function of the soma position) |

The 49 columns are not a partition of the brain and do not sum to the axon length; MY is the ancestor of IO, NTS,
PARN and IRN. 10 of the 4,448 training profiles are zero in all 49 regions; they stay in training and scoring.
Every method's regional prediction is clipped at zero from below (the same rule for all).

## 6. The four tasks

Each neuron is predicted once, by states fitted on other brain samples.

| task | hidden (to be estimated) | observed (inputs) |
|---|---|---|
| **T1** | regional profile, all 49 regions | soma position + source label |
| **T2** | soma position (AP, DV, ML) | the measured regional profile, 49 regions + source label |
| **T3** | source label | soma position + the measured regional profile |
| **T4** | the 33 regions outside the hypothalamus (MRN ... OLF) | the 16 measured hypothalamic regions (ARH ... ME) + soma position + source label |

T1 is the only hidden-axon forecast: nothing derived from the neuron's axon is an input. T2, T3 and T4 are
**declared** as taking a measured profile; the API accepts a measured profile as input in these three conditionals
only, and the axon adapter rejects it always.

T4's block was fixed without outcomes: the first 16 published columns are exactly the descendants of HY among the
49 (checked against the local CCF ontology), and MY with its four descendants lies entirely in the hidden block,
so no parent sits on one side with its child on the other.

T3 rule for labels without training support: a neuron whose true label has no training sample in its fold cannot
be scored by any method; it is excluded from accuracy and log-loss for all methods and counted (Crh and Trh have
one training sample each, so their neurons fall here).

## 7. The learned engine

**What is unchanged and what this layer adds.** The model is the unchanged `JointModel` of
`analysis/brm_joint.py`: a latent-class mixture with diagonal Gaussian factors for continuous fields and tables for
categorical fields, fitted by EM on observed entries. Its file is hash-checked at import and never edited. Two
things around it are this layer's own: (a) a batch reader that computes the model's conditionals for many rows at
once (it reproduces `JointModel.infer()` to 1e-9, tested); (b) the **presence/level ("hurdle") encoding** of a
zero-inflated region and the mean and quantiles that are read from it. The encoding only changes how a region is
presented to the same model (two fields instead of one); the formulas for its mean and quantiles are written in
this layer, not in the model, and are checked against the model's own sampler (tested). Whether it is used is
decided on inner folds; setting `representations` to `["log"]` removes it.

One fitted state answers all four tasks: evidence fields condition the class probabilities, hidden fields are
read from the conditional mixture (mean, median, central 90% interval; class probabilities for the label).

Two settings are chosen on inner folds, per task and outer fold:

* **number of components K** in {1, 2, 4, 8, 12, 16} (16 is the model's limit). **With K = 1 the engine has learned
  nothing from its inputs:** its prediction is the neuron-weighted training mean (label target: the
  neuron-weighted label frequencies), up to the model's small prior smoothing. The report says for every task in
  how many folds K = 1 was chosen;
* **representation of the regional profile**, both using the same unchanged model:
  * `log`: one Gaussian field per region on log10(1 + length): 53 fields (3 soma, 1 label, 49 regions);
  * `hurdle`: per region a categorical presence field (length above zero or not) and a Gaussian level field that
    is observed only when the length is above zero and missing otherwise: 102 fields. Zero is then an outcome
    with its own probability. The conditional mean is the sum over classes of P(class) x P(present | class) x
    mean level; quantiles come from the zero-inflated mixture and are exactly zero where the mass at zero reaches
    them.

Fixed constructor settings: `n_init = 5`, `max_iter = 300`, `tol = 1e-7`, `variance_floor = 1e-3`,
`dirichlet_alpha = 1.1` (the last two are the library defaults). Seeds are a hash of
(salt, representation, stage, outer fold, inner fold, K). The model standardises every continuous field with
training statistics itself; there is no other scaling choice for the engine. A fit that stops at the iteration
limit is counted and reported, never extended.

Point estimates: the conditional **mean** (regional values clipped at zero) is the engine's estimate for every
headline number and verdict. The conditional **median** is reported as a second, declared summary and scored on
absolute error only as a secondary line.

## 8. Baselines and references

All use the training rows of the fold only.

Always reported:

* **MEAN**: the training mean of the target (T3: the label frequencies). It has ONE declared setting, chosen on
  inner folds like k and lambda: **weighting** in {`sample`: every training brain sample counts once; `neuron`:
  every training neuron counts once}. Reason: the two are estimates of the same thing; with very uneven numbers of
  neurons per sample the sample-weighted one is the noisier, and a one-component engine (which IS the
  neuron-weighted mean) was "better than MEAN with an interval excluding zero" in most no-signal worlds while MEAN
  could only weight samples (review; test
  `test_a_one_component_engine_is_the_neuron_weighted_mean_and_does_not_beat_the_tuned_mean`).
* **IDENTITY**: the source label's training mean, with the same weighting choice (its own, chosen on inner
  folds); global fallback (MEAN with that weighting) for a label with no training sample. Not applicable to T3
  (the label is the target).
* **COPY**: the value of the nearest training neuron of another brain sample in the observed space. Distance is
  Euclidean: T1 on the soma in micrometres; T2 on the 49 log-lengths; T3 and T4 on all observed coordinates
  standardised with the fold's training mean and SD (mixed units). Ties go to the first training neuron in neuron-ID
  order. T3: the copied label gets probability (1 + prior) / 2, the others prior / 2, where the prior is the
  sample-weighted label frequency (fixed; it is not the tuned MEAN).
* **COPY, same source label**: the same among training neurons with the query's label; fallback to COPY when the
  label has no training neuron. Not applicable to T3.

Two further simple references, tuned on inner folds:

* **KNN**: average over the k nearest training neurons of other samples; k in {5, 10, 20, 40, 80}; pool in
  {all, same source label (fallback: all)}; distance scaling in {raw, standardised} for T1 and T2, standardised for
  T3 and T4. T3: label counts plus the prior, divided by k + 1.
* **RIDGE**: ridge regression of the target on the standardised observed coordinates and the one-hot label; penalty
  lambda x n with lambda in {0.001, 0.01, 0.1, 1, 10}; T1 also chooses between linear and quadratic soma terms.
  T3: multinomial logistic regression with a ridge penalty, lambda in {0.0001, 0.001, 0.01, 0.1}, L-BFGS, at most
  500 iterations (whether it stopped by itself is recorded).

MEAN and IDENTITY also give an empirical central 90% interval (quantiles of the training values, of the label's
rows for IDENTITY; the same for both weightings) so that the engine's interval has something to be compared with.

## 9. Splits

Whole brain samples, by a salted hash; salt `SYNARCH-recon-2026-11-10`, written before the folds were first
computed and not to be changed.

* **Outer: 5 folds.** Inside each source label the samples are ordered by SHA-256(salt:sample); the labels are
  taken in sorted order of the table's spelling and the samples are dealt in turn. Fold sizes differ by at most
  one sample and every label is spread as evenly as its sample count allows.
* **Inner: 3 folds** inside each outer-training part, the same rule with salt `<salt>:inner:<outer fold>`.

Result on the 101 training samples (inputs only; SHA-256 of the fold lists
`032d4b796b09881911c6fcf052453be72d90bd2d07164f547e2c33bee0fe0d86`): 21 / 20 / 20 / 20 / 20 samples and
955 / 762 / 863 / 795 / 1,073 neurons. Avp and Vip (2 samples each) are split over folds 0 and 1; Crh (1 sample)
is in fold 2 and Trh (1 sample) in fold 4, so their neurons are predicted without any training sample of their
label: the engine then marginalises the label, IDENTITY and same-label COPY fall back to their global rule, and
T3 cannot score them.

## 10. Selection

Every choice is made on the inner folds of the outer-training part only: engine (representation, K), MEAN
(weighting), IDENTITY (weighting), KNN (k, pool, scaling), RIDGE (lambda, degree). The criterion is the task's
primary metric (section 11), averaged per brain sample inside each inner validation fold and then over the three
inner folds. Ties go to the earlier setting in the declared order (log before hurdle, smaller K; `sample` before
`neuron`; the other grids in the order above); for T3 a tie in accuracy is broken by lower log-loss. The chosen
setting is refitted on the whole outer-training part and predicts the held-out samples once. No outer result
changes a setting.

A fit that raises a numerical error is recorded; a setting with a failed inner fit is not eligible; if the chosen
setting fails on the outer-training part, the next-best eligible setting is used and listed as a deviation.

Size: 5 x 3 x 12 = 180 inner engine fits, at most 20 outer fits, at most 4 all-training fits.

## 11. Metrics

Per neuron, then the mean over the neurons of a brain sample, then the mean over brain samples.

* **T1, T4**: mean squared error (**primary**) and mean absolute error over the target regions on log10(1 + length);
  Pearson correlation of the predicted and measured profile across the target regions (undefined for a constant
  measured profile: excluded from correlation only, with counts; a constant prediction scores 0).
* **T2**: Euclidean error in micrometres: mean (**primary**) and median over a sample's neurons; absolute error
  per axis.
* **T3**: accuracy (**primary**) and log-loss (probability floor 1e-12).
* **Intervals**: for every inferred 90% interval, empirical coverage (share of target coordinates inside, per
  neuron, per sample, over samples), mean width, and coverage per coordinate. T3: coverage and mean size of the
  smallest label set holding 90% probability.

## 12. Comparisons, verdict, status, multiplicity

For the engine against each comparator and metric: the paired brain-sample difference, its mean, a percentile
bootstrap interval over brain samples (10,000 resamples, fixed seeds, 95%), the number of samples favouring the
engine, and the **relative gain** (the mean advantage divided by the comparator's own mean). **All intervals are
descriptive**: folds share training data, labels share batches, and the data were seen before. A paired difference
smaller than 1e-12 relative to the metric is an exact tie (two predictors that agree to rounding error are not
better than one another).

Per task the report says whether the engine beat **every, some or none** of (a) the always-reported baselines and
(b) all comparators, on the primary metric, in three ways that are shown side by side:

* *any amount*: the engine's mean is better, however little;
* *with interval*: also, the descriptive 95% interval of the paired difference excludes zero;
* *clear*: also, the relative gain over that comparator is at least the **declared minimum of 2%**
  (`min_relative_gain = 0.02`).

Why a minimum: when two predictors are almost the same the paired interval is extremely narrow, and a difference
of a hundredth of a percent falls outside zero. In no-signal synthetic worlds the interval rule alone called a
task supported with gains between 0.005% and 1.6%; the planted worlds that deserved it showed gains of 5% to 13%
(review, and section 17). The number was fixed before any real outcome was scored.

**Status registered for the task's conditional in the reconstruction API. This is the only decision taken from the
benchmark:**

* `exploratory_supported` only if the engine is better than **every** comparator (the four baselines and the two
  references) on the primary metric, **each** descriptive 95% interval of the paired difference excludes zero,
  **and** the relative gain over each comparator is at least 2% (equivalently: at least 2% over the best one);
* otherwise `no_better_than_baseline`. The API then returns, by the declared default policy, the answer of the
  best comparator, flagged `no_better_than_baseline` and `baseline_level`; with the other policy it refuses.

Conditionals that were not benchmarked are `not_validated` (refused by default). Nothing can become `validated`
here: independent confirmation does not exist.

**Multiplicity.** Four tasks, one primary metric each: four applications of the status rule, not adjusted for
one another. Every other interval in the report (other metrics, single comparators, source labels: about ninety)
is an unadjusted description and decides nothing. How often the rule calls at least one of the four tasks
supported when there is nothing to learn was measured on no-signal synthetic worlds with the frozen pipeline
(`nullrate` command; outputs `data/interim/synarch_recon/synthetic/nullrate_tiny.json` and `nullrate_real.json`):

* **Small shape (32 training samples, about 450 neurons, 8 regions; the complete small grid, K in {1, 2, 4})**: 200 worlds, 800 task results. **Worlds with at least one task called supported: 0 of 200.** With the interval rule alone, without the minimum gain, it would have been 6 of 200 (task results: T1 3, T2 0, T3 0, T4 3; their gains over the best comparator: T1: 0.30%, 0.30%, 0.44%; T4: 0.26%, 0.35%, 1.00%).
  Largest gain over the best comparator in any of these worlds: T1 0.58%, T2 0.54%, T3 13.47%, T4 1.00%. The engine beat every always-reported baseline *by any amount* in 99 (T1), 63 (T2), 90 (T3), 87 (T4) of 200 worlds, and *clearly* in 0 (T1), 0 (T2), 0 (T3), 0 (T4).
* **Real shape (101 training samples, 4,448 neurons, 49 regions, the real label counts; every declared setting except that the engine's grid stops at K = 4, because the three largest values take most of a run's time)**: 8 worlds, 32 task results. **Worlds with at least one task called supported: 0 of 8.** With the interval rule alone, without the minimum gain, it would have been 0 of 8 (task results: T1 0, T2 0, T3 0, T4 0; their gains over the best comparator: none).
  Largest gain over the best comparator in any of these worlds: T1 0.00%, T2 0.00%, T3 7.23%, T4 0.00%. The engine beat every always-reported baseline *by any amount* in 4 (T1), 3 (T2), 5 (T3), 4 (T4) of 8 worlds, and *clearly* in 0 (T1), 0 (T2), 0 (T3), 0 (T4).
* One further no-signal world at the real shape with the complete grid (K up to 16; section 17): no task supported.

What these numbers say and do not say. Zero observed is not zero: with 200 worlds the chance of at least one false 'supported' among the four tasks is below about 1.5% at the small shape (95% bound, rule of three). The real shape was simulated 8 times with a shortened grid and once with the complete grid, which bounds little by itself; what it adds is the size of chance differences there. For T1, T2 and T4 the largest gain over the best comparator in any no-signal world was 1.00% at the small shape and 0.00% at the real shape, below the minimum of 2%: for these three tasks the minimum gain is what removes the narrow-interval artefact. For T3 the minimum gain protects nothing, at either shape: accuracy near chance level moves by several percent of itself by chance (largest gain in a no-signal world: 13.5% at the small shape, 7.2% at the real shape); there it was the requirement that all six intervals exclude zero that held, in each of the 209 worlds. 'By any amount' the engine beats every always-reported baseline in about half of the no-signal worlds: a comparison of means without a margin says nothing, which is why the verdict of a task is the clear comparison.

## 13. When is an interval's coverage flagged as checked

The flag `coverage_validated` on an answer means "the empirical coverage passed the rule below in this
exploratory benchmark". It is not a claim of calibration.

* T2 (continuous target): |coverage - 0.90| <= 0.05 (two-sided).
* T1, T4 (the truth has a mass at exactly zero, so a central interval cannot hit its level): coverage **at least**
  0.85 and mean width no larger than that of the MEAN baseline's unconditional interval. This rule is one-sided:
  an interval that covers 95% where 90% was asked passes it (the no-signal synthetic runs do exactly that), so
  for T1 and T4 the flag says "covers at least as often as stated and is no wider than the unconditional
  interval", nothing more. The share of non-zero truths inside their interval is reported as a description; it
  is not a calibration target (a central interval is exactly [0, 0] where zero is very likely).
* T3: coverage of the 90% label set **at least** 0.85 (one-sided, as a label set cannot hit a level exactly).

The flag applies to answers from measured inputs inside the training support, at level 0.90, for units of brain
samples outside the state. The API drops it otherwise (unseen label, input outside the training range, inferred
or unverified input, evidence that could not condition, unknown or seen brain sample, another level).

## 14. Breakdown by source label

For each source label: brain samples, neurons, every method's metrics, and whether the engine was better than
each comparator. Paired bootstrap intervals are given only for labels with at least five brain samples; smaller
groups are marked "descriptive only".

Also reported: neurons whose label had no training sample in their fold, neurons with an input outside the range
of their fold's training samples, the settings chosen per fold (K = 1 counted separately) and how often K sat at a
grid edge, fits that stopped at the iteration limit, failed fits and deviations.

## 15. Freeze, the single run, outputs

```
python3 -I -m pytest tests/test_synarch_recon.py -q
python3 -I analysis/synarch_recon_benchmark.py inputs                      # input audit (training samples only)
python3 -I analysis/synarch_recon_benchmark.py freeze --reviewed "who reviewed, and where"
python3 -I analysis/synarch_recon_benchmark.py run --workers 6             # once
python3 -I analysis/synarch_recon_benchmark.py status
```

* `freeze` reads inputs only. It writes the pre-run declaration (exclusive create): tasks, comparators, grids
  (the weighting of MEAN and IDENTITY among them), the minimum gain, the salt, the fold lists, the settings of the
  state reader (support margin 0.05, 60 quantile steps, state format), hashes of code, tests, this protocol,
  inputs and the data contract, library versions. (A synthetic world's declaration pins the code and the tests
  but not this file: the protocol binds the real cohort.)
* `run` refuses unless the declaration matches the current files exactly, writes an exclusive start marker, and
  refuses a second start. It saves every out-of-fold prediction and state before any metric is computed
  (`predictions/`, `states/`, `fits.json`, `selection.json`), then scores from the saved predictions
  (`metrics.json`, `evidence.json`, each written once) and writes `run_complete.json`, which binds by SHA-256 the
  declaration, the prediction manifest, the fit records, the selection, the metrics and the evidence. The
  benchmark is then complete and stands. `run` prints the SHA-256 of `run_complete.json`.
* Three steps follow in the same command; they read no new data, cannot change a metric, and can be repeated by
  hand (`serving`, `demo`, `report`): (1) one all-training state per task, with the engine setting chosen most
  often over the five folds (ties: earlier in the declared order); it is never scored, it is what the API uses
  for units outside the training partition, and its evidence is this benchmark; `serving.json` records each
  state's fingerprint and the SHA-256 of the `run_complete.json` it belongs to; a numerical failure there is
  recorded and changes nothing else; (2) the API demonstration (`api_demo.json`) on one predetermined unit, the
  first neuron of outer fold 0 in ID order, answered by states that never saw its sample; (3) the report.
* The reconstruction API over a finished run (`build_reconstructor`) serves nothing unless: the declaration still
  matches code and inputs; every file bound by `run_complete.json` has its recorded hash; every state file holds
  exactly the state whose fingerprint and training groups were recorded when it was fitted, made with the
  declared reader settings; `serving.json` belongs to this `run_complete.json`.
* The results are independent of the number of workers (tested).
* All outputs are under `data/interim/synarch_recon/`. Nothing is committed, downloaded or uploaded.

## 16. If something breaks

* `abort --reason "..."` sets aside a freeze that was never run, or a start that broke for a technical reason.
  It is **refused as soon as any outcome-derived file exists**: the predictions or their manifest, `fits.json`
  (it holds the inner-fold scores), `selection.json`, metrics, evidence or the completion record. From then on
  outcomes have been used and the run stands. The declaration, the start marker and fitted but unscored states
  move to `aborted/<time>/`, are listed in the next declaration and in the report, and nothing is deleted. What
  the abort record says about existing files is computed, not assumed.
* A run that was interrupted after its predictions were saved is finished, not restarted: `score` computes the
  metrics from the saved predictions and writes the completion record.
* **After a STARTED run was set aside, the next freeze must be the same benchmark.** Tasks, comparators, grids,
  metrics settings, minimum gain, salt, folds, inputs and reader settings are compared with the declaration that
  was set aside, and any difference is refused. Only code and protocol text may differ (the repair of the
  recorded defect); the new declaration records the old declaration's SHA-256 and the list of changed files, and
  the report says so. A freeze that was withdrawn before any run binds nothing.
* A different design is a new benchmark. A second completed run of any kind is reported next to the first, never
  instead of it.
* `metrics.json` and `evidence.json` are written once. `score` on a finished run is a re-score: it refuses
  predictions whose hash changed, always uses the frozen settings, writes to `rescores/<n>/` and changes nothing
  the API or the report reads. If the code that scores or reports differs from the frozen code, the report says so
  in its first lines.

## 17. Synthetic checks before the freeze

Same code path, same file formats, synthetic cohorts (`synthetic` command; outputs under
`data/interim/synarch_recon/synthetic/`; the builder's runs of the draft code are kept in
`data/interim/synarch_recon/synthetic_builder_draft/`):

* a world where soma, label and profile share latent classes with zero inflation, and a world where they are
  independent, at a small size and at the real size (101 training samples, 4,448 neurons, 49 regions, the real
  label counts, a region that is always zero, all-zero profiles, a duplicated soma, labels with one sample);
* in the independent world every task must end `no_better_than_baseline`;
* unit tests (71, all passing) cover the pairing rule, statuses and lineage, refusals, sealing canaries, folds,
  baselines on a toy set, metric aggregation, equality with `JointModel.infer()`, recovery of planted
  conditionals, the axon adapter, the freeze and the single run, and the three review blockers of section 19
  (one focused test per item, each shown to fail without its fix).

Results of the synthetic runs with the final code (seed 0; primary metric of the engine against the best
comparator, with the relative gain; full tables in each world's `out/results.md`). "supported" is
`exploratory_supported`, "not" is `no_better_than_baseline`:

| synthetic world | T1 (log MSE) | T2 (mean error, um) | T3 (accuracy) | T4 (log MSE) |
|---|---|---|---|---|
| shared classes, small (32 samples, 452 neurons) | supported: 0.451 vs RIDGE 0.478 (+5.63%) | supported: 245 vs KNN 283 (+13.46%) | not: 0.528 vs KNN 0.542 (-2.52%) | supported: 0.393 vs RIDGE 0.426 (+7.77%) |
| independent, small (32 samples, 452 neurons) | not: 1.820 vs MEAN 1.823 (+0.17%) | not: 711 vs RIDGE 712 (+0.18%) | not: 0.312 vs MEAN 0.323 (-3.33%) | not: 1.602 vs MEAN 1.591 (-0.69%) |
| shared classes, real size (101 samples, 4,448 neurons) | supported: 0.570 vs KNN 0.606 (+5.99%) | not: 206 vs KNN 192 (-7.39%) | not: 0.438 vs KNN 0.447 (-1.94%) | not: 0.575 vs KNN 0.584 (+1.56%) |
| independent, real size (101 samples, 4,448 neurons) | not: 1.9339 vs RIDGE 1.9338 (-0.01%) | not: 727 vs RIDGE 726 (-0.10%) | not: 0.1315 vs MEAN 0.1313 (+0.16%) | not: 1.9057 vs RIDGE 1.9054 (-0.02%) |

* In the independent worlds no task is supported. The engine's mean is nominally the best in some cells, by
  fractions of a percent: the minimum gain and the interval rule of section 12 are what prevent a false
  "supported".
* In the shared-classes world at the real size the engine chose (outer folds): T1: hurdle K=16 in 5 folds; T2: hurdle K=12 in 1 fold, hurdle K=16 in 4 folds; T3: hurdle K=12 in 2 folds, hurdle K=16 in 3 folds; T4: hurdle K=16 in 5 folds.
  K sits at or near the model's limit of 16, so the cap of 16 components matters, and a "not" on real data
  may partly be that cap. T4 there shows the status rule at work: the engine is ahead of every comparator, but
  against the best one (KNN) its gain is 1.56% and the interval of the difference includes zero; it is not
  called supported.
* In the independent world at the real size K = 1 was chosen in 5 of 5 folds (T1), 3 of 5 folds (T2), 1 of 5 folds (T3), 5 of 5 folds (T4).
* Run time at the real size with 4 workers on a busy machine: 5 min 16 s for the shared-classes world (187 engine fits, 0 at the iteration limit), 17 min 5 s for the
  independent one (194 engine fits, 27 at the iteration limit: EM converges slowly on noise; these fits are counted, never extended).
  One fit on 3,467 rows, one process: `log` (53 fields) 1.2 / 2.7 / 3.4 / 5.3 / 5.0 / 7.2 s for K = 1 / 2 / 4 / 8 / 12 / 16;
  `hurdle` (102 fields) 4.1 / 4.8 / 4.7 / 7.1 / 12.9 / 22.0 s (`timing` command, the same busy machine).
* The results do not depend on the number of workers: the two small worlds were run with 1 and with 3 workers;
  metrics, every saved prediction array and every selection are identical
  (`data/interim/synarch_recon/synthetic/workers_check.txt`).
* No-signal simulations of the status rule: section 12 (`nullrate_tiny.json`, `nullrate_real.json`).
* These runs used `analysis/synarch_recon.py` `5ad109dd...`, `analysis/synarch_recon_modules.py` `ed0d90fa...`,
  `analysis/synarch_recon_benchmark.py` `c899dd6a...`, `tests/test_synarch_recon.py` `a2c7285d...` (full hashes in
  `synthetic/code_at_run.sha256` and in each synthetic `pre_run_declaration.json`). A synthetic declaration
  does not pin this file (section 15), so filling in this section after the runs changes nothing they bind.

## 18. Decisions taken after the reviews, and what the next reviewer may still change before the freeze

Taken (each is one setting in `CONFIG`, declared and frozen with the rest):

1. MEAN and IDENTITY choose their weighting (brain sample or neuron) on inner folds (section 8).
2. The status needs a relative gain of at least 2% over every comparator, besides the six intervals (section 12).
3. The two tuned references (KNN, RIDGE) stay among the comparators the engine must beat.
4. The `hurdle` representation stays in the engine's grid; it is this layer's addition around the unchanged model
   (section 7). Setting `representations` to `["log"]` removes it (53 fields only).
5. T3 and T4 stay; removing either from `CONFIG["tasks"]` removes it from the run.
6. The coverage flag for zero-inflated targets is one-sided plus a width condition, and is worded as such
   (section 13).
7. The fine soma label is loaded but is not an input anywhere.

Still a choice for the reviewer of this revision: the value of the minimum gain (2%); whether T3 and T4 are worth
their place among four decisions.

## 19. Rules of the reconstruction API that every answer goes through

These rules are code, with one focused regression test per item (`tests/test_synarch_recon.py`, section
"Regression: review blockers K1 K2 K3"); each test fails on a copy of the code with that fix switched off and
passes on the code as it is (`data/interim/synarch_recon/fix_evidence/`).

**K1. What may feed an answer about a neuron whose axon is hidden.** The axon is hidden in every request that is
not one of the three declared profile-input tasks (T2, T3, T4). One rule (`hidden_axon_problem`) is applied to
every record of such a request, whether the answering module uses it or not:

* an inferred record, **whatever its variable name** (a soma position or a source label can be inferred from the
  hidden profile as well as anything else), is admitted only if (i) this layer itself issued exactly that record
  (it is in the layer's ledger of issued values; any change of unit, value, scale statement, producing state or
  lineage makes it another record), (ii) it was issued for this neuron and brain sample, and (iii) for that
  inference and every inference behind it, transitively: the producing state is registered now, with the training
  groups the layer recorded when it issued the value; those groups do not contain the query's brain sample; the
  brain sample is known;
* when the answer itself is axon-derived (a **hidden-axon forecast**: T1, every undeclared conditional with an
  axon-derived target, and the axon adapter), a measured record must be one of: soma position, source label, fine
  soma label, hemisphere, reconstruction type (`reports/axon_dataset_plan.md`, section 8), and so must every
  measured ancestor of an inferred record; the neuron's own measured profile or trace is refused, and so is
  everything that descends from it;
* anything else is refused: missing, unregistered or unverifiable provenance, another status.

In the three declared profile-input tasks the measured profile is an input by declaration, so an inferred input
is flagged there instead of refused: `conditioned_on_inferred_input` always, `inferred_input_unverified` when this
layer did not issue it, `group_in_training_state` when a state behind it holds the unit's brain sample. In every
task, nothing may descend from the requested property itself.

A direct answer by a state that holds other neurons of the unit's brain sample is returned with the notice
`group_in_training_state` and without the coverage flag (a neuron that is itself in the state is refused). Used
as an INPUT of a hidden-axon request, such an answer is refused by the rule above.

Identity: a Jiao neuron's brain sample is the prefix of its ID, and a request that names another one is refused;
identifiers with surrounding blanks are refused; a sample spelled in another case is still that sample. A value
inferred for one neuron cannot be made a record of another. Comparators that copy donors exclude the query's own
brain sample by the training rows' spelling.

The axon adapter declares `axon_arbor` with status **`blocked: independent validation pending`**. Nothing passed
to it unblocks it: not a backend, not an evidence object. It imports no axon engine.

**Limit, stated plainly:** records are not signed. A value that somebody re-enters by hand as a "measured" record,
or a neuron that is presented under a made-up ID, cannot be told from the real thing by any check in this layer.
The guarantee covers measured records as the loaders make them, inferred records exactly as `reconstruct()` issued
them in the running process, and the IDs as given. The ledger is written by `reconstruct()` only and lives in
memory: an inferred value that was saved and reloaded in another process is not admitted where the axon is
hidden; it has to be inferred again.

**K2. Scale and frame.** Every record that holds a value states its scale in one way: a required value tag
(`identity`: the registered measurement scale and units; `log10_1p`: already on the modelling scale) and the one
text registered for that tag; and it states the registered frame of its variable. A missing, unknown or
contradictory statement, other units, another or a missing frame are refused: by the request API, by the planner,
by the record store and by the axon adapter. A value reaches the modelling scale by exactly one conversion, in
one function (`model_values`): 99 um becomes 2.0 once; an estimate of 2.0 made into a record and used as
evidence conditions at 2.0. Every estimate a module returns is held to the same statement before it is issued.
An inferred input is still flagged (`conditioned_on_inferred_input`: its uncertainty is not propagated) and its
answer never carries the coverage flag. An inferred input computed only from measurements that the same request
supplies adds no information and is ignored, and said to be; an inferred value of a zero-inflated region is
refused as evidence under the presence/level encoding (it is a mean over "absent" and "present").

**K3. Identity of a saved state.** A fitted state is read-only, and its fingerprint is computed from exactly
what is saved, in an encoding that cannot be made to collide, and is never cached. It binds every item that can
change an answer, a refusal or a notice:

1. the state format and the layer's schema name;
2. the SHA-256 of the unchanged joint-model source;
3. the reader's settings: the **support margin** (how far outside the training range a continuous input may lie
   before it is refused or marginalised: 5% of the range) and the bisection steps of every interval quantile;
4. the layout: namespace, variables, and every field (name, variable, coordinate, kind, modelling transform) in
   order;
5. the encoding of every continuous field (the representation);
6. per model field: type, supported or not, block, the standardisation constants (training mean and scale), the
   categories and their order;
7. the fit settings: components requested, seed, restarts, iteration limit, tolerance, variance floor, Dirichlet
   prior, representation;
8. per block: name, fields in order, rows, number of components, the convergence record;
9. every parameter array (weights, means, variances, tables): values, shape, type, and which block, kind and
   field it belongs to;
10. the **support ranges** of the continuous fields;
11. the training units and the training groups, as lists.

Loading refuses any file that is not exactly a state as it is saved (one changed item, an unknown or missing
entry, a malformed or missing range, a non-canonical text). And a state is never loaded on the file's own word:
the hash stored in a file only shows that the file agrees with itself, and whoever edits the file can rewrite it.
The caller must name the fingerprint it expects (the one recorded in `selection.json` or `serving.json` when the
state was fitted), and any other state, however consistent, is refused. A module also refuses a state that is no
longer the one it was registered with, and keeps a read-only copy of the comparator rows.

Not part of a state, and so not in its fingerprint: the policy and the interval level (run-time settings,
returned with every result); the task definitions and the benchmark evidence (bound by the pre-run declaration
and by `run_complete.json`, section 15); the comparator bank, whose own hash binds its rows, the comparator that
answered and every setting of it. What file hashes cannot detect: somebody who rewrites every file of a finished
run consistently, `run_complete.json` included. The anchor is the SHA-256 of `run_complete.json`, which `run`
prints and the results report states in its first lines.

## 20. What changed after the reviews (2026-11-11)

Sources: the other agent's review (`reports/synarch_recon_codex_review.md`, three blockers), a full review of the
builder's draft (twelve findings), and three attacks on the builder's first fixes. All changes were developed and
tested on synthetic data.

* K1, K2, K3: completed as in section 19 (the first fixes covered the reviewed cases only).
* Baselines: weighting of MEAN and IDENTITY chosen on inner folds (section 8).
* Status rule: minimum relative gain of 2%; multiplicity stated and measured (section 12).
* Run-once: abort refused once outcomes were used; a restart keeps the design; metrics written once; the files of
  a finished run are bound and verified before anything is served (sections 15, 16).
* API: a declared task answers in its own coordinate order, so names and values cannot be swapped; a
  baseline-level answer is guarded by, and reports, the comparator bank that made it; a reference value attaches
  only to the unit's own group; one unit has one group; a variable is requested once.
* Wording: sections 7 and 13.
