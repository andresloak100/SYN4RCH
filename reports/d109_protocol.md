# D-109 protocol: predicting regional axon length of held-out hypothalamic neurons (Jiao et al. 2025)

Written 2026-11-10, revised the same day after three independent reviews (leakage, protocol against code,
statistics), before any real D-109 outcome is scored. Implementation: `analysis/d109_jiao.py`. Inputs and their
audit: `reports/d109_jiao_audit.md`. Conventions and sources: `reports/d109_jiao_conventions.md`.

**Status is not recorded in this file.** Approval is the line `VERDICT: APPROVED TO RUN` in
`reports/d109_review.md`, which also names the SHA-256 of this file; this file is not edited after review.
Until then no D-109 score exists and `run` refuses (section 14). All development ran on synthetic outcomes or on
outcomes whose whole rows were permuted across neurons with a fixed seed. Section 16 lists what has been seen
or computed on real data outside this protocol.

## 1. Question and scope

Can the regional axon length of a hypothalamic neuron from an animal the model has never seen be predicted
from its soma position, and does the animal's genetic tag improve that prediction?

Not claimed, whatever the result:

* **Not 3D paths.** The outcome is axon length per brain region, a table, not a trajectory.
* **Not connectivity.** Axon length in a region is not a count of synapses or of partners.
* **Not causal.** A tag is a Cre line (or the orexin promoter virus). Every animal carries exactly one tag, so
  the tag stands for everything the animals of a line share: the labelled cell population, but also injection
  target, virus, batch and date. A tag benefit is a statement about prediction only.
* **Not genetic transfer.** A model that has never seen a tag receives no information about it; the unseen-tag
  split measures transfer by position only.
* **Not a bound on the information in the tag.** The primary comparison measures the net gain of one model
  (extra-trees with a tag one-hot). A negative or small result says that this model did not gain, not that the
  tag carries little (section 12).
* **Not gene expression.** No neuron has a measured expression profile.
* Male mice only; one laboratory; soma coordinates as published, which are inferred to be mirrored to one
  hemisphere (section 2); laterality is not available.

## 2. Data

7,180 neurons, 166 animals (ID prefix = fMOST Mouse ID = one brain), 16 tags, one tag per animal.

* Inputs: `All_soma_infor_HY_neurons.csv`, columns `file, geno, location, AP, DV, ML` only; and
  `geno_location_labels.csv`. The columns `type`, `cluster` and `Class` are derived from the axons and are never
  read. `all_neuron_soma_infor.xlsx` is not read by the script.
* Position: AP, DV, ML exactly as published. **Inferred, not stated for this file** (conventions report,
  sections 2 and 3): micrometres in the CCFv3 10 um frame with AP = x, DV = y, ML = z (resolved from the authors'
  code for sibling tables with the same columns; for this file the authors' notebook plots the three columns
  unscaled against their CCF structure meshes); mirrored to one side (documented mirroring rule; all ML below
  5,700; no sentence says this table is the mirrored one). The models use the three numbers as given, so the
  analysis does not depend on these inferences; the description of the result does. Original hemisphere is not
  in the open tables.
* Outcome: `all_neurons_projection_length_fig1.csv`, 49 region columns in the published order.
* The animal is taken to be the ID prefix. No source states that no animal was imaged under two IDs
  (unresolved; inputs only: no two somata of different animals lie within 1 um).
* Li et al. 2024 is not used (section 13).

## 3. Target, units, regions, all-zero profiles

* **Target.** `y[i, r] = log10(1 + L[i, r])`, with `L` the published axon length of neuron `i` in region `r`.
  **Unit: micrometres is indicated, not resolved.** The table carries no unit; the authors' axis labels say
  "(um)" for other length tables, and their package computes lengths as distances in SWC coordinates, but which
  tool made this table is unresolved. The `+ 1` offset therefore **assumes micrometres** (it is 1 um; if the unit
  were millimetres the offset would compress the scale). No normalisation by the neuron's total: total length
  is part of what is predicted.
* **Regions.** The 49 columns exactly as published, equal weight, no re-aggregation. Known and unresolved:
  MY is the CCF ancestor of IO, NTS, PARN and IRN and may count them twice; HPF, OLF and Isocortex are parents
  whose summation rule is not documented; whether columns sum both hemispheres is not stated. Because of this
  the columns are treated as 49 published measurements, not as a partition of the brain. Sensitivity (reported,
  not decisive): the primary comparison **scored on the 48 columns without MY, models unchanged** (the same
  49-output models; nothing is refitted).
* **All-zero profiles.** 13 neurons are zero in all 49 regions; no source explains them. They **stay in training
  and in test** for the primary analysis (their target is the zero vector, and the primary metric is defined for
  it). `run` refuses if their number is not 13. Sensitivity: the primary comparison recomputed without them in
  the test set (they stay in training). They are identified inside `run` by the outcome-only rule "row sum = 0".

## 4. Nucleus vocabulary

One authoritative vocabulary: the **authors' fine label**, the suffix of the published ID (59 CCF acronyms). It
equalled the xlsx `location` column for all neurons when checked by hand on 2026-11-10; the script does not
read the xlsx and does not assert this. It is the authors' assignment; its procedure is not fully documented,
and it is not known whether it was made with knowledge of the reconstruction. No nucleus is derived from
coordinates (no verified annotation volume in the repository).

Crosswalk: fine label -> hold-out key = the fine label, replaced by `MBO` or `ZI` when it lies in that CCF
subtree (`data/raw/atlas/structure_graph_1.json`). The script asserts this equals the authors' 50-label column
for every neuron. The 16-label column (`others`) is not used: `others` is a remainder of 35 labels, not a place.

* As a model input the nucleus is used **only by the nucleus-mean baseline** (section 7). Learned models get
  coordinates, not labels.
* As hold-out key (S2): the 50-label key. `HY` (1,316 somata in hypothalamic voxels without a nucleus) and
  fibre-tract labels (136 somata) are never a hold-out fold; they stay in training when their animal is not purged.

## 5. Splits

Written by `python3 -I analysis/d109_jiao.py splits` to `data/interim/d109_splits.json` from IDs, positions,
labels and tags only; seed 20261110. `run` regenerates the file and refuses on any difference. SHA-256 of the
split file: `606dafac82774a77dc873f3aad65a986b21b5f9d55eb5142d0d865f89d972947` (unchanged by the review fixes).

**S1, tag benefit within represented lines.** Whole-animal 5-fold. Within each tag (sorted by name) the animals
(sorted by ID) are shuffled with a seeded generator and dealt round-robin to folds, with one running offset
across tags. Every tag is in the training set of every fold (asserted; the smallest tags have 2 animals). Every
neuron is tested once. Folds: 1,319 / 1,597 / 1,385 / 1,440 / 1,439 neurons, 34 / 33 / 33 / 33 / 33 animals.
The 166 animals are very unevenly spread over tags: six tags have 16 to 22 animals each (112 animals), three
have 9 to 11 (29), three have 4 to 6 (15), four have 2 or 3 (10). **The primary mean is therefore mostly a
statement about the six large lines.** A test animal of a 2- or 3-animal tag is predicted from one or two
tag-mates, so there the tag arm learns one animal's effect; no per-tag claim is made for tags with fewer than
5 animals. The benefit is printed by tag-size group (2-3, 4-8, 9-15, 16 or more animals), descriptively.

**S2, unseen nucleus label.** One fold per hold-out key that is under HY in the CCF, is not `HY` itself, and has
at least 30 neurons from at least 5 animals: 23 nuclei (AHN ARH AVP DMH LHA LPO MBO MPN MPO PH PMd PMv PSTN PVH
PVHd PVp PeF RCH SBPV STN TU VMH ZI), 5,576 test neurons. Test = every neuron with that key. Training = every
neuron of every animal that has **no** test neuron: all records of any animal with a test neuron are purged,
whatever their nucleus. Training sizes range from 2,462 neurons / 78 animals (LHA) to 6,859 / 160 (STN). 27
labels are not folds (listed by `splits`).
**What S2 is and is not.** It is transfer to a nucleus *label* that no training neuron carries, with every
animal that has a neuron of that label removed. It is **not** transfer to an unsampled region of space: somata
labelled only `HY`, fibre-tract labels and neighbouring or sibling labels (for example PVHd when PVH is held
out) from non-test animals stay in training and can lie next to or inside the held-out nucleus (inputs only:
28 pairs of somata from different animals lie within 10 um of each other). A position model can score well on
S2 from near neighbours with another label. In 8 of the 23 folds the purge also removes every animal of some
tag, so some test animals have a tag absent from training (AHN 5, ARH 5, LHA 6, MPN 2, PVH 8, PeF 2, RCH 2,
SBPV 2 animals); for them a with-tag model is its position-only fallback. The count is printed per fold.

**S3, unseen tag.** Leave one tag out: test = all animals of the tag, training = all others; 16 folds, 2 to 22
animals each. Reported as **position transfer to an unseen line**, never as genetic transfer. With-tag models
are not fitted or printed under S3: by the fallback they would be position-only models at another tuned
setting, and any difference would be tuning noise.

In every fold of every split no animal is on both sides (asserted by `splits`, `check`, `dry` and `run`).

## 6. Preprocessing and tuning

* Every mean, scale, category list and group mean is computed from the training rows of the fit at hand. A
  predictor is a function of (training positions, tags, labels, outcomes; test positions, tags, labels) and
  nothing else. `check` verifies that a test neuron's prediction does not change when other test neurons are
  removed or when test outcomes are corrupted.
* **Nested whole-animal tuning** for every hyperparameter: inside each outer training set, 3 whole-animal folds
  made by the same dealing rule (seeded by split and fold name). Inner loss = the primary metric inside the
  training set: every training animal is an inner test animal exactly once; its error is the mean over its
  neurons; the loss is the mean over all training animals. The value with the smallest loss wins; ties go to the
  first value in the grid. The model is then refitted on the whole outer training set. A tag may be missing from
  an inner training set; the fallbacks apply.
* Grids: k of k-nearest in {5, 10, 20, 40, 80}; ridge lambda in {0.1, 1, 10, 100, 1000}; trees
  `min_samples_leaf` in {3, 10, 30, 100}; network L2 penalty in {1e-4, 1e-2, 1}. (The value 100 for the trees
  was added after review, on synthetic evidence only: larger leaves reduce the penalty the tag one-hot pays for
  fitting training-animal effects.)
* **Grid-edge rule.** The number of folds in which the chosen value is the first or last of its grid is printed
  for every model. Grids are not changed after the outcomes are opened. A wider grid may be run afterwards only
  as a labelled deviation that does not replace the first result.

## 7. Models

Baselines:

1. `global`: training mean profile.
2. `nucleus`: mean profile of training neurons with the same fine label. Back-off when fewer than 5: for a label
   under HY, each CCF ancestor in turn up to HY, pooling all training neurons whose label lies in that subtree;
   then the global mean. `HY` and fibre-tract labels: own label, then global mean.
3. `knn`: mean profile of the k nearest training somata (Euclidean, in the published coordinates). Training
   neurons are always from other animals.
4. `knn_tag`: the same among training neurons of the test neuron's tag; if that tag has fewer than k training
   neurons (or none), `knn`.
5. `ridge`: ridge regression on the 9 standardised terms of a degree-2 polynomial of standardised AP, DV, ML;
   intercept not penalised.
6. `ridge_tag`: plus a centred one-hot of the tag, same penalty.

Learned models, declared in advance, two:

7. `et` / `et_tag`: extremely randomised trees, one multi-output forest for the 49 regions (scikit-learn 1.8.0
   `ExtraTreesRegressor`, 100 trees, `max_features=1.0`, fixed seed), on raw AP, DV, ML, without / with the tag
   one-hot.
8. `mlp` / `mlp_tag`: one multi-output network with two shared hidden layers of 64 units (scikit-learn
   `MLPRegressor`, Adam, batch 128, 60 epochs, no early stopping, fixed seed), on standardised coordinates
   without / with the tag one-hot; outcomes standardised per region on the training set.

Library versions at writing: Python 3.12.4, numpy 1.26.4, scipy 1.13.1, pandas 2.2.2, scikit-learn 1.8.0. They
are written into the run marker and the results; `run` refuses if scikit-learn is not 1.8.0 or if any version
differs from the one `check` passed with.

**Fallback for an unseen tag** (ridge, trees, network): a test neuron whose tag is absent from the training set
is predicted by the position-only model of the same kind and setting, fitted on the same training set.
**Fallback for an unseen nucleus:** the back-off chain of baseline 2; no other model uses labels.

## 8. Metric and the single primary comparison

* Per-neuron error: mean over the 49 regions of `(prediction - y)^2`.
* **Mouse-level aggregation.** S1: the mean error of an animal's test neurons, then the mean over the 166
  animals, each animal weighing the same. S2 and S3: the same within each fold (nucleus or tag), then the mean
  over folds, each fold weighing the same. **One aggregation per split is used everywhere** (skill table and
  comparisons).
* **Primary comparison (one):** split S1, `et_tag` against `et`. Statistic: mean over the 166 animals of
  `d_a = error_a(et) - error_a(et_tag)`; relative benefit = that mean divided by the mean of `error_a(et)`.
  **Estimand:** the net gain of this model from the tag one-hot. Under "animals differ, tags carry nothing" its
  expected value is negative, because the one-hot lets trees fit training-animal effects (reviewer's reduced
  simulation: about -3% to -12%).
* Reported next to the primary, descriptive and never decisive: `ridge_tag` vs `ridge` and `knn_tag` vs `knn`
  on S1 (in simulation the ridge pair realised far more of a planted tag effect than the trees did).
* Also reported per model and split, descriptively: skill = 1 - error / error of `global`; per-neuron Pearson r
  across the 49 regions, aggregated the same way (undefined for a constant observed profile, e.g. all-zero: those
  neurons are left out of r only; a constant prediction scores 0).

## 9. Uncertainty

* Unit = animal (S1), or fold (S2: nucleus, n = 23; S3: tag, n = 16), where a fold's value is the mean over its
  test animals. Per-fold values are listed with their number of animals.
* Intervals: percentile bootstrap over units, 10,000 resamples; 95% for every comparison and also 99% for the
  relative difference. p value: one-sided sign-flip over units, 10,000 random flips, with the +1 correction.
* Effective n: S1 166 animals; S2 23 nuclei (6 to 88 animals each; animals recur across nuclei, so folds are not
  independent); S3 16 tags (2 to 22 animals each). Animals of one training fold share fitted models and animals
  of one tag share whatever the tag stands for, so units are not independent.
* **The intervals are approximate and can be too narrow.** In the statistics reviewer's reduced simulation the
  95% animal bootstrap covered the truth in 92-96% of replicates with independent animal effects, 89% when one
  animal-level offset shifts all regions together (the most plausible real case), and 60% with very
  heavy-tailed animal effects. For this reason "not shown" uses the 99% interval (section 12), and the five
  animals with the largest `|d_a|` are printed together with the primary comparison without them.

## 10. Tag-shuffle control

Tags are permuted **between animals**: each animal keeps all its neurons and receives the tag of another animal.
Primary version: within 4 spatial strata of animals (soma-centroid AP at or below / above the median over
animals, then centroid DV at or below / above the median within each half; 41 or 42 animals each), 200
permutations. Secondary: the same without strata, 100 permutations.

**The whole procedure is repeated on each shuffled assignment**, so the null is built exactly as the observed
statistic: the S1 folds are dealt again from the shuffled tags by the S1 rule (every tag is then in every
training fold, as in the observed arrangement, and no test animal falls back to the position-only model), and
both arms (`et`, whose inner tuning folds depend on the tags, and `et_tag`) are tuned and fitted again. `check`
asserts that the rule applied to the real tags reproduces S1 and the observed statistic.
p = (1 + number of permuted statistics >= observed) / (permutations + 1). Also reported: null mean, null SD,
and (observed - null mean) / null SD.

* What it is: a **one-sided guard, not a calibrated test.** Before the folds were dealt again the reviewers
  found the control valid but conservative (under a simulated null it never rejected, and the observed value
  lay about 1.9 null SDs below the null mean). The cause was not isolated; fixed folds were one candidate, the
  spatial compactness of real tags is another. Dealing the folds again did not remove it: in the one synthetic
  no-tag data set of `check` (39 permutations) the observed value still lay 3.4 null SDs below the null mean,
  so real tags, being spatially compact, let the trees fit training-animal effects more than shuffled tags do.
  The control therefore stays conservative under "animals differ, tags carry nothing"; its error rates have not
  been calibrated. In practice the 2% minimum effect binds before this control does.
* What it can show: that the real assignment of tags to animals predicts new animals better than random
  assignments that keep animals whole and roughly keep where each tag was sampled.
* What it cannot show: that the benefit is genetic or causal. Tag is confounded with animal and with everything
  a line's animals share (injection target, virus, batch). A batch effect aligned with the tag would pass every
  condition, this control included. The strata are coarse (tags are sampled in different places, and finer
  strata would leave almost nothing to permute), so a small p can also mean that the tag carries finer
  positional or targeting information than the position model extracted.
* Inputs-only description of the confound with batch, for which animal-ID order is the only available proxy:
  along the ID order the 166 animals form 57 runs of equal tag (a random assignment gives about 152, minimum 140
  in 2,000 draws; perfect blocking would give 16). Tags are strongly clustered in ID order; several small tags
  are single blocks (Agrp 202812-202816, Trh 212201-212202, Crh, Vip). `splits` prints the table.

## 11. Controls on synthetic outcomes (`check`), and multiplicity

Synthetic outcomes are generated from inputs and a seed only: a smooth random function of position, plus a
random tag-by-region effect, plus a random animal-by-region effect, plus noise.

* Planted tag effect (tag SD 0.7 of noise SD): the primary rule must return SHOWN.
* No tag effect, with position and animal effects: the primary rule must not return SHOWN.
* Moderate planted tag effect (tag SD 0.3): the benefit realised by `et_tag` and by `ridge_tag` is printed.
  **Not a gate.**
* Leak detector: outcomes that depend on the animal only; in S1, S2 and S3 no model may beat `global` by more
  than 1%. Positive controls of the detector: S2 without the animal purge, and neuron-level folds
  (animals on both sides), must each break that 1% criterion. The detector is weak by nature: somata of
  different animals are spatially mixed, so even neuron-level folds give position models only about 3% skill on
  animal-only outcomes. The guarantee against leakage is the set of assertions on the split file (section 5);
  the detector is a supplementary test.
* Outcome rows permuted across neurons: no skill, no tag benefit.

`check` uses 39 shuffles per planted condition (smallest possible p 0.025) and one data set per condition: it is
a functional test, not a calibration of error rates, and **says nothing about power**. Reviewer's estimate from
a reduced simulation (12 regions, 30 trees, no nested tuning; not the full pipeline): the rule has high power
only when the realised benefit is near +3% or more, and `et_tag` realised about +0.2% of a planted effect worth
+5% to an oracle (tag SD 0.3) and +3.6% of one worth +13% (tag SD 0.5). The full pipeline agrees on the one
moderate data set of `check` (tag SD 0.3): `et_tag` realised +0.08% (95% interval -0.33% to +0.49%) where
`ridge_tag` realised +3.83% (+3.20% to +4.44%). A moderate real tag effect can therefore end as "not shown for
this model" or "insufficient evidence". **The primary instrument is known in advance to be insensitive to
moderate tag effects; whether to keep it as primary is an open decision for the supervisor before the freeze
(recorded in the review).**

Multiplicity: one primary comparison, no correction. Secondary family, Holm at 0.05 over seven sign-flip p values:
H0 S1 `et` vs `global`; H1 S1 `et` vs `nucleus`; H2 S1 `et` vs `knn`; H3 S1 `mlp_tag` vs `mlp`; H4 S2 `et` vs
`global`; H5 S2 `et_tag` vs `et`; H6 S3 `et` vs `global`. A secondary comparison is claimed only if its adjusted
p is below 0.05 **and** its relative difference is at least 2%. H1 and H2 are claimed only if H0 is. **H3 and H5
are tag comparisons without a shuffle control: they are claimed only if the primary result is SHOWN, and then
only as consistency of the primary finding; they can never support a tag-benefit statement by themselves.** H5
partly measures the unseen-tag fallback (section 5). Reason for the 2% floor, seen on permuted outcomes: with no
signal at all, a model that overfits less beats one that overfits more by a fraction of a percent in almost
every animal, and the sign-flip p is then tiny. A p value between two models is therefore not evidence that
either predicts anything. Everything else is descriptive.

## 12. Reading the primary result

Let `rel` be the relative benefit, `[lo, hi]` its 95% interval and `hi99` the upper end of its 99% interval.

* **Shown:** `rel` >= 2%, the 95% interval of the difference excludes 0, sign-flip p < 0.025, stratified
  tag-shuffle p < 0.05, and `check` passed. Wording: "the tag improves prediction of regional axon length for
  new animals of represented lines"; predictive, not causal; mostly about the six large lines.
* **Not shown for this model:** `hi99` < 2%. Wording: "this extra-trees model did not gain 2% from the tag
  one-hot". **This does not bound the information in the tag**: the model realises only part of a real tag
  effect and pays a penalty under no effect (sections 8 and 11). The descriptive `ridge_tag` and `knn_tag`
  comparisons are reported next to it and are not decisive either way.
* **Insufficient evidence:** anything else.

The 2% threshold is a declared minimum of practical interest, not derived from data.

Declared sensitivity analyses of the primary comparison (reported, never decisive): without the all-zero test
profiles; scored on 48 columns without MY; without the five animals with the largest `|d_a|`; without the 32
animals of the SWC sample (section 16); and the breakdown by tag-size group.

## 13. Li et al. 2024

Not used in D-109. Verified at the animal level: Li and Jiao share brains 195037 and 195038 (106 Jiao neurons);
"106" is that derived count, no published statement of overlap was found, and neuron-level identity is unresolved.
Any later use must drop both brains whole; only brain 18696 is a different animal, same laboratory, one animal.

## 14. Freeze, review, run

Order: (1) independent review of this protocol, the script and the split file, recorded in
`reports/d109_review.md`. To approve, the review must contain, each on its own line, exactly once:

    VERDICT: APPROVED TO RUN
    REVIEWED SCRIPT SHA-256 (PINNED block blanked): <64 hex>
    REVIEWED PROTOCOL SHA-256: <64 hex>
    REVIEWED SPLITS SHA-256: <64 hex>

(`hashes` prints the three lines). The script hash is that of the file with the contents of its `PINNED` block
removed, so pasting `PINNED` is the only change to the script that the run accepts after review; any other
change to the script, this protocol or the split rule needs a new review. (2) `hashes`, and its `PINNED` output
(SHA-256 of the three input tables, the CCF structure graph, this protocol, the split file and the review)
pasted into the script. (3) Commit the script, this protocol, the conventions report and the review. (4) `run`,
once. `check` need not be repeated after the paste: its record carries the blanked script hash.

`run` refuses unless: `PINNED` has all seven keys and every hash matches; the review has the verdict line and
its three hashes equal the script (blanked), this protocol and the split file; the script, this protocol, the
conventions report and the review are the committed versions with no uncommitted change; scikit-learn is
1.8.0; the split file is byte-identical to a fresh `splits`; `check` has passed with this script (blanked
hash), split file and library versions; the number of all-zero profiles is 13; and
`data/interim/d109_run_opened.json` does not exist. That marker is written before the outcome table is read,
is printed in the run output and stored in the results file, and must be copied into the committed results
report. The code is frozen by the commit recorded in it.

**Reruns.** The first completed run is the result. `D109_RERUN=1` is only for a crash or a bug; it is reported
as a deviation, writes to new file names (`d109_run_repeat<n>_*`) and never overwrites; `run` also refuses if a
results file of the same name exists.

**Limits of the guard, stated plainly.** (a) It is in-process and protects the scripted path only: the raw
table is an ordinary file, and code that imports the module could call the reader directly. No other code may
import `analysis/d109_jiao.py`, read the table, or set `_RUN_OPEN` before the run; reviewers search the
repository for `projection_length`, `_read_outcome_table` and `_RUN_OPEN` just before `run`. (b) The `check`
record and the run marker live in `data/interim/`, which is not version controlled and can be edited or deleted
by hand; "once" therefore holds on the honest path only, and the copy of the marker in the committed results
report is what makes a silent repeat visible. (c) `data/interim/d109_splits.json`, the `check` record, the logs
and the copies of the source documents are likewise on this machine only. They are reproducible: the split
file from the script (hash above), the documents by `sources --fetch` (19 hashes in the script). Whether to
force-add the small split file is the owner's decision.

Changing anything above after the outcomes are opened is a deviation and must be reported as one.

## 15. Open points that limit interpretation

Listed in `reports/d109_jiao_conventions.md`, section 9. The ones that touch this protocol: the length unit
(micrometres indicated, not resolved); axes, units and mirroring of the soma table (inferred for this file);
both-hemisphere sums not stated; MY / HPF / OLF / Isocortex summation; the 13 all-zero profiles; the label
assignment procedure and whether it used the reconstruction; the tool that made the table; no original
hemisphere; animal = ID prefix not stated in so many words. An outcome-only comparison with the 385-column
table would settle three of them; it has not been run and needs the supervisor's approval. These are to be
repeated next to the S1 and S2 results.

## 16. Disclosures: what was seen or computed before the freeze

1. **D-110 and D-111 (separate decisions, another outcome).** Before this protocol was frozen, 37 Jiao SWC
   reconstructions (32 animals; `data/raw/jiao2025/swc/MANIFEST.tsv`) were downloaded and viewed, and D-110
   and D-111 scored soma-conditioned predictions of **3D axon geometry** on 16 of them (14 animals:
   201230_079, 201232_040, 201232_056, 201894_001, 202170_016, 202559_043, 202617_025, 202620_016, 202824_062,
   210007_142, 210104_007, 210796_002, 211353_070, 211889_096, 212202_011, 212202_026), with whole-animal
   hold-out and no tag (`reports/d110_results.md`, `reports/d111_real_results.md`). That is a statistic
   relating soma position to real axon data of D-109 neurons. It concerns 3D paths, not the 49-region length
   table, it uses no tag and it does not touch the primary comparison. No D-109 design choice used it; the
   writer read the two result summaries only to write this paragraph, after the design was fixed, and the
   changes made afterwards are those the reviewers asked for. Whether this is acceptable under the rule "no
   statistic relating an input to the real axon outcomes before review" is the supervisor's decision.
   Declared sensitivity: the primary comparison without the 32 animals of the 37 neurons (134 animals).
2. **Outcome-only numbers seen.** PERMUTED printouts contain marginals of the real table that involve no input
   (error of the global mean, mean correlation with the global mean profile). The committed audit holds
   input-only and outcome-only marginals. Published summary tables (Nature Supplementary Tables 3, 4, 6) were
   displayed once by a verifier; nothing was computed from them.
3. **Protocol rules added after seeing PERMUTED or SYNTHETIC output, never real input-outcome relations:** the
   2% floor and the H0 prerequisite of the secondary family (permuted output); and, after the reviews
   (synthetic simulations): tree leaf size 100, the 99% interval for "not shown" and its wording, the gate on
   H3 / H5, folds dealt again inside the tag shuffle, the grid-edge and rerun rules, pooled inner loss, fold
   aggregation for S2 / S3, no with-tag rows under S3, the declared sensitivity analyses.
4. For the Li 2024 check the first 500 bytes (soma line) of 6 Li SWC files were read; no axon data.
5. **A separate analysis ("BRM", `analysis/brm_real.py`, `reports/brm_real_results.md`; not part of D-109, by
   another author) scored real predictions of this same 49-region table on 2026-11-10**, after the D-109 design
   and its review fixes were written and while its final `check` was running (its marker
   `data/interim/brm_real/RUN_STARTED.json` is dated 19:28 local time; the D-109 script was last changed at
   18:36). Its protocol says it uses the same 7,180 neurons, soma position plus the source label, and the D-109
   S1 split file. So a statistic relating soma position and tag to the real D-109 outcomes on the D-109 primary
   split **exists before the D-109 review is complete**, against the rule this protocol was written under. The
   D-109 writer read only the first sentence of that report (its model had lower held-out error than copying
   rules and simple averages) and none of its numbers; no D-109 choice used it. Consequence: from that moment
   D-109 is no longer the first look at these outcomes. Anyone who has read the BRM results must not change this
   protocol; whether D-109 can still be called pre-registered, and whether it needs another reviewer who has not
   seen them, is the supervisor's decision.
