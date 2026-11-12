# Axon completion v2: one error-informed development follow-up

STATUS: FINAL — one fixed error-informed TRAIN-development follow-up.
ROOT_CLEARANCE: APPROVED — run only after the cache-coverage plan passes and source/tests/protocol are committed.

Owned new files: `analysis/axon_completion_v2.py`, `tests/test_axon_completion_v2.py`, this protocol.
Outputs: `data/interim/axon_completion_v2/`. No v1, contract, frozen result or raw file is edited.

This is **error-informed DEVELOPMENT on spent TRAIN**, not independent confirmation.
The completed v1 result is preserved (SHA-256
`bc76687b0bdd99063dca2aa82a089bf2522cec21afdd268c8fed5bd5773ba428`).
Its 300-um synthetic-masking completion task is unchanged; there is no whole-neuron,
gene-to-circuit, synaptic or viewer-promotion claim. Old validation/test and external
MouseLight traces remain closed.

## Fixed hypothesis and unchanged components

The v1 coarse donor summaries may omit route identity, and its dual-head rule did not
constrain cable error. The ONE v2 candidate adds a fit-only donor geometry representation
and a third cable-risk head. This tests a bundle; no gain will be attributed separately
to either change. There is no architecture, landmark-count, target or threshold sweep.

Reuse the exact v1 run frozen at `b9ecc2b`, its `RUN_STARTED.json` source pins and completed
receipts. The following stay identical:

* Observation300um, metric50um/200um, original coordinates, append-only observation.
* Shared frontier-graft decoder including its shared-trunk deduplication and donor-cable guard.
* Eligibility, donor rules, 16-soma/8-absolute/8-centred candidate union and all controls.
* All5outer/4inner group memberships, copied from the pinned `folds.json`.
* Exact per-fold fallback and metadata-fallback choices from pinned v1 `fold_info.json`.
  No inner refit or fallback reselection is necessary for this fixed follow-up.
* Sample-balanced row weights and fixed HGB settings:200iterations, learning rate0.05,
  max31leaves, min50rows/leaf, random seed20261113. The same configuration is used for each head.

## Fit-only route representation

For EACH outer fitting pool, use the independently verified native D through the frozen
FitDoor, with both endpoint groups checked. No query hidden geometry or query-to-landmark
distance is a prediction input.

1. Compute the sample-balanced fit medoid: each fitting sample has equal total mass,
   neurons inside a sample equal mass; minimize weighted D, ties by donor ID.
2. Start the landmark bank with that medoid. Repeatedly maximize the minimum D to selected
   landmarks, among neurons from samples not yet represented. Ties use ascending donor ID.
3. Choose32landmarks; cap only when there are fewer than32fitting samples.
4. Append donor `log1p(D(d,landmark)/200)` and donor-minus-fallback values to the unchanged
   v1 prefix Featurizer:64new columns for the actual outer pools.
5. For a fitting query, set BOTH columns for any landmark from the query's own brain
   sample to NaN. HGB handles missing values. This prevents donor-to-own-query landmark
   geometry from being a fitting shortcut. Outer-eval query samples are outside the fit
   landmark bank, so their landmark columns are unmasked.

All vocabularies, centrality, landmarks and donor-feature reads stay in the current fit
animals. Landmark IDs are provenance of a learned representation, never query-ID features.
Only their donor distances become model inputs. The cached F matrix is never used.

## Three fixed targets and decision

Fit only on VALID fitting labels from the immutable v1 `labels_fit.jsonl`, through LabelDoor:

* `gd = log1p(distance_fb/200) - log1p(distance_d/200)` (unchanged).
* `gf = F1_d - F1_fb` (unchanged).
* `gc = abs(extra_fb-hidden)/hidden - abs(extra_d-hidden)/hidden`.

Missing/nonfinite cable fields in a supposedly valid fitting label are a hard error, not a
new selective exclusion. Invalid labels retain v1's counted training-row handling.

Switch only when candidate differs from fallback and predicted `gd>0`, `gf>0`, `gc>=0`.
Among eligible candidates maximize `gd+gf`, ID tie. Otherwise retain the EXACT fallback.
Nonfinite predictions fall back. The observation and every emitted geometry operator are
unchanged, so label meaning remains identical.

## Controls and gate

Every required v1 control remains. Add the frozen v1 `completion_selector` as the required
`completion_selector_v1` control. Fixed controls, metadata-only and diagnostic observation-only
reuse their identical v1 forecast choices and hashes.

Recreate `switch_random` on the v2 selector's switched queries, using the same seed recipe.
Recreate `permuted_labels` by fitting the same new three-head pipeline with one shared
permutation applied to all three fitting targets. No old switch sets are reused.

Before fitting each outer fold, reconstruct its evaluation locals and assert every v1
fallback-set donor exactly matches the frozen forecast. Assert that old selector,
metadata-only, random and permuted donors belong to the same candidate union.
Full-pool copying controls may legitimately lie outside that union.

The v1 gate is unchanged and also applies against `completion_selector_v1`:

* Zero invalid, empty, duplicate or cable-guard-failed candidate forecasts.
* Against EVERY required control: distance delta<=-10um and F1 delta>=+0.010;
  paired sample95%intervals must respectively exclude0 below/above.
* An invalid/unevaluable required control makes the joint result FAIL/inconclusive.
* Inflated share (`extra>2*hidden`) <= every required control's share+0.02.
* Sample-equal mean absolute relative cable error no higher than adaptive fallback.
* Positive switch rate. Observation-only stays a reported failure diagnostic, not a required control.

Any failure is reported as failure. A development pass still requires fresh independent data.

## Cache reuse, sealing and runtime

All reused v1 files are pinned explicitly in the new source, and original source pins are
checked against v1 `RUN_STARTED.json`. v1 results and receipts are checked before/after the run.

Fitting loads only `labels_fit.jsonl`; per-fold doors enforce both endpoints. A label may
exist globally because a neuron fits another fold, but the current fit cannot read its
outer-eval labels. Each outer fold performs a refused outer-label probe.

Save per-fold choices/checkpoints, then decode the new selector/random/permuted arms from
public observations. Save observation SHA and emitted completion SHA for every method and
seal `forecasts.json` BEFORE admitting `labels_score.jsonl` into the evaluator.
The decode worker, like v1, loads the canonical query arbor for acquisition and calls
`observe`, then passes ONLY the public prefix to the decoder. This is a code-guarded
acquisition handoff, not physical absence of hidden query geometry on disk.

After the seal, reuse the exact original fit and score caches. Validate each admitted
selected score's completion hash against the new emitted hash. Compute ONLY query/donor
pairs missing from both original caches, and only after sealing. No185k-label rebuild,
no9.5million-pair audit, no full-matrix serialization. A bounded64-donor LRU inside each
worker avoids repeated donor loading/cutting without changing the decoder.

`plan` does no fitting or biological geometry/score computation. It checks fitting-cache
coverage and gives a conservative13,104new-evaluation-pair upper bound (three new arms
times4,368queries); exact missing coverage is recorded only after choices are sealed.

Run continuously with <=4workers and all four BLAS thread variables1 using
`/opt/anaconda3/bin/python -I`. Reuse same-code checkpoints if interrupted; changed code,
protocol or inputs refuse resume. A completed/failed run is never repeated under this version.

## Before fitting

Focused synthetic tests must verify deterministic fit-only landmarks, forbidden held-out
landmark reads, hidden-outer-cache invariance, gc sign/negative-risk veto, exact fallback,
valid-label row weighting, the required frozen-v1 comparison and post-seal cache admission.
One fixed-HGB toy problem must recover planted routes from permitted prefix and donor
landmark features on12new toy queries and reject attractive negative-cable-gain donors.
This is software-capability validation only, not biological accuracy.

ROOT reviews the concrete protocol/tests. Only after explicit clearance is this document
marked FINAL/APPROVED and committed with source/tests. The executable refuses biological
fitting while this document remains DRAFT/PENDING.
