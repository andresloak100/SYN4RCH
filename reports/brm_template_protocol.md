# Conditional template selection with local soma attachment

Declared 2026-11-10 in response to the request to diagnose and improve the
failed regional-to-3D bridge. Existing biological runs remain frozen. The
replacement is built and tested **synthetically only** while Agent 1 prepares
a larger, independently split geometry dataset. The previously exposed 16
neurons may inform diagnostic calculations, not another model-performance
claim. Oracle selectors using their hidden measurements are diagnostic only.

## Question and scope

Can a conditional distribution over measured training arbors select a useful
whole template without damaging its distal anatomical coordinates? This
replaces the three-axis RMS decoder and nearest-regional-mean selector. It
does not alter the joint regional probability model, invent neuron-specific
genes, average unmatched branches or create a synaptic circuit map.

The existing regional biological result is already qualified: animal-held-out
log-MSE and profile correlation improved, but MAE did not. Regional length
columns do not supply a sequence of traversed structures or an ordered path.
Neither regional success nor a synthetic geometry success proves accurate
biological pathway reconstruction.

## Conditional distribution over donors

Each eligible measured donor j has an animal ID, source label, measured
regional profile and canonical 50-micrometre arbor. Query inputs contain only
known soma position, source label and identity/provenance. The query's true
profile, path, length, branch count and RMS are not inference arguments.

Use one query-safe regional state for the query posterior pi and all donor
full-evidence responsibilities r. The query and its entire animal must be
absent from that state's regional fitting and all geometry donor records.
Reuse the checked frozen posterior arithmetic; do not silently substitute a
donor's different out-of-fold state.

Define a_j = 1 / number of eligible donor neurons from animal(j), A = number
of donor animals, m_k = sum_j a_j r_jk, with fixed shrinkage lambda = 1:

```
v_jk = (a_j r_jk + lambda a_j/A) / (m_k + lambda)
w_j(q) = sum_k pi_qk v_jk
```

Each v column and w sum to one. The prior is the animal-balanced distribution
over whole arbors, not a mean of geometrically unmatched points. A query
with no useful conditional evidence should not force a novel deformation.
Record weights, component support and source provenance.

For uniform guidance use w_j = a_j/A. For conditional guidance use the above
w. Select the measured template i minimizing sum_j w_j D(i,j), where D is
the unchanged symmetric arc-weighted 3D nearest-sample distance between
soma-centred training arbors. This is a discrete geometric-risk decision,
not Euclidean closeness to a 49-dimensional regional average. Record the
corresponding weighted F1 risk as a diagnostic/control; it is not a hidden
query F1. Ties follow ascending neuron ID.

## Five policies and anatomical attachment

Keep exactly five policies:

1. Existing animal-balanced geometric medoid, translated as in D-110.
2. Uniform geometric medoid with a 500-micrometre local attachment.
3. Conditional geometric medoid with that 500-micrometre attachment.
4. Uniform geometric medoid with a 1,000-micrometre local attachment.
5. Conditional geometric medoid with that 1,000-micrometre attachment.

Policy 1 is an exact fallback, including its canonical coordinates and
selection rule. Its global translation can move distal targets, so matching
uniform locally attached controls separate anatomical anchoring from useful
regional conditioning. No broad hyperparameter search is planned.

For a chosen template let d_i be its original root-to-node path length,
delta = query soma minus donor soma, t_i = min(d_i/L, 1), and
u_i = 1 - 3 t_i^2 + 2 t_i^3. Move node i to x_i + delta u_i. The soma moves
to the query; EVERY existing vertex with d_i >= L stays exactly at its
measured atlas coordinate. Preserve all original branches and vertices,
subdividing enlarged edges to at most 50 micrometres without chord smoothing.

The derivative magnitude is at most 1.5/L along the tree. Require
norm(delta) <= L/3, giving a displacement difference no greater than half an
original edge's length. Verify each positive-length edge's length ratio is
within [0.5, 1.5] and its dot product with the original edge is positive.
Unsupported or distorted attachments return the EXACT baseline for that
query, with an explicit reason; no query may be omitted from evaluation.
These checks bound local edge distortion; they do not prove absence of
global self-intersections, valid obstacles or correct anatomical routes.

## Validation and outcome boundaries

Policy selection uses a separate whole-animal validation set only. A
nonbaseline policy must have strictly lower animal-mean distance AND higher
animal-mean F1 than the exact baseline on the complete declared validation
set, including per-query fallbacks. Among qualifying policies choose lowest
distance, then highest F1, then fixed policy order. Otherwise retain baseline.
This is a selection rule, not a statistical guarantee on unseen animals.

Conditional policies additionally require improvement over the uniform
attached control at the same radius on both validation measures before they
are eligible. A test claim specifically about joint guidance also requires
improvement over that matched control on both held-out measures. An
attachment-only improvement must be labelled as such. Record selection and
per-query fallback reasons. Never choose per-query templates from true
query geometry or promote an oracle to a deployable model.

Regional and geometry fitting animals, validation animals and test animals
must be disjoint where each outcome is withheld. Regional states used for
validation must exclude ALL validation AND final-test animals. An outer
query-safe state can still contain an inner-validation animal's outcomes;
the old five S1 states alone do not validate nested geometry tuning. Refuse
invalid membership/proofs rather than quietly falling back. Any future
nested tuning needs upstream states excluding both levels, or a single
declared train/validation/test partition with regional fits from training
animals only. Agent 1's larger dataset must resolve this before real runs.

## Synthetic verification and reporting

New files only: `analysis/brm_template.py`, its focused tests, a synthetic
benchmark, and new protocol/review/results files. No import may fit a model
or read biological data. Pin the unchanged joint, regional and geometry
dependencies. Synthetic artifacts live in a separate
`data/interim/brm_template_synthetic/` bundle with source/protocol/seed hashes.

Before the benchmark, independently review the implementation and fixtures.
Check normalized animal-balanced donor weights, exact baseline identity,
component and donor permutation behavior, duplicate-animal invariance,
query-outcome omission, unseen tags, distal anchor preservation, edge support
bounds, topology/vertex preservation and rejection of unsafe states/splits.

The user's additional required positive control has two distinct measured
synthetic donors and a known appropriate donor B. The query's permitted
label, through the synthetic regional state and training-only geometry,
must select B. Regenerating the forecast after changing the external hidden
query axon and true regional profile must preserve the donor, weights and
coordinates bit for bit. External outcome loaders are forbidden during
this test. Test geometry and truth-based oracle selections are never inputs
or candidate policies of the replacement.

Use three declared synthetic worlds with independent train, validation and
test animals: informative molecular/source labels connected to distinct
routes; no observable route information; and a misleading association that
fails in validation. Fit the unchanged joint model on synthetic training
rows only. Do not read or import any old real model state for these worlds.
Declare generation, seeds, counts and pass conditions in the benchmark
before running; retain all five policies, validation decisions and held-out
scores. A fixture with informative labels should exercise conditional
benefit; uninformative/misleading validation should preserve fallback. Do
not tune generator difficulty from held-out benchmark performance.

Save selected policy and every test forecast before opening test truth for
scoring. Report all five methods, selected-policy and baseline outcomes,
animal-mean distance/F1, query support/fallbacks and split/state hashes.
An independent replay must verify predictions, validation selection and
metrics from saved outputs. A passing synthetic fixture demonstrates the
implementation's ability under its declared assumptions, not biological
improvement. The next biological test requires Agent 1's new dataset,
verified animal independence, a separately frozen protocol and one test run.
The reviewed algorithm, candidate family and validation rule must be
hash-frozen after synthetic validation and independent audit, before fitting
or evaluating the replacement on Agent 1's biological cohort. Its training
parameters and validation-selected policy must then be pinned before the
independent test outcomes are opened.
