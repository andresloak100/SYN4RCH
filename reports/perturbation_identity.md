# SYN4RCH perturbation identity report

**Under `lab` grouping, 33.9% of placed records have a perturbation identity that occurs in exactly one group: when that group is held out, nothing in training shares the manipulation, so a fold-averaged accuracy over them measures the label prior and not generalisation.**

Every number below is counted from the input named in the Run table. Nothing is estimated, imputed or projected, and no target is inferred from a targeting strategy. SCHEMA.md 5b: a cross-validated number must be reported alongside the fraction of test rows whose perturbation identity was unseen in training. This report is that fraction.

## Run

| field | value |
|---|---|
| generated (UTC) | 2026-11-06T17:34:32Z |
| input | `data/curated/records.jsonl` |
| input present | yes |
| scope | all records read |
| records analysed | 508 |
| perturbation identity | `perturbation_class / perturbation_polarity / perturbation_target` |
| target resolution state | targets partially resolved |
| records with a resolved `target_id` | 170 |
| records with `target_raw` only | 0 |
| records whose class needs a target | 175 (of which 170 have one) |
| load errors | 0 |

> **What state the data was in.** Targets partially resolved. A molecular record's identity is its target gene; with the target null, every knockout in the corpus shares one identity and every overexpression another, so such records are counted as `target_unresolved` and are **excluded from both** the learnable and the unlearnable bucket rather than counted as learnable. Reading them as learnable is the one serious misinterpretation available here, and it is the artefact of a field that has not been filled in yet.

## What a perturbation identity is, and when it is missing

An identity is `(perturbation_class, perturbation_polarity, perturbation_target)` — `syn4rch.qc.keys.perturbation_key`. The three travel together because a polarity means different physical things in different classes: `decrease` is less firing under `circuit_activity` and less of a gene product under `molecular`, and nothing in this project is allowed to key on polarity alone.

A class whose identity *is* its target (`molecular`, `pharmacological`) and whose target is null has no identity to place. Those records are counted on their own line, never folded into either bucket:

- `resolved` — as specific as the class allows.
- `target_unresolved` — the class needs a target and there is none. For `molecular` that is a missing field, being filled in now. For `pharmacological` the identity is the *drug*, and `perturbation_key` has no ligand slot at all, so those records cannot reach `resolved` today and this is a limitation of the key, not of the curation.
- `unresolvable` — no `method`, so no class and no axis.

For `circuit_activity` and `structural` a null target is correct, not missing: the manipulation acts on the population, and the population is a separate component of the claim key.

## Identity reach, per declared grouping

Groupings are `syn4rch.models.task.GROUPINGS`, reused rather than restated: SCHEMA.md 5 makes that registry the single place a grouping name becomes legitimate. No grouping is omitted, including the leakage-prone one, because a grouping that produces an inconvenient number is the one most worth printing.

*Placed* records are those with a resolved identity and a group. *Single-group* records are those whose identity occurs in exactly one group — under leave-one-group-out, the fold that tests them is the fold that removed every record sharing their manipulation.

| grouping | leakage-prone | groups | placed records | single-group records | unseen-identity fraction | identity not placeable | ungroupable |
|---|---|---|---|---|---|---|---|
| `lab` | no | 68 | 502 | 170 | 33.9% | 6 | 0 |
| `population` | no | 44 | 502 | 122 | 24.3% | 6 | 0 |
| `population_outcome` | no | 182 | 502 | 0 | 0.0% | 6 | 0 |
| `record` | yes | 508 | 502 | 0 | 0.0% | 6 | 0 |
| `study_id` | no | 75 | 502 | 170 | 33.9% | 6 | 0 |

Every fraction above is stated with its grouping, as SCHEMA.md 5 requires. `IdentityCoverage` has no `__float__`, so the bare number cannot be formatted without naming the field that carries the grouping beside it.

## By perturbation class — where the contrast is

This is the breakdown the analysis exists for. `circuit_activity` and `molecular` are not two flavours of the same thing; they behave completely differently under cross-study grouping, and the reason is what each class's identity *is*.

### Grouping `lab`

| class | records | identity resolved | distinct identities | single-group | multi-group | unseen-identity fraction | not placeable |
|---|---|---|---|---|---|---|---|
| `<null>` | 1 | 0 | 0 | 0 | 0 | not measurable | 1 |
| `circuit_activity` | 314 | 314 | 3 | 0 | 314 | 0.0% | 0 |
| `molecular` | 170 | 170 | 22 | 170 | 0 | 100.0% | 0 |
| `pharmacological` | 5 | 0 | 0 | 0 | 0 | not measurable | 5 |
| `structural` | 18 | 18 | 1 | 0 | 18 | 0.0% | 0 |

### Grouping `population`

| class | records | identity resolved | distinct identities | single-group | multi-group | unseen-identity fraction | not placeable |
|---|---|---|---|---|---|---|---|
| `<null>` | 1 | 0 | 0 | 0 | 0 | not measurable | 1 |
| `circuit_activity` | 314 | 314 | 3 | 0 | 314 | 0.0% | 0 |
| `molecular` | 170 | 170 | 22 | 122 | 48 | 71.8% | 0 |
| `pharmacological` | 5 | 0 | 0 | 0 | 0 | not measurable | 5 |
| `structural` | 18 | 18 | 1 | 0 | 18 | 0.0% | 0 |

### Grouping `population_outcome`

| class | records | identity resolved | distinct identities | single-group | multi-group | unseen-identity fraction | not placeable |
|---|---|---|---|---|---|---|---|
| `<null>` | 1 | 0 | 0 | 0 | 0 | not measurable | 1 |
| `circuit_activity` | 314 | 314 | 3 | 0 | 314 | 0.0% | 0 |
| `molecular` | 170 | 170 | 22 | 0 | 170 | 0.0% | 0 |
| `pharmacological` | 5 | 0 | 0 | 0 | 0 | not measurable | 5 |
| `structural` | 18 | 18 | 1 | 0 | 18 | 0.0% | 0 |

### Grouping `study_id`

| class | records | identity resolved | distinct identities | single-group | multi-group | unseen-identity fraction | not placeable |
|---|---|---|---|---|---|---|---|
| `<null>` | 1 | 0 | 0 | 0 | 0 | not measurable | 1 |
| `circuit_activity` | 314 | 314 | 3 | 0 | 314 | 0.0% | 0 |
| `molecular` | 170 | 170 | 22 | 170 | 0 | 100.0% | 0 |
| `pharmacological` | 5 | 0 | 0 | 0 | 0 | not measurable | 5 |
| `structural` | 18 | 18 | 1 | 0 | 18 | 0.0% | 0 |

## The same contrast at cell level — and why this is the honest one

A *cell* is `(population, perturbation identity, outcome)`: the claim key without its condition slots. The identity tables above ask whether a **perturbation identity** recurs across groups, and for `circuit_activity` that question is close to vacuous. `perturbation_key` is `(class, polarity, target)` and a circuit record has no target, so the whole circuit half of this corpus holds a handful of identities and none of them is confined to one group. That is a property of a coarse key, not of the literature.

What actually tells one circuit manipulation from another is the **population** it was applied to, and the population lives in the cell. So the question a held-out record cares about — *is there anything specific here to have learned from?* — is the one below. See `docs/DECISIONS.md` D-048, which records this correction and the target it fixes.

One table below is 0 by construction and is printed anyway. A cell contains the population, so a cell can never span two `population` groups, and the same holds for `population_outcome`. Those rows measure nothing; dropping them would be tidier and would also be the habit of dropping the grouping whose number is inconvenient, so every declared grouping is shown and this sentence says which to ignore.

### Grouping `lab`

| class | records | cells | cells spanning >=2 groups | share |
|---|---|---|---|---|
| `circuit_activity` | 314 | 156 | 26 | 17% |
| `molecular` | 170 | 132 | 0 | 0% |
| `pharmacological` | 5 | 3 | 0 | 0% |
| `structural` | 18 | 13 | 0 | 0% |
| `unclassified` | 1 | 1 | 0 | 0% |

### Grouping `population`

| class | records | cells | cells spanning >=2 groups | share |
|---|---|---|---|---|
| `circuit_activity` | 314 | 156 | 0 | 0% |
| `molecular` | 170 | 132 | 0 | 0% |
| `pharmacological` | 5 | 3 | 0 | 0% |
| `structural` | 18 | 13 | 0 | 0% |
| `unclassified` | 1 | 1 | 0 | 0% |

### Grouping `population_outcome`

| class | records | cells | cells spanning >=2 groups | share |
|---|---|---|---|---|
| `circuit_activity` | 314 | 156 | 0 | 0% |
| `molecular` | 170 | 132 | 0 | 0% |
| `pharmacological` | 5 | 3 | 0 | 0% |
| `structural` | 18 | 13 | 0 | 0% |
| `unclassified` | 1 | 1 | 0 | 0% |

### Grouping `record`

| class | records | cells | cells spanning >=2 groups | share |
|---|---|---|---|---|
| `circuit_activity` | 314 | 156 | 49 | 31% |
| `molecular` | 170 | 132 | 30 | 23% |
| `pharmacological` | 5 | 3 | 1 | 33% |
| `structural` | 18 | 13 | 3 | 23% |
| `unclassified` | 1 | 1 | 0 | 0% |

### Grouping `study_id`

| class | records | cells | cells spanning >=2 groups | share |
|---|---|---|---|---|
| `circuit_activity` | 314 | 156 | 26 | 17% |
| `molecular` | 170 | 132 | 0 | 0% |
| `pharmacological` | 5 | 3 | 0 | 0% |
| `structural` | 18 | 13 | 0 | 0% |
| `unclassified` | 1 | 1 | 0 | 0% |

## Why: the structure of how this field publishes

The numbers above are not a sample-size problem. They are a description of publication practice, and that description is measurable whether or not the targets have been resolved yet.

| class | records | studies | labs | labs with >1 study | targeting strategies | strategies in >1 study | one study per lab |
|---|---|---|---|---|---|---|---|
| `<null>` | 1 | 1 | 1 | 0 | 0 | 0 | yes |
| `circuit_activity` | 314 | 57 | 50 | 4 | 89 | 0 | no |
| `molecular` | 170 | 16 | 16 | 0 | 28 | 0 | yes |
| `pharmacological` | 5 | 2 | 2 | 0 | 2 | 0 | yes |
| `structural` | 18 | 2 | 2 | 0 | 1 | 0 | yes |

`<null>` rests on a single study, so for it the unseen-identity fraction is not yet a statement about how the field publishes: it is a statement that the corpus holds one paper.

For `molecular`, `pharmacological`, `structural`, every study comes from a different laboratory and no laboratory contributed two: the study-to-lab map is one-to-one. Holding out a lab therefore holds out exactly one study, so `study_id` and `lab` grouping must give these classes the same unseen-identity fraction, now and after the targets are resolved. The two groupings cannot be played off against each other to make the number smaller.

## What this corpus can and cannot support

In plain terms, and without a model having been fitted:

- **`<null>` (1 records): one laboratory, so there is no cross-study question yet.** The whole class is 1 study from 1 laboratory. Held out under `lab`, there is no other lab left. This is not the one-gene-one-lab structure described below; it is a class the corpus has barely begun to cover.
- **`circuit_activity` (314 records): can support cross-study prediction, as far as identity goes.** Every manipulation here recurs across laboratories, so a held-out record's perturbation identity is one the training set has seen. That is a necessary condition for generalisation, not a sufficient one — it says a model *could* have learned something, not that it did.
- **`molecular` (170 records): cannot support cross-study prediction at all.** Every one of these records has a manipulation that appears in exactly one laboratory. Held out, there is nothing in training that shares it. No amount of modelling fixes this and no additional paper from the same groups fixes it either: the only remedy is a second laboratory repeating the same manipulation.
- **`pharmacological` (5 records): cannot be assessed from this corpus today.** Not one of these records has an identity specific enough to look for in a training set, because the field that would carry it is null. This is not a finding that they are predictable; it is a finding that the question cannot be put yet.
- **`structural` (18 records): can support cross-study prediction, as far as identity goes.** Every manipulation here recurs across laboratories, so a held-out record's perturbation identity is one the training set has seen. That is a necessary condition for generalisation, not a sufficient one — it says a model *could* have learned something, not that it did.

The reason is not a shortage of papers. It is that this field publishes one laboratory, one gene, one paper: a group establishes a gene's role in a population and the next group picks a different gene. Circuit manipulations escape this because their identity is a method everyone shares — a chemogenetic activation is the same manipulation in every laboratory that performs one — while a molecular manipulation's identity is the gene, and the gene is what distinguishes one laboratory's programme from another's.

Where a class is spread over several laboratories and its unseen-identity fraction is still 100%, the honest statement is the strong one: **that class of records can never be predicted from this corpus, under any grouping that controls leakage, until a second laboratory repeats the same manipulation.** Not "with more data", and not "with a better model". With a replication. A class that rests on one laboratory is a weaker version of the same statement: there is nothing to hold out yet.

## What this report does not claim

- It does not fit, score or cross-validate a model. It measures a property of the data that holds before any model exists.
- A low unseen-identity fraction is a necessary condition for generalisation, never a sufficient one. It says a manipulation recurred across groups, not that anything about it was learnable.
- It does not resolve targets, and it does not guess one from `targeting_strategy`. The strategy string is counted (how many are shared across studies) but never parsed into a gene: that is `target_raw`/`target_id`'s job and inferring it here would be the inference SCHEMA.md 2 rule 4 forbids.
- It cannot rule out two studies sharing one animal cohort (SCHEMA.md 5 item 3): the schema records no cohort identity.
- `pharmacological` records cannot reach `resolved` under the current `perturbation_key`, because the drug that defines them has no slot in it. They are counted as `target_unresolved` and that is a limitation of the key, which this report states rather than works around.

## The other scope, for comparison

The analysis above ran over **all records read**. The same measurement over **usable observations only (syn4rch.feasibility.coverage, schema v0.2.0 predicate)** is below, so that the choice of scope cannot be a choice of number.

| grouping | records | placed | single-group | unseen-identity fraction |
|---|---|---|---|---|
| `lab` | 492 | 492 | 170 | 34.6% |
| `population` | 492 | 492 | 122 | 24.8% |
| `population_outcome` | 492 | 492 | 0 | 0.0% |
| `record` | 492 | 492 | 0 | 0.0% |
| `study_id` | 492 | 492 | 170 | 34.6% |
