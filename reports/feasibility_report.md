# SYN4RCH feasibility report (Phase 2)

**Verdict: 1 of 93 (population, outcome) task(s) are feasible for predictive modelling (1 marginal, 91 insufficient), from 347 usable record(s) from 39 studies and 39 identified lab(s).**

Every number below is counted from the input named in the Run table. Nothing is estimated, imputed or projected. An `insufficient` verdict is a measurement of the corpus, not a failure of this module.

## Run

| field | value |
|---|---|
| generated (UTC) | 2026-11-04T04:23:38Z |
| input | `data/curated/records.jsonl` |
| input present | yes |
| records read | 362 |
| records usable as observations | 347 |
| records held out of the matrix | 15 |
| usability predicate | schema v0.2.0 |
| claim key | `population + perturbation_class + perturbation_polarity + perturbation_target + outcome_id + sex + diet + age_group + challenge` (schema v0.2.0, 9 components) |
| distinct studies | 39 |
| distinct labs (last author + affiliation) | 39 |
| usable records with no lab identity | 0 |
| usable records stating none of the 4 conditions | 146 |
| load errors | 0 |

### Why records were held out of the matrix

Counted by reason, so nothing is filtered silently. A record may appear under more than one reason. Being held out here says nothing about its admissibility to `data/curated/`, which is `syn4rch.schema.validate`'s decision.

| reason | records | meaning |
|---|---|---|
| PERTURBATION_POLARITY_UNCLEAR | 15 | perturbation_polarity is 'unclear': the manipulation's sign is not derivable even along its own class's axis, so it cannot serve as a model input; kept in the dataset, counted here, never coerced (schema v0.2.0 predicate) |

## 1. How many neuronal populations have usable intervention data?

**15.** A population counts when at least one record names it, names a perturbation `method`, and carries a `perturbation_polarity` that is neither absent nor `unclear` — a perturbation of unstated sign states nothing an input could encode. Under schema v0.1.0 this predicate read `valence` instead, which held out every molecular record; the difference is measured in the section on the key change. Of these populations, 10 are unresolved and grouped by their raw name under a `raw:` prefix (no vocabulary match; SCHEMA.md 3.2 forbids force-mapping them).

| population | records | studies | outcomes measured |
|---|---|---|---|
| `arc_agrp` | 166 | 19 | 19 |
| `arc_pomc` | 98 | 19 | 18 |
| `raw:type 1 agrp neurons` | 17 | 1 | 12 |
| `raw:arc` | 14 | 1 | 9 |
| `pvh_sim1` | 11 | 1 | 8 |
| `raw:mbh` | 11 | 1 | 8 |
| `vmh_sf1` | 9 | 1 | 3 |
| `raw:pomc neurons` | 7 | 2 | 6 |
| `raw:lepr neurons` | 5 | 1 | 3 |
| `arc_kisspeptin` | 2 | 1 | 1 |
| `raw:npy1rpvh neurons` | 2 | 1 | 2 |
| `raw:ntsth → pvn neurons` | 2 | 1 | 1 |
| `raw:drnsert → pvn projection` | 1 | 1 | 1 |
| `raw:neurons derived from pomc-expressing progenitors in the arc including pomc- and agrp-positive neurons` | 1 | 1 | 1 |
| `raw:trpv1-like receptor-expressing arc pomc neurons` | 1 | 1 | 1 |

## 2. How many physiological outcomes have been measured?

**25.** An outcome counts when at least one record resolves it to an `outcome_id` and states a `direction` other than `not_reported`. `no_change` counts: SCHEMA.md 3.5 keeps null results, and dropping them is the leading source of direction bias.

| outcome | tier | records | studies | populations measured |
|---|---|---|---|---|
| `food_intake` | A | 98 | 30 | 10 |
| `body_weight` | A | 39 | 22 | 7 |
| `energy_expenditure` | A | 27 | 12 | 7 |
| `fat_mass` | A | 22 | 12 | 7 |
| `glucose_tolerance` | A | 21 | 12 | 6 |
| `insulin_sensitivity` | A | 21 | 11 | 6 |
| `blood_glucose` | A | 17 | 7 | 3 |
| `bat_thermogenesis` | A | 16 | 6 | 6 |
| `rer` | A | 14 | 6 | 4 |
| `locomotor_activity` | A | 13 | 5 | 4 |
| `plasma_insulin` | A | 9 | 8 | 3 |
| `hepatic_lipid` | A | 8 | 5 | 4 |
| `core_temperature` | A | 7 | 4 | 5 |
| `corticosterone` | A | 7 | 4 | 3 |
| `sympathetic_nerve_activity` | A | 6 | 3 | 3 |
| `blood_pressure` | A | 5 | 2 | 2 |
| `lean_mass` | A | 4 | 4 | 4 |
| `gonadotropins` | A | 3 | 1 | 2 |
| `plasma_glucagon` | A | 3 | 2 | 1 |
| `hepatic_glucose_production` | A | 2 | 2 | 1 |
| `estrous_cyclicity` | A | 1 | 1 | 1 |
| `heart_rate` | A | 1 | 1 | 1 |
| `osmolality_avp` | A | 1 | 1 | 1 |
| `thyroid_hormone` | A | 1 | 1 | 1 |
| `wat_lipolysis` | A | 1 | 1 | 1 |

## 3. How many neuron-outcome relationships have independent replication?

The unit is the claim `syn4rch.qc.keys.claim_key`, which under schema v0.2.0 is `(population, perturbation_class, perturbation_polarity, perturbation_target, outcome_id, sex, diet, age_group, challenge)`. SCHEMA.md 1 forbids aggregating across perturbations, so a circuit activation and a gene knockout on the same pair are two claims, not two observations of one; and the condition is in the key for the same reason one step further out, because a fasted male cohort and a fed female one were never comparable.

**11 of 237 claim(s) are replicated across two or more labs**; a further 0 are replicated across studies within a single lab and are discounted (SCHEMA.md 5.2 names repeat publication by one group as the largest leakage risk here); 191 claim(s) rest on a single observation.

| field | count |
|---|---|
| claims (strict v0.2.0 claim key) | 237 |
| single observation | 191 |
| replicated across >= 2 studies | 11 |
| replicated across >= 2 labs | 11 |
| same lab only (discounted) | 0 |
| claims containing records with no lab identity | 0 |
| claims matched on 'both unknown' conditions | 80 |

### Every replication figure, by minimum condition specificity

`condition_specificity` counts how many of the 4 condition slots (`sex`, `diet`, `age_group`, `challenge`) a record actually states. The **minimum** over a cell's records is what qualifies the cell, because one record stating nothing is enough to make the whole match a match on `unknown`. A row at specificity 0 says the claim's records were grouped together partly because their conditions were not known, not because they were the same.

| min condition specificity | claims | replicated >= 2 studies | replicated >= 2 labs | agreement defined | mean agreement |
|---|---|---|---|---|---|
| 0 of 4 | 80 | 8 | 8 | 8 | 0.50 |
| 1 of 4 | 119 | 3 | 3 | 2 | 1.00 |
| 2 of 4 | 33 | 0 | 0 | 0 | undefined |
| 3 of 4 | 5 | 0 | 0 | 0 | undefined |
| 4 of 4 | 0 | 0 | 0 | 0 | undefined |

The pooled mean agreement across all specificities is 0.60. It is printed last and deliberately: a mean taken across cells that matched on nothing and cells that matched on three stated conditions describes neither, and the row-by-row breakdown above is the figure to read.

No claim in this corpus reaches specificity 4 of 4: not one record states all four conditions. Every replication figure above is therefore qualified, and none of them describes two experiments known to have run under the same conditions.

## 4. How often do independent studies agree on effect direction?

Agreement is the fraction of independent studies stating the modal direction, with a 95% Wilson score interval (chosen because the Wald interval has width zero at a unanimous 2 of 2, and the standard library has no incomplete beta for Clopper-Pearson; `syn4rch.feasibility.replication` documents this).

**Agreement is undefined, not 100%, for a claim with one independent observation**: 227 of 237 claim(s) have no second independent study, and a defined agreement exists for only 10.

Of the 10 claim(s) with a defined agreement, 2 are unanimous and 9 are contested (independent studies state incompatible directions, by the QC definition in docs/DECISIONS.md D-011). 80 of them matched on conditions that were unknown rather than known to be the same; the mean agreement is broken down by specificity in section 3 rather than pooled here. Nothing was resolved, averaged or dropped.

The `min cond` column is the minimum condition specificity of the cell, out of 4. It is not a footnote: an agreement at `0/4` was computed over studies that matched on their conditions being unknown.

| population | perturbation (class / polarity / target) | outcome | condition | min cond | studies | labs | modal direction | agreement | 95% CI | conflicts |
|---|---|---|---|---|---|---|---|---|---|---|
| `arc_agrp` | circuit_activity / increase / <null> | `food_intake` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown | 0/4 | 7 | 7 | increase | 0.71 | [0.36, 0.92] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | circuit_activity / decrease / <null> | `food_intake` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown | 0/4 | 5 | 5 | decrease | 0.80 | [0.38, 0.96] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_pomc` | circuit_activity / increase / <null> | `body_weight` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown | 0/4 | 4 | 4 | no_change | 0.50 | [0.15, 0.85] | QC_CONFLICT_EFFECT_VS_NULL, QC_CONFLICT_OPPOSED |
| `arc_pomc` | circuit_activity / increase / <null> | `blood_glucose` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown | 0/4 | 3 | 3 | increase (tied) | 0.33 | [0.06, 0.79] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_pomc` | circuit_activity / decrease / <null> | `food_intake` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown | 0/4 | 3 | 3 | increase (tied) | 0.33 | [0.06, 0.79] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_pomc` | circuit_activity / increase / <null> | `food_intake` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown | 0/4 | 3 | 3 | decrease (tied) | 0.33 | [0.06, 0.79] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | circuit_activity / decrease / <null> | `body_weight` | sex=female, diet=unknown, age_group=unknown, challenge=unknown | 1/4 | 2 | 2 | no_change | 1.00 | [0.34, 1.00] | - |
| `arc_agrp` | circuit_activity / increase / <null> | `food_intake` | sex=male, diet=unknown, age_group=unknown, challenge=unknown | 1/4 | 2 | 2 | increase | 1.00 | [0.34, 1.00] | - |
| `arc_agrp` | circuit_activity / decrease / <null> | `insulin_sensitivity` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown | 0/4 | 2 | 2 | increase (tied) | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_pomc` | circuit_activity / increase / <null> | `corticosterone` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown | 0/4 | 2 | 2 | increase (tied) | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |

### Modal direction per claim

| modal direction | claims |
|---|---|
| decrease | 63 |
| increase | 54 |
| no_change | 99 |
| undetermined | 21 |

## 5. Which outcomes have sufficient evidence for predictive modelling?

**Feasible for at least one population and valence: `food_intake`.**

Verdicts are per claim-keyed task (the nine-component schema v0.2.0 key) and then rolled up to the `(population, outcome)` pairs the Phase 2 question is phrased in; a pair is as feasible as its best claim and the rollup names which perturbation and which condition that was. No threshold in `gates.py` changed at v0.2.0 — only the unit they apply to, which is finer, so fewer cells can clear a gate that counts studies and labs. The section on the key change measures that cost.

| verdict | tasks | (population, outcome) pairs |
|---|---|---|
| feasible | 2 | 1 |
| marginal | 2 | 1 |
| insufficient | 233 | 91 |

### Tasks reaching marginal or better

| population | perturbation | outcome | condition | min cond | verdict | records | studies | labs | agreement | reasons | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `arc_agrp` | circuit_activity / increase / <null> | `food_intake` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown | 0/4 | feasible | 19 | 7 | 7 | 0.71 | `MEETS_ALL_FEASIBLE_GATES` | `CONDITIONS_ALL_UNKNOWN`, `WITHIN_STUDY_DISAGREEMENT`, `DIRECTION_DISAGREEMENT_PRESENT` |
| `arc_agrp` | circuit_activity / decrease / <null> | `food_intake` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown | 0/4 | feasible | 8 | 5 | 5 | 0.80 | `MEETS_ALL_FEASIBLE_GATES` | `CONDITIONS_ALL_UNKNOWN`, `WITHIN_STUDY_DISAGREEMENT`, `DIRECTION_DISAGREEMENT_PRESENT` |
| `arc_agrp` | circuit_activity / decrease / <null> | `body_weight` | sex=female, diet=unknown, age_group=unknown, challenge=unknown | 1/4 | marginal | 4 | 2 | 2 | 1.00 | `STUDIES_BELOW_FEASIBLE` | `CONDITIONS_PARTIALLY_UNKNOWN` |
| `arc_agrp` | circuit_activity / increase / <null> | `food_intake` | sex=male, diet=unknown, age_group=unknown, challenge=unknown | 1/4 | marginal | 4 | 2 | 2 | 1.00 | `STUDIES_BELOW_FEASIBLE` | `CONDITIONS_PARTIALLY_UNKNOWN` |

### Verdicts by minimum condition specificity

A `feasible` verdict whose cell matched on unknown conditions is still a `feasible` verdict — there is no specificity gate, because inventing a threshold here would be inventing a knob (`MIN_CONDITION_SPECIFICITY_GATE` is `None` on purpose). It is counted separately instead, so the qualification cannot be lost between this table and the headline count.

| verdict | 0/4 | 1/4 | 2/4 | 3/4 | 4/4 |
|---|---|---|---|---|---|
| feasible | 2 | 0 | 0 | 0 | 0 |
| marginal | 0 | 2 | 0 | 0 | 0 |
| insufficient | 78 | 117 | 33 | 5 | 0 |

### Why tasks did not pass

| reason | tasks | meaning |
|---|---|---|
| `AGREEMENT_UNDEFINED` | 227 | fewer than 2 studies state a direction, so agreement is undefined (not 100%) |
| `DIRECTION_CONTESTED` | 6 | independent studies do not agree on a direction (agreement below 0.60); reported as contested, not resolved |
| `MEETS_ALL_FEASIBLE_GATES` | 2 | every gate passed |
| `NO_INDEPENDENT_REPLICATION` | 226 | fewer than 2 distinct studies |
| `SINGLE_OBSERVATION` | 191 | exactly one record: nothing to replicate against |
| `STUDIES_BELOW_FEASIBLE` | 2 | fewer than 3 distinct studies |

### Promotion gates

Every threshold is a named constant in `src/syn4rch/feasibility/gates.py`, each with the reason that number was chosen. They are deliberately pessimistic: lab-level independence is required for `feasible`, agreement at one observation is undefined and undefined never passes, and the Wilson lower bound must clear the chance rate, so unanimity among very few studies lands at `marginal`.

| gate | value |
|---|---|
| MIN_STUDIES_MARGINAL | 2 |
| MIN_STUDIES_FEASIBLE | 3 |
| MIN_LABS_MARGINAL | 1 |
| MIN_LABS_FEASIBLE | 2 |
| MIN_RECORDS_FEASIBLE | 3 |
| MIN_DIRECTION_AGREEMENT_MARGINAL | 0.6 |
| MIN_DIRECTION_AGREEMENT_FEASIBLE | 0.7 |
| CHANCE_DIRECTION_RATE | 0.3333 |
| REQUIRE_AGREEMENT_CI_ABOVE_CHANCE | True |
| EXCLUDE_TIER_B_AS_TARGET | True |
| MIN_CONDITION_SPECIFICITY_GATE | None |

## What schema v0.2.0 changed, measured both ways

`claim_key` was `(population, valence, outcome_id)`. It is now `(population, perturbation_class, perturbation_polarity, perturbation_target, outcome_id, sex, diet, age_group, challenge)`. `syn4rch.qc.keys.claim_key_pooled` preserves the old key exactly, for one purpose: without it this change would be unfalsifiable. Both columns below are computed over the same usable records, so the only thing varying between them is the key.

| quantity | pooled key (v0.1.0) | strict key (v0.2.0) |
|---|---|---|
| records keyed | 347 | 347 |
| cells | 124 | 237 |
| directly opposed cells (increase and decrease in one cell) | 18 | 2 |
| cells carrying a QC conflict code across studies | 32 | 9 |
| multi-study cells | 35 | 11 |
| multi-lab cells | 35 | 11 |
| cells with a defined agreement | 34 | 10 |
| single-observation cells | 62 | 191 |

### What the split did, cell by cell

| quantity | count |
|---|---|
| pooled cells the strict key splits in two or more | 44 |
| largest number of strict cells out of one pooled cell | 10 |
| opposed pooled cells whose refinements are all internally consistent (the contradiction was heterogeneity) | 16 |
| opposed pooled cells still opposed after the split (a real head-to-head disagreement) | 2 |
| multi-study pooled cells whose every refinement has one study (replication evidence given up) | 26 |
| — of those, lost because the refinements are different perturbations | 18 |
| — of those, lost because two refinements state *different* conditions | 1 |
| — of those, lost **only** because one paper states a condition the other leaves unknown | 7 |

**What it bought:** 16 of the 18 head-to-head `increase` vs `decrease` cells under the old key were not contradictions. A finer key separated them, which means they had been a male cohort against a female one, a chow-fed cohort against a high-fat one, or a gene knockout against a gene overexpression — experiments that were never comparable, recorded as laboratories disagreeing. 2 remain, and those are the real disagreements.

**What it cost:** 26 of the 35 multi-study cells stopped being multi-study — 74.3% of the replication evidence the old key appeared to have. Most of it was never evidence: 18 of the 26 were records describing different *perturbations* (a knockout of one gene against a knockout of another), and 1 were two different stated conditions. Pooling those and counting the papers is the defect v0.2.0 exists to fix. The loss is still real in the sense that matters for modelling: there are now fewer cells that can clear a gate counting studies and labs.

**What cannot be adjudicated:** the remaining 7 lost cell(s) were split *only* because one paper states a condition the other leaves unknown. `unknown` pools "the curator did not ask" with "the paper does not say", and neither can match a stated value: a paper reporting male mice and a paper reporting no sex at all are two keys, so they stop replicating each other even if both ran males. **Nothing in this corpus can say whether they did.** The alternative — assuming the unstated condition matched — is the inference SCHEMA.md 2 rule 4 forbids, so this cost is paid and reported rather than argued away. It is counted on its own line and never as heterogeneity resolved.

### Gate arithmetic under both keys

The same gates, the same records, two keys. The pooled column is **not a result** — a pooled cell may mix two perturbations or two cohorts, so a verdict over it describes no single experiment. It is here so the cost of the correct key is visible rather than absorbed.

| verdict | tasks (pooled v0.1.0) | tasks (strict v0.2.0) | pairs (pooled) | pairs (strict) |
|---|---|---|---|---|
| feasible | 1 | 2 | 1 | 1 |
| marginal | 9 | 2 | 8 | 1 |
| insufficient | 114 | 233 | 84 | 91 |

### The usability predicate also changed

Usability used to require a `valence` that was neither absent nor `unclear`. Under v0.2.0 it requires a `perturbation_polarity` that is neither absent nor `unclear`. This is a corrected predicate, not a loosened threshold: a conditional knockout's `valence` is honestly `unclear` (SCHEMA.md 3.4 and 6 forbid coercing it) while its *manipulation* has a perfectly well defined sign, which is `decrease` on the amount of the target gene product. The old predicate tested a field that had stopped carrying the information. Every number in `gates.py` is unchanged.

| quantity | count |
|---|---|
| records read | 362 |
| usable under the v0.1.0 predicate | 177 |
| usable under the v0.2.0 predicate | 347 |
| admitted by v0.2.0 only | 170 |
| admitted by v0.1.0 only | 0 |

The records the change admits, by perturbation class — a predicate change that *admits* records is the kind that most needs showing, so it is broken down rather than stated as one number:

| perturbation class | records admitted by v0.2.0 only |
|---|---|
| `molecular` | 170 |

## Coverage matrix (population x outcome)

The matrix is built dense and the zeros are kept: which relationships nobody has measured is half the Phase 2 answer, and it disappears if only the populated cells are stored.

| field | value |
|---|---|
| matrix shape | 15 x 25 |
| cells | 375 |
| cells with at least one record | 93 |
| cells with no records | 282 |
| cells with two or more labs | 26 |
| sparsity (fraction of cells empty) | 75.2% |

The dense grid is not drawn: 15 x 25 exceeds the 12 x 12 printable limit. 282 cell(s) are empty and are listed in the JSON output by omission from `populated_cells`.

### Populated cells

| population | outcome | records | studies | labs |
|---|---|---|---|---|
| `arc_agrp` | `food_intake` | 58 | 17 | 17 |
| `arc_pomc` | `food_intake` | 27 | 14 | 14 |
| `arc_pomc` | `body_weight` | 16 | 12 | 12 |
| `arc_agrp` | `body_weight` | 15 | 9 | 9 |
| `arc_agrp` | `glucose_tolerance` | 13 | 7 | 7 |
| `arc_agrp` | `energy_expenditure` | 11 | 6 | 6 |
| `arc_agrp` | `insulin_sensitivity` | 10 | 6 | 6 |
| `arc_agrp` | `locomotor_activity` | 10 | 4 | 4 |
| `arc_pomc` | `fat_mass` | 10 | 5 | 5 |
| `arc_agrp` | `rer` | 9 | 5 | 5 |
| `arc_agrp` | `blood_glucose` | 8 | 3 | 3 |
| `arc_pomc` | `blood_glucose` | 8 | 4 | 4 |
| `arc_agrp` | `fat_mass` | 7 | 5 | 5 |
| `arc_pomc` | `insulin_sensitivity` | 7 | 6 | 6 |
| `arc_agrp` | `plasma_insulin` | 6 | 5 | 5 |
| `arc_agrp` | `bat_thermogenesis` | 5 | 3 | 3 |
| `arc_pomc` | `corticosterone` | 4 | 2 | 2 |
| `arc_pomc` | `energy_expenditure` | 4 | 3 | 3 |
| `arc_pomc` | `glucose_tolerance` | 4 | 4 | 4 |
| `arc_pomc` | `hepatic_lipid` | 4 | 2 | 2 |
| `raw:arc` | `bat_thermogenesis` | 4 | 1 | 1 |
| `raw:type 1 agrp neurons` | `energy_expenditure` | 4 | 1 | 1 |
| `vmh_sf1` | `energy_expenditure` | 4 | 1 | 1 |
| `arc_agrp` | `plasma_glucagon` | 3 | 2 | 2 |
| `arc_agrp` | `sympathetic_nerve_activity` | 3 | 2 | 2 |
| `arc_pomc` | `blood_pressure` | 3 | 1 | 1 |
| `pvh_sim1` | `food_intake` | 3 | 1 | 1 |
| `vmh_sf1` | `rer` | 3 | 1 | 1 |
| `arc_agrp` | `corticosterone` | 2 | 1 | 1 |
| `arc_agrp` | `hepatic_lipid` | 2 | 2 | 2 |
| `arc_kisspeptin` | `gonadotropins` | 2 | 1 | 1 |
| `arc_pomc` | `bat_thermogenesis` | 2 | 1 | 1 |
| `arc_pomc` | `hepatic_glucose_production` | 2 | 2 | 2 |
| `arc_pomc` | `plasma_insulin` | 2 | 2 | 2 |
| `pvh_sim1` | `body_weight` | 2 | 1 | 1 |
| `raw:arc` | `body_weight` | 2 | 1 | 1 |
| `raw:arc` | `food_intake` | 2 | 1 | 1 |
| `raw:lepr neurons` | `bat_thermogenesis` | 2 | 1 | 1 |
| `raw:lepr neurons` | `sympathetic_nerve_activity` | 2 | 1 | 1 |
| `raw:mbh` | `bat_thermogenesis` | 2 | 1 | 1 |

... 53 more row(s) not shown.

## Leakage-aware splits available on this corpus

SCHEMA.md 5: a cross-validated number reported without naming its grouping is an invalid result in this project. `syn4rch.feasibility.splits` therefore returns folds only inside a labelled fold set, and `attach_metric` raises if a caller tries to bind a number to folds without the label. Shared animal cohorts (SCHEMA.md 5 item 3) are not controlled by any grouping here: the schema records no cohort identity.

**The claim key change does not reach this table.** All three fold generators key on `study_key`, `lab_key` and `population_key`; none of them reads `claim_key`, a perturbation field or a condition field, so fold membership is identical before and after schema v0.2.0. That is checked rather than assumed: `tests/test_feasibility_splits.py` holds records that differ only in `perturbation_class`, `perturbation_polarity`, `target_id` and the four conditions, and asserts the folds come out unchanged.

Grouping is also not the whole of leakage control. SCHEMA.md 5b: a split can keep a lab off both sides and still hand the model a test row whose perturbation identity never appeared in training, which is unlearnable rather than merely hard. That fraction is measured in `reports/perturbation_identity.md` and is required next to any cross-validated number from this corpus.

| grouping | groups | usable folds | cross-validation possible | records excluded as ungroupable |
|---|---|---|---|---|
| `by_study` | 39 | 39 | yes | 0 |
| `by_lab` | 39 | 39 | yes | 0 |
| `leave_one_population_out` | 15 | 15 | yes | 0 |

## What this report does not claim

- It does not estimate how many records curation will eventually yield. Every number above is a count of what is in the input file at the time shown.
- It does not resolve, average, prefer or drop either side of a disagreement. Contested claims are reported as contested.
- It does not treat a single observation as agreement, and it does not treat repeated publication by one lab as independent replication.
- It cannot rule out two studies sharing one animal cohort (SCHEMA.md 5 item 3): the schema records no cohort identity.
- A `feasible` verdict says a prediction task can be posed and cross-validated with the grouping named, not that a model will work. It does not say the held-out rows were learnable: see `reports/perturbation_identity.md` (SCHEMA.md 5b).
- It does not claim that two records in one cell ran under the same conditions. Where the minimum condition specificity is below 4, they were grouped partly because the conditions were unknown. No claim in this corpus states all 4.
- The pooled (v0.1.0) key appears only in the section on the key change. It is never a verdict: it pools experiments that were never comparable, which is the defect v0.2.0 exists to fix.
- Tier B outcomes are excluded as prediction targets by `EXCLUDE_TIER_B_AS_TARGET` (SCHEMA.md 6: Tier B is an annotation layer, not a training target).
