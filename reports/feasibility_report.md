# SYN4RCH feasibility report (Phase 2)

**Verdict: 1 of 265 (population, outcome) task(s) are feasible for predictive modelling (3 marginal, 261 insufficient), from 678 usable record(s) from 97 studies and 85 identified lab(s).**

Every number below is counted from the input named in the Run table. Nothing is estimated, imputed or projected. An `insufficient` verdict is a measurement of the corpus, not a failure of this module.

## Run

| field | value |
|---|---|
| generated (UTC) | 2026-11-07T05:20:37Z |
| input | `data/curated/records.jsonl` |
| input present | yes |
| records read | 709 |
| records usable as observations | 678 |
| records held out of the matrix | 31 |
| usability predicate | schema v0.2.0 |
| claim key | `population + perturbation_class + perturbation_polarity + perturbation_target + outcome_id + sex + diet + age_group + challenge + background_genotype` (schema v0.2.0, 10 components) |
| distinct studies | 97 |
| distinct labs (last author, D-056) | 85 |
| usable records with no lab identity | 0 |
| usable records stating none of the 5 conditions | 262 |
| load errors | 0 |

### Why records were held out of the matrix

Counted by reason, so nothing is filtered silently. A record may appear under more than one reason. Being held out here says nothing about its admissibility to `data/curated/`, which is `syn4rch.schema.validate`'s decision.

| reason | records | meaning |
|---|---|---|
| NO_PERTURBATION_METHOD | 2 | method is absent: the intervention is unspecified |
| NO_PERTURBATION_POLARITY | 2 | perturbation_polarity is absent and cannot be derived from method: the manipulation has no sign on any axis (schema v0.2.0 predicate) |
| PERTURBATION_POLARITY_UNCLEAR | 29 | perturbation_polarity is 'unclear': the manipulation's sign is not derivable even along its own class's axis, so it cannot serve as a model input; kept in the dataset, counted here, never coerced (schema v0.2.0 predicate) |

## 1. How many neuronal populations have usable intervention data?

**77.** A population counts when at least one record names it, names a perturbation `method`, and carries a `perturbation_polarity` that is neither absent nor `unclear` — a perturbation of unstated sign states nothing an input could encode. Under schema v0.1.0 this predicate read `valence` instead, which held out every molecular record; the difference is measured in the section on the key change. Of these populations, 49 are unresolved and grouped by their raw name under a `raw:` prefix (no vocabulary match; SCHEMA.md 3.2 forbids force-mapping them).

| population | records | studies | outcomes measured |
|---|---|---|---|
| `arc_agrp` | 178 | 25 | 19 |
| `arc_pomc` | 102 | 21 | 18 |
| `lha_orexin` | 52 | 9 | 14 |
| `lha_mch` | 43 | 5 | 9 |
| `vmh_sf1` | 22 | 2 | 6 |
| `raw:npy+/lepr+ neurons` | 20 | 1 | 10 |
| `raw:type 1 agrp neurons` | 17 | 1 | 12 |
| `pvh_crh` | 15 | 6 | 10 |
| `pvh_sim1` | 15 | 3 | 8 |
| `raw:arc` | 14 | 1 | 9 |
| `raw:dmm crh neurons` | 14 | 1 | 4 |
| `raw:mbh` | 11 | 1 | 8 |
| `raw:npy−/lepr+ neurons` | 10 | 1 | 10 |
| `pvh_trh` | 8 | 1 | 6 |
| `lha_vgat` | 7 | 1 | 2 |
| `raw:dmh-located trh neurons (dmhtrh)` | 7 | 1 | 7 |
| `raw:pomc neurons` | 7 | 2 | 6 |
| `raw:trh neurons in the mpa (mpatrh)` | 7 | 1 | 7 |
| `raw:type-i nnos neurons` | 7 | 1 | 6 |
| `raw:sumvglut2 neurons` | 6 | 1 | 5 |
| `poa_galanin` | 5 | 2 | 3 |
| `pvh_oxt` | 5 | 1 | 4 |
| `raw:lepr neurons` | 5 | 1 | 3 |
| `raw:lha leprb neurons` | 5 | 1 | 3 |
| `raw:vtavgat neurons->raw:lh` | 4 | 1 | 1 |
| `raw:zbtb16arc neurons` | 4 | 1 | 4 |
| `raw:zbtb16pvh neurons` | 4 | 1 | 4 |
| `pvh_crh->lha` | 3 | 1 | 1 |
| `raw:arcghr+ neurons` | 3 | 1 | 3 |
| `raw:central 5-ht neurons` | 3 | 1 | 3 |
| `raw:dmhvglut2 neurons` | 3 | 1 | 1 |
| `raw:nac d1r neurons` | 3 | 1 | 2 |
| `raw:sumnos1` | 3 | 1 | 3 |
| `raw:sumnos1/vglut2 neurons` | 3 | 1 | 2 |
| `raw:vtagad67+` | 3 | 1 | 2 |
| `scn_vip` | 3 | 2 | 2 |
| `tmn_histamine` | 3 | 2 | 2 |
| `arc_agrp->raw:bnst` | 2 | 2 | 1 |
| `arc_agrp->raw:lha` | 2 | 2 | 1 |
| `arc_agrp->raw:pvh` | 2 | 2 | 1 |

... 37 more row(s) not shown.

## 2. How many physiological outcomes have been measured?

**37.** An outcome counts when at least one record resolves it to an `outcome_id` and states a `direction` other than `not_reported`. `no_change` counts: SCHEMA.md 3.5 keeps null results, and dropping them is the leading source of direction bias.

| outcome | tier | records | studies | populations measured |
|---|---|---|---|---|
| `food_intake` | A | 148 | 49 | 29 |
| `locomotor_activity` | A | 47 | 22 | 20 |
| `body_weight` | A | 44 | 26 | 12 |
| `energy_expenditure` | A | 41 | 17 | 15 |
| `wakefulness` | A | 36 | 15 | 21 |
| `sleep` | A | 28 | 16 | 11 |
| `nrem_sleep` | A | 26 | 15 | 14 |
| `bat_thermogenesis` | A | 25 | 12 | 13 |
| `insulin_sensitivity` | A | 25 | 13 | 7 |
| `fat_mass` | A | 24 | 13 | 9 |
| `rem_sleep` | A | 24 | 11 | 12 |
| `glucose_tolerance` | A | 23 | 13 | 8 |
| `core_temperature` | A | 22 | 11 | 15 |
| `blood_glucose` | A | 21 | 9 | 6 |
| `rer` | A | 21 | 8 | 9 |
| `sleep_fragmentation` | A | 21 | 6 | 6 |
| `corticosterone` | A | 20 | 7 | 6 |
| `plasma_insulin` | A | 14 | 11 | 7 |
| `sleep_latency` | A | 10 | 3 | 5 |
| `hepatic_lipid` | A | 9 | 6 | 5 |
| `blood_pressure` | A | 7 | 3 | 3 |
| `lean_mass` | A | 6 | 5 | 6 |
| `sympathetic_nerve_activity` | A | 6 | 3 | 3 |
| `heart_rate` | A | 5 | 4 | 4 |
| `thyroid_hormone` | A | 4 | 2 | 4 |
| `gonadotropins` | A | 3 | 1 | 2 |
| `plasma_glucagon` | A | 3 | 2 | 1 |
| `bone_density` | A | 2 | 1 | 2 |
| `estrous_cyclicity` | A | 2 | 2 | 2 |
| `gut_motility` | A | 2 | 1 | 1 |
| `hepatic_glucose_production` | A | 2 | 2 | 1 |
| `inflammatory_markers_other` | A | 2 | 1 | 1 |
| `il6` | A | 1 | 1 | 1 |
| `osmolality_avp` | A | 1 | 1 | 1 |
| `tnf` | A | 1 | 1 | 1 |
| `wat_lipolysis` | A | 1 | 1 | 1 |
| `water_intake` | A | 1 | 1 | 1 |

## 3. How many neuron-outcome relationships have independent replication?

The unit is the claim `syn4rch.qc.keys.claim_key`, which under schema v0.2.0 is `(population, perturbation_class, perturbation_polarity, perturbation_target, outcome_id, sex, diet, age_group, challenge, background_genotype)`. SCHEMA.md 1 forbids aggregating across perturbations, so a circuit activation and a gene knockout on the same pair are two claims, not two observations of one; and the condition is in the key for the same reason one step further out, because a fasted male cohort and a fed female one were never comparable.

**12 of 517 claim(s) are replicated across two or more labs**; a further 0 are replicated across studies within a single lab and are discounted (SCHEMA.md 5.2 names repeat publication by one group as the largest leakage risk here); 420 claim(s) rest on a single observation.

| field | count |
|---|---|
| claims (strict v0.2.0 claim key) | 517 |
| single observation | 420 |
| replicated across >= 2 studies | 12 |
| replicated across >= 2 labs | 12 |
| same lab only (discounted) | 0 |
| claims containing records with no lab identity | 0 |
| claims matched on 'both unknown' conditions | 192 |

### Every replication figure, by minimum condition specificity

`condition_specificity` counts how many of the 5 condition slots (`sex`, `diet`, `age_group`, `challenge`, `background_genotype`) a record actually states. The **minimum** over a cell's records is what qualifies the cell, because one record stating nothing is enough to make the whole match a match on `unknown`. A row at specificity 0 says the claim's records were grouped together partly because their conditions were not known, not because they were the same.

| min condition specificity | claims | replicated >= 2 studies | replicated >= 2 labs | agreement defined | mean agreement |
|---|---|---|---|---|---|
| 0 of 5 | 192 | 9 | 9 | 9 | 0.69 |
| 1 of 5 | 137 | 2 | 2 | 2 | 1.00 |
| 2 of 5 | 172 | 1 | 1 | 1 | 0.50 |
| 3 of 5 | 16 | 0 | 0 | 0 | undefined |
| 4 of 5 | 0 | 0 | 0 | 0 | undefined |
| 5 of 5 | 0 | 0 | 0 | 0 | undefined |

The pooled mean agreement across all specificities is 0.72. It is printed last and deliberately: a mean taken across cells that matched on nothing and cells that matched on three stated conditions describes neither, and the row-by-row breakdown above is the figure to read.

No claim in this corpus reaches specificity 5 of 5: not one record states every condition. Every replication figure above is therefore qualified, and none of them describes two experiments known to have run under the same conditions.

## 4. How often do independent studies agree on effect direction?

Agreement is the fraction of independent studies stating the modal direction, with a 95% Wilson score interval (chosen because the Wald interval has width zero at a unanimous 2 of 2, and the standard library has no incomplete beta for Clopper-Pearson; `syn4rch.feasibility.replication` documents this).

**Agreement is undefined, not 100%, for a claim with one independent observation**: 505 of 517 claim(s) have no second independent study, and a defined agreement exists for only 12.

Of the 12 claim(s) with a defined agreement, 5 are unanimous and 7 are contested (independent studies state incompatible directions, by the QC definition in docs/DECISIONS.md D-011). 192 of them matched on conditions that were unknown rather than known to be the same; the mean agreement is broken down by specificity in section 3 rather than pooled here. Nothing was resolved, averaged or dropped.

The `min cond` column is the minimum condition specificity of the cell, out of 5. It is not a footnote: an agreement at `0/5` was computed over studies that matched on their conditions being unknown.

| population | perturbation (class / polarity / target) | outcome | condition | min cond | studies | labs | modal direction | agreement | 95% CI | conflicts |
|---|---|---|---|---|---|---|---|---|---|---|
| `lha_orexin` | circuit_activity / increase / <null> | `sleep` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | 4 | 4 | no_change | 0.50 | [0.15, 0.85] | QC_CONFLICT_EFFECT_VS_NULL, QC_CONFLICT_OPPOSED |
| `arc_agrp` | circuit_activity / increase / <null> | `food_intake` | sex=unknown, diet=unknown, age_group=unknown, challenge=ad_libitum, background_genotype=unknown | 1/5 | 3 | 3 | increase | 1.00 | [0.44, 1.00] | - |
| `lha_orexin` | circuit_activity / increase / <null> | `wakefulness` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | 3 | 3 | no_change | 0.67 | [0.21, 0.94] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | circuit_activity / increase / <null> | `food_intake` | sex=male, diet=chow, age_group=unknown, challenge=unknown, background_genotype=unknown | 2/5 | 2 | 2 | increase | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | circuit_activity / increase / <null> | `food_intake` | sex=unknown, diet=chow, age_group=unknown, challenge=unknown, background_genotype=unknown | 1/5 | 2 | 2 | increase | 1.00 | [0.34, 1.00] | - |
| `arc_agrp` | circuit_activity / increase / <null> | `food_intake` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | 2 | 2 | increase | 1.00 | [0.34, 1.00] | - |
| `lha_mch` | circuit_activity / increase / <null> | `nrem_sleep` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | 2 | 2 | increase (tied) | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |
| `lha_mch` | circuit_activity / increase / <null> | `rem_sleep` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | 2 | 2 | increase | 1.00 | [0.34, 1.00] | - |
| `lha_orexin` | circuit_activity / increase / <null> | `locomotor_activity` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | 2 | 2 | increase (tied) | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |
| `lha_orexin` | circuit_activity / decrease / <null> | `nrem_sleep` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | 2 | 2 | decrease (tied) | 0.50 | [0.09, 0.91] | QC_CONFLICT_OPPOSED |
| `lha_orexin` | circuit_activity / increase / <null> | `nrem_sleep` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | 2 | 2 | decrease | 1.00 | [0.34, 1.00] | - |
| `scn_vip` | circuit_activity / increase / <null> | `sleep` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | 2 | 2 | increase (tied) | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |

### Modal direction per claim

| modal direction | claims |
|---|---|
| decrease | 115 |
| increase | 172 |
| no_change | 184 |
| undetermined | 46 |

## 5. Which outcomes have sufficient evidence for predictive modelling?

**Feasible for at least one population and valence: `food_intake`.**

Verdicts are per claim-keyed task (the full claim key) and then rolled up to the `(population, outcome)` pairs the Phase 2 question is phrased in; a pair is as feasible as its best claim and the rollup names which perturbation and which condition that was. No threshold in `gates.py` changed at v0.2.0 — only the unit they apply to, which is finer, so fewer cells can clear a gate that counts studies and labs. The section on the key change measures that cost.

| verdict | tasks | (population, outcome) pairs |
|---|---|---|
| feasible | 1 | 1 |
| marginal | 5 | 3 |
| insufficient | 511 | 261 |

### Tasks reaching marginal or better

| population | perturbation | outcome | condition | min cond | verdict | records | studies | labs | agreement | reasons | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `arc_agrp` | circuit_activity / increase / <null> | `food_intake` | sex=unknown, diet=unknown, age_group=unknown, challenge=ad_libitum, background_genotype=unknown | 1/5 | feasible | 3 | 3 | 3 | 1.00 | `MEETS_ALL_FEASIBLE_GATES` | `CONDITIONS_PARTIALLY_UNKNOWN` |
| `lha_orexin` | circuit_activity / increase / <null> | `wakefulness` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | marginal | 3 | 3 | 3 | 0.67 | `AGREEMENT_BELOW_FEASIBLE`, `AGREEMENT_CI_INCLUDES_CHANCE` | `CONDITIONS_ALL_UNKNOWN`, `DIRECTION_DISAGREEMENT_PRESENT` |
| `arc_agrp` | circuit_activity / increase / <null> | `food_intake` | sex=unknown, diet=chow, age_group=unknown, challenge=unknown, background_genotype=unknown | 1/5 | marginal | 4 | 2 | 2 | 1.00 | `STUDIES_BELOW_FEASIBLE` | `CONDITIONS_PARTIALLY_UNKNOWN` |
| `arc_agrp` | circuit_activity / increase / <null> | `food_intake` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | marginal | 2 | 2 | 2 | 1.00 | `STUDIES_BELOW_FEASIBLE`, `RECORDS_BELOW_FEASIBLE` | `CONDITIONS_ALL_UNKNOWN` |
| `lha_mch` | circuit_activity / increase / <null> | `rem_sleep` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | marginal | 3 | 2 | 2 | 1.00 | `STUDIES_BELOW_FEASIBLE` | `CONDITIONS_ALL_UNKNOWN` |
| `lha_orexin` | circuit_activity / increase / <null> | `nrem_sleep` | sex=unknown, diet=unknown, age_group=unknown, challenge=unknown, background_genotype=unknown | 0/5 | marginal | 2 | 2 | 2 | 1.00 | `STUDIES_BELOW_FEASIBLE`, `RECORDS_BELOW_FEASIBLE` | `CONDITIONS_ALL_UNKNOWN` |

### Verdicts by minimum condition specificity

A `feasible` verdict whose cell matched on unknown conditions is still a `feasible` verdict — there is no specificity gate, because inventing a threshold here would be inventing a knob (`MIN_CONDITION_SPECIFICITY_GATE` is `None` on purpose). It is counted separately instead, so the qualification cannot be lost between this table and the headline count.

| verdict | 0/5 | 1/5 | 2/5 | 3/5 | 4/5 | 5/5 |
|---|---|---|---|---|---|---|
| feasible | 0 | 1 | 0 | 0 | 0 | 0 |
| marginal | 4 | 1 | 0 | 0 | 0 | 0 |
| insufficient | 188 | 135 | 172 | 16 | 0 | 0 |

### Why tasks did not pass

| reason | tasks | meaning |
|---|---|---|
| `AGREEMENT_BELOW_FEASIBLE` | 1 | direction agreement below 0.70 |
| `AGREEMENT_CI_INCLUDES_CHANCE` | 1 | the Wilson 95% lower bound on agreement does not clear the 0.33 chance rate |
| `AGREEMENT_UNDEFINED` | 505 | fewer than 2 studies state a direction, so agreement is undefined (not 100%) |
| `DIRECTION_CONTESTED` | 6 | independent studies do not agree on a direction (agreement below 0.60); reported as contested, not resolved |
| `MEETS_ALL_FEASIBLE_GATES` | 1 | every gate passed |
| `NO_INDEPENDENT_REPLICATION` | 505 | fewer than 2 distinct studies |
| `RECORDS_BELOW_FEASIBLE` | 2 | fewer than 3 records |
| `SINGLE_OBSERVATION` | 420 | exactly one record: nothing to replicate against |
| `STUDIES_BELOW_FEASIBLE` | 4 | fewer than 3 distinct studies |

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
| records keyed | 678 | 678 |
| cells | 328 | 517 |
| directly opposed cells (increase and decrease in one cell) | 29 | 12 |
| cells carrying a QC conflict code across studies | 43 | 7 |
| multi-study cells | 53 | 12 |
| multi-lab cells | 53 | 12 |
| cells with a defined agreement | 51 | 12 |
| single-observation cells | 211 | 420 |

### What the split did, cell by cell

| quantity | count |
|---|---|
| pooled cells the strict key splits in two or more | 74 |
| largest number of strict cells out of one pooled cell | 14 |
| opposed pooled cells whose refinements are all internally consistent (the contradiction was heterogeneity) | 18 |
| opposed pooled cells still opposed after the split (a real head-to-head disagreement) | 11 |
| multi-study pooled cells whose every refinement has one study (replication evidence given up) | 44 |
| — of those, lost because the refinements are different perturbations | 18 |
| — of those, lost because two refinements state *different* conditions | 9 |
| — of those, lost **only** because one paper states a condition the other leaves unknown | 17 |

**What it bought:** 18 of the 29 head-to-head `increase` vs `decrease` cells under the old key were not contradictions. A finer key separated them, which means they had been a male cohort against a female one, a chow-fed cohort against a high-fat one, or a gene knockout against a gene overexpression — experiments that were never comparable, recorded as laboratories disagreeing. 11 remain, and those are the real disagreements.

**What it cost:** 44 of the 53 multi-study cells stopped being multi-study — 83.0% of the replication evidence the old key appeared to have. Most of it was never evidence: 18 of the 44 were records describing different *perturbations* (a knockout of one gene against a knockout of another), and 9 were two different stated conditions. Pooling those and counting the papers is the defect v0.2.0 exists to fix. The loss is still real in the sense that matters for modelling: there are now fewer cells that can clear a gate counting studies and labs.

**What cannot be adjudicated:** the remaining 17 lost cell(s) were split *only* because one paper states a condition the other leaves unknown. `unknown` pools "the curator did not ask" with "the paper does not say", and neither can match a stated value: a paper reporting male mice and a paper reporting no sex at all are two keys, so they stop replicating each other even if both ran males. **Nothing in this corpus can say whether they did.** The alternative — assuming the unstated condition matched — is the inference SCHEMA.md 2 rule 4 forbids, so this cost is paid and reported rather than argued away. It is counted on its own line and never as heterogeneity resolved.

### Gate arithmetic under both keys

The same gates, the same records, two keys. The pooled column is **not a result** — a pooled cell may mix two perturbations or two cohorts, so a verdict over it describes no single experiment. It is here so the cost of the correct key is visible rather than absorbed.

| verdict | tasks (pooled v0.1.0) | tasks (strict v0.2.0) | pairs (pooled) | pairs (strict) |
|---|---|---|---|---|
| feasible | 2 | 1 | 2 | 1 |
| marginal | 16 | 5 | 15 | 3 |
| insufficient | 310 | 511 | 248 | 261 |

### The usability predicate also changed

Usability used to require a `valence` that was neither absent nor `unclear`. Under v0.2.0 it requires a `perturbation_polarity` that is neither absent nor `unclear`. This is a corrected predicate, not a loosened threshold: a conditional knockout's `valence` is honestly `unclear` (SCHEMA.md 3.4 and 6 forbid coercing it) while its *manipulation* has a perfectly well defined sign, which is `decrease` on the amount of the target gene product. The old predicate tested a field that had stopped carrying the information. Every number in `gates.py` is unchanged.

| quantity | count |
|---|---|
| records read | 709 |
| usable under the v0.1.0 predicate | 506 |
| usable under the v0.2.0 predicate | 678 |
| admitted by v0.2.0 only | 172 |
| admitted by v0.1.0 only | 0 |

The records the change admits, by perturbation class — a predicate change that *admits* records is the kind that most needs showing, so it is broken down rather than stated as one number:

| perturbation class | records admitted by v0.2.0 only |
|---|---|
| `molecular` | 172 |

## Coverage matrix (population x outcome)

The matrix is built dense and the zeros are kept: which relationships nobody has measured is half the Phase 2 answer, and it disappears if only the populated cells are stored.

| field | value |
|---|---|
| matrix shape | 77 x 37 |
| cells | 2849 |
| cells with at least one record | 265 |
| cells with no records | 2584 |
| cells with two or more labs | 48 |
| sparsity (fraction of cells empty) | 90.7% |

The dense grid is not drawn: 77 x 37 exceeds the 12 x 12 printable limit. 2584 cell(s) are empty and are listed in the JSON output by omission from `populated_cells`.

### Populated cells

| population | outcome | records | studies | labs |
|---|---|---|---|---|
| `arc_agrp` | `food_intake` | 62 | 23 | 20 |
| `arc_pomc` | `food_intake` | 29 | 15 | 15 |
| `arc_pomc` | `body_weight` | 16 | 12 | 12 |
| `arc_agrp` | `body_weight` | 15 | 9 | 8 |
| `arc_agrp` | `locomotor_activity` | 14 | 5 | 5 |
| `arc_agrp` | `energy_expenditure` | 13 | 7 | 6 |
| `arc_agrp` | `glucose_tolerance` | 13 | 7 | 7 |
| `lha_orexin` | `sleep` | 13 | 8 | 6 |
| `arc_agrp` | `insulin_sensitivity` | 12 | 7 | 6 |
| `arc_pomc` | `fat_mass` | 10 | 5 | 5 |
| `arc_agrp` | `rer` | 9 | 5 | 5 |
| `lha_mch` | `rem_sleep` | 9 | 2 | 2 |
| `arc_agrp` | `blood_glucose` | 8 | 3 | 3 |
| `arc_pomc` | `blood_glucose` | 8 | 4 | 4 |
| `raw:dmm crh neurons` | `sleep_fragmentation` | 8 | 1 | 1 |
| `vmh_sf1` | `food_intake` | 8 | 2 | 2 |
| `arc_agrp` | `fat_mass` | 7 | 5 | 5 |
| `arc_pomc` | `insulin_sensitivity` | 7 | 6 | 6 |
| `lha_orexin` | `corticosterone` | 7 | 1 | 1 |
| `pvh_sim1` | `food_intake` | 7 | 3 | 3 |
| `arc_agrp` | `plasma_insulin` | 6 | 5 | 5 |
| `lha_mch` | `food_intake` | 6 | 2 | 2 |
| `lha_mch` | `nrem_sleep` | 6 | 3 | 3 |
| `lha_mch` | `sleep_fragmentation` | 6 | 2 | 2 |
| `lha_mch` | `sleep_latency` | 6 | 1 | 1 |
| `lha_mch` | `wakefulness` | 6 | 2 | 2 |
| `lha_orexin` | `wakefulness` | 6 | 4 | 3 |
| `arc_agrp` | `bat_thermogenesis` | 5 | 3 | 3 |
| `lha_orexin` | `locomotor_activity` | 5 | 2 | 2 |
| `lha_orexin` | `nrem_sleep` | 5 | 3 | 2 |
| `lha_vgat` | `food_intake` | 5 | 1 | 1 |
| `arc_pomc` | `corticosterone` | 4 | 2 | 2 |
| `arc_pomc` | `energy_expenditure` | 4 | 3 | 3 |
| `arc_pomc` | `glucose_tolerance` | 4 | 4 | 4 |
| `arc_pomc` | `hepatic_lipid` | 4 | 2 | 2 |
| `raw:arc` | `bat_thermogenesis` | 4 | 1 | 1 |
| `raw:type 1 agrp neurons` | `energy_expenditure` | 4 | 1 | 1 |
| `raw:vtavgat neurons->raw:lh` | `locomotor_activity` | 4 | 1 | 1 |
| `vmh_sf1` | `energy_expenditure` | 4 | 1 | 1 |
| `arc_agrp` | `plasma_glucagon` | 3 | 2 | 2 |

... 225 more row(s) not shown.

## Leakage-aware splits available on this corpus

SCHEMA.md 5: a cross-validated number reported without naming its grouping is an invalid result in this project. `syn4rch.feasibility.splits` therefore returns folds only inside a labelled fold set, and `attach_metric` raises if a caller tries to bind a number to folds without the label. Shared animal cohorts (SCHEMA.md 5 item 3) are not controlled by any grouping here: the schema records no cohort identity.

**The claim key change does not reach this table.** All three fold generators key on `study_key`, `lab_key` and `population_key`; none of them reads `claim_key`, a perturbation field or a condition field, so fold membership is identical before and after schema v0.2.0. That is checked rather than assumed: `tests/test_feasibility_splits.py` holds records that differ only in `perturbation_class`, `perturbation_polarity`, `target_id` and the every condition, and asserts the folds come out unchanged.

Grouping is also not the whole of leakage control. SCHEMA.md 5b: a split can keep a lab off both sides and still hand the model a test row whose perturbation identity never appeared in training, which is unlearnable rather than merely hard. That fraction is measured in `reports/perturbation_identity.md` and is required next to any cross-validated number from this corpus.

| grouping | groups | usable folds | cross-validation possible | records excluded as ungroupable |
|---|---|---|---|---|
| `by_study` | 97 | 97 | yes | 0 |
| `by_lab` | 85 | 85 | yes | 0 |
| `leave_one_population_out` | 77 | 77 | yes | 0 |

## What this report does not claim

- It does not estimate how many records curation will eventually yield. Every number above is a count of what is in the input file at the time shown.
- It does not resolve, average, prefer or drop either side of a disagreement. Contested claims are reported as contested.
- It does not treat a single observation as agreement, and it does not treat repeated publication by one lab as independent replication.
- It cannot rule out two studies sharing one animal cohort (SCHEMA.md 5 item 3): the schema records no cohort identity.
- A `feasible` verdict says a prediction task can be posed and cross-validated with the grouping named, not that a model will work. It does not say the held-out rows were learnable: see `reports/perturbation_identity.md` (SCHEMA.md 5b).
- It does not claim that two records in one cell ran under the same conditions. Where the minimum condition specificity is below 5, they were grouped partly because the conditions were unknown. No claim in this corpus states all 5.
- The pooled (v0.1.0) key appears only in the section on the key change. It is never a verdict: it pools experiments that were never comparable, which is the defect v0.2.0 exists to fix.
- Tier B outcomes are excluded as prediction targets by `EXCLUDE_TIER_B_AS_TARGET` (SCHEMA.md 6: Tier B is an annotation layer, not a training target).
