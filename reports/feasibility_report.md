# SYN4RCH feasibility report (Phase 2)

**Verdict: 1 of 46 (population, outcome) task(s) are feasible for predictive modelling (2 marginal, 43 insufficient), from 177 usable record(s) from 25 studies and 25 identified lab(s).**

Every number below is counted from the input named in the Run table. Nothing is estimated, imputed or projected. An `insufficient` verdict is a measurement of the corpus, not a failure of this module.

## Run

| field | value |
|---|---|
| generated (UTC) | 2026-11-04T00:53:38Z |
| input | `data/curated/records.jsonl` |
| input present | yes |
| records read | 362 |
| records usable as observations | 177 |
| records held out of the matrix | 185 |
| distinct studies | 25 |
| distinct labs (last author + affiliation) | 25 |
| usable records with no lab identity | 0 |
| load errors | 0 |

### Why records were held out of the matrix

Counted by reason, so nothing is filtered silently. A record may appear under more than one reason. Being held out here says nothing about its admissibility to `data/curated/`, which is `syn4rch.schema.validate`'s decision.

| reason | records | meaning |
|---|---|---|
| VALENCE_UNCLEAR | 185 | valence is 'unclear' -- a valid schema answer (SCHEMA.md 3.4) that states no direction of perturbation, so it cannot serve as a model input; kept in the dataset, counted here, never coerced |

## 1. How many neuronal populations have usable intervention data?

**10.** A population counts when at least one record names it, names a perturbation `method`, and carries a `valence` that is neither absent nor `unclear` — a perturbation of unstated direction states nothing an input could encode. Of these, 6 are unresolved populations grouped by their raw name under a `raw:` prefix (no vocabulary match; SCHEMA.md 3.2 forbids force-mapping them).

| population | records | studies | outcomes measured |
|---|---|---|---|
| `arc_agrp` | 101 | 13 | 14 |
| `arc_pomc` | 41 | 11 | 10 |
| `raw:type 1 agrp neurons` | 17 | 1 | 12 |
| `vmh_sf1` | 9 | 1 | 3 |
| `arc_kisspeptin` | 2 | 1 | 1 |
| `raw:npy1rpvh neurons` | 2 | 1 | 2 |
| `raw:ntsth → pvn neurons` | 2 | 1 | 1 |
| `raw:drnsert → pvn projection` | 1 | 1 | 1 |
| `raw:neurons derived from pomc-expressing progenitors in the arc including pomc- and agrp-positive neurons` | 1 | 1 | 1 |
| `raw:trpv1-like receptor-expressing arc pomc neurons` | 1 | 1 | 1 |

## 2. How many physiological outcomes have been measured?

**21.** An outcome counts when at least one record resolves it to an `outcome_id` and states a `direction` other than `not_reported`. `no_change` counts: SCHEMA.md 3.5 keeps null results, and dropping them is the leading source of direction bias.

| outcome | tier | records | studies | populations measured |
|---|---|---|---|---|
| `food_intake` | A | 74 | 21 | 8 |
| `blood_glucose` | A | 14 | 5 | 3 |
| `energy_expenditure` | A | 13 | 4 | 3 |
| `body_weight` | A | 12 | 7 | 3 |
| `insulin_sensitivity` | A | 11 | 4 | 3 |
| `rer` | A | 9 | 3 | 3 |
| `corticosterone` | A | 7 | 4 | 3 |
| `glucose_tolerance` | A | 7 | 3 | 3 |
| `plasma_insulin` | A | 6 | 5 | 3 |
| `locomotor_activity` | A | 4 | 1 | 1 |
| `sympathetic_nerve_activity` | A | 4 | 2 | 2 |
| `plasma_glucagon` | A | 3 | 2 | 1 |
| `bat_thermogenesis` | A | 2 | 2 | 2 |
| `blood_pressure` | A | 2 | 1 | 1 |
| `core_temperature` | A | 2 | 1 | 1 |
| `gonadotropins` | A | 2 | 1 | 1 |
| `fat_mass` | A | 1 | 1 | 1 |
| `heart_rate` | A | 1 | 1 | 1 |
| `hepatic_glucose_production` | A | 1 | 1 | 1 |
| `hepatic_lipid` | A | 1 | 1 | 1 |
| `osmolality_avp` | A | 1 | 1 | 1 |

## 3. How many neuron-outcome relationships have independent replication?

The unit is the claim `(population, valence, outcome)`: SCHEMA.md 1 forbids aggregating across valences, so an activation result and an inhibition result on the same pair are two claims, not two observations of one.

**17 of 59 claim(s) are replicated across two or more labs**; a further 0 are replicated across studies within a single lab and are discounted (SCHEMA.md 5.2 names repeat publication by one group as the largest leakage risk here); 28 claim(s) rest on a single observation.

| field | count |
|---|---|
| claims (population x valence x outcome) | 59 |
| single observation | 28 |
| replicated across >= 2 studies | 17 |
| replicated across >= 2 labs | 17 |
| same lab only (discounted) | 0 |
| claims containing records with no lab identity | 0 |

## 4. How often do independent studies agree on effect direction?

Agreement is the fraction of independent studies stating the modal direction, with a 95% Wilson score interval (chosen because the Wald interval has width zero at a unanimous 2 of 2, and the standard library has no incomplete beta for Clopper-Pearson; `syn4rch.feasibility.replication` documents this).

**Agreement is undefined, not 100%, for a claim with one independent observation**: 42 of 59 claim(s) have no second independent study, and a defined agreement exists for only 17.

Of the 17 claim(s) with a defined agreement, 1 are unanimous and 16 are contested (independent studies state incompatible directions, by the QC definition in docs/DECISIONS.md D-011). Mean agreement where defined: 0.51. Nothing was resolved, averaged or dropped.

| population | valence | outcome | studies | labs | modal direction | agreement | 95% CI | conflicts |
|---|---|---|---|---|---|---|---|---|
| `arc_agrp` | activation | `food_intake` | 9 | 9 | increase | 0.67 | [0.35, 0.88] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | inhibition | `food_intake` | 7 | 7 | decrease | 0.57 | [0.25, 0.84] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_pomc` | activation | `food_intake` | 5 | 5 | decrease (tied) | 0.20 | [0.04, 0.62] | QC_CONFLICT_EFFECT_VS_NULL, QC_CONFLICT_OPPOSED |
| `arc_pomc` | activation | `body_weight` | 4 | 4 | no_change | 0.50 | [0.15, 0.85] | QC_CONFLICT_EFFECT_VS_NULL, QC_CONFLICT_OPPOSED |
| `arc_pomc` | inhibition | `food_intake` | 4 | 4 | increase (tied) | 0.25 | [0.05, 0.70] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | inhibition | `body_weight` | 3 | 3 | no_change | 1.00 | [0.44, 1.00] | - |
| `arc_agrp` | activation | `insulin_sensitivity` | 3 | 3 | decrease | 0.67 | [0.21, 0.94] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_pomc` | activation | `blood_glucose` | 3 | 3 | increase (tied) | 0.33 | [0.06, 0.79] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | activation | `blood_glucose` | 2 | 2 | no_change | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | activation | `energy_expenditure` | 2 | 2 | no_change | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | activation | `glucose_tolerance` | 2 | 2 | no_change | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | inhibition | `insulin_sensitivity` | 2 | 2 | increase (tied) | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | activation | `plasma_glucagon` | 2 | 2 | increase (tied) | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | activation | `plasma_insulin` | 2 | 2 | increase (tied) | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | activation | `rer` | 2 | 2 | increase | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_agrp` | activation | `sympathetic_nerve_activity` | 2 | 2 | decrease | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |
| `arc_pomc` | activation | `corticosterone` | 2 | 2 | increase (tied) | 0.50 | [0.09, 0.91] | QC_CONFLICT_EFFECT_VS_NULL |

### Modal direction per claim

| modal direction | claims |
|---|---|
| decrease | 15 |
| increase | 15 |
| no_change | 22 |
| undetermined | 7 |

## 5. Which outcomes have sufficient evidence for predictive modelling?

**Feasible for at least one population and valence: `body_weight`.**

Verdicts are per task `(population, valence, outcome)` and then rolled up to the `(population, outcome)` pairs the Phase 2 question is phrased in; a pair is as feasible as its best valence and the rollup names which one.

| verdict | tasks | (population, outcome) pairs |
|---|---|---|
| feasible | 1 | 1 |
| marginal | 2 | 2 |
| insufficient | 56 | 43 |

### Tasks reaching marginal or better

| population | valence | outcome | verdict | records | studies | labs | agreement | reasons |
|---|---|---|---|---|---|---|---|---|
| `arc_agrp` | inhibition | `body_weight` | feasible | 5 | 3 | 3 | 1.00 | `MEETS_ALL_FEASIBLE_GATES` |
| `arc_agrp` | activation | `food_intake` | marginal | 24 | 9 | 9 | 0.67 | `AGREEMENT_BELOW_FEASIBLE` |
| `arc_agrp` | activation | `insulin_sensitivity` | marginal | 6 | 3 | 3 | 0.67 | `AGREEMENT_BELOW_FEASIBLE`, `AGREEMENT_CI_INCLUDES_CHANCE` |

### Why tasks did not pass

| reason | tasks | meaning |
|---|---|---|
| `AGREEMENT_BELOW_FEASIBLE` | 2 | direction agreement below 0.70 |
| `AGREEMENT_CI_INCLUDES_CHANCE` | 1 | the Wilson 95% lower bound on agreement does not clear the 0.33 chance rate |
| `AGREEMENT_UNDEFINED` | 42 | fewer than 2 studies state a direction, so agreement is undefined (not 100%) |
| `DIRECTION_CONTESTED` | 14 | independent studies do not agree on a direction (agreement below 0.60); reported as contested, not resolved |
| `MEETS_ALL_FEASIBLE_GATES` | 1 | every gate passed |
| `NO_INDEPENDENT_REPLICATION` | 42 | fewer than 2 distinct studies |
| `SINGLE_OBSERVATION` | 28 | exactly one record: nothing to replicate against |

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

## Coverage matrix (population x outcome)

The matrix is built dense and the zeros are kept: which relationships nobody has measured is half the Phase 2 answer, and it disappears if only the populated cells are stored.

| field | value |
|---|---|
| matrix shape | 10 x 21 |
| cells | 210 |
| cells with at least one record | 46 |
| cells with no records | 164 |
| cells with two or more labs | 15 |
| sparsity (fraction of cells empty) | 78.1% |

The dense grid is not drawn: 10 x 21 exceeds the 12 x 12 printable limit. 164 cell(s) are empty and are listed in the JSON output by omission from `populated_cells`.

### Populated cells

| population | outcome | records | studies | labs |
|---|---|---|---|---|
| `arc_agrp` | `food_intake` | 47 | 13 | 13 |
| `arc_pomc` | `food_intake` | 19 | 8 | 8 |
| `arc_agrp` | `insulin_sensitivity` | 8 | 4 | 4 |
| `arc_agrp` | `blood_glucose` | 7 | 2 | 2 |
| `arc_agrp` | `body_weight` | 6 | 3 | 3 |
| `arc_pomc` | `blood_glucose` | 6 | 3 | 3 |
| `arc_agrp` | `energy_expenditure` | 5 | 2 | 2 |
| `arc_agrp` | `glucose_tolerance` | 5 | 2 | 2 |
| `arc_agrp` | `rer` | 5 | 2 | 2 |
| `arc_pomc` | `body_weight` | 5 | 4 | 4 |
| `arc_agrp` | `locomotor_activity` | 4 | 1 | 1 |
| `arc_agrp` | `plasma_insulin` | 4 | 3 | 3 |
| `arc_pomc` | `corticosterone` | 4 | 2 | 2 |
| `raw:type 1 agrp neurons` | `energy_expenditure` | 4 | 1 | 1 |
| `vmh_sf1` | `energy_expenditure` | 4 | 1 | 1 |
| `arc_agrp` | `plasma_glucagon` | 3 | 2 | 2 |
| `arc_agrp` | `sympathetic_nerve_activity` | 3 | 2 | 2 |
| `vmh_sf1` | `rer` | 3 | 1 | 1 |
| `arc_agrp` | `corticosterone` | 2 | 1 | 1 |
| `arc_kisspeptin` | `gonadotropins` | 2 | 1 | 1 |
| `arc_pomc` | `insulin_sensitivity` | 2 | 2 | 2 |
| `raw:ntsth → pvn neurons` | `food_intake` | 2 | 1 | 1 |
| `raw:type 1 agrp neurons` | `blood_pressure` | 2 | 1 | 1 |
| `raw:type 1 agrp neurons` | `core_temperature` | 2 | 1 | 1 |
| `vmh_sf1` | `food_intake` | 2 | 1 | 1 |
| `arc_agrp` | `bat_thermogenesis` | 1 | 1 | 1 |
| `arc_agrp` | `hepatic_lipid` | 1 | 1 | 1 |
| `arc_pomc` | `glucose_tolerance` | 1 | 1 | 1 |
| `arc_pomc` | `hepatic_glucose_production` | 1 | 1 | 1 |
| `arc_pomc` | `plasma_insulin` | 1 | 1 | 1 |
| `arc_pomc` | `rer` | 1 | 1 | 1 |
| `arc_pomc` | `sympathetic_nerve_activity` | 1 | 1 | 1 |
| `raw:drnsert → pvn projection` | `food_intake` | 1 | 1 | 1 |
| `raw:neurons derived from pomc-expressing progenitors in the arc including pomc- and agrp-positive neurons` | `food_intake` | 1 | 1 | 1 |
| `raw:npy1rpvh neurons` | `glucose_tolerance` | 1 | 1 | 1 |
| `raw:npy1rpvh neurons` | `insulin_sensitivity` | 1 | 1 | 1 |
| `raw:trpv1-like receptor-expressing arc pomc neurons` | `food_intake` | 1 | 1 | 1 |
| `raw:type 1 agrp neurons` | `bat_thermogenesis` | 1 | 1 | 1 |
| `raw:type 1 agrp neurons` | `blood_glucose` | 1 | 1 | 1 |
| `raw:type 1 agrp neurons` | `body_weight` | 1 | 1 | 1 |

... 6 more row(s) not shown.

## Leakage-aware splits available on this corpus

SCHEMA.md 5: a cross-validated number reported without naming its grouping is an invalid result in this project. `syn4rch.feasibility.splits` therefore returns folds only inside a labelled fold set, and `attach_metric` raises if a caller tries to bind a number to folds without the label. Shared animal cohorts (SCHEMA.md 5 item 3) are not controlled by any grouping here: the schema records no cohort identity.

| grouping | groups | usable folds | cross-validation possible | records excluded as ungroupable |
|---|---|---|---|---|
| `by_study` | 25 | 25 | yes | 0 |
| `by_lab` | 25 | 25 | yes | 0 |
| `leave_one_population_out` | 10 | 10 | yes | 0 |

## What this report does not claim

- It does not estimate how many records curation will eventually yield. Every number above is a count of what is in the input file at the time shown.
- It does not resolve, average, prefer or drop either side of a disagreement. Contested claims are reported as contested.
- It does not treat a single observation as agreement, and it does not treat repeated publication by one lab as independent replication.
- It cannot rule out two studies sharing one animal cohort (SCHEMA.md 5 item 3): the schema records no cohort identity.
- A `feasible` verdict says a prediction task can be posed and cross-validated with the grouping named, not that a model will work.
- Tier B outcomes are excluded as prediction targets by `EXCLUDE_TIER_B_AS_TARGET` (SCHEMA.md 6: Tier B is an annotation layer, not a training target).
