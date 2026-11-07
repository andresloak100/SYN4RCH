# Condition fill (top-up), schema v0.4.0

Generated 2026-11-07T05:58:03Z by `python3 -m syn4rch.curate.fill_conditions run`.

Read 709 records from 102 papers. The pass reads `sex`, `diet` and `age_group` from each paper's Methods, and only for the dimensions that paper holds constant; it reads `challenge` from each record's own passage. Nothing is inferred: every value below is backed by a verbatim sentence and a locator stored on the record.

## What the papers could and could not settle

A dimension a paper holds constant can be applied to every record from that paper. A dimension that varies within the paper cannot be applied at all, because the Methods do not say which cohort an individual figure panel came from. `varies` is therefore not a failure of this pass; it is the paper declining to answer the question at the scope the question was asked.

| dimension | constant | varies | absent |
|---|---|---|---|
| `sex` | 47 | 33 | 22 |
| `diet` | 29 | 34 | 39 |
| `age_group` | 13 | 11 | 78 |
| `challenge` | 0 | 54 | 48 |
| `background_genotype` | 0 | 26 | 76 |

## What changed, per dimension

| dimension | filled | corrected | confirmed | conflicted |
|---|---|---|---|---|
| `sex` | 169 | 0 | 121 | 0 |
| `diet` | 85 | 0 | 104 | 0 |
| `age_group` | 45 | 0 | 36 | 0 |
| `challenge` | 1 | 1 | 29 | 6 |
| `background_genotype` | 0 | 0 | 0 | 0 |

`filled` was null. `corrected` said `not_reported`, which an extraction working from one passage could assert and the Methods then contradict. `confirmed` already held the value and gained only its evidence, so its `record_id` did not change. `conflicted` holds a different positive value and was **left exactly as curated**; those records are listed below for a human to settle.

## Before and after

* `sex`
  * before: null 494, male 121, female 43, not_reported 27, both 24
  * after:  null 325, male 251, both 62, female 44, not_reported 27
* `diet`
  * before: null 325, chow 152, not_reported 124, high_fat 107, other 1
  * after:  null 240, chow 237, not_reported 124, high_fat 107, other 1
* `age_group`
  * before: null 589, not_reported 77, young_adult 42, juvenile 1
  * after:  null 544, young_adult 87, not_reported 77, juvenile 1
* `challenge`
  * before: null 573, not_reported 41, other 34, fasted 21, ad_libitum 18, refed 14, cold_exposure 4, restraint_stress 4
  * after:  null 572, not_reported 40, other 34, fasted 22, ad_libitum 19, refed 14, cold_exposure 4, restraint_stress 4
* `background_genotype`
  * before: null 616, not_reported 37, npy:loss_of_function:cell_type_restricted 16, lep:loss_of_function:global 15, wild_type 15, mixed 7, asip:gain_of_function:global 2, trpv1:loss_of_function:global 1
  * after:  null 616, not_reported 37, npy:loss_of_function:cell_type_restricted 16, lep:loss_of_function:global 15, wild_type 15, mixed 7, asip:gain_of_function:global 2, trpv1:loss_of_function:global 1

## Which rule read which value

| rule | values |
|---|---|
| `methods_single_sex` | 238 |
| `methods_single_diet` | 189 |
| `methods_single_age_band` | 81 |
| `methods_both_sexes` | 52 |
| `quote_names_the_state` | 21 |
| `outcome_names_the_state` | 9 |
| `quote_names_the_state_with_a_competing_token` | 1 |

This table is the audit handle. A rule that turns out to be wrong is wrong for every value it produced, and this says how many that is and therefore how much of the dataset a correction has to touch.

## Identity

All four condition fields are `record_id` components (D-026), so a filled value re-identifies its record. 219 of 709 ids changed; 490 did not, those being the records whose values were confirmed rather than changed.

No two records collapsed onto one id.

## Evidence

479 of 709 records now carry at least one condition evidence entry, 743 entries in all. An entry means a sentence was found in the retrieved text and checked verbatim; it is not the same statement as a value of `not_reported`, which asserts that no such sentence exists, and the two are never added together.

## What the extractor refused, and why (469)

A refusal is a reading the rules declined rather than a sentence they failed to find. Each is named, so the set of records a rule affects can be re-examined without re-reading the corpus.

| refusal | records |
|---|---|
| `refused_methods_cannot_state_a_measurement_state` | 92 |
| `refused_methods_cannot_establish_one_background` | 76 |
| `refused_housing_statement_is_not_a_challenge` | 38 |
| `refused_protocol_start_age_is_not_measurement_age` | 37 |
| `refused_assay_implies_its_own_challenge` | 34 |
| `refused_figure_scoped_cohort_sentence` | 27 |
| `refused_methods_name_more_than_one_background` | 26 |
| `refused_age_band_requires_mouse` | 18 |
| `refused_disjunctive_diet_is_two_arms` | 16 |
| `refused_quote_not_at_its_locator` | 15 |
| `refused_quote_names_more_than_one_state` | 14 |
| `refused_unrecognised_diet_term` | 12 |
| `refused_age_not_constant_after_a_refusal` | 10 |
| `refused_straddling_age_range` | 10 |
| `refused_ambiguous_sex_sentence` | 9 |
| `refused_sex_not_constant_after_a_refusal` | 7 |
| `refused_second_sex_reported_outside_the_methods` | 6 |
| `refused_disjunctive_sex_is_unresolved` | 5 |
| `refused_restricted_feeding_has_no_enum_member` | 5 |
| `refused_rat_age_band` | 4 |
| `refused_tissue_preparation_age` | 4 |
| `refused_in_vitro_age_statement` | 2 |
| `refused_qualifier_inverts_the_state` | 1 |
| `refused_two_states_in_one_outcome` | 1 |

## Disagreements left unresolved (6 records)

Each of these holds a curated value that the paper's own text appears to contradict. Nothing was changed. A human has to read the paper and decide, and until then the record says what the original curation said.

* `cf668711e4369657` (doi:10.1073/pnas.1802237115)
  * challenge: curated 'refed', paper says 'fasted' (passage, Results para 5) -- left as curated
* `e527341d2a14a370` (doi:10.1038/s41467-026-75650-0)
  * challenge: curated 'refed', paper says 'fasted' (passage, Results para 13) -- left as curated
* `d206091ce714a8f6` (doi:10.1038/s41467-018-05711-6)
  * challenge: curated 'refed', paper says 'fasted' (passage, Results para 5) -- left as curated
* `bf6a25067c559d29` (doi:10.1016/j.molmet.2022.101455)
  * challenge: curated 'refed', paper says 'fasted' (passage, Results para 9) -- left as curated
* `e1e102a40065cdef` (doi:10.1016/j.celrep.2017.11.089)
  * challenge: curated 'refed', paper says 'fasted' (passage, Results para 6) -- left as curated
* `3f9e8f6e7b671752` (doi:10.1016/j.celrep.2017.11.089)
  * challenge: curated 'refed', paper says 'fasted' (passage, Results para 8) -- left as curated

## Does every quote sit where it says it does?

743 of 743 evidence entries were resolved to the passage their locator names and the quote found there. 0 failed.

This is a stricter check than the extractor's own gate, which asks only whether a quote occurs somewhere in the paper. A real quote carrying the wrong address sends an auditor to the wrong paragraph, and an audit trail that misdirects is worse than none, because it looks sound.

## Hard findings

* before: none
* after:  none

A pass that introduces a hard finding has written something the schema forbids and must not be committed.

