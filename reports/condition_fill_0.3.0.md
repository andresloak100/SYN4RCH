# Condition fill, schema v0.2.0 -> v0.3.0

Generated 2026-11-04T20:32:30Z by `python3 -m syn4rch.curate.fill_conditions run`.

Read 362 records from 42 papers. The pass reads `sex`, `diet` and `age_group` from each paper's Methods, and only for the dimensions that paper holds constant; it reads `challenge` from each record's own passage. Nothing is inferred: every value below is backed by a verbatim sentence and a locator stored on the record.

## What the papers could and could not settle

A dimension a paper holds constant can be applied to every record from that paper. A dimension that varies within the paper cannot be applied at all, because the Methods do not say which cohort an individual figure panel came from. `varies` is therefore not a failure of this pass; it is the paper declining to answer the question at the scope the question was asked.

| dimension | constant | varies | absent |
|---|---|---|---|
| `sex` | 23 | 7 | 12 |
| `diet` | 14 | 22 | 6 |
| `age_group` | 5 | 3 | 34 |
| `challenge` | 0 | 31 | 11 |

## What changed, per dimension

| dimension | filled | corrected | confirmed | conflicted |
|---|---|---|---|---|
| `sex` | 112 | 7 | 55 | 2 |
| `diet` | 86 | 11 | 7 | 0 |
| `age_group` | 36 | 0 | 0 | 0 |
| `challenge` | 25 | 0 | 0 | 0 |

`filled` was null. `corrected` said `not_reported`, which an extraction working from one passage could assert and the Methods then contradict. `confirmed` already held the value and gained only its evidence, so its `record_id` did not change. `conflicted` holds a different positive value and was **left exactly as curated**; those records are listed below for a human to settle.

## Before and after

* `sex`
  * before: null 224, male 67, female 38, not_reported 22, both 11
  * after:  male 175, null 112, female 38, both 22, not_reported 15
* `diet`
  * before: null 197, high_fat 89, chow 50, not_reported 26
  * after:  chow 147, null 111, high_fat 89, not_reported 15
* `age_group`
  * before: null 286, not_reported 69, young_adult 6, juvenile 1
  * after:  null 250, not_reported 69, young_adult 42, juvenile 1
* `challenge`
  * before: null 362
  * after:  null 337, fasted 16, ad_libitum 6, cold_exposure 2, refed 1

## Which rule read which value

| rule | values |
|---|---|
| `methods_single_sex` | 160 |
| `methods_single_diet` | 104 |
| `methods_single_age_band` | 36 |
| `quote_names_the_state` | 16 |
| `methods_both_sexes` | 14 |
| `outcome_names_the_state` | 8 |
| `quote_names_the_state_with_a_competing_token` | 1 |

This table is the audit handle. A rule that turns out to be wrong is wrong for every value it produced, and this says how many that is and therefore how much of the dataset a correction has to touch.

## Identity

All four condition fields are `record_id` components (D-026), so a filled value re-identifies its record. 227 of 362 ids changed; 135 did not, those being the records whose values were confirmed rather than changed.

No two records collapsed onto one id.

## Evidence

252 of 362 records now carry at least one condition evidence entry, 339 entries in all. An entry means a sentence was found in the retrieved text and checked verbatim; it is not the same statement as a value of `not_reported`, which asserts that no such sentence exists, and the two are never added together.

## What the extractor refused, and why (164)

A refusal is a reading the rules declined rather than a sentence they failed to find. Each is named, so the set of records a rule affects can be re-examined without re-reading the corpus.

| refusal | records |
|---|---|
| `refused_methods_cannot_state_a_measurement_state` | 40 |
| `refused_assay_implies_its_own_challenge` | 26 |
| `refused_protocol_start_age_is_not_measurement_age` | 18 |
| `refused_quote_not_at_its_locator` | 15 |
| `refused_disjunctive_diet_is_two_arms` | 12 |
| `refused_housing_statement_is_not_a_challenge` | 9 |
| `refused_age_band_requires_mouse` | 7 |
| `refused_figure_scoped_cohort_sentence` | 7 |
| `refused_unrecognised_diet_term` | 6 |
| `refused_straddling_age_range` | 5 |
| `refused_ambiguous_sex_sentence` | 4 |
| `refused_quote_names_more_than_one_state` | 4 |
| `refused_age_not_constant_after_a_refusal` | 3 |
| `refused_restricted_feeding_has_no_enum_member` | 3 |
| `refused_disjunctive_sex_is_unresolved` | 1 |
| `refused_qualifier_inverts_the_state` | 1 |
| `refused_rat_age_band` | 1 |
| `refused_tissue_preparation_age` | 1 |
| `refused_two_states_in_one_outcome` | 1 |

## Disagreements left unresolved (2 records)

Each of these holds a curated value that the paper's own text appears to contradict. Nothing was changed. A human has to read the paper and decide, and until then the record says what the original curation said.

* `1f2c0189ace296a2` (doi:10.1007/s00125-013-3140-5)
  * sex: curated 'female', paper says 'male' (methods, Methods para 3) -- left as curated
* `6153f87b799a9b76` (doi:10.1007/s00125-013-3140-5)
  * sex: curated 'female', paper says 'male' (methods, Methods para 3) -- left as curated

## Does every quote sit where it says it does?

339 of 339 evidence entries were resolved to the passage their locator names and the quote found there. 0 failed.

This is a stricter check than the extractor's own gate, which asks only whether a quote occurs somewhere in the paper. A real quote carrying the wrong address sends an auditor to the wrong paragraph, and an audit trail that misdirects is worse than none, because it looks sound.

## Hard findings

* before: S006_SCHEMA_VERSION_MISMATCH 362
* after:  none

A pass that introduces a hard finding has written something the schema forbids and must not be committed.

