# SYN4RCH feasibility report (Phase 2)

**Verdict: insufficient evidence — `data/curated/records.jsonl` does not exist yet (0 records), so none of the five Phase 2 questions can be answered and no prediction task is modellable.**

Every number below is counted from the input named in the Run table. Nothing is estimated, imputed or projected. An `insufficient` verdict is a measurement of the corpus, not a failure of this module.

## Run

| field | value |
|---|---|
| generated (UTC) | 2026-11-03T20:14:21Z |
| input | `data/curated/records.jsonl` |
| input present | no (treated as 0 records) |
| records read | 0 |
| records usable as observations | 0 |
| records held out of the matrix | 0 |
| distinct studies | 0 |
| distinct labs (last author + affiliation) | 0 |
| usable records with no lab identity | 0 |
| load errors | 0 |

## The five Phase 2 questions

**None of the five can be answered: the input file does not exist yet, so there are 0 usable observations.** No table below is computed over more than zero records, and no derived statistic is printed for them: a sparsity, an agreement rate or a mean over an empty corpus would be arithmetic on nothing.

| question | measured answer |
|---|---|
| 1. populations with usable intervention data | 0 |
| 2. physiological outcomes measured | 0 |
| 3. relationships with independent replication | 0 |
| 4. direction agreement among independent studies | undefined (no claim has two independent studies) |
| 5. outcomes with sufficient evidence for modelling | none |

The promotion gates these answers were judged against are listed below, so the verdict can be re-checked when the corpus grows.

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

## What this report does not claim

- It does not estimate how many records curation will eventually yield. Every number above is a count of what is in the input file at the time shown.
- It does not resolve, average, prefer or drop either side of a disagreement. Contested claims are reported as contested.
- It does not treat a single observation as agreement, and it does not treat repeated publication by one lab as independent replication.
- It cannot rule out two studies sharing one animal cohort (SCHEMA.md 5 item 3): the schema records no cohort identity.
- A `feasible` verdict says a prediction task can be posed and cross-validated with the grouping named, not that a model will work.
- Tier B outcomes are excluded as prediction targets by `EXCLUDE_TIER_B_AS_TARGET` (SCHEMA.md 6: Tier B is an annotation layer, not a training target).
