# Condition corrections (D-053)

`fill_conditions` fills an empty condition dimension and refuses to overwrite one that is not empty (D-039). That refusal is right, and it left seven records this pass corrects: six storing `fasted` for a fast-refeed intake test, and one whose `challenge` came from a comparator rather than from the animals measured.

## Run

| field | value |
|---|---|
| generated (UTC) | 2026-11-06T17:33:33Z |
| input | `data/curated/records.jsonl` |
| schema version | 0.3.0 |
| records | 508 |
| corrections declared | 7 |
| applied | 7 |
| refused | 0 |
| record ids changed | 7 |
| condition quotes re-checked at their locator | 402 |
| locator failures | 0 |

Outcomes: `applied` 7.

## What changed

| old id | new id | dimension | from | to | outcome |
|---|---|---|---|---|---|
| `0397b7caf40ad06b` | `cf668711e4369657` | challenge | `fasted` | `refed` | applied |
| `087960c4021bf6f2` | `e527341d2a14a370` | challenge | `fasted` | `refed` | applied |
| `41fca07890c238c8` | `d206091ce714a8f6` | challenge | `fasted` | `refed` | applied |
| `fcea91b0669ad33f` | `bf6a25067c559d29` | challenge | `fasted` | `refed` | applied |
| `93c2e8f0c13a30ed` | `e1e102a40065cdef` | challenge | `fasted` | `refed` | applied |
| `fdb78695221a531e` | `3f9e8f6e7b671752` | challenge | `fasted` | `refed` | applied |
| `c1b98d981dba36bb` | `1068b48a9cf28019` | challenge | `fasted` | `null` | applied |

### Why, one by one

* `0397b7caf40ad06b` — food intake is the readout, so the animals had food at the moment of measurement; 'in fasted mice' describes the state they were in before it was returned
* `087960c4021bf6f2` — deprive, return food, measure intake: the measurement is of the refeeding
* `41fca07890c238c8` — fasting-induced hyperphagia is eating, which cannot be measured in an animal that still has no food
* `fcea91b0669ad33f` — deprive, return food, measure intake: the measurement is of the refeeding
* `93c2e8f0c13a30ed` — feeding after a fast is the refed state, not the fasted one
* `fdb78695221a531e` — 'in fasted mice' names the state before food was returned; the readout is the feeding that followed
* `c1b98d981dba36bb` — 'fasted mice' is the comparator the effect is measured against, not the state of the animals measured; the quote supports no challenge for this record

## The `challenge` distribution, before and after

| value | before | after | delta |
|---|---|---|---|
| `ad_libitum` | 18 | 18 | +0 |
| `cold_exposure` | 4 | 4 | +0 |
| `fasted` | 28 | 21 | -7 |
| `null` | 443 | 444 | +1 |
| `other` | 7 | 7 | +0 |
| `refed` | 8 | 14 | +6 |

## What this pass does not do

* It applies no judgement of its own. Every change is a row of `CORRECTIONS` in `syn4rch.curate.correct_conditions`, naming the record, the value that had to be there, the value to put there, the quote fragment that had to match and the decision that authorises it.
* It does not re-read the papers for new conditions. That is `fill_conditions`, which only fills what is empty.
* It does not touch any other field. `challenge` is a `RECORD_ID_FIELDS` component, so the identity moves with it and the id map is how every external reference follows.

