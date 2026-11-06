# Background genotype fill (D-060, result in D-062)

Schema v0.4.0 declared `background_genotype` and left it null on all 508 records. This pass fills it from a declared table, reviewable in a diff, produced by eight readers — one per study carrying an adversarially confirmed background genotype — each deciding every record of its study under the arm rule.

## Run

| field | value |
|---|---|
| generated (UTC) | 2026-11-06T21:08:12Z |
| input | `data/curated/records.jsonl` |
| assignments table | `data/curated/background_genotype_assignments.json` |
| records | 508 |
| **filled** | **77** |
| refused | 0 |
| evidence entries added | 51 |
| quotes re-checked at their locator | 453 |
| locator failures | 0 |
| **record ids unchanged** | **508 of 508** |
| identity map required | no |

Outcomes: `filled` 77.

## The dimension, before and after

| value | before | after | delta |
|---|---|---|---|
| `mixed` | 0 | 5 | +5 |
| `not_reported` | 0 | 26 | +26 |
| `null` | 508 | 431 | -77 |
| `segment(s)` | 0 | 31 | +31 |
| `wild_type` | 0 | 15 | +15 |

`null` is not a gap to be closed. D-060 defines it as "the curator has not addressed genotype here", and `condition_key` pools it with `not_reported` into `unknown`. A record absent from the table keeps it.

## What `not_reported` and `wild_type` mean here, and why they differ

`not_reported` asserts that the paper was read for this record and establishes no background manipulation. It carries **no evidence entry**, because an absence cannot be quoted (A009, A011), and this pass refuses a `not_reported` row that arrives with a quote rather than dropping it silently.

`wild_type` is the opposite: a *positive, quoted* establishment that the animals carry nothing beyond the recorded perturbation. It does **not** pool with `unknown`, so an unstated background and an established wild-type background never share a claim cell (D-044). A background strain is not sufficient for it — a strain says what the animals are, not that nothing further was done to them.

`mixed` is a positive value for a record whose own quote reports a joint result across two or more genotypes. It is a fact about the record rather than an absence, and the record must **not** be split into one per genotype: deriving two observations from a sentence reporting one is what D-044 refuses in general.

## Hard findings, before and after

| when | hard findings |
|---|---|
| before | none |
| after | none |

## What this pass does not do

* It decides nothing. Every value is a row of the declared table, carrying its own quote, locator and reasoning.
* It does not overwrite. A record whose `background_genotype` is already set is refused, not updated — that is `correct_conditions`'s job and D-039's rule.
* It re-identifies nothing. `background_genotype` is not a `RECORD_ID_FIELDS` component (D-060), the pass checks that per record, and it refuses rather than writing an id map it does not have.

