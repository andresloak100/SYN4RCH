# Schema migration v0.3.0 -> v0.4.0

v0.4.0 adds one field, `background_genotype`, as the fifth condition dimension (D-060). That is the whole change to the data: the field is **added, not filled**. Every record begins in the `None` state, which D-060 defines as "the curator has not addressed genotype here" and which `condition_key` pools with `not_reported` into `unknown`. Filling it is a separate curation pass.

## Run

| field | value |
|---|---|
| generated (UTC) | 2026-11-06T20:44:29Z |
| input | `../../../../var/folders/_5/rynzp98n41x8_zqylycjclv80000gn/T/tmp8kcv97si/records.jsonl` |
| from / to | 0.3.0 -> 0.4.0 |
| records | 1 |
| versions in | 0.3.0 x1 |
| already at 0.4.0 (no-op) | 0 |
| **record ids unchanged** | **1 of 1** |
| record ids changed | 0 |
| **identity map required** | **no** |

## No identity map, and why that is the finding

This is the first migration in this project that leaves no identity map. `background_genotype` is deliberately **not** a component of `RECORD_ID_FIELDS`, so adding it re-identifies nothing — every `record_id` below is byte-identical to its input, checked per record rather than assumed, and the pass refuses outright if one moves.

The three passes before this one each rewrote identities and left a map that every external reference must still be resolved through: `record_id_map_0.1.0_to_0.2.0.json`, `record_id_map_0.2.0_to_0.3.0.json`, `record_id_map_d053.json`. That debt is permanent. Avoiding a fourth was a design criterion in D-060, not a happy accident: the variant that would have put the genotype in the hash was scored against this one and rejected, and the measured comparison is in that entry.

## What changed, field by field

| field | records changed |
|---|---|
| `schema_version` | 1 |

`schema_version` is the stamp, and on a clean run it is the only row. `background_genotype` does **not** appear, which is correct and worth stating: the dataclass declares it with a default of `None`, so reading a v0.3.0 record already materialises it as null and there is no value for the migration to change. The field is added to the *schema*, not written to the *data*. Any other field in this table would be a defect.

## Hard findings, before and after

| when | hard findings |
|---|---|
| before | {'S006_SCHEMA_VERSION_MISMATCH': 1} |
| after | none |

A migration that introduces a hard finding has broken something. The two rows are printed rather than compared in prose so a reader can see which code it would have been.

## What this pass does not do

* It does not curate. `background_genotype` is null on every record afterwards, and the 71 records D-060 identifies as carrying an establishable background genotype are curated by a separate pass.
* It does not touch any other field, and the table above is the evidence rather than the claim.
* It does not renumber anything. `FROM_VERSION` and `TO_VERSION` are pinned literals rather than `SCHEMA_VERSION`, so this report still recomputes after the schema moves again (D-046).

