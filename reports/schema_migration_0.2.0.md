# Schema migration: curated records v0.1.0 -> v0.2.0

Recomputed 2026-11-04T20:36:28Z by `python3 -m syn4rch.curate.migrate report` from the committed dataset. No record was written.

Every number below is recomputed from `data/curated/records.0.2.0.jsonl.bak` by `python3 -m syn4rch.curate.migrate report` (`docs/CURATION_PROTOCOL.md` 6). Nothing here is estimated.

## 1. What moved

| | |
|---|---|
| records in | 362 |
| records out | 362 |
| records dropped or added | 0 |
| `record_id` changed | 362 |
| `record_id` unchanged | 0 |
| `conflicts_with` references rewritten | 0 of 0 |
| `conflicts_with` references not in the id map | 0 |
| distinct `record_id` | 362 |
| `record_id` consistent with its own fields | 362 of 362 |
| id map | `data/curated/record_id_map_0.1.0_to_0.2.0.json` |
| pre-migration file kept at | `data/curated/records.0.1.0.jsonl.bak` |

`RECORD_ID_FIELDS` grew from 7 components to 12, so **every id in the corpus changed**. That is the whole reason for the version bump. The five components added were `target_raw`, `sex`, `diet`, `age_group`, `challenge`.

The map from old id to new id is complete and one-to-one over the input:
every v0.1.0 id appears exactly once as a key, and no two of them share a
new id. Two records colliding on one new id would have merged two
observations into one, so the migration refuses to write at all rather
than writing the collision and mentioning it here.

Schema versions now present in the file: 0.2.0=362.

## 2. The perturbation axis

Filled by rule from fields the curators already wrote:
`perturbation_class` from `method`, `perturbation_polarity` from
`(method, valence)`. Both are pure functions
(`syn4rch.schema.rules`), so nothing here is a judgement.

### `perturbation_class`

| class | records | share |
|---|---|---|
| `circuit_activity` | 170 | 47.0% |
| `molecular` | 170 | 47.0% |
| `structural` | 17 | 4.7% |
| `pharmacological` | 5 | 1.4% |

### `perturbation_polarity`

| polarity | records | share |
|---|---|---|
| `increase` | 181 | 50.0% |
| `decrease` | 149 | 41.2% |
| `remove` | 17 | 4.7% |
| `unclear` | 15 | 4.1% |

### Polarity within class

Polarity read without its class is meaningless (SCHEMA.md 3.4b):
`decrease` is less firing under `circuit_activity` and less of the target
gene product under `molecular`.

| class | polarity | records |
|---|---|---|
| `circuit_activity` | `increase` | 100 |
| `circuit_activity` | `decrease` | 60 |
| `circuit_activity` | `unclear` | 10 |
| `molecular` | `decrease` | 89 |
| `molecular` | `increase` | 81 |
| `structural` | `remove` | 17 |
| `pharmacological` | `unclear` | 5 |

For reference, the `method` and `valence` these were derived from: method chemogenetic=144, conditional_knockout=89, conditional_overexpression=81, optogenetic=26, genetic_ablation=17, local_pharmacology=5; valence unclear=185, activation=100, inhibition=60, ablation=17.

## 3. Molecular targets

The target was never missing from this dataset; it was curated inside the free text of `targeting_strategy` and sometimes `actuator`. This migration lifted it into `target_raw` / `target_id` against a new vocabulary, `data/vocab/targets.json` (17 entries, 27 matchable strings). **It did not read any paper.** A record's target is only as good as its `targeting_strategy`: a curator's wrong gene there becomes a confidently resolved wrong gene here, and nothing in this pipeline would notice.

| | |
|---|---|
| eligible records (`molecular` or `pharmacological`) | 175 |
| resolved | 170 |
| unresolved | 5 |
| of which, no synonym matched | 5 |
| of which, ambiguous (two or more genes matched) | 0 |
| already curated, left untouched | 0 |
| ineligible records carrying a target (must be 0; `S008`) | 0 |

Re-running the resolver over the committed file reproduces the stored `target_id` for 170 of 170 resolved records, with 0 disagreement(s) and 0 ambiguity(ies). A disagreement would mean the stored target no longer follows from the text it was taken from.

### Ambiguous matches, for a human to curate

None: no record's curated text names two different genes of the
vocabulary. Had one, it would be listed here by `record_id` and left
unresolved; an ambiguous match is never guessed.

### Eligible but unresolved

| record_id | method | targeting_strategy (truncated) |
|---|---|---|
| `54b55771e8efacb8` | `local_pharmacology` | capsaicin directly infused into the ARC |
| `8bafb2b422520641` | `local_pharmacology` | intra-ARC microinjection (30 pmol, 0.2 uL) |
| `a6252beb7bf4bbf1` | `local_pharmacology` | intra-ARC microinjection (30 pmol, 0.2 uL) |
| `8c66db83a4be42d7` | `local_pharmacology` | intra-ARC microinjection (30 pmol, 0.2 uL) |
| `bb183415b6f72867` | `local_pharmacology` | intra-ARC microinjection (30 pmol, 0.2 uL) |

`data/vocab/targets.json` is a vocabulary of *genes*, so a locally infused drug resolves to nothing in it. Naming the receptor a drug is thought to act on would be inference about pharmacology, not resolution of a curated string, and this code does not do it. A pharmacological record's compound stays where the curator put it, in `ligand`.

### Targets found

| target_id | records | studies | labs |
|---|---|---|---|
| `trp53` | 34 | 1 | 1 |
| `sh2b1` | 20 | 1 | 1 |
| `sirt1` | 18 | 1 | 1 |
| `trim31` | 14 | 1 | 1 |
| `htr2c` | 11 | 1 | 1 |
| `bcl6` | 10 | 1 | 1 |
| `fam172a` | 9 | 1 | 1 |
| `smad7` | 9 | 1 | 1 |
| `prkce` | 8 | 1 | 1 |
| `npy` | 7 | 1 | 1 |
| `raf1` | 7 | 1 | 1 |
| `rptor` | 7 | 1 | 1 |
| `trsp` | 7 | 1 | 1 |
| `kiss1r` | 3 | 1 | 1 |
| `il17ra` | 2 | 1 | 1 |
| `manf` | 2 | 1 | 1 |
| `socs3` | 2 | 1 | 1 |

## 4. The gene is the study is the lab

**The 17 target genes in this corpus come from 16 studies and 16 laboratories, and 17 of the 17 appear in exactly one study while 17 appear in exactly one laboratory: the gene, the study and the lab are very nearly the same variable.** A model told the gene has therefore been told the study and the lab, so a cross-validated number that does not group on study or lab is measuring how well 17 single-gene experiments can be memorised, not whether a molecular perturbation predicts a physiological outcome.

This is the single biggest structural fact about the molecular half of the corpus, and the reason `target_id` had to exist: without it, "molecular perturbation of this population" pools every one of those experiments into one cell, and the cell then looks like a replicated finding. With it, the collinearity is visible instead of hidden -- it is not removed. Nothing a schema can do makes one laboratory's single gene generalise.

Studies contributing more than one gene:

* `doi:10.1016/j.isci.2024.110259`: `il17ra`, `socs3`

Leakage control is specified in SCHEMA.md 5: group on `study_id` and on
lab identity. Under this corpus those groupings also group the gene, so
leave-one-study-out on molecular records is, in effect,
leave-one-gene-out, and no split of this dataset can test whether a
molecular perturbation of one gene predicts another.

## 5. The condition dimension was not filled, and could not be

**The migration was not permitted to fill `challenge`, `challenge_detail`, `condition_source` or `condition_locator`, and did not fill them.** Those four require reading the paper, and this code read no paper. It did not write `not_reported` either: per SCHEMA.md 3.3b that means a curator checked the paper and it does not say, so writing it from a migration would record a curation act that never happened (SCHEMA.md 2 rule 4). They are left `null` and counted here.

| field | set | null | of which `not_reported` |
|---|---|---|---|
| `sex` | 138 | 224 | 22 |
| `diet` | 165 | 197 | 26 |
| `age_group` | 76 | 286 | 69 |
| `challenge` | 0 | 362 | 0 |
| `challenge_detail` | 0 | 362 | 0 |
| `condition_source` | 0 | 362 | 0 |
| `condition_locator` | 0 | 362 | 0 |

1448 of the 1448 condition-field slots the v0.2.0 schema added are unfilled (362 records x 4 fields). `sex`, `diet` and `age_group` are curated where known and were not touched.

One consequence is worth stating in advance, because it will look like a defect when it happens. `challenge` is a `record_id` component, so the curation pass that fills it will change the id of every record it touches, exactly as this migration changed all 362. That is the identity model working as designed -- a fasted cohort and an ad-libitum one are two observations -- but it means condition curation is an id-changing act and needs the same map-and-remap treatment this module implements, not an in-place edit.

The consequence is in the dataset's own warnings: `W015` fires on the 153 records for which none of `(sex, diet, age_group, challenge)` is known, and those records can only ever be compared with other condition-blind records. Closing that gap is curation work, not migration work, and no amount of migration will do it.

## 6. Validation, before and after

The migration refuses to write if the result would carry a hard validation finding the input did not already have, and the two tables below are how that is checked. `S005` (id inconsistent) and `S006` (version mismatch) are expected on every *input* record, because at v0.2.0 the seven-component ids no longer hash and the old version string no longer matches: they are the artefacts the migration exists to remove, and they are excluded when the two sides are compared.

### After the migration, recomputed from the committed dataset

| finding | severity | records |
|---|---|---|
| `S006_SCHEMA_VERSION_MISMATCH` | HARD | 362 |
| `W001_NO_CONTROL_TYPE` | WARNING | 152 |
| `W002_NO_SAMPLE_SIZE` | WARNING | 352 |
| `W003_POPULATION_ID_NULL` | WARNING | 67 |
| `W005_UNREVIEWED_LLM_EXTRACTION` | WARNING | 362 |
| `W006_NO_EFFECT_SIZE` | WARNING | 191 |
| `W008_SIGNIFICANCE_WITHOUT_TEST` | WARNING | 84 |
| `W010_NO_LAB_IDENTITY` | WARNING | 10 |
| `W015_CONDITION_UNSPECIFIED` | WARNING | 153 |
| `W016_POLARITY_MAY_BE_INVERTED_BY_ALLELE` | WARNING | 8 |
| `W017_CONDITION_WITHOUT_EVIDENCE` | WARNING | 209 |

362 of 362 records carry no hard finding.

### Before the migration, recomputed from the preserved input

Counted by running `syn4rch.schema.validate` over the v0.1.0 file at the v0.2.0 code's standards. Checkable: the pre-migration file is kept, and `migrate report` re-reads it.

| finding | severity | records |
|---|---|---|
| `S005_RECORD_ID_INCONSISTENT` | HARD | 362 |
| `S006_SCHEMA_VERSION_MISMATCH` | HARD | 362 |
| `W001_NO_CONTROL_TYPE` | WARNING | 152 |
| `W002_NO_SAMPLE_SIZE` | WARNING | 352 |
| `W003_POPULATION_ID_NULL` | WARNING | 67 |
| `W005_UNREVIEWED_LLM_EXTRACTION` | WARNING | 362 |
| `W006_NO_EFFECT_SIZE` | WARNING | 191 |
| `W008_SIGNIFICANCE_WITHOUT_TEST` | WARNING | 84 |
| `W010_NO_LAB_IDENTITY` | WARNING | 10 |
| `W014_MOLECULAR_TARGET_UNIDENTIFIED` | WARNING | 170 |
| `W015_CONDITION_UNSPECIFIED` | WARNING | 153 |
| `W016_POLARITY_MAY_BE_INVERTED_BY_ALLELE` | WARNING | 8 |
| `W017_CONDITION_WITHOUT_EVIDENCE` | WARNING | 209 |

362 of 362 input records carried no hard finding at all, and 362 of 362 do now.

Setting `S005` and `S006` aside, the input carried no hard finding and the migrated dataset carries no hard finding. The migration therefore introduced no hard finding, which is the condition `migrate run` checks before it will write anything.

`W014` (molecular perturbation with no identified target) is the one warning this migration was able to retire, and it did so by resolving a field the curators had already written, not by learning anything new about any paper.

**This report describes the dataset as it stood at v0.2.0.** The v0.3.0 condition fill has since rewritten `data/curated/records.jsonl` and re-identified 227 of its 362 records, so `migrate report` can no longer recompute these numbers from the committed dataset: it refuses, naming the version mismatch. Recomputing them needs the preserved v0.2.0 state (`data/curated/records.0.2.0.jsonl.bak`, untracked like every backup), and `data/curated/record_id_map_0.2.0_to_0.3.0.json` carries any id in this report forward to the record it is now. This is the ordinary consequence of a second migration, and it is stated here rather than left for a reader to discover from an error message.

**Two notes on reading the counts above.** They are what the *current* validator sees in these two files, because a report whose numbers cannot be recomputed from the committed dataset does not belong in this repository (`docs/CURATION_PROTOCOL.md` 6). So a rule introduced after this migration appears in its table: `W017`, added by schema v0.3.0, counts condition values that carry no verbatim sentence, and it describes a gap this migration never claimed to close.

And the "carried no hard finding at all" figure sets `S005` and `S006` aside, as every other before/after comparison here does. An earlier rendering did not, which made it read **0 of 362** for the input -- a number produced entirely by a 0.1.0 file being validated by newer code, and one that contradicted the sentence immediately below it. The migration neither introduced nor removed a hard finding, and both columns now say so.

## 7. What the migration filled, and what it refused to

| v0.2.0 field | filled by the migration? | how |
|---|---|---|
| `perturbation_class` | yes | total function of `method` |
| `perturbation_polarity` | yes | function of `(method, valence)` |
| `target_raw` | only when exactly one gene occurs in the curated text | resolution of `targeting_strategy` / `actuator` |
| `target_id` | as `target_raw` | `data/vocab/targets.json` |
| `challenge` | **no** | requires reading the paper |
| `challenge_detail` | **no** | requires reading the paper |
| `condition_source` | **no** | requires reading the paper |
| `condition_locator` | **no** | requires reading the paper |
| `schema_version` | yes | set to `0.2.0` |
| `record_id` | yes | recomputed over 12 components |
| `conflicts_with` | yes | rewritten through the id map |
| everything else | no | untouched |

`docs/DECISIONS.md` D-017 makes the gate the only writer of the curated file. This migration is the one other writer, and a narrower one: it creates no record, admits none, drops none, and writes no value that a curator did not already write or that a published rule does not derive from one. D-017 should be amended to name the exception.

