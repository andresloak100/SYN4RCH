# SYN4RCH Tier A Record Schema (v0.1.0)

This document is the **binding contract**. Code, extraction prompts and QC checks
must conform to it. Changes require a version bump and an entry in `docs/DECISIONS.md`.

## 1. Unit of record

One record = **one observation**:

    (neuronal population) x (perturbation) x (physiological readout) x (study)

A single paper yields many records. A single figure panel usually yields one record.
Records are never aggregated across populations, valences, or outcomes.

## 2. Anti-fabrication rules (non-negotiable)

1. `support_quote` is **mandatory** and must be text copied verbatim from the source.
2. The stated `direction` must be derivable from `support_quote` alone.
3. `source_locator` must name where it came from (e.g. "Fig 3C", "Results para 4", "Table 1").
4. If a field is not reported, it is `null`. Never inferred, never defaulted, never guessed.
5. No record may be created for a paper whose text was not actually retrieved.
6. `extraction_method` must honestly record whether a human reviewed the record.

A record failing any of these is invalid and must not enter `data/curated/`.

## 3. Fields

### 3.1 Identity and provenance
| field | type | notes |
|---|---|---|
| `record_id` | str | deterministic hash of (study_id, population_raw, perturbation, outcome_raw, source_locator) |
| `study_id` | str | canonical: `doi:...`; fallback `pmid:...` |
| `pmid` / `pmcid` / `doi` | str? | all that are known |
| `title`, `journal`, `year` | str/int | |
| `first_author`, `last_author`, `last_author_affiliation` | str? | **required for lab-level leakage control** |
| `retrieved_from` | enum | `pubmed`, `europepmc_oa`, `europepmc_abstract`, `crossref`, `manual` |
| `full_text_available` | bool | abstract-only records are marked and down-weighted |
| `schema_version` | str | `0.1.0` |
| `extraction_method` | enum | `human`, `llm`, `llm_human_reviewed` |
| `extracted_by`, `extracted_at` | str | |
| `review_status` | enum | `unreviewed`, `accepted`, `rejected`, `needs_full_text` |

### 3.2 Neuronal population
| field | type | notes |
|---|---|---|
| `population_raw` | str | exactly as the paper names it |
| `population_id` | str? | id from `vocab/populations.yaml`; `null` if no confident match |
| `marker_genes` | list[str] | e.g. `["Agrp","Npy"]`; HGNC/MGI symbols |
| `region_raw` | str | as written |
| `region_id` | str? | Allen Brain Atlas / HypoMap identifier where resolvable |
| `neurotransmitter` | list[enum] | `GABA`, `glutamate`, `dopamine`, `serotonin`, `histamine`, `peptidergic`, `mixed`, `unknown` — only if the paper or the vocab entry states it |
| `projection_specific` | bool | true if the manipulation targeted a projection, not the soma population |
| `projection_target_raw` / `projection_target_id` | str? | |

Populations not in the vocabulary are kept with `population_id = null` and surfaced
in a report for vocabulary extension. They are **not** force-mapped to the nearest name.

### 3.3 Subject
`species` (`mus_musculus`, `rattus_norvegicus`, `homo_sapiens`, other),
`strain`, `genetic_background`, `sex` (`male`,`female`,`both`,`not_reported`),
`age_value` + `age_unit` (`days`,`weeks`,`months`,`years`), `age_group`
(`juvenile`,`young_adult`,`middle_aged`,`aged`,`not_reported`),
`diet` (`chow`,`high_fat`,`other`,`not_reported`), `n_animals`.

`age_group` is assigned only by explicit rule from `age_value` (documented in code),
never by impression. Mouse: <8 wk juvenile, 8-26 wk young_adult, 27-64 wk middle_aged,
>=65 wk aged.

### 3.4 Perturbation
| field | type | notes |
|---|---|---|
| `method` | enum | `optogenetic`, `chemogenetic`, `electrical_stimulation`, `lesion`, `genetic_ablation`, `conditional_knockout`, `conditional_overexpression`, `local_pharmacology`, `other` |
| `actuator` | str? | `ChR2`, `hM3Dq`, `hM4Di`, `KORD`, `eNpHR`, `ArchT`, `Casp3`, `DTA`, `TeTx`, ... |
| `valence` | enum | `activation`, `inhibition`, `ablation`, `unclear` |
| `targeting_strategy` | str | Cre line + viral vector, knock-in, etc. |
| `acuity` | enum | `acute` (<24 h), `subchronic` (1-7 d), `chronic` (>7 d), `not_reported` |
| `duration_value`/`duration_unit` | | |
| `stim_frequency_hz`, `stim_pulse_width_ms`, `stim_duty_cycle`, `stim_power_mw`, `light_wavelength_nm` | float? | optogenetics only |
| `ligand`, `ligand_dose_value`, `ligand_dose_unit`, `ligand_route` | | chemogenetics only |
| `laterality` | enum | `unilateral`, `bilateral`, `not_reported` |
| `time_to_readout_value`/`_unit` | | delay between perturbation and measurement |

`valence = unclear` is a valid and expected answer (e.g. mixed Gq/Gi, lesion with
unclear net circuit effect). It must not be coerced.

### 3.5 Readout and effect
| field | type | notes |
|---|---|---|
| `outcome_raw` | str | as written |
| `outcome_id` | str | from `vocab/outcomes.yaml`; required |
| `outcome_tier` | enum | `A` (physiological) or `B` (aging/lifespan/healthspan) |
| `tissue_raw`/`tissue_id` | str? | measurement compartment |
| `assay` | str | e.g. "tail-vein glucometer", "indirect calorimetry", "ELISA" |
| `direction` | enum | `increase`, `decrease`, `no_change`, `non_monotonic`, `not_reported` |
| `effect_size_value`, `effect_size_type`, `effect_size_units` | | `percent_change`, `absolute_difference`, `fold_change`, `cohens_d`, `auc_difference` |
| `dispersion_type`, `dispersion_value`, `ci_low`, `ci_high` | | `sem`, `sd`, `ci95` |
| `n_treatment`, `n_control` | int? | |
| `statistical_test`, `p_value`, `significant` | | `significant` only if the paper states a threshold comparison |
| `is_primary_outcome` | bool | distinguishes headline claim from incidental measurement |

`no_change` is recorded explicitly. Dropping null results is a leading source of
direction bias and is forbidden.

### 3.6 Controls and evidence grading
| field | type | notes |
|---|---|---|
| `control_types` | list[enum] | `opsin_free`, `reporter_only`, `vehicle`, `ligand_in_nonexpressing`, `light_off`, `contralateral`, `cre_negative`, `wildtype`, `none_reported` |
| `has_actuator_control` | bool | the specific control that makes a chemo/optogenetic claim interpretable |
| `evidence_class` | enum | `causal_interventional`, `functional_association`, `anatomical` — Tier A records are `causal_interventional` by construction; the other two belong to the knowledge graph, not this table |
| `confidence` | enum | `high`, `moderate`, `low` — assigned by documented rule from (full text available, actuator control present, n reported, effect size reported), never by vibe |
| `conflicts_with` | list[record_id] | populated by QC, not by the extractor |
| `support_quote` | str | **mandatory verbatim** |
| `source_locator` | str | **mandatory** |
| `extraction_notes` | str? | ambiguities, caveats |

## 4. Storage layout and immutability

    data/
      raw/        harvest dumps, append-only, never edited by hand
      interim/    screening decisions, parsed texts
      curated/    records.jsonl  <- the dataset; version controlled
      inferred/   anything derived (KG edges, imputations)
      predictions/ model outputs only
      vocab/      controlled vocabularies

These four layers never mix. A value in `curated/` must trace to `raw/`.
A value in `inferred/` or `predictions/` may never be cited as evidence.

## 5. Leakage control (defines the benchmark)

Splits must be grouped so that none of the following cross the train/test boundary:

1. `study_id` — the obvious one.
2. **Lab identity** (`last_author` + affiliation) — the same group publishing the
   same circuit repeatedly is the single largest leakage risk in this corpus.
3. **Shared cohort** — multiple papers on one animal cohort; flagged manually.
4. **Population x outcome cell** — for generalization claims, evaluate
   leave-one-population-out as well as leave-one-study-out.

Reporting a single cross-validated number without naming the grouping is
considered an invalid result in this project.

## 6. What this schema deliberately refuses to represent

- Quantitative dose-response curves (insufficient reporting across the corpus).
- Inferred anatomical wiring presented as measured (belongs to the KG with its own class).
- Any "rejuvenation score". Tier B is an annotation layer, not a training target.
