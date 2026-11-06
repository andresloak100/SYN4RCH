# SYN4RCH Tier A Record Schema (v0.3.0)

This document is the **binding contract**. Code, extraction prompts and QC checks
must conform to it. Changes require a version bump and an entry in `docs/DECISIONS.md`.

**v0.3.0 gave the condition dimension evidence, and nothing else.** Every filled
`sex`, `diet`, `age_group` and `challenge` now carries a `condition_evidence`
entry: a verbatim sentence, the passage it sits in, and the name of the rule that
read it (3.3b, D-037). The record-level `condition_source`/`condition_locator`
pair it replaces could describe one provenance where a record genuinely has two —
the animals' sex is stated once in Methods, the acute state of a measurement only
by the passage the record came from — and so described neither; it was null on all
362 records of the first corpus.

The claim key, the perturbation axis and the `record_id` *input set* are unchanged
from v0.2.0. The condition **values** changed, and because all four are hash
components (D-026), the curation pass that filled them re-identified every record
it touched; `data/curated/record_id_map_0.2.0_to_0.3.0.json` is the map and
`reports/condition_fill_0.3.0.md` is what it did, refusal by refusal. Four hard
rules (A008-A011) make an entry auditable or inadmissible, and one warning (W017)
counts the values still carrying no sentence rather than rejecting them (D-042).

**v0.2.0 changed two things and every `record_id` in the corpus.** It added the
condition dimension (3.3b) and split the perturbation axis (3.4b). The reason was
measured rather than assumed: under v0.1.0 the first 362 curated records produced
36 cross-study claims holding a conflict, 12 of them in direct `increase` vs
`decrease` opposition. Under v0.2.0 those become 9 and 1. They had been a male
cohort compared against a female one, a chow-fed cohort against a high-fat-fed one,
or a gene knockout against a gene overexpression, pooled into one claim by a key too
coarse to tell them apart.

The same change cuts claims reported by more than one study from 39 to 11, and that
cost is the other half of the finding: most of the apparent replication in this
corpus was an artifact of the same coarseness that manufactured the conflicts.
See `docs/DECISIONS.md` D-023 to D-028.

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
| `record_id` | str | deterministic hash; the twelve components are named in `RECORD_ID_FIELDS` and fixed by D-002 + D-026 |
| `study_id` | str | canonical: `doi:...`; fallback `pmid:...` |
| `pmid` / `pmcid` / `doi` | str? | all that are known |
| `title`, `journal`, `year` | str/int | |
| `first_author`, `last_author`, `last_author_affiliation` | str? | `last_author` is **the lab identity** used for leakage control (D-056: a laboratory is a person). The affiliation is free text and is *not* part of that identity; its one job is the homonym check `QC_POSSIBLE_LAB_HOMONYM`. Missing author → `W010`; missing affiliation → `W018` |
| `retrieved_from` | enum | `pubmed`, `europepmc_oa`, `europepmc_abstract`, `crossref`, `manual` |
| `full_text_available` | bool | abstract-only records are marked and down-weighted |
| `schema_version` | str | `0.2.0` |
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

### 3.3b Condition (v0.2.0)

| field | type | notes |
|---|---|---|
| `challenge` | enum? | `ad_libitum`, `fasted`, `refed`, `cold_exposure`, `thermoneutral`, `glucose_load`, `insulin_load`, `high_fat_challenge`, `restraint_stress`, `exercise`, `other`, `not_reported` |
| `challenge_detail` | str? | the paper's own wording, e.g. "16 h overnight fast" |
| `condition_evidence` | list[obj] | one entry per **filled** dimension (v0.3.0) |

Each `condition_evidence` entry carries `dimension` (one of `sex`, `diet`,
`age_group`, `challenge`), `value` (the enum value the entry supports), `source`
(`passage` or `methods`), `locator`, `quote` and `rule`. `locator` and `quote`
are both mandatory, and `value` must equal the field it describes; A008, A009,
A010 and A011 enforce those four things separately, because a missing locator, a
missing quote, a malformed entry and an entry that contradicts its own field are
four different mistakes with four different fixes.

A dimension recorded as `not_reported` carries **no** entry. That value asserts
an absence, and an absence has no sentence to quote. It is the one condition
value in the schema that is a claim about the paper rather than a reading of it,
which is also why it is the one value a wider-scope reading is allowed to
replace (D-039).

`rule` names the extraction rule that produced the entry. It exists so that a
rule later found to be wrong can be traced to every value it produced, and the
correction applied once to all of them rather than record by record. Nothing
reads it to make a decision.

`challenge` is the acute physiological state at the moment of measurement, and is
distinct from `diet`, which is the maintenance diet over weeks. In this literature
the challenge frequently decides the sign of the result: AgRP activation raises food
intake in a fed mouse and raises it much further in a fasted one, and a glucose
reading means opposite things before and after a glucose load.

**The condition dimension is `(sex, diet, age_group, challenge)`.** Two records whose
condition dimension differs are **not comparable**: they may disagree because the
biology disagrees, or because the experiments were never the same experiment. The
dataset cannot distinguish those, so it refuses to.

`condition_evidence` exists because rule 4 of section 2 and the way papers are
written are in direct tension. A result lives in Results; the diet and the sex
live in Methods. Forbidding the Methods route left 197 of the first 362 records
with no diet and 224 with no sex, which made the condition gap unfixable by
curating harder. `source = methods` permits the Methods route **and makes it
auditable**: the entry must name the Methods passage and quote the sentence, and
the verbatim gate applies to that quote as it does to any other. It is never a
licence to infer: a condition that no passage of the paper states is still
`null`.

The provenance is per dimension rather than per record (v0.3.0, D-037) because a
single record genuinely draws its dimensions from different places. The sex of
the animals is stated once in Methods and holds for every figure in the paper;
the acute state at the moment of measurement is stated, when it is stated at
all, by the passage the record was read from. A single record-level
`condition_source` could describe only one of those, so it described neither,
and it was null on all 362 records of the first corpus.

**What a Methods reading may and may not settle** (D-038). A dimension a paper
holds *constant* across every experiment it reports can be read from Methods and
applied to each of that paper's records. A dimension that *varies* within the
paper cannot be applied at all: a paper comparing a chow group against a
high-fat group states both diets in Methods, and nothing in Methods says which
group any individual figure panel belongs to. `varies` is not a failure to
extract. It is the paper declining to answer the question at the scope the
question was asked, and the honest record of that is `null`.

`not_reported` and `null` are different statements and both are kept. `null` means
the curator did not address the question; `not_reported` means the paper was checked
and does not say. For *grouping* they are pooled, because both mean the condition
could be anything (`syn4rch.qc.keys.condition_key`), and `condition_specificity`
reports how often that pooling happened.

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

### 3.4b Perturbation class, polarity and target (v0.2.0)

| field | type | notes |
|---|---|---|
| `perturbation_class` | enum? | `circuit_activity`, `molecular`, `structural`, `pharmacological`, `other`. Derived from `method`; admits no override |
| `perturbation_polarity` | enum? | `increase`, `decrease`, `remove`, `bidirectional`, `unclear`. The sign **along the class's own axis** |
| `target_raw` | str? | the paper's own wording for the molecular target, e.g. `Raptor` |
| `target_id` | str? | resolved against `data/vocab/targets.json` |

`valence` answers one question: did the manipulation drive this population's firing
up or down? For a conditional knockout of a receptor gene the honest answer is
`unclear`, and section 3.4 rightly forbids coercing it. Under v0.1.0 that left 185
of 362 records with no signed axis at all, 170 of them molecular, and so half the
dataset could not be a model input.

v0.2.0 stops asking `valence` to carry the molecular manipulations. `valence` keeps
its meaning exactly; `perturbation_class` and `perturbation_polarity` carry the rest.

**Polarity is meaningless without its class.** `decrease` is less firing under
`circuit_activity` and less of the target gene product under `molecular`. Those are
different physical interventions, and pooling them is the same sign error that
D-019 was written against. The pairing is enforced structurally:
`syn4rch.qc.keys.perturbation_key` returns the class, the polarity and the target as
one value, and nothing in this project may key on polarity alone.

`target_raw` / `target_id` follow the raw-then-resolved pattern already used for
populations and outcomes. They matter more than their size suggests. The 170
molecular records cover **17 distinct target genes across 16 studies and 16
laboratories, and no gene is manipulated by more than one laboratory.** Without a
target field, "molecular perturbation of AgRP neurons" reads as 89 observations of
one thing when it is seventeen unrelated single-gene experiments. A molecular record with no identified
target is admitted but warned about (`W014`), because it pools with every other
molecular perturbation of the same population and that is almost always wrong.

Only `molecular` and `pharmacological` perturbations may carry a target. A target on
a circuit or structural perturbation is refused (`S008`): it corrupts the
perturbation key, which is a component of `record_id`.

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
2. **Lab identity** (the normalised `last_author`) — the same group publishing the
   same circuit repeatedly is the single largest leakage risk in this corpus.
   This was `last_author` **+ affiliation** until D-056, and that key did not
   work: the affiliation is free text, the same department is spelled
   differently in every paper, and the corpus held 75 studies and 75 "labs", so
   `by_lab` was numerically identical to `by_study` while four laboratories had
   published eleven papers between them. A laboratory is a person. The failure
   mode of the looser key is a homonym, which `QC_POSSIBLE_LAB_HOMONYM` detects
   from the affiliations the key no longer uses.
3. **Shared cohort** — multiple papers on one animal cohort; flagged manually.
4. **Population x outcome cell** — for generalization claims, evaluate
   leave-one-population-out as well as leave-one-study-out.

Reporting a single cross-validated number without naming the grouping is
considered an invalid result in this project.

### 5b. Perturbation identity (v0.2.0)

Grouping prevents the same study or lab appearing on both sides of a split. It does
not establish that anything on the test side was *learnable* from the train side.

Every one of the 17 molecular target genes in this corpus appears in exactly one
study and one laboratory. Under `by_study` or `by_lab` grouping, a held-out molecular
record's perturbation identity therefore never appears in training. There is nothing
for a model to have learned, and an accuracy figure averaged over such folds measures
the label prior rather than generalisation.

A cross-validated number must therefore be reported alongside the fraction of test
rows whose perturbation identity was unseen in training. `syn4rch.feasibility.identity`
computes it. A high unseen fraction is not a modelling failure to be tuned away; it is
a statement about how this field publishes — one lab, one gene, one paper — and the
only remedy is a second laboratory repeating the same manipulation.

**For circuit records the operative version of this is cell depth, not identity
(D-048, D-057).** `perturbation_key` is `(class, polarity, target)` and a circuit
record has no target, so the whole circuit half holds about three identities and
none is confined to one group — the unseen-identity fraction is 0% and measures
nothing. What distinguishes one circuit manipulation from another is the
population, which lives in the cell. D-057 measured the direction result by how
many laboratories had reached the test row's own cell: at zero the model is
*identical* to the floor by construction, at one the difference spans zero, and
the whole of the established effect is in cells two or more laboratories had
reached. A cross-validated number on circuit records should therefore be read
with that stratification, which `syn4rch.models.generalise` reports.

## 6. What this schema deliberately refuses to represent

- Quantitative dose-response curves (insufficient reporting across the corpus).
- Inferred anatomical wiring presented as measured (belongs to the KG with its own class).
- Any "rejuvenation score". Tier B is an annotation layer, not a training target.
- A `direction` for the assay number where that number runs opposite to the quantity
  `outcome_id` names (D-019). Falling glucose in an insulin tolerance test is
  *rising* insulin sensitivity.
- A condition the paper does not state. `condition_evidence` widens where a
  condition may be *read from*; it does not permit one to be inferred (3.3b).
- A net circuit `valence` for a molecular perturbation. `unclear` is the answer, and
  `perturbation_polarity` carries the sign of what was actually done (3.4b).
- A `challenge` taken from a husbandry sentence. "Mice were allowed ad libitum
  access to water and standard chow" describes how animals were kept, not the
  state any particular measurement was made in, and the same paper measures
  glucose during a glucose load and food intake after a fast. Reading it as a
  challenge would stamp `ad_libitum` across a whole paper and report the result
  as near-complete coverage (D-040).
- A `challenge` implied by the assay that produced the outcome. A glucose
  tolerance test is a glucose load by definition, so deriving `glucose_load`
  from `outcome_id = glucose_tolerance` restates what the outcome already says:
  it discriminates nothing in any key, while raising the apparent condition
  coverage and re-identifying every record it touches. The schema wants the
  acute state of measurements whose assay does *not* fix it -- food intake,
  body weight, blood glucose -- and refuses the rest (D-041).
