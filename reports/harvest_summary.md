# SYN4RCH harvest summary

Every number below is computed from `data/interim/candidates.jsonl` and
`docs/QUERIES.md`, both written by this run. Nothing is estimated.

**Caveat that governs the whole table.** This is a title/abstract screen. It
cannot tell you how many *extractable* (population x perturbation x outcome)
records exist -- only full-text curation can. These are candidate papers, an
UPPER BOUND on usable records' source papers, and the `include` count in
particular will shrink once full texts are read.

## Run

| field | value |
|---|---|
| started (UTC) | 2026-11-03T23:33:26Z |
| finished (UTC) | 2026-11-03T23:55:48Z |
| queries issued | 404 |
| queries ok / errored | 404 / 0 |
| queries truncated by retmax (250 / broad 1000) | 22 |
| unique papers after dedupe | 9653 |
| duplicate merge events | 7803 |
| cached HTTP responses in sqlite | 1667 |
| cache hits / misses this run | 693 / 632 |

## Screening

| decision | papers |
|---|---|
| include | 5288 |
| uncertain | 1567 |
| exclude | 2798 |
| **total** | **9653** |

### Exclusion reasons (top)

| reason | papers |
|---|---|
| `publication_type` | 1579 |
| `missing_modality_site_readout` | 391 |
| `missing_site_readout` | 385 |
| `missing_modality_site` | 233 |
| `missing_modality_readout` | 210 |

### Flags (recorded, never used to exclude)

| flag | papers |
|---|---|
| `abstract_missing` | 206 |
| `human_observational` | 35 |
| `possibly_in_vitro_only` | 4 |

## Physiological readout mentions

| field | value |
|---|---|
| candidate papers whose title/abstract mentions >=1 vocabulary readout | 7692 |
| papers with an abstract at all | 9447 |

Per outcome id (papers mentioning it; a paper can mention many):

| outcome_id | papers |
|---|---|
| `blood_glucose` | 450 |
| `glucose_tolerance` | 196 |
| `insulin_sensitivity` | 328 |
| `plasma_insulin` | 231 |
| `plasma_glucagon` | 31 |
| `hepatic_glucose_production` | 66 |
| `corticosterone` | 787 |
| `gh_secretion` | 335 |
| `igf1` | 62 |
| `thyroid_hormone` | 294 |
| `gonadotropins` | 943 |
| `food_intake` | 2224 |
| `energy_expenditure` | 606 |
| `rer` | 60 |
| `body_weight` | 1805 |
| `fat_mass` | 362 |
| `lean_mass` | 43 |
| `bat_thermogenesis` | 184 |
| `core_temperature` | 794 |
| `heart_rate` | 672 |
| `blood_pressure` | 1047 |
| `sympathetic_nerve_activity` | 395 |
| `parasympathetic_tone` | 30 |
| `wat_lipolysis` | 53 |
| `hepatic_lipid` | 59 |
| `bone_density` | 44 |
| `il6` | 147 |
| `tnf` | 183 |
| `inflammatory_markers_other` | 409 |
| `immune_cell_composition` | 31 |
| `gut_motility` | 110 |
| `sleep` | 1210 |
| `locomotor_activity` | 686 |
| `water_intake` | 456 |
| `osmolality_avp` | 143 |
| `estrous_cyclicity` | 310 |
| `median_lifespan` | 34 |
| `maximum_lifespan` | 0 |
| `frailty_index` | 1 |
| `grip_strength` | 7 |
| `cognitive_decline` | 260 |
| `tissue_aging_markers` | 1 |
| `senescence_markers` | 29 |

## Unique papers per population vocabulary entry

`matched_in_text` = the screener found one of that entry's terms in the
title/abstract/MeSH. `found_by_its_queries` = the paper came back from a query
built from that entry. Zeros are kept, and they are the point of the table.

| population_id | curation_confidence | matched_in_text | of which `include` | found_by_its_queries |
|---|---|---|---|---|
| `arc_agrp` | high | 321 | 208 | 445 |
| `lha_orexin` | high | 257 | 161 | 398 |
| `arc_pomc` | high | 225 | 148 | 392 |
| `arc_kisspeptin` | medium | 216 | 80 | 278 |
| `pvh_crh` | high | 171 | 91 | 313 |
| `pvh_oxt` | high | 161 | 61 | 273 |
| `lha_mch` | medium | 102 | 57 | 270 |
| `pvh_avp` | high | 89 | 40 | 284 |
| `tmn_histamine` | high | 82 | 45 | 271 |
| `scn_vip` | high | 58 | 17 | 252 |
| `vmh_sf1` | high | 46 | 23 | 251 |
| `pvh_trh` | medium | 35 | 19 | 218 |
| `arc_ghrh` | medium | 29 | 17 | 110 |
| `avpv_kisspeptin` | medium | 28 | 15 | 223 |
| `pvh_mc4r` | medium | 26 | 6 | 194 |
| `arc_th` | medium | 25 | 13 | 219 |
| `lha_vgat` | high | 25 | 16 | 256 |
| `zi_gaba` | high | 22 | 8 | 250 |
| `pvh_sim1` | medium | 19 | 8 | 249 |
| `pvh_nos1` | low | 18 | 6 | 253 |
| `poa_warm_sensitive` | low | 18 | 10 | 167 |
| `dmh_lepr` | medium | 14 | 7 | 168 |
| `poa_galanin` | medium | 14 | 5 | 169 |
| `lha_vglut2` | high | 12 | 7 | 248 |
| `arc_sst` | low | 9 | 3 | 242 |
| `vmh_bdnf` | medium | 4 | 1 | 95 |
| `scn_avp` | high | 4 | 2 | 59 |
| `mbh_ikkb` | low | 2 | 0 | 2 |
| `arc_tbx3` | low | 0 | 0 | 2 |
| `dmh_bdnf` | low | 0 | 0 | 1 |

Populations with **zero** text matches in the whole corpus: `arc_tbx3`, `dmh_bdnf`

## Full text availability

| field | papers |
|---|---|
| open access per Europe PMC | 3480 |
| not open access per Europe PMC | 6173 |
| OA status unknown (lookup budget not reached) | 0 |
| has a PMCID | 5527 |
| OA full text XML actually downloaded to data/raw/fulltext/ | 974 |
| full-text retrieval attempted / succeeded | 992 / 974 |
| abstract only (no full text retrieved) | 8679 |

Open-access lookups: 5258 papers in 106 Europe PMC batch queries, 5238 resolved.

## Bibliographic completeness (SCHEMA.md 3.1)

| field present | papers |
|---|---|
| DOI | 9264 |
| title+journal+year+first_author+last_author | 9463 |
| last_author_affiliation (needed for lab-level leakage control) | 5199 |
| Crossref backfills applied | 65 |

## Known weaknesses of this search

* retmax (250, broad 1000) truncated 22 of 404 successful queries, so the corpus is a
  relevance-ranked sample of the matrix, not its closure.
* Phrase matching on title/abstract misses papers that name a population only
  in the methods section, which is common for Cre-line targeting.
* The screener's `site` criterion is satisfied by a region term alone, so
  region-level manipulations that cannot yield a population-specific record
  still screen in.
* Ambiguous acronyms are excluded from matching, which costs recall for papers
  that only ever write `PVN` or `ARC` -- but the alternative is false matches on
  English words.
* **Europe PMC (the `E1`/`E3` queries) searches FULL TEXT, not just the
  abstract**, so those queries are far less precise than the PubMed `[tiab]`
  queries: a paper that mentions a population once in its methods matches. This
  is why `found_by_its_queries` is much larger than `matched_in_text` for most
  populations. `matched_in_text` is the number to trust; `found_by_its_queries`
  is only a record of which query surfaced the paper.
* PubMed silently splits a quoted phrase it does not hold in its phrase index
  into separate terms, which loosens some blocks in ways the log cannot show.
* The `include` bucket is inflated by papers where the modality, the site and
  the readout are all mentioned but are not connected to each other -- e.g. a
  pharmacology paper citing optogenetic work in its introduction.
* No non-English literature, no books/theses, no conference abstracts.

