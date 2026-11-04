# SYN4RCH extraction yield assessment (pre-curation)

**Every figure here is an upper bound, not a dataset.** An extraction unit is
a passage in which a neuronal population, a perturbation method and a
physiological readout co-occur. Co-occurrence is not a finding
(`docs/CURATION_PROTOCOL.md` 2). The number of real records is known only
after papers are read and proposals pass the verbatim gate.

## Run

| field | value |
|---|---|
| units file | `data/interim/extraction_units.jsonl` |
| units file written (UTC) | 2026-11-04T00:18:35Z |
| extraction units | 8860 |
| papers yielding >=1 unit | 489 |
| units per paper, median / max | 10 / 168 |
| lab identity available | yes, from the harvest candidates |

## Where the units sit, and why it decides the yield

| section | units | can state the paper's own measured result? |
|---|---|---|
| Discussion | 2943 | almost never: interprets, often other work |
| Results | 2254 | yes, with a citable locator |
| Other | 1362 | unclassified |
| Abstract | 1032 | yes, but with no figure-level locator |
| Introduction | 978 | almost never: cites other work |
| Methods | 291 | yes, for perturbation detail |

Only **2545 of 8860** units (29%) sit in `Results` or `Methods`, the sections where a measured
result with a direction can be quoted against a locator. A further 1032 sit in
the abstract, which states the paper's own claim but gives no figure to cite.
The Introduction and Discussion units are the largest block and the least
usable: the protocol requires the reader to discard a passage that cites or
interprets somebody else's experiment.

## Ceiling on the number of records

| bound | all sections | Results/Methods only |
|---|---|---|
| distinct (population, method, outcome) triples | 1858 | 945 |
| distinct study x triple (**hard ceiling on records**) | 4850 | 1772 |
| papers that could contribute | 489 | 278 |

At most **1772** records can exist from these full texts, and realistically far
fewer: a passage naming all three terms frequently reports no direction for
that particular triple, or reports it for a different population in the same
nucleus. Treat 1772 as the number that cannot be exceeded, not as a forecast.

## Could any prediction task clear the feasibility gate?

`syn4rch.feasibility.gates` calls a task feasible at >=3 independent studies
across >=2 distinct labs. Applying those thresholds to units rather than to
records gives the most optimistic possible answer:

| (population, outcome) cells, Results/Methods units | count | share |
|---|---|---|
| with any extractable unit | 287 | 100% |
| >=3 studies and >=2 labs (gate, upper bound) | 90 | 31% |
| >=2 studies | 135 | 47% |
| exactly 1 study (agreement undefined) | 152 | 53% |

Lab identity resolved for 278 of 278 contributing studies (0 without).

### Cells with the most independent studies

Study and lab counts are nearly equal throughout, which means this corpus is
not one group publishing repeatedly: the independence that SCHEMA.md 5 calls
the project's largest leakage risk is genuinely present at the top of the
distribution.

| population | outcome | studies | labs |
|---|---|---|---|
| `arc_agrp` | `food_intake` | 70 | 70 |
| `arc_pomc` | `food_intake` | 50 | 50 |
| `arc_agrp` | `body_weight` | 46 | 46 |
| `arc_pomc` | `body_weight` | 42 | 42 |
| `lha_orexin` | `sleep` | 26 | 24 |
| `arc_agrp` | `energy_expenditure` | 21 | 21 |
| `arc_agrp` | `locomotor_activity` | 18 | 18 |
| `arc_pomc` | `energy_expenditure` | 18 | 18 |
| `arc_agrp` | `fat_mass` | 17 | 17 |
| `arc_agrp` | `insulin_sensitivity` | 15 | 15 |
| `arc_agrp` | `glucose_tolerance` | 14 | 14 |
| `arc_pomc` | `fat_mass` | 14 | 14 |
| `arc_pomc` | `glucose_tolerance` | 13 | 13 |
| `arc_pomc` | `insulin_sensitivity` | 13 | 13 |
| `arc_agrp` | `blood_glucose` | 12 | 12 |
| `pvh_sim1` | `food_intake` | 12 | 10 |
| `arc_agrp` | `lean_mass` | 11 | 11 |
| `arc_agrp` | `rer` | 11 | 11 |
| `pvh_crh` | `corticosterone` | 11 | 11 |
| `arc_pomc` | `bat_thermogenesis` | 10 | 10 |

### Populations by extractable units

| population | units | papers |
|---|---|---|
| `arc_agrp` | 718 | 76 |
| `arc_pomc` | 490 | 66 |
| `lha_orexin` | 177 | 37 |
| `pvh_crh` | 157 | 26 |
| `vmh_sf1` | 123 | 15 |
| `pvh_sim1` | 94 | 12 |
| `lha_vgat` | 94 | 12 |
| `lha_mch` | 92 | 13 |
| `zi_gaba` | 89 | 9 |
| `lha_vglut2` | 87 | 13 |
| `pvh_oxt` | 84 | 17 |
| `arc_kisspeptin` | 50 | 14 |
| `pvh_avp` | 36 | 5 |
| `pvh_mc4r` | 35 | 5 |
| `pvh_trh` | 35 | 4 |
| `tmn_histamine` | 31 | 8 |
| `dmh_lepr` | 27 | 5 |
| `pvh_nos1` | 25 | 6 |
| `scn_vip` | 20 | 8 |
| `avpv_kisspeptin` | 19 | 6 |

## What this does and does not establish

**Establishes.** The corpus is not uniformly thin. A narrow set of cells,
concentrated in arcuate AgRP and POMC neurons against feeding and body-weight
readouts, has dozens of independent laboratories behind it. If curation
recovers even a modest fraction of those passages as records, direction-of-
effect prediction becomes testable for those cells.

**Does not establish.** That any record exists. That a passage reports a
direction. That independent studies agree. 53% of cells still rest on a single
study, where agreement is undefined rather than perfect, so a model trained
on the dense corner will not generalise to the sparse remainder, and
leave-one-population-out evaluation should be expected to fail.

**Known bias.** Only open-access full texts are here. Paywalled papers
cannot yield records at all, so the dense corner is partly a map of which
laboratories publish open access.
