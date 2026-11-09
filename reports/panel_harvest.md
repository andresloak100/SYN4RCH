# Panel harvest (D-097)

Produced by `python3 analysis/panel_harvest.py`. One Europe PMC query per circuit and readout: the circuit and the readout in the title or abstract, a circuit-activation method anywhere, open access, not a review. `new` counts returned papers with a PMCID that no curated record cites. Nothing here has been read.

## New open-access papers per circuit and readout

| readout | arc_agrp | arc_kisspeptin | arc_pomc | lha_mch | lha_orexin | lha_vgat | poa_galanin | pvh_crh | pvh_oxt | scn_vip | tmn_histamine | vmh_sf1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| food_intake | 74 | 2 | 48 | 9 | 8 | 4 | 0 | 5 | 5 | 0 | 1 | 4 |
| corticosterone | 3 | 0 | 3 | 0 | 1 | 0 | 0 | 16 | 1 | 1 | 0 | 0 |
| wakefulness | 1 | 0 | 0 | 14 | 41 | 2 | 2 | 6 | 7 | 6 | 12 | 0 |
| locomotor_activity | 3 | 1 | 5 | 2 | 17 | 1 | 0 | 0 | 2 | 6 | 2 | 1 |
| body_weight | 34 | 0 | 18 | 5 | 4 | 3 | 0 | 0 | 1 | 0 | 0 | 5 |
| blood_glucose | 8 | 0 | 4 | 0 | 2 | 0 | 0 | 2 | 1 | 0 | 0 | 6 |
| energy_expenditure | 15 | 0 | 10 | 2 | 4 | 0 | 0 | 0 | 3 | 0 | 0 | 5 |
| core_temperature | 6 | 0 | 2 | 0 | 2 | 1 | 1 | 2 | 3 | 1 | 0 | 2 |
| heart_rate | 2 | 0 | 2 | 0 | 3 | 0 | 0 | 2 | 4 | 1 | 0 | 0 |

351 distinct papers returned, 286 not cited by any curated record.

## Every query

### PANEL-arc_agrp-food_intake

hits 102, returned 100, with PMCID 99, new 74

```
(TITLE:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre") OR ABSTRACT:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia") OR ABSTRACT:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_agrp-corticosterone

hits 9, returned 9, with PMCID 9, new 3

```
(TITLE:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre") OR ABSTRACT:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone") OR ABSTRACT:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_agrp-wakefulness

hits 2, returned 2, with PMCID 1, new 1

```
(TITLE:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre") OR ABSTRACT:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep") OR ABSTRACT:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_agrp-locomotor_activity

hits 5, returned 5, with PMCID 5, new 3

```
(TITLE:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre") OR ABSTRACT:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity") OR ABSTRACT:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_agrp-body_weight

hits 41, returned 41, with PMCID 40, new 34

```
(TITLE:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre") OR ABSTRACT:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body weight" OR "body mass" OR "weight gain" OR "weight loss") OR ABSTRACT:("body weight" OR "body mass" OR "weight gain" OR "weight loss")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_agrp-blood_glucose

hits 10, returned 10, with PMCID 10, new 8

```
(TITLE:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre") OR ABSTRACT:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia") OR ABSTRACT:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_agrp-energy_expenditure

hits 20, returned 20, with PMCID 20, new 15

```
(TITLE:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre") OR ABSTRACT:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry") OR ABSTRACT:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_agrp-core_temperature

hits 9, returned 9, with PMCID 9, new 6

```
(TITLE:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre") OR ABSTRACT:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis") OR ABSTRACT:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_agrp-heart_rate

hits 3, returned 3, with PMCID 3, new 2

```
(TITLE:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre") OR ABSTRACT:("ARC AgRP/NPY neurons" OR "AgRP neurons" OR "AgRP/NPY neurons" OR "NPY/AgRP neurons" OR "agouti-related peptide neurons" OR "AgRP-expressing neurons" OR "arcuate AgRP neurons" OR "orexigenic AgRP neurons" OR "Agrp-IRES-Cre" OR "Npy-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular") OR ABSTRACT:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_kisspeptin-food_intake

hits 2, returned 2, with PMCID 2, new 2

```
(TITLE:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre") OR ABSTRACT:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia") OR ABSTRACT:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_kisspeptin-corticosterone

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre") OR ABSTRACT:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone") OR ABSTRACT:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_kisspeptin-wakefulness

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre") OR ABSTRACT:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep") OR ABSTRACT:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_kisspeptin-locomotor_activity

hits 1, returned 1, with PMCID 1, new 1

```
(TITLE:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre") OR ABSTRACT:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity") OR ABSTRACT:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_kisspeptin-body_weight

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre") OR ABSTRACT:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body weight" OR "body mass" OR "weight gain" OR "weight loss") OR ABSTRACT:("body weight" OR "body mass" OR "weight gain" OR "weight loss")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_kisspeptin-blood_glucose

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre") OR ABSTRACT:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia") OR ABSTRACT:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_kisspeptin-energy_expenditure

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre") OR ABSTRACT:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry") OR ABSTRACT:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_kisspeptin-core_temperature

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre") OR ABSTRACT:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis") OR ABSTRACT:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_kisspeptin-heart_rate

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre") OR ABSTRACT:("ARC kisspeptin (KNDy) neurons" OR "KNDy neurons" OR "kisspeptin neurons" OR "arcuate kisspeptin neurons" OR "Kiss1 neurons" OR "kisspeptin/neurokinin B/dynorphin neurons" OR "neurokinin B neurons" OR "NKB neurons" OR "Kiss1-Cre" OR "Tac2-Cre" OR "Pdyn-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular") OR ABSTRACT:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_pomc-food_intake

hits 59, returned 59, with PMCID 58, new 48

```
(TITLE:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2") OR ABSTRACT:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia") OR ABSTRACT:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_pomc-corticosterone

hits 4, returned 4, with PMCID 4, new 3

```
(TITLE:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2") OR ABSTRACT:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone") OR ABSTRACT:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_pomc-wakefulness

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2") OR ABSTRACT:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep") OR ABSTRACT:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_pomc-locomotor_activity

hits 5, returned 5, with PMCID 5, new 5

```
(TITLE:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2") OR ABSTRACT:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity") OR ABSTRACT:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_pomc-body_weight

hits 21, returned 21, with PMCID 20, new 18

```
(TITLE:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2") OR ABSTRACT:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body weight" OR "body mass" OR "weight gain" OR "weight loss") OR ABSTRACT:("body weight" OR "body mass" OR "weight gain" OR "weight loss")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_pomc-blood_glucose

hits 5, returned 5, with PMCID 5, new 4

```
(TITLE:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2") OR ABSTRACT:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia") OR ABSTRACT:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_pomc-energy_expenditure

hits 13, returned 13, with PMCID 13, new 10

```
(TITLE:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2") OR ABSTRACT:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry") OR ABSTRACT:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_pomc-core_temperature

hits 4, returned 4, with PMCID 4, new 2

```
(TITLE:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2") OR ABSTRACT:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis") OR ABSTRACT:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-arc_pomc-heart_rate

hits 3, returned 3, with PMCID 2, new 2

```
(TITLE:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2") OR ABSTRACT:("ARC POMC neurons" OR "POMC neurons" OR "proopiomelanocortin neurons" OR "POMC-expressing neurons" OR "arcuate POMC neurons" OR "POMC/CART neurons" OR "anorexigenic POMC neurons" OR "alpha-MSH neurons" OR "Pomc-Cre" OR "Pomc-CreERT2")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular") OR ABSTRACT:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_mch-food_intake

hits 9, returned 9, with PMCID 9, new 9

```
(TITLE:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre") OR ABSTRACT:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia") OR ABSTRACT:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_mch-corticosterone

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre") OR ABSTRACT:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone") OR ABSTRACT:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_mch-wakefulness

hits 18, returned 18, with PMCID 18, new 14

```
(TITLE:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre") OR ABSTRACT:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep") OR ABSTRACT:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_mch-locomotor_activity

hits 2, returned 2, with PMCID 2, new 2

```
(TITLE:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre") OR ABSTRACT:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity") OR ABSTRACT:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_mch-body_weight

hits 5, returned 5, with PMCID 5, new 5

```
(TITLE:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre") OR ABSTRACT:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body weight" OR "body mass" OR "weight gain" OR "weight loss") OR ABSTRACT:("body weight" OR "body mass" OR "weight gain" OR "weight loss")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_mch-blood_glucose

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre") OR ABSTRACT:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia") OR ABSTRACT:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_mch-energy_expenditure

hits 2, returned 2, with PMCID 2, new 2

```
(TITLE:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre") OR ABSTRACT:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry") OR ABSTRACT:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_mch-core_temperature

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre") OR ABSTRACT:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis") OR ABSTRACT:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_mch-heart_rate

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre") OR ABSTRACT:("LHA MCH neurons" OR "MCH neurons" OR "melanin-concentrating hormone neurons" OR "Pmch neurons" OR "MCH-expressing neurons" OR "lateral hypothalamic MCH neurons" OR "Pmch-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular") OR ABSTRACT:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_orexin-food_intake

hits 9, returned 9, with PMCID 9, new 8

```
(TITLE:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre") OR ABSTRACT:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia") OR ABSTRACT:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_orexin-corticosterone

hits 1, returned 1, with PMCID 1, new 1

```
(TITLE:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre") OR ABSTRACT:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone") OR ABSTRACT:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_orexin-wakefulness

hits 56, returned 56, with PMCID 56, new 41

```
(TITLE:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre") OR ABSTRACT:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep") OR ABSTRACT:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_orexin-locomotor_activity

hits 21, returned 21, with PMCID 21, new 17

```
(TITLE:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre") OR ABSTRACT:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity") OR ABSTRACT:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_orexin-body_weight

hits 5, returned 5, with PMCID 5, new 4

```
(TITLE:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre") OR ABSTRACT:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body weight" OR "body mass" OR "weight gain" OR "weight loss") OR ABSTRACT:("body weight" OR "body mass" OR "weight gain" OR "weight loss")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_orexin-blood_glucose

hits 3, returned 3, with PMCID 2, new 2

```
(TITLE:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre") OR ABSTRACT:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia") OR ABSTRACT:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_orexin-energy_expenditure

hits 5, returned 5, with PMCID 5, new 4

```
(TITLE:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre") OR ABSTRACT:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry") OR ABSTRACT:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_orexin-core_temperature

hits 2, returned 2, with PMCID 2, new 2

```
(TITLE:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre") OR ABSTRACT:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis") OR ABSTRACT:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_orexin-heart_rate

hits 3, returned 3, with PMCID 3, new 3

```
(TITLE:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre") OR ABSTRACT:("LHA orexin/hypocretin neurons" OR "orexin neurons" OR "hypocretin neurons" OR "orexinergic neurons" OR "orexin/hypocretin neurons" OR "Hcrt neurons" OR "lateral hypothalamic orexin neurons" OR "Hcrt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular") OR ABSTRACT:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_vgat-food_intake

hits 4, returned 4, with PMCID 4, new 4

```
(TITLE:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre") OR ABSTRACT:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia") OR ABSTRACT:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_vgat-corticosterone

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre") OR ABSTRACT:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone") OR ABSTRACT:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_vgat-wakefulness

hits 2, returned 2, with PMCID 2, new 2

```
(TITLE:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre") OR ABSTRACT:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep") OR ABSTRACT:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_vgat-locomotor_activity

hits 1, returned 1, with PMCID 1, new 1

```
(TITLE:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre") OR ABSTRACT:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity") OR ABSTRACT:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_vgat-body_weight

hits 3, returned 3, with PMCID 3, new 3

```
(TITLE:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre") OR ABSTRACT:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body weight" OR "body mass" OR "weight gain" OR "weight loss") OR ABSTRACT:("body weight" OR "body mass" OR "weight gain" OR "weight loss")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_vgat-blood_glucose

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre") OR ABSTRACT:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia") OR ABSTRACT:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_vgat-energy_expenditure

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre") OR ABSTRACT:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry") OR ABSTRACT:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_vgat-core_temperature

hits 1, returned 1, with PMCID 1, new 1

```
(TITLE:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre") OR ABSTRACT:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis") OR ABSTRACT:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-lha_vgat-heart_rate

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre") OR ABSTRACT:("LHA GABAergic (Vgat) neurons" OR "LH GABA neurons" OR "LHA GABAergic neurons" OR "LH Vgat neurons" OR "lateral hypothalamic GABAergic neurons" OR "LHA Vgat neurons" OR "GABAergic lateral hypothalamic neurons" OR "Vgat-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular") OR ABSTRACT:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-poa_galanin-food_intake

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre") OR ABSTRACT:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia") OR ABSTRACT:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-poa_galanin-corticosterone

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre") OR ABSTRACT:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone") OR ABSTRACT:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-poa_galanin-wakefulness

hits 4, returned 4, with PMCID 3, new 2

```
(TITLE:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre") OR ABSTRACT:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep") OR ABSTRACT:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-poa_galanin-locomotor_activity

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre") OR ABSTRACT:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity") OR ABSTRACT:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-poa_galanin-body_weight

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre") OR ABSTRACT:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body weight" OR "body mass" OR "weight gain" OR "weight loss") OR ABSTRACT:("body weight" OR "body mass" OR "weight gain" OR "weight loss")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-poa_galanin-blood_glucose

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre") OR ABSTRACT:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia") OR ABSTRACT:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-poa_galanin-energy_expenditure

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre") OR ABSTRACT:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry") OR ABSTRACT:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-poa_galanin-core_temperature

hits 2, returned 2, with PMCID 2, new 1

```
(TITLE:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre") OR ABSTRACT:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis") OR ABSTRACT:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-poa_galanin-heart_rate

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre") OR ABSTRACT:("Preoptic galanin neurons" OR "galanin neurons" OR "POA galanin neurons" OR "MPOA galanin neurons" OR "galaninergic preoptic neurons" OR "VLPO galanin neurons" OR "Gal-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular") OR ABSTRACT:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_crh-food_intake

hits 5, returned 5, with PMCID 5, new 5

```
(TITLE:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre") OR ABSTRACT:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia") OR ABSTRACT:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_crh-corticosterone

hits 18, returned 18, with PMCID 18, new 16

```
(TITLE:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre") OR ABSTRACT:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone") OR ABSTRACT:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_crh-wakefulness

hits 10, returned 10, with PMCID 10, new 6

```
(TITLE:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre") OR ABSTRACT:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep") OR ABSTRACT:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_crh-locomotor_activity

hits 2, returned 2, with PMCID 2, new 0

```
(TITLE:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre") OR ABSTRACT:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity") OR ABSTRACT:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_crh-body_weight

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre") OR ABSTRACT:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body weight" OR "body mass" OR "weight gain" OR "weight loss") OR ABSTRACT:("body weight" OR "body mass" OR "weight gain" OR "weight loss")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_crh-blood_glucose

hits 2, returned 2, with PMCID 2, new 2

```
(TITLE:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre") OR ABSTRACT:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia") OR ABSTRACT:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_crh-energy_expenditure

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre") OR ABSTRACT:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry") OR ABSTRACT:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_crh-core_temperature

hits 2, returned 2, with PMCID 2, new 2

```
(TITLE:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre") OR ABSTRACT:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis") OR ABSTRACT:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_crh-heart_rate

hits 2, returned 2, with PMCID 2, new 2

```
(TITLE:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre") OR ABSTRACT:("PVH CRH neurons" OR "CRH neurons" OR "CRF neurons" OR "corticotropin-releasing hormone neurons" OR "corticotropin-releasing factor neurons" OR "PVN CRH neurons" OR "PVH CRF neurons" OR "CRH-expressing PVN neurons" OR "Crh-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular") OR ABSTRACT:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_oxt-food_intake

hits 8, returned 8, with PMCID 7, new 5

```
(TITLE:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre") OR ABSTRACT:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia") OR ABSTRACT:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_oxt-corticosterone

hits 1, returned 1, with PMCID 1, new 1

```
(TITLE:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre") OR ABSTRACT:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone") OR ABSTRACT:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_oxt-wakefulness

hits 7, returned 7, with PMCID 7, new 7

```
(TITLE:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre") OR ABSTRACT:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep") OR ABSTRACT:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_oxt-locomotor_activity

hits 2, returned 2, with PMCID 2, new 2

```
(TITLE:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre") OR ABSTRACT:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity") OR ABSTRACT:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_oxt-body_weight

hits 1, returned 1, with PMCID 1, new 1

```
(TITLE:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre") OR ABSTRACT:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body weight" OR "body mass" OR "weight gain" OR "weight loss") OR ABSTRACT:("body weight" OR "body mass" OR "weight gain" OR "weight loss")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_oxt-blood_glucose

hits 1, returned 1, with PMCID 1, new 1

```
(TITLE:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre") OR ABSTRACT:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia") OR ABSTRACT:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_oxt-energy_expenditure

hits 3, returned 3, with PMCID 3, new 3

```
(TITLE:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre") OR ABSTRACT:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry") OR ABSTRACT:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_oxt-core_temperature

hits 3, returned 3, with PMCID 3, new 3

```
(TITLE:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre") OR ABSTRACT:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis") OR ABSTRACT:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-pvh_oxt-heart_rate

hits 5, returned 5, with PMCID 5, new 4

```
(TITLE:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre") OR ABSTRACT:("PVH oxytocin neurons" OR "oxytocin neurons" OR "OXT neurons" OR "PVN oxytocin neurons" OR "oxytocinergic neurons" OR "Oxt-expressing neurons" OR "magnocellular oxytocin neurons" OR "Oxt-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular") OR ABSTRACT:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-scn_vip-food_intake

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre") OR ABSTRACT:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia") OR ABSTRACT:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-scn_vip-corticosterone

hits 1, returned 1, with PMCID 1, new 1

```
(TITLE:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre") OR ABSTRACT:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone") OR ABSTRACT:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-scn_vip-wakefulness

hits 6, returned 6, with PMCID 6, new 6

```
(TITLE:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre") OR ABSTRACT:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep") OR ABSTRACT:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-scn_vip-locomotor_activity

hits 7, returned 7, with PMCID 7, new 6

```
(TITLE:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre") OR ABSTRACT:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity") OR ABSTRACT:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-scn_vip-body_weight

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre") OR ABSTRACT:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body weight" OR "body mass" OR "weight gain" OR "weight loss") OR ABSTRACT:("body weight" OR "body mass" OR "weight gain" OR "weight loss")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-scn_vip-blood_glucose

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre") OR ABSTRACT:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia") OR ABSTRACT:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-scn_vip-energy_expenditure

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre") OR ABSTRACT:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry") OR ABSTRACT:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-scn_vip-core_temperature

hits 1, returned 1, with PMCID 1, new 1

```
(TITLE:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre") OR ABSTRACT:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis") OR ABSTRACT:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-scn_vip-heart_rate

hits 1, returned 1, with PMCID 1, new 1

```
(TITLE:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre") OR ABSTRACT:("SCN VIP neurons" OR "VIP neurons" OR "vasoactive intestinal peptide neurons" OR "VIPergic SCN neurons" OR "SCN core VIP neurons" OR "Vip-expressing SCN neurons" OR "Vip-IRES-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular") OR ABSTRACT:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-tmn_histamine-food_intake

hits 2, returned 2, with PMCID 1, new 1

```
(TITLE:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre") OR ABSTRACT:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia") OR ABSTRACT:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-tmn_histamine-corticosterone

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre") OR ABSTRACT:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone") OR ABSTRACT:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-tmn_histamine-wakefulness

hits 17, returned 17, with PMCID 16, new 12

```
(TITLE:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre") OR ABSTRACT:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep") OR ABSTRACT:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-tmn_histamine-locomotor_activity

hits 5, returned 5, with PMCID 4, new 2

```
(TITLE:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre") OR ABSTRACT:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity") OR ABSTRACT:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-tmn_histamine-body_weight

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre") OR ABSTRACT:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body weight" OR "body mass" OR "weight gain" OR "weight loss") OR ABSTRACT:("body weight" OR "body mass" OR "weight gain" OR "weight loss")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-tmn_histamine-blood_glucose

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre") OR ABSTRACT:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia") OR ABSTRACT:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-tmn_histamine-energy_expenditure

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre") OR ABSTRACT:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry") OR ABSTRACT:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-tmn_histamine-core_temperature

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre") OR ABSTRACT:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis") OR ABSTRACT:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-tmn_histamine-heart_rate

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre") OR ABSTRACT:("TMN histaminergic neurons" OR "histaminergic neurons" OR "TMN histamine neurons" OR "tuberomammillary histamine neurons" OR "histamine neurons" OR "HDC neurons" OR "histidine decarboxylase neurons" OR "TMN neurons" OR "Hdc-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular") OR ABSTRACT:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-vmh_sf1-food_intake

hits 6, returned 6, with PMCID 6, new 4

```
(TITLE:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre") OR ABSTRACT:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia") OR ABSTRACT:("food intake" OR "feeding" OR "food consumption" OR "caloric intake" OR "hyperphagia" OR "hypophagia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-vmh_sf1-corticosterone

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre") OR ABSTRACT:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone") OR ABSTRACT:("corticosterone" OR "HPA axis" OR "glucocorticoid" OR "adrenocorticotropic hormone")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-vmh_sf1-wakefulness

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre") OR ABSTRACT:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep") OR ABSTRACT:("wakefulness" OR "arousal" OR "sleep" OR "NREM" OR "REM sleep")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-vmh_sf1-locomotor_activity

hits 2, returned 2, with PMCID 2, new 1

```
(TITLE:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre") OR ABSTRACT:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity") OR ABSTRACT:("locomotor activity" OR "locomotion" OR "physical activity" OR "spontaneous activity" OR "wheel running" OR "home cage activity")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-vmh_sf1-body_weight

hits 5, returned 5, with PMCID 5, new 5

```
(TITLE:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre") OR ABSTRACT:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body weight" OR "body mass" OR "weight gain" OR "weight loss") OR ABSTRACT:("body weight" OR "body mass" OR "weight gain" OR "weight loss")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-vmh_sf1-blood_glucose

hits 6, returned 6, with PMCID 6, new 6

```
(TITLE:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre") OR ABSTRACT:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia") OR ABSTRACT:("blood glucose" OR "glycemia" OR "glycaemia" OR "plasma glucose" OR "hyperglycemia" OR "hypoglycemia")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-vmh_sf1-energy_expenditure

hits 6, returned 6, with PMCID 6, new 5

```
(TITLE:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre") OR ABSTRACT:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry") OR ABSTRACT:("energy expenditure" OR "oxygen consumption" OR "metabolic rate" OR "indirect calorimetry")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-vmh_sf1-core_temperature

hits 2, returned 2, with PMCID 2, new 2

```
(TITLE:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre") OR ABSTRACT:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis") OR ABSTRACT:("body temperature" OR "core temperature" OR "hypothermia" OR "hyperthermia" OR "thermogenesis")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

### PANEL-vmh_sf1-heart_rate

hits 0, returned 0, with PMCID 0, new 0

```
(TITLE:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre") OR ABSTRACT:("VMH SF1 (Nr5a1) neurons" OR "SF1 neurons" OR "SF-1 neurons" OR "steroidogenic factor 1 neurons" OR "VMH SF1 neurons" OR "Nr5a1 neurons" OR "SF1-expressing VMH neurons" OR "Sf1-Cre")) AND ("chemogenetic" OR "chemogenetically" OR "DREADD" OR "DREADDs" OR "hM3Dq" OR "Gq-DREADD" OR "clozapine-N-oxide" OR "optogenetic" OR "optogenetically" OR "channelrhodopsin" OR "ChR2" OR "photostimulation" OR "Chrimson" OR "optogenetic activation" OR "chemogenetic activation") AND (TITLE:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular") OR ABSTRACT:("heart rate" OR "bradycardia" OR "tachycardia" OR "blood pressure" OR "cardiovascular")) AND OPEN_ACCESS:y AND NOT PUB_TYPE:"review"
```

