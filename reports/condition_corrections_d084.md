# Condition corrections (D-053)

`fill_conditions` fills an empty condition dimension and refuses to overwrite one that is not empty (D-039). That refusal is right, and it left seven records this pass corrects: six storing `fasted` for a fast-refeed intake test, and one whose `challenge` came from a comparator rather than from the animals measured.

## Run

| field | value |
|---|---|
| generated (UTC) | 2026-11-07T05:55:07Z |
| input | `data/curated/records.jsonl` |
| schema version | 0.4.0 |
| records | 709 |
| corrections declared | 62 |
| applied | 55 |
| refused | 7 |
| record ids changed | 55 |
| condition quotes re-checked at their locator | 442 |
| locator failures | 0 |

Outcomes: `applied` 55, `refused_record_not_found` 7.

## What changed

| old id | new id | dimension | from | to | outcome |
|---|---|---|---|---|---|
| `0397b7caf40ad06b` | `-` | challenge | `fasted` | `refed` | refused_record_not_found |
| `087960c4021bf6f2` | `-` | challenge | `fasted` | `refed` | refused_record_not_found |
| `41fca07890c238c8` | `-` | challenge | `fasted` | `refed` | refused_record_not_found |
| `fcea91b0669ad33f` | `-` | challenge | `fasted` | `refed` | refused_record_not_found |
| `93c2e8f0c13a30ed` | `-` | challenge | `fasted` | `refed` | refused_record_not_found |
| `fdb78695221a531e` | `-` | challenge | `fasted` | `refed` | refused_record_not_found |
| `c1b98d981dba36bb` | `-` | challenge | `fasted` | `null` | refused_record_not_found |
| `0597926b69553b2e` | `b3572306da643fbb` | sex | `male` | `null` | applied |
| `207836afe1bbdc2a` | `68c50534d04b4154` | sex | `male` | `null` | applied |
| `53727385b18e0ef6` | `ca0a7f8b8bb43dcf` | sex | `male` | `null` | applied |
| `59ad77826998af41` | `d5e9059338deabec` | sex | `male` | `null` | applied |
| `5ae957cceb49c472` | `3b88335eb733a89c` | sex | `male` | `null` | applied |
| `7edc37f497db2baa` | `2ca0deb2d7a41243` | sex | `male` | `null` | applied |
| `c090086821e1b7e6` | `5315b16c4af5fe8a` | sex | `male` | `null` | applied |
| `d38e944a2ead9f12` | `ad6611a292421be2` | sex | `male` | `null` | applied |
| `e740c67c9602d534` | `2c59d96f59c401e6` | sex | `male` | `null` | applied |
| `18bb2a91c838002a` | `2a8176784670ac66` | sex | `male` | `null` | applied |
| `23b123614a0837b2` | `f4e0a8418300c0e3` | sex | `male` | `null` | applied |
| `2d45f698f9d5b39b` | `c68e9efa5cb6cedd` | sex | `male` | `null` | applied |
| `306eda2851dfc8e7` | `9b44303b48dff2ed` | sex | `male` | `null` | applied |
| `3fb9a23d967289cc` | `1809184b652c944a` | sex | `male` | `null` | applied |
| `47d621b6a5848641` | `2d444ed6957e901f` | sex | `male` | `null` | applied |
| `53c2504388fd8567` | `f9bb373bf8857d17` | sex | `male` | `null` | applied |
| `559b55f5f1c9550d` | `c821e00978b5f8ed` | sex | `male` | `null` | applied |
| `6da604c021094c11` | `3dcf4a10394a4216` | sex | `male` | `null` | applied |
| `7d2be54e91040f3c` | `71540be821b1d3dd` | sex | `male` | `null` | applied |
| `92ec5f0baa4f646a` | `4f29b194d93bee8e` | sex | `male` | `null` | applied |
| `b372933b164eb9f8` | `c5d091c8b1421c67` | sex | `male` | `null` | applied |
| `b42212bc92a1f9a7` | `2018978e6712d946` | sex | `male` | `null` | applied |
| `bbdfc6990941e272` | `87f794a2388a8e38` | sex | `male` | `null` | applied |
| `c3dc855ab1017089` | `16488bddeeb60bf4` | sex | `male` | `null` | applied |
| `d247141bef42765c` | `2c8993f35a1ae97e` | sex | `male` | `null` | applied |
| `d255954c15990b02` | `a69dd1102d12b1df` | sex | `male` | `null` | applied |
| `d427db2d93c11b97` | `7d12bba264d12086` | sex | `male` | `null` | applied |
| `da2e00a4a6298e7f` | `af1d8fa41ac611dd` | sex | `male` | `null` | applied |
| `df66f47098b985df` | `d87060e22af877cf` | sex | `male` | `null` | applied |
| `e17d05d7d557ca5d` | `07ea7fb1070ec3e0` | sex | `male` | `null` | applied |
| `e256eb5de8f78800` | `0ee8b30a2bd190f7` | sex | `male` | `null` | applied |
| `1b30e4473d0d7202` | `6d35d8224751a11e` | sex | `male` | `null` | applied |
| `23fe1e3dc459a0e4` | `929739b50a199b98` | sex | `male` | `null` | applied |
| `71eb042f9f9f2cb7` | `c3a80b16ff6988af` | sex | `male` | `null` | applied |
| `80874b755076808a` | `be9a38dc45a5c38f` | sex | `male` | `null` | applied |
| `9a678ea47cbabe39` | `556a43db01598a8c` | sex | `male` | `null` | applied |
| `d113f0f9d665c1b4` | `f20ed9fc2790008e` | sex | `male` | `null` | applied |
| `f83bc46df17e1ae1` | `21f8ba4ef6384562` | sex | `male` | `null` | applied |
| `2a73f93ea6b13cca` | `7dadc4fd8f0e891a` | sex | `male` | `null` | applied |
| `2df1fee7fd87178a` | `b7957da14cf3768d` | sex | `male` | `null` | applied |
| `33148b7354b90b52` | `40555b61b292d72c` | sex | `male` | `null` | applied |
| `43db6565053e145d` | `8d411fe5fdcde119` | sex | `male` | `null` | applied |
| `5a0778e270e63af0` | `93d72dde7676b3c0` | sex | `male` | `null` | applied |
| `5b48a0d35bc5b378` | `d2376427c0ad4cab` | sex | `male` | `null` | applied |
| `65d21ec042d52cbf` | `094b0f81b04bf1cd` | sex | `male` | `null` | applied |
| `7f997d4c565bb5fe` | `e587d291576af288` | sex | `male` | `null` | applied |
| `b45891b14274922d` | `1669f22136825ec9` | sex | `male` | `null` | applied |
| `cc0375883b50d448` | `ad33fa6476e31b0c` | sex | `male` | `null` | applied |
| `d1414d8588de21e3` | `b3138b7d62b93f41` | sex | `male` | `null` | applied |
| `ddec4c842f1573e9` | `4fe8fa65be5abc6f` | sex | `male` | `null` | applied |
| `f12cdc2916b65836` | `02ae26a5215f1d51` | sex | `male` | `null` | applied |
| `f91902e5fb544751` | `95f4f6c310f1c69f` | sex | `male` | `null` | applied |
| `fb098e47a9945beb` | `8b3deb59ba9256f7` | sex | `male` | `null` | applied |
| `fbf940efe71f13da` | `f9d2df73f54aff86` | sex | `male` | `null` | applied |
| `58eebeb60a276869` | `7e9045d6a02d73fa` | sex | `male` | `null` | applied |

### Why, one by one

* `0397b7caf40ad06b` — food intake is the readout, so the animals had food at the moment of measurement; 'in fasted mice' describes the state they were in before it was returned
* `087960c4021bf6f2` — deprive, return food, measure intake: the measurement is of the refeeding
* `41fca07890c238c8` — fasting-induced hyperphagia is eating, which cannot be measured in an animal that still has no food
* `fcea91b0669ad33f` — deprive, return food, measure intake: the measurement is of the refeeding
* `93c2e8f0c13a30ed` — feeding after a fast is the refed state, not the fasted one
* `fdb78695221a531e` — 'in fasted mice' names the state before food was returned; the readout is the feeding that followed
* `c1b98d981dba36bb` — 'fasted mice' is the comparator the effect is measured against, not the state of the animals measured; the quote supports no challenge for this record
* `0597926b69553b2e` — the Results report the same phenotypes in female Fam172a-loxP and PFKO mice, so sex varies and the Methods purchase sentence does not scope to this record 
* `207836afe1bbdc2a` — the Results report the same phenotypes in female Fam172a-loxP and PFKO mice, so sex varies and the Methods purchase sentence does not scope to this record 
* `53727385b18e0ef6` — the Results report the same phenotypes in female Fam172a-loxP and PFKO mice, so sex varies and the Methods purchase sentence does not scope to this record 
* `59ad77826998af41` — the Results report the same phenotypes in female Fam172a-loxP and PFKO mice, so sex varies and the Methods purchase sentence does not scope to this record 
* `5ae957cceb49c472` — the Results report the same phenotypes in female Fam172a-loxP and PFKO mice, so sex varies and the Methods purchase sentence does not scope to this record 
* `7edc37f497db2baa` — the Results report the same phenotypes in female Fam172a-loxP and PFKO mice, so sex varies and the Methods purchase sentence does not scope to this record 
* `c090086821e1b7e6` — the Results report the same phenotypes in female Fam172a-loxP and PFKO mice, so sex varies and the Methods purchase sentence does not scope to this record 
* `d38e944a2ead9f12` — the Results report the same phenotypes in female Fam172a-loxP and PFKO mice, so sex varies and the Methods purchase sentence does not scope to this record 
* `e740c67c9602d534` — the Results report the same phenotypes in female Fam172a-loxP and PFKO mice, so sex varies and the Methods purchase sentence does not scope to this record 
* `18bb2a91c838002a` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `23b123614a0837b2` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `2d45f698f9d5b39b` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `306eda2851dfc8e7` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `3fb9a23d967289cc` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `47d621b6a5848641` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `53c2504388fd8567` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `559b55f5f1c9550d` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `6da604c021094c11` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `7d2be54e91040f3c` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `92ec5f0baa4f646a` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `b372933b164eb9f8` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `b42212bc92a1f9a7` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `bbdfc6990941e272` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `c3dc855ab1017089` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `d247141bef42765c` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `d255954c15990b02` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `d427db2d93c11b97` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `da2e00a4a6298e7f` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `df66f47098b985df` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `e17d05d7d557ca5d` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `e256eb5de8f78800` — the Methods sentence is an anaesthesia procedure and the paper reports female cohorts in its Results ('improved glucose tolerance in female AGRP-hM3q'), so sex is not constant and this record's own passage does not state it 
* `1b30e4473d0d7202` — the paper carries out metabolic studies on female POMC-G12D mice in its Results, so sex varies; this record's own passage does not state it 
* `23fe1e3dc459a0e4` — the Methods value is inadmissible as for the rest of this paper. NOTE: this record's own passage DOES state its cohort -- 'POMC-G12D mice (males with free access to regular chow)' -- so sex=male is re-sourceable from the passage under 3b. It is removed rather than re-sourced because this module keeps the evidence sentence and changes only the value by design, and a re-source replaces the sentence. Recorded in D-084 as deferred re-curation, not as a lost value 
* `71eb042f9f9f2cb7` — the paper carries out metabolic studies on female POMC-G12D mice in its Results, so sex varies; this record's own passage does not state it 
* `80874b755076808a` — the paper carries out metabolic studies on female POMC-G12D mice in its Results, so sex varies; this record's own passage does not state it 
* `9a678ea47cbabe39` — the paper carries out metabolic studies on female POMC-G12D mice in its Results, so sex varies; this record's own passage does not state it 
* `d113f0f9d665c1b4` — the paper carries out metabolic studies on female POMC-G12D mice in its Results, so sex varies; this record's own passage does not state it 
* `f83bc46df17e1ae1` — the paper carries out metabolic studies on female POMC-G12D mice in its Results, so sex varies; this record's own passage does not state it 
* `2a73f93ea6b13cca` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `2df1fee7fd87178a` — this record's own passage names BOTH sexes ('Male Agrp-Sw ... whereas female Agrp-Sh'), which CURATION_PROTOCOL 3b makes a refusal rather than a pick, so it cannot be re-sourced from the passage either 
* `33148b7354b90b52` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `43db6565053e145d` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `5a0778e270e63af0` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `5b48a0d35bc5b378` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `65d21ec042d52cbf` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `7f997d4c565bb5fe` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `b45891b14274922d` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `cc0375883b50d448` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `d1414d8588de21e3` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `ddec4c842f1573e9` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `f12cdc2916b65836` — same as 2df1fee7fd87178a: the passage names the male and the female arm in one sentence, so no single value is readable from it 
* `f91902e5fb544751` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `fb098e47a9945beb` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `fbf940efe71f13da` — the paper is explicitly sexually dimorphic -- the dataset's own records 1f2c0189ace296a2 and 6153f87b799a9b76 carry sex=female from Results para 8 -- so the Methods purchase sentence cannot speak for every figure 
* `58eebeb60a276869` — the slice recording reported here was made in neurons 'from both male and female mice', so the paper does not hold sex constant 

## The `challenge` distribution, before and after

| value | before | after | delta |
|---|---|---|---|
| `ad_libitum` | 18 | 18 | +0 |
| `cold_exposure` | 4 | 4 | +0 |
| `fasted` | 21 | 21 | +0 |
| `not_reported` | 41 | 41 | +0 |
| `null` | 573 | 573 | +0 |
| `other` | 34 | 34 | +0 |
| `refed` | 14 | 14 | +0 |
| `restraint_stress` | 4 | 4 | +0 |

## What this pass does not do

* It applies no judgement of its own. Every change is a row of `CORRECTIONS` in `syn4rch.curate.correct_conditions`, naming the record, the value that had to be there, the value to put there, the quote fragment that had to match and the decision that authorises it.
* It does not re-read the papers for new conditions. That is `fill_conditions`, which only fills what is empty.
* It does not touch any other field. `challenge` is a `RECORD_ID_FIELDS` component, so the identity moves with it and the id map is how every external reference follows.

