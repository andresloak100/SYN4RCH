# Independent-animal cohort for the real 3D axon benchmark: the locked plan

**Locked 2026-11-10 by the commit that adds this file.** Made from metadata only: no trace file and no regional
outcome table was opened (the script blocks both), and no result of the other agent's experiments (`brm_*`,
`d110_*`, `d111_*`) was read before the lock. The one prior awareness is the sentence recorded in
`reports/d109_jiao_audit.md` (BRM "had lower held-out error than copying rules and simple averages").

Everything below is reproduced by `python3 -I analysis/axon_dataset_audit.py`. It writes
`reports/axon_dataset_manifest.tsv` (one row per neuron, SHA-256 `649f829a592cd3d13bcda5cffdbd049e759568d38bacb144f6a61220610a1ca4`)
and `reports/axon_dataset_animals.tsv` (one row per animal, SHA-256 `11a0680bfcb53b204e7224a1f982174bbf11b25fa6aedb3e9476ed465f69df9f`);
those two tables are not version controlled (repository rule for `reports/`), so the animal table is also printed in the
appendix here. Source notes with quotations and addresses: `reports/axon_dataset_sources.md` (more Jiao traces) and
`reports/axon_dataset_other_sources.md` (other datasets). Each was written by one agent and has **not** been
independently reviewed; this plan relies on them only for access, licence and archive sizes, and re-checked their
per-animal counts against the index (166 of 166 animals agree).

## 1. Source, units and identifiers

* **Data:** Jiao et al. 2025 (Nat Neurosci, doi:10.1038/s41593-025-01919-0): 7,180 reconstructed hypothalamic neurons
  from 166 mice.
* **Neuron ID** `<animal>_<n>` (for example `201232_056`); **animal** = the ID prefix (the portal's own sample ID).
* **Source label:** the animal's Cre line (or the orexin promoter virus). One label per animal, 16 labels. It is a
  genetic label, not measured gene expression.
* **Soma label:** the fine label of the paper's soma table (59 labels; 1,316 neurons are only `HY`, 136 are on
  fibre-tract labels). The portal index uses 62 acronyms and differs on 42 neurons; it is kept in the local manifest only.

## 2. Conventions (from sources; details in `reports/swc_access.md`, `reports/d109_jiao_conventions.md`)

* **Trace format:** SWC, columns `id type AP DV ML radius parent`, micrometres, Allen CCFv3 (2017) on the 10 um grid;
  +AP posterior, +DV ventral, +ML toward the animal's right; midline ML 5,700.
* **Hemisphere:** portal and data-centre files keep the true hemisphere (index: Left 3,071, Right 4,109). The paper's
  tables and the Zenodo files are mirrored to the left by ML' = 11,400 - ML (whole neuron).
* **Benchmark frame:** raw files are stored unchanged. Modelling and scoring use the soma-ipsilateral frame, i.e. a
  right-hemisphere neuron is mirrored by the rule above, which keeps ipsilateral and contralateral relations.
* **Node types:** 1 soma, 2 axon, 3 dendrite, 0 undefined; 10, 11, 12 are undocumented single marker nodes. A file
  codes its axon as type 2 **or** as type 0, never both; which one is known only after the file is read. Of the 37
  traces on disk, 16 use type 2 (exactly the 16 development neurons) and 21 use type 0; no animal has both.
* **Axon rule for the benchmark (the authors' own, from their length code):** axon = all cable that is not dendrite
  (types 3, 4) and not the soma; a marker node takes the type of its nearest non-marker ancestor. The file's coding
  (2 or 0) is recorded at ingestion.
* **Reconstruction type (index):** `Axon_and_dendrite` 6,871, `Axon_only` 309.
* **Registration:** fMOST, Fast Neurite Tracer with human checks, landmark-based ANTs registration to CCFv3. The
  accuracy figure was not read. Three of 37 inspected somas lie 9 to 33 um outside their labelled nucleus mesh.

## 3. Classes and exact counts

| class | neurons | animals |
|---|---:|---:|
| development (the other agent's 16 neurons; stay development data) | 16 | 14 |
| inspected before this plan, not development (viewed and validated, never fitted or scored) | 21 | 18 |
| new neurons of the 32 already-inspected animals | 1,575 | (same 32) |
| neurons of genuinely new animals (no trace of the animal ever inspected) | 5,568 | 134 |
| **all** | **7,180** | **166** |

The 16 development IDs are those recorded in `reports/d109_protocol.md` section 16; the other agent's own manifest
has not been read and should confirm them. Two animals are also in Li et al. 2024: 195037 (new, 58 neurons) and
195038 (inspected, 48 neurons).

## 4. The split: whole animals, fixed by a hash

**Rule.**

1. Never held out (training only): the 32 inspected animals (14 of them development animals), the Li-shared new
   animal 195037, and the 10 new animals with fewer than 5 neurons (26 neurons).
2. The other 123 new animals are **eligible**. Inside each source label they are ordered by
   SHA-256(`SYN4RCH-axon-cohort-2026-11-10:` + animal ID). With m eligible animals in a label: m >= 10 gives
   round(0.30 m) test and round(0.20 m) validation animals; 4 <= m <= 9 the same with at least one of each; m = 2 or 3
   gives one test and one validation; m = 1 gives one test. The first animals in hash order are test, the next
   validation, the rest training.
3. The salt was written before the split was first computed, the script was run once, and the result is not to be
   re-drawn.

**Result.**

| part | animals | neurons | Left / Right | Axon_only | soma labels (only `HY` or a tract) | data-centre archives (exact bytes) |
|---|---:|---:|---|---:|---|---:|
| training | 101 | 4,448 | 1,871 / 2,577 | 174 | 57 (900) | 4,037,546,018 |
| validation | 25 | 1,185 | 464 / 721 | 24 | 54 (237) | 1,201,150,823 |
| **test** | **40** | **1,547** | 736 / 811 | 111 | 46 (315) | 1,421,163,124 |
| all | 166 | 7,180 | 3,071 / 4,109 | 309 | 59 (1,452) | 6,659,859,965 |

Training pool: 14 development animals (887 neurons), 18 inspected animals (725), the Li-shared new animal (58), 10
small new animals (26), and 58 new animals placed by hash (2,752).

Animals / neurons per source label:

| label | training | validation | test |
|---|---|---|---|
| Adcyap1 | 12 / 709 | 4 / 140 | 6 / 199 |
| Sst | 13 / 813 | 4 / 188 | 5 / 194 |
| Penk | 10 / 304 | 3 / 157 | 5 / 188 |
| Nts | 9 / 542 | 3 / 265 | 4 / 232 |
| Pdyn | 11 / 445 | 2 / 123 | 4 / 287 |
| Tac2 | 10 / 371 | 3 / 131 | 4 / 161 |
| Oxt | 8 / 290 | 1 / 9 | 2 / 49 |
| Pomc | 6 / 238 | 1 / 20 | 2 / 46 |
| Tac1 | 6 / 386 | 1 / 76 | 2 / 64 |
| Orexin | 4 / 69 | 1 / 19 | 1 / 15 |
| Agrp | 3 / 87 | 1 / 10 | 1 / 5 |
| Crh | 1 / 16 | 1 / 47 | 1 / 13 |
| Pmch | 3 / 107 | 0 | 1 / 9 |
| Avp | 2 / 13 | 0 | 1 / 28 |
| Trh | 1 / 42 | 0 | 1 / 57 |
| Vip | 2 / 16 | 0 | 0 |

**Donor coverage (inputs only).** Distance from a held-out soma to the nearest training soma, mirrored frame: test
median 56 um, 90th percentile 101 um, maximum 222 um; with the same source label median 112 um, 90th percentile
251 um, maximum 619 um. Validation: 59 / 103 / 531 um, and 117 / 279 / 1,373 um with the same label. So, once the
training traces are available, nearly every held-out neuron has a measured neighbour from another animal within
about 0.1 mm.

**What the split cannot guarantee.** Source labels are clustered in animal-ID order (batches). Ordering by hash mixes
a label's animals over the three parts, but no batch record exists, so a held-out animal may share a batch with a
training animal.

## 5. Overlap exclusions

* No validation or test animal has any previously inspected trace; all 37 traces on disk belong to training animals.
* The 16 development neurons and their 14 animals are training only.
* Both Li-shared animals are training only; if Li et al. 2024 is ever added, brains 195037 and 195038 are dropped
  from it whole.
* MouseLight neurons are another laboratory's animals: no overlap with Jiao animals.

## 6. What is on disk, and what the benchmark needs

**On disk now:** 37 traces, all of training animals: the 16 development neurons (axon coded type 2) and 21 others
(axon coded type 0). **No validation or test trace is on disk.** The regional table of the paper (axon length in 49
regions for all 7,180 neurons) is on disk; rows of validation and test animals are sealed (section 7).

**Needed:** the traces of all three parts. Exact archive bytes are in the table of section 4; every animal's archive
size is in the appendix.

| stage | what | data-centre archives | as single portal files |
|---|---|---|---|
| A | training + validation (develop, choose and freeze a method) | 126 zips, 5,238,696,841 bytes | 5,596 files not on disk, about 6 to 8 GB |
| B | test, fetched only after the freeze | 40 zips, 1,421,163,124 bytes | 1,547 files, about 1.5 to 2.2 GB |

Fetching the test traces only after the method is frozen is the strongest seal and costs nothing.

**Routes** (`reports/axon_dataset_sources.md`, A2 to A7):

* **Data centre (BSDC, doi:10.12412/BSDC.1667278800.20003):** one zip per animal, licence CC BY-NC-ND 4.0. Its
  download buttons need a logged-in account (registration with e-mail confirmation, picture code, terms box).
  Only the user can create the account and accept the terms.
* **Portal, per-neuron address:** answered without login for the 23 files fetched on 2026-11-10, but the portal's own
  download button needs a login, and the portal shows only a copyright line. Using it for thousands of files is the
  user's decision.
* **Licence consequence:** private, non-commercial analysis and model fitting are covered. Mirrored, simplified or
  otherwise derived traces cannot be shared. Whether model weights fitted on the traces count as adapted material is
  unresolved.

## 7. Sealing

* Traces are stored per part: `data/raw/jiao2025/traces/{train,validation,test}/<animal>/`. The fetch script places a
  file by its animal's part in the manifest and records bytes and SHA-256.
* The engine reads traces and regional rows only through loaders that take the manifest and a stage. Stage `develop`
  can open training animals only; `validate` adds validation animals **as queries**; `test` adds test animals as
  queries and requires the frozen protocol, pinned hashes and a committed tree, and writes an opened marker first.
* Regional outcome rows of validation and test animals are never loaded into any fit. They are used only to score
  those animals, at their stage.
* Canaries: a planted leak (a query animal put into a state's animal list) must make the engine refuse; replacing a
  held-out animal's traces and regional rows by random values must leave its predictions byte-identical.

## 8. Query-safe contract for any regional model and geometry decoder

A **regional state** is anything fitted on region-level or trace data of a set of animals: parameters, latent factors,
cluster assignments, normalisation statistics, tuned settings, donor weights.

1. Every regional state is fitted on **training animals only**. Validation and test animals are in no state, at any
   time, including the states that produce donor responsibilities.
2. The prediction for a query neuron of animal a comes from states that do not contain a. For validation and test
   queries that is the training state. For training animals that need their own regional prediction as a decoder
   input, it is an out-of-fold state by animal.
3. Donors are traces of training animals other than a. Donor responsibilities are computed from the query's
   permitted inputs and a state that does not contain a.
4. **Permitted query inputs:** soma position, source label, hemisphere, soma label, reconstruction type, and
   properties inferred from them by models that satisfy rule 1. **Never inputs:** the query's own axon, dendrite, or
   regional profile; these are used for scoring only.
5. Settings are chosen on validation animals only. A held-out set is scored once; nothing is tuned afterwards.
6. Audit trail: a state manifest (state ID, hash of the sorted animal list, code hash) and a prediction log (query,
   state ID, hash of the donor animal set), with an automatic check that the query's animal is in neither.

## 9. Mechanical inclusion rules after fetching (fixed now, independent of any model)

1. Tree: node 1 is the soma (type 1, parent -1); every parent exists; one connected tree; no cycle.
2. Soma within 5 um of the table's soma after the frame rule.
3. Soma hemisphere equals the index.
4. Node types within {0, 1, 2, 3, 4, 10, 11, 12}; the axon is coded as type 2 or type 0, not both.
5. `Axon_only` files have no dendrite node; `Axon_and_dendrite` files have at least one, otherwise the neuron is
   excluded from scoring (its dendrite could be counted as axon).
6. Axon cable of at least 1,000 um.
7. Cable outside the atlas box (AP 0 to 13,200, DV 0 to 8,000, ML 0 to 11,400) is clipped for scoring, as in the
   authors' own measurements, and the clipped length is logged.

Every exclusion is logged with its rule and applies to all methods alike; counts are reported per part.

## 10. What a result on this cohort can and cannot say

* The unit is the animal: 40 test animals. A method must beat copying a representative measured neuron on the
  animal-mean geometric distance **and** on F1, with a matched uniform-attachment control.
* Labels with enough test animals for a per-label statement: Adcyap1 6, Sst 5, Penk 5, Nts 4, Pdyn 4, Tac2 4.
  Oxt, Pomc and Tac1 have 2; Orexin, Agrp, Crh, Pmch, Avp and Trh have 1; Vip has none.
* For Avp and Trh the training pool holds very little of the label (13 neurons in 2 animals; 42 neurons in 1 animal).
* 111 of the 1,547 test neurons are axon-only reconstructions.
* A label effect is confounded with the animal and its batch (one label per animal).
* Other real datasets (`reports/axon_dataset_other_sources.md`): MouseLight has 97 hypothalamic neurons from 22
  samples of another laboratory (21 CC BY 4.0, about 16.5 MB; the rest CC BY-NC), usable later as an independent
  check; Li et al. 2024 adds one new animal (18696, 115 oxytocin neurons; archive CC BY 4.0 at the data centre).

## 11. Open points

* The development list must be confirmed by the other agent's manifest.
* Axon coding (type 2 or 0) of the 134 animals without a trace on disk is unknown until fetched.
* Whether portal, data-centre and paper versions of a trace are identical; meaning of node types 10 to 12.
* Whether the data centre's terms extend to the portal.
* The source notes were not independently reviewed.

**Decision needed before any fetch:** which route (data-centre account, or the portal's per-neuron address), and the
go-ahead for stage A.

## Appendix: the 166 animals

`status`: development = one of the 14 development animals; inspected = a trace was inspected, not development;
new = no trace ever inspected.

| animal | source label | neurons | part | status | archive bytes |
|---|---|---:|---|---|---:|
| 192710 | Pomc | 38 | train | inspected | 23,594,115 |
| 192711 | Pomc | 20 | validation | new | 12,444,164 |
| 192712 | Pomc | 23 | test | new | 11,887,298 |
| 192713 | Pomc | 23 | test | new | 19,497,792 |
| 192714 | Pomc | 11 | train | new | 9,337,842 |
| 195037 | Oxt | 58 | train | new | 17,202,219 |
| 195038 | Oxt | 48 | train | inspected | 9,671,750 |
| 195039 | Oxt | 62 | train | inspected | 17,990,544 |
| 201221 | Nts | 33 | train | inspected | 43,108,652 |
| 201223 | Nts | 86 | validation | new | 103,850,382 |
| 201224 | Nts | 110 | validation | new | 121,709,427 |
| 201225 | Nts | 92 | test | new | 115,615,958 |
| 201226 | Nts | 76 | test | new | 81,025,055 |
| 201227 | Nts | 37 | train | new | 41,423,457 |
| 201228 | Nts | 51 | test | new | 54,096,137 |
| 201229 | Tac1 | 118 | train | new | 109,574,143 |
| 201230 | Tac1 | 68 | train | development | 73,747,532 |
| 201232 | Tac1 | 54 | train | development | 61,188,318 |
| 201578 | Nts | 50 | train | new | 39,862,840 |
| 201580 | Nts | 81 | train | new | 115,334,088 |
| 201789 | Pdyn | 47 | test | new | 26,140,674 |
| 201793 | Pdyn | 56 | train | new | 49,595,158 |
| 201795 | Pdyn | 36 | validation | new | 24,171,071 |
| 201796 | Tac1 | 70 | train | new | 102,726,529 |
| 201797 | Tac1 | 76 | validation | new | 97,310,012 |
| 201892 | Adcyap1 | 20 | test | new | 20,154,606 |
| 201893 | Adcyap1 | 59 | train | new | 50,806,942 |
| 201894 | Adcyap1 | 44 | train | development | 55,884,214 |
| 201895 | Tac2 | 26 | train | new | 24,187,024 |
| 201896 | Tac2 | 13 | train | new | 11,284,434 |
| 201922 | Nts | 60 | train | inspected | 63,233,576 |
| 201924 | Nts | 135 | train | new | 139,313,803 |
| 201925 | Nts | 69 | validation | new | 71,168,204 |
| 201926 | Penk | 19 | train | new | 11,160,131 |
| 201929 | Sst | 30 | test | new | 26,819,602 |
| 202166 | Penk | 45 | test | new | 22,391,940 |
| 202167 | Penk | 10 | test | new | 5,809,249 |
| 202168 | Penk | 24 | train | new | 31,653,156 |
| 202170 | Penk | 71 | train | development | 54,065,799 |
| 202171 | Penk | 84 | validation | new | 78,957,744 |
| 202172 | Sst | 11 | validation | new | 7,868,564 |
| 202174 | Sst | 47 | train | new | 37,335,572 |
| 202175 | Sst | 9 | train | new | 4,841,898 |
| 202275 | Tac1 | 23 | test | new | 33,407,691 |
| 202276 | Tac1 | 10 | train | new | 12,269,510 |
| 202455 | Penk | 54 | test | new | 71,911,971 |
| 202456 | Penk | 25 | train | new | 38,133,795 |
| 202457 | Pomc | 53 | train | inspected | 46,854,028 |
| 202458 | Pomc | 51 | train | new | 47,871,066 |
| 202459 | Pomc | 45 | train | new | 31,570,175 |
| 202460 | Pomc | 40 | train | new | 31,360,063 |
| 202557 | Sst | 33 | validation | new | 35,455,930 |
| 202558 | Sst | 36 | test | new | 37,852,901 |
| 202559 | Tac2 | 45 | train | development | 49,978,425 |
| 202560 | Tac2 | 30 | test | new | 34,722,695 |
| 202561 | Tac2 | 41 | test | new | 49,209,712 |
| 202562 | Tac2 | 73 | train | new | 66,363,798 |
| 202564 | Tac2 | 18 | test | new | 23,233,627 |
| 202570 | Penk | 14 | train | new | 18,664,951 |
| 202571 | Pdyn | 140 | test | new | 150,210,256 |
| 202572 | Pdyn | 87 | validation | new | 93,713,310 |
| 202600 | Oxt | 5 | train | new | 4,396,799 |
| 202601 | Oxt | 110 | train | new | 32,722,147 |
| 202602 | Oxt | 28 | test | new | 22,566,361 |
| 202603 | Oxt | 21 | test | new | 22,216,956 |
| 202614 | Avp | 3 | train | new | 1,010,189 |
| 202616 | Sst | 29 | train | new | 47,535,219 |
| 202617 | Sst | 109 | train | development | 98,223,108 |
| 202618 | Sst | 68 | train | inspected | 63,228,169 |
| 202619 | Sst | 110 | train | new | 98,961,255 |
| 202620 | Sst | 74 | train | development | 70,767,989 |
| 202662 | Orexin | 2 | train | new | 380,730 |
| 202792 | Pdyn | 22 | train | new | 19,965,205 |
| 202793 | Pdyn | 15 | train | new | 13,743,965 |
| 202812 | Agrp | 28 | train | new | 14,415,787 |
| 202813 | Agrp | 32 | train | inspected | 13,861,264 |
| 202814 | Agrp | 5 | test | new | 2,819,613 |
| 202815 | Agrp | 27 | train | new | 16,662,955 |
| 202816 | Agrp | 10 | validation | new | 4,926,588 |
| 202817 | Oxt | 2 | train | new | 599,741 |
| 202818 | Oxt | 2 | train | new | 2,307,943 |
| 202820 | Oxt | 3 | train | new | 1,144,409 |
| 202822 | Vip | 3 | train | new | 2,817,383 |
| 202824 | Vip | 13 | train | development | 12,325,923 |
| 202825 | Avp | 10 | train | inspected | 5,626,970 |
| 202826 | Avp | 28 | test | new | 8,575,659 |
| 210004 | Penk | 1 | train | new | 1,154,575 |
| 210006 | Penk | 37 | train | new | 35,547,739 |
| 210007 | Pmch | 93 | train | development | 88,110,802 |
| 210016 | Penk | 48 | train | new | 44,706,964 |
| 210073 | Penk | 42 | train | inspected | 54,404,178 |
| 210074 | Penk | 20 | test | new | 17,486,882 |
| 210104 | Orexin | 29 | train | development | 39,524,540 |
| 210105 | Orexin | 12 | train | inspected | 7,793,167 |
| 210245 | Oxt | 9 | validation | new | 4,792,575 |
| 210247 | Pdyn | 42 | train | new | 42,910,264 |
| 210660 | Penk | 23 | train | new | 35,204,226 |
| 210672 | Adcyap1 | 19 | test | new | 19,221,230 |
| 210673 | Adcyap1 | 39 | validation | new | 38,470,583 |
| 210674 | Adcyap1 | 13 | test | new | 11,015,395 |
| 210757 | Pdyn | 24 | train | new | 15,844,569 |
| 210758 | Pdyn | 4 | train | new | 2,914,827 |
| 210796 | Pdyn | 3 | train | development | 949,450 |
| 210797 | Pdyn | 62 | test | new | 51,654,586 |
| 210951 | Pmch | 10 | train | inspected | 6,680,386 |
| 210952 | Pmch | 4 | train | new | 2,135,725 |
| 210953 | Pmch | 9 | test | new | 9,586,934 |
| 210975 | Tac1 | 66 | train | new | 75,337,892 |
| 210989 | Orexin | 26 | train | new | 27,567,620 |
| 210990 | Orexin | 15 | test | new | 14,641,227 |
| 210991 | Orexin | 19 | validation | new | 31,560,939 |
| 210992 | Adcyap1 | 44 | train | new | 50,453,136 |
| 210993 | Adcyap1 | 35 | test | new | 33,976,484 |
| 210994 | Adcyap1 | 34 | train | new | 33,846,702 |
| 210995 | Adcyap1 | 48 | validation | new | 45,744,441 |
| 210996 | Adcyap1 | 10 | train | new | 6,400,583 |
| 210997 | Adcyap1 | 52 | train | new | 56,835,390 |
| 211028 | Sst | 38 | train | new | 30,656,809 |
| 211266 | Pdyn | 55 | train | inspected | 56,871,819 |
| 211267 | Pdyn | 38 | test | new | 25,304,860 |
| 211268 | Pdyn | 63 | train | new | 39,298,606 |
| 211269 | Tac2 | 72 | test | new | 63,355,037 |
| 211270 | Tac2 | 27 | train | new | 15,901,868 |
| 211271 | Adcyap1 | 38 | train | new | 22,944,261 |
| 211272 | Adcyap1 | 22 | validation | new | 20,378,614 |
| 211273 | Adcyap1 | 31 | test | new | 15,353,159 |
| 211275 | Adcyap1 | 49 | train | new | 45,206,314 |
| 211276 | Adcyap1 | 84 | train | new | 68,391,180 |
| 211353 | Adcyap1 | 121 | train | development | 99,595,067 |
| 211354 | Adcyap1 | 81 | test | new | 70,212,098 |
| 211355 | Adcyap1 | 48 | train | new | 51,430,646 |
| 211392 | Adcyap1 | 31 | validation | new | 33,300,265 |
| 211393 | Adcyap1 | 126 | train | new | 96,312,981 |
| 211458 | Tac1 | 41 | test | new | 25,348,727 |
| 211460 | Nts | 86 | train | new | 97,250,327 |
| 211461 | Nts | 47 | train | inspected | 41,790,845 |
| 211462 | Nts | 13 | train | new | 9,498,340 |
| 211463 | Nts | 13 | test | new | 12,995,232 |
| 211464 | Tac2 | 33 | validation | new | 41,338,902 |
| 211465 | Tac2 | 27 | train | inspected | 45,109,551 |
| 211491 | Tac2 | 7 | validation | new | 9,915,331 |
| 211493 | Tac2 | 111 | train | new | 128,152,981 |
| 211494 | Penk | 59 | test | new | 62,286,305 |
| 211525 | Penk | 58 | validation | new | 55,702,889 |
| 211526 | Penk | 15 | validation | new | 15,179,929 |
| 211551 | Sst | 72 | test | new | 61,941,419 |
| 211552 | Sst | 64 | train | new | 56,585,742 |
| 211553 | Sst | 125 | train | new | 85,041,220 |
| 211554 | Sst | 124 | validation | new | 101,155,927 |
| 211555 | Sst | 59 | train | new | 58,673,105 |
| 211556 | Tac2 | 2 | train | new | 2,066,230 |
| 211557 | Tac2 | 23 | train | new | 23,912,941 |
| 211558 | Tac2 | 91 | validation | new | 113,837,208 |
| 211559 | Tac2 | 24 | train | new | 32,464,912 |
| 211595 | Sst | 72 | train | inspected | 53,249,473 |
| 211596 | Sst | 42 | test | new | 33,691,860 |
| 211597 | Sst | 20 | validation | new | 13,767,249 |
| 211600 | Sst | 9 | train | new | 6,729,756 |
| 211601 | Sst | 14 | test | new | 8,141,590 |
| 211889 | Pdyn | 121 | train | development | 98,543,556 |
| 211890 | Pdyn | 40 | train | inspected | 47,189,944 |
| 212201 | Trh | 57 | test | new | 41,664,998 |
| 212202 | Trh | 42 | train | development | 28,977,766 |
| 212247 | Crh | 13 | test | new | 3,119,348 |
| 212248 | Crh | 16 | train | inspected | 5,490,346 |
| 220102 | Crh | 47 | validation | new | 24,430,575 |
