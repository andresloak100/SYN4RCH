# Other real single-neuron axon datasets: metadata audit

Task `other_datasets` of the dataset audit for a larger real 3D benchmark. Written 2026-11-10 (local time;
web reads logged 2026-11-11 01:10 to 01:48 UTC). This is an inventory, not a model run.

**Rules followed.** Metadata and documentation only. No SWC or other geometry file was requested, not even
its headers. No download button was clicked, no login, no account, no terms accepted. From every portal index
only the fields `file, id, acronym, name, type_name, type_array, parent_type_name, hemisphere,
reconstruction_type` were kept; projectome type and class fields were never printed, counted or used. Nothing
in the repository was changed except this new file. The other agent's files and the sibling audit files were
not opened.

**Disclosures.** (1) The public neuron indexes of the eight other mouse projects on the portal and seven
animal tables were fetched (about 14.2 MB of JSON in all; the `pvh_oxt` pair was already in the session
cache); they are lists of file names, soma labels and animal records. (2) One large listing was fetched:
NeuroXiv's `NeuroXiv_metadata.csv` (19.4 MB, 18 columns: identifiers, source, atlas, soma position,
hemisphere, soma region, arbor flags; no tracing nodes, no projection measurements). (3) The MouseLight
search service was asked for all neurons; the answer holds identifiers and one soma record per tracing, no
tracing nodes. (4) Two generic keywords ("projectome", "hypothal") were typed into the data centre's public
catalogue search. (5) Soma coordinates of non-Jiao neurons were used only to tell axes and hemisphere apart.

**Status tags.** VERIFIED = read from a source listed in section 10 (hash given) or computed from one.
INFERRED = reasoning from verified facts, stated by no source. UNVERIFIED = not checked, or impossible to
check under the rules.

## 1. In plain words

- **A second real set of hypothalamic neurons from mice Jiao never used is on the same portal**: 266 neurons
  from 5 mice (Ding et al. 2025, "somatosensory ascending pathways"). No genetic label; lateral hypothalamus
  only; same institute and pipeline as Jiao.
- **MouseLight (Janelia, another laboratory) has 97 hypothalamic neurons from 22 mice**, all with axon and
  dendrite. Best independence, but few cells per mouse, another file convention and mixed licences.
- **The oxytocin set (Li et al. 2024) adds one new mouse** (18696, 115 neurons). Its archive carries CC BY 4.0
  at the data centre. Its other two mice are Jiao animals.
- **The portal's "LHA" project is not new data**: it is exactly Jiao's 916 LHA neurons, file for file.
- Nothing else found has downloadable hypothalamic neurons with long axons. The cortex, hippocampus, amygdala
  and spinal sets (about 31,000 more neurons on the portal) could only be generic training geometry.
- The data centre (BSDC) needs an account to download. The portal serves single files without one (shown
  earlier for Jiao; not tested for the new set).
- Proposed external pool, before anything is fetched: **478 neurons from 28 mice outside Jiao's 166**
  (266 / 5 + 97 / 22 + 115 / 1). Section 7 lists exactly what a fetch would be.

## 2. Summary table

HY = soma label is `HY` or a descendant of it in `data/raw/atlas/structure_graph_1.json` (91 acronyms; ZI, STN
and the mammillary body are inside). Use: **T** extra training geometry, **X** independent test, **N** neither.
All counts VERIFIED from the indexes named in section 10 unless marked.

| dataset | neurons / animals | soma in HY: neurons / animals | label per neuron | animals shared with Jiao | licence shown | use |
|---|---|---|---|---|---|---|
| Ding 2025, portal `spcd` | 2,249 / 79 | **266 / 5** | none | 0 | BSDC: CC BY-NC-ND 4.0; portal: copyright line | X (the 266), T (rest) |
| Li 2024, portal `pvh_oxt` | 264 / 3 | 232 / 3 (mouse 18696: 101 of 115) | Oxt (all) | **2** (195037, 195038) | BSDC: CC BY 4.0; portal: copyright line | X, small (18696 only) |
| Ji 2025, portal `lha` | 916 / 89 | 916 / 89 | 12 lines | **89 (all); every file is a Jiao file** | as Jiao | N (duplicate) |
| CeA, portal `cea` (Dou 2024 + a 2026 set) | 488 / 20 | 0 | SST, Crh, PKC, MOR | 0 | BSDC: CC BY-NC-ND 4.0 | T |
| Gao 2022 / 2023, portal `pfc` | 6,357 / 161 | 0 | none | 0 | BSDC: CC BY-NC-ND 4.0 | T |
| Qiu 2024, portal `hipp` | 10,100 / 240 | 0 | none | 0 | BSDC: CC BY-NC-ND 4.0 | T |
| Whole cortex, portal `whole-cortex` | 18,621 / 382 (includes the 6,357 PFC) | 0 | wild type or 9 Cre lines | 0 | BSDC: CC BY-NC-ND 4.0 | T |
| Retina, portal `rbm` | 2,133 / 95 | 0 (retina) | WT, Opn4-cre, Jam2-cre | 0 | not released | N |
| MouseLight (Winnubst 2019) | 1,653 (1,641 listed) / 55 | **97 / 22** | none in the public index | 0 (other laboratory) | CC BY-NC 4.0 or CC BY 4.0, by record | X, T |
| SEU-ALLEN / BICCN fMOST (Peng 2021) | 1,876 / 39 | 0 by source label | not in the listing read | 0 | BIL: CC BY-SA 4.0; NeuroXiv copy: CC BY-NC | T |
| NeuroXiv (aggregator) | 175,149 rows, 19,406 with axon | MouseLight 55 to 70; SEU-ALLEN 27 border cells; ION 0 | none | none in its listing | CC BY-NC, no redistribution | a route, not a source |
| NeuroMorpho.Org | 1,037 mouse HY records, 13 archives | 66 MouseLight with complete axon; rest are slices | varies | none | CC BY 4.0 | N |
| Lin 2023 (histamine neurons) | 60 | 60 (TMN) | HDC-CreERT2 | UNVERIFIED | no deposit found | N |
| Senol 2025 (tuberal Sst) | 186 | 186 | SST-Cre | UNVERIFIED (Jiao's group did the imaging) | no deposit found | N |

## 3. The other projects on mouse.digital-brain.cn

VERIFIED: the public project list (`api/v1/public/projects`, 7,821 bytes, same hash as the session cache)
names 12 projects. Mouse: prefrontal cortex (Gao 2022), hippocampus (Qiu 2024), CeA (Dou 2024), PVH OXT
(Li 2024), hypothalamus (Jiao 2025), LHA (Ji 2025), SPCD (Ding 2025), whole cortex (Gao 2025), retina (2026).
Others: two zebrafish projects and macaque PFC. Each mouse project has a neuron index
`info/mouse/<project>/mouse.neuron.info.json` and most have an animal table.

VERIFIED: across the seven mouse projects other than `hy` and `lha` there are 819 distinct animal IDs; the
only ones also in Jiao's 166 are 195037 and 195038 (through `pvh_oxt`). None of the 37 inspected neurons and
none of the 16 development neurons is in any project other than `lha`.

### 3.1 Ding et al. 2025: somatosensory ascending pathways (`spcd`)

- **Citation.** VERIFIED (CrossRef): Ding W-Q, Song W, Shi X, ... Sun Y-G. Single-neuron projectome reveals
  organization of somatosensory ascending pathways in the mouse brain. Neuron 113:2083-2101.e5 (2025),
  doi:10.1016/j.neuron.2025.03.007. Data record: BSDC, doi:10.12412/BSDC.1737020663.20001,
  https://cstr.cn/33145.11.BSDC.1737015588.1879450714607951874.
- **Counts.** VERIFIED (index): 2,249 neurons, 79 animals: 785 spinal projection neurons (label `Spcd`,
  53 animals) and 1,464 brain "central relay" neurons (26 animals). The record's abstract says 785 and 1,463.
- **Hypothalamus.** VERIFIED: 266 relay neurons have a HY-subtree soma label, all from 5 animals:

| animal | HY-soma neurons (of all its neurons) | labels | hemisphere field | reconstruction |
|---|---|---|---|---|
| 210726 | 77 (155) | LHA 28, HY 27, ZI 16, PSTN 3, PH 2, PST 1 | right 46, left 31 | axon and dendrite 77 |
| 210661 | 73 (169) | HY 28, LHA 26, PSTN 7, ZI 7, PST 4, TU 1 | left 73 | axon and dendrite 71, axon only 2 |
| 210662 | 58 (98) | LHA 26, HY 18, PSTN 10, ZI 4 | left 58 | axon and dendrite 58 |
| 211986 | 36 (62) | ZI 15, HY 11, LHA 6, PeF 2, MM 1, PSTN 1 | left 36 | axon and dendrite 36 |
| 210727 | 22 (43) | HY 11, LHA 7, ZI 2, PSTN 1, TM 1 | right 15, left 7 | axon and dendrite 22 |

  Totals: HY 95, LHA 93, ZI 44, PSTN 22, PST 5, PH 2, PeF 2, TU 1, TM 1, MM 1. Declared neighbours of the
  hypothalamus in the same set: PAG 99, VTA 1 (a list declared in this audit, not measured).
- **Animals.** VERIFIED (animal table `spcd/SPCD.json`, 75 rows): the five are C57BL/6J, laboratory field
  "孙衍刚" (Sun Yan-Gang), injected at two sites, "LHA;PO" (210661, 210662) or "LHA;VPL" (210726, 210727,
  211986), with AAV2/9-hEF1a-DIO-EYFP; age 65 to 68 days (211986: 179). Sex: the field is `0`; its meaning is
  UNVERIFIED. None of the five IDs is in the Jiao index, in the 32 inspected animals or in any other portal
  project.
- **What these cells are.** VERIFIED (abstract): relay neurons that receive direct spinal input. INFERRED: a
  Cre-dependent reporter in wild-type mice means the Cre came across the synapse from the spinal cord, so the
  sample is chosen by input and by injection site, not by peptide. It covers the lateral hypothalamus, zona
  incerta and parasubthalamic region only.
- **Label.** VERIFIED: `type_name` is empty and `type_array` is `unknown` for all 2,249. No source label.
- **IDs.** VERIFIED: `<animal>_<number>`, file `<animal>_<number>.swc`, as in Jiao.
- **Frame.** VERIFIED: the soma labels include PeF, PoT, Eth and PIL. INFERRED (atlas history not re-read
  here): those structures entered the Allen ontology with the 2017 CCFv3 annotation, so the labels were
  assigned with that version. INFERRED: micrometres in CCFv3 2017, axis order AP, DV, ML, true hemisphere, as
  in Jiao's portal files: same viewer, same file specification on the portal's Download page, both
  hemispheres present in the index. UNVERIFIED until a file is opened. The registration method is
  UNVERIFIED: the paper's Methods are behind a paywall (the page offers only "Get full text access"). The
  institute's published protocol (Ren et al. 2025, Bio-protocol, doi:10.21769/BioProtoc.5300) describes two
  independent tracers plus a third who merges, a quality score of at least 90, and ANTs registration of
  region masks to CCFv3; INFERRED, not stated, that this set followed it.
- **Node types.** VERIFIED (portal Download page): 0 undefined, 1 soma, 2 axon, 3 dendrite; seven columns;
  parent -1 at the soma. UNVERIFIED whether these files, like Jiao's, sometimes write axon as type 0 or carry
  types 10 to 12.
- **Licence and access.** VERIFIED: the BSDC record shows "License Agreement CC BY-NC-ND 4.0", badge
  "Open Access", 10 files, 5.86 GB; the geometry is `/somatosensory _data/data/swc/brain.zip` (1,086,478 KB)
  and `.../swc/spcd.zip` (1,288,724 KB). VERIFIED (portal code): the viewer fetches single files from
  `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/spcd/spcd_neuron_download/<animal>_<number>.swc`
  (the code builds it from the base, `/info/mouse/` and the project part, exactly as for Jiao's `hy` files).
  UNVERIFIED: whether that address answers without login (the Jiao address does). The portal shows only
  "Copyright © 2023 Digital Brain. All Rights Reserved."
- **Independence.** VERIFIED (CrossRef author lists): Sun Yan-Gang, last author here, is also an author of
  Jiao 2025; Xiaofei Wang, who deposited the Jiao files, is an author here. Same institute, same imaging and
  tracing platform, overlapping authors: new animals, not a new laboratory.
- **Verdict: independent test set (X) for the 266**, from genuinely new animals, with three limits: same
  institute and pipeline as Jiao, no label, lateral regions only. The other 1,983 neurons: training only.

Exact IDs of the 266 (number after the animal; `*` = axon only):

```
210661 HY 001 004 006 007 008 021 032 034 035 036 040 041 071* 083 084 085 086 087 089 091 092 096 097 098 111 115 117 224
210661 LHA 005 014 015 016 029 031 033 037 038 039 082 093 094 095 099 100 101 105 181 182 183 184 185 186 193 225
210661 PSTN 003 009 010 011 012 013 107 | ZI 106 118 151 152 161* 162 195 | PST 088 104 194 226 | TU 102
210662 LHA 024 025 026 027 028 031 032 033 034 035 036 041 043 046 050 051 052 055 056 057 060 061 089 090 091 097
210662 HY 002 003 005 014 042 044 045 047 048 049 053 062 079 085 094 095 096 101
210662 PSTN 029 030 074 075 076 077 078 086 092 093 | ZI 006 084 088 104
210726 LHA 001 003 004 010 011 014 018 020 021 022 023 024 025 026 028 029 030 034 036 053 054 146 147 148 161 163 170 171
210726 HY 006 007 008 013 015 017 027 031 032 035 049 055 064 090 092 095 101 102 125 126 129 138 139 151 152 156 157
210726 ZI 009 043 044 046 047 048 051 052 056 098 118 119 124 162 167 173 | PSTN 002 037 160 | PH 016 127 | PST 033
210727 HY 007 009 012 013 016 025 026 035 036 044 060 | LHA 014 015 037 038 039 040 041 | ZI 008 046 | PSTN 010 | TM 017
211986 ZI 019 020 022 023 029 030 039 044 045 046 048 050 067 088 095 | HY 025 027 032 040 043 049 073 082 106 112 113
211986 LHA 038 070 079 083 085 104 | PeF 080 081 | MM 041 | PSTN 109
```

### 3.2 Li et al. 2024: PVH oxytocin neurons (`pvh_oxt`)

- **Citation.** VERIFIED: Li H, Jiang T, An S, Xu M, ... Yu X. Neuron 112:1081-1099.e7 (2024),
  doi:10.1016/j.neuron.2023.12.022. The article is now open archive and was read in a browser.
- **Counts.** VERIFIED (index and paper): 264 neurons, all Oxt, all `Axon_only`: 115 from 18696, 65 from
  195037, 84 from 195038. Soma labels: PVH 149, HY 59, TH 28, PVHd 23, RE 4, SBPV 1.
- **Overlap.** VERIFIED: 195037 and 195038 are Jiao animals (58 and 48 Jiao neurons); 195038 is one of the 32
  inspected animals (Jiao neuron 195038_041), neither is a development animal. 18696 is not a Jiao animal.
  No file stem equals a Jiao stem (Li stems have three parts, e.g. `18696_000_L`; the `_L` is on all 264
  files, Left and Right alike, so it is not a hemisphere code).
- **Are the Li and Jiao cells in the two shared mice the same cells?** Still UNVERIFIED. New hint, INFERRED:
  the portal's animal table gives 123 traced cells for 195037 and 133 for 195038; Li has 65 and 84, Jiao 58
  and 48, and 65 + 58 = 123, 84 + 48 = 132. That fits two complementary sets better than one shared set.
- **Sex conflict.** VERIFIED: the paper describes 18696 as a lactating female and 195037 and 195038 as virgin
  females. Jiao 2025 is titled "in male mice" and its Reporting Summary says male only
  (`reports/d109_jiao_conventions.md`). The two shared animals cannot be both. UNVERIFIED which is right.
- **Frame.** VERIFIED (paper): CCFv3; features drawn in Amira on the 10 µm counterstain channel, files warped
  with the deformation field; the authors say registration was weaker around the ventricles and below the
  third ventricle, and that the registration in the paper differs slightly from the website's. This agrees
  with the 130 to 350 µm soma shifts already noted for the shared brains. Both hemispheres appear in the
  index (18696: Left 108, Right 7), so the portal copy is not mirrored (INFERRED).
- **Tracing.** VERIFIED (paper): manual, FNT software, two independent tracers and a third to merge; dendrites
  only for 195038. The portal index marks all 264 as axon only, so the served files may lack those dendrites
  (UNVERIFIED).
- **Licence and access.** VERIFIED: the paper's Key Resources Table gives the deposit
  doi:10.12412/BSDC.1704691586.20001. That BSDC record ("Paraventricular nucleus oxytocin neurons", Xiang Yu,
  2024-01-03) shows **"License Agreement CC BY 4.0"**, one file `/PVH.rar` (39,628 KB), badge "Open Access".
  The portal copy has only the copyright line. Single files (portal code):
  `.../info/mouse/pvh_oxt/pvh_oxt_neuron_download/<stem>.swc`. The earlier D-109 check read the first 500
  bytes of six Li files by range request (`reports/d109_jiao_conventions.md`); INFERRED: through this address
  and without login.
- **Verdict: a small independent check (X) for mouse 18696 only**: one animal, one label, one nucleus, axon
  only, same institute, another registration. The 149 neurons of 195037 and 195038 are extra traces in two
  Jiao animals and must follow those animals wherever the split puts them.

Exact IDs of mouse 18696 (number between `18696_` and `_L`):

```
PVH  000 005 006 008 011 013 015 016 017 018 024 025 026 030 033 036 037 040 045 046 048 049 050 060 061 070 072 073 074 077 085 086 087 088 090 091 092 093 094 095 096 098 100 103 104 106 109 110 112 113 115 118 121 122 123 124 125 127 129 130 132 133 134 135 136 138 142
HY   004 007 009 012 014 019 022 023 027 034 035 052 065 078 080 081 089 097 101 102 105 107 108 111 114 116 128 131 137 141
TH   020 031 043 044 062 064 069 071 075 084 126 139 | PVHd 032 059 066 | RE 021 063 | SBPV 140
```

### 3.3 Ji et al. 2025: "Mouse LHA" (`lha`)

VERIFIED: the paper is Ji C, Zhang Y, Lin Z, ... Nat Commun (2025), doi:10.1038/s41467-025-67133-5; it uses
916 sparsely labelled LHA neurons of 12 lines and points to the portal. The `lha` index has 916 records, and
all 916 file stems, soma labels, lines, hemispheres and reconstruction types equal the Jiao index's 916
LHA-soma records; the portal serves them from the Jiao address. It holds 5 of the 37 inspected neurons and 2
of the 16 development neurons because they are those same files. **Verdict: neither; no new neuron, no new
animal.**

### 3.4 Central amygdala (`cea`)

VERIFIED: 488 neurons, 20 animals, soma CEAc 227, CEAm 155, CEAl 106; labels SST 206 (8 animals), Crh 119
(4), PKC 96 (4), MOR 67 (4); 438 axon and dendrite, 43 axon only, 6 dendrite only, 1 undefined. The 67 MOR
neurons are Dou Y-N et al., Natl Sci Rev 11:nwae195 (2024), doi:10.1093/nsr/nwae195 (article CC BY 4.0);
the BSDC record says they were registered to CCFv3 and is CC BY-NC-ND 4.0 (one zip, 21,597 KB). The other 421
are the BSDC record "Connectome of SST, CRH, and PKC neurons in the central amygdala of mice" (2026-05-06,
421 files, 430.19 MB, CC BY-NC-ND 4.0, no article listed). No animal shared with Jiao. **Verdict: training
only.** It is the nearest in kind to the hypothalamus (subcortical, peptide-labelled, long descending axons).

### 3.5 Prefrontal cortex, hippocampus, whole cortex, retina

| project | citation (VERIFIED, CrossRef) | neurons / animals | reconstruction | BSDC record (VERIFIED) |
|---|---|---|---|---|
| `pfc` | Gao L et al., Nat Neurosci 25:515 (2022), doi:10.1038/s41593-022-01041-5; dendrites: 26:1111 (2023), doi:10.1038/s41593-023-01339-y | 6,357 / 161 | axon only 4,437; axon and dendrite 1,920 | doi:10.12412/BSDC.1690164952.20001, one zip 11,136,078 KB; older record 163 zips, 11.63 GB |
| `hipp` | Qiu S et al., Science 383:eadj9198 (2024), doi:10.1126/science.adj9198 | 10,100 / 240 | axon only 5,863; axon and dendrite 4,237 | doi:10.12412/BSDC.1667278800.20001, `/hip.zip` 4,737,786 KB |
| `whole-cortex` | none listed ("Gao et al., 2025", paper link "TBD") | 18,621 / 382 = the 6,357 PFC + 12,264 | axon only 4,850; axon and dendrite 13,771 | doi:10.12412/BSDC.1747279998.20001, `/swc_in_ccf.zip` 6,803,774 KB (the 12,264 only) |
| `rbm` | none ("TBD et al, 2026") | 2,133 / 95, soma label `R` (retina) | axon and dendrite | "coming soon" |

All CC BY-NC-ND 4.0 at BSDC. No hypothalamic soma, no animal shared with Jiao. Labels: none for `pfc` and
`hipp`; `whole-cortex` has wild type 8,760 and nine Cre lines (Wfs1 2,335, Rbp4 2,102, PlexinD1 1,519, Cux2
1,073, Tlx3 875, Tle4 679, Ntsr1 611, Rorb 559, Scnn1a 108). **Verdict: training only**, and of doubtful value
for hypothalamic shapes: these are cortical and hippocampal projection neurons.

### 3.6 Other records at the data centre

VERIFIED (BSDC catalogue searches "projectome", 22 hits, and "hypothal", 17 hits): the only hypothalamic
single-neuron records are the two Jiao records, the Li record (plus a second Li record holding figures only)
and the Ding record. Other single-neuron
records: nucleus accumbens (Open Access, 1.78 GB), cortical PT neurons (Open Access, 2.95 GB), S1 layer 2/3
(Open Access); and, not downloadable, VTA dopamine, BLA cholecystokinin, parafascicular thalamus and sensory
cortex ("Download by project members only"), mediodorsal thalamus ("Embargoed Access"). The portal code also
names unlisted projects `bla, cea2, cm, md, nac, pf, pvt, rt, vpl, vta`; they are not in the public list and
were not probed.

## 4. Janelia MouseLight

- **Citation.** VERIFIED: Winnubst J, ... Chandrashekar J. Cell 179:268-281 (2019),
  doi:10.1016/j.cell.2019.07.042 (free at PMC6754285).
- **Counts.** VERIFIED (public NeuronBrowser search service, version 1.6.4): `neuronCount` 1,653; a
  whole-brain soma query returned 1,641 neurons from 55 samples. One sample is one brain (the project's
  Resources page: about 100 neurons per sample).
- **Hypothalamus.** VERIFIED: **97 neurons from 22 samples** have a CCFv3 soma compartment in the HY subtree:
  LHA 29, ZI 15, HY 13, mammillary body 17 (Mmd 9, Mmme 7, Mml 1), PH 6, AHN 3, FF 3, STN 2, and one each of
  MPO, LZ, PST, PSTN, PeF, PVH, PVi, SBPV, VMH. Without ZI, FF and STN that is 77, the figure in
  `reports/conversion_map.md`. By the legacy CCFv2.5 assignment it is 95 neurons; 22 neurons change side of
  the HY boundary between the two atlases; 108 are in HY by at least one of three assignments. All 97 have
  both an axon and a dendrite tracing.
- **Animals.** VERIFIED: per sample, 17, 10, 9, 8, 8, 7, 5, 5, 4, 4, 3, then 2 or 1. VERIFIED (2019 paper):
  adult females, mostly wild-type C57BL/6 with sparse viral labelling; one Sim1-Cre brain (2018-08-01) and one
  Rorb-Cre brain. UNVERIFIED: sex and strain of the samples whose records were all issued after that paper
  (samples 31, 34, 37, 40 and 43: 14 of the 97 neurons); which public sample number is the Sim1 brain, so
  which neurons carry a Cre label (the public index exposes only a sample number; NeuroMorpho lists all 66 of
  its MouseLight hypothalamic neurons as C57BL/6).
- **IDs.** VERIFIED: neuron `AA0001`-style (20 newer ones `AL`), one DOI per neuron
  (`10.25378/janelia.<n>`), sample number 1 to 56.
- **Frame.** VERIFIED: the browser serves two frames, "CCFv3" and "CCFv2.5 (ML legacy)". VERIFIED (paper):
  registered to the Allen 10 µm average template with an affine step and landmark thin-plate splines in 3D
  Slicer; error typically under 100 µm; 110 ± 60 µm between repeated registrations. INFERRED from the soma
  coordinates returned by the CCFv3 query (x 1,241 to 10,124; y 463 to 7,192; z 2,045 to 13,164): micrometres
  with x = left-right, y = dorso-ventral, z = antero-posterior, so Jiao's (AP, DV, ML) is MouseLight's
  (z, y, x). Both sides present (hypothalamic somata: 76 with x < 5,700, 21 above), so not mirrored.
- **Node types.** VERIFIED (service constants): 0 path, 1 soma, 2 axon, 3 (basal) dendrite, 4 apical
  dendrite, 5 branch point, 6 end point; axon and dendrite are separate tracings. UNVERIFIED how an exported
  SWC writes them. VERIFIED (paper): every neuron traced by two annotators; consensus expected above 90% of
  the axon.
- **Licence, with the conflict.** VERIFIED: the NeuronBrowser terms page says the data "are available under
  the CC-BY NC license" and asks for the DOI of each neuron used. The per-neuron figshare records differ by
  year: of the 97, **76 are "CC BY-NC 4.0"** (2017 and 2019 records, which hold only a link) and **21 are
  "CC BY 4.0"** (2024 records, which hold `<id>.swc` and `<id>.json`). The AWS open-data entry for MouseLight
  says "Creative Commons Attribution 4.0 International". INFERRED for the whole index by DOI range: about
  1,096 older records, 526 from 2024, 19 with another DOI form.
- **Access.** VERIFIED: the figshare records and the browser are public pages with no login step; three
  routes exist: files attached to the 2024 figshare records; the NeuronBrowser's export for any neuron; a bulk
  archive "released 2022.02.28" that predates the 2024 additions (size UNVERIFIED). INFERRED, not tried: the
  files themselves download without login.
- **Verdict: independent test set (X)**, the only one from another laboratory, pipeline, microscope and
  registration; 22 animals but at most 17 neurons each; opposite sex to Jiao; no label. Also training (T) for
  its other 1,544 listed neurons. The price is a conversion step that must be checked on the first file.

Exact IDs of the 97 (`BY` = CC BY 4.0 record; all others CC BY-NC 4.0):

```
sample 5:  AA0080 PH, AA0081 PH, AA0082 HY, AA0084 HY, AA0086 HY, AA0087 HY, AA0090 VMH, AA0092 SBPV, AA0093 AHN
sample 7:  AA0128 ZI                      sample 8:  AA0141 AHN, AA0142 PVi
sample 9:  AA0152 AHN, AA0153 PVH, AA0163 HY, AA0167 PeF
sample 10: AA0174 ZI, AA0175 ZI, AA0177 ZI
sample 11: AA0207 LHA, AA0212 HY, AA0213 PH, AA0214 PH, AA0216 HY, AA0217 HY, AA0218 HY, AA0219 HY, AA0220 PH, AA0222 PH
sample 16: AA0330 ZI                      sample 18: AA0389 HY, AA0428 LHA
sample 19: AA0443 ZI, AA0449 HY, AA0486 ZI, AA0488 ZI          sample 20: AA0436 LHA, AA0437 LHA
sample 23: AA0537 LHA, AA0542 LHA, AA0555 LHA, AA0557 LHA, AA0558 LHA
sample 25: AA0606 LHA, AA0607 LHA, AA0608 LHA, AA0609 LHA, AA0614 ZI, AA0615 ZI, AA0616 LHA
sample 26: AA0730 LHA                     sample 27: AA0760 LHA, AA0761 LHA
sample 28: AA0768 ZI, AA0779 LHA, AA0808 ZI, AA0904 ZI, AA0917 LHA
sample 29: AA0930 Mmd, AA0933 Mmd, AA0934 Mmd, AA0935 Mmme, AA0936 Mmme, AA0937 Mmd, AA0938 PSTN, AA0939 Mml, AA0940 Mmd,
           AA1034 Mmd, AA1036 Mmme, AA1039 Mmd, AA1040 Mmd, AA1041 Mmd, AA1048 Mmme, AA1067 Mmme, AA1081 ZI
sample 30: AA1100 LHA, AA1101 PST BY, AA1107 HY BY, AA1120 LHA BY, AA1150 LHA BY, AA1153 LHA BY, AA1154 LHA BY, AA1175 STN BY
sample 31: AA1188 LHA BY, AA1198 Mmme BY, AA1200 LHA BY, AA1207 LZ BY, AA1209 FF BY, AA1218 LHA BY, AA1220 STN BY, AA1226 LHA BY
sample 34: AA1316 ZI BY                   sample 37: AA1412 LHA BY, AA1419 Mmme BY
sample 40: AA1462 FF BY, AA1464 FF BY     sample 43: AA1526 MPO BY
```

Note: Jiao's own Zenodo archive holds 18 MouseLight files (`AA*.swc`, `reports/swc_access.md`). Which 18 was
not read here (the listing is outside this task's inputs). They keep MouseLight's licence, not Zenodo's.

## 5. Other sources

### 5.1 SEU-ALLEN / BICCN fMOST (Peng et al. 2021)

VERIFIED: Peng H, Xie P, Liu L, ... Zeng H. Nature 598:174-181 (2021), doi:10.1038/s41586-021-03941-1
(article CC BY 4.0): 1,741 full reconstructions, extended to 1,876 (Liu Y et al., Nat Commun 15 (2024),
doi:10.1038/s41467-024-54745-6). Brain Image Library doi:10.35077/g.25: each cell has an "original" and a
"registered" SWC in the Allen CCF; the library's stated licence is "Creative Commons Attribution-ShareAlike
4.0 International". In NeuroXiv's listing the 1,876 come from 39 brains, none of whose IDs is a portal or
Jiao animal. **Hypothalamus: none by the source's own soma label.** 27 neurons have a NeuroXiv soma voxel in
ZI 18, LHA 4, FF 3 or HY 2 (7 brains), but the source labels them thalamic (VPM, VPL and their parvicellular
parts, VM, MG, RT, VAL, PoT): border cells, not hypothalamic neurons. **Verdict: training only.**

### 5.2 NeuroXiv

VERIFIED: Jiang S, Wang L, Yun Z, Chen H, ... Peng H. Nat Methods 22:1195-1198 (2025),
doi:10.1038/s41592-025-02687-2; https://neuroxiv.org, files at https://download.neuroxiv.org. Its metadata
listing (19,401,529 bytes, 18 metadata columns, no measurements) has 175,149 rows, all tagged CCFv3:
SEU-ALLEN 157,619 (1,876 with axon, 155,743 local dendrites), ION 16,380 (PFC and hippocampus, 401 brains),
MouseLight 1,150. IDs look like `ION_full_17099_002_CCFv3`, `MouseLight_full_AA0001_CCFv3`. Soma position is
given in 25 µm units with a Left / Right column. Hypothalamic somata among neurons with an axon: MouseLight
55 by coordinates (69 by source label, 70 either), SEU-ALLEN 27 (see 5.1), ION 0.

VERIFIED: the file list now also has **`ION_Hypo_full_SWC_CCFv3.zip`, dated 2025-12-31**, which the README
and the listing do not describe. INFERRED: a re-registered copy of the Jiao set. Terms (README): "available
under the CC-BY NC license"; among the prohibited uses: "Distribute, sublicense, or sell the Dataset or any
derivative works". Zip sizes UNVERIFIED (not requested). UNVERIFIED: how NeuroXiv resamples or re-registers files.
**Verdict: a convenient single-frame route to MouseLight, SEU-ALLEN and ION, not a new source**; its
MouseLight copy is older and smaller than the source (1,150 against 1,653).

### 5.3 NeuroMorpho.Org

VERIFIED (public search service, release 8.7 of 2026-10-29): 1,037 mouse records with brain region
hypothalamus, in 13 archives. Only 338 have an axon domain, and only the 66 MouseLight ones are whole-brain
reconstructions with a complete axon. The rest are slice or thick-section work with local or cut axons (Chen_Dong 217 PVH
oxytocin, Rancillac 66 preoptic, Chee 34 MCH, Cai 11 POMC, Burdakov 33, Stagkourakis_Broberger 10) or have no
axon at all (Xiong 335 microglia, Dias_Lima 185, Peca 46, Wang_Zhou 30); Ruiyu and Peng add 2 each. There is
no Jiao archive. Licence: "Creative Commons Attribution 4.0 International". VERIFIED (FAQ): standardised
files are moved so the soma sits at the origin and rotated to principal axes, so they are no longer in an
atlas frame; programmatic access to the original files is "Temporarily not accessible". **Verdict: neither.**

### 5.4 Hypothalamic reconstructions in papers, with no public files found

- Lin W et al., PNAS 120 (2023), doi:10.1073/pnas.2216231120. VERIFIED: 60 histamine neurons of
  the tuberomammillary nucleus (HDC-CreERT2), 52 with long projections, fMOST, registered to CCF v3.0, saved
  as SWC; availability statement: all data are in the article and supplement. No repository found.
- Senol E et al., Nat Commun (2025), doi:10.1038/s41467-025-60585-9. VERIFIED: 186 tuberal Sst neurons, fMOST
  and FNT, imaging contributed by Jiao's group; only source data provided. The Jiao index has 133 Sst
  neurons labelled TU (11 animals). UNVERIFIED whether the 186 include them.
- Zhang B et al., Neuron 109:331-346 (2021), doi:10.1016/j.neuron.2020.10.032: magnocellular oxytocin system.
  Abstract only; no single-neuron file deposit found (UNVERIFIED beyond the abstract).

**Verdict for all three: neither**, unless the authors share files. Search limits: Europe PMC and CrossRef
keyword searches, the BSDC catalogue and NeuroMorpho; a dataset described only in a supplement could be
missed.

## 6. Conventions side by side

| | Jiao (portal) | Ding `spcd` | Li `pvh_oxt` | MouseLight | NeuroXiv copies |
|---|---|---|---|---|---|
| frame | CCFv3 2017, 10 µm grid (VERIFIED earlier) | 2017-style labels (VERIFIED); frame INFERRED same as Jiao | CCFv3 (VERIFIED); website copy and paper differ slightly | CCFv3 and a legacy CCFv2.5 (VERIFIED) | "CCFv3" tag (VERIFIED) |
| units, axes | µm; AP, DV, ML | INFERRED same | INFERRED same | µm; x = ML, y = DV, z = AP (INFERRED) | soma in 25 µm units; files UNVERIFIED |
| hemisphere | true side; midline ML 5,700 | both in index; INFERRED true side | both in index; INFERRED true side | true side (VERIFIED from somata) | Left / Right column |
| registration | UNVERIFIED (Methods paywalled) | UNVERIFIED; protocol: ANTs on region masks | Amira features, deformation field (VERIFIED) | affine + landmark splines (VERIFIED) | UNVERIFIED |
| node types | 0, 1, 2, 3 documented; axon as 2 or as 0; 10 to 12 undocumented | documented 0, 1, 2, 3; files UNVERIFIED | UNVERIFIED; index says axon only | 0 to 6 in the service; export UNVERIFIED | UNVERIFIED |
| full axon | axon in all 7,180; dendrite in 6,871 | axon in all 266; dendrite in 264 | axon only | axon and dendrite, two annotators | as source |
| label | 16 lines | none | Oxt | none public | none |
| sex | male (title) | UNVERIFIED | female (VERIFIED) | female (VERIFIED) | - |
| file | `<animal>_<n>.swc` | `<animal>_<n>.swc` | `<animal>_<n>_L.swc` | `<id>.swc` or JSON | zip per source |

## 7. What a fetch would be (nothing was fetched)

| step | files | bytes | route | login | terms that bind |
|---|---|---|---|---|---|
| A. Ding, the 266 | 266 single files | UNVERIFIED; about 260 MB by analogy with the 23 Jiao portal files already on disk (mean 992,448 bytes) | portal address in 3.1 | UNVERIFIED (INFERRED none) | portal copyright line; BSDC CC BY-NC-ND 4.0 |
| A'. Ding, all relay neurons | 1 zip (1,464 neurons) | 1,086,478 KB as displayed | BSDC | **account needed** | CC BY-NC-ND 4.0 |
| B. MouseLight, the 21 CC BY records | 21 SWC | **16,540,071** (VERIFIED, figshare metadata) | figshare record of each DOI | none | CC BY 4.0; cite each DOI |
| B'. MouseLight, the other 76 | 76 SWC | UNVERIFIED; about 58 MB at the mean of the 25 sizes seen (767,016 bytes) | NeuronBrowser export | none | CC BY-NC 4.0; cite each DOI |
| C. Li, mouse 18696 | 115 single files, or `/PVH.rar` (all 264) | single files UNVERIFIED; archive 39,628 KB | portal address in 3.2, or BSDC | portal: none (INFERRED, see 3.2); BSDC **account needed** | portal copyright line; BSDC CC BY 4.0 |
| D. Training sets, only if wanted | CeA 488; `spcd` rest 1,983; cortex 18,621; hippocampus 10,100 | BSDC zips: 0.45, 2.4, 17.9 and 4.7 GB as displayed | BSDC, or portal single files | BSDC account | CC BY-NC-ND 4.0 |

VERIFIED from the data centre's page code, not by clicking: both download handlers first call a helper that
returns true only when a stored access token exists and otherwise sends the browser to the login page. So
**BSDC downloads need a registered account even for "Open Access" records.** Its terms of use say data may
be browsed, downloaded or used subject to the owner's licence, forbid excessive load and getting around
access limits, and say the IP address of each view and download is recorded.

**What the user must do.**
1. Say which of A, B, B', C are wanted. A, B and C are the minimum for an external test; D is optional.
2. For any BSDC archive: create a BSDC account and log in personally (not done here, and not something this
   audit may do). Not needed if single files from the portal are enough.
3. Accept the licence limits: keep raw files local; no redistribution of converted or re-framed files under
   the ND and NC terms; cite each MouseLight DOI; credit BSDC and the authors.
4. Optional, by e-mail to the authors: the licence of portal-served copies; sex and registration for the Ding
   set; whether the Li and Jiao cells in 195037 and 195038 are distinct; which MouseLight sample is Sim1-Cre.

## 8. What each could be used for, and the split

**Order of value as an independent test.**
1. **Ding, 266 neurons / 5 mice.** Same frame and file family as Jiao (INFERRED), so almost no conversion.
   Weakness: same institute, pipeline and atlas registration habits, so it tests new animals, not a new
   laboratory. No label: it can score a location-only model and a model run with its label input withheld,
   not the label's contribution.
2. **MouseLight, 97 neurons / 22 mice.** The only cross-laboratory check. Needs axis reordering, a node-type
   map and a frame check on the first file. Registration differences of the order of 100 µm are expected and
   are part of what it tests.
3. **Li, mouse 18696, 115 neurons.** One animal, oxytocin, axon only: a sanity check for one label.

**Proposed handling in a whole-animal split.**
- All 28 external animals are test-only by default: never in training, never donors. If the supervisor wants
  some for tuning, split by animal: for example Ding 210661 and 210726 for validation, the other three for
  test; MouseLight by sample.
- 195037 and 195038 are Jiao animals. Li's 149 neurons from them go wherever the Jiao split puts those two
  animals, on both sides, or are dropped. 195038 is already an inspected animal.
- `lha` adds nothing and must not be counted twice.
- **Query safety.** A regional state built from Jiao animals contains no external animal, so external
  queries are safe by construction, provided no external geometry enters training. If section 7 D sets are
  ever used for training, the rule is the same as inside Jiao: no state, and no decoder donor pool, may
  contain a neuron of the animal being queried; `spcd` relay neurons from the five hypothalamic animals
  (527 neurons in all, 266 of them hypothalamic) must then be held out together.
- Before any scoring: fix the eligibility rule (HY-subtree label; axon present; complete), the conversion per
  source and the metric, and check the frame on files with no hypothalamic outcome attached to the decision.

## 9. New facts for earlier notes, and open points

New or corrected (all VERIFIED unless marked):
- BSDC downloads need an account (read from the site's page code, not tried). `reports/swc_access.md` and
  `reports/d109_jiao_conventions.md` left this open.
- The Li 2024 archive is CC BY 4.0 at BSDC. `reports/conversion_map.md` says "none stated".
- A second Jiao record exists at BSDC: "Projectome of hypothalamic peptidergic neurons in male mice",
  2025-07-09, one file `/HY_swc_all.zip` (3,088,813 KB, 2.88 GB), CC BY-NC-ND 4.0,
  https://cstr.cn/33145.11.BSDC.1752068512.1942934383845957634; and NeuroXiv has `ION_Hypo_full_SWC_CCFv3.zip`.
- The portal's LHA project equals Jiao's LHA subset exactly.
- Sex of 195037 and 195038 conflicts between Li 2024 (female) and Jiao 2025 (male).
- MouseLight now lists 1,653 neurons (NeuroXiv's README still says 1,200), and its licence is not uniform.

Open (UNVERIFIED):
1. Frame, axes, node types and file sizes of the Ding files; whether their portal address needs a login.
2. Registration method and sex for the Ding set (Methods paywalled).
3. Whether Li and Jiao cells in the two shared mice are distinct.
4. Which of the 97 MouseLight neurons come from the Sim1-Cre brain; the frame of the figshare SWC copies.
5. Why the MouseLight query returned 1,641 of 1,653.
6. Sizes and processing of the NeuroXiv archives; whether `ION_Hypo` is the full Jiao set.
7. What licence applies to files taken from the portal rather than from BSDC.
8. Which anatomical side is "Left" in each frame.

## 10. Sources relied on

Copies are only in the scratch directory
`/private/tmp/claude-501/-Users-andresarriaga-Documents-SYN4RCH/6fc8191b-130f-4b76-b6a8-69340161f99e/scratchpad/bench/axon_audit_other_datasets/web/`
(every request, including failures, is in `.../axon_audit_other_datasets/fetch_log.tsv`). The cached portal
viewer code read for the download addresses is `.../scratchpad/li/app.js` (6,482,711 bytes, SHA-256
`58502a6d7075576a50b6e99a501e8e7daed63d696560e52495aee64dc4df8099`).

| local copy | bytes | SHA-256 | fetched (UTC) | URL |
|---|---|---|---|---|
| `portal_projects.json` | 7,821 | `7154a66174ea955c05766310ea07f599260d48a87cc68eac1ca4177788680b82` | 2026-11-11T01:10:49Z | `https://mouse.digital-brain.cn/projectome/2/srv/api/v1/public/projects` |
| `portal_pfc.neuron.info.json` | 1,655,767 | `29dc94374f9eb87bb6722c247af0712d54c906229feb11ebc3793083231a9469` | 2026-11-11T01:11:41Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/mouse.neuron.info.json` |
| `portal_hipp.neuron.info.json` | 1,823,150 | `472418549780fd2fd7c6eb6164a93ac541babc05fdf004921523bb20f0c209f8` | 2026-11-11T01:11:54Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/hipp/mouse.neuron.info.json` |
| `portal_cea.neuron.info.json` | 190,309 | `78282e2645d266f65d62a815d539990c7483e65dfd7534163a0c42b3a022bd96` | 2026-11-11T01:12:04Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/cea/mouse.neuron.info.json` |
| `portal_pvh_oxt.neuron.info.json` | 52,388 | `61203bb4056330349c260cf0668035a7a9b043c42561113fd8f8b89785d55b86` | 2026-11-11T01:12:09Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/pvh_oxt/mouse.neuron.info.json` |
| `portal_lha.neuron.info.json` | 253,755 | `186d66593a708e3a11828629e22b9c0b8d9c639ddb6f5c5b83abc1b55b668e27` | 2026-11-11T01:12:13Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/lha/mouse.neuron.info.json` |
| `portal_spcd.neuron.info.json` | 548,951 | `fd7c63ed2f9c736747146eb8dcdc94f8fc0bc5e9d274a43b7da5ec13e9854927` | 2026-11-11T01:12:18Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/spcd/mouse.neuron.info.json` |
| `portal_cortex.neuron.info.json` | 7,499,161 | `229251317a07f8853bf0c84e1a4721d25be5af0ea0604b600f9ee17cca9a93d5` | 2026-11-11T01:12:24Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/cortex/mouse.neuron.info.json` |
| `portal_rbm.neuron.info.json` | 1,197,653 | `7317bba899ce79de078063126eb8339a9eefc2c8dd72a3183af3cb00468bfd5b` | 2026-11-11T01:12:45Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/rbm/mouse.neuron.info.json` |
| `portal_sample_PFC.json` | 137,231 | `7f7681da2bba09cc53dc3e1b15a65d3b7c6511726f9fe1683f05a40aee2c0fce` | 2026-11-11T01:12:52Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/sampleInfo/PFC.json` |
| `portal_sample_HIPP.json` | 210,418 | `db9a6fdecfa5f0627ecbb8a0feafb26c61a06b703d482b71d93915d72c84786c` | 2026-11-11T01:12:56Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/sampleInfo/HIPP.json` |
| `portal_sample_cea.json` | 21,892 | `4416b3f23124b84c0b584bb3972fe6978aceb90bc20938dc02f28b1c65561ae6` | 2026-11-11T01:13:01Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/cea/sampleInfo/PVH_Oxt.json` |
| `portal_sample_pvh_oxt.json` | 7,418 | `a15fec89f5a8a85062c4e4967e53346509eddc9f15b8fe0e78e12999521562ed` | 2026-11-11T01:13:05Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/pvh_oxt/sampleInfo/PVH_Oxt.json` |
| `portal_sample_spcd.json` | 54,715 | `b9819742cd11c236867fac4929f77c253b2fad09ccf08264cde1ec799321f263` | 2026-11-11T01:13:08Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/spcd/SPCD.json` |
| `portal_sample_cortex.json` | 416,406 | `f7d7bba85f91f815ebca925f06f34025599acae6589a224fbccac8155ac5a1fb` | 2026-11-11T01:13:12Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/cortex/sample.info.json` |
| `portal_sample_rbm.json` | 93,274 | `87e0fbd5169f4ce77349c6e71113efa95c8509d95f841de335f419d520a71c3e` | 2026-11-11T01:13:18Z | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/rbm/RBM.json` |
| `portal_chunk_1.js` (Download page: dataset table, file specification, copyright line) | 148,908 | `dd54ab673039290231a4c2968bf49b907275ff052fc9605226eb3844300818d3` | 2026-11-11T01:16:50Z | `https://mouse.digital-brain.cn/projectome/js/1.js` |
| `bsdc_app.js` (login helper) | 64,914 | `1c6a5748223d6ab26d8a736834bdd7d8dfaa841b4c5ccdddc034bcbae7febfad` | 2026-11-11T01:17:49Z | `https://www.braindatacenter.cn/datacenter/web/js/app.833b21ae.js` |
| `bsdc_chunk_details_a.js` (dataset page: download handlers) | 142,649 | `d8d7f3b41a9d96e4b7f8fb7dbe3d9d4e1a28ab4ad91931e8346dfe5b22bc7fc3` | 2026-11-11T01:43:23Z | `https://www.braindatacenter.cn/datacenter/web/js/chunk-ea9c37ce.b0c5bbc6.js` |
| `bsdc_records_rendered_2026-11-10.txt` | 15,899 | `0ba7147164a647d2d1f924e8a65a2f024966af2f1411cbb7369138211da94b4a` | 2026-11-11 | transcription of BSDC record, terms and catalogue pages read in a browser; URLs inside |
| `cell_pages_rendered_2026-11-10.txt` | 3,395 | `a23d56812ffed1bba764f594d74d64cf418aae12ba027a1d40c6d8cde04a5198` | 2026-11-11 | paraphrased notes of the Li 2024 and Ding 2025 publisher pages read in a browser; URLs inside |
| `doi_handle_hipp.json` | 485 | `d956ff2e9523622d7641841dab46303d03df8d00215e75912f6d78d6019aa216` | 2026-11-11T01:18:08Z | `https://doi.org/api/handles/10.12412/BSDC.1667278800.20001` |
| `doi_handle_cortex.json` | 485 | `fb652f9df2a6ec3816b63135ac2bbedecee2757c467557da7d7f7bc4811de791` | 2026-11-11T01:18:10Z | `https://doi.org/api/handles/10.12412/BSDC.1747279998.20001` |
| `crossref_ding2025.json` | 32,950 | `eef020a947cdfddb95cb9b4fff1b973327c0056a0156d0ecdfc7a025e2b10dc2` | 2026-11-11T01:24:25Z | `https://api.crossref.org/works/10.1016/j.neuron.2025.03.007` |
| `crossref_dou2024.json` | 32,501 | `8b6f66a7321a168b663da7d1bb93fd1a0dff7d456a097faec481d12356d2fd2b` | 2026-11-11T01:24:30Z | `https://api.crossref.org/works/10.1093/nsr/nwae195` |
| `crossref_li2024.json` | 37,014 | `704bbe359d908bfe18418310b1908713cc25a25c1c9038b9f495484a4bd8b111` | 2026-11-11T01:24:34Z | `https://api.crossref.org/works/10.1016/j.neuron.2023.12.022` |
| `crossref_qiu2024.json` | 26,441 | `40108fb0bebd5098bbfe3c2f2b35b725d3af9d5e9d09d39b0a1bb44cd70d6c7e` | 2026-11-11T01:24:38Z | `https://api.crossref.org/works/10.1126/science.adj9198` |
| `crossref_gao2022.json` | 27,437 | `05eef65c7b0189bfe9d1c62c94e0a452025b0cd7f4a7fda565c718f98bcd012a` | 2026-11-11T01:24:42Z | `https://api.crossref.org/works/10.1038/s41593-022-01041-5` |
| `crossref_gao2023.json` | 23,219 | `6770bc3a9586113466df54b4f1cbbe27474e65375453595ed5bd0fc1bbb40a51` | 2026-11-11T01:24:47Z | `https://api.crossref.org/works/10.1038/s41593-023-01339-y` |
| `crossref_jiao2025.json` | 36,556 | `c4a521c75a2a384fb014c09a4f98062bee4322b4b172a91871ea77a3626b8719` | 2026-11-11T01:24:51Z | `https://api.crossref.org/works/10.1038/s41593-025-01919-0` |
| `crossref_peng2021.json` | 40,714 | `f7c76c21c1d11b6756bbef9e7f87a34b3b7d5112bf5bfb2bb4858dc2f7e72e61` | 2026-11-11T01:35:34Z | `https://api.crossref.org/works/10.1038/s41586-021-03941-1` |
| `crossref_q_neuroxiv.json` | 6,005 | `9713f1772f26444af4241080cbf3cdf6aff097816b644f71b92f4f34f2c528b0` | 2026-11-11T01:35:28Z | CrossRef title search for the NeuroXiv paper |
| `crossref_q_liu2024.json` | 18,180 | `3e3c44a99c148f7b8de57ad58f5082793a1f4fc6652aaea5e69f584201bd1bfc` | 2026-11-11T01:35:39Z | CrossRef title search for Liu 2024 |
| `crossref_q_hns.json` | 6,103 | `5821c58c9520eaa2ca8caaad216d933e4baf543330ef29cf1041bf4bbcd8182e` | 2026-11-11T01:35:45Z | CrossRef title search for Zhang 2021 |
| `crossref_q_hist.json` | 8,598 | `27b1e6b58de2cf6105acc97336ff9eda9e91854da015c8622fd4a40458c789cc` | 2026-11-11T01:35:51Z | CrossRef title search for Lin 2023 |
| `epmc_ding2025.json` | 17,747 | `0d65b61995db4f49437590e867e51bc2d514e9f179d158cbd878daa966a35608` | 2026-11-11T01:25:12Z | Europe PMC record of doi:10.1016/j.neuron.2025.03.007 |
| `epmc_s1.json` | 21,475 | `c7e16ef84f646c355d5e3ca8c0e3a818773625daef010409f83c23fa9e74cba8` | 2026-11-11T01:40:37Z | Europe PMC search: single-neuron projectome, hypothalamus |
| `pmc12800304_efetch.xml` (Ji 2025) | 244,544 | `0beb5be7a7b05ea32ae5c55644dc268b2c674c39f59e7564c31a1b3703d0c7d7` | 2026-11-11T01:41:02Z | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=12800304&rettype=xml` |
| `pmc12067299_efetch.xml` (Ren 2025 protocol) | 60,472 | `c2b1a62bbafca950ac047bf63eb605e146b2ee8f795fc4d833eaa2dc8fe27d0b` | 2026-11-11T01:41:03Z | `...efetch.fcgi?db=pmc&id=12067299&rettype=xml` |
| `pmc12216382_efetch.xml` (Senol 2025) | 190,426 | `e52480911456bd98e0ca0dc0259b2d364dff8aaf44ced463089f74fc0cf2b141` | 2026-11-11T01:41:05Z | `...efetch.fcgi?db=pmc&id=12216382&rettype=xml` |
| `pmc10083611_efetch.xml` (Lin 2023) | 176,492 | `1a5040744cc5cafa611f3538b59bd5cedf3738596bd50359d42bac7027ef1001` | 2026-11-11T01:36:17Z | `...efetch.fcgi?db=pmc&id=10083611&rettype=xml` |
| `pmc6754285_efetch.xml` (Winnubst 2019) | 195,667 | `bfd98d8e0797d0e842625707818db5d20adce420a91de90f6858cb61606b6a35` | 2026-11-11T01:32:56Z | `...efetch.fcgi?db=pmc&id=6754285&rettype=xml` |
| `ml_bundle.js` (NeuronBrowser page code) | 2,039,911 | `903e1dcae93c6a4be22ce97b20f63d8e4a99cfa319cbf3f6be518a76177651c9` | 2026-11-11T01:26:35Z | `https://ml-neuronbrowser.janelia.org/bundle.js` |
| `ml_constants.json` | 279,086 | `b208e58b94a78dc5f80f2b87a8ddb4c1b916758c75f297b289e13188a0514757` | 2026-11-11T01:27:16Z | query to `https://ml-neuronbrowser.janelia.org/graphql`; body saved as `ml_q_constants.json` |
| `ml_all_ccf30.json` (identifiers and one soma record per tracing; no tracing nodes) | 1,389,723 | `6fae0d984a6f88de8348e5eb726311cecbaa2f69c990a9f97c80156fe71c0f6f` | 2026-11-11T01:27:49Z | same service; body saved as `ml_q_all.json` |
| `figshare_ml/<neuron>.json` (108 records) | 511,998 | `ea25f19de0a225346f3bffa26c6f256089a439de4a7412750a879dd8d3d43e6a` (over the sorted list of name and file hash) | 2026-11-11 | `https://api.figshare.com/v2/articles/<id>` |
| `janelia_ml_neuronbrowser.html` (terms) | 62,583 | `82c128ae4622f6647c6b410865d696b2f72175285ebe01cc2b499bcb1a175c72` | 2026-11-11T01:31:57Z | `https://www.janelia.org/project-team/mouselight/neuronbrowser` |
| `janelia_ml_resources.html` | 56,322 | `86a4dae638887340681af861830059bbca54f8060744a14972d79eb4d973027d` | 2026-11-11T01:32:06Z | `https://www.janelia.org/project-team/mouselight/resources` |
| `aws_registry_mouselight.html` | 8,014 | `445a862a22b1fb2b8b7f5c0dc2d1ed82fb524fdfce0e62fb70de4efe793c531e` | 2026-11-11T01:32:14Z | `https://registry.opendata.aws/janelia-mouselight/` |
| `neuroxiv_download.html` | 6,859 | `424897eed6ccc95ab1dc252bb1a8e318e05fdaa5ab8606613c44f9a092bd32aa` | 2026-11-11T01:33:47Z | `https://download.neuroxiv.org/` |
| `neuroxiv_list_files.json` | 1,594 | `05391fdbe16f3b712fd92464d40225a1c6b165007ca736d6f6d11be1699bab39` | 2026-11-11T01:33:59Z | `https://d36ajqhpoeuszk.cloudfront.net/list-files` |
| `neuroxiv_README.md` | 9,393 | `a0b8c830ecb43c6c467eb2fd4549f06647dbdfa58538f44e402cf113434dcc43` | 2026-11-11T01:34:05Z | `https://d36ajqhpoeuszk.cloudfront.net/README.md` |
| `neuroxiv_metadata.csv` | 19,401,529 | `67616c1157ebcc6e61c4bbaf2920977a88115f1eca4c374dd063b1cc888445e2` | 2026-11-11T01:34:38Z | `https://d36ajqhpoeuszk.cloudfront.net/NeuroXiv_metadata.csv` |
| `nmo_mouse_hy_p0.json`, `_p1`, `_p2` | 891,372; 808,132; 56,507 | `3d659a2112c7948210fa0baf92c807b0e9ec531ef029b8321adc54b1f041a4e9`; `1e44c38212b73eed34b10730edc166aa75720478e337d92a752893f774abf617`; `828b4b67857171b867f7610c0f083a24fdc7148350ed32e39725c7da5fea5298` | 2026-11-11T01:36:43Z to 01:36:53Z | `https://neuromorpho.org/api/neuron/select?q=species:mouse&fq=brain_region:hypothalamus&size=500&page=0` (1, 2) |
| `nmo_useterm.html` | 16,466 | `1217b07d558eb8ae86d50e0ef07547c71ee9bcba97ac9609621e6cea00d04423` | 2026-11-11T01:37:37Z | `https://neuromorpho.org/useterm.jsp` |
| `nmo_faq.html` | 65,898 | `807c5e7cde8c8438e9a98b6e4277dd9458ca8ecc3dd115e5706a636b02ea3d81` | 2026-11-11T01:37:46Z | `https://neuromorpho.org/myfaq.jsp` |
| `bil_g25.html` | 138,120 | `5ddff81354c11573efa10ed78b38648abd0c18c89f7622e3aa30899086f98713` | 2026-11-11T01:37:55Z | `https://doi.org/10.35077/g.25` |
| `bil_citation.html` | 10,125 | `8b977541fb5a88fc0abd2b6967a40d984a9708ba21ab885d07b140e23515c26e` | 2026-11-11T01:38:16Z | `https://www.brainimagelibrary.org/citation.html` |

Read in a browser, so not hashable as bytes (transcriptions hashed above): the BSDC record pages
`https://www.braindatacenter.cn/datacenter/web/#/dataSet/details?id=<id>` for ids 1879450714607951874 (Ding),
1742375820989538306 (Li), 1722177494323273729 (Li, figures), 1942934383845957634 (Jiao, 2025),
1922847323835588609 (whole cortex), 1755143583550119938 (hippocampus), 1585453592638009346 and
1681926304704638978 (PFC), 1925514407162388481 and 2051907638291587073 (CeA), 2000812574081802241 (nucleus
accumbens), 1923208705126948866 (PT neurons); its terms page `#/termsForUsage`; its catalogue `#/dataSet`
with the keywords "projectome" and "hypothal"; and the publisher pages
`https://www.cell.com/neuron/fulltext/S0896-6273(23)01010-3` and `.../S0896-6273(25)00179-5`.

Local inputs used: `data/raw/jiao2025/swc/portal_index/mouse.neuron.info.json` (allowed fields only),
`data/raw/jiao2025/swc/MANIFEST.tsv` (IDs and byte counts only), `data/raw/atlas/structure_graph_1.json`,
and the reports `swc_access.md`, `swc_validation.md`, `d109_jiao_conventions.md`, `conversion_map.md`,
`paired_data_search.md`, `d109_protocol.md` (section 16).

## 11. How to reproduce

Scripts are in `.../axon_audit_other_datasets/scripts/` and run as `python3 -I <script>`:
`portal_inventory.py` and `portal_details.py` (portal indexes, overlap with Jiao, the 37 and the 16),
`cross_ids.py` (animal IDs across sources), `ml_inventory.py`, `ml_figshare_summary.py`, `ml_sizes.py`
(MouseLight), `neuroxiv_inventory.py`, `nmo_inventory.py`, `id_lists.py` (the ID lists above, also saved as
`.../axon_audit_other_datasets/candidate_external_ids.json`), `sources_table.py` (re-checks every byte count
and hash in section 10 against the log). `fetch_doc.sh` is the only fetcher; it refuses any address ending
in a geometry or archive extension.
