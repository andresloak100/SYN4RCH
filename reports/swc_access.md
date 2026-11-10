# How to get individual Jiao 2025 SWC files publicly

Checked 2026-11-10. Read-only: no account, no login, no terms dialog accepted, nothing downloaded in bulk.
Source paper: Jiao et al. 2025, Nat Neurosci, doi:10.1038/s41593-025-01919-0 (7,180 neurons, 166 brains).

## Short answer

- **RESOLVED: single SWC files are public.** The portal serves each neuron as a plain-text file at
  `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/hy/hy_neuron_download/<brain>_<n>.swc`.
  No login, no token, no terms dialog. This confirms `conversion_map.md` ("single SWC files without login").
  The scout's "login probable" applies to the BSDC record, not to the portal.
- **RESOLVED: the portal frame and the Zenodo frame are different.** Portal files keep the true hemisphere.
  The authors' own code mirrors every right-hemisphere neuron to the left (ML' = 11400 - ML) before making
  the Zenodo files and the paper tables. Keep both frames apart and label every file with its frame.
- **Licences differ by source** (see Terms). Zenodo: CC BY 4.0. BSDC: now shows CC BY-NC-ND 4.0.
  Portal: copyright line only.
- **Recommended sample: 47 portal SWCs, about 52 MB.** List at the end.

## 1. Portal access

### Evidence

| what | finding |
|---|---|
| Viewer | `https://mouse.digital-brain.cn/projectome/hy` (title "Mouse Projectome Atlas - Digital Brain CEBSIT(ION)"); Vue app in `/projectome/js/app.js` |
| Download code | `DataOperationPanel.vue` (in app.js), function `downloadSWC`: base `https://mouse.digital-brain.cn/projectome/2/srv/` + `/info/mouse/` + `hy/hy_neuron_download/<file>.swc`, fetched with a plain `axios.get`. The same route serves project `lha`. |
| Login gate | None on the download path. `onDownloadNeurons` refuses only selections of 10,000 neurons or more ("Please select 10000 or fewer neurons to download"). It then fetches in batches of 100, 4 at a time, and zips them in the browser. Login appears only for saving groups and animations ("Please login to see your data"). A global axios interceptor sends `Authorization: Bearer <localStorage access_token>` (null when logged out); the server answers 200 anyway. |
| Index | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/hy/mouse.neuron.info.json`, HTTP 200, `application/json`, 1,980,365 bytes, Last-Modified Thu, 07 Dec 2023 09:47:56 GMT, SHA-256 `fdf1fa2eff099bc755e5218dc6437e772f12fe0a934c01a4072d846d9b3f4da3` (same bytes as the session cache `scratchpad/li/hy.neuron.info.json`; `hy.neuron.info.json` on the server returns 404, so that is only the local name). Listed in `/projectome/2/srv/info/projects.json` under "Mouse Hypothalamus, Jiao et al. 2025". |
| Index content | `neuron_data`: 7,180 rows, one per file, no duplicate files, 166 brains. Fields: `file`, `acronym`/`name` (soma nucleus), `type_name` (peptide line, 16 values), `hemisphere` (Left 3,071 / Right 4,109), `reconstruction_type` (Axon_and_dendrite 6,871 / Axon_only 309), `type`/`class` (projectome type id). |
| SWC fetch (from the user's probe, `/private/tmp/synarch_jiao_geometry_probe/probe.json`) | e.g. `192710_002.swc`: HTTP 200, `nginx/1.16.1`, `text/plain; charset=utf-8`, Content-Length 1,090,525, `Accept-Ranges: bytes`, Last-Modified Tue, 27 Aug 2024 03:19:04 GMT, SHA-256 `8873e27f89d858c36483467bfd4dd14dc5a9318072d7148d953484def72f5255`. No new SWC was fetched for this report. |
| robots.txt | None (the path returns the site's HTML page). |
| Portal "Data Download" page | `https://mouse.digital-brain.cn/projectome/download`. For "Single-neuron datasets for mouse hypothalamus" (2024-06-12) it links only to the BSDC record. |

### Probe files (read-only; hashes recomputed here)

| file | frame | bytes | SHA-256 |
|---|---|---|---|
| 192710_002.swc (portal) | portal, true hemisphere (Left) | 1,090,525 | 8873e27f89d858c36483467bfd4dd14dc5a9318072d7148d953484def72f5255 |
| 195038_041.swc (portal) | portal, true hemisphere (Left) | 530,638 | 8d6f299313ed1a015329bb21706274dd30ebc6d142eaf224042ef5b7ec3b55d5 |
| 201232_056.swc (portal) | portal, true hemisphere (Right) | 380,151 | 8a6bb67e6acb61937f96957204a06f298e9634b200d77e635b0b13781bdd78b3 |
| zenodo_195038_041_Oxt_PVH.swc | Zenodo, mirrored-left (Left, so unchanged) | 540,214 | 712527f83796a6a7ddc0101d4680f1feb714f6b5e3e69b60032acd23ec37e114 |
| zenodo_201232_056_Tac1_MMm.swc | Zenodo, mirrored-left (was Right) | 387,204 | d74c0b02bb1a70871f7d9aef6f39742f1c7b79604d2b9a96331b1fb1980b262b |

The two Oxt files have the same values. They differ in bytes only because the line endings differ (the probe found them equal once newlines were normalised).

### BSDC record

`https://www.braindatacenter.cn/datacenter/web/#/dataSet/details?id=1800731789745954817`
(DOI 10.12412/BSDC.1667278800.20003). It holds 166 files, one zip per brain, 6.20 GB in all. The listed zips are 9–122 MB
each (for example `/192710.zip` 23,594 KB and `/201224.zip` 121,709 KB). They are **not** 6.2–6.7 GB each, as the scout reported.
The page has Download / "Copy Download URL" / "Batch Download" buttons. None was clicked.
UNRESOLVED: whether those downloads need a login. They are not needed, because the portal serves the same neurons one by one.
The authors' notebook reads `<brain>/swc_allen_space/<n>.swc` from these zips, so a zip holds the Allen-space SWCs in the original hemisphere (see section 3).

## 2. Terms (exact text)

**Portal** (footer of `https://mouse.digital-brain.cn/projectome/download`; the viewer page shows no terms):
"Copyright © 2023 Digital Brain. All Rights Reserved." No licence and no terms of use. No dialog appears before downloading.
What this allows: nothing beyond viewing and the download the site offers. Treat the files as all rights reserved:
local analysis and private display with attribution; no redistribution of the raw SWCs.

**BSDC record page**, section "Copyright & License":
"data users should clearly indicate the source of the data and the author of the data in the research results generated by using the data … For re-posting (second or multiple releases) data, the author must also indicate the source of the original data."
Suggested acknowledgement: "The dataset is provided by Brain Science Data Center, Chinese Academy of Sciences (https://braindatacenter.cn/)."
Citation: Xiaofei Wang (2024), Single-neuron datasets for mouse hypothalamus, BSDC, https://doi.org/10.12412/BSDC.1667278800.20003.
**Details panel: "License Agreement CC BY-NC-ND 4.0".** This contradicts `conversion_map.md` line 38 ("BSDC licence field is empty … 'CC BY-NC-ND 4.0' was wrong"). That line should be corrected: the field shows CC BY-NC-ND 4.0 today.
Under CC BY-NC-ND: viewing, non-commercial analysis and sharing of the unmodified data with attribution are allowed. Redistributing adapted versions (for example converted or re-framed SWCs) is not. Figures of the analysis are generally fine with attribution, but the NC term rules out commercial use.
UNRESOLVED: the BSDC "Terms of use" footer page did not open from the record page, so its text is not quoted.

**Zenodo 14588094** (`https://zenodo.org/records/14588094`): "License: Creative Commons Attribution 4.0 International". The API gives `"license": {"id": "cc-by-4.0"}`.
CC BY 4.0 allows viewing, analysis and redistribution of the files and derived figures, including modified versions, with attribution.
This covers only the 36 SWCs in the zip (`Supplementary Fig/Supplementary Fig4/swc/`). 18 are Jiao neurons; 18 are MouseLight (`AA*.swc`). The MouseLight files are not Jiao neurons: keep them out of Jiao counts and give them their own MouseLight source and licence.

**Practical rule:** keep the raw portal files local (data/raw, not committed, not republished). Show them only on private pages with attribution to Jiao et al. 2025, CEBSIT/Digital Brain and BSDC. Only the 18 Zenodo Jiao SWCs may be redistributed (CC BY 4.0).

## 3. Coordinate conventions

| item | value | source |
|---|---|---|
| Columns | `id type x y z radius parent`, 7 per line, space-separated, no header lines | portal Download page "swc format specification"; probe files have 7 columns on every row and no `#` lines |
| Axis order | x = AP, y = DV, z = ML | authors' `two_D_plots.py` (Zenodo): `dfneuron.columns = ["ID","type","AP","DV","ML","R","father"]` |
| Units / space | micrometres in Allen CCFv3 (2017), placed on the 10 µm grid (voxel = coordinate / 10) | `two_D_plots.py` divides by 10 against `average_template_10_8bit.nrrd`; `Fig 6a.ipynb` loads `annotation_10.nrrd` from `mouse_ccf/annotation/ccf_2017/`. The portal viewer also loads `average_template_10.nrrd` / `annotation_10.nrrd`. The probe files' ranges fit CCF (AP 4,968–10,556, DV 2,443–7,117, ML 2,875–5,877). |
| Midline / hemisphere | midline ML = 5,700 µm (570 voxels). In the portal index, Left means soma ML < 5,700 (195038_041: soma z 5,519, Left; 201232_056: soma z 5,877, Right; 192710_002: soma z 4,796, Left) | `two_D_plots.py` uses 570 as the voxel split; three-neuron check |
| Mirroring | **Portal and BSDC Allen-space SWCs keep the original hemisphere.** For the Zenodo SWCs and the paper tables, the authors mirrored every neuron whose soma has ML > 5,700: `if float(first_line.split(" ")[4])>5700: records[4] = 11400 - float(records[4])`, giving the output folder `swc_mirrow`. | Zenodo `Main Figs/fig_for_plot_neurons/fig_for_plot_neurons.ipynb` (SHA-256 f0c5b0ec3a8edc18607a79d8e65db9714e3447b1b21c57506b5538a6d16f83b8). This matches the probe: ML_portal + ML_zenodo = 11,400 at all 6,875 nodes of 201232_056. `All_soma_infor_HY_neurons.csv` (AP, DV, ML) is in the mirrored frame. |
| Soma | node 1, type 1, parent -1; one root per file | spec plus all probe files |
| Type codes | 0 undefined, 1 soma, 2 axon, 3 dendrite (documented). Also seen: 10 and 12, undocumented (192710_002 has one node of each; 201232_056 has 4 nodes of type 10) | portal spec; probe counts |
| Type use differs between files | 192710_002: 18,570 type 0 and 733 type 3. 195038_041: 9,575 type 0 and no 2 or 3 (Axon_only). 201232_056: 5,746 type 2 and 1,124 type 3. So axon nodes appear as 0 in some exports and as 2 in others. | probe |
| Radius | column 6 ranges from 0.0039 to 55.9 µm. These are not usable diameters. | probe |

Display rules that follow:
- Keep type 0 as "type 0 (unlabelled)" and types 10 and 12 as "type 10/12 (undocumented)". Do not relabel them.
- Draw with a constant line width, and say that the radii are not used.
- Label every file and figure "portal frame (true hemisphere)" or "Zenodo/paper frame (mirrored to left, ML < 5,700)".

UNRESOLVED:
- The paper's Methods are paywalled (no PMC copy; Europe PMC returns PMID 40140607 with no open full text). The reconstruction software, the registration method and an official statement of what types 10 and 12 mean are therefore not quoted from the paper. The conventions above come from the authors' code and the portal spec.
- Whether the BSDC zips carry the same type labels as the portal.

## 4. Recommended small public sample

47 portal SWCs, estimated at about 52 MB. All come from the portal route above, in the portal frame (true hemisphere).
Fetch them one at a time, a few seconds apart, with User-Agent `syn4rch-analysis/0.1`, and record URL, bytes and SHA-256 for each.

**A. The 18 Jiao neurons that also have a CC BY 4.0 Zenodo copy (25.7 MB, sizes known).**
Fetch the portal copy of each. This checks the frames on 18 neurons (12 Right, 6 Left): the Right ones should satisfy ML_p + ML_z = 11,400, and the Left ones should be identical.
195038_041 Oxt PVH L (Axon_only); 201232_040 Tac1 HY R; 201232_056 Tac1 MMm R; 201894_001 Adcyap1 PSTN R; 201922_073 Nts MPO R; 202559_043 Tac2 HY R; 202617_025 Sst ZI L; 202618_143 Sst LHA R; 202620_016 Sst PeF R; 210796_002 Pdyn VMH L; 211266_057 Pdyn HY R (Axon_only); 211353_070 Adcyap1 AHN L; 211461_026 Nts PH R; 211465_137 Tac2 ZI L; 211595_037 Sst HY L; 211595_063 Sst HY R; 211889_096 Pdyn HY R; 211890_005 Pdyn VMH R.

**B. Fill in the peptide lines and hemispheres (19 neurons, about 18 MB estimated).**
Each is Axon_and_dendrite, in its line's typical nucleus where possible, and the median-length candidate there.
192710_002 Pomc TU L (already probed); 202457_066 Pomc ARH R; 202813_035 Agrp ARH L; 202813_009 Agrp ARH R; 202825_072 Avp PVH L; 202825_065 Avp PVH R; 212248_071 Crh PVH R (Crh has no Left neurons); 201221_022 Nts LHA L; 210105_001 Orexin DMH L; 210104_007 Orexin LHA R; 195039_009 Oxt PVH R; 210073_044 Penk VMH L; 202170_016 Penk VMH R; 210007_142 Pmch LHA L; 210951_004 Pmch LHA R; 201230_079 Tac1 VMH L; 212202_026 Trh PVHd L; 212202_011 Trh PVHd R; 202824_062 Vip MPN R (Vip has no Left neurons).

**C. One dense nucleus for the intertwined-arbor view (10 ARH neurons, about 9 MB estimated).**
Pomc and Agrp, 5 Left and 5 Right, Axon_and_dendrite, chosen around the median length:
192710_050, 202458_008, 192710_035 (Pomc L); 202813_029, 202815_029 (Agrp L); 202460_045, 202460_053, 192711_018 (Pomc R); 202813_023, 202813_011 (Agrp R).

Coverage: all 16 peptide lines; both hemispheres for 14 of them (Crh and Vip have only Right neurons in the index); 16 nuclei; 2 Axon_only neurons.
Sizes for B and C are estimated at 40.5 bytes per µm of the projection length in `all_neurons_projection_length_fig1.csv`, calibrated on the 18 Zenodo files (range 30–140 bytes/µm). The true total is known only after fetching; stop if it passes 100 MB.
The selection was made with `scratchpad/scripts/select.py` from the portal index and that table.
