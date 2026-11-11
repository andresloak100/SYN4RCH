# Axon dataset sources for a larger real 3D benchmark

Dataset audit: metadata and source documentation only. Written 2026-11-10; documents read 2026-11-11, 01:12 to 01:40 UTC.
This file holds Part A only.

**Rules followed.** No SWC file, no archive of traces and no bulk table was requested, opened or downloaded; not even
a size check was made on a trace file or a trace archive. No account, no login, no picture code, no terms box. The
outcome tables and the outcome-derived fields of the index and the soma tables were not read. What was read: public
pages, site code, licence texts, public metadata listings and the authors' code. Copies, with URL, time, bytes and
SHA-256, are in the scratch directory named in A10.

Three things to know about how this was done. (1) The Zenodo archive mixes code, figures, tables and 36 SWCs; only
its file list and five code files were read, by range request, and notebook outputs were discarded. (2) The two portal
index files contain outcome-derived fields; they were parsed with a filter that drops those fields before anything is
kept. (3) The published Methods were read from a copy on a third-party web page (A10).

**Marks.** **V** = verified, source given. **I** = inferred. **U** = unverified.

**Quotations.** Only short labels, interface strings and openly licensed code are quoted. Longer passages are
paraphrased and located (URL, SHA-256, search key), so the exact wording can be re-read.

## Part A: additional Jiao traces

Jiao et al. 2025, Nat Neurosci 28:1073-1088, doi:10.1038/s41593-025-01919-0. This part extends
`reports/swc_access.md` and corrects some of its statements (A9). GB means 10^9 bytes.

### A0. Short answer

- **How much is left.** 7,143 of the 7,180 traces have never been inspected. 5,568 of them are in 134 animals of which
  no trace was ever inspected (genuinely new animals). The other 1,575 are new neurons in the 32 animals already touched.
- **Where they are.** Two official places: the portal (one file per neuron) and the BSDC record (one zip per animal;
  names and exact sizes are public). Zenodo and the authors' GitHub hold only the 18 neurons already on disk.
  NeuroMorpho has none. NeuroXiv lists a hypothalamus archive but does not describe it.
- **Access.** Both official places put their download buttons behind a login with a BSDC account. The portal's
  per-neuron address answered without login on 2026-11-10, but that is not the way the site offers the files.
- **Size of the new-animal set.** 134 zips, 5,222,228,699 bytes (exact, compressed). As single portal files: 5,568
  requests, about 6.4 to 7.1 GB (estimate; range 5.2 to 8.2).
- **Licence.** BSDC: CC BY-NC-ND 4.0. Portal: a copyright line and no licence. Private, non-commercial analysis and
  model fitting are covered by the BSDC licence. Sharing mirrored, resampled, simplified or converted traces is not.
- **Axon rule.** The authors' own code counts as axon every piece of cable that is not dendrite (node types 3 and 4).
  A file codes its axon as type 2 or as type 0; which one is known only after the file is fetched. The 16 development
  neurons are exactly the 16 type-2 files among the 37 inspected.
- **What the user would have to do.** For the official routes: register a BSDC account (name, organisation, phone and
  address fields), confirm the e-mail, log in with a picture code and tick the Terms of Service box.

### A1. Inventory recomputed

From `data/raw/jiao2025/swc/portal_index/mouse.neuron.info.json` (allowed fields only) and
`data/raw/jiao2025/swc/MANIFEST.tsv`. All **V** unless marked.

| quantity | value |
|---|---|
| neurons in the index | 7,180 (7,180 distinct file names, all `<6 digits>_<3 digits>.swc`) |
| animals (file name before the underscore) | 166 |
| source labels | 16; every animal has exactly one |
| hemisphere in the index | Left 3,071 / Right 4,109 |
| `reconstruction_type` | Axon_and_dendrite 6,871 / Axon_only 309 |
| previously inspected | 37 neurons from 32 animals (41 files) |
| development set | 16 neurons from 14 animals, all inside the 37 |
| genuinely new animals | 134 animals, 5,568 neurons (Left 2,382 / Right 3,186; Axon_and_dendrite 5,344 / Axon_only 224) |
| new neurons in animals already touched | 1,575 (871 in the 14 development animals, 704 in the 18 other inspected animals) |
| new animals per label | Adcyap1 20, Sst 18, Penk 16, Tac2 15, Pdyn 13, Nts 13, Oxt 9, Pomc 7, Tac1 7, Orexin 4, Agrp 4, Pmch 2, Avp 2, Crh 2, Vip 1, Trh 1 |
| new neurons per label | Nts 899, Adcyap1 883, Sst 872, Pdyn 636, Tac2 591, Penk 536, Tac1 404, Oxt 238, Pomc 213, Agrp 70, Orexin 62, Crh 60, Trh 57, Avp 31, Pmch 13, Vip 3 |
| size of a new animal | 1 to 140 neurons, median 33.5; 10 animals under 5; 66 with neurons on both sides, 68 on one side |
| Axon_only neurons in new animals | 224, in 45 of the 134 animals; 159 of them are in 7 animals where at least half the neurons are Axon_only (192713 and 201228 are entirely Axon_only) |
| animals shared with Li 2024 | 195037 (58 neurons, new) and 195038 (48 neurons, inspected) |

Notes.

- The development list is the one printed in `reports/d109_protocol.md`, section 16. The other agent's own manifest
  was not read, so that it is that agent's actual list is **U**.
- Animal = file prefix is also the portal's own rule: its filter takes the sample ID as the part of the file name
  before the underscore (**V**: portal `app.js`, search key `sampleIDs.includes`).
- Soma label. The index acronym and the `location` of `all_neuron_soma_infor.xlsx` agree for 7,138 neurons and differ
  for 42. Three acronyms occur only in the index (BST, V3, VM; one neuron each), which is why the index has 62 labels
  and the table 59 (**V**: computed from the allowed columns). The paper's label is the xlsx one.
- Label name. The index says `Orexin` where the xlsx says `Wt(orexin)` (103 neurons) (**V**).
- The portal project "Mouse LHA (Ji et al, 2025)" is not extra data. Its index lists 916 files; all 916 are among the
  7,180 (89 animals, all with soma label LHA) and they are served from the same folder (**V**).

### A2. Routes

| | route | what it serves | open without login? | new neurons it adds |
|---|---|---|---|---|
| R1 | portal, direct file address | one SWC per neuron, true hemisphere | answered without login on 2026-11-10; the site's own button needs login | all 7,143 |
| R2 | portal, export in the browser | the same files, zipped in the browser in lots of 100 | no: login | all 7,143 |
| R3 | portal, "Data Download" page | a link to the BSDC record | no: the link opens only after login | none by itself |
| R4 | BSDC record | 166 zips, one per animal | no: the buttons lead to the login page | all 7,143 |
| R5 | Zenodo 14588094, authors' GitHub | 18 Jiao SWCs, mirrored | yes | none (already on disk) |
| R6 | NeuroXiv download server | one undescribed archive `ION_Hypo_full_SWC_CCFv3.zip` | yes (public link) | unknown |
| R7 | NeuroMorpho.org | nothing from this paper | - | none |
| R8 | the authors | anything else, on request to the corresponding authors (paper's Code availability statement) | - | - |

**R1. Portal, direct file address.**

- Address: `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/hy/hy_neuron_download/<brain>_<nnn>.swc`.
  The site builds exactly this address for the projects `hy` and `lha` (**V**: `app.js`, function `downloadSWC`).
- One file per neuron. The prefrontal project offers a registered and an original file per neuron (`_reg`, `_orig`);
  the hypothalamus project offers one (**V**: same function).
- The file: plain text, LF line ends, 7 columns, no header lines, `text/plain`, `Accept-Ranges: bytes`. All 21
  recorded `Last-Modified` dates are 27 Aug 2024 (**V**: `reports/swc_access.md`, `MANIFEST.tsv`,
  `reports/swc_validation.md`). The index is dated 7 Dec 2023.
- Access as observed: HTTP 200 without login for 23 files on 2026-11-10 (3 by the user's probe, 20 by the prototype
  build) (**V**: `MANIFEST.tsv`). Not tested again today, because no SWC was requested.
- Access as offered by the site: the viewer's Download button runs `downloadFunc`, which goes on only
  `if (this.userInfo)` and otherwise opens the dialog "Please login first." (**V**: `app.js`, component
  `NeuronData.vue`). The paper describes the portal as a place where "registered users" can explore, visualise and
  download (**V**: published text, third-party copy, see A10). So a direct fetch without login goes around the
  interface. Nothing written on the portal forbids it, and nothing licenses it.
- No size, checksum or version is published per file. There is no `robots.txt` (**V**).

**R2. Portal, export in the browser.**

- After login: choose neurons (the filter takes sample IDs, labels, hemisphere and reconstruction type) and press
  Download. The browser fetches the R1 addresses 4 at a time and saves zips of 100 files, built in the browser and
  named `neurons_<b>_<B>.zip` (**V**: `app.js`, functions `onDownloadNeurons`, `downloadSWC`).
- Limit: 10,000 or more neurons are refused with "Please select 10000 or fewer neurons to download." (**V**).
  7,180 neurons give 72 zips; 5,568 give 56.
- Every downloaded address is reported to the site's Matomo statistics (**V**: `_paq`, `trackLink`).

**R3. Portal, "Data Download" page** (`https://mouse.digital-brain.cn/projectome/download`).

- Lists "Single-neuron datasets for mouse hypothalamus", 2024-06-12, pointing to the BSDC record. The link icon
  works only for a logged-in user; otherwise the login dialog opens (**V**: `1.js`, function `openTarget`).
- Holds the SWC column description (A4) and the copyright footer.

**R4. BSDC record.**

- `https://www.braindatacenter.cn/datacenter/web/#/dataSet/details?id=1800731789745954817`. The DOI
  10.12412/BSDC.1667278800.20003 resolves to it (**V**: DOI handle record). The paper's Data availability statement
  names this record and the portal, and nothing else (**V**: nature.com article page).
- Record (**V**: record metadata): title "Single-neuron datasets for mouse hypothalamus"; description "Over 7000
  single neurons."; released 2024-06-12; 166 files; licence code 1 = CC BY-NC-ND 4.0; copyright holder Xiaofei Wang;
  2,380 page views and 105,014 file downloads so far.
- Files (**V**): 166 zips named `<brain>.zip`. The 166 names are exactly the 166 animals of the index. Sizes are given
  in bytes: 6,659,859,965 in all (6.66 GB; the page's "6.20 GB" is the binary unit). Smallest 380,730 (202662, 2
  neurons), median 33,353,978, largest 150,210,256 (202571, 140 neurons). No checksums. The whole list is in the appendix.
- What a zip contains is not documented anywhere read (**U**). Inferred (**I**):
  - a folder per brain with `swc_allen_space/<nnn>.swc`: the authors' notebook reads that path from what it calls the
    BSDC download;
  - more than one copy of each neuron: the zips hold 927,557 compressed bytes per neuron, about 2.5 to 2.8 times what
    one portal-style SWC per neuron would need at the compression measured on the 18 Zenodo SWCs (3.26 to 3.51, median
    3.42); the authors' client code names a `swc_raw` folder beside the registered one.
- Access (**V**: site code, module `cf45` and functions `downloadfile`, `copyUrl`, `streamDownload`): the page shows
  the file list and the buttons Download, "Copy Download URL" and "Batch Download" to every visitor (access category
  3, `showFile` true). Using any of them without a stored login token sends the browser to
  `https://www.braindatacenter.cn/iam/admin/#/login`. The page carries no "login required" text. No
  request-to-the-author step applies to this record (**I**: that step is coded for another access category).
- Whether the file addresses themselves answer without a token was not tested (**U**), and should not be tested
  without the user: the site's terms forbid going around its access measures (A5).
- The site records each download against the account, and gives data owners a view of downloads by account
  (**V**: names in the site's API map).
- The page shows a general tip for folders with several archives: download all of them before unpacking. The zips
  here are named per animal, so a subset is presumably complete in itself (**I**).

**R5. Zenodo 14588094 and the authors' GitHub.**

- Zenodo: CC BY 4.0; one archive of 918,302,741 bytes with 1,153 members. 36 are SWCs, all in
  `Supplementary Fig/Supplementary Fig4/swc/`: 18 Jiao neurons (23,725,631 bytes) and 18 MouseLight files (**V**:
  central directory only).
- GitHub `xulab2022/Single-Neuron-Projectome-of-Hypothalamus` (named in the paper's Reporting Summary), branch
  `master`: the same 36 SWC names; no licence file (**V**: file tree).
- Nothing new in either.

**R6. NeuroXiv.**

- The download server lists `ION_Hypo_full_SWC_CCFv3.zip`, dated 2025-12-31 04:17:41, with no description (**V**:
  file listing of `https://download.neuroxiv.org`).
- Its README (last changed 2026-02-12) describes two ION sets only, prefrontal cortex and hippocampus. The portal home
  page still shows totals without a hypothalamus set (19,406 axons = 1,876 + 16,380 + 1,150) (**V**). So the content
  of the archive, its neuron count, its source version and its frame are undocumented (**U**). From the name alone: a
  NeuroXiv-processed copy of ION hypothalamus traces (**I**).
- NeuroXiv's documented processing of ION data (**V**: preprint Methods; **U** that it was applied to this archive):
  node types re-coded to 1 soma, 2 axon, 3 basal dendrite, 4 apical dendrite; files kept only if they are a single
  tree with one type-1 root and consistent types; mapped with mBrainAligner; **resampled so that parent and child
  nodes are 1 um apart**; renamed `<resource>_<full/local>_<brainid>_<neuronid>_..._<atlas>`. Such files would not
  hold the authors' own nodes.
- Size: the listing gives names and dates only. Not measured, because that needs a request on the archive (**U**).
- Terms and licence chain: see A5.

**R7. NeuroMorpho.org.** Not there (**V**). No record cites PMID 40140607 or the DOI; the same query form returns 23
records for a control paper. The whole site has 1,037 mouse hypothalamus records, far fewer than 7,180. No archive is
named Jiao or Xu (archive names alone could not exclude a deposit under another name; the first two facts do).

**R8. Other things checked.**

- The authors' Python client `neuron-vis` (`IONData.py`) fetches SWCs from an institute server at an internal address
  (10.10.48.110). The public host name left in a comment (bap.cebsit.ac.cn) refused the connection. Not a route (**V**).
- Older host `mouse.braindatacenter.cn`: expired certificate. Not a route (**V**).
- Europe PMC has no open full text. The Author Correction of 2 Dec 2025 only fixes the spelling of one author's name (**V**).

### A3. Frame, units, registration and node types, by route

| | R1/R2 portal file | R4 BSDC zip | R5 Zenodo / GitHub | R6 NeuroXiv archive |
|---|---|---|---|---|
| hemisphere | true side (**V** on 23 files) | true side (**I**: the authors mirror these files themselves) | mirrored to ML < 5,700 (**V** on 18) | **U** |
| columns, units | id, type, AP, DV, ML, radius, parent; micrometres (**V**) | same (**I**) | same (**V**) | micrometres (**I**) |
| atlas | Allen CCFv3, 2017 annotation, 10 um grid, 13,200 x 8,000 x 11,400 um (**V**: authors' code) | same (**I**) | same (**V**) | CCFv3 (from the name) |
| nodes | as traced | as traced (**I**) | as traced | resampled to 1 um, if the standard processing was used (**U**) |
| type codes | 0, 1, 2, 3 documented; 10, 11, 12 present and undocumented (**V**) | **U** | as the portal file (**V** on 4 pairs) | 1 to 4, if the standard processing was used (**U**) |
| dendrites | present in Axon_and_dendrite files, absent in Axon_only (**V** on 37, **I** for the rest) | **U** | as labelled (**V**) | **U** |

Facts common to all routes.

- Midline ML = 5,700. Mirroring rule: when the soma has ML > 5,700 the whole neuron gets ML' = 11,400 - ML (**V**:
  Zenodo notebook `fig_for_plot_neurons.ipynb`; pyswcloader `swc_preprocess`).
- Side names. The index calls ML < 5,700 "Left". pyswcloader's comment says somas are mirrored "to the left brain".
  neuron-vis adds cable with z > 5,700 to `projectionRight` (**V**). That this is the animal's anatomical left is not
  stated in anything read (**U**).
- Registration (**V**: published Methods and preprint Methods agree): the whole brain, with its traced neurons, was
  registered to CCFv3 with the group's earlier method (Gao et al. 2022). Landmark regions were segmented from the
  cytoarchitecture channel; a diffeomorphic (symmetric normalisation) registration was computed in ANTs; the transform
  was applied to the traced points. Gao 2022's figure legend describes an affine step followed by a non-rigid step.
- Accuracy statement (**V** for the method, **U** for the number): the authors checked 4 brains and 3 structures
  (fimbria, median eminence, interpeduncular nucleus) and measured the distance between each structure's centre of mass
  after registration and in the atlas. The result is shown only in a figure (preprint Extended Data Fig. 2a; in the
  final paper Supplementary Fig. 1b, by its main-text reference). No number is in any text read. The Supplementary
  Information PDF (20.7 MB) was not opened.
- Evidence already on file (`reports/swc_validation.md`): somas are 0 to 1.92 um from the paper's table; 3 PVH somas
  lie 9 to 33 um outside the PVH mesh; 5 of 37 neurons run past the caudal end of the atlas (AP > 13,192).
- Truncation. The paper itself notes that some descending axons may continue beyond the imaged volume into the spinal
  cord (**V**: published text, third-party copy). The authors' code moves coordinates beyond the atlas box to its edge
  before measuring (**V**: pyswcloader `swc.py`, lines 39-47), so their tables do not treat that cable as the raw
  files do.
- Imaging and tracing (**V**): fMOST at 0.32 x 0.32 x 1 um; tracing with the in-house Fast Neurite Tracer, each step
  confirmed by a person; every neuron checked by a second person; random quality checks per brain; tracers blind to
  soma location and label; neurons with tracing errors, or with a soma outside the hypothalamus, were excluded.
- The radius column is not a diameter (`reports/swc_access.md`).

### A4. Axon-type convention

What the sources say.

1. Portal Download page (**V**: `1.js`): the type column is introduced with "eg: 0,1, 2, 3..." and four codes are
   named: "0: undefined", "1: soma", "2: axon", "3: dendrite". Nothing on units, axes or frame.
2. pyswcloader README (**V**): the SWC format follows the standard specification (neuronland.org link).
3. Which tool made the paper's numbers (**V**): the published Methods say soma location, axon length and terminal
   number per region were computed with pyswcloader (third-party copy); the public Code availability statement says
   this library was used for the feature calculations. The 2023 preprint named another tool (NeuronView); the
   published text supersedes it.
4. pyswcloader's rule (**V**, commit bc64d46):
   - `total_length` and `projection_length` add a segment only when
     `_is_axon not in [3, 4] and _parent_is_axon not in [3, 4]`, the two variables being the node's type and its
     parent's type (`swc.py` lines 70-72; `projection.py` lines 20-22).
   - So **axon = all cable that is not dendrite**. Types 0 and 2, the first segment out of the soma, and the marker
     types 10, 11 and 12 all count. A segment that touches a type-3 or type-4 node does not.
   - Terminals are the leaf nodes whose type is not 3 or 4 (`projection.py` line 116).
   - Before measuring, the neuron is mirrored to z < 5,700 and clamped to the atlas box; extra roots are dropped and
     extra type-1 nodes are treated as type 0 (`swc_preprocess`).
5. neuron-vis, the authors' viewer (**V**: `SwcLoader.py` lines 116-138): dendrites are the edges whose second point
   has type 3; the axon starts at the first edge that leaves the soma, or a type-3 point, towards a point that is not
   type 3. The same rule: what is not dendrite is axon; an axon may leave from a dendrite.
6. The paper's Methods do not say how dendrites were traced or labelled. One figure legend colours soma, axon and
   dendrites separately (**V**: published text, third-party copy).
7. `reconstruction_type` has no written definition in anything read. The portal's filter knows four values
   (`axon_only`, `axon_and_dendrite`, `dendrite_only`, `undefined`); the hypothalamus index uses two (**V**). Its
   meaning is taken from the names and from the 37 inspected neurons (**I**).

What the 37 inspected files show (**V**: `reports/swc_validation.md`, the index).

- All 35 Axon_and_dendrite files contain type-3 nodes; both Axon_only files contain none.
- Each file codes its non-dendritic cable as type 2 or as type 0, never both: 16 neurons use type 2, 21 use type 0.
- **The 16 type-2 neurons are exactly the 16 development neurons.**
- No animal has both codings: 14 animals are type 2, 18 are type 0. The 5 animals with two inspected neurons agree
  within the animal. Coding by animal or by export batch is therefore likely (**I**).

Not knowable for a file until it is fetched.

- Whether its axon is coded 2 or 0 (probably the same for its whole animal, **I**).
- How many 10 / 11 / 12 marker nodes it has. Their meaning is not documented in anything read.
- Whether it really has, or lacks, dendrite nodes as the index says.
- Whether cable leaves the atlas volume caudally.
- Its size, node count and node spacing.
- Whether the portal copy, the BSDC copy and the version used in the paper are the same (one portal soma differs
  from the paper's table by 1.92 um).
- Whether it has a single root (the authors' code guards against extra roots).

Rule proposed for declaration before any fetch: use the authors' rule (axon = segments whose two end nodes are not
type 3 or 4), report type-0 and type-2 animals separately, and never relabel raw files. A stricter rule (explicit
type 2 only) would remove an unknown share of animals: 18 of the 32 inspected animals are type 0, but they were not
chosen at random (**U** for the 134).

### A5. Licence and terms

The texts (**V**).

- **Portal.** Footer: "Copyright © 2023 Digital Brain. All Rights Reserved." No licence, no terms of use, no dialog
  (`1.js`; `reports/swc_access.md`). Its login and registration use the BSDC account system (`app.js`).
- **BSDC record.** "License Agreement": CC BY-NC-ND 4.0 (record value `license: 1`; dictionary entry 1 is
  "CC BY-NC-ND 4.0"). The attribution and acknowledgement text is quoted in `reports/swc_access.md`. The record also
  carries a restricted-sharing note that reads "unpublished" in Chinese; the page shows that note only for records
  that need an access request (**I**: not shown for this record).
- **BSDC site terms** (Chinese; paraphrased from `bsdc_terms_of_use.json`): using the site means accepting the terms.
  Data may be browsed, downloaded and used, subject to the data owner's licence. Services other than data sharing are
  for non-commercial use. Forbidden: putting an excessive load on the site, interfering with it, and going around
  measures that prevent or limit access. The site logs browsing, and the IP address of each download, for statistics.
  Registration details must be true.
- **CC BY-NC-ND 4.0 legal code**, section 2(a)(1): the material may be reproduced and shared for non-commercial
  purposes, and the licensee may "produce and reproduce, but not Share, Adapted Material for NonCommercial purposes
  only". Section 3 repeats that there is no permission to share adapted material.
- **Zenodo.** `cc-by-4.0` (record metadata).
- **NeuroXiv.** CC BY-NC (paper's Data availability statement and README). The README allows research, educational
  and personal use, and forbids commercial use and distributing the dataset or derivative works; continued use counts
  as acceptance of changed terms. The source record is CC BY-NC-ND, which gives no permission to share a re-coded,
  resampled copy; whether NeuroXiv has the authors' permission is not stated (**U**).
- **Code.** pyswcloader: MIT in the PyPI metadata, no licence file on GitHub. neuron-vis: its README ends with
  "GFDL v1.2+", the repository field says not specified. Zenodo code: CC BY 4.0.

What each allows. This is a reading of the texts above, not legal advice.

| use | portal copy | BSDC copy | Zenodo 18 | NeuroXiv copy |
|---|---|---|---|---|
| private analysis | no licence; treat as private research use only (**U**) | yes, non-commercial (**V**) | yes | yes, academic (**V**) |
| fitting models, kept private | as above (**U**) | yes: adapted material may be produced, not shared (**V**) | yes | yes |
| sharing model weights | no basis (**U**) | unsettled whether weights are adapted material; ask the authors (**U**) | yes, with attribution | not if they count as derivative work (**U**) |
| sharing derived geometry (mirrored, resampled, simplified, completed) | no | no (ND) (**V**) | yes, with attribution | no (README) |
| sharing unmodified files | no | yes, non-commercial, with attribution (**V**) | yes | README says no |
| figures and summary numbers | numbers yes; trace images privately or with permission (**I**) | numbers yes; a drawn trace is arguably adapted material: private, or ask (**I**) | yes | as BSDC (**I**) |

### A6. Download needs

| set | animals | neurons | BSDC zips: requests, bytes (exact) | portal files: requests, bytes (estimate) |
|---|---:|---:|---|---|
| everything | 166 | 7,180 | 166 zips, 6,659,859,965 | 7,180 files; 8.2 GB by E1 (6.7 to 9.9) or 9.0 GB by E2 (7.6 to 10.4) |
| genuinely new animals | 134 | 5,568 | 134 zips, 5,222,228,699 | 5,568 files; 6.4 GB by E1 (5.2 to 7.7) or 7.1 GB by E2 (6.0 to 8.2) |
| inspected animals, not development | 18 | 725 (704 not inspected) | 18 zips, 605,748,777 | 704 files; about 0.8 GB |
| development animals | 14 | 887 (871 not inspected) | 14 zips, 831,882,489 | 871 files; about 1.0 to 1.1 GB |

A zip always holds the whole animal, including neurons already inspected. Of the 7,180 portal files, 23 are already on
disk in the portal's frame; 14 more neurons are on disk only as mirrored Zenodo copies.

How the portal estimate was made. No file was opened.

- The `bytes` column of `MANIFEST.tsv`: 41 files, 46,551,940 bytes; mean 1,135,413; median 838,646; SD 699,917;
  range 142,096 to 3,055,432.
- These are 37 neurons. Put in the portal's form (23 portal files as they are; 14 Zenodo-only files less one byte per
  line for the carriage return, node counts from `reports/swc_validation.md`): mean 1,147,248 bytes; median 929,633;
  SD 688,812; quartiles 629,975 and 1,575,691; range 142,096 to 3,002,962. Bootstrap 95% interval of the mean: 0.94 to
  1.37 MB.
- E1 = that mean times the number of neurons. E2 is a cross-check that uses each animal's exact zip size: the 37
  files hold 1.35 times the zip bytes per neuron of their animals (interval 1.14 to 1.57; correlation 0.58), so an
  animal's portal bytes are put at 1.35 times its zip bytes.
- The 37 were not drawn at random (18 are the authors' showcase neurons, 19 were picked near a median length), so the
  bias is unknown. Planning figure: about 6 to 8 GB for the new animals and 8 to 10 GB for everything.
- On the wire the files may be smaller if the server compresses text. Not tested (**U**).
- Time at one request every 3 seconds: 4.6 hours for 5,568 files, 6.0 hours for 7,180, plus transfer.
- Portal export in the browser: 56 zips of 100 for the new animals, 72 for everything.
- BSDC unpacked size is not documented. At the compression measured on the Zenodo SWCs (3.4 times), the 134 zips
  would unpack to about 18 GB and all 166 to about 23 GB (**I**).
- NeuroXiv: 1 request; size not listed (**U**).

Per animal: the appendix gives, for each of the 166 animals, its label, neuron count, hemisphere counts, Axon_only
count, status, exact zip bytes and the E2 estimate of its portal bytes.

Integrity: neither source publishes checksums. Record size and SHA-256 at download. BSDC sizes can be checked against
the appendix.

### A7. What the user would have to do personally

| route | registration | login | terms | notes |
|---|---|---|---|---|
| R1 portal, direct address | none | none | none shown | the user decides whether to use an address the site offers only behind login; keep it slow and private |
| R2 portal, browser export | BSDC account | yes | Terms of Service box at login | runs in the user's own browser |
| R4 BSDC zips | BSDC account | yes | Terms of Service box at login | 134 or 166 downloads, or the page's Batch Download |
| R5 Zenodo, GitHub | none | none | none | nothing new there |
| R6 NeuroXiv archive | none | none | use counts as acceptance (README) | one public link; content undescribed; licence chain unclear |
| R8 the authors | - | - | - | an e-mail to the corresponding authors is the only way to get permission to share derived traces or models |

The BSDC account, as coded in the forms (**V**: registration and login form code; the live forms were not opened).

- Two steps: basic information, then a detailed profile.
- Fields marked required in the form's rules: e-mail, password (8 to 24 characters, letters and digits), password
  again, sex, real name, nickname, telephone, organisation, title, academic title, department, profession, city,
  address, postcode.
- The account must be confirmed from a link sent by e-mail before login works.
- Login: account, password, a 4-character picture code, and a box saying the Terms of Service were read and accepted.
- Alternative: login with a CSTCloud passport (escience.cn).
- Not known: whether accounts are approved by hand, and whether foreign phone numbers and addresses are accepted (**U**).

Creating the account, solving the picture code and accepting the terms are the user's own steps; an agent must not do
them.

### A8. Small probes that would settle the open points

None was done. Each needs explicit approval.

1. Zip contents and unpacked size: after the user logs in, read one zip's file list, or unpack one small animal
   (202662: 380,730 bytes, 2 neurons).
2. Axon coding of each new animal (type 2 or type 0): the first and last kilobyte of one portal file per animal (134
   range requests of about 2 KB) would probably show it. This is geometry access.
3. NeuroXiv archive: one size request and one read of its file list would show how many files it holds and whether
   they are resampled.
4. The registration accuracy number: it is in a figure of the free Supplementary Information PDF (20.7 MB).
5. Portal copy against BSDC copy: compare one animal once both are in hand.

### A9. Corrections to earlier notes

These files were not edited.

1. `reports/swc_access.md`, portal table, "Login gate: None on the download path". The server does not check, but the
   interface does: the Download button works only for a logged-in user (R1).
2. `reports/swc_access.md`, BSDC section. The zips are 0.38 to 150.2 MB, not 9 to 122 MB. "6.20 GB" is the binary
   unit: 6,659,859,965 bytes.
3. `reports/swc_access.md`, "UNRESOLVED: whether those downloads need a login". Resolved for the site's buttons:
   they do (R4).
4. `reports/swc_access.md`: the project list is at `/projectome/2/srv/api/v1/public/projects`;
   `/projectome/2/srv/info/projects.json` returns 404.
5. `reports/swc_prototype.md`, "about 6.2 GB" for the full set. That is the compressed BSDC total. As single files
   the set is estimated at 8 to 10 GB.
6. `reports/conversion_map.md`, "the NeuroXiv mirror is CC BY-NC 4.0". NeuroXiv documents no hypothalamus set; only
   an undescribed archive name exists (R6).
7. `reports/d109_jiao_conventions.md`, open point b (which tool made the length table): the published Methods name
   pyswcloader. Open point h (name of the ML < 5,700 side): the authors' code calls it left.

### A10. Sources

Copies are in
`/private/tmp/claude-501/-Users-andresarriaga-Documents-SYN4RCH/6fc8191b-130f-4b76-b6a8-69340161f99e/scratchpad/bench/axon_audit_jiao_access/`
(`web/`, `zenodo/`, with `LOG.tsv` in each; scripts in `scripts/`). Nothing was written to the repository except this
file. Requests used a generic client name and no credentials. The BSDC open endpoints were sent the two public
constants that the site's own page sends for every visitor (`accessId`, `appId`); one of them is a read-only search
call sent as POST. nature.com, biorxiv.org and gitee.com needed cookies held in memory to pass their redirects;
nothing was stored.

| what | URL | read (UTC) | bytes | SHA-256 |
|---|---|---|---:|---|
| Portal site code `app.js`: per-neuron download, batch export, sample-ID filter, login calls | `https://mouse.digital-brain.cn/projectome/js/app.js` | 2026-11-11 01:12:56 | 6,482,711 | `58502a6d7075576a50b6e99a501e8e7daed63d696560e52495aee64dc4df8099` |
| Portal `1.js`: Data Download page (SWC format text, dataset links, login gate on the links, copyright footer) | `https://mouse.digital-brain.cn/projectome/js/1.js` | 2026-11-11 01:13:55 | 148,908 | `dd54ab673039290231a4c2968bf49b907275ff052fc9605226eb3844300818d3` |
| Portal public project list | `https://mouse.digital-brain.cn/projectome/2/srv/api/v1/public/projects` | 2026-11-11 01:12:59 | 7,821 | `7154a66174ea955c05766310ea07f599260d48a87cc68eac1ca4177788680b82` |
| Portal index, project `hy` (same bytes as the repository copy; only allowed fields parsed) | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/hy/mouse.neuron.info.json` | 2026-11-11 01:15:02 | 1,980,365 | `fdf1fa2eff099bc755e5218dc6437e772f12fe0a934c01a4072d846d9b3f4da3` |
| Portal index, project `lha` (only allowed fields parsed) | `https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/lha/mouse.neuron.info.json` | 2026-11-11 01:14:55 | 253,755 | `186d66593a708e3a11828629e22b9c0b8d9c639ddb6f5c5b83abc1b55b668e27` |
| Portal `robots.txt` request (returns the app page: no robots file) | `https://mouse.digital-brain.cn/robots.txt` | 2026-11-11 01:33:53 | 1,296 | `87d1fb8e4676e922c8e1be0856f088504096eeae425af42ae8d3aeba51f3fab0` |
| BSDC site code `app.833b21ae.js`: API map, login guard (module `cf45`), download helper, interface strings | `https://www.braindatacenter.cn/datacenter/web/js/app.833b21ae.js` | 2026-11-11 01:15:30 | 64,914 | `1c6a5748223d6ab26d8a736834bdd7d8dfaa841b4c5ccdddc034bcbae7febfad` |
| BSDC dataset page code `chunk-ea9c37ce.b0c5bbc6.js`: access states, Download / Copy URL / Batch Download | `https://www.braindatacenter.cn/datacenter/web/js/chunk-ea9c37ce.b0c5bbc6.js` | 2026-11-11 01:15:49 | 142,649 | `d8d7f3b41a9d96e4b7f8fb7dbe3d9d4e1a28ab4ad91931e8346dfe5b22bc7fc3` |
| BSDC record metadata, with the 166 zip names and byte sizes (open endpoint the page calls; sent the web app's public `accessId` / `appId` constants, no cookie, no token) | `https://www.braindatacenter.cn/api/collect/open/web/datasetInfo/getDatasetInfoById?datasetInfoId=1800731789745954817` | 2026-11-11 01:17:35 | 44,574 | `5091c29deff9cfc9b89d5c7e9dcd40b72802256657c494eef33276e4fcf45fe1` |
| BSDC licence dictionary (code 1 = CC BY-NC-ND 4.0) | `https://www.braindatacenter.cn/api/cms/dict/open/get/childInfos/datacenter_dataset_license` | 2026-11-11 01:18:26 | 4,628 | `26b0cc48f56402d999df0953306dfcf7ab6f74428142428ddaa96ca85f7a39c0` |
| BSDC site Terms of use, Chinese text (read-only search call made by the `#/termsForUsage` page) | `https://www.braindatacenter.cn/api/cms/open/searchArticles` | 2026-11-11 01:18:50 | 11,188 | `41560d9959314a1303247a115b8f3507bc3889c879546c1942476cebd27b7912` |
| BSDC project description | `https://www.braindatacenter.cn/api/collect/open/web/projectInfo/getProjectWebDetail?projectId=1585507260749971457` | 2026-11-11 01:35:31 | 1,633 | `3ade72e00de12056cceeeee49799b882a7e4c7400f68abd51d1c449f2ad7c9f4` |
| BSDC registration form code (read only; nothing submitted) | `https://www.braindatacenter.cn/iam/admin/js/119.beb04638.js` | 2026-11-11 01:19:31 | 21,969 | `06be998a7756bc467888cd8d0df08bd7cf9800c7e7d4d9df498f451d8af4f2bb` |
| BSDC login form code (read only; nothing submitted) | `https://www.braindatacenter.cn/iam/admin/js/504.59103c90.js` | 2026-11-11 01:19:32 | 20,178 | `5e964a871361b057256c5259c079899f095c86f550e3a498e0a920dc06c1cf8f` |
| DOI handle record of 10.12412/BSDC.1667278800.20003 | `https://doi.org/api/handles/10.12412/BSDC.1667278800.20003` | 2026-11-11 01:32:22 | 475 | `5e5f8d2f002ca084ccc6bb386bc7c9659a53ed05d3d95846b32624fbe7f7d452` |
| Nature Neuroscience article page (public sections: Data availability, Code availability, Rights, Change history) | `https://www.nature.com/articles/s41593-025-01919-0` | 2026-11-11 01:20:59 | 595,947 | `78bf73b5b34289cfb5f890e4f3fe9c80e74d3bb9a32c0996d0d52ccf3a27b0b1` |
| Author Correction page | `https://www.nature.com/articles/s41593-025-02184-x` | 2026-11-11 01:21:03 | 369,175 | `e006a3b299199fda879756d231889a2329e3c305723ff685d27a15145b88cd1c` |
| Europe PMC record (no open full text) | `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:10.1038/s41593-025-01919-0&resultType=core&format=json` | 2026-11-11 01:23:02 | 27,024 | `8dddf3626f3d485a69849725c3b3729a8bf6caa427a9a5c9108463d414d7ac21` |
| bioRxiv metadata of the preprint (one version, CC BY-NC-ND) | `https://api.biorxiv.org/details/biorxiv/10.1101/2023.05.25.542241` | 2026-11-11 01:22:01 | 2,651 | `b1c71d5bffe3b29d50e7f6954e713be8717ed1ddac88f00bee66bf33700b53dc` |
| bioRxiv preprint full text (Methods read) | `https://www.biorxiv.org/content/early/2023/05/25/2023.05.25.542241.source.xml` | 2026-11-11 01:22:10 | 174,590 | `a5833a7cfdb5a9a78754eedd8ad903980001d5c04908dd4ee4d70984d50903c5` |
| Gao et al. 2022 article page (public sections; the method the paper refers to) | `https://www.nature.com/articles/s41593-022-01041-5` | 2026-11-11 01:25:02 | 484,360 | `1ee3d94866ce6408713497dd20cabeed7ccb48b0d502bb45ee8272e2fae970bd` |
| PyPI metadata of pyswcloader | `https://pypi.org/pypi/pyswcloader/json` | 2026-11-11 01:25:23 | 5,481 | `699a3fa6529369d21fab579b19c954cddc308526c3e50fc3954f8efe4f3dd845` |
| GitHub metadata of txgxxx/pyswcloader | `https://api.github.com/repos/txgxxx/pyswcloader` | 2026-11-11 01:25:24 | 5,065 | `50416fe4ae8a48fa4c986851c0173467e113fa5a2c7dfd503fd5df76d46d4471` |
| pyswcloader README at commit bc64d46 | `https://raw.githubusercontent.com/txgxxx/pyswcloader/bc64d464976dd4fad88619b64d91acd5999526fd/README.md` | 2026-11-11 01:25:33 | 1,048 | `69a6721089d03adb9d3d62867f01b5fd53220eb6c2a7381c2558912378d88f5a` |
| pyswcloader `swc.py` at commit bc64d46 | `https://raw.githubusercontent.com/txgxxx/pyswcloader/bc64d464976dd4fad88619b64d91acd5999526fd/pyswcloader/swc.py` | 2026-11-11 01:25:33 | 7,419 | `5bbd41f890979edfd809dbd91224448fea00df5be5144d0350d10df17ee41c22` |
| pyswcloader `projection.py` at commit bc64d46 | `https://raw.githubusercontent.com/txgxxx/pyswcloader/bc64d464976dd4fad88619b64d91acd5999526fd/pyswcloader/projection.py` | 2026-11-11 01:25:34 | 8,946 | `814107e57600324f4c11a8f3d0f9de60054dc158c4c4b8bc23951276fb740a80` |
| pyswcloader `brain.py` at commit bc64d46 | `https://raw.githubusercontent.com/txgxxx/pyswcloader/bc64d464976dd4fad88619b64d91acd5999526fd/pyswcloader/brain.py` | 2026-11-11 01:25:35 | 2,385 | `f1d324a41549e6789845b3cc2908b7d69a6b5fd2b68a02cc63d742699323ce11` |
| Zenodo record 14588094 metadata | `https://zenodo.org/api/records/14588094` | 2026-11-11 01:26:26 | 5,498 | `9771c8d05c3b1f0103252157432d5efe99caa4bafb10443f52ed30132731cdba` |
| Zenodo archive: central directory only (member names and sizes) | `https://zenodo.org/api/records/14588094/files/Code%20for%20Hypothalamic%20projectome.zip/content` | 2026-11-11 01:27:26 | 143,430 | `3268fda3a5f4caf7f6114c46c299235f4fc219939985851367f263c590aee59b` |
| Zenodo archive member `Main Figs/fig_for_plot_neurons/fig_for_plot_neurons.ipynb` (hash of the member; only cell sources kept) | `https://zenodo.org/api/records/14588094/files/Code%20for%20Hypothalamic%20projectome.zip/content` | 2026-11-11 01:27:27 | 35,262 | `f0c5b0ec3a8edc18607a79d8e65db9714e3447b1b21c57506b5538a6d16f83b8` |
| Zenodo archive member `Extended Data Fig/Extended Data Fig4&5_axon_showcase/two_D_plots.py` | `https://zenodo.org/api/records/14588094/files/Code%20for%20Hypothalamic%20projectome.zip/content` | 2026-11-11 01:27:28 | 30,513 | `047568a3cd47ef9d813f7afab0a0cd21280f5eff3b9673fcb7f24fe84a4a62b1` |
| gitee page of bigduduwxf/neuron-vis (text summary served to non-browser clients) | `https://gitee.com/bigduduwxf/neuron-vis` | 2026-11-11 01:27:53 | 4,977 | `e4b53b53399ddf1809e229794fa474e446d795ee50f833c4cdf3e00b612bb096` |
| neuron-vis `neuronVis/IONData.py` (branch master) | `https://gitee.com/bigduduwxf/neuron-vis/raw/master/neuronVis/IONData.py` | 2026-11-11 01:28:14 | 13,385 | `8b8e200a3366bf03837a861a83a2bfbd6186efcd782a30c569b6ab91e2130988` |
| neuron-vis `neuronVis/SwcLoader.py` (branch master) | `https://gitee.com/bigduduwxf/neuron-vis/raw/master/neuronVis/SwcLoader.py` | 2026-11-11 01:28:17 | 8,575 | `8d46c8b38d07371a05d5236556a10e38fa6d379643580da4fe20f8a014e1f927` |
| neuron-vis `neuronVis/Project.py` (branch master) | `https://gitee.com/bigduduwxf/neuron-vis/raw/master/neuronVis/Project.py` | 2026-11-11 01:28:23 | 2,052 | `d330f27765b8296c4039f4ecf161a5d4b561411f671c3f8dcca450aaed766e75` |
| GitHub metadata of xulab2022/Single-Neuron-Projectome-of-Hypothalamus | `https://api.github.com/repos/xulab2022/Single-Neuron-Projectome-of-Hypothalamus` | 2026-11-11 01:39:31 | 6,474 | `d5a34c50d2abf7e2104d79539555aee682032e08eb761a3d517a5ef078e27e97` |
| File tree of that repository, branch master (names and sizes only) | `https://api.github.com/repos/xulab2022/Single-Neuron-Projectome-of-Hypothalamus/git/trees/master?recursive=1` | 2026-11-11 01:39:59 | 391,377 | `b88b506fdc9e779ec32ad5fe56b66545b81f2ae81f04b654389eb091dcd24de2` |
| NeuroXiv paper page (public sections: Data and Code availability) | `https://www.nature.com/articles/s41592-025-02687-2` | 2026-11-11 01:29:21 | 331,768 | `7284e58cc5d9c55ec1b6b385953fff3c0f48e4b42c2beded6da1fb81a9f6dbc3` |
| NeuroXiv download page | `https://download.neuroxiv.org/` | 2026-11-11 01:29:29 | 6,859 | `424897eed6ccc95ab1dc252bb1a8e318e05fdaa5ab8606613c44f9a092bd32aa` |
| NeuroXiv download server file listing (names and dates; no sizes) | `https://d36ajqhpoeuszk.cloudfront.net/list-files` | 2026-11-11 01:29:39 | 1,594 | `05391fdbe16f3b712fd92464d40225a1c6b165007ca736d6f6d11be1699bab39` |
| NeuroXiv download server README (datasets, terms of use) | `https://d36ajqhpoeuszk.cloudfront.net/README.md` | 2026-11-11 01:29:52 | 9,393 | `a0b8c830ecb43c6c467eb2fd4549f06647dbdfa58538f44e402cf113434dcc43` |
| NeuroXiv portal home page (counts shown) | `https://neuroxiv.org/` | 2026-11-11 01:29:30 | 18,223 | `f93331a949964ed83801309e8083df90ec465e8232b55a35bc9d185f10a8f0a9` |
| bioRxiv metadata of the NeuroXiv preprint (CC BY-NC-ND) | `https://api.biorxiv.org/details/biorxiv/10.1101/2024.09.22.614319` | 2026-11-11 01:33:03 | 1,929 | `56b09fb687ef5edcfc520278be4fa2c0aa01fbf37ce7e55a88e2fab61f7c94f3` |
| NeuroXiv preprint full text (Methods on standardisation read) | `https://www.biorxiv.org/content/early/2024/09/24/2024.09.22.614319.source.xml` | 2026-11-11 01:33:11 | 175,078 | `c114665d92a4c2845da54a68eeff8234f5fab41029e34c87e49164298e4ac977` |
| NeuroMorpho archive names (1,081) | `https://neuromorpho.org/api/neuron/fields/archive?size=2000` | 2026-11-11 01:31:29 | 12,422 | `f31069431465c4dbad394a887a0cee11c4f5311f1a82fb0bd23a307fd2c4e9fa` |
| NeuroMorpho count of mouse hypothalamus records (size=1 query) | `https://neuromorpho.org/api/neuron/select?q=brain_region:hypothalamus&fq=species:mouse&size=1` | 2026-11-11 01:31:50 | 1,581 | `e4f0051bb5a886d37c66f0907282a61bc25fee03394201affe66e88fc5e8ddad` |
| NeuroMorpho records citing PMID 40140607 (none; HTTP 404 = no match) | `https://neuromorpho.org/api/neuron/select?q=reference_pmid:40140607&size=1` | 2026-11-11 01:31:41 | 104 | `bcc4b8b6c1a03f3d1a786308f9b808ab1fbe83617a7f38497cfa8aa3667a5acf` |
| NeuroMorpho records citing the article DOI (none) | `https://neuromorpho.org/api/neuron/select?q=reference_doi:10.1038/s41593-025-01919-0&size=1` | 2026-11-11 01:31:41 | 104 | `bcc4b8b6c1a03f3d1a786308f9b808ab1fbe83617a7f38497cfa8aa3667a5acf` |
| NeuroMorpho positive control for the same query form (23 records) | `https://neuromorpho.org/api/neuron/select?q=reference_pmid:30455464&size=1` | 2026-11-11 01:32:10 | 1,618 | `aad098a6604ceb14d3f887951ac32721bd35506d608058c2c8b2f79bb6667c1b` |
| CC BY-NC-ND 4.0 legal code | `https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode.en` | 2026-11-11 01:37:32 | 49,896 | `e4acba200a8b8139aebcffffca4a7fcc41b357d7a5acbe9e639b6e960cbb97cb` |

Other sources.

- Published paper, Methods and Reporting Summary (pages 17 to 20 and 32 to 34 of a 34-page publisher PDF). The article
  is behind the publisher's paywall and Europe PMC has no open copy. The text was read from a copy on a third-party
  laboratory web page:
  `https://tan-lab-siat.github.io/assets/pdf/Zhuolei Jiao,Xiaohong Xu et al., nature neuroscience 2025 Projectome-based characterization of hypothalamic peptidergic neurons in male mice, Nature Neuroscience 2025.pdf`,
  8,375,264 bytes, SHA-256 `0e88ad7e885f4d2bcfd9241404f4b9afc1be25f20fea440f941703d76c8ba57a`. It was not compared
  with the publisher's copy. The fetch tool saved the PDF outside the repository, at
  `/Users/andresarriaga/.claude/projects/-Users-andresarriaga-Documents-SYN4RCH/6fc8191b-130f-4b76-b6a8-69340161f99e/tool-results/webfetch-1794360214498-e54wtk.pdf`;
  it is copyrighted and should stay private. Only the Methods and Reporting Summary pages were kept as text. A keyword
  search also showed a few sentences of the Results; none is used here except the remark on spinal cord truncation
  and the "registered users" wording.
- Local documents: `reports/swc_access.md`, `reports/swc_validation.md`, `reports/swc_prototype.md`,
  `reports/d109_jiao_conventions.md`, `reports/d109_jiao_audit.md`, `reports/d109_protocol.md` (section 16),
  `reports/conversion_map.md`, `data/raw/jiao2025/swc/MANIFEST.tsv`.
- Not reachable: `http://bap.cebsit.ac.cn/neuronbrowser/api/user/selectNeurons?id=192710` (connection refused);
  `https://mouse.braindatacenter.cn/` (expired certificate).
- Queries that returned nothing and are not evidence: the NeuroMorpho literature endpoint returned zero rows for every
  PMID tried, including a control, so only the neuron endpoint (with its control) is used.

### Appendix: the 166 animals

Status: `new` = no trace of the animal was ever inspected; `inspected` = at least one trace inspected, not in the
development set; `development` = one of the 14 development animals. Zip bytes are exact (BSDC record). Portal bytes
are the E2 estimate (1.35 times the zip bytes), not a measurement.

| animal | source label | neurons | Left | Right | Axon_only | status | BSDC zip bytes (exact) | portal bytes (estimate E2) |
|---|---|---:|---:|---:|---:|---|---:|---:|
| 192710 | Pomc | 38 | 37 | 1 | 2 | inspected | 23,594,115 | 31,865,373 |
| 192711 | Pomc | 20 | 13 | 7 | 0 | new | 12,444,164 | 16,806,645 |
| 192712 | Pomc | 23 | 14 | 9 | 2 | new | 11,887,298 | 16,054,562 |
| 192713 | Pomc | 23 | 13 | 10 | 23 | new | 19,497,792 | 26,333,024 |
| 192714 | Pomc | 11 | 5 | 6 | 8 | new | 9,337,842 | 12,611,357 |
| 195037 | Oxt | 58 | 0 | 58 | 53 | new (also in Li 2024) | 17,202,219 | 23,232,705 |
| 195038 | Oxt | 48 | 22 | 26 | 47 | inspected (also in Li 2024) | 9,671,750 | 13,062,322 |
| 195039 | Oxt | 62 | 1 | 61 | 23 | inspected | 17,990,544 | 24,297,389 |
| 201221 | Nts | 33 | 33 | 0 | 0 | inspected | 43,108,652 | 58,221,012 |
| 201223 | Nts | 86 | 40 | 46 | 0 | new | 103,850,382 | 140,256,633 |
| 201224 | Nts | 110 | 48 | 62 | 0 | new | 121,709,427 | 164,376,424 |
| 201225 | Nts | 92 | 92 | 0 | 1 | new | 115,615,958 | 156,146,801 |
| 201226 | Nts | 76 | 0 | 76 | 0 | new | 81,025,055 | 109,429,558 |
| 201227 | Nts | 37 | 0 | 37 | 0 | new | 41,423,457 | 55,945,048 |
| 201228 | Nts | 51 | 0 | 51 | 51 | new | 54,096,137 | 73,060,319 |
| 201229 | Tac1 | 118 | 36 | 82 | 0 | new | 109,574,143 | 147,986,941 |
| 201230 | Tac1 | 68 | 36 | 32 | 0 | development | 73,747,532 | 99,600,794 |
| 201232 | Tac1 | 54 | 24 | 30 | 0 | development | 61,188,318 | 82,638,766 |
| 201578 | Nts | 50 | 0 | 50 | 0 | new | 39,862,840 | 53,837,334 |
| 201580 | Nts | 81 | 0 | 81 | 1 | new | 115,334,088 | 155,766,118 |
| 201789 | Pdyn | 47 | 31 | 16 | 0 | new | 26,140,674 | 35,304,665 |
| 201793 | Pdyn | 56 | 43 | 13 | 2 | new | 49,595,158 | 66,981,457 |
| 201795 | Pdyn | 36 | 31 | 5 | 1 | new | 24,171,071 | 32,644,589 |
| 201796 | Tac1 | 70 | 70 | 0 | 0 | new | 102,726,529 | 138,738,797 |
| 201797 | Tac1 | 76 | 0 | 76 | 1 | new | 97,310,012 | 131,423,442 |
| 201892 | Adcyap1 | 20 | 0 | 20 | 0 | new | 20,154,606 | 27,220,094 |
| 201893 | Adcyap1 | 59 | 0 | 59 | 1 | new | 50,806,942 | 68,618,049 |
| 201894 | Adcyap1 | 44 | 0 | 44 | 2 | development | 55,884,214 | 75,475,232 |
| 201895 | Tac2 | 26 | 26 | 0 | 0 | new | 24,187,024 | 32,666,134 |
| 201896 | Tac2 | 13 | 13 | 0 | 0 | new | 11,284,434 | 15,240,355 |
| 201922 | Nts | 60 | 0 | 60 | 0 | inspected | 63,233,576 | 85,401,019 |
| 201924 | Nts | 135 | 0 | 135 | 3 | new | 139,313,803 | 188,152,268 |
| 201925 | Nts | 69 | 2 | 67 | 0 | new | 71,168,204 | 96,117,245 |
| 201926 | Penk | 19 | 0 | 19 | 0 | new | 11,160,131 | 15,072,476 |
| 201929 | Sst | 30 | 28 | 2 | 0 | new | 26,819,602 | 36,221,601 |
| 202166 | Penk | 45 | 45 | 0 | 0 | new | 22,391,940 | 30,241,758 |
| 202167 | Penk | 10 | 0 | 10 | 0 | new | 5,809,249 | 7,845,765 |
| 202168 | Penk | 24 | 0 | 24 | 0 | new | 31,653,156 | 42,749,627 |
| 202170 | Penk | 71 | 0 | 71 | 0 | development | 54,065,799 | 73,019,345 |
| 202171 | Penk | 84 | 0 | 84 | 0 | new | 78,957,744 | 106,637,521 |
| 202172 | Sst | 11 | 11 | 0 | 0 | new | 7,868,564 | 10,627,003 |
| 202174 | Sst | 47 | 47 | 0 | 0 | new | 37,335,572 | 50,424,096 |
| 202175 | Sst | 9 | 6 | 3 | 0 | new | 4,841,898 | 6,539,295 |
| 202275 | Tac1 | 23 | 0 | 23 | 0 | new | 33,407,691 | 45,119,239 |
| 202276 | Tac1 | 10 | 0 | 10 | 0 | new | 12,269,510 | 16,570,764 |
| 202455 | Penk | 54 | 15 | 39 | 1 | new | 71,911,971 | 97,121,751 |
| 202456 | Penk | 25 | 0 | 25 | 2 | new | 38,133,795 | 51,502,147 |
| 202457 | Pomc | 53 | 35 | 18 | 0 | inspected | 46,854,028 | 63,279,384 |
| 202458 | Pomc | 51 | 46 | 5 | 0 | new | 47,871,066 | 64,652,959 |
| 202459 | Pomc | 45 | 35 | 10 | 0 | new | 31,570,175 | 42,637,556 |
| 202460 | Pomc | 40 | 14 | 26 | 1 | new | 31,360,063 | 42,353,786 |
| 202557 | Sst | 33 | 30 | 3 | 0 | new | 35,455,930 | 47,885,518 |
| 202558 | Sst | 36 | 35 | 1 | 0 | new | 37,852,901 | 51,122,782 |
| 202559 | Tac2 | 45 | 0 | 45 | 0 | development | 49,978,425 | 67,499,083 |
| 202560 | Tac2 | 30 | 0 | 30 | 0 | new | 34,722,695 | 46,895,237 |
| 202561 | Tac2 | 41 | 0 | 41 | 0 | new | 49,209,712 | 66,460,887 |
| 202562 | Tac2 | 73 | 73 | 0 | 0 | new | 66,363,798 | 89,628,585 |
| 202564 | Tac2 | 18 | 18 | 0 | 0 | new | 23,233,627 | 31,378,510 |
| 202570 | Penk | 14 | 5 | 9 | 0 | new | 18,664,951 | 25,208,219 |
| 202571 | Pdyn | 140 | 38 | 102 | 1 | new | 150,210,256 | 202,868,630 |
| 202572 | Pdyn | 87 | 21 | 66 | 0 | new | 93,713,310 | 126,565,864 |
| 202600 | Oxt | 5 | 2 | 3 | 0 | new | 4,396,799 | 5,938,160 |
| 202601 | Oxt | 110 | 54 | 56 | 0 | new | 32,722,147 | 44,193,368 |
| 202602 | Oxt | 28 | 0 | 28 | 1 | new | 22,566,361 | 30,477,325 |
| 202603 | Oxt | 21 | 0 | 21 | 5 | new | 22,216,956 | 30,005,431 |
| 202614 | Avp | 3 | 1 | 2 | 0 | new | 1,010,189 | 1,364,325 |
| 202616 | Sst | 29 | 18 | 11 | 1 | new | 47,535,219 | 64,199,376 |
| 202617 | Sst | 109 | 51 | 58 | 0 | development | 98,223,108 | 132,656,636 |
| 202618 | Sst | 68 | 46 | 22 | 1 | inspected | 63,228,169 | 85,393,716 |
| 202619 | Sst | 110 | 47 | 63 | 0 | new | 98,961,255 | 133,653,552 |
| 202620 | Sst | 74 | 45 | 29 | 0 | development | 70,767,989 | 95,576,729 |
| 202662 | Orexin | 2 | 0 | 2 | 0 | new | 380,730 | 514,200 |
| 202792 | Pdyn | 22 | 0 | 22 | 0 | new | 19,965,205 | 26,964,296 |
| 202793 | Pdyn | 15 | 9 | 6 | 0 | new | 13,743,965 | 18,562,110 |
| 202812 | Agrp | 28 | 17 | 11 | 3 | new | 14,415,787 | 19,469,449 |
| 202813 | Agrp | 32 | 11 | 21 | 1 | inspected | 13,861,264 | 18,720,530 |
| 202814 | Agrp | 5 | 4 | 1 | 2 | new | 2,819,613 | 3,808,069 |
| 202815 | Agrp | 27 | 9 | 18 | 1 | new | 16,662,955 | 22,504,394 |
| 202816 | Agrp | 10 | 3 | 7 | 1 | new | 4,926,588 | 6,653,675 |
| 202817 | Oxt | 2 | 0 | 2 | 0 | new | 599,741 | 809,989 |
| 202818 | Oxt | 2 | 0 | 2 | 0 | new | 2,307,943 | 3,117,026 |
| 202820 | Oxt | 3 | 0 | 3 | 0 | new | 1,144,409 | 1,545,598 |
| 202822 | Vip | 3 | 0 | 3 | 1 | new | 2,817,383 | 3,805,057 |
| 202824 | Vip | 13 | 0 | 13 | 0 | development | 12,325,923 | 16,646,953 |
| 202825 | Avp | 10 | 2 | 8 | 0 | inspected | 5,626,970 | 7,599,586 |
| 202826 | Avp | 28 | 14 | 14 | 9 | new | 8,575,659 | 11,581,980 |
| 210004 | Penk | 1 | 0 | 1 | 0 | new | 1,154,575 | 1,559,328 |
| 210006 | Penk | 37 | 14 | 23 | 1 | new | 35,547,739 | 48,009,512 |
| 210007 | Pmch | 93 | 43 | 50 | 0 | development | 88,110,802 | 118,999,316 |
| 210016 | Penk | 48 | 47 | 1 | 1 | new | 44,706,964 | 60,379,636 |
| 210073 | Penk | 42 | 42 | 0 | 0 | inspected | 54,404,178 | 73,476,348 |
| 210074 | Penk | 20 | 20 | 0 | 0 | new | 17,486,882 | 23,617,161 |
| 210104 | Orexin | 29 | 0 | 29 | 0 | development | 39,524,540 | 53,380,438 |
| 210105 | Orexin | 12 | 12 | 0 | 0 | inspected | 7,793,167 | 10,525,174 |
| 210245 | Oxt | 9 | 0 | 9 | 0 | new | 4,792,575 | 6,472,681 |
| 210247 | Pdyn | 42 | 42 | 0 | 0 | new | 42,910,264 | 57,953,077 |
| 210660 | Penk | 23 | 0 | 23 | 2 | new | 35,204,226 | 47,545,576 |
| 210672 | Adcyap1 | 19 | 16 | 3 | 0 | new | 19,221,230 | 25,959,510 |
| 210673 | Adcyap1 | 39 | 35 | 4 | 0 | new | 38,470,583 | 51,957,001 |
| 210674 | Adcyap1 | 13 | 13 | 0 | 1 | new | 11,015,395 | 14,877,001 |
| 210757 | Pdyn | 24 | 21 | 3 | 0 | new | 15,844,569 | 21,399,111 |
| 210758 | Pdyn | 4 | 3 | 1 | 0 | new | 2,914,827 | 3,936,662 |
| 210796 | Pdyn | 3 | 3 | 0 | 0 | development | 949,450 | 1,282,293 |
| 210797 | Pdyn | 62 | 62 | 0 | 1 | new | 51,654,586 | 69,762,847 |
| 210951 | Pmch | 10 | 0 | 10 | 0 | inspected | 6,680,386 | 9,022,292 |
| 210952 | Pmch | 4 | 0 | 4 | 0 | new | 2,135,725 | 2,884,434 |
| 210953 | Pmch | 9 | 0 | 9 | 0 | new | 9,586,934 | 12,947,772 |
| 210975 | Tac1 | 66 | 66 | 0 | 0 | new | 75,337,892 | 101,748,678 |
| 210989 | Orexin | 26 | 0 | 26 | 1 | new | 27,567,620 | 37,231,847 |
| 210990 | Orexin | 15 | 0 | 15 | 10 | new | 14,641,227 | 19,773,920 |
| 210991 | Orexin | 19 | 0 | 19 | 13 | new | 31,560,939 | 42,625,082 |
| 210992 | Adcyap1 | 44 | 0 | 44 | 1 | new | 50,453,136 | 68,140,211 |
| 210993 | Adcyap1 | 35 | 0 | 35 | 0 | new | 33,976,484 | 45,887,431 |
| 210994 | Adcyap1 | 34 | 0 | 34 | 0 | new | 33,846,702 | 45,712,152 |
| 210995 | Adcyap1 | 48 | 42 | 6 | 5 | new | 45,744,441 | 61,780,815 |
| 210996 | Adcyap1 | 10 | 10 | 0 | 0 | new | 6,400,583 | 8,644,400 |
| 210997 | Adcyap1 | 52 | 44 | 8 | 0 | new | 56,835,390 | 76,759,856 |
| 211028 | Sst | 38 | 3 | 35 | 0 | new | 30,656,809 | 41,403,996 |
| 211266 | Pdyn | 55 | 0 | 55 | 3 | inspected | 56,871,819 | 76,809,056 |
| 211267 | Pdyn | 38 | 27 | 11 | 0 | new | 25,304,860 | 34,175,844 |
| 211268 | Pdyn | 63 | 31 | 32 | 1 | new | 39,298,606 | 53,075,300 |
| 211269 | Tac2 | 72 | 71 | 1 | 0 | new | 63,355,037 | 85,565,060 |
| 211270 | Tac2 | 27 | 27 | 0 | 0 | new | 15,901,868 | 21,476,497 |
| 211271 | Adcyap1 | 38 | 5 | 33 | 0 | new | 22,944,261 | 30,987,703 |
| 211272 | Adcyap1 | 22 | 2 | 20 | 0 | new | 20,378,614 | 27,522,631 |
| 211273 | Adcyap1 | 31 | 15 | 16 | 0 | new | 15,353,159 | 20,735,431 |
| 211275 | Adcyap1 | 49 | 0 | 49 | 1 | new | 45,206,314 | 61,054,040 |
| 211276 | Adcyap1 | 84 | 0 | 84 | 0 | new | 68,391,180 | 92,366,695 |
| 211353 | Adcyap1 | 121 | 121 | 0 | 1 | development | 99,595,067 | 134,509,556 |
| 211354 | Adcyap1 | 81 | 80 | 1 | 0 | new | 70,212,098 | 94,825,963 |
| 211355 | Adcyap1 | 48 | 38 | 10 | 1 | new | 51,430,646 | 69,460,401 |
| 211392 | Adcyap1 | 31 | 31 | 0 | 0 | new | 33,300,265 | 44,974,154 |
| 211393 | Adcyap1 | 126 | 126 | 0 | 0 | new | 96,312,981 | 130,076,887 |
| 211458 | Tac1 | 41 | 2 | 39 | 0 | new | 25,348,727 | 34,235,089 |
| 211460 | Nts | 86 | 0 | 86 | 0 | new | 97,250,327 | 131,342,833 |
| 211461 | Nts | 47 | 3 | 44 | 3 | inspected | 41,790,845 | 56,441,229 |
| 211462 | Nts | 13 | 2 | 11 | 1 | new | 9,498,340 | 12,828,120 |
| 211463 | Nts | 13 | 0 | 13 | 0 | new | 12,995,232 | 17,550,898 |
| 211464 | Tac2 | 33 | 32 | 1 | 0 | new | 41,338,902 | 55,830,851 |
| 211465 | Tac2 | 27 | 27 | 0 | 0 | inspected | 45,109,551 | 60,923,355 |
| 211491 | Tac2 | 7 | 1 | 6 | 0 | new | 9,915,331 | 13,391,293 |
| 211493 | Tac2 | 111 | 57 | 54 | 0 | new | 128,152,981 | 173,078,859 |
| 211494 | Penk | 59 | 19 | 40 | 0 | new | 62,286,305 | 84,121,668 |
| 211525 | Penk | 58 | 21 | 37 | 2 | new | 55,702,889 | 75,230,341 |
| 211526 | Penk | 15 | 8 | 7 | 0 | new | 15,179,929 | 20,501,472 |
| 211551 | Sst | 72 | 0 | 72 | 2 | new | 61,941,419 | 83,655,878 |
| 211552 | Sst | 64 | 11 | 53 | 0 | new | 56,585,742 | 76,422,691 |
| 211553 | Sst | 125 | 0 | 125 | 0 | new | 85,041,220 | 114,853,647 |
| 211554 | Sst | 124 | 27 | 97 | 0 | new | 101,155,927 | 136,617,598 |
| 211555 | Sst | 59 | 30 | 29 | 1 | new | 58,673,105 | 79,241,809 |
| 211556 | Tac2 | 2 | 0 | 2 | 1 | new | 2,066,230 | 2,790,577 |
| 211557 | Tac2 | 23 | 15 | 8 | 0 | new | 23,912,941 | 32,295,968 |
| 211558 | Tac2 | 91 | 66 | 25 | 0 | new | 113,837,208 | 153,744,485 |
| 211559 | Tac2 | 24 | 14 | 10 | 0 | new | 32,464,912 | 43,845,956 |
| 211595 | Sst | 72 | 26 | 46 | 0 | inspected | 53,249,473 | 71,916,844 |
| 211596 | Sst | 42 | 16 | 26 | 0 | new | 33,691,860 | 45,503,028 |
| 211597 | Sst | 20 | 0 | 20 | 0 | new | 13,767,249 | 18,593,557 |
| 211600 | Sst | 9 | 0 | 9 | 0 | new | 6,729,756 | 9,088,969 |
| 211601 | Sst | 14 | 14 | 0 | 0 | new | 8,141,590 | 10,995,742 |
| 211889 | Pdyn | 121 | 43 | 78 | 0 | development | 98,543,556 | 133,089,422 |
| 211890 | Pdyn | 40 | 11 | 29 | 0 | inspected | 47,189,944 | 63,733,060 |
| 212201 | Trh | 57 | 34 | 23 | 0 | new | 41,664,998 | 56,271,265 |
| 212202 | Trh | 42 | 15 | 27 | 0 | development | 28,977,766 | 39,136,340 |
| 212247 | Crh | 13 | 0 | 13 | 1 | new | 3,119,348 | 4,212,880 |
| 212248 | Crh | 16 | 0 | 16 | 2 | inspected | 5,490,346 | 7,415,066 |
| 220102 | Crh | 47 | 0 | 47 | 1 | new | 24,430,575 | 32,995,066 |
