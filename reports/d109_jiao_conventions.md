# D-109: conventions of the Jiao 2025 tables, from source documentation

Written 2026-11-10, before any model is fitted to real outcomes. Companion to `reports/d109_jiao_audit.md`
(input audit) and `reports/d109_protocol.md` (protocol). **No statistic relating an input to the axon-length
outcome was computed for this file.** Every point is marked RESOLVED (with the document that settles it),
INFERRED / INDICATED (evidence that does not settle it for the file actually used) or UNRESOLVED. Table ranges
were never used as proof of a convention. Revised 2026-11-10 after independent review: several points first
marked RESOLVED were relabelled because the evidence concerned sibling tables or was not saved.

Two disclosures. (1) While previewing Nature Supplementary Tables 3, 4 and 6 the verifier displayed their first
rows; they are the paper's own published summaries, nothing was computed from them and they are not used.
(2) For the Li 2024 check the first 500 bytes (the soma line) of 6 Li SWC files were read by range request;
no axon data and no bulk download.

## 1. Sources, with full URLs and SHA-256

The data tables themselves (fetched by `analysis/fetch_jiao.py`) are listed with URLs, zip member paths and
hashes in `reports/d109_jiao_audit.md`. The documentation below is fetched and verified by
`python3 -I analysis/d109_jiao.py sources --fetch [dir]` (URLs and hashes are in the `SOURCES` table of that
script); `python3 -I analysis/d109_jiao.py sources` verifies the local copies in `data/interim/d109_docs/`
(not version controlled). A fresh fetch on 2026-11-10 reproduced all 19 hashes.

Zenodo record 14588094 (doi:10.5281/zenodo.14588094, CC BY 4.0), archive
`https://zenodo.org/records/14588094/files/Code%20for%20Hypothalamic%20projectome.zip` (918,302,741 bytes),
members read by HTTP range request:

| zip member | SHA-256 |
|---|---|
| `Main Figs/fig_for_plot_neurons/fig_for_plot_neurons.ipynb` | f0c5b0ec3a8edc18607a79d8e65db9714e3447b1b21c57506b5538a6d16f83b8 |
| `Main Figs/fig5b&5f&Extended Data Fig8a/fig5b&5f&Extended Data Fig8a.ipynb` | e112254b81a48ad0a42b62046f35d5e7a82f69ad45d83cf4f2ffea0ba3798b16 |
| `Main Figs/fig6a/Fig 6a.ipynb` | 0108a2881ec24b7e8636eee5f266da3e6f275380607fb9dc922c36ba4019a357 |
| `Main Figs/fig1b/brain_areas.xlsx` | 38c22fbac2df9c29136fcd0de2fa7a9198c2ccd4126df10a93b3bc1f4765da3a |
| `Main Figs/fig1b/fig 1b.ipynb` | 99f9aca8f60943e371f2e6e486bf12502c3d8e68d6c1f0c2765bfac19e35cb09 |
| `Main Figs/fig1b/visualization.py` | 77c7e74bb1d347e58e1e2550bd695c7ce0a9cf1051f4b1c43b1e19066058d6fa |
| `Main Figs/fig1f&Supplementary fig3/fig_1f&Supplementary_fig3.R` | 52be573fe531f4ab208e25d37de3c89822137ce2a7a6e5a3bd1cf6450550ad16 |
| `Extended Data Fig/Extended Data Fig10a/Extended Data Fig10a.ipynb` | 28e8a8edb9fb657cd4e0afb9323fd9287d569bc6194d855a6c87d40651f15e67 |
| `Main Figs/fig5i-j/fig 5i-j.ipynb` | f092342bc8ed915b5c285d424b48edb01377d1db388e85a7bc149e4885294e27 |
| `Main Figs/fig1f&Supplementary fig3/brain_area_sort.csv` | 96fac9c9079d10cf73e77c00a8f090ea7fd08a9135e0a74ac4b9e642539e7fea |
| `Extended Data Fig/Extended Data Fig2b/Allen_volumne.csv` | a7801789f221d1d87c75893b546a2ac81ce017bcfe1a03923b882014c8ce3df8 |

Also read from the archive: every notebook, `.py` and `.R` file, and the header row only of
`Extended Data Fig/Extended Data Fig4a&fig5a/All_neurons_projection_infor&type.csv` (7,963,560 bytes; the
"385-column table"). The member listing is saved as `data/interim/d109_docs/zenodo/zip_listing.txt`
(sha256 a417e659a5eb10dc604244d8c779e1efbce4022951880a23c66468545da5d465).

pyswcloader (`https://github.com/txgxxx/pyswcloader`, commit bc64d464976dd4fad88619b64d91acd5999526fd, no
licence file; its README cites the paper), files at
`https://raw.githubusercontent.com/txgxxx/pyswcloader/bc64d464976dd4fad88619b64d91acd5999526fd/pyswcloader/<file>`:

| file | SHA-256 |
|---|---|
| swc.py | 5bbd41f890979edfd809dbd91224448fea00df5be5144d0350d10df17ee41c22 |
| projection.py | 814107e57600324f4c11a8f3d0f9de60054dc158c4c4b8bc23951276fb740a80 |
| brain.py | f1d324a41549e6789845b3cc2908b7d69a6b5fd2b68a02cc63d742699323ce11 |

Nature supplementary files (free to read; no licence stated):

| file | URL | SHA-256 |
|---|---|---|
| Reporting Summary | `https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-025-01919-0/MediaObjects/41593_2025_1919_MOESM2_ESM.pdf` | 2dd3d2ec841601dd6b0f6b4fd73e4230d8b702fabeb04e01707082334143012d |
| Supplementary Tables 1-6 | `https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-025-01919-0/MediaObjects/41593_2025_1919_MOESM3_ESM.xlsx` | 893d614b9490c83c3469b176594526716d11055db0461b804b5341ba37b4359b |

Portal index files (base `https://mouse.digital-brain.cn/projectome/2/srv/`; no licence stated):

| path | bytes | Last-Modified | SHA-256 |
|---|---|---|---|
| `info/mouse/hy/mouse.neuron.info.json` | 1,980,365 | 2023-12-07 | fdf1fa2eff099bc755e5218dc6437e772f12fe0a934c01a4072d846d9b3f4da3 |
| `info/mouse/pvh_oxt/mouse.neuron.info.json` | 52,388 | 2024-06-18 | 61203bb4056330349c260cf0668035a7a9b043c42561113fd8f8b89785d55b86 |
| `info/mouse/pvh_oxt/sampleInfo/PVH_Oxt.json` | 7,418 | 2024-06-17 | a15fec89f5a8a85062c4e4967e53346509eddc9f15b8fe0e78e12999521562ed |
| `api/v1/public/projects` (dynamic, may change) | 7,821 | - | 7154a66174ea955c05766310ea07f599260d48a87cc68eac1ca4177788680b82 |

Read in a browser, not hashable: bioRxiv 10.1101/2023.05.25.542241v1 full-text Methods; BSDC dataset page
`https://www.braindatacenter.cn/datacenter/web/#/dataSet/details?id=1800731789745954817`
(doi:10.12412/BSDC.1667278800.20003) and the site's Terms of use (`#/termsForUsage`); portal main view
`https://mouse.digital-brain.cn/projectome/hy`.

**Not read.** The final Nature Methods (paywalled; Europe PMC lists no open copy) and the Supplementary
Information PDF (MOESM1, 20.7 MB). The bioRxiv quotations below came through a summarising fetch tool and
**have not been re-checked word for word**; the points that rest on them are marked.

## 2. Axes and units of AP / DV / ML: RESOLVED for sibling tables, INFERRED for the input table

AP = CCFv3 x, DV = y, ML = z, in micrometres, in the 10 um `ccf_2017` frame (13,200 x 8,000 x 11,400 um).

* Sibling tables (same column names, other files): `fig5b&5f&Extended Data Fig8a.ipynb` loads
  `ccf_2017/annotation_10.nrrd`, takes the coronal slice `Image_slice` computed from the AP values divided by 10,
  and overlays `ax.scatter(dfplot.ML/10,dfplot.DV/10,...)`; `Fig 6a.ipynb` does the same. These notebooks plot
  `orexin_soma_infor_type.xlsx` and `Agrp_Pomc_ARH_soma_infor.xlsx`, **not** the file D-109 uses.
* The input table itself (`All_soma_infor_HY_neurons.csv`): `Main Figs/fig1b/fig 1b.ipynb` reads it and passes
  `df[["AP","DV","ML"]].values` to `generate_go_scatters`, which plots columns 0, 1, 2 as x, y, z, unscaled, in
  the same 3D scene as the authors' structure meshes (`visualization.py`, `generate_go_ba`, mesh vertices used
  as they are). So the three columns are in the coordinate frame of those meshes, in the order x, y, z. That the
  meshes are CCFv3 meshes in micrometres is not stated in the archive. **Inferred, not resolved, for this file.**
* pyswcloader `swc.py`: `swc_preprocess(..., flip=True, dimension=[13200, 8000, 11400])`; `brain.py`:
  `annotataion[int(coords[0]/resolution), int(coords[1]/resolution), int(coords[2]/resolution)]`.

The D-109 models use the three columns as given (trees on raw values, the others on standardised values), so
no modelling step depends on the unit or the axis names.

UNRESOLVED, and without effect on modelling: which anatomical side z < 5,700 is (pyswcloader calls it "the left
side"; the mesh folder `fig 1b` uses is named `HY_obj_right_part`; Allen's convention was not re-read). The code
uses AP 5,472 um as its bregma reference; not otherwise verified.

## 3. Mirroring: rule RESOLVED; that the soma table is mirrored is INFERRED; original hemisphere NOT in the open tables

* Rule (authors' notebook `fig_for_plot_neurons.ipynb`, with the note "we recommend you mirrow the neurons in
  the same brain side using the following code"): if the first SWC node has z > 5700, every node gets
  `records[4] = 11400 - float(records[4])`; otherwise the file is copied unchanged. pyswcloader `swc.py`:
  `if flip==True and float(data.loc[1, 'z'])>(11400/2)` then `data.z = dimension[2] - data.z`.
* The whole neuron is reflected (soma and axon), so ipsilateral / contralateral relations are kept.
* bioRxiv Methods, "Neuron exclusion", reports that all neurons were mirrored to the same hemisphere before
  analysis (quotation not re-checked word for word).
* **Inference, not a statement about the file:** that the published soma table is the mirrored one. It rests on
  the documented rule plus the fact that all 7,180 ML values are below 5,700; the mirroring notebook reads the
  soma CSV only for IDs and mirrors raw SWC files, and says nothing about the table. Supporting evidence, not
  documentation: Li 2024's SWC soma lines for the two shared brains carry ML on both sides of 5,700.
* Original hemisphere per neuron: absent from every Zenodo table read. It exists only in the BSDC SWC files
  (the notebook reads raw `<brain>\swc_allen_space\<n>.swc` and mirrors them itself) and in the portal index
  (Left / Right). **Bilateral soma coordinates therefore cannot be preserved from the current inputs.** The
  protocol uses the published (mirrored) coordinates and says so; laterality is Stage 2 material.

## 4. The axon-length table

* Units: **micrometres INDICATED, not resolved.** The table carries no unit. Evidence, none of it about this
  file: (i) pyswcloader `projection.py` computes lengths as Euclidean node-to-parent distances (`math.dist`) in
  SWC coordinates, which are micrometres in the frame above, but which tool made this table is unresolved (next
  point); (ii) the authors label axes "Total axon length(um)" in `Extended Data Fig10a.ipynb` and
  "Axon length(um)" in `fig 5i-j.ipynb`, which plot other tables (`non_hypothalamic_area_top10.xlsx`,
  `orexin_data_for_violin_plot.xlsx`); (iii) `fig_1f&Supplementary_fig3.R`, the only code that reads this
  table, draws it with a colour scale from 0 to 2,000 and no unit. The protocol's `log10(1 + L)` assumes
  micrometres and says so.
* Per-region rule in pyswcloader `projection_length`: each node takes the `annotation_10` voxel label; a segment
  inside one label adds its whole length; a segment crossing two labels is split half and half; dendrite nodes
  (SWC type 3 or 4) and label 0 are skipped.
* UNRESOLVED: which tool made `all_neurons_projection_length_fig1.csv`. bioRxiv Methods name "NeuronView",
  the Nature Reporting Summary lists pyswcloader.
* UNRESOLVED: whether the 49 columns sum both hemispheres. Strongly indicated (the default function sums both
  sides; side-specific columns elsewhere are named `_ipsi` / `_contra`), never stated.
* The 49 columns are all CCF acronyms: 25 leaves (ARH PMv TU TMv MPO PH LHA LPO ME SNr VTA DR IO PARN IRN PRNr
  PCG PA MS SI ACB LH PVT PT SPFm) and 24 parents (MM AHN VMH DMH PVH MPN ZI MRN PAG SCm MY NTS PB BMA BLA BST
  MEA CEA LS MD AM HPF Isocortex OLF). No custom merges. They are a subset: the full table has about 383 region
  columns, fibre tracts and ventricles included. No rule for choosing the 49 is documented; `brain_area_sort.csv`
  only groups them for display.
* Hierarchy: exactly one ancestor / descendant conflict among the 49: **MY is the ancestor of IO, NTS, PARN and
  IRN.** UNRESOLVED: whether the MY column is the whole medulla (it would then count those four twice) or only
  voxels labelled MY directly; the same question for HPF and OLF; how Isocortex was summed. A hint only: the
  385-column header lists MY, HPF and OLF next to their children and has no Isocortex column.
* UNRESOLVED: why 13 neurons are zero in all 49 regions. No source mentions them. Possible, undocumented: their
  axon lies wholly in regions outside the 49.
* **Proposed, not run, needs the supervisor's approval:** an outcome-only comparison of the 49-column table with
  the 385-column table (no inputs involved) would settle the MY / HPF / OLF / Isocortex and all-zero questions.

## 5. Location labels: three vocabularies (inputs only)

* Fine label = suffix of the published ID (`<brain>_<n>_<tag>_<location>`; Supplementary Table 5 legend gives
  this ID format). It equals the `location` column of `all_neuron_soma_infor.xlsx` for all 7,180 neurons
  (checked 2026-11-10). 59 values, all CCF acronyms.
* `location` in `All_soma_infor_HY_neurons.csv`: 50 values. It is the fine label with the CCF subtree of MBO
  pooled (TMv 79, MMl 75, SUM 56, MMm 55, LM 37, MMme 27, MMd 22, MMp 4, TMd 3 = 358) and FF (17) pooled into
  its CCF parent ZI. `analysis/d109_jiao.py` asserts this rule reproduces the column exactly.
* `location` in `geno_location_labels.csv`: the 15 most frequent labels, everything else "others" (2,143).
  "others" is a remainder, not a place: HY 1,316 (soma in voxels labelled hypothalamus with no nucleus), PeF 120,
  PSTN 81, PMd 76, PVp 73, RCH 61, fx 52, AVP 49, SBPV 47, nst 32, STN 32, VMPO 25, AVPV 20, ME 14, ADP 13, PS 12,
  SO 12, PVi 11, VLPO 11, och 11, mtt 10, mfb 9, PVa 7, PD 6, MEPO 6, PST 6, opt 5, sup 5, pm 4, sm 4, SCH 4,
  ASO 3, mtg 2, int 2, PVpo 2. 136 of these somata carry fibre-tract labels.
* UNRESOLVED: how the fine label was assigned (voxel lookup or manual correction, before or after mirroring).
  bioRxiv says the soma location was calculated from the CCFv3 annotation file; Supplementary Table 2 says some
  names follow Paxinos and Franklin or manual demarcation. No code in the archive makes the 15 + others collapse.
* A coordinate-derived nucleus assignment is **not** made in D-109: no annotation volume is in the repository and
  the label procedure is unresolved. The authors' fine label is used and described as such.

## 6. Brains and animals: RESOLVED

* Nature Supplementary Table 1 ("Details of brain samples"): 166 rows, 166 unique "fMOST Mouse ID", equal to the
  166 ID prefixes; per-brain neuron counts match for all 166 (sum 7,180); the marker matches the tag for all
  (the table says "Orexin" where the data say "Wt(orexin)", 6 brains).
* Reporting Summary: all tracing data came from 166 samples; each peptidergic type from at least two mice.
* One row = one mouse ID = one brain. UNRESOLVED only in the strict sense that no sentence says no animal was
  imaged twice.
* Sex: male only (Reporting Summary). Age at injection 53 to 150 days (Supplementary Table 1, all 166 filled);
  53 days is below the "over 8 weeks" of the preprint. Strain: C57BL/6J wild type or Cre lines; the table gives
  the marker per brain, not a background column.
* Input-side covariates in that table, not used by D-109: age, expression time, target region, injection
  coordinates (many bilateral), viruses, titres, volume, QC score. They are line-level experimental covariates
  that the tag also stands for (protocol, section 10).
* Reporting Summary: tracing staff were blind to soma location and peptidergic type.

## 7. Licence of the full SWC set: as displayed 2026-11-10, CONTRADICTS earlier notes

* The BSDC page (title "Single-neuron datasets for mouse hypothalamus", author Xiaofei Wang, released
  2024-06-12, 6.20 GB, 166 files, one zip per brain) displays **"License Agreement CC BY-NC-ND 4.0"**, requires
  naming the source and the author, and suggests the acknowledgement "The dataset is provided by Brain Science
  Data Center, Chinese Academy of Sciences (https://braindatacenter.cn/)."
* `reports/conversion_map.md` line 38 now shows CC BY-NC-ND 4.0 for the BSDC record (corrected by its owner,
  citing `reports/swc_access.md`). **`reports/d109_jiao_audit.md` (lines 25-26) still says the licence field is
  empty, and its lines 41-47 still call units, axes and mirroring unresolved without pointing here; it needs a
  dated pointer to this file** (not edited here: outside this task's files). NoDerivatives matters for
  redistributing anything derived from the SWCs.
* Site Terms of use: data may be browsed, downloaded and used under the owner's licence; no excessive load, no
  bypassing of access restrictions. Whether download needs a login was not tested (no account, no download).
* The portal main view shows no licence text. The Zenodo tables used by D-109 are CC BY 4.0.
* Stage 2 (the SWC set) stays a read-only question of access and terms: no registration, no bulk download.

## 8. Li et al. 2024 (Neuron, doi:10.1016/j.neuron.2023.12.022) and Jiao: overlap verified at the animal level

* **No primary statement of the overlap was found.** Li 2024 and the Jiao Methods are paywalled; the preprint
  and the public Nature text are silent. The "106 are also in Jiao" of `reports/conversion_map.md` line 39 is
  therefore a **derived count, not a published one**: it equals the number of Jiao neurons in the two brains
  both datasets share (58 + 48).
* Portal indexes (hashes in section 1): Li has 264 neurons, all Oxt, from brains 18696 (115), 195037 (65) and
  195038 (84). The Jiao index has 7,180 files whose `<brain>_<n>` stems equal the 7,180 local IDs exactly.
* **Shared animals: brains 195037 and 195038** (Jiao 58 + 48 = 106 neurons; Li 65 + 84 = 149). Brain 18696 is
  not in Jiao; Jiao's 195039 is not in Li.
* Neuron level: matching `<brain>_<n>` between the two indexes gives 68 pairs, not 106, and a matching number is
  **not proven to be the same cell** (soma coordinates differ by roughly 130 to 350 um in AP with a consistent
  sign per brain, i.e. a different registration; hemisphere labels disagree in 10 of 68; soma acronyms agree in
  37). Renumbered duplicates among the unmatched IDs cannot be excluded. UNRESOLVED.
  Number-matched Jiao IDs, brain 195037 (50): 001 002 003 005 006 007 008 009 010 011 012 013 014 015 018 019
  020 021 022 023 024 025 026 027 029 030 031 032 033 034 035 036 037 038 040 043 044 045 046 048 049 050 051 052
  053 054 055 056 058 059. Brain 195038 (18): 001 003 004 005 006 007 008 010 011 012 013 016 017 018 021 024 034
  040. (Full IDs: `data/interim/d109_docs/portal/shared.json`.)
* **Rule.** Any use of Li as a check must drop brains 195037 and 195038 whole, on both sides. Only brain 18696
  is a different animal, and it is not an independent pipeline (same institute, imaging group and portal,
  overlapping authors), one animal, one tag (Oxt), one area (PVH), axon only, with a registration and region
  annotation that visibly differ from Jiao's. **D-109 does not use Li at all** (protocol, section 13); it could
  at most become a small same-laboratory check after its regional lengths are shown comparable.
* Not stated on the portal: sex coding and licence of the Li data.

## 9. Open points carried into the protocol

a. Final Nature Methods not read; preprint quotations not re-checked word for word.
b. Tool that produced the 49-column table.
c. Whether the 49 columns sum both hemispheres.
d. MY / HPF / OLF whole-structure or residual; Isocortex sum; rule for choosing the 49.
e. Why 13 profiles are all zero.
f. Soma-label assignment procedure.
g. Original hemisphere per neuron (Stage 2 only).
h. Anatomical name of the z < 5,700 side.
i. No explicit statement that no animal was imaged twice.
j. Portal terms beyond the main view; whether BSDC download needs a login.
k. Neuron-level identity between Li and Jiao; no primary statement of the overlap.
l. Owner edits still needed outside this task's files: `reports/d109_jiao_audit.md` lines 25-26 (licence) and
   41-47 (conventions) need a dated pointer to this file; `reports/conversion_map.md` line 39 still gives "106
   are also in Jiao" as fact and should read "two shared brains (106 Jiao neurons); neuron identity unresolved".
m. Length unit of the 49-column table: micrometres indicated, not resolved (section 4).
n. Axes and units of the soma table used: inferred from sibling tables and the fig 1b notebook (section 2);
   that the table is the mirrored one: inferred (section 3).
o. Whether the authors' soma label was assigned with knowledge of the reconstruction (section 5).
