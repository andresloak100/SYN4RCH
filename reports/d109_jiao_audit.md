# D-109 input audit: Jiao et al. 2025 region-level tables

Audit of 2026-11-10, before any model is fitted or scored. Fetched by `analysis/fetch_jiao.py` into
`data/raw/jiao2025/` (not version controlled). Nothing below relates an input to the axon outcome.

## Sources

Paper: Jiao Z, Gao T, Wang X, et al. Projectome-based characterization of hypothalamic peptidergic neurons
in male mice. Nat Neurosci 28(5):1073-1088 (2025), doi:10.1038/s41593-025-01919-0. Not open access;
preprint bioRxiv 10.1101/2023.05.25.542241.

Zenodo record 14588094 (doi:10.5281/zenodo.14588094, CC BY 4.0), file
`https://zenodo.org/records/14588094/files/Code%20for%20Hypothalamic%20projectome.zip` (918,302,741 bytes,
md5 d5b6528efb4f1bd22a285e331d9a1645 per the record; members read by range request, archive not downloaded).

| local file | bytes | SHA-256 | zip member |
|---|---|---|---|
| all_neurons_projection_length_fig1.csv | 2,283,305 | 2072bc20550607c195be0a42cee9fcbe0e0ab5aaebfaee0a1016b2a43d35e63c | `Main Figs/fig1f&Supplementary fig3/all_neurons_projection_length_fig1.csv` |
| all_neurons_terminal_terminal_number_Supplementary_fig3_sort.csv | 789,174 | 479f4ba922109737eb3fd0fe92e347d810c67ef4aa9dfc7999db74270386b399 | `Main Figs/fig1f&Supplementary fig3/all_neurons_terminal_terminal_number_Supplementary_fig3_sort.csv` |
| All_soma_infor_HY_neurons.csv | 708,172 | cc01d048bf5a7c9c9f46869fb750f2d67a98ca934903d2e2e71c1a67b1d97cb9 | `Extended Data Fig/Extended Data Fig2a/All_soma_infor_HY_neurons.csv` |
| geno_location_labels.csv | 221,087 | 8260a279767f4d0d2e28b7b426061f1fa0ea17158b65109537738ffbf784c3da | `Main Figs/fig1f&Supplementary fig3/geno_location_labels.csv` |
| all_neuron_soma_infor.xlsx | 521,398 | 7b7ae75b94688619dbb367d5f9e71d8b157615624932375a00e80c796f21f162 | `Main Figs/fig2b&fig2g/all_neuron_soma_infor.xlsx` |
| nature_source_data_fig2.xlsx | 521,398 | 7b7ae75b94688619dbb367d5f9e71d8b157615624932375a00e80c796f21f162 | Nature Source Data Fig. 2: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-025-01919-0/MediaObjects/41593_2025_1919_MOESM5_ESM.xlsx` (no licence stated; byte-identical to the Zenodo xlsx) |

Licence of the full SWC set (not fetched): `reports/conversion_map.md` records that the BSDC licence field is
empty and asks only for attribution; an earlier search's "CC BY-NC-ND 4.0" is recorded there as wrong. Unresolved.

## Findings

* **Neuron IDs.** 7,180 unique IDs (`<brain>_<n>_<tag>_<location>`), identical sets in the projection, soma (csv
  and xlsx) and label tables. The terminal-count table has only a 0..7179 row index and "sort" in its name: its
  row order is unverified and it is not used.
* **Animals.** Brain ID = ID prefix: 166 brains, 1 to 140 neurons each (median 37.5).
* **Tag by animal.** Every brain carries exactly one genetic tag: tag and animal are fully confounded.
  Brains / neurons per tag: Sst 22 / 1,195; Adcyap1 22 / 1,048; Nts 16 / 1,039; Pdyn 17 / 855; Tac2 17 / 663;
  Penk 18 / 649; Tac1 9 / 526; Oxt 11 / 348; Pomc 9 / 304; Pmch 4 / 116; Wt(orexin) 6 / 103; Agrp 5 / 102;
  Trh 2 / 99; Crh 3 / 76; Avp 3 / 41; Vip 2 / 16.
* **Genetic label, not expression.** The tag is a Cre line (or the orexin promoter virus). No neuron has a
  measured gene-expression profile.
* **Soma coordinates.** Columns AP, DV, ML; ranges AP 4,621-8,381, DV 4,659-7,288, ML 3,683-5,699. Consistent
  with micrometres in the 10 um CCFv3 frame (midline 5,700), but **units, axis order and frame are not yet
  confirmed from source documentation.** All 7,180 ML values are below 5,700. This is compatible with somata
  mirrored to one hemisphere, but **does not prove it**; the portal index lists Left 3,071 / Right 4,109.
  UNRESOLVED: whether the tables are mirrored, and by what rule. Bilateral soma coordinates are not present
  in these tables.
* **Outcome table.** Axon length for 49 named regions per neuron (ARH ... Isocortex, OLF), no ipsilateral /
  contralateral split, units not stated in the table (UNRESOLVED). No missing values; 72% zeros; 13 neurons
  are zero in every region; median row total 31,525.
* **Location label.** 16 values; 2,143 neurons are "others". The location in the ID string equals the table's
  label for 64.9% of neurons (the ID keeps finer names). Region ontology of the 49 columns and of the
  location labels is to be fixed in the protocol.
* **Outcome-derived columns.** `cluster` / `type` (31 projectome types) and `Class` (Class1 2,906, Class2
  4,274) are derived from the axons and are excluded from model inputs.
* **Duplicates.** Two neurons share exact soma coordinates with another neuron.
* **Overlap with Li et al. 2024** (Neuron, doi:10.1016/j.neuron.2023.12.022; `reports/conversion_map.md`
  lines 38-39: 264 Oxt neurons from 3 brains, 106 also in Jiao). Not yet verified here. Li must not be used
  as independent validation until the shared neurons are identified and removed.

## Addendum (independent audit, 2026-11-10): three nucleus vocabularies that must not be mixed

Reproduced on the label columns only:

* `all_neuron_soma_infor.xlsx` (= Nature Source Data Fig. 2): 59 labels, the finest. Includes 1,316 `HY`
  (unassigned), fibre-tract labels (fx 52, nst 32, och 11, mtt 10, mfb 9, opt 5, sup 5, pm 4, sm 4, mtg 2,
  int 2) and the fine mammillary / tuberomammillary labels.
* `All_soma_infor_HY_neurons.csv`: 50 labels. Differs from the xlsx on **375 neurons**: 358 collapsed to
  `MBO` (TMv 79, MMl 75, SUM 56, MMm 55, LM 37, MMme 27, MMd 22, MMp 4, TMd 3) and 17 `FF` relabelled `ZI`.
* `geno_location_labels.csv`: 16 labels. Differs from the xlsx on 2,518 neurons, 2,143 of them collapsed
  to `others`.

**Binding requirement for the protocol.** Choose ONE authoritative nucleus vocabulary and document a
crosswalk to the other two and to the CCF ontology (`data/raw/atlas/structure_graph_1.json`); never mix
them silently. Default: keep the fine xlsx label, unless source documentation supports a different
aggregation. Any coordinate-derived nucleus assignment must be verified against a named annotation volume
and registration convention (which CCF annotation version and resolution, which axis order and units)
before it is used; until then it is labelled unverified and is not used as an input or hold-out key.
The protocol must say how `HY` (unassigned) and fibre-tract labels are treated in the unseen-nucleus split.

## Addendum: outside exposure and a parallel test (recorded 2026-11-10)

* A separate agent (Codex) is running its own real-data test of a joint mathematical model against copying
  on the same Jiao tables (`analysis/brm_real.py`, `reports/brm_real_protocol.md`). It reads D-109's S1 split
  file read-only and does not edit D-109 files.
* During that agent's schema inspection, **one raw regional axon-length profile (one neuron's row) was
  displayed**. Reported: no statistic or design choice used its values. Recorded here as outside exposure.
* D-109's protocol choices (metric, splits, baselines, models, controls, readings) must be fixed
  independently of any BRM result. If BRM results become known before D-109 is frozen, D-109's protocol
  must say so and list them under "Already seen".

## Addendum: what later verification changed, and status (2026-11-10, evening)

This supersedes the lines above marked unresolved where it says so. Details and sources:
`reports/d109_jiao_conventions.md`, `reports/swc_access.md`, `reports/swc_validation.md`.

* **Licences.** The BSDC record's details panel shows "License Agreement CC BY-NC-ND 4.0"; the portal shows
  only a copyright line ("All Rights Reserved"); the Zenodo record is CC BY 4.0. The note above that the BSDC
  field is empty is superseded.
* **Soma coordinates.** AP = CCF x, DV = y, ML = z, in micrometres on the 10 um ccf_2017 frame: shown by the
  authors' code for sibling tables, inferred for the table used here. **Mirroring:** the authors reflect every
  neuron whose soma has ML > 5,700 by ML' = 11,400 - ML (their notebook), and this was verified node by node
  on SWC pairs: portal files keep the true hemisphere, Zenodo files are mirrored. That the published soma
  table is the mirrored one is an inference from that rule and the ML range. The original hemisphere exists
  only in the portal index and the SWC files.
* **Outcome table.** Length unit: micrometres indicated, not resolved. Both hemispheres summed: indicated, not
  stated. Of the 49 columns 25 are CCF leaves and 24 are parents; `MY` is an ancestor of `IO`, `NTS`, `PARN`
  and `IRN`, so a double count is possible and unresolved. The 13 all-zero profiles are unexplained.
* **Li et al. 2024.** Verified at the animal level: Li and Jiao share brains 195037 and 195038 (106 Jiao
  neurons, 149 Li neurons); Li's third brain, 18696, is not in Jiao. "106" is that derived count; no published
  statement of the overlap was found, and neuron-level identity is unresolved (68 ID-number matches, not
  proven to be the same cells). Any use of Li must drop both shared brains whole. D-109 does not use Li.
* **Exposure before the freeze.** (i) BRM (`analysis/brm_real.py`, another author, requested separately by the
  user) scored predictions of this same 49-region table on D-109's S1 split on 2026-11-10 (run marker 19:28
  local), after the D-109 design and review fixes were written and before D-109 was frozen. (ii) D-110 and
  D-111 (another author) scored 3D axon geometry on 16 of the 37 downloaded SWC neurons. What the D-109
  session knows of the BRM result: one sentence, that its model had lower held-out error than copying rules
  and simple averages; no number. No D-109 choice used any of this. **D-109 is therefore not the first look
  at these outcomes; if it is run, it is reported as a pre-specified but non-blind check.**
* **Status.** Protocol (`reports/d109_protocol.md`), implementation (`analysis/d109_jiao.py`) and conventions
  (`reports/d109_jiao_conventions.md`) are written and were reviewed three ways (leakage, protocol against
  code, statistics on synthetic data). **Not frozen:** the review fixes need a confirming review, and the
  choice of primary comparison is open (the tree model recovers almost none of a planted moderate tag effect
  that ridge recovers). `run` refuses. No D-109 score exists.
* **Parked (2026-11-10, user instruction).** D-109 is not to be frozen or run. The three reviews and the
  writer's response are preserved in `reports/d109_parked.md`.

## Scope

Regional axon length is not a 3D path and not synaptic connectivity. Any tag effect is predictive, not causal.
