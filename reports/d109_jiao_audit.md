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

## Scope

Regional axon length is not a 3D path and not synaptic connectivity. Any tag effect is predictive, not causal.
