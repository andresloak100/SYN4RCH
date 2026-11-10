# Measured 3D prototype: individual hypothalamic neurons (Jiao et al. 2025)

Built 2026-11-10. Nothing was committed. No account, login or terms dialog was used. Every trace drawn is a measured reconstruction: no arbor was synthesised, completed or added.

Source paper: Jiao et al. 2025, Nat Neurosci, doi:10.1038/s41593-025-01919-0 (7,180 traced neurons, 166 brains).
How the files can be reached publicly and on what terms: `reports/swc_access.md`. Every per-neuron check: `reports/swc_validation.md`.

## What was built

| step | file | output |
|---|---|---|
| fetch | `analysis/swc_fetch.py` | `data/raw/jiao2025/swc/{zenodo,portal,portal_index}/` plus `MANIFEST.tsv` (file, neuron, bytes, full SHA-256, source URL or zip member, frame, licence, how obtained) |
| validate | `analysis/swc_validate.py` | `reports/swc_validation.md`, `data/interim/swc_validation.json` |
| render | `analysis/swc_render.py` | `reports/swc_viewer.html` (8.4 MB, private, untracked), `data/interim/swc_view_dorsal.png`, `data/interim/swc_view_sagittal.png`, `data/interim/swc_render.json` |

All three scripts run with `python3 -I`. They can be re-run: the fetch skips any file that is already present with its manifest hash.

## Sample

- **41 raw SWC files, 46.55 MB, for 37 distinct Jiao neurons.** A further 1.98 MB file is the public portal index.
  - **Zenodo 14588094 (CC BY 4.0):** all 18 Jiao SWCs of Supplementary Fig. 4, read out of the 918 MB zip by HTTP range requests. These files are in the mirrored-left frame. The 18 MouseLight SWCs in the same zip folder were not fetched, so no MouseLight neuron is in any count or view.
  - **Portal `mouse.digital-brain.cn/.../hy_neuron_download/<id>.swc` (copyright line only, "Copyright © 2023 Digital Brain. All Rights Reserved."; the BSDC record of the same data shows CC BY-NC-ND 4.0):** 23 files in the true-hemisphere frame.
    - 20 were new GETs, made one at a time 3 s apart with User-Agent `syn4rch-analysis/0.1`.
    - 3 are byte-identical copies of the user's probe files; their SHA-256 matched the recorded hashes.
  - **4 neurons have both a portal file and a Zenodo file** (195038_041, 201232_056, 202618_143, 211461_026).
- **Coverage:**
  - Peptide lines: all 16. Oxt 2, Tac1 3, Adcyap1 2, Nts 3, Tac2 2, Sst 5, Pdyn 4, Pomc 2, Agrp 2, Avp 2, Crh 1, Orexin 2, Penk 2, Pmch 2, Trh 2, Vip 1.
  - Soma hemisphere (true): 15 Left, 22 Right.
  - Soma labels: 16 (AHN, ARH, DMH, HY, LHA, MMm, MPN, MPO, PH, PSTN, PVH, PVHd, PeF, TU, VMH, ZI).
  - Crh and Vip have only Right-hemisphere neurons in the whole index.
- **Cable:** 748,329 nodes in the 37 drawn files.
  - Type 2 plus type 0 cable: 1,813 mm.
  - Type 3 (dendrite): 102 mm.
  - Lengths are sums of straight segments, in micrometres, in Allen CCFv3 2017.
- **Licence handling:**
  - The Zenodo files may be redistributed with attribution.
  - The portal files stay local (`data/raw/` is gitignored) and are shown only on private pages with attribution.
  - `reports/swc_viewer.html` embeds simplified and reflected portal geometry. It is gitignored, was not published, and is for local viewing only: sharing it, even privately, would distribute all-rights-reserved geometry, and the BSDC CC BY-NC-ND 4.0 terms do not allow adapted versions.
  - `reports/swc_viewer_zenodo_only.html` (`python3 -I analysis/swc_render.py --share`) is the only build meant for sharing. It contains only the 18 Zenodo CC BY 4.0 neurons (drawn from their Zenodo files, also for the 4 paired neurons) and the Allen meshes, with attribution and a statement of the changes (simplified, and Right neurons reflected for the true-hemisphere frame). It was not published; it is gitignored like every `reports/*.html`.

## What the viewer shows (`reports/swc_viewer.html`)

- **Traces:** real traces in CCFv3 micrometres. The scene axes are x = animal's left (−ML), y = dorsal (−DV), z = anterior (−AP). This is right-handed like the data, so the brain is not mirror-imaged (the first build had determinant −1 and showed a mirror image; fixed). An orientation marker (A/P, L/R, D/V) in the top-left corner follows the camera. Each neuron has its own colour. Checkboxes select neurons, grouped by peptide line, with all/none buttons for each line and overall.
- **Frame selector:**
  - *True hemisphere (portal frame)* is the default. This is the measured bilateral frame.
  - *Mirrored-left (Zenodo/paper frame)* is the alternative.
  - Where a neuron is drawn reflected for the chosen frame (ML' = 11400 − ML, applied to the whole neuron at display time), its list entry is marked ⇆ and its provenance panel says so.
  - Default view: 14 Zenodo-only neurons, 9 of them Right, so 9 neurons are drawn reflected back to their true hemisphere. The reflection rule is the authors' own (`fig_for_plot_neurons.ipynb`). It was verified exact on all nodes of the 3 Right-hemisphere pairs; the Left-hemisphere pair (195038_041) is identical in both files.
- **Context:**
  - Root (brain) mesh, transparent.
  - HY mesh.
  - Optionally, 18 nucleus meshes. All of these are drawn on both hemispheres, as stored in CCF.
  - The midline plane at ML = 5,700 µm.
  - Soma spheres (fixed 28 µm radius markers, not measured sizes).
- **Styling by raw SWC type, with nothing relabelled:**
  - Dendrite (type 3): opaque.
  - Axon (type 2): translucent.
  - **Type 0 (unlabelled):** translucent and dashed.
  - Other (types 10, 11 and 12, undocumented): points. The edge into each marker node is drawn in the style of its nearest non-marker ancestor's type, so there are no gaps (the first build left these 236 edges, 548.8 µm, undrawn).
  - Each class has its own toggle.
- **Provenance panel (click a neuron ID):**
  - Peptide line, soma label from the xlsx (and the portal index acronym), hemisphere and reconstruction type.
  - Each source file with its frame, bytes, SHA-256 prefix, URL or zip member, and licence: "CC BY 4.0, Jiao et al. 2025, Zenodo 14588094" or the portal copyright line.
  - Node types, lengths, branch points and tips, the number of nodes drawn, the maximum display displacement, and the result of the portal versus Zenodo comparison.
- **Visible note:** "Measured traces only; no synthesised arbors."
- **Light and dark themes** follow the system setting, with a toggle; both were checked in headless Chrome. Preset views:
  - *dorsal*: from above, anterior at the top, the animal's left on the left of the screen;
  - *left lateral*: from the animal's left, anterior to the left, dorsal up;
  - *frontal (anterior)*: from in front, dorsal up, the animal's left on the right of the screen (the first build put this camera behind the brain; fixed).
  - These were checked in headless Chrome by projecting test points (ML 4,000 vs 7,400; AP 4,000 vs 9,000; DV 3,000 vs 7,000) through each preset camera, and by screenshots: olfactory bulbs at the top in the dorsal view and on the left in the lateral view. The frame switch reflected 9 neurons in the true frame and 13 portal Right neurons in the mirrored frame, as expected.
  - Mouse orbit, pan and zoom are written inline. Clicking, dragging and the provenance panel were not tested by hand.
  - The only external script is three.js r128 from cdnjs, now with a subresource-integrity hash (sha512 from the cdnjs API, `api.cdnjs.com/libraries/three.js/r128?fields=sri`) and `crossorigin="anonymous"`. Headless Chrome loaded it with the hash in place, so the hash matches the served file. Opening the page still makes one request to cdnjs.
- **Static views:** `data/interim/swc_view_dorsal.png` and `data/interim/swc_view_sagittal.png`, made with matplotlib.
  - They use the true-hemisphere frame and draw every node.
  - The plot limits include the descending axons beyond the caudal edge of the atlas.

## Display simplifications (raw files are never changed)

1. **Viewer line decimation.**
   - Rule: within each unbranched run, Ramer–Douglas–Peucker with a 1.0 µm tolerance.
   - Always kept: the soma, every branch point, every tip, and every node where the type changes. Type styling is therefore exact.
   - Result: 228,397 of 748,329 nodes are drawn (30.5%). The largest measured distance of a dropped node from the drawn line is 1.000 µm. The PNGs are not decimated.
2. **Meshes are not decimated.**
3. **No radii.** Lines are 1 px wide, and SWC radii are not used: they range from 0.004 to 56 µm and are not diameters. Somas are fixed 28 µm-radius marker spheres, not measured sizes.
4. **Marker nodes.** The 236 undocumented type 10/11/12 nodes are drawn as points. The 236 edges into them (548.8 µm in all, longest 9.8 µm) are drawn in the style of the nearest non-marker ancestor's type; their raw type is unchanged in the files.
5. **Reflection.** In each frame, neurons whose file is in the other frame are mirrored as whole neurons (ML' = 11400 − ML) and marked ⇆.
6. **One trace per neuron.** For the 4 paired neurons the viewer draws the portal file. In the mirrored frame it reflects that file, which equals the Zenodo file to within 0.000 µm on every node.
7. **Colour carries identity only.** Thirty-seven hues spaced by the golden angle cannot all be told apart. Identity is also given by the list, the swatches and the provenance panel.
8. **PNG styling.** Thick lines are type 3; thin lines are type 2 or type 0 together; an edge into a marker node takes its parent's style, as in the viewer. Context is drawn as scattered mesh vertices.

## Validation summary (full list: `reports/swc_validation.md`, 106 anomaly lines)

- **Tree checks (all 41 files):**
  - Exactly one root, at id 1 with type 1, and one soma node.
  - Sequential ids, all parents present, no cycles, one connected component.
- **Hemisphere:**
  - The soma side matches the portal index for all 23 portal files.
  - All 18 Zenodo somas are at ML < 5,700, as the mirrored frame requires.
  - The frame check on the 4 pairs:
    - The Left neuron is identical.
    - The 3 Right neurons satisfy ML_p + ML_z = 11,400 at every node, with identical AP, DV, ids, types and parents.
- **Soma position against `all_neuron_soma_infor.xlsx`** (a mirrored-left table): 0.00 µm for 39 files after the expected reflection. The two exceptions:
  - portal 195039_009: 0.66 µm.
  - portal 201230_079: 1.92 µm.
- **Soma inside the meshes:** each soma test is exact ray parity against every triangle along 15 fixed random directions (inside = at least 8 of 15), with the exact point-to-triangle distance to the surface.
  - All 41 somas are inside the root mesh (15/15 rays each) and inside the HY mesh (15/15, except Sst 202620_016, 14/15).
  - A nucleus test against the labelled nucleus's own mesh was possible for 9 labels (PVH, ARH, VMH, DMH, LHA, ZI, AHN, PH, MPO). The labels HY, TU, MMm, PSTN, PeF, PVHd and MPN have no local mesh, so those somas were tested against HY only.
  - Every testable soma lies inside its labelled nucleus except three portal PVH neurons, which are outside the PVH mesh (0 of 15 rays) by 33.3 µm (202825_072 Avp), 25.3 µm (202825_065 Avp) and 8.9 µm (212248_071 Crh) from the mesh surface. The first report gave 35–43 µm; those were nearest-vertex distances, an upper bound. This is a small disagreement at the PVH border, of the order of a few cell diameters, not of a nucleus.
- **Node types:**
  - Each file uses either type 2 or type 0 for its non-dendritic cable, never both: 16 drawn neurons use type 2 and 21 use type 0.
  - Type 0 appears in neurons that the index calls Axon_and_dendrite and in Axon_only ones. It is reported as "type 0 (unlabelled)" and was not relabelled.
  - Undocumented types 10, 11 and 12 (type 11 was not in the access report) occur as 1–28 single nodes inside type 0 or 2 arbors (one sits inside a type 3 arbor). Every off-soma type change in the files enters or leaves one of these marker nodes; there is no change directly between types 0, 2 and 3.
  - The index says Axon_only for 195038_041 and 211266_057, and both files have no type 3 nodes. Consistent.
- **Nodes outside the root mesh:** 6 portal neurons (2,203 nodes; majority vote of three grid parity tests with rays along AP, DV and ML, instead of the first build's single AP ray; 69 nodes in all files have a split vote).
  - Five are descending axons that run past the caudal end of the CCF volume (root AP max 13,192 µm) to AP up to 14,207 µm: Avp 202825_072, Orexin 210104_007, Pmch 210951_004, and Trh 212202_026 and 212202_011. The traces continue into the spinal cord, outside the atlas.
  - Oxt 195039_009 has 10 nodes at most 5 µm outside the brain surface near the median eminence, within the 14 µm grid resolution. The first report's "41 µm" was a nearest-vertex distance; this is not a meaningful excursion.
- **Other anomalies:**
  - 5 segments are longer than 50 µm (the longest is 138 µm, in 210796_002).
  - The xlsx genotype "Wt(orexin)" differs from the portal's "Orexin" (2 neurons). This is naming only.
  - The meshes are not perfectly closed: HY has 89 open edges, AVPV 6, TM 4 and root 6. Soma tests are therefore 15-direction exact votes, and node tests a 3-axis vote on a 20 µm grid, so nodes within about 14 µm of the root surface can fall either way.
  - Line endings: all 18 Zenodo files use CRLF and all 23 portal files LF. This explains the byte difference between the two copies of a paired neuron (for example 540,214 vs 530,638 bytes for 195038_041, 9,576 lines); the parsed content, including radii, is identical.
  - Portal 201230_079 differs from the xlsx soma by 1.92 µm. The portal export may be a slightly different trace version from the paper's; when more portal files are added, compare whole arbors, not only somas, where a Zenodo copy exists.

## Actual geometry blockers for a full-scale model

1. **Access.**
   - Per-file portal access works, but the full 7,180-neuron set (about 6.2 GB) would be a bulk download. The hard rules here forbid that without an explicit decision.
   - The BSDC zips may need a login (UNRESOLVED).
2. **Licence.**
   - The portal has no licence, only all rights reserved. The BSDC record now shows CC BY-NC-ND 4.0 (`reports/swc_access.md`); `reports/conversion_map.md` line 38 has been corrected to say so. ND rules out redistributing re-framed, simplified or converted SWCs, including inside a shared viewer.
   - Only the 18 Zenodo SWCs are CC BY 4.0.
   - A shared or public model must therefore ship code and derived summaries, not portal geometry.
3. **Mirroring.**
   - The paper tables and Zenodo files are mirrored to the left; the portal and BSDC keep the true hemisphere.
   - Any mix must carry a per-file frame label. Laterality (for example the 0.96 contralateral fraction of Pomc 202457_066) exists only in the true frame.
   - `data/interim/d108_3d_brain.npz` is folded onto one hemisphere, so it can only be compared in the mirrored frame, or after unfolding.
4. **Registration error.**
   - Three PVH somas sit 9–33 µm outside the PVH mesh surface (exact distance), and every other testable soma is inside its labelled nucleus.
   - The paper's registration method and its error are behind the paywall (UNRESOLVED).
   - This is well below nucleus scale (about 100–300 µm), but matters for border assignments and for scoring at 50 µm.
5. **Truncation and atlas edge.**
   - Five of 37 neurons leave the CCF volume caudally, and their distal axons are not represented in any atlas mesh or annotation.
   - Whether the imaging stopped at a fixed level, so that spinal arbors are cut, is not documented here.
   - Complete arbors cannot be assumed for descending lines (Avp, Orexin, Pmch, Trh).
6. **Ambiguous labels.**
   - Type 0 versus type 2 differs between exports.
   - Types 10, 11 and 12 are undefined, with no official definition quoted.
   - An axon/dendrite split for type 0 files needs a documented rule, for example "every non-type-3 child of the soma". That rule has not been adopted.
7. **Sample size per line.**
   - Here there are 1–5 neurons per line.
   - In the full index, lines range from hundreds of neurons down to small counts, and Crh and Vip have no Left neurons.
   - Many lines come from few animals, so whole-animal hold-out leaves few training animals for some lines.
8. **Molecular profiles.** Neurons carry only a Cre-line label. There is no per-neuron transcriptome, so a link to HypoMap types is by line and soma position only.
9. **Radii unusable.** Calibre cannot be modelled from these files.

## What a held-out trajectory model must beat (nothing is scored here)

**Hold-out rule.** Hold out whole animals: all neurons from one brain ID (the first 6 digits of the neuron ID) are kept out together, so a neuron is never predicted from a sibling in the same brain. Report results per peptide line and for both hemispheres. Score in the true-hemisphere frame; if a baseline needs the folded `d108_3d_brain.npz`, fold both the prediction and the truth the same way.

**Baselines** (each predicts a held-out neuron's axon from its soma position and peptide line):

| baseline | how it works |
|---|---|
| B1, nearest measured neuron | Copy the axon of the training neuron (from other animals) with the nearest soma, translated to the held-out soma. |
| B2, same-line nearest | B1 restricted to the same peptide line. |
| B3, same-line, same-nucleus medoid | The training arbor with the smallest mean distance to the other same-line, same-nucleus arbors. |
| B4, straight line | A line from the soma to the centroid of each target region of the line; target regions are taken from training neurons of the same line. |
| B5, tract-weighted shortest path | Shortest path from the soma to the same targets through the 0.2 mm voxel brain in `data/interim/d108_3d_brain.npz` (its `mask`, `tract` and `fields` cost maps), with the costs fixed beforehand. |
| B6, mirrored copy | For the frame check only: the same-line nearest neuron, reflected to the held-out soma's hemisphere. |

**Metrics.** Use type 2 plus type 0 cable as "axon" only under a rule that is declared beforehand.

- **Distance:** bidirectional mean nearest-point distance (µm) between predicted and measured axon, resampled to 5 µm. Report each direction and the symmetric mean; also report the 95th percentile.
- **Region occupancy overlap:** Jaccard and Dice of the voxels occupied, at 50, 100, 200 and 400 µm.
- **Terminal-region F1:** CCF regions that contain axon tips (at least one tip, or at least 100 µm of cable), at the structure level used in the D-109 benchmark.
- **Laterality:** the error in the contralateral cable fraction.
- **Uncertainty:** report each metric with a bootstrap over held-out animals.

A model counts as better only if it beats B2 and B5 on the bidirectional distance and on terminal F1, under the whole-animal hold-out.

## Files

- Scripts:
  - `analysis/swc_fetch.py`
  - `analysis/swc_validate.py`
  - `analysis/swc_render.py`
- Raw data:
  - `data/raw/jiao2025/swc/MANIFEST.tsv`
  - `data/raw/jiao2025/swc/zenodo/*.swc` (18)
  - `data/raw/jiao2025/swc/portal/*.swc` (23)
  - `data/raw/jiao2025/swc/portal_index/mouse.neuron.info.json`
- Reports:
  - `reports/swc_validation.md`
  - `reports/swc_viewer.html` (local only; contains portal geometry)
  - `reports/swc_viewer_zenodo_only.html` (shareable build: 18 CC BY 4.0 Zenodo neurons only)
- Derived:
  - `data/interim/swc_validation.json`
  - `data/interim/swc_render.json`
  - `data/interim/swc_view_dorsal.png`
  - `data/interim/swc_view_sagittal.png`

## Review fixes (second pass, 2026-11-10)

Two reviews (geometry and viewer) checked the first build. Every verified problem was fixed; none was rejected.

| problem | fix |
|---|---|
| Scene was a mirror image (transform determinant −1) | `world.scale(-1,-1,-1)`; determinant now +1; orientation marker added; presets checked by projection and screenshots |
| "Frontal" camera was behind the brain; dorsal view had posterior at top; wrong code comment | Presets redefined (dorsal: anterior up, left on left; left lateral; frontal from anterior); comment corrected |
| 236 edges into type 10/11/12 marker nodes (548.8 µm) not drawn | Drawn in the nearest non-marker ancestor's style; disclosed in the simplification text |
| Fixed 28 µm soma spheres not disclosed | Added to the simplification text |
| three.js tag without integrity hash | sha512 SRI from the cdnjs API plus `crossorigin`; verified to load in headless Chrome |
| Viewer and conversion_map licence text outdated | Portal: copyright line only; BSDC: CC BY-NC-ND 4.0. Full viewer marked local-only; `--share` build with only the 18 CC BY 4.0 Zenodo neurons; `conversion_map.md` line 38 corrected |
| PVH distances were nearest-vertex upper bounds (35–43 µm) | Exact point-to-triangle distances: 33.3, 25.3, 8.9 µm |
| Single-axis parity inside test fragile on open meshes | Somas: exact 15-direction ray vote; nodes: 3-axis grid majority vote; votes reported |
| CRLF vs LF byte difference undocumented | Line endings counted per file and explained in `swc_validation.md` |
| 201230_079 soma 1.92 µm from the xlsx | Already an anomaly; added the advice to compare whole arbors when more portal files arrive |

Consequences of the stricter tests: nodes outside the root mesh changed from 2,159 to 2,203 (same 6 neurons, same caudal interpretation); the 195039_009 excursion shrank from "41 µm" to at most 5 µm; one soma HY vote is 14/15 (202620_016, still inside).
