# Mouse hypothalamus: what is measured together, and the conversions it allows

The reconstruction (D-098) is built from **conversions**. Each one learns to turn a
property we have into a property we need, from a dataset where **both are measured in the
same cells**, and is graded on cells it never saw against a location-only baseline. A chain
of conversions carries its uncertainty forward. Nothing that is measured nowhere is filled in.

Compiled 2026-11-09 from five independent searches and a final pass that merged them and
re-opened the records for Epi-retro-seq, Rossi, Kim, Li 2026, Jiao and Timonidis (workflow
`matched-data-search-hypothalamus`). Every other record was opened by the search that found it.
Counts come from deposited metadata where it was tabulated, otherwise from the paper's text.
A licence is given only where one was read; "none stated" means the record names none.

## What exists

`G` genes (transcriptome), `M` DNA methylome, `T` genetic tag (Cre or promoter line), `L`
location, `P` projection target(s), `A` complete axon (single-neuron morphology), `E`
electrophysiology.

### Same-neuron: identity and projection measured in the same cells

| dataset | where | measured together | cells | access |
|---|---|---|---|---|
| Kim et al. 2019, Cell, doi:10.1016/j.cell.2019.09.020 | VMHvl | G + P: one of 4 targets (dPAG 491, lPAG 395, MPOA 317, AHNa 138) | 1,341 of 4,574 sequenced; 16 male mice | Mendeley Data 10.17632/ypx3sw2f7c, CC BY 4.0; metadata 0.7 MB, matrix 188 MB |
| Li et al. 2026, Cell Reports, doi:10.1016/j.celrep.2025.116904 | PVH | G + P: spinal cord or parabrachial (NTS-projecting libraries listed); MERFISH in separate mice places the clusters | per target only in Table S12, not verified; the paper gives percentages per cluster | Zenodo 10.5281/zenodo.15983704, CC BY 4.0, 13.7 GB in all (projection objects 43 MB and 200 MB; includes HypoMap PVH cells integrated with theirs); GEO GSE303256 |
| Zhou et al. 2023, Nature, doi:10.1038/s41586-023-06823-w (Epi-retro-seq) | whole HY as one pooled dissection; includes some thalamic-type cells | M + P: one of 10 targets, from 17 experiments (AMY male only; PAL and VTA female only); Allen subclass per cell (47 subclasses) | 1,538 passing QC (the paper's figure uses 1,572) | GEO GSE230782, none stated; NeMO copy CC BY 4.0; labels in Supp. Table 2 (7 MB); methylome 60 GB as one archive, about 200 GB as per-cell files |
| Rossi et al. 2021, Neuron, doi:10.1016/j.neuron.2021.09.020 | LHA glutamate | G + P: LHb 257 or VTA 1,595, read from reporter transcripts; no label table is deposited and the calling threshold is not stated | 1,852 | GEO GSE169176, 183 MB |
| Lewis et al. 2020, Neuron, doi:10.1016/j.neuron.2020.10.002 | PVN oxytocin | G + P: takes up blood-borne Fluoro-Gold (projects outside the blood-brain barrier) 91, or not 55 | 146 (the GEO text says 215 magnocellular, the deposited table 125: unresolved) | GEO GSE147092, small |
| Winter et al. 2023, Nature, doi:10.1038/s41586-023-06817-8 | HY among 8 dissected regions | G + P: spinal cord, cervical or lumbar | HY count not verified (2.9% of spinal-projecting neurons by imaging) | GEO GSE247602; NeMO; some of its cells are in Li 2026's spinal object |
| Beine et al. 2022, J Neurosci, doi:10.1523/JNEUROSCI.1197-22.2022 | PVH and LHA clusters among 14 populations | G + P: lumbar cord | not verified | GEO GSE212409; some of its cells are in Li 2026's spinal object |
| Chung et al. 2017, Nature, doi:10.1038/nature22350 | POA | G + P: all to TMN, so no contrast between targets | 76 | GEO GSE79108, small |
| Calvigioni et al. 2023, Nat Neurosci, doi:10.1038/s41593-023-01367-8 | LHA | G + E + dendritic morphology; P: all to LHb, so no contrast between targets | 163 | BioStudies S-BSST1069 |

### Single-neuron projectomes: complete axons with a genetic label

| dataset | measured together | neurons | access |
|---|---|---|---|
| Jiao et al. 2025, Nat Neurosci, doi:10.1038/s41593-025-01919-0 | T (15 Cre lines plus an orexin promoter virus) + soma L in CCFv3 + A | 7,180 from 166 male mice | digital-brain.cn portal (single SWC files without login); BSDC doi:10.12412/BSDC.1667278800.20003, 6.2 GB; cable-length matrix 13 MB. **No licence stated**: the BSDC licence field is empty and asks only for attribution (one search's "CC BY-NC-ND 4.0" was wrong); the NeuroXiv mirror is CC BY-NC 4.0 |
| Li et al. 2024, Neuron, doi:10.1016/j.neuron.2023.12.022 | T (Oxt) + L + A; median eminence only (177) or widespread (87) | 264 from 3 brains; 106 are also in Jiao | same portal; none stated |
| Janelia MouseLight | L + A (9 with a Cre label) | 77 in classic HY (97 with any HY soma) | NeuronBrowser; 2024 figshare items CC BY 4.0, 2017-2019 items CC BY-NC 4.0 |

### References that place cell types and carry labels between atlases

| resource | gives | access |
|---|---|---|
| HypoMap (Steuernagel et al. 2022) | G for 465 HY cell types, 384,925 cells | data/raw/hypomap (in use) |
| Allen Brain Cell Atlas (Yao et al. 2023) | G and the whole-brain taxonomy (5,322 clusters); HY 10x files 7.3 and 2.9 GB | CC BY-NC 4.0 |
| Zhuang MERFISH (Zhang et al. 2023) | 1,122 genes; about 2.6 million cells with CCF coordinates | CC BY 4.0 |
| MapMyCells | maps any cells onto the Allen taxonomy; no HypoMap-to-Allen crosswalk is published | brain-map.org |
| Allen Mouse Brain Connectivity Atlas | P per Cre-line population (D-099's data) | public API |
| Allen Mouse Brain Atlas ISH | bulk G per brain structure: 3,318 coronal datasets, 3,190 genes | public API |

**Nothing at scale measures a full transcriptome and a complete axon in the same
hypothalamic neuron, and in every same-neuron dataset above each cell has exactly one
target.** No hypothalamic neuron has a measured vector of several targets alongside its
genes. The chain has to bridge that gap.

Looked at and not usable as training data: Berkhout et al.'s integrated PVN atlas (its
projection labels were transferred by similarity, so testing genes against them would be
circular); MERGE-seq (the hypothalamus is only a target); Connect-seq (input side, not output).

## Conversion paths

| conversion | trained on | graded how | ready |
|---|---|---|---|
| identity -> projection target across HY | Epi-retro-seq (10 targets, Allen subclass per cell), with subclass -> RNA fingerprint and subclass -> 3D position from the Allen atlas | D-099's statistic over subclass pairs, location partialled | labels now (7 MB); Allen 10x hypothalamus files about 10 GB; methylome (60 GB) not needed |
| identity -> projection target, location fixed by design | Kim 2019 VMHvl (4 targets); Rossi 2021 LHA (2); Li 2026 PVH (2-3); Lewis 2020 PVN (2) | held-out animals | now; under 1 GB |
| tag + exact position -> single-neuron projection | Jiao 2025 | held-out neurons and held-out lines, against position alone | now (13 MB matrix) |
| target area's genes -> how strongly it is hit (relational) | Allen ISH + connectivity, Timonidis 2020's design | spatially blocked folds and a location baseline, which the original lacked | now; about 4,700 API calls |
| genes -> 3D position of every cell type | HypoMap mapped onto the MERFISH cells (MapMyCells or Tangram) | held-out genes and cells | now; several GB, not all sizes checked |
| HypoMap type -> tag -> projection | marker genes linking HypoMap types to Jiao's lines | held-out lines | after the two above |

## Order

1. **D-099's question again, with identity and target measured in the same nucleus:
   Epi-retro-seq.** Removes the Cre-line matching that limited D-100, and covers many nuclei
   and 10 targets. Each cell already carries an Allen subclass name, so two conversions attach
   without a crosswalk: subclass to RNA fingerprint (Allen 10x hypothalamus data) and subclass
   to 3D position (Allen MERFISH cells in the hypothalamus). Unit: Allen subclasses with at
   least 20 neurons, less the thalamic-type RE-Xi Nox4 Glut (21 subclasses, 1,182 neurons,
   210 pairs, by the final pass's count). Projection profile: each subclass's enrichment in
   each target's experiments relative to its share of all HY neurons. Then D-099's statistic,
   location partialled. Limits: methylation, not RNA, decides each cell's subclass; position
   is known per subclass, not per cell; three targets come from one sex only.
2. **Identity to target inside one nucleus: Kim 2019 VMHvl.** Per-cell RNA, four targets,
   under 200 MB; a learned classifier graded on cells it never saw. The paper already reports
   that most VMHvl types contain cells labelled from more than one target, so this is a
   replication with a stricter test. 16 male mice; 3 gave cells labelled from both MPOA and
   dPAG, the rest one target each, so target and animal are largely entangled and grading
   holds out whole animals.
3. **Tag and position to single-neuron projection: Jiao 2025.** Tests, neuron by neuron,
   whether the genetic tag adds to exact soma position, and supplies the axon layer of the 3D
   map. No licence is stated: analysis only, nothing derived from it redistributed, until the
   authors confirm terms.
4. **Genes to 3D position** for all 465 HypoMap types: the spatial base of the map.

## Corrections this search made

* **Timonidis et al. 2020 predicted from the target area's genes, not the injection site's.**
  For each injection, the projection strength into each of 428 target structures was predicted
  from that target's bulk ISH profile (3,318 datasets), with folds over target structures, not
  spatially blocked, and no location baseline. Ridge median R-squared 0.54 over 1,397
  experiments. The median AUROC of 0.89 is, by one search's reading of the code, the maximum
  over 61 cutoffs chosen on the test predictions, which would make it optimistic; the final
  pass did not re-check this. D-099 described the features as injection-site
  expression; that was wrong. Its conclusion still holds for a different reason: a target's
  expression profile also encodes where the target is, and nothing in the design separates
  the two. The result supports the relational idea (target genes), not a cell's own genes
  predicting its wiring.
* **HypoMap under-represents PVH subtypes.** Li et al. 2026 integrated HypoMap's 5,119 Sim1+
  "PVH" neurons with their own: seven Sim1+ clusters were almost entirely their cells, and
  many HypoMap PVH cells express peri-PVH markers (Cabp7, Onecut3, Gsc). D-099 placed
  Cre-line tracing on HypoMap types, so part of its null may be HypoMap's resolution in PVH.
