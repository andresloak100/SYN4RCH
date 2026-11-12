# Paired-data search and frozen exploratory transfer test

2026-11-12. BARseq2 A1 is admitted for one descriptive, two-mouse test of captured projection destinations. No biological score had been computed when the admission and protocol were frozen.

## Admitted data

The public [Sun et al. BARseq2 deposit, version 2](https://data.mendeley.com/datasets/jnx89bmv4s/2) is CC BY 4.0. Its compact A1 object pairs measured marker counts, image-derived soma depth and observed target barcode counts within cells. [The paper](https://doi.org/10.1038/s41593-021-00842-4) describes this paired assay; original processing scripts and sample sheets establish the animal crosswalk.

The frozen cohort contains **598 cells: 304 from XC95 and 294 from XC97**, with **23 measured genes** and **11 anatomical targets assayed in both animals**. The separate 593-row excitatory-only animal-index file is not used to identify the larger matrix. Composite animal/local-cell keys prevent collisions between mice. The pooled source matches both original experiment objects after the documented target-column reorder.

QC retains cells with a raw projection count above 2 somewhere in the source panel and at least 5 Slc17a7 or Gad1 spots. Each target's outcome is raw target count above 2. Lower counts mean assay nondetection, not proven biological absence. The olfactory-bulb control is excluded as a target. Five inhibitory cells remain in this published paired cohort.

Anatomy is **one measured cortical-depth dimension**, normalized using the author's top-distance / cortical-thickness convention. It is not absolute XYZ and does not control tangential position. Technical controls use measured panel-total and unused-code counts; panel total also contains biology.

ROOT extracted the four input tables without fitting. An independent comparison checked every retained gene value (13,754), every label (6,578), depth, technical covariates and all 22 mouse/target mask entries against the original MAT data. All 19 compact source-download hashes match the provider's hashes. Source and extraction pins are recorded in [the admission](../docs/paired_projection_barseq2_admission.json).

## Frozen comparison

Train on one animal and evaluate the other, then reverse direction. All slices and cells of an animal stay together. The sole candidate is the fixed gene-plus-depth L2 logistic model. Controls include frequency, depth-only logistic and 15-neighbour averaging, molecular-only averaging, copying one depth-nearest donor, and copying or averaging donors under the candidate's exact gene/depth representation. Every donor is from the training mouse. All methods score identical cells and targets.

The [protocol](../docs/paired_projection_real_protocol.json), code and admission are fixed before the single run `barseq2_a1_598_logo_v1`. ROOT reproduced the synthetic software results exactly and **121 focused reference/runner tests pass**. Runtime checks bind the animal mapping, inputs, outcomes, folds, configuration, runner, imported reference and protocol. An exclusive started record prevents silent reruns, and keyed out-of-fold predictions permit independent score replay.

Results will report animal-equal Brier score and log loss, each animal and target, and the target-equal aggregate. Two animals do not support a population-confirmatory win or a significance claim. Several targets have only one to five captured cells in a mouse. The cohort is conditional on capture somewhere and cannot establish all-neuron probabilities. This is auditory cortex, so any gain does not establish hypothalamic transfer, complete axon reconstruction or synaptic circuits. No tuning or recalibration follows the result.

## Other sources checked

| Source | Qualification for this test |
|---|---|
| Wang et al. CeA Projection-FISH | Public exact-cell table, two mice and five shared targets. Shared XYZ origin, units and inter-animal transforms were not established; therefore no strong absolute-position comparison was admitted. |
| Insula Retro-MERFISH 2026, GSE329718 | Seven tracing mice, but each target was assayed in only one mouse. Sections cannot provide independent same-target animals. The large image/object archive was not downloaded. |
| Axonal BARseq2 VISp 2026 | Promising paired measurements from two mice. Synapse lists a scripts archive but anonymous download is denied; no compact paired table or complete assay mask was recovered. |
| MOp Retro-MERFISH and Fan thalamic BARseq | Checked public access did not supply an executable paired table with the required metadata. |

Detailed metadata and provenance remain in `data/interim/paired_projection_search_round4/`. Only the previously approved Rossi issue was posted; it had no reply at the check in this work period. Closed VMHvl, Jiao and completion evaluations remain closed, and no viewer output is promoted by this protocol.
