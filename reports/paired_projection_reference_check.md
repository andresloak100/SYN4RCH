# Paired projection reference benchmark: software validation

2026-11-12. **Software checks pass; biological prediction is NOT_EVALUATED.**

Claude implemented one fixed L2 logistic model using measured endogenous markers and measured position, with frequency, technical, location-neighbour, molecular-neighbour and Naive Bayes controls. The [protocol](../docs/paired_projection_reference_protocol.json) keeps captured-target prediction separate from the exclusive two-target endpoint. Neither endpoint establishes an axon's complete geometry or a synaptic connection.

ROOT reviewed the implementation with an independent reviewer and reproduced its synthetic run. All **68 focused tests pass**. The planted-data run passes all software checks and 21 deliberate refusal cases. ROOT's independently reproduced macro losses match Claude's exactly; hand-calculated animal-equal losses also match. These numbers describe simulated data only.

Review fixes enforce non-synthetic refusal before data access, assay-mask and endpoint-group checks before expression parsing, exact probability-vector shape and range, predictor-name/role checks, and nonblank animal identifiers. Canonical cell-ID ordering removes a small row-order dependency in the molecular-neighbour control. Shuffling all four input tables now leaves every method's predictions and metrics identical. Fitted preprocessing remains unchanged when hidden evaluation values are altered.

Frozen review hashes:

- Source: `9ae3dd30a7fad55c6b7d232e8faad4ec678bf406c105faa61cade544039723a3`
- Protocol: `c9df12ffffcb696be2b94b57d99a873fb0e1f81d3dc0f8e9c95cb95d57aa7336`
- ROOT receipt: `data/interim/paired_projection_reference/root_review/check_final.json`, SHA-256 `dc58afc48ddd4e696eb80b1801a378b145c0088a84ed73d048e976ff79cadb9c`

## Biological data admission

No real-data fit was run. The [SPIDER primary paper](https://doi.org/10.1093/nsr/nwag004) describes its main 124,829-cell spatial projectome across 36 PFC slices in one mouse. The public object's matching cell/slice counts support its association with that cohort, but no authoritative cell-to-animal crosswalk was recovered. The paper also reports a repeat in another mouse; that observation does not establish a released multi-animal paired table or an unseen-animal prediction test. Slices and hemispheres must not stand in for animals.

The separate three-mouse sequencing experiment covers different injected targets; its sample count cannot be assigned to the spatial object. Per-animal target coverage, assay/QC masks, repeat-animal paired measurements and enough independent animals remain necessary. LHA is a destination of these PFC cells, so this resource concerns cortical input to hypothalamus rather than hypothalamic-source circuit reconstruction.

The checked supplemental archive transfer did not complete, and the separate spatial-data listing returned an access error. Neither failure establishes that the missing metadata does not exist. The [approved Rossi data request](https://github.com/stuberlab/Rossi-et-al-2021/issues/1) is posted once; no reply was present at the last check in this work period.

The source/admission record is in `data/interim/paired_projection_search_round3/spider_independence_check.json` (SHA-256 `64be93efe89a0e0af8ef0259284f6c91260ba25b9a29803ba188950e43bda909`) and the committed search summary. No improvement on biological prediction, new circuit discovery or viewer promotion is claimed. Existing closed biological tests remain closed.
