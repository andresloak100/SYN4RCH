# Published algorithms for neuronal projection prediction

Synarch's next reference model should be regularized logistic regression on genuinely measured genes and position, with a strong nonlinear position control. Naive Bayes and XGBoost provide published comparators. Shared latent factors are a later option when the measurements support them. These are proposed adaptations; no model below has been run or validated by this literature review.

The current public SPIDER spatial object supplies a measured 32-gene panel and captured projection-barcode fields. Animal identities and target-by-animal assay masks remain unresolved. Literature examples do not remove that requirement, and destination prediction does not establish complete axon paths or synaptic connections.

## Direct prediction references

| Reference | Published algorithm | Relevant use | Validation boundary |
|---|---|---|---|
| [Wang et al., eLife 2023](https://doi.org/10.7554/eLife.84262), central amygdala | L2 logistic regression from 29 measured marker genes; recursive feature elimination | An interpretable small-panel reference, close in feature count to our spatial panel | Main scores use repeated stratified cell folds. Selected marker predictions were also checked in a third animal; this is narrower than validating the complete classifier across many animals. |
| [Xu et al., eLife 2024](https://doi.org/10.7554/eLife.85419), MERGE-seq | Target-specific Naive Bayes from endogenous RNA | A simple molecular comparator for captured target labels | Random cell splitting and repeated subsampling do not establish animal-held-out generalisation. Feature selection must be confined to fitting groups in Synarch. |
| [Zheng et al., National Science Review 2026](https://doi.org/10.1093/nsr/nwag004), SPIDER-Seq | XGBoost binary classifiers from transcriptomic PCs and XYZ | A nonlinear comparator; relevant to PFC inputs to LHA | The inspected workflow uses cell splitting. Comparisons with MERFISH spatial patterns are not same-cell projection ground truth. Its full-transcriptome inputs differ from our measured 32-gene panel. |

The amygdala example demonstrates that a small measured gene panel can carry projection-label information in a defined circuit. It does not establish hypothalamic transfer, exclusivity of a neuron's destinations, or a causal genetic wiring rule. Its unlabelled class is an assay status, not proof of no projection. [Wang et al.](https://doi.org/10.7554/eLife.84262)

## Spatial controls and shared representations

**BARseq2:** [Sun et al., Nature Neuroscience 2021](https://doi.org/10.1038/s41593-021-00842-4) combines NMF projection modules, Gaussian-process regression from laminar depth, and gene-related regression. The most useful immediate lesson is a strong nonlinear anatomy control: a weak linear position model could overstate a molecular gain. The inspected projection-variance analysis uses neuron-level cross-validation. Any output-module basis must be fitted inside the training groups, and a target-derived module cannot become an unknown neuron's input.

**Sparse reduced-rank regression:** [Sorensen et al., Nature 2026](https://doi.org/10.1038/s41586-026-10424-8) learns gene-linked local phenotypic factors and uses factors plus source location in projection models. RNA/local phenotypes and whole-neuron morphology come from separate cohorts, so these are not individual RNA/full-axon pairs. The checked projection script withholds cells, not animals. This is a reference for shared representations, not a demonstrated hypothalamic path model. With only our measured small panel, a direct predictor is the simpler first comparison.

## Unlabelled cells and detection

[Elkan and Noto, KDD 2008](https://doi.org/10.1145/1401890.1401920) describes positive-unlabelled learning. Its probability correction assumes labelled positives are randomly selected from the positive population. Viral capture may depend on cell properties, tracer and injection site, so Synarch cannot assume that condition holds. This is an optional framework for modelling detection, not an automatic way to turn barcode nondetection into anatomical negatives. The [authors' paper](https://cseweb.ucsd.edu/~elkan/posonly.pdf) explains the assumption.

In particular, the current SPIDER zero-capture records must not silently become verified non-projecting cells. Until the assay mask and detection meaning are established, distinguish a captured-label probability from a biological projection probability.

## Proposed Synarch reference pipeline

This is a design reference, not a frozen experiment or permission to fit the current data.

1. Establish exact cell identities, disjoint animal or justified pool identities, assay masks, detection rules and independent evaluation power.
2. Separate fitting, calibration and evaluation animals before feature selection or scaling. An animal's slices and cells stay together.
3. Use measured endogenous genes and actual source coordinates. Exclude reporter/vector genes, target-derived classes and hidden outcomes. Keep inferred coordinates or transcriptomes explicitly separate.
4. Fit one regularized logistic candidate. Compare frequency/constant, technical, anatomy-only and nearest-neighbour controls; retain Naive Bayes and XGBoost as references for a prospectively specified comparison.
5. Use a nonlinear location control appropriate to the sample size, such as location nearest-neighbours. Exact Gaussian processes are not a practical default for the full spatial object.
6. Assess animal-averaged Brier score and log loss, plus separate reliability/calibration diagnostics and discrimination. Choose margins and the familywise rule before evaluation. Calibrate outside the evaluation animals.
7. Advance to an independent geometry experiment only if destination inference provides a reproducible gain. No new 3D prediction or discovery follows from this reference list.

Language-independent outline:

```text
admit exact measured pairs, biological groups and target-assay masks
reserve evaluation groups
inside fitting groups: select/scale permitted inputs and choose regularization
fit probability(target capture | measured genes, measured position)
seal predictions before revealing evaluation outcomes
compare with strong controls on the same groups and proper scores
report the complete gate, calibration, failures and assay limitations
```

[Friedman, Hastie and Tibshirani, JSS 2010](https://www.jstatsoft.org/article/view/v033i01) provides the penalized-regression algorithm and coordinate-descent implementation reference. The [official glmnet vignette](https://cran.r-project.org/web/packages/glmnet/vignettes/glmnet.pdf) supports explicit fold assignments. [GroupKFold](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.GroupKFold.html) documents non-overlapping group folds. Neither supplies missing animal identities.

[Probability calibration documentation](https://scikit-learn.org/stable/modules/calibration.html) also distinguishes proper scores from calibration alone: lower Brier loss does not by itself prove better calibration. Avoid reporting accuracy alone for an imbalanced captured-label endpoint.

## Code and data references

| Method | Reference resource | Reuse status checked |
|---|---|---|
| Regularized logistic regression | [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html); [scikit-learn](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html) | glmnet lists GPL-2. Use a version-pinned implementation rather than copying a paper's environment. |
| MERGE-seq | [Author analysis repository](https://github.com/MichaelPeibo/MERGE-seq-analysis/tree/revised) | No repository licence established in the check. Reference the method; do not assume unrestricted source redistribution. |
| BARseq2 | [Processed data and analysis deposit](https://data.mendeley.com/datasets/jnx89bmv4s/2) | Deposit lists CC BY 4.0. |
| SPIDER-Seq | [Pinned prediction script](https://github.com/ZhengTiger/SPIDER-Seq/blob/047eb5aaa15087c241570745bd6df07882b0dd78/Main%20Figures/Figure7.R); [public spatial resource](https://huggingface.co/spaces/TigerZheng/SPIDER-web/tree/main/data) | Author analysis code is MIT. Check processed-data terms separately. |
| Sorensen et al. | [Pinned projection script](https://github.com/AllenInstitute/exc_vis_manuscript/blob/7f187ac1efdc0d81b18b5d539fa602b3bae1a550/scripts/glm_lf_subclass_loo.R); [licence](https://github.com/AllenInstitute/exc_vis_manuscript/blob/main/LICENSE) | Repository licence contains a noncommercial restriction; it is not ordinary BSD/MIT. |
| Central amygdala | [EASI-FISH data deposit](https://doi.org/10.25378/janelia.21171373) | Data reference established; a standalone prediction-code licence was not established in this review. |

Checked 12 November 2026. Only original method summaries and pseudocode are included here. No third-party algorithm source, model weights or research datasets are republished, and no author code was executed.
