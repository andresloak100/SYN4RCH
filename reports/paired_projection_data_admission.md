# Paired molecular -> destination data admission audit (bounded, read-only)

**Verdict: REJECT.** No source was admitted from the checked public sources. Rossi 2021 is rejected, and the one
backup screened (Luo 2025) is also rejected. There is no fitting, outcome association or analysis of assay values.
The next step waits for ROOT's review.

**Corrected after ROOT's review.** This version was re-audited from the cached sources only: no new fetches, and
all 10 fetch hashes were re-verified. It changes the wording on the RNA pools, the HCR availability and Luo's
per-nucleus labels. The missing-data checklist is in `docs/paired_projection_data_request.txt`, unsent. ROOT's draft
contract is `data/interim/paired_projection_v1/test_contract_draft.json`.

**Scope:**
* Only small public metadata and text were read, each file at most 50 MB.
* Nothing was downloaded from the expression matrices, the reads, the GSA raw data, MouseLight or Jiao.
* There were no accounts and no author contact.
* Excluded on instruction:
  * Kim 2019 VMHvl, the closed pilot;
  * Li 2026 PVH, which ROOT has already marked unqualified.
* Code: `analysis/paired_projection_data_admission.py` (SHA-256 `acb52e10…`).
  * The `cached` mode re-audits from `raw/` after checking each fetch hash.
  * Only per-cell HCR file names count as found HCR data: HALO `*_quant.xlsx` exports or `HCR_data.csv`. RNA DEG or
    Seurat CSVs never count.
  * Finding such a file can only move a source to PENDING. It never verifies cell pairs.
* Tests: `tests/test_paired_projection_data_admission.py` (`93626f48…`), 9 passing. They include a fixture where an RNA
  DEG CSV is present and the verdict is still REJECT.
* Output: `data/interim/paired_projection_data_admission/admission.json` (`0f0316be…`), with the raw copies in
  `raw/`. Earlier versions were `f32aa182…` and `e10d1522…`.
* ROOT's admission-only contract draft is `docs/paired_projection_test_contract.json` (`1dec9a86…`). It authorizes
  no evaluation.
* The data request stays unsent until the user authorizes it.

## 1. Rossi et al. 2021, LHA Vglut2 -> LHb / VTA: REJECT

**Source:** Neuron 109:3823 (2021), DOI 10.1016/j.neuron.2021.09.020, PMC8812999.
* **Article licence:** a PMC author manuscript, available for text mining under fair use. It is not CC.
* **GEO:** no data licence stated.
* **GitHub repo:** no licence file.

**What is public:**

| item | content |
|---|---|
| GEO GSE169176 | 3 samples, GSM5179460–62, "LHA sample pool 1–3", each "pooled ... from LHA tissue of two animals (male and female)". Only filtered 10x h5 matrices (56.7, 66.6 and 68.3 MB; RAW.tar 191.7 MB). Not downloaded: over the cap, and they hold assay values. |
| SRA SRP311270 / BioSamples | 3 runs, one per pool. No animal, sex or subject fields. |
| Pooling, as each source puts it | **GEO:** three samples, each a two-animal pool (male and female). **Paper:** "tissue from all animals was ... pooled" before three reactions, which is ambiguous between all six together and within each pool. Both statements are kept. |
| Paper, scRNA methods | 7 mice (4 male, 3 female). One was used for histology. "Tissue from all animals was ... pooled" before 3 Chromium reactions. Labelled neurons: eYFP (LHb) 257, tdTomato (VTA) 1,595. |
| Paper, HCR methods | 3 mice (1 male, 2 female), separate from the RNA mice. 9 genes in 3 rounds: 7 markers (Pax6, Sostdc1, Vgat, Vglut2, Nptx2, Pdyn, Pitx2) plus eYFP and tdTomato. 6–12 sections per mouse, HALO puncta per cell, thresholds stated. |
| Data availability | Only GEO, plus "All code ... on the Stuber lab Github". No HCR deposit. |
| GitHub `stuberlab/Rossi-et-al-2021` @ `41e5b21f` | 10 code files only (notebooks and R scripts); no data files (`.xlsx`, `.csv`, `.h5` and so on). |
| HCR notebook (`Fig3_S2_HCR analysis_Final.ipynb`) | Reads Dropbox HALO exports `AG09_C1/C2`, `AG10_A2`, `AG11_C1/C2_quant.xlsx`, i.e. mice AG09, AG10, AG11. Derives mouse, slide, slice and hemisphere from the sheet names. Hand-assigns AP coordinates per slice. XY comes from cell bounding boxes. Printed: 92,866 cells, 816 YFP+ and 857 tdTomato+. Its SVM uses random 50% cell splits, not splits by animal. |
| Supplement PDF | PMC returned an empty body to a script. Europe PMC says the article is not open access. Not read. |

**Requested items:**

| item | public? |
|---|---|
| Per-cell HCR counts | **Not found** in the checked sources |
| Exact cell -> animal / section join | **Not found**: the sheet-name scheme is visible in the code, but the tables were not found |
| Measured coordinates | **Not found**: bounding boxes and AP assignments are only in the HCR tables |
| Tracer labels per cell | RNA: in the h5 matrices (not opened). HCR: **not found** |
| Per-cell animal ID in the RNA data | **Not found**: there are 3 two-animal pools, and the IDs are not recoverable from the checked sources |
| Assays | stated |
| Animal counts | 6 RNA mice in 3 pools + 3 HCR mice |

Checked sources: GEO, SRA, BioSample, the PMC text and the authors' GitHub. Not retrieved: the PMC supplement PDF
and the publisher's supplementary files. So "not found" means not found in the checked sources. It is not a claim
that the data are non-public everywhere.

**Hard blockers:**
1. The HCR per-cell data (tables, joins, coordinates and tracer calls) were not found in the checked public sources.
   Author contact is not allowed here. So no external RNA -> HCR animals are available for testing now.
2. The RNA data come from three two-animal pools in three reactions.
   * The pools may count as biologically independent units, but only if animal-disjointness is confirmed.
   * They never give six recovered animal IDs.
3. Tracer identity is confounded with destination in both assays:
   * eYFP (hSyn) always marks LHb, and tdTomato (CAG) always marks VTA;
   * only the separate ephys cohort was counterbalanced.
4. Even if the HCR tables were released:
   * 3 mice and 7 shared genes make it feasibility or development only, not a powered confirmatory test;
   * it would also need a locked RNA -> HCR transform, calibrated outside the HCR test animals.
5. There is no novelty: the authors already classified eYFP vs tdTomato with an SVM in both assays.

## 2. Backup screened: Luo et al. 2025, SuM -> dorsal / ventral DG: REJECT

**Source:** Neuron (2025), DOI 10.1016/j.neuron.2025.09.006, PMC12581905, text-mining licence. Data are in GSA
CRA019056; no raw data were downloaded. This confirms ROOT's check.
* **Patch-seq:** the Methods say "Two cells were collected in one tube as a sample." The 22 profiles therefore cannot
  be taken as individually paired cells.
* **snRNA-seq, labels:** eGFP was added to the reference, and nuclei with an eGFP count above 0 are called positive.
  * So label-positive nuclei can pair RNA with a captured label.
  * The checked text reports only four positives, all from vDG mice, and states no dDG count.
  * The remaining nuclei are unknown (not captured), not negatives.
* **snRNA-seq, animals:** the Methods say "Tissues from 5 mice of the same group were pooled together." The pool IDs,
  the number of pools (replication) and the sample unit are unresolved, so independent evaluation animals are not
  established.
* **Endpoint:** dDG vs vDG is a different endpoint from the Epi-retro-seq list. It is not a hard blocker for a new
  task.

**Not screened:** this audit was bounded to one backup, and Li 2026 PVH was not re-checked. In the earlier search
notes, Calvigioni 2023 and Chung 2017 have a single target with all cells positive, so they cannot give a
destination contrast.

## Source receipts (`raw/`, fetched 2026-11-12)

| file | URL | bytes | SHA-256 |
|---|---|---:|---|
| rossi_geo_soft | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169176&targ=all&form=text&view=brief | 5,800,338 | `a18083db7c1ce82cec10c84416a52e07e5aa816d1ec2a65d38d7602b9e942081` |
| rossi_geo_filelist | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169176/suppl/filelist.txt | 346 | `572ca19777ecaa0a84ece6ee06c0f0cb1f19ee18de848e622869c9080c899c8d` |
| rossi_sra_runinfo | https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=SRP311270 | 1,936 | `acf315cbf08fc3cd869e29e9cb6cc043c06132520e1fb973e657e6f033f7f3e1` |
| rossi_biosample | eutils efetch biosample SAMN18350787,85,84 | 4,771 | `23bb0e98ce266e8af84513b8ac6a00a66a4f2371564ce08ffeabe721bcc39f11` |
| rossi_pmc_xml | eutils efetch pmc 8812999 | 194,235 | `9515f8cc1a2b5b8481229473bb52681f535cac9328a52b6c676f93089385acd3` |
| rossi_epmc_supp | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8812999/supplementaryFiles | 296 | `2c4a397e867b0f2fbd525c70f8fed35f4b51df7a841773f723c43b7471a10167` |
| rossi_github_tree | api.github.com repos/stuberlab/Rossi-et-al-2021 tree @41e5b21f | 3,205 | `3766c3314fc4e7dbef0cbdb320f7ca60dc6e31115b382ef2fef467aed2c3aa4c` |
| rossi_github_repo | api.github.com repos/stuberlab/Rossi-et-al-2021 | 6,433 | `3dcb8a907198334c93bcdda84842d73a66e2fe855c3d5616b500a7165e63b8bc` |
| rossi_hcr_notebook | raw.githubusercontent.com …/41e5b21f…/Fig3_S2_HCR analysis_Final.ipynb | 526,068 | `554a2357bd59b1a983cb29b32dc8042eec45c3030373411d608b911119fbc4e6` |
| luo_pmc_xml | eutils efetch pmc 12581905 | 211,821 | `a0a7a9d6ebf73d811c0fc68ad8314af0ce77fccfcf01aebde36807021f27a288` |

The GEO and GitHub API responses can change over time, so a re-fetch may give different hashes. The values above
are from this run.

## ROOT verification

**Result: PASS.** ROOT independently checked the following and found them all correct:
* all 10 cached source hashes and the code pin;
* its own parse of GEO, which gives three two-animal pools;
* that the HCR tables are missing from the pinned repo;
* Luo's actual methods, including the eGFP-positive statements;
* that the contract copy is unchanged;
* that the old test and the v1/v2 result hashes are unchanged.

ROOT ran the 9 tests itself, and all passed.

**Receipt:** `data/interim/paired_projection_data_admission_root/check.json`, SHA-256
`82a55d6554ce9cf609e868f5fb7f1067518ef5cfc651b37c4ba92d5def0612e9`.

**Contract:** it stays ADMISSION_ONLY. No model run is justified.

**Public request draft:** ROOT's draft GitHub issue is `docs/paired_projection_public_data_request.json`, SHA-256
`009e5be6a2f717531c68e524da489132863f672f760841b251d681f433cf211b`.
* Status: UNSENT_REQUIRES_USER_AUTHORIZATION.
* It has not been posted, and no one has been contacted.
* If the user approves a post, ROOT will make it separately.
