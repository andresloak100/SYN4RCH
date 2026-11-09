# D-101 `run` output

The pre-registered test, run once by `python3 analysis/d101_programmes.py run`.

```
commit 51aa9e657f405a023744996e2b6322c52499f98f; every input matches its pinned SHA-256

FAMILY: D-099 all genes (p = 0.2057, already run) and the two programmes; Holm-adjusted one-sided p < 0.05 passes
  all_genes              (D-099)                                        Holm 0.6171
  guidance_adhesion      r = -0.057  null 95th +0.239  p = 0.6522  Holm 1.0000  -> not significant
                         transfer gain +0.0207  null 95th +0.0488  p = 0.2496  -> not tested
  transcription_factors  r = -0.090  null 95th +0.193  p = 0.7858  Holm 1.0000  -> not significant
                         transfer gain -0.0087  null 95th +0.0466  p = 0.6258  -> not tested

SECONDARY (no decision)
  guidance_adhesion      within region, location NOT partialled  r = +0.039  p = 0.3888
                         correlation with the all-gene similarity over the eligible pairs +0.869
  transcription_factors  within region, location NOT partialled  r = +0.033  p = 0.3738
                         correlation with the all-gene similarity over the eligible pairs +0.818
```
