# D-099 `run` output

The pre-registered test, run once by `python3 analysis/d099_gene_wiring.py run`.

```
commit 657368f6266da316f207ad71fe6cc1ecf5d4afc0; every input matches its pinned SHA-256

DECISION 1 - molecular identity predicts wiring within a region, location partialled out
  86 eligible pairs: r = +0.124   null median -0.005, 95th +0.258   one-sided p = 0.2057  -> not significant
DECISION 2 - transfer by nearest molecular neighbour beats a random same-region neighbour
  mean wiring-cosine gain +0.0466   null median +0.0008, 95th +0.0488   one-sided p = 0.0590  -> not tested (decision 1 did not pass)

SECONDARY (no decision)
  s1 genes, all pairs, location included        r = +0.105  p = 0.0589  (1240 pairs)
  s2 location only, all pairs: closeness/wiring r = +0.331  p = 0.0001
     closeness/genes                             r = +0.265  p = 0.0001
  s4 tier A/B only (relabelled among A/B)      r = +0.416  p = 0.1972  (17 eligible pairs)
  s5 all targets incl. the hypothalamus         r = +0.031  p = 0.4033  (836 structures)
  s7 within region, location NOT partialled     r = +0.213  p = 0.0801
  s6 each region left out of decision 1:
     without (Pre)Mammillary region                               r = +0.165  p = 0.1749
     without Anterior hypothalamic nucleus                        r = +0.114  p = 0.2359
     without Arcuate hypothalamic nucleus                         r = +0.207  p = 0.1019
     without Dorsomedial nucleus of the hypothalamus              r = +0.128  p = 0.2089
     without Lateral hypothalamic area                            r = +0.125  p = 0.2309
     without Lateral preoptic area                                r = +0.124  p = 0.2154
     without Medial preoptic area                                 r = +0.119  p = 0.2259
     without Paraventricular hypothalamic nucleus                 r = +0.046  p = 0.3733
     without Periventricular hypothalamic nucleus, posterior part r = +0.124  p = 0.2154
     without Suprachiasmatic nucleus                              r = +0.124  p = 0.2154
     without Ventromedial hypothalamic nucleus                    r = +0.140  p = 0.2104
     without Zona incerta                                         r = +0.105  p = 0.2534

per region (diagnostic): unadjusted correlation over eligible pairs
  (Pre)Mammillary region                                  r = -0.045  (19 pairs)
  Anterior hypothalamic nucleus                           2 pair(s)
  Arcuate hypothalamic nucleus                            r = -0.049  (17 pairs)
  Dorsomedial nucleus of the hypothalamus                 r = +0.488  (3 pairs)
  Lateral hypothalamic area                               r = +0.357  (10 pairs)
  Lateral preoptic area                                   1 pair(s)
  Medial preoptic area                                    r = -0.579  (6 pairs)
  Paraventricular hypothalamic nucleus                    r = +0.392  (17 pairs)
  Suprachiasmatic nucleus                                 1 pair(s)
  Ventromedial hypothalamic nucleus                       r = +0.859  (5 pairs)
  Zona incerta                                            r = +0.588  (5 pairs)
```
