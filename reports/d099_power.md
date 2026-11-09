# D-099 power

Produced by `python3 analysis/d099_gene_wiring.py power` before the test was run: the real gene
fingerprints, pair structure and injection sites with synthetic wiring. The last block also reads
the measured wiring against the injection sites, with no gene similarity involved. No statistic
relating genes to the measured wiring is computed.

```
populations 73; wiring structures (grey matter outside HY, ancestors removed) 663
same-region pairs 225: 86 eligible; 139 excluded (98 share a placed cell type, the rest by cross-expression)
   19  (Pre)Mammillary region
   17  Arcuate hypothalamic nucleus
   17  Paraventricular hypothalamic nucleus
   10  Lateral hypothalamic area
    6  Medial preoptic area
    5  Ventromedial hypothalamic nucleus
    5  Zona incerta
    3  Dorsomedial nucleus of the hypothalamus
    2  Anterior hypothalamic nucleus
    1  Suprachiasmatic nucleus
    1  Lateral preoptic area
rates over 400 synthetic datasets; primary = location-adjusted, one-sided 0.05; transfer = primary AND nearest-neighbour gain, each one-sided 0.05
  gene            0.00   mean adjusted r +0.006   primary 0.060   unadjusted 0.052   transfer 0.020
  gene            0.25   mean adjusted r +0.049   primary 0.107   unadjusted 0.128   transfer 0.030
  gene            0.50   mean adjusted r +0.181   primary 0.395   unadjusted 0.427   transfer 0.200
  gene            0.75   mean adjusted r +0.377   primary 0.880   unadjusted 0.902   transfer 0.650
  gene            1.00   mean adjusted r +0.517   primary 0.993   unadjusted 0.993   transfer 0.882
  position        0.50   mean adjusted r +0.015   primary 0.062   unadjusted 0.185   transfer 0.035
  position        1.00   mean adjusted r +0.069   primary 0.115   unadjusted 0.583   transfer 0.102
  cells           1.00   mean adjusted r -0.036   primary 0.025   unadjusted 0.022   transfer 0.010
  position+cells  1.00   mean adjusted r +0.031   primary 0.095   unadjusted 0.228   transfer 0.065
  leaf            1.00   mean adjusted r +0.110   primary 0.138   unadjusted 0.152   transfer 0.092
  leaf            2.00   mean adjusted r +0.152   primary 0.250   unadjusted 0.247   transfer 0.172
  leaf            3.00   mean adjusted r +0.155   primary 0.270   unadjusted 0.282   transfer 0.203
real wiring (no genes): within-region r with injection closeness +0.316; location covariates explain 14.7% of within-region wiring similarity
  position matched to the real data (strength 0.6, explains 14.0%): false pass, primary 0.085, transfer 0.045
python3 analysis/d099_gene_wiring.py power  42.52s user 4.15s system 98% cpu 47.216 total
```
