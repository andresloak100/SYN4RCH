# D-095 `--power` output

Produced by `python3 analysis/d095_pair_test.py --power` before the test was run. It
never computes the observed statistic.

```
inputs (SHA-256):
  genes_npz            5b103fecaacc8d814966db2e2f638e120c7de2d907f29e4755ee44af5155a496 (1 file)
  hy_injections        031aa47061ec084737e5ef0128e3027ccf1d9f1dba86864be308352f55ae9fc0 (1 file)
  circuit_experiments  b2ec3ef17fcd589a9b34c3028730ab778a63b8ee16a06161db669b6c2d3e8f8b (1 file)
  projections          aca145fb21d6867db9c526b0a00840c7adbf589f3dab6115e33efab274768692 (43 files)
  structure_graph      b5f0d024d1df09ee18aef15093f55240f4a70b3924ee9689add612dadff80386 (1 file)
  records              1638e5ae81bda8cf53cffe4fa89eaddd7e2673e2613c02bde6a899792a93d85d (1 file)
  cluster_mapping      66b343397db8f928535832b31e48260ce67b8dadfcb83566e9b123e8a6f38a33 (1 file)
  populations          6c870a9203332e83edbda0cbaefe7b918288bb7d62746e892618d064abbe5a26 (1 file)

circuit records: 490 usable; 137 molecular; 19 projection-specific
rows: 471 somatic with a signed manipulation; 433 in the primary analysis
  set aside, by background and circuit: asip:gain_of_function:global on arc_agrp 1, asip:gain_of_function:global on arc_pomc 1, lep:loss_of_function:global on arc_agrp 8, lep:loss_of_function:global on arc_pomc 4, mixed on arc_agrp 5, mixed on lha_orexin 2, npy:loss_of_function:cell_type_restricted on arc_agrp 16; non_monotonic 1
  condition slots known in the primary rows: sex 50%, diet 59%, age_group 13%, challenge 15%, background_genotype 3%; mean 1.39 of 5
  molecular cells measured in two or more circuits: 9 (sirt1, trp53); compared after the independence rule: 0
  genes: D-091/D-093 key scored 102 findings, 24 valence-unclear, 7 holding a projection-specific record; corrected key scores 75, of which 70 have exactly the record set of an originally scored finding
  connectivity: D-091/D-093 key scored 82 findings, 16 valence-unclear, 5 holding a projection-specific record; corrected key scores 65, of which 60 have exactly the record set of an originally scored finding

connectivity variant, ancestors of HY removed: Spearman 1.000 against D-093's feature over 66 pairs
connectivity variant, grey matter only: Spearman 0.975 against D-093's feature over 66 pairs

genes
  circuits with findings      : 15, of which 14 in a comparison
  pairs eligible              : 102
  pairs with a comparison     : 40
  comparisons                 : 73   (22 distinct readout x manipulation cells)
  comparisons per pair        : median 1, max 4
  by manipulation class       : {'circuit_activity': 73}
  involving lha_orexin/arc_agrp: 48 of 73
  with a rat record on a side : 9 of 73
  circuit appearances         : lha_orexin 31, arc_agrp 21, arc_pomc 16, pvh_crh 14, tmn_histamine 11, pvh_trh 10, lha_mch 9, pvh_sim1 8, lha_vglut2 6, poa_galanin 6, lha_vgat 5, pvh_oxt 5, scn_vip 3, arc_ghrh 1
  most similar compared pairs : lha_vgat-pvh_sim1 rank 101/102 (1); arc_agrp-arc_pomc rank 97/102 (4); lha_vglut2-poa_galanin rank 95/102 (1); arc_agrp-arc_ghrh rank 94/102 (1); arc_pomc-lha_vgat rank 91/102 (1); pvh_crh-pvh_trh rank 90/102 (1); arc_pomc-pvh_sim1 rank 89/102 (2); arc_agrp-pvh_trh rank 86/102 (2)
  sensitivity comparisons     : a 66, b 60, c1 42, c2 52, e 74, f 84, g 55
  nearest-neighbour re-score: 199 findings, 75 scored, 48 with 2+ candidate neighbours
  pair test: two-sided 0.025, power counts passes in the modelled direction only. nearest neighbour: one-sided 0.05 (D-091's rule) / 0.0125
  size -- no link, 2000 datasets each:
    tau 0.0   pair 0.019   nearest 0.051 / 0.012
    tau 1.0   pair 0.025   nearest 0.049 / 0.017
  latent model -- similarity sets the correlation between circuits:
    tau 0.0 w 0.25   pair 0.02   nearest 0.08 / 0.02
    tau 0.0 w 0.50   pair 0.00   nearest 0.09 / 0.02
    tau 0.0 w 1.00   pair 0.03   nearest 0.12 / 0.03
    tau 1.0 w 0.25   pair 0.01   nearest 0.05 / 0.01
    tau 1.0 w 0.50   pair 0.01   nearest 0.07 / 0.00
    tau 1.0 w 1.00   pair 0.03   nearest 0.12 / 0.03
  family model -- a circuit does what its similarity family does, with probability pi:
    average k=n/4 sizes [12, 1, 1, 1]            align +0.65   pi 0.5   pair 0.07   nearest 0.09 / 0.02
    average k=n/4 sizes [12, 1, 1, 1]            align +0.65   pi 0.8   pair 0.22   nearest 0.37 / 0.17
    average k=n/4 sizes [12, 1, 1, 1]            align +0.65   pi 1.0   pair 0.68   nearest 0.92 / 0.57
    average k=n/3 sizes [10, 2, 1, 1, 1]         align +0.69   pi 0.5   pair 0.11   nearest 0.12 / 0.05
    average k=n/3 sizes [10, 2, 1, 1, 1]         align +0.69   pi 0.8   pair 0.36   nearest 0.28 / 0.08
    average k=n/3 sizes [10, 2, 1, 1, 1]         align +0.69   pi 1.0   pair 0.68   nearest 0.63 / 0.30
    average k=n/2 sizes [5, 4, 1, 1, 1, 1, 1, 1] align +0.53   pi 0.5   pair 0.02   nearest 0.15 / 0.05
    average k=n/2 sizes [5, 4, 1, 1, 1, 1, 1, 1] align +0.53   pi 0.8   pair 0.13   nearest 0.28 / 0.10
    average k=n/2 sizes [5, 4, 1, 1, 1, 1, 1, 1] align +0.53   pi 1.0   pair 0.20   nearest 0.52 / 0.28
    ward    k=n/4 sizes [5, 4, 3, 3]             align -0.05 * pi 0.5   pair 0.01   nearest 0.09 / 0.02
    ward    k=n/4 sizes [5, 4, 3, 3]             align -0.05 * pi 0.8   pair 0.02   nearest 0.14 / 0.05
    ward    k=n/4 sizes [5, 4, 3, 3]             align -0.05 * pi 1.0   pair 0.00   nearest 0.21 / 0.03
    ward    k=n/3 sizes [4, 3, 3, 3, 2]          align +0.15 * pi 0.5   pair 0.02   nearest 0.10 / 0.04
    ward    k=n/3 sizes [4, 3, 3, 3, 2]          align +0.15 * pi 0.8   pair 0.03   nearest 0.17 / 0.07
    ward    k=n/3 sizes [4, 3, 3, 3, 2]          align +0.15 * pi 1.0   pair 0.03   nearest 0.17 / 0.05
    ward    k=n/2 sizes [4, 3, 2, 2, 1, 1, 1, 1] align +0.49   pi 0.5   pair 0.03   nearest 0.10 / 0.05
    ward    k=n/2 sizes [4, 3, 2, 2, 1, 1, 1, 1] align +0.49   pi 0.8   pair 0.07   nearest 0.25 / 0.10
    ward    k=n/2 sizes [4, 3, 2, 2, 1, 1, 1, 1] align +0.49   pi 1.0   pair 0.16   nearest 0.37 / 0.18
    (* same-family tracks similarity rank below 0.3 on the compared pairs: not a model of function following similarity)
  independent comparisons -- agreement changes by delta from least to most similar pair:
    delta -0.8   pair 0.91
    delta -0.4   pair 0.28
    delta +0.2   pair 0.10
    delta +0.4   pair 0.32
    delta +0.6   pair 0.61
    delta +0.8   pair 0.92

connectivity
  circuits with findings      : 12, of which 11 in a comparison
  pairs eligible              : 66
  pairs with a comparison     : 29
  comparisons                 : 59   (18 distinct readout x manipulation cells)
  comparisons per pair        : median 2, max 4
  by manipulation class       : {'circuit_activity': 59}
  involving lha_orexin/arc_agrp: 41 of 59
  with a rat record on a side : 8 of 59
  circuit appearances         : lha_orexin 28, arc_agrp 17, arc_pomc 13, pvh_crh 13, vmh_sf1 11, tmn_histamine 11, lha_mch 8, lha_vgat 5, pvh_oxt 5, poa_galanin 5, scn_vip 2
  most similar compared pairs : lha_orexin-lha_vgat rank 64/66 (2); arc_pomc-vmh_sf1 rank 63/66 (2); lha_mch-lha_orexin rank 60/66 (4); lha_orexin-tmn_histamine rank 58/66 (4); arc_agrp-arc_pomc rank 57/66 (4); arc_pomc-lha_orexin rank 56/66 (3); arc_pomc-lha_mch rank 54/66 (1); lha_orexin-poa_galanin rank 49/66 (2)
  sensitivity comparisons     : a 52, b 49, c1 31, c2 42, e 60, f 70
  nearest-neighbour re-score: 184 findings, 65 scored, 40 with 2+ candidate neighbours
  pair test: two-sided 0.025, power counts passes in the modelled direction only. nearest neighbour: one-sided 0.05 (D-091's rule) / 0.0125
  size -- no link, 2000 datasets each:
    tau 0.0   pair 0.020   nearest 0.070 / 0.021
    tau 1.0   pair 0.029   nearest 0.060 / 0.024
  latent model -- similarity sets the correlation between circuits:
    tau 0.0 w 0.25   pair 0.02   nearest 0.09 / 0.04
    tau 0.0 w 0.50   pair 0.02   nearest 0.10 / 0.05
    tau 0.0 w 1.00   pair 0.03   nearest 0.16 / 0.05
    tau 1.0 w 0.25   pair 0.01   nearest 0.09 / 0.03
    tau 1.0 w 0.50   pair 0.01   nearest 0.11 / 0.03
    tau 1.0 w 1.00   pair 0.03   nearest 0.15 / 0.07
  family model -- a circuit does what its similarity family does, with probability pi:
    average k=n/4 sizes [9, 2, 1]                align +0.70   pi 0.5   pair 0.08   nearest 0.11 / 0.03
    average k=n/4 sizes [9, 2, 1]                align +0.70   pi 0.8   pair 0.30   nearest 0.33 / 0.15
    average k=n/4 sizes [9, 2, 1]                align +0.70   pi 1.0   pair 0.80   nearest 0.92 / 0.60
    average k=n/3 sizes [5, 4, 2, 1]             align +0.63   pi 0.5   pair 0.04   nearest 0.19 / 0.06
    average k=n/3 sizes [5, 4, 2, 1]             align +0.63   pi 0.8   pair 0.28   nearest 0.54 / 0.34
    average k=n/3 sizes [5, 4, 2, 1]             align +0.63   pi 1.0   pair 0.61   nearest 0.92 / 0.72
    average k=n/2 sizes [4, 4, 1, 1, 1, 1]       align +0.61   pi 0.5   pair 0.04   nearest 0.17 / 0.07
    average k=n/2 sizes [4, 4, 1, 1, 1, 1]       align +0.61   pi 0.8   pair 0.13   nearest 0.44 / 0.22
    average k=n/2 sizes [4, 4, 1, 1, 1, 1]       align +0.61   pi 1.0   pair 0.31   nearest 0.79 / 0.52
    ward    k=n/4 sizes [5, 4, 3]                align +0.54   pi 0.5   pair 0.07   nearest 0.20 / 0.11
    ward    k=n/4 sizes [5, 4, 3]                align +0.54   pi 0.8   pair 0.16   nearest 0.55 / 0.30
    ward    k=n/4 sizes [5, 4, 3]                align +0.54   pi 1.0   pair 0.37   nearest 0.87 / 0.61
    ward    k=n/3 sizes [5, 4, 2, 1]             align +0.63   pi 0.5   pair 0.08   nearest 0.24 / 0.12
    ward    k=n/3 sizes [5, 4, 2, 1]             align +0.63   pi 0.8   pair 0.25   nearest 0.53 / 0.32
    ward    k=n/3 sizes [5, 4, 2, 1]             align +0.63   pi 1.0   pair 0.57   nearest 0.89 / 0.69
    ward    k=n/2 sizes [4, 2, 2, 2, 1, 1]       align +0.37   pi 0.5   pair 0.01   nearest 0.10 / 0.04
    ward    k=n/2 sizes [4, 2, 2, 2, 1, 1]       align +0.37   pi 0.8   pair 0.05   nearest 0.13 / 0.07
    ward    k=n/2 sizes [4, 2, 2, 2, 1, 1]       align +0.37   pi 1.0   pair 0.02   nearest 0.20 / 0.08
    (* same-family tracks similarity rank below 0.3 on the compared pairs: not a model of function following similarity)
  independent comparisons -- agreement changes by delta from least to most similar pair:
    delta -0.8   pair 0.80
    delta -0.4   pair 0.24
    delta +0.2   pair 0.06
    delta +0.4   pair 0.17
    delta +0.6   pair 0.45
    delta +0.8   pair 0.82

data target: family model pi 0.8, every eligible pair compared on the same 5 readouts (power in the modelled direction; align in brackets)
  genes, the corpus's own 15 circuits (510 comparisons): average n/4 0.44 (+0.61)  average n/3 0.70 (+0.66)  average n/2 0.68 (+0.56)  ward n/4 0.18 (+0.21)  ward n/3 0.24 (+0.34)  ward n/2 0.29 (+0.45)
  connectivity, the corpus's own 12 circuits (330 comparisons): average n/4 0.74 (+0.62)  average n/3 0.87 (+0.65)  average n/2 0.75 (+0.62)  ward n/4 0.82 (+0.56)  ward n/3 0.84 (+0.65)  ward n/2 0.47 (+0.53)
  random HypoMap cell types, 15 circuits (525 comparisons): average n/4 0.90 (+0.74)  average n/3 0.86 (+0.70)  average n/2 0.50 (+0.55)  ward n/4 0.75 (+0.67)  ward n/3 0.71 (+0.62)  ward n/2 0.49 (+0.50)
  random HypoMap cell types, 30 circuits (2175 comparisons): average n/4 0.99 (+0.65)  average n/3 0.95 (+0.59)  average n/2 0.79 (+0.44)  ward n/4 0.85 (+0.54)  ward n/3 0.90 (+0.49)  ward n/2 0.61 (+0.39)
  random HypoMap cell types, 60 circuits (8850 comparisons): average n/4 1.00 (+0.56)  average n/3 0.98 (+0.47)  average n/2 0.84 (+0.34)  ward n/4 0.99 (+0.43)  ward n/3 0.87 (+0.38)  ward n/2 0.70 (+0.29)
  mean pairwise gene cosine: the corpus's 15 circuits 0.255; 15 random cell types 0.016
```
