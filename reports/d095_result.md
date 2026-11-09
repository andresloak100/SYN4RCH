# D-095 `--run` output

The pre-registered test, run once by `python3 analysis/d095_pair_test.py --run`.

```
commit a71a6e511573ec669b42baa4367f3e3348672daf; every input matches its pinned SHA-256

PRIMARY (the only decisions: two-sided p < 0.025)
  genes                                n=73  r=-0.045  null 2.5/50/97.5 = -0.267/+0.003/+0.259  p=0.7381  -> not significant
                                       agreement: more-similar half 0.424 (33), less-similar half 0.425 (40), overall 0.425
  connectivity                         n=59  r=-0.011  null 2.5/50/97.5 = -0.289/+0.001/+0.301  p=0.9441  -> not significant
                                       agreement: more-similar half 0.351 (37), less-similar half 0.364 (22), overall 0.356

SENSITIVITY (no decision; read by the rule fixed in D-095)
 genes
  a  laboratory excluded               n=66  r=-0.065  null 2.5/50/97.5 = -0.286/+0.005/+0.277  p=0.6525
                                       agreement: more-similar half 0.414 (29), less-similar half 0.432 (37), overall 0.424
  b  same-nucleus pairs out            n=60  r=+0.041  null 2.5/50/97.5 = -0.306/+0.003/+0.312  p=0.8129
                                       agreement: more-similar half 0.519 (27), less-similar half 0.455 (33), overall 0.483
  c1 without lha_orexin                n=42  r=-0.046  null 2.5/50/97.5 = -0.339/+0.002/+0.340  p=0.7985
                                       agreement: more-similar half 0.524 (21), less-similar half 0.429 (21), overall 0.476
  c2 without arc_agrp                  n=52  r=-0.060  null 2.5/50/97.5 = -0.318/+0.002/+0.316  p=0.7517
                                       agreement: more-similar half 0.370 (27), less-similar half 0.520 (25), overall 0.442
  e  rat records dropped               n=74  r=-0.012  null 2.5/50/97.5 = -0.263/+0.002/+0.264  p=0.9233
                                       agreement: more-similar half 0.424 (33), less-similar half 0.415 (41), overall 0.419
  f  mutant and mixed backgrounds in   n=84  r=-0.029  null 2.5/50/97.5 = -0.254/+0.001/+0.251  p=0.8111
                                       agreement: more-similar half 0.425 (40), less-similar half 0.455 (44), overall 0.440
  g  tier-C broad markers out          n=55  r=-0.123  null 2.5/50/97.5 = -0.286/+0.004/+0.284  p=0.4192
                                       agreement: more-similar half 0.379 (29), less-similar half 0.462 (26), overall 0.418
  d  both-no-change dropped            n=70  r=-0.045  null 2.5/50/97.5 = -0.290/+0.004/+0.280  p=0.7621
                                       agreement: more-similar half 0.406 (32), less-similar half 0.395 (38), overall 0.400
 connectivity
  a  laboratory excluded               n=52  r=+0.025  null 2.5/50/97.5 = -0.320/+0.001/+0.325  p=0.9031
                                       agreement: more-similar half 0.344 (32), less-similar half 0.350 (20), overall 0.346
  b  same-nucleus pairs out            n=49  r=+0.049  null 2.5/50/97.5 = -0.326/+0.002/+0.338  p=0.8065
                                       agreement: more-similar half 0.379 (29), less-similar half 0.400 (20), overall 0.388
  c1 without lha_orexin                n=31  r=-0.074  null 2.5/50/97.5 = -0.370/-0.001/+0.384  p=0.7519
                                       agreement: more-similar half 0.333 (15), less-similar half 0.438 (16), overall 0.387
  c2 without arc_agrp                  n=42  r=+0.020  null 2.5/50/97.5 = -0.370/+0.003/+0.381  p=0.9377
                                       agreement: more-similar half 0.357 (28), less-similar half 0.357 (14), overall 0.357
  e  rat records dropped               n=60  r=-0.067  null 2.5/50/97.5 = -0.290/+0.000/+0.301  p=0.6817
                                       agreement: more-similar half 0.333 (39), less-similar half 0.381 (21), overall 0.350
  f  mutant and mixed backgrounds in   n=70  r=+0.006  null 2.5/50/97.5 = -0.291/+0.000/+0.298  p=0.9743
                                       agreement: more-similar half 0.372 (43), less-similar half 0.370 (27), overall 0.371
  d  both-no-change dropped            n=57  r=-0.053  null 2.5/50/97.5 = -0.312/+0.002/+0.320  p=0.7591
                                       agreement: more-similar half 0.314 (35), less-similar half 0.364 (22), overall 0.333


SECONDARY (no decision): both features on the 11 circuits both place that have findings, pairs eligible under both
  genes                                n=48  r=-0.196  null 2.5/50/97.5 = -0.314/+0.003/+0.303  p=0.2442
                                       agreement: more-similar half 0.333 (24), less-similar half 0.458 (24), overall 0.396
  connectivity                         n=48  r=-0.106  null 2.5/50/97.5 = -0.328/+0.002/+0.319  p=0.5635
                                       agreement: more-similar half 0.355 (31), less-similar half 0.471 (17), overall 0.396

CORRECTION (no decision): D-091 / D-093's nearest-neighbour test on the original and the corrected key
  genes (D-092), original key        circuits=15 findings=102 acc=0.422  null med=0.412 95th=0.461  p=0.381
  genes (D-092), corrected key       circuits=15 findings=75  acc=0.440  null med=0.413 95th=0.467  p=0.297
  connectivity (D-094), original key circuits=12 findings=82  acc=0.427  null med=0.378 95th=0.427  p=0.057
  connectivity (D-094), corrected key circuits=12 findings=65  acc=0.400  null med=0.338 95th=0.400  p=0.085
  genes, shared circuits, original key circuits=11 findings=79  acc=0.418  null med=0.392 95th=0.443  p=0.313
  genes, shared circuits, corrected key circuits=11 findings=62  acc=0.387  null med=0.371 95th=0.435  p=0.397
  connectivity, shared circuits, original key circuits=11 findings=79  acc=0.418  null med=0.392 95th=0.443  p=0.313
  connectivity, shared circuits, corrected key circuits=11 findings=62  acc=0.387  null med=0.371 95th=0.435  p=0.397

every compared pair, by genes similarity:
    lha_vgat        pvh_sim1        sim +0.733   agree 0 of 1
    arc_agrp        arc_pomc        sim +0.544   agree 1 of 4
    lha_vglut2      poa_galanin     sim +0.510   agree 1 of 1
    arc_agrp        arc_ghrh        sim +0.475   agree 0 of 1
    arc_pomc        lha_vgat        sim +0.463   agree 0 of 1
    pvh_crh         pvh_trh         sim +0.439   agree 0 of 1
    arc_pomc        pvh_sim1        sim +0.433   agree 2 of 2
    arc_agrp        pvh_trh         sim +0.396   agree 2 of 2
    arc_pomc        pvh_trh         sim +0.372   agree 1 of 2
    arc_agrp        lha_vgat        sim +0.365   agree 1 of 1
    arc_agrp        pvh_sim1        sim +0.364   agree 0 of 2
    lha_vglut2      tmn_histamine   sim +0.361   agree 1 of 1
    arc_pomc        lha_vglut2      sim +0.361   agree 1 of 1
    arc_agrp        pvh_crh         sim +0.337   agree 1 of 1
    lha_orexin      pvh_trh         sim +0.331   agree 2 of 3
    lha_orexin      pvh_sim1        sim +0.314   agree 0 of 2
    arc_agrp        lha_vglut2      sim +0.310   agree 0 of 1
    lha_orexin      poa_galanin     sim +0.303   agree 0 of 2
    lha_orexin      lha_vglut2      sim +0.280   agree 0 of 1
    arc_pomc        lha_orexin      sim +0.260   agree 1 of 3
    lha_mch         lha_orexin      sim +0.257   agree 0 of 4
    lha_orexin      lha_vgat        sim +0.254   agree 1 of 2
    lha_vglut2      pvh_crh         sim +0.233   agree 1 of 1
    arc_pomc        pvh_oxt         sim +0.225   agree 0 of 1
    poa_galanin     pvh_crh         sim +0.223   agree 1 of 1
    lha_orexin      pvh_crh         sim +0.211   agree 2 of 4
    pvh_crh         scn_vip         sim +0.209   agree 0 of 1
    arc_agrp        pvh_oxt         sim +0.202   agree 0 of 3
    arc_agrp        lha_orexin      sim +0.191   agree 2 of 4
    lha_orexin      pvh_oxt         sim +0.179   agree 0 of 1
    lha_mch         pvh_sim1        sim +0.175   agree 0 of 1
    arc_pomc        pvh_crh         sim +0.163   agree 1 of 1
    lha_orexin      tmn_histamine   sim +0.150   agree 3 of 4
    arc_pomc        lha_mch         sim +0.143   agree 0 of 1
    pvh_trh         scn_vip         sim +0.133   agree 0 of 1
    lha_orexin      scn_vip         sim +0.093   agree 0 of 1
    pvh_crh         tmn_histamine   sim +0.081   agree 3 of 4
    lha_mch         pvh_trh         sim +0.045   agree 1 of 1
    arc_agrp        lha_mch         sim +0.044   agree 1 of 2
    poa_galanin     tmn_histamine   sim -0.097   agree 1 of 2

every compared pair, by connectivity similarity:
    lha_orexin      lha_vgat        sim +0.755   agree 1 of 2
    arc_pomc        vmh_sf1         sim +0.748   agree 2 of 2
    lha_mch         lha_orexin      sim +0.705   agree 0 of 4
    lha_orexin      tmn_histamine   sim +0.651   agree 3 of 4
    arc_agrp        arc_pomc        sim +0.651   agree 1 of 4
    arc_pomc        lha_orexin      sim +0.636   agree 1 of 3
    arc_pomc        lha_mch         sim +0.613   agree 0 of 1
    lha_orexin      poa_galanin     sim +0.556   agree 0 of 2
    lha_mch         vmh_sf1         sim +0.537   agree 0 of 1
    arc_agrp        lha_orexin      sim +0.520   agree 2 of 4
    lha_orexin      vmh_sf1         sim +0.510   agree 0 of 3
    lha_orexin      pvh_crh         sim +0.488   agree 2 of 4
    poa_galanin     tmn_histamine   sim +0.487   agree 1 of 2
    lha_orexin      pvh_oxt         sim +0.460   agree 0 of 1
    arc_pomc        lha_vgat        sim +0.457   agree 0 of 1
    tmn_histamine   vmh_sf1         sim +0.450   agree 0 of 1
    arc_agrp        vmh_sf1         sim +0.441   agree 0 of 2
    arc_pomc        pvh_crh         sim +0.440   agree 1 of 1
    arc_pomc        pvh_oxt         sim +0.436   agree 0 of 1
    lha_vgat        vmh_sf1         sim +0.424   agree 0 of 1
    arc_agrp        lha_mch         sim +0.412   agree 1 of 2
    pvh_crh         vmh_sf1         sim +0.402   agree 0 of 1
    poa_galanin     pvh_crh         sim +0.382   agree 1 of 1
    arc_agrp        lha_vgat        sim +0.371   agree 1 of 1
    arc_agrp        pvh_oxt         sim +0.361   agree 0 of 3
    pvh_crh         tmn_histamine   sim +0.354   agree 3 of 4
    arc_agrp        pvh_crh         sim +0.302   agree 1 of 1
    pvh_crh         scn_vip         sim +0.231   agree 0 of 1
    lha_orexin      scn_vip         sim +0.213   agree 0 of 1

```
