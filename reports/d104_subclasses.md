# D-104 `tables` output

Cell types, positions and fingerprints for D-104, by `PYTHONPATH=src python3 analysis/d104_epiretro.py tables`. Neurons per experiment and subclass are label counts; no gene-to-projection statistic is computed.

```
passing hypothalamic neurons 1538; subclasses with >= 20 of them: 22; compared (>= 50% of MERFISH cells in HY): 19, 1112 neurons, 958 in the 14 primary experiments
experiments (target, sex, neurons of the compared subclasses, astrocyte share of all sorted nuclei):
  Jf6H7J8K9L   MOB  female   27    45% of 99
  Jm6H7J8K9L   MOB  male     45    20% of 90
  Nm6H7J8K9L   AMY  male     25    62% of 111
  Pf6H7J8K9L   STR  female   45    43% of 132
  Pm6H7J8K9L   STR  male     43    39% of 133
  Qf6H7J8K9L   PAL  female   25    63% of 119
  Rf6H7J8K9L   TH   female   78    28% of 151
  Rm6H7J8K9L   TH   male     97    11% of 161
  Tf6H7J8K9L   SC   female  152     1% of 177
  Tm6H7J8K9L   SC   male     94     3% of 109
  Uf6H7J8K9L   VTA  female  104     1% of 144
  Vf6H7J8K9L   P    female   92    16% of 123
  Vm6H7J8K9L   P    male     43     9% of 57
  Wf6H7J8K9L   MY   female   41    39% of 79
  Wm6H7J8K9L   MY   male     81     5% of 92
  Zf6H7J8K9L   PFC  female   37    37% of 118
  Zm6H7J8K9L   PFC  male     83     1% of 155
astrocytes among all 2542 sorted hypothalamic nuclei: 933 (37%)
subclass                             neurons primary clean MERFISH HY HY share  zone  cells by zone PVZ/PVR/MEZ/LZ/unassigned
AHN Onecut3 Gaba                          26      22     1       1198      96%   MEZ  38/141/447/281/291
AHN-RCH-LHA Otp Fezf1 Glut                29      21    23       2420      93%   MEZ  52/108/1278/391/591
AHN-SBPV-PVHd Pdrm12 Gaba                 35      29    27       2128      98%   MEZ  50/226/1228/390/234
AVPV-MEPO-SFO Tbr1 Glut                   36      34    13       1060      77%   PVR  0/939/4/21/96
BST-MPN Six3 Nrgn Gaba                    31      26    27       1644      56%   MEZ  59/167/977/65/376
DMH Hmx2 Gaba                             42      26    14       2766      96%   PVR  175/887/218/662/824
DMH-LHA Gsx1 Gaba                         23      17    14       2087      95%    LZ  250/399/120/846/472
DMH-LHA Vgll2 Glut                       170     164   148        583      95%    LZ  12/89/39/294/149
LHA Barhl2 Glut                           58      46    43        940      80%    LZ  28/23/13/505/371
LHA-AHN-PVH Otp Trh Glut                  31      28    23       1408      84%   PVR  139/366/288/332/283
MPN-MPO-PVpo Hmx2 Glut                    21      15    17       1482      93%   PVR  7/834/476/46/119
NDB-SI-MA-STRv Lhx8 Gaba                  27      24    27        709       7%    LZ  0/136/10/301/262  (not hypothalamic)
PH-ant-LHA Otp Bsx Glut                   59      55    24       1814      93%   MEZ  24/381/472/470/467
PVH-SO-PVa Otp Glut                       55      39    33       1810      88%   PVZ  1070/134/90/182/334
PVHd-DMH Lhx6 Gaba                        31      25     1       1016      96%    LZ  14/67/281/339/315
PVHd-SBPV Six3 Prox1 Gaba                 26      17    13       1325      95%   MEZ  60/198/446/168/453
PVpo-VMPO-MPN Hmx2 Gaba                  101      87    63       3234      86%   PVR  65/1695/862/163/449
RE-Xi Nox4 Glut                          139     125   138        142       4%    LZ  7/0/8/69/58  (not hypothalamic)
RT-ZI Gnb3 Gaba                           43      36    33        259       3%    LZ  0/3/6/222/28  (not hypothalamic)
STN-PSTN Pitx2 Glut                       92      92    78       2329      84%    LZ  0/6/205/2004/114
VMH Nr5a1 Glut                           135     119   125       4032      99%   MEZ  254/29/2720/427/602
ZI Pax6 Gaba                             111      96    75       4909      81%    LZ  46/11/324/3616/912
of the compared subclasses' neurons, 12.9% carry a label other than their methylation cluster's majority; 215 of 266 clusters are clean (>= 90% one subclass)
largest share of methylation clusters two compared subclasses have in common: 0.20
MERFISH panel genes found in the 10x data: 500; hypothalamic parcels: 45
```
