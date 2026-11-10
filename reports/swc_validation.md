# Jiao 2025 SWC sample: validation

Written by `analysis/swc_validate.py` (python3 -I). Raw files in `data/raw/jiao2025/swc/` are read, never changed; hashes and sources are in `data/raw/jiao2025/swc/MANIFEST.tsv`. Machine-readable copy: `data/interim/swc_validation.json`.

- 41 files, 37 distinct Jiao neurons: 18 Zenodo files (CC BY 4.0, mirrored-left frame) and 23 portal files (copyright line only, all rights reserved; BSDC record of the same data: CC BY-NC-ND 4.0; true-hemisphere frame); 4 neurons are in both. No MouseLight file is included.
- Units: micrometres, Allen CCFv3 2017; SWC x, y, z = AP, DV, ML. Midline ML = 5,700 um. Segment type = child node type. 'type0' = type 0 (unlabelled), 'other' = undocumented types 10, 11 and 12; neither is relabelled.
- Lengths are summed straight segments between consecutive nodes. Branch points: non-root nodes with 2 or more children (soma children listed separately). Tips: non-root nodes with no children.
- Soma mesh tests (root, HY, nucleus): exact ray parity against every triangle along 15 fixed random directions (seed 0); inside = at least 8 of 15 rays odd. The vote is reported; any non-unanimous vote is listed as an anomaly. Distances to a mesh are exact point-to-triangle distances to the surface (the nearest-vertex distance, an upper bound, is also kept in the JSON).
- Node tests (every node against the root mesh): majority vote of three grid parity tests with rays along AP, DV and ML (20 um column grid), because a single axis is fragile where a mesh has holes. A node within about 14 um of the surface can fall either way; nodes where the three axes disagree are counted per file (`root_vote_split_nodes` in the JSON).
- Line endings: all Zenodo files use CRLF and all portal files LF (counted per file below). This alone explains the byte difference between the two copies of one neuron (one extra byte per line); the parsed content is identical.
- Contralateral fraction: share of the segment length on the side of the midline opposite the soma, in the file's own frame (a reflection does not change it).

## Mesh check (data/raw/atlas/meshes)

Coordinate frame: OBJ vertices are AP, DV, ML in micrometres, CCFv3, both hemispheres (root ML 486 to 10,891; HY ML 3,158 to 8,219, centred on 5,700), the same frame as the portal SWCs.

| mesh | id | vertices | faces | open edges | non-manifold edges | grid columns with odd hits (AP/DV/ML rays) | min (AP, DV, ML) | max |
|---|---|---|---|---|---|---|---|---|
| AHN | 88 | 713 | 1418 | 0 | 0 | 0/0/0 | [5551.0, 5777.0, 4779.0] | [6654.0, 6805.0, 6606.0] |
| ARH | 223 | 349 | 694 | 0 | 0 | 0/0/0 | [6368.0, 6732.0, 5195.0] | [7440.0, 7265.0, 6181.0] |
| AVPV | 272 | 262 | 516 | 6 | 0 | 7/20/7 | [4588.0, 6072.0, 5266.0] | [5208.0, 6833.0, 6100.0] |
| DMH | 830 | 477 | 946 | 0 | 0 | 0/0/0 | [6678.0, 5740.0, 4996.0] | [7469.0, 6717.0, 6389.0] |
| HY | 1097 | 6217 | 12329 | 89 | 0 | 87/100/94 | [4510.0, 3108.0, 3158.0] | [8464.0, 7428.0, 8219.0] |
| LHA | 194 | 1672 | 3336 | 0 | 0 | 0/0/0 | [5784.0, 5329.0, 3921.0] | [7842.0, 6599.0, 7460.0] |
| LPO | 226 | 779 | 1550 | 0 | 0 | 0/0/0 | [4549.0, 5525.0, 4179.0] | [5790.0, 6670.0, 7181.0] |
| ME | 10671 | 157 | 310 | 0 | 0 | 0/0/0 | [6680.0, 7146.0, 5405.0] | [7603.0, 7422.0, 5964.0] |
| MEPO | 452 | 137 | 270 | 0 | 0 | 0/0/0 | [4872.0, 4535.0, 5615.0] | [5353.0, 6334.0, 5765.0] |
| MPO | 523 | 893 | 1782 | 0 | 0 | 0/0/0 | [4697.0, 5742.0, 4759.0] | [5812.0, 6742.0, 6604.0] |
| PH | 946 | 723 | 1442 | 0 | 0 | 0/0/0 | [6998.0, 4755.0, 4899.0] | [8020.0, 6314.0, 6466.0] |
| PMv | 1004 | 292 | 576 | 0 | 0 | 0/0/0 | [7328.0, 6360.0, 4765.0] | [7820.0, 6877.0, 6597.0] |
| PVH | 38 | 250 | 496 | 0 | 0 | 0/0/0 | [5672.0, 5323.0, 5304.0] | [6432.0, 6025.0, 6074.0] |
| SCH | 286 | 138 | 272 | 0 | 0 | 0/0/0 | [5667.0, 6678.0, 5397.0] | [6155.0, 6993.0, 5990.0] |
| SO | 390 | 170 | 332 | 0 | 0 | 0/0/0 | [5576.0, 6534.0, 4174.0] | [6151.0, 6941.0, 7193.0] |
| TM | 557 | 327 | 644 | 4 | 0 | 5/5/1 | [7352.0, 6304.0, 4403.0] | [7739.0, 6938.0, 6957.0] |
| VMH | 693 | 554 | 1100 | 0 | 0 | 0/0/0 | [6333.0, 6306.0, 4873.0] | [7293.0, 7063.0, 6491.0] |
| VMPO | 576073699 | 129 | 250 | 0 | 0 | 0/0/0 | [5242.0, 6678.0, 5099.0] | [5556.0, 6945.0, 6267.0] |
| ZI | 797 | 1927 | 3846 | 0 | 0 | 0/0/0 | [5872.0, 4099.0, 3165.0] | [8204.0, 5571.0, 8216.0] |
| root | 997 | 49324 | 98638 | 6 | 0 | 0/0/0 | [-17.0, 134.0, 486.0] | [13192.0, 7564.0, 10891.0] |

## Per file

| file | frame | line | soma label (xlsx) | hemisph. (index) | soma side | nodes | types | axon um (type 2) | type0 um | dendrite um | branch pts | tips | soma children | tree ok | soma vs xlsx um | in root | in HY | nucleus test | contra axon / type0 | nodes outside root |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| zenodo/195038_041_Oxt_PVH.swc | zenodo | Oxt | PVH | Left | Left | 9576 | 0:9575 1:1 | 0 | 24057 | 0 | 13 | 14 | 1 | yes | 0.00 | yes 15/15 | yes 15/15 | PVH: in 15/15 | - / 0.000 | 0 |
| zenodo/201232_040_Tac1_HY.swc | zenodo | Tac1 | HY | Right | Left | 36312 | 1:1 2:35441 3:867 10:3 | 91879 | 0 | 2075 | 96 | 102 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | HY: in 15/15 | 0.051 / - | 0 |
| zenodo/201232_056_Tac1_MMm.swc | zenodo | Tac1 | MMm | Right | Left | 6875 | 1:1 2:5746 3:1124 10:4 | 14493 | 0 | 2970 | 73 | 78 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | HY (ancestor of MMm): in 15/15 | 0.000 / - | 0 |
| zenodo/201894_001_Adcyap1_PSTN.swc | zenodo | Adcyap1 | PSTN | Right | Left | 35511 | 1:1 2:34067 3:1438 10:5 | 81600 | 0 | 3599 | 804 | 810 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | HY (ancestor of PSTN): in 15/15 | 0.000 / - | 0 |
| zenodo/201922_073_Nts_MPO.swc | zenodo | Nts | MPO | Right | Left | 52470 | 0:50958 1:1 3:1487 10:20 11:4 | 0 | 126068 | 4655 | 1009 | 1049 | 5 | yes | 0.00 | yes 15/15 | yes 15/15 | MPO: in 15/15 | - / 0.121 | 0 |
| zenodo/202559_043_Tac2_HY.swc | zenodo | Tac2 | HY | Right | Left | 12515 | 1:1 2:11650 3:842 10:22 | 29978 | 0 | 2315 | 243 | 247 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | HY: in 15/15 | 0.050 / - | 0 |
| zenodo/202617_025_Sst_ZI.swc | zenodo | Sst | ZI | Left | Left | 11248 | 1:1 2:10250 3:992 10:1 11:1 12:3 | 24489 | 0 | 2560 | 156 | 160 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | ZI: in 15/15 | 0.000 / - | 0 |
| zenodo/202618_143_Sst_LHA.swc | zenodo | Sst | LHA | Right | Left | 12574 | 0:10959 1:1 3:1607 10:5 11:2 | 0 | 29599 | 4007 | 206 | 212 | 6 | yes | 0.00 | yes 15/15 | yes 15/15 | LHA: in 15/15 | - / 0.238 | 0 |
| zenodo/202620_016_Sst_PeF.swc | zenodo | Sst | PeF | Right | Left | 35079 | 1:1 2:32200 3:2872 10:2 11:1 12:3 | 79254 | 0 | 7795 | 695 | 700 | 3 | yes | 0.00 | yes 15/15 | yes 14/15 | HY (ancestor of PeF): in 14/15 | 0.034 / - | 0 |
| zenodo/210796_002_Pdyn_VMH.swc | zenodo | Pdyn | VMH | Left | Left | 14681 | 1:1 2:13581 3:1092 10:2 11:1 12:4 | 38290 | 0 | 3230 | 175 | 177 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | VMH: in 15/15 | 0.171 / - | 0 |
| zenodo/211266_057_Pdyn_HY.swc | zenodo | Pdyn | HY | Right | Left | 22255 | 0:22252 1:1 12:2 | 0 | 58329 | 0 | 488 | 500 | 1 | yes | 0.00 | yes 15/15 | yes 15/15 | HY: in 15/15 | - / 0.970 | 0 |
| zenodo/211353_070_Adcyap1_AHN.swc | zenodo | Adcyap1 | AHN | Left | Left | 25519 | 1:1 2:24403 3:1087 10:28 | 64711 | 0 | 3307 | 332 | 334 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | AHN: in 15/15 | 0.010 / - | 0 |
| zenodo/211461_026_Nts_PH.swc | zenodo | Nts | PH | Right | Left | 36243 | 0:34773 1:1 3:1468 12:1 | 0 | 88264 | 3720 | 849 | 876 | 5 | yes | 0.00 | yes 15/15 | yes 15/15 | PH: in 15/15 | - / 0.000 | 0 |
| zenodo/211465_137_Tac2_ZI.swc | zenodo | Tac2 | ZI | Left | Left | 21276 | 0:20350 1:1 3:921 10:4 | 0 | 49595 | 2334 | 463 | 468 | 5 | yes | 0.00 | yes 15/15 | yes 15/15 | ZI: in 15/15 | - / 0.032 | 0 |
| zenodo/211595_037_Sst_HY.swc | zenodo | Sst | HY | Left | Left | 12304 | 0:11097 1:1 3:1206 | 0 | 30877 | 3194 | 127 | 133 | 6 | yes | 0.00 | yes 15/15 | yes 15/15 | HY: in 15/15 | - / 0.000 | 0 |
| zenodo/211595_063_Sst_HY.swc | zenodo | Sst | HY | Right | Left | 14414 | 0:13421 1:1 3:991 10:1 | 0 | 37622 | 2858 | 246 | 253 | 4 | yes | 0.00 | yes 15/15 | yes 15/15 | HY: in 15/15 | - / 0.032 | 0 |
| zenodo/211889_096_Pdyn_HY.swc | zenodo | Pdyn | HY | Right | Left | 9701 | 1:1 2:9475 3:223 12:2 | 23957 | 0 | 584 | 67 | 69 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | HY: in 15/15 | 0.073 / - | 0 |
| zenodo/211890_005_Pdyn_VMH.swc | zenodo | Pdyn | VMH | Right | Left | 40816 | 0:39738 1:1 3:1075 10:2 | 0 | 99695 | 3353 | 525 | 531 | 6 | yes | 0.00 | yes 15/15 | yes 15/15 | VMH: in 15/15 | - / 0.222 | 0 |
| portal/192710_002.swc | portal | Pomc | TU | Left | Left | 19306 | 0:18570 1:1 3:733 10:1 12:1 | 0 | 49058 | 2080 | 241 | 269 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | HY (ancestor of TU): in 15/15 | - / 0.000 | 0 |
| portal/195038_041.swc | portal | Oxt | PVH | Left | Left | 9576 | 0:9575 1:1 | 0 | 24057 | 0 | 13 | 14 | 1 | yes | 0.00 | yes 15/15 | yes 15/15 | PVH: in 15/15 | - / 0.000 | 0 |
| portal/201232_056.swc | portal | Tac1 | MMm | Right | Right | 6875 | 1:1 2:5746 3:1124 10:4 | 14493 | 0 | 2970 | 73 | 78 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | HY (ancestor of MMm): in 15/15 | 0.000 / - | 0 |
| portal/202457_066.swc | portal | Pomc | ARH | Right | Right | 12967 | 0:11783 1:1 3:1182 10:1 | 0 | 28770 | 3395 | 92 | 100 | 5 | yes | 0.00 | yes 15/15 | yes 15/15 | ARH: in 15/15 | - / 0.960 | 0 |
| portal/202813_035.swc | portal | Agrp | ARH | Left | Left | 11814 | 0:11458 1:1 3:354 10:1 | 0 | 29004 | 914 | 133 | 136 | 3 | yes | 0.00 | yes 15/15 | yes 15/15 | ARH: in 15/15 | - / 0.003 | 0 |
| portal/202813_009.swc | portal | Agrp | ARH | Right | Right | 11384 | 0:11042 1:1 3:340 10:1 | 0 | 30735 | 915 | 44 | 47 | 3 | yes | 0.00 | yes 15/15 | yes 15/15 | ARH: in 15/15 | - / 0.000 | 0 |
| portal/202825_072.swc | portal | Avp | PVH | Left | Left | 8429 | 0:7817 1:1 3:611 | 0 | 19060 | 1635 | 24 | 28 | 4 | yes | 0.00 | yes 15/15 | yes 15/15 | PVH: OUT (33.3 um from surface) 0/15 | - / 0.000 | 522 |
| portal/202825_065.swc | portal | Avp | PVH | Right | Right | 8478 | 0:7866 1:1 3:611 | 0 | 19762 | 1781 | 65 | 68 | 3 | yes | 0.00 | yes 15/15 | yes 15/15 | PVH: OUT (25.3 um from surface) 0/15 | - / 0.000 | 0 |
| portal/212248_071.swc | portal | Crh | PVH | Right | Right | 3337 | 0:2730 1:1 3:606 | 0 | 6998 | 1657 | 45 | 48 | 3 | yes | 0.00 | yes 15/15 | yes 15/15 | PVH: OUT (8.9 um from surface) 0/15 | - / 0.160 | 0 |
| portal/201221_022.swc | portal | Nts | LHA | Left | Left | 27132 | 0:26326 1:1 3:784 10:11 11:7 12:3 | 0 | 64113 | 2313 | 422 | 429 | 3 | yes | 0.00 | yes 15/15 | yes 15/15 | LHA: in 15/15 | - / 0.118 | 0 |
| portal/210105_001.swc | portal | Orexin | DMH | Left | Left | 20489 | 0:19863 1:1 3:624 10:1 | 0 | 53949 | 1807 | 71 | 77 | 4 | yes | 0.00 | yes 15/15 | yes 15/15 | DMH: in 15/15 | - / 0.319 | 0 |
| portal/210104_007.swc | portal | Orexin | LHA | Right | Right | 26984 | 1:1 2:26225 3:755 10:2 11:1 | 67405 | 0 | 2036 | 64 | 68 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | LHA: in 15/15 | 0.000 / - | 378 |
| portal/195039_009.swc | portal | Oxt | PVH | Right | Right | 2593 | 0:1952 1:1 3:640 | 0 | 5222 | 1500 | 11 | 14 | 3 | yes | 0.66 | yes 15/15 | yes 15/15 | PVH: in 15/15 | - / 0.036 | 10 |
| portal/210073_044.swc | portal | Penk | VMH | Left | Left | 33255 | 0:32343 1:1 3:911 | 0 | 79133 | 2557 | 611 | 616 | 5 | yes | 0.00 | yes 15/15 | yes 15/15 | VMH: in 15/15 | - / 0.025 | 0 |
| portal/202170_016.swc | portal | Penk | VMH | Right | Right | 23989 | 1:1 2:23385 3:591 10:8 11:3 12:1 | 59786 | 0 | 1583 | 510 | 525 | 3 | yes | 0.00 | yes 15/15 | yes 15/15 | VMH: in 15/15 | 0.498 / - | 0 |
| portal/210007_142.swc | portal | Pmch | LHA | Left | Left | 13225 | 1:1 2:11427 3:1794 12:3 | 29475 | 0 | 4774 | 281 | 283 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | LHA: in 15/15 | 0.021 / - | 0 |
| portal/210951_004.swc | portal | Pmch | LHA | Right | Right | 16520 | 0:13980 1:1 3:2532 10:7 | 0 | 36876 | 6560 | 162 | 164 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | LHA: in 15/15 | - / 0.073 | 728 |
| portal/201230_079.swc | portal | Tac1 | VMH | Left | Left | 38246 | 1:1 2:36911 3:1306 10:28 | 92904 | 0 | 3621 | 414 | 416 | 2 | yes | 1.92 | yes 15/15 | yes 15/15 | VMH: in 15/15 | 0.063 / - | 0 |
| portal/212202_026.swc | portal | Trh | PVHd | Left | Left | 22542 | 1:1 2:21750 3:786 12:5 | 55821 | 0 | 2078 | 334 | 336 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | HY (ancestor of PVHd): in 15/15 | 0.000 / - | 364 |
| portal/212202_011.swc | portal | Trh | PVHd | Right | Right | 9952 | 1:1 2:8223 3:1727 10:1 | 22128 | 0 | 5135 | 151 | 154 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | HY (ancestor of PVHd): in 15/15 | 0.000 / - | 201 |
| portal/202824_062.swc | portal | Vip | MPN | Right | Right | 28318 | 1:1 2:27040 3:1250 10:5 11:13 12:9 | 69943 | 0 | 3547 | 412 | 414 | 2 | yes | 0.00 | yes 15/15 | yes 15/15 | HY (ancestor of MPN): in 15/15 | 0.135 / - | 0 |
| portal/211461_026.swc | portal | Nts | PH | Right | Right | 36243 | 0:34773 1:1 3:1468 12:1 | 0 | 88264 | 3720 | 849 | 876 | 5 | yes | 0.00 | yes 15/15 | yes 15/15 | PH: in 15/15 | - / 0.000 | 0 |
| portal/202618_143.swc | portal | Sst | LHA | Right | Right | 12574 | 0:10959 1:1 3:1607 10:5 11:2 | 0 | 29599 | 4007 | 206 | 212 | 6 | yes | 0.00 | yes 15/15 | yes 15/15 | LHA: in 15/15 | - / 0.238 | 0 |

Bounding boxes (um, AP DV ML) and maximum segment length:

| file | min | max | max segment um | radius range (not used) | line ending | root-test vote-split nodes |
|---|---|---|---|---|---|---|
| zenodo/195038_041_Oxt_PVH.swc | [4967.7, 2442.6, 3813.9] | [10556.2, 7014.1, 5536.2] | 22.9 | 0.0039 to 55.9 | CRLF (9576) | 0 |
| zenodo/201232_040_Tac1_HY.swc | [3343.7, 424.1, 2397.4] | [8973.9, 6900.1, 6095.7] | 38.5 | 0.0039 to 48.5 | CRLF (36312) | 0 |
| zenodo/201232_056_Tac1_MMm.swc | [5990.0, 4168.1, 4764.6] | [10091.8, 6255.5, 5706.4] | 17.5 | 0.0039 to 50.6 | CRLF (6875) | 0 |
| zenodo/201894_001_Adcyap1_PSTN.swc | [7194.8, 3568.6, 3567.5] | [12832.2, 6418.1, 5639.6] | 25.7 | 0.0039 to 55.8 | CRLF (35511) | 0 |
| zenodo/201922_073_Nts_MPO.swc | [4429.0, 2928.3, 4218.2] | [10856.9, 7117.3, 6314.3] | 27.3 | 0.0039 to 67.3 | CRLF (52470) | 0 |
| zenodo/202559_043_Tac2_HY.swc | [5163.4, 5560.3, 4181.0] | [7849.1, 6958.3, 6352.6] | 10.9 | 0.0039 to 25.7 | CRLF (12515) | 0 |
| zenodo/202617_025_Sst_ZI.swc | [5902.3, 4939.2, 4282.3] | [7285.7, 5757.1, 5637.9] | 19.2 | 0.0039 to 59.0 | CRLF (11248) | 0 |
| zenodo/202618_143_Sst_LHA.swc | [6148.4, 3432.2, 3854.7] | [8228.3, 6001.7, 6123.7] | 16.8 | 0.0039 to 37.6 | CRLF (12574) | 0 |
| zenodo/202620_016_Sst_PeF.swc | [4038.9, 3029.0, 4151.7] | [7857.3, 6417.6, 6133.8] | 40.2 | 0.0039 to 60.3 | CRLF (35079) | 0 |
| zenodo/210796_002_Pdyn_VMH.swc | [5558.9, 1599.6, 4481.2] | [10172.4, 7053.2, 6325.4] | 138.0 | 0.0039 to 55.9 | CRLF (14681) | 0 |
| zenodo/211266_057_Pdyn_HY.swc | [6522.2, 3323.8, 4897.3] | [9724.0, 6064.3, 7773.3] | 40.5 | 0.0039 to 39.1 | CRLF (22255) | 0 |
| zenodo/211353_070_Adcyap1_AHN.swc | [4661.0, 4687.4, 2814.1] | [7985.2, 7210.2, 5989.6] | 29.4 | 0.0039 to 64.2 | CRLF (25519) | 0 |
| zenodo/211461_026_Nts_PH.swc | [6566.9, 2277.2, 3717.6] | [11875.1, 7309.4, 5710.8] | 30.3 | 0.0039 to 17.3 | CRLF (36243) | 0 |
| zenodo/211465_137_Tac2_ZI.swc | [6491.5, 1353.7, 3739.5] | [10150.9, 5511.2, 5921.6] | 18.9 | 0.0039 to 40.1 | CRLF (21276) | 0 |
| zenodo/211595_037_Sst_HY.swc | [6781.8, 3233.9, 3770.8] | [9994.1, 5775.7, 5553.6] | 53.2 | 0.0039 to 61.0 | CRLF (12304) | 0 |
| zenodo/211595_063_Sst_HY.swc | [5455.8, 2911.3, 4579.6] | [7900.4, 7136.2, 6234.7] | 34.3 | 0.0039 to 68.2 | CRLF (14414) | 0 |
| zenodo/211889_096_Pdyn_HY.swc | [5698.7, 3190.9, 4064.3] | [7448.9, 6189.2, 6254.7] | 27.9 | 0.0039 to 37.3 | CRLF (9701) | 0 |
| zenodo/211890_005_Pdyn_VMH.swc | [5264.6, 1649.0, 3133.2] | [10427.8, 7153.2, 7114.8] | 29.2 | 0.0039 to 70.1 | CRLF (40816) | 0 |
| portal/192710_002.swc | [5737.9, 5146.4, 2874.9] | [7637.9, 7116.5, 5628.7] | 11.4 | 0.0039 to 5.4 | LF (19306) | 0 |
| portal/195038_041.swc | [4967.7, 2442.6, 3813.9] | [10556.2, 7014.1, 5536.2] | 22.9 | 0.0039 to 55.9 | LF (9576) | 0 |
| portal/201232_056.swc | [5990.0, 4168.1, 5693.6] | [10091.8, 6255.5, 6635.4] | 17.5 | 0.0039 to 50.6 | LF (6875) | 0 |
| portal/202457_066.swc | [5463.1, 3401.6, 2894.9] | [7567.7, 7282.0, 6282.1] | 12.5 | 0.0039 to 2.7 | LF (12967) | 0 |
| portal/202813_035.swc | [4246.8, 3075.2, 4424.9] | [10132.0, 7245.2, 5706.5] | 16.6 | 0.0039 to 37.9 | LF (11814) | 0 |
| portal/202813_009.swc | [5013.0, 5350.9, 5728.2] | [7803.7, 7158.6, 7000.0] | 28.2 | 0.0039 to 21.1 | LF (11384) | 0 |
| portal/202825_072.swc | [6195.3, 3695.8, 3156.7] | [14206.9, 7435.2, 5809.0] | 30.2 | 0.0039 to 62.6 | LF (8429) | 3 |
| portal/202825_065.swc | [6150.9, 3738.7, 5734.2] | [10468.2, 6514.2, 7801.9] | 42.6 | 0.0039 to 55.7 | LF (8478) | 0 |
| portal/212248_071.swc | [5838.7, 5274.0, 5584.6] | [7669.0, 7324.4, 7206.3] | 22.4 | 0.0039 to 48.0 | LF (3337) | 0 |
| portal/201221_022.swc | [3509.7, 2521.1, 3741.1] | [10929.1, 6956.1, 6521.6] | 19.7 | 0.0039 to 63.6 | LF (27132) | 0 |
| portal/210105_001.swc | [4742.1, 3077.9, 4273.4] | [10972.3, 7228.4, 7014.1] | 77.1 | 0.0039 to 65.8 | LF (20489) | 0 |
| portal/210104_007.swc | [1614.5, 1139.9, 6075.3] | [13390.9, 7485.1, 8556.0] | 34.1 | 0.0039 to 65.2 | LF (26984) | 16 |
| portal/195039_009.swc | [5979.2, 5434.8, 5631.6] | [7676.2, 7277.6, 7155.9] | 8.4 | 0.0039 to 8.2 | LF (2593) | 13 |
| portal/210073_044.swc | [5195.1, 2218.0, 2999.3] | [10444.6, 7206.7, 5989.0] | 23.4 | 0.0039 to 51.9 | LF (33255) | 0 |
| portal/202170_016.swc | [4899.4, 5173.8, 2711.4] | [7793.0, 7239.1, 7512.0] | 25.1 | 0.0039 to 88.9 | LF (23989) | 0 |
| portal/210007_142.swc | [4851.5, 4753.5, 3738.1] | [7999.8, 6999.3, 5874.4] | 37.8 | 0.0039 to 59.9 | LF (13225) | 0 |
| portal/210951_004.swc | [3947.2, 2273.5, 5223.9] | [13398.8, 7194.5, 7421.9] | 44.1 | 0.0039 to 61.9 | LF (16520) | 22 |
| portal/201230_079.swc | [5138.3, 1648.2, 2831.4] | [10341.0, 7092.4, 6392.8] | 20.8 | 0.0039 to 39.7 | LF (38246) | 0 |
| portal/212202_026.swc | [5541.6, 2816.4, 3617.2] | [13439.9, 6572.6, 5773.5] | 25.4 | 0.0039 to 30.3 | LF (22542) | 14 |
| portal/212202_011.swc | [6192.0, 3764.5, 6111.3] | [13538.1, 7370.7, 7657.0] | 26.6 | 0.0039 to 92.5 | LF (9952) | 1 |
| portal/202824_062.swc | [4881.8, 1488.3, 4908.8] | [10512.7, 7198.0, 6947.3] | 24.7 | 0.0039 to 38.7 | LF (28318) | 0 |
| portal/211461_026.swc | [6566.9, 2277.2, 5689.2] | [11875.1, 7309.4, 7682.4] | 30.3 | 0.0039 to 17.3 | LF (36243) | 0 |
| portal/202618_143.swc | [6148.4, 3432.2, 5276.3] | [8228.3, 6001.7, 7545.3] | 16.8 | 0.0039 to 37.6 | LF (12574) | 0 |

## Frame check: neurons with both a portal and a Zenodo file

| neuron | index hemisphere | nodes | ids/types/parents identical | max abs dAP, dDV | max abs dML | max abs (ML_p + ML_z - 11400) | relation | bytes portal / Zenodo | byte difference |
|---|---|---|---|---|---|---|---|---|---|
| 195038_041 | Left | 9576 | True | 0.0 | 0.0 | 3772.26 | identical | 530638 / 540214 | 9576 |
| 201232_056 | Right | 6875 | True | 0.0 | 1870.7000000000007 | 0.0 | reflected ML'=11400-ML | 380151 / 387204 | 7053 |
| 202618_143 | Right | 12574 | True | 0.0 | 3690.5200000000004 | 0.0 | reflected ML'=11400-ML | 702662 / 718251 | 15589 |
| 211461_026 | Right | 36243 | True | 0.0 | 3964.76 | 0.0 | reflected ML'=11400-ML | 2066493 / 2113999 | 47506 |

The byte difference is the Zenodo file's CRLF line endings (one extra byte per line, see the line-ending column above), not a content change; radii are also identical.

## Every anomaly (106)

- `zenodo/195038_041_Oxt_PVH.swc`: type 0 (unlabelled) nodes: 9575 (24057 um); index reconstruction_type Axon_only
- `zenodo/201232_040_Tac1_HY.swc`: undocumented node types {'10': 3} (kept as 'other')
- `zenodo/201232_040_Tac1_HY.swc`: 6 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/201232_056_Tac1_MMm.swc`: undocumented node types {'10': 4} (kept as 'other')
- `zenodo/201232_056_Tac1_MMm.swc`: 6 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/201894_001_Adcyap1_PSTN.swc`: undocumented node types {'10': 5} (kept as 'other')
- `zenodo/201894_001_Adcyap1_PSTN.swc`: 9 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/201922_073_Nts_MPO.swc`: undocumented node types {'10': 20, '11': 4} (kept as 'other')
- `zenodo/201922_073_Nts_MPO.swc`: type 0 (unlabelled) nodes: 50958 (126068 um); index reconstruction_type Axon_and_dendrite
- `zenodo/201922_073_Nts_MPO.swc`: 43 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/202559_043_Tac2_HY.swc`: undocumented node types {'10': 22} (kept as 'other')
- `zenodo/202559_043_Tac2_HY.swc`: 43 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/202617_025_Sst_ZI.swc`: undocumented node types {'10': 1, '11': 1, '12': 3} (kept as 'other')
- `zenodo/202617_025_Sst_ZI.swc`: 8 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/202618_143_Sst_LHA.swc`: undocumented node types {'10': 5, '11': 2} (kept as 'other')
- `zenodo/202618_143_Sst_LHA.swc`: type 0 (unlabelled) nodes: 10959 (29599 um); index reconstruction_type Axon_and_dendrite
- `zenodo/202618_143_Sst_LHA.swc`: 15 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/202620_016_Sst_PeF.swc`: undocumented node types {'10': 2, '11': 1, '12': 3} (kept as 'other')
- `zenodo/202620_016_Sst_PeF.swc`: 9 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/202620_016_Sst_PeF.swc`: soma HY test not unanimous: 14 of 15 rays inside
- `zenodo/210796_002_Pdyn_VMH.swc`: undocumented node types {'10': 2, '11': 1, '12': 4} (kept as 'other')
- `zenodo/210796_002_Pdyn_VMH.swc`: 9 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/210796_002_Pdyn_VMH.swc`: 2 segments longer than 50 um (max 138 um)
- `zenodo/211266_057_Pdyn_HY.swc`: undocumented node types {'12': 2} (kept as 'other')
- `zenodo/211266_057_Pdyn_HY.swc`: type 0 (unlabelled) nodes: 22252 (58329 um); index reconstruction_type Axon_only
- `zenodo/211266_057_Pdyn_HY.swc`: 2 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/211353_070_Adcyap1_AHN.swc`: undocumented node types {'10': 28} (kept as 'other')
- `zenodo/211353_070_Adcyap1_AHN.swc`: 47 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/211461_026_Nts_PH.swc`: undocumented node types {'12': 1} (kept as 'other')
- `zenodo/211461_026_Nts_PH.swc`: type 0 (unlabelled) nodes: 34773 (88264 um); index reconstruction_type Axon_and_dendrite
- `zenodo/211461_026_Nts_PH.swc`: 1 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/211465_137_Tac2_ZI.swc`: undocumented node types {'10': 4} (kept as 'other')
- `zenodo/211465_137_Tac2_ZI.swc`: type 0 (unlabelled) nodes: 20350 (49595 um); index reconstruction_type Axon_and_dendrite
- `zenodo/211465_137_Tac2_ZI.swc`: 8 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/211595_037_Sst_HY.swc`: type 0 (unlabelled) nodes: 11097 (30877 um); index reconstruction_type Axon_and_dendrite
- `zenodo/211595_037_Sst_HY.swc`: 2 segments longer than 50 um (max 53 um)
- `zenodo/211595_063_Sst_HY.swc`: undocumented node types {'10': 1} (kept as 'other')
- `zenodo/211595_063_Sst_HY.swc`: type 0 (unlabelled) nodes: 13421 (37622 um); index reconstruction_type Axon_and_dendrite
- `zenodo/211595_063_Sst_HY.swc`: 1 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/211889_096_Pdyn_HY.swc`: undocumented node types {'12': 2} (kept as 'other')
- `zenodo/211889_096_Pdyn_HY.swc`: 3 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `zenodo/211890_005_Pdyn_VMH.swc`: undocumented node types {'10': 2} (kept as 'other')
- `zenodo/211890_005_Pdyn_VMH.swc`: type 0 (unlabelled) nodes: 39738 (99695 um); index reconstruction_type Axon_and_dendrite
- `zenodo/211890_005_Pdyn_VMH.swc`: 4 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/192710_002.swc`: undocumented node types {'10': 1, '12': 1} (kept as 'other')
- `portal/192710_002.swc`: type 0 (unlabelled) nodes: 18570 (49058 um); index reconstruction_type Axon_and_dendrite
- `portal/192710_002.swc`: 6 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/195038_041.swc`: type 0 (unlabelled) nodes: 9575 (24057 um); index reconstruction_type Axon_only
- `portal/201232_056.swc`: undocumented node types {'10': 4} (kept as 'other')
- `portal/201232_056.swc`: 6 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/202457_066.swc`: undocumented node types {'10': 1} (kept as 'other')
- `portal/202457_066.swc`: type 0 (unlabelled) nodes: 11783 (28770 um); index reconstruction_type Axon_and_dendrite
- `portal/202457_066.swc`: 1 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/202813_035.swc`: undocumented node types {'10': 1} (kept as 'other')
- `portal/202813_035.swc`: type 0 (unlabelled) nodes: 11458 (29004 um); index reconstruction_type Axon_and_dendrite
- `portal/202813_035.swc`: 2 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/202813_009.swc`: undocumented node types {'10': 1} (kept as 'other')
- `portal/202813_009.swc`: type 0 (unlabelled) nodes: 11042 (30735 um); index reconstruction_type Axon_and_dendrite
- `portal/202813_009.swc`: 2 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/202825_072.swc`: type 0 (unlabelled) nodes: 7817 (19060 um); index reconstruction_type Axon_and_dendrite
- `portal/202825_072.swc`: soma outside the PVH mesh for label PVH: 33.3 um from the mesh surface (exact point-to-triangle; nearest vertex 43 um); 0 of 15 rays inside
- `portal/202825_072.swc`: 522 nodes outside the root mesh (3-axis majority vote; max 1073 um from the root surface; types {'0': 522}; AP 13164 to 14207, beyond the caudal end of the CCF root mesh (AP max 13192): descending axon continuing into the spinal cord, outside the atlas volume)
- `portal/202825_065.swc`: type 0 (unlabelled) nodes: 7866 (19762 um); index reconstruction_type Axon_and_dendrite
- `portal/202825_065.swc`: soma outside the PVH mesh for label PVH: 25.3 um from the mesh surface (exact point-to-triangle; nearest vertex 38 um); 0 of 15 rays inside
- `portal/212248_071.swc`: type 0 (unlabelled) nodes: 2730 (6998 um); index reconstruction_type Axon_and_dendrite
- `portal/212248_071.swc`: soma outside the PVH mesh for label PVH: 8.9 um from the mesh surface (exact point-to-triangle; nearest vertex 35 um); 0 of 15 rays inside
- `portal/201221_022.swc`: undocumented node types {'10': 11, '11': 7, '12': 3} (kept as 'other')
- `portal/201221_022.swc`: type 0 (unlabelled) nodes: 26326 (64113 um); index reconstruction_type Axon_and_dendrite
- `portal/201221_022.swc`: 39 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/210105_001.swc`: undocumented node types {'10': 1} (kept as 'other')
- `portal/210105_001.swc`: type 0 (unlabelled) nodes: 19863 (53949 um); index reconstruction_type Axon_and_dendrite
- `portal/210105_001.swc`: 1 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/210105_001.swc`: 3 segments longer than 50 um (max 77 um)
- `portal/210105_001.swc`: xlsx geno Wt(orexin) != portal index type_name Orexin
- `portal/210104_007.swc`: undocumented node types {'10': 2, '11': 1} (kept as 'other')
- `portal/210104_007.swc`: 3 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/210104_007.swc`: 378 nodes outside the root mesh (3-axis majority vote; max 299 um from the root surface; types {'2': 376, '10': 1, '11': 1}; AP 13071 to 13391, beyond the caudal end of the CCF root mesh (AP max 13192): descending axon continuing into the spinal cord, outside the atlas volume)
- `portal/210104_007.swc`: xlsx geno Wt(orexin) != portal index type_name Orexin
- `portal/195039_009.swc`: type 0 (unlabelled) nodes: 1952 (5222 um); index reconstruction_type Axon_and_dendrite
- `portal/195039_009.swc`: 10 nodes outside the root mesh (3-axis majority vote; max 5 um from the root surface; types {'0': 10}; median (AP, DV, ML) [7668.0, 7135.0, 5648.0])
- `portal/210073_044.swc`: type 0 (unlabelled) nodes: 32343 (79133 um); index reconstruction_type Axon_and_dendrite
- `portal/202170_016.swc`: undocumented node types {'10': 8, '11': 3, '12': 1} (kept as 'other')
- `portal/202170_016.swc`: 24 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/210007_142.swc`: undocumented node types {'12': 3} (kept as 'other')
- `portal/210007_142.swc`: 5 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/210951_004.swc`: undocumented node types {'10': 7} (kept as 'other')
- `portal/210951_004.swc`: type 0 (unlabelled) nodes: 13980 (36876 um); index reconstruction_type Axon_and_dendrite
- `portal/210951_004.swc`: 11 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/210951_004.swc`: 728 nodes outside the root mesh (3-axis majority vote; max 207 um from the root surface; types {'0': 728}; AP 13184 to 13399, beyond the caudal end of the CCF root mesh (AP max 13192): descending axon continuing into the spinal cord, outside the atlas volume)
- `portal/201230_079.swc`: undocumented node types {'10': 28} (kept as 'other')
- `portal/201230_079.swc`: 50 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/201230_079.swc`: soma differs from the xlsx soma by 1.92 um (after the frame's expected reflection)
- `portal/212202_026.swc`: undocumented node types {'12': 5} (kept as 'other')
- `portal/212202_026.swc`: 9 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/212202_026.swc`: 364 nodes outside the root mesh (3-axis majority vote; max 248 um from the root surface; types {'2': 364}; AP 13183 to 13440, beyond the caudal end of the CCF root mesh (AP max 13192): descending axon continuing into the spinal cord, outside the atlas volume)
- `portal/212202_011.swc`: undocumented node types {'10': 1} (kept as 'other')
- `portal/212202_011.swc`: 1 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/212202_011.swc`: 201 nodes outside the root mesh (3-axis majority vote; max 358 um from the root surface; types {'2': 201}; AP 13123 to 13538, beyond the caudal end of the CCF root mesh (AP max 13192): descending axon continuing into the spinal cord, outside the atlas volume)
- `portal/202824_062.swc`: undocumented node types {'10': 5, '11': 13, '12': 9} (kept as 'other')
- `portal/202824_062.swc`: 53 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/211461_026.swc`: undocumented node types {'12': 1} (kept as 'other')
- `portal/211461_026.swc`: type 0 (unlabelled) nodes: 34773 (88264 um); index reconstruction_type Axon_and_dendrite
- `portal/211461_026.swc`: 1 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
- `portal/202618_143.swc`: undocumented node types {'10': 5, '11': 2} (kept as 'other')
- `portal/202618_143.swc`: type 0 (unlabelled) nodes: 10959 (29599 um); index reconstruction_type Axon_and_dendrite
- `portal/202618_143.swc`: 15 parent-child edges change type away from the soma; 0 of them between documented types 0/2/3 (the rest enter or leave an undocumented type 10/11/12 node inside an arbor)
