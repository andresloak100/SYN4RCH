# analysis/

One-off analyses behind decision-log entries, kept so every number in
`docs/DECISIONS.md` can be rerun from the repository. Not part of the package:
these scripts may use numpy, scipy and anndata (D-001 binds `schema/` and `qc/`).

| script | produces | entry |
|---|---|---|
| `extract_hypomap_means.py` | per-cell-type mean expression from HypoMap | D-091 input |
| `fetch_connectivity.py` | Allen tracing experiments matched to circuits, projection densities | D-093 input |
| `d095_pair_test.py` | the pair-level test, its power (`reports/d095_power.md`), and the D-092 / D-094 correction | D-095 |
| `panel_harvest.py` | open-access papers per circuit and readout (`reports/panel_harvest.md`) | D-097 |
| `panel_select.py` | the panel's five readouts and the papers to read (`reports/panel_selection.md`) | D-097 |
| `panel_stage.py` | full text and staged units for the panel's papers (paused, D-098) | D-097 |
| `d099_gene_wiring.py` | traced cell types placed in HypoMap, and the gene-to-wiring gate test (`reports/d099_*.md`) | D-099 |

Raw inputs live under `data/raw/`, which is not version controlled; each script
says what it reads and caches.
