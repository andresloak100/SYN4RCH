# analysis/

One-off analyses behind decision-log entries, kept so every number in
`docs/DECISIONS.md` can be rerun from the repository. Not part of the package:
these scripts may use numpy, scipy and anndata (D-001 binds `schema/` and `qc/`).

| script | produces | entry |
|---|---|---|
| `extract_hypomap_means.py` | per-cell-type mean expression from HypoMap | D-091 input |
| `fetch_connectivity.py` | Allen tracing experiments matched to circuits, projection densities | D-093 input |
| `d095_pair_test.py` | the pair-level test, its power (`reports/d095_power.md`), and the D-092 / D-094 correction | D-095 |

Raw inputs live under `data/raw/`, which is not version controlled; each script
says what it reads and caches.
