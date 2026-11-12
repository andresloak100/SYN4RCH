"""Extract the prospectively selected BARseq2 A1 598-cell cohort; never fit a model.

Source: Sun et al. 2021, Mendeley jnx89bmv4s/2, CC BY 4.0.
The animal crosswalk and target coverage were independently checked against the
deposited pooling scripts and sample sheets. Imported MAT files are data only.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/interim/paired_projection_search_round4/barseq2/A1pooleddata.mat'
OUT = ROOT / 'data/interim/paired_projection_real_benchmark/barseq2_a1_inputs'
SOURCE_SHA = 'd4b1ac270b90c15ccb21a5ce06e9baf02eb38eebd52388c1067147be8d8ffa31'
GENES = ['Cdh4', 'Cdh6', 'Cdh8', 'Cdh9', 'Cdh13', 'Cdh18', 'Pcdh20', 'Pcdh17',
         'Pcdh8', 'Pcdh11x', 'Cdh10', 'Cdh11', 'Cdh12', 'Cdh15', 'Cdh20',
         'Cdh22', 'Cdh24', 'Pcdh19', 'Pcdh7', 'Pcdh9', 'Slc17a7', 'Gad1', 'Slc30a3']
TARGETS = ['OFC', 'M', 'Rstr', 'SS', 'Cstr', 'Amyg', 'VisIp', 'VisC', 'AudC', 'Thal', 'Tect']
ANIMALS = {1: 'XC95', 2: 'XC97'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build_tables(source=SOURCE):
    if sha(source) != SOURCE_SHA:
        raise ValueError('source MAT differs from the pinned public deposit')
    d = loadmat(source, simplify_cells=True)
    if list(d['genelabels']) != GENES or list(d['labels']) != ['OB', *TARGETS]:
        raise ValueError('source gene/target order changed')
    if list(d['genelabelsctrl']) != ['Cdh7', 'Cdh2', 'Unused-1', 'Unused-2']:
        raise ValueError('control-column order changed')
    g = np.asarray(d['allbcexpmat'], dtype=float)
    r = np.asarray(d['projraw'], dtype=float)
    dep = np.asarray(d['allbcdepths'], dtype=float)
    ctrl = np.asarray(d['allbcexpmatctrl'], dtype=float)
    brain = np.asarray(d['brainidx']).astype(int)
    local = np.asarray(d['allbcid'], dtype=float)
    if g.shape != (1283, 23) or r.shape != (1283, 12) or dep.shape != (1283, 2) or ctrl.shape != (1283, 4):
        raise ValueError('source shapes changed')
    if not np.array_equal(brain, np.r_[np.ones(702, dtype=int), np.full(581, 2, dtype=int)]):
        raise ValueError('source mouse mapping changed')
    for a in (g, r, dep, ctrl, local):
        if not np.isfinite(a).all() or (a < 0).any():
            raise ValueError('nonfinite or negative source value')
    if (dep.sum(axis=1) <= 0).any():
        raise ValueError('undefined cortical depth')
    selected = (r.max(axis=1) > 2) & (g[:, [GENES.index('Slc17a7'), GENES.index('Gad1')]].max(axis=1) >= 5)
    counts = {ANIMALS[b]: int((selected & (brain == b)).sum()) for b in ANIMALS}
    if counts != {'XC95': 304, 'XC97': 294}:
        raise ValueError(f'frozen 598-cell cohort changed: {counts}')
    rows = np.flatnonzero(selected)
    group = np.array([ANIMALS[b] for b in brain[rows]])
    ids = np.array([f'{ANIMALS[brain[i]]}:local:{local[i]:.17g}' for i in rows])
    if len(set(ids)) != len(ids):
        raise ValueError('nonunique composite cell IDs')
    cells = pd.DataFrame({'cell_id': ids, 'group_id': group,
                          'depth_norm1200_um': 1200 * dep[rows, 0] / dep[rows].sum(axis=1),
                          'log_panel_total': np.log1p(g[rows].sum(axis=1)),
                          'log_unused_spots': np.log1p(ctrl[rows, 2:].sum(axis=1))})
    genes = pd.DataFrame(g[rows], columns=GENES)
    genes.insert(0, 'cell_id', ids)
    labels = pd.DataFrame([{'cell_id': cid, 'target': tg, 'captured': int(r[i, j] > 2), 'qc_pass': 1}
                           for cid, i in zip(ids, rows) for j, tg in enumerate(['OB', *TARGETS]) if j > 0])
    mask = pd.DataFrame([{'group_id': a, 'target': t, 'assayed': 1} for a in ANIMALS.values() for t in TARGETS])
    crosswalk = pd.DataFrame({'cell_id': ids, 'group_id': group, 'source_row_0based': rows,
                              'source_local_id': local[rows], 'source_brainidx': brain[rows]})
    return {'cells': cells, 'genes': genes, 'labels': labels, 'mask': mask}, crosswalk


def extract(out=OUT):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    tables, crosswalk = build_tables()
    files = {}
    for k, frame in {**tables, 'crosswalk': crosswalk}.items():
        p = out / f'{k}.csv'
        frame.to_csv(p, index=False, float_format='%.17g')
        files[k] = {'path': p.name, 'sha256': sha(p), 'rows': len(frame)}
    ids = tables['cells']['cell_id'].tolist()
    cohort = {'n': len(ids), 'sha256': hashlib.sha256('\n'.join(sorted(ids)).encode()).hexdigest()}
    count_rows = []
    for animal in ANIMALS.values():
        cell_ids = set(tables['cells'].loc[tables['cells'].group_id == animal, 'cell_id'])
        for target in TARGETS:
            y = tables['labels'].loc[tables['labels'].cell_id.isin(cell_ids) &
                                     (tables['labels'].target == target), 'captured']
            if not (0 < y.sum() < len(y)):
                raise ValueError(f'{animal}/{target}: one-class endpoint')
            count_rows.append({'animal': animal, 'target': target, 'n': len(y), 'captured': int(y.sum())})
    rec = {'schema': 'synarch_barseq2_a1_extraction_v1', 'source_path': str(SOURCE.relative_to(ROOT)),
           'source_sha256': SOURCE_SHA, 'extractor_sha256': sha(__file__), 'files': files,
           'genes': GENES, 'targets': TARGETS, 'animal_counts': {'XC95': 304, 'XC97': 294},
           'frozen_cohort': {t: cohort for t in TARGETS}, 'target_counts': count_rows,
           'sampling': '598 marker-QC cells conditional on some source projraw > 2; includes five inhibitory cells; not an all-neuron census',
           'capture': 'raw target molecule count > 2; <=2 is nondetection, not biological absence',
           'position': 'one measured laminar dimension normalized to a 1200-um cortical-thickness reference; not absolute XYZ; does not control tangential anatomy',
           'technical': 'log_panel_total is log1p(sum of 23 measured marker counts), a capture-size proxy that also contains biology; log_unused_spots is log1p of two unused-code counts, an assay-background proxy; no target barcode counts used as predictors',
           'models_fitted': False}
    p = out / 'extraction_receipt.json'
    p.write_text(json.dumps(rec, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'receipt': str(p), 'sha256': sha(p), 'cells': len(ids), 'animals': rec['animal_counts']}))
    return rec


if __name__ == '__main__':
    extract()
