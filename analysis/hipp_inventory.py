"""Structural inventory of the NeuroXiv ION hippocampus SWC zip (CC BY-NC; lab use only, no redistribution).

Reads straight from the zip; never extracts. Per file: field/id/root checks, coordinate magnitude, 50 um subdivision
budget, axon cable and extent, soma, joined to NeuroXiv metadata. Also checks animal-ID overlap with the earlier
hypothalamic sources so the "never seen" claim is literal.

  /opt/anaconda3/bin/python -I -u analysis/hipp_inventory.py
"""
from __future__ import annotations

import csv
import io
import json
import zipfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ZIP = ROOT / 'data/raw/hippocampus_neuroxiv/ION_Hipp_full_SWC_CCFv3.zip'
META = ROOT / 'data/raw/axon_pretrain_v1/metadata/neuroxiv_metadata.csv'
OUT = ROOT / 'data/interim/hippocampus'
MAX_ABS = 1e6
MAX_POINTS = 5_000_000
ATLAS = np.array([13200., 8000., 11400.])


def inspect(name):
    rec = {'file': name.split('/')[-1], 'neuron_id': name.split('/')[-1][:-4]}
    rec['animal'] = rec['neuron_id'].split('_')[2]
    fails = []
    with zipfile.ZipFile(ZIP) as z:
        text = z.read(name).decode('utf-8', 'replace')
    rows = []
    for ln in text.splitlines():
        s = ln.strip()
        if not s or s.startswith('#'):
            continue
        p = s.split()
        if len(p) < 7:
            fails.append('fields')
            break
        rows.append(p[:7])
    try:
        a = np.asarray(rows, float)
    except ValueError:
        a = None
        fails.append('non_numeric')
    if a is None or len(a) == 0 or 'fields' in fails:
        rec.update(ok=False, failed=';'.join(sorted(set(fails)) or ['empty']))
        return rec
    if not np.isfinite(a).all():
        fails.append('nonfinite')
    ids, typ, xyz, par = a[:, 0].astype(np.int64), a[:, 1].astype(np.int64), a[:, 2:5], a[:, 6].astype(np.int64)
    rec['nodes'] = len(a)
    if len(np.unique(ids)) != len(ids):
        fails.append('dup_ids')
    roots = int((par == -1).sum())
    rec['roots'] = roots
    rec['soma_rows'] = int((typ == 1).sum())
    pos = {int(i): k for k, i in enumerate(ids)}
    pidx = np.array([pos.get(int(q), -2) if q != -1 else -1 for q in par])
    if (pidx == -2).any():
        fails.append('missing_parent')
    if roots != 1:
        fails.append('multi_root' if roots > 1 else 'no_root')
    rec['max_abs_um'] = float(np.abs(xyz).max()) if np.isfinite(xyz).all() else None
    if rec['max_abs_um'] is None or rec['max_abs_um'] > MAX_ABS:
        fails.append('coordinate_magnitude')
    rec['types'] = ','.join(str(t) for t in sorted(set(typ.tolist())))
    child = np.flatnonzero(pidx >= 0)
    if len(child) and 'coordinate_magnitude' not in fails:
        L = np.linalg.norm(xyz[child] - xyz[pidx[child]], axis=1)
        pts = int(np.ceil(L / 50.).sum())
        rec['subdivision_points'] = pts
        if pts > MAX_POINTS:
            fails.append('subdivision_budget')
        ax = typ[child] == 2
        rec['axon_cable_um'] = float(L[ax].sum())
        rec['total_cable_um'] = float(L.sum())
        rec['max_edge_um'] = float(L.max())
        r0 = int(np.flatnonzero(par == -1)[0]) if roots >= 1 else 0
        soma = xyz[r0]
        rec['soma_x'], rec['soma_y'], rec['soma_z'] = map(float, soma)
        rec['soma_type'] = int(typ[r0])
        axn = xyz[typ == 2]
        rec['axon_extent_um'] = float(np.linalg.norm(axn - soma, axis=1).max()) if len(axn) else 0.
        rec['frac_outside_atlas'] = float(np.mean(((xyz < 0) | (xyz > ATLAS)).any(1)))
        if rec['axon_cable_um'] <= 0:
            fails.append('no_axon')
    rec['ok'] = not fails
    rec['failed'] = ';'.join(sorted(set(fails)))
    return rec


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(ZIP) as z:
        names = [i.filename for i in z.infolist() if i.filename.endswith('.swc')]
    print(len(names), 'swc files', flush=True)
    recs = []
    with ProcessPoolExecutor(6) as ex:
        for k, r in enumerate(ex.map(inspect, names, chunksize=16), 1):
            recs.append(r)
            if k % 1000 == 0:
                print(f'  {k}/{len(names)}', flush=True)
    meta = {m['Morphology ID']: m for m in csv.DictReader(open(META))}
    for r in recs:
        m = meta.get(r['neuron_id'])
        r['in_neuroxiv_metadata'] = m is not None
        r['soma_region'] = m['Soma region (based on soma coordinates in ARA brain structures)'] if m else ''
        r['soma_region_source'] = m['Soma region (provide by data source, manually confirmed by human experts)'] if m else ''
        r['hemisphere'] = m['Hemisphere'] if m else ''
    keys = sorted({k for r in recs for k in r})
    with open(OUT / 'inventory.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        w.writerows(sorted(recs, key=lambda r: r['neuron_id']))

    # animal overlap with earlier hypothalamic sources
    prior = set()
    man = ROOT / 'reports/axon_dataset_manifest.tsv'
    if man.exists():
        rows = list(csv.DictReader(open(man), delimiter='\t'))
        col = next((c for c in rows[0] if 'animal' in c.lower() or 'brain' in c.lower() or 'sample' in c.lower()), None)
        prior |= {str(x[col]) for x in rows} if col else set()
    for idx in ('spcd_index.json', 'pvh_oxt_index.json'):
        d = json.loads((ROOT / 'data/raw/axon_pretrain_v1/metadata' / idx).read_text())
        nd = d['neuron_data']
        for n in (nd.values() if isinstance(nd, dict) else nd):
            prior.add(str(n['file']).split('_')[0])
    ok = [r for r in recs if r['ok']]
    animals = Counter(r['animal'] for r in ok)
    overlap = sorted(set(animals) & prior)
    summary = {
        'zip_sha256': (ZIP.parent / (ZIP.name + '.sha256')).read_text().split()[0],
        'files': len(recs), 'structurally_ok': len(ok),
        'failures': Counter(r['failed'] for r in recs if not r['ok']).most_common(),
        'animals_ok': len(animals),
        'neurons_per_animal_p10_p50_p90_max': [int(np.quantile(list(animals.values()), q)) for q in (.1, .5, .9)] + [max(animals.values())],
        'soma_regions_top': Counter(r['soma_region'] for r in ok).most_common(15),
        'axon_cable_um_p10_p50_p90': np.quantile([r['axon_cable_um'] for r in ok], [.1, .5, .9]).tolist(),
        'axon_extent_um_p10_p50_p90': np.quantile([r['axon_extent_um'] for r in ok], [.1, .5, .9]).tolist(),
        'max_abs_um_ok': max(r['max_abs_um'] for r in ok),
        'max_edge_um_ok_p99_max': [float(np.quantile([r['max_edge_um'] for r in ok], .99)), max(r['max_edge_um'] for r in ok)],
        'hemisphere': Counter(r['hemisphere'] for r in ok).most_common(),
        'missing_from_metadata': sum(not r['in_neuroxiv_metadata'] for r in recs),
        'prior_hypothalamic_animal_ids_checked': len(prior),
        'animal_overlap_with_prior_sources': overlap,
    }
    (OUT / 'inventory_summary.json').write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
