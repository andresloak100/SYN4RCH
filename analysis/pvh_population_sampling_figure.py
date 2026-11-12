"""Export a static sampling-evidence figure from the completed PVH atlas.

Counts describe collected nuclei and author-inferred identities, not projection
prevalence, strength or absent connections. No new fitting or inference occurs.
"""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / 'data/interim/pvh_population_pilot/pvh_population_atlas_v1'
OUTPUT = ROOT / 'reports/figures/pvh_population_sampling_evidence.png'


def main():
    source = ATLAS / 'atlas_counts.csv'
    with source.open(newline='') as stream:
        rows = [r for r in csv.DictReader(stream) if r['grouping'] == 'author_inferred_type']
    labels = sorted({r['group'] for r in rows})
    assert len(labels) == 12, 'Unexpected inferred-type vocabulary'
    index = {label: i for i, label in enumerate(labels)}
    matrix = np.zeros((len(labels), 2), int)
    seen = set()
    for row in rows:
        key = (row['group'], row['positive_collection_target'])
        assert key not in seen, 'Duplicate type/collection key'
        seen.add(key)
        matrix[index[key[0]], ('PB', 'SC').index(key[1])] = int(row['n_nuclei'])
    assert matrix.sum(axis=0).tolist() == [598, 627]

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
    fig = plt.figure(figsize=(8.8, 8.3), facecolor='white')
    ax = fig.add_axes([.31, .23, .43, .58])
    scale = ax.imshow(matrix, cmap='Blues', vmin=0, vmax=int(matrix.max()), aspect='auto')
    ax.set_xticks([0, 1], ['PB\n598 collected', 'Spinal cord\n627 collected'])
    ax.xaxis.tick_top()
    ax.tick_params(axis='both', which='both', length=0, pad=9)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_ylabel('Author-inferred molecular identity', labelpad=13)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for row in range(len(labels)):
        for col in range(2):
            value = int(matrix[row, col])
            ax.text(col, row, str(value) if value else '—', ha='center', va='center', fontsize=12,
                    color='white' if value > matrix.max() * .57 else '#172532')
    ax.set_xticks(np.arange(-.5, 2), minor=True)
    ax.set_yticks(np.arange(-.5, len(labels)), minor=True)
    ax.grid(which='minor', color='white', linewidth=2)
    cax = fig.add_axes([.79, .25, .022, .54])
    bar = fig.colorbar(scale, cax=cax)
    bar.set_label('Sampled nuclei (count)', labelpad=10)
    bar.outline.set_visible(False)
    fig.text(.08, .94, 'PVH population sampling evidence', fontsize=19, weight='bold', color='#172532')
    fig.text(.08, .90, 'Positive retrograde collection groups · existing published data', fontsize=11, color='#465666')
    fig.text(.08, .13, '1,225 original-study nuclei; molecular identities are inferred by the study authors.',
             fontsize=10, color='#273745')
    fig.text(.08, .099, 'A dash means none were sampled. It does not establish an absent projection.',
             fontsize=10, color='#273745')
    fig.text(.08, .068, 'Counts are not connection strength or prevalence. No individual spatial positions are shown.',
             fontsize=9.5, color='#273745')
    fig.text(.08, .032, 'Source: Li et al., Cell Reports (2026), DOI 10.1016/j.celrep.2025.116904',
             fontsize=9, color='#465666')
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=160, facecolor='white')
    plt.close(fig)
    receipt = {'source': str(source.relative_to(ROOT)),
               'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
               'output': str(OUTPUT.relative_to(ROOT)),
               'output_sha256': hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
               'rows': labels, 'columns': ['PB collection', 'SC collection'], 'sampled_counts': matrix.tolist(),
               'scope': 'Sampling counts by AUTHOR_INFERRED identity, not connectivity strength or prevalence.'}
    path = ROOT / 'data/interim/pvh_population_pilot/sampling_figure_receipt.json'
    path.write_text(json.dumps(receipt, indent=2) + '\n')
    print(OUTPUT)


if __name__ == '__main__':
    main()
