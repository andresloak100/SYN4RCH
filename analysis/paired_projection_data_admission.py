"""Bounded data-admission audit for genuinely paired hypothalamic molecular -> destination data.

Read-only. Fetches only small public metadata/text (each <= 50 MB): GEO series text and file lists,
SRA run info, BioSample records, PMC full-text XML, and the authors' public GitHub tree/notebook.
It never downloads expression matrices or raw reads, never reads assay values, never fits a model
and never relates molecules to projection outcomes. Output: admission.json with source hashes,
extracted facts and an ADMIT / REJECT / PENDING verdict with hard blockers.

Usage: /opt/anaconda3/bin/python -I analysis/paired_projection_data_admission.py run|cached
(cached = no network; re-audit the saved raw files after verifying their fetch hashes)
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/interim/paired_projection_data_admission'
RAW = OUT / 'raw'          # untrusted downloads only; no code lives here
MAX_BYTES = 50 * 1024 * 1024
UA = 'research-metadata-audit'

ROSSI_COMMIT = '41e5b21f7578fc6e9ef6d01c41f6671751a66a1c'
SOURCES = {
    'rossi_geo_soft': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169176&targ=all&form=text&view=brief',
    'rossi_geo_filelist': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169176/suppl/filelist.txt',
    'rossi_sra_runinfo': 'https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=SRP311270',
    'rossi_biosample': ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample'
                        '&id=SAMN18350787,SAMN18350785,SAMN18350784'),
    'rossi_pmc_xml': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=8812999&rettype=xml',
    'rossi_epmc_supp': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8812999/supplementaryFiles',
    'rossi_github_tree': ('https://api.github.com/repos/stuberlab/Rossi-et-al-2021/git/trees/'
                          + ROSSI_COMMIT + '?recursive=1'),
    'rossi_github_repo': 'https://api.github.com/repos/stuberlab/Rossi-et-al-2021',
    'rossi_hcr_notebook': ('https://raw.githubusercontent.com/stuberlab/Rossi-et-al-2021/' + ROSSI_COMMIT
                           + '/Fig3_S2_HCR%20analysis_Final.ipynb'),
    'luo_pmc_xml': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=12581905&rettype=xml',
}


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def fetch(name, url):
    req = urllib.request.Request(url, headers={'User-Agent': UA})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = r.read(MAX_BYTES + 1)
    if len(data) > MAX_BYTES:
        raise ValueError(f'{name}: exceeds 50 MB cap')
    RAW.mkdir(parents=True, exist_ok=True)
    (RAW / name).write_bytes(data)
    return {'name': name, 'url': url, 'bytes': len(data), 'sha256': sha256_bytes(data),
            'fetched_utc': datetime.now(timezone.utc).isoformat()}


def plain_text(xml):
    t = re.sub(r'<[^>]+>', ' ', xml)
    return re.sub(r'\s+', ' ', t).replace('&lt;', '<').replace('&gt;', '>').replace('&amp;', '&')


def quote(text, phrase, after=300):
    i = text.find(phrase)
    return None if i < 0 else text[i:i + after]


def parse_geo_soft(soft):
    samples, cur = [], None
    for line in soft.splitlines():
        if line.startswith('^SAMPLE'):
            cur = {'gsm': line.split('=', 1)[1].strip()}
            samples.append(cur)
        elif cur is not None and '=' in line:
            k, v = (s.strip() for s in line.split('=', 1))
            if k in ('!Sample_title', '!Sample_description'):
                cur[k[8:].lower()] = v
            elif k.startswith('!Sample_supplementary_file'):
                cur.setdefault('files', []).append(v.rsplit('/', 1)[-1])
            elif k == '!Sample_characteristics_ch1':
                cur.setdefault('characteristics', []).append(v)
    series_files = [l.split('=', 1)[1].strip().rsplit('/', 1)[-1]
                    for l in soft.splitlines() if l.startswith('!Series_supplementary_file')]
    return {'samples': samples, 'series_files': series_files}


def parse_filelist(txt):
    rows = []
    for line in txt.splitlines()[1:]:
        p = line.split('\t')
        if len(p) >= 5:
            rows.append({'kind': p[0], 'name': p[1], 'bytes': int(p[3]), 'type': p[4]})
    return rows


def parse_notebook(nb_json):
    """Extract only input file names, metadata derivations and printed counts (no values)."""
    nb = json.loads(nb_json)
    code = '\n'.join(''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'code')
    outs = []
    for c in nb['cells']:
        for o in c.get('outputs', []):
            if 'text' in o:
                outs.append(''.join(o['text']))
    out = '\n'.join(outs)
    inputs = sorted(set(re.findall(r"'([^'\\/]+\.(?:xlsx|csv))'", code)))
    dirs = sorted(set(re.findall(r'"(C:\\\\[^"]+)"', code)))
    total = re.findall(r'Total Cells:\s+(\d+)', out)
    yfp = re.findall(r'Number of YFP cells:\s+(\d+)', out)
    tom = re.findall(r'Number of Tomato cells:\s+(\d+)', out)
    mice = sorted({f.split('_')[0] for f in inputs if re.match(r'AG\d+_', f)})
    return {
        'input_files': inputs, 'input_dirs_private': dirs, 'hcr_mouse_ids_from_filenames': mice,
        'printed_total_cells': [int(x) for x in total], 'printed_yfp_cells': [int(x) for x in yfp],
        'printed_tomato_cells': [int(x) for x in tom],
        'derives_mouse_slide_slice_hemisphere_from_sheet_names': 'MouseSlideSlice' in code,
        'ap_coord_hand_assigned_per_slice': "'AP coord'" in code,
        'xy_from_bounding_box_columns': 'Xcenter' in code,
        'svm_split_random_50pct_not_by_animal': 'test_data=.5' in code,
    }



CHECKED_FOR_HCR = ['GEO GSE169176 series/sample text and suppl filelist', 'SRA SRP311270 run info',
                   'BioSamples SAMN18350784/85/87', 'PMC8812999 full-text XML (incl. Data and Code Availability)',
                   'GitHub stuberlab/Rossi-et-al-2021 tree @' + ROSSI_COMMIT[:8]]
UNRETRIEVED = ['PMC supplement NIHMS1771318-supplement-1.pdf (empty body to a script; Europe PMC: '
               'not open access)', 'publisher (Cell Press) supplementary files (not attempted)']


HCR_CANDIDATE = re.compile(r'^(?:AG\d+_[A-Z]\d+_quant\.xlsx|.*_quant\.xlsx|HCR_data\.csv)$')


def is_hcr_candidate(name):
    """True only for per-cell HCR table names (HALO *_quant.xlsx or HCR_data.csv)."""
    return bool(HCR_CANDIDATE.match(name))


def audit(raw):
    """raw: dict name -> text. Returns the facts and the verdicts."""
    geo = parse_geo_soft(raw['rossi_geo_soft'])
    files = parse_filelist(raw['rossi_geo_filelist'])
    tree = json.loads(raw['rossi_github_tree'])
    repo = json.loads(raw['rossi_github_repo'])
    nb = parse_notebook(raw['rossi_hcr_notebook'])
    rt = plain_text(raw['rossi_pmc_xml'])
    lt = plain_text(raw['luo_pmc_xml'])
    repo_files = [x['path'] for x in tree.get('tree', []) if x['type'] == 'blob']
    data_like = [p for p in repo_files if re.search(r'\.(xlsx|csv|tsv|h5|h5ad|rds|mat|npy|zip)$', p, re.I)]
    # Only per-cell HCR candidates count: HALO '*_quant.xlsx' exports or HCR_data.csv. RNA DEG/Seurat CSVs never
    # count. Presence alone can only move the verdict to PENDING; it never verifies cell pairs.
    hcr_candidates = sorted(f for f in repo_files if is_hcr_candidate(Path(f).name))
    hcr_found = bool(hcr_candidates)
    two_animal_pools = sum('two animals' in s.get('description', '') for s in geo['samples'])
    rossi = {
        'citation': 'Rossi MA et al. Neuron 109:3823-3837.e6 (2021). DOI 10.1016/j.neuron.2021.09.020; PMC8812999',
        'article_licence': 'PMC author manuscript; "available for text mining ... fair use" (not CC)',
        'geo_licence': 'none stated by GEO',
        'github_licence': repo.get('license'),
        'geo_samples': geo['samples'], 'geo_series_files': geo['series_files'], 'geo_filelist': files,
        'github_repo_files': repo_files, 'github_data_like_files': data_like,
        'hcr_notebook': nb,
        'hcr_input_tables_found_in_checked_sources': hcr_found,
        'hcr_candidate_files_found': hcr_candidates,
        'checked_sources_for_hcr_inputs': CHECKED_FOR_HCR, 'unretrieved_sources': UNRETRIEVED,
        'epmc_supp_response': raw['rossi_epmc_supp'][:300],
        'quotes': {k: quote(rt, p, n) for k, (p, n) in {
            'data_availability': ('Data and Code Availability', 420),
            'scrna_animals': ('Single-Cell sequencing C57BL/6J mice', 330),
            'scrna_pooling': ('Following digestion, tissue from all animals', 160),
            'scrna_reactions': ('Approximately 16,000-18,000 pooled cells', 160),
            'scrna_labels': ('We detected eYFP (n=257 neurons)', 160),
            'hcr_animals': ('HCR C57BL/6J mice', 360),
            'hcr_thresholds': ('Puncta filters were as follows', 140),
            'svm_split': ('Training was performed on a random 50%', 90),
            'ephys_counterbalanced': ('viruses counterbalanced between brain regions', 60),
        }.items()},
        'facts': {
            'rna_mice_sequenced': 6, 'rna_reactions': 3,
            'rna_geo_pools_of_two_animals': two_animal_pools,
            'rna_pooling_sources': {
                'GEO': 'each of 3 samples: "pooled live single-cell suspension generated from LHA tissue of '
                       'two animals (male and female)"',
                'paper': '"tissue from all animals was ... pooled" before 3 reactions (ambiguous: all six '
                         'together, or within each two-animal pool)'},
            'rna_pools_independence': 'three two-animal pools; biologically independent only if animal-'
                                      'disjointness is confirmed; per-cell animal IDs are not recoverable '
                                      'from checked sources (never six recovered animal IDs)',
            'rna_labelled_neurons_reported': {'eYFP_LHb': 257, 'tdTomato_VTA': 1595},
            'hcr_mice': 3, 'hcr_genes': 9, 'hcr_marker_genes_shared_with_rna': 7,
            'hcr_per_cell_counts_found_in_checked_sources': hcr_found,
            'fluorophore_target_counterbalanced_rna': False, 'fluorophore_target_counterbalanced_hcr': False,
            'tracer_assignment': 'retroAAV2-hSyn-eYFP -> LHb; retroAAV2-CAG-tdTomato -> VTA (both assays)',
        },
    }
    rossi['hard_blockers'] = [
        'Per-cell HCR counts, cell->mouse/slide/section/hemisphere joins, XY positions and tracer calls were '
        'NOT FOUND in the checked public sources (GEO, SRA, BioSample, PMC text, authors\' GitHub). The '
        'authors\' notebook reads Dropbox HALO exports ' + ', '.join(nb['input_files'][:5]) + '. The PMC '
        'supplement could not be retrieved, so this is not a claim that the data are universally non-public.',
        'Author contact is not permitted in this audit, so the missing items cannot be requested here.',
        'RNA: three two-animal pools (GEO) in three reactions; no per-cell animal ID in checked sources. '
        'Pools may count as independent units only if animal-disjointness is confirmed.',
        'Tracer identity is confounded with destination in both assays (eYFP/hSyn always LHb, '
        'tdTomato/CAG always VTA); only the separate ephys cohort was counterbalanced.',
        'Even if the HCR tables were obtained, 3 HCR mice with 7 shared marker genes are feasibility/'
        'development only, not a powered confirmatory test; an RNA->HCR transform would have to be locked '
        'and calibrated outside the HCR test animals.',
    ]
    rossi['verdict'] = 'REJECT' if not hcr_found else 'PENDING'
    rossi['verdict_scope'] = 'not admitted from the checked public sources; re-screen only if the missing items appear'
    luo = {
        'citation': 'Luo YJ et al. Neuron (2025). DOI 10.1016/j.neuron.2025.09.006; PMC12581905',
        'article_licence': 'PMC: "available for text mining ... fair use" (not CC)',
        'data': 'GSA CRA019056 (no raw download)',
        'quotes': {
            'patchseq_two_cells_per_tube': quote(lt, 'Two cells were collected in one tube as a sample', 60),
            'patchseq_results_22': quote(lt, 'yielding 22 high-quality profiles', 60),
            'snrna_pooling': quote(lt, 'Tissues from 5 mice of the same group were pooled together', 70),
            'snrna_egfp_reference': quote(lt, 'The eGFP gene was added to GRCm39', 140),
            'snrna_egfp_threshold': quote(lt, 'Cells with eGFP count values greater than 0', 110),
            'snrna_egfp_positives': quote(lt, 'Among the four eGFP-labeled cells recovered from SuM-vDG mice', 120),
            'data_availability': quote(lt, 'Single-nucleus RNA-seq and patch-seq data in this study were deposited', 200),
        },
        'facts': {
            'snrna_label_detected_per_nucleus': True,
            'snrna_label_rule': 'eGFP count > 0 in a reference with eGFP added',
            'snrna_vdg_positive_nuclei_reported': 4,
            'snrna_ddg_positive_nuclei_reported': 'not stated in checked text',
            'snrna_label_negative_meaning': 'unknown (not captured), not a verified absence',
            'snrna_mice_per_pool': 5,
            'unresolved': ['pool IDs and pool-to-library mapping', 'number of pools per group (replication)',
                           'the sample unit of the 22 Patch-seq profiles', 'per-nucleus eGFP status for dDG mice'],
            'endpoint_note': 'dDG vs vDG is an endpoint difference from the Epi-retro-seq list, not a hard '
                             'blocker for a new task',
        },
        'hard_blockers': [
            'Patch-seq Methods: "Two cells were collected in one tube as a sample", so the 22 profiles cannot '
            'be taken as individually paired cells.',
            'snRNA-seq: label-positive nuclei can pair RNA with a captured label (eGFP > 0), but only four '
            'positives (all from vDG mice) are reported in the checked text; the remaining nuclei are '
            'unknown, not negatives.',
            'snRNA-seq: 5 mice pooled per group; pool IDs, replication and the sample unit are unresolved in '
            'the checked text, so independent evaluation animals are not established.',
        ],
        'verdict': 'REJECT',
    }
    return {'rossi2021_lha': rossi, 'luo2025_sum': luo}


def write(result, sources, mode):
    path = OUT / 'admission.json'
    result['_meta'] = {
        'script_sha256': sha256_bytes(Path(__file__).read_bytes()),
        'mode': mode, 'finished_utc': datetime.now(timezone.utc).isoformat(),
        'sources': sources,
        'not_downloaded': ['GSE169176_RAW.tar (191,692,800 B) and its three h5 matrices (57-68 MB each): '
                           'over the 50 MB cap and assay values',
                           'SRA reads SRR14001610-12', 'GSA CRA019056 raw data'] + UNRETRIEVED,
        'excluded_by_instruction': ['Kim 2019 VMHvl (closed pilot)', 'Li 2026 PVH (already unqualified)',
                                    'Jiao 2025', 'MouseLight'],
        'overall': 'REJECT (no source admitted)',
    }
    path.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    print('wrote', path, sha256_bytes(path.read_bytes()))
    return result


def load_cached(sources):
    """Re-read cached raw files; refuse if any hash differs from its fetch receipt."""
    raw = {}
    for s in sources:
        b = (RAW / s['name']).read_bytes()
        if sha256_bytes(b) != s['sha256']:
            raise ValueError(f"{s['name']}: cached hash differs from fetch receipt")
        raw[s['name']] = b.decode('utf-8', errors='replace')
    return raw


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    receipts = [fetch(n, u) for n, u in SOURCES.items()]
    return write(audit(load_cached(receipts)), receipts, 'fetch')


def cached():
    """Correction pass: no network; reuse the original fetch receipts and verify every hash."""
    prior = json.loads((OUT / 'admission.json').read_text())
    receipts = prior['_meta']['sources']
    return write(audit(load_cached(receipts)), receipts, 'cached')


if __name__ == '__main__':
    if sys.argv[1:] == ['run']:
        run()
    elif sys.argv[1:] == ['cached']:
        cached()
    else:
        sys.exit(__doc__)
