"""Verify completed native TRAIN audit against the original root pre-audit pin.

No independent pair replay, model fitting, scoring or held-out geometry access.
Only invoke after the full receipt exists AND root authorizes this verifier.
Use the same four thread variables (=1) as the audit; runtime is fail-closed.
"""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
PIN = ROOT / 'data/interim/axon_trace_download_audit/referee_cache_pre_audit_pin.json'
PIN_SHA256 = '94bb1746886a0dee9d959e1ba7fa644355882ce445565b321ee8ee6f0a5fd2f5'
THREAD_KEYS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')


def require(condition, message):
    if not condition: raise ValueError(message)


def read(path): return json.loads(Path(path).read_text())
def sha(path):
    with Path(path).open('rb') as handle: return hashlib.file_digest(handle, 'sha256').hexdigest()
def canonical(value): return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n'
def digest(value): return hashlib.sha256(canonical(value).encode()).hexdigest()


def under(base, relative):
    require(isinstance(relative, str) and not Path(relative).is_absolute(), 'Expected relative artifact path')
    target = (Path(base) / relative).resolve()
    require(target.is_relative_to(Path(base).resolve()), 'Artifact escaped pinned bundle')
    return target


def check_receipt(receipt, n, pairs):
    require(receipt.get('schema') == 'referee_full_distance_audit_v1' and receipt.get('status') == 'complete', 'No COMPLETE full audit receipt')
    require(receipt.get('scope') == 'ALL_upper_triangle_D_bitwise;F_hash_and_structure_only', 'Wrong independent audit scope')
    require(receipt.get('pairs_verified') == pairs == receipt.get('pairs_expected'), 'Incomplete pair verification')
    require(receipt.get('row_pairs_verified') == [n - i - 1 for i in range(n)], 'Wrong complete upper-triangle row counts')
    require(receipt.get('bitwise_mismatches') == 0, 'Distance replay mismatch recorded')
    body = {k: v for k, v in receipt.items() if k != 'receipt_sha256'}
    require(receipt.get('receipt_sha256') == digest(body), 'Receipt content checksum mismatch')
    return receipt['binding']


def check_coverage(mask, n):
    import numpy as np
    require(mask.shape == (n, n) and mask.dtype == np.dtype('uint8'), 'Wrong independent coverage shape/type')
    require(not mask.diagonal().any(), 'Independent coverage diagonal is nonzero')
    for start in range(0, n, 128):
        end = min(start + 128, n); block = mask[start:end]
        expected = np.ones(block.shape, np.uint8); expected[np.arange(end - start), np.arange(start, end)] = 0
        require(np.array_equal(block, expected), 'Independent full coverage contains an unverified/invalid entry')
        require(np.array_equal(block, mask[:, start:end].T), 'Independent coverage is asymmetric')
    row_upper = [int(mask[i, i + 1:].sum()) for i in range(n)]
    require(row_upper == [n - i - 1 for i in range(n)], 'Independent coverage row count mismatch')
    return row_upper


def verify(pin_path=PIN, expected_pin_sha=PIN_SHA256):
    require(sha(pin_path) == expected_pin_sha, 'Original pre-audit pin changed')
    pin = read(pin_path); base = Path(pin['bundle']).resolve()
    require(pin['status'] == 'ROOT_VERIFIED_FRESH_PRE_AUDIT_PIN', 'Wrong original pin status')
    n, pairs = pin['training_neurons'], pin['expected_upper_triangle_pairs']
    require(pairs == n * (n - 1) // 2, 'Original pinned pair count is inconsistent')
    declaration_path = under(base, 'draft_declaration.json')
    require(sha(declaration_path) == pin['files_sha256']['draft_declaration.json'], 'Original draft changed')
    declaration = read(declaration_path); descriptor = declaration['training_distance_cache']
    receipt_path = under(base, descriptor['audit_path'])
    require(receipt_path.is_file(), 'Full receipt is absent; do not run acceptance early')
    receipt_file_sha = sha(receipt_path)
    receipt = read(receipt_path); binding = check_receipt(receipt, n, pairs)
    # Every original byte pin is checked before importing a loader or opening a
    # TRAIN tree. Held-out rows below are metadata only, never followed to files.
    artifact_hashes = {p: sha(under(base, p)) for p in pin['files_sha256']}
    source_hashes = {p: sha(under(ROOT, p)) for p in pin['source_sha256']}
    require(artifact_hashes == pin['files_sha256'], 'Pinned original cache/declaration artifact changed')
    require(source_hashes == pin['source_sha256'], 'Pinned independent source changed')
    groups = {k: {r['animal'] for r in rows} for k, rows in declaration['splits'].items()}
    require(set(groups) == {'train', 'validation', 'test'}, 'Exact split names required')
    require(all(not groups[a] & groups[b] for a, b in (('train', 'validation'), ('train', 'test'), ('validation', 'test'))), 'Split animals overlap')
    require({k: len(v) for k, v in groups.items()} == pin['brain_sample_counts'], 'Pinned animal counts changed')
    rows = sorted(declaration['splits']['train'], key=lambda r: r['id'])
    ids = [r['id'] for r in rows]
    require(len(ids) == len(set(ids)) == n, 'TRAIN donor list duplicate/incomplete')
    expected_paths = {r['arbor_path']: r['arbor_sha256'] for r in rows}
    require(len(expected_paths) == n and expected_paths == pin['training_arbors'], 'TRAIN declaration differs from pinned NPZ listing')
    require(all(r['arbor_path'] == f"training/{r['id']}.npz" for r in rows), 'Noncanonical TRAIN geometry path in donor list')
    require(not (base / 'training').is_symlink() and all(not (base / p).is_symlink() for p in expected_paths), 'TRAIN file alias refused before opening geometry')
    actual_listing = {str(p.relative_to(base)) for p in (base / 'training').rglob('*.npz')}
    require(actual_listing == set(expected_paths), 'Extra/missing TRAIN NPZ artifact')
    current_arbors = {p: sha(under(base, p)) for p in expected_paths}
    require(current_arbors == pin['training_arbors'], 'Pinned TRAIN NPZ bytes changed')
    listing = ''.join(f'{h}  {p}\n' for p, h in sorted(current_arbors.items()))
    require(hashlib.sha256(listing.encode()).hexdigest() == pin['training_arbor_listing_sha256'], 'Pinned TRAIN path/checksum listing digest differs')
    index = read(under(base, descriptor['index_path']))
    require(index['ids'] == ids and index['n'] == n and index['pairs'] == pairs, 'Cache is not the exact declared TRAIN pool')
    require(index['ids_sha256'] == hashlib.sha256('\n'.join(ids).encode()).hexdigest(), 'Sorted TRAIN ID digest differs')
    expected_arbors = [{'id': r['id'], 'animal': r['animal'], 'npz_sha256': r['arbor_sha256'], 'arrays_sha256': h}
                       for r, h in zip(rows, index['arbor_sha256'])]
    require(len(index['arbor_sha256']) == n and binding['training_arbors'] == expected_arbors, 'Receipt TRAIN semantic/file identities differ')
    require(binding['n'] == n and binding['pairs'] == pairs and binding['ids_sha256'] == index['ids_sha256']
            and binding['library_fingerprint'] == index['library_fingerprint'], 'Receipt library identity differs')
    keys = ('index_path', 'complete_path', 'distance_path', 'f1_path', 'mask_path')
    require(binding['files_sha256'] == {key: pin['files_sha256'][descriptor[key]] for key in keys}, 'Receipt does not bind original cache bytes')
    require(binding['independent_source_sha256'] == {'helper': source_hashes['analysis/prediction_distance_cache.py'],
                                                   'geometry_math': source_hashes['analysis/brm_geometry_audit.py']}, 'Receipt independent source differs')
    for name in ('analysis/d110_geometry.py', 'analysis/brm_template.py'):
        require(binding['producer_metric'][name] == source_hashes[name], 'Producer metric source differs')
    require(binding['producer_metric'] == index['metric'], 'Receipt/native metric definition differs')
    progress_path = receipt_path.with_name(receipt_path.name + '.progress.json')
    coverage_path = receipt_path.with_name(receipt_path.name + '.verified.npy')
    progress_file_sha, coverage_file_sha = sha(progress_path), sha(coverage_path)
    progress = read(progress_path)
    require(progress['binding'] == binding and progress['coverage_sha256'] == coverage_file_sha, 'Final checkpoint binding/mask changed')
    import numpy as np
    import scipy
    coverage = np.load(coverage_path, mmap_mode='r', allow_pickle=False)
    require(check_coverage(coverage, n) == receipt['row_pairs_verified'], 'Receipt/coverage row mismatch')
    runtime = {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__, 'byteorder': sys.byteorder,
               'platform': platform.platform(), 'machine': platform.machine(), 'numpy_build': getattr(np.__config__, 'CONFIG', {}),
               'thread_environment': {k: os.environ.get(k) for k in THREAD_KEYS}}
    require(runtime == binding['runtime'], 'Current verification/loader runtime differs from audit')
    helper_path = under(ROOT, 'analysis/prediction_distance_cache.py')
    spec = importlib.util.spec_from_file_location('post_audit_original_loader', helper_path)
    helper = importlib.util.module_from_spec(spec); spec.loader.exec_module(helper)
    require(all(path in declaration['method_artifacts'] for path in helper.descriptor_paths(descriptor)), 'Original draft did not declare every cache artifact')
    donors = [{**r, 'arbor': helper.G.load_tree(under(base, r['arbor_path']))} for r in rows]
    # This validates actual TRAIN arrays/file hashes, native matrix values,
    # symmetry and the complete content/runtime-bound receipt. No D replay.
    distance = helper.load_cache(base, descriptor, donors)
    require(distance.shape == (n, n) and distance.dtype == np.dtype('float64') and not distance.flags.writeable, 'Loader did not return immutable native D')
    # Boundaries must remain stable through the final acceptance call.
    require(sha(pin_path) == expected_pin_sha and sha(receipt_path) == receipt_file_sha, 'Pin/receipt changed during acceptance')
    require(sha(progress_path) == progress_file_sha and sha(coverage_path) == coverage_file_sha, 'Final checkpoint changed during acceptance')
    require(all(sha(under(base, p)) == h for p, h in artifact_hashes.items()), 'Original cache artifacts changed during acceptance')
    require(all(sha(under(ROOT, p)) == h for p, h in source_hashes.items()), 'Original source changed during acceptance')
    return {'status': 'PASS', 'scope': 'independent_original_pin_receipt_coverage_and_actual_TRAIN_loader_acceptance',
            'pre_audit_pin_sha256': expected_pin_sha, 'bundle': str(base), 'training_neurons': n,
            'training_animals': len(groups['train']), 'pairs_verified': pairs, 'bitwise_mismatches': 0,
            'original_artifacts_rehashed': len(artifact_hashes), 'original_sources_rehashed': len(source_hashes),
            'original_training_npz_rehashed': len(current_arbors), 'training_listing_sha256_original': pin['training_arbor_listing_sha256'],
            'training_listing_sha256_compact_json': digest(current_arbors), 'receipt_sha256': receipt_file_sha,
            'receipt_content_sha256': receipt['receipt_sha256'], 'coverage_sha256': coverage_file_sha,
            'progress_sha256': progress_file_sha, 'actual_frozen_loader_accepted': True,
            'source_sha256': sha(__file__), 'runtime': runtime, 'utc': datetime.now(timezone.utc).isoformat(),
            'limitations': ['No second all-pair replay: verifies the full original receipt and exact frozen loader acceptance.',
                'D independently audited ALL upper-triangle pairs bitwise; F hash/structure only, not independently metric-audited.',
                'Only TRAIN geometry opened; validation/test rows inspected as split metadata only.',
                'Draft declaration remains a draft; this does not freeze predictions, select models, or authorize truth scoring.']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'Refuse to overwrite an existing verification record')
    result = verify()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as handle: handle.write(canonical(result))
    print(canonical({k: result[k] for k in ('status', 'training_neurons', 'pairs_verified', 'actual_frozen_loader_accepted')}), end='')
