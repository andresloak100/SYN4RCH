"""Synthetic receipt/coverage/refusal checks; never open the real bundle."""
from pathlib import Path
import importlib.util
import sys

import numpy as np
import pytest

P = Path(__file__).resolve().parents[1] / 'analysis/referee_cache_post_audit_verify.py'
spec = importlib.util.spec_from_file_location('post_audit_verifier_tested', P)
V = importlib.util.module_from_spec(spec); sys.modules[spec.name] = V; spec.loader.exec_module(V)


def receipt(n=4):
    body = {'schema': 'referee_full_distance_audit_v1', 'status': 'complete',
            'scope': 'ALL_upper_triangle_D_bitwise;F_hash_and_structure_only',
            'pairs_verified': n * (n - 1) // 2, 'pairs_expected': n * (n - 1) // 2,
            'row_pairs_verified': [n - i - 1 for i in range(n)], 'bitwise_mismatches': 0,
            'binding': {'synthetic': True}}
    return {**body, 'receipt_sha256': V.digest(body)}


def test_exact_complete_receipt_and_mask_pass():
    r = receipt(); assert V.check_receipt(r, 4, 6) == {'synthetic': True}
    mask = np.ones((4, 4), np.uint8); np.fill_diagonal(mask, 0)
    assert V.check_coverage(mask, 4) == [3, 2, 1, 0]


@pytest.mark.parametrize('field,value', [('pairs_verified', 5), ('row_pairs_verified', [2, 2, 2, 0]),
                                       ('bitwise_mismatches', 1), ('status', 'paused'),
                                       ('scope', 'sample_only')])
def test_semantically_wrong_but_rehashed_receipt_rejected(field, value):
    r = receipt(); r[field] = value; r['receipt_sha256'] = V.digest({k: v for k, v in r.items() if k != 'receipt_sha256'})
    with pytest.raises(ValueError): V.check_receipt(r, 4, 6)


def test_receipt_content_tampering_rejected():
    r = receipt(); r['binding'] = {'tampered': True}
    with pytest.raises(ValueError, match='checksum'): V.check_receipt(r, 4, 6)


@pytest.mark.parametrize('mode', ['missing_pair', 'diagonal', 'invalid_value', 'wrong_dtype'])
def test_corrupt_or_partial_mask_never_looks_complete(mode):
    a = np.ones((4, 4), np.uint8); np.fill_diagonal(a, 0)
    if mode == 'missing_pair': a[0, 3] = a[3, 0] = 0
    if mode == 'diagonal': a[1, 1] = 1
    if mode == 'invalid_value': a[0, 1] = a[1, 0] = 2
    if mode == 'wrong_dtype': a = a.astype(np.float64)
    with pytest.raises(ValueError): V.check_coverage(a, 4)


def test_changed_pre_pin_refused_before_any_bundle_access(tmp_path):
    p = tmp_path / 'pin.json'; p.write_text('{}')
    with pytest.raises(ValueError, match='pre-audit pin changed'): V.verify(p)


def test_missing_receipt_refused_before_training_files_are_opened(tmp_path):
    base = tmp_path / 'bundle'; base.mkdir()
    declaration = base / 'draft_declaration.json'
    declaration.write_text(V.canonical({'training_distance_cache': {'audit_path': 'cache/full_audit.json'}}))
    p = tmp_path / 'pin.json'
    p.write_text(V.canonical({'status': 'ROOT_VERIFIED_FRESH_PRE_AUDIT_PIN', 'bundle': str(base),
                             'training_neurons': 2, 'expected_upper_triangle_pairs': 1,
                             'files_sha256': {'draft_declaration.json': V.sha(declaration)}}))
    with pytest.raises(ValueError, match='receipt is absent'): V.verify(p, V.sha(p))
    assert not (base / 'training').exists()


@pytest.mark.parametrize('path', ['../outside', '/absolute'])
def test_artifact_escape_rejected(tmp_path, path):
    with pytest.raises(ValueError): V.under(tmp_path, path)


def test_tiny_synthetic_receipt_reaches_actual_loader_without_heldout_files(tmp_path):
    # Reuse a synthetic three-donor cache builder, never the real bundle.
    p = Path(__file__).with_name('test_prediction_distance_cache.py')
    s = importlib.util.spec_from_file_location('post_pin_synthetic_cache_fixture', p)
    T = importlib.util.module_from_spec(s); sys.modules[s.name] = T; s.loader.exec_module(T)
    base = tmp_path / 'bundle'; base.mkdir()
    descriptor, donors, _ = T.fixture(base)
    (base / 'train').rename(base / 'training')
    for donor in donors: donor['arbor_path'] = donor['arbor_path'].replace('train/', 'training/')
    rows = [{k: v for k, v in d.items() if k != 'arbor'} for d in donors]
    declaration = {'training_distance_cache': descriptor,
                   'method_artifacts': {p: 'pending_audit' for p in T.C.descriptor_paths(descriptor)},
                   'splits': {'train': rows, 'validation': [{'animal': 'VAL', 'arbor_path': 'NEVER_OPEN'}],
                              'test': [{'animal': 'TEST', 'arbor_path': 'NEVER_OPEN'}]}}
    (base / 'draft_declaration.json').write_text(V.canonical(declaration))
    artifacts = [descriptor[k] for k in ('index_path', 'complete_path', 'distance_path', 'f1_path', 'mask_path')] + ['draft_declaration.json']
    source_paths = ('analysis/prediction_distance_cache.py', 'analysis/brm_geometry_audit.py',
                    'analysis/d110_geometry.py', 'analysis/brm_template.py', 'analysis/prediction_independent_evaluator.py')
    arbors = {r['arbor_path']: r['arbor_sha256'] for r in rows}
    listing = ''.join(f'{h}  {p}\n' for p, h in sorted(arbors.items()))
    pin = {'status': 'ROOT_VERIFIED_FRESH_PRE_AUDIT_PIN', 'bundle': str(base), 'training_neurons': 3,
           'expected_upper_triangle_pairs': 3, 'brain_sample_counts': {'train': 2, 'validation': 1, 'test': 1},
           'files_sha256': {p: V.sha(base / p) for p in artifacts},
           'source_sha256': {p: V.sha(V.ROOT / p) for p in source_paths}, 'training_arbors': arbors,
           'training_arbor_listing_sha256': V.hashlib.sha256(listing.encode()).hexdigest()}
    path = tmp_path / 'pre_pin.json'; path.write_text(V.canonical(pin)); expected_sha = V.sha(path)
    T.C.audit_cache(base, descriptor, donors, workers=1)
    result = V.verify(path, expected_sha)
    assert result['status'] == 'PASS' and result['actual_frozen_loader_accepted']
    assert result['original_training_npz_rehashed'] == 3 and result['pairs_verified'] == 3
    assert not (base / 'NEVER_OPEN').exists()
