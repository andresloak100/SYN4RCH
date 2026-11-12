import csv
import json

import numpy as np
import pytest

from analysis import pvh_population_pilot_replay as R


def test_pool_then_fold_scoring_is_not_pooled_cell_scoring():
    cohort = [dict(object='PB', cell_id='a', pool='p1', target='PB'),
              dict(object='SC', cell_id='b', pool='p2', target='SC'),
              dict(object='SC', cell_id='c', pool='p2', target='SC'),
              dict(object='PB', cell_id='d', pool='p3', target='PB'),
              dict(object='SC', cell_id='e', pool='p4', target='SC')]
    forecasts = {}
    for i, r in enumerate(cohort):
        fold = 'f1' if i < 3 else 'f2'
        for method in R.ALPHAS:
            forecasts[fold, r['object'], r['cell_id'], method] = {'p': [.1, .1, .1, .8, .2][i]}
    folds = {'f1': {'test': ['p1', 'p2']}, 'f2': {'test': ['p3', 'p4']}}
    result = R.scores(forecasts, cohort, folds)['train_prior']
    assert result['brier_macro'] == pytest.approx(((.81+.01)/2+.04)/2)
    assert result['log_loss_macro'] == pytest.approx(((-np.log(.1)-np.log(.9))/2-np.log(.8))/2)


def test_reporter_aliases_do_not_remove_native_genes():
    for name in ['ERCC-0001', 'ERCC_002', 'DsRed', 'EGFP', 'BFP', 'mCherry', 'Cre', 'barcode_1', 'HTO']:
        assert R.excluded_gene(name)
    for name in ['Creb1', 'Creld1', 'Gfpt1', 'Gfpt2', 'Sim1']:
        assert not R.excluded_gene(name)


def test_completion_hashes_are_checked_before_artifact_loading(tmp_path):
    protocol = tmp_path / 'protocol.json'; protocol.write_text('{}')
    folder = tmp_path / 'run'; folder.mkdir()
    (folder / 'protocol_snapshot.json').write_text('{}')
    (folder / 'predictions.csv').write_text('original')
    (folder / 'STARTED.json').write_text(json.dumps({
        'protocol_snapshot_sha256': R.sha(protocol),
        'source_sha256': R.sha(R.ROOT / 'analysis/pvh_population_pilot.py'),
        'extractor_sha256': R.sha(R.ROOT / 'analysis/pvh_population_extract.R')}))
    complete = {p.name: R.sha(p) for p in folder.iterdir()}
    (folder / 'RUN_COMPLETE.json').write_text(json.dumps({'files_sha256': complete}))
    R.verify_run(folder, protocol)
    (folder / 'predictions.csv').write_text('tampered')
    with pytest.raises(ValueError, match='output hash'):
        R.verify_run(folder, protocol)


def test_near_rejects_nonfinite_or_wrong_shape():
    R.near([1., 2.], [1., 2.+1e-14], 'fixture')
    with pytest.raises(ValueError):
        R.near([1., 2.], [np.nan, 2.], 'fixture')
    with pytest.raises(ValueError):
        R.near([1., 2.], [[1., 2.]], 'fixture')


def test_historical_atlas_keeps_its_own_source_and_protocol(tmp_path, monkeypatch):
    import hashlib
    old = {'analysis/pvh_population_pilot.py': b'old primary',
           'analysis/pvh_population_extract.R': b'old extractor',
           'docs/pvh_population_pilot_protocol.json': b'{"status":"old frozen"}'}
    monkeypatch.setattr(R, 'historical_blob', lambda commit, name: old[name])
    folder = tmp_path / 'atlas'; folder.mkdir()
    protocol = tmp_path / 'current.json'; protocol.write_text('{"status":"new frozen"}')
    (folder / 'protocol_snapshot.json').write_bytes(old['docs/pvh_population_pilot_protocol.json'])
    (folder / 'STARTED.json').write_text(json.dumps({
        'source_sha256': hashlib.sha256(old['analysis/pvh_population_pilot.py']).hexdigest(),
        'extractor_sha256': hashlib.sha256(old['analysis/pvh_population_extract.R']).hexdigest(),
        'protocol_snapshot_sha256': R.sha(folder / 'protocol_snapshot.json')}))
    (folder / 'RUN_COMPLETE.json').write_text(json.dumps({'files_sha256': {
        name: R.sha(folder / name) for name in ['STARTED.json', 'protocol_snapshot.json']}}))
    R.verify_run(folder, protocol, historical_commit='old-commit')
    old['analysis/pvh_population_pilot.py'] = b'changed old source'
    with pytest.raises(ValueError, match='historical atlas'):
        R.verify_run(folder, protocol, historical_commit='old-commit')


def test_forecast_key_must_match_its_source_collection_and_pool(tmp_path):
    cohort = [{'object': 'PB', 'cell_id': 'x', 'pool': 'p2', 'sex': 'male', 'target': 'PB'}]
    protocol = {'folds': {'fold_one': {'train': ['p1'], 'test': ['p2']}}}
    path = tmp_path / 'predictions.csv'
    with path.open('w', newline='') as stream:
        columns = ['fold', 'object', 'cell_id', 'pool', 'sex', 'collection_target', 'method', 'probability', 'raw_probability']
        writer = csv.DictWriter(stream, fieldnames=columns); writer.writeheader()
        writer.writerow(dict(fold='fold_one', object='PB', cell_id='x', pool='WRONG', sex='male',
                             collection_target='PB', method='train_prior', probability='.5', raw_probability='.5'))
    with pytest.raises(ValueError, match='forecast identity'):
        R.verify_fit(tmp_path, cohort, np.ones((1, 1)), ['G1'], protocol, {})
