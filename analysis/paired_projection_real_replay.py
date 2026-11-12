"""Independently replay saved BARseq2 scores; never fit or rerun a predictor."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data/interim/paired_projection_real_benchmark'
RUN = BASE / 'barseq2_a1_598_logo_v1'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def replay():
    manifest_path = BASE / 'barseq2_a1_inputs/manifest.json'
    manifest = json.loads(manifest_path.read_text())
    complete = json.loads((RUN / 'RUN_COMPLETE.json').read_text())
    summary = json.loads((RUN / 'result.json').read_text())
    for key, path in [('manifest', manifest_path), ('result', RUN / 'result.json'),
                      ('predictions', RUN / 'predictions.csv'), ('started', RUN / 'STARTED.json')]:
        assert sha(path) == complete[key + '_sha256'], key
    admission_path = Path(manifest['admission']['record_path'])
    assert sha(admission_path) == complete['admission_record_sha256'] == manifest['admission']['record_sha256']
    admission = json.loads(admission_path.read_text())
    assert manifest['run_id'] == complete['run_id'] == admission['run_id'] == RUN.name
    assert complete['runner_sha256'] == sha(ROOT / 'analysis/paired_projection_real_benchmark.py')
    assert complete['reference_sha256'] == sha(ROOT / 'analysis/paired_projection_reference.py')
    assert complete['protocol_sha256'] == sha(ROOT / 'docs/paired_projection_real_protocol.json')
    tables = {}
    for key, pin in manifest['files'].items():
        p = manifest_path.parent / pin['path']
        assert sha(p) == pin['sha256'] == complete['file_sha256'][key] == admission['file_sha256'][key]
        tables[key] = pd.read_csv(p)
    cells = tables['cells'].sort_values('cell_id').reset_index(drop=True)
    truth = tables['labels'][['cell_id', 'target', 'captured']].rename(columns={'target': 'endpoint'})
    pred = pd.read_csv(RUN / 'predictions.csv')
    assert not pred.duplicated(['cell_id', 'endpoint', 'method']).any()
    joined = pred.merge(truth, on=['cell_id', 'endpoint'], how='left', validate='many_to_one', suffixes=('', '_source'))
    assert len(joined) == len(pred) and (joined.captured == joined.captured_source).all()
    cg = dict(zip(cells.cell_id, cells.group_id))
    assert all(cg[c] == a for c, a in zip(pred.cell_id, pred.group_id))
    group_fold = {g: i for i, g in enumerate(sorted(cells.group_id.unique()))}
    assert all(group_fold[g] == f for g, f in zip(pred.group_id, pred.fold))
    methods, targets, animals = summary['methods'], manifest['targets'], sorted(cells.group_id.unique())
    assert set(pred.method) == set(methods) and set(pred.endpoint) == set(targets)
    assert len(pred) == len(cells) * len(targets) * len(methods)
    assert np.isfinite(pred.probability).all() and pred.probability.between(0, 1).all()
    errors, scores, by_animal = [], {}, {}
    for method in methods:
        ep_scores = []
        by_animal[method] = {a: [] for a in animals}
        for target in targets:
            animal_scores = []
            for animal in animals:
                sub = pred[(pred.method == method) & (pred.endpoint == target) & (pred.group_id == animal)]
                assert set(sub.cell_id) == set(cells.loc[cells.group_id == animal, 'cell_id'])
                p, y = sub.probability.to_numpy(), sub.captured.to_numpy()
                pc = np.clip(p, 1e-6, 1 - 1e-6)
                s = np.array([np.mean((p-y)**2), -np.mean(y*np.log(pc)+(1-y)*np.log(1-pc))])
                animal_scores.append(s)
                by_animal[method][animal].append(s)
                reported = summary['per_endpoint'][target]['methods'][method]['per_animal'][animal]
                errors.extend(abs(s - [reported['brier'], reported['log_loss']]))
            em = np.mean(animal_scores, axis=0)
            ep_scores.append(em)
            reported = summary['per_endpoint'][target]['methods'][method]
            errors.extend(abs(em - [reported['brier_animal_equal'], reported['log_loss_animal_equal']]))
        macro = np.mean(ep_scores, axis=0)
        reported = summary['macro']['methods'][method]
        errors.extend(abs(macro - [reported['brier_target_equal'], reported['log_loss_target_equal']]))
        scores[method] = dict(brier=float(macro[0]), log_loss=float(macro[1]))
        for animal in animals:
            s = np.mean(by_animal[method][animal], axis=0)
            r = reported['per_animal_macro'][animal]
            errors.extend(abs(s - [r['brier'], r['log_loss']]))
    # Reconstruct only deterministic copying controls from permitted training inputs.
    # This is an audit of saved predictions, not another predictor fit or a new test.
    genes = np.log1p(tables['genes'].set_index('cell_id').loc[cells.cell_id, manifest['endogenous_genes']].to_numpy())
    depth = cells[manifest['position_columns']].to_numpy()
    labels = truth.pivot(index='cell_id', columns='endpoint', values='captured').loc[cells.cell_id, targets].to_numpy()
    prediction_lookup = pred.set_index(['cell_id', 'endpoint', 'method']).probability
    copy_checks = 0
    for heldout in animals:
        tr, ev = cells.group_id.to_numpy() != heldout, cells.group_id.to_numpy() == heldout
        for method, inputs in [('copy_1nn_location', depth), ('copy_1nn_matched', np.hstack([genes, depth]))]:
            mu, sd = inputs[tr].mean(axis=0), inputs[tr].std(axis=0)
            sd[sd == 0] = 1
            a, b = (inputs[tr]-mu)/sd, (inputs[ev]-mu)/sd
            dist = np.sum((b[:, None, :] - a[None, :, :])**2, axis=2)
            for qi, cid in enumerate(cells.cell_id[ev]):
                tied = np.flatnonzero(np.isclose(dist[qi], dist[qi].min(), rtol=1e-12, atol=1e-12))
                for ti, target in enumerate(targets):
                    p = prediction_lookup.loc[(cid, target, method)]
                    assert p in labels[tr][tied, ti], (cid, target, method, 'not a nearest training donor label')
                    copy_checks += 1
    max_error = float(max(errors))
    assert max_error < 1e-12, max_error
    record = {'schema': 'synarch_barseq2_saved_prediction_replay_v1', 'decision': 'PASS',
              'replay_code_sha256': sha(__file__), 'result_sha256': sha(RUN / 'result.json'),
              'predictions_sha256': sha(RUN / 'predictions.csv'), 'complete_sha256': sha(RUN / 'RUN_COMPLETE.json'),
              'prediction_rows': len(pred), 'cells': len(cells), 'targets': len(targets), 'animals': animals,
              'same_query_ids_labels_groups_folds_verified': True, 'copy_control_predictions_checked': copy_checks,
              'max_absolute_score_difference': max_error, 'macro_replayed': scores,
              'predictor_fits_performed_by_this_check': 0,
              'scope': 'Replay of the one saved exploratory run; no refitting, new model, selection, p-value or promotion.'}
    p = RUN / 'ROOT_REPLAY.json'
    with p.open('x') as f:
        json.dump(record, f, indent=2, sort_keys=True)
    print(json.dumps({'decision': 'PASS', 'prediction_rows': len(pred), 'copy_checks': copy_checks,
                      'max_abs_error': max_error, 'receipt': str(p), 'sha256': sha(p), 'macro': scores}))
    return record


if __name__ == '__main__':
    replay()
