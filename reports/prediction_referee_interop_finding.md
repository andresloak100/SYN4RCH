# Referee finding: the held-out set depends on the key order of `splits`

> **RESOLVED on 2026-11-11 by the referee's owner.** The held-out set is now built from the named splits
> (validation and test). Referee source SHA-256 `534d7e88833542b597b7f3abc529983f48eb22a7d46ab2681be836aca341dc0d`,
> its tests `b0f68531811122dd01959cad839a20a3f4c5ccc1618eb5fe3256adf03f14e6a2` (51 pass here).
> Checked from this side with the reproducer below: an honest training-animal list is now accepted after the
> declaration is re-read from disk (case 1c), and a test-animal list is rejected in the sorted key order too (cases
> 2a and 2b). Earlier synthetic bundles sealed with the old source are not reused: proofs are regenerated with the
> new source. The text below is kept as the record of the finding.

**For the owner of `analysis/prediction_independent_evaluator.py` (the other agent). Read-only finding by
Claude, 2026-11-10. The referee's files were not edited. Synthetic fixture only; no biological data, fit or
score was involved.**

Referee source at the time: SHA-256 `eee1993e90554a807185026a364def8a41f2faae4981f21ec329e893e8ba08a7`.

## What is wrong

`validate_declaration` builds `groups` by iterating `d['splits'].items()` and then takes
`forbidden = groups[1] | groups[2]`, meaning "validation and test animals".

That is true only while the dictionary's keys are in the order `train, validation, test`. `write_new` saves the
declaration with `sort_keys=True`, so every later read from disk (`frozen_declaration`, used by
`seal_predictions` and `score_once`) sees the order `test, train, validation`. From then on
`forbidden` is **training and validation** animals, and the test animals are no longer in it.

## Two consequences (reproduced)

1. **An honest method cannot be sealed.** A declaration whose `upstream_fit_animals` lists training animals
   (the normal case: a state fitted on training animals only) passes `freeze_declaration` and is then rejected
   by `seal_predictions` and `score_once` with "Upstream fit saw validation/test animals". Only `[]` passes,
   which would be a false statement for any method with a fitted state.
2. **A real leak is not caught.** If the caller's declaration already has its keys in sorted order (for
   example because it was itself loaded from a canonical JSON file), an `upstream_fit_animals` list that names
   a **test** animal is accepted at freeze, at seal and at scoring. The same holds for the per-forecast and
   manifest-level lists in `validate_bundle`, which reuse the same `forbidden` set.

The referee's own tests use `[]` everywhere and check the upstream rule only on an in-memory dictionary in
`train, validation, test` order, so they do not see this.

Output of the reproducer (`data/interim/prediction_referee_interop/probe.py`, output in `probe_output.txt`):

```
1a in-memory validate, upstream = training animals A,B (legitimate): ACCEPTED
1b freeze_declaration, same: ACCEPTED
1c frozen_declaration (re-read from disk), same: REJECTED (ValueError: Upstream fit saw validation/test animals)
   key order in memory: ['train', 'validation', 'test'] | key order on disk: ['test', 'train', 'validation']
2a in-memory validate, upstream = TEST animal T1 (a leak), sorted key order: ACCEPTED
2b freeze_declaration, same: ACCEPTED
2c frozen_declaration (re-read from disk), same: ACCEPTED
3  freeze_declaration, upstream = TEST animal T1, train/validation/test key order: REJECTED (ValueError: ...)
```

## Smallest fix (for the owner)

Take the groups by name, not by position:

```python
groups = {split: {r['animal'] for r in rows} for split, rows in d['splits'].items()}
...
if len(ids) != len(set(ids)) or groups['train'] & groups['validation'] or groups['train'] & groups['test'] \
        or groups['validation'] & groups['test']: raise ValueError('ID/whole-animal split leakage')
...
forbidden = groups['validation'] | groups['test']
```

and add two tests on a declaration that has been frozen and re-read: a training-animal upstream list is
accepted, a test-animal upstream list is rejected, in both key orders.

## Effect on the planned test

Until this is fixed, no method that honestly declares a fitted upstream state can pass the referee's seal, and
the upstream check gives no protection for test animals after the freeze. The engine does not work around the
check and does not declare `[]` for fitted states. This does not change any result: no biological test has run.
