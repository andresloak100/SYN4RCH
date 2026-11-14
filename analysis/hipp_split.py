"""Freeze a by-mouse train/val/test split of the ION hippocampus set. Metadata only: animal IDs and a fixed salt.

Rank = sha256(SALT + animal). Lowest 48 mice -> test, next 38 -> val, rest -> train. No geometry, region or
performance quantity is used. Test mice are scored exactly once, after the protocol is frozen.
"""
import csv, hashlib, json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INV = ROOT / 'data/interim/hippocampus/inventory.csv'
OUT = ROOT / 'data/interim/hippocampus/split.json'
SALT = 'synarch-hipp-split-2026-11-14'
N_TEST, N_VAL = 48, 38

rows = list(csv.DictReader(open(INV)))
animals = sorted({r['animal'] for r in rows})
rank = sorted(animals, key=lambda a: hashlib.sha256((SALT + a).encode()).hexdigest())
part = {a: 'test' for a in rank[:N_TEST]}
part.update({a: 'val' for a in rank[N_TEST:N_TEST + N_VAL]})
part.update({a: 'train' for a in rank[N_TEST + N_VAL:]})
n = Counter(part[r['animal']] for r in rows)
split = {'schema': 'synarch_hipp_split_v1', 'salt': SALT, 'rule': __doc__.strip(),
         'created_utc': datetime.now(timezone.utc).isoformat(),
         'inventory_sha256': hashlib.sha256(INV.read_bytes()).hexdigest(),
         'animals': {k: sorted(a for a in part if part[a] == k) for k in ('train', 'val', 'test')},
         'neurons': dict(n)}
with OUT.open('x') as f:
    json.dump(split, f, indent=1)
print({k: len(v) for k, v in split['animals'].items()}, dict(n))
print('split sha256', hashlib.sha256(OUT.read_bytes()).hexdigest())
