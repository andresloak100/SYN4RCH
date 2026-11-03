"""Read the two evidence layers the viewer shows, without embellishing them.

* Literature coverage: papers whose title/abstract/MeSH mention a population
  (harvest screen). A count of mentions, NOT evidence of an effect.
* Causal perturbation records: rows of data/curated/records.jsonl.
"""

from __future__ import annotations

import datetime as _dt
import json
import os
import re
from typing import Dict, Optional


def _utc_mtime(path: str) -> str:
    ts = os.path.getmtime(path)
    return _dt.datetime.fromtimestamp(ts, _dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def literature_from_candidates(path: str) -> Dict:
    counts: Dict[str, Dict[str, int]] = {}
    total = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            total += 1
            s = rec.get("screen") or {}
            for pid in s.get("matched_populations") or []:
                c = counts.setdefault(pid, {"matched_in_text": 0, "include": 0})
                c["matched_in_text"] += 1
                if s.get("decision") == "include":
                    c["include"] += 1
    return {
        "source": os.path.relpath(path),
        "source_kind": "candidates.jsonl",
        "run_timestamp": _utc_mtime(path),
        "run_timestamp_kind": "file modified (UTC)",
        "total_papers": total,
        "counts": counts,
    }


_ROW = re.compile(r"^\|\s*`([a-z0-9_]+)`\s*\|\s*(\w+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*$")


def literature_from_summary(path: str) -> Dict:
    """Parse the per-population table of reports/harvest_summary.md."""
    counts: Dict[str, Dict[str, int]] = {}
    finished: Optional[str] = None
    started: Optional[str] = None
    total: Optional[int] = None
    in_table = False
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            m = re.match(r"^\|\s*finished \(UTC\)\s*\|\s*(\S+)\s*\|", line)
            if m:
                finished = m.group(1)
            m = re.match(r"^\|\s*started \(UTC\)\s*\|\s*(\S+)\s*\|", line)
            if m:
                started = m.group(1)
            m = re.match(r"^\|\s*unique papers after dedupe\s*\|\s*(\d+)\s*\|", line)
            if m:
                total = int(m.group(1))
            if line.startswith("| population_id |"):
                in_table = True
                continue
            if in_table:
                m = _ROW.match(line)
                if m:
                    counts[m.group(1)] = {"matched_in_text": int(m.group(3)),
                                          "include": int(m.group(4))}
                elif not line.startswith("|---") and not line.startswith("|-"):
                    if counts:
                        in_table = False
    if not counts:
        raise ValueError("no per-population table found in %s" % path)
    return {
        "source": os.path.relpath(path),
        "source_kind": "harvest_summary.md",
        "run_timestamp": finished,
        "run_started": started,
        "run_timestamp_kind": "harvest run finished (UTC)",
        "total_papers": total,
        "counts": counts,
    }


def literature(candidates_path: str, summary_path: str) -> Dict:
    if os.path.exists(candidates_path):
        return literature_from_candidates(candidates_path)
    lit = literature_from_summary(summary_path)
    lit["note"] = ("%s not present at build time; counts read from the "
                   "summary report written by the same harvest run" % os.path.relpath(candidates_path))
    return lit


def causal_records(path: str) -> Dict:
    if not os.path.exists(path):
        return {"source": os.path.relpath(path), "exists": False, "n_records": 0,
                "by_population": {},
                "status": "curation in progress; no curated records exist yet"}
    by_pop: Dict[str, int] = {}
    n = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            n += 1
            pid = rec.get("population_id") or (rec.get("population") or {}).get("population_id")
            if pid:
                by_pop[pid] = by_pop.get(pid, 0) + 1
    return {"source": os.path.relpath(path), "exists": True, "n_records": n,
            "by_population": by_pop, "built_from_mtime": _utc_mtime(path),
            "status": "%d curated record(s)" % n}
