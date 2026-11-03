"""Crossref client -- DOI metadata backfill ONLY.

Crossref is never used to discover papers for this corpus (its subject coverage
is not selective enough to screen on). It is used to fill gaps in bibliographic
fields that PubMed or Europe PMC left empty: journal, year, first/last author and
the last author's affiliation, which SCHEMA.md 3.1 requires for lab-level leakage
control.

Backfill never overwrites a value that already exists, and every field it adds is
listed in the record's ``crossref_backfilled`` key so provenance stays visible.
No email address ("mailto") is sent, even though Crossref's polite pool asks for
one; we stay in the anonymous pool instead.
"""

from __future__ import annotations

import json
import urllib.parse
from typing import Dict, List, Optional

from .cache import Cache, RawLog, FetchError, fetch

BASE = "https://api.crossref.org/works"
SOURCE = "crossref"

_BACKFILLABLE = ("journal", "year", "first_author", "last_author",
                 "last_author_affiliation", "title")


def _name(author: Dict) -> Optional[str]:
    given = (author.get("given") or "").strip()
    family = (author.get("family") or "").strip()
    if given and family:
        return "%s %s" % (given, family)
    if family:
        return family
    name = (author.get("name") or "").strip()
    return name or None


def _affiliation(author: Dict) -> Optional[str]:
    affs = author.get("affiliation") or []
    for aff in affs:
        name = (aff.get("name") or "").strip()
        if name:
            return name
    return None


def lookup(doi: str, *, cache: Optional[Cache] = None,
           raw_log: Optional[RawLog] = None) -> Optional[Dict]:
    """Fetch and flatten Crossref metadata for one DOI, or None."""
    if not doi:
        return None
    url = "%s/%s" % (BASE, urllib.parse.quote(doi, safe=""))
    try:
        body, _cached = fetch(url, source=SOURCE, query=doi, cache=cache,
                              raw_log=raw_log, ext="json")
    except FetchError:
        return None
    try:
        msg = json.loads(body).get("message") or {}
    except ValueError:
        return None

    titles: List[str] = msg.get("title") or []
    containers: List[str] = msg.get("container-title") or []
    authors: List[Dict] = msg.get("author") or []

    year = None
    for key in ("published-print", "published-online", "published", "issued", "created"):
        parts = (msg.get(key) or {}).get("date-parts") or []
        if parts and parts[0] and parts[0][0]:
            try:
                year = int(parts[0][0])
            except (TypeError, ValueError):
                year = None
            if year:
                break

    return {
        "doi": (msg.get("DOI") or doi).lower(),
        "title": titles[0] if titles else None,
        "journal": containers[0] if containers else None,
        "year": year,
        "first_author": _name(authors[0]) if authors else None,
        "last_author": _name(authors[-1]) if authors else None,
        "last_author_affiliation": _affiliation(authors[-1]) if authors else None,
        "crossref_type": msg.get("type"),
    }


def backfill(record: Dict, *, cache: Optional[Cache] = None,
             raw_log: Optional[RawLog] = None) -> Dict:
    """Fill only the empty bibliographic fields of ``record`` from Crossref."""
    doi = record.get("doi")
    missing = [f for f in _BACKFILLABLE if not record.get(f)]
    if not doi or not missing:
        return record
    meta = lookup(doi, cache=cache, raw_log=raw_log)
    if not meta:
        return record
    filled = []
    for field in missing:
        val = meta.get(field)
        if val:
            record[field] = val
            filled.append(field)
    if filled:
        record["crossref_backfilled"] = filled
        if record.get("crossref_type"):
            pass
        record.setdefault("crossref_type", meta.get("crossref_type"))
    return record
