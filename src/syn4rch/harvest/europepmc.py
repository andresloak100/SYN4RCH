"""Europe PMC REST client: search, open-access detection, OA full-text retrieval.

Europe PMC is used for two things PubMed cannot do:

* it indexes preprints and some non-MEDLINE content, so it widens the search;
* it states open-access status (``isOpenAccess``) and serves OA full text as XML,
  which is what decides ``full_text_available`` in docs/SCHEMA.md 3.1.

Full texts are stored verbatim under ``data/raw/fulltext/<PMCID>.xml``. No email
address is sent.
"""

from __future__ import annotations

import json
import os
import urllib.parse
from typing import Dict, List, Optional, Tuple

from .cache import Cache, RawLog, FetchError, data_path, fetch

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
SOURCE = "europepmc"
POST_THRESHOLD = 1800


def _norm_record(item: Dict) -> Dict:
    doi = item.get("doi")
    pmcid = item.get("pmcid")
    is_oa = (item.get("isOpenAccess") or "").upper() == "Y"
    in_epmc = (item.get("inEPMC") or "").upper() == "Y"
    year = item.get("pubYear")
    try:
        year_i: Optional[int] = int(year) if year else None
    except (TypeError, ValueError):
        year_i = None
    author_string = item.get("authorString") or ""
    authors = [a.strip() for a in author_string.split(",") if a.strip()]
    pub_types = item.get("pubTypeList", {}).get("pubType", []) if isinstance(item.get("pubTypeList"), dict) else []
    if isinstance(pub_types, str):
        pub_types = [pub_types]
    return {
        "pmid": item.get("pmid"),
        "pmcid": pmcid.upper() if pmcid else None,
        "doi": doi.lower() if doi else None,
        "title": item.get("title"),
        "abstract": item.get("abstractText"),
        "journal": (item.get("journalInfo", {}) or {}).get("journal", {}).get("title")
                   or item.get("journalTitle"),
        "year": year_i,
        "first_author": authors[0] if authors else None,
        "last_author": authors[-1] if authors else None,
        "last_author_affiliation": None,  # Europe PMC search does not return affiliations
        "publication_types": [p for p in pub_types if p],
        "mesh_terms": [],
        "n_authors": len(authors),
        "is_open_access": is_oa,
        "in_epmc": in_epmc,
        "source_db": item.get("source"),
        "epmc_id": item.get("id"),
        "retrieved_from": "europepmc_oa" if is_oa else "europepmc_abstract",
    }


def search(query: str, *, page_size: int = 100, cache: Optional[Cache] = None,
           raw_log: Optional[RawLog] = None) -> Tuple[List[Dict], int]:
    """Search Europe PMC. Returns ``(records, hit_count)``.

    ``hit_count`` is Europe PMC's own total for the query; ``records`` is at most
    ``page_size`` long (one page -- no silent deep pagination).
    """
    params = [
        ("query", query),
        ("format", "json"),
        ("pageSize", str(min(page_size, 1000))),
        ("resultType", "core"),
    ]
    endpoint = BASE + "/search"
    url = endpoint + "?" + urllib.parse.urlencode(params)
    post_data = urllib.parse.urlencode(params).encode("utf-8") if len(url) > POST_THRESHOLD else None
    body, _cached = fetch(url, source=SOURCE, query=query, cache=cache,
                          raw_log=raw_log, ext="json",
                          post_url=endpoint, post_data=post_data)
    try:
        doc = json.loads(body)
    except ValueError as exc:
        raise FetchError(url, None, "unparseable JSON: %s" % exc)
    hit_count = int(doc.get("hitCount") or 0)
    items = (doc.get("resultList") or {}).get("result") or []
    return [_norm_record(it) for it in items], hit_count


def oa_lookup(pmid: Optional[str] = None, doi: Optional[str] = None, *,
              cache: Optional[Cache] = None,
              raw_log: Optional[RawLog] = None) -> Optional[Dict]:
    """Look one article up by PMID (preferred) or DOI to learn its OA status."""
    if pmid:
        query = "EXT_ID:%s AND SRC:MED" % pmid
    elif doi:
        query = 'DOI:"%s"' % doi
    else:
        return None
    try:
        recs, _ = search(query, page_size=1, cache=cache, raw_log=raw_log)
    except FetchError:
        return None
    return recs[0] if recs else None


def fulltext_path(pmcid: str) -> str:
    return data_path("raw", "fulltext", "%s.xml" % pmcid.upper())


def fetch_fulltext(pmcid: str, *, cache: Optional[Cache] = None,
                   raw_log: Optional[RawLog] = None) -> Optional[str]:
    """Retrieve OA full text XML for a PMCID; returns the local path or None.

    Not every OA record has a full-text XML endpoint; a 404 is a real answer and
    is returned as None rather than retried or faked.
    """
    pmcid = pmcid.upper()
    path = fulltext_path(pmcid)
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    url = "%s/%s/fullTextXML" % (BASE, pmcid)
    try:
        body, _cached = fetch(url, source="europepmc_fulltext", query=pmcid,
                              cache=cache, raw_log=raw_log, ext="xml")
    except FetchError:
        return None
    if not body or body.lstrip().startswith("<error") or len(body) < 500:
        return None
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(body)
    return path
