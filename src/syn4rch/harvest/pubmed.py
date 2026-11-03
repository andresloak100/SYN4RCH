"""PubMed E-utilities client (ESearch + EFetch), XML parsed with xml.etree.

No email address or api_key is sent. NCBI's anonymous rate ceiling is respected
by the throttle in :mod:`syn4rch.harvest.cache` (0.5 s minimum spacing).

Fields extracted per record, matching docs/SCHEMA.md 3.1:
    pmid, pmcid, doi, title, abstract, journal, year, first_author, last_author,
    last_author_affiliation, publication_types, mesh_terms
Anything absent in the XML is ``None`` or an empty list. Nothing is inferred.
"""

from __future__ import annotations

import urllib.parse
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Tuple

from .cache import Cache, RawLog, fetch

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SOURCE = "pubmed"

# URLs longer than this are sent as POST; NCBI answers 414 to long GETs and the
# vocabulary-generated query blocks routinely run to several kB.
POST_THRESHOLD = 1800


def esearch(term: str, *, retmax: int = 200, mindate: Optional[str] = None,
            maxdate: Optional[str] = None, cache: Optional[Cache] = None,
            raw_log: Optional[RawLog] = None) -> Tuple[List[str], int]:
    """Run ESearch. Returns ``(pmids, total_count)``.

    ``total_count`` is PubMed's own Count for the query, which may exceed
    ``len(pmids)`` when retmax truncates. Both numbers are reported downstream.
    """
    params = [
        ("db", "pubmed"),
        ("term", term),
        ("retmax", str(retmax)),
        ("retmode", "xml"),
        ("sort", "relevance"),
    ]
    if mindate:
        params.append(("mindate", mindate))
    if maxdate:
        params.append(("maxdate", maxdate))
    if mindate or maxdate:
        params.append(("datetype", "pdat"))
    endpoint = EUTILS + "/esearch.fcgi"
    url = endpoint + "?" + urllib.parse.urlencode(params)
    post_data = urllib.parse.urlencode(params).encode("utf-8") if len(url) > POST_THRESHOLD else None
    body, _cached = fetch(url, source=SOURCE, query=term, cache=cache,
                          raw_log=raw_log, ext="xml",
                          post_url=endpoint, post_data=post_data)
    root = ET.fromstring(body)
    count_el = root.find("Count")
    total = int(count_el.text) if (count_el is not None and count_el.text and count_el.text.isdigit()) else 0
    pmids = [el.text for el in root.findall("./IdList/Id") if el.text]
    return pmids, total


def _text(el: Optional[ET.Element]) -> Optional[str]:
    if el is None:
        return None
    s = "".join(el.itertext()).strip()
    return s or None


def _year_of(article: ET.Element) -> Optional[int]:
    for path in ("./Journal/JournalIssue/PubDate/Year",
                 "./Journal/JournalIssue/PubDate/MedlineDate"):
        el = article.find(path)
        txt = _text(el)
        if txt:
            for token in txt.replace("-", " ").split():
                if len(token) == 4 and token.isdigit():
                    return int(token)
    return None


def _author_name(auth: ET.Element) -> Optional[str]:
    collective = _text(auth.find("CollectiveName"))
    if collective:
        return collective
    last = _text(auth.find("LastName"))
    initials = _text(auth.find("Initials"))
    fore = _text(auth.find("ForeName"))
    if last and fore:
        return "%s %s" % (fore, last)
    if last and initials:
        return "%s %s" % (initials, last)
    return last


def _affiliation(auth: ET.Element) -> Optional[str]:
    affs = [_text(a) for a in auth.findall("./AffiliationInfo/Affiliation")]
    affs = [a for a in affs if a]
    return affs[0] if affs else None


def parse_efetch(xml_text: str) -> List[Dict]:
    """Parse an EFetch PubmedArticleSet into records."""
    out: List[Dict] = []
    root = ET.fromstring(xml_text)
    for art in list(root.findall("./PubmedArticle")) + list(root.findall("./PubmedBookArticle")):
        cite = art.find("./MedlineCitation") or art.find("./BookDocument")
        if cite is None:
            continue
        pmid = _text(cite.find("PMID"))
        article = cite.find("./Article")

        doi = None
        pmcid = None
        for aid in art.findall("./PubmedData/ArticleIdList/ArticleId"):
            idtype = (aid.get("IdType") or "").lower()
            val = _text(aid)
            if not val:
                continue
            if idtype == "doi" and doi is None:
                doi = val.lower()
            elif idtype == "pmc" and pmcid is None:
                pmcid = val.upper()
        if doi is None and article is not None:
            for eid in article.findall("./ELocationID"):
                if (eid.get("EIdType") or "").lower() == "doi":
                    doi = (_text(eid) or "").lower() or None
                    break

        title = _text(article.find("ArticleTitle")) if article is not None else None
        if title is None:
            title = _text(cite.find("./ArticleTitle"))

        abstract_parts: List[str] = []
        if article is not None:
            for ab in article.findall("./Abstract/AbstractText"):
                label = ab.get("Label")
                txt = _text(ab)
                if not txt:
                    continue
                abstract_parts.append("%s: %s" % (label, txt) if label else txt)
        abstract = "\n".join(abstract_parts) or None

        journal = None
        year = None
        if article is not None:
            journal = (_text(article.find("./Journal/Title"))
                       or _text(article.find("./Journal/ISOAbbreviation")))
            year = _year_of(article)

        authors: List[ET.Element] = []
        if article is not None:
            authors = [a for a in article.findall("./AuthorList/Author")]
        first_author = _author_name(authors[0]) if authors else None
        last_author = _author_name(authors[-1]) if len(authors) >= 1 else None
        last_author_affiliation = _affiliation(authors[-1]) if authors else None
        if last_author_affiliation is None and authors:
            # fall back to any listed affiliation, recording nothing if none exist
            for auth in reversed(authors):
                aff = _affiliation(auth)
                if aff:
                    last_author_affiliation = aff
                    break

        pub_types = [t for t in (_text(x) for x in (
            article.findall("./PublicationTypeList/PublicationType") if article is not None else []
        )) if t]

        mesh: List[str] = []
        for mh in cite.findall("./MeshHeadingList/MeshHeading"):
            name = _text(mh.find("DescriptorName"))
            if name:
                mesh.append(name)

        out.append({
            "pmid": pmid,
            "pmcid": pmcid,
            "doi": doi,
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": year,
            "first_author": first_author,
            "last_author": last_author,
            "last_author_affiliation": last_author_affiliation,
            "publication_types": pub_types,
            "mesh_terms": mesh,
            "n_authors": len(authors),
            "retrieved_from": "pubmed",
        })
    return out


def efetch(pmids: List[str], *, batch_size: int = 100, cache: Optional[Cache] = None,
           raw_log: Optional[RawLog] = None) -> List[Dict]:
    """EFetch full records for a list of PMIDs, in batches."""
    records: List[Dict] = []
    for i in range(0, len(pmids), batch_size):
        chunk = [p for p in pmids[i:i + batch_size] if p]
        if not chunk:
            continue
        query = "efetch:" + ",".join(chunk)
        params = [
            ("db", "pubmed"),
            ("id", ",".join(chunk)),
            ("retmode", "xml"),
        ]
        endpoint = EUTILS + "/efetch.fcgi"
        url = endpoint + "?" + urllib.parse.urlencode(params)
        post_data = urllib.parse.urlencode(params).encode("utf-8") if len(url) > POST_THRESHOLD else None
        body, _cached = fetch(url, source=SOURCE, query=query, cache=cache,
                              raw_log=raw_log, ext="xml",
                              post_url=endpoint, post_data=post_data)
        try:
            records.extend(parse_efetch(body))
        except ET.ParseError as exc:
            # A malformed batch is reported, not silently dropped.
            records.append({
                "pmid": None, "pmcid": None, "doi": None,
                "title": None, "abstract": None, "journal": None, "year": None,
                "first_author": None, "last_author": None,
                "last_author_affiliation": None,
                "publication_types": [], "mesh_terms": [], "n_authors": 0,
                "retrieved_from": "pubmed",
                "_parse_error": "%s for pmids %s" % (exc, ",".join(chunk)),
            })
    return records
