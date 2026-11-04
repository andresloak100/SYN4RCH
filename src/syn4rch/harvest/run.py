"""Harvest runner.

    python3 -m syn4rch.harvest.run --max-queries 404 [--resume]

Pipeline: build the query matrix from the vocabularies -> issue each query
(PubMed ESearch / Europe PMC search) -> EFetch full bibliographic records ->
merge Europe PMC open-access status -> Crossref backfill for missing
bibliographic fields -> dedupe by DOI, then PMID, then normalised title ->
rule-based screening -> write

    data/interim/candidates.jsonl
    docs/QUERIES.md             (every query issued, with its result counts)
    reports/harvest_summary.md  (measured counts only)

Every number in the summary is computed from the records actually written.
Nothing is estimated. Queries that error out are recorded as errors with their
message and counted, not dropped.

Standard library only; Python 3.10 compatible.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from typing import Dict, List, Optional, Tuple

from ..schema.record import SCHEMA_VERSION
from . import crossref, europepmc, pubmed, queries as qmod
from .cache import Cache, FetchError, RawLog, data_path, repo_root
from .screen import screen_record
from ..vocab import load

STATE_PATH = data_path("interim", "harvest_state.json")
CANDIDATES_PATH = data_path("interim", "candidates.jsonl")
QUERIES_MD = os.path.join(repo_root(), "docs", "QUERIES.md")
SUMMARY_MD = os.path.join(repo_root(), "reports", "harvest_summary.md")

_NONALNUM = re.compile(r"[^a-z0-9]+")


def norm_title(title: Optional[str]) -> str:
    return _NONALNUM.sub(" ", (title or "").lower()).strip()


# --------------------------------------------------------------------------- #
# record store
# --------------------------------------------------------------------------- #

class Store:
    """Keyed record store that merges duplicates (DOI -> PMID -> title)."""

    def __init__(self) -> None:
        self.records: Dict[str, Dict] = {}
        self.by_doi: Dict[str, str] = {}
        self.by_pmid: Dict[str, str] = {}
        self.by_title: Dict[str, str] = {}
        self.merge_events: List[Tuple[str, str]] = []

    def _find(self, rec: Dict) -> Optional[str]:
        doi = (rec.get("doi") or "").lower() or None
        pmid = rec.get("pmid") or None
        title = norm_title(rec.get("title"))
        if doi and doi in self.by_doi:
            return self.by_doi[doi]
        if pmid and pmid in self.by_pmid:
            return self.by_pmid[pmid]
        if title and len(title) > 25 and title in self.by_title:
            return self.by_title[title]
        return None

    def _reindex(self, key: str, rec: Dict) -> None:
        doi = (rec.get("doi") or "").lower() or None
        pmid = rec.get("pmid") or None
        title = norm_title(rec.get("title"))
        if doi:
            self.by_doi[doi] = key
        if pmid:
            self.by_pmid[pmid] = key
        if title and len(title) > 25:
            self.by_title[title] = key

    def add(self, rec: Dict, *, qid: str, population: Optional[str],
            tier: str) -> str:
        existing_key = self._find(rec)
        if existing_key is None:
            key = ("doi:" + rec["doi"]) if rec.get("doi") else (
                ("pmid:" + rec["pmid"]) if rec.get("pmid") else
                ("title:" + norm_title(rec.get("title"))[:120]))
            rec = dict(rec)
            rec["study_id"] = key
            rec["found_by_queries"] = []
            rec["found_by_tiers"] = []
            rec["found_for_populations"] = []
            self.records[key] = rec
            self._reindex(key, rec)
        else:
            key = existing_key
            cur = self.records[key]
            for field, val in rec.items():
                if field in ("found_by_queries", "found_by_tiers",
                             "found_for_populations", "study_id"):
                    continue
                if not cur.get(field) and val:
                    cur[field] = val
            # Europe PMC OA flags must win when they are positive
            if rec.get("is_open_access"):
                cur["is_open_access"] = True
            if rec.get("pmcid") and not cur.get("pmcid"):
                cur["pmcid"] = rec["pmcid"]
            self.merge_events.append((key, rec.get("retrieved_from") or "?"))
            self._reindex(key, cur)
        cur = self.records[key]
        if qid not in cur["found_by_queries"]:
            cur["found_by_queries"].append(qid)
        if tier not in cur["found_by_tiers"]:
            cur["found_by_tiers"].append(tier)
        if population and population not in cur["found_for_populations"]:
            cur["found_for_populations"].append(population)
        return key


# --------------------------------------------------------------------------- #
# state (for --resume)
# --------------------------------------------------------------------------- #

def load_state() -> Dict:
    if os.path.exists(STATE_PATH):
        with open(STATE_PATH, "r", encoding="utf-8") as fh:
            return json.load(fh)
    return {"query_log": [], "done_qids": []}


def save_state(state: Dict) -> None:
    """Write the resume state. Never allowed to kill a run in progress."""
    os.makedirs(os.path.dirname(STATE_PATH), exist_ok=True)
    tmp = "%s.%d.tmp" % (STATE_PATH, os.getpid())
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(state, fh, indent=1, sort_keys=True)
        os.replace(tmp, STATE_PATH)
    except OSError as exc:
        sys.stderr.write("  warning: could not save resume state: %s\n" % exc)
        try:
            if os.path.exists(tmp):
                os.remove(tmp)
        except OSError:
            pass


# --------------------------------------------------------------------------- #
# main stages
# --------------------------------------------------------------------------- #

BROAD_TIERS = ("T3_region_only", "T4_modality_only", "E3_region_epmc")


def issue_queries(matrix: List[Dict], *, retmax: int, broad_retmax: int,
                  cache: Cache, raw_log: RawLog,
                  state: Dict, store: Store, verbose: bool = True,
                  prior_log: Optional[List[Dict]] = None,
                  prior_pmids: Optional[Dict[str, List[Dict]]] = None
                  ) -> Tuple[List[Dict], Dict[str, List[Dict]]]:
    """Issue every query in ``matrix``; return the query log entries.

    Broad queries (region-only and modality-only tiers) get ``broad_retmax``
    instead of ``retmax``: they are the net for populations missing from the
    vocabulary and have hit counts in the thousands, so capping them at the same
    retmax as a population-specific query would throw that net away.
    """
    log: List[Dict] = []
    done = set(state.get("done_qids") or [])
    prior_log = list(prior_log or [])
    pubmed_pmids: Dict[str, List[Dict]] = dict(prior_pmids or {})

    for i, q in enumerate(matrix, 1):
        qid = q["qid"]
        entry = {
            "qid": qid, "tier": q["tier"], "engine": q["engine"],
            "population": q["population"], "cluster": q["cluster"],
            "term": q["term"],
            "issued_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
        rm = broad_retmax if q["tier"] in BROAD_TIERS else retmax
        entry["retmax_used"] = rm
        try:
            if q["engine"] == "pubmed":
                pmids, total = pubmed.esearch(q["term"], retmax=rm,
                                              cache=cache, raw_log=raw_log)
                entry["engine_total_hits"] = total
                entry["ids_retrieved"] = len(pmids)
                entry["truncated_by_retmax"] = total > len(pmids)
                for pmid in pmids:
                    pubmed_pmids.setdefault(pmid, []).append(q)
            else:
                recs, hits = europepmc.search(q["term"], page_size=rm,
                                              cache=cache, raw_log=raw_log)
                entry["engine_total_hits"] = hits
                entry["ids_retrieved"] = len(recs)
                entry["truncated_by_retmax"] = hits > len(recs)
                for rec in recs:
                    store.add(rec, qid=qid, population=q["population"], tier=q["tier"])
            entry["status"] = "ok"
        except FetchError as exc:
            entry["status"] = "error"
            entry["error"] = str(exc)
            entry["engine_total_hits"] = None
            entry["ids_retrieved"] = 0
        log.append(entry)
        done.add(qid)
        if verbose and (i % 20 == 0 or i == len(matrix)):
            sys.stderr.write("  query %d/%d (%s) unique pmids so far: %d\n"
                             % (i, len(matrix), qid, len(pubmed_pmids)))
            sys.stderr.flush()
        if i % 25 == 0:
            state["done_qids"] = sorted(done)
            state["query_log"] = prior_log + log
            state["pubmed_pmids"] = {k: [x["qid"] for x in v]
                                     for k, v in pubmed_pmids.items()}
            save_state(state)

    state["done_qids"] = sorted(done)
    state["pubmed_pmids"] = {k: [q["qid"] for q in v] for k, v in pubmed_pmids.items()}
    state["query_log"] = prior_log + log
    save_state(state)
    return prior_log + log, pubmed_pmids


def fetch_pubmed_records(pubmed_pmids: Dict[str, List[Dict]], *, cache: Cache,
                         raw_log: RawLog, store: Store, verbose: bool = True) -> int:
    pmids = sorted(pubmed_pmids.keys())
    n = 0
    batch = 100
    for i in range(0, len(pmids), batch):
        chunk = pmids[i:i + batch]
        try:
            recs = pubmed.efetch(chunk, batch_size=batch, cache=cache, raw_log=raw_log)
        except FetchError as exc:
            sys.stderr.write("  efetch failed for %d pmids: %s\n" % (len(chunk), exc))
            continue
        for rec in recs:
            if rec.get("_parse_error"):
                sys.stderr.write("  parse error: %s\n" % rec["_parse_error"])
                continue
            pmid = rec.get("pmid")
            sources = pubmed_pmids.get(pmid or "", [])
            if not sources:
                store.add(rec, qid="(unattributed)", population=None, tier="unknown")
            for q in sources:
                store.add(rec, qid=q["qid"], population=q["population"], tier=q["tier"])
            n += 1
        if verbose:
            sys.stderr.write("  efetched %d/%d pmids\n" % (min(i + batch, len(pmids)), len(pmids)))
            sys.stderr.flush()
    return n


def annotate_open_access(store: Store, *, cache: Cache, raw_log: RawLog,
                         max_lookups: int, verbose: bool = True) -> Dict[str, int]:
    """Batch Europe PMC EXT_ID lookups to learn OA status for PubMed-only records."""
    todo = [r for r in store.records.values()
            if r.get("is_open_access") is None and r.get("pmid")]
    todo.sort(key=lambda r: r.get("pmid") or "")
    stats = {"looked_up": 0, "resolved": 0, "batches": 0}
    chunk_size = 50
    budget = max_lookups
    for i in range(0, len(todo), chunk_size):
        if budget <= 0:
            break
        chunk = todo[i:i + chunk_size]
        budget -= len(chunk)
        query = " OR ".join("EXT_ID:%s" % r["pmid"] for r in chunk)
        try:
            recs, _hits = europepmc.search("(%s) AND SRC:MED" % query,
                                           page_size=chunk_size,
                                           cache=cache, raw_log=raw_log)
        except FetchError as exc:
            sys.stderr.write("  epmc OA batch failed: %s\n" % exc)
            continue
        stats["batches"] += 1
        stats["looked_up"] += len(chunk)
        found = {r.get("pmid"): r for r in recs if r.get("pmid")}
        for rec in chunk:
            hit = found.get(rec["pmid"])
            if hit is None:
                rec["is_open_access"] = False
                rec["oa_lookup"] = "not_in_europepmc"
                continue
            rec["is_open_access"] = bool(hit.get("is_open_access"))
            rec["in_epmc"] = hit.get("in_epmc")
            if hit.get("pmcid") and not rec.get("pmcid"):
                rec["pmcid"] = hit["pmcid"]
            if hit.get("doi") and not rec.get("doi"):
                rec["doi"] = hit["doi"]
            rec["oa_lookup"] = "europepmc"
            stats["resolved"] += 1
        if verbose:
            sys.stderr.write("  OA lookups: %d/%d\n" % (min(i + chunk_size, len(todo)), len(todo)))
            sys.stderr.flush()
    for rec in store.records.values():
        if rec.get("is_open_access") is None:
            rec["oa_lookup"] = "not_attempted"
    return stats


def backfill_crossref(store: Store, *, cache: Cache, raw_log: RawLog,
                      max_lookups: int) -> int:
    needed = [r for r in store.records.values()
              if r.get("doi") and not all(r.get(f) for f in
                                          ("journal", "year", "last_author"))]
    needed.sort(key=lambda r: r.get("doi") or "")
    n = 0
    for rec in needed[:max_lookups]:
        before = dict(rec)
        crossref.backfill(rec, cache=cache, raw_log=raw_log)
        if rec.get("crossref_backfilled"):
            n += 1
        elif before is not rec:
            pass
    return n


def retrieve_fulltext(store: Store, *, cache: Cache, raw_log: RawLog,
                      max_files: int, verbose: bool = True) -> Dict[str, int]:
    """Retrieve OA full text for screened-in records, highest score first."""
    cands = [r for r in store.records.values()
             if r.get("pmcid") and r.get("is_open_access")
             and (r.get("screen") or {}).get("decision") == "include"]
    cands.sort(key=lambda r: -((r.get("screen") or {}).get("score") or 0))
    stats = {"attempted": 0, "retrieved": 0}
    for rec in cands[:max_files]:
        stats["attempted"] += 1
        path = europepmc.fetch_fulltext(rec["pmcid"], cache=cache, raw_log=raw_log)
        if path:
            rec["full_text_path"] = os.path.relpath(path, repo_root())
            rec["full_text_retrieved"] = True
            stats["retrieved"] += 1
        else:
            rec["full_text_retrieved"] = False
        if verbose and stats["attempted"] % 20 == 0:
            sys.stderr.write("  fulltext %d/%d\n" % (stats["attempted"], min(len(cands), max_files)))
            sys.stderr.flush()
    return stats


# --------------------------------------------------------------------------- #
# outputs
# --------------------------------------------------------------------------- #

def write_candidates(store: Store) -> int:
    os.makedirs(os.path.dirname(CANDIDATES_PATH), exist_ok=True)
    rows = sorted(store.records.values(),
                  key=lambda r: (-(r.get("screen") or {}).get("score", 0),
                                 r.get("study_id") or ""))
    with open(CANDIDATES_PATH, "w", encoding="utf-8") as fh:
        for rec in rows:
            out = dict(rec)
            out["full_text_available"] = bool(out.get("full_text_retrieved"))
            out["schema_version"] = SCHEMA_VERSION
            fh.write(json.dumps(out, sort_keys=True) + "\n")
    return len(rows)


def write_queries_md(log: List[Dict], *, run_started: str, run_finished: str,
                     args: argparse.Namespace, matrix_size: int) -> None:
    os.makedirs(os.path.dirname(QUERIES_MD), exist_ok=True)
    ok = [e for e in log if e.get("status") == "ok"]
    errs = [e for e in log if e.get("status") != "ok"]
    trunc = [e for e in ok if e.get("truncated_by_retmax")]
    lines: List[str] = []
    lines.append("# SYN4RCH query log")
    lines.append("")
    lines.append("Generated by `python3 -m syn4rch.harvest.run`. Every query below was")
    lines.append("actually issued; `engine_total_hits` is the search engine's own total for")
    lines.append("that query and `ids_retrieved` is how many ids this run took from it")
    lines.append("(capped by `--retmax`). Where they differ the query is marked truncated,")
    lines.append("and the corpus is a relevance-ranked sample of that query, not all of it.")
    lines.append("")
    lines.append("## Run")
    lines.append("")
    lines.append("| field | value |")
    lines.append("|---|---|")
    lines.append("| started (UTC) | %s |" % run_started)
    lines.append("| finished (UTC) | %s |" % run_finished)
    lines.append("| queries in matrix | %d |" % matrix_size)
    lines.append("| queries issued | %d |" % len(log))
    lines.append("| queries ok | %d |" % len(ok))
    lines.append("| queries errored | %d |" % len(errs))
    lines.append("| queries truncated by retmax | %d |" % len(trunc))
    lines.append("| --retmax (population-specific queries) | %d |" % args.retmax)
    lines.append("| --broad-retmax (region-only / modality-only) | %d |" % args.broad_retmax)
    lines.append("| --max-queries | %s |" % args.max_queries)
    lines.append("")
    lines.append("## Curation decisions recorded here")
    lines.append("")
    lines.append("(`docs/DECISIONS.md` is owned by another agent, so harvest/vocabulary")
    lines.append("decisions are logged in this file instead.)")
    lines.append("")
    lines.append("1. **All `allen_id` values in `data/vocab/regions.json` are `null`.** The")
    lines.append("   curator did not have verified Allen Mouse Brain Atlas structure ids to")
    lines.append("   hand and would have had to guess the integers. A wrong atlas id is worse")
    lines.append("   than a missing one, so all 23 are null and must be filled from a")
    lines.append("   downloaded copy of the Allen structure graph.")
    lines.append("2. **`higher_is` is `null` for all but four outcomes.** Direction of benefit")
    lines.append("   is context-dependent for essentially every Tier A physiological variable.")
    lines.append("   Only median lifespan, maximum lifespan, grip strength (better) and")
    lines.append("   frailty index (worse) carry a value.")
    lines.append("3. **Cre lines are listed only where the curator was confident the driver")
    lines.append("   exists under approximately that name.** GHRH, TRH, Tbx3, VMH/DMH BDNF and")
    lines.append("   the QRFP preoptic population therefore have empty `cre_lines` with the")
    lines.append("   uncertainty stated in `notes`, rather than plausible-looking invented")
    lines.append("   allele names. No JAX stock numbers appear anywhere.")
    lines.append("4. **`resolve()` has no fuzzy matching.** Exact and curated-synonym only, so")
    lines.append("   `None` is frequent and correct. Ambiguous acronyms (`LH`, `ARC`, `ME`,")
    lines.append("   `ZI`, `PH`, `GH`, `T3`/`T4`, ...) are excluded from the screener's")
    lines.append("   matchers entirely; see `AMBIGUOUS_ACRONYMS` in `screen.py`.")
    lines.append("5. **Query term blocks are built round-robin across vocabulary entries**, not")
    lines.append("   as a flat truncated list. A flat list cut at N silently searched only for")
    lines.append("   opto/chemo methods and glucose outcomes; this was found and fixed before")
    lines.append("   the real run.")
    lines.append("6. **Queries are POSTed to E-utilities when the URL exceeds ~1.8 kB.** The")
    lines.append("   vocabulary-generated blocks run to several kB and NCBI answers HTTP 414")
    lines.append("   to a GET. The canonical expanded URL is still what gets logged and cached.")
    lines.append("7. **`data/raw/` keeps every API response verbatim**, plus a sqlite response")
    lines.append("   cache in `data/interim/`. Together these run to hundreds of MB after a")
    lines.append("   full harvest. Provenance is mandatory (SCHEMA.md section 4) so they are")
    lines.append("   not pruned here; `.gitignore` is owned by another agent, so whether they")
    lines.append("   are tracked is that file's decision, not this pipeline's.")
    lines.append("8. **`data/raw/fulltext/` can hold more files than the summary's full-text")
    lines.append("   count.** Earlier runs with a different `--retmax` selected a different")
    lines.append("   top-N of screened-in papers; the summary counts only records in the")
    lines.append("   current `candidates.jsonl`, while the directory is append-only.")
    lines.append("9. **No email address is sent to any API** (no `email=`, no `tool=`, no")
    lines.append("   Crossref `mailto`). We stay in the anonymous pools and inside a 2 req/s")
    lines.append("   NCBI limit instead of inventing contact details.")
    lines.append("")
    lines.append("## Queries issued")
    lines.append("")
    lines.append("| # | qid | tier | engine | population | cluster | engine total hits | retmax | ids retrieved | truncated | status |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for i, e in enumerate(log, 1):
        lines.append("| %d | `%s` | %s | %s | %s | %s | %s | %s | %d | %s | %s |" % (
            i, e["qid"], e["tier"], e["engine"], e.get("population") or "-",
            e.get("cluster") or "-",
            "-" if e.get("engine_total_hits") is None else e["engine_total_hits"],
            e.get("retmax_used") or "-",
            e.get("ids_retrieved") or 0,
            "yes" if e.get("truncated_by_retmax") else "no",
            e.get("status")))
    if errs:
        lines.append("")
        lines.append("### Errors")
        lines.append("")
        for e in errs:
            lines.append("* `%s`: %s" % (e["qid"], e.get("error")))
    lines.append("")
    lines.append("## Exact query strings")
    lines.append("")
    lines.append("Verbatim, so any query can be re-issued by hand.")
    lines.append("")
    for e in log:
        lines.append("### `%s` (%s, %s)" % (e["qid"], e["tier"], e["engine"]))
        lines.append("")
        lines.append("```")
        lines.append(e["term"])
        lines.append("```")
        lines.append("")
    with open(QUERIES_MD, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def summarise(store: Store, log: List[Dict], *, run_started: str, run_finished: str,
              args: argparse.Namespace, oa_stats: Dict, ft_stats: Dict,
              crossref_n: int, cache: Cache) -> str:
    recs = list(store.records.values())
    pops = load("populations")

    decisions = {"include": 0, "uncertain": 0, "exclude": 0}
    reasons: Dict[str, int] = {}
    excl_reasons: Dict[str, int] = {}
    flags: Dict[str, int] = {}
    with_readout = 0
    outcome_counts: Dict[str, int] = {}
    for r in recs:
        s = r.get("screen") or {}
        d = s.get("decision")
        if d in decisions:
            decisions[d] += 1
        reasons[s.get("primary_reason") or "?"] = reasons.get(s.get("primary_reason") or "?", 0) + 1
        if d == "exclude":
            key = s.get("primary_reason") or "?"
            excl_reasons[key] = excl_reasons.get(key, 0) + 1
        for f in s.get("flags") or []:
            flags[f] = flags.get(f, 0) + 1
        if s.get("matched_outcomes"):
            with_readout += 1
            for oid in s["matched_outcomes"]:
                outcome_counts[oid] = outcome_counts.get(oid, 0) + 1

    # population distribution: papers whose SCREEN matched that population
    pop_screen: Dict[str, int] = {pid: 0 for pid in pops}
    pop_screen_incl: Dict[str, int] = {pid: 0 for pid in pops}
    pop_query: Dict[str, int] = {pid: 0 for pid in pops}
    for r in recs:
        s = r.get("screen") or {}
        for pid in s.get("matched_populations") or []:
            pop_screen[pid] = pop_screen.get(pid, 0) + 1
            if s.get("decision") == "include":
                pop_screen_incl[pid] = pop_screen_incl.get(pid, 0) + 1
        for pid in r.get("found_for_populations") or []:
            if pid in pop_query:
                pop_query[pid] += 1

    oa_yes = sum(1 for r in recs if r.get("is_open_access"))
    oa_no = sum(1 for r in recs if r.get("is_open_access") is False)
    oa_unknown = len(recs) - oa_yes - oa_no
    with_pmcid = sum(1 for r in recs if r.get("pmcid"))
    ft_local = sum(1 for r in recs if r.get("full_text_retrieved"))
    biblio_complete = sum(1 for r in recs if all(r.get(f) for f in
                          ("title", "journal", "year", "first_author", "last_author")))
    with_aff = sum(1 for r in recs if r.get("last_author_affiliation"))
    with_doi = sum(1 for r in recs if r.get("doi"))
    with_abstract = sum(1 for r in recs if r.get("abstract"))

    ok_q = [e for e in log if e.get("status") == "ok"]
    trunc_q = [e for e in ok_q if e.get("truncated_by_retmax")]

    L: List[str] = []
    L.append("# SYN4RCH harvest summary")
    L.append("")
    L.append("Every number below is computed from `data/interim/candidates.jsonl` and")
    L.append("`docs/QUERIES.md`, both written by this run. Nothing is estimated.")
    L.append("")
    L.append("**Caveat that governs the whole table.** This is a title/abstract screen. It")
    L.append("cannot tell you how many *extractable* (population x perturbation x outcome)")
    L.append("records exist -- only full-text curation can. These are candidate papers, an")
    L.append("UPPER BOUND on usable records' source papers, and the `include` count in")
    L.append("particular will shrink once full texts are read.")
    L.append("")
    L.append("## Run")
    L.append("")
    L.append("| field | value |")
    L.append("|---|---|")
    L.append("| started (UTC) | %s |" % run_started)
    L.append("| finished (UTC) | %s |" % run_finished)
    L.append("| queries issued | %d |" % len(log))
    L.append("| queries ok / errored | %d / %d |" % (len(ok_q), len(log) - len(ok_q)))
    L.append("| queries truncated by retmax (%d / broad %d) | %d |" % (
        args.retmax, args.broad_retmax, len(trunc_q)))
    L.append("| unique papers after dedupe | %d |" % len(recs))
    L.append("| duplicate merge events | %d |" % len(store.merge_events))
    L.append("| cached HTTP responses in sqlite | %d |" % cache.stats()["rows"])
    L.append("| cache hits / misses this run | %d / %d |" % (cache.hits, cache.misses))
    L.append("")
    L.append("## Screening")
    L.append("")
    L.append("| decision | papers |")
    L.append("|---|---|")
    for d in ("include", "uncertain", "exclude"):
        L.append("| %s | %d |" % (d, decisions[d]))
    L.append("| **total** | **%d** |" % len(recs))
    L.append("")
    L.append("### Exclusion reasons (top)")
    L.append("")
    L.append("| reason | papers |")
    L.append("|---|---|")
    for key, n in sorted(excl_reasons.items(), key=lambda kv: -kv[1]):
        L.append("| `%s` | %d |" % (key, n))
    L.append("")
    L.append("### Flags (recorded, never used to exclude)")
    L.append("")
    if flags:
        L.append("| flag | papers |")
        L.append("|---|---|")
        for key, n in sorted(flags.items(), key=lambda kv: -kv[1]):
            L.append("| `%s` | %d |" % (key, n))
    else:
        L.append("None.")
    L.append("")
    L.append("## Physiological readout mentions")
    L.append("")
    L.append("| field | value |")
    L.append("|---|---|")
    L.append("| candidate papers whose title/abstract mentions >=1 vocabulary readout | %d |" % with_readout)
    L.append("| papers with an abstract at all | %d |" % with_abstract)
    L.append("")
    L.append("Per outcome id (papers mentioning it; a paper can mention many):")
    L.append("")
    L.append("| outcome_id | papers |")
    L.append("|---|---|")
    outs = load("outcomes")
    for oid in outs:
        L.append("| `%s` | %d |" % (oid, outcome_counts.get(oid, 0)))
    L.append("")
    L.append("## Unique papers per population vocabulary entry")
    L.append("")
    L.append("`matched_in_text` = the screener found one of that entry's terms in the")
    L.append("title/abstract/MeSH. `found_by_its_queries` = the paper came back from a query")
    L.append("built from that entry. Zeros are kept, and they are the point of the table.")
    L.append("")
    L.append("| population_id | curation_confidence | matched_in_text | of which `include` | found_by_its_queries |")
    L.append("|---|---|---|---|---|")
    for pid, entry in sorted(pops.items(), key=lambda kv: -pop_screen.get(kv[0], 0)):
        L.append("| `%s` | %s | %d | %d | %d |" % (
            pid, entry.get("curation_confidence"), pop_screen.get(pid, 0),
            pop_screen_incl.get(pid, 0), pop_query.get(pid, 0)))
    zero = [p for p in pops if pop_screen.get(p, 0) == 0]
    L.append("")
    L.append("Populations with **zero** text matches in the whole corpus: %s" %
             (", ".join("`%s`" % p for p in zero) if zero else "none"))
    L.append("")
    L.append("## Full text availability")
    L.append("")
    L.append("| field | papers |")
    L.append("|---|---|")
    L.append("| open access per Europe PMC | %d |" % oa_yes)
    L.append("| not open access per Europe PMC | %d |" % oa_no)
    L.append("| OA status unknown (lookup budget not reached) | %d |" % oa_unknown)
    L.append("| has a PMCID | %d |" % with_pmcid)
    L.append("| OA full text XML actually downloaded to data/raw/fulltext/ | %d |" % ft_local)
    L.append("| full-text retrieval attempted / succeeded | %d / %d |" % (
        ft_stats.get("attempted", 0), ft_stats.get("retrieved", 0)))
    L.append("| abstract only (no full text retrieved) | %d |" % (len(recs) - ft_local))
    L.append("")
    L.append("Open-access lookups: %d papers in %d Europe PMC batch queries, %d resolved." % (
        oa_stats.get("looked_up", 0), oa_stats.get("batches", 0), oa_stats.get("resolved", 0)))
    L.append("")
    L.append("## Bibliographic completeness (SCHEMA.md 3.1)")
    L.append("")
    L.append("| field present | papers |")
    L.append("|---|---|")
    L.append("| DOI | %d |" % with_doi)
    L.append("| title+journal+year+first_author+last_author | %d |" % biblio_complete)
    L.append("| last_author_affiliation (needed for lab-level leakage control) | %d |" % with_aff)
    L.append("| Crossref backfills applied | %d |" % crossref_n)
    L.append("")
    L.append("## Known weaknesses of this search")
    L.append("")
    L.append("* retmax (%d, broad %d) truncated %d of %d successful queries, so the corpus is a" % (
        args.retmax, args.broad_retmax, len(trunc_q), len(ok_q)))
    L.append("  relevance-ranked sample of the matrix, not its closure.")
    L.append("* Phrase matching on title/abstract misses papers that name a population only")
    L.append("  in the methods section, which is common for Cre-line targeting.")
    L.append("* The screener's `site` criterion is satisfied by a region term alone, so")
    L.append("  region-level manipulations that cannot yield a population-specific record")
    L.append("  still screen in.")
    L.append("* Ambiguous acronyms are excluded from matching, which costs recall for papers")
    L.append("  that only ever write `PVN` or `ARC` -- but the alternative is false matches on")
    L.append("  English words.")
    L.append("* **Europe PMC (the `E1`/`E3` queries) searches FULL TEXT, not just the")
    L.append("  abstract**, so those queries are far less precise than the PubMed `[tiab]`")
    L.append("  queries: a paper that mentions a population once in its methods matches. This")
    L.append("  is why `found_by_its_queries` is much larger than `matched_in_text` for most")
    L.append("  populations. `matched_in_text` is the number to trust; `found_by_its_queries`")
    L.append("  is only a record of which query surfaced the paper.")
    L.append("* PubMed silently splits a quoted phrase it does not hold in its phrase index")
    L.append("  into separate terms, which loosens some blocks in ways the log cannot show.")
    L.append("* The `include` bucket is inflated by papers where the modality, the site and")
    L.append("  the readout are all mentioned but are not connected to each other -- e.g. a")
    L.append("  pharmacology paper citing optogenetic work in its introduction.")
    L.append("* No non-English literature, no books/theses, no conference abstracts.")
    L.append("")
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="python3 -m syn4rch.harvest.run")
    ap.add_argument("--max-queries", type=int, default=50,
                    help="number of queries from the matrix to issue (priority order)")
    ap.add_argument("--retmax", type=int, default=50,
                    help="max ids taken per population-specific query")
    ap.add_argument("--broad-retmax", type=int, default=400,
                    help="max ids taken per region-only / modality-only query")
    ap.add_argument("--resume", action="store_true",
                    help="skip queries already recorded as done in data/interim/harvest_state.json")
    ap.add_argument("--no-epmc", action="store_true", help="skip Europe PMC queries")
    ap.add_argument("--max-oa-lookups", type=int, default=20000)
    ap.add_argument("--max-crossref", type=int, default=150)
    ap.add_argument("--max-fulltext", type=int, default=60)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    verbose = not args.quiet

    run_started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    t0 = time.time()

    cache = Cache()
    raw_log = RawLog()
    store = Store()
    state = load_state() if args.resume else {"query_log": [], "done_qids": []}

    matrix = qmod.build_matrix(include_epmc=not args.no_epmc)
    matrix_size = len(matrix)
    selected = matrix[:args.max_queries]
    prior_log: List[Dict] = []
    prior_pmids: Dict[str, List[Dict]] = {}
    if args.resume:
        done = set(state.get("done_qids") or [])
        by_qid = {q["qid"]: q for q in matrix}
        selected_qids = set(q["qid"] for q in selected)
        # Europe PMC results live only in memory, so those queries are always
        # re-issued; the sqlite cache makes that free.
        skippable = set(qid for qid in done
                        if (by_qid.get(qid) or {}).get("engine") == "pubmed")
        prior_log = [e for e in (state.get("query_log") or [])
                     if e.get("qid") in skippable and e.get("qid") in selected_qids]
        for pmid, qids in (state.get("pubmed_pmids") or {}).items():
            qs = [by_qid[q] for q in qids if q in by_qid and q in skippable]
            if qs:
                prior_pmids[pmid] = qs
        selected = [q for q in selected if q["qid"] not in skippable]
        if verbose:
            sys.stderr.write("resume: %d queries restored from state (%d pmids), "
                             "%d still to issue\n"
                             % (len(prior_log), len(prior_pmids), len(selected)))

    if verbose:
        sys.stderr.write("matrix: %d queries, issuing %d (retmax=%d)\n"
                         % (matrix_size, len(selected), args.retmax))

    log, pubmed_pmids = issue_queries(selected, retmax=args.retmax,
                                      broad_retmax=args.broad_retmax, cache=cache,
                                      raw_log=raw_log, state=state, store=store,
                                      verbose=verbose, prior_log=prior_log,
                                      prior_pmids=prior_pmids)
    if verbose:
        sys.stderr.write("efetching %d unique pmids\n" % len(pubmed_pmids))
    fetch_pubmed_records(pubmed_pmids, cache=cache, raw_log=raw_log, store=store,
                         verbose=verbose)

    if verbose:
        sys.stderr.write("open-access lookups\n")
    oa_stats = annotate_open_access(store, cache=cache, raw_log=raw_log,
                                    max_lookups=args.max_oa_lookups, verbose=verbose)

    if verbose:
        sys.stderr.write("crossref backfill\n")
    crossref_n = backfill_crossref(store, cache=cache, raw_log=raw_log,
                                   max_lookups=args.max_crossref)

    if verbose:
        sys.stderr.write("screening %d records\n" % len(store.records))
    for rec in store.records.values():
        rec["screen"] = screen_record(rec)

    if verbose:
        sys.stderr.write("retrieving OA full text (max %d)\n" % args.max_fulltext)
    ft_stats = retrieve_fulltext(store, cache=cache, raw_log=raw_log,
                                 max_files=args.max_fulltext, verbose=verbose)

    n = write_candidates(store)
    run_finished = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    write_queries_md(log, run_started=run_started, run_finished=run_finished,
                     args=args, matrix_size=matrix_size)
    os.makedirs(os.path.dirname(SUMMARY_MD), exist_ok=True)
    with open(SUMMARY_MD, "w", encoding="utf-8") as fh:
        fh.write(summarise(store, log, run_started=run_started,
                           run_finished=run_finished, args=args, oa_stats=oa_stats,
                           ft_stats=ft_stats, crossref_n=crossref_n, cache=cache))

    state["last_run"] = {"started": run_started, "finished": run_finished,
                         "queries": len(log), "records": n}
    save_state(state)
    cache.close()
    if verbose:
        sys.stderr.write("done: %d candidate papers in %.1f min\n"
                         % (n, (time.time() - t0) / 60.0))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
