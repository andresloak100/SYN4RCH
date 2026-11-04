"""Select and package a reproducible curation pilot.

    python3 -m syn4rch.curate.pilot --populations arc_agrp,arc_pomc --papers 50 --batches 5

"Which papers were curated, and why those?" must have an exact answer, so
selection is deterministic and documented rather than ad hoc.

Selection rule
--------------
1. Eligible papers are those with at least one extraction unit in a
   result-bearing section (Results or Methods) for one of the requested
   populations. Introduction and Discussion units never make a paper eligible:
   `docs/CURATION_PROTOCOL.md` 2 requires the reader to discard them.
2. Eligible papers are split into strata by which requested populations they
   cover, and papers are drawn from each stratum in proportion to its size, so
   the pilot's mix matches the corpus's.
3. Within a stratum, order is by `sha256(seed | study_id)`. This is an
   arbitrary but fixed order: it is not correlated with how many units, how
   many outcomes or how recent a paper is, so the measured records-per-paper
   yield is an unbiased estimate for that stratum. Sorting by unit count
   instead would select the richest papers and inflate the yield, which is the
   one number this pilot exists to measure.

What a batch contains
---------------------
Every result-bearing unit of each selected paper, including units for
populations other than the requested ones: once a passage is being read, any
observation it reports is admissible, and discarding them would waste reading.
Each paper also carries a `record_stub` with its bibliographic and provenance
fields already filled from the harvest, so a curator never retypes an
identifier. Fields a curator must supply are listed in `fields_to_complete`.

Standard library only (DECISIONS.md D-001). Reads data/, writes data/interim/.
"""

from __future__ import annotations

import argparse
import collections
import datetime as _dt
import hashlib
import json
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

from ..schema.rules import perturbation_class_from_method
from ..vocab.core import load as vocab_load

RESULT_SECTIONS = ("Results", "Methods")

# Filled from the harvest; a curator must never retype these.
STUB_FROM_HARVEST = ("study_id", "pmid", "pmcid", "doi", "title", "journal", "year",
                     "first_author", "last_author", "last_author_affiliation",
                     "retrieved_from", "full_text_available", "schema_version")

#: Fields the stub fills beyond :data:`STUB_FROM_HARVEST` -- provenance and the
#: Tier A constant. A curator must not retype these either.
STUB_CONSTANTS = ("extraction_method", "review_status", "evidence_class")

#: Fields no curator may fill, because something else owns them.
#:
#: ``record_id`` is computed from the record (SCHEMA.md 3.1).
#: ``perturbation_class`` and ``perturbation_polarity`` are derived by rule from
#: ``method`` and ``valence`` at admission (D-025); a curator writing them by
#: hand is writing a guess where a total function exists.
#: ``conflicts_with`` is written by QC, never by a reader.
#: ``extracted_by`` and ``extracted_at`` identify the run, and the runner sets
#: them.
NOT_THE_CURATORS_JOB = (
    "record_id",
    "perturbation_class",
    "perturbation_polarity",
    "conflicts_with",
    "extracted_by",
    "extracted_at",
)


def fields_to_complete() -> Tuple[str, ...]:
    """Every schema field a curator must supply, derived from the schema.

    This list used to be written out by hand and drifted: it still described
    schema v0.1.0 after v0.2.0 added the perturbation axis and v0.3.0 added the
    condition dimension, so ``challenge``, ``challenge_detail``, ``age_group``,
    ``condition_evidence``, the five ``stim_*`` fields, ``region_id``,
    ``projection_target_id``, ``tissue_id``, ``ci_low``/``ci_high``,
    ``genetic_background`` and ``confidence`` were all missing from it. Three
    separate readers hit that in one pass: each was told to record the condition
    dimension and handed a field list that did not mention it, and two of them
    put a stated stimulation frequency and a resolvable projection target into
    ``extraction_notes`` because the list offered nowhere else to put them.

    Deriving it means the list cannot drift again, and the partition is total:
    every field of :class:`~syn4rch.schema.record.TierARecord` is in the stub,
    in :data:`NOT_THE_CURATORS_JOB`, or here. A test asserts that.
    """
    from ..schema.record import TierARecord

    taken = set(STUB_FROM_HARVEST) | set(STUB_CONSTANTS) | set(NOT_THE_CURATORS_JOB)
    return tuple(f for f in TierARecord.field_names() if f not in taken)


def _repo_root() -> str:
    try:
        from ..harvest.cache import repo_root

        return repo_root()
    except Exception:
        return os.path.abspath(os.getcwd())


def _default(*parts: str) -> str:
    return os.path.join(_repo_root(), *parts)


def _order_key(seed: str, study_id: str) -> str:
    return hashlib.sha256(("%s|%s" % (seed, study_id)).encode("utf-8")).hexdigest()


def load_units(path: str) -> List[Dict]:
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def load_papers(path: str) -> Dict[str, Dict]:
    papers: Dict[str, Dict] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            sid = rec.get("study_id")
            if sid:
                papers[sid] = rec
    return papers


#: ``--populations`` value meaning "every population the units resolve".
ALL_POPULATIONS = "all"


def eligible(units: Sequence[Dict], populations: Sequence[str]) -> Dict[str, Dict]:
    """``study_id -> {'covers': set, 'units': [...]}`` for result-bearing units.

    ``populations == [ALL_POPULATIONS]`` accepts every resolved population
    instead of a named list. That is not a convenience: restricting a selection
    to the populations the corpus already covers most would guarantee that new
    papers land in existing cells, which is the very thing D-048's secondary
    metric is supposed to *measure*. A metric satisfied by the selection rule
    measures the selection rule.
    """
    want = set(populations)
    take_all = ALL_POPULATIONS in want
    result_units = [u for u in units if u.get("section_label") in RESULT_SECTIONS]
    by_study: Dict[str, List[Dict]] = collections.defaultdict(list)
    for u in result_units:
        by_study[u["study_id"]].append(u)
    out: Dict[str, Dict] = {}
    for sid, us in by_study.items():
        present = {u.get("candidate_population_id") for u in us if u.get("candidate_population_id")}
        covers = present if take_all else (present & want)
        if covers:
            out[sid] = {"covers": covers, "units": us}
    return out


def method_class_of(method_id: Optional[str]) -> Optional[str]:
    """Perturbation class of a *screening* method id, by the project's own chain.

    The units carry a finer vocabulary than the schema does -- ``opto_chr2``,
    ``chemo_hm3dq`` -- so the class is reached in two documented steps:
    ``data/vocab/methods.json`` maps the screening id to a schema ``method``,
    and :func:`syn4rch.schema.rules.perturbation_class_from_method` maps that to
    a class. Nothing here invents a third vocabulary, so "circuit" in a
    selection means exactly what it means in the claim key (SCHEMA.md 3.4b).
    """
    if not method_id:
        return None
    try:
        entry = vocab_load("methods").get(method_id)
    except Exception:  # pragma: no cover - vocab absent in unit tests
        return None
    if not entry:
        return None
    cls = perturbation_class_from_method(entry.get("method"))
    return getattr(cls, "value", cls)


def filter_units_by_class(units: Sequence[Dict], method_class: str) -> List[Dict]:
    """Units whose candidate method belongs to one perturbation class."""
    return [u for u in units
            if method_class_of(u.get("candidate_method_id")) == method_class]


def drop_curated_studies(
    elig: Dict[str, Dict], curated_path: str
) -> Tuple[Dict[str, Dict], int]:
    """Remove studies that already have records. Returns ``(kept, dropped)``.

    A pilot that re-selects a curated paper spends a reader on work already
    done, and -- worse -- would propose records whose ``record_id`` collides
    with one already in the dataset.
    """
    curated = set()
    with open(curated_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                sid = json.loads(line).get("study_id")
                if sid:
                    curated.add(sid)
    kept = {s: v for s, v in elig.items() if s not in curated}
    return kept, len(elig) - len(kept)


def select(elig: Dict[str, Dict], *, n_papers: int, populations: Sequence[str],
           seed: str) -> Dict:
    """Proportional stratified selection, deterministic in ``seed``."""
    strata: Dict[str, List[str]] = collections.defaultdict(list)
    for sid, info in elig.items():
        key = "+".join(sorted(info["covers"]))
        strata[key].append(sid)
    for key in strata:
        strata[key].sort(key=lambda s: _order_key(seed, s))

    total = len(elig)
    n_papers = min(n_papers, total)
    # Largest-remainder apportionment, so the quotas sum exactly to n_papers.
    quotas: Dict[str, int] = {}
    remainders = []
    for key, sids in strata.items():
        exact = n_papers * len(sids) / total
        quotas[key] = int(exact)
        remainders.append((exact - int(exact), key))
    short = n_papers - sum(quotas.values())
    for _, key in sorted(remainders, key=lambda t: (-t[0], t[1]))[:short]:
        quotas[key] += 1

    chosen: List[str] = []
    for key in sorted(strata):
        chosen.extend(strata[key][:quotas[key]])
    chosen.sort(key=lambda s: _order_key(seed, s))
    return {
        "chosen": chosen,
        "strata_sizes": {k: len(v) for k, v in strata.items()},
        "quotas": quotas,
        "eligible_total": total,
    }


def build_stub(paper: Dict, schema_version: str) -> Dict:
    stub = {k: paper.get(k) for k in STUB_FROM_HARVEST}
    stub["schema_version"] = schema_version
    stub["full_text_available"] = bool(paper.get("full_text_available"))
    if not stub.get("retrieved_from"):
        stub["retrieved_from"] = "europepmc_oa" if stub["full_text_available"] else "pubmed"
    stub["extraction_method"] = "llm"
    stub["review_status"] = "unreviewed"
    stub["evidence_class"] = "causal_interventional"
    return stub


def package(chosen: Sequence[str], elig: Dict[str, Dict], papers: Dict[str, Dict],
            *, n_batches: int, schema_version: str) -> List[List[Dict]]:
    """Round-robin papers into batches so each gets a similar reading load."""
    rows: List[Dict] = []
    for sid in chosen:
        info = elig[sid]
        paper = papers.get(sid) or {"study_id": sid}
        # One entry per distinct passage; the triples that pointed at it are
        # listed as leads, because a passage may yield several records or none.
        passages: Dict[str, Dict] = {}
        for u in sorted(info["units"], key=lambda x: (-x.get("score", 0), x.get("locator") or "")):
            loc = u.get("locator")
            p = passages.get(loc)
            if p is None:
                p = passages[loc] = {"source_locator": loc, "section": u.get("section_label"),
                                     "passage_text": u.get("passage_text"), "leads": []}
            p["leads"].append({
                "candidate_population_id": u.get("candidate_population_id"),
                "candidate_method_id": u.get("candidate_method_id"),
                "candidate_outcome_id": u.get("candidate_outcome_id"),
                "evidence_terms": u.get("evidence_terms"),
                "unit_id": u.get("unit_id"),
            })
        rows.append({
            "study_id": sid,
            "covers_requested_populations": sorted(info["covers"]),
            "full_text_path": paper.get("full_text_path"),
            "record_stub": build_stub(paper, schema_version),
            "n_passages": len(passages),
            "passages": [passages[k] for k in sorted(passages)],
            "fields_to_complete": list(fields_to_complete()),
            "_reminder": ("A passage is a place worth reading, not a record. Many "
                          "passages yield none. Never propose a record whose "
                          "direction does not follow from its own support_quote."),
        })
    batches: List[List[Dict]] = [[] for _ in range(max(1, n_batches))]
    for i, row in enumerate(rows):
        batches[i % len(batches)].append(row)
    return [b for b in batches if b]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python3 -m syn4rch.curate.pilot",
                                 description=__doc__.split("\n\n")[0])
    ap.add_argument("--units", default=_default("data", "interim", "extraction_units.jsonl"))
    ap.add_argument("--candidates", default=_default("data", "interim", "candidates.jsonl"))
    ap.add_argument("--out-dir", default=_default("data", "interim", "pilot"))
    ap.add_argument("--populations", default="arc_agrp,arc_pomc",
                    help="comma-separated population ids, or 'all'")
    ap.add_argument("--papers", type=int, default=50)
    ap.add_argument("--batches", type=int, default=5)
    ap.add_argument("--seed", default="syn4rch-pilot-1")
    ap.add_argument("--method-class", default=None,
                    help="keep only units whose candidate method belongs to "
                         "this perturbation class (e.g. circuit_activity), "
                         "resolved through data/vocab/methods.json and "
                         "syn4rch.schema.rules.perturbation_class_from_method")
    ap.add_argument("--exclude-curated", default=None,
                    help="drop studies that already appear in this curated JSONL")
    args = ap.parse_args(argv)

    for path in (args.units, args.candidates):
        if not os.path.exists(path):
            sys.stderr.write("missing %s\n" % path)
            return 1
    # No try/except fallback here. A hardcoded version that silently stands in
    # for the real one is how the stub writer came to label records with a
    # version the schema package had already moved past.
    from ..schema.record import SCHEMA_VERSION

    populations = [p.strip() for p in args.populations.split(",") if p.strip()]
    units = load_units(args.units)
    papers = load_papers(args.candidates)
    units_in = len(units)
    if args.method_class:
        units = filter_units_by_class(units, args.method_class)
        if not units:
            sys.stderr.write("no units of class %s\n" % args.method_class)
            return 1
    elig = eligible(units, populations)
    if not elig:
        sys.stderr.write("no eligible papers for %s\n" % ", ".join(populations))
        return 1
    dropped_curated = 0
    if args.exclude_curated:
        if not os.path.exists(args.exclude_curated):
            sys.stderr.write("missing %s\n" % args.exclude_curated)
            return 1
        elig, dropped_curated = drop_curated_studies(elig, args.exclude_curated)
        if not elig:
            sys.stderr.write("every eligible paper is already curated\n")
            return 1
    sel = select(elig, n_papers=args.papers, populations=populations, seed=args.seed)
    batches = package(sel["chosen"], elig, papers,
                      n_batches=args.batches, schema_version=SCHEMA_VERSION)

    os.makedirs(args.out_dir, exist_ok=True)
    written = []
    for i, batch in enumerate(batches, 1):
        path = os.path.join(args.out_dir, "batch_%02d.jsonl" % i)
        with open(path, "w", encoding="utf-8") as fh:
            for row in batch:
                fh.write(json.dumps(row, sort_keys=True) + "\n")
        written.append({"batch": i, "path": os.path.relpath(path, _repo_root()),
                        "papers": len(batch),
                        "passages": sum(r["n_passages"] for r in batch),
                        "chars": sum(len(p["passage_text"] or "")
                                     for r in batch for p in r["passages"])})

    manifest = {
        "built_at": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "selection_rule": "proportional stratified; order = sha256(seed|study_id)",
        "seed": args.seed,
        "requested_populations": populations,
        "method_class_filter": args.method_class,
        "units_before_method_filter": units_in,
        "units_after_method_filter": len(units),
        "excluded_curated_from": (
            os.path.relpath(args.exclude_curated, _repo_root())
            if args.exclude_curated else None
        ),
        "studies_dropped_as_already_curated": dropped_curated,
        "result_sections": list(RESULT_SECTIONS),
        "eligible_papers": sel["eligible_total"],
        "strata_sizes": sel["strata_sizes"],
        "quotas": sel["quotas"],
        "papers_selected": len(sel["chosen"]),
        "selected_study_ids": sel["chosen"],
        "batches": written,
        "extraction_method": "llm",
        "review_status": "unreviewed",
        "note": ("Counts of passages are reading load, not an expected record "
                 "count. Yield is measured by the gate, not predicted here."),
    }
    mpath = os.path.join(args.out_dir, "manifest.json")
    with open(mpath, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=1)

    print("eligible %d papers; strata %s" % (sel["eligible_total"], sel["strata_sizes"]))
    print("selected %d papers into %d batches (seed %r)"
          % (len(sel["chosen"]), len(written), args.seed))
    for b in written:
        print("  %s  %2d papers  %3d passages  %6d chars" %
              (b["path"], b["papers"], b["passages"], b["chars"]))
    print("manifest %s" % os.path.relpath(mpath, _repo_root()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
