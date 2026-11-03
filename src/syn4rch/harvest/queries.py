"""Build the query matrix programmatically from the controlled vocabularies.

Four tiers, all generated from ``data/vocab/*.json`` so that extending a
vocabulary extends the search:

* **T1 population-wide** -- (population synonyms + Cre lines) x (all perturbation
  modality terms) x (all physiological readout terms). One query per population.
* **T2 population x readout cluster** -- the same, but with the readout block
  narrowed to one cluster (metabolic, autonomic, endocrine, ...). This is where
  most of the recall comes from, because a single huge OR block gets dominated by
  whatever PubMed's relevance ranking likes.
* **T3 region-only** -- (hypothalamic region synonyms, no population terms) x
  modality x readout cluster. This is the net for populations that are NOT in the
  vocabulary; hits here that no T1/T2 query returned are exactly the
  vocabulary-extension candidates.
* **T4 modality-only** -- modality x hypothalamus, with no readout constraint,
  plus the Tier B aging block. Deliberately low-precision.

Every query object carries the tier, the population/cluster it came from and the
engine it is meant for, and ``run.py`` writes each issued query with its result
count to ``docs/QUERIES.md``.
"""

from __future__ import annotations

import re
from typing import Dict, Iterable, List, Optional

from ..vocab import load

# --------------------------------------------------------------------------- #
# readout clusters (hand-assigned grouping of outcome ids; every outcome id in
# outcomes.json must appear exactly once -- checked by check_clusters())
# --------------------------------------------------------------------------- #

OUTCOME_CLUSTERS: Dict[str, List[str]] = {
    "glucose_metabolism": [
        "blood_glucose", "glucose_tolerance", "insulin_sensitivity",
        "plasma_insulin", "plasma_glucagon", "hepatic_glucose_production",
    ],
    "energy_balance": [
        "food_intake", "energy_expenditure", "rer", "body_weight",
        "fat_mass", "lean_mass", "wat_lipolysis", "hepatic_lipid",
    ],
    "thermoregulation": [
        "bat_thermogenesis", "core_temperature",
    ],
    "cardiovascular_autonomic": [
        "heart_rate", "blood_pressure", "sympathetic_nerve_activity",
        "parasympathetic_tone",
    ],
    "neuroendocrine": [
        "corticosterone", "gh_secretion", "igf1", "thyroid_hormone",
    ],
    "reproduction": [
        "gonadotropins", "estrous_cyclicity",
    ],
    "fluid_balance": [
        "water_intake", "osmolality_avp",
    ],
    "inflammation_immune": [
        "il6", "tnf", "inflammatory_markers_other", "immune_cell_composition",
    ],
    "behaviour_state": [
        "sleep", "locomotor_activity", "gut_motility",
    ],
    "musculoskeletal": [
        "bone_density",
    ],
    "aging_tier_b": [
        "median_lifespan", "maximum_lifespan", "frailty_index", "grip_strength",
        "cognitive_decline", "tissue_aging_markers", "senescence_markers",
    ],
}

# Terms that are too generic or too ambiguous to put in a PubMed phrase block.
# "LH" and "lipid oxidation" style collisions are the reason this exists.
_TERM_STOPLIST = {
    "lh", "gh secretion", "fi", "cort", "hr", "map", "ffa", "me", "zi", "ph",
    "son", "hy", "mbh", "arc", "arh", "vmh", "vmn", "dmh", "pvh", "pvn", "lha",
    "poa", "mpa", "mpo", "mpn", "vmpo", "mnpo", "lpo", "scn", "sch", "tmn",
    "pmv", "avpv", "rp3v", "pefa", "paov", "pav", "flox", "flox/flox", "cko",
    "dta", "dtr", "acr", "arch", "nphr", "ox neurons", "t3", "t4", "free t3",
    "free t4", "p16", "p21", "il 6", "tnf", "crp", "ccl2", "cxcl1", "il 10",
    "gal neurons", "sna", "rsna", "msna", "rq", "rer", "ee", "vo2", "vco2",
    "bmd", "bv/tv", "ztd", "wsn", "tcore", "a12 dopamine neurons", "sasp",
    "ma", "fsh",
}

_SAFE = re.compile(r"^[0-9a-zA-Z][0-9a-zA-Z \-/\.\(\)'+,]*$")


def _usable(term: str) -> bool:
    t = (term or "").strip()
    if len(t) < 3:
        return False
    if t.lower() in _TERM_STOPLIST:
        return False
    if not _SAFE.match(t):
        return False
    # single short acronyms are dropped; multi-word phrases are kept
    if " " not in t and len(t) <= 4:
        return False
    return True


def _dedupe(terms: Iterable[str]) -> List[str]:
    seen = set()
    out = []
    for t in terms:
        key = t.lower().strip()
        if key and key not in seen:
            seen.add(key)
            out.append(t.strip())
    return out


# --------------------------------------------------------------------------- #
# term blocks
# --------------------------------------------------------------------------- #

def population_terms(pop_id: str, limit: int = 12) -> List[str]:
    pops = load("populations")
    entry = pops[pop_id]
    terms = [entry.get("label") or ""]
    terms += list(entry.get("synonyms") or [])
    terms += list(entry.get("cre_lines") or [])
    return [t for t in _dedupe(terms) if _usable(t)][:limit]


def modality_terms(per_entry: int = 6, limit: int = 140) -> List[str]:
    """Modality block, taking up to ``per_entry`` terms from EVERY method entry.

    Round-robin rather than a flat truncated list: a flat list cut at N silently
    drops whole method families (ablation, lesion, knockout, pharmacology all sit
    at the end of methods.json) and the query then only finds opto/chemo papers.
    """
    meths = load("methods")
    buckets: List[List[str]] = []
    for entry in meths.values():
        cand: List[str] = []
        if entry.get("actuator"):
            cand.append(str(entry["actuator"]))
        cand.extend(str(s) for s in (entry.get("synonyms") or []))
        cand = [t for t in _dedupe(cand) if _usable(t)]
        buckets.append(cand[:per_entry])
    terms: List[str] = []
    for i in range(per_entry):
        for bucket in buckets:
            if i < len(bucket):
                terms.append(bucket[i])
    return _dedupe(terms)[:limit]


def outcome_terms(cluster: Optional[str] = None, *, per_entry: int = 8,
                  limit: int = 150) -> List[str]:
    """Readout block. Round-robin across outcome entries for the same reason as
    :func:`modality_terms`: a flat truncated list would cover only the first few
    outcomes in the file and silently search for nothing else."""
    outs = load("outcomes")
    ids = list(outs.keys()) if cluster is None else OUTCOME_CLUSTERS[cluster]
    if cluster is None:
        per_entry = min(per_entry, 3)
    buckets: List[List[str]] = []
    for oid in ids:
        entry = outs[oid]
        cand = [entry.get("label") or ""] + [str(s) for s in (entry.get("synonyms") or [])]
        cand = [t for t in _dedupe(cand) if _usable(t)]
        buckets.append(cand[:per_entry])
    terms: List[str] = []
    for i in range(per_entry):
        for bucket in buckets:
            if i < len(bucket):
                terms.append(bucket[i])
    return _dedupe(terms)[:limit]


def region_terms(limit: int = 30) -> List[str]:
    regions = load("regions")
    terms: List[str] = []
    for entry in regions.values():
        terms.append(entry.get("label") or "")
        terms.extend(str(s) for s in (entry.get("synonyms") or []))
    return [t for t in _dedupe(terms) if _usable(t)][:limit]


def check_clusters() -> None:
    """Every outcome id must be in exactly one cluster."""
    outs = set(load("outcomes").keys())
    seen: Dict[str, str] = {}
    for cluster, ids in OUTCOME_CLUSTERS.items():
        for oid in ids:
            if oid not in outs:
                raise ValueError("cluster %s references unknown outcome %r" % (cluster, oid))
            if oid in seen:
                raise ValueError("outcome %r in two clusters: %s, %s" % (oid, seen[oid], cluster))
            seen[oid] = cluster
    missing = sorted(outs - set(seen))
    if missing:
        raise ValueError("outcomes not assigned to any cluster: %s" % ", ".join(missing))


# --------------------------------------------------------------------------- #
# query string construction
# --------------------------------------------------------------------------- #

def pubmed_block(terms: List[str], field: str = "tiab") -> str:
    return "(" + " OR ".join('"%s"[%s]' % (t, field) for t in terms) + ")"


def epmc_block(terms: List[str]) -> str:
    return "(" + " OR ".join('"%s"' % t for t in terms) + ")"


ANIMAL_FILTER_PUBMED = '(mice[tiab] OR mouse[tiab] OR rat[tiab] OR rats[tiab] OR murine[tiab] OR "animals"[MeSH Terms])'


class Query(dict):
    """A query to be issued. A dict so it serialises straight to JSON."""


def _mk(qid: str, tier: str, engine: str, term: str, *, population: Optional[str],
        cluster: Optional[str]) -> Query:
    return Query(qid=qid, tier=tier, engine=engine, term=term,
                 population=population, cluster=cluster)


def build_matrix(*, include_epmc: bool = True,
                 clusters: Optional[List[str]] = None) -> List[Query]:
    """Return the full ordered query matrix. Order is the priority order:
    T1, then T2, then T3, then T4, so ``--max-queries`` truncates sensibly."""
    check_clusters()
    pops = load("populations")
    cluster_names = clusters or list(OUTCOME_CLUSTERS.keys())

    mod = modality_terms()
    mod_pm = pubmed_block(mod)
    mod_epmc = epmc_block(mod[:30])
    all_out_pm = pubmed_block(outcome_terms(None))
    regions_pm = pubmed_block(region_terms())
    regions_epmc = epmc_block(region_terms(limit=20))

    out: List[Query] = []

    # T1: one broad query per population
    for pid in pops:
        pterms = population_terms(pid)
        if not pterms:
            continue
        term = "%s AND %s AND %s" % (pubmed_block(pterms), mod_pm, all_out_pm)
        out.append(_mk("T1-%s" % pid, "T1_population_wide", "pubmed", term,
                       population=pid, cluster=None))

    # T2: population x readout cluster
    for pid in pops:
        pterms = population_terms(pid)
        if not pterms:
            continue
        for cname in cluster_names:
            oterms = outcome_terms(cname)
            if not oterms:
                continue
            term = "%s AND %s AND %s" % (pubmed_block(pterms), mod_pm, pubmed_block(oterms))
            out.append(_mk("T2-%s-%s" % (pid, cname), "T2_population_x_cluster",
                           "pubmed", term, population=pid, cluster=cname))

    # T3: region-only x readout cluster (catches populations missing from vocab)
    for cname in cluster_names:
        oterms = outcome_terms(cname)
        term = "%s AND %s AND %s AND %s" % (regions_pm, mod_pm, pubmed_block(oterms),
                                            ANIMAL_FILTER_PUBMED)
        out.append(_mk("T3-region-%s" % cname, "T3_region_only", "pubmed", term,
                       population=None, cluster=cname))

    # T4: modality-only, no readout constraint
    out.append(_mk("T4-modality-hypothalamus", "T4_modality_only", "pubmed",
                   "%s AND %s AND %s" % (regions_pm, mod_pm, ANIMAL_FILTER_PUBMED),
                   population=None, cluster=None))
    out.append(_mk("T4-aging-hypothalamus", "T4_modality_only", "pubmed",
                   "%s AND %s AND %s" % (regions_pm, mod_pm,
                                         pubmed_block(outcome_terms("aging_tier_b"))),
                   population=None, cluster=None))

    if include_epmc:
        # Europe PMC: population-level sweep (catches preprints and non-MEDLINE)
        for pid in pops:
            pterms = population_terms(pid, limit=8)
            if not pterms:
                continue
            term = "%s AND %s" % (epmc_block(pterms), mod_epmc)
            out.append(_mk("E1-%s" % pid, "E1_population_epmc", "europepmc", term,
                           population=pid, cluster=None))
        out.append(_mk("E3-region", "E3_region_epmc", "europepmc",
                       "%s AND %s" % (regions_epmc, mod_epmc),
                       population=None, cluster=None))

    return out


if __name__ == "__main__":  # pragma: no cover
    qs = build_matrix()
    print("%d queries" % len(qs))
    from collections import Counter
    for tier, n in Counter(q["tier"] for q in qs).items():
        print("  %-26s %d" % (tier, n))
    print(qs[0]["term"][:400])
