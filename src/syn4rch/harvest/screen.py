"""Transparent rule-based title/abstract screener.

Produces, for every record, a decision in {include, uncertain, exclude}, an
integer score, and an explicit list of REASONS. Nothing is ever silently
dropped: a record the rules cannot place lands in ``uncertain``, which is a real
bucket meant to be used liberally, not a failure mode.

What it can and cannot do
-------------------------
It answers: does this title/abstract mention (a) a perturbation modality,
(b) a hypothalamic site, and (c) a physiological readout -- and is it a review.
It cannot tell you whether a usable (population x perturbation x outcome)
record can be EXTRACTED from the paper. That needs full text. Screening counts
are therefore an upper bound on candidate papers.

Matching rules
--------------
* Long terms are matched case-insensitively on a normalised copy of the text,
  with word boundaries.
* Short all-capital acronyms (<=5 chars, e.g. PVN, ARC, ChR2, DTA) are matched
  CASE-SENSITIVELY on the original text, because "arc" is an English word and
  "LH" is both a nucleus and a hormone. Acronyms on an explicit ambiguity
  blacklist are not matched at all.
* All terms come from the vocabularies; nothing is hard-coded except the
  publication-type exclusion list, the animal/human cue lists and the ambiguity
  blacklist, each of which is spelled out below.
"""

from __future__ import annotations

import re
from typing import Dict, Iterable, List, Optional, Set, Tuple

from ..vocab import load
from ..vocab.core import meta

# Publication types that cannot contain a primary perturbation experiment.
EXCLUDED_PUBTYPES = {
    "review", "systematic review", "meta-analysis", "comment", "editorial",
    "news", "published erratum", "retraction of publication", "retracted publication",
    "biography", "historical article", "newspaper article", "patient education handout",
    "address", "congress", "lecture", "interview", "autobiography", "bibliography",
    "video-audio media", "practice guideline", "guideline", "scoping review",
    "narrative review", "preprint",
}

# Acronyms that are too ambiguous to screen on even case-sensitively.
AMBIGUOUS_ACRONYMS = {
    "LH", "ARC", "ME", "ZI", "PH", "MAP", "HR", "EE", "FI", "RQ", "RER",
    "CORT", "T3", "T4", "SON", "HY", "ACR", "DBS", "MSNA", "C21", "DCZ",
    "AHA", "PV", "PaV", "GH", "FSH", "CRP", "SNA",
}

ANIMAL_CUES = (
    "mice", "mouse", "rat", "rats", "murine", "mus musculus", "rattus",
    "c57bl", "rodent", "rodents", "animals", "zebrafish", "macaque", "sheep",
    "ewes", "pigs", "porcine", "ovine", "hamster", "gerbil", "primate",
)
HUMAN_CUES = (
    "patients", "patient", "participants", "volunteers", "human subjects",
    "men and women", "cohort", "cross-sectional", "observational study",
    "questionnaire", "biobank", "nhanes", "clinical trial", "subjects were recruited",
)
IN_VITRO_CUES = (
    "cell line", "hek293", "cultured neurons", "primary culture", "organoid",
    "in vitro only", "cell-free",
)

_WS = re.compile(r"\s+")


def _norm_text(s: str) -> str:
    s = (s or "").replace("‐", "-").replace("‑", "-").replace("–", "-")
    s = s.replace("α", "alpha").replace("β", "beta").replace("κ", "kappa")
    return _WS.sub(" ", s.lower())


def _is_acronym(term: str) -> bool:
    t = term.strip()
    if " " in t or len(t) > 5:
        return False
    letters = [c for c in t if c.isalpha()]
    return bool(letters) and all(c.isupper() for c in letters)


def _boundary_pattern(term: str) -> str:
    esc = re.escape(term)
    esc = esc.replace(r"\ ", r"[\s\-]+")
    return r"(?<![0-9A-Za-z])" + esc + r"(?![0-9A-Za-z])"


class TermMatcher:
    """Compiled alternations for one category of vocabulary terms."""

    def __init__(self, name: str, terms_by_id: Dict[str, Iterable[str]]) -> None:
        self.name = name
        self._ci: List[Tuple[str, "re.Pattern[str]"]] = []
        self._cs: List[Tuple[str, "re.Pattern[str]"]] = []
        for vid, terms in terms_by_id.items():
            ci_terms, cs_terms = [], []
            for term in terms:
                t = (term or "").strip()
                if len(t) < 2:
                    continue
                if _is_acronym(t):
                    if t in AMBIGUOUS_ACRONYMS:
                        continue
                    cs_terms.append(t)
                else:
                    ci_terms.append(_norm_text(t))
            if ci_terms:
                self._ci.append((vid, re.compile(
                    "|".join(_boundary_pattern(t) for t in sorted(set(ci_terms), key=len, reverse=True))
                )))
            if cs_terms:
                self._cs.append((vid, re.compile(
                    "|".join(_boundary_pattern(t) for t in sorted(set(cs_terms), key=len, reverse=True))
                )))

    def match(self, raw_text: str, norm_text: str) -> List[str]:
        hits: Set[str] = set()
        for vid, rx in self._ci:
            if rx.search(norm_text):
                hits.add(vid)
        for vid, rx in self._cs:
            if rx.search(raw_text):
                hits.add(vid)
        return sorted(hits)


def _pop_terms() -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for pid, entry in load("populations").items():
        terms = [entry.get("label") or ""] + list(entry.get("synonyms") or [])
        terms += list(entry.get("cre_lines") or [])
        out[pid] = [t for t in terms if t]
    return out


def _region_terms() -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for rid, entry in load("regions").items():
        terms = [entry.get("label") or ""] + list(entry.get("synonyms") or [])
        out[rid] = [t for t in terms if t]
    return out


def _method_terms() -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for mid, entry in load("methods").items():
        terms = [entry.get("actuator") or ""] + list(entry.get("synonyms") or [])
        out[mid] = [t for t in terms if t]
    return out


def _outcome_terms() -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for oid, entry in load("outcomes").items():
        terms = [entry.get("label") or ""] + list(entry.get("synonyms") or [])
        out[oid] = [t for t in terms if t]
    return out


#: Matcher over the readouts `outcomes.json` declares out of scope (D-063), so
#: a passage mentioning one can be flagged rather than silently proposed under
#: whatever in-scope id it also happens to mention.
_OUT_OF_SCOPE: Optional[TermMatcher] = None


def out_of_scope_matcher() -> TermMatcher:
    """One matcher over every declared out-of-scope readout, by category."""
    global _OUT_OF_SCOPE
    if _OUT_OF_SCOPE is None:
        block = meta("outcomes").get("out_of_scope_readouts") or {}
        terms = {
            category: list(values)
            for category, values in block.items()
            if not category.startswith("_") and isinstance(values, list)
        }
        _OUT_OF_SCOPE = TermMatcher("out_of_scope", terms)
    return _OUT_OF_SCOPE


def out_of_scope_readouts(raw_text: str) -> List[str]:
    """Categories of out-of-scope readout this passage mentions.

    D-063 put neural activity out of scope and `tests/test_vocab_scope.py`
    guards the *vocabulary* -- no entry may list an out-of-scope term as a
    synonym. Nothing guarded the *matcher*, and it did not need to: a sentence
    reporting "EEG delta power during wakefulness increased" resolves to
    `wakefulness` through an incidental word, so an out-of-scope readout
    arrives under a legitimate id and no rule is broken on the way in. 224 of
    946 proposed sleep units carry one and 133 still resolve to an in-scope
    state (D-070).

    This does not suppress the id, because the passage may genuinely report
    both an in-scope state and an out-of-scope measure in one sentence, and
    deciding which the record is about is a reading judgement the vocabulary
    refuses to make for a curator (D-021). It names what is there, so the
    judgement is prompted rather than skipped.
    """
    return out_of_scope_matcher().match(raw_text, _norm_text(raw_text))


_MATCHERS: Optional[Dict[str, TermMatcher]] = None


def matchers() -> Dict[str, TermMatcher]:
    global _MATCHERS
    if _MATCHERS is None:
        _MATCHERS = {
            "populations": TermMatcher("populations", _pop_terms()),
            "regions": TermMatcher("regions", _region_terms()),
            "methods": TermMatcher("methods", _method_terms()),
            "outcomes": TermMatcher("outcomes", _outcome_terms()),
        }
    return _MATCHERS


# Score weights, stated here so a reader can recompute any score by hand.
WEIGHTS = {
    "modality": 3,
    "site": 3,
    "readout": 3,
    "named_population": 2,
    "animal_cue": 1,
    "multiple_readouts": 1,
    "excluded_pubtype": -8,
    "human_only": -3,
    "in_vitro_only": -2,
    "no_abstract": -2,
}


def screen_record(rec: Dict) -> Dict:
    """Screen one bibliographic record. Returns the screening block to attach."""
    title = rec.get("title") or ""
    abstract = rec.get("abstract") or ""
    mesh = " ; ".join(rec.get("mesh_terms") or [])
    raw = "\n".join([title, abstract, mesh])
    norm = _norm_text(raw)

    m = matchers()
    pops = m["populations"].match(raw, norm)
    regions = m["regions"].match(raw, norm)
    methods = m["methods"].match(raw, norm)
    outcomes = m["outcomes"].match(raw, norm)

    pubtypes = [str(p).strip().lower() for p in (rec.get("publication_types") or [])]
    excluded_types = sorted(set(pubtypes) & EXCLUDED_PUBTYPES)
    title_norm = _norm_text(title)
    if not excluded_types and re.search(r"\b(a|systematic|narrative|scoping|mini)?\s*review\b", title_norm):
        excluded_types = ["review (inferred from title)"]

    animal = [c for c in ANIMAL_CUES if c in norm]
    human = [c for c in HUMAN_CUES if c in norm]
    in_vitro = [c for c in IN_VITRO_CUES if c in norm]

    reasons: List[str] = []
    flags: List[str] = []
    score = 0

    has_modality = bool(methods)
    has_site = bool(pops or regions)
    has_readout = bool(outcomes)

    if has_modality:
        score += WEIGHTS["modality"]
        reasons.append("modality mentioned: %s" % ", ".join(methods[:4]))
    else:
        reasons.append("no perturbation modality term found")
    if has_site:
        score += WEIGHTS["site"]
        if pops:
            score += WEIGHTS["named_population"]
            reasons.append("population mentioned: %s" % ", ".join(pops[:4]))
        else:
            reasons.append("hypothalamic region mentioned but no vocabulary population: %s"
                           % ", ".join(regions[:4]))
    else:
        reasons.append("no hypothalamic region or population term found")
    if has_readout:
        score += WEIGHTS["readout"]
        reasons.append("physiological readout mentioned: %s" % ", ".join(outcomes[:5]))
        if len(outcomes) > 1:
            score += WEIGHTS["multiple_readouts"]
    else:
        reasons.append("no physiological readout term found")

    if animal:
        score += WEIGHTS["animal_cue"]
    if not abstract:
        score += WEIGHTS["no_abstract"]
        reasons.append("abstract not available at screening time")
        flags.append("abstract_missing")
    if human and not animal:
        score += WEIGHTS["human_only"]
        flags.append("human_observational")
        reasons.append("human study cues (%s) with no animal cue: flagged, not excluded"
                       % ", ".join(human[:3]))
    if in_vitro and not animal:
        score += WEIGHTS["in_vitro_only"]
        flags.append("possibly_in_vitro_only")
    if excluded_types:
        score += WEIGHTS["excluded_pubtype"]

    # --- decision ---------------------------------------------------------
    n_criteria = sum([has_modality, has_site, has_readout])
    if excluded_types:
        decision = "exclude"
        reasons.insert(0, "excluded publication type: %s" % ", ".join(excluded_types))
        primary = "publication_type"
    elif n_criteria == 3:
        if not abstract:
            decision = "uncertain"
            primary = "all_three_but_no_abstract"
        else:
            decision = "include"
            primary = "modality+site+readout"
    elif n_criteria == 2:
        decision = "uncertain"
        missing = [n for n, ok in (("modality", has_modality), ("site", has_site),
                                   ("readout", has_readout)) if not ok]
        primary = "missing_" + "_".join(missing)
    else:
        decision = "exclude"
        missing = [n for n, ok in (("modality", has_modality), ("site", has_site),
                                   ("readout", has_readout)) if not ok]
        primary = "missing_" + "_".join(missing)
        reasons.insert(0, "fewer than two of the three required criteria present")

    return {
        "decision": decision,
        "score": score,
        "primary_reason": primary,
        "reasons": reasons,
        "flags": flags,
        "matched_populations": pops,
        "matched_regions": regions,
        "matched_methods": methods,
        "matched_outcomes": outcomes,
        "n_criteria_met": n_criteria,
        "screener_version": "0.1.0",
        "screened_on": "title+abstract+mesh",
    }


def screen_all(records: Iterable[Dict]) -> List[Dict]:
    out = []
    for rec in records:
        rec = dict(rec)
        rec["screen"] = screen_record(rec)
        out.append(rec)
    return out
