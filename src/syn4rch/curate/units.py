"""Extraction unit candidates: places worth reading. Stage 5 of the protocol.

**An extraction unit candidate is not a record and must never be stored as
one.** It asserts nothing about biology. It says only this: inside one passage
of one paper, a population term, a perturbation-method term and an outcome
term from the controlled vocabularies all occur. Co-occurrence is not a
finding (``docs/CURATION_PROTOCOL.md`` 2). The passage may be an introduction
sentence citing someone else's work, a methods sentence with no result in it,
or a discussion sentence speculating about a mechanism, and the reader at
stage 6 must throw all three away. A low units-to-records yield is the
protocol working, not a bug.

The type is called :class:`ExtractionUnitCandidate` for that reason, and it is
built so that it cannot be mistaken for a
:class:`~syn4rch.schema.record.TierARecord`:

* it has no ``direction``, no ``valence``, no ``effect_size_value`` and no
  ``confidence`` -- every field a record needs in order to make a claim is
  absent, because every one of them requires reading;
* it has no ``support_quote``. It carries ``passage_text``, the whole passage
  verbatim. A quote is something a reader *selects* out of that;
* :meth:`ExtractionUnitCandidate.to_dict` stamps every serialised row with
  ``"_kind": "extraction_unit_candidate"`` and ``"_not_a_record": true``, so a
  row that leaks into a records file is visible to ``grep``;
* ``candidate_population_id`` and friends are named ``candidate_`` because a
  vocabulary term appearing in a passage is not yet that passage's subject.

Matching reuses :class:`syn4rch.harvest.screen.TermMatcher` unchanged -- the
same boundary-aware patterns, the same case-sensitivity rule for short
acronyms, the same ``AMBIGUOUS_ACRONYMS`` stoplist. There is deliberately no
second matcher in this project: two matchers would eventually disagree, and
then the screening counts and the curation counts would describe different
corpora.

Localisation, and what it costs
-------------------------------
All three matches must fall inside the *same passage*. A Results paragraph
reading "stimulation increased food intake", whose population was named two
paragraphs earlier, produces no unit. That is a recall cost, accepted
knowingly: a unit whose population comes from elsewhere in the paper would
invite a reader to write a ``support_quote`` that does not contain its own
subject, which is how a plausible record with no basis gets made.
:func:`paper_term_summary` reports what a paper mentions overall, so the
recall gap can be measured instead of guessed.

Python 3.10 compatible, standard library only.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# Imported, never reimplemented: the screener's matchers are the project's
# single definition of "this vocabulary term occurs in this text".
from ..harvest.screen import TermMatcher, matchers
from ..harvest.screen import _norm_text as normalize_for_matching
from .jats import ParsedPaper
from .passages import (
    KIND_FIGURE_CAPTION,
    KIND_PARAGRAPH,
    KIND_TABLE_CAPTION,
    KIND_TITLE,
    Passage,
    split_passages,
)

__all__ = [
    "ExtractionUnitCandidate",
    "PassageMatches",
    "SCORE_WEIGHTS",
    "UNIT_ID_LENGTH",
    "match_passage",
    "score_unit",
    "units_from_paper",
    "units_from_papers",
    "paper_term_summary",
    "matched_terms",
]

#: Length of the hex ``unit_id``. Units are addresses, not dataset rows, so a
#: shorter digest than ``record_id``'s would still be fine; it matches for
#: readability.
UNIT_ID_LENGTH = 16

_HASH_SEP = "\x1f"
_DIGIT = re.compile(r"\d")

#: Lexical cues that a passage reports a *change* rather than a procedure.
#: Used for RANKING ONLY. Nothing here is ever written into a record: the
#: direction of an effect is a reading judgement (SCHEMA.md 2.2), and a word
#: list cannot make it. A passage containing "unchanged" and a passage
#: containing "increased" score identically, deliberately -- both are worth
#: reading, and ``no_change`` results must not be ranked below positive ones
#: (SCHEMA.md 3.5).
CHANGE_CUES = (
    "increase",
    "increased",
    "increases",
    "decrease",
    "decreased",
    "decreases",
    "reduced",
    "reduction",
    "elevated",
    "lowered",
    "higher",
    "lower",
    "greater",
    "attenuated",
    "abolished",
    "unchanged",
    "no change",
    "no significant",
    "did not alter",
    "did not change",
    "comparable",
    "indistinguishable",
    "impaired",
    "improved",
    "suppressed",
    "enhanced",
    "blunted",
    "rescued",
)

#: Ranking weights. Every score in ``extraction_units.jsonl`` is recomputable
#: by hand from these numbers plus the ``score_breakdown`` the row carries.
#:
#: The ordering these produce is a reading *queue*, nothing more. A unit with
#: score 2 is as admissible a source for a record as one with score 11; the
#: score only says which passage a reader should open first.
SCORE_WEIGHTS = {
    # Where the passage sits. Effects are reported in results and in float
    # captions; methods sections describe procedures, and introductions and
    # discussions describe other people's findings and one's own
    # interpretations.
    "section_results": 4,
    "section_caption": 3,
    "section_abstract": 2,
    "section_methods": -2,
    "section_introduction": -2,
    "section_discussion": -1,
    "section_title": -1,
    # What the passage contains.
    "has_digit": 1,
    "change_cue": 2,
    # Coherence between the matched terms.
    "region_matched": 1,
    "region_is_population_nucleus": 1,
    # Passage shape. A passage too short to contain a subject, a manipulation
    # and a readout is usually a heading or a fragment.
    "very_short_passage": -2,
    # A passage matching several populations or several outcomes is more
    # likely to be a summary sentence than a single observation.
    "ambiguous_population_set": -1,
    "ambiguous_outcome_set": -1,
}

#: A passage shorter than this (characters) earns ``very_short_passage``.
SHORT_PASSAGE_CHARS = 120


# --------------------------------------------------------------------------- #
# matching one passage
# --------------------------------------------------------------------------- #


def matched_terms(matcher: TermMatcher, raw_text: str, norm_text: str) -> Dict[str, List[str]]:
    """Surface strings that triggered each id in ``matcher``, for the reader.

    Reads the patterns :class:`~syn4rch.harvest.screen.TermMatcher` already
    compiled rather than compiling anything new, so a term matches here if and
    only if it matches in screening.

    The strings returned are **match evidence**, not field values. A
    case-insensitive hit is reported as it appears in the normalised copy of
    the text (lower-cased, dashes folded), so it is not necessarily verbatim.
    ``population_raw`` must be copied out of ``passage_text`` by the reader --
    SCHEMA.md 3.2 says "exactly as the paper names it", and this is not that.
    """
    out: Dict[str, List[str]] = {}
    for vid, pattern in getattr(matcher, "_ci", ()):
        hits = [m.group(0) for m in pattern.finditer(norm_text)]
        if hits:
            out.setdefault(vid, []).extend(hits)
    for vid, pattern in getattr(matcher, "_cs", ()):
        hits = [m.group(0) for m in pattern.finditer(raw_text)]
        if hits:
            out.setdefault(vid, []).extend(hits)
    return {vid: sorted(dict.fromkeys(hits)) for vid, hits in sorted(out.items())}


@dataclass(frozen=True)
class PassageMatches:
    """Vocabulary terms found inside one passage. Still not a finding."""

    passage: Passage
    populations: Tuple[str, ...] = ()
    regions: Tuple[str, ...] = ()
    methods: Tuple[str, ...] = ()
    outcomes: Tuple[str, ...] = ()
    terms: Dict[str, Dict[str, List[str]]] = field(default_factory=dict)

    @property
    def is_candidate(self) -> bool:
        """True when a population, a method and an outcome all co-occur here."""
        return bool(self.populations and self.methods and self.outcomes)


def match_passage(passage: Passage, matcher_set: Optional[Dict[str, TermMatcher]] = None) -> PassageMatches:
    """Run the screener's matchers over one passage."""
    ms = matcher_set if matcher_set is not None else matchers()
    raw = passage.text
    norm = normalize_for_matching(raw)
    kinds = ("populations", "regions", "methods", "outcomes")
    ids = {kind: tuple(ms[kind].match(raw, norm)) for kind in kinds}
    terms = {
        kind: matched_terms(ms[kind], raw, norm) for kind in kinds if ids[kind]
    }
    return PassageMatches(
        passage=passage,
        populations=ids["populations"],
        regions=ids["regions"],
        methods=ids["methods"],
        outcomes=ids["outcomes"],
        terms=terms,
    )


# --------------------------------------------------------------------------- #
# the candidate
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ExtractionUnitCandidate:
    """One (population, method, outcome) co-occurrence inside one passage.

    **Not a record.** See the module docstring. Nothing in this class states a
    direction, an effect size, a valence or a confidence, because none of
    those can be derived from co-occurrence, and a field that exists gets
    filled in.

    The fields a reader at stage 6 needs are ``passage_text`` (the source,
    verbatim) and ``locator`` (what a resulting record's ``source_locator``
    must say, unmodified -- ``docs/CURATION_PROTOCOL.md`` 3.3). The
    ``candidate_*`` ids are hints about which vocabulary entries might apply.
    """

    unit_id: str
    study_id: Optional[str]
    pmcid: Optional[str]
    pmid: Optional[str]
    doi: Optional[str]
    locator: str
    section: str
    section_label: str
    passage_kind: str
    passage_text: str
    candidate_population_id: str
    candidate_method_id: str
    candidate_outcome_id: str
    candidate_region_ids: Tuple[str, ...] = ()
    evidence_terms: Dict[str, List[str]] = field(default_factory=dict)
    score: int = 0
    score_breakdown: Dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, object]:
        """JSON-serialisable row, stamped as what it is.

        ``_kind`` and ``_not_a_record`` are written on every row so that a
        unit which ends up in the wrong file is findable with ``grep`` instead
        of being indistinguishable from curated evidence.
        """
        return {
            "_kind": "extraction_unit_candidate",
            "_not_a_record": True,
            "unit_id": self.unit_id,
            "study_id": self.study_id,
            "pmcid": self.pmcid,
            "pmid": self.pmid,
            "doi": self.doi,
            "locator": self.locator,
            "section": self.section,
            "section_label": self.section_label,
            "passage_kind": self.passage_kind,
            "passage_text": self.passage_text,
            "candidate_population_id": self.candidate_population_id,
            "candidate_method_id": self.candidate_method_id,
            "candidate_outcome_id": self.candidate_outcome_id,
            "candidate_region_ids": list(self.candidate_region_ids),
            "evidence_terms": {k: list(v) for k, v in self.evidence_terms.items()},
            "score": self.score,
            "score_breakdown": dict(self.score_breakdown),
        }

    def sort_key(self) -> Tuple:
        """Total order: score descending, then identity, for determinism."""
        return (
            -self.score,
            str(self.study_id or ""),
            self.locator,
            self.candidate_population_id,
            self.candidate_method_id,
            self.candidate_outcome_id,
        )


def _unit_id(
    study_id: Optional[str],
    locator: str,
    population_id: str,
    method_id: str,
    outcome_id: str,
) -> str:
    payload = _HASH_SEP.join(
        [str(study_id or ""), locator, population_id, method_id, outcome_id]
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:UNIT_ID_LENGTH]


# --------------------------------------------------------------------------- #
# scoring
# --------------------------------------------------------------------------- #

_SECTION_WEIGHT_KEY = {
    "Results": "section_results",
    "Abstract": "section_abstract",
    "Methods": "section_methods",
    "Introduction": "section_introduction",
    "Discussion": "section_discussion",
}


def score_unit(
    passage: Passage,
    population_id: str,
    region_ids: Sequence[str],
    population_nucleus: Optional[str],
    n_populations: int,
    n_outcomes: int,
) -> Tuple[int, Dict[str, int]]:
    """Score one candidate. Returns ``(score, breakdown)``.

    The rule is the sum of the applicable :data:`SCORE_WEIGHTS` entries, and
    the breakdown names every one that applied, so any score in
    ``extraction_units.jsonl`` is checkable by hand. The score ranks a reading
    queue; it is not a measure of evidence, and no record field is derived
    from it.
    """
    breakdown: Dict[str, int] = {}

    def add(key: str) -> None:
        breakdown[key] = SCORE_WEIGHTS[key]

    if passage.kind in (KIND_FIGURE_CAPTION, KIND_TABLE_CAPTION):
        add("section_caption")
    elif passage.kind == KIND_TITLE:
        add("section_title")
    else:
        key = _SECTION_WEIGHT_KEY.get(passage.section_label)
        if key:
            add(key)

    if _DIGIT.search(passage.text):
        add("has_digit")
    lowered = passage.text.lower()
    if any(cue in lowered for cue in CHANGE_CUES):
        add("change_cue")

    if region_ids:
        add("region_matched")
        if population_nucleus and population_nucleus in region_ids:
            add("region_is_population_nucleus")

    if len(passage.text) < SHORT_PASSAGE_CHARS:
        add("very_short_passage")
    if n_populations > 2:
        add("ambiguous_population_set")
    if n_outcomes > 3:
        add("ambiguous_outcome_set")

    return sum(breakdown.values()), breakdown


# --------------------------------------------------------------------------- #
# building
# --------------------------------------------------------------------------- #


def _population_nucleus(population_id: str) -> Optional[str]:
    try:
        from ..vocab import load

        entry = load("populations").get(population_id) or {}
    except Exception:
        return None
    nucleus = entry.get("nucleus")
    return str(nucleus) if nucleus else None


def units_from_paper(
    paper: ParsedPaper,
    study_id: Optional[str] = None,
    passages: Optional[Sequence[Passage]] = None,
    matcher_set: Optional[Dict[str, TermMatcher]] = None,
    max_per_passage: int = 0,
) -> List[ExtractionUnitCandidate]:
    """Every localized (population, method, outcome) co-occurrence in ``paper``.

    Returns candidates sorted by :meth:`ExtractionUnitCandidate.sort_key`
    (score descending, then a stable identity tie-break). ``max_per_passage``
    of ``0`` means no cap; a positive value keeps the highest-scoring N per
    passage, which is a *reading budget*, so the caller is responsible for
    reporting how many were dropped (``syn4rch.curate.run`` does).
    """
    if passages is None:
        passages = split_passages(paper)
    resolved_study = study_id or _default_study_id(paper)
    out: List[ExtractionUnitCandidate] = []
    for passage in passages:
        matches = match_passage(passage, matcher_set=matcher_set)
        if not matches.is_candidate:
            continue
        per_passage: List[ExtractionUnitCandidate] = []
        for population_id in matches.populations:
            nucleus = _population_nucleus(population_id)
            score, breakdown = score_unit(
                passage,
                population_id,
                matches.regions,
                nucleus,
                len(matches.populations),
                len(matches.outcomes),
            )
            for method_id in matches.methods:
                for outcome_id in matches.outcomes:
                    evidence = {
                        "populations": matches.terms.get("populations", {}).get(
                            population_id, []
                        ),
                        "methods": matches.terms.get("methods", {}).get(method_id, []),
                        "outcomes": matches.terms.get("outcomes", {}).get(
                            outcome_id, []
                        ),
                        "regions": sorted(
                            {
                                term
                                for rid in matches.regions
                                for term in matches.terms.get("regions", {}).get(rid, [])
                            }
                        ),
                    }
                    per_passage.append(
                        ExtractionUnitCandidate(
                            unit_id=_unit_id(
                                resolved_study,
                                passage.locator,
                                population_id,
                                method_id,
                                outcome_id,
                            ),
                            study_id=resolved_study,
                            pmcid=paper.pmcid,
                            pmid=paper.pmid,
                            doi=paper.doi,
                            locator=passage.locator,
                            section=passage.section,
                            section_label=passage.section_label,
                            passage_kind=passage.kind,
                            passage_text=passage.text,
                            candidate_population_id=population_id,
                            candidate_method_id=method_id,
                            candidate_outcome_id=outcome_id,
                            candidate_region_ids=tuple(matches.regions),
                            evidence_terms=evidence,
                            score=score,
                            score_breakdown=breakdown,
                        )
                    )
        per_passage.sort(key=lambda unit: unit.sort_key())
        if max_per_passage > 0:
            per_passage = per_passage[:max_per_passage]
        out.extend(per_passage)
    out.sort(key=lambda unit: unit.sort_key())
    return out


def units_from_papers(
    papers: Iterable[Tuple[Optional[str], ParsedPaper]],
    max_per_passage: int = 0,
) -> List[ExtractionUnitCandidate]:
    """Units for several ``(study_id, paper)`` pairs, globally sorted."""
    shared = matchers()
    out: List[ExtractionUnitCandidate] = []
    for study_id, paper in papers:
        out.extend(
            units_from_paper(
                paper,
                study_id=study_id,
                matcher_set=shared,
                max_per_passage=max_per_passage,
            )
        )
    out.sort(key=lambda unit: unit.sort_key())
    return out


def _default_study_id(paper: ParsedPaper) -> Optional[str]:
    """``study_id`` implied by a paper's own ids (SCHEMA.md 3.1 order)."""
    if paper.doi:
        return "doi:%s" % paper.doi
    if paper.pmid:
        return "pmid:%s" % paper.pmid
    if paper.pmcid:
        return "pmcid:%s" % paper.pmcid
    return None


def paper_term_summary(
    paper: ParsedPaper,
    passages: Optional[Sequence[Passage]] = None,
    matcher_set: Optional[Dict[str, TermMatcher]] = None,
) -> Dict[str, List[str]]:
    """What the paper mentions anywhere, regardless of co-occurrence.

    Exists so the recall cost of localisation can be *measured*: a paper that
    mentions a population, a method and an outcome somewhere but yields no
    unit is a paper whose claim is spread across paragraphs. Counting those is
    how the gap stops being a guess. These lists are never used to build a
    unit.
    """
    if passages is None:
        passages = split_passages(paper)
    ms = matcher_set if matcher_set is not None else matchers()
    acc: Dict[str, set] = {k: set() for k in ("populations", "regions", "methods", "outcomes")}
    for passage in passages:
        raw = passage.text
        norm = normalize_for_matching(raw)
        for kind in acc:
            acc[kind].update(ms[kind].match(raw, norm))
    return {kind: sorted(values) for kind, values in acc.items()}
