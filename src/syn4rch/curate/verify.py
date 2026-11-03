"""The verbatim gate. Stage 7 of docs/CURATION_PROTOCOL.md.

This is the one check in the project that a plausible-sounding fabrication
cannot talk its way past, and it works by being mechanical and small.

What it does
------------
``verify_support(record, paper)`` answers two questions about one proposed
record, and only these two:

1. Does ``record.support_quote`` occur as a substring of this paper's
   retrieved text, under exactly one normalisation?
2. Does ``record.source_locator`` name a passage that exists in this paper?

The normalisation is fixed and is the whole of it
(:func:`normalize_for_match`):

* runs of whitespace collapse to a single space;
* unicode dashes, quotation marks and non-breaking/typographic spaces fold to
  their ASCII equivalents;
* leading and trailing whitespace is stripped.

**Nothing else.** No case folding. No stemming or lemmatisation. No dropped,
added or reordered words. No punctuation stripping. No edit distance, token
overlap, embedding similarity or any other approximate match. ``"Food intake"``
does not verify against ``"food intake"``, and it is not supposed to: a reader
who could not be bothered to copy the capitalisation did not copy the
sentence.

The same normalisation is applied to both sides of every comparison, so the
check is symmetric and has no preferred direction.

What it refuses to do
---------------------
It never repairs anything. A failed quote is not re-searched with a looser
rule, not trimmed to its longest matching prefix, not re-quoted from the
nearest passage, and not written back to the record in any form. ``gate()``
returns rejections; it returns no corrected records, because a gate that
repairs its inputs is a formality (``docs/CURATION_PROTOCOL.md`` 4). Records
are treated as immutable throughout this module.

What it cannot catch -- stated plainly, because a reader of the dataset is
entitled to know (``docs/CURATION_PROTOCOL.md`` 5, DECISIONS D-015):

* It cannot tell whether the quote is **representative**. A real sentence
  cherry-picked to suit a conclusion passes.
* It cannot tell whether the stated ``direction``, ``valence``, effect size,
  sample size or any other field actually **follows from** the quote. A
  verbatim quote attached to the opposite direction passes.
* It cannot tell whether the quote is the paper's own claim. An introduction
  sentence describing somebody else's result is verbatim text in this paper.
* It cannot see any source that is not on disk: supplementary files, figure
  images, and the body of a paywalled paper are all invisible, so a locator
  naming one is rejected rather than checked.
* It cannot detect a quote assembled from two different papers' sentences
  beyond checking each record against the paper its own ``study_id`` names.

``extraction_method`` and ``review_status`` remain the only honest record of
whether a human made the judgements this gate cannot.

Python 3.10 compatible, standard library only.
"""

from __future__ import annotations

import sys
import unicodedata
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from ..schema.record import TierARecord
from ..schema.validate import ValidationError, hard_findings, validate
from .jats import ParsedPaper
from .passages import Passage, PassageIndex, build_index, parse_locator

__all__ = [
    "CODES",
    "FORBIDDEN_MODULES",
    "NORMALISATION",
    "forbidden_imports",
    "zs_coverage_is_complete",
    "Rejection",
    "SupportVerification",
    "VerificationFinding",
    "normalize_for_match",
    "quote_occurs_in",
    "verify_support",
    "gate",
    "index_papers",
    "find_paper",
]


# --------------------------------------------------------------------------- #
# the one normalisation
# --------------------------------------------------------------------------- #

#: Code points folded to ASCII ``-``. HYPHEN, NON-BREAKING HYPHEN, FIGURE
#: DASH, EN DASH, EM DASH, HORIZONTAL BAR, MINUS SIGN, the two small/fullwidth
#: forms and HYPHEN BULLET. Publisher XML and PDF extractors disagree about
#: these constantly and the disagreement carries no meaning.
#:
#: Every set below is written as escapes, not as the characters themselves, so
#: that the normalisation can be audited by reading this file and cannot be
#: changed by a copy-paste that silently swaps one dash for another.
_DASHES = "\u2010\u2011\u2012\u2013\u2014\u2015\u2212\ufe58\ufe63\uff0d\u2043"

#: Code points folded to ASCII ``'`` (the single-quote and prime family).
_SINGLE_QUOTES = "\u2018\u2019\u201a\u201b\u2032\u2035\u00b4\u02b9\u02bc\u2039\u203a"

#: Code points folded to ASCII ``"`` (the double-quote, double-prime and
#: guillemet family).
_DOUBLE_QUOTES = "\u201c\u201d\u201e\u201f\u2033\u2036\u00ab\u00bb\u02ba\u3003"

#: Code points folded to ASCII space: every Unicode space separator (general
#: category Zs -- which is where NO-BREAK SPACE U+00A0, NARROW NO-BREAK SPACE
#: U+202F and FIGURE SPACE U+2007 live), plus LINE SEPARATOR U+2028 and
#: PARAGRAPH SEPARATOR U+2029. Folding these is what "non-breaking spaces fold
#: to ASCII" means. Because whitespace runs then collapse, the fold is
#: unobservable except that it stops a non-breaking space from being a
#: different character than an ordinary one.
#:
#: Listed explicitly rather than derived by scanning all of Unicode, so that
#: the set is auditable and importing this module is cheap.
#: :func:`zs_coverage_is_complete` lets the test suite confirm the list is
#: still the whole of category Zs.
_SPACE_SEPARATORS = "\u0020\u00a0\u1680\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029"


def zs_coverage_is_complete() -> bool:
    """True when :data:`_SPACE_SEPARATORS` still covers all of category Zs.

    Called by the test suite. A future Unicode version adding a space
    separator would make this False, which is a finding about the
    normalisation, not about the test.
    """
    zs = {
        chr(cp) for cp in range(0x110000) if unicodedata.category(chr(cp)) == "Zs"
    }
    return zs <= set(_SPACE_SEPARATORS)


_TRANSLATION = {}
for _ch in _DASHES:
    _TRANSLATION[ord(_ch)] = "-"
for _ch in _SINGLE_QUOTES:
    _TRANSLATION[ord(_ch)] = "'"
for _ch in _DOUBLE_QUOTES:
    _TRANSLATION[ord(_ch)] = '"'
for _ch in _SPACE_SEPARATORS:
    _TRANSLATION[ord(_ch)] = " "
del _ch

#: One-line description of the normalisation, for reports.
NORMALISATION = (
    "collapse whitespace runs to one space; fold unicode dashes, quotation "
    "marks and non-breaking spaces to ASCII; strip the ends. Nothing else: "
    "no case folding, no stemming, no dropped words, no fuzzy matching."
)


def normalize_for_match(text: Optional[str]) -> str:
    """The gate's only normalisation. See :data:`NORMALISATION`.

    Applied identically to the quote and to the source text. Case is
    preserved. Word content is preserved. Punctuation other than the dash and
    quote families is preserved exactly.
    """
    if text is None:
        return ""
    folded = str(text).translate(_TRANSLATION)
    return " ".join(folded.split())


def quote_occurs_in(quote: str, text: str) -> bool:
    """True if ``quote`` occurs verbatim in ``text`` under the normalisation.

    A plain substring test on two normalised strings. An empty quote returns
    ``False``: the empty string is a substring of everything, which would make
    a missing quote verify.
    """
    needle = normalize_for_match(quote)
    if not needle:
        return False
    return needle in normalize_for_match(text)


# --------------------------------------------------------------------------- #
# findings
# --------------------------------------------------------------------------- #

#: Every code this module can emit. ``code -> (blocking, meaning)``.
#:
#: ``blocking`` is what decides admission. A blocking code means the record
#: does not enter ``data/curated/`` -- there is no override, no threshold and
#: no "accept with a note" path.
CODES: Dict[str, Tuple[bool, str]] = {
    "C001_TEXT_NOT_RETRIEVED": (
        True,
        "No retrieved text is held for this record's study, so there is "
        "nothing to check the quote against (SCHEMA.md rule 2.5).",
    ),
    "C002_MISSING_SUPPORT_QUOTE": (
        True,
        "support_quote is absent or whitespace-only (SCHEMA.md rule 2.1).",
    ),
    "C003_MISSING_SOURCE_LOCATOR": (
        True,
        "source_locator is absent or whitespace-only (SCHEMA.md rule 2.3).",
    ),
    "C004_QUOTE_NOT_VERBATIM": (
        True,
        "support_quote does not occur in the paper's retrieved text under the "
        "documented normalisation (SCHEMA.md rule 2.1).",
    ),
    "C005_LOCATOR_UNPARSEABLE": (
        True,
        "source_locator is not in a form this project can resolve to a "
        "passage (SCHEMA.md rule 2.3).",
    ),
    "C006_LOCATOR_NOT_IN_PAPER": (
        True,
        "source_locator parses but names a passage this paper does not have "
        "(SCHEMA.md rule 2.3).",
    ),
    "C007_SCHEMA_HARD": (
        True,
        "syn4rch.schema.validate reported a hard finding; the record is "
        "inadmissible regardless of its quote.",
    ),
    "C008_UNREADABLE_PROPOSAL": (
        True,
        "The proposed record could not be read as a Tier A record at all.",
    ),
    "W101_QUOTE_NOT_IN_LOCATED_PASSAGE": (
        False,
        "The quote is verbatim somewhere in this paper but not in the passage "
        "the locator names. Reported, not blocking: paragraph numbering is "
        "parser-dependent, and CURATION_PROTOCOL.md 4 defines the gate as "
        "'substring of the paper's retrieved text' plus 'the locator names a "
        "passage that exists'. Blocking by default (DECISIONS.md D-018); pass strict_locator_passage=False to report only.",
    ),
    "W102_LOCATOR_IS_WHOLE_SECTION": (
        False,
        "source_locator names a whole section rather than a paragraph, "
        "figure or table. Admissible but a weak address.",
    ),
    "W103_LOCATOR_PANEL_SUFFIX": (
        False,
        "source_locator carries a figure panel suffix (e.g. 'Fig 3C'). "
        "Panels are not separately addressable in text, so it resolved to "
        "that figure's caption; the figure number matched exactly.",
    ),
    "W104_FULL_TEXT_FLAG_MISMATCH": (
        False,
        "full_text_available is not true although retrieved text was found "
        "for this study: the record misdescribes its own provenance.",
    ),
    "W105_PAPER_PARSE_WARNINGS": (
        False,
        "The source document did not parse cleanly, so the text checked "
        "against may be incomplete (a truncated harvest file, for example).",
    ),
    "W106_QUOTE_VERY_SHORT": (
        False,
        "The quote is short enough that occurring in the paper is weak "
        "evidence of anything.",
    ),
    "W107_LOCATOR_NUMBERED_BY_ORDER": (
        False,
        "The figure or table this locator resolved to carried no label, so "
        "it is numbered by document order and the number may not be the "
        "publisher's.",
    ),
}

#: Quotes at or below this many normalised characters earn
#: ``W106_QUOTE_VERY_SHORT``. Chosen so that a short clause ("food intake was
#: unchanged") clears it and a bare noun phrase does not.
SHORT_QUOTE_CHARS = 25


@dataclass(frozen=True)
class VerificationFinding:
    """One structured finding about one proposed record."""

    code: str
    blocking: bool
    message: str
    detail: Dict[str, object] = field(default_factory=dict)

    def __str__(self) -> str:
        return "[%s] %s: %s" % (
            "BLOCKING" if self.blocking else "note",
            self.code,
            self.message,
        )

    def to_dict(self) -> Dict[str, object]:
        return {
            "code": self.code,
            "blocking": self.blocking,
            "message": self.message,
            "detail": dict(self.detail),
        }


def _finding(code: str, message: Optional[str] = None, **detail) -> VerificationFinding:
    blocking, default = CODES[code]
    return VerificationFinding(
        code=code,
        blocking=blocking,
        message=message if message is not None else default,
        detail=detail,
    )


@dataclass(frozen=True)
class SupportVerification:
    """The structured result of checking one record against one paper."""

    record_id: Optional[str]
    study_id: Optional[str]
    source_locator: Optional[str]
    findings: Tuple[VerificationFinding, ...] = ()
    #: Locators of the passages in which the quote occurs verbatim.
    quote_found_in: Tuple[str, ...] = ()
    #: Locators the ``source_locator`` resolved to.
    located_passages: Tuple[str, ...] = ()
    #: How the locator resolved; see :meth:`PassageIndex.resolve`.
    locator_resolution: str = "not_checked"
    #: Which identifier matched the paper, or ``None`` when none did.
    paper_key: Optional[str] = None

    @property
    def blocking(self) -> Tuple[VerificationFinding, ...]:
        return tuple(f for f in self.findings if f.blocking)

    @property
    def notes(self) -> Tuple[VerificationFinding, ...]:
        return tuple(f for f in self.findings if not f.blocking)

    @property
    def ok(self) -> bool:
        """True when nothing blocking was found."""
        return not self.blocking

    @property
    def quote_is_verbatim(self) -> bool:
        return bool(self.quote_found_in)

    @property
    def locator_exists(self) -> bool:
        return bool(self.located_passages)

    def codes(self) -> List[str]:
        return [f.code for f in self.findings]

    def to_dict(self) -> Dict[str, object]:
        return {
            "record_id": self.record_id,
            "study_id": self.study_id,
            "source_locator": self.source_locator,
            "ok": self.ok,
            "quote_is_verbatim": self.quote_is_verbatim,
            "locator_exists": self.locator_exists,
            "locator_resolution": self.locator_resolution,
            "paper_key": self.paper_key,
            "quote_found_in": list(self.quote_found_in),
            "located_passages": list(self.located_passages),
            "findings": [f.to_dict() for f in self.findings],
        }


# --------------------------------------------------------------------------- #
# paper lookup
# --------------------------------------------------------------------------- #

PaperSource = Union[ParsedPaper, PassageIndex, None]


def _blank(value) -> bool:
    return value is None or not str(value).strip()


def _probe_keys(record: TierARecord) -> List[str]:
    """Identifier strings to look a record's paper up by, in priority order."""
    out: List[str] = []
    for prefix, value in (
        (None, record.study_id),
        ("doi", record.doi),
        ("pmid", record.pmid),
        ("pmcid", record.pmcid),
    ):
        if _blank(value):
            continue
        bare = str(value).strip().lower()
        out.append(bare)
        if ":" in bare:
            out.append(bare.split(":", 1)[1])
        elif prefix:
            out.append("%s:%s" % (prefix, bare))
    seen = set()
    unique = []
    for key in out:
        if key and key not in seen:
            seen.add(key)
            unique.append(key)
    return unique


def index_papers(
    papers: Union[Mapping[str, ParsedPaper], Iterable[ParsedPaper]],
) -> Dict[str, PassageIndex]:
    """Build a ``{identifier -> PassageIndex}`` map for :func:`gate`.

    Accepts either a mapping (whose keys are added as identifiers as given,
    lower-cased) or a bare iterable of papers (identified by their own
    ``doi``/``pmid``/``pmcid``, bare and prefixed). Each paper is split into
    passages exactly once.
    """
    out: Dict[str, PassageIndex] = {}

    def register(keys: Sequence[str], index: PassageIndex) -> None:
        for key in keys:
            if key:
                out.setdefault(str(key).strip().lower(), index)

    if isinstance(papers, Mapping):
        for key, paper in papers.items():
            index = paper if isinstance(paper, PassageIndex) else build_index(paper)
            keys = [key]
            if isinstance(paper, ParsedPaper):
                keys.extend(paper.study_keys())
            register(keys, index)
        return out

    for paper in papers:
        if isinstance(paper, PassageIndex):  # pragma: no cover - defensive
            continue
        register(paper.study_keys(), build_index(paper))
    return out


def find_paper(
    record: TierARecord, index: Mapping[str, PassageIndex]
) -> Tuple[Optional[PassageIndex], Optional[str]]:
    """The indexed paper this record's identifiers name, and which one matched.

    Exact identifier matching only. A record whose ``study_id`` matches no
    retrieved paper is rejected under SCHEMA.md rule 2.5 rather than checked
    against the nearest available paper, which would be the most damaging
    possible form of helpfulness.
    """
    for key in _probe_keys(record):
        hit = index.get(key)
        if hit is not None:
            return hit, key
    return None, None


# --------------------------------------------------------------------------- #
# the check
# --------------------------------------------------------------------------- #


def _as_index(paper: PaperSource) -> Optional[PassageIndex]:
    if paper is None:
        return None
    if isinstance(paper, PassageIndex):
        return paper
    return build_index(paper)


def verify_support(
    record: TierARecord,
    paper: PaperSource,
    strict_locator_passage: bool = True,
) -> SupportVerification:
    """Check one proposed record's quote and locator against one paper.

    ``paper`` may be a :class:`~syn4rch.curate.jats.ParsedPaper`, a
    pre-built :class:`~syn4rch.curate.passages.PassageIndex`, or ``None`` for
    "no text was retrieved for this study" -- which is itself a verdict
    (``C001_TEXT_NOT_RETRIEVED``), not an error.

    The record is never modified. Returns a :class:`SupportVerification`;
    ``.ok`` is the verdict and ``.findings`` says why.
    """
    findings: List[VerificationFinding] = []
    quote = record.support_quote
    locator = record.source_locator
    index = _as_index(paper)

    if isinstance(paper, ParsedPaper) and paper.parse_warnings:
        findings.append(
            _finding(
                "W105_PAPER_PARSE_WARNINGS",
                parse_warnings=list(paper.parse_warnings)[:5],
                parse_warning_count=len(paper.parse_warnings),
            )
        )

    if index is None or not index.passages:
        findings.append(
            _finding(
                "C001_TEXT_NOT_RETRIEVED",
                detail_reason=(
                    "no paper supplied"
                    if index is None
                    else "the paper parsed to zero passages"
                ),
                study_id=record.study_id,
            )
        )
        return SupportVerification(
            record_id=record.record_id,
            study_id=record.study_id,
            source_locator=locator,
            findings=tuple(findings),
        )

    if _blank(quote):
        findings.append(_finding("C002_MISSING_SUPPORT_QUOTE"))
    if _blank(locator):
        findings.append(_finding("C003_MISSING_SOURCE_LOCATOR"))

    # --- 1. is the quote verbatim in this paper? --------------------------
    needle = normalize_for_match(quote)
    containing: List[Passage] = []
    if needle:
        containing = [
            passage
            for passage in index.passages
            if needle in normalize_for_match(passage.text)
        ]
        if not containing:
            findings.append(
                _finding(
                    "C004_QUOTE_NOT_VERBATIM",
                    normalised_quote=needle[:300],
                    normalised_quote_chars=len(needle),
                    passages_searched=len(index.passages),
                    normalisation=NORMALISATION,
                )
            )
        elif len(needle) <= SHORT_QUOTE_CHARS:
            findings.append(
                _finding("W106_QUOTE_VERY_SHORT", normalised_quote_chars=len(needle))
            )

    # --- 2. does the locator name a passage that exists? ------------------
    located: List[Passage] = []
    resolution = "not_checked"
    if not _blank(locator):
        located, resolution = index.resolve(locator)
        if resolution == "unparseable":
            findings.append(
                _finding(
                    "C005_LOCATOR_UNPARSEABLE",
                    source_locator=locator,
                    example_locators=index.locators()[:8],
                )
            )
        elif not located:
            findings.append(
                _finding(
                    "C006_LOCATOR_NOT_IN_PAPER",
                    source_locator=locator,
                    parsed_as=(
                        parse_locator(locator).describe()
                        if parse_locator(locator)
                        else None
                    ),
                    available_locators=index.locators()[:40],
                    available_locator_count=len(index.passages),
                )
            )
        else:
            if resolution == "panel_suffix":
                findings.append(
                    _finding("W103_LOCATOR_PANEL_SUFFIX", source_locator=locator)
                )
            if resolution == "section_whole":
                findings.append(
                    _finding(
                        "W102_LOCATOR_IS_WHOLE_SECTION",
                        source_locator=locator,
                        paragraphs_in_section=len(located),
                    )
                )
            if any(p.numbered_by_order for p in located):
                findings.append(
                    _finding("W107_LOCATOR_NUMBERED_BY_ORDER", source_locator=locator)
                )

    # --- 3. do the two agree? --------------------------------------------
    if containing and located:
        located_set = {p.locator for p in located}
        if not any(p.locator in located_set for p in containing):
            note = _finding(
                "W101_QUOTE_NOT_IN_LOCATED_PASSAGE",
                source_locator=locator,
                quote_found_in=[p.locator for p in containing],
            )
            if strict_locator_passage:
                note = VerificationFinding(
                    code=note.code,
                    blocking=True,
                    message=note.message,
                    detail=note.detail,
                )
            findings.append(note)

    if record.full_text_available is not True:
        findings.append(
            _finding(
                "W104_FULL_TEXT_FLAG_MISMATCH",
                full_text_available=record.full_text_available,
            )
        )

    return SupportVerification(
        record_id=record.record_id,
        study_id=record.study_id,
        source_locator=locator,
        findings=tuple(findings),
        quote_found_in=tuple(p.locator for p in containing),
        located_passages=tuple(p.locator for p in located),
        locator_resolution=resolution,
        paper_key=None,
    )


# --------------------------------------------------------------------------- #
# the gate
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Rejection:
    """One rejected proposal, with a machine-readable reason.

    ``codes`` is the full set of blocking codes; ``code`` is the first, which
    is the one a report groups by. ``record`` is the proposal exactly as it
    arrived -- unmodified, because nothing in this module repairs a record.
    """

    record: TierARecord
    codes: Tuple[str, ...]
    reasons: Tuple[str, ...]
    verification: Optional[SupportVerification] = None
    schema_findings: Tuple[ValidationError, ...] = ()
    detail: Dict[str, object] = field(default_factory=dict)

    @property
    def code(self) -> str:
        return self.codes[0] if self.codes else "C008_UNREADABLE_PROPOSAL"

    @property
    def record_id(self) -> Optional[str]:
        return getattr(self.record, "record_id", None)

    def to_dict(self) -> Dict[str, object]:
        return {
            "record_id": self.record_id,
            "study_id": getattr(self.record, "study_id", None),
            "source_locator": getattr(self.record, "source_locator", None),
            "code": self.code,
            "codes": list(self.codes),
            "reasons": list(self.reasons),
            "schema_findings": [f.to_dict() for f in self.schema_findings],
            "verification": self.verification.to_dict() if self.verification else None,
            "detail": dict(self.detail),
        }


def gate(
    records: Iterable[TierARecord],
    papers: Union[Mapping[str, ParsedPaper], Iterable[ParsedPaper], None] = None,
    strict_locator_passage: bool = True,
    run_schema_validate: bool = True,
) -> Tuple[List[TierARecord], List[Rejection]]:
    """Admit or reject each proposal. Returns ``(admitted, rejected)``.

    A record is admitted only when ``syn4rch.schema.validate`` reports no hard
    finding **and** :func:`verify_support` reports nothing blocking. Admitted
    records are the same objects that were passed in: the gate selects, it
    does not transform, and in particular it never edits a quote, a locator or
    any other field to make a proposal pass.

    ``papers`` may be a mapping of identifier to paper, an iterable of papers,
    or ``None``. ``None`` means no text is held at all, so every proposal is
    rejected under SCHEMA.md rule 2.5 -- which is the correct answer, not a
    configuration mistake.
    """
    index = index_papers(papers) if papers is not None else {}
    admitted: List[TierARecord] = []
    rejected: List[Rejection] = []

    for record in records:
        codes: List[str] = []
        reasons: List[str] = []
        schema_hard: Tuple[ValidationError, ...] = ()

        if run_schema_validate:
            schema_hard = tuple(hard_findings(validate(record)))
            if schema_hard:
                codes.append("C007_SCHEMA_HARD")
                reasons.append(
                    "%s: %s"
                    % (
                        CODES["C007_SCHEMA_HARD"][1],
                        "; ".join("%s (%s)" % (f.code, f.field) for f in schema_hard),
                    )
                )

        paper_index, paper_key = find_paper(record, index)
        verification = verify_support(
            record, paper_index, strict_locator_passage=strict_locator_passage
        )
        if paper_key is not None:
            verification = SupportVerification(
                record_id=verification.record_id,
                study_id=verification.study_id,
                source_locator=verification.source_locator,
                findings=verification.findings,
                quote_found_in=verification.quote_found_in,
                located_passages=verification.located_passages,
                locator_resolution=verification.locator_resolution,
                paper_key=paper_key,
            )
        for finding in verification.blocking:
            codes.append(finding.code)
            reasons.append(finding.message)

        if codes:
            rejected.append(
                Rejection(
                    record=record,
                    codes=tuple(codes),
                    reasons=tuple(reasons),
                    verification=verification,
                    schema_findings=schema_hard,
                    detail={"paper_key": paper_key},
                )
            )
        else:
            admitted.append(record)

    return admitted, rejected


# Guard against a future edit quietly importing a similarity library into the
# one module whose value depends on not having one. Standard library only
# (DECISIONS D-001) and, specifically, nothing that can approximate.
FORBIDDEN_MODULES = (
    "difflib",
    "rapidfuzz",
    "fuzzywuzzy",
    "Levenshtein",
    "jellyfish",
    "regex",
)


def forbidden_imports() -> List[str]:
    """Approximate-matching modules imported by this module's own source.

    Called by the test suite, not at import time. A non-empty result means
    somebody has put fuzzy matching inside the gate, which would end the one
    guarantee the gate makes.
    """
    import ast
    import inspect

    try:
        source = inspect.getsource(sys.modules[__name__])
    except (OSError, TypeError):  # pragma: no cover - source is always there
        return []
    found: List[str] = []
    for node in ast.walk(ast.parse(source)):
        names: List[str] = []
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            names = [node.module or ""]
        for name in names:
            root = name.split(".")[0]
            if root in FORBIDDEN_MODULES and root not in found:
                found.append(root)
    return found
