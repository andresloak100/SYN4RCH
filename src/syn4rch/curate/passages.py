"""Passages: the addressable units of a retrieved paper. Stage 4 of the protocol.

A :class:`Passage` is a span of text from one paper together with a
``source_locator`` in the style SCHEMA.md 3.6 requires -- ``"Results para 4"``,
``"Fig 3C caption"``, ``"Methods para 2"``, ``"Table 1 caption"``. The locator
is the only thing that lets a reader of ``data/curated/records.jsonl`` go back
to the sentence a record came from, so it has two obligations:

**Stable.** The same XML always yields the same locators. Paragraph numbering
runs over a section's whole subtree in document order (see
:class:`syn4rch.curate.jats.Section`), so it does not move when a publisher
nests a heading one level deeper. Figure numbering uses the publisher's own
label (``"Figure 3"`` -> ``Fig 3``), falling back to document order only when
there is no label, in which case the fallback is recorded.

**Unique within a paper.** Two passages never share a locator; a collision
(two sections with the same heading, two figures labelled "Figure 3") gets a
``" (2)"`` suffix, and that suffix is stripped again when a locator is parsed,
so a record citing either form still resolves.

This module owns the locator *grammar*: it both writes locators
(:func:`split_passages`) and reads them (:func:`parse_locator`). Keeping the
two in one file is deliberate -- a locator this project writes must be a
locator this project can resolve, and that invariant is testable only when
both halves are defined together.

Python 3.10 compatible, standard library only.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .jats import ParsedPaper, classify_section

__all__ = [
    "Passage",
    "LocatorRef",
    "PassageIndex",
    "split_passages",
    "parse_locator",
    "build_index",
    "KIND_TITLE",
    "KIND_PARAGRAPH",
    "KIND_FIGURE_CAPTION",
    "KIND_TABLE_CAPTION",
]

KIND_TITLE = "title"
KIND_PARAGRAPH = "paragraph"
KIND_FIGURE_CAPTION = "figure_caption"
KIND_TABLE_CAPTION = "table_caption"

_WS = re.compile(r"\s+")
#: The dedup suffix :func:`split_passages` appends on a locator collision.
_DEDUP_SUFFIX = re.compile(r"\s*\((\d+)\)\s*$")
#: Trailing words that say "this is a caption" and carry no address.
_CAPTION_WORD = re.compile(r"\s*\b(caption|legend)\b\s*$", re.IGNORECASE)

_FIGURE_LOCATOR = re.compile(
    r"^(?:fig|figure|figs|figures)\.?\s*"
    r"(?P<num>\d+|[ivxlcdm]+)"
    r"(?P<panel>\s*[a-z](?:\s*[-,–—/]\s*[a-z])*)?$",
    re.IGNORECASE,
)
_TABLE_LOCATOR = re.compile(
    r"^(?:table|tbl|tab)\.?\s*(?P<num>\d+|[ivxlcdm]+)"
    r"(?P<panel>\s*[a-z](?:\s*[-,–—/]\s*[a-z])*)?$",
    re.IGNORECASE,
)
_SECTION_PARA_LOCATOR = re.compile(
    r"^(?P<sec>.+?)[,:]?\s*(?:para|paragraph|par)\.?\s*(?P<num>\d+)$",
    re.IGNORECASE,
)


def _collapse(text: Optional[str]) -> str:
    if not text:
        return ""
    return _WS.sub(" ", str(text)).strip()


def _norm_key(text: Optional[str]) -> str:
    """Normalised comparison key for a locator component.

    Lower-cases, collapses whitespace, and strips edge punctuation. This is a
    *locator* normalisation only; it is never applied to quoted text, where
    case is load-bearing (see :mod:`syn4rch.curate.verify`).
    """
    return _collapse(text).lower().strip(" .,:;()[]-")


# --------------------------------------------------------------------------- #
# data model
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Passage:
    """One addressable span of a paper's retrieved text.

    ``locator`` is what a record's ``source_locator`` must say.
    ``section`` is the display name the locator uses (``"Results"``, or a
    non-standard section's own heading). ``text`` is the passage verbatim,
    with whitespace runs collapsed -- the same collapse the verbatim gate
    applies to both sides of its comparison, so a quote copied from here
    verifies against here.
    """

    locator: str
    section: str
    text: str
    kind: str = KIND_PARAGRAPH
    #: Canonical section label from :data:`syn4rch.curate.jats.SECTION_LABELS`.
    section_label: str = "Other"
    #: Paragraph number within the section, or the float's number/ordinal.
    ordinal: Optional[int] = None
    #: True when a figure/table locator had to fall back to document order
    #: because the float carried no usable label.
    numbered_by_order: bool = False

    def to_dict(self) -> Dict[str, object]:
        return {
            "locator": self.locator,
            "section": self.section,
            "section_label": self.section_label,
            "kind": self.kind,
            "ordinal": self.ordinal,
            "numbered_by_order": self.numbered_by_order,
            "text": self.text,
        }


@dataclass(frozen=True)
class LocatorRef:
    """A parsed ``source_locator``: what it addresses, in resolvable terms.

    ``kind`` is one of ``"title"``, ``"section_para"``, ``"section"``,
    ``"figure"``, ``"table"``. ``panel`` holds a figure panel suffix
    (``"Fig 3C"`` -> number ``"3"``, panel ``"c"``) which is recorded but
    never required to match anything: a figure's panels are not separately
    addressable in text, so ``"Fig 3C caption"`` resolves to figure 3's
    caption. The figure *number* must match exactly.
    """

    kind: str
    name: str = ""
    number: Optional[int] = None
    float_number: str = ""
    panel: Optional[str] = None
    raw: str = ""

    def describe(self) -> str:
        if self.kind == "title":
            return "the article title"
        if self.kind == "section_para":
            return "paragraph %s of section %r" % (self.number, self.name)
        if self.kind == "section":
            return "section %r as a whole" % self.name
        if self.kind == "figure":
            return "the caption of figure %s" % self.float_number
        if self.kind == "table":
            return "the caption of table %s" % self.float_number
        return self.raw  # pragma: no cover - defensive


# --------------------------------------------------------------------------- #
# writing locators
# --------------------------------------------------------------------------- #


def _section_display_names(paper: ParsedPaper) -> List[str]:
    """Display name per section, disambiguated, in paper order.

    A canonical section uses its label. An ``Other`` section uses its own
    heading, because that is the only honest address for it; with no heading it
    becomes ``"Section <k>"`` by position.
    """
    names: List[str] = []
    for index, section in enumerate(paper.sections, start=1):
        if section.label != "Other":
            name = section.label
        elif section.title:
            name = _collapse(section.title)
        else:
            name = "Section %d" % index
        names.append(name)
    return _disambiguate(names)


def _disambiguate(names: Sequence[str]) -> List[str]:
    seen: Dict[str, int] = {}
    out: List[str] = []
    for name in names:
        key = _norm_key(name)
        seen[key] = seen.get(key, 0) + 1
        out.append(name if seen[key] == 1 else "%s (%d)" % (name, seen[key]))
    return out


def split_passages(paper: ParsedPaper, include_title: bool = True) -> List[Passage]:
    """Every addressable passage of ``paper``, in document order.

    Deterministic: the same :class:`~syn4rch.curate.jats.ParsedPaper` always
    yields the same list, and every ``locator`` is unique within the result.
    An empty paper yields an empty list rather than a placeholder, because a
    placeholder is a passage a quote could be claimed against.
    """
    passages: List[Passage] = []
    used: Dict[str, int] = {}

    def add(locator: str, section: str, text: str, **kwargs) -> None:
        text = _collapse(text)
        if not text:
            return
        key = _norm_key(locator)
        used[key] = used.get(key, 0) + 1
        final = locator if used[key] == 1 else "%s (%d)" % (locator, used[key])
        passages.append(Passage(final, section, text, **kwargs))

    if include_title and paper.title:
        add(
            "Title",
            "Title",
            paper.title,
            kind=KIND_TITLE,
            section_label="Other",
            ordinal=None,
        )

    for section, display in zip(paper.sections, _section_display_names(paper)):
        for number, paragraph in enumerate(section.paragraphs, start=1):
            add(
                "%s para %d" % (display, number),
                display,
                paragraph,
                kind=KIND_PARAGRAPH,
                section_label=section.label,
                ordinal=number,
            )

    for caption in paper.figures:
        number = caption.number or str(caption.ordinal)
        add(
            "Fig %s caption" % number,
            "Figure",
            caption.text,
            kind=KIND_FIGURE_CAPTION,
            section_label="Other",
            ordinal=caption.ordinal,
            numbered_by_order=caption.number is None,
        )
    for caption in paper.tables:
        number = caption.number or str(caption.ordinal)
        add(
            "Table %s caption" % number,
            "Table",
            caption.text,
            kind=KIND_TABLE_CAPTION,
            section_label="Other",
            ordinal=caption.ordinal,
            numbered_by_order=caption.number is None,
        )

    return passages


# --------------------------------------------------------------------------- #
# reading locators
# --------------------------------------------------------------------------- #


def parse_locator(locator: Optional[str]) -> Optional[LocatorRef]:
    """Parse a ``source_locator`` string, or return ``None`` if it is not one.

    Accepted forms (case-insensitively, with ``"caption"``/``"legend"`` and a
    ``" (2)"`` dedup suffix optional throughout)::

        Title
        Results para 4            Methods, paragraph 2      Discussion par 1
        Fig 3                     Figure 3 caption          Fig 3C-E caption
        Table 1                   Table II caption
        Results                   Abstract

    ``None`` is returned for anything else -- ``"Supplementary Figure 2"``,
    ``"page 7"``, ``"personal communication"`` -- and the gate rejects a
    record whose locator does not parse. That is the intended behaviour: a
    locator naming something outside the retrieved text is not an address this
    dataset can honour (SCHEMA.md rule 2.5).
    """
    raw = _collapse(locator)
    if not raw:
        return None
    text = _DEDUP_SUFFIX.sub("", raw)
    text = _CAPTION_WORD.sub("", text).strip(" .,:;")
    if not text:
        return None

    if _norm_key(text) in ("title", "article title", "paper title"):
        return LocatorRef(kind="title", name="Title", raw=raw)

    match = _FIGURE_LOCATOR.match(text)
    if match:
        return LocatorRef(
            kind="figure",
            float_number=_float_token(match.group("num")),
            panel=_norm_key(match.group("panel")) or None,
            raw=raw,
        )
    match = _TABLE_LOCATOR.match(text)
    if match:
        return LocatorRef(
            kind="table",
            float_number=_float_token(match.group("num")),
            panel=_norm_key(match.group("panel")) or None,
            raw=raw,
        )
    match = _SECTION_PARA_LOCATOR.match(text)
    if match:
        name = _collapse(match.group("sec"))
        if name:
            return LocatorRef(
                kind="section_para",
                name=name,
                number=int(match.group("num")),
                raw=raw,
            )
        return None

    # A bare section name. Allowed, and weaker than a paragraph address: the
    # quote check still has to pass, but the reader has not said where in the
    # section to look. Reported, not rejected.
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9 &/'’.,\-]*", text):
        return LocatorRef(kind="section", name=text, raw=raw)
    return None


def _float_token(token: str) -> str:
    token = _collapse(token)
    return token if token.isdigit() else token.upper()


# --------------------------------------------------------------------------- #
# resolution
# --------------------------------------------------------------------------- #


@dataclass
class PassageIndex:
    """Passages of one paper, indexed for locator resolution.

    Built once per paper so that gating a few thousand records does not
    re-split the same document a few thousand times.
    """

    passages: List[Passage] = field(default_factory=list)
    _by_locator: Dict[str, List[Passage]] = field(default_factory=dict, repr=False)
    _by_section: Dict[str, List[Passage]] = field(default_factory=dict, repr=False)
    _by_float: Dict[Tuple[str, str], List[Passage]] = field(
        default_factory=dict, repr=False
    )

    def __post_init__(self) -> None:
        for passage in self.passages:
            self._by_locator.setdefault(_norm_key(passage.locator), []).append(passage)
            if passage.kind == KIND_PARAGRAPH:
                self._by_section.setdefault(_norm_key(passage.section), []).append(
                    passage
                )
            elif passage.kind in (KIND_FIGURE_CAPTION, KIND_TABLE_CAPTION):
                kind = "figure" if passage.kind == KIND_FIGURE_CAPTION else "table"
                number = _float_number_of(passage)
                self._by_float.setdefault((kind, number.lower()), []).append(passage)

    # --- lookup -----------------------------------------------------------

    def locators(self) -> List[str]:
        return [p.locator for p in self.passages]

    def section_names(self) -> List[str]:
        names: List[str] = []
        for passage in self.passages:
            if passage.kind == KIND_PARAGRAPH and passage.section not in names:
                names.append(passage.section)
        return names

    def resolve(self, locator: Optional[str]) -> Tuple[List[Passage], str]:
        """Resolve a locator to the passages it addresses.

        Returns ``(passages, how)``. ``how`` is a short machine-readable word
        naming the route taken, so a report can show that a locator matched
        only via a documented relaxation:

        ``"exact"``
            The locator string is one this paper's passages carry.
        ``"section_name"``
            A section-paragraph locator whose section name is a synonym of a
            section this paper has (``"Materials and Methods para 2"`` ->
            ``"Methods para 2"``), resolved with the same cue table
            :func:`syn4rch.curate.jats.classify_section` uses. Not fuzzy
            matching: the paragraph number must be exact.
        ``"panel_suffix"``
            A figure/table locator carrying a panel letter (``"Fig 3C"``).
            Panels are not separately addressable in text, so this resolves to
            that float's caption. The number still has to match exactly.
        ``"section_whole"``
            A bare section name; every paragraph of that section.
        ``"unparseable"`` / ``"not_found"``
            Nothing was addressed.
        """
        ref = parse_locator(locator)
        if ref is None:
            return [], "unparseable"
        exact = self._by_locator.get(_norm_key(locator))
        if exact:
            return list(exact), "exact"
        stripped = _norm_key(_DEDUP_SUFFIX.sub("", _collapse(locator)))
        exact = self._by_locator.get(stripped)
        if exact:
            return list(exact), "exact"

        if ref.kind == "title":
            hits = [p for p in self.passages if p.kind == KIND_TITLE]
            return hits, "exact" if hits else "not_found"

        if ref.kind in ("figure", "table"):
            hits = self._by_float.get((ref.kind, ref.float_number.lower()), [])
            if hits:
                return list(hits), "panel_suffix" if ref.panel else "exact"
            return [], "not_found"

        if ref.kind == "section_para":
            for passage in self._matching_section(ref.name):
                if passage.ordinal == ref.number:
                    return [passage], "section_name"
            return [], "not_found"

        if ref.kind == "section":
            hits = self._matching_section(ref.name)
            return hits, "section_whole" if hits else "not_found"

        return [], "not_found"  # pragma: no cover - defensive

    # --- helpers ----------------------------------------------------------

    def _matching_section(self, name: str) -> List[Passage]:
        key = _norm_key(name)
        direct = self._by_section.get(key)
        if direct:
            return list(direct)
        # Synonym route: the name classifies to the same canonical label as a
        # section this paper has. Uses the published cue table, so the set of
        # accepted synonyms is the documented one, not an open-ended guess.
        label = classify_section(None, name)
        if label == "Other":
            return []
        hits = [
            p
            for p in self.passages
            if p.kind == KIND_PARAGRAPH and p.section_label == label
        ]
        return hits


def _float_number_of(passage: Passage) -> str:
    """The number a figure/table passage is addressed by."""
    match = re.match(
        r"^(?:fig|table)\s*(\S+)", _collapse(passage.locator), re.IGNORECASE
    )
    if match:
        return _DEDUP_SUFFIX.sub("", match.group(1)).strip()
    return str(passage.ordinal or "")


def build_index(paper: ParsedPaper, include_title: bool = True) -> PassageIndex:
    """Split ``paper`` and index the result for locator resolution."""
    return PassageIndex(passages=split_passages(paper, include_title=include_title))
