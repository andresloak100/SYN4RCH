"""JATS XML -> :class:`ParsedPaper`. Stage 3/4 of docs/CURATION_PROTOCOL.md.

Reads the ``fullTextXML`` documents that ``syn4rch.harvest.europepmc`` writes to
``data/raw/fulltext/<PMCID>.xml`` and turns one into an addressable object:
bibliographic identity, ordered sections, figure captions and table captions.

Three properties this module is built around
--------------------------------------------

**It never raises on bad input.** ``data/raw/fulltext/`` is written by a
concurrent harvest process, so a file on disk may be a few hundred bytes of a
document whose closing tags have not been written yet. A truncated or malformed
file is a *normal* input here, not an error: :func:`parse_jats` recovers every
element the parser managed to see, records what went wrong in
``ParsedPaper.parse_warnings``, and returns. A caller that wants to know whether
the parse was clean reads that list; a caller that does not gets partial text
rather than a traceback.

**It never invents text.** Inline markup (``italic``, ``sup``, ``sub``,
``xref``, ``bold``, ...) is flattened by concatenating character data in
document order, which is exactly what the published sentence says. Nothing is
reordered, substituted or paraphrased. The only transformation applied to
character data is that runs of whitespace collapse to a single space, which is
the same normalisation the verbatim gate uses (``verify.normalize_for_match``),
so a quote copied out of a passage still verifies against it.

**Inline citation markers stay in the text.** ``<xref ref-type="bibr">12</xref>``
after a word flattens to ``...food intake12``, because that is what the XML
says. Dropping it would mean the gate verifies quotes against a string that is
not the source. The consequence, which is deliberate: a quote must be copied
from the passage text this module produces (as
``docs/CURATION_PROTOCOL.md`` 3.1 already requires), not retyped from a PDF.

Python 3.10 compatible, standard library only.
"""

from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

__all__ = [
    "SECTION_LABELS",
    "Caption",
    "ParsedPaper",
    "Section",
    "parse_jats",
    "parse_jats_file",
    "flatten_text",
]

# --------------------------------------------------------------------------- #
# canonical section labels
# --------------------------------------------------------------------------- #

#: The canonical section labels this module emits, in the order a paper has
#: them. ``"Other"`` is a real answer: a section whose heading does not say
#: what it is keeps its own title and is labelled ``Other`` rather than being
#: guessed into one of the five.
SECTION_LABELS = (
    "Abstract",
    "Introduction",
    "Methods",
    "Results",
    "Discussion",
    "Other",
)

#: Lower-cased cues matched against ``sec-type`` and the section's own title,
#: in this order. First hit wins, which is why ``results`` precedes
#: ``discussion``: a "Results and Discussion" section is labelled ``Results``,
#: because that is where its measurements are reported.
_SECTION_CUES: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    (
        "Introduction",
        ("intro", "introduction", "background"),
    ),
    (
        "Methods",
        (
            "methods",
            "method",
            "materials and methods",
            "methods and materials",
            "material and methods",
            "materials & methods",
            "experimental procedures",
            "experimental procedure",
            "experimental section",
            "star methods",
            "methodology",
            "subjects and methods",
            "patients and methods",
            "online methods",
        ),
    ),
    (
        "Results",
        ("results", "result", "results and discussion", "findings"),
    ),
    (
        "Discussion",
        (
            "discussion",
            "discussion and conclusions",
            "conclusion",
            "conclusions",
            "concluding remarks",
            "general discussion",
        ),
    ),
)

#: Tags whose character data belongs to a *different* flow than the paragraph
#: that contains them, so a space is forced on either side when flattening.
#: Everything not listed here (``italic``, ``sup``, ``sub``, ``xref``,
#: ``bold``, ``sc``, ``monospace``, ``styled-content``, ``named-content``,
#: inline formulae, ...) is inline and concatenated with no added character.
_BLOCK_TAGS = frozenset(
    {
        "p",
        "sec",
        "title",
        "label",
        "caption",
        "list",
        "list-item",
        "def-list",
        "def-item",
        "term",
        "def",
        "disp-quote",
        "disp-formula",
        "table",
        "table-wrap",
        "tbody",
        "thead",
        "tr",
        "td",
        "th",
        "fig",
        "graphic",
        "boxed-text",
        "break",
        "array",
        "abstract",
        "ack",
        "statement",
        "supplementary-material",
    }
)

#: Tags dropped wholesale when flattening running text: their content is
#: apparatus rather than the author's sentence. ``xref`` and ``ext-link`` are
#: deliberately *not* here -- their character data is part of the sentence
#: (see the module docstring).
_DROP_TAGS = frozenset(
    {
        "table-wrap-foot",
        "fn",
        "fn-group",
        "object-id",
        "processing-meta",
    }
)

#: Tags that are *floats*: they carry their own caption and must not be
#: swallowed into the running text of the section that contains them.
_FLOAT_TAGS = frozenset({"fig", "fig-group", "table-wrap", "table-wrap-group"})

_WS = re.compile(r"\s+")


# --------------------------------------------------------------------------- #
# data model
# --------------------------------------------------------------------------- #


@dataclass
class Section:
    """One top-level section of a paper.

    ``label`` is one of :data:`SECTION_LABELS`. ``title`` is the paper's own
    heading, verbatim, or ``None`` when the section has none. ``paragraphs``
    holds every ``<p>`` in the section's subtree in document order, including
    those inside subsections, and *excluding* figure and table captions, which
    are addressed separately. Flattening subsections is what makes
    ``"Results para 4"`` a stable address: paragraph numbering does not depend
    on how deeply a publisher chose to nest its headings.
    """

    label: str
    title: Optional[str]
    paragraphs: List[str] = field(default_factory=list)
    #: Headings of nested subsections, in order, for a reader's orientation.
    subsection_titles: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        return {
            "label": self.label,
            "title": self.title,
            "paragraphs": list(self.paragraphs),
            "subsection_titles": list(self.subsection_titles),
        }


@dataclass
class Caption:
    """A figure or table caption.

    ``label`` is the publisher's own label flattened verbatim (``"Figure 3"``,
    ``"Fig. 3"``, ``"Table 1"``); ``number`` is the integer or roman numeral
    pulled out of it, or ``None`` when the label carries none. ``ordinal`` is
    the 1-based position of this float among floats of its kind in document
    order, and is what addressing falls back to when ``number`` is ``None``.
    """

    kind: str  # "figure" | "table"
    label: Optional[str]
    number: Optional[str]
    ordinal: int
    title: Optional[str]
    text: str

    def to_dict(self) -> Dict[str, object]:
        return {
            "kind": self.kind,
            "label": self.label,
            "number": self.number,
            "ordinal": self.ordinal,
            "title": self.title,
            "text": self.text,
        }


@dataclass
class ParsedPaper:
    """Everything this module could read out of one JATS document.

    Every field is independently optional. A ``ParsedPaper`` with no title, no
    sections and three ``parse_warnings`` is the correct result for a file the
    harvest was still writing, and callers must handle it: see
    :func:`syn4rch.curate.verify.verify_support`, which treats a paper with no
    passages as "text not retrieved" under SCHEMA.md rule 2.5.
    """

    pmcid: Optional[str] = None
    pmid: Optional[str] = None
    doi: Optional[str] = None
    title: Optional[str] = None
    journal: Optional[str] = None
    year: Optional[int] = None
    first_author: Optional[str] = None
    last_author: Optional[str] = None
    last_author_affiliation: Optional[str] = None
    sections: List[Section] = field(default_factory=list)
    figures: List[Caption] = field(default_factory=list)
    tables: List[Caption] = field(default_factory=list)
    parse_warnings: List[str] = field(default_factory=list)
    #: Where the XML came from, when it came from a file.
    source_path: Optional[str] = None

    # --- convenience ------------------------------------------------------

    @property
    def is_empty(self) -> bool:
        """True when no readable text was recovered at all.

        The title alone does not count as text: a paper known only by its
        title has no passage a quote could be checked against.
        """
        return not any(s.paragraphs for s in self.sections) and not (
            self.figures or self.tables
        )

    def study_keys(self) -> List[str]:
        """Identifier strings this paper can be looked up by.

        Includes the bare ids and the ``doi:``/``pmid:``/``pmcid:`` prefixed
        forms that SCHEMA.md 3.1 uses for ``study_id``. Lower-cased; the
        caller is expected to lower-case its probe too.
        """
        out: List[str] = []
        for prefix, value in (
            ("doi", self.doi),
            ("pmid", self.pmid),
            ("pmcid", self.pmcid),
            ("pmc", self.pmcid),
        ):
            if not value:
                continue
            bare = str(value).strip().lower()
            if not bare:
                continue
            out.append(bare)
            if not bare.startswith(prefix + ":"):
                out.append("%s:%s" % (prefix, bare))
        seen = set()
        unique = []
        for key in out:
            if key not in seen:
                seen.add(key)
                unique.append(key)
        return unique

    def section(self, label: str) -> Optional[Section]:
        """First section with this canonical label, or ``None``."""
        for sec in self.sections:
            if sec.label == label:
                return sec
        return None

    def to_dict(self) -> Dict[str, object]:
        return {
            "pmcid": self.pmcid,
            "pmid": self.pmid,
            "doi": self.doi,
            "title": self.title,
            "journal": self.journal,
            "year": self.year,
            "first_author": self.first_author,
            "last_author": self.last_author,
            "last_author_affiliation": self.last_author_affiliation,
            "sections": [s.to_dict() for s in self.sections],
            "figures": [c.to_dict() for c in self.figures],
            "tables": [c.to_dict() for c in self.tables],
            "parse_warnings": list(self.parse_warnings),
            "source_path": self.source_path,
        }


# --------------------------------------------------------------------------- #
# namespace-agnostic element access
# --------------------------------------------------------------------------- #


def _local(tag) -> str:
    """Local name of an element tag, with any XML namespace stripped.

    Europe PMC serves JATS in no namespace, but publisher-supplied and
    archive-converted files sometimes carry one (and ``xlink:`` always appears
    on attributes). Comparing local names means both parse identically instead
    of one silently yielding an empty paper.
    """
    if not isinstance(tag, str):
        return ""
    if tag.startswith("{"):
        return tag.split("}", 1)[1]
    return tag


def _iter_children(elem, *names: str):
    wanted = set(names)
    for child in list(elem):
        if _local(child.tag) in wanted:
            yield child


def _find_child(elem, *names: str):
    for child in _iter_children(elem, *names):
        return child
    return None


def _iter_descendants(elem, *names: str):
    """Every descendant (and ``elem`` itself) whose local name matches."""
    wanted = set(names)
    stack = [elem]
    while stack:
        node = stack.pop()
        if _local(node.tag) in wanted:
            yield node
        stack.extend(reversed(list(node)))


def _attr(elem, name: str) -> Optional[str]:
    """Attribute lookup that ignores namespaces on the attribute name."""
    if elem is None:
        return None
    value = elem.get(name)
    if value is not None:
        return value
    for key, val in elem.attrib.items():
        if _local(key) == name:
            return val
    return None


# --------------------------------------------------------------------------- #
# text flattening
# --------------------------------------------------------------------------- #


def _collapse(text: str) -> str:
    return _WS.sub(" ", text).strip()


def _flatten(elem, skip_tags=frozenset(), drop_floats: bool = False) -> str:
    parts: List[str] = []

    def walk(node, is_root: bool) -> None:
        tag = _local(node.tag)
        dropped = tag in _DROP_TAGS or tag in skip_tags
        if drop_floats and tag in _FLOAT_TAGS and not is_root:
            dropped = True
        if dropped:
            # The element's own content goes, but its tail is character data
            # belonging to the parent's flow and must be kept.
            if not is_root and node.tail:
                parts.append(node.tail)
            return
        block = tag in _BLOCK_TAGS
        if block and not is_root:
            parts.append(" ")
        if node.text:
            parts.append(node.text)
        for child in list(node):
            walk(child, False)
        if block and not is_root:
            parts.append(" ")
        if not is_root and node.tail:
            parts.append(node.tail)

    if elem is None:
        return ""
    walk(elem, True)
    return _collapse("".join(parts))


def flatten_text(elem, skip_tags=frozenset(), drop_floats: bool = False) -> str:
    """Flatten one element's subtree to plain text.

    Character data is concatenated in document order. Inline markup adds
    nothing; a block-level child (another ``<p>``, a list item, a table cell)
    is surrounded by a space so that two blocks cannot be glued into one word.
    Whitespace runs then collapse to a single space.

    ``drop_floats`` removes nested ``<fig>`` / ``<table-wrap>`` subtrees, which
    is what section bodies want: a float carries its own caption and is
    addressed on its own.
    """
    return _flatten(elem, skip_tags=skip_tags, drop_floats=drop_floats)


def _text_or_none(elem, **kwargs) -> Optional[str]:
    if elem is None:
        return None
    text = _flatten(elem, **kwargs)
    return text or None


# --------------------------------------------------------------------------- #
# tolerant XML parsing
# --------------------------------------------------------------------------- #


def _salvage(text: str, warnings: List[str]):
    """Recover as much of a malformed/truncated document as the parser saw.

    Uses :class:`xml.etree.ElementTree.XMLPullParser`: it hands out each
    element at its *start* tag and keeps filling that same object as parsing
    proceeds, so when the feed dies partway through, every element already
    opened is still reachable from the root with whatever text and children
    were read. That is precisely the "the harvest is still writing this file"
    case. Returns the root element, or ``None`` if not even one start tag was
    readable.
    """
    parser = ET.XMLPullParser(events=("start",))
    root = None
    try:
        parser.feed(text)
    except ET.ParseError as exc:
        warnings.append("recovery: feed stopped at %s" % exc)
    except Exception as exc:  # pragma: no cover - defensive
        warnings.append("recovery: feed raised %s: %s" % (type(exc).__name__, exc))
    try:
        parser.close()
    except Exception:
        # Expected for a truncated document: unclosed elements.
        pass
    events = parser.read_events()
    while True:
        try:
            _event, elem = next(events)
        except StopIteration:
            break
        except Exception as exc:
            warnings.append(
                "recovery: event stream stopped (%s: %s)" % (type(exc).__name__, exc)
            )
            break
        if root is None:
            root = elem
    if root is None:
        warnings.append("recovery: no element could be read")
    return root


def _tag_scan(prefix: str):
    """Open-tag stack of a prefix of an XML document, as the source spells it.

    Returns ``(stack, cut)``: the element names left open, outermost first, and
    the offset to truncate at so that the result does not end inside a tag.

    A deliberately small scanner, not an XML parser. It exists only to work out
    which closing tags a truncated document is missing, and the names it
    returns are the *source* names (``j:p``, not ``{ns}p``) because they are
    going to be written back as text. If it gets anything wrong, the repair it
    feeds fails to parse and the caller falls back to
    :func:`_salvage`, so being wrong is safe. Comments, CDATA sections,
    processing instructions and declarations are skipped.
    """
    stack: List[str] = []
    index = 0
    size = len(prefix)
    while True:
        start = prefix.find("<", index)
        if start < 0:
            return stack, size
        if prefix.startswith("<!--", start):
            end = prefix.find("-->", start + 4)
            if end < 0:
                return stack, start
            index = end + 3
            continue
        if prefix.startswith("<![CDATA[", start):
            end = prefix.find("]]>", start + 9)
            if end < 0:
                return stack, start
            index = end + 3
            continue
        if prefix.startswith("<?", start):
            end = prefix.find("?>", start + 2)
            if end < 0:
                return stack, start
            index = end + 2
            continue
        if prefix.startswith("<!", start):
            # A declaration; skip a balanced internal subset, then to '>'.
            cursor = start + 2
            depth = 0
            while cursor < size:
                char = prefix[cursor]
                if char == "[":
                    depth += 1
                elif char == "]":
                    depth -= 1
                elif char == ">" and depth <= 0:
                    break
                cursor += 1
            if cursor >= size:
                return stack, start
            index = cursor + 1
            continue
        cursor = start + 1
        quote = ""
        while cursor < size:
            char = prefix[cursor]
            if quote:
                if char == quote:
                    quote = ""
            elif char in "\"'":
                quote = char
            elif char == ">":
                break
            cursor += 1
        if cursor >= size:
            # Truncated inside a tag: drop the partial tag.
            return stack, start
        body = prefix[start + 1 : cursor]
        index = cursor + 1
        if body.startswith("/"):
            name = body[1:].strip()
            if name in stack:
                while stack:
                    if stack.pop() == name:
                        break
            continue
        if body.rstrip().endswith("/"):
            continue
        parts = body.split()
        if parts:
            stack.append(parts[0])
    return stack, size  # pragma: no cover - unreachable


def _offset_of(text: str, position) -> int:
    """Character offset of an :class:`ET.ParseError` ``(line, column)``."""
    if not position:
        return len(text)
    line, column = position
    lines = text.splitlines(keepends=True)
    if line < 1:
        return 0
    if line > len(lines):
        return len(text)
    offset = sum(len(l) for l in lines[: line - 1])
    return min(len(text), offset + max(0, column))


def _repair(text: str, exc, warnings: List[str]):
    """Re-parse the document up to the failure, closing what was left open.

    This is the recovery that matters for ``data/raw/fulltext/``. The raw
    salvage in :func:`_salvage` loses the character data immediately before
    the break, because the tree builder has not flushed it when the parse
    dies; supplying the missing closing tags makes that last sentence --
    often a whole Results paragraph -- parse normally.

    It repairs *markup*, never text: not one character of character data is
    added, removed or changed, so a quote verified against the result is still
    a quote from the file on disk. Returns the root element, or ``None`` when
    the repair did not parse (the caller then falls back to the raw salvage).
    """
    prefix = text[: _offset_of(text, getattr(exc, "position", None))]
    stack, cut = _tag_scan(prefix)
    prefix = prefix[:cut]
    if not stack or not prefix.strip():
        return None
    patched = prefix + "".join("</%s>" % name for name in reversed(stack))
    try:
        root = ET.fromstring(patched)
    except ET.ParseError:
        return None
    except Exception:  # pragma: no cover - defensive
        return None
    warnings.append(
        "recovered the readable prefix by closing %d unclosed element(s) "
        "(%s); no character data was altered"
        % (len(stack), ", ".join(reversed(stack)))
    )
    return root


def _parse_root(text: str, warnings: List[str]):
    if text is None:
        warnings.append("input was None")
        return None
    if not text.strip():
        warnings.append("input was empty")
        return None
    try:
        return ET.fromstring(text)
    except ET.ParseError as exc:
        warnings.append("XML is not well-formed (%s); parsed what was readable" % exc)
        repaired = _repair(text, exc, warnings)
        if repaired is not None:
            return repaired
    except Exception as exc:  # pragma: no cover - defensive
        warnings.append(
            "XML parse raised %s: %s; parsed what was readable"
            % (type(exc).__name__, exc)
        )
    return _salvage(text, warnings)


# --------------------------------------------------------------------------- #
# front matter
# --------------------------------------------------------------------------- #

_PMCID = re.compile(r"PMC\d+", re.IGNORECASE)
_YEAR = re.compile(r"\b(1[6-9]\d{2}|20\d{2}|21\d{2})\b")
_LABEL_NUMBER = re.compile(
    r"(?:fig(?:ure)?|table|tbl|scheme|chart)\.?\s*([0-9]+|[IVXLCDM]+)\b",
    re.IGNORECASE,
)
_TRAILING_NUMBER = re.compile(r"([0-9]+|[IVXLCDM]+)\s*$")


def _collect_ids(root, paper: ParsedPaper) -> None:
    for node in _iter_descendants(root, "article-id"):
        kind = (_attr(node, "pub-id-type") or "").strip().lower()
        value = _collapse(node.text or "")
        if not value:
            continue
        if kind in ("pmcid", "pmc") and paper.pmcid is None:
            match = _PMCID.search(value)
            if match:
                paper.pmcid = match.group(0).upper()
            else:
                digits = re.sub(r"\D", "", value)
                paper.pmcid = ("PMC" + digits) if digits else None
        elif kind == "pmid" and paper.pmid is None:
            paper.pmid = value
        elif kind == "doi" and paper.doi is None:
            paper.doi = value.lower()


def _author_name(contrib) -> Optional[str]:
    name = _find_child(contrib, "name", "string-name", "name-alternatives")
    if name is None:
        collab = _find_child(contrib, "collab")
        return _text_or_none(collab)
    if _local(name.tag) == "name-alternatives":
        name = _find_child(name, "name", "string-name") or name
    surname = _text_or_none(_find_child(name, "surname"))
    given = _text_or_none(_find_child(name, "given-names"))
    if surname and given:
        return "%s %s" % (surname, given)
    if surname:
        return surname
    flat = _flatten(name)
    return flat or None


def _affiliation_text(aff) -> Optional[str]:
    """Flatten an ``<aff>``, dropping its presentational label/sup marker."""
    return _text_or_none(aff, skip_tags=frozenset({"label", "sup"}))


def _collect_authors(root, paper: ParsedPaper) -> None:
    front = _find_child(root, "front")
    scope = front if front is not None else root
    contribs: List = []
    for group in _iter_descendants(scope, "contrib-group"):
        for contrib in _iter_children(group, "contrib"):
            kind = (_attr(contrib, "contrib-type") or "author").strip().lower()
            if kind and kind != "author":
                continue
            contribs.append(contrib)
    if not contribs:
        return
    names = [(_author_name(c), c) for c in contribs]
    named = [(n, c) for n, c in names if n]
    if not named:
        return
    paper.first_author = named[0][0]
    paper.last_author = named[-1][0]

    # Affiliation of the last author, resolved through the author's own
    # <xref ref-type="aff"> when present. No nearest-neighbour guessing: a
    # wrong affiliation silently breaks the lab-level leakage control that
    # SCHEMA.md 5 calls this corpus's largest risk.
    affs: Dict[str, str] = {}
    all_affs: List[str] = []
    for aff in _iter_descendants(scope, "aff"):
        text = _affiliation_text(aff)
        if not text:
            continue
        all_affs.append(text)
        aff_id = _attr(aff, "id")
        if aff_id:
            affs.setdefault(aff_id, text)

    last_contrib = named[-1][1]
    rids: List[str] = []
    for xref in _iter_descendants(last_contrib, "xref"):
        if (_attr(xref, "ref-type") or "").strip().lower() != "aff":
            continue
        rid = _attr(xref, "rid") or ""
        rids.extend(r for r in rid.split() if r)
    resolved = [affs[r] for r in rids if r in affs]
    if resolved:
        paper.last_author_affiliation = "; ".join(dict.fromkeys(resolved))
    elif len(all_affs) == 1:
        # Exactly one affiliation in the whole front matter: it is the only
        # thing it can be. More than one and unlinked -> None, not a guess.
        paper.last_author_affiliation = all_affs[0]
    elif all_affs:
        paper.parse_warnings.append(
            "last author has no <xref ref-type='aff'> and the paper lists %d "
            "affiliations; last_author_affiliation left null" % len(all_affs)
        )


def _collect_bibliographic(root, paper: ParsedPaper) -> None:
    front = _find_child(root, "front")
    scope = front if front is not None else root

    for group in _iter_descendants(scope, "title-group"):
        title = _find_child(group, "article-title")
        if title is not None:
            paper.title = _text_or_none(title)
            break
    if paper.title is None:
        for title in _iter_descendants(scope, "article-title"):
            paper.title = _text_or_none(title)
            break

    for node in _iter_descendants(scope, "journal-title"):
        paper.journal = _text_or_none(node)
        break
    if paper.journal is None:
        for node in _iter_descendants(scope, "journal-id"):
            paper.journal = _text_or_none(node)
            break

    # Year: prefer an electronic/print publication date over "accepted" or
    # "received" history dates, which is why <history> is not searched.
    preferred = ("epub", "ppub", "collection", "pub", "epub-ppub", "")
    candidates: List[Tuple[int, int]] = []
    for date in _iter_descendants(scope, "pub-date"):
        kind = (_attr(date, "pub-type") or _attr(date, "date-type") or "").strip().lower()
        year_text = _text_or_none(_find_child(date, "year")) or ""
        match = _YEAR.search(year_text)
        if not match:
            continue
        rank = preferred.index(kind) if kind in preferred else len(preferred)
        candidates.append((rank, int(match.group(1))))
    if candidates:
        candidates.sort(key=lambda pair: (pair[0], pair[1]))
        paper.year = candidates[0][1]


# --------------------------------------------------------------------------- #
# sections
# --------------------------------------------------------------------------- #


def classify_section(sec_type: Optional[str], title: Optional[str]) -> str:
    """Canonical label for one section, from its ``sec-type`` and heading.

    ``sec-type`` is consulted first because it is the publisher's own
    assertion; the heading text is the fallback. Matching is on the normalised
    string and on its first ``|``-separated component (JATS allows compound
    types such as ``materials|methods``). Anything unrecognised is ``Other``:
    guessing a section's identity would put a figure legend in "Results" and
    make every locator built on it wrong.
    """
    probes: List[str] = []
    for raw in (sec_type, title):
        if not raw:
            continue
        norm = _collapse(str(raw)).lower().strip(" .:;-")
        if not norm:
            continue
        probes.append(norm)
        if "|" in norm:
            probes.extend(part.strip() for part in norm.split("|") if part.strip())
    for probe in probes:
        for label, cues in _SECTION_CUES:
            if probe in cues:
                return label
    # Substring fallback, applied only to the heading, and only as whole
    # words, so "Results" matches "3. Results" and "Methods" matches
    # "Materials and Methods" without "Methodological notes" matching.
    for probe in probes:
        words = set(re.findall(r"[a-z]+", probe))
        for label, cues in _SECTION_CUES:
            for cue in cues:
                cue_words = set(re.findall(r"[a-z]+", cue))
                if cue_words and cue_words <= words:
                    return label
    return "Other"


def _section_paragraphs(sec) -> Tuple[List[str], List[str]]:
    """Paragraphs and subsection headings of one section, in document order.

    Walks the subtree once. ``<p>`` elements inside floats
    (``<fig>``, ``<table-wrap>``) are skipped: those are captions and are
    addressed as captions.
    """
    paragraphs: List[str] = []
    subtitles: List[str] = []

    def walk(node, depth: int) -> None:
        for child in list(node):
            tag = _local(child.tag)
            if tag in _FLOAT_TAGS or tag in _DROP_TAGS:
                continue
            if tag == "p":
                text = _flatten(child, drop_floats=True)
                if text:
                    paragraphs.append(text)
            elif tag == "title":
                if depth > 0:
                    text = _flatten(child)
                    if text:
                        subtitles.append(text)
            elif tag == "sec":
                walk(child, depth + 1)
            elif tag in ("list", "def-list", "disp-quote", "boxed-text", "statement"):
                text = _flatten(child, drop_floats=True)
                if text:
                    paragraphs.append(text)
            else:
                walk(child, depth)

    walk(sec, 0)
    return paragraphs, subtitles


def _collect_abstract(root, paper: ParsedPaper) -> None:
    front = _find_child(root, "front")
    scope = front if front is not None else root
    for abstract in _iter_descendants(scope, "abstract"):
        kind = (_attr(abstract, "abstract-type") or "").strip().lower()
        if kind in ("graphical", "teaser", "video", "key-points", "author-highlights"):
            continue
        paragraphs, subtitles = _section_paragraphs(abstract)
        if not paragraphs:
            flat = _flatten(abstract, drop_floats=True, skip_tags=frozenset({"title"}))
            if flat:
                paragraphs = [flat]
        if not paragraphs:
            continue
        paper.sections.append(
            Section(
                label="Abstract",
                title=_text_or_none(_find_child(abstract, "title")),
                paragraphs=paragraphs,
                subsection_titles=subtitles,
            )
        )
        return
    paper.parse_warnings.append("no abstract found")


def _collect_body_sections(root, paper: ParsedPaper) -> None:
    body = _find_child(root, "body")
    if body is None:
        paper.parse_warnings.append("no <body> element found")
        return
    top_secs = list(_iter_children(body, "sec"))
    if not top_secs:
        # Some documents put paragraphs straight under <body> with no <sec>.
        paragraphs, _ = _section_paragraphs(body)
        if paragraphs:
            paper.sections.append(
                Section(label="Other", title=None, paragraphs=paragraphs)
            )
            paper.parse_warnings.append(
                "<body> has no <sec>; its paragraphs were kept as one "
                "'Other' section"
            )
        else:
            paper.parse_warnings.append("<body> contains no paragraphs")
        return
    for sec in top_secs:
        title = _text_or_none(_find_child(sec, "title"))
        label = classify_section(_attr(sec, "sec-type"), title)
        paragraphs, subtitles = _section_paragraphs(sec)
        if not paragraphs and not title:
            continue
        paper.sections.append(
            Section(
                label=label,
                title=title,
                paragraphs=paragraphs,
                subsection_titles=subtitles,
            )
        )


# --------------------------------------------------------------------------- #
# floats
# --------------------------------------------------------------------------- #


def _number_token(token: str) -> str:
    return token if token.isdigit() else token.upper()


def _label_number(label: Optional[str], caption_title: Optional[str]) -> Optional[str]:
    """The figure/table number named by a label, or ``None``.

    Read from the label first (``"Figure 3"`` -> ``"3"``, ``"Table II"`` ->
    ``"II"``), then from the caption's own title for publishers that put the
    label there. A bare numeric label (``"3"``) is accepted. Nothing is
    inferred from document order here: that is what ``Caption.ordinal`` is
    for, and the two are kept apart so a locator can say which one it used.
    """
    for text in (label, caption_title):
        if not text:
            continue
        match = _LABEL_NUMBER.search(text)
        if match:
            return _number_token(match.group(1))
    if label:
        match = _TRAILING_NUMBER.search(_collapse(label))
        if match:
            return _number_token(match.group(1))
    return None


def _collect_floats(root, paper: ParsedPaper) -> None:
    fig_ordinal = 0
    table_ordinal = 0
    for node in _iter_descendants(root, "fig", "table-wrap"):
        tag = _local(node.tag)
        kind = "figure" if tag == "fig" else "table"
        label = _text_or_none(_find_child(node, "label"))
        caption = _find_child(node, "caption")
        caption_title = _text_or_none(_find_child(caption, "title")) if caption is not None else None
        if caption is not None:
            body_text = _flatten(caption, skip_tags=frozenset({"title"}))
        else:
            body_text = ""
        pieces = [piece for piece in (label, caption_title, body_text) if piece]
        text = _collapse(" ".join(pieces))
        if not text:
            continue
        if kind == "figure":
            fig_ordinal += 1
            ordinal = fig_ordinal
        else:
            table_ordinal += 1
            ordinal = table_ordinal
        entry = Caption(
            kind=kind,
            label=label,
            number=_label_number(label, caption_title),
            ordinal=ordinal,
            title=caption_title,
            text=text,
        )
        if kind == "figure":
            paper.figures.append(entry)
        else:
            paper.tables.append(entry)


# --------------------------------------------------------------------------- #
# entry points
# --------------------------------------------------------------------------- #


def parse_jats(xml_text: str, source_path: Optional[str] = None) -> ParsedPaper:
    """Parse JATS ``fullTextXML`` into a :class:`ParsedPaper`.

    **Never raises.** Malformed, truncated, empty and non-XML input all return
    a ``ParsedPaper`` whose ``parse_warnings`` says what happened and whose
    other fields hold whatever was readable.
    """
    warnings: List[str] = []
    paper = ParsedPaper(parse_warnings=warnings, source_path=source_path)
    root = _parse_root(xml_text, warnings)
    if root is None:
        return paper

    root_tag = _local(root.tag)
    if root_tag not in ("article", "pmc-articleset", "articles", "result"):
        warnings.append("unexpected root element <%s>" % root_tag)
    if root_tag != "article":
        inner = None
        for node in _iter_descendants(root, "article"):
            inner = node
            break
        if inner is not None:
            root = inner

    for step, func in (
        ("ids", _collect_ids),
        ("bibliographic metadata", _collect_bibliographic),
        ("authors", _collect_authors),
        ("abstract", _collect_abstract),
        ("body sections", _collect_body_sections),
        ("figure and table captions", _collect_floats),
    ):
        try:
            func(root, paper)
        except Exception as exc:  # pragma: no cover - last-resort isolation
            # One broken subtree must not cost the caller the rest of the
            # document; the loss is recorded rather than hidden.
            warnings.append(
                "while reading %s: %s: %s" % (step, type(exc).__name__, exc)
            )

    if paper.is_empty:
        warnings.append("no readable passages were recovered")
    return paper


def parse_jats_file(path: str) -> ParsedPaper:
    """Read and parse one file. Never raises; I/O failures become warnings.

    A file whose name is a PMCID (how ``harvest.europepmc`` writes them)
    supplies ``pmcid`` when the XML itself does not carry one. That is a fact
    about where the bytes came from, not an inference about their content.
    """
    text = ""
    warnings: List[str] = []
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            text = handle.read()
    except OSError as exc:
        paper = ParsedPaper(source_path=path)
        paper.parse_warnings.append("could not read %s: %s" % (path, exc))
        return paper
    paper = parse_jats(text, source_path=path)
    paper.parse_warnings[:0] = warnings
    if paper.pmcid is None:
        stem = os.path.splitext(os.path.basename(path))[0]
        match = _PMCID.fullmatch(stem.strip())
        if match:
            paper.pmcid = match.group(0).upper()
    return paper
