"""Curation: retrieved open-access full text -> schema-valid Tier A records.

Implements stages 4-7 of ``docs/CURATION_PROTOCOL.md``. The ordering of this
package is the ordering of the protocol, and each module reads only the one
above it:

:mod:`syn4rch.curate.jats`
    Parse a Europe PMC ``fullTextXML`` document into a
    :class:`~syn4rch.curate.jats.ParsedPaper`. Never raises: a file the
    harvest is still writing yields a partial paper and a populated
    ``parse_warnings``.

:mod:`syn4rch.curate.passages`
    Split a paper into :class:`~syn4rch.curate.passages.Passage` objects, each
    carrying a ``source_locator`` in the style SCHEMA.md 3.6 requires
    (``"Results para 4"``, ``"Fig 3C caption"``). Owns the locator grammar in
    both directions: it writes locators and it resolves them.

:mod:`syn4rch.curate.units`
    Build :class:`~syn4rch.curate.units.ExtractionUnitCandidate` objects --
    localized (population, method, outcome) co-occurrences. **These are
    candidates for a human or LLM reader. They are not records and assert
    nothing.** Matching reuses :class:`syn4rch.harvest.screen.TermMatcher`
    unchanged.

:mod:`syn4rch.curate.conditions`
    Read the condition dimension -- ``sex``, ``diet``, ``age_group``,
    ``challenge`` -- out of a paper's own text, with a verbatim quote and a
    locator for every value returned. Owns the distinction between a dimension
    a paper holds *constant* (admissible for every record from it) and one that
    *varies* within the paper (admissible for none), and the named refusals
    that keep a confident wrong value from being produced where a refusal
    belongs (D-038, D-040, D-041, D-044).

:mod:`syn4rch.curate.fill_conditions`
    ``python3 -m syn4rch.curate.fill_conditions run``: applies the above to the
    curated table, attaches the evidence, re-identifies what it changed and
    leaves an old-to-new map. The second writer to touch the dataset (D-031).

:mod:`syn4rch.curate.verify`
    The verbatim gate. ``verify_support`` confirms a record's
    ``support_quote`` occurs in the paper's retrieved text under one fixed
    normalisation and that its ``source_locator`` names a passage that exists.
    ``gate`` returns ``(admitted, rejected)``. Nothing here repairs a record.

:mod:`syn4rch.curate.run`
    ``python3 -m syn4rch.curate.run units`` and ``... gate``.

What this package does *not* contain, deliberately: the stage that reads a
paper and proposes a record (stage 6). That stage needs reading comprehension,
is therefore the only stage that can fabricate, and is kept outside the
deterministic code so that nothing here can be mistaken for having done it.

The Tier A schema is defined once, in :mod:`syn4rch.schema`. This package
imports it and never restates it.

Python 3.10 compatible, standard library only. No network.
"""

from .jats import (  # noqa: F401
    SECTION_LABELS,
    Caption,
    ParsedPaper,
    Section,
    classify_section,
    flatten_text,
    parse_jats,
    parse_jats_file,
)
from .passages import (  # noqa: F401
    KIND_FIGURE_CAPTION,
    KIND_PARAGRAPH,
    KIND_TABLE_CAPTION,
    KIND_TITLE,
    LocatorRef,
    Passage,
    PassageIndex,
    build_index,
    parse_locator,
    split_passages,
)
from .units import (  # noqa: F401
    SCORE_WEIGHTS,
    ExtractionUnitCandidate,
    PassageMatches,
    match_passage,
    paper_term_summary,
    score_unit,
    units_from_paper,
    units_from_papers,
)
from .verify import (  # noqa: F401
    CODES,
    NORMALISATION,
    Rejection,
    SupportVerification,
    VerificationFinding,
    find_paper,
    gate,
    index_papers,
    normalize_for_match,
    quote_occurs_in,
    verify_support,
)

from .conditions import (  # noqa: F401
    DIMENSIONS,
    REFUSALS,
    RULES,
    VERDICTS,
    Extracted,
    PaperConditions,
    paper_text,
    read_paper_conditions,
    read_record_conditions,
    sentences,
    verify_against_paper,
)

__all__ = [
    # jats
    "SECTION_LABELS",
    "Caption",
    "ParsedPaper",
    "Section",
    "classify_section",
    "flatten_text",
    "parse_jats",
    "parse_jats_file",
    # passages
    "KIND_FIGURE_CAPTION",
    "KIND_PARAGRAPH",
    "KIND_TABLE_CAPTION",
    "KIND_TITLE",
    "LocatorRef",
    "Passage",
    "PassageIndex",
    "build_index",
    "parse_locator",
    "split_passages",
    # conditions
    "DIMENSIONS",
    "REFUSALS",
    "RULES",
    "VERDICTS",
    "Extracted",
    "PaperConditions",
    "paper_text",
    "read_paper_conditions",
    "read_record_conditions",
    "sentences",
    "verify_against_paper",
    # units
    "SCORE_WEIGHTS",
    "ExtractionUnitCandidate",
    "PassageMatches",
    "match_passage",
    "paper_term_summary",
    "score_unit",
    "units_from_paper",
    "units_from_papers",
    # verify
    "CODES",
    "NORMALISATION",
    "Rejection",
    "SupportVerification",
    "VerificationFinding",
    "find_paper",
    "gate",
    "index_papers",
    "normalize_for_match",
    "quote_occurs_in",
    "verify_support",
]
