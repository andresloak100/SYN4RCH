"""Reading the condition dimension out of a paper's own text.

The condition dimension is ``(sex, diet, age_group, challenge)`` (SCHEMA.md
3.3b). This module proposes values for it, and every value it proposes carries
a verbatim sentence and a locator. It proposes nothing else: there is no
default, no fallback and no "most likely" value anywhere below. A dimension
the text does not state comes back as ``absent``, which is a result.

The one distinction this module is built around
-----------------------------------------------

A dimension a paper holds **constant** across all of its experiments can be
read from the Methods and applied to every record from that paper: if the only
animals in the paper are male, then the animals behind any one of its
measurements were male, and the Methods sentence saying so is evidence for
every record.

A dimension that **varies** within the paper cannot be assigned from the
Methods at all. A paper with a chow group and a high-fat group states both
diets in its Methods and states nowhere in its Methods which cohort any
individual record came from. Reading "high-fat diet" out of such a Methods and
stamping it on every record of the paper would be fabrication: the sentence is
real, the attribution is invented. For a varying dimension the only admissible
evidence is the record's own passage, which is why
:func:`read_record_conditions` exists and why it is not given the paper.

So :func:`read_paper_conditions` returns a *verdict* per dimension --
``"constant"``, ``"varies"`` or ``"absent"`` -- and an extraction only for the
constant ones. ``varies`` is not a failure. It is the module declining to
answer a question the Methods cannot answer, and on this corpus it is the
answer for ``diet`` in every diet-comparison study and for ``challenge`` in
very nearly all of them.

Why a tolerance test's challenge is not extracted
-------------------------------------------------

``GTT``, ``ITT``, ``ipGTT``, ``glucose tolerance test`` and friends are assay
names, and the acute state of a tolerance test is a function of the outcome
itself: every record whose ``outcome_id`` is a GTT readout was taken under a
glucose load, by definition of the assay. Writing ``glucose_load`` into
``challenge`` for those records therefore discriminates nothing that
``outcome_id`` does not already discriminate, while adding one more filled
condition field per record and so inflating the apparent condition coverage.
Coverage bought that way measures the extractor, not the corpus. The refusal
is named ``refused_assay_implies_its_own_challenge`` and it is reported, not
silent.

Why ``challenge`` is never read from the Methods
------------------------------------------------

The first version of this module had a Methods rule for ``challenge`` that
fired when the Methods named exactly one acute state. It refused husbandry
sentences by name, and then routed around its own refusal: in
``doi:10.1038/s41467-018-05711-6`` it read "Ad libitum fed rats were
anesthetized by an i.p. injection of ketamine" -- a *surgical preparation*
sentence, not husbandry -- and reported ``challenge = ad_libitum`` as a
constant of the paper. The animals were ad libitum fed when they were
anaesthetised for stereotaxic injection, which says nothing about the state
any animal was in during any measurement.

The general rule is the one that survives: **a Methods section cannot say
which acute state a given measurement was made in**, whatever the sentence's
wording. Husbandry ("ad libitum access to water and standard chow") is one
case of it; surgical preparation is another; a Methods that describes a
fasting protocol does not say which of the paper's measurements were taken
under it. So no Methods sentence is admissible evidence for ``challenge``, the
paper-level verdict is always ``absent`` or ``varies``, and the refusal is
reported as ``refused_methods_cannot_state_a_measurement_state``. The housing
case keeps its own name,
``refused_housing_statement_is_not_a_challenge``, because stamping
``ad_libitum`` on every record of such a paper would be fabrication at scale
and would look like a 100% coverage success.

``challenge`` is therefore read only from the record's own fields, by
:func:`read_record_conditions`.

Why an age at the start of a protocol is not the animals' age
-------------------------------------------------------------

"C57BL/6 mice were fed with HFD for 12 consecutive weeks (starting at 5 weeks
of age)" states an age of 5 weeks and describes animals that were about 17
weeks old when they were measured. ``age_group`` is the animal's band *when
the readout was taken*, so reading ``juvenile`` out of that sentence is wrong
in substance and wrong in the direction that matters. An age statement whose
sentence also carries a protocol-start or protocol-duration cue is refused:
``refused_protocol_start_age_is_not_measurement_age``.

The same reading refuses an age stated at acquisition ("12 weeks on arrival"),
at surgery ("7-9-month-old mice were anesthetized and placed in a stereotaxic
frame"), at weaning, at an acclimation or housing transition, at perfusion or
slice preparation (``refused_tissue_preparation_age``), or for one named assay
("For acute cold exposure, 8-week-old mice..."). None of those is the age at
the readout, and the gap between them and the readout is weeks the paper does
not state. This is expensive -- it takes this corpus from 14 ``age_group``
constants to 5 -- and the five that remain are the sentences that say what the
whole paper did ("All mice were aged 10-18 weeks at the time when the
experiments were performed"). See "Why refusing costs less than guessing"
below for why that is the right side of the trade.

The asymmetry with ``sex`` is deliberate: the *same* surgery sentence is
admissible evidence for sex and inadmissible for age, because an animal's sex
does not change between surgery and readout and its age does.

A refusal is also never a licence to promote what survived it. The same paper
states "8-10-week-old ... male Sprague Dawley rats and C57/BL6 WT mice", which
the rat rule skips; the 5-week statement was then the only one left and looked
unanimous. If **any** age statement in the Methods was refused for any reason,
the verdict is not ``constant``: the paper states an age the rule cannot use,
which is a reason to decline the dimension rather than to fall back on the
remainder.

Why refusing costs less than guessing
-------------------------------------

A wrong condition value is worse than a missing one, and the two costs are not
symmetric. The condition dimension decides **comparability** (SCHEMA.md 3.3b):
two records whose condition differs are declared incomparable. A ``null``
pools with every other unknown (``syn4rch.qc.keys.condition_key``), so two
records that should match still match. A *wrong* value actively splits them --
it manufactures a non-comparability that no experiment supports -- and because
``sex``, ``diet``, ``age_group`` and ``challenge`` are all ``record_id``
components (``RECORD_ID_FIELDS``), it does so under a new identity.

So the cost of refusing is coverage, and the cost of accepting is fabricated
distinctions plus a re-identified record. Five well-founded ``age_group``
constants beat fourteen with nine bad ones, and that is the trade every
refusal below is making.

Every rule here fails closed
----------------------------

A rule that cannot classify something must refuse, never guess, because every
leak this module has had was a fail-*open* one: an unlisted diet abbreviation
made one arm of a two-arm study invisible and the other arm was promoted to a
constant; an unlisted compound ("pair-fed", which is the *negation* of ad
libitum access) matched on the ``fed`` inside it; an unhandled qualifier
("below thermoneutrality", which is a cold challenge) matched the term it
inverts. In each case the rule produced a confident wrong value where "I do
not know" was the available and correct answer.

So the diet rule scans for diet-*shaped* tokens independently of its own
synonym table and refuses on any it cannot classify; the state table excludes
restricted-feeding compounds and qualifier-inverted terms by name; and a
refusal never promotes whatever survived it.

The one exception is in :func:`read_record_conditions`, where a state-like
token excluded by a trap alongside a *different* token that matched yields the
value with the rule ``..._with_a_competing_token`` and a note naming both.
``challenge`` is the scarcest dimension in this corpus and the surviving match
is usually right, so the uncertainty is published for a human to check rather
than the value dropped or the competition hidden.

**Known limitation, stated rather than papered over.** The record-level state
table (:data:`_STATE_TABLE`) has no fail-closed scanner, and the diet rule
does. An acute state this corpus names in words the table does not list -- a
hyperinsulinaemic clamp, "glucose-loaded", a hypoglycaemic challenge -- yields
no extraction *silently*, where an unclassified diet-shaped term yields a
refusal. That is the same shape as the unrecognised-diet-abbreviation defect,
one step short of being fixed: fixing it needs a vocabulary of what a
state-shaped phrase even looks like, and guessing that vocabulary here would
be the fabrication this module exists to prevent. The consequence is that a
missing ``challenge`` means either "the record says nothing" or "the record
says something this table cannot read", and the two are not distinguishable
from the output.

The verbatim gate
-----------------

:func:`verify_against_paper` re-checks every proposed quote against the
paper's own retrieved text with :func:`syn4rch.curate.verify.quote_occurs_in`,
and :func:`read_paper_conditions` drops any extraction that fails. This is
mechanical: a rule below could be outright buggy and still could not emit a
sentence the paper does not contain.

Nothing here writes a record, a file or a ``not_reported``. ``not_reported``
asserts an absence and an absence cannot be quoted (see
:class:`syn4rch.schema.record.ConditionEvidence`); where a rule would produce
it, this module produces no extraction instead.

Python 3.10 compatible, standard library only (D-001).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from ..schema.enums import AgeGroup, AgeUnit, Challenge, Diet, Sex, Species
from ..schema.rules import age_in_weeks, assign_age_group
from .jats import ParsedPaper
from .verify import quote_occurs_in

__all__ = [
    "DIMENSIONS",
    "VERDICTS",
    "RULES",
    "REFUSALS",
    "Extracted",
    "PaperConditions",
    "paper_text",
    "inadmissible_methods_conditions",
    "read_paper_conditions",
    "read_record_conditions",
    "sentences",
    "verify_against_paper",
]


#: The dimensions this module reads, in a locally fixed iteration order.
#:
#: Equal to :data:`syn4rch.schema.record.CONDITION_DIMENSIONS` (D-075).
#:
#: It was a **subset** from D-060 until D-075, excluding ``background_genotype``
#: on the argument that deciding whether an allele is the recorded perturbation
#: or the substrate it was applied to is a per-record judgement this module
#: cannot make. That argument is sound and it is about **extraction**. Being out
#: of this tuple also took the dimension out of the **verdict**, and so out of
#: :func:`inadmissible_methods_conditions` -- which needs only the verdict. One
#: tuple was carrying two decisions, and the second one was never argued.
#:
#: So genotype is in, with no Methods extraction rule, exactly as ``challenge``
#: is (see :func:`_read_genotype`). A dimension this module cannot extract is
#: still a dimension it can refuse, and refusing is the whole of what 3b asks.
DIMENSIONS = ("sex", "diet", "age_group", "challenge", "background_genotype")

#: The three verdicts a dimension can earn. See the module docstring.
VERDICTS = ("constant", "varies", "absent")

#: Every rule that can produce an extraction. A value in the corpus whose
#: ``rule`` is not in this tuple did not come from this module.
RULES = (
    "methods_single_sex",
    "methods_both_sexes",
    "methods_single_diet",
    "methods_single_age_band",
    # There is deliberately no Methods rule for ``challenge``. See the module
    # docstring: no Methods sentence can say which acute state a given
    # measurement was made in.
    "outcome_names_the_state",
    "quote_names_the_state",
    # D-044's one exception: a value kept with its uncertainty published
    # rather than refused. See :func:`read_record_conditions`.
    "outcome_names_the_state_with_a_competing_token",
    "quote_names_the_state_with_a_competing_token",
)

#: Every named refusal. These are reported in ``notes`` rather than being
#: resolved into a value, because the point of each one is that the honest
#: answer is "no extraction" and a silent refusal is indistinguishable from a
#: rule that never ran.
REFUSALS = (
    "refused_methods_cannot_state_a_measurement_state",
    "refused_figure_scoped_cohort_sentence",
    "refused_disjunctive_sex_is_unresolved",
    "refused_second_sex_reported_outside_the_methods",
    "refused_ambiguous_sex_sentence",
    "refused_unrecognised_diet_term",
    "refused_disjunctive_diet_is_two_arms",
    "refused_housing_statement_is_not_a_challenge",
    "refused_assay_implies_its_own_challenge",
    "refused_two_states_in_one_outcome",
    "refused_quote_names_more_than_one_state",
    "refused_restricted_feeding_has_no_enum_member",
    "refused_qualifier_inverts_the_state",
    "refused_fed_inside_a_diet_phrase",
    "refused_sex_not_constant_after_a_refusal",
    "refused_rat_age_band",
    "refused_age_band_requires_mouse",
    "refused_straddling_age_range",
    "refused_in_vitro_age_statement",
    "refused_protocol_start_age_is_not_measurement_age",
    "refused_tissue_preparation_age",
    "refused_age_not_constant_after_a_refusal",
    "refused_record_has_no_locator",
    "dropped_quote_not_verbatim",
)


# --------------------------------------------------------------------------- #
# public data model
# --------------------------------------------------------------------------- #


@dataclass
class Extracted:
    """One proposed condition value, with the evidence that produced it.

    Maps one-to-one onto :class:`syn4rch.schema.record.ConditionEvidence`:
    ``source`` is a :class:`~syn4rch.schema.enums.ConditionSource` value
    string and ``value`` is the target enum's value string, so a caller can
    build the evidence entry without translating anything.
    """

    dimension: str
    value: str
    source: str
    locator: str
    quote: str
    rule: str


@dataclass
class PaperConditions:
    """What one paper's Methods can and cannot say about its own conditions."""

    study_id: str
    pmcid: str
    #: dimension -> "constant" | "varies" | "absent"
    verdicts: Dict[str, str]
    #: dimension -> Extracted, only for dimensions whose verdict is "constant"
    constants: Dict[str, "Extracted"]
    #: dimension -> the refusal names recorded while reading it (D-079).
    #:
    #: Attribution is by *which reader was called*, never by parsing a refusal
    #: name. `REFUSALS` is a flat tuple with no dimension in it, and two names
    #: were already missing from it, so a hand-authored name->dimension map
    #: would be one more declared thing that drifts from the thing it
    #: describes. The caller already knows the dimension; this records what it
    #: knew.
    refusals: Dict[str, List[str]] = field(default_factory=dict)
    #: what the Methods said about the animals, for the age rule
    species: Optional[str] = None
    age_value: Optional[float] = None
    age_unit: Optional[str] = None
    notes: List[str] = field(default_factory=list)


# --------------------------------------------------------------------------- #
# sentence splitting
# --------------------------------------------------------------------------- #

#: Tokens that end in a period without ending a sentence. Deliberately short:
#: a missed split merges two sentences, which costs precision in the quote but
#: cannot invent a value, while an over-eager split can cut a quote in half
#: and is therefore the worse failure.
_ABBREVIATIONS = frozenset(
    """
    e.g i.e etc vs cf ca approx incl resp fig figs eq eqs ref refs no nos
    al dr prof mr mrs ms st jr sr inc ltd co dept univ
    i.p s.c i.v i.c.v p.o o.n
    min sec hr hrs h wk wks mo yr yrs
    """.split()
)

#: End-of-sentence candidate: terminal punctuation, optional closing quotes or
#: brackets, then whitespace.
_SENTENCE_BREAK = re.compile(r"([.!?])([\"')\]’”]*)(\s+)")

#: Word immediately before a candidate period.
_TRAILING_WORD = re.compile(r"([A-Za-z](?:[A-Za-z.]*[A-Za-z])?|[A-Za-z])\.?$")


def _is_abbreviation(text: str, dot: int) -> bool:
    """True when the period at ``dot`` belongs to an abbreviation, not a stop."""
    match = _TRAILING_WORD.search(text[:dot])
    if match is None:
        return False
    token = match.group(1).rstrip(".").lower()
    if not token:
        return False
    if token in _ABBREVIATIONS:
        return True
    # A single letter before a period is an initial ("Smith, J. et al."), and
    # "C57BL/6J." is not an initial because the token is longer than one char.
    return len(token) == 1


def sentences(text: str) -> List[str]:
    """Split ``text`` into sentences, conservatively.

    Every returned sentence is a contiguous substring of ``text`` with only
    surrounding whitespace removed. Nothing is reordered, substituted or
    re-punctuated, so a sentence returned here still verifies verbatim against
    the text it came from (the gate's normalisation collapses whitespace and
    strips the ends, which is exactly the difference this function can
    introduce).

    A period is treated as a sentence end only when it is followed by
    whitespace and then something that can start a sentence, and only when the
    word before it is not a known abbreviation or a single-letter initial. So
    ``"p < 0.01"``, ``"8.5 weeks"``, ``"e.g. mice"`` and ``"Smith, J. et al."``
    do not split; ``"...were used. Mice were..."`` does.
    """
    if not text:
        return []
    body = str(text)
    out: List[str] = []
    start = 0
    for match in _SENTENCE_BREAK.finditer(body):
        if match.start(1) < start:
            continue
        if _is_abbreviation(body, match.start(1)):
            continue
        nxt = body[match.end(3) : match.end(3) + 1]
        if nxt and not (
            nxt.isupper() or nxt.isdigit() or nxt in "(“‘\"'["
        ):
            continue
        piece = body[start : match.end(2)].strip()
        if piece:
            out.append(piece)
        start = match.end(3)
    tail = body[start:].strip()
    if tail:
        out.append(tail)
    return out


# --------------------------------------------------------------------------- #
# shared vocabulary
# --------------------------------------------------------------------------- #

#: Dash code points that appear interchangeably in publisher XML. Ranges in
#: this corpus use ASCII hyphens, en dashes and em dashes without distinction,
#: so every pattern below accepts any of them.
#: The ASCII hyphen is written last so that every ``[...]`` class built from
#: this string ends with it, where a hyphen is a literal and cannot open a
#: character range.
_DASHES = "‐‑‒–—―−-"
_DASH = "[" + _DASHES + "]"
#: Separator inside a multi-word phrase: whitespace, a dash, or both.
_GAP = "[\\s" + _DASHES + "]+"


def _phrase(text: str) -> str:
    """Regex source for a phrase whose internal gaps may be spaces or dashes."""
    return r"\b" + _GAP.join(re.escape(part) for part in text.split()) + r"\b"


def _any(*patterns: str) -> "re.Pattern[str]":
    return re.compile("|".join(patterns), re.IGNORECASE)


#: Numbers this corpus spells out. "Twelve-week-old C57BL/6 mice" and
#: "Three-week-old female breeders" are both real wordings, and reading the
#: word "twelve" as 12 is parsing, not inference: the sentence says twelve.
#: Stopping at twenty is where the forms stop appearing.
_WORD_NUMERALS = {
    "one": 1.0, "two": 2.0, "three": 3.0, "four": 4.0, "five": 5.0,
    "six": 6.0, "seven": 7.0, "eight": 8.0, "nine": 9.0, "ten": 10.0,
    "eleven": 11.0, "twelve": 12.0, "thirteen": 13.0, "fourteen": 14.0,
    "fifteen": 15.0, "sixteen": 16.0, "seventeen": 17.0, "eighteen": 18.0,
    "nineteen": 19.0, "twenty": 20.0,
}

#: Longest first, so that "seventeen" is read before the "seven" inside it.
_WORD_NUMERAL_ALTERNATION = "|".join(
    sorted(_WORD_NUMERALS, key=len, reverse=True)
)

#: A number, in digits or in words.
_NUM = r"(?:\d+(?:\.\d+)?|" + _WORD_NUMERAL_ALTERNATION + ")"


def _as_number(text: str):
    """``float`` for a digit string or a spelled-out numeral, else ``None``."""
    if text is None:
        return None
    try:
        return float(text)
    except (TypeError, ValueError):
        return _WORD_NUMERALS.get(str(text).strip().lower())


#: Sentences without one of these do not describe the experimental animals, and
#: a sentence that does not describe the animals cannot state their sex or age.
_ANIMAL_CUE = _any(
    r"\bmice\b", r"\bmouse\b", r"\brats?\b", r"\banimals\b",
    r"C57BL", r"\bSprague\b", r"\bWistar\b",
)

#: Mouse cues and rat cues, kept apart because ``assign_age_group`` is
#: mouse-only by design (SCHEMA.md 3.3) and refuses to extrapolate.
_MOUSE_CUE = _any(r"\bmice\b", r"\bmouse\b", r"C57BL")
_RAT_CUE = _any(
    r"\brats?\b", r"\bSprague\b", r"\bWistar\b", r"\bLong" + _DASH + r"?\s?Evans\b"
)

#: A colony's female breeders say nothing about the animals that were measured,
#: so a sentence about breeding casts no sex vote. ``dam`` and ``pup`` are
#: matched as whole words: ``dam`` as a substring also matches "damage" and
#: ``pup`` also matches "pupil", and excluding a real cohort sentence on that
#: basis would be its own silent error.
#:
#: ``litter`` is matched as a whole word and **not** as a substring, because
#: "littermate" is a cohort word, not a breeding word: "male littermates served
#: as controls" is exactly the kind of sentence that states the cohort's sex,
#: and excluding it as breeding would discard real evidence.
_BREEDING_CUE = _any(
    r"breed", r"\bdams?\b", r"pregnan", r"\blitters?\b", r"\bpups?\b",
    r"sperm", r"oocyte", r"embryo",
)

#: Cues that an age names the *start* or the *span* of a protocol rather than
#: the animals' age when they were measured. "starting at 5 weeks of age" in a
#: sentence about a twelve-week diet describes a 17-week-old animal at readout,
#: and ``age_group`` is the band at readout. Any of these in the sentence
#: refuses the statement; the refusal costs coverage on sentences that happen
#: to mention a duration, which is the right direction to be wrong in.
_PROTOCOL_AGE_CUE = _any(
    r"\bstarting\b",
    r"\bstart(?:ed|s|ing)?\s+(?:at|from|on)\b",
    r"\bbeginning\b",
    r"\bbegan\b",
    r"\bbegins\b",
    r"\bonwards?\b",
    r"\bthereafter\b",
    r"\buntil\b",
    r"\bfollowed by\b",
    r"\bfor\s+(?:a\s+(?:further|total)\s+(?:of\s+)?)?" + _NUM + r"\s*"
    r"(?:more\s+|additional\s+|consecutive\s+)?(?:days?|weeks?|months?|years?)\b",
    r"\bfrom\s+" + _NUM + r"\s*(?:days?|weeks?|months?|years?)\s+of\s+age\b",
    r"\bover\s+(?:a\s+period|" + _NUM + r"\s*(?:days?|weeks?|months?))\b",
    r"\bduring\s+" + _NUM + r"\s*(?:days?|weeks?|months?)\b",
    r"\bat\s+the\s+age\s+of\b",
    r"\bwean(?:ing|ed)\b",
    # Acquisition. An age on the invoice is the age on arrival, and the
    # readout comes after the acclimation that the next sentence describes.
    r"\b(?:obtained|purchased|received|sourced|supplied|bred)\s+(?:from|by)\b",
    r"\barrived\b", r"\binitial\s+body\s+weight\b",
    # Surgery. The age at stereotaxic injection or implantation is the age
    # weeks before the readout, while the virus expresses.
    r"\ban(?:a)?esthet(?:i[sz]ed|ic|ised|ized)\b", r"\bstereotax",
    r"\bwere\s+(?:implanted|operated|injected\s+bilaterally|cannulated)\b",
    r"\bsurgery\b", r"\bsurgical",
    # A housing or acclimation transition stated with an age is the age at
    # which the protocol began.
    r"\bseparated\s+into\b", r"\b(?:singly|individually)\s+housed\b",
    r"\bsingle\s+cages?\b", r"\bacclimat", r"\bhabituat",
)

#: Terminal-procedure and tissue-preparation cues. An age stated in a sentence
#: about perfusion, slice preparation, decapitation or tissue collection is the
#: age of the animals used for *that* procedure, which is routinely a separate
#: cohort from the one a given record was measured in -- the same shape as the
#: surgical-preparation sentence that made the Methods ``challenge`` rule
#: untenable, and as the protocol-start age above.
#:
#: Found while auditing ``doi:10.1016/j.celrep.2021.109997``, whose only age
#: statement is "3- to 5-week-old male or female mice were deeply
#: anesthetized ... and trans-cardially perfused". ``juvenile`` would have been
#: stamped on every record of that paper, including its adult recordings.
_TISSUE_PREPARATION_CUE = _any(
    r"\bperfus", r"\btrans" + _DASH + r"?cardial", r"\bdecapitat",
    r"\bslices?\s+were\s+(?:prepared|cut|obtained)\b",
    r"\bslice\s+preparation\b", r"\bsacrific", r"\beuthani[sz]",
    r"\bbrains?\s+were\s+(?:removed|collected|extracted|dissected)\b",
    r"\btissues?\s+(?:was|were)\s+(?:collected|harvested|fixed)\b",
    r"\bcryostat\b", r"\bimmunohistochem", r"\bin situ hybridi",
    r"\bpost" + _DASH + r"?mortem\b",
)

#: An age given for cultured cells is not the animals' age.
_IN_VITRO_CUE = _any(
    r"\bin vitro\b", r"\bcultur", r"\bplat(?:ed|ing)\b", r"\bDIV\b",
    r"\bcell lines?\b", r"\bpassage[sd]?\b", r"\bdissociated\b",
)


#: Assay classes: a qualifier naming ONE kind of experiment among several.
#:
#: "All behavioural experiments ... were done on male mice" scopes itself to a
#: subset; "All experiments were performed on male mice" does not. A
#: word-presence test cannot tell them apart, because "All animal experiments"
#: and "All the mouse experiments in this study" carry qualifiers too and those
#: qualify the *subject* rather than selecting a subset of assays. The
#: distinction is carried by this lexicon, so ``animal``, ``mouse``, ``mice``,
#: ``rodent`` and ``this`` are deliberately absent from it.
#:
#: **Unlike the genotype guard (D-075), a miss here changes the outcome.** That
#: guard was safe to leave modest because both of its verdicts refused; an
#: assay class this lexicon does not know leaves a wrong ``constant`` standing.
#: So the lexicon is generous and the cost of that is measured in D-078 rather
#: than assumed -- 974 papers re-read before and after.
_ASSAY_CLASS = (
    r"behaviou?ral|electrophysiolog\w*|ephys|optogenetic|chemogenetic"
    r"|photometr\w*|imaging|histolog\w*|immunohisto\w*|patch[\s" + _DASHES + r"]?clamp"
    r"|slice|calcium|fibre|fiber|EEG|EMG|telemetr\w*|multichannel"
    r"|microdialysis|tracing|anatomical|vitro|vivo|qPCR|western"
    r"|perfusion|surgic\w*|feeding|locomotor"
)

#: The nouns an assay class attaches to.
_ASSAY_NOUN = (
    r"experiments?|recordings?|studies|assays?|measurements?"
    r"|sessions?|tests?|manipulations?"
)

#: A figure, panel, cohort or experiment reference scopes a sentence to a
#: *subset* of the paper, so it cannot vote on a paper-level verdict.
#: "Experiments (Figure 3) were performed in female mice only" says the
#: opposite of "this paper used female mice": it says there is a female-only
#: subset, and says nothing about the sex used in the paper's other figures.
_SCOPED_TO_A_SUBSET = _any(
    r"\bfig(?:ure|s|\.)?\s*\d",
    r"\bsupplementary\s+fig",
    r"\bextended\s+data\s+fig",
    r"\bpanels?\s*[A-Za-z0-9]",
    r"\bfor\s+(?:the\s+)?experiments?\s+(?:in|of|shown)\b",
    r"\bexperiments?\s+(?:in|shown\s+in|of)\s+fig",
    # "For the diet-induced obesity experiments, ...", "For acute cold
    # exposure, ..." -- a sentence that opens by naming one assay describes
    # that assay's cohort, not the paper's.
    r"\bfor\s+the\s+[\w" + _DASHES + r"]+(?:\s+[\w" + _DASHES + r"]+){0,2}"
    r"\s+experiments?\b",
    r"\bfor\s+(?:acute|chronic)\s+\w+",
    r"\bin\s+(?:this|these|that|those|one|a\s+separate|a\s+subset\s+of|the\s+second)"
    r"\s+(?:cohort|cohorts|group|groups|subset|experiment|experiments|"
    r"set\s+of\s+experiments)\b",
    r"\bcohort\s*\d",
    # A direct subset marker: "performed in female mice only".
    r"\b(?:males?|females?)\s+(?:mice|rats|animals)\s+only\b",
    # D-078. A universal quantifier over ONE assay class is not the paper:
    # "All behavioural experiments conducted to target Hcrt neurons were done
    # on ... male mice" left 20 of 27 records -- corticosterone, blood
    # pressure, heart rate -- outside the sentence that spoke for them.
    r"\b(?:all|every)\s+(?:the\s+)?(?:" + _ASSAY_CLASS + r")[\w"
    + _DASHES + r"]*\s+(?:" + _ASSAY_NOUN + r")\b",
    # "For chemogenetic experiments, only male mice were used"; "were used for
    # electrophysiological recordings". The form above this one required the
    # literal "the" and the literal word "experiments", so both missed.
    r"\b(?:for|in|during)\s+(?:the\s+|our\s+)?(?:" + _ASSAY_CLASS
    + r")[\w" + _DASHES + r"]*\s+(?:" + _ASSAY_NOUN + r")\b",
    # The enumerated-assay subject: "Optogenetic manipulations, fiber
    # photometry, in vivo multichannel recording, and RNA-inference
    # experiments were performed in male ... mice".
    r"^\s*(?:" + _ASSAY_CLASS + r")[\w" + _DASHES + r"]*\s+(?:"
    + _ASSAY_NOUN + r")\b",
    # An explicit denial of universality. "In most cases, behavioural
    # experiments were with male mice" says there are other cases.
    r"\bin\s+most\s+cases\b",
    r"\bmost\s+(?:of\s+the\s+)?(?:" + _ASSAY_NOUN + r")\b",
    # "...were used for most experiments, except for optogenetic experiments
    # that used males only" -- the exception names the second cohort.
    r"\bexcept\s+(?:for\s+)?[^.;]{0,80}?\b(?:only|alone)\b",
    # "Sixteen of the Orexin-tTA mice were reared with chow that included
    # doxycycline" -- a counted subset of the colony, not the colony.
    r"\b(?:\d+|" + "|".join(_WORD_NUMERALS) + r")\s+of\s+"
    r"(?:the\s+|these\s+|those\s+)?[\w" + _DASHES + r"]*\s*"
    r"(?:mice|mouse|rats?|animals)\b",
)


# --------------------------------------------------------------------------- #
# sex
# --------------------------------------------------------------------------- #

_MALE = re.compile(r"\bmales?\b", re.IGNORECASE)
_FEMALE = re.compile(r"\bfemales?\b", re.IGNORECASE)
#: An explicit conjunction of sexes. ``mixed-sex`` joins "both sexes" rather
#: than "either sex": it asserts that both were used, which is what
#: :class:`Sex.both` means.
_BOTH_SEXES = re.compile(
    r"\bboth\s+sexes\b|\bboth\s+genders\b|\bmixed[\s" + _DASHES + r"]sex\b",
    re.IGNORECASE,
)

#: "mice of either sex" leaves an individual record's sex unresolved, and a
#: ``constant`` verdict claims that no second value exists.
#:
#: The phrase contains neither a ``male`` nor a ``female`` token, so
#: :data:`_SEX_DISJUNCTION` never saw it, and a sentence naming a male cohort
#: *and* an either-sex cohort -- "Male mice at 10-16 weeks were used for
#: behavioral experiments, whereas mice of either sex at 4-8 weeks were used
#: for in vitro patch-clamp experiments" -- voted ``male`` with no refusal.
#: The sentence stored as the claim that no second value exists stated the
#: second value verbatim (D-078).
_EITHER_SEX = re.compile(
    r"\b(?:either|any)\s+(?:sex|sexes|gender|genders)\b"
    r"|\bsex(?:es)?\s+(?:was|were)\s+not\s+"
    r"(?:specified|determined|recorded|reported)\b",
    re.IGNORECASE,
)

_SEX_TOKEN = r"(?:males?|females?)"
#: An optional parenthetical between the two sex words, so that
#: "male (n = 5) and female (n = 5) mice" reads as the conjunction it is.
_PARENTHETICAL = r"(?:\([^)]*\)\s*)?"

#: Only an explicit **conjunction** may vote ``both``.
_SEX_CONJUNCTION = re.compile(
    r"\b" + _SEX_TOKEN + r"\b\s*" + _PARENTHETICAL + r"and\s+"
    + _PARENTHETICAL + r"\b" + _SEX_TOKEN + r"\b",
    re.IGNORECASE,
)
#: A **disjunction** states that either sex was used, which for any individual
#: record leaves the sex unresolved. "male or female mice were anesthetized"
#: is not ``both``: ``both`` asserts that both were used, and the sentence
#: asserts only that one of them was.
_SEX_DISJUNCTION = re.compile(
    r"\b" + _SEX_TOKEN + r"\b\s*" + _PARENTHETICAL + r"or\s+"
    + _PARENTHETICAL + r"\b" + _SEX_TOKEN + r"\b",
    re.IGNORECASE,
)


def _sex_vote(sentence: str) -> Tuple[Optional[Sex], Optional[str]]:
    """``(vote, refusal)`` for one Methods sentence.

    Exactly one of the two is ever set. ``(None, None)`` means the sentence
    says nothing about the cohort's sex; a refusal name means it says something
    this rule will not turn into a value.

    ``both`` requires an explicit conjunction ("male and female", "males and
    females", "both sexes"). Two sex words joined by "or", or two sex words
    with neither construction, leave the sex unresolved and are refused: a
    ``constant`` verdict claims that no second value exists in the paper, and
    neither sentence establishes that.
    """
    if not _ANIMAL_CUE.search(sentence):
        return None, None
    if _BREEDING_CUE.search(sentence):
        # A colony's female breeders say nothing about the animals measured,
        # so this is not a cohort sentence at all -- not a refusal.
        return None, None
    # Before the token test below, which cannot see a sexless phrasing.
    if _EITHER_SEX.search(sentence):
        return None, "refused_disjunctive_sex_is_unresolved"
    male = bool(_MALE.search(sentence))
    female = bool(_FEMALE.search(sentence))
    if not (male or female or _BOTH_SEXES.search(sentence)):
        return None, None
    if _SCOPED_TO_A_SUBSET.search(sentence):
        return None, "refused_figure_scoped_cohort_sentence"
    if _BOTH_SEXES.search(sentence):
        return Sex.both, None
    if male and female:
        if _SEX_CONJUNCTION.search(sentence):
            return Sex.both, None
        if _SEX_DISJUNCTION.search(sentence):
            return None, "refused_disjunctive_sex_is_unresolved"
        return None, "refused_ambiguous_sex_sentence"
    return (Sex.male if male else Sex.female), None


# --------------------------------------------------------------------------- #
# diet
# --------------------------------------------------------------------------- #

#: Modifiers that are determiners or back-references: they point at a diet
#: named somewhere else, so they are not evidence of an unknown diet type.
_DIET_MODIFIER_DETERMINER = frozenset(
    """
    a an the this that these those their its each same both either neither
    any all other another respective individual by between on in with of
    """.split()
)

#: Modifiers naming physical form or handling rather than composition.
_DIET_MODIFIER_PHYSICAL_FORM = frozenset(
    """
    pelleted pellet pelletised pelletized powdered powder liquid solid
    semi-solid semisolid granular paste mash extruded irradiated autoclaved
    sterile sterilised sterilized free inaccessible accessible available
    fresh dry dried wet gel gelled
    """.split()
)

#: Modifiers naming the species or purpose a feed is sold for. These are the
#: case-insensitive form of the product-name rule below: "PicoLab Rodent Diet
#: 20" and "a rodent diet (Harlan)" are the same kind of statement, and only
#: one of them is capitalised. None of these words says anything about
#: composition.
_DIET_MODIFIER_CATEGORY = frozenset(
    """
    rodent rodents mouse mice murine rat rats animal laboratory lab
    maintenance breeder breeding growth commercial
    """.split()
)

#: Modifiers naming a composition class. These are the real signal: a
#: composition the synonym table does not classify may be a second arm.
#: Listed for the record and for the report; the decision below does not
#: depend on the list being complete, because an unrecognised *lower-case*
#: modifier is treated as a signal either way.
_DIET_MODIFIER_COMPOSITION = frozenset(
    """
    fat sugar sucrose fructose glucose carbohydrate carbohydrates protein
    cholesterol calorie caloric energy western cafeteria ketogenic choline
    methionine fibre fiber lard obesogenic
    """.split()
)

#: The synonym table. ``chow`` on its own covers "standard chow", "normal
#: chow", "regular chow" and "normal chow diet"; the ``(standard|normal|
#: regular) diet`` branch covers the chow-free wordings. The abbreviations are
#: matched case-**sensitively**: lower-case "hfd" is not a word this
#: literature uses, and a case-insensitive ``\bSCD\b`` or ``\bND\b`` would
#: fire on unrelated acronyms.
#:
#: This table is **incomplete and always will be**, which is the whole reason
#: :func:`_diet_occurrences` also scans for diet-*shaped* tokens it cannot
#: classify. ``ND`` was missing from it once, one arm of a two-arm study went
#: invisible, and the surviving arm was promoted to a constant for every
#: record of that paper.
#: Words that may sit between "standard/normal/regular" and "diet" without
#: changing what is being said: "a regular rodent diet" and "a standard
#: pelleted diet" are both "regular diet" with a category or a physical form
#: inserted. Built from the modifier sets below so the two cannot drift apart.
_DIET_INTERVENING = r"(?:\s+(?:%s))?" % "|".join(
    sorted(_DIET_MODIFIER_CATEGORY | _DIET_MODIFIER_PHYSICAL_FORM)
)

_DIET_CLASSIFIERS = (
    (
        Diet.chow,
        _any(
            _phrase("chow"),
            r"\b(?:standard|normal|regular)" + _DIET_INTERVENING + r"\s+diet\b",
        ),
    ),
    (
        Diet.high_fat,
        _any(
            r"\bhigh" + _GAP + r"fat\s+diet\b",
            _phrase("western diet"),
            r"\b60\s*%\s*(?:of\s*)?(?:kcal|kJ|calories)?\s*(?:from\s+)?fat\b",
        ),
    ),
    (Diet.chow, re.compile(r"\b(?:NCD|SCD)\b")),
    (Diet.high_fat, re.compile(r"\bHFD\b")),
)

#: ``ND`` means "normal diet" in this literature and also means "not
#: determined" in plenty of other prose, so it is classified only in a
#: sentence that is talking about food at all. Out of that context it stays
#: unclassified and earns the refusal, which is the safe direction.
_DIET_ND = re.compile(r"\bND\b")
_DIET_CONTEXT = _any(
    r"\bdiets?\b", r"\bchow\b", r"\bfed\b", r"\bfeeding\b", r"\bkcal\b",
    r"\bfat\b", r"\bHFD\b", r"\bcalorie", r"\bnutrition",
)

#: Diet-*shaped* tokens, scanned independently of the synonym table above. An
#: occurrence of one of these that the table does not classify means the
#: lexicon is incomplete, and an incomplete lexicon must produce a refusal.
_DIET_SHAPES = (
    re.compile(r"\bchow\b", re.IGNORECASE),
    re.compile(r"\b(?:HFD|HFHS|ND|NCD|SCD|LFD|CD|RD|WD)\b"),
)

#: Every occurrence of the word "diet", with whatever modifies it. The
#: modifier may be attached by a space, by a dash, or not at all
#: ("LabDiet"). Whether the occurrence is a *signal* is then decided by
#: :func:`_diet_modifier_is_a_signal`, not by the bare fact that the word
#: appeared.
#:
#: The first version of this scanner was ``\b\w*\s?diet\b`` -- any word before
#: "diet" -- and it refused on 11 of the 42 papers over terms that said
#: nothing about composition: product names ("PicoLab Rodent Diet 20",
#: "LabDiet 5021"), English back-references ("this diet", "each diet"), and
#: physical form ("pelleted diet", "liquid diet"). Those are how a paper
#: *identifies* or *refers back to* its chow, not a second arm, and refusing
#: on them cost 11 real ``chow`` constants.
_DIET_WORD = re.compile(
    r"(?:(?P<mod>[\w%" + _DASHES + r"]+)[\s" + _DASHES + r"]*)?"
    r"(?P<diet>[Dd][Ii][Ee][Tt])\b"
)


#: The food words the stimulus guard screens, and the abbreviations it did not.
_FOOD_WORD = r"(?:chow|food|diets?)"
#: Matched case-sensitively; see the note in :data:`_DIET_AS_STIMULUS`.
_FOOD_ABBR = r"(?:HFD|HFHS|NCD|SCD|LFD|ND|CD|RD|WD)"

#: Food as a *stimulus* is not a maintenance diet. "response to chow" in a
#: neural-recording sentence is chow presented to the animal during the
#: recording; the paper may well be chow-fed, but this sentence is not
#: evidence of it. This is the ``diet`` analogue of the husbandry trap, and
#: unlike a refusal it removes the occurrence rather than declining the
#: dimension: a stimulus mention is no evidence of a *second* diet either.
_DIET_AS_STIMULUS = (
    _any(
        r"\b(?:response|responses|responsive|responding|respond(?:ed|s)?)\s+to\s+"
        r"(?:the\s+)?(?:\w+[\s" + _DASHES + r"]+){0,2}?(?:chow|food|diet)\b",
    ),
    _any(
        r"\b" + _FOOD_WORD + r"\s+(?:presentation|delivery|exposure)\b",
        r"\b" + _FOOD_WORD + r"\s+was\s+(?:presented|delivered|offered|introduced)\b",
        r"\bupon\s+(?:chow|food)\b",
        r"\bpresentation\s+of\s+(?:\w+[\s" + _DASHES + r"]+){0,2}?(?:chow|food)\b",
        r"\b(?:presented|delivered|offered)\s+(?:with\s+)?"
        r"(?:\w+[\s" + _DASHES + r"]+){0,2}?(?:chow|food)\b",
        # D-078. The food token BEFORE the verb: "The food presented in the
        # feeding experimental sessions was standard laboratory chow".
        r"\b" + _FOOD_WORD + r"\s+(?:that\s+was\s+)?"
        r"(?:presented|delivered|offered|introduced)\b",
        # "2-min exposure to salient stimuli including ... HFD (in a small
        # petri dish)". The guard's own cue list had "chow exposure" and not
        # "exposure to chow".
        r"\b(?:exposure|exposed)\s+to\s+(?:\w+[\s" + _DASHES + r"]+){0,6}?"
        + _FOOD_WORD + r"\b",
        r"\bstimuli\s+includ\w*\s+(?:\w+[\s" + _DASHES + r"]+){0,8}?"
        + _FOOD_WORD + r"\b",
        r"\bin\s+a\s+(?:small\s+)?petri\s+dish\b",
    ),
    # Case-SENSITIVE, because the diet abbreviations are: `\bCD\b` and
    # `\bND\b` under IGNORECASE would match English words. This is the half
    # the guard was missing entirely -- its food alternation was
    # `(?:chow|food|diet)` while `_DIET_SHAPES` also knows the abbreviations,
    # so a stimulus written "HFD" was invisible to it and wrote `high_fat`
    # onto 9 records from a behavioural assay (D-078).
    re.compile(
        r"(?:\b(?:[Ee]xposure|[Ee]xposed)\s+to\s+(?:\w+[\s" + _DASHES
        + r"]+){0,6}?" + _FOOD_ABBR + r"\b"
        r"|\bstimuli\s+includ\w*\s+(?:\w+[\s" + _DASHES + r"]+){0,8}?"
        + _FOOD_ABBR + r"\b"
        r"|\b" + _FOOD_ABBR + r"\s*\(\s*in\s+a\s+(?:small\s+)?petri"
        r"|\b" + _FOOD_ABBR + r"\s+(?:presentation|delivery|exposure)\b"
        r"|\b" + _FOOD_ABBR + r"\s+(?:(?:that\s+)?was\s+)?"
        r"(?:presented|delivered|offered|introduced)\b)"
    ),
)

_EITHER = re.compile(r"\beither\b", re.IGNORECASE)
_OR = re.compile(r"\bor\b", re.IGNORECASE)


@dataclass(frozen=True)
class _DietHit:
    """One diet-shaped occurrence. ``value`` is ``None`` when unclassified."""

    start: int
    end: int
    value: Optional[Diet]
    text: str


def _overlaps(start: int, end: int, spans) -> bool:
    return any(not (end <= s or start >= e) for s, e in spans)


def _diet_modifier_is_a_signal(modifier, diet_token, at_sentence_start) -> bool:
    """Does this "diet" occurrence mean the lexicon is missing a diet *type*?

    The question is about the modifier, not about the word "diet". An
    unclassified term should mean "the synonym table does not know this diet
    type", never "the word diet appeared next to something".

    Not a signal:

    * no modifier at all -- a bare "diet" refers to a diet named elsewhere;
    * a determiner or back-reference ("this diet", "each diet", "by diet");
    * physical form or handling ("pelleted diet", "irradiated diet");
    * a product name, recognised as a Capitalised or CamelCase modifier
      ("Rodent Diet", "LabDiet", "Special Diet") or a capitalised "Diet"
      mid-sentence. A brand is how a paper identifies its chow;
    * the species or purpose a feed is sold for ("rodent diet", "maintenance
      diet"), which is the same statement as a brand with the capitals
      removed.

    A signal: an ALL-CAPS abbreviation ("HFHS diet"), anything containing a
    digit or a percentage (a product code or a composition, "2920x diet",
    "60% fat diet"), and any other lower-case modifier. That last clause is
    what keeps the rule fail-closed: a composition word nobody listed is
    still a composition word.
    """
    if not modifier:
        return False
    token = modifier.strip(_DASHES + "-_")
    if not token:
        return False
    lowered = token.lower()
    if lowered in _DIET_MODIFIER_DETERMINER:
        return False
    if lowered in _DIET_MODIFIER_PHYSICAL_FORM:
        return False
    if lowered in _DIET_MODIFIER_CATEGORY:
        return False
    letters = [ch for ch in token if ch.isalpha()]
    if letters and all(ch.isupper() for ch in letters) and len(letters) > 1:
        return True  # HFHS, HF, NC: an abbreviation, not a brand
    if any(ch.isdigit() for ch in token) or "%" in token:
        return True  # 2920x, 60%: a product code or a composition
    if token[:1].isupper():
        return False  # Rodent, Special, Lab: a product name
    if diet_token[:1].isupper() and not at_sentence_start:
        return False  # "... 20 (Diet) ...": part of a product name
    return True


def _diet_occurrences(sentence: str) -> List[_DietHit]:
    """Every diet-shaped occurrence in one sentence, classified or not.

    Classification runs first and its spans are recorded, so a diet-shaped
    span that overlaps a classified one is already accounted for: "high fat
    diet" classifies as ``high_fat`` and the bare shape "fat diet" inside it
    is not then reported as an unknown term.
    """
    stimulus = [
        (match.start(), match.end())
        for pattern in _DIET_AS_STIMULUS
        for match in pattern.finditer(sentence)
    ]
    out: List[_DietHit] = []
    taken: List[Tuple[int, int]] = []

    classified: List[Tuple[int, int, Diet, str]] = []
    for value, pattern in _DIET_CLASSIFIERS:
        for match in pattern.finditer(sentence):
            classified.append((match.start(), match.end(), value, match.group(0)))
    if _DIET_CONTEXT.search(sentence):
        for match in _DIET_ND.finditer(sentence):
            classified.append(
                (match.start(), match.end(), Diet.chow, match.group(0))
            )
    for start, end, value, text in classified:
        if _overlaps(start, end, stimulus):
            continue
        out.append(_DietHit(start, end, value, text))
        taken.append((start, end))

    # The word "diet" first, so that its span (modifier included) is taken
    # before the abbreviation scanner can report the same thing twice.
    for match in _DIET_WORD.finditer(sentence):
        start, end = match.start(), match.end()
        if _overlaps(start, end, stimulus) or _overlaps(start, end, taken):
            continue
        taken.append((start, end))
        if _diet_modifier_is_a_signal(
            match.group("mod"), match.group("diet"), start == 0
        ):
            out.append(_DietHit(start, end, None, match.group(0)))

    for pattern in _DIET_SHAPES:
        for match in pattern.finditer(sentence):
            start, end = match.start(), match.end()
            if _overlaps(start, end, stimulus) or _overlaps(start, end, taken):
                continue
            taken.append((start, end))
            out.append(_DietHit(start, end, None, match.group(0)))

    out.sort(key=lambda hit: hit.start)
    return out


def _diet_disjunction(sentence: str, occurrences) -> bool:
    """True when this sentence states its diet as a two-arm disjunction.

    "fed either a ND ... or a HFD" and "a standard diet or a chow diet" are a
    paper telling you it has two arms, **whatever the arms are called**, which
    is why the construction is detected and not only the terms: the arm whose
    name the lexicon does not know is exactly the arm that goes invisible.

    Position matters, and it is what separates a diet disjunction from an
    unrelated one. "Mice fed standard chow received either saline or CNO" has
    a diet term and an either/or, but the diet term sits *before* the
    "either": the disjunction is about the drug. So "either" counts only when
    a diet term follows it, and a bare "or" counts only when diet terms sit on
    both sides of it.
    """
    if not occurrences:
        return False
    starts = [hit.start for hit in occurrences]
    ends = [hit.end for hit in occurrences]
    for match in _EITHER.finditer(sentence):
        if any(start >= match.end() for start in starts) and _OR.search(
            sentence, match.end()
        ):
            return True
    for match in _OR.finditer(sentence):
        before = any(0 <= match.start() - end <= 60 for end in ends)
        after = any(0 <= start - match.end() <= 60 for start in starts)
        if before and after:
            return True
    return False


# --------------------------------------------------------------------------- #
# age
# --------------------------------------------------------------------------- #

_AGE_UNITS = {
    "day": AgeUnit.days,
    "days": AgeUnit.days,
    "week": AgeUnit.weeks,
    "weeks": AgeUnit.weeks,
    "month": AgeUnit.months,
    "months": AgeUnit.months,
    "year": AgeUnit.years,
    "years": AgeUnit.years,
}

_UNIT = r"days?|weeks?|months?|years?"
#: A dash or "to" spans a range ("8-10 weeks"); "and" enumerates two cohorts
#: ("15 and 78 weeks old"), which is two statements, not one span.
_RANGE_SEP = "(?:" + _DASH + r"|\bto\b)"
_LIST_SEP = r"(?:\band\b|&)"
_SEP = "(?P<sep>" + _RANGE_SEP + "|" + _LIST_SEP + ")"
_SECOND = "(?:\\s*" + _SEP + "\\s*(?P<b>" + _NUM + "))?"

#: The three shapes an age statement takes in this corpus. ``_NUM`` accepts a
#: spelled-out numeral as well as digits, so "Twelve-week-old mice" reads the
#: same as "12-week-old mice".
_AGE_PATTERNS = (
    # "8-10-week-old", "6-7 week old", "4 months old", "3-weeks-old",
    # "Twelve-week-old"
    re.compile(
        r"\b(?P<a>" + _NUM + ")" + _SECOND + "[\\s" + _DASHES + "]*"
        r"(?P<unit>" + _UNIT + ")[\\s" + _DASHES + r"]*old\b",
        re.IGNORECASE,
    ),
    # "aged 10-12 weeks", "aged 9 weeks", "aged 6 to 8 weeks"
    re.compile(
        r"\baged\s+(?P<a>" + _NUM + ")" + _SECOND
        + r"\s*(?P<unit>" + _UNIT + r")\b",
        re.IGNORECASE,
    ),
    # "12 to 16 weeks of age", "8-10 weeks of age"
    re.compile(
        r"\b(?P<a>" + _NUM + ")" + _SECOND
        + r"\s*(?P<unit>" + _UNIT + r")\s+of\s+age\b",
        re.IGNORECASE,
    ),
)


@dataclass(frozen=True)
class _AgeStatement:
    """One age the text states, as stated: a scalar or a two-ended range."""

    values: Tuple[float, ...]
    unit: AgeUnit
    is_range: bool
    raw: str

    @property
    def is_scalar(self) -> bool:
        return not self.is_range and len(self.values) == 1


def _parse_age_match(match: "re.Match[str]") -> Optional[List[_AgeStatement]]:
    unit = _AGE_UNITS.get(match.group("unit").lower())
    if unit is None:  # pragma: no cover - the pattern cannot produce one
        return None
    first = _as_number(match.group("a"))
    if first is None:  # pragma: no cover - the pattern cannot produce one
        return None
    raw = match.group(0)
    second_text = match.groupdict().get("b")
    if second_text is None:
        return [_AgeStatement((first,), unit, False, raw)]
    second = _as_number(second_text)
    if second is None:  # pragma: no cover
        return None
    separator = (match.group("sep") or "").strip().lower()
    if separator in ("and", "&"):
        # Two cohorts, not one span: "15 and 78 weeks old" states two ages.
        return [
            _AgeStatement((first,), unit, False, raw),
            _AgeStatement((second,), unit, False, raw),
        ]
    return [_AgeStatement((first, second), unit, True, raw)]


def _age_statements(sentence: str) -> List[_AgeStatement]:
    """Age statements in one sentence, non-overlapping, in order of position."""
    accepted: List[Tuple[int, int, List[_AgeStatement]]] = []
    candidates: List[Tuple[int, int, List[_AgeStatement]]] = []
    for pattern in _AGE_PATTERNS:
        for match in pattern.finditer(sentence):
            parsed = _parse_age_match(match)
            if parsed:
                candidates.append((match.start(), match.end(), parsed))
    # Longest span first, so "aged 8-10 weeks of age" is read once.
    candidates.sort(key=lambda item: (item[0] - item[1], item[0]))
    for start, end, parsed in candidates:
        if any(not (end <= s or start >= e) for s, e, _ in accepted):
            continue
        accepted.append((start, end, parsed))
    accepted.sort(key=lambda item: item[0])
    out: List[_AgeStatement] = []
    for _, _, parsed in accepted:
        out.extend(parsed)
    return out


def _mouse_band(value: float, unit: AgeUnit) -> Optional[AgeGroup]:
    """The mouse band for one scalar age, or ``None`` when there is not one.

    The age is put into weeks by :func:`syn4rch.schema.rules.age_in_weeks` and
    banded by :func:`syn4rch.schema.rules.assign_age_group`; neither the
    conversion factors nor the band boundaries are restated here.
    ``not_reported`` is ``assign_age_group``'s refusal, so it comes back as
    ``None``: never as a value to store.
    """
    if age_in_weeks(value, unit.value) is None:
        return None
    band = assign_age_group(value, unit.value, Species.mus_musculus.value)
    if band is AgeGroup.not_reported:
        return None
    return band


def _band_of(statement: _AgeStatement) -> Optional[AgeGroup]:
    """The band this statement supports, or ``None`` when it supports none.

    A range is assigned a band only when **both** ends fall in it. ``6-12
    weeks`` straddles ``juvenile`` and ``young_adult``, so it supports neither,
    and picking either end (or the midpoint) would be inventing a cohort the
    paper did not describe.
    """
    bands = [_mouse_band(value, statement.unit) for value in statement.values]
    if not bands or any(band is None for band in bands):
        return None
    if len(set(bands)) != 1:
        return None
    return bands[0]


# --------------------------------------------------------------------------- #
# challenge
# --------------------------------------------------------------------------- #

#: The explicit phrase table. Not loose token matching: ``fed`` is a word that
#: appears inside "chow-fed" and "HFD-fed", where it names a diet and not a
#: measurement state, and a token matcher cannot tell the difference.
_STATE_TABLE = (
    (
        Challenge.fasted,
        (
            "fasted", "fasting", "after a fast", "following an overnight fast",
            "overnight fast", "fasting-induced", "fast-induced",
            "fasted-induced",
        ),
    ),
    (
        Challenge.ad_libitum,
        (
            "fed state", "under fed conditions", "in the fed state",
            "ad libitum fed", "ad lib fed", "fed blood glucose", "fed",
        ),
    ),
    (Challenge.refed, ("refed", "re-fed", "refeeding", "re-feeding")),
    (
        Challenge.cold_exposure,
        ("cold exposure", "cold-exposed", "cold challenge"),
    ),
    (Challenge.thermoneutral, ("thermoneutral", "thermoneutrality")),
    (Challenge.restraint_stress, ("restraint stress",)),
    (Challenge.exercise, ("treadmill", "exercise", "exercised")),
)

#: ``(regex, value, phrase)`` sorted longest phrase first, so that a long
#: phrase consumes its span before a short one inside it is considered:
#: "ad libitum fed" is read once as ``ad_libitum``, not also as bare ``fed``,
#: and "re-fed" is ``refed`` rather than a hyphen-prefixed ``fed``.
_STATE_MATCHERS = tuple(
    sorted(
        (
            (re.compile(_phrase(phrase), re.IGNORECASE), value, phrase)
            for value, phrases in _STATE_TABLE
            for phrase in phrases
        ),
        key=lambda item: -len(item[2]),
    )
)

#: Words that make a following ``fed`` a diet rather than a state. "chow-fed",
#: "HFD-fed", "high-fat-fed" and "Western diet-fed" describe what the animals
#: ate for weeks, which is ``diet``; they say nothing about whether the
#: measurement was taken in the fed state.
_FED_IS_A_DIET = re.compile(
    r"(?:chow|HFD|NCD|SCD|LFD|diets?|fat|western|standard|normal|regular)"
    "[\\s" + _DASHES + "]*$",
    re.IGNORECASE,
)

#: The same trap from the other side. "mice were fed a high-fat diet" and
#: "animals were fed standard chow" put the diet *after* the word, so a
#: backward-looking guard alone would read them as a fed state. Applied to the
#: text immediately following a bare ``fed``.
_FED_INTRODUCES_A_DIET = re.compile(
    r"^\s*(?:with\s+|on\s+)?(?:an?|the)?\s*(?:ad libitum\s+)?"
    r"(?:\d+\s*%\s*(?:kcal\s*)?(?:from\s+)?fat\s+)?"
    r"(?:high[\s" + _DASHES + r"]*fat|western|standard|normal|regular"
    r"|control|chow|HFD|NCD|SCD|LFD)",
    re.IGNORECASE,
)

#: Husbandry wordings for unrestricted food access. Used **only** by the
#: Methods route, and only in order to recognise a housing statement so that
#: it can be refused: "ad libitum access to water and chow" names no
#: measurement, so it must not reach ``challenge``. It is deliberately not in
#: :data:`_STATE_TABLE`, which is the table for a record's own wording.
_AD_LIBITUM_HUSBANDRY = _any(
    _phrase("ad libitum"),
    _phrase("ad lib"),
    _phrase("freely available"),
    _phrase("free access"),
)

#: Assay names whose challenge is a function of the outcome itself. See the
#: module docstring.
_ASSAY_NAME = _any(
    r"\b(?:ip|i\.p\.|o|op)?GTT\b",
    r"\b(?:ip|i\.p\.)?ITT\b",
    r"\bglucose tolerance\b",
    r"\binsulin tolerance\b",
    r"\btolerance test\b",
)

#: Restricted-feeding protocols. Pair-feeding gives an animal a ration matched
#: to another group's intake: it is the deliberate *negation* of ad libitum
#: access, and it reached ``ad_libitum`` through the bare ``fed`` inside the
#: compound "pair-fed". :class:`~syn4rch.schema.enums.Challenge` has no member
#: for restricted feeding, so these are **refused**. Mapping them to ``other``
#: would pool a restricted-feeding experiment with every unrelated state and
#: make the pooling invisible.
#:
#: Matched before the state table, and the span is removed from consideration,
#: so the inner ``fed`` cannot match either.
_RESTRICTED_FEEDING = _any(
    _phrase("pair fed"),
    _phrase("pair feeding"),
    _phrase("food restricted"),
    _phrase("food restriction"),
    _phrase("calorie restricted"),
    _phrase("calorie restriction"),
    _phrase("caloric restriction"),
    _phrase("restricted fed"),
    _phrase("restricted feeding"),
    _phrase("time restricted"),
)

#: Qualifiers that invert the state term they govern. "temperatures below
#: thermoneutrality" is a *cold* challenge, and the term that matched was
#: ``thermoneutral``. The qualifier must govern the term directly -- no
#: intervening words -- so "did not suppress refeeding" is untouched while
#: "not fasted" and "below thermoneutrality" are refused.
#:
#: This is deliberately a refusal and not a conversion. "Below
#: thermoneutrality" is a cold challenge to a physiologist and an *inference*
#: to this codebase, and an inference is ``null`` (SCHEMA.md 2 rule 4). A
#: human curates it.
_INVERTING_QUALIFIER = re.compile(
    r"\b(?:below|above|outside(?:\s+of)?|near|sub|other\s+than|rather\s+than"
    r"|instead\s+of|not)"
    r"[\s" + _DASHES + r"]*$",
    re.IGNORECASE,
)

#: Husbandry vocabulary. An ``ad_libitum`` phrase in a sentence carrying one of
#: these describes how the animals were kept, not how a measurement was made.
_HOUSING_CUE = _any(
    r"\baccess\b", r"\bhoused?\b", r"\bhousing\b", r"\bcages?\b",
    r"\bvivarium\b", r"\bbedding\b", r"\bmaintained\b", r"\bwater\b",
    r"\blight" + _DASH + r"dark\b", r"\blight/dark\b", r"\bambient\b",
    r"\b\d+\s*h\s*(?:light|dark)\b", r"\btemperature of\b",
)


@dataclass(frozen=True)
class _StateHit:
    start: int
    end: int
    value: Challenge
    phrase: str


@dataclass(frozen=True)
class _StateScan:
    """What a state scan found, including what it refused to read.

    ``excluded`` is the point of this object. Every trap in this module is a
    place where a state-like token was *seen and not used*, and a trap that
    reports nothing is a fail-open rule: an unlisted compound ("pair-fed") or
    an unhandled qualifier ("below thermoneutrality") would otherwise produce a
    confident wrong value, and an excluded token sitting next to a matched one
    would let the matched one look unopposed.
    """

    hits: Tuple[_StateHit, ...] = ()
    #: ``refusal name -> the token(s) that earned it``, in order of discovery.
    excluded: Tuple[Tuple[str, str], ...] = ()

    @property
    def values(self):
        return {hit.value for hit in self.hits}


def _scan_states(text: str) -> _StateScan:
    """Scan ``text`` for state phrases, recording matches *and* exclusions."""
    if not text:
        return _StateScan()
    hits: List[_StateHit] = []
    excluded: List[Tuple[str, str]] = []
    taken: List[Tuple[int, int]] = []

    # Restricted feeding first: the span is consumed so that the ``fed``
    # inside "pair-fed" cannot be read as a fed state.
    for match in _RESTRICTED_FEEDING.finditer(text):
        taken.append((match.start(), match.end()))
        excluded.append(
            ("refused_restricted_feeding_has_no_enum_member", match.group(0))
        )

    for pattern, value, phrase in _STATE_MATCHERS:
        for match in pattern.finditer(text):
            start, end = match.start(), match.end()
            if _overlaps(start, end, taken):
                continue
            if phrase == "fed" and (
                _FED_IS_A_DIET.search(text[:start])
                or _FED_INTRODUCES_A_DIET.match(text[end : end + 60])
            ):
                taken.append((start, end))
                excluded.append(
                    ("refused_fed_inside_a_diet_phrase", text[max(0, start - 12) : end])
                )
                continue
            if _INVERTING_QUALIFIER.search(text[:start]):
                taken.append((start, end))
                excluded.append(
                    (
                        "refused_qualifier_inverts_the_state",
                        text[max(0, start - 16) : end],
                    )
                )
                continue
            taken.append((start, end))
            hits.append(_StateHit(start, end, value, match.group(0)))

    hits.sort(key=lambda hit: hit.start)
    return _StateScan(hits=tuple(hits), excluded=tuple(excluded))


def _state_hits(text: str) -> List[_StateHit]:
    """Non-overlapping state-phrase occurrences in ``text``, in order."""
    return list(_scan_states(text).hits)


def _methods_state_hits(text: str) -> List[_StateHit]:
    """State-phrase occurrences in a Methods sentence.

    :func:`_state_hits` plus the husbandry wordings of
    :data:`_AD_LIBITUM_HUSBANDRY`, which exist here so that a housing
    statement can be *seen* and then refused rather than silently missed.
    """
    hits = _state_hits(text)
    taken = [(hit.start, hit.end) for hit in hits]
    for match in _AD_LIBITUM_HUSBANDRY.finditer(text):
        start, end = match.start(), match.end()
        if _overlaps(start, end, taken):
            continue
        taken.append((start, end))
        hits.append(_StateHit(start, end, Challenge.ad_libitum, match.group(0)))
    hits.sort(key=lambda hit: hit.start)
    return hits


def _names_an_assay(text: str) -> bool:
    return bool(text) and bool(_ASSAY_NAME.search(text))


# --------------------------------------------------------------------------- #
# the verbatim gate
# --------------------------------------------------------------------------- #


def paper_text(paper: ParsedPaper) -> str:
    """Everything this paper's retrieved text contains, as one string.

    Title, section headings, every paragraph, and every figure and table
    caption. This is the haystack the gate checks a quote against; it is
    assembled rather than cached so that it cannot drift from the paper.
    """
    parts: List[str] = []
    if paper.title:
        parts.append(paper.title)
    for section in paper.sections:
        if section.title:
            parts.append(section.title)
        parts.extend(section.paragraphs)
    for caption in list(paper.figures) + list(paper.tables):
        if caption.title:
            parts.append(caption.title)
        if caption.text:
            parts.append(caption.text)
    return "\n".join(parts)


def verify_against_paper(extracted: "Extracted", paper: ParsedPaper) -> bool:
    """True when ``extracted.quote`` occurs verbatim in this paper's text.

    The check is :func:`syn4rch.curate.verify.quote_occurs_in` under the
    project's single normalisation: no case folding, no fuzzy matching, no
    repair. An empty quote is not verbatim.
    """
    return quote_occurs_in(extracted.quote, paper_text(paper))


# --------------------------------------------------------------------------- #
# reading one paper's Methods
# --------------------------------------------------------------------------- #


def _methods_units(paper: ParsedPaper) -> List[Tuple[str, str]]:
    """``(sentence, locator)`` for every sentence of the Methods section.

    The locator is the project's convention for paragraph *N* (1-based) of a
    section -- ``"Methods para 3"`` -- and matches
    :func:`syn4rch.curate.passages.parse_locator`. Paragraph numbering comes
    from ``Section.paragraphs``, which is already flattened across
    subsections, so it does not depend on how the publisher nested headings.
    """
    section = paper.section("Methods")
    if section is None:
        return []
    out: List[Tuple[str, str]] = []
    for index, paragraph in enumerate(section.paragraphs, start=1):
        locator = "Methods para %d" % index
        for sentence in sentences(paragraph):
            out.append((sentence, locator))
    return out


def _claim_units(paper: ParsedPaper) -> List[Tuple[str, str]]:
    """``(sentence, locator)`` for the text where a paper reports its own
    measurements: the Results section, and the figure and table captions.

    **Not the Discussion**, and that is the whole design (D-082). Every shape a
    second-cohort check must not catch lives there -- a citation of another
    lab's female study, a limitation ("the absence of female cohorts may limit
    the generalizability"), and future work. Measured over the 42 papers in
    this corpus holding a single-sex ``constant``: Results and captions give
    6 real defects and **0** false positives; adding the Discussion keeps the
    6, adds 3 false positives and costs 8 further entries. Whole-paper is
    identical to Results + Discussion + captions here, because the Abstract and
    Introduction contribute nothing for sex.

    So the exculpation is section scope, not a lexicon. The one false positive
    scope cannot reach -- a pregnant dam reported in the Results -- is already
    caught by :data:`_BREEDING_CUE`, which this module has had all along.
    """
    out: List[Tuple[str, str]] = []
    section = paper.section("Results")
    if section is not None:
        for index, paragraph in enumerate(section.paragraphs, start=1):
            locator = "Results para %d" % index
            for sentence in sentences(paragraph):
                out.append((sentence, locator))
    for caption in list(paper.figures or []) + list(paper.tables or []):
        locator = caption.label or ("%s %d" % (caption.kind, caption.ordinal))
        for sentence in sentences(caption.text or ""):
            out.append((sentence, locator))
    return out


def _second_sex_in_claims(claims, value):
    """The first claim sentence reporting a sex other than ``value``.

    ``(sentence, locator)`` or ``None``. Only ``male`` and ``female`` can be
    contradicted this way: ``both`` already names both, so there is no second
    value for a claim sentence to add.

    This is the one dimension where widening the haystack does not widen the
    *vocabulary*, which is why it is fixable at all. D-075 retired the Methods
    rule for ``background_genotype`` because genetic background is an open set
    and absence can never be established; sex has two tokens, so scanning more
    text still exhausts the vocabulary.
    """
    if value not in (Sex.male, Sex.female):
        return None
    # `==`, not `is`. `Extracted.value` is a plain `str`, and `Sex` is a
    # str-enum, so `value in (Sex.male, Sex.female)` is True while
    # `value is Sex.male` is False -- which selected the pattern for the SAME
    # sex and fired on every single-sex paper. Caught by the four papers whose
    # constants this check had to preserve; without them it would have shipped.
    other = _FEMALE if value == Sex.male else _MALE
    for sentence, locator in claims:
        if _BREEDING_CUE.search(sentence):
            continue
        if other.search(sentence):
            return sentence, locator
    return None


def _methods_species(units: List[Tuple[str, str]]) -> Tuple[Optional[str], bool, bool]:
    """``(species, has_mouse_cue, has_rat_cue)`` from the Methods animal cues.

    ``species`` is filled only when the cues are unambiguous. A Methods naming
    both mice and rats leaves it ``None``: the paper used both, and which one
    a given sentence is about is decided per sentence, not per paper.
    """
    mouse = any(_MOUSE_CUE.search(sentence) for sentence, _ in units)
    rat = any(_RAT_CUE.search(sentence) for sentence, _ in units)
    if mouse and not rat:
        return Species.mus_musculus.value, mouse, rat
    if rat and not mouse:
        return Species.rattus_norvegicus.value, mouse, rat
    return None, mouse, rat


#: Refusals that merely *block* ``constant``: the sentence is set aside and the
#: verdict comes from whatever else voted, except that ``constant`` is no
#: longer available. The rest force ``varies`` outright, because they are
#: statements that a second value exists.
_BLOCKING_ONLY = frozenset({"refused_figure_scoped_cohort_sentence"})


def _read_sex(units, verdicts, constants, notes) -> None:
    votes: List[Tuple[Sex, str, str]] = []
    refusals: List[str] = []
    for sentence, locator in units:
        vote, refusal = _sex_vote(sentence)
        if vote is not None:
            votes.append((vote, sentence, locator))
        elif refusal is not None and refusal not in refusals:
            refusals.append(refusal)

    for refusal in refusals:
        notes.append(
            {
                "refused_figure_scoped_cohort_sentence": (
                    "refused_figure_scoped_cohort_sentence: a sentence naming "
                    "a sex is scoped to a figure, panel, cohort or named set "
                    "of experiments, so it describes a subset. "
                    '"Experiments (Figure 3) were performed in female mice '
                    'only" states that a female-only subset exists and says '
                    "nothing about the sex used in the paper's other figures."
                ),
                "refused_disjunctive_sex_is_unresolved": (
                    "refused_disjunctive_sex_is_unresolved: the Methods say "
                    '"male or female", which states that either was used, not '
                    "that both were. For any individual record the sex is "
                    "unresolved, and unresolved is not a value."
                ),
                "refused_ambiguous_sex_sentence": (
                    "refused_ambiguous_sex_sentence: a sentence names both "
                    "sex words with neither an explicit conjunction nor an "
                    "explicit disjunction, so it does not establish that both "
                    "sexes were used."
                ),
            }[refusal]
        )

    forces_varies = [name for name in refusals if name not in _BLOCKING_ONLY]
    distinct = {vote for vote, _, _ in votes}

    if forces_varies or len(distinct) > 1:
        verdicts["sex"] = "varies" if (votes or forces_varies) else "absent"
        if len(distinct) > 1:
            notes.append(
                "sex varies: the Methods vote for %s, and the Methods do not "
                "say which cohort any one record came from"
                % ", ".join(sorted(v.value for v in distinct))
            )
        return
    if not votes:
        verdicts["sex"] = "absent"
        return
    if refusals:
        # A blocking refusal: something named a sex and was set aside, so the
        # absence of a second value was not established.
        verdicts["sex"] = "varies"
        notes.append(
            "refused_sex_not_constant_after_a_refusal: a vote for %s survived "
            "but a sex-naming sentence was set aside (%s), so the absence of a "
            "second value was not established."
            % (
                ", ".join(sorted(v.value for v in distinct)),
                ", ".join(refusals),
            )
        )
        return
    value = next(iter(distinct))
    _, sentence, locator = votes[0]
    verdicts["sex"] = "constant"
    constants["sex"] = Extracted(
        dimension="sex",
        value=value.value,
        source="methods",
        locator=locator,
        quote=sentence,
        rule="methods_both_sexes" if value is Sex.both else "methods_single_sex",
    )


def _read_diet(units, verdicts, constants, notes) -> None:
    classes: List[Diet] = []
    first: Dict[Diet, Tuple[str, str]] = {}
    unrecognised: List[str] = []
    disjunction: List[str] = []
    scoped = False

    for sentence, locator in units:
        occurrences = _diet_occurrences(sentence)
        if not occurrences:
            continue
        if _SCOPED_TO_A_SUBSET.search(sentence):
            scoped = True
            continue
        if _diet_disjunction(sentence, occurrences) and sentence not in disjunction:
            disjunction.append(sentence)
        for hit in occurrences:
            if hit.value is None:
                if hit.text not in unrecognised:
                    unrecognised.append(hit.text)
            elif hit.value not in classes:
                classes.append(hit.value)
                first[hit.value] = (sentence, locator)

    if unrecognised:
        notes.append(
            "refused_unrecognised_diet_term: the Methods contain diet-shaped "
            "term(s) the synonym table does not classify (%s). An "
            "unclassified term means the lexicon is incomplete, and an "
            "incomplete lexicon must produce a refusal: the invisible term "
            "may be the second arm of a two-arm study, which is exactly how "
            "a diet-comparison paper once came out as a constant."
            % ", ".join(repr(term) for term in unrecognised)
        )
    if disjunction:
        notes.append(
            "refused_disjunctive_diet_is_two_arms: the Methods state a diet "
            "as a disjunction (%r). A paper writing \"either X or Y\" is "
            "telling you it has two arms, whatever the arms are called."
            % disjunction[0]
        )
    if scoped:
        notes.append(
            "refused_figure_scoped_cohort_sentence: a sentence naming a diet "
            "is scoped to a figure, panel or named set of experiments, so it "
            "describes a subset and cannot vote at paper level."
        )

    if unrecognised or disjunction or len(classes) > 1:
        verdicts["diet"] = "varies"
        if len(classes) > 1:
            notes.append(
                "diet varies: the Methods name %s. This is a diet-comparison "
                "study and its Methods cannot say which group a record is "
                "from; diet must come from the record's own passage or stay "
                "null." % " and ".join(value.value for value in classes)
            )
        return
    if not classes:
        verdicts["diet"] = "absent"
        return
    if scoped:
        verdicts["diet"] = "varies"
        return
    value = classes[0]
    sentence, locator = first[value]
    verdicts["diet"] = "constant"
    constants["diet"] = Extracted(
        dimension="diet",
        value=value.value,
        source="methods",
        locator=locator,
        quote=sentence,
        rule="methods_single_diet",
    )


#: The refusals the age rule can record, with the note each one writes.
_AGE_REFUSAL_NOTES = {
    "refused_rat_age_band": (
        "refused_rat_age_band: an age statement appeared in a sentence naming "
        "rats; the mouse bands of SCHEMA.md 3.3 do not transfer, so the "
        "statement was skipped rather than converted."
    ),
    "refused_in_vitro_age_statement": (
        "refused_in_vitro_age_statement: an age statement appeared in a "
        "sentence about cultured cells, where the age is the culture's and "
        "not the animals'."
    ),
    "refused_protocol_start_age_is_not_measurement_age": (
        "refused_protocol_start_age_is_not_measurement_age: an age statement "
        "appeared in a sentence that also states when a protocol started or "
        "how long it ran, so the age is the age at the start and not the age "
        "at the readout. age_group is the band the animal was in when it was "
        "measured."
    ),
    "refused_straddling_age_range": (
        "refused_straddling_age_range: an age range's two ends fall in "
        "different bands, so it supports neither. Picking an end or a "
        "midpoint would invent a cohort the paper did not describe."
    ),
    "refused_tissue_preparation_age": (
        "refused_tissue_preparation_age: an age statement appeared in a "
        "sentence about perfusion, slice preparation or tissue collection, so "
        "it is the age of the animals used for that procedure and not "
        "necessarily the age of the cohort any record was measured in."
    ),
    "refused_figure_scoped_cohort_sentence": (
        "refused_figure_scoped_cohort_sentence: an age statement is scoped to "
        "a figure, panel, cohort or named set of experiments, so it describes "
        "a subset and cannot settle the paper-level band."
    ),
}


def _read_age(units, verdicts, constants, notes, species):
    """Returns ``(age_value, age_unit)`` alongside filling the verdict.

    Three preconditions, each of which produces ``absent`` or ``varies``
    rather than a band:

    1. the Methods must establish the species as mouse, because
       :func:`syn4rch.schema.rules.assign_age_group` is mouse-only by design
       and that refusal must not be bypassed by handing it a hard-coded
       ``"mus_musculus"`` for a paper whose own species is unresolved;
    2. every age statement in the Methods must be usable -- one refused
       statement means the paper states an age the rule cannot use, and the
       verdict is then not ``constant`` however unanimous the remainder looks;
    3. the usable statements must agree on one band.
    """
    found: List[Tuple[_AgeStatement, str, str]] = []
    refused: List[str] = []
    stated_any = False

    def refuse(name: str) -> None:
        if name not in refused:
            refused.append(name)

    for sentence, locator in units:
        if not _ANIMAL_CUE.search(sentence):
            continue
        statements = _age_statements(sentence)
        if not statements:
            continue
        stated_any = True
        if _RAT_CUE.search(sentence):
            refuse("refused_rat_age_band")
            continue
        if _IN_VITRO_CUE.search(sentence):
            refuse("refused_in_vitro_age_statement")
            continue
        if _PROTOCOL_AGE_CUE.search(sentence):
            refuse("refused_protocol_start_age_is_not_measurement_age")
            continue
        if _TISSUE_PREPARATION_CUE.search(sentence):
            refuse("refused_tissue_preparation_age")
            continue
        if _SCOPED_TO_A_SUBSET.search(sentence):
            refuse("refused_figure_scoped_cohort_sentence")
            continue
        for statement in statements:
            if _band_of(statement) is None:
                refuse("refused_straddling_age_range")
                continue
            found.append((statement, sentence, locator))

    for name in refused:
        notes.append(_AGE_REFUSAL_NOTES[name])

    # Precondition 1: the species must be established as mouse.
    if species != Species.mus_musculus.value:
        verdicts["age_group"] = "absent"
        if stated_any:
            notes.append(
                "refused_age_band_requires_mouse: an age was stated but the "
                "Methods do not establish the species as mus_musculus "
                "(species=%r). assign_age_group refuses every other species "
                "by design and the refusal is not bypassed here." % (species,)
            )
        return None, None

    if not found:
        verdicts["age_group"] = "absent"
        return None, None

    # Precondition 2: a refusal is not a licence to promote the remainder.
    if refused:
        verdicts["age_group"] = "varies"
        notes.append(
            "refused_age_not_constant_after_a_refusal: %d age statement(s) "
            "survived (%s) but at least one was refused (%s). The paper "
            "states an age this rule cannot use, which is a reason to decline "
            "the dimension, not to fall back on whatever survived."
            % (
                len(found),
                ", ".join(
                    sorted({_band_of(s).value for s, _, _ in found})
                ),
                ", ".join(refused),
            )
        )
        return None, None

    bands = {_band_of(statement) for statement, _, _ in found}
    if len(bands) != 1:
        verdicts["age_group"] = "varies"
        notes.append(
            "age_group varies: the Methods state ages in %s"
            % ", ".join(sorted(band.value for band in bands))
        )
        return None, None

    band = next(iter(bands))
    statement, sentence, locator = found[0]
    verdicts["age_group"] = "constant"
    constants["age_group"] = Extracted(
        dimension="age_group",
        value=band.value,
        source="methods",
        locator=locator,
        quote=sentence,
        rule="methods_single_age_band",
    )

    # age_value/age_unit only for a single scalar age every statement agrees
    # on. A range states a span, not a value, and storing either end as the
    # age would be a fabrication the band assignment does not need.
    scalars = {
        (statement.values[0], statement.unit)
        for statement, _, _ in found
        if statement.is_scalar
    }
    if len(scalars) == 1 and all(statement.is_scalar for statement, _, _ in found):
        value, unit = next(iter(scalars))
        return float(value), unit.value
    return None, None


def _read_challenge(units, verdicts, constants, notes) -> None:
    """Record what the Methods say about acute state, and refuse all of it.

    There is no extraction path here and no rule name, by design. The verdict
    is ``varies`` when any Methods sentence names a state outside husbandry --
    the paper does measure in named states, and the Methods cannot say which
    measurement was made in which -- and ``absent`` when nothing but husbandry
    or nothing at all is named. ``constant`` is unreachable; see the module
    docstring for the surgical-preparation sentence that made the earlier
    Methods rule untenable.
    """
    measured: Dict[Challenge, Tuple[str, str]] = {}
    anywhere: Dict[Challenge, Tuple[str, str]] = {}
    housing: Dict[Challenge, Tuple[str, str]] = {}

    for sentence, locator in units:
        hits = _methods_state_hits(sentence)
        if not hits:
            continue
        is_housing = bool(_HOUSING_CUE.search(sentence))
        for hit in hits:
            anywhere.setdefault(hit.value, (sentence, locator))
            target = housing if is_housing else measured
            target.setdefault(hit.value, (sentence, locator))

    if housing and not measured:
        notes.append(
            "refused_housing_statement_is_not_a_challenge: the Methods "
            "describe husbandry (%s). Husbandry is not the state a "
            "measurement was taken in: the same paper measures glucose "
            "during a tolerance test and food intake after a fast. Stamping "
            "this on every record would be fabrication at scale and would "
            "look like full coverage."
            % ", ".join(sorted(value.value for value in housing))
        )

    if not anywhere:
        verdicts["challenge"] = "absent"
        return

    notes.append(
        "refused_methods_cannot_state_a_measurement_state: the Methods name "
        "%s, but no Methods sentence says which acute state any given "
        "measurement was made in -- a husbandry sentence does not, and "
        "neither does a surgical-preparation or protocol sentence. challenge "
        "is read only from a record's own outcome_raw and support_quote."
        % ", ".join(sorted(value.value for value in anywhere))
    )
    verdicts["challenge"] = "varies" if measured else "absent"


#: Strain and background vocabulary, used only to *describe* a refusal.
#:
#: Deliberately not a classifier of genetic background: see
#: :func:`_read_genotype` for why no detector here can change an outcome.
_BACKGROUND_TOKEN = re.compile(
    r"\b(?:"
    r"C57BL/6[A-Z]?|BALB/?c|129[A-Za-z/]*|FVB(?:/[A-Za-z]+)?|DBA/\d"
    r"|CD\s?-?\s?1\b|ICR\b|Swiss\s+Webster|NOD\b|NMRI\b"
    r")",
    re.IGNORECASE,
)


def _read_genotype(units, verdicts, notes) -> None:
    """Record what the Methods say about genetic background, and refuse it all.

    There is no extraction path here and no rule name, as for ``challenge``.
    ``constant`` is unreachable, and for a stronger reason than challenge's.
    The module docstring's standing rule is that ``constant`` is "a claim that
    no second value exists in the paper", returned only when the absence of a
    second value was *established*. Sex has two tokens and diet has a handful,
    so their absence can be established by exhausting the vocabulary. Genetic
    background has no closed vocabulary -- strain names, alleles, floxed lines
    and Cre drivers are open sets -- so a second value can never be shown
    absent, and promoting the only one detected is exactly the failure that
    cost this module three defects (see the module docstring).

    **Both verdicts this function can return refuse.** ``varies`` and
    ``absent`` are equally inadmissible under 3b, so how well the tokens below
    match changes the *explanation* a curator reads and never the decision.
    That is what makes an open-vocabulary dimension safe to gate: the gate does
    not depend on recognising the vocabulary.

    The verdicts are therefore each a true statement rather than a guess:
    ``varies`` where two or more distinct backgrounds are named, which the
    paper itself has then established; ``absent`` otherwise, meaning this
    module read nothing it can act on -- not that the paper is uniform.
    """
    seen: Dict[str, str] = {}
    for sentence, locator in units:
        for match in _BACKGROUND_TOKEN.finditer(sentence):
            token = match.group(0)
            seen.setdefault(token.upper().replace(" ", ""), token)

    if len(seen) >= 2:
        notes.append(
            "refused_methods_name_more_than_one_background: the Methods name "
            "%s. The paper itself has established that a second background "
            "exists, so no Methods sentence says which one any given record's "
            "animals carried. background_genotype is read only from a "
            "record's own passage, under the arm rule (D-060)."
            % ", ".join(sorted(seen.values()))
        )
        verdicts["background_genotype"] = "varies"
        return

    notes.append(
        "refused_methods_cannot_establish_one_background: genetic background "
        "has no closed vocabulary, so this module cannot show that a second "
        "background is absent, and a single detected background is not a "
        "paper-wide one%s. 'constant' is unreachable for this dimension; the "
        "value belongs to the record's own passage under the arm rule (D-060)."
        % ("" if not seen else " (detected: %s)" % ", ".join(sorted(seen.values())))
    )
    verdicts["background_genotype"] = "absent"


def read_paper_conditions(
    paper: ParsedPaper, study_id: str, pmcid: str
) -> PaperConditions:
    """Read the condition dimension out of one paper's Methods.

    Returns a verdict for each of :data:`DIMENSIONS` and an extraction only
    for the dimensions whose verdict is ``"constant"``. Every extraction is
    re-checked with :func:`verify_against_paper`; one whose quote does not
    occur in the paper is dropped, its verdict demoted to ``"absent"``, and
    the drop recorded in ``notes``. A ``"constant"`` verdict therefore always
    has an entry in ``constants``, and no entry in ``constants`` carries a
    quote this paper does not contain.

    ``verdicts["challenge"]`` and ``verdicts["background_genotype"]`` are
    never ``"constant"``: there is no Methods rule for either. Both are read
    from a record's own passage instead, and both are refused here so that
    :func:`inadmissible_methods_conditions` can say so.

    **A ``constant`` verdict is a claim that no second value exists in the
    paper's retrieved text**, so it is returned only when the absence of a
    second value was established -- never when a second value merely was not detected. Every
    refusal, every set-aside sentence and every unrecognised token is a reason
    to decline the dimension, not a reason to trust what is left. Three
    defects in this module's first version were all this one mistake: a rat
    age statement was skipped and the remaining statement looked unanimous; a
    figure-scoped sex sentence was excluded and an unscoped one looked
    unopposed; an unrecognised diet abbreviation made one arm of a two-arm
    study invisible and the other arm was promoted. All three failed *open*.
    """
    verdicts: Dict[str, str] = {dimension: "absent" for dimension in DIMENSIONS}
    constants: Dict[str, Extracted] = {}
    notes: List[str] = []

    units = _methods_units(paper)
    if not units:
        notes.append(
            "no Methods section with paragraphs was found, so no dimension "
            "could be read from the Methods."
        )
        return PaperConditions(
            study_id=study_id,
            pmcid=pmcid,
            verdicts=verdicts,
            constants=constants,
            notes=notes,
        )

    species, has_mouse, has_rat = _methods_species(units)
    if has_mouse and has_rat:
        notes.append(
            "the Methods name both mice and rats, so species was left null "
            "and each age statement was judged by its own sentence."
        )

    # Each reader is given its own notes list so that whatever it records is
    # attributed to the dimension it was called for. Nothing is parsed.
    refusals: Dict[str, List[str]] = {}

    def _run(dimension, reader, *args):
        own: List[str] = []
        result = reader(units, verdicts, constants, own, *args) if reader is not _read_genotype \
            else reader(units, verdicts, own)
        recorded = [n.split(":", 1)[0] for n in own if n.split(":", 1)[0] in REFUSALS]
        if recorded:
            refusals[dimension] = recorded
        notes.extend(own)
        return result

    _run("sex", _read_sex)
    _run("diet", _read_diet)
    age_value, age_unit = _run("age_group", _read_age, species)
    _run("challenge", _read_challenge)
    _run("background_genotype", _read_genotype)

    # D-082. The Methods' scope is not the paper's. `constant` claims that no
    # second value exists; a Methods section can name one cohort cleanly while
    # the Results report the other, and six papers in this corpus do exactly
    # that. PMC3940852's stored evidence is "C57Bl/6, db/+ and db/db male mice
    # were purchased from CLEA Japan" while its Results read "Male ... mice
    # showed a lean phenotype, whereas female ... mice tended toward increased
    # body weight". No quote screen can find that, because the Methods sentence
    # is clean.
    if verdicts.get("sex") == "constant":
        found = _second_sex_in_claims(
            _claim_units(paper), getattr(constants.get("sex"), "value", None)
        )
        if found is not None:
            sentence, locator = found
            constants.pop("sex", None)
            verdicts["sex"] = "varies"
            notes.append(
                "refused_second_sex_reported_outside_the_methods: the Methods "
                "name one sex, but this paper reports a measurement in the "
                "other at %s -- %r. A `constant` verdict claims no second "
                "value exists in the retrieved text, and the Methods cannot "
                "establish that about the Results."
                % (locator, sentence[:160])
            )
            refusals.setdefault("sex", []).append(
                "refused_second_sex_reported_outside_the_methods"
            )

    # The verbatim gate. Mechanical, and applied to every extraction however
    # it was produced: a buggy rule above cannot emit a sentence the paper
    # does not contain.
    for dimension in DIMENSIONS:
        extracted = constants.get(dimension)
        if extracted is None:
            continue
        if verify_against_paper(extracted, paper):
            continue
        del constants[dimension]
        verdicts[dimension] = "absent"
        notes.append(
            "dropped_quote_not_verbatim: the %s extraction (rule %s, locator "
            "%s) quoted a sentence that does not occur in this paper's "
            "retrieved text, so it was dropped."
            % (dimension, extracted.rule, extracted.locator)
        )
        refusals.setdefault(dimension, []).append("dropped_quote_not_verbatim")
        if dimension == "age_group":
            age_value, age_unit = None, None

    return PaperConditions(
        study_id=study_id,
        pmcid=pmcid,
        verdicts=verdicts,
        constants=constants,
        refusals=refusals,
        species=species,
        age_value=age_value,
        age_unit=age_unit,
        notes=notes,
    )


# --------------------------------------------------------------------------- #
# reading one record's own passage
# --------------------------------------------------------------------------- #


#: Exclusions that stop the record-level route outright, with the clause each
#: one writes into its note. Both caught a token that names a state this
#: schema cannot hold, so there is nothing to fall back to.
_STOPPING_STATE_REFUSALS = {
    "refused_restricted_feeding_has_no_enum_member": (
        "describes a restricted-feeding protocol -- the negation of ad "
        "libitum access -- for which Challenge has no member; mapping it to "
        "'other' would pool it with every unrelated state"
    ),
    "refused_qualifier_inverts_the_state": (
        "carries a relational or negating qualifier that inverts the state "
        "term it governs; converting it would be an inference, and an "
        "inference is null (SCHEMA.md 2 rule 4)"
    ),
}


def _text_of(value) -> str:
    return "" if value is None else str(value)


def _sentence_carrying(text: str, position: int) -> str:
    """The sentence of ``text`` that contains the character at ``position``.

    Falls back to the whole (stripped) text when the split does not cover the
    position, which cannot drop evidence: the whole text is still verbatim.
    """
    cursor = 0
    for sentence in sentences(text):
        found = text.find(sentence, cursor)
        if found < 0:  # pragma: no cover - sentences are substrings
            continue
        if found <= position < found + len(sentence):
            return sentence
        cursor = found + len(sentence)
    return text.strip()


#: Refusals that do not unsettle a `constant` on the same dimension.
#:
#: These two ARE the demotion -- they are recorded when a prior refusal already
#: forced the verdict away from `constant`, so a verdict that is still
#: `constant` cannot have been produced by them. Treating them as blocking
#: would refuse every paper whose reader did its job.
_REFUSALS_THAT_DO_NOT_UNSETTLE = frozenset(
    {
        "refused_sex_not_constant_after_a_refusal",
        "refused_age_not_constant_after_a_refusal",
    }
)


def inadmissible_methods_conditions(record, paper: "PaperConditions"):
    """Methods-sourced conditions this paper does not hold constant (3b).

    Returns ``[(dimension, value, verdict), ...]``, empty when every
    Methods-sourced entry is admissible. A record carrying one of these is
    inadmissible: the Methods named two cohorts on that dimension, so they say
    the values exist and nothing about which one produced this figure.

    **This closes a rule that was stated and unenforced.** ``A008``-``A011``
    check that a ``condition_evidence`` entry is auditable -- locator present,
    quote present, value matching the field. None of them asks whether the
    value was allowed to be read from the Methods at all. So a record could
    carry ``sex='male'`` sourced from the Methods of a paper that studied both
    sexes, with a perfect locator and a verbatim quote, and pass every check in
    the project.

    The hole was never exercised while conditions came from
    :mod:`syn4rch.curate.fill_conditions`, which asks
    :func:`read_paper_conditions` first and declines on anything but
    ``constant`` -- and the 508 records curated that way hold 320
    Methods-sourced entries with zero violations. It was exercised the moment a
    reader filled conditions by hand: the first 30 papers of the D-072 pass
    produced 272 Methods-sourced entries of which **115 were on dimensions this
    classifier refuses**. See D-074.

    Dimensions outside :data:`DIMENSIONS` are skipped rather than failed,
    because a dimension with no verdict has nothing to check.

    ``background_genotype`` was such a dimension until D-075, and the six
    Methods-sourced genotype entries in the corpus went unexamined in both
    directions as a result. It is now in :data:`DIMENSIONS` with no Methods
    rule, so its verdict is never ``constant`` and every Methods-sourced
    genotype entry is refused. That is a reading and not a gap: the refusal
    says the Methods cannot pin a record to one arm, which is the same thing
    ``challenge``'s refusal says, and the stored quotes bear it out -- each of
    the six names two cohorts.
    """
    out = []
    for entry in getattr(record, "condition_evidence", None) or []:
        source = entry.get("source") if isinstance(entry, dict) else getattr(entry, "source", None)
        dimension = entry.get("dimension") if isinstance(entry, dict) else getattr(entry, "dimension", None)
        value = entry.get("value") if isinstance(entry, dict) else getattr(entry, "value", None)
        source = getattr(source, "value", source)
        dimension = getattr(dimension, "value", dimension)
        if str(source) != "methods" or dimension not in DIMENSIONS:
            continue
        verdict = (paper.verdicts or {}).get(dimension)
        if verdict != "constant":
            out.append((dimension, getattr(value, "value", value), verdict))
            continue
        # D-079: `constant` is not enough. The verdict claims no second value
        # exists; a refusal on the same dimension is the reader saying it could
        # not establish that. D-044 asks for a value that was ESTABLISHED, and
        # a constant standing beside its own dimension's refusal was merely
        # undetected -- the ten defects of D-077/D-078 were all that state.
        # The demotion already exists for sex and age
        # (`refused_sex_not_constant_after_a_refusal`) and fires only when a
        # refusal was *detected*, so this is the backstop for the next
        # detector gap rather than a second copy of the same rule.
        recorded = (getattr(paper, "refusals", None) or {}).get(dimension) or []
        blocking = [name for name in recorded if name not in _REFUSALS_THAT_DO_NOT_UNSETTLE]
        if blocking:
            out.append(
                (dimension, getattr(value, "value", value),
                 "constant beside %s" % ", ".join(sorted(set(blocking))))
            )
    return out


def read_record_conditions(record, notes: Optional[List[str]] = None):
    """Condition values readable from one record's **own** fields.

    Reads ``outcome_raw``, ``support_quote`` and ``source_locator`` and
    nothing else. It is deliberately not given the paper: this is the route
    for dimensions that vary within a paper, where the Methods cannot say
    which cohort a record came from and only the record's own passage can.

    Only ``challenge`` is read here. The acute state is the dimension the
    curator transcribes into ``outcome_raw`` ("EE in the fasted state"), and
    it is the dimension a Results sentence actually pins down.

    Two sources, in order:

    * ``outcome_raw`` (rule ``outcome_names_the_state``). The quote is
      ``outcome_raw`` itself, because it is the paper's own wording of the
      measurement, transcribed verbatim by the curator.
    * ``support_quote`` (rule ``quote_names_the_state``), accepted only when
      exactly one state phrase occurs in the whole quote. A quote contrasting
      fed with fasted animals cannot say which cohort this record is. The
      returned quote is the sentence within ``support_quote`` carrying the
      phrase, which is a substring of it and so verifies verbatim wherever
      ``support_quote`` does.

    The locator is always the record's own ``source_locator``.

    Two traps stop a route outright, because the token they caught names a
    state the schema cannot hold:
    ``refused_restricted_feeding_has_no_enum_member`` (pair-feeding is the
    negation of ad libitum access and :class:`Challenge` has no member for it)
    and ``refused_qualifier_inverts_the_state`` ("below thermoneutrality" is a
    cold challenge to a physiologist and an inference here).

    One trap does **not** refuse. When a state-like token was excluded and a
    *different* token matched, the extraction is kept, the rule becomes
    ``..._with_a_competing_token``, and a note names both tokens. This is the
    single place in the module where a refusal does not decline the value:
    ``challenge`` is the scarcest dimension in this corpus, the surviving
    match is usually right, and the honest handling of "usually right" is to
    publish the uncertainty for a human to check rather than to drop the
    value or to hide the competition.

    Pass ``notes`` to collect all of this: this function returns a list of
    extractions, so a refusal is otherwise indistinguishable from a record
    that simply says nothing, and the refusals are the part worth reading.
    """
    out_notes = notes if notes is not None else []
    outcome = _text_of(getattr(record, "outcome_raw", None)).strip()
    quote = _text_of(getattr(record, "support_quote", None))
    locator = _text_of(getattr(record, "source_locator", None)).strip()
    record_id = getattr(record, "record_id", None)

    if _names_an_assay(outcome):
        out_notes.append(
            "refused_assay_implies_its_own_challenge: outcome %r names an "
            "assay, and a tolerance test's acute state is a function of the "
            "outcome itself. Storing it would discriminate nothing that "
            "outcome_id does not already discriminate, while inflating the "
            "apparent condition coverage. (record %s)" % (outcome, record_id)
        )
        return []

    for source_text, field_name, base_rule in (
        (outcome, "outcome_raw", "outcome_names_the_state"),
        (quote, "support_quote", "quote_names_the_state"),
    ):
        if not source_text.strip():
            continue
        scan = _scan_states(source_text)

        if field_name == "outcome_raw" and len(scan.values) > 1:
            out_notes.append(
                "refused_two_states_in_one_outcome: outcome %r names %s. A "
                "record whose own outcome names two conditions is a curation "
                "defect worth surfacing, not a coin flip. (record %s)"
                % (
                    source_text,
                    ", ".join(sorted(v.value for v in scan.values)),
                    record_id,
                )
            )
            return []
        if field_name == "support_quote" and len(scan.hits) > 1:
            out_notes.append(
                "refused_quote_names_more_than_one_state: the support_quote "
                "carries %d state phrases (%s), so it cannot say which "
                "cohort this record is. (record %s)"
                % (
                    len(scan.hits),
                    ", ".join(hit.phrase for hit in scan.hits),
                    record_id,
                )
            )
            return []

        if not scan.hits:
            stopped = [
                (name, token)
                for name, token in scan.excluded
                if name in _STOPPING_STATE_REFUSALS
            ]
            for name, token in stopped:
                out_notes.append(
                    "%s: %s names %r, which %s. No value is stored. "
                    "(record %s)"
                    % (
                        name,
                        field_name,
                        token,
                        _STOPPING_STATE_REFUSALS[name],
                        record_id,
                    )
                )
            if stopped:
                return []
            continue

        hit = scan.hits[0]
        if not locator:
            out_notes.append(
                "refused_record_has_no_locator: %s names a state but the "
                "record has no source_locator to attribute it to. (record %s)"
                % (field_name, record_id)
            )
            return []

        rule = base_rule
        if scan.excluded:
            rule = base_rule + "_with_a_competing_token"
            out_notes.append(
                "%s: %s matched %r (-> %s) while %s was excluded by %s. The "
                "value is kept and flagged rather than refused, so that this "
                "small set can be checked by hand. (record %s)"
                % (
                    rule,
                    field_name,
                    hit.phrase,
                    hit.value.value,
                    ", ".join(repr(token) for _, token in scan.excluded),
                    ", ".join(sorted({name for name, _ in scan.excluded})),
                    record_id,
                )
            )

        return [
            Extracted(
                dimension="challenge",
                value=hit.value.value,
                source="passage",
                locator=locator,
                quote=(
                    source_text
                    if field_name == "outcome_raw"
                    else _sentence_carrying(source_text, hit.start)
                ),
                rule=rule,
            )
        ]

    return []
