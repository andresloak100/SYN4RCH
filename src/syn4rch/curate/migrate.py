"""Schema migration of ``data/curated/records.jsonl``: v0.1.0 -> v0.2.0.

    python3 -m syn4rch.curate.migrate run
    python3 -m syn4rch.curate.migrate report

``run``
    Reads the curated file at v0.1.0, fills the v0.2.0 fields it is *allowed*
    to fill, recomputes every ``record_id``, rewrites the file atomically,
    writes the complete old-id -> new-id map to
    ``data/curated/record_id_map_0.1.0_to_0.2.0.json`` and the report to
    ``reports/schema_migration_0.2.0.md``. The pre-migration file is kept at
    ``data/curated/records.0.1.0.jsonl.bak``.

``report``
    Recomputes every number in that report from the committed dataset and
    re-renders it, so that no figure in it has to be taken on trust
    (``docs/CURATION_PROTOCOL.md`` 6).

What a migration is allowed to do
---------------------------------
A migration is code. It cannot read a paper, so it may only compute what the
curators' own fields already determine:

* ``perturbation_class``    -- a total function of ``method``
  (:func:`syn4rch.schema.rules.perturbation_class_from_method`).
* ``perturbation_polarity`` -- a function of ``(method, valence)``
  (:func:`syn4rch.schema.rules.perturbation_polarity_from`).
* ``schema_version``        -- set to :data:`TO_VERSION`, the fixed version
  this migration targets, which is not necessarily the version the schema
  package implements today.
* ``record_id``             -- recomputed, because ``RECORD_ID_FIELDS`` grew
  from seven components to twelve. **Every id in the corpus changes**, so the
  old-id -> new-id map is written out and used to rewrite every
  ``conflicts_with`` reference. ``conflicts_with`` is empty on all 362 records
  today and the rewrite is therefore a no-op in practice; it is implemented and
  tested because the next QC run will populate it, and a migration that
  silently orphaned those references would be a trap discovered much later.
* ``target_raw`` / ``target_id`` -- only under the strict rule below.

And what it must not do
-----------------------
``challenge``, ``challenge_detail``, ``condition_source`` and
``condition_locator`` stay ``null``. They require reading the paper, and this
code has read no paper. In particular it does **not** write ``not_reported``:
per SCHEMA.md 3.3b that means a curator checked the paper and it does not say,
and writing it here would record a curation act that never happened. The
unfilled fields are counted and published instead.

Two of those four no longer exist. Schema v0.3.0 retired the record-level
``condition_source``/``condition_locator`` pair -- null on all 362 records --
and replaced it with ``condition_evidence``, one entry per filled condition
dimension, each carrying its own source, locator and quote (D-037). The reason
is the one this migration ran into from the other side: a record's sex is
stated in the paper's Methods and its acute metabolic state, if stated at all,
in the passage the record was read from, so one scalar pair could never have
described both. This module still names the pair, because it targets v0.2.0 and
that is what v0.2.0 had; it never constructs the fields.

``sex``, ``diet`` and ``age_group`` are already curated where known and are not
touched.

Target resolution is resolution, not reading
--------------------------------------------
This is the one place where the migration looks at prose, so the limit is worth
stating plainly. The molecular target was never absent from the dataset; it was
curated inside the free text of ``targeting_strategy`` (e.g. ``"AgRPp53 KO (p53
deleted in AgRP neurons); wild type (WT) controls"``) and sometimes
``actuator``. This module lifts it into a field of its own, exactly as the
project already lifts ``population_raw`` to ``population_id`` and
``outcome_raw`` to ``outcome_id``. It does not verify the gene against the
paper. **A record's target is therefore only as good as its
``targeting_strategy``**: if a curator wrote the wrong gene there, this code
will faithfully resolve the wrong gene, and nothing here would notice.

The rule, which is deliberately unhelpful rather than clever:

1. Only ``molecular`` and ``pharmacological`` records are eligible (SCHEMA.md
   3.4b, enforced as ``S008``). A chemogenetic activation has no molecular
   target and keeps ``null``.
2. A synonym from ``data/vocab/targets.json`` must occur in
   ``targeting_strategy`` or ``actuator``, matched case-insensitively as a
   substring with no word boundary, because genotype strings concatenate the
   symbol with the driver and the floxed allele (``AgRPp53``, ``Sh2b1f/f``,
   ``Htr2cfl/Y``, ``TrspAgrpKO``). The only guard is that neither neighbouring
   character may be a digit, so ``Raf1`` cannot be found inside ``Raf12``.
3. Exactly one *entry* of the vocabulary must match. Zero matches: both fields
   stay ``null``, counted. Two or more entries: both stay ``null``, counted,
   and the record is listed by id in the report for a human to curate. An
   ambiguous match is never guessed. Several synonyms of the *same* entry are
   not an ambiguity -- a study that writes both ``SH2B1`` and ``Sh2b1`` names
   one gene.
4. ``target_raw`` is the matched string exactly as it appears in the text;
   ``target_id`` is the vocabulary id.
5. A record that already carries a target is left alone. Curation outranks
   resolution, and that is also what makes a second run a no-op.

Relationship to D-017
---------------------
``docs/DECISIONS.md`` D-017 says the gate in :mod:`syn4rch.curate.run` is the
only writer of ``data/curated/records.jsonl``. A migration is the one other
thing that may hold that pen, and it is a narrower act: it creates no record,
admits no record, drops no record, and writes no value that a curator did not
already write or that a published rule does not derive from one. The guard
rails are that it refuses to start unless every input line declares the
version it migrates *from*, and refuses to write if the result would carry a
hard validation finding the input did not already have. D-017 should be
amended to name this exception; that is the dataset owner's call and not this
module's.

Exit status:

    0  migrated (or reported) and nothing needs a human
    1  migrated, but something needs a human: an ambiguous target match, or a
       ``conflicts_with`` reference that is not in the id map
    2  refused: wrong input version, unreadable line, an id collision, or a
       new hard validation finding. Nothing is written.

A refusal writes nothing at all. Re-running is safe: the write is a temp file
plus :func:`os.replace`, the pre-migration file is preserved, and ``run``
declines to migrate a file that is not at the version it migrates from rather
than guessing what to do with it.

Python 3.10 compatible, standard library only (``docs/DECISIONS.md`` D-001).
No network.
"""

from __future__ import annotations

import argparse
import datetime
import difflib
import hashlib
import json
import os
import re
import shutil
import stat
import sys
import tempfile
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ..schema.enums import PerturbationClass
from ..schema.io import JsonlError, dumps_record, loads_record
from ..schema.record import (
    RECORD_ID_FIELDS,
    RECORD_ID_LENGTH,
    TierARecord,
    normalize_for_hash,
)
# ``SCHEMA_VERSION`` is deliberately *not* imported here: see TO_VERSION.
from ..schema.rules import perturbation_class_from_method, perturbation_polarity_from
from ..schema.validate import Severity, validate

__all__ = [
    "FROM_VERSION",
    "TO_VERSION",
    "RECORD_ID_FIELDS_ADDED_IN_0_2_0",
    "RECORD_ID_FIELDS_0_1_0",
    "CONDITION_FIELDS_THE_MIGRATION_MAY_NOT_FILL",
    "CONDITION_FIELDS_STILL_IN_THE_SCHEMA",
    "TARGET_BEARING_CLASSES",
    "TARGET_SEARCH_FIELDS",
    "MIN_SYNONYM_CHARS",
    "RESOLVED",
    "NO_MATCH",
    "AMBIGUOUS",
    "NOT_ELIGIBLE",
    "ALREADY_SET",
    "TargetMatch",
    "TargetResolution",
    "MigrationResult",
    "SynonymIndex",
    "load_targets",
    "build_synonym_index",
    "fold",
    "occurrences",
    "resolve_target",
    "compute_record_id_0_1_0",
    "migrate_record",
    "migrate_records",
    "remap_conflicts",
    "dataset_stats",
    "render_report_md",
    "build_parser",
    "cmd_run",
    "cmd_report",
    "main",
]


#: The version this migration reads. A file that does not declare it is
#: refused, not guessed at.
FROM_VERSION = "0.1.0"

#: The version this migration writes. A **literal**, deliberately, and not
#: ``SCHEMA_VERSION``.
#:
#: A migration names a fixed pair of versions by nature: this module fills the
#: fields v0.2.0 added, recomputes the twelve-component v0.2.0 ``record_id``,
#: and wrote ``reports/schema_migration_0.2.0.md``,
#: ``data/curated/record_id_map_0.1.0_to_0.2.0.json`` and
#: ``data/curated/records.0.1.0.jsonl.bak``. Reading the target from the schema
#: package made it move on its own at the v0.3.0 bump, which would have
#: stamped records ``0.3.0`` while filling nothing v0.3.0 added -- a file
#: claiming a schema whose fields it does not carry. Each later schema version
#: needs its own migration module; this one stays what it is.
TO_VERSION = "0.2.0"

#: The five ``record_id`` components v0.2.0 added (D-026). Named here so the
#: v0.1.0 id function can be reconstructed by subtraction rather than by a
#: second hand-maintained copy of the old tuple.
RECORD_ID_FIELDS_ADDED_IN_0_2_0 = (
    "target_raw",
    "sex",
    "diet",
    "age_group",
    "challenge",
)

#: The seven v0.1.0 ``record_id`` components, in v0.1.0 order.
RECORD_ID_FIELDS_0_1_0 = tuple(
    name for name in RECORD_ID_FIELDS if name not in RECORD_ID_FIELDS_ADDED_IN_0_2_0
)

#: Condition fields the migration is forbidden to fill: each one requires
#: reading the paper (SCHEMA.md 3.3b). ``not_reported`` is forbidden too -- it
#: is a claim that the paper was checked.
#:
#: Named as the **v0.2.0** schema named them, which is the schema this module
#: migrates to. ``condition_source`` and ``condition_locator`` existed then,
#: were left null on all 362 records, and were retired in v0.3.0 in favour of
#: per-dimension ``condition_evidence`` (D-037) for exactly the reason this
#: migration could not fill them: one record-level pair cannot carry the
#: provenance of a sex read from Methods and a challenge read from the
#: record's own passage. They stay in this tuple because the four-slot count
#: published in ``reports/schema_migration_0.2.0.md`` is a statement about the
#: v0.2.0 schema, and rewriting history to two slots would make that report
#: unreproducible. Nothing here constructs them: see
#: :data:`CONDITION_FIELDS_STILL_IN_THE_SCHEMA`.
CONDITION_FIELDS_THE_MIGRATION_MAY_NOT_FILL = (
    "challenge",
    "challenge_detail",
    "condition_source",
    "condition_locator",
)

#: The subset of the above that is still a field of :class:`TierARecord`.
#: Every code path that reads or writes a record goes through this; the full
#: tuple is only ever counted and printed. Derived rather than listed, so that
#: retiring or restoring a condition field cannot leave the two disagreeing.
CONDITION_FIELDS_STILL_IN_THE_SCHEMA = tuple(
    name
    for name in CONDITION_FIELDS_THE_MIGRATION_MAY_NOT_FILL
    if name in TierARecord.field_names()
)

#: Perturbation classes that may carry a molecular target (SCHEMA.md 3.4b).
#: Must agree with ``S008_TARGET_ON_NON_MOLECULAR`` in
#: :mod:`syn4rch.schema.validate`; ``tests/test_curate_migrate.py`` asserts the
#: agreement behaviourally rather than by importing a private name.
TARGET_BEARING_CLASSES = frozenset(
    {PerturbationClass.molecular, PerturbationClass.pharmacological}
)

#: Curated fields searched for a target, in priority order.
TARGET_SEARCH_FIELDS = ("targeting_strategy", "actuator")

#: A synonym shorter than this is refused rather than matched: matching has no
#: word boundary (see the module docstring), and a one- or two-character
#: substring would fire inside half the genotype strings in the corpus.
MIN_SYNONYM_CHARS = 3

_HASH_SEP = "\x1f"

_DASHES = "-‐‑‒–—―−"
_SPACES = "\t\n\r\x0b\x0c    　"

#: Character-for-character fold: every dash-like character becomes ``-`` and
#: every space-like character becomes a space. **Length-preserving on purpose**
#: -- a match's offsets in the folded string are offsets in the original, so
#: ``target_raw`` can be sliced out of the text exactly as the curator wrote
#: it. Case is handled by ``re.IGNORECASE`` for the same reason: lowercasing
#: can change a string's length for some Unicode characters.
_FOLD = str.maketrans(
    {ch: "-" for ch in _DASHES} | {ch: " " for ch in _SPACES}
)


# --------------------------------------------------------------------------- #
# paths (mirrors syn4rch.curate.run so this module stands on its own)
# --------------------------------------------------------------------------- #


def _repo_root() -> str:
    try:
        from ..harvest.cache import repo_root

        return repo_root()
    except Exception:  # pragma: no cover - defensive
        return os.path.abspath(os.getcwd())


def _default(*parts: str) -> str:
    return os.path.join(_repo_root(), *parts)


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# --------------------------------------------------------------------------- #
# the v0.1.0 record_id, reconstructed
# --------------------------------------------------------------------------- #


def compute_record_id_0_1_0(record) -> str:
    """The ``record_id`` v0.1.0 would have given ``record``.

    Same construction as :func:`syn4rch.schema.record.compute_record_id` over
    the seven :data:`RECORD_ID_FIELDS_0_1_0` components. It exists for one
    purpose: to *check* that the input file really was written by v0.1.0. If a
    stored id does not reproduce, the file was hand-edited after its ids were
    computed, which the migration reports rather than papering over.

    Accepts a :class:`~syn4rch.schema.record.TierARecord` or a plain dict.
    """
    getter = (
        record.get
        if isinstance(record, dict)
        else (lambda name: getattr(record, name, None))
    )
    payload = _HASH_SEP.join(
        normalize_for_hash(getter(name)) for name in RECORD_ID_FIELDS_0_1_0
    )
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    return digest[:RECORD_ID_LENGTH]


# --------------------------------------------------------------------------- #
# target resolution
# --------------------------------------------------------------------------- #


def fold(text: Optional[str]) -> str:
    """Fold dash-like and space-like characters. Length-preserving."""
    if text is None:
        return ""
    return str(text).translate(_FOLD)


@dataclass(frozen=True)
class SynonymIndex:
    """Every ``(target_id, synonym)`` pair of the targets vocabulary.

    Longest synonym first, so that where a text names both ``Trp53`` and
    ``p53`` the more specific string is the one ``target_raw`` records.
    """

    pairs: Tuple[Tuple[str, str], ...]

    @property
    def ids(self) -> Tuple[str, ...]:
        seen: List[str] = []
        for target_id, _ in self.pairs:
            if target_id not in seen:
                seen.append(target_id)
        return tuple(seen)

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self.pairs)


def load_targets(path: Optional[str] = None) -> Dict[str, dict]:
    """Load the targets vocabulary as ``{id: entry}``.

    With no ``path`` this is ``syn4rch.vocab.load("targets")``, i.e.
    ``data/vocab/targets.json`` with the project's invariants enforced. A
    ``path`` reads that file instead, which is what ``--targets`` is for: a
    curator trying a vocabulary out, and a test that must not depend on the
    committed one.
    """
    if path is None:
        from ..vocab import load  # local: the vocab package is not always needed

        return load("targets")
    with open(path, "r", encoding="utf-8") as handle:
        document = json.load(handle)
    entries = document.get("entries")
    if not isinstance(entries, dict):
        raise ValueError("%s: a top-level 'entries' object is required" % path)
    for key, entry in entries.items():
        if not isinstance(entry, dict):
            raise ValueError("%s: entry %r is not an object" % (path, key))
        if entry.get("id") != key:
            raise ValueError(
                "%s: entry key %r != id %r" % (path, key, entry.get("id"))
            )
    return entries


def build_synonym_index(entries: Optional[Dict[str, dict]] = None) -> SynonymIndex:
    """Build the synonym index from ``data/vocab/targets.json``.

    ``entries`` matches ``syn4rch.vocab.load("targets")`` and is injectable so
    a test need not depend on the real vocabulary. The entry id itself is
    included as a matchable string: a paper writing ``Rptor`` should resolve
    even though the corpus only ever writes ``Raptor``.

    Raises :class:`ValueError` on an entry with no synonyms, on a synonym
    shorter than :data:`MIN_SYNONYM_CHARS`, and on a synonym that is a
    substring of another entry's synonym. The last one is the dangerous case:
    substring matching could not tell such a pair apart, so it is refused at
    load time instead of producing a silently wrong target later.
    """
    if entries is None:
        entries = load_targets()

    pairs: List[Tuple[str, str]] = []
    for target_id in sorted(entries):
        entry = entries[target_id] or {}
        synonyms = list(entry.get("synonyms") or [])
        if not synonyms:
            raise ValueError(
                "targets/%s has no synonyms; a target that cannot be matched "
                "in text would silently never resolve" % target_id
            )
        seen = set()
        for synonym in [target_id] + [str(s) for s in synonyms]:
            text = str(synonym).strip()
            key = fold(text).lower()
            if not key or key in seen:
                continue
            if len(text) < MIN_SYNONYM_CHARS:
                raise ValueError(
                    "targets/%s: synonym %r is shorter than %d characters; "
                    "matching has no word boundary, so it is refused rather "
                    "than allowed to fire inside unrelated genotypes"
                    % (target_id, text, MIN_SYNONYM_CHARS)
                )
            seen.add(key)
            pairs.append((target_id, text))

    by_key = [(fold(s).lower(), tid, s) for tid, s in pairs]
    for key, tid, synonym in by_key:
        for other_key, other_tid, other_synonym in by_key:
            if other_tid == tid:
                continue
            if key in other_key:
                raise ValueError(
                    "targets: synonym %r of %r is a substring of synonym %r of "
                    "%r; one occurrence would resolve to two genes"
                    % (synonym, tid, other_synonym, other_tid)
                )

    pairs.sort(key=lambda pair: (-len(pair[1]), pair[0], pair[1]))
    return SynonymIndex(pairs=tuple(pairs))


def occurrences(text: Optional[str], synonym: str) -> List[Tuple[int, int]]:
    """Every ``(start, end)`` at which ``synonym`` occurs in ``text``.

    Case-insensitive, dash- and space-folded, and **without a word boundary**:
    the symbol is routinely glued to the Cre driver and the floxed allele
    (``AgRPp53``, ``Sh2b1f/f``, ``Htr2cfl/Y``). The single guard is that
    neither neighbouring character may be a digit, which is what keeps
    ``Raf1`` out of ``Raf12``. Offsets index ``text`` itself, because the fold
    is length-preserving.
    """
    hay = fold(text)
    needle = fold(synonym)
    if not hay or len(needle) < MIN_SYNONYM_CHARS:
        return []
    out: List[Tuple[int, int]] = []
    for match in re.finditer(re.escape(needle), hay, re.IGNORECASE):
        start, end = match.span()
        if start > 0 and hay[start - 1].isdigit():
            continue
        if end < len(hay) and hay[end].isdigit():
            continue
        out.append((start, end))
    return out


@dataclass(frozen=True)
class TargetMatch:
    """One occurrence of one vocabulary synonym in one curated field."""

    target_id: str
    raw: str
    synonym: str
    source_field: str
    start: int


#: Resolution outcomes, all of them named rather than signalled by ``None``.
RESOLVED = "resolved"
NO_MATCH = "no_match"
AMBIGUOUS = "ambiguous"
NOT_ELIGIBLE = "not_eligible"
ALREADY_SET = "already_set"


@dataclass(frozen=True)
class TargetResolution:
    """What the resolver concluded about one record."""

    status: str
    target_raw: Optional[str] = None
    target_id: Optional[str] = None
    candidates: Tuple[str, ...] = ()
    matches: Tuple[TargetMatch, ...] = ()

    @property
    def is_resolved(self) -> bool:
        return self.status == RESOLVED


def resolve_target(record: TierARecord, index: SynonymIndex) -> TargetResolution:
    """Resolve ``record``'s molecular target from its own curated text.

    See the module docstring for the rule and for what this does and does not
    prove. Never raises, never guesses, and never writes to ``record``.
    """
    if not _blank(record.target_raw) or not _blank(record.target_id):
        return TargetResolution(
            status=ALREADY_SET,
            target_raw=record.target_raw,
            target_id=record.target_id,
        )

    cls = perturbation_class_from_method(record.method)
    if cls not in TARGET_BEARING_CLASSES:
        return TargetResolution(status=NOT_ELIGIBLE)

    matches: List[TargetMatch] = []
    for position, field_name in enumerate(TARGET_SEARCH_FIELDS):
        text = getattr(record, field_name, None)
        if _blank(text):
            continue
        for target_id, synonym in index.pairs:
            for start, end in occurrences(text, synonym):
                matches.append(
                    TargetMatch(
                        target_id=target_id,
                        raw=str(text)[start:end],
                        synonym=synonym,
                        source_field=field_name,
                        start=position * 10 ** 6 + start,
                    )
                )

    found = sorted({match.target_id for match in matches})
    if not found:
        return TargetResolution(status=NO_MATCH)
    if len(found) > 1:
        return TargetResolution(
            status=AMBIGUOUS,
            candidates=tuple(found),
            matches=tuple(matches),
        )

    target_id = found[0]
    # Longest match wins, then the earlier field, then the earlier position.
    best = min(
        (m for m in matches if m.target_id == target_id),
        key=lambda m: (-len(m.raw), m.start),
    )
    return TargetResolution(
        status=RESOLVED,
        target_raw=best.raw,
        target_id=target_id,
        candidates=(target_id,),
        matches=tuple(matches),
    )


def _blank(value) -> bool:
    return value is None or not str(value).strip()


# --------------------------------------------------------------------------- #
# migrating
# --------------------------------------------------------------------------- #


@dataclass
class MigrationResult:
    """Everything one migration run produced, and everything it counted."""

    records: List[TierARecord] = field(default_factory=list)
    id_map: Dict[str, str] = field(default_factory=dict)
    counts: Dict[str, int] = field(default_factory=dict)
    resolutions: Dict[str, TargetResolution] = field(default_factory=dict)
    ambiguous: List[Tuple[str, Tuple[str, ...], str]] = field(default_factory=list)
    unresolved: List[Tuple[str, str, str]] = field(default_factory=list)
    unmapped_conflicts: List[Tuple[str, str]] = field(default_factory=list)
    stale_input_ids: List[str] = field(default_factory=list)
    id_collisions: List[Tuple[str, List[str]]] = field(default_factory=list)

    @property
    def needs_a_human(self) -> bool:
        return bool(self.ambiguous or self.unmapped_conflicts)


def migrate_record(
    record: TierARecord, index: SynonymIndex
) -> Tuple[TierARecord, TargetResolution]:
    """Migrate one record. Returns ``(new_record, target_resolution)``.

    Pure: ``record`` is not modified. ``conflicts_with`` is left as it is here
    and rewritten by :func:`remap_conflicts` once the whole map is known.

    Idempotent by construction. ``perturbation_class`` and
    ``perturbation_polarity`` are recomputed from ``method`` and ``valence``
    every time, which cannot lose information because they are pure functions
    of those fields; the target is only filled when empty; the condition
    fields are never touched. Applying this to its own output changes nothing.
    """
    data = record.to_dict()
    data["schema_version"] = TO_VERSION

    cls = perturbation_class_from_method(record.method)
    polarity = perturbation_polarity_from(record.method, record.valence)
    data["perturbation_class"] = cls.value if cls is not None else None
    data["perturbation_polarity"] = polarity.value if polarity is not None else None

    resolution = resolve_target(record, index)
    if resolution.is_resolved:
        data["target_raw"] = resolution.target_raw
        data["target_id"] = resolution.target_id

    for name in CONDITION_FIELDS_STILL_IN_THE_SCHEMA:
        # Explicit, not incidental: the value the curators left is the value
        # that stays, and for all 362 records of the first corpus that is
        # ``None``. Writing ``not_reported`` here would claim a paper was
        # checked (SCHEMA.md 3.3b).
        data[name] = getattr(record, name)

    data["record_id"] = None  # recomputed by TierARecord.__post_init__
    migrated = TierARecord.from_dict(data)
    return migrated, resolution


def remap_conflicts(
    records: Sequence[TierARecord], id_map: Dict[str, str]
) -> List[Tuple[str, str]]:
    """Rewrite every ``conflicts_with`` reference through ``id_map``, in place.

    Returns the ``(record_id, reference)`` pairs that were **not** in the map
    and were therefore left untouched. An unmapped reference is reported, never
    dropped: dropping it would turn a dangling pointer into a silent absence of
    conflict, which is the more expensive of the two failures.
    """
    unmapped: List[Tuple[str, str]] = []
    for record in records:
        if not record.conflicts_with:
            continue
        rewritten: List[str] = []
        for reference in record.conflicts_with:
            if reference in id_map:
                rewritten.append(id_map[reference])
            else:
                rewritten.append(reference)
                unmapped.append((record.record_id or "", reference))
        record.conflicts_with = rewritten
    return unmapped


def migrate_records(
    records: Sequence[TierARecord], index: SynonymIndex
) -> MigrationResult:
    """Migrate a whole corpus: fields, ids, the id map and ``conflicts_with``.

    Counts everything it did. Detects, rather than tolerates, the two ways this
    could corrupt the dataset: an input id that is not the one v0.1.0 would
    have computed (so the file was edited after its ids were written), and two
    records colliding on one new id (so two observations would become one).
    """
    result = MigrationResult()
    counts = {
        "records_in": len(records),
        "records_out": 0,
        "ids_changed": 0,
        "ids_unchanged": 0,
        "schema_version_set": 0,
        "perturbation_class_filled": 0,
        "perturbation_class_null": 0,
        "perturbation_polarity_filled": 0,
        "perturbation_polarity_null": 0,
        "target_eligible": 0,
        "target_resolved": 0,
        "target_unresolved_no_match": 0,
        "target_ambiguous": 0,
        "target_not_eligible": 0,
        "target_already_set": 0,
        "conflicts_with_references_in": 0,
        "conflicts_with_references_remapped": 0,
        "conflicts_with_references_unmapped": 0,
        "records_already_at_target_version": 0,
        "input_ids_not_reproducible_under_0_1_0": 0,
    }

    by_new_id: Dict[str, List[str]] = {}
    for record in records:
        old_id = record.record_id or ""
        if record.schema_version == TO_VERSION:
            counts["records_already_at_target_version"] += 1
        elif (
            record.schema_version == FROM_VERSION
            and compute_record_id_0_1_0(record) != old_id
        ):
            counts["input_ids_not_reproducible_under_0_1_0"] += 1
            result.stale_input_ids.append(old_id)

        counts["conflicts_with_references_in"] += len(record.conflicts_with)

        migrated, resolution = migrate_record(record, index)
        result.records.append(migrated)
        result.resolutions[old_id] = resolution

        counts["schema_version_set"] += 1
        if migrated.perturbation_class is None:
            counts["perturbation_class_null"] += 1
        else:
            counts["perturbation_class_filled"] += 1
        if migrated.perturbation_polarity is None:
            counts["perturbation_polarity_null"] += 1
        else:
            counts["perturbation_polarity_filled"] += 1

        if resolution.status == NOT_ELIGIBLE:
            counts["target_not_eligible"] += 1
        else:
            if resolution.status != ALREADY_SET:
                counts["target_eligible"] += 1
            if resolution.status == RESOLVED:
                counts["target_resolved"] += 1
            elif resolution.status == NO_MATCH:
                counts["target_unresolved_no_match"] += 1
                result.unresolved.append(
                    (
                        migrated.record_id or "",
                        migrated.method.value if migrated.method else "null",
                        (migrated.targeting_strategy or "")[:120],
                    )
                )
            elif resolution.status == AMBIGUOUS:
                counts["target_ambiguous"] += 1
                result.ambiguous.append(
                    (
                        migrated.record_id or "",
                        resolution.candidates,
                        (migrated.targeting_strategy or "")[:120],
                    )
                )
            elif resolution.status == ALREADY_SET:
                counts["target_already_set"] += 1

        new_id = migrated.record_id or ""
        result.id_map[old_id] = new_id
        by_new_id.setdefault(new_id, []).append(old_id)
        if new_id == old_id:
            counts["ids_unchanged"] += 1
        else:
            counts["ids_changed"] += 1

    result.id_collisions = [
        (new_id, olds) for new_id, olds in sorted(by_new_id.items()) if len(olds) > 1
    ]

    result.unmapped_conflicts = remap_conflicts(result.records, result.id_map)
    counts["conflicts_with_references_unmapped"] = len(result.unmapped_conflicts)
    counts["conflicts_with_references_remapped"] = (
        counts["conflicts_with_references_in"]
        - counts["conflicts_with_references_unmapped"]
    )
    counts["records_out"] = len(result.records)
    result.counts = counts
    return result


# --------------------------------------------------------------------------- #
# counting the dataset (every number the report prints)
# --------------------------------------------------------------------------- #


def _dist(values: Iterable[Optional[str]]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for value in values:
        key = "null" if value is None else str(getattr(value, "value", value))
        out[key] = out.get(key, 0) + 1
    return out


def dataset_stats(
    records: Sequence[TierARecord], index: SynonymIndex
) -> Dict[str, object]:
    """Recompute, from records alone, every quantity the report states.

    This is what makes the report checkable (``docs/CURATION_PROTOCOL.md`` 6):
    the ``report`` subcommand calls this on the committed dataset and must
    reproduce the numbers ``run`` published. It also re-runs the target
    resolver and compares it against the stored ``target_id``, so a hand-edited
    target shows up as a disagreement rather than being assumed correct.
    """
    stats: Dict[str, object] = {}
    stats["records"] = len(records)
    stats["schema_versions"] = _dist(r.schema_version for r in records)
    stats["method"] = _dist(r.method for r in records)
    stats["valence"] = _dist(r.valence for r in records)
    stats["perturbation_class"] = _dist(r.perturbation_class for r in records)
    stats["perturbation_polarity"] = _dist(r.perturbation_polarity for r in records)

    cross: Dict[str, Dict[str, int]] = {}
    for record in records:
        cls = "null" if record.perturbation_class is None else record.perturbation_class.value
        pol = (
            "null"
            if record.perturbation_polarity is None
            else record.perturbation_polarity.value
        )
        cross.setdefault(cls, {})
        cross[cls][pol] = cross[cls].get(pol, 0) + 1
    stats["class_by_polarity"] = cross

    # --- targets --------------------------------------------------------
    eligible = [
        r
        for r in records
        if perturbation_class_from_method(r.method) in TARGET_BEARING_CLASSES
    ]
    stats["target_eligible_records"] = len(eligible)
    stats["target_resolved_records"] = sum(1 for r in eligible if not _blank(r.target_id))
    stats["target_unresolved_records"] = sum(1 for r in eligible if _blank(r.target_id))
    stats["target_on_ineligible_records"] = sum(
        1
        for r in records
        if perturbation_class_from_method(r.method) not in TARGET_BEARING_CLASSES
        and (not _blank(r.target_id) or not _blank(r.target_raw))
    )
    stats["target_id"] = _dist(r.target_id for r in eligible)

    recheck = {"agree": 0, "disagree": [], "ambiguous": [], "no_match": 0}
    for record in eligible:
        probe = TierARecord.from_dict(dict(record.to_dict(), target_raw=None, target_id=None))
        resolution = resolve_target(probe, index)
        if resolution.status == AMBIGUOUS:
            recheck["ambiguous"].append(
                (
                    record.record_id or "",
                    resolution.candidates,
                    (record.targeting_strategy or "")[:120],
                )
            )
            continue
        if resolution.status == NO_MATCH:
            recheck["no_match"] += 1
            if not _blank(record.target_id):
                recheck["disagree"].append(
                    (record.record_id or "", record.target_id, "no_match")
                )
            continue
        if resolution.target_id == record.target_id:
            recheck["agree"] += 1
        else:
            recheck["disagree"].append(
                (record.record_id or "", record.target_id, resolution.target_id)
            )
    stats["target_recheck"] = recheck

    # --- gene / study / lab collinearity --------------------------------
    genes: Dict[str, Dict[str, object]] = {}
    for record in eligible:
        if _blank(record.target_id):
            continue
        slot = genes.setdefault(
            record.target_id, {"records": 0, "studies": set(), "labs": set()}
        )
        slot["records"] = int(slot["records"]) + 1
        slot["studies"].add(record.study_id or "")
        slot["labs"].add(
            (record.last_author or "", record.last_author_affiliation or "")
        )
    stats["genes"] = {
        gene: {
            "records": slot["records"],
            "studies": len(slot["studies"]),
            "labs": len(slot["labs"]),
            "study_ids": sorted(slot["studies"]),
        }
        for gene, slot in sorted(
            genes.items(), key=lambda kv: (-int(kv[1]["records"]), kv[0])
        )
    }
    all_studies = set()
    all_labs = set()
    for slot in genes.values():
        all_studies |= slot["studies"]
        all_labs |= slot["labs"]
    stats["genes_total"] = len(genes)
    stats["gene_studies_total"] = len(all_studies)
    stats["gene_labs_total"] = len(all_labs)
    stats["genes_in_exactly_one_study"] = sum(
        1 for slot in genes.values() if len(slot["studies"]) == 1
    )
    stats["genes_in_exactly_one_lab"] = sum(
        1 for slot in genes.values() if len(slot["labs"]) == 1
    )
    studies_to_genes: Dict[str, set] = {}
    for gene, slot in genes.items():
        for study in slot["studies"]:
            studies_to_genes.setdefault(study, set()).add(gene)
    stats["studies_with_more_than_one_gene"] = {
        study: sorted(gs) for study, gs in sorted(studies_to_genes.items()) if len(gs) > 1
    }

    # --- conditions ------------------------------------------------------
    condition = {}
    for name in ("sex", "diet", "age_group") + CONDITION_FIELDS_THE_MIGRATION_MAY_NOT_FILL:
        # ``getattr`` with a default, for the two fields v0.3.0 retired: they
        # were null on every record of the v0.2.0 corpus, so counting them as
        # null is what the v0.2.0 report said and still true of the data. The
        # alternative -- dropping them from the count -- would silently change
        # the published slot total.
        values = [getattr(r, name, None) for r in records]
        condition[name] = {
            "set": sum(1 for v in values if v is not None),
            "null": sum(1 for v in values if v is None),
            "not_reported": sum(
                1 for v in values if getattr(v, "value", v) == "not_reported"
            ),
        }
    stats["condition_fields"] = condition
    stats["condition_field_slots_unfilled"] = sum(
        condition[name]["null"]
        for name in CONDITION_FIELDS_THE_MIGRATION_MAY_NOT_FILL
    )

    # --- validation ------------------------------------------------------
    hard: Dict[str, int] = {}
    warn: Dict[str, int] = {}
    clean = 0
    for record in records:
        findings = validate(record)
        # Version artefacts are set aside here for the same reason they are set
        # aside in every before/after comparison (:func:`_non_artefact`): this
        # module writes a file at a pinned TO_VERSION, so once the schema moved
        # past it, S005/S006 fire on every output record and say only that the
        # validator is newer than the file. Counting them here made the figure
        # "N of M records carry no hard finding" collapse to zero, and
        # contradict the sentence printed immediately after it, which sets the
        # same two codes aside.
        if not any(
            f.severity is Severity.HARD and f.code not in _VERSION_ARTEFACT_CODES
            for f in findings
        ):
            clean += 1
        for finding in findings:
            bucket = hard if finding.severity is Severity.HARD else warn
            bucket[finding.code] = bucket.get(finding.code, 0) + 1
    stats["hard_findings"] = hard
    stats["warning_findings"] = warn
    #: Records with no hard finding **other than** a version artefact.
    stats["records_hard_clean"] = clean

    stats["conflicts_with_references"] = sum(len(r.conflicts_with) for r in records)
    stats["record_ids_unique"] = len({r.record_id for r in records})
    stats["record_ids_consistent"] = sum(1 for r in records if r.id_is_consistent())
    return stats


# --------------------------------------------------------------------------- #
# the report
# --------------------------------------------------------------------------- #


def _table(rows: Dict[str, int], header: Tuple[str, str], total: Optional[int] = None):
    out = ["| %s | %s | %s |" % (header[0], header[1], "share")]
    out.append("|---|---|---|")
    denominator = total if total else max(1, sum(rows.values()))
    for key, count in sorted(rows.items(), key=lambda kv: (-kv[1], kv[0])):
        out.append("| `%s` | %d | %.1f%% |" % (key, count, 100.0 * count / denominator))
    return out


def render_report_md(
    stats: Dict[str, object],
    *,
    counts: Optional[Dict[str, int]] = None,
    baseline: Optional[Dict[str, object]] = None,
    result: Optional[MigrationResult] = None,
    records_path: str,
    id_map_path: str,
    backup_path: Optional[str],
    index: SynonymIndex,
    generated: str,
    recomputed_only: bool = False,
    event_complete: bool = True,
) -> str:
    """Render ``reports/schema_migration_0.2.0.md``.

    Every number comes from ``stats`` (recomputed from the dataset by
    :func:`dataset_stats`) or from ``counts`` (counted during the run and
    re-derivable from the committed id map and the preserved input). Nothing is
    estimated and nothing is carried over from a previous run.

    ``event_complete`` is ``False`` when the caller could only see the
    committed dataset and not the migration event -- then the few rows that
    only the event can answer are omitted rather than guessed, which makes the
    render a line-subsequence of the full report and so still checkable
    against it.
    """
    L: List[str] = []
    A = L.append
    A("# Schema migration: curated records v0.1.0 -> v0.2.0")
    A("")
    if recomputed_only:
        A(
            "Recomputed %s by `python3 -m syn4rch.curate.migrate report` from the "
            "committed dataset. No record was written." % generated
        )
    else:
        A("Generated %s by `python3 -m syn4rch.curate.migrate run`." % generated)
    A("")
    A(
        "Every number below is recomputed from `%s` by "
        "`python3 -m syn4rch.curate.migrate report` "
        "(`docs/CURATION_PROTOCOL.md` 6). Nothing here is estimated."
        % os.path.relpath(records_path, _repo_root())
    )
    A("")

    # --- what moved -------------------------------------------------------
    A("## 1. What moved")
    A("")
    A("| | |")
    A("|---|---|")
    if counts:
        A("| records in | %d |" % counts["records_in"])
        A("| records out | %d |" % counts["records_out"])
        A("| records dropped or added | %d |" % (counts["records_in"] - counts["records_out"]))
        A("| `record_id` changed | %d |" % counts["ids_changed"])
        A("| `record_id` unchanged | %d |" % counts["ids_unchanged"])
        A(
            "| `conflicts_with` references rewritten | %d of %d |"
            % (
                counts["conflicts_with_references_remapped"],
                counts["conflicts_with_references_in"],
            )
        )
        A(
            "| `conflicts_with` references not in the id map | %d |"
            % counts["conflicts_with_references_unmapped"]
        )
    else:
        A("| records in the committed dataset | %d |" % stats["records"])
        A("| `conflicts_with` references | %d |" % stats["conflicts_with_references"])
    A("| distinct `record_id` | %d |" % stats["record_ids_unique"])
    A(
        "| `record_id` consistent with its own fields | %d of %d |"
        % (stats["record_ids_consistent"], stats["records"])
    )
    A("| id map | `%s` |" % os.path.relpath(id_map_path, _repo_root()))
    if backup_path:
        A("| pre-migration file kept at | `%s` |" % os.path.relpath(backup_path, _repo_root()))
    A("")
    A(
        "`RECORD_ID_FIELDS` grew from %d components to %d, so **every id in the "
        "corpus changed**. That is the whole reason for the version bump. The "
        "five components added were %s."
        % (
            len(RECORD_ID_FIELDS_0_1_0),
            len(RECORD_ID_FIELDS),
            ", ".join("`%s`" % f for f in RECORD_ID_FIELDS_ADDED_IN_0_2_0),
        )
    )
    A("")
    A("The map from old id to new id is complete and one-to-one over the input:")
    A("every v0.1.0 id appears exactly once as a key, and no two of them share a")
    A("new id. Two records colliding on one new id would have merged two")
    A("observations into one, so the migration refuses to write at all rather")
    A("than writing the collision and mentioning it here.")
    A("")
    A("Schema versions now present in the file: %s." % _inline(stats["schema_versions"]))
    A("")

    # --- the perturbation axis -------------------------------------------
    A("## 2. The perturbation axis")
    A("")
    A("Filled by rule from fields the curators already wrote:")
    A("`perturbation_class` from `method`, `perturbation_polarity` from")
    A("`(method, valence)`. Both are pure functions")
    A("(`syn4rch.schema.rules`), so nothing here is a judgement.")
    A("")
    A("### `perturbation_class`")
    A("")
    L.extend(_table(stats["perturbation_class"], ("class", "records"), int(stats["records"])))
    A("")
    A("### `perturbation_polarity`")
    A("")
    L.extend(_table(stats["perturbation_polarity"], ("polarity", "records"), int(stats["records"])))
    A("")
    A("### Polarity within class")
    A("")
    A("Polarity read without its class is meaningless (SCHEMA.md 3.4b):")
    A("`decrease` is less firing under `circuit_activity` and less of the target")
    A("gene product under `molecular`.")
    A("")
    A("| class | polarity | records |")
    A("|---|---|---|")
    cross = stats["class_by_polarity"]
    for cls in sorted(cross, key=lambda c: -sum(cross[c].values())):
        for pol, count in sorted(cross[cls].items(), key=lambda kv: (-kv[1], kv[0])):
            A("| `%s` | `%s` | %d |" % (cls, pol, count))
    A("")
    A(
        "For reference, the `method` and `valence` these were derived from: "
        "method %s; valence %s."
        % (_inline(stats["method"]), _inline(stats["valence"]))
    )
    A("")

    # --- targets ----------------------------------------------------------
    A("## 3. Molecular targets")
    A("")
    A(
        "The target was never missing from this dataset; it was curated inside "
        "the free text of `targeting_strategy` and sometimes `actuator`. This "
        "migration lifted it into `target_raw` / `target_id` against a new "
        "vocabulary, `data/vocab/targets.json` (%d entries, %d matchable "
        "strings). **It did not read any paper.** A record's target is only as "
        "good as its `targeting_strategy`: a curator's wrong gene there becomes "
        "a confidently resolved wrong gene here, and nothing in this pipeline "
        "would notice."
        % (len(index.ids), len(index))
    )
    A("")
    A("| | |")
    A("|---|---|")
    A(
        "| eligible records (`molecular` or `pharmacological`) | %d |"
        % stats["target_eligible_records"]
    )
    A("| resolved | %d |" % stats["target_resolved_records"])
    A("| unresolved | %d |" % stats["target_unresolved_records"])
    if counts and event_complete:
        # Only when the migration event itself was observed (a ``run``, or a
        # ``report`` that re-read the preserved input). From the committed
        # dataset alone, an unresolved record does not say which of the two
        # reasons left it unresolved.
        A("| of which, no synonym matched | %d |" % counts["target_unresolved_no_match"])
        A("| of which, ambiguous (two or more genes matched) | %d |" % counts["target_ambiguous"])
        A("| already curated, left untouched | %d |" % counts["target_already_set"])
    A(
        "| ineligible records carrying a target (must be 0; `S008`) | %d |"
        % stats["target_on_ineligible_records"]
    )
    A("")
    recheck = stats["target_recheck"]
    A(
        "Re-running the resolver over the committed file reproduces the stored "
        "`target_id` for %d of %d resolved records, with %d disagreement(s) and "
        "%d ambiguity(ies). A disagreement would mean the stored target no "
        "longer follows from the text it was taken from."
        % (
            recheck["agree"],
            stats["target_resolved_records"],
            len(recheck["disagree"]),
            len(recheck["ambiguous"]),
        )
    )
    A("")
    if result is not None and result.ambiguous:
        A("### Ambiguous matches, for a human to curate")
        A("")
        A("Two or more genes occur in one record's curated text. The migration")
        A("left `target_raw` and `target_id` null rather than picking one.")
        A("")
        A("| record_id | candidates | targeting_strategy (truncated) |")
        A("|---|---|---|")
        for record_id, candidates, text in result.ambiguous:
            A(
                "| `%s` | %s | %s |"
                % (
                    record_id,
                    ", ".join("`%s`" % c for c in candidates),
                    str(text).replace("|", "/"),
                )
            )
        A("")
    else:
        A("### Ambiguous matches, for a human to curate")
        A("")
        A("None: no record's curated text names two different genes of the")
        A("vocabulary. Had one, it would be listed here by `record_id` and left")
        A("unresolved; an ambiguous match is never guessed.")
        A("")
    if result is not None and result.unresolved:
        A("### Eligible but unresolved")
        A("")
        A("| record_id | method | targeting_strategy (truncated) |")
        A("|---|---|---|")
        for record_id, method, text in result.unresolved:
            A("| `%s` | `%s` | %s |" % (record_id, method, text.replace("|", "/")))
        A("")
        A(
            "`data/vocab/targets.json` is a vocabulary of *genes*, so a locally "
            "infused drug resolves to nothing in it. Naming the receptor a drug "
            "is thought to act on would be inference about pharmacology, not "
            "resolution of a curated string, and this code does not do it. A "
            "pharmacological record's compound stays where the curator put it, "
            "in `ligand`."
        )
        A("")
    A("### Targets found")
    A("")
    A("| target_id | records | studies | labs |")
    A("|---|---|---|---|")
    for gene, slot in stats["genes"].items():
        A("| `%s` | %d | %d | %d |" % (gene, slot["records"], slot["studies"], slot["labs"]))
    A("")

    # --- collinearity -----------------------------------------------------
    A("## 4. The gene is the study is the lab")
    A("")
    A(
        "**The %d target genes in this corpus come from %d studies and %d "
        "laboratories, and %d of the %d appear in exactly one study while %d "
        "appear in exactly one laboratory: the gene, the study and the lab are "
        "very nearly the same variable.** A model told the gene has therefore "
        "been told the study and the lab, so a cross-validated number that "
        "does not group on study or lab is measuring how well %d single-gene "
        "experiments can be memorised, not whether a molecular perturbation "
        "predicts a physiological outcome."
        % (
            stats["genes_total"],
            stats["gene_studies_total"],
            stats["gene_labs_total"],
            stats["genes_in_exactly_one_study"],
            stats["genes_total"],
            stats["genes_in_exactly_one_lab"],
            stats["genes_total"],
        )
    )
    A("")
    A(
        "This is the single biggest structural fact about the molecular half of "
        "the corpus, and the reason `target_id` had to exist: without it, "
        '"molecular perturbation of this population" pools every one of those '
        "experiments into one cell, and the cell then looks like a replicated "
        "finding. With it, the collinearity is visible instead of hidden -- it "
        "is not removed. Nothing a schema can do makes one laboratory's single "
        "gene generalise."
    )
    A("")
    multi = stats["studies_with_more_than_one_gene"]
    if multi:
        A("Studies contributing more than one gene:")
        A("")
        for study, gs in multi.items():
            A("* `%s`: %s" % (study, ", ".join("`%s`" % g for g in gs)))
    else:
        A("No study in the corpus perturbs more than one gene.")
    A("")
    A("Leakage control is specified in SCHEMA.md 5: group on `study_id` and on")
    A("lab identity. Under this corpus those groupings also group the gene, so")
    A("leave-one-study-out on molecular records is, in effect,")
    A("leave-one-gene-out, and no split of this dataset can test whether a")
    A("molecular perturbation of one gene predicts another.")
    A("")

    # --- conditions -------------------------------------------------------
    A("## 5. The condition dimension was not filled, and could not be")
    A("")
    A(
        "**The migration was not permitted to fill `challenge`, "
        "`challenge_detail`, `condition_source` or `condition_locator`, and did "
        "not fill them.** Those four require reading the paper, and this code "
        "read no paper. It did not write `not_reported` either: per SCHEMA.md "
        "3.3b that means a curator checked the paper and it does not say, so "
        "writing it from a migration would record a curation act that never "
        "happened (SCHEMA.md 2 rule 4). They are left `null` and counted here."
    )
    A("")
    A("| field | set | null | of which `not_reported` |")
    A("|---|---|---|---|")
    for name, slot in stats["condition_fields"].items():
        A(
            "| `%s` | %d | %d | %d |"
            % (name, slot["set"], slot["null"], slot["not_reported"])
        )
    A("")
    A(
        "%d of the %d condition-field slots the v0.2.0 schema added are "
        "unfilled (%d records x %d fields). `sex`, `diet` and `age_group` are "
        "curated where known and were not touched."
        % (
            stats["condition_field_slots_unfilled"],
            len(CONDITION_FIELDS_THE_MIGRATION_MAY_NOT_FILL) * int(stats["records"]),
            stats["records"],
            len(CONDITION_FIELDS_THE_MIGRATION_MAY_NOT_FILL),
        )
    )
    A("")
    A(
        "One consequence is worth stating in advance, because it will look like "
        "a defect when it happens. `challenge` is a `record_id` component, so "
        "the curation pass that fills it will change the id of every record it "
        "touches, exactly as this migration changed all %d. That is the "
        "identity model working as designed -- a fasted cohort and an "
        "ad-libitum one are two observations -- but it means condition "
        "curation is an id-changing act and needs the same map-and-remap "
        "treatment this module implements, not an in-place edit."
        % stats["records"]
    )
    A("")
    A(
        "The consequence is in the dataset's own warnings: `W015` fires on the "
        "%d records for which none of `(sex, diet, age_group, challenge)` is "
        "known, and those records can only ever be compared with other "
        "condition-blind records. Closing that gap is curation work, not "
        "migration work, and no amount of migration will do it."
        % stats["warning_findings"].get("W015_CONDITION_UNSPECIFIED", 0)
    )
    A("")

    # --- validation -------------------------------------------------------
    A("## 6. Validation, before and after")
    A("")
    A(
        "The migration refuses to write if the result would carry a hard "
        "validation finding the input did not already have, and the two tables "
        "below are how that is checked. `S005` (id inconsistent) and `S006` "
        "(version mismatch) are expected on every *input* record, because at "
        "v0.2.0 the seven-component ids no longer hash and the old version "
        "string no longer matches: they are the artefacts the migration exists "
        "to remove, and they are excluded when the two sides are compared."
    )
    A("")
    A("### After the migration, recomputed from the committed dataset")
    A("")
    A("| finding | severity | records |")
    A("|---|---|---|")
    for code, count in sorted(stats["hard_findings"].items()):
        A("| `%s` | HARD | %d |" % (code, count))
    for code, count in sorted(stats["warning_findings"].items()):
        A("| `%s` | WARNING | %d |" % (code, count))
    A("")
    A(
        "%d of %d records carry no hard finding."
        % (stats["records_hard_clean"], stats["records"])
    )
    A("")
    if baseline is not None:
        A("### Before the migration, recomputed from the preserved input")
        A("")
        A(
            "Counted by running `syn4rch.schema.validate` over the v%s file at "
            "the v%s code's standards. Checkable: the pre-migration file is "
            "kept, and `migrate report` re-reads it." % (FROM_VERSION, TO_VERSION)
        )
        A("")
        A("| finding | severity | records |")
        A("|---|---|---|")
        for code, count in sorted(baseline["hard_findings"].items()):
            A("| `%s` | HARD | %d |" % (code, count))
        for code, count in sorted(baseline["warning_findings"].items()):
            A("| `%s` | WARNING | %d |" % (code, count))
        A("")
        before_other = _non_artefact(baseline["hard_findings"])
        after_other = _non_artefact(stats["hard_findings"])
        A(
            "%d of %d input records carried no hard finding at all, and %d of "
            "%d do now."
            % (
                baseline["records_hard_clean"],
                baseline["records"],
                stats["records_hard_clean"],
                stats["records"],
            )
        )
        A("")
        A(
            "Setting `S005` and `S006` aside, the input carried %s and the "
            "migrated dataset carries %s. %s"
            % (
                _findings_phrase(before_other),
                _findings_phrase(after_other),
                (
                    "The migration therefore introduced no hard finding, which "
                    "is the condition `migrate run` checks before it will "
                    "write anything."
                    if not _positive_deltas(before_other, after_other)
                    else "The migration introduced at least one hard finding, "
                    "and `migrate run` was told explicitly to write anyway "
                    "(`--allow-new-findings`)."
                ),
            )
        )
        A("")
    A(
        "`W014` (molecular perturbation with no identified target) is the one "
        "warning this migration was able to retire, and it did so by resolving "
        "a field the curators had already written, not by learning anything new "
        "about any paper."
    )
    A("")
    A(
        "**This report describes the dataset as it stood at v0.2.0.** The "
        "v0.3.0 condition fill has since rewritten `data/curated/records.jsonl` "
        "and re-identified 227 of its 362 records, so `migrate report` can no "
        "longer recompute these numbers from the committed dataset: it refuses, "
        "naming the version mismatch. Recomputing them needs the preserved "
        "v0.2.0 state (`data/curated/records.0.2.0.jsonl.bak`, untracked like "
        "every backup), and `data/curated/record_id_map_0.2.0_to_0.3.0.json` "
        "carries any id in this report forward to the record it is now. This is "
        "the ordinary consequence of a second migration, and it is stated here "
        "rather than left for a reader to discover from an error message."
    )
    A("")
    A(
        "**Two notes on reading the counts above.** They are what the *current* "
        "validator sees in these two files, because a report whose numbers "
        "cannot be recomputed from the committed dataset does not belong in "
        "this repository (`docs/CURATION_PROTOCOL.md` 6). So a rule introduced "
        "after this migration appears in its table: `W017`, added by schema "
        "v0.3.0, counts condition values that carry no verbatim sentence, and "
        "it describes a gap this migration never claimed to close."
    )
    A("")
    A(
        "And the \"carried no hard finding at all\" figure sets `S005` and "
        "`S006` aside, as every other before/after comparison here does. An "
        "earlier rendering did not, which made it read **0 of 362** for the "
        "input -- a number produced entirely by a 0.1.0 file being validated "
        "by newer code, and one that contradicted the sentence immediately "
        "below it. The migration neither introduced nor removed a hard "
        "finding, and both columns now say so."
    )
    A("")

    # --- what was filled --------------------------------------------------
    A("## 7. What the migration filled, and what it refused to")
    A("")
    A("| v0.2.0 field | filled by the migration? | how |")
    A("|---|---|---|")
    A("| `perturbation_class` | yes | total function of `method` |")
    A("| `perturbation_polarity` | yes | function of `(method, valence)` |")
    A("| `target_raw` | only when exactly one gene occurs in the curated text | resolution of `targeting_strategy` / `actuator` |")
    A("| `target_id` | as `target_raw` | `data/vocab/targets.json` |")
    A("| `challenge` | **no** | requires reading the paper |")
    A("| `challenge_detail` | **no** | requires reading the paper |")
    A("| `condition_source` | **no** | requires reading the paper |")
    A("| `condition_locator` | **no** | requires reading the paper |")
    A("| `schema_version` | yes | set to `%s` |" % TO_VERSION)
    A("| `record_id` | yes | recomputed over 12 components |")
    A("| `conflicts_with` | yes | rewritten through the id map |")
    A("| everything else | no | untouched |")
    A("")
    A(
        "`docs/DECISIONS.md` D-017 makes the gate the only writer of the curated "
        "file. This migration is the one other writer, and a narrower one: it "
        "creates no record, admits none, drops none, and writes no value that a "
        "curator did not already write or that a published rule does not derive "
        "from one. D-017 should be amended to name the exception."
    )
    A("")
    if result is not None and result.stale_input_ids:
        A("## 8. Input ids that did not reproduce under v0.1.0")
        A("")
        A(
            "%d input record(s) carried an id that the v0.1.0 hash of their own "
            "fields does not reproduce, which means the file was edited after "
            "its ids were written. Listed so that the edit can be found:"
            % len(result.stale_input_ids)
        )
        A("")
        for record_id in result.stale_input_ids:
            A("* `%s`" % record_id)
        A("")
    return "\n".join(L) + "\n"


#: Hard findings that a version difference alone explains, on either side of
#: the migration: every v0.1.0 input record carries both when read by code
#: implementing a later schema, and every v0.2.0 output record carries
#: ``S006`` when read by code implementing a later schema still -- which, from
#: v0.3.0 onwards, is all of it. Excluded when the before and after validation
#: states are compared, because counting them as "fixed" would flatter the
#: migration and counting them as "introduced" would block a migration that
#: did nothing wrong.
_VERSION_ARTEFACT_CODES = (
    "S005_RECORD_ID_INCONSISTENT",
    "S006_SCHEMA_VERSION_MISMATCH",
)


def _non_artefact(findings: Dict[str, int]) -> Dict[str, int]:
    """``findings`` without the codes a version difference alone explains.

    Applied to **both** sides of every before/after comparison. It used to be
    applied to the input side only, which was invisible while
    :data:`TO_VERSION` tracked the schema package -- the output then declared
    the current version and carried no ``S006`` to subtract. With the target
    pinned to v0.2.0 the one-sided filter counted the output's version marker
    as a finding the migration had introduced, and ``migrate run`` refused to
    write at all.
    """
    return {
        code: count
        for code, count in findings.items()
        if code not in _VERSION_ARTEFACT_CODES
    }


def _findings_phrase(findings: Dict[str, int]) -> str:
    if not findings:
        return "no hard finding"
    return ", ".join(
        "%d x `%s`" % (count, code) for code, count in sorted(findings.items())
    )


def _inline(dist: Dict[str, int]) -> str:
    return ", ".join(
        "%s=%d" % (key, count)
        for key, count in sorted(dist.items(), key=lambda kv: (-kv[1], kv[0]))
    )


# --------------------------------------------------------------------------- #
# reading, writing
# --------------------------------------------------------------------------- #


#: Fields that existed at 0.2.0 and no longer exist in the schema (D-037).
#:
#: This module reads files written at :data:`TO_VERSION`, and the current
#: ``TierARecord`` refuses an unknown key by design, so without this a 0.2.0
#: file became unreadable the moment 0.3.0 removed a field -- and with it went
#: the ability to recompute ``reports/schema_migration_0.2.0.md``.
#: `docs/CURATION_PROTOCOL.md` 6 says a number that cannot be recomputed does
#: not belong in a report, so the reader is taught what was retired rather than
#: the report being left unverifiable.
#:
#: A retired field holding a **non-null** value is still an error: both were
#: null on all 362 records, which is the only reason retiring them was free.
FIELDS_RETIRED_AFTER_0_2_0 = ("condition_source", "condition_locator")


def strip_retired_fields(data: dict) -> List[str]:
    """Drop retired keys from a decoded record. Returns any that held a value.

    Targeted rather than ``strict=False``, which would also have swallowed a
    misspelled field name -- and a silently dropped key is how a column quietly
    stops being curated.
    """
    kept: List[str] = []
    for name in FIELDS_RETIRED_AFTER_0_2_0:
        if name in data:
            value = data.pop(name)
            if value is not None:
                kept.append("%s=%r" % (name, value))
    return kept


def read_records(
    path: str, expect_version: Optional[str]
) -> Tuple[List[TierARecord], List[Tuple[int, str]], Dict[str, int]]:
    """Read a curated JSONL. Returns ``(records, errors, versions_seen)``.

    ``expect_version`` is checked here rather than inside
    :func:`syn4rch.schema.io.loads_record` so that the error message can name
    every version found instead of only the first offending line.

    Keys retired after 0.2.0 are stripped (:func:`strip_retired_fields`) so
    that a file this module wrote stays readable by the code that reads it.
    """
    records: List[TierARecord] = []
    errors: List[Tuple[int, str]] = []
    versions: Dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                decoded = json.loads(line)
                if isinstance(decoded, dict):
                    held = strip_retired_fields(decoded)
                    if held:
                        errors.append(
                            (
                                number,
                                "refusing to drop a retired field that holds a "
                                "value: %s" % ", ".join(held),
                            )
                        )
                        continue
                    line = json.dumps(decoded)
                record = loads_record(
                    line, expect_version=None, path=path, line_number=number
                )
            except JsonlError as exc:
                errors.append((number, str(exc)))
                continue
            except Exception as exc:  # pragma: no cover - defensive
                errors.append((number, "%s: %s" % (type(exc).__name__, exc)))
                continue
            key = str(record.schema_version)
            versions[key] = versions.get(key, 0) + 1
            records.append(record)
    return records, errors, versions


#: Mode for a file this module creates. ``tempfile`` makes its files 0600 and
#: :func:`os.replace` keeps the temp file's mode, so without this a migrated
#: dataset would come out readable only by the user who ran the migration.
#: An existing file's own mode wins.
_NEW_FILE_MODE = 0o644


def _atomic_write_text(path: str, text: str) -> None:
    """Write ``text`` to ``path`` via a temp file in the same directory.

    The rename is :func:`os.replace`, which is atomic on POSIX, so a reader
    sees either the old file or the new one and never a half-written dataset.
    The replacement keeps the mode the file already had (:data:`_NEW_FILE_MODE`
    for a new one): a dataset that came out 0600 because of how it was written
    would be a migration artefact nobody asked for.
    """
    parent = os.path.dirname(os.path.abspath(path)) or "."
    os.makedirs(parent, exist_ok=True)
    try:
        mode = stat.S_IMODE(os.stat(path).st_mode)
    except OSError:
        mode = _NEW_FILE_MODE
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        dir=parent,
        prefix=os.path.basename(path) + ".",
        suffix=".tmp",
        delete=False,
    )
    try:
        with handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(handle.name, mode)
        os.replace(handle.name, path)
    except BaseException:
        try:
            os.unlink(handle.name)
        except OSError:  # pragma: no cover - best effort
            pass
        raise


def _records_text(records: Iterable[TierARecord]) -> str:
    return "".join(dumps_record(record) + "\n" for record in records)


def _id_map_document(
    result: MigrationResult, *, records_path: str, generated: str
) -> dict:
    return {
        "_meta": {
            "from_schema_version": FROM_VERSION,
            "to_schema_version": TO_VERSION,
            "generated": generated,
            "generated_by": "python3 -m syn4rch.curate.migrate run",
            "source": os.path.relpath(records_path, _repo_root()),
            "records": result.counts["records_in"],
            "ids_changed": result.counts["ids_changed"],
            "ids_unchanged": result.counts["ids_unchanged"],
            "record_id_fields_0_1_0": list(RECORD_ID_FIELDS_0_1_0),
            "record_id_fields_0_2_0": list(RECORD_ID_FIELDS),
            "how_to_read": (
                "map[old_id] = new_id. Complete and one-to-one over the input: "
                "every v0.1.0 record_id is a key exactly once and no two keys "
                "share a value. Use it to follow any external reference to a "
                "v0.1.0 record_id, including conflicts_with lists in older "
                "QC output, which this migration rewrote in the dataset itself."
            ),
        },
        "map": dict(sorted(result.id_map.items())),
    }


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #


def _positive_deltas(before: Dict[str, int], after: Dict[str, int]) -> Dict[str, int]:
    out = {}
    for code, count in after.items():
        if count > before.get(code, 0):
            out[code] = count - before.get(code, 0)
    return out


def cmd_run(args) -> int:
    generated = _now()
    if not os.path.isfile(args.records):
        sys.stderr.write("error: no such file: %s\n" % args.records)
        return 2

    records, errors, versions = read_records(args.records, FROM_VERSION)
    if errors:
        sys.stderr.write(
            "error: %d line(s) of %s could not be read; refusing to migrate a "
            "file this code cannot read in full\n" % (len(errors), args.records)
        )
        for number, message in errors:
            sys.stderr.write("  line %d: %s\n" % (number, message))
        return 2
    if not records:
        sys.stderr.write("error: %s holds no records\n" % args.records)
        return 2
    if set(versions) != {args.from_version}:
        sys.stderr.write(
            "error: refusing to migrate. This migration reads schema_version "
            "%r only; %s declares %s. Nothing was written.\n"
            % (
                args.from_version,
                args.records,
                ", ".join(
                    "%s (%d record%s)" % (v, n, "" if n == 1 else "s")
                    for v, n in sorted(versions.items())
                ),
            )
        )
        if set(versions) == {TO_VERSION}:
            sys.stderr.write(
                "       This file is already at %s. A second run would be a "
                "no-op; nothing to do.\n" % TO_VERSION
            )
        return 2

    try:
        index = build_synonym_index(load_targets(args.targets))
    except (OSError, ValueError) as exc:
        sys.stderr.write("error: targets vocabulary: %s\n" % exc)
        return 2
    baseline = dataset_stats(records, index)
    result = migrate_records(records, index)

    if result.id_collisions:
        sys.stderr.write(
            "error: %d new record_id(s) would be shared by more than one "
            "record, which would merge distinct observations. Nothing was "
            "written.\n" % len(result.id_collisions)
        )
        for new_id, olds in result.id_collisions:
            sys.stderr.write("  %s <- %s\n" % (new_id, ", ".join(olds)))
        return 2
    if result.stale_input_ids and not args.allow_stale_input_ids:
        sys.stderr.write(
            "error: %d input record(s) carry a record_id that the v0.1.0 hash "
            "of their own fields does not reproduce, so the file was edited "
            "after its ids were written and the old ids cannot be trusted as a "
            "key. Nothing was written. Pass --allow-stale-input-ids to migrate "
            "anyway (the map will then carry those ids as they are stored).\n"
            % len(result.stale_input_ids)
        )
        for record_id in result.stale_input_ids[:20]:
            sys.stderr.write("  %s\n" % record_id)
        return 2

    stats = dataset_stats(result.records, index)
    new_hard = _positive_deltas(
        _non_artefact(baseline["hard_findings"]),
        _non_artefact(stats["hard_findings"]),
    )
    if new_hard and not args.allow_new_findings:
        sys.stderr.write(
            "error: the migrated records would carry hard validation finding(s) "
            "the input did not have. Nothing was written.\n"
        )
        for code, delta in sorted(new_hard.items()):
            sys.stderr.write("  %s: +%d\n" % (code, delta))
        return 2

    report = render_report_md(
        stats,
        counts=result.counts,
        baseline=baseline,
        result=result,
        records_path=args.records,
        id_map_path=args.id_map,
        backup_path=None if args.no_backup else args.backup,
        index=index,
        generated=generated,
    )

    if not args.dry_run:
        if not args.no_backup:
            status = _make_backup(args.records, args.backup)
            if status is not None:
                sys.stderr.write("error: %s\n" % status)
                return 2
        _atomic_write_text(args.records, _records_text(result.records))
        _atomic_write_text(
            args.id_map,
            json.dumps(
                _id_map_document(result, records_path=args.records, generated=generated),
                ensure_ascii=False,
                indent=1,
            )
            + "\n",
        )
        _atomic_write_text(args.report, report)

    _print_run_summary(args, result, baseline, stats, index)
    return 1 if result.needs_a_human else 0


def _make_backup(records_path: str, backup_path: str) -> Optional[str]:
    """Preserve the pre-migration file. Returns an error message, or ``None``.

    An existing backup is never overwritten unless it is byte-identical to the
    file about to be replaced: it is the only copy of the v0.1.0 dataset, and
    overwriting it with a half-migrated file would make the migration
    irreversible.
    """
    if os.path.exists(backup_path):
        if _same_bytes(records_path, backup_path):
            return None
        return (
            "%s already exists and differs from %s. That backup is the only "
            "copy of the pre-migration dataset, so it is not overwritten. "
            "Nothing was written. Move it aside deliberately, or pass "
            "--no-backup if you have preserved it elsewhere."
            % (backup_path, records_path)
        )
    parent = os.path.dirname(os.path.abspath(backup_path)) or "."
    os.makedirs(parent, exist_ok=True)
    shutil.copy2(records_path, backup_path)
    return None


def _same_bytes(left: str, right: str) -> bool:
    try:
        with open(left, "rb") as a, open(right, "rb") as b:
            return a.read() == b.read()
    except OSError:  # pragma: no cover - defensive
        return False


def _print_run_summary(args, result, baseline, stats, index) -> None:
    counts = result.counts
    print("schema migration: curated records %s -> %s" % (args.from_version, TO_VERSION))
    print("  records:      %s" % args.records)
    suffix = " (not written: --dry-run)" if args.dry_run else ""
    print("  id map:       %s%s" % (args.id_map, suffix))
    print("  report:       %s%s" % (args.report, suffix))
    if args.no_backup:
        print("  backup:       none (--no-backup)")
    else:
        print("  backup:       %s%s" % (args.backup, suffix))
    print("  targets vocab: %d entries, %d matchable strings" % (len(index.ids), len(index)))
    print("")
    for key in (
        "records_in",
        "records_out",
        "ids_changed",
        "ids_unchanged",
        "perturbation_class_filled",
        "perturbation_class_null",
        "perturbation_polarity_filled",
        "perturbation_polarity_null",
        "target_eligible",
        "target_resolved",
        "target_unresolved_no_match",
        "target_ambiguous",
        "target_already_set",
        "target_not_eligible",
        "conflicts_with_references_in",
        "conflicts_with_references_remapped",
        "conflicts_with_references_unmapped",
        "records_already_at_target_version",
        "input_ids_not_reproducible_under_0_1_0",
    ):
        print("  %-48s %d" % (key, counts[key]))
    print("")
    print("  perturbation_class (counted):")
    for key, count in sorted(stats["perturbation_class"].items(), key=lambda kv: (-kv[1], kv[0])):
        print("    %-44s %d" % (key, count))
    print("  perturbation_polarity (counted):")
    for key, count in sorted(stats["perturbation_polarity"].items(), key=lambda kv: (-kv[1], kv[0])):
        print("    %-44s %d" % (key, count))
    print("  target_id (counted):")
    for key, count in sorted(stats["target_id"].items(), key=lambda kv: (-kv[1], kv[0])):
        print("    %-44s %d" % (key, count))
    print("")
    print(
        "  %-48s %d of %d"
        % ("records_hard_clean", stats["records_hard_clean"], stats["records"])
    )
    if stats["hard_findings"]:
        print("  hard findings remaining (counted):")
        for code, count in sorted(stats["hard_findings"].items()):
            print("    %-44s %d" % (code, count))
    print("  warning findings (before -> after, counted):")
    codes = sorted(set(baseline["warning_findings"]) | set(stats["warning_findings"]))
    for code in codes:
        print(
            "    %-44s %d -> %d"
            % (code, baseline["warning_findings"].get(code, 0), stats["warning_findings"].get(code, 0))
        )
    print("")
    print("  condition fields the migration may not fill (left null):")
    for name in CONDITION_FIELDS_THE_MIGRATION_MAY_NOT_FILL:
        slot = stats["condition_fields"][name]
        print("    %-44s set=%d null=%d" % (name, slot["set"], slot["null"]))
    print("")
    print(
        "  %d genes / %d studies / %d labs; %d gene(s) in exactly one study, "
        "%d in exactly one lab."
        % (
            stats["genes_total"],
            stats["gene_studies_total"],
            stats["gene_labs_total"],
            stats["genes_in_exactly_one_study"],
            stats["genes_in_exactly_one_lab"],
        )
    )
    if result.ambiguous:
        print("")
        print("  ambiguous target matches, left null for a human (counted):")
        for record_id, candidates, _text in result.ambiguous:
            print("    %s  %s" % (record_id, ", ".join(candidates)))
    if result.unmapped_conflicts:
        print("")
        print("  conflicts_with references not in the id map (left as they are):")
        for record_id, reference in result.unmapped_conflicts:
            print("    %s -> %s" % (record_id, reference))


def cmd_report(args) -> int:
    """Recompute the report from the committed dataset (and the kept input).

    Two modes, and the command says which one it is in:

    *complete* -- the preserved pre-migration file is present, so the whole
    migration is re-run in memory from it. Every number in the published
    report is then reproduced, including the before/after validation tables,
    and the re-migration's output is compared record by record against the
    committed dataset. That is as close to an audit as this project can get.

    *dataset only* -- the preserved input is gone (it is not version
    controlled). Everything the committed ``records.jsonl`` and the committed
    id map determine is still recomputed; the handful of rows that only the
    migration event can answer are omitted rather than invented, so the render
    remains a line-subsequence of the published report and is compared as one.
    """
    generated = _now()
    if not os.path.isfile(args.records):
        sys.stderr.write("error: no such file: %s\n" % args.records)
        return 2
    records, errors, versions = read_records(args.records, TO_VERSION)
    if errors:
        for number, message in errors:
            sys.stderr.write("unreadable line %d: %s\n" % (number, message))
        return 2
    if set(versions) != {TO_VERSION}:
        sys.stderr.write(
            "error: report is for schema_version %r; %s declares %s. Run the "
            "migration first.\n"
            % (TO_VERSION, args.records, ", ".join(sorted(versions)))
        )
        return 2

    try:
        index = build_synonym_index(load_targets(args.targets))
    except (OSError, ValueError) as exc:
        sys.stderr.write("error: targets vocabulary: %s\n" % exc)
        return 2
    stats = dataset_stats(records, index)

    counts: Optional[Dict[str, int]] = None
    baseline: Optional[Dict[str, object]] = None
    result: Optional[MigrationResult] = None
    event_complete = False
    reproduces: Optional[bool] = None
    backup_path: Optional[str] = None

    if args.backup and os.path.isfile(args.backup):
        kept, kept_errors, kept_versions = read_records(args.backup, FROM_VERSION)
        if not kept_errors and set(kept_versions) == {args.from_version}:
            baseline = dataset_stats(kept, index)
            result = migrate_records(kept, index)
            counts = result.counts
            event_complete = True
            backup_path = args.backup
            reproduces = [r.to_dict() for r in result.records] == [
                r.to_dict() for r in records
            ]

    if counts is None:
        # No usable preserved input: fall back to the committed id map, which
        # is version controlled and answers the identity questions on its own.
        if os.path.isfile(args.id_map):
            with open(args.id_map, "r", encoding="utf-8") as handle:
                document = json.load(handle)
            mapping = document.get("map") or {}
            counts = {
                "records_in": len(mapping),
                "records_out": len(records),
                "ids_changed": sum(1 for old, new in mapping.items() if old != new),
                "ids_unchanged": sum(1 for old, new in mapping.items() if old == new),
                "conflicts_with_references_in": stats["conflicts_with_references"],
                "conflicts_with_references_remapped": stats[
                    "conflicts_with_references"
                ],
                "conflicts_with_references_unmapped": 0,
            }
        result = MigrationResult(ambiguous=list(stats["target_recheck"]["ambiguous"]))

    report = render_report_md(
        stats,
        counts=counts,
        baseline=baseline,
        result=result,
        records_path=args.records,
        id_map_path=args.id_map,
        backup_path=backup_path,
        index=index,
        generated=generated,
        recomputed_only=True,
        event_complete=event_complete,
    )

    if args.out and not args.dry_run:
        _atomic_write_text(args.out, report)
    if args.print_report:
        print(report)
    else:
        print("recomputed from %s" % args.records)
        print(
            "  mode                     %s"
            % (
                "complete (re-ran the migration from %s)" % args.backup
                if event_complete
                else "dataset only (no preserved pre-migration file)"
            )
        )
        print("  records                  %d" % stats["records"])
        print("  record_ids unique        %d" % stats["record_ids_unique"])
        print("  record_ids consistent    %d" % stats["record_ids_consistent"])
        print("  records_hard_clean       %d" % stats["records_hard_clean"])
        print(
            "  targets resolved         %d of %d eligible"
            % (stats["target_resolved_records"], stats["target_eligible_records"])
        )
        print(
            "  genes / studies / labs   %d / %d / %d"
            % (
                stats["genes_total"],
                stats["gene_studies_total"],
                stats["gene_labs_total"],
            )
        )
        print(
            "  condition slots unfilled %d" % stats["condition_field_slots_unfilled"]
        )
        if reproduces is not None:
            print(
                "  re-migrating the preserved input reproduces the committed "
                "dataset: %s" % ("yes" if reproduces else "NO")
            )
        if args.out:
            print("  report written to        %s" % args.out)

    if reproduces is False:
        print(
            "  the committed dataset is not what this code produces from the "
            "preserved input; it was edited afterwards."
        )
        return 1

    if args.compare and os.path.isfile(args.compare):
        with open(args.compare, "r", encoding="utf-8") as handle:
            published = _report_body(handle.read())
        recomputed = _report_body(report)
        if event_complete:
            agrees = published == recomputed
            how = "the published report matches the recomputation"
        else:
            agrees = _is_subsequence(recomputed.splitlines(), published.splitlines())
            how = (
                "every line the recomputation could derive matches the "
                "published report, in order"
            )
        if agrees:
            print("  %s." % how)
        else:
            print("  published report DIFFERS from the recomputation:")
            for line in list(
                difflib.unified_diff(
                    published.splitlines(),
                    recomputed.splitlines(),
                    fromfile=args.compare,
                    tofile="recomputed",
                    lineterm="",
                )
            )[:60]:
                print("    %s" % line)
            return 1
    return 0


def _is_subsequence(needle: Sequence[str], haystack: Sequence[str]) -> bool:
    """True if every line of ``needle`` occurs in ``haystack``, in order."""
    iterator = iter(haystack)
    return all(any(line == other for other in iterator) for line in needle)


def _report_body(text: str) -> str:
    """The report minus its generation timestamp line, for comparison."""
    return "\n".join(
        line
        for line in text.splitlines()
        if not line.startswith("Generated ") and not line.startswith("Recomputed ")
    )


# --------------------------------------------------------------------------- #
# argument parsing
# --------------------------------------------------------------------------- #


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python3 -m syn4rch.curate.migrate",
        description=(
            "Migrate data/curated/records.jsonl from schema v0.1.0 to v0.2.0: "
            "derive the perturbation axis, resolve the molecular target out of "
            "the curated targeting text, recompute every record_id, and leave "
            "the condition dimension null because filling it would require "
            "reading the papers."
        ),
    )
    sub = parser.add_subparsers(dest="command")

    run = sub.add_parser(
        "run",
        help="migrate the curated file in place (atomically, with a backup)",
        description=(
            "Refuses to start unless every input record declares the version "
            "it migrates from, and refuses to write if the result would carry "
            "a hard validation finding the input did not have."
        ),
    )
    run.add_argument("--records", default=None, metavar="PATH",
                     help="curated records (default: data/curated/records.jsonl)")
    run.add_argument("--id-map", default=None, metavar="PATH",
                     help="old-id -> new-id map (default: "
                          "data/curated/record_id_map_0.1.0_to_0.2.0.json)")
    run.add_argument("--report", default=None, metavar="PATH",
                     help="migration report (default: reports/schema_migration_0.2.0.md)")
    run.add_argument("--backup", default=None, metavar="PATH",
                     help="pre-migration copy (default: "
                          "data/curated/records.0.1.0.jsonl.bak)")
    run.add_argument("--no-backup", action="store_true",
                     help="do not keep a pre-migration copy. Only for a file you "
                          "have preserved elsewhere; the migration is otherwise "
                          "irreversible")
    run.add_argument("--from-version", default=FROM_VERSION, metavar="V",
                     help="the version to migrate from; an input declaring "
                          "anything else is refused (default: %s)" % FROM_VERSION)
    run.add_argument("--allow-stale-input-ids", action="store_true",
                     help="migrate even though some input record_id is not the "
                          "one v0.1.0 would have computed from its own fields")
    run.add_argument("--allow-new-findings", action="store_true",
                     help="write even though the migrated records carry a hard "
                          "validation finding the input did not have")
    run.add_argument("--targets", default=None, metavar="PATH",
                     help="targets vocabulary to resolve against (default: "
                          "data/vocab/targets.json via syn4rch.vocab)")
    run.add_argument("--dry-run", action="store_true",
                     help="count, validate and print; write nothing")
    run.set_defaults(func=cmd_run)

    report = sub.add_parser(
        "report",
        help="recompute every number in the migration report from the dataset",
        description=(
            "Reads the curated file at %s and re-derives the report. "
            "docs/CURATION_PROTOCOL.md 6: a number that cannot be recomputed "
            "from the committed dataset does not belong in a report." % TO_VERSION
        ),
    )
    report.add_argument("--records", default=None, metavar="PATH")
    report.add_argument("--id-map", default=None, metavar="PATH")
    report.add_argument("--out", default=None, metavar="PATH",
                        help="rewrite the report here (default: write nothing)")
    report.add_argument("--compare", default=None, metavar="PATH",
                        help="compare the recomputation against a published "
                             "report and exit non-zero if they differ "
                             "(default: reports/schema_migration_0.2.0.md)")
    report.add_argument("--targets", default=None, metavar="PATH")
    report.add_argument("--backup", default=None, metavar="PATH",
                        help="the preserved pre-migration file; when present "
                             "the whole migration is re-run from it and the "
                             "before/after tables are recomputed too "
                             "(default: data/curated/records.%s.jsonl.bak)"
                             % FROM_VERSION)
    report.add_argument("--from-version", default=FROM_VERSION, metavar="V")
    report.add_argument("--print-report", action="store_true",
                        help="print the whole recomputed report to stdout")
    report.add_argument("--dry-run", action="store_true")
    report.set_defaults(func=cmd_report)
    return parser


def _apply_defaults(args) -> None:
    if args.records is None:
        args.records = _default("data", "curated", "records.jsonl")
    if args.id_map is None:
        args.id_map = _default(
            "data", "curated", "record_id_map_%s_to_%s.json" % (FROM_VERSION, TO_VERSION)
        )
    if args.command == "run":
        if args.report is None:
            args.report = _default("reports", "schema_migration_%s.md" % TO_VERSION)
        if args.backup is None:
            args.backup = _default(
                "data", "curated", "records.%s.jsonl.bak" % FROM_VERSION
            )
    elif args.command == "report":
        if args.compare is None:
            args.compare = _default("reports", "schema_migration_%s.md" % TO_VERSION)
        if args.backup is None:
            args.backup = _default(
                "data", "curated", "records.%s.jsonl.bak" % FROM_VERSION
            )


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        return 2
    _apply_defaults(args)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
