"""SYN4RCH controlled vocabularies.

Public interface (the rest of the project codes against exactly this)::

    from syn4rch.vocab import load, resolve, report_unmatched

    load(kind)                  -> dict mapping id -> entry dict
    resolve(kind, raw)          -> str | None   (exact + curated-synonym only)
    report_unmatched(kind, xs)  -> list[str]    (distinct strings that did not resolve)

`resolve` performs NO fuzzy matching of any kind: no edit distance, no substring
containment, no prefix scoring. A string resolves only if, after normalisation, it
equals an entry id, an entry label, or one of that entry's curated synonyms.
Returning None is normal and expected; it is how unknown populations stay honest
(see docs/SCHEMA.md 3.2: population_id is null when there is no confident match).

Python 3.10 compatible, standard library only.
"""

from .core import (  # noqa: F401
    KINDS,
    VocabError,
    VocabCollisionError,
    load,
    load_all,
    resolve,
    resolve_all,
    report_unmatched,
    normalise,
    selfcheck,
    vocab_dir,
)

__all__ = [
    "KINDS",
    "VocabError",
    "VocabCollisionError",
    "load",
    "load_all",
    "resolve",
    "resolve_all",
    "report_unmatched",
    "normalise",
    "selfcheck",
    "vocab_dir",
]
