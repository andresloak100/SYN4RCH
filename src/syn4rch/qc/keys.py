"""Grouping keys shared by the QC checks.

Keeping them in one place means "same observation", "same lab" and "same
claim" mean exactly one thing across duplicate, conflict and
non-independence detection.

Python 3.10 compatible, standard library only.
"""

from typing import Optional, Tuple

from ..schema.record import normalize_for_hash

__all__ = [
    "norm",
    "population_key",
    "claim_key",
    "lab_key",
    "observation_key",
    "study_key",
]


def norm(value) -> str:
    """Case-folded, whitespace-collapsed, NFKC form. ``None`` -> ``""``."""
    return normalize_for_hash(value).casefold()


def study_key(record) -> str:
    return norm(record.study_id)


def population_key(record) -> str:
    """Identity of the population for grouping.

    Prefers the resolved ``population_id``. When it is ``None`` (SCHEMA.md 3.2
    keeps such records rather than force-mapping them) the paper's own wording
    is used, prefixed ``raw:`` so an unresolved population can never silently
    group with a resolved one. See docs/DECISIONS.md (D-010).
    """
    if record.population_id:
        return norm(record.population_id)
    return "raw:" + norm(record.population_raw)


def claim_key(record) -> Tuple[str, str, str]:
    """``(population, valence, outcome_id)`` -- the claim under test.

    Two records with the same claim key make a statement about the same thing
    and are therefore comparable for conflict.
    """
    return (
        population_key(record),
        norm(record.valence),
        norm(record.outcome_id),
    )


def lab_key(record) -> Optional[Tuple[str, str]]:
    """``(last_author, affiliation)`` or ``None`` if lab identity is unknown.

    SCHEMA.md 5 names lab identity as the largest leakage risk in this corpus.
    ``None`` means "cannot be grouped": such records are skipped by
    non-independence detection rather than lumped into one anonymous lab.
    """
    author = norm(record.last_author)
    if not author:
        return None
    return (author, norm(record.last_author_affiliation))


def observation_key(record) -> Tuple[str, str, str, str]:
    """``(study, population, valence, outcome_id)`` plus locator is used by
    duplicate detection; this is its study-scoped part."""
    return (study_key(record),) + claim_key(record)
