"""Grouping keys shared by the QC checks.

Keeping them in one place means "same observation", "same lab", "same
perturbation" and "same claim" mean exactly one thing across duplicate,
conflict and non-independence detection.

What changed in schema v0.2.0
-----------------------------
:func:`claim_key` used to be ``(population, valence, outcome_id)``. Measured
against the first 362 curated records, that key put 18 cells into direct
``increase`` vs ``decrease`` opposition. The v0.2.0 key leaves 2: **16 of the 18
apparent contradictions in this dataset were never contradictions.** They were a
male cohort compared against a female one, a chow-fed cohort against a
high-fat-fed one, or one gene's knockout against another gene's -- experiments
that were never comparable, recorded as laboratories disagreeing.

The cost is in the same measurement and is larger than the benefit looks.
Cells reported by more than one study fall from 39 to 11. Of the 26 lost
(measured over usable records by :mod:`syn4rch.feasibility.replication`), 18 go
to a genuine perturbation difference and 1 to two differently stated
conditions -- both correct separations of work that was never comparable. The
remaining **7 are lost only because one paper states a condition the other
leaves unknown**: ``unknown`` cannot match a stated value, so a paper reporting
male mice stops replicating a paper that reports no sex, even if both ran
males. Nothing in this corpus can adjudicate that, and assuming they matched is
what SCHEMA.md 2 rule 4 forbids. That 7 is reported on its own line and never
as heterogeneity resolved.

``valence`` has left the claim key, replaced by :func:`perturbation_key`, which
carries the perturbation's class, its signed polarity and its molecular target
together. See docs/DECISIONS.md (D-023 .. D-026).

Python 3.10 compatible, standard library only.
"""

from typing import Optional, Tuple

from ..schema.enums import Challenge, Sex, Diet, AgeGroup
from ..schema.record import normalize_for_hash
from ..schema.rules import (
    perturbation_class_from_method,
    perturbation_polarity_from,
)

__all__ = [
    "UNKNOWN",
    "NULL",
    "CONDITION_FIELDS",
    "CLAIM_KEY_FIELDS",
    "POOLED_CLAIM_KEY_FIELDS",
    "norm",
    "population_key",
    "perturbation_key",
    "condition_key",
    "condition_specificity",
    "claim_key",
    "claim_key_pooled",
    "lab_key",
    "observation_key",
    "study_key",
]

#: Token for a condition that is not known, whether because the curator did
#: not address it (``None``) or because the paper was checked and does not say
#: (``not_reported``). The two are different statements about the *curation*
#: and both are preserved in the record, but they are the same statement about
#: the *experiment*: the condition could be anything. Grouping them apart
#: would mean two records that each fail to state a diet are held to be
#: incomparable, which would make a conflict undetectable in the 197 records of
#: 362 that carry no diet -- and would equally make replication undetectable
#: there. They are therefore pooled, and :func:`condition_specificity` is how a
#: report says that the pooling happened. See docs/DECISIONS.md (D-023).
UNKNOWN = "unknown"

#: Token for an absent non-condition component, kept visible so that a key
#: printed in a report cannot be mistaken for a populated one.
NULL = "<null>"

#: The condition dimension, in the order :func:`condition_key` returns it.
CONDITION_FIELDS: Tuple[str, ...] = ("sex", "diet", "age_group", "challenge")

#: Component names of :func:`claim_key`, in order, for labelling reports.
CLAIM_KEY_FIELDS: Tuple[str, ...] = (
    "population",
    "perturbation_class",
    "perturbation_polarity",
    "perturbation_target",
    "outcome_id",
) + CONDITION_FIELDS

#: Component names of :func:`claim_key_pooled` -- the schema v0.1.0 key.
POOLED_CLAIM_KEY_FIELDS: Tuple[str, ...] = ("population", "valence", "outcome_id")

#: Condition enum members that mean "the paper does not say". Pooled with
#: ``None`` by :func:`condition_key`; see :data:`UNKNOWN`.
_NOT_REPORTED = frozenset(
    {
        Sex.not_reported.value,
        Diet.not_reported.value,
        AgeGroup.not_reported.value,
        Challenge.not_reported.value,
    }
)


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


def _target_component(record) -> str:
    """Identity of the molecular target, or :data:`NULL`.

    Same resolved-then-raw pattern as :func:`population_key`, and the same
    reason for the ``raw:`` prefix: an unresolved gene symbol must not pool
    with a resolved one.
    """
    target_id = getattr(record, "target_id", None)
    if target_id:
        return norm(target_id)
    target_raw = getattr(record, "target_raw", None)
    if target_raw:
        return "raw:" + norm(target_raw)
    return NULL


def perturbation_key(record) -> Tuple[str, str, str]:
    """``(class, polarity, target)`` -- what was done to the population.

    The three travel together by construction, because a polarity means
    different physical things in different classes: ``decrease`` is less firing
    under ``circuit_activity`` and less of a gene product under ``molecular``.
    Nothing in this project may key on polarity alone, and this function is why
    it is awkward to try. See :class:`syn4rch.schema.enums.PerturbationPolarity`.

    ``class`` and ``polarity`` prefer the stored fields and fall back to
    :mod:`syn4rch.schema.rules`, so a record written before schema v0.2.0 keys
    identically to the same record after migration.
    """
    cls = getattr(record, "perturbation_class", None)
    if cls is None:
        cls = perturbation_class_from_method(record.method)
    pol = getattr(record, "perturbation_polarity", None)
    if pol is None:
        pol = perturbation_polarity_from(record.method, record.valence)
    return (
        norm(cls) or NULL,
        norm(pol) or NULL,
        _target_component(record),
    )


def condition_key(record) -> Tuple[str, str, str, str]:
    """``(sex, diet, age_group, challenge)`` -- the condition the work ran under.

    Any slot that is not known reads :data:`UNKNOWN`. Two records whose
    condition keys differ are **not comparable**: they may disagree because the
    biology disagrees, or because one was a fasted male on a high-fat diet and
    the other a fed female on chow. The dataset cannot tell those apart, so
    this key refuses to.
    """
    out = []
    for name in CONDITION_FIELDS:
        value = norm(getattr(record, name, None))
        out.append(UNKNOWN if (not value or value in _NOT_REPORTED) else value)
    return tuple(out)  # type: ignore[return-value]


def condition_specificity(record) -> int:
    """How many of the four condition slots are actually known (0..4).

    A conflict found between two records of specificity 0 is a different kind
    of statement from one found between two records of specificity 4, and a
    report that does not distinguish them is overclaiming. ``0`` means the
    comparison matched on "both unknown", not on "both the same".
    """
    return sum(1 for slot in condition_key(record) if slot != UNKNOWN)


def claim_key(record):
    """The claim under test: population, perturbation, outcome, condition.

    Two records with the same claim key make a statement about the same thing
    under the same conditions and are therefore comparable for conflict and for
    replication. Components are named by :data:`CLAIM_KEY_FIELDS`.

    ``valence`` is deliberately absent. For a circuit manipulation it is
    restated by ``perturbation_polarity``; for a molecular one it is honestly
    ``unclear`` and carries no information. Keeping it would have fragmented
    the key without separating anything.
    """
    return (
        (population_key(record),)
        + perturbation_key(record)
        + (norm(record.outcome_id),)
        + condition_key(record)
    )


def claim_key_pooled(record) -> Tuple[str, str, str]:
    """``(population, valence, outcome_id)`` -- the schema v0.1.0 claim key.

    Kept, and kept exact, for one purpose: the difference between a count under
    this key and the same count under :func:`claim_key` *is* the measurement of
    how much apparent disagreement in the corpus was condition heterogeneity.
    Deleting it would have made the v0.2.0 change unfalsifiable.

    It must not be used to decide whether two records conflict. It pools
    experiments that were never comparable, which is the defect v0.2.0 exists
    to fix.
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


def observation_key(record):
    """``(study,) + claim_key`` -- duplicate detection's study-scoped part.

    Duplicate detection adds the locator. Note that under v0.2.0 this is
    *stricter* than it was: two records in one study at one locator reporting
    one readout for a male cohort and a female cohort are two observations, not
    a duplicated one, and are no longer flagged as a duplicate.
    """
    return (study_key(record),) + claim_key(record)
