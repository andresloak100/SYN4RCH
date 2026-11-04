"""Conflict detection between studies.

A conflict is two or more *different studies* making the same claim -- same
:func:`syn4rch.qc.keys.claim_key` -- with incompatible ``direction``.

What schema v0.2.0 changed here
-------------------------------
The claim key used to be ``(population, valence, outcome_id)``. It now carries
the perturbation's class, polarity and molecular target, and the condition the
experiment ran under. Measured on the first 362 curated records, that takes the
number of directly opposed cells from 18 to 2. Sixteen of the eighteen
"laboratories disagree" findings this module used to emit were a male cohort
compared against a female one, a chow-fed cohort against a high-fat-fed one, or
a gene knockout against a gene overexpression -- comparisons that were never
valid. They were not conflicts, and reporting them as conflicts was the single
largest distortion in the dataset.

:func:`key_comparison` keeps the old count computable, because the difference
between the two numbers is the only evidence that the change did anything.

A conflict is also now **qualified by how much of the condition is known**. Two
records that each fail to state sex, diet, age and challenge match each other on
"both unknown", not on "both the same", and a disagreement between them may
still be condition heterogeneity that neither paper reported. Every finding
carries ``condition_specificity`` (0..4) so that a report cannot print those two
situations as one number. See docs/DECISIONS.md (D-023, D-027).

**Conflicts are data.** The literature disagrees about hypothalamic
perturbation outcomes, and that disagreement is one of the most informative
signals in this corpus: a claim that replicates across labs and a claim that
flips between labs are not the same evidence. Therefore this module:

* records the conflict on **both** sides, in ``conflicts_with`` (SCHEMA.md 3.6
  says this field is populated by QC, not by the extractor);
* **never** drops, prefers, down-ranks, reconciles or averages either side;
* **never** rewrites ``direction``;
* reports at severity ``INFO``, because the dataset is correct and the
  literature is in tension.

Conflicts *within* one study are not conflicts here: they are either
different conditions (and the claim key is too coarse) or an extraction error,
and both are handled by duplicate/locator checks.

Python 3.10 compatible, standard library only.
"""

from typing import Dict, List, Optional, Tuple

from ..schema.enums import Direction
from .findings import QCFinding, QCSeverity
from .keys import (
    CLAIM_KEY_FIELDS,
    claim_key,
    claim_key_pooled,
    condition_specificity,
    study_key,
)

__all__ = [
    "OPPOSED_PAIRS",
    "directions_conflict",
    "describe_claim",
    "find_conflicts",
    "annotate_conflicts",
    "key_comparison",
]

#: Directly opposed directions.
OPPOSED_PAIRS = frozenset(
    {
        frozenset({Direction.increase, Direction.decrease}),
    }
)

#: An asserted effect against an asserted null for the same claim.
_EFFECT_DIRECTIONS = frozenset({Direction.increase, Direction.decrease})


def directions_conflict(a, b) -> Optional[str]:
    """Return a conflict code for directions ``a`` vs ``b``, else ``None``.

    Two rules only, both narrow on purpose (docs/DECISIONS.md D-011):

    * ``increase`` vs ``decrease`` -> ``QC_CONFLICT_OPPOSED``
    * ``increase``/``decrease`` vs ``no_change`` -> ``QC_CONFLICT_EFFECT_VS_NULL``

    ``not_reported`` conflicts with nothing (there is no claim to contradict).
    ``non_monotonic`` is not treated as conflicting with a monotonic direction:
    a non-monotonic response can look like either arm depending on the dose or
    time window sampled, so calling that a contradiction would be an
    interpretation, not an observation.
    """
    if a is None or b is None:
        return None
    if a is Direction.not_reported or b is Direction.not_reported:
        return None
    if frozenset({a, b}) in OPPOSED_PAIRS:
        return "QC_CONFLICT_OPPOSED"
    if (a in _EFFECT_DIRECTIONS and b is Direction.no_change) or (
        b in _EFFECT_DIRECTIONS and a is Direction.no_change
    ):
        return "QC_CONFLICT_EFFECT_VS_NULL"
    return None


def describe_claim(claim) -> str:
    """Render a claim key as ``field=value`` pairs, labelled by name.

    The key has nine components since v0.2.0. Printing it as a bare tuple in a
    finding would leave a reader guessing which slot was the diet, so the names
    travel with the values.
    """
    return ", ".join(
        "%s=%s" % (name, value)
        for name, value in zip(CLAIM_KEY_FIELDS, claim)
    )


def _condition_caveat(a, b) -> str:
    """One sentence on what the condition match actually established.

    When neither record states any condition, the two were matched on "both
    unknown". That is a weaker basis for calling something a contradiction than
    a matched sex and diet, and the finding says so rather than leaving the
    reader to assume the stronger reading.
    """
    spec = min(condition_specificity(a), condition_specificity(b))
    if spec == 0:
        return (
            "Neither record states sex, diet, age group or challenge, so the "
            "two were matched on 'condition unknown' rather than on matching "
            "conditions; this disagreement may still be condition "
            "heterogeneity that neither paper reported."
        )
    if spec < len(CLAIM_KEY_FIELDS) - 5:
        return (
            "Conditions match on %d of 4 slots; the rest are unknown in at "
            "least one of the two records." % spec
        )
    return "Conditions match on all four slots."


def _pairs_in_claim(group):
    """Yield ``(code, record_a, record_b)`` for conflicting cross-study pairs."""
    for i in range(len(group)):
        for j in range(i + 1, len(group)):
            a, b = group[i], group[j]
            if study_key(a) == study_key(b):
                continue
            code = directions_conflict(a.direction, b.direction)
            if code is not None:
                yield code, a, b


def find_conflicts(records) -> List[QCFinding]:
    """Report every cross-study direction conflict. Records are not modified."""
    by_claim: Dict[Tuple[str, ...], List] = {}
    for record in records:
        if not record.outcome_id:
            continue  # missing outcome_id is a HARD schema finding
        by_claim.setdefault(claim_key(record), []).append(record)

    out: List[QCFinding] = []
    for claim, group in sorted(by_claim.items()):
        for code, a, b in _pairs_in_claim(group):
            out.append(
                QCFinding(
                    code=code,
                    severity=QCSeverity.INFO,
                    message=(
                        "Claim (%s) is reported as %s by study %s (%s) and as "
                        "%s by study %s (%s). %s Recorded on both sides; not "
                        "resolved. Conflicts are data."
                        % (
                            describe_claim(claim),
                            a.direction.value if a.direction else "null",
                            a.study_id,
                            a.record_id,
                            b.direction.value if b.direction else "null",
                            b.study_id,
                            b.record_id,
                            _condition_caveat(a, b),
                        )
                    ),
                    record_ids=sorted([a.record_id, b.record_id]),
                    field_name="direction",
                    details={
                        "claim": list(claim),
                        "claim_fields": list(CLAIM_KEY_FIELDS),
                        "condition_specificity": min(
                            condition_specificity(a), condition_specificity(b)
                        ),
                        "sides": [
                            {
                                "record_id": a.record_id,
                                "study_id": a.study_id,
                                "direction": a.direction.value
                                if a.direction
                                else None,
                                "support_quote": a.support_quote,
                                "source_locator": a.source_locator,
                            },
                            {
                                "record_id": b.record_id,
                                "study_id": b.study_id,
                                "direction": b.direction.value
                                if b.direction
                                else None,
                                "support_quote": b.support_quote,
                                "source_locator": b.source_locator,
                            },
                        ],
                    },
                )
            )
    return out


def annotate_conflicts(records) -> List[QCFinding]:
    """Populate ``conflicts_with`` on both sides of every conflict, in place.

    Idempotent: ids are added to a sorted set, never duplicated. Existing ids
    are preserved -- a human may have recorded a conflict this narrow rule
    cannot see, and QC does not overwrite that. Returns the findings.
    """
    findings = find_conflicts(records)
    by_id: Dict[str, List] = {}
    for record in records:
        by_id.setdefault(record.record_id, []).append(record)

    for finding in findings:
        left, right = finding.record_ids[0], finding.record_ids[1]
        for this_id, other_id in ((left, right), (right, left)):
            for record in by_id.get(this_id, []):
                if other_id != record.record_id and other_id not in record.conflicts_with:
                    record.conflicts_with = sorted(
                        set(record.conflicts_with) | {other_id}
                    )
    return findings


def key_comparison(records) -> dict:
    """Conflict counts under the v0.2.0 claim key and under the v0.1.0 one.

    The v0.1.0 key (:func:`syn4rch.qc.keys.claim_key_pooled`) pools experiments
    that were never comparable. It is computed here for exactly one reason: the
    gap between the two counts is the measurement of how much of this corpus's
    apparent disagreement was condition and perturbation heterogeneity rather
    than laboratories contradicting each other. Without the comparison the
    v0.2.0 change would be an assertion instead of a result.

    The ``pooled`` numbers must never be used to decide that two records
    conflict. They are a baseline to be reported against, nothing more.
    """
    out = {}
    for label, keyfn in (("strict", claim_key), ("pooled", claim_key_pooled)):
        groups: Dict[Tuple[str, ...], List] = {}
        for record in records:
            if not record.outcome_id:
                continue
            groups.setdefault(tuple(keyfn(record)), []).append(record)
        codes: Dict[str, int] = {}
        cells_with_conflict = 0
        for group in groups.values():
            found = {code for code, _, _ in _pairs_in_claim(group)}
            for code in found:
                codes[code] = codes.get(code, 0) + 1
            if found:
                cells_with_conflict += 1
        out[label] = {
            "n_cells": len(groups),
            "n_cells_with_conflict": cells_with_conflict,
            "cells_by_code": dict(sorted(codes.items())),
        }
    out["key_fields"] = {
        "strict": list(CLAIM_KEY_FIELDS),
        "pooled": ["population", "valence", "outcome_id"],
    }
    out["note"] = (
        "A cell counted under 'pooled' but not under 'strict' is not a "
        "resolved disagreement. It is a comparison that should never have "
        "been made."
    )
    return out
