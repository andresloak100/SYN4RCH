"""Conflict detection between studies.

A conflict is two or more *different studies* making the same claim --
same (population, valence, outcome_id) -- with incompatible ``direction``.

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
from .keys import claim_key, study_key

__all__ = [
    "OPPOSED_PAIRS",
    "directions_conflict",
    "find_conflicts",
    "annotate_conflicts",
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
    by_claim: Dict[Tuple[str, str, str], List] = {}
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
                        "Claim (population=%s, valence=%s, outcome=%s) is "
                        "reported as %s by study %s (%s) and as %s by study %s "
                        "(%s). Recorded on both sides; not resolved. Conflicts "
                        "are data."
                        % (
                            claim[0],
                            claim[1],
                            claim[2],
                            a.direction.value if a.direction else "null",
                            a.study_id,
                            a.record_id,
                            b.direction.value if b.direction else "null",
                            b.study_id,
                            b.record_id,
                        )
                    ),
                    record_ids=sorted([a.record_id, b.record_id]),
                    field_name="direction",
                    details={
                        "claim": list(claim),
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
