"""Duplicate, collision and non-independence detection.

Three distinct things, deliberately not conflated:

1. **Exact duplicates** -- byte-identical records (same ``record_id``, same
   content). A harvesting or merge bug. WARNING.
2. **Locator collisions** -- same study, same ``source_locator``, same claim,
   but *different* content. Either one figure panel was extracted twice with
   divergent values, or the locator is too coarse to identify the
   observation. WARNING; needs a human.
3. **Possible non-independence** -- different studies from the *same lab*
   reporting the same (population, valence, outcome, direction). Not an
   error: it is normal and expected that a group follows up its own work.
   Reported as INFO so that leakage control (SCHEMA.md 5) can group these
   records, and never as something to delete.

Nothing here removes, merges or edits a record.

Python 3.10 compatible, standard library only.
"""

from typing import Dict, List, Tuple

from .findings import QCFinding, QCSeverity
from .keys import claim_key, lab_key, norm, study_key

__all__ = [
    "find_exact_duplicates",
    "find_locator_collisions",
    "find_possible_non_independence",
    "run_duplicate_checks",
]


def _content_signature(record) -> str:
    """Record content with ``conflicts_with`` excluded.

    ``conflicts_with`` is populated by QC itself (SCHEMA.md 3.6), so including
    it would make duplicate detection depend on whether conflict annotation
    had already run.
    """
    import json

    data = record.to_dict()
    data.pop("conflicts_with", None)
    return json.dumps(data, sort_keys=True, ensure_ascii=False)


def find_exact_duplicates(records) -> List[QCFinding]:
    """Records sharing a ``record_id``, split by whether content agrees."""
    by_id: Dict[str, List] = {}
    for record in records:
        by_id.setdefault(record.record_id, []).append(record)

    out: List[QCFinding] = []
    for record_id, group in sorted(by_id.items(), key=lambda kv: str(kv[0])):
        if len(group) < 2:
            continue
        signatures = {_content_signature(r) for r in group}
        ids = [r.record_id for r in group]
        if len(signatures) == 1:
            out.append(
                QCFinding(
                    code="QC_DUP_EXACT",
                    severity=QCSeverity.WARNING,
                    message=(
                        "%d byte-identical copies of record %s; keep one."
                        % (len(group), record_id)
                    ),
                    record_ids=ids,
                    field_name="record_id",
                    details={"copies": len(group)},
                )
            )
        else:
            out.append(
                QCFinding(
                    code="QC_DUP_ID_CONTENT_DIVERGENT",
                    severity=QCSeverity.WARNING,
                    message=(
                        "%d records share record_id %s but their content "
                        "differs; the identifying fields are identical while "
                        "the extracted values are not. Needs a human."
                        % (len(group), record_id)
                    ),
                    record_ids=ids,
                    field_name="record_id",
                    details={"variants": len(signatures)},
                )
            )
    return out


def find_locator_collisions(records) -> List[QCFinding]:
    """Same study + same ``source_locator`` + same claim, different content."""
    by_slot: Dict[Tuple[str, str, Tuple[str, str, str]], List] = {}
    for record in records:
        locator = norm(record.source_locator)
        if not locator:
            continue  # absence is a HARD schema finding, not a collision
        by_slot.setdefault(
            (study_key(record), locator, claim_key(record)), []
        ).append(record)

    out: List[QCFinding] = []
    for (study, locator, claim), group in sorted(
        by_slot.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])
    ):
        if len(group) < 2:
            continue
        if len({r.record_id for r in group}) == 1:
            continue  # already reported by find_exact_duplicates
        out.append(
            QCFinding(
                code="QC_LOCATOR_COLLISION",
                severity=QCSeverity.WARNING,
                message=(
                    "%d distinct records in study %s cite the same locator "
                    "%r for the same claim %s. Either the panel was extracted "
                    "twice or the locator is too coarse."
                    % (len(group), study, locator, "/".join(claim))
                ),
                record_ids=[r.record_id for r in group],
                field_name="source_locator",
                details={"study_id": study, "source_locator": locator},
            )
        )
    return out


def find_possible_non_independence(records) -> List[QCFinding]:
    """Same lab, different studies, same (population, valence, outcome, direction)."""
    buckets: Dict[Tuple, List] = {}
    for record in records:
        lab = lab_key(record)
        if lab is None:
            continue
        buckets.setdefault(lab + claim_key(record) + (norm(record.direction),), []).append(
            record
        )

    out: List[QCFinding] = []
    for key, group in sorted(buckets.items(), key=lambda kv: tuple(str(p) for p in kv[0])):
        studies = sorted({study_key(r) for r in group if r.study_id})
        if len(studies) < 2:
            continue
        out.append(
            QCFinding(
                code="QC_POSSIBLE_NON_INDEPENDENCE",
                severity=QCSeverity.INFO,
                message=(
                    "Lab %r reports the same claim (%s, valence=%s, outcome=%s, "
                    "direction=%s) in %d separate studies (%s). These records "
                    "are plausibly not independent evidence; group them in any "
                    "train/test split (SCHEMA.md 5). This is not an error and "
                    "nothing should be removed."
                    % (
                        key[0],
                        key[2],
                        key[3],
                        key[4],
                        key[5],
                        len(studies),
                        ", ".join(studies),
                    )
                ),
                record_ids=sorted(r.record_id for r in group),
                details={
                    "lab": list(key[:2]),
                    "studies": studies,
                    "n_records": len(group),
                },
            )
        )
    return out


def run_duplicate_checks(records) -> List[QCFinding]:
    return (
        find_exact_duplicates(records)
        + find_locator_collisions(records)
        + find_possible_non_independence(records)
    )
