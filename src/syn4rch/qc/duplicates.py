"""Duplicate, collision and non-independence detection.

Three distinct things, deliberately not conflated:

1. **Exact duplicates** -- byte-identical records (same ``record_id``, same
   content). A harvesting or merge bug. WARNING.
2. **Locator collisions** -- same study, same ``source_locator``, same claim,
   but *different* content. Either one figure panel was extracted twice with
   divergent values, or the locator is too coarse to identify the
   observation. WARNING; needs a human.
3. **Possible non-independence** -- different studies from the *same lab*
   reporting the same claim and the same direction. Not an error: it is normal
   and expected that a group follows up its own work. Reported as INFO so that
   leakage control (SCHEMA.md 5) can group these records, and never as
   something to delete.
4. **Lab identity collisions** -- one last author appearing under affiliations
   with nothing in common, i.e. the homonym that
   :func:`syn4rch.qc.keys.lab_key` can no longer rule out now that it is the
   author alone (D-056). WARNING, because merging two researchers into one lab
   costs test rows and manufactures a replication.

Check 3 returned **nothing on this corpus until D-056**, and not because the
corpus is clean: lab identity included the free-text affiliation, so four labs
that had published eleven papers between them counted as eleven labs and no
lab ever appeared in two studies. A detector that cannot fire is worse than an
absent one, because its silence reads as a result.

Nothing here removes, merges or edits a record.

Python 3.10 compatible, standard library only.
"""

from typing import Dict, List, Tuple

from .findings import QCFinding, QCSeverity
from .keys import CLAIM_KEY_FIELDS, claim_key, lab_affiliation, lab_key, norm, study_key

__all__ = [
    "AFFILIATION_STOPWORDS",
    "affiliation_tokens",
    "find_exact_duplicates",
    "find_lab_identity_collisions",
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
    """Same lab, different studies, same claim and same direction.

    The lab is :func:`syn4rch.qc.keys.lab_key`, the normalised last author
    (D-056). While that key also carried the affiliation this function could
    not fire at all, because no two papers in this corpus spell one department
    identically.
    """
    buckets: Dict[Tuple, List] = {}
    for record in records:
        lab = lab_key(record)
        if lab is None:
            continue
        key = (lab,) + claim_key(record) + (norm(record.direction),)
        buckets.setdefault(key, []).append(record)

    out: List[QCFinding] = []
    for key, group in sorted(buckets.items(), key=lambda kv: tuple(str(p) for p in kv[0])):
        studies = sorted({study_key(r) for r in group if r.study_id})
        if len(studies) < 2:
            continue
        # The claim components are named from CLAIM_KEY_FIELDS rather than
        # positionally. They were printed positionally until D-056, under the
        # wrong names: claim_key stopped carrying `valence` when v0.2.0 added
        # the perturbation axis, so the message had been labelling the
        # perturbation class "valence", the polarity "outcome" and the target
        # "direction" ever since. Nothing read it, because the check never
        # fired.
        claim = dict(zip(CLAIM_KEY_FIELDS, key[1 : 1 + len(CLAIM_KEY_FIELDS)]))
        direction = key[1 + len(CLAIM_KEY_FIELDS)]
        out.append(
            QCFinding(
                code="QC_POSSIBLE_NON_INDEPENDENCE",
                severity=QCSeverity.INFO,
                message=(
                    "Lab %r reports the same claim (%s) with direction=%s in "
                    "%d separate studies (%s). These records are plausibly not "
                    "independent evidence; group them in any train/test split "
                    "(SCHEMA.md 5 item 2). This is not an error and nothing "
                    "should be removed."
                    % (
                        key[0],
                        ", ".join(
                            "%s=%s" % (name, claim[name])
                            for name in CLAIM_KEY_FIELDS
                        ),
                        direction,
                        len(studies),
                        ", ".join(studies),
                    )
                ),
                record_ids=sorted(r.record_id for r in group),
                details={
                    "lab": key[0],
                    "claim": claim,
                    "direction": direction,
                    "studies": studies,
                    "n_records": len(group),
                },
            )
        )
    return out


#: Affiliation words that identify no particular institution. Two unrelated
#: laboratories both sit in a "department" at a "university" in a country, so a
#: shared token from this list is not evidence that two affiliations name one
#: place. Kept as data, and deliberately short: the check is meant to fire on a
#: genuine homonym, not to adjudicate how two addresses were typed.
AFFILIATION_STOPWORDS: frozenset = frozenset(
    """
    academy affiliated biology biomedical cancer center centre central chinese
    city clinic clinical college department departments division faculty
    graduate group health hospital institute institutes institution
    international laboratory laboratories medical medicine national
    neuroscience program programme research school science sciences state
    system systems technology unit university universities
    australia austria belgium brazil canada china denmark england finland
    france germany india ireland israel italy japan korea mexico netherlands
    norway poland portugal russia scotland singapore spain sweden switzerland
    taiwan turkey wales
    """.split()
)

#: Shortest affiliation token that can identify a place. Three letters and
#: under are postcodes, state abbreviations and prepositions.
MIN_AFFILIATION_TOKEN = 4


def affiliation_tokens(text) -> frozenset:
    """Substantive words of an affiliation: the place-identifying tokens.

    Alphanumeric runs of at least :data:`MIN_AFFILIATION_TOKEN` characters,
    minus :data:`AFFILIATION_STOPWORDS`. Street numbers and postcodes survive
    and that is fine -- they are highly specific when they match and simply
    absent when they do not.
    """
    out = set()
    word = []
    for char in norm(text):
        if char.isalnum():
            word.append(char)
            continue
        if word:
            out.add("".join(word))
            word = []
    if word:
        out.add("".join(word))
    return frozenset(
        token
        for token in out
        if len(token) >= MIN_AFFILIATION_TOKEN and token not in AFFILIATION_STOPWORDS
    )


def find_lab_identity_collisions(records) -> List[QCFinding]:
    """One last author whose affiliations have no substantive token in common.

    The failure mode D-056 accepted when it made :func:`lab_key` the author
    alone: two unrelated researchers whose normalised surname and initials
    coincide merge into one laboratory, which removes test rows from every
    leave-one-lab-out split and reports their independent results as one lab
    replicating itself.

    The test is deliberately weak, and warns rather than erring. An affiliation
    may be absent, abbreviated to an institute's initials, or written for a
    different campus of one university, and none of those is a homonym. A blank
    affiliation is skipped entirely: it is an absence of evidence about the lab,
    not evidence of two labs. An author with one affiliation, or none, cannot
    trip it.
    """
    by_author: Dict[str, Dict[str, List]] = {}
    for record in records:
        author = lab_key(record)
        if author is None:
            continue
        affiliation = lab_affiliation(record)
        if not affiliation:
            continue
        by_author.setdefault(author, {}).setdefault(affiliation, []).append(record)

    out: List[QCFinding] = []
    for author, affiliations in sorted(by_author.items()):
        if len(affiliations) < 2:
            continue
        token_sets = {
            affiliation: affiliation_tokens(affiliation)
            for affiliation in affiliations
        }
        if any(not tokens for tokens in token_sets.values()):
            # An affiliation with no substantive token at all says nothing
            # either way, so this author is not evidence of a collision.
            continue
        shared = frozenset.intersection(*token_sets.values())
        if shared:
            continue
        group = [r for records_ in affiliations.values() for r in records_]
        studies = sorted({study_key(r) for r in group if r.study_id})
        out.append(
            QCFinding(
                code="QC_POSSIBLE_LAB_HOMONYM",
                severity=QCSeverity.WARNING,
                message=(
                    "Last author %r appears under %d affiliations with no "
                    "substantive word in common, across %d studies (%s). "
                    "lab_key is the author alone (D-056), so these are being "
                    "grouped as one laboratory; if they are two researchers "
                    "who share a name, every leave-one-lab-out split is "
                    "dropping rows it should be testing and the replication "
                    "counts are overstated. Affiliations: %s"
                    % (
                        author,
                        len(affiliations),
                        len(studies),
                        ", ".join(studies),
                        " || ".join(sorted(affiliations)),
                    )
                ),
                record_ids=sorted(r.record_id for r in group),
                field_name="last_author_affiliation",
                details={
                    "lab": author,
                    "affiliations": sorted(affiliations),
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
        + find_lab_identity_collisions(records)
    )
