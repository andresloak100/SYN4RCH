"""Leakage-aware fold generators, and a guard against unlabelled metrics.

Implements SCHEMA.md 5. Three groupings, each leave-one-group-out:

* :func:`by_study` -- item 1, ``study_id``.
* :func:`by_lab` -- item 2, the normalised ``last_author`` (D-056). The one
  SCHEMA.md calls the largest leakage risk in this corpus.
* :func:`leave_one_population_out` -- item 4, the generalisation question:
  can the model say anything about a population it has never seen?

Item 3 of SCHEMA.md 5 (shared animal cohorts across papers) is **not**
implemented and cannot be: nothing in the schema records a cohort identity, so
a split here can group by study and by lab but cannot promise that two papers
did not share mice. Every fold set says so in :attr:`FoldSet.note`, rather than
leaving a reader to assume the guarantee is complete. This matches the
deferral already recorded in docs/DECISIONS.md D-013.

**Records that cannot be grouped are dropped from the split entirely** -- not
assigned to train, not assigned to test -- and listed in
:attr:`FoldSet.ungroupable_record_ids`. A record with no lab identity
(docs/DECISIONS.md D-009) might belong to any lab, including the held-out one;
putting it in train is a possible leak and putting it in test is a possible
duplicate of train. Excluding it costs data and states the cost. It is never
lumped with other unknown-lab records into one anonymous "lab", which would
invent a group.

**The grouping label is not optional.** SCHEMA.md 5: "Reporting a single
cross-validated number without naming the grouping is considered an invalid
result in this project." So this module does not hand out bare folds. Folds
arrive inside a :class:`FoldSet` that carries its ``grouping``, and the only
way to attach a measurement to them is :func:`attach_metric`, which raises
:class:`UnlabelledResultError` unless the grouping travels with the number.
:class:`GroupedMetric` deliberately defines no ``__float__``: a caller that
wants the naked number has to reach for ``.value`` and will be reading the
``grouping`` field in the same line.

Python 3.10 compatible, standard library only.
"""

from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

from ..qc.keys import lab_key, population_key, study_key
from .coverage import UNRESOLVED_POPULATION_PREFIX

__all__ = [
    "GROUPINGS",
    "Fold",
    "FoldSet",
    "GroupedMetric",
    "UnlabelledResultError",
    "all_fold_sets",
    "attach_metric",
    "by_lab",
    "by_study",
    "leave_one_population_out",
    "require_grouping",
]

#: The only grouping labels this project recognises. A metric labelled with
#: anything else is refused: an unrecognised label is as uninformative as none.
GROUPINGS: Tuple[str, ...] = ("by_study", "by_lab", "leave_one_population_out")

_SHARED_COHORT_CAVEAT = (
    "Shared animal cohorts across papers (SCHEMA.md 5 item 3) are NOT "
    "controlled: the schema records no cohort identity, so this split cannot "
    "rule out two studies reporting one cohort."
)

_NOTES: Dict[str, str] = {
    "by_study": (
        "Leave-one-study-out on study_id (SCHEMA.md 5 item 1). Does not "
        "control for one lab publishing the same circuit repeatedly, which "
        "SCHEMA.md 5 item 2 calls the largest leakage risk here: use by_lab "
        "for that. " + _SHARED_COHORT_CAVEAT
    ),
    "by_lab": (
        "Leave-one-lab-out on the normalised last_author (SCHEMA.md 5 item 2; "
        "docs/DECISIONS.md D-056 -- the affiliation was part of this key until "
        "then, and made it identical to by_study). Records with no lab "
        "identity are excluded from every fold and listed; they are not "
        "pooled into an anonymous lab (docs/DECISIONS.md D-009). "
        + _SHARED_COHORT_CAVEAT
    ),
    "leave_one_population_out": (
        "Leave-one-population-out (SCHEMA.md 5 item 4), grouping on "
        "population_id or, when unresolved, 'raw:' + population_raw "
        "(docs/DECISIONS.md D-010). Tests generalisation to an unseen "
        "population; it does NOT also hold out the study or the lab, so a "
        "held-out population studied by a lab present in train remains a "
        "leakage path. " + _SHARED_COHORT_CAVEAT
    ),
}


class UnlabelledResultError(ValueError):
    """A cross-validated number was offered without its grouping.

    SCHEMA.md 5 declares such a number an invalid result, so this is an error
    and not a warning.
    """


@dataclass(frozen=True)
class Fold:
    """One leave-one-group-out fold, as record ids.

    Record ids rather than records: a fold is an index into the dataset, and
    copying records into folds invites a caller to mutate one copy.
    """

    index: int
    grouping: str
    held_out: str
    train_record_ids: Tuple[str, ...] = ()
    test_record_ids: Tuple[str, ...] = ()

    @property
    def n_train(self) -> int:
        return len(self.train_record_ids)

    @property
    def n_test(self) -> int:
        return len(self.test_record_ids)

    @property
    def is_degenerate(self) -> bool:
        """True if either side is empty: nothing can be fitted or scored."""
        return self.n_train == 0 or self.n_test == 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "index": self.index,
            "grouping": self.grouping,
            "held_out": self.held_out,
            "n_train": self.n_train,
            "n_test": self.n_test,
            "train_record_ids": list(self.train_record_ids),
            "test_record_ids": list(self.test_record_ids),
        }


@dataclass(frozen=True)
class FoldSet:
    """Folds plus the grouping label that makes them interpretable."""

    grouping: str
    folds: Tuple[Fold, ...] = ()
    groups: Tuple[str, ...] = ()
    ungroupable_record_ids: Tuple[str, ...] = ()
    note: str = ""

    def __post_init__(self) -> None:
        if self.grouping not in GROUPINGS:
            raise UnlabelledResultError(
                "grouping %r is not one of %s; an unrecognised grouping label "
                "is as uninformative as none (SCHEMA.md 5)"
                % (self.grouping, ", ".join(GROUPINGS))
            )

    @property
    def n_folds(self) -> int:
        return len(self.folds)

    @property
    def n_groups(self) -> int:
        return len(self.groups)

    @property
    def n_ungroupable(self) -> int:
        return len(self.ungroupable_record_ids)

    @property
    def n_usable_folds(self) -> int:
        return sum(1 for fold in self.folds if not fold.is_degenerate)

    @property
    def is_usable(self) -> bool:
        """True if at least two folds can actually be trained and scored.

        One group means leave-one-out has nothing to leave out; a single
        degenerate fold is not cross-validation.
        """
        return self.n_usable_folds >= 2

    def to_dict(self) -> Dict[str, Any]:
        return {
            "grouping": self.grouping,
            "n_groups": self.n_groups,
            "n_folds": self.n_folds,
            "n_usable_folds": self.n_usable_folds,
            "is_usable": self.is_usable,
            "groups": list(self.groups),
            "n_ungroupable": self.n_ungroupable,
            "ungroupable_record_ids": list(self.ungroupable_record_ids),
            "note": self.note,
            "folds": [f.to_dict() for f in self.folds],
        }


def _leave_one_group_out(
    records: Iterable,
    grouping: str,
    key_of: Callable[[Any], Optional[str]],
) -> FoldSet:
    records = list(records)
    grouped: Dict[str, List[str]] = {}
    ungroupable: List[str] = []
    for record in records:
        key = key_of(record)
        rid = record.record_id or ""
        if key is None:
            ungroupable.append(rid)
            continue
        grouped.setdefault(key, []).append(rid)

    groups = tuple(sorted(grouped))
    folds: List[Fold] = []
    for index, held_out in enumerate(groups):
        test = tuple(sorted(grouped[held_out]))
        train = tuple(
            sorted(
                rid
                for group, ids in grouped.items()
                if group != held_out
                for rid in ids
            )
        )
        folds.append(
            Fold(
                index=index,
                grouping=grouping,
                held_out=held_out,
                train_record_ids=train,
                test_record_ids=test,
            )
        )
    return FoldSet(
        grouping=grouping,
        folds=tuple(folds),
        groups=groups,
        ungroupable_record_ids=tuple(sorted(ungroupable)),
        note=_NOTES[grouping],
    )


def _study_group(record) -> Optional[str]:
    key = study_key(record)
    return key or None


def _lab_group(record) -> Optional[str]:
    # D-056: the lab key is the normalised last author, and is already the
    # group label. Joining an affiliation onto it is what made this split
    # indistinguishable from by_study.
    return lab_key(record)


def _population_group(record) -> Optional[str]:
    key = population_key(record)
    if not key or key == UNRESOLVED_POPULATION_PREFIX:
        return None
    return key


def by_study(records: Iterable) -> FoldSet:
    """Leave-one-study-out folds (SCHEMA.md 5 item 1)."""
    return _leave_one_group_out(records, "by_study", _study_group)


def by_lab(records: Iterable) -> FoldSet:
    """Leave-one-lab-out folds (SCHEMA.md 5 item 2).

    Lab identity is ``last_author`` + ``last_author_affiliation``. Records
    without it are excluded from every fold and listed.
    """
    return _leave_one_group_out(records, "by_lab", _lab_group)


def leave_one_population_out(records: Iterable) -> FoldSet:
    """Leave-one-population-out folds (SCHEMA.md 5 item 4)."""
    return _leave_one_group_out(
        records, "leave_one_population_out", _population_group
    )


def all_fold_sets(records: Iterable) -> Dict[str, FoldSet]:
    """Every grouping, keyed by its label. ``records`` is iterated once each."""
    records = list(records)
    return {
        "by_study": by_study(records),
        "by_lab": by_lab(records),
        "leave_one_population_out": leave_one_population_out(records),
    }


def require_grouping(folds: Any) -> str:
    """Return the grouping label of ``folds``, or raise.

    Raises :class:`UnlabelledResultError` if ``folds`` is a bare sequence of
    folds, carries no ``grouping``, or carries one this project does not
    recognise.
    """
    if isinstance(folds, (list, tuple, set, frozenset)):
        raise UnlabelledResultError(
            "a bare sequence of folds carries no grouping label; pass the "
            "FoldSet returned by by_study / by_lab / "
            "leave_one_population_out. SCHEMA.md 5: a cross-validated number "
            "without its grouping is an invalid result."
        )
    grouping = getattr(folds, "grouping", None)
    if not grouping:
        raise UnlabelledResultError(
            "%s carries no grouping label; SCHEMA.md 5 makes an unlabelled "
            "cross-validated number an invalid result"
            % type(folds).__name__
        )
    if grouping not in GROUPINGS:
        raise UnlabelledResultError(
            "grouping %r is not one of %s" % (grouping, ", ".join(GROUPINGS))
        )
    if getattr(folds, "folds", None) is None:
        raise UnlabelledResultError(
            "%s has a grouping label but no folds; there is nothing the "
            "number could have been cross-validated over" % type(folds).__name__
        )
    return grouping


@dataclass(frozen=True)
class GroupedMetric:
    """A number that cannot be stated without the grouping it was computed on.

    No ``__float__`` and no ``__index__`` on purpose: there is no implicit
    conversion that would let this be formatted as a bare number by accident.
    """

    name: str
    value: float
    grouping: str
    n_folds: int
    n_usable_folds: int
    n_groups: int
    n_ungroupable: int = 0
    note: str = ""

    def __post_init__(self) -> None:
        if self.grouping not in GROUPINGS:
            raise UnlabelledResultError(
                "grouping %r is not one of %s" % (self.grouping, ", ".join(GROUPINGS))
            )
        if not self.name:
            raise UnlabelledResultError("a metric must be named")

    def __str__(self) -> str:
        text = "%s = %s [grouping=%s, folds=%d usable of %d, groups=%d" % (
            self.name,
            ("%.4g" % self.value) if isinstance(self.value, float) else self.value,
            self.grouping,
            self.n_usable_folds,
            self.n_folds,
            self.n_groups,
        )
        if self.n_ungroupable:
            text += ", %d record(s) excluded as ungroupable" % self.n_ungroupable
        return text + "]"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "value": self.value,
            "grouping": self.grouping,
            "n_folds": self.n_folds,
            "n_usable_folds": self.n_usable_folds,
            "n_groups": self.n_groups,
            "n_ungroupable": self.n_ungroupable,
            "note": self.note,
        }


def attach_metric(folds: Any, name: str, value: float) -> GroupedMetric:
    """Bind a measurement to the fold set it was computed over.

    The only supported way to report a cross-validated number in this project.
    Raises :class:`UnlabelledResultError` when the grouping label does not come
    with the folds -- including when a caller passes a plain list of folds, an
    object with an empty or unknown ``grouping``, or ``None``.
    """
    grouping = require_grouping(folds)
    return GroupedMetric(
        name=name,
        value=value,
        grouping=grouping,
        n_folds=getattr(folds, "n_folds", len(getattr(folds, "folds", ()) or ())),
        n_usable_folds=getattr(folds, "n_usable_folds", 0),
        n_groups=getattr(folds, "n_groups", 0),
        n_ungroupable=getattr(folds, "n_ungroupable", 0),
        note=getattr(folds, "note", "") or "",
    )
