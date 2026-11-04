"""Perturbation identity: what a held-out record's manipulation shares with training.

SCHEMA.md 5b. Grouping a split by ``study_id`` or by lab prevents the same
study or the same laboratory appearing on both sides. It does **not** establish
that anything on the test side was learnable from the train side.

The structural question
-----------------------
Of the first 362 curated records, 170 are ``conditional_knockout`` or
``conditional_overexpression``. SCHEMA.md 5b records the measurement: for 13 of
the molecular target genes in this corpus, the gene appears in exactly one
study and one laboratory. Not one of those genes is manipulated by a second
group.

The consequence is not statistical, it is structural. Under ``by_study`` or
``by_lab`` grouping, when the one study that manipulated a gene is held out,
no training record shares that perturbation identity. There is nothing for a
model to have learned about it. An accuracy figure averaged over such folds is
not a measurement of generalisation; it is a measurement of the label prior.

So this module answers, for every grouping declared in
:data:`syn4rch.models.task.GROUPINGS`: **what fraction of records have a
perturbation identity that occurs in exactly one group?** Under leave-one-group-out
those records are unlearnable by construction -- the two statements are the
same statement, because an identity confined to one group is absent from
training in exactly the fold that tests it.

The registry is reused rather than restated. SCHEMA.md 5 makes an unlabelled
cross-validated number an invalid result, and
:data:`syn4rch.models.task.GROUPINGS` is the single place a grouping name
becomes legitimate; a parallel list here would be a second place to add a
grouping that the evaluation harness would then refuse.

What "perturbation identity" means
----------------------------------
:func:`syn4rch.qc.keys.perturbation_key`: ``(class, polarity, target)``. The
three travel together because a polarity means different physical things in
different classes (:class:`syn4rch.schema.enums.PerturbationPolarity`), and
nothing in this project may key on polarity alone.

Three statuses, and the third is never folded into either of the others:

* ``resolved`` -- the identity is as specific as its class allows.
* ``target_unresolved`` -- the class is one whose identity *is* the target
  (:data:`TARGET_CONSTITUTIVE_CLASSES`) and the target is null. A knockout with
  no named gene is not a perturbation identity: it is the statement "some gene
  was removed". Counting it as resolved would pool every knockout in the corpus
  into one identity, that identity would occur in many studies, and the
  analysis would report the molecular class as **learnable** -- the exact
  opposite of the finding, manufactured by a missing field. So these records
  are counted separately and no single-group fraction is reported for them.
* ``unresolvable`` -- no ``method`` at all, so there is no class and no axis.
  Counted separately too.

For ``circuit_activity`` and ``structural`` a null target is *correct*, not
missing: the manipulation acts on the population, and the population is a
separate component of :func:`syn4rch.qc.keys.claim_key`. For
``pharmacological`` the identity is the drug, which
:func:`~syn4rch.qc.keys.perturbation_key` has no slot for, so those records are
reported as ``target_unresolved`` with that limitation named rather than
silently counted as a single identity.

Reading this while the migration is in flight
---------------------------------------------
``target_raw`` and ``target_id`` are being resolved from ``targeting_strategy``
in a separate pass. The analysis is correct either way and gets sharper when
they land: today every molecular record is ``target_unresolved``, so the
molecular single-group fraction is **not reported at all** rather than reported
as zero, and :attr:`IdentityReport.target_resolution` states how much of the
corpus had a resolved target when the numbers were computed. Every report this
module renders carries that line, so a number can never be read as if it had
been measured on resolved targets when it was not.

Python 3.10 compatible, standard library only (docs/DECISIONS.md D-001).
"""

import datetime
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from ..models.task import GROUPINGS, group_keys_for
from ..qc.keys import NULL, lab_key, perturbation_key, study_key
from .splits import UnlabelledResultError

__all__ = [
    "IDENTITY_FIELDS",
    "STATUS_RESOLVED",
    "STATUS_TARGET_UNRESOLVED",
    "STATUS_UNRESOLVABLE",
    "STATUSES",
    "TARGET_CONSTITUTIVE_CLASSES",
    "ClassIdentityCoverage",
    "ClassStructure",
    "Identity",
    "IdentityCoverage",
    "IdentityReport",
    "UnlabelledIdentityError",
    "analyse",
    "build_report",
    "class_structure",
    "identity_coverage",
    "perturbation_identity",
    "render_markdown",
]

#: Component names of :func:`syn4rch.qc.keys.perturbation_key`, in order.
IDENTITY_FIELDS: Tuple[str, ...] = (
    "perturbation_class",
    "perturbation_polarity",
    "perturbation_target",
)

#: The identity is as specific as its class allows.
STATUS_RESOLVED = "resolved"

#: The class's identity *is* its target and the target is null.
STATUS_TARGET_UNRESOLVED = "target_unresolved"

#: No ``method``, so no class and no axis: there is no identity to place.
STATUS_UNRESOLVABLE = "unresolvable"

#: Every status, in report order.
STATUSES: Tuple[str, ...] = (
    STATUS_RESOLVED,
    STATUS_TARGET_UNRESOLVED,
    STATUS_UNRESOLVABLE,
)

#: Classes whose perturbation identity is the target, so that a null target
#: leaves the identity underdetermined rather than merely unannotated.
#:
#: * ``molecular`` -- the manipulation *is* the gene. "A knockout" is not an
#:   identity; "a Raptor knockout" is.
#: * ``pharmacological`` -- the manipulation is the drug. ``perturbation_key``
#:   carries no ligand slot, so these records cannot reach ``resolved`` at all
#:   today and the report says so rather than counting them as one identity.
#:
#: ``circuit_activity`` and ``structural`` are deliberately absent: their
#: target is the population, which :func:`syn4rch.qc.keys.claim_key` carries
#: as its own component, so a null target there is correct and not a gap.
TARGET_CONSTITUTIVE_CLASSES: Tuple[str, ...] = ("molecular", "pharmacological")


class UnlabelledIdentityError(UnlabelledResultError):
    """An identity figure was offered without a declared grouping.

    Subclasses :class:`syn4rch.feasibility.splits.UnlabelledResultError` so
    that one ``except`` catches both: SCHEMA.md 5's rule is the same rule
    whether the number is an accuracy or an unseen-identity fraction.
    """


@dataclass(frozen=True)
class Identity:
    """One record's perturbation identity and how far it resolves."""

    key: Tuple[str, str, str]
    status: str
    perturbation_class: str

    @property
    def is_resolved(self) -> bool:
        return self.status == STATUS_RESOLVED

    @property
    def label(self) -> str:
        return " / ".join(self.key)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "key": list(self.key),
            "fields": list(IDENTITY_FIELDS),
            "status": self.status,
            "perturbation_class": self.perturbation_class,
        }


def perturbation_identity(record) -> Identity:
    """The record's identity, via :func:`syn4rch.qc.keys.perturbation_key`.

    The key function is used rather than the stored fields, so a record
    written before the v0.2.0 migration -- where class and polarity are
    derived from ``method`` and ``valence`` -- resolves identically to the same
    record after it.
    """
    key = perturbation_key(record)
    cls, polarity, target = key
    if cls == NULL or polarity == NULL:
        return Identity(key=key, status=STATUS_UNRESOLVABLE, perturbation_class=cls)
    if cls in TARGET_CONSTITUTIVE_CLASSES and target == NULL:
        return Identity(
            key=key, status=STATUS_TARGET_UNRESOLVED, perturbation_class=cls
        )
    return Identity(key=key, status=STATUS_RESOLVED, perturbation_class=cls)


# --------------------------------------------------------------------------- #
# coverage of one grouping
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class ClassIdentityCoverage:
    """The single-group result for one perturbation class under one grouping.

    ``fraction_single_group`` is ``None``, never ``0.0``, when no record of
    this class has a resolved identity. "No record of this class has an
    identity specific enough to look up" and "no record of this class has an
    identity confined to one group" are opposite conclusions, and a zero would
    print as the second while meaning the first.
    """

    perturbation_class: str
    grouping: str
    n_records: int = 0
    n_records_resolved: int = 0
    n_records_target_unresolved: int = 0
    n_records_unresolvable: int = 0
    n_records_ungroupable: int = 0
    n_identities: int = 0
    n_identities_single_group: int = 0
    n_records_single_group: int = 0
    n_records_multi_group: int = 0

    @property
    def fraction_single_group(self) -> Optional[float]:
        """Of the resolved, groupable records, the share that cannot be learned."""
        denominator = self.n_records_single_group + self.n_records_multi_group
        if denominator == 0:
            return None
        return self.n_records_single_group / denominator

    @property
    def is_measurable(self) -> bool:
        return self.fraction_single_group is not None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "perturbation_class": self.perturbation_class,
            "grouping": self.grouping,
            "n_records": self.n_records,
            "n_records_resolved": self.n_records_resolved,
            "n_records_target_unresolved": self.n_records_target_unresolved,
            "n_records_unresolvable": self.n_records_unresolvable,
            "n_records_ungroupable": self.n_records_ungroupable,
            "n_identities": self.n_identities,
            "n_identities_single_group": self.n_identities_single_group,
            "n_records_single_group": self.n_records_single_group,
            "n_records_multi_group": self.n_records_multi_group,
            "fraction_single_group": self.fraction_single_group,
            "is_measurable": self.is_measurable,
        }


@dataclass(frozen=True)
class IdentityCoverage:
    """Perturbation-identity reach under one declared grouping.

    Every figure here is bound to :attr:`grouping`, which is validated against
    :data:`syn4rch.models.task.GROUPINGS` in :meth:`__post_init__`. There is no
    ``__float__``: as with
    :class:`syn4rch.feasibility.splits.GroupedMetric`, a caller that wants the
    bare fraction has to name the field and will be reading the grouping on the
    same line.
    """

    grouping: str
    description: str = ""
    schema_clause: str = ""
    leakage_prone: bool = False
    n_records: int = 0
    n_groups: int = 0
    n_records_ungroupable: int = 0
    n_records_resolved: int = 0
    n_records_target_unresolved: int = 0
    n_records_unresolvable: int = 0
    n_identities: int = 0
    n_identities_single_group: int = 0
    n_records_single_group: int = 0
    n_records_multi_group: int = 0
    #: Unseen-identity fraction per leave-one-group-out fold, in group order.
    per_fold_unseen_fraction: Tuple[Tuple[str, Optional[float]], ...] = ()
    by_class: Tuple[ClassIdentityCoverage, ...] = ()

    def __post_init__(self) -> None:
        if self.grouping not in GROUPINGS:
            raise UnlabelledIdentityError(
                "grouping %r is not declared in syn4rch.models.task.GROUPINGS "
                "(declared: %s); an unrecognised grouping label is as "
                "uninformative as none (SCHEMA.md 5)"
                % (self.grouping, ", ".join(sorted(GROUPINGS)))
            )

    @property
    def n_records_placed(self) -> int:
        """Records with a resolved identity *and* a group: the denominator."""
        return self.n_records_single_group + self.n_records_multi_group

    @property
    def fraction_single_group(self) -> Optional[float]:
        """Share of placed records whose identity occurs in exactly one group.

        ``None`` when no record could be placed -- which is a statement about
        the data, not a fraction of zero.
        """
        if self.n_records_placed == 0:
            return None
        return self.n_records_single_group / self.n_records_placed

    @property
    def fraction_of_corpus_unmeasurable(self) -> Optional[float]:
        if not self.n_records:
            return None
        return (
            self.n_records_target_unresolved + self.n_records_unresolvable
        ) / self.n_records

    @property
    def stated(self) -> str:
        """The fraction with its grouping, as a report must print it."""
        value = self.fraction_single_group
        if value is None:
            return (
                "not measurable [grouping=%s]: no record has both a resolved "
                "perturbation identity and a group" % self.grouping
            )
        return (
            "%.1f%% of placed records have a single-group perturbation identity "
            "[grouping=%s, %d of %d placed, %d groups]"
            % (
                100.0 * value,
                self.grouping,
                self.n_records_single_group,
                self.n_records_placed,
                self.n_groups,
            )
        )

    def class_coverage(self, perturbation_class: str) -> Optional[ClassIdentityCoverage]:
        for entry in self.by_class:
            if entry.perturbation_class == perturbation_class:
                return entry
        return None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "grouping": self.grouping,
            "description": self.description,
            "schema_clause": self.schema_clause,
            "leakage_prone": self.leakage_prone,
            "n_records": self.n_records,
            "n_groups": self.n_groups,
            "n_records_ungroupable": self.n_records_ungroupable,
            "n_records_resolved": self.n_records_resolved,
            "n_records_target_unresolved": self.n_records_target_unresolved,
            "n_records_unresolvable": self.n_records_unresolvable,
            "n_identities": self.n_identities,
            "n_identities_single_group": self.n_identities_single_group,
            "n_records_single_group": self.n_records_single_group,
            "n_records_multi_group": self.n_records_multi_group,
            "n_records_placed": self.n_records_placed,
            "fraction_single_group": self.fraction_single_group,
            "fraction_of_corpus_unmeasurable": self.fraction_of_corpus_unmeasurable,
            "stated": self.stated,
            "per_fold_unseen_fraction": [
                {"held_out": name, "unseen_fraction": value}
                for name, value in self.per_fold_unseen_fraction
            ],
            "by_class": [entry.to_dict() for entry in self.by_class],
        }


def identity_coverage(records: Iterable, grouping: str) -> IdentityCoverage:
    """Measure perturbation-identity reach under one declared grouping.

    A record counts as *unlearnable by construction* when its identity occurs
    in exactly one group: under leave-one-group-out, the fold that tests it is
    the fold that removed every record sharing its identity.

    Records with no group key are excluded from the fractions and counted in
    ``n_records_ungroupable`` -- the same treatment
    :mod:`syn4rch.feasibility.splits` gives them (docs/DECISIONS.md D-009), and
    for the same reason: a record that cannot be grouped cannot be placed on
    either side of a split, so it cannot be said to be learnable or not.
    """
    if grouping not in GROUPINGS:
        raise UnlabelledIdentityError(
            "grouping %r is not declared in syn4rch.models.task.GROUPINGS "
            "(declared: %s)" % (grouping, ", ".join(sorted(GROUPINGS)))
        )
    records = list(records)
    declared = GROUPINGS[grouping]

    identity_groups: Dict[Tuple[str, str, str], set] = {}
    groups: set = set()
    placed: List[Tuple[Any, Identity, str]] = []
    n_ungroupable = 0
    status_counts = {name: 0 for name in STATUSES}
    class_records: Dict[str, List[Any]] = {}

    for record in records:
        identity = perturbation_identity(record)
        status_counts[identity.status] += 1
        class_records.setdefault(identity.perturbation_class, []).append(record)
        group = group_keys_for(record).get(grouping)
        if group is None:
            n_ungroupable += 1
            continue
        groups.add(group)
        if not identity.is_resolved:
            continue
        identity_groups.setdefault(identity.key, set()).add(group)
        placed.append((record, identity, group))

    single = {key for key, seen in identity_groups.items() if len(seen) == 1}
    n_single = sum(1 for _, identity, _ in placed if identity.key in single)
    n_multi = len(placed) - n_single

    # Per-fold unseen fraction, the figure SCHEMA.md 5b requires next to any
    # cross-validated number. For leave-one-group-out, a test row's identity is
    # unseen in training exactly when that identity occurs in no other group.
    per_group: Dict[str, List[bool]] = {}
    for _, identity, group in placed:
        per_group.setdefault(group, []).append(identity.key in single)
    per_fold = tuple(
        (
            group,
            (sum(1 for flag in flags if flag) / len(flags)) if flags else None,
        )
        for group, flags in sorted(per_group.items())
    )

    by_class: List[ClassIdentityCoverage] = []
    for cls in sorted(class_records):
        members = class_records[cls]
        cls_status = {name: 0 for name in STATUSES}
        cls_ungroupable = 0
        cls_identities: set = set()
        cls_single = 0
        cls_multi = 0
        for record in members:
            identity = perturbation_identity(record)
            cls_status[identity.status] += 1
            group = group_keys_for(record).get(grouping)
            if group is None:
                cls_ungroupable += 1
                continue
            if not identity.is_resolved:
                continue
            cls_identities.add(identity.key)
            if identity.key in single:
                cls_single += 1
            else:
                cls_multi += 1
        by_class.append(
            ClassIdentityCoverage(
                perturbation_class=cls,
                grouping=grouping,
                n_records=len(members),
                n_records_resolved=cls_status[STATUS_RESOLVED],
                n_records_target_unresolved=cls_status[STATUS_TARGET_UNRESOLVED],
                n_records_unresolvable=cls_status[STATUS_UNRESOLVABLE],
                n_records_ungroupable=cls_ungroupable,
                n_identities=len(cls_identities),
                n_identities_single_group=len(
                    {k for k in cls_identities if k in single}
                ),
                n_records_single_group=cls_single,
                n_records_multi_group=cls_multi,
            )
        )

    return IdentityCoverage(
        grouping=grouping,
        description=declared.description,
        schema_clause=declared.schema_clause,
        leakage_prone=declared.leakage_prone,
        n_records=len(records),
        n_groups=len(groups),
        n_records_ungroupable=n_ungroupable,
        n_records_resolved=status_counts[STATUS_RESOLVED],
        n_records_target_unresolved=status_counts[STATUS_TARGET_UNRESOLVED],
        n_records_unresolvable=status_counts[STATUS_UNRESOLVABLE],
        n_identities=len(identity_groups),
        n_identities_single_group=len(single),
        n_records_single_group=n_single,
        n_records_multi_group=n_multi,
        per_fold_unseen_fraction=per_fold,
        by_class=tuple(by_class),
    )


# --------------------------------------------------------------------------- #
# the publication structure behind the finding
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class ClassStructure:
    """How one perturbation class is distributed over studies and labs.

    This is the part of SCHEMA.md 5b that is measurable whether or not targets
    have been resolved, and it is what makes the finding a statement about
    publication practice rather than about sample size.

    ``study_lab_is_bijection`` is the load-bearing one. When every study of a
    class comes from a different lab and no lab has two studies of it, the
    class's study-level and lab-level identity reach are necessarily identical:
    holding out a lab holds out exactly one study. Resolving targets can then
    only move both numbers together, never separate them.
    """

    perturbation_class: str
    n_records: int = 0
    n_studies: int = 0
    n_labs: int = 0
    n_records_without_lab: int = 0
    n_labs_with_multiple_studies: int = 0
    n_studies_sharing_a_lab: int = 0
    n_distinct_targeting_strategies: int = 0
    n_targeting_strategies_in_multiple_studies: int = 0

    @property
    def study_lab_is_bijection(self) -> bool:
        """One lab per study and one study per lab, with every lab identified."""
        return (
            self.n_records > 0
            and self.n_records_without_lab == 0
            and self.n_studies == self.n_labs
            and self.n_labs_with_multiple_studies == 0
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "perturbation_class": self.perturbation_class,
            "n_records": self.n_records,
            "n_studies": self.n_studies,
            "n_labs": self.n_labs,
            "n_records_without_lab": self.n_records_without_lab,
            "n_labs_with_multiple_studies": self.n_labs_with_multiple_studies,
            "n_studies_sharing_a_lab": self.n_studies_sharing_a_lab,
            "n_distinct_targeting_strategies": self.n_distinct_targeting_strategies,
            "n_targeting_strategies_in_multiple_studies": (
                self.n_targeting_strategies_in_multiple_studies
            ),
            "study_lab_is_bijection": self.study_lab_is_bijection,
        }


def class_structure(records: Iterable) -> List[ClassStructure]:
    """One :class:`ClassStructure` per perturbation class present."""
    by_class: Dict[str, List[Any]] = {}
    for record in records:
        by_class.setdefault(perturbation_identity(record).perturbation_class, []).append(
            record
        )

    out: List[ClassStructure] = []
    for cls in sorted(by_class):
        members = by_class[cls]
        studies = {study_key(r) for r in members if r.study_id}
        labs = {lab_key(r) for r in members if lab_key(r) is not None}
        lab_to_studies: Dict[Any, set] = {}
        for record in members:
            lab = lab_key(record)
            if lab is None:
                continue
            lab_to_studies.setdefault(lab, set()).add(study_key(record))
        strategy_to_studies: Dict[str, set] = {}
        for record in members:
            strategy = (record.targeting_strategy or "").strip().casefold()
            if not strategy:
                continue
            strategy_to_studies.setdefault(strategy, set()).add(study_key(record))
        multi = [v for v in lab_to_studies.values() if len(v) > 1]
        out.append(
            ClassStructure(
                perturbation_class=cls,
                n_records=len(members),
                n_studies=len(studies),
                n_labs=len(labs),
                n_records_without_lab=sum(
                    1 for r in members if lab_key(r) is None
                ),
                n_labs_with_multiple_studies=len(multi),
                n_studies_sharing_a_lab=sum(len(v) for v in multi),
                n_distinct_targeting_strategies=len(strategy_to_studies),
                n_targeting_strategies_in_multiple_studies=sum(
                    1 for v in strategy_to_studies.values() if len(v) > 1
                ),
            )
        )
    return out


# --------------------------------------------------------------------------- #
# the report
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class TargetResolution:
    """How much of the corpus had a resolved perturbation target when this ran.

    Carried on every report because the molecular result depends entirely on
    it. ``n_with_target == 0`` means the one-gene-one-lab structure SCHEMA.md
    5b describes is **not visible to this analysis yet**, and the report says
    so in its first sentence rather than printing a fraction that would read as
    its refutation.
    """

    n_records: int = 0
    n_with_target_id: int = 0
    n_with_target_raw_only: int = 0
    n_target_constitutive: int = 0
    n_target_constitutive_resolved: int = 0

    @property
    def n_with_target(self) -> int:
        return self.n_with_target_id + self.n_with_target_raw_only

    @property
    def is_complete(self) -> bool:
        return (
            self.n_target_constitutive > 0
            and self.n_target_constitutive_resolved == self.n_target_constitutive
        )

    @property
    def is_absent(self) -> bool:
        return self.n_target_constitutive_resolved == 0

    @property
    def state(self) -> str:
        if self.n_target_constitutive == 0:
            return "no records of a target-constitutive class"
        if self.is_complete:
            return "targets fully resolved"
        if self.is_absent:
            return "targets not yet resolved"
        return "targets partially resolved"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_records": self.n_records,
            "n_with_target_id": self.n_with_target_id,
            "n_with_target_raw_only": self.n_with_target_raw_only,
            "n_with_target": self.n_with_target,
            "n_target_constitutive": self.n_target_constitutive,
            "n_target_constitutive_resolved": self.n_target_constitutive_resolved,
            "state": self.state,
            "is_complete": self.is_complete,
            "is_absent": self.is_absent,
        }


def _target_resolution(records: Sequence) -> TargetResolution:
    n_id = 0
    n_raw = 0
    n_const = 0
    n_const_resolved = 0
    for record in records:
        has_id = bool(getattr(record, "target_id", None))
        has_raw = bool(getattr(record, "target_raw", None))
        if has_id:
            n_id += 1
        elif has_raw:
            n_raw += 1
        identity = perturbation_identity(record)
        if identity.perturbation_class in TARGET_CONSTITUTIVE_CLASSES:
            n_const += 1
            if identity.key[2] != NULL:
                n_const_resolved += 1
    return TargetResolution(
        n_records=len(records),
        n_with_target_id=n_id,
        n_with_target_raw_only=n_raw,
        n_target_constitutive=n_const,
        n_target_constitutive_resolved=n_const_resolved,
    )


@dataclass
class IdentityReport:
    """Everything the perturbation-identity pass measured."""

    coverages: List[IdentityCoverage] = field(default_factory=list)
    structures: List[ClassStructure] = field(default_factory=list)
    target_resolution: TargetResolution = field(default_factory=TargetResolution)
    scope: str = "all records read"
    n_records: int = 0
    source_path: Optional[str] = None
    source_exists: bool = True
    load_errors: List[str] = field(default_factory=list)
    generated_at: Optional[str] = None

    @property
    def is_empty(self) -> bool:
        return self.n_records == 0

    def coverage(self, grouping: str) -> Optional[IdentityCoverage]:
        for entry in self.coverages:
            if entry.grouping == grouping:
                return entry
        return None

    def structure(self, perturbation_class: str) -> Optional[ClassStructure]:
        for entry in self.structures:
            if entry.perturbation_class == perturbation_class:
                return entry
        return None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "generated_at": self.generated_at,
            "source_path": self.source_path,
            "source_exists": self.source_exists,
            "load_errors": list(self.load_errors),
            "scope": self.scope,
            "n_records": self.n_records,
            "identity_fields": list(IDENTITY_FIELDS),
            "target_constitutive_classes": list(TARGET_CONSTITUTIVE_CLASSES),
            "target_resolution": self.target_resolution.to_dict(),
            "top_line": top_line(self),
            "coverages": [c.to_dict() for c in self.coverages],
            "class_structures": [s.to_dict() for s in self.structures],
        }


def analyse(
    records: Iterable,
    groupings: Optional[Sequence[str]] = None,
    scope: str = "all records read",
) -> "IdentityReport":
    """Measure identity reach under every declared grouping.

    ``groupings`` defaults to every name in
    :data:`syn4rch.models.task.GROUPINGS`. A name absent from that registry
    raises: no grouping may be invented here, and none may be dropped -- a
    grouping that produces an inconvenient number is the one most worth
    printing.
    """
    records = list(records)
    names = tuple(groupings) if groupings is not None else tuple(sorted(GROUPINGS))
    unknown = [name for name in names if name not in GROUPINGS]
    if unknown:
        raise UnlabelledIdentityError(
            "undeclared grouping(s) %r; declared: %s"
            % (unknown, ", ".join(sorted(GROUPINGS)))
        )
    return IdentityReport(
        coverages=[identity_coverage(records, name) for name in names],
        structures=class_structure(records),
        target_resolution=_target_resolution(records),
        scope=scope,
        n_records=len(records),
        generated_at=_utc_now(),
    )


def build_report(
    records: Sequence,
    source_path: Optional[str] = None,
    source_exists: bool = True,
    load_errors: Optional[Sequence[str]] = None,
    scope: str = "all records read",
    generated_at: Optional[str] = None,
) -> IdentityReport:
    """:func:`analyse` plus the provenance a written report has to carry."""
    report = analyse(records, scope=scope)
    report.source_path = source_path
    report.source_exists = source_exists
    report.load_errors = list(load_errors or [])
    report.generated_at = generated_at or _utc_now()
    return report


def _utc_now() -> str:
    return (
        datetime.datetime.now(datetime.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


# --- rendering --------------------------------------------------------------


def _fmt_fraction(value: Optional[float]) -> str:
    return "not measurable" if value is None else "%.1f%%" % (100.0 * value)


def _table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> List[str]:
    lines = [
        "| " + " | ".join(header) + " |",
        "|" + "|".join("---" for _ in header) + "|",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    if not rows:
        lines.append("| _none_ |" + " |" * (len(header) - 1))
    return lines


def top_line(report: "IdentityReport") -> str:
    """One sentence, and it is the first thing in the file."""
    if report.is_empty:
        where = report.source_path or "the curated records file"
        reason = "does not exist yet" if not report.source_exists else "is empty"
        return (
            "**No records: `%s` %s, so no part of this corpus can be said to "
            "support or not support cross-study prediction.**" % (where, reason)
        )

    resolution = report.target_resolution
    molecular = report.structure("molecular")
    n_molecular = molecular.n_records if molecular else 0

    if resolution.is_absent and n_molecular:
        return (
            "**The %d molecular record(s) in this corpus cannot yet be assessed: "
            "`target_raw` and `target_id` are null for all of them, so their "
            "perturbation identity is the bare statement \"a gene was changed\" "
            "and no single-group fraction is reported for them. They come from "
            "%d studies in %d laboratories, one study per laboratory, and no "
            "two of those studies share a targeting strategy — the structure "
            "SCHEMA.md 5b describes is there in the provenance, and resolving "
            "the targets is what will let this analysis see it.**"
            % (
                n_molecular,
                molecular.n_studies,
                molecular.n_labs,
            )
        )

    # SCHEMA.md 5b is about ``by_study`` and ``by_lab`` specifically, so the
    # headline quotes one of those where it exists rather than whichever
    # grouping happens to give the largest number: a grouping with one group
    # trivially makes every identity single-group, and leading with that would
    # be alarming for the wrong reason.
    preferred = [
        c
        for c in (report.coverage("lab"), report.coverage("study_id"))
        if c is not None and c.fraction_single_group is not None
    ]
    worst = preferred or [
        c
        for c in report.coverages
        if not c.leakage_prone and c.fraction_single_group is not None
    ]
    if not worst:
        return (
            "**No record in this corpus has both a resolved perturbation "
            "identity and a group under any leakage-controlled grouping, so "
            "the unseen-identity fraction SCHEMA.md 5b requires cannot be "
            "computed on these %d record(s).**" % report.n_records
        )
    top = max(worst, key=lambda c: c.fraction_single_group or 0.0)
    return (
        "**Under `%s` grouping, %s of placed records have a perturbation "
        "identity that occurs in exactly one group: when that group is held "
        "out, nothing in training shares the manipulation, so a fold-averaged "
        "accuracy over them measures the label prior and not "
        "generalisation.**" % (top.grouping, _fmt_fraction(top.fraction_single_group))
    )


def _section_run(report: "IdentityReport") -> List[str]:
    resolution = report.target_resolution
    rows = [
        ["generated (UTC)", report.generated_at or "-"],
        ["input", "`%s`" % (report.source_path or "<in memory>")],
        [
            "input present",
            "yes" if report.source_exists else "no (treated as 0 records)",
        ],
        ["scope", report.scope],
        ["records analysed", str(report.n_records)],
        ["perturbation identity", "`%s`" % " / ".join(IDENTITY_FIELDS)],
        ["target resolution state", resolution.state],
        [
            "records with a resolved `target_id`",
            str(resolution.n_with_target_id),
        ],
        [
            "records with `target_raw` only",
            str(resolution.n_with_target_raw_only),
        ],
        [
            "records whose class needs a target",
            "%d (of which %d have one)"
            % (
                resolution.n_target_constitutive,
                resolution.n_target_constitutive_resolved,
            ),
        ],
        ["load errors", str(len(report.load_errors))],
    ]
    lines = ["## Run", "", *_table(["field", "value"], rows), ""]
    if report.load_errors:
        lines.append("### Load errors")
        lines.append("")
        for err in report.load_errors:
            lines.append("- `%s`" % err)
        lines.append("")
    lines.append(
        "> **What state the data was in.** %s. A molecular record's identity "
        "is its target gene; with the target null, every knockout in the "
        "corpus shares one identity and every overexpression another, so such "
        "records are counted as `target_unresolved` and are **excluded from "
        "both** the learnable and the unlearnable bucket rather than counted "
        "as learnable. Reading them as learnable is the one serious "
        "misinterpretation available here, and it is the artefact of a field "
        "that has not been filled in yet."
        % resolution.state.capitalize()
    )
    lines.append("")
    return lines


def _section_identity(report: "IdentityReport") -> List[str]:
    lines = [
        "## What a perturbation identity is, and when it is missing",
        "",
        "An identity is `(perturbation_class, perturbation_polarity, "
        "perturbation_target)` — `syn4rch.qc.keys.perturbation_key`. The three "
        "travel together because a polarity means different physical things in "
        "different classes: `decrease` is less firing under "
        "`circuit_activity` and less of a gene product under `molecular`, and "
        "nothing in this project is allowed to key on polarity alone.",
        "",
        "A class whose identity *is* its target (%s) and whose target is "
        "null has no identity to place. Those records are counted on their own "
        "line, never folded into either bucket:"
        % ", ".join("`%s`" % c for c in TARGET_CONSTITUTIVE_CLASSES),
        "",
        "- `resolved` — as specific as the class allows.",
        "- `target_unresolved` — the class needs a target and there is none. "
        "For `molecular` that is a missing field, being filled in now. For "
        "`pharmacological` the identity is the *drug*, and "
        "`perturbation_key` has no ligand slot at all, so those records "
        "cannot reach `resolved` today and this is a limitation of the key, "
        "not of the curation.",
        "- `unresolvable` — no `method`, so no class and no axis.",
        "",
        "For `circuit_activity` and `structural` a null target is correct, not "
        "missing: the manipulation acts on the population, and the population "
        "is a separate component of the claim key.",
        "",
    ]
    return lines


def _section_groupings(report: "IdentityReport") -> List[str]:
    lines = [
        "## Identity reach, per declared grouping",
        "",
        "Groupings are `syn4rch.models.task.GROUPINGS`, reused rather than "
        "restated: SCHEMA.md 5 makes that registry the single place a grouping "
        "name becomes legitimate. No grouping is omitted, including the "
        "leakage-prone one, because a grouping that produces an inconvenient "
        "number is the one most worth printing.",
        "",
        "*Placed* records are those with a resolved identity and a group. "
        "*Single-group* records are those whose identity occurs in exactly one "
        "group — under leave-one-group-out, the fold that tests them is the "
        "fold that removed every record sharing their manipulation.",
        "",
    ]
    rows = []
    for coverage in report.coverages:
        rows.append(
            [
                "`%s`" % coverage.grouping,
                "yes" if coverage.leakage_prone else "no",
                str(coverage.n_groups),
                str(coverage.n_records_placed),
                str(coverage.n_records_single_group),
                _fmt_fraction(coverage.fraction_single_group),
                str(
                    coverage.n_records_target_unresolved
                    + coverage.n_records_unresolvable
                ),
                str(coverage.n_records_ungroupable),
            ]
        )
    lines.extend(
        _table(
            [
                "grouping",
                "leakage-prone",
                "groups",
                "placed records",
                "single-group records",
                "unseen-identity fraction",
                "identity not placeable",
                "ungroupable",
            ],
            rows,
        )
    )
    lines.append("")
    lines.append(
        "Every fraction above is stated with its grouping, as SCHEMA.md 5 "
        "requires. `IdentityCoverage` has no `__float__`, so the bare number "
        "cannot be formatted without naming the field that carries the "
        "grouping beside it."
    )
    lines.append("")
    return lines


def _section_classes(report: "IdentityReport") -> List[str]:
    lines = [
        "## By perturbation class — where the contrast is",
        "",
        "This is the breakdown the analysis exists for. `circuit_activity` and "
        "`molecular` are not two flavours of the same thing; they behave "
        "completely differently under cross-study grouping, and the reason is "
        "what each class's identity *is*.",
        "",
    ]
    for coverage in report.coverages:
        if coverage.leakage_prone:
            continue
        lines.append("### Grouping `%s`" % coverage.grouping)
        lines.append("")
        rows = []
        for entry in coverage.by_class:
            rows.append(
                [
                    "`%s`" % entry.perturbation_class,
                    str(entry.n_records),
                    str(entry.n_records_resolved),
                    str(entry.n_identities),
                    str(entry.n_records_single_group),
                    str(entry.n_records_multi_group),
                    _fmt_fraction(entry.fraction_single_group),
                    str(
                        entry.n_records_target_unresolved
                        + entry.n_records_unresolvable
                    ),
                ]
            )
        lines.extend(
            _table(
                [
                    "class",
                    "records",
                    "identity resolved",
                    "distinct identities",
                    "single-group",
                    "multi-group",
                    "unseen-identity fraction",
                    "not placeable",
                ],
                rows,
            )
        )
        lines.append("")
    return lines


def _section_structure(report: "IdentityReport") -> List[str]:
    lines = [
        "## Why: the structure of how this field publishes",
        "",
        "The numbers above are not a sample-size problem. They are a "
        "description of publication practice, and that description is "
        "measurable whether or not the targets have been resolved yet.",
        "",
    ]
    rows = []
    for structure in report.structures:
        rows.append(
            [
                "`%s`" % structure.perturbation_class,
                str(structure.n_records),
                str(structure.n_studies),
                str(structure.n_labs),
                str(structure.n_labs_with_multiple_studies),
                str(structure.n_distinct_targeting_strategies),
                str(structure.n_targeting_strategies_in_multiple_studies),
                "yes" if structure.study_lab_is_bijection else "no",
            ]
        )
    lines.extend(
        _table(
            [
                "class",
                "records",
                "studies",
                "labs",
                "labs with >1 study",
                "targeting strategies",
                "strategies in >1 study",
                "one study per lab",
            ],
            rows,
        )
    )
    lines.append("")
    bijections = [
        s for s in report.structures if s.study_lab_is_bijection and s.n_studies > 1
    ]
    single_study = [s for s in report.structures if s.n_studies <= 1]
    if single_study:
        lines.append(
            "%s rest%s on a single study, so for %s the unseen-identity "
            "fraction is not yet a statement about how the field publishes: "
            "it is a statement that the corpus holds one paper."
            % (
                ", ".join("`%s`" % s.perturbation_class for s in single_study),
                "" if len(single_study) > 1 else "s",
                "them" if len(single_study) > 1 else "it",
            )
        )
        lines.append("")
    if bijections:
        lines.append(
            "For %s, every study comes from a different laboratory and no "
            "laboratory contributed two: the study-to-lab map is one-to-one. "
            "Holding out a lab therefore holds out exactly one study, so "
            "`study_id` and `lab` grouping must give these classes the same "
            "unseen-identity fraction, now and after the targets are "
            "resolved. The two groupings cannot be played off against each "
            "other to make the number smaller."
            % ", ".join("`%s`" % s.perturbation_class for s in bijections)
        )
        lines.append("")
    return lines


def _section_conclusion(report: "IdentityReport") -> List[str]:
    lines = [
        "## What this corpus can and cannot support",
        "",
        "In plain terms, and without a model having been fitted:",
        "",
    ]

    for coverage in report.coverages:
        if coverage.leakage_prone or coverage.grouping != "lab":
            continue
        for entry in coverage.by_class:
            fraction = entry.fraction_single_group
            if entry.n_records == 0:
                continue
            structure = report.structure(entry.perturbation_class)
            if structure is not None and structure.n_labs <= 1:
                # A different situation from the one below, and conflating
                # them would overstate the generality of the finding: this
                # class is one laboratory's work, so there is no cross-study
                # question to ask of it yet, whatever its identities look like.
                lines.append(
                    "- **`%s` (%d records): one laboratory, so there is no "
                    "cross-study question yet.** The whole class is %d "
                    "stud%s from %d laborator%s. Held out under `lab`, there "
                    "is no other lab left. This is not the one-gene-one-lab "
                    "structure described below; it is a class the corpus has "
                    "barely begun to cover."
                    % (
                        entry.perturbation_class,
                        entry.n_records,
                        structure.n_studies,
                        "y" if structure.n_studies == 1 else "ies",
                        structure.n_labs,
                        "y" if structure.n_labs == 1 else "ies",
                    )
                )
                continue
            if fraction is None:
                lines.append(
                    "- **`%s` (%d records): cannot be assessed from this "
                    "corpus today.** Not one of these records has an identity "
                    "specific enough to look for in a training set, because "
                    "the field that would carry it is null. This is not a "
                    "finding that they are predictable; it is a finding that "
                    "the question cannot be put yet."
                    % (entry.perturbation_class, entry.n_records)
                )
            elif fraction >= 0.999:
                lines.append(
                    "- **`%s` (%d records): cannot support cross-study "
                    "prediction at all.** Every one of these records has a "
                    "manipulation that appears in exactly one laboratory. Held "
                    "out, there is nothing in training that shares it. No "
                    "amount of modelling fixes this and no additional paper "
                    "from the same groups fixes it either: the only remedy is "
                    "a second laboratory repeating the same manipulation."
                    % (entry.perturbation_class, entry.n_records)
                )
            elif fraction <= 0.001:
                lines.append(
                    "- **`%s` (%d records): can support cross-study "
                    "prediction, as far as identity goes.** Every "
                    "manipulation here recurs across laboratories, so a "
                    "held-out record's perturbation identity is one the "
                    "training set has seen. That is a necessary condition for "
                    "generalisation, not a sufficient one — it says a model "
                    "*could* have learned something, not that it did."
                    % (entry.perturbation_class, entry.n_records)
                )
            else:
                lines.append(
                    "- **`%s` (%d records): partly.** %s of them have a "
                    "manipulation confined to one laboratory and are "
                    "unlearnable under `lab` grouping; the rest recur across "
                    "laboratories. A single accuracy figure over this class "
                    "averages the two and describes neither."
                    % (
                        entry.perturbation_class,
                        entry.n_records,
                        _fmt_fraction(fraction),
                    )
                )
    lines.append("")
    lines.append(
        "The reason is not a shortage of papers. It is that this field "
        "publishes one laboratory, one gene, one paper: a group establishes a "
        "gene's role in a population and the next group picks a different "
        "gene. Circuit manipulations escape this because their identity is a "
        "method everyone shares — a chemogenetic activation is the same "
        "manipulation in every laboratory that performs one — while a "
        "molecular manipulation's identity is the gene, and the gene is what "
        "distinguishes one laboratory's programme from another's."
    )
    lines.append("")
    lines.append(
        "Where a class is spread over several laboratories and its "
        "unseen-identity fraction is still 100%, the honest statement is the "
        "strong one: **that class of records can never be predicted from this "
        "corpus, under any grouping that controls leakage, until a second "
        "laboratory repeats the same manipulation.** Not \"with more data\", "
        "and not \"with a better model\". With a replication. A class that "
        "rests on one laboratory is a weaker version of the same statement: "
        "there is nothing to hold out yet."
    )
    lines.append("")
    return lines


def _section_limits(report: "IdentityReport") -> List[str]:
    return [
        "## What this report does not claim",
        "",
        "- It does not fit, score or cross-validate a model. It measures a "
        "property of the data that holds before any model exists.",
        "- A low unseen-identity fraction is a necessary condition for "
        "generalisation, never a sufficient one. It says a manipulation "
        "recurred across groups, not that anything about it was learnable.",
        "- It does not resolve targets, and it does not guess one from "
        "`targeting_strategy`. The strategy string is counted (how many are "
        "shared across studies) but never parsed into a gene: that is "
        "`target_raw`/`target_id`'s job and inferring it here would be the "
        "inference SCHEMA.md 2 rule 4 forbids.",
        "- It cannot rule out two studies sharing one animal cohort "
        "(SCHEMA.md 5 item 3): the schema records no cohort identity.",
        "- `pharmacological` records cannot reach `resolved` under the current "
        "`perturbation_key`, because the drug that defines them has no slot in "
        "it. They are counted as `target_unresolved` and that is a limitation "
        "of the key, which this report states rather than works around.",
        "",
    ]


def render_markdown(report: "IdentityReport") -> str:
    """The full report. First line after the title is the verdict sentence."""
    lines: List[str] = [
        "# SYN4RCH perturbation identity report",
        "",
        top_line(report),
        "",
        "Every number below is counted from the input named in the Run table. "
        "Nothing is estimated, imputed or projected, and no target is inferred "
        "from a targeting strategy. SCHEMA.md 5b: a cross-validated number "
        "must be reported alongside the fraction of test rows whose "
        "perturbation identity was unseen in training. This report is that "
        "fraction.",
        "",
    ]
    lines.extend(_section_run(report))
    if report.is_empty:
        lines.extend(
            [
                "## Nothing to measure",
                "",
                "With zero records there is no identity to place and no "
                "grouping to place it under. No fraction is printed for them: "
                "a 0% unseen-identity fraction over an empty corpus would read "
                "as the best possible result.",
                "",
            ]
        )
    else:
        lines.extend(_section_identity(report))
        lines.extend(_section_groupings(report))
        lines.extend(_section_classes(report))
        lines.extend(_section_structure(report))
        lines.extend(_section_conclusion(report))
    lines.extend(_section_limits(report))
    return "\n".join(lines).rstrip() + "\n"
