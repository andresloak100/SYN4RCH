"""The direction-of-effect prediction task, and the rules that define its data.

    given (population_id, valence, outcome_id, species, acuity)
    predict direction in {increase, decrease, no_change}

This is deliberately a three-way *classification* of the sign of an effect, not
a quantitative simulation. ``docs/SCHEMA.md`` 6 refuses to represent
dose-response curves at all (the corpus does not report them consistently), and
refuses any aggregate "rejuvenation score". A sign is the strongest claim the
dataset can support, and whether even that is predictable is the open question
this module exists to answer.

What this module will not do
----------------------------
**It never silently drops a record.** Two of the five ``Direction`` values
cannot be turned into a label of this task:

* ``not_reported`` -- the paper measured something and did not say which way it
  went. Guessing would fabricate the observation.
* ``non_monotonic`` -- a real, recorded response shape with no single sign.
  ``docs/DECISIONS.md`` D-011 already refuses to treat it as contradicting a
  monotonic direction; collapsing it into one here would be the same error.

Both are **excluded and counted**, with the reason attached, in
:class:`ExclusionLedger`. A dataset whose exclusions are not reported is a
dataset whose label distribution cannot be interpreted: dropping
``not_reported`` rows quietly is how a benchmark comes to look better than the
literature it was built from.

**It never trains on Tier B.** ``docs/SCHEMA.md`` 6: "Tier B is an annotation
layer, not a training target." Tier B outcomes (lifespan, frailty, grip
strength) are therefore excluded by default and counted under their own reason,
not merged into some other bucket. ``exclude_tier_b=False`` exists only so a
caller can *measure* how much Tier B data there is; using it to train is a
schema violation, and the ledger records that the switch was flipped.

Grouping
--------
Every grouping this task supports is declared in :data:`GROUPINGS` and the
group keys are computed here, at dataset-build time, from
``syn4rch.qc.keys`` -- the same functions QC uses, so that "same lab" and
"same population" mean one thing across the project. ``docs/SCHEMA.md`` 5
makes an unlabelled cross-validated number an invalid result, so the grouping
registry is the single place a grouping name becomes legitimate, and
``syn4rch.models.evaluate`` refuses any name that is not in it.

Python 3.10 compatible, standard library only.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from ..qc.keys import lab_key, norm, population_key, study_key
from ..schema.enums import Direction, OutcomeTier

__all__ = [
    "LABELS",
    "LABELLABLE_DIRECTIONS",
    "UNLABELLABLE_DIRECTIONS",
    "FEATURE_FIELDS",
    "NULL_FEATURE",
    "ExclusionReason",
    "EXCLUSION_REASONS",
    "ExclusionLedger",
    "Features",
    "Example",
    "Grouping",
    "GROUPINGS",
    "DEFAULT_GROUPINGS",
    "DirectionTask",
    "TASK",
    "TaskDataset",
    "build_dataset",
    "load_dataset",
    "default_records_path",
    "repo_root",
]

#: The three labels of this task, in a fixed order. Every probability vector,
#: confusion matrix and per-class table in this package uses this order, so
#: that two reports can be compared column by column.
LABELS: Tuple[str, ...] = ("increase", "decrease", "no_change")

#: ``Direction`` values that map onto :data:`LABELS`.
LABELLABLE_DIRECTIONS: Tuple[Direction, ...] = (
    Direction.increase,
    Direction.decrease,
    Direction.no_change,
)

#: ``Direction`` values that are valid records but not labels of this task.
UNLABELLABLE_DIRECTIONS: Tuple[Direction, ...] = (
    Direction.not_reported,
    Direction.non_monotonic,
)

#: The feature tuple, in the order the task statement gives it.
FEATURE_FIELDS: Tuple[str, ...] = (
    "population_id",
    "valence",
    "outcome_id",
    "species",
    "acuity",
)

#: Stand-in for an absent feature value inside a composite cell key. It is a
#: visible token rather than an empty string so that a cell key printed in a
#: report cannot be mistaken for a populated one.
NULL_FEATURE = "<null>"


# --------------------------------------------------------------------------- #
# exclusions
# --------------------------------------------------------------------------- #

class ExclusionReason:
    """Reason codes for records that cannot become examples of this task.

    Plain string constants rather than an ``Enum``: these codes are written
    into reports and JSON, and are compared against strings in tests.
    """

    direction_absent = "direction_absent"
    direction_not_reported = "direction_not_reported"
    direction_non_monotonic = "direction_non_monotonic"
    outcome_id_absent = "outcome_id_absent"
    outcome_tier_b = "outcome_tier_b"
    population_absent = "population_absent"
    valence_absent = "valence_absent"


#: Reason code -> (one-line explanation, the clause that justifies excluding it).
EXCLUSION_REASONS: Dict[str, Tuple[str, str]] = {
    ExclusionReason.direction_absent: (
        "direction is null: the record states no effect at all",
        "SCHEMA.md 2 rule 4 (a field not reported is null) + 3.5",
    ),
    ExclusionReason.direction_not_reported: (
        "direction is 'not_reported': measured, but the paper does not say "
        "which way it went, so no label exists to predict",
        "SCHEMA.md 3.5",
    ),
    ExclusionReason.direction_non_monotonic: (
        "direction is 'non_monotonic': a recorded response shape with no "
        "single sign; collapsing it to one would be an interpretation",
        "SCHEMA.md 3.5 + DECISIONS.md D-011",
    ),
    ExclusionReason.outcome_id_absent: (
        "outcome_id is null, so the readout cannot be identified and the "
        "feature tuple is incomplete",
        "SCHEMA.md 3.5 (outcome_id required)",
    ),
    ExclusionReason.outcome_tier_b: (
        "outcome is Tier B (aging/lifespan/healthspan), which this project "
        "refuses to use as a training target",
        "SCHEMA.md 6",
    ),
    ExclusionReason.population_absent: (
        "neither population_id nor population_raw is present, so the "
        "population cannot be identified even by the paper's own wording",
        "SCHEMA.md 3.2 + DECISIONS.md D-010",
    ),
    ExclusionReason.valence_absent: (
        "valence is null; note that valence 'unclear' is NOT excluded, it is "
        "a legitimate feature value (SCHEMA.md 3.4 forbids coercing it)",
        "SCHEMA.md 3.4",
    ),
}

#: Order in which reasons are tested. The first match is the record's *primary*
#: reason (what the headline ledger counts); every reason that applied is also
#: counted separately, so a record excluded for two reasons shows up in both
#: tallies and neither is hidden behind the other.
_REASON_ORDER: Tuple[str, ...] = (
    ExclusionReason.direction_absent,
    ExclusionReason.direction_not_reported,
    ExclusionReason.direction_non_monotonic,
    ExclusionReason.outcome_tier_b,
    ExclusionReason.outcome_id_absent,
    ExclusionReason.population_absent,
    ExclusionReason.valence_absent,
)


@dataclass
class ExclusionLedger:
    """Counted, reasoned account of every input record that is not an example.

    ``primary`` partitions the excluded records (each counted exactly once).
    ``co_occurring`` counts every reason that applied to any record, so it can
    exceed ``n_excluded``. ``record_ids`` keeps the identifiers, capped, so a
    curator can go and look.
    """

    n_input: int = 0
    n_examples: int = 0
    primary: Dict[str, int] = field(default_factory=dict)
    co_occurring: Dict[str, int] = field(default_factory=dict)
    record_ids: Dict[str, List[str]] = field(default_factory=dict)
    max_ids_per_reason: int = 25
    tier_b_excluded: bool = True

    @property
    def n_excluded(self) -> int:
        return sum(self.primary.values())

    def add(self, record_id: Optional[str], reasons: Sequence[str]) -> None:
        if not reasons:
            raise ValueError("ExclusionLedger.add called with no reason")
        unknown = [r for r in reasons if r not in EXCLUSION_REASONS]
        if unknown:
            raise ValueError("undeclared exclusion reason(s): %r" % (unknown,))
        primary = reasons[0]
        self.primary[primary] = self.primary.get(primary, 0) + 1
        for reason in reasons:
            self.co_occurring[reason] = self.co_occurring.get(reason, 0) + 1
        bucket = self.record_ids.setdefault(primary, [])
        if record_id and len(bucket) < self.max_ids_per_reason:
            bucket.append(record_id)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_input_records": self.n_input,
            "n_examples": self.n_examples,
            "n_excluded": self.n_excluded,
            "excluded_fraction": (
                None if self.n_input == 0 else self.n_excluded / self.n_input
            ),
            "tier_b_excluded": self.tier_b_excluded,
            "primary_reason_counts": dict(sorted(self.primary.items())),
            "co_occurring_reason_counts": dict(sorted(self.co_occurring.items())),
            "reason_explanations": {
                code: {"why": why, "clause": clause}
                for code, (why, clause) in EXCLUSION_REASONS.items()
            },
            "example_record_ids": {
                code: list(ids) for code, ids in sorted(self.record_ids.items())
            },
        }


# --------------------------------------------------------------------------- #
# features, labels, examples
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Features:
    """The feature tuple of one example.

    ``population_id`` is the *resolved* vocabulary id and is ``None`` when the
    paper's population did not resolve (SCHEMA.md 3.2 forbids force-mapping
    it). ``population_key`` is the grouping key from ``syn4rch.qc.keys``, which
    falls back to ``"raw:" + normalised population_raw`` (DECISIONS.md D-010)
    so two unresolved populations never collide into one cell.

    The two are kept apart on purpose: the frequency baselines key on
    ``population_key`` (an unresolved population is still a population), while
    :class:`~syn4rch.models.baselines.NeurotransmitterSign` needs
    ``population_id`` because a transmitter can only be read out of the
    vocabulary, and must abstain when there is no id to look up.
    """

    population_id: Optional[str]
    population_key: str
    valence: Optional[str]
    outcome_id: Optional[str]
    species: Optional[str]
    acuity: Optional[str]

    def cell(self, *fields_: str) -> Tuple[str, ...]:
        """A hashable cell key over the named fields, nulls made visible."""
        out: List[str] = []
        for name in fields_:
            if name not in ("population_key",) and name not in FEATURE_FIELDS:
                raise KeyError("not a feature of this task: %r" % name)
            value = getattr(self, name)
            out.append(NULL_FEATURE if value in (None, "") else str(value))
        return tuple(out)

    def to_dict(self) -> Dict[str, Optional[str]]:
        return {
            "population_id": self.population_id,
            "population_key": self.population_key,
            "valence": self.valence,
            "outcome_id": self.outcome_id,
            "species": self.species,
            "acuity": self.acuity,
        }


@dataclass(frozen=True)
class Example:
    """One labelled example: features, label, and its group keys.

    The group keys are computed once, here, so that the evaluation harness
    never needs the record back and therefore cannot quietly group by something
    that is not declared in :data:`GROUPINGS`.
    """

    record_id: Optional[str]
    features: Features
    label: str
    group_keys: Mapping[str, Optional[str]]

    def group(self, grouping: str) -> Optional[str]:
        if grouping not in GROUPINGS:
            raise KeyError(
                "undeclared grouping %r; declared: %s"
                % (grouping, ", ".join(sorted(GROUPINGS)))
            )
        return self.group_keys.get(grouping)


def extract_features(record) -> Features:
    """Feature tuple of a ``TierARecord``. Never invents a value."""
    return Features(
        population_id=(record.population_id or None),
        population_key=population_key(record),
        valence=(str(record.valence) if record.valence is not None else None),
        outcome_id=(record.outcome_id or None),
        species=(str(record.species) if record.species is not None else None),
        acuity=(str(record.acuity) if record.acuity is not None else None),
    )


def extract_label(record) -> Optional[str]:
    """The task label, or ``None`` when the record carries no label.

    ``None`` means "this record cannot be labelled", which is a different
    statement from any of the three labels and must not be collapsed into one.
    """
    direction = record.direction
    if direction is None:
        return None
    if direction in LABELLABLE_DIRECTIONS:
        return str(direction)
    return None


def exclusion_reasons(record, exclude_tier_b: bool = True) -> List[str]:
    """Every reason this record cannot be an example, in :data:`_REASON_ORDER`.

    Empty list means the record is usable. ``valence == 'unclear'`` is usable:
    SCHEMA.md 3.4 calls it "a valid and expected answer" and forbids coercing
    it, so it travels as a feature value and the baselines may abstain on it.
    """
    found: List[str] = []
    direction = record.direction
    if direction is None:
        found.append(ExclusionReason.direction_absent)
    elif direction is Direction.not_reported:
        found.append(ExclusionReason.direction_not_reported)
    elif direction is Direction.non_monotonic:
        found.append(ExclusionReason.direction_non_monotonic)

    if exclude_tier_b and record.outcome_tier is OutcomeTier.B:
        found.append(ExclusionReason.outcome_tier_b)
    if not record.outcome_id:
        found.append(ExclusionReason.outcome_id_absent)
    if not record.population_id and not norm(record.population_raw):
        found.append(ExclusionReason.population_absent)
    if record.valence is None:
        found.append(ExclusionReason.valence_absent)

    order = {code: i for i, code in enumerate(_REASON_ORDER)}
    return sorted(set(found), key=lambda c: order[c])


# --------------------------------------------------------------------------- #
# groupings (SCHEMA.md 5)
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Grouping(object):
    """A declared way to group rows so that a split does not leak.

    ``leakage_prone`` marks a grouping that SCHEMA.md 5 warns against relying
    on. It is kept available because a number computed under it is a useful
    *upper bound* -- and because the honest way to show that row-level
    cross-validation flatters a model is to print it next to the grouped
    number, labelled, rather than to omit it.
    """

    name: str
    description: str
    schema_clause: str
    leakage_prone: bool = False
    may_be_ungroupable: bool = False
    #: The name of the equivalent generator in ``syn4rch.feasibility.splits``,
    #: or ``None`` when that module has no counterpart. Used only to *look*
    #: for a fold generator there; a generator found this way is still
    #: validated against the actual data before it is used.
    feasibility_name: Optional[str] = None


def _lab_group(record) -> Optional[str]:
    key = lab_key(record)
    if key is None:
        return None
    return "%s | %s" % key


def _population_outcome_group(record) -> Optional[str]:
    outcome = record.outcome_id
    if not outcome:
        return None
    return "%s x %s" % (population_key(record), norm(outcome))


#: Every grouping the evaluation harness will accept. A name absent from here
#: cannot be used, which is how SCHEMA.md 5's "naming the grouping" is made
#: structural rather than a convention someone remembers.
GROUPINGS: Dict[str, Grouping] = {
    "study_id": Grouping(
        name="study_id",
        description="no study appears in both train and test",
        schema_clause="SCHEMA.md 5 item 1",
        feasibility_name="by_study",
    ),
    "lab": Grouping(
        name="lab",
        description=(
            "no (last_author + affiliation) appears in both train and test; "
            "records with no lab identity cannot be grouped and are reported, "
            "not quietly pooled into one anonymous lab"
        ),
        schema_clause="SCHEMA.md 5 item 2 + DECISIONS.md D-009",
        may_be_ungroupable=True,
        feasibility_name="by_lab",
    ),
    "population": Grouping(
        name="population",
        description="leave-one-population-out: generalisation to a new population",
        schema_clause="SCHEMA.md 5 item 4",
        feasibility_name="leave_one_population_out",
    ),
    "population_outcome": Grouping(
        name="population_outcome",
        description=(
            "leave-one-(population x outcome)-cell-out: generalisation to a "
            "claim never seen, which is the condition a hypothesis generator "
            "actually operates under"
        ),
        schema_clause="SCHEMA.md 5 item 4",
        may_be_ungroupable=True,
    ),
    "record": Grouping(
        name="record",
        description=(
            "each row is its own group, i.e. ordinary row-level "
            "cross-validation. Leakage-prone by construction: the same study "
            "and the same lab appear on both sides"
        ),
        schema_clause="SCHEMA.md 5 (the thing it forbids reporting unlabelled)",
        leakage_prone=True,
    ),
}

#: Groupings the CLI reports by default, in report order.
DEFAULT_GROUPINGS: Tuple[str, ...] = (
    "study_id",
    "lab",
    "population",
    "population_outcome",
    "record",
)

_GROUP_KEY_FNS = {
    "study_id": lambda r: (study_key(r) or None),
    "lab": _lab_group,
    "population": lambda r: (population_key(r) or None),
    "population_outcome": _population_outcome_group,
    "record": lambda r: (r.record_id or None),
}


def group_keys_for(record) -> Dict[str, Optional[str]]:
    """All declared group keys for one record. ``None`` = cannot be grouped."""
    return {name: _GROUP_KEY_FNS[name](record) for name in GROUPINGS}


# --------------------------------------------------------------------------- #
# the task object
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class DirectionTask:
    """The task definition, as data.

    Exists so that a report, a saved prediction file and a test all quote the
    same statement of what was being predicted, rather than three paraphrases.
    """

    name: str = "direction_of_effect"
    version: str = "0.1.0"
    statement: str = (
        "given (population_id, valence, outcome_id, species, acuity), "
        "predict direction in {increase, decrease, no_change}"
    )
    labels: Tuple[str, ...] = LABELS
    feature_fields: Tuple[str, ...] = FEATURE_FIELDS
    is_classification: bool = True
    refuses: Tuple[str, ...] = (
        "quantitative effect size prediction (SCHEMA.md 6)",
        "any aggregate 'rejuvenation score' (SCHEMA.md 6)",
        "Tier B outcomes as a training target (SCHEMA.md 6)",
        "a cross-validated number without its grouping (SCHEMA.md 5)",
    )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "task": self.name,
            "task_version": self.version,
            "statement": self.statement,
            "labels": list(self.labels),
            "feature_fields": list(self.feature_fields),
            "refuses": list(self.refuses),
        }


#: The singleton task definition.
TASK = DirectionTask()


# --------------------------------------------------------------------------- #
# dataset
# --------------------------------------------------------------------------- #

@dataclass
class TaskDataset:
    """Labelled examples plus the counted account of what was left out."""

    examples: List[Example] = field(default_factory=list)
    #: The ``TierARecord`` each example came from, aligned index-for-index with
    #: :attr:`examples`. Kept so that a fold generator which takes *records*
    #: -- as ``syn4rch.feasibility.splits`` does -- can be handed the real
    #: records rather than a reconstruction of them, and so that a fold it
    #: returns can be mapped back to example indices without guesswork.
    records: List[Any] = field(default_factory=list)
    ledger: ExclusionLedger = field(default_factory=ExclusionLedger)
    source: Optional[str] = None
    load_errors: List[str] = field(default_factory=list)
    task: DirectionTask = TASK

    def __len__(self) -> int:
        return len(self.examples)

    @property
    def is_empty(self) -> bool:
        return not self.examples

    def labels(self) -> List[str]:
        return [e.label for e in self.examples]

    def features(self) -> List[Features]:
        return [e.features for e in self.examples]

    def label_counts(self) -> Dict[str, int]:
        counts = {label: 0 for label in LABELS}
        for example in self.examples:
            counts[example.label] += 1
        return counts

    def group_values(self, grouping: str) -> List[Optional[str]]:
        """Group key per example, aligned with :attr:`examples`."""
        if grouping not in GROUPINGS:
            raise KeyError(
                "undeclared grouping %r; declared: %s"
                % (grouping, ", ".join(sorted(GROUPINGS)))
            )
        return [e.group_keys.get(grouping) for e in self.examples]

    def grouping_coverage(self, grouping: str) -> Dict[str, Any]:
        """How many rows and distinct groups a grouping can actually cover."""
        values = self.group_values(grouping)
        groupable = [v for v in values if v is not None]
        sizes: Dict[str, int] = {}
        for value in groupable:
            sizes[value] = sizes.get(value, 0) + 1
        return {
            "grouping": grouping,
            "n_rows_total": len(values),
            "n_rows_groupable": len(groupable),
            "n_rows_ungroupable": len(values) - len(groupable),
            "n_groups": len(sizes),
            "largest_group_rows": max(sizes.values()) if sizes else 0,
            "singleton_groups": sum(1 for n in sizes.values() if n == 1),
            "leakage_prone": GROUPINGS[grouping].leakage_prone,
        }

    def subset(self, indices: Iterable[int]) -> List[Example]:
        return [self.examples[i] for i in indices]

    def record_index(self) -> Dict[str, int]:
        """``record_id -> example index``.

        A record id that occurs twice among the examples would make this map
        ambiguous, so it raises: QC detects exact duplicates
        (``syn4rch.qc.duplicates``) and a modelling layer that silently kept
        one of two identical rows would hide that finding.
        """
        out: Dict[str, int] = {}
        for index, example in enumerate(self.examples):
            rid = example.record_id
            if not rid:
                continue
            if rid in out:
                raise ValueError(
                    "duplicate record_id %r among the examples; run "
                    "syn4rch.qc first" % rid
                )
            out[rid] = index
        return out

    def describe(self) -> Dict[str, Any]:
        return {
            "task": self.task.to_dict(),
            "source": self.source,
            "n_examples": len(self.examples),
            "label_counts": self.label_counts(),
            "exclusions": self.ledger.to_dict(),
            "groupings": {
                name: self.grouping_coverage(name) for name in GROUPINGS
            },
            "load_errors": list(self.load_errors),
        }


def build_dataset(
    records: Iterable[Any],
    *,
    exclude_tier_b: bool = True,
    source: Optional[str] = None,
) -> TaskDataset:
    """Turn ``TierARecord`` objects into a :class:`TaskDataset`.

    Every input record either becomes an example or appears in the ledger with
    at least one declared reason. There is no third outcome; that invariant is
    asserted below rather than documented, because a counted exclusion that
    silently becomes an uncounted one is exactly the failure this ledger
    exists to prevent.
    """
    ledger = ExclusionLedger(tier_b_excluded=exclude_tier_b)
    examples: List[Example] = []
    kept: List[Any] = []
    for record in records:
        ledger.n_input += 1
        reasons = exclusion_reasons(record, exclude_tier_b=exclude_tier_b)
        label = extract_label(record)
        if reasons or label is None:
            if not reasons:
                # Defensive: a label of None with no reason would be a silent
                # drop. Attribute it rather than lose it.
                reasons = [ExclusionReason.direction_absent]
            ledger.add(getattr(record, "record_id", None), reasons)
            continue
        kept.append(record)
        examples.append(
            Example(
                record_id=getattr(record, "record_id", None),
                features=extract_features(record),
                label=label,
                group_keys=group_keys_for(record),
            )
        )
    ledger.n_examples = len(examples)
    if ledger.n_examples + ledger.n_excluded != ledger.n_input:
        raise AssertionError(
            "exclusion ledger does not balance: %d examples + %d excluded "
            "!= %d input records" % (
                ledger.n_examples, ledger.n_excluded, ledger.n_input
            )
        )
    return TaskDataset(
        examples=examples, records=kept, ledger=ledger, source=source
    )


# --------------------------------------------------------------------------- #
# locating and loading the curated table
# --------------------------------------------------------------------------- #

def repo_root() -> Optional[str]:
    """Directory holding ``data/`` and ``pyproject.toml``, or ``None``.

    Walks up from this file, as ``syn4rch.vocab.core.vocab_dir`` does. Returns
    ``None`` instead of guessing, so a caller in an unusual layout gets a clear
    "I do not know where the dataset is" rather than a path that does not exist.
    """
    import os

    current = os.path.abspath(os.path.dirname(__file__))
    for _ in range(8):
        if os.path.isdir(os.path.join(current, "data")) and os.path.isfile(
            os.path.join(current, "pyproject.toml")
        ):
            return current
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent
    return None


def default_records_path() -> Optional[str]:
    """``<repo>/data/curated/records.jsonl``, whether or not it exists."""
    import os

    root = repo_root()
    if root is None:
        return None
    return os.path.join(root, "data", "curated", "records.jsonl")


def load_dataset(
    path: Optional[str] = None,
    *,
    exclude_tier_b: bool = True,
    any_version: bool = False,
) -> TaskDataset:
    """Read the curated JSONL and build the dataset.

    A missing or empty file is **not** an error: it yields an empty dataset
    whose ``source`` records the path that was looked at. The caller (the CLI,
    the report) is then responsible for saying so. Reading one unparseable line
    stops the read and records the message in ``load_errors``; the records read
    before it are kept, because a partial dataset that says it is partial is
    more useful than an exception.
    """
    import os

    from ..schema.io import JsonlError, iter_jsonl

    if path is None:
        path = default_records_path()
    if path is None:
        dataset = build_dataset([], exclude_tier_b=exclude_tier_b, source=None)
        dataset.load_errors.append(
            "could not locate the repository root, so data/curated/records.jsonl "
            "could not be found; pass an explicit path"
        )
        return dataset

    if not os.path.isfile(path):
        dataset = build_dataset([], exclude_tier_b=exclude_tier_b, source=path)
        dataset.load_errors.append("no such file: %s" % path)
        return dataset

    records: List[Any] = []
    errors: List[str] = []
    try:
        for record in iter_jsonl(
            path, expect_version=None if any_version else "0.1.0"
        ):
            records.append(record)
    except JsonlError as exc:
        errors.append(str(exc))

    dataset = build_dataset(records, exclude_tier_b=exclude_tier_b, source=path)
    dataset.load_errors.extend(errors)
    return dataset
