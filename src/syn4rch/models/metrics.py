"""Scalar metric functions. Standard library only, deterministic.

Everything here returns a **bare number or ``None``** -- never a dict. That is
deliberate. ``docs/SCHEMA.md`` 5 makes a cross-validated number reported
without its grouping an invalid result, so the only objects in this package
with a ``to_dict()`` are the grouping-labelled containers in
``syn4rch.models.evaluate``. If these functions returned ready-to-serialise
dicts, an unlabelled metric would be one ``json.dump`` away; as it is, getting
a number out of this module and into a file requires passing through the
labelled wrapper.

``None`` always means "this metric is not defined on this input", and the
caller is expected to say why rather than substitute a value. Returning 0.0 for
an undefined accuracy, or 0.5 for an undefined AUC, is how an empty evaluation
comes to look like a mediocre one.

No optional library is used; see ``_compat.py`` for why.

Python 3.10 compatible.
"""

import math
import random
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from .task import LABELS

__all__ = [
    "PredPair",
    "scored_pairs",
    "accuracy",
    "confusion_matrix",
    "per_class_prf",
    "macro_f1",
    "brier_score",
    "ReliabilityBin",
    "reliability_bins",
    "expected_calibration_error",
    "roc_auc_ovr",
    "roc_auc_ovr_macro",
    "scores_are_constant",
    "bootstrap_group_ci",
    "mcnemar_exact_p",
    "percentile",
]

#: (true label, predicted label or None) pairs.
PredPair = Tuple[str, Optional[str]]


def scored_pairs(
    y_true: Sequence[str], y_pred: Sequence[Optional[str]]
) -> List[Tuple[str, str]]:
    """Pairs where the model actually answered.

    An abstention (``None``) is removed rather than counted as an error.
    Scoring an abstention as wrong would punish a model for saying "I do not
    know", which on this corpus is frequently the correct answer; scoring it as
    right would be absurd. The cost of removing it is that accuracy is then
    computed on a subset, so :func:`~syn4rch.models.evaluate` reports coverage
    beside every accuracy and the report prints the two together.
    """
    if len(y_true) != len(y_pred):
        raise ValueError(
            "y_true and y_pred differ in length (%d vs %d)" % (len(y_true), len(y_pred))
        )
    return [(t, p) for t, p in zip(y_true, y_pred) if p is not None]


def accuracy(y_true: Sequence[str], y_pred: Sequence[Optional[str]]) -> Optional[float]:
    """Accuracy over answered examples, or ``None`` if there are none."""
    pairs = scored_pairs(y_true, y_pred)
    if not pairs:
        return None
    return sum(1 for t, p in pairs if t == p) / len(pairs)


def confusion_matrix(
    y_true: Sequence[str],
    y_pred: Sequence[Optional[str]],
    labels: Sequence[str] = LABELS,
) -> Dict[str, Dict[str, int]]:
    """``{true_label: {predicted_label: n}}``, with an ``abstain`` column."""
    out = {
        t: {p: 0 for p in list(labels) + ["abstain"]} for t in labels
    }
    for true, pred in zip(y_true, y_pred):
        row = out.setdefault(true, {p: 0 for p in list(labels) + ["abstain"]})
        key = pred if pred is not None else "abstain"
        row[key] = row.get(key, 0) + 1
    return out


def per_class_prf(
    y_true: Sequence[str],
    y_pred: Sequence[Optional[str]],
    labels: Sequence[str] = LABELS,
) -> Dict[str, Dict[str, Optional[float]]]:
    """Precision, recall, F1 and support per class.

    Precision is ``None`` (not 0.0) when the class was never predicted, and
    recall is ``None`` when the class never occurs in the truth. A class that
    was never predicted has no precision; calling it zero would silently drag
    a macro average down for a reason the average cannot express.
    """
    out: Dict[str, Dict[str, Optional[float]]] = {}
    for label in labels:
        tp = sum(1 for t, p in zip(y_true, y_pred) if p == label and t == label)
        fp = sum(1 for t, p in zip(y_true, y_pred) if p == label and t != label)
        support = sum(1 for t in y_true if t == label)
        fn = support - tp
        precision = None if (tp + fp) == 0 else tp / (tp + fp)
        recall = None if support == 0 else tp / support
        if precision is None or recall is None or (precision + recall) == 0:
            f1: Optional[float] = (
                None if (precision is None or recall is None) else 0.0
            )
        else:
            f1 = 2 * precision * recall / (precision + recall)
        out[label] = {
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "support": support,
            "tp": tp,
            "fp": fp,
            "fn": fn,
        }
    return out


def macro_f1(
    y_true: Sequence[str],
    y_pred: Sequence[Optional[str]],
    labels: Sequence[str] = LABELS,
) -> Optional[float]:
    """Unweighted mean F1 over classes **present in the truth**.

    A class with no support is excluded from the average rather than entered as
    zero: on a fold where ``no_change`` never occurs, counting it as zero would
    cap macro-F1 at 0.67 for reasons that belong to the split, not the model.
    A class that is present in the truth but never predicted contributes 0,
    which is the honest penalty for never predicting it.
    """
    stats = per_class_prf(y_true, y_pred, labels)
    present = [label for label in labels if stats[label]["support"] > 0]
    if not present:
        return None
    total = 0.0
    for label in present:
        value = stats[label]["f1"]
        total += 0.0 if value is None else value
    return total / len(present)


def brier_score(
    y_true: Sequence[str],
    probas: Sequence[Optional[Dict[str, float]]],
    labels: Sequence[str] = LABELS,
) -> Optional[float]:
    """Multiclass Brier score: mean over examples of sum_k (p_k - y_k)^2.

    Range 0 (perfect) to 2 (confidently wrong on a one-hot target). Examples
    with no probability vector are skipped, not imputed; ``None`` comes back if
    none of them had one.
    """
    total = 0.0
    n = 0
    for true, proba in zip(y_true, probas):
        if proba is None:
            continue
        n += 1
        for label in labels:
            target = 1.0 if label == true else 0.0
            total += (float(proba.get(label, 0.0)) - target) ** 2
    if n == 0:
        return None
    return total / n


class ReliabilityBin(object):
    """One confidence bin of a reliability diagram.

    ``meaningful`` is False when the bin holds fewer than ``min_n`` examples.
    A reliability diagram with two examples in a bin is a diagram that says
    nothing about calibration in that range, and plotting the point anyway is
    the standard way calibration gets over-read.
    """

    __slots__ = ("low", "high", "n", "mean_confidence", "observed_accuracy", "min_n")

    def __init__(self, low, high, n, mean_confidence, observed_accuracy, min_n):
        self.low = low
        self.high = high
        self.n = n
        self.mean_confidence = mean_confidence
        self.observed_accuracy = observed_accuracy
        self.min_n = min_n

    @property
    def meaningful(self) -> bool:
        return self.n >= self.min_n

    @property
    def gap(self) -> Optional[float]:
        if self.mean_confidence is None or self.observed_accuracy is None:
            return None
        return self.observed_accuracy - self.mean_confidence

    def as_row(self) -> Dict[str, object]:
        """Plain values for a labelled container to embed. Not a metric dict:
        it carries no number that stands alone without its grouping."""
        return {
            "bin_low": self.low,
            "bin_high": self.high,
            "n": self.n,
            "mean_confidence": self.mean_confidence,
            "observed_accuracy": self.observed_accuracy,
            "gap": self.gap,
            "meaningful": self.meaningful,
            "min_n_for_meaning": self.min_n,
        }


def reliability_bins(
    y_true: Sequence[str],
    probas: Sequence[Optional[Dict[str, float]]],
    n_bins: int = 10,
    min_n: int = 10,
) -> List[ReliabilityBin]:
    """Confidence-calibration bins on the probability of the predicted class.

    For each example the predicted class is ``argmax p``, its confidence is
    ``max p``, and the bin's observed accuracy is the fraction of its examples
    whose predicted class was right. Equal-width bins on [0, 1]; the last bin
    is closed on the right so confidence exactly 1.0 lands in it.
    """
    if n_bins < 1:
        raise ValueError("n_bins must be >= 1")
    buckets: List[List[Tuple[float, bool]]] = [[] for _ in range(n_bins)]
    for true, proba in zip(y_true, probas):
        if not proba:
            continue
        # argmax with a deterministic tie-break (highest probability, then the
        # LABELS order) so a tied vector does not bin differently run to run.
        best_confidence, predicted = max(
            (proba.get(label, 0.0), -index, label)
            for index, label in enumerate(LABELS)
        )[0::2]
        confidence = float(best_confidence)
        index = min(int(confidence * n_bins), n_bins - 1)
        buckets[index].append((confidence, predicted == true))
    out: List[ReliabilityBin] = []
    for index, bucket in enumerate(buckets):
        low = index / n_bins
        high = (index + 1) / n_bins
        if not bucket:
            out.append(ReliabilityBin(low, high, 0, None, None, min_n))
            continue
        mean_confidence = sum(c for c, _ in bucket) / len(bucket)
        observed = sum(1 for _, correct in bucket if correct) / len(bucket)
        out.append(
            ReliabilityBin(low, high, len(bucket), mean_confidence, observed, min_n)
        )
    return out


def expected_calibration_error(bins: Sequence[ReliabilityBin]) -> Optional[float]:
    """Sample-weighted mean |accuracy - confidence| over non-empty bins."""
    total = sum(b.n for b in bins)
    if total == 0:
        return None
    error = 0.0
    for b in bins:
        if b.n == 0 or b.gap is None:
            continue
        error += (b.n / total) * abs(b.gap)
    return error


def scores_are_constant(
    probas: Sequence[Optional[Dict[str, float]]], label: str, tol: float = 1e-12
) -> bool:
    """True when every example got the same score for ``label``.

    A constant score cannot rank anything, so an AUC computed from it is 0.5 by
    arithmetic rather than by measurement. :class:`GlobalMajority` is exactly
    this case, and printing 0.5 for it would suggest a chance-level ranker
    where there is in fact no ranker.
    """
    values = [p.get(label, 0.0) for p in probas if p]
    if len(values) < 2:
        return True
    first = values[0]
    return all(abs(v - first) <= tol for v in values)


def roc_auc_ovr(
    y_true: Sequence[str],
    probas: Sequence[Optional[Dict[str, float]]],
    label: str,
) -> Optional[float]:
    """One-vs-rest ROC AUC for one class, by the rank (Mann-Whitney) identity.

    Ties get average ranks, so a partially-tied score set gets the correct
    0.5-credit rather than an order-dependent answer. ``None`` when the class
    has no positives or no negatives in this fold, or when the scores are
    constant: in each case the number would be an artefact.
    """
    pairs = [
        (float(p.get(label, 0.0)), t == label)
        for t, p in zip(y_true, probas)
        if p
    ]
    if not pairs:
        return None
    n_pos = sum(1 for _, positive in pairs if positive)
    n_neg = len(pairs) - n_pos
    if n_pos == 0 or n_neg == 0:
        return None
    if scores_are_constant([p for p in probas if p], label):
        return None
    ordered = sorted(range(len(pairs)), key=lambda i: pairs[i][0])
    ranks = [0.0] * len(pairs)
    i = 0
    while i < len(ordered):
        j = i
        while (
            j + 1 < len(ordered)
            and pairs[ordered[j + 1]][0] == pairs[ordered[i]][0]
        ):
            j += 1
        average = (i + j + 2) / 2.0  # 1-based average rank of the tie block
        for k in range(i, j + 1):
            ranks[ordered[k]] = average
        i = j + 1
    rank_sum_pos = sum(r for r, (_, positive) in zip(ranks, pairs) if positive)
    return (rank_sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def roc_auc_ovr_macro(
    y_true: Sequence[str],
    probas: Sequence[Optional[Dict[str, float]]],
    labels: Sequence[str] = LABELS,
) -> Tuple[Optional[float], Dict[str, Optional[float]]]:
    """(macro AUC over the classes where it is defined, per-class AUCs)."""
    per_class = {label: roc_auc_ovr(y_true, probas, label) for label in labels}
    defined = [v for v in per_class.values() if v is not None]
    if not defined:
        return None, per_class
    return sum(defined) / len(defined), per_class


# --------------------------------------------------------------------------- #
# bootstrap over groups
# --------------------------------------------------------------------------- #

def percentile(values: Sequence[float], q: float) -> Optional[float]:
    """Linear-interpolation percentile, ``q`` in [0, 1]. Documented so that an
    interval printed by this project can be reproduced without guessing which
    of the nine conventional definitions was used."""
    if not values:
        return None
    if not 0.0 <= q <= 1.0:
        raise ValueError("q must be in [0, 1]")
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = q * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def bootstrap_group_ci(
    metric_fn: Callable[[Sequence[int]], Optional[float]],
    groups: Sequence[Sequence[int]],
    n_boot: int = 1000,
    alpha: float = 0.05,
    seed: int = 20240101,
) -> Dict[str, object]:
    """Percentile bootstrap confidence interval, **resampling whole groups**.

    ``groups`` is a partition of row indices into independent blocks -- one
    study, one lab, one population, depending on the grouping. Each bootstrap
    replicate draws ``len(groups)`` groups *with replacement* and concatenates
    their rows, so a group drawn twice contributes all of its rows twice and a
    group not drawn contributes none.

    This is the cluster bootstrap, and the difference from resampling rows is
    the entire point. ``docs/SCHEMA.md`` 5 names lab identity as the largest
    source of non-independence in this corpus: one lab contributing forty rows
    about one circuit is forty correlated observations, not forty independent
    ones. A row-level bootstrap treats them as independent and returns an
    interval that is too narrow -- by roughly the square root of the average
    group size -- which makes a baseline look more precisely measured than the
    evidence allows. Resampling groups propagates the uncertainty of *which
    labs happened to study this*, which is the uncertainty a reader of the
    benchmark actually cares about.

    The returned dict is intentionally not a metric dict: it holds an interval
    and the bookkeeping to interpret it, and it is embedded by a
    grouping-labelled container in ``syn4rch.models.evaluate`` rather than
    serialised on its own.
    """
    if n_boot < 1:
        raise ValueError("n_boot must be >= 1")
    blocks = [list(g) for g in groups if len(g) > 0]
    n_groups = len(blocks)
    n_rows = sum(len(b) for b in blocks)
    result: Dict[str, object] = {
        "method": "percentile cluster bootstrap, resampled over groups",
        "resampling_unit": "group",
        "n_groups": n_groups,
        "n_rows": n_rows,
        "n_bootstrap_requested": n_boot,
        "n_bootstrap_usable": 0,
        "n_bootstrap_undefined": 0,
        "alpha": alpha,
        "seed": seed,
        "ci_low": None,
        "ci_high": None,
        "bootstrap_median": None,
    }
    if n_groups == 0:
        result["note"] = "no groups, so no interval"
        return result
    if n_groups == 1:
        result["note"] = (
            "a single group cannot be resampled: every replicate is the same "
            "data, so no interval is reported"
        )
        return result

    rng = random.Random(seed)
    values: List[float] = []
    undefined = 0
    for _ in range(n_boot):
        indices: List[int] = []
        for _ in range(n_groups):
            indices.extend(blocks[rng.randrange(n_groups)])
        value = metric_fn(indices)
        if value is None:
            undefined += 1
        else:
            values.append(float(value))
    result["n_bootstrap_usable"] = len(values)
    result["n_bootstrap_undefined"] = undefined
    if not values:
        result["note"] = (
            "the metric was undefined in every replicate, so no interval "
            "is reported"
        )
        return result
    result["ci_low"] = percentile(values, alpha / 2.0)
    result["ci_high"] = percentile(values, 1.0 - alpha / 2.0)
    result["bootstrap_median"] = percentile(values, 0.5)
    # Meaningful only when the metric is a *difference*: the share of
    # replicates in which the candidate did not gain. Reported for every
    # metric because computing it is free and because for a difference it is
    # the number a reader wants beside the interval; for a plain accuracy it is
    # the share of replicates at or below zero, which is ordinarily 0 and
    # harmless.
    result["fraction_at_or_below_zero"] = (
        sum(1 for v in values if v <= 0.0) / len(values)
    )
    return result


def mcnemar_exact_p(b: int, c: int) -> Optional[float]:
    """Two-sided exact McNemar p-value from the two discordant counts.

    ``b`` and ``c`` are the rows where exactly one of the two predictors was
    right. Under the null that neither is better, each discordant row is a fair
    coin, so the p-value is the two-sided binomial tail of ``min(b, c)``
    successes in ``b + c`` trials. Exact rather than the chi-square
    approximation because the counts here are tens, not thousands.

    ``None`` when there are no discordant rows: two predictors that agree
    everywhere give no evidence either way, which is not the same as a
    p-value of 1 for a test that was actually run.

    This is a **row-level** test and it treats rows as independent, which in
    this corpus they are not (SCHEMA.md 5). It is reported *beside* the
    group-clustered bootstrap, never instead of it: when the two disagree, the
    honest reading is that the effect is suggestive and not established, and
    the disagreement is itself the finding.
    """
    n = b + c
    if n <= 0:
        return None
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / (2.0 ** n)
    return min(1.0, 2.0 * tail)


def index_groups(group_values: Sequence[Optional[str]]) -> List[List[int]]:
    """Row indices blocked by group key. Rows with key ``None`` are omitted.

    Omitting them is deliberate and must be reported by the caller: a record
    whose lab is unknown (DECISIONS.md D-009) cannot be grouped, and pooling
    every such record into one anonymous "unknown lab" group would invent a
    collaboration between unrelated papers.
    """
    buckets: Dict[str, List[int]] = {}
    for index, value in enumerate(group_values):
        if value is None:
            continue
        buckets.setdefault(value, []).append(index)
    return [buckets[key] for key in sorted(buckets)]
