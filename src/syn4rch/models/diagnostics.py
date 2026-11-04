"""Diagnostics that ask whether a feature can carry signal at all.

A benchmark table answers "how well did this baseline do". It does not answer
"was there anything there to do well on", and the two are different questions
whose answers are easy to confuse. This module computes the second kind, from
the dataset alone, before any model is fitted:

* :func:`molecular_axis_mirroring` -- the falsification test for the schema
  v0.2.0 molecular axis. If "loss of a gene product in this population"
  predicts the sign of a readout, then a loss-of-function record and a
  gain-of-function record in the same ``(population, outcome)`` cell should
  tend to point *opposite ways*. Where they do not, the axis is not carrying a
  sign in this corpus, and the honest thing is a counted statement of that,
  not a search for a grouping under which the number looks better.

* :func:`perturbation_identity_novelty` -- SCHEMA.md 5b: a cross-validated
  number must be reported alongside the fraction of test rows whose
  perturbation identity never appeared in training. Grouping stops a study
  leaking across a split; it does not establish that anything on the test side
  was learnable from the train side.

* :func:`cell_structure` -- what each cell definition does to the data. A
  finer cell is a smaller cell. In-sample purity rises mechanically as cells
  shrink, so purity is reported *next to* the singleton rate, which is the
  thing that makes a pure cell worthless out of fold.

Nothing here fits a model, resamples, or produces a score. Everything is a
count, and a count of what was in the table.

Python 3.10 compatible, standard library only.
"""

from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .task import (
    CONDITION_FEATURE_FIELDS,
    LABELS,
    NULL_FEATURE,
    PERTURBATION_AXIS,
    Example,
    TaskDataset,
)

__all__ = [
    "LOSS_OF_FUNCTION",
    "GAIN_OF_FUNCTION",
    "MOLECULAR",
    "CELL_DEFINITIONS",
    "modal_label",
    "molecular_axis_mirroring",
    "perturbation_identity_novelty",
    "cell_structure",
    "dataset_diagnostics",
]

#: ``perturbation_polarity`` of a molecular loss of function (knockout).
LOSS_OF_FUNCTION = "decrease"

#: ``perturbation_polarity`` of a molecular gain of function (overexpression).
GAIN_OF_FUNCTION = "increase"

#: The ``perturbation_class`` whose axis is the amount of a gene product.
MOLECULAR = "molecular"

#: Cell definitions compared by :func:`cell_structure`, in the order a report
#: should read them: the v0.1.0 cell, then each v0.2.0 refinement of it.
CELL_DEFINITIONS: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("population_outcome", ("population_key", "outcome_id")),
    (
        "population_axis_outcome",
        ("population_key", PERTURBATION_AXIS, "outcome_id"),
    ),
    (
        "population_axis_outcome_condition",
        ("population_key", PERTURBATION_AXIS, "outcome_id")
        + CONDITION_FEATURE_FIELDS,
    ),
)

#: The two labels that assert a direction of change. ``no_change`` is a real,
#: recorded result (SCHEMA.md 3.5 forbids dropping it) but it is not a sign, so
#: a pair involving it can neither mirror nor fail to mirror.
_SIGNED_LABELS = ("increase", "decrease")


def modal_label(counts: Mapping[str, int]) -> Optional[str]:
    """Modal label, ties broken on :data:`~syn4rch.models.task.LABELS` order.

    The same documented, deterministic, arbitrary rule the baselines use, for
    the same reason: a tie means the data does not choose, and resolving it by
    file order would make a diagnostic depend on how the file was written.
    """
    best: Optional[str] = None
    best_count = -1
    for label in LABELS:
        count = counts.get(label, 0)
        if count > best_count:
            best, best_count = label, count
    return None if best_count <= 0 else best


def _counts(labels: Iterable[str]) -> Dict[str, int]:
    out = {label: 0 for label in LABELS}
    for label in labels:
        out[label] = out.get(label, 0) + 1
    return out


def _nonzero(counts: Mapping[str, int]) -> Dict[str, int]:
    return {label: counts[label] for label in LABELS if counts.get(label)}


def _mirroring(
    molecular: Sequence[Example],
    key_fields: Tuple[str, ...],
    name: str,
    description: str,
) -> Dict[str, Any]:
    """One mirroring analysis, at the cell granularity ``key_fields`` defines.

    Three per-cell verdicts:

    ``mirrored``
        the two sides' modal directions are opposite (one ``increase``, one
        ``decrease``) -- what the axis predicts;
    ``same``
        the two sides' modal directions are the same assertive direction --
        what the axis predicts against;
    ``no_change_involved``
        at least one side's mode is ``no_change``, so the pair makes no
        statement about a sign either way.

    Record pairs are counted the same way, because a cell with four records on
    one side and one on the other is weaker evidence than its modal comparison
    suggests.

    ``n_distinct_targets``, ``n_distinct_studies`` and ``n_distinct_labs`` per
    cell are what decides whether a verdict means anything. A cell whose two
    sides remove and add *different genes* was never a mirror pair; a cell
    whose records all come from one laboratory tests the polarity only within
    that laboratory's own report, which is also how a coherent narrative looks.
    """
    cells: Dict[Tuple[str, ...], Dict[str, List[Example]]] = {}
    for example in molecular:
        polarity = example.features.perturbation_polarity or NULL_FEATURE
        if polarity not in (LOSS_OF_FUNCTION, GAIN_OF_FUNCTION):
            continue
        key = example.features.cell(*key_fields)
        side = "loss_of_function" if polarity == LOSS_OF_FUNCTION else "gain_of_function"
        cells.setdefault(key, {}).setdefault(side, []).append(example)

    both = {
        key: sides
        for key, sides in cells.items()
        if sides.get("loss_of_function") and sides.get("gain_of_function")
    }

    rows: List[Dict[str, Any]] = []
    verdicts = {"mirrored": 0, "same": 0, "no_change_involved": 0}
    pairs = {"opposite_sign": 0, "same_sign": 0, "no_change_involved": 0}
    by_provenance: Dict[str, Dict[str, int]] = {
        "single_study": {"mirrored": 0, "same": 0, "no_change_involved": 0},
        "multiple_studies": {"mirrored": 0, "same": 0, "no_change_involved": 0},
    }
    n_multi_lab = 0
    for key, sides in sorted(both.items()):
        loss = sides["loss_of_function"]
        gain = sides["gain_of_function"]
        loss_counts = _counts(e.label for e in loss)
        gain_counts = _counts(e.label for e in gain)
        loss_mode = modal_label(loss_counts)
        gain_mode = modal_label(gain_counts)
        if loss_mode in _SIGNED_LABELS and gain_mode in _SIGNED_LABELS:
            verdict = "mirrored" if loss_mode != gain_mode else "same"
        else:
            verdict = "no_change_involved"
        verdicts[verdict] += 1
        for left in loss:
            for right in gain:
                if left.label in _SIGNED_LABELS and right.label in _SIGNED_LABELS:
                    pairs[
                        "opposite_sign" if left.label != right.label else "same_sign"
                    ] += 1
                else:
                    pairs["no_change_involved"] += 1
        targets = {
            e.features.perturbation_target
            for e in loss + gain
            if e.features.perturbation_target
        }
        studies = {e.group_keys.get("study_id") for e in loss + gain}
        labs = {
            e.group_keys.get("lab")
            for e in loss + gain
            if e.group_keys.get("lab") is not None
        }
        if len(labs) > 1:
            n_multi_lab += 1
        provenance = "single_study" if len(studies) == 1 else "multiple_studies"
        by_provenance[provenance][verdict] += 1
        rows.append(
            {
                "cell": list(key),
                "cell_provenance": provenance,
                "loss_of_function": _nonzero(loss_counts),
                "gain_of_function": _nonzero(gain_counts),
                "loss_modal_direction": loss_mode,
                "gain_modal_direction": gain_mode,
                "verdict": verdict,
                "n_records": len(loss) + len(gain),
                "n_distinct_targets": len(targets),
                "n_distinct_studies": len(studies),
                "n_distinct_labs": len(labs),
            }
        )

    signed_pairs = pairs["opposite_sign"] + pairs["same_sign"]
    return {
        "level": name,
        "cell_fields": list(key_fields),
        "description": description,
        "n_cells": len(cells),
        "n_cells_with_both_polarities": len(both),
        "n_cells_with_both_polarities_and_two_or_more_labs": n_multi_lab,
        "cell_verdicts": verdicts,
        # The same verdicts split by whether the cell's records all come from
        # one study. At gene level a single-study cell is one laboratory
        # manipulating one gene both ways -- the cleanest test available here,
        # and also the one most exposed to a paper's internal coherence. A
        # multi-study cell at (population, outcome) level pools genes, so a
        # verdict there is a statement about the pooling, not the polarity.
        "cell_verdicts_by_provenance": by_provenance,
        "record_pairs": pairs,
        "fraction_of_signed_pairs_mirrored": (
            None if signed_pairs == 0 else pairs["opposite_sign"] / signed_pairs
        ),
        "cells": rows,
    }


def molecular_axis_mirroring(dataset: TaskDataset) -> Dict[str, Any]:
    """Does the molecular polarity predict the sign of the readout?

    The axis says: loss of a gene product and gain of the same gene product
    should move a readout opposite ways. The test is run at two granularities,
    because they answer different questions and only one of them is a test of
    the polarity at all:

    ``by_population_outcome``
        the cell as it was before ``target_raw`` existed: every molecular
        record sharing a population and a readout, whatever gene it concerns.
        Reported because it is what a model keyed on ``(population, class,
        polarity, outcome)`` actually sees when the target is unknown, and
        because the difference between the two levels *is* the measurement of
        what pooling genes does.

    ``by_population_outcome_target``
        the honest test: the same gene, in the same population, on the same
        readout, manipulated both ways.

    ``n_cells_with_both_polarities_and_two_or_more_labs`` is the number that
    decides whether any of this can support prediction. A mirrored pair inside
    one paper is consistent with the polarity carrying a sign and equally
    consistent with a paper reporting a coherent story; only an independent
    laboratory reproducing the pair distinguishes them, and only a cell seen in
    training can transfer to a held-out study at all (SCHEMA.md 5b).
    """
    molecular = [
        e
        for e in dataset.examples
        if (e.features.perturbation_class or "") == MOLECULAR
    ]
    n_known_target = sum(1 for e in molecular if e.features.perturbation_target)
    return {
        "n_molecular_examples": len(molecular),
        "n_examples_with_known_target": n_known_target,
        "by_population_outcome": _mirroring(
            molecular,
            ("population_key", "outcome_id"),
            "by_population_outcome",
            "gene-pooled: all molecular records sharing a population and a "
            "readout, whatever gene was manipulated",
        ),
        "by_population_outcome_target": _mirroring(
            molecular,
            ("population_key", "outcome_id", "perturbation_target"),
            "by_population_outcome_target",
            "gene-specific: the same gene in the same population on the same "
            "readout, manipulated both ways",
        ),
        "what_would_falsify_this": (
            "the axis predicts that loss and gain of function on one gene "
            "point opposite ways. Mirrored cells and opposite-signed record "
            "pairs are evidence for it; 'same' cells and same-signed pairs are "
            "evidence against; pairs involving no_change are evidence of "
            "neither. At the gene-pooled level a cell with several distinct "
            "targets tests nothing about the polarity, because its two sides "
            "concern different genes -- and that, rather than the polarity "
            "being wrong, is what a negative result there means."
        ),
    }


def _identity(example: Example) -> Tuple[str, ...]:
    """The perturbation identity of one example: its axis cell key.

    ``(class, polarity, target)`` from
    :func:`syn4rch.qc.keys.perturbation_key`, via
    :meth:`syn4rch.models.task.Features.cell`, so the identity used here is the
    identity the baselines key on and not a second notion of it.
    """
    return example.features.cell(PERTURBATION_AXIS)


def perturbation_identity_novelty(
    dataset: TaskDataset,
    folds: Sequence[Tuple[Sequence[int], Sequence[int]]],
) -> Dict[str, Any]:
    """Fraction of test rows whose perturbation identity was unseen in training.

    Required beside every cross-validated number by ``docs/SCHEMA.md`` 5b: for
    13 of the molecular target genes in this corpus the gene appears in exactly
    one study and one laboratory, so under ``by_study`` or ``by_lab`` grouping a
    held-out molecular record's perturbation never appeared in training, and an
    accuracy averaged over such folds measures the label prior rather than
    generalisation.

    ``identity_is_coarse`` is the caveat that decides how to read the number,
    and it is about **molecular** records only. A molecular perturbation whose
    gene is not recorded has an incomplete identity: it degenerates to
    ``(molecular, polarity)``, which is seen in training almost always, so the
    unseen fraction is then a lower bound and the real novelty is higher. A
    circuit, structural or pharmacological manipulation has no gene to name, so
    ``(class, polarity)`` *is* its identity at the granularity this schema
    records, and its contribution to the fraction is not an underestimate.
    """
    n_test = 0
    n_unseen = 0
    per_fold: List[Dict[str, Any]] = []
    for train_idx, test_idx in folds:
        seen = {_identity(dataset.examples[i]) for i in train_idx}
        unseen = sum(
            1 for i in test_idx if _identity(dataset.examples[i]) not in seen
        )
        n_test += len(test_idx)
        n_unseen += unseen
        per_fold.append(
            {
                "n_test_rows": len(test_idx),
                "n_unseen_identity": unseen,
                "fraction_unseen": (
                    None if not test_idx else unseen / len(test_idx)
                ),
            }
        )
    n_known_target = sum(
        1 for e in dataset.examples if e.features.perturbation_target
    )
    molecular = [
        e
        for e in dataset.examples
        if (e.features.perturbation_class or "") == MOLECULAR
    ]
    n_molecular_without_target = sum(
        1 for e in molecular if not e.features.perturbation_target
    )
    return {
        "identity": "(perturbation_class, perturbation_polarity, perturbation_target)",
        "identity_source": "syn4rch.qc.keys.perturbation_key",
        "n_test_rows": n_test,
        "n_test_rows_with_unseen_identity": n_unseen,
        "fraction_unseen": (None if n_test == 0 else n_unseen / n_test),
        "n_examples_with_known_target": n_known_target,
        "n_molecular_examples": len(molecular),
        "n_molecular_without_target": n_molecular_without_target,
        "identity_is_coarse": n_molecular_without_target > 0,
        "note": (
            "SCHEMA.md 5b. %d of %d molecular examples name the gene that was "
            "manipulated.%s A non-molecular perturbation has no gene to name, "
            "so (class, polarity) is its identity at the granularity this "
            "schema records."
            % (
                len(molecular) - n_molecular_without_target,
                len(molecular),
                ""
                if n_molecular_without_target == 0
                else " For the %d that do not, the identity degenerates to "
                "(molecular, polarity) and this fraction is a lower bound on "
                "how novel the held-out perturbations really were."
                % n_molecular_without_target,
            )
        ),
        "per_fold": per_fold,
    }


def cell_structure(dataset: TaskDataset) -> Dict[str, Any]:
    """What each cell definition does to the data: size, singletons, purity.

    ``in_sample_purity`` is the share of rows that agree with their own cell's
    modal label. It is the **ceiling** an exact-cell lookup could reach if it
    were allowed to memorise the data it is scored on, and it rises
    mechanically as cells shrink: a cell of one row is always 100% pure. It is
    therefore reported next to ``fraction_singleton_cells``, which is the share
    of cells that hold a single row and so can never contribute a prediction to
    a held-out row of the same cell. A refinement whose purity rises while its
    singleton share rises faster has fragmented the data, not explained it.
    """
    rows: List[Dict[str, Any]] = []
    n_examples = len(dataset.examples)
    for name, fields_ in CELL_DEFINITIONS:
        cells: Dict[Tuple[str, ...], Dict[str, int]] = {}
        for example in dataset.examples:
            key = example.features.cell(*fields_)
            cell = cells.setdefault(key, {label: 0 for label in LABELS})
            cell[example.label] = cell.get(example.label, 0) + 1
        sizes = [sum(counts.values()) for counts in cells.values()]
        modal_total = sum(
            max(counts.get(label, 0) for label in LABELS) for counts in cells.values()
        )
        singletons = sum(1 for size in sizes if size == 1)
        rows.append(
            {
                "cell": name,
                "fields": list(fields_),
                "n_cells": len(cells),
                "n_rows": n_examples,
                "mean_rows_per_cell": (
                    None if not cells else n_examples / len(cells)
                ),
                "largest_cell_rows": max(sizes) if sizes else 0,
                "n_singleton_cells": singletons,
                "fraction_singleton_cells": (
                    None if not cells else singletons / len(cells)
                ),
                "n_rows_in_singleton_cells": singletons,
                "in_sample_purity": (
                    None if n_examples == 0 else modal_total / n_examples
                ),
            }
        )
    return {
        "definitions": rows,
        "note": (
            "in_sample_purity is a ceiling measured on the rows it was "
            "computed from, not an accuracy; read it against "
            "fraction_singleton_cells, because a cell holding one row is "
            "trivially pure and contributes nothing to a held-out row."
        ),
    }


def dataset_diagnostics(dataset: TaskDataset) -> Dict[str, Any]:
    """Every dataset-level diagnostic, for the benchmark result and the report."""
    return {
        "molecular_axis_mirroring": molecular_axis_mirroring(dataset),
        "cell_structure": cell_structure(dataset),
    }
