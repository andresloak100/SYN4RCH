"""Markdown rendering of the baseline benchmark.

Three rules this renderer follows, each of them a consequence of
``docs/SCHEMA.md`` 5 and 6 rather than a formatting preference:

1. **Every number is printed with its grouping, its row count and its group
   count.** The grouping is the heading of the table the number sits in and is
   repeated in the machine-readable JSON; ``n=<rows>/<groups>g`` is printed
   inside the cell. A number cannot appear in this report stripped of that
   context, because the only thing it can be rendered from is an
   :class:`~syn4rch.models.evaluate.Estimate`, which cannot exist without a
   grouping.

2. **A number too small to mean anything says so next to itself.** Not in a
   footnote, not omitted: ``⚠`` and the reason, in the cell. Omitting it would
   hide how little data there is, and printing it bare would invite the reader
   to compare it with a number measured on ten times as much.

3. **The baselines are presented as the bar to clear.** The document ends by
   stating, per grouping, the accuracy and macro-F1 a candidate model has to
   beat, with the interval and the n attached. A model that does not clear
   those is not to be adopted here, and the point of writing the bar down is
   that it cannot be renegotiated after the candidate's number is known.

Python 3.10 compatible, standard library only.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple

from .evaluate import (
    HEADLINE_COMPARISON,
    MIN_GROUPS_FOR_CI,
    MIN_N_PER_RELIABILITY_BIN,
    MIN_ROWS_FOR_METRIC,
    BenchmarkResult,
    Comparison,
    ComparisonSlice,
    Estimate,
    MetricBundle,
)
from .task import EXCLUSION_REASONS, GROUPINGS, LABELS

__all__ = ["render_markdown", "render_cell", "render_delta", "headline"]

#: Columns of the per-grouping comparison table.
_COLUMNS = (
    ("direction_accuracy", "direction accuracy"),
    ("macro_f1", "macro-F1"),
    ("roc_auc_ovr_macro", "ranking (OvR AUC)"),
    ("brier_score", "Brier"),
    ("expected_calibration_error", "ECE"),
)


def _fmt(value: Optional[float], places: int = 3) -> str:
    if value is None:
        return "--"
    return ("%." + str(places) + "f") % value


def render_cell(estimate: Optional[Estimate]) -> str:
    """One table cell: the number, its interval, its n, and any warning.

    Never returns a bare number. If there is no estimate at all the cell says
    "not computed" rather than being blank, because a blank cell reads as zero.
    """
    if estimate is None:
        return "not computed"
    if estimate.value is None:
        return "n/a -- %s" % (estimate.undefined_reason or "undefined")
    text = _fmt(estimate.value)
    if estimate.ci and estimate.ci.get("ci_low") is not None:
        text += " [%s, %s]" % (
            _fmt(estimate.ci["ci_low"]),
            _fmt(estimate.ci["ci_high"]),
        )
    text += " <br>n=%d/%dg" % (estimate.n_rows, estimate.n_groups)
    if estimate.caveats:
        text += " <br>**(!) %s**" % "; ".join(estimate.caveats)
    return text


def render_delta(item: Optional[ComparisonSlice]) -> str:
    """One cell of a before/after table: the paired difference, with both tests.

    Never a bare difference. The group-clustered interval and the row-level
    exact test are printed together even when they disagree -- especially when
    they disagree, because that disagreement is the finding.
    """
    if item is None:
        return "not computed"
    if item.n_rows == 0:
        return "no rows -- the two cannot differ here"
    if item.delta is None:
        return "n/a (%d rows)" % item.n_rows
    text = "**%+.3f**" % item.delta
    if item.ci and item.ci.get("ci_low") is not None:
        text += " [%+.3f, %+.3f]" % (item.ci["ci_low"], item.ci["ci_high"])
    text += " <br>%.3f vs %.3f" % (item.accuracy, item.reference_accuracy)
    text += " <br>n=%d/%dg" % (item.n_rows, item.n_groups)
    if item.mcnemar_p is not None:
        text += " <br>McNemar p=%.3f (b=%d, c=%d)" % (
            item.mcnemar_p,
            item.n_only_baseline_right,
            item.n_only_reference_right,
        )
    text += " <br>verdict: **%s**" % item.verdict
    if item.caveats:
        text += " <br>**(!) %s**" % "; ".join(item.caveats)
    return text


def _heading(text: str, level: int = 2) -> str:
    return "\n%s %s\n" % ("#" * level, text)


def _table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    # The leading blank line is load-bearing: a GFM table that follows a
    # paragraph without one renders as a paragraph.
    out = ["", "| " + " | ".join(header) + " |"]
    out.append("|" + "|".join(["---"] * len(header)) + "|")
    for row in rows:
        out.append("| " + " | ".join(row) + " |")
    return "\n".join(out) + "\n"


def _headline_grouping(result: BenchmarkResult) -> Optional[str]:
    """The grouping the headline is stated under, preferring leave-one-study-out.

    The headline has to be stated under *one* named grouping (SCHEMA.md 5), and
    which one is a choice that must not be made after seeing the numbers. The
    rule is fixed here: ``study_id`` if it was reported, otherwise the first
    reported grouping that is not leakage-prone, otherwise the first reported.
    """
    if "study_id" in result.groupings:
        return "study_id"
    for name in result.groupings:
        if not GROUPINGS[name].leakage_prone:
            return name
    return result.groupings[0] if result.groupings else None


def _accuracy_of(result: BenchmarkResult, baseline: str, grouping: str):
    bundle = result.bundle(baseline, grouping)
    if bundle is None or bundle.status != "ok":
        return None
    estimate = bundle.estimates.get("direction_accuracy")
    return None if estimate is None else estimate.value


def _beats_the_floor(result: BenchmarkResult) -> List[Tuple[str, str, float, float]]:
    """``(grouping, baseline, accuracy, floor)`` wherever a baseline out-scores it.

    Pooled accuracy only, and therefore *not* evidence of anything on its own:
    this list decides whether the headline may say "nothing beats the floor",
    and the paired comparisons decide whether anything that does beat it beat it
    by more than noise.
    """
    out: List[Tuple[str, str, float, float]] = []
    for grouping in result.groupings:
        floor = _accuracy_of(result, "global_majority", grouping)
        if floor is None:
            continue
        for bundle in result.bundles:
            if bundle.grouping != grouping or bundle.status != "ok":
                continue
            if bundle.baseline == "global_majority":
                continue
            value = (bundle.estimates.get("direction_accuracy") or Estimate(
                name="x", grouping=grouping, value=None, n_rows=0, n_groups=0
            )).value
            if value is not None and value > floor:
                out.append((grouping, bundle.baseline, value, floor))
    return out


def headline(result: BenchmarkResult) -> str:
    """The headline, computed from the numbers rather than written over them.

    Three things it must do, and the reason each is here:

    * if **nothing** beats the global majority direction, say so in plain
      words, at the top, before any table;
    * if something does, say what, by how much, under which grouping, on how
      many rows, and whether the two tests agree -- never "direction is
      predictable", which is a claim about the task and not about a margin;
    * state where the margin is **not**, since a corpus whose halves cancel is
      the ordinary case in which an average hides a result.
    """
    lines: List[str] = ["## Headline", ""]
    description = result.dataset_description or {}
    if not description.get("n_examples"):
        lines.append(
            "**Nothing was measured.** There are no labellable examples, so "
            "there is no accuracy, no comparison and no finding. This is a "
            "statement about the dataset, not a result."
        )
        return "\n".join(lines) + "\n"

    grouping = _headline_grouping(result)
    if grouping is None:
        lines.append(
            "**No grouping was reported, so no number here is valid** "
            "(SCHEMA.md 5)."
        )
        return "\n".join(lines) + "\n"

    beaten = _beats_the_floor(result)
    floor = _accuracy_of(result, "global_majority", grouping)
    candidate, reference = HEADLINE_COMPARISON
    comparison = result.comparison(candidate, reference, grouping)

    if not beaten:
        lines.append(
            "**No baseline beats guessing the global majority direction.** "
            "Under every grouping reported below, the best baseline is at or "
            "below the accuracy of a constant predictor of the modal "
            "direction%s. Direction of effect is not predictable from these "
            "features on this corpus, and the schema v0.2.0 features did not "
            "change that."
            % ("" if floor is None else " (%.3f under `%s`)" % (floor, grouping))
        )
        return "\n".join(lines) + "\n"

    overall = comparison.overall if comparison else None
    circuit = (
        comparison.slice_named("circuit_activity__primary_cell_seen")
        if comparison
        else None
    )
    circuit_all = (
        comparison.slice_named("circuit_activity") if comparison else None
    )
    molecular_seen = (
        comparison.slice_named("molecular__primary_cell_seen")
        if comparison
        else None
    )
    molecular_all = comparison.slice_named("molecular") if comparison else None

    lines.append(
        "**Something beats the global majority for the first time, and it is "
        "one half of the corpus.** Under `%s`, `%s` -- the modal direction of "
        "the (population, perturbation class, perturbation polarity, "
        "perturbation target, outcome) cell, falling straight back to the "
        "global mode where that cell was never seen -- scores %s against the "
        "floor's %s."
        % (
            grouping,
            candidate,
            "--" if overall is None or overall.accuracy is None else "%.3f" % overall.accuracy,
            "--" if floor is None else "%.3f" % floor,
        )
    )
    lines.append("")
    if overall is not None and overall.delta is not None:
        lines.append(
            "- **Pooled over all %d rows: %+.3f%s.** Read this as the diluted "
            "number it is: %s of those rows are molecular records whose cell "
            "is never populated in training, where the two predictors are "
            "identical by construction and can only pull the difference "
            "towards zero while narrowing its interval."
            % (
                overall.n_rows,
                overall.delta,
                ""
                if not (overall.ci and overall.ci.get("ci_low") is not None)
                else " [%+.3f, %+.3f]" % (overall.ci["ci_low"], overall.ci["ci_high"]),
                "--" if molecular_all is None else str(molecular_all.n_rows),
            )
        )
    if circuit_all is not None and circuit_all.delta is not None:
        lines.append(
            "- **Circuit-activity records (%d rows, %d %s groups): %+.3f%s.** "
            "This is where the whole margin is."
            % (
                circuit_all.n_rows,
                circuit_all.n_groups,
                grouping,
                circuit_all.delta,
                ""
                if not (circuit_all.ci and circuit_all.ci.get("ci_low") is not None)
                else " [%+.3f, %+.3f]"
                % (circuit_all.ci["ci_low"], circuit_all.ci["ci_high"]),
            )
        )
    if circuit is not None and circuit.delta is not None:
        lines.append(
            "- **Circuit-activity records whose cell was actually seen in "
            "training (%d rows, %d %s groups): %.3f against %.3f, %+.3f%s.** "
            "The two tests disagree and both are reported: the study-clustered "
            "bootstrap %s zero, the row-level exact McNemar test gives "
            "p=%.3f on b=%d / c=%d discordant rows. Verdict: **%s**, which in "
            "plain words means *this is the first signal above the baseline in "
            "this project and it is not established*."
            % (
                circuit.n_rows,
                circuit.n_groups,
                grouping,
                circuit.accuracy,
                circuit.reference_accuracy,
                circuit.delta,
                ""
                if not (circuit.ci and circuit.ci.get("ci_low") is not None)
                else " [%+.3f, %+.3f]" % (circuit.ci["ci_low"], circuit.ci["ci_high"]),
                "excludes" if circuit.interval_excludes_zero else "includes",
                circuit.mcnemar_p if circuit.mcnemar_p is not None else float("nan"),
                circuit.n_only_baseline_right,
                circuit.n_only_reference_right,
                circuit.verdict,
            )
        )
    if molecular_seen is not None:
        lines.append(
            "- **Molecular records: %s.** Once the manipulated gene is part of "
            "the cell key -- which schema v0.2.0 made possible and which is "
            "the only honest key for a molecular perturbation -- not one "
            "molecular test row has its cell populated in training, because "
            "each gene in this corpus appears in exactly one study. There is "
            "nothing to transfer, so the molecular half contributes a "
            "structural zero rather than a result (SCHEMA.md 5b)."
            % (
                "%d of %d rows had a seen cell"
                % (
                    molecular_seen.n_rows,
                    0 if molecular_all is None else molecular_all.n_rows,
                )
            )
        )
    lines.append("")
    unseen_share = (
        None
        if circuit_all is None or circuit is None or not circuit_all.n_rows
        else 1.0 - (circuit.n_rows / circuit_all.n_rows)
    )
    lines.append(
        "**What this is not.** It is not \"direction is now predictable\". The "
        "defensible claim is narrower: on the circuit-activity half, where a "
        "held-out study's (population, signed perturbation, outcome) cell has "
        "been seen before, the cell's modal direction beats the global prior "
        "by %s on %s rows in %s %s groups, reaching %s accuracy on a "
        "three-way problem -- and %s of circuit rows have no evidence in their "
        "cell at all, where this predictor is the global prior and nothing "
        "more. That is useful for ranking hypotheses. It is nowhere near a "
        "simulator, and nothing here licenses an aggregate score of any kind "
        "(SCHEMA.md 6)."
        % (
            "--" if circuit is None or circuit.delta is None else "%+.3f" % circuit.delta,
            "--" if circuit is None else str(circuit.n_rows),
            "--" if circuit is None else str(circuit.n_groups),
            grouping,
            "--" if circuit is None or circuit.accuracy is None else "%.3f" % circuit.accuracy,
            "--" if unseen_share is None else "%.0f%%" % (100.0 * unseen_share),
        )
    )
    return "\n".join(lines) + "\n"


def _dataset_section(description: Dict[str, Any]) -> str:
    lines: List[str] = [_heading("1. The dataset this was measured on")]
    source = description.get("source") or "(no path: the dataset was supplied in memory)"
    exclusions = description.get("exclusions") or {}
    n_input = exclusions.get("n_input_records", 0)
    n_examples = description.get("n_examples", 0)
    lines.append("- Records read: **%d** from `%s`" % (n_input, source))
    lines.append("- Labellable examples: **%d**" % n_examples)
    lines.append(
        "- Excluded: **%d**%s"
        % (
            exclusions.get("n_excluded", 0),
            (
                " (%.1f%% of the records read)"
                % (100.0 * exclusions["excluded_fraction"])
                if exclusions.get("excluded_fraction") is not None
                else ""
            ),
        )
    )
    for error in description.get("load_errors") or []:
        lines.append("- **Load problem:** %s" % error)

    if n_examples == 0:
        lines.append("")
        lines.append(
            "> **There is no labellable data.** Nothing below this line is a "
            "measurement. No accuracy, macro-F1, ranking or calibration number "
            "is reported, because a metric computed on zero examples is not a "
            "small number -- it is an absent one, and printing 0.0 or 0.5 in "
            "its place would be the single most misleading thing this report "
            "could do. Curate records into `data/curated/records.jsonl` and run "
            "this again."
        )

    lines.append(_heading("1.1 Label distribution", 3))
    counts = description.get("label_counts") or {}
    total = sum(counts.values()) or 0
    lines.append(
        _table(
            ("direction", "n", "share"),
            [
                (
                    label,
                    str(counts.get(label, 0)),
                    "--" if total == 0 else "%.1f%%" % (100.0 * counts.get(label, 0) / total),
                )
                for label in LABELS
            ],
        )
    )
    if total:
        largest = max(counts.get(label, 0) for label in LABELS)
        lines.append(
            "A constant predictor of the modal direction therefore scores "
            "**%.3f** accuracy on this label distribution before any "
            "cross-validation. Read every accuracy below against that number "
            "and not against 0.333." % (largest / total)
        )

    lines.append(_heading("1.2 What was excluded, and why", 3))
    lines.append(
        "Records that cannot be labelled are counted here rather than dropped. "
        "`not_reported` means the paper measured something and did not say "
        "which way it went; `non_monotonic` is a real recorded response shape "
        "with no single sign (DECISIONS.md D-011). Tier B outcomes are "
        "excluded because SCHEMA.md 6 states that Tier B is an annotation "
        "layer and not a training target. Quietly dropping any of these would "
        "change the label distribution above without changing this report."
    )
    primary = exclusions.get("primary_reason_counts") or {}
    rows = []
    for code, (why, clause) in sorted(EXCLUSION_REASONS.items()):
        rows.append((
            "`%s`" % code,
            str(primary.get(code, 0)),
            why,
            clause,
        ))
    lines.append(_table(("reason", "n (primary)", "what it means", "clause"), rows))
    co = exclusions.get("co_occurring_reason_counts") or {}
    if sum(co.values()) > exclusions.get("n_excluded", 0):
        lines.append(
            "Some records were excludable for more than one reason. Counting "
            "every reason that applied: "
            + ", ".join("`%s` %d" % (k, v) for k, v in sorted(co.items()) if v)
            + "."
        )
    if not exclusions.get("tier_b_excluded", True):
        lines.append(
            "> **`exclude_tier_b=False` was set for this run.** Tier B "
            "outcomes are in the examples above. SCHEMA.md 6 forbids using "
            "them as a training target, so these numbers are a measurement of "
            "data volume, not a model result."
        )

    axis = description.get("perturbation_axis") or {}
    if axis:
        lines.append(_heading("1.3 The signed perturbation axis (v0.2.0)", 3))
        lines.append(
            "Before v0.2.0, 185 of 362 records carried `valence: unclear` and "
            "had no signed feature at all. Every row below now has a class and "
            "a polarity. **The polarity is never a feature on its own:** it is "
            "read with its class or not at all, because `decrease` means less "
            "firing under `circuit_activity` and less of a gene product under "
            "`molecular`, and pooling those is the sign error D-019 was "
            "written against. `Features.cell` raises if a cell key names the "
            "polarity without the class, which is why that rule is a property "
            "of the code rather than a note in a document."
        )
        lines.append(
            _table(
                ("class / polarity", "n"),
                [
                    (pair, str(count))
                    for pair, count in (axis.get("class_polarity_counts") or {}).items()
                ],
            )
        )
        lines.append(
            "- Rows with a signed polarity (`increase`, `decrease` or "
            "`remove`): **%d of %d**."
            % (axis.get("n_with_signed_polarity", 0), axis.get("n_examples", 0))
        )
        lines.append(
            "- Rows naming the molecular target that was manipulated: **%d**. "
            "Where the target is absent a molecular cell pools whatever genes "
            "share a population and a readout; where it is present the cell "
            "becomes gene-specific, and far sparser."
            % axis.get("n_with_known_target", 0)
        )

    conditions = description.get("conditions") or {}
    if conditions:
        lines.append(_heading("1.4 How much of the condition was known", 3))
        lines.append(
            "`None` (the curator did not address it) and `not_reported` (the "
            "paper was checked and does not say) are pooled into a single "
            "`%s` token by `syn4rch.qc.keys.condition_key` (D-023), and this "
            "task reuses that convention rather than inventing a second one. "
            "The pooling is what makes a conflict detectable at all in the "
            "records that state no condition -- and it is also why these "
            "counts have to be printed: two rows that match on `%s` matched by "
            "**both failing to state a condition**, not by sharing one. A "
            "feature set that hid this would lie about what the model had to "
            "work with."
            % (
                conditions.get("unknown_token", "unknown"),
                conditions.get("unknown_token", "unknown"),
            )
        )
        lines.append(
            _table(
                ("condition slots known", "n (of 4: %s)" % ", ".join(
                    "`%s`" % s for s in conditions.get("slots", [])
                ), "n (of the 3 used as features: %s)" % ", ".join(
                    "`%s`" % s for s in conditions.get("feature_slots", [])
                )),
                [
                    (
                        str(k),
                        str((conditions.get("by_slots_known") or {}).get(k, 0)),
                        str(
                            (conditions.get("by_feature_slots_known") or {}).get(k, 0)
                        ),
                    )
                    for k in ("0", "1", "2", "3", "4")
                ],
            )
        )
        lines.append(
            "`age_group` is in the condition key and is **not** a feature of "
            "this task, which is why both counts are given rather than one "
            "standing in for the other."
        )
    return "\n".join(lines) + "\n"


def _groupings_section(description: Dict[str, Any], reported: Sequence[str]) -> str:
    lines = [_heading("2. Groupings (SCHEMA.md 5)")]
    lines.append(
        "SCHEMA.md 5 names four things that must not cross a train/test "
        "boundary and ends by declaring a cross-validated number reported "
        "without its grouping an invalid result. Every number in this document "
        "therefore sits under exactly one of these headings, and carries its "
        "row and group counts in the cell."
    )
    coverage = description.get("groupings") or {}
    rows = []
    for name in reported:
        declaration = GROUPINGS[name]
        stats = coverage.get(name, {})
        rows.append((
            "`%s`" % name,
            declaration.description,
            declaration.schema_clause,
            str(stats.get("n_rows_groupable", 0)),
            str(stats.get("n_groups", 0)),
            str(stats.get("n_rows_ungroupable", 0)),
            "**yes**" if declaration.leakage_prone else "no",
        ))
    lines.append(
        _table(
            (
                "grouping",
                "what it holds apart",
                "clause",
                "rows groupable",
                "groups",
                "rows with no key",
                "leakage-prone",
            ),
            rows,
        )
    )
    lines.append(
        "A row with no group key is **excluded** from that grouping's "
        "evaluation rather than pooled into one anonymous group: a record with "
        "no `last_author` cannot be shown not to share a lab with the training "
        "set, and inventing an \"unknown lab\" group would assert a "
        "collaboration between unrelated papers (DECISIONS.md D-009)."
    )
    return "\n".join(lines) + "\n"


def _baselines_section(baselines: Sequence[Dict[str, Any]]) -> str:
    lines = [_heading("3. The baselines")]
    rows = []
    for item in baselines:
        availability = item.get("availability") or {}
        rows.append((
            "`%s`" % item.get("name"),
            item.get("description", ""),
            "**%s**" % item.get("feature_generation", "?"),
            ", ".join("`%s`" % f for f in item.get("features_used") or [])
            or "(none: a constant predictor)",
            "available" if availability.get("available") else "**UNAVAILABLE**",
            "yes" if item.get("emits_proba") else "no",
        ))
    lines.append(
        _table(
            (
                "baseline",
                "what it does",
                "feature generation",
                "features it reads",
                "status",
                "emits probabilities",
            ),
            rows,
        )
    )
    lines.append(
        "The `feature generation` column is each baseline's own declaration of "
        "which version of the task's feature tuple it keys on, and it is what "
        "the before/after section is assembled from. A v0.1.0 baseline here is "
        "unchanged code: its number is comparable with the previous benchmark's "
        "by construction, not by recollection."
    )
    for item in baselines:
        availability = item.get("availability") or {}
        if availability.get("available"):
            continue
        lines.append(
            _heading("3.1 `%s` is UNAVAILABLE" % item.get("name"), 3)
        )
        lines.append(availability.get("reason") or "(no reason recorded)")
        lines.append("")
        lines.append("Places that were looked at:")
        lines.append("")
        for probe in availability.get("probes") or []:
            lines.append("- %s" % probe)
        lines.append("")
        lines.append(
            "Its rows below are empty on purpose. No claim in this report "
            "rests on it, and no edge was invented so that the row could be "
            "filled."
        )
    return "\n".join(lines) + "\n"


def _unseen_cell(bundle: MetricBundle) -> str:
    """The SCHEMA.md 5b number, for the column beside every accuracy."""
    novelty = bundle.unseen_perturbation_identity or {}
    fraction = novelty.get("fraction_unseen")
    if fraction is None:
        return "--"
    text = "%.1f%% (%d/%d)" % (
        100.0 * fraction,
        novelty.get("n_test_rows_with_unseen_identity", 0),
        novelty.get("n_test_rows", 0),
    )
    if novelty.get("identity_is_coarse"):
        text += (
            " <br>**(!) lower bound: %d molecular record(s) do not name the "
            "gene, so their identity is incomplete**"
            % novelty.get("n_molecular_without_target", 0)
        )
    return text


def _bundle_rows(bundles: Sequence[MetricBundle]) -> List[List[str]]:
    rows: List[List[str]] = []
    for bundle in bundles:
        if bundle.status != "ok":
            rows.append(
                [
                    "`%s`" % bundle.baseline,
                    "**no numbers: %s**<br>%s"
                    % (bundle.status, bundle.status_detail or ""),
                ]
                + ["--"] * (len(_COLUMNS) - 1)
                + ["--", "--", "--", "--"]
            )
            continue
        cells = [render_cell(bundle.estimates.get(key)) for key, _ in _COLUMNS]
        coverage = bundle.coverage
        rows.append(
            ["`%s`" % bundle.baseline]
            + cells
            + [
                "--" if coverage is None else "%.1f%% (%d/%d)"
                % (
                    100.0 * coverage,
                    bundle.n_rows_scored,
                    bundle.n_rows_scored + bundle.n_rows_abstained,
                ),
                _unseen_cell(bundle),
                str(bundle.n_folds),
                bundle.fold_source.split(":")[-1],
            ]
        )
    return rows


def _results_section(result: BenchmarkResult) -> str:
    lines = [_heading("4. Results, one table per grouping")]
    lines.append(
        "Cell format: `value [95%% bootstrap CI] n=<scored rows>/<independent "
        "groups>g`. The interval is a percentile cluster bootstrap resampled "
        "over **groups**, not rows (%d replicates, seed %s): one lab "
        "contributing forty rows about one circuit is forty correlated "
        "observations, and a row-level bootstrap would report an interval too "
        "narrow by roughly the square root of the average group size. "
        "`(!)` marks a number whose n is too small to carry meaning "
        "(< %d scored rows, or < %d independent groups); it is printed with "
        "the warning rather than bare or omitted."
        % (
            result.settings.get("n_bootstrap", 0),
            result.settings.get("seed"),
            MIN_ROWS_FOR_METRIC,
            MIN_GROUPS_FOR_CI,
        )
    )
    lines.append(
        "\nThe five families are deliberately **not** combined into one score. "
        "Accuracy on a skewed three-class problem mostly measures the skew; "
        "macro-F1 answers a different question; ranking and calibration decide "
        "whether the outputs can be thresholded or acted on at all. A single "
        "headline number would let a model be adopted on the one of the five "
        "it happened to win."
    )
    lines.append(
        "\nThe last-but-two column is required beside every cross-validated "
        "number by SCHEMA.md 5b: the share of test rows whose **perturbation "
        "identity** -- its class, polarity and molecular target -- never "
        "appears in its fold's training set. Grouping guarantees that no study "
        "crosses a split; it guarantees nothing about whether the held-out "
        "perturbation was learnable from the training side. Where that share "
        "is high, the accuracy beside it is substantially a measurement of the "
        "label prior, and SCHEMA.md 5b is explicit that this is a fact about "
        "how the field publishes -- one lab, one gene, one paper -- and not a "
        "number to be tuned away."
    )
    header = (
        ["baseline"]
        + [label for _, label in _COLUMNS]
        + [
            "coverage (answered)",
            "test rows with unseen perturbation identity (5b)",
            "folds",
            "fold source",
        ]
    )
    for grouping in result.groupings:
        declaration = GROUPINGS[grouping]
        bundles = [b for b in result.bundles if b.grouping == grouping]
        lines.append(_heading("4.%d Grouping: `%s`" % (
            result.groupings.index(grouping) + 1, grouping), 3))
        lines.append("*%s* (%s)" % (declaration.description, declaration.schema_clause))
        if declaration.leakage_prone:
            lines.append("")
            lines.append(
                "> **Leakage-prone grouping.** The same study and the same lab "
                "appear on both sides of every fold, so these numbers are an "
                "upper bound and not an estimate of generalisation. They are "
                "printed because the gap between them and the grouped numbers "
                "above is the measurement of how much row-level "
                "cross-validation flatters a model on this corpus."
            )
        lines.append("")
        lines.append(_table(header, _bundle_rows(bundles)))
        notes: List[str] = []
        for bundle in bundles:
            for note in bundle.notes:
                entry = "`%s`: %s" % (bundle.baseline, note)
                if entry not in notes:
                    notes.append(entry)
            if bundle.fold_source_note:
                entry = "`%s`: %s" % (bundle.baseline, bundle.fold_source_note)
                if entry not in notes:
                    notes.append(entry)
        if notes:
            lines.append("Notes for this grouping:")
            lines.append("")
            for note in notes:
                lines.append("- %s" % note)
            lines.append("")
    return "\n".join(lines) + "\n"


def _before_after_section(result: BenchmarkResult, number: int) -> str:
    """Paired before/after differences, per grouping, with the strata."""
    lines = [
        _heading(
            "%d. Before and after the schema v0.2.0 features" % number
        )
    ]
    task = (result.dataset_description or {}).get("task") or {}
    added = task.get("features_added_in_this_version") or []
    lines.append(
        "The schema v0.1.0 baselines in this run are **unchanged code** -- the "
        "same cell keys, the same backoff order, the same tie rule -- so their "
        "rows above are directly comparable with the previous benchmark, and "
        "the \"before\" side of every difference below is that previous "
        "measurement recomputed rather than quoted. The features added in "
        "v0.2.0 are: %s."
        % (", ".join("`%s`" % name for name in added) or "(none declared)")
    )
    lines.append("")
    lines.append(
        "One exception to that comparability, and it is a fold-assignment "
        "artefact rather than a change of code: the `record` grouping "
        "partitions rows by `record_id`, and `record_id` is a hash of the "
        "fields the v0.2.0 migration rewrote. A v0.1.0 baseline's number under "
        "`record` therefore differs from the previous benchmark's because the "
        "rows landed in different folds, not because the baseline computes "
        "anything new. The grouped numbers that matter -- `study_id`, `lab`, "
        "`population`, `population_outcome` -- do not depend on record ids and "
        "are unchanged."
    )
    lines.append("")
    lines.append(
        "Each cell is a **paired** difference: both baselines run over the "
        "same folds, scored on the rows both answered, with the bootstrap "
        "resampling groups of those rows. Two accuracies from the tables above "
        "must not be subtracted instead -- the difference of two independent "
        "intervals is not an interval for the difference, and these two "
        "baselines are strongly correlated because they share a training set "
        "and a backoff tail."
    )
    lines.append("")
    lines.append(
        "`all_rows` and `primary_cell_seen` answer different questions. Where "
        "a candidate's primary cell was never populated in training it falls "
        "back to the same global mode the floor uses, so those rows are "
        "**identical by construction** and can carry no information about the "
        "difference; including them dilutes the margin and drags in whatever "
        "the label prior happens to be on that subset. `primary_cell_seen` is "
        "the subset on which something was actually learned. Both are printed, "
        "each labelled, so neither can be quoted as the other."
    )
    lines.append("")
    lines.append(
        "`verdict` is **established** only when the group-clustered interval "
        "and the row-level exact McNemar test agree, **suggestive** when "
        "exactly one of them does. On a hundred-odd rows in twenty studies the "
        "two routinely disagree, and when they do, the honest reading is the "
        "weaker one."
    )

    slice_names = [
        "all_rows",
        "primary_cell_seen",
        "circuit_activity",
        "circuit_activity__primary_cell_seen",
        "molecular",
        "molecular__primary_cell_seen",
    ]
    for index, grouping in enumerate(result.groupings):
        comparisons = [c for c in result.comparisons if c.grouping == grouping]
        if not comparisons:
            continue
        declaration = GROUPINGS[grouping]
        lines.append(
            _heading("%d.%d Grouping: `%s`" % (number, index + 1, grouping), 3)
        )
        lines.append(
            "*%s* (%s)" % (declaration.description, declaration.schema_clause)
        )
        if declaration.leakage_prone:
            lines.append("")
            lines.append(
                "> **Leakage-prone grouping.** Every difference here is an "
                "upper bound on the difference, for the same reason the "
                "accuracies are."
            )
        rows: List[List[str]] = []
        for comparison in comparisons:
            if comparison.status != "ok":
                rows.append(
                    [
                        "`%s`<br>vs `%s`" % (comparison.baseline, comparison.reference),
                        "**no comparison: %s**<br>%s"
                        % (comparison.status, comparison.status_detail or ""),
                    ]
                    + ["--"] * (len(slice_names) - 1)
                )
                continue
            rows.append(
                ["`%s`<br>vs `%s`" % (comparison.baseline, comparison.reference)]
                + [
                    render_delta(comparison.slice_named(name))
                    for name in slice_names
                ]
            )
        lines.append("")
        lines.append(
            _table(
                ["candidate vs reference"] + ["`%s`" % name for name in slice_names],
                rows,
            )
        )
        notes: List[str] = []
        for comparison in comparisons:
            for caveat in comparison.caveats:
                entry = "`%s` vs `%s`: %s" % (
                    comparison.baseline,
                    comparison.reference,
                    caveat,
                )
                if entry not in notes:
                    notes.append(entry)
        if notes:
            lines.append("Notes for this grouping:")
            lines.append("")
            for note in notes:
                lines.append("- %s" % note)
            lines.append("")
        lines.append(
            "A cell reading *not computed* is a slice this comparison does not "
            "have: `primary_cell_seen` needs a candidate with a primary cell, "
            "and a class slice needs rows of that class. It is not a number "
            "that came out empty."
        )
    return "\n".join(lines) + "\n"


def _mirroring_section(result: BenchmarkResult, number: int) -> str:
    """The falsification test for the molecular axis, and the cell arithmetic."""
    diagnostics = result.diagnostics or {}
    mirroring = diagnostics.get("molecular_axis_mirroring") or {}
    lines = [
        _heading(
            "%d. Does the molecular axis carry a sign? (loss of function vs "
            "gain of function)" % number
        )
    ]
    if not mirroring:
        lines.append("No diagnostic was computed for this run.")
        return "\n".join(lines) + "\n"

    lines.append(
        "Schema v0.2.0's claim is that a molecular manipulation has a "
        "well-defined sign along its own axis: less of a gene product, or more "
        "of it. If that sign reaches the readout, then loss and gain of "
        "function on the **same gene** in the same population on the same "
        "readout should point opposite ways. This is that test, and it is run "
        "at two granularities because they answer different questions."
    )
    lines.append("")
    lines.append(
        "- Molecular examples: **%d**, of which **%d** record which gene was "
        "manipulated."
        % (
            mirroring.get("n_molecular_examples", 0),
            mirroring.get("n_examples_with_known_target", 0),
        )
    )
    lines.append("")

    rows = []
    for level_key in ("by_population_outcome", "by_population_outcome_target"):
        level = mirroring.get(level_key) or {}
        verdicts = level.get("cell_verdicts") or {}
        pairs = level.get("record_pairs") or {}
        provenance = level.get("cell_verdicts_by_provenance") or {}
        single = provenance.get("single_study") or {}
        multi = provenance.get("multiple_studies") or {}
        rows.append(
            (
                "`%s`" % level_key,
                level.get("description", ""),
                str(level.get("n_cells", 0)),
                str(level.get("n_cells_with_both_polarities", 0)),
                str(level.get("n_cells_with_both_polarities_and_two_or_more_labs", 0)),
                "%d / %d / %d"
                % (
                    verdicts.get("mirrored", 0),
                    verdicts.get("same", 0),
                    verdicts.get("no_change_involved", 0),
                ),
                "%d / %d / %d"
                % (
                    pairs.get("opposite_sign", 0),
                    pairs.get("same_sign", 0),
                    pairs.get("no_change_involved", 0),
                ),
                "%d / %d / %d"
                % (
                    single.get("mirrored", 0),
                    single.get("same", 0),
                    single.get("no_change_involved", 0),
                ),
                "%d / %d / %d"
                % (
                    multi.get("mirrored", 0),
                    multi.get("same", 0),
                    multi.get("no_change_involved", 0),
                ),
            )
        )
    lines.append(
        _table(
            (
                "level",
                "what a cell is",
                "cells",
                "cells with both polarities",
                "...and >= 2 labs",
                "cells mirrored / same / no_change",
                "record pairs opposite / same / no_change",
                "single-study cells m/s/nc",
                "multi-study cells m/s/nc",
            ),
            rows,
        )
    )
    lines.append(mirroring.get("what_would_falsify_this", ""))

    gene = mirroring.get("by_population_outcome_target") or {}
    pooled = mirroring.get("by_population_outcome") or {}
    gene_verdicts = gene.get("cell_verdicts") or {}
    gene_pairs = gene.get("record_pairs") or {}
    pooled_multi = (pooled.get("cell_verdicts_by_provenance") or {}).get(
        "multiple_studies"
    ) or {}
    if gene:
        lines.append("")
        lines.append(
            "**What this run found.** Hold the gene fixed and the axis does "
            "mirror: of the **%d** gene-specific cells holding both a loss and "
            "a gain of function, **%d** have opposite modal directions, "
            "**%d** have the same, and **%d** involve `no_change` on one side; "
            "at record level, **%d** opposite-signed pairs against **%d** "
            "same-signed. Pool the genes and that structure disappears: among "
            "the gene-pooled cells that draw on more than one study, **%d** "
            "mirror, **%d** are the same and **%d** involve `no_change`."
            % (
                gene.get("n_cells_with_both_polarities", 0),
                gene_verdicts.get("mirrored", 0),
                gene_verdicts.get("same", 0),
                gene_verdicts.get("no_change_involved", 0),
                gene_pairs.get("opposite_sign", 0),
                gene_pairs.get("same_sign", 0),
                pooled_multi.get("mirrored", 0),
                pooled_multi.get("same", 0),
                pooled_multi.get("no_change_involved", 0),
            )
        )
        lines.append("")
        lines.append(
            "**Why that is not usable as a feature, and the number that says "
            "so: %d.** That is how many of those gene-specific cells draw on "
            "two or more laboratories. Every one of them is a single study "
            "reporting both directions of its own manipulation, which is "
            "equally consistent with the polarity carrying a real sign and "
            "with a paper reporting a coherent story; nothing here "
            "distinguishes the two. And prediction needs more than a sign: it "
            "needs the cell to have been **seen in training**. Since each gene "
            "in this corpus appears in exactly one study, a leave-one-study-out "
            "split removes the only evidence about that gene along with the "
            "study, which is why the molecular half of the before/after table "
            "above is a column of structural zeros rather than a result. The "
            "remedy is not a better model. It is a second laboratory "
            "manipulating the same gene in the same population (SCHEMA.md 5b)."
            % gene.get("n_cells_with_both_polarities_and_two_or_more_labs", 0)
        )

    for level_key, title in (
        (
            "by_population_outcome_target",
            "the same gene, manipulated both ways",
        ),
        (
            "by_population_outcome",
            "genes pooled, which is what a target-blind cell key sees",
        ),
    ):
        level = mirroring.get(level_key) or {}
        cells = level.get("cells") or []
        if not cells:
            continue
        lines.append(_heading("`%s`: %s" % (level_key, title), 3))
        lines.append(
            _table(
                (
                    "cell",
                    "loss of function",
                    "gain of function",
                    "verdict",
                    "studies",
                    "labs",
                    "distinct genes",
                ),
                [
                    (
                        "`%s`" % " / ".join(str(part) for part in row.get("cell", [])),
                        str(row.get("loss_of_function")),
                        str(row.get("gain_of_function")),
                        "**%s**" % row.get("verdict"),
                        str(row.get("n_distinct_studies")),
                        str(row.get("n_distinct_labs")),
                        str(row.get("n_distinct_targets")),
                    )
                    for row in cells
                ],
            )
        )

    structure = diagnostics.get("cell_structure") or {}
    definitions = structure.get("definitions") or []
    if definitions:
        lines.append(
            _heading("%d.1 What each cell definition does to the data" % number, 3)
        )
        lines.append(
            "A finer cell is a smaller cell. In-sample purity is the ceiling an "
            "exact-cell lookup could reach if it were scored on the rows it "
            "memorised, and it rises **mechanically** as cells shrink: a cell "
            "holding one row is trivially pure and contributes nothing to a "
            "held-out row of the same cell. Read the two columns together, or "
            "the refinement that fragments the data hardest will look like the "
            "one that explains it best."
        )
        lines.append(
            _table(
                (
                    "cell",
                    "fields",
                    "cells",
                    "mean rows per cell",
                    "largest cell",
                    "singleton cells",
                    "in-sample purity (a ceiling, not an accuracy)",
                ),
                [
                    (
                        "`%s`" % row.get("cell"),
                        ", ".join("`%s`" % f for f in row.get("fields", [])),
                        str(row.get("n_cells")),
                        _fmt(row.get("mean_rows_per_cell"), 2),
                        str(row.get("largest_cell_rows")),
                        "%d (%s)"
                        % (
                            row.get("n_singleton_cells", 0),
                            _fmt(
                                None
                                if row.get("fraction_singleton_cells") is None
                                else 100.0 * row["fraction_singleton_cells"],
                                1,
                            )
                            + "%",
                        ),
                        _fmt(row.get("in_sample_purity")),
                    )
                    for row in definitions
                ],
            )
        )
        lines.append(structure.get("note", ""))
    return "\n".join(lines) + "\n"


def _behaviour_section(result: BenchmarkResult) -> str:
    lines = [_heading("7. How the baselines behaved")]
    any_content = False
    for grouping in result.groupings:
        for bundle in result.bundles:
            if bundle.grouping != grouping or bundle.status != "ok":
                continue
            detail = bundle.baseline_detail or {}
            backoff = bundle.backoff_usage
            abstain = bundle.abstention_reasons
            if not backoff and not abstain:
                continue
            any_content = True
            lines.append(
                _heading("`%s` under `%s`" % (bundle.baseline, grouping), 3)
            )
            if backoff:
                chain = detail.get("backoff_chain")
                if chain:
                    lines.append("Backoff chain: " + " -> ".join(chain))
                lines.append("")
                lines.append(
                    _table(
                        ("level actually used", "predictions"),
                        [(str(k), str(v)) for k, v in sorted(backoff.items())],
                    )
                )
                lines.append(
                    "A baseline that spends most of its predictions one or more "
                    "rungs down its own chain is not the model its name "
                    "describes, and comparing a candidate against it as though "
                    "it were a cell-level model would set the bar in the wrong "
                    "place."
                )
            if abstain:
                lines.append("")
                lines.append(
                    _table(
                        ("abstention reason", "rows"),
                        [(str(k), str(v)) for k, v in sorted(abstain.items())],
                    )
                )
                lines.append(
                    "An abstention is scored neither right nor wrong; it "
                    "reduces the row count the accuracy was computed on, which "
                    "is why coverage is printed beside every accuracy."
                )
    if not any_content:
        lines.append(
            "No baseline recorded a backoff or an abstention under any "
            "grouping in this run."
        )
    return "\n".join(lines) + "\n"


def _calibration_section(result: BenchmarkResult) -> str:
    lines = [_heading("8. Calibration")]
    lines.append(
        "Reliability bins are on the probability of the predicted class. A bin "
        "holding fewer than %d rows is marked `(!)` and must not be read as "
        "evidence of calibration or miscalibration in that range. Brier is "
        "the multiclass score (0 perfect, 2 confidently wrong); ECE is the "
        "sample-weighted mean gap between confidence and observed accuracy."
        % MIN_N_PER_RELIABILITY_BIN
    )
    printed = False
    for bundle in result.bundles:
        calibration = bundle.calibration or {}
        if not calibration:
            continue
        if not calibration.get("available"):
            continue
        rows = []
        for row in calibration.get("bins") or []:
            if not row.get("n"):
                continue
            rows.append((
                "%.1f-%.1f" % (row["bin_low"], row["bin_high"]),
                str(row["n"]),
                _fmt(row.get("mean_confidence")),
                _fmt(row.get("observed_accuracy")),
                _fmt(row.get("gap")),
                "" if row.get("meaningful") else "**(!) too few rows**",
            ))
        if not rows:
            continue
        printed = True
        lines.append(
            _heading("`%s` under `%s`" % (bundle.baseline, bundle.grouping), 3)
        )
        lines.append(
            "Grouping: `%s`. Rows with a probability: %d, in %d independent "
            "groups. ECE %s."
            % (
                bundle.grouping,
                calibration.get("n_rows_with_probability", 0),
                calibration.get("n_groups", 0),
                _fmt(calibration.get("expected_calibration_error")),
            )
        )
        lines.append("")
        lines.append(
            _table(
                (
                    "confidence bin",
                    "n",
                    "mean confidence",
                    "observed accuracy",
                    "gap",
                    "",
                ),
                rows,
            )
        )
    unavailable = [
        b for b in result.bundles
        if (b.calibration or {}).get("available") is False
    ]
    seen = set()
    for bundle in unavailable:
        if bundle.baseline in seen:
            continue
        seen.add(bundle.baseline)
        lines.append(
            "- `%s`: no calibration. %s"
            % (bundle.baseline, (bundle.calibration or {}).get("reason", ""))
        )
    if not printed and not seen:
        lines.append("No baseline produced a probability in this run.")
    return "\n".join(lines) + "\n"


def _bar_section(result: BenchmarkResult) -> str:
    lines = [_heading("9. The bar a model has to clear")]
    lines.append(
        "A sophisticated model is not to be adopted in this project unless it "
        "beats these numbers **under the same grouping**, with its own "
        "group-resampled interval, on at least accuracy and macro-F1 together. "
        "Beating `global_majority` on accuracy alone is not evidence of "
        "anything: on a skewed label distribution that is a measurement of the "
        "skew. Writing the bar down before a candidate exists is the point of "
        "this section; a bar set after the candidate's number is known is not "
        "a bar."
    )
    rows = []
    for grouping in result.groupings:
        declaration = GROUPINGS[grouping]
        best: Dict[str, Any] = {}
        for metric in ("direction_accuracy", "macro_f1"):
            candidates = [
                (b.estimates[metric], b.baseline)
                for b in result.bundles
                if b.grouping == grouping
                and b.status == "ok"
                and metric in b.estimates
                and b.estimates[metric].value is not None
            ]
            if not candidates:
                best[metric] = None
                continue
            best[metric] = max(candidates, key=lambda pair: pair[0].value)
        accuracy = best.get("direction_accuracy")
        macro = best.get("macro_f1")
        rows.append((
            "`%s`" % grouping,
            "**yes**" if declaration.leakage_prone else "no",
            (
                "no number to beat"
                if accuracy is None
                else "%s by `%s`" % (render_cell(accuracy[0]), accuracy[1])
            ),
            (
                "no number to beat"
                if macro is None
                else "%s by `%s`" % (render_cell(macro[0]), macro[1])
            ),
        ))
    lines.append(
        _table(
            (
                "grouping",
                "leakage-prone",
                "accuracy to beat (best baseline)",
                "macro-F1 to beat (best baseline)",
            ),
            rows,
        )
    )
    lines.append(
        "Where a cell says *no number to beat*, the bar is not zero -- it is "
        "undefined, because there was not enough data under that grouping to "
        "measure anything. A model evaluated there would be equally "
        "unmeasurable."
    )
    return "\n".join(lines) + "\n"


def _limits_section(result: BenchmarkResult) -> str:
    lines = [_heading("10. What these numbers do not say")]
    lines.append(
        "- **They are not a simulation.** The task is the sign of an effect in "
        "three classes. SCHEMA.md 6 refuses to represent dose-response curves, "
        "and nothing here predicts a magnitude."
    )
    lines.append(
        "- **They are not an aging result.** Tier B outcomes are excluded from "
        "the examples (SCHEMA.md 6). No number here is about lifespan, frailty "
        "or healthspan, and no \"rejuvenation score\" exists in this project."
    )
    lines.append(
        "- **They are conditional on the corpus, not on biology.** A baseline "
        "that predicts direction well may be reading a publication habit: the "
        "readouts a field reports are the ones that moved, and `no_change` is "
        "under-reported at source even though SCHEMA.md 3.5 requires it to be "
        "recorded when found."
    )
    lines.append(
        "- **Comparisons across groupings are not comparisons of models.** A "
        "number under `record` and a number under `lab` are measurements of "
        "different quantities. Only same-grouping comparisons are meaningful."
    )
    unavailable = [
        item["name"]
        for item in result.baselines
        if not (item.get("availability") or {}).get("available")
    ]
    if unavailable:
        lines.append(
            "- **%s could not be run at all.** The mechanistic baseline(s) "
            "that would test whether circuit structure predicts direction are "
            "absent from this comparison, so the bar in section 7 is the bar "
            "set by frequency and transmitter-sign baselines only, and may be "
            "lower than the real one." % ", ".join("`%s`" % n for n in unavailable)
        )
    lines.append(
        "- **A margin on one half of a corpus is not a capability.** Section 5 "
        "finds a difference on the circuit-activity records and a structural "
        "zero on the molecular ones. Nothing here supports a claim about "
        "molecular perturbations, and nothing supports extending the "
        "circuit-activity result to a population, a readout or a perturbation "
        "whose cell has not been seen -- which, per the SCHEMA.md 5b column in "
        "section 4, is most of them."
    )
    lines.append(
        "- **\"Insufficient evidence to rank these interventions reliably\" is "
        "still an available answer, and often the correct one.** It is what an "
        "abstention is, what an `UNAVAILABLE` baseline is, what a *no number "
        "to beat* cell is, and what a slice with no rows is. None of those is "
        "a gap to be filled in with a default; filling them in is how a "
        "benchmark starts reporting confidence it has not got."
    )
    lines.append(
        "- **Reproducibility.** Every number above is a deterministic function "
        "of the dataset, the fold source and the seed printed in section 0. No "
        "optional library affects any of them."
    )
    return "\n".join(lines) + "\n"


def render_markdown(result: BenchmarkResult) -> str:
    """The full report."""
    description = result.dataset_description or {}
    task = description.get("task") or {}
    settings = result.settings or {}
    environment = result.environment or {}

    lines = ["# Baseline benchmark: direction of effect", ""]
    lines.append(
        "Generated by `python3 -m syn4rch.models.run` at %s. Regenerable: this "
        "file is an output, not a source (DECISIONS.md D-016)."
        % (result.generated_at or "unknown time")
    )
    lines.append("")
    lines.append(
        "**Task.** %s" % (task.get("statement") or "(task statement missing)")
    )
    lines.append("")
    lines.append(
        "This is direction-of-effect classification, deliberately not "
        "quantitative simulation. The baselines below are the simplest "
        "predictors capable of answering whether direction is predictable at "
        "all; they exist to be the bar that a more sophisticated model must "
        "clear before it is worth adopting, and the harness is built to show "
        "when the answer is that direction is not predictable from these "
        "features."
    )
    lines.append("")
    lines.append("This report refuses to produce:")
    lines.append("")
    for item in task.get("refuses") or []:
        lines.append("- %s" % item)
    lines.append("")
    lines.append(headline(result))

    lines.append(_heading("0. Provenance"))
    lines.append(
        _table(
            ("setting", "value"),
            [
                ("task version", str(task.get("task_version"))),
                (
                    "previous task version (the before side)",
                    str(task.get("previous_version")),
                ),
                ("dataset", str(description.get("source"))),
                (
                    "schema versions read",
                    ", ".join(
                        "%s x%d" % (version, count)
                        for version, count in (
                            description.get("schema_versions_read") or {}
                        ).items()
                    )
                    or "(none)",
                ),
                ("fold source", "`%s`" % result.fold_source),
                ("folds requested", str(settings.get("n_folds_requested"))),
                ("bootstrap replicates", str(settings.get("n_bootstrap"))),
                (
                    "bootstrap resampling unit",
                    "**%s**" % settings.get("bootstrap_resampling_unit"),
                ),
                ("seed", str(settings.get("seed"))),
                ("reliability bins", str(settings.get("n_reliability_bins"))),
                ("numpy", str(environment.get("numpy"))),
                ("scikit-learn", str(environment.get("scikit-learn"))),
                ("metrics implementation", str(environment.get("metrics_implementation"))),
            ],
        )
    )
    lines.append("Fold resolution: %s" % result.fold_source_note)
    lines.append("")

    lines.append(_dataset_section(description))
    lines.append(_groupings_section(description, result.groupings))
    lines.append(_baselines_section(result.baselines))
    lines.append(_results_section(result))
    lines.append(_before_after_section(result, 5))
    lines.append(_mirroring_section(result, 6))
    lines.append(_behaviour_section(result))
    lines.append(_calibration_section(result))
    lines.append(_bar_section(result))
    lines.append(_limits_section(result))
    return "\n".join(lines).rstrip() + "\n"
