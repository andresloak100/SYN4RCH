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

from typing import Any, Dict, List, Optional, Sequence

from .evaluate import (
    MIN_GROUPS_FOR_CI,
    MIN_N_PER_RELIABILITY_BIN,
    MIN_ROWS_FOR_METRIC,
    BenchmarkResult,
    Estimate,
    MetricBundle,
)
from .task import EXCLUSION_REASONS, GROUPINGS, LABELS

__all__ = ["render_markdown", "render_cell"]

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
            "available" if availability.get("available") else "**UNAVAILABLE**",
            "yes" if item.get("emits_proba") else "no",
        ))
    lines.append(
        _table(("baseline", "what it does", "status", "emits probabilities"), rows)
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
                + ["--", "--", "--"]
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
    header = (
        ["baseline"]
        + [label for _, label in _COLUMNS]
        + ["coverage (answered)", "folds", "fold source"]
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


def _behaviour_section(result: BenchmarkResult) -> str:
    lines = [_heading("5. How the baselines behaved")]
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
    lines = [_heading("6. Calibration")]
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
    lines = [_heading("7. The bar a model has to clear")]
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
    lines = [_heading("8. What these numbers do not say")]
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
        "quantitative simulation. The five baselines below are the simplest "
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

    lines.append(_heading("0. Provenance"))
    lines.append(
        _table(
            ("setting", "value"),
            [
                ("task version", str(task.get("task_version"))),
                ("dataset", str(description.get("source"))),
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
    lines.append(_behaviour_section(result))
    lines.append(_calibration_section(result))
    lines.append(_bar_section(result))
    lines.append(_limits_section(result))
    return "\n".join(lines).rstrip() + "\n"
