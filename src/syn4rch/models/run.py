"""Baseline benchmark command line entry point.

    python3 -m syn4rch.models.run
    python3 -m syn4rch.models.run --records data/curated/records.jsonl \
        --out reports/baseline_benchmark.md --json reports/baseline_benchmark.json

Writes ``reports/baseline_benchmark.md`` by default (tracked, per
DECISIONS.md D-016, because a measured claim with no committed provenance is
not reproducible by a reader who clones this repository).

Exit status:

    0  the report was written -- including when there was no data to measure,
       which is a result and is stated in the report
    1  the report could not be written (e.g. the output path is unwritable)
    2  usage error

**No data is not an error.** A missing or empty ``data/curated/records.jsonl``
produces a report that says so, in those words, and contains no metrics. The
alternative -- a non-zero exit, or a table of zeroes -- would either make an
empty dataset indistinguishable from a broken one, or make it look like a
measured result.

Nothing here writes to ``data/`` (SCHEMA.md 4). Python 3.10 compatible,
standard library only; ``numpy`` and ``scikit-learn`` are detected for the
provenance block and are not required.

Python 3.10 compatible.
"""

import argparse
import json
import os
import sys
from typing import Optional, Sequence

from .evaluate import run_benchmark
from .report import render_markdown
from .task import DEFAULT_GROUPINGS, GROUPINGS, default_records_path, load_dataset, repo_root

__all__ = ["main", "build_parser", "default_report_path"]


def default_report_path() -> str:
    root = repo_root()
    if root is None:
        return os.path.join("reports", "baseline_benchmark.md")
    return os.path.join(root, "reports", "baseline_benchmark.md")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python3 -m syn4rch.models.run",
        description=(
            "Benchmark the direction-of-effect baselines and write "
            "reports/baseline_benchmark.md. Every number is reported under a "
            "named grouping (SCHEMA.md 5)."
        ),
    )
    parser.add_argument(
        "--records",
        default=None,
        metavar="PATH",
        help="curated JSONL (default: <repo>/data/curated/records.jsonl)",
    )
    parser.add_argument(
        "--out",
        default=None,
        metavar="PATH",
        help="Markdown report path (default: <repo>/reports/baseline_benchmark.md)",
    )
    parser.add_argument(
        "--json",
        dest="json_path",
        default=None,
        metavar="PATH",
        help="also write the full result as JSON",
    )
    parser.add_argument(
        "--groupings",
        default=",".join(DEFAULT_GROUPINGS),
        metavar="A,B,C",
        help="comma-separated declared groupings (default: %s)"
        % ",".join(DEFAULT_GROUPINGS),
    )
    parser.add_argument(
        "--folds",
        type=int,
        default=5,
        metavar="N",
        help="folds requested when this module generates them (default: 5); "
        "ignored by a leave-one-group-out generator",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=1000,
        metavar="N",
        help="bootstrap replicates, resampled over groups (default: 1000)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20240101,
        metavar="N",
        help="bootstrap seed (default: 20240101); every number is a "
        "deterministic function of it",
    )
    parser.add_argument(
        "--reliability-bins",
        type=int,
        default=10,
        metavar="N",
        help="reliability diagram bins (default: 10)",
    )
    parser.add_argument(
        "--include-tier-b",
        action="store_true",
        help="keep Tier B outcomes among the examples. SCHEMA.md 6 forbids "
        "using Tier B as a training target, so this is for measuring how much "
        "Tier B data exists, not for producing a model result; the report "
        "states that the switch was set.",
    )
    parser.add_argument(
        "--any-version",
        action="store_true",
        help="read records whose schema_version differs from this code's "
        "(migration inspection only)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="do not print the summary to stdout",
    )
    return parser


def _write(path: str, text: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    groupings = [g.strip() for g in args.groupings.split(",") if g.strip()]
    undeclared = [g for g in groupings if g not in GROUPINGS]
    if undeclared:
        sys.stderr.write(
            "error: undeclared grouping(s) %s. A grouping must be registered "
            "in syn4rch.models.task.GROUPINGS so that its leakage properties "
            "are documented (SCHEMA.md 5). Declared: %s\n"
            % (", ".join(undeclared), ", ".join(sorted(GROUPINGS)))
        )
        return 2
    if not groupings:
        sys.stderr.write(
            "error: no groupings requested. SCHEMA.md 5 does not permit a "
            "cross-validated number without one.\n"
        )
        return 2
    if args.folds < 2:
        sys.stderr.write("error: --folds must be >= 2\n")
        return 2
    if args.bootstrap < 1:
        sys.stderr.write("error: --bootstrap must be >= 1\n")
        return 2

    records_path = args.records or default_records_path()
    dataset = load_dataset(
        records_path,
        exclude_tier_b=not args.include_tier_b,
        any_version=args.any_version,
    )

    result = run_benchmark(
        dataset,
        groupings=groupings,
        n_folds=args.folds,
        n_bootstrap=args.bootstrap,
        seed=args.seed,
        n_reliability_bins=args.reliability_bins,
    )

    out_path = args.out or default_report_path()
    try:
        _write(out_path, render_markdown(result))
    except OSError as exc:
        sys.stderr.write("error: could not write %s: %s\n" % (out_path, exc))
        return 1

    if args.json_path:
        try:
            parent = os.path.dirname(os.path.abspath(args.json_path))
            if parent:
                os.makedirs(parent, exist_ok=True)
            with open(args.json_path, "w", encoding="utf-8") as handle:
                json.dump(result.to_dict(), handle, indent=2, ensure_ascii=False)
                handle.write("\n")
        except OSError as exc:
            sys.stderr.write("error: could not write %s: %s\n" % (args.json_path, exc))
            return 1

    if not args.quiet:
        sys.stdout.write(_summary(dataset, result, out_path, args.json_path))
    return 0


def _summary(dataset, result, out_path, json_path) -> str:
    lines = []
    lines.append("records read:          %d" % dataset.ledger.n_input)
    lines.append("labellable examples:   %d" % len(dataset))
    lines.append("excluded (counted):    %d" % dataset.ledger.n_excluded)
    for error in dataset.load_errors:
        lines.append("load problem:          %s" % error)
    if dataset.is_empty:
        lines.append("")
        lines.append(
            "No labellable examples, so no metric was computed. The report "
            "says so rather than printing zeroes; an empty dataset is a "
            "result, not a failure."
        )
    lines.append("")
    lines.append("fold source:           %s" % result.fold_source)
    for item in result.baselines:
        availability = item.get("availability") or {}
        if not availability.get("available"):
            lines.append(
                "UNAVAILABLE baseline:  %s (%s)"
                % (item.get("name"), (availability.get("reason") or "").split(".")[0])
            )
    if not dataset.is_empty:
        lines.append("")
        for grouping in result.groupings:
            for bundle in result.bundles:
                if bundle.grouping != grouping or bundle.status != "ok":
                    continue
                estimate = bundle.estimates.get("direction_accuracy")
                if estimate is None:
                    continue
                lines.append(
                    "  %-28s %-20s %s"
                    % (bundle.baseline, grouping, estimate.render())
                )
    lines.append("")
    lines.append("wrote %s" % out_path)
    if json_path:
        lines.append("wrote %s" % json_path)
    return "\n".join(lines) + "\n"


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
