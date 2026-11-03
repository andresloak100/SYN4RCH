"""QC command line entry point.

    python3 -m syn4rch.qc.run data/curated/records.jsonl [--json reports/qc.json]

Exit status:

    0  no hard schema violations (warnings and conflicts do not fail the run)
    1  at least one hard schema violation, or a line that could not be read
    2  usage error (missing file, bad arguments)

Conflicts and non-independence never affect the exit code: they are findings
about the literature, not defects in the dataset.

Python 3.10 compatible, standard library only.
"""

import argparse
import json
import os
import sys

from ..schema.io import JsonlError, iter_jsonl
from .report import render_text, run_qc

__all__ = ["main", "build_parser"]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python3 -m syn4rch.qc.run",
        description="Quality control over a SYN4RCH Tier A records.jsonl file.",
    )
    parser.add_argument("path", help="path to records.jsonl")
    parser.add_argument(
        "--json",
        dest="json_path",
        default=None,
        metavar="PATH",
        help="also write the full report as JSON (e.g. reports/qc.json)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        metavar="N",
        help="show at most N findings per section in the text report",
    )
    parser.add_argument(
        "--no-annotate",
        action="store_true",
        help="do not populate conflicts_with in memory (report only)",
    )
    parser.add_argument(
        "--any-version",
        action="store_true",
        help="read records whose schema_version differs from this code's "
        "(for migration inspection only; off by default)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="suppress the text report (useful with --json)",
    )
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if not os.path.isfile(args.path):
        sys.stderr.write("error: no such file: %s\n" % args.path)
        return 2

    records = []
    load_errors = []
    try:
        for record in iter_jsonl(
            args.path, expect_version=None if args.any_version else "0.1.0"
        ):
            records.append(record)
    except JsonlError as exc:
        load_errors.append(str(exc))

    report = run_qc(records, annotate=not args.no_annotate)
    report.load_errors = load_errors

    if not args.quiet:
        sys.stdout.write(render_text(report, path=args.path, limit=args.limit) + "\n")

    if args.json_path:
        parent = os.path.dirname(os.path.abspath(args.json_path))
        if parent:
            os.makedirs(parent, exist_ok=True)
        with open(args.json_path, "w", encoding="utf-8") as handle:
            json.dump(report.to_dict(), handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        if not args.quiet:
            sys.stdout.write("wrote %s\n" % args.json_path)

    return 0 if report.ok else 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
