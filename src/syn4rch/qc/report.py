"""Assembly and rendering of the QC report.

``run_qc(records)`` runs every check and returns a :class:`QCReport`. The
report distinguishes, at the top level and in the exit code:

* **hard schema violations** -- records that may not enter ``data/curated/``;
* **warnings** -- things to look at;
* **info** -- properties of the literature (conflicts, non-independence).

``annotate`` controls whether ``conflicts_with`` is written onto the records
in memory. The caller decides whether to persist them; this module never
writes to ``data/``.

Python 3.10 compatible, standard library only.
"""

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

from ..schema.validate import Severity as SchemaSeverity
from ..schema.validate import validate
from .comparability import comparability_coverage, find_missing_comparability_fields, unclear_valence_count
from .conflicts import annotate_conflicts, find_conflicts, key_comparison
from .keys import condition_specificity
from .duplicates import run_duplicate_checks
from .findings import QCFinding, QCSeverity
from .units import check_units

__all__ = ["QCReport", "run_qc", "render_text"]


def _schema_findings(records) -> List[QCFinding]:
    out: List[QCFinding] = []
    for record in records:
        for finding in validate(record):
            out.append(
                QCFinding(
                    code=finding.code,
                    severity=(
                        QCSeverity.HARD
                        if finding.severity is SchemaSeverity.HARD
                        else QCSeverity.WARNING
                    ),
                    message=finding.message,
                    record_ids=[record.record_id] if record.record_id else [],
                    field_name=finding.field,
                )
            )
    return out


@dataclass
class QCReport:
    """Everything the QC pass found, plus a few counts for the header."""

    n_records: int = 0
    findings: List[QCFinding] = field(default_factory=list)
    stats: Dict[str, Any] = field(default_factory=dict)
    load_errors: List[str] = field(default_factory=list)

    def by_severity(self, severity: QCSeverity) -> List[QCFinding]:
        return [f for f in self.findings if f.severity is severity]

    @property
    def hard(self) -> List[QCFinding]:
        return self.by_severity(QCSeverity.HARD)

    @property
    def warnings(self) -> List[QCFinding]:
        return self.by_severity(QCSeverity.WARNING)

    @property
    def info(self) -> List[QCFinding]:
        return self.by_severity(QCSeverity.INFO)

    @property
    def blocking_record_ids(self) -> List[str]:
        ids = set()
        for finding in self.hard:
            ids.update(finding.record_ids)
        return sorted(ids)

    @property
    def ok(self) -> bool:
        """True if nothing blocks admission to ``data/curated/``."""
        return not self.hard and not self.load_errors

    def counts_by_code(self) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for finding in self.findings:
            counts[finding.code] = counts.get(finding.code, 0) + 1
        return dict(sorted(counts.items()))

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_records": self.n_records,
            "ok": self.ok,
            "load_errors": list(self.load_errors),
            "counts": {
                "hard": len(self.hard),
                "warning": len(self.warnings),
                "info": len(self.info),
            },
            "counts_by_code": self.counts_by_code(),
            "stats": dict(self.stats),
            "blocking_record_ids": self.blocking_record_ids,
            "findings": [f.to_dict() for f in self.findings],
        }


def run_qc(
    records,
    annotate: bool = True,
    vocab_loader: Optional[Callable[[str], Dict]] = None,
) -> QCReport:
    """Run every QC check over ``records``.

    ``annotate=True`` populates ``conflicts_with`` on both sides of each
    conflict (in memory). No check ever edits any other field, drops a record
    or reconciles a disagreement.
    """
    records = list(records)
    report = QCReport(n_records=len(records))

    report.findings.extend(_schema_findings(records))
    report.findings.extend(run_duplicate_checks(records))
    report.findings.extend(check_units(records, loader=vocab_loader))
    report.findings.extend(find_missing_comparability_fields(records))
    conflicts = annotate_conflicts(records) if annotate else find_conflicts(records)
    report.findings.extend(conflicts)

    coverage = comparability_coverage(records)
    report.stats = {
        "n_studies": len({r.study_id for r in records if r.study_id}),
        "n_labs": len(
            {
                (r.last_author or "", r.last_author_affiliation or "")
                for r in records
                if r.last_author
            }
        ),
        "n_populations_unresolved": sum(1 for r in records if r.population_id is None),
        "n_abstract_only": sum(1 for r in records if r.full_text_available is not True),
        "n_valence_unclear": unclear_valence_count(records),
        "n_conflicting_records": len(
            {rid for f in conflicts for rid in f.record_ids}
        ),
        "comparability_coverage": coverage,
        "direction_counts": _direction_counts(records),
        "claim_key_comparison": key_comparison(records),
        "condition_specificity": _specificity_counts(records),
    }
    return report


def _specificity_counts(records) -> Dict[str, int]:
    """How many records know 0, 1, 2, 3 or 4 of their condition slots.

    Reported because a conflict or a replication found between two records of
    specificity 0 matched on "both unknown", not on "both the same", and a
    report that prints one number for both is overclaiming (DECISIONS D-023).
    """
    counts = {str(n): 0 for n in range(5)}
    for record in records:
        counts[str(condition_specificity(record))] += 1
    return counts


def _direction_counts(records) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for record in records:
        key = record.direction.value if record.direction else "null"
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def _section(title, findings, limit=None) -> List[str]:
    lines = ["%s (%d)" % (title, len(findings)), "-" * 72]
    if not findings:
        lines.append("  none")
    shown = findings if limit is None else findings[:limit]
    for finding in shown:
        ids = ", ".join(finding.record_ids) if finding.record_ids else "-"
        lines.append("  %s  [%s]" % (finding.code, ids))
        lines.append("      %s" % finding.message)
    if limit is not None and len(findings) > limit:
        lines.append("  ... %d more" % (len(findings) - limit))
    lines.append("")
    return lines


def render_text(report: QCReport, path=None, limit=None) -> str:
    """Human-readable report."""
    lines = [
        "=" * 72,
        "SYN4RCH QC report",
        "=" * 72,
        "input              : %s" % (path if path else "<in memory>"),
        "records read       : %d" % report.n_records,
        "studies            : %s" % report.stats.get("n_studies", "-"),
        "labs identified    : %s" % report.stats.get("n_labs", "-"),
        "hard violations    : %d  (block admission to data/curated/)"
        % len(report.hard),
        "warnings           : %d" % len(report.warnings),
        "info (data, not defects): %d" % len(report.info),
        "",
    ]

    if report.load_errors:
        lines.append("LOAD ERRORS (%d)" % len(report.load_errors))
        lines.append("-" * 72)
        for err in report.load_errors:
            lines.append("  %s" % err)
        lines.append("")

    coverage = report.stats.get("comparability_coverage", {})
    if coverage:
        total = coverage.get("_total_records", report.n_records) or 1
        lines.append("Comparability fields reported (not the condition dimension)")
        lines.append("-" * 72)
        for name in ("control_types", "n", "species", "valence"):
            if name in coverage:
                lines.append(
                    "  %-14s %4d/%d reported" % (name, coverage[name], total)
                )
        lines.append("")

    directions = report.stats.get("direction_counts", {})
    if directions:
        lines.append("Reported directions")
        lines.append("-" * 72)
        for key, value in directions.items():
            lines.append("  %-14s %4d" % (key, value))
        lines.append(
            "  (no_change is kept deliberately; dropping null results biases "
            "direction -- SCHEMA.md 3.5)"
        )
        lines.append("")

    comparison = report.stats.get("claim_key_comparison") or {}
    if comparison.get("strict") and comparison.get("pooled"):
        strict, pooled = comparison["strict"], comparison["pooled"]
        lines.append("What counts as the same claim (schema v0.2.0 vs v0.1.0)")
        lines.append("-" * 72)
        lines.append("  %-34s %8s %8s" % ("", "v0.2.0", "v0.1.0"))
        lines.append(
            "  %-34s %8d %8d" % ("distinct claims", strict["n_cells"], pooled["n_cells"])
        )
        lines.append(
            "  %-34s %8d %8d"
            % (
                "claims with a cross-study conflict",
                strict["n_cells_with_conflict"],
                pooled["n_cells_with_conflict"],
            )
        )
        for code in sorted(
            set(strict["cells_by_code"]) | set(pooled["cells_by_code"])
        ):
            lines.append(
                "    %-32s %8d %8d"
                % (
                    code.replace("QC_CONFLICT_", "").lower(),
                    strict["cells_by_code"].get(code, 0),
                    pooled["cells_by_code"].get(code, 0),
                )
            )
        lines.append(
            "  v0.1.0 keyed on (population, valence, outcome) and so pooled a male"
        )
        lines.append(
            "  cohort with a female one, and a gene knockout with an overexpression."
        )
        lines.append(
            "  A claim counted under v0.1.0 but not under v0.2.0 is not a resolved"
        )
        lines.append(
            "  disagreement; it is a comparison that should never have been made."
        )
        lines.append("")

    specificity = report.stats.get("condition_specificity") or {}
    if specificity:
        lines.append("Condition specificity (slots known of sex/diet/age/challenge)")
        lines.append("-" * 72)
        for slots in sorted(specificity):
            lines.append("  %s of 4 known   %4d records" % (slots, specificity[slots]))
        blind = specificity.get("0", 0)
        if blind:
            lines.append(
                "  %d record(s) know none, so they can only ever be compared with"
                % blind
            )
            lines.append(
                "  other condition-blind records, and any agreement or disagreement"
            )
            lines.append(
                "  found that way may be condition heterogeneity instead."
            )
        lines.append("")

    lines.extend(
        _section("HARD schema violations -- must not enter data/curated/", report.hard, limit)
    )
    lines.extend(_section("WARNINGS -- review, not blocking", report.warnings, limit))
    lines.extend(
        _section(
            "INFO -- conflicts and possible non-independence (these are DATA; "
            "nothing was resolved, averaged or dropped)",
            report.info,
            limit,
        )
    )

    lines.append("Findings by code")
    lines.append("-" * 72)
    for code, count in report.counts_by_code().items():
        lines.append("  %-38s %4d" % (code, count))
    lines.append("")
    lines.append(
        "VERDICT: %s"
        % (
            "no hard schema violations"
            if report.ok
            else "%d hard violation(s) across %d record(s)"
            % (len(report.hard), len(report.blocking_record_ids))
        )
    )
    return "\n".join(lines)
