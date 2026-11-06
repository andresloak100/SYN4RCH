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

from ..schema.record import CONDITION_DIMENSIONS
from ..schema.validate import Severity as SchemaSeverity
from ..schema.validate import unevidenced_conditions, validate
from .comparability import comparability_coverage, find_missing_comparability_fields, unclear_valence_count
from .conflicts import annotate_conflicts, find_conflicts, key_comparison
from .duplicates import run_duplicate_checks
from .keys import CONDITION_FIELDS, NULL, condition_specificity, lab_key
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
        # lab_key, not a second definition of lab identity. This counted
        # (last_author, last_author_affiliation) inline until D-056 and so
        # kept reporting 75 labs for 75 studies after the shared key was
        # corrected -- the exact duplication keys.py exists to prevent.
        "n_labs": len({lab_key(r) for r in records if lab_key(r) is not None}),
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
        "condition_evidence": _condition_evidence_counts(records),
    }
    return report


def _specificity_counts(records) -> Dict[str, int]:
    """How many records know 0, 1, ... or all of their condition slots.

    Reported because a conflict or a replication found between two records of
    specificity 0 matched on "both unknown", not on "both the same", and a
    report that prints one number for both is overclaiming (DECISIONS D-023).

    The number of buckets is ``len(CONDITION_FIELDS) + 1``, never a literal:
    it was ``range(5)`` for four slots until D-060 added a fifth, which left
    the top bucket uncounted while the printed denominator said five.
    """
    counts = {str(n): 0 for n in range(len(CONDITION_FIELDS) + 1)}
    for record in records:
        counts[str(condition_specificity(record))] += 1
    return counts


#: The one condition value that asserts an absence. It carries no evidence
#: entry by construction, so it is neither evidenced nor a gap. Spelled here
#: rather than imported because :mod:`syn4rch.qc.keys` pools it with ``None``
#: for key-building, and this section needs the two apart.
_NOT_REPORTED = "not_reported"


def _condition_value(record, dimension: str):
    """One condition dimension's own value, as a plain string or ``None``."""
    raw = getattr(record, dimension, None)
    return getattr(raw, "value", raw)


def _condition_evidence_counts(records) -> Dict[str, Any]:
    """What the condition dimension can be *checked* against (schema v0.3.0).

    Distinct from :func:`_specificity_counts`, which counts how many
    dimensions are *known*. A dimension can be known and unquoted: all 362
    records curated before v0.3.0 are, which is the gap W017 counts.

    Three counts here are deliberately kept apart rather than summed:

    * ``evidenced`` -- a verbatim sentence was found and recorded;
    * ``not_reported`` -- the paper was read and states nothing, so there is
      no sentence to quote and no evidence entry by construction;
    * ``null`` -- the dimension was never addressed.

    Adding the first two into a "covered" share would read as provenance for
    an assertion that by definition has none.

    ``by_rule`` carries ``entries`` *and* ``records`` per rule because one
    record can be touched by one rule on three dimensions: the entry count
    says how often a rule fired, and only the record count says how much of
    the dataset a wrong rule would affect.

    ``study_reach`` is here because the per-dimension shares cannot tell
    evidence drawn from two papers from evidence spread across forty.
    """
    records = list(records)
    dimensions: Dict[str, Dict[str, int]] = {
        name: {"positive": 0, "evidenced": 0, "not_reported": 0, "null": 0}
        for name in CONDITION_DIMENSIONS
    }
    by_source: Dict[str, int] = {"passage": 0, "methods": 0}
    rules: Dict[str, Dict[str, int]] = {}
    n_entries = 0
    gap_records = 0
    gap_values = 0
    reached_studies = set()
    evidenced_without_study = 0

    for record in records:
        entries = list(record.condition_evidence or [])
        # The gap is read from the rule that raises W017 rather than recomputed
        # here. Two implementations of one rule agree until the day one of them
        # is corrected, and this section's whole claim is that its numbers are
        # the validator's numbers.
        missing = set(unevidenced_conditions(record))
        record_gap = 0
        record_evidenced = 0
        for name in CONDITION_DIMENSIONS:
            row = dimensions[name]
            value = _condition_value(record, name)
            if value is None:
                row["null"] += 1
            elif value == _NOT_REPORTED:
                row["not_reported"] += 1
            else:
                row["positive"] += 1
                if name in missing:
                    record_gap += 1
                else:
                    row["evidenced"] += 1
                    record_evidenced += 1
        if record_gap:
            gap_records += 1
            gap_values += record_gap
        if record_evidenced:
            if record.study_id:
                reached_studies.add(record.study_id)
            else:
                evidenced_without_study += 1

        seen_rules = set()
        for entry in entries:
            n_entries += 1
            source = getattr(entry, "source", None)
            token = getattr(source, "value", source) or NULL
            by_source[token] = by_source.get(token, 0) + 1
            rule = getattr(entry, "rule", None) or NULL
            tally = rules.setdefault(rule, {"entries": 0, "records": 0})
            tally["entries"] += 1
            if rule not in seen_rules:
                seen_rules.add(rule)
                tally["records"] += 1

    return {
        "_total_records": len(records),
        "n_entries": n_entries,
        "dimensions": dimensions,
        "by_source": by_source,
        "by_rule": dict(
            sorted(
                rules.items(),
                key=lambda kv: (-kv[1]["entries"], -kv[1]["records"], kv[0]),
            )
        ),
        "gap": {"records": gap_records, "values": gap_values},
        "study_reach": {
            "n_studies": len({r.study_id for r in records if r.study_id}),
            "n_studies_with_evidenced_value": len(reached_studies),
            "n_evidenced_records_without_study_id": evidenced_without_study,
        },
    }


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


def _render_condition_evidence(stats: Dict[str, Any]) -> List[str]:
    """The condition-evidence section, as report lines.

    Prints counts only. There is deliberately no single figure of merit for a
    record's or a dataset's condition provenance: the four dimensions are
    evidenced through different routes, and a dataset that is fully evidenced
    on sex and not at all on challenge is not 50% of anything.
    """
    if not stats or not stats.get("_total_records"):
        return []
    dimensions = stats.get("dimensions") or {}
    gap = stats.get("gap") or {}
    reach = stats.get("study_reach") or {}

    lines = [
        "Condition evidence (what is quoted, not what is known)",
        "-" * 72,
        "  %-12s %8s %9s %12s %5s"
        % ("dimension", "positive", "evidenced", "not_reported", "null"),
    ]
    for name in CONDITION_DIMENSIONS:
        row = dimensions.get(name) or {}
        lines.append(
            "  %-12s %8d %9d %12d %5d"
            % (
                name,
                row.get("positive", 0),
                row.get("evidenced", 0),
                row.get("not_reported", 0),
                row.get("null", 0),
            )
        )

    by_source = stats.get("by_source") or {}
    lines.append("  evidence entries by source (%d in all)" % stats.get("n_entries", 0))
    ordered = ["passage", "methods"] + sorted(set(by_source) - {"passage", "methods"})
    for token in ordered:
        lines.append("    %-12s %8d" % (token, by_source.get(token, 0)))

    by_rule = stats.get("by_rule") or {}
    lines.append("  evidence entries by rule (most frequent first -- the audit handle)")
    if not by_rule:
        lines.append("    none")
    for rule, tally in by_rule.items():
        lines.append(
            "    %-22s %6d entries  %5d record(s)"
            % (rule, tally.get("entries", 0), tally.get("records", 0))
        )

    lines.append(
        "  positive values with no evidence entry: %d, on %d record(s) "
        "(the W017 population)" % (gap.get("values", 0), gap.get("records", 0))
    )
    lines.append(
        "  studies with at least one evidenced value: %d of %d"
        % (
            reach.get("n_studies_with_evidenced_value", 0),
            reach.get("n_studies", 0),
        )
    )
    orphans = reach.get("n_evidenced_records_without_study_id", 0)
    if orphans:
        lines.append(
            "    (%d record(s) carry evidence but no study_id and so cannot be"
            % orphans
        )
        lines.append("    counted into that reach figure)")
    lines.append(
        "  Reach is reported because no single coverage figure can tell evidence"
    )
    lines.append("  drawn from two papers from evidence spread across forty.")
    lines.append(
        "  An evidence entry means a verbatim sentence was found and checked, while"
    )
    lines.append(
        "  a not_reported value asserts that the paper carries no such sentence and"
    )
    lines.append(
        '  so has none by construction -- the two must never be added into one'
    )
    lines.append('  "covered" figure.')
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
        lines.append(
            "Condition specificity (slots known of %s)"
            % "/".join(CONDITION_FIELDS)
        )
        lines.append("-" * 72)
        for slots in sorted(specificity):
            lines.append(
                "  %s of %d known   %4d records"
                % (slots, len(CONDITION_FIELDS), specificity[slots])
            )
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

    lines.extend(_render_condition_evidence(report.stats.get("condition_evidence") or {}))

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
