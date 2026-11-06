"""Assembly and Markdown rendering of the Phase 2 feasibility report.

The report answers the five Phase 2 questions in order, from counted records
only. Three rules govern it:

1. **It opens with the verdict.** :func:`top_line` is one sentence and it is
   the first thing in the file. A reader who stops there must not be left with
   a more optimistic impression than a reader who finishes.
1b. **It shows the before and the after of the schema v0.2.0 key change, not
   only the after.** Two sections exist for nothing else: the claim-key
   comparison (:func:`syn4rch.feasibility.replication.compare_claim_keys`) and
   the usability-predicate comparison
   (:func:`syn4rch.feasibility.coverage.compare_usability`). A reader has to be
   able to see what the change bought -- apparent contradictions that were
   condition or perturbation heterogeneity -- and what it cost -- replication
   evidence given up by a stricter key. Printing only the new numbers would
   make the change unfalsifiable.
1c. **No replication or agreement figure appears without the minimum condition
   specificity of its cell.** The condition key pools "not stated" with "not
   asked", so a cell can match on *both unknown* rather than on *both the
   same*; a claim replicated across three fully specified studies and one
   replicated across three studies that specify nothing are different evidence,
   and one number for both overclaims. Every table below carries the
   specificity column and every count is also broken down by it.
2. **A small corpus is described as small, in that sentence.** Below the floors
   below, the top line says the dataset is too small and the tables are
   labelled a census of what exists rather than statistics about a field.
   Degenerate numbers are suppressed, not formatted: a 1x1 matrix has no
   meaningful sparsity, and "0% empty" on one cell would read as complete
   coverage.
3. **Nothing is estimated.** Every figure is a count, a ratio of counts, or a
   confidence interval on such a ratio with its method named. No imputation,
   no projection of how many records curation "should" yield.

Python 3.10 compatible, standard library only.
"""

import datetime
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from ..qc.keys import CONDITION_FIELDS
from . import gates as gates_mod
from .coverage import (
    EXCLUSION_CODES,
    CoverageMatrix,
    UsabilityComparison,
    build_matrix,
    compare_usability,
    exclusion_note,
    exclusion_reasons,
)
from .gates import PairVerdict, TaskVerdict, pair_verdicts, summarize_verdicts
from .replication import (
    CONFIDENCE_LEVEL,
    KEY_POOLED,
    KEY_STRICT,
    N_CONDITION_SLOTS,
    KeyComparison,
    ReplicationCell,
    build_replication_cells,
    compare_claim_keys,
    summarize_replication,
)
from .splits import FoldSet, all_fold_sets

__all__ = [
    "MAX_GRID_DIM",
    "MAX_TABLE_ROWS",
    "MIN_MATRIX_DIM_FOR_SPARSITY",
    "SMALL_CORPUS_RECORDS",
    "SMALL_CORPUS_STUDIES",
    "FeasibilityReport",
    "build_report",
    "render_markdown",
    "top_line",
]

#: Below this many usable records the corpus is described as too small in the
#: top line. 30 = 10 x :data:`syn4rch.feasibility.gates.MIN_STUDIES_FEASIBLE`:
#: beneath it, not even ten tasks could clear the study gate in the best case,
#: so no distributional claim in this report would have ten points behind it.
SMALL_CORPUS_RECORDS = 30

#: Below this many distinct studies the corpus is described as too small.
#: Ten, because :data:`syn4rch.feasibility.gates.MIN_LABS_FEASIBLE` is 2 and
#: each study appears in one lab: with fewer than ten studies the corpus is a
#: reading list, and a matrix sparsity computed from it measures which papers
#: have been curated so far, not which relationships are unstudied.
SMALL_CORPUS_STUDIES = 10

#: Sparsity is printed only when both axes have at least this many entries. A
#: 1xN matrix's emptiness is an artefact of its own axes being drawn from the
#: records that exist.
MIN_MATRIX_DIM_FOR_SPARSITY = 2

#: The dense grid (with its explicit zeros) is drawn only up to this many rows
#: and columns; beyond it the grid is unreadable in Markdown and the zeros are
#: reported as a count instead.
MAX_GRID_DIM = 12

#: Row cap for the long tables; the remainder is reported as "... N more".
MAX_TABLE_ROWS = 40


@dataclass
class FeasibilityReport:
    """Everything the feasibility pass measured."""

    coverage: CoverageMatrix = field(default_factory=CoverageMatrix)
    cells: List[ReplicationCell] = field(default_factory=list)
    verdicts: List[TaskVerdict] = field(default_factory=list)
    pairs: List[PairVerdict] = field(default_factory=list)
    fold_sets: Dict[str, FoldSet] = field(default_factory=dict)
    replication_summary: Dict[str, Any] = field(default_factory=dict)
    verdict_summary: Dict[str, Any] = field(default_factory=dict)
    #: The v0.1.0 and v0.2.0 claim keys measured on the same records. The
    #: before/after of the schema change, and the reason
    #: :func:`syn4rch.qc.keys.claim_key_pooled` is kept.
    key_comparison: Optional[KeyComparison] = None
    #: The v0.1.0 and v0.2.0 usability predicates on the same records.
    usability_comparison: Optional[UsabilityComparison] = None
    #: Verdicts under the pooled key, for the gate-arithmetic comparison only.
    #: Never reported as the current verdict set.
    pooled_verdict_summary: Dict[str, Any] = field(default_factory=dict)
    pooled_pairs: List[PairVerdict] = field(default_factory=list)
    source_path: Optional[str] = None
    source_exists: bool = True
    load_errors: List[str] = field(default_factory=list)
    generated_at: Optional[str] = None

    # --- shape of the corpus ---------------------------------------------

    @property
    def n_records_read(self) -> int:
        return self.coverage.n_records_read

    @property
    def n_records_usable(self) -> int:
        return self.coverage.n_records_usable

    @property
    def is_empty(self) -> bool:
        return self.n_records_read == 0

    @property
    def has_usable_records(self) -> bool:
        return self.n_records_usable > 0

    @property
    def is_small_corpus(self) -> bool:
        """True if the corpus is too small for the tables to be statistics."""
        return (
            self.n_records_usable < SMALL_CORPUS_RECORDS
            or self.coverage.n_studies < SMALL_CORPUS_STUDIES
        )

    @property
    def n_feasible(self) -> int:
        return sum(1 for p in self.pairs if p.verdict == gates_mod.VERDICT_FEASIBLE)

    @property
    def n_marginal(self) -> int:
        return sum(1 for p in self.pairs if p.verdict == gates_mod.VERDICT_MARGINAL)

    @property
    def n_insufficient(self) -> int:
        return sum(
            1 for p in self.pairs if p.verdict == gates_mod.VERDICT_INSUFFICIENT
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "generated_at": self.generated_at,
            "source_path": self.source_path,
            "source_exists": self.source_exists,
            "load_errors": list(self.load_errors),
            "top_line": top_line(self),
            "is_empty": self.is_empty,
            "is_small_corpus": self.is_small_corpus,
            "coverage": self.coverage.to_dict(),
            "key_comparison": (
                self.key_comparison.to_dict() if self.key_comparison else None
            ),
            "usability_comparison": (
                self.usability_comparison.to_dict()
                if self.usability_comparison
                else None
            ),
            "pooled_verdict_summary": dict(self.pooled_verdict_summary),
            "replication_summary": dict(self.replication_summary),
            "verdict_summary": dict(self.verdict_summary),
            "replication_cells": [c.to_dict() for c in self.cells],
            "task_verdicts": [v.to_dict() for v in self.verdicts],
            "pair_verdicts": [p.to_dict() for p in self.pairs],
            "fold_sets": {
                name: {
                    k: v
                    for k, v in fs.to_dict().items()
                    # The per-fold record id lists are large and are an index
                    # into the dataset, not a finding; the shape is what the
                    # report is about.
                    if k != "folds"
                }
                for name, fs in self.fold_sets.items()
            },
        }


def build_report(
    records: Sequence,
    source_path: Optional[str] = None,
    source_exists: bool = True,
    load_errors: Optional[Sequence[str]] = None,
    generated_at: Optional[str] = None,
) -> FeasibilityReport:
    """Measure ``records`` and assemble the report.

    Safe on an empty sequence: every section degrades to an explicit zero and
    :func:`top_line` says the dataset is empty.
    """
    records = list(records)
    coverage = build_matrix(records)
    # Replication is computed over the usable records only, for the same
    # reason the matrix axes are: a record with no stated direction cannot
    # agree or disagree with anything. The held-out records stay counted in
    # coverage.excluded_counts.
    usable = [r for r in records if not _is_excluded(r)]
    cells = build_replication_cells(usable, KEY_STRICT)
    verdicts = gates_mod.evaluate_cells(cells)

    # The before/after. Both comparisons run over the *same* usable set as the
    # headline numbers so that the only thing varying between the two columns
    # is the key, and over all records for the usability predicate, where the
    # record set is itself what differs.
    pooled_verdicts = gates_mod.evaluate_records(usable, KEY_POOLED)

    return FeasibilityReport(
        coverage=coverage,
        cells=cells,
        verdicts=verdicts,
        pairs=pair_verdicts(verdicts),
        fold_sets=all_fold_sets(usable),
        replication_summary=summarize_replication(cells),
        verdict_summary=summarize_verdicts(verdicts),
        key_comparison=compare_claim_keys(usable),
        usability_comparison=compare_usability(records),
        pooled_verdict_summary=summarize_verdicts(pooled_verdicts),
        pooled_pairs=pair_verdicts(pooled_verdicts),
        source_path=source_path,
        source_exists=source_exists,
        load_errors=list(load_errors or []),
        generated_at=generated_at or _utc_now(),
    )


def _is_excluded(record) -> bool:
    return bool(exclusion_reasons(record))


def _utc_now() -> str:
    return (
        datetime.datetime.now(datetime.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


# --- the top line -----------------------------------------------------------


def top_line(report: "FeasibilityReport") -> str:
    """One sentence. The first thing in the file, and the whole verdict.

    It never leads with a statistic computed over a corpus too small to
    support it.
    """
    path = report.source_path or "the curated records file"
    if report.load_errors and report.n_records_read == 0:
        return (
            "**Verdict: not computed — `%s` could not be read (%d load error(s)), "
            "so no feasibility question is answered below.**" % (path, len(report.load_errors))
        )
    if report.n_records_read == 0:
        reason = "does not exist yet" if not report.source_exists else "is empty"
        return (
            "**Verdict: insufficient evidence — `%s` %s (0 records), so none of "
            "the five Phase 2 questions can be answered and no prediction task "
            "is modellable.**" % (path, reason)
        )
    if not report.has_usable_records:
        return (
            "**Verdict: insufficient evidence — %d record(s) were read but none "
            "carries a usable (population, intervention, measured outcome) "
            "triple, so no population x outcome task exists to assess.**"
            % report.n_records_read
        )

    shape = "%d usable record(s) from %d stud%s and %d identified lab(s)" % (
        report.n_records_usable,
        report.coverage.n_studies,
        "y" if report.coverage.n_studies == 1 else "ies",
        report.coverage.n_labs,
    )
    n_pairs = len(report.pairs)
    if report.is_small_corpus:
        if report.n_feasible:
            return (
                "**Verdict: insufficient evidence overall — the corpus is %s, "
                "below the floor of %d records and %d studies at which these "
                "tables become statistics, and only %d of %d (population, "
                "outcome) task(s) clear the promotion gates.**"
                % (
                    shape,
                    SMALL_CORPUS_RECORDS,
                    SMALL_CORPUS_STUDIES,
                    report.n_feasible,
                    n_pairs,
                )
            )
        return (
            "**Verdict: insufficient evidence — the corpus is %s, below the "
            "floor of %d records and %d studies at which these tables become "
            "statistics, and no (population, outcome) task clears the "
            "promotion gates; the counts below are a census of what has been "
            "curated, not a description of the field.**"
            % (shape, SMALL_CORPUS_RECORDS, SMALL_CORPUS_STUDIES)
        )
    if report.n_feasible == 0:
        return (
            "**Verdict: insufficient evidence — none of the %d (population, "
            "outcome) task(s) in %s clears the promotion gates (%d marginal, "
            "%d insufficient).**"
            % (n_pairs, shape, report.n_marginal, report.n_insufficient)
        )
    return (
        "**Verdict: %d of %d (population, outcome) task(s) are feasible for "
        "predictive modelling (%d marginal, %d insufficient), from %s.**"
        % (
            report.n_feasible,
            n_pairs,
            report.n_marginal,
            report.n_insufficient,
            shape,
        )
    )


# --- rendering helpers ------------------------------------------------------


def _fmt_ratio(value: Optional[float]) -> str:
    return "undefined" if value is None else "%.2f" % value


def _fmt_pct(value: Optional[float]) -> str:
    return "not meaningful" if value is None else "%.1f%%" % (100.0 * value)


def _fmt_ci(ci: Optional[Tuple[float, float]]) -> str:
    if ci is None:
        return "undefined"
    return "[%.2f, %.2f]" % (ci[0], ci[1])


def _table(header: Sequence[str], rows: Sequence[Sequence[str]], cap: int = MAX_TABLE_ROWS) -> List[str]:
    lines = [
        "| " + " | ".join(header) + " |",
        "|" + "|".join("---" for _ in header) + "|",
    ]
    for row in rows[:cap]:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    if len(rows) > cap:
        lines.append("")
        lines.append("... %d more row(s) not shown." % (len(rows) - cap))
    if not rows:
        lines.append("| _none_ |" + " |" * (len(header) - 1))
    return lines


def _studies_phrase(n: int) -> str:
    return "%d study" % n if n == 1 else "%d studies" % n


def _census_caveat(report: "FeasibilityReport") -> List[str]:
    if not report.is_small_corpus or not report.has_usable_records:
        return []
    return [
        "> **Read the tables below as a census, not as statistics.** The corpus "
        "holds %d usable record(s) from %s, below the floors of %d records and "
        "%d studies this report requires before a proportion computed over it "
        "describes anything but the curation queue."
        % (
            report.n_records_usable,
            _studies_phrase(report.coverage.n_studies),
            SMALL_CORPUS_RECORDS,
            SMALL_CORPUS_STUDIES,
        ),
        "",
    ]


# --- sections ---------------------------------------------------------------


def _section_run(report: "FeasibilityReport") -> List[str]:
    coverage = report.coverage
    rows = [
        ["generated (UTC)", report.generated_at or "-"],
        ["input", "`%s`" % (report.source_path or "<in memory>")],
        [
            "input present",
            "yes" if report.source_exists else "no (treated as 0 records)",
        ],
        ["records read", str(coverage.n_records_read)],
        ["records usable as observations", str(coverage.n_records_usable)],
        ["records held out of the matrix", str(coverage.n_records_excluded)],
        ["usability predicate", "schema v%s" % coverage.key_version],
        [
            "claim key",
            "`%s` (schema v0.2.0, %d components)"
            % (
                " + ".join(report.replication_summary.get("claim_fields", ())),
                len(report.replication_summary.get("claim_fields", ())),
            ),
        ],
        ["distinct studies", str(coverage.n_studies)],
        ["distinct labs (last author, D-056)", str(coverage.n_labs)],
        ["usable records with no lab identity", str(coverage.n_records_without_lab)],
        [
            "usable records stating none of the %d conditions" % N_CONDITION_SLOTS,
            str(coverage.n_records_condition_unspecified),
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
    if coverage.excluded_counts:
        lines.append("### Why records were held out of the matrix")
        lines.append("")
        lines.append(
            "Counted by reason, so nothing is filtered silently. A record may "
            "appear under more than one reason. Being held out here says "
            "nothing about its admissibility to `data/curated/`, which is "
            "`syn4rch.schema.validate`'s decision."
        )
        lines.append("")
        lines.extend(
            _table(
                ["reason", "records", "meaning"],
                [
                    [code, str(count), exclusion_note(code)]
                    for code, count in coverage.excluded_counts.items()
                ],
                cap=len(EXCLUSION_CODES),
            )
        )
        lines.append("")
    return lines


def _section_populations(report: "FeasibilityReport") -> List[str]:
    coverage = report.coverage
    lines = [
        "## 1. How many neuronal populations have usable intervention data?",
        "",
        "**%d.** A population counts when at least one record names it, names a "
        "perturbation `method`, and carries a `perturbation_polarity` that is "
        "neither absent nor `unclear` — a perturbation of unstated sign states "
        "nothing an input could encode. Under schema v0.1.0 this predicate "
        "read `valence` instead, which held out every molecular record; the "
        "difference is measured in the section on the key change. Of these "
        "populations, %d are unresolved and grouped by their raw name under a "
        "`raw:` prefix (no vocabulary match; SCHEMA.md 3.2 forbids "
        "force-mapping them)."
        % (coverage.n_populations, coverage.n_populations_unresolved),
        "",
    ]
    rows = [
        [
            "`%s`" % pop,
            str(coverage.records_per_population.get(pop, 0)),
            str(coverage.studies_per_population.get(pop, 0)),
            str(sum(1 for c in coverage.row(pop) if not c.is_empty)),
        ]
        for pop in sorted(
            coverage.populations,
            key=lambda p: (-coverage.records_per_population.get(p, 0), p),
        )
    ]
    lines.extend(
        _table(["population", "records", "studies", "outcomes measured"], rows)
    )
    lines.append("")
    return lines


def _section_outcomes(report: "FeasibilityReport") -> List[str]:
    coverage = report.coverage
    lines = [
        "## 2. How many physiological outcomes have been measured?",
        "",
        "**%d.** An outcome counts when at least one record resolves it to an "
        "`outcome_id` and states a `direction` other than `not_reported`. "
        "`no_change` counts: SCHEMA.md 3.5 keeps null results, and dropping "
        "them is the leading source of direction bias."
        % coverage.n_outcomes,
        "",
    ]
    rows = [
        [
            "`%s`" % out,
            ", ".join(coverage.outcome_tiers.get(out, ())) or "-",
            str(coverage.records_per_outcome.get(out, 0)),
            str(coverage.studies_per_outcome.get(out, 0)),
            str(
                sum(
                    1
                    for pop in coverage.populations
                    if not coverage.cell(pop, out).is_empty
                )
            ),
        ]
        for out in sorted(
            coverage.outcomes,
            key=lambda o: (-coverage.records_per_outcome.get(o, 0), o),
        )
    ]
    lines.extend(
        _table(
            ["outcome", "tier", "records", "studies", "populations measured"], rows
        )
    )
    lines.append("")
    return lines


def _section_matrix(report: "FeasibilityReport") -> List[str]:
    coverage = report.coverage
    lines = [
        "## Coverage matrix (population x outcome)",
        "",
        "The matrix is built dense and the zeros are kept: which relationships "
        "nobody has measured is half the Phase 2 answer, and it disappears if "
        "only the populated cells are stored.",
        "",
    ]
    sparsity_printable = (
        coverage.n_populations >= MIN_MATRIX_DIM_FOR_SPARSITY
        and coverage.n_outcomes >= MIN_MATRIX_DIM_FOR_SPARSITY
    )
    rows = [
        ["matrix shape", "%d x %d" % (coverage.n_populations, coverage.n_outcomes)],
        ["cells", str(coverage.n_cells_total)],
        ["cells with at least one record", str(coverage.n_cells_populated)],
        ["cells with no records", str(coverage.n_cells_empty)],
        [
            "cells with two or more labs",
            str(coverage.n_cells_multi_lab),
        ],
        [
            "sparsity (fraction of cells empty)",
            _fmt_pct(coverage.sparsity)
            if sparsity_printable
            else "not meaningful at %d x %d (needs at least %d x %d)"
            % (
                coverage.n_populations,
                coverage.n_outcomes,
                MIN_MATRIX_DIM_FOR_SPARSITY,
                MIN_MATRIX_DIM_FOR_SPARSITY,
            ),
        ],
    ]
    lines.extend(_table(["field", "value"], rows))
    lines.append("")

    if (
        coverage.n_populations
        and coverage.n_outcomes
        and coverage.n_populations <= MAX_GRID_DIM
        and coverage.n_outcomes <= MAX_GRID_DIM
    ):
        lines.append("### Record counts per cell (zeros shown)")
        lines.append("")
        header = ["population \\ outcome"] + ["`%s`" % o for o in coverage.outcomes]
        grid_rows = [
            ["`%s`" % pop]
            + [str(coverage.cell(pop, out).n_records) for out in coverage.outcomes]
            for pop in coverage.populations
        ]
        lines.extend(_table(header, grid_rows, cap=MAX_GRID_DIM))
        lines.append("")
    elif coverage.n_cells_total:
        lines.append(
            "The dense grid is not drawn: %d x %d exceeds the %d x %d printable "
            "limit. %d cell(s) are empty and are listed in the JSON output by "
            "omission from `populated_cells`."
            % (
                coverage.n_populations,
                coverage.n_outcomes,
                MAX_GRID_DIM,
                MAX_GRID_DIM,
                coverage.n_cells_empty,
            )
        )
        lines.append("")

    lines.append("### Populated cells")
    lines.append("")
    lines.extend(
        _table(
            ["population", "outcome", "records", "studies", "labs"],
            [
                [
                    "`%s`" % c.population,
                    "`%s`" % c.outcome,
                    str(c.n_records),
                    str(c.n_studies),
                    str(c.n_labs),
                ]
                for c in coverage.populated_cells()
            ],
        )
    )
    lines.append("")
    return lines


def _section_replication(report: "FeasibilityReport") -> List[str]:
    summary = report.replication_summary
    fields = summary.get("claim_fields", ())
    lines = [
        "## 3. How many neuron-outcome relationships have independent "
        "replication?",
        "",
        "The unit is the claim `syn4rch.qc.keys.claim_key`, which under schema "
        "v0.2.0 is `(%s)`. SCHEMA.md 1 forbids aggregating across "
        "perturbations, so a circuit activation and a gene knockout on the "
        "same pair are two claims, not two observations of one; and the "
        "condition is in the key for the same reason one step further out, "
        "because a fasted male cohort and a fed female one were never "
        "comparable." % ", ".join(fields),
        "",
        "**%d of %d claim(s) are replicated across two or more labs**; a "
        "further %d are replicated across studies within a single lab and are "
        "discounted (SCHEMA.md 5.2 names repeat publication by one group as "
        "the largest leakage risk here); %d claim(s) rest on a single "
        "observation."
        % (
            summary.get("n_cells_replicated_across_labs", 0),
            summary.get("n_cells", 0),
            summary.get("n_cells_same_lab_only", 0),
            summary.get("n_cells_single_observation", 0),
        ),
        "",
    ]
    rows = [
        ["claims (strict v0.2.0 claim key)", summary.get("n_cells", 0)],
        ["single observation", summary.get("n_cells_single_observation", 0)],
        [
            "replicated across >= 2 studies",
            summary.get("n_cells_replicated_across_studies", 0),
        ],
        [
            "replicated across >= 2 labs",
            summary.get("n_cells_replicated_across_labs", 0),
        ],
        ["same lab only (discounted)", summary.get("n_cells_same_lab_only", 0)],
        [
            "claims containing records with no lab identity",
            summary.get("n_cells_with_unknown_lab_records", 0),
        ],
        [
            "claims matched on 'both unknown' conditions",
            summary.get("n_cells_matched_on_unknown_conditions", 0),
        ],
    ]
    lines.extend(_table(["field", "count"], [[r[0], str(r[1])] for r in rows]))
    lines.append("")
    lines.extend(_specificity_block(report))
    return lines


def _specificity_block(report: "FeasibilityReport") -> List[str]:
    """Replication counts broken down by minimum condition specificity.

    The breakdown, not a footnote under a pooled number. A replication across
    three studies that each state sex, diet and age is different evidence from
    one across three studies that state nothing, and the second matched on
    "both unknown" rather than on "both the same".
    """
    summary = report.replication_summary
    all_cells = summary.get("cells_by_min_condition_specificity", {})
    studies = summary.get(
        "replicated_across_studies_by_min_condition_specificity", {}
    )
    labs = summary.get("replicated_across_labs_by_min_condition_specificity", {})
    defined = summary.get("agreement_defined_by_min_condition_specificity", {})
    means = summary.get("mean_agreement_by_min_condition_specificity", {})

    lines = [
        "### Every replication figure, by minimum condition specificity",
        "",
        "`condition_specificity` counts how many of the %d condition slots "
        "(%s) a record actually states. "
        "The **minimum** over a cell's records is what qualifies the cell, "
        "because one record stating nothing is enough to make the whole match "
        "a match on `unknown`. A row at specificity 0 says the claim's records "
        "were grouped together partly because their conditions were not known, "
        "not because they were the same."
        % (
            N_CONDITION_SLOTS,
            ", ".join("`%s`" % f for f in CONDITION_FIELDS),
        ),
        "",
    ]
    rows = []
    for level in range(N_CONDITION_SLOTS + 1):
        key = str(level)
        mean = means.get(key)
        rows.append(
            [
                "%d of %d" % (level, N_CONDITION_SLOTS),
                str(all_cells.get(key, 0)),
                str(studies.get(key, 0)),
                str(labs.get(key, 0)),
                str(defined.get(key, 0)),
                "undefined" if mean is None else "%.2f" % mean,
            ]
        )
    lines.extend(
        _table(
            [
                "min condition specificity",
                "claims",
                "replicated >= 2 studies",
                "replicated >= 2 labs",
                "agreement defined",
                "mean agreement",
            ],
            rows,
            cap=N_CONDITION_SLOTS + 1,
        )
    )
    lines.append("")
    pooled_mean = summary.get("mean_agreement_where_defined")
    lines.append(
        "The pooled mean agreement across all specificities is %s. It is "
        "printed last and deliberately: a mean taken across cells that matched "
        "on nothing and cells that matched on three stated conditions "
        "describes neither, and the row-by-row breakdown above is the figure "
        "to read." % _fmt_ratio(pooled_mean)
    )
    lines.append("")
    if all_cells.get(str(N_CONDITION_SLOTS), 0) == 0:
        lines.append(
            "No claim in this corpus reaches specificity %d of %d: not one "
            "record states every condition. Every replication figure above "
            "is therefore qualified, and none of them describes two "
            "experiments known to have run under the same conditions."
            % (N_CONDITION_SLOTS, N_CONDITION_SLOTS)
        )
        lines.append("")
    return lines


def _section_agreement(report: "FeasibilityReport") -> List[str]:
    summary = report.replication_summary
    defined = summary.get("n_cells_agreement_defined", 0)
    lines = [
        "## 4. How often do independent studies agree on effect direction?",
        "",
        "Agreement is the fraction of independent studies stating the modal "
        "direction, with a %d%% Wilson score interval (chosen because the Wald "
        "interval has width zero at a unanimous 2 of 2, and the standard "
        "library has no incomplete beta for Clopper-Pearson; "
        "`syn4rch.feasibility.replication` documents this)."
        % int(round(100 * CONFIDENCE_LEVEL)),
        "",
        "**Agreement is undefined, not 100%%, for a claim with one independent "
        "observation**: %d of %d claim(s) have no second independent study, "
        "and a defined agreement exists for only %d."
        % (
            summary.get("n_cells_agreement_undefined", 0),
            summary.get("n_cells", 0),
            defined,
        ),
        "",
    ]
    if defined:
        lines.append(
            "Of the %d claim(s) with a defined agreement, %d are unanimous and "
            "%d are contested (independent studies state incompatible "
            "directions, by the QC definition in docs/DECISIONS.md D-011). "
            "%d of them matched on conditions that were unknown rather than "
            "known to be the same; the mean agreement is broken down by "
            "specificity in section 3 rather than pooled here. Nothing was "
            "resolved, averaged or dropped."
            % (
                defined,
                summary.get("n_cells_unanimous_where_defined", 0),
                summary.get("n_cells_contested", 0),
                summary.get("n_cells_matched_on_unknown_conditions", 0),
            )
        )
        lines.append("")
        replicated = [c for c in report.cells if c.agreement_is_defined]
        lines.append(
            "The `min cond` column is the minimum condition specificity of the "
            "cell, out of %d. It is not a footnote: an agreement at `0/%d` was "
            "computed over studies that matched on their conditions being "
            "unknown." % (N_CONDITION_SLOTS, N_CONDITION_SLOTS)
        )
        lines.append("")
        lines.extend(
            _table(
                [
                    "population",
                    "perturbation (class / polarity / target)",
                    "outcome",
                    "condition",
                    "min cond",
                    "studies",
                    "labs",
                    "modal direction",
                    "agreement",
                    "%d%% CI" % int(round(100 * CONFIDENCE_LEVEL)),
                    "conflicts",
                ],
                [
                    [
                        "`%s`" % c.population,
                        c.perturbation_label,
                        "`%s`" % c.outcome,
                        c.condition_label,
                        "%d/%d" % (c.min_condition_specificity, N_CONDITION_SLOTS),
                        str(c.n_studies),
                        str(c.n_labs),
                        (c.modal_direction or "-")
                        + (" (tied)" if c.modal_direction_is_tied else ""),
                        _fmt_ratio(c.agreement_studies_value),
                        _fmt_ci(c.agreement_studies_ci),
                        ", ".join(c.conflict_codes) or "-",
                    ]
                    for c in sorted(
                        replicated,
                        key=lambda c: (
                            -c.n_labs,
                            -c.n_studies,
                            -c.min_condition_specificity,
                            c.population,
                            c.outcome,
                        ),
                    )
                ],
            )
        )
        lines.append("")
    else:
        lines.append(
            "No claim in this corpus has two independent studies stating a "
            "direction, so there is no agreement rate to report. This is the "
            "measured answer to question 4, not a gap in the computation."
        )
        lines.append("")

    directions = summary.get("modal_direction_counts", {})
    if directions:
        lines.append("### Modal direction per claim")
        lines.append("")
        lines.extend(
            _table(
                ["modal direction", "claims"],
                [[k, str(v)] for k, v in directions.items()],
            )
        )
        lines.append("")
    return lines


def _section_gates(report: "FeasibilityReport") -> List[str]:
    summary = report.verdict_summary
    by_verdict = summary.get("by_verdict", {})
    lines = [
        "## 5. Which outcomes have sufficient evidence for predictive "
        "modelling?",
        "",
        "**%s**"
        % (
            "None."
            if not summary.get("any_feasible")
            else "Feasible for at least one population and valence: "
            + ", ".join("`%s`" % o for o in summary.get("outcomes_feasible", []))
            + "."
        ),
        "",
        "Verdicts are per claim-keyed task (the full claim key) and then "
        "rolled up to the `(population, outcome)` pairs the "
        "Phase 2 question is phrased in; a pair is as feasible as its best "
        "claim and the rollup names which perturbation and which condition "
        "that was. No threshold in `gates.py` changed at v0.2.0 — only the "
        "unit they apply to, which is finer, so fewer cells can clear a gate "
        "that counts studies and labs. The section on the key change measures "
        "that cost.",
        "",
    ]
    lines.extend(
        _table(
            ["verdict", "tasks", "(population, outcome) pairs"],
            [
                [
                    "feasible",
                    str(by_verdict.get("feasible", 0)),
                    str(report.n_feasible),
                ],
                [
                    "marginal",
                    str(by_verdict.get("marginal", 0)),
                    str(report.n_marginal),
                ],
                [
                    "insufficient",
                    str(by_verdict.get("insufficient", 0)),
                    str(report.n_insufficient),
                ],
            ],
        )
    )
    lines.append("")

    promoted = [
        v
        for v in report.verdicts
        if v.verdict in (gates_mod.VERDICT_FEASIBLE, gates_mod.VERDICT_MARGINAL)
    ]
    if promoted:
        lines.append("### Tasks reaching marginal or better")
        lines.append("")
        lines.extend(
            _table(
                [
                    "population",
                    "perturbation",
                    "outcome",
                    "condition",
                    "min cond",
                    "verdict",
                    "records",
                    "studies",
                    "labs",
                    "agreement",
                    "reasons",
                    "notes",
                ],
                [
                    [
                        "`%s`" % v.population,
                        v.perturbation_label,
                        "`%s`" % v.outcome,
                        v.condition_label,
                        "%d/%d" % (v.min_condition_specificity, N_CONDITION_SLOTS),
                        v.verdict,
                        str(v.numbers.get("n_records", 0)),
                        str(v.numbers.get("n_studies", 0)),
                        str(v.numbers.get("n_labs", 0)),
                        _fmt_ratio(v.numbers.get("agreement_studies")),
                        ", ".join("`%s`" % r for r in v.reasons),
                        ", ".join("`%s`" % n for n in v.notes) or "-",
                    ]
                    for v in sorted(
                        promoted,
                        key=lambda v: (
                            -v.rank,
                            -v.numbers.get("n_labs", 0),
                            -v.min_condition_specificity,
                            v.population,
                            v.outcome,
                        ),
                    )
                ],
            )
        )
        lines.append("")
        summary_spec = report.verdict_summary.get(
            "by_verdict_and_min_condition_specificity", {}
        )
        if summary_spec:
            lines.append("### Verdicts by minimum condition specificity")
            lines.append("")
            lines.append(
                "A `feasible` verdict whose cell matched on unknown conditions "
                "is still a `feasible` verdict — there is no specificity gate, "
                "because inventing a threshold here would be inventing a knob "
                "(`MIN_CONDITION_SPECIFICITY_GATE` is `None` on purpose). It "
                "is counted separately instead, so the qualification cannot be "
                "lost between this table and the headline count."
            )
            lines.append("")
            lines.extend(
                _table(
                    ["verdict"]
                    + [
                        "%d/%d" % (k, N_CONDITION_SLOTS)
                        for k in range(N_CONDITION_SLOTS + 1)
                    ],
                    [
                        [name]
                        + [
                            str(summary_spec.get(name, {}).get(str(k), 0))
                            for k in range(N_CONDITION_SLOTS + 1)
                        ]
                        for name in reversed(gates_mod.VERDICTS)
                    ],
                    cap=len(gates_mod.VERDICTS),
                )
            )
            lines.append("")
    else:
        lines.append(
            "No task reaches `marginal`. Every task in this corpus is "
            "`insufficient`, which is a legitimate Phase 2 conclusion: it says "
            "the corpus does not yet support a prediction task, and it is "
            "reached from counts rather than from an opinion about the "
            "literature."
        )
        lines.append("")

    by_reason = summary.get("by_reason", {})
    if by_reason:
        lines.append("### Why tasks did not pass")
        lines.append("")
        lines.extend(
            _table(
                ["reason", "tasks", "meaning"],
                [
                    [
                        "`%s`" % code,
                        str(count),
                        gates_mod.REASON_NOTES.get(code, "-"),
                    ]
                    for code, count in by_reason.items()
                ],
                cap=len(gates_mod.REASON_NOTES) + 1,
            )
        )
        lines.append("")

    lines.append("### Promotion gates")
    lines.append("")
    lines.append(
        "Every threshold is a named constant in "
        "`src/syn4rch/feasibility/gates.py`, each with the reason that number "
        "was chosen. They are deliberately pessimistic: lab-level "
        "independence is required for `feasible`, agreement at one "
        "observation is undefined and undefined never passes, and the Wilson "
        "lower bound must clear the chance rate, so unanimity among very few "
        "studies lands at `marginal`."
    )
    lines.append("")
    lines.extend(
        _table(
            ["gate", "value"],
            [
                [k, ("%.4g" % v) if isinstance(v, float) else str(v)]
                for k, v in gates_mod.THRESHOLDS.items()
            ],
            cap=len(gates_mod.THRESHOLDS),
        )
    )
    lines.append("")
    return lines


def _section_key_change(report: "FeasibilityReport") -> List[str]:
    """The before and after of the schema v0.2.0 claim key change.

    Kept as its own section, and never reduced to the "after" column, because
    the difference between the two keys *is* the measurement of how much
    apparent disagreement in this corpus was condition and perturbation
    heterogeneity. That is why ``syn4rch.qc.keys.claim_key_pooled`` still
    exists.
    """
    comparison = report.key_comparison
    usability = report.usability_comparison
    lines = [
        "## What schema v0.2.0 changed, measured both ways",
        "",
        "`claim_key` was `(population, valence, outcome_id)`. It is now "
        "`(population, perturbation_class, perturbation_polarity, "
        "perturbation_target, outcome_id, sex, diet, age_group, challenge)`. "
        "`syn4rch.qc.keys.claim_key_pooled` preserves the old key exactly, for "
        "one purpose: without it this change would be unfalsifiable. Both "
        "columns below are computed over the same usable records, so the only "
        "thing varying between them is the key.",
        "",
    ]
    if comparison is None:
        lines.append("_No records, so there is nothing to compare._")
        lines.append("")
        return lines

    pooled = comparison.pooled
    strict = comparison.strict
    rows = [
        ["records keyed", str(pooled.n_records), str(strict.n_records)],
        ["cells", str(pooled.n_cells), str(strict.n_cells)],
        [
            "directly opposed cells (increase and decrease in one cell)",
            str(pooled.n_cells_directly_opposed),
            str(strict.n_cells_directly_opposed),
        ],
        [
            "cells carrying a QC conflict code across studies",
            str(pooled.n_cells_with_conflict_codes),
            str(strict.n_cells_with_conflict_codes),
        ],
        [
            "multi-study cells",
            str(pooled.n_cells_multi_study),
            str(strict.n_cells_multi_study),
        ],
        [
            "multi-lab cells",
            str(pooled.n_cells_multi_lab),
            str(strict.n_cells_multi_lab),
        ],
        [
            "cells with a defined agreement",
            str(pooled.n_cells_agreement_defined),
            str(strict.n_cells_agreement_defined),
        ],
        [
            "single-observation cells",
            str(pooled.n_cells_single_observation),
            str(strict.n_cells_single_observation),
        ],
    ]
    lines.extend(
        _table(
            ["quantity", "pooled key (v0.1.0)", "strict key (v0.2.0)"],
            rows,
            cap=len(rows),
        )
    )
    lines.append("")
    lines.append("### What the split did, cell by cell")
    lines.append("")
    rows = [
        [
            "pooled cells the strict key splits in two or more",
            str(comparison.n_pooled_cells_split),
        ],
        [
            "largest number of strict cells out of one pooled cell",
            str(comparison.max_strict_cells_per_pooled_cell),
        ],
        [
            "opposed pooled cells whose refinements are all internally "
            "consistent (the contradiction was heterogeneity)",
            str(comparison.n_opposed_cells_resolved_by_split),
        ],
        [
            "opposed pooled cells still opposed after the split (a real "
            "head-to-head disagreement)",
            str(comparison.n_opposed_cells_surviving),
        ],
        [
            "multi-study pooled cells whose every refinement has one study "
            "(replication evidence given up)",
            str(comparison.n_multi_study_cells_lost),
        ],
        [
            "— of those, lost because the refinements are different "
            "perturbations",
            str(comparison.n_lost_to_perturbation_split),
        ],
        [
            "— of those, lost because two refinements state *different* "
            "conditions",
            str(comparison.n_lost_to_stated_condition_difference),
        ],
        [
            "— of those, lost **only** because one paper states a condition "
            "the other leaves unknown",
            str(comparison.n_lost_to_unknown_against_stated),
        ],
    ]
    lines.extend(_table(["quantity", "count"], rows, cap=len(rows)))
    lines.append("")
    bought = comparison.n_opposed_cells_resolved_by_split
    cost = comparison.n_multi_study_cells_lost
    lines.append(
        "**What it bought:** %d of the %d head-to-head `increase` vs "
        "`decrease` cells under the old key were not contradictions. A finer "
        "key separated them, which means they had been a male cohort against a "
        "female one, a chow-fed cohort against a high-fat one, or a gene "
        "knockout against a gene overexpression — experiments that were never "
        "comparable, recorded as laboratories disagreeing. %d remain, and "
        "those are the real disagreements."
        % (bought, pooled.n_cells_directly_opposed, comparison.n_opposed_cells_surviving)
    )
    lines.append("")
    lines.append(
        "**What it cost:** %d of the %d multi-study cells stopped being "
        "multi-study — %s of the replication evidence the old key appeared to "
        "have. Most of it was never evidence: %d of the %d were records "
        "describing different *perturbations* (a knockout of one gene against "
        "a knockout of another), and %d were two different stated conditions. "
        "Pooling those and counting the papers is the defect v0.2.0 exists to "
        "fix. The loss is still real in the sense that matters for modelling: "
        "there are now fewer cells that can clear a gate counting studies and "
        "labs."
        % (
            cost,
            pooled.n_cells_multi_study,
            _fmt_pct(
                None
                if not pooled.n_cells_multi_study
                else cost / pooled.n_cells_multi_study
            ),
            comparison.n_lost_to_perturbation_split,
            cost,
            comparison.n_lost_to_stated_condition_difference,
        )
    )
    lines.append("")
    if comparison.cost_that_cannot_be_adjudicated:
        lines.append(
            "**What cannot be adjudicated:** the remaining %d lost cell(s) "
            "were split *only* because one paper states a condition the other "
            "leaves unknown. `unknown` pools \"the curator did not ask\" with "
            "\"the paper does not say\", and neither can match a stated value: "
            "a paper reporting male mice and a paper reporting no sex at all "
            "are two keys, so they stop replicating each other even if both "
            "ran males. **Nothing in this corpus can say whether they did.** "
            "The alternative — assuming the unstated condition matched — is "
            "the inference SCHEMA.md 2 rule 4 forbids, so this cost is paid "
            "and reported rather than argued away. It is counted on its own "
            "line and never as heterogeneity resolved."
            % comparison.cost_that_cannot_be_adjudicated
        )
        lines.append("")

    pooled_summary = report.pooled_verdict_summary
    if pooled_summary:
        lines.append("### Gate arithmetic under both keys")
        lines.append("")
        lines.append(
            "The same gates, the same records, two keys. The pooled column is "
            "**not a result** — a pooled cell may mix two perturbations or two "
            "cohorts, so a verdict over it describes no single experiment. It "
            "is here so the cost of the correct key is visible rather than "
            "absorbed."
        )
        lines.append("")
        pooled_by = pooled_summary.get("by_verdict", {})
        strict_by = report.verdict_summary.get("by_verdict", {})
        pooled_pair_counts = _pair_verdict_counts(report.pooled_pairs)
        strict_pair_counts = _pair_verdict_counts(report.pairs)
        lines.extend(
            _table(
                [
                    "verdict",
                    "tasks (pooled v0.1.0)",
                    "tasks (strict v0.2.0)",
                    "pairs (pooled)",
                    "pairs (strict)",
                ],
                [
                    [
                        name,
                        str(pooled_by.get(name, 0)),
                        str(strict_by.get(name, 0)),
                        str(pooled_pair_counts.get(name, 0)),
                        str(strict_pair_counts.get(name, 0)),
                    ]
                    for name in reversed(gates_mod.VERDICTS)
                ],
                cap=len(gates_mod.VERDICTS),
            )
        )
        lines.append("")

    if usability is not None:
        lines.append("### The usability predicate also changed")
        lines.append("")
        lines.append(
            "Usability used to require a `valence` that was neither absent nor "
            "`unclear`. Under v0.2.0 it requires a `perturbation_polarity` that "
            "is neither absent nor `unclear`. This is a corrected predicate, "
            "not a loosened threshold: a conditional knockout's `valence` is "
            "honestly `unclear` (SCHEMA.md 3.4 and 6 forbid coercing it) while "
            "its *manipulation* has a perfectly well defined sign, which is "
            "`decrease` on the amount of the target gene product. The old "
            "predicate tested a field that had stopped carrying the "
            "information. Every number in `gates.py` is unchanged."
        )
        lines.append("")
        rows = [
            ["records read", str(usability.n_records)],
            ["usable under the v0.1.0 predicate", str(usability.n_usable_v1)],
            ["usable under the v0.2.0 predicate", str(usability.n_usable_v2)],
            [
                "admitted by v0.2.0 only",
                str(usability.n_admitted_by_v2_only),
            ],
            [
                "admitted by v0.1.0 only",
                str(usability.n_admitted_by_v1_only),
            ],
        ]
        lines.extend(_table(["quantity", "count"], rows, cap=len(rows)))
        lines.append("")
        if usability.admitted_by_v2_only_by_class:
            lines.append(
                "The records the change admits, by perturbation class — a "
                "predicate change that *admits* records is the kind that most "
                "needs showing, so it is broken down rather than stated as one "
                "number:"
            )
            lines.append("")
            lines.extend(
                _table(
                    ["perturbation class", "records admitted by v0.2.0 only"],
                    [
                        ["`%s`" % cls, str(count)]
                        for cls, count in usability.admitted_by_v2_only_by_class.items()
                    ],
                    cap=len(usability.admitted_by_v2_only_by_class) + 1,
                )
            )
            lines.append("")
    return lines


def _pair_verdict_counts(pairs: Sequence[PairVerdict]) -> Dict[str, int]:
    counts: Dict[str, int] = {name: 0 for name in gates_mod.VERDICTS}
    for pair in pairs:
        counts[pair.verdict] = counts.get(pair.verdict, 0) + 1
    return counts


def _section_splits(report: "FeasibilityReport") -> List[str]:
    lines = [
        "## Leakage-aware splits available on this corpus",
        "",
        "SCHEMA.md 5: a cross-validated number reported without naming its "
        "grouping is an invalid result in this project. "
        "`syn4rch.feasibility.splits` therefore returns folds only inside a "
        "labelled fold set, and `attach_metric` raises if a caller tries to "
        "bind a number to folds without the label. Shared animal cohorts "
        "(SCHEMA.md 5 item 3) are not controlled by any grouping here: the "
        "schema records no cohort identity.",
        "",
        "**The claim key change does not reach this table.** All three fold "
        "generators key on `study_key`, `lab_key` and `population_key`; none "
        "of them reads `claim_key`, a perturbation field or a condition field, "
        "so fold membership is identical before and after schema v0.2.0. That "
        "is checked rather than assumed: "
        "`tests/test_feasibility_splits.py` holds records that differ only in "
        "`perturbation_class`, `perturbation_polarity`, `target_id` and the "
        "every condition, and asserts the folds come out unchanged.",
        "",
        "Grouping is also not the whole of leakage control. SCHEMA.md 5b: a "
        "split can keep a lab off both sides and still hand the model a test "
        "row whose perturbation identity never appeared in training, which is "
        "unlearnable rather than merely hard. That fraction is measured in "
        "`reports/perturbation_identity.md` and is required next to any "
        "cross-validated number from this corpus.",
        "",
    ]
    rows = []
    for name, fold_set in report.fold_sets.items():
        rows.append(
            [
                "`%s`" % name,
                str(fold_set.n_groups),
                str(fold_set.n_usable_folds),
                "yes" if fold_set.is_usable else "no",
                str(fold_set.n_ungroupable),
            ]
        )
    lines.extend(
        _table(
            [
                "grouping",
                "groups",
                "usable folds",
                "cross-validation possible",
                "records excluded as ungroupable",
            ],
            rows,
            cap=len(report.fold_sets) + 1,
        )
    )
    lines.append("")
    unusable = [n for n, fs in report.fold_sets.items() if not fs.is_usable]
    if unusable:
        lines.append(
            "Cross-validation is impossible under %s on this corpus: fewer "
            "than two folds have both a train and a test side."
            % ", ".join("`%s`" % n for n in unusable)
        )
        lines.append("")
    return lines


def _section_limits(report: "FeasibilityReport") -> List[str]:
    return [
        "## What this report does not claim",
        "",
        "- It does not estimate how many records curation will eventually "
        "yield. Every number above is a count of what is in the input file at "
        "the time shown.",
        "- It does not resolve, average, prefer or drop either side of a "
        "disagreement. Contested claims are reported as contested.",
        "- It does not treat a single observation as agreement, and it does "
        "not treat repeated publication by one lab as independent replication.",
        "- It cannot rule out two studies sharing one animal cohort "
        "(SCHEMA.md 5 item 3): the schema records no cohort identity.",
        "- A `feasible` verdict says a prediction task can be posed and "
        "cross-validated with the grouping named, not that a model will work. "
        "It does not say the held-out rows were learnable: see "
        "`reports/perturbation_identity.md` (SCHEMA.md 5b).",
        "- It does not claim that two records in one cell ran under the same "
        "conditions. Where the minimum condition specificity is below %d, they "
        "were grouped partly because the conditions were unknown. No claim in "
        "this corpus states all %d."
        % (N_CONDITION_SLOTS, N_CONDITION_SLOTS),
        "- The pooled (v0.1.0) key appears only in the section on the key "
        "change. It is never a verdict: it pools experiments that were never "
        "comparable, which is the defect v0.2.0 exists to fix.",
        "- Tier B outcomes are excluded as prediction targets by "
        "`EXCLUDE_TIER_B_AS_TARGET` (SCHEMA.md 6: Tier B is an annotation "
        "layer, not a training target).",
        "",
    ]


def _section_nothing_to_measure(report: "FeasibilityReport") -> List[str]:
    """Replaces sections 1-5 when no record is usable.

    Printing "0 of 0 claims agree" and a sparsity over a 0x0 matrix would be
    arithmetic performed on nothing. The five questions are answered, in
    words, with the single fact that decides all five.
    """
    if report.n_records_read == 0:
        cause = (
            "the input file does not exist yet"
            if not report.source_exists
            else "the input file is empty"
        )
    else:
        cause = (
            "all %d record(s) read were held out of the matrix (see the Run "
            "table for the reasons)" % report.n_records_read
        )
    return [
        "## The five Phase 2 questions",
        "",
        "**None of the five can be answered: %s, so there are 0 usable "
        "observations.** No table below is computed over more than zero "
        "records, and no derived statistic is printed for them: a sparsity, an "
        "agreement rate or a mean over an empty corpus would be arithmetic on "
        "nothing." % cause,
        "",
        *_table(
            ["question", "measured answer"],
            [
                ["1. populations with usable intervention data", "0"],
                ["2. physiological outcomes measured", "0"],
                ["3. relationships with independent replication", "0"],
                [
                    "4. direction agreement among independent studies",
                    "undefined (no claim has two independent studies)",
                ],
                [
                    "5. outcomes with sufficient evidence for modelling",
                    "none",
                ],
            ],
            cap=5,
        ),
        "",
        "The promotion gates these answers were judged against are listed "
        "below, so the verdict can be re-checked when the corpus grows.",
        "",
        *_table(
            ["gate", "value"],
            [
                [k, ("%.4g" % v) if isinstance(v, float) else str(v)]
                for k, v in gates_mod.THRESHOLDS.items()
            ],
            cap=len(gates_mod.THRESHOLDS),
        ),
        "",
    ]


def render_markdown(report: "FeasibilityReport") -> str:
    """The full report. First line after the title is the verdict sentence."""
    lines: List[str] = [
        "# SYN4RCH feasibility report (Phase 2)",
        "",
        top_line(report),
        "",
        "Every number below is counted from the input named in the Run table. "
        "Nothing is estimated, imputed or projected. An `insufficient` verdict "
        "is a measurement of the corpus, not a failure of this module.",
        "",
    ]
    lines.extend(_census_caveat(report))
    lines.extend(_section_run(report))
    if report.has_usable_records:
        lines.extend(_section_populations(report))
        lines.extend(_section_outcomes(report))
        lines.extend(_section_replication(report))
        lines.extend(_section_agreement(report))
        lines.extend(_section_gates(report))
        lines.extend(_section_key_change(report))
        lines.extend(_section_matrix(report))
        lines.extend(_section_splits(report))
    else:
        lines.extend(_section_nothing_to_measure(report))
    lines.extend(_section_limits(report))
    text = "\n".join(lines).rstrip() + "\n"
    return text
