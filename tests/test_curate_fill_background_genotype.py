"""The background-genotype fill (D-060, applied in D-062).

This writer fills a dimension that is null everywhere and must never overwrite
one that is not (D-039), so most of what follows tests its refusals. The grammar
tests matter as much as the behaviour tests: a value that is nearly right is the
one that enters a claim key and is never looked at again.
"""

import _path  # noqa: F401
import io
import json
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

import synthetic
from syn4rch.curate.fill_background_genotype import (
    DIMENSION,
    EFFECTS,
    OUTCOMES,
    SCOPES,
    SENTINELS,
    Assignment,
    GrammarError,
    IdentityMovedError,
    apply_assignments,
    build_parser,
    fill_rule,
    main,
    parse_value,
    plan_assignments,
    read_table,
    render_report_md,
)
from syn4rch.schema.record import RECORD_ID_FIELDS, TierARecord
from syn4rch.schema.validate import Severity, validate

GENES = ["asip", "lep", "npy", "slc17a6"]
QUOTE = "Synthetic quote: these animals additionally lacked Npy, by 25%."


def record(**overrides):
    return synthetic.good_record(support_quote=QUOTE, **overrides)


def assignment(record_id, value="npy:loss_of_function:global", **kw):
    data = dict(
        record_id=record_id,
        study_id="synthetic:test-001",
        value=value,
        quote=QUOTE if value not in ("not_reported",) else "",
        locator="Results para 1" if value not in ("not_reported",) else "",
        source="passage",
        confidence="high",
        reasoning="synthetic",
    )
    data.update(kw)
    return Assignment(**data)


class GrammarTest(unittest.TestCase):
    def test_the_three_sentinels_parse_to_no_segments(self):
        for value in SENTINELS:
            self.assertEqual(parse_value(value, GENES), [])

    def test_a_well_formed_segment_parses(self):
        self.assertEqual(
            parse_value("npy:loss_of_function:global", GENES),
            [("npy", "loss_of_function", "global")],
        )

    def test_null_is_refused_because_it_is_the_absence_of_a_row(self):
        with self.assertRaises(GrammarError) as c:
            parse_value("null", GENES)
        self.assertIn("omit the record", str(c.exception))

    def test_an_unknown_gene_is_refused_not_prefixed(self):
        """Unlike population_key's raw: fallback: the vocabulary was extended
        for exactly the genes this pass needs, so an unknown one means the table
        and the vocabulary disagree."""
        with self.assertRaises(GrammarError) as c:
            parse_value("notagene:loss_of_function:global", GENES)
        self.assertIn("does not resolve", str(c.exception))

    def test_an_unknown_effect_or_scope_is_refused(self):
        for bad in ("npy:nonsense:global", "npy:loss_of_function:nowhere"):
            with self.assertRaises(GrammarError):
                parse_value(bad, GENES)

    def test_an_uppercase_gene_is_refused(self):
        with self.assertRaises(GrammarError):
            parse_value("Npy:loss_of_function:global", GENES)

    def test_a_malformed_segment_is_refused(self):
        for bad in ("npy", "npy:global", "npy:loss_of_function:global:extra", ""):
            with self.assertRaises(GrammarError):
                parse_value(bad, GENES)

    def test_segments_must_be_sorted_and_deduplicated(self):
        ok = "lep:loss_of_function:global+npy:loss_of_function:global"
        self.assertEqual(len(parse_value(ok, GENES)), 2)
        with self.assertRaises(GrammarError) as c:
            parse_value("npy:loss_of_function:global+lep:loss_of_function:global", GENES)
        self.assertIn("sorted", str(c.exception))
        with self.assertRaises(GrammarError):
            parse_value("npy:loss_of_function:global+npy:loss_of_function:global", GENES)

    def test_the_declared_enums_are_what_the_grammar_accepts(self):
        for effect in EFFECTS:
            parse_value("npy:%s:global" % effect, GENES)
        for scope in SCOPES:
            parse_value("npy:loss_of_function:%s" % scope, GENES)


class RefusalTest(unittest.TestCase):
    def test_a_missing_record_is_refused(self):
        got = plan_assignments([record()], [assignment("not-an-id")], GENES)
        self.assertEqual([o.outcome for o in got], ["refused_record_not_found"])

    def test_an_existing_value_is_never_overwritten(self):
        """D-039. This is a FILL; overwriting is correct_conditions' job."""
        rec = record(background_genotype="wild_type")
        got = plan_assignments([rec], [assignment(rec.record_id)], GENES)
        self.assertEqual(got[0].outcome, "refused_value_already_set")
        self.assertEqual(got[0].stored_before, "wild_type")

    def test_a_positive_value_without_a_quote_is_refused(self):
        rec = record()
        got = plan_assignments(
            [rec], [assignment(rec.record_id, quote="", locator="")], GENES
        )
        self.assertEqual(got[0].outcome, "refused_missing_quote")

    def test_not_reported_carrying_a_quote_is_refused(self):
        """An absence cannot be quoted (A009, A011), and the pass refuses the
        row rather than silently dropping the quote. All 26 not_reported rows
        of the real pass arrived like this, carrying ARM evidence."""
        rec = record()
        got = plan_assignments(
            [rec],
            [assignment(rec.record_id, value="not_reported", quote=QUOTE,
                        locator="Results para 1")],
            GENES,
        )
        self.assertEqual(got[0].outcome, "refused_quote_on_an_absence")

    def test_an_unknown_gene_in_a_row_is_its_own_outcome(self):
        rec = record()
        got = plan_assignments(
            [rec], [assignment(rec.record_id, value="zzz:loss_of_function:global")], GENES
        )
        self.assertEqual(got[0].outcome, "refused_unknown_gene")

    def test_two_rows_for_one_record_are_refused(self):
        rec = record()
        got = plan_assignments(
            [rec],
            [assignment(rec.record_id), assignment(rec.record_id, value="wild_type")],
            GENES,
        )
        self.assertEqual([o.outcome for o in got], ["filled", "refused_duplicate_row"])

    def test_every_outcome_name_is_declared(self):
        rec = record()
        for o in plan_assignments([rec], [assignment(rec.record_id)], GENES):
            self.assertIn(o.outcome, OUTCOMES)


class ApplyTest(unittest.TestCase):
    def test_the_value_and_an_evidence_entry_are_both_written(self):
        rec = record()
        result = apply_assignments(
            [rec], [assignment(rec.record_id)], known_genes=GENES
        )
        after = result.records[0]
        self.assertEqual(after.background_genotype, "npy:loss_of_function:global")
        entry = [e for e in after.condition_evidence if e.dimension == DIMENSION]
        self.assertEqual(len(entry), 1)
        self.assertEqual(entry[0].value, "npy:loss_of_function:global")
        self.assertEqual(entry[0].quote, QUOTE)
        self.assertTrue(entry[0].rule.startswith("read_under_the_arm_rule"))
        self.assertEqual(result.evidence_added, 1)

    def test_not_reported_gets_no_evidence_entry(self):
        rec = record()
        result = apply_assignments(
            [rec],
            [assignment(rec.record_id, value="not_reported", quote="", locator="")],
            known_genes=GENES,
        )
        after = result.records[0]
        self.assertEqual(after.background_genotype, "not_reported")
        self.assertEqual(
            [e for e in after.condition_evidence if e.dimension == DIMENSION], []
        )
        self.assertEqual(result.evidence_added, 0)

    def test_a_filled_record_earns_no_hard_finding(self):
        rec = record()
        result = apply_assignments(
            [rec], [assignment(rec.record_id)], known_genes=GENES
        )
        hard = [f for f in validate(result.records[0]) if f.severity is Severity.HARD]
        self.assertEqual(hard, [])

    def test_wild_type_and_mixed_are_positive_and_carry_evidence(self):
        for value in ("wild_type", "mixed"):
            rec = record()
            result = apply_assignments(
                [rec], [assignment(rec.record_id, value=value)], known_genes=GENES
            )
            entry = [
                e for e in result.records[0].condition_evidence
                if e.dimension == DIMENSION
            ]
            self.assertEqual(len(entry), 1, value)
            self.assertEqual(entry[0].value, value)

    def test_no_record_id_moves(self):
        rec = record()
        before = rec.record_id
        result = apply_assignments(
            [rec], [assignment(rec.record_id)], known_genes=GENES
        )
        self.assertEqual(result.records[0].record_id, before)
        self.assertEqual(result.ids_unchanged, 1)
        self.assertIs(result.to_dict()["record_id_map_required"], False)

    def test_the_dimension_is_not_a_hash_component(self):
        """If this ever fails, the no-id-map property of D-060 is void and this
        pass must grow a map before it may write."""
        self.assertNotIn(DIMENSION, RECORD_ID_FIELDS)

    def test_a_moved_identity_refuses(self):
        import syn4rch.curate.fill_background_genotype as M

        class Moved(TierARecord):
            @classmethod
            def from_dict(cls, data):
                r = TierARecord.from_dict(data)
                if data.get(DIMENSION):
                    r.record_id = "0" * 16
                return r

        rec = record()
        real = M.TierARecord
        try:
            M.TierARecord = Moved
            with self.assertRaises(IdentityMovedError):
                apply_assignments([rec], [assignment(rec.record_id)], known_genes=GENES)
        finally:
            M.TierARecord = real

    def test_the_input_records_are_not_mutated(self):
        rec = record()
        snapshot = rec.to_dict()
        apply_assignments([rec], [assignment(rec.record_id)], known_genes=GENES)
        self.assertEqual(rec.to_dict(), snapshot)

    def test_running_it_twice_reports_already_filled(self):
        rec = record()
        once = apply_assignments([rec], [assignment(rec.record_id)], known_genes=GENES)
        twice = apply_assignments(
            once.records, [assignment(rec.record_id)], known_genes=GENES
        )
        self.assertEqual([o.outcome for o in twice.outcomes], ["already_filled"])
        self.assertEqual(twice.evidence_added, 0)

    def test_the_rule_names_the_decision(self):
        self.assertEqual(fill_rule("D-062"), "read_under_the_arm_rule_d062")


class TableAndReportTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_the_committed_table_reads_and_is_internally_consistent(self):
        rows, meta = read_table("data/curated/background_genotype_assignments.json")
        self.assertTrue(rows)
        self.assertEqual(meta.get("decision"), "D-062")
        for row in rows:
            if row.value == "not_reported":
                self.assertEqual(row.quote, "", row.record_id)
            else:
                self.assertTrue(row.quote.strip(), row.record_id)
                self.assertTrue(row.locator.strip(), row.record_id)

    def test_every_committed_value_parses(self):
        rows, _ = read_table("data/curated/background_genotype_assignments.json")
        genes = sorted(json.load(open("data/vocab/targets.json"))["entries"])
        for row in rows:
            parse_value(row.value, genes)

    def test_the_report_explains_the_three_sentinels(self):
        rec = record()
        before = [TierARecord.from_dict(rec.to_dict())]
        result = apply_assignments(
            [rec], [assignment(rec.record_id)], known_genes=GENES
        )
        text = render_report_md(
            result, before, generated="2026-01-01T00:00:00Z",
            records_path="r.jsonl", table_path="t.json",
        )
        self.assertIn("not_reported", text)
        self.assertIn("wild_type", text)
        self.assertIn("mixed", text)
        self.assertIn("cannot be quoted", text)
        self.assertIn("not a gap to be closed", text)


class CLITest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.records = os.path.join(self.tmp, "records.jsonl")
        self.table = os.path.join(self.tmp, "table.json")
        self.report = os.path.join(self.tmp, "report.md")
        self.rec = record()
        from syn4rch.schema.io import write_jsonl

        write_jsonl(self.records, [self.rec])
        json.dump(
            {
                "_meta": {"decision": "D-062"},
                "assignments": [assignment(self.rec.record_id).to_dict()],
            },
            open(self.table, "w"),
        )

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = main(argv)
        return code, out.getvalue(), err.getvalue()

    def test_dry_run_writes_nothing(self):
        unchanged = open(self.records, encoding="utf-8").read()
        code, out, err = self._run(
            [self.records, "--table", self.table, "--report", self.report,
             "--fulltext-dir", self.tmp, "--dry-run"]
        )
        self.assertEqual(code, 0, err)
        self.assertIn("dry run", out)
        self.assertEqual(open(self.records, encoding="utf-8").read(), unchanged)

    def test_it_writes_the_dataset_a_backup_and_a_report(self):
        code, _, err = self._run(
            [self.records, "--table", self.table, "--report", self.report,
             "--fulltext-dir", self.tmp, "--quiet"]
        )
        self.assertEqual(code, 0, err)
        rows = [json.loads(l) for l in open(self.records, encoding="utf-8") if l.strip()]
        self.assertEqual(rows[0][DIMENSION], "npy:loss_of_function:global")
        self.assertTrue(os.path.exists(self.records + ".pre-background-genotype.bak"))
        self.assertIn("Background genotype fill", open(self.report, encoding="utf-8").read())

    def test_a_refusal_blocks_the_write_unless_allowed(self):
        json.dump(
            {
                "_meta": {},
                "assignments": [assignment("not-an-id").to_dict()],
            },
            open(self.table, "w"),
        )
        unchanged = open(self.records, encoding="utf-8").read()
        code, _, err = self._run(
            [self.records, "--table", self.table, "--report", self.report,
             "--fulltext-dir", self.tmp]
        )
        self.assertEqual(code, 1)
        self.assertIn("refused", err)
        self.assertEqual(open(self.records, encoding="utf-8").read(), unchanged)

    def test_an_empty_table_is_refused(self):
        json.dump({"_meta": {}, "assignments": []}, open(self.table, "w"))
        code, _, err = self._run(
            [self.records, "--table", self.table, "--report", self.report]
        )
        self.assertEqual(code, 2)
        self.assertIn("no assignments", err)

    def test_the_parser_declares_its_defaults(self):
        args = build_parser().parse_args([])
        self.assertIsNone(args.records)
        self.assertEqual(args.decision, "D-062")
        self.assertFalse(args.allow_refusals)


class CommittedCorpusAdmissibilityTest(unittest.TestCase):
    """D-075: this pass let the arm rule read a Methods paragraph.

    The arm rule decides whether an allele is the recorded perturbation or the
    substrate it was applied to **for this record's arm**. A Methods paragraph
    describes every arm, so it cannot pin a record to one -- and six committed
    records carry a genotype whose evidence is sourced ``methods`` with the
    rule ``llm_read_the_passage_under_the_arm_rule``. The rule name says a
    passage was read; the source says which section that passage was in.

    C009 could not see them until ``background_genotype`` entered
    ``conditions.DIMENSIONS`` (D-075). The number is recorded in the decision
    entry rather than here, because a test that must be edited when the data
    is corrected is a test that argues against correcting it. What is asserted
    here is the invariant that holds before and after: no Methods-sourced
    genotype entry is admissible, because ``constant`` is unreachable for this
    dimension.
    """

    def _corpus(self):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        path = os.path.join(root, "data", "curated", "records.jsonl")
        if not os.path.exists(path):  # pragma: no cover
            self.skipTest("no committed corpus")
        with open(path, encoding="utf-8") as handle:
            return [json.loads(l) for l in handle if l.strip()]

    def test_no_methods_sourced_genotype_entry_is_admissible(self):
        from syn4rch.curate.conditions import inadmissible_methods_conditions

        class _R:
            def __init__(self, evidence):
                self.condition_evidence = evidence

        class _P:
            verdicts = {"background_genotype": "varies"}

        for row in self._corpus():
            methods_entries = [
                e for e in row.get("condition_evidence") or []
                if e.get("dimension") == "background_genotype"
                and e.get("source") == "methods"
            ]
            if not methods_entries:
                continue
            self.assertTrue(
                inadmissible_methods_conditions(_R(methods_entries), _P()),
                row.get("record_id"),
            )

    def test_the_passage_route_carries_the_overwhelming_majority(self):
        """If this inverts, the fill has started preferring the Methods."""
        sources = {}
        for row in self._corpus():
            for e in row.get("condition_evidence") or []:
                if e.get("dimension") == "background_genotype":
                    sources[e.get("source")] = sources.get(e.get("source"), 0) + 1
        self.assertGreater(
            sources.get("passage", 0), 5 * sources.get("methods", 0),
            "the arm rule belongs to a record's own passage (D-060): %s" % sources,
        )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
