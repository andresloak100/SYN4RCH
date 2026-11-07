"""The condition corrector (D-053, built in D-058).

`fill_conditions` fills an empty dimension and refuses to overwrite a curated
one. This module overwrites curated values, so the tests here are mostly about
what it refuses to do: the two locks (the declared `from` value and the stored
quote), the re-identification, and the fact that nothing else moves.
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
from syn4rch.curate.correct_conditions import (
    CORRECTIONS,
    CORRECTION_RULE,
    Correction,
    CorrectionCollision,
    OUTCOMES,
    apply_corrections,
    build_parser,
    main,
    plan_corrections,
    render_report_md,
)
from syn4rch.schema.enums import Challenge
from syn4rch.schema.io import read_jsonl, write_jsonl
from syn4rch.schema.record import CONDITION_DIMENSIONS, RECORD_ID_FIELDS, TierARecord
from syn4rch.schema.validate import Severity, validate

QUOTE = "Synthetic quote: mice were fasted overnight and intake rose by 25%."


def record(**overrides):
    data = dict(
        challenge="fasted",
        challenge_detail="fasted overnight (fabricated)",
        support_quote=QUOTE,
    )
    data.update(overrides)
    return synthetic.good_record(**data)


def with_challenge_quote(rec, quote=QUOTE):
    """Point the record's `challenge` evidence at ``quote``."""
    for entry in rec.condition_evidence or []:
        if entry.dimension == "challenge":
            entry.quote = quote
    return rec


def correction(**overrides):
    data = dict(
        record_id="",
        dimension="challenge",
        from_value="fasted",
        to_value="refed",
        quote_contains="fasted overnight",
        reason="synthetic reason",
        decision="D-053",
    )
    data.update(overrides)
    return Correction(**data)


def outcomes_of(records, corrections):
    return [o.outcome for o in plan_corrections(records, corrections)]


class DeclaredCorrectionsTest(unittest.TestCase):
    """The table itself, which is the thing a reviewer reads in a diff."""

    def test_every_correction_names_a_declared_condition_dimension(self):
        for item in CORRECTIONS:
            self.assertIn(item.dimension, CONDITION_DIMENSIONS, item.record_id)

    def test_every_value_is_legal_for_its_own_dimension_or_a_removal(self):
        """Per dimension, from the schema's own map.

        This read `{c.value for c in Challenge}` while the table held only
        `challenge` corrections, so the first `sex` correction (D-084) failed a
        test whose intent it satisfied. The enum comes from
        `record._ENUM_FIELDS` rather than a dict written here, so a dimension
        added to the schema cannot be checked against the wrong enum.
        """
        from syn4rch.schema.record import _ENUM_FIELDS

        for item in CORRECTIONS:
            enum = _ENUM_FIELDS[item.dimension]
            values = {member.value for member in enum}
            self.assertIn(item.from_value, values, "%s %s" % (item.record_id, item.dimension))
            if item.to_value is not None:
                self.assertIn(item.to_value, values, item.record_id)

    def test_the_table_covers_more_than_one_dimension(self):
        """Pins that the generalisation above is exercised, not merely written."""
        self.assertGreater(len({c.dimension for c in CORRECTIONS}), 1)

    def test_no_correction_is_a_no_op(self):
        for item in CORRECTIONS:
            self.assertNotEqual(item.from_value, item.to_value, item.record_id)

    def test_one_correction_per_record_and_dimension(self):
        keys = [(c.record_id, c.dimension) for c in CORRECTIONS]
        self.assertEqual(len(keys), len(set(keys)))

    def test_every_correction_carries_a_reason_a_decision_and_a_quote_lock(self):
        for item in CORRECTIONS:
            self.assertTrue(item.reason.strip(), item.record_id)
            self.assertTrue(item.decision.strip(), item.record_id)
            self.assertTrue(item.quote_contains.strip(), item.record_id)

    def test_the_corrected_dimension_is_a_record_id_component(self):
        """If it were not, no re-identification would be needed and this tool
        would not need an id map. The test pins the reason the map exists."""
        for item in CORRECTIONS:
            self.assertIn(item.dimension, RECORD_ID_FIELDS, item.dimension)


class RefusalTest(unittest.TestCase):
    def test_a_record_that_is_not_there_is_refused(self):
        rows = [with_challenge_quote(record())]
        self.assertEqual(
            outcomes_of(rows, [correction(record_id="not-an-id")]),
            ["refused_record_not_found"],
        )

    def test_a_dimension_that_is_not_a_condition_is_refused(self):
        rec = with_challenge_quote(record())
        self.assertEqual(
            outcomes_of(
                [rec], [correction(record_id=rec.record_id, dimension="direction")]
            ),
            ["refused_dimension_not_a_condition"],
        )

    def test_a_value_other_than_the_declared_from_is_refused(self):
        """The lock that matters most: never overwrite whatever happens to be there."""
        rec = with_challenge_quote(record(challenge="ad_libitum"))
        got = plan_corrections([rec], [correction(record_id=rec.record_id)])
        self.assertEqual(got[0].outcome, "refused_value_is_not_the_declared_from")
        self.assertEqual(got[0].stored_value, "ad_libitum")

    def test_a_value_with_no_evidence_entry_is_refused(self):
        rec = record()
        rec.condition_evidence = [
            e for e in rec.condition_evidence or [] if e.dimension != "challenge"
        ]
        self.assertEqual(
            outcomes_of([rec], [correction(record_id=rec.record_id)]),
            ["refused_no_evidence_entry"],
        )

    def test_a_quote_that_does_not_match_is_refused(self):
        """The second lock: a correction is addressed to the evidence too.

        A record_id hashes twelve fields and not the quote, so it survives
        re-curation of the sentence behind a value. Without this check a
        correction written about one sentence could apply to another.
        """
        rec = with_challenge_quote(
            record(), quote="Synthetic quote: something else rose by 25%."
        )
        self.assertEqual(
            outcomes_of([rec], [correction(record_id=rec.record_id)]),
            ["refused_quote_does_not_match"],
        )

    def test_two_corrections_for_one_dimension_are_refused(self):
        rec = with_challenge_quote(record())
        got = outcomes_of(
            [rec],
            [
                correction(record_id=rec.record_id),
                correction(record_id=rec.record_id, to_value="ad_libitum"),
            ],
        )
        self.assertEqual(got, ["applied", "refused_duplicate_correction"])

    def test_every_outcome_name_is_declared(self):
        rec = with_challenge_quote(record())
        got = plan_corrections([rec], [correction(record_id=rec.record_id)])
        for item in got:
            self.assertIn(item.outcome, OUTCOMES)


class ApplyTest(unittest.TestCase):
    def test_the_value_changes_and_the_id_is_recomputed(self):
        rec = with_challenge_quote(record())
        old = rec.record_id
        result = apply_corrections([rec], [correction(record_id=old)])
        after = result.records[0]
        self.assertEqual(after.challenge, Challenge.refed)
        self.assertNotEqual(after.record_id, old)
        self.assertEqual(after.record_id, after.compute_id())
        self.assertTrue(after.id_is_consistent())

    def test_the_evidence_keeps_its_sentence_and_changes_its_value_and_rule(self):
        rec = with_challenge_quote(record())
        before = next(
            e for e in rec.condition_evidence if e.dimension == "challenge"
        )
        result = apply_corrections([rec], [correction(record_id=rec.record_id)])
        after = next(
            e
            for e in result.records[0].condition_evidence
            if e.dimension == "challenge"
        )
        self.assertEqual(after.value, "refed")
        self.assertEqual(after.quote, before.quote)
        self.assertEqual(after.locator, before.locator)
        self.assertEqual(after.source, before.source)
        self.assertTrue(after.rule.startswith(CORRECTION_RULE))
        self.assertIn("d053", after.rule)

    def test_a_removal_drops_the_evidence_entry(self):
        """A quote left on a null dimension would assert a value that is not
        stored, which the validator refuses in both directions (A009, A011)."""
        rec = with_challenge_quote(record())
        result = apply_corrections(
            [rec], [correction(record_id=rec.record_id, to_value=None)]
        )
        after = result.records[0]
        self.assertIsNone(after.challenge)
        self.assertEqual(
            [e.dimension for e in after.condition_evidence if e.dimension == "challenge"],
            [],
        )
        self.assertEqual(
            [f for f in validate(after) if f.severity is Severity.HARD], []
        )

    def test_nothing_else_about_the_record_moves(self):
        rec = with_challenge_quote(record())
        before = rec.to_dict()
        after = apply_corrections(
            [rec], [correction(record_id=rec.record_id)]
        ).records[0].to_dict()
        moved = {
            key
            for key in set(before) | set(after)
            if before.get(key) != after.get(key)
        }
        self.assertEqual(moved, {"record_id", "challenge", "condition_evidence"})

    def test_the_input_records_are_not_mutated(self):
        rec = with_challenge_quote(record())
        before = rec.to_dict()
        apply_corrections([rec], [correction(record_id=rec.record_id)])
        self.assertEqual(rec.to_dict(), before)

    def test_applying_twice_is_a_no_op(self):
        rec = with_challenge_quote(record())
        first = apply_corrections([rec], [correction(record_id=rec.record_id)])
        again = apply_corrections(
            first.records, [correction(record_id=rec.record_id)]
        )
        # The record has been re-identified, so the correction no longer
        # addresses any record in the dataset.
        self.assertEqual([o.outcome for o in again.outcomes], ["refused_record_not_found"])
        self.assertEqual(again.ids_changed, 0)

    def test_a_correction_already_satisfied_reports_already_applied(self):
        rec = with_challenge_quote(record(challenge="refed", challenge_detail="refed"))
        got = plan_corrections([rec], [correction(record_id=rec.record_id)])
        self.assertEqual(got[0].outcome, "already_applied")

    def test_the_id_map_is_complete_and_bijective(self):
        rows = [
            with_challenge_quote(record()),
            with_challenge_quote(
                record(
                    study_id="synthetic:test-077",
                    doi="synthetic:test-077",
                    source_locator="Fig 7A (synthetic)",
                )
            ),
        ]
        result = apply_corrections(rows, [correction(record_id=rows[0].record_id)])
        self.assertEqual(set(result.id_map), {r.record_id for r in rows})
        self.assertEqual(len(set(result.id_map.values())), len(result.id_map))
        self.assertEqual(result.ids_changed, 1)
        self.assertEqual(result.id_map[rows[1].record_id], rows[1].record_id)

    def test_conflicts_with_is_rewritten_through_the_map(self):
        target = with_challenge_quote(record())
        other = with_challenge_quote(
            record(
                study_id="synthetic:test-078",
                doi="synthetic:test-078",
                source_locator="Fig 8A (synthetic)",
                challenge="ad_libitum",
                challenge_detail="ad libitum (fabricated)",
            )
        )
        other.conflicts_with = [target.record_id]
        result = apply_corrections(
            [target, other], [correction(record_id=target.record_id)]
        )
        new_target = result.records[0].record_id
        self.assertEqual(result.records[1].conflicts_with, [new_target])
        self.assertEqual(result.unmapped_conflicts, [])

    def test_a_dangling_conflict_reference_is_reported_not_dropped(self):
        rec = with_challenge_quote(record())
        rec.conflicts_with = ["0000000000000000"]
        result = apply_corrections([rec], [correction(record_id=rec.record_id)])
        self.assertEqual(result.records[0].conflicts_with, ["0000000000000000"])
        self.assertEqual(len(result.unmapped_conflicts), 1)

    def test_a_correction_that_would_collide_two_records_refuses(self):
        """Two records differing only in the corrected condition collapse into
        one id. Keeping one of them would delete an observation."""
        already = with_challenge_quote(
            record(challenge="refed", challenge_detail="fasted overnight (fabricated)")
        )
        target = with_challenge_quote(record())
        with self.assertRaises(CorrectionCollision):
            apply_corrections([target, already], [correction(record_id=target.record_id)])


class ReportTest(unittest.TestCase):
    def test_the_report_states_the_counts_and_the_reasons(self):
        rec = with_challenge_quote(record())
        before = [TierARecord.from_dict(rec.to_dict())]
        result = apply_corrections([rec], [correction(record_id=rec.record_id)])
        text = render_report_md(
            result, before, generated="2026-01-01T00:00:00Z", records_path="x.jsonl"
        )
        self.assertIn("Condition corrections (D-053)", text)
        self.assertIn("synthetic reason", text)
        self.assertIn("`refed`", text)
        self.assertIn("What this pass does not do", text)

    def test_locator_failures_are_printed_when_present(self):
        rec = with_challenge_quote(record())
        result = apply_corrections([rec], [correction(record_id=rec.record_id)])
        result.locator_failures = [
            {
                "record_id": "x",
                "dimension": "challenge",
                "locator": "Results para 1",
                "why": "quote does not occur in the passage it names",
            }
        ]
        text = render_report_md(
            result, [rec], generated="2026-01-01T00:00:00Z", records_path="x.jsonl"
        )
        self.assertIn("no longer verify at their locator", text)


class CLITest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.path = os.path.join(self.tmp, "records.jsonl")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = main(argv)
        return code, out.getvalue(), err.getvalue()

    def test_show_lists_every_declared_correction(self):
        code, out, _ = self._run(["show"])
        self.assertEqual(code, 0)
        for item in CORRECTIONS:
            self.assertIn(item.record_id, out)

    def test_show_json_is_machine_readable(self):
        code, out, _ = self._run(["show", "--json"])
        self.assertEqual(code, 0)
        self.assertEqual(len(json.loads(out)), len(CORRECTIONS))

    def test_an_empty_dataset_is_refused(self):
        open(self.path, "w", encoding="utf-8").close()
        code, _, err = self._run(["run", self.path])
        self.assertEqual(code, 2)
        self.assertIn("holds no records", err)

    def test_a_refused_correction_blocks_the_write(self):
        write_jsonl(self.path, [with_challenge_quote(record())])
        unchanged = open(self.path, encoding="utf-8").read()
        code, _, err = self._run(["run", self.path, "--fulltext-dir", self.tmp])
        self.assertEqual(code, 1)
        self.assertIn("refused", err)
        self.assertEqual(open(self.path, encoding="utf-8").read(), unchanged)

    def test_dry_run_writes_nothing(self):
        rec = with_challenge_quote(record())
        write_jsonl(self.path, [rec])
        unchanged = open(self.path, encoding="utf-8").read()
        code, out, _ = self._run(
            [
                "run", self.path,
                "--fulltext-dir", self.tmp,
                "--dry-run",
                "--allow-refusals",
            ]
        )
        self.assertEqual(code, 0)
        self.assertIn("dry run", out)
        self.assertEqual(open(self.path, encoding="utf-8").read(), unchanged)

    def test_the_parser_requires_a_subcommand(self):
        with self.assertRaises(SystemExit):
            build_parser().parse_args([])


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
