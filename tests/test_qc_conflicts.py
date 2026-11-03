"""Conflict detection: it must find the conflict and leave it unresolved."""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.qc import QCSeverity, annotate_conflicts, directions_conflict, find_conflicts
from syn4rch.schema import Direction


def planted_pair():
    """Two synthetic studies that disagree about the same synthetic claim."""
    up = synthetic.good_record(
        direction="increase",
        support_quote="Synthetic quote: the readout increased by 25% after activation.",
    )
    down = synthetic.other_lab_record(
        direction="decrease",
        effect_size_value=-30.0,
        support_quote="Synthetic quote: the readout decreased by 30% after activation.",
    )
    return up, down


class TestDirectionsConflict(unittest.TestCase):
    def test_increase_vs_decrease_is_opposed(self):
        self.assertEqual(
            directions_conflict(Direction.increase, Direction.decrease),
            "QC_CONFLICT_OPPOSED",
        )
        self.assertEqual(
            directions_conflict(Direction.decrease, Direction.increase),
            "QC_CONFLICT_OPPOSED",
        )

    def test_effect_vs_null(self):
        for effect in (Direction.increase, Direction.decrease):
            self.assertEqual(
                directions_conflict(effect, Direction.no_change),
                "QC_CONFLICT_EFFECT_VS_NULL",
            )
            self.assertEqual(
                directions_conflict(Direction.no_change, effect),
                "QC_CONFLICT_EFFECT_VS_NULL",
            )

    def test_agreement_is_not_a_conflict(self):
        for direction in Direction:
            self.assertIsNone(directions_conflict(direction, direction))

    def test_not_reported_conflicts_with_nothing(self):
        for direction in Direction:
            self.assertIsNone(directions_conflict(Direction.not_reported, direction))
            self.assertIsNone(directions_conflict(direction, Direction.not_reported))

    def test_non_monotonic_is_not_called_a_contradiction(self):
        for direction in (Direction.increase, Direction.decrease, Direction.no_change):
            self.assertIsNone(directions_conflict(Direction.non_monotonic, direction))

    def test_null_direction_is_not_a_conflict(self):
        self.assertIsNone(directions_conflict(None, Direction.increase))


class TestFindConflicts(unittest.TestCase):
    def test_planted_conflict_is_found(self):
        up, down = planted_pair()
        findings = find_conflicts([up, down])
        self.assertEqual([f.code for f in findings], ["QC_CONFLICT_OPPOSED"])
        finding = findings[0]
        self.assertEqual(finding.severity, QCSeverity.INFO)
        self.assertEqual(sorted(finding.record_ids), sorted([up.record_id, down.record_id]))
        self.assertIn("Conflicts are data", finding.message)

    def test_both_sides_are_reported_with_their_evidence(self):
        up, down = planted_pair()
        sides = find_conflicts([up, down])[0].details["sides"]
        self.assertEqual(len(sides), 2)
        self.assertEqual(
            {s["direction"] for s in sides}, {"increase", "decrease"}
        )
        for side in sides:
            self.assertTrue(side["support_quote"])
            self.assertTrue(side["source_locator"])
            self.assertTrue(side["study_id"].startswith("synthetic:"))

    def test_effect_vs_null_conflict_is_found(self):
        up = synthetic.good_record()
        null = synthetic.other_lab_record(
            direction="no_change",
            effect_size_value=None,
            effect_size_type=None,
            effect_size_units=None,
            support_quote="Synthetic quote: no significant difference was observed.",
            confidence="moderate",
        )
        self.assertEqual(
            [f.code for f in find_conflicts([up, null])],
            ["QC_CONFLICT_EFFECT_VS_NULL"],
        )

    def test_agreeing_studies_produce_no_conflict(self):
        up = synthetic.good_record()
        also_up = synthetic.other_lab_record()
        self.assertEqual(find_conflicts([up, also_up]), [])

    def test_within_study_disagreement_is_not_a_cross_study_conflict(self):
        a = synthetic.good_record(direction="increase")
        b = synthetic.good_record(
            direction="decrease",
            source_locator="Fig S1B (synthetic)",
            effect_size_value=-25.0,
            support_quote="Synthetic quote: the readout decreased by 25%.",
        )
        self.assertEqual(find_conflicts([a, b]), [])

    def test_different_claims_do_not_conflict(self):
        up = synthetic.good_record()
        for override in (
            {"outcome_id": "syn:outcome_2", "outcome_raw": "synthetic readout two"},
            {"valence": "inhibition"},
            {"population_id": "syn:pop_beta", "population_raw": "synthetic population beta"},
        ):
            down = synthetic.other_lab_record(
                direction="decrease",
                effect_size_value=-30.0,
                support_quote="Synthetic quote: the readout decreased by 30%.",
                **override
            )
            self.assertEqual(find_conflicts([up, down]), [], override)

    def test_three_way_disagreement_reports_each_pair(self):
        up, down = planted_pair()
        third = synthetic.same_lab_other_study_record(
            direction="decrease",
            effect_size_value=-10.0,
            support_quote="Synthetic quote: the readout decreased by 10%.",
        )
        findings = find_conflicts([up, down, third])
        self.assertEqual(len(findings), 2)  # up/down and up/third
        self.assertTrue(all(f.code == "QC_CONFLICT_OPPOSED" for f in findings))

    def test_records_missing_an_outcome_id_are_skipped(self):
        up, down = planted_pair()
        up.outcome_id = None
        self.assertEqual(find_conflicts([up, down]), [])


class TestAnnotateConflictsNeverResolves(unittest.TestCase):
    def test_conflicts_with_is_populated_on_both_sides(self):
        up, down = planted_pair()
        annotate_conflicts([up, down])
        self.assertEqual(up.conflicts_with, [down.record_id])
        self.assertEqual(down.conflicts_with, [up.record_id])

    def test_nothing_is_dropped_averaged_or_rewritten(self):
        up, down = planted_pair()
        before = [
            (r.direction, r.effect_size_value, r.support_quote, r.study_id, r.record_id)
            for r in (up, down)
        ]
        records = [up, down]
        annotate_conflicts(records)
        self.assertEqual(len(records), 2)
        after = [
            (r.direction, r.effect_size_value, r.support_quote, r.study_id, r.record_id)
            for r in (up, down)
        ]
        self.assertEqual(before, after)
        self.assertIs(up.direction, Direction.increase)
        self.assertIs(down.direction, Direction.decrease)

    def test_annotation_is_idempotent(self):
        up, down = planted_pair()
        for _ in range(3):
            annotate_conflicts([up, down])
        self.assertEqual(up.conflicts_with, [down.record_id])
        self.assertEqual(down.conflicts_with, [up.record_id])

    def test_pre_existing_manual_annotations_are_preserved(self):
        up, down = planted_pair()
        up.conflicts_with = ["ffffffffffffffff"]
        annotate_conflicts([up, down])
        self.assertEqual(
            up.conflicts_with, sorted(["ffffffffffffffff", down.record_id])
        )

    def test_a_record_never_conflicts_with_itself(self):
        up, down = planted_pair()
        annotate_conflicts([up, down])
        self.assertNotIn(up.record_id, up.conflicts_with)

    def test_three_way_annotation_lists_every_counterpart(self):
        up, down = planted_pair()
        third = synthetic.same_lab_other_study_record(
            direction="decrease",
            effect_size_value=-10.0,
            support_quote="Synthetic quote: the readout decreased by 10%.",
        )
        annotate_conflicts([up, down, third])
        self.assertEqual(up.conflicts_with, sorted([down.record_id, third.record_id]))
        self.assertEqual(down.conflicts_with, [up.record_id])
        self.assertEqual(third.conflicts_with, [up.record_id])

    def test_conflicts_survive_a_jsonl_round_trip(self):
        import os
        import shutil
        import tempfile

        from syn4rch.schema import read_jsonl, write_jsonl

        up, down = planted_pair()
        annotate_conflicts([up, down])
        tmp = tempfile.mkdtemp(prefix="syn4rch-test-")
        try:
            path = os.path.join(tmp, "records.jsonl")
            write_jsonl(path, [up, down])
            loaded = read_jsonl(path)
            self.assertEqual(loaded[0].conflicts_with, [down.record_id])
            self.assertEqual(loaded[1].conflicts_with, [up.record_id])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
