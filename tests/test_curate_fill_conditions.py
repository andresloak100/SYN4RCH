"""The condition fill decides, and refuses, exactly what it says it does.

These tests drive :func:`plan_record` and :func:`fill_records` with
hand-built :class:`PaperConditions`, so that the decision logic is tested
apart from the text extraction that feeds it. The extractor has its own
tests in ``test_curate_conditions.py``; what is at stake here is precedence,
what counts as a correction, and what the pass refuses to touch.
"""

import _path  # noqa: F401
import json
import os
import tempfile
import unittest

import synthetic
from syn4rch.curate.conditions import Extracted, PaperConditions
from syn4rch.curate.fill_conditions import (
    REFUSED_NOT_AT_LOCATOR,
    CONFIRMED,
    CONFLICTED,
    CORRECTED,
    FILLED,
    TO_VERSION,
    fill_records,
    fill_stats,
    plan_record,
    render_report_md,
    verify_plan,
)
from syn4rch.schema.record import TierARecord


def extracted(dimension, value, source="methods", locator="Methods para 1", rule="r"):
    return Extracted(
        dimension=dimension,
        value=value,
        source=source,
        locator=locator,
        quote="Synthetic quote: %s was %s." % (dimension, value),
        rule=rule,
    )


def paper(**constants):
    verdicts = {d: "constant" for d in constants}
    return PaperConditions(
        study_id="synthetic:study-1",
        pmcid="PMC0000001",
        verdicts=verdicts,
        constants={d: extracted(d, v) for d, v in constants.items()},
    )


def record(**overrides):
    """A record with no condition evidence and nothing filled by default."""
    data = dict(
        sex=None,
        diet=None,
        age_group=None,
        challenge=None,
        challenge_detail=None,
        condition_evidence=[],
        age_value=None,
        age_unit=None,
        species=None,
    )
    data.update(overrides)
    return synthetic.good_record(**data)


class TestTheFourOutcomes(unittest.TestCase):
    def test_a_null_field_is_filled_and_the_id_changes(self):
        r = record(sex=None)
        plan = plan_record(r, paper(sex="male"))
        self.assertEqual([o.outcome for o in plan.outcomes], [FILLED])
        self.assertEqual(plan.updates["sex"], "male")
        self.assertTrue(plan.changes_identity)

    def test_not_reported_is_corrected_because_it_was_asserted_too_narrowly(self):
        """``not_reported`` claims the paper is silent. The Methods are not."""
        r = record(sex="not_reported")
        plan = plan_record(r, paper(sex="female"))
        self.assertEqual([o.outcome for o in plan.outcomes], [CORRECTED])
        self.assertEqual(plan.updates["sex"], "female")
        self.assertTrue(plan.changes_identity)

    def test_an_agreeing_value_is_confirmed_and_keeps_its_id(self):
        r = record(sex="male")
        plan = plan_record(r, paper(sex="male"))
        self.assertEqual([o.outcome for o in plan.outcomes], [CONFIRMED])
        self.assertEqual(plan.updates, {})
        self.assertFalse(plan.changes_identity)
        self.assertEqual(len(plan.evidence), 1)
        self.assertEqual(plan.evidence[0].dimension, "sex")

    def test_a_disagreeing_value_is_reported_and_left_alone(self):
        r = record(sex="male")
        plan = plan_record(r, paper(sex="female"))
        self.assertEqual([o.outcome for o in plan.outcomes], [CONFLICTED])
        self.assertEqual(plan.updates, {})
        self.assertEqual(plan.evidence, [])
        self.assertTrue(plan.notes)
        self.assertIn("left as curated", plan.notes[0])


class TestPrecedence(unittest.TestCase):
    def test_a_paper_level_value_never_overwrites_a_curated_one(self):
        """The Methods describe the study; the passage describes this panel."""
        r = record(diet="chow")
        plan = plan_record(r, paper(diet="high_fat"))
        self.assertEqual(plan.updates, {})
        self.assertEqual([o.outcome for o in plan.outcomes], [CONFLICTED])

    def test_challenge_is_never_taken_from_the_methods(self):
        """An acute state is a property of a measurement, not of a study."""
        p = paper(sex="male")
        p.verdicts["challenge"] = "constant"
        p.constants["challenge"] = extracted("challenge", "fasted")
        plan = plan_record(record(), p)
        dims = [o.dimension for o in plan.outcomes]
        self.assertIn("sex", dims)
        self.assertNotIn("challenge", dims)
        self.assertNotIn("challenge", plan.updates)

    def test_a_dimension_the_paper_does_not_settle_is_left_untouched(self):
        p = PaperConditions(
            study_id="synthetic:study-1",
            pmcid="PMC0000001",
            verdicts={"sex": "varies", "diet": "absent"},
            constants={},
        )
        plan = plan_record(record(), p)
        self.assertEqual(plan.outcomes, [])
        self.assertEqual(plan.updates, {})

    def test_a_record_with_no_paper_at_all_is_left_untouched(self):
        plan = plan_record(record(), None)
        self.assertEqual(plan.updates, {})
        self.assertEqual(plan.evidence, [])


class TestEvidenceIsNotDuplicated(unittest.TestCase):
    def test_a_dimension_that_already_has_evidence_gains_none(self):
        r = record(
            sex="male",
            condition_evidence=[
                {
                    "dimension": "sex",
                    "value": "male",
                    "source": "passage",
                    "locator": "Results para 2",
                    "quote": "Synthetic quote: male mice.",
                    "rule": "earlier_pass",
                }
            ],
        )
        plan = plan_record(r, paper(sex="male"))
        self.assertEqual([o.outcome for o in plan.outcomes], [CONFIRMED])
        self.assertEqual(plan.evidence, [])


class TestAgeRidesAlongWithItsBand(unittest.TestCase):
    def test_a_scalar_age_and_species_are_filled_beside_the_band(self):
        p = paper(age_group="young_adult")
        p.age_value = 12.0
        p.age_unit = "weeks"
        p.species = "mus_musculus"
        plan = plan_record(record(), p)
        self.assertEqual(plan.updates["age_group"], "young_adult")
        self.assertEqual(plan.updates["age_value"], 12.0)
        self.assertEqual(plan.updates["age_unit"], "weeks")
        self.assertEqual(plan.updates["species"], "mus_musculus")

    def test_a_curated_age_value_is_not_overwritten(self):
        p = paper(age_group="young_adult")
        p.age_value = 12.0
        p.age_unit = "weeks"
        plan = plan_record(record(age_value=9.0, age_unit="weeks"), p)
        self.assertNotIn("age_value", plan.updates)


class TestFillRecordsOverADataset(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def test_every_record_is_rewritten_to_the_new_version(self):
        rows = [record(outcome_raw="synthetic readout %d" % i) for i in range(3)]
        result = fill_records(rows, self.dir)
        self.assertEqual(len(result.records), 3)
        for r in result.records:
            self.assertEqual(r.schema_version, TO_VERSION)

    def test_the_id_map_is_complete_over_the_input(self):
        rows = [record(outcome_raw="synthetic readout %d" % i) for i in range(4)]
        result = fill_records(rows, self.dir)
        self.assertEqual(
            sorted(result.id_map), sorted(r.record_id for r in rows)
        )

    def test_a_paper_with_no_retrieved_text_is_named_not_guessed(self):
        rows = [record()]
        result = fill_records(rows, self.dir)
        self.assertEqual(result.missing_fulltext, [rows[0].study_id])
        self.assertEqual(result.counts["papers"], 0)


class TestIdentityCollisions(unittest.TestCase):
    """Filling a gap can merge two records that differed only by that gap."""

    def test_two_records_differing_only_in_an_unfilled_sex_collide(self):
        filled_in = record(sex="male")
        gap = record(sex=None)
        # Same study, population, perturbation, outcome and locator: the only
        # hash component that differs is the one about to be filled.
        self.assertNotEqual(filled_in.record_id, gap.record_id)

        rows = [filled_in, gap]
        result = fill_records(rows, tempfile.mkdtemp())
        # With no full text there is nothing to fill, so no collision yet.
        self.assertEqual(result.collisions, [])

        # Now force the fill by planning directly and applying it.
        plan = plan_record(gap, paper(sex="male"))
        data = gap.to_dict()
        data.update(plan.updates)
        data["record_id"] = None
        merged = TierARecord.from_dict(data)
        self.assertEqual(merged.record_id, filled_in.record_id)


class TestReportAndStats(unittest.TestCase):
    def test_the_report_states_the_four_outcomes_and_names_the_rules(self):
        rows = [record(sex="not_reported"), record(sex="male", outcome_raw="o2")]
        result = fill_records(rows, tempfile.mkdtemp())
        # No full text, so nothing fired; the report must still be honest.
        stats = fill_stats(rows, result.records, result)
        text = render_report_md(stats, generated="2026-01-01T00:00:00Z")
        self.assertIn("No rule fired. Nothing was filled.", text)
        self.assertIn("filled", text)
        self.assertIn("conflicted", text)
        self.assertIn("left exactly as curated", text)

    def test_stats_are_recomputed_from_the_records_not_the_plans(self):
        rows = [record(sex=None)]
        result = fill_records(rows, tempfile.mkdtemp())
        stats = fill_stats(rows, result.records, result)
        self.assertEqual(stats["dimensions"]["sex"]["before"], {"null": 1})
        self.assertEqual(stats["dimensions"]["sex"]["after"], {"null": 1})

    def test_a_pass_that_introduces_a_hard_finding_is_visible_in_the_stats(self):
        rows = [record()]
        result = fill_records(rows, tempfile.mkdtemp())
        stats = fill_stats(rows, result.records, result)
        self.assertEqual(stats["hard_findings_before"], stats["hard_findings_after"])


if __name__ == "__main__":
    unittest.main()


class TestThePassageRouteOutranksTheMethods(unittest.TestCase):
    """A passage describes this measurement; the Methods describe the study."""

    def test_a_passage_value_wins_the_dimension_outright(self):
        r = record(sex=None)
        passage = extracted(
            "sex", "female", source="passage", locator="Results para 4",
            rule="quote_names_the_sex",
        )
        plan = plan_record(r, paper(sex="male"), record_extractions=[passage])
        self.assertEqual(plan.updates["sex"], "female")
        self.assertEqual(len(plan.evidence), 1)
        self.assertEqual(plan.evidence[0].source.value, "passage")
        self.assertEqual([o.rule for o in plan.outcomes], ["quote_names_the_sex"])

    def test_the_methods_still_fill_a_dimension_the_passage_is_silent_on(self):
        r = record(sex=None, diet=None)
        passage = extracted("sex", "female", source="passage", locator="Results para 4")
        plan = plan_record(
            r, paper(sex="male", diet="chow"), record_extractions=[passage]
        )
        self.assertEqual(plan.updates["sex"], "female")
        self.assertEqual(plan.updates["diet"], "chow")

    def test_a_passage_value_contradicting_a_curated_one_is_still_not_applied(self):
        """Precedence decides which reading is used, not whether to overwrite."""
        r = record(sex="male")
        passage = extracted("sex", "female", source="passage", locator="Results para 4")
        plan = plan_record(r, None, record_extractions=[passage])
        self.assertEqual(plan.updates, {})
        self.assertEqual([o.outcome for o in plan.outcomes], [CONFLICTED])


class _FakeIndex:
    """Resolves one locator to one passage. Enough to drive verify_plan."""

    class _P:
        def __init__(self, text):
            self.text = text

    def __init__(self, mapping):
        self._mapping = mapping

    def resolve(self, locator):
        if locator in self._mapping:
            return [self._P(self._mapping[locator])], "exact"
        return [], "none"


class TestEvidenceMustSitWhereItSaysItDoes(unittest.TestCase):
    """A real quote at the wrong address misdirects an auditor."""

    def test_a_quote_present_in_its_passage_is_kept(self):
        r = record(sex=None)
        plan = plan_record(
            r, None,
            record_extractions=[
                extracted("sex", "male", source="passage", locator="Results para 2")
            ],
        )
        index = _FakeIndex({"Results para 2": "Synthetic quote: sex was male. More text."})
        self.assertEqual(verify_plan(plan, index), 0)
        self.assertEqual(len(plan.evidence), 1)
        self.assertEqual(plan.updates["sex"], "male")

    def test_a_quote_absent_from_its_passage_loses_the_value_too(self):
        r = record(sex=None)
        plan = plan_record(
            r, None,
            record_extractions=[
                extracted("sex", "male", source="passage", locator="Results para 2")
            ],
        )
        index = _FakeIndex({"Results para 2": "A paragraph that says nothing of the sort."})
        self.assertEqual(verify_plan(plan, index), 1)
        self.assertEqual(plan.evidence, [])
        self.assertNotIn("sex", plan.updates)
        self.assertEqual(plan.outcomes, [])
        self.assertTrue(
            any(n.startswith(REFUSED_NOT_AT_LOCATOR) for n in plan.notes),
            plan.notes,
        )

    def test_an_unresolvable_locator_is_a_refusal_not_a_pass(self):
        r = record(sex=None)
        plan = plan_record(
            r, None,
            record_extractions=[
                extracted("sex", "male", source="passage", locator="Nowhere para 9")
            ],
        )
        self.assertEqual(verify_plan(plan, _FakeIndex({})), 1)
        self.assertEqual(plan.evidence, [])

    def test_no_paper_means_no_verification_and_so_no_evidence(self):
        """An unverifiable claim is not a verified one."""
        r = record(sex=None)
        plan = plan_record(
            r, None,
            record_extractions=[
                extracted("sex", "male", source="passage", locator="Results para 2")
            ],
        )
        self.assertEqual(verify_plan(plan, None), 1)
        self.assertEqual(plan.updates, {})

    def test_confirmed_values_survive_losing_their_evidence(self):
        """Dropping evidence for a value already curated leaves the value."""
        r = record(sex="male")
        plan = plan_record(
            r, None,
            record_extractions=[
                extracted("sex", "male", source="passage", locator="Results para 2")
            ],
        )
        verify_plan(plan, _FakeIndex({"Results para 2": "nothing relevant here"}))
        self.assertEqual(plan.evidence, [])
        self.assertEqual(plan.updates, {})
        self.assertEqual(r.sex.value, "male")


class TestTheTopUpIsIdempotent(unittest.TestCase):
    """A later curation pass tops up new papers without disturbing old ones."""

    def test_a_second_pass_over_filled_records_changes_no_identity(self):
        evidence = [{
            "dimension": "sex", "value": "male", "source": "methods",
            "locator": "Methods para 1", "quote": "Synthetic quote: male mice.",
            "rule": "methods_single_sex",
        }]
        r = record(sex="male", condition_evidence=evidence)
        before = r.record_id
        plan = plan_record(r, paper(sex="male"))
        self.assertEqual([o.outcome for o in plan.outcomes], [CONFIRMED])
        self.assertEqual(plan.updates, {})
        self.assertEqual(plan.evidence, [], "a confirmed value gains no second entry")
        self.assertFalse(plan.changes_identity)
        self.assertEqual(r.record_id, before)

    def test_both_input_versions_are_declared_and_the_written_one_is_fixed(self):
        from syn4rch.curate.fill_conditions import ACCEPTED_INPUT_VERSIONS, FROM_VERSION
        self.assertIn(FROM_VERSION, ACCEPTED_INPUT_VERSIONS)
        self.assertIn(TO_VERSION, ACCEPTED_INPUT_VERSIONS)
        self.assertEqual(TO_VERSION, "0.3.0")
