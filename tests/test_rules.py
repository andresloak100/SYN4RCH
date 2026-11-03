"""The two derived fields: age_group (SCHEMA.md 3.3) and confidence (3.6)."""

import _path  # noqa: F401
import itertools
import unittest

from syn4rch.schema import AgeGroup, Confidence, Species, age_in_weeks
from syn4rch.schema import assign_age_group as age_group
from syn4rch.schema import assign_confidence as confidence

MOUSE = "mus_musculus"


class TestAgeGroupBoundaries(unittest.TestCase):
    """The exact boundaries of SCHEMA.md 3.3, nothing more."""

    def test_juvenile_to_young_adult_boundary(self):
        self.assertIs(age_group(7, "weeks", MOUSE), AgeGroup.juvenile)
        self.assertIs(age_group(7.999, "weeks", MOUSE), AgeGroup.juvenile)
        self.assertIs(age_group(8, "weeks", MOUSE), AgeGroup.young_adult)

    def test_young_adult_to_middle_aged_boundary(self):
        self.assertIs(age_group(26, "weeks", MOUSE), AgeGroup.young_adult)
        self.assertIs(age_group(26.5, "weeks", MOUSE), AgeGroup.young_adult)
        self.assertIs(age_group(27, "weeks", MOUSE), AgeGroup.middle_aged)

    def test_middle_aged_to_aged_boundary(self):
        self.assertIs(age_group(64, "weeks", MOUSE), AgeGroup.middle_aged)
        self.assertIs(age_group(64.5, "weeks", MOUSE), AgeGroup.middle_aged)
        self.assertIs(age_group(65, "weeks", MOUSE), AgeGroup.aged)
        self.assertIs(age_group(120, "weeks", MOUSE), AgeGroup.aged)

    def test_zero_is_juvenile(self):
        self.assertIs(age_group(0, "weeks", MOUSE), AgeGroup.juvenile)


class TestAgeUnitConversion(unittest.TestCase):
    def test_days_weeks_months_years(self):
        self.assertAlmostEqual(age_in_weeks(56, "days"), 8.0)
        self.assertAlmostEqual(age_in_weeks(8, "weeks"), 8.0)
        self.assertAlmostEqual(age_in_weeks(12, "months"), 52.178571, places=4)
        self.assertAlmostEqual(age_in_weeks(1, "years"), 52.178571, places=4)

    def test_conversions_land_in_the_right_band(self):
        self.assertIs(age_group(56, "days", MOUSE), AgeGroup.young_adult)   # 8 wk
        self.assertIs(age_group(55, "days", MOUSE), AgeGroup.juvenile)      # 7.86 wk
        self.assertIs(age_group(3, "months", MOUSE), AgeGroup.young_adult)  # 13.0 wk
        self.assertIs(age_group(18, "months", MOUSE), AgeGroup.aged)        # 78.3 wk
        self.assertIs(age_group(2, "years", MOUSE), AgeGroup.aged)

    def test_unreported_and_nonsense_inputs(self):
        self.assertIsNone(age_in_weeks(None, "weeks"))
        self.assertIsNone(age_in_weeks(8, None))
        self.assertIsNone(age_in_weeks(8, "fortnights"))
        self.assertIsNone(age_in_weeks(-1, "weeks"))
        self.assertIsNone(age_in_weeks("old", "weeks"))
        for args in ((None, "weeks"), (8, None), (8, "fortnights"), (-1, "weeks")):
            self.assertIs(age_group(args[0], args[1], MOUSE), AgeGroup.not_reported)


class TestAgeGroupSpeciesRestriction(unittest.TestCase):
    """Mouse boundaries must not be extrapolated to other species."""

    def test_non_mouse_species_return_not_reported(self):
        for species in (
            Species.rattus_norvegicus,
            Species.homo_sapiens,
            Species.other,
            "rattus_norvegicus",
            "homo_sapiens",
        ):
            self.assertIs(
                age_group(65, "weeks", species),
                AgeGroup.not_reported,
                "species %s must not reuse mouse boundaries" % species,
            )

    def test_missing_or_unknown_species(self):
        self.assertIs(age_group(12, "weeks", None), AgeGroup.not_reported)
        self.assertIs(age_group(12, "weeks", "canis_familiaris"), AgeGroup.not_reported)

    def test_docstring_states_the_restriction(self):
        doc = age_group.__doc__ or ""
        self.assertIn("mouse-specific", doc)


class TestConfidenceTruthTable(unittest.TestCase):
    """All 16 input combinations are pinned here.

    Rule: low if no full text, or if fewer than 2 of the other three hold;
    high if all four hold; moderate otherwise.
    """

    EXPECTED = {
        # (full_text, actuator_control, n_reported, effect_size)
        (False, False, False, False): Confidence.low,
        (False, False, False, True): Confidence.low,
        (False, False, True, False): Confidence.low,
        (False, False, True, True): Confidence.low,
        (False, True, False, False): Confidence.low,
        (False, True, False, True): Confidence.low,
        (False, True, True, False): Confidence.low,
        (False, True, True, True): Confidence.low,
        (True, False, False, False): Confidence.low,
        (True, False, False, True): Confidence.low,
        (True, False, True, False): Confidence.low,
        (True, False, True, True): Confidence.moderate,
        (True, True, False, False): Confidence.low,
        (True, True, False, True): Confidence.moderate,
        (True, True, True, False): Confidence.moderate,
        (True, True, True, True): Confidence.high,
    }

    def test_table_is_complete(self):
        self.assertEqual(
            set(self.EXPECTED), set(itertools.product((False, True), repeat=4))
        )

    def test_every_row(self):
        for inputs, expected in sorted(self.EXPECTED.items()):
            self.assertIs(
                confidence(*inputs),
                expected,
                "confidence%s should be %s" % (inputs, expected.value),
            )

    def test_none_counts_as_not_met(self):
        self.assertIs(confidence(None, True, True, True), Confidence.low)
        self.assertIs(confidence(True, None, None, None), Confidence.low)
        self.assertIs(confidence(True, None, True, True), Confidence.moderate)

    def test_truthy_non_true_values_do_not_count(self):
        """Only strict ``True`` counts: 1, "yes" and [] must not pass."""
        self.assertIs(confidence(1, 1, 1, 1), Confidence.low)
        self.assertIs(confidence("yes", "yes", "yes", "yes"), Confidence.low)

    def test_high_requires_all_four(self):
        highs = [k for k, v in self.EXPECTED.items() if v is Confidence.high]
        self.assertEqual(highs, [(True, True, True, True)])


class TestActuatorControlHelper(unittest.TestCase):
    def test_recognises_actuator_isolating_controls(self):
        from syn4rch.schema import has_actuator_control_from_controls as helper

        self.assertTrue(helper(["opsin_free"]))
        self.assertTrue(helper(["vehicle", "cre_negative"]))
        self.assertFalse(helper(["vehicle"]))
        self.assertFalse(helper(["none_reported"]))
        self.assertFalse(helper([]))
        self.assertFalse(helper(None))
        self.assertFalse(helper(["not_a_control_type"]))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
