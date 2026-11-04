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



class PerturbationClassTest(unittest.TestCase):
    """``perturbation_class`` is a total function of ``method`` (v0.2.0)."""

    def test_table_is_total_over_the_method_vocabulary(self):
        from syn4rch.schema.enums import PerturbationMethod
        from syn4rch.schema.rules import METHOD_CLASS

        self.assertEqual(set(METHOD_CLASS), set(PerturbationMethod))

    def test_every_method_yields_a_class(self):
        from syn4rch.schema.enums import PerturbationMethod
        from syn4rch.schema.rules import perturbation_class_from_method

        for method in PerturbationMethod:
            self.assertIsNotNone(perturbation_class_from_method(method))

    def test_the_split_that_matters(self):
        from syn4rch.schema.enums import PerturbationClass
        from syn4rch.schema.rules import perturbation_class_from_method

        self.assertIs(
            perturbation_class_from_method("chemogenetic"),
            PerturbationClass.circuit_activity,
        )
        self.assertIs(
            perturbation_class_from_method("conditional_knockout"),
            PerturbationClass.molecular,
        )

    def test_absent_or_unknown_method_gives_none_not_other(self):
        from syn4rch.schema.rules import perturbation_class_from_method

        self.assertIsNone(perturbation_class_from_method(None))
        # 'other' is a curator's explicit declaration, never a fallback for an
        # unrecognised string: writing it would hide a vocabulary gap.
        self.assertIsNone(perturbation_class_from_method("not_a_method"))

    def test_string_and_enum_inputs_agree(self):
        from syn4rch.schema.enums import PerturbationMethod
        from syn4rch.schema.rules import perturbation_class_from_method

        self.assertIs(
            perturbation_class_from_method("optogenetic"),
            perturbation_class_from_method(PerturbationMethod.optogenetic),
        )


class PerturbationPolarityTest(unittest.TestCase):
    """``perturbation_polarity`` is the sign along the class's own axis."""

    def polarity(self, method, valence):
        from syn4rch.schema.rules import perturbation_polarity_from

        value = perturbation_polarity_from(method, valence)
        return None if value is None else value.value

    def test_circuit_manipulation_restates_valence(self):
        self.assertEqual(self.polarity("chemogenetic", "activation"), "increase")
        self.assertEqual(self.polarity("chemogenetic", "inhibition"), "decrease")
        self.assertEqual(self.polarity("optogenetic", "ablation"), "remove")
        self.assertEqual(self.polarity("chemogenetic", "unclear"), "unclear")

    def test_a_circuit_manipulation_without_valence_stays_unclear(self):
        # The same virus can excite or silence; the method does not say which.
        self.assertEqual(self.polarity("chemogenetic", None), "unclear")

    def test_molecular_polarity_comes_from_the_method_and_ignores_valence(self):
        """The whole point of the v0.2.0 split.

        A knockout's ``valence`` is honestly ``unclear`` and must stay so, but
        the *manipulation* has a perfectly well defined sign.
        """
        for valence in (None, "unclear", "activation", "inhibition"):
            self.assertEqual(
                self.polarity("conditional_knockout", valence), "decrease"
            )
            self.assertEqual(
                self.polarity("conditional_overexpression", valence), "increase"
            )

    def test_structural_removal_is_remove_whatever_the_valence(self):
        for valence in (None, "unclear", "activation"):
            self.assertEqual(self.polarity("genetic_ablation", valence), "remove")
            self.assertEqual(self.polarity("lesion", valence), "remove")

    def test_pharmacology_reads_the_valence(self):
        self.assertEqual(
            self.polarity("local_pharmacology", "inhibition"), "decrease"
        )
        self.assertEqual(self.polarity("local_pharmacology", None), "unclear")

    def test_other_is_unclear_not_a_guess(self):
        self.assertEqual(self.polarity("other", "activation"), "unclear")

    def test_absent_method_gives_none(self):
        self.assertIsNone(self.polarity(None, "activation"))
        self.assertIsNone(self.polarity("not_a_method", "activation"))

    def test_unrecognised_valence_does_not_crash_or_guess(self):
        self.assertEqual(self.polarity("chemogenetic", "sideways"), "unclear")

    def test_every_method_and_valence_pair_is_defined(self):
        from syn4rch.schema.enums import PerturbationMethod, Valence

        for method in PerturbationMethod:
            for valence in list(Valence) + [None]:
                self.assertIsNotNone(self.polarity(method, valence))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
