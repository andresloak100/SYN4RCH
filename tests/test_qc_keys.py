"""The grouping keys, and what schema v0.2.0 changed about "the same claim".

These tests exist because the claim key decides what counts as a contradiction
in this corpus. Getting it wrong does not look like a bug: it looks like the
literature disagreeing with itself.
"""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.qc import keys as K
from syn4rch.schema.enums import PerturbationClass, PerturbationPolarity


class PopulationKeyTest(unittest.TestCase):
    def test_resolved_population_uses_the_id(self):
        r = synthetic.good_record(population_id="syn:pop_alpha")
        self.assertEqual(K.population_key(r), "syn:pop_alpha")

    def test_unresolved_population_is_prefixed_and_never_pools(self):
        unresolved = synthetic.good_record(population_id=None)
        self.assertTrue(K.population_key(unresolved).startswith("raw:"))
        resolved = synthetic.good_record(
            population_id="synthetic population alpha"
        )
        self.assertNotEqual(K.population_key(unresolved), K.population_key(resolved))


class PerturbationKeyTest(unittest.TestCase):
    def test_class_polarity_and_target_travel_together(self):
        r = synthetic.molecular_record()
        self.assertEqual(
            K.perturbation_key(r), ("molecular", "decrease", "syn:synthgene1")
        )

    def test_knockout_and_overexpression_do_not_share_a_key(self):
        """The 170-record defect: both were valence 'unclear' under v0.1.0."""
        ko = synthetic.molecular_record(
            method="conditional_knockout", perturbation_polarity="decrease"
        )
        oe = synthetic.molecular_record(
            method="conditional_overexpression", perturbation_polarity="increase"
        )
        self.assertEqual(ko.valence, oe.valence)  # both honestly 'unclear'
        self.assertNotEqual(K.perturbation_key(ko), K.perturbation_key(oe))
        # ... and under the old key they were indistinguishable:
        self.assertEqual(K.claim_key_pooled(ko), K.claim_key_pooled(oe))

    def test_molecular_decrease_never_pools_with_circuit_decrease(self):
        """'decrease' means a gene product in one case and firing in the other."""
        molecular = synthetic.molecular_record()
        circuit = synthetic.good_record(
            valence="inhibition", perturbation_polarity="decrease"
        )
        self.assertEqual(
            K.perturbation_key(molecular)[1], K.perturbation_key(circuit)[1]
        )
        self.assertNotEqual(
            K.perturbation_key(molecular)[0], K.perturbation_key(circuit)[0]
        )
        self.assertNotEqual(K.claim_key(molecular), K.claim_key(circuit))

    def test_two_genes_do_not_pool(self):
        a = synthetic.molecular_record(target_id="syn:a", target_raw="SynA")
        b = synthetic.molecular_record(target_id="syn:b", target_raw="SynB")
        self.assertNotEqual(K.perturbation_key(a), K.perturbation_key(b))

    def test_unresolved_target_is_prefixed_and_never_pools(self):
        raw_only = synthetic.molecular_record(target_id=None, target_raw="SynA")
        self.assertTrue(K.perturbation_key(raw_only)[2].startswith("raw:"))
        resolved = synthetic.molecular_record(target_id="syna", target_raw="SynA")
        self.assertNotEqual(
            K.perturbation_key(raw_only), K.perturbation_key(resolved)
        )

    def test_absent_target_reads_null_not_empty(self):
        r = synthetic.good_record()
        self.assertEqual(K.perturbation_key(r)[2], K.NULL)

    def test_stored_fields_and_derivation_agree(self):
        """A pre-v0.2.0 record keys identically to the same record migrated."""
        stored = synthetic.good_record()
        underived = synthetic.good_record(
            perturbation_class=None, perturbation_polarity=None
        )
        self.assertEqual(K.perturbation_key(stored), K.perturbation_key(underived))

    def test_stored_field_wins_over_derivation(self):
        r = synthetic.good_record(
            perturbation_class=PerturbationClass.other,
            perturbation_polarity=PerturbationPolarity.bidirectional,
        )
        self.assertEqual(K.perturbation_key(r)[:2], ("other", "bidirectional"))


class ConditionKeyTest(unittest.TestCase):
    def test_known_conditions_appear_in_order(self):
        r = synthetic.good_record(
            sex="female", diet="high_fat", age_group="aged", challenge="fasted"
        )
        self.assertEqual(
            K.condition_key(r), ("female", "high_fat", "aged", "fasted")
        )
        self.assertEqual(K.CONDITION_FIELDS,
                         ("sex", "diet", "age_group", "challenge"))

    def test_none_and_not_reported_both_read_unknown(self):
        absent = synthetic.good_record(sex=None)
        checked = synthetic.good_record(sex="not_reported")
        self.assertEqual(K.condition_key(absent)[0], K.UNKNOWN)
        self.assertEqual(K.condition_key(checked)[0], K.UNKNOWN)
        self.assertEqual(K.condition_key(absent), K.condition_key(checked))

    def test_the_record_still_distinguishes_them(self):
        """Pooling is a grouping choice, not a loss of the stored distinction."""
        absent = synthetic.good_record(sex=None)
        checked = synthetic.good_record(sex="not_reported")
        self.assertIsNone(absent.sex)
        self.assertEqual(checked.sex.value, "not_reported")

    def test_differing_conditions_are_not_comparable(self):
        male = synthetic.good_record(sex="male")
        female = synthetic.good_record(sex="female")
        self.assertNotEqual(K.claim_key(male), K.claim_key(female))
        # Under the v0.1.0 key they were one claim, and a disagreement between
        # them was recorded as laboratories contradicting each other.
        self.assertEqual(K.claim_key_pooled(male), K.claim_key_pooled(female))

    def test_specificity_counts_known_slots(self):
        self.assertEqual(K.condition_specificity(synthetic.good_record()), 4)
        blind = synthetic.good_record(
            sex=None, diet=None, age_group=None, challenge=None
        )
        self.assertEqual(K.condition_specificity(blind), 0)
        partial = synthetic.good_record(
            sex="male", diet=None, age_group="not_reported", challenge=None
        )
        self.assertEqual(K.condition_specificity(partial), 1)


class ClaimKeyTest(unittest.TestCase):
    def test_component_names_match_the_key_length(self):
        key = K.claim_key(synthetic.good_record())
        self.assertEqual(len(key), len(K.CLAIM_KEY_FIELDS))

    def test_valence_is_absent_from_the_strict_key(self):
        self.assertNotIn("valence", K.CLAIM_KEY_FIELDS)
        self.assertIn("valence", K.POOLED_CLAIM_KEY_FIELDS)

    def test_pooled_key_is_exactly_the_v0_1_0_key(self):
        r = synthetic.good_record()
        self.assertEqual(
            K.claim_key_pooled(r),
            (K.norm(r.population_id), K.norm(r.valence), K.norm(r.outcome_id)),
        )

    def test_same_experiment_same_key(self):
        a = synthetic.good_record()
        b = synthetic.other_lab_record()
        self.assertEqual(K.claim_key(a), K.claim_key(b))

    def test_observation_key_is_study_scoped(self):
        a = synthetic.good_record()
        b = synthetic.other_lab_record()
        self.assertNotEqual(K.observation_key(a), K.observation_key(b))
        self.assertEqual(K.observation_key(a)[1:], K.claim_key(a))


class LabKeyTest(unittest.TestCase):
    def test_missing_last_author_cannot_be_grouped(self):
        self.assertIsNone(K.lab_key(synthetic.good_record(last_author=None)))

    def test_same_lab_different_study(self):
        a = synthetic.good_record()
        b = synthetic.same_lab_other_study_record()
        self.assertEqual(K.lab_key(a), K.lab_key(b))
        self.assertNotEqual(K.study_key(a), K.study_key(b))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
