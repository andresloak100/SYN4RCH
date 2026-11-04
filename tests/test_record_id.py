"""record_id must be deterministic, portable, and sensitive to its inputs."""

import _path  # noqa: F401
import inspect
import unittest

import synthetic
from syn4rch.schema import RECORD_ID_FIELDS, TierARecord, compute_record_id
from syn4rch.schema.record import RECORD_ID_LENGTH, normalize_for_hash

#: Golden digest. Pinned so that a change to the hash input set, the
#: separator, the normalisation or the digest length fails loudly instead of
#: silently re-identifying the whole corpus. Recomputing this value is a
#: schema change and needs a version bump.
#:
#: It was last recomputed for schema v0.2.0, which added ``target_raw`` and the
#: four condition components. Every ``record_id`` in the corpus changed, which
#: is the whole reason v0.2.0 is a version bump and not a patch; the migration
#: at ``syn4rch.curate.migrate`` carries the old-to-new map.
GOLDEN_ARGS = (
    "synthetic:test-001",
    "synthetic population alpha",
    "chemogenetic",
    "activation",
    "synthetic-actuator",
    None,                       # target_raw: a chemogenetic record has none
    "male",
    "chow",
    "young_adult",
    "ad_libitum",
    "synthetic readout one",
    "Fig S1A (synthetic)",
)
GOLDEN_ID = "d456f8ae94c75b93"

#: Twelve components joined by one separator.
EXPECTED_N_COMPONENTS = 12


class TestRecordIdContract(unittest.TestCase):
    """The declared field list, the function signature and the fixture agree.

    Three places state the hash input set, and a change to one of them that
    misses the others is exactly the kind of silent re-identification the
    golden digest exists to prevent. These tests make the three agree by
    assertion rather than by somebody remembering.
    """

    def test_signature_matches_declared_fields(self):
        params = tuple(inspect.signature(compute_record_id).parameters)
        self.assertEqual(params, RECORD_ID_FIELDS)

    def test_component_count_is_pinned(self):
        self.assertEqual(len(RECORD_ID_FIELDS), EXPECTED_N_COMPONENTS)

    def test_every_hashed_field_is_a_real_record_field(self):
        names = set(TierARecord.field_names())
        for field in RECORD_ID_FIELDS:
            self.assertIn(field, names, "%r is hashed but is not a field" % field)

    def test_golden_args_are_the_fixture_values(self):
        record = synthetic.good_record()
        actual = tuple(
            getattr(record, name) for name in RECORD_ID_FIELDS
        )
        normalised_actual = tuple(normalize_for_hash(v) for v in actual)
        normalised_golden = tuple(normalize_for_hash(v) for v in GOLDEN_ARGS)
        self.assertEqual(normalised_actual, normalised_golden)


class TestRecordIdDeterminism(unittest.TestCase):
    def test_golden_value_is_stable(self):
        self.assertEqual(compute_record_id(*GOLDEN_ARGS), GOLDEN_ID)

    def test_length_and_alphabet(self):
        rid = compute_record_id(*GOLDEN_ARGS)
        self.assertEqual(len(rid), RECORD_ID_LENGTH)
        self.assertTrue(all(c in "0123456789abcdef" for c in rid))

    def test_same_inputs_same_id_across_instances(self):
        a = synthetic.good_record()
        b = synthetic.good_record()
        self.assertEqual(a.record_id, b.record_id)
        self.assertEqual(a.record_id, GOLDEN_ID)

    def test_repeated_calls_are_identical(self):
        record = synthetic.good_record()
        self.assertEqual({record.compute_id() for _ in range(50)}, {record.record_id})

    def test_does_not_depend_on_python_hash_randomisation(self):
        """A subprocess with a different PYTHONHASHSEED must agree."""
        import os
        import subprocess
        import sys

        src = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"
        )
        code = (
            "import sys; sys.path.insert(0, %r);"
            "from syn4rch.schema import compute_record_id;"
            "print(compute_record_id(*%r))" % (src, GOLDEN_ARGS)
        )
        env = dict(os.environ, PYTHONHASHSEED="12345")
        out = subprocess.check_output([sys.executable, "-c", code], env=env)
        self.assertEqual(out.decode().strip(), GOLDEN_ID)

    def test_enum_and_string_inputs_agree(self):
        """Passing enums or their string values must hash identically."""
        from syn4rch.schema import (
            AgeGroup,
            Challenge,
            Diet,
            PerturbationMethod,
            Sex,
            Valence,
        )

        with_enums = compute_record_id(
            GOLDEN_ARGS[0],
            GOLDEN_ARGS[1],
            PerturbationMethod.chemogenetic,
            Valence.activation,
            GOLDEN_ARGS[4],
            GOLDEN_ARGS[5],
            Sex.male,
            Diet.chow,
            AgeGroup.young_adult,
            Challenge.ad_libitum,
            GOLDEN_ARGS[10],
            GOLDEN_ARGS[11],
        )
        self.assertEqual(with_enums, GOLDEN_ID)


class TestRecordIdSensitivity(unittest.TestCase):
    #: One perturbation per hashed field, as a record kwarg.
    #:
    #: ``target_raw`` on a chemogenetic record is not an admissible record
    #: (``S008_TARGET_ON_NON_MOLECULAR`` refuses it) but it is a legitimate
    #: hash input to test: these tests pin the *hash*, and the gate into
    #: ``data/curated/`` is ``validate``'s job, tested separately.
    FIELD_PERTURBATIONS = {
        "study_id": {"study_id": "synthetic:test-999"},
        "population_raw": {"population_raw": "synthetic population beta"},
        "method": {"method": "optogenetic"},
        "valence": {"valence": "inhibition"},
        "actuator": {"actuator": "synthetic-actuator-2"},
        "target_raw": {"target_raw": "Synthgene9"},
        "sex": {"sex": "female"},
        "diet": {"diet": "high_fat"},
        "age_group": {"age_group": "aged"},
        "challenge": {"challenge": "fasted"},
        "outcome_raw": {"outcome_raw": "synthetic readout two"},
        "source_locator": {"source_locator": "Fig S9Z (synthetic)"},
    }

    def test_every_hashed_field_changes_the_id(self):
        baseline = synthetic.good_record().record_id
        self.assertEqual(
            sorted(self.FIELD_PERTURBATIONS), sorted(RECORD_ID_FIELDS)
        )
        seen = {baseline}
        for name, override in self.FIELD_PERTURBATIONS.items():
            rid = synthetic.good_record(**override).record_id
            self.assertNotEqual(
                rid, baseline, "record_id ignored a change to %s" % name
            )
            self.assertNotIn(rid, seen, "collision after changing %s" % name)
            seen.add(rid)

    def test_the_v0_2_0_collision_class_is_closed(self):
        """Two cohorts differing only in a condition must not share an id.

        This is the defect v0.2.0 closed. Curation batch 04 produced three
        record_id collisions because a male and a female cohort, reported in
        one passage for one readout, hashed to the same value under v0.1.0.
        """
        base = synthetic.good_record(sex="male")
        other = synthetic.good_record(sex="female")
        self.assertNotEqual(base.record_id, other.record_id)
        chow = synthetic.good_record(diet="chow")
        hfd = synthetic.good_record(diet="high_fat")
        self.assertNotEqual(chow.record_id, hfd.record_id)
        fed = synthetic.good_record(challenge="ad_libitum")
        fasted = synthetic.good_record(challenge="fasted")
        self.assertNotEqual(fed.record_id, fasted.record_id)

    def test_two_genes_in_one_study_do_not_collide(self):
        """A knockout of one gene and of another are different observations."""
        a = synthetic.molecular_record(target_raw="SynA", target_id="syn:a")
        b = synthetic.molecular_record(target_raw="SynB", target_id="syn:b")
        self.assertNotEqual(a.record_id, b.record_id)

    def test_non_hashed_field_does_not_change_the_id(self):
        baseline = synthetic.good_record().record_id
        for override in (
            {"direction": "decrease"},
            {"effect_size_value": 99.0},
            {"n_animals": 99},
            {"targeting_strategy": "a completely different synthetic strategy"},
            {"extraction_notes": "different notes"},
            # Derived from method and valence, both already hashed, so adding
            # them to the hash would buy no discrimination (DECISIONS D-026).
            {"perturbation_class": "other"},
            {"perturbation_polarity": "unclear"},
            # target_id is the resolved form of target_raw, which is hashed.
            {"target_id": "syn:something_else"},
            {"challenge_detail": "different fabricated detail"},
            # Schema v0.3.0. ``condition_evidence`` says *why* a condition
            # holds the value it holds; the value itself is hashed and the
            # provenance is not. A record whose sex was read from Methods and
            # the same record whose sex was later traced to a figure legend
            # are one observation, so re-sourcing a condition must leave the
            # id alone -- otherwise filling in provenance would re-identify
            # the corpus and orphan every conflicts_with reference to it.
            # (This case replaces the retired condition_source scalar.)
            {
                "condition_evidence": [
                    {
                        "dimension": "sex",
                        "value": "male",
                        "source": "passage",
                        "locator": "Fig S1A legend (synthetic)",
                        "quote": "Synthetic quote: male mice, as above.",
                        "rule": "synthetic-other",
                    }
                ]
            },
            {"condition_evidence": []},
        ):
            self.assertEqual(
                synthetic.good_record(**override).record_id,
                baseline,
                "record_id should not depend on %s" % list(override),
            )

    def test_component_boundaries_cannot_be_forged(self):
        """Shifting text across a component boundary must change the id."""
        blanks = [""] * (EXPECTED_N_COMPONENTS - 2)
        a = compute_record_id("a", "b", *blanks)
        b = compute_record_id("ab", "", *blanks)
        self.assertNotEqual(a, b)

    def test_explicit_record_id_is_preserved(self):
        record = TierARecord(study_id="synthetic:test-001", record_id="deadbeefdeadbeef")
        self.assertEqual(record.record_id, "deadbeefdeadbeef")
        self.assertFalse(record.id_is_consistent())


class TestHashNormalisation(unittest.TestCase):
    """Normalisation behaviour, exercised through the full component list."""

    @staticmethod
    def _id(**overrides):
        args = dict(zip(RECORD_ID_FIELDS, ("s", "p", "m", "v", "a", "t",
                                           "sx", "d", "ag", "ch", "o", "l")))
        args.update(overrides)
        return compute_record_id(*(args[name] for name in RECORD_ID_FIELDS))

    def test_whitespace_is_collapsed(self):
        self.assertEqual(normalize_for_hash("  Fig   3C \n"), "Fig 3C")
        self.assertEqual(
            self._id(source_locator=" Fig  3C "),
            self._id(source_locator="Fig 3C"),
        )

    def test_unicode_is_nfkc_normalised(self):
        composed = "café"            # é as one code point
        decomposed = "café"         # e + combining acute
        self.assertNotEqual(composed, decomposed)
        self.assertEqual(
            self._id(population_raw=composed),
            self._id(population_raw=decomposed),
        )

    def test_case_is_preserved(self):
        self.assertNotEqual(
            self._id(population_raw="AgRP-like-synthetic"),
            self._id(population_raw="agrp-like-synthetic"),
        )

    def test_none_and_empty_string_agree(self):
        nones = {name: None for name in RECORD_ID_FIELDS}
        empties = {name: "" for name in RECORD_ID_FIELDS}
        nones["study_id"] = empties["study_id"] = "s"
        self.assertEqual(self._id(**nones), self._id(**empties))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
