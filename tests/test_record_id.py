"""record_id must be deterministic, portable, and sensitive to its inputs."""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.schema import RECORD_ID_FIELDS, TierARecord, compute_record_id
from syn4rch.schema.record import RECORD_ID_LENGTH, normalize_for_hash

#: Golden digest. Pinned so that a change to the hash input set, the
#: separator, the normalisation or the digest length fails loudly instead of
#: silently re-identifying the whole corpus. Recomputing this value is a
#: schema change and needs a version bump.
GOLDEN_ARGS = (
    "synthetic:test-001",
    "synthetic population alpha",
    "chemogenetic",
    "activation",
    "synthetic-actuator",
    "synthetic readout one",
    "Fig S1A (synthetic)",
)
GOLDEN_ID = "95a8a980e8695a08"


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
        from syn4rch.schema import PerturbationMethod, Valence

        with_enums = compute_record_id(
            GOLDEN_ARGS[0],
            GOLDEN_ARGS[1],
            PerturbationMethod.chemogenetic,
            Valence.activation,
            GOLDEN_ARGS[4],
            GOLDEN_ARGS[5],
            GOLDEN_ARGS[6],
        )
        self.assertEqual(with_enums, GOLDEN_ID)


class TestRecordIdSensitivity(unittest.TestCase):
    #: One perturbation per hashed field, as a record kwarg.
    FIELD_PERTURBATIONS = {
        "study_id": {"study_id": "synthetic:test-999"},
        "population_raw": {"population_raw": "synthetic population beta"},
        "method": {"method": "optogenetic"},
        "valence": {"valence": "inhibition"},
        "actuator": {"actuator": "synthetic-actuator-2"},
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

    def test_non_hashed_field_does_not_change_the_id(self):
        baseline = synthetic.good_record().record_id
        for override in (
            {"direction": "decrease"},
            {"effect_size_value": 99.0},
            {"n_animals": 99},
            {"targeting_strategy": "a completely different synthetic strategy"},
            {"extraction_notes": "different notes"},
        ):
            self.assertEqual(
                synthetic.good_record(**override).record_id,
                baseline,
                "record_id should not depend on %s" % list(override),
            )

    def test_component_boundaries_cannot_be_forged(self):
        """Shifting text across a component boundary must change the id."""
        a = compute_record_id("a", "b", "", "", "", "", "")
        b = compute_record_id("ab", "", "", "", "", "", "")
        self.assertNotEqual(a, b)

    def test_explicit_record_id_is_preserved(self):
        record = TierARecord(study_id="synthetic:test-001", record_id="deadbeefdeadbeef")
        self.assertEqual(record.record_id, "deadbeefdeadbeef")
        self.assertFalse(record.id_is_consistent())


class TestHashNormalisation(unittest.TestCase):
    def test_whitespace_is_collapsed(self):
        self.assertEqual(normalize_for_hash("  Fig   3C \n"), "Fig 3C")
        self.assertEqual(
            compute_record_id("s", "p", "m", "v", "a", "o", " Fig  3C "),
            compute_record_id("s", "p", "m", "v", "a", "o", "Fig 3C"),
        )

    def test_unicode_is_nfkc_normalised(self):
        composed = "café"          # é as one code point
        decomposed = "café"       # e + combining acute
        self.assertNotEqual(composed, decomposed)
        self.assertEqual(
            compute_record_id("s", composed, "m", "v", "a", "o", "l"),
            compute_record_id("s", decomposed, "m", "v", "a", "o", "l"),
        )

    def test_case_is_preserved(self):
        self.assertNotEqual(
            compute_record_id("s", "AgRP-like-synthetic", "m", "v", "a", "o", "l"),
            compute_record_id("s", "agrp-like-synthetic", "m", "v", "a", "o", "l"),
        )

    def test_none_and_empty_string_agree(self):
        self.assertEqual(
            compute_record_id("s", None, None, None, None, None, None),
            compute_record_id("s", "", "", "", "", "", ""),
        )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
