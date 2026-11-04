"""Polarity annotations on the outcomes vocabulary (DECISIONS.md D-019/D-021).

These guard the trap that produced the first real defect in the dataset: a
synonym matched the words in a paper, so a curator filed an inverted proxy at
face value and the direction came out backwards.
"""

import _path  # noqa: F401
import json
import os
import unittest

from syn4rch import vocab

OUTCOMES = os.path.join(vocab.vocab_dir(), "outcomes.json")

# The traps found so far. New entries may be added; these must not regress.
MUST_BE_INVERTED = {
    "insulin_sensitivity": {"insulin resistance", "HOMA-IR"},
    "glucose_tolerance": {"glucose excursion"},
    "sleep": {"sleep latency"},
    # Renamed from `cognitive_decline` at v0.2.0 (D-029). Its two synonyms were
    # classified as *signed* while the id named a quantity and its negation at
    # once; once the id names performance only, they are plainly inverted.
    "cognitive_performance": {"cognitive decline", "memory impairment",
                              "escape latency"},
}
MUST_BE_SIGNED = {
    "blood_glucose": {"hyperglycemia", "hypoglycemia"},
    "food_intake": {"hyperphagia", "hypophagia"},
    "core_temperature": {"hypothermia", "hyperthermia"},
    "blood_pressure": {"hypertension", "hypotension"},
    "locomotor_activity": {"hyperactivity", "hypoactivity"},
    "glucose_tolerance": {"impaired glucose tolerance"},
}


class PolarityAnnotationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(OUTCOMES, encoding="utf-8") as fh:
            cls.doc = json.load(fh)
        cls.entries = cls.doc["entries"]

    def test_known_inverted_synonyms_are_annotated(self):
        for oid, expected in MUST_BE_INVERTED.items():
            got = set(self.entries[oid].get("inverted_synonyms") or [])
            self.assertTrue(expected <= got, "%s missing %s" % (oid, expected - got))

    def test_known_signed_synonyms_are_annotated(self):
        for oid, expected in MUST_BE_SIGNED.items():
            got = set(self.entries[oid].get("signed_synonyms") or [])
            self.assertTrue(expected <= got, "%s missing %s" % (oid, expected - got))

    def test_every_annotated_synonym_really_is_a_synonym(self):
        # An annotation naming a string the entry does not list would never fire.
        for oid, entry in self.entries.items():
            syns = set(entry.get("synonyms") or [])
            for field in ("inverted_synonyms", "signed_synonyms"):
                for s in entry.get(field) or []:
                    self.assertIn(s, syns, "%s.%s lists %r, not in synonyms" % (oid, field, s))

    def test_a_synonym_is_never_both_inverted_and_signed(self):
        for oid, entry in self.entries.items():
            inv = set(entry.get("inverted_synonyms") or [])
            sig = set(entry.get("signed_synonyms") or [])
            self.assertEqual(inv & sig, set(), "%s: %s in both" % (oid, inv & sig))

    def test_polarity_policy_is_recorded_and_forbids_automatic_rewriting(self):
        policy = " ".join(self.doc["_meta"]["polarity_policy"]).lower()
        self.assertIn("d-019", policy)
        self.assertIn("opposite", policy)
        self.assertIn("never be used to rewrite a direction automatically", policy)

    def test_itt_and_gtt_entries_point_raw_glucose_elsewhere(self):
        # D-019 step 1: the literally measured quantity wins when it exists.
        for oid in ("insulin_sensitivity", "glucose_tolerance"):
            notes = (self.entries[oid].get("notes") or "").lower()
            self.assertIn("blood_glucose", notes, oid)
            self.assertIn("d-019", notes, oid)

    def test_the_ambiguous_outcome_id_was_renamed(self):
        """D-029. The old id named a quantity and its negation at once."""
        self.assertNotIn("cognitive_decline", self.entries)
        self.assertIn("cognitive_performance", self.entries)
        entry = self.entries["cognitive_performance"]
        self.assertEqual(entry["label"], "Cognitive performance")
        notes = (entry.get("notes") or "").lower()
        self.assertIn("d-029", notes)
        self.assertIn("increase = better", notes)
        # The wording that used to be filed as merely "signed" is inverted:
        # a paper reporting less decline is reporting more performance.
        self.assertNotIn("signed_synonyms", entry)

    def test_a_renamed_id_still_resolves_the_old_wording(self):
        """The rename must not lose a paper's own words."""
        for wording in ("cognitive decline", "memory impairment",
                        "cognitive performance", "Morris water maze"):
            self.assertEqual(
                vocab.resolve("outcomes", wording), "cognitive_performance",
                wording,
            )

    def test_the_rename_is_recorded_in_the_changelog(self):
        """The rename must stay recorded as later entries are appended.

        This used to assert ``version == "0.2.0"``, which pinned the version
        the rename happened at and broke on the next legitimate bump (0.3.0,
        D-049) while the guarantee it is named for -- that the rename is
        written down -- was never in danger. The version is now checked for
        shape only; what matters is that the changelog still carries the
        rename, whatever has been appended after it.
        """
        changelog = " ".join(self.doc["_meta"].get("changelog") or []).lower()
        self.assertIn("cognitive_decline", changelog)
        self.assertIn("cognitive_performance", changelog)
        version = self.doc["_meta"]["version"]
        self.assertRegex(version, r"^\d+\.\d+\.\d+$")

    def test_annotations_do_not_disturb_synonym_resolution(self):
        # The new fields are metadata; matching still works off `synonyms`.
        self.assertEqual(vocab.resolve("outcomes", "insulin resistance"), "insulin_sensitivity")
        self.assertEqual(vocab.resolve("outcomes", "hyperphagia"), "food_intake")
        self.assertEqual(vocab.resolve("outcomes", "blood glucose"), "blood_glucose")

    def test_vocabulary_still_loads_and_selfchecks(self):
        self.assertEqual(len(vocab.load("outcomes")), len(self.entries))
        vocab.selfcheck()


if __name__ == "__main__":
    unittest.main()
