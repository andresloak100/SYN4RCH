"""What a Tier A outcome is, made checkable (D-063).

The boundary existed only as one clause of D-054 -- "Tier A is systemic
physiology" -- until a reader asked whether EEG helps when it records the cortex
rather than the hypothalamus. The vocabulary now states the rule, and these
tests are what stop the next curator reaching for "theta power" on judgement.

The rule is NOT central versus peripheral. `inflammatory_markers_other`
deliberately admits microglial activation, Iba1, GFAP and gliosis because the
hypothalamic-inflammation literature reports them as the systemic-inflammation
readout, and a central/peripheral rule would have retroactively invalidated an
entry that was argued for. The rule is **state of the organism versus activity
of neurons**, and one of these tests pins the exception so a later sweep cannot
quietly remove it.
"""

import json
import os
import unittest

from _path import *  # noqa: F401,F403  (puts src/ on sys.path)

VOCAB = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "vocab", "outcomes.json",
)
SYNONYM_FIELDS = ("synonyms", "inverted_synonyms", "signed_synonyms")


def load():
    with open(VOCAB, encoding="utf-8") as handle:
        return json.load(handle)


def out_of_scope_terms(vocab):
    out = set()
    for key, terms in vocab["_meta"]["out_of_scope_readouts"].items():
        if key.startswith("_"):
            continue
        out |= {t.strip().lower() for t in terms}
    return out


class ScopePolicyTest(unittest.TestCase):
    def test_the_policy_exists_and_is_prose_a_curator_can_act_on(self):
        meta = load()["_meta"]
        self.assertIn("scope_policy", meta)
        self.assertGreaterEqual(len(meta["scope_policy"]), 3)
        joined = " ".join(meta["scope_policy"]).lower()
        for needle in ("state of the organism", "activity of neurons", "proxy"):
            self.assertIn(needle, joined)

    def test_the_policy_cites_the_schema_rather_than_a_preference(self):
        """The unit of record is what makes neural activity a perturbation."""
        joined = " ".join(load()["_meta"]["scope_policy"])
        self.assertIn("SCHEMA.md 1", joined)

    def test_the_policy_says_the_list_is_not_a_licence(self):
        """Absence from the exclusion list must not read as inclusion."""
        joined = " ".join(load()["_meta"]["scope_policy"]).lower()
        self.assertIn("not exhaustive", joined)
        self.assertIn("decision entry", joined)


class OutOfScopeTest(unittest.TestCase):
    def test_no_entry_claims_an_out_of_scope_readout_as_a_synonym(self):
        """The guard. A failure here means the vocabulary grew past its scope.

        Checked across every synonym field, because an out-of-scope readout
        smuggled in as an `inverted_synonym` would be curated just as readily as
        one smuggled in as a plain synonym.
        """
        vocab = load()
        banned = out_of_scope_terms(vocab)
        self.assertGreater(len(banned), 20)
        offences = [
            (outcome_id, field, term)
            for outcome_id, entry in vocab["entries"].items()
            for field in SYNONYM_FIELDS
            for term in (entry.get(field) or [])
            if term.strip().lower() in banned
        ]
        self.assertEqual(offences, [])

    def test_eeg_delta_power_is_gone_from_sleep(self):
        """The one deletion D-063 made, pinned by name.

        It invited a cortical spectral measure into the same claim cell as a
        sleep duration -- two quantities that are not even the same kind of
        thing, which is worse than the NREM/REM pooling the split addresses.
        """
        sleep = load()["entries"]["sleep"]
        for field in SYNONYM_FIELDS:
            self.assertNotIn("EEG delta power", sleep.get(field) or [])
        self.assertIn("eeg delta power", out_of_scope_terms(load()))

    def test_the_state_readouts_are_still_in_the_vocabulary_somewhere(self):
        """D-063 removed a spectral measure, not the sleep states themselves.

        Asserted across the whole file rather than against `sleep`, because
        D-064 then moved each state to an id of its own. The invariant that
        matters is that no state readout was lost between the two passes --
        excluding a cortical measure and relocating a vigilance state are
        different operations and only one of them happened here.
        """
        reachable = {
            term.lower()
            for entry in load()["entries"].values()
            for field in SYNONYM_FIELDS
            for term in (entry.get(field) or [])
        }
        for kept in (
            "nrem sleep",
            "rem sleep",
            "slow wave sleep",
            "wakefulness",
            "sleep latency",
            "sleep fragmentation",
            "total sleep time",
        ):
            self.assertIn(kept, reachable, kept)

    def test_the_behavioural_refusals_of_d054_are_listed_not_remembered(self):
        banned = out_of_scope_terms(load())
        for term in ("anxiety", "impulsivity", "progressive ratio break point"):
            self.assertIn(term, banned, term)

    def test_the_declared_categories_name_their_justification(self):
        note = load()["_meta"]["out_of_scope_readouts"]["_note"]
        self.assertIn("D-054", note)


class GlialExceptionTest(unittest.TestCase):
    """The documented exception, pinned so a later scope sweep cannot drop it."""

    def test_central_glial_markers_remain_in_scope(self):
        entry = load()["entries"]["inflammatory_markers_other"]
        synonyms = {s.lower() for s in entry["synonyms"]}
        for marker in ("microglial activation", "iba1", "gfap", "gliosis", "neuroinflammation"):
            self.assertIn(marker, synonyms, marker)

    def test_and_the_entry_still_says_why(self):
        """An undocumented exception is indistinguishable from an oversight."""
        notes = (load()["entries"]["inflammatory_markers_other"].get("notes") or "").lower()
        self.assertIn("central glial markers", notes)
        self.assertIn("systemic", notes)


class SleepSplitTest(unittest.TestCase):
    """The split D-064 made, and the pooling it ended.

    One `sleep` id held NREM, REM, wakefulness, latency and fragmentation as
    synonyms. NREM and REM are not complements and move independently, so a
    paper's NREM time and its wakefulness time landed in one claim cell -- the
    D-050 / D-052 defect class, a key too coarse to tell two experiments apart.
    """

    STATE_IDS = (
        "nrem_sleep",
        "rem_sleep",
        "wakefulness",
        "sleep_latency",
        "sleep_fragmentation",
    )

    def test_the_five_state_ids_exist_and_are_tier_a(self):
        entries = load()["entries"]
        for outcome_id in self.STATE_IDS:
            self.assertIn(outcome_id, entries)
            self.assertEqual(entries[outcome_id]["tier"], "A", outcome_id)

    def test_sleep_no_longer_claims_the_states_it_pooled(self):
        """Mandatory, not stylistic: vocab selfcheck forbids one normalised
        string mapping to two ids, so a term in both would fail to load."""
        sleep = load()["entries"]["sleep"]
        everything = {
            t.lower() for field in SYNONYM_FIELDS for t in (sleep.get(field) or [])
        }
        for moved in (
            "nrem sleep",
            "rem sleep",
            "slow wave sleep",
            "wakefulness",
            "sleep latency",
            "sleep fragmentation",
        ):
            self.assertNotIn(moved, everything, moved)

    def test_sleep_keeps_the_composite_readouts(self):
        """A sentence reporting two states jointly cannot be assigned to one."""
        sleep = load()["entries"]["sleep"]
        synonyms = {s.lower() for s in sleep["synonyms"]}
        for kept in ("sleep", "total sleep time", "sleep architecture", "sleep-wake amounts"):
            self.assertIn(kept, synonyms, kept)

    def test_bare_arousal_is_a_synonym_of_nothing(self):
        """Ambiguous between the waking state and a count of awakenings.

        The vocabulary annotates traps rather than resolving them (D-021), so an
        ambiguous term resolves to no id and forces a reading.
        """
        for entry in load()["entries"].values():
            for field in SYNONYM_FIELDS:
                self.assertNotIn(
                    "arousal",
                    [t.lower() for t in (entry.get(field) or [])],
                    entry["id"],
                )

    def test_the_unambiguous_arousal_senses_do_resolve(self):
        fragmentation = load()["entries"]["sleep_fragmentation"]
        synonyms = {s.lower() for s in fragmentation["synonyms"]}
        self.assertIn("arousal index", synonyms)
        self.assertIn("arousals", synonyms)

    def test_each_new_id_states_which_way_its_direction_runs(self):
        """These are the ids a curator will flip a sign on if left unstated."""
        entries = load()["entries"]
        for outcome_id in ("wakefulness", "sleep_latency", "sleep_fragmentation"):
            notes = (entries[outcome_id].get("notes") or "").lower()
            self.assertIn("direction", notes, outcome_id)
            self.assertIn("increase", notes, outcome_id)

    def test_fragmentation_warns_that_its_statistics_run_opposite_ways(self):
        """Bout count and bout duration move opposite ways for one change."""
        notes = load()["entries"]["sleep_fragmentation"]["notes"].lower()
        self.assertIn("opposite", notes)

    def test_nrem_and_rem_are_declared_non_complementary(self):
        notes = load()["entries"]["nrem_sleep"]["notes"].lower()
        self.assertIn("not a complement", notes)


class VersionTest(unittest.TestCase):
    def test_the_version_and_changelog_moved_together(self):
        meta = load()["_meta"]
        self.assertEqual(meta["version"], "0.5.1")
        self.assertTrue(
            any(entry.startswith("0.5.1:") for entry in meta["changelog"]),
            "a version bump with no changelog entry is an unexplained change",
        )

    def test_each_changelog_entry_cites_its_decision(self):
        changelog = load()["_meta"]["changelog"]
        for version, decision in (("0.4.0", "D-063"), ("0.5.0", "D-064"), ("0.5.1", "D-070")):
            entry = next(e for e in changelog if e.startswith(version + ":"))
            self.assertIn(decision, entry)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
