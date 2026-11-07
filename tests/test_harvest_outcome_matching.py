"""Outcome matching after the sleep split (D-070).

Two defects found while costing the EEG curation pass, both invisible until the
matcher was run over the passages it would actually be used on:

1. "NREM and REM sleep both decreased" matched `rem_sleep` and missed
   `nrem_sleep`, because the matcher is a word-boundary alternation and neither
   `NREM sleep` nor `non-REM sleep` occurs contiguously in that phrase. 23 of
   946 proposed sleep units.
2. "EEG delta power during wakefulness increased" resolves to `wakefulness`
   through an incidental word, so a readout D-063 put out of scope arrives
   under a legitimate id. 142 of 946. `tests/test_vocab_scope.py` guards the
   vocabulary and could not catch this: no rule is broken on the way in.
"""

import json
import os
import unittest

from _path import *  # noqa: F401,F403  (puts src/ on sys.path)

from syn4rch.harvest.screen import (
    _norm_text,
    matchers,
    out_of_scope_readouts,
)
from syn4rch.vocab.core import load, meta

STATES = ("nrem_sleep", "rem_sleep", "wakefulness", "sleep_latency", "sleep_fragmentation")


def ids(text):
    return set(matchers()["outcomes"].match(text, _norm_text(text)))


class JointPhrasingTest(unittest.TestCase):
    """The commonest composite phrasing in this literature, and it lost a half."""

    def test_nrem_and_rem_matches_both_states(self):
        for phrase in (
            "NREM and REM sleep both decreased",
            "non-REM and REM sleep decreased",
            "REM and NREM sleep were unchanged",
        ):
            got = ids(phrase)
            self.assertIn("nrem_sleep", got, phrase)
            self.assertIn("rem_sleep", got, phrase)

    def test_the_composite_id_is_still_offered_alongside(self):
        """Not a defect: under D-064 a joint report belongs on `sleep`.

        The matcher proposes candidates and the curator decides which the record
        is about. Offering `sleep` beside the two states is what lets a reader
        follow D-064's refusal rule; suppressing it would make that choice for
        them.
        """
        self.assertIn("sleep", ids("NREM and REM sleep both decreased"))

    def test_each_state_alone_still_resolves(self):
        self.assertIn("nrem_sleep", ids("NREM sleep time was shorter"))
        self.assertIn("rem_sleep", ids("REM sleep was longer"))
        self.assertIn("rem_sleep", ids("rapid eye movement sleep rose"))

    def test_the_elided_forms_are_in_the_vocabulary_not_the_code(self):
        """The fix is data, so it is reviewable in a diff."""
        synonyms = load("outcomes")["nrem_sleep"]["synonyms"]
        self.assertIn("NREM", synonyms)
        self.assertIn("non-REM", synonyms)


class BareRemTest(unittest.TestCase):
    """Bare `REM` is listed, and the overlap it brings was measured first.

    I first withheld it, reasoning only from the false positive: `\\bREM\\b`
    matches inside `non-REM sleep`. That was the wrong frame for a *candidate*
    proposer. A false candidate is reviewable by the curator who sees it; a
    missed one is invisible. Measured before accepting: 25 units mention
    non-REM and all 25 also discuss REM in its own right, so the spurious match
    never fires alone on this corpus, while withholding the synonym lost REM on
    12 units outright.
    """

    def test_bare_rem_is_a_synonym_so_the_reversed_ordering_resolves(self):
        self.assertIn("REM", load("outcomes")["rem_sleep"]["synonyms"])
        got = ids("REM and NREM sleep were unchanged")
        self.assertIn("rem_sleep", got)
        self.assertIn("nrem_sleep", got)

    def test_the_non_rem_overlap_is_pinned_rather_than_patched(self):
        """Left as a known, measured overlap instead of a negative lookbehind.

        A lookbehind in `_boundary_pattern` would apply to every vocabulary,
        which is a wide change to fix a narrow trap that costs nothing here.
        This pins the behaviour so a future change to it is deliberate.
        """
        got = ids("non-REM sleep was longer")
        self.assertIn("nrem_sleep", got)
        self.assertIn("rem_sleep", got)

    def test_a_plain_nrem_passage_names_nrem(self):
        self.assertIn("nrem_sleep", ids("NREM sleep time was shorter"))


class OutOfScopeDetectorTest(unittest.TestCase):
    def test_it_names_the_category_of_an_out_of_scope_readout(self):
        self.assertEqual(
            out_of_scope_readouts("EEG delta power during wakefulness increased"),
            ["neural_activity"],
        )
        self.assertEqual(
            out_of_scope_readouts("anxiety-like behaviour in the open field"),
            ["behavioural_affective"],
        )

    def test_it_is_silent_on_a_clean_state_passage(self):
        self.assertEqual(out_of_scope_readouts("NREM sleep time was shorter"), [])

    def test_it_catches_the_readouts_d063_excluded_by_name(self):
        for text, category in (
            ("sleep spindles were denser", "neural_activity"),
            ("c-Fos counts rose", "neural_activity"),
            ("the local field potential shifted", "neural_activity"),
            ("firing rate increased", "neural_activity"),
            ("progressive ratio break point fell", "behavioural_affective"),
        ):
            self.assertIn(category, out_of_scope_readouts(text), text)

    def test_it_does_not_suppress_the_in_scope_id(self):
        """Deliberate. A sentence can report both, and choosing is a reading.

        D-021: the vocabulary annotates traps rather than resolving them,
        because the fix is a judgement. The detector prompts the judgement; it
        does not make it.
        """
        text = "EEG delta power during wakefulness increased"
        self.assertIn("wakefulness", ids(text))
        self.assertTrue(out_of_scope_readouts(text))

    def test_every_declared_term_is_matchable_not_one_per_category(self):
        """Sampling a category proves it is wired, not that its members are.

        This checked the first term of each category until the genotype session
        observed that both of our near-miss audits had looked in one direction
        only -- it then found `arc_ghrh` losing a cluster at 0.583 under a 0.60
        cut, a removal-direction failure neither audit could have caught. The
        same one-sidedness was here.

        It matters because `TermMatcher` silently drops terms shorter than two
        characters and acronyms listed in `AMBIGUOUS_ACRONYMS` -- which holds
        `ARC`, `LH`, `ME`, `ZI` and `HY`, because they collide with populations
        and regions. So a term could be declared out of scope and undetectable
        at once, which would make it one more declared thing nothing reads.

        Measured when written: 39 declared, 0 undetectable.
        """
        block = meta("outcomes")["out_of_scope_readouts"]
        categories = [k for k in block if not k.startswith("_")]
        self.assertTrue(categories)
        undetectable = [
            (category, term)
            for category in categories
            for term in block[category]
            if category not in out_of_scope_readouts(term)
        ]
        self.assertEqual(
            undetectable,
            [],
            "declared out-of-scope term(s) the matcher cannot detect; "
            "TermMatcher drops sub-2-character terms and ambiguous acronyms",
        )

    def test_the_declared_list_is_large_enough_to_be_a_rule(self):
        block = meta("outcomes")["out_of_scope_readouts"]
        total = sum(len(v) for k, v in block.items() if not k.startswith("_"))
        self.assertGreaterEqual(total, 30)


class MetaAccessorTest(unittest.TestCase):
    def test_meta_exposes_the_policy_blocks_load_hides(self):
        """Added so a consumer does not reimplement the loader's path resolution."""
        block = meta("outcomes")
        self.assertIn("out_of_scope_readouts", block)
        self.assertIn("scope_policy", block)
        self.assertEqual(block["version"], "0.5.2")

    def test_meta_is_empty_rather_than_raising_for_a_file_without_one(self):
        for kind in ("regions", "methods", "populations", "targets"):
            self.assertIsInstance(meta(kind), dict)




class UnitsArtefactFreshnessTest(unittest.TestCase):
    """Every outcome id in the units file must exist in the vocabulary (D-072).

    The artefact carried `cognitive_decline` on 111 units for the life of the
    project -- an id renamed to `cognitive_performance` at outcomes 0.2.0 by
    D-029. Regenerating replaced all 111 with the new id, so no readout was
    ever lost; the rename simply never reached the artefact and nothing checked.

    Same class as the stale committed signaling graph: an invariant about a
    generated file that no generator can enforce after the fact, because the
    generator is not what went stale.
    """

    UNITS = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "data", "interim", "extraction_units.jsonl",
    )

    def test_no_unit_carries_an_outcome_id_the_vocabulary_lacks(self):
        if not os.path.exists(self.UNITS):
            self.skipTest("data/interim is untracked machine output (D-016)")
        known = set(load("outcomes"))
        seen, stale = set(), set()
        with open(self.UNITS, encoding="utf-8") as handle:
            for line in handle:
                outcome = json.loads(line).get("candidate_outcome_id")
                if not outcome:
                    continue
                seen.add(outcome)
                if outcome not in known:
                    stale.add(outcome)
        self.assertTrue(seen, "no outcome ids in the units file at all")
        self.assertEqual(
            sorted(stale),
            [],
            "units carry outcome id(s) the vocabulary no longer has; "
            "regenerate rather than edit the artefact",
        )

    def test_the_split_ids_reached_the_artefact(self):
        """The regeneration is what makes the EEG pass safe; pin that it ran."""
        if not os.path.exists(self.UNITS):
            self.skipTest("data/interim is untracked machine output (D-016)")
        with open(self.UNITS, encoding="utf-8") as handle:
            present = {
                json.loads(line).get("candidate_outcome_id") for line in handle
            }
        for state in ("nrem_sleep", "rem_sleep", "wakefulness"):
            self.assertIn(state, present, state)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
