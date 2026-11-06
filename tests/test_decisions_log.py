"""Invariants of the append-only decision log.

`docs/DECISIONS.md` is the project's citable record: code comments, reports and
other entries refer to decisions by number, so the number is an identifier and
not a position. These tests assert the properties that make a citation
resolvable, and deliberately do **not** assert one that looks obvious.

Both failures these tests exist to prevent happened, within a day of each other,
between two sessions working the same checkout:

* a code comment cited `D-061` while `DECISIONS.md` ended at `D-060`, inverting
  the project's own rule that code cites a decision and the decision explains
  the code;
* two entries sat uncommitted in the file at once, so a commit staging it by
  path took both, and the numbers came out non-monotonic.
"""

import _path  # noqa: F401
import os
import re
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOG = os.path.join(ROOT, "docs", "DECISIONS.md")

#: A heading: the identifier, then a title. Both are required -- a numbered
#: entry with no title cannot be found by a reader who only half-remembers it.
HEADING = re.compile(r"^## D-(\d{3}) - (.+)$", re.M)

#: A citation anywhere else. Bounded on the right so "D-0481" is not read as
#: "D-048", and case-sensitive because the project writes them one way.
CITATION = re.compile(r"\bD-(\d{3})\b")

#: Where a citation may appear. Reports are excluded: they are regenerated
#: outputs, so a stale citation in one is a signal to regenerate rather than a
#: defect in the log.
CITING_TREES = ("src", "docs", "tests")


def log_text() -> str:
    with open(LOG, encoding="utf-8") as handle:
        return handle.read()


def numbers():
    return [int(n) for n, _ in HEADING.findall(log_text())]


def citations():
    """``{number: [paths]}`` for every D-NNN cited outside the log itself."""
    found = {}
    for tree in CITING_TREES:
        base = os.path.join(ROOT, tree)
        for dirpath, _dirs, files in os.walk(base):
            if "__pycache__" in dirpath:
                continue
            for name in files:
                if not name.endswith((".py", ".md", ".json")):
                    continue
                path = os.path.join(dirpath, name)
                if os.path.abspath(path) == os.path.abspath(LOG):
                    continue
                try:
                    with open(path, encoding="utf-8") as handle:
                        text = handle.read()
                except (OSError, UnicodeDecodeError):  # pragma: no cover
                    continue
                for n in CITATION.findall(text):
                    found.setdefault(int(n), []).append(
                        os.path.relpath(path, ROOT)
                    )
    return found


class HeadingTest(unittest.TestCase):
    def test_the_log_has_entries(self):
        self.assertGreater(len(numbers()), 50)

    def test_every_entry_has_a_number_and_a_title(self):
        raw = [
            line for line in log_text().splitlines()
            if line.startswith("## D-")
        ]
        self.assertEqual(
            len(raw),
            len(HEADING.findall(log_text())),
            "an entry heading does not match '## D-NNN - title'",
        )


class NumberingTest(unittest.TestCase):
    def test_no_number_is_claimed_twice(self):
        """Two entries under one number make every citation to it ambiguous."""
        seen = numbers()
        dupes = sorted({n for n in seen if seen.count(n) > 1})
        self.assertEqual(dupes, [], "duplicate decision numbers: %s" % dupes)

    def test_there_are_no_gaps(self):
        """A gap means an entry was lost, or a number was burned silently.

        Either is worth knowing: a reader who cannot find D-0NN should be able
        to conclude it never existed rather than wonder.
        """
        seen = set(numbers())
        missing = sorted(set(range(min(seen), max(seen) + 1)) - seen)
        self.assertEqual(missing, [], "missing decision numbers: %s" % missing)

    def test_the_numbers_are_deliberately_NOT_required_to_ascend(self):
        """The file order is chronological; the numbers are identifiers.

        Two sessions writing concurrently can append in an order that does not
        match their reserved numbers -- which is what happened to D-061 and
        D-062. The append order is the truthful one: D-062 *was* written first.
        Reordering the body to make the numbers ascend would rewrite an
        append-only file to make it look tidier, which is the worse trade.

        So this test asserts only that the inversion is bounded. If the log ever
        becomes substantially unordered, that is a coordination failure rather
        than a cosmetic one, and this catches it.
        """
        seen = numbers()
        inversions = [(a, b) for a, b in zip(seen, seen[1:]) if b < a]
        self.assertLessEqual(
            len(inversions),
            1,
            "more than one out-of-order adjacency: %s" % inversions,
        )


class CitationTest(unittest.TestCase):
    def test_every_cited_decision_exists(self):
        """The failure this test was written for.

        A code comment cited D-061 while the log ended at D-060. The project's
        rule is that code cites a decision and the decision explains the code;
        a citation to an entry that does not exist inverts it, and the reader
        who follows the pointer finds nothing.
        """
        existing = set(numbers())
        cited = citations()
        dangling = {
            n: sorted(set(paths))
            for n, paths in cited.items()
            if n not in existing
        }
        self.assertEqual(
            dangling,
            {},
            "citations to decisions that do not exist: %s" % dangling,
        )

    def test_the_citation_pattern_does_not_over_match(self):
        self.assertEqual(CITATION.findall("see D-048 and D-056."), ["048", "056"])
        self.assertEqual(CITATION.findall("D-0481"), [])
        self.assertEqual(CITATION.findall("A008-A011"), [])

    def test_most_entries_are_cited_somewhere(self):
        """Not all of them: an entry may record a decision no code depends on.

        A very low rate would mean the log had drifted away from the code it is
        supposed to explain, so this is a floor rather than a requirement.
        """
        existing = set(numbers())
        cited = set(citations()) & existing
        self.assertGreater(
            len(cited),
            len(existing) // 3,
            "only %d of %d entries are cited anywhere" % (len(cited), len(existing)),
        )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
