"""Evidence layers: literature counts are read, never invented; missing records read as zero."""

import _path  # noqa: F401
import json
import os
import tempfile
import unittest

from syn4rch.atlas import evidence

SUMMARY = """# x
| started (UTC) | 2026-01-01T00:00:00Z |
| finished (UTC) | 2026-01-01T00:02:00Z |
| unique papers after dedupe | 12 |

| population_id | curation_confidence | matched_in_text | of which `include` | found_by_its_queries |
|---|---|---|---|---|
| `arc_agrp` | high | 7 | 5 | 9 |
| `arc_tbx3` | low | 0 | 0 | 2 |

Populations with zero text matches
"""


class EvidenceTest(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def p(self, name):
        return os.path.join(self.d, name)

    def test_summary_table(self):
        with open(self.p("s.md"), "w") as fh:
            fh.write(SUMMARY)
        lit = evidence.literature(self.p("missing.jsonl"), self.p("s.md"))
        self.assertEqual(lit["source_kind"], "harvest_summary.md")
        self.assertEqual(lit["run_timestamp"], "2026-01-01T00:02:00Z")
        self.assertEqual(lit["counts"]["arc_agrp"], {"matched_in_text": 7, "include": 5})
        self.assertEqual(lit["counts"]["arc_tbx3"]["matched_in_text"], 0)
        self.assertEqual(lit["total_papers"], 12)
        self.assertIn("note", lit)

    def test_candidates_preferred_when_present(self):
        with open(self.p("c.jsonl"), "w") as fh:
            fh.write(json.dumps({"screen": {"decision": "include", "matched_populations": ["a", "b"]}}) + "\n")
            fh.write(json.dumps({"screen": {"decision": "exclude", "matched_populations": ["a"]}}) + "\n")
        lit = evidence.literature(self.p("c.jsonl"), self.p("nope.md"))
        self.assertEqual(lit["counts"]["a"], {"matched_in_text": 2, "include": 1})
        self.assertEqual(lit["total_papers"], 2)

    def test_summary_without_table_raises(self):
        with open(self.p("e.md"), "w") as fh:
            fh.write("# empty\n")
        with self.assertRaises(ValueError):
            evidence.literature_from_summary(self.p("e.md"))

    def test_absent_records_are_zero(self):
        c = evidence.causal_records(self.p("records.jsonl"))
        self.assertFalse(c["exists"])
        self.assertEqual(c["n_records"], 0)


if __name__ == "__main__":
    unittest.main()
