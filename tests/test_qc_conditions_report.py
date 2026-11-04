"""The condition-evidence section of the QC report (schema v0.3.0).

What this section must answer, and what these tests pin:

* how many condition values are *evidenced* -- a verbatim sentence was found
  and recorded -- as opposed to merely *known*, which is what
  ``condition_specificity`` already counts;
* which extraction rule produced them, and on how many records, so that a
  rule found to be wrong can be costed before it is fixed;
* how far across the corpus's studies the evidence reaches.

The three states a dimension can be in are counted apart on purpose.
``not_reported`` is an assertion that the paper says nothing, so it has no
sentence to quote and no evidence entry by construction; adding it to the
evidenced count would manufacture provenance. Several tests below exist only
to keep that sum from ever being taken.
"""

import _path  # noqa: F401
import json
import unittest

import synthetic
from syn4rch.qc import render_text, run_qc
from syn4rch.qc.keys import NULL

#: Title line of the section under test, as rendered.
TITLE = "Condition evidence (what is quoted, not what is known)"

#: Title line of the *other* condition section, which must keep working.
SPECIFICITY_TITLE = "Condition specificity"


def evidence(dimension, value, source="methods", rule="synthetic"):
    """One fabricated condition evidence entry. Nothing here is from a paper."""
    return {
        "dimension": dimension,
        "value": value,
        "source": source,
        "locator": "Methods para 1" if source == "methods" else "Results para 1",
        "quote": "Synthetic quote: fabricated sentence for %s." % dimension,
        "rule": rule,
    }


def section_lines(text):
    """The rendered section, from its title to the blank line that ends it."""
    lines = text.split("\n")
    start = lines.index(TITLE)
    out = []
    for line in lines[start:]:
        if not line.strip():
            break
        out.append(line)
    return out


class ConditionEvidenceCountsTest(unittest.TestCase):
    def setUp(self):
        self.loader = synthetic.stub_vocab_loader()

    def stats(self, records):
        return run_qc(records, vocab_loader=self.loader).stats["condition_evidence"]

    def test_a_fully_evidenced_dataset_has_no_gap(self):
        stats = self.stats([synthetic.good_record(), synthetic.other_lab_record()])
        for name in ("sex", "diet", "age_group", "challenge"):
            self.assertEqual(stats["dimensions"][name]["positive"], 2, name)
            self.assertEqual(stats["dimensions"][name]["evidenced"], 2, name)
        self.assertEqual(stats["gap"], {"records": 0, "values": 0})
        self.assertEqual(stats["n_entries"], 8)

    def test_a_dataset_with_no_evidence_is_entirely_gap(self):
        """The state of all 362 records curated before v0.3.0."""
        records = [
            synthetic.good_record(condition_evidence=[]),
            synthetic.other_lab_record(condition_evidence=[]),
        ]
        stats = self.stats(records)
        for name in ("sex", "diet", "age_group", "challenge"):
            self.assertEqual(stats["dimensions"][name]["positive"], 2, name)
            self.assertEqual(stats["dimensions"][name]["evidenced"], 0, name)
        self.assertEqual(stats["gap"], {"records": 2, "values": 8})
        self.assertEqual(stats["n_entries"], 0)
        self.assertEqual(stats["by_rule"], {})

    def test_the_gap_record_count_is_the_w017_population(self):
        """Two numbers for one thing must not be allowed to drift apart."""
        records = [
            synthetic.good_record(condition_evidence=[]),
            synthetic.other_lab_record(condition_evidence=[]),
            synthetic.same_lab_other_study_record(),  # fully evidenced
        ]
        report = run_qc(records, vocab_loader=self.loader)
        w017 = [
            f for f in report.findings
            if f.code == "W017_CONDITION_WITHOUT_EVIDENCE"
        ]
        self.assertEqual(
            report.stats["condition_evidence"]["gap"]["records"], len(w017)
        )
        self.assertEqual(len(w017), 2)

    def test_not_reported_is_neither_evidenced_nor_a_gap(self):
        """An absence cannot be quoted, so it is not a missing quote either."""
        record = synthetic.good_record(
            diet="not_reported",
            condition_evidence=[
                evidence("sex", "male"),
                evidence("age_group", "young_adult"),
                evidence("challenge", "ad_libitum", source="passage"),
            ],
        )
        stats = self.stats([record])
        diet = stats["dimensions"]["diet"]
        self.assertEqual(diet["not_reported"], 1)
        self.assertEqual(diet["positive"], 0)
        self.assertEqual(diet["evidenced"], 0)
        self.assertEqual(stats["gap"], {"records": 0, "values": 0})

    def test_null_is_neither_evidenced_nor_a_gap(self):
        record = synthetic.good_record(
            challenge=None,
            challenge_detail=None,
            condition_evidence=[
                evidence("sex", "male"),
                evidence("diet", "chow"),
                evidence("age_group", "young_adult"),
            ],
        )
        stats = self.stats([record])
        challenge = stats["dimensions"]["challenge"]
        self.assertEqual(challenge["null"], 1)
        self.assertEqual(challenge["positive"], 0)
        self.assertEqual(challenge["evidenced"], 0)
        self.assertEqual(stats["gap"], {"records": 0, "values": 0})

    def test_every_dimension_row_accounts_for_every_record(self):
        records = [
            synthetic.good_record(),
            synthetic.other_lab_record(diet="not_reported"),
            synthetic.same_lab_other_study_record(
                sex=None, challenge=None, challenge_detail=None
            ),
        ]
        stats = self.stats(records)
        for name, row in stats["dimensions"].items():
            self.assertEqual(
                row["positive"] + row["not_reported"] + row["null"],
                stats["_total_records"],
                name,
            )
            self.assertLessEqual(row["evidenced"], row["positive"], name)

    def test_an_entry_for_an_unfilled_dimension_is_not_a_covered_value(self):
        """Provenance pointing at a value that is not there (schema A011)."""
        record = synthetic.good_record(
            diet=None,
            condition_evidence=[
                evidence("sex", "male"),
                evidence("diet", "chow"),  # the field is null
                evidence("age_group", "young_adult"),
                evidence("challenge", "ad_libitum", source="passage"),
            ],
        )
        stats = self.stats([record])
        self.assertEqual(stats["dimensions"]["diet"]["null"], 1)
        self.assertEqual(stats["dimensions"]["diet"]["evidenced"], 0)
        self.assertEqual(stats["n_entries"], 4, "the stray entry is still counted")

    def test_sources_are_counted_apart(self):
        stats = self.stats([synthetic.good_record()])
        self.assertEqual(stats["by_source"]["methods"], 3)
        self.assertEqual(stats["by_source"]["passage"], 1)

    def test_the_rule_table_counts_entries_and_records_and_orders_by_frequency(self):
        records = [
            synthetic.good_record(
                condition_evidence=[
                    evidence("sex", "male", rule="methods_header"),
                    evidence("diet", "chow", rule="methods_header"),
                    evidence("age_group", "young_adult", rule="methods_header"),
                    evidence(
                        "challenge", "ad_libitum", source="passage", rule="passage_scan"
                    ),
                ]
            ),
            synthetic.other_lab_record(
                condition_evidence=[
                    evidence("sex", "male", source="passage", rule="passage_scan"),
                    evidence("diet", "chow", source="passage", rule="passage_scan"),
                    evidence(
                        "age_group", "young_adult", source="passage", rule="passage_scan"
                    ),
                    evidence(
                        "challenge", "ad_libitum", source="passage", rule="passage_scan"
                    ),
                ]
            ),
        ]
        by_rule = self.stats(records)["by_rule"]
        self.assertEqual(list(by_rule), ["passage_scan", "methods_header"])
        self.assertEqual(by_rule["passage_scan"], {"entries": 5, "records": 2})
        # Three entries on one record is one record's worth of exposure, which
        # is the number an audit of a wrong rule needs.
        self.assertEqual(by_rule["methods_header"], {"entries": 3, "records": 1})

    def test_an_entry_with_no_rule_is_visible_rather_than_dropped(self):
        record = synthetic.good_record(
            condition_evidence=[evidence("sex", "male", rule=None)],
            diet=None,
            age_group=None,
            age_value=None,
            age_unit=None,
            challenge=None,
            challenge_detail=None,
        )
        self.assertEqual(self.stats([record])["by_rule"], {NULL: {"entries": 1, "records": 1}})

    def test_study_reach_distinguishes_concentration_from_spread(self):
        """Same coverage counts, two different corpora."""
        concentrated = [
            synthetic.good_record(),
            synthetic.good_record(source_locator="Fig S1B (synthetic)"),
            synthetic.other_lab_record(condition_evidence=[]),
        ]
        spread = [
            synthetic.good_record(),
            synthetic.other_lab_record(),
            synthetic.same_lab_other_study_record(condition_evidence=[]),
        ]
        a, b = self.stats(concentrated), self.stats(spread)
        self.assertEqual(
            a["dimensions"], b["dimensions"], "the two must be indistinguishable here"
        )
        self.assertEqual(a["gap"], b["gap"])
        self.assertEqual(a["study_reach"]["n_studies_with_evidenced_value"], 1)
        self.assertEqual(a["study_reach"]["n_studies"], 2)
        self.assertEqual(b["study_reach"]["n_studies_with_evidenced_value"], 2)
        self.assertEqual(b["study_reach"]["n_studies"], 3)

    def test_evidence_on_a_record_with_no_study_id_is_not_silently_lost(self):
        records = [synthetic.good_record(study_id=None)]
        reach = self.stats(records)["study_reach"]
        self.assertEqual(reach["n_studies"], 0)
        self.assertEqual(reach["n_studies_with_evidenced_value"], 0)
        self.assertEqual(reach["n_evidenced_records_without_study_id"], 1)

    def test_stats_survive_a_json_round_trip_with_the_rule_order_intact(self):
        records = [
            synthetic.good_record(
                condition_evidence=[
                    evidence("sex", "male", rule="methods_header"),
                    evidence("diet", "chow", rule="methods_header"),
                    evidence("age_group", "young_adult", rule="other_rule"),
                    evidence(
                        "challenge", "ad_libitum", source="passage", rule="other_rule"
                    ),
                ]
            ),
            synthetic.other_lab_record(
                condition_evidence=[evidence("sex", "male", rule="methods_header")]
            ),
        ]
        report = run_qc(records, vocab_loader=self.loader)
        reloaded = json.loads(json.dumps(report.to_dict()))
        stats = reloaded["stats"]["condition_evidence"]
        self.assertEqual(list(stats["by_rule"]), ["methods_header", "other_rule"])
        self.assertEqual(stats["by_rule"]["methods_header"]["records"], 2)


class ConditionEvidenceSectionTest(unittest.TestCase):
    def setUp(self):
        self.loader = synthetic.stub_vocab_loader()

    def render(self, records):
        return render_text(run_qc(records, vocab_loader=self.loader))

    def test_the_section_refuses_to_add_not_reported_to_evidenced(self):
        text = self.render([synthetic.good_record(diet="not_reported")])
        self.assertIn(TITLE, text)
        joined = " ".join(section_lines(text))
        self.assertIn("verbatim sentence was found and checked", joined)
        self.assertIn("carries no such sentence", joined)
        self.assertIn("must never be added", joined)

    def test_the_gap_is_stated_as_values_and_as_records(self):
        records = [
            synthetic.good_record(condition_evidence=[]),
            synthetic.other_lab_record(condition_evidence=[]),
        ]
        joined = " ".join(section_lines(self.render(records)))
        self.assertIn("positive values with no evidence entry: 8, on 2 record(s)", joined)
        self.assertIn("W017", joined)

    def test_reach_is_rendered_with_its_denominator(self):
        records = [synthetic.good_record(), synthetic.other_lab_record(condition_evidence=[])]
        joined = " ".join(section_lines(self.render(records)))
        self.assertIn("studies with at least one evidenced value: 1 of 2", joined)

    def test_the_rule_table_is_rendered_as_the_audit_handle(self):
        joined = " ".join(section_lines(self.render([synthetic.good_record()])))
        self.assertIn("audit handle", joined)
        self.assertIn("synthetic", joined)

    def test_the_section_reports_no_score_or_share(self):
        """Counts and named gaps only; no figure of merit for a corpus."""
        lines = section_lines(self.render([synthetic.good_record()]))
        for line in lines:
            self.assertNotIn("%", line, line)
            self.assertNotIn("score", line.lower(), line)

    def test_the_specificity_section_still_stands_on_its_own(self):
        """Known and evidenced are different questions and stay apart."""
        blind = synthetic.good_record(
            sex=None, diet=None, age_group=None, challenge=None,
            age_value=None, age_unit=None, challenge_detail=None,
            condition_evidence=[],
        )
        report = run_qc([synthetic.good_record(), blind], vocab_loader=self.loader)
        text = render_text(report)
        self.assertIn(SPECIFICITY_TITLE, text)
        self.assertIn(TITLE, text)
        self.assertEqual(report.stats["condition_specificity"]["4"], 1)
        self.assertEqual(report.stats["condition_specificity"]["0"], 1)
        # The condition-blind record is not a gap: it has nothing to quote.
        self.assertEqual(
            report.stats["condition_evidence"]["gap"], {"records": 0, "values": 0}
        )

    def test_an_empty_dataset_renders_no_section_rather_than_a_table_of_zeros(self):
        text = render_text(run_qc([], vocab_loader=self.loader))
        self.assertNotIn(TITLE, text)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
