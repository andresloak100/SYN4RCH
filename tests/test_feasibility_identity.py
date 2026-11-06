"""Perturbation identity: what a held-out record shares with training.

SCHEMA.md 5b. Grouping keeps a study or a lab off both sides of a split; it
does not establish that anything on the test side was learnable from the train
side. These tests pin the three things that are easy to get wrong:

1. an identity whose target is missing is counted **separately**, never as
   learnable -- otherwise a null field would read as the refutation of the
   finding;
2. the grouping travels with every fraction;
3. the registry is ``syn4rch.models.task.GROUPINGS`` and nothing may be
   invented or dropped.
"""

import _path  # noqa: F401
import contextlib
import io
import json
import os
import shutil
import tempfile
import unittest

import synthetic
from syn4rch.feasibility.identity import (
    ClassCellReach,
    IDENTITY_FIELDS,
    cell_reach,
    STATUS_RESOLVED,
    STATUS_TARGET_UNRESOLVED,
    STATUS_UNRESOLVABLE,
    TARGET_CONSTITUTIVE_CLASSES,
    UnlabelledIdentityError,
    analyse,
    class_structure,
    identity_coverage,
    perturbation_identity,
    render_markdown,
    top_line,
)
from syn4rch.feasibility.identity_run import (
    DEFAULT_RECORDS,
    DEFAULT_REPORT,
    SCOPES,
)
from syn4rch.feasibility.identity_run import main as identity_main
from syn4rch.feasibility.splits import UnlabelledResultError
from syn4rch.models.task import GROUPINGS
from syn4rch.schema import write_jsonl


def _mol(target, study, author, affiliation="Institute, Somewhere", **overrides):
    """A fabricated molecular record: one gene, one study, one lab."""
    data = dict(
        target_raw=target,
        target_id="syn:%s" % target.lower(),
        study_id="synthetic:%s" % study,
        doi="synthetic:%s" % study,
        last_author=author,
        last_author_affiliation=affiliation,
    )
    data.update(overrides)
    return synthetic.molecular_record(**data)


def _circuit(study, author, affiliation="Institute, Somewhere", **overrides):
    """A fabricated circuit record: the same manipulation everyone performs."""
    data = dict(
        study_id="synthetic:%s" % study,
        doi="synthetic:%s" % study,
        last_author=author,
        last_author_affiliation=affiliation,
    )
    data.update(overrides)
    return synthetic.good_record(**data)


class TestWhatAnIdentityIs(unittest.TestCase):
    def test_it_is_the_perturbation_key_class_polarity_target(self):
        self.assertEqual(
            IDENTITY_FIELDS,
            (
                "perturbation_class",
                "perturbation_polarity",
                "perturbation_target",
            ),
        )
        identity = perturbation_identity(synthetic.molecular_record())
        self.assertEqual(
            identity.key, ("molecular", "decrease", "syn:synthgene1")
        )
        self.assertEqual(identity.status, STATUS_RESOLVED)

    def test_a_circuit_record_with_no_target_is_resolved_not_missing(self):
        # The manipulation acts on the population, and the population is a
        # separate component of the claim key.
        identity = perturbation_identity(synthetic.good_record())
        self.assertEqual(identity.perturbation_class, "circuit_activity")
        self.assertEqual(identity.key[2], "<null>")
        self.assertEqual(identity.status, STATUS_RESOLVED)
        self.assertTrue(identity.is_resolved)

    def test_a_molecular_record_with_no_target_is_underdetermined(self):
        identity = perturbation_identity(
            synthetic.molecular_record(target_raw=None, target_id=None)
        )
        self.assertEqual(identity.status, STATUS_TARGET_UNRESOLVED)
        self.assertFalse(identity.is_resolved)
        self.assertIn("molecular", TARGET_CONSTITUTIVE_CLASSES)

    def test_a_record_with_no_method_is_unresolvable(self):
        identity = perturbation_identity(synthetic.good_record(method=None))
        self.assertEqual(identity.status, STATUS_UNRESOLVABLE)

    def test_an_unresolved_target_prefers_the_raw_symbol_over_nothing(self):
        identity = perturbation_identity(
            synthetic.molecular_record(target_id=None)
        )
        self.assertEqual(identity.status, STATUS_RESOLVED)
        # The 'raw:' prefix keeps it from colliding with a resolved id.
        self.assertTrue(identity.key[2].startswith("raw:"))


class TestUnresolvedTargetsAreNeverCountedAsLearnable(unittest.TestCase):
    """The one serious misreading available here, refused structurally.

    With ``target`` null, every knockout in a corpus shares one identity, that
    identity occurs in many studies, and a naive count would report the
    molecular class as *learnable* -- the exact opposite of the finding,
    manufactured by a field that has not been filled in.
    """

    def setUp(self):
        self.records = [
            synthetic.molecular_record(
                target_raw=None,
                target_id=None,
                study_id="synthetic:m1",
                doi="synthetic:m1",
                last_author="One, A.",
            ),
            synthetic.molecular_record(
                target_raw=None,
                target_id=None,
                study_id="synthetic:m2",
                doi="synthetic:m2",
                last_author="Two, B.",
            ),
        ]

    def test_they_are_in_neither_bucket(self):
        coverage = identity_coverage(self.records, "lab")
        self.assertEqual(coverage.n_records, 2)
        self.assertEqual(coverage.n_records_target_unresolved, 2)
        self.assertEqual(coverage.n_records_resolved, 0)
        self.assertEqual(coverage.n_records_single_group, 0)
        self.assertEqual(coverage.n_records_multi_group, 0)
        self.assertEqual(coverage.n_records_placed, 0)

    def test_the_fraction_is_not_measurable_rather_than_zero(self):
        coverage = identity_coverage(self.records, "lab")
        # 0.0 would print as "none of them is unlearnable", which is the
        # opposite of what a missing field licenses.
        self.assertIsNone(coverage.fraction_single_group)
        self.assertIn("not measurable", coverage.stated)
        self.assertEqual(coverage.fraction_of_corpus_unmeasurable, 1.0)

    def test_resolving_the_target_is_what_reveals_the_structure(self):
        before = identity_coverage(self.records, "lab")
        self.assertIsNone(before.fraction_single_group)
        after = identity_coverage(
            [
                _mol("Gene1", "m1", "One, A."),
                _mol("Gene2", "m2", "Two, B.", "Other Institute, Elsewhere"),
            ],
            "lab",
        )
        self.assertEqual(after.n_records_placed, 2)
        self.assertEqual(after.fraction_single_group, 1.0)

    def test_a_pharmacological_record_cannot_reach_resolved_today(self):
        # Its identity is the drug and perturbation_key has no ligand slot.
        # Stated as a limitation of the key, not counted as one identity.
        record = synthetic.good_record(
            method="local_pharmacology", valence="activation"
        )
        identity = perturbation_identity(record)
        self.assertEqual(identity.perturbation_class, "pharmacological")
        self.assertEqual(identity.status, STATUS_TARGET_UNRESOLVED)
        self.assertIn("pharmacological", TARGET_CONSTITUTIVE_CLASSES)


class TestTheFindingItself(unittest.TestCase):
    """One lab, one gene, one paper -- against a shared circuit method."""

    def _corpus(self):
        return [
            # Three genes, three studies, three labs: a perfect one-to-one map.
            _mol("GeneA", "mol-a", "Alpha, A.", "Institute Alpha, A-town"),
            _mol("GeneB", "mol-b", "Beta, B.", "Institute Beta, B-town"),
            _mol("GeneC", "mol-c", "Gamma, C.", "Institute Gamma, C-town"),
            # Three labs performing the same chemogenetic activation.
            _circuit("cir-a", "Delta, D.", "Institute Delta, D-town"),
            _circuit("cir-b", "Epsilon, E.", "Institute Epsilon, E-town"),
            _circuit("cir-c", "Zeta, F.", "Institute Zeta, F-town"),
        ]

    def test_molecular_records_are_unlearnable_under_study_and_lab(self):
        for grouping in ("study_id", "lab"):
            coverage = identity_coverage(self._corpus(), grouping)
            molecular = coverage.class_coverage("molecular")
            self.assertEqual(molecular.n_records, 3, grouping)
            self.assertEqual(molecular.n_identities, 3, grouping)
            self.assertEqual(molecular.n_identities_single_group, 3, grouping)
            self.assertEqual(molecular.fraction_single_group, 1.0, grouping)

    def test_circuit_records_are_learnable_in_identity_terms(self):
        for grouping in ("study_id", "lab"):
            coverage = identity_coverage(self._corpus(), grouping)
            circuit = coverage.class_coverage("circuit_activity")
            self.assertEqual(circuit.n_records, 3, grouping)
            # One shared manipulation, seen in every group.
            self.assertEqual(circuit.n_identities, 1, grouping)
            self.assertEqual(circuit.n_identities_single_group, 0, grouping)
            self.assertEqual(circuit.fraction_single_group, 0.0, grouping)

    def test_that_contrast_is_the_point_and_is_not_averaged_away(self):
        coverage = identity_coverage(self._corpus(), "lab")
        # The corpus-level number is halfway between two classes that behave
        # completely differently, which is why the breakdown exists.
        self.assertEqual(coverage.fraction_single_group, 0.5)
        fractions = {
            entry.perturbation_class: entry.fraction_single_group
            for entry in coverage.by_class
        }
        self.assertEqual(fractions, {"circuit_activity": 0.0, "molecular": 1.0})

    def test_the_per_fold_unseen_fraction_is_one_in_the_molecular_folds(self):
        coverage = identity_coverage(self._corpus(), "study_id")
        per_fold = dict(coverage.per_fold_unseen_fraction)
        self.assertEqual(per_fold["synthetic:mol-a"], 1.0)
        self.assertEqual(per_fold["synthetic:cir-a"], 0.0)

    def test_a_second_lab_repeating_the_manipulation_is_the_only_remedy(self):
        corpus = self._corpus()
        corpus.append(
            _mol("GeneA", "mol-a-repeat", "Theta, G.", "Institute Theta, G-town")
        )
        coverage = identity_coverage(corpus, "lab")
        molecular = coverage.class_coverage("molecular")
        self.assertEqual(molecular.n_records, 4)
        # GeneA now occurs in two labs, so its two records become learnable.
        self.assertEqual(molecular.n_records_single_group, 2)
        self.assertEqual(molecular.n_records_multi_group, 2)
        self.assertEqual(molecular.fraction_single_group, 0.5)

    def test_another_paper_from_the_same_lab_changes_nothing(self):
        corpus = self._corpus()
        corpus.append(
            _mol(
                "GeneA",
                "mol-a-again",
                "Alpha, A.",
                "Institute Alpha, A-town",
                source_locator="Fig 7B (synthetic)",
            )
        )
        coverage = identity_coverage(corpus, "lab")
        molecular = coverage.class_coverage("molecular")
        self.assertEqual(molecular.n_records, 4)
        self.assertEqual(molecular.fraction_single_group, 1.0)


class TestEveryFractionNamesItsGrouping(unittest.TestCase):
    def test_an_undeclared_grouping_is_refused(self):
        with self.assertRaises(UnlabelledIdentityError):
            identity_coverage([synthetic.good_record()], "by_vibes")

    def test_it_is_refused_by_the_same_error_splits_uses(self):
        # SCHEMA.md 5's rule is one rule, whether the number is an accuracy or
        # an unseen-identity fraction.
        with self.assertRaises(UnlabelledResultError):
            identity_coverage([synthetic.good_record()], "by_vibes")

    def test_analyse_refuses_an_undeclared_grouping_too(self):
        with self.assertRaises(UnlabelledIdentityError):
            analyse([synthetic.good_record()], groupings=["study_id", "nope"])

    def test_the_stated_form_always_carries_the_grouping(self):
        coverage = identity_coverage([synthetic.good_record()], "study_id")
        self.assertIn("grouping=study_id", coverage.stated)

    def test_there_is_no_implicit_conversion_to_a_bare_fraction(self):
        coverage = identity_coverage([synthetic.good_record()], "study_id")
        with self.assertRaises(TypeError):
            float(coverage)

    def test_every_declared_grouping_is_measured_and_none_is_dropped(self):
        report = analyse([synthetic.good_record(), synthetic.molecular_record()])
        measured = sorted(c.grouping for c in report.coverages)
        self.assertEqual(measured, sorted(GROUPINGS))

    def test_the_leakage_prone_grouping_is_reported_not_hidden(self):
        report = analyse([synthetic.good_record()])
        record_level = report.coverage("record")
        self.assertIsNotNone(record_level)
        self.assertTrue(record_level.leakage_prone)


class TestUngroupableRecords(unittest.TestCase):
    def test_a_record_with_no_lab_is_excluded_and_counted(self):
        records = [
            _mol("GeneA", "mol-a", "Alpha, A."),
            synthetic.molecular_record(
                target_raw="GeneB",
                target_id="syn:geneb",
                study_id="synthetic:mol-b",
                doi="synthetic:mol-b",
                last_author=None,
                last_author_affiliation=None,
            ),
        ]
        coverage = identity_coverage(records, "lab")
        self.assertEqual(coverage.n_records_ungroupable, 1)
        self.assertEqual(coverage.n_records_placed, 1)
        # It is not pooled into an anonymous lab (DECISIONS.md D-009), and it
        # is not silently dropped either.
        self.assertEqual(coverage.n_groups, 1)


class TestPublicationStructure(unittest.TestCase):
    def test_one_study_per_lab_is_detected_as_a_bijection(self):
        structures = {
            s.perturbation_class: s
            for s in class_structure(
                [
                    _mol("GeneA", "mol-a", "Alpha, A.", "Institute Alpha, A"),
                    _mol("GeneB", "mol-b", "Beta, B.", "Institute Beta, B"),
                ]
            )
        }
        molecular = structures["molecular"]
        self.assertEqual(molecular.n_studies, 2)
        self.assertEqual(molecular.n_labs, 2)
        self.assertEqual(molecular.n_labs_with_multiple_studies, 0)
        self.assertTrue(molecular.study_lab_is_bijection)

    def test_two_studies_from_one_lab_break_the_bijection(self):
        structures = {
            s.perturbation_class: s
            for s in class_structure(
                [
                    _mol("GeneA", "mol-a", "Alpha, A.", "Institute Alpha, A"),
                    _mol("GeneB", "mol-b", "Alpha, A.", "Institute Alpha, A"),
                ]
            )
        }
        molecular = structures["molecular"]
        self.assertEqual(molecular.n_studies, 2)
        self.assertEqual(molecular.n_labs, 1)
        self.assertEqual(molecular.n_labs_with_multiple_studies, 1)
        self.assertFalse(molecular.study_lab_is_bijection)

    def test_an_unidentified_lab_prevents_the_bijection_claim(self):
        structures = {
            s.perturbation_class: s
            for s in class_structure(
                [
                    _mol("GeneA", "mol-a", "Alpha, A.", "Institute Alpha, A"),
                    synthetic.molecular_record(
                        target_raw="GeneB",
                        target_id="syn:geneb",
                        study_id="synthetic:mol-b",
                        doi="synthetic:mol-b",
                        last_author=None,
                        last_author_affiliation=None,
                    ),
                ]
            )
        }
        self.assertFalse(structures["molecular"].study_lab_is_bijection)

    def test_the_targeting_strategy_is_counted_and_never_parsed(self):
        structures = {
            s.perturbation_class: s
            for s in class_structure(
                [
                    _mol("GeneA", "mol-a", "Alpha, A.", "Institute Alpha, A"),
                    _mol("GeneB", "mol-b", "Beta, B.", "Institute Beta, B"),
                ]
            )
        }
        molecular = structures["molecular"]
        # Both fixtures share one strategy string, so it occurs in two studies.
        self.assertEqual(molecular.n_distinct_targeting_strategies, 1)
        self.assertEqual(
            molecular.n_targeting_strategies_in_multiple_studies, 1
        )
        # And it was never turned into a gene: the identity came from
        # target_id, not from the strategy text.
        for record in (_mol("GeneA", "mol-a", "Alpha, A."),):
            self.assertEqual(
                perturbation_identity(record).key[2], "syn:genea"
            )


class TestEmptyAndDegenerate(unittest.TestCase):
    def test_no_records_measures_nothing_and_says_so(self):
        report = analyse([])
        self.assertTrue(report.is_empty)
        for coverage in report.coverages:
            self.assertIsNone(coverage.fraction_single_group)
            self.assertEqual(coverage.n_records, 0)
        self.assertIn("No records", top_line(report))

    def test_an_empty_report_renders_without_a_zero_percent_headline(self):
        text = render_markdown(analyse([]))
        self.assertIn("Nothing to measure", text)
        self.assertNotIn("0.0% of placed records", text)

    def test_a_rendered_report_states_the_target_resolution_state(self):
        unresolved = render_markdown(
            analyse(
                [synthetic.molecular_record(target_raw=None, target_id=None)]
            )
        )
        self.assertIn("targets not yet resolved", unresolved.lower())
        resolved = render_markdown(analyse([synthetic.molecular_record()]))
        self.assertIn("targets fully resolved", resolved.lower())

    def test_the_report_is_deterministic_and_serialisable(self):
        records = [
            _mol("GeneA", "mol-a", "Alpha, A.", "Institute Alpha, A"),
            _circuit("cir-a", "Delta, D.", "Institute Delta, D"),
        ]
        at = "2001-01-01T00:00:00Z"
        first = analyse(records)
        second = analyse(list(reversed(records)))
        first.generated_at = second.generated_at = at
        self.assertEqual(first.to_dict(), second.to_dict())


class TestCli(unittest.TestCase):
    """The report must be exactly reproducible by the CLI."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="syn4rch-identity-")
        self.out = os.path.join(self.tmp, "perturbation_identity.md")
        self.records = os.path.join(self.tmp, "records.jsonl")
        write_jsonl(
            self.records,
            [
                _mol("GeneA", "mol-a", "Alpha, A.", "Institute Alpha, A"),
                _mol("GeneB", "mol-b", "Beta, B.", "Institute Beta, B"),
                _circuit("cir-a", "Delta, D.", "Institute Delta, D"),
                _circuit("cir-b", "Epsilon, E.", "Institute Epsilon, E"),
            ],
        )

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, *argv):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = identity_main(list(argv))
        return code, buffer.getvalue()

    def test_a_missing_input_is_not_an_error(self):
        code, _ = self._run(
            "--records", os.path.join(self.tmp, "nope.jsonl"), "--out", self.out
        )
        self.assertEqual(code, 0)
        with open(self.out) as handle:
            self.assertIn("No records", handle.read())

    def test_an_unlearnable_class_exits_zero(self):
        # A whole class being unpredictable is the success case, not a failure.
        code, stdout = self._run("--records", self.records, "--out", self.out)
        self.assertEqual(code, 0)
        # The headline quotes a grouping SCHEMA.md 5b names, not whichever
        # grouping happens to give the biggest number.
        self.assertIn("`lab`", stdout)
        with open(self.out) as handle:
            text = handle.read()
        self.assertIn("100.0%", text)
        self.assertIn("cannot support cross-study prediction", text)

    def test_the_contrast_between_the_classes_is_in_the_written_report(self):
        self._run("--records", self.records, "--out", self.out)
        with open(self.out) as handle:
            text = handle.read()
        self.assertIn("circuit_activity", text)
        self.assertIn("molecular", text)
        self.assertIn("one laboratory, one gene, one paper", text)
        self.assertIn("a second laboratory repeating the same manipulation", text)

    def test_both_scopes_appear_whichever_one_leads(self):
        for scope in SCOPES:
            self._run(
                "--records", self.records, "--out", self.out, "--scope", scope
            )
            with open(self.out) as handle:
                text = handle.read()
            self.assertIn("The other scope, for comparison", text, scope)

    def test_json_output_carries_both_scopes_and_every_grouping(self):
        json_path = os.path.join(self.tmp, "identity.json")
        code, _ = self._run(
            "--records", self.records, "--out", self.out, "--json", json_path
        )
        self.assertEqual(code, 0)
        with open(json_path) as handle:
            data = json.load(handle)
        self.assertEqual(data["primary_scope"], "all")
        self.assertEqual(data["secondary_scope"], "usable")
        measured = sorted(c["grouping"] for c in data["primary"]["coverages"])
        self.assertEqual(measured, sorted(GROUPINGS))
        for coverage in data["primary"]["coverages"]:
            # No fraction is ever serialised without its grouping.
            self.assertIn("grouping", coverage)
            self.assertIn("fraction_single_group", coverage)

    def test_stdout_mode_writes_no_file(self):
        code, text = self._run("--records", self.records, "--stdout")
        self.assertEqual(code, 0)
        self.assertFalse(os.path.exists(self.out))
        self.assertIn("# SYN4RCH perturbation identity report", text)

    def test_quiet_suppresses_stdout_but_still_writes_the_report(self):
        code, text = self._run(
            "--records", self.records, "--out", self.out, "--quiet"
        )
        self.assertEqual(code, 0)
        self.assertEqual(text, "")
        self.assertTrue(os.path.exists(self.out))

    def test_it_is_reproducible(self):
        self._run("--records", self.records, "--out", self.out, "--quiet")
        with open(self.out) as handle:
            first = handle.read()
        second_path = os.path.join(self.tmp, "again.md")
        self._run("--records", self.records, "--out", second_path, "--quiet")
        with open(second_path) as handle:
            second = handle.read()

        def strip_timestamp(text):
            return "\n".join(
                line for line in text.splitlines() if "generated (UTC)" not in line
            )

        self.assertEqual(strip_timestamp(first), strip_timestamp(second))

    def test_it_does_not_modify_its_input(self):
        with open(self.records, "rb") as handle:
            before = handle.read()
        self._run("--records", self.records, "--out", self.out, "--quiet")
        with open(self.records, "rb") as handle:
            self.assertEqual(handle.read(), before)

    def test_the_default_paths_are_the_documented_ones(self):
        self.assertEqual(DEFAULT_RECORDS, os.path.join("data", "curated", "records.jsonl"))
        self.assertEqual(
            DEFAULT_REPORT, os.path.join("reports", "perturbation_identity.md")
        )

    def test_an_unwritable_output_path_is_a_usage_error(self):
        code, _ = self._run(
            "--records",
            self.records,
            "--out",
            os.path.join(self.records, "nested", "out.md"),
        )
        self.assertEqual(code, 2)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()


class CellReachTest(unittest.TestCase):
    """Cell-level reach: the level at which the class contrast is real."""

    def test_two_labs_sharing_a_cell_make_it_multi_group(self):
        rows = [
            synthetic.good_record(population_raw="p1", population_id="syn:p1",
                                  outcome_id="food_intake", outcome_raw="o"),
            synthetic.other_lab_record(population_raw="p1", population_id="syn:p1",
                                       outcome_id="food_intake", outcome_raw="o"),
        ]
        got = {r.perturbation_class: r for r in cell_reach(rows, "lab")}
        circuit = got["circuit_activity"]
        self.assertEqual(circuit.n_cells, 1)
        self.assertEqual(circuit.n_cells_multi_group, 1)
        self.assertEqual(circuit.fraction_cells_multi_group, 1.0)

    def test_one_lab_per_cell_makes_none_multi_group(self):
        rows = [
            synthetic.good_record(population_raw="p1", population_id="syn:p1",
                                  outcome_id="food_intake", outcome_raw="o"),
            synthetic.other_lab_record(population_raw="p2", population_id="syn:p2",
                                       outcome_id="body_weight", outcome_raw="o2"),
        ]
        got = {r.perturbation_class: r for r in cell_reach(rows, "lab")}
        circuit = got["circuit_activity"]
        self.assertEqual(circuit.n_cells, 2)
        self.assertEqual(circuit.n_cells_multi_group, 0)

    def test_records_in_multi_group_cells_are_counted(self):
        """The field was declared, defaulted to 0, and never populated (D-056).

        `to_dict` serialised it, so a reader of the JSON was told that zero
        records lived in the 26 multi-group circuit cells of the real corpus.
        A number that is always 0 is worse than an absent one.
        """
        rows = [
            synthetic.good_record(population_raw="p1", population_id="syn:p1",
                                  outcome_id="food_intake", outcome_raw="o"),
            synthetic.other_lab_record(population_raw="p1", population_id="syn:p1",
                                       outcome_id="food_intake", outcome_raw="o"),
            # A second cell, one lab only: its rows must not be counted.
            synthetic.other_lab_record(population_raw="p2", population_id="syn:p2",
                                       outcome_id="body_weight", outcome_raw="o2",
                                       source_locator="Fig 9Z (synthetic)"),
        ]
        circuit = {r.perturbation_class: r for r in cell_reach(rows, "lab")}[
            "circuit_activity"
        ]
        self.assertEqual(circuit.n_records, 3)
        self.assertEqual(circuit.n_cells, 2)
        self.assertEqual(circuit.n_cells_multi_group, 1)
        self.assertEqual(circuit.n_records_in_multi_group_cells, 2)
        self.assertEqual(
            circuit.to_dict()["n_records_in_multi_group_cells"], 2
        )

    def test_no_multi_group_cell_means_no_records_in_one(self):
        rows = [
            synthetic.good_record(population_raw="p1", population_id="syn:p1",
                                  outcome_id="food_intake", outcome_raw="o"),
            synthetic.other_lab_record(population_raw="p2", population_id="syn:p2",
                                       outcome_id="body_weight", outcome_raw="o2"),
        ]
        circuit = {r.perturbation_class: r for r in cell_reach(rows, "lab")}[
            "circuit_activity"
        ]
        self.assertEqual(circuit.n_records_in_multi_group_cells, 0)

    def test_the_share_is_none_not_zero_when_there_are_no_cells(self):
        """"No cells" and "no cell reaches two groups" are opposite findings."""
        row = ClassCellReach(perturbation_class="molecular", grouping="lab")
        self.assertIsNone(row.fraction_cells_multi_group)

    def test_an_undeclared_grouping_raises_rather_than_defaulting(self):
        with self.assertRaises(UnlabelledIdentityError):
            cell_reach([synthetic.good_record()], "by_vibes")

    def test_the_report_carries_cell_reach_for_every_grouping(self):
        report = analyse([synthetic.good_record()])
        groupings = {r.grouping for r in report.cell_reaches}
        self.assertEqual(groupings, {c.grouping for c in report.coverages})
