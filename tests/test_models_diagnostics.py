"""Dataset diagnostics: the molecular mirroring test, cells, and the headline.

These tests are about whether a *finding* can be produced, including a negative
one. Each fixture below is built so that the right answer is known in advance:
a corpus where loss and gain of function mirror, one where they contradict, one
where the cell pools two genes, and one where the readout never moves.
"""

import _path  # noqa: F401
import unittest

import synthetic
from syn4rch.models.diagnostics import (
    CELL_DEFINITIONS,
    cell_structure,
    dataset_diagnostics,
    modal_label,
    molecular_axis_mirroring,
    perturbation_identity_novelty,
)
from syn4rch.models.evaluate import run_benchmark
from syn4rch.models.report import headline, render_markdown
from syn4rch.models.task import build_dataset


def molecular(
    index,
    polarity="decrease",
    direction="increase",
    target="Gene1",
    outcome="body_weight",
    population="arc_pomc",
    study=None,
    lab=None,
):
    """One FABRICATED molecular record with the stated sign and result."""
    study = study if study is not None else "synthetic:m%d" % index
    lab = lab if lab is not None else "Lab%d" % index
    return synthetic.molecular_record(
        study_id=study,
        doi=study,
        last_author=lab,
        last_author_affiliation="Inst of %s" % lab,
        population_id=population,
        outcome_id=outcome,
        method=(
            "conditional_knockout"
            if polarity == "decrease"
            else "conditional_overexpression"
        ),
        perturbation_polarity=polarity,
        target_raw=target,
        target_id="syn:%s" % target.lower(),
        direction=direction,
        source_locator="Fig %d" % index,
    )


class TestModalLabel(unittest.TestCase):
    def test_ties_break_on_the_label_order_the_baselines_use(self):
        self.assertEqual(modal_label({"increase": 1, "decrease": 1}), "increase")
        self.assertEqual(modal_label({"decrease": 2, "increase": 1}), "decrease")
        self.assertIsNone(modal_label({}))
        self.assertIsNone(modal_label({"increase": 0}))


class TestMolecularMirroring(unittest.TestCase):
    def mirroring(self, records):
        return molecular_axis_mirroring(build_dataset(records))

    def test_a_gene_manipulated_both_ways_in_opposite_directions_mirrors(self):
        result = self.mirroring(
            [
                molecular(1, polarity="decrease", direction="increase"),
                molecular(2, polarity="increase", direction="decrease"),
            ]
        )
        gene = result["by_population_outcome_target"]
        self.assertEqual(gene["n_cells_with_both_polarities"], 1)
        self.assertEqual(gene["cell_verdicts"]["mirrored"], 1)
        self.assertEqual(gene["record_pairs"]["opposite_sign"], 1)
        self.assertEqual(gene["fraction_of_signed_pairs_mirrored"], 1.0)

    def test_a_gene_whose_two_directions_agree_is_counted_against_the_axis(self):
        result = self.mirroring(
            [
                molecular(1, polarity="decrease", direction="increase"),
                molecular(2, polarity="increase", direction="increase"),
            ]
        )
        gene = result["by_population_outcome_target"]
        self.assertEqual(gene["cell_verdicts"]["same"], 1)
        self.assertEqual(gene["cell_verdicts"]["mirrored"], 0)
        self.assertEqual(gene["record_pairs"]["same_sign"], 1)
        self.assertEqual(gene["fraction_of_signed_pairs_mirrored"], 0.0)

    def test_no_change_is_neither_evidence_for_nor_against(self):
        result = self.mirroring(
            [
                molecular(1, polarity="decrease", direction="no_change"),
                molecular(2, polarity="increase", direction="decrease"),
            ]
        )
        gene = result["by_population_outcome_target"]
        self.assertEqual(gene["cell_verdicts"]["no_change_involved"], 1)
        self.assertEqual(gene["record_pairs"]["no_change_involved"], 1)
        self.assertIsNone(gene["fraction_of_signed_pairs_mirrored"])

    def test_pooling_two_genes_is_reported_as_pooling_not_as_a_contradiction(self):
        # Two different genes, each manipulated one way. The gene-pooled cell
        # sees a both-polarity cell; the gene-specific view sees two cells with
        # one polarity each, i.e. no test at all.
        result = self.mirroring(
            [
                molecular(1, polarity="decrease", direction="increase", target="GeneA"),
                molecular(2, polarity="increase", direction="increase", target="GeneB"),
            ]
        )
        pooled = result["by_population_outcome"]
        gene = result["by_population_outcome_target"]
        self.assertEqual(pooled["n_cells_with_both_polarities"], 1)
        self.assertEqual(pooled["cells"][0]["n_distinct_targets"], 2)
        self.assertEqual(gene["n_cells_with_both_polarities"], 0)
        self.assertEqual(gene["n_cells"], 2)

    def test_single_study_and_multi_study_cells_are_counted_separately(self):
        records = [
            molecular(1, polarity="decrease", direction="increase",
                      study="synthetic:one", lab="LabOne"),
            molecular(2, polarity="increase", direction="decrease",
                      study="synthetic:one", lab="LabOne"),
        ]
        gene = self.mirroring(records)["by_population_outcome_target"]
        self.assertEqual(
            gene["cell_verdicts_by_provenance"]["single_study"]["mirrored"], 1
        )
        self.assertEqual(
            gene["cell_verdicts_by_provenance"]["multiple_studies"]["mirrored"], 0
        )
        self.assertEqual(
            gene["n_cells_with_both_polarities_and_two_or_more_labs"], 0
        )

    def test_two_labs_on_one_gene_are_counted_as_the_replication_they_are(self):
        records = [
            molecular(1, polarity="decrease", direction="increase", lab="LabA"),
            molecular(2, polarity="increase", direction="decrease", lab="LabB"),
        ]
        gene = self.mirroring(records)["by_population_outcome_target"]
        self.assertEqual(
            gene["n_cells_with_both_polarities_and_two_or_more_labs"], 1
        )

    def test_circuit_records_are_not_part_of_the_molecular_test(self):
        result = self.mirroring(
            [
                synthetic.good_record(source_locator="c1"),
                synthetic.good_record(valence="inhibition", source_locator="c2"),
            ]
        )
        self.assertEqual(result["n_molecular_examples"], 0)
        self.assertEqual(result["by_population_outcome"]["n_cells"], 0)

    def test_a_record_with_no_target_still_counts_in_the_pooled_view(self):
        result = self.mirroring(
            [
                molecular(1, polarity="decrease", direction="increase"),
                synthetic.molecular_record(
                    study_id="synthetic:nt", doi="synthetic:nt",
                    method="conditional_overexpression",
                    perturbation_polarity="increase",
                    target_raw=None, target_id=None,
                    population_id="arc_pomc", outcome_id="body_weight",
                    direction="decrease", source_locator="Fig NT",
                ),
            ]
        )
        self.assertEqual(result["n_examples_with_known_target"], 1)
        self.assertEqual(
            result["by_population_outcome"]["n_cells_with_both_polarities"], 1
        )
        # The two differ in the gene-specific view: one cell has a named gene,
        # the other the <null> token, so they are not the same cell.
        self.assertEqual(
            result["by_population_outcome_target"]["n_cells_with_both_polarities"], 0
        )


class TestCellStructure(unittest.TestCase):
    def test_every_declared_cell_definition_is_measured(self):
        structure = cell_structure(build_dataset([synthetic.good_record()]))
        self.assertEqual(
            [row["cell"] for row in structure["definitions"]],
            [name for name, _ in CELL_DEFINITIONS],
        )

    def test_a_finer_cell_fragments_and_its_purity_rises_mechanically(self):
        records = [
            synthetic.good_record(
                population_id="arc_agrp",
                outcome_id="food_intake",
                valence=("activation", "inhibition")[i % 2],
                perturbation_polarity=("increase", "decrease")[i % 2],
                direction=("increase", "decrease")[i % 2],
                study_id="synthetic:f%d" % i,
                doi="synthetic:f%d" % i,
                source_locator="Fig %d" % i,
            )
            for i in range(8)
        ]
        rows = {r["cell"]: r for r in cell_structure(build_dataset(records))["definitions"]}
        coarse = rows["population_outcome"]
        fine = rows["population_axis_outcome"]
        self.assertEqual(coarse["n_cells"], 1)
        self.assertEqual(fine["n_cells"], 2)
        self.assertLess(coarse["in_sample_purity"], fine["in_sample_purity"])
        self.assertEqual(fine["in_sample_purity"], 1.0)

    def test_purity_is_printed_next_to_the_singleton_share(self):
        records = [
            synthetic.good_record(
                population_id="syn:p%d" % i, outcome_id="syn:o%d" % i,
                source_locator="Fig %d" % i,
            )
            for i in range(5)
        ]
        row = cell_structure(build_dataset(records))["definitions"][0]
        self.assertEqual(row["in_sample_purity"], 1.0)
        self.assertEqual(row["fraction_singleton_cells"], 1.0)


class TestIdentityNovelty(unittest.TestCase):
    def test_folds_with_no_shared_identity_report_a_fraction_of_one(self):
        records = [
            molecular(i, target="Gene%d" % i, direction="increase")
            for i in range(4)
        ]
        dataset = build_dataset(records)
        folds = [((0, 1), (2, 3)), ((2, 3), (0, 1))]
        novelty = perturbation_identity_novelty(dataset, folds)
        self.assertEqual(novelty["fraction_unseen"], 1.0)
        self.assertEqual(novelty["n_test_rows"], 4)
        self.assertEqual(len(novelty["per_fold"]), 2)

    def test_a_shared_identity_is_not_counted_as_unseen(self):
        records = [
            molecular(i, target="Shared", direction="increase") for i in range(4)
        ]
        dataset = build_dataset(records)
        novelty = perturbation_identity_novelty(dataset, [((0, 1), (2, 3))])
        self.assertEqual(novelty["fraction_unseen"], 0.0)

    def test_a_molecular_record_with_no_gene_makes_it_a_declared_lower_bound(self):
        records = [
            synthetic.molecular_record(
                study_id="synthetic:u%d" % i, doi="synthetic:u%d" % i,
                target_raw=None, target_id=None,
                source_locator="Fig %d" % i,
            )
            for i in range(4)
        ]
        dataset = build_dataset(records)
        novelty = perturbation_identity_novelty(dataset, [((0, 1), (2, 3))])
        self.assertTrue(novelty["identity_is_coarse"])
        self.assertEqual(novelty["n_molecular_without_target"], 4)
        self.assertIn("lower bound", novelty["note"])

    def test_a_circuit_record_has_no_gene_to_name_so_is_not_coarse(self):
        # A chemogenetic manipulation has no molecular target by construction.
        # Calling its identity incomplete would invent a missing field.
        records = [
            synthetic.good_record(
                study_id="synthetic:c%d" % i, doi="synthetic:c%d" % i,
                source_locator="Fig %d" % i,
            )
            for i in range(4)
        ]
        dataset = build_dataset(records)
        novelty = perturbation_identity_novelty(dataset, [((0, 1), (2, 3))])
        self.assertFalse(novelty["identity_is_coarse"])
        self.assertEqual(novelty["n_molecular_examples"], 0)
        self.assertNotIn("lower bound", novelty["note"])


class TestHeadline(unittest.TestCase):
    def test_an_empty_dataset_produces_a_statement_that_nothing_was_measured(self):
        result = run_benchmark(build_dataset([]), n_bootstrap=10)
        text = headline(result)
        self.assertIn("Nothing was measured", text)
        self.assertNotIn("beats", text)

    def test_when_nothing_beats_the_floor_the_headline_says_so_in_plain_words(self):
        # Every row the same direction: no cell can out-score the global mode.
        records = [
            synthetic.good_record(
                study_id="synthetic:n%d" % i, doi="synthetic:n%d" % i,
                population_id=["arc_agrp", "arc_pomc"][i % 2],
                outcome_id="food_intake",
                direction="increase",
                source_locator="Fig %d" % i,
            )
            for i in range(12)
        ]
        result = run_benchmark(
            build_dataset(records), groupings=("study_id",), n_bootstrap=10
        )
        text = headline(result)
        self.assertIn(
            "No baseline beats guessing the global majority direction", text
        )
        self.assertIn("not predictable", text)

    def test_the_headline_names_its_grouping_and_refuses_to_overclaim(self):
        records = []
        for i in range(24):
            polarity = ("increase", "decrease")[i % 2]
            records.append(
                synthetic.good_record(
                    study_id="synthetic:h%d" % (i % 8),
                    doi="synthetic:h%d" % (i % 8),
                    population_id="arc_agrp",
                    outcome_id="food_intake",
                    valence="activation" if polarity == "increase" else "inhibition",
                    perturbation_polarity=polarity,
                    direction="increase" if polarity == "increase" else "decrease",
                    source_locator="Fig %d" % i,
                )
            )
        result = run_benchmark(
            build_dataset(records), groupings=("study_id",), n_bootstrap=50
        )
        text = headline(result)
        self.assertIn("study_id", text)
        self.assertIn("It is not \"direction is now predictable\"", text)
        self.assertNotIn("direction is now predictable.**", text)

    def test_the_report_carries_the_headline_and_both_new_sections(self):
        result = run_benchmark(
            build_dataset([synthetic.good_record(), synthetic.other_lab_record()]
                          + [synthetic.molecular_record(source_locator="m%d" % i)
                             for i in range(2)]),
            groupings=("study_id",),
            n_bootstrap=10,
        )
        text = render_markdown(result)
        self.assertIn("## Headline", text)
        self.assertIn("Before and after the schema v0.2.0 features", text)
        self.assertIn("Does the molecular axis carry a sign?", text)
        self.assertIn("The bar a model has to clear", text)
        self.assertIn("Insufficient evidence to rank these interventions", text)


class TestDatasetDiagnostics(unittest.TestCase):
    def test_both_diagnostics_are_present_and_name_what_they_measure(self):
        out = dataset_diagnostics(build_dataset([synthetic.molecular_record()]))
        self.assertIn("molecular_axis_mirroring", out)
        self.assertIn("cell_structure", out)
        self.assertIn(
            "what_would_falsify_this", out["molecular_axis_mirroring"]
        )

    def test_diagnostics_on_an_empty_dataset_are_empty_not_absent(self):
        out = dataset_diagnostics(build_dataset([]))
        self.assertEqual(out["molecular_axis_mirroring"]["n_molecular_examples"], 0)
        for row in out["cell_structure"]["definitions"]:
            self.assertEqual(row["n_cells"], 0)
            self.assertIsNone(row["in_sample_purity"])


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
