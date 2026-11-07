"""The similarity feature, and the pre-registered choices it must not drift from.

D-059 requires the gene set and `k` to be fixed in code before the first
evaluation, and D-066 fixes them. Several tests here exist only to make a later
change to those constants fail loudly, because changing one after a delta is
visible is the move the pre-registration exists to prevent.
"""

import _path  # noqa: F401
import unittest

from syn4rch.atlas.similarity import (
    BOOTSTRAP_REPLICATES,
    BOOTSTRAP_SEED,
    CENTRE_GENES,
    STABILITY_MIN,
    bootstrap_gene_draws,
    centre_genes,
    is_rankable,
    nearest_stability,
    EXCLUDED_PREFIXES,
    gene_variances,
    select_from_variance,
    K_EXPLORATORY,
    K_NEIGHBOURS,
    LIVE_GENE_FLOOR,
    N_GENES,
    PRIMARY_TIERS,
    SECONDARY_TIERS,
    NeighbourChoice,
    PopulationVector,
    choose_neighbour,
    cosine,
    jaccard,
    select_genes,
    similarity_matrix,
)


def vec(pid, values, clusters=("c1",), tier="A", cells=100, leaves=None):
    return PopulationVector(population_id=pid, tier=tier,
                            values=tuple(values), clusters=tuple(clusters),
                            n_cells=cells,
                            leaf_clusters=tuple(leaves or ()))


class PreRegisteredConstantsTest(unittest.TestCase):
    """If one of these fails, check D-066 before changing the test."""

    def test_the_gene_count_is_the_declared_one(self):
        self.assertEqual(N_GENES, 2000)

    def test_k_is_one_for_the_primary(self):
        self.assertEqual(K_NEIGHBOURS, 1)

    def test_the_exploratory_k_is_separate_and_cannot_be_the_primary(self):
        self.assertEqual(K_EXPLORATORY, 3)
        self.assertNotEqual(K_EXPLORATORY, K_NEIGHBOURS)

    def test_only_tier_A_is_admitted_to_the_primary(self):
        self.assertEqual(PRIMARY_TIERS, ("A",))
        self.assertEqual(SECONDARY_TIERS, ("A", "B", "C"))

    def test_technical_covariates_are_excluded_from_selection(self):
        self.assertEqual(EXCLUDED_PREFIXES, ("mt-", "Rps", "Rpl"))

    def test_genes_are_centred_per_D067(self):
        self.assertIs(CENTRE_GENES, True)

    def test_the_tie_guard_thresholds_are_the_declared_ones(self):
        self.assertEqual(BOOTSTRAP_REPLICATES, 200)
        self.assertEqual(STABILITY_MIN, 0.95)


class CentringTest(unittest.TestCase):
    def test_each_gene_loses_its_mean_across_populations(self):
        out = centre_genes([[1.0, 10.0], [3.0, 20.0]])
        self.assertEqual(out, ((-1.0, -5.0), (1.0, 5.0)))

    def test_centring_breaks_the_near_degeneracy_it_was_added_for(self):
        """D-067: two neurons sharing a large baseline are ~1.0 apart under a
        raw cosine, which is what made a 0.0051 top-1 margin decide every row.
        """
        a = [10.0, 10.0, 1.0]
        b = [10.0, 10.0, 0.0]
        self.assertGreater(cosine(a, b), 0.98)
        ca, cb = centre_genes([a, b])
        self.assertLess(cosine(ca, cb), 0.0)

    def test_an_empty_matrix_centres_to_nothing(self):
        self.assertEqual(centre_genes([]), ())

    def test_a_single_population_centres_to_zero(self):
        """One population is its own gene mean, so nothing survives -- and the
        cosine machinery must return 0.0 rather than divide by zero."""
        out = centre_genes([[1.0, 2.0]])
        self.assertEqual(out, ((0.0, 0.0),))
        self.assertEqual(cosine(out[0], out[0]), 0.0)


class BootstrapTest(unittest.TestCase):
    def test_the_draws_are_deterministic_for_a_fixed_seed(self):
        """A tie guard whose verdict moved between runs would decide which rows
        are reported on nothing at all."""
        self.assertEqual(bootstrap_gene_draws(50, replicates=4, seed=1),
                         bootstrap_gene_draws(50, replicates=4, seed=1))

    def test_a_different_seed_gives_different_draws(self):
        self.assertNotEqual(bootstrap_gene_draws(50, replicates=4, seed=1),
                            bootstrap_gene_draws(50, replicates=4, seed=2))

    def test_each_draw_is_the_full_width_with_replacement(self):
        draws = bootstrap_gene_draws(10, replicates=3, seed=BOOTSTRAP_SEED)
        self.assertEqual(len(draws), 3)
        for d in draws:
            self.assertEqual(len(d), 10)
            self.assertTrue(all(0 <= i < 10 for i in d))

    def test_the_default_seed_is_the_declared_one(self):
        self.assertEqual(BOOTSTRAP_SEED, 20670)


class StabilityTest(unittest.TestCase):
    def test_an_unshakeable_nearest_neighbour_scores_one(self):
        def nearest(_genes):
            return {"a": "b", "b": "a"}
        out = nearest_stability(["a", "b"], nearest, 10, replicates=20)
        self.assertEqual(out, {"a": 1.0, "b": 1.0})

    def test_a_neighbour_that_never_reappears_scores_zero(self):
        calls = {"n": 0}

        def nearest(_genes):
            calls["n"] += 1
            return {"a": "b"} if calls["n"] == 1 else {"a": "c"}
        out = nearest_stability(["a"], nearest, 10, replicates=20)
        self.assertEqual(out["a"], 0.0)

    def test_the_full_gene_set_sets_the_reference_not_a_replicate(self):
        seen = []

        def nearest(genes):
            seen.append(genes)
            return {"a": "b"}
        nearest_stability(["a"], nearest, 5, replicates=2)
        self.assertEqual(seen[0], (0, 1, 2, 3, 4))

    def test_a_row_with_no_neighbour_at_all_does_not_count_as_agreement(self):
        """None == None must not read as a stable ranking."""
        def nearest(_genes):
            return {"a": None}
        self.assertEqual(nearest_stability(["a"], nearest, 5, replicates=10)["a"],
                         0.0)


class RankabilityTest(unittest.TestCase):
    def test_a_stable_row_is_rankable(self):
        self.assertTrue(is_rankable(1.0))
        self.assertTrue(is_rankable(STABILITY_MIN))

    def test_an_unstable_row_is_not_and_that_is_the_honest_answer(self):
        """The mechanism for "insufficient evidence to rank these reliably"."""
        self.assertFalse(is_rankable(STABILITY_MIN - 0.01))
        self.assertFalse(is_rankable(0.0))


class GeneSelectionTest(unittest.TestCase):
    def test_genes_are_chosen_by_variance_across_clusters(self):
        means = [[0.0, 5.0, 1.0], [10.0, 5.0, 1.2]]
        names = ["varies", "flat", "almost_flat"]
        self.assertEqual(select_genes(means, names, n_genes=1), (0,))

    def test_a_gene_high_everywhere_is_not_selected(self):
        """`Malat1` tops nearly every cluster in this atlas. High is not
        informative; different is."""
        means = [[9.0, 0.0, 0.1], [9.0, 4.0, 0.1]]
        names = ["high_and_flat", "varies", "low_and_flat"]
        self.assertEqual(select_genes(means, names, n_genes=1), (1,))

    def test_mitochondrial_and_ribosomal_genes_cannot_be_selected(self):
        means = [[0.0, 0.0], [99.0, 1.0]]
        names = ["mt-Cytb", "real"]
        self.assertEqual(select_genes(means, names, n_genes=2), (1,))

    def test_the_result_is_sorted_and_deterministic(self):
        means = [[0.0, 0.0, 0.0], [3.0, 1.0, 2.0]]
        names = ["a", "b", "c"]
        picked = select_genes(means, names, n_genes=2)
        self.assertEqual(picked, tuple(sorted(picked)))
        self.assertEqual(picked, select_genes(means, names, n_genes=2))

    def test_ties_break_on_index_so_two_runs_agree(self):
        means = [[0.0, 0.0], [1.0, 1.0]]
        self.assertEqual(select_genes(means, ["a", "b"], n_genes=1), (0,))

    def test_no_clusters_selects_nothing(self):
        self.assertEqual(select_genes([], ["a"], n_genes=5), ())

    def test_asking_for_more_genes_than_exist_returns_all_eligible(self):
        means = [[0.0, 1.0], [2.0, 3.0]]
        self.assertEqual(select_genes(means, ["a", "b"], n_genes=99), (0, 1))

    def test_the_policy_is_separate_from_the_arithmetic_and_they_agree(self):
        """The real run computes variance with numpy over 465 x 51,216 values.

        Only the arithmetic may differ; the policy applied to a variance vector
        must be the single declared one, so a fast path cannot quietly change
        which genes are chosen.
        """
        means = [[0.0, 5.0, 1.0, 9.0], [10.0, 5.0, 1.2, 9.0]]
        names = ["varies", "flat", "almost_flat", "mt-Cytb"]
        variances = gene_variances(means, len(names))
        self.assertEqual(
            select_genes(means, names, n_genes=2),
            select_from_variance(variances, names, n_genes=2))

    def test_variance_is_the_population_form_not_the_sample_form(self):
        """Divisor n, not n-1. Fixed so a different convention in a fast path
        shows up as a test failure rather than a slightly different gene set."""
        self.assertAlmostEqual(gene_variances([[0.0], [2.0]], 1)[0], 1.0)


class CosineTest(unittest.TestCase):
    def test_identical_vectors_are_one(self):
        self.assertAlmostEqual(cosine([1.0, 2.0], [1.0, 2.0]), 1.0)

    def test_scale_does_not_matter(self):
        self.assertAlmostEqual(cosine([1.0, 2.0], [10.0, 20.0]), 1.0)

    def test_orthogonal_vectors_are_zero(self):
        self.assertAlmostEqual(cosine([1.0, 0.0], [0.0, 1.0]), 0.0)

    def test_an_empty_vector_is_zero_rather_than_a_crash(self):
        """A population whose every selected gene is silent must not take the
        run down; it must show up as a zero-similarity row in the log."""
        self.assertEqual(cosine([0.0, 0.0], [1.0, 1.0]), 0.0)


class JaccardTest(unittest.TestCase):
    def test_disjoint_cluster_sets_are_zero(self):
        self.assertEqual(jaccard(["a"], ["b"]), 0.0)

    def test_identical_sets_are_one(self):
        self.assertEqual(jaccard(["a", "b"], ["b", "a"]), 1.0)

    def test_a_subset_is_the_ratio_of_sizes(self):
        self.assertAlmostEqual(jaccard(["a"], ["a", "b"]), 0.5)

    def test_two_empty_sets_are_zero_not_one(self):
        self.assertEqual(jaccard([], []), 0.0)


class VectorDiagnosticsTest(unittest.TestCase):
    def test_live_genes_counts_only_what_clears_the_floor(self):
        v = vec("p", [0.0, LIVE_GENE_FLOOR, LIVE_GENE_FLOOR + 0.01, 5.0])
        self.assertEqual(v.live_genes, 2)

    def test_a_population_with_no_live_genes_reports_zero(self):
        """The case this diagnostic exists for: a miss on such a population is
        not evidence about similarity, and without the count it looks the same
        as any other miss."""
        self.assertEqual(vec("p", [0.0, 0.0]).live_genes, 0)

    def test_the_norm_is_reported_alongside(self):
        self.assertAlmostEqual(vec("p", [3.0, 4.0]).norm, 5.0)


class SimilarityMatrixTest(unittest.TestCase):
    def test_the_matrix_is_symmetric(self):
        sims = similarity_matrix([vec("a", [1.0, 0.0]), vec("b", [0.0, 1.0])])
        self.assertEqual(sims[("a", "b")], sims[("b", "a")])

    def test_a_population_is_not_compared_with_itself(self):
        sims = similarity_matrix([vec("a", [1.0]), vec("b", [1.0])])
        self.assertNotIn(("a", "a"), sims)


class NeighbourChoiceTest(unittest.TestCase):
    def setUp(self):
        self.held = vec("held", [1.0, 0.0, 0.0], clusters=["shared", "own"])

    def test_the_nearest_non_overlapping_neighbour_wins(self):
        near = vec("near", [0.9, 0.1, 0.0], clusters=["elsewhere"])
        far = vec("far", [0.0, 0.0, 1.0], clusters=["other"])
        sims = {("held", "near"): 0.95, ("held", "far"): 0.05}
        c = choose_neighbour(self.held, [far, near], sims)
        self.assertEqual(c.chosen, "near")
        self.assertAlmostEqual(c.similarity, 0.95)
        self.assertFalse(c.fell_back)

    def test_an_overlapping_neighbour_is_skipped_for_the_next_nearest(self):
        """D-066. The overlapping pair's similarity is real but correlated by
        construction, so it must not supply a prediction."""
        overlapping = vec("overlap", [1.0, 0.0, 0.0], clusters=["shared"])
        clean = vec("clean", [0.5, 0.5, 0.0], clusters=["elsewhere"])
        sims = {("held", "overlap"): 0.99, ("held", "clean"): 0.40}
        c = choose_neighbour(self.held, [overlapping, clean], sims)
        self.assertEqual(c.chosen, "clean")
        self.assertEqual(c.skipped_for_overlap, ("overlap",))
        self.assertEqual(c.jaccard_with_chosen, 0.0)

    def test_when_every_neighbour_overlaps_the_row_falls_back(self):
        """Counted unpredicted, never counted wrong."""
        a = vec("a", [1.0, 0.0, 0.0], clusters=["shared"])
        b = vec("b", [1.0, 0.0, 0.0], clusters=["own"])
        c = choose_neighbour(self.held, [a, b], {("held", "a"): 0.9,
                                                 ("held", "b"): 0.8})
        self.assertTrue(c.fell_back)
        self.assertIsNone(c.chosen)
        self.assertEqual(sorted(c.skipped_for_overlap), ["a", "b"])

    def test_no_candidates_at_all_falls_back_and_records_zero_available(self):
        c = choose_neighbour(self.held, [], {})
        self.assertTrue(c.fell_back)
        self.assertEqual(c.n_neighbours_available, 0)

    def test_disqualification_can_be_turned_off_only_explicitly(self):
        """The secondary reported alongside the primary needs the un-gated
        figure; the default must be the pre-registered behaviour."""
        overlapping = vec("overlap", [1.0, 0.0, 0.0], clusters=["shared"])
        sims = {("held", "overlap"): 0.99}
        gated = choose_neighbour(self.held, [overlapping], sims)
        ungated = choose_neighbour(self.held, [overlapping], sims,
                                   disqualify_overlap=False)
        self.assertTrue(gated.fell_back)
        self.assertEqual(ungated.chosen, "overlap")
        self.assertAlmostEqual(ungated.jaccard_with_chosen, 0.5)

    def test_disqualification_sees_overlap_across_levels(self):
        """D-066 disqualifies an overlapping neighbour, and before this fix it
        could only see overlap that shared a cluster NAME -- so two populations
        mapped at different levels sharing most of their cells passed the test.
        """
        held = vec("held", [1.0, 0.0], clusters=["C66-46"], leaves=["L1", "L2"])
        same_cells = vec("other", [0.9, 0.1], clusters=["C185-115"],
                         leaves=["L1", "L3"])
        sims = {("held", "other"): 0.99}
        c = choose_neighbour(held, [same_cells], sims)
        self.assertTrue(c.fell_back, "overlap across levels was not detected")
        self.assertEqual(c.skipped_for_overlap, ("other",))

    def test_overlap_key_falls_back_to_names_without_leaves(self):
        held = vec("held", [1.0], clusters=["shared"])
        self.assertEqual(held.overlap_key, ("shared",))
        withleaves = vec("h2", [1.0], clusters=["shared"], leaves=["L1"])
        self.assertEqual(withleaves.overlap_key, ("L1",))

    def test_every_diagnostic_d059_and_d066_require_is_populated(self):
        near = vec("near", [0.9, 0.1, 0.0], clusters=["elsewhere"])
        c = choose_neighbour(self.held, [near], {("held", "near"): 0.95})
        self.assertEqual(c.n_neighbours_available, 1)       # D-059 (a)
        self.assertEqual(c.chosen, "near")                  # D-059 (b)
        self.assertIsNotNone(c.similarity)                  # D-059 (b)
        self.assertIsNotNone(c.jaccard_with_chosen)         # D-066
        self.assertGreater(c.live_genes_self, 0)            # D-066 amendment
        self.assertGreater(c.live_genes_chosen, 0)

    def test_the_full_ranking_is_logged_not_just_the_winner(self):
        a = vec("a", [1.0, 0.0, 0.0], clusters=["x"])
        b = vec("b", [0.0, 1.0, 0.0], clusters=["y"])
        c = choose_neighbour(self.held, [a, b], {("held", "a"): 0.9,
                                                 ("held", "b"): 0.2})
        self.assertEqual([p for p, _s in c.ranked], ["a", "b"])

    def test_ties_in_similarity_break_on_population_id(self):
        a = vec("aaa", [1.0, 0.0, 0.0], clusters=["x"])
        b = vec("bbb", [1.0, 0.0, 0.0], clusters=["y"])
        sims = {("held", "aaa"): 0.5, ("held", "bbb"): 0.5}
        self.assertEqual(choose_neighbour(self.held, [b, a], sims).chosen, "aaa")

    def test_a_missing_similarity_is_treated_as_zero_not_an_error(self):
        a = vec("a", [1.0], clusters=["x"])
        self.assertEqual(choose_neighbour(self.held, [a], {}).chosen, "a")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
