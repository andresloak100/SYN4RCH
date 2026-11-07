"""The population -> cluster rule, and the four defects it was written against.

Every test here names a mistake that was actually made while this mapping was
first built by hand, and which the rule now cannot repeat.
"""

import _path  # noqa: F401
import unittest

from syn4rch.atlas.population_clusters import (
    GOOD_RANK,
    MARKER_FLOOR,
    MIN_CONTRAST,
    MIN_REGION_SHARE,
    OUTCOMES,
    TIERS,
    Assignment,
    Candidate,
    assign_population,
    nesting,
)


def cand(level, cluster, mean, *, region=100.0, rank=1, named=True):
    return Candidate(level=level, cluster=cluster, marker_mean=mean,
                     region_share=region, marker_rank_in_cluster=rank,
                     name_carries_marker=named)


class OutcomeVocabularyTest(unittest.TestCase):
    def test_every_outcome_the_rule_can_return_is_declared(self):
        seen = set()
        seen.add(assign_population("p", None, []).outcome)
        seen.add(assign_population("p", "Gene", []).outcome)
        seen.add(assign_population(
            "p", "Gene", [cand("C66", "a", 9.0, region=10.0)]).outcome)
        seen.add(assign_population(
            "p", "Gene", [cand("C66", "a", 0.01)]).outcome)
        seen.add(assign_population(
            "p", "Gene", [cand("C66", "a", 3.0), cand("C66", "b", 0.01),
                          cand("C66", "c", 0.01)]).outcome)
        self.assertEqual(seen, set(OUTCOMES))

    def test_a_mapped_assignment_is_not_refused_and_a_refusal_carries_a_reason(self):
        ok = assign_population("p", "Gene",
                               [cand("C66", "a", 3.0), cand("C66", "b", 0.01),
                                cand("C66", "c", 0.01)])
        self.assertFalse(ok.refused)
        self.assertIsNone(ok.reason)
        bad = assign_population("p", "Gene", [cand("C66", "a", 0.01)])
        self.assertTrue(bad.refused)
        self.assertTrue(bad.reason)


class FloorTest(unittest.TestCase):
    def test_a_receptor_that_never_clears_the_floor_is_refused(self):
        """Lepr peaks at 0.074 and Mc4r at 0.055 across their whole region.

        Both are real populations with real Cre lines; neither has enough
        transcript to mark a cluster. Refusing is the honest outcome, and the
        reason has to say where the best attempt landed.
        """
        a = assign_population("dmh_lepr", "Lepr", [
            cand("C185", "x", 0.074), cand("C465", "y", 0.060)])
        self.assertEqual(a.outcome, "refused_below_floor")
        self.assertIn("0.074", a.reason)
        self.assertIn("C185", a.reason)

    def test_the_floor_sits_between_the_two_observed_populations_of_markers(self):
        self.assertLess(0.391, MARKER_FLOOR)   # Nr5a1, the best TF
        self.assertLess(MARKER_FLOOR, 1.793)   # Kiss1, the weakest clean one


class LevelChoiceTest(unittest.TestCase):
    def test_the_level_is_chosen_by_contrast_not_by_fineness(self):
        """`pvh_crh` is the case: no PVH-majority Crh cluster exists at C185.

        A finer level always wins on raw magnitude because it averages fewer
        cells, so choosing "the finest that passes" would fragment every
        population. Contrast is what identifies the resolution at which the
        marker actually separates a subset of the region.
        """
        a = assign_population("p", "Gene", [
            # coarse: the marker is high but so is everything around it
            cand("C66", "coarse", 2.0),
            cand("C66", "coarse_other_1", 1.9),
            cand("C66", "coarse_other_2", 1.8),
            cand("C66", "coarse_other_3", 1.7),
            # fine: lower in absolute terms, but it stands alone
            cand("C465", "fine", 1.0),
            cand("C465", "fine_other_1", 0.02),
            cand("C465", "fine_other_2", 0.02),
            cand("C465", "fine_other_3", 0.02),
        ])
        self.assertEqual(a.level, "C465")
        self.assertEqual(a.clusters, ("fine",))

    def test_a_level_whose_best_is_below_the_floor_cannot_win_on_contrast(self):
        a = assign_population("p", "Gene", [
            cand("C66", "passes", 3.0), cand("C66", "passes_other", 0.2),
            cand("C66", "passes_other_2", 0.2),
            # enormous contrast, but 0.1 does not mark a cell type
            cand("C465", "tiny", 0.1), cand("C465", "tiny_other", 0.0001),
            cand("C465", "tiny_other_2", 0.0001),
        ])
        self.assertEqual(a.level, "C66")


class RegionTest(unittest.TestCase):
    def test_region_filters_candidates_and_is_never_the_evidence(self):
        """Every mapping an earlier adversarial pass refuted had a consistent
        region. Anatomy identifies the nucleus, not the cell type."""
        a = assign_population("p", "Gene", [
            cand("C185", "in_region", 3.0, region=100.0),
            cand("C185", "out_of_region", 99.0, region=20.0),
        ])
        self.assertEqual(a.clusters, ("in_region",))

    def test_nothing_in_the_region_is_its_own_outcome(self):
        a = assign_population("p", "Gene",
                              [cand("C185", "x", 9.0, region=10.0)])
        self.assertEqual(a.outcome, "refused_no_region_match")

    def test_the_share_boundary_is_inclusive(self):
        a = assign_population("p", "Gene", [
            cand("C185", "exactly_at", 3.0, region=MIN_REGION_SHARE),
            cand("C185", "other", 0.1, region=MIN_REGION_SHARE),
            cand("C185", "other_2", 0.1, region=MIN_REGION_SHARE)])
        self.assertEqual(a.clusters[0], "exactly_at")


class KeepFractionTest(unittest.TestCase):
    def test_a_population_spanning_several_clusters_keeps_them_all(self):
        a = assign_population("p", "Gene", [
            cand("C465", "a", 4.0), cand("C465", "b", 3.0),
            cand("C465", "c", 2.1), cand("C465", "far_below", 0.3),
            cand("C465", "far_below_2", 0.3)])
        self.assertEqual(a.clusters, ("a", "b", "c"))

    def test_the_kept_clusters_come_back_strongest_first(self):
        a = assign_population("p", "Gene", [
            cand("C465", "weak", 2.1), cand("C465", "strong", 4.0),
            cand("C465", "middle", 3.0), cand("C465", "out", 0.1),
            cand("C465", "out_2", 0.1)])
        self.assertEqual(a.clusters, ("strong", "middle", "weak"))


class ContrastTest(unittest.TestCase):
    def test_contrast_is_scale_free_so_ambient_rna_cannot_inflate_it(self):
        """The defect that made an absolute coverage check useless.

        Ambient RNA lifts every cluster's mean together. A ratio is unchanged
        by that; an absolute threshold or a mass share is not.
        """
        plain = assign_population("p", "Gene", [
            cand("C66", "a", 3.0), cand("C66", "b", 0.3), cand("C66", "c", 0.3)])
        lifted = assign_population("p", "Gene", [
            cand("C66", "a", 30.0), cand("C66", "b", 3.0), cand("C66", "c", 3.0)])
        self.assertEqual(plain.clusters, ("a",))
        self.assertAlmostEqual(plain.contrast, lifted.contrast)

    def test_a_pan_marker_that_marks_the_whole_region_lands_in_tier_C(self):
        """`zi_gaba` on Slc32a1. Every cluster in the zona incerta is
        GABAergic, so the marker keeps all of them and picks nothing out."""
        a = assign_population("zi_gaba", "Slc32a1", [
            cand("C465", "a", 1.75, named=False), cand("C465", "b", 1.3, named=False),
            cand("C465", "c", 1.2, named=False), cand("C465", "d", 1.1, named=False)])
        self.assertEqual(len(a.clusters), 4)
        self.assertEqual(a.outcome, "mapped")
        self.assertLess(a.contrast, MIN_CONTRAST)
        self.assertEqual(a.tier, "C")

    def test_keeping_every_candidate_means_the_marker_did_not_discriminate(self):
        """The pan-marker case has nothing to stand apart from.

        Contrast is 1.0 rather than infinity: a marker that keeps the whole
        region has not picked a population out of it, and that is the opposite
        of confidence.
        """
        a = assign_population("p", "Gene",
                              [cand("C185", "only", 3.0, named=False)])
        self.assertEqual(a.contrast, 1.0)
        self.assertEqual(a.tier, "C")

    def test_contrast_measures_the_kept_set_against_what_was_excluded(self):
        """Not best-against-all-others, which counts a population's own
        siblings against it. `arc_agrp` spans five C465 clusters."""
        a = assign_population("p", "Gene", [
            cand("C465", "sib_a", 4.0), cand("C465", "sib_b", 3.0),
            cand("C465", "sib_c", 2.1),
            cand("C465", "out_a", 0.1), cand("C465", "out_b", 0.1)])
        self.assertEqual(len(a.clusters), 3)
        self.assertAlmostEqual(a.contrast, 30.0)   # median(3.0) / median(0.1)


class TierTest(unittest.TestCase):
    def test_tier_A_requires_the_atlas_authors_to_agree_on_every_cluster(self):
        a = assign_population("p", "Gene", [
            cand("C465", "Gene.Thing-1", 4.0, named=True),
            cand("C465", "Gene.Thing-2", 3.0, named=True),
            cand("C465", "other", 0.05, named=False),
            cand("C465", "other_2", 0.05, named=False)])
        self.assertEqual(a.tier, "A")
        self.assertEqual(a.named_clusters, 2)

    def test_partial_name_support_is_tier_B(self):
        a = assign_population("p", "Gene", [
            cand("C465", "Gene.Thing", 4.0, named=True),
            cand("C465", "Unrelated", 3.0, named=False),
            cand("C465", "floor", 0.05, named=False),
            cand("C465", "floor_2", 0.05, named=False)])
        self.assertEqual(a.tier, "B")

    def test_no_name_support_but_a_top_ranked_marker_is_tier_B(self):
        """`pvh_avp`: the clusters are named Oxt/Th/Sim1, yet Avp is the single
        highest-expressed gene in all three. Magnocellular PVH really does
        co-express both, so this is ambiguity to record, not an error."""
        a = assign_population("pvh_avp", "Avp", [
            cand("C465", "Oxt.Ebf3.Caprin2.GLU-6", 7.679, rank=1, named=False),
            cand("C465", "other", 0.5, rank=900, named=False),
            cand("C465", "other_2", 0.5, rank=900, named=False)])
        self.assertEqual(a.tier, "B")
        self.assertEqual(a.best_marker_rank, 1)

    def test_neither_names_nor_rank_is_tier_C(self):
        a = assign_population("p", "Gene", [
            cand("C465", "x", 4.0, rank=GOOD_RANK + 1, named=False),
            cand("C465", "y", 0.05, rank=5000, named=False),
            cand("C465", "z", 0.05, rank=5000, named=False)])
        self.assertEqual(a.tier, "C")

    def test_full_name_support_survives_low_contrast_as_tier_B(self):
        """`poa_galanin` has Gal in all four cluster names and contrast 9.01.

        The atlas authors naming every cluster for the marker is evidence this
        rule did not generate, so a round-number contrast threshold must not
        discard it -- an earlier gate filed Gal beside `zi_gaba`, whose marker
        keeps the whole zona incerta and picks nothing out of it at all.
        """
        a = assign_population("p", "Gene", [
            cand("C465", "Gene.a", 2.0, named=True),
            cand("C465", "Gene.b", 1.9, named=True),
            cand("C465", "out", 0.3, named=False)])
        self.assertLess(a.contrast, MIN_CONTRAST)
        self.assertEqual(a.tier, "B")

    def test_name_support_cannot_reach_tier_A_without_contrast(self):
        """Tier A needs both lines of evidence, so B is the ceiling here."""
        a = assign_population("p", "Gene", [
            cand("C465", "Gene.a", 2.0, named=True),
            cand("C465", "Gene.b", 1.9, named=True),
            cand("C465", "out", 0.3, named=False)])
        self.assertNotEqual(a.tier, "A")

    def test_rank_is_recorded_but_never_gates_a_mapping(self):
        """Crh ranks 544th inside the cluster HypoMap itself calls
        `Crh.Tent5a.GLU-2`. Rank tracks transcript abundance, which differs by
        orders of magnitude between a neuropeptide and a transcription factor,
        so gating on it would silently select for abundance."""
        a = assign_population("pvh_crh", "Crh", [
            cand("C465", "Crh.Tent5a.GLU-2", 0.859, rank=544, named=True),
            cand("C465", "other", 0.008, rank=9000, named=False),
            cand("C465", "other_2", 0.008, rank=9000, named=False)])
        self.assertEqual(a.outcome, "mapped")
        self.assertEqual(a.tier, "A")
        self.assertEqual(a.best_marker_rank, 544)

    def test_every_tier_returned_is_declared(self):
        self.assertIn(assign_population(
            "p", "Gene", [cand("C66", "Gene.x", 3.0), cand("C66", "y", 0.1),
                          cand("C66", "z", 0.1)]).tier, TIERS)


class NestingTest(unittest.TestCase):
    def test_overlap_across_LEVELS_is_seen_though_no_name_is_shared(self):
        """The bug an adversarial pass found and these tests had not.

        `nesting` compared cluster NAME strings, so two populations mapped at
        different levels were invisible to it however many cells they shared.
        In the real mapping it reported 1 of 11 cell-level overlapping pairs;
        `arc_agrp` (C66) and `arc_tbx3` (C185) share 4663 of `arc_agrp`'s 5244
        cells -- 88.9% -- and not one cluster name.
        """
        coarse = Assignment("arc_agrp", "mapped", level="C66",
                            clusters=("C66-46: Agrp.GABA-4",),
                            leaf_clusters=("C465-286", "C465-287", "C465-288"))
        finer = Assignment("arc_tbx3", "mapped", level="C185",
                           clusters=("C185-115: Npy.Agrp.GABA-4",),
                           leaf_clusters=("C465-287", "C465-288", "C465-999"))
        self.assertEqual(
            {c.split(": ")[0] for c in coarse.clusters}
            & {c.split(": ")[0] for c in finer.clusters},
            set(), "the fixture must share no cluster name, or it proves nothing")
        self.assertEqual(nesting([coarse, finer]),
                         (("arc_agrp", "arc_tbx3", "partial"),))

    def test_leaves_are_preferred_over_names_when_both_are_present(self):
        """Names can agree while cells do not, at different levels."""
        a = Assignment("a", "mapped", clusters=("shared-name",),
                       leaf_clusters=("leaf-1",))
        b = Assignment("b", "mapped", clusters=("shared-name",),
                       leaf_clusters=("leaf-2",))
        self.assertEqual(nesting([a, b]), ())

    def test_names_are_still_used_when_no_leaves_are_supplied(self):
        a = Assignment("a", "mapped", clusters=("x", "y"))
        b = Assignment("b", "mapped", clusters=("y",))
        self.assertEqual(nesting([a, b]), (("a", "b", "superset"),))

    def test_a_subset_relation_is_reported_in_both_directions_once(self):
        """`arc_pomc` inside `arc_tbx3` is the biology, not a bug -- but two
        populations sharing clusters have expression vectors correlated by
        construction, so a downstream similarity ranking must be able to see it."""
        small = Assignment("arc_pomc", "mapped", clusters=("a", "b"))
        big = Assignment("arc_tbx3", "mapped", clusters=("a", "b", "c"))
        self.assertEqual(nesting([small, big]), (("arc_pomc", "arc_tbx3", "subset"),))

    def test_partial_overlap_is_distinguished_from_nesting(self):
        one = Assignment("p", "mapped", clusters=("a", "b"))
        two = Assignment("q", "mapped", clusters=("b", "c"))
        self.assertEqual(nesting([one, two]), (("p", "q", "partial"),))

    def test_disjoint_sets_report_nothing(self):
        one = Assignment("p", "mapped", clusters=("a",))
        two = Assignment("q", "mapped", clusters=("b",))
        self.assertEqual(nesting([one, two]), ())

    def test_refused_populations_are_not_compared(self):
        one = Assignment("p", "mapped", clusters=("a",))
        two = Assignment("q", "refused_below_floor", clusters=())
        self.assertEqual(nesting([one, two]), ())


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
