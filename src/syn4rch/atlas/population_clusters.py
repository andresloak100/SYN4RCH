"""Assign a declared population to a set of HypoMap transcriptomic clusters.

This is a **Tier B annotation layer**: it says which clusters in a public mouse
atlas most plausibly correspond to a population the vocabulary already declares.
It is never a training target, and a mapping here is not evidence that any
intervention did anything.

The decision logic lives here as a pure function over a prepared candidate
table, so it can be tested without the 4 GB atlas; the extraction that builds
that table is the caller's job (see `scripts/` and D-065).

Four things this module exists to prevent, each of which actually happened
while the first version of this mapping was made by hand:

* **An absolute coverage fraction is not a validity check.** "What share of the
  atlas's marker-positive cells does this cluster set hold?" rejected 12 of 12
  mappings, because droplet data spreads a highly-expressed neuropeptide into
  tens of thousands of glial droplets as ambient RNA. Pmch read as "positive"
  in 58,480 of 384,925 cells; 15% of the hypothalamus is not MCH neurons. Mass
  share fails the same way for the same reason: ambient contributes little per
  cell but there are vastly more contaminated cells, so glia dominate the
  total. The usable quantity is **contrast** -- the best candidate's mean over
  the median of the other candidates in the same region -- which is scale-free
  and so cannot be inflated by ambient.

* **A rule pinned to one cluster level manufactures refusals.** `pvh_crh` has
  no PVH-majority Crh cluster at C185 at all; it resolves only at C286/C465.
  So the level is chosen, per population, as the one where contrast is
  greatest -- not the finest, which with less averaging always wins on raw
  magnitude and would fragment every population.

* **Within-cluster rank measures abundance, not specificity.** Crh ranks 544th
  of 51,216 genes inside the very cluster HypoMap's own authors named
  `Crh.Tent5a.GLU-2`. Neuropeptides (Pmch, Oxt, Avp, Hcrt, Pomc) are orders of
  magnitude more abundant than the transcription factors and receptors
  (Tbx3, Sim1, Nos1, Crh, Mc4r) that define other populations, so a rank
  threshold silently selects for abundance. Rank is recorded, never gated on.

* **Anatomy identifies a nucleus, not a cell type.** Every one of the three
  mappings refuted by an earlier adversarial pass had a consistent region.
  Region is therefore a *filter* on candidates, never the evidence.

The tier is the output a caller should actually read. Tier A means HypoMap's own
cluster names carry a marker the vocabulary declares for that population -- the
atlas authors' judgement, independent of every threshold chosen here.
"""

from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Sequence, Tuple

#: A candidate's region must be at least this much of the cluster, in percent.
#: Nearly half of HypoMap's clusters carry a substantial `NA` region share, so
#: this filter is weaker than it looks and is not treated as evidence.
MIN_REGION_SHARE = 50.0

#: A marker mean below this does not mark a cell type. Calibrated against the
#: observed split: cleanly-marked neuropeptide populations reach 1.8-6.3, while
#: receptor- and kinase-defined ones (Lepr 0.074, Mc4r 0.055, Ikbkb 0.328) and
#: transcription factors (Nr5a1 0.391) never approach it.
MARKER_FLOOR = 0.5

#: Keep every candidate within this fraction of the best one, so a population
#: that genuinely spans several clusters is not cut down to its single peak.
KEEP_FRACTION = 0.5

#: Contrast at or above this is required for tier A or B.
MIN_CONTRAST = 10.0

#: Marker rank good enough to earn tier B without any name support.
GOOD_RANK = 50

OUTCOMES = (
    "mapped",
    "refused_no_marker_in_atlas",
    "refused_no_region_match",
    "refused_below_floor",
)

TIERS = ("A", "B", "C")


@dataclass(frozen=True)
class Candidate:
    """One cluster considered for one population, at one level."""

    level: str
    cluster: str
    marker_mean: float
    region_share: float
    #: Rank of the marker among all genes inside this cluster. Recorded for the
    #: reader; never gated on, because it tracks abundance (see module docstring).
    marker_rank_in_cluster: int
    #: Does HypoMap's own name for this cluster carry a marker the vocabulary
    #: declares for this population? The atlas authors' judgement, not ours.
    name_carries_marker: bool


@dataclass(frozen=True)
class Assignment:
    population_id: str
    outcome: str
    marker: Optional[str] = None
    level: Optional[str] = None
    clusters: Tuple[str, ...] = ()
    contrast: Optional[float] = None
    tier: Optional[str] = None
    named_clusters: int = 0
    best_marker_rank: Optional[int] = None
    reason: Optional[str] = None

    @property
    def refused(self) -> bool:
        return self.outcome != "mapped"


def _median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    mid = len(ordered) // 2
    return (ordered[mid] if len(ordered) % 2
            else (ordered[mid - 1] + ordered[mid]) / 2.0)


def _contrast(kept: Sequence[float], excluded: Sequence[float]) -> float:
    """How far the chosen set stands above the candidates left out of it.

    Scale-free on purpose: ambient RNA raises every cluster's mean together, so
    a ratio survives it where an absolute or mass-share figure does not.

    Measured kept-vs-excluded rather than best-vs-all-others, because a
    population that genuinely spans several clusters has its own siblings in
    "all others", which drags the comparison up and penalises it for being
    real. `arc_agrp` occupies five C465 clusters; `scn_avp` seven.

    Keeping *every* candidate in the region is the pan-marker case -- Slc32a1
    across the zona incerta, where every cluster is GABAergic -- so there is
    nothing to stand apart from and the contrast is 1.0, not infinity.
    """
    if not excluded:
        return 1.0
    below = _median(excluded)
    if below <= 0:
        return float("inf")
    return _median(kept) / below


def _tier(named: int, total: int, contrast: float, best_rank: int) -> str:
    """Combine two independent lines of support; do not let either gate the other.

    Name support and contrast are different kinds of evidence: one is the atlas
    authors' judgement, the other is this rule's own measurement. An earlier
    version gated everything behind contrast, which threw away full name
    support on a round number -- `poa_galanin` has Gal in all four of its
    cluster names and contrast 9.01, and was filed next to `zi_gaba`, whose
    marker keeps the whole region and picks nothing out at all.

    Tier boundaries are convenience labels over continuous quantities, so a
    population sitting a hair either side of `MIN_CONTRAST` is not meaningfully
    different from its neighbour; callers that need the real figure should read
    `contrast` and `named_clusters` rather than the letter.
    """
    full_names = bool(total) and named == total
    if full_names and contrast >= MIN_CONTRAST:
        return "A"
    if full_names:
        return "B"                      # independent corroboration on its own
    if contrast >= MIN_CONTRAST and (named >= 1 or best_rank <= GOOD_RANK):
        return "B"
    return "C"


def assign_population(
    population_id: str,
    marker: Optional[str],
    candidates: Iterable[Candidate],
    *,
    min_region_share: float = MIN_REGION_SHARE,
    floor: float = MARKER_FLOOR,
    keep_fraction: float = KEEP_FRACTION,
) -> Assignment:
    """Choose a level and a cluster set for one population, or refuse.

    `candidates` may span several levels; the level with the greatest contrast
    wins among those whose best candidate clears `floor`.
    """
    if not marker:
        return Assignment(population_id, "refused_no_marker_in_atlas",
                          reason="no declared marker gene is present in the atlas")

    by_level: Dict[str, list] = {}
    for c in candidates:
        if c.region_share >= min_region_share:
            by_level.setdefault(c.level, []).append(c)
    if not by_level:
        return Assignment(population_id, "refused_no_region_match", marker=marker,
                          reason="no cluster reaches %.0f%% of the declared region"
                                 % min_region_share)

    scored = []
    for level, group in by_level.items():
        means = [c.marker_mean for c in group]
        top = max(means)
        kept = [c for c, m in zip(group, means) if m >= keep_fraction * top]
        left = [m for c, m in zip(group, means)
                if m < keep_fraction * top]
        contrast = _contrast([c.marker_mean for c in kept], left)
        scored.append((level, group, top, kept, contrast))

    passing = [s for s in scored if s[2] >= floor]
    if not passing:
        level, group, top, _kept, _c = max(scored, key=lambda s: s[2])
        best = max(group, key=lambda c: c.marker_mean)
        return Assignment(
            population_id, "refused_below_floor", marker=marker,
            reason="no level reaches the %.2f floor (best %.3f at %s, cluster %s)"
                   % (floor, top, level, best.cluster))

    level, group, top, keep, contrast = max(passing, key=lambda s: s[4])
    keep = sorted(keep, key=lambda c: -c.marker_mean)
    named = sum(1 for c in keep if c.name_carries_marker)
    best_rank = min(c.marker_rank_in_cluster for c in keep)
    return Assignment(
        population_id, "mapped", marker=marker, level=level,
        clusters=tuple(c.cluster for c in keep),
        contrast=contrast,
        tier=_tier(named, len(keep), contrast, best_rank),
        named_clusters=named, best_marker_rank=best_rank)


def nesting(assignments: Iterable[Assignment]) -> Tuple[Tuple[str, str, str], ...]:
    """Pairs of populations whose cluster sets intersect.

    Nesting is reported, not forbidden: `arc_pomc` inside `arc_tbx3` is the
    biology -- Tbx3 is expressed in POMC neurons. But two populations sharing
    clusters will have expression vectors correlated **by construction**, so
    anything downstream that ranks populations by transcriptomic similarity has
    to be able to see the overlap rather than rediscover it as a finding.
    """
    sets = {a.population_id: frozenset(a.clusters)
            for a in assignments if not a.refused}
    out = []
    for a in sorted(sets):
        for b in sorted(sets):
            if a >= b or not (sets[a] & sets[b]):
                continue
            if sets[a] <= sets[b]:
                rel = "subset"
            elif sets[b] <= sets[a]:
                rel = "superset"
            else:
                rel = "partial"
            out.append((a, b, rel))
    return tuple(out)
