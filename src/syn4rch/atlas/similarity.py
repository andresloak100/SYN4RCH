"""Population expression vectors and cosine similarity, fixed per D-066.

Everything a result could be tuned against is a module constant declared before
the first evaluation, because D-059 requires the gene set and `k` to be fixed in
code and run once. The decision logic is pure functions over arrays so it can be
tested without the 4 GB atlas.

**Tier B throughout, and never a training target.** These vectors come from an
external observational atlas and are used to generalise a direction label that
still comes only from interventional experiments.

Three diagnostics are computed and returned whether or not anything asks for
them. Each is cheap now and unrecoverable afterwards, and each separates two
failures that would otherwise look identical in a delta:

* `live_genes` -- how much of the 2000-gene vector is actually informative for a
  population. A cosine of 0.8 computed over 40 live genes is not the same claim
  as 0.8 over 900, and without this a miss on a one-record population cannot be
  told apart from "this population never had a usable vector".
* `jaccard` on the cluster sets of a (held-out, neighbour) pair -- D-066
  disqualifies an overlapping neighbour, and the figure is logged so the rule
  can be seen to have bitten or not.
* `n_neighbours_available` and the nearest neighbour's identity -- required by
  D-059's second amendment.

A note on what the gene set is NOT selected by. A second session argued the
selection is biased toward well-represented regions, and specifically that
HypoMap underrepresents the VMH. Measured, that premise is false: at C465 the
VMH holds 45 clusters and 33,694 cells (8.8% of the atlas), effectively tied
with the arcuate's 45 and 37,667 (9.8%). The thin region among the populations
that carry records is the **PVH**, with 16 clusters and 8,710 cells (2.3%).
This matters beyond bookkeeping: `vmh_sf1` is refused for `Nr5a1` peaking at
0.391, and it cannot be blamed on sampling when the VMH is the second-best
sampled named region in the atlas. The cause is transcript abundance, which is
the point of D-065's "a Cre line is not a marker". The diagnostic below is
worth having regardless of which account of the bias is right, because it
measures the vector directly rather than reasoning about its provenance.
"""

from dataclasses import dataclass, field
from math import sqrt
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

#: Number of genes retained, by variance of per-cluster mean across every
#: cluster at the finest level. Declared in D-066. An unrestricted cosine over
#: all 51,216 genes measures library composition rather than cell identity:
#: `Malat1`, `Meg3` and `Snhg11` top nearly every cluster in this atlas.
N_GENES = 2000

#: Technical covariates rather than identity, excluded before selection.
EXCLUDED_PREFIXES = ("mt-", "Rps", "Rpl")

#: Neighbours consulted for a prediction. One, because with at most eight
#: neighbours available -- several contributing a single record -- a larger `k`
#: averages over populations that have no record in the relevant cell anyway,
#: and `k = 1` makes the nearest-neighbour log identical to the cause of the
#: prediction. D-066 forecloses letting any other `k` touch the primary.
K_NEIGHBOURS = 1

#: Exploratory only; labelled as such and cannot change the primary verdict.
K_EXPLORATORY = 3

#: A gene counts as live for a population if its mean exceeds this. A fifth of
#: D-065's marker floor: low enough to count real but modest expression, high
#: enough to exclude the ambient trickle that makes absolute thresholds
#: meaningless at the single-cell level.
LIVE_GENE_FLOOR = 0.1

#: D-067: each gene is centred across the populations being compared before the
#: cosine. An un-centred cosine over non-negative vectors sharing a large
#: component is always near 1 -- measured on the real pool, every pair fell
#: between 0.893 and 0.974 and the gap between a row's first and second nearest
#: neighbour was 0.0051. With `k = 1` the whole prediction turned on that gap.
CENTRE_GENES = True

#: D-067's tie guard. A row is predicted only if its nearest neighbour survives
#: resampling the gene set: the same population must come first in at least
#: `STABILITY_MIN` of `BOOTSTRAP_REPLICATES` replicates. A row that fails is
#: reported unranked -- unpredicted, never wrong. This is the mechanism for
#: answering "insufficient evidence to rank these reliably", which this
#: experiment had no implementation of before.
BOOTSTRAP_REPLICATES = 200
STABILITY_MIN = 0.95

#: Fixed so the pre-registered run is reproducible by anyone re-running it.
BOOTSTRAP_SEED = 20670

#: Mapping tiers admitted to the primary (D-066). The diagnostics below are
#: reported for every population and **do not** gate admission -- measuring a
#: population's vector quality and then excluding it on what we saw would be
#: the post-hoc move the pre-registration exists to prevent.
PRIMARY_TIERS = ("A",)
SECONDARY_TIERS = ("A", "B", "C")


@dataclass(frozen=True)
class PopulationVector:
    population_id: str
    tier: str
    #: Pseudobulk mean over exactly the cells the mapping assigns, so the vector
    #: is comparable across populations assigned at different cluster levels.
    values: Tuple[float, ...]
    clusters: Tuple[str, ...]
    n_cells: int
    #: Finest-level clusters, for the overlap test. Two populations mapped at
    #: different levels share cells without sharing a cluster name, so a
    #: name-based disqualification is blind to them (see
    #: `population_clusters.nesting`). D-066 disqualifies an overlapping
    #: neighbour, and it can only do that on what this field carries.
    leaf_clusters: Tuple[str, ...] = ()

    @property
    def overlap_key(self) -> Tuple[str, ...]:
        return self.leaf_clusters or self.clusters

    @property
    def live_genes(self) -> int:
        return sum(1 for v in self.values if v > LIVE_GENE_FLOOR)

    @property
    def norm(self) -> float:
        return sqrt(sum(v * v for v in self.values))


@dataclass
class NeighbourChoice:
    """Why one test row got the prediction it got."""

    population_id: str
    n_neighbours_available: int
    chosen: Optional[str] = None
    similarity: Optional[float] = None
    jaccard_with_chosen: float = 0.0
    live_genes_self: int = 0
    live_genes_chosen: int = 0
    skipped_for_overlap: Tuple[str, ...] = ()
    fell_back: bool = False
    ranked: List[Tuple[str, float]] = field(default_factory=list)


def select_from_variance(
    variances: Sequence[float],
    gene_names: Sequence[str],
    *,
    n_genes: int = N_GENES,
) -> Tuple[int, ...]:
    """The declared selection policy, given each gene's variance.

    Split out from the arithmetic so there is exactly one implementation of the
    *policy* -- exclude technical covariates, take the top `n_genes`, break ties
    on index for determinism -- while the variance itself can be computed with
    whatever is fast. The caller computing variance differently cannot change
    which genes the policy would pick from it.
    """
    scored = [
        (variances[i], i) for i, name in enumerate(gene_names)
        if not str(name).startswith(EXCLUDED_PREFIXES)
    ]
    scored.sort(key=lambda pair: (-pair[0], pair[1]))
    return tuple(sorted(i for _v, i in scored[:n_genes]))


def gene_variances(
    cluster_means: Sequence[Sequence[float]],
    n_genes_total: int,
) -> Tuple[float, ...]:
    """Population variance of each gene's per-cluster mean. Reference version."""
    n = len(cluster_means)
    if not n:
        return tuple(0.0 for _ in range(n_genes_total))
    out = []
    for i in range(n_genes_total):
        total = 0.0
        for row in cluster_means:
            total += row[i]
        mean = total / n
        var = 0.0
        for row in cluster_means:
            d = row[i] - mean
            var += d * d
        out.append(var / n)
    return tuple(out)


def select_genes(
    cluster_means: Sequence[Sequence[float]],
    gene_names: Sequence[str],
    *,
    n_genes: int = N_GENES,
) -> Tuple[int, ...]:
    """Indices of the `n_genes` highest-variance genes across all clusters.

    Depends only on the atlas: not on which populations carry records, not on
    the mapping, and not on any delta.
    """
    if not len(cluster_means):
        return ()
    variances = gene_variances(cluster_means, len(gene_names))
    return select_from_variance(variances, gene_names, n_genes=n_genes)


def cosine(a: Sequence[float], b: Sequence[float]) -> float:
    num = sum(x * y for x, y in zip(a, b))
    na = sqrt(sum(x * x for x in a))
    nb = sqrt(sum(y * y for y in b))
    if na == 0.0 or nb == 0.0:
        return 0.0
    return num / (na * nb)


def jaccard(a: Iterable[str], b: Iterable[str]) -> float:
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return 0.0
    return len(sa & sb) / len(sa | sb)


def similarity_matrix(
    vectors: Sequence[PopulationVector],
) -> Dict[Tuple[str, str], float]:
    out = {}
    for i, u in enumerate(vectors):
        for v in vectors[i + 1:]:
            s = cosine(u.values, v.values)
            out[(u.population_id, v.population_id)] = s
            out[(v.population_id, u.population_id)] = s
    return out


def choose_neighbour(
    held_out: PopulationVector,
    candidates: Sequence[PopulationVector],
    sims: Dict[Tuple[str, str], float],
    *,
    disqualify_overlap: bool = True,
) -> NeighbourChoice:
    """The nearest training population that does not share a cluster.

    D-066: an overlapping neighbour is skipped and the next nearest used,
    because two populations sharing clusters have vectors correlated by
    construction. If every candidate overlaps, the row falls back -- counted
    unpredicted, never counted wrong.
    """
    ranked = sorted(
        ((c, sims.get((held_out.population_id, c.population_id), 0.0))
         for c in candidates),
        key=lambda pair: (-pair[1], pair[0].population_id),
    )
    choice = NeighbourChoice(
        population_id=held_out.population_id,
        n_neighbours_available=len(candidates),
        live_genes_self=held_out.live_genes,
        ranked=[(c.population_id, s) for c, s in ranked],
    )
    skipped = []
    for cand, sim in ranked:
        overlap = jaccard(held_out.overlap_key, cand.overlap_key)
        if disqualify_overlap and overlap > 0.0:
            skipped.append(cand.population_id)
            continue
        choice.chosen = cand.population_id
        choice.similarity = sim
        choice.jaccard_with_chosen = overlap
        choice.live_genes_chosen = cand.live_genes
        choice.skipped_for_overlap = tuple(skipped)
        return choice
    choice.skipped_for_overlap = tuple(skipped)
    choice.fell_back = True
    return choice


def centre_genes(
    matrix: Sequence[Sequence[float]],
) -> Tuple[Tuple[float, ...], ...]:
    """Subtract each gene's mean across the populations being compared.

    D-067. This is the whole of the metric change: it removes the shared
    baseline that every neuron contributes, leaving what differs between them.
    """
    if not matrix:
        return ()
    n_pop = len(matrix)
    width = len(matrix[0])
    means = [sum(row[g] for row in matrix) / n_pop for g in range(width)]
    return tuple(tuple(row[g] - means[g] for g in range(width))
                 for row in matrix)


def bootstrap_gene_draws(
    n_genes: int,
    *,
    replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = BOOTSTRAP_SEED,
) -> Tuple[Tuple[int, ...], ...]:
    """Gene indices resampled with replacement, deterministically.

    Seeded rather than left to chance so that the single pre-registered run is
    reproducible -- a tie guard whose verdict changes between runs would decide
    which rows are reported on nothing.
    """
    import random

    rng = random.Random(seed)
    return tuple(
        tuple(rng.randrange(n_genes) for _ in range(n_genes))
        for _ in range(replicates)
    )


def nearest_stability(
    population_ids: Sequence[str],
    nearest_for_genes,
    n_genes: int,
    *,
    replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = BOOTSTRAP_SEED,
) -> Dict[str, float]:
    """Fraction of replicates in which each population's nearest is unchanged.

    `nearest_for_genes(gene_indices) -> {population_id: nearest_id}` supplies
    the arithmetic, so the policy here is the single declared one however the
    similarity is computed -- the same split as `select_from_variance`.
    """
    full = nearest_for_genes(tuple(range(n_genes)))
    agree = {pid: 0 for pid in population_ids}
    draws = bootstrap_gene_draws(n_genes, replicates=replicates, seed=seed)
    for draw in draws:
        got = nearest_for_genes(draw)
        for pid in population_ids:
            if got.get(pid) is not None and got.get(pid) == full.get(pid):
                agree[pid] += 1
    return {pid: agree[pid] / replicates for pid in population_ids}


def is_rankable(stability: float, *, minimum: float = STABILITY_MIN) -> bool:
    """Whether a row may be predicted at all (D-067).

    Below the threshold the honest output is "insufficient evidence to rank",
    not a neighbour the feature cannot defend.
    """
    return stability >= minimum
