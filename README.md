# SYN4RCH

A rigorously curated dataset of **hypothalamic neuronal perturbation
experiments and their systemic physiological consequences**, plus the schema
and quality-control code that defines it.

One record = one observation:

    (neuronal population) x (perturbation) x (physiological readout) x (study)

## What this is, and what it is not

SYN4RCH is a **hypothesis-generation platform**. It is not a validated
simulator of anything, and in particular it is not a rejuvenation model. It
does not score interventions, rank them, or predict what a perturbation would
do in an animal that has not been measured.

The first milestone is deliberately modest and falsifiable:

> **Determine whether the outcomes of hypothalamic neuronal interventions are
> predictable at all from the evidence that already exists in the published
> literature.**

"The available evidence is insufficient to predict these outcomes" is a
legitimate result of that work, and a publishable one. It is a finding about
the state of the field, and arriving at it is a success for this project, not
a failure. The schema and QC layers are built so that this answer *can* be
reached honestly: null results are recorded, conflicts between studies are
preserved, unmappable populations are kept unmapped, and every stored value
carries a verbatim quote and a locator so a human can go and check it.

Design commitments that follow from that, enforced in code:

* **Nothing is inferred into a field.** If a paper does not report it, the
  value is `null`. Never defaulted, never guessed (`docs/SCHEMA.md` 2).
* **Every value carries its evidence.** A record without a `support_quote` or
  a `source_locator` cannot be stored.
* **Null results are kept.** `direction = no_change` is recorded explicitly;
  dropping null results is a leading source of direction bias and is
  forbidden.
* **Conflicts are data.** When two studies disagree about the same claim, QC
  records the conflict on both sides and *resolves nothing*. Nothing in this
  codebase averages, prefers or discards either side.
* **`unclear` is a valid answer.** A perturbation with no clear net valence is
  recorded as `unclear` rather than coerced into activation or inhibition.
* **Leakage control defines the benchmark.** A cross-validated number reported
  without naming its grouping (study, lab, population) is treated as an
  invalid result (`docs/SCHEMA.md` 5).

`docs/SCHEMA.md` is the binding contract for the data. `docs/DECISIONS.md` is
an append-only log of every judgement call made while implementing it,
including the places where SCHEMA.md was underspecified and how that was
resolved.

## Where this stands

The first milestone has a partial answer, and it is not the one this project
started with.

The dataset holds **362 Tier A records from 42 studies**, curated from a corpus
of 9,653 screened papers of which 974 full texts were retrievable. Every record
carries a verbatim quote that was mechanically verified against the retrieved
source; 100% of them are still `review_status: unreviewed`, which the reports
state rather than hide.

**Direction of effect is predictable in one half of the corpus and not in the
other.** Leave-one-study-out, frequency majority, two independent
implementations agreeing:

| subset | rows with evidence in their cell | studies | guess the mode | cell's modal direction |
|---|---|---|---|---|
| circuit manipulation (opto-/chemogenetic) | 116 | 20 | 0.397 | **0.517** |
| molecular manipulation (knockout / overexpression) | **0 of 170** | — | — | nothing to predict from |

Pooled, nothing beats guessing, which is what every earlier benchmark here
reported. Split, the circuit half gains 12 points (study-clustered bootstrap:
95% interval [+0.008, +0.271]; row-level McNemar p=0.156, so suggestive rather
than established).

The molecular half is not *bad* at prediction, it has nothing to predict from:
with the manipulated gene in the cell key, not one of its 170 rows has its cell
populated in training, because each gene appears in exactly one study.

The same schema change also cut the number of claims reported by more than one
study from 39 to 11. That is not a loss of data; it is the discovery that most
of the apparent replication here came from pooling a male cohort with a female
one and one gene with another. This corpus holds eleven independently
replicated claims.

Two things explain the predictability split, and both are now schema-level facts
rather than suspicions:

* **The molecular records are 17 genes across 16 studies and 16 laboratories,
  and no gene is manipulated by more than one laboratory.** So under study- or
  lab-grouped validation a held-out molecular record's perturbation was never
  seen in training, and there is nothing to have learned. The axis is at least
  internally consistent: where the *same* gene was both knocked out and
  overexpressed, 6 of 9 cells mirror cleanly and none contradicts. But all 9 of
  those cells are one laboratory and one study, so that consistency is equally
  consistent with a paper telling a coherent story, and nothing here separates
  the two. `syn4rch.feasibility.identity` measures the reach;
  `docs/DECISIONS.md` D-027 records both halves.
* **Until schema v0.2.0 there was no field that could tell the two halves
  apart.** All 170 molecular records carried `valence: unclear`, which was
  indistinguishable from a circuit manipulation whose net effect the paper did
  not resolve. v0.2.0 did not create the circuit signal — on that subset the old
  and new keys score identically — it made the signal *findable* (D-025, D-028).

One more honest result, since it goes against the change: **conditioning the
model on sex, diet and metabolic state makes prediction slightly worse**, by
−0.006 pooled. Four fifths of the conditioned cells hold a single record, so the
split fragments the evidence faster than it sharpens it. The condition dimension
is kept because it is what removed the false contradictions, not because it
predicts better (`docs/DECISIONS.md` D-036).

## What the condition curation pass changed (schema v0.3.0)

D-036's verdict above was measured on a field that was empty. A curation pass has
now read all 42 papers and filled what their text supports, with a verbatim
sentence and a locator on every value (`reports/condition_fill_0.3.0.md`):

| dimension | null before | null after | values carrying evidence |
|---|---|---|---|
| `sex` | 224 | 112 | 174 of 235 |
| `diet` | 197 | 111 | 104 of 236 |
| `age_group` | 286 | 250 | 36 of 43 |
| `challenge` | **362** | 337 | 25 of 25 |

**The verdict did not change.** Conditioned against axis-only, under `study_id`:
−0.006 [−0.021, +0.006], 0.479 against 0.485, McNemar p=0.754, `no_difference` —
the same −0.006 as on empty conditions, and *worse* (−0.061) on the rows whose
conditioned cell was seen in training. Real condition information makes the
conditioned cell finer, not better informed: 203 → 268 cells, 69.5% → 81.3%
singletons. Three schema versions of work on this dimension have not made it
useful for prediction.

**What it did change is the quality of the replication evidence.** Before the
pass this corpus held 11 claims reported by more than one study, **8 of which
matched at condition specificity 0** — grouped together partly because nobody
knew their conditions. After it there are 5, and **none matches on nothing
known**. Records knowing no condition slot at all fell from 153 to 46. Cross-study
conflicting claims fell from 9 to 4 and the last directly-opposed pair to 0.

**And a prediction written down beforehand failed.** 7 multi-study claims had been
lost *only* because one paper stated a condition another omitted, and curating the
conditions was supposed to recover some. It made it worse: 7 became 10. Filling
the conditions a paper reports makes it *less* comparable with a paper reporting
none, since `unknown` matches no stated value. Partial condition curation
increases unadjudicable non-comparability; only curating both sides of a pair
reduces it (`docs/DECISIONS.md` D-047).

No claim in the corpus reaches specificity 4 of 4, so **no replication figure here
describes two experiments known to have run under the same conditions** — the same
caveat as before the pass, still true after it.

## What is being tested next, written down in advance

The replication reach of the two halves, at the level of a
`(population, perturbation, outcome)` cell:

| class | rows | labs | cells | cells with ≥2 labs |
|---|---|---|---|---|
| circuit | 170 | 26 | 56 | **17** |
| molecular | 170 | 16 | 132 | **0** |

290 uncurated circuit papers already have retrieved full text. The plan is to
curate 60 of them and evaluate **once**: the circuit direction comparison must
reach verdict `established` (today `suggestive`, +0.121 [+0.008, +0.271],
McNemar p=0.156), and multi-lab circuit cells must rise from 17 to ≥40.

The target, the stopping point and the results that would abandon the plan are
fixed in `docs/DECISIONS.md` D-048, committed before the first new paper is read,
because a target chosen after seeing the data is not a target. **"Adding 60
circuit papers did not make direction of effect predictable" is a result**, and
reporting it is what this project is for.

So "the evidence is insufficient" was an artifact of a missing field for half
the corpus, and remains the right answer for the other half. What to curate next
changed accordingly: more optogenetic and chemogenetic papers, whose
perturbation identity recurs across laboratories by construction, rather than
more single-gene papers.

**What is not claimed.** 0.52 on a three-way problem, over 116 rows and 20
studies, with a row-level paired test that does not reach significance
(McNemar p=0.156) while the study-level bootstrap excludes zero. That is a first
signal above baseline, useful for ranking hypotheses. It is not a predictor of
what a perturbation would do in an animal nobody has measured, and a third of
circuit rows still have no evidence in their cell at all.

## The four data layers never mix

    data/
      raw/          harvest dumps, append-only, never hand-edited   [untracked]
      interim/      screening decisions, parsed texts               [untracked]
      curated/      records.jsonl  <- THE dataset                   [tracked]
                    record_id_map_*.json  old-to-new ids per migration [tracked]
      inferred/     anything derived (KG edges, imputations)        [untracked]
      predictions/  model outputs only                              [untracked]
      vocab/        controlled vocabularies (JSON)                  [tracked]

Two rules govern the layers:

1. **Every value in `curated/` must trace to `raw/`** -- via `study_id`,
   `source_locator` and `support_quote`.
2. **Nothing in `inferred/` or `predictions/` may ever be cited as
   evidence.** A model output is not an observation. Keeping these in separate
   directories is what makes that distinction auditable rather than a matter
   of good intentions.

`curated/` and `vocab/` are version controlled; the other layers are
regenerable and are not (see `.gitignore`).

## Requirements

Python **3.10 or newer**, and nothing else. The schema and QC layers use the
**standard library only** -- no pydantic, no pandas, no pyyaml -- so they run
on any interpreter, with no install step and no network. See
`docs/DECISIONS.md` D-001 (stdlib) and D-006 (why vocabularies are JSON rather
than the YAML that SCHEMA.md mentions).

## Running the QC CLI

From the repository root, with no install:

    PYTHONPATH=src python3 -m syn4rch.qc.run data/curated/records.jsonl

Or writing the full report as JSON as well:

    PYTHONPATH=src python3 -m syn4rch.qc.run data/curated/records.jsonl \
        --json reports/qc.json

Useful flags: `--limit N` (truncate each section of the text report),
`--quiet` (JSON only), `--no-annotate` (report conflicts without writing
`conflicts_with` in memory), `--any-version` (inspect records written under a
different `schema_version`; off by default).

Exit status:

| code | meaning |
|---|---|
| `0` | no hard schema violations -- warnings and conflicts do not fail the run |
| `1` | at least one hard schema violation, or a line that could not be read |
| `2` | usage error (missing file, bad arguments) |

Conflicts and possible non-independence never affect the exit code. An exit
code that fired on disagreement in the literature would teach whoever runs it
to suppress disagreement.

The report separates three things, and so should any downstream use of it:

* **HARD** -- the record may not enter `data/curated/`.
* **WARNING** -- a curator should look (missing control, missing n, duplicate,
  unit mismatch, a derived field disagreeing with its rule).
* **INFO** -- a property of the literature, correctly recorded: cross-study
  conflicts, and same-lab repeats that are plausibly not independent evidence.

To try it on the fabricated fixtures:

    PYTHONPATH=src python3 -m syn4rch.qc.run tests/fixtures/records_clean.jsonl

## Running the tests

    python3 -m unittest discover -s tests

No network, no install, no fixtures outside the repository. The suite stubs
the controlled-vocabulary interface, so it passes whether or not
`data/vocab/*.json` has been populated yet.

## Synthetic data only, everywhere

Every example record in this repository is **fabricated**, and every
fabricated `study_id` carries the `synthetic:` prefix (`synthetic:test-001`)
so that it can never be mistaken for evidence in a file, a report or a `grep`.
`tests/fixtures/README.md` repeats this warning next to the files themselves.

No identifier, quote, number, author or biological claim anywhere in the code,
the tests or this README refers to a real publication. Where a realistic
example would have been clearer, a synthetic one is used instead.

## Layout

    docs/SCHEMA.md         binding data contract (v0.3.0)
    docs/DECISIONS.md      append-only decision log
    src/syn4rch/schema/    enums, TierARecord, derived rules, validation, JSONL IO
    src/syn4rch/qc/        duplicates, conflicts, units, conditions, report, CLI
    src/syn4rch/vocab/     controlled-vocabulary loader and resolver
    src/syn4rch/harvest/   literature retrieval
    src/syn4rch/curate/    extraction, curation, conditions, migrations
    src/syn4rch/feasibility/  can this question be answered from this corpus?
    src/syn4rch/models/    modelling experiments (never a source of evidence)
    src/syn4rch/kg/        knowledge graph (association and anatomy classes)
    tests/                 unittest suite; fabricated fixtures only
    reports/               generated QC reports (untracked)
