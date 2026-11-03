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

## The four data layers never mix

    data/
      raw/          harvest dumps, append-only, never hand-edited   [untracked]
      interim/      screening decisions, parsed texts               [untracked]
      curated/      records.jsonl  <- THE dataset                   [tracked]
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

    docs/SCHEMA.md         binding data contract (v0.1.0)
    docs/DECISIONS.md      append-only decision log
    src/syn4rch/schema/    enums, TierARecord, derived rules, validation, JSONL IO
    src/syn4rch/qc/        duplicates, conflicts, units, conditions, report, CLI
    src/syn4rch/vocab/     controlled-vocabulary loader and resolver
    src/syn4rch/harvest/   literature retrieval
    src/syn4rch/curate/    extraction and curation tooling
    src/syn4rch/feasibility/  can this question be answered from this corpus?
    src/syn4rch/models/    modelling experiments (never a source of evidence)
    src/syn4rch/kg/        knowledge graph (association and anatomy classes)
    tests/                 unittest suite; fabricated fixtures only
    reports/               generated QC reports (untracked)
