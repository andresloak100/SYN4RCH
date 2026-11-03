# SYN4RCH curation protocol (v0.1.0)

How a published experiment becomes a Tier A record. `docs/SCHEMA.md` defines
*what* a record is; this document defines *how* one is allowed to come into
existence. Where the two appear to disagree, SCHEMA.md wins and the
disagreement is a defect to be logged in `docs/DECISIONS.md`.

The protocol exists for one reason. The single greatest threat to this dataset
is not a missing field or a wrong unit. It is a plausible-sounding record that
no paper actually supports. Every rule below is written against that threat.

---

## 1. Stages

Each stage reads the layer above it and writes only its own output. No stage
edits an earlier layer.

| # | stage | input | output |
|---|---|---|---|
| 1 | harvest | PubMed / Europe PMC / Crossref | `data/raw/`, `data/interim/candidates.jsonl` |
| 2 | screen | candidates | screening decision on each candidate |
| 3 | full-text retrieval | open-access candidates | `data/raw/fulltext/*.xml` |
| 4 | passage segmentation | retrieved full text | addressable passages with locators |
| 5 | unit building | passages + vocabularies | `data/interim/extraction_units.jsonl` |
| 6 | record proposal | one extraction unit | a *proposed* record |
| 7 | verbatim gate | proposed records + full text | admitted / rejected |
| 8 | schema validation | admitted records | hard findings and warnings |
| 9 | QC | validated records | duplicates, unit checks, conflicts |
| 10 | dataset | everything that survived | `data/curated/records.jsonl` |

Stages 1–5 and 7–9 are deterministic code. Stage 6 is the only stage that
requires reading comprehension, and it is therefore the only stage that can
fabricate. Stages 7–9 exist to catch it when it does.

## 2. An extraction unit is not a record

Stage 5 emits units: a passage of text in which a population term, a
perturbation-method term and an outcome term happen to co-occur. A unit is a
*place worth reading*. It asserts nothing.

Co-occurrence is not a finding. A unit may be an introduction sentence citing
somebody else's work, a methods sentence with no result in it, or a discussion
sentence speculating about a mechanism. The reader at stage 6 must discard all
three. The expected yield from units to records is low, and a low yield is the
protocol working.

## 3. Stage 6: proposing a record

The reader — human or model — is given one unit: its passage text and its
locator. From that unit alone the reader proposes at most one record.

**Required of every proposal.**

1. `support_quote` is copied, character for character, from the passage.
2. The `direction` must follow from `support_quote` read on its own, with no
   other sentence and no background knowledge. If it does not, there is no
   record.
3. `source_locator` is the unit's locator, unmodified.
4. Any field the passage does not report is `null`. Not a default, not a
   typical value, not a value from elsewhere in the paper.

**Forbidden.**

- Quoting a figure, panel or table the reader cannot read as text.
- Merging two passages into one record, or one passage into two records that
  restate the same observation.
- Resolving `valence` to `activation` or `inhibition` when the passage leaves
  the net circuit effect open. `unclear` is a correct answer and must not be
  coerced (SCHEMA.md 3.4).
- Dropping a null result. `no_change` is recorded explicitly, and omitting it
  is the single most effective way to put a direction bias into this dataset
  (SCHEMA.md 3.5).
- Proposing any record for a paper whose text was not retrieved (SCHEMA.md
  rule 2.5). Abstract-only papers are carried as candidates, never as records.

**Honest labelling.** `extraction_method` records what actually happened:
`llm` for an unreviewed model proposal, `llm_human_reviewed` only once a human
has read the passage and the record side by side, `human` for direct curation.
A model may never write `llm_human_reviewed`. A record whose
`extraction_method` is `llm` carries `review_status = unreviewed`, and the
dataset reports how much of itself is in that state.

## 4. Stage 7: the verbatim gate

The gate is mechanical and its verdict is final.

`support_quote` must occur as a substring of the paper's retrieved text under
exactly one normalisation: runs of whitespace collapse to a single space, and
unicode dashes, quotation marks and non-breaking spaces fold to their ASCII
equivalents. Nothing else. No case folding, no stemming, no dropped words, no
fuzzy or approximate matching. The locator must name a passage that exists in
that paper.

A proposal that fails is **rejected and recorded as rejected**. It is never
repaired, never re-quoted, never nudged into passing. Repairing a failed quote
would convert the one check that cannot be talked around into a formality.

Rejections are counted and published in `reports/curation_rejects.md`. A
rising rejection rate is information about the reader, and suppressing it
would hide the only signal the protocol has.

## 5. What this protocol does not guarantee

Stated plainly, because a reader of the dataset is entitled to know.

- The gate proves a quote is **present** in the source. It cannot prove the
  quote is **representative** of the paper, nor that the reader did not select
  the one sentence that suits a conclusion.
- The gate cannot check that the stated `direction` genuinely follows from the
  quote. That remains a reading judgement, and `extraction_method` plus
  `review_status` are the only honest record of whether a human made it.
- Effect sizes, dispersions and sample sizes transcribed from text are not
  recomputed from the underlying data. They are as right as the paper is.
- Nothing detects two papers reporting one animal cohort
  (`docs/DECISIONS.md` D-013). Lab-level grouping is the available proxy and
  it is weaker than cohort-level grouping.
- A paper behind a paywall is invisible to this dataset. The corpus is
  therefore biased toward open-access publishing, and that bias is not
  corrected anywhere.

## 6. Reproducibility obligation

Every record traces to a locator in a retrieved file in `data/raw/`
(SCHEMA.md 4). Any report that cites a count must be regenerable from
`data/curated/records.jsonl` by a committed CLI. A number that cannot be
recomputed from the committed dataset does not belong in a report.
