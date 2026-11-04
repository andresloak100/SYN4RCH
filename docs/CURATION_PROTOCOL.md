# SYN4RCH curation protocol (v0.3.0)

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

The reader — human or model — is given one passage: its text and its locator.
From that passage alone, and nothing else, the reader proposes every record it
supports. A passage reporting three readouts for one manipulation yields three
records, because SCHEMA.md 1 makes one observation the unit of record and
forbids aggregating across outcomes (`docs/DECISIONS.md` D-020). A passage
reporting no result yields none.

**Required of every proposal.**

1. `support_quote` is copied, character for character, from the passage.
2. The `direction` must follow from `support_quote` read on its own, with no
   other sentence and no background knowledge. If it does not, there is no
   record.
3. `source_locator` is the unit's locator, unmodified.
4. Any field the passage does not report is `null`. Not a default, not a
   typical value, not a value from elsewhere in the paper. The one controlled
   exception is the condition dimension, covered in 3b below.
5. `perturbation_class` and `perturbation_polarity` are **derived, never
   judged**. Take them from `syn4rch.schema.rules`, or leave them null and let
   the gate fill them. A hand-written `perturbation_class` that disagrees with
   `method` is refused outright (`S007`); a hand-written
   `perturbation_polarity` that disagrees with the rule is admitted with a
   warning (`W013`), because for a locally infused drug a curator who has read
   the pharmacology may genuinely know the sign when `valence` does not. Where
   you do override, say why in `extraction_notes`.
6. For a `molecular` or `pharmacological` perturbation, `target_raw` is the
   paper's own name for what was manipulated and `target_id` resolves it
   against `data/vocab/targets.json`. Neither may be guessed from context: an
   ambiguous match is `null`, and `null` is reported (`W014`) rather than
   filled in.

**Forbidden.**

- Quoting a figure, panel or table the reader cannot read as text.
- Merging two passages into one record, or one passage into two records that
  restate the same observation.
- Resolving `valence` to `activation` or `inhibition` when the passage leaves
  the net circuit effect open. `unclear` is a correct answer and must not be
  coerced (SCHEMA.md 3.4).
- Recording a `direction` that describes the assay number rather than the
  quantity `outcome_id` names, where the two have opposite polarity. Falling
  glucose in an insulin tolerance test is *rising* insulin sensitivity
  (`docs/DECISIONS.md` D-019).
- Dropping a null result. `no_change` is recorded explicitly, and omitting it
  is the single most effective way to put a direction bias into this dataset
  (SCHEMA.md 3.5).
- Proposing any record for a paper whose text was not retrieved (SCHEMA.md
  rule 2.5). Abstract-only papers are carried as candidates, never as records.
- Setting `target_raw` or `target_id` on a circuit or structural perturbation.
  A chemogenetic actuator is an `actuator`, not a molecular target, and the
  confusion corrupts `record_id` (SCHEMA.md 3.4b, refused as `S008`).
- Writing a condition the paper does not state anywhere, including
  `not_reported` when nothing was checked. `not_reported` is an assertion that
  the paper was read and is silent; a migration or a reader who did not look
  writes `null`.
- Letting the derived polarity stand silently on a **dominant-negative or
  catalytically dead allele**. Overexpressing one is `increase` by method and a
  loss of function in fact, so name the allele in `targeting_strategy` and say
  so in `extraction_notes`. `W016` catches the wordings it knows; it is not a
  substitute for the reader noticing (`docs/DECISIONS.md` D-032).
- Reading a `perturbation_polarity` as though it described the readout. The
  polarity and the `direction` share the words `increase` and `decrease` and
  mean different things: the polarity is what was done inside the brain, the
  direction is what the body did (`docs/DECISIONS.md` D-033).
- Applying a Methods-stated condition to a paper's records when the Methods
  describe **two** cohorts on that dimension. Both diets of a diet-comparison
  study appear in its Methods, and neither is the answer for any particular
  figure (`docs/DECISIONS.md` D-038).
- Reading a husbandry sentence as a `challenge`, or an assay name as its own
  `challenge` (`docs/DECISIONS.md` D-040, D-041). Both are set out in 3b.
- Writing a condition value with no `condition_evidence` entry behind it. The
  value may be right and a reader still cannot check it without re-reading the
  paper; `W017` counts how much of the corpus is in that state.

## 3b. Stage 6: the condition dimension

The condition dimension is `(sex, diet, age_group, challenge)` and it is the
one place where a value may come from outside the record's own passage.

Papers put the result in Results and the sex and diet in Methods. Reading the
passage alone therefore guaranteed the gap this version exists to close: 197 of
the first 362 records carried no diet and 224 no sex, and the curators were
right every single time. So a condition **may** be taken from the paper's
Methods, on three conditions, all carried by a `condition_evidence` entry on the
record (v0.3.0):

1. `source = methods`, stated, not implied.
2. `locator` names the Methods passage it came from. Its absence is a HARD
   finding (`A008`) and the record is refused.
3. `quote` is the sentence itself, verbatim. Its absence is a HARD finding
   (`A009`). An entry asserting that evidence exists while withholding it is
   worse than a `null`, because it reads as checked.

The entry also names the `dimension` it supports and repeats the `value`, and
that value must equal the field (`A011`). One entry per filled dimension: a
record may take its sex from Methods and its challenge from its own passage, and
each says so separately.

Everything else about rule 4 stands. A condition that no passage of the paper
states is `null`. `condition_evidence` widens where a value may be read *from*;
it does not widen what may be written down. A reader who infers a diet from a
figure legend's y-axis, from a strain's typical husbandry, or from what the lab
usually does has fabricated it.

**A Methods value applies only where the paper holds the dimension constant.**
If the Methods describe one cohort on a dimension -- "male C57BL/6 mice" -- that
value is true of every figure in the paper and may be written on every record
from it. If they describe two -- a chow group and a high-fat group -- the Methods
have told you both diets exist and nothing about which group produced Figure 3C.
Then the dimension stays `null` on every record whose own passage does not say,
however obvious the answer feels. See D-038.

`challenge` is the acute state at the moment of measurement — fasted, refed,
glucose-loaded, cold-exposed — and is not `diet`, which is the maintenance diet
over weeks. In this literature the challenge frequently decides the sign of the
result, which is why it is now a field and not a note.

**`challenge` is almost never a Methods value.** It is a property of a
measurement, and a Methods section describing a fast, a glucose load and a
metabolic cage run does not say which one Figure 3C reports. Two specific
readings are forbidden outright:

* **Husbandry is not a challenge.** "Mice were allowed ad libitum access to
  water and standard chow" says how the animals were kept, not the state any
  measurement was made in. Nearly every paper here contains such a sentence, so
  accepting it would fill `challenge` across most of the corpus and report the
  result as near-complete coverage. It is the most damaging mistake available in
  this stage: wrong in the flattering direction, wrong on hundreds of records at
  once, and re-identifying every one of them. See D-040.
* **An assay does not supply its own challenge.** A glucose tolerance test is a
  glucose load by definition, so writing `glucose_load` on a
  `glucose_tolerance` record restates the outcome column. It separates nothing
  that `outcome_id` does not already separate, while raising the apparent
  coverage. See D-041.

What `challenge` is *for* is the measurement whose assay does not fix its state:
food intake, body weight, energy expenditure, blood glucose. "Food intake in the
fed state" and "food intake after a fast" are two different experiments, and
telling them apart is the entire value of the field. Those come from the
record's own passage, where the paper says it.

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
- The standing verbatim gate (`syn4rch.curate.verify.gate`) reads
  `source_locator` and `support_quote` only. It does not re-check the
  `locator` and `quote` on a `condition_evidence` entry, so a condition
  evidence entry written or edited outside the condition fill is not caught by
  the gate every record passes through.
  The condition fill checks its own output and checks it harder than the gate
  would: `syn4rch.curate.fill_conditions.check_evidence_locators` resolves each
  entry's locator to the passage it names and confirms the quote occurs *in
  that passage*, not merely somewhere in the paper, and refuses to write if any
  entry fails. So every entry in the dataset today has been checked for
  substance; what is missing is a standing check that keeps that true for
  entries added by some other route (`docs/DECISIONS.md` D-024, D-037).
- Nothing detects a polarity inversion mechanically. The outcomes vocabulary
  annotates the known traps (D-021) and the perturbation key keeps molecular
  and circuit signs apart (D-025), but a reader who states the wrong direction
  for an un-annotated assay will pass every check (D-019).
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
