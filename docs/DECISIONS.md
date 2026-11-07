# Decision log

Append-only. Each entry records a choice that a reader of the code could
otherwise mistake for an accident. Never edit or delete an entry; supersede
it with a new one that names the entry it replaces.

Entries D-001 .. D-013 were written while implementing the schema and QC core
against `docs/SCHEMA.md` v0.1.0.

---

## D-001 - Standard library only for the schema and QC layers

**Decision.** `src/syn4rch/schema/**` and `src/syn4rch/qc/**` import nothing
outside the Python standard library, and target Python 3.10 syntax.

**Why.** This layer is the definition of the dataset. It must run on any
interpreter, on a laptop with no network, years from now, without resolving a
dependency tree. `pydantic`, `pandas` and `pyyaml` would each be convenient
and each would make the dataset unreadable without them.

**Cost.** Validation, serialisation and the unit checks are hand-written. No
dataframe conveniences in QC.

---

## D-002 - `record_id` hashes a named subset of the perturbation

**Decision.** `record_id = sha256(...)[:16]` over, in order:
`study_id`, `population_raw`, `method`, `valence`, `actuator`,
`outcome_raw`, `source_locator` -- joined by U+001F, NFKC-normalised,
whitespace-collapsed, case preserved.

**Why.** SCHEMA.md 3.1 specifies the hash over
`(study_id, population_raw, perturbation, outcome_raw, source_locator)`, but
`perturbation` is a *section* of the schema (3.4), not a field. The three
chosen perturbation fields are categorical and stable. Free-text fields such
as `targeting_strategy` are deliberately excluded: re-wording one during
curation would silently re-identify the record and orphan every
`conflicts_with` reference to it.

**Defect in SCHEMA.md.** 3.1 should name the perturbation fields explicitly.
Recorded here rather than changed in the doc.

**Pinned.** `tests/test_record_id.py` holds a golden digest. Changing the
input set, the separator, the normalisation or the truncation length changes
every id in the corpus and is a schema version bump.

---

## D-003 - `age_group` band edges are read as half-open intervals

**Decision.** Mouse: `<8` juvenile, `8 <= w < 27` young_adult,
`27 <= w < 65` middle_aged, `>= 65` aged.

**Why.** SCHEMA.md 3.3 gives the bands in whole weeks ("8-26", "27-64",
">=65"), which leaves 26.1-26.9 wk and 64.1-64.9 wk undefined. The half-open
reading is integer-identical to the table and defines the gaps. Boundary
values 7/8, 26/27 and 64/65 wk are pinned in `tests/test_rules.py`.

---

## D-004 - Age unit conversion uses the mean Gregorian year

**Decision.** `days/7`, `weeks*1`, `months*(365.25/12/7)`, `years*(365.25/7)`.

**Why.** "4 weeks per month" would place 6 months at 24 wk (young_adult) when
it is 26.1 wk, still young_adult but within 1 wk of the middle_aged edge; the
error is not negligible at band edges. A single documented factor is better
than two approximations.

---

## D-005 - The `confidence` mapping is this project's choice, not SCHEMA.md's

**Decision.** With `S = has_actuator_control + n_reported + effect_size_reported`
(strict `True` only):
`low` if not `full_text_available` **or** `S <= 1`;
`high` if `full_text_available` and `S == 3`;
`moderate` otherwise.

**Why.** SCHEMA.md 3.6 names the four inputs and requires a documented rule
but gives no mapping. Without full text a record cannot be re-checked at all,
so it cannot be better than `low` regardless of the other three. All 16 rows
are pinned in `tests/test_rules.py`.

**Defect in SCHEMA.md.** 3.6 should contain this table. Until it does, this
entry is the specification.

---

## D-006 - Deviation: controlled vocabularies are JSON, not YAML

**Decision.** Vocabularies live in `data/vocab/populations.json`,
`outcomes.json`, `methods.json`, `regions.json`.

**Why.** SCHEMA.md 3.2/3.5 refer to `vocab/populations.yaml` and
`vocab/outcomes.yaml`. Parsing YAML requires `pyyaml`, which contradicts
D-001. JSON is in the standard library and loses nothing the vocabularies
need.

**Status.** Deviation from the binding document, accepted knowingly. SCHEMA.md
should be amended to say JSON at the next version bump.

---

## D-007 - Added field: `assay_units`

**Decision.** `TierARecord` carries an optional `assay_units` alongside
`assay`.

**Why.** The QC layer is required to check measured units against the `units`
field of the `outcomes` vocabulary entry, and SCHEMA.md 3.5 provides nowhere
to record the unit a value was measured in: `effect_size_units` describes the
*effect*, which for a percent change or Cohen's d is not the readout's unit.
The field is optional and additive; existing records without it stay valid.

**Status.** Extension beyond SCHEMA.md v0.1.0; needs ratifying at the next
version bump.

---

## D-008 - Enumerated time units; `title`/`journal`/`year` treated as nullable

**Decision.** `duration_unit` and `time_to_readout_unit` use a `TimeUnit` enum
(`seconds`..`months`). `title`, `journal` and `year` are nullable, with a
WARNING (`W012_MISSING_BIBLIOGRAPHY`) when absent.

**Why.** SCHEMA.md 3.4 leaves the duration units unenumerated, which would
admit "wks", "weeks" and "w" in one column. SCHEMA.md 3.1 does not mark the
bibliographic fields optional, but harvest metadata is routinely incomplete,
and rule 4 ("if a field is not reported, it is null") outranks a table's
formatting. Making them hard requirements would reject records whose evidence
is intact.

---

## D-009 - Lab identity is a WARNING, not a hard failure

**Decision.** Missing `last_author` / `last_author_affiliation` emits
`W010_NO_LAB_IDENTITY`, not a hard failure.

> **Amended by D-056 and D-057.** Lab identity is the `last_author` alone, so
> `W010` now fires only when the author is missing and the new
> `W018_NO_LAB_AFFILIATION` covers a missing affiliation. The severity decision
> below is unchanged; only which field triggers which warning.

**Why.** SCHEMA.md 3.1 types them `str?` yet annotates them "required for
lab-level leakage control" -- an internal inconsistency. Blocking admission
would discard usable evidence; ignoring it would silently break the leakage
control that SCHEMA.md 5 calls the project's largest risk. A counted warning
plus a report statistic keeps the problem visible. Any split that cannot group
a record's lab must say so, per SCHEMA.md 5.

---

## D-010 - Unresolved populations group by raw name, under a `raw:` prefix

**Decision.** QC grouping keys use `population_id` when present, else
`"raw:" + normalised population_raw`.

**Why.** SCHEMA.md 3.2 keeps unmapped populations with `population_id = null`
and forbids force-mapping them. Grouping all nulls together would invent
claims ("these two labs disagree") between unrelated populations; skipping
them entirely would hide real duplicates. The prefix guarantees an unresolved
population never groups with a resolved one.

---

## D-011 - Only two direction pairs count as a conflict

**Decision.** `increase` vs `decrease` -> `QC_CONFLICT_OPPOSED`.
`increase`/`decrease` vs `no_change` -> `QC_CONFLICT_EFFECT_VS_NULL`.
`not_reported` conflicts with nothing. `non_monotonic` is not treated as
contradicting a monotonic direction.

**Why.** A non-monotonic response can present as either arm depending on the
dose or time window sampled, so calling that a contradiction would be an
interpretation rather than an observation. Conflicts are recorded on both
sides in `conflicts_with` and reported at `INFO`: they are properties of the
literature, not defects, and nothing in this codebase resolves, prefers,
averages or drops either side.

---

## D-012 - Hard vs warning, and what the CLI's exit code means

**Decision.** Hard findings: the seven SCHEMA.md 2 rules (`A*`), structural
requirements (`S*`) and value sanity (`V*`). Warnings: `W*` plus QC
duplicates, unit mismatches and missing conditions. Conflicts and possible
non-independence are `INFO`. `python3 -m syn4rch.qc.run` exits non-zero only
on hard schema violations or an unreadable line.

**Why.** An exit code that fires on disagreement in the literature would
train whoever runs it to suppress disagreement.

---

## D-013 - Deferred, knowingly

* **Non-mouse `age_group`.** Rats and humans return `not_reported`; the mouse
  boundaries do not transfer and no species-specific rule has been authored.
* **Shared-cohort detection** (SCHEMA.md 5 item 3) is manual; nothing detects
  two papers sharing one animal cohort. Non-independence detection only
  catches same-lab repeats of the same claim.
* **Population x outcome leave-one-out grouping** (SCHEMA.md 5 item 4) has no
  code yet; the keys it needs exist in `src/syn4rch/qc/keys.py`.
* **Numeric unit conversion.** Unit checking is textual, with a small alias
  table. It never converts a value and never rewrites a record.
* **Semantic verification of `support_quote`.** Code can check that a quote
  exists and contains a digit when a number is claimed; it cannot check that
  the quote is verbatim from the source, nor that the stated `direction`
  really follows from it. That remains a human obligation, and
  `extraction_method` is the only honest record of whether it was met.
* **`significant` threshold checking** is a warning only; whether the paper
  actually stated a threshold comparison cannot be verified from the record.

---

## D-014 - Unit matching is textual, with optional parenthetical qualifiers

**Decision.** `outcomes` vocabulary entries annotate units in prose (e.g.
`"degC (BAT temperature)"`, `"Ucp1 mRNA fold change"`). When comparing a
record's units against them, the parenthetical qualifier is optional: both the
full string and the part before `(` are accepted. A small alias table folds
obvious spelling variants (`mg/dL` ~ `mg per dL`, `percent` ~ `%`, `degC` ~
`°C`).

**Why.** Without this, nearly every annotated vocabulary unit would raise a
spurious `QC_UNIT_MISMATCH`, and a check that is always wrong gets ignored.

**Limits, knowingly accepted.** The comparison remains purely textual. It
cannot tell `mg/dL` from `mmol/L` (a real dimensional difference needing a
molar mass), it will not catch a value recorded in the wrong one of two
declared units, and it never converts or rewrites anything. Vocabulary units
that are free prose (`"UCP1 protein a.u."`) are effectively unverifiable.
A dimensional unit system is out of scope for v0.1.0.

---

## D-015 - Verbatim `support_quote` verification is now mechanised

**Decision.** `src/syn4rch/curate/verify.py` checks that `support_quote`
occurs as a substring of the paper's retrieved full text, under one
normalisation: whitespace runs collapse to a single space, and unicode dashes,
quotation marks and non-breaking spaces fold to ASCII. No case folding, no
stemming, no fuzzy matching. A proposal that fails is rejected and logged, and
is never repaired.

**Supersedes.** The bullet "Semantic verification of `support_quote`" in
D-013, but only in part. D-013 said code "cannot check that the quote is
verbatim from the source". That was true when no source text was held locally.
It is no longer true: stage 3 of `docs/CURATION_PROTOCOL.md` writes
open-access full text to `data/raw/fulltext/`, so the source string is on disk
and the check is a normalised substring test.

**Still deferred, unchanged from D-013.** Code cannot check that the stated
`direction` follows from the quote, nor that the quote is representative of
the paper rather than the one sentence that suits a conclusion. Those remain
reading judgements, and `extraction_method` plus `review_status` are the only
honest record of whether a human made them.

**Scope limit.** The check applies only to papers whose full text was
retrieved. Abstract-only candidates cannot yield records at all (SCHEMA.md
rule 2.5), so there is no unverifiable remainder — a record either has a
locally held source to check against or it does not exist.

---

## D-016 - Markdown reports are version controlled

**Decision.** `.gitignore` tracks `reports/*.md` while continuing to ignore
everything else written to `reports/`.

**Why.** `data/raw/` and `data/interim/` are deliberately untracked (they are
large, and reproducible from the harvest). A reader who clones this repository
therefore holds the curated dataset and the code, but not the harvest. If the
reports were also untracked, every measured claim this project makes about its
own coverage would be absent from the clone and unreproducible by the reader.
Tracking the reports keeps the claim and its provenance together.

**Cost, knowingly accepted.** A committed report can fall out of date with
respect to a regenerated dataset. Mitigation: every report states the run
timestamp and the inputs it was computed from, and is regenerable by a
committed CLI. A report that cannot be regenerated is a defect.

---

## D-017 - Records are proposed and then gated, never written directly

**Decision.** Nothing writes to `data/curated/records.jsonl` except the gate
in `src/syn4rch/curate/run.py`. Extraction produces *proposals* in
`data/interim/`; admission is a separate, deterministic step that runs the
verbatim check (D-015), then `syn4rch.schema.validate`, then QC.

**Why.** The stage that reads papers is the only stage capable of fabricating
a record, and it is the stage that cannot be made deterministic. Separating
proposal from admission means the fallible stage never holds the pen on the
dataset. Collapsing the two — letting an extractor append directly to the
curated file — would leave no point at which a fabricated record could be
caught, and the anti-fabrication rules of SCHEMA.md 2 would be advisory.

**Consequence.** The rejection report is part of the dataset's documentation,
not a debugging artefact. The ratio of proposals to admitted records is a
published statistic.

---

## D-018 - A quote must sit in the passage its locator names

**Decision.** The curation gate rejects a record whose `support_quote` is
verbatim in the paper but not in the passage named by `source_locator`
(`W101_QUOTE_NOT_IN_LOCATED_PASSAGE` blocks by default). `--lenient-passage`
downgrades it to a note, for hand-curated records only.

**Why.** CURATION_PROTOCOL.md 3.3 requires the locator to be the extraction
unit's locator, unmodified, and units are built from the same passage
segmentation the gate uses, so in the normal pipeline the two can never
disagree by accident. When they do disagree, the likeliest cause is a reader
pairing a direction from one passage with a quote from another, which is the
fabrication pattern the gate exists to catch. The lenient mode exists because
a human counting paragraphs in a PDF may number them differently from the
parser; that is a reason to allow an override, not to make the override the
default.

**Also recorded, not yet acted on** (raised while building the gate):
supplementary-material locators (e.g. "Fig S1A") can never be admitted,
because supplements are not in the retrieved XML; and `record_id` (D-002)
does not include the quote, so two proposals differing only in quote share an
id in the rejection report.

---

## D-019 - `direction` describes the outcome, not the assay number

**Decision.** `direction` must describe the quantity named by `outcome_id`,
never the raw assay number, whenever the two have opposite polarity. Two steps,
in order:

1. Prefer the `outcome_id` the passage literally measured, when the vocabulary
   has one. "A significant decrease in blood glucose levels in the ITT" is a
   `blood_glucose` measurement: record `blood_glucose` / `decrease`.
2. Only when the vocabulary has no entry for the measured quantity, use the
   interpreted outcome and invert the direction to match it. "Ameliorated
   insulin resistance" under `insulin_sensitivity` is `increase`. "BAT lipid
   droplet size decreased" under `bat_thermogenesis` is `increase`. The
   inversion and the measured quantity go in `extraction_notes`.

If the inversion cannot be read from the `support_quote` alone, there is no
record. A clean result is not discarded merely because an inversion is involved.

**Why.** Found in the first pilot batch, where "a significant decrease in blood
glucose levels in the ITT" had been recorded as `insulin_sensitivity` /
`decrease`. Falling glucose in an insulin tolerance test is *rising* insulin
sensitivity, so the sign was inverted. Three records in one batch of 41 were
affected. A second batch met the same situation and discarded the result
instead, so two curators produced opposite treatments of one case: worse than
either choice alone, because the inconsistency is invisible once the records
are pooled.

A sign error of this kind is the most damaging defect this dataset can carry.
It does not look like missing data or a failed check. It trains a model
backwards, and every downstream accuracy figure would be computed against
corrupted labels while looking entirely healthy.

**The vocabulary invites the error.** `insulin_sensitivity` lists `ITT` among
its synonyms, `glucose_tolerance` lists `GTT`, and `bat_thermogenesis` admits
lipid-droplet measures. In all three the assay's number moves opposite to the
named quantity. That a synonym matches the words in the text does not make it
the correct `outcome_id`.

**Not yet done.** Nothing mechanically detects a polarity inversion; the check
is a reading judgement, like `direction` itself. Marking the hazardous
vocabulary entries with an explicit polarity flag, and warning in QC when such
an entry is used, is the obvious next step and is not implemented.

---

## D-020 - One passage may yield several records

**Decision.** A passage yields every record it supports, not at most one.

**Why.** `docs/CURATION_PROTOCOL.md` 3 originally said the reader "proposes at
most one record" from a unit, which contradicts SCHEMA.md 1: one record is one
(population x perturbation x readout x study) observation, and aggregating
across outcomes is forbidden. A passage reporting food intake, body weight and
energy expenditure for one manipulation describes three observations. Forcing
them into one record would either drop two or aggregate across readouts.

**Resolution.** SCHEMA.md wins, per the protocol's own preamble, and the
protocol's wording has been corrected. What stays forbidden is the reverse:
splitting one observation into two records that restate it, and merging two
passages into one record.

---

## D-021 - Polarity hazards are annotated in the outcomes vocabulary

**Decision.** `data/vocab/outcomes.json` carries two new per-entry fields and a
`_meta.polarity_policy`. Vocabulary version 0.1.1.

* `inverted_synonyms` — the quantity these words name moves **opposite** to the
  quantity the entry is named for. `insulin_sensitivity` lists "insulin
  resistance" and "HOMA-IR"; `glucose_tolerance` lists "glucose excursion";
  `sleep` lists "sleep latency".
* `signed_synonyms` — the words already carry a direction or name a clinical
  state, so `direction` must describe the *change* in the entry's quantity, not
  the presence of the state. "Reduced hyperglycemia" is `blood_glucose`
  **decrease**. 17 such synonyms across 10 entries.

`insulin_sensitivity` and `glucose_tolerance` additionally carry a note sending
raw ITT/GTT glucose values to `blood_glucose` (D-019 step 1).

**Why.** D-019 recorded that nothing made the polarity trap visible. It was
invisible in exactly the wrong place: at the moment a curator matches the words
in a paper against a synonym list. Four of five pilot batches hit it, and two
reached opposite conclusions on the same case. Annotating the entry puts the
warning where the mistake is made.

**`higher_is` is a different field and was not reused.** It encodes whether a
higher value is better or worse *for the organism* and is deliberately null for
39 of 43 entries. Polarity is about whether the assay number tracks the named
quantity, which is a measurement question, not a value judgement. Conflating
them would have made both useless.

**Deliberately not automated.** The policy states that these lists must never
be used to rewrite a direction automatically, and that they are not exhaustive.
They record the traps found so far. A QC warning when `outcome_raw` matches an
annotated synonym would be useful and is not implemented.

**Found a second defect, not fixed.** The id `cognitive_decline` and its label
"Cognitive performance / decline" name both a performance and a decline, so
`direction` is ambiguous on its face. Its note now rules that direction
describes cognitive *performance*, and records that the id should be renamed at
the next version bump.

---

## D-022 - Rulings the first pilot forced, collected

Five curators hit the same ambiguities independently. Settled as follows; none
required rework, so these record what was accepted.

**Region-level manipulations are kept, with `population_id = null`.** An
intra-nucleus injection of a ubiquitous-promoter vector, or an i.c.v. drug, is
not a defined neuronal population, which SCHEMA.md 1 makes the unit of record.
They are nonetheless causal interventional evidence, 29 records across two
batches. They are admitted with `population_raw` as the paper names it and
`population_id` null, which is sufficient because `qc.keys.population_key`
(D-010) prefixes an unresolved population `raw:`, and `models.task` and
`feasibility.splits` both use that function. A region-level record can
therefore never pool into a cell-type-specific cell, nor be used as evidence
for one. A `raw:` group is not a cell-type-specific claim and must not be read
as one.

**A condition qualifier may go in `outcome_raw`.** `record_id` (D-002) hashes no
diet, sex, timepoint or fed/fasted field, so two genuine observations differing
only in condition collide. Curators disambiguated with the paper's own wording,
e.g. "food intake (chow diet)", "weight gain, female cohort", 13 records in one
batch alone. Accepted: the alternative is dropping one of the two, and the
collisions fall disproportionately on null results, which rule 6 forbids
dropping. A condition field belongs in the next schema version.

**`significant` is set from a bare "significantly".** The paper asserting
significance *is* the threshold comparison SCHEMA.md 3.5 asks for. That it does
not name the test is already reported as `W008`, so the information is not lost.

**`has_actuator_control` is null, not false, when the passage does not mention
controls.** `false` would assert the paper has no such control, which the
passage does not support.

**`is_primary_outcome` is left null.** No passage states it; filling it would be
curator impression.

**Co-perturbation and dual-population designs yield no record.** No field holds
a second simultaneous manipulation, so a direction would silently attach to the
wrong cause. Discarded with the reason in `extraction_notes`, not dropped
silently.

**A reported non-significant trend is `no_change` with `significant = false`**,
and the trend quoted in `extraction_notes`. `decrease` is unsupported and
dropping it is forbidden.

**`non_monotonic` covers a biphasic response over dose or over time.** The enum
named neither.

---

## D-023 - The claim key carries the condition; unknown pools with unknown

**Decision.** Schema v0.2.0 adds a condition dimension
`(sex, diet, age_group, challenge)` to the schema and to
`syn4rch.qc.keys.claim_key`. Two records whose condition dimension differs are
not comparable, so they can neither conflict with nor replicate each other.

Within the key, a condition that is `null` and one that is `not_reported` are
pooled into a single `unknown` token. Both are preserved as stored values,
because they are different statements about the *curation*: `null` means the
curator did not address the question, `not_reported` means the paper was
checked and does not say. They are the same statement about the *experiment* —
the condition could be anything — and the key answers the experimental
question.

**Why.** D-022 recorded that curators were disambiguating condition collisions
by hand, inside `outcome_raw`, with the paper's own wording. That was the right
local decision and the wrong structural one. Measured on the first 362 records:

| claim key | cells | conflicting claims | directly opposed | multi-study claims |
|---|---|---|---|---|
| v0.1.0 `(population, valence, outcome_id)` | 127 | 36 | 12 | 39 |
| v0.2.0 (this decision + D-025 + D-026) | 249 | 9 | 1 | 11 |

Eleven of twelve head-to-head `increase` vs `decrease` contradictions in this
dataset were never contradictions, and 27 of 36 conflicting claims of any kind
were not conflicts.

**Revised before commit, and the revision matters.** An earlier draft of this
entry reported 229 cells, 20 conflicts and multi-study cells falling 39 to 29,
and concluded that "most replication evidence survives the change". That was
measured before `target_raw` and `target_id` were resolved (D-026), while every
gene knockout still pooled with every other. With the targets resolved the
figures above are the real ones, and the conclusion inverts: **multi-study
claims fall from 39 to 11, so roughly three quarters of the apparent cross-study
replication in this corpus does not survive the change.**

That is the finding, not a side effect. Under a key that pools a male cohort
with a female one and one gene with another, this dataset looked like it held 39
independently replicated claims. It holds 11. The replication was an artifact of
the same coarseness that manufactured the conflicts, and both numbers moved for
one reason.

**Why not keep every unknown apart.** It would mean two records that each fail
to state a diet are held incomparable, making both conflict *and* replication
undetectable across the 197 of 362 records that carry no diet. That discards
real information: two labs independently reporting that AgRP activation raises
food intake is replication, whether or not either states the diet.

**The cost, stated.** A cell matched on "both unknown" is weaker evidence than
one matched on a known sex and diet. `condition_specificity` (0..4) is reported
with every conflict and every replication figure so the two cannot be printed
as one number. Across the corpus: 153 records know no condition slot, 161 know
one, 43 know two, 5 know three, and **none knows all four.**

**Pinned.** `tests/test_qc_keys.py`.

---

## D-024 - A condition may be read from Methods, with a locator

**Decision.** `condition_source` ∈ {`passage`, `methods`, `not_reported`} and
`condition_locator`. When `condition_source = methods`, `condition_locator` is
mandatory and its absence is a HARD finding
(`A008_CONDITION_FROM_METHODS_WITHOUT_LOCATOR`).

**Why.** SCHEMA.md 2 rule 4 requires a field the record's passage does not
report to be `null`, and papers put the result in Results and the diet in
Methods. Obeying rule 4 literally therefore guaranteed the condition gap: 197
of 362 records with no diet, 224 with no sex. No amount of curating harder
fixes that, because the curators were right every time.

The conditions are a property of the experiment, not of the sentence. Reading
them from Methods is legitimate. Reading them from Methods *without saying so*
is not, because it would make rule 4 unenforceable everywhere else.

**What this does not permit.** Inference. A condition no passage of the paper
states is still `null`. `condition_source` widens where a condition may be read
*from*; it does not widen what may be written down. The verbatim gate applies
to `condition_locator` exactly as it does to `source_locator`.

**Not yet implemented.** The verbatim gate does not currently re-check
`condition_locator` (`syn4rch.curate.verify` reads `source_locator` only), so
today the obligation is documented and validated for presence but not for
substance. Recorded here as a gap rather than left to be discovered.

---

## D-025 - `perturbation_class` and `perturbation_polarity` split the axis

**Decision.** Two derived fields, both functions of `method` (and, for two
classes, of `valence`):

* `perturbation_class` ∈ {`circuit_activity`, `molecular`, `structural`,
  `pharmacological`, `other`}. A total function of `method`
  (`syn4rch.schema.rules.METHOD_CLASS`), admitting **no** override: a stored
  value disagreeing with the rule is a HARD finding
  (`S007_PERTURBATION_CLASS_MISMATCH`).
* `perturbation_polarity` ∈ {`increase`, `decrease`, `remove`,
  `bidirectional`, `unclear`} — the sign of the manipulation along its own
  class's axis. A stored value disagreeing with the rule is a WARNING
  (`W013`), not a refusal, because for a locally infused drug a curator who has
  read the pharmacology may know the sign when `valence` does not.

`valence` keeps its v0.1.0 meaning untouched and `unclear` remains a correct
answer that must not be coerced.

**Why.** `valence` was being asked to answer two questions at once. For a
conditional knockout of a receptor gene, "did firing go up or down" honestly
has no answer, so `valence` was `unclear` — correctly — and the record had no
signed axis. That accounted for 185 of 362 records, 170 of them molecular, so
roughly half the dataset could not be a model input. The manipulation itself,
meanwhile, has a perfectly well defined sign: a knockout removes a gene
product and an overexpression adds one.

**Polarity is meaningless without its class.** `decrease` is less firing under
`circuit_activity` and less gene product under `molecular`. Pooling those is
the D-019 error in a new place. Enforced structurally rather than by
convention: `syn4rch.qc.keys.perturbation_key` returns
`(class, polarity, target)` as one value, and nothing in this project keys on
polarity alone. `tests/test_qc_keys.py` pins it.

**Honest about what it buys.** This makes half the corpus *representable*. It
does not make it *predictable*, and D-027 records why not.

---

## D-026 - `record_id` gains the target and the condition (supersedes part of D-002)

**Decision.** `RECORD_ID_FIELDS` grows from 7 components to 12, adding
`target_raw`, `sex`, `diet`, `age_group`, `challenge`. Every `record_id` in the
corpus changes. `syn4rch.curate.migrate` carries the old-to-new map and rewrites
`conflicts_with` through it.

**Why.** D-022 reported three `record_id` collisions in curation batch 04. The
cause was structural: a male cohort and a female cohort, reported in one passage
for one readout, were two observations that hashed to one identity. The
condition-qualifier-in-`outcome_raw` convention worked around it; this fixes it.

**Free text in the hash, against D-002's own rule.** D-002 excluded free-text
fields because re-wording one during curation would silently re-identify the
record. `target_raw` is admitted anyway: a knockout of *Raptor* and a knockout
of *p53* in the same neurons, same readout, same locator genuinely collide
otherwise, and a certain collision is worse than a hypothetical
re-identification. A gene symbol is also about as stable as free text gets.
`target_id` is *not* hashed, because it is the resolved form of a component
already in the hash.

**Still excluded, deliberately.** `targeting_strategy` and every other prose
field. `perturbation_class` and `perturbation_polarity`, which are functions of
`method` and `valence` and so add no discrimination. `time_to_readout_value`,
because a float in a hash is a portability hazard and the D-022 qualifier
convention already separates timepoints through `outcome_raw`.

**Pinned.** `tests/test_record_id.py` holds the new golden digest, asserts that
the function signature matches `RECORD_ID_FIELDS`, and asserts that the batch-04
collision class is closed.

---

## D-027 - Perturbation identity, and why v0.2.0 does not make direction predictable

**Decision.** A cross-validated number must be reported alongside the fraction
of test rows whose perturbation identity (`perturbation_key`) never appears in
the training side of the same split. `syn4rch.feasibility.identity` computes it.
A high unseen fraction is not a modelling deficiency to be tuned away; it is
reported as a property of the corpus.

**Why.** All 170 molecular records resolved a target (D-026). They cover **17
distinct target genes across 16 studies and 16 laboratories.** No gene is
manipulated by more than one study or more than one laboratory; one study carries
two genes (`il17ra` and `socs3`). Under `by_study` or `by_lab` grouping a
held-out molecular record's perturbation identity is therefore *always* unseen in
training. There is nothing for a model to have learned, and an accuracy figure
averaged over those folds measures the label prior.

**The axis looks sound, and the look cannot be corroborated.** Pooled across
genes, the 17 `(population, outcome)` cells holding both a loss-of-function and
a gain-of-function record mostly do not mirror — `arc_pomc / hepatic_lipid`
reads loss {increase 1, decrease 1} against gain {decrease 1, increase 1} —
which looks like noise until the target is in the key. Resolve the gene and
there are 9 `(population, outcome, target)` cells where the *same* gene was both
knocked out and overexpressed. Six mirror cleanly (`bcl6` on food intake, BAT
thermogenesis and hepatic lipid; `fam172a` on food intake and energy
expenditure; `trim31` on hepatic lipid), one agrees in both arms as a null, two
have no clear modal direction, and **none contradicts.**

**Revised before commit.** An earlier draft read that as "`perturbation_polarity`
is measuring something real". It does not establish that. **All 9 of those cells
are a single laboratory and a single study**, verified twice: every mirrored pair
is one paper reporting both arms of its own experiment. That is exactly as
consistent with a real biological sign as it is with a paper presenting a
coherent story, and no analysis available here separates the two. Zero
gene-level bidirectional comparisons in this corpus span two laboratories.

So the honest statement is narrower: the molecular axis is *internally
consistent*, which is a necessary condition for it to be measuring anything, and
its identity never repeats across laboratories, so neither its validity nor its
transferability can be tested. Nothing learned about one gene transfers to
another, and nothing in the corpus checks whether a gene's reported sign would
replicate.

**What follows.** The schema fix was still worth making: it closed an id
collision class, removed eleven of twelve false contradictions, and made the
corpus's real structure visible. But it does not unlock prediction for the
molecular half, and curating more molecular papers will not either. What would:
a second laboratory manipulating a gene a first laboratory already manipulated.
Until then the honest answer for that half of the corpus is the refusal
SCHEMA.md 6 provides for — insufficient evidence to rank these interventions
reliably — and the circuit-activity half is where a model has any chance.

**Standing instruction.** Do not present v0.2.0 as a predictive improvement
unless a benchmark shows one. A schema fix that does not help is a result.

---

## D-028 - Direction is predictable in the circuit half, and v0.2.0 is why we can see it

**The measurement.** Leave-one-study-out, frequency majority with a
global-majority backoff, deterministic tie-break in the order
(increase, decrease, no_change). 361 labellable records, 42 studies.

**Two implementations agree.** The figures below were first measured by a
script sharing no code with `syn4rch.models` beyond `syn4rch.qc.keys`,
specifically so that the committed benchmark had something to be checked
against rather than only agreeing with itself. `syn4rch.models` then reproduced
the circuit-half figures exactly: 116 rows, 20 study groups, 0.397 against
0.517, +0.121, bootstrap interval excluding zero, McNemar p=0.156 on b=49 and
c=35 discordant rows.

They differed on the molecular half, and the committed implementation was
right. See the note under the table.

Pooled over the whole corpus — which is what every earlier benchmark did:

| baseline | accuracy | on rows whose cell was seen |
|---|---|---|
| global majority | 0.446 | 0.446 |
| v0.1.0 key `(population, valence, outcome)` | 0.432 | 0.444 (n=241) |
| v0.2.0 key `(population, class, polarity, outcome)` | 0.457 | 0.465 (n=215) |

Nothing. Now split by `perturbation_class`:

| subset | rows | studies | global majority | v0.2.0 key | difference |
|---|---|---|---|---|---|
| `circuit_activity`, cell seen | 116 | 20 | 0.397 | **0.517** | **+0.121** |
| `molecular`, cell seen, target **in** the key | **0** | — | — | — | nothing to measure |
| `molecular`, cell seen, target **omitted** | 99 | 15 | 0.505 | 0.404 | −0.101 |

Study-clustered bootstrap over the studies: the committed benchmark
(`reports/baseline_benchmark.md`, 1000 replicates) gives 95% interval
[+0.008, +0.271]; the independent script (4000 replicates) gives median +0.123,
interval [+0.011, +0.256], with 1.7% of resamples showing no gain. The
interval's lower bound sits a hair above zero either way, which is the honest
shape of this result.

**The molecular half is a structural zero, not a negative result.** With the
manipulated gene in the cell key — the only honest key for a molecular
perturbation, and the one v0.2.0 made expressible — **not one of the 170
molecular rows has its cell populated in the training side of any study-grouped
fold**, because each gene appears in exactly one study. The predictor is the
global prior for every one of them, so there is no difference to measure.

The −0.101 on the second row is what you get by *omitting* the target, which
pools `trp53` with `rptor` with `sirt1` and asks the model to learn from the
pool. It is worse than guessing because the pool is noise. An earlier draft of
this entry reported that figure as the molecular result; it is better read as
the cost of the pooling the target field removes.

Either way the two halves cancel in the pooled number, which is the whole
reason the corpus looked unpredictable.

**Why rows whose cell was unseen are excluded from the comparison.** On such a
row both baselines predict the global majority, so their predictions are
identical by construction and the row cannot contribute to the difference.
Including them dilutes the effect and smuggles in a composition artifact: the 53
unseen circuit rows happen to score 0.528 under both baselines. Both figures are
reported above; the "cell seen" column is the one that measures learning.

**What v0.2.0 did and did not do.** It did **not** create the signal. On the
circuit subset the v0.1.0 key and the v0.2.0 key score *identically* (0.517 on
the same 116 rows), because for a circuit manipulation `perturbation_polarity`
is a restatement of `valence`, which v0.1.0 already had.

What v0.2.0 created is the ability to say which records belong in that subset.
All 170 molecular records carried `valence: unclear`, and so did 10 circuit
records whose net effect the paper did not resolve. No field distinguished "a
gene was manipulated" from "the net circuit effect is unresolved", so every
analysis pooled them, and the molecular half's anti-signal cancelled the circuit
half's signal. `perturbation_class` is what makes the split expressible.

So the earlier conclusion — that direction of effect is not predictable from
this corpus — was an artifact of a missing field, for half the corpus. It
remains correct for the other half, and D-027 explains why: 17 genes, 16
studies, 16 labs, and no gene repeated across laboratories.

**What may and may not be claimed.** The defensible claim is that this is the
first signal above baseline in this project and that it is confined to the
circuit-activity half. Not "direction is predictable":

* 116 rows, 20 studies, and the bootstrap interval's lower bound is +0.011.
* McNemar on the same rows gives b=49, c=35, two-sided exact p=0.156, so the
  row-level paired test does **not** reach significance. Two tests pointing
  different ways is the honest summary: suggestive, not established.
* 0.52 on a three-way problem is useful for ranking hypotheses and nowhere near
  a simulator.
* A third of circuit rows (53 of 169) still have no evidence in their cell at
  all, and for those the system has nothing to say but the prior.
* **The circuit half is learnable partly because its perturbation vocabulary is
  coarse.** `syn4rch.feasibility.identity` measures 3 distinct perturbation
  identities across all 170 circuit records (increase, decrease, unclear) and 22
  across the 170 molecular ones. A chemogenetic activation is "the same
  manipulation" in every laboratory only because the identity does not record
  which receptor, which promoter, which dose or which titre. So the 0% unseen
  fraction for `circuit_activity` is a statement about the resolution of the
  key as much as about the biology, and the modest 0.517 accuracy is what that
  coarseness buys. Sharpening the circuit identity — an actuator or dose
  component — would raise the unseen fraction and very likely lower the
  accuracy. That is a trade to make deliberately, with the numbers in view, not
  by accident.

**Consequences for what to curate next.** Curating more molecular papers will
not help (D-027). Curating more optogenetic and chemogenetic papers plausibly
will, because their perturbation identity — a population and an activation or
inhibition — recurs across laboratories by construction. That is a change of
plan from "curate the remaining 924 full texts", and it is the first time this
project has had a measured reason to prefer one paper over another.

---

## D-029 - `cognitive_decline` renamed to `cognitive_performance`

**Decision.** The outcome id `cognitive_decline` becomes
`cognitive_performance` (outcomes vocabulary 0.2.0). `direction` on it
describes the change in cognitive **performance**: `increase` = better. Its
two annotated synonyms move from `signed_synonyms` to `inverted_synonyms`, and
`escape latency` is added to both the synonym list and the inverted list.

**Why.** D-021 found this defect and deferred it: the id and its label named
*both* a performance and a decline, so `direction` was ambiguous before a
curator had read anything. "Cognitive decline decreased" and "cognitive
performance decreased" are opposite claims that the old id could not tell
apart. That is the D-019 failure mode built into a vocabulary entry rather
than stumbled into by a reader.

The reclassification matters as much as the rename. Under the old id, "cognitive
decline" was *signed* — a word carrying its own direction. Under an id that
names performance only, the same word is *inverted*: a paper reporting less
decline is reporting more performance, and a reader who copies the direction
across will record it backwards. `escape latency` is the same trap in assay
form, since a shorter latency is better memory.

**Cost: none.** No curated record used the id — the 362 records cover 25
outcomes and this was not among them — so the rename orphans nothing. Three
references existed in code and tests and all three were updated. Doing it now
is the point: D-021 said "at the next version bump", and this is that bump.
Deferring it again would have meant carrying a known polarity trap into
whatever gets curated next.

**`higher_is` stays null**, and not for want of an opinion. The measures are
not comparable across assays, and this is a Tier B outcome, which SCHEMA.md 6
forbids as a training target in any case.

**Pinned.** `tests/test_vocab_polarity.py` asserts the old id is gone, the new
one resolves every old wording, and the changelog records the rename.

---

## D-030 - `qc.conditions` renamed to `qc.comparability`

**Decision.** The module `syn4rch/qc/conditions.py` becomes
`qc/comparability.py`, and with it:

| was | is |
|---|---|
| `REQUIRED_CONDITIONS` | `REQUIRED_COMPARABILITY_FIELDS` |
| `missing_conditions` | `missing_comparability_fields` |
| `find_missing_conditions` | `find_missing_comparability_fields` |
| `condition_coverage` | `comparability_coverage` |
| `QC_MISSING_CONDITIONS` | `QC_MISSING_COMPARABILITY_FIELDS` |

No aliases are kept. Nothing outside `syn4rch.qc` imported these, and a
deprecation shim would preserve the ambiguity this decision exists to remove.

**Why.** The module checked four methodological fields — a control, a group
size, a species, a perturbation valence — and called them "conditions", stating
that they "are needed before a record can be compared with another record".
D-023 then introduced the *experimental condition dimension*,
`(sex, diet, age_group, challenge)`, which is the thing that actually decides
comparability. Two records agreeing on all four fields in this module are still
not comparable if one is a fasted male on a high-fat diet and the other a fed
female on chow.

So after v0.2.0 the codebase had two unrelated meanings of "condition", one of
them attached to a false claim about comparability. This project has twice had
to undo a name that meant two things at once — D-019 on `direction`, D-029 on
`cognitive_decline` — and both times the cost was a record filed backwards.
Introducing a third was not acceptable, particularly since this one was
introduced by the same change that created the ambiguity.

**The QC code changed too**, deliberately. `QC_MISSING_CONDITIONS` fired 359
times on the first 362 records and was read — including by me — as evidence of
the condition gap that D-023 addresses. It was never measuring that. A stable
code is worth keeping only while it means what its name says.

**The finding message now says what it is not**, naming the condition dimension
explicitly, so the same misreading cannot be made from the report alone.

---

## D-031 - The migration is a second writer to `data/curated/` (amends D-017)

**Decision.** D-017 says records are proposed and then gated, and that the gate
at `syn4rch.curate.run` is the only thing that writes
`data/curated/records.jsonl`. `syn4rch.curate.migrate` is a second writer. The
exception is granted, narrowly, on these conditions, all of which it meets:

1. It **creates no record, admits no record and drops no record.** Count in
   equals count out, asserted in code, and a run that cannot satisfy that writes
   nothing.
2. It **fills no field a published rule does not derive.** `perturbation_class`
   and `perturbation_polarity` come from `syn4rch.schema.rules`; `target_raw` and
   `target_id` are resolved from the already-curated `targeting_strategy` and
   `actuator` against a committed vocabulary; nothing else is written except
   `schema_version` and the recomputed `record_id`.
3. It **refuses to write if validation got worse.** The hard-finding count per
   code must not exceed its pre-migration baseline, with `S005` and `S006`
   discounted because those two *are* the version change.
4. It is **reversible and re-checkable.** The pre-migration file is preserved at
   `data/curated/records.0.1.0.jsonl.bak`, the old-to-new id map is committed,
   and `migrate report` re-runs the whole migration in memory from the backup and
   asserts that it reproduces the committed dataset byte for byte.
5. It is **idempotent.** A second run exits non-zero and leaves the file
   unchanged rather than migrating twice.

**Why an exception at all.** A version bump that changes every `record_id`
cannot go through the gate: the gate reads *proposals*, and these are already
admitted records. The alternative was re-curating 362 records from their
passages, which would have re-run the only fallible stage in the pipeline for no
reason and risked changing data that is correct.

**What is not granted.** No future writer gets this by precedent. Each one needs
its own entry and must meet the five conditions above or better. In particular,
the pass that fills the condition dimension will re-identify every record it
touches, because `challenge` is a `record_id` component (D-026) and is null on
all 362 records today. That pass is a *curation* pass, not a migration: it reads
papers, so it goes through the gate, and it needs its own old-to-new id map.

---

## D-032 - A dominant-negative allele inverts the derived polarity (W016)

**The defect.** `perturbation_polarity` is derived from `method`, so
`conditional_overexpression` yields `increase`. Overexpressing a
dominant-negative or catalytically dead allele *reduces* the target's function,
so the derived polarity describes the direction of **expression** while the
experiment is a loss of function. This is the D-019 sign error one level below
the place `perturbation_polarity` was introduced to fix it.

**Found in the data, not hypothesised.** The `sirt1` records (18, one study)
include both wild-type Sirt1 overexpression and the deacetylase-dead
**Sirt1-H355Y** dominant negative. Both resolve to `target_id: sirt1`, both
derive `polarity: increase`, and so a gain of function and a loss of function
sat in one cell. Eight records carry the marker. Six of those eight are worse
than mis-signed: their `targeting_strategy` names *both* arms
("Rosa26 Sirt1-WT (Agrp-Sw) or Rosa26 Sirt1-H355Y"), so the record does not
identify which allele it measured.

**Decision.** A new WARNING, `W016_POLARITY_MAY_BE_INVERTED_BY_ALLELE`, fires
when a `molecular` record's targeting text matches a committed pattern of
function-reducing allele wordings (dominant negative, kinase/deacetylase dead,
catalytically inactive, `H355Y`, `K<n>R`, `D<n>A/N`, inactive mutant).

**It never rewrites the polarity**, for the same reason the outcomes
vocabulary's polarity annotations never rewrite a direction (D-021): a
mechanical guess at a sign is how an inversion gets into data in the first
place. The stored value is quoted in the warning so a curator can see what is
being questioned.

**Deliberately narrow.** The pattern does not match "mutant", "transgenic" or
"floxed", which would flag most of the corpus and therefore mean nothing. A
test pins both the wordings that must fire and the wordings that must not.

**Not closed.** The eight records keep their derived polarity and their
warning. Fixing them properly needs either an allele field or distinct target
ids per allele, and either is a schema change that should wait until a curator
has re-read that paper and settled which of the six ambiguous records describes
which arm. Recorded here so that the warning is not mistaken for a fix. Until
then those eight records should be excluded from any molecular analysis that
depends on the sign, and `W016` is how such an analysis can find them.

---

## D-033 - `axis` and `direction` use the same two words and are never printed bare

**The hazard.** `perturbation_polarity` and `direction` share a value
vocabulary. `increase` on the polarity means *more firing, or more of the
target gene product*. `increase` on the direction means *the readout went up*.
They are different statements about different ends of the experiment, and in a
report, a tooltip or a prompt they look identical.

This is D-019 arriving in a new place. D-019 was a reader recording the assay
number's direction instead of the outcome's; this would be a reader, or a
model, or a page, treating a manipulation's sign as a result's sign. The
v0.2.0 split created the hazard by giving the manipulation its own signed
field, which was the right thing to do and does not make the hazard less real.

**Decision.** Anything that renders a polarity renders it as a phrase, never as
the bare enum value. `syn4rch.atlas.signaling` carries `axis_phrase` on every
edge — "more firing", "less of the target gene product", "population or gene
product removed", "sign of the manipulation not established" — and the viewer
prints `PERTURBATION less firing · circuit activity (axis = decrease)`, with the
raw value shown only in parentheses after the words. The JSON carries a
top-level `axis_trap_note` stating the distinction, and the page renders it.

The stored values do not change. This is a presentation rule, because the
alternative — renaming the polarity members to something like `more`/`less` —
would make `perturbation_polarity` and `direction` harder to compare in code
where comparing them is legitimate, and would not stop a report printing the
bare word anyway.

**Where this is not yet enforced.** Only the atlas layer carries
`axis_phrase`. A future report or prompt that prints a polarity must add its
own phrasing; nothing mechanical prevents a bare `increase` from reaching a
reader elsewhere. Recorded as a gap, not as a solved problem.

---

## D-034 - `contested` is agreement below 0.7 over at least 3 records

**Decision.** An edge in the signalling graph is `contested` when
`agreement < 0.7` and `n_records >= 3`, where agreement is the share of records
reporting the modal direction. A second candidate rule — "two or more
directions each supported by two or more independent studies" — was measured
and **rejected as the contested rule**, and kept as a separately named
`minority_replicated` flag.

**Why rejected.** It flags only a couple more edges, and one of them is
`arc_agrp → food_intake` under more firing: 24 records, 9 studies, 9 labs, 21
reporting an increase and 3 reporting no change, agreement 0.875. That is the
most replicated claim in this corpus and the best-established result in the
field. Drawing it as contested would present agreement as conflict, which is
the opposite error from the one the flag exists to prevent and a worse one: a
viewer who sees the field's firmest finding marked "disputed" learns to ignore
the marking.

**Why the signal is kept anyway.** "Three independent studies report no change
while twenty-one report an increase" is real information about an edge, and
suppressing it to keep the picture clean would be the same sin in the other
direction. It is reported in the provenance readout under its own name, and it
is never drawn as conflict.

**Both rules are published, not buried.** `contested_rule` and
`minority_replicated_rule` are top-level fields in
`viewer/data/signaling.json`, stated on the page, and pinned in
`tests/test_atlas_signaling.py`. A reader can check what the word means without
reading the code.

---

## D-035 - The usability predicate reads the polarity, not the valence

**Decision.** `syn4rch.feasibility.coverage.has_usable_intervention` used to
require a `valence` that was neither null nor `unclear`. It now reads
`perturbation_polarity`. Usable records go from **177 to 347 of 362**, the 170
added being exactly the molecular ones, with 15 still excluded as
`PERTURBATION_POLARITY_UNCLEAR`.

**Why this is a correction and not a loosened threshold.** The old predicate
asked a sensible question — "does this record say which way the manipulation
went?" — of a field that, after D-025, no longer answers it. A conditional
knockout's `valence` is honestly `unclear` and SCHEMA.md 3.4 forbids coercing
it, while the manipulation itself has a perfectly well defined sign that now
lives in `perturbation_polarity`. Leaving the predicate on `valence` would have
meant the usability gate excluded precisely the records v0.2.0 was written to
make representable. No threshold constant changed; all ten in
`feasibility.gates.THRESHOLDS` are pinned literally by test.

**It is the change in this version that most needs watching**, because it is the
one that makes a published number larger, so it is instrumented rather than
asserted: `key_version="0.1.0"` reproduces the old count of 177 exactly,
`compare_usability` reports both, and the report prints the admission broken
down by perturbation class so the +170 cannot be mistaken for a general
improvement.

**"Usable" is not "learnable", and the two must be read together.** These 347
records have the fields an example of the direction task needs. The identity
analysis (D-027) says 100% of the 170 molecular ones are unlearnable under
study or lab grouping. Both are true and neither cancels the other: a record
can be well formed and still have nothing in the training set that shares its
manipulation. Reporting usability alone would overstate the corpus; hiding the
molecular records at the usability stage would have concealed *why* it is
overstated. They are published side by side.

**Fewer cells pass the gates, and that is the expected direction.** Strict key:
237 tasks give 2 feasible, 2 marginal, 233 insufficient; 93 (population,
outcome) pairs give 1 feasible, 1 marginal, 91 insufficient. The pooled key on
the same records gave 8 marginal pairs where the strict key gives 1, so the
correct key costs 7 marginal pairs.

**Both surviving `feasible` verdicts sit at condition specificity 0.** They
carry `CONDITIONS_ALL_UNKNOWN` and are counted separately.
`MIN_CONDITION_SPECIFICITY_GATE = None` is published as a constant so that the
absence of a specificity gate is a decision on the record rather than an
oversight. Inventing a threshold there would be inventing the knob `gates.py`
exists to warn against.

---

## D-036 - The condition dimension fixes conflict detection and does not help prediction

**The measurement.** Adding the condition to the model's cell key makes
prediction slightly *worse*, not better. Paired against the axis-only cell,
under `study_id`: −0.006 [−0.015, +0.000] pooled and −0.028 on the rows whose
cell was seen. The conditioned rung reaches only 88 of 361 predictions, so most
rows fall through it to a coarser one.

Why, in one table (`reports/baseline_benchmark.md`):

| cell | cells | mean rows per cell | singleton cells | in-sample purity |
|---|---|---|---|---|
| population x outcome | 95 | 3.80 | 50.5% | 0.632 |
| + perturbation axis | 203 | 1.78 | 69.5% | 0.809 |
| + condition | 249 | 1.45 | 81.1% | 0.875 |

In-sample purity rises monotonically and means nothing: it is a ceiling
measured on the rows it was computed from, and a cell holding one row is
trivially pure while contributing nothing to a held-out row of the same cell.
Four fifths of the conditioned cells hold exactly one record. The refinement
that fragments the data hardest looks like the one that explains it best, and
the two columns are published side by side so that it cannot.

**This does not retract D-023.** The condition dimension does what it was added
for. It removed 16 of 18 false contradictions and it is the only reason the
remaining 9 conflicting claims can be believed. Conflict detection asks "are
these two records comparable?", which the condition answers. Prediction asks
"have I seen enough records like this one?", which the condition makes harder by
construction, because every slot it adds divides the evidence.

**Both are true at once and neither is the headline.** A schema change can be
right for the dataset's integrity and useless for its predictive power, and
reporting only the half that flatters the change would be the error this project
exists to avoid. The condition dimension is kept, in the claim key and in
`record_id`, on the integrity argument alone.

**What would change the verdict.** Conditions curated rather than inferred.
`challenge` is null on all 362 records and 153 know no condition slot at all, so
today the conditioned cell is mostly splitting `unknown` from `unknown`. The
measurement above is of a condition dimension that exists in the schema and not
yet in the data, and it should be re-run once a curation pass has filled it.
Until then "conditioning hurts" means "conditioning on absent data hurts",
which is not a finding about biology.

**Also measured, and reported separately rather than blended.** Under
`study_id` the axis improves ranking (one-vs-rest AUC 0.567 to 0.639) and makes
calibration distinctly worse (expected calibration error 0.001 to 0.193),
because cell frequencies drawn from cells of one or two rows are overconfident.
Accuracy, ranking, calibration and uncertainty are four numbers and this project
does not combine them.

---

## D-037 - Condition provenance is per dimension, not per record

**Supersedes D-024**, which required a `condition_locator` on a record whose
`condition_source` was `methods`. The requirement was right and the shape it
was attached to was wrong; the locator requirement survives, per entry and
whatever the source (A008).

**Decision.** v0.3.0 removes the record-level `condition_source` and
`condition_locator` fields and replaces them with `condition_evidence`, a list
holding one entry per *filled* condition dimension. Each entry carries the
dimension, the value it supports, the source class (`passage` or `methods`), a
locator, a verbatim quote, and the name of the rule that read it.

**Why.** The scalar pair could describe one provenance, and a single record
genuinely has two. The sex of the animals is stated once in a paper's Methods
and holds for every figure in it. The acute metabolic state at the moment of
measurement is stated, when it is stated at all, by the passage the record was
read from. A record filling `sex` from Methods and `challenge` from its own
passage cannot say so with one `condition_source`, so the honest options were a
field that describes half the record or no field at all. In practice the second
happened: `condition_source` was null on all 362 records of the first corpus,
and the condition gap went unfixed partly because the schema had no way to
record a fix.

**Cost.** The schema gains its first list-of-objects field, so `to_dict`,
`from_dict` and the coercion path in `record.py` grew a case. Four validation
rules replaced one (A008, A009, A010, A011), because a missing locator, a
missing quote, a malformed entry and an entry contradicting its own field are
four different mistakes with four different fixes, and a report that collapses
them tells a curator only that something is wrong with the condition.

**Nothing was lost.** Both removed fields were null on every record in the
corpus. This change was free in data and would have stopped being free the
moment one record used them.

---

## D-038 - A paper may settle a dimension it holds constant, and may not settle one that varies

**Decision.** The condition fill reads `sex`, `diet` and `age_group` from a
paper's Methods and applies the value to every record from that paper **only
when the paper holds that dimension constant across every experiment it
reports**. Where the Methods describe two cohorts on that dimension, the fill
assigns nothing and the field stays `null`.

**Why.** This is the whole argument for why a Methods-sourced condition is
admissible at all. A paper stating "male C57BL/6 mice were used" has told you
the sex of the animals behind every figure in it. A paper comparing a chow group
against a high-fat group has told you that both diets exist in the study and
nothing whatever about which group produced Figure 3C. Treating the second like
the first would stamp one of the two diets across every record of a
diet-comparison paper -- and diet-comparison papers are most of this corpus's
high-fat records, so the error would land precisely where it does most damage.

**`varies` is not an extraction failure.** It is the paper declining to answer
the question at the scope the question was asked, and `null` is the correct
record of that. A pass that reported `varies` as a gap to be closed by better
parsing would be measuring the wrong thing: the information is not in the
Methods to be found.

**Cost.** The fill cannot reach the dimensions that matter most for conflict
detection in exactly the papers where those dimensions were deliberately
manipulated. Closing those requires reading each figure panel's own cohort,
which is a human curation task, not a parsing task.

---

## D-039 - `not_reported` asserted from one passage may be corrected by the Methods; a curated value may not be overwritten

**Decision.** The fill applies a Methods-sourced value over a stored
`not_reported`, and records that as a *correction*. It never applies one over a
stored positive value: where the two disagree the disagreement is reported and
the curated value stands.

**Why.** `not_reported` is a claim about the paper -- that it was checked and
does not say (SCHEMA.md 3.3b). The extraction that produced the first corpus
could see one passage, so a record reading `sex: not_reported` means "this
passage does not say", which at paper scope is a different and often false
statement. One record's `extraction_notes` says so outright: *"Species not
stated in this passage, so left null."* Correcting those is not overruling a
curator's judgement; it is widening the scope the judgement was made at.

A stored *positive* value is the opposite case. It was read from the record's
own passage, which describes this measurement, while the Methods describe the
study. The passage is the more specific statement and it wins. A paper that used
both sexes and a record that says `male` are not in conflict at all -- the
record is simply more precise than the paper-level reading.

**And genuine disagreements are not resolved.** A stored `male` against a
Methods section that describes only females is a real contradiction, and this
project does not have a rule for picking a side. Those records are listed in
`reports/condition_fill_0.3.0.md` with both values and the locator, and they
keep saying what the original curation said until a human reads the paper.

---

## D-040 - No Methods sentence can state a measurement's challenge

**Decision.** "Mice were allowed ad libitum access to water and standard chow",
and every sentence like it, may not be read as `challenge = ad_libitum`. The
rule that would do so is implemented as a named refusal,
`refused_housing_statement_is_not_a_challenge`.

**Why.** It describes how animals were kept between experiments. The same paper
measures glucose during a glucose load, food intake after an overnight fast, and
energy expenditure in a metabolic cage, and the housing sentence is true of the
animals in all three while describing the state of measurement in none. Nearly
every paper in this corpus contains such a sentence, so accepting it would have
filled `challenge` on most of the 362 records, taken the dimension's coverage
from zero to near-total, and made the condition dimension look solved.

**This is the most dangerous available mistake in this pass**, which is why it
is a tested refusal rather than a matter of care. It is wrong in the direction
that produces an impressive number, it is wrong on hundreds of records at once,
and because `challenge` is a `record_id` component it would have re-identified
the whole corpus on the strength of a sentence about cage husbandry.

**Revised before commit, and widened.** The first version of this rule refused
husbandry sentences specifically. Testing the extractor against
`doi:10.1038/s41467-018-05711-6` showed the rule being routed around by a
sentence that is not husbandry at all: *"Ad libitum fed rats were anesthetized
by an i.p. injection of ketamine..."* -- a stereotaxic surgery preparation. It
was read as a paper-level `challenge` of `ad_libitum`, which would have been
stamped on all 34 of that paper's records on the strength of the animals'
feeding state at the moment they were anaesthetised.

Patching the phrase list would have left the next wording to be found. The
correct rule is the general one: **no sentence in a Methods section is
admissible evidence for `challenge`**, because the Methods describe the study's
procedures and the field describes the state of one measurement. The Methods
route for this dimension is therefore removed rather than filtered, and the
refusal is named `refused_methods_cannot_state_a_measurement_state`. `challenge`
has exactly one admissible source: the record's own passage, where the paper
says which state *this* readout was taken in.

---

## D-041 - An assay does not supply its own challenge

**Decision.** `challenge` is not derived from the outcome. `outcome_id =
glucose_tolerance` does not yield `glucose_load`, and `GTT`, `ITT`, `glucose
tolerance test` and `ipGTT` in a record's `outcome_raw` yield no challenge.
Named refusal: `refused_assay_implies_its_own_challenge`.

**Why.** It is true and useless. A glucose tolerance test is a glucose load by
definition, so the derived value is a function of a key component the record
already carries. In the claim key it discriminates nothing, since any two
records it separates were already separated by `outcome_id`. In the model's cell
key it adds a feature with no information. What it would do is raise the
reported condition coverage by the 43 records whose outcome is a tolerance test,
and re-identify every one of them, in exchange for nothing.

**What the dimension is actually for.** The acute state of measurements whose
assay does *not* fix it: food intake, body weight, energy expenditure, blood
glucose. "Food intake in the fed state" and "food intake after a fast" are two
different experiments that `outcome_id` alone cannot separate, and separating
them is the entire value of the field. Those come from the record's own passage.

**Consequence for the coverage figure.** `challenge` coverage after this pass is
low, and it is low honestly. A reader comparing it against a figure that
included the assay-derived values would conclude this pass did worse, when the
difference is that the other figure counts restatements of the outcome column.

---

## D-042 - The unevidenced values already in the corpus are counted, not rejected

**Decision.** A condition dimension holding a positive value with no evidence
entry raises `W017_CONDITION_WITHOUT_EVIDENCE`, a warning. It is not a hard
finding, so the 362 records curated before v0.3.0 remain admissible.

**Why.** Making it hard would reject the entire dataset at the moment the rule
was introduced, which is not a standard being enforced but a dataset being
deleted. The values are not worthless: they were read from a passage by an
extraction pass whose output was quote-verified, and they are what the first
corpus was measured on. What is missing is the sentence that lets a reader check
one without re-reading the paper.

**What the warning is for.** Counting. The QC report says how many condition
*values* and how many *records* are in that state, so the remaining gap is a
number in a report rather than an impression. A later pass closes it by
attaching evidence to values it confirms, and the count is how progress is
visible.

**Why not retro-fit evidence by searching for a sentence that agrees.** Because
that is confirmation bias mechanised: searching a paper for a sentence matching
a value you already hold will find one in a paper that mentions both sexes, and
the resulting entry would document agreement that was never tested. The fill
derives the value from the paper first and compares second, which is the only
order in which a disagreement can surface.

---

---

## D-043 - An age at the start of a protocol is not the age at measurement

**Decision.** An age statement is refused when its sentence also carries a
protocol-start or protocol-duration cue -- "starting at", "beginning at", "from
N weeks of age", "for N consecutive weeks". Named refusal:
`refused_protocol_start_age_is_not_measurement_age`. Separately, if **any** age
statement in a paper's Methods is refused for any reason, the paper's `age_group`
verdict may not be `constant`.

**Why.** Both halves came from one sentence in
`doi:10.1038/s41467-018-05711-6`: *"C57BL/6 mice were fed with HFD for 12
consecutive weeks (starting at 5 weeks of age)."* The first draft of the
extractor read 5 weeks, converted it, and returned `age_group = juvenile`. Those
animals were five weeks old when their diet began and roughly seventeen weeks
old when anything was measured in them. `juvenile` is not a slightly imprecise
answer, it is the wrong band, and it is wrong about the only thing the band is
for: the animal's state when the readout was taken.

**The second half is the subtler bug.** That paper's first Methods paragraph
says "8-10-week-old ... male Sprague Dawley rats and C57/BL6 WT mice", which the
rat refusal (D-003's mouse-only boundaries) correctly skipped. With that
statement gone the 5-week statement stood alone, looked unanimous, and was
promoted to `constant`. A refusal is evidence that the paper states an age the
rule cannot use; treating the survivors as consensus turns a declined reading
into a confident one. That is cherry-picking implemented as a loop, and it is
invisible in the output -- the verdict says `constant` and names a real quote.

**Also.** No age band is assigned unless the paper's species is established as
mouse. The same paper returned `species = None`, because its Methods name both
mice and rats, while still receiving a mouse band. `assign_age_group` refuses
every non-mouse by design (D-003), and passing it a hard-coded
`"mus_musculus"` for a paper whose species is unresolved bypasses that refusal
rather than satisfying it.

---

## D-044 - A `constant` verdict must be established, not merely undetected

**Decision.** The condition extractor may return `constant` for a dimension only
when the absence of a second value has been *established*. Every refused
sentence, every skipped statement and every term the vocabulary does not
recognise forces the verdict away from `constant`. Specifically: an unclassified
diet-shaped token (`\b\w*\s?diet\b`, `chow`, or an abbreviation of the `HFD` /
`ND` / `NCD` / `LFD` shape) yields `varies` and names the token; a disjunction
("fed either a ND or a HFD") yields `varies` whatever its arms are called; a
cohort sentence scoped to one figure or one procedure does not vote at paper
level; and a disjunctive sex ("male or female mice") votes nothing rather than
`both`.

**Why.** Every one of these was a real wrong answer in the first draft, found by
reading all 42 papers' extracted values against the sentence each came from. The
worst was `doi:10.1016/j.isci.2024.110259`:

> "Mice were fed either a ND (18% fat ...) or a HFD (45% fat ..."

a plain diet-comparison study, returned as `constant high_fat`. The cause was
that the synonym table I wrote listed `NCD`, `SCD` and `normal diet` but not the
bare `ND`, so one arm of the comparison was simply invisible, the other stood
unopposed, and the rule concluded the paper held diet constant. All four of that
paper's records would have been stamped `high_fat` on the strength of a missing
dictionary entry.

**The shape of the bug is what matters, not the instance.** `constant` is a
claim about what the paper does *not* contain, and a claim about absence drawn
from a lexicon is only as good as the lexicon. Adding `ND` fixes this paper and
leaves the next abbreviation to be found by whoever reads the dataset. So the
rule fails closed: anything diet-shaped that the table cannot classify is
evidence that the table is incomplete, and an incomplete table must produce a
refusal.

The same reasoning retired three other readings. A sentence saying "Experiments
(Figure 3) were performed in female mice only" was read as a paper-level
`female`, when the word "only" beside a figure reference says precisely that the
paper has other figures with other animals. A slice-preparation sentence reading
"male or female mice were deeply anesthetized" was read as `both`, when a
disjunction leaves any individual record's sex unresolved. And "animals with the
expected ... response to chow" was read as `diet = chow`, when the chow in that
sentence is a stimulus presented during a recording.

**The contrast with D-040 is the useful part.** A stereotaxic surgery sentence
*is* admissible evidence for `sex`, because an animal's sex does not change
between surgery and readout, and inadmissible for `challenge`, because its
feeding state does. What makes a sentence evidence is not where it sits in the
Methods but whether the property it states is invariant between that sentence
and the measurement. `sex` is invariant, `diet` is invariant over weeks and so
survives a housing sentence, `age_group` is invariant only if the sentence
describes the animals at the experiment rather than at the start of a protocol
(D-043), and `challenge` is invariant over nothing at all (D-040).

**One marked exception, where the value is kept and the doubt is published.**
Where a record's quote carries a state token that a trap rule excluded *and* a
different token that matched, the extraction is returned rather than refused,
under the distinct rule name
`quote_names_the_state_with_a_competing_token`, with both tokens named in the
notes. The case that prompted it:

> "chow-fed Htr2cSim1KO mice exhibited an anorectic response to a sub-anorectic
> dose of dFen (1 mg/kg) that otherwise did not suppress refeeding (during the
> first 2 h) in wild-type littermates"

`chow-fed` is correctly excluded as a diet compound, which leaves `refeeding`
unopposed, which is the promotion-after-a-refusal shape this entry otherwise
forbids. `refed` is probably right, and `challenge` is the scarcest dimension in
this corpus, so refusing it costs more than publishing the doubt. The rule name
is what makes that choice visible: these extractions are counted on their own
line in the fill report, so the set a human has to re-read is small and named
rather than mixed into the rest.

---

## D-045 - The age band's evidence is its quote plus the rule that read it

**Decision.** All four condition dimensions carry their own evidence entry,
`age_group` included. An `age_group` entry quotes the sentence stating the age
and names the rule that converted it (`EVIDENCE_DIMENSIONS`, the identity
mapping in `EVIDENCE_FOR_CONDITION`). An alternative was built and reverted:
entries keyed on `age_value` instead, so that the quote and the value it
supports would be the same fact.

**The objection that prompted it is a good one.** A reviewer reading the
fixtures pointed out that an entry pairing *"mice were 12 weeks old at
surgery"* with the value `young_adult` offers a quote that does not contain the
value. A009 asks for the sentence a value came from and A011 asks the entry to
agree with its field, and for a value no paper writes, those look
irreconcilable.

**Why it was reverted anyway.** Two measurements.

Half of this corpus's age statements are ranges -- "8-10-week-old",
"6-24 weeks", "2-3 months" -- 16 ranges against 15 scalars across the 42
papers. A range yields a band when both of its ends fall in one band, and
yields no number to put in `age_value`. Keying the evidence on `age_value`
would therefore have left about half the filled bands with no entry at all,
flagged by `W017` as unevidenced, while the sentence that settles them sat in
the paper unquoted. That is a worse account of the evidence, not a better one.

And the objection proves too much. **No condition value is literally present
in its own quote.** `diet = chow` comes from "animals received ad libitum
access to standard rodent chow (ssniff V1554)", and `sex = male` from "C57Bl/6,
db/+ and db/db male mice were purchased from CLEA Japan". Every entry in this
schema pairs a sentence with the schema's term for what the sentence says. The
mapping is the curation.

**What is genuinely different about the age, and what follows.** Its mapping is
*parametric* rather than lexical: arithmetic against a boundary table that
could be revised without a single quote becoming wrong (D-003 records that
those boundaries were a judgement call). So an `age_group` entry is the one
whose audit needs both halves -- the `quote` for what the paper stated, and
`rule` for how the band was reached -- where the other three are checkable from
the quote alone. That is the reason `rule` is on the entry at all.

**A limitation this surfaced, which is not fixed.** `W009` checks a stored
`age_group` against `assign_age_group(age_value, age_unit, species)`, and for a
range-stated age `age_value` is null, so the rule returns `not_reported` and
the check cannot run. **A band derived from a range is not mechanically
verifiable by anything in this repository**, whichever field its evidence hangs
on. It is verifiable by a reader, from the quote and the rule name, which is
why both are stored. Making it machine-checkable needs a representation for a
stated range, and that is a schema change with no demand behind it yet.

**On the process.** This entry records a change that was made, tested green,
reviewed, rebuilt the other way, measured, and reverted. The measurement is
what decided it, and it was cheap: counting scalars against ranges took one
command. The first version was reverted because it was right about the
principle and wrong about this corpus.

---

## D-046 - A version artefact is set aside on both sides, in every comparison

**Decision.** `S005` and `S006` -- the findings that fire when a file's
`schema_version` is not the validator's -- are excluded from **both** sides of
every before/after comparison in `syn4rch.curate.migrate`, including the
`records_hard_clean` figure, via one shared `_non_artefact` helper. And
`read_records` strips the fields retired after 0.2.0 (`condition_source`,
`condition_locator`) so that a file this module wrote stays readable by the
code that reads it.

**Why, and what it was hiding.** Pinning the 0.1.0 -> 0.2.0 migration's target
to the literal `0.2.0` (so it would not drift when the schema moved to 0.3.0)
exposed two latent faults at once.

The first stopped the migration from writing anything: `cmd_run`'s gate
filtered the artefacts out of the *baseline* and not out of the result, so once
the output declared a version older than the validator, every output record
tripped `S006` and the gate read it as `+362` hard findings the migration had
introduced. The report and the gate disagreed about what they were comparing,
and the disagreement was invisible for as long as the target version tracked
the schema package.

The second was a published number. `reports/schema_migration_0.2.0.md` read
*"0 of 362 input records carried no hard finding at all, and 362 of 362 do
now"*, which flatters the migration and contradicts the sentence immediately
after it, which sets `S005` and `S006` aside and states that neither file
carries a hard finding. The `0` was produced entirely by a 0.1.0 file being
validated by newer code. Both columns now read 362, and the regenerated report
says what the correction was.

**And the report stays recomputable.** Once 0.3.0 removed two fields, 0.3.0
code could no longer read a 0.2.0 file at all, which silently cost this
repository the ability to check the published migration numbers --
`docs/CURATION_PROTOCOL.md` 6 forbids exactly that. Teaching the reader which
keys were retired restores it. The strip is targeted rather than
`strict=False`, which would also have swallowed a misspelled field name, and a
retired field found holding a non-null value is still refused: both were null
on all 362 records, and that is the only reason retiring them was free.

**One consequence to read carefully.** Recomputation runs under current code,
so a rule introduced *after* a migration appears in that migration's table.
`W017`, added by v0.3.0, counts 209 records with a condition value carrying no
verbatim sentence. It describes a gap the 0.2.0 migration never claimed to
close, and the report now says so rather than leaving a reader to wonder how a
2026 migration failed a rule that did not exist.

---

## D-047 - What filling the conditions measured, including the prediction that failed

**The pass.** `syn4rch.curate.fill_conditions` read the 42 curated papers and
wrote 362 records at v0.3.0. Per dimension — filled / corrected / confirmed /
conflicted: `sex` 112 / 7 / 55 / 2, `diet` 86 / 11 / 7 / 0, `age_group` 36 / 0 /
0 / 0, `challenge` 25 / 0 / 0 / 0. 227 of 362 `record_id`s changed. 339 evidence
entries on 252 records, 314 from Methods and 25 from a record's own passage; all
339 resolved to the passage their locator names with the quote found there.
Per dimension: 174 for `sex`, 104 for `diet`, 36 for `age_group`, 25 for
`challenge`.
165 named refusals across 13 rules. Nulls: `sex` 224 → 112, `diet` 197 → 111,
`age_group` 286 → 250, `challenge` 362 → 337.

**What it bought, and it is not the headline I expected.** Cross-study
conflicting claims fell from 9 to 4 and the last directly-opposed pair from 1 to
0 (QC, 362 records; the feasibility report's 2 → 1 counts the same change over
its 347 usable records). But the decomposition matters more than the count. Of
the 9 conflicts, **4 split because two records state genuinely different
conditions** and **5 split because one record states a condition the other
leaves unknown**. The second kind is not a contradiction resolved; it is a
comparison withdrawn, and this corpus cannot say whether the two experiments
agreed.

The real gain is in *replication quality*. Before the fill the corpus held 11
multi-study claims, **8 of which matched at condition specificity 0** — grouped
together partly because nobody knew their conditions. After it there are 5, and
**none of them matches on nothing known**: 4 sit at 1 of 4 stated slots and 1 at
2 of 4. Records knowing no condition slot at all fell from 153 to 46. So the
count halved and what remains is qualified, which is the trade this dimension
exists to make.

**The prediction I wrote down beforehand, and it failed.** `reports/` recorded
that 7 multi-study cells had been lost *only* because one paper stated a
condition the other omitted, and I expected a curation pass to let some of them
be adjudicated. It went the other way: **7 became 10**. Filling the conditions a
paper reports makes that paper *less* comparable with a paper that reports none,
because `unknown` matches no stated value. Partial condition curation therefore
increases unadjudicable non-comparability, and only curating *both* sides of a
pair can reduce it. The honest summary of the cost line is that 32 of 35
multi-study pooled cells gave up their replication — 18 to different
perturbations, 4 to genuinely different conditions, and 10 to this asymmetry.

**And D-036 stands, now measured on data instead of absence.** D-036 said
conditioning the model hurts slightly and that the measurement should be re-run
once the conditions were curated, because it was then "of a condition dimension
that exists in the schema and not yet in the data". Re-run: conditioned against
axis-only, under `study_id`, **−0.006 [−0.021, +0.006]**, 0.479 against 0.485,
McNemar p=0.754, verdict `no_difference`. The pooled figure is the same −0.006 it
was on empty conditions, and on the rows whose conditioned cell was seen in
training it is now **−0.061**, worse than the −0.028 it was before.

The mechanism is the one D-036 named, and filling the dimension sharpened it.
The cell refinement table now reads 95 → 203 → **268** cells, singletons 50.5% →
69.5% → **81.3%**, in-sample purity 0.632 → 0.809 → **0.898**. Real condition
information makes the conditioned cell *finer*, so more of it is singletons, so
it generalises worse. Four fifths of conditioned cells hold one record.

**One thing the re-run did establish**, reported separately because it is a
different question: on the 49 rows whose conditioned cell was seen in training,
the conditioned baseline scores 0.673 against 0.327 for guessing the mode —
**+0.347 [+0.077, +0.695]**, McNemar p=0.006, verdict `established`. Where a
conditioned cell exists it is strongly informative. It just loses to the
unconditioned axis cell on those same rows (0.673 against 0.735), because the
coarser cell has more evidence in it. Both facts are true and neither is the
headline.

The circuit-activity signal is untouched: 0.517 against 0.397, +0.121 [+0.008,
+0.271], McNemar p=0.156, still `suggestive` (D-028).

**So the condition dimension is still kept on the integrity argument alone**, and
that argument is now stronger than it was: it is what took 8 of 11 replications
off a match on mutual ignorance. It remains useless for prediction, and three
schema versions of work on it have not made it otherwise.

**What would change the verdict now.** Not more of this. `challenge` reaches 25
records because the only admissible source for it is the record's own passage
(D-040), and `age_group` reaches 43 because most stated ages are at surgery, at
vendor arrival or at the start of a diet (D-043). No claim in the corpus reaches
specificity 4 of 4, so no replication figure here describes two experiments known
to have matched on all four conditions — the same sentence this project wrote
before the pass, still true after it. Closing that needs the condition recorded
per figure panel by a human reading the paper, which is curation, not parsing.

---

## D-048 - Pre-registration: what curating more circuit papers must achieve, and what would abandon the plan

**Written before any new paper is curated**, and committed first so that the
result can fail. D-028 concluded that curating more optogenetic and chemogenetic
papers was the first measured reason this project had to prefer one paper over
another. This entry fixes what "it worked" means, because a target chosen after
seeing the data is not a target.

**The hypothesis.** A circuit manipulation is specified by a population and a
sign — "activate AgRP neurons" — and that specification recurs across
laboratories. A molecular manipulation is specified by a gene, and in this
corpus no gene is manipulated by more than one laboratory. So adding circuit
papers should add *replicated cells*, where adding molecular papers adds only
singletons (D-027).

**The baseline, measured today on the 362 committed records:**

| class | rows | studies | labs | cells | cells with >=2 labs |
|---|---|---|---|---|---|
| circuit_activity | 170 | 26 | 26 | 56 | **17** |
| molecular | 170 | 16 | 16 | 132 | **0** |
| structural | 17 | 1 | 1 | 12 | 0 |
| pharmacological | 5 | 2 | 2 | 3 | 0 |

Cell here is `(population, perturbation, outcome)`. Under the full claim key,
which also carries the condition, the whole corpus holds 5 multi-lab claims
(D-047) — conditions split those 17 circuit cells down to a handful, and that is
expected rather than a new problem.

**A correction to my own earlier framing.** I had said circuit perturbation
identity "recurs across laboratories by construction". It does, and for a weaker
reason than the phrase implies: `perturbation_key` is
`(class, polarity, target)` and a circuit record has no target, so the entire
circuit half of the corpus holds **3** distinct identities and 0% of rows have an
unseen identity in leave-one-study-out. That is a property of a coarse key, not
of the biology. The population — which is what actually distinguishes one
circuit manipulation from another — lives in the cell key, not the perturbation
key. So the claim is about *cells*, and the table above is the honest version of
it. Any metric defined on perturbation identity alone would have been trivially
satisfied and would have measured nothing.

**The material.** 290 papers with a circuit-manipulation extraction unit have
retrieved full text and are not yet curated, carrying 4,436 candidate units.
The plan is not limited by what is on disk.

**The target.** Curate **60** uncurated circuit papers, then evaluate **once**:

1. *Primary.* The circuit-subset direction comparison
   (`perturbation_axis_cell_majority` against `global_majority`, leave-one-study-out)
   reaches verdict `established` — the study-clustered bootstrap interval excludes
   zero **and** the row-level exact McNemar test reaches p < 0.05. Today it is
   `suggestive`: +0.121 [+0.008, +0.271] with p=0.156.
2. *Secondary.* Circuit cells with >=2 laboratories rises from 17 to **>=40**,
   i.e. the added papers land in existing cells rather than founding new
   singletons.

**Evaluated once, at 60 papers.** Not after each batch. A suggestive result
tested repeatedly becomes significant by chance, and this project has one
suggestive result it would be easy to flatter. One mechanism check is permitted
at 30 papers, on the secondary metric only: if multi-lab circuit cells have not
risen at all by then, the hypothesis is already wrong and curation stops.

**What would abandon the plan.** Any of:

* multi-lab circuit cells flat at 30 papers (the mechanism check above);
* the secondary metric reached while the primary is not — the papers replicate
  each other and direction is still unpredictable, which would say the signal
  in D-028 was a 116-row artefact;
* the primary reached only after the planned stopping point is passed.

**And the result that would be reported either way.** "Adding 60 circuit papers
did not make direction of effect predictable" is a publishable finding about the
state of this literature, and is the answer this project was built to be able to
give (README, milestone 1). It is not a reason to curate a 61st paper.

**One thing this does not fix.** `perturbation_key` still has no slot for a
ligand or a population, so pharmacological identity remains unmeasurable and
circuit identity remains coarse. Changing it is a schema version and is out of
scope here; it is recorded so that the 3-identity figure above is not mistaken
for a finding.

---

## D-049 - Thirteen polarity traps annotated, and the eight records that could have hit them checked

**Found by a curator, not by a test.** A reader working the circuit pass hit this
sentence in `doi:10.1038/s41467-018-08239-x`:

> "an increase in RER … a decrease in fat utilization … an increase in
> carbohydrate utilization"

One calorimetry measurement, stated three ways. `data/vocab/outcomes.json` listed
`fat oxidation` and `lipid oxidation` as **plain** synonyms of `rer` with no
`inverted_synonyms` entry — but RER is VCO2/VO2, so burning more fat *lowers* it.
A curator following the synonym list would enter `rer: decrease` from a sentence
that supports `rer: increase`, and nothing in the pipeline would object: the
quote is verbatim, the locator resolves, the direction is stated in the text.

**The trap was not one entry.** Sweeping every entry against its own id using the
vocabulary's own rule — *an increase in the term implies a decrease in this
entry's quantity* — turned up 13 unannotated terms across 9 entries:

| entry | newly annotated as inverted |
|---|---|
| `rer` | fat oxidation, lipid oxidation |
| `sleep` | wakefulness, arousal, insomnia, sleep fragmentation |
| `blood_glucose` | hypoglycemia |
| `blood_pressure` | hypotension |
| `body_weight` | weight loss |
| `core_temperature` | hypothermia |
| `food_intake` | hypophagia |
| `locomotor_activity` | hypoactivity |
| `osmolality_avp` | hyponatremia |

The `hyper-` halves of those pairs are *not* inverted and are left alone;
`substrate utilization`, `fuel preference`, `glucose clearance` and
`hyperinsulinemic-euglycemic clamp` name no direction at all and are neither.

**And the records that could already have hit them were checked, one by one.**
Eight records in the corpus have an `outcome_raw` matching an inverted synonym of
their own `outcome_id`. **All eight are correct.** Three of them matched `weight
loss`, which was *not* annotated at the time they were curated, and the readers
inverted it anyway: "activation of POMC^ARC reduced fasting-induced body weight
loss" is stored as `body_weight: increase`. So this annotation is **preventive,
not corrective** — nothing was rewritten, because nothing was wrong.

That distinction is worth stating rather than leaving implied. An annotation pass
that silently "fixed" eight directions would have been indistinguishable in the
diff from one that corrupted eight directions, and the only way to tell is to
check each against its quote and say so.

**Annotation only.** No code reads `inverted_synonyms` (D-021: the vocabulary
annotates the known traps, it does not resolve them, because the fix is a reading
judgement). No id changed, no stored direction changed, no key changed, so no
number in any report moves. The outcomes vocabulary goes to 0.3.0 with a
changelog entry.

**One defect found and not fixed.** `osmolality_avp` bundles plasma osmolality
with AVP, and the two can move in opposite directions — AVP rises in SIADH while
osmolality falls. `hyponatremia` is annotated inverted with respect to
osmolality, the quantity it names, and the entry's notes now ask a record here to
say in `extraction_notes` which of the two its quote measures. Splitting the id
is a vocabulary change with no demand behind it yet: one record in the corpus
uses it.

---

## D-050 - A projection-specific manipulation is a different manipulation

**Decision.** `population_key` now distinguishes a projection-targeted
manipulation from a somatic one: `arc_agrp->pvh`, not `arc_agrp`. The
`projection_specific` flag decides it, and a record that declares itself
projection-specific while naming no target keys as `...->unknown` rather than
rejoining the somata — the same reasoning as the `raw:` prefix on an unresolved
population (D-010).

**Why, and found by a curator mid-pass.** A reader working the circuit batches
pointed out that AgRP→PVH, AgRP→LHA and AgRP→BNST all key as `arc_agrp`, so
three separately targeted projections pool into one claim. That costs twice, in
exactly the two directions this project has already paid for once (D-023):

* two laboratories manipulating *different* projections read as two
  laboratories replicating one claim — manufactured replication;
* a somatic result disagreeing with a terminal result reads as a contradiction —
  manufactured conflict. The same reader's batch contains MCH somata
  (`food_intake` increase) and MCH→ME terminals (`food_intake` no_change), which
  would have pooled.

**It was latent, and measured before being fixed.** On the committed 362-record
corpus the change moves cells 274 → 277 and moves **nothing else**: cross-study
conflicting claims stay at 4, directly opposed at 0, multi-study claims at 5. The
one pooled cell was a single study, a single laboratory and four `increase`
records, so it cost nothing yet. It would have started costing the moment the
new circuit papers landed, which is why the fix went in before they were
admitted rather than after.

**And it does not move the pre-registered target.** D-048's secondary metric is
circuit cells with two or more laboratories. Under the corrected key the
baseline is 59 circuit cells of which **17** are multi-laboratory — the same 17
as before. Only the denominator moved. The target "17 to >=40" therefore stands
exactly as written, which matters: a key change that moved the baseline would
have been goalpost-moving however well argued, and the way to tell the
difference is to measure both and publish the pair.

**`record_id` is unaffected.** It hashes `population_raw`, not
`population_key`, so no record was re-identified by this.

---

## D-051 - The 30-paper mechanism check: not flat, and not on track

**The check D-048 permitted.** Thirty of the sixty circuit papers are curated:
96 records admitted from 30 papers, 0 rejected, 3 of those papers declining to
yield anything on scope and 7 more yielding nothing from the passages supplied.
The corpus is 458 records from 61 studies.

| | records | cells | cells with >=2 labs |
|---|---|---|---|
| circuit, before | 170 | 59 | **17** |
| circuit, at 30 papers | 265 | 116 | **23** |
| molecular, unchanged | 170 | 132 | **0** |

**The abandon condition does not fire.** Multi-laboratory circuit cells rose,
so the hypothesis that circuit manipulations recur across laboratories is not
refuted, and curation continues to 60 papers.

**But the rate is not on track, and saying so now is the point of saying it at
all.** Cells grew by 57 while multi-laboratory cells grew by 6, so 51 of the 57
new cells are single-laboratory and the *share* of multi-laboratory cells fell
from 29% to 20%. Thirty more papers at this rate land at about 29, against a
target of 40. **On the evidence so far the secondary metric will miss.**

That extrapolation is recorded before the second half is curated so that neither
outcome can be rationalised afterwards. If the final figure comes in near 29,
this entry predicted it. If it comes in at 40 or above, the second thirty papers
were materially different from the first thirty and that difference is the thing
to explain.

**Why the new papers founded cells instead of joining them.** The selection took
all populations rather than the two the corpus already covers most (D-038's
reasoning applied to selection: restricting would have guaranteed the metric).
The 30 papers brought orexin, CRH, VIP, SF1, MCH, LPO and GHR populations, and a
cell needs a *second laboratory on the same population, perturbation and
outcome*. Breadth is what this corpus is short of and also what makes replication
rarer, and both halves of that are the finding rather than a flaw in the pass.

**The primary metric was not measured, and `reports/baseline_benchmark.md` is
deliberately stale.** It still describes the 362-record corpus. D-048 fixes one
evaluation at 60 papers, and regenerating the benchmark now would be the second
look that turns a suggestive result significant by chance. A stale report in a
repository that requires recomputable numbers is a cost; taking it knowingly,
and writing down why, is cheaper than the alternative.

---

## D-052 - The one new contradiction pools two genotypes, and is left standing

**What appeared.** Admitting the 30 circuit papers took cross-study conflicting
claims from 4 to 6 and produced the corpus's first directly-opposed pair since
the condition fill removed the last one (D-047):

| | direction | study | lab |
|---|---|---|---|
| `0c27e142` | decrease | `molmet.2022.101455` | Herzog |
| `f2ba156b` | no_change | `molmet.2022.101455` | Herzog |
| `ad4a6ccf` | **increase** | `s41467-020-14291-3` | Brüning |

Cell: `arc_agrp`, chemogenetic activation, `locomotor_activity`, every condition
slot unknown.

**It is partly real and partly an artefact, and the artefact is nameable.** Both
papers manipulate AgRP neurons *lacking NPY*, and `arc_agrp` does not record
that. Herzog's record is explicitly "AGRP neurons that lack Npy"; Brüning's
quote covers `NPY^wt/wt` and `NPY^Δ/Δ` together. So the cell pools an
NPY-intact population with an NPY-null one — the same defect class as the
projection pooling fixed in D-050, one level further in: **a key too coarse to
tell two experiments apart manufactures a contradiction.** The curator who
curated it said so at the time, and recommended a genotype dimension rather than
a `population_raw` string.

**Not fixed, deliberately.** Adding genotype to `population_key` would split
cells across the whole corpus and move D-048's pre-registered baseline, which is
the one thing a mid-experiment change must not do (D-050 explains why that
mattered there, and the same argument binds here). The conflict therefore stands
as recorded, which is also what SCHEMA.md requires: conflicts are data, nothing
is resolved, and both sides keep their quote. The caveat is written here so that
a reader meeting this cell in a QC report knows what it is.

**A second, smaller miss in the same cell, and the system behaved correctly.**
The two Herzog records come from one sentence — "a significant reduction in
locomotor activity in the fed but not [fasted]" — and are stored with
`challenge` null on both, so a fed result and a fasted result share a cell and
read as a within-study disagreement. The condition dimension exists for exactly
this. The v0.3.0 extractor did not fill them, and was right not to: a quote
naming two states triggers
`refused_quote_names_more_than_one_state` rather than picking one (D-044). The
fix is a human reading those two records, and the refusal is what makes the
small set findable.

---

## D-053 - A fast-refeed intake test is `refed`, and six records disagree

**The convention, decided here because the corpus has been inconsistent about it
since the condition fill.** `challenge` is the acute physiological state **at the
moment of measurement** (SCHEMA.md 3.3b). So:

* `fasted` — the measurement was taken while the animal was still without food:
  fasting glucose, fasting-induced lipolysis, energy expenditure in the fasted
  state, plasma hormones in a fasted animal.
* `refed` — food was returned and the measurement is of what happened next. **A
  food-intake measurement after a fast is this**, whatever the paper calls it:
  the animal is eating, and "the state at the moment of measurement" is
  refeeding.

**Why it needed deciding.** A curator asked, and the corpus answered both ways.
Eleven food-intake records carry `fasted` or `refed`, and the same design —
deprive, return food, measure intake — appears under both labels:

| labelled `fasted` | labelled `refed` |
|---|---|
| "food intake following an overnight fast" | "food intake upon refeeding" |
| "food intake after a fast" | "subsequent food intake ... in a refeeding paradigm" |
| "feeding after a fast" | "food intake", detail "concurrent with re-feeding after a 24 h fast" |
| "fasted-induced hyperphagia" | |
| "Stimulation of Arc POMC neurons in fasted mice ... reduction in food intake" | |
| "activation of SF1 neurons in fasted mice suppressed feeding" | |

Those last two read as "measured in a fasted animal", but an animal whose food
intake is being measured has food. They are the same fast-refeed test as the
right-hand column.

**The cost of leaving it.** This is the condition dimension manufacturing
non-comparability out of vocabulary rather than biology: a fast-refeed intake
test from one laboratory will not match the same test from another because one
curator wrote `fasted` and the other `refed`. That is the exact failure the
dimension was added to prevent (D-023), arriving by the back door.

**One probable misread found in the same pass**, recorded rather than quietly
fixed: `c1b98d981dba36bb` carries `challenge: fasted` from "AgRP chemogenetic
stimulation led to voracious total food intake, **comparable to fasted mice**".
"Fasted mice" there is the comparator the effect is measured against, not the
state of the animals measured. The condition should probably be null.

**Deferred, and why.** Six records need `fasted` changed to `refed` and one
needs a condition removed. `challenge` is a `record_id` component, so this is a
re-identification, and `fill_conditions` fills what is missing — it has no
corrector, by design, because overwriting a curated value is the thing it
refuses to do (D-039). Building one is the follow-up.

> **Built and applied in D-058.** All seven corrected, 0 refused. The
> inconsistency turned out to be inventing a replication rather than hiding
> one: the only claim that changed status was a two-study cell in which both
> records carried a `fasted` label neither deserved.

It is deliberately **not** done during the D-048 experiment. It does not touch
that measurement — `cell_reach` keys on population, perturbation and outcome and
carries no condition slot — but it does move the claim-key conflict and
replication counts, and changing those while a pre-registered pass is half
curated is how a result becomes unfalsifiable. The convention is fixed now so
that the second thirty papers are curated against it; the existing six are
corrected afterwards.

---

## D-054 - Readouts the circuit pass had to decline for want of an id, listed not added

**Decision.** The Tier A outcome vocabulary is **not** extended during this pass.
The readouts below were genuine interventional observations that three readers
declined because `outcome_id` had nowhere to put them. They are recorded here so
the next pass can extend the vocabulary deliberately, which is this project's
rule for an unmappable value (SCHEMA.md 3.2 applies the same discipline to
populations).

**Physiological readouts with no Tier A id:**

| readout | paper | note |
|---|---|---|
| bone mineral **content** | `molmet.2023.101790` | `bone_density` exists with density-only units. One sentence reports BMC up and BMD unchanged; mapping both to `bone_density` would put an increase and a no-change in one cell, so only the BMD null was recorded |
| skeletal-muscle thermogenesis | `biom14070821` | `bat_thermogenesis` is BAT-specific |
| serum CRH | `elife.90191` | `corticosterone` is not CRH |
| leptin sensitivity | `neuron.2020.04.020` | `insulin_sensitivity` has no leptin analogue. Would also need a convention for *interaction* records: the finding is that a manipulation changes leptin's effect |
| median eminence barrier permeability | `neuron.2020.04.020` | |
| skin / back temperature, and the BAT-minus-back differential | `molmet.2023.101790` | the differential is that paper's primary thermogenesis statistic |
| mineral apposition rate | `molmet.2023.101790` | |
| progressive-ratio break point | `nature11270` | motivation to work for food |

**Populations with no id**, kept unmapped with `population_raw` filled, as
SCHEMA.md 3.2 requires: LPO Vgat (lateral preoptic — `region_id` resolves to
`lpo`, there is no population entry), zebrafish Hcrt (the paper places them in
the hypothalamus, but `lha_orexin` is an LHA population a zebrafish does not
have), ARC GHR+ (the paper's own point is that GHR and GHRH neurons dissociate),
and the intersectional Arc Lepr/Npy subsets.

**Deliberate exclusions, not gaps.** Behavioural-affective readouts — anxiety,
impulsivity, novelty-seeking, alcohol intake, flavour preference, operant lever
pressing, open-field centre time — map to nothing because Tier A is systemic
*physiology*. Listing them here stops the next reader proposing them as
vocabulary gaps.

**Why not just add the ids.** Each is a scope decision about what this dataset
is for, and several carry a second question with them: BMC needs units and a
relationship to BMD, leptin sensitivity needs an interaction convention the
schema does not have, and the BAT-minus-back differential is a derived quantity
rather than a measurement. Adding eight ids in the middle of a pre-registered
pass, on the judgement of the agent doing the pass, is how a vocabulary becomes
whatever the last curator wanted.

---

## D-055 - The D-048 evaluation: primary established, secondary missed, and a baseline I should have re-measured

**Run once, at 60 papers, as pre-registered.** 146 records admitted from the 60
papers (96 from the first thirty, 50 from the second), 0 rejected by the gate.
Corpus: **508 records, 75 studies**.

### Primary: MET

Circuit subset, `perturbation_axis_cell_majority` against `global_majority`,
leave-one-study-out, rows whose cell was populated in training:

| corpus | code | effect | accuracy | n / groups | McNemar | verdict |
|---|---|---|---|---|---|---|
| 362 records | pre-D-050 key | +0.121 [+0.008, +0.271] | 0.517 vs 0.397 | 116 / 20 | p=0.156 | suggestive |
| 362 records | **final code** | +0.097 [−0.019, +0.245] | 0.504 vs 0.407 | 113 / 20 | p=0.266 | **no_difference** |
| **508 records** | **final code** | **+0.208 [+0.089, +0.354]** | **0.558 vs 0.351** | **154 / 35** | **p=0.004** | **established** |

The target was absolute — "reaches verdict `established`" — and it is reached:
the study-clustered bootstrap excludes zero and the row-level exact McNemar test
gives p=0.004, which is the first time those two have agreed on anything in this
project.

**And a mistake of mine, which the middle row exists to show.** D-050 claimed the
projection-key fix "does not move the pre-registered target". I verified that for
the *secondary* metric — 17 multi-lab cells before and after — and did not check
the primary. The primary baseline did move, from `suggestive` to
`no_difference`. So the figure D-048 recorded as the baseline was stale from the
moment I fixed the key, and I did not notice for a day.

It does not change whether the target was met, because the target names an
absolute verdict rather than a difference from baseline. It does change the size
of the claimed improvement, in my favour, which is exactly why it has to be in
the table: a control that flatters the result is still a control, and the way to
keep it honest is to publish it next to the number it flatters. Had I quoted
"+0.121 → +0.208" without re-measuring, the comparison would have mixed two
versions of the code.

### Secondary: MISSED

**26** circuit cells span two or more laboratories, against a target of **40**.
The corpus went 59 cells / 17 multi-lab → 116 / 23 → **156 / 26**.

**D-051 called this.** At the 30-paper checkpoint I wrote that the rate implied
"near 29, against a target of 40" and that "on the evidence so far the secondary
metric will miss". It landed at 26.

### What the split result means, which is the part worth reading

The primary improved substantially while the secondary missed badly, and those
are not in tension once the two are stated precisely. **The signal did not come
from broad cross-laboratory replication. It came from more rows landing in the
few cells that already recur.** 26 cells carry the 154 rows the primary is
measured on; the other 130 circuit cells hold one laboratory each.

So the mechanism in D-048's hypothesis — "circuit manipulations recur across
laboratories, so adding circuit papers adds replicated cells" — is only weakly
supported. Adding 60 papers added 97 cells and 9 multi-laboratory ones. What it
did add is *depth* in the handful of populations the field keeps returning to,
and depth is what the direction model needed.

**The result therefore rests on a narrow base and the next finding to chase is
whether it generalises.** 0.558 on a three-way problem over 154 rows in 35
studies is established, not accurate, and it is established on the subset where
a cell was seen — a third of circuit rows still have no evidence in their cell
at all.

> **Chased, in D-056 (pre-registration) and D-057 (result).** It generalises
> across laboratories — `established` at the same +0.208 under a lab grouping
> that, until D-056, had silently been identical to `study_id`. It does not
> generalise across cell depth: the whole effect lives in the 88 rows whose
> cell two or more laboratories had already reached, while the 66 rows with one
> prior laboratory give +0.106 spanning zero. The "26 cells carry the 154 rows"
> sentence below is now computed by `cell_reach` rather than by hand.

### Also, and separately: the first feasible task

`syn4rch.feasibility` now reports **1 of 180 (population, outcome) tasks
feasible**, where every previous run said none. It is `arc_agrp` / circuit
activation / `food_intake`: 3 studies, 3 laboratories, agreement 1.00, clearing
every gate.

Two things keep that in proportion. It is the most reproduced result in this
literature — AgRP activation makes a mouse eat — so finding it is a sanity check
on the pipeline rather than a discovery. And both of its claim-keyed variants
carry `CONDITIONS_ALL_UNKNOWN` or `CONDITIONS_PARTIALLY_UNKNOWN`: they matched on
0 and 1 of 4 stated condition slots, so even the feasible task is a match on
mutual ignorance about three of its four conditions.

### No abandon condition fired

D-048 named three. Multi-lab cells were not flat at 30 papers. The primary was
not reached after the stopping point — it was reached at it. And the one that
would have been damning, "the secondary reached while the primary is not", is the
mirror image of what happened: the metric that would have said the signal was a
116-row artefact is the one that came in strong.

---

## D-056 - Pre-registration: does the established direction result generalise, and a lab identity that has never worked

**Written before any of it is measured**, and committed before the code that
measures it, for the same reason D-048 was: D-055 ends by saying the result
"rests on a narrow base and the next finding to chase is whether it
generalises", and a generalisation test designed after seeing its own answer is
not a test.

**This does not reopen D-048.** That pre-registration named leave-one-study-out
and is closed; its primary was met and stays reported as met. What follows is a
*harder* question asked of the same result, and nothing here can retroactively
change whether D-048's target was reached.

### The defect found while setting the test up

`lab_key` is `(last_author, last_author_affiliation)`. Affiliation is free text
harvested from PubMed, and the same laboratory writes it differently in every
paper. The corpus has **75 studies and 75 "labs"**, so leave-one-lab-out has
been numerically identical to leave-one-study-out in every report this project
has published — and SCHEMA.md 5 item 2, which names lab identity "the largest
leakage risk in this corpus", has therefore never been controlled for.

Four last authors contributed more than one study. All four are unmistakably one
laboratory each:

| last author | studies | what splits them |
|---|---|---|
| Brüning JC. | 5 | "max planck institute for metabolism research, department of neuronal control of metabolism, gleueler str. 50" / "...gleueler strasse 50" / department and institute in the other order / email appended / "neuronal control of metabolism group" |
| Herzog H. | 2 | "neuroscience division, garvan institute of medical research, st vincent's hospital, darlinghurst, nsw 2010" vs "garvan institute of medical research, nsw" |
| Sweeney P. | 2 | field order, and one carries `sweenp@illinois.edu` |
| Withers DJ. | 2 | "mrc london institute of medical sciences; institute of clinical sciences, imperial college london" vs "institute of clinical sciences, imperial college" |

Every difference is typographic. **121 of the 314 circuit records — 39% — belong
to those four laboratories**, which is precisely the fraction on which the
leakage control was supposed to bite.

### The correction: a laboratory is a person

`lab_key` becomes the normalised `last_author` alone. The affiliation adds no
discrimination that can be trusted and in this corpus it only ever splits a lab
into pieces, never joins two.

The failure mode of the looser key is a **homonym**: two unrelated researchers
whose normalised surname and initials coincide would merge into one lab, which
costs test rows and manufactures a replication. That is made detectable rather
than assumed away — a new warning fires when one last author appears under
affiliations sharing no substantive token. All four labs above share a strong
token (`max planck institute for metabolism research`, `garvan institute of
medical research`, `university of illinois urbana-champaign`,
`institute of clinical sciences`), so the warning is silent on this corpus and
would fire on a genuine collision.

`record_id` does not hash the author or the affiliation (`RECORD_ID_FIELDS`), so
this re-identifies nothing.

### The six predictions

Each is a number, written now:

1. **Grouping.** Lab groups fall 75 → **68**, with 0 records ungroupable.
2. **The primary, cross-lab.** The headline comparison
   (`perturbation_axis_cell_majority` vs `global_majority`, circuit subset, rows
   whose cell was populated in training) under the **corrected `lab` grouping**
   keeps verdict `established`, with delta ≥ **+0.10** — at least half of the
   +0.208 it shows under leave-one-study-out.
   *If it comes in below `established`, D-055's primary is in part one
   laboratory agreeing with itself, and the README headline and D-055 must say
   so.*
3. **The secondary, corrected.** Circuit cells spanning ≥2 laboratories fall
   from 26 to somewhere in **20–25**. D-048's secondary metric was already
   missed at 26 against a target of 40; the corrected figure is the one to
   quote, and it can only be worse.
4. **The one feasible task survives.** `arc_agrp` / circuit activation /
   `food_intake` is reported feasible on 3 studies and 3 laboratories. Predict
   its three studies have three distinct last authors and it stays feasible.
   *If two of them are one lab, the first feasible task this project ever found
   was an artefact of lab over-splitting, and feasibility returns to 0 of 180.*
5. **Depth.** Among cell-seen circuit rows, split by how many distinct
   laboratories contributed training rows to the test row's own cell: exactly 1,
   or 2 and above. Predict the delta is larger at ≥2 laboratories than at 1, and
   **positive at both**. If depth 1 is flat, the model needs a cell to have been
   *replicated*, not merely observed once, and that is a sharper statement of
   what it can do than the headline is.
6. **The verdict rule is calibrated.** Permuting labels **within** study
   clusters (preserving each study's label multiset, destroying the association
   between a cell and a direction across studies) and re-running the whole
   comparison 200 times reaches `established` in **fewer than 5%** of
   permutations.
   *If it reaches it more often, the "clustered interval and McNemar agree" rule
   manufactures significance at n=154 in 35 clusters, and every `established`
   verdict in this project is suspect — including D-055's.*

### Three things already measured, read out rather than re-run

These are in the committed `reports/baseline_benchmark.json` as of c21ce6d and
are **not** predictions. They are the existing answer to "does it generalise",
and they are unflattering:

* **Unseen population: nothing.** Under leave-one-population-out the headline
  candidate's delta against the floor is exactly **0.000** on every subset, and
  the cell-seen slice holds **0 rows**. Holding out a population empties the
  cell key of which it is a component, so the model reduces *exactly* to the
  global mode. Its only mechanism is cell lookup.
* **No transfer to the molecular half.** Under leave-one-study-out the molecular
  subset has **0** cell-seen rows and delta exactly 0.000: every molecular cell
  is a single study, so holding that study out empties it (D-027, D-048).
* **43% of the cell-seen set is within-study.** Row-level cross-validation gives
  268 cell-seen circuit rows; leave-one-study-out gives 154. The other 114 rows
  are cell-seen only via other records of their *own* study.

### What a reader should take from the shape of this entry

The generalisation question has a known negative answer on two of its three
faces before anything new is run. The open face is whether the within-cell
signal is cross-laboratory, and the reason that face was open is a key that has
silently agreed with `study_id` for the life of the project.


---

## D-057 - The generalisation test: cross-laboratory yes, unreplicated cells no, and one prediction wrong

**The D-056 probes, run once.** Five of six predictions held. The one that
failed is stated first, because a pre-registration whose misses are buried at
the bottom is a press release.

### Prediction 3 was wrong

I predicted that correcting lab identity would cut circuit cells spanning two or
more laboratories from 26 to "somewhere in 20-25", and that the corrected figure
"can only be worse". **It is 26.** Not one of the 26 multi-study circuit cells
was a laboratory agreeing with itself: the four repeat laboratories' cells either
are not multi-study cells, or carry a third laboratory as well.

So D-048's secondary metric missed at 26 against a target of 40, and it missed at
26 under both definitions of a laboratory. My reasoning — 39% of circuit records
sit in four repeat labs, therefore some of the replication must be theirs — was
an inference from the right number to the wrong conclusion, and the measurement
says so.

### Prediction 1: lab groups 75 -> 68, nothing ungroupable. Held

68 groups, 0 records ungroupable, and the largest group grows from 34 rows to 54
(Herzog's two papers). `syn4rch.qc.report` had a **fourth** copy of the lab
identity, computed inline as `(last_author, last_author_affiliation)`, and kept
printing "labs identified: 75" after the shared key was corrected. It now calls
`lab_key`, which is the duplication `qc/keys.py` exists to prevent.

### Prediction 2: the signal is cross-laboratory. Held

| grouping | delta | clustered 95% CI | accuracy | rows / groups | McNemar | verdict |
|---|---|---|---|---|---|---|
| `study_id` | +0.208 | [+0.089, +0.354] | 0.558 vs 0.351 | 154 / 35 | p=0.004 | **established** |
| `lab` (corrected) | +0.208 | [+0.086, +0.385] | 0.558 vs 0.351 | 154 / 32 | p=0.004 | **established** |

The threshold was `established` with delta >= +0.10. It is `established` at
+0.208, so D-055's primary is not one laboratory agreeing with itself.

**And the accuracies being equal to four decimals is a coincidence, not an
identity.** Merging the four repeat laboratories into single folds changed
exactly **two** predictions: one row with truth `decrease` went
`no_change` -> `decrease` (wrong to right) and one with truth `no_change` went
`no_change` -> `decrease` (right to wrong). They cancel. Reporting "the
predictions are unchanged" would have been a stronger and false claim; what
changed is the clustering, 35 groups to 32, which widens the interval.

### Prediction 5: it needs a replicated cell, not merely a seen one. Held, with a caveat on my own wording

Circuit rows by how many distinct **laboratories** had reached the test row's own
cell in training:

| stratum | rows | delta | clustered 95% CI | accuracy | McNemar | verdict |
|---|---|---|---|---|---|---|
| `depth_0` | 159 | **+0.000** | [+0.000, +0.000] | 0.459 vs 0.459 | n/a | no_difference |
| `depth_1` | 66 | +0.106 | [−0.113, +0.421] | 0.485 vs 0.379 | p=0.349 | no_difference |
| `depth_2_plus` | 88 | **+0.284** | [+0.108, +0.484] | 0.614 vs 0.330 | p=0.006 | **established** |

`depth_0` is exactly 0.000 with zero discordant rows, which is the instrument
check: where the cell was never populated the candidate *is* the floor, and a
non-zero value there would have meant this baseline is not what the project says
it is. It is reported rather than assumed.

I predicted the delta would be "larger at >=2 laboratories than at 1, and
positive at both". Literally that held: +0.284 against +0.106, both positive. But
+0.106 on 66 rows in 18 groups spans zero in both tests, so the honest reading is
**one prior laboratory in a cell buys nothing this corpus can detect, and the
whole of the established effect lives in the 88 rows whose cell two or more
laboratories had reached.** That is a sharper and smaller claim than the
headline, and it is the claim the system can actually support.

### Prediction 6: the verdict rule is calibrated. Held, and the pre-registered 200 permutations were not enough to say so

Labels permuted **within** study, so every study keeps its own label multiset and
the cell structure is untouched; what is destroyed is the association between a
cell and a direction across studies.

| permutations | reached `established` | rate | null deltas >= observed | empirical p |
|---|---|---|---|---|
| 200 (pre-registered) | 6 | 3.0% | 0 of 200 | 0.005 |
| 2000 (precision) | 43 | 2.15% | 0 of 2000 | 0.0005 |

The rule fires on null data 2.15% of the time, under the 5% threshold. **But
6/200 has a binomial 95% upper bound of about 6.4%, which crosses the very
threshold the probe was testing** — the pre-registered permutation count was too
small to answer its own question, and I only noticed when the answer was close
to the line. The 2000-permutation run is more of the same quantity at the same
threshold, not a new metric, and it settles it at about [1.6%, 2.9%].

**The null is not centred on zero, and that is the most useful thing this probe
produced.** Its median delta is +0.039 and its 95th percentile +0.104, because
permuting within study preserves each study's label multiset: a cell one study
dominates keeps that study's labels, so cell lookup still recovers study-level
label structure. Anyone reading "delta > 0" as evidence that a cell carries
direction information is comparing against the wrong origin. The observed +0.208
exceeds every one of 2000 permuted deltas (max +0.195).

### Prediction 4: the one feasible task survives. Held

`arc_agrp` / circuit activation / `food_intake` is still 3 studies, **3
laboratories**, agreement 1.00, both claim-keyed variants still passing every
gate. Its three studies have three distinct last authors, so the first feasible
task this project ever found is not an artefact of lab over-splitting. The cell
it sits in did lose laboratories: `arc_agrp` / `food_intake` goes from 23
studies / 23 labs to 23 studies / **20** labs.

### What the three probes say together

The result generalises **across laboratories** and not **across cells**. Those
are different axes and only the first was open:

* across laboratories: yes, with the same effect size and a wider interval;
* across cell depth: no — the effect is entirely in cells two or more
  laboratories had already reached;
* to an unseen population: no, and exactly so. D-056 read this out of the
  committed benchmark rather than re-running it: the candidate's delta against
  the floor is **0.000** under leave-one-population-out. The population is a
  component of the cell, so holding it out empties the cell and the model
  reduces to the global mode. Its only mechanism is cell lookup;
* to the molecular half: no, for the same reason at a different key — 0 of 170
  molecular rows have a populated cell under leave-one-study-out.

So the system predicts the direction of a circuit manipulation when at least two
laboratories have already published that exact manipulation and readout. That is
a real and measurable capability, it is not a model of the hypothalamus, and the
distance between those two sentences is the project's position.

### Three defects found while running it, each fixed as a class

**1. `n_records_in_multi_group_cells` was declared, serialised and never
populated.** It defaulted to 0 and `cell_reach` never passed it, so anything
reading the identity JSON was told that **0** records lived in the 26 multi-group
circuit cells. It now reports **154** — which independently confirms D-055's
assertion that "26 cells carry the 154 rows the primary is measured on". That
sentence was true when I wrote it and I had computed it by hand; the code now
computes it.

**2. `W010_NO_LAB_IDENTITY` contradicted the key it exists to protect.** It
fired when the author *or* the affiliation was blank, so after D-056 made lab
identity the author alone, **17 records that group perfectly well were being told
they have no lab identity**. W010 now tracks `lab_key` exactly, and a new
`W018_NO_LAB_AFFILIATION` covers the real cost of a missing affiliation: the
homonym check cannot look at that author. Without the second warning, a corpus
with no affiliations at all would leave `QC_POSSIBLE_LAB_HOMONYM` silently
vacuous, which is the failure D-056 was written about. Total warnings are
unchanged at 2714: the same 17 records, under the code that describes them.

**3. `find_possible_non_independence` printed its claim components under the
wrong names.** `claim_key` stopped carrying `valence` when v0.2.0 added the
perturbation axis, and the message had since been labelling the perturbation
class "valence", the polarity "outcome" and the target "direction". Nothing
noticed because the check never fired. The components are now named from
`CLAIM_KEY_FIELDS`.

### And a correction to D-056's own framing

D-056 said the non-independence detector "returned nothing on a corpus in which
four labs had published eleven papers", implying the affiliation was why. The
cause was right and the implication was wrong: **after the fix it still returns
nothing.** `claim_key` has nine components including all four conditions, and no
laboratory in this corpus repeats a claim that specific. Under the coarser
`(population, outcome, direction)` there are six lab-level repeats, among them
Brüning reporting `arc_agrp` -> `food_intake` increase in three separate studies.

So the detector keys on something **finer than the thing labs actually repeat**,
and finer than the cell the model learns on. That is a real gap and it is left
open deliberately: the claim key is the right granularity for conflict and
replication, the cell is the right granularity for leakage, and they are two
questions. The leakage question is now answered directly by probe 1 rather than
inferred from a detector — which is the better answer anyway.

### What moved in the committed reports

Every `by_lab` number in `reports/baseline_benchmark.md` changed, and every
`study_id` number is bit-identical to the version committed at c21ce6d. That
equality is also the control on a refactor made for this work:
`evaluate.paired_slice` was extracted from a closure so the depth strata could
use the same interval and the same test as the headline slice, and the way to
show an extraction changed nothing is to show the numbers it does not touch
standing still.

One substantive change beyond the group counts: `neurotransmitter_sign` under
`by_lab` falls from `established` to `suggestive`. It is *worse* than the floor
by 0.162, and the corrected grouping makes that harm less certain, not more.


---

## D-058 - The D-053 corrections applied, and the inconsistency was inventing replication rather than hiding it

**The follow-up D-053 deferred.** `fill_conditions` fills an empty condition
dimension and refuses to overwrite one that is not empty (D-039), so the six
records storing `fasted` for a fast-refeed intake test, and the one whose
`challenge` came from a comparator, could not be fixed by any tool in the
project. `syn4rch.curate.correct_conditions` is that tool.

### The corrector decides nothing

Every change is a row of `CORRECTIONS`: the record, the dimension, **the value
that must currently be stored**, the value to put there, a fragment of the
stored evidence quote that must match, the reason, and the decision that
authorises it. A correction the table does not contain cannot be applied, and
the table is reviewed in a diff rather than inferred at runtime.

Two locks have to open, and the second one matters more than it looks:

* `from_value` is not a wildcard. A corrector that overwrites whatever it finds
  cannot be reviewed, and overwriting a curated value is precisely the thing
  D-039 refuses to do casually.
* `quote_contains` addresses the correction to the **evidence**, not only to the
  record. A `record_id` hashes twelve fields and not the quote, so it survives
  re-curation of the sentence behind a value; without this check a correction
  written about one sentence could silently apply to another.

After the corrections, every condition quote in the whole dataset is re-checked
against the passage its locator names — the locator-strict check from the fill,
not the weaker "occurs somewhere in the paper". Any failure, any refusal, or any
hard schema violation blocks the write.

A correction keeps its quote, locator and source: the same sentence is still the
evidence, only how it is read has changed. Only the value moves, and `rule`
becomes `corrected_by_decision_d053` so the entry says a convention decided it
rather than claiming an extractor read it that way. Where a value is **removed**
the evidence entry is deleted, because an entry asserting a value that is no
longer stored is what A009 and A011 refuse in both directions.

### What was applied

7 of 7, 0 refused. `challenge` is a `RECORD_ID_FIELDS` component, so all seven
were re-identified; the id map is complete over the input, `conflicts_with` was
rewritten through it, and a collision would have refused the run.

| old id | new id | from | to |
|---|---|---|---|
| `0397b7caf40ad06b` | `cf668711e4369657` | fasted | refed |
| `087960c4021bf6f2` | `e527341d2a14a370` | fasted | refed |
| `41fca07890c238c8` | `d206091ce714a8f6` | fasted | refed |
| `fcea91b0669ad33f` | `bf6a25067c559d29` | fasted | refed |
| `93c2e8f0c13a30ed` | `e1e102a40065cdef` | fasted | refed |
| `fdb78695221a531e` | `3f9e8f6e7b671752` | fasted | refed |
| `c1b98d981dba36bb` | `1068b48a9cf28019` | fasted | **null** |

`challenge` goes `fasted` 28 → 21, `refed` 8 → 14, null 443 → 444. Diffing the
whole dataset through the id map, the only fields that changed anywhere are
`challenge` and `condition_evidence`, on exactly 7 records. No direction moved.

### D-053's reasoning was right about the mechanism and wrong about its sign

D-053 argued the inconsistency was "the condition dimension manufacturing
non-comparability out of vocabulary rather than biology: a fast-refeed intake
test from one laboratory will not match the same test from another because one
curator wrote `fasted` and the other `refed`." The expectation was therefore that
correcting it would **recover** matches.

It did the opposite, and the opposite is the better outcome. Claims reported by
more than one study go **10 → 9** and cross-study conflicts go **9 → 8** — and
both are the *same pair of records*:

| | direction | study | quote |
|---|---|---|---|
| `c1b98d981dba36bb` | increase | `molmet.2024.102021` | "voracious total food intake, **comparable to fasted mice**" |
| `fcea91b0669ad33f` | no_change | `molmet.2022.101455` | "**food intake after a fast**" |

Cell: `arc_agrp`, circuit activation, `food_intake`, `challenge: fasted`,
everything else unknown. That cell was simultaneously counted as a replicated
claim and reported as a cross-study disagreement — which is what a claim key
does when two non-comparable experiments land in it. **Neither record deserved
its `fasted` label.** One is a fast-refeed intake test; the other took the label
from the comparator the effect was measured against. Correcting both dissolved a
false replication and a false contradiction in one move, about the single most
reproduced result in this literature.

So in this corpus the vocabulary inconsistency was **inventing comparability**,
not suppressing it. The mechanism D-053 named is real — a key too coarse or too
loosely applied to tell two experiments apart (D-050, D-052) — and its sign here
runs the other way.

**The six `refed` corrections created no new matches at all.** Each landed in a
cell with no cell-mates: there is no second laboratory in this corpus running the
same fast-refeed intake test with the other three conditions matching. The
convention is still worth having for the records curated after it (D-053 fixed it
before the second thirty papers for that reason); it bought no replication here.

### And the one feasible task got stronger

The nulled record moved into the all-unknown `arc_agrp` / activation /
`food_intake` cell, which is the feasible task of D-055. It goes from 4 records /
3 studies / 3 labs to **5 records / 4 studies / 4 labs**, agreement still 1.00
with the interval tightening from [0.44, 1.00] to [0.51, 1.00]. The spurious
`fasted` variant — 2 studies, agreement 0.50, tied modal direction, carrying
`QC_CONFLICT_EFFECT_VS_NULL` — is gone from the report entirely.

That is the shape of this pass: the replication *count* fell by one and the
replication *evidence* improved, because the thing removed was noise and the
record removed from it went where it belonged.

### What did not move, as D-053 predicted

The direction-of-effect result is untouched: +0.208 [+0.089, +0.354], p=0.004,
`established`, 154 rows in 35 studies, and identically under the corrected lab
grouping. D-053 said so in advance — "it does not touch that measurement:
`cell_reach` keys on population, perturbation and outcome and carries no
condition slot" — and that prediction held exactly, which is why the corrections
were safe to make after the pre-registered evaluation closed rather than before.

The one benchmark number that did move is the condition-conditioned baseline,
which is the only one that reads `challenge`: `perturbation_axis_condition_majority`
against the floor goes from +0.361 on 61 rows to **+0.400 on 60 rows**, still
`established`, and still `no_difference` against the unconditioned axis baseline.
Conditioning continues to fragment more than it helps (D-036), on one row fewer.

### The cost of this pass, stated

A re-identification is not free. Seven ids in the published record are now
superseded, and `data/curated/record_id_map_d053.json` is the third such map in
the project's history (0.1.0 → 0.2.0, 0.2.0 → 0.3.0, and this one). Anything
citing one of those seven ids — including QC output written before this pass —
must be resolved through the map. The alternative was leaving seven records
knowingly wrong, which is worse, but the map is a permanent obligation rather
than a one-off.


---

## D-059 - Pre-registration: can transcriptomic similarity predict an unseen population's direction?

**Written before the atlas is opened**, and committed before the code that uses
it. D-057 measured the sharpest failure this project has: under
leave-one-population-out the headline candidate's delta against the floor is
**exactly +0.0000** on every subset, because the population is a component of
the cell key and holding it out empties the cell. The model's only mechanism is
cell lookup, so on a population it has never seen it reduces to the global mode
by construction. This entry fixes what "transcriptomic similarity fixes that"
would have to mean.

### The honest bar is not zero

Beating +0.0000 would be trivial, and the benchmark already shows why. Under the
same leave-one-population-out grouping, on the circuit subset:

| comparison | delta | verdict |
|---|---|---|
| `perturbation_axis_cell_majority` vs `global_majority` | **+0.0000** | no_difference |
| `perturbation_axis_majority` vs `global_majority` | +0.0543 | no_difference |
| **`outcome_majority` vs `global_majority`** | **+0.0735** | **established** |

`outcome_majority` ignores the population entirely and still generalises,
because `outcome_id` is not what the grouping holds out: the base rate of a
readout is portable. So any population-similarity model must beat **the baseline
that does not look at the population at all**. That is the comparison this
entry pre-registers, and choosing the weaker one would have manufactured a
result.

### The data, named before it is used

**HypoMap** (Steuernagel et al., *Nature Metabolism* 2022,
`doi:10.1038/s42255-022-00657-y`), the unified single-cell atlas of the murine
hypothalamus: 384,925 cells, one 3.8 GB H5AD from CELLxGENE collection
`d86517f0-fa7e-4266-b82e-a521350d6d36`. Downloaded to `data/raw/hypomap/`, which
is untracked machine input (D-016).

It is the right atlas because it is hypothalamus-specific, mouse, integrated
across the same literature this corpus is drawn from, and its clusters are
annotated with the marker genes our `data/vocab/populations.json` already uses.

**Why the cheap version was abandoned, with the measurement that killed it.**
The obvious first attempt is similarity from the marker genes already in the
vocabulary. It cannot work: markers are *chosen* to be unique to a population.
Across the 30 vocabulary entries there are 38 marker slots holding 33 distinct
genes, and **6 of 435 population pairs share any marker at all** — `Kiss1` in
ARC and AVPV, `Bdnf` in DMH/VMH/POA, `Avp` in PVH and SCN, `Slc32a1` in LHA and
ZI. A similarity matrix built from them is the identity matrix, which carries no
information about which populations resemble which. Real expression profiles over
thousands of genes are the only version of this that can carry signal.

### The feature, fixed now

1. Per-cluster mean log-normalised expression from HypoMap, over a gene set
   declared in code (not chosen after seeing results).
2. A **mapping** from this project's `population_id` values to HypoMap clusters,
   curated with the cluster's own annotation and marker genes as the evidence,
   recorded as data and reviewable in a diff — the `CORRECTIONS`-table pattern
   from D-058, for the same reason: a mapping inferred at runtime cannot be
   reviewed.
3. Cosine similarity between population mean-expression vectors.
4. A new baseline, `population_similarity_majority`: for a test row whose own
   `(population, axis, outcome)` cell is empty, take the modal direction of the
   same `(axis, outcome)` cell in the **k transcriptionally nearest populations
   that were seen in training**, weighted by similarity. Where no neighbour has
   that cell, fall back exactly as the existing baselines do, so the comparison
   isolates the similarity step and nothing else.

### The predictions

1. **Primary.** `population_similarity_majority` against **`outcome_majority`**,
   grouping `population`, circuit subset, reaches delta > 0 with verdict
   `established` — both the population-clustered bootstrap interval excluding
   zero and the row-level exact McNemar test under 0.05.
2. **Secondary.** It does not damage the result D-055 established: against
   `global_majority` under leave-one-**study**-out on the cell-seen circuit
   rows, it stays `established` and its delta stays at or above +0.15 (the
   current figure is +0.208).
3. **Mapping coverage.** At least **20** of the 44 `population_id` values in the
   corpus map to a HypoMap cluster. Below that the experiment does not run, and
   the finding is that this corpus's populations are not the atlas's.

> **Prediction 3 was void as written, and is withdrawn as a gate rather than
> repaired into one (amended before the primary was measured; the amendment was
> raised independently by a second session reading this entry).**
>
> There are not 44 `population_id` values. 44 is the count of `population_key`s,
> which is a different thing: 32 soma keys (13 resolved, 19 `raw:`) plus 12
> projection keys (8 on a resolved soma, 4 on a `raw:` one). The corpus holds
> **13 distinct `population_id` values**, on 380 of 508 records. The other 128
> records carry no id at all, so their 23 keys have nothing to map *from*
> (D-010) — they are a stated ceiling on coverage, not a shortfall. Mapping all
> 13 would therefore still leave 13 < 20 and the experiment could never have run
> by its own terms. (An earlier breakdown of mine in this session said
> "13 + 12 + 23" of 44, which sums to 48: that classifier double-counted the 4
> keys that are both projection-targeted and unresolved. The split above does
> not.)
>
> **Why it is withdrawn instead of re-set.** The honest denominator is the 13
> resolved ids — and all 13 carry circuit records, so the circuit-subset
> denominator the primary needs is also 13. But I had already measured coverage
> against them while validating the atlas ingest, before noticing the error.
> **13 of 13 markers are present in HypoMap**, and at the C185 level 10 of 13
> have a cluster whose own name corroborates the marker (`arc_agrp` ->
> `Npy.Agrp.GABA-4`, `arc_pomc` -> `Anxa2.Pomc.GLU-5`, `lha_orexin` ->
> `Hcrt.Rfx4.GLU-4`, `scn_vip` -> `Vip.Vipr2.GABA-2`, `lha_mch` -> `Pmch.GLU-7`,
> plus `pvh_trh`, `pvh_crh`, `arc_kisspeptin`, `arc_ghrh`, `tmn_histamine`).
> Three are weak and need curated judgement rather than name matching:
> `vmh_sf1` (Nr5a1 peaks at 0.34; HypoMap is ARC/PVH/LHA-weighted and the VMH is
> underrepresented), `pvh_sim1` (Sim1 at 0.76), and `lha_vglut2` (`Slc17a6` is
> pan-glutamatergic, so "top cluster" is not a meaningful question for it).
>
> > **Correction, D-068: "the VMH is underrepresented" is false, and it was my
> > claim.** Measured at C465, the VMH holds 45 clusters and 33,694 cells, 8.8%
> > of the atlas, against the arcuate's 45 and 37,667 (9.8%) — effectively tied.
> > The thin region among populations carrying records is the **PVH**: 16
> > clusters, 8,710 cells, 2.3%. I asserted it here without measuring it, a
> > second session took the phrase from this entry and built an argument on it,
> > and it stood for two days as the explanation for why `vmh_sf1` maps weakly.
> > It is not the explanation. With sampling retired, that refusal rests on
> > transcript abundance alone — which makes D-065's "a Cre line is not a
> > marker" a stronger claim than it was while a confound was available.
>
> A threshold whose answer I already know is not a gate, and setting one now at a
> level I know passes would be the goalpost move D-050 and D-056 both argue
> against — in the most flattering direction available. So prediction 3 becomes a
> **reported measurement**: coverage is 13/13 on markers and 10/13 on
> name-corroborated clusters, the mapping proceeds, and the honest statement of
> its limit is that **a quarter of the corpus (128 records, 23 keys) is
> unmappable by construction**. Zebrafish Hcrt, which D-054 lists as
> deliberately unmapped and which could not map to a murine atlas in any case,
> has `population_id: None` and is already outside the 13 — no explicit
> exclusion is needed.
>
> **Predictions 1 and 2 are untouched and still unmeasured.** The primary
> (similarity against `outcome_majority`, not against the floor) and the
> secondary (the established leave-one-study-out result must survive) are the
> ones that carry the experiment, and neither has been looked at.

> **Second amendment, also before the primary is measured: two per-row
> diagnostics are now REQUIRED, because without them a miss is uninterpretable.**
> Raised by a second session, and it is right. The run must log, for every test
> row, (a) how many mapped neighbour populations were available to it and
> (b) the identity and cosine similarity of the nearest neighbour that actually
> supplied its prediction. Neither is recoverable after the fact.
>
> The reason is the shape of the corpus, measured now rather than discovered
> later. Of the 314 circuit records, 245 carry a resolved `population_id`, and
> they are extremely unevenly weighted:
>
> | population | circuit records |
> |---|---|
> | `arc_agrp` | 123 |
> | `arc_pomc` | 52 |
> | `vmh_sf1` | 22 |
> | `lha_orexin` | 16 |
> | `pvh_trh` | 8 |
> | `pvh_crh` | 7 |
> | `pvh_sim1` | 4 |
> | `lha_vglut2`, `scn_vip`, `lha_mch` | 3 each |
> | `arc_kisspeptin` | 2 |
> | `arc_ghrh`, `tmn_histamine` | 1 each |
>
> **175 of those 245 — 71% — are `arc_agrp` plus `arc_pomc`**, and 6 of the 13
> ids contribute three circuit records or fewer. Three of the 13 map weakly
> (`vmh_sf1` at a 0.34 Nr5a1 peak, and it is the third-largest circuit
> population; `pvh_sim1` at 0.76; `lha_vglut2`, whose `Slc17a6` is
> pan-glutamatergic so a top-cluster question is ill-posed). So about 10
> populations map cleanly, and under leave-one-population-out a held-out
> population has **at most 9 neighbours**, several of them contributing one to
> three records.
>
> Without the two diagnostics, a miss cannot be separated from "there were never
> enough mapped populations to tell" — and that ambiguity would make the whole
> experiment unfalsifiable in the direction it is most likely to go.
>
> **It also sharpens the prediction, in a way that makes it easier to refute.**
> For a held-out `arc_agrp` the nearest mapped neighbour is almost certainly
> `arc_pomc`: same nucleus, large shared transcriptome, and opposite effects on
> food intake. That is the exact pair this entry already names as the reason
> expression similarity should fail to transfer a sign. So the primary is, in
> large part, a test of that single pair, and the per-row log will show whether
> the misses are concentrated there. If the AgRP/POMC rows miss while the rest
> hold, the finding is specific and mechanistic rather than a flat negative.

**My prediction, recorded so it can be wrong: the primary will miss.** I expect
transcriptomic similarity to be a poor predictor of the *sign* of a systemic
response. Two populations can be transcriptionally adjacent and physiologically
opposed — ARC AgRP and ARC POMC neurons are neighbours in the same nucleus,
share a great deal of their transcriptome, and drive food intake in opposite
directions. If that is right, the result is a negative finding worth having:
expression similarity does not transfer direction of effect, and the unseen-
population failure is not fixable this way.

### What would abandon it

* fewer than 20 populations mapped (above);
* the secondary violated — a model that generalises to new populations by
  giving up the one established result is not an improvement;
* any need to choose the gene set, `k`, or the mapping *after* seeing a
  delta. All three are fixed in code before the first evaluation, and the
  evaluation is run **once**.

### What this is not

It is not a claim that expression predicts physiology, and nothing here becomes
a training target in the sense SCHEMA.md 6 forbids. It is a **feature** derived
from an external observational atlas, used to generalise a direction label that
still comes only from interventional experiments. If it works it is a
computational hypothesis generator; if it fails that is a measurement of how far
transcriptomic identity carries, and either way it is reported separately from
the demonstrated relationships in the corpus.


---

## D-060 - The genotype dimension: a fifth condition, not a population suffix, and the two risks were the wrong way round

**A measured design decision, not a pre-registration.** The quantity that
decides between the options — what each one costs the established result — was
*measured* for both before choosing, which is the standard D-050 set: "a key
change that moved the baseline would have been goalpost-moving however well
argued, and the way to tell the difference is to measure both and publish the
pair." Nothing is predicted here. D-052 declined this fix because D-048 was open;
D-048 closed at D-055, so the reason no longer binds.

### What the corpus actually contains

Seventy-five of 75 papers read for genetic manipulations beyond the Cre driver
and the actuator, every entry carrying a verbatim quote and a locator. Then every
claimed *background* genotype re-checked by an independent skeptic instructed to
refute it.

| | n |
|---|---|
| manipulations found, all roles | 665 |
| claimed as background (not the recorded perturbation) | 26 |
| **CONFIRMED on adversarial re-check** | **12** |
| refuted | 13 |
| uncertain | 1 |
| **quote failures** | **0** |

Zero quote failures: every surveyed quote was verbatim at the passage its locator
named. The survey's *citation* was perfect and its *judgement* over-reported
background by 50%. Three of the thirteen refutations, to show they are
substantive rather than pedantic:

* `Sh2b1` — the gene that paper is *testing*. Its "global knockout" appears
  experimentally in exactly one place, as an **antibody negative control**; the
  quoted sentence itself contains the words "(negative control)".
* `Prkce` — likewise the tested gene. The claimed "PKCεfl/Δ background" is the
  allelic configuration *of the perturbation*, plus a germline-deletion breeding
  artefact. The skeptic also found that paper's Methods and Results contradicting
  each other on the experimental genotype.
* `Lepr` (db/db) — a separately purchased cohort carrying no Cre and no Sirt1
  allele, used for one descriptive NAD+ readout. "Recording it as background
  would tell an auditor that some Sirt1-overexpression arm was run on a
  leptin-receptor-null background, which the paper never did."

**And the survey asked the question at the wrong granularity, which is the
finding that settles the design.** `is_the_recorded_perturbation` was asked once
per *paper*. Herzog `doi:10.1016/j.molmet.2022.101455`, 23 records, splits exactly
along it:

| records | method | class | `target_raw` | `population_raw` |
|---|---|---|---|---|
| 7 | `conditional_knockout` | molecular | **`Npy`** | "AGRP neurons" |
| 16 | `chemogenetic` | circuit_activity | `null` | **"Npy-depleted AGRP neurons"** / "AGRP neurons that lack Npy" |

Within one study the same allele is the recorded perturbation for 7 records and
unrecorded background for 16. Having nowhere to put the background, the curator
put it in `population_raw` as free text — which `population_key` then discards,
which is precisely how D-052's contradiction was manufactured. So the dimension
must be **per record**, and it must distinguish "perturbation on gene X" from
"background carrying gene X". The establishable set is **71 records in 8
studies**, not the 157 in 16 the raw survey implied.

### The decision

A **fifth condition dimension**, `background_genotype`, a controlled scalar
holding `gene_id:effect:scope` segments joined by `+` and sorted, resolved
against `data/vocab/targets.json` — whose `id_policy` already demands the
lowercased official MGI symbol ("trp53 rather than p53, rptor rather than raptor,
kiss1r rather than gpr54") and already resolves 18 of the 48 surveyed genes. The
existing free-text `genetic_background` field becomes the declared RAW half and
the new field the RESOLVED half, the pattern this schema already uses for
`population_raw`/`population_id` and `target_raw`/`target_id`. No existing value
is rewritten.

It enters `CONDITION_FIELDS` and therefore `claim_key` and `condition_key`. It
does **not** enter `RECORD_ID_FIELDS` and does **not** enter
`CONDITION_FEATURE_FIELDS`.

> **Correction: the `age_group` precedent cited here was half wrong, and the
> half that mattered pointed the other way.** A second session checked the
> tuple. `age_group` **is** a `RECORD_ID_FIELDS` component — as are `sex`,
> `diet` and `challenge`, so **all four pre-existing condition dimensions are
> hashed**. The precedent supports excluding genotype from the model's feature
> set, where `age_group` is correctly absent, and is an argument for
> *inclusion* in the hash. The sentence below leaned on both halves and was
> entitled to only one.
>
> **The decision does not change, and here is the argument that actually
> carries it.** Three grounds, none of them a precedent:
>
> 1. A fourth identity map is permanent debt borne by every external reference
>    to this dataset. Three already exist, one of them added two days ago by
>    D-058.
> 2. **The pooled-genotype record is answered by design, not by a split.** The
>    live case is `ad4a6ccfc070e39c`, whose quote reports one joint result for
>    an NPY-intact and an NPY-null arm: "a slight increase in locomotor
>    activity in both hM3DGq-AgRP; NPY-wt/wt, as well as in ... NPY-Δ/Δ-mice".
>    D-060 gives that record the `mixed` value. Splitting it would derive two
>    observations from a sentence reporting one, which this project refuses
>    everywhere else — D-044 makes a quote naming two states a refusal rather
>    than a pick, and D-052 records exactly that refusal for `challenge` on two
>    records of the same study.
> 3. Therefore no curation action the schema permits requires the genotype in
>    the hash. Where a paper reports the two arms separately, the records differ
>    in `outcome_raw` or `source_locator` and no collision arises.
>
> **The collision is real, and ground 3 holds only while something enforces
> it.** Verified: splitting `ad4a6ccfc070e39c` into its two arms yields the
> identical id `ad4a6ccfc070e39c` for both halves, because they share all
> twelve hash components. `syn4rch.curate.run` now refuses a batch in which two
> admitted records share a `record_id`, printing which fields they differ in.
> Its previous check compared admitted ids only against the existing file and
> could not see a batch colliding with itself; such a batch was admitted and
> raised much later in `TaskDataset.record_index`, at the wrong place and time
> for a curator to read. `correct_conditions` already refused on collision
> (`CorrectionCollision`); the curation gate now does too, and a test asserts
> the two-records-one-id property directly so that genotype silently entering
> the hash would fail loudly rather than void this entry.

That the exclusion from the feature set follows `age_group` is sound: `age_group`
is in the condition key and in `condition_specificity` (4 slots) and is not a
model feature (3 slots), and `task.py` already documents that both counts are
printed "rather than one being quietly substituted for the other".
Genotype makes those 5 and 3. Comparability and prediction are different
questions — D-036's finding — and this dimension is on the comparability side by
construction.

Absence is fail-closed per D-044, in five states. `None` (the curator has not
addressed it; all 508 records start here), `not_reported` (the paper was read and
states none), `wild_type` (a *positive, quoted* establishment that the animals
carry nothing beyond the recorded perturbation), `mixed` (a positive value for a
record whose own quote reports two genotypes jointly — Brüning's case, which is a
fact about the record rather than an absence), and one or more segments.
`condition_key` pools `None` with `not_reported` into `unknown`, inheriting D-023
with no code change, and does **not** pool `wild_type` with `unknown`. An
unstated background and an established wild-type background therefore never share
a cell. `strain: C57BL/6J`, populated on 8 records, is **not** sufficient to
write `wild_type`: a strain says what the animals are, not that nothing further
was done to them.

### The two risks were the wrong way round, and both options were measured

| | (a) condition dimension, unhashed | (b) `population_key` suffix |
|---|---|---|
| `record_id`s re-identified | **0** | 0 |
| model primary cells | 305 → **305** | 305 → 314 |
| cell strings whose identity changes | **none** | **all 305** |
| circuit records in a splitting cell | **0** | **77** |
| the 154-row established slice | **untouched** | perturbed |

Option (b) — the cheap-looking one, and the one the task framed as safe because
it re-identifies nothing — is safe for `record_id` and **unsafe for the
pre-registered baseline**, moving 77 of 314 circuit records in exactly the slice
the established +0.208 is measured on. Option (a) threatens `record_id` only if
it is hashed, which this design does not do. The panel scored (b) last on every
lens, 14 of 30.

### Measured consequences of the decision, on the 71 establishable records

| | before | after |
|---|---|---|
| claim cells | 406 | **413** |
| model primary cells | 305 | **305** |
| singleton model cells | 223 | **223** |
| cross-study conflicts | 8 | **4** |
| `QC_CONFLICT_OPPOSED` | **1** | **0** |

**D-052's opposed pair dissolves.** It was exactly
`['0c27e14215860fb0', 'ad4a6ccfc070e39c']` — Herzog decrease against Brüning
increase. They separate because Herzog deleted `Npy` only in AgRP neurons
(`cell_type_restricted`) while Brüning used a whole-animal knockout (`global`),
and those are genuinely different experiments. Three of the seven effect-vs-null
conflicts go with it. Zero model cells split, so the direction result cannot
move — that is a proof from the cell key, not a prediction.

### Two corrections to the process that produced this

**The judges' labels could not be trusted, and the tally they implied was
inverted.** Three judges scored four proposals on separate lenses. All three
described the proposals correctly and two of three **swapped the numbers of
proposals 3 and 4**, so aggregating on the number gave 25/30 to the wrong design.
Re-aggregating on the descriptions gives the arm-scoped design 25/30 with two
judge wins, and the scalar variant 18. Had the numeric tally been trusted, the
wrong design would have been built. Corroboration that the arm-scoped winner is
right came from the verification, independently: **all 12 confirmed background
genotypes returned `applies_to_all_animals: false`** — only some arms — so arm
scoping is what the data demands rather than a refinement.

**The winning design's central engineering claim was false, and checking it took
one command.** It claimed the whole change was "four tuple appends" with
`models/task.py` untouched. But `src/syn4rch/models/task.py:555` reads

```python
sex, diet, _age_group, challenge = condition_key(record)
```

a fixed four-tuple unpack. A fifth condition field raises `ValueError: too many
values to unpack` there and takes the entire model layer with it. The real change
is five edits, the fifth being to read those three features by name rather than
positionally. Every number in the table above was measured with that fix applied;
without it nothing runs at all.


---

## D-062 - The genotype fill: D-052's contradiction dissolves, and four of eight studies had no record in the arm

**The curation pass D-060 designed.** Eight readers, one per study carrying an
adversarially confirmed background genotype, each deciding every record of its
study under the arm rule. 78 assignments, **zero positive values arriving
without a verbatim quote**, 69 of 78 at high confidence.

| | |
|---|---|
| filled | **77** |
| left `null` (arm not establishable) | 1 |
| refused | 0 |
| evidence entries added | 51 |
| record ids unchanged | **508 of 508** |
| identity map required | no |

### D-052 is resolved

`QC_CONFLICT_OPPOSED` goes **1 → 0**. The pair was
`['0c27e14215860fb0', 'ad4a6ccfc070e39c']` and it separates like this:

| record | direction | `background_genotype` |
|---|---|---|
| Herzog `0c27e14215860fb0` | decrease | `npy:loss_of_function:cell_type_restricted` |
| Brüning `ad4a6ccfc070e39c` | increase | **`mixed`** |

Herzog deleted `Npy` only in AgRP neurons. Brüning's record reports a *joint*
result across an NPY-intact and an NPY-null arm — "a slight increase in
locomotor activity in both hM3DGq-AgRP; NPY-wt/wt, as well as in ...
NPY-Δ/Δ-mice" — so it is comparable with neither of them. Two different claim
cells, and the contradiction was never a contradiction.

That `mixed` is doing the work here matters. D-060 invented it for exactly this
record, and the alternative — splitting it into one record per genotype — would
have derived two observations from a sentence reporting one, which D-044 refuses
in general and which the `record_id` collision demonstrated concretely (both
halves hash identically).

### What else moved

| | before | after |
|---|---|---|
| claim cells | 406 | **416** |
| cross-study conflicts | 8 | **4** |
| `QC_CONFLICT_OPPOSED` | 1 | **0** |
| `QC_LOCATOR_COLLISION` | 53 | **52** |
| records stating no condition at all | 139 | **117** |
| total QC findings | 2723 | 2718 |

The locator collision that resolved is a *correct* resolution and was not
predicted. Collisions group on (study, locator, claim), so two records at one
locator that differ only in genotype had been reported as one observation
extracted twice. They are two experiments, and the finer key says so.

**The established result is untouched**, as D-060 proved from the cell key
rather than predicted: +0.208 [+0.089, +0.354], p=0.004, 154 rows / 35 studies,
`established` — and identically under the corrected `lab` grouping. Zero model
cells split, because `background_genotype` is not in
`PERTURBATION_AXIS_CELL_BACKOFF`.

**And the cost, which is real.** One `(population, outcome)` pair dropped from
`marginal` to `insufficient`, because the dimension split its best claim in two
and neither half clears a gate counting studies. The one `feasible` task
survives. This is the same trade D-023 made for the other four conditions: a
finer key buys comparability and spends replication.

### Four of eight studies had no record in the carrying arm, which is the finding

Four readers returned `is_really_background: false`, and **not one of them
because the segment was bogus or was really the perturbation.** In every case
the confirmed background genotype is genuine for that paper and no *curated
record* comes from the arm that carries it:

* `advs.202403245` — the A53T-SNCA background belongs to the alpha-synuclein
  arm; the single curated record is from the MPTP cohort.
* `cmet.2023.03.019` — the NPY-null arms are Figures 2 and 6B; all 8 curated
  records are NPY-intact.
* `s41467-017-02534-9` — the four Pdx1-driven segments sit in arms neither
  curated record is from.
* `nn.2739` — 2 of 7 records are in the A^y arm, 5 are not.

So the genotype dimension is sparser than the survey implied for a reason that
is about *which results were curated*, not about the papers. 31 records carry a
segment, 15 an established `wild_type`, 5 `mixed`, 26 `not_reported`, and 431
remain `null`. A reader must not take the four `false` flags as verdicts against
the segments; the segments are confirmed, and the arms they belong to are
mostly uncurated.

### A refusal that was right, and a promotion I declined

All 26 `not_reported` rows arrived **carrying a quote**, and the writer refused
every one of them. It was right to. The quotes are *arm* evidence — "MPTP-induced
weight loss", "we crossed our conditional Npy knockout line", "normoglycaemic
POMC-cre mice" — and none of them asserts that nothing further was done to the
animals. A009 and A011 forbid an evidence entry on a value that asserts an
absence, and an absence cannot be quoted. The quotes are retained in the
assignments table as `arm_evidence_quote` rather than dropped, because they are
why the reader chose the arm and a reviewer needs them; they simply do not enter
`condition_evidence`.

**And 13 of those 26 are arguably `wild_type`, which I deliberately did not
promote.** They come from `s41598-019-49295-7` and quote "normoglycaemic
AgRP-ires-cre" or "normoglycaemic POMC-cre" mice in a paper whose other arm is
`ob/ob` — naming the complete allele complement, which is the exact criterion
another reader used to assign `wild_type` in `cmet.2023.03.019`. Promoting them
would be the better-looking result: `wild_type` is a positive value that does
not pool with `unknown`, so it buys comparability that `not_reported` does not.
It would also be me overriding a curation judgement on my own reading of a
truncated quote, which is the thing this project's whole structure exists to
prevent. They stay at `not_reported`, flagged in the table for a later pass that
can read the papers properly.

### What the table is

`data/curated/background_genotype_assignments.json`, tracked and reviewable in a
diff, one row per record with its value, quote, locator, source, confidence and
the reader's reasoning. `syn4rch.curate.fill_background_genotype` applies it and
decides nothing: it fills only where the value is `null` (D-039 — overwriting is
`correct_conditions`' job), validates every segment against
`data/vocab/targets.json`, attaches a `ConditionEvidence` entry to every
positive value and refuses one on an absence, re-checks every quote against the
passage its locator names, and asserts per record that no `record_id` moved.


---

## D-061 - Europe PMC returned MeSH and affiliations all along; the parser discarded them

**The defect.** `syn4rch.harvest.europepmc.search` has always requested
`resultType=core` — the Europe PMC result type that carries `meshHeadingList`
and author affiliations. `_norm_record` hardcoded `"mesh_terms": []` and
`"last_author_affiliation": None`, the second under a comment asserting *"Europe
PMC search does not return affiliations"*. Meanwhile `screen.py:197` joins MeSH
into the text it matches against and `screen.py:302` stamps every candidate
`screened_on: "title+abstract+mesh"`.

So the screener was told it had read MeSH, and had not. Verified cross-tab over
all 9,653 candidates:

| source | MeSH present | MeSH empty |
|---|---|---|
| `pubmed` | 5,020 | 238 |
| `europepmc_oa` | 325 | **2,519** |
| `europepmc_abstract` | 321 | **1,230** |

Same family as D-056: a check that cannot fire, whose silence reads as a result.

**Both fields were already on disk, under names no reader guesses.** MeSH at
`meshHeadingList.meshHeading[].descriptorName`; the affiliation per author at
`authorList.author[].authorAffiliationDetailsList`, with a result-level
`affiliation` as fallback. My own first probe missed the affiliation because it
looked for a key called `affiliation` inside each author. Re-parsing the 137
cached search responses — no network — takes MeSH coverage from 0% to **85%**
and last-author affiliation from 0% to **91%**.

Qualifiers are dropped on purpose: the screener matches descriptor names, and
`Hypothalamus/metabolism` would match neither `Hypothalamus` nor `metabolism`.

### The measurement, and the wrong version of it I ran first

**Published because the first number was wrong and I quoted it before checking.**

My first experiment rebuilt each candidate from the Europe PMC cache and
compared its screen against the *stored* candidate. That measures the difference
between two **sources**, not the difference the fix makes: a stored row merged
from the PubMed path carries pubtypes and abstract text the Europe PMC row does
not. It reported 64 changed decisions and 31 newly includable papers.

`PMC7617324` is what exposed it — `exclude` → `include` with an empty
gained-terms set, which is impossible if a new match caused it. Its stored row
carries a **preprint** pubtype worth −8, so it scored 5 and stored as `exclude`;
the Europe PMC row has no pubtype at all, scores 11, and screens as `include`.
Nothing to do with MeSH. I had said at the time that "an unexplained transition
in a 64-row result is the row that turns out to be the bug", and it was.

**The correct experiment holds the record fixed and varies only the one field the
fix populates** — screen each cached record twice, once with `mesh_terms`
emptied exactly as the old parser produced it:

| | |
|---|---|
| decisions changed by MeSH alone | **98** |
| `uncertain` → `include` | 65 |
| `exclude` → `uncertain` | 33 |
| newly includable papers | **65** (25 open access, 36 with a PMCID) |
| transitions with no gained term | **0** |

Every change is in the permissive direction, which is what a missing matcher
input predicts. The corrected figure is *larger* than the wrong one, so the
error was not self-serving — but it was still wrong, and the entry carries both
because a reader meeting 64/31 in this conversation's history needs to know
which number survived.

What MeSH contributes is mostly the **anatomical site**, which a title often
omits: `hy` 22, `optogenetic_unspecified` 11, `gonadotropins` 9,
`corticosterone` 9, `sleep` 8, `blood_glucose` 7, `food_intake` 6, `pvh` 5.

### Two magnitudes, and which one means what

"~3,750 candidates carried a provenance string asserting evidence never
fetched" and "98 decisions and 65 papers changed" are both true and differ by a
factor of forty. The first is the size of the **provenance defect** — rows whose
`screened_on` lied. The second is the size of the **consequence**. Letting the
first imply the second would be the mirror of D-058, where a replication count
fell by one while the replication evidence improved; those read as opposites
until you say which number measures what.

### The boundary: this does not reach the curated dataset

**Found by the session curating the genotype dimension, and verified here before
being accepted.** Curated records do not get their affiliation from the search
API at all. `curate/jats.py:748-767` sets `last_author_affiliation` from the
full text, resolving the last author's `<xref ref-type="aff">` and falling back
to the sole affiliation when the front matter lists exactly one — refusing rather
than guessing when several are unlinked. Measured on the 508 records:

| harvest path | affiliation present | missing |
|---|---|---|
| `europepmc_oa` | 350 | 17 |
| `pubmed` | 141 | 0 |

So the 17 gaps are papers whose JATS listed several unlinked affiliations and
correctly got none, not papers the API failed. **D-056's 68 lab groups do not
move**, and I withdraw the claim I made that the homonym check would gain
affiliations on this path — it already had them.

This fix therefore matters at the **candidate and screening stage** and nowhere
else. Said plainly here so the entry cannot be read as a claim about the
dataset.

### What this does not do

* The committed `data/interim/candidates.jsonl` still carries
  `screened_on: "title+abstract+mesh"` on rows screened without it. The parser
  fix makes future screens honest; those rows stay wrong until the harvest is
  replayed. That file is untracked (D-016), so replaying is cheap, and it is a
  separate decision from this one.
* The 65 newly includable papers are **not curated**. 36 have a PMCID so their
  full text is retrievable; none is in the corpus. Whether they are worth
  curating is a question for the next selection pass, not an implication of this
  entry.
* Nothing here touches the schema, the keys, the vocabularies or the dataset.


---

## D-063 - What a Tier A readout is: a state of the organism, not the activity of neurons

**Prompted by a reader asking the right question.** Told that 91 uncurated EEG
papers in the retrieved corpus were the largest available seam, they asked
whether EEG scans help *"as they are on the cortex"*. The project had never
written down an answer. The rule existed as a single clause of D-054 — "Tier A
is systemic physiology" — invoked to decline anxiety, impulsivity and operant
lever-pressing, and nowhere else.

**The boundary matters because this is a hypothalamic dataset measured in the
cortex.** Stimulate AgRP neurons, record EEG, and the electrode is reporting
cortical activity, not what the animal's body did. Whether that is an outcome at
all decides how much of the EEG seam is admissible, and the answer turned out to
be about half of it.

### The rule I first proposed was wrong

I had been saying "Tier A readouts measure the body, so cortical measures are
out". Reading the vocabulary properly kills that. `inflammatory_markers_other`
deliberately admits `microglial activation`, `Iba1`, `GFAP`, `gliosis` and
`neuroinflammation`, and its notes say why: the hypothalamic-inflammation
literature reports them *as* the systemic-inflammation readout. A
central-versus-peripheral rule would have retroactively invalidated an entry
that was argued for, which is the kind of tidy rule that is worse than no rule.

### The rule, stated against the schema rather than against taste

**A Tier A outcome is a state of the organism, or a tissue measurement the
literature uses as a declared proxy for one. It is not the activity of neurons.**

The basis is the unit of record. SCHEMA.md 1 makes it
`(population) x (perturbation) x (physiological readout)`. Neural activity sits
on the **perturbation** side of that product: it is how a manipulation is
verified to have worked, not what the manipulation did to the animal. Curating
it as an outcome would put the same quantity on both sides of one record.

That line sorts every case cleanly:

| readout | verdict | why |
|---|---|---|
| `sleep` (min/h of NREM, REM, wake) | **in** | a whole-animal state; EEG is the instrument, no more the readout than a calorimeter is |
| glial markers under `inflammatory_markers_other` | **in** | a tissue state the field reports as the systemic one, with the proxy stated |
| EEG delta/theta power, spindles, oscillations | **out** | cortical activity |
| c-Fos, firing rate, LFP, calcium imaging, patch clamp | **out** | neural activity, and usually the manipulation's own verification |

**The instrument is not the test.** This is the part the question exposed: all 36
Tier A outcomes measure the body, and two of them — `sympathetic_nerve_activity`
and `parasympathetic_tone` — *are* electrophysiology. Electrodes were never
disqualifying. Peripheral efferent outflow is the body's autonomic command
signal; cortical voltage is another part of the brain.

### What changed, which is one deletion

`EEG delta power` comes out of `sleep.synonyms`. It invited a cortical spectral
measure to be curated against an id whose units are `min/h` and `% of time` —
putting a spectral quantity in the same claim cell as a sleep duration, which is
worse than the NREM/REM pooling D-064 addresses, because those at least measure
the same kind of thing.

`data/vocab/outcomes.json` goes to 0.4.0 with two new `_meta` blocks:
`scope_policy`, the rule above, and `out_of_scope_readouts`, 39 declared terms in
two categories — `neural_activity`, excluded by this entry, and
`behavioural_affective`, which D-054 had already declined and which was until now
only remembered.

`tests/test_vocab_scope.py` asserts no entry claims any of those 39 as a synonym
of **any** kind, because a readout smuggled in as an `inverted_synonym` would be
curated just as readily as a plain one. It also pins the glial exception *and its
justifying note*, so a later scope sweep cannot quietly delete the documented
case and leave a tidier rule behind.

### Preventive, not corrective

**No curated record violates this.** Zero of 508 records have an `outcome_raw`
that is an out-of-scope readout, and zero even contain `delta power`,
`spectral`, `firing rate`, `c-Fos`, `LFP` or `calcium`. Nothing was rewritten and
no direction changed — the same shape as D-049, where thirteen polarity traps
were annotated and all eight records at risk turned out already correct. Saying
so matters: a scope pass that silently dropped records would be
indistinguishable in the diff from one that corrupted them.

### What it costs, measured

Of the 638 proposed extraction units in the EEG seam that map to `sleep`:

| | units | |
|---|---|---|
| state amounts only | **334** (52%) | admissible |
| cortical/spectral only | 31 (5%) | excluded by this entry |
| both in one passage | 151 (24%) | a curator must separate the sentence |
| neither | 122 (19%) | methods and introduction passages |

So the rule removes 5% of that seam outright and complicates another quarter.
The seam remains the best thing available and is about half the size I first
told the reader it was, which is a direct consequence of their question.

### What this entry does not do

* It does not extend the vocabulary. The `sleep` id is still one id; splitting it
  into state amounts is D-064, which this entry is the premise for.
* It does not settle the eight readouts D-054 declined for want of an id. Those
  are a different question — no home rather than out of scope — and D-054's rule
  stands: adding ids is a decision entry, never a curator's judgement mid-pass.
* `out_of_scope_readouts` is not a licence. A readout absent from the list is not
  thereby in scope, and the policy says so in those words.


---

## D-064 - Splitting `sleep` into the five quantities it was pooling, and the row the primary metric correctly lost

**The premise is D-063**, which fixed what a Tier A readout is. This entry acts
on it. One `sleep` id held `NREM sleep`, `REM sleep`, `slow wave sleep`,
`wakefulness`, `sleep latency` and `sleep fragmentation` as synonyms, with units
`min/h`. NREM and REM are not complements and move independently under one
manipulation, so a paper's NREM time and its wakefulness time landed in the same
claim cell — the defect class of D-050 (projections) and D-052 (genotypes), a key
too coarse to tell two experiments apart.

`data/vocab/outcomes.json` goes to 0.5.0 with five new Tier A ids —
`nrem_sleep`, `rem_sleep`, `wakefulness`, `sleep_latency`,
`sleep_fragmentation` — and `sleep` re-scoped to the composite readout that does
not resolve to one state.

**The moves were forced, not stylistic.** `syn4rch.vocab.core.selfcheck` refuses
a normalised string mapping to two ids in one file, so each state term had to
leave `sleep` to live in its own entry. `wakefulness` and `sleep latency` stop
being `inverted_synonyms` and become quantities in their own right, which takes
a hand-flip away from the curator: `direction` on `wakefulness` now describes
wakefulness, and more wake is `increase`.

**Bare `arousal` is deliberately a synonym of nothing.** It names the waking
state in one literature and a count of awakenings in the other. `arousal index`
and `arousals` resolve to `sleep_fragmentation`; bare `arousal` resolves to
nothing and forces a reading, which is what D-021 says this vocabulary does with
a trap it cannot settle.

### The curators had already found this defect, and declined data over it

The strongest evidence for the split is not a measurement. It is three
`extraction_notes` written during the circuit pass, each recording observations
the reader **refused to curate** because the vocabulary could not hold them:

> "The same paragraph reports a significant REM-sleep *decrease* and a longer REM
> latency, which map to this same outcome_id in the opposite direction; they are
> deliberately not proposed, because an increase and a decrease in one cell would
> be a contradiction manufactured by the vocabulary having no NREM/REM
> distinction."

Five such observations across three records are now recordable: a REM decrease
(`Results para 8`), a non-significant REM increase and a shorter NREM latency
(`para 10`), and a REM decrease and a longer REM latency (`para 11`), all in
`doi:10.1371/journal.pone.0020360`. **They are listed, not curated** — adding the
slot and filling it are separate named passes, following D-060 and D-062.

**And an independent confirmation from a different paper.** A trial curating five
deposit-bearing EEG papers found `doi:10.7554/eLife.69909` reporting **seven**
sleep readouts in one paragraph for one population, one perturbation and one
condition: wakefulness +140%, NREM −81.6%, REM −94.5%, NREM in the following dark
period unchanged, low delta power up, high delta power down, NREM power density
unchanged. Under the old id that single cell takes `increase`, `decrease` and
`no_change` simultaneously, and the conflicts detector cannot tell it from a real
contradiction. Two of those seven are delta power, which D-063 had already put
out of scope.

### What moved, measured before and after

D-023's standard: a finer key **buys comparability and spends replication**, and
both halves are the result.

| | before | after |
|---|---|---|
| sleep-family records | 10 | 10 (7 `sleep`, 3 `nrem_sleep`) |
| claim cells they occupy | 7 | **8** |
| cells with >=2 studies | 1 | 1 |
| cells with >=2 labs | 1 | 1 |
| conflicts involving them | 1 | 1 |
| corpus claim cells | 416 | **417** |
| corpus conflicts | 4 | 4 |

**Only 3 of 10 records moved**, and the refusal rule is why: a record is moved
only when its `outcome_raw` names one state *and* its stored direction is already
correct under the new id. Flipping a direction is curation, not a vocabulary
change, and needs the quote re-read. **Under-moving and saying so beats flipping
a stored direction inside a vocabulary pass.**

**The rule as I first wrote it was keyed on the wrong thing**, and the session
curating the genotype dimension corrected it. I had said "refuse to move any
record whose `outcome_raw` matches an inverted synonym". The right statement is
"refuse to move any record whose stored direction was **derived by inversion**",
and `c9d907bbb8f604ea` is why: its `outcome_raw` is "wakefulness and REM sleep
time in the dark period", naming one inverted term and one plain one, and its
`increase` is a synthesis of two raw measurements moving in opposite directions
that happen to agree about sleep. A substring rule catches that record for the
wrong reason; a rule keyed on the derivation catches it for the right one, and
the record is irreducibly a `sleep` record that can move to no state id at all.

Seven stayed: four because the quote names two states jointly or a whole-cycle
composite (`sleep-wake amounts`), three because the direction was derived by
inversion from a wakefulness or arousal measurement.

**The three ambiguous ones were checked against their quotes before being frozen
in place**, because a refusal rule applied to an already-wrong record preserves
the error rather than preventing it. `19e2a8a1115e2cc7` — "the wakefulness time
was significantly **longer**", stored `sleep: decrease`. `c9d907bbb8f604ea` —
"wakefulness **decreased** accompanied by REM sleep **increasing**", stored
`sleep: increase`, both halves agreeing. `bd0fa9320548db4c` — arousals terminate
sleep, stored `sleep: decrease` per D-019. All three correctly inverted, with the
reasoning already in `extraction_notes`. So D-049's finding — that every record
whose `outcome_raw` matches an inverted synonym of its own `outcome_id` is
correctly stored — now holds at 508 records, having been established at 362.

That 7-of-10 figure is itself a finding: the existing sleep records are mostly
composites *because* there was one id to aim at. A reader given `nrem_sleep` and
`rem_sleep` records two observations where a reader given `sleep` records one.
So the split's value is almost entirely prospective — the 334 admissible
state-amount units in the EEG seam — and almost nothing retrospectively.

`outcome_id` is not a `RECORD_ID_FIELDS` component, so **no record was
re-identified**: 508 of 508 ids unchanged, and across the whole dataset the only
fields that changed anywhere are `outcome_id` and `extraction_notes`, on those 3
records.

### The primary metric lost one row, and it was right to

`reports/generalisation.md` moves from **+0.208 [+0.089, +0.354], 154 rows / 35
studies, p=0.004** to **+0.203 [+0.086, +0.347], 153 rows / 35 studies,
p=0.005**. Still `established` under both groupings, and still clear of D-059's
secondary threshold of +0.15. But the D-055 and D-057 baseline is now +0.203, and
that is this entry's doing, not the genotype pass's.

The row is `75c57af3d1474438`, identified by diffing the cell-seen set against
the pre-split dataset in git — exactly one left, none entered. Its old
`lha_orexin` / circuit-activation / `sleep` cell had three members:

| | direction | `outcome_raw` | lab |
|---|---|---|---|
| other study | decrease | NREM and REM sleep time (wakefulness time rose) | Niu J. |
| own study | decrease | NREM sleep time | Sakurai |
| own study | decrease | wakefulness time | Sakurai |

So its only cross-study evidence was an **NREM+REM composite** standing in for an
NREM-only measurement. After the split the composite stays on `sleep` and the
NREM record moves to `nrem_sleep`, so the link is severed and the record's cell
holds only its own study — which leave-one-study-out empties.

**The row was comparable only because the id was too coarse.** The metric got
smaller and more honest in the same move, and that is the shape of every key
refinement in this project: D-050 removed manufactured replication, D-058 removed
a false match and a false conflict in one pair, and this removes a cross-study
link that was comparing two different quantities.

### And the feasibility verdict paid the same price

`1 of 181 tasks feasible (0 marginal, 180 insufficient)`, against
`1 of 180 feasible (1 marginal, 178 insufficient)` before. The marginal task is
gone and the task count rose by one, because the split created a new
(population, outcome) pair and took records from an existing one. The genotype
pass cost a marginal task the same way, which the session running it warned me
about before I measured — so this is the second instance of a predicted cost
arriving on schedule, and D-023's framing has now survived two tests.

### What this entry does not do

* It does not curate the five declined observations, nor the 334 EEG units. The
  ids exist and are nearly empty by design, the way `background_genotype` was
  added by D-060 and filled by D-062.
* It does not move the seven records that stayed. Three of them need a direction
  flipped under a new id, which is a correction pass needing the quotes re-read.
* `atlas/signaling.py` needed five `OUTCOME_SYSTEM` entries to accept the new
  ids, all mapping to `behaviour` as `sleep` already does. That file belongs to
  another session this pass; it made the change rather than me.


## D-065 - Mapping populations to atlas clusters: the coverage invariant that rejected everything, and why ambient RNA made it do that

**A second session proposed three structural fixes** after an adversarial pass
refuted three of thirteen hand-made population → HypoMap cluster mappings. All
three were aimed at real gaps. One of them, applied literally, rejected 12 of 12
mappings — and that failure is the most useful thing in this entry.

`data/curated/population_cluster_mapping.json` now holds the mapping for all 30
declared populations, built by `syn4rch.atlas.population_clusters` from the
atlas (`doi:10.1038/s42255-022-00657-y`, 384,925 mouse cells × 51,216 genes,
asserted on `NCBITaxon:10090`). **Tier B, and never a training target.** A
mapping here says a cluster plausibly corresponds to a population the
vocabulary already declares; it is not evidence that any intervention did
anything.

### The proposed coverage invariant failed, and the reason generalises

The proposal: for each mapped population, record what fraction of the atlas's
marker-positive cells the chosen cluster set captures. It targets the one defect
adversarial verification *cannot* catch — `scn_vip` was confirmed while holding
roughly 44% of the high-Vip SCN cells, because the mapping was not wrong, only
incomplete.

Measured, it rejected everything. 12 of 12 populations held under half their
marker-positive cells; `lha_mch` held 1%. The denominator is the problem: 58,480
of 384,925 cells read as Pmch-positive. **15% of the hypothalamus is not MCH
neurons.** Droplet dissociation leaves a soup of transcripts from the most
highly-expressing cells, and a highly-expressed neuropeptide is detectable in
nearly every droplet. An absolute threshold on log-normalised counts cannot tell
ambient contamination from expression.

Transcript *mass* share fails for the same reason, which is less obvious. Mass
looked robust — ambient contributes little per cell — but there are vastly more
contaminated cells, so the product dominates. Ranked by Pmch mass, the top three
clusters in the atlas are astrocytes, oligodendrocytes and astrocytes again;
`Pmch.GLU-7`, with a mean of 6.316 against their 0.392, comes fourth.

**A metric that fails every case is measuring itself, not the cases.** What
survives is a scale-free ratio: `contrast`, the median marker mean of the chosen
clusters over the median of the region candidates left out. Ambient lifts every
cluster together, so a ratio is unchanged by it.

### Three refusals that were mine, not the data's

* The region table invented `Tuberomammillary`. The atlas does not use it — the
  histaminergic cluster is filed under `(Pre)Mammillary region` — so
  `tmn_histamine` refused for want of a label, and now maps at 118× contrast.
* The rule was pinned to one cluster level. `pvh_crh` has no PVH-majority Crh
  cluster at C185 at all; it resolves only at C286/C465. The level is now chosen
  per population by greatest contrast — not the finest, which averages fewer
  cells, always wins on raw magnitude, and would fragment every population.
* Within-cluster rank was treated as specificity. Crh ranks **544th of 51,216
  genes** inside the cluster HypoMap's own authors named `Crh.Tent5a.GLU-2`.
  Rank tracks abundance, which differs by orders of magnitude between a
  neuropeptide and a transcription factor, so gating on it silently selects for
  abundance. Rank is recorded and never gated on.

A fourth was caught by the tests rather than the data: contrast was first
measured best-against-all-others, which counts a population's own siblings
against it and penalises it for genuinely spanning several clusters —
`arc_agrp` occupies five C465 clusters, `scn_avp` seven.

### What the atlas authors' own names decide

The tier is the output a caller should read, and tier A requires the marker to
appear in HypoMap's name for *every* chosen cluster: an expert judgement this
rule did not generate. Name support and contrast are combined rather than
gated, because an earlier version discarded full name support on a round
number — `poa_galanin` has Gal in all four of its cluster names and contrast
9.01, and was filed beside `zi_gaba`, whose pan-GABAergic marker keeps the whole
zona incerta and picks nothing out of it.

**11 tier A, 4 tier B, 8 tier C, 7 refused.** The refusals form a regularity
worth stating: every population defined by a neuropeptide maps at 17–140×
contrast, and every one defined by a receptor or a kinase fails the floor —
`Lepr` peaks at 0.074 across the DMH, `Mc4r` at 0.055 across the PVH, `Ikbkb`
at 0.328. `Nr5a1` reaches 0.391 and `vmh_sf1` is refused. These are real
populations with real Cre lines; the driver gene simply does not carry enough
transcript to mark a cluster. **A Cre line is not a marker**, and the mapping
says so rather than guessing.

### The rule was declared after looking, and the file says so

`declared_after_inspecting_the_atlas: true` is recorded in the artefact. This is
not pre-registration in the sense of D-048/D-056 and must not be cited as if it
were: the parameters were chosen knowing what the atlas contains. What the rule
buys is reproducibility and uniform application to all 30 populations, including
the ones where it returns an answer nobody wanted.

### One nesting survives, and D-059 has to be able to see it

`arc_kisspeptin` is a subset of `arc_tbx3`. That is the biology — Tbx3 is
expressed in the ARC Kiss1/Tac2 neurons — and it is reported rather than
forbidden. But two populations sharing clusters have expression vectors
correlated **by construction**, so a similarity ranking would find them alike
for a reason that is an artefact of this mapping rather than a fact about
transcriptomes. The pre-registered per-row diagnostics of D-059 must therefore
carry the cluster-set overlap of each nearest-neighbour pair, or the experiment
can rediscover its own input as a result.

The third proposal — assert region composition sums to ~100% — passes on all
1002 clusters. The composition data was always sound; the truncation that
manufactured one of the three refutations was in the consumer, which kept only
the top region per cluster and discarded a 44% PVH share.

## D-066 - Third amendment to D-059, written while it provably cannot change the answer: the free parameters close, and overlapping populations are disqualified rather than diagnosed

**D-059 says the gene set, `k` and the mapping must all be fixed in code before
the first evaluation.** D-065 fixed the mapping. This entry fixes the other two,
plus two questions D-065 raised and could not settle alone. Nothing here has
been measured against an outcome; the primary and secondary of D-059 remain
unmeasured.

### Why this amendment is legitimate, stated as a checkable fact

D-065 found one surviving nesting: `arc_kisspeptin` is a subset of `arc_tbx3`.
Two populations sharing clusters have expression vectors correlated **by
construction**, so a similarity ranking could rediscover this mapping's own
input as a finding.

**`arc_tbx3` holds zero records.** It can neither lend a direction as a training
neighbour nor appear as a test row, so the overlapping pair cannot arise in this
evaluation at all. That is what makes now the only moment this rule is free: a
decision taken while it provably cannot affect the result costs nothing, and
after the run it is unavailable at any price. The second session raised this and
it is the right reading — changing a pre-registration before running is still a
change, and doing it for a pair that cannot occur buys nothing while spending
the protocol's integrity.

### Overlapping populations are disqualified, not merely logged

**The rule: if a test row's nearest training neighbour shares one or more
clusters with the held-out population, that neighbour is skipped and the next
nearest is used. If every available neighbour overlaps, the row falls back
exactly as the existing baselines do and is counted as unpredicted, not wrong.**

The argument is not biological — the nesting is real and the similarity genuine.
It is about which choice leaves the primary interpretable **whichever way it
lands**:

* Admit overlapping pairs and the primary passes, and you must then ask whether
  it passed only on them. If it did, you need a re-run excluding them — a second
  look at a pre-registered metric, which D-048 forbids. That branch has no
  honest exit.
* Disqualify them and a pass is a pass on the clean set, which is the claim
  worth making. A failure is a failure on the clean set, and whatever the
  overlapping pairs would have contributed becomes a labelled secondary.

Both branches of the second are interpretable; one branch of the first is not.

**The diagnostic is additive and logged regardless.** For every test row the run
records the cluster-set Jaccard overlap between the held-out population and the
neighbour that supplied its prediction, alongside the two diagnostics D-059's
second amendment already requires. On this corpus the figure will be 0.0 on
every row, which is the point: it is unrecoverable afterwards, and it goes live
the moment any curation gives `arc_tbx3` a record.

### The gene set, fixed

The **2000 genes with the highest variance of per-cluster mean across all 465
C465 clusters**, excluding mitochondrial (`mt-`) and ribosomal (`Rps`, `Rpl`)
genes as technical covariates rather than identity.

Chosen this way because it depends only on the atlas — not on which populations
carry records, not on the mapping, and not on any delta. Cosine similarity over
all 51,216 genes would be dominated by what is high in every cluster alike:
`Malat1`, `Meg3` and `Snhg11` are the top genes of nearly every cluster in this
atlas, so an unrestricted cosine measures library composition rather than cell
identity.

A **population's vector** is the cell-count-weighted mean of its assigned
clusters' means — the pseudobulk mean over exactly the cells the mapping
assigns. This is level-independent on purpose, because D-065 assigns different
populations at different levels (C66 through C465) and the vectors must be
comparable.

### `k = 1`, and no other `k` can touch the primary

With at most eight neighbours available, several contributing one to three
records, a larger `k` would average over populations that have no record in the
relevant cell in any case. `k = 1` also makes the per-row nearest-neighbour log
*identical to the cause of the prediction*, which is what makes a miss
diagnosable rather than merely disappointing.

`k = 3` with similarity weighting will be computed and reported as an
exploratory figure, **labelled as such, and it cannot change the primary
verdict.** Recorded here so that reporting the better of the two afterwards is
foreclosed.

### The primary runs on tier A only, and that choice costs something

Mapping tier (D-065) gates admission: **the primary admits tier A mappings
only**, because D-065 states in terms that tier C mappings must not be treated
as established, and `lha_vglut2`'s pan-glutamatergic `Slc17a6` makes its vector
a blur over most excitatory neurons rather than a population.

| | circuit records | populations | neighbours available |
|---|---|---|---|
| resolved `population_id` | 245 | 13 | — |
| **primary (tier A)** | **213** | **9** | **at most 8** |
| secondary (tier A+B+C) | 245 | 12 | at most 11 |

Excluded from the primary: `vmh_sf1` 22 records (refused — `Nr5a1` peaks at
0.391 and never clears the marker floor), `pvh_sim1` 4 (tier C), `lha_vglut2` 3
(tier C), `scn_vip` 3 (tier B, and the known ~44% coverage gap).

**The cost, stated rather than buried: purity makes the pool *more* dominated by
one pair, not less.** `arc_agrp` (123) plus `arc_pomc` (52) are 175 of 245
resolved circuit records, 71% — but 175 of 213 in the tier A primary, **82%**.
Six of the nine admitted populations contribute three records or fewer. So the
primary is, in large part, a test of a single pair, exactly as D-059's second
amendment already warned; restricting to tier A sharpens that rather than
relieving it. The tier A+B+C secondary is reported **regardless of the primary's
outcome**, so the information in B and C is neither lost nor available to rescue
a failed primary.

### What each outcome means, pre-committed so none of them is unfalsifiable

D-059 records the prediction that the primary will miss, and the mechanism:
ARC AgRP and ARC POMC neurons are transcriptionally adjacent and drive food
intake in opposite directions. The per-row log separates two failures that would
otherwise look alike:

* **Uninformative.** Delta ≈ 0 and the nearest-neighbour log shows neighbours
  effectively unrelated to the held-out population, or most rows falling back
  for want of a neighbour with the cell. The finding is that the pool was too
  small to tell, not that similarity fails.
* **Actively misleading.** Delta ≤ 0 with the misses concentrated on the
  `arc_agrp`/`arc_pomc` rows, each naming the other as its nearest neighbour at
  high cosine similarity. The finding is specific and mechanistic: expression
  similarity transfers identity but not the sign of a systemic response.

The second is the more valuable result, and it is the one that serves the
standing requirement that this system be able to answer "insufficient evidence
to rank these interventions reliably" — a feature that is confidently wrong in a
nameable way is worse than no feature, and this is how we would know.

**If the primary passes**, the per-row breakdown for `arc_agrp` and `arc_pomc`
is reported in the same table, not separately. Should it pass while those rows
miss, the pass rests on populations contributing one to three records each and
must be reported in those words.

### The rule is pre-registered, not the list

The figures above are computed against the mapping committed at `91b3a28`, and
an adversarial pass over the eleven tier A claims was in flight when this was
written. If it refutes one, that population leaves the primary pool **by the
rule already stated here**, not by a new decision — which is the reason this
entry fixes "tier A is admitted" rather than naming nine populations. The pool
table is therefore a measurement as of that commit, and the run must print its
own pool alongside the result so the two can be compared.

A refutation moving a population out is not an amendment and needs none. A
change to *which tiers are admitted* would be, and is foreclosed.

### Unchanged

Predictions 1 and 2 of D-059 stand as written and are still unmeasured. The
evaluation runs **once**. This remains a Tier B feature over an external
observational atlas, generalising a direction label that comes only from
interventional experiments; D-065's mapping is `never_a_training_target: true`
and nothing here alters that.

## D-067 - The pre-registered similarity metric cannot rank: a 0.0051 margin, the tie guard that belongs here, and the last amendment before the run

**Measured before any delta, and it is a property of the feature, not a result.**
D-059 specifies cosine similarity between population mean-expression vectors.
Built as specified over the gene set D-066 fixed, on the 9-population tier A
primary pool, it is very nearly degenerate:

| metric | min | median | max | spread | median top-1 margin |
|---|---|---|---|---|---|
| **cosine, raw (as pre-registered)** | 0.893 | 0.942 | 0.974 | 0.081 | **0.0051** |
| Pearson across genes | 0.749 | 0.872 | 0.943 | 0.194 | 0.0097 |
| cosine after centring each gene | −0.367 | −0.138 | 0.231 | 0.598 | **0.0605** |
| cosine on z-scored genes | −0.375 | −0.127 | 0.173 | 0.549 | 0.0533 |

Every pair of populations sits between 0.893 and 0.974, and **the gap between a
row's first and second nearest neighbour is about half a percent.** With `k = 1`
the whole prediction turns on that gap. This is the defect a second session
first raised as a `tied` outcome, which I answered by saying `vmh_sf1`'s refusal
had removed the only ranking it applied to. That was wrong — it applies here,
to the similarity ranking, which is where a margin inside its own uncertainty
actually decides something.

The cause is arithmetic rather than biological: cosine between non-negative
vectors that share a large component is always near 1. All these populations are
neurons and share most of their transcriptome, so the shared part dominates the
inner product. Centring each gene across populations removes it.

### The gene set was suspected first, tested, and exonerated

A second session argued the gene set is selected indirectly by well-represented
regions, citing the claim that HypoMap underrepresents the VMH. **Measured, that
premise is false — and the premise was mine, not theirs.** They took the phrase
verbatim from D-059, which is my entry and asserted it without measuring it, so
the error originated here and propagated outward. Their argument from it was
their own and was worth acting on; attributing the bad premise to them, as an
earlier wording of this paragraph did, moved my error onto a peer. At C465 the VMH holds 45 clusters and 33,694 cells, 8.8% of
the atlas, effectively tied with the arcuate's 45 and 37,667 (9.8%). The thin
region among populations carrying records is the **PVH**: 16 clusters, 8,710
cells, 2.3%. This is worth more than a correction — `vmh_sf1` is refused for
`Nr5a1` peaking at 0.391, and that cannot be blamed on sampling when the VMH is
the second-best-sampled named region in the atlas. The cause is transcript
abundance, which is exactly D-065's "a Cre line is not a marker".

I then tested the obvious fix anyway, since the top genes by variance across all
465 clusters are `Plp1`, `Mbp`, `Apoe` and `Cst3` — myelin and glial markers,
i.e. the gene set was partly measuring the neuron/glia axis while every
population under test is a neuron. Restricting variance to the 320 neuronal
clusters changes the gene set completely (top genes become `Avp`, `Sst`, `Gad2`,
`Dlk1`) and changes the similarity structure **not at all**: spread 0.081 →
0.088, median 0.942 → 0.941, and eight of nine nearest neighbours identical.

So **the gene set stands exactly as D-066 fixed it.** A pre-registered parameter
that survives a genuine attempt to improve it is worth more than one quietly
replaced, and recording the failed amendment is the point.

### Two changes, both before any delta

**1. The metric becomes cosine after centring each gene across the populations
being compared.** Declared reason: an un-centred cosine over non-negative data
is a known-degenerate choice, and D-066's own stated purpose for restricting the
gene set — measure cell identity, not what is high in everything alike — is the
same argument one level down. Gene centring is chosen over z-scoring because it
removes the shared baseline without reweighting every gene by its variance
across nine populations, which with nine points is itself noisy.

The selection among the four rows above used **spread and top-1 margin only**.
No direction label enters either quantity, so this cannot tilt the primary
toward similarity working. It can only decide whether the feature is able to
rank at all.

**2. A per-row tie guard, which is what the `tied` proposal was for.** A row is
predicted only if its nearest neighbour is **stable**: resampling the 2000 genes
with replacement, B = 200, the same population must be nearest in **at least
95%** of replicates. A row failing that is reported **unranked** — counted as
unpredicted, exactly like D-066's overlap fallback, and never counted wrong.

This is the mechanism for the standing requirement that this system be able to
answer *"insufficient evidence to rank these interventions reliably."* Until now
that requirement had no implementation in this experiment; a `k = 1` ranking on
a 0.0051 margin would have asserted a neighbour it could not defend.

### The stopping rule, which matters more than either change

Counting D-059's three amendments, this is the fifth time this experiment has
been re-specified before being run. Every one was forced by a measurement of the
feature and none by a delta, which is the distinction that keeps them
legitimate — but the pattern itself is the risk, because a sufficiently
determined amender can reach any result without ever looking at one.

**So: this is the last amendment. The evaluation runs next.** If the tie guard
leaves too few ranked rows to measure, the finding is *"transcriptomic
similarity cannot rank these populations reliably, and the experiment reports
insufficient evidence"* — which is a real answer, is reported as the result, and
is **not** grounds for a sixth re-specification. No further change to the gene
set, the metric, `k`, the tie guard, the mapping, or the admitted tiers will be
made before the primary is measured, and any change after it is a new experiment
under a new number with the first result reported alongside.

### Unchanged

D-059's predictions 1 and 2 stand and are unmeasured. D-066's overlap
disqualification, `k = 1`, the exploratory `k = 3` that cannot touch the primary,
tier A admission, and every logged diagnostic stand as written. This remains a
Tier B feature that is never a training target.

## D-068 - The result: transcriptomic similarity cannot rank these populations, and D-059's primary is unevaluable by its own criterion

**This is the outcome D-067's stopping rule named in advance, and it is
reported as the result rather than re-specified.** No delta was computed.

### The measurement

Under leave-one-population-out on the 9-population tier A pool, with genes
centred on the training fold alone, each held-out population's nearest
non-overlapping neighbour was tested for stability against resampling the 2000
genes (B = 200, seed 20670, threshold 95%):

| held out | circuit records | nearest | stability | rankable |
|---|---|---|---|---|
| `arc_agrp` | 123 | `arc_pomc` | **52.0%** | no |
| `arc_pomc` | 52 | `arc_ghrh` | 78.0% | no |
| `lha_orexin` | 16 | `pvh_trh` | 79.0% | no |
| `pvh_trh` | 8 | `pvh_crh` | **99.5%** | yes |
| `pvh_crh` | 7 | `pvh_trh` | **100.0%** | yes |
| `lha_mch` | 3 | `lha_orexin` | 79.5% | no |
| `arc_kisspeptin` | 2 | `arc_agrp` | 38.0% | no |
| `arc_ghrh` | 1 | `arc_pomc` | 83.5% | no |
| `tmn_histamine` | 1 | `arc_pomc` | 68.0% | no |

**2 of 9 populations are rankable, holding 15 of 213 circuit records.** The
population carrying 123 of those 213 records has a nearest neighbour that
survives resampling barely half the time.

The figures above use fold-wise centring. An earlier pass centred genes across
all nine populations including the held-out one — no label leakage, since atlas
expression is not the label, but the wrong protocol for leave-one-population-out.
Correcting it moves individual stabilities (`lha_mch` 93.0% → 79.5%,
`lha_orexin` 51.5% → 79.0%) and **changes the verdict not at all**: 2 of 9
either way. The conclusion does not rest on the error.

### Why the primary is unevaluable rather than failed

D-059's prediction 1 requires **both** a population-clustered bootstrap interval
excluding zero **and** a row-level exact McNemar test under 0.05. Two populations
supply every ranked row, and a clustered bootstrap over two clusters admits
three distinct resamples; the interval cannot be formed. So the primary is not a
miss, a null, or a delta near zero — **it cannot be computed by the criterion
that was pre-registered for it**, and reporting any delta from 15 rows across
two mutually-nearest populations would assert a measurement this design
explicitly forbids.

Calling that a failure would be as wrong as calling it a pass. D-059 expected a
miss and recorded the mechanism; this is a third outcome neither entry
anticipated, and the honest name for it is the one the standing requirement
already had: **insufficient evidence to rank these reliably.**

### What was nevertheless learned, and it is not nothing

**The feature carries anatomical signal and no usable subtype signal.** Both
stable pairs are same-nucleus and mutually nearest: `pvh_trh` ↔ `pvh_crh` at
99.5/100%. The near-misses are also same-nucleus — `lha_mch` → `lha_orexin` at
79.5%. Every cross-nucleus choice is unstable. So population mean expression
locates a nucleus reliably and distinguishes cell types within and between them
too weakly to act on, which is the same boundary D-065 found from the other
side: anatomy identifies a nucleus, not a cell type.

**D-059's mechanistic prediction is corroborated as a tendency while being too
weak to test.** It predicted in writing that a held-out `arc_agrp`'s nearest
neighbour would be `arc_pomc` — same nucleus, shared transcriptome, opposite
effects on food intake. It is, in 52% of replicates, the modal choice. But 52%
is a coin flip, so whether similarity transfers the *sign* of an effect remains
untested. The predicted pair is real; the ranking that would have let us
interrogate it is not stable enough to use.

**The negative result is about this corpus, not about the idea.** Nine
populations, six of them contributing three records or fewer, cannot support a
nearest-neighbour experiment whatever the metric. The binding constraint is
curation depth, not the atlas: HypoMap supplied usable vectors for all 23 mapped
populations (1460–1679 live genes of 2000, median 1640), and the diagnostic a
second session asked for found no population with a thin vector. What is missing
is records, in populations other than `arc_agrp` and `arc_pomc`.

### What this forecloses, and what it does not

D-059's prediction 2 — that the similarity model must not damage the leave-one-
study-out result D-055 established — is **vacuously satisfied and therefore not
evidence of anything**: a model that declines to rank 93% of rows cannot damage
a result, and it is recorded as vacuous rather than passed.

Per D-067, no sixth re-specification follows. The route forward is **more
curated records in under-represented populations**, not a different metric: the
next version of this experiment is a new entry under a new number, run only when
the pool can support a clustered interval, and it must report this result
alongside its own.

The mapping of D-065 stands and remains Tier B, never a training target. Nothing
here is evidence that any intervention does anything.

> **Sensitivity, measured after the fact: the conclusion is robust and the table
> is not.** A second session's finding in D-070 — that for a tool which
> *proposes* candidates, a false one is reviewable while a missed one is
> invisible — prompted an audit of D-065's keep threshold, which records what it
> kept and nothing about what it nearly kept.
>
> > **The asymmetry does not transfer, and this paragraph first said it did.**
> > That argument holds because a curator reads every unit a matcher proposes.
> > Nothing reads this mapping before it becomes a mean-expression vector, so
> > here a wrong **inclusion** is the expensive error and a near-miss exclusion
> > the cheap one — the reverse. Stated as it was, it argues for loosening this
> > threshold, which would be wrong. The session that supplied the asymmetry
> > corrected its direction for this module, and the general form is: **prefer
> > the error a reviewer will see, and that depends on whether a reviewer
> > exists.** See D-073, which acts on it. Audited, 7 of 23 mapped populations have
> a region candidate within 0.35–0.50 of their best, now recorded as
> `near_miss_excluded` in the artefact. One of them bears directly on this entry.
>
> `pvh_trh` excluded **`C185-12: Ebf1.Trh.GLU-2`** — Trh-named, 100% PVH, mean
> 1.234 — at 0.442 of its best against a 0.50 cut, missing by 0.058. And
> `pvh_trh` is one of the only two rankable populations above.
>
> Adding that single cluster to `pvh_trh` takes the result from **2 of 9 rankable
> to 0 of 9**: `pvh_trh` → `pvh_crh` falls from 99.5% to 88.0% and `pvh_crh` →
> `pvh_trh` from 100% to 77.0%, so the one stable pair dissolves and no
> population clears the guard.
>
> **So "2 of 9" rests on a threshold decision, and it is fragile in one
> direction only — toward zero.** The finding this entry reports, that
> transcriptomic similarity cannot rank these populations reliably, holds under
> both: at 2 of 9 the clustered interval cannot be formed, and at 0 of 9 there is
> nothing to form it from. The conclusion survives; the specific table should not
> be quoted without this note. It also means the `pvh_crh` caveat I flagged
> elsewhere was understated — a refutation is not needed to reach zero, a 0.058
> shift in one threshold does it.
>
> The mapping is **not** changed. D-067 forecloses altering it after the primary
> is measured, and a revision would be a new experiment under a new number
> reporting this result alongside. This is a robustness note, not a
> re-specification.

## D-069 - The re-run trigger, declared while no EEG record exists and the answer cannot be known

**The user has approved the EEG curation pass, conditional on the matcher fixes
landing.** That pass is projected to add 150–170 records across `lha_orexin`,
`pvh_crh`, `lha_mch`, `tmn_histamine` and `pvh_oxt` — the five populations whose
thinness made D-068's primary unevaluable.

D-068 said the experiment re-runs "when the pool can support a clustered
interval" and **did not define what that means.** That is an undefined threshold
that would otherwise be judged after the records exist, by someone who can see
whether the experiment would pass — the exact failure D-048, D-050 and D-056 were
each written about. So it is defined here, now, while the pass has produced
nothing and the answer is unknowable. A second session made this argument about
the overlap rule in D-066 and it was right: a decision taken while it provably
cannot affect the result is the cheapest this project can buy, and it is
unavailable afterwards at any price.

### The trigger

The similarity experiment re-runs **once**, under a new number, only when all
three hold on the circuit subset:

1. **At least 8 populations rankable** under D-067's tie guard — nearest
   neighbour stable in ≥95% of 200 gene-resampled replicates, fold-wise centred.
2. **No single population holds more than 50% of ranked rows.** D-068's deeper
   defect was not the count but the concentration: `arc_agrp` carried 123 of 213
   records and had a coin-flip neighbour, so the primary was close to a test of
   one pair whatever the statistics said.
3. **At least 5 rankable populations hold ≥5 records each**, so the clustered
   interval rests on several populations that can each move it rather than on
   one that cannot.

Checked on the pool and the feature only. No delta is computed until all three
pass.

### What a near-miss means, stated before one can happen

**A near-miss is not grounds for relaxing any of the three.** If the pass
delivers 7 rankable populations, or 8 with one holding 55% of rows, the finding
is that the corpus still cannot support the experiment — reported in those
words, with the numbers, as D-068 was. Loosening a threshold to admit a pool
that just failed it is the goalpost move, and it is more tempting after a
multi-day curation effort than it was before one, which is precisely why the
threshold is fixed while nobody has spent that effort yet.

Everything D-066 and D-067 fixed carries over unchanged: the 2000-gene set,
gene-centred cosine, `k = 1`, the exploratory `k = 3` that cannot touch the
primary, tier A admission, overlap disqualification, the tie guard and all logged
diagnostics. A change to any of them makes it a different experiment under a
different number, not a re-run of this one.

### A leakage control the approval creates

The curation that follows is **sampling, not labelling** — choosing which papers
to read is not choosing what they say, and the records remain traceable to
verbatim quotes under the existing discipline. But a session now knows which
populations need records for the experiment to become evaluable, and knows the
similarity structure, because I sent them the nearest-neighbour table.

So: **curation for this pass must not consult the similarity matrix, the
nearest-neighbour table, or the stability figures.** Direction labels come from
the passage and nothing else. The populations to target are chosen by record
count, which is public and has no bearing on what any record says; which
direction a record reports must never be informed by what the feature expects.
Recorded because the opportunity for this particular leak did not exist before
today, and the person best placed to avoid it is the one who can no longer
un-know the table.

### What this does not promise

Nothing here predicts the re-run will pass, and D-059's expectation that
similarity fails to transfer the sign of an effect stands unmeasured. The three
conditions make the question *askable*; they do not make the answer likely
either way. D-065's mapping remains Tier B and never a training target.

---

## D-070 - The matcher lost half of every joint sleep statement, and my scope guard was narrower than the rule it guarded

**Found while costing the EEG curation pass, by running the matcher over the
passages it would actually be used on.** Not by a test, and not by inspection —
both defects are invisible in the code and in the vocabulary, and appear only in
the output.

The session curating the genotype dimension flagged the prerequisite: all 946
proposed sleep units carry `candidate_outcome_id: "sleep"` and predate D-064, so
a curator accepting the candidate id on an NREM measurement would land it in the
pooled cell and nothing would error, because `sleep` is a legitimate id. True,
and regenerating the units was going to fix it. Then I ran the live matcher over
those 946 passages.

### Defect 1: each ordering silently lost its elided half

`TermMatcher` is a word-boundary alternation over synonyms. So:

| passage | matched | missed | units |
|---|---|---|---|
| "NREM and REM sleep both decreased" | `rem_sleep` | **`nrem_sleep`** | 47 |
| "REM and NREM sleep were unchanged" | `nrem_sleep` | **`rem_sleep`** | 12 |

Neither `NREM sleep` nor `non-REM sleep` occurs contiguously in the first, and
`REM sleep` does not in the second. The commonest composite phrasing in this
literature, and the elided half was dropped without trace.

The fix is **data, not code**: `NREM` and `non-REM` added to `nrem_sleep`, `REM`
to `rem_sleep`. It is reviewable in a diff, and it needed no change to the
matcher — which matters, because a change to `_boundary_pattern` applies to
every vocabulary in the project.

**And I got this wrong once before getting it right.** I first withheld bare
`REM`, reasoning that `\bREM\b` fires inside `non-REM sleep` and would claim
every NREM passage. That was the wrong frame. This matcher **proposes
candidates a curator accepts or rejects**, so a false candidate is reviewable by
the person who sees it and a missed one is invisible to everybody. Withholding
the synonym to avoid a false positive lost REM on 12 units outright.

Then I measured the cost I had been avoiding: 25 units mention non-REM, and
**all 25 also discuss REM in its own right**, so the spurious match never fires
alone on this corpus. The precision cost is bounded and visible; the recall cost
was neither. Left as a measured, pinned overlap rather than patched with a
negative lookbehind — the D-049 pattern of annotating a trap instead of
resolving it automatically.

### Defect 2: an out-of-scope readout arrives under a legitimate id

"EEG delta power during wakefulness increased" resolves to `wakefulness`,
through an incidental word. **133 of 946 units are in this state**, and 224
carry an out-of-scope readout at all.

**That number is a correction, and the way it was wrong is the point.** I first
measured 142 with a regex I hand-rolled for the purpose; the other session
measured 125 with one it hand-rolled independently. Neither was right. Running
the detector that reads the vocabulary's own `out_of_scope_readouts` gives 224
carrying and **133 still resolving to an in-scope state**, split
`neural_activity` 179 and `behavioural_affective` 52 — and both sessions
reproduce it exactly, because there is now one detector instead of two. Two
hand-rolled matchers disagreeing in opposite directions around the right answer
is the signature this project already knows: `syn4rch.curate.units` imports the
screener's matchers rather than reimplementing them for exactly this reason, and
I reimplemented one anyway to measure a defect about reimplemented matchers.

D-063 put neural activity out of scope and `tests/test_vocab_scope.py` guards
the *vocabulary*: no entry may list an out-of-scope term as a synonym of any
kind. That guard could not catch this and was never going to. **No rule is
broken on the way in** — the id is legitimate, the synonym is legitimate, and
the passage simply contains both an in-scope state and an out-of-scope measure.
My guard was narrower than the rule it was written for, and the gap is exactly
one indirection wide.

`syn4rch.harvest.screen.out_of_scope_readouts` now names the categories a
passage mentions. It **does not suppress the id**, deliberately: a sentence can
report both, and deciding which the record is about is a reading judgement the
vocabulary refuses to make for a curator (D-021). It prompts the judgement
rather than making it.

### A third thing I reported as a defect and was wrong about

I told the other session that 535 units resolving to a state id *and* still
`sleep` was a defect, "the pooled id stays on the menu for exactly the passages
it should not". That is backwards. Under D-064 a joint report *belongs* on
`sleep`, and the proposer offering all three ids is what lets a reader follow
that rule. Suppressing `sleep` would make the choice for them. Corrected here
rather than left standing, because the count was quoted.

### One change outside my own paths

`syn4rch.vocab.core.meta(kind)` returns a vocabulary's `_meta` block. `load`
deliberately exposes only `entries`, and the alternative was reimplementing this
module's path resolution inside `harvest/` to read the same file a second way —
which is how two views of one vocabulary come to disagree. `vocab/` is assigned
to neither session under the current ownership split; flagged to the other
session rather than taken silently.

### What this does not claim

* It does not re-propose anything. The units are regenerated against outcomes
  0.5.1; no record changed, no direction moved, and the dataset is untouched.
* It does not make the EEG pass safe on its own. The 133 out-of-scope passages
  are now *flagged*, not excluded, and a curator still has to read them.
* `outcomes.json` goes to 0.5.1 — a synonym addition, not a new id.


## D-071 - Adversarial verification refuted none of the eleven tier A mappings, and found a bug in the tooling doing the verifying

**44 agents, three independent lenses per claim then an adversarial judge
defaulting to REFUTED: 11 of 11 judged, 0 refuted, 0 tier A name claims false,
11 carrying caveats.** No mapping was confirmed clean, and none fell.

That is worth stating against the baseline: an earlier pass over the
**hand-made** mapping refuted 3 of 13. The same attack on mappings produced by
D-065's declared rule refuted nothing. The rule is the difference, not the
attacker.

### The finding that matters is about my code, not the mappings

`nesting()` compared **cluster name strings**. Two populations mapped at
different levels share cells without sharing a name, so it was structurally
blind to them. Recomputed from the atlas's own per-cell labels, with cell counts
reconciling exactly for all 23 mapped populations:

**11 cell-level overlapping pairs. One declared. Ten invisible, and all ten
share no cluster name at all.**

`arc_agrp` (C66) and `arc_tbx3` (C185) overlap on **4663 of `arc_agrp`'s 5244
cells — 88.9%** — and intersect on nothing a string comparison can see.
`arc_ghrh` is wholly inside `arc_th`; `pvh_trh` wholly inside `pvh_sim1`.

The docstring of the broken function said downstream similarity ranking "has to
be able to see the overlap rather than rediscover it as a finding". It could
not. **My own tests passed throughout** — `test_no_unrecorded_overlap_exists`
compared names too, so it asserted the same blindness it was written to prevent.
An adversarial agent found this by resolving every mapping to its finest-level
leaves and reconciling cell counts, which is a check I never wrote.

### D-068 survives, and that is luck rather than design

D-066 disqualifies an overlapping neighbour, and D-068's run logged a Jaccard of
0.00 on every row. Those zeros were **correct**: the nine tier A populations
holding records are cell-level disjoint, verified. So the published result stands
unchanged and needs no re-run.

But the mechanism protecting it was not working, and it would have failed on the
tier A+B+C secondary D-066 requires reported regardless of the primary:
**`pvh_trh` lies entirely inside `pvh_sim1`** — 2151 cells, 100% of `pvh_trh` —
and `pvh_trh` is one of the only two rankable populations. With tier C admitted,
a name-based test would have let a population supply a prediction for a
population it wholly contains. That secondary has not been run, so the fix lands
before the defect could reach a result.

### The fix

`Assignment` and `PopulationVector` carry `leaf_clusters`, the finest-level
clusters a mapping resolves to; `nesting()` and D-066's disqualification compare
those. Names remain the fallback for same-level comparisons, where names and
cells agree. The artefact now records `leaf_clusters` per entry and declares all
**11** pairs. Tests assert the cross-level case explicitly, and the fixture is
required to share no cluster name — otherwise it would stop demonstrating what
name matching missed.

### Two things the judges got right that their own examiners got wrong

Recorded because it is evidence the panel's structure did work:

* Two lenses reported that the ARC Npy+/Agrp− branch was "silently dropped".
  The judge checked and it is not: `C286-179: Ptprk.Npy.Sst.GABA-4` (1433 cells,
  100% ARC, Sst 4.37 against Agrp 0.18) is mapped as `arc_sst`, descends from
  `C66-47` and not `C66-46`, and shares zero cells with `arc_agrp`. The real
  issue is narrower — `Npy`-resolved evidence in the ARC is ambiguous between
  two entries, and neither records it.
* Two lenses called `best_marker_rank: 5` wrong because `Npy` ranks 2. The judge
  read the module: it is the best rank of **the chosen marker across kept
  clusters**, not across declared markers. Not a defect.

### A caveat with no fix available here

`arc_agrp`'s mapped cluster is chemically heterogeneous: `C465-290:
Lepr.Agrp.GABA-4`, 581 of 5244 cells (11.1%), has summed mitochondrial mean
**3.60** against 18.29/18.45/21.25/21.62 for its four siblings and an atlas
median of 20.72 — a single-nucleus prep among whole-cell ones. It is still 100%
ARC with Agrp 2.23, so these are AgRP neurons on different chemistry, but the
headline mean of 3.4254 is a blend across library preparations. Recorded rather
than corrected: dropping the cluster would change the mapping after the primary
was measured, which D-067 forecloses.

The mapping is otherwise unchanged. Tier B, never a training target.

---

## D-072 - Pre-registration: the EEG curation pass, sized on the regenerated units, with the leakage control the approval created

**Written before any paper is read**, and committed before the briefs. The user
authorised the pass after the D-070 matcher fixes landed; this fixes what it
will curate, what would count as working, and what a reader is allowed to know.

### What regenerating the units produced

`syn4rch.curate.run units` re-segmented all 974 retrieved papers against
outcomes 0.5.1. Units **8,860 -> 9,836**, and the sleep bucket stopped being one
bucket:

| outcome id | before | after |
|---|---|---|
| `sleep` | 946 | 767 |
| `wakefulness` | 0 | **416** |
| `rem_sleep` | 0 | **356** |
| `nrem_sleep` | 0 | **268** |
| `sleep_fragmentation` | 0 | **84** |
| `sleep_latency` | 0 | **31** |
| total | 946 | **1,922** |

The total more than doubles because a passage reporting NREM, REM and wake now
yields three units where it yielded one. Those are three readouts, and the 767
remaining on `sleep` are genuine composites.

**290 units were lost, and both causes are correct.** 179 were `sleep` units
matching only on bare `arousal`, which D-070 deliberately made a synonym of
nothing — three papers lost every unit they had, and reading their passages
confirms the verdict: all three are about GCaMP signals, circuit tracing and
microcircuit anatomy, and none reports a vigilance state. An ambiguous word was
carrying them.

**And 111 were stale for two vocabulary versions.** They carried
`cognitive_decline`, an id renamed to `cognitive_performance` at outcomes 0.2.0
by D-029. The new file has `cognitive_performance` 111 — the same count, so no
readout was lost and the rename simply had never reached the artefact. **The
units file had been carrying an outcome id the vocabulary did not have, for the
life of the project, and nothing checked.** That is the stale-artefact class
again: the fourth instance after `lab_key`'s affiliation (D-056), `screened_on`
(D-061) and the signaling graph. An invariant — every outcome id in the units
file must exist in the vocabulary — is cheap and now tested.

### The seam, re-measured

| population | tier | records now | units before | units after |
|---|---|---|---|---|
| `lha_orexin` | A | 16 | 306 | **560** |
| `lha_mch` | A | 3 | 129 | **301** |
| `pvh_crh` | A | 7 | 146 | **269** |
| `tmn_histamine` | A | 1 | 58 | **131** |
| `pvh_oxt` | A | 0 | 29 | 33 |
| `arc_agrp` | A | 123 | 8 | 8 |
| `arc_pomc` | A | 52 | 4 | 4 |

Seam total **945 -> 1,724**. Tier A units outside `arc_agrp` and `arc_pomc`:
**668 -> 1,294**, against **38** records in those populations today. The two
dominant populations still draw 12 units between them, so the complementarity
survived the matcher fix and strengthened.

91 papers, 59 with at least one proposed unit, 42 with ten or more.

### The leakage control, which the approval created

Knowing which populations the similarity experiment needs is new as of today,
and the session that found it raised this before I could act on it. Choosing
which papers to read by record count is **sampling** and is legitimate. What a
record *reports* must come from the passage alone. So:

1. Readers receive the passage, the locator, the vocabularies and the protocol.
   They do **not** receive the similarity matrix, the nearest-neighbour table,
   the stability figures, or any statement of which populations are wanted.
2. The targeted populations are named **here**, in this entry, and nowhere a
   reader will see them. That makes the sampling auditable and the readers
   blind at once.
3. No population name appears in a brief as a thing to look for. A reader told
   "we need `pvh_crh` records" will find `pvh_crh` records, and every one of
   them will be individually defensible, which is why no downstream
   quote-level check catches it.

**Targeted, for the record:** `lha_orexin`, `lha_mch`, `pvh_crh`,
`tmn_histamine`, `pvh_oxt` — the five tier A populations the seam reaches that
hold 27 records between them.

### The predictions

1. **Volume.** At least **150 records** admitted from the 91 papers. The circuit
   pass yielded 146 from 60.
2. **Depth where it is missing.** At least **6** tier A populations holding 5 or
   more records, against 3 today (`lha_orexin` 16, `pvh_trh` 8, `pvh_crh` 7).
3. **Concentration.** `arc_agrp` + `arc_pomc` fall below **65%** of tier A
   circuit records, from 82%.
4. **Replication.** Circuit cells spanning two or more laboratories rise by at
   least **15**, from 26. The seam brings 68 laboratories absent from the
   corpus.

**And the prediction I expect to fail, recorded because D-051 taught me to.**
Number 4. D-064 split `sleep` into five ids, so an NREM record now only matches
another NREM record — the comparability the split bought is paid for in
fragmentation, and this pass is where the bill arrives. Two laboratories both
reporting "sleep decreased" used to share a cell; now one reporting NREM and one
reporting wakefulness do not. I expect volume and depth to come in and
replication to disappoint, and if that happens the finding is that **splitting
an outcome id and counting replication across it are in tension**, which is
D-023's trade appearing a second time at a larger scale.

### What would abandon the pass

* Fewer than 60 records after the first 30 papers — the yield is not there.
* Any reader brief found to name a targeted population, which voids the leakage
  control for that reader's output.
* `pvh_crh` refuted by the adversarial mapping pass *and* the pass yielding
  fewer than 10 `pvh_crh` records: the second-largest seam would then be
  supporting neither question.

### Named limitations, before anyone finds them

* **My screen is conservative and this entry does not fix it.** Requiring EEG
  mentions plus intervention plus hypothalamic terms misses 50 studies that
  carry sleep units, found by the other session using a different filter. The
  limitation is about the screen's history rather than this pass: the same
  threshold has run on every seam this project has ever selected.
* **133 of the sleep passages carry an out-of-scope readout under an in-scope
  id** and are flagged, not excluded (D-070). A curator still has to read them.
* **`pvh_crh` is unsettled.** The adversarial pass on the eleven tier A mappings
  has not reported, and `pvh_crh` is both the second-largest seam here and the
  only stable neighbour of `pvh_trh` in D-068.
* The pass cannot satisfy D-069's re-run trigger by itself. That requires 8
  rankable populations, none over 50% of ranked rows, and 5 rankable with 5+
  records — and rankability is a property of the mapping and the tie guard, not
  of record counts. This pass moves one input of three.



## D-073 - A tier that cannot survive a 0.058 shift: the threshold-stability veto, and the result it would have given from the start

**D-067 forecloses changing the mapping or the admitted tiers, and says any
change after the primary is measured is a new experiment under a new number with
the first result reported alongside.** This is that entry. D-068's conclusion is
unchanged; its table is superseded.

### What was wrong

`pvh_trh` held tier A while excluding `C185-12: Ebf1.Trh.GLU-2` — Trh-named,
100% PVH, 3034 cells — at 0.442 of its best against a 0.50 cut. D-068's
sensitivity note already recorded that adding that one cluster takes the result
from 2 of 9 rankable to 0 of 9. **A tier whose membership cannot survive a 0.058
shift in an arbitrary cut was never load-bearing**, and recording the near miss
fixed the reporting while leaving the tier intact.

A second session measured the exposure and it is concentrated rather than
diffuse: of 19 near misses, 14 sit on tier C populations and 4 on tier B — **one
is on tier A**, and it is the one that collapses the table. They also noted that
0.442 is sixth of nineteen in a band running 0.361 to 0.491, so it is mid-pack:
if 0.442 is "nearly in", so are five others.

### The rule

**Tier A additionally requires the kept set to be invariant to `keep_fraction`
across a declared band, `[0.40, 0.60]`.** Stability is a veto on tier A, not
another axis of support: a mapping whose very membership depends on where the
cut falls is not established, whatever its names and contrast say. It cannot
push a population below B.

Both directions count, and the second matters more. A candidate that would
**enter** at the low end is a cluster the vector may be wrongly missing. A member
that would **leave** at the high end is one the vector may be wrongly carrying —
and that is the expensive error here, because nothing reviews this mapping before
it becomes a mean. This is the corrected asymmetry: prefer the error a reviewer
will see, and there is no reviewer in this path.

The band is declared rather than fitted. A threshold chosen to capture 0.442
would be calibrated to the case under discussion; a stability requirement over a
band is the same discipline as D-067's tie guard one level up, and it does not
fire on a single case: **7 of 23 mapped populations are threshold-unstable.**

### What it changes

| | |
|---|---|
| `arc_ghrh` | A → B, destabilised by `C465-285: Gsx1.Ghrh.GABA-3` at 0.583 (would leave) |
| `pvh_trh` | A → B, destabilised by `C185-12: Ebf1.Trh.GLU-2` at 0.442 (would enter) |

Tiers are now A 9, B 6, C 8, refused 7. **`arc_ghrh` was caught by the rule and
not by anyone's attention** — the near-miss audit only looked at the inclusion
direction, so a symmetric rule found a second tier A failure that no one had
named. Both mappings are unchanged and both survived D-071's adversarial pass;
only the tier label moves.

### Re-measured, and the conclusion strengthens

The primary pool falls from 9 populations and 213 circuit records to **7 and
204**. Re-running D-067's tie guard on it, with D-071's cell-level overlap test:

**0 of 7 rankable. Not one population has a stable nearest neighbour.**

| held out | records | nearest | stability |
|---|---|---|---|
| `arc_agrp` | 123 | `arc_pomc` | 58.0% |
| `arc_pomc` | 52 | `arc_agrp` | 55.5% |
| `lha_orexin` | 16 | `pvh_crh` | 58.5% |
| `pvh_crh` | 7 | `arc_kisspeptin` | 48.0% |
| `lha_mch` | 3 | `lha_orexin` | 64.5% |
| `arc_kisspeptin` | 2 | `arc_pomc` | 57.0% |
| `tmn_histamine` | 1 | `arc_pomc` | 70.0% |

D-068 reported 2 of 9 and "insufficient evidence to rank". The corrected answer
is the same finding with nothing left to qualify it: there is no clustered
interval to form because there is no ranked row at all. The other session
predicted exactly this — that the honest table was 1 of 9 or 0 of 9 from the
start — before it was computed.

**One thing is sharper than before.** `arc_agrp` and `arc_pomc` now name each
other reciprocally, at 58.0% and 55.5%. D-059 predicted in writing that this
pair would be mutually nearest and that their opposite effects on food intake
were the reason expression similarity should fail to transfer a sign. The pairing
is confirmed in both directions and is a coin flip in both, so the mechanism is
corroborated and remains untestable for want of a stable ranking.

### Scope

No parameter of D-066 or D-067 moves: gene set, centred cosine, `k = 1`, the
exploratory `k = 3`, tier A admission, overlap disqualification and every logged
diagnostic stand. D-069's re-run trigger stands and now reads against corrected
tiers, which is why this had to land before the EEG pass produces a record
rather than after. Tier B, never a training target.
