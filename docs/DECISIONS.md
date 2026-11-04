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
