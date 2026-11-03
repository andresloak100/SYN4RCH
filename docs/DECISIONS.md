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
