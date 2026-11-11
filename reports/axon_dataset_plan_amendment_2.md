# Amendment 2 to the locked cohort plan: the access route is decided; a neuron the source does not serve

**Dated 2026-11-11, before any Stage A trace is on disk (the directory `data/raw/jiao2025/traces/` does not exist yet)
and before any validation or test trace has been opened by anyone.** It amends `reports/axon_dataset_plan.md` (locked
by commit c96b167; amendment 1 in c659eac). It does **not** change the animal split, the classes, the counts, the
sealing rules, the query-safe contract, inclusion rules 1 to 6 or amendment 1 (no clipping).

## 1. Access route (plan section 6 left this to the user)

On 2026-11-11 the user said: "download the missing neuron traces so we can run the real test against copying."
The other agent relayed it and took ownership of the download.

* **Route:** the portal's public per-neuron address, no account, no terms accepted on the user's behalf.
* **Stage A now:** training and validation only: 5,596 files not yet on disk (4,411 training, 1,185 validation), about
  6 to 8 GB. Raw bytes are stored unchanged, with a hash and a receipt per file.
* **Stage B later:** the 1,547 test traces stay unfetched until the method is frozen. No test file may be on disk
  before the freeze; the freeze records that.
* **Owner:** the other agent's `analysis/axon_trace_download.py` with its manifests and receipts. The engine starts no
  download of its own. What ingestion needs from the final listing is in `reports/axon_engine_ingestion_fields.md`.
* The licence position of plan section 6 is unchanged: private, non-commercial analysis; traces and anything derived
  from them stay local.

The earlier blocker "waiting for the user's choice of route" is closed. What remains before a real test is the
download itself, then ingestion, validation and the freeze.

## 2. New mechanical rule 8: a neuron the source does not serve

Plan section 9 had no rule for a trace that cannot be obtained, so one missing file would have blocked a whole stage.

> **8. Unavailable at the source.** A neuron for which the source serves no trace (the download listing records it as
> `unavailable_at_source`, with the reason) is excluded, logged under rule 8 and counted per part. Never silently absent.

* The rule is decided by the source alone. It never looks at a model, a forecast or a score, and it applies to every
  method alike.
* A training neuron under rule 8 is simply not a donor. A validation or test neuron under rule 8 is listed in the
  referee declaration's exclusions before any forecast is sealed.
* Reporting duty: the counts under rule 8 are given per part and per source label. If a whole held-out brain sample is
  lost, or more than 5% of a held-out part, the report says in its first lines that the scored cohort is smaller than
  the locked one.

## 3. Why this is safe to decide now

Nothing was fitted, selected or scored on validation or test animals, and none of their traces is on disk. Both points
are fixed before the first Stage A file is ingested.
