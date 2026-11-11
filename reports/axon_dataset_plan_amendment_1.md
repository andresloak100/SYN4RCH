# Amendment 1 to the locked cohort plan: no clipping; complete measured trees are scored

**Dated 2026-11-10, before any validation or test trace exists on disk and before any biological scoring.**
It amends `reports/axon_dataset_plan.md` (locked by commit c96b167). It does **not** change the animal split,
the classes, the counts, the sealing rules or the query-safe contract.

## What changes

Section 9, rule 7 of the plan said: "Cable outside the atlas box (AP 0 to 13,200, DV 0 to 8,000, ML 0 to 11,400) is
clipped for scoring, as in the authors' own measurements, and the clipped length is logged."

It is replaced by:

> **7. No clipping.** A trace is ingested, used as a donor and scored as the complete measured tree. Truth,
> copying controls and candidate forecasts live in the same full coordinate domain (CCFv3 micrometres, mirrored-left
> frame), including cable that lies outside the atlas box. The length of cable outside the box is logged per neuron
> and reported per part; it is never removed, moved or reconnected.

## Why

* A read-only review by the other agent (`reports/axon_engine_codex_review.md`, with synthetic counterexamples in
  `data/interim/axon_engine_codex_review/clipping_counterexamples.json`) showed that the engine draft's clipping
  joined an exit point to a later re-entry point with a straight chord, adding cable that was never measured
  (400 um inside the box became 1,400 um), and dropped an edge whose two ends were outside the box although it
  crossed the box (13,300 um became 100 um). It also clipped truth and donors but not translated or attached
  forecasts.
* A clipped arbor is in general a forest. The forecast and metric format of this project and of the independent
  referee is one rooted tree with edges of at most 50 um. A forest cannot be stored in it without inventing
  connections, and inventing cable is not acceptable.
* The authors clamp coordinates to the atlas box only to look up region labels for their length table. The geometric
  measures used here (nearest-sample distance, precision and recall at 200 um) need no atlas lookup and are defined
  for any coordinates, so nothing requires clipping.
* On the 37 inspected traces, cable outside the box was rare: five descending axons continue past the caudal end of
  the atlas (to AP 14,207) and one neuron has three nodes about 5 um outside the ventral surface.

## Consequences

* The scoring domain is the complete measured tree, for every method alike. This is stated in the engine protocol.
* The engine's clipping function is removed, and the two counterexamples become regression tests: ingestion must
  return exactly the measured cable.
* Rule 6 (at least 1,000 um of axon cable) now refers to the complete tree.
* Nothing was fitted, selected or scored on validation or test animals before this amendment; no trace of theirs is on
  disk.
