# FABRICATED TEST DATA -- NEVER EVIDENCE

Every file in this directory is **synthetic test data, invented to exercise
code**. None of it was extracted from any publication.

* No `study_id`, `doi`, `pmid`, title, author, affiliation, journal,
  population, marker gene, assay, quote, number or direction here refers to
  anything real. Each one was made up.
* Every record carries a `study_id` beginning with `synthetic:` so that a
  fabricated record can never be mistaken for a curated one, in a file, in a
  report, or in a `grep`.
* These files must never be copied into `data/curated/`, cited, aggregated,
  used to train or evaluate a model, or quoted as a finding about biology.
* The fixtures deliberately contain *defective* records -- missing quotes,
  unsupported numbers, planted conflicts -- because their purpose is to make
  the QC layer fail in the right places.

If you are looking for real data, there is none in this repository yet.
