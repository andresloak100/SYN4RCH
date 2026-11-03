"""Literature harvest pipeline for SYN4RCH.

Modules
-------
``cache``       sqlite response cache, polite HTTP client, append-only raw log
``pubmed``      E-utilities ESearch/EFetch, XML parsed with xml.etree
``europepmc``   Europe PMC search, open-access detection, OA full-text retrieval
``crossref``    DOI metadata backfill only
``queries``     query matrix built programmatically from data/vocab/*.json
``screen``      transparent rule-based screener (include/uncertain/exclude + reasons)
``run``         CLI: ``python3 -m syn4rch.harvest.run --max-queries N [--resume]``

No third-party dependencies. No email address is sent to any API.
"""
