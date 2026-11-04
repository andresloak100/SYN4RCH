# Curation gate: rejections

Generated 2026-11-04T21:36:18Z by `python3 -m syn4rch.curate.run gate`.

| input | value |
|---|---|
| proposals | `data/interim/proposals_circuit/all.jsonl` |
| retrieved full text | `/Users/andresarriaga/Documents/SYN4RCH/data/raw/fulltext` |
| full-text files parsed | 974 |
| of those, with parse warnings | 39 |
| of those, yielding no passages | 0 |
| proposals read | 96 |
| proposal lines unreadable | 0 |
| admitted | 96 |
| rejected | 0 |
| locator-passage agreement | blocking (default) |

## The check

`support_quote` must occur as a substring of the paper's retrieved
text under exactly one normalisation: collapse whitespace runs to one space; fold unicode dashes, quotation marks and non-breaking spaces to ASCII; strip the ends. Nothing else: no case folding, no stemming, no dropped words, no fuzzy matching.

`source_locator` must name a passage that exists in that paper.

A failed proposal is rejected and recorded. It is never repaired,
re-quoted or nudged into passing (docs/CURATION_PROTOCOL.md 4).

## Rejections by reason

None. Every proposal read was admitted.

## Every rejection

None.

