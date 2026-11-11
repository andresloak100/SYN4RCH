# What the engine's ingestion needs from the trace download

**For the owner of `analysis/axon_trace_download.py` (the other agent). Written 2026-11-11 by Claude, before any
Stage A trace is consumed. The downloader and its manifests and receipts are not touched from this side.**

The engine reads traces only through one listing file. It never scans directories and never opens a file that is
not listed. If the downloader's final manifest has these fields under other names, say so and the listing below will be
written from it by a small converter on the engine's side; the raw files stay where they are, byte for byte.

## Files

* One file per neuron, raw bytes as served, never rewritten, renamed inside, mirrored or resampled.
* Path, relative to `data/raw/jiao2025/traces/`: `<part>/<brain sample>/<neuron>.swc`, for example
  `validation/192711/192711_005.swc`. This is section 7 of the locked plan. `<part>` must be the part the locked split
  gives that neuron (`reports/axon_dataset_manifest.tsv`, column `part`); the engine refuses a file stored under another
  part. Stage A holds `train` and `validation` only. **No file of a test brain sample may exist on disk before the freeze.**
* The 37 traces already in `data/raw/jiao2025/swc/` (all training) stay as they are. A neuron may have a second copy
  from the portal; the engine prefers the true-hemisphere copy and checks the two against each other.

## Listing: `data/raw/jiao2025/traces/MANIFEST.tsv` (tab-separated, one header line)

| column | required | meaning |
|---|---|---|
| `file` | yes | path relative to `data/raw/jiao2025/traces/`, ending in `.swc`, first component = part |
| `neuron` | yes | cohort neuron ID exactly as in the cohort manifest: `<brain sample>_<number>`, for example `192711_005` |
| `bytes` | yes | size of the stored file |
| `sha256` | yes | SHA-256 of the stored bytes; the engine re-hashes on open and refuses a mismatch |
| `source` | yes | the exact address the file came from |
| `frame` | yes | `portal_true_hemisphere` for portal files (they keep the true hemisphere; the engine mirrors right-hemisphere neurons itself). A neuron the source does not serve gets one row with `frame` = `unavailable_at_source` and no file |
| `licence` | wanted | the notice the source shows |
| `how` | wanted | fetch method, time (UTC) and HTTP status |

Rows whose `file` does not end in `.swc` are ignored, so receipts may share the listing.

## Things only the downloader can tell the engine

1. **Every neuron accounted for.** For each of the 5,633 Stage A neurons (4,448 training, 1,185 validation; 37 training
   traces are already on disk): either a file row or an `unavailable_at_source` row with the reason (HTTP status,
   empty body, HTML instead of a trace). Never silently absent.
2. **Completeness marker.** A way to know the download is finished and the listing is final (for example a final
   manifest SHA-256 in a receipt). The engine pins the listing's SHA-256, so a listing that still changes cannot be used.
3. **No test file.** A statement in the final receipt that no address of a test brain sample was requested.
4. **Duplicates.** If the same neuron was fetched twice with different bytes, both hashes and which one is listed.

## What the engine checks itself after that (no action needed from the downloader)

Parsing, node types (axon coded as type 2 or type 0), one soma, hemisphere against the portal index, soma against the
paper's table (after the frame rule), tree integrity, at least 1,000 um of axon. These are the mechanical rules of plan
section 9; none looks at a model or a score. Excluded neurons are logged with the rule number and counted per part.
