"""Loader, normaliser and strict resolver for the SYN4RCH vocabularies.

Design rules (deliberate, do not "improve" without a DECISIONS entry):

* Matching is EXACT on a normalised string. Normalisation only removes
  presentational noise (case, whitespace, hyphen/underscore style, surrounding
  punctuation, a few Greek letters and the diacritics in "oestrous"). It never
  removes words and never truncates.
* No fuzzy matching. "POMC-adjacent neurons" must resolve to None, not to
  ``arc_pomc``. Substring and edit-distance heuristics are the standard way a
  curated vocabulary silently starts fabricating annotations, so they are absent.
* A synonym that would map to two different ids is a curation bug, not something
  to arbitrate at runtime: :func:`selfcheck` raises.

Standard library only; valid syntax for Python 3.10.
"""

from __future__ import annotations

import json
import os
import re
import unicodedata
from typing import Dict, Iterable, List, Optional

KINDS = ("populations", "outcomes", "methods", "regions")


class VocabError(Exception):
    """Base class for vocabulary problems."""


class VocabCollisionError(VocabError):
    """Raised when an id is duplicated or a synonym maps to more than one id."""


# --------------------------------------------------------------------------- #
# locating data/vocab
# --------------------------------------------------------------------------- #

_ENV_VAR = "SYN4RCH_VOCAB_DIR"


def vocab_dir() -> str:
    """Absolute path of the directory holding the vocabulary JSON files.

    Resolution order: ``$SYN4RCH_VOCAB_DIR``, then ``<repo>/data/vocab`` found by
    walking up from this file. Raises if neither exists, rather than guessing.
    """
    env = os.environ.get(_ENV_VAR)
    if env:
        if not os.path.isdir(env):
            raise VocabError("%s=%r is not a directory" % (_ENV_VAR, env))
        return os.path.abspath(env)

    here = os.path.abspath(os.path.dirname(__file__))
    cur = here
    for _ in range(8):
        cand = os.path.join(cur, "data", "vocab")
        if os.path.isdir(cand):
            return cand
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    raise VocabError(
        "could not locate data/vocab above %s; set %s" % (here, _ENV_VAR)
    )


# --------------------------------------------------------------------------- #
# normalisation
# --------------------------------------------------------------------------- #

_GREEK = {
    "α": "alpha",
    "Α": "alpha",
    "β": "beta",
    "Β": "beta",
    "γ": "gamma",
    "Γ": "gamma",
    "δ": "delta",
    "Δ": "delta",
    "κ": "kappa",
    "Κ": "kappa",
    "μ": "micro",
    "µ": "micro",
    "σ": "sigma",
    "ω": "omega",
}

# hyphen-like, space-like and underscore characters all become a single space
_SEPARATORS = re.compile(
    "["
    "-‐‑‒–—―−"  # hyphen/dash family
    "_\\s"
    "]+"
)
_STRIP_EDGE = re.compile(r"^[^0-9a-z]+|[^0-9a-z\)]+$")
_SLASH_PAD = re.compile(r"\s*/\s*")


def normalise(raw: str) -> str:
    """Normalise a string for exact matching.

    Lowercases, maps Greek letters to their Latin names, strips diacritics,
    turns every hyphen/dash/underscore/whitespace run into one space, tightens
    whitespace around slashes, and trims leading/trailing punctuation.

    Word content is never changed, so two strings only ever collide if they are
    the same words in the same order.
    """
    if raw is None:
        return ""
    s = str(raw)
    s = "".join(_GREEK.get(ch, ch) for ch in s)
    s = unicodedata.normalize("NFKD", s)
    s = "".join(ch for ch in s if not unicodedata.combining(ch))
    s = s.lower()
    s = _SEPARATORS.sub(" ", s)
    s = _SLASH_PAD.sub("/", s)
    s = s.strip()
    s = _STRIP_EDGE.sub("", s)
    return s.strip()


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #

_CACHE: Dict[str, Dict[str, dict]] = {}
_INDEX: Dict[str, Dict[str, str]] = {}


def _check_kind(kind: str) -> str:
    if kind not in KINDS:
        raise VocabError("unknown vocabulary kind %r; expected one of %s" % (kind, ", ".join(KINDS)))
    return kind


def _read_file(kind: str) -> Dict[str, dict]:
    path = os.path.join(vocab_dir(), "%s.json" % kind)
    if not os.path.isfile(path):
        raise VocabError("missing vocabulary file: %s" % path)
    with open(path, "r", encoding="utf-8") as fh:
        doc = json.load(fh)
    entries = doc.get("entries")
    if not isinstance(entries, dict):
        raise VocabError("%s: top-level 'entries' object is required" % path)
    out: Dict[str, dict] = {}
    for key, entry in entries.items():
        if not isinstance(entry, dict):
            raise VocabError("%s: entry %r is not an object" % (path, key))
        eid = entry.get("id")
        if not eid:
            raise VocabError("%s: entry %r has no 'id'" % (path, key))
        if eid != key:
            raise VocabError("%s: entry key %r != id %r" % (path, key, eid))
        if eid in out:
            raise VocabCollisionError("%s: duplicate id %r" % (path, eid))
        out[eid] = entry
    return out


def _build_index(kind: str, entries: Dict[str, dict]) -> Dict[str, str]:
    """Map normalised id/label/synonym -> entry id, raising on cross-entry clashes."""
    index: Dict[str, str] = {}
    for eid, entry in entries.items():
        keys: List[str] = [eid]
        label = entry.get("label")
        if label:
            keys.append(str(label))
        for syn in entry.get("synonyms") or []:
            keys.append(str(syn))
        # methods entries are also addressable by their actuator string
        act = entry.get("actuator")
        if act:
            keys.append(str(act))
        seen_here = set()
        for key in keys:
            nk = normalise(key)
            if not nk or nk in seen_here:
                continue
            seen_here.add(nk)
            prev = index.get(nk)
            if prev is not None and prev != eid:
                raise VocabCollisionError(
                    "%s: string %r maps to both %r and %r" % (kind, key, prev, eid)
                )
            index[nk] = eid
    return index


def load(kind: str) -> Dict[str, dict]:
    """Return ``{id: entry}`` for one vocabulary. Cached; result is shared, treat as read-only."""
    _check_kind(kind)
    if kind not in _CACHE:
        entries = _read_file(kind)
        _INDEX[kind] = _build_index(kind, entries)
        _CACHE[kind] = entries
    return _CACHE[kind]


def load_all() -> Dict[str, Dict[str, dict]]:
    """Load every vocabulary."""
    return {k: load(k) for k in KINDS}


def _index(kind: str) -> Dict[str, str]:
    load(kind)
    return _INDEX[kind]


# --------------------------------------------------------------------------- #
# resolving
# --------------------------------------------------------------------------- #

def resolve(kind: str, raw: str) -> Optional[str]:
    """Resolve ``raw`` to a vocabulary id, or ``None`` when there is no exact match.

    Exact and curated-synonym matching only, on the normalised string. There is
    no fuzzy fallback and there never should be: ``None`` is a legitimate,
    frequent and informative answer.
    """
    _check_kind(kind)
    if raw is None:
        return None
    nk = normalise(raw)
    if not nk:
        return None
    return _index(kind).get(nk)


def resolve_all(raw: str) -> Dict[str, Optional[str]]:
    """Try every vocabulary; useful for triage of unknown strings."""
    return {k: resolve(k, raw) for k in KINDS}


def report_unmatched(kind: str, strings: Iterable[str]) -> List[str]:
    """Distinct input strings that did not resolve, for vocabulary-extension reports.

    Order of first appearance is preserved; de-duplication is on the normalised
    form but the originally-written string is what comes back, because that is
    what a curator needs to read.
    """
    _check_kind(kind)
    seen = set()
    out: List[str] = []
    for raw in strings:
        nk = normalise(raw)
        if not nk or nk in seen:
            continue
        if resolve(kind, raw) is None:
            seen.add(nk)
            out.append(raw if isinstance(raw, str) else str(raw))
        else:
            seen.add(nk)
    return out


# --------------------------------------------------------------------------- #
# self-check
# --------------------------------------------------------------------------- #

def selfcheck(verbose: bool = False) -> Dict[str, dict]:
    """Load every vocabulary, enforce invariants, and return a summary.

    Invariants enforced (by raising):
      * entry key == entry id, ids unique within a file
      * no normalised string maps to two ids within a file
      * populations: 'nucleus' must name an id in regions.json
      * outcomes: tier in {A,B}; higher_is in {better,worse,None}
      * methods: canonical_valence in {activation,inhibition,ablation,unclear}
      * populations: curation_confidence in {high,medium,low}
    """
    summary: Dict[str, dict] = {}
    regions = load("regions")
    for kind in KINDS:
        entries = load(kind)
        summary[kind] = {
            "entries": len(entries),
            "index_keys": len(_index(kind)),
        }

    pops = load("populations")
    conf_counts = {"high": 0, "medium": 0, "low": 0}
    for eid, entry in pops.items():
        nuc = entry.get("nucleus")
        if nuc not in regions:
            raise VocabError(
                "populations/%s: nucleus %r is not an id in regions.json" % (eid, nuc)
            )
        conf = entry.get("curation_confidence")
        if conf not in conf_counts:
            raise VocabError(
                "populations/%s: curation_confidence %r not in high/medium/low" % (eid, conf)
            )
        conf_counts[conf] += 1
        nts = entry.get("neurotransmitter")
        if not isinstance(nts, list) or not nts:
            raise VocabError("populations/%s: neurotransmitter must be a non-empty list" % eid)
        allowed_nt = {
            "GABA", "glutamate", "dopamine", "serotonin", "histamine",
            "peptidergic", "mixed", "unknown",
        }
        bad = [n for n in nts if n not in allowed_nt]
        if bad:
            raise VocabError("populations/%s: neurotransmitter values %r not in SCHEMA enum" % (eid, bad))
    summary["populations"]["curation_confidence"] = conf_counts

    outs = load("outcomes")
    tiers = {"A": 0, "B": 0}
    higher_is_counts = {"better": 0, "worse": 0, "null": 0}
    for eid, entry in outs.items():
        tier = entry.get("tier")
        if tier not in tiers:
            raise VocabError("outcomes/%s: tier %r not in A/B" % (eid, tier))
        tiers[tier] += 1
        hi = entry.get("higher_is", "MISSING")
        if hi == "MISSING":
            raise VocabError("outcomes/%s: higher_is is required (use null deliberately)" % eid)
        if hi is None:
            higher_is_counts["null"] += 1
        elif hi in ("better", "worse"):
            higher_is_counts[hi] += 1
        else:
            raise VocabError("outcomes/%s: higher_is %r not in better/worse/null" % (eid, hi))
    summary["outcomes"]["tiers"] = tiers
    summary["outcomes"]["higher_is"] = higher_is_counts

    meths = load("methods")
    valences = {"activation": 0, "inhibition": 0, "ablation": 0, "unclear": 0}
    for eid, entry in meths.items():
        val = entry.get("canonical_valence")
        if val not in valences:
            raise VocabError("methods/%s: canonical_valence %r invalid" % (eid, val))
        valences[val] += 1
    summary["methods"]["canonical_valence"] = valences

    allen_known = sum(1 for e in regions.values() if e.get("allen_id") is not None)
    summary["regions"]["allen_id_set"] = allen_known
    summary["regions"]["allen_id_null"] = len(regions) - allen_known

    if verbose:
        print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


if __name__ == "__main__":  # pragma: no cover
    selfcheck(verbose=True)
