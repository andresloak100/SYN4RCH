"""On-disk cache, polite HTTP client and append-only raw provenance log.

Three jobs:

1. :class:`Cache` -- a sqlite table keyed by (source, url) holding the response
   body, so a re-run costs no network.
2. :func:`fetch` -- urllib with a fixed User-Agent, a per-host rate limiter and
   retry-with-backoff on 429/5xx. **No email address is ever sent**: not in the
   URL, not in a header. NCBI asks for a tool/email pair for higher limits; we
   stay anonymous and inside the anonymous limit instead of inventing contact
   details.
3. The raw log -- every response actually fetched from the network is written
   byte-for-byte under ``data/raw/<source>/`` and a line is appended to
   ``data/raw/manifest.jsonl`` recording the UTC timestamp, the source, the exact
   query string that produced it, the full URL, the HTTP status, the sha256 and
   the file path. Append-only: nothing in ``data/raw`` is ever rewritten.

Standard library only; Python 3.10 compatible.
"""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import time
import urllib.error
import http.client
import urllib.request
from typing import Dict, Optional, Tuple

USER_AGENT = "SYN4RCH/0.1 (research tool)"

# NCBI's documented anonymous ceiling is 3 requests/second; the project brief
# caps us at 2/s, so we use 0.5 s minimum spacing for NCBI and 0.34 s elsewhere.
_MIN_INTERVAL: Dict[str, float] = {
    "eutils.ncbi.nlm.nih.gov": 0.5,
    "www.ebi.ac.uk": 0.34,
    "api.crossref.org": 0.34,
}
_DEFAULT_INTERVAL = 0.5

_last_call: Dict[str, float] = {}


def repo_root() -> str:
    here = os.path.abspath(os.path.dirname(__file__))
    cur = here
    for _ in range(8):
        if os.path.isdir(os.path.join(cur, "data")) and os.path.isdir(os.path.join(cur, "docs")):
            return cur
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    raise RuntimeError("could not locate repo root above %s" % here)


def data_path(*parts: str) -> str:
    return os.path.join(repo_root(), "data", *parts)


def _ensure_dir(path: str) -> None:
    if path and not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)


def _host_of(url: str) -> str:
    try:
        return url.split("://", 1)[1].split("/", 1)[0].lower()
    except IndexError:
        return url


def _throttle(url: str) -> None:
    host = _host_of(url)
    interval = _MIN_INTERVAL.get(host, _DEFAULT_INTERVAL)
    prev = _last_call.get(host)
    now = time.time()
    if prev is not None:
        wait = interval - (now - prev)
        if wait > 0:
            time.sleep(wait)
    _last_call[host] = time.time()


class RawLog:
    """Append-only provenance log under data/raw/."""

    def __init__(self, root: Optional[str] = None) -> None:
        self.root = root or data_path("raw")
        _ensure_dir(self.root)
        self.manifest = os.path.join(self.root, "manifest.jsonl")

    def record(self, source: str, query: str, url: str, status: int, body: bytes,
               ext: str = "txt") -> str:
        sha = hashlib.sha256(body).hexdigest()
        subdir = os.path.join(self.root, source)
        _ensure_dir(subdir)
        path = os.path.join(subdir, "%s.%s" % (sha, ext))
        if not os.path.exists(path):
            with open(path, "wb") as fh:
                fh.write(body)
        line = {
            "retrieved_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "source": source,
            "query": query,
            "url": url,
            "http_status": status,
            "bytes": len(body),
            "sha256": sha,
            "path": os.path.relpath(path, repo_root()),
        }
        with open(self.manifest, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(line, sort_keys=True) + "\n")
        return path


class Cache:
    """sqlite response cache keyed by (source, url)."""

    def __init__(self, db_path: Optional[str] = None) -> None:
        self.db_path = db_path or data_path("interim", "harvest_cache.sqlite")
        _ensure_dir(os.path.dirname(self.db_path))
        self.conn = sqlite3.connect(self.db_path)
        self.conn.execute(
            "CREATE TABLE IF NOT EXISTS responses ("
            " source TEXT NOT NULL,"
            " url TEXT NOT NULL,"
            " query TEXT,"
            " http_status INTEGER,"
            " body TEXT,"
            " sha256 TEXT,"
            " fetched_at_utc TEXT,"
            " PRIMARY KEY (source, url))"
        )
        self.conn.commit()
        self.hits = 0
        self.misses = 0

    def get(self, source: str, url: str) -> Optional[str]:
        cur = self.conn.execute(
            "SELECT body FROM responses WHERE source=? AND url=?", (source, url)
        )
        row = cur.fetchone()
        if row is None:
            return None
        return row[0]

    def put(self, source: str, url: str, query: str, status: int, body: str) -> None:
        self.conn.execute(
            "INSERT OR REPLACE INTO responses"
            " (source,url,query,http_status,body,sha256,fetched_at_utc)"
            " VALUES (?,?,?,?,?,?,?)",
            (
                source, url, query, status, body,
                hashlib.sha256(body.encode("utf-8")).hexdigest(),
                time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            ),
        )
        self.conn.commit()

    def stats(self) -> Dict[str, int]:
        cur = self.conn.execute("SELECT COUNT(*) FROM responses")
        return {"rows": cur.fetchone()[0], "hits": self.hits, "misses": self.misses}

    def close(self) -> None:
        self.conn.close()


class FetchError(Exception):
    def __init__(self, url: str, status: Optional[int], message: str) -> None:
        super().__init__("%s (status=%s) %s" % (url, status, message))
        self.url = url
        self.status = status


def fetch(url: str, *, source: str, query: str, cache: Optional[Cache] = None,
          raw_log: Optional[RawLog] = None, ext: str = "txt",
          max_retries: int = 4, timeout: int = 60,
          use_cache: bool = True, post_url: Optional[str] = None,
          post_data: Optional[bytes] = None) -> Tuple[str, bool]:
    """GET ``url``, returning ``(body_text, from_cache)``.

    Retries on 429 and 5xx with exponential backoff (1, 2, 4, 8 s) plus the
    server's Retry-After when it gives one. Never attaches an email address.

    When ``post_data`` is given the request is sent as a POST to ``post_url``
    (needed because the vocabulary-generated queries exceed the ~8 kB URL limit
    and NCBI answers 414 to a GET). ``url`` then stays the canonical, fully
    expanded request URL and remains the cache key and the logged provenance
    URL, so a POSTed query is still exactly reproducible from the raw log.
    """
    if cache is not None and use_cache:
        hit = cache.get(source, url)
        if hit is not None:
            cache.hits += 1
            return hit, True
        cache.misses += 1

    headers = {"User-Agent": USER_AGENT, "Accept": "*/*"}
    if post_data is not None:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
        req = urllib.request.Request(post_url or url.split("?", 1)[0],
                                     data=post_data, headers=headers, method="POST")
    else:
        req = urllib.request.Request(url, headers=headers)

    last_err: Optional[str] = None
    last_status: Optional[int] = None
    for attempt in range(max_retries + 1):
        _throttle(url)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                status = getattr(resp, "status", 200) or 200
                body = resp.read()
        except urllib.error.HTTPError as exc:  # noqa: PERF203
            status = exc.code
            try:
                body = exc.read()
            except Exception:
                body = b""
            if status in (429, 500, 502, 503, 504) and attempt < max_retries:
                retry_after = exc.headers.get("Retry-After") if exc.headers else None
                wait = 2.0 ** attempt
                if retry_after:
                    try:
                        wait = max(wait, float(retry_after))
                    except ValueError:
                        pass
                last_status, last_err = status, "retrying after HTTP %s" % status
                time.sleep(wait)
                continue
            raise FetchError(url, status, body[:400].decode("utf-8", "replace"))
        except urllib.error.URLError as exc:
            if attempt < max_retries:
                last_status, last_err = None, str(exc.reason)
                time.sleep(2.0 ** attempt)
                continue
            raise FetchError(url, None, str(exc.reason))
        except TimeoutError as exc:
            if attempt < max_retries:
                last_status, last_err = None, "timeout: %s" % exc
                time.sleep(2.0 ** attempt)
                continue
            raise FetchError(url, None, "timeout")
        except (http.client.HTTPException, ConnectionError, OSError) as exc:
            # A server that drops the connection mid-response (e.g.
            # http.client.RemoteDisconnected) is a transient network failure,
            # exactly like a URLError. It must be retried and, if it persists,
            # surfaced as a FetchError the caller already handles. Letting it
            # escape aborted a whole harvest run during full-text retrieval and
            # discarded the in-memory record store before candidates were written.
            if attempt < max_retries:
                last_status, last_err = None, "%s: %s" % (type(exc).__name__, exc)
                time.sleep(2.0 ** attempt)
                continue
            raise FetchError(url, None, "%s: %s" % (type(exc).__name__, exc))

        if raw_log is not None:
            raw_log.record(source, query, url, status, body, ext=ext)
        text = body.decode("utf-8", "replace")
        if cache is not None:
            cache.put(source, url, query, status, text)
        return text, False

    raise FetchError(url, last_status, last_err or "exhausted retries")
