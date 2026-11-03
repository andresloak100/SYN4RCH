"""JSONL persistence for ``data/curated/records.jsonl``.

One record per line, UTF-8, ``\\n``-terminated, keys in schema declaration
order. Round-tripping is exact: ``read_jsonl(write_jsonl(rs)) == rs``.

Version handling is strict by design. A file written under a different
``schema_version`` is **refused**, not coerced: loading v0.2 data with v0.1
field semantics is how a dataset acquires values nobody can account for.
Migration is an explicit, separate step.

Python 3.10 compatible, standard library only.
"""

import json
import os
from typing import Iterable, Iterator, List, Optional

from .record import SCHEMA_VERSION, SchemaCoercionError, TierARecord

__all__ = [
    "JsonlError",
    "SchemaVersionError",
    "dumps_record",
    "loads_record",
    "iter_jsonl",
    "read_jsonl",
    "write_jsonl",
    "append_jsonl",
]


class JsonlError(ValueError):
    """A line of a JSONL file could not be read as a Tier A record."""

    def __init__(self, message, path=None, line_number=None):
        self.path = path
        self.line_number = line_number
        where = ""
        if path is not None or line_number is not None:
            where = " (%s:%s)" % (path if path else "<stream>", line_number)
        super().__init__(message + where)


class SchemaVersionError(JsonlError):
    """A line declares a ``schema_version`` this code does not implement."""


def dumps_record(record: TierARecord) -> str:
    """Serialise one record to a single JSON line (no trailing newline)."""
    return json.dumps(record.to_dict(), ensure_ascii=False, sort_keys=False)


def loads_record(
    line: str,
    expect_version: Optional[str] = SCHEMA_VERSION,
    strict_fields: bool = True,
    path=None,
    line_number=None,
) -> TierARecord:
    """Parse one JSONL line.

    ``expect_version`` of ``None`` disables the version check (for migration
    tooling only, which must then handle the differences explicitly).
    """
    try:
        data = json.loads(line)
    except json.JSONDecodeError as exc:
        raise JsonlError("invalid JSON: %s" % exc.msg, path, line_number)
    if not isinstance(data, dict):
        raise JsonlError(
            "expected a JSON object, got %s" % type(data).__name__, path, line_number
        )
    found = data.get("schema_version")
    if expect_version is not None and found != expect_version:
        raise SchemaVersionError(
            "refusing to load schema_version=%r; this code implements %r. "
            "Migrate explicitly rather than loading across versions."
            % (found, expect_version),
            path,
            line_number,
        )
    try:
        return TierARecord.from_dict(data, strict=strict_fields)
    except SchemaCoercionError as exc:
        raise JsonlError(str(exc), path, line_number)


def iter_jsonl(
    path,
    expect_version: Optional[str] = SCHEMA_VERSION,
    strict_fields: bool = True,
) -> Iterator[TierARecord]:
    """Yield records from a JSONL file. Blank lines are skipped.

    Raises :class:`JsonlError` / :class:`SchemaVersionError` on the first bad
    line, naming the line number.
    """
    with open(path, "r", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            yield loads_record(
                line,
                expect_version=expect_version,
                strict_fields=strict_fields,
                path=path,
                line_number=number,
            )


def read_jsonl(
    path,
    expect_version: Optional[str] = SCHEMA_VERSION,
    strict_fields: bool = True,
) -> List[TierARecord]:
    return list(iter_jsonl(path, expect_version=expect_version, strict_fields=strict_fields))


def _write(path, records: Iterable[TierARecord], mode: str) -> int:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    count = 0
    with open(path, mode, encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(dumps_record(record))
            handle.write("\n")
            count += 1
    return count


def write_jsonl(path, records: Iterable[TierARecord]) -> int:
    """Write (replacing) ``path``. Returns the number of records written."""
    return _write(path, records, "w")


def append_jsonl(path, records: Iterable[TierARecord]) -> int:
    """Append to ``path``, creating it if needed (``data/raw`` discipline)."""
    return _write(path, records, "a")
