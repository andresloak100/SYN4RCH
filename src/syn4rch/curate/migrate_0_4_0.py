"""Migrate the curated dataset from schema v0.3.0 to v0.4.0 (D-060).

v0.4.0 adds one field, ``background_genotype``, as the fifth condition
dimension. That is the whole change to the data.

**This is the first migration in this project that needs no identity map, and
the reason is the point of D-060.** ``background_genotype`` is deliberately not
a component of :data:`~syn4rch.schema.record.RECORD_ID_FIELDS`, so adding it
re-identifies nothing: every ``record_id`` in the output is byte-identical to
its input. The three migrations before this one each rewrote identities and left
a map behind that every external reference must still be resolved through
(``record_id_map_0.1.0_to_0.2.0.json``, ``record_id_map_0.2.0_to_0.3.0.json``,
``record_id_map_d053.json``). This one leaves no such debt, and that was a
design criterion rather than a happy accident: the option that would have put
the genotype in the hash was scored and rejected, and the measured comparison is
in D-060.

So the migration asserts the property rather than assuming it. If any
``record_id`` moves, the run refuses and writes nothing — a silent
re-identification is exactly the failure an id map exists to paper over, and
this pass has no map to offer.

What it does
------------
1. Reads the input as raw JSON, because
   :func:`syn4rch.schema.io.read_jsonl` correctly refuses to load a v0.3.0 file
   against v0.4.0 code. Reading across versions is this module's job and
   nothing else's.
2. Leaves ``background_genotype`` null on every record. The field is *added*,
   not *filled*: all 508 records begin in the ``None`` state, which D-060
   defines as "the curator has not addressed genotype here". Filling it is a
   separate curation pass, for the same reason ``fill_conditions`` was separate
   from the v0.2.0 migration — a migration that also curated would make the two
   impossible to review apart.
3. Sets ``schema_version`` to :data:`TO_VERSION`.
4. Changes nothing else, and checks that claim field by field.

``FROM_VERSION`` and ``TO_VERSION`` are pinned literals, not
:data:`~syn4rch.schema.record.SCHEMA_VERSION`. D-046: a migration whose target
follows the current version stops being able to recompute the report it already
published.

Python 3.10 compatible, standard library only.
"""

import argparse
import datetime
import json
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from ..schema.record import (
    RECORD_ID_FIELDS,
    SCHEMA_VERSION,
    TierARecord,
)
from ..schema.validate import Severity, validate
from .migrate import _atomic_write_text, _records_text

__all__ = [
    "ADDED_FIELDS",
    "ACCEPTED_INPUT_VERSIONS",
    "FROM_VERSION",
    "TO_VERSION",
    "MigrationResult",
    "IdentityMovedError",
    "build_parser",
    "main",
    "migrate_records",
    "read_v0_3_0_records",
    "render_report_md",
]

#: Pinned, not :data:`SCHEMA_VERSION` (D-046).
FROM_VERSION = "0.3.0"
TO_VERSION = "0.4.0"

#: Inputs this pass will read. ``0.4.0`` is accepted so a second run is a
#: no-op rather than an error: a migration that cannot be run twice is one
#: an operator is afraid to re-run after an interruption.
ACCEPTED_INPUT_VERSIONS: Tuple[str, ...] = ("0.3.0", "0.4.0")

#: The only fields v0.4.0 adds. Named as data so the report and the test agree.
ADDED_FIELDS: Tuple[str, ...] = ("background_genotype",)


class IdentityMovedError(RuntimeError):
    """A ``record_id`` changed, which this migration has no map to explain."""


def read_v0_3_0_records(
    path: str,
) -> Tuple[List[TierARecord], List[Tuple[int, str]], Dict[str, int]]:
    """Read a curated JSONL at :data:`ACCEPTED_INPUT_VERSIONS`.

    Returns ``(records, errors, versions_seen)``, matching
    :func:`syn4rch.curate.migrate.read_records` and
    :func:`syn4rch.curate.fill_conditions.read_v0_2_0_records`.

    A record that already carries a non-null ``background_genotype`` is an
    **error**, not something to keep. The field does not exist before v0.4.0, so
    a value in it means the file has been through a pass this one does not know
    about, and migrating it would silently endorse that pass's work.
    """
    records: List[TierARecord] = []
    errors: List[Tuple[int, str]] = []
    versions: Dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                data = json.loads(line)
            except ValueError as exc:
                errors.append((number, "not JSON: %s" % exc))
                continue
            if not isinstance(data, dict):
                errors.append((number, "not a JSON object"))
                continue
            declared = str(data.get("schema_version"))
            if declared not in ACCEPTED_INPUT_VERSIONS:
                errors.append(
                    (
                        number,
                        "schema_version=%r is not one of %s"
                        % (declared, ", ".join(ACCEPTED_INPUT_VERSIONS)),
                    )
                )
                continue
            premature = [
                "%s=%r" % (name, data[name])
                for name in ADDED_FIELDS
                if declared == FROM_VERSION and data.get(name) is not None
            ]
            if premature:
                errors.append(
                    (
                        number,
                        "a v%s record already carries %s; these fields do not "
                        "exist before v%s, so this file has been through a pass "
                        "this migration does not know about"
                        % (FROM_VERSION, ", ".join(premature), TO_VERSION),
                    )
                )
                continue
            try:
                record = TierARecord.from_dict(data)
            except Exception as exc:
                errors.append((number, "%s: %s" % (type(exc).__name__, exc)))
                continue
            versions[declared] = versions.get(declared, 0) + 1
            records.append(record)
    return records, errors, versions


@dataclass
class MigrationResult:
    """The migrated records and everything a reader needs to check the pass."""

    records: List[TierARecord] = field(default_factory=list)
    versions_in: Dict[str, int] = field(default_factory=dict)
    n_already_migrated: int = 0
    fields_changed: Dict[str, int] = field(default_factory=dict)
    ids_unchanged: int = 0
    hard_findings_before: Dict[str, int] = field(default_factory=dict)
    hard_findings_after: Dict[str, int] = field(default_factory=dict)

    @property
    def n_records(self) -> int:
        return len(self.records)

    def to_dict(self) -> Dict[str, object]:
        return {
            "from_schema_version": FROM_VERSION,
            "to_schema_version": TO_VERSION,
            "added_fields": list(ADDED_FIELDS),
            "n_records": self.n_records,
            "versions_in": dict(sorted(self.versions_in.items())),
            "n_already_migrated": self.n_already_migrated,
            "fields_changed": dict(sorted(self.fields_changed.items())),
            "ids_unchanged": self.ids_unchanged,
            "ids_changed": self.n_records - self.ids_unchanged,
            "record_id_map_required": False,
            "hard_findings_before": dict(sorted(self.hard_findings_before.items())),
            "hard_findings_after": dict(sorted(self.hard_findings_after.items())),
        }


def _hard(records: Sequence[TierARecord]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for record in records:
        for finding in validate(record):
            if finding.severity is Severity.HARD:
                out[finding.code] = out.get(finding.code, 0) + 1
    return out


def migrate_records(records: Sequence[TierARecord]) -> MigrationResult:
    """Stamp v0.4.0 onto every record and prove nothing else moved.

    ``records`` is not mutated: the result holds fresh objects, so a caller can
    diff the before against the after.
    """
    result = MigrationResult()
    before = [TierARecord.from_dict(r.to_dict()) for r in records]
    result.hard_findings_before = _hard(before)

    for original in before:
        was = original.to_dict()
        if str(original.schema_version) == TO_VERSION:
            result.n_already_migrated += 1
        data = dict(was)
        data["schema_version"] = TO_VERSION
        for name in ADDED_FIELDS:
            data.setdefault(name, None)
        migrated = TierARecord.from_dict(data)

        # The whole safety claim of this migration, checked field by field
        # rather than asserted in a docstring.
        now = migrated.to_dict()
        for key in sorted(set(was) | set(now)):
            if was.get(key) != now.get(key):
                result.fields_changed[key] = result.fields_changed.get(key, 0) + 1

        if migrated.record_id != original.record_id:
            raise IdentityMovedError(
                "record_id moved from %r to %r. v%s adds %s, which is not a "
                "component of RECORD_ID_FIELDS (%s), so no identity may change "
                "and this pass has no id map to resolve one that did."
                % (
                    original.record_id,
                    migrated.record_id,
                    TO_VERSION,
                    ", ".join(ADDED_FIELDS),
                    ", ".join(RECORD_ID_FIELDS),
                )
            )
        result.ids_unchanged += 1
        result.records.append(migrated)

    result.hard_findings_after = _hard(result.records)
    return result


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #

def _table(header: Sequence[str], rows: Sequence[Sequence[object]]) -> List[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    for row in rows:
        out.append("| " + " | ".join(str(c) for c in row) + " |")
    return out


def render_report_md(
    result: MigrationResult, *, generated: str, records_path: str
) -> str:
    data = result.to_dict()
    lines: List[str] = []
    add = lines.append
    add("# Schema migration v%s -> v%s" % (FROM_VERSION, TO_VERSION))
    add("")
    add(
        "v%s adds one field, `background_genotype`, as the fifth condition "
        "dimension (D-060). That is the whole change to the data: the field is "
        "**added, not filled**. Every record begins in the `None` state, which "
        "D-060 defines as \"the curator has not addressed genotype here\" and "
        "which `condition_key` pools with `not_reported` into `unknown`. "
        "Filling it is a separate curation pass." % TO_VERSION
    )
    add("")
    add("## Run")
    add("")
    lines.extend(
        _table(
            ["field", "value"],
            [
                ["generated (UTC)", generated],
                ["input", "`%s`" % records_path],
                ["from / to", "%s -> %s" % (FROM_VERSION, TO_VERSION)],
                ["records", data["n_records"]],
                ["versions in", ", ".join("%s x%d" % (k, v) for k, v in data["versions_in"].items())],
                ["already at %s (no-op)" % TO_VERSION, data["n_already_migrated"]],
                ["**record ids unchanged**", "**%s of %s**" % (data["ids_unchanged"], data["n_records"])],
                ["record ids changed", data["ids_changed"]],
                ["**identity map required**", "**no**"],
            ],
        )
    )
    add("")
    add("## No identity map, and why that is the finding")
    add("")
    add(
        "This is the first migration in this project that leaves no identity "
        "map. `background_genotype` is deliberately **not** a component of "
        "`RECORD_ID_FIELDS`, so adding it re-identifies nothing — every "
        "`record_id` below is byte-identical to its input, checked per record "
        "rather than assumed, and the pass refuses outright if one moves."
    )
    add("")
    add(
        "The three passes before this one each rewrote identities and left a "
        "map that every external reference must still be resolved through: "
        "`record_id_map_0.1.0_to_0.2.0.json`, "
        "`record_id_map_0.2.0_to_0.3.0.json`, `record_id_map_d053.json`. That "
        "debt is permanent. Avoiding a fourth was a design criterion in D-060, "
        "not a happy accident: the variant that would have put the genotype in "
        "the hash was scored against this one and rejected, and the measured "
        "comparison is in that entry."
    )
    add("")
    add("## What changed, field by field")
    add("")
    changed = data["fields_changed"]
    if changed:
        lines.extend(
            _table(
                ["field", "records changed"],
                [[("`%s`" % k), v] for k, v in changed.items()],
            )
        )
        add("")
        add(
            "`schema_version` is the stamp, and on a clean run it is the only "
            "row. `background_genotype` does **not** appear, which is correct "
            "and worth stating: the dataclass declares it with a default of "
            "`None`, so reading a v%s record already materialises it as null "
            "and there is no value for the migration to change. The field is "
            "added to the *schema*, not written to the *data*. Any other field "
            "in this table would be a defect." % FROM_VERSION
        )
    else:
        add("Nothing changed: every input record was already at v%s." % TO_VERSION)
    add("")
    add("## Hard findings, before and after")
    add("")
    lines.extend(
        _table(
            ["when", "hard findings"],
            [
                ["before", data["hard_findings_before"] or "none"],
                ["after", data["hard_findings_after"] or "none"],
            ],
        )
    )
    add("")
    add(
        "A migration that introduces a hard finding has broken something. The "
        "two rows are printed rather than compared in prose so a reader can see "
        "which code it would have been."
    )
    add("")
    add("## What this pass does not do")
    add("")
    add(
        "* It does not curate. `background_genotype` is null on every record "
        "afterwards, and the 71 records D-060 identifies as carrying an "
        "establishable background genotype are curated by a separate pass."
    )
    add(
        "* It does not touch any other field, and the table above is the "
        "evidence rather than the claim."
    )
    add(
        "* It does not renumber anything. `FROM_VERSION` and `TO_VERSION` are "
        "pinned literals rather than `SCHEMA_VERSION`, so this report still "
        "recomputes after the schema moves again (D-046)."
    )
    add("")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _repo_root() -> str:
    here = os.path.abspath(__file__)
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(here))))


def _default(*parts: str) -> str:
    return os.path.join(_repo_root(), *parts)


def _now() -> str:
    return (
        datetime.datetime.now(datetime.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python3 -m syn4rch.curate.migrate_0_4_0",
        description=(
            "Migrate the curated dataset from schema v%s to v%s: add "
            "background_genotype as the fifth condition dimension, null on "
            "every record, and prove no record_id moved."
            % (FROM_VERSION, TO_VERSION)
        ),
    )
    parser.add_argument("records", nargs="?", help="curated JSONL to rewrite")
    parser.add_argument("--backup", default=None)
    parser.add_argument("--report", default=None)
    parser.add_argument("--json", dest="json_path", default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    records_path = args.records or _default("data", "curated", "records.jsonl")

    try:
        records, errors, versions = read_v0_3_0_records(records_path)
    except OSError as exc:
        sys.stderr.write("error: could not read %s: %s\n" % (records_path, exc))
        return 2
    if errors:
        sys.stderr.write(
            "error: %d line(s) could not be read; nothing written\n" % len(errors)
        )
        for number, why in errors[:10]:
            sys.stderr.write("  line %d: %s\n" % (number, why))
        return 1
    if not records:
        sys.stderr.write("error: %s holds no records\n" % records_path)
        return 2

    try:
        result = migrate_records(records)
    except IdentityMovedError as exc:
        sys.stderr.write("error: %s\n" % exc)
        return 1
    result.versions_in = versions

    introduced = {
        code: n
        for code, n in result.hard_findings_after.items()
        if n > result.hard_findings_before.get(code, 0)
    }
    if introduced:
        sys.stderr.write(
            "error: the migration would introduce hard finding(s) %s; nothing "
            "written\n" % introduced
        )
        return 1

    generated = _now()
    relative = os.path.relpath(os.path.abspath(records_path), _repo_root())
    report = render_report_md(result, generated=generated, records_path=relative)
    report_path = args.report or _default(
        "reports", "schema_migration_%s.md" % TO_VERSION
    )

    if not args.quiet:
        sys.stdout.write(
            "records %d | ids unchanged %d | fields changed %s | no id map needed\n"
            % (
                result.n_records,
                result.ids_unchanged,
                dict(sorted(result.fields_changed.items())),
            )
        )
    if args.dry_run:
        sys.stdout.write("dry run: nothing written\n")
        return 0

    backup_path = args.backup or (records_path + ".pre-%s.bak" % TO_VERSION)
    _atomic_write_text(backup_path, _records_text(records))
    _atomic_write_text(records_path, _records_text(result.records))
    _atomic_write_text(report_path, report)
    if args.json_path:
        _atomic_write_text(
            args.json_path, json.dumps(result.to_dict(), indent=2) + "\n"
        )
    if not args.quiet:
        sys.stdout.write(
            "wrote %s\nbackup %s\nreport %s\n"
            % (records_path, backup_path, report_path)
        )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
