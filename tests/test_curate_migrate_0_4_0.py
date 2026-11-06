"""The v0.3.0 -> v0.4.0 migration (D-060).

The property this migration exists to guarantee is that **no identity moves**,
because `background_genotype` is not a `RECORD_ID_FIELDS` component and this
pass therefore leaves no id map. Most of what follows tests that it refuses
rather than proceeding when that guarantee does not hold.
"""

import _path  # noqa: F401
import io
import json
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

import synthetic
from syn4rch.curate.migrate_0_4_0 import (
    ACCEPTED_INPUT_VERSIONS,
    ADDED_FIELDS,
    FROM_VERSION,
    TO_VERSION,
    IdentityMovedError,
    build_parser,
    main,
    migrate_records,
    read_v0_3_0_records,
    render_report_md,
)
from syn4rch.schema.record import RECORD_ID_FIELDS, SCHEMA_VERSION, TierARecord


def as_v0_3_0(record=None, **overrides) -> dict:
    """A record dict as a v0.3.0 file would hold it: old stamp, no new field."""
    data = (record or synthetic.good_record()).to_dict()
    data["schema_version"] = FROM_VERSION
    for name in ADDED_FIELDS:
        data.pop(name, None)
    data.update(overrides)
    return data


def write(path, dicts):
    with open(path, "w", encoding="utf-8") as handle:
        for d in dicts:
            handle.write(json.dumps(d, ensure_ascii=False) + "\n")


class VersionsTest(unittest.TestCase):
    def test_the_versions_are_pinned_literals_not_the_current_schema(self):
        """D-046: a target that follows SCHEMA_VERSION cannot recompute its own
        published report once the schema moves again."""
        self.assertEqual(FROM_VERSION, "0.3.0")
        self.assertEqual(TO_VERSION, "0.4.0")
        self.assertIn(TO_VERSION, ACCEPTED_INPUT_VERSIONS)
        self.assertIn(FROM_VERSION, ACCEPTED_INPUT_VERSIONS)

    def test_the_added_field_is_not_a_record_id_component(self):
        """The reason this migration needs no id map, asserted rather than
        trusted: if the field ever entered the hash, every test below is void."""
        for name in ADDED_FIELDS:
            self.assertNotIn(name, RECORD_ID_FIELDS)


class MigrateTest(unittest.TestCase):
    def test_only_the_stamp_changes(self):
        original = TierARecord.from_dict(as_v0_3_0())
        result = migrate_records([original])
        self.assertEqual(result.n_records, 1)
        self.assertEqual(result.fields_changed, {"schema_version": 1})
        self.assertEqual(str(result.records[0].schema_version), TO_VERSION)

    def test_the_added_field_is_null_not_filled(self):
        """A migration that also curated would make the two impossible to
        review apart. All records begin in the `None` state (D-060)."""
        result = migrate_records([TierARecord.from_dict(as_v0_3_0())])
        for name in ADDED_FIELDS:
            self.assertIsNone(getattr(result.records[0], name))

    def test_every_record_id_is_unchanged(self):
        rows = [
            TierARecord.from_dict(as_v0_3_0()),
            TierARecord.from_dict(
                as_v0_3_0(synthetic.other_lab_record())
            ),
            TierARecord.from_dict(as_v0_3_0(synthetic.molecular_record())),
        ]
        before = [r.record_id for r in rows]
        result = migrate_records(rows)
        self.assertEqual([r.record_id for r in result.records], before)
        self.assertEqual(result.ids_unchanged, len(rows))
        self.assertEqual(result.to_dict()["ids_changed"], 0)

    def test_no_identity_map_is_required(self):
        result = migrate_records([TierARecord.from_dict(as_v0_3_0())])
        self.assertIs(result.to_dict()["record_id_map_required"], False)

    def test_a_moved_identity_refuses_rather_than_proceeding(self):
        """There is no map to explain one, so the pass must stop.

        `RECORD_ID_FIELDS` cannot be monkeypatched to simulate this --
        `compute_record_id` has a fixed twelve-argument signature and raises
        `TypeError` on a thirteenth, which is itself a useful guard. So the
        identity is moved at the seam `migrate_records` actually compares: the
        record it rebuilds from the dict.
        """
        import syn4rch.curate.migrate_0_4_0 as M

        class Moved(TierARecord):
            @classmethod
            def from_dict(cls, data):
                record = TierARecord.from_dict(data)
                # Only the MIGRATED side, which is the only dict already
                # stamped to the target. migrate_records rebuilds the before
                # copy through this same call, so moving both would leave them
                # equal and prove nothing.
                if str(data.get("schema_version")) == TO_VERSION:
                    record.record_id = "0" * 16
                return record

        real = M.TierARecord
        try:
            M.TierARecord = Moved
            with self.assertRaises(IdentityMovedError) as caught:
                migrate_records([TierARecord.from_dict(as_v0_3_0())])
            self.assertIn("record_id moved", str(caught.exception))
            self.assertIn("no id map", str(caught.exception))
        finally:
            M.TierARecord = real

    def test_running_it_twice_is_a_no_op(self):
        once = migrate_records([TierARecord.from_dict(as_v0_3_0())])
        twice = migrate_records(once.records)
        self.assertEqual(twice.fields_changed, {})
        self.assertEqual(twice.n_already_migrated, 1)
        self.assertEqual(
            [r.record_id for r in twice.records],
            [r.record_id for r in once.records],
        )

    def test_it_resolves_the_version_mismatch_rather_than_introducing_anything(self):
        """The before-state legitimately carries a hard finding.

        A v0.3.0 record read under v0.4.0 code earns
        `S006_SCHEMA_VERSION_MISMATCH`, which is the whole reason the loader
        refuses the file. The migration's job is to *resolve* it, and the test
        asserts both halves: the finding is there before and gone after, and
        nothing new appears.
        """
        result = migrate_records([TierARecord.from_dict(as_v0_3_0())])
        self.assertEqual(
            result.hard_findings_before, {"S006_SCHEMA_VERSION_MISMATCH": 1}
        )
        self.assertEqual(result.hard_findings_after, {})

    def test_the_input_records_are_not_mutated(self):
        original = TierARecord.from_dict(as_v0_3_0())
        snapshot = original.to_dict()
        migrate_records([original])
        self.assertEqual(original.to_dict(), snapshot)


class ReaderTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.path = os.path.join(self.tmp, "records.jsonl")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_it_reads_across_the_version_the_loader_refuses(self):
        """`read_jsonl` correctly refuses a v0.3.0 file against v0.4.0 code.
        Reading across versions is this module's job and nothing else's."""
        from syn4rch.schema.io import SchemaVersionError, read_jsonl

        write(self.path, [as_v0_3_0()])
        with self.assertRaises(SchemaVersionError):
            list(read_jsonl(self.path))
        records, errors, versions = read_v0_3_0_records(self.path)
        self.assertEqual(errors, [])
        self.assertEqual(len(records), 1)
        self.assertEqual(versions, {FROM_VERSION: 1})

    def test_a_v0_4_0_input_is_accepted_so_a_rerun_is_safe(self):
        write(self.path, [synthetic.good_record().to_dict()])
        records, errors, versions = read_v0_3_0_records(self.path)
        self.assertEqual(errors, [])
        self.assertEqual(versions, {SCHEMA_VERSION: 1})

    def test_an_unaccepted_version_is_an_error_not_a_coercion(self):
        write(self.path, [as_v0_3_0(schema_version="0.1.0")])
        records, errors, _ = read_v0_3_0_records(self.path)
        self.assertEqual(records, [])
        self.assertEqual(len(errors), 1)
        self.assertIn("not one of", errors[0][1])

    def test_a_premature_value_in_the_new_field_is_refused(self):
        """The field does not exist before v0.4.0, so a value in it means the
        file went through a pass this migration does not know about."""
        data = as_v0_3_0()
        data["background_genotype"] = "npy:loss_of_function:global"
        write(self.path, [data])
        records, errors, _ = read_v0_3_0_records(self.path)
        self.assertEqual(records, [])
        self.assertEqual(len(errors), 1)
        self.assertIn("does not know about", errors[0][1])

    def test_bad_json_is_reported_with_its_line_number(self):
        with open(self.path, "w", encoding="utf-8") as handle:
            handle.write("{not json\n")
        records, errors, _ = read_v0_3_0_records(self.path)
        self.assertEqual(records, [])
        self.assertEqual(errors[0][0], 1)


class ReportTest(unittest.TestCase):
    def test_the_report_states_that_no_id_map_is_needed_and_why(self):
        result = migrate_records([TierARecord.from_dict(as_v0_3_0())])
        text = render_report_md(
            result, generated="2026-01-01T00:00:00Z", records_path="x.jsonl"
        )
        self.assertIn("No identity map", text)
        self.assertIn("RECORD_ID_FIELDS", text)
        self.assertIn("record_id_map_d053.json", text)
        self.assertIn("added, not filled", text)

    def test_the_report_explains_why_the_new_field_is_absent_from_the_diff(self):
        result = migrate_records([TierARecord.from_dict(as_v0_3_0())])
        text = render_report_md(
            result, generated="2026-01-01T00:00:00Z", records_path="x.jsonl"
        )
        self.assertIn("does **not** appear", text)
        self.assertIn("added to the *schema*, not written to the *data*", text)


class CLITest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.path = os.path.join(self.tmp, "records.jsonl")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = main(argv)
        return code, out.getvalue(), err.getvalue()

    def test_dry_run_writes_nothing(self):
        write(self.path, [as_v0_3_0()])
        unchanged = open(self.path, encoding="utf-8").read()
        code, out, err = self._run([self.path, "--dry-run"])
        self.assertEqual(code, 0, err)
        self.assertIn("dry run", out)
        self.assertIn("no id map needed", out)
        self.assertEqual(open(self.path, encoding="utf-8").read(), unchanged)

    def test_it_writes_the_dataset_a_backup_and_a_report(self):
        write(self.path, [as_v0_3_0()])
        report = os.path.join(self.tmp, "migration.md")
        js = os.path.join(self.tmp, "migration.json")
        code, _, err = self._run(
            [self.path, "--report", report, "--json", js, "--quiet"]
        )
        self.assertEqual(code, 0, err)
        rows = [json.loads(l) for l in open(self.path, encoding="utf-8") if l.strip()]
        self.assertEqual([r["schema_version"] for r in rows], [TO_VERSION])
        self.assertTrue(os.path.exists(self.path + ".pre-%s.bak" % TO_VERSION))
        self.assertIn("No identity map", open(report, encoding="utf-8").read())
        self.assertIs(json.load(open(js))["record_id_map_required"], False)

    def test_a_backup_preserves_the_pre_migration_file(self):
        write(self.path, [as_v0_3_0()])
        before = open(self.path, encoding="utf-8").read()
        self._run([self.path, "--quiet"])
        backup = open(self.path + ".pre-%s.bak" % TO_VERSION, encoding="utf-8").read()
        self.assertIn(FROM_VERSION, backup)
        self.assertEqual(
            json.loads(before)["record_id"], json.loads(backup)["record_id"]
        )

    def test_an_unreadable_line_blocks_the_whole_write(self):
        with open(self.path, "w", encoding="utf-8") as handle:
            handle.write(json.dumps(as_v0_3_0()) + "\n{not json\n")
        unchanged = open(self.path, encoding="utf-8").read()
        code, _, err = self._run([self.path])
        self.assertEqual(code, 1)
        self.assertIn("could not be read", err)
        self.assertEqual(open(self.path, encoding="utf-8").read(), unchanged)

    def test_an_empty_file_is_refused(self):
        open(self.path, "w", encoding="utf-8").close()
        code, _, err = self._run([self.path])
        self.assertEqual(code, 2)
        self.assertIn("no records", err)

    def test_the_parser_defaults_to_the_curated_dataset(self):
        args = build_parser().parse_args([])
        self.assertIsNone(args.records)
        self.assertFalse(args.dry_run)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
