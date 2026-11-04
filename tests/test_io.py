"""JSONL round-tripping and version refusal."""

import _path  # noqa: F401
import json
import os
import shutil
import tempfile
import unittest

import synthetic
from syn4rch.schema import (
    SCHEMA_VERSION,
    JsonlError,
    SchemaVersionError,
    TierARecord,
    dumps_record,
    loads_record,
    read_jsonl,
    write_jsonl,
)
from syn4rch.schema.record import SchemaCoercionError


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="syn4rch-test-")
        self.path = os.path.join(self.tmp, "records.jsonl")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)


class TestRoundTrip(TempDirCase):
    def test_write_then_read_is_lossless(self):
        original = [
            synthetic.good_record(),
            synthetic.other_lab_record(),
            synthetic.good_record(
                direction="not_reported",
                effect_size_value=None,
                effect_size_type=None,
                effect_size_units=None,
                source_locator="Table S1 (synthetic)",
                confidence="moderate",
            ),
        ]
        self.assertEqual(write_jsonl(self.path, original), 3)
        loaded = read_jsonl(self.path)
        self.assertEqual(len(loaded), 3)
        for before, after in zip(original, loaded):
            self.assertEqual(before, after)
            self.assertEqual(before.to_dict(), after.to_dict())
            self.assertEqual(before.record_id, after.record_id)

    def test_enums_and_lists_survive(self):
        record = synthetic.good_record(
            neurotransmitter=["GABA", "peptidergic"],
            marker_genes=["Synth1", "Synth2"],
            control_types=["vehicle", "ligand_in_nonexpressing"],
            conflicts_with=["aaaaaaaaaaaaaaaa"],
        )
        back = loads_record(dumps_record(record))
        self.assertEqual([n.value for n in back.neurotransmitter], ["GABA", "peptidergic"])
        self.assertEqual(back.marker_genes, ["Synth1", "Synth2"])
        self.assertEqual([c.value for c in back.control_types],
                         ["vehicle", "ligand_in_nonexpressing"])
        self.assertEqual(back.conflicts_with, ["aaaaaaaaaaaaaaaa"])

    def test_one_json_object_per_line(self):
        write_jsonl(self.path, [synthetic.good_record(), synthetic.other_lab_record()])
        with open(self.path, "r", encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        self.assertEqual(len(lines), 2)
        for line in lines:
            self.assertIsInstance(json.loads(line), dict)

    def test_blank_lines_are_skipped(self):
        with open(self.path, "w", encoding="utf-8") as handle:
            handle.write(dumps_record(synthetic.good_record()) + "\n")
            handle.write("\n   \n")
            handle.write(dumps_record(synthetic.other_lab_record()) + "\n")
        self.assertEqual(len(read_jsonl(self.path)), 2)

    def test_non_ascii_is_preserved(self):
        record = synthetic.good_record(
            last_author_affiliation="Institut für synthetische Daten, Überall"
        )
        back = loads_record(dumps_record(record))
        self.assertEqual(back.last_author_affiliation, record.last_author_affiliation)


class TestVersionRefusal(TempDirCase):
    def _line_with_version(self, version):
        data = synthetic.good_record().to_dict()
        data["schema_version"] = version
        return json.dumps(data)

    #: A version that is deliberately not the current one. Spelled against
    #: SCHEMA_VERSION rather than as a literal: the literal that used to sit
    #: here became the *correct* version at the v0.2.0 bump, and the test
    #: silently stopped checking anything.
    FOREIGN = SCHEMA_VERSION + "-not-this"

    def test_foreign_version_is_refused(self):
        for version in (self.FOREIGN, "1.0.0", "0.1", None, ""):
            with self.assertRaises(SchemaVersionError) as ctx:
                loads_record(self._line_with_version(version))
            self.assertIn("refusing to load", str(ctx.exception))

    def test_refusal_names_the_line_number(self):
        with open(self.path, "w", encoding="utf-8") as handle:
            handle.write(dumps_record(synthetic.good_record()) + "\n")
            handle.write(self._line_with_version(self.FOREIGN) + "\n")
        with self.assertRaises(SchemaVersionError) as ctx:
            read_jsonl(self.path)
        self.assertIn(":2", str(ctx.exception))

    def test_current_version_is_accepted(self):
        record = loads_record(self._line_with_version(SCHEMA_VERSION))
        self.assertEqual(record.schema_version, SCHEMA_VERSION)

    def test_version_check_can_be_disabled_explicitly(self):
        record = loads_record(
            self._line_with_version(self.FOREIGN), expect_version=None
        )
        self.assertEqual(record.schema_version, self.FOREIGN)


class TestMalformedInput(TempDirCase):
    def test_invalid_json_names_the_line(self):
        with open(self.path, "w", encoding="utf-8") as handle:
            handle.write(dumps_record(synthetic.good_record()) + "\n")
            handle.write("{not json\n")
        with self.assertRaises(JsonlError) as ctx:
            read_jsonl(self.path)
        self.assertIn(":2", str(ctx.exception))

    def test_json_array_is_rejected(self):
        with self.assertRaises(JsonlError):
            loads_record("[1, 2, 3]")

    def test_unknown_field_is_refused_by_default(self):
        data = synthetic.good_record().to_dict()
        data["rejuvenation_score"] = 0.9
        with self.assertRaises(JsonlError):
            loads_record(json.dumps(data))
        record = loads_record(json.dumps(data), strict_fields=False)
        self.assertFalse(hasattr(record, "rejuvenation_score"))

    def test_unknown_enum_value_is_refused(self):
        data = synthetic.good_record().to_dict()
        data["direction"] = "improved"
        with self.assertRaises(JsonlError) as ctx:
            loads_record(json.dumps(data))
        self.assertIn("direction", str(ctx.exception))

    def test_bare_string_for_a_list_field_is_refused(self):
        with self.assertRaises(SchemaCoercionError):
            TierARecord(study_id="synthetic:test-001", marker_genes="Synth1")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
