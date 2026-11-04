"""Schema migration v0.1.0 -> v0.2.0: what it fills, and what it must not.

Every record here is fabricated (``tests/synthetic.py``) and every path is a
temp path. Nothing in this module reads or writes ``data/`` or ``reports/``,
with one deliberate exception: one test loads the committed
``data/vocab/targets.json`` to check that the real vocabulary can be indexed at
all. Everything else resolves against a fixture vocabulary, so a later edit to
the real one cannot quietly change what these tests prove.
"""

import _path  # noqa: F401
import io
import json
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout, redirect_stderr

import synthetic
from syn4rch.curate.migrate import (
    AMBIGUOUS,
    CONDITION_FIELDS_THE_MIGRATION_MAY_NOT_FILL,
    FROM_VERSION,
    NO_MATCH,
    NOT_ELIGIBLE,
    RECORD_ID_FIELDS_0_1_0,
    RESOLVED,
    TO_VERSION,
    build_synonym_index,
    compute_record_id_0_1_0,
    dataset_stats,
    load_targets,
    main,
    migrate_record,
    migrate_records,
    occurrences,
    remap_conflicts,
    render_report_md,
    resolve_target,
)
from syn4rch.schema.record import RECORD_ID_FIELDS, TierARecord
from syn4rch.schema.validate import Severity, validate

#: A fabricated targets vocabulary. Two real-looking gene symbols are enough to
#: exercise resolution, ambiguity and the no-match path; ``fakegene`` exists so
#: that no assertion here depends on the committed vocabulary's contents.
FIXTURE_TARGETS = {
    "trp53": {
        "id": "trp53",
        "label": "Synthetic p53 entry",
        "synonyms": ["Trp53", "p53"],
        "notes": "Fixture entry; the official symbol choice is tested elsewhere.",
    },
    "rptor": {
        "id": "rptor",
        "label": "Synthetic Raptor entry",
        "synonyms": ["Rptor", "Raptor"],
        "notes": "Fixture entry.",
    },
    "fakegene": {
        "id": "fakegene",
        "label": "Fabricated gene that is in no paper",
        "synonyms": ["Fakegene"],
        "notes": "Fixture entry; exists only so these tests cannot pass by "
        "accident against the real vocabulary.",
    },
}


def fixture_index():
    return build_synonym_index(FIXTURE_TARGETS)


def v010_record(**overrides) -> TierARecord:
    """A fabricated record as schema v0.1.0 would have stored it.

    The v0.2.0 fields are cleared and ``record_id`` is the seven-component
    v0.1.0 hash, so this is what the migration actually meets on disk.
    """
    record = synthetic.good_record(**overrides)
    data = record.to_dict()
    data["schema_version"] = FROM_VERSION
    for name in (
        "perturbation_class",
        "perturbation_polarity",
        "target_raw",
        "target_id",
    ) + CONDITION_FIELDS_THE_MIGRATION_MAY_NOT_FILL:
        data[name] = None
    data["record_id"] = compute_record_id_0_1_0(data)
    return TierARecord.from_dict(data)


def read_text(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()


def read_bytes(path):
    with open(path, "rb") as handle:
        return handle.read()


def read_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def write_jsonl_text(path, records):
    from syn4rch.schema.io import dumps_record

    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(dumps_record(record) + "\n")


class DerivedFieldsTest(unittest.TestCase):
    """What the migration is allowed to compute, it computes by rule."""

    def test_molecular_knockout_gets_class_and_polarity(self):
        record = v010_record(
            method="conditional_knockout",
            valence="unclear",
            actuator=None,
            targeting_strategy="synthetic floxed Trp53 x synthetic-Cre",
        )
        self.assertIsNone(record.perturbation_class)
        migrated, _ = migrate_record(record, fixture_index())
        self.assertEqual(migrated.perturbation_class.value, "molecular")
        self.assertEqual(migrated.perturbation_polarity.value, "decrease")
        self.assertEqual(migrated.schema_version, TO_VERSION)

    def test_overexpression_is_a_gain_of_function(self):
        migrated, _ = migrate_record(
            v010_record(method="conditional_overexpression", valence="unclear"),
            fixture_index(),
        )
        self.assertEqual(migrated.perturbation_class.value, "molecular")
        self.assertEqual(migrated.perturbation_polarity.value, "increase")

    def test_circuit_polarity_restates_valence(self):
        for valence, polarity in (
            ("activation", "increase"),
            ("inhibition", "decrease"),
            ("unclear", "unclear"),
        ):
            migrated, _ = migrate_record(
                v010_record(method="chemogenetic", valence=valence), fixture_index()
            )
            self.assertEqual(migrated.perturbation_class.value, "circuit_activity")
            self.assertEqual(migrated.perturbation_polarity.value, polarity)

    def test_structural_is_remove(self):
        migrated, _ = migrate_record(
            v010_record(method="genetic_ablation", valence="ablation"), fixture_index()
        )
        self.assertEqual(migrated.perturbation_class.value, "structural")
        self.assertEqual(migrated.perturbation_polarity.value, "remove")

    def test_migrated_records_carry_no_new_hard_finding(self):
        record = v010_record(method="conditional_knockout", valence="unclear")
        before = {f.code for f in validate(record) if f.severity is Severity.HARD}
        self.assertEqual(
            before, {"S005_RECORD_ID_INCONSISTENT", "S006_SCHEMA_VERSION_MISMATCH"}
        )
        migrated, _ = migrate_record(record, fixture_index())
        self.assertEqual(
            [f.code for f in validate(migrated) if f.severity is Severity.HARD], []
        )


class ConditionFieldsTest(unittest.TestCase):
    """The four condition fields are never invented. This is the whole point."""

    def test_condition_fields_stay_null(self):
        record = v010_record()
        migrated, _ = migrate_record(record, fixture_index())
        for name in CONDITION_FIELDS_THE_MIGRATION_MAY_NOT_FILL:
            self.assertIsNone(
                getattr(migrated, name),
                "%s was filled by the migration; it requires reading the paper"
                % name,
            )

    def test_not_reported_is_not_written_either(self):
        # ``not_reported`` is a claim that a curator checked the paper
        # (SCHEMA.md 3.3b). The migration checked nothing, so it may not say so.
        migrated, _ = migrate_record(v010_record(), fixture_index())
        for name in CONDITION_FIELDS_THE_MIGRATION_MAY_NOT_FILL:
            value = getattr(migrated, name)
            self.assertNotEqual(getattr(value, "value", value), "not_reported")

    def test_a_curated_condition_value_is_preserved(self):
        record = v010_record()
        data = record.to_dict()
        data["challenge"] = "fasted"
        data["challenge_detail"] = "16 h synthetic fast"
        data["condition_source"] = "methods"
        data["condition_locator"] = "Methods para 2 (synthetic)"
        data["record_id"] = compute_record_id_0_1_0(data)
        migrated, _ = migrate_record(TierARecord.from_dict(data), fixture_index())
        self.assertEqual(migrated.challenge.value, "fasted")
        self.assertEqual(migrated.condition_source.value, "methods")
        self.assertEqual(migrated.condition_locator, "Methods para 2 (synthetic)")

    def test_curated_subject_fields_are_untouched(self):
        record = v010_record(sex="female", diet="high_fat")
        migrated, _ = migrate_record(record, fixture_index())
        self.assertEqual(migrated.sex.value, "female")
        self.assertEqual(migrated.diet.value, "high_fat")
        self.assertEqual(migrated.age_group.value, record.age_group.value)


class TargetResolutionTest(unittest.TestCase):
    def test_one_gene_in_targeting_strategy_resolves(self):
        record = v010_record(
            method="conditional_knockout",
            valence="unclear",
            actuator=None,
            targeting_strategy="SynPOPp53 KO (p53 deleted in synthetic neurons)",
        )
        resolution = resolve_target(record, fixture_index())
        self.assertEqual(resolution.status, RESOLVED)
        self.assertEqual(resolution.target_id, "trp53")
        # The raw string is what the text said, not the vocabulary's spelling.
        self.assertEqual(resolution.target_raw, "p53")

    def test_the_gene_may_come_from_actuator(self):
        record = v010_record(
            method="conditional_overexpression",
            valence="unclear",
            actuator="Raptor",
            targeting_strategy="synthetic AAV into a synthetic nucleus",
        )
        migrated, resolution = migrate_record(record, fixture_index())
        self.assertEqual(resolution.status, RESOLVED)
        self.assertEqual(migrated.target_id, "rptor")
        self.assertEqual(migrated.target_raw, "Raptor")

    def test_two_genes_resolve_to_null_and_are_reported(self):
        record = v010_record(
            method="conditional_knockout",
            valence="unclear",
            actuator=None,
            targeting_strategy="synthetic double mutant: p53 and Raptor both floxed",
        )
        resolution = resolve_target(record, fixture_index())
        self.assertEqual(resolution.status, AMBIGUOUS)
        self.assertEqual(resolution.candidates, ("rptor", "trp53"))
        migrated, _ = migrate_record(record, fixture_index())
        self.assertIsNone(migrated.target_raw)
        self.assertIsNone(migrated.target_id)

        result = migrate_records([record], fixture_index())
        self.assertEqual(result.counts["target_ambiguous"], 1)
        self.assertEqual(result.counts["target_resolved"], 0)
        self.assertEqual(len(result.ambiguous), 1)
        self.assertEqual(result.ambiguous[0][0], result.records[0].record_id)
        self.assertTrue(result.needs_a_human)
        report = render_report_md(
            dataset_stats(result.records, fixture_index()),
            counts=result.counts,
            result=result,
            records_path="records.jsonl",
            id_map_path="map.json",
            backup_path=None,
            index=fixture_index(),
            generated="TEST",
        )
        self.assertIn("Ambiguous matches", report)
        self.assertIn(result.records[0].record_id, report)
        self.assertIn("`trp53`", report)

    def test_several_synonyms_of_one_gene_are_not_an_ambiguity(self):
        record = v010_record(
            method="conditional_knockout",
            valence="unclear",
            actuator="p53",
            targeting_strategy="floxed Trp53 (p53) in synthetic neurons",
        )
        resolution = resolve_target(record, fixture_index())
        self.assertEqual(resolution.status, RESOLVED)
        self.assertEqual(resolution.target_id, "trp53")
        # Longest match wins, so target_raw is the more specific spelling.
        self.assertEqual(resolution.target_raw, "Trp53")

    def test_no_match_leaves_it_null(self):
        record = v010_record(
            method="local_pharmacology",
            valence="unclear",
            actuator=None,
            targeting_strategy="synthetic drug infused into a synthetic nucleus",
        )
        resolution = resolve_target(record, fixture_index())
        self.assertEqual(resolution.status, NO_MATCH)
        migrated, _ = migrate_record(record, fixture_index())
        self.assertIsNone(migrated.target_id)

    def test_a_circuit_perturbation_may_not_carry_a_target(self):
        # A chemogenetic activation has no molecular target, even when the
        # curated text happens to name a gene.
        for method in ("chemogenetic", "optogenetic", "genetic_ablation", "lesion"):
            record = v010_record(
                method=method,
                valence="activation" if "genetic" not in method else "ablation",
                actuator=None,
                targeting_strategy="synthetic p53-reporter line, Raptor-naive",
            )
            resolution = resolve_target(record, fixture_index())
            self.assertEqual(resolution.status, NOT_ELIGIBLE, method)
            migrated, _ = migrate_record(record, fixture_index())
            self.assertIsNone(migrated.target_id, method)
            self.assertIsNone(migrated.target_raw, method)
            self.assertNotIn(
                "S008_TARGET_ON_NON_MOLECULAR",
                [f.code for f in validate(migrated)],
                method,
            )

    def test_an_already_curated_target_is_left_alone(self):
        record = v010_record(method="conditional_knockout", valence="unclear")
        data = record.to_dict()
        data["target_raw"] = "a curator's own spelling"
        data["target_id"] = "rptor"
        data["targeting_strategy"] = "text that names p53 instead"
        data["record_id"] = compute_record_id_0_1_0(data)
        migrated, resolution = migrate_record(
            TierARecord.from_dict(data), fixture_index()
        )
        self.assertEqual(resolution.status, "already_set")
        self.assertEqual(migrated.target_raw, "a curator's own spelling")
        self.assertEqual(migrated.target_id, "rptor")

    def test_matching_has_no_word_boundary_but_guards_digits(self):
        # Genotype strings glue the symbol to the driver and the allele.
        self.assertTrue(occurrences("AgRPp53 KO; WT controls", "p53"))
        self.assertTrue(occurrences("Sh2b1f/f;POMC-Cre", "Sh2b1"))
        self.assertTrue(occurrences("Htr2cfl/Y", "Htr2c"))
        self.assertTrue(occurrences("TrspAgrpKO", "Trsp"))
        # ...but a digit on either side means a different symbol.
        self.assertEqual(occurrences("Raf12 knockout", "Raf1"), [])
        self.assertEqual(occurrences("floxed 1p53x", "p53"), [])
        # Offsets index the original string, so target_raw can be sliced out.
        text = "bilateral AAV9-CAG-SH2B1 into the synthetic nucleus"
        start, end = occurrences(text, "Sh2b1")[0]
        self.assertEqual(text[start:end], "SH2B1")

    def test_dashes_and_spaces_fold_before_matching(self):
        self.assertTrue(occurrences("AgRP–Raptor mice", "Raptor"))
        self.assertTrue(occurrences("PKC epsilon floxed", "PKC epsilon"))


class SynonymIndexTest(unittest.TestCase):
    def test_cross_entry_substring_is_refused(self):
        with self.assertRaises(ValueError):
            build_synonym_index(
                {
                    "aaa": {"id": "aaa", "synonyms": ["Sirt1"]},
                    "bbb": {"id": "bbb", "synonyms": ["Sirt1a"]},
                }
            )

    def test_a_two_character_synonym_is_refused(self):
        with self.assertRaises(ValueError):
            build_synonym_index({"xy": {"id": "xy", "synonyms": ["xy"]}})

    def test_an_entry_with_no_synonyms_is_refused(self):
        with self.assertRaises(ValueError):
            build_synonym_index({"abc": {"id": "abc", "synonyms": []}})

    def test_the_entry_id_is_itself_matchable(self):
        index = fixture_index()
        self.assertIn(("rptor", "rptor"), index.pairs)

    def test_the_committed_vocabulary_indexes(self):
        # The one test here that reads data/: the real file must be usable.
        index = build_synonym_index(load_targets())
        self.assertGreaterEqual(len(index.ids), 1)
        for target_id, synonym in index.pairs:
            self.assertTrue(synonym.strip())
            self.assertEqual(target_id, target_id.lower())


class IdMapTest(unittest.TestCase):
    def setUp(self):
        self.records = [
            v010_record(
                method="conditional_knockout",
                valence="unclear",
                actuator=None,
                targeting_strategy="synthetic p53 floxed line",
                source_locator="Fig S1A (synthetic)",
            ),
            v010_record(
                method="chemogenetic",
                valence="activation",
                source_locator="Fig S1B (synthetic)",
            ),
            synthetic.other_lab_record(),
        ]
        data = self.records[2].to_dict()
        data["schema_version"] = FROM_VERSION
        for name in ("perturbation_class", "perturbation_polarity", "target_raw", "target_id"):
            data[name] = None
        data["record_id"] = compute_record_id_0_1_0(data)
        self.records[2] = TierARecord.from_dict(data)

    def test_the_map_is_complete_and_bijective_over_the_input(self):
        result = migrate_records(self.records, fixture_index())
        old_ids = [r.record_id for r in self.records]
        self.assertEqual(sorted(result.id_map), sorted(old_ids))
        self.assertEqual(len(set(result.id_map.values())), len(old_ids))
        self.assertEqual(
            [result.id_map[o] for o in old_ids],
            [r.record_id for r in result.records],
        )
        self.assertEqual(result.id_collisions, [])
        self.assertEqual(result.counts["records_in"], len(old_ids))
        self.assertEqual(result.counts["records_out"], len(old_ids))

    def test_every_id_changes_because_the_hash_input_grew(self):
        result = migrate_records(self.records, fixture_index())
        self.assertEqual(result.counts["ids_changed"], len(self.records))
        self.assertEqual(result.counts["ids_unchanged"], 0)
        self.assertEqual(len(RECORD_ID_FIELDS_0_1_0), 7)
        self.assertEqual(len(RECORD_ID_FIELDS), 12)

    def test_a_field_new_to_the_key_changes_the_id_but_not_the_old_one(self):
        base = v010_record(sex="male")
        other = v010_record(sex="female")
        self.assertEqual(base.record_id, other.record_id)  # sex was not in the v0.1.0 key
        index = fixture_index()
        self.assertNotEqual(
            migrate_record(base, index)[0].record_id,
            migrate_record(other, index)[0].record_id,
        )

    def test_conflicts_with_is_remapped(self):
        a, b = self.records[0], self.records[1]
        data = a.to_dict()
        data["conflicts_with"] = [b.record_id]
        data["record_id"] = compute_record_id_0_1_0(data)
        a = TierARecord.from_dict(data)
        result = migrate_records([a, b], fixture_index())
        new_a, new_b = result.records
        self.assertEqual(new_a.conflicts_with, [new_b.record_id])
        self.assertNotEqual(new_a.conflicts_with, [b.record_id])
        self.assertEqual(result.counts["conflicts_with_references_in"], 1)
        self.assertEqual(result.counts["conflicts_with_references_remapped"], 1)
        self.assertEqual(result.counts["conflicts_with_references_unmapped"], 0)
        self.assertFalse(result.needs_a_human)

    def test_an_unmappable_reference_is_kept_and_reported(self):
        a = self.records[0]
        data = a.to_dict()
        data["conflicts_with"] = ["deadbeefdeadbeef"]
        data["record_id"] = compute_record_id_0_1_0(data)
        result = migrate_records([TierARecord.from_dict(data)], fixture_index())
        self.assertEqual(result.records[0].conflicts_with, ["deadbeefdeadbeef"])
        self.assertEqual(result.counts["conflicts_with_references_unmapped"], 1)
        self.assertTrue(result.needs_a_human)

    def test_remap_conflicts_rewrites_in_place(self):
        record = synthetic.good_record(conflicts_with=["old1", "old2"])
        unmapped = remap_conflicts([record], {"old1": "new1"})
        self.assertEqual(record.conflicts_with, ["new1", "old2"])
        self.assertEqual(unmapped, [(record.record_id, "old2")])


class IdempotenceTest(unittest.TestCase):
    def test_migrating_twice_changes_nothing(self):
        index = fixture_index()
        records = [
            v010_record(
                method="conditional_knockout",
                valence="unclear",
                actuator=None,
                targeting_strategy="synthetic Raptor floxed line",
            ),
            v010_record(method="chemogenetic", valence="inhibition"),
        ]
        once = migrate_records(records, index)
        twice = migrate_records(once.records, index)
        self.assertEqual(
            [r.to_dict() for r in once.records], [r.to_dict() for r in twice.records]
        )
        self.assertEqual(twice.counts["ids_changed"], 0)
        self.assertEqual(twice.counts["ids_unchanged"], len(records))
        self.assertEqual(
            twice.counts["records_already_at_target_version"], len(records)
        )
        # The second pass's map is the identity, so composing the two maps is
        # still the migration's own map.
        self.assertEqual(
            {old: twice.id_map[new] for old, new in once.id_map.items()},
            once.id_map,
        )

    def test_a_record_already_at_0_2_0_is_not_double_migrated(self):
        already = synthetic.good_record(
            method="conditional_knockout",
            valence="unclear",
            actuator=None,
            targeting_strategy="synthetic p53 line",
            target_raw="p53",
            target_id="trp53",
            perturbation_class="molecular",
            perturbation_polarity="decrease",
        )
        migrated, resolution = migrate_record(already, fixture_index())
        self.assertEqual(resolution.status, "already_set")
        self.assertEqual(migrated.to_dict(), already.to_dict())
        self.assertEqual(migrated.record_id, already.record_id)


class CliTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.records_path = os.path.join(self.tmp, "curated", "records.jsonl")
        os.makedirs(os.path.dirname(self.records_path))
        self.id_map = os.path.join(self.tmp, "curated", "id_map.json")
        self.backup = os.path.join(self.tmp, "curated", "records.0.1.0.jsonl.bak")
        self.report = os.path.join(self.tmp, "reports", "migration.md")
        self.targets = os.path.join(self.tmp, "targets.json")
        with open(self.targets, "w", encoding="utf-8") as handle:
            json.dump({"_meta": {"vocab": "targets"}, "entries": FIXTURE_TARGETS}, handle)
        self.corpus = [
            v010_record(
                method="conditional_knockout",
                valence="unclear",
                actuator=None,
                targeting_strategy="SynPOPp53 KO (p53 deleted in synthetic neurons)",
                source_locator="Fig S1A (synthetic)",
            ),
            v010_record(
                method="conditional_overexpression",
                valence="unclear",
                actuator="Raptor",
                source_locator="Fig S1B (synthetic)",
            ),
            v010_record(
                method="chemogenetic",
                valence="activation",
                source_locator="Fig S1C (synthetic)",
            ),
            v010_record(
                method="conditional_knockout",
                valence="unclear",
                actuator=None,
                targeting_strategy="synthetic double mutant: p53 and Raptor floxed",
                source_locator="Fig S1D (synthetic)",
            ),
        ]
        write_jsonl_text(self.records_path, self.corpus)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def run_cli(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = main(argv)
        return code, out.getvalue(), err.getvalue()

    def base_argv(self, command="run"):
        argv = [
            command,
            "--records", self.records_path,
            "--id-map", self.id_map,
            "--targets", self.targets,
        ]
        argv += ["--backup", self.backup]
        if command == "run":
            argv += ["--report", self.report]
        return argv

    def read_out(self):
        return [
            json.loads(line)
            for line in read_text(self.records_path).splitlines()
            if line.strip()
        ]

    def test_run_migrates_and_writes_everything(self):
        before = read_bytes(self.records_path)
        code, out, err = self.run_cli(self.base_argv())
        # Exit 1, not 0: one record's target is ambiguous and needs a human.
        self.assertEqual(code, 1, err)
        rows = self.read_out()
        self.assertEqual(len(rows), len(self.corpus))
        self.assertEqual({r["schema_version"] for r in rows}, {TO_VERSION})
        self.assertEqual(len(rows[0]), len(TierARecord.field_names()))

        from syn4rch.schema.io import read_jsonl

        for record in read_jsonl(self.records_path):
            self.assertEqual(
                [f.code for f in validate(record) if f.severity is Severity.HARD],
                [],
                record.record_id,
            )

        document = read_json(self.id_map)
        self.assertEqual(document["_meta"]["from_schema_version"], FROM_VERSION)
        self.assertEqual(document["_meta"]["to_schema_version"], TO_VERSION)
        self.assertEqual(
            sorted(document["map"]), sorted(r.record_id for r in self.corpus)
        )
        self.assertEqual(
            sorted(document["map"].values()), sorted(r["record_id"] for r in rows)
        )
        self.assertEqual(len(set(document["map"].values())), len(rows))

        self.assertEqual(read_bytes(self.backup), before)
        report = read_text(self.report)
        self.assertIn("v0.1.0 -> v0.2.0", report)
        self.assertIn("was not permitted to fill", report)
        self.assertIn("Ambiguous matches", report)
        self.assertIn("target_resolved", out)

    def test_the_file_is_replaced_atomically_and_no_temp_file_is_left(self):
        self.run_cli(self.base_argv())
        leftovers = [
            name
            for name in os.listdir(os.path.dirname(self.records_path))
            if name.endswith(".tmp")
        ]
        self.assertEqual(leftovers, [])

    def test_run_refuses_a_file_that_is_not_at_0_1_0(self):
        write_jsonl_text(
            self.records_path,
            [synthetic.good_record()],  # already at TO_VERSION
        )
        before = read_bytes(self.records_path)
        code, _out, err = self.run_cli(self.base_argv())
        self.assertEqual(code, 2)
        self.assertIn("refusing to migrate", err)
        self.assertIn(TO_VERSION, err)
        self.assertEqual(read_bytes(self.records_path), before)
        self.assertFalse(os.path.exists(self.id_map))
        self.assertFalse(os.path.exists(self.report))
        self.assertFalse(os.path.exists(self.backup))

    def test_run_refuses_an_unknown_version_rather_than_guessing(self):
        data = self.corpus[0].to_dict()
        data["schema_version"] = "0.3.0"
        with open(self.records_path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(json.dumps(data) + "\n")
        code, _out, err = self.run_cli(self.base_argv())
        self.assertEqual(code, 2)
        self.assertIn("0.3.0", err)

    def test_run_refuses_a_mixed_version_file(self):
        mixed = list(self.corpus[:2]) + [synthetic.good_record()]
        write_jsonl_text(self.records_path, mixed)
        code, _out, err = self.run_cli(self.base_argv())
        self.assertEqual(code, 2)
        self.assertIn("refusing to migrate", err)

    def test_run_refuses_an_unreadable_line(self):
        with open(self.records_path, "a", encoding="utf-8") as handle:
            handle.write("{not json\n")
        code, _out, err = self.run_cli(self.base_argv())
        self.assertEqual(code, 2)
        self.assertIn("could not be read", err)

    def test_a_second_run_is_refused_and_changes_nothing(self):
        code, _out, _err = self.run_cli(self.base_argv())
        self.assertEqual(code, 1)
        after_first = read_bytes(self.records_path)
        backup_after_first = read_bytes(self.backup)
        code, _out, err = self.run_cli(self.base_argv())
        self.assertEqual(code, 2)
        self.assertIn("already at %s" % TO_VERSION, err)
        self.assertEqual(read_bytes(self.records_path), after_first)
        self.assertEqual(read_bytes(self.backup), backup_after_first)

    def test_an_existing_different_backup_is_not_overwritten(self):
        with open(self.backup, "w", encoding="utf-8") as handle:
            handle.write("an older pre-migration file\n")
        before = read_bytes(self.records_path)
        code, _out, err = self.run_cli(self.base_argv())
        self.assertEqual(code, 2)
        self.assertIn("not overwritten", err)
        self.assertEqual(read_bytes(self.records_path), before)
        self.assertEqual(
            read_text(self.backup), "an older pre-migration file\n"
        )

    def test_run_refuses_a_stale_input_id(self):
        data = self.corpus[0].to_dict()
        data["record_id"] = "0000000000000000"
        with open(self.records_path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(json.dumps(data) + "\n")
        code, _out, err = self.run_cli(self.base_argv())
        self.assertEqual(code, 2)
        self.assertIn("edited after its ids were written", err)
        self.assertFalse(os.path.exists(self.backup))
        code, _out, _err = self.run_cli(self.base_argv() + ["--allow-stale-input-ids"])
        self.assertIn(code, (0, 1))
        document = read_json(self.id_map)
        self.assertIn("0000000000000000", document["map"])

    def test_dry_run_writes_nothing(self):
        before = read_bytes(self.records_path)
        code, out, _err = self.run_cli(self.base_argv() + ["--dry-run"])
        self.assertIn(code, (0, 1))
        self.assertEqual(read_bytes(self.records_path), before)
        self.assertFalse(os.path.exists(self.id_map))
        self.assertFalse(os.path.exists(self.report))
        self.assertFalse(os.path.exists(self.backup))
        self.assertIn("not written: --dry-run", out)

    def test_missing_file_is_a_usage_error(self):
        code, _out, _err = self.run_cli(
            ["run", "--records", os.path.join(self.tmp, "nope.jsonl"), "--dry-run"]
        )
        self.assertEqual(code, 2)

    def test_report_recomputes_the_published_numbers(self):
        self.run_cli(self.base_argv())
        code, out, err = self.run_cli(
            self.base_argv("report") + ["--compare", self.report]
        )
        self.assertEqual(code, 0, err + out)
        self.assertIn("matches the recomputation", out)
        self.assertIn("reproduces the committed dataset: yes", out)

    def test_report_refuses_a_file_that_is_not_at_0_2_0(self):
        code, _out, err = self.run_cli(self.base_argv("report"))
        self.assertEqual(code, 2)
        self.assertIn("Run the migration first", err)

    def test_report_notices_a_dataset_that_no_longer_matches(self):
        self.run_cli(self.base_argv())
        rows = self.read_out()
        rows.pop()
        with open(self.records_path, "w", encoding="utf-8", newline="\n") as handle:
            for row in rows:
                handle.write(json.dumps(row) + "\n")
        code, out, _err = self.run_cli(
            self.base_argv("report") + ["--compare", self.report]
        )
        self.assertEqual(code, 1)
        self.assertIn("reproduces the committed dataset: NO", out)
        self.assertIn("edited afterwards", out)

    def test_report_without_the_preserved_input_checks_what_it_can(self):
        self.run_cli(self.base_argv())
        os.unlink(self.backup)
        code, out, err = self.run_cli(
            self.base_argv("report") + ["--compare", self.report]
        )
        self.assertEqual(code, 0, err + out)
        self.assertIn("dataset only", out)
        self.assertIn("matches the published report", out)


class ReportContentTest(unittest.TestCase):
    def test_the_report_states_what_was_refused_and_why(self):
        index = fixture_index()
        result = migrate_records(
            [
                v010_record(
                    method="conditional_knockout",
                    valence="unclear",
                    actuator=None,
                    targeting_strategy="synthetic p53 floxed line",
                ),
                v010_record(method="chemogenetic", valence="activation"),
            ],
            index,
        )
        stats = dataset_stats(result.records, index)
        report = render_report_md(
            stats,
            counts=result.counts,
            baseline=None,
            result=result,
            records_path="data/curated/records.jsonl",
            id_map_path="data/curated/map.json",
            backup_path="data/curated/records.0.1.0.jsonl.bak",
            index=index,
            generated="TEST",
        )
        for expected in (
            "records in",
            "record_id` changed",
            "perturbation_class",
            "perturbation_polarity",
            "condition dimension was not filled",
            "was not permitted to fill",
            "`not_reported`",
            "The gene is the study is the lab",
            "conflicts_with",
        ):
            self.assertIn(expected, report, expected)
        self.assertEqual(stats["condition_field_slots_unfilled"], 2 * 4)


if __name__ == "__main__":
    unittest.main()
