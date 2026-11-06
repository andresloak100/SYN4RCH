"""The condition reader: what the Methods can say, and what only a record can.

FABRICATED TEST DATA - NEVER EVIDENCE.

Every XML document below is hand-written in this file and every identifier in
it is a deliberately unissuable form (``PMC0000...``, ``10.0000/...``), in the
style of ``tests/fixtures/jats/``; every ``study_id`` carries the
``synthetic:`` prefix of ``tests/synthetic.py``, so that nothing here can be
mistaken for evidence in a file, a report or a grep. Nothing under ``data/``
is read: the suite must run with no corpus present.

Some sentences are *paraphrases of the shape* of real Methods sentences that
the extractor got wrong (a figure-scoped sex statement, a two-arm diet
disjunction, an age at the start of a protocol). They are reconstructions of
the construction, not transcriptions of anybody's prose, and they assert
nothing about any paper.
"""

import _path  # noqa: F401
import unittest
from unittest import mock
from xml.sax.saxutils import escape

from syn4rch.curate import conditions as C
from syn4rch.curate.jats import parse_jats

STUDY = "synthetic:cond-001"
PMCID = "PMC0000009"


# --------------------------------------------------------------------------- #
# synthetic JATS
# --------------------------------------------------------------------------- #


def jats(methods, results=("Invented result sentence.",), title=None):
    """A minimal fabricated JATS article with the given Methods paragraphs."""
    def paras(items):
        return "\n".join("      <p>%s</p>" % escape(p) for p in items)

    return """<?xml version="1.0" encoding="UTF-8"?>
<article article-type="research-article">
  <front>
    <article-meta>
      <article-id pub-id-type="pmcid">%s</article-id>
      <article-id pub-id-type="doi">10.0000/synthetic.0009</article-id>
      <title-group><article-title>%s</article-title></title-group>
      <pub-date pub-type="epub"><year>2001</year></pub-date>
      <abstract><p>Invented synthetic abstract.</p></abstract>
    </article-meta>
  </front>
  <body>
    <sec sec-type="materials|methods">
      <title>Materials and Methods</title>
%s
    </sec>
    <sec sec-type="results">
      <title>Results</title>
%s
    </sec>
  </body>
</article>
""" % (
        PMCID,
        escape(title or "Synthetic fixture for the condition reader"),
        paras(methods),
        paras(results),
    )


def read(methods, **kwargs):
    paper = parse_jats(jats(methods, **kwargs))
    return C.read_paper_conditions(paper, STUDY, PMCID)


def has_note(result, name):
    return any(note.startswith(name) for note in result.notes)


class Record:
    """A stand-in for the three fields ``read_record_conditions`` reads."""

    def __init__(self, outcome_raw=None, support_quote=None,
                 source_locator="Results para 1",
                 record_id="synthetic:rec-0001"):
        self.outcome_raw = outcome_raw
        self.support_quote = support_quote
        self.source_locator = source_locator
        self.record_id = record_id


# --------------------------------------------------------------------------- #
# the module's own vocabulary
# --------------------------------------------------------------------------- #


class TestModuleVocabulary(unittest.TestCase):
    def test_dimensions_are_a_subset_of_the_schema(self):
        """D-060: a SUBSET, and deliberately so.

        This asserted equality until D-060 added ``background_genotype``. This
        module is the paper-reading extractor for the conditions a Methods
        section states about a cohort; deciding whether an allele is the
        recorded perturbation or the substrate it was applied to is a
        per-record judgement it cannot make, so genotype is curated by its own
        pass. The subset relation is the true claim, and it is asserted rather
        than left as a comment so a dimension cannot silently go unextracted
        without this test naming it.
        """
        from syn4rch.schema.record import CONDITION_DIMENSIONS

        self.assertTrue(set(C.DIMENSIONS) <= set(CONDITION_DIMENSIONS))
        self.assertEqual(
            tuple(d for d in CONDITION_DIMENSIONS if d in set(C.DIMENSIONS)),
            C.DIMENSIONS,
            "the extractor's order must follow the schema's",
        )
        self.assertEqual(
            sorted(set(CONDITION_DIMENSIONS) - set(C.DIMENSIONS)),
            ["background_genotype"],
            "only background_genotype is extracted by another pass (D-060)",
        )

    def test_every_rule_name_is_declared(self):
        # No rule may be emitted that RULES does not name: a value whose rule
        # is unlisted cannot be traced back to the rule that made it.
        self.assertIn("methods_single_sex", C.RULES)
        self.assertIn("quote_names_the_state_with_a_competing_token", C.RULES)
        self.assertNotIn("methods_single_challenge", C.RULES)

    def test_no_methods_rule_for_challenge(self):
        self.assertFalse(any(rule.endswith("_challenge") for rule in C.RULES))


# --------------------------------------------------------------------------- #
# sentences
# --------------------------------------------------------------------------- #


class TestSentences(unittest.TestCase):
    def test_splits_on_a_full_stop_before_a_capital(self):
        self.assertEqual(
            C.sentences("Mice were used. Rats were not."),
            ["Mice were used.", "Rats were not."],
        )

    def test_does_not_split_a_decimal_or_an_abbreviation(self):
        self.assertEqual(len(C.sentences("Significance was p < 0.01 throughout.")), 1)
        self.assertEqual(len(C.sentences("Drugs were given i.p. 30 min before.")), 1)
        self.assertEqual(len(C.sentences("Tracers e.g. CTB were injected.")), 1)
        self.assertEqual(len(C.sentences("As reported by Smith, J. et al.")), 1)

    def test_every_sentence_is_a_verbatim_substring(self):
        text = "Male mice were used. Glucose was measured. Data are mean."
        for sentence in C.sentences(text):
            self.assertIn(sentence, text)

    def test_empty_input(self):
        self.assertEqual(C.sentences(""), [])
        self.assertEqual(C.sentences(None), [])


# --------------------------------------------------------------------------- #
# sex
# --------------------------------------------------------------------------- #


class TestSexFromMethods(unittest.TestCase):
    def test_constant_single_sex(self):
        result = read([
            "All experiments used male C57BL/6J mice.",
            "Mice were group housed.",
        ])
        self.assertEqual(result.verdicts["sex"], "constant")
        extracted = result.constants["sex"]
        self.assertEqual(extracted.value, "male")
        self.assertEqual(extracted.rule, "methods_single_sex")
        self.assertEqual(extracted.source, "methods")
        self.assertEqual(extracted.locator, "Methods para 1")
        self.assertEqual(extracted.quote,
                         "All experiments used male C57BL/6J mice.")
        self.assertEqual(extracted.dimension, "sex")

    def test_constant_both_sexes_by_conjunction(self):
        for sentence in (
            "Male and female C57BL/6J mice were used throughout.",
            "Males and females were used in equal numbers in all mouse cohorts.",
            "Mice of both sexes were studied.",
            "Female and male mice were used throughout.",
            "Male (n = 5) and female (n = 5) mice were used throughout.",
        ):
            with self.subTest(sentence=sentence):
                result = read([sentence])
                self.assertEqual(result.verdicts["sex"], "constant", sentence)
                self.assertEqual(result.constants["sex"].value, "both")
                self.assertEqual(result.constants["sex"].rule,
                                 "methods_both_sexes")

    def test_varies_when_two_sentences_disagree(self):
        result = read([
            "Male C57BL/6J mice were used for the calorimetry measurements.",
            "Female C57BL/6J mice were used for the behavioural measurements.",
        ])
        self.assertEqual(result.verdicts["sex"], "varies")
        self.assertNotIn("sex", result.constants)

    def test_absent_when_no_sentence_names_a_sex(self):
        result = read(["Procedures were approved by the institutional committee."])
        self.assertEqual(result.verdicts["sex"], "absent")
        self.assertNotIn("sex", result.constants)

    def test_breeding_sentences_do_not_vote(self):
        result = read([
            "Female mice used as breeders were maintained in the colony.",
            "Male C57BL/6J mice were used for all experiments.",
        ])
        self.assertEqual(result.verdicts["sex"], "constant")
        self.assertEqual(result.constants["sex"].value, "male")

    # --- defect 5: a sentence scoped to one figure does not vote ---------- #

    def test_figure_scoped_sentence_does_not_vote_and_leaves_absent(self):
        result = read(["Experiments (Figure 3) were performed in female mice only."])
        self.assertEqual(result.verdicts["sex"], "absent")
        self.assertNotIn("sex", result.constants)
        self.assertTrue(has_note(result, "refused_figure_scoped_cohort_sentence"))

    def test_figure_scoped_sentence_blocks_constant_for_the_rest(self):
        # The scoped sentence says a female-only subset exists, so the absence
        # of a second value is not established and 'male' must not be promoted.
        result = read([
            "Male C57BL/6J mice were used.",
            "Experiments (Figure 3) were performed in female mice only.",
        ])
        self.assertEqual(result.verdicts["sex"], "varies")
        self.assertNotIn("sex", result.constants)
        self.assertTrue(
            has_note(result, "refused_sex_not_constant_after_a_refusal")
        )

    # --- defect 6: "male or female" is not both --------------------------- #

    def test_disjunctive_sex_is_varies_not_both(self):
        result = read([
            "Briefly, 3- to 5-week-old male or female mice were deeply "
            "anesthetized with an i.p. injection of chloral hydrate."
        ])
        self.assertEqual(result.verdicts["sex"], "varies")
        self.assertNotIn("sex", result.constants)
        self.assertTrue(
            has_note(result, "refused_disjunctive_sex_is_unresolved")
        )

    def test_two_sex_words_with_no_explicit_construction_is_refused(self):
        result = read([
            "Mice were assigned to groups; male littermates served as "
            "controls for the female measurements reported here."
        ])
        self.assertEqual(result.verdicts["sex"], "varies")
        self.assertNotIn("sex", result.constants)


# --------------------------------------------------------------------------- #
# diet
# --------------------------------------------------------------------------- #


class TestDietFromMethods(unittest.TestCase):
    def test_constant_chow(self):
        result = read(["Mice were maintained on standard chow."])
        self.assertEqual(result.verdicts["diet"], "constant")
        extracted = result.constants["diet"]
        self.assertEqual(extracted.value, "chow")
        self.assertEqual(extracted.rule, "methods_single_diet")
        self.assertEqual(extracted.locator, "Methods para 1")

    def test_constant_high_fat(self):
        result = read([
            "Mice received a high-fat diet (HFD) for the duration of the study."
        ])
        self.assertEqual(result.verdicts["diet"], "constant")
        self.assertEqual(result.constants["diet"].value, "high_fat")

    def test_housing_sentence_is_a_legitimate_source_for_diet(self):
        # Husbandry is inadmissible for challenge and admissible for diet: the
        # maintenance diet is exactly what a housing sentence states.
        result = read([
            "Mice were allowed ad libitum access to water and standard chow."
        ])
        self.assertEqual(result.verdicts["diet"], "constant")
        self.assertEqual(result.constants["diet"].value, "chow")

    def test_varies_when_both_classes_are_named(self):
        result = read([
            "Mice were maintained on standard chow.",
            "A second cohort received a high-fat diet (60% kcal fat).",
        ])
        self.assertEqual(result.verdicts["diet"], "varies")
        self.assertNotIn("diet", result.constants)

    def test_absent_when_no_diet_is_named(self):
        result = read(["Male C57BL/6J mice were used for all experiments."])
        self.assertEqual(result.verdicts["diet"], "absent")

    # --- defect 7: chow as a stimulus is not a maintenance diet ---------- #

    def test_chow_as_a_recording_stimulus_is_not_a_diet(self):
        result = read([
            "Only animals with the expected inhibitory (AgRP) or excitatory "
            "(POMC) response to chow were included for additional analysis."
        ])
        self.assertEqual(result.verdicts["diet"], "absent")
        self.assertNotIn("diet", result.constants)

    def test_chow_presentation_is_not_a_diet(self):
        for sentence in (
            "Chow was presented after the baseline period.",
            "Firing rates were quantified upon chow.",
            "Responses to chow presentation were recorded.",
        ):
            with self.subTest(sentence=sentence):
                result = read([sentence])
                self.assertEqual(result.verdicts["diet"], "absent", sentence)

    # --- defect 8: an unrecognised diet term is not evidence of absence -- #

    def test_an_unrecognised_term_must_not_promote_the_recognised_one(self):
        # One arm is classified and the other is an abbreviation the synonym
        # table does not know. The unknown arm must not go invisible and leave
        # the known one looking unopposed.
        result = read([
            "Mice were maintained on standard chow; a second cohort "
            "received a CD for eight weeks."
        ])
        self.assertEqual(result.verdicts["diet"], "varies")
        self.assertNotIn("diet", result.constants)
        self.assertTrue(has_note(result, "refused_unrecognised_diet_term"))

    def test_a_product_name_is_not_an_unrecognised_diet_type(self):
        # A brand is how a paper identifies its chow, not a second arm.
        for sentence in (
            "Mice had ad libitum access to mouse chow (PicoLab Rodent Diet "
            "20, 5053).",
            "Mice were maintained on standard chow (LabDiet 5021).",
            "Mice were maintained on standard chow (Special Diet Services).",
            "Mice were maintained on standard chow (TestDiet 58Y1).",
        ):
            with self.subTest(sentence=sentence):
                result = read([sentence])
                self.assertEqual(result.verdicts["diet"], "constant", sentence)
                self.assertEqual(result.constants["diet"].value, "chow")
                self.assertFalse(
                    has_note(result, "refused_unrecognised_diet_term"),
                    sentence,
                )

    def test_a_back_reference_is_not_an_unrecognised_diet_type(self):
        for sentence in (
            "Mice were maintained on standard chow; body weight was "
            "measured on this diet.",
            "Mice were maintained on standard chow and each diet was "
            "weighed.",
            "Mice were maintained on standard chow and the diet was replaced "
            "weekly.",
            "Mice were maintained on standard chow; body weight differed "
            "between diet groups.",
        ):
            with self.subTest(sentence=sentence):
                result = read([sentence])
                self.assertEqual(result.verdicts["diet"], "constant", sentence)
                self.assertEqual(result.constants["diet"].value, "chow")

    def test_physical_form_and_category_are_not_diet_types(self):
        for sentence in (
            "Mice were maintained on standard chow as a pelleted diet.",
            "Mice were maintained on standard chow; an irradiated diet was "
            "used.",
            "Mice were maintained on standard chow and a liquid diet was "
            "available.",
        ):
            with self.subTest(sentence=sentence):
                result = read([sentence])
                self.assertEqual(result.verdicts["diet"], "constant", sentence)

    def test_a_category_word_between_regular_and_diet_still_classifies(self):
        result = read([
            "A regular rodent diet (2.99 kcal/g; 9.4% kcal derived from fat) "
            "was provided to the mice."
        ])
        self.assertEqual(result.verdicts["diet"], "constant")
        self.assertEqual(result.constants["diet"].value, "chow")

    def test_a_composition_modifier_the_table_cannot_classify_refuses(self):
        for sentence in (
            "Mice were fed a HFHS diet for twelve weeks.",
            "Mice received 2920x diet from the vendor.",
            "Mice received a control diet.",
            "Mice were fed a soy-free diet.",
        ):
            with self.subTest(sentence=sentence):
                result = read([sentence])
                self.assertEqual(result.verdicts["diet"], "varies", sentence)
                self.assertNotIn("diet", result.constants)
                self.assertTrue(
                    has_note(result, "refused_unrecognised_diet_term"),
                    sentence,
                )

    def test_unrecognised_term_alone_forces_varies(self):
        # No classifiable term at all, but something diet-shaped: the lexicon
        # is incomplete and an incomplete lexicon must refuse.
        result = read(["Mice were fed the indicated diet from weaning."])
        self.assertEqual(result.verdicts["diet"], "varies")
        self.assertTrue(has_note(result, "refused_unrecognised_diet_term"))

    def test_nd_and_hfd_two_arm_study_is_varies(self):
        result = read([
            "Mice were fed either a ND (18% fat, 58% carbohydrates and 24% "
            "proteins) or a HFD (45% fat) for twelve weeks."
        ])
        self.assertEqual(result.verdicts["diet"], "varies")
        self.assertNotIn("diet", result.constants)

    def test_disjunction_forces_varies_on_its_own(self):
        # Both arms happen to classify to the same class here, so only the
        # construction can tell you the paper has two arms.
        result = read([
            "Mice were fed either a standard diet or a normal chow diet."
        ])
        self.assertEqual(result.verdicts["diet"], "varies")
        self.assertTrue(
            has_note(result, "refused_disjunctive_diet_is_two_arms")
        )

    def test_a_disjunction_about_something_else_does_not_force_varies(self):
        # The either/or is about the drug, not the diet: the diet term sits
        # before the "either". Detecting the construction must not mean
        # detecting every construction.
        result = read([
            "Mice fed standard chow received either saline or CNO."
        ])
        self.assertEqual(result.verdicts["diet"], "constant")
        self.assertEqual(result.constants["diet"].value, "chow")

    def test_figure_scoped_diet_sentence_does_not_vote(self):
        result = read([
            "For the experiments in Figure 4, mice received a high-fat diet."
        ])
        self.assertEqual(result.verdicts["diet"], "absent")
        self.assertTrue(has_note(result, "refused_figure_scoped_cohort_sentence"))


# --------------------------------------------------------------------------- #
# age_group
# --------------------------------------------------------------------------- #


class TestAgeFromMethods(unittest.TestCase):
    def test_constant_band_from_a_range(self):
        result = read(["Male C57BL/6J mice aged 8-10 weeks were used."])
        self.assertEqual(result.verdicts["age_group"], "constant")
        extracted = result.constants["age_group"]
        self.assertEqual(extracted.value, "young_adult")
        self.assertEqual(extracted.rule, "methods_single_age_band")
        self.assertEqual(result.species, "mus_musculus")
        # A range states a span, not a value.
        self.assertIsNone(result.age_value)
        self.assertIsNone(result.age_unit)

    def test_en_dash_and_em_dash_ranges_read_the_same(self):
        bands = []
        for dash in ("-", "–", "—"):
            result = read([
                "Male C57BL/6J mice aged 8%s10 weeks were used." % dash
            ])
            bands.append(result.constants["age_group"].value)
        self.assertEqual(bands, ["young_adult"] * 3)

    def test_the_three_statement_shapes(self):
        for sentence in (
            "Male C57BL/6J mice 8-10-week-old were used.",
            "Male C57BL/6J mice aged 8–10 weeks were used.",
            "Male C57BL/6J mice 8-10 weeks of age were used.",
            "Male C57BL/6J mice 4 months old were used.",
            "Male C57BL/6J mice 12 to 16 weeks of age were used.",
        ):
            with self.subTest(sentence=sentence):
                result = read([sentence])
                self.assertEqual(result.verdicts["age_group"], "constant",
                                 sentence)
                self.assertEqual(result.constants["age_group"].value,
                                 "young_adult")

    def test_spelled_out_numerals_read_as_numbers(self):
        # Reading the word "twelve" as 12 is parsing, not inference: the
        # sentence says twelve.
        for sentence, band in (
            ("Twelve-week-old C57BL/6J mice were used.", "young_adult"),
            ("Three-week-old C57BL/6J mice were used.", "juvenile"),
            ("Eight to ten week old C57BL/6J mice were used.",
             "young_adult"),
            ("C57BL/6J mice aged twenty weeks were used.", "young_adult"),
            ("C57BL/6J mice four months old were used.", "young_adult"),
        ):
            with self.subTest(sentence=sentence):
                result = read([sentence])
                self.assertEqual(result.verdicts["age_group"], "constant",
                                 sentence)
                self.assertEqual(result.constants["age_group"].value, band)

    def test_a_spelled_out_numeral_and_its_digits_read_the_same(self):
        words = read(["Twelve-week-old C57BL/6J mice were used."])
        digits = read(["12-week-old C57BL/6J mice were used."])
        self.assertEqual(words.constants["age_group"].value,
                         digits.constants["age_group"].value)
        self.assertEqual(words.age_value, digits.age_value)
        self.assertEqual(words.age_unit, digits.age_unit)

    def test_a_spelled_out_protocol_duration_still_refuses(self):
        result = read([
            "C57BL/6J mice were fed the diet for twelve consecutive weeks "
            "starting at five weeks of age."
        ])
        self.assertEqual(result.verdicts["age_group"], "absent")
        self.assertTrue(
            has_note(result,
                     "refused_protocol_start_age_is_not_measurement_age")
        )

    def test_a_spelled_out_numeral_inside_a_longer_one(self):
        # "seventeen" must not read as the "seven" inside it.
        result = read(["Seventeen-week-old C57BL/6J mice were used."])
        self.assertEqual(result.age_value, 17.0)

    def test_scalar_age_sets_age_value_and_unit(self):
        result = read(["Male C57BL/6J mice were 12 weeks of age when tested."])
        self.assertEqual(result.verdicts["age_group"], "constant")
        self.assertEqual(result.age_value, 12.0)
        self.assertEqual(result.age_unit, "weeks")

    def test_straddling_range_yields_nothing(self):
        # 6 weeks is juvenile and 12 weeks is young_adult, so the range
        # supports neither band.
        result = read(["Male C57BL/6J mice aged 6–12 weeks were used."])
        self.assertEqual(result.verdicts["age_group"], "absent")
        self.assertNotIn("age_group", result.constants)
        self.assertTrue(has_note(result, "refused_straddling_age_range"))

    def test_rat_sentence_is_skipped_not_extrapolated(self):
        result = read([
            "Adult male Sprague-Dawley rats (8-10 weeks old) were used.",
        ])
        self.assertEqual(result.verdicts["age_group"], "absent")
        self.assertNotIn("age_group", result.constants)
        self.assertTrue(has_note(result, "refused_rat_age_band"))

    def test_varies_when_two_bands_are_stated(self):
        result = read([
            "Male C57BL/6J mice aged 8-10 weeks were used for cohort one.",
            "Male C57BL/6J mice aged 70-80 weeks were used for cohort two.",
        ])
        self.assertEqual(result.verdicts["age_group"], "varies")
        self.assertNotIn("age_group", result.constants)

    def test_absent_when_no_age_is_stated(self):
        result = read(["Male C57BL/6J mice were used for all experiments."])
        self.assertEqual(result.verdicts["age_group"], "absent")

    def test_two_ages_joined_by_and_are_two_cohorts(self):
        result = read([
            "Male C57BL/6J mice aged 15 and 78 weeks old were compared."
        ])
        self.assertEqual(result.verdicts["age_group"], "varies")
        self.assertIsNone(result.age_value)

    def test_in_vitro_age_is_not_the_animals_age(self):
        result = read([
            "Mouse cortical neurons were recorded aged 13–17 days after "
            "plating in culture."
        ])
        self.assertEqual(result.verdicts["age_group"], "absent")
        self.assertTrue(has_note(result, "refused_in_vitro_age_statement"))

    # --- defect 2: a protocol-start age is not a measurement age --------- #

    def test_protocol_start_age_is_refused(self):
        result = read([
            "C57BL/6 mice were fed with HFD for 12 consecutive weeks "
            "(starting at 5 weeks of age)."
        ])
        self.assertEqual(result.verdicts["age_group"], "absent")
        self.assertNotIn("age_group", result.constants)
        self.assertIsNone(result.age_value)
        self.assertTrue(
            has_note(result,
                     "refused_protocol_start_age_is_not_measurement_age")
        )

    def test_protocol_duration_cue_alone_refuses(self):
        result = read([
            "Male C57BL/6J mice aged 8 weeks were maintained on the diet for "
            "12 weeks before measurement."
        ])
        self.assertEqual(result.verdicts["age_group"], "absent")
        self.assertTrue(
            has_note(result,
                     "refused_protocol_start_age_is_not_measurement_age")
        )

    def test_an_age_stated_at_surgery_is_refused(self):
        # The age at stereotaxic injection is the age weeks before the
        # readout, while the virus expresses. Same class as the protocol-start
        # age: an age that is not the age at the measurement.
        result = read([
            "7-9-month-old male Agrp-cre mice were anesthetized with "
            "isoflurane and placed in a stereotaxic frame."
        ])
        self.assertEqual(result.verdicts["age_group"], "absent")
        self.assertTrue(
            has_note(result,
                     "refused_protocol_start_age_is_not_measurement_age")
        )

    def test_the_same_surgery_sentence_is_still_valid_for_sex(self):
        # The asymmetry is the point: an animal's sex does not change between
        # surgery and readout, and its age does.
        result = read([
            "7-9-month-old male Agrp-cre mice were anesthetized with "
            "isoflurane and placed in a stereotaxic frame."
        ])
        self.assertEqual(result.verdicts["sex"], "constant")
        self.assertEqual(result.constants["sex"].value, "male")

    def test_an_age_on_arrival_is_refused(self):
        for sentence in (
            "Adult male C57BL/6J mice aged 12 weeks on arrival were obtained "
            "from the supplier.",
            "Male C57BL/6 mice (6-7 weeks old with initial body weight of "
            "20-22 g) were purchased from the vendor.",
            "Male C57BL/6J mice were acclimatised for one week at 8 weeks of "
            "age.",
        ):
            with self.subTest(sentence=sentence):
                result = read([sentence])
                self.assertEqual(result.verdicts["age_group"], "absent",
                                 sentence)

    def test_an_age_at_weaning_is_refused(self):
        result = read([
            "The high-fat diet was fed to mice after weaning at 3 weeks of age."
        ])
        self.assertEqual(result.verdicts["age_group"], "absent")

    def test_an_age_stated_for_one_named_assay_is_refused(self):
        for sentence in (
            "For acute cold exposure, 8-week-old mice were placed at 4 C "
            "for 6 h.",
            "For the diet-induced obesity experiments, 4-week-old mice were "
            "used.",
        ):
            with self.subTest(sentence=sentence):
                result = read([sentence])
                self.assertEqual(result.verdicts["age_group"], "absent",
                                 sentence)

    def test_a_perfusion_age_is_refused(self):
        result = read([
            "Briefly, 3- to 5-week-old mice were deeply anesthetized and "
            "trans-cardially perfused."
        ])
        self.assertEqual(result.verdicts["age_group"], "absent")
        self.assertNotIn("age_group", result.constants)

    def test_a_paper_level_age_statement_still_reads(self):
        # The refusals above must not swallow the sentences that do state the
        # age of the animals the paper measured.
        for sentence, band in (
            ("All experiments were performed in adult mice, aged 12-20 weeks.",
             "young_adult"),
            ("All mice were aged 10-18 weeks at the time when the "
             "experiments were performed.", "young_adult"),
            ("Mice aged 2 to 3 months were used.", "young_adult"),
        ):
            with self.subTest(sentence=sentence):
                result = read([sentence])
                self.assertEqual(result.verdicts["age_group"], "constant",
                                 sentence)
                self.assertEqual(result.constants["age_group"].value, band)

    # --- defect 3: a refusal does not promote the remainder -------------- #

    def test_a_refused_statement_blocks_constant_for_the_survivors(self):
        result = read([
            "Male Sprague Dawley rats 8-10-week-old were used for the "
            "electrophysiology.",
            "Male C57BL/6J mice aged 12 weeks were used for the calorimetry.",
        ])
        self.assertEqual(result.verdicts["age_group"], "absent")
        self.assertNotIn("age_group", result.constants)

    def test_mouse_only_paper_with_a_refusal_is_varies_not_constant(self):
        result = read([
            "Male C57BL/6J mice aged 6–12 weeks were used for the "
            "electrophysiology.",
            "Male C57BL/6J mice aged 12 weeks were used for the calorimetry.",
        ])
        self.assertEqual(result.verdicts["age_group"], "varies")
        self.assertNotIn("age_group", result.constants)
        self.assertTrue(
            has_note(result, "refused_age_not_constant_after_a_refusal")
        )

    # --- defect 4: the mouse precondition is explicit -------------------- #

    def test_no_band_when_species_is_not_established_as_mouse(self):
        # Both mice and rats are named, so species is unresolved and no mouse
        # band may be assigned, however clean the sentence looks.
        result = read([
            "Male Sprague Dawley rats and C57BL/6 mice were used.",
            "Mice were 12 weeks of age when tested.",
        ])
        self.assertIsNone(result.species)
        self.assertEqual(result.verdicts["age_group"], "absent")
        self.assertNotIn("age_group", result.constants)
        self.assertIsNone(result.age_value)
        self.assertTrue(has_note(result, "refused_age_band_requires_mouse"))

    def test_rat_only_paper_gets_no_band(self):
        result = read([
            "Male Wistar rats were used.",
            "Rats were 12 weeks of age when tested.",
        ])
        self.assertEqual(result.species, "rattus_norvegicus")
        self.assertEqual(result.verdicts["age_group"], "absent")
        self.assertIsNone(result.age_value)

    def test_assign_age_group_refusal_is_never_stored_as_a_value(self):
        for result in (
            read(["Male Wistar rats were 12 weeks of age."]),
            read(["Male C57BL/6J mice aged 6–12 weeks were used."]),
        ):
            self.assertNotIn("age_group", result.constants)
            for extracted in result.constants.values():
                self.assertNotEqual(extracted.value, "not_reported")


# --------------------------------------------------------------------------- #
# challenge, at paper level
# --------------------------------------------------------------------------- #


class TestChallengeFromMethods(unittest.TestCase):
    def test_housing_statement_is_refused(self):
        result = read([
            "Mice were allowed ad libitum access to water and standard chow."
        ])
        self.assertEqual(result.verdicts["challenge"], "absent")
        self.assertNotIn("challenge", result.constants)
        self.assertTrue(
            has_note(result,
                     "refused_housing_statement_is_not_a_challenge")
        )

    def test_a_surgical_preparation_sentence_is_refused(self):
        # Defect 1: this routed around the housing refusal because it is not
        # husbandry. The animals were ad libitum fed when they were
        # anaesthetised, which says nothing about any measurement.
        result = read([
            "Ad libitum fed rats were anesthetized by an i.p. injection of "
            "ketamine before stereotaxic injection."
        ])
        self.assertNotEqual(result.verdicts["challenge"], "constant")
        self.assertNotIn("challenge", result.constants)
        self.assertTrue(
            has_note(result,
                     "refused_methods_cannot_state_a_measurement_state")
        )

    def test_a_single_named_state_is_varies_never_constant(self):
        result = read([
            "Food intake was measured in ad libitum fed mice over 24 h."
        ])
        self.assertEqual(result.verdicts["challenge"], "varies")
        self.assertNotIn("challenge", result.constants)

    def test_varies_when_the_methods_name_several_states(self):
        result = read([
            "Mice were allowed ad libitum access to water and standard chow.",
            "Blood glucose was measured after an overnight fast.",
        ])
        self.assertEqual(result.verdicts["challenge"], "varies")
        self.assertNotIn("challenge", result.constants)

    def test_absent_when_no_state_is_named(self):
        result = read(["Male C57BL/6J mice were used for all experiments."])
        self.assertEqual(result.verdicts["challenge"], "absent")

    def test_challenge_is_never_constant(self):
        for methods in (
            ["Blood glucose was measured after an overnight fast."],
            ["Mice were fasted overnight before every measurement."],
            ["Energy expenditure was measured during cold exposure."],
            ["Mice were allowed ad libitum access to food and water."],
        ):
            with self.subTest(methods=methods):
                result = read(methods)
                self.assertNotEqual(result.verdicts["challenge"], "constant")
                self.assertNotIn("challenge", result.constants)


# --------------------------------------------------------------------------- #
# challenge, from the record itself
# --------------------------------------------------------------------------- #


class TestRecordChallenge(unittest.TestCase):
    def one(self, **kwargs):
        notes = []
        extracted = C.read_record_conditions(Record(**kwargs), notes)
        return extracted, notes

    def test_outcome_names_the_state(self):
        for outcome, expected in (
            ("EE in the fasted state", "fasted"),
            ("food intake under fed conditions", "ad_libitum"),
            ("food intake following an overnight fast", "fasted"),
            ("fasting-induced lipolysis", "fasted"),
            ("locomotor activity in the fed state", "ad_libitum"),
            ("RQ under fasting conditions", "fasted"),
            ("fasted-induced hyperphagia", "fasted"),
            ("food intake after refeeding", "refed"),
            ("EE during cold exposure", "cold_exposure"),
            ("BAT temperature at thermoneutrality", "thermoneutral"),
            ("corticosterone after restraint stress", "restraint_stress"),
            ("VO2 during treadmill running", "exercise"),
        ):
            with self.subTest(outcome=outcome):
                extracted, _ = self.one(outcome_raw=outcome)
                self.assertEqual(len(extracted), 1, outcome)
                self.assertEqual(extracted[0].value, expected)
                self.assertEqual(extracted[0].dimension, "challenge")
                self.assertEqual(extracted[0].source, "passage")
                self.assertEqual(extracted[0].rule, "outcome_names_the_state")
                # The quote is outcome_raw itself, verbatim.
                self.assertEqual(extracted[0].quote, outcome)
                self.assertEqual(extracted[0].locator, "Results para 1")

    def test_fed_inside_a_diet_phrase_is_not_a_fed_state(self):
        for outcome in (
            "RER (HFD-fed animals)",
            "food intake in chow-fed mice",
            "body weight in high-fat-fed mice",
            "glucose tolerance after 1 week of fat feeding",
            "body weight in mice fed a high-fat diet",
            "body weight in mice fed standard chow",
            "food intake (Western diet-fed)",
        ):
            with self.subTest(outcome=outcome):
                extracted, _ = self.one(outcome_raw=outcome)
                self.assertEqual(extracted, [], outcome)

    def test_two_states_in_one_outcome_is_a_refusal(self):
        extracted, notes = self.one(outcome_raw="fed and fasting serum insulin")
        self.assertEqual(extracted, [])
        self.assertTrue(
            any(n.startswith("refused_two_states_in_one_outcome") for n in notes)
        )

    def test_assay_names_yield_no_challenge(self):
        for outcome in (
            "GTT",
            "ITT",
            "glucose tolerance",
            "insulin tolerance test",
            "glucose tolerance test",
            "insulin excursions during the ipGTT",
            "blood glucose levels in the ITT",
        ):
            with self.subTest(outcome=outcome):
                extracted, notes = self.one(
                    outcome_raw=outcome,
                    support_quote="Glucose fell in fasted mice.",
                )
                self.assertEqual(extracted, [], outcome)
                self.assertTrue(
                    any(
                        n.startswith(
                            "refused_assay_implies_its_own_challenge"
                        )
                        for n in notes
                    )
                )

    def test_quote_names_the_state_when_exactly_one_phrase_occurs(self):
        quote = "Food intake was increased in male mice fasted for 24 h."
        extracted, _ = self.one(outcome_raw="food intake", support_quote=quote)
        self.assertEqual(len(extracted), 1)
        self.assertEqual(extracted[0].value, "fasted")
        self.assertEqual(extracted[0].rule, "quote_names_the_state")
        self.assertEqual(extracted[0].quote, quote)

    def test_quote_sentence_is_the_one_carrying_the_phrase(self):
        quote = ("Body weight was unchanged. Food intake rose in mice fasted "
                 "for 24 h.")
        extracted, _ = self.one(outcome_raw="food intake", support_quote=quote)
        self.assertEqual(
            extracted[0].quote,
            "Food intake rose in mice fasted for 24 h.",
        )
        self.assertIn(extracted[0].quote, quote)

    def test_quote_with_two_state_phrases_is_refused(self):
        quote = ("Food intake rose in fasted mice but not in ad libitum fed "
                 "mice.")
        extracted, notes = self.one(outcome_raw="food intake",
                                    support_quote=quote)
        self.assertEqual(extracted, [])
        self.assertTrue(
            any(
                n.startswith("refused_quote_names_more_than_one_state")
                for n in notes
            )
        )

    def test_outcome_wins_over_the_quote(self):
        extracted, _ = self.one(
            outcome_raw="EE in the fasted state",
            support_quote="Energy expenditure rose during cold exposure.",
        )
        self.assertEqual(len(extracted), 1)
        self.assertEqual(extracted[0].value, "fasted")
        self.assertEqual(extracted[0].rule, "outcome_names_the_state")

    def test_no_locator_is_a_refusal(self):
        extracted, notes = self.one(outcome_raw="EE in the fasted state",
                                    source_locator=None)
        self.assertEqual(extracted, [])
        self.assertTrue(
            any(n.startswith("refused_record_has_no_locator") for n in notes)
        )

    def test_nothing_is_read_for_the_other_three_dimensions(self):
        extracted, _ = self.one(
            outcome_raw="food intake",
            support_quote="Food intake rose in male chow-fed mice aged 10 "
                          "weeks that had been fasted.",
        )
        self.assertEqual({e.dimension for e in extracted}, {"challenge"})

    # --- defect 9: pair-feeding has no enum member ----------------------- #

    def test_pair_fed_is_refused_not_mapped_to_ad_libitum(self):
        quote = ("Importantly, mice injected with the inhibitory actuator "
                 "maintained lower glucose levels following CNO treatment in "
                 "pair-fed mice.")
        extracted, notes = self.one(outcome_raw="blood glucose",
                                    support_quote=quote)
        self.assertEqual(extracted, [])
        self.assertTrue(
            any(
                n.startswith(
                    "refused_restricted_feeding_has_no_enum_member"
                )
                for n in notes
            )
        )

    def test_every_restricted_feeding_wording_is_refused(self):
        for phrase in (
            "pair-fed", "pair fed", "food-restricted", "food restricted",
            "calorie-restricted", "calorie restricted", "restricted-fed",
            "time-restricted",
        ):
            with self.subTest(phrase=phrase):
                extracted, _ = self.one(
                    outcome_raw="body weight",
                    support_quote="Body weight fell in %s mice." % phrase,
                )
                self.assertEqual(extracted, [], phrase)

    def test_restricted_feeding_is_never_mapped_to_other(self):
        extracted, _ = self.one(
            outcome_raw="body weight",
            support_quote="Body weight fell in pair-fed mice.",
        )
        self.assertEqual([e.value for e in extracted], [])

    # --- defect 10: a qualifier inverts the state it governs ------------- #

    def test_below_thermoneutrality_is_refused(self):
        quote = ("After ablation, however, mice exhibited reduced TC defense "
                 "when exposed to temperatures below thermoneutrality.")
        extracted, notes = self.one(outcome_raw="core temperature",
                                    support_quote=quote)
        self.assertEqual(extracted, [])
        self.assertTrue(
            any(
                n.startswith("refused_qualifier_inverts_the_state")
                for n in notes
            )
        )

    def test_below_thermoneutrality_is_not_converted_to_cold_exposure(self):
        extracted, _ = self.one(
            outcome_raw="core temperature",
            support_quote="Mice were held below thermoneutrality.",
        )
        self.assertEqual([e.value for e in extracted], [])

    def test_other_inverting_qualifiers(self):
        for quote in (
            "Measurements were made in mice that were not fasted.",
            "Housing was above thermoneutrality throughout.",
            "Animals were kept at sub-thermoneutral temperatures.",
            "Mice were studied outside thermoneutrality.",
        ):
            with self.subTest(quote=quote):
                extracted, _ = self.one(outcome_raw="core temperature",
                                        support_quote=quote)
                self.assertEqual(extracted, [], quote)

    def test_a_qualifier_must_govern_the_term_directly(self):
        # "did not suppress refeeding" has a word between the qualifier and
        # the term, so the qualifier does not govern it.
        extracted, _ = self.one(
            outcome_raw="food intake",
            support_quote="The dose did not suppress refeeding in littermates.",
        )
        self.assertEqual(len(extracted), 1)
        self.assertEqual(extracted[0].value, "refed")

    # --- defect 11: flag a competing token, do not refuse it ------------- #

    def test_a_competing_excluded_token_is_flagged_not_refused(self):
        quote = ("chow-fed knockout mice exhibited an anorectic response to a "
                 "sub-anorectic dose of the agonist that otherwise did not "
                 "suppress refeeding in wild-type littermates.")
        extracted, notes = self.one(outcome_raw="food intake",
                                    support_quote=quote)
        self.assertEqual(len(extracted), 1)
        self.assertEqual(extracted[0].value, "refed")
        self.assertEqual(
            extracted[0].rule,
            "quote_names_the_state_with_a_competing_token",
        )
        flag = [
            n for n in notes
            if n.startswith("quote_names_the_state_with_a_competing_token")
        ]
        self.assertEqual(len(flag), 1)
        self.assertIn("refeeding", flag[0])
        self.assertIn("chow-fed", flag[0])

    def test_a_clean_quote_keeps_the_plain_rule_name(self):
        extracted, _ = self.one(
            outcome_raw="food intake",
            support_quote="Food intake rose in mice fasted for 24 h.",
        )
        self.assertEqual(extracted[0].rule, "quote_names_the_state")

    def test_the_outcome_route_flags_a_competing_token_too(self):
        extracted, notes = self.one(outcome_raw="RER (HFD-fed, fasted)")
        self.assertEqual(len(extracted), 1)
        self.assertEqual(extracted[0].value, "fasted")
        self.assertEqual(
            extracted[0].rule,
            "outcome_names_the_state_with_a_competing_token",
        )


# --------------------------------------------------------------------------- #
# the verbatim gate
# --------------------------------------------------------------------------- #


class TestVerbatimGate(unittest.TestCase):
    def test_verify_accepts_a_sentence_from_the_paper(self):
        paper = parse_jats(jats(["Male C57BL/6J mice were used."]))
        extracted = C.Extracted(
            dimension="sex", value="male", source="methods",
            locator="Methods para 1",
            quote="Male C57BL/6J mice were used.",
            rule="methods_single_sex",
        )
        self.assertTrue(C.verify_against_paper(extracted, paper))

    def test_verify_rejects_a_sentence_the_paper_does_not_contain(self):
        paper = parse_jats(jats(["Male C57BL/6J mice were used."]))
        for quote in (
            "Female C57BL/6J mice were used.",
            "male C57BL/6J mice were used.",  # case is not folded
            "",
        ):
            with self.subTest(quote=quote):
                extracted = C.Extracted(
                    dimension="sex", value="male", source="methods",
                    locator="Methods para 1", quote=quote,
                    rule="methods_single_sex",
                )
                self.assertFalse(C.verify_against_paper(extracted, paper))

    def test_every_extraction_a_real_read_produces_is_verbatim(self):
        paper = parse_jats(jats([
            "Male C57BL/6J mice aged 8-10 weeks were maintained on standard "
            "chow.",
        ]))
        result = C.read_paper_conditions(paper, STUDY, PMCID)
        self.assertTrue(result.constants)
        for extracted in result.constants.values():
            self.assertTrue(C.verify_against_paper(extracted, paper))

    def test_a_fabricated_quote_is_dropped_and_the_verdict_demoted(self):
        # The gate must hold even when a rule above it is wrong, so the rule
        # is made wrong: ``sentences`` is replaced by one that returns a
        # sentence the paper does not contain. Nothing reaches ``constants``.
        paper = parse_jats(jats(["Male C57BL/6J mice were used."]))
        fabricated = "Male mice were used, in a sentence the paper lacks."
        with mock.patch.object(C, "sentences", lambda text: [fabricated]):
            result = C.read_paper_conditions(paper, STUDY, PMCID)
        self.assertEqual(result.verdicts["sex"], "absent")
        self.assertNotIn("sex", result.constants)
        self.assertTrue(has_note(result, "dropped_quote_not_verbatim"))

    def test_a_dropped_age_extraction_clears_age_value(self):
        paper = parse_jats(jats(["Male C57BL/6J mice were 12 weeks of age."]))
        fabricated = "Male C57BL/6J mice were 12 weeks of age, allegedly."
        with mock.patch.object(C, "sentences", lambda text: [fabricated]):
            result = C.read_paper_conditions(paper, STUDY, PMCID)
        self.assertEqual(result.verdicts["age_group"], "absent")
        self.assertIsNone(result.age_value)
        self.assertIsNone(result.age_unit)


# --------------------------------------------------------------------------- #
# the shape of the result
# --------------------------------------------------------------------------- #


class TestPaperConditionsShape(unittest.TestCase):
    def test_all_four_dimensions_get_a_verdict(self):
        result = read(["Male C57BL/6J mice were used."])
        self.assertEqual(set(result.verdicts), set(C.DIMENSIONS))
        for verdict in result.verdicts.values():
            self.assertIn(verdict, C.VERDICTS)

    def test_constants_only_holds_constant_dimensions(self):
        result = read([
            "Male C57BL/6J mice aged 8-10 weeks were maintained on standard "
            "chow.",
            "Blood glucose was measured after an overnight fast.",
        ])
        for dimension, extracted in result.constants.items():
            self.assertEqual(result.verdicts[dimension], "constant")
            self.assertEqual(extracted.dimension, dimension)
        for dimension, verdict in result.verdicts.items():
            if verdict == "constant":
                self.assertIn(dimension, result.constants)

    def test_identity_is_carried_through(self):
        result = read(["Male C57BL/6J mice were used."])
        self.assertEqual(result.study_id, STUDY)
        self.assertEqual(result.pmcid, PMCID)

    def test_a_paper_with_no_methods_reads_nothing(self):
        xml = """<?xml version="1.0"?><article><front><article-meta>
        <title-group><article-title>Synthetic</article-title></title-group>
        </article-meta></front><body><sec sec-type="results">
        <title>Results</title><p>Male mice were used and food intake rose.</p>
        </sec></body></article>"""
        result = C.read_paper_conditions(parse_jats(xml), STUDY, PMCID)
        self.assertEqual(set(result.verdicts.values()), {"absent"})
        self.assertEqual(result.constants, {})
        self.assertTrue(result.notes)

    def test_every_extraction_source_is_a_condition_source_value(self):
        from syn4rch.schema.enums import ConditionSource

        allowed = {member.value for member in ConditionSource}
        result = read([
            "Male C57BL/6J mice aged 8-10 weeks were maintained on standard "
            "chow."
        ])
        for extracted in result.constants.values():
            self.assertIn(extracted.source, allowed)
            self.assertIn(extracted.rule, C.RULES)

    def test_extractions_convert_to_condition_evidence(self):
        from syn4rch.schema.record import ConditionEvidence

        result = read(["Male C57BL/6J mice were used."])
        extracted = result.constants["sex"]
        evidence = ConditionEvidence(
            dimension=extracted.dimension,
            value=extracted.value,
            source=extracted.source,
            locator=extracted.locator,
            quote=extracted.quote,
            rule=extracted.rule,
        )
        self.assertEqual(evidence.source.value, "methods")
        self.assertEqual(evidence.value, "male")

    def test_locators_parse_as_project_locators(self):
        from syn4rch.curate.passages import parse_locator

        result = read([
            "Procedures were approved by the committee.",
            "Male C57BL/6J mice were maintained on standard chow.",
        ])
        self.assertTrue(result.constants)
        for extracted in result.constants.values():
            reference = parse_locator(extracted.locator)
            self.assertIsNotNone(reference, extracted.locator)
            self.assertEqual(reference.kind, "section_para")
            self.assertEqual(reference.name, "Methods")
            self.assertEqual(reference.number, 2)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
