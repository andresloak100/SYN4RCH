"""Controlled enumerations for the SYN4RCH Tier A record schema.

Every enum here mirrors a vocabulary fixed by ``docs/SCHEMA.md`` (v0.2.0).
The *string values* are the contract: they are what is written to
``data/curated/records.jsonl`` and what extraction prompts must emit. Member
names are an implementation convenience and may be renamed; values may not.

Python 3.10 compatible, standard library only.
"""

from enum import Enum

__all__ = [
    "RetrievedFrom",
    "ExtractionMethod",
    "ReviewStatus",
    "Neurotransmitter",
    "Species",
    "Sex",
    "AgeUnit",
    "AgeGroup",
    "Diet",
    "Challenge",
    "ConditionSource",
    "PerturbationMethod",
    "Valence",
    "PerturbationClass",
    "PerturbationPolarity",
    "Acuity",
    "Laterality",
    "OutcomeTier",
    "Direction",
    "EffectSizeType",
    "DispersionType",
    "ControlType",
    "EvidenceClass",
    "Confidence",
    "TimeUnit",
]


class _StrEnum(str, Enum):
    """``str``-valued enum.

    Subclassing ``str`` keeps ``json.dumps`` and f-strings honest without
    needing ``enum.StrEnum`` (3.11+). ``__str__`` is pinned to the value so
    that reports do not print ``Direction.increase``.
    """

    def __str__(self) -> str:  # pragma: no cover - trivial
        return str(self.value)


# --- 3.1 identity and provenance -------------------------------------------

class RetrievedFrom(_StrEnum):
    pubmed = "pubmed"
    europepmc_oa = "europepmc_oa"
    europepmc_abstract = "europepmc_abstract"
    crossref = "crossref"
    manual = "manual"


class ExtractionMethod(_StrEnum):
    human = "human"
    llm = "llm"
    llm_human_reviewed = "llm_human_reviewed"


class ReviewStatus(_StrEnum):
    unreviewed = "unreviewed"
    accepted = "accepted"
    rejected = "rejected"
    needs_full_text = "needs_full_text"


# --- 3.2 neuronal population ----------------------------------------------

class Neurotransmitter(_StrEnum):
    GABA = "GABA"
    glutamate = "glutamate"
    dopamine = "dopamine"
    serotonin = "serotonin"
    histamine = "histamine"
    peptidergic = "peptidergic"
    mixed = "mixed"
    unknown = "unknown"


# --- 3.3 subject -----------------------------------------------------------

class Species(_StrEnum):
    mus_musculus = "mus_musculus"
    rattus_norvegicus = "rattus_norvegicus"
    homo_sapiens = "homo_sapiens"
    other = "other"


class Sex(_StrEnum):
    male = "male"
    female = "female"
    both = "both"
    not_reported = "not_reported"


class AgeUnit(_StrEnum):
    days = "days"
    weeks = "weeks"
    months = "months"
    years = "years"


class AgeGroup(_StrEnum):
    juvenile = "juvenile"
    young_adult = "young_adult"
    middle_aged = "middle_aged"
    aged = "aged"
    not_reported = "not_reported"


class Diet(_StrEnum):
    chow = "chow"
    high_fat = "high_fat"
    other = "other"
    not_reported = "not_reported"


class Challenge(_StrEnum):
    """The physiological state the readout was taken under (v0.2.0).

    Distinct from :class:`Diet`, which is the *maintenance* diet over weeks.
    ``Challenge`` is the acute state at the moment of measurement, and in this
    literature it frequently decides the sign of the result: AgRP activation
    raises food intake in a fed mouse and raises it far more in a fasted one,
    and a glucose reading means opposite things before and after a glucose
    load. Before v0.2.0 the schema had no field for it, so two incomparable
    experiments could occupy one prediction cell and their disagreement was
    recorded as a contradiction between laboratories.

    ``not_reported`` and ``None`` are different statements. ``None`` means the
    curator did not address the question; ``not_reported`` means the paper was
    checked and does not say. See docs/DECISIONS.md (D-023).
    """

    ad_libitum = "ad_libitum"
    fasted = "fasted"
    refed = "refed"
    cold_exposure = "cold_exposure"
    thermoneutral = "thermoneutral"
    glucose_load = "glucose_load"
    insulin_load = "insulin_load"
    high_fat_challenge = "high_fat_challenge"
    restraint_stress = "restraint_stress"
    exercise = "exercise"
    other = "other"
    not_reported = "not_reported"


class ConditionSource(_StrEnum):
    """Where one condition dimension of a record came from (v0.2.0).

    SCHEMA.md 2 rule 4 requires that a field the record's own passage does not
    report be ``null``. Conditions are a property of the *experiment*, not of
    the sentence, and in practice they are stated in Methods while the result
    is stated in Results. Forbidding the Methods route left 197 of 362 records
    with no diet and 224 with no sex, which is what made the condition gap
    unfixable by curating harder.

    ``methods`` permits the Methods route and makes it auditable. ``passage``
    means the record's own quoted passage states the value. See
    docs/DECISIONS.md (D-024).

    v0.3.0 moved this from a record-level field to
    :class:`~syn4rch.schema.record.ConditionEvidence`, one entry per filled
    dimension, because a record's sex and its acute metabolic state routinely
    arrive by different routes and one scalar could not say so (D-037). The
    locator moved with it and is now required on **every** entry, whatever its
    source (``A008``): the sentence stating the diet is not the sentence the
    effect was read from even when both are in Results, so ``source_locator``
    does not locate it. The verbatim gate applies to an entry's locator as it
    does to any other.
    """

    passage = "passage"
    methods = "methods"
    not_reported = "not_reported"


# --- 3.4 perturbation ------------------------------------------------------

class PerturbationMethod(_StrEnum):
    optogenetic = "optogenetic"
    chemogenetic = "chemogenetic"
    electrical_stimulation = "electrical_stimulation"
    lesion = "lesion"
    genetic_ablation = "genetic_ablation"
    conditional_knockout = "conditional_knockout"
    conditional_overexpression = "conditional_overexpression"
    local_pharmacology = "local_pharmacology"
    other = "other"


class Valence(_StrEnum):
    """Net effect of the perturbation on this population's *activity*.

    ``unclear`` is a correct answer and SCHEMA.md 3.4 forbids coercing it. It
    is what a conditional knockout of a receptor gene honestly deserves: the
    manipulation changes a gene product, and whether firing went up or down is
    usually not measured. v0.2.0 therefore stops asking ``valence`` to carry
    the molecular manipulations as well, and gives them
    :class:`PerturbationClass` and :class:`PerturbationPolarity` instead.
    """

    activation = "activation"
    inhibition = "inhibition"
    ablation = "ablation"
    unclear = "unclear"


class PerturbationClass(_StrEnum):
    """What *kind* of thing was done to the population (v0.2.0).

    A function of ``method``; see
    :func:`syn4rch.schema.rules.perturbation_class_from_method`. It exists as
    a stored field, rather than being recomputed at every use, because it
    qualifies the meaning of :class:`PerturbationPolarity` and the two must
    never travel apart.

    * ``circuit_activity`` -- firing was directly driven or silenced
      (optogenetic, chemogenetic, electrical).
    * ``molecular`` -- a gene product inside the population was removed or
      added (conditional knockout, conditional overexpression). The effect on
      firing is a separate, usually unmeasured question.
    * ``structural`` -- the population itself was destroyed (lesion, genetic
      ablation).
    * ``pharmacological`` -- a drug was applied locally; the sign depends on
      the drug.
    * ``other`` -- declared rather than guessed.
    """

    circuit_activity = "circuit_activity"
    molecular = "molecular"
    structural = "structural"
    pharmacological = "pharmacological"
    other = "other"


class PerturbationPolarity(_StrEnum):
    """The sign of the manipulation **along its own class's axis** (v0.2.0).

    Read this with :class:`PerturbationClass` or not at all. ``decrease``
    means *less firing* for ``circuit_activity`` and *less of the target gene
    product* for ``molecular``. Those are different physical interventions and
    pooling them would be the same sign error that
    ``docs/DECISIONS.md`` D-019 was written against. The pairing is enforced
    structurally: :func:`syn4rch.qc.keys.perturbation_key` returns the class,
    the polarity and the target together, and nothing in this project is
    permitted to key on polarity alone.

    * ``increase`` -- activation, or gain of function.
    * ``decrease`` -- inhibition, or loss of function.
    * ``remove`` -- the population or the gene product is gone outright.
    * ``bidirectional`` -- the experiment ran both ways and the record covers
      both; it is not a midpoint.
    * ``unclear`` -- not derivable. Never a default; see
      :func:`syn4rch.schema.rules.perturbation_polarity_from`.
    """

    increase = "increase"
    decrease = "decrease"
    remove = "remove"
    bidirectional = "bidirectional"
    unclear = "unclear"


class Acuity(_StrEnum):
    acute = "acute"
    subchronic = "subchronic"
    chronic = "chronic"
    not_reported = "not_reported"


class Laterality(_StrEnum):
    unilateral = "unilateral"
    bilateral = "bilateral"
    not_reported = "not_reported"


class TimeUnit(_StrEnum):
    """Unit for ``duration_*`` and ``time_to_readout_*``.

    SCHEMA.md leaves these unenumerated; see docs/DECISIONS.md (D-006).
    """

    seconds = "seconds"
    minutes = "minutes"
    hours = "hours"
    days = "days"
    weeks = "weeks"
    months = "months"


# --- 3.5 readout and effect ------------------------------------------------

class OutcomeTier(_StrEnum):
    A = "A"
    B = "B"


class Direction(_StrEnum):
    increase = "increase"
    decrease = "decrease"
    no_change = "no_change"
    non_monotonic = "non_monotonic"
    not_reported = "not_reported"


class EffectSizeType(_StrEnum):
    percent_change = "percent_change"
    absolute_difference = "absolute_difference"
    fold_change = "fold_change"
    cohens_d = "cohens_d"
    auc_difference = "auc_difference"


class DispersionType(_StrEnum):
    sem = "sem"
    sd = "sd"
    ci95 = "ci95"


# --- 3.6 controls and evidence grading ------------------------------------

class ControlType(_StrEnum):
    opsin_free = "opsin_free"
    reporter_only = "reporter_only"
    vehicle = "vehicle"
    ligand_in_nonexpressing = "ligand_in_nonexpressing"
    light_off = "light_off"
    contralateral = "contralateral"
    cre_negative = "cre_negative"
    wildtype = "wildtype"
    none_reported = "none_reported"


class EvidenceClass(_StrEnum):
    causal_interventional = "causal_interventional"
    functional_association = "functional_association"
    anatomical = "anatomical"


class Confidence(_StrEnum):
    high = "high"
    moderate = "moderate"
    low = "low"


#: Control types that make a chemo-/optogenetic causal claim interpretable.
#: Used by :func:`syn4rch.schema.rules.has_actuator_control_from_controls`.
ACTUATOR_CONTROL_TYPES = frozenset(
    {
        ControlType.opsin_free,
        ControlType.reporter_only,
        ControlType.ligand_in_nonexpressing,
        ControlType.cre_negative,
        ControlType.light_off,
    }
)

#: Directions that assert an observed effect (as opposed to "not reported").
ASSERTIVE_DIRECTIONS = frozenset(
    {
        Direction.increase,
        Direction.decrease,
        Direction.no_change,
        Direction.non_monotonic,
    }
)
