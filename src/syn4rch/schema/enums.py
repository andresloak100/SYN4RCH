"""Controlled enumerations for the SYN4RCH Tier A record schema.

Every enum here mirrors a vocabulary fixed by ``docs/SCHEMA.md`` (v0.1.0).
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
    "PerturbationMethod",
    "Valence",
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
    activation = "activation"
    inhibition = "inhibition"
    ablation = "ablation"
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
