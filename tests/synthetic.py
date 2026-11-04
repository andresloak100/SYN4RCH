"""Builders for FABRICATED test records.

Everything produced here is invented for the purpose of exercising code. No
value corresponds to any real paper, author, population, measurement or
finding. Every ``study_id`` carries the ``synthetic:`` prefix precisely so
that a synthetic record can never be mistaken for evidence, in a file, in a
report or in a grep.

These helpers must never be imported by anything under ``src/``.
"""

import _path  # noqa: F401  (sys.path bootstrap)

from syn4rch.schema import TierARecord
from syn4rch.schema.rules import (
    perturbation_class_from_method,
    perturbation_polarity_from,
)

#: Prefix every fabricated study identifier carries.
SYNTHETIC_PREFIX = "synthetic:"


def good_record(**overrides) -> TierARecord:
    """A fabricated record that satisfies every rule, hard and soft.

    ``validate(good_record())`` returns an empty list: no hard findings and no
    warnings. Tests mutate one field at a time from this baseline so that a
    single assertion pins a single rule.
    """
    data = dict(
        study_id="synthetic:test-001",
        doi="synthetic:test-001",
        pmid=None,
        title="Synthetic test record for schema validation",
        journal="Journal of Synthetic Test Data",
        year=2001,
        first_author="Testcase, A.",
        last_author="Fixture, B.",
        last_author_affiliation="Institute of Synthetic Data, Nowhere",
        retrieved_from="manual",
        full_text_available=True,
        extraction_method="llm_human_reviewed",
        extracted_by="unit-test",
        extracted_at="2001-01-01T00:00:00Z",
        review_status="accepted",
        population_raw="synthetic population alpha",
        population_id="syn:pop_alpha",
        marker_genes=["Synth1"],
        region_raw="synthetic region",
        region_id="syn:region_1",
        neurotransmitter=["unknown"],
        projection_specific=False,
        species="mus_musculus",
        strain="synthetic strain",
        sex="male",
        age_value=12.0,
        age_unit="weeks",
        age_group="young_adult",
        diet="chow",
        n_animals=16,
        # condition dimension, schema v0.2.0. Fully specified on purpose, so
        # that condition_specificity(good_record()) == 4 and a test that cares
        # about a *missing* condition has to say so explicitly.
        challenge="ad_libitum",
        challenge_detail="fed ad libitum (fabricated)",
        condition_source="passage",
        condition_locator=None,
        method="chemogenetic",
        actuator="synthetic-actuator",
        valence="activation",
        # circuit_activity + increase, which is what the rule derives from
        # (chemogenetic, activation). A chemogenetic record has no molecular
        # target, so target_raw/target_id stay null.
        perturbation_class="circuit_activity",
        perturbation_polarity="increase",
        target_raw=None,
        target_id=None,
        targeting_strategy="synthetic Cre line + synthetic vector",
        acuity="acute",
        ligand="synthetic ligand",
        ligand_dose_value=1.0,
        ligand_dose_unit="mg/kg",
        ligand_route="ip",
        laterality="bilateral",
        time_to_readout_value=1.0,
        time_to_readout_unit="hours",
        outcome_raw="synthetic readout one",
        outcome_id="syn:outcome_1",
        outcome_tier="A",
        tissue_raw="synthetic tissue",
        assay="synthetic assay",
        assay_units="syn_unit",
        direction="increase",
        effect_size_value=25.0,
        effect_size_type="percent_change",
        effect_size_units="%",
        dispersion_type="sem",
        dispersion_value=3.0,
        n_treatment=8,
        n_control=8,
        statistical_test="synthetic test",
        p_value=0.01,
        significant=True,
        is_primary_outcome=True,
        control_types=["vehicle", "ligand_in_nonexpressing"],
        has_actuator_control=True,
        evidence_class="causal_interventional",
        confidence="high",
        support_quote=(
            "Synthetic quote: the readout increased by 25% in the treated "
            "group (fabricated text, not from any publication)."
        ),
        source_locator="Fig S1A (synthetic)",
        extraction_notes="Fabricated fixture; never evidence.",
    )
    data.update(overrides)
    data.pop("record_id", None)
    # ``perturbation_class`` and ``perturbation_polarity`` are derived fields
    # (schema v0.2.0). Re-derive them from whatever method and valence the
    # caller ended up with, unless the caller set them deliberately. Without
    # this, ``good_record(valence="inhibition")`` would carry a stale
    # ``perturbation_polarity`` of "increase": the record would trip W013, and
    # every test that changes a valence to make a *different claim* would
    # silently keep the original claim key.
    if "perturbation_class" not in overrides:
        derived = perturbation_class_from_method(data.get("method"))
        data["perturbation_class"] = derived
    if "perturbation_polarity" not in overrides:
        data["perturbation_polarity"] = perturbation_polarity_from(
            data.get("method"), data.get("valence")
        )
    return TierARecord(**data)


def molecular_record(**overrides) -> TierARecord:
    """A fabricated *molecular* perturbation: loss of function on a target gene.

    Exists because schema v0.2.0 made molecular and circuit perturbations two
    different things, and a fixture that is only ever chemogenetic cannot
    exercise the difference. ``valence`` stays ``unclear`` here and that is
    correct, not an omission: a gene knockout's net effect on firing is
    genuinely unknown, and SCHEMA.md 3.4 forbids coercing it.
    """
    data = dict(
        study_id="synthetic:test-mol-001",
        doi="synthetic:test-mol-001",
        title="Synthetic molecular perturbation record",
        method="conditional_knockout",
        valence="unclear",
        perturbation_class="molecular",
        perturbation_polarity="decrease",
        target_raw="Synthgene1",
        target_id="syn:synthgene1",
        actuator=None,
        ligand=None,
        ligand_dose_value=None,
        ligand_dose_unit=None,
        ligand_route=None,
        control_types=["cre_negative", "wildtype"],
        targeting_strategy="floxed Synthgene1 crossed to a synthetic Cre line",
        source_locator="Fig S4D (synthetic)",
    )
    data.update(overrides)
    return good_record(**data)


def other_lab_record(**overrides) -> TierARecord:
    """A fabricated record from a *different* synthetic study and lab."""
    data = dict(
        study_id="synthetic:test-002",
        doi="synthetic:test-002",
        title="Second synthetic test record",
        last_author="Other, C.",
        last_author_affiliation="Second Institute of Synthetic Data, Elsewhere",
        source_locator="Fig S2B (synthetic)",
    )
    data.update(overrides)
    return good_record(**data)


def same_lab_other_study_record(**overrides) -> TierARecord:
    """Different synthetic study, same synthetic lab (non-independence case)."""
    data = dict(
        study_id="synthetic:test-003",
        doi="synthetic:test-003",
        title="Third synthetic test record, same synthetic lab",
        source_locator="Fig S3C (synthetic)",
    )
    data.update(overrides)
    return good_record(**data)


def stub_vocab_loader(outcomes=None):
    """A stand-in for ``syn4rch.vocab.load`` for unit tests.

    The real vocabulary files are owned by another component and may not
    exist; the QC suite must not depend on them.
    """
    outcomes = (
        {"syn:outcome_1": {"label": "synthetic readout one", "units": "syn_unit"}}
        if outcomes is None
        else outcomes
    )

    def load(kind):
        if kind == "outcomes":
            return outcomes
        return {}

    return load
