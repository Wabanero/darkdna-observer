"""Operational primitive ontology for Prompt 1.

Prompt 1 is sequence-first. It emits candidate proxy labels only, never
confirmed temporal, dynamical, teleological, or active-inference primitives.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PrimitiveClass:
    ontology_key: str
    candidate_name: str
    confirmed_name: str
    required_input_level: str
    allowed_in_prompt1: bool
    requires_dynamic_data: bool
    description: str
    prompt1_interpretation: str
    forbidden_interpretation: str
    suggested_prompt2_view: str


PRIMITIVES: dict[str, PrimitiveClass] = {
    "fractal_scaffold": PrimitiveClass(
        "fractal_scaffold",
        "fractal_scaffold_candidate",
        "fractal_scaffold",
        "sequence_plus_matched_nulls_and_optional_physical_assay",
        True,
        False,
        "Scale-persistent sequence texture, numeric-walk structure, or compression/fractal-like summaries.",
        "Sequence architecture compatible with scale-constraint or folding-substrate hypotheses.",
        "Do not interpret as proof of holographic genome physics or true mathematical self-similarity.",
        "multiscale perturbation / folding validation view",
    ),
    "constraint_grammar_region": PrimitiveClass(
        "constraint_grammar_region",
        "constraint_grammar_region_candidate",
        "constraint_grammar_region",
        "sequence_plus_matched_nulls_and_sequence_scramble_assay",
        True,
        False,
        "Sequence grammar, spacing, transition surprise, or token recurrence that remains anomalous after controls.",
        "Sequence architecture compatible with grammar-constraint hypotheses.",
        "Do not interpret as a decoded semantic language or teleological program.",
        "grammar scramble / tiling perturbation view",
    ),
    "non_B_DNA_physical_susceptibility": PrimitiveClass(
        "non_B_DNA_physical_susceptibility",
        "non_B_DNA_physical_susceptibility_candidate",
        "physical_susceptibility_domain",
        "sequence_plus_optional_physical_assay",
        True,
        False,
        "Sequence-derived physical-susceptibility proxy based on G-richness, G4 potential, oxidation-prone contexts, G-skewness, non-B-DNA propensity and charge-transfer-like sequence contexts.",
        "Sequence-derived physical-susceptibility proxy based on G-richness, G4 potential, oxidation-prone contexts, G-skewness, non-B-DNA propensity and charge-transfer-like sequence contexts.",
        "Do not claim quantum susceptibility, actual quantum effects, or charge-transfer dynamics from sequence proxies alone.",
        "physical validation / charge-oxidation susceptibility view",
    ),
    "replication_instability": PrimitiveClass(
        "replication_instability",
        "replication_instability_candidate",
        "replication_instability_region",
        "replication_stress_or_fork_assay",
        True,
        False,
        "Repeat, palindrome, tract, skew, or fork-texture proxies suggesting replication-stress susceptibility.",
        "Sequence architecture compatible with replication-instability hypotheses.",
        "Do not claim replication instability without replication assays.",
        "replication fork texture / stress assay view",
    ),
    "periodic_spacing_grammar": PrimitiveClass(
        "periodic_spacing_grammar",
        "periodic_spacing_grammar_candidate",
        "periodic_spacing_grammar",
        "sequence_plus_matched_nulls_and_phase_scramble_assay",
        True,
        False,
        "Helical, nucleosome-scale, or spacing periodicity that remains anomalous after composition controls.",
        "Sequence architecture compatible with periodic spacing or phasing hypotheses.",
        "Do not claim pulse decoding, frequency response, or chromatin motion from static sequence.",
        "temporal pulse / live-locus motion view after an intermediate molecular bridge",
    ),
    "asymmetric_repeat_architecture": PrimitiveClass(
        "asymmetric_repeat_architecture",
        "asymmetric_repeat_architecture_candidate",
        "asymmetric_repeat_architecture",
        "sequence_plus_matched_nulls_and_orientation_or_repeat_assay",
        True,
        False,
        "Left/right composition asymmetry, nested repeats, or recurrent-k-mer orientation remaining after controls.",
        "Sequence architecture compatible with oriented or nested repeat hypotheses.",
        "Do not claim hysteresis, memory, or history-dependent response from static sequence.",
        "timecourse recovery / hysteresis view after an intermediate molecular bridge",
    ),
    "negative_space_element": PrimitiveClass(
        "negative_space_element",
        "negative_space_element_candidate",
        "negative_space_element",
        "sequence_plus_matched_nulls_and_rescue_or_scramble_assay",
        True,
        False,
        "Structured absence, depleted words, deserts, or feature voids that remain anomalous after controls.",
        "Sequence-derived negative-space substrate candidate.",
        "Do not interpret absence as function without matched-null and perturbation evidence.",
        "negative-space rescue / scramble assay view",
    ),
    "sequence_regime_boundary": PrimitiveClass(
        "sequence_regime_boundary",
        "sequence_regime_boundary_candidate",
        "sequence_regime_boundary",
        "sequence_plus_boundary_validation_assay",
        True,
        False,
        "Candidate boundary between intrinsic sequence regimes.",
        "Sequence-derived regime-boundary candidate.",
        "Do not claim enhancer, insulator, or chromatin boundary function without validation.",
        "boundary disruption / insulation validation view",
    ),
    "TE_grammar_node": PrimitiveClass(
        "TE_grammar_node",
        "TE_grammar_node_candidate",
        "TE_grammar_node",
        "sequence_plus_TE_annotation_and_perturbation_assay",
        True,
        False,
        "TE-derived or TE-mosaic architecture that remains anomalous after controlling for simple TE overlap.",
        "Sequence and annotation-derived TE grammar candidate.",
        "Do not claim TE exaptation or regulatory function from overlap alone.",
        "TE grammar reconstruction / perturbation view",
    ),
    "unexplained_dark_anomaly": PrimitiveClass(
        "unexplained_dark_anomaly",
        "unexplained_dark_anomaly_candidate",
        "unexplained_dark_anomaly",
        "sequence_plus_matched_nulls_and_artifact_review",
        True,
        False,
        "High residual sequence anomaly without one dominant candidate-specific explanation.",
        "Sequence-derived residual anomaly candidate requiring controls and validation.",
        "Do not over-interpret; prioritize controls and artifact review.",
        "matched-null review / orthogonal validation view",
    ),
}


# Retired Mode E identity labels. Static sequence can motivate a later assay
# hypothesis, but these names are not architecture classes and must not appear
# as primary primitive labels.
LEGACY_PRIMITIVE_ALIASES: dict[str, str] = {
    "periodic_spacing_grammar_candidate_score": "periodic_spacing_grammar",
    "asymmetric_repeat_architecture_candidate_score": "asymmetric_repeat_architecture",
    "hysteresis": "asymmetric_repeat_architecture",
    "hysteresis_candidate": "asymmetric_repeat_architecture",
    "hysteresis_element": "asymmetric_repeat_architecture",
    "hysteresis_candidate_score": "asymmetric_repeat_architecture",
    "resonant_pulse_decoder": "periodic_spacing_grammar",
    "resonant_pulse_decoder_candidate": "periodic_spacing_grammar",
    "resonant_pulse_decoder_candidate_score": "periodic_spacing_grammar",
    "chromatin_motion_oscillator": "periodic_spacing_grammar",
    "chromatin_motion_oscillator_candidate": "periodic_spacing_grammar",
    "chromatin_motion_oscillator_candidate_score": "periodic_spacing_grammar",
    "decoherence_boundary": "sequence_regime_boundary",
    "decoherence_boundary_candidate": "sequence_regime_boundary",
    "decoherence_boundary_candidate_score": "sequence_regime_boundary",
    "possibility_gate": "sequence_regime_boundary",
    "possibility_gate_candidate": "sequence_regime_boundary",
    "possibility_gate_candidate_score": "sequence_regime_boundary",
    "criticality_tuner": "sequence_regime_boundary",
    "criticality_tuner_candidate": "sequence_regime_boundary",
    "criticality_tuner_candidate_score": "sequence_regime_boundary",
    "quantum_susceptible_domain": "non_B_DNA_physical_susceptibility",
    "quantum_susceptible_domain_candidate": "non_B_DNA_physical_susceptibility",
    "quantum_susceptible_domain_candidate_score": "non_B_DNA_physical_susceptibility",
}

MODE_E_IDENTITY_LABELS = frozenset(
    {
        "hysteresis_candidate",
        "resonant_pulse_decoder_candidate",
        "possibility_gate_candidate",
        "criticality_tuner_candidate",
        "chromatin_motion_oscillator_candidate",
        "decoherence_boundary_candidate",
    }
)


_ALIASES: dict[str, str] = {}
for key, primitive in PRIMITIVES.items():
    _ALIASES[key] = key
    _ALIASES[primitive.candidate_name] = key
    _ALIASES[primitive.confirmed_name] = key
    _ALIASES[f"{primitive.candidate_name}_score"] = key
_ALIASES.update(LEGACY_PRIMITIVE_ALIASES)


def primitive_names() -> list[str]:
    return [primitive.candidate_name for primitive in PRIMITIVES.values()]


def get_primitive(name: str) -> PrimitiveClass:
    key = _ALIASES.get(name, "unexplained_dark_anomaly")
    return PRIMITIVES[key]
