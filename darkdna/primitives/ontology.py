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
    "quantum_susceptible_domain": PrimitiveClass(
        "quantum_susceptible_domain",
        "quantum_susceptible_domain_candidate",
        "physical_susceptibility_domain",
        "sequence_plus_optional_physical_assay",
        True,
        False,
        "Sequence-derived physical-susceptibility proxy based on G-richness, G4 potential, oxidation-prone contexts, G-skewness, non-B-DNA propensity and charge-transfer-like sequence contexts.",
        "Sequence-derived physical-susceptibility proxy based on G-richness, G4 potential, oxidation-prone contexts, G-skewness, non-B-DNA propensity and charge-transfer-like sequence contexts.",
        "Do not claim actual quantum effects.",
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
    "chromatin_motion_oscillator": PrimitiveClass(
        "chromatin_motion_oscillator",
        "chromatin_motion_oscillator_candidate",
        "chromatin_motion_oscillator",
        "live_locus_motion_or_spatial_dynamics",
        True,
        True,
        "Static sequence architecture that motivates locus-motion or spatial-dynamics testing.",
        "Sequence architecture compatible with possible spatial-dynamics or motion-response hypotheses.",
        "Do not claim oscillatory chromatin behavior from static sequence.",
        "live-locus motion / spatial dynamics view",
    ),
    "decoherence_boundary": PrimitiveClass(
        "decoherence_boundary",
        "decoherence_boundary_candidate",
        "decoherence_boundary",
        "single_cell_variance_or_noise_propagation_assay",
        True,
        True,
        "Entropy/noise or sequence-regime boundary proxy that motivates noise-propagation tests.",
        "Sequence architecture compatible with possible variance/noise boundary behavior.",
        "Do not claim physical decoherence or confirmed noise insulation from static sequence.",
        "single-cell variance/noise propagation view",
    ),
    "hysteresis": PrimitiveClass(
        "hysteresis",
        "hysteresis_candidate",
        "hysteresis_element",
        "perturbation_timecourse_or_recovery",
        True,
        True,
        "Asymmetry, repeat architecture, and non-B propensity compatible with history-dependent perturbation tests.",
        "Sequence architecture compatible with possible memory-like or metastable behavior.",
        "Do not claim actual hysteresis or memory from static sequence.",
        "timecourse_recovery / hysteresis view",
    ),
    "resonant_pulse_decoder": PrimitiveClass(
        "resonant_pulse_decoder",
        "resonant_pulse_decoder_candidate",
        "resonant_pulse_decoder",
        "temporal_pulse_assay_or_timecourse",
        True,
        True,
        "Periodic, phase, or spacing grammar compatible with temporal pulse perturbation tests.",
        "Periodic/phase/spacing grammar compatible with pulse-decoding hypotheses.",
        "Do not claim frequency decoding from sequence alone.",
        "temporal pulse / timecourse view",
    ),
    "possibility_gate": PrimitiveClass(
        "possibility_gate",
        "possibility_gate_candidate",
        "possibility_gate",
        "state_transition_graph_or_pseudotime",
        True,
        True,
        "Boundary or constraint-like sequence architecture compatible with future state-transition tests.",
        "Sequence architecture compatible with boundary/constraint-like behavior.",
        "Do not claim future-state bias or reachable-state modulation from static sequence.",
        "state_transition / constructor view",
    ),
    "criticality_tuner": PrimitiveClass(
        "criticality_tuner",
        "criticality_tuner_candidate",
        "criticality_tuner",
        "dose_gradient_timecourse_or_pseudotime",
        True,
        True,
        "Sequence-regime boundary or entropy/compression transition compatible with threshold-like hypotheses.",
        "Sequence-regime boundary or entropy/compression transition compatible with threshold-like hypotheses.",
        "Do not claim biological criticality or transition threshold shift from sequence alone.",
        "criticality / dose-gradient view",
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


_ALIASES: dict[str, str] = {}
for key, primitive in PRIMITIVES.items():
    _ALIASES[key] = key
    _ALIASES[primitive.candidate_name] = key
    _ALIASES[primitive.confirmed_name] = key


def primitive_names() -> list[str]:
    return [primitive.candidate_name for primitive in PRIMITIVES.values()]


def get_primitive(name: str) -> PrimitiveClass:
    key = _ALIASES.get(name, "unexplained_dark_anomaly")
    return PRIMITIVES[key]
