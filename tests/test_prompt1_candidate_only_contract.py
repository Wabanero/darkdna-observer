import pandas as pd

from darkdna.primitives.ontology import primitive_names
from darkdna.views.primitive_scores import score_primitives


def test_prompt1_candidate_only_labels():
    forbidden_confirmed = {
        "hysteresis_element",
        "resonant_pulse_decoder",
        "possibility_gate",
        "criticality_tuner",
        "attractor_gate",
        "future_state_biaser",
        "trajectory_constraint",
        "active_inference_prior_region",
        "genomic_reservoir",
        "contextual_operator",
    }
    allowed_candidates = {
        "hysteresis_candidate",
        "resonant_pulse_decoder_candidate",
        "possibility_gate_candidate",
        "criticality_tuner_candidate",
        "chromatin_motion_oscillator_candidate",
        "replication_instability_candidate",
        "non_B_DNA_physical_susceptibility_candidate",
        "sequence_regime_boundary_candidate",
        "negative_space_element_candidate",
        "TE_grammar_node_candidate",
        "fractal_scaffold_candidate",
        "constraint_grammar_region_candidate",
        "decoherence_boundary_candidate",
        "unexplained_dark_anomaly_candidate",
    }
    names = set(primitive_names())
    assert forbidden_confirmed.isdisjoint(names)
    assert allowed_candidates.issubset(names)
    assert "quantum_susceptible_domain_candidate" not in names


def test_sequence_only_commands_do_not_emit_prompt2_dynamic_scores():
    scores = score_primitives(pd.DataFrame({"region_id": ["r1"], "fractal_score": [1.0]}))
    forbidden_scores = {
        "future_state_prediction_gain",
        "transition_bias_score",
        "reachable_state_delta",
        "state_space_coverage_change",
        "path_dependence_score",
        "transition_threshold_shift_score",
        "state_stabilization_score",
        "active_inference_prior_score",
        "teleology_operational_score",
    }
    assert forbidden_scores.isdisjoint(scores.columns)
