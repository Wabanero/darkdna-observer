import pandas as pd

from darkdna.primitives.ontology import get_primitive, primitive_names
from darkdna.views.primitive_scores import score_primitives


def test_prompt1_outputs_candidate_labels_not_confirmed_dynamic_labels():
    names = primitive_names()
    assert "hysteresis_candidate" in names
    assert "resonant_pulse_decoder_candidate" in names
    assert "possibility_gate_candidate" in names
    assert "criticality_tuner_candidate" in names
    assert "hysteresis_element" not in names
    assert "resonant_pulse_decoder" not in names
    assert "possibility_gate" not in names
    assert "criticality_tuner" not in names


def test_sequence_only_scores_do_not_produce_prompt2_dynamic_scores():
    features = pd.DataFrame({"region_id": ["r1"], "fractal_score": [1.0], "boundary_condition_candidate_score": [0.5]})
    scores = score_primitives(features)
    forbidden = {
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
    assert forbidden.isdisjoint(scores.columns)
    assert "possibility_gate_candidate_score" in scores.columns
    assert "criticality_tuner_candidate_score" in scores.columns


def test_ontology_distinguishes_candidate_and_confirmed_names():
    hysteresis = get_primitive("hysteresis_candidate")
    assert hysteresis.candidate_name == "hysteresis_candidate"
    assert hysteresis.confirmed_name == "hysteresis_element"
    assert hysteresis.allowed_in_prompt1 is True
    assert hysteresis.requires_dynamic_data is True
    assert "static sequence" in hysteresis.forbidden_interpretation

    quantum = get_primitive("quantum_susceptible_domain_candidate")
    assert quantum.confirmed_name == "physical_susceptibility_domain"
    assert quantum.requires_dynamic_data is False
    assert "Do not claim actual quantum effects" in quantum.forbidden_interpretation
