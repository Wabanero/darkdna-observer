import pandas as pd

from darkdna.primitives.ontology import get_primitive, primitive_names
from darkdna.views.primitive_scores import primitive_score_manifest, score_primitives


def test_prompt1_outputs_architecture_classes_not_mode_e_identities():
    names = primitive_names()
    assert "asymmetric_repeat_architecture_candidate" in names
    assert "periodic_spacing_grammar_candidate" in names
    assert "sequence_regime_boundary_candidate" in names
    assert "hysteresis_candidate" not in names
    assert "resonant_pulse_decoder_candidate" not in names
    assert "possibility_gate_candidate" not in names
    assert "criticality_tuner_candidate" not in names
    assert "chromatin_motion_oscillator_candidate" not in names
    assert "decoherence_boundary_candidate" not in names
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
    assert "asymmetric_repeat_architecture_candidate_score" in scores.columns
    assert "periodic_spacing_grammar_candidate_score" in scores.columns
    assert "hysteresis_candidate_score" not in scores.columns
    assert "possibility_gate_candidate_score" not in scores.columns
    assert "criticality_tuner_candidate_score" not in scores.columns
    assert scores.iloc[0]["asymmetric_repeat_architecture_candidate_score_weighting_scheme"] == "covariance_aware_robust_cohort_standardization"
    assert "not_null_calibrated" in scores.iloc[0]["asymmetric_repeat_architecture_candidate_score_calibration_status"]
    assert scores.iloc[0]["asymmetric_repeat_architecture_candidate_score_empirical_p_value_status"] == "unavailable"


def test_primitive_score_manifest_marks_composites_as_screening_views():
    manifest = primitive_score_manifest()
    assert manifest["score_status"] == "covariance_aware_cohort_standardized_screening_score"
    assert manifest["interpretation_order"] == [
        "measured_feature_profile",
        "statistical_anomaly_after_explicit_controls_and_nulls",
        "post_hoc_mechanistic_hypothesis",
    ]
    assert "asymmetric_repeat_architecture_candidate_score" in manifest["components"]
    assert "hysteresis_candidate_score" not in manifest["components"]
    assert "correlation_and_double_counting_audit" in manifest["required_validation_before_mechanistic_use"]
    assert manifest["empirical_p_value_policy"] == "unavailable_without_explicit_null_distribution"


def test_ontology_distinguishes_candidate_and_confirmed_names():
    architecture = get_primitive("hysteresis_candidate")
    assert architecture.candidate_name == "asymmetric_repeat_architecture_candidate"
    assert architecture.confirmed_name == "asymmetric_repeat_architecture"
    assert architecture.allowed_in_prompt1 is True
    assert architecture.requires_dynamic_data is False
    assert "hysteresis" in architecture.forbidden_interpretation

    physical = get_primitive("non_B_DNA_physical_susceptibility_candidate")
    assert physical.confirmed_name == "physical_susceptibility_domain"
    assert physical.requires_dynamic_data is False
    assert "Do not claim quantum susceptibility" in physical.forbidden_interpretation

    legacy = get_primitive("quantum_susceptible_domain_candidate")
    assert legacy.candidate_name == "non_B_DNA_physical_susceptibility_candidate"

    boundary = get_primitive("decoherence_boundary_candidate")
    assert boundary.candidate_name == "sequence_regime_boundary_candidate"
    periodic = get_primitive("resonant_pulse_decoder_candidate")
    assert periodic.candidate_name == "periodic_spacing_grammar_candidate"
