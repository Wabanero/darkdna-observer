import pandas as pd

from darkdna.residuals.null_models import build_matched_null_models


def test_matched_null_model_schema():
    scores = pd.DataFrame(
        {
            "region_id": ["r1", "r2", "r3"],
            "fractal_scaffold_candidate_score": [3.0, 1.0, 2.0],
            "constraint_grammar_region_candidate_score": [1.0, 2.0, 3.0],
            "quantum_susceptible_domain_candidate_score": [0.1, 0.2, 0.3],
            "replication_instability_candidate_score": [0.1, 0.2, 0.3],
            "decoherence_boundary_candidate_score": [0.1, 0.2, 0.3],
            "resonant_pulse_decoder_candidate_score": [0.1, 0.2, 0.3],
            "hysteresis_candidate_score": [0.1, 0.2, 0.3],
            "possibility_gate_candidate_score": [0.1, 0.2, 0.3],
            "criticality_tuner_candidate_score": [0.1, 0.2, 0.3],
            "chromatin_motion_oscillator_candidate_score": [0.1, 0.2, 0.3],
            "negative_space_element_candidate_score": [0.1, 0.2, 0.3],
            "sequence_regime_boundary_candidate_score": [0.1, 0.2, 0.3],
            "TE_grammar_node_candidate_score": [0.1, 0.2, 0.3],
            "unexplained_dark_anomaly_candidate_score": [0.5, 0.6, 0.7],
        }
    )
    features = pd.DataFrame({"region_id": ["r1", "r2", "r3"], "gc_content": [0.4, 0.5, 0.6], "length": [100, 100, 100]})
    nulls = build_matched_null_models(scores, features, n_controls=2)
    assert {"null_model_id", "region_id", "primitive", "null_zscore", "empirical_p_value"}.issubset(nulls.columns)
    assert not nulls.empty
