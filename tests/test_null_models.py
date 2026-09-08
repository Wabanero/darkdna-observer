import pandas as pd

from darkdna.residuals.null_models import build_matched_null_models, null_panel_status


def test_matched_null_model_schema():
    scores = pd.DataFrame(
        {
            "region_id": ["r1", "r2", "r3"],
            "fractal_scaffold_candidate_score": [3.0, 1.0, 2.0],
            "constraint_grammar_region_candidate_score": [1.0, 2.0, 3.0],
            "non_B_DNA_physical_susceptibility_candidate_score": [0.1, 0.2, 0.3],
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
    assert {
        "null_model_id",
        "region_id",
        "primitive",
        "null_zscore",
        "empirical_p_value",
        "null_panel_status",
        "missing_or_partial_null_models",
    }.issubset(nulls.columns)
    assert not nulls.empty
    assert "dinucleotide_preserving_shuffle" in nulls.iloc[0]["missing_or_partial_null_models"]
    assert nulls["null_panel_status"].eq("partial_null_panel_not_for_promotion").all()
    assert (nulls["sequence_null_model_count"] == 0).all()


def test_null_panel_status_is_not_single_zscore_sufficient():
    status = null_panel_status()
    assert status["status"] == "insufficient_single_matched_null_until_complementary_nulls_pass"
    assert "matched_controls_v1" in status["implemented_null_models"]
    assert "syntenic_ortholog_controls" in status["missing_or_partial_null_models"]
