import pandas as pd

from darkdna.reports.region_cards import make_region_cards


def test_region_card_contains_assay_blueprint_fields():
    windows = pd.DataFrame(
        {
            "region_id": ["r1"],
            "chrom": ["scaffold_A"],
            "start": [10],
            "end": [210],
            "window_size": [200],
            "parent_region_id": [None],
            "child_region_ids": [""],
            "artifact_risk_flags": [""],
        }
    )
    labels = pd.DataFrame({"region_id": ["r1"], "primitive_class": ["negative_space_element_candidate"], "primitive_confidence": [0.8], "top_supporting_features": ["depleted_kmer_score"]})
    residuals = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive": ["negative_space_element_candidate_score"],
            "observed_score": [2.0],
            "residual_zscore": [3.0],
            "matched_null_zscore": [2.5],
            "empirical_p_value": [0.01],
            "classical_explanation_fraction": [0.2],
            "covariates_used": ["gc_content"],
        }
    )
    cards = make_region_cards(windows, labels, residuals)
    assert cards[0]["key_interaction_test"].startswith("negative_space_rescue_effect")
    assert cards[0]["candidate_only"] is True
    assert cards[0]["forbidden_interpretation"]
    assert cards[0]["recommended_primitive_assay"] == "Negative-Space Rescue/Scramble Assay"


def test_region_card_uses_primitive_specific_key_tests():
    windows = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "chrom": ["scaffold_A", "scaffold_A"],
            "start": [10, 300],
            "end": [210, 500],
            "window_size": [200, 200],
            "parent_region_id": [None, None],
            "child_region_ids": ["", ""],
            "artifact_risk_flags": ["", ""],
        }
    )
    labels = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "primitive_class": ["fractal_scaffold_candidate", "quantum_susceptible_domain_candidate"],
            "primitive_confidence": [0.8, 0.8],
            "top_supporting_features": ["fractal_score", "G4_susceptibility_proxy"],
        }
    )
    residuals = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "primitive": ["fractal_scaffold_candidate_score", "quantum_susceptible_domain_candidate_score"],
            "observed_score": [2.0, 2.0],
            "residual_zscore": [3.0, 3.0],
            "matched_null_zscore": [2.5, 2.5],
            "empirical_p_value": [0.01, 0.01],
            "classical_explanation_fraction": [0.2, 0.2],
            "covariates_used": ["gc_content", "gc_content"],
        }
    )
    cards = make_region_cards(windows, labels, residuals)
    tests_by_primitive = {card["primitive_class"]: card["key_interaction_test"] for card in cards}
    assert tests_by_primitive["fractal_scaffold_candidate"].startswith("folding_scale_effect")
    assert tests_by_primitive["quantum_susceptible_domain_candidate"].startswith("physical_susceptibility_effect")
    assert len(set(tests_by_primitive.values())) == 2


def test_no_call_and_unexplained_do_not_use_generic_key_test():
    windows = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "chrom": ["scaffold_A", "scaffold_A"],
            "start": [10, 300],
            "end": [210, 500],
            "window_size": [200, 200],
            "parent_region_id": [None, None],
            "child_region_ids": ["", ""],
            "artifact_risk_flags": ["", ""],
        }
    )
    labels = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "primitive_class": ["no_call", "unexplained_dark_anomaly_candidate"],
            "primitive_confidence": [0.0, 0.8],
            "top_supporting_features": ["", ""],
        }
    )
    residuals = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "primitive": ["unexplained_dark_anomaly_candidate_score", "unexplained_dark_anomaly_candidate_score"],
            "observed_score": [0.1, 2.0],
            "residual_zscore": [0.1, 3.0],
            "matched_null_zscore": [0.1, 2.5],
            "empirical_p_value": [0.9, 0.01],
            "classical_explanation_fraction": [0.8, 0.2],
            "covariates_used": ["gc_content", "gc_content"],
        }
    )
    cards = make_region_cards(windows, labels, residuals)
    tests_by_primitive = {card["primitive_class"]: card["key_interaction_test"] for card in cards}
    assert tests_by_primitive["no_call"].startswith("no_call_review")
    assert tests_by_primitive["unexplained_dark_anomaly_candidate"].startswith("dark_anomaly_effect")
    assert "Native_treatment" not in " ".join(tests_by_primitive.values())
