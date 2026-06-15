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
    assert cards[0]["key_interaction_test"].startswith("effect =")
    assert cards[0]["candidate_only"] is True
    assert cards[0]["forbidden_interpretation"]
    assert cards[0]["recommended_primitive_assay"] == "Negative-Space Rescue/Scramble Assay"
