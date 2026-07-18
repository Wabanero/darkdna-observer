import pandas as pd

from darkdna.primitives.labeler import assign_primitive_labels
from darkdna.views.primitive_scores import score_primitives


def test_primitive_score_generation_and_assignment():
    features = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "fractal_score": [5.0, 0.1],
            "scale_persistence_score": [1.0, 0.0],
            "compression_anomaly_score": [1.0, 0.0],
        }
    )
    scores = score_primitives(features)
    assert "fractal_scaffold_candidate_score" in scores.columns
    assert "fractal_scaffold_score" not in scores.columns
    residuals = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive": ["fractal_scaffold_candidate_score"],
            "observed_score": [5.0],
            "predicted_classical_score": [0.0],
            "residual_score": [5.0],
            "residual_zscore": [3.0],
            "matched_null_zscore": [3.0],
            "empirical_p_value": [0.01],
            "classical_explanation_fraction": [0.1],
            "covariates_used": ["gc_content"],
        }
    )
    labels = assign_primitive_labels(residuals, features=features)
    assert labels.iloc[0]["primitive_class"] == "fractal_scaffold_candidate"
    assert labels.iloc[0]["candidate_promotion_status"] == "screening_only_legacy_null_metadata_unavailable"


def test_primitive_promotion_requires_agreement_across_severe_null_panel():
    residuals = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "primitive": ["fractal_scaffold_candidate_score"] * 2,
            "residual_zscore": [3.0, 3.0],
            "matched_null_zscore": [3.0, 3.0],
            "null_panel_status": ["severe_null_panel_available"] * 2,
            "null_model_count": [4, 4],
            "null_model_agreement": [1.0, 0.5],
            "null_model_conflict": [False, True],
        }
    )

    labels = assign_primitive_labels(residuals)

    status = labels.set_index("region_id")["candidate_promotion_status"].to_dict()
    assert status["r1"] == "eligible_for_candidate_promotion"
    assert status["r2"] == "screening_only_conflicting_null_families"
    assert labels.set_index("region_id").loc["r1", "survives_severe_null_panel"]
    assert not labels.set_index("region_id").loc["r2", "survives_severe_null_panel"]


def test_primitive_labeler_handles_unavailable_null_count():
    residuals = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive": ["fractal_scaffold_candidate_score"],
            "residual_zscore": [3.0],
            "matched_null_zscore": [float("nan")],
            "null_model_count": [float("nan")],
        }
    )

    label = assign_primitive_labels(residuals).iloc[0]

    assert label["null_model_count"] == 0
    assert label["candidate_promotion_status"] == "screening_only_legacy_null_metadata_unavailable"
