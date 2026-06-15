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
