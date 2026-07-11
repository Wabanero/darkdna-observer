import pandas as pd

from darkdna.residuals.classical_covariates import prepare_classical_covariates
from darkdna.residuals.residual_model import residualize_scores


def test_residualization_output_schema_and_leakage_guard():
    scores = pd.DataFrame(
        {
            "region_id": ["r1", "r2", "r3", "r4"],
            "fractal_scaffold_candidate_score": [1.0, 2.0, 3.0, 4.0],
            "unexplained_dark_anomaly_candidate_score": [1.0, 1.5, 2.5, 3.5],
        }
    )
    covariates = pd.DataFrame({"region_id": ["r1", "r2", "r3", "r4"], "gc_content": [0.4, 0.5, 0.6, 0.7], "fractal_scaffold_candidate_score": [99, 99, 99, 99]})
    residuals, summary = residualize_scores(scores, covariates, method="linear")
    assert {"observed_score", "predicted_classical_score", "residual_zscore", "classical_explanation_fraction"}.issubset(residuals.columns)
    fractal_summary = summary[summary["primitive"] == "fractal_scaffold_candidate_score"].iloc[0]
    assert "fractal_scaffold_candidate_score" not in fractal_summary["covariates_used"]
    assert "unexplained_dark_anomaly_candidate_score" not in fractal_summary["covariates_used"]
    assert fractal_summary["covariates_used"] == "gc_content"


def test_prepare_classical_covariates_excludes_anomaly_features():
    windows = pd.DataFrame(
        {
            "region_id": ["r1"],
            "chrom": ["chr1"],
            "gc_content": [0.42],
            "n_fraction": [0.0],
            "overlaps_promoter": [False],
            "mappability": [0.98],
            "window_size": [1000],
        }
    )
    features = pd.DataFrame(
        {
            "region_id": ["r1"],
            "simple_repeat_fraction": [0.12],
            "fractal_score": [0.9],
            "hysteresis_candidate_score": [0.8],
        }
    )
    covariates = prepare_classical_covariates(windows, features)
    assert {"gc_content", "simple_repeat_fraction", "mappability", "window_size"}.issubset(covariates.columns)
    assert "fractal_score" not in covariates.columns
    assert "hysteresis_candidate_score" not in covariates.columns
