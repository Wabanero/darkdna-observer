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
    assert {"observed_score", "predicted_classical_score", "residual_zscore", "classical_model_global_r2", "classical_explanation_fraction"}.issubset(residuals.columns)
    assert residuals["empirical_p_value"].isna().all()
    assert set(residuals["empirical_p_value_status"]) == {"unavailable"}
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


def test_residualization_uses_blocked_cv_groups_without_covariate_leakage():
    scores = pd.DataFrame(
        {
            "region_id": [f"r{i}" for i in range(6)],
            "fractal_scaffold_candidate_score": [1.0, 1.2, 2.0, 2.2, 3.0, 3.2],
        }
    )
    covariates = pd.DataFrame(
        {
            "region_id": [f"r{i}" for i in range(6)],
            "block_id": ["b1", "b1", "b2", "b2", "b3", "b3"],
            "gc_content": [0.40, 0.42, 0.50, 0.52, 0.60, 0.62],
            "start": [0, 500, 100000, 100500, 200000, 200500],
            "end": [1000, 1500, 101000, 101500, 201000, 201500],
        }
    )

    residuals, summary = residualize_scores(scores, covariates, method="linear")

    assert {"cv_strategy", "cv_group_col", "cv_group_id"}.issubset(residuals.columns)
    assert set(residuals["cv_group_col"]) == {"block_id"}
    fractal_summary = summary[summary["primitive"] == "fractal_scaffold_candidate_score"].iloc[0]
    assert fractal_summary["cv_strategy"] == "leave_one_block_id_out"
    assert "block_id" not in fractal_summary["covariates_used"]
    assert "start" not in fractal_summary["covariates_used"]
    assert "end" not in fractal_summary["covariates_used"]
