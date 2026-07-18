import numpy as np
import pandas as pd

from darkdna.residuals.residual_model import residualize_scores


def test_residuals_include_blocked_conditional_quantile_and_conformal_calibration():
    x = np.linspace(0.0, 1.0, 30)
    scores = pd.DataFrame(
        {
            "region_id": [f"r{i}" for i in range(30)],
            "fractal_scaffold_candidate_score": 1.5 * x + np.sin(np.arange(30)) * (0.05 + x * 0.2),
        }
    )
    covariates = pd.DataFrame(
        {
            "region_id": scores["region_id"],
            "block_id": ["b1"] * 10 + ["b2"] * 10 + ["b3"] * 10,
            "gc_content": x,
        }
    )

    residuals, summary = residualize_scores(scores, covariates)

    assert set(residuals["conditional_variance_status"].str.split(":").str[0]) == {"available"}
    assert residuals["conditional_residual_zscore"].notna().all()
    assert residuals["quantile_residual"].notna().all()
    assert residuals["conformal_interval_lower"].notna().all()
    assert residuals["empirical_p_value"].isna().all()
    assert residuals["classical_model_global_r2"].nunique() == 1
    assert summary.iloc[0]["block_bootstrap_status"].startswith("available")

