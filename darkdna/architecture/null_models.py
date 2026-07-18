"""Block-aware nulls for Mode B quantity scores."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from darkdna.nulls.calibration import infer_genomic_blocks
from darkdna.utils.stats import empirical_p_value


def build_architecture_nulls(
    scores: pd.DataFrame,
    features: pd.DataFrame,
    controls: pd.DataFrame | None = None,
    *,
    block_size_bp: int = 100_000,
    minimum_independent_blocks: int = 5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = features.merge(scores[["region_id", "quantity_anomaly_score"]], on="region_id", how="inner")
    frame["calibration_block_id"] = infer_genomic_blocks(frame, block_size_bp)
    rows = []
    residuals = []
    for _, row in frame.iterrows():
        matched_controls = frame.loc[frame["calibration_block_id"].astype(str) != str(row["calibration_block_id"])]
        values = pd.to_numeric(matched_controls["quantity_anomaly_score"], errors="coerce").dropna().to_numpy(dtype=float)
        observed = float(row["quantity_anomaly_score"])
        mean = float(np.mean(values)) if values.size else math.nan
        std = float(np.std(values, ddof=1)) if values.size > 1 else math.nan
        zscore = (observed - mean) / std if np.isfinite(std) and std > 0 else math.nan
        blocks = int(matched_controls["calibration_block_id"].nunique())
        status = "available_block_calibrated" if blocks >= minimum_independent_blocks and np.isfinite(zscore) else "partial_not_for_promotion"
        p_value = empirical_p_value(observed, values, higher=True) if values.size else math.nan
        rows.append(
            {
                "region_id": str(row["region_id"]),
                "null_model_id": "matched_interval_independent_blocks",
                "observed_quantity_score": observed,
                "null_mean": mean,
                "null_std": std,
                "architecture_null_zscore": zscore,
                "architecture_null_empirical_p": p_value,
                "null_sample_size": int(values.size),
                "independent_block_count": blocks,
                "architecture_null_status": status,
            }
        )
        residuals.append(
            {
                "region_id": str(row["region_id"]),
                "observed_quantity_score": observed,
                "predicted_matched_quantity_score": mean,
                "architecture_residual": observed - mean if np.isfinite(mean) else math.nan,
                "architecture_residual_zscore": zscore,
                "residual_status": status,
            }
        )
    if controls is not None and not controls.empty:
        for region_id, group in controls.groupby("region_id", sort=False):
            native = group.loc[group["transformation"] == "native", "prediction"]
            observed = float(native.iloc[0]) if not native.empty else math.nan
            for family, model_id in [
                ("identity", "equal_length_replacement_panel"),
                ("length", "length_titration_panel"),
                ("copy_number", "copy_number_titration_panel"),
            ]:
                subset = group.loc[group["control_family"] == family]
                values = pd.to_numeric(subset["prediction"], errors="coerce").dropna().to_numpy(dtype=float)
                rows.append(
                    {
                        "region_id": str(region_id),
                        "null_model_id": model_id,
                        "observed_quantity_score": observed,
                        "null_mean": float(np.mean(values)) if values.size else math.nan,
                        "null_std": float(np.std(values, ddof=1)) if values.size > 1 else math.nan,
                        "architecture_null_zscore": math.nan,
                        "architecture_null_empirical_p": math.nan,
                        "null_sample_size": int(values.size),
                        "independent_block_count": 0,
                        "architecture_null_status": "available_model_perturbation_not_exchangeable_null" if values.size else "unavailable",
                        "architecture_null_reason": "A deterministic perturbation panel is not treated as an exchangeable null distribution; no p-value is emitted.",
                    }
                )
    return pd.DataFrame(rows), pd.DataFrame(residuals)
