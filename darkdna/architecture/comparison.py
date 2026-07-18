"""Mode A versus Mode B comparison without collapsing the axes."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


def _mode_a_summary(mode_a: pd.DataFrame) -> pd.DataFrame:
    if mode_a is None or mode_a.empty:
        return pd.DataFrame(columns=["region_id", "sequence_specific_score"])
    frame = mode_a.copy()
    if "region_id" not in frame.columns:
        if "locus_id" in frame.columns:
            frame["region_id"] = frame["locus_id"].astype(str)
        else:
            return pd.DataFrame(columns=["region_id", "sequence_specific_score"])
    if "residual_zscore" in frame.columns:
        values = pd.to_numeric(frame["residual_zscore"], errors="coerce").abs()
        frame = frame.assign(_value=values).groupby("region_id", as_index=False)["_value"].max()
        frame = frame.rename(columns={"_value": "sequence_specific_score"})
    elif "support_score" in frame.columns:
        frame = frame[["region_id", "support_score"]].rename(columns={"support_score": "sequence_specific_score"})
    elif "sequence_identity_sensitivity" in frame.columns:
        frame = frame[["region_id", "sequence_identity_sensitivity"]].rename(columns={"sequence_identity_sensitivity": "sequence_specific_score"})
    else:
        frame = frame[["region_id"]].assign(sequence_specific_score=math.nan)
    return frame


def compare_sequence_vs_quantity(mode_a: pd.DataFrame, mode_b: pd.DataFrame) -> pd.DataFrame:
    sequence = _mode_a_summary(mode_a)
    quantity_columns = [
        column
        for column in [
            "region_id",
            "quantity_anomaly_score",
            "sequence_identity_sensitivity",
            "length_sensitivity",
            "copy_number_sensitivity",
            "spacing_null_zscore",
            "sequence_indifference_score",
            "sequence_quantity_interaction_score",
            "sequence_indifferent_candidate",
            "dominant_mode",
            "mode_b_score_status",
        ]
        if column in mode_b.columns
    ]
    merged = mode_b[quantity_columns].merge(sequence, on="region_id", how="left")
    merged["mode_a_available"] = pd.to_numeric(merged["sequence_specific_score"], errors="coerce").notna()
    merged["mode_b_available"] = pd.to_numeric(merged.get("quantity_anomaly_score"), errors="coerce").notna()
    merged["comparison_status"] = np.where(
        merged["mode_a_available"] & merged["mode_b_available"],
        "both_axes_available_kept_separate",
        "partial_axis_availability",
    )
    merged["comparison_caveat"] = "Mode A and Mode B are separate evidence axes; no universal DarkDNA score is calculated."
    return merged
