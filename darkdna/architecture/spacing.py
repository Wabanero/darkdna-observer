"""Anchor-distance and spacing-constraint features."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


def compute_spacing_features(intervals: pd.DataFrame, anchors: pd.DataFrame | None = None) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    raw_distances: list[float] = []
    provisional: list[dict[str, object]] = []
    for _, interval in intervals.iterrows():
        same = anchors.loc[anchors["chrom"].astype(str) == str(interval["chrom"])].copy() if anchors is not None and not anchors.empty and {"chrom", "start", "end"}.issubset(anchors.columns) else pd.DataFrame()
        left = same.loc[pd.to_numeric(same["end"], errors="coerce") <= int(interval["start"])].sort_values("end").tail(1) if not same.empty else pd.DataFrame()
        right = same.loc[pd.to_numeric(same["start"], errors="coerce") >= int(interval["end"])].sort_values("start").head(1) if not same.empty else pd.DataFrame()
        if left.empty or right.empty:
            left_distance = right_distance = anchor_distance = math.nan
            status = "unavailable_missing_flanking_anchors"
        else:
            left_distance = int(interval["start"]) - int(left.iloc[0]["end"])
            right_distance = int(right.iloc[0]["start"]) - int(interval["end"])
            anchor_distance = int(right.iloc[0]["start"]) - int(left.iloc[0]["end"])
            raw_distances.append(float(anchor_distance))
            status = "available_observed_spacing"
        provisional.append(
            {
                "region_id": str(interval["region_id"]),
                "left_anchor_distance": left_distance,
                "right_anchor_distance": right_distance,
                "anchor_to_anchor_distance": anchor_distance,
                "spacing_conservation": math.nan,
                "orientation_conservation": math.nan,
                "sequence_turnover": math.nan,
                "spacing_status": status,
                "spacing_reason": "Comparative spacing conservation requires homologous anchor sets across assemblies or species.",
            }
        )
    mean = float(np.mean(raw_distances)) if raw_distances else math.nan
    std = float(np.std(raw_distances, ddof=1)) if len(raw_distances) > 1 else math.nan
    for record in provisional:
        value = float(record["anchor_to_anchor_distance"])
        record["spacing_null_zscore"] = (value - mean) / std if np.isfinite(value) and np.isfinite(std) and std > 0 else math.nan
        rows.append(record)
    return pd.DataFrame(rows)
