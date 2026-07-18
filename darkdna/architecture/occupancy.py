"""Physical-occupancy interval summaries."""

from __future__ import annotations

import math

import pandas as pd

from darkdna.architecture.amount_features import _overlap_bp


def compute_occupancy_features(intervals: pd.DataFrame, occupancy: pd.DataFrame | None = None) -> pd.DataFrame:
    rows = []
    for _, interval in intervals.iterrows():
        if occupancy is None or occupancy.empty or not {"chrom", "start", "end"}.issubset(occupancy.columns):
            rows.append({"region_id": str(interval["region_id"]), "physical_occupancy": math.nan, "occupancy_status": "unavailable_missing_track"})
            continue
        same = occupancy.loc[occupancy["chrom"].astype(str) == str(interval["chrom"])].copy()
        same = same.loc[(same["end"] > int(interval["start"])) & (same["start"] < int(interval["end"]))]
        value_column = "value" if "value" in same.columns else None
        weights = [_overlap_bp(int(interval["start"]), int(interval["end"]), int(row.start), int(row.end)) for row in same.itertuples()]
        if value_column and sum(weights) > 0:
            value = sum(weight * float(row[value_column]) for weight, (_, row) in zip(weights, same.iterrows())) / sum(weights)
        else:
            value = sum(weights) / max(1, int(interval["end"]) - int(interval["start"])) if weights else math.nan
        rows.append({"region_id": str(interval["region_id"]), "physical_occupancy": value, "occupancy_status": "available" if weights else "unavailable_no_overlap"})
    return pd.DataFrame(rows)
