"""Bulk heterochromatin amount summaries."""

from __future__ import annotations

import math

import pandas as pd

from darkdna.architecture.amount_features import _overlap_bp


def compute_heterochromatin_mass(intervals: pd.DataFrame, heterochromatin: pd.DataFrame | None = None) -> pd.DataFrame:
    rows = []
    for _, interval in intervals.iterrows():
        if heterochromatin is None or heterochromatin.empty:
            rows.append({"region_id": str(interval["region_id"]), "heterochromatic_bp": math.nan, "heterochromatin_mass_status": "unavailable_missing_track"})
            continue
        same = heterochromatin.loc[heterochromatin["chrom"].astype(str) == str(interval["chrom"])]
        bp = sum(_overlap_bp(int(interval["start"]), int(interval["end"]), int(row.start), int(row.end)) for row in same.itertuples())
        rows.append({"region_id": str(interval["region_id"]), "heterochromatic_bp": int(bp), "heterochromatin_mass_status": "available"})
    return pd.DataFrame(rows)
