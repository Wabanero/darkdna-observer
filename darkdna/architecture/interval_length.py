"""Interval-length summaries."""

from __future__ import annotations

import pandas as pd


def compute_interval_lengths(intervals: pd.DataFrame) -> pd.DataFrame:
    output = intervals[["region_id", "chrom", "start", "end"]].copy()
    output["interval_length"] = pd.to_numeric(output["end"], errors="coerce") - pd.to_numeric(output["start"], errors="coerce")
    output["interval_length_status"] = "available_observed_coordinates"
    return output
