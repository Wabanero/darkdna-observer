"""Presence/absence variation features."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


def _to_long(table: pd.DataFrame) -> pd.DataFrame:
    if {"region_id", "sample_id", "present"}.issubset(table.columns):
        return table[["region_id", "sample_id", "present"]].copy()
    if "region_id" in table.columns:
        value_columns = [column for column in table.columns if column != "region_id"]
        return table.melt(id_vars="region_id", value_vars=value_columns, var_name="sample_id", value_name="present")
    return pd.DataFrame(columns=["region_id", "sample_id", "present"])


def compute_presence_absence_features(intervals: pd.DataFrame, table: pd.DataFrame | None = None) -> pd.DataFrame:
    long = _to_long(table) if table is not None and not table.empty else pd.DataFrame()
    rows = []
    for region_id in intervals["region_id"].astype(str):
        subset = long.loc[long["region_id"].astype(str) == region_id] if not long.empty else pd.DataFrame()
        values = pd.to_numeric(subset.get("present", pd.Series(dtype=float)), errors="coerce").dropna().to_numpy(dtype=float)
        if values.size:
            present = values > 0
            frequency = float(np.mean(present))
            status, reason = "available", "Presence/absence frequency from supplied samples."
        else:
            frequency = math.nan
            status, reason = "unavailable_missing_presence_absence_input", "No presence/absence observations were supplied for this interval."
        rows.append(
            {
                "region_id": region_id,
                "presence_absence_frequency": frequency,
                "deletion_frequency_from_presence_absence": 1.0 - frequency if np.isfinite(frequency) else math.nan,
                "presence_absence_sample_size": int(values.size),
                "presence_absence_status": status,
                "presence_absence_reason": reason,
            }
        )
    return pd.DataFrame(rows)
