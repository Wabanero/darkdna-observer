"""Masking helpers for windows."""

from __future__ import annotations

import pandas as pd

from darkdna.io.bed import overlaps_any


def filter_overlapping(windows: pd.DataFrame, intervals: pd.DataFrame, keep_flag: bool = False) -> pd.DataFrame:
    out = windows.copy()
    mask = [overlaps_any(row.chrom, row.start, row.end, intervals) for row in out.itertuples()]
    out["masked_overlap"] = mask
    return out if keep_flag else out.loc[~pd.Series(mask, index=out.index)].copy()


def mark_overlap(windows: pd.DataFrame, intervals: pd.DataFrame, column: str) -> pd.DataFrame:
    out = windows.copy()
    out[column] = [overlaps_any(row.chrom, row.start, row.end, intervals) for row in out.itertuples()]
    return out
