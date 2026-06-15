"""BED, bedGraph, and interval helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


BED_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]


def read_bed(path: str | Path | None, names: list[str] | None = None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame(columns=names or BED_COLUMNS)
    rows: list[list[str]] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            rows.append(line.rstrip("\n").split("\t"))
    if not rows:
        return pd.DataFrame(columns=names or BED_COLUMNS)
    ncols = max(len(row) for row in rows)
    base = names or BED_COLUMNS + [f"field_{idx}" for idx in range(len(BED_COLUMNS), ncols)]
    cols = base[:ncols]
    padded = [row + [None] * (ncols - len(row)) for row in rows]
    df = pd.DataFrame(padded, columns=cols)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    if "strand" not in df.columns:
        df["strand"] = "."
    return df


def read_bedgraph(path: str | Path | None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame(columns=["chrom", "start", "end", "value"])
    df = read_bed(path, names=["chrom", "start", "end", "value"])
    if not df.empty:
        df["value"] = pd.to_numeric(df["value"], errors="coerce")
    return df


def write_bed(df: pd.DataFrame, path: str | Path, columns: list[str] | None = None) -> None:
    cols = columns or [c for c in ["chrom", "start", "end", "region_id", "score", "strand"] if c in df.columns]
    with Path(path).open("w", encoding="utf-8") as handle:
        for row in df[cols].itertuples(index=False, name=None):
            handle.write("\t".join("" if pd.isna(v) else str(v) for v in row) + "\n")


def write_bedgraph(df: pd.DataFrame, path: str | Path, value_col: str) -> None:
    with Path(path).open("w", encoding="utf-8") as handle:
        for row in df.itertuples():
            value = getattr(row, value_col, np.nan)
            if pd.isna(value):
                value = 0
            handle.write(f"{row.chrom}\t{int(row.start)}\t{int(row.end)}\t{float(value):.6g}\n")


def interval_overlap_bp(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    return max(0, min(int(end_a), int(end_b)) - max(int(start_a), int(start_b)))


def overlaps_any(chrom: str, start: int, end: int, intervals: pd.DataFrame) -> bool:
    if intervals is None or intervals.empty:
        return False
    subset = intervals[intervals["chrom"].astype(str) == str(chrom)]
    for row in subset.itertuples():
        if interval_overlap_bp(start, end, int(row.start), int(row.end)) > 0:
            return True
    return False


def overlap_fraction(chrom: str, start: int, end: int, intervals: pd.DataFrame) -> float:
    if intervals is None or intervals.empty:
        return 0.0
    length = max(1, int(end) - int(start))
    subset = intervals[intervals["chrom"].astype(str) == str(chrom)]
    overlap = 0
    for row in subset.itertuples():
        overlap += interval_overlap_bp(start, end, int(row.start), int(row.end))
    return min(1.0, overlap / length)


def weighted_interval_mean(chrom: str, start: int, end: int, intervals: pd.DataFrame, value_col: str = "value") -> float:
    if intervals is None or intervals.empty or value_col not in intervals.columns:
        return np.nan
    length = max(1, int(end) - int(start))
    total = 0.0
    weight = 0
    subset = intervals[intervals["chrom"].astype(str) == str(chrom)]
    for row in subset.itertuples():
        bp = interval_overlap_bp(start, end, int(row.start), int(row.end))
        if bp > 0:
            total += bp * float(getattr(row, value_col))
            weight += bp
    if weight == 0:
        return np.nan
    return float(total / length)


def nearest_interval_distance(chrom: str, start: int, end: int, intervals: pd.DataFrame, point_cols: tuple[str, str] = ("start", "end")) -> float:
    if intervals is None or intervals.empty:
        return np.nan
    center = (int(start) + int(end)) // 2
    subset = intervals[intervals["chrom"].astype(str) == str(chrom)]
    distances = []
    for row in subset.itertuples():
        a = int(getattr(row, point_cols[0]))
        b = int(getattr(row, point_cols[1]))
        if interval_overlap_bp(start, end, a, b) > 0:
            distances.append(0)
        else:
            distances.append(min(abs(center - a), abs(center - b)))
    return float(min(distances)) if distances else np.nan


def collect_overlapping_values(
    chrom: str,
    start: int,
    end: int,
    intervals: pd.DataFrame,
    value_col: str,
    max_values: int = 5,
) -> list[str]:
    if intervals is None or intervals.empty or value_col not in intervals.columns:
        return []
    values: list[str] = []
    subset = intervals[intervals["chrom"].astype(str) == str(chrom)]
    for row in subset.itertuples():
        if interval_overlap_bp(start, end, int(row.start), int(row.end)) > 0:
            value = getattr(row, value_col)
            if pd.notna(value) and str(value) not in values:
                values.append(str(value))
        if len(values) >= max_values:
            break
    return values


def concatenate_flags(flags: Iterable[str]) -> str:
    unique = [flag for flag in dict.fromkeys(flags) if flag]
    return ";".join(unique)
