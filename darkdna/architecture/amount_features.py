"""Amount and interval-burden features for sequence-indifferent architecture."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


def _overlap_bp(start: int, end: int, other_start: int, other_end: int) -> int:
    return max(0, min(end, other_end) - max(start, other_start))


def _overlaps(row: pd.Series, intervals: pd.DataFrame | None) -> pd.DataFrame:
    if intervals is None or intervals.empty or not {"chrom", "start", "end"}.issubset(intervals.columns):
        return pd.DataFrame()
    same = intervals.loc[intervals["chrom"].astype(str) == str(row["chrom"])].copy()
    if same.empty:
        return same
    mask = (pd.to_numeric(same["end"], errors="coerce") > int(row["start"])) & (
        pd.to_numeric(same["start"], errors="coerce") < int(row["end"])
    )
    return same.loc[mask]


def compute_amount_features(
    intervals: pd.DataFrame,
    *,
    repeats: pd.DataFrame | None = None,
    heterochromatin: pd.DataFrame | None = None,
    replication_domains: pd.DataFrame | None = None,
    genome_sizes: dict[str, int] | None = None,
) -> pd.DataFrame:
    if not {"chrom", "start", "end"}.issubset(intervals.columns):
        raise ValueError("Mode B intervals require chrom, start, and end")
    frame = intervals.copy().reset_index(drop=True)
    if "region_id" not in frame.columns:
        if "locus_id" in frame.columns:
            frame["region_id"] = frame["locus_id"].astype(str)
        else:
            frame["region_id"] = frame.apply(lambda row: f"{row['chrom']}:{int(row['start'])}-{int(row['end'])}", axis=1)
    inferred_sizes = frame.groupby(frame["chrom"].astype(str))["end"].max().astype(int).to_dict()
    sizes = {**inferred_sizes, **(genome_sizes or {})}
    genome_bp = sum(max(0, int(value)) for value in sizes.values()) or 1
    family_counts = frame.get("family", pd.Series("", index=frame.index)).fillna("").astype(str).value_counts()
    rows: list[dict[str, object]] = []
    for _, row in frame.iterrows():
        start, end = int(row["start"]), int(row["end"])
        length = max(0, end - start)
        repeat_rows = _overlaps(row, repeats)
        repeat_bp = sum(
            _overlap_bp(start, end, int(item.start), int(item.end)) for item in repeat_rows.itertuples()
        ) if not repeat_rows.empty else 0
        family_column = next((column for column in ("family", "TE_family", "repeat_family") if column in repeat_rows.columns), None)
        overlapping_families = sorted(set(repeat_rows[family_column].dropna().astype(str))) if family_column else []
        if family_column and repeats is not None:
            burdens = (
                repeats.assign(_length=pd.to_numeric(repeats["end"], errors="coerce") - pd.to_numeric(repeats["start"], errors="coerce"))
                .groupby(family_column)["_length"]
                .sum()
            )
            repeat_family_total_bp = float(max((burdens.get(family, 0.0) for family in overlapping_families), default=math.nan))
        else:
            repeat_family_total_bp = math.nan
        hetero_rows = _overlaps(row, heterochromatin)
        hetero_bp = sum(_overlap_bp(start, end, int(item.start), int(item.end)) for item in hetero_rows.itertuples()) if not hetero_rows.empty else 0
        replication_rows = _overlaps(row, replication_domains)
        replication_length = (
            float(max(int(item.end) - int(item.start) for item in replication_rows.itertuples()))
            if not replication_rows.empty
            else math.nan
        )
        same_chrom = frame.loc[frame["chrom"].astype(str) == str(row["chrom"])].copy()
        centres = (pd.to_numeric(same_chrom["start"]) + pd.to_numeric(same_chrom["end"])) / 2.0
        centre = (start + end) / 2.0
        distances = (centres - centre).abs()
        nonself = distances.loc[distances > 0]
        nearest_spacing = float(nonself.min()) if not nonself.empty else math.nan
        adjacent = int((nonself <= max(length * 2, 1)).sum())
        family = str(row.get("family", "") or "")
        rows.append(
            {
                "region_id": str(row["region_id"]),
                "source_locus_id": str(row.get("locus_id", "") or ""),
                "source_representative_region_id": str(row.get("representative_region_id", "") or ""),
                "chrom": str(row["chrom"]),
                "start": start,
                "end": end,
                "interval_length": length,
                "locus_length": length,
                "repeat_array_length": int(repeat_bp),
                "repeat_family_total_bp": repeat_family_total_bp,
                "local_repeat_fraction": repeat_bp / length if length else math.nan,
                "chromosome_fraction_occupied": length / max(1, int(sizes.get(str(row["chrom"]), end))),
                "genome_fraction_occupied": length / genome_bp,
                "heterochromatin_overlap": hetero_bp / length if length and heterochromatin is not None else math.nan,
                "heterochromatin_status": "available" if heterochromatin is not None and not heterochromatin.empty else "unavailable_missing_track",
                "replication_domain_length": replication_length,
                "replication_domain_status": "available" if replication_domains is not None and not replication_domains.empty else "unavailable_missing_track",
                "candidate_family_burden": int(family_counts.get(family, 1)) if family else math.nan,
                "adjacent_copy_count": adjacent,
                "inter_copy_spacing": nearest_spacing,
                "overlapping_repeat_families": ";".join(overlapping_families),
                "amount_feature_status": "available_observed_interval_features",
                "amount_feature_caveat": "Amount features are descriptive until controlled quantity perturbations or comparative constraints are supplied.",
            }
        )
    return pd.DataFrame(rows)
