"""Classical covariates and artifact-risk helpers."""

from __future__ import annotations

import pandas as pd


BLOCK_CV_SIZE_BP = 100_000

CLASSICAL_COVARIATES = [
    "gc_content",
    "n_fraction",
    "CpG_density",
    "simple_repeat_fraction",
    "low_complexity_mask_fraction",
    "TE_overlap_fraction",
    "overlaps_TE",
    "overlaps_cCRE",
    "overlaps_enhancer",
    "overlaps_promoter",
    "overlaps_exon",
    "overlaps_intron",
    "overlaps_utr",
    "distance_to_nearest_tss",
    "overlaps_blacklist",
    "mappability",
    "overlaps_assembly_gap",
    "overlaps_segmental_duplication",
    "scaffold_edge_distance",
    "window_size",
    "local_TE_density",
    "local_GC_background",
    "local_CpG_background",
]


def _block_id(row: pd.Series, block_size: int = BLOCK_CV_SIZE_BP) -> str:
    chrom = str(row.get("chrom", "unknown"))
    try:
        start = int(row.get("start", 0) or 0)
    except Exception:
        start = 0
    return f"{chrom}:block_{start // block_size}"


def build_classical_covariates(windows: pd.DataFrame, features: pd.DataFrame | None = None) -> pd.DataFrame:
    base = windows.copy()
    if features is not None and not features.empty:
        base = base.merge(features, on="region_id", how="left", suffixes=("", "_feature"))
    if "TE_overlap_fraction" not in base.columns and "overlaps_TE" in base.columns:
        base["TE_overlap_fraction"] = base["overlaps_TE"].astype(float)
    if "local_TE_density" not in base.columns:
        base["local_TE_density"] = base.get("TE_overlap_fraction", 0.0)
    if "local_GC_background" not in base.columns:
        base["local_GC_background"] = base.get("gc_content", base.get("gc_content_feature", 0.0))
    if "local_CpG_background" not in base.columns:
        base["local_CpG_background"] = base.get("CpG_density", 0.0)
    cols = ["region_id"]
    for metadata_col in ["chrom", "start", "end"]:
        if metadata_col in base.columns:
            cols.append(metadata_col)
    if "chrom" in base.columns and "start" in base.columns:
        base["block_id"] = base.apply(_block_id, axis=1)
        base["chromosome_cv_fold"] = base["chrom"].astype(str)
        cols.extend(["block_id", "chromosome_cv_fold"])
    for col in CLASSICAL_COVARIATES:
        if col not in base.columns:
            base[col] = 0.0
        cols.append(col)
    return base[cols]


def artifact_flags_to_list(flags: str | float | None) -> list[str]:
    if flags is None or pd.isna(flags) or flags == "":
        return []
    return [flag for flag in str(flags).split(";") if flag]


def artifact_risk_score(flags: str | float | None) -> float:
    items = artifact_flags_to_list(flags)
    return min(1.0, len(items) / 5.0)
