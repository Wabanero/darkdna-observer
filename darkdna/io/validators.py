"""Input validation helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


REGION_COLUMNS = [
    "region_id",
    "chrom",
    "start",
    "end",
    "window_size",
    "scale_level",
    "parent_region_id",
    "child_region_ids",
    "is_dark",
    "overlaps_exon",
    "overlaps_promoter",
    "overlaps_intron",
    "overlaps_utr",
    "overlaps_TE",
    "TE_family",
    "overlaps_cCRE",
    "overlaps_enhancer",
    "overlaps_blacklist",
    "overlaps_assembly_gap",
    "overlaps_segmental_duplication",
    "nearest_gene",
    "distance_to_nearest_tss",
    "gc_content",
    "n_fraction",
    "mappability",
    "low_complexity_mask_fraction",
    "scaffold_edge_distance",
    "artifact_risk_flags",
]


def require_path(path: str | Path | None, label: str) -> Path:
    if not path:
        raise ValueError(f"{label} is required")
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{label} not found: {p}")
    return p


def ensure_region_schema(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in REGION_COLUMNS:
        if col not in out.columns:
            out[col] = None
    return out[REGION_COLUMNS + [c for c in out.columns if c not in REGION_COLUMNS]]
