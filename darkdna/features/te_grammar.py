"""Transposable-element grammar features."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from darkdna.io.bed import collect_overlapping_values, interval_overlap_bp, overlap_fraction
from darkdna.utils.stats import shannon_entropy


def compute_te_grammar_for_region(chrom: str, start: int, end: int, te_annotations: pd.DataFrame | None) -> dict[str, float | str | bool]:
    if te_annotations is None or te_annotations.empty:
        return {
            "TE_overlap_fraction": 0.0,
            "TE_family": "",
            "TE_class": "",
            "TE_superfamily": "",
            "TE_nesting_score": math.nan,
            "TE_fragmentation_score": math.nan,
            "TE_orientation_architecture": "",
            "TE_family_mosaic_score": math.nan,
            "TE_boundary_score": math.nan,
            "distance_to_nearest_TE_boundary": math.nan,
            "LTR_internal_discordance": math.nan,
            "TE_derived_candidate_flag": False,
            "TE_age_divergence_proxy": math.nan,
            "local_TE_density": 0.0,
            "local_TE_diversity": 0.0,
            "TE_strand_orientation_entropy": math.nan,
        }
    subset = te_annotations[te_annotations["chrom"].astype(str) == str(chrom)].copy()
    overlaps = []
    for row in subset.itertuples():
        bp = interval_overlap_bp(start, end, int(row.start), int(row.end))
        if bp > 0:
            overlaps.append((row, bp))
    length = max(1, end - start)
    overlap_bp = sum(bp for _, bp in overlaps)
    families = [str(getattr(row, "family", "")) for row, _ in overlaps if pd.notna(getattr(row, "family", None))]
    classes = [str(getattr(row, "class", "")) for row, _ in overlaps if pd.notna(getattr(row, "class", None))]
    supers = [str(getattr(row, "superfamily", "")) for row, _ in overlaps if pd.notna(getattr(row, "superfamily", None))]
    strands = [str(getattr(row, "strand", ".")) for row, _ in overlaps if pd.notna(getattr(row, "strand", None))]
    boundaries = []
    for row in subset.itertuples():
        boundaries.extend([abs(start - int(row.end)), abs(end - int(row.start)), abs(start - int(row.start)), abs(end - int(row.end))])
    divergence_vals = []
    for row, _ in overlaps:
        try:
            divergence_vals.append(float(getattr(row, "divergence")))
        except Exception:
            pass
    return {
        "TE_overlap_fraction": overlap_bp / length,
        "TE_family": ";".join(dict.fromkeys(families)),
        "TE_class": ";".join(dict.fromkeys(classes)),
        "TE_superfamily": ";".join(dict.fromkeys(supers)),
        "TE_nesting_score": max(0.0, len(overlaps) - len(set(families))) / max(1, len(overlaps)) if overlaps else 0.0,
        "TE_fragmentation_score": len(overlaps) / max(1, overlap_bp / 100.0) if overlaps else 0.0,
        "TE_orientation_architecture": ";".join(strands),
        "TE_family_mosaic_score": len(set(families)) / max(1, len(families)) if families else 0.0,
        "TE_boundary_score": 1.0 / (1.0 + min(boundaries)) if boundaries else 0.0,
        "distance_to_nearest_TE_boundary": float(min(boundaries)) if boundaries else math.nan,
        "LTR_internal_discordance": math.nan,
        "TE_derived_candidate_flag": overlap_bp / length > 0.1,
        "TE_age_divergence_proxy": float(np.nanmean(divergence_vals)) if divergence_vals else math.nan,
        "local_TE_density": overlap_fraction(chrom, start - length, end + length, te_annotations),
        "local_TE_diversity": len(set(families)),
        "TE_strand_orientation_entropy": shannon_entropy(strands) if strands else 0.0,
    }


def annotate_te_grammar(windows: pd.DataFrame, te_annotations: pd.DataFrame | None) -> pd.DataFrame:
    rows = []
    for row in windows.itertuples():
        rows.append({"region_id": row.region_id, **compute_te_grammar_for_region(row.chrom, int(row.start), int(row.end), te_annotations)})
    return pd.DataFrame(rows)
