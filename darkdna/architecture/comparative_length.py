"""Length conservation despite sequence turnover."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


def compute_length_conservation(intervals: pd.DataFrame, syntenic: pd.DataFrame | None = None) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for region_id in intervals["region_id"].astype(str):
        subset = syntenic.loc[syntenic["region_id"].astype(str) == region_id].copy() if syntenic is not None and not syntenic.empty and "region_id" in syntenic.columns else pd.DataFrame()
        if subset.empty:
            rows.append(
                {
                    "region_id": region_id,
                    "mean_interval_length": math.nan,
                    "interval_length_cv": math.nan,
                    "normalized_length_variance": math.nan,
                    "sequence_identity": math.nan,
                    "alignment_coverage": math.nan,
                    "length_conservation_score": math.nan,
                    "sequence_conservation_score": math.nan,
                    "length_minus_sequence_conservation": math.nan,
                    "synteny_confidence": math.nan,
                    "assembly_confidence": math.nan,
                    "length_conservation_status": "unavailable_missing_syntenic_intervals",
                    "length_conservation_reason": "Poor or absent alignment is not evidence of functional sequence turnover.",
                }
            )
            continue
        if "interval_length" in subset.columns:
            lengths = pd.to_numeric(subset["interval_length"], errors="coerce").dropna().to_numpy(dtype=float)
        else:
            lengths = (pd.to_numeric(subset["end"], errors="coerce") - pd.to_numeric(subset["start"], errors="coerce")).dropna().to_numpy(dtype=float)
        mean_length = float(np.mean(lengths)) if lengths.size else math.nan
        cv = float(np.std(lengths, ddof=1) / mean_length) if lengths.size > 1 and mean_length > 0 else math.nan
        normalized_variance = float(np.var(lengths, ddof=1) / (mean_length**2)) if lengths.size > 1 and mean_length > 0 else math.nan
        identity = float(pd.to_numeric(subset.get("sequence_identity", pd.Series(dtype=float)), errors="coerce").mean())
        coverage = float(pd.to_numeric(subset.get("alignment_coverage", pd.Series(dtype=float)), errors="coerce").mean())
        synteny_confidence = float(pd.to_numeric(subset.get("synteny_confidence", pd.Series(dtype=float)), errors="coerce").mean())
        assembly_confidence = float(pd.to_numeric(subset.get("assembly_confidence", pd.Series(dtype=float)), errors="coerce").mean())
        length_score = 1.0 / (1.0 + cv) if np.isfinite(cv) else math.nan
        sequence_score = identity if np.isfinite(identity) else math.nan
        confident = np.isfinite(synteny_confidence) and synteny_confidence >= 0.5 and np.isfinite(coverage) and coverage >= 0.5
        rows.append(
            {
                "region_id": region_id,
                "mean_interval_length": mean_length,
                "interval_length_cv": cv,
                "normalized_length_variance": normalized_variance,
                "sequence_identity": identity,
                "alignment_coverage": coverage,
                "length_conservation_score": length_score if confident else math.nan,
                "sequence_conservation_score": sequence_score if confident else math.nan,
                "length_minus_sequence_conservation": length_score - sequence_score if confident and np.isfinite(sequence_score) else math.nan,
                "synteny_confidence": synteny_confidence,
                "assembly_confidence": assembly_confidence,
                "length_conservation_status": "available_confident_synteny" if confident else "unavailable_low_synteny_or_alignment_confidence",
                "length_conservation_reason": "Length conservation is interpreted only with adequate synteny and alignment coverage.",
            }
        )
    return pd.DataFrame(rows)
