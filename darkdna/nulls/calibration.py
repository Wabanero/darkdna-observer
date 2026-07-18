"""Block-aware severe-null calibration for candidate score tables."""

from __future__ import annotations

import math
from collections.abc import Iterable

import numpy as np
import pandas as pd

from darkdna.nulls.registry import assess_null_availability
from darkdna.utils.stats import empirical_p_value


MATCH_COLUMNS: dict[str, tuple[str, ...]] = {
    "same_length_gc_matched": ("length", "window_size", "gc_content"),
    "local_genomic_matched": ("start",),
    "gene_tss_distance_matched": ("distance_to_nearest_tss",),
    "te_family_matched": ("TE_family",),
    "te_subfamily_matched": ("TE_subfamily",),
    "te_age_matched": ("TE_age", "TE_divergence"),
    "repeat_density_matched": ("local_TE_density", "simple_repeat_fraction", "TE_overlap_fraction"),
    "chromatin_compartment_matched": ("chromatin_compartment",),
    "replication_timing_matched": ("replication_timing",),
    "recombination_matched": ("recombination_rate",),
    "mutation_rate_matched": ("mutation_rate",),
    "damage_environment_matched": ("damage_environment",),
    "mappability_matched": ("mappability",),
    "assembly_confidence_matched": ("assembly_confidence",),
    "syntenic_ortholog": ("synteny_group",),
    "population_frequency_matched": ("population_frequency",),
    "copy_number_matched": ("copy_number",),
    "presence_absence_matched": ("presence_absence_frequency",),
}


EXACT_MATCH_COLUMNS = {"TE_family", "TE_subfamily", "chromatin_compartment", "synteny_group", "damage_environment"}


def infer_genomic_blocks(features: pd.DataFrame, block_size_bp: int = 100_000) -> pd.Series:
    if {"chrom", "start"}.issubset(features.columns):
        starts = pd.to_numeric(features["start"], errors="coerce").fillna(0).astype(int)
        return features["chrom"].astype(str) + ":" + (starts // max(1, block_size_bp)).astype(str)
    return pd.Series([f"independent_row:{idx}" for idx in range(len(features))], index=features.index)


def _score_columns(scores: pd.DataFrame) -> list[str]:
    return [
        column
        for column in scores.columns
        if column != "region_id" and column.endswith("_score") and pd.api.types.is_numeric_dtype(scores[column])
    ]


def _usable_match_columns(frame: pd.DataFrame, model_id: str) -> list[str]:
    candidates = MATCH_COLUMNS.get(model_id, ())
    if model_id == "same_length_gc_matched":
        length = "length" if "length" in frame.columns else "window_size" if "window_size" in frame.columns else None
        return [column for column in [length, "gc_content"] if column and column in frame.columns]
    return [column for column in candidates if column in frame.columns]


def _matched_control_indices(
    frame: pd.DataFrame,
    row_index: int,
    model_id: str,
    n_controls: int,
) -> list[int]:
    row = frame.iloc[row_index]
    controls = frame.loc[frame.index != frame.index[row_index]].copy()
    if "calibration_block_id" in frame.columns:
        controls = controls.loc[controls["calibration_block_id"].astype(str) != str(row["calibration_block_id"])]
    if controls.empty:
        return []
    if model_id == "local_genomic_matched" and "chrom" in frame.columns:
        same_chrom = controls.loc[controls["chrom"].astype(str) == str(row.get("chrom", ""))]
        if not same_chrom.empty:
            controls = same_chrom
    distance = pd.Series(0.0, index=controls.index)
    used = 0
    for column in _usable_match_columns(frame, model_id):
        if column in EXACT_MATCH_COLUMNS:
            exact = controls[column].fillna("").astype(str) == str(row.get(column, ""))
            if exact.any():
                controls = controls.loc[exact]
                distance = distance.loc[controls.index]
                used += 1
            continue
        values = pd.to_numeric(controls[column], errors="coerce")
        observed = pd.to_numeric(pd.Series([row.get(column)]), errors="coerce").iloc[0]
        scale = float(values.std(ddof=1))
        if np.isfinite(observed) and np.isfinite(scale) and scale > 0:
            distance = distance.loc[controls.index] + ((values - observed).abs() / scale).fillna(0.0)
            used += 1
    if used == 0 or controls.empty:
        return []
    return distance.sort_values().head(max(1, n_controls)).index.tolist()


def _precompute_control_indices(frame: pd.DataFrame, model_id: str, n_controls: int) -> list[list[int]]:
    """Vectorized pairwise matching, excluding the focal genomic block."""

    columns = _usable_match_columns(frame, model_id)
    size = len(frame)
    if not columns or size == 0:
        return [[] for _ in range(size)]
    distance = np.zeros((size, size), dtype=float)
    allowed = np.ones((size, size), dtype=bool)
    blocks = frame["calibration_block_id"].astype(str).to_numpy()
    allowed &= blocks[:, None] != blocks[None, :]
    np.fill_diagonal(allowed, False)
    if model_id == "local_genomic_matched" and "chrom" in frame.columns:
        chrom = frame["chrom"].astype(str).to_numpy()
        same_chrom = chrom[:, None] == chrom[None, :]
        for index in range(size):
            if np.any(allowed[index] & same_chrom[index]):
                allowed[index] &= same_chrom[index]
    for column in columns:
        if column in EXACT_MATCH_COLUMNS:
            values = frame[column].fillna("").astype(str).to_numpy()
            same = values[:, None] == values[None, :]
            nonempty = values != ""
            for index in range(size):
                if nonempty[index] and np.any(allowed[index] & same[index]):
                    allowed[index] &= same[index]
            continue
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
        scale = float(np.nanstd(values, ddof=1))
        if np.isfinite(scale) and scale > 0:
            pair_distance = np.abs(values[:, None] - values[None, :]) / scale
            pair_distance[~np.isfinite(pair_distance)] = 0.0
            distance += pair_distance
    controls: list[list[int]] = []
    for index in range(size):
        candidates = np.flatnonzero(allowed[index])
        if candidates.size:
            order = candidates[np.argsort(distance[index, candidates], kind="stable")]
            controls.append(order[: max(1, n_controls)].tolist())
        else:
            controls.append([])
    return controls


def build_severe_null_panel(
    scores: pd.DataFrame,
    features: pd.DataFrame,
    *,
    n_controls: int = 25,
    block_size_bp: int = 100_000,
    minimum_independent_blocks: int = 5,
    agreement_z_threshold: float = 2.0,
    score_columns: Iterable[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if "region_id" not in scores.columns or "region_id" not in features.columns:
        raise ValueError("scores and features must contain region_id")
    frame = features.drop_duplicates("region_id").merge(scores, on="region_id", how="inner", suffixes=("", "_score_table"))
    frame = frame.reset_index(drop=True)
    frame["calibration_block_id"] = infer_genomic_blocks(frame, block_size_bp)
    primitives = list(score_columns) if score_columns is not None else _score_columns(scores)
    availability = assess_null_availability(frame.columns)
    registry_by_id = {str(item["null_model_id"]): item for item in availability}
    detailed_rows: list[dict[str, object]] = []
    for model_id in MATCH_COLUMNS:
        spec = registry_by_id.get(model_id, {})
        columns = _usable_match_columns(frame, model_id)
        model_available = bool(spec.get("available")) or bool(columns)
        if not model_available:
            continue
        control_map = _precompute_control_indices(frame, model_id, n_controls)
        for row_index, row in frame.iterrows():
            control_indices = control_map[row_index]
            independent_blocks = frame.loc[control_indices, "calibration_block_id"].nunique() if control_indices else 0
            for primitive in primitives:
                observed = pd.to_numeric(pd.Series([row.get(primitive)]), errors="coerce").iloc[0]
                values = pd.to_numeric(frame.loc[control_indices, primitive], errors="coerce").dropna().to_numpy(dtype=float) if primitive in frame.columns else np.array([], dtype=float)
                mean = float(np.mean(values)) if values.size else math.nan
                std = float(np.std(values, ddof=1)) if values.size > 1 else math.nan
                zscore = float((observed - mean) / std) if np.isfinite(observed) and np.isfinite(std) and std > 0 else math.nan
                p_value = empirical_p_value(float(observed), values, higher=True) if np.isfinite(observed) and values.size else math.nan
                if values.size < 2:
                    status = "unavailable_insufficient_controls"
                    reason = "Fewer than two controls remained after block-aware matching."
                elif independent_blocks < minimum_independent_blocks:
                    status = "partial_insufficient_independent_blocks"
                    reason = f"Only {independent_blocks} independent blocks; {minimum_independent_blocks} required for promotion."
                elif not np.isfinite(zscore):
                    status = "partial_zero_null_variance"
                    reason = "The matched null distribution had zero or undefined variance."
                else:
                    status = "available_block_calibrated"
                    reason = "Empirical calibration used controls from independent genomic blocks."
                detailed_rows.append(
                    {
                        "region_id": str(row["region_id"]),
                        "primitive": primitive,
                        "null_model_id": model_id,
                        "primitive_score": float(observed) if np.isfinite(observed) else math.nan,
                        "null_mean": mean,
                        "null_std": std,
                        "null_zscore": zscore,
                        "empirical_p_value": float(p_value),
                        "null_sample_size": int(values.size),
                        "independent_block_count": int(independent_blocks),
                        "calibration_block_id": str(row["calibration_block_id"]),
                        "matched_features_used": ",".join(columns),
                        "null_status": status,
                        "null_reason": reason,
                    }
                )
    details = pd.DataFrame(detailed_rows)
    summaries: list[dict[str, object]] = []
    all_ids = [str(item["null_model_id"]) for item in availability]
    for (region_id, primitive), group in details.groupby(["region_id", "primitive"], sort=False):
        usable = group.loc[group["null_sample_size"] >= 2].copy()
        calibrated = usable.loc[usable["null_zscore"].map(np.isfinite)]
        available_ids = usable["null_model_id"].astype(str).tolist()
        missing_ids = [model_id for model_id in all_ids if model_id not in available_ids]
        zvalues = calibrated["null_zscore"].to_numpy(dtype=float)
        pvalues = calibrated["empirical_p_value"].to_numpy(dtype=float)
        survival = zvalues >= agreement_z_threshold
        agreement = float(np.mean(survival)) if survival.size else math.nan
        conflict = bool(survival.any() and (~survival).any()) if survival.size else False
        independent_blocks = int(usable["independent_block_count"].max()) if not usable.empty else 0
        if len(available_ids) >= 3 and independent_blocks >= minimum_independent_blocks:
            panel_status = "severe_null_panel_available"
        elif available_ids:
            panel_status = "partial_null_panel_not_for_promotion"
        else:
            panel_status = "null_panel_unavailable"
        summaries.append(
            {
                "region_id": region_id,
                "primitive": primitive,
                "null_model_id": "severe_null_panel_conservative_aggregate",
                "primitive_score": float(group["primitive_score"].iloc[0]),
                "null_zscore": float(np.min(zvalues)) if zvalues.size else math.nan,
                "empirical_p_value": float(np.max(pvalues[np.isfinite(pvalues)])) if np.isfinite(pvalues).any() else math.nan,
                "empirical_p_value_status": "available_explicit_named_null_panel" if zvalues.size else "unavailable",
                "empirical_p_value_reason": "Conservative maximum empirical p-value across named calibrated nulls." if zvalues.size else "No calibrated named null distribution was available.",
                "empirical_p_value_tail": "right_tail_high_score",
                "null_panel_status": panel_status,
                "available_null_models": ",".join(available_ids),
                "missing_null_models": ",".join(missing_ids),
                "missing_or_partial_null_models": ",".join(missing_ids),
                "null_model_count": int(len(available_ids)),
                "null_model_agreement": agreement,
                "null_model_conflict": conflict,
                "independent_block_count": independent_blocks,
            }
        )
    return pd.DataFrame(summaries), details
