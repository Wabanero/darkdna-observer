"""Copy-number features with grouped validation when phenotype data exist."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from darkdna.architecture.amount_features import _overlap_bp


def _sample_copy_number(interval: pd.Series, tracks: pd.DataFrame) -> list[dict[str, object]]:
    same = tracks.loc[tracks["chrom"].astype(str) == str(interval["chrom"])].copy()
    same = same.loc[(pd.to_numeric(same["end"]) > int(interval["start"])) & (pd.to_numeric(same["start"]) < int(interval["end"]))]
    if same.empty:
        return []
    value_column = "copy_number" if "copy_number" in same.columns else "value"
    sample_column = "sample_id" if "sample_id" in same.columns else None
    rows = []
    for sample_id, group in same.groupby(sample_column, dropna=False) if sample_column else [("reference", same)]:
        weights = np.asarray([
            _overlap_bp(int(interval["start"]), int(interval["end"]), int(row.start), int(row.end))
            for row in group.itertuples()
        ], dtype=float)
        values = pd.to_numeric(group[value_column], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(values) & (weights > 0)
        if valid.any():
            rows.append({"sample_id": str(sample_id), "copy_number": float(np.average(values[valid], weights=weights[valid]))})
    return rows


def _grouped_dosage_effect(values: pd.DataFrame, phenotype: pd.DataFrame | None) -> tuple[float, float, str, str]:
    if phenotype is None or phenotype.empty or not {"sample_id", "phenotype"}.issubset(phenotype.columns):
        return math.nan, math.nan, "unavailable", "A phenotype table with sample_id and phenotype is required."
    merged = values.merge(phenotype, on="sample_id", how="inner")
    group_column = next((column for column in ("group_id", "haplotype", "strain", "accession", "family_id") if column in merged.columns), None)
    if group_column is None or merged[group_column].nunique() < 2:
        return math.nan, math.nan, "unavailable_grouped_cv", "At least two independent groups are required; random sample splitting is not used."
    slopes, errors = [], []
    for held_out in merged[group_column].dropna().unique():
        train = merged.loc[merged[group_column] != held_out]
        test = merged.loc[merged[group_column] == held_out]
        x_train = pd.to_numeric(train["copy_number"], errors="coerce").to_numpy(dtype=float)
        y_train = pd.to_numeric(train["phenotype"], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(x_train) & np.isfinite(y_train)
        if valid.sum() < 2 or np.std(x_train[valid]) == 0:
            continue
        x_valid = x_train[valid]
        y_valid = y_train[valid]
        x_mean = float(np.mean(x_valid))
        y_mean = float(np.mean(y_valid))
        denominator = float(np.sum((x_valid - x_mean) ** 2))
        if denominator == 0:
            continue
        slope = float(np.sum((x_valid - x_mean) * (y_valid - y_mean)) / denominator)
        intercept = y_mean - slope * x_mean
        x_test = pd.to_numeric(test["copy_number"], errors="coerce").to_numpy(dtype=float)
        y_test = pd.to_numeric(test["phenotype"], errors="coerce").to_numpy(dtype=float)
        test_valid = np.isfinite(x_test) & np.isfinite(y_test)
        if test_valid.any():
            slopes.append(float(slope))
            errors.extend((y_test[test_valid] - (intercept + slope * x_test[test_valid])).tolist())
    if not slopes:
        return math.nan, math.nan, "unavailable_grouped_cv", "No held-out group had an estimable training slope."
    baseline = float(pd.to_numeric(merged["phenotype"], errors="coerce").var())
    mse = float(np.mean(np.square(errors))) if errors else math.nan
    cv_score = 1.0 - mse / baseline if np.isfinite(baseline) and baseline > 0 and np.isfinite(mse) else math.nan
    return float(np.mean(slopes)), cv_score, "available_grouped_cv", f"Leave-one-{group_column}-out validation."


def compute_copy_number_features(
    intervals: pd.DataFrame,
    copy_number_tracks: pd.DataFrame | None = None,
    phenotype: pd.DataFrame | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, interval in intervals.iterrows():
        region_id = str(interval["region_id"])
        if copy_number_tracks is None or copy_number_tracks.empty or not {"chrom", "start", "end"}.issubset(copy_number_tracks.columns):
            rows.append(
                {
                    "region_id": region_id,
                    "sample_specific_copy_number": "",
                    "copy_number_mean": math.nan,
                    "copy_number_variance": math.nan,
                    "deletion_frequency": math.nan,
                    "duplication_frequency": math.nan,
                    "copy_number_effect_size": math.nan,
                    "copy_number_cv_score": math.nan,
                    "dosage_response_shape": "unavailable",
                    "copy_number_status": "unavailable_missing_copy_number_input",
                    "copy_number_reason": "No bedGraph/CNV copy-number table was supplied.",
                }
            )
            continue
        per_sample = pd.DataFrame(_sample_copy_number(interval, copy_number_tracks))
        if per_sample.empty:
            values = np.array([], dtype=float)
        else:
            values = per_sample["copy_number"].to_numpy(dtype=float)
        effect, cv_score, cv_status, cv_reason = _grouped_dosage_effect(per_sample, phenotype)
        rows.append(
            {
                "region_id": region_id,
                "sample_specific_copy_number": ";".join(f"{row.sample_id}:{row.copy_number:.6g}" for row in per_sample.itertuples()),
                "copy_number_mean": float(np.mean(values)) if values.size else math.nan,
                "copy_number_variance": float(np.var(values, ddof=1)) if values.size > 1 else math.nan,
                "deletion_frequency": float(np.mean(values < 1.5)) if values.size else math.nan,
                "duplication_frequency": float(np.mean(values > 2.5)) if values.size else math.nan,
                "copy_number_effect_size": effect,
                "copy_number_cv_score": cv_score,
                "dosage_response_shape": "linear_screen" if np.isfinite(effect) else "unavailable",
                "copy_number_status": "available" if values.size else "unavailable_no_overlapping_measurements",
                "copy_number_cv_status": cv_status,
                "copy_number_reason": cv_reason,
                "copy_number_sample_size": int(values.size),
            }
        )
    return pd.DataFrame(rows)
