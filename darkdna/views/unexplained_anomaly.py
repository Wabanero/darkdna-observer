"""Block-aware multivariate outlier screening for unexplained anomalies."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


def _robust_location_scale(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    values = frame.to_numpy(dtype=float)
    location = np.nanmedian(values, axis=0)
    mad = np.nanmedian(np.abs(values - location), axis=0) / 0.6745
    standard = np.nanstd(values, axis=0, ddof=1)
    scale = np.where(np.isfinite(mad) & (mad > 1e-9), mad, standard)
    scale = np.where(np.isfinite(scale) & (scale > 1e-9), scale, 1.0)
    return location, scale


def _mahalanobis_from_training(train: pd.DataFrame, test: pd.DataFrame) -> np.ndarray:
    location, scale = _robust_location_scale(train)
    train_z = (train.to_numpy(dtype=float) - location) / scale
    test_z = (test.to_numpy(dtype=float) - location) / scale
    train_z = np.nan_to_num(train_z, nan=0.0, posinf=0.0, neginf=0.0)
    test_z = np.nan_to_num(test_z, nan=0.0, posinf=0.0, neginf=0.0)
    if train_z.shape[1] == 1:
        return np.abs(test_z[:, 0])
    feature_count = train_z.shape[1]
    means = [float(np.mean(train_z[:, index])) for index in range(feature_count)]
    denominator = max(1, len(train_z) - 1)
    covariance = [[0.0 for _ in range(feature_count)] for _ in range(feature_count)]
    for left in range(feature_count):
        for right in range(feature_count):
            value = sum(
                (float(row[left]) - means[left]) * (float(row[right]) - means[right])
                for row in train_z
            ) / denominator
            covariance[left][right] = 0.8 * value + (0.2 if left == right else 0.0)
    inverse = _gauss_jordan_inverse(covariance)
    distances: list[float] = []
    for row in test_z:
        vector = [float(value) for value in row]
        transformed = [sum(inverse[left][right] * vector[right] for right in range(feature_count)) for left in range(feature_count)]
        squared = sum(vector[index] * transformed[index] for index in range(feature_count))
        distances.append(math.sqrt(max(0.0, squared)))
    return np.asarray(distances, dtype=float)


def _gauss_jordan_inverse(matrix: list[list[float]]) -> list[list[float]]:
    """Invert a small shrinkage covariance matrix without native LAPACK calls."""

    size = len(matrix)
    augmented = [
        [float(matrix[row][column]) for column in range(size)]
        + [1.0 if row == column else 0.0 for column in range(size)]
        for row in range(size)
    ]
    for column in range(size):
        pivot_row = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot_row][column]) < 1e-12:
            augmented[pivot_row][column] += 1e-6
        augmented[column], augmented[pivot_row] = augmented[pivot_row], augmented[column]
        pivot = augmented[column][column]
        augmented[column] = [value / pivot for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                augmented[row][index] - factor * augmented[column][index]
                for index in range(2 * size)
            ]
    return [row[size:] for row in augmented]


def cross_fitted_unexplained_outlierness(
    feature_matrix: pd.DataFrame,
    groups: pd.Series | None,
    *,
    minimum_training_rows: int = 8,
) -> tuple[pd.Series, dict[str, object]]:
    """Score held-out multivariate outlierness without random row leakage.

    ``groups`` must encode chromosomes or non-overlapping genomic blocks.  The
    function refuses to random-split rows because overlapping windows are not
    independent evidence.
    """

    numeric = feature_matrix.apply(pd.to_numeric, errors="coerce")
    informative = [
        column
        for column in numeric.columns
        if int(numeric[column].notna().sum()) >= 3 and float(numeric[column].std(skipna=True) or 0.0) > 1e-12
    ]
    output = pd.Series(np.nan, index=feature_matrix.index, dtype=float)
    if len(informative) < 2:
        return output, {
            "status": "unavailable",
            "method": "cross_fitted_shrinkage_robust_mahalanobis",
            "reason": "At least two non-degenerate feature axes are required.",
            "feature_count": len(informative),
        }
    if groups is None:
        return output, {
            "status": "unavailable",
            "method": "cross_fitted_shrinkage_robust_mahalanobis",
            "reason": "No chromosome or non-overlapping genomic-block grouping was supplied; random row splitting is forbidden.",
            "feature_count": len(informative),
        }
    group_values = groups.reindex(feature_matrix.index).fillna("missing").astype(str)
    if group_values.nunique() < 2:
        return output, {
            "status": "unavailable",
            "method": "cross_fitted_shrinkage_robust_mahalanobis",
            "reason": "At least two independent chromosome/block groups are required.",
            "feature_count": len(informative),
        }
    used_groups = 0
    matrix = numeric[informative]
    for group in group_values.drop_duplicates():
        test_mask = group_values == group
        train_mask = ~test_mask
        if int(train_mask.sum()) < minimum_training_rows:
            continue
        train = matrix.loc[train_mask]
        test = matrix.loc[test_mask]
        distances = _mahalanobis_from_training(train, test)
        # Calibrate distance only against the corresponding training data.
        train_distances = _mahalanobis_from_training(train, train)
        median = float(np.median(train_distances))
        mad = float(np.median(np.abs(train_distances - median)) / 0.6745)
        if not np.isfinite(mad) or mad <= 1e-9:
            continue
        output.loc[test_mask] = (distances - median) / mad
        used_groups += 1
    status = "available" if int(output.notna().sum()) == len(output) else "partial"
    if int(output.notna().sum()) == 0:
        status = "unavailable"
    return output, {
        "status": status,
        "method": "cross_fitted_shrinkage_robust_mahalanobis",
        "reason": (
            "Held-out distances were calibrated against training-block distances."
            if status == "available"
            else "Some or all folds lacked enough independent training rows or non-degenerate calibration distances."
        ),
        "feature_count": len(informative),
        "features": informative,
        "group_count": int(group_values.nunique()),
        "groups_scored": int(used_groups),
        "score_interpretation": "held_out_multivariate_outlier_screen_not_null_significance_or_function",
    }


def infer_block_groups(features: pd.DataFrame, *, block_size: int = 100_000) -> pd.Series | None:
    if "block_id" in features.columns and features["block_id"].fillna("").astype(str).nunique() >= 2:
        return features["block_id"].astype(str)
    if {"chrom", "start"}.issubset(features.columns):
        chrom = features["chrom"].fillna("missing").astype(str)
        start = pd.to_numeric(features["start"], errors="coerce")
        blocks = chrom + ":block_" + (start.fillna(-1).astype(int) // block_size).astype(str)
        if blocks.nunique() >= 2:
            return blocks
    if "chrom" in features.columns and features["chrom"].fillna("").astype(str).nunique() >= 2:
        return features["chrom"].astype(str)
    return None
