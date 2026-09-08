"""Statistical helpers shared by feature scoring and residualization."""

from __future__ import annotations

import math
from collections.abc import Iterable

import numpy as np


EPS = 1e-12


def safe_divide(numerator: float, denominator: float, default: float = 0.0) -> float:
    if denominator is None or abs(float(denominator)) < EPS:
        return default
    return float(numerator) / float(denominator)


def shannon_entropy(values: Iterable[str] | dict[str, float]) -> float:
    if isinstance(values, dict):
        counts = np.array(list(values.values()), dtype=float)
    else:
        seq = list(values)
        if not seq:
            return 0.0
        _, counts = np.unique(seq, return_counts=True)
        counts = counts.astype(float)
    total = counts.sum()
    if total <= 0:
        return 0.0
    probs = counts / total
    probs = probs[probs > 0]
    return float(-(probs * np.log2(probs)).sum())


def zscore(value: float, values: Iterable[float]) -> float:
    arr = np.array([v for v in values if np.isfinite(v)], dtype=float)
    if arr.size < 2:
        return 0.0
    std = float(arr.std(ddof=1))
    if std < EPS:
        return 0.0
    return float((value - arr.mean()) / std)


def robust_zscore(value: float, values: Iterable[float]) -> float:
    arr = np.array([v for v in values if np.isfinite(v)], dtype=float)
    if arr.size < 2:
        return 0.0
    median = float(np.median(arr))
    mad = float(np.median(np.abs(arr - median)))
    if mad < EPS:
        return zscore(value, arr)
    return float(0.6745 * (value - median) / mad)


def empirical_p_value(value: float, null_values: Iterable[float], higher: bool = True) -> float:
    arr = np.array([v for v in null_values if np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return math.nan
    if higher:
        extreme = (arr >= value).sum()
    else:
        extreme = (arr <= value).sum()
    return float((extreme + 1) / (arr.size + 1))


def minmax01(value: float, lo: float, hi: float) -> float:
    if not np.isfinite(value) or hi <= lo:
        return 0.0
    return float(np.clip((value - lo) / (hi - lo), 0.0, 1.0))


def optional_float(value: object, default: float = math.nan) -> float:
    """Coerce a scalar to a finite float, otherwise NA.

    Missing, non-numeric, and non-finite values stay unavailable. Measured
    zeros are preserved. Callers must not treat the default as evidence.
    """

    if value is None:
        return float(default)
    if isinstance(value, (bytes, str)) and str(value).strip() == "":
        return float(default)
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float(default)
    return number if np.isfinite(number) else float(default)


def optional_row_float(row: dict, *keys: str) -> float:
    """Return the first finite numeric value among *keys*, else NA."""

    getter = getattr(row, "get", None)
    for key in keys:
        if isinstance(row, dict) and key not in row:
            continue
        raw = getter(key) if getter is not None else None
        number = optional_float(raw)
        if np.isfinite(number):
            return number
    return math.nan


def finite_mean(values: Iterable[float]) -> float:
    finite = [float(value) for value in values if np.isfinite(value)]
    return float(np.mean(finite)) if finite else math.nan


def finite_max(values: Iterable[float]) -> float:
    finite = [float(value) for value in values if np.isfinite(value)]
    return float(np.max(finite)) if finite else math.nan


def robust_scale_series(values: Iterable[float]) -> np.ndarray:
    arr = np.array(list(values), dtype=float)
    if arr.size == 0:
        return arr
    finite = np.isfinite(arr)
    out = np.zeros_like(arr, dtype=float)
    if finite.sum() < 2:
        return out
    median = np.nanmedian(arr[finite])
    mad = np.nanmedian(np.abs(arr[finite] - median))
    scale = mad / 0.6745 if mad > EPS else np.nanstd(arr[finite])
    if not np.isfinite(scale) or scale < EPS:
        return out
    out[finite] = (arr[finite] - median) / scale
    return out


def row_score(row: dict, columns: list[str]) -> float:
    """Mean of available columns. Missing columns stay out of the mean; all-missing is NA."""

    return finite_mean(optional_row_float(row, column) for column in columns)
