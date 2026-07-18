"""Validated multiscale texture screens with explicit scaling diagnostics.

The module intentionally avoids claiming that a sequence is fractal.  It only
reports a multiscale texture screening score when a scaling interval, fit
diagnostics, same-composition surrogates, and window-shift stability can be
computed.  Historical ``fractal_*`` fields remain as deprecated aliases.
"""

from __future__ import annotations

import json
import math
import random

import numpy as np
import pandas as pd

from darkdna.features.sequence import NUMERIC_WALK_MAPPINGS, clean_sequence, kmer_counts
from darkdna.utils.optional_deps import optional_import
from darkdna.utils.stats import shannon_entropy


MULTISCALE_METHOD_VERSION = "multiscale_texture_v2"


def simple_slope(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2 or y.size < 2:
        return math.nan
    centered = x - x.mean()
    denominator = float(np.sum(centered**2))
    if denominator <= 0:
        return math.nan
    return float(np.sum(centered * (y - y.mean())) / denominator)


def log_spaced_scales(
    length: int,
    *,
    min_scale: int = 8,
    max_scale_fraction: float = 0.25,
    n_scales: int = 10,
) -> np.ndarray:
    maximum = int(length * max_scale_fraction)
    if maximum < min_scale:
        return np.array([], dtype=int)
    raw = np.geomspace(min_scale, maximum, num=max(2, n_scales))
    scales = sorted({int(round(value)) for value in raw if min_scale <= value <= maximum})
    return np.asarray([scale for scale in scales if length // scale >= 4], dtype=int)


def _segment_variances(profile: np.ndarray, scale: int) -> np.ndarray:
    count = len(profile) // scale
    if count < 2:
        return np.array([], dtype=float)
    x = np.arange(scale, dtype=float)
    x_centered = x - float(x.mean())
    denominator = float(np.sum(x_centered**2))
    variance_blocks: list[np.ndarray] = []
    for source in (profile, profile[::-1]):
        segments = np.asarray(source[: count * scale], dtype=float).reshape(count, scale)
        means = segments.mean(axis=1)
        slopes = ((segments - means[:, None]) * x_centered[None, :]).sum(axis=1) / denominator
        trends = slopes[:, None] * x[None, :] + (means - slopes * float(x.mean()))[:, None]
        variances = np.mean((segments - trends) ** 2, axis=1)
        variance_blocks.append(variances[np.isfinite(variances) & (variances > 0)])
    return np.concatenate(variance_blocks) if variance_blocks else np.array([], dtype=float)


def dfa_fluctuations(values: np.ndarray, scales: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return np.array([], dtype=float)
    profile = np.cumsum(values - np.mean(values))
    fluctuations = []
    for scale in scales:
        variances = _segment_variances(profile, int(scale))
        fluctuations.append(float(np.sqrt(variances.mean())) if variances.size else math.nan)
    return np.asarray(fluctuations, dtype=float)


def _linear_fit_diagnostics(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    slope = simple_slope(x, y)
    if not np.isfinite(slope):
        return {"slope": math.nan, "intercept": math.nan, "r2": math.nan, "rmse": math.nan, "ci_low": math.nan, "ci_high": math.nan}
    intercept = float(y.mean() - slope * x.mean())
    predicted = slope * x + intercept
    residuals = y - predicted
    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = math.nan if ss_tot <= 0 else float(1.0 - ss_res / ss_tot)
    rmse = float(np.sqrt(np.mean(residuals**2)))
    denominator = float(np.sum((x - x.mean()) ** 2))
    if len(x) > 2 and denominator > 0:
        standard_error = math.sqrt(max(0.0, ss_res) / (len(x) - 2) / denominator)
        ci_low, ci_high = slope - 1.96 * standard_error, slope + 1.96 * standard_error
    else:
        ci_low = ci_high = math.nan
    return {"slope": float(slope), "intercept": intercept, "r2": r2, "rmse": rmse, "ci_low": float(ci_low), "ci_high": float(ci_high)}


def identify_scaling_region(
    scales: np.ndarray,
    fluctuations: np.ndarray,
    *,
    minimum_points: int = 4,
    minimum_r2: float = 0.90,
) -> dict[str, float | int | str]:
    valid = np.isfinite(fluctuations) & (fluctuations > 0) & (scales > 0)
    scale_values = scales[valid]
    fluctuation_values = fluctuations[valid]
    if len(scale_values) < minimum_points:
        return {
            "status": "unavailable",
            "reason": f"Only {len(scale_values)} valid DFA scales; at least {minimum_points} are required.",
            "n_points": int(len(scale_values)),
        }
    best: dict[str, float | int | str] | None = None
    for start in range(0, len(scale_values) - minimum_points + 1):
        for stop in range(start + minimum_points, len(scale_values) + 1):
            x = np.log(scale_values[start:stop].astype(float))
            y = np.log(fluctuation_values[start:stop])
            diagnostics = _linear_fit_diagnostics(x, y)
            r2 = float(diagnostics["r2"])
            if not np.isfinite(r2):
                continue
            candidate = {
                **diagnostics,
                "start_scale": int(scale_values[start]),
                "end_scale": int(scale_values[stop - 1]),
                "n_points": int(stop - start),
            }
            # Prefer longer intervals, then better fits.  This makes a perfect
            # four-point subrange less attractive than a stable broad interval.
            merit = r2 + 0.015 * (stop - start)
            if best is None or merit > float(best["merit"]):
                best = {**candidate, "merit": merit}
    if best is None:
        return {"status": "unavailable", "reason": "No finite log-log scaling fit was found.", "n_points": int(len(scale_values))}
    best.pop("merit", None)
    if float(best["r2"]) < minimum_r2:
        return {**best, "status": "unavailable", "reason": f"Best scaling fit R2={best['r2']:.3f} is below {minimum_r2:.2f}."}
    return {**best, "status": "available", "reason": "A diagnostic scaling interval met the configured fit threshold."}


def dfa_diagnostics(values: np.ndarray, **scale_options: object) -> dict[str, object]:
    scales = log_spaced_scales(len(values), **scale_options)
    fluctuations = dfa_fluctuations(values, scales)
    fit = identify_scaling_region(scales, fluctuations)
    return {
        **fit,
        "scales": scales.tolist(),
        "fluctuations": [float(value) if np.isfinite(value) else None for value in fluctuations],
    }


def dfa_like_exponent(values: np.ndarray) -> float:
    """Deprecated compatibility wrapper for the previous three-scale estimator."""

    diagnostics = dfa_diagnostics(np.asarray(values, dtype=float))
    return float(diagnostics.get("slope", math.nan))


def multifractal_dfa(values: np.ndarray, scales: np.ndarray, q_values: tuple[int, ...] = (-4, -2, 0, 2, 4)) -> dict[str, object]:
    if len(values) < 512 or len(scales) < 5:
        return {
            "status": "unavailable",
            "reason": "Multifractal DFA requires at least 512 bases and five valid scales.",
            "spectrum_width": math.nan,
            "generalized_hurst": {},
        }
    profile = np.cumsum(np.asarray(values, dtype=float) - np.mean(values))
    exponents: dict[str, float] = {}
    for q in q_values:
        fq: list[float] = []
        valid_scales: list[int] = []
        for scale in scales:
            variances = _segment_variances(profile, int(scale))
            if variances.size == 0:
                continue
            if q == 0:
                value = float(np.exp(0.5 * np.mean(np.log(variances))))
            else:
                value = float(np.mean(variances ** (q / 2.0)) ** (1.0 / q))
            if value > 0 and np.isfinite(value):
                fq.append(value)
                valid_scales.append(int(scale))
        if len(fq) >= 4:
            fit = _linear_fit_diagnostics(np.log(np.asarray(valid_scales, dtype=float)), np.log(np.asarray(fq, dtype=float)))
            if np.isfinite(fit["slope"]):
                exponents[str(q)] = float(fit["slope"])
    if len(exponents) < 3:
        return {
            "status": "unavailable",
            "reason": "Too few generalized Hurst exponents were estimable.",
            "spectrum_width": math.nan,
            "generalized_hurst": exponents,
        }
    width = float(max(exponents.values()) - min(exponents.values()))
    return {"status": "available", "reason": "Generalized exponents were estimated across q values.", "spectrum_width": width, "generalized_hurst": exponents}


def fourier_spectrum_summaries(values: np.ndarray) -> dict[str, float]:
    if values.size < 8:
        return {"spectrum_slope": math.nan, "low_frequency_power": math.nan, "spectral_entropy": math.nan}
    spectrum = np.abs(np.fft.rfft(values - values.mean())) ** 2
    powers = spectrum[1:]
    if powers.size == 0 or powers.sum() <= 0:
        return {"spectrum_slope": math.nan, "low_frequency_power": math.nan, "spectral_entropy": math.nan}
    frequencies = np.arange(1, powers.size + 1, dtype=float)
    positive = powers > 0
    slope = simple_slope(np.log(frequencies[positive]), np.log(powers[positive])) if positive.sum() > 2 else math.nan
    probabilities = powers / powers.sum()
    entropy = float(-(probabilities[probabilities > 0] * np.log2(probabilities[probabilities > 0])).sum())
    return {
        "spectrum_slope": float(slope),
        "low_frequency_power": float(powers[: max(1, powers.size // 10)].sum() / powers.sum()),
        "spectral_entropy": entropy,
    }


def wavelet_spectrum_summaries(values: np.ndarray) -> dict[str, float | str]:
    pywt, warning = optional_import("pywt")
    if pywt is None:
        return {
            "wavelet_spectrum_status": "unavailable",
            "wavelet_spectrum_reason": warning or "PyWavelets is unavailable.",
            "wavelet_energy_total": math.nan,
            "wavelet_energy_entropy": math.nan,
            "wavelet_logscale_slope": math.nan,
        }
    if values.size < 32:
        return {
            "wavelet_spectrum_status": "unavailable",
            "wavelet_spectrum_reason": "At least 32 values are required for a multilevel wavelet spectrum.",
            "wavelet_energy_total": math.nan,
            "wavelet_energy_entropy": math.nan,
            "wavelet_logscale_slope": math.nan,
        }
    level = min(6, int(np.log2(values.size)) - 2)
    coefficients = pywt.wavedec(values, "db2", mode="periodization", level=level)
    detail = coefficients[1:]
    energies = np.asarray([float(np.mean(coefficient**2)) for coefficient in detail], dtype=float)
    total = float(energies.sum())
    probabilities = energies / total if total > 0 else np.full(len(energies), math.nan)
    entropy = float(-(probabilities[probabilities > 0] * np.log2(probabilities[probabilities > 0])).sum()) if total > 0 else math.nan
    positive = energies > 0
    slope = simple_slope(np.arange(1, len(energies) + 1, dtype=float)[positive], np.log2(energies[positive])) if positive.sum() >= 3 else math.nan
    return {
        "wavelet_spectrum_status": "available",
        "wavelet_spectrum_reason": "Computed a db2 multilevel detail-energy spectrum.",
        "wavelet_energy_total": total,
        "wavelet_energy_entropy": entropy,
        "wavelet_logscale_slope": float(slope),
    }


def wavelet_energy(values: np.ndarray) -> dict[str, float | str]:
    """Compatibility alias for the former wavelet helper."""

    return wavelet_spectrum_summaries(values)


def chaos_game_features(seq: str, k: int) -> dict[str, float]:
    seq = clean_sequence(seq)
    if len(seq) < k:
        return {f"chaos_game_k{k}_occupancy": math.nan, f"chaos_game_k{k}_entropy": math.nan}
    counts = kmer_counts(seq, k)
    return {
        f"chaos_game_k{k}_occupancy": float(len(counts) / (4**k)),
        f"chaos_game_k{k}_entropy": shannon_entropy(counts),
    }


def _mapping_steps(seq: str, mapping: dict[str, float]) -> np.ndarray:
    return np.asarray([mapping.get(base, 0.0) for base in seq], dtype=float)


def _surrogate_dfa_test(values: np.ndarray, observed_slope: float, *, n_surrogates: int, seed: int) -> dict[str, float | int | str]:
    if not np.isfinite(observed_slope) or n_surrogates < 4:
        return {"status": "unavailable", "reason": "A valid observed slope and at least four surrogates are required.", "zscore": math.nan, "empirical_p_value": math.nan, "n_valid": 0}
    rng = np.random.default_rng(seed)
    slopes: list[float] = []
    for _ in range(n_surrogates):
        shuffled = rng.permutation(values)
        diagnostic = dfa_diagnostics(shuffled)
        slope = float(diagnostic.get("slope", math.nan))
        if diagnostic.get("status") == "available" and np.isfinite(slope):
            slopes.append(slope)
    if len(slopes) < 4:
        return {"status": "unavailable", "reason": "Fewer than four surrogate sequences had valid scaling fits.", "zscore": math.nan, "empirical_p_value": math.nan, "n_valid": len(slopes)}
    null = np.asarray(slopes, dtype=float)
    standard_deviation = float(null.std(ddof=1))
    zscore = (observed_slope - float(null.mean())) / standard_deviation if standard_deviation > 1e-12 else math.nan
    p_value = float((np.sum(np.abs(null - null.mean()) >= abs(observed_slope - null.mean())) + 1) / (len(null) + 1))
    return {
        "status": "available" if np.isfinite(zscore) else "unavailable",
        "reason": "Mononucleotide-composition-preserving permutation test." if np.isfinite(zscore) else "Surrogate slope variance was zero.",
        "zscore": float(zscore),
        "empirical_p_value": p_value if np.isfinite(zscore) else math.nan,
        "n_valid": int(len(null)),
        "null_mean": float(null.mean()),
        "null_std": standard_deviation,
    }


def _window_shift_stability(values: np.ndarray) -> dict[str, float | int | str]:
    slopes: list[float] = []
    max_shift = max(1, len(values) // 50)
    for offset in (0, max_shift, 2 * max_shift, 3 * max_shift):
        shifted = values[offset:]
        diagnostics = dfa_diagnostics(shifted)
        slope = float(diagnostics.get("slope", math.nan))
        if diagnostics.get("status") == "available" and np.isfinite(slope):
            slopes.append(slope)
    if len(slopes) < 3:
        return {"status": "unavailable", "reason": "Fewer than three shifted windows had valid scaling fits.", "slope_sd": math.nan, "n_valid": len(slopes)}
    return {"status": "available", "reason": "DFA slope stability across four nearby start offsets.", "slope_sd": float(np.std(slopes, ddof=1)), "n_valid": len(slopes)}


def compute_scale_fractal_features(
    seq: str,
    *,
    n_surrogates: int = 4,
    seed: int = 13,
) -> dict[str, float | int | str]:
    seq = clean_sequence(seq)
    mapping_diagnostics: dict[str, dict[str, object]] = {}
    primary_name = "strong_weak_H_bond_numeric_walk"
    for name, definition in NUMERIC_WALK_MAPPINGS.items():
        steps = _mapping_steps(seq, definition["mapping"])
        mapping_diagnostics[name] = dfa_diagnostics(steps)

    primary = mapping_diagnostics[primary_name]
    primary_steps = _mapping_steps(seq, NUMERIC_WALK_MAPPINGS[primary_name]["mapping"])
    primary_slope = float(primary.get("slope", math.nan))
    surrogate = _surrogate_dfa_test(primary_steps, primary_slope, n_surrogates=n_surrogates, seed=seed)
    stability = _window_shift_stability(primary_steps)
    scales = np.asarray(primary.get("scales", []), dtype=int)
    multifractal = multifractal_dfa(primary_steps, scales)
    fourier = fourier_spectrum_summaries(primary_steps)
    wavelet = wavelet_spectrum_summaries(primary_steps)

    if primary.get("status") == "available" and surrogate.get("status") == "available" and stability.get("status") == "available":
        score = abs(float(surrogate["zscore"])) * max(0.0, float(primary.get("r2", 0.0))) / (1.0 + float(stability["slope_sd"]))
        status = "available"
        reason = "Valid scaling interval, surrogate calibration, and window-shift diagnostics were available."
    else:
        score = math.nan
        status = "unavailable"
        missing = [
            name
            for name, diagnostic in (("scaling_interval", primary), ("surrogate_test", surrogate), ("window_shift_stability", stability))
            if diagnostic.get("status") != "available"
        ]
        reason = "Required multiscale diagnostics unavailable: " + ", ".join(missing)

    chaos: dict[str, float] = {}
    for k in (3, 4, 5, 6):
        chaos.update(chaos_game_features(seq, k))

    deprecation = (
        "fractal_score is a deprecated alias for multiscale_texture_screening_score. "
        "The screen does not establish fractality or multifractality."
    )
    return {
        "multiscale_method": MULTISCALE_METHOD_VERSION,
        "multiscale_texture_status": status,
        "multiscale_texture_reason": reason,
        "multiscale_texture_screening_score": float(score),
        "multiscale_primary_mapping": primary_name,
        "multiscale_mapping_diagnostics": json.dumps(mapping_diagnostics, sort_keys=True),
        "DFA_exponent": primary_slope,
        "DFA_scaling_r2": float(primary.get("r2", math.nan)),
        "DFA_scaling_residual_rmse": float(primary.get("rmse", math.nan)),
        "DFA_exponent_ci_low": float(primary.get("ci_low", math.nan)),
        "DFA_exponent_ci_high": float(primary.get("ci_high", math.nan)),
        "DFA_scaling_start": float(primary.get("start_scale", math.nan)),
        "DFA_scaling_end": float(primary.get("end_scale", math.nan)),
        "DFA_valid_scale_count": int(primary.get("n_points", 0)),
        "DFA_surrogate_zscore": float(surrogate.get("zscore", math.nan)),
        "DFA_surrogate_empirical_p_value": float(surrogate.get("empirical_p_value", math.nan)),
        "DFA_surrogate_status": str(surrogate.get("status", "unavailable")),
        "DFA_surrogate_count": int(surrogate.get("n_valid", 0)),
        "DFA_window_shift_slope_sd": float(stability.get("slope_sd", math.nan)),
        "DFA_window_shift_status": str(stability.get("status", "unavailable")),
        "multifractal_DFA_status": str(multifractal["status"]),
        "multifractal_DFA_reason": str(multifractal["reason"]),
        "multifractal_spectrum_width": float(multifractal["spectrum_width"]),
        "multifractal_generalized_hurst": json.dumps(multifractal["generalized_hurst"], sort_keys=True),
        "fractal_score": float(score),
        "fractal_score_status": "deprecated_alias",
        "fractal_score_deprecation_warning": deprecation,
        "DFA_like_exponent_estimator": primary_slope,
        "Hurst_like_exponent_estimator": primary_slope,
        "DFA_like_exponent_deprecation_warning": "Use DFA_exponent with DFA scaling diagnostics; it is not a Hurst estimate without validated assumptions.",
        "numeric_walk_generation": json.dumps(np.cumsum(primary_steps)[:200].round(4).tolist()),
        "scale_persistence_score": math.nan,
        "parent_child_window_consistency": math.nan,
        "scale_persistence_deprecation_warning": "Parent-child similarity is not scale persistence; use multiscale_parent_child_similarity_screen.",
        "cross_scale_feature_correlation": math.nan,
        "scale_breakpoint_detection": math.nan,
        "micro_to_macro_gain": math.nan,
        "renormalization_profile": math.nan,
        "renormalization_profile_status": "unavailable_not_implemented_without_validated_renormalization_model",
        **fourier,
        **wavelet,
        **chaos,
    }


def add_multiscale_profiles(feature_table: pd.DataFrame, windows: pd.DataFrame | None = None) -> pd.DataFrame:
    out = feature_table.copy()
    if windows is not None and "region_id" in windows.columns:
        context_columns = [column for column in ("region_id", "parent_region_id", "child_region_ids", "window_size") if column in windows.columns]
        missing_context = [column for column in context_columns if column != "region_id" and column not in out.columns]
        if missing_context:
            out = out.merge(windows[["region_id", *missing_context]], on="region_id", how="left")
    score_column = "multiscale_texture_screening_score"
    if score_column not in out.columns:
        return out
    out["multiscale_parent_child_similarity_screen"] = math.nan
    out["multiscale_parent_child_similarity_status"] = "unavailable_no_related_valid_scores"
    by_region = out.set_index("region_id")
    for index, row in out.iterrows():
        related: list[float] = []
        parent = row.get("parent_region_id")
        if isinstance(parent, str) and parent in by_region.index:
            parent_value = pd.to_numeric(pd.Series([by_region.loc[parent, score_column]]), errors="coerce").iloc[0]
            if np.isfinite(parent_value):
                related.append(float(parent_value))
        children = str(row.get("child_region_ids") or "").split(",")
        for child in children:
            if child and child in by_region.index:
                child_value = pd.to_numeric(pd.Series([by_region.loc[child, score_column]]), errors="coerce").iloc[0]
                if np.isfinite(child_value):
                    related.append(float(child_value))
        score = pd.to_numeric(pd.Series([row.get(score_column)]), errors="coerce").iloc[0]
        if related and np.isfinite(score):
            similarities = [1.0 / (1.0 + abs(float(score) - value)) for value in related]
            similarity = float(np.mean(similarities))
            out.loc[index, "multiscale_parent_child_similarity_screen"] = similarity
            out.loc[index, "multiscale_parent_child_similarity_status"] = "available_descriptive_not_scaling_evidence"
            out.loc[index, "micro_to_macro_gain"] = float(score) - float(np.mean(related))
    out["scale_persistence_score"] = out["multiscale_parent_child_similarity_screen"]
    out["parent_child_window_consistency"] = out["multiscale_parent_child_similarity_screen"]
    return out
