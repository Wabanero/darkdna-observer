"""Scale/fractal view features."""

from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd

from darkdna.utils.optional_deps import optional_import
from darkdna.utils.stats import safe_divide
from darkdna.features.sequence import clean_sequence, compression_ratio, numeric_walk


def simple_slope(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2 or y.size < 2:
        return 0.0
    x_centered = x - x.mean()
    denom = float(np.sum(x_centered**2))
    if denom == 0:
        return 0.0
    return float(np.sum(x_centered * (y - y.mean())) / denom)


def dfa_like_exponent(values: np.ndarray) -> float:
    if values.size < 16:
        return math.nan
    scales = np.array([s for s in [4, 8, 16, 32, 64, 128] if s < values.size // 2])
    if scales.size < 2:
        return math.nan
    fluctuations = []
    x = np.arange(values.size)
    profile = values - np.mean(values)
    for scale in scales:
        chunks = values.size // scale
        rms = []
        for idx in range(chunks):
            y = profile[idx * scale : (idx + 1) * scale]
            xx = x[:scale]
            slope = simple_slope(xx, y)
            intercept = float(y.mean() - slope * xx.mean())
            trend = slope * xx + intercept
            rms.append(np.sqrt(np.mean((y - trend) ** 2)))
        fluctuations.append(np.mean(rms))
    if any(f <= 0 for f in fluctuations):
        return math.nan
    slope = simple_slope(np.log(scales), np.log(fluctuations))
    return float(slope)


def fourier_spectrum_summaries(values: np.ndarray) -> dict[str, float]:
    if values.size < 4:
        return {"spectrum_slope": math.nan, "low_frequency_power": 0.0, "spectral_entropy": 0.0}
    spectrum = np.abs(np.fft.rfft(values - values.mean())) ** 2
    powers = spectrum[1:]
    if powers.size == 0 or powers.sum() <= 0:
        return {"spectrum_slope": math.nan, "low_frequency_power": 0.0, "spectral_entropy": 0.0}
    freqs = np.arange(1, powers.size + 1)
    positive = powers > 0
    slope = simple_slope(np.log(freqs[positive]), np.log(powers[positive])) if positive.sum() > 1 else math.nan
    probs = powers / powers.sum()
    entropy = float(-(probs[probs > 0] * np.log2(probs[probs > 0])).sum())
    return {
        "spectrum_slope": float(slope),
        "low_frequency_power": float(powers[: max(1, powers.size // 10)].sum() / powers.sum()),
        "spectral_entropy": entropy,
    }


def wavelet_energy(values: np.ndarray) -> dict[str, float]:
    pywt, warning = optional_import("pywt")
    if pywt is None or values.size < 8:
        return {"wavelet_energy_total": math.nan, "wavelet_energy_entropy": math.nan}
    coeffs = pywt.wavedec(values, "db1", mode="periodization", level=min(4, int(np.log2(values.size))))
    energies = np.array([float(np.sum(c**2)) for c in coeffs])
    total = energies.sum()
    probs = energies / total if total else energies
    entropy = float(-(probs[probs > 0] * np.log2(probs[probs > 0])).sum()) if total else math.nan
    return {"wavelet_energy_total": float(total), "wavelet_energy_entropy": entropy}


def chaos_game_features(seq: str, k: int) -> dict[str, float]:
    seq = clean_sequence(seq)
    if len(seq) < k:
        return {f"chaos_game_k{k}_occupancy": 0.0, f"chaos_game_k{k}_entropy": 0.0}
    from darkdna.features.sequence import kmer_counts
    from darkdna.utils.stats import shannon_entropy

    counts = kmer_counts(seq, k)
    occupancy = len(counts) / (4**k)
    return {f"chaos_game_k{k}_occupancy": float(occupancy), f"chaos_game_k{k}_entropy": shannon_entropy(counts)}


def compute_scale_fractal_features(seq: str) -> dict[str, float | str]:
    seq = clean_sequence(seq)
    walk = numeric_walk(seq, {"G": 1, "C": 1, "A": -1, "T": -1})
    fourier = fourier_spectrum_summaries(walk)
    wavelet = wavelet_energy(walk)
    dfa = dfa_like_exponent(walk)
    chaos = {}
    for k in [3, 4, 5, 6]:
        chaos.update(chaos_game_features(seq, k))
    fractal_score = float(np.nanmean([abs(fourier.get("spectrum_slope", 0.0)), fourier["low_frequency_power"], wavelet.get("wavelet_energy_entropy", 0.0) or 0.0, dfa if np.isfinite(dfa) else 0.0]))
    compression_anomaly = compression_ratio(seq) - 1.0
    return {
        "numeric_walk_generation": json.dumps(walk[:200].round(4).tolist()),
        "DFA_like_exponent_estimator": dfa,
        "Hurst_like_exponent_estimator": dfa,
        "fractal_score": fractal_score,
        "compression_anomaly_score": float(compression_anomaly),
        "scale_persistence_score": 0.0,
        "parent_child_window_consistency": 0.0,
        "cross_scale_feature_correlation": math.nan,
        "scale_breakpoint_detection": 0.0,
        "micro_to_macro_gain": 0.0,
        "renormalization_profile": json.dumps({}),
        **fourier,
        **wavelet,
        **chaos,
    }


def add_multiscale_profiles(feature_table: pd.DataFrame, windows: pd.DataFrame | None = None) -> pd.DataFrame:
    out = feature_table.copy()
    if windows is not None and "region_id" in windows.columns:
        context = windows[["region_id", "parent_region_id", "child_region_ids", "window_size"]].copy()
        out = out.merge(context, on="region_id", how="left")
    if "fractal_score" not in out.columns:
        return out
    out["scale_persistence_score"] = 0.0
    out["parent_child_window_consistency"] = 0.0
    by_region = out.set_index("region_id")
    for idx, row in out.iterrows():
        related = []
        parent = row.get("parent_region_id")
        if isinstance(parent, str) and parent in by_region.index:
            related.append(float(by_region.loc[parent, "fractal_score"]))
        children = str(row.get("child_region_ids") or "").split(",")
        for child in children:
            if child and child in by_region.index:
                related.append(float(by_region.loc[child, "fractal_score"]))
        if related:
            score = float(row["fractal_score"])
            diffs = [1.0 / (1.0 + abs(score - r)) for r in related]
            out.loc[idx, "scale_persistence_score"] = float(np.mean(diffs))
            out.loc[idx, "parent_child_window_consistency"] = float(np.mean(diffs))
            out.loc[idx, "micro_to_macro_gain"] = score - float(np.mean(related))
    return out
