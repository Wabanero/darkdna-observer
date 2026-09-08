"""Entropy/noise sequence-first view."""

from __future__ import annotations

from darkdna.utils.stats import finite_mean, optional_row_float


def compute_entropy_noise_view(row: dict) -> dict[str, float]:
    local_entropy_variance = optional_row_float(row, "local_entropy_cliffs")
    entropy_cliff = optional_row_float(row, "entropy_boundary_score")
    low_complexity_boundary = optional_row_float(row, "repeat_boundary_score")
    compression = optional_row_float(row, "compression_anomaly_score", "compression_boundary_score")
    entropy_asym = optional_row_float(row, "left_right_entropy_asymmetry")
    void = optional_row_float(row, "local_feature_void_score")
    decoherence = finite_mean(
        [entropy_cliff, low_complexity_boundary, compression, entropy_asym, void]
    )
    return {
        "local_entropy_variance_across_subwindows": local_entropy_variance,
        "entropy_cliff_score": entropy_cliff,
        "low_complexity_boundary_score": low_complexity_boundary,
        "entropy_noise_compression_anomaly_score": compression,
        "entropy_asymmetry_score": entropy_asym,
        "feature_void_score": void,
        "entropy_noise_boundary_score": decoherence,
    }
