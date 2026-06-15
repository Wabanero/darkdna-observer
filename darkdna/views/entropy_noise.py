"""Entropy/noise sequence-first view."""

from __future__ import annotations

import numpy as np


def compute_entropy_noise_view(row: dict) -> dict[str, float]:
    local_entropy_variance = float(row.get("local_entropy_cliffs", 0.0) or 0.0)
    entropy_cliff = float(row.get("entropy_boundary_score", 0.0) or 0.0)
    low_complexity_boundary = float(row.get("repeat_boundary_score", 0.0) or 0.0)
    compression = float(row.get("compression_anomaly_score", row.get("compression_boundary_score", 0.0)) or 0.0)
    entropy_asym = float(row.get("left_right_entropy_asymmetry", 0.0) or 0.0)
    void = float(row.get("local_feature_void_score", 0.0) or 0.0)
    decoherence = float(np.mean([entropy_cliff, low_complexity_boundary, compression, entropy_asym, void]))
    return {
        "local_entropy_variance_across_subwindows": local_entropy_variance,
        "entropy_cliff_score": entropy_cliff,
        "low_complexity_boundary_score": low_complexity_boundary,
        "entropy_noise_compression_anomaly_score": compression,
        "entropy_asymmetry_score": entropy_asym,
        "feature_void_score": void,
        "decoherence_boundary_candidate_score": decoherence,
    }
