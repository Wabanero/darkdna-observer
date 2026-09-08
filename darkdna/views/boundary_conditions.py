"""Boundary-condition primitive view."""

from __future__ import annotations

from darkdna.utils.stats import finite_mean, optional_row_float


def compute_boundary_condition_view(row: dict) -> dict[str, float]:
    sequence_boundary = optional_row_float(row, "left_right_regime_difference_score")
    entropy_transition = optional_row_float(row, "entropy_boundary_score")
    repeat_unique = optional_row_float(row, "repeat_boundary_score")
    te_boundary = optional_row_float(row, "TE_boundary_score")
    negative_boundary = optional_row_float(row, "negative_space_boundary_score")
    compression = optional_row_float(row, "compression_anomaly_score", "compression_boundary_score")
    entropy_asym = optional_row_float(row, "left_right_entropy_asymmetry")
    candidate = finite_mean(
        [
            sequence_boundary,
            entropy_transition,
            repeat_unique,
            te_boundary,
            negative_boundary,
            compression,
            entropy_asym,
        ]
    )
    return {
        "sequence_regime_boundary_score": sequence_boundary,
        "entropy_transition_score": entropy_transition,
        "repeat_to_unique_transition_score": repeat_unique,
        "TE_boundary_behavior_score": te_boundary,
        "negative_space_boundary_score_view": negative_boundary,
        "compression_transition_score": compression,
        "entropy_asymmetry_boundary_score": entropy_asym,
        "boundary_condition_candidate_score": candidate,
    }
