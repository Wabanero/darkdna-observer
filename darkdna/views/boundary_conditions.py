"""Boundary-condition primitive view."""

from __future__ import annotations

import numpy as np


def compute_boundary_condition_view(row: dict) -> dict[str, float]:
    sequence_boundary = float(row.get("left_right_regime_difference_score", 0.0) or 0.0)
    entropy_transition = float(row.get("entropy_boundary_score", 0.0) or 0.0)
    repeat_unique = float(row.get("repeat_boundary_score", 0.0) or 0.0)
    te_boundary = float(row.get("TE_boundary_score", 0.0) or 0.0)
    negative_boundary = float(row.get("negative_space_boundary_score", 0.0) or 0.0)
    candidate = float(np.mean([sequence_boundary, entropy_transition, repeat_unique, te_boundary, negative_boundary]))
    return {
        "sequence_regime_boundary_score": sequence_boundary,
        "entropy_transition_score": entropy_transition,
        "repeat_to_unique_transition_score": repeat_unique,
        "TE_boundary_behavior_score": te_boundary,
        "negative_space_boundary_score_view": negative_boundary,
        "boundary_condition_candidate_score": candidate,
    }
