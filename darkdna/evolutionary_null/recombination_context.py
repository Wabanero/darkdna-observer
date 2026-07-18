"""Recombination-context rate modifiers."""

from __future__ import annotations


def recombination_rate_modifier(value: float | None, median: float | None = None) -> tuple[float, str]:
    if value is None or median is None or median <= 0:
        return 1.0, "unavailable_no_recombination_track"
    return max(0.25, min(4.0, float(value) / float(median))), "available_relative_modifier"
