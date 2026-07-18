"""Regional neutral-rate modifiers."""

from __future__ import annotations


def regional_rate_modifier(*values: float | None) -> tuple[float, str]:
    finite = [float(value) for value in values if value is not None]
    if not finite:
        return 1.0, "unavailable_generic_rate"
    return max(0.1, min(10.0, sum(finite) / len(finite))), "available_relative_rate"
