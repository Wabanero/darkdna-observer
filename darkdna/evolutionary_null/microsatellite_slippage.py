"""Microsatellite slippage proxy preserving total interval length."""

from __future__ import annotations

import random


def apply_slippage_proxy(seq: str, rate: float, *, seed: int) -> str:
    if len(seq) < 12 or rate <= 0:
        return seq
    rng = random.Random(seed)
    values = list(seq)
    for start in range(0, len(values) - 6, 6):
        if rng.random() < min(rate, 1.0):
            period = rng.choice((1, 2, 3))
            motif = values[start : start + period]
            values[start : start + period] = motif[::-1]
    return "".join(values)
