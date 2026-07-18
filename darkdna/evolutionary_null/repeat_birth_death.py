"""Repeat birth/death proxy for neutral sequence generation."""

from __future__ import annotations

import random


def apply_repeat_turnover(seq: str, rate: float, *, seed: int) -> str:
    if len(seq) < 20 or rate <= 0:
        return seq
    rng = random.Random(seed)
    values = list(seq)
    events = int(len(values) * min(rate, 0.1) / 5)
    for _ in range(events):
        start = rng.randrange(0, len(values) - 5)
        motif = values[start : start + rng.choice((2, 3, 4))]
        target = rng.randrange(0, len(values) - len(motif) + 1)
        values[target : target + len(motif)] = motif
    return "".join(values)
