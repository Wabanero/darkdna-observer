"""CpG-deamination component for neutral surrogates."""

from __future__ import annotations

import random


def apply_cpg_deamination(seq: str, probability: float, *, seed: int) -> str:
    rng = random.Random(seed)
    values = list(seq.upper())
    for index in range(len(values) - 1):
        if values[index] == "C" and values[index + 1] == "G" and rng.random() < max(0.0, probability):
            values[index] = "T"
    return "".join(values)
