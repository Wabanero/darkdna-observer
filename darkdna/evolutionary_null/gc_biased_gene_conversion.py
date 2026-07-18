"""GC-biased gene-conversion perturbation for neutral surrogates."""

from __future__ import annotations

import random


def apply_gc_bias(seq: str, strength: float, *, seed: int) -> str:
    rng = random.Random(seed)
    output = []
    probability = min(1.0, max(0.0, abs(strength)))
    for base in seq.upper():
        if strength > 0 and base in "AT" and rng.random() < probability:
            output.append("G" if base == "A" else "C")
        elif strength < 0 and base in "GC" and rng.random() < probability:
            output.append("A" if base == "G" else "T")
        else:
            output.append(base)
    return "".join(output)
