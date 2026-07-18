"""Length-preserving proxy for regional indel turnover."""

from __future__ import annotations

import random


def apply_length_preserving_indel_turnover(seq: str, rate: float, *, seed: int) -> str:
    """Move short sequence blocks while preserving benchmark length.

    True indels require ancestral alignments. This transformation is explicitly a
    turnover proxy for fixed-coordinate benchmarks.
    """

    if len(seq) < 8 or rate <= 0:
        return seq
    rng = random.Random(seed)
    values = list(seq)
    moves = int(len(values) * min(rate, 0.25))
    for _ in range(moves):
        source = rng.randrange(len(values))
        base = values.pop(source)
        values.insert(rng.randrange(len(values) + 1), base)
    return "".join(values)
