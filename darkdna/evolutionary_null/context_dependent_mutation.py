"""Reference-conditioned first-order sequence generator."""

from __future__ import annotations

import random
from collections import Counter, defaultdict


def fit_transition_counts(sequences: list[str]) -> dict[str, dict[str, int]]:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for sequence in sequences:
        clean = "".join(base for base in sequence.upper() if base in "ACGT")
        for left, right in zip(clean, clean[1:]):
            counts[left][right] += 1
    return {left: dict(rights) for left, rights in counts.items()}


def generate_markov_sequence(
    length: int,
    transitions: dict[str, dict[str, int]],
    base_frequencies: dict[str, float],
    *,
    seed: int,
) -> str:
    if length <= 0:
        return ""
    rng = random.Random(seed)
    bases = list("ACGT")
    base_weights = [max(0.0, float(base_frequencies.get(base, 0.25))) for base in bases]
    current = rng.choices(bases, weights=base_weights, k=1)[0]
    output = [current]
    for _ in range(length - 1):
        options = transitions.get(current, {})
        if options:
            next_bases = sorted(options)
            weights = [max(0, int(options[base])) for base in next_bases]
            current = rng.choices(next_bases, weights=weights, k=1)[0]
        else:
            current = rng.choices(bases, weights=base_weights, k=1)[0]
        output.append(current)
    return "".join(output)
