"""Sequence-derived non-B-DNA propensity proxies.

These scores are sequence-derived susceptibility proxies only. They do not
claim experimentally confirmed structures.
"""

from __future__ import annotations

import numpy as np
from collections import Counter

from .repeats import g_tract_density
from .sequence import clean_sequence


def longest_alternating(seq: str, groups: tuple[set[str], set[str]]) -> int:
    best = 0
    current = 0
    last_group = None
    for base in seq:
        group = 0 if base in groups[0] else 1 if base in groups[1] else None
        if group is None:
            current = 0
            last_group = None
        elif last_group is None or group != last_group:
            current += 1
            last_group = group
        else:
            current = 1
            last_group = group
        best = max(best, current)
    return best


def repeat_density(seq: str, mode: str = "direct", k: int = 6) -> float:
    seq = clean_sequence(seq)
    if len(seq) < 2 * k:
        return 0.0
    comp = str.maketrans("ACGT", "TGCA")
    counts = Counter(seq[i : i + k] for i in range(len(seq) - k + 1) if "N" not in seq[i : i + k])
    total = sum(counts.values())
    if total == 0:
        return 0.0
    hits = 0
    for token, count in counts.items():
        if mode == "direct":
            if count > 1:
                hits += count
        elif mode == "inverted":
            if counts.get(token.translate(comp)[::-1], 0) > 0:
                hits += count
        elif mode == "mirror":
            if counts.get(token[::-1], 0) > 0:
                hits += count
    return hits / total


def compute_nonb_dna_features(seq: str) -> dict[str, float]:
    seq = clean_sequence(seq)
    length = max(1, len(seq))
    g = seq.count("G")
    c = seq.count("C")
    gc = g + c
    g_skew = (g - c) / max(1, gc)
    g4 = g_tract_density(seq, min_run=3) * (1 + max(0.0, g_skew))
    z_alt = longest_alternating(seq, ({"G", "C"}, {"A", "T"})) / length
    purine_pyr_alt = longest_alternating(seq, ({"A", "G"}, {"C", "T"})) / length
    a_phased = sum(1 for i in range(0, max(0, len(seq) - 4), 10) if seq[i : i + 4].count("A") >= 3) / max(1, length / 10)
    inverted = repeat_density(seq, "inverted")
    direct = repeat_density(seq, "direct")
    mirror = repeat_density(seq, "mirror")
    triplex = (purine_pyr_alt + max(0.0, g_skew) + seq.count("AGG") / length) / 3
    cruciform = (inverted + mirror) / 2
    gc_skew = (g - c) / max(1, gc)
    c_skew = (c - g) / max(1, gc)
    rloop = np.mean([max(0.0, gc_skew), max(0.0, g_skew), max(0.0, -c_skew), g / length, abs(g_skew)])
    aggregate = float(np.mean([g4, z_alt, a_phased, inverted, direct, mirror, triplex, cruciform, rloop]))
    return {
        "G4_susceptibility_proxy": float(g4),
        "Z_DNA_propensity_proxy": float(z_alt),
        "A_phased_tract_score": float(a_phased),
        "inverted_repeat_density": float(inverted),
        "direct_repeat_density": float(direct),
        "mirror_repeat_density": float(mirror),
        "triplex_H_DNA_proxy": float(triplex),
        "cruciform_forming_potential": float(cruciform),
        "R_loop_forming_potential": float(rloop),
        "non_B_DNA_aggregate_score": aggregate,
    }
