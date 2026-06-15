"""Negative-space features: structured absence and depleted tokens."""

from __future__ import annotations

import itertools

import numpy as np

from .repeats import g_tract_density, simple_repeat_fraction
from .sequence import clean_sequence, cpg_density, kmer_counts, kmer_frequencies


def depleted_kmer_score(seq: str, k: int = 4) -> float:
    seq = clean_sequence(seq)
    possible = ["".join(chars) for chars in itertools.product("ACGT", repeat=k)]
    seen = set(kmer_counts(seq, k))
    expected_possible = min(len(possible), max(1, len(seq) - k + 1))
    missing = len(set(possible[:]) - seen)
    return missing / max(1, len(possible)) if len(seq) >= k else 1.0


def local_void(seq: str, block: int = 50) -> float:
    seq = clean_sequence(seq)
    if len(seq) < block:
        block = max(5, len(seq) // 4 or len(seq))
    scores = []
    for idx in range(0, len(seq), block):
        chunk = seq[idx : idx + block]
        if chunk:
            scores.append(depleted_kmer_score(chunk, k=3))
    return float(np.max(scores)) if scores else 0.0


def compute_negative_space_features(seq: str) -> dict[str, float]:
    seq = clean_sequence(seq)
    freqs = kmer_frequencies(seq, 4)
    low_interaction_desert = 1.0 - min(1.0, sum(v for k, v in freqs.items() if k.count("G") >= 2 or k == k[::-1]))
    motif_desert = max(0.0, 1.0 - len(freqs) / max(1, min(4**4, len(seq))))
    repeat_desert = 1.0 - simple_repeat_fraction(seq)
    cpg_desert = 1.0 - min(1.0, cpg_density(seq) * 20)
    gtract_desert = 1.0 - min(1.0, g_tract_density(seq) * 20)
    depleted = depleted_kmer_score(seq, 4)
    void = local_void(seq)
    unexpected = float(np.mean([depleted, motif_desert, repeat_desert, cpg_desert, gtract_desert, low_interaction_desert]))
    return {
        "depleted_kmer_score": depleted,
        "forbidden_word_score": depleted,
        "motif_desert_score": motif_desert,
        "repeat_desert_score": repeat_desert,
        "CpG_desert_score": cpg_desert,
        "G_tract_desert_score": gtract_desert,
        "low_interaction_token_desert_score": low_interaction_desert,
        "unexpected_silence_score": unexpected,
        "local_feature_void_score": void,
        "negative_space_boundary_score": abs(void - unexpected),
    }
