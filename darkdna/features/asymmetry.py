"""Directional and strand asymmetry features."""

from __future__ import annotations

import json

import numpy as np

from .sequence import clean_sequence, compression_ratio, kmer_counts


def skew(a: int, b: int) -> float:
    denom = a + b
    return 0.0 if denom == 0 else (a - b) / denom


def entropy(seq: str) -> float:
    from darkdna.utils.stats import shannon_entropy

    return shannon_entropy([b for b in seq if b in "ACGT"])


def compute_asymmetry_features(seq: str) -> dict[str, float | str]:
    seq = clean_sequence(seq)
    mid = len(seq) // 2
    left, right = seq[:mid], seq[mid:]
    counts = {base: seq.count(base) for base in "ACGT"}
    left_gc = (left.count("G") + left.count("C")) / max(1, sum(left.count(b) for b in "ACGT"))
    right_gc = (right.count("G") + right.count("C")) / max(1, sum(right.count(b) for b in "ACGT"))
    strand_asym = {}
    comp = str.maketrans("ACGT", "TGCA")
    counts4 = kmer_counts(seq, 4)
    for token, count in counts4.items():
        rc = token.translate(comp)[::-1]
        denom = count + counts4.get(rc, 0)
        if denom:
            strand_asym[token] = abs(count - counts4.get(rc, 0)) / denom
    return {
        "GC_skew": skew(counts["G"], counts["C"]),
        "AT_skew": skew(counts["A"], counts["T"]),
        "G_skew": skew(seq[:mid].count("G"), seq[mid:].count("G")),
        "C_skew": skew(seq[:mid].count("C"), seq[mid:].count("C")),
        "purine_skew": skew(counts["A"] + counts["G"], counts["C"] + counts["T"]),
        "pyrimidine_skew": skew(counts["C"], counts["T"]),
        "left_right_GC_asymmetry": abs(left_gc - right_gc),
        "left_right_entropy_asymmetry": abs(entropy(left) - entropy(right)),
        "left_right_compression_asymmetry": abs(compression_ratio(left or "N") - compression_ratio(right or "N")),
        "kmer_strand_asymmetry": float(np.mean(list(strand_asym.values()))) if strand_asym else 0.0,
        "directional_repeat_asymmetry": abs(left.count("G") - right.count("G")) / max(1, seq.count("G")),
        "orientation_bias_of_recurrent_kmers": float(np.mean(list(strand_asym.values()))) if strand_asym else 0.0,
        "kmer_strand_asymmetry_profile": json.dumps(dict(sorted(strand_asym.items())[:20])),
    }
