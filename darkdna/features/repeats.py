"""Repeat and low-complexity sequence features."""

from __future__ import annotations

import re
from collections import Counter

import numpy as np


def homopolymer_runs(seq: str) -> list[int]:
    seq = seq.upper()
    runs: list[int] = []
    if not seq:
        return runs
    current = seq[0]
    length = 1
    for base in seq[1:]:
        if base == current:
            length += 1
        else:
            runs.append(length)
            current = base
            length = 1
    runs.append(length)
    return runs


def simple_repeat_fraction(seq: str, min_period: int = 1, max_period: int = 6, min_copies: int = 3) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    covered = np.zeros(len(seq), dtype=bool)
    for period in range(min_period, max_period + 1):
        pattern = re.compile(f"([ACGT]{{{period}}})(?:\\1){{{min_copies - 1},}}")
        for match in pattern.finditer(seq):
            covered[match.start() : match.end()] = True
    return float(covered.mean())


def periodicity_proxy(seq: str, max_period: int = 20) -> float:
    seq = seq.upper()
    if len(seq) < 4:
        return 0.0
    scores = []
    for period in range(1, min(max_period, len(seq) // 2) + 1):
        matches = sum(1 for i in range(len(seq) - period) if seq[i] == seq[i + period] and seq[i] in "ACGT")
        scores.append(matches / max(1, len(seq) - period))
    return float(max(scores) if scores else 0.0)


def palindrome_density(seq: str, min_len: int = 4, max_len: int = 12) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    comp = str.maketrans("ACGTN", "TGCAN")
    hits = 0
    possible = 0
    for width in range(min_len, min(max_len, len(seq)) + 1):
        for idx in range(0, len(seq) - width + 1):
            token = seq[idx : idx + width]
            possible += 1
            if token == token.translate(comp)[::-1]:
                hits += 1
    return hits / max(1, possible)


def g_tract_density(seq: str, min_run: int = 3) -> float:
    runs = homopolymer_runs(seq.upper().replace("A", " ").replace("C", " ").replace("T", " ").replace("N", " "))
    # The replacement trick leaves spaces in the run list, so use a direct scan.
    seq = seq.upper()
    total = 0
    run = 0
    for base in seq:
        if base == "G":
            run += 1
        else:
            if run >= min_run:
                total += run
            run = 0
    if run >= min_run:
        total += run
    return total / max(1, len(seq))


def recurrent_kmer_orientation_bias(seq: str, k: int = 4) -> float:
    seq = seq.upper()
    if len(seq) < k:
        return 0.0
    comp = str.maketrans("ACGT", "TGCA")
    counts = Counter(seq[i : i + k] for i in range(len(seq) - k + 1) if "N" not in seq[i : i + k])
    diffs = []
    for kmer, count in counts.items():
        rc = kmer.translate(comp)[::-1]
        if kmer <= rc:
            denom = count + counts.get(rc, 0)
            if denom:
                diffs.append(abs(count - counts.get(rc, 0)) / denom)
    return float(np.mean(diffs)) if diffs else 0.0
