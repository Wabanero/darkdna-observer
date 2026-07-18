"""Deterministic sequence transformations used by null panels and benchmarks."""

from __future__ import annotations

import random
from collections import Counter, defaultdict


COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(COMPLEMENT)[::-1]


def mononucleotide_shuffle(seq: str, seed: int = 13) -> str:
    values = list(seq.upper())
    random.Random(seed).shuffle(values)
    return "".join(values)


def dinucleotide_preserving_shuffle(seq: str, seed: int = 13) -> str:
    seq = seq.upper()
    if len(seq) < 3:
        return mononucleotide_shuffle(seq, seed)
    rng = random.Random(seed)
    edges: dict[str, list[str]] = defaultdict(list)
    for left, right in zip(seq, seq[1:]):
        edges[left].append(right)
    for values in edges.values():
        rng.shuffle(values)
    current = seq[0]
    output = [current]
    bases = list(seq) or list("ACGT")
    for _ in range(len(seq) - 1):
        current = edges[current].pop() if edges.get(current) else rng.choice(bases)
        output.append(current)
    return "".join(output)


def kmer_preserving_shuffle(seq: str, k: int = 3, seed: int = 13) -> str:
    seq = seq.upper()
    if len(seq) <= k:
        return mononucleotide_shuffle(seq, seed)
    blocks = [seq[start : start + k] for start in range(0, len(seq), max(1, k))]
    random.Random(seed).shuffle(blocks)
    return "".join(blocks)[: len(seq)]


def local_mononucleotide_shuffle(seq: str, seed: int = 13, block_size: int = 1_000) -> str:
    rng = random.Random(seed)
    blocks = []
    for start in range(0, len(seq), max(1, block_size)):
        block = list(seq[start : start + block_size].upper())
        rng.shuffle(block)
        blocks.append("".join(block))
    return "".join(blocks)


def markov_chain_surrogate(seq: str, seed: int = 13) -> str:
    seq = seq.upper()
    if len(seq) < 2:
        return seq
    rng = random.Random(seed)
    transitions: dict[str, list[str]] = defaultdict(list)
    for left, right in zip(seq, seq[1:]):
        transitions[left].append(right)
    bases = [base for base in seq if base in "ACGTN"] or list("ACGT")
    current = rng.choice(bases)
    output = [current]
    for _ in range(len(seq) - 1):
        choices = transitions.get(current) or bases
        current = rng.choice(choices)
        output.append(current)
    return "".join(output)


def synthetic_equal_composition(seq: str, seed: int = 13) -> str:
    counts = Counter(seq.upper())
    values = [base for base in sorted(counts) for _ in range(counts[base])]
    random.Random(seed + 7_919).shuffle(values)
    return "".join(values)


def transform_sequence(
    seq: str,
    method: str,
    *,
    seed: int = 13,
    k: int = 3,
    block_size: int = 1_000,
) -> str:
    seq = seq.upper()
    methods = {
        "native": lambda: seq,
        "reversed_sequence": lambda: seq[::-1],
        "whole_genome_reverse": lambda: seq[::-1],
        "reverse_complement": lambda: reverse_complement(seq),
        "mononucleotide_preserving": lambda: mononucleotide_shuffle(seq, seed),
        "global_mononucleotide_shuffle": lambda: mononucleotide_shuffle(seq, seed),
        "local_mononucleotide_shuffle": lambda: local_mononucleotide_shuffle(seq, seed, block_size),
        "dinucleotide_preserving": lambda: dinucleotide_preserving_shuffle(seq, seed),
        "dinucleotide_preserving_shuffle": lambda: dinucleotide_preserving_shuffle(seq, seed),
        "global_dinucleotide_shuffle": lambda: dinucleotide_preserving_shuffle(seq, seed),
        "local_dinucleotide_shuffle": lambda: "".join(
            dinucleotide_preserving_shuffle(seq[start : start + block_size], seed + start)
            for start in range(0, len(seq), max(1, block_size))
        ),
        "kmer_preserving": lambda: kmer_preserving_shuffle(seq, k, seed),
        "kmer_preserving_shuffle": lambda: kmer_preserving_shuffle(seq, k, seed),
        "markov_chain_surrogate": lambda: markov_chain_surrogate(seq, seed),
        "synthetic_equal_composition": lambda: synthetic_equal_composition(seq, seed),
    }
    if method not in methods:
        raise ValueError(f"Unknown sequence null transformation: {method}")
    transformed = methods[method]()
    if len(transformed) != len(seq):
        raise RuntimeError(f"Null transformation {method} changed sequence length unexpectedly")
    return transformed
