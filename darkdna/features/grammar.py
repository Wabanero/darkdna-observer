"""Sequence grammar features independent of motif databases."""

from __future__ import annotations

import json
import random
from collections import Counter, defaultdict

import numpy as np

from darkdna.utils.stats import shannon_entropy
from .sequence import clean_sequence, compression_ratio, kmer_counts, kmer_frequencies


def autocorrelation_periodicity(seq: str, mapping: dict[str, float] | None = None, max_lag: int = 200) -> dict[str, float]:
    seq = clean_sequence(seq)
    if len(seq) < 4:
        return {"best_lag": 0.0, "best_score": 0.0, "lag_10_score": 0.0, "lag_147_score": 0.0}
    mapping = mapping or {"G": 1.0, "C": 1.0, "A": -1.0, "T": -1.0, "N": 0.0}
    x = np.array([mapping.get(base, 0.0) for base in seq], dtype=float)
    x = x - x.mean()
    scores = {}
    for lag in range(1, min(max_lag, len(seq) - 1) + 1):
        denom = np.linalg.norm(x[:-lag]) * np.linalg.norm(x[lag:])
        scores[lag] = float(np.dot(x[:-lag], x[lag:]) / denom) if denom else 0.0
    best_lag, best_score = max(scores.items(), key=lambda item: abs(item[1])) if scores else (0, 0.0)
    return {
        "best_lag": float(best_lag),
        "best_score": float(abs(best_score)),
        "lag_10_score": float(abs(scores.get(10, 0.0))),
        "lag_147_score": float(abs(scores.get(147, 0.0))),
    }


def fourier_summary(seq: str) -> dict[str, float]:
    seq = clean_sequence(seq)
    if len(seq) < 4:
        return {"fourier_peak_frequency": 0.0, "fourier_peak_power": 0.0}
    x = np.array([1.0 if base in "GC" else -1.0 if base in "AT" else 0.0 for base in seq])
    spectrum = np.abs(np.fft.rfft(x - x.mean())) ** 2
    if spectrum.size <= 1:
        return {"fourier_peak_frequency": 0.0, "fourier_peak_power": 0.0}
    idx = int(np.argmax(spectrum[1:]) + 1)
    return {"fourier_peak_frequency": float(idx / len(seq)), "fourier_peak_power": float(spectrum[idx] / max(1.0, spectrum.sum()))}


def shuffled_forbidden_word_score(seq: str, k: int = 4, shuffles: int = 20, seed: int = 13) -> float:
    seq = clean_sequence(seq)
    observed = len(kmer_counts(seq, k))
    if len(seq) < k:
        return 0.0
    rng = random.Random(seed)
    expected = []
    bases = list(seq)
    for _ in range(shuffles):
        rng.shuffle(bases)
        expected.append(len(kmer_counts("".join(bases), k)))
    mean = float(np.mean(expected)) if expected else observed
    return max(0.0, (mean - observed) / max(1.0, mean))


def entropy_cliffs(seq: str, block: int = 50) -> float:
    seq = clean_sequence(seq)
    if len(seq) < block * 2:
        block = max(5, len(seq) // 4)
    entropies = [shannon_entropy(seq[i : i + block]) for i in range(0, len(seq), block) if seq[i : i + block]]
    if len(entropies) < 2:
        return 0.0
    return float(np.max(np.abs(np.diff(entropies))))


def kmer_adjacency(seq: str, k: int = 4) -> tuple[int, int, float]:
    seq = clean_sequence(seq)
    edges: dict[tuple[str, str], int] = defaultdict(int)
    nodes = set()
    for i in range(len(seq) - k):
        a = seq[i : i + k]
        b = seq[i + 1 : i + 1 + k]
        if "N" in a or "N" in b:
            continue
        nodes.add(a)
        nodes.add(b)
        edges[(a, b)] += 1
    entropy = shannon_entropy(edges)
    return len(nodes), len(edges), entropy


def markov_surprise(seq: str, k: int = 2) -> float:
    seq = clean_sequence(seq)
    if len(seq) < k + 2:
        return 0.0
    context_counts = Counter(seq[i : i + k] for i in range(len(seq) - k) if "N" not in seq[i : i + k])
    trans_counts = Counter(seq[i : i + k + 1] for i in range(len(seq) - k) if "N" not in seq[i : i + k + 1])
    surprises = []
    for token, count in trans_counts.items():
        context = token[:-1]
        prob = count / max(1, context_counts[context])
        surprises.append(-np.log2(max(prob, 1e-12)))
    return float(np.mean(surprises)) if surprises else 0.0


def distant_mutual_information(seq: str, lag: int | None = None) -> float:
    seq = clean_sequence(seq)
    if len(seq) < 20:
        return 0.0
    lag = lag or max(5, len(seq) // 4)
    pairs = [(seq[i], seq[i + lag]) for i in range(len(seq) - lag) if seq[i] in "ACGT" and seq[i + lag] in "ACGT"]
    if not pairs:
        return 0.0
    joint = Counter(pairs)
    left = Counter(a for a, _ in pairs)
    right = Counter(b for _, b in pairs)
    total = len(pairs)
    mi = 0.0
    for (a, b), count in joint.items():
        pxy = count / total
        px = left[a] / total
        py = right[b] / total
        mi += pxy * np.log2(pxy / max(px * py, 1e-12))
    return float(mi)


def compute_grammar_features(seq: str) -> dict[str, float | str]:
    seq = clean_sequence(seq)
    ac = autocorrelation_periodicity(seq)
    fft = fourier_summary(seq)
    nodes, edges, graph_entropy = kmer_adjacency(seq)
    freqs = kmer_frequencies(seq, 4)
    recurrent = [k for k, v in freqs.items() if v >= 0.02]
    nested = sum(seq.count(token * 2) for token in recurrent[:20]) / max(1, len(seq))
    compression_disc = abs(compression_ratio(seq[: len(seq) // 2] or "N") - compression_ratio(seq[len(seq) // 2 :] or "N"))
    return {
        "spacing_periodicity_autocorrelation": ac["best_score"],
        "spacing_periodicity_best_lag": ac["best_lag"],
        "spacing_periodicity_fourier_power": fft["fourier_peak_power"],
        "spacing_periodicity_fourier_frequency": fft["fourier_peak_frequency"],
        "phase_periodicity_around_10bp": ac["lag_10_score"],
        "nucleosome_scale_periodicity_around_147bp": ac["lag_147_score"],
        "forbidden_word_depletion_enrichment": shuffled_forbidden_word_score(seq),
        "orientation_asymmetry": abs(seq.count("A") - seq.count("T")) / max(1, seq.count("A") + seq.count("T")),
        "palindrome_break_score": 1.0 - min(1.0, abs(seq.count("AT") - seq.count("TA")) / max(1, len(seq))),
        "local_entropy_cliffs": entropy_cliffs(seq),
        "compression_discontinuity": compression_disc,
        "motif_like_token_recurrence": len(recurrent) / max(1, len(freqs)),
        "kmer_adjacency_graph_nodes": float(nodes),
        "kmer_adjacency_graph_edges": float(edges),
        "motif_token_cooccurrence_graph_edges": float(edges),
        "recursive_block_detection_proxy": nested,
        "nested_repeat_architecture_score": nested,
        "grammar_entropy": graph_entropy,
        "Markov_order_anomaly": markov_surprise(seq, k=2),
        "n_gram_transition_surprise": markov_surprise(seq, k=3),
        "long_range_dependency_proxy": distant_mutual_information(seq),
        "mutual_information_between_distant_sequence_positions": distant_mutual_information(seq),
        "frequent_kmer_tokens": json.dumps(recurrent[:25]),
    }
