"""Sequence-regime boundary features."""

from __future__ import annotations

import numpy as np

from darkdna.utils.optional_deps import optional_import
from darkdna.utils.stats import shannon_entropy
from .repeats import simple_repeat_fraction
from .sequence import clean_sequence, compression_ratio, cpg_density, gc_content, kmer_frequencies


def left_right(seq: str) -> tuple[str, str]:
    seq = clean_sequence(seq)
    mid = len(seq) // 2
    return seq[:mid], seq[mid:]


def distribution_shift(left: dict[str, float], right: dict[str, float]) -> float:
    keys = set(left) | set(right)
    return float(sum(abs(left.get(k, 0.0) - right.get(k, 0.0)) for k in keys) / 2)


def sliding_transition(seq: str, block: int = 50) -> float:
    seq = clean_sequence(seq)
    if len(seq) < 20:
        return 0.0
    block = min(block, max(5, len(seq) // 4))
    values = []
    for idx in range(0, max(1, len(seq) - block + 1), block):
        chunk = seq[idx : idx + block]
        values.append([gc_content(chunk), shannon_entropy(chunk), simple_repeat_fraction(chunk)])
    arr = np.array(values, dtype=float)
    if len(arr) < 2:
        return 0.0
    return float(np.nanmax(np.linalg.norm(np.diff(arr, axis=0), axis=1)))


def rupture_breakpoint_score(seq: str) -> float:
    ruptures, warning = optional_import("ruptures")
    if ruptures is None:
        return sliding_transition(seq)
    try:  # pragma: no cover - optional dependency not always installed.
        x = np.array([[1.0 if base in "GC" else -1.0 if base in "AT" else 0.0] for base in clean_sequence(seq)])
        if len(x) < 20:
            return 0.0
        algo = ruptures.Pelt(model="rbf").fit(x)
        breakpoints = algo.predict(pen=3)
        if not breakpoints:
            return 0.0
        center = len(seq) / 2
        return float(max(0.0, 1.0 - min(abs(bp - center) for bp in breakpoints) / max(1.0, center)))
    except Exception:
        return sliding_transition(seq)


def compute_boundary_features(seq: str) -> dict[str, float]:
    left, right = left_right(seq)
    entropy_boundary = abs(shannon_entropy(left) - shannon_entropy(right))
    gc_boundary = abs(gc_content(left) - gc_content(right))
    cpg_boundary = abs(cpg_density(left) - cpg_density(right))
    repeat_boundary = abs(simple_repeat_fraction(left) - simple_repeat_fraction(right))
    kshift = distribution_shift(kmer_frequencies(left, 4), kmer_frequencies(right, 4))
    comp_boundary = abs(compression_ratio(left or "N") - compression_ratio(right or "N"))
    walks = np.cumsum([1 if b in "GC" else -1 if b in "AT" else 0 for b in clean_sequence(seq)])
    if len(walks) >= 4:
        direction_change = abs(np.mean(np.diff(walks[: len(walks) // 2])) - np.mean(np.diff(walks[len(walks) // 2 :])))
    else:
        direction_change = 0.0
    segmentation = rupture_breakpoint_score(seq)
    regime = float(np.mean([entropy_boundary, gc_boundary, cpg_boundary, repeat_boundary, kshift, comp_boundary]))
    local_transition = sliding_transition(seq)
    return {
        "entropy_boundary_score": float(entropy_boundary),
        "GC_boundary_score": float(gc_boundary),
        "CpG_boundary_score": float(cpg_boundary),
        "repeat_boundary_score": float(repeat_boundary),
        "kmer_distribution_shift_score": float(kshift),
        "compression_boundary_score": float(comp_boundary),
        "numeric_walk_direction_change_score": float(direction_change),
        "segmentation_breakpoint_score": float(segmentation),
        "left_right_regime_difference_score": regime,
        "local_feature_transition_score": float(local_transition),
    }
