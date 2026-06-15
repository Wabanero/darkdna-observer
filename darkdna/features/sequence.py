"""Core sequence feature extraction."""

from __future__ import annotations

import bz2
import gzip
import json
import lzma
from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from darkdna.io.fasta import read_fasta
from darkdna.utils.stats import safe_divide, shannon_entropy
from .repeats import g_tract_density, homopolymer_runs, palindrome_density, periodicity_proxy, simple_repeat_fraction


DNA = "ACGT"


def clean_sequence(seq: str) -> str:
    return "".join(base if base in "ACGTN" else "N" for base in seq.upper())


def base_counts(seq: str) -> Counter:
    return Counter(clean_sequence(seq))


def gc_content(seq: str) -> float:
    counts = base_counts(seq)
    usable = counts["A"] + counts["C"] + counts["G"] + counts["T"]
    return safe_divide(counts["G"] + counts["C"], usable, default=np.nan)


def at_content(seq: str) -> float:
    counts = base_counts(seq)
    usable = counts["A"] + counts["C"] + counts["G"] + counts["T"]
    return safe_divide(counts["A"] + counts["T"], usable, default=np.nan)


def n_fraction(seq: str) -> float:
    seq = clean_sequence(seq)
    return safe_divide(seq.count("N"), len(seq), default=1.0)


def cpg_density(seq: str) -> float:
    seq = clean_sequence(seq)
    return safe_divide(sum(1 for i in range(len(seq) - 1) if seq[i : i + 2] == "CG"), max(1, len(seq) - 1))


def cpg_observed_expected(seq: str) -> float:
    seq = clean_sequence(seq)
    usable = len([b for b in seq if b in DNA])
    if usable == 0:
        return np.nan
    obs = sum(1 for i in range(len(seq) - 1) if seq[i : i + 2] == "CG")
    c = seq.count("C")
    g = seq.count("G")
    expected = (c * g) / max(1, usable)
    return safe_divide(obs, expected, default=0.0)


def kmer_counts(seq: str, k: int) -> Counter:
    seq = clean_sequence(seq)
    return Counter(seq[i : i + k] for i in range(len(seq) - k + 1) if "N" not in seq[i : i + k])


def kmer_frequencies(seq: str, k: int) -> dict[str, float]:
    counts = kmer_counts(seq, k)
    total = sum(counts.values())
    if total == 0:
        return {}
    return {kmer: count / total for kmer, count in sorted(counts.items())}


def kmer_entropy(seq: str, k: int) -> float:
    return shannon_entropy(kmer_counts(seq, k))


def block_entropy_profile(seq: str, block_size: int = 50) -> list[float]:
    seq = clean_sequence(seq)
    if not seq:
        return []
    if len(seq) < block_size:
        block_size = max(1, len(seq) // 4 or len(seq))
    return [shannon_entropy(seq[idx : idx + block_size]) for idx in range(0, len(seq), block_size)]


def lempel_ziv_complexity(seq: str) -> float:
    seq = clean_sequence(seq)
    if not seq:
        return 0.0
    dictionary: set[str] = set()
    token = ""
    complexity = 0
    for char in seq:
        token += char
        if token not in dictionary:
            dictionary.add(token)
            complexity += 1
            token = ""
    return safe_divide(complexity, len(seq), default=0.0)


def compression_ratio(seq: str, method: str = "gzip") -> float:
    raw = clean_sequence(seq).encode("ascii")
    if not raw:
        return 0.0
    if method == "gzip":
        compressed = gzip.compress(raw)
    elif method == "bz2":
        compressed = bz2.compress(raw)
    elif method == "lzma":
        compressed = lzma.compress(raw)
    else:
        raise ValueError(f"Unknown compression method: {method}")
    return len(compressed) / len(raw)


def normalized_compression_distance(seq: str, controls: Iterable[str] | None = None) -> float:
    controls = list(controls or [])
    if not controls:
        return np.nan
    x = clean_sequence(seq).encode("ascii")
    cx = len(gzip.compress(x))
    distances = []
    for control in controls:
        y = clean_sequence(control).encode("ascii")
        cy = len(gzip.compress(y))
        cxy = len(gzip.compress(x + y))
        distances.append((cxy - min(cx, cy)) / max(cx, cy, 1))
    return float(np.mean(distances)) if distances else np.nan


def numeric_walk(seq: str, mapping: dict[str, float]) -> np.ndarray:
    seq = clean_sequence(seq)
    steps = np.array([mapping.get(base, 0.0) for base in seq], dtype=float)
    return np.cumsum(steps)


def walk_summary(seq: str, mapping: dict[str, float]) -> dict[str, float]:
    walk = numeric_walk(seq, mapping)
    if walk.size == 0:
        return {"final": 0.0, "range": 0.0, "std": 0.0}
    return {"final": float(walk[-1]), "range": float(walk.max() - walk.min()), "std": float(walk.std())}


def compute_sequence_features(seq: str, controls: Iterable[str] | None = None) -> dict[str, object]:
    seq = clean_sequence(seq)
    counts = base_counts(seq)
    usable = counts["A"] + counts["C"] + counts["G"] + counts["T"]
    features: dict[str, object] = {
        "length": len(seq),
        "gc_content": gc_content(seq),
        "at_content": at_content(seq),
        "n_fraction": n_fraction(seq),
        "CpG_density": cpg_density(seq),
        "observed_expected_CpG": cpg_observed_expected(seq),
        "Shannon_entropy": shannon_entropy([base for base in seq if base in "ACGT"]),
        "block_entropy_profile": json.dumps(block_entropy_profile(seq)),
        "Lempel_Ziv_complexity": lempel_ziv_complexity(seq),
        "gzip_compression_ratio": compression_ratio(seq, "gzip"),
        "bz2_compression_ratio": compression_ratio(seq, "bz2"),
        "lzma_compression_ratio": compression_ratio(seq, "lzma"),
        "normalized_compression_distance_to_matched_controls": normalized_compression_distance(seq, controls),
        "simple_repeat_fraction": simple_repeat_fraction(seq),
        "tandem_repeat_like_periodicity_proxy": periodicity_proxy(seq),
        "palindrome_density": palindrome_density(seq),
        "G_tract_density": g_tract_density(seq),
        "predicted_G_quadruplex_proxy_score": g_tract_density(seq, min_run=3) * (1.0 + safe_divide(counts["G"], usable, 0.0)),
    }
    for base in "ACGTN":
        features[f"mono_freq_{base}"] = safe_divide(counts[base], len(seq), default=0.0)
    for k in [2, 3, 4, 5, 6]:
        freqs = kmer_frequencies(seq, k)
        features[f"k{k}_frequencies"] = json.dumps(freqs, sort_keys=True)
        features[f"k{k}_mer_entropy"] = kmer_entropy(seq, k)
    features["di_nucleotide_frequencies"] = features["k2_frequencies"]
    features["tri_nucleotide_frequencies"] = features["k3_frequencies"]
    runs = homopolymer_runs(seq)
    features["homopolymer_run_max"] = max(runs) if runs else 0
    features["homopolymer_run_mean"] = float(np.mean(runs)) if runs else 0.0
    features["homopolymer_run_p95"] = float(np.percentile(runs, 95)) if runs else 0.0
    for prefix, mapping in {
        "purine_pyrimidine_numeric_walk": {"A": 1, "G": 1, "C": -1, "T": -1},
        "strong_weak_H_bond_numeric_walk": {"G": 1, "C": 1, "A": -1, "T": -1},
        "GC_AT_numeric_walk": {"G": 1, "C": 1, "A": -1, "T": -1},
    }.items():
        summary = walk_summary(seq, mapping)
        for key, value in summary.items():
            features[f"{prefix}_{key}"] = value
    return features


def compute_all_sequence_features(seq: str, controls: Iterable[str] | None = None) -> dict[str, object]:
    features = compute_sequence_features(seq, controls=controls)
    from .asymmetry import compute_asymmetry_features
    from .boundaries import compute_boundary_features
    from .grammar import compute_grammar_features
    from .negative_space import compute_negative_space_features
    from .nonb_dna import compute_nonb_dna_features
    from .physical_shape import compute_physical_shape_features
    from darkdna.views.scale_fractal import compute_scale_fractal_features

    for fn in [
        compute_nonb_dna_features,
        compute_physical_shape_features,
        compute_asymmetry_features,
        compute_grammar_features,
        compute_boundary_features,
        compute_negative_space_features,
        compute_scale_fractal_features,
    ]:
        features.update(fn(seq))
    return features


def extract_features_for_windows(windows: pd.DataFrame, fasta: str | Path) -> pd.DataFrame:
    genome = read_fasta(fasta)
    rows = []
    for row in windows.itertuples():
        seq = genome.get(str(row.chrom), "")[int(row.start) : int(row.end)]
        features = compute_all_sequence_features(seq)
        rows.append({"region_id": row.region_id, "chrom": row.chrom, "start": row.start, "end": row.end, **features})
    return pd.DataFrame(rows)


def write_sequence_features(features: pd.DataFrame, outdir: str | Path) -> Path:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "sequence_features.parquet"
    features.to_parquet(path, index=False)
    return path
