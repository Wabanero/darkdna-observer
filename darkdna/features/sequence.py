"""Core sequence feature extraction."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from darkdna.io.fasta import read_fasta
from darkdna.utils.progress import ProgressReporter, progress_message
from darkdna.utils.stats import safe_divide, shannon_entropy
from .compression import calibrate_compression, compressed_size
from .repeats import g_tract_density, homopolymer_runs, palindrome_density, periodicity_proxy, simple_repeat_fraction


DNA = "ACGT"


NUMERIC_WALK_MAPPINGS: dict[str, dict[str, object]] = {
    "purine_pyrimidine_numeric_walk": {
        "rationale": "Contrasts purines (A/G) with pyrimidines (C/T).",
        "mapping": {"A": 1.0, "G": 1.0, "C": -1.0, "T": -1.0},
        "orientation_sensitivity": "orientation invariant; sign reverses under complement",
        "reverse_complement_behaviour": "reversal and sign inversion of step order",
        "null_expectation": "zero expected step for balanced purine/pyrimidine composition",
    },
    "strong_weak_H_bond_numeric_walk": {
        "rationale": "Contrasts three-hydrogen-bond G/C with two-hydrogen-bond A/T bases.",
        "mapping": {"G": 1.0, "C": 1.0, "A": -1.0, "T": -1.0},
        "orientation_sensitivity": "orientation invariant",
        "reverse_complement_behaviour": "same step values in reversed order",
        "null_expectation": "zero expected step at 50 percent GC",
    },
    "amino_keto_numeric_walk": {
        "rationale": "Contrasts amino bases (A/C) with keto bases (G/T).",
        "mapping": {"A": 1.0, "C": 1.0, "G": -1.0, "T": -1.0},
        "orientation_sensitivity": "orientation sensitive",
        "reverse_complement_behaviour": "same signs in reversed order after reverse complement",
        "null_expectation": "zero expected step for balanced amino/keto composition",
    },
}

NUMERIC_WALK_ALIASES = {
    "GC_AT_numeric_walk": "strong_weak_H_bond_numeric_walk",
}


def numeric_walk_mapping_registry() -> list[dict[str, object]]:
    """Return canonical, non-duplicated walk mappings and their assumptions."""

    return [
        {
            "name": name,
            "biological_or_mathematical_rationale": definition["rationale"],
            "mapping_table": definition["mapping"],
            "orientation_sensitivity": definition["orientation_sensitivity"],
            "reverse_complement_behaviour": definition["reverse_complement_behaviour"],
            "null_expectation": definition["null_expectation"],
        }
        for name, definition in NUMERIC_WALK_MAPPINGS.items()
    ]


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
    """Deprecated raw compression ratio retained as a migration alias.

    Use the header-corrected and same-length-null-calibrated fields returned by
    :func:`compute_sequence_features` for inference.
    """

    cleaned = clean_sequence(seq)
    if not cleaned:
        return 0.0
    return compressed_size(cleaned, method) / len(cleaned)


def normalized_compression_distance(seq: str, controls: Iterable[str] | None = None) -> float:
    controls = list(controls or [])
    if not controls:
        return np.nan
    x = clean_sequence(seq)
    cx = compressed_size(x, "gzip")
    distances = []
    for control in controls:
        y = clean_sequence(control)
        cy = compressed_size(y, "gzip")
        cxy = compressed_size(x + y, "gzip")
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
    compression = calibrate_compression(seq, n_surrogates=2, seed=13)
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
        "raw_compression_ratio_deprecation_warning": (
            "Raw compressor ratios are retained only as compatibility aliases; do not compare them across window sizes. "
            "Use header-corrected, same-length-null-calibrated compression fields."
        ),
        "simple_repeat_fraction": simple_repeat_fraction(seq),
        "tandem_repeat_like_periodicity_proxy": periodicity_proxy(seq),
        "palindrome_density": palindrome_density(seq),
        "G_tract_density": g_tract_density(seq),
        "predicted_G_quadruplex_proxy_score": g_tract_density(seq, min_run=3) * (1.0 + safe_divide(counts["G"], usable, 0.0)),
        **compression,
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
    features["numeric_walk_mapping_registry"] = json.dumps(numeric_walk_mapping_registry(), sort_keys=True)
    for prefix, definition in NUMERIC_WALK_MAPPINGS.items():
        mapping = definition["mapping"]
        summary = walk_summary(seq, mapping)
        for key, value in summary.items():
            features[f"{prefix}_{key}"] = value
    for alias, canonical in NUMERIC_WALK_ALIASES.items():
        for key in ("final", "range", "std"):
            features[f"{alias}_{key}"] = features[f"{canonical}_{key}"]
        features[f"{alias}_deprecation_warning"] = (
            f"{alias} duplicated {canonical} and is a deprecated compatibility alias. "
            f"Use {canonical}."
        )
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


def extract_features_for_windows(windows: pd.DataFrame, fasta: str | Path, *, progress: bool = False) -> pd.DataFrame:
    if progress:
        progress_message("extract-features", "loading FASTA sequence")
    genome = read_fasta(fasta)
    rows = []
    reporter = ProgressReporter("extract-features", total=len(windows)) if progress else None
    if reporter:
        reporter.start("computing sequence features")
    for idx, row in enumerate(windows.itertuples(), start=1):
        seq = genome.get(str(row.chrom), "")[int(row.start) : int(row.end)]
        features = compute_all_sequence_features(seq)
        rows.append({"region_id": row.region_id, "chrom": row.chrom, "start": row.start, "end": row.end, **features})
        if reporter:
            reporter.update(idx, message=str(row.region_id))
    if reporter:
        reporter.finish()
    return pd.DataFrame(rows)


def write_sequence_features(features: pd.DataFrame, outdir: str | Path) -> Path:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "sequence_features.parquet"
    features.to_parquet(path, index=False)
    return path
