"""Length-aware compression screens calibrated against explicit sequence nulls.

Compression is used here as a statistical-structure screen.  It is not
"functional information", and compressor byte counts are not comparable across
window sizes unless the header overhead and a same-length null are reported.
"""

from __future__ import annotations

import bz2
import gzip
import lzma
import math
import random
from collections import defaultdict
from functools import lru_cache

import numpy as np


COMPRESSOR_VERSION = "stdlib-python:gzip-level1,bz2-level1,lzma-preset0"
COMPRESSORS = ("gzip", "bz2", "lzma")


def _compress(payload: bytes, method: str) -> bytes:
    if method == "gzip":
        return gzip.compress(payload, compresslevel=1, mtime=0)
    if method == "bz2":
        return bz2.compress(payload, compresslevel=1)
    if method == "lzma":
        return lzma.compress(payload, preset=0)
    raise ValueError(f"Unknown compression method: {method}")


def compressed_size(seq: str, method: str = "gzip") -> int:
    """Return the raw compressed byte count for a sequence."""

    return len(_compress(seq.upper().encode("ascii", errors="ignore"), method))


@lru_cache(maxsize=None)
def compressor_header_bytes(method: str = "gzip") -> int:
    """Estimate fixed stream overhead from the compressor's empty payload."""

    return compressed_size("", method)


def header_corrected_compressed_size(seq: str, method: str = "gzip") -> int:
    """Return payload bytes after subtracting the empty-stream overhead."""

    return max(0, compressed_size(seq, method) - compressor_header_bytes(method))


def header_corrected_compression_ratio(seq: str, method: str = "gzip") -> float:
    if not seq:
        return math.nan
    return header_corrected_compressed_size(seq, method) / len(seq)


def _mononucleotide_shuffle(seq: str, rng: random.Random) -> str:
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def _dinucleotide_shuffle(seq: str, rng: random.Random) -> str:
    """Generate an Eulerian surrogate preserving dinucleotide counts when possible."""

    if len(seq) < 3:
        return _mononucleotide_shuffle(seq, rng)
    edges: dict[str, list[str]] = defaultdict(list)
    for left, right in zip(seq, seq[1:]):
        edges[left].append(right)
    for values in edges.values():
        rng.shuffle(values)
    current = seq[0]
    shuffled = [current]
    for _ in range(len(seq) - 1):
        choices = edges.get(current, [])
        if choices:
            current = choices.pop()
        else:
            remaining = [base for values in edges.values() for base in values]
            if not remaining:
                return _mononucleotide_shuffle(seq, rng)
            current = rng.choice(remaining)
        shuffled.append(current)
    candidate = "".join(shuffled)
    # A malformed fallback must never be advertised as count preserving.
    original = sorted(zip(seq, seq[1:]))
    observed = sorted(zip(candidate, candidate[1:]))
    return candidate if original == observed else _mononucleotide_shuffle(seq, rng)


def _kmer_block_shuffle(seq: str, rng: random.Random, k: int = 3) -> str:
    """Permute non-overlapping k-mer blocks while preserving length and block multiset."""

    if len(seq) <= k:
        return _mononucleotide_shuffle(seq, rng)
    blocks = [seq[index : index + k] for index in range(0, len(seq), k)]
    rng.shuffle(blocks)
    return "".join(blocks)[: len(seq)]


def _null_zscore(observed: float, null_values: list[float]) -> tuple[float, float, float]:
    finite = np.asarray([value for value in null_values if np.isfinite(value)], dtype=float)
    if finite.size < 2:
        return math.nan, math.nan, math.nan
    mean = float(finite.mean())
    std = float(finite.std(ddof=1))
    if not np.isfinite(std) or std <= 1e-12:
        return math.nan, mean, std
    # Positive values mean more compressible than the null family.
    return float((mean - observed) / std), mean, std


def calibrate_compression(
    seq: str,
    *,
    n_surrogates: int = 8,
    seed: int = 13,
    k: int = 3,
) -> dict[str, float | int | str]:
    """Calibrate corrected compression sizes against same-length null families.

    The returned z-scores are descriptive null-standardized statistics, not
    probabilities of function.  Degenerate or too-small null distributions are
    represented as unavailable with ``NaN`` values and an explicit reason.
    """

    cleaned = "".join(base for base in seq.upper() if base in "ACGTN")
    result: dict[str, float | int | str] = {
        "compression_window_size_bp": len(cleaned),
        "compression_null_replicates": int(max(0, n_surrogates)),
        "compression_null_seed": int(seed),
        "compression_calibration_method": "header_corrected_same_length_sequence_surrogates_v1",
        "compression_predictor_version": COMPRESSOR_VERSION,
        "compression_interpretation": "statistical_structure_not_functional_information",
    }
    if len(cleaned) < 12 or n_surrogates < 2:
        reason = "At least 12 bases and two null surrogates are required for compression calibration."
        result.update({"compression_calibration_status": "unavailable", "compression_calibration_reason": reason})
        for method in COMPRESSORS:
            result[f"{method}_header_corrected_ratio"] = math.nan
            result[f"{method}_dinucleotide_null_compression_zscore"] = math.nan
            result[f"{method}_kmer_null_compression_zscore"] = math.nan
        result["multiple_compressor_agreement"] = math.nan
        return result

    rng = random.Random(seed)
    dinucleotide_nulls = [_dinucleotide_shuffle(cleaned, rng) for _ in range(n_surrogates)]
    kmer_nulls = [_kmer_block_shuffle(cleaned, rng, k=k) for _ in range(n_surrogates)]
    zscores: list[float] = []
    unavailable: list[str] = []
    for method in COMPRESSORS:
        observed = header_corrected_compression_ratio(cleaned, method)
        dinuc_values = [header_corrected_compression_ratio(control, method) for control in dinucleotide_nulls]
        kmer_values = [header_corrected_compression_ratio(control, method) for control in kmer_nulls]
        dinuc_z, dinuc_mean, dinuc_std = _null_zscore(observed, dinuc_values)
        kmer_z, kmer_mean, kmer_std = _null_zscore(observed, kmer_values)
        result.update(
            {
                f"{method}_raw_compression_ratio": compressed_size(cleaned, method) / len(cleaned),
                f"{method}_compressor_header_bytes": compressor_header_bytes(method),
                f"{method}_header_corrected_ratio": observed,
                f"{method}_dinucleotide_null_mean": dinuc_mean,
                f"{method}_dinucleotide_null_std": dinuc_std,
                f"{method}_dinucleotide_null_compression_zscore": dinuc_z,
                f"{method}_kmer_null_mean": kmer_mean,
                f"{method}_kmer_null_std": kmer_std,
                f"{method}_kmer_null_compression_zscore": kmer_z,
            }
        )
        for family, value in (("dinucleotide", dinuc_z), ("kmer", kmer_z)):
            if np.isfinite(value):
                zscores.append(float(value))
            else:
                unavailable.append(f"{method}:{family}")

    if zscores:
        signs = np.sign(np.asarray(zscores, dtype=float))
        majority = 1.0 if float(np.sum(signs >= 0)) >= len(signs) / 2 else -1.0
        agreement = float(np.mean(signs == majority))
    else:
        agreement = math.nan
    result["multiple_compressor_agreement"] = agreement
    low_replicate = n_surrogates < 8
    result["compression_calibration_status"] = "available" if not unavailable and not low_replicate else "partial"
    reasons: list[str] = []
    if unavailable:
        reasons.append("Degenerate null variance for: " + ", ".join(unavailable))
    if low_replicate:
        reasons.append(f"Only {n_surrogates} surrogate replicates; increase this count for inferential use.")
    result["compression_calibration_reason"] = (
        "; ".join(reasons)
        if reasons
        else "All requested compressor/null combinations had non-degenerate same-length null distributions."
    )
    return result
