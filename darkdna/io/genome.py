"""Genome metadata helpers."""

from __future__ import annotations

from pathlib import Path

from .fasta import genome_sizes_from_fasta, read_chrom_sizes


def load_genome_sizes(fasta: str | Path | None = None, chrom_sizes: str | Path | None = None) -> dict[str, int]:
    if chrom_sizes:
        return read_chrom_sizes(chrom_sizes)
    if fasta:
        return genome_sizes_from_fasta(fasta)
    raise ValueError("Provide either genome FASTA or chrom sizes.")


def scaffold_edge_distance(start: int, end: int, chrom_size: int) -> int:
    return int(min(max(0, start), max(0, chrom_size - end)))


def is_unplaced_or_unlocalized(chrom: str) -> bool:
    text = chrom.lower()
    tokens = ["unplaced", "unlocalized", "random", "un_", "chrunknown", "scaffold", "contig"]
    return any(token in text for token in tokens)
