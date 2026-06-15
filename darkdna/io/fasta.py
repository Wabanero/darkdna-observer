"""FASTA readers and sequence extraction.

The core reader is deliberately small and pure-Python. Optional FASTA packages
can accelerate production use, but the MVP should work anywhere.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


def read_fasta(path: str | Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records.setdefault(current, [])
            elif current is not None:
                records[current].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def write_fasta(records: dict[str, str], path: str | Path, width: int = 80) -> None:
    with Path(path).open("w", encoding="utf-8") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n")
            for idx in range(0, len(seq), width):
                handle.write(seq[idx : idx + width] + "\n")


def read_chrom_sizes(path: str | Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                parts = line.split()
            sizes[parts[0]] = int(parts[1])
    return sizes


def write_chrom_sizes(sizes: dict[str, int], path: str | Path) -> None:
    with Path(path).open("w", encoding="utf-8") as handle:
        for chrom, size in sizes.items():
            handle.write(f"{chrom}\t{size}\n")


def genome_sizes_from_fasta(path: str | Path) -> dict[str, int]:
    return {chrom: len(seq) for chrom, seq in read_fasta(path).items()}


def fetch_sequence(genome: dict[str, str], chrom: str, start: int, end: int) -> str:
    seq = genome.get(chrom, "")
    start = max(0, int(start))
    end = min(len(seq), int(end))
    if end <= start:
        return ""
    return seq[start:end].upper()


def add_sequences_to_windows(windows: pd.DataFrame, fasta: str | Path) -> pd.DataFrame:
    genome = read_fasta(fasta)
    out = windows.copy()
    out["sequence"] = [
        fetch_sequence(genome, str(row.chrom), int(row.start), int(row.end)) for row in out.itertuples()
    ]
    return out


def write_sequences_for_windows(windows: pd.DataFrame, fasta: str | Path, output_fasta: str | Path) -> None:
    genome = read_fasta(fasta)
    records = {
        str(row.region_id): fetch_sequence(genome, str(row.chrom), int(row.start), int(row.end))
        for row in windows.itertuples()
    }
    write_fasta(records, output_fasta)


def iter_fasta_windows(fasta: str | Path, windows: pd.DataFrame) -> Iterable[tuple[str, str]]:
    genome = read_fasta(fasta)
    for row in windows.itertuples():
        yield str(row.region_id), fetch_sequence(genome, str(row.chrom), int(row.start), int(row.end))
