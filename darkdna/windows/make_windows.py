"""Dark/noncoding window generation."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from darkdna.io.bed import (
    collect_overlapping_values,
    concatenate_flags,
    overlap_fraction,
    overlaps_any,
    read_bed,
    read_bedgraph,
    weighted_interval_mean,
    write_bed,
)
from darkdna.io.fasta import read_fasta
from darkdna.io.genome import is_unplaced_or_unlocalized, load_genome_sizes, scaffold_edge_distance
from darkdna.io.gff import annotation_tables, read_te_annotation
from darkdna.io.validators import ensure_region_schema
from darkdna.utils.config import DEFAULT_ARTIFACT_THRESHOLDS
from .multiscale import assign_multiscale_context, scale_level_for_size


def simple_gc(seq: str) -> float:
    clean = [b for b in seq.upper() if b in "ACGT"]
    if not clean:
        return np.nan
    return float((clean.count("G") + clean.count("C")) / len(clean))


def n_fraction(seq: str) -> float:
    if not seq:
        return 1.0
    return float(seq.upper().count("N") / len(seq))


def low_complexity_fraction(seq: str) -> float:
    if not seq:
        return 1.0
    seq = seq.upper()
    runs = 0
    current = ""
    length = 0
    for base in seq:
        if base == current:
            length += 1
        else:
            if length >= 5:
                runs += length
            current = base
            length = 1
    if length >= 5:
        runs += length
    return float(runs / len(seq))


def generate_windows_for_chrom(chrom: str, size: int, window_size: int, step: int) -> Iterable[dict]:
    if size <= 0:
        return
    for start in range(0, max(1, size - window_size + 1), max(1, step)):
        end = min(size, start + window_size)
        if end <= start:
            continue
        yield {"chrom": chrom, "start": start, "end": end, "window_size": window_size}
    if size < window_size:
        yield {"chrom": chrom, "start": 0, "end": size, "window_size": window_size}
    elif (size - window_size) % max(1, step) != 0:
        yield {"chrom": chrom, "start": max(0, size - window_size), "end": size, "window_size": window_size}


def load_optional_bed(path: str | Path | None) -> pd.DataFrame:
    return read_bed(path) if path else pd.DataFrame(columns=["chrom", "start", "end"])


def annotate_windows(
    windows: pd.DataFrame,
    genome: dict[str, str] | None,
    chrom_sizes: dict[str, int],
    annotation_path: str | Path | None = None,
    blacklist_path: str | Path | None = None,
    te_annotation_path: str | Path | None = None,
    ccre_path: str | Path | None = None,
    enhancer_path: str | Path | None = None,
    promoter_path: str | Path | None = None,
    mappability_path: str | Path | None = None,
    assembly_gaps_path: str | Path | None = None,
    segmental_duplication_path: str | Path | None = None,
    centromere_telomere_path: str | Path | None = None,
    promoter_bp: int = 1000,
    artifact_thresholds: dict[str, float] | None = None,
) -> pd.DataFrame:
    thresholds = {**DEFAULT_ARTIFACT_THRESHOLDS, **(artifact_thresholds or {})}
    annotations = annotation_tables(annotation_path, promoter_bp=promoter_bp)
    blacklist = load_optional_bed(blacklist_path)
    te = read_te_annotation(te_annotation_path)
    ccre = load_optional_bed(ccre_path)
    enhancer = load_optional_bed(enhancer_path)
    explicit_promoter = load_optional_bed(promoter_path)
    assembly_gaps = load_optional_bed(assembly_gaps_path)
    segdups = load_optional_bed(segmental_duplication_path)
    centelo = load_optional_bed(centromere_telomere_path)
    mappability = read_bedgraph(mappability_path) if mappability_path else pd.DataFrame(columns=["chrom", "start", "end", "value"])

    rows = []
    for row in windows.itertuples():
        chrom, start, end = str(row.chrom), int(row.start), int(row.end)
        seq = genome.get(chrom, "")[start:end] if genome is not None else ""
        gc = simple_gc(seq) if seq else np.nan
        nfrac = n_fraction(seq) if seq else np.nan
        low_complexity = low_complexity_fraction(seq) if seq else np.nan
        chrom_size = chrom_sizes.get(chrom, end)
        edge = scaffold_edge_distance(start, end, chrom_size)
        mapp = weighted_interval_mean(chrom, start, end, mappability) if not mappability.empty else np.nan

        overlaps_exon = overlaps_any(chrom, start, end, annotations["exons"])
        overlaps_promoter = overlaps_any(chrom, start, end, annotations["promoters"]) or overlaps_any(chrom, start, end, explicit_promoter)
        overlaps_intron = overlaps_any(chrom, start, end, annotations["introns"])
        overlaps_utr = overlaps_any(chrom, start, end, annotations["utrs"])
        overlaps_te = overlaps_any(chrom, start, end, te)
        overlaps_blacklist = overlaps_any(chrom, start, end, blacklist)
        overlaps_gap = overlaps_any(chrom, start, end, assembly_gaps)
        overlaps_segdup = overlaps_any(chrom, start, end, segdups)
        overlaps_centelo = overlaps_any(chrom, start, end, centelo)

        tss = annotations["tss"]
        nearest_gene = ""
        distance_to_nearest_tss = np.nan
        if not tss.empty:
            same = tss[tss["chrom"].astype(str) == chrom]
            center = (start + end) // 2
            if not same.empty:
                distances = (same["start"].astype(int) - center).abs()
                idx = distances.idxmin()
                nearest_gene = str(same.loc[idx, "name"])
                distance_to_nearest_tss = float(distances.loc[idx])

        flags = []
        if np.isfinite(nfrac) and nfrac >= thresholds["high_n_fraction"]:
            flags.append("high_N_fraction")
        if np.isfinite(mapp) and mapp <= thresholds["low_mappability"]:
            flags.append("low_mappability")
        if overlaps_blacklist:
            flags.append("blacklist_overlap")
        if overlaps_gap:
            flags.append("assembly_gap_overlap")
        if overlaps_segdup:
            flags.append("segmental_duplication_overlap")
        if overlaps_centelo:
            flags.append("centromeric_telomeric_or_proximal_repeat_context")
        if np.isfinite(low_complexity) and low_complexity >= thresholds["extreme_low_complexity"]:
            flags.append("extreme_low_complexity")
        repeat_frac = overlap_fraction(chrom, start, end, te)
        if repeat_frac >= thresholds["extreme_repeat_density"]:
            flags.append("extreme_repeat_density")
        if edge <= int(thresholds["scaffold_edge_bp"]):
            flags.append("scaffold_edge_proximity")
        usable = 1.0 - (nfrac if np.isfinite(nfrac) else 0.0)
        if usable < thresholds["very_short_usable_fraction"]:
            flags.append("very_short_usable_sequence_after_masking")
        if is_unplaced_or_unlocalized(chrom):
            flags.append("unplaced_or_unlocalized_contig")

        rows.append(
            {
                **row._asdict(),
                "region_id": f"{chrom}:{start}-{end}:{int(row.window_size)}",
                "scale_level": scale_level_for_size(int(row.window_size), sorted(windows["window_size"].unique())),
                "parent_region_id": None,
                "child_region_ids": "",
                "is_dark": not overlaps_exon,
                "overlaps_exon": bool(overlaps_exon),
                "overlaps_promoter": bool(overlaps_promoter),
                "overlaps_intron": bool(overlaps_intron),
                "overlaps_utr": bool(overlaps_utr),
                "overlaps_TE": bool(overlaps_te),
                "TE_family": ";".join(collect_overlapping_values(chrom, start, end, te, "family")),
                "overlaps_cCRE": overlaps_any(chrom, start, end, ccre),
                "overlaps_enhancer": overlaps_any(chrom, start, end, enhancer),
                "overlaps_blacklist": bool(overlaps_blacklist),
                "overlaps_assembly_gap": bool(overlaps_gap),
                "overlaps_segmental_duplication": bool(overlaps_segdup),
                "nearest_gene": nearest_gene,
                "distance_to_nearest_tss": distance_to_nearest_tss,
                "gc_content": gc,
                "n_fraction": nfrac,
                "mappability": mapp,
                "low_complexity_mask_fraction": low_complexity,
                "scaffold_edge_distance": edge,
                "artifact_risk_flags": concatenate_flags(flags),
            }
        )
    out = pd.DataFrame(rows)
    out = assign_multiscale_context(out)
    return ensure_region_schema(out)


def make_dark_windows(
    fasta: str | Path | None = None,
    chrom_sizes_path: str | Path | None = None,
    window_sizes: list[int] | None = None,
    step_fraction: float = 0.5,
    exclude_coding_exons: bool = True,
    **annotation_kwargs,
) -> pd.DataFrame:
    sizes = load_genome_sizes(fasta=fasta, chrom_sizes=chrom_sizes_path)
    genome = read_fasta(fasta) if fasta else None
    window_sizes = window_sizes or [200, 1000, 5000, 10000, 50000]
    rows: list[dict] = []
    for chrom, chrom_size in sizes.items():
        for window_size in window_sizes:
            step = max(1, int(window_size * step_fraction))
            rows.extend(generate_windows_for_chrom(chrom, chrom_size, int(window_size), step))
    windows = pd.DataFrame(rows)
    windows = annotate_windows(windows, genome=genome, chrom_sizes=sizes, **annotation_kwargs)
    if exclude_coding_exons:
        windows = windows.loc[~windows["overlaps_exon"].astype(bool)].copy()
    return ensure_region_schema(windows.reset_index(drop=True))


def write_windows(windows: pd.DataFrame, outdir: str | Path) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    parquet = out / "dark_windows.parquet"
    bed = out / "dark_windows.bed"
    windows.to_parquet(parquet, index=False)
    write_bed(windows, bed, columns=["chrom", "start", "end", "region_id"])
    return {"parquet": parquet, "bed": bed}
