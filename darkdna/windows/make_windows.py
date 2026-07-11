"""Dark/noncoding window generation."""

from __future__ import annotations

import re
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
from darkdna.utils.progress import ProgressReporter, progress_message
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


def sequence_stat_prefixes(seq: str) -> dict[str, np.ndarray]:
    seq = seq.upper()
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    low_complexity = np.zeros(len(arr), dtype=np.uint8)
    for match in re.finditer(r"([ACGTN])\1{4,}", seq):
        low_complexity[match.start() : match.end()] = 1

    def prefix(mask: np.ndarray) -> np.ndarray:
        out = np.empty(len(mask) + 1, dtype=np.int32)
        out[0] = 0
        out[1:] = np.cumsum(mask, dtype=np.int32)
        return out

    return {
        "gc": prefix((arr == ord("G")) | (arr == ord("C"))),
        "n": prefix(arr == ord("N")),
        "low_complexity": prefix(low_complexity),
    }


def sequence_stats_from_prefixes(prefixes: dict[str, np.ndarray] | None, start: int, end: int) -> tuple[float, float, float]:
    if prefixes is None:
        return np.nan, np.nan, np.nan
    max_len = len(prefixes["gc"]) - 1
    start = max(0, min(int(start), max_len))
    end = max(start, min(int(end), max_len))
    length = max(1, end - start)
    n_count = int(prefixes["n"][end] - prefixes["n"][start])
    gc_count = int(prefixes["gc"][end] - prefixes["gc"][start])
    low_count = int(prefixes["low_complexity"][end] - prefixes["low_complexity"][start])
    usable = max(0, length - n_count)
    gc = float(gc_count / usable) if usable else np.nan
    return gc, float(n_count / length), float(low_count / length)


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


def build_interval_index(df: pd.DataFrame, value_col: str | None = None) -> dict[str, dict[str, object]]:
    if df is None or df.empty or not {"chrom", "start", "end"}.issubset(df.columns):
        return {}
    work = df.copy()
    work["_start_num"] = pd.to_numeric(work["start"], errors="coerce")
    work["_end_num"] = pd.to_numeric(work["end"], errors="coerce")
    work = work.dropna(subset=["_start_num", "_end_num"])
    work["_start_num"] = work["_start_num"].astype(int)
    work["_end_num"] = work["_end_num"].astype(int)
    index: dict[str, dict[str, object]] = {}
    for chrom, group in work.sort_values(["chrom", "_start_num", "_end_num"]).groupby("chrom", sort=False):
        group = group.reset_index(drop=True)
        record: dict[str, object] = {
            "df": group,
            "starts": group["_start_num"].to_numpy(dtype=int),
            "ends": group["_end_num"].to_numpy(dtype=int),
        }
        record["max_end_prefix"] = np.maximum.accumulate(record["ends"])  # type: ignore[arg-type]
        if value_col and value_col in group.columns:
            record["values"] = pd.to_numeric(group[value_col], errors="coerce").fillna(0.0).to_numpy(dtype=float)
        index[str(chrom)] = record
    return index


def _indexed_positions(index: dict[str, dict[str, object]], chrom: str, start: int, end: int) -> tuple[dict[str, object] | None, np.ndarray]:
    record = index.get(str(chrom))
    if record is None:
        return None, np.array([], dtype=int)
    starts = record["starts"]  # type: ignore[index]
    ends = record["ends"]  # type: ignore[index]
    cutoff = int(np.searchsorted(starts, int(end), side="left"))
    if cutoff <= 0:
        return record, np.array([], dtype=int)
    mask = ends[:cutoff] > int(start)
    return record, np.flatnonzero(mask)


def indexed_overlaps_any(index: dict[str, dict[str, object]], chrom: str, start: int, end: int) -> bool:
    record = index.get(str(chrom))
    if record is None:
        return False
    starts = record["starts"]  # type: ignore[index]
    max_end_prefix = record["max_end_prefix"]  # type: ignore[index]
    cutoff = int(np.searchsorted(starts, int(end), side="left"))
    if cutoff <= 0:
        return False
    return bool(max_end_prefix[cutoff - 1] > int(start))


def indexed_overlap_fraction(index: dict[str, dict[str, object]], chrom: str, start: int, end: int) -> float:
    record, positions = _indexed_positions(index, chrom, start, end)
    if record is None or not positions.size:
        return 0.0
    starts = record["starts"][positions]  # type: ignore[index]
    ends = record["ends"][positions]  # type: ignore[index]
    overlaps = np.minimum(ends, int(end)) - np.maximum(starts, int(start))
    return min(1.0, float(np.clip(overlaps, 0, None).sum()) / max(1, int(end) - int(start)))


def indexed_weighted_interval_mean(index: dict[str, dict[str, object]], chrom: str, start: int, end: int) -> float:
    record, positions = _indexed_positions(index, chrom, start, end)
    if record is None or not positions.size or "values" not in record:
        return np.nan
    starts = record["starts"][positions]  # type: ignore[index]
    ends = record["ends"][positions]  # type: ignore[index]
    values = record["values"][positions]  # type: ignore[index]
    overlaps = np.clip(np.minimum(ends, int(end)) - np.maximum(starts, int(start)), 0, None)
    if float(overlaps.sum()) == 0.0:
        return np.nan
    return float((overlaps * values).sum() / max(1, int(end) - int(start)))


def indexed_collect_values(
    index: dict[str, dict[str, object]],
    chrom: str,
    start: int,
    end: int,
    value_col: str,
    max_values: int = 5,
) -> list[str]:
    record, positions = _indexed_positions(index, chrom, start, end)
    if record is None or not positions.size:
        return []
    df = record["df"]  # type: ignore[assignment]
    if value_col not in df.columns:  # type: ignore[union-attr]
        return []
    values: list[str] = []
    for value in df.iloc[positions][value_col].dropna():  # type: ignore[union-attr]
        text = str(value)
        if text and text not in values:
            values.append(text)
        if len(values) >= max_values:
            break
    return values


def nearest_tss_from_index(index: dict[str, dict[str, object]], chrom: str, start: int, end: int) -> tuple[str, float]:
    record = index.get(str(chrom))
    if record is None:
        return "", np.nan
    starts = record["starts"]  # type: ignore[index]
    if len(starts) == 0:
        return "", np.nan
    center = (int(start) + int(end)) // 2
    pos = int(np.searchsorted(starts, center, side="left"))
    candidates = [idx for idx in (pos - 1, pos, pos + 1) if 0 <= idx < len(starts)]
    if not candidates:
        return "", np.nan
    distances = [abs(int(starts[idx]) - center) for idx in candidates]
    best_pos = candidates[int(np.argmin(distances))]
    df = record["df"]  # type: ignore[assignment]
    name = ""
    if "name" in df.columns:  # type: ignore[union-attr]
        name = str(df.iloc[best_pos].get("name", ""))  # type: ignore[union-attr]
    return name, float(abs(int(starts[best_pos]) - center))


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
    progress: bool = False,
) -> pd.DataFrame:
    thresholds = {**DEFAULT_ARTIFACT_THRESHOLDS, **(artifact_thresholds or {})}
    if progress:
        progress_message("make-windows", "loading annotation tracks")
    annotations = annotation_tables(annotation_path, promoter_bp=promoter_bp, progress=progress)
    blacklist = load_optional_bed(blacklist_path)
    te = read_te_annotation(te_annotation_path, progress=progress)
    ccre = load_optional_bed(ccre_path)
    enhancer = load_optional_bed(enhancer_path)
    explicit_promoter = load_optional_bed(promoter_path)
    assembly_gaps = load_optional_bed(assembly_gaps_path)
    segdups = load_optional_bed(segmental_duplication_path)
    centelo = load_optional_bed(centromere_telomere_path)
    mappability = read_bedgraph(mappability_path) if mappability_path else pd.DataFrame(columns=["chrom", "start", "end", "value"])
    if progress:
        progress_message("make-windows", "indexing annotation tracks")
    annotation_indexes = {key: build_interval_index(value) for key, value in annotations.items()}
    blacklist_index = build_interval_index(blacklist)
    te_index = build_interval_index(te)
    ccre_index = build_interval_index(ccre)
    enhancer_index = build_interval_index(enhancer)
    explicit_promoter_index = build_interval_index(explicit_promoter)
    assembly_gaps_index = build_interval_index(assembly_gaps)
    segdups_index = build_interval_index(segdups)
    centelo_index = build_interval_index(centelo)
    mappability_index = build_interval_index(mappability, value_col="value")

    rows = []
    current_chrom = None
    current_prefixes: dict[str, np.ndarray] | None = None
    reporter = ProgressReporter("make-windows annotate", total=len(windows)) if progress else None
    if reporter:
        reporter.start("annotating windows")
    for idx, row in enumerate(windows.itertuples(), start=1):
        chrom, start, end = str(row.chrom), int(row.start), int(row.end)
        if genome is not None and chrom != current_chrom:
            current_chrom = chrom
            current_prefixes = sequence_stat_prefixes(genome.get(chrom, ""))
            if progress:
                progress_message("make-windows", f"prepared sequence stats for {chrom}")
        gc, nfrac, low_complexity = sequence_stats_from_prefixes(current_prefixes, start, end)
        chrom_size = chrom_sizes.get(chrom, end)
        edge = scaffold_edge_distance(start, end, chrom_size)
        mapp = indexed_weighted_interval_mean(mappability_index, chrom, start, end) if mappability_index else np.nan

        overlaps_exon = indexed_overlaps_any(annotation_indexes["exons"], chrom, start, end)
        overlaps_promoter = indexed_overlaps_any(annotation_indexes["promoters"], chrom, start, end) or indexed_overlaps_any(explicit_promoter_index, chrom, start, end)
        overlaps_intron = indexed_overlaps_any(annotation_indexes["introns"], chrom, start, end)
        overlaps_utr = indexed_overlaps_any(annotation_indexes["utrs"], chrom, start, end)
        overlaps_te = indexed_overlaps_any(te_index, chrom, start, end)
        overlaps_blacklist = indexed_overlaps_any(blacklist_index, chrom, start, end)
        overlaps_gap = indexed_overlaps_any(assembly_gaps_index, chrom, start, end)
        overlaps_segdup = indexed_overlaps_any(segdups_index, chrom, start, end)
        overlaps_centelo = indexed_overlaps_any(centelo_index, chrom, start, end)

        nearest_gene, distance_to_nearest_tss = nearest_tss_from_index(annotation_indexes["tss"], chrom, start, end)

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
        repeat_frac = indexed_overlap_fraction(te_index, chrom, start, end)
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
                "TE_family": ";".join(indexed_collect_values(te_index, chrom, start, end, "family")),
                "overlaps_cCRE": indexed_overlaps_any(ccre_index, chrom, start, end),
                "overlaps_enhancer": indexed_overlaps_any(enhancer_index, chrom, start, end),
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
        if reporter:
            reporter.update(idx, message=f"{chrom}:{start}-{end}")
    out = pd.DataFrame(rows)
    if reporter:
        reporter.finish()
    if progress:
        progress_message("make-windows", "assigning multiscale parent/child context")
    out = assign_multiscale_context(out)
    return ensure_region_schema(out)


def make_dark_windows(
    fasta: str | Path | None = None,
    chrom_sizes_path: str | Path | None = None,
    window_sizes: list[int] | None = None,
    step_fraction: float = 0.5,
    exclude_coding_exons: bool = True,
    progress: bool = False,
    **annotation_kwargs,
) -> pd.DataFrame:
    if progress:
        progress_message("make-windows", "loading genome sizes")
    sizes = load_genome_sizes(fasta=fasta, chrom_sizes=chrom_sizes_path)
    mappability_path = annotation_kwargs.get("mappability_path")
    if mappability_path:
        mappability = read_bedgraph(mappability_path)
        if not mappability.empty:
            for chrom, end in mappability.groupby("chrom")["end"].max().items():
                sizes.setdefault(str(chrom), int(end))
    if progress:
        total_bp = sum(sizes.values())
        progress_message("make-windows", f"loaded {len(sizes)} sequences ({total_bp:,} bp)")
        progress_message("make-windows", "loading FASTA sequence")
    genome = read_fasta(fasta) if fasta else None
    window_sizes = window_sizes or [200, 1000, 5000, 10000, 50000]
    rows: list[dict] = []
    reporter = ProgressReporter("make-windows generate", total=len(sizes) * len(window_sizes)) if progress else None
    if reporter:
        reporter.start("generating raw windows")
    for chrom, chrom_size in sizes.items():
        for window_size in window_sizes:
            step = max(1, int(window_size * step_fraction))
            before = len(rows)
            rows.extend(generate_windows_for_chrom(chrom, chrom_size, int(window_size), step))
            if reporter:
                reporter.step(message=f"{chrom} size={window_size} added={len(rows) - before}")
    if reporter:
        reporter.finish(f"raw_windows={len(rows)}")
    if progress:
        progress_message("make-windows", f"raw windows generated: {len(rows):,}")
    windows = pd.DataFrame(rows)
    windows = annotate_windows(windows, genome=genome, chrom_sizes=sizes, progress=progress, **annotation_kwargs)
    if exclude_coding_exons:
        before = len(windows)
        windows = windows.loc[~windows["overlaps_exon"].astype(bool)].copy()
        if progress:
            progress_message("make-windows", f"excluded coding-exon windows: {before - len(windows):,}; retained={len(windows):,}")
    return ensure_region_schema(windows.reset_index(drop=True))


def write_windows(windows: pd.DataFrame, outdir: str | Path) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    parquet = out / "dark_windows.parquet"
    bed = out / "dark_windows.bed"
    windows.to_parquet(parquet, index=False)
    write_bed(windows, bed, columns=["chrom", "start", "end", "region_id"])
    return {"parquet": parquet, "bed": bed}
