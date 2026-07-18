"""Native-versus-randomized default-state benchmark.

The benchmark measures dependence on the chosen null definition. Native excess
is reported as statistical structure and is never converted into functional
probability or selected-function evidence.
"""

from __future__ import annotations

import math
import random
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from darkdna.evolutionary_null import fit_evolutionary_null, simulate_neutral_sequence
from darkdna.features.sequence import compute_all_sequence_features
from darkdna.nulls.sequence import transform_sequence
from darkdna.nulls.sequence import reverse_complement
from darkdna.views.primitive_scores import score_primitives


DEFAULT_METHODS = (
    "native",
    "whole_genome_reverse",
    "reverse_complement",
    "global_mononucleotide_shuffle",
    "local_mononucleotide_shuffle",
    "global_dinucleotide_shuffle",
    "local_dinucleotide_shuffle",
    "kmer_preserving_shuffle",
    "repeat_only_shuffle",
    "non_repeat_only_shuffle",
    "TE_orientation_reversal",
    "evolutionary_process_generated",
)


def _selected_windows(windows: pd.DataFrame, max_windows: int) -> pd.DataFrame:
    if windows.empty:
        return windows.copy()
    order = windows.sort_values([column for column in ["chrom", "start", "end", "region_id"] if column in windows.columns])
    if len(order) <= max_windows:
        return order.reset_index(drop=True)
    indices = np.linspace(0, len(order) - 1, max_windows, dtype=int)
    return order.iloc[indices].reset_index(drop=True)


def _transform_genomes(
    genome: dict[str, str],
    methods: tuple[str, ...],
    *,
    seed: int,
    kmer_size: int,
    block_size: int,
    repeat_intervals: pd.DataFrame | None,
    te_annotations: pd.DataFrame | None,
) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    transformed: dict[str, dict[str, str]] = {"native": {chrom: sequence.upper() for chrom, sequence in genome.items()}}
    unavailable: dict[str, str] = {}
    evolutionary_model = fit_evolutionary_null(list(genome.values()))
    for method_index, method in enumerate(methods):
        if method == "native" or method.startswith("local_"):
            continue
        if method in {"repeat_only_shuffle", "non_repeat_only_shuffle"}:
            if repeat_intervals is None or repeat_intervals.empty:
                unavailable[method] = "A repeat mask with interval-resolved repeat/non-repeat sequence is required."
                continue
            transformed[method] = {}
            for chrom_index, (chrom, sequence) in enumerate(genome.items()):
                values = list(sequence.upper())
                mask = np.zeros(len(values), dtype=bool)
                same = repeat_intervals.loc[repeat_intervals["chrom"].astype(str) == str(chrom)]
                for interval in same.itertuples():
                    mask[max(0, int(interval.start)) : min(len(values), int(interval.end))] = True
                selected = mask if method == "repeat_only_shuffle" else ~mask
                positions = np.flatnonzero(selected)
                selected_values = [values[position] for position in positions]
                random.Random(seed + method_index * 10_007 + chrom_index).shuffle(selected_values)
                for position, base in zip(positions, selected_values):
                    values[int(position)] = base
                transformed[method][chrom] = "".join(values)
            continue
        if method == "TE_orientation_reversal":
            annotations = te_annotations if te_annotations is not None else repeat_intervals
            if annotations is None or annotations.empty:
                unavailable[method] = "A TE annotation with copy coordinates is required."
                continue
            transformed[method] = {}
            for chrom, sequence in genome.items():
                values = list(sequence.upper())
                same = annotations.loc[annotations["chrom"].astype(str) == str(chrom)].sort_values("start")
                for interval in same.itertuples():
                    start, end = max(0, int(interval.start)), min(len(values), int(interval.end))
                    values[start:end] = list(reverse_complement("".join(values[start:end])))
                transformed[method][chrom] = "".join(values)
            continue
        if method == "evolutionary_process_generated":
            transformed[method] = {
                chrom: simulate_neutral_sequence(evolutionary_model, len(sequence), seed=seed + method_index * 10_007 + chrom_index)
                for chrom_index, (chrom, sequence) in enumerate(genome.items())
            }
            continue
        transformed[method] = {
            chrom: transform_sequence(
                sequence,
                method,
                seed=seed + method_index * 10_007 + chrom_index,
                k=kmer_size,
                block_size=block_size,
            )
            for chrom_index, (chrom, sequence) in enumerate(genome.items())
        }
    return transformed, unavailable


def _feature_rows(
    genome: dict[str, str],
    windows: pd.DataFrame,
    methods: tuple[str, ...],
    transformed: dict[str, dict[str, str]],
    unavailable: dict[str, str],
    *,
    seed: int,
    kmer_size: int,
    block_size: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for method_index, method in enumerate(methods):
        if method in unavailable:
            rows.append(
                {
                    "benchmark_region_id": "",
                    "region_id": "",
                    "null_method": method,
                    "status": "unavailable",
                    "reason": unavailable[method],
                }
            )
            continue
        for row_index, row in windows.iterrows():
            chrom = str(row["chrom"])
            start, end = int(row["start"]), int(row["end"])
            native = genome.get(chrom, "")[start:end]
            if method.startswith("local_"):
                sequence = transform_sequence(
                    native,
                    method,
                    seed=seed + method_index * 10_007 + row_index,
                    k=kmer_size,
                    block_size=block_size,
                )
            else:
                sequence = transformed[method].get(chrom, "")[start:end]
            region_id = str(row.get("region_id", f"{chrom}:{start}-{end}"))
            rows.append(
                {
                    "benchmark_region_id": f"{method}|{region_id}",
                    "region_id": region_id,
                    "null_method": method,
                    "status": "available",
                    "reason": "",
                    **compute_all_sequence_features(sequence),
                }
            )
    return pd.DataFrame(rows)


def _shift_table(table: pd.DataFrame, id_columns: set[str]) -> pd.DataFrame:
    available = table.loc[table["status"] == "available"].copy() if "status" in table.columns else table.copy()
    if available.empty:
        return pd.DataFrame()
    numeric = [
        column
        for column in available.select_dtypes(include=[np.number]).columns
        if column not in id_columns and not column.endswith("_p_value")
    ]
    rows: list[dict[str, object]] = []
    native = available.loc[available["null_method"] == "native"]
    for method, group in available.groupby("null_method", sort=False):
        for feature in numeric:
            native_values = pd.to_numeric(native[feature], errors="coerce").dropna().to_numpy(dtype=float)
            method_values = pd.to_numeric(group[feature], errors="coerce").dropna().to_numpy(dtype=float)
            native_mean = float(np.mean(native_values)) if native_values.size else math.nan
            method_mean = float(np.mean(method_values)) if method_values.size else math.nan
            pooled = float(np.std(np.concatenate([native_values, method_values]), ddof=1)) if native_values.size + method_values.size > 2 else math.nan
            rows.append(
                {
                    "null_method": method,
                    "feature": feature,
                    "native_mean": native_mean,
                    "null_mean": method_mean,
                    "native_minus_null": native_mean - method_mean if np.isfinite(native_mean) and np.isfinite(method_mean) else math.nan,
                    "standardized_shift": (native_mean - method_mean) / pooled if np.isfinite(pooled) and pooled > 0 else math.nan,
                    "native_n": int(native_values.size),
                    "null_n": int(method_values.size),
                    "interpretation": "statistical_structure_shift_not_selected_function",
                }
            )
    return pd.DataFrame(rows)


def benchmark_default_state(
    genome: dict[str, str],
    windows: pd.DataFrame,
    *,
    seed: int = 13,
    methods: tuple[str, ...] = DEFAULT_METHODS,
    max_windows: int = 256,
    local_block_size: int = 1_000,
    kmer_size: int = 3,
    foundation_predictors: dict[str, Callable[[str], float]] | None = None,
    configured_foundation_models: dict[str, str] | None = None,
    repeat_intervals: pd.DataFrame | None = None,
    te_annotations: pd.DataFrame | None = None,
) -> dict[str, pd.DataFrame]:
    selected = _selected_windows(windows, max_windows)
    transformed, unavailable = _transform_genomes(
        genome,
        methods,
        seed=seed,
        kmer_size=kmer_size,
        block_size=local_block_size,
        repeat_intervals=repeat_intervals,
        te_annotations=te_annotations,
    )
    features = _feature_rows(
        genome,
        selected,
        methods,
        transformed,
        unavailable,
        seed=seed,
        kmer_size=kmer_size,
        block_size=local_block_size,
    )
    available_features = features.loc[features["status"] == "available"].copy()
    primitive_input = available_features.drop(columns=["region_id"]).rename(columns={"benchmark_region_id": "region_id"})
    primitive_scores = score_primitives(primitive_input) if not primitive_input.empty else pd.DataFrame()
    if not primitive_scores.empty:
        primitive_scores["null_method"] = primitive_scores["region_id"].str.split("|", n=1).str[0]
        primitive_scores["source_region_id"] = primitive_scores["region_id"].str.split("|", n=1).str[1]
    feature_shift = _shift_table(features, {"length"})
    primitive_shift = _shift_table(primitive_scores, set()) if not primitive_scores.empty else pd.DataFrame()
    sensitivity = (
        feature_shift.loc[feature_shift["null_method"] != "native"]
        .groupby("feature", as_index=False)["standardized_shift"]
        .agg(null_method_shift_min="min", null_method_shift_max="max", null_method_shift_std="std")
        if not feature_shift.empty
        else pd.DataFrame()
    )
    model_rows: list[dict[str, object]] = []
    predictors = foundation_predictors or {}
    configured = configured_foundation_models or {}
    for model in sorted(set(predictors) | set(configured) | {"Puffin_or_Puffin-D", "Enformer", "Borzoi", "AlphaGenome"}):
        if model not in predictors:
            model_rows.append(
                {
                    "model": model,
                    "status": "unavailable",
                    "value": math.nan,
                    "reason": "No local callable adapter was supplied; models are never downloaded automatically.",
                    "local_path": configured.get(model, ""),
                }
            )
            continue
        for chrom, sequence in genome.items():
            model_rows.append(
                {
                    "model": model,
                    "status": "available_local_adapter",
                    "chrom": chrom,
                    "value": float(predictors[model](sequence)),
                    "reason": "User-supplied local adapter inference.",
                    "local_path": configured.get(model, ""),
                }
            )
    return {
        "feature_rows": features,
        "primitive_rows": primitive_scores,
        "feature_shift": feature_shift,
        "primitive_shift": primitive_shift,
        "model_background": pd.DataFrame(model_rows),
        "method_sensitivity": sensitivity,
    }


def write_default_state_benchmark(results: dict[str, pd.DataFrame], outdir: str | Path) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        "feature_shift": out / "native_vs_null_feature_shift.parquet",
        "primitive_shift": out / "native_vs_null_primitive_shift.parquet",
        "model_background": out / "transcription_initiation_background.parquet",
        "method_sensitivity": out / "null_method_sensitivity.parquet",
        "html": out / "default_state_benchmark.html",
    }
    results["feature_shift"].to_parquet(paths["feature_shift"], index=False)
    results["primitive_shift"].to_parquet(paths["primitive_shift"], index=False)
    results["model_background"].to_parquet(paths["model_background"], index=False)
    results["method_sensitivity"].to_parquet(paths["method_sensitivity"], index=False)
    feature_html = results["feature_shift"].head(200).to_html(index=False, escape=True) if not results["feature_shift"].empty else "<p>No calibrated feature shifts.</p>"
    primitive_html = results["primitive_shift"].head(200).to_html(index=False, escape=True) if not results["primitive_shift"].empty else "<p>No primitive shifts.</p>"
    paths["html"].write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>Default-state benchmark</title></head><body>"
        "<h1>Native versus randomized default-state benchmark</h1>"
        "<p>Native-versus-null excess is statistical structure, not proof of selected function. Results depend on null definition.</p>"
        f"<h2>Feature shifts</h2>{feature_html}<h2>Primitive shifts</h2>{primitive_html}</body></html>",
        encoding="utf-8",
    )
    return paths
