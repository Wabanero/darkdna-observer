"""Matched null model generation and summaries.

The current implementation provides a first matched-control null. It also
emits a registry of required complementary nulls so downstream reports do not
mistake one z-score for a severe null panel.
"""

from __future__ import annotations

import json
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.utils.progress import ProgressReporter
from darkdna.utils.stats import empirical_p_value
from darkdna.views.primitive_scores import PRIMITIVE_SCORE_COLUMNS


NULL_MODEL_REGISTRY: list[dict[str, str]] = [
    {
        "null_model_id": "matched_controls_v1",
        "status": "implemented",
        "description": "Nearby/similar genomic controls matched on available covariates.",
    },
    {
        "null_model_id": "same_length_gc_matched",
        "status": "partially_available_via_matched_controls",
        "description": "Controls with similar length/window size and GC content.",
    },
    {
        "null_model_id": "dinucleotide_preserving_shuffle",
        "status": "helper_available_not_pipeline_integrated",
        "description": "Synthetic sequence shuffle preserving much dinucleotide texture.",
    },
    {
        "null_model_id": "kmer_preserving_shuffle",
        "status": "helper_available_not_pipeline_integrated",
        "description": "Synthetic sequence shuffle preserving local k-mer blocks.",
    },
    {
        "null_model_id": "te_family_age_matched",
        "status": "requires_tracks",
        "description": "Controls matched on TE family and TE age/divergence annotations.",
    },
    {
        "null_model_id": "chromatin_compartment_matched",
        "status": "requires_tracks",
        "description": "Controls matched on A/B compartment or comparable chromatin state.",
    },
    {
        "null_model_id": "replication_timing_matched",
        "status": "requires_tracks",
        "description": "Controls matched on replication timing signal.",
    },
    {
        "null_model_id": "recombination_mutation_environment_matched",
        "status": "requires_tracks",
        "description": "Controls matched on recombination, mutation-rate, or damage environment.",
    },
    {
        "null_model_id": "gene_tss_distance_matched",
        "status": "implemented_if_distance_to_nearest_tss_available",
        "description": "Controls matched on gene/TSS proximity.",
    },
    {
        "null_model_id": "nearby_genomic_controls",
        "status": "partially_available_via_matched_controls",
        "description": "Local genomic controls rather than genome-wide random controls.",
    },
    {
        "null_model_id": "syntenic_ortholog_controls",
        "status": "not_implemented",
        "description": "Orthologous or syntenic sequence controls across assemblies/species.",
    },
    {
        "null_model_id": "population_frequency_controls",
        "status": "not_implemented",
        "description": "Controls matched on allele frequency, presence/absence, or copy number.",
    },
    {
        "null_model_id": "reversed_or_synthetic_sequences",
        "status": "not_pipeline_integrated",
        "description": "Reverse, scrambled, or synthetic controls for transformation robustness.",
    },
]


def null_model_registry() -> list[dict[str, str]]:
    return [dict(item) for item in NULL_MODEL_REGISTRY]


def null_panel_status() -> dict:
    implemented = [item["null_model_id"] for item in NULL_MODEL_REGISTRY if item["status"] == "implemented"]
    not_fully_implemented = [
        item["null_model_id"]
        for item in NULL_MODEL_REGISTRY
        if item["status"] != "implemented"
    ]
    return {
        "status": "insufficient_single_matched_null_until_complementary_nulls_pass",
        "implemented_null_models": implemented,
        "missing_or_partial_null_models": not_fully_implemented,
    }


def gc_matched_shuffle(seq: str, seed: int = 13) -> str:
    rng = random.Random(seed)
    chars = list(seq.upper())
    rng.shuffle(chars)
    return "".join(chars)


def dinucleotide_preserving_shuffle(seq: str, seed: int = 13) -> str:
    """Lightweight Eulerian-ish fallback preserving much dinucleotide texture."""

    seq = seq.upper()
    if len(seq) < 3:
        return gc_matched_shuffle(seq, seed)
    rng = random.Random(seed)
    edges: dict[str, list[str]] = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        edges[a].append(b)
    for values in edges.values():
        rng.shuffle(values)
    current = seq[0]
    out = [current]
    for _ in range(len(seq) - 1):
        if edges.get(current):
            nxt = edges[current].pop()
        else:
            nxt = rng.choice(list("ACGTN"))
        out.append(nxt)
        current = nxt
    return "".join(out)


def kmer_preserving_shuffle(seq: str, k: int = 3, seed: int = 13) -> str:
    seq = seq.upper()
    if len(seq) <= k:
        return gc_matched_shuffle(seq, seed)
    rng = random.Random(seed)
    blocks = [seq[i : i + k] for i in range(0, len(seq), k)]
    rng.shuffle(blocks)
    return "".join(blocks)[: len(seq)]


def matched_feature_columns(features: pd.DataFrame) -> list[str]:
    candidates = [
        "gc_content",
        "length",
        "window_size",
        "mappability",
        "simple_repeat_fraction",
        "CpG_density",
        "TE_overlap_fraction",
        "TE_family",
        "TE_age",
        "TE_divergence",
        "local_TE_density",
        "distance_to_nearest_tss",
        "chromatin_compartment",
        "replication_timing",
        "recombination_rate",
        "mutation_rate",
        "copy_number",
        "presence_absence_frequency",
        "population_frequency",
        "chrom",
    ]
    return [col for col in candidates if col in features.columns]


def _numeric_distance(row: pd.Series, controls: pd.DataFrame, cols: list[str]) -> pd.Series:
    dist = pd.Series(0.0, index=controls.index)
    for col in cols:
        if col == "chrom":
            dist += (controls[col].astype(str) != str(row[col])).astype(float) * 0.25
            continue
        values = pd.to_numeric(controls[col], errors="coerce")
        value = pd.to_numeric(pd.Series([row[col]]), errors="coerce").iloc[0]
        scale = values.std()
        if pd.isna(value) or pd.isna(scale) or scale == 0:
            continue
        dist += ((values - value).abs() / scale).fillna(0.0)
    return dist


def select_matched_controls(features: pd.DataFrame, region_id: str, n: int = 25) -> pd.DataFrame:
    if features.empty:
        return features
    row = features.loc[features["region_id"] == region_id]
    if row.empty:
        return features.head(0)
    row = row.iloc[0]
    controls = features.loc[features["region_id"] != region_id].copy()
    if controls.empty:
        return controls
    cols = matched_feature_columns(features)
    dist = _numeric_distance(row, controls, cols)
    return controls.loc[dist.sort_values().head(n).index]


def build_matched_null_models(
    scores: pd.DataFrame,
    features: pd.DataFrame | None = None,
    n_controls: int = 25,
    *,
    progress: bool = False,
) -> pd.DataFrame:
    features = features if features is not None and not features.empty else scores[["region_id"]].copy()
    if "region_id" not in features.columns:
        features = scores[["region_id"]].copy()
    merged_scores = scores.set_index("region_id")
    rows = []
    available_score_cols = [col for col in PRIMITIVE_SCORE_COLUMNS if col in scores.columns]
    region_ids = scores["region_id"].astype(str).tolist()
    reporter = ProgressReporter("build-null-models", total=len(region_ids)) if progress else None
    if reporter:
        reporter.start(f"matching controls n_controls={n_controls}")
    panel_status = null_panel_status()
    registry_ids = [item["null_model_id"] for item in NULL_MODEL_REGISTRY]
    missing_or_partial = panel_status["missing_or_partial_null_models"]
    for idx, region_id in enumerate(region_ids, start=1):
        controls = select_matched_controls(features, region_id, n=n_controls)
        if controls.empty:
            controls = features.loc[features["region_id"].astype(str) != region_id]
        if controls.empty:
            controls = features
        control_ids = [rid for rid in controls["region_id"].astype(str).tolist() if rid in merged_scores.index]
        for primitive in available_score_cols:
            observed = float(merged_scores.loc[region_id, primitive])
            null_values = merged_scores.loc[control_ids, primitive].astype(float).to_numpy() if control_ids else scores[primitive].astype(float).to_numpy()
            null_values = null_values[np.isfinite(null_values)]
            null_mean = float(np.mean(null_values)) if null_values.size else np.nan
            null_std = float(np.std(null_values, ddof=1)) if null_values.size > 1 else 0.0
            null_z = 0.0 if null_std == 0 or not np.isfinite(null_std) else (observed - null_mean) / null_std
            rows.append(
                {
                    "null_model_id": "matched_controls_v1",
                    "region_id": region_id,
                    "primitive": primitive,
                    "primitive_score": observed,
                    "null_mean": null_mean,
                    "null_std": null_std,
                    "null_zscore": float(null_z),
                    "empirical_p_value": empirical_p_value(observed, null_values, higher=True),
                    "matched_features_used": ",".join(matched_feature_columns(features)),
                    "null_panel_status": panel_status["status"],
                    "available_null_models": ",".join(registry_ids),
                    "missing_or_partial_null_models": ",".join(missing_or_partial),
                }
            )
        if reporter:
            reporter.update(idx, message=region_id)
    if reporter:
        reporter.finish()
    return pd.DataFrame(rows)


def write_matched_nulls(nulls: pd.DataFrame, outdir: str | Path) -> Path:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "matched_nulls.parquet"
    nulls.to_parquet(path, index=False)
    alias = out / "null_model_summary.parquet"
    nulls.to_parquet(alias, index=False)
    registry_path = out / "null_model_registry.json"
    registry_path.write_text(json.dumps(null_model_registry(), indent=2), encoding="utf-8")
    return path
