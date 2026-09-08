"""Matched genomic controls plus candidate-level sequence-transform nulls.

Genomic matching is a covariate control. Reverse, dinucleotide, k-mer, and
evolutionary-process families test whether the exact interval sequence is
unusual and enter promotion when they are fully calibrated.
"""

from __future__ import annotations

import json
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.nulls.calibration import build_severe_null_panel
from darkdna.nulls.registry import assess_null_availability, null_model_registry as severe_null_model_registry
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
        "status": "candidate_level_when_sequence_available",
        "description": "Synthetic sequence shuffle preserving much dinucleotide texture.",
    },
    {
        "null_model_id": "kmer_preserving_shuffle",
        "status": "candidate_level_when_sequence_available",
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
        "status": "candidate_level_when_sequence_available",
        "description": "Reverse, scrambled, or synthetic controls for transformation robustness.",
    },
]


def null_model_registry() -> list[dict[str, str]]:
    severe = [dict(item) for item in severe_null_model_registry()]
    severe.append(
        {
            "null_model_id": "matched_controls_v1",
            "family": "legacy_alias",
            "execution_mode": "matched_table",
            "required_columns": [],
            "required_input": "",
            "description": "Deprecated compatibility alias for the earlier single matched-control implementation.",
            "status": "deprecated_alias",
        }
    )
    return severe


def null_panel_status() -> dict:
    assessed = assess_null_availability(())
    implemented = [
        "matched_controls_v1",
        "dinucleotide_preserving_shuffle",
        "kmer_preserving_shuffle",
        "mononucleotide_preserving",
        "markov_chain_surrogate",
        "reversed_sequence",
        "reverse_complement",
        "synthetic_equal_composition",
        "evolutionary_process_generated",
    ]
    not_fully_implemented = [str(item["null_model_id"]) for item in assessed if not item["available"]]
    if "syntenic_ortholog_controls" not in not_fully_implemented:
        not_fully_implemented.append("syntenic_ortholog_controls")
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
    block_size_bp: int = 100_000,
    minimum_independent_blocks: int = 5,
    agreement_z_threshold: float = 2.0,
    sequences: dict[str, str] | None = None,
    n_sequence_surrogates: int = 8,
    seed: int = 13,
    kmer_size: int = 3,
    progress: bool = False,
) -> pd.DataFrame:
    features = features if features is not None and not features.empty else scores[["region_id"]].copy()
    if "region_id" not in features.columns:
        features = scores[["region_id"]].copy()
    reporter = ProgressReporter("build-null-models", total=1) if progress else None
    if reporter:
        reporter.start(f"severe block-aware panel n_controls={n_controls} sequence_surrogates={n_sequence_surrogates}")
    summary, details = build_severe_null_panel(
        scores,
        features,
        n_controls=n_controls,
        block_size_bp=block_size_bp,
        minimum_independent_blocks=minimum_independent_blocks,
        agreement_z_threshold=agreement_z_threshold,
        score_columns=[column for column in PRIMITIVE_SCORE_COLUMNS if column in scores.columns],
        sequences=sequences,
        n_sequence_surrogates=n_sequence_surrogates,
        seed=seed,
        kmer_size=kmer_size,
    )
    summary.attrs["null_details"] = details
    if reporter:
        reporter.finish(f"summary_rows={len(summary)} detail_rows={len(details)}")
    return summary


def write_matched_nulls(nulls: pd.DataFrame, outdir: str | Path) -> Path:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "matched_nulls.parquet"
    details = nulls.attrs.get("null_details")
    serializable_nulls = nulls.copy()
    serializable_nulls.attrs = {}
    serializable_nulls.to_parquet(path, index=False)
    alias = out / "null_model_summary.parquet"
    serializable_nulls.to_parquet(alias, index=False)
    if isinstance(details, pd.DataFrame):
        details.to_parquet(out / "severe_null_details.parquet", index=False)
    registry_path = out / "null_model_registry.json"
    registry_path.write_text(json.dumps(null_model_registry(), indent=2), encoding="utf-8")
    return path
