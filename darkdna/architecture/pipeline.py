"""Mode B orchestration and output contracts."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from darkdna.architecture.amount_features import compute_amount_features
from darkdna.architecture.comparative_length import compute_length_conservation
from darkdna.architecture.comparison import compare_sequence_vs_quantity
from darkdna.architecture.copy_number import compute_copy_number_features
from darkdna.architecture.heterochromatin_mass import compute_heterochromatin_mass
from darkdna.architecture.null_models import build_architecture_nulls
from darkdna.architecture.occupancy import compute_occupancy_features
from darkdna.architecture.presence_absence import compute_presence_absence_features
from darkdna.architecture.scoring import architecture_score_manifest, score_architecture_features
from darkdna.architecture.sequence_indifference import Predictor, evaluate_sequence_indifference
from darkdna.architecture.spacing import compute_spacing_features


def _normalise_intervals(intervals: pd.DataFrame) -> pd.DataFrame:
    frame = intervals.copy()
    if "region_id" not in frame.columns:
        if "locus_id" in frame.columns:
            frame["region_id"] = frame["locus_id"].astype(str)
        else:
            frame["region_id"] = frame.apply(lambda row: f"{row['chrom']}:{int(row['start'])}-{int(row['end'])}", axis=1)
    return frame


def run_sequence_indifferent_architecture(
    intervals: pd.DataFrame,
    *,
    sequences: dict[str, str],
    mode_a: pd.DataFrame | None = None,
    repeats: pd.DataFrame | None = None,
    copy_number: pd.DataFrame | None = None,
    presence_absence: pd.DataFrame | None = None,
    syntenic_intervals: pd.DataFrame | None = None,
    anchors: pd.DataFrame | None = None,
    occupancy: pd.DataFrame | None = None,
    heterochromatin: pd.DataFrame | None = None,
    replication_domains: pd.DataFrame | None = None,
    phenotype: pd.DataFrame | None = None,
    genome_sizes: dict[str, int] | None = None,
    predictor: Predictor | None = None,
    seed: int = 13,
    kmer_size: int = 3,
    block_size_bp: int = 100_000,
    minimum_independent_blocks: int = 5,
) -> dict[str, pd.DataFrame]:
    frame = _normalise_intervals(intervals)
    amount = compute_amount_features(
        frame,
        repeats=repeats,
        heterochromatin=heterochromatin,
        replication_domains=replication_domains,
        genome_sizes=genome_sizes,
    )
    copy_features = compute_copy_number_features(frame, copy_number, phenotype)
    presence_features = compute_presence_absence_features(frame, presence_absence)
    length_features = compute_length_conservation(frame, syntenic_intervals)
    spacing_features = compute_spacing_features(frame, anchors)
    occupancy_features = compute_occupancy_features(frame, occupancy)
    heterochromatin_features = compute_heterochromatin_mass(frame, heterochromatin)
    feature_tables = [amount, copy_features, presence_features, length_features, spacing_features, occupancy_features, heterochromatin_features]
    architecture_features = feature_tables[0]
    for table in feature_tables[1:]:
        architecture_features = architecture_features.merge(table, on="region_id", how="left")
    if "artifact_risk_flags" in frame.columns:
        architecture_features = architecture_features.merge(frame[["region_id", "artifact_risk_flags"]], on="region_id", how="left")
    pool = [sequence for sequence in sequences.values() if sequence]
    summaries, controls = [], []
    spacing_lookup = spacing_features.set_index("region_id") if not spacing_features.empty else pd.DataFrame()
    for index, region_id in enumerate(frame["region_id"].astype(str)):
        sequence = sequences.get(region_id, "")
        if not sequence:
            summaries.append(
                {
                    "region_id": region_id,
                    "sequence_identity_sensitivity": float("nan"),
                    "length_sensitivity": float("nan"),
                    "copy_number_sensitivity": float("nan"),
                    "sequence_indifference_score": float("nan"),
                    "sequence_quantity_interaction_score": float("nan"),
                    "sequence_indifference_status": "unavailable_missing_sequence",
                    "sequence_indifference_reason": "No FASTA sequence was available for this interval.",
                }
            )
            continue
        spacing_value = None
        if isinstance(spacing_lookup, pd.DataFrame) and region_id in spacing_lookup.index:
            spacing_value = pd.to_numeric(pd.Series([spacing_lookup.loc[region_id].get("anchor_to_anchor_distance")]), errors="coerce").iloc[0]
        summary, region_controls = evaluate_sequence_indifference(
            region_id,
            sequence,
            predictor=predictor,
            replacement_pool=pool,
            spacing=float(spacing_value) if pd.notna(spacing_value) else None,
            seed=seed + index * 10_007,
            k=kmer_size,
        )
        summaries.append(summary)
        controls.extend(region_controls)
    sensitivity = pd.DataFrame(summaries)
    control_table = pd.DataFrame(controls)
    candidates = score_architecture_features(architecture_features, sensitivity)
    nulls, residuals = build_architecture_nulls(
        candidates,
        architecture_features,
        control_table,
        block_size_bp=block_size_bp,
        minimum_independent_blocks=minimum_independent_blocks,
    )
    matched_nulls = nulls.loc[nulls["null_model_id"] == "matched_interval_independent_blocks"]
    candidates = candidates.merge(
        matched_nulls[["region_id", "architecture_null_zscore", "architecture_null_empirical_p", "architecture_null_status"]],
        on="region_id",
        how="left",
    )
    candidates["promotion_status"] = candidates["architecture_null_status"].map(
        lambda status: "eligible_for_candidate_promotion" if status == "available_block_calibrated" else "screening_only_not_for_promotion"
    )
    comparison = compare_sequence_vs_quantity(mode_a if mode_a is not None else pd.DataFrame(), candidates)
    return {
        "architecture_features": architecture_features,
        "copy_number_features": copy_features,
        "presence_absence_features": presence_features,
        "spacing_features": spacing_features,
        "length_conservation": length_features,
        "sequence_indifference_controls": control_table,
        "sequence_indifference_summary": sensitivity,
        "sequence_vs_quantity_scores": comparison,
        "architecture_nulls": nulls,
        "architecture_residuals": residuals,
        "architecture_candidates": candidates,
    }


def write_architecture_outputs(results: dict[str, pd.DataFrame], outdir: str | Path) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    names = {
        "architecture_features": "architecture_features.parquet",
        "copy_number_features": "copy_number_features.parquet",
        "presence_absence_features": "presence_absence_features.parquet",
        "spacing_features": "spacing_features.parquet",
        "length_conservation": "length_conservation.parquet",
        "sequence_indifference_controls": "sequence_indifference_controls.parquet",
        "sequence_indifference_summary": "sequence_indifference_summary.parquet",
        "sequence_vs_quantity_scores": "sequence_vs_quantity_scores.parquet",
        "architecture_nulls": "architecture_nulls.parquet",
        "architecture_residuals": "architecture_residuals.parquet",
        "architecture_candidates": "architecture_candidates.parquet",
    }
    paths: dict[str, Path] = {}
    for key, filename in names.items():
        path = out / filename
        results.get(key, pd.DataFrame()).to_parquet(path, index=False)
        paths[key] = path
    manifest = out / "architecture_score_manifest.json"
    manifest.write_text(json.dumps(architecture_score_manifest(), indent=2), encoding="utf-8")
    paths["manifest"] = manifest
    return paths
