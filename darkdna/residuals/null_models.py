"""Matched null model generation and summaries."""

from __future__ import annotations

import random
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.utils.stats import empirical_p_value
from darkdna.views.primitive_scores import PRIMITIVE_SCORE_COLUMNS


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
        "local_TE_density",
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
) -> pd.DataFrame:
    features = features if features is not None and not features.empty else scores[["region_id"]].copy()
    if "region_id" not in features.columns:
        features = scores[["region_id"]].copy()
    merged_scores = scores.set_index("region_id")
    rows = []
    available_score_cols = [col for col in PRIMITIVE_SCORE_COLUMNS if col in scores.columns]
    for region_id in scores["region_id"].astype(str):
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
                }
            )
    return pd.DataFrame(rows)


def write_matched_nulls(nulls: pd.DataFrame, outdir: str | Path) -> Path:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "matched_nulls.parquet"
    nulls.to_parquet(path, index=False)
    alias = out / "null_model_summary.parquet"
    nulls.to_parquet(alias, index=False)
    return path
