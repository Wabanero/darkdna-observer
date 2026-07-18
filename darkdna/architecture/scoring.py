"""Conservative Mode B scoring and operational candidate classification."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.architecture.schemas import EVOLUTIONARY_CAVEAT, MODEL_CAVEAT


QUANTITY_COLUMNS = (
    "interval_length",
    "repeat_array_length",
    "repeat_family_total_bp",
    "local_repeat_fraction",
    "heterochromatin_overlap",
    "replication_domain_length",
    "candidate_family_burden",
    "adjacent_copy_count",
    "inter_copy_spacing",
    "copy_number_variance",
    "presence_absence_frequency",
    "length_conservation_score",
    "spacing_null_zscore",
    "physical_occupancy",
)


def _robust_abs_z(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    median = float(numeric.median())
    mad = float((numeric - median).abs().median())
    if not np.isfinite(mad) or mad == 0:
        return pd.Series(np.nan, index=values.index)
    return (numeric - median).abs() / (1.4826 * mad)


def _candidate_label(row: pd.Series) -> tuple[str, str]:
    artifact = str(row.get("artifact_risk_flags", "") or "")
    if artifact:
        return "artifact_compatible_candidate", "artifact"
    identity = float(row.get("sequence_identity_sensitivity", math.nan))
    length = float(row.get("length_sensitivity", math.nan))
    copy = float(row.get("copy_number_sensitivity", math.nan))
    spacing = abs(float(row.get("spacing_null_zscore", math.nan))) / 3.0 if np.isfinite(float(row.get("spacing_null_zscore", math.nan))) else math.nan
    hetero = float(row.get("heterochromatin_overlap", math.nan))
    occupancy = float(row.get("physical_occupancy", math.nan))
    quantity_values = [value for value in (length, copy, spacing) if np.isfinite(value)]
    quantity = max(quantity_values) if quantity_values else math.nan
    high_identity = np.isfinite(identity) and identity >= 0.20
    high_quantity = np.isfinite(quantity) and quantity >= 0.20
    if high_identity and high_quantity:
        return "mixed_sequence_quantity_candidate", "mixed"
    if high_identity and not high_quantity:
        return "sequence_specific_architecture_candidate", "Mode_A"
    if not high_identity and high_quantity:
        if np.isfinite(copy) and copy == quantity:
            return "copy_number_constraint_candidate", "Mode_B"
        if np.isfinite(spacing) and spacing == quantity:
            return "spacing_constraint_candidate", "Mode_B"
        return "length_constrained_interval_candidate", "Mode_B"
    if np.isfinite(hetero) and hetero >= 0.5:
        return "bulk_heterochromatin_candidate", "Mode_B"
    if np.isfinite(occupancy) and occupancy > 0:
        return "occupancy_dependent_region_candidate", "Mode_B"
    return "unresolved_architecture_candidate", "unresolved"


def score_architecture_features(features: pd.DataFrame, sensitivities: pd.DataFrame) -> pd.DataFrame:
    frame = features.merge(sensitivities, on="region_id", how="left")
    z_columns = []
    for column in QUANTITY_COLUMNS:
        if column in frame.columns and pd.to_numeric(frame[column], errors="coerce").notna().sum() >= 3:
            z_column = f"{column}_descriptive_robust_abs_z"
            frame[z_column] = _robust_abs_z(frame[column])
            z_columns.append(z_column)
    frame["quantity_anomaly_score"] = frame[z_columns].max(axis=1, skipna=True) if z_columns else math.nan
    labels = frame.apply(_candidate_label, axis=1)
    frame["sequence_indifferent_candidate"] = [label for label, _ in labels]
    frame["dominant_mode"] = [mode for _, mode in labels]
    frame["mode_b_score_status"] = "descriptive_and_model_based_screen_not_causal_probability"
    frame["candidate_only"] = True
    frame["allowed_interpretation"] = "Prioritizes intervals for equal-length replacement, length titration, copy-number titration, or spacer assays."
    frame["forbidden_interpretation"] = "Do not infer selected function, adaptive origin, or biological causality from this score."
    frame["model_based_caveat"] = MODEL_CAVEAT
    frame["evolutionary_caveat"] = EVOLUTIONARY_CAVEAT
    return frame


def architecture_score_manifest() -> dict[str, object]:
    return {
        "mode": "Mode_B_sequence_indifferent_architecture",
        "score_status": "screening_axis_not_probability",
        "quantity_columns": list(QUANTITY_COLUMNS),
        "classification_thresholds": {"identity_sensitivity": 0.20, "quantity_sensitivity": 0.20},
        "calibration": "Block-aware architecture nulls are reported separately; descriptive robust z-scores are not p-values.",
        "caveats": [MODEL_CAVEAT, EVOLUTIONARY_CAVEAT],
    }


def write_architecture_manifest(path: str | Path) -> Path:
    target = Path(path)
    target.write_text(json.dumps(architecture_score_manifest(), indent=2), encoding="utf-8")
    return target
