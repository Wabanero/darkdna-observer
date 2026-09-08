"""Assign operational primitive labels from residual and matched-null evidence.

Missing z-scores stay NA. A class is assigned only when a finite residual or
matched-null z-score meets its threshold. When several classes survive, every
survivor is emitted; the labeler does not pick a single winner.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.features.classical import artifact_risk_score
from darkdna.primitives.ontology import get_primitive
from darkdna.utils.progress import ProgressReporter
from darkdna.utils.stats import optional_float
from .assay_recommender import recommend_assay


UNEXPLAINED = "unexplained_dark_anomaly_candidate"

DOMINANCE_FEATURES = {
    "fractal_scaffold_candidate": ["multiscale_texture_screening_score", "DFA_surrogate_zscore", "multiscale_parent_child_similarity_screen"],
    "constraint_grammar_region_candidate": ["grammar_entropy", "Markov_order_anomaly", "motif_like_token_recurrence"],
    "non_B_DNA_physical_susceptibility_candidate": ["G4_sequence_potential", "Z_DNA_sequence_potential", "R_loop_susceptibility_sequence_potential", "charge_oxidation_susceptibility_score"],
    "replication_instability_candidate": ["fork_texture_score", "simple_repeat_fraction", "palindrome_density"],
    "periodic_spacing_grammar_candidate": ["phase_periodicity_around_10bp", "spacing_periodicity_fourier_power", "spacing_periodicity_autocorrelation"],
    "asymmetric_repeat_architecture_candidate": ["left_right_GC_asymmetry", "nested_repeat_architecture_score", "orientation_bias_of_recurrent_kmers"],
    "negative_space_element_candidate": ["depleted_kmer_score", "unexpected_silence_score"],
    "sequence_regime_boundary_candidate": ["boundary_condition_candidate_score", "left_right_regime_difference_score"],
    "TE_grammar_node_candidate": ["TE_family_mosaic_score", "TE_boundary_score", "TE_overlap_fraction"],
}


def _finite_number(value: object) -> float:
    return optional_float(value)


def _survives_thresholds(rz: float, nz: float, residual_threshold: float, matched_null_threshold: float) -> bool:
    residual_ok = np.isfinite(rz) and rz >= residual_threshold
    null_ok = np.isfinite(nz) and nz >= matched_null_threshold
    return bool(residual_ok or null_ok)


def _rank_score(rz: float, nz: float) -> float:
    if not np.isfinite(rz) and not np.isfinite(nz):
        return math.nan
    score = 0.0
    if np.isfinite(rz):
        score += rz
    if np.isfinite(nz):
        score += 0.5 * nz
    return score


def _priority(rz: float, nz: float, artifact_score: float, severe_null_support: bool) -> float:
    parts: list[float] = []
    if np.isfinite(rz):
        parts.append(max(rz, 0.0))
    if np.isfinite(nz):
        parts.append(max(nz, 0.0))
    if not parts:
        return math.nan
    priority = min(1.0, sum(parts) / 8.0) * (1.0 - 0.35 * artifact_score)
    if not severe_null_support:
        priority *= 0.5
    return float(priority)


def _supporting_features(feature_row: pd.Series | None, primitive: str) -> list[str]:
    if feature_row is None:
        return []
    features = DOMINANCE_FEATURES.get(primitive, [])
    scored = []
    for name in features:
        if name in feature_row.index:
            number = _finite_number(feature_row[name])
            if np.isfinite(number):
                scored.append((name, number))
    return [name for name, _ in sorted(scored, key=lambda item: abs(item[1]), reverse=True)]


def _promotion_fields(record: pd.Series, nz: float, matched_null_threshold: float, minimum_null_models_for_promotion: int, minimum_null_agreement_for_promotion: float) -> dict:
    null_count_value = record.get("null_model_count", np.nan)
    null_count = int(null_count_value) if pd.notna(null_count_value) else 0
    agreement_value = record.get("null_model_agreement", np.nan)
    null_agreement = _finite_number(agreement_value)
    conflict_value = record.get("null_model_conflict", False)
    null_conflict = bool(conflict_value) if pd.notna(conflict_value) else False
    null_panel = str(record.get("null_panel_status", "") or "")
    severe_null_support = bool(
        null_panel == "severe_null_panel_available"
        and null_count >= minimum_null_models_for_promotion
        and np.isfinite(null_agreement)
        and null_agreement >= minimum_null_agreement_for_promotion
        and not null_conflict
        and np.isfinite(nz)
        and nz >= matched_null_threshold
    )
    if severe_null_support:
        promotion_status = "eligible_for_candidate_promotion"
    elif not null_panel:
        promotion_status = "screening_only_legacy_null_metadata_unavailable"
    elif null_conflict:
        promotion_status = "screening_only_conflicting_null_families"
    else:
        promotion_status = "screening_only_insufficient_severe_null_support"
    return {
        "null_count": null_count,
        "null_agreement": null_agreement,
        "null_conflict": null_conflict,
        "null_panel": null_panel,
        "severe_null_support": severe_null_support,
        "promotion_status": promotion_status,
    }


def _label_row(
    *,
    region_id: object,
    primitive: str,
    score_name: object,
    rz: float,
    nz: float,
    empirical_p: float,
    empirical_p_status: str,
    promotion: dict,
    priority: float,
    supporting: list[str],
    flags: object,
    labeling_status: str,
    competing: list[str],
) -> dict:
    exclusive = labeling_status == "single_surviving_class" and len(competing) == 1
    return {
        "region_id": region_id,
        "primitive_class": primitive,
        "primitive_score_name": score_name,
        "primitive_priority": priority,
        "primitive_priority_status": "uncalibrated_ranking_priority_not_probability",
        "primitive_confidence": priority,
        "primitive_confidence_deprecation_warning": (
            "Deprecated alias for primitive_priority; this 0-1 ranking heuristic is not calibrated confidence or probability."
        ),
        "residual_zscore": rz,
        "matched_null_zscore": nz,
        "empirical_p_value": empirical_p,
        "empirical_p_value_status": empirical_p_status,
        "null_panel_status": promotion["null_panel"],
        "null_model_count": promotion["null_count"],
        "null_model_agreement": promotion["null_agreement"],
        "null_model_conflict": promotion["null_conflict"],
        "survives_severe_null_panel": promotion["severe_null_support"],
        "candidate_promotion_status": promotion["promotion_status"],
        "top_supporting_features": ";".join(supporting),
        "artifact_risk_flags": flags,
        "recommended_assay": recommend_assay(primitive).get("assay", ""),
        "labeling_status": labeling_status,
        "competing_primitive_classes": ";".join(competing),
        "competing_primitive_count": int(len(competing)),
        "is_exclusive_label": bool(exclusive),
    }


def assign_primitive_labels(
    residuals: pd.DataFrame,
    features: pd.DataFrame | None = None,
    windows: pd.DataFrame | None = None,
    residual_threshold: float = 2.0,
    matched_null_threshold: float = 2.0,
    minimum_null_models_for_promotion: int = 3,
    minimum_null_agreement_for_promotion: float = 0.5,
    *,
    progress: bool = False,
) -> pd.DataFrame:
    feature_lookup = features.set_index("region_id") if features is not None and not features.empty else pd.DataFrame()
    window_lookup = windows.set_index("region_id") if windows is not None and not windows.empty else pd.DataFrame()
    rows = []
    reporter = ProgressReporter("infer-primitives", total=int(residuals["region_id"].nunique())) if progress else None
    if reporter:
        reporter.start("assigning primitive labels")
    for idx, (region_id, group) in enumerate(residuals.groupby("region_id"), start=1):
        feature_row = feature_lookup.loc[region_id] if not feature_lookup.empty and region_id in feature_lookup.index else None
        window_row = window_lookup.loc[region_id] if not window_lookup.empty and region_id in window_lookup.index else None
        flags = window_row.get("artifact_risk_flags", "") if window_row is not None else ""
        artifact_score = artifact_risk_score(flags)
        survivors: list[dict] = []
        metadata_fallback: dict | None = None
        fallback_rank = -math.inf
        for _, record in group.iterrows():
            primitive = get_primitive(str(record["primitive"])).candidate_name
            rz = _finite_number(record.get("residual_zscore"))
            nz = _finite_number(record.get("matched_null_zscore"))
            empirical_p = _finite_number(record.get("empirical_p_value"))
            promotion = _promotion_fields(
                record,
                nz,
                matched_null_threshold,
                minimum_null_models_for_promotion,
                minimum_null_agreement_for_promotion,
            )
            candidate = {
                "primitive": primitive,
                "score_name": record["primitive"],
                "rz": rz,
                "nz": nz,
                "empirical_p": empirical_p,
                "empirical_p_status": str(record.get("empirical_p_value_status", "")),
                "promotion": promotion,
                "priority": _priority(rz, nz, artifact_score, promotion["severe_null_support"]),
                "supporting": _supporting_features(feature_row, primitive),
                "rank": _rank_score(rz, nz),
            }
            rank = candidate["rank"]
            if np.isfinite(rank) and rank > fallback_rank:
                fallback_rank = rank
                metadata_fallback = candidate
            elif metadata_fallback is None:
                metadata_fallback = candidate
            if _survives_thresholds(rz, nz, residual_threshold, matched_null_threshold):
                survivors.append(candidate)
        collapsed: dict[str, dict] = {}
        for item in survivors:
            current = collapsed.get(item["primitive"])
            if current is None or (np.isfinite(item["rank"]) and (not np.isfinite(current["rank"]) or item["rank"] > current["rank"])):
                collapsed[item["primitive"]] = item
        survivors = list(collapsed.values())
        specific = [item for item in survivors if item["primitive"] != UNEXPLAINED]
        if specific:
            selected = sorted(specific, key=lambda item: (-item["rank"] if np.isfinite(item["rank"]) else math.inf, item["primitive"]))
            labeling_status = "single_surviving_class" if len(selected) == 1 else "competing_hypotheses_not_exclusive"
            competing = [item["primitive"] for item in selected]
        elif survivors:
            selected = [item for item in survivors if item["primitive"] == UNEXPLAINED]
            labeling_status = "unexplained_residual_without_dominant_class"
            competing = [UNEXPLAINED]
        else:
            fallback = metadata_fallback or {
                "primitive": "no_call",
                "score_name": "",
                "rz": math.nan,
                "nz": math.nan,
                "empirical_p": math.nan,
                "empirical_p_status": "",
                "promotion": {
                    "null_count": 0,
                    "null_agreement": math.nan,
                    "null_conflict": False,
                    "null_panel": "",
                    "severe_null_support": False,
                    "promotion_status": "screening_only_legacy_null_metadata_unavailable",
                },
                "priority": math.nan,
                "supporting": [],
            }
            selected = [fallback]
            labeling_status = "no_call"
            competing = []
            selected[0] = {**fallback, "primitive": "no_call", "supporting": []}
        for item in selected:
            rows.append(
                _label_row(
                    region_id=region_id,
                    primitive=item["primitive"] if labeling_status != "no_call" else "no_call",
                    score_name=item["score_name"],
                    rz=item["rz"],
                    nz=item["nz"],
                    empirical_p=item["empirical_p"],
                    empirical_p_status=item["empirical_p_status"],
                    promotion=item["promotion"],
                    priority=item["priority"] if labeling_status != "no_call" else math.nan,
                    supporting=item["supporting"] if labeling_status != "no_call" else [],
                    flags=flags,
                    labeling_status=labeling_status,
                    competing=competing,
                )
            )
        if reporter:
            reporter.update(idx, message=str(region_id))
    if reporter:
        reporter.finish()
    return pd.DataFrame(rows)


def write_primitive_labels(labels: pd.DataFrame, outdir: str | Path) -> Path:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "primitive_labels.parquet"
    labels.to_parquet(path, index=False)
    alias = out / "candidate_primitives.parquet"
    labels.to_parquet(alias, index=False)
    return path
