"""Assign operational primitive labels from residual and matched-null evidence."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.features.classical import artifact_risk_score
from darkdna.utils.progress import ProgressReporter
from .assay_recommender import recommend_assay


SCORE_TO_PRIMITIVE = {
    "fractal_scaffold_candidate_score": "fractal_scaffold_candidate",
    "constraint_grammar_region_candidate_score": "constraint_grammar_region_candidate",
    "non_B_DNA_physical_susceptibility_candidate_score": "non_B_DNA_physical_susceptibility_candidate",
    "quantum_susceptible_domain_candidate_score": "non_B_DNA_physical_susceptibility_candidate",
    "replication_instability_candidate_score": "replication_instability_candidate",
    "decoherence_boundary_candidate_score": "decoherence_boundary_candidate",
    "resonant_pulse_decoder_candidate_score": "resonant_pulse_decoder_candidate",
    "hysteresis_candidate_score": "hysteresis_candidate",
    "possibility_gate_candidate_score": "possibility_gate_candidate",
    "criticality_tuner_candidate_score": "criticality_tuner_candidate",
    "chromatin_motion_oscillator_candidate_score": "chromatin_motion_oscillator_candidate",
    "negative_space_element_candidate_score": "negative_space_element_candidate",
    "sequence_regime_boundary_candidate_score": "sequence_regime_boundary_candidate",
    "TE_grammar_node_candidate_score": "TE_grammar_node_candidate",
    "unexplained_dark_anomaly_candidate_score": "unexplained_dark_anomaly_candidate",
}

DOMINANCE_FEATURES = {
    "fractal_scaffold_candidate": ["multiscale_texture_screening_score", "DFA_surrogate_zscore", "multiscale_parent_child_similarity_screen"],
    "constraint_grammar_region_candidate": ["grammar_entropy", "Markov_order_anomaly", "motif_like_token_recurrence"],
    "non_B_DNA_physical_susceptibility_candidate": ["G4_sequence_potential", "Z_DNA_sequence_potential", "R_loop_susceptibility_sequence_potential", "charge_oxidation_susceptibility_score"],
    "replication_instability_candidate": ["fork_texture_score", "simple_repeat_fraction", "palindrome_density"],
    "decoherence_boundary_candidate": ["decoherence_boundary_candidate_score", "entropy_boundary_score"],
    "resonant_pulse_decoder_candidate": ["phase_periodicity_around_10bp", "spacing_periodicity_fourier_power"],
    "hysteresis_candidate": ["left_right_GC_asymmetry", "nested_repeat_architecture_score", "G4_sequence_potential", "Z_DNA_sequence_potential"],
    "possibility_gate_candidate": ["boundary_condition_candidate_score", "negative_space_boundary_score", "forbidden_word_depletion_enrichment"],
    "criticality_tuner_candidate": ["entropy_boundary_score", "compression_boundary_score", "local_feature_transition_score"],
    "chromatin_motion_oscillator_candidate": ["spacing_periodicity_autocorrelation", "DNA_bendability_proxy", "left_right_entropy_asymmetry"],
    "negative_space_element_candidate": ["depleted_kmer_score", "unexpected_silence_score"],
    "sequence_regime_boundary_candidate": ["boundary_condition_candidate_score", "left_right_regime_difference_score"],
    "TE_grammar_node_candidate": ["TE_family_mosaic_score", "TE_boundary_score", "TE_overlap_fraction"],
}


def _best_residual_for_region(residual_rows: pd.DataFrame) -> pd.Series:
    ranked = residual_rows.assign(rank_score=residual_rows["residual_zscore"].fillna(0) + residual_rows["matched_null_zscore"].fillna(0) * 0.5)
    return ranked.sort_values("rank_score", ascending=False).iloc[0]


def _supporting_features(feature_row: pd.Series | None, primitive: str) -> list[str]:
    if feature_row is None:
        return []
    features = DOMINANCE_FEATURES.get(primitive, [])
    scored = []
    for name in features:
        if name in feature_row.index:
            try:
                scored.append((name, float(feature_row[name])))
            except Exception:
                scored.append((name, 0.0))
    return [name for name, _ in sorted(scored, key=lambda item: abs(item[1]), reverse=True) if name]


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
        best = _best_residual_for_region(group)
        primitive = SCORE_TO_PRIMITIVE.get(best["primitive"], "unexplained_dark_anomaly_candidate")
        rz_value = best.get("residual_zscore", np.nan)
        rz = float(rz_value) if pd.notna(rz_value) and np.isfinite(float(rz_value)) else 0.0
        nz = float(best.get("matched_null_zscore", 0.0) if pd.notna(best.get("matched_null_zscore", np.nan)) else 0.0)
        null_count_value = best.get("null_model_count", 0)
        null_count = int(null_count_value) if pd.notna(null_count_value) else 0
        agreement_value = best.get("null_model_agreement", np.nan)
        null_agreement = float(agreement_value) if pd.notna(agreement_value) else np.nan
        conflict_value = best.get("null_model_conflict", False)
        null_conflict = bool(conflict_value) if pd.notna(conflict_value) else False
        null_panel = str(best.get("null_panel_status", "") or "")
        severe_null_support = bool(
            null_panel == "severe_null_panel_available"
            and null_count >= minimum_null_models_for_promotion
            and np.isfinite(null_agreement)
            and null_agreement >= minimum_null_agreement_for_promotion
            and not null_conflict
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
        if rz < residual_threshold and nz < matched_null_threshold:
            primitive = "unexplained_dark_anomaly_candidate" if float(group["residual_zscore"].max()) >= residual_threshold else "no_call"
        feature_row = feature_lookup.loc[region_id] if not feature_lookup.empty and region_id in feature_lookup.index else None
        window_row = window_lookup.loc[region_id] if not window_lookup.empty and region_id in window_lookup.index else None
        flags = window_row.get("artifact_risk_flags", "") if window_row is not None else ""
        artifact_score = artifact_risk_score(flags)
        priority = max(0.0, min(1.0, (max(rz, 0.0) + max(nz, 0.0)) / 8.0)) * (1.0 - 0.35 * artifact_score)
        if not severe_null_support:
            priority *= 0.5
        supporting = _supporting_features(feature_row, primitive)
        rows.append(
            {
                "region_id": region_id,
                "primitive_class": primitive,
                "primitive_score_name": best["primitive"],
                "primitive_priority": float(priority),
                "primitive_priority_status": "uncalibrated_ranking_priority_not_probability",
                "primitive_confidence": float(priority),
                "primitive_confidence_deprecation_warning": (
                    "Deprecated alias for primitive_priority; this 0-1 ranking heuristic is not calibrated confidence or probability."
                ),
                "residual_zscore": rz,
                "matched_null_zscore": nz,
                "empirical_p_value": (
                    float(best.get("empirical_p_value"))
                    if pd.notna(best.get("empirical_p_value", np.nan))
                    else np.nan
                ),
                "empirical_p_value_status": str(best.get("empirical_p_value_status", "")),
                "null_panel_status": null_panel,
                "null_model_count": null_count,
                "null_model_agreement": null_agreement,
                "null_model_conflict": null_conflict,
                "survives_severe_null_panel": severe_null_support,
                "candidate_promotion_status": promotion_status,
                "top_supporting_features": ";".join(supporting),
                "artifact_risk_flags": flags,
                "recommended_assay": recommend_assay(primitive).get("assay", ""),
            }
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
