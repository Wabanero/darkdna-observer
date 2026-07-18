"""Covariance-aware primitive screening views.

Candidate scores are cohort-standardized screening statistics, not
probabilities and not null significance.  Historical equal-weight means are
retained only in ``*_legacy_screening_composite`` columns.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.utils.progress import ProgressReporter
from darkdna.utils.stats import robust_scale_series
from .boundary_conditions import compute_boundary_condition_view
from .entropy_noise import compute_entropy_noise_view
from .physical_susceptibility import compute_physical_susceptibility_view
from .scale_fractal import compute_scale_fractal_features
from .unexplained_anomaly import cross_fitted_unexplained_outlierness, infer_block_groups


PRIMITIVE_SCORE_COLUMNS = [
    "fractal_scaffold_candidate_score",
    "constraint_grammar_region_candidate_score",
    "non_B_DNA_physical_susceptibility_candidate_score",
    "replication_instability_candidate_score",
    "decoherence_boundary_candidate_score",
    "resonant_pulse_decoder_candidate_score",
    "hysteresis_candidate_score",
    "possibility_gate_candidate_score",
    "criticality_tuner_candidate_score",
    "chromatin_motion_oscillator_candidate_score",
    "negative_space_element_candidate_score",
    "sequence_regime_boundary_candidate_score",
    "TE_grammar_node_candidate_score",
    "unexplained_dark_anomaly_candidate_score",
]


PRIMITIVE_SCORE_COMPONENTS: dict[str, list[str]] = {
    "fractal_scaffold_candidate_score": [
        "multiscale_texture_screening_score",
        "multiscale_parent_child_similarity_screen",
        "gzip_dinucleotide_null_compression_zscore",
        "gzip_kmer_null_compression_zscore",
    ],
    "constraint_grammar_region_candidate_score": [
        "grammar_entropy",
        "forbidden_word_depletion_enrichment",
        "motif_like_token_recurrence",
        "Markov_order_anomaly",
        "long_range_dependency_proxy",
    ],
    "non_B_DNA_physical_susceptibility_candidate_score": [
        "G4_sequence_potential",
        "i_motif_sequence_potential",
        "Z_DNA_sequence_potential",
        "triplex_H_DNA_sequence_potential",
        "cruciform_sequence_potential",
        "slipped_DNA_sequence_potential",
        "R_loop_susceptibility_sequence_potential",
    ],
    "replication_instability_candidate_score": ["fork_texture_score"],
    "decoherence_boundary_candidate_score": ["decoherence_boundary_candidate_score"],
    "resonant_pulse_decoder_candidate_score": [
        "phase_periodicity_around_10bp",
        "nucleosome_scale_periodicity_around_147bp",
        "spacing_periodicity_fourier_power",
    ],
    "hysteresis_candidate_score": [
        "abs(left_right_GC_asymmetry)",
        "nested_repeat_architecture_score",
        "max_nonb_sequence_potential",
        "orientation_bias_of_recurrent_kmers",
    ],
    "possibility_gate_candidate_score": [
        "boundary_condition_candidate_score",
        "negative_space_element_candidate_score",
        "forbidden_word_depletion_enrichment",
    ],
    "criticality_tuner_candidate_score": [
        "boundary_condition_candidate_score",
        "entropy_boundary_score",
        "compression_boundary_score",
    ],
    "chromatin_motion_oscillator_candidate_score": [
        "spacing_periodicity_autocorrelation",
        "legacy_dinucleotide_bendability_screen",
        "left_right_entropy_asymmetry",
    ],
    "negative_space_element_candidate_score": [
        "depleted_kmer_score",
        "unexpected_silence_score",
        "local_feature_void_score",
    ],
    "sequence_regime_boundary_candidate_score": ["boundary_condition_candidate_score"],
    "TE_grammar_node_candidate_score": [
        "TE_overlap_fraction",
        "TE_family_mosaic_score",
        "TE_boundary_score",
    ],
    "unexplained_dark_anomaly_candidate_score": [
        "cross_fitted_multivariate_outlierness_after_known_screening_axes",
    ],
}


FEATURE_ORIENTATION_AUDIT: dict[str, dict[str, str | int]] = {
    feature: {
        "orientation": 1,
        "meaning": "higher values increase the corresponding candidate-screen priority after robust cohort standardization",
    }
    for components in PRIMITIVE_SCORE_COMPONENTS.values()
    for feature in components
    if not feature.startswith("cross_fitted_")
}


FEATURE_ALIASES = {
    "multiscale_texture_screening_score": "fractal_score",
    "multiscale_parent_child_similarity_screen": "scale_persistence_score",
    "legacy_dinucleotide_bendability_screen": "DNA_bendability_proxy",
    "G4_sequence_potential": "G4_susceptibility_proxy",
    "Z_DNA_sequence_potential": "Z_DNA_propensity_proxy",
    "triplex_H_DNA_sequence_potential": "triplex_H_DNA_proxy",
    "cruciform_sequence_potential": "cruciform_forming_potential",
    "R_loop_susceptibility_sequence_potential": "R_loop_forming_potential",
}


SCREENING_COMPOSITE_CAVEAT = (
    "Primitive score columns are covariance-aware, robustly cohort-standardized screening views. "
    "They are not probabilities, not p-values, not selected mechanistic measurements, and not evidence of function. "
    "An explicit null panel and held-out validation remain necessary."
)


REQUIRED_COMPOSITE_VALIDATION = [
    "feature_orientation_audit",
    "feature_specific_null_calibration",
    "correlation_and_double_counting_audit",
    "covariance_aware_evidence_combination",
    "cross_validated_weighting_when_labels_exist",
    "test_stability_across_organisms",
    "test_robustness_to_sequence_transformations",
    "reproduce_on_held_out_chromosomes_or_blocks",
]


def primitive_score_manifest() -> dict:
    return {
        "score_status": "covariance_aware_cohort_standardized_screening_score",
        "legacy_score_status": "deprecated_equal_weight_legacy_screening_composite",
        "calibration_status": "cohort_standardized_not_null_calibrated_not_probability",
        "caveat": SCREENING_COMPOSITE_CAVEAT,
        "interpretation_order": [
            "measured_feature_profile",
            "statistical_anomaly_after_explicit_controls_and_nulls",
            "post_hoc_mechanistic_hypothesis",
        ],
        "required_validation_before_mechanistic_use": REQUIRED_COMPOSITE_VALIDATION,
        "components": PRIMITIVE_SCORE_COMPONENTS,
        "feature_orientation_audit": FEATURE_ORIENTATION_AUDIT,
        "combination_method": "robust_component_zscores_combined_with_observed_correlation_covariance_denominator",
        "unexplained_anomaly_method": "cross_fitted_shrinkage_robust_mahalanobis_on_known_screening_axes",
        "empirical_p_value_policy": "unavailable_without_explicit_null_distribution",
    }


def value(row: dict, key: str, default: float = math.nan) -> float:
    candidate_keys = [key]
    if key in FEATURE_ALIASES:
        candidate_keys.append(FEATURE_ALIASES[key])
    for candidate in candidate_keys:
        try:
            raw = row.get(candidate, math.nan)
            if raw is None or pd.isna(raw):
                continue
            return float(raw)
        except Exception:
            continue
    return default


def _finite_mean(values: list[float]) -> float:
    finite = np.asarray([number for number in values if np.isfinite(number)], dtype=float)
    return float(finite.mean()) if finite.size else math.nan


def _finite_max(values: list[float]) -> float:
    finite = np.asarray([number for number in values if np.isfinite(number)], dtype=float)
    return float(finite.max()) if finite.size else math.nan


def _enrich_row(row: dict) -> dict:
    enriched = dict(row)
    if "multiscale_texture_screening_score" not in enriched and "fractal_score" not in enriched and "sequence" in enriched:
        enriched.update(compute_scale_fractal_features(str(enriched["sequence"])))
    enriched.update(compute_physical_susceptibility_view(enriched))
    enriched.update(compute_entropy_noise_view(enriched))
    enriched.update(compute_boundary_condition_view(enriched))
    enriched["max_nonb_sequence_potential"] = _finite_max(
        [
            value(enriched, "G4_sequence_potential"),
            value(enriched, "i_motif_sequence_potential"),
            value(enriched, "Z_DNA_sequence_potential"),
            value(enriched, "triplex_H_DNA_sequence_potential"),
            value(enriched, "cruciform_sequence_potential"),
            value(enriched, "slipped_DNA_sequence_potential"),
            value(enriched, "R_loop_susceptibility_sequence_potential"),
        ]
    )
    return enriched


def primitive_scores_for_row(row: dict) -> dict[str, float]:
    """Return deprecated row-local legacy composites for compatibility.

    Cohort-standardized canonical scores are created by :func:`score_primitives`.
    The unexplained-anomaly score is deliberately unavailable here because it
    requires held-out multivariate calibration.
    """

    enriched = _enrich_row(row)
    physical_values = [value(enriched, feature) for feature in PRIMITIVE_SCORE_COMPONENTS["non_B_DNA_physical_susceptibility_candidate_score"]]
    scores = {
        "fractal_scaffold_candidate_score": _finite_mean([value(enriched, feature) for feature in PRIMITIVE_SCORE_COMPONENTS["fractal_scaffold_candidate_score"]]),
        "constraint_grammar_region_candidate_score": _finite_mean([value(enriched, feature) for feature in PRIMITIVE_SCORE_COMPONENTS["constraint_grammar_region_candidate_score"]]),
        "non_B_DNA_physical_susceptibility_candidate_score": _finite_max(physical_values),
        "replication_instability_candidate_score": value(enriched, "fork_texture_score"),
        "decoherence_boundary_candidate_score": value(enriched, "decoherence_boundary_candidate_score"),
        "resonant_pulse_decoder_candidate_score": _finite_mean([value(enriched, feature) for feature in PRIMITIVE_SCORE_COMPONENTS["resonant_pulse_decoder_candidate_score"]]),
        "hysteresis_candidate_score": _finite_mean(
            [
                abs(value(enriched, "left_right_GC_asymmetry")),
                value(enriched, "nested_repeat_architecture_score"),
                value(enriched, "max_nonb_sequence_potential"),
                value(enriched, "orientation_bias_of_recurrent_kmers"),
            ]
        ),
        "possibility_gate_candidate_score": _finite_mean([value(enriched, feature) for feature in PRIMITIVE_SCORE_COMPONENTS["possibility_gate_candidate_score"]]),
        "criticality_tuner_candidate_score": _finite_mean([value(enriched, feature) for feature in PRIMITIVE_SCORE_COMPONENTS["criticality_tuner_candidate_score"]]),
        "chromatin_motion_oscillator_candidate_score": _finite_mean([value(enriched, feature) for feature in PRIMITIVE_SCORE_COMPONENTS["chromatin_motion_oscillator_candidate_score"]]),
        "negative_space_element_candidate_score": _finite_mean([value(enriched, feature) for feature in PRIMITIVE_SCORE_COMPONENTS["negative_space_element_candidate_score"]]),
        "sequence_regime_boundary_candidate_score": value(enriched, "boundary_condition_candidate_score"),
        "TE_grammar_node_candidate_score": _finite_mean([value(enriched, feature) for feature in PRIMITIVE_SCORE_COMPONENTS["TE_grammar_node_candidate_score"]]),
        "unexplained_dark_anomaly_candidate_score": math.nan,
    }
    return {name: float(number) for name, number in scores.items()}


def _component_series(enriched: pd.DataFrame, expression: str) -> pd.Series:
    absolute = expression.startswith("abs(") and expression.endswith(")")
    key = expression[4:-1] if absolute else expression
    candidate = key if key in enriched.columns else FEATURE_ALIASES.get(key)
    if candidate is None or candidate not in enriched.columns:
        return pd.Series(np.nan, index=enriched.index, dtype=float)
    series = pd.to_numeric(enriched[candidate], errors="coerce")
    return series.abs() if absolute else series


def _robust_component_zscore(series: pd.Series) -> tuple[pd.Series, str]:
    finite = series.notna() & np.isfinite(series)
    output = pd.Series(np.nan, index=series.index, dtype=float)
    if int(finite.sum()) < 2:
        return output, "unavailable_fewer_than_two_finite_values"
    values = series.loc[finite].to_numpy(dtype=float)
    scaled = robust_scale_series(values)
    if not np.isfinite(scaled).all() or float(np.nanstd(scaled)) <= 1e-12:
        return output, "unavailable_degenerate_cohort_distribution"
    output.loc[finite] = scaled
    return output, "available_robust_cohort_standardization"


def _covariance_aware_score(enriched: pd.DataFrame, components: list[str]) -> tuple[pd.Series, dict[str, object]]:
    zscores: dict[str, pd.Series] = {}
    statuses: dict[str, str] = {}
    for component in components:
        zscore, status = _robust_component_zscore(_component_series(enriched, component))
        statuses[component] = status
        if zscore.notna().any():
            zscores[component] = zscore
    output = pd.Series(np.nan, index=enriched.index, dtype=float)
    if not zscores:
        return output, {
            "status": "unavailable",
            "reason": "No component had a non-degenerate cohort distribution.",
            "component_status": statuses,
            "effective_components": [],
            "high_correlation_pairs": [],
        }
    matrix = pd.DataFrame(zscores)
    correlation = matrix.corr(min_periods=2).fillna(0.0).copy()
    for component in correlation.columns:
        correlation.loc[component, component] = 1.0
    high_pairs: list[dict[str, object]] = []
    for left_index, left in enumerate(correlation.columns):
        for right in correlation.columns[left_index + 1 :]:
            coefficient = float(correlation.loc[left, right])
            if abs(coefficient) >= 0.85:
                high_pairs.append({"left": left, "right": right, "correlation": coefficient})
    for index, row in matrix.iterrows():
        available = row.dropna()
        if available.empty:
            continue
        subcorrelation = correlation.loc[available.index, available.index].to_numpy(dtype=float)
        denominator_squared = float(np.ones(len(available)) @ subcorrelation @ np.ones(len(available)))
        if denominator_squared <= 1e-12:
            continue
        output.loc[index] = float(available.sum() / math.sqrt(denominator_squared))
    return output, {
        "status": "available" if output.notna().all() else "partial",
        "reason": "Robust component z-scores combined using their observed correlation denominator.",
        "component_status": statuses,
        "effective_components": list(zscores),
        "high_correlation_pairs": high_pairs,
    }


def score_primitives(features: pd.DataFrame, *, progress: bool = False) -> pd.DataFrame:
    records = features.to_dict(orient="records")
    reporter = ProgressReporter("score-primitives", total=len(records)) if progress else None
    if reporter:
        reporter.start("auditing and standardizing primitive candidate screens")
    enriched_rows: list[dict] = []
    legacy_rows: list[dict] = []
    for index, row in enumerate(records, start=1):
        enriched = _enrich_row(row)
        enriched_rows.append(enriched)
        identity = {key: row.get(key) for key in ("region_id", "chrom", "start", "end", "block_id") if key in row}
        legacy_rows.append({**identity, **primitive_scores_for_row(enriched)})
        if reporter:
            reporter.update(index, message=str(identity.get("region_id", "")))
    enriched = pd.DataFrame(enriched_rows)
    out = pd.DataFrame([{key: row.get(key) for key in ("region_id", "chrom", "start", "end", "block_id") if key in row} for row in records])
    legacy = pd.DataFrame(legacy_rows)
    audits: dict[str, dict[str, object]] = {}
    known_score_matrix = pd.DataFrame(index=out.index)
    for column in PRIMITIVE_SCORE_COLUMNS:
        if column == "unexplained_dark_anomaly_candidate_score":
            continue
        score, audit = _covariance_aware_score(enriched, PRIMITIVE_SCORE_COMPONENTS[column])
        out[column] = score
        out[f"{column}_legacy_screening_composite"] = pd.to_numeric(legacy.get(column, np.nan), errors="coerce")
        known_score_matrix[column] = score
        audits[column] = audit

    groups = infer_block_groups(enriched)
    unexplained, unexplained_audit = cross_fitted_unexplained_outlierness(known_score_matrix, groups)
    unexplained_column = "unexplained_dark_anomaly_candidate_score"
    out[unexplained_column] = unexplained
    out[f"{unexplained_column}_legacy_screening_composite"] = math.nan
    audits[unexplained_column] = unexplained_audit

    for column in PRIMITIVE_SCORE_COLUMNS:
        robust, robust_status = _robust_component_zscore(pd.to_numeric(out[column], errors="coerce"))
        out[f"{column}_robust_zscore"] = robust
        out[f"{column}_empirical_p_value"] = math.nan
        out[f"{column}_empirical_p_value_status"] = "unavailable"
        out[f"{column}_empirical_p_value_reason"] = "No explicit null distribution is defined at primitive-score construction; observed genomic ranks are not null significance."
        out[f"{column}_component_features"] = ";".join(PRIMITIVE_SCORE_COMPONENTS.get(column, []))
        out[f"{column}_weighting_scheme"] = "covariance_aware_robust_cohort_standardization"
        out[f"{column}_calibration_status"] = "cohort_standardized_not_null_calibrated_not_probability"
        out[f"{column}_robust_zscore_status"] = robust_status
        out[f"{column}_correlation_audit"] = json.dumps(audits[column], sort_keys=True)
    if reporter:
        reporter.finish()
    return out


def write_primitive_scores(scores: pd.DataFrame, outdir: str | Path) -> Path:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "primitive_scores.parquet"
    scores.to_parquet(path, index=False)
    manifest_path = out / "primitive_score_manifest.json"
    manifest_path.write_text(json.dumps(primitive_score_manifest(), indent=2), encoding="utf-8")
    return path
