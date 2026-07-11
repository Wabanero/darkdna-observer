"""Primitive-specific screening composite composition.

These columns are compatibility screening views, not calibrated probabilities
and not direct mechanistic measurements.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.utils.progress import ProgressReporter
from darkdna.utils.stats import empirical_p_value, robust_scale_series
from .boundary_conditions import compute_boundary_condition_view
from .entropy_noise import compute_entropy_noise_view
from .physical_susceptibility import compute_physical_susceptibility_view
from .scale_fractal import compute_scale_fractal_features


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
        "fractal_score",
        "scale_persistence_score",
        "compression_anomaly_score",
    ],
    "constraint_grammar_region_candidate_score": [
        "grammar_entropy",
        "forbidden_word_depletion_enrichment",
        "motif_like_token_recurrence",
        "Markov_order_anomaly",
        "long_range_dependency_proxy",
    ],
    "non_B_DNA_physical_susceptibility_candidate_score": ["nonB_physical_susceptibility_score"],
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
        "non_B_DNA_aggregate_score",
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
        "DNA_bendability_proxy",
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
        "mean_of_other_screening_composites",
    ],
}

SCREENING_COMPOSITE_CAVEAT = (
    "Primitive score columns are equal-weight screening composites for ranking and review. "
    "They are not calibrated probabilities, not selected mechanistic measurements, and not "
    "evidence for the primitive hypothesis without residualization, severe null models, "
    "held-out calibration, and validation."
)

REQUIRED_COMPOSITE_VALIDATION = [
    "justify_or_learn_weights",
    "test_stability_across_organisms",
    "test_robustness_to_sequence_transformations",
    "audit_double_counting_between_correlated_features",
    "calibrate_probabilistically",
    "reproduce_on_held_out_chromosomes_or_blocks",
]


def primitive_score_manifest() -> dict:
    return {
        "score_status": "uncalibrated_equal_weight_screening_composite",
        "caveat": SCREENING_COMPOSITE_CAVEAT,
        "interpretation_order": [
            "measured_feature_profile",
            "statistical_anomaly_after_controls_and_nulls",
            "post_hoc_mechanistic_hypothesis",
        ],
        "required_validation_before_mechanistic_use": REQUIRED_COMPOSITE_VALIDATION,
        "components": PRIMITIVE_SCORE_COMPONENTS,
    }


def value(row: dict, key: str, default: float = 0.0) -> float:
    try:
        v = row.get(key, default)
        if v is None or pd.isna(v):
            return default
        return float(v)
    except Exception:
        return default


def primitive_scores_for_row(row: dict) -> dict[str, float]:
    enriched = dict(row)
    if "fractal_score" not in enriched and "sequence" in enriched:
        enriched.update(compute_scale_fractal_features(str(enriched["sequence"])))
    enriched.update(compute_physical_susceptibility_view(enriched))
    enriched.update(compute_entropy_noise_view(enriched))
    enriched.update(compute_boundary_condition_view(enriched))

    fractal = np.mean([value(enriched, "fractal_score"), value(enriched, "scale_persistence_score"), value(enriched, "compression_anomaly_score")])
    grammar = np.mean(
        [
            value(enriched, "grammar_entropy"),
            value(enriched, "forbidden_word_depletion_enrichment"),
            value(enriched, "motif_like_token_recurrence"),
            value(enriched, "Markov_order_anomaly"),
            value(enriched, "long_range_dependency_proxy"),
        ]
    )
    physical_susceptibility = value(enriched, "nonB_physical_susceptibility_score")
    replication = value(enriched, "fork_texture_score")
    decoherence = value(enriched, "decoherence_boundary_candidate_score")
    resonant = np.mean(
        [
            value(enriched, "phase_periodicity_around_10bp"),
            value(enriched, "nucleosome_scale_periodicity_around_147bp"),
            value(enriched, "spacing_periodicity_fourier_power"),
        ]
    )
    hysteresis = np.mean(
        [
            abs(value(enriched, "left_right_GC_asymmetry")),
            value(enriched, "nested_repeat_architecture_score"),
            value(enriched, "non_B_DNA_aggregate_score"),
            value(enriched, "orientation_bias_of_recurrent_kmers"),
        ]
    )
    negative = np.mean(
        [
            value(enriched, "depleted_kmer_score"),
            value(enriched, "unexpected_silence_score"),
            value(enriched, "local_feature_void_score"),
        ]
    )
    boundary = value(enriched, "boundary_condition_candidate_score")
    te = np.mean([value(enriched, "TE_overlap_fraction"), value(enriched, "TE_family_mosaic_score"), value(enriched, "TE_boundary_score")])
    # Prompt 1 can only derive static sequence-compatible candidate proxies for
    # dynamic labels. It does not infer future states, thresholds, memory, or
    # trajectories without Prompt 2 dynamic data.
    possibility = np.mean([boundary, negative, value(enriched, "forbidden_word_depletion_enrichment")])
    criticality = np.mean([boundary, value(enriched, "entropy_boundary_score"), value(enriched, "compression_boundary_score")])
    chromatin_motion = np.mean(
        [
            value(enriched, "spacing_periodicity_autocorrelation"),
            value(enriched, "DNA_bendability_proxy"),
            value(enriched, "left_right_entropy_asymmetry"),
        ]
    )
    scores = {
        "fractal_scaffold_candidate_score": float(fractal),
        "constraint_grammar_region_candidate_score": float(grammar),
        "non_B_DNA_physical_susceptibility_candidate_score": float(physical_susceptibility),
        "replication_instability_candidate_score": float(replication),
        "decoherence_boundary_candidate_score": float(decoherence),
        "resonant_pulse_decoder_candidate_score": float(resonant),
        "hysteresis_candidate_score": float(hysteresis),
        "possibility_gate_candidate_score": float(possibility),
        "criticality_tuner_candidate_score": float(criticality),
        "chromatin_motion_oscillator_candidate_score": float(chromatin_motion),
        "negative_space_element_candidate_score": float(negative),
        "sequence_regime_boundary_candidate_score": float(boundary),
        "TE_grammar_node_candidate_score": float(te),
    }
    scores["unexplained_dark_anomaly_candidate_score"] = float(np.nanmean(list(scores.values())))
    return scores


def score_primitives(features: pd.DataFrame, *, progress: bool = False) -> pd.DataFrame:
    rows = []
    records = features.to_dict(orient="records")
    reporter = ProgressReporter("score-primitives", total=len(records)) if progress else None
    if reporter:
        reporter.start("scoring primitive candidates")
    for idx, row in enumerate(records, start=1):
        identity = {k: row.get(k) for k in ["region_id", "chrom", "start", "end"] if k in row}
        scores = primitive_scores_for_row(row)
        rows.append({**identity, **scores})
        if reporter:
            reporter.update(idx, message=str(identity.get("region_id", "")))
    if reporter:
        reporter.finish()
    out = pd.DataFrame(rows)
    for col in PRIMITIVE_SCORE_COLUMNS:
        if col not in out.columns:
            out[col] = 0.0
        out[f"{col}_robust_zscore"] = robust_scale_series(out[col].fillna(0.0))
        values = out[col].fillna(0.0).to_numpy()
        out[f"{col}_empirical_p_value"] = [empirical_p_value(v, values, higher=True) for v in values]
        out[f"{col}_component_features"] = ";".join(PRIMITIVE_SCORE_COMPONENTS.get(col, []))
        out[f"{col}_weighting_scheme"] = "equal_weight_mean_screening_composite"
        out[f"{col}_calibration_status"] = "uncalibrated_requires_held_out_and_null_panel_calibration"
    return out


def write_primitive_scores(scores: pd.DataFrame, outdir: str | Path) -> Path:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "primitive_scores.parquet"
    scores.to_parquet(path, index=False)
    manifest_path = out / "primitive_score_manifest.json"
    manifest_path.write_text(json.dumps(primitive_score_manifest(), indent=2), encoding="utf-8")
    return path
