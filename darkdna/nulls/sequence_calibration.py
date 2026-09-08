"""Candidate-level sequence-transform and evolutionary-process nulls.

These families test whether a primitive screen is unusual for the exact
interval sequence. Genomic matching is a covariate control and is not a
substitute for this counterfactual.
"""

from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np
import pandas as pd

from darkdna.evolutionary_null import fit_evolutionary_null, simulate_neutral_sequence
from darkdna.features.sequence import compute_all_sequence_features
from darkdna.nulls.sequence import transform_sequence
from darkdna.utils.stats import empirical_p_value
from darkdna.views.primitive_scores import PRIMITIVE_SCORE_COLUMNS, primitive_scores_for_row


STOCHASTIC_SEQUENCE_NULLS: dict[str, str] = {
    "mononucleotide_preserving": "mononucleotide_preserving",
    "dinucleotide_preserving_shuffle": "dinucleotide_preserving_shuffle",
    "kmer_preserving_shuffle": "kmer_preserving_shuffle",
    "markov_chain_surrogate": "markov_chain_surrogate",
    "synthetic_equal_composition": "synthetic_equal_composition",
}
STOCHASTIC_SEED_OFFSET = {model_id: index * 97 for index, model_id in enumerate(STOCHASTIC_SEQUENCE_NULLS)}

PAIRED_ORIENTATION_NULLS: dict[str, str] = {
    "reversed_sequence": "reversed_sequence",
    "reverse_complement": "reverse_complement",
}

GENERATIVE_SEQUENCE_NULLS: dict[str, str] = {
    "evolutionary_process_generated": "evolutionary_process_generated",
}

SEQUENCE_NULL_MODEL_IDS: tuple[str, ...] = tuple(
    list(STOCHASTIC_SEQUENCE_NULLS) + list(PAIRED_ORIENTATION_NULLS) + list(GENERATIVE_SEQUENCE_NULLS)
)

UNEXPLAINED = "unexplained_dark_anomaly_candidate_score"


def sequences_from_feature_table(features: pd.DataFrame) -> dict[str, str]:
    if features is None or features.empty or "region_id" not in features.columns:
        return {}
    if "sequence" not in features.columns:
        return {}
    sequences: dict[str, str] = {}
    for row in features.itertuples(index=False):
        sequence = str(getattr(row, "sequence", "") or "")
        if sequence:
            sequences[str(row.region_id)] = sequence
    return sequences


def _score_sequence(sequence: str) -> dict[str, float]:
    return primitive_scores_for_row(compute_all_sequence_features(sequence))


def _detail_row(
    *,
    region_id: str,
    primitive: str,
    model_id: str,
    observed: float,
    values: np.ndarray,
    calibration_block_id: str,
    status: str,
    reason: str,
    matched_features_used: str,
) -> dict[str, object]:
    finite = values[np.isfinite(values)] if values.size else np.array([], dtype=float)
    mean = float(np.mean(finite)) if finite.size else math.nan
    std = float(np.std(finite, ddof=1)) if finite.size > 1 else math.nan
    zscore = float((observed - mean) / std) if np.isfinite(observed) and np.isfinite(std) and std > 0 else math.nan
    p_value = empirical_p_value(float(observed), finite, higher=True) if np.isfinite(observed) and finite.size else math.nan
    return {
        "region_id": region_id,
        "primitive": primitive,
        "null_model_id": model_id,
        "primitive_score": float(observed) if np.isfinite(observed) else math.nan,
        "null_mean": mean,
        "null_std": std,
        "null_zscore": zscore,
        "empirical_p_value": float(p_value) if np.isfinite(p_value) else math.nan,
        "null_sample_size": int(finite.size),
        "independent_block_count": int(finite.size),
        "calibration_block_id": calibration_block_id,
        "matched_features_used": matched_features_used,
        "null_status": status,
        "null_reason": reason,
        "null_execution_mode": "sequence_transform" if model_id != "evolutionary_process_generated" else "generative",
    }


def build_sequence_transform_null_details(
    scores: pd.DataFrame,
    features: pd.DataFrame,
    sequences: Mapping[str, str],
    primitives: list[str],
    *,
    n_surrogates: int = 8,
    seed: int = 13,
    kmer_size: int = 3,
    include_evolutionary: bool = True,
) -> pd.DataFrame:
    """Score named sequence nulls per region using transforms of the focal sequence."""

    if not sequences:
        return pd.DataFrame()
    native_lookup = features.drop_duplicates("region_id").set_index("region_id") if "region_id" in features.columns else pd.DataFrame()
    evolutionary_model = fit_evolutionary_null(list(sequences.values())) if include_evolutionary else None
    n_surrogates = max(2, int(n_surrogates))
    rows: list[dict[str, object]] = []
    region_ids = [str(region_id) for region_id in scores["region_id"].drop_duplicates()] if "region_id" in scores.columns else list(sequences)
    for region_index, region_id in enumerate(region_ids):
        sequence = sequences.get(region_id, "")
        calibration_block = str(native_lookup.loc[region_id, "calibration_block_id"]) if not native_lookup.empty and region_id in native_lookup.index and "calibration_block_id" in native_lookup.columns else f"sequence:{region_id}"
        if not sequence or not any(base in sequence.upper() for base in "ACGT"):
            for model_id in SEQUENCE_NULL_MODEL_IDS:
                for primitive in primitives:
                    rows.append(
                        _detail_row(
                            region_id=region_id,
                            primitive=primitive,
                            model_id=model_id,
                            observed=math.nan,
                            values=np.array([], dtype=float),
                            calibration_block_id=calibration_block,
                            status="unavailable_missing_sequence",
                            reason="No ACGT sequence was available for this interval.",
                            matched_features_used="sequence",
                        )
                    )
            continue
        # Score the native interval with the same row-local screens as its
        # transforms. Cohort-standardized Mode A scores live on a different
        # scale and must not enter this counterfactual.
        native_scores = _score_sequence(sequence)

        def observed_value(primitive: str) -> float:
            if primitive in native_scores and np.isfinite(native_scores[primitive]):
                return float(native_scores[primitive])
            return math.nan

        for model_id, method in STOCHASTIC_SEQUENCE_NULLS.items():
            surrogate_scores: dict[str, list[float]] = {primitive: [] for primitive in primitives}
            for replicate in range(n_surrogates):
                transformed = transform_sequence(
                    sequence,
                    method,
                    seed=seed + region_index * 100_003 + replicate * 1_013 + STOCHASTIC_SEED_OFFSET[model_id],
                    k=kmer_size,
                )
                scored = _score_sequence(transformed)
                for primitive in primitives:
                    surrogate_scores[primitive].append(float(scored.get(primitive, math.nan)))
            for primitive in primitives:
                if primitive == UNEXPLAINED:
                    rows.append(
                        _detail_row(
                            region_id=region_id,
                            primitive=primitive,
                            model_id=model_id,
                            observed=math.nan,
                            values=np.array([], dtype=float),
                            calibration_block_id=calibration_block,
                            status="unavailable_not_a_sequence_null",
                            reason="Held-out multivariate outlierness is not calibrated by sequence transformation.",
                            matched_features_used="sequence",
                        )
                    )
                    continue
                values = np.asarray(surrogate_scores[primitive], dtype=float)
                finite = values[np.isfinite(values)]
                observed = observed_value(primitive)
                if finite.size < 2:
                    status, reason = "unavailable_insufficient_controls", "Fewer than two finite sequence-transform scores were produced."
                elif not np.isfinite(observed):
                    status, reason = "unavailable_observed_score_na", "The native sequence screen is unavailable, so the sequence null cannot be evaluated."
                elif float(np.std(finite, ddof=1)) <= 0:
                    status, reason = "partial_zero_null_variance", "The sequence-transform null distribution had zero or undefined variance."
                else:
                    status, reason = "available_sequence_calibrated", "Empirical calibration used transformations of the focal interval sequence."
                rows.append(
                    _detail_row(
                        region_id=region_id,
                        primitive=primitive,
                        model_id=model_id,
                        observed=observed,
                        values=values,
                        calibration_block_id=calibration_block,
                        status=status,
                        reason=reason,
                        matched_features_used="sequence",
                    )
                )

        for model_id, method in PAIRED_ORIENTATION_NULLS.items():
            transformed = transform_sequence(sequence, method, seed=seed)
            scored = _score_sequence(transformed)
            for primitive in primitives:
                if primitive == UNEXPLAINED:
                    rows.append(
                        _detail_row(
                            region_id=region_id,
                            primitive=primitive,
                            model_id=model_id,
                            observed=math.nan,
                            values=np.array([], dtype=float),
                            calibration_block_id=calibration_block,
                            status="unavailable_not_a_sequence_null",
                            reason="Held-out multivariate outlierness is not calibrated by sequence transformation.",
                            matched_features_used="sequence",
                        )
                    )
                    continue
                observed = observed_value(primitive)
                transformed_value = float(scored.get(primitive, math.nan))
                values = np.asarray([transformed_value], dtype=float)
                rows.append(
                    _detail_row(
                        region_id=region_id,
                        primitive=primitive,
                        model_id=model_id,
                        observed=observed,
                        values=values,
                        calibration_block_id=calibration_block,
                        status="partial_paired_orientation_transform",
                        reason="Reverse/reverse-complement is a paired orientation control, not a resampled null distribution.",
                        matched_features_used="sequence",
                    )
                )

        if evolutionary_model is not None:
            surrogate_scores = {primitive: [] for primitive in primitives}
            for replicate in range(n_surrogates):
                generated = simulate_neutral_sequence(
                    evolutionary_model,
                    len(sequence),
                    seed=seed + region_index * 17_903 + replicate,
                )
                scored = _score_sequence(generated)
                for primitive in primitives:
                    surrogate_scores[primitive].append(float(scored.get(primitive, math.nan)))
            for primitive in primitives:
                if primitive == UNEXPLAINED:
                    rows.append(
                        _detail_row(
                            region_id=region_id,
                            primitive=primitive,
                            model_id="evolutionary_process_generated",
                            observed=math.nan,
                            values=np.array([], dtype=float),
                            calibration_block_id=calibration_block,
                            status="unavailable_not_a_sequence_null",
                            reason="Held-out multivariate outlierness is not calibrated by sequence transformation.",
                            matched_features_used="sequence,evolutionary_model",
                        )
                    )
                    continue
                values = np.asarray(surrogate_scores[primitive], dtype=float)
                finite = values[np.isfinite(values)]
                observed = observed_value(primitive)
                if finite.size < 2:
                    status, reason = "unavailable_insufficient_controls", "Fewer than two evolutionary-process surrogates were produced."
                elif not np.isfinite(observed):
                    status, reason = "unavailable_observed_score_na", "The native sequence screen is unavailable, so the evolutionary-process null cannot be evaluated."
                elif float(np.std(finite, ddof=1)) <= 0:
                    status, reason = "partial_zero_null_variance", "The evolutionary-process null distribution had zero or undefined variance."
                else:
                    status, reason = (
                        "available_sequence_calibrated",
                        "Reference-conditioned evolutionary-process surrogate of the focal interval; not an ancestral reconstruction.",
                    )
                rows.append(
                    _detail_row(
                        region_id=region_id,
                        primitive=primitive,
                        model_id="evolutionary_process_generated",
                        observed=observed,
                        values=values,
                        calibration_block_id=calibration_block,
                        status=status,
                        reason=reason,
                        matched_features_used="sequence,evolutionary_model",
                    )
                )
    return pd.DataFrame(rows)
