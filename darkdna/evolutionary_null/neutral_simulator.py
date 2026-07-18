"""Reference-conditioned neutral simulator and empirical feature calibration."""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.evolutionary_null.context_dependent_mutation import fit_transition_counts, generate_markov_sequence
from darkdna.evolutionary_null.cpg_deamination import apply_cpg_deamination
from darkdna.evolutionary_null.gc_biased_gene_conversion import apply_gc_bias
from darkdna.evolutionary_null.indel_model import apply_length_preserving_indel_turnover
from darkdna.evolutionary_null.microsatellite_slippage import apply_slippage_proxy
from darkdna.evolutionary_null.mutation_spectrum import MutationSpectrum, mutation_spectrum_from_table
from darkdna.evolutionary_null.repeat_birth_death import apply_repeat_turnover
from darkdna.features.repeats import simple_repeat_fraction
from darkdna.features.sequence import cpg_density, gc_content, lempel_ziv_complexity
from darkdna.utils.stats import empirical_p_value, shannon_entropy


@dataclass(frozen=True)
class EvolutionaryNullModel:
    base_frequencies: dict[str, float]
    transition_counts: dict[str, dict[str, int]]
    mutation_spectrum: MutationSpectrum
    cpg_event_probability: float = 0.01
    gc_bgc_strength: float = 0.0
    indel_turnover_rate: float = 0.005
    slippage_rate: float = 0.005
    repeat_turnover_rate: float = 0.002
    calibration_status: str = "reference_conditioned_generic_processes"
    limitation: str = "Not an ancestral reconstruction and not organism-specific without supplied mutation data."


def fit_evolutionary_null(
    sequences: list[str],
    mutation_table: pd.DataFrame | None = None,
) -> EvolutionaryNullModel:
    counts = Counter(base for sequence in sequences for base in sequence.upper() if base in "ACGT")
    total = sum(counts.values()) or 1
    spectrum = mutation_spectrum_from_table(mutation_table)
    status = "organism_conditioned_partial" if spectrum.calibration_status == "organism_conditioned_partial" else "reference_conditioned_generic_processes"
    return EvolutionaryNullModel(
        base_frequencies={base: counts[base] / total for base in "ACGT"},
        transition_counts=fit_transition_counts(sequences),
        mutation_spectrum=spectrum,
        calibration_status=status,
    )


def simulate_neutral_sequence(model: EvolutionaryNullModel, length: int, *, seed: int) -> str:
    sequence = generate_markov_sequence(length, model.transition_counts, model.base_frequencies, seed=seed)
    sequence = apply_cpg_deamination(sequence, model.cpg_event_probability, seed=seed + 1)
    sequence = apply_gc_bias(sequence, model.gc_bgc_strength, seed=seed + 2)
    sequence = apply_length_preserving_indel_turnover(sequence, model.indel_turnover_rate, seed=seed + 3)
    sequence = apply_slippage_proxy(sequence, model.slippage_rate, seed=seed + 4)
    sequence = apply_repeat_turnover(sequence, model.repeat_turnover_rate, seed=seed + 5)
    return sequence


def _features(sequence: str) -> dict[str, float]:
    clean = "".join(base for base in sequence.upper() if base in "ACGT")
    return {
        "gc_content": gc_content(clean),
        "CpG_density": cpg_density(clean),
        "Shannon_entropy": shannon_entropy(clean),
        "Lempel_Ziv_complexity": lempel_ziv_complexity(clean),
        "simple_repeat_fraction": simple_repeat_fraction(clean),
    }


def build_evolutionary_null_scores(
    sequences: dict[str, str],
    *,
    n_surrogates: int = 25,
    seed: int = 13,
    mutation_table: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, EvolutionaryNullModel]:
    model = fit_evolutionary_null(list(sequences.values()), mutation_table)
    rows: list[dict[str, object]] = []
    for region_index, (region_id, sequence) in enumerate(sequences.items()):
        observed = _features(sequence)
        null_features = [
            _features(simulate_neutral_sequence(model, len(sequence), seed=seed + region_index * 100_003 + replicate))
            for replicate in range(max(1, n_surrogates))
        ]
        for feature, value in observed.items():
            null_values = np.asarray([record[feature] for record in null_features], dtype=float)
            mean = float(np.mean(null_values))
            std = float(np.std(null_values, ddof=1)) if null_values.size > 1 else math.nan
            zscore = float((value - mean) / std) if np.isfinite(std) and std > 0 else math.nan
            p_value = empirical_p_value(value, null_values, higher=True)
            rows.append(
                {
                    "region_id": region_id,
                    "feature": feature,
                    "observed_value": float(value),
                    "evolutionary_null_mean": mean,
                    "evolutionary_null_std": std,
                    "evolutionary_process_null_zscore": zscore,
                    "evolutionary_null_empirical_p": float(p_value),
                    "null_sample_size": int(null_values.size),
                    "calibration_status": model.calibration_status,
                    "limitation": model.limitation,
                }
            )
    return pd.DataFrame(rows), model


def write_evolutionary_null_outputs(
    scores: pd.DataFrame,
    model: EvolutionaryNullModel,
    outdir: str | Path,
) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    scores_path = out / "evolutionary_null_scores.parquet"
    manifest_path = out / "evolutionary_model_manifest.json"
    scores.to_parquet(scores_path, index=False)
    payload = asdict(model)
    payload["mutation_spectrum"] = asdict(model.mutation_spectrum)
    payload["interpretation"] = "Neutral-process baseline only; excess does not prove selected function."
    manifest_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return {"scores": scores_path, "manifest": manifest_path}
