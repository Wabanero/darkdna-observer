"""Counterfactual sequence/quantity transformations for Mode B screening."""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np

from darkdna.features.repeats import simple_repeat_fraction
from darkdna.features.sequence import gc_content, lempel_ziv_complexity
from darkdna.nulls.sequence import transform_sequence
from darkdna.utils.stats import shannon_entropy


Predictor = Callable[[str, int, float, float | None], float]


def default_architecture_screen_predictor(
    sequence: str,
    length: int,
    copy_number: float,
    spacing: float | None,
) -> float:
    """Transparent screening proxy used only when no user model is supplied."""

    clean = "".join(base for base in sequence.upper() if base in "ACGT")
    if not clean or length <= 0 or copy_number <= 0:
        return 0.0
    entropy_structure = 1.0 - min(1.0, shannon_entropy(clean) / 2.0)
    repeat_structure = simple_repeat_fraction(clean)
    gc_extremeness = abs(gc_content(clean) - 0.5) * 2.0
    lz = min(1.0, lempel_ziv_complexity(clean))
    sequence_component = float(np.mean([entropy_structure, repeat_structure, gc_extremeness, 1.0 - lz]))
    quantity_component = 0.10 * math.log2(max(2, length)) + 0.20 * math.log2(1.0 + copy_number)
    spacing_component = 0.0 if spacing is None or not np.isfinite(spacing) else 0.02 * math.log2(1.0 + max(0.0, spacing))
    return sequence_component + quantity_component + spacing_component


def _resize_composition_matched(sequence: str, target_length: int, seed: int) -> str:
    if target_length <= 0:
        return ""
    if not sequence:
        return "N" * target_length
    repeated = (sequence * (target_length // len(sequence) + 1))[:target_length]
    return transform_sequence(repeated, "mononucleotide_preserving", seed=seed)


def _resize_native_pattern(sequence: str, target_length: int) -> str:
    """Resize while sampling across the native interval instead of one biased prefix."""

    if target_length <= 0 or not sequence:
        return ""
    if target_length >= len(sequence):
        copies, remainder = divmod(target_length, len(sequence))
        return sequence * copies + _resize_native_pattern(sequence, remainder)
    block_count = min(10, target_length)
    pieces = []
    assigned = 0
    for block in range(block_count):
        source_start = int(round(block * len(sequence) / block_count))
        source_end = int(round((block + 1) * len(sequence) / block_count))
        take = int(round((block + 1) * target_length / block_count)) - assigned
        assigned += take
        source = sequence[source_start:source_end]
        pieces.append((source * (take // max(1, len(source)) + 1))[:take])
    return "".join(pieces)[:target_length]


def _closest_replacement(sequence: str, pool: list[str], *, by_repeat: bool) -> str | None:
    candidates = [candidate for candidate in pool if candidate and candidate != sequence]
    if not candidates:
        return None
    if by_repeat:
        target = simple_repeat_fraction(sequence)
        return min(candidates, key=lambda candidate: abs(simple_repeat_fraction(candidate) - target))
    target = gc_content(sequence)
    return min(candidates, key=lambda candidate: abs(gc_content(candidate) - target))


def generate_sequence_indifference_controls(
    sequence: str,
    *,
    replacement_pool: list[str] | None = None,
    seed: int = 13,
    k: int = 3,
) -> list[dict[str, object]]:
    sequence = sequence.upper()
    pool = replacement_pool or []
    gc_replacement = _closest_replacement(sequence, pool, by_repeat=False)
    repeat_replacement = _closest_replacement(sequence, pool, by_repeat=True)
    rows: list[dict[str, object]] = [
        {"transformation": "native", "sequence": sequence, "length_factor": 1.0, "copy_number": 1.0, "control_family": "native"},
        {"transformation": "reverse", "sequence": transform_sequence(sequence, "reversed_sequence", seed=seed), "length_factor": 1.0, "copy_number": 1.0, "control_family": "orientation"},
        {"transformation": "reverse_complement", "sequence": transform_sequence(sequence, "reverse_complement", seed=seed), "length_factor": 1.0, "copy_number": 1.0, "control_family": "orientation"},
        {"transformation": "mononucleotide_shuffle", "sequence": transform_sequence(sequence, "mononucleotide_preserving", seed=seed), "length_factor": 1.0, "copy_number": 1.0, "control_family": "identity"},
        {"transformation": "dinucleotide_shuffle", "sequence": transform_sequence(sequence, "dinucleotide_preserving_shuffle", seed=seed), "length_factor": 1.0, "copy_number": 1.0, "control_family": "identity"},
        {"transformation": "kmer_preserving_shuffle", "sequence": transform_sequence(sequence, "kmer_preserving_shuffle", seed=seed, k=k), "length_factor": 1.0, "copy_number": 1.0, "control_family": "identity"},
        {"transformation": "GC_matched_unrelated_equal_length", "sequence": _resize_composition_matched(gc_replacement or sequence, len(sequence), seed + 1), "length_factor": 1.0, "copy_number": 1.0, "control_family": "identity"},
        {"transformation": "repeat_matched_unrelated_equal_length", "sequence": _resize_composition_matched(repeat_replacement or sequence, len(sequence), seed + 2), "length_factor": 1.0, "copy_number": 1.0, "control_family": "identity"},
        {"transformation": "equal_length_replacement", "sequence": _resize_composition_matched(gc_replacement or sequence, len(sequence), seed + 3), "length_factor": 1.0, "copy_number": 1.0, "control_family": "identity"},
    ]
    for factor in (0.25, 0.5, 0.75, 1.25, 1.5, 2.0):
        rows.append(
            {
                "transformation": f"length_titration_{factor:g}x",
                "sequence": _resize_native_pattern(sequence, max(1, int(round(len(sequence) * factor)))),
                "length_factor": factor,
                "copy_number": 1.0,
                "control_family": "length",
            }
        )
    for copies in (0.0, 1.0, 2.0, 4.0, 8.0):
        rows.append(
            {
                "transformation": f"copy_number_titration_{copies:g}",
                "sequence": sequence,
                "length_factor": 1.0,
                "copy_number": copies,
                "control_family": "copy_number",
            }
        )
    rows.extend(
        [
            {"transformation": "length_2x_copy_1", "sequence": _resize_native_pattern(sequence, len(sequence) * 2), "length_factor": 2.0, "copy_number": 1.0, "control_family": "factorial"},
            {"transformation": "length_1x_copy_2", "sequence": sequence, "length_factor": 1.0, "copy_number": 2.0, "control_family": "factorial"},
            {"transformation": "length_2x_copy_2", "sequence": _resize_native_pattern(sequence, len(sequence) * 2), "length_factor": 2.0, "copy_number": 2.0, "control_family": "factorial"},
        ]
    )
    return rows


def evaluate_sequence_indifference(
    region_id: str,
    sequence: str,
    *,
    predictor: Predictor | None = None,
    replacement_pool: list[str] | None = None,
    spacing: float | None = None,
    seed: int = 13,
    k: int = 3,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    model = predictor or default_architecture_screen_predictor
    predictor_id = "user_supplied_predictor" if predictor is not None else "intrinsic_plus_quantity_screen_v1"
    controls = generate_sequence_indifference_controls(sequence, replacement_pool=replacement_pool, seed=seed, k=k)
    for record in controls:
        transformed = str(record.pop("sequence"))
        length = len(transformed)
        score = float(model(transformed, length, float(record["copy_number"]), spacing))
        record.update(
            {
                "region_id": region_id,
                "transformed_length": length,
                "prediction": score,
                "predictor_id": predictor_id,
                "evidence_scope": "model_based_perturbation_not_biological_causality",
            }
        )
    native = next(record for record in controls if record["transformation"] == "native")
    denominator = max(1.0, abs(float(native["prediction"])))

    def sensitivity(family: str) -> float:
        values = [float(record["prediction"]) for record in controls if record["control_family"] == family]
        return (max(values) - min(values)) / denominator if len(values) > 1 else math.nan

    identity_values = [abs(float(record["prediction"]) - float(native["prediction"])) / denominator for record in controls if record["control_family"] == "identity"]
    identity = float(np.median(identity_values)) if identity_values else math.nan
    length_sensitivity = sensitivity("length")
    copy_sensitivity = sensitivity("copy_number")
    orientation = sensitivity("orientation")
    composition_records = [record for record in controls if record["transformation"] in {"mononucleotide_shuffle", "GC_matched_unrelated_equal_length"}]
    composition = float(np.mean([abs(float(record["prediction"]) - float(native["prediction"])) / denominator for record in composition_records]))
    repeat_record = next(record for record in controls if record["transformation"] == "repeat_matched_unrelated_equal_length")
    repeat_sensitivity = abs(float(repeat_record["prediction"]) - float(native["prediction"])) / denominator
    factorial = {record["transformation"]: float(record["prediction"]) for record in controls if record["control_family"] == "factorial"}
    interaction = abs(
        factorial["length_2x_copy_2"]
        - factorial["length_2x_copy_1"]
        - factorial["length_1x_copy_2"]
        + float(native["prediction"])
    ) / denominator
    quantity = max(length_sensitivity, copy_sensitivity)
    indifference = max(0.0, quantity) * max(0.0, 1.0 - min(1.0, identity))
    summary = {
        "region_id": region_id,
        "predictor_id": predictor_id,
        "sequence_identity_sensitivity": identity,
        "length_sensitivity": length_sensitivity,
        "copy_number_sensitivity": copy_sensitivity,
        "composition_sensitivity": composition,
        "repeat_fraction_sensitivity": repeat_sensitivity,
        "orientation_sensitivity": orientation,
        "sequence_indifference_score": indifference,
        "sequence_quantity_interaction_score": interaction,
        "sequence_indifference_status": "available_model_based_screen",
        "sequence_indifference_reason": "Equal-length identity controls and length/copy titrations were evaluated with the named predictor.",
        "sequence_indifference_caveat": "Model-based perturbation evidence is not biological causality or selected-function evidence.",
    }
    return summary, controls
