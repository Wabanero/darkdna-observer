import numpy as np

from darkdna.evolutionary_null import (
    build_evolutionary_null_scores,
    fit_evolutionary_null,
    simulate_neutral_sequence,
)


def test_evolutionary_null_is_deterministic_length_preserving_and_explicitly_calibrated():
    sequences = {"r1": "ACGTCG" * 20, "r2": "ATATGC" * 20}
    model = fit_evolutionary_null(list(sequences.values()))
    first = simulate_neutral_sequence(model, 121, seed=17)
    second = simulate_neutral_sequence(model, 121, seed=17)
    assert first == second
    assert len(first) == 121
    assert model.calibration_status == "reference_conditioned_generic_processes"
    scores, fitted = build_evolutionary_null_scores(sequences, n_surrogates=7, seed=19)
    assert {
        "evolutionary_process_null_zscore",
        "evolutionary_null_empirical_p",
        "calibration_status",
        "limitation",
    }.issubset(scores.columns)
    assert scores["evolutionary_null_empirical_p"].between(0, 1).all()
    assert fitted.limitation
