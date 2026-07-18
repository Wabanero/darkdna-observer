import math

from darkdna.views.scale_fractal import compute_scale_fractal_features


def test_multiscale_texture_reports_diagnostics_and_deprecated_alias():
    sequence = ("ACGTGGCCATAT" * 180)[:2000]
    result = compute_scale_fractal_features(sequence, n_surrogates=6, seed=7)

    assert result["multiscale_method"] == "multiscale_texture_v2"
    assert result["multiscale_texture_status"] in {"available", "unavailable"}
    assert result["fractal_score_status"] == "deprecated_alias"
    assert result["fractal_score"] == result["multiscale_texture_screening_score"] or (
        math.isnan(result["fractal_score"]) and math.isnan(result["multiscale_texture_screening_score"])
    )
    assert result["DFA_valid_scale_count"] >= 0
    assert result["DFA_window_shift_status"] in {"available", "unavailable"}
    assert "strong_weak_H_bond_numeric_walk" in result["multiscale_mapping_diagnostics"]


def test_short_sequence_is_explicitly_unavailable_not_zero_placeholder():
    result = compute_scale_fractal_features("ACGTACGT", n_surrogates=6)

    assert result["multiscale_texture_status"] == "unavailable"
    assert math.isnan(result["multiscale_texture_screening_score"])
    assert math.isnan(result["scale_persistence_score"])
    assert result["renormalization_profile_status"].startswith("unavailable")

