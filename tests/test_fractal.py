from darkdna.views.scale_fractal import compute_scale_fractal_features


def test_fractal_scale_features_have_fallbacks():
    features = compute_scale_fractal_features("ACGTACGTGGCC" * 40)
    assert "DFA_like_exponent_estimator" in features
    assert "chaos_game_k3_occupancy" in features
    assert features["fractal_score"] >= 0
