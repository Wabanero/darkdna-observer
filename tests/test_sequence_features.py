from darkdna.features.sequence import compute_all_sequence_features, compute_sequence_features


def test_gc_entropy_lz_and_compression_features():
    features = compute_sequence_features("ACGTACGTNN")
    assert features["length"] == 10
    assert round(features["gc_content"], 3) == 0.5
    assert features["Shannon_entropy"] > 1.9
    assert features["Lempel_Ziv_complexity"] > 0
    assert features["gzip_compression_ratio"] > 0
    assert "k4_mer_entropy" in features


def test_all_sequence_features_include_view_inputs():
    features = compute_all_sequence_features("GGGAGGGCGGGTTACGTACGTACGT")
    assert "G4_susceptibility_proxy" in features
    assert "grammar_entropy" in features
    assert "fractal_score" in features
