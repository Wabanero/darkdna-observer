from darkdna.features.asymmetry import compute_asymmetry_features


def test_asymmetry_and_skew_features():
    features = compute_asymmetry_features("GGGGCCCCAAAATTTT" * 5)
    assert "GC_skew" in features
    assert "left_right_GC_asymmetry" in features
    assert "kmer_strand_asymmetry" in features
