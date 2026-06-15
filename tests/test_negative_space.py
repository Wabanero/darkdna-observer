from darkdna.features.negative_space import compute_negative_space_features


def test_negative_space_depleted_kmer_features():
    features = compute_negative_space_features("ATATTAATTA" * 30)
    assert features["depleted_kmer_score"] > 0
    assert features["unexpected_silence_score"] > 0
    assert "negative_space_boundary_score" in features
