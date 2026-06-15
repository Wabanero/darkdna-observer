from darkdna.features.boundaries import compute_boundary_features


def test_boundary_features_detect_left_right_shift():
    features = compute_boundary_features("A" * 100 + "GCGCGCGC" * 12)
    assert features["GC_boundary_score"] > 0.5
    assert features["entropy_boundary_score"] >= 0
    assert features["segmentation_breakpoint_score"] >= 0
