from darkdna.features.nonb_dna import compute_nonb_dna_features


def test_nonb_dna_proxy_scores_are_conservative_numeric_outputs():
    features = compute_nonb_dna_features("GGGAGGGCGGGTT" * 10)
    assert features["G4_susceptibility_proxy"] > 0
    assert features["R_loop_forming_potential"] >= 0
    assert features["non_B_DNA_aggregate_score"] >= 0
