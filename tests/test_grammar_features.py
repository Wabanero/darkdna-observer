from darkdna.features.grammar import compute_grammar_features


def test_grammar_graph_and_periodicity_features():
    features = compute_grammar_features("AAAAACCCCC" * 40)
    assert features["kmer_adjacency_graph_nodes"] > 0
    assert features["kmer_adjacency_graph_edges"] > 0
    assert "phase_periodicity_around_10bp" in features
    assert "mutual_information_between_distant_sequence_positions" in features
