from darkdna.features.te_grammar import compute_te_grammar_for_region
from darkdna.io.gff import read_te_annotation
from darkdna.toy_data import make_toy_data


def test_te_grammar_features_from_gff3(tmp_path):
    paths = make_toy_data(tmp_path)
    te = read_te_annotation(paths["te_gff3"])
    features = compute_te_grammar_for_region("scaffold_A", 6700, 7200, te)
    assert features["TE_overlap_fraction"] > 0.5
    assert features["TE_family_mosaic_score"] > 0
    assert features["TE_derived_candidate_flag"] is True
