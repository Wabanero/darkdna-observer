import pandas as pd

from darkdna.nulls.calibration import build_severe_null_panel
from darkdna.nulls.sequence_calibration import SEQUENCE_NULL_MODEL_IDS, build_sequence_transform_null_details


STOCHASTIC_IDS = {
    "mononucleotide_preserving",
    "dinucleotide_preserving_shuffle",
    "kmer_preserving_shuffle",
    "markov_chain_surrogate",
    "synthetic_equal_composition",
}
PRIMITIVE = "constraint_grammar_region_candidate_score"
SEQUENCE = "ATGCGGATCCATGCGGATCCATGCGGATCCTAGCTAGCTAGCTAGCTAGCATGCGGATCCATGCGGATCCTAGCTAGCTA"


def _tables():
    scores = pd.DataFrame({"region_id": ["r1"], PRIMITIVE: [1.4]})
    features = pd.DataFrame(
        {
            "region_id": ["r1"],
            "chrom": ["chr1"],
            "start": [0],
            "length": [len(SEQUENCE)],
            "gc_content": [0.5],
            "sequence": [SEQUENCE],
        }
    )
    return scores, features


def test_sequence_transforms_calibrate_candidate_nulls():
    scores, features = _tables()
    details = build_sequence_transform_null_details(
        scores,
        features,
        {"r1": SEQUENCE},
        [PRIMITIVE],
        n_surrogates=4,
        seed=13,
        include_evolutionary=True,
    )
    by_model = details.set_index("null_model_id")["null_status"].to_dict()
    for model_id in STOCHASTIC_IDS:
        assert by_model[model_id] == "available_sequence_calibrated"
        assert int(details.loc[details["null_model_id"] == model_id, "null_sample_size"].iloc[0]) >= 2
    assert by_model["evolutionary_process_generated"] == "available_sequence_calibrated"
    assert by_model["reversed_sequence"] == "partial_paired_orientation_transform"
    assert by_model["reverse_complement"] == "partial_paired_orientation_transform"


def test_sequence_column_on_feature_table_is_enough_for_candidate_calibration():
    scores, features = _tables()
    summary, details = build_severe_null_panel(
        scores,
        features,
        score_columns=[PRIMITIVE],
        n_controls=2,
        minimum_independent_blocks=5,
        n_sequence_surrogates=4,
        seed=13,
    )
    assert not details.empty
    assert set(SEQUENCE_NULL_MODEL_IDS).issubset(set(details["null_model_id"]))
    row = summary.iloc[0]
    assert "dinucleotide_preserving_shuffle" not in str(row["missing_or_partial_null_models"])
    assert int(row["sequence_null_model_count"]) >= 1
    assert row["null_panel_status"] == "severe_null_panel_available"


def test_missing_sequence_does_not_promote_from_genomic_matching_alone():
    scores = pd.DataFrame({"region_id": ["r1", "r2", "r3"], PRIMITIVE: [3.0, 1.0, 2.0]})
    features = pd.DataFrame(
        {
            "region_id": ["r1", "r2", "r3"],
            "gc_content": [0.4, 0.5, 0.6],
            "length": [100, 100, 100],
        }
    )
    summary, details = build_severe_null_panel(
        scores,
        features,
        score_columns=[PRIMITIVE],
        n_controls=2,
        minimum_independent_blocks=5,
    )
    assert details["null_execution_mode"].eq("matched_table").all()
    assert "dinucleotide_preserving_shuffle" in summary.iloc[0]["missing_or_partial_null_models"]
    assert int(summary.iloc[0]["sequence_null_model_count"]) == 0
    assert summary["null_panel_status"].eq("partial_null_panel_not_for_promotion").all()
