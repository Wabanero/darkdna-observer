import pandas as pd

from darkdna.nulls.calibration import build_severe_null_panel
from darkdna.nulls.registry import assess_null_availability, null_model_registry


def test_severe_null_registry_is_complete_and_candidate_summary_is_conservative():
    registry = null_model_registry()
    ids = {row["null_model_id"] for row in registry}
    assert {
        "same_length_gc_matched",
        "dinucleotide_preserving_shuffle",
        "te_subfamily_matched",
        "assembly_confidence_matched",
        "evolutionary_process_generated",
    }.issubset(ids)
    availability = assess_null_availability(
        ["chrom", "start", "length", "gc_content", "mappability"], sequence_available=True
    )
    by_id = {row["null_model_id"]: row for row in availability}
    assert by_id["mononucleotide_preserving"]["available"]
    assert not by_id["te_subfamily_matched"]["available"]

    scores = pd.DataFrame({"region_id": [f"r{i}" for i in range(8)], "candidate_score": list(range(8))})
    features = pd.DataFrame(
        {
            "region_id": [f"r{i}" for i in range(8)],
            "chrom": ["chr1"] * 8,
            "start": [i * 100_000 for i in range(8)],
            "length": [100] * 8,
            "gc_content": [0.4 + i * 0.01 for i in range(8)],
            "mappability": [0.9] * 8,
        }
    )
    summary, details = build_severe_null_panel(
        scores,
        features,
        score_columns=["candidate_score"],
        block_size_bp=100_000,
        minimum_independent_blocks=3,
        n_controls=6,
    )
    assert not details.empty
    assert {"available_null_models", "missing_null_models", "null_model_count", "null_model_agreement", "null_model_conflict"}.issubset(summary.columns)
    assert summary["null_model_id"].eq("severe_null_panel_conservative_aggregate").all()
