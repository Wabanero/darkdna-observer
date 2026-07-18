import pandas as pd

from darkdna.architecture.pipeline import run_sequence_indifferent_architecture, write_architecture_outputs


def test_mode_b_pipeline_writes_separate_axes_and_explicit_missing_inputs(tmp_path):
    intervals = pd.DataFrame(
        {
            "region_id": ["r1", "r2", "r3", "r4", "r5", "r6"],
            "chrom": ["chr1"] * 6,
            "start": [0, 100_000, 200_000, 300_000, 400_000, 500_000],
            "end": [100, 100_120, 200_140, 300_160, 400_180, 500_200],
            "artifact_risk_flags": [""] * 6,
        }
    )
    sequences = {row.region_id: ("ACGT" * 60)[: row.end - row.start] for row in intervals.itertuples()}
    results = run_sequence_indifferent_architecture(
        intervals,
        sequences=sequences,
        block_size_bp=100_000,
        minimum_independent_blocks=3,
    )
    candidates = results["architecture_candidates"]
    assert {
        "sequence_indifferent_candidate",
        "dominant_mode",
        "sequence_identity_sensitivity",
        "length_sensitivity",
        "copy_number_sensitivity",
        "promotion_status",
    }.issubset(candidates.columns)
    assert results["copy_number_features"]["copy_number_status"].eq("unavailable_missing_copy_number_input").all()
    comparison = results["sequence_vs_quantity_scores"]
    assert comparison["comparison_caveat"].str.contains("separate evidence axes").all()
    assert {
        "matched_interval_independent_blocks",
        "equal_length_replacement_panel",
        "length_titration_panel",
        "copy_number_titration_panel",
    }.issubset(set(results["architecture_nulls"]["null_model_id"]))
    transformed = results["architecture_nulls"].query("null_model_id != 'matched_interval_independent_blocks'")
    assert transformed["architecture_null_empirical_p"].isna().all()
    paths = write_architecture_outputs(results, tmp_path)
    assert all(path.exists() for path in paths.values())
