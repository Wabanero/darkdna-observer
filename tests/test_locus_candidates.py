import pandas as pd

from darkdna.reports.locus_candidates import (
    block_bootstrap_locus_summary,
    merge_candidate_loci,
    write_candidate_locus_outputs,
)


def test_overlapping_windows_merge_to_one_locus(tmp_path):
    windows = pd.DataFrame(
        {
            "region_id": ["r1", "r2", "r3", "r4"],
            "chrom": ["chrI", "chrI", "chrI", "chrI"],
            "start": [0, 500, 0, 9000],
            "end": [1000, 1500, 5000, 10000],
            "window_size": [1000, 1000, 5000, 1000],
            "artifact_risk_flags": ["", "", "", ""],
        }
    )
    labels = pd.DataFrame(
        {
            "region_id": ["r1", "r2", "r3", "r4"],
            "primitive_class": [
                "negative_space_element_candidate",
                "negative_space_element_candidate",
                "negative_space_element_candidate",
                "no_call",
            ],
            "primitive_score_name": [
                "negative_space_element_candidate_score",
                "negative_space_element_candidate_score",
                "negative_space_element_candidate_score",
                "unexplained_dark_anomaly_candidate_score",
            ],
            "primitive_confidence": [0.8, 0.7, 0.6, 0.0],
            "residual_zscore": [3.0, 2.5, 2.2, 0.1],
            "matched_null_zscore": [2.0, 1.5, 1.1, 0.0],
            "empirical_p_value": [0.01, 0.02, 0.03, 0.9],
        }
    )
    residuals = pd.DataFrame(
        {
            "region_id": ["r1", "r2", "r3", "r4"],
            "primitive": [
                "negative_space_element_candidate_score",
                "negative_space_element_candidate_score",
                "negative_space_element_candidate_score",
                "unexplained_dark_anomaly_candidate_score",
            ],
            "observed_score": [1.0, 0.8, 0.7, 0.1],
            "predicted_classical_score": [0.1, 0.1, 0.1, 0.1],
            "residual_score": [0.9, 0.7, 0.6, 0.0],
            "classical_explanation_fraction": [0.1, 0.2, 0.2, 1.0],
        }
    )

    loci = merge_candidate_loci(windows, labels, residuals, block_size=5000)

    assert len(loci) == 1
    locus = loci.iloc[0]
    assert locus["start"] == 0
    assert locus["end"] == 5000
    assert locus["n_windows"] == 3
    assert locus["window_sizes"] == "1000;5000"
    assert locus["scale_validation_status"] == "cross_scale_supported"
    assert locus["locus_effective_test_count"] == 2
    assert locus["locus_empirical_p_value"] == 0.02
    assert locus["global_bh_q_value"] == 0.02

    summary = block_bootstrap_locus_summary(loci, block_size=5000, n_bootstrap=20)
    assert summary.iloc[0]["block_bootstrap_status"] == "insufficient_independent_blocks"

    paths = write_candidate_locus_outputs(windows, labels, residuals, tmp_path, block_size=5000)
    assert paths["candidate_loci_parquet"].exists()
    assert paths["candidate_loci_tsv"].exists()
    assert paths["candidate_loci_bed"].exists()
    assert paths["candidate_loci_block_bootstrap_tsv"].exists()
