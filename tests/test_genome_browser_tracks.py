import pandas as pd

from darkdna.reports.genome_browser_tracks import make_tracks


def test_genome_browser_tracks_are_written(tmp_path):
    windows = pd.DataFrame({"region_id": ["r1"], "chrom": ["scaffold_A"], "start": [0], "end": [200], "artifact_risk_flags": ["high_N_fraction"]})
    labels = pd.DataFrame({"region_id": ["r1"], "primitive_class": ["negative_space_element_candidate"], "primitive_confidence": [0.9]})
    residuals = pd.DataFrame({"region_id": ["r1"], "primitive": ["unexplained_dark_anomaly_candidate_score"], "residual_zscore": [4.0]})
    paths = make_tracks(windows, labels, residuals, tmp_path)
    assert paths["all_residual_scores"].exists()
    assert paths["negative_space_element_candidate"].exists()
    assert paths["artifact_risk_flags"].exists()
