import pandas as pd

from darkdna.architecture.spacing import compute_spacing_features


def test_spacing_reports_anchor_distance_not_intervening_sequence_function():
    intervals = pd.DataFrame(
        {
            "region_id": ["r1", "r2", "r3"],
            "chrom": ["chr1"] * 3,
            "start": [100, 400, 800],
            "end": [200, 500, 900],
        }
    )
    anchors = pd.DataFrame(
        {
            "chrom": ["chr1"] * 4,
            "start": [0, 250, 600, 1_000],
            "end": [50, 300, 650, 1_050],
        }
    )
    result = compute_spacing_features(intervals, anchors)
    assert result.loc[0, "left_anchor_distance"] == 50
    assert result.loc[0, "right_anchor_distance"] == 50
    assert result["spacing_null_zscore"].notna().all()
    assert "Comparative spacing" in result.loc[0, "spacing_reason"]
