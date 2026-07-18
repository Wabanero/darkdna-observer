import pandas as pd

from darkdna.architecture.amount_features import compute_amount_features


def test_amount_features_measure_bp_burden_without_calling_it_function():
    intervals = pd.DataFrame({"region_id": ["r1"], "chrom": ["chr1"], "start": [100], "end": [200]})
    repeats = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [90, 150],
            "end": [120, 230],
            "family": ["A", "B"],
        }
    )
    heterochromatin = pd.DataFrame({"chrom": ["chr1"], "start": [180], "end": [250]})
    features = compute_amount_features(
        intervals,
        repeats=repeats,
        heterochromatin=heterochromatin,
        genome_sizes={"chr1": 1_000},
    )
    row = features.iloc[0]
    assert row["interval_length"] == 100
    assert row["repeat_array_length"] == 70
    assert row["local_repeat_fraction"] == 0.7
    assert row["heterochromatin_overlap"] == 0.2
    assert "descriptive" in row["amount_feature_caveat"]
