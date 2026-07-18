import pandas as pd
import pytest

from darkdna.architecture.presence_absence import compute_presence_absence_features


def test_presence_absence_supports_wide_accession_matrix():
    intervals = pd.DataFrame({"region_id": ["r1", "r2"]})
    matrix = pd.DataFrame({"region_id": ["r1", "r2"], "a1": [1, 0], "a2": [1, 1], "a3": [0, 0]})
    features = compute_presence_absence_features(intervals, matrix).set_index("region_id")
    assert features.loc["r1", "presence_absence_frequency"] == 2 / 3
    assert features.loc["r2", "deletion_frequency_from_presence_absence"] == pytest.approx(2 / 3)
    assert features.loc["r1", "presence_absence_sample_size"] == 3
