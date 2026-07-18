import numpy as np
import pandas as pd

from darkdna.architecture.comparative_length import compute_length_conservation


def test_length_conservation_requires_confident_synteny_not_poor_alignment():
    intervals = pd.DataFrame({"region_id": ["r1", "r2"]})
    syntenic = pd.DataFrame(
        {
            "region_id": ["r1", "r1", "r1", "r2"],
            "interval_length": [100, 102, 98, 100],
            "sequence_identity": [0.2, 0.25, 0.22, 0.1],
            "alignment_coverage": [0.9, 0.9, 0.8, 0.1],
            "synteny_confidence": [0.9, 0.9, 0.8, 0.1],
            "assembly_confidence": [0.9, 0.9, 0.9, 0.2],
        }
    )
    result = compute_length_conservation(intervals, syntenic).set_index("region_id")
    assert result.loc["r1", "length_conservation_status"] == "available_confident_synteny"
    assert result.loc["r1", "length_minus_sequence_conservation"] > 0
    assert result.loc["r2", "length_conservation_status"] == "unavailable_low_synteny_or_alignment_confidence"
    assert np.isnan(result.loc["r2", "length_conservation_score"])
