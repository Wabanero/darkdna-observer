import json
import math

from darkdna.features.nonb_dna import compute_nonb_dna_features, compute_nonb_structure_records


def test_nonb_records_are_structure_specific_level_one_evidence():
    sequence = "TTTGGGAGGGAGGGAGGGTTTCCCCACCCCACCCCACCCC"
    records = compute_nonb_structure_records(sequence, region_id="r1")

    types = {record["structure_type"] for record in records}
    assert "G4" in types
    assert "i_motif" in types
    for record in records:
        assert record["region_id"] == "r1"
        assert record["evidence_level"] == "level_1_sequence_potential"
        assert record["formation_context_status"].startswith("unavailable")
        assert record["experimental_support"].startswith("unavailable")
        assert 0 <= record["start_offset"] < record["end_offset"] <= len(sequence)

    features = compute_nonb_dna_features(sequence)
    assert json.loads(features["nonb_structure_records_json"])
    assert math.isnan(features["non_B_DNA_aggregate_score"])

