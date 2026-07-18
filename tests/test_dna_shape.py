import json
import math

from darkdna.features.physical_shape import compute_physical_shape_features


def test_validated_shape_fields_are_records_and_missing_predictors_are_na():
    result = compute_physical_shape_features("ACGT" * 20)
    records = json.loads(result["dna_shape_records_json"])

    assert {record["feature"] for record in records} >= {
        "minor_groove_width",
        "roll",
        "propeller_twist",
        "helix_twist",
        "curvature",
        "bendability",
    }
    assert all(
        {"predictor", "predictor_version", "raw_value", "normalization", "availability_status", "warning"}.issubset(record)
        for record in records
    )
    assert math.isnan(result["DNA_stiffness_proxy"])
    assert result["DNA_stiffness_proxy_deprecation_warning"].startswith("Unavailable")
    assert "legacy" in result["legacy_dinucleotide_shape_warning"].lower() or "deprecated" in result["legacy_dinucleotide_shape_warning"].lower()

