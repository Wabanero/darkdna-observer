import math

from darkdna.features.compression import (
    calibrate_compression,
    compressed_size,
    compressor_header_bytes,
    header_corrected_compressed_size,
)
from darkdna.features.sequence import numeric_walk_mapping_registry


def test_compression_subtracts_headers_and_uses_same_length_nulls():
    sequence = "ACGT" * 60
    result = calibrate_compression(sequence, n_surrogates=6, seed=11)

    assert result["compression_window_size_bp"] == len(sequence)
    assert result["compression_null_replicates"] == 6
    assert result["compression_calibration_status"] in {"available", "partial"}
    for method in ("gzip", "bz2", "lzma"):
        assert header_corrected_compressed_size(sequence, method) == max(
            0, compressed_size(sequence, method) - compressor_header_bytes(method)
        )
        assert f"{method}_dinucleotide_null_compression_zscore" in result
        assert f"{method}_kmer_null_compression_zscore" in result
    assert not math.isnan(result["multiple_compressor_agreement"])


def test_numeric_walk_registry_has_unique_canonical_mapping_tables():
    registry = numeric_walk_mapping_registry()
    mappings = [tuple(sorted(item["mapping_table"].items())) for item in registry]

    assert len(mappings) == len(set(mappings))
    assert {item["name"] for item in registry} == {
        "purine_pyrimidine_numeric_walk",
        "strong_weak_H_bond_numeric_walk",
        "amino_keto_numeric_walk",
    }

