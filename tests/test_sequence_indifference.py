from darkdna.architecture.sequence_indifference import evaluate_sequence_indifference


def _motif_predictor(sequence, length, copy_number, spacing):
    return sequence.count("AAAA") / max(1, length)


def _length_predictor(sequence, length, copy_number, spacing):
    return length / 100.0


def _copy_predictor(sequence, length, copy_number, spacing):
    return copy_number


def test_sequence_specific_length_specific_and_copy_specific_synthetic_mechanisms_separate():
    sequence = "AAAA" * 25 + "CGCG" * 25
    replacement_pool = ["ACGT" * 50, "GGCC" * 50]
    motif, _ = evaluate_sequence_indifference("motif", sequence, predictor=_motif_predictor, replacement_pool=replacement_pool)
    length, _ = evaluate_sequence_indifference("length", sequence, predictor=_length_predictor, replacement_pool=replacement_pool)
    copy, controls = evaluate_sequence_indifference("copy", sequence, predictor=_copy_predictor, replacement_pool=replacement_pool)
    assert motif["sequence_identity_sensitivity"] > motif["length_sensitivity"]
    assert length["length_sensitivity"] > length["sequence_identity_sensitivity"]
    assert copy["copy_number_sensitivity"] > copy["sequence_identity_sensitivity"]
    assert {row["copy_number"] for row in controls if row["control_family"] == "copy_number"} == {0.0, 1.0, 2.0, 4.0, 8.0}
    assert "not biological causality" in copy["sequence_indifference_caveat"]
