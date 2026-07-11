import json

import pandas as pd

from darkdna.reports.html_report import generate_html_report


def test_html_report_generation(tmp_path):
    windows = pd.DataFrame({"region_id": ["r1"], "chrom": ["scaffold_A"], "start": [0], "end": [200], "window_size": [200], "artifact_risk_flags": [""]})
    labels = pd.DataFrame({"region_id": ["r1"], "primitive_class": ["negative_space_element_candidate"], "primitive_confidence": [0.9]})
    residuals = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "primitive": ["negative_space_element_candidate_score", "TE_grammar_node_candidate_score"],
            "observed_score": [1.0, 0.0],
            "predicted_classical_score": [0.2, 0.0],
            "residual_score": [0.8, 0.0],
            "residual_zscore": [3.0, 0.0],
            "matched_null_zscore": [2.0, 0.0],
        }
    )
    cards = [
        {
            "region_id": "r1",
            "primitive_class": "negative_space_element_candidate",
            "coordinates": "scaffold_A:0-200",
            "primitive_confidence": 0.9,
            "candidate_only": True,
            "candidate_statement": "This is a sequence-derived candidate hypothesis, not a confirmed biological primitive.",
            "observed_feature_evidence": {"supporting_features": ["depleted_kmer_score"]},
            "primitive_hypothesis": {"hypothesis_statement": "Sequence-derived negative-space substrate candidate."},
            "terminology_scope": {"dark_operational_use": "Dark is operational."},
            "assembly_pangenome_context": {"caveat": "Reference-scoped candidate."},
            "score_methodology": {
                "score_status": "uncalibrated_equal_weight_screening_composite",
                "caveat": "Screening composite.",
            },
            "null_model_panel": {
                "status": "insufficient_single_matched_null_until_complementary_nulls_pass",
                "caveat": "Single null is insufficient.",
            },
            "recommended_primitive_assay": "Negative-Space Rescue/Scramble Assay",
            "key_interaction_test": "effect = test",
            "interpretation_caveat": "hypothesis",
        }
    ]
    path = generate_html_report(windows, labels, residuals, cards, tmp_path)
    assert path.exists()
    html = path.read_text(encoding="utf-8")
    assert "hypothesis-generating" in html
    assert "Locus-Level Candidate Evidence" in html
    assert "Overlapping multiscale windows are dependent observations" in html
    assert "Observed feature evidence" in html
    assert "Primitive hypothesis" in html
    assert "Assembly context" in html
    assert "Terminology" in html
    assert "Score methodology" in html
    assert "Null model panel" in html
    assert "Prompt" not in html
    assert (tmp_path / "multipanel_summary.svg").exists()
    assert (tmp_path / "classical_control_multipanel.svg").exists()
    assert (tmp_path / "multipanel_summary.svg").read_text(encoding="utf-8").count("Candidate Summary Multipanel") == 1
    assert (tmp_path / "classical_control_multipanel.svg").read_text(encoding="utf-8").count("Classical Explanation Removal Multipanel") == 1
    assert "No active TE-grammar candidates" in html
    assert "TE_grammar_node_candidate_score</td>" not in html
