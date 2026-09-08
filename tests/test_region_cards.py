import pandas as pd

from darkdna.reports.region_cards import make_region_cards


def test_region_card_contains_assay_blueprint_fields():
    windows = pd.DataFrame(
        {
            "region_id": ["r1"],
            "chrom": ["scaffold_A"],
            "start": [10],
            "end": [210],
            "window_size": [200],
            "parent_region_id": [None],
            "child_region_ids": [""],
            "artifact_risk_flags": [""],
        }
    )
    labels = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive_class": ["negative_space_element_candidate"],
            "primitive_score_name": ["negative_space_element_candidate_score"],
            "primitive_confidence": [0.8],
            "top_supporting_features": ["depleted_kmer_score"],
        }
    )
    residuals = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive": ["negative_space_element_candidate_score"],
            "observed_score": [2.0],
            "residual_zscore": [3.0],
            "matched_null_zscore": [2.5],
            "empirical_p_value": [0.01],
            "classical_explanation_fraction": [0.2],
            "covariates_used": ["gc_content"],
        }
    )
    cards = make_region_cards(windows, labels, residuals)
    assert cards[0]["key_interaction_test"].startswith("negative_space_rescue_effect")
    assert cards[0]["candidate_only"] is True
    assert cards[0]["forbidden_interpretation"]
    assert cards[0]["recommended_primitive_assay"] == "Negative-Space Rescue/Scramble Assay"
    assert cards[0]["observed_feature_evidence"]["supporting_features"] == ["depleted_kmer_score"]
    assert cards[0]["primitive_hypothesis"]["candidate_label"] == "negative_space_element_candidate"
    assert cards[0]["mechanistic_bridge"]["bridge_status"] == "hypothesized_bridge_requires_validation"
    assert cards[0]["mechanistic_bridge"]["direct_primitive_validation_allowed"] is False
    assert "exploratory_only" in cards[0]["assay_validation_scope"]
    assert "endogenous_CRISPRi_or_CRISPRa" in cards[0]["causal_validation_hierarchy"]
    assert "genomic position" in cards[0]["native_context_caveat"]
    assert "not observed molecular properties" in cards[0]["feature_hypothesis_boundary"]
    assert cards[0]["labeling_status"] in {"", "single_surviving_class", "competing_hypotheses_not_exclusive", "no_call", "unexplained_residual_without_dominant_class"}
    assert "Dark is an operational project term" in cards[0]["terminology_scope"]["dark_operational_use"]
    assert cards[0]["assembly_pangenome_context"]["current_scope"] == "reference_based_window"
    assert cards[0]["score_methodology"]["score_status"] == "covariance_aware_cohort_standardized_screening_score"
    assert "depleted_kmer_score" in cards[0]["score_methodology"]["component_features"]
    assert cards[0]["null_model_panel"]["status"] == "insufficient_single_matched_null_until_complementary_nulls_pass"
    assert cards[0]["negative_evidence"]["candidate_status"] == "insufficient_evidence"
    assert "does not establish" in cards[0]["evolutionary_interpretation_caveat"]


def test_region_card_uses_primitive_specific_key_tests():
    windows = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "chrom": ["scaffold_A", "scaffold_A"],
            "start": [10, 300],
            "end": [210, 500],
            "window_size": [200, 200],
            "parent_region_id": [None, None],
            "child_region_ids": ["", ""],
            "artifact_risk_flags": ["", ""],
        }
    )
    labels = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "primitive_class": ["fractal_scaffold_candidate", "non_B_DNA_physical_susceptibility_candidate"],
            "primitive_confidence": [0.8, 0.8],
            "top_supporting_features": ["fractal_score", "G4_susceptibility_proxy"],
        }
    )
    residuals = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "primitive": ["fractal_scaffold_candidate_score", "non_B_DNA_physical_susceptibility_candidate_score"],
            "observed_score": [2.0, 2.0],
            "residual_zscore": [3.0, 3.0],
            "matched_null_zscore": [2.5, 2.5],
            "empirical_p_value": [0.01, 0.01],
            "classical_explanation_fraction": [0.2, 0.2],
            "covariates_used": ["gc_content", "gc_content"],
        }
    )
    cards = make_region_cards(windows, labels, residuals)
    tests_by_primitive = {card["primitive_class"]: card["key_interaction_test"] for card in cards}
    bridge_by_primitive = {card["primitive_class"]: card["mechanistic_bridge"] for card in cards}
    assert tests_by_primitive["fractal_scaffold_candidate"].startswith("folding_scale_effect")
    assert tests_by_primitive["non_B_DNA_physical_susceptibility_candidate"].startswith("physical_susceptibility_effect")
    assert len(set(tests_by_primitive.values())) == 2
    assert "polymer physics or coarse-grained DNA simulations" in bridge_by_primitive["fractal_scaffold_candidate"]["candidate_intermediate_processes"]
    assert "test GC-, dinucleotide-, and k-mer-preserved controls" in bridge_by_primitive["fractal_scaffold_candidate"]["required_bridge_evidence"]
    assert "quantum" not in bridge_by_primitive["non_B_DNA_physical_susceptibility_candidate"]["measured_feature"].lower()


def test_no_call_and_unexplained_do_not_use_generic_key_test():
    windows = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "chrom": ["scaffold_A", "scaffold_A"],
            "start": [10, 300],
            "end": [210, 500],
            "window_size": [200, 200],
            "parent_region_id": [None, None],
            "child_region_ids": ["", ""],
            "artifact_risk_flags": ["", ""],
        }
    )
    labels = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "primitive_class": ["no_call", "unexplained_dark_anomaly_candidate"],
            "primitive_confidence": [0.0, 0.8],
            "top_supporting_features": ["", ""],
        }
    )
    residuals = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "primitive": ["unexplained_dark_anomaly_candidate_score", "unexplained_dark_anomaly_candidate_score"],
            "observed_score": [0.1, 2.0],
            "residual_zscore": [0.1, 3.0],
            "matched_null_zscore": [0.1, 2.5],
            "empirical_p_value": [0.9, 0.01],
            "classical_explanation_fraction": [0.8, 0.2],
            "covariates_used": ["gc_content", "gc_content"],
        }
    )
    cards = make_region_cards(windows, labels, residuals)
    tests_by_primitive = {card["primitive_class"]: card["key_interaction_test"] for card in cards}
    assert tests_by_primitive["no_call"].startswith("no_call_review")
    assert tests_by_primitive["unexplained_dark_anomaly_candidate"].startswith("dark_anomaly_effect")
    assert "Native_treatment" not in " ".join(tests_by_primitive.values())


def test_region_cards_keep_observed_features_separate_from_hypothesis():
    windows = pd.DataFrame(
        {
            "region_id": ["r1"],
            "chrom": ["scaffold_A"],
            "start": [10],
            "end": [210],
            "window_size": [200],
            "parent_region_id": [None],
            "child_region_ids": [""],
            "artifact_risk_flags": ["low_mappability"],
            "mappability": [0.2],
            "n_fraction": [0.01],
            "low_complexity_mask_fraction": [0.1],
            "scaffold_edge_distance": [5000],
            "overlaps_assembly_gap": [False],
            "overlaps_segmental_duplication": [True],
        }
    )
    labels = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive_class": ["sequence_regime_boundary_candidate"],
            "primitive_score_name": ["sequence_regime_boundary_candidate_score"],
            "primitive_confidence": [0.7],
            "top_supporting_features": ["entropy_boundary_score;gc_left_right_delta"],
        }
    )
    residuals = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive": ["sequence_regime_boundary_candidate_score"],
            "observed_score": [1.8],
            "residual_zscore": [2.7],
            "matched_null_zscore": [2.3],
            "empirical_p_value": [0.02],
            "classical_explanation_fraction": [0.25],
            "covariates_used": ["gc_content;mappability"],
            "null_panel_status": ["insufficient_single_matched_null_until_complementary_nulls_pass"],
            "available_null_models": ["matched_controls_v1"],
            "missing_or_partial_null_models": ["dinucleotide_preserving_shuffle,kmer_preserving_shuffle"],
        }
    )

    card = make_region_cards(windows, labels, residuals)[0]

    observed = card["observed_feature_evidence"]
    hypothesis = card["primitive_hypothesis"]
    assert observed["supporting_features"] == ["entropy_boundary_score", "gc_left_right_delta"]
    assert observed["artifact_risk_flags"] == "low_mappability"
    assert card["assembly_pangenome_context"]["available_reference_context"]["mappability"] == 0.2
    assert card["assembly_pangenome_context"]["available_reference_context"]["overlaps_segmental_duplication"] is True
    assert "presence_absence_variation" in card["assembly_pangenome_context"]["missing_first_class_inputs"]
    assert card["score_methodology"]["component_features"] == ["boundary_condition_candidate_score"]
    assert "dinucleotide_preserving_shuffle" in card["null_model_panel"]["missing_or_partial_null_models"]
    assert "kmer_preserving_shuffle" in card["null_model_panel"]["missing_or_partial_null_models"]
    assert hypothesis["candidate_label"] == "sequence_regime_boundary_candidate"
    assert hypothesis["hypothesis_statement"] != observed["supporting_features"][0]
    assert card["candidate_only"] is True


def test_periodic_spacing_candidate_requires_intermediate_mechanistic_bridge():
    windows = pd.DataFrame(
        {
            "region_id": ["r1"],
            "chrom": ["scaffold_A"],
            "start": [10],
            "end": [210],
            "window_size": [200],
            "parent_region_id": [None],
            "child_region_ids": [""],
            "artifact_risk_flags": [""],
        }
    )
    labels = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive_class": ["periodic_spacing_grammar_candidate"],
            "primitive_confidence": [0.8],
            "top_supporting_features": ["phase_periodicity_around_10bp"],
        }
    )
    residuals = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive": ["periodic_spacing_grammar_candidate_score"],
            "observed_score": [2.0],
            "residual_zscore": [3.0],
            "matched_null_zscore": [2.5],
            "empirical_p_value": [0.01],
            "classical_explanation_fraction": [0.2],
            "covariates_used": ["gc_content"],
        }
    )

    card = make_region_cards(windows, labels, residuals)[0]
    bridge = card["mechanistic_bridge"]

    assert "10 bp/147 bp periodicity" in bridge["measured_feature"]
    assert "spacing-dependent molecular readout" in bridge["proposed_dynamic_phenotype"]
    assert "nucleosome phasing" in bridge["candidate_intermediate_processes"]
    assert "TF cooperative binding at phased spacing" in bridge["candidate_intermediate_processes"]
    assert "show periodicity affects an intermediate molecular process" in bridge["required_bridge_evidence"]
    assert "pulse decoding" in " ".join(bridge["required_bridge_evidence"])
    assert bridge["direct_primitive_validation_allowed"] is False
    assert "exploratory_only" in card["assay_validation_scope"]


def test_missing_confidence_is_not_treated_as_low_confidence():
    windows = pd.DataFrame(
        {
            "region_id": ["r1"],
            "chrom": ["scaffold_A"],
            "start": [10],
            "end": [210],
            "window_size": [200],
            "parent_region_id": [None],
            "child_region_ids": [""],
            "artifact_risk_flags": [""],
        }
    )
    labels = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive_class": ["negative_space_element_candidate"],
            "primitive_confidence": [float("nan")],
            "labeling_status": ["single_surviving_class"],
            "competing_primitive_classes": ["negative_space_element_candidate"],
            "competing_primitive_count": [1],
            "is_exclusive_label": [True],
            "top_supporting_features": ["depleted_kmer_score"],
        }
    )
    residuals = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive": ["negative_space_element_candidate_score"],
            "observed_score": [2.0],
            "residual_zscore": [3.0],
            "matched_null_zscore": [2.5],
            "empirical_p_value": [0.01],
            "classical_explanation_fraction": [0.2],
            "covariates_used": ["gc_content"],
        }
    )

    card = make_region_cards(windows, labels, residuals)[0]
    assert card["primitive_confidence"] is None
    assert "low_confidence" not in card["conflicting_features"]
    assert card["is_exclusive_label"] is True
    assert card["competing_primitive_classes"] == ["negative_space_element_candidate"]
