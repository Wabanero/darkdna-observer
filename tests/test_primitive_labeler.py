import math

import pandas as pd

from darkdna.primitives.labeler import assign_primitive_labels
from darkdna.views.primitive_scores import score_primitives


def test_primitive_score_generation_and_assignment():
    features = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "fractal_score": [5.0, 0.1],
            "scale_persistence_score": [1.0, 0.0],
            "compression_anomaly_score": [1.0, 0.0],
        }
    )
    scores = score_primitives(features)
    assert "fractal_scaffold_candidate_score" in scores.columns
    assert "fractal_scaffold_score" not in scores.columns
    residuals = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive": ["fractal_scaffold_candidate_score"],
            "observed_score": [5.0],
            "predicted_classical_score": [0.0],
            "residual_score": [5.0],
            "residual_zscore": [3.0],
            "matched_null_zscore": [3.0],
            "empirical_p_value": [0.01],
            "classical_explanation_fraction": [0.1],
            "covariates_used": ["gc_content"],
        }
    )
    labels = assign_primitive_labels(residuals, features=features)
    assert labels.iloc[0]["primitive_class"] == "fractal_scaffold_candidate"
    assert labels.iloc[0]["candidate_promotion_status"] == "screening_only_legacy_null_metadata_unavailable"
    assert labels.iloc[0]["labeling_status"] == "single_surviving_class"
    assert labels.iloc[0]["is_exclusive_label"]
    assert labels.iloc[0]["competing_primitive_classes"] == "fractal_scaffold_candidate"


def test_primitive_promotion_requires_agreement_across_severe_null_panel():
    residuals = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "primitive": ["fractal_scaffold_candidate_score"] * 2,
            "residual_zscore": [3.0, 3.0],
            "matched_null_zscore": [3.0, 3.0],
            "null_panel_status": ["severe_null_panel_available"] * 2,
            "null_model_count": [4, 4],
            "null_model_agreement": [1.0, 0.5],
            "null_model_conflict": [False, True],
        }
    )

    labels = assign_primitive_labels(residuals)

    status = labels.set_index("region_id")["candidate_promotion_status"].to_dict()
    assert status["r1"] == "eligible_for_candidate_promotion"
    assert status["r2"] == "screening_only_conflicting_null_families"
    assert labels.set_index("region_id").loc["r1", "survives_severe_null_panel"]
    assert not labels.set_index("region_id").loc["r2", "survives_severe_null_panel"]


def test_primitive_labeler_handles_unavailable_null_count():
    residuals = pd.DataFrame(
        {
            "region_id": ["r1"],
            "primitive": ["fractal_scaffold_candidate_score"],
            "residual_zscore": [3.0],
            "matched_null_zscore": [float("nan")],
            "null_model_count": [float("nan")],
        }
    )

    label = assign_primitive_labels(residuals).iloc[0]

    assert label["null_model_count"] == 0
    assert label["candidate_promotion_status"] == "screening_only_legacy_null_metadata_unavailable"
    assert math.isnan(label["matched_null_zscore"])
    assert label["residual_zscore"] == 3.0


def test_labeler_emits_every_surviving_class_instead_of_a_winner():
    residuals = pd.DataFrame(
        {
            "region_id": ["r1", "r1", "r1"],
            "primitive": [
                "constraint_grammar_region_candidate_score",
                "sequence_regime_boundary_candidate_score",
                "unexplained_dark_anomaly_candidate_score",
            ],
            "residual_zscore": [3.2, 2.8, 4.0],
            "matched_null_zscore": [2.1, 2.4, 2.2],
        }
    )

    labels = assign_primitive_labels(residuals)

    assert set(labels["primitive_class"]) == {
        "constraint_grammar_region_candidate",
        "sequence_regime_boundary_candidate",
    }
    assert (labels["labeling_status"] == "competing_hypotheses_not_exclusive").all()
    assert (labels["is_exclusive_label"] == False).all()
    assert labels["competing_primitive_count"].eq(2).all()
    competing = set(labels.iloc[0]["competing_primitive_classes"].split(";"))
    assert competing == {
        "constraint_grammar_region_candidate",
        "sequence_regime_boundary_candidate",
    }


def test_labeler_keeps_unavailable_zscores_as_na_and_does_not_treat_them_as_zero():
    residuals = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "primitive": ["fractal_scaffold_candidate_score", "fractal_scaffold_candidate_score"],
            "residual_zscore": [float("nan"), 0.4],
            "matched_null_zscore": [float("nan"), 0.3],
        }
    )

    labels = assign_primitive_labels(residuals).set_index("region_id")

    assert labels.loc["r1", "primitive_class"] == "no_call"
    assert labels.loc["r1", "labeling_status"] == "no_call"
    assert math.isnan(labels.loc["r1", "residual_zscore"])
    assert math.isnan(labels.loc["r1", "matched_null_zscore"])
    assert math.isnan(labels.loc["r1", "primitive_priority"])
    assert labels.loc["r2", "primitive_class"] == "no_call"
    assert labels.loc["r2", "residual_zscore"] == 0.4
    assert labels.loc["r2", "matched_null_zscore"] == 0.3
