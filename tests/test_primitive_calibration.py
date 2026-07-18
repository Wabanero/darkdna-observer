import json

import numpy as np
import pandas as pd

from darkdna.views.primitive_scores import score_primitives


def test_primitive_scores_are_covariance_aware_and_do_not_invent_p_values():
    values = np.linspace(0.1, 2.0, 12)
    features = pd.DataFrame(
        {
            "region_id": [f"r{i}" for i in range(12)],
            "chrom": ["chr1"] * 6 + ["chr2"] * 6,
            "start": list(range(0, 600_000, 50_000)),
            "grammar_entropy": values,
            "forbidden_word_depletion_enrichment": values * 2.0,
            "motif_like_token_recurrence": values[::-1],
            "Markov_order_anomaly": values**2,
            "long_range_dependency_proxy": values + 0.3,
        }
    )
    scores = score_primitives(features)
    column = "constraint_grammar_region_candidate_score"

    assert scores[column].notna().all()
    assert f"{column}_legacy_screening_composite" in scores
    assert scores[f"{column}_empirical_p_value"].isna().all()
    assert set(scores[f"{column}_empirical_p_value_status"]) == {"unavailable"}
    audit = json.loads(scores.iloc[0][f"{column}_correlation_audit"])
    assert audit["effective_components"]
    assert audit["high_correlation_pairs"]
