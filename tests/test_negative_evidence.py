from darkdna.validation.negative_evidence import evaluate_negative_evidence


def test_decisive_negative_evidence_rejects_or_downgrades_candidate():
    result = evaluate_negative_evidence(
        {
            "survives_severe_nulls": False,
            "assembly_artifact_supported": True,
            "mechanistic_bridge_plausible": True,
        }
    )

    assert result["candidate_status"] == "candidate_not_supported"
    assert result["decision"] == "reject_or_downgrade"
    assert result["decisive_negative_count"] >= 2


def test_missing_controls_are_insufficient_not_positive_support():
    result = evaluate_negative_evidence({})

    assert result["candidate_status"] == "insufficient_evidence"
    assert result["decision"] == "hold_without_escalation"

