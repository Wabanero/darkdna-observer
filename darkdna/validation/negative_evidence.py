"""First-class negative evidence for conservative candidate downgrading."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class NegativeEvidenceRecord:
    criterion: str
    status: str
    decisive: bool
    observed_value: str
    interpretation: str
    source: str


DECISIVE_STATUSES = {"failed", "artifact_supported", "explained_by_null", "negative_sufficient_power"}


def _record(
    criterion: str,
    status: str,
    observed_value: Any,
    interpretation: str,
    source: str,
    *,
    decisive: bool | None = None,
) -> NegativeEvidenceRecord:
    is_decisive = status in DECISIVE_STATUSES if decisive is None else decisive
    return NegativeEvidenceRecord(
        criterion=criterion,
        status=status,
        decisive=bool(is_decisive),
        observed_value=str(observed_value),
        interpretation=interpretation,
        source=source,
    )


def evaluate_negative_evidence(evidence: dict[str, Any]) -> dict[str, Any]:
    """Evaluate supplied negative controls without treating missing data as failure."""

    records: list[NegativeEvidenceRecord] = []
    null_survival = evidence.get("survives_severe_nulls")
    if null_survival is False:
        records.append(_record("severe_null_panel", "explained_by_null", False, "The signal disappears under the severe null panel.", "null_models"))
    elif null_survival is True:
        records.append(_record("severe_null_panel", "passed", True, "The signal survives the configured severe null panel; this is anomaly evidence, not function.", "null_models", decisive=False))
    else:
        records.append(_record("severe_null_panel", "unavailable", "NA", "A complete severe null panel is unavailable; absence of evidence is not negative evidence.", "null_models", decisive=False))

    if evidence.get("te_family_age_explains_signal") is True:
        records.append(_record("TE_family_and_age", "explained_by_null", True, "TE family/subfamily and age explain the candidate signal.", "TE_annotation"))
    if evidence.get("window_shift_stable") is False:
        records.append(_record("nearby_window_shift_stability", "failed", False, "The signal is unstable under nearby window shifts.", "multiscale_diagnostics"))
    if evidence.get("assembly_artifact_supported") is True:
        records.append(_record("assembly_representation", "artifact_supported", True, "Assembly gap/collapse or repeat representation can explain the signal.", "artifact_flags"))
    if evidence.get("survives_mutation_environment_null") is False:
        records.append(_record("organism_specific_mutation_null", "explained_by_null", False, "The signal disappears under organism-specific mutation controls.", "evolutionary_null"))
    if evidence.get("mechanistic_bridge_plausible") is False:
        records.append(_record("mechanistic_bridge", "failed", False, "No plausible intermediate molecular process has been specified.", "region_card"))
    if evidence.get("endogenous_perturbation_negative_sufficient_power") is True:
        records.append(_record("endogenous_perturbation", "negative_sufficient_power", True, "A sufficiently powered endogenous perturbation was negative.", "causal_experiment"))
    if evidence.get("matched_rescue_negative") is True:
        records.append(_record("matched_rescue", "negative_sufficient_power", True, "The matched rescue failed with adequate power.", "causal_experiment"))
    if evidence.get("independent_replication_failed") is True:
        records.append(_record("independent_replication", "failed", True, "Independent replication failed.", "replication_study"))

    decisive = [record for record in records if record.decisive]
    unavailable = [record for record in records if record.status == "unavailable"]
    if decisive:
        candidate_status = "candidate_not_supported"
        decision = "reject_or_downgrade"
    elif unavailable:
        candidate_status = "insufficient_evidence"
        decision = "hold_without_escalation"
    else:
        candidate_status = "unresolved"
        decision = "eligible_for_next_predeclared_validation_stage"
    return {
        "candidate_status": candidate_status,
        "decision": decision,
        "decisive_negative_count": len(decisive),
        "records": [asdict(record) for record in records],
        "caveat": "Negative controls can reject or downgrade a candidate; missing controls do not count as positive support.",
    }


def candidate_negative_evidence(
    *,
    artifact_risk_flags: str = "",
    null_panel_status: str = "",
    survives_severe_nulls: bool | None = None,
    window_shift_status: str = "",
    bridge_status: str = "",
) -> dict[str, Any]:
    flags = {flag.strip() for flag in str(artifact_risk_flags).replace(",", ";").split(";") if flag.strip()}
    assembly_flags = {
        "assembly_gap",
        "assembly_gap_overlap",
        "overlaps_assembly_gap",
        "collapsed_repeat",
        "assembly_collapse",
        "very_high_n_fraction",
    }
    legacy_complete_null = null_panel_status in {"complete_severe_null_panel_passed", "available_multiple_severe_nulls"}
    if survives_severe_nulls is None and legacy_complete_null:
        survives_severe_nulls = True
    bridge_plausible = bridge_status not in {"bridge_missing_feature_audit_required", "missing", ""}
    window_stable: bool | None
    if window_shift_status == "available":
        window_stable = True
    elif window_shift_status in {"failed", "unstable"}:
        window_stable = False
    else:
        window_stable = None
    return evaluate_negative_evidence(
        {
            "survives_severe_nulls": survives_severe_nulls,
            "window_shift_stable": window_stable,
            "assembly_artifact_supported": bool(flags.intersection(assembly_flags)),
            "mechanistic_bridge_plausible": bridge_plausible,
        }
    )


def write_negative_evidence(rows: list[dict[str, Any]], outdir: str | Path) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    json_path = out / "negative_evidence.json"
    parquet_path = out / "negative_evidence.parquet"
    json_path.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    flat: list[dict[str, Any]] = []
    for row in rows:
        region_id = row.get("region_id", "")
        for record in row.get("records", []):
            flat.append({"region_id": region_id, "candidate_status": row.get("candidate_status"), "decision": row.get("decision"), **record})
    pd.DataFrame(flat).to_parquet(parquet_path, index=False)
    return {"json": json_path, "parquet": parquet_path}
