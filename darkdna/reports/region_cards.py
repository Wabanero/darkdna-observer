"""Region card generation."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.primitives.assay_recommender import recommend_assay
from darkdna.primitives.ontology import get_primitive


def top_scores_for_region(residuals: pd.DataFrame, region_id: str, n: int = 5) -> list[dict]:
    subset = residuals[residuals["region_id"].astype(str) == str(region_id)].copy()
    if subset.empty:
        return []
    subset["abs_residual"] = subset["residual_zscore"].abs()
    cols = ["primitive", "observed_score", "residual_zscore", "matched_null_zscore", "empirical_p_value", "classical_explanation_fraction"]
    return subset.sort_values("abs_residual", ascending=False)[cols].head(n).to_dict(orient="records")


def matched_null_summary(residuals: pd.DataFrame, region_id: str) -> dict:
    top = top_scores_for_region(residuals, region_id, n=1)
    if not top:
        return {}
    return {
        "top_primitive": top[0]["primitive"],
        "matched_null_zscore": top[0].get("matched_null_zscore"),
        "empirical_p_value": top[0].get("empirical_p_value"),
    }


def residualization_summary(residuals: pd.DataFrame, region_id: str) -> dict:
    top = top_scores_for_region(residuals, region_id, n=1)
    if not top:
        return {}
    return {
        "observed_score": top[0].get("observed_score"),
        "residual_zscore": top[0].get("residual_zscore"),
        "classical_explanation_fraction": top[0].get("classical_explanation_fraction"),
    }


def build_region_card(
    window: pd.Series,
    label: pd.Series,
    residuals: pd.DataFrame,
    features: pd.DataFrame | None = None,
) -> dict:
    region_id = str(window["region_id"])
    primitive = str(label.get("primitive_class", "unexplained_dark_anomaly_candidate"))
    assay = recommend_assay(primitive)
    ontology = get_primitive(primitive)
    feature_row = None
    if features is not None and not features.empty and region_id in set(features["region_id"].astype(str)):
        feature_row = features.loc[features["region_id"].astype(str) == region_id].iloc[0]
    supporting = str(label.get("top_supporting_features", "")).split(";") if pd.notna(label.get("top_supporting_features", "")) else []
    conflicting = []
    if str(window.get("artifact_risk_flags", "")):
        conflicting.append("artifact_risk_flags_present")
    if float(label.get("primitive_confidence", 0.0) or 0.0) < 0.4:
        conflicting.append("low_confidence")
    controlled = ""
    subset = residuals[residuals["region_id"].astype(str) == region_id]
    if not subset.empty:
        controlled = str(subset.iloc[0].get("covariates_used", ""))
    card = {
        "region_id": region_id,
        "coordinates": f"{window.get('chrom')}:{int(window.get('start'))}-{int(window.get('end'))}",
        "window_size": int(window.get("window_size", int(window.get("end")) - int(window.get("start")))),
        "parent_child_multiscale_context": {
            "parent_region_id": window.get("parent_region_id"),
            "child_region_ids": window.get("child_region_ids"),
        },
        "primitive_class": primitive,
        "candidate_only": True,
        "candidate_statement": "This is a sequence-derived candidate hypothesis, not a confirmed biological primitive.",
        "confirmed_name": ontology.confirmed_name,
        "requires_dynamic_validation": bool(ontology.requires_dynamic_data),
        "required_validation_data": ontology.required_input_level,
        "prompt1_allowed_interpretation": ontology.prompt1_interpretation,
        "suggested_prompt2_view": ontology.suggested_prompt2_view,
        "primitive_confidence": float(label.get("primitive_confidence", 0.0) or 0.0),
        "top_scores": top_scores_for_region(residuals, region_id),
        "top_supporting_features": [f for f in supporting if f],
        "conflicting_features": conflicting,
        "artifact_risk_flags": str(window.get("artifact_risk_flags", "")),
        "matched_null_summary": matched_null_summary(residuals, region_id),
        "residualization_summary": residualization_summary(residuals, region_id),
        "classical_covariates_controlled": controlled,
        "nearest_classical_explanation": "Review TE/repeat/GC/mappability/gene-proximity covariates before prioritization.",
        "why_not_classical": "Candidate remains prioritized by residual and matched-null evidence after available classical covariate controls.",
        "predicted_hidden_property": ontology.description,
        "what_standard_assays_may_miss": "Static annotation or gene-centric assays may miss sequence-intrinsic perturbation interactions.",
        "recommended_primitive_assay": assay.get("assay"),
        "recommended_classical_validation_assay": assay.get("classical_validation"),
        "control_sequence_design": assay.get("sequence_controls"),
        "treatment_perturbation_design": assay.get("perturbations"),
        "key_interaction_test": assay.get("key_interaction_test"),
        "expected_positive_result": "Native sequence shows a larger perturbation-dependent readout shift than matched control sequences.",
        "expected_negative_result": "Native and matched control sequences respond similarly under the perturbation.",
        "allowed_interpretation": assay.get("allowed_interpretation"),
        "forbidden_interpretation": ontology.forbidden_interpretation,
        "assay_feasibility": "MVP estimate: medium; refine with organism, cell type, and available assays.",
        "assay_cost_complexity": "MVP estimate: moderate to high depending on readout.",
        "recommended_next_step": "Inspect matched-null diagnostics, artifact flags, and design native/control perturbation assay.",
        "interpretation_caveat": "This card is an assay-generating hypothesis, not a functional annotation.",
    }
    if assay.get("temporal_interaction_test"):
        card["temporal_interaction_test"] = assay["temporal_interaction_test"]
    return card


def make_region_cards(
    windows: pd.DataFrame,
    labels: pd.DataFrame,
    residuals: pd.DataFrame,
    features: pd.DataFrame | None = None,
    top_n: int | None = None,
) -> list[dict]:
    merged = labels.merge(windows, on="region_id", how="left", suffixes=("_label", ""))
    merged = merged.sort_values("primitive_confidence", ascending=False)
    if top_n:
        merged = merged.head(top_n)
    cards = []
    for row in merged.iterrows():
        record = row[1]
        card = build_region_card(record, record, residuals, features=features)
        cards.append(card)
    return cards


def write_region_cards(cards: list[dict], outdir: str | Path) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    json_path = out / "region_cards.json"
    json_path.write_text(json.dumps(cards, indent=2, default=str), encoding="utf-8")
    table_path = out / "region_cards.tsv"
    flat = pd.DataFrame(
        [
            {
                "region_id": c["region_id"],
                "coordinates": c["coordinates"],
                "primitive_class": c["primitive_class"],
                "primitive_confidence": c["primitive_confidence"],
                "artifact_risk_flags": c["artifact_risk_flags"],
                "recommended_primitive_assay": c["recommended_primitive_assay"],
            }
            for c in cards
        ]
    )
    flat.to_csv(table_path, sep="\t", index=False)
    return {"json": json_path, "tsv": table_path}
