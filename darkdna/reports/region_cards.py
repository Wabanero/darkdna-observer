"""Region card generation."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.primitives.assay_recommender import recommend_assay
from darkdna.primitives.ontology import get_primitive
from darkdna.residuals.null_models import null_panel_status
from darkdna.utils.progress import ProgressReporter
from darkdna.validation.negative_evidence import candidate_negative_evidence
from darkdna.views.primitive_scores import primitive_score_manifest


TERMINOLOGY_SCOPE = {
    "studied_sequence_scope": "unannotated_and_noncoding_genomic_sequence_architectures",
    "dark_operational_use": (
        "Dark is an operational project term; it does not imply function, "
        "selected-effect status, assembly absence, or assembly difficulty."
    ),
    "function_caveat": (
        "Biochemical activity such as transcription, accessibility, or protein binding "
        "is not sufficient evidence of selected biological function."
    ),
}

CAUSAL_VALIDATION_HIERARCHY = [
    "in_silico_mutagenesis_or_sequence_model_sanity_check",
    "MPRA_or_STARR_seq_or_synthetic_reporter_library",
    "endogenous_CRISPRi_or_CRISPRa",
    "small_deletion_base_editing_or_prime_editing",
    "single_cell_perturbation_RNA_ATAC_or_multiome",
    "knockin_or_knockout_model_when_justified",
    "phenotype_under_challenge_ageing_stress_differentiation_or_development",
]


def _finite_float(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def _bool_field(row: pd.Series, field: str) -> bool:
    value = row.get(field, False)
    return bool(value) if pd.notna(value) else False


def assembly_pangenome_context(window: pd.Series) -> dict:
    """Summarize current reference-representation evidence and missing pangenome layers."""

    return {
        "current_scope": "reference_based_window",
        "available_reference_context": {
            "mappability": _finite_float(window.get("mappability")),
            "n_fraction": _finite_float(window.get("n_fraction")),
            "low_complexity_mask_fraction": _finite_float(window.get("low_complexity_mask_fraction")),
            "scaffold_edge_distance": _finite_float(window.get("scaffold_edge_distance")),
            "overlaps_assembly_gap": _bool_field(window, "overlaps_assembly_gap"),
            "overlaps_segmental_duplication": _bool_field(window, "overlaps_segmental_duplication"),
            "artifact_risk_flags": str(window.get("artifact_risk_flags", "")),
        },
        "missing_first_class_inputs": [
            "assembly_confidence",
            "copy_number_variation",
            "presence_absence_variation",
            "strain_or_accession_specificity",
            "repeat_array_completeness",
            "graph_or_pangenome_coordinates",
            "assembly_liftover",
        ],
        "caveat": (
            "Candidate interpretation is reference-scoped until assembly confidence, "
            "copy-number, presence/absence, and pangenome context are available."
        ),
    }


def top_scores_for_region(residuals: pd.DataFrame, region_id: str, n: int = 5) -> list[dict]:
    subset = residuals[residuals["region_id"].astype(str) == str(region_id)].copy()
    if subset.empty:
        return []
    subset["abs_residual"] = subset["residual_zscore"].abs()
    cols = [
        column
        for column in [
            "primitive",
            "observed_score",
            "residual_zscore",
            "residual_zscore_method",
            "matched_null_zscore",
            "empirical_p_value",
            "empirical_p_value_status",
            "classical_model_global_r2",
            "classical_explanation_fraction",
        ]
        if column in subset.columns
    ]
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
        "residual_zscore_method": top[0].get("residual_zscore_method"),
        "classical_model_global_r2": top[0].get("classical_model_global_r2", top[0].get("classical_explanation_fraction")),
        "classical_explanation_fraction": top[0].get("classical_explanation_fraction"),
        "migration_warning": "classical_explanation_fraction is a deprecated alias for the model-level classical_model_global_r2.",
    }


def score_methodology_summary(label: pd.Series) -> dict:
    manifest = primitive_score_manifest()
    primitive_score_name = str(label.get("primitive_score_name", ""))
    return {
        "primitive_score_name": primitive_score_name,
        "score_status": manifest["score_status"],
        "component_features": manifest["components"].get(primitive_score_name, []),
        "weighting_scheme": "equal_weight_mean_screening_composite",
        "caveat": manifest["caveat"],
        "interpretation_order": manifest["interpretation_order"],
        "required_validation_before_mechanistic_use": manifest["required_validation_before_mechanistic_use"],
    }


def null_model_panel_summary(residuals: pd.DataFrame, region_id: str) -> dict:
    subset = residuals[residuals["region_id"].astype(str) == str(region_id)]
    fallback = null_panel_status()
    caveat = "Promotion requires agreement across several appropriate null families; survival of one convenient null is insufficient."
    if subset.empty:
        return {
            "status": fallback["status"],
            "available_null_models": fallback["implemented_null_models"],
            "missing_or_partial_null_models": fallback["missing_or_partial_null_models"],
            "caveat": caveat,
        }
    row = subset.iloc[0]
    available = str(row.get("available_null_models", "") or "")
    missing = str(row.get("missing_or_partial_null_models", "") or "")
    return {
        "status": str(row.get("null_panel_status", "") or fallback["status"]),
        "available_null_models": [item for item in available.split(",") if item]
        or fallback["implemented_null_models"],
        "null_model_count": int(row.get("null_model_count", 0) or 0),
        "null_model_agreement": _finite_float(row.get("null_model_agreement")),
        "null_model_conflict": _bool_field(row, "null_model_conflict"),
        "missing_or_partial_null_models": [item for item in missing.split(",") if item]
        or fallback["missing_or_partial_null_models"],
        "caveat": caveat,
    }


def build_region_card(
    window: pd.Series,
    label: pd.Series,
    residuals: pd.DataFrame,
    features: pd.DataFrame | None = None,
    architecture_candidates: pd.DataFrame | None = None,
) -> dict:
    region_id = str(window["region_id"])
    primitive = str(label.get("primitive_class", "unexplained_dark_anomaly_candidate"))
    assay = recommend_assay(primitive)
    ontology = get_primitive(primitive)
    supporting = (
        str(label.get("top_supporting_features", "")).split(";")
        if pd.notna(label.get("top_supporting_features", ""))
        else []
    )
    conflicting = []
    if str(window.get("artifact_risk_flags", "")):
        conflicting.append("artifact_risk_flags_present")
    if float(label.get("primitive_confidence", 0.0) or 0.0) < 0.4:
        conflicting.append("low_confidence")
    controlled = ""
    subset = residuals[residuals["region_id"].astype(str) == region_id]
    if not subset.empty:
        controlled = str(subset.iloc[0].get("covariates_used", ""))
    top_scores = top_scores_for_region(residuals, region_id)
    null_summary = matched_null_summary(residuals, region_id)
    resid_summary = residualization_summary(residuals, region_id)
    supporting_features = [f for f in supporting if f]
    artifact_risk_flags = str(window.get("artifact_risk_flags", ""))
    observed_feature_evidence = {
        "supporting_features": supporting_features,
        "conflicting_features": conflicting,
        "top_scores": top_scores,
        "matched_null_summary": null_summary,
        "residualization_summary": resid_summary,
        "artifact_risk_flags": artifact_risk_flags,
    }
    primitive_hypothesis = {
        "candidate_label": primitive,
        "hypothesis_statement": ontology.prompt1_interpretation,
        "hypothesis_rationale": ontology.description,
        "required_validation_data": ontology.required_input_level,
        "requires_dynamic_validation": bool(ontology.requires_dynamic_data),
        "recommended_assay": assay.get("assay"),
        "forbidden_interpretation": ontology.forbidden_interpretation,
    }
    mechanistic_bridge = assay.get("mechanistic_bridge", {})
    feature_row = None
    if features is not None and not features.empty and "region_id" in features.columns:
        feature_subset = features[features["region_id"].astype(str) == region_id]
        if not feature_subset.empty:
            feature_row = feature_subset.iloc[0]
    window_shift_status = str(feature_row.get("DFA_window_shift_status", "")) if feature_row is not None else ""
    negative_evidence = candidate_negative_evidence(
        artifact_risk_flags=artifact_risk_flags,
        null_panel_status=str(null_model_panel_summary(residuals, region_id).get("status", "")),
        survives_severe_nulls=(
            bool(label.get("survives_severe_null_panel"))
            if pd.notna(label.get("survives_severe_null_panel", np.nan))
            else None
        ),
        window_shift_status=window_shift_status,
        bridge_status=str(mechanistic_bridge.get("bridge_status", "")),
    )
    architecture_row = None
    if architecture_candidates is not None and not architecture_candidates.empty:
        direct = architecture_candidates.loc[architecture_candidates["region_id"].astype(str) == region_id]
        if direct.empty and "source_representative_region_id" in architecture_candidates.columns:
            direct = architecture_candidates.loc[
                architecture_candidates["source_representative_region_id"].astype(str) == region_id
            ]
        if not direct.empty:
            architecture_row = direct.iloc[0]
    analysis_modes = ["Mode_A_sequence_specific"]
    mode_b_summary = {
        "status": "unavailable_not_run",
        "reason": "Mode B was disabled or no matching interval-level result was available.",
    }
    if architecture_row is not None:
        analysis_modes.append("Mode_B_sequence_indifferent")
        mode_b_summary = {
            "status": str(architecture_row.get("mode_b_score_status", "available")),
            "sequence_indifferent_candidate": str(architecture_row.get("sequence_indifferent_candidate", "unresolved_architecture_candidate")),
            "dominant_mode": str(architecture_row.get("dominant_mode", "unresolved")),
            "sequence_identity_sensitivity": _finite_float(architecture_row.get("sequence_identity_sensitivity")),
            "length_sensitivity": _finite_float(architecture_row.get("length_sensitivity")),
            "copy_number_sensitivity": _finite_float(architecture_row.get("copy_number_sensitivity")),
            "spacing_sensitivity": _finite_float(architecture_row.get("spacing_null_zscore")),
            "sequence_indifference_score": _finite_float(architecture_row.get("sequence_indifference_score")),
            "architecture_null_zscore": _finite_float(architecture_row.get("architecture_null_zscore")),
            "architecture_null_status": str(architecture_row.get("architecture_null_status", "")),
            "promotion_status": str(architecture_row.get("promotion_status", "screening_only_not_for_promotion")),
            "evidence_scope": "model_based_perturbation_and_quantity_screen_not_biological_causality",
        }
    card = {
        "region_id": region_id,
        "coordinates": f"{window.get('chrom')}:{int(window.get('start'))}-{int(window.get('end'))}",
        "window_size": int(window.get("window_size", int(window.get("end")) - int(window.get("start")))),
        "parent_child_multiscale_context": {
            "parent_region_id": window.get("parent_region_id"),
            "child_region_ids": window.get("child_region_ids"),
        },
        "primitive_class": primitive,
        "analysis_modes": analysis_modes,
        "sequence_specific_candidate": primitive,
        "sequence_indifferent_architecture": mode_b_summary,
        "sequence_indifferent_candidate": mode_b_summary.get("sequence_indifferent_candidate", "unavailable"),
        "dominant_mode": mode_b_summary.get("dominant_mode", "Mode_A"),
        "sequence_identity_sensitivity": mode_b_summary.get("sequence_identity_sensitivity"),
        "length_sensitivity": mode_b_summary.get("length_sensitivity"),
        "copy_number_sensitivity": mode_b_summary.get("copy_number_sensitivity"),
        "spacing_sensitivity": mode_b_summary.get("spacing_sensitivity"),
        "candidate_only": True,
        "candidate_statement": "This is a sequence-derived candidate hypothesis, not a confirmed biological primitive.",
        "terminology_scope": TERMINOLOGY_SCOPE,
        "assembly_pangenome_context": assembly_pangenome_context(window),
        "score_methodology": score_methodology_summary(label),
        "null_model_panel": null_model_panel_summary(residuals, region_id),
        "negative_evidence": negative_evidence,
        "candidate_support_status": negative_evidence["candidate_status"],
        "candidate_promotion_status": str(
            label.get("candidate_promotion_status", "screening_only_legacy_null_metadata_unavailable")
        ),
        "survives_severe_null_panel": (
            bool(label.get("survives_severe_null_panel"))
            if pd.notna(label.get("survives_severe_null_panel", np.nan))
            else None
        ),
        "confirmed_name": ontology.confirmed_name,
        "requires_dynamic_validation": bool(ontology.requires_dynamic_data),
        "required_validation_data": ontology.required_input_level,
        "prompt1_allowed_interpretation": ontology.prompt1_interpretation,
        "observed_feature_evidence": observed_feature_evidence,
        "primitive_hypothesis": primitive_hypothesis,
        "mechanistic_bridge": mechanistic_bridge,
        "assay_validation_scope": mechanistic_bridge.get("assay_scope_if_bridge_missing", ""),
        "causal_validation_hierarchy": CAUSAL_VALIDATION_HIERARCHY,
        "native_context_caveat": (
            "Plasmid or synthetic reporter assays lose genomic position, chromatin, 3D contacts, "
            "replication timing, allele context, nearby TE context, and nuclear compartment."
        ),
        "feature_hypothesis_boundary": (
            "Observed features are measured sequence/statistical evidence; "
            "primitive labels are assay-generating hypotheses, not observed molecular properties."
        ),
        "suggested_prompt2_view": ontology.suggested_prompt2_view,
        "primitive_confidence": float(label.get("primitive_confidence", 0.0) or 0.0),
        "top_scores": top_scores,
        "top_supporting_features": supporting_features,
        "conflicting_features": conflicting,
        "artifact_risk_flags": artifact_risk_flags,
        "matched_null_summary": null_summary,
        "residualization_summary": resid_summary,
        "classical_covariates_controlled": controlled,
        "nearest_classical_explanation": "Review TE/repeat/GC/mappability/gene-proximity covariates before prioritization.",
        "why_not_classical": (
            "Candidate remains prioritized after classical covariate control and agrees across the configured severe null panel."
            if bool(label.get("survives_severe_null_panel", False))
            else "Candidate remains a screening hypothesis after classical covariate control; severe null support is incomplete or conflicting."
        ),
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
        "recommended_next_step": "Inspect matched-null diagnostics, artifact flags, mechanistic bridge evidence, and native/control perturbation design.",
        "interpretation_caveat": (
            "This card is an assay-generating hypothesis, not a functional annotation. "
            "If the mechanistic bridge is unvalidated, the assay is exploratory rather than direct primitive validation."
        ),
        "evolutionary_interpretation_caveat": (
            "A causal or quantity-dependent effect does not establish that the element originated or is maintained by selection for that effect."
        ),
    }
    if assay.get("temporal_interaction_test"):
        card["temporal_interaction_test"] = assay["temporal_interaction_test"]
    return card


def make_region_cards(
    windows: pd.DataFrame,
    labels: pd.DataFrame,
    residuals: pd.DataFrame,
    features: pd.DataFrame | None = None,
    architecture_candidates: pd.DataFrame | None = None,
    top_n: int | None = None,
    *,
    progress: bool = False,
) -> list[dict]:
    merged = labels.merge(windows, on="region_id", how="left", suffixes=("_label", ""))
    merged = merged.sort_values("primitive_confidence", ascending=False)
    if top_n:
        merged = merged.head(top_n)
    cards = []
    reporter = ProgressReporter("make-region-cards", total=len(merged)) if progress else None
    if reporter:
        reporter.start("building assay cards")
    for idx, row in enumerate(merged.iterrows(), start=1):
        record = row[1]
        card = build_region_card(
            record,
            record,
            residuals,
            features=features,
            architecture_candidates=architecture_candidates,
        )
        cards.append(card)
        if reporter:
            reporter.update(idx, message=str(record.get("region_id", "")))
    if reporter:
        reporter.finish()
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
                "analysis_modes": ";".join(c.get("analysis_modes", [])),
                "sequence_indifferent_candidate": c.get("sequence_indifferent_candidate", "unavailable"),
                "dominant_mode": c.get("dominant_mode", "Mode_A"),
                "primitive_confidence": c["primitive_confidence"],
                "artifact_risk_flags": c["artifact_risk_flags"],
                "feature_hypothesis_boundary": c["feature_hypothesis_boundary"],
                "assembly_pangenome_caveat": c["assembly_pangenome_context"]["caveat"],
                "score_methodology_caveat": c["score_methodology"]["caveat"],
                "null_model_panel_status": c["null_model_panel"]["status"],
                "candidate_support_status": c["candidate_support_status"],
                "negative_evidence_decision": c["negative_evidence"]["decision"],
                "mechanistic_bridge_status": c["mechanistic_bridge"].get("bridge_status", ""),
                "assay_validation_scope": c["assay_validation_scope"],
                "native_context_caveat": c["native_context_caveat"],
                "recommended_primitive_assay": c["recommended_primitive_assay"],
            }
            for c in cards
        ]
    )
    flat.to_csv(table_path, sep="\t", index=False)
    return {"json": json_path, "tsv": table_path}
