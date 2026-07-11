"""HTML report generation."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from jinja2 import Template

from .locus_candidates import block_bootstrap_locus_summary, merge_candidate_loci
from .plots import write_basic_plots


REPORT_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{{ project_name }} - DarkDNA Observer Report</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 32px; color: #1d2528; background: #fbfbf8; }
    h1, h2 { color: #223c3b; }
    table { border-collapse: collapse; width: 100%; margin: 12px 0 24px; background: white; }
    th, td { border: 1px solid #d7ded8; padding: 6px 8px; font-size: 13px; vertical-align: top; }
    th { background: #e7efea; text-align: left; }
    .caveat { padding: 12px; background: #fff5d9; border-left: 4px solid #c6912f; }
    .card { border: 1px solid #d7ded8; border-radius: 6px; padding: 12px; margin: 12px 0; background: white; }
    img { max-width: 100%; height: auto; display: block; margin: 12px 0 20px; }
  </style>
</head>
<body>
  <h1>{{ project_name }}</h1>
  <p class="caveat">This report is hypothesis-generating. Sequence-only labels are candidate hypotheses, not confirmed biological primitives.</p>

  <h2>Project Summary</h2>
  <p>Input windows: {{ window_count }}. Primitive labels: {{ primitive_count }}. Residual rows: {{ residual_count }}.</p>

  <h2>Input Summary</h2>
  <table><tbody>
    {% for key, value in input_summary.items() %}<tr><th>{{ key }}</th><td>{{ value }}</td></tr>{% endfor %}
  </tbody></table>

  <h2>Run Provenance Summary</h2>
  <pre>{{ provenance }}</pre>

  <h2>Window Counts</h2>
  {{ window_counts_html }}

  <h2>Primitive Counts</h2>
  {{ primitive_counts_html }}

  <h2>Artifact-Risk Summary</h2>
  {{ artifact_summary_html }}

  <h2>Locus-Level Candidate Evidence</h2>
  <p class="caveat">Overlapping multiscale windows are dependent observations. Candidate counts should be read from merged loci, not from raw window rows.</p>
  {{ locus_candidates_html }}

  <h2>Block Bootstrap Summary</h2>
  {{ block_bootstrap_html }}

  <h2>Top Candidate Windows</h2>
  {{ top_candidates_html }}

  <h2>Top Residual Anomalies</h2>
  {{ top_residuals_html }}

  <h2>Top Negative-Space Candidates</h2>
  {{ negative_space_html }}

  <h2>Top Boundary-Condition Candidates</h2>
  {{ boundary_html }}

  <h2>Top TE-Grammar Candidates</h2>
  {{ te_html }}

  <h2>Primitive-Specific Plots</h2>
  {% for _, path in plots.items() %}<img src="{{ path.name }}" alt="{{ path.name }}">{% endfor %}

  <h2>Top Region Cards</h2>
  {% for card in cards %}
    <div class="card">
      <h3>{{ card.region_id }} - {{ card.primitive_class }}</h3>
      <p><strong>Candidate-only:</strong> {{ card.candidate_only }}. {{ card.candidate_statement }}</p>
      {% if card.observed_feature_evidence %}
      <p><strong>Observed feature evidence:</strong> {{ card.observed_feature_evidence.supporting_features|join(", ") or "No dominant supporting feature recorded" }}</p>
      {% endif %}
      {% if card.primitive_hypothesis %}
      <p><strong>Primitive hypothesis:</strong> {{ card.primitive_hypothesis.hypothesis_statement }}</p>
      {% endif %}
      {% if card.terminology_scope %}
      <p><strong>Terminology:</strong> {{ card.terminology_scope.dark_operational_use }}</p>
      {% endif %}
      {% if card.assembly_pangenome_context %}
      <p><strong>Assembly context:</strong> {{ card.assembly_pangenome_context.caveat }}</p>
      {% endif %}
      {% if card.score_methodology %}
      <p><strong>Score methodology:</strong> {{ card.score_methodology.score_status }}. {{ card.score_methodology.caveat }}</p>
      {% endif %}
      {% if card.null_model_panel %}
      <p><strong>Null model panel:</strong> {{ card.null_model_panel.status }}. {{ card.null_model_panel.caveat }}</p>
      {% endif %}
      {% if card.mechanistic_bridge %}
      <p><strong>Mechanistic bridge:</strong> {{ card.mechanistic_bridge.measured_feature }} &rarr; {{ card.mechanistic_bridge.proposed_dynamic_phenotype }}</p>
      <p><strong>Bridge status:</strong> {{ card.mechanistic_bridge.bridge_status }}. {{ card.mechanistic_bridge.assay_scope_if_bridge_missing }}</p>
      {% endif %}
      {% if card.native_context_caveat %}
      <p><strong>Native-context caveat:</strong> {{ card.native_context_caveat }}</p>
      {% endif %}
      <p><strong>Coordinates:</strong> {{ card.coordinates }}; <strong>Confidence:</strong> {{ "%.3f"|format(card.primitive_confidence) }}</p>
      <p><strong>Assay:</strong> {{ card.recommended_primitive_assay }}</p>
      <p><strong>Key test:</strong> {{ card.key_interaction_test }}</p>
      <p><strong>Caveat:</strong> {{ card.interpretation_caveat }}</p>
    </div>
  {% endfor %}
</body>
</html>
"""


def _table_html(df: pd.DataFrame, max_rows: int = 20, empty_message: str = "No rows.") -> str:
    if df is None or df.empty:
        return f"<p>{empty_message}</p>"
    return df.head(max_rows).to_html(index=False, escape=True)


def _active_evidence_rows(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return df.copy() if df is not None else pd.DataFrame()
    score_cols = [col for col in ["observed_score", "residual_score", "residual_zscore", "matched_null_zscore"] if col in df.columns]
    if not score_cols:
        return df
    active = pd.Series(False, index=df.index)
    for col in score_cols:
        values = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
        active = active | (values.abs() > 1e-12)
    return df.loc[active].copy()


def _report_locus_columns(loci: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "locus_id",
        "primitive_class",
        "chrom",
        "start",
        "end",
        "n_windows",
        "window_sizes",
        "scale_validation_status",
        "representative_region_id",
        "support_score",
        "locus_empirical_p_value",
        "global_bh_q_value",
        "primitive_bh_q_value",
        "block_id",
    ]
    return loci[[col for col in columns if col in loci.columns]].copy()


def generate_html_report(
    windows: pd.DataFrame,
    labels: pd.DataFrame,
    residuals: pd.DataFrame,
    cards: list[dict],
    outdir: str | Path,
    project_name: str = "darkdna-observer",
    input_summary: dict | None = None,
) -> Path:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    plots = write_basic_plots(residuals, out, labels=labels, windows=windows)
    window_counts = windows.groupby("window_size").size().reset_index(name="count") if "window_size" in windows.columns else pd.DataFrame()
    primitive_counts = labels.groupby("primitive_class").size().reset_index(name="count") if "primitive_class" in labels.columns else pd.DataFrame()
    artifact_summary = (
        windows["artifact_risk_flags"].fillna("").str.get_dummies(sep=";").sum().reset_index().rename(columns={"index": "flag", 0: "count"})
        if "artifact_risk_flags" in windows.columns
        else pd.DataFrame()
    )
    loci = merge_candidate_loci(windows, labels, residuals)
    block_bootstrap = block_bootstrap_locus_summary(loci)
    top_candidates = labels.sort_values("primitive_confidence", ascending=False).head(20)
    top_residuals = _active_evidence_rows(residuals).sort_values("residual_zscore", ascending=False).head(20)
    negative = _active_evidence_rows(residuals[residuals["primitive"] == "negative_space_element_candidate_score"]).sort_values("residual_zscore", ascending=False).head(20)
    boundary = _active_evidence_rows(residuals[residuals["primitive"] == "sequence_regime_boundary_candidate_score"]).sort_values("residual_zscore", ascending=False).head(20)
    te = _active_evidence_rows(residuals[residuals["primitive"] == "TE_grammar_node_candidate_score"]).sort_values("residual_zscore", ascending=False).head(20)
    provenance_path = out / "run_metadata.json"
    provenance = provenance_path.read_text(encoding="utf-8") if provenance_path.exists() else "Provenance not found in report directory."
    html = Template(REPORT_TEMPLATE).render(
        project_name=project_name,
        input_summary=input_summary or {},
        provenance=provenance,
        window_count=len(windows),
        primitive_count=len(labels),
        residual_count=len(residuals),
        window_counts_html=_table_html(window_counts),
        primitive_counts_html=_table_html(primitive_counts),
        artifact_summary_html=_table_html(artifact_summary),
        locus_candidates_html=_table_html(
            _report_locus_columns(loci),
            empty_message="No merged candidate loci passed the candidate filters.",
        ),
        block_bootstrap_html=_table_html(
            block_bootstrap,
            empty_message="No candidate loci available for block bootstrap summary.",
        ),
        top_candidates_html=_table_html(top_candidates),
        top_residuals_html=_table_html(top_residuals),
        negative_space_html=_table_html(negative),
        boundary_html=_table_html(boundary),
        te_html=_table_html(te, empty_message="No active TE-grammar candidates. TE annotation may be absent, or all TE grammar scores are zero after classical controls and matched-null comparison."),
        plots=plots,
        cards=cards[:25],
    )
    path = out / "darkdna_report.html"
    path.write_text(html, encoding="utf-8")
    return path
