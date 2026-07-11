"""HTML report generation."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from jinja2 import Template

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

  <h2>Top Candidate Regions</h2>
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
      <p><strong>Coordinates:</strong> {{ card.coordinates }}; <strong>Confidence:</strong> {{ "%.3f"|format(card.primitive_confidence) }}</p>
      <p><strong>Assay:</strong> {{ card.recommended_primitive_assay }}</p>
      <p><strong>Key test:</strong> {{ card.key_interaction_test }}</p>
      <p><strong>Caveat:</strong> {{ card.interpretation_caveat }}</p>
    </div>
  {% endfor %}
</body>
</html>
"""


def _table_html(df: pd.DataFrame, max_rows: int = 20) -> str:
    if df is None or df.empty:
        return "<p>No rows.</p>"
    return df.head(max_rows).to_html(index=False, escape=True)


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
    top_candidates = labels.sort_values("primitive_confidence", ascending=False).head(20)
    top_residuals = residuals.sort_values("residual_zscore", ascending=False).head(20)
    negative = residuals[residuals["primitive"] == "negative_space_element_candidate_score"].sort_values("residual_zscore", ascending=False).head(20)
    boundary = residuals[residuals["primitive"] == "sequence_regime_boundary_candidate_score"].sort_values("residual_zscore", ascending=False).head(20)
    te = residuals[residuals["primitive"] == "TE_grammar_node_candidate_score"].sort_values("residual_zscore", ascending=False).head(20)
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
        top_candidates_html=_table_html(top_candidates),
        top_residuals_html=_table_html(top_residuals),
        negative_space_html=_table_html(negative),
        boundary_html=_table_html(boundary),
        te_html=_table_html(te),
        plots=plots,
        cards=cards[:25],
    )
    path = out / "darkdna_report.html"
    path.write_text(html, encoding="utf-8")
    return path
