from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
import yaml
from typer.testing import CliRunner

from darkdna.cli import app


CONFIRMED_FORBIDDEN_LABELS = {
    "hysteresis_element",
    "resonant_pulse_decoder",
    "possibility_gate",
    "criticality_tuner",
    "attractor_gate",
    "future_state_biaser",
    "trajectory_constraint",
    "active_inference_prior_region",
    "genomic_reservoir",
    "contextual_operator",
}


def write_temp_config(source: str | Path, tmp_path: Path, output_name: str) -> Path:
    source = Path(source)
    data = yaml.safe_load(source.read_text(encoding="utf-8"))
    data["output_dir"] = str(tmp_path / output_name)
    path = tmp_path / source.name
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return path


def run_cli_pipeline(config: Path) -> Path:
    runner = CliRunner()
    commands = [
        ["make-windows", "--config", str(config)],
        ["extract-sequence-features", "--config", str(config)],
        ["score-primitives", "--config", str(config)],
        ["build-null-models", "--config", str(config)],
        ["residualize", "--config", str(config)],
        ["infer-primitives", "--config", str(config)],
        ["make-region-cards", "--config", str(config)],
        ["report", "--config", str(config)],
        ["make-tracks", "--config", str(config)],
    ]
    data = yaml.safe_load(config.read_text(encoding="utf-8"))
    outdir = Path(data["output_dir"])
    for command in commands:
        result = runner.invoke(app, command)
        assert result.exit_code == 0, result.output
    return outdir


def assert_pipeline_outputs(outdir: Path) -> None:
    required = [
        "dark_windows.bed",
        "dark_windows.parquet",
        "sequence_features.parquet",
        "primitive_scores.parquet",
        "primitive_score_manifest.json",
        "classical_covariates.parquet",
        "null_model_summary.parquet",
        "null_model_registry.json",
        "residual_scores.parquet",
        "candidate_primitives.parquet",
        "candidate_loci.parquet",
        "candidate_loci.tsv",
        "candidate_loci.bed",
        "candidate_loci_block_bootstrap.tsv",
        "region_cards.json",
        "darkdna_report.html",
        "all_residual_scores.bedGraph",
        "primitive_labels.bed",
    ]
    for name in required:
        assert (outdir / name).exists(), name

    windows = pd.read_parquet(outdir / "dark_windows.parquet")
    labels = pd.read_parquet(outdir / "candidate_primitives.parquet")
    loci = pd.read_parquet(outdir / "candidate_loci.parquet")
    residuals = pd.read_parquet(outdir / "residual_scores.parquet")
    nulls = pd.read_parquet(outdir / "null_model_summary.parquet")
    cards = json.loads((outdir / "region_cards.json").read_text(encoding="utf-8"))

    assert not windows.empty
    assert windows["region_id"].is_unique
    assert (windows["start"] >= 0).all()
    assert (windows["end"] > windows["start"]).all()
    assert (windows["window_size"] > 0).all()
    assert "artifact_risk_flags" in windows.columns
    assert "matched_null_zscore" in residuals.columns
    assert {"null_model_id", "region_id", "null_zscore"}.issubset(nulls.columns)
    assert "primitive_class" in labels.columns
    assert {
        "locus_id",
        "primitive_class",
        "n_windows",
        "locus_empirical_p_value",
        "global_bh_q_value",
        "block_id",
        "scale_validation_status",
    }.issubset(loci.columns)
    assert CONFIRMED_FORBIDDEN_LABELS.isdisjoint(set(labels["primitive_class"].dropna()))
    for label in labels["primitive_class"].dropna().unique():
        assert label == "no_call" or label.endswith("_candidate")
    assert cards
    assert all(card.get("candidate_only") is True for card in cards)
    assert all(card.get("forbidden_interpretation") for card in cards)
    assert all(card.get("observed_feature_evidence") for card in cards)
    assert all(card.get("primitive_hypothesis") for card in cards)
    assert all(card.get("mechanistic_bridge") for card in cards)
    assert all("assay_scope_if_bridge_missing" in card.get("mechanistic_bridge", {}) for card in cards)
    assert all(card.get("causal_validation_hierarchy") for card in cards)
    assert all(card.get("native_context_caveat") for card in cards)
    assert all(card.get("feature_hypothesis_boundary") for card in cards)
    assert all(card.get("terminology_scope") for card in cards)
    assert all(card.get("assembly_pangenome_context") for card in cards)
    assert all(card.get("score_methodology") for card in cards)
    assert all(card.get("null_model_panel") for card in cards)
    assert all(card.get("required_validation_data") for card in cards)
    assert all(card.get("suggested_prompt2_view") for card in cards)
