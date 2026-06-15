from pathlib import Path

import pandas as pd
import pytest

from tests.integration_helpers import assert_pipeline_outputs, run_cli_pipeline, write_temp_config


SCORE_BY_CANDIDATE = {
    "quantum_susceptible_domain_candidate": "quantum_susceptible_domain_candidate_score",
    "resonant_pulse_decoder_candidate": "resonant_pulse_decoder_candidate_score",
    "sequence_regime_boundary_candidate": "sequence_regime_boundary_candidate_score",
    "negative_space_element_candidate": "negative_space_element_candidate_score",
    "TE_grammar_node_candidate": "TE_grammar_node_candidate_score",
}


def overlaps(a_start, a_end, b_start, b_end):
    return min(a_end, b_end) > max(a_start, b_start)


@pytest.mark.integration
def test_toy_prompt1_pipeline_outputs_and_expected_candidates(tmp_path):
    config = write_temp_config("configs/test_toy.yaml", tmp_path, "toy_run")
    outdir = run_cli_pipeline(config)
    assert_pipeline_outputs(outdir)

    windows = pd.read_parquet(outdir / "dark_windows.parquet")
    scores = pd.read_parquet(outdir / "primitive_scores.parquet")
    labels = pd.read_parquet(outdir / "candidate_primitives.parquet")
    expected = pd.read_csv("data/toy/expected_candidates.tsv", sep="\t")
    scored = windows.merge(scores, on=["region_id", "chrom", "start", "end"], how="inner")
    labeled = windows.merge(labels, on="region_id", how="left")

    for candidate_type, score_col in SCORE_BY_CANDIDATE.items():
        row = expected[expected["expected_candidate_type"] == candidate_type].iloc[0]
        mask = (scored["chrom"] == row["chrom"]) & scored.apply(lambda r: overlaps(r["start"], r["end"], row["start"], row["end"]), axis=1)
        overlapping = scored[mask]
        assert not overlapping.empty, candidate_type
        label_mask = (labeled["chrom"] == row["chrom"]) & labeled.apply(lambda r: overlaps(r["start"], r["end"], row["start"], row["end"]), axis=1)
        has_label = candidate_type in set(labeled.loc[label_mask, "primitive_class"].dropna())
        high_score = overlapping[score_col].max() >= scored[score_col].quantile(0.70)
        assert has_label or high_score, candidate_type
