"""Genome browser BED/bedGraph outputs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from darkdna.io.bed import write_bed, write_bedgraph
from darkdna.utils.progress import ProgressReporter, progress_message
from .locus_candidates import write_candidate_locus_outputs


PRIMITIVE_TRACKS = {
    "fractal_scaffold_candidate": "fractal_scaffold_candidates.bed",
    "constraint_grammar_region_candidate": "constraint_grammar_candidates.bed",
    "quantum_susceptible_domain_candidate": "physical_susceptibility_candidates.bed",
    "replication_instability_candidate": "replication_instability_candidates.bed",
    "decoherence_boundary_candidate": "decoherence_boundary_candidates.bed",
    "resonant_pulse_decoder_candidate": "resonant_pulse_decoder_candidates.bed",
    "hysteresis_candidate": "hysteresis_candidates.bed",
    "negative_space_element_candidate": "negative_space_candidates.bed",
    "sequence_regime_boundary_candidate": "sequence_regime_boundary_candidates.bed",
    "TE_grammar_node_candidate": "TE_grammar_node_candidates.bed",
}


def make_tracks(
    windows: pd.DataFrame,
    labels: pd.DataFrame,
    residuals: pd.DataFrame,
    outdir: str | Path,
    *,
    progress: bool = False,
) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    if progress:
        progress_message("make-tracks", f"writing tracks to {out}")
    coords = windows[["region_id", "chrom", "start", "end", "artifact_risk_flags"]].copy()
    global_residual = residuals[residuals["primitive"] == "unexplained_dark_anomaly_candidate_score"][["region_id", "residual_zscore"]]
    scored = coords.merge(global_residual, on="region_id", how="left")
    scored["residual_zscore"] = scored["residual_zscore"].fillna(0.0)
    paths: dict[str, Path] = {}
    bedgraph_path = out / "all_residual_scores.bedGraph"
    write_bedgraph(scored.rename(columns={"residual_zscore": "score"}), bedgraph_path, "score")
    paths["all_residual_scores"] = bedgraph_path

    label_coords = coords.merge(labels, on="region_id", how="inner")
    primitive_bed = out / "primitive_labels.bed"
    label_coords["name"] = label_coords["primitive_class"] + "|" + label_coords["region_id"].astype(str)
    write_bed(label_coords, primitive_bed, columns=["chrom", "start", "end", "name"])
    paths["primitive_labels"] = primitive_bed

    locus_paths = write_candidate_locus_outputs(windows, labels, residuals, out)
    paths["candidate_loci"] = locus_paths["candidate_loci_bed"]

    reporter = ProgressReporter("make-tracks", total=len(PRIMITIVE_TRACKS)) if progress else None
    if reporter:
        reporter.start("writing primitive BED tracks")
    for idx, (primitive, filename) in enumerate(PRIMITIVE_TRACKS.items(), start=1):
        subset = label_coords[label_coords["primitive_class"] == primitive].copy()
        path = out / filename
        if subset.empty:
            path.write_text("", encoding="utf-8")
        else:
            write_bed(subset, path, columns=["chrom", "start", "end", "name"])
        paths[primitive] = path
        if reporter:
            reporter.update(idx, message=primitive)
    if reporter:
        reporter.finish()

    artifact_path = out / "artifact_risk_flags.bed"
    artifact = coords[coords["artifact_risk_flags"].fillna("") != ""].copy()
    artifact["name"] = artifact["artifact_risk_flags"]
    if artifact.empty:
        artifact_path.write_text("", encoding="utf-8")
    else:
        write_bed(artifact, artifact_path, columns=["chrom", "start", "end", "name"])
    paths["artifact_risk_flags"] = artifact_path
    return paths
