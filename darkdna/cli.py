"""Command line interface for darkdna-observer."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Optional

import pandas as pd
import typer

from darkdna.features.sequence import extract_features_for_windows, write_sequence_features
from darkdna.features.te_grammar import annotate_te_grammar
from darkdna.io.gff import read_te_annotation
from darkdna.primitives.labeler import assign_primitive_labels, write_primitive_labels
from darkdna.reports.genome_browser_tracks import make_tracks as write_tracks
from darkdna.reports.html_report import generate_html_report
from darkdna.reports.region_cards import make_region_cards, write_region_cards
from darkdna.residuals.null_models import build_matched_null_models, write_matched_nulls
from darkdna.residuals.residual_model import residualize_scores, write_residual_outputs
from darkdna.toy_data import make_toy_data
from darkdna.utils.config import load_config, resolve_config_path, write_config_snapshot
from darkdna.utils.logging import get_logger
from darkdna.utils.provenance import write_provenance
from darkdna.views.multiscale_profiles import compute_multiscale_profiles
from darkdna.views.primitive_scores import score_primitives, write_primitive_scores
from darkdna.windows.make_windows import make_dark_windows, write_windows


app = typer.Typer(help="DarkDNA-Observer sequence-first dark/noncoding DNA hypothesis generator.")
logger = get_logger()


def read_table(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if p.suffix == ".parquet":
        return pd.read_parquet(p)
    if p.suffix in {".tsv", ".bed", ".bedGraph"}:
        return pd.read_csv(p, sep="\t")
    if p.suffix == ".csv":
        return pd.read_csv(p)
    return pd.read_table(p)


def cfg_path(config_path: Path | None, cfg, field: str, explicit: Path | None = None) -> Path | None:
    return explicit or resolve_config_path(config_path, getattr(cfg, field, None))


def cfg_outdir(config_path: Path | None, cfg, explicit: Path | None = None) -> Path:
    resolved = explicit or resolve_config_path(config_path, getattr(cfg, "output_dir", None))
    return resolved or Path("darkdna_run")


def require_configured_path(value: Path | None, label: str) -> Path:
    if value is None:
        raise typer.BadParameter(f"{label} is required unless supplied by --config")
    return value


@app.command()
def init(outdir: Path = typer.Option(Path("."), help="Directory to initialize.")) -> None:
    """Create example config files."""

    outdir.mkdir(parents=True, exist_ok=True)
    cfg = load_config(None)
    config_dir = outdir / "configs"
    config_dir.mkdir(exist_ok=True)
    write_config_snapshot(cfg, config_dir / "example_sequence_first.yaml")
    schema_src = Path(__file__).resolve().parents[1] / "configs" / "schema.yaml"
    if schema_src.exists():
        shutil.copyfile(schema_src, config_dir / "schema.yaml")
    logger.info("Initialized darkdna-observer config files in %s", outdir)


@app.command("make-windows")
def make_windows_cmd(
    fasta: Optional[Path] = typer.Option(None, help="Genome FASTA."),
    chrom_sizes: Optional[Path] = typer.Option(None, help="Chrom sizes file."),
    annotation: Optional[Path] = typer.Option(None, help="GTF/GFF3 gene annotation."),
    blacklist: Optional[Path] = typer.Option(None, help="Blacklist BED."),
    te_annotation: Optional[Path] = typer.Option(None, help="TE BED/GFF3 annotation."),
    ccre: Optional[Path] = typer.Option(None, help="cCRE BED."),
    enhancer: Optional[Path] = typer.Option(None, help="Enhancer BED."),
    promoter: Optional[Path] = typer.Option(None, help="Promoter BED."),
    mappability: Optional[Path] = typer.Option(None, help="Mappability bigWig/bedGraph/BED."),
    assembly_gaps: Optional[Path] = typer.Option(None, help="Assembly gaps BED."),
    segmental_duplication: Optional[Path] = typer.Option(None, help="Segmental duplication BED."),
    centromere_telomere: Optional[Path] = typer.Option(None, help="Centromere/telomere BED."),
    outdir: Optional[Path] = typer.Option(None, help="Output directory."),
    config: Optional[Path] = typer.Option(None, help="YAML config."),
    window_sizes: Optional[str] = typer.Option(None, help="Comma-separated window sizes."),
    step_fraction: Optional[float] = typer.Option(None, help="Step as fraction of window size."),
    exclude_coding_exons: Optional[bool] = typer.Option(None, help="Exclude coding exon-overlapping windows."),
) -> None:
    cfg = load_config(
        config,
        {
            "window_sizes": [int(x) for x in window_sizes.split(",")] if window_sizes else None,
            "step_fraction": step_fraction,
            "exclude_coding_exons": exclude_coding_exons,
        },
    )
    fasta = cfg_path(config, cfg, "fasta", fasta)
    chrom_sizes = cfg_path(config, cfg, "chrom_sizes", chrom_sizes)
    annotation = cfg_path(config, cfg, "annotation", annotation)
    blacklist = cfg_path(config, cfg, "blacklist", blacklist)
    te_annotation = cfg_path(config, cfg, "te_annotation", te_annotation)
    ccre = cfg_path(config, cfg, "ccre", ccre)
    enhancer = cfg_path(config, cfg, "enhancer", enhancer)
    promoter = cfg_path(config, cfg, "promoter", promoter)
    mappability = cfg_path(config, cfg, "mappability", mappability)
    assembly_gaps = cfg_path(config, cfg, "assembly_gaps", assembly_gaps)
    segmental_duplication = cfg_path(config, cfg, "segmental_duplication", segmental_duplication)
    centromere_telomere = cfg_path(config, cfg, "centromere_telomere", centromere_telomere)
    outdir = cfg_outdir(config, cfg, outdir)
    windows = make_dark_windows(
        fasta=fasta,
        chrom_sizes_path=chrom_sizes,
        window_sizes=cfg.window_sizes,
        step_fraction=cfg.step_fraction,
        exclude_coding_exons=cfg.exclude_coding_exons,
        annotation_path=annotation,
        blacklist_path=blacklist,
        te_annotation_path=te_annotation,
        ccre_path=ccre,
        enhancer_path=enhancer,
        promoter_path=promoter,
        mappability_path=mappability,
        assembly_gaps_path=assembly_gaps,
        segmental_duplication_path=segmental_duplication,
        centromere_telomere_path=centromere_telomere,
        promoter_bp=cfg.promoter_bp,
        artifact_thresholds=cfg.artifact_thresholds,
    )
    paths = write_windows(windows, outdir)
    write_provenance(outdir, "darkdna make-windows", cfg, [fasta, chrom_sizes, annotation, blacklist, te_annotation])
    logger.info("Wrote %d windows to %s and %s", len(windows), paths["bed"], paths["parquet"])


@app.command("extract-sequence-features")
def extract_sequence_features_cmd(
    fasta: Optional[Path] = typer.Option(None, help="Genome FASTA."),
    windows: Optional[Path] = typer.Option(None, help="Window table parquet/tsv."),
    outdir: Optional[Path] = typer.Option(None, help="Output directory."),
    te_annotation: Optional[Path] = typer.Option(None, help="Optional TE annotation for TE grammar features."),
    config: Optional[Path] = typer.Option(None, help="YAML config."),
) -> None:
    cfg = load_config(config)
    outdir = cfg_outdir(config, cfg, outdir)
    fasta = require_configured_path(cfg_path(config, cfg, "fasta", fasta), "FASTA")
    windows = windows or outdir / "dark_windows.parquet"
    te_annotation = cfg_path(config, cfg, "te_annotation", te_annotation)
    window_table = read_table(windows)
    features = extract_features_for_windows(window_table, fasta)
    if te_annotation:
        te_features = annotate_te_grammar(window_table, read_te_annotation(te_annotation))
        features = features.merge(te_features, on="region_id", how="left")
    features = compute_multiscale_profiles(features, windows=window_table)
    path = write_sequence_features(features, outdir)
    write_provenance(outdir, "darkdna extract-sequence-features", cfg, [fasta, windows, te_annotation])
    logger.info("Wrote sequence features for %d regions to %s", len(features), path)


@app.command("score-primitives")
def score_primitives_cmd(
    features: Optional[Path] = typer.Option(None, help="Sequence feature table."),
    outdir: Optional[Path] = typer.Option(None, help="Output directory."),
    config: Optional[Path] = typer.Option(None, help="YAML config."),
) -> None:
    cfg = load_config(config)
    outdir = cfg_outdir(config, cfg, outdir)
    features = features or outdir / "sequence_features.parquet"
    feature_table = read_table(features)
    scores = score_primitives(feature_table)
    path = write_primitive_scores(scores, outdir)
    write_provenance(outdir, "darkdna score-primitives", cfg, [features])
    logger.info("Wrote primitive scores to %s", path)


@app.command("build-null-models")
def build_null_models_cmd(
    scores: Optional[Path] = typer.Option(None, help="Primitive score table."),
    features: Optional[Path] = typer.Option(None, help="Feature/covariate table."),
    outdir: Optional[Path] = typer.Option(None, help="Output directory."),
    n_controls: Optional[int] = typer.Option(None, help="Matched controls per region."),
    config: Optional[Path] = typer.Option(None, help="YAML config."),
) -> None:
    cfg = load_config(config)
    outdir = cfg_outdir(config, cfg, outdir)
    scores = scores or outdir / "primitive_scores.parquet"
    features = features or outdir / "sequence_features.parquet"
    n_controls = n_controls if n_controls is not None else cfg.n_null
    score_table = read_table(scores)
    feature_table = read_table(features)
    nulls = build_matched_null_models(score_table, feature_table, n_controls=n_controls)
    path = write_matched_nulls(nulls, outdir)
    write_provenance(outdir, "darkdna build-null-models", cfg, [scores, features])
    logger.info("Wrote matched null summaries to %s", path)


@app.command()
def residualize(
    scores: Optional[Path] = typer.Option(None, help="Primitive score table."),
    covariates: Optional[Path] = typer.Option(None, help="Classical covariate/feature table."),
    nulls: Optional[Path] = typer.Option(None, help="Matched null summary table."),
    outdir: Optional[Path] = typer.Option(None, help="Output directory."),
    method: str = typer.Option("linear", help="linear, robust_linear, random_forest, gradient_boosting, xgboost, or lightgbm."),
    config: Optional[Path] = typer.Option(None, help="YAML config."),
) -> None:
    cfg = load_config(config)
    outdir = cfg_outdir(config, cfg, outdir)
    scores = scores or outdir / "primitive_scores.parquet"
    covariates = covariates or outdir / "sequence_features.parquet"
    nulls = nulls or outdir / "matched_nulls.parquet"
    score_table = read_table(scores)
    covariate_table = read_table(covariates)
    null_table = read_table(nulls) if nulls else pd.DataFrame()
    residuals, summary = residualize_scores(score_table, covariate_table, null_table, method=method)
    paths = write_residual_outputs(residuals, summary, outdir)
    write_provenance(outdir, "darkdna residualize", cfg, [scores, covariates, nulls])
    logger.info("Wrote residual scores to %s", paths["residuals"])


@app.command("infer-primitives")
def infer_primitives_cmd(
    residuals: Optional[Path] = typer.Option(None, help="Residual score table."),
    features: Optional[Path] = typer.Option(None, help="Feature table."),
    windows: Optional[Path] = typer.Option(None, help="Window table."),
    outdir: Optional[Path] = typer.Option(None, help="Output directory."),
    residual_threshold: float = typer.Option(2.0, help="Residual z-score threshold."),
    matched_null_threshold: float = typer.Option(2.0, help="Matched null z-score threshold."),
    config: Optional[Path] = typer.Option(None, help="YAML config."),
) -> None:
    cfg = load_config(config)
    outdir = cfg_outdir(config, cfg, outdir)
    residuals = residuals or outdir / "residual_scores.parquet"
    features = features or outdir / "sequence_features.parquet"
    windows = windows or outdir / "dark_windows.parquet"
    residual_table = read_table(residuals)
    feature_table = read_table(features) if features else pd.DataFrame()
    window_table = read_table(windows) if windows else pd.DataFrame()
    labels = assign_primitive_labels(
        residual_table,
        features=feature_table,
        windows=window_table,
        residual_threshold=residual_threshold,
        matched_null_threshold=matched_null_threshold,
    )
    path = write_primitive_labels(labels, outdir)
    write_provenance(outdir, "darkdna infer-primitives", cfg, [residuals, features, windows])
    logger.info("Wrote primitive labels to %s", path)


@app.command("make-region-cards")
def make_region_cards_cmd(
    windows: Optional[Path] = typer.Option(None, help="Window table."),
    labels: Optional[Path] = typer.Option(None, help="Primitive labels table."),
    residuals: Optional[Path] = typer.Option(None, help="Residual score table."),
    features: Optional[Path] = typer.Option(None, help="Feature table."),
    outdir: Optional[Path] = typer.Option(None, help="Output directory."),
    top_n: Optional[int] = typer.Option(None, help="Limit number of cards."),
    config: Optional[Path] = typer.Option(None, help="YAML config."),
) -> None:
    cfg = load_config(config)
    outdir = cfg_outdir(config, cfg, outdir)
    windows = windows or outdir / "dark_windows.parquet"
    labels = labels or outdir / "primitive_labels.parquet"
    residuals = residuals or outdir / "residual_scores.parquet"
    features = features or outdir / "sequence_features.parquet"
    window_table = read_table(windows)
    label_table = read_table(labels)
    residual_table = read_table(residuals)
    feature_table = read_table(features) if features else pd.DataFrame()
    cards = make_region_cards(window_table, label_table, residual_table, feature_table, top_n=top_n)
    paths = write_region_cards(cards, outdir)
    write_provenance(outdir, "darkdna make-region-cards", cfg, [windows, labels, residuals, features])
    logger.info("Wrote %d region cards to %s", len(cards), paths["json"])


@app.command()
def report(
    windows: Optional[Path] = typer.Option(None, help="Window table."),
    labels: Optional[Path] = typer.Option(None, help="Primitive labels table."),
    residuals: Optional[Path] = typer.Option(None, help="Residual score table."),
    cards: Optional[Path] = typer.Option(None, help="Region cards JSON."),
    outdir: Optional[Path] = typer.Option(None, help="Output directory."),
    project_name: str = typer.Option("darkdna-observer", help="Report project name."),
    config: Optional[Path] = typer.Option(None, help="YAML config."),
) -> None:
    cfg = load_config(config)
    outdir = cfg_outdir(config, cfg, outdir)
    windows = windows or outdir / "dark_windows.parquet"
    labels = labels or outdir / "primitive_labels.parquet"
    residuals = residuals or outdir / "residual_scores.parquet"
    cards = cards or outdir / "region_cards.json"
    project_name = cfg.project_name if config and project_name == "darkdna-observer" else project_name
    window_table = read_table(windows)
    label_table = read_table(labels)
    residual_table = read_table(residuals)
    card_list = json.loads(Path(cards).read_text(encoding="utf-8"))
    path = generate_html_report(window_table, label_table, residual_table, card_list, outdir, project_name=project_name)
    write_provenance(outdir, "darkdna report", cfg, [windows, labels, residuals, cards])
    logger.info("Wrote HTML report to %s", path)


@app.command("make-tracks")
def make_tracks_cmd(
    windows: Optional[Path] = typer.Option(None, help="Window table."),
    labels: Optional[Path] = typer.Option(None, help="Primitive labels table."),
    residuals: Optional[Path] = typer.Option(None, help="Residual score table."),
    outdir: Optional[Path] = typer.Option(None, help="Output directory."),
    config: Optional[Path] = typer.Option(None, help="YAML config."),
) -> None:
    cfg = load_config(config)
    outdir = cfg_outdir(config, cfg, outdir)
    windows = windows or outdir / "dark_windows.parquet"
    labels = labels or outdir / "primitive_labels.parquet"
    residuals = residuals or outdir / "residual_scores.parquet"
    paths = write_tracks(read_table(windows), read_table(labels), read_table(residuals), outdir)
    write_provenance(outdir, "darkdna make-tracks", cfg, [windows, labels, residuals])
    logger.info("Wrote %d genome browser tracks to %s", len(paths), outdir)


@app.command("make-toy-data")
def make_toy_data_cmd(
    outdir: Path = typer.Option(Path("toy_darkdna"), "--out", "--outdir", help="Toy-data output directory."),
    seed: int = typer.Option(42, help="Random seed."),
) -> None:
    paths = make_toy_data(outdir, seed=seed)
    write_provenance(outdir, "darkdna make-toy-data", load_config(None), [])
    logger.info("Wrote toy data to %s", outdir)
    for key, path in paths.items():
        logger.info("%s: %s", key, path)


if __name__ == "__main__":
    app()
