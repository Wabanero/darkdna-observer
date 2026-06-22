<p align="center">
  <img src="assets/darkdna-observer-logo.png" alt="DarkDNA-Observer logo" width="360">
</p>

# DarkDNA-Observer

DarkDNA-Observer is a sequence-first, view-first scientific Python MVP for discovering candidate dark/noncoding DNA windows whose intrinsic sequence architecture looks unusual after classical covariate controls and matched null models.

The package is designed to produce assay-generating hypotheses, not functional annotations. The main object is a genomic window, not a gene.

## Scientific Scope

DarkDNA-Observer can detect dark/noncoding windows with residual sequence, scale, entropy, grammar, boundary, physical-susceptibility, non-B-DNA, negative-space, or TE-grammar proxy anomalies that are not trivially explained by available classical genomic covariates.

It does not prove:

- teleology
- mystical function
- holographic genome physics
- true future-state bias
- actual quantum DNA effects
- experimentally confirmed non-B-DNA structures
- confirmed temporal or dynamical primitive behavior from sequence alone

Primitive names are operational candidate labels. All high-priority candidates should be interpreted together with matched null models, residualization summaries, and artifact-risk flags.

## Candidate Proxies vs Confirmed Dynamic Primitives

It can detect candidate substrates with intrinsic sequence architectures compatible with memory, rhythm, boundary, physical susceptibility, negative-space, TE grammar, or scale-constraint hypotheses.

It cannot infer true future-state bias, hysteresis, active inference, transition thresholds, possibility gates, or teleological behavior from static sequence alone. Those require dynamic data: time-course, perturbation, pseudotime, state-transition graphs, single-cell data, recovery experiments, live-locus motion, or dose gradients.

Candidate labels:

- `fractal_scaffold_candidate`
- `constraint_grammar_region_candidate`
- `quantum_susceptible_domain_candidate`
- `replication_instability_candidate`
- `chromatin_motion_oscillator_candidate`
- `decoherence_boundary_candidate`
- `hysteresis_candidate`
- `resonant_pulse_decoder_candidate`
- `possibility_gate_candidate`
- `criticality_tuner_candidate`
- `negative_space_element_candidate`
- `sequence_regime_boundary_candidate`
- `TE_grammar_node_candidate`
- `unexplained_dark_anomaly_candidate`

Confirmed dynamic names are reserved for later workflows with appropriate validation data.

## Pipeline

```text
genome + annotation
-> dark window generation
-> intrinsic sequence feature extraction
-> primitive/view candidate scoring
-> classical covariate annotation
-> matched null models
-> residual anomaly scoring
-> candidate primitive assignment
-> artifact-risk annotation
-> region cards
-> primitive-specific assay recommendation
-> HTML report and genome browser tracks
```

This is not a gene-centric pipeline. Protein-coding exons can be excluded, but promoters, introns, TEs, cCREs, enhancers, lncRNAs, and other classical annotations are retained as covariates and labels unless explicitly filtered.

## Installation

DarkDNA-Observer targets Python 3.11 or newer.

```bash
pip install -e .
```

Core dependencies include `typer`, `pandas`, `numpy`, `scipy`, `scikit-learn`, `pyarrow`, `jinja2`, `matplotlib`, `networkx`, and `PyYAML`. Optional dependencies extend feature coverage but are not required for the MVP path.

Missing optional libraries are reported and skipped gracefully. The package should continue running when optional physical-shape, fractal, R/Bioconductor, bigWig, or accelerated FASTA libraries are absent.

## Quick Start

Generate toy data:

```bash
darkdna make-toy-data --outdir demo/toy
```

Generate windows:

```bash
darkdna make-windows \
  --fasta demo/toy/toy.fa \
  --annotation demo/toy/toy.gff3 \
  --te-annotation demo/toy/toy_te.gff3 \
  --blacklist demo/toy/toy_blacklist.bed \
  --mappability demo/toy/toy_mappability.bedGraph \
  --outdir demo/run
```

Run the sequence-first scoring pipeline:

```bash
darkdna extract-sequence-features \
  --fasta demo/toy/toy.fa \
  --windows demo/run/dark_windows.parquet \
  --te-annotation demo/toy/toy_te.gff3 \
  --outdir demo/run

darkdna score-primitives \
  --features demo/run/sequence_features.parquet \
  --outdir demo/run

darkdna build-null-models \
  --scores demo/run/primitive_scores.parquet \
  --features demo/run/sequence_features.parquet \
  --outdir demo/run

darkdna residualize \
  --scores demo/run/primitive_scores.parquet \
  --covariates demo/run/sequence_features.parquet \
  --nulls demo/run/matched_nulls.parquet \
  --outdir demo/run

darkdna infer-primitives \
  --residuals demo/run/residual_scores.parquet \
  --features demo/run/sequence_features.parquet \
  --windows demo/run/dark_windows.parquet \
  --outdir demo/run
```

Create cards, report, and tracks:

```bash
darkdna make-region-cards \
  --windows demo/run/dark_windows.parquet \
  --labels demo/run/primitive_labels.parquet \
  --residuals demo/run/residual_scores.parquet \
  --features demo/run/sequence_features.parquet \
  --outdir demo/run

darkdna report \
  --windows demo/run/dark_windows.parquet \
  --labels demo/run/primitive_labels.parquet \
  --residuals demo/run/residual_scores.parquet \
  --cards demo/run/region_cards.json \
  --outdir demo/run

darkdna make-tracks \
  --windows demo/run/dark_windows.parquet \
  --labels demo/run/primitive_labels.parquet \
  --residuals demo/run/residual_scores.parquet \
  --outdir demo/run
```

## CLI Commands

- `darkdna init`: write example configuration files.
- `darkdna make-toy-data`: create a deterministic toy genome, annotations, tracks, expected anomalies, and toy windows.
- `darkdna make-windows`: generate multiscale dark/noncoding genomic windows.
- `darkdna extract-sequence-features`: compute intrinsic sequence, grammar, boundary, negative-space, non-B-DNA, asymmetry, physical-shape proxy, TE grammar, and scale/fractal features.
- `darkdna score-primitives`: combine features into Prompt 1 candidate scores.
- `darkdna build-null-models`: build matched null summaries.
- `darkdna residualize`: estimate residual anomaly after classical covariate controls.
- `darkdna infer-primitives`: assign candidate labels from residual and matched-null evidence.
- `darkdna make-region-cards`: create JSON/TSV candidate region cards and assay blueprints.
- `darkdna report`: generate an HTML report and diagnostic plots.
- `darkdna make-tracks`: generate BED/bedGraph files for genome browsers.

## Key Outputs

Window generation writes:

- `dark_windows.bed`
- `dark_windows.parquet`

Feature and score commands write:

- `sequence_features.parquet`
- `primitive_scores.parquet`
- `matched_nulls.parquet`
- `residual_scores.parquet`
- `residualization_summary.json`
- `primitive_labels.parquet`

Reporting writes:

- `region_cards.json`
- `region_cards.tsv`
- `darkdna_report.html`
- residual diagnostic plots
- genome browser tracks

Genome browser tracks include:

- `all_residual_scores.bedGraph`
- `primitive_labels.bed`
- `fractal_scaffold_candidates.bed`
- `constraint_grammar_candidates.bed`
- `physical_susceptibility_candidates.bed`
- `replication_instability_candidates.bed`
- `decoherence_boundary_candidates.bed`
- `resonant_pulse_decoder_candidates.bed`
- `hysteresis_candidates.bed`
- `negative_space_candidates.bed`
- `sequence_regime_boundary_candidates.bed`
- `TE_grammar_node_candidates.bed`
- `artifact_risk_flags.bed`

## Region Cards

Every region card states:

```text
This is a sequence-derived candidate hypothesis, not a confirmed biological primitive.
```

Cards include:

- candidate-only status
- coordinates and multiscale context
- candidate primitive label and confidence
- top scores and supporting features
- conflicting features and artifact-risk flags
- matched-null and residualization summaries
- classical covariates controlled
- allowed and forbidden interpretation
- required validation data
- suggested view
- recommended primitive assay
- recommended classical validation assay
- control sequence design
- treatment or perturbation design
- expected positive and negative results

All assay blueprints include the required sequence-by-treatment interaction:

```text
effect = (Native_treatment - Native_control) - (ControlSequence_treatment - ControlSequence_control)
```

Temporal or memory-like candidates require sequence-by-treatment-by-time/history validation before any confirmed dynamic interpretation.

## Supported Inputs

Required:

- genome FASTA or chrom sizes

Optional:

- GFF3/GTF/BED gene annotation
- blacklist BED
- TE BED/GFF3
- cCRE/enhancer/promoter BED
- mappability bigWig, bedGraph, or BED
- assembly gaps BED
- segmental duplication BED
- centromere/telomere BED

The MVP supports scaffold and contig genomes, plant and non-model genomes, missing gene names, missing transcript biotypes, and non-GENCODE GFF3 attributes. It does not assume human chromosome names such as `chr1`.

## Development

Run tests with:

```bash
pytest
```

In this workspace, the available Python executable is:

```powershell
C:\Users\User\miniconda3\envs\genetichyper\python.exe -m pytest -q
```

## License

This project is **source-available**, but it is **not open source**.

Copyright © 2026 Filippo Bergeretti. All rights reserved.

The code, documentation, examples, assets, logos, genomic feature definitions, scoring concepts, sequence-analysis logic, and associated research ideas are publicly visible for evaluation and portfolio-review purposes only.

You may not copy, modify, redistribute, repackage, publish, sublicense, use commercially, or create derivative works from this project without explicit written permission.

Scientific note: this project treats "quantum-susceptible" sequence regions conservatively as physical-susceptibility proxies based on sequence composition and context. It does not claim to demonstrate actual quantum effects in genomic sequences.

See [`LICENSE`](./LICENSE) for the full terms.

## License

MIT
