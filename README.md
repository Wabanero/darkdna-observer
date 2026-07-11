<p align="center">
  <img src="assets/darkdna-observer-logo.png" alt="DarkDNA-Observer logo" width="360">
</p>

# DarkDNA-Observer

DarkDNA-Observer is a sequence-first analysis toolkit for finding unusual
dark/noncoding DNA windows. It turns a genome plus optional annotations into
candidate regions, primitive scores, residual anomaly scores, validation cards,
HTML reports, and genome-browser tracks.

The project is not a gene annotation pipeline and it is not a functional-claim
engine. Its job is narrower: rank genomic windows whose intrinsic sequence
architecture is unusual after classical covariates, matched null models, and
artifact flags have been considered.

## What It Does

DarkDNA-Observer starts from genomic windows rather than genes. It can exclude
protein-coding exons, while retaining promoters, introns, transposable
elements, enhancers, cCREs, low-mappability regions, and other annotations as
covariates or risk labels.

For each window it computes sequence-derived views:

- base composition, CpG, entropy, k-mer spectra, compression, numeric walks
- repeat, palindrome, homopolymer, periodicity, and fork-texture proxies
- grammar features such as transition surprise, spacing periodicity, recurrent
  k-mers, forbidden-word depletion, and long-range dependency proxies
- left/right asymmetry and sequence-regime boundary features
- non-B-DNA susceptibility proxies such as G4, Z-DNA, triplex/H-DNA,
  cruciform, R-loop, and A-tract signals
- negative-space features that measure structured absence, deserts, and local
  feature voids
- transposable-element grammar features when TE annotations are available
- multiscale/fractal-like summaries and scale persistence

The outputs are candidate hypotheses for prioritization and assay design. A
high score means "this window deserves controlled follow-up", not "this window
has proven function".

## Primitive Properties

Primitive labels are operational classes. They group windows by the type of
sequence evidence that dominates after residualization and matched-null review.

| Primitive candidate | Main evidence | What it means in this tool | Required follow-up |
| --- | --- | --- | --- |
| `fractal_scaffold_candidate` | fractal score, scale persistence, compression anomaly | Scale-persistent sequence texture or numeric-walk structure. | Folding, compaction, or multiscale perturbation assays with shuffled controls. |
| `constraint_grammar_region_candidate` | grammar entropy, forbidden-word depletion, motif-like recurrence, Markov anomaly | Sequence grammar or spacing structure that is not explained by simple composition. | Grammar scramble, tiling, motif/spacing perturbation assays. |
| `quantum_susceptible_domain_candidate` | G4, non-B-DNA aggregate, charge/oxidation susceptibility proxies | Physical-susceptibility candidate based on G-rich and non-B-prone sequence contexts. | G4/oxidation/nanopore or other physical validation. No claim of actual quantum effects. |
| `replication_instability_candidate` | fork texture, repeats, palindromes, skew, non-B propensity | Sequence architecture compatible with replication-stress susceptibility. | Replication timing, fork-pausing, or stress/recovery assays. |
| `chromatin_motion_oscillator_candidate` | spacing autocorrelation, bendability proxy, entropy asymmetry | Static sequence candidate for locus-motion or spatial-dynamics testing. | Live-locus motion or spatial-dynamics assays. |
| `decoherence_boundary_candidate` | entropy cliffs, boundary scores, compression changes, feature voids | Noise/variance boundary candidate from sequence-regime transitions. | Single-cell variance, reporter covariance, or noise-propagation assays. |
| `hysteresis_candidate` | GC asymmetry, nested repeats, non-B propensity, recurrent-kmer orientation | Candidate for history-dependent behavior testing. | Perturbation, recovery, and time/history interaction experiments. |
| `resonant_pulse_decoder_candidate` | 10 bp/147 bp periodicity, Fourier/autocorrelation spacing power | Periodic or phase-like sequence grammar candidate. | Pulse-frequency or time-course perturbation assays. |
| `possibility_gate_candidate` | boundary condition score, negative-space boundary, forbidden-word depletion | Boundary/constraint-like region that may alter reachable state tests. | State-transition, pseudotime, or fate-mapping validation. |
| `criticality_tuner_candidate` | sequence boundary, entropy boundary, compression boundary | Candidate threshold-like sequence transition. | Dose-gradient, perturbation-threshold, or recovery-rate assays. |
| `negative_space_element_candidate` | depleted k-mers, motif deserts, repeat/CpG/G-tract deserts, local voids | Structured absence that remains anomalous after controls. | Rescue/scramble assays that reinsert or disrupt missing tokens. |
| `sequence_regime_boundary_candidate` | left/right regime difference, entropy/GC/CpG/repeat/compression shifts | Candidate boundary between intrinsic sequence regimes. | Boundary smoothing, disruption, insulation, or accessibility assays. |
| `TE_grammar_node_candidate` | TE overlap, TE mosaic score, TE boundary score, TE orientation entropy | TE-derived or TE-mosaic sequence architecture beyond simple overlap. | TE-order/orientation perturbation and TE-derived regulatory comparison. |
| `unexplained_dark_anomaly_candidate` | high residual anomaly without one dominant class | A prioritized dark window that needs review before interpretation. | Matched-null review, artifact checks, and orthogonal validation. |

Every primitive score is paired with robust z-scores, empirical p-values,
matched-null summaries, residual anomaly scores, and artifact-risk flags.

## Interpretation Rules

DarkDNA-Observer deliberately stays conservative.

It can say:

- this window is sequence-derived and candidate-only
- this primitive score is high relative to matched controls
- this residual remains high after available classical covariates
- these features support or conflict with the candidate label
- these assays would test the hypothesis

It cannot say:

- the region has confirmed biological function
- a static sequence proves memory, oscillation, threshold behavior, state bias,
  or any other dynamic mechanism
- a physical-susceptibility proxy proves actual quantum effects
- a grammar-like pattern is a semantic code or a designed program

Confirmed dynamic interpretations require dynamic data such as time-course,
perturbation, recovery, single-cell state transitions, pseudotime, live-locus
motion, or dose-gradient experiments.

## Included Data

The repository includes enough data to run analyses immediately.

| Dataset | Path | Purpose |
| --- | --- | --- |
| Toy synthetic fixture | `data/toy/` | Fast deterministic smoke/integration data with known artificial candidate intervals. |
| Yeast R64 chrI subset | `data/reference/yeast_R64_chrI/` | Small real-genome technical fixture for *Saccharomyces cerevisiae* chromosome I. |
| Arabidopsis TAIR10 chr1 1-2 Mb subset | `data/reference/arabidopsis_TAIR10_chr1_2Mb/` | Small plant/non-model fixture for scaffold/GFF3 compatibility. |

The real fixtures are subset FASTA/GFF3 files plus chrom-size/checksum metadata.
They are intentionally small and versioned. Full source genomes are not kept in
the repository.

## Installation

DarkDNA-Observer targets Python 3.11 or newer.

Recommended Conda setup:

```powershell
cd C:\Users\User\Documents\darkdna-observer
conda env create --file environment.yaml
conda activate darkdna-observer
```

The `environment.yaml` file creates a dedicated `darkdna-observer` environment
from this repository's runtime and development requirements. It installs the
package in editable mode with `pytest` for local validation.

To update an existing environment after dependency changes:

```powershell
conda env update --file environment.yaml --prune
conda activate darkdna-observer
```

If this Windows Conda install raises `CERTIFICATE_VERIFY_FAILED` while solving
the environment, use the one-command workaround instead of changing global Conda
settings:

```powershell
conda env create --insecure --file environment.yaml
```

Without Conda:

```bash
python -m pip install -e ".[dev]"
```

Core dependencies are declared in `pyproject.toml` and mirrored in
`environment.yaml`: `typer`, `pydantic`, `pandas`, `numpy`, `scipy`,
`scikit-learn`, `pyarrow`, `jinja2`, `matplotlib`, `networkx`, and `PyYAML`.

Optional dependencies extend feature coverage for interval operations, FASTA
speedups, bigWig tracks, physical-shape modelling, scale/fractal summaries, and
boosted residual models. Missing optional libraries are skipped gracefully.

## Quick Start With Included Toy Data

Run the full pipeline against the committed toy fixture:

```bash
darkdna run --config configs/test_toy.yaml
```

On Windows PowerShell, to keep a live console log and save the same output:

```powershell
python -m darkdna.cli run --config configs/test_toy.yaml --outdir results/run_progress_toy |
  Tee-Object -FilePath results/run_logs/run_progress_toy.log
```

If you also need to capture native `stderr` without PowerShell formatting it as
`NativeCommandError`, let `cmd` merge the streams before `Tee-Object`:

```powershell
cmd /c "python -m darkdna.cli run --config configs/test_toy.yaml --outdir results/run_progress_toy 2>&1" |
  Tee-Object -FilePath results/run_logs/run_progress_toy.log
```

This is equivalent to running each stage manually:

```bash
darkdna make-windows --config configs/test_toy.yaml
darkdna extract-sequence-features --config configs/test_toy.yaml
darkdna score-primitives --config configs/test_toy.yaml
darkdna build-null-models --config configs/test_toy.yaml
darkdna residualize --config configs/test_toy.yaml
darkdna infer-primitives --config configs/test_toy.yaml
darkdna make-region-cards --config configs/test_toy.yaml
darkdna report --config configs/test_toy.yaml
darkdna make-tracks --config configs/test_toy.yaml
```

Results are written to `results/test_toy/`.

## Quick Start With Included Real Fixtures

Yeast:

```bash
darkdna run --config configs/test_yeast_R64_chrI.yaml
```

Or stage-by-stage:

```bash
darkdna make-windows --config configs/test_yeast_R64_chrI.yaml
darkdna extract-sequence-features --config configs/test_yeast_R64_chrI.yaml
darkdna score-primitives --config configs/test_yeast_R64_chrI.yaml
darkdna build-null-models --config configs/test_yeast_R64_chrI.yaml
darkdna residualize --config configs/test_yeast_R64_chrI.yaml
darkdna infer-primitives --config configs/test_yeast_R64_chrI.yaml
darkdna make-region-cards --config configs/test_yeast_R64_chrI.yaml
darkdna report --config configs/test_yeast_R64_chrI.yaml
darkdna make-tracks --config configs/test_yeast_R64_chrI.yaml
```

Arabidopsis:

```bash
darkdna run --config configs/test_arabidopsis_TAIR10_chr1_2Mb.yaml
```

Or stage-by-stage:

```bash
darkdna make-windows --config configs/test_arabidopsis_TAIR10_chr1_2Mb.yaml
darkdna extract-sequence-features --config configs/test_arabidopsis_TAIR10_chr1_2Mb.yaml
darkdna score-primitives --config configs/test_arabidopsis_TAIR10_chr1_2Mb.yaml
darkdna build-null-models --config configs/test_arabidopsis_TAIR10_chr1_2Mb.yaml
darkdna residualize --config configs/test_arabidopsis_TAIR10_chr1_2Mb.yaml
darkdna infer-primitives --config configs/test_arabidopsis_TAIR10_chr1_2Mb.yaml
darkdna make-region-cards --config configs/test_arabidopsis_TAIR10_chr1_2Mb.yaml
darkdna report --config configs/test_arabidopsis_TAIR10_chr1_2Mb.yaml
darkdna make-tracks --config configs/test_arabidopsis_TAIR10_chr1_2Mb.yaml
```

Results are written under `results/`.

## Regenerating Reference Fixtures

The committed reference fixtures can be regenerated from public Ensembl Genomes
sources:

```bash
python scripts/prepare_test_references.py \
  --dataset yeast_R64_chrI \
  --download \
  --force \
  --out data/reference/yeast_R64_chrI

python scripts/prepare_test_references.py \
  --dataset arabidopsis_TAIR10_chr1_2Mb \
  --download \
  --force \
  --out data/reference/arabidopsis_TAIR10_chr1_2Mb
```

The script downloads temporary full-source files, writes only the configured
subsets, and records checksums.

## CLI Commands

- `darkdna init`: write example configuration files.
- `darkdna run`: execute the complete config-driven pipeline with visible console progress.
- `darkdna make-toy-data`: create deterministic synthetic test data.
- `darkdna make-windows`: generate multiscale dark/noncoding genomic windows.
- `darkdna extract-sequence-features`: compute intrinsic sequence and annotation-derived features.
- `darkdna score-primitives`: combine features into candidate primitive scores.
- `darkdna build-null-models`: build matched null summaries.
- `darkdna residualize`: estimate residual anomaly after classical covariate controls.
- `darkdna infer-primitives`: assign candidate labels from residual and matched-null evidence.
- `darkdna make-region-cards`: create JSON/TSV candidate cards and assay blueprints.
- `darkdna report`: generate an HTML report and diagnostic plots.
- `darkdna make-tracks`: generate BED/bedGraph files for genome browsers.

## Key Outputs

Window generation writes:

- `dark_windows.bed`
- `dark_windows.parquet`

Feature, score, and residual commands write:

- `sequence_features.parquet`
- `primitive_scores.parquet`
- `matched_nulls.parquet`
- `null_model_summary.parquet`
- `residual_scores.parquet`
- `residualization_summary.json`
- `primitive_labels.parquet`
- `candidate_primitives.parquet`

Reporting writes:

- `region_cards.json`
- `region_cards.tsv`
- `darkdna_report.html`
- residual diagnostic plots
- genome-browser BED/bedGraph tracks

## Region Cards

Region cards are designed for review and experiment planning. Each card
contains:

- candidate-only status and interpretation caveat
- coordinates and multiscale context
- primitive class and confidence
- top supporting and conflicting features
- residual and matched-null evidence
- artifact-risk flags
- controlled covariates
- required validation data
- recommended primitive-specific assay
- recommended classical validation assay
- control sequence design
- perturbation design
- expected positive and negative outcomes

The required interaction test for assay blueprints is:

```text
effect = (Native_treatment - Native_control) - (ControlSequence_treatment - ControlSequence_control)
```

Temporal or memory-like candidates require sequence-by-treatment-by-time/history
validation before any dynamic interpretation.

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

The MVP supports scaffold and contig genomes, plant and non-model genomes,
missing gene names, missing transcript biotypes, and non-GENCODE GFF3
attributes. It does not assume human chromosome names such as `chr1`.

## What Is Still Missing

- Better public real-data fixtures with curated TE, mappability, blacklist, gap,
  and centromere/telomere tracks for each organism.
- Continuous-integration jobs that run the real-fixture pipelines, not only the
  fast unit paths.
- Larger validation notebooks or benchmark reports comparing candidate recovery
  across multiple organisms.
- A cleaner public vocabulary for dynamic follow-up views; some internal config
  keys still keep legacy names for compatibility.
- Packaging polish: release metadata, citation file, changelog, and documented
  optional dependency groups per workflow.

## Development

Run tests with:

```powershell
conda activate darkdna-observer
pytest
```

## License

This project is **source-available**, but it is **not open source**.

Copyright (c) 2026 Filippo Bergeretti. All rights reserved.

The code, documentation, examples, assets, logos, genomic feature definitions, scoring concepts, sequence-analysis logic, and associated research ideas are publicly visible for evaluation and portfolio-review purposes only.

You may not copy, modify, redistribute, repackage, publish, sublicense, use commercially, or create derivative works from this project without explicit written permission.

Scientific note: this project treats "quantum-susceptible" sequence regions conservatively as physical-susceptibility proxies based on sequence composition and context. It does not claim to demonstrate actual quantum effects in genomic sequences.

See [`LICENSE`](./LICENSE) for the full terms.

