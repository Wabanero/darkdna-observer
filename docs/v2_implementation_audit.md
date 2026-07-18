# DarkDNA-Observer v2 implementation audit

Audit date: 2026-07-18  
Audited baseline commit: `0ee4bcd` (`Add candidate card visuals`)  
Scope: repository state before the v2 Phase 1 corrections in this change set.

## Baseline verification

The fast suite completed with `34 passed, 3 deselected in 17.79s` using Python from the `darkdna-observer` Conda environment. The three deselected tests are the toy, yeast, and Arabidopsis pipeline integrations. An all-tests run reached the large Arabidopsis fixture and timed out after 121.8 seconds; 14 tests had passed and no failure had been reported before the timeout. This is not recorded as a full integration pass.

The worktree already contained a user modification to `configs/test_arabidopsis_TAIR10_chr1_2Mb.yaml`. This implementation does not modify or revert that file.

## Existing architecture and public surface

The repository is a compact Python package with a Typer CLI, Pydantic configuration (plus a dataclass fallback), pandas/Parquet output contracts, sequence feature modules, candidate primitive views, matched genomic controls, blocked residualization, locus merging, region cards, HTML reports, browser tracks, deterministic toy data, and small yeast/Arabidopsis fixtures.

Existing public CLI commands are:

- `init`
- `make-windows`
- `extract-sequence-features`
- `score-primitives`
- `build-null-models`
- `residualize`
- `infer-primitives`
- `make-region-cards`
- `report`
- `make-tracks`
- `run`
- `make-toy-data`

The default `run` path executes nine stages from window generation through browser tracks. Phase 1 retains this sequence-first path and does not add Mode B or the later v2 commands.

Current output contracts include `dark_windows`, `sequence_features`, `primitive_scores`, `matched_nulls`, `classical_covariates`, `residual_scores`, primitive labels, merged candidate loci, block-bootstrap summaries, region cards, HTML/plot files, and browser tracks. Provenance includes config snapshots, checksums, software versions, command logs, and optional-dependency status.

Optional dependencies are grouped for interval tools, FASTA readers, tracks, physical modelling, scale analysis, boosted models, and reports. Missing optional packages generally do not stop the sequence-first workflow.

## What is scientifically adequate

- Dark, noncoding, unannotated, assembly-difficult, and evolutionary-junk terminology is kept separate.
- Primitive labels are described as candidate hypotheses rather than confirmed mechanisms.
- The retired quantum label maps to the conservative physical-susceptibility candidate.
- Static sequence outputs do not include the prohibited Prompt 2 dynamic measurements.
- Windows carry artifact flags for mappability, gaps, segmental duplications, sequence ambiguity, repeat density, low complexity, scaffold edges, and related reference risks.
- Classical covariates are assembled separately from anomaly views, reducing direct target leakage.
- Residual predictions use chromosome or genomic-block holdouts when those groups are available.
- Overlapping discovery windows are merged into locus-level evidence before FDR and counting.
- Region cards separate observed features, candidate hypotheses, mechanistic bridges, allowed interpretations, forbidden interpretations, and factorial assay contrasts.
- The null registry already states that one convenient matched-control family is not a severe null panel.
- The software records missing pangenome, copy-number, presence/absence, and assembly-completeness inputs instead of silently assuming reference completeness.

## MVP heuristics found in the baseline

- The multiscale module used three fixed DFA-like scales, one GC/AT walk, a direct mean of unrelated Fourier/wavelet/DFA values, and several zero or empty placeholder outputs.
- Parent-child score similarity was called scale persistence even though it was only a descriptive cross-window similarity.
- Raw gzip/bz2/lzma ratios retained fixed compressor overhead and were compared without same-length null calibration.
- `strong_weak_H_bond_numeric_walk` and `GC_AT_numeric_walk` used the same mapping under two names.
- DNA shape used a small dinucleotide screen; stiffness was defined as `2 - bendability`; predictor-backed fields stayed unavailable even when rpy2 was detected.
- Non-B-DNA proxies for unlike structures were averaged into `non_B_DNA_aggregate_score`.
- Primitive scores were equal-weight means of heterogeneous, correlated, differently oriented features.
- `unexplained_dark_anomaly_candidate_score` was the mean of the other primitive scores and therefore was not unexplained multivariate outlierness.
- Primitive-score “empirical p-values” were ranks among observed genomic windows, not tails under a named null.
- Residual significance relied on one unconditional residual standard deviation when conditional variance was not modelled.
- `classical_explanation_fraction` was one model-level R2 repeated on every row, despite its per-region-sounding name.
- `primitive_confidence` was an uncalibrated 0-1 prioritization heuristic.

## Misleading names and Phase 1 decisions

| Baseline field | Phase 1 canonical field or decision | Compatibility policy |
| --- | --- | --- |
| `fractal_score` | `multiscale_texture_screening_score` | Deprecated alias retained with an explicit warning; unavailable is `NA` when diagnostics fail. |
| `DFA_like_exponent_estimator` / `Hurst_like_exponent_estimator` | `DFA_exponent` plus scaling interval, R2, residual, CI, surrogate, and shift diagnostics | Deprecated aliases retained. |
| `scale_persistence_score` | `multiscale_parent_child_similarity_screen` | Deprecated alias retained; explicitly not scaling evidence. |
| raw compression ratios | header-corrected ratios and same-length dinucleotide/k-mer null z-scores | Raw ratios remain compatibility diagnostics with a warning. |
| duplicate `GC_AT_numeric_walk` | canonical `strong_weak_H_bond_numeric_walk` | Alias values retained with a deprecation warning; the registry contains one canonical mapping. |
| predictor-like shape placeholders | predictor records with value, version, normalization, status, and warning | Missing validated predictors produce `NA`; legacy dinucleotide screens are separately named. |
| `DNA_stiffness_proxy` | retired | Alias is `NA` with the reason that the direct transformation lacked support. |
| `non_B_DNA_aggregate_score` | structure-specific Level-1 records and fields | Aggregate is `NA` with a deprecation warning. |
| equal-weight primitive score | covariance-aware robust cohort-standardized screen | Equal-weight value retained only as `*_legacy_screening_composite`. |
| observed-rank primitive p-value | no p-value without explicit null | Value is `NA`, status `unavailable`, and reason is human-readable. |
| mean-based unexplained anomaly | block-aware cross-fitted robust multivariate outlierness | No row-local fallback mean is produced. |
| `classical_explanation_fraction` | `classical_model_global_r2` | Deprecated alias retained with a warning. |
| `primitive_confidence` | `primitive_priority` | Deprecated alias retained; status says it is not confidence or probability. |

## What can be extended safely

- The Typer command/stage structure can accept optional later-mode commands without changing the default Mode A run.
- Pydantic config can add nested mode settings with disabled defaults and explicit skip records.
- Parquet tables and JSON manifests can add columns without removing legacy aliases.
- Region cards already have natural extension points for evidence tensors, life history, Mode B, conformation, RNA mechanisms, stopping rules, and integrated evidence.
- The existing blocked CV, locus merging, block bootstrap, and null registry are suitable foundations for the Phase 2 severe-null work.
- Assembly and artifact fields can be extended with pangenome confidence and coordinate portability without changing sequence feature APIs.

## Data required by each v2 mode

| Mode | Required or strongly preferred data | Missing-data behavior |
| --- | --- | --- |
| A — sequence-specific architecture | FASTA; windows; composition, repeat, annotation, mappability, and assembly controls; explicit sequence and genomic nulls | Sequence features can run from FASTA, but strong anomaly claims remain unavailable when diagnostic/null requirements fail. |
| B — sequence-indifferent architecture | interval sources, length/copy-number/presence-absence/spacing data, matched replacements or controlled perturbations, phenotype table for causal analyses | Entire optional stage must be skipped with `status=unavailable`, `value=NA`, and a reason. |
| C — conformation and molecular state | sequence potential plus torsion, transcription, ions, temperature, methylation, nucleosomes, replication, readers/resolvers, and ideally observed-structure assays | Sequence alone yields Level 1 only; Level 2/3 remain unavailable. |
| D — transcription/RNA/translation | CAGE/RAMPAGE or nascent RNA, processed/localized RNA, stability, RBP/structure data, small-RNA evidence, Ribo-seq and proteomics, plus causal contrasts | Correlation never auto-assigns a mechanism; absent tracks skip the relevant evidence layer. |
| E — dynamic perturbation | time/pseudotime, treatment and recovery, dose or pulse design, ordered perturbations, matched controls, replication | Static inputs cannot confirm hysteresis, pulse decoding, criticality, state gating, motion, or noise insulation. |
| F — evolutionary history | synteny, phylogeny, population variation, TE copy age/subfamily/phylogeny, deletion and mutation context, effective population size, and explicit evolutionary nulls | Current role, maintenance, origin, selfishness, deleterious burden, and exaptation remain unresolved when these data are absent. |

## What current inputs cannot establish

The current sequence-first inputs cannot establish selected biological function, organism fitness effects, maintenance selection, adaptive origin, exaptation, selfish origin, deleterious persistence, transcription-process function, RNA-product function, hidden translation, actual non-B structure formation, charge transfer, quantum effects, physical hysteresis, pulse decoding, criticality, future-state bias, live chromatin motion, pangenome representativeness, or sequence-indifferent causality.

Sequence anomaly, entropy, compression, motif presence, predicted shape, conservation, transcription, accessibility, binding, association, and reporter activity are individually insufficient for those claims.

## Phase 2 and Phase 3 implementation

Phase 2 adds a complete severe-null framework: a registry for every required
family, input-dependent availability, block-aware matched calibration,
candidate-level agreement/conflict fields, deterministic sequence transforms, a
reference-conditioned evolutionary-process simulator, and a native-versus-null
benchmark whose result remains conditional on null definition. Foundation models
are local-adapter-only and are never downloaded.

Phase 3 adds first-class Mode B packages for amount, interval length, repeat
burden, copy number, presence/absence, comparative length, spacing, occupancy,
heterochromatic mass, sequence-indifference transformations, block-aware
architecture nulls, candidate classification, Mode A/Mode B comparison, output
manifests, CLI commands, region-card fields, and report sections. Mode B is
disabled by default. Missing optional tracks generate explicit unavailable
records rather than zero placeholders.

## Current phase boundary

The Function Evidence Tensor and Element Life-History model (Phase 4), copy-level
TE and conformation graph (Phase 5), RNA/hidden translation (Phase 6), dynamic
validation (Phase 7), and Onion/comparative multi-genome benchmark (Phase 8)
remain deferred. The default-state native-versus-randomized component originally
listed under the broader comparative roadmap is complete in Phase 2.

## Phase 1 verification

Final verification on the Phase 1 implementation:

- non-integration suite: `47 passed, 3 deselected in 18.24s`;
- toy end-to-end pipeline: passed in `170.2s`;
- yeast R64 chrI end-to-end pipeline: passed in `29.2s`;
- Arabidopsis TAIR10 chr1 2 Mb end-to-end pipeline: passed in `155.8s` using the denser pre-existing user configuration;
- `python -m compileall -q darkdna`: passed.

All tests are offline and no model or dataset download is performed.
