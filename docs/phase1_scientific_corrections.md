# Phase 1 scientific corrections

## Multiscale texture

`compute_scale_fractal_features` now uses configurable log-spaced DFA scales, automatic contiguous scaling-region selection, fit R2 and residual RMSE, slope confidence intervals, composition-preserving surrogate tests, three non-duplicated sequence mappings, multifractal DFA when length permits, wavelet-spectrum summaries when PyWavelets is installed, and nearby-window shift stability. The canonical output is `multiscale_texture_screening_score`.

If a required diagnostic fails, the score is `NA`. The deprecated `fractal_score` alias carries the same value and an explicit warning; it is not a fractal claim.

## Compression

Compression features subtract the empty-stream byte overhead for gzip, bz2, and lzma. Each compressor is calibrated separately against deterministic same-length dinucleotide and k-mer-block surrogate families. Outputs include null means, null standard deviations, z-scores, null counts, seed, window size, calibration status, and multi-compressor directional agreement.

The default feature-extraction path uses two replicates per null family to keep whole-genome screening tractable and therefore marks calibration `partial`; callers must increase `n_surrogates` for inferential use. gzip/bz2 use level 1 and lzma uses preset 0, with those settings recorded in the predictor version.

Positive compression z-scores mean the sequence compressed more than that null family. They do not mean functional information.

## Numeric walks

Canonical mappings are purine/pyrimidine, strong/weak hydrogen bond, and amino/keto. The mapping registry records rationale, table, orientation sensitivity, reverse-complement behavior, and null expectation. `GC_AT_numeric_walk` is a deprecated alias of the strong/weak walk and is no longer a separately defined mapping.

## DNA shape

Every strong DNA-shape field has predictor, predictor version, raw value, normalization, availability status, and warning. No external model is downloaded. Missing DNAshapeR/dnacurve or another validated model produces `NA`.

The historical dinucleotide table is retained only as `legacy_dinucleotide_*_screen`. The unsupported `stiffness = 2 - bendability` transformation is retired.

## Non-B structure

The module emits structure records for G4, i-motif, Z-DNA, triplex/H-DNA, cruciform, hairpin, slipped DNA, R-loop susceptibility, and A-tract curvature. Each record has strand, offsets, predictor/version, raw score, calibrated score, context status, experimental support, prediction agreement, and evidence level.

`non_B_DNA_aggregate_score` is `NA`. Unlike structures are not averaged into a mechanistic score.

## Primitive calibration

Canonical primitive scores are robustly standardized across the analysis cohort and combined using an observed-correlation covariance denominator. Each score carries its component list, weighting method, calibration status, robust-z status, and correlation audit. The old row-local means are available only as `*_legacy_screening_composite`.

Primitive construction emits no empirical p-value because it has no explicit null. Matched-control p-values are created later by the named null model.

## Unexplained anomaly

The old mean of primitive scores is retired. The new screen uses shrinkage robust Mahalanobis distance fitted and calibrated outside the held-out chromosome or non-overlapping genomic block. Without two groups, enough training rows, and two informative axes, the result is unavailable.

## Residual calibration

The canonical global metric is `classical_model_global_r2`, which may be negative under blocked cross-validation. Per-row residual outputs identify their z-score method and include robust unconditional diagnostics, cross-fitted conditional scale/z, held-out quantile residuals, conformal intervals, heteroscedasticity diagnostics, and block-bootstrap summary intervals where the data support them.

No residual p-value is synthesized from the observed score distribution.

## Negative evidence

Region cards and standalone negative-evidence outputs distinguish decisive failure, explicit negative evidence, missing controls, and unresolved evidence. Assembly artifacts, severe-null explanation, unstable window shifts, absent bridges, powered negative perturbations, failed rescue, and failed replication can reject or downgrade a candidate.
