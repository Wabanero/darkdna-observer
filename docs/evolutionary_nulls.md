# Evolutionary and severe nulls

## Interpretation boundary

DarkDNA-Observer tests observed structure against named null families. A signal
that exceeds a null is an anomaly relative to that null; it is not selected
function, adaptive origin, or proof that the null is the true evolutionary
history.

## Registry and availability

`null_model_registry.json` contains composition, local-genomic, gene-distance,
TE family/subfamily/age, repeat-density, chromatin, replication, recombination,
mutation, damage, mappability, assembly, synteny, population, copy-number,
presence/absence, orientation, synthetic-composition, and evolutionary-process
families. Each family declares its execution mode and required inputs.

For every candidate and primitive, the summary reports:

- `available_null_models` and `missing_null_models`;
- `null_model_count`;
- `null_model_agreement` and `null_model_conflict`;
- `null_panel_status`;
- the conservative minimum survival z-score and maximum empirical p-value among
  calibrated named nulls.

Controls from the same genomic block are excluded. A panel with too few
independent blocks remains `partial_null_panel_not_for_promotion`.

## Evolutionary-process surrogate

`darkdna build-evolutionary-nulls` fits base frequencies and first-order
transitions from the supplied reference and supports a user mutation-spectrum
table. The simulator includes explicit CpG, GC-bias, indel-turnover,
microsatellite, and repeat-turnover components while preserving benchmark
interval length.

Without organism-specific mutation data the status is
`reference_conditioned_generic_processes`. The model is not an ancestral
reconstruction. TE insertion remains unavailable unless a local TE library and
insertion-preference model are supplied; no model is downloaded.

Outputs are `evolutionary_null_scores.parquet` and
`evolutionary_model_manifest.json`.
