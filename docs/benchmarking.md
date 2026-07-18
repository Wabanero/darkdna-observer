# Native-versus-randomized default-state benchmark

Run:

```bash
darkdna benchmark-default-state --config config.yaml
```

The benchmark compares native sequence with whole-genome reversal, reverse
complement, global and local mononucleotide and dinucleotide shuffles,
k-mer-preserving shuffles, and a reference-conditioned evolutionary-process
genome. Repeat-only, non-repeat-only, and TE-orientation controls are marked
unavailable unless their required annotations are supplied.

The key output is not a binary decision. It is the shift under each null and the
sensitivity of the conclusion to null definition.

Outputs:

- `native_vs_null_feature_shift.parquet`;
- `native_vs_null_primitive_shift.parquet`;
- `transcription_initiation_background.parquet`;
- `null_method_sensitivity.parquet`;
- `default_state_benchmark.html`.

Puffin/Puffin-D, Enformer, Borzoi, AlphaGenome, or future models require a local
callable adapter. A configured path without an adapter is reported as
unavailable. Analysis never downloads a model or makes a network call.

Native excess is statistical structure and may reflect selection, mutation
bias, repeat history, or other processes. It is not evidence that every active
or anomalous interval is functional.
