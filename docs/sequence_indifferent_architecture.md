# Mode B — sequence-indifferent architecture

Mode B evaluates amount, length, copy number, presence/absence, spacing,
occupancy, heterochromatic mass, and replication-domain dimensions separately
from Mode A sequence anomaly.

The central counterfactual is whether a modelled effect remains stable after an
equal-length sequence replacement while changing under length or copy-number
titration. This is model-based perturbation evidence, not biological causality.

## Inputs

The only mandatory inputs are genomic intervals and their FASTA sequence.
Optional first-class inputs are repeat intervals, CNV/bedGraph copy number,
wide or long presence/absence matrices, syntenic interval tables, anchor sets,
occupancy, heterochromatin, replication domains, and phenotype tables.

When optional data are absent, the corresponding value is `NA`, the status is
`unavailable`, and a human-readable reason is emitted. Missing evidence is not
converted to zero.

Copy-number phenotype models require grouped validation. Random splits among
related samples, accessions, strains, haplotypes, overlapping loci, or family
copies are not used.

## Transformations

The control table includes native, reverse, reverse complement,
mononucleotide/dinucleotide/k-mer shuffles, GC- and repeat-matched equal-length
replacement, length titration at 0.25×–2×, copy-number titration at 0–8 copies,
and a length-by-copy factorial contrast.

Reported axes include sequence-identity, length, copy-number, composition,
repeat-fraction, and orientation sensitivity plus a sequence-by-quantity
interaction. Candidate labels are operational screening labels and include
sequence-specific, sequence-indifferent, mixed, length, copy-number, spacing,
heterochromatin, occupancy, artifact-compatible, and unresolved states.

Mode A and Mode B remain separate in `sequence_vs_quantity_scores.parquet`.
There is no universal DarkDNA score.

> A causal or quantity-dependent effect does not establish that the element
> originated or is maintained by selection for that effect.
