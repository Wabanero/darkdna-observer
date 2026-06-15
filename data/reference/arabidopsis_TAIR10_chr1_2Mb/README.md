# arabidopsis_TAIR10_chr1_2Mb

Optional Level C plant/non-model real-genome integration fixture.

This directory intentionally does not include full Arabidopsis reference downloads by default. Generate the small Chr1 2 Mb subset with:

```bash
python scripts/prepare_test_references.py --dataset arabidopsis_TAIR10_chr1_2Mb --no-download --source-dir /path/to/local/source --out data/reference/arabidopsis_TAIR10_chr1_2Mb
```

Normal pytest skips the Arabidopsis integration test unless `genome.fa`, `chrom.sizes`, and `genes.gff3` exist.
