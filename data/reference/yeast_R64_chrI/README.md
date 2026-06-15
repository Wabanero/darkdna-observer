# yeast_R64_chrI

Optional Level B real-genome technical integration fixture.

This directory intentionally does not include full yeast reference downloads by default. Generate the small chromosome-I subset with:

```bash
python scripts/prepare_test_references.py --dataset yeast_R64_chrI --no-download --source-dir /path/to/local/source --out data/reference/yeast_R64_chrI
```

Normal pytest skips the yeast integration test unless `genome.fa`, `chrom.sizes`, and `genes.gff3` exist.
