from pathlib import Path

import pytest

from tests.integration_helpers import assert_pipeline_outputs, run_cli_pipeline, write_temp_config


@pytest.mark.integration
def test_arabidopsis_reference_pipeline_skips_or_runs(tmp_path):
    root = Path("data/reference/arabidopsis_TAIR10_chr1_2Mb")
    required = [root / "genome.fa", root / "chrom.sizes", root / "genes.gff3"]
    if not all(path.exists() for path in required):
        pytest.skip("Optional arabidopsis_TAIR10_chr1_2Mb reference fixture is absent; run scripts/prepare_test_references.py to enable.")
    config = write_temp_config("configs/test_arabidopsis_TAIR10_chr1_2Mb.yaml", tmp_path, "arabidopsis_run")
    outdir = run_cli_pipeline(config)
    assert_pipeline_outputs(outdir)
