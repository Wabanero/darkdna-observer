from pathlib import Path

import pytest

from tests.integration_helpers import assert_pipeline_outputs, run_cli_pipeline, write_temp_config


@pytest.mark.integration
def test_yeast_reference_pipeline_skips_or_runs(tmp_path):
    root = Path("data/reference/yeast_R64_chrI")
    required = [root / "genome.fa", root / "chrom.sizes", root / "genes.gff3"]
    if not all(path.exists() for path in required):
        pytest.skip("Optional yeast_R64_chrI reference fixture is absent; run scripts/prepare_test_references.py to enable.")
    config = write_temp_config("configs/test_yeast_R64_chrI.yaml", tmp_path, "yeast_run")
    outdir = run_cli_pipeline(config)
    assert_pipeline_outputs(outdir)
