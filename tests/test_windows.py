from darkdna.io.fasta import read_fasta
from darkdna.toy_data import make_toy_data
from darkdna.windows.make_windows import make_dark_windows


def test_window_generation_excludes_exons_and_keeps_covariates(tmp_path):
    paths = make_toy_data(tmp_path)
    windows = make_dark_windows(
        fasta=paths["fasta"],
        annotation_path=paths["annotation_gff3"],
        te_annotation_path=paths["te_gff3"],
        blacklist_path=paths["blacklist"],
        mappability_path=paths["mappability"],
        window_sizes=[200],
        exclude_coding_exons=True,
    )
    assert not windows.empty
    assert "region_id" in windows.columns
    assert not windows["overlaps_exon"].any()
    assert windows["overlaps_TE"].any()
    assert "artifact_risk_flags" in windows.columns


def test_fasta_extraction_and_scaffold_support(tmp_path):
    paths = make_toy_data(tmp_path)
    genome = read_fasta(paths["fasta"])
    assert "toy_chr1" in genome
    assert "toy_chr2" in genome
    assert "toy_scaffoldA" in genome
    assert len(genome["toy_chr1"]) == 100_000
