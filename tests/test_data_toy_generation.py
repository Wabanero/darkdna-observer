import pandas as pd

from darkdna.io.fasta import read_fasta
from darkdna.toy_data import make_toy_data


def test_level_a_toy_generation_schema(tmp_path):
    paths = make_toy_data(tmp_path / "toy", seed=42)
    expected_files = [
        "toy_genome.fa",
        "toy_genome.fa.fai",
        "toy.chrom.sizes",
        "toy_annotation.gff3",
        "toy_te_annotation.gff3",
        "toy_blacklist.bed",
        "toy_mappability.bedGraph",
        "expected_candidates.tsv",
        "README.md",
    ]
    for name in expected_files:
        assert (tmp_path / "toy" / name).exists(), name
    genome = read_fasta(paths["fasta"])
    assert {name: len(seq) for name, seq in genome.items()} == {
        "toy_chr1": 100_000,
        "toy_chr2": 80_000,
        "toy_scaffoldA": 50_000,
    }
    expected = pd.read_csv(tmp_path / "toy" / "expected_candidates.tsv", sep="\t")
    assert set(expected.columns) == {
        "region_label",
        "chrom",
        "start",
        "end",
        "expected_candidate_type",
        "expected_high_score_family",
    }
    assert len(expected) >= 9
