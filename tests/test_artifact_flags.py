from darkdna.toy_data import make_toy_data
from darkdna.windows.make_windows import make_dark_windows


def test_artifact_flags_for_unplaced_low_mappability_contig(tmp_path):
    paths = make_toy_data(tmp_path)
    windows = make_dark_windows(
        fasta=paths["fasta"],
        blacklist_path=paths["blacklist"],
        mappability_path=paths["mappability"],
        window_sizes=[200],
        exclude_coding_exons=False,
    )
    flagged = windows[windows["chrom"] == "contig_unplaced_01"]["artifact_risk_flags"].str.cat(sep=";")
    assert "low_mappability" in flagged
    assert "unplaced_or_unlocalized_contig" in flagged
