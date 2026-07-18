import numpy as np
import pandas as pd

from darkdna.architecture.copy_number import compute_copy_number_features
from darkdna.cli import read_table


def test_copy_number_uses_grouped_cv_and_reports_sample_variation():
    intervals = pd.DataFrame({"region_id": ["r1"], "chrom": ["chr1"], "start": [0], "end": [100]})
    tracks = pd.DataFrame(
        {
            "chrom": ["chr1"] * 4,
            "start": [0] * 4,
            "end": [100] * 4,
            "sample_id": ["s1", "s2", "s3", "s4"],
            "copy_number": [1.0, 2.0, 3.0, 4.0],
        }
    )
    phenotype = pd.DataFrame(
        {
            "sample_id": ["s1", "s2", "s3", "s4"],
            "phenotype": [2.0, 4.0, 6.0, 8.0],
            "group_id": ["g1", "g1", "g2", "g2"],
        }
    )
    row = compute_copy_number_features(intervals, tracks, phenotype).iloc[0]
    assert row["copy_number_sample_size"] == 4
    assert row["copy_number_variance"] > 0
    assert row["copy_number_cv_status"] == "available_grouped_cv"
    assert np.isfinite(row["copy_number_effect_size"])


def test_structural_variant_vcf_cn_field_is_read_as_copy_number(tmp_path):
    path = tmp_path / "cnv.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        "chr1\t101\tcnv1\tN\t<CNV>\t.\tPASS\tEND=200;SVTYPE=CNV\tGT:CN\t0/1:1\t0/1:4\n",
        encoding="utf-8",
    )
    table = read_table(path)
    assert table["start"].tolist() == [100, 100]
    assert table["end"].tolist() == [200, 200]
    assert table.set_index("sample_id").loc["s2", "copy_number"] == 4
