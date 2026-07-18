import pandas as pd

from darkdna.benchmarks.default_state import benchmark_default_state, write_default_state_benchmark


def test_default_state_benchmark_reports_null_method_dependence_and_skips_missing_masks(tmp_path):
    genome = {"chr1": ("ACGT" * 80) + ("AAAAACCCCC" * 20)}
    windows = pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "chrom": ["chr1", "chr1"],
            "start": [0, 160],
            "end": [120, 280],
        }
    )
    results = benchmark_default_state(
        genome,
        windows,
        seed=7,
        max_windows=2,
        methods=(
            "native",
            "whole_genome_reverse",
            "global_mononucleotide_shuffle",
            "evolutionary_process_generated",
            "repeat_only_shuffle",
        ),
    )
    unavailable = results["feature_rows"].query("null_method == 'repeat_only_shuffle'")
    assert unavailable.iloc[0]["status"] == "unavailable"
    assert set(results["feature_shift"]["null_method"]) >= {
        "native",
        "whole_genome_reverse",
        "global_mononucleotide_shuffle",
        "evolutionary_process_generated",
    }
    assert (results["feature_shift"]["interpretation"] == "statistical_structure_shift_not_selected_function").all()
    paths = write_default_state_benchmark(results, tmp_path)
    assert all(path.exists() for path in paths.values())


def test_default_state_benchmark_uses_supplied_repeat_and_te_intervals():
    genome = {"chr1": "ACGT" * 100}
    windows = pd.DataFrame({"region_id": ["r1"], "chrom": ["chr1"], "start": [0], "end": [120]})
    repeats = pd.DataFrame({"chrom": ["chr1"], "start": [20], "end": [80], "strand": ["+"]})
    results = benchmark_default_state(
        genome,
        windows,
        methods=("native", "repeat_only_shuffle", "non_repeat_only_shuffle", "TE_orientation_reversal"),
        repeat_intervals=repeats,
        te_annotations=repeats,
    )
    rows = results["feature_rows"]
    assert set(rows["null_method"]) == {
        "native",
        "repeat_only_shuffle",
        "non_repeat_only_shuffle",
        "TE_orientation_reversal",
    }
    assert rows["status"].eq("available").all()
