from typer.testing import CliRunner

from darkdna.cli import app


def test_existing_sequence_first_commands_remain_public():
    result = CliRunner().invoke(app, ["--help"])

    assert result.exit_code == 0
    for command in (
        "make-windows",
        "extract-sequence-features",
        "score-primitives",
        "build-null-models",
        "residualize",
        "infer-primitives",
        "make-region-cards",
        "report",
        "make-tracks",
        "run",
    ):
        assert command in result.output


def test_phase2_and_phase3_commands_are_public_without_removing_mode_a():
    result = CliRunner().invoke(app, ["--help"])
    assert result.exit_code == 0
    for command in (
        "build-evolutionary-nulls",
        "benchmark-default-state",
        "extract-architecture-features",
        "score-sequence-indifferent-architecture",
        "compare-sequence-vs-quantity",
        "infer-architecture-candidates",
    ):
        assert command in result.output
