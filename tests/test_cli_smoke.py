from typer.testing import CliRunner

from darkdna.cli import app


def test_cli_smoke_init_and_toy_data(tmp_path):
    runner = CliRunner()
    result = runner.invoke(app, ["init", "--outdir", str(tmp_path / "init")])
    assert result.exit_code == 0
    result = runner.invoke(app, ["make-toy-data", "--outdir", str(tmp_path / "toy")])
    assert result.exit_code == 0
    assert (tmp_path / "toy" / "toy.fa").exists()
    assert (tmp_path / "toy" / "dark_windows.parquet").exists()
