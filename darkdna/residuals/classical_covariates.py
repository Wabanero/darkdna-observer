"""Classical covariate table assembly."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from darkdna.features.classical import CLASSICAL_COVARIATES, build_classical_covariates


def prepare_classical_covariates(windows: pd.DataFrame, features: pd.DataFrame | None = None) -> pd.DataFrame:
    return build_classical_covariates(windows, features)


def write_classical_covariates(covariates: pd.DataFrame, outdir: str | Path) -> Path:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "classical_covariates.parquet"
    covariates.to_parquet(path, index=False)
    return path
