"""Classical covariate table assembly."""

from __future__ import annotations

import pandas as pd

from darkdna.features.classical import CLASSICAL_COVARIATES, build_classical_covariates


def prepare_classical_covariates(windows: pd.DataFrame, features: pd.DataFrame | None = None) -> pd.DataFrame:
    return build_classical_covariates(windows, features)
