"""Table-level multiscale profile helpers."""

from __future__ import annotations

import pandas as pd

from .scale_fractal import add_multiscale_profiles


def compute_multiscale_profiles(features: pd.DataFrame, windows: pd.DataFrame | None = None) -> pd.DataFrame:
    return add_multiscale_profiles(features, windows=windows)
