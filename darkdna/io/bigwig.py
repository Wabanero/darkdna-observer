"""bigWig/bedGraph access with graceful fallback."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .bed import read_bedgraph, weighted_interval_mean
from darkdna.utils.optional_deps import optional_import


def mean_signal(path: str | Path | None, chrom: str, start: int, end: int) -> float:
    if not path:
        return np.nan
    p = Path(path)
    if p.suffix.lower() in {".bedgraph", ".bed", ".bg"}:
        return weighted_interval_mean(chrom, start, end, read_bedgraph(p))
    pybigwig, warning = optional_import("pyBigWig")
    if pybigwig is None:
        return np.nan
    try:  # pragma: no cover - pyBigWig availability varies.
        bw = pybigwig.open(str(p))
        value = bw.stats(chrom, int(start), int(end), type="mean")[0]
        bw.close()
        return float(value) if value is not None else np.nan
    except Exception:
        return np.nan


def annotate_mean_signal(windows: pd.DataFrame, path: str | Path | None, column: str = "mappability") -> pd.DataFrame:
    out = windows.copy()
    if not path:
        out[column] = np.nan
        return out
    p = Path(path)
    if p.suffix.lower() in {".bedgraph", ".bed", ".bg"}:
        intervals = read_bedgraph(p)
        out[column] = [weighted_interval_mean(row.chrom, row.start, row.end, intervals) for row in out.itertuples()]
        return out
    out[column] = [mean_signal(p, row.chrom, row.start, row.end) for row in out.itertuples()]
    return out
