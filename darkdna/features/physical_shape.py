"""Physical shape and curvature proxies."""

from __future__ import annotations

import math

import numpy as np

from darkdna.utils.optional_deps import optional_import
from .sequence import clean_sequence


DINUC_BENDABILITY = {
    "AA": 0.92,
    "AT": 0.82,
    "TA": 1.05,
    "CA": 1.10,
    "GT": 1.10,
    "CT": 1.02,
    "GA": 1.02,
    "CG": 0.78,
    "GC": 1.20,
    "GG": 0.95,
    "CC": 0.95,
    "AC": 0.96,
    "AG": 0.98,
    "TC": 0.98,
    "TG": 0.96,
    "TT": 0.92,
}

DINUC_STIFFNESS = {k: 2.0 - v for k, v in DINUC_BENDABILITY.items()}


def dinuc_values(seq: str, table: dict[str, float]) -> list[float]:
    seq = clean_sequence(seq)
    return [table[seq[i : i + 2]] for i in range(len(seq) - 1) if seq[i : i + 2] in table]


def compute_physical_shape_features(seq: str) -> dict[str, float | str]:
    bend = dinuc_values(seq, DINUC_BENDABILITY)
    stiff = dinuc_values(seq, DINUC_STIFFNESS)
    bendability = float(np.mean(bend)) if bend else math.nan
    curvature = float(np.std(bend)) if bend else math.nan
    stiffness = float(np.mean(stiff)) if stiff else math.nan

    rpy2, warning = optional_import("rpy2")
    shape_values = {
        "minor_groove_width_proxy": math.nan,
        "helix_twist_proxy": math.nan,
        "propeller_twist_proxy": math.nan,
        "roll_proxy": math.nan,
        "physical_shape_warning": warning or "",
    }
    if rpy2 is not None:  # pragma: no cover - R/DNAshapeR rarely available in CI.
        shape_values["physical_shape_warning"] = "DNAshapeR integration detected but not invoked by MVP fallback."
    return {
        "DNA_bendability_proxy": bendability,
        "DNA_curvature_proxy": curvature,
        "DNA_stiffness_proxy": stiffness,
        **shape_values,
    }
