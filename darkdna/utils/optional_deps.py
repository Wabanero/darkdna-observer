"""Optional dependency helpers.

The MVP treats physical-shape, fractal, interval, browser-track, and R-based
integrations as feature coverage boosters. Missing packages should produce
warnings and NaN/fallback values, never crash the main sequence-first path.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from typing import Any


OPTIONAL_DEPENDENCIES = {
    "Bio": "FASTA parsing fallback is used when Biopython is unavailable.",
    "pyfaidx": "Pure-Python FASTA reader is used when pyfaidx is unavailable.",
    "pyfastx": "Pure-Python FASTA reader is used when pyfastx is unavailable.",
    "skbio": "Native lightweight sequence summaries are used when scikit-bio is unavailable.",
    "pyranges": "Pure-Python interval overlap is used when pyranges is unavailable.",
    "pybedtools": "Pure-Python interval overlap is used when pybedtools is unavailable.",
    "pyBigWig": "BED/bedGraph fallbacks are used when pyBigWig is unavailable.",
    "pywt": "Wavelet summaries are skipped when PyWavelets is unavailable.",
    "nolds": "DFA/Hurst fallback estimators are used when nolds is unavailable.",
    "antropy": "Native entropy estimates are used when antropy is unavailable.",
    "dnacurve": "Simple dinucleotide physical-shape proxies are used when dnacurve is unavailable.",
    "rpy2": "DNAshapeR-backed shape proxies are skipped when rpy2 is unavailable.",
    "xgboost": "sklearn models are used when xgboost is unavailable.",
    "lightgbm": "sklearn models are used when lightgbm is unavailable.",
    "statsmodels": "sklearn/NumPy regressions are used when statsmodels is unavailable.",
}


def optional_import(module_name: str) -> tuple[Any | None, str | None]:
    """Import an optional module and return ``(module, warning)``."""

    try:
        return importlib.import_module(module_name), None
    except Exception as exc:  # pragma: no cover - exact optional failures vary.
        reason = OPTIONAL_DEPENDENCIES.get(module_name, "Optional feature skipped.")
        return None, f"{module_name} unavailable ({exc}); {reason}"


def is_available(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


def optional_dependency_report() -> dict[str, dict[str, Any]]:
    report: dict[str, dict[str, Any]] = {}
    for name, fallback in OPTIONAL_DEPENDENCIES.items():
        available = is_available(name)
        version = None
        if available:
            mod, _ = optional_import(name)
            version = getattr(mod, "__version__", None) if mod is not None else None
        report[name] = {"available": available, "version": version, "fallback": fallback}
    return report


def write_optional_dependency_report(path: str | Path) -> None:
    Path(path).write_text(json.dumps(optional_dependency_report(), indent=2), encoding="utf-8")
