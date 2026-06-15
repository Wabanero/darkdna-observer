"""Diagnostic plots."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def write_basic_plots(residuals: pd.DataFrame, outdir: str | Path) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 4))
        ax.hist(residuals["residual_zscore"].dropna(), bins=30, color="#426b69", edgecolor="white")
        ax.set_xlabel("Residual z-score")
        ax.set_ylabel("Count")
        ax.set_title("Residual Score Histogram")
        path = out / "residual_score_histogram.png"
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
        paths["residual_score_histogram"] = path

        fig, ax = plt.subplots(figsize=(5, 4))
        ax.scatter(residuals["predicted_classical_score"], residuals["observed_score"], s=12, alpha=0.6, color="#7b4b73")
        ax.set_xlabel("Predicted classical score")
        ax.set_ylabel("Observed score")
        ax.set_title("Observed vs Predicted Classical Score")
        path = out / "observed_vs_predicted_classical_score.png"
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
        paths["observed_vs_predicted_classical_score"] = path
    except Exception:
        pass
    return paths
