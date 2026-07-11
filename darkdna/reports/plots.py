"""Dependency-light diagnostic plots."""

from __future__ import annotations

from html import escape
from pathlib import Path

import numpy as np
import pandas as pd


def _finite(values: pd.Series) -> np.ndarray:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    return arr[np.isfinite(arr)]


def _svg_frame(width: int, height: int, title: str, body: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n'
        '<rect width="100%" height="100%" fill="#fbfbf8"/>\n'
        f'<text x="24" y="28" font-family="Arial" font-size="16" fill="#223c3b">{escape(title)}</text>\n'
        f"{body}\n"
        "</svg>\n"
    )


def _write_histogram(values: np.ndarray, path: Path, title: str) -> None:
    width, height = 620, 420
    left, right, top, bottom = 54, 24, 52, 44
    plot_w = width - left - right
    plot_h = height - top - bottom
    if values.size == 0:
        body = '<text x="54" y="96" font-family="Arial" font-size="13" fill="#555">No finite values.</text>'
        path.write_text(_svg_frame(width, height, title, body), encoding="utf-8")
        return
    counts, edges = np.histogram(values, bins=min(30, max(1, values.size)))
    max_count = max(1, int(counts.max()))
    parts = [
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#65736f"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#65736f"/>',
    ]
    bin_w = plot_w / max(1, len(counts))
    for idx, count in enumerate(counts):
        bar_h = (count / max_count) * plot_h
        x = left + idx * bin_w + 1
        y = top + plot_h - bar_h
        parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(1.0, bin_w - 2):.2f}" height="{bar_h:.2f}" fill="#426b69"/>')
    parts.extend(
        [
            f'<text x="{left}" y="{height - 14}" font-family="Arial" font-size="12" fill="#333">Residual z-score</text>',
            f'<text x="8" y="{top + 14}" font-family="Arial" font-size="12" fill="#333">Count</text>',
            f'<text x="{left}" y="{top + plot_h + 18}" font-family="Arial" font-size="11" fill="#555">{edges[0]:.2g}</text>',
            f'<text x="{left + plot_w - 36}" y="{top + plot_h + 18}" font-family="Arial" font-size="11" fill="#555">{edges[-1]:.2g}</text>',
        ]
    )
    path.write_text(_svg_frame(width, height, title, "\n".join(parts)), encoding="utf-8")


def _write_scatter(x_values: np.ndarray, y_values: np.ndarray, path: Path, title: str) -> None:
    width, height = 620, 420
    left, right, top, bottom = 64, 26, 52, 48
    plot_w = width - left - right
    plot_h = height - top - bottom
    mask = np.isfinite(x_values) & np.isfinite(y_values)
    x_values = x_values[mask]
    y_values = y_values[mask]
    if x_values.size == 0:
        body = '<text x="64" y="96" font-family="Arial" font-size="13" fill="#555">No finite paired values.</text>'
        path.write_text(_svg_frame(width, height, title, body), encoding="utf-8")
        return
    xmin, xmax = float(x_values.min()), float(x_values.max())
    ymin, ymax = float(y_values.min()), float(y_values.max())
    if xmin == xmax:
        xmin -= 1.0
        xmax += 1.0
    if ymin == ymax:
        ymin -= 1.0
        ymax += 1.0
    parts = [
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#65736f"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#65736f"/>',
    ]
    for x_val, y_val in zip(x_values, y_values):
        x = left + ((float(x_val) - xmin) / (xmax - xmin)) * plot_w
        y = top + plot_h - ((float(y_val) - ymin) / (ymax - ymin)) * plot_h
        parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3" fill="#7b4b73" opacity="0.65"/>')
    parts.extend(
        [
            f'<text x="{left}" y="{height - 14}" font-family="Arial" font-size="12" fill="#333">Predicted classical score</text>',
            f'<text x="8" y="{top + 14}" font-family="Arial" font-size="12" fill="#333">Observed</text>',
        ]
    )
    path.write_text(_svg_frame(width, height, title, "\n".join(parts)), encoding="utf-8")


def write_basic_plots(residuals: pd.DataFrame, outdir: str | Path) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    hist_path = out / "residual_score_histogram.svg"
    _write_histogram(_finite(residuals.get("residual_zscore", pd.Series(dtype=float))), hist_path, "Residual Score Histogram")
    paths["residual_score_histogram"] = hist_path

    scatter_path = out / "observed_vs_predicted_classical_score.svg"
    predicted = pd.to_numeric(residuals.get("predicted_classical_score", pd.Series(dtype=float)), errors="coerce").to_numpy(dtype=float)
    observed = pd.to_numeric(residuals.get("observed_score", pd.Series(dtype=float)), errors="coerce").to_numpy(dtype=float)
    _write_scatter(predicted, observed, scatter_path, "Observed vs Predicted Classical Score")
    paths["observed_vs_predicted_classical_score"] = scatter_path
    return paths
