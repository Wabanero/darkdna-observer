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


def _scale(value: float, vmin: float, vmax: float, low: float, high: float) -> float:
    if not np.isfinite(value) or vmax == vmin:
        return (low + high) / 2
    return low + ((value - vmin) / (vmax - vmin)) * (high - low)


def _primitive_label(value: object) -> str:
    text = str(value).replace("_candidate_score", "").replace("_candidate", "")
    return text.replace("_", " ")


def _write_multipanel_summary(
    residuals: pd.DataFrame,
    labels: pd.DataFrame | None,
    windows: pd.DataFrame | None,
    path: Path,
    title: str = "Candidate Summary Multipanel",
) -> None:
    width, height = 1120, 780
    margin = 32
    gap = 28
    panel_w = (width - (2 * margin) - gap) / 2
    panel_h = (height - 92 - gap) / 2
    panels = [
        (margin, 70, "Candidate counts"),
        (margin + panel_w + gap, 70, "Top candidate confidence"),
        (margin, 70 + panel_h + gap, "Residual z-score distribution"),
        (margin + panel_w + gap, 70 + panel_h + gap, "Genomic residual landscape"),
    ]
    parts = [
        f'<text x="{margin}" y="34" font-family="Arial" font-size="22" fill="#223c3b">{escape(title)}</text>',
        f'<text x="{margin}" y="56" font-family="Arial" font-size="12" fill="#555">Sequence-derived candidate overview for prioritization and follow-up validation.</text>',
    ]
    for x, y, panel_title in panels:
        parts.extend(
            [
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{panel_w:.1f}" height="{panel_h:.1f}" fill="white" stroke="#d7ded8" rx="4"/>',
                f'<text x="{x + 14:.1f}" y="{y + 24:.1f}" font-family="Arial" font-size="14" fill="#223c3b">{escape(panel_title)}</text>',
            ]
        )

    labels = labels if labels is not None else pd.DataFrame()
    windows = windows if windows is not None else pd.DataFrame()

    # Panel 1: primitive label counts.
    x0, y0, _ = panels[0]
    if labels.empty or "primitive_class" not in labels.columns:
        parts.append(f'<text x="{x0 + 18:.1f}" y="{y0 + 58:.1f}" font-family="Arial" font-size="12" fill="#555">No candidate labels available.</text>')
    else:
        counts = labels["primitive_class"].astype(str).value_counts().head(8)
        max_count = max(1, int(counts.max()))
        bar_left = x0 + 190
        bar_top = y0 + 46
        row_h = max(22, (panel_h - 72) / max(1, len(counts)))
        max_bar_w = panel_w - 230
        for idx, (primitive, count) in enumerate(counts.items()):
            y = bar_top + idx * row_h
            label = escape(_primitive_label(primitive)[:26])
            bar_w = (float(count) / max_count) * max_bar_w
            parts.append(f'<text x="{x0 + 18:.1f}" y="{y + 14:.1f}" font-family="Arial" font-size="11" fill="#333">{label}</text>')
            parts.append(f'<rect x="{bar_left:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="14" fill="#426b69"/>')
            parts.append(f'<text x="{bar_left + bar_w + 6:.1f}" y="{y + 12:.1f}" font-family="Arial" font-size="11" fill="#333">{int(count)}</text>')

    # Panel 2: top confidence candidates.
    x0, y0, _ = panels[1]
    needed = {"region_id", "primitive_class", "primitive_confidence"}
    if labels.empty or not needed.issubset(labels.columns):
        parts.append(f'<text x="{x0 + 18:.1f}" y="{y0 + 58:.1f}" font-family="Arial" font-size="12" fill="#555">No confidence values available.</text>')
    else:
        top = labels.copy()
        top["primitive_confidence"] = pd.to_numeric(top["primitive_confidence"], errors="coerce")
        top = top.dropna(subset=["primitive_confidence"]).sort_values("primitive_confidence", ascending=False).head(8)
        max_conf = max(1e-9, float(top["primitive_confidence"].max())) if not top.empty else 1.0
        bar_left = x0 + 200
        bar_top = y0 + 46
        row_h = max(22, (panel_h - 72) / max(1, len(top)))
        max_bar_w = panel_w - 245
        for idx, row in enumerate(top.itertuples(index=False)):
            y = bar_top + idx * row_h
            primitive = escape(_primitive_label(getattr(row, "primitive_class"))[:24])
            confidence = float(getattr(row, "primitive_confidence"))
            bar_w = (confidence / max_conf) * max_bar_w
            parts.append(f'<text x="{x0 + 18:.1f}" y="{y + 14:.1f}" font-family="Arial" font-size="11" fill="#333">{primitive}</text>')
            parts.append(f'<rect x="{bar_left:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="14" fill="#7b4b73"/>')
            parts.append(f'<text x="{bar_left + bar_w + 6:.1f}" y="{y + 12:.1f}" font-family="Arial" font-size="11" fill="#333">{confidence:.2f}</text>')

    # Panel 3: residual histogram.
    x0, y0, _ = panels[2]
    values = _finite(residuals.get("residual_zscore", pd.Series(dtype=float)))
    plot_left = x0 + 52
    plot_top = y0 + 48
    plot_w = panel_w - 82
    plot_h = panel_h - 88
    if values.size == 0:
        parts.append(f'<text x="{x0 + 18:.1f}" y="{y0 + 58:.1f}" font-family="Arial" font-size="12" fill="#555">No finite residuals available.</text>')
    else:
        counts, edges = np.histogram(values, bins=min(24, max(1, values.size)))
        max_count = max(1, int(counts.max()))
        bin_w = plot_w / max(1, len(counts))
        parts.append(f'<line x1="{plot_left:.1f}" y1="{plot_top + plot_h:.1f}" x2="{plot_left + plot_w:.1f}" y2="{plot_top + plot_h:.1f}" stroke="#65736f"/>')
        parts.append(f'<line x1="{plot_left:.1f}" y1="{plot_top:.1f}" x2="{plot_left:.1f}" y2="{plot_top + plot_h:.1f}" stroke="#65736f"/>')
        for idx, count in enumerate(counts):
            bar_h = (float(count) / max_count) * plot_h
            x = plot_left + idx * bin_w + 1
            y = plot_top + plot_h - bar_h
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(1.0, bin_w - 2):.1f}" height="{bar_h:.1f}" fill="#426b69" opacity="0.85"/>')
        parts.append(f'<text x="{plot_left:.1f}" y="{y0 + panel_h - 14:.1f}" font-family="Arial" font-size="11" fill="#555">{edges[0]:.2f}</text>')
        parts.append(f'<text x="{plot_left + plot_w - 30:.1f}" y="{y0 + panel_h - 14:.1f}" font-family="Arial" font-size="11" fill="#555">{edges[-1]:.2f}</text>')

    # Panel 4: strongest residual by region along genomic position.
    x0, y0, _ = panels[3]
    if residuals.empty or windows.empty or not {"region_id", "residual_zscore"}.issubset(residuals.columns) or not {"region_id", "start", "end"}.issubset(windows.columns):
        parts.append(f'<text x="{x0 + 18:.1f}" y="{y0 + 58:.1f}" font-family="Arial" font-size="12" fill="#555">No genomic residual positions available.</text>')
    else:
        strongest = residuals.copy()
        strongest["residual_zscore"] = pd.to_numeric(strongest["residual_zscore"], errors="coerce")
        strongest["abs_z"] = strongest["residual_zscore"].abs()
        strongest = strongest.dropna(subset=["abs_z"])
        if strongest.empty:
            parts.append(f'<text x="{x0 + 18:.1f}" y="{y0 + 58:.1f}" font-family="Arial" font-size="12" fill="#555">No finite genomic residuals available.</text>')
        else:
            strongest = strongest.loc[strongest.groupby("region_id")["abs_z"].idxmax()]
            coords = windows[["region_id", "start", "end"]].drop_duplicates("region_id")
            strongest = strongest.merge(coords, on="region_id", how="left").dropna(subset=["start", "end"])
            strongest["midpoint"] = (pd.to_numeric(strongest["start"], errors="coerce") + pd.to_numeric(strongest["end"], errors="coerce")) / 2
            strongest = strongest.dropna(subset=["midpoint"])
            plot_left = x0 + 54
            plot_top = y0 + 48
            plot_w = panel_w - 86
            plot_h = panel_h - 92
            if strongest.empty:
                parts.append(f'<text x="{x0 + 18:.1f}" y="{y0 + 58:.1f}" font-family="Arial" font-size="12" fill="#555">No finite genomic residuals available.</text>')
            else:
                xmin, xmax = float(strongest["midpoint"].min()), float(strongest["midpoint"].max())
                zmax = max(1.0, float(strongest["residual_zscore"].abs().max()))
                parts.append(f'<line x1="{plot_left:.1f}" y1="{plot_top + plot_h / 2:.1f}" x2="{plot_left + plot_w:.1f}" y2="{plot_top + plot_h / 2:.1f}" stroke="#c8d1cc"/>')
                parts.append(f'<line x1="{plot_left:.1f}" y1="{plot_top:.1f}" x2="{plot_left:.1f}" y2="{plot_top + plot_h:.1f}" stroke="#65736f"/>')
                parts.append(f'<line x1="{plot_left:.1f}" y1="{plot_top + plot_h:.1f}" x2="{plot_left + plot_w:.1f}" y2="{plot_top + plot_h:.1f}" stroke="#65736f"/>')
                for row in strongest.itertuples(index=False):
                    x = _scale(float(getattr(row, "midpoint")), xmin, xmax, plot_left, plot_left + plot_w)
                    z = float(getattr(row, "residual_zscore"))
                    y = _scale(z, -zmax, zmax, plot_top + plot_h, plot_top)
                    color = "#7b4b73" if z >= 0 else "#426b69"
                    radius = 2.5 + min(4.0, abs(z) * 0.55)
                    parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius:.1f}" fill="{color}" opacity="0.72"/>')
                parts.append(f'<text x="{plot_left:.1f}" y="{y0 + panel_h - 14:.1f}" font-family="Arial" font-size="11" fill="#555">{xmin / 1000:.1f} kb</text>')
                parts.append(f'<text x="{plot_left + plot_w - 48:.1f}" y="{y0 + panel_h - 14:.1f}" font-family="Arial" font-size="11" fill="#555">{xmax / 1000:.1f} kb</text>')
                parts.append(f'<text x="{x0 + 12:.1f}" y="{plot_top + 10:.1f}" font-family="Arial" font-size="11" fill="#555">+z</text>')
                parts.append(f'<text x="{x0 + 12:.1f}" y="{plot_top + plot_h:.1f}" font-family="Arial" font-size="11" fill="#555">-z</text>')

    path.write_text(_svg_frame(width, height, title, "\n".join(parts)), encoding="utf-8")


def write_basic_plots(
    residuals: pd.DataFrame,
    outdir: str | Path,
    labels: pd.DataFrame | None = None,
    windows: pd.DataFrame | None = None,
) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    multipanel_path = out / "multipanel_summary.svg"
    _write_multipanel_summary(residuals, labels, windows, multipanel_path)
    paths["multipanel_summary"] = multipanel_path

    hist_path = out / "residual_score_histogram.svg"
    _write_histogram(_finite(residuals.get("residual_zscore", pd.Series(dtype=float))), hist_path, "Residual Score Histogram")
    paths["residual_score_histogram"] = hist_path

    scatter_path = out / "observed_vs_predicted_classical_score.svg"
    predicted = pd.to_numeric(residuals.get("predicted_classical_score", pd.Series(dtype=float)), errors="coerce").to_numpy(dtype=float)
    observed = pd.to_numeric(residuals.get("observed_score", pd.Series(dtype=float)), errors="coerce").to_numpy(dtype=float)
    _write_scatter(predicted, observed, scatter_path, "Observed vs Predicted Classical Score")
    paths["observed_vs_predicted_classical_score"] = scatter_path
    return paths
