"""Dependency-light diagnostic plots."""

from __future__ import annotations

from html import escape
from pathlib import Path

import numpy as np
import pandas as pd


def _finite(values: pd.Series) -> np.ndarray:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    return arr[np.isfinite(arr)]


def _svg_frame(width: int, height: int, title: str, body: str, *, show_title: bool = True) -> str:
    title_text = f'<text x="24" y="28" font-family="Arial" font-size="16" fill="#223c3b">{escape(title)}</text>\n' if show_title else ""
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n'
        '<rect width="100%" height="100%" fill="#fbfbf8"/>\n'
        f"{title_text}"
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


def _add_x_ticks(
    parts: list[str],
    ticks: list[float],
    vmin: float,
    vmax: float,
    left: float,
    baseline: float,
    width: float,
    *,
    label_y: float | None = None,
    suffix: str = "",
) -> None:
    label_y = baseline + 16 if label_y is None else label_y
    for tick in ticks:
        x = _scale(tick, vmin, vmax, left, left + width)
        parts.append(f'<line x1="{x:.1f}" y1="{baseline:.1f}" x2="{x:.1f}" y2="{baseline + 4:.1f}" stroke="#65736f"/>')
        parts.append(f'<text x="{x - 12:.1f}" y="{label_y:.1f}" font-family="Arial" font-size="10" fill="#555">{tick:g}{suffix}</text>')


def _add_y_ticks(
    parts: list[str],
    ticks: list[float],
    vmin: float,
    vmax: float,
    left: float,
    top: float,
    height: float,
    *,
    label_x: float | None = None,
    suffix: str = "",
) -> None:
    label_x = left - 34 if label_x is None else label_x
    for tick in ticks:
        y = _scale(tick, vmin, vmax, top + height, top)
        parts.append(f'<line x1="{left - 4:.1f}" y1="{y:.1f}" x2="{left:.1f}" y2="{y:.1f}" stroke="#65736f"/>')
        parts.append(f'<text x="{label_x:.1f}" y="{y + 3:.1f}" font-family="Arial" font-size="10" fill="#555">{tick:g}{suffix}</text>')


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

    path.write_text(_svg_frame(width, height, title, "\n".join(parts), show_title=False), encoding="utf-8")


def _write_classical_control_multipanel(
    residuals: pd.DataFrame,
    labels: pd.DataFrame | None,
    path: Path,
    title: str = "Classical Explanation Removal Multipanel",
) -> None:
    width, height = 1120, 780
    margin = 32
    gap = 28
    panel_w = (width - (2 * margin) - gap) / 2
    panel_h = (height - 92 - gap) / 2
    panels = [
        (margin, 70, "Classical explanation fraction"),
        (margin + panel_w + gap, 70, "Observed vs predicted classical score"),
        (margin, 70 + panel_h + gap, "Residual vs matched-null evidence"),
        (margin + panel_w + gap, 70 + panel_h + gap, "Top post-control candidates"),
    ]
    parts = [
        f'<text x="{margin}" y="34" font-family="Arial" font-size="22" fill="#223c3b">{escape(title)}</text>',
        f'<text x="{margin}" y="56" font-family="Arial" font-size="12" fill="#555">How much primitive signal is predicted by classical covariates, and what remains after subtraction.</text>',
        '<rect x="690" y="21" width="10" height="10" fill="#9b5f3f"/>',
        '<text x="705" y="30" font-family="Arial" font-size="10" fill="#555">mostly classical</text>',
        '<circle cx="812" cy="26" r="4" fill="#7b4b73" opacity="0.72"/>',
        '<text x="823" y="30" font-family="Arial" font-size="10" fill="#555">strong post-control</text>',
        '<line x1="948" y1="26" x2="964" y2="26" stroke="#c6912f" stroke-dasharray="3 4"/>',
        '<text x="970" y="30" font-family="Arial" font-size="10" fill="#555">review threshold</text>',
    ]
    for x, y, panel_title in panels:
        parts.extend(
            [
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{panel_w:.1f}" height="{panel_h:.1f}" fill="white" stroke="#d7ded8" rx="4"/>',
                f'<text x="{x + 14:.1f}" y="{y + 24:.1f}" font-family="Arial" font-size="14" fill="#223c3b">{escape(panel_title)}</text>',
            ]
        )

    labels = labels if labels is not None else pd.DataFrame()

    # Panel 1: R2-like classical explanation fraction by primitive.
    x0, y0, _ = panels[0]
    r2_column = "classical_model_global_r2" if "classical_model_global_r2" in residuals.columns else "classical_explanation_fraction"
    if residuals.empty or "primitive" not in residuals.columns or r2_column not in residuals.columns:
        parts.append(f'<text x="{x0 + 18:.1f}" y="{y0 + 58:.1f}" font-family="Arial" font-size="12" fill="#555">No classical explanation fractions available.</text>')
    else:
        explained = residuals[["primitive", r2_column]].drop_duplicates("primitive").copy()
        explained[r2_column] = pd.to_numeric(explained[r2_column], errors="coerce").fillna(0.0)
        explained = explained.sort_values(r2_column, ascending=False).head(8)
        bar_left = x0 + 210
        bar_top = y0 + 46
        row_h = max(22, (panel_h - 74) / max(1, len(explained)))
        max_bar_w = panel_w - 255
        for idx, row in enumerate(explained.itertuples(index=False)):
            y = bar_top + idx * row_h
            value = float(getattr(row, r2_column))
            label = escape(_primitive_label(getattr(row, "primitive"))[:28])
            bar_w = max(0.0, min(1.0, value)) * max_bar_w
            color = "#9b5f3f" if value >= 0.7 else "#426b69"
            parts.append(f'<text x="{x0 + 18:.1f}" y="{y + 14:.1f}" font-family="Arial" font-size="11" fill="#333">{label}</text>')
            parts.append(f'<rect x="{bar_left:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="14" fill="{color}"/>')
            parts.append(f'<text x="{bar_left + bar_w + 6:.1f}" y="{y + 12:.1f}" font-family="Arial" font-size="11" fill="#333">{value:.2f}</text>')
        axis_y = y0 + panel_h - 30
        threshold_x = bar_left + 0.7 * max_bar_w
        parts.append(f'<line x1="{bar_left:.1f}" y1="{axis_y:.1f}" x2="{bar_left + max_bar_w:.1f}" y2="{axis_y:.1f}" stroke="#65736f"/>')
        _add_x_ticks(parts, [0.0, 0.5, 0.7, 1.0], 0.0, 1.0, bar_left, axis_y, max_bar_w)
        parts.append(f'<line x1="{threshold_x:.1f}" y1="{y0 + 42:.1f}" x2="{threshold_x:.1f}" y2="{axis_y:.1f}" stroke="#c6912f" stroke-dasharray="3 4"/>')
        parts.append(f'<text x="{threshold_x + 5:.1f}" y="{axis_y - 8:.1f}" font-family="Arial" font-size="10" fill="#8a6a1f">0.70 mostly classical</text>')

    # Panel 2: observed primitive scores against classical predictions.
    x0, y0, _ = panels[1]
    needed = {"observed_score", "predicted_classical_score"}
    if residuals.empty or not needed.issubset(residuals.columns):
        parts.append(f'<text x="{x0 + 18:.1f}" y="{y0 + 58:.1f}" font-family="Arial" font-size="12" fill="#555">No observed/predicted score pairs available.</text>')
    else:
        x_values = pd.to_numeric(residuals["predicted_classical_score"], errors="coerce").to_numpy(dtype=float)
        y_values = pd.to_numeric(residuals["observed_score"], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(x_values) & np.isfinite(y_values)
        x_values = x_values[mask]
        y_values = y_values[mask]
        plot_left = x0 + 60
        plot_top = y0 + 48
        plot_w = panel_w - 88
        plot_h = panel_h - 90
        if x_values.size == 0:
            parts.append(f'<text x="{x0 + 18:.1f}" y="{y0 + 58:.1f}" font-family="Arial" font-size="12" fill="#555">No finite observed/predicted values.</text>')
        else:
            xmin, xmax = float(np.min(x_values)), float(np.max(x_values))
            ymin, ymax = float(np.min(y_values)), float(np.max(y_values))
            if xmin == xmax:
                xmin -= 1.0
                xmax += 1.0
            if ymin == ymax:
                ymin -= 1.0
                ymax += 1.0
            diag_min = max(xmin, ymin)
            diag_max = min(xmax, ymax)
            parts.append(f'<line x1="{plot_left:.1f}" y1="{plot_top + plot_h:.1f}" x2="{plot_left + plot_w:.1f}" y2="{plot_top + plot_h:.1f}" stroke="#65736f"/>')
            parts.append(f'<line x1="{plot_left:.1f}" y1="{plot_top:.1f}" x2="{plot_left:.1f}" y2="{plot_top + plot_h:.1f}" stroke="#65736f"/>')
            if diag_min < diag_max:
                x1 = _scale(diag_min, xmin, xmax, plot_left, plot_left + plot_w)
                y1 = _scale(diag_min, ymin, ymax, plot_top + plot_h, plot_top)
                x2 = _scale(diag_max, xmin, xmax, plot_left, plot_left + plot_w)
                y2 = _scale(diag_max, ymin, ymax, plot_top + plot_h, plot_top)
                parts.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#c6912f" stroke-dasharray="4 4"/>')
                parts.append(f'<text x="{x2 - 106:.1f}" y="{y2 - 6:.1f}" font-family="Arial" font-size="10" fill="#8a6a1f">y=x classical prediction</text>')
            stride = max(1, int(np.ceil(x_values.size / 450)))
            for x_val, y_val in zip(x_values[::stride], y_values[::stride]):
                x = _scale(float(x_val), xmin, xmax, plot_left, plot_left + plot_w)
                y = _scale(float(y_val), ymin, ymax, plot_top + plot_h, plot_top)
                parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.6" fill="#7b4b73" opacity="0.55"/>')
            x_mid = (xmin + xmax) / 2
            y_mid = (ymin + ymax) / 2
            _add_x_ticks(parts, [xmin, x_mid, xmax], xmin, xmax, plot_left, plot_top + plot_h, plot_w)
            _add_y_ticks(parts, [ymin, y_mid, ymax], ymin, ymax, plot_left, plot_top, plot_h)
            parts.append(f'<text x="{plot_left:.1f}" y="{y0 + panel_h - 14:.1f}" font-family="Arial" font-size="11" fill="#555">predicted primitive score, unitless</text>')
            parts.append(f'<text x="{x0 + 12:.1f}" y="{plot_top + 12:.1f}" font-family="Arial" font-size="11" fill="#555">observed score, unitless</text>')
            parts.append(f'<circle cx="{x0 + panel_w - 168:.1f}" cy="{y0 + 52:.1f}" r="4" fill="#7b4b73" opacity="0.55"/>')
            parts.append(f'<text x="{x0 + panel_w - 157:.1f}" y="{y0 + 56:.1f}" font-family="Arial" font-size="10" fill="#555">primitive score row</text>')
            parts.append(f'<line x1="{x0 + panel_w - 172:.1f}" y1="{y0 + 68:.1f}" x2="{x0 + panel_w - 162:.1f}" y2="{y0 + 68:.1f}" stroke="#c6912f" stroke-dasharray="4 4"/>')
            parts.append(f'<text x="{x0 + panel_w - 157:.1f}" y="{y0 + 72:.1f}" font-family="Arial" font-size="10" fill="#555">classical fit line</text>')

    # Panel 3: residual z-score versus matched-null z-score.
    x0, y0, _ = panels[2]
    needed = {"residual_zscore", "matched_null_zscore"}
    if residuals.empty or not needed.issubset(residuals.columns):
        parts.append(f'<text x="{x0 + 18:.1f}" y="{y0 + 58:.1f}" font-family="Arial" font-size="12" fill="#555">No matched-null evidence available.</text>')
    else:
        x_values = pd.to_numeric(residuals["residual_zscore"], errors="coerce").to_numpy(dtype=float)
        y_values = pd.to_numeric(residuals["matched_null_zscore"], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(x_values) & np.isfinite(y_values)
        x_values = x_values[mask]
        y_values = y_values[mask]
        plot_left = x0 + 60
        plot_top = y0 + 48
        plot_w = panel_w - 88
        plot_h = panel_h - 92
        if x_values.size == 0:
            parts.append(f'<text x="{x0 + 18:.1f}" y="{y0 + 58:.1f}" font-family="Arial" font-size="12" fill="#555">No finite residual/null z-score pairs.</text>')
        else:
            zmax = max(2.5, float(np.nanmax(np.abs(np.r_[x_values, y_values]))))
            parts.append(f'<line x1="{plot_left:.1f}" y1="{plot_top + plot_h / 2:.1f}" x2="{plot_left + plot_w:.1f}" y2="{plot_top + plot_h / 2:.1f}" stroke="#c8d1cc"/>')
            parts.append(f'<line x1="{plot_left + plot_w / 2:.1f}" y1="{plot_top:.1f}" x2="{plot_left + plot_w / 2:.1f}" y2="{plot_top + plot_h:.1f}" stroke="#c8d1cc"/>')
            for threshold in (-2.0, 2.0):
                x = _scale(threshold, -zmax, zmax, plot_left, plot_left + plot_w)
                y = _scale(threshold, -zmax, zmax, plot_top + plot_h, plot_top)
                parts.append(f'<line x1="{x:.1f}" y1="{plot_top:.1f}" x2="{x:.1f}" y2="{plot_top + plot_h:.1f}" stroke="#c6912f" stroke-dasharray="3 4" opacity="0.65"/>')
                parts.append(f'<line x1="{plot_left:.1f}" y1="{y:.1f}" x2="{plot_left + plot_w:.1f}" y2="{y:.1f}" stroke="#c6912f" stroke-dasharray="3 4" opacity="0.65"/>')
            x_threshold = _scale(2.0, -zmax, zmax, plot_left, plot_left + plot_w)
            y_threshold = _scale(2.0, -zmax, zmax, plot_top + plot_h, plot_top)
            parts.append(f'<text x="{x_threshold + 5:.1f}" y="{plot_top + 12:.1f}" font-family="Arial" font-size="10" fill="#8a6a1f">residual z = 2</text>')
            parts.append(f'<text x="{plot_left + 5:.1f}" y="{y_threshold - 6:.1f}" font-family="Arial" font-size="10" fill="#8a6a1f">null z = 2</text>')
            parts.append(f'<line x1="{plot_left:.1f}" y1="{plot_top + plot_h:.1f}" x2="{plot_left + plot_w:.1f}" y2="{plot_top + plot_h:.1f}" stroke="#65736f"/>')
            parts.append(f'<line x1="{plot_left:.1f}" y1="{plot_top:.1f}" x2="{plot_left:.1f}" y2="{plot_top + plot_h:.1f}" stroke="#65736f"/>')
            stride = max(1, int(np.ceil(x_values.size / 450)))
            for x_val, y_val in zip(x_values[::stride], y_values[::stride]):
                x = _scale(float(x_val), -zmax, zmax, plot_left, plot_left + plot_w)
                y = _scale(float(y_val), -zmax, zmax, plot_top + plot_h, plot_top)
                strong = x_val >= 2.0 and y_val >= 2.0
                color = "#7b4b73" if strong else "#426b69"
                parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{3.8 if strong else 2.6}" fill="{color}" opacity="0.62"/>')
            tick_max = round(zmax, 1)
            _add_x_ticks(parts, [-tick_max, -2.0, 0.0, 2.0, tick_max], -zmax, zmax, plot_left, plot_top + plot_h, plot_w)
            _add_y_ticks(parts, [-tick_max, -2.0, 0.0, 2.0, tick_max], -zmax, zmax, plot_left, plot_top, plot_h)
            parts.append(f'<text x="{plot_left:.1f}" y="{y0 + panel_h - 14:.1f}" font-family="Arial" font-size="11" fill="#555">residual z-score after classical controls</text>')
            parts.append(f'<text x="{x0 + 12:.1f}" y="{plot_top + 12:.1f}" font-family="Arial" font-size="11" fill="#555">null z</text>')
            parts.append(f'<circle cx="{x0 + panel_w - 170:.1f}" cy="{y0 + 52:.1f}" r="4" fill="#7b4b73" opacity="0.62"/>')
            parts.append(f'<text x="{x0 + panel_w - 158:.1f}" y="{y0 + 56:.1f}" font-family="Arial" font-size="10" fill="#555">residual and null z >= 2</text>')
            parts.append(f'<circle cx="{x0 + panel_w - 170:.1f}" cy="{y0 + 68:.1f}" r="4" fill="#426b69" opacity="0.62"/>')
            parts.append(f'<text x="{x0 + panel_w - 158:.1f}" y="{y0 + 72:.1f}" font-family="Arial" font-size="10" fill="#555">background/review point</text>')

    # Panel 4: strongest post-control candidate rows.
    x0, y0, _ = panels[3]
    if residuals.empty or not {"region_id", "primitive", "residual_zscore"}.issubset(residuals.columns):
        parts.append(f'<text x="{x0 + 18:.1f}" y="{y0 + 58:.1f}" font-family="Arial" font-size="12" fill="#555">No residual candidate rows available.</text>')
    else:
        top = residuals.copy()
        top["residual_zscore"] = pd.to_numeric(top["residual_zscore"], errors="coerce")
        if "matched_null_zscore" in top.columns:
            top["matched_null_zscore"] = pd.to_numeric(top["matched_null_zscore"], errors="coerce")
        else:
            top["matched_null_zscore"] = np.nan
        residual_part = top["residual_zscore"].clip(lower=0.0)
        null_part = top["matched_null_zscore"].clip(lower=0.0)
        top["support_score"] = residual_part.fillna(0.0) + 0.5 * null_part.fillna(0.0)
        top = top.dropna(subset=["residual_zscore"]).sort_values("support_score", ascending=False).head(8)
        if not labels.empty and {"region_id", "primitive_class"}.issubset(labels.columns):
            label_cols = ["region_id", "primitive_class"]
            if "primitive_score_name" in labels.columns:
                label_cols.append("primitive_score_name")
                top = top.merge(
                    labels[label_cols].drop_duplicates(label_cols),
                    left_on=["region_id", "primitive"],
                    right_on=["region_id", "primitive_score_name"],
                    how="left",
                )
            else:
                top = top.merge(labels[label_cols], on="region_id", how="left")
        max_score = max(1e-9, float(top["support_score"].max())) if not top.empty else 1.0
        bar_left = x0 + 218
        bar_top = y0 + 46
        row_h = max(22, (panel_h - 74) / max(1, len(top)))
        max_bar_w = panel_w - 265
        for idx, row in enumerate(top.itertuples(index=False)):
            y = bar_top + idx * row_h
            primitive = getattr(row, "primitive_class", None) or getattr(row, "primitive")
            label = escape(_primitive_label(primitive)[:28])
            score = float(getattr(row, "support_score"))
            rz = float(getattr(row, "residual_zscore"))
            nz = float(getattr(row, "matched_null_zscore"))
            bar_w = (score / max_score) * max_bar_w
            parts.append(f'<text x="{x0 + 18:.1f}" y="{y + 14:.1f}" font-family="Arial" font-size="11" fill="#333">{label}</text>')
            parts.append(f'<rect x="{bar_left:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="14" fill="#7b4b73"/>')
            parts.append(f'<text x="{bar_left + bar_w + 6:.1f}" y="{y + 12:.1f}" font-family="Arial" font-size="11" fill="#333">z {rz:.1f} / n {nz:.1f}</text>')
        axis_y = y0 + panel_h - 30
        parts.append(f'<line x1="{bar_left:.1f}" y1="{axis_y:.1f}" x2="{bar_left + max_bar_w:.1f}" y2="{axis_y:.1f}" stroke="#65736f"/>')
        _add_x_ticks(parts, [0.0, round(max_score / 2, 1), round(max_score, 1)], 0.0, max_score, bar_left, axis_y, max_bar_w)
        if max_score >= 3.0:
            threshold_x = _scale(3.0, 0.0, max_score, bar_left, bar_left + max_bar_w)
            parts.append(f'<line x1="{threshold_x:.1f}" y1="{y0 + 42:.1f}" x2="{threshold_x:.1f}" y2="{axis_y:.1f}" stroke="#c6912f" stroke-dasharray="3 4"/>')
            parts.append(f'<text x="{threshold_x + 5:.1f}" y="{axis_y - 8:.1f}" font-family="Arial" font-size="10" fill="#8a6a1f">review guide 3</text>')
        parts.append(f'<text x="{bar_left:.1f}" y="{y0 + panel_h - 12:.1f}" font-family="Arial" font-size="11" fill="#555">support score = residual z + 0.5 x matched-null z</text>')

    path.write_text(_svg_frame(width, height, title, "\n".join(parts), show_title=False), encoding="utf-8")


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

    classical_path = out / "classical_control_multipanel.svg"
    _write_classical_control_multipanel(residuals, labels, classical_path)
    paths["classical_control_multipanel"] = classical_path

    hist_path = out / "residual_score_histogram.svg"
    _write_histogram(_finite(residuals.get("residual_zscore", pd.Series(dtype=float))), hist_path, "Residual Score Histogram")
    paths["residual_score_histogram"] = hist_path

    scatter_path = out / "observed_vs_predicted_classical_score.svg"
    predicted = pd.to_numeric(residuals.get("predicted_classical_score", pd.Series(dtype=float)), errors="coerce").to_numpy(dtype=float)
    observed = pd.to_numeric(residuals.get("observed_score", pd.Series(dtype=float)), errors="coerce").to_numpy(dtype=float)
    _write_scatter(predicted, observed, scatter_path, "Observed vs Predicted Classical Score")
    paths["observed_vs_predicted_classical_score"] = scatter_path
    return paths
