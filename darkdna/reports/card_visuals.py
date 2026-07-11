"""Per-card visual summaries and standalone card bundles."""

from __future__ import annotations

import json
import re
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.io.bed import write_bed


PALETTE = [
    "#426b69",
    "#7b4b73",
    "#9b5f3f",
    "#3f6b8f",
    "#7a6b2c",
    "#6b5f95",
    "#587842",
    "#8a4f5f",
]


STANDALONE_CARD_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{title}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 32px; color: #1d2528; background: #fbfbf8; }}
    main {{ max-width: 1120px; }}
    h1, h2 {{ color: #223c3b; }}
    a {{ color: #315f8f; }}
    .caveat {{ padding: 12px; background: #fff5d9; border-left: 4px solid #c6912f; }}
    .card-section {{ border: 1px solid #d7ded8; border-radius: 6px; padding: 12px; margin: 14px 0; background: white; }}
    .card-visual svg {{ max-width: 100%; height: auto; display: block; margin: 10px 0 18px; }}
    table {{ border-collapse: collapse; width: 100%; margin: 12px 0; background: white; }}
    th, td {{ border: 1px solid #d7ded8; padding: 6px 8px; font-size: 13px; vertical-align: top; text-align: left; }}
    th {{ background: #e7efea; width: 220px; }}
  </style>
</head>
<body>
<main>
  <p><a href="../../darkdna_report.html">Back to report</a></p>
  <h1>{heading}</h1>
  <div class="card-visual">{visual_html}</div>
  <p class="caveat">{candidate_statement}</p>
  <div class="card-section">
    <h2>Candidate Interpretation</h2>
    <table><tbody>
      {summary_rows}
    </tbody></table>
  </div>
  <div class="card-section">
    <h2>Card Files</h2>
    <table><tbody>
      <tr><th>Focal BED</th><td><a href="region.bed">region.bed</a></td></tr>
      <tr><th>Plotted span BED</th><td><a href="view_window.bed">view_window.bed</a></td></tr>
      <tr><th>Visible primitive intervals</th><td><a href="primitive_intervals.bed">primitive_intervals.bed</a></td></tr>
      <tr><th>Visible candidate loci</th><td><a href="candidate_loci.bed">candidate_loci.bed</a></td></tr>
      <tr><th>Card JSON</th><td><a href="card.json">card.json</a></td></tr>
    </tbody></table>
  </div>
</main>
</body>
</html>
"""


def augment_cards_with_visuals(
    cards: list[dict],
    windows: pd.DataFrame,
    labels: pd.DataFrame,
    residuals: pd.DataFrame,
    loci: pd.DataFrame,
    outdir: str | Path,
) -> list[dict]:
    """Attach SVG visuals to cards and write standalone card bundles."""

    out = Path(outdir)
    cards_dir = out / "cards"
    cards_dir.mkdir(parents=True, exist_ok=True)
    used_slugs: set[str] = set()
    augmented: list[dict] = []
    for card in cards:
        copy = dict(card)
        focus = _focus_interval(copy, windows)
        primitive = str(copy.get("primitive_class", "candidate"))
        slug = _unique_slug(_slug(f"{focus['name']}_{primitive}"), used_slugs)
        folder = cards_dir / slug
        view = _view_interval(focus, windows)
        intervals = _primitive_intervals_for_view(windows, labels, residuals, view)
        view_loci = _loci_for_view(loci, view)
        visual_html = build_card_visual(copy, windows, labels, residuals, loci, focus=focus, view=view, intervals=intervals, view_loci=view_loci)
        _write_card_bundle(folder, copy, focus, view, intervals, view_loci, visual_html)
        copy["card_visual_html"] = visual_html
        copy["card_href"] = f"cards/{slug}/index.html"
        copy["card_region_bed_href"] = f"cards/{slug}/region.bed"
        copy["card_primitive_intervals_bed_href"] = f"cards/{slug}/primitive_intervals.bed"
        augmented.append(copy)
    return augmented


def build_card_visual(
    card: dict,
    windows: pd.DataFrame,
    labels: pd.DataFrame,
    residuals: pd.DataFrame,
    loci: pd.DataFrame,
    *,
    focus: dict | None = None,
    view: dict | None = None,
    intervals: pd.DataFrame | None = None,
    view_loci: pd.DataFrame | None = None,
) -> str:
    """Return an inline SVG for genome position plus primitive behavior."""

    focus = focus or _focus_interval(card, windows)
    view = view or _view_interval(focus, windows)
    intervals = intervals if intervals is not None else _primitive_intervals_for_view(windows, labels, residuals, view)
    view_loci = view_loci if view_loci is not None else _loci_for_view(loci, view)

    primitive = str(card.get("primitive_class", "candidate"))
    lanes = _lane_order(intervals, primitive)
    lane_count = max(1, len(lanes))
    width = 1040
    left = 230
    right = 34
    plot_w = width - left - right
    y_global = 54
    y_local = 116
    y_loci = 152
    y_lanes = 194
    lane_h = 26
    model_y = y_lanes + lane_count * lane_h + 54
    height = model_y + 178
    focus_x = _scale(focus["start"], view["start"], view["end"], left, plot_w)
    focus_w = max(2.0, _scale(focus["end"], view["start"], view["end"], left, plot_w) - focus_x)
    genome_x = _scale((focus["start"] + focus["end"]) / 2, 0, view["chrom_size"], left, plot_w)

    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{escape(_visual_alt(card, focus, intervals))}">',
        "<title>Candidate card genomic context and primitive hypothesis</title>",
        f'<rect width="{width}" height="{height}" fill="#fbfbf8"/>',
        '<text x="24" y="28" font-family="Arial" font-size="15" fill="#223c3b">Genome context</text>',
        '<text x="24" y="45" font-family="Arial" font-size="11" fill="#555">observed chromosome span and local candidate neighborhood</text>',
        f'<text x="24" y="{y_global + 4}" font-family="Arial" font-size="12" fill="#333">{escape(str(focus["chrom"]))}</text>',
        f'<line x1="{left}" y1="{y_global}" x2="{left + plot_w}" y2="{y_global}" stroke="#65736f" stroke-width="2"/>',
        f'<circle cx="{genome_x:.1f}" cy="{y_global}" r="6" fill="#c6912f" stroke="#8a6a1f"/>',
        f'<text x="{left}" y="{y_global + 22}" font-family="Arial" font-size="10" fill="#555">0 bp</text>',
        f'<text x="{left + plot_w - 88}" y="{y_global + 22}" font-family="Arial" font-size="10" fill="#555">{_fmt_bp(view["chrom_size"])} observed</text>',
        f'<text x="24" y="{y_local + 4}" font-family="Arial" font-size="12" fill="#333">local view</text>',
        f'<rect x="{left}" y="{y_local - 8}" width="{plot_w}" height="18" fill="#eef3f0" stroke="#d7ded8"/>',
        f'<rect x="{focus_x:.1f}" y="{y_local - 12}" width="{focus_w:.1f}" height="26" fill="#fff5d9" stroke="#c6912f"/>',
        f'<text x="{min(left + plot_w - 120, focus_x + focus_w + 6):.1f}" y="{y_local - 17}" font-family="Arial" font-size="10" fill="#8a6a1f">focal window</text>',
    ]
    _append_axis_ticks(parts, view["start"], view["end"], left, y_local + 16, plot_w)
    parts.extend(
        [
            f'<text x="24" y="{y_loci + 4}" font-family="Arial" font-size="12" fill="#333">merged loci</text>',
            f'<line x1="{left}" y1="{y_loci}" x2="{left + plot_w}" y2="{y_loci}" stroke="#d7ded8"/>',
        ]
    )
    if view_loci.empty:
        parts.append(f'<text x="{left}" y="{y_loci + 4}" font-family="Arial" font-size="11" fill="#777">No merged candidate locus overlaps this plotted span.</text>')
    else:
        for row in view_loci.itertuples(index=False):
            x, w = _scaled_interval(row.start, row.end, view, left, plot_w)
            name = escape(str(getattr(row, "primitive_class", "candidate")))
            score = _fmt_float(getattr(row, "max_primitive_confidence", np.nan), 3)
            parts.append(
                f'<rect x="{x:.1f}" y="{y_loci - 8}" width="{w:.1f}" height="16" fill="none" stroke="#7b4b73" stroke-width="1.5">'
                f"<title>{name}; confidence {score}</title></rect>"
            )
    parts.append(f'<text x="24" y="{y_lanes - 13}" font-family="Arial" font-size="15" fill="#223c3b">Primitive signals in this region</text>')
    if not lanes:
        parts.append(f'<text x="{left}" y="{y_lanes + 8}" font-family="Arial" font-size="12" fill="#555">No non-no_call primitive intervals overlap this view.</text>')
    else:
        for lane_idx, lane in enumerate(lanes):
            y = y_lanes + lane_idx * lane_h
            color = _color_for_primitive(lane)
            is_focal_lane = lane == primitive
            if is_focal_lane:
                parts.append(f'<rect x="{left - 8}" y="{y - 13}" width="{plot_w + 16}" height="{lane_h - 3}" fill="#fff5d9" opacity="0.55"/>')
            parts.append(f'<text x="24" y="{y + 4}" font-family="Arial" font-size="11" fill="#333">{escape(_primitive_label(lane, 27))}</text>')
            parts.append(f'<line x1="{left}" y1="{y}" x2="{left + plot_w}" y2="{y}" stroke="#e2e7e3"/>')
            subset = intervals[intervals["primitive_class"].astype(str) == lane]
            for row in subset.itertuples(index=False):
                x, w = _scaled_interval(row.start, row.end, view, left, plot_w)
                confidence = _finite_float(getattr(row, "primitive_confidence", 0.0), default=0.0)
                opacity = 0.24 + 0.66 * min(1.0, max(0.0, confidence))
                is_focal = str(getattr(row, "region_id", "")) == str(focus["name"])
                stroke = "#1d2528" if is_focal else color
                stroke_w = 2 if is_focal else 0.8
                tooltip = _interval_tooltip(row)
                parts.append(
                    f'<rect x="{x:.1f}" y="{y - 8}" width="{w:.1f}" height="16" fill="{color}" fill-opacity="{opacity:.2f}" '
                    f'stroke="{stroke}" stroke-width="{stroke_w}"><title>{escape(tooltip)}</title></rect>'
                )
    parts.extend(_primitive_model_panel(card, primitive, model_y, width))
    parts.append("</svg>")
    return "\n".join(parts)


def _write_card_bundle(
    folder: Path,
    card: dict,
    focus: dict,
    view: dict,
    intervals: pd.DataFrame,
    view_loci: pd.DataFrame,
    visual_html: str,
) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "card.json").write_text(json.dumps(card, indent=2, default=str), encoding="utf-8")
    _write_bed_files(folder, card, focus, view, intervals, view_loci)
    html = STANDALONE_CARD_TEMPLATE.format(
        title=escape(f"{focus['name']} - {card.get('primitive_class', 'candidate')}"),
        heading=escape(f"{focus['name']} - {card.get('primitive_class', 'candidate')}"),
        visual_html=visual_html,
        candidate_statement=escape(str(card.get("candidate_statement", "This is a sequence-derived candidate hypothesis."))),
        summary_rows=_summary_rows(card),
    )
    (folder / "index.html").write_text(html, encoding="utf-8")


def _write_bed_files(folder: Path, card: dict, focus: dict, view: dict, intervals: pd.DataFrame, view_loci: pd.DataFrame) -> None:
    confidence = int(round(1000 * _finite_float(card.get("primitive_confidence", 0.0), default=0.0)))
    region = pd.DataFrame(
        [
            {
                "chrom": focus["chrom"],
                "start": focus["start"],
                "end": focus["end"],
                "name": f"{card.get('primitive_class', 'candidate')}|{focus['name']}",
                "score": max(0, min(1000, confidence)),
                "strand": ".",
            }
        ]
    )
    write_bed(region, folder / "region.bed", columns=["chrom", "start", "end", "name", "score", "strand"])
    view_row = pd.DataFrame(
        [{"chrom": view["chrom"], "start": view["start"], "end": view["end"], "name": f"view|{focus['name']}", "score": 0, "strand": "."}]
    )
    write_bed(view_row, folder / "view_window.bed", columns=["chrom", "start", "end", "name", "score", "strand"])

    primitive_bed = intervals.copy()
    if primitive_bed.empty:
        (folder / "primitive_intervals.bed").write_text("", encoding="utf-8")
    else:
        primitive_bed["name"] = primitive_bed["primitive_class"].astype(str) + "|" + primitive_bed["region_id"].astype(str)
        primitive_bed["score"] = (
            pd.to_numeric(primitive_bed.get("primitive_confidence", 0.0), errors="coerce").fillna(0.0).clip(0.0, 1.0) * 1000
        ).round().astype(int)
        primitive_bed["strand"] = "."
        write_bed(primitive_bed, folder / "primitive_intervals.bed", columns=["chrom", "start", "end", "name", "score", "strand"])

    loci_bed = view_loci.copy()
    if loci_bed.empty:
        (folder / "candidate_loci.bed").write_text("", encoding="utf-8")
    else:
        loci_bed["name"] = loci_bed.get("locus_id", loci_bed["primitive_class"].astype(str))
        loci_bed["score"] = (
            pd.to_numeric(loci_bed.get("max_primitive_confidence", 0.0), errors="coerce").fillna(0.0).clip(0.0, 1.0) * 1000
        ).round().astype(int)
        loci_bed["strand"] = "."
        write_bed(loci_bed, folder / "candidate_loci.bed", columns=["chrom", "start", "end", "name", "score", "strand"])


def _focus_interval(card: dict, windows: pd.DataFrame) -> dict:
    region_id = str(card.get("region_id", "candidate"))
    if windows is not None and not windows.empty and "region_id" in windows.columns:
        subset = windows[windows["region_id"].astype(str) == region_id]
        if not subset.empty:
            row = subset.iloc[0]
            return {
                "chrom": str(row.get("chrom", "region")),
                "start": int(row.get("start", 0)),
                "end": int(row.get("end", 1)),
                "name": region_id,
            }
    parsed = _parse_coordinates(str(card.get("coordinates", "")))
    if parsed:
        chrom, start, end = parsed
        return {"chrom": chrom, "start": start, "end": end, "name": region_id}
    return {"chrom": "region", "start": 0, "end": max(1, int(card.get("window_size", 1) or 1)), "name": region_id}


def _parse_coordinates(value: str) -> tuple[str, int, int] | None:
    match = re.match(r"^(.+):(\d+)-(\d+)$", value.strip())
    if not match:
        return None
    return match.group(1), int(match.group(2)), int(match.group(3))


def _view_interval(focus: dict, windows: pd.DataFrame) -> dict:
    chrom = str(focus["chrom"])
    chrom_size = max(int(focus["end"]), _chrom_size(chrom, windows))
    length = max(1, int(focus["end"]) - int(focus["start"]))
    flank = max(5000, length * 4)
    start = max(0, int(focus["start"]) - flank)
    end = min(chrom_size, int(focus["end"]) + flank)
    if end <= start:
        end = start + length
    return {"chrom": chrom, "start": start, "end": end, "chrom_size": chrom_size}


def _chrom_size(chrom: str, windows: pd.DataFrame) -> int:
    if windows is None or windows.empty or not {"chrom", "end"}.issubset(windows.columns):
        return 1
    subset = windows[windows["chrom"].astype(str) == str(chrom)]
    if subset.empty:
        return 1
    return int(pd.to_numeric(subset["end"], errors="coerce").fillna(0).max())


def _primitive_intervals_for_view(windows: pd.DataFrame, labels: pd.DataFrame, residuals: pd.DataFrame, view: dict) -> pd.DataFrame:
    needed = {"region_id", "chrom", "start", "end"}
    if windows is None or labels is None or windows.empty or labels.empty or not needed.issubset(windows.columns) or "region_id" not in labels.columns:
        return pd.DataFrame(columns=["region_id", "chrom", "start", "end", "primitive_class", "primitive_confidence"])
    window_cols = [col for col in ["region_id", "chrom", "start", "end", "window_size", "artifact_risk_flags"] if col in windows.columns]
    label_cols = [col for col in ["region_id", "primitive_class", "primitive_score_name", "primitive_confidence", "top_supporting_features"] if col in labels.columns]
    merged = windows[window_cols].merge(labels[label_cols], on="region_id", how="inner")
    if merged.empty or "primitive_class" not in merged.columns:
        return pd.DataFrame(columns=merged.columns)
    merged["start"] = pd.to_numeric(merged["start"], errors="coerce")
    merged["end"] = pd.to_numeric(merged["end"], errors="coerce")
    merged = merged.dropna(subset=["start", "end"])
    merged["start"] = merged["start"].astype(int)
    merged["end"] = merged["end"].astype(int)
    merged = merged[
        (merged["chrom"].astype(str) == str(view["chrom"]))
        & (merged["start"] < int(view["end"]))
        & (merged["end"] > int(view["start"]))
        & (merged["primitive_class"].fillna("").astype(str).ne(""))
        & (merged["primitive_class"].fillna("").astype(str).ne("no_call"))
    ].copy()
    if merged.empty:
        return merged
    if "primitive_score_name" not in merged.columns:
        merged["primitive_score_name"] = merged["primitive_class"].astype(str) + "_score"
    if residuals is not None and not residuals.empty and {"region_id", "primitive"}.issubset(residuals.columns):
        res_cols = [col for col in ["region_id", "primitive", "residual_zscore", "matched_null_zscore", "observed_score", "empirical_p_value"] if col in residuals.columns]
        merged = merged.merge(
            residuals[res_cols],
            left_on=["region_id", "primitive_score_name"],
            right_on=["region_id", "primitive"],
            how="left",
        )
    merged["primitive_confidence"] = pd.to_numeric(merged.get("primitive_confidence", 0.0), errors="coerce").fillna(0.0)
    return merged.sort_values(["primitive_class", "start", "end"], kind="mergesort")


def _loci_for_view(loci: pd.DataFrame, view: dict) -> pd.DataFrame:
    if loci is None or loci.empty or not {"chrom", "start", "end"}.issubset(loci.columns):
        return pd.DataFrame(columns=[] if loci is None else loci.columns)
    subset = loci.copy()
    subset["start"] = pd.to_numeric(subset["start"], errors="coerce")
    subset["end"] = pd.to_numeric(subset["end"], errors="coerce")
    subset = subset.dropna(subset=["start", "end"])
    subset["start"] = subset["start"].astype(int)
    subset["end"] = subset["end"].astype(int)
    return subset[
        (subset["chrom"].astype(str) == str(view["chrom"]))
        & (subset["start"] < int(view["end"]))
        & (subset["end"] > int(view["start"]))
    ].copy()


def _lane_order(intervals: pd.DataFrame, focal_primitive: str) -> list[str]:
    if intervals is None or intervals.empty or "primitive_class" not in intervals.columns:
        return []
    counts = intervals.groupby("primitive_class")["primitive_confidence"].max().sort_values(ascending=False)
    lanes = [str(value) for value in counts.index]
    if focal_primitive in lanes:
        lanes.remove(focal_primitive)
        lanes.insert(0, focal_primitive)
    return lanes


def _primitive_model_panel(card: dict, primitive: str, y: float, width: int) -> list[str]:
    model = _model_text(card)
    x0 = 24
    x1 = 360
    x2 = 690
    parts = [
        f'<line x1="24" y1="{y - 24}" x2="{width - 24}" y2="{y - 24}" stroke="#d7ded8"/>',
        f'<text x="24" y="{y}" font-family="Arial" font-size="15" fill="#223c3b">Hypothesized primitive behavior</text>',
        f'<text x="24" y="{y + 18}" font-family="Arial" font-size="11" fill="#555">sequence proxy -> candidate behavior -> validation readout</text>',
        f'<text x="{x0}" y="{y + 50}" font-family="Arial" font-size="12" fill="#333">Measured proxy</text>',
        f'<text x="{x1}" y="{y + 50}" font-family="Arial" font-size="12" fill="#333">Primitive sketch</text>',
        f'<text x="{x2}" y="{y + 50}" font-family="Arial" font-size="12" fill="#333">Behavior to test</text>',
        f'<line x1="{x0 + 210}" y1="{y + 94}" x2="{x1 - 24}" y2="{y + 94}" stroke="#65736f" marker-end="url(#arrowhead)"/>',
        f'<line x1="{x1 + 240}" y1="{y + 94}" x2="{x2 - 24}" y2="{y + 94}" stroke="#65736f" marker-end="url(#arrowhead)"/>',
        '<defs><marker id="arrowhead" markerWidth="8" markerHeight="8" refX="7" refY="3.5" orient="auto"><polygon points="0 0, 8 3.5, 0 7" fill="#65736f"/></marker></defs>',
    ]
    parts.extend(_wrapped_svg_text(model["feature"], x0, y + 72, 44, max_lines=3, fill="#333"))
    parts.extend(_primitive_symbol(primitive, x1, y + 64))
    parts.extend(_wrapped_svg_text(model["behavior"], x2, y + 72, 43, max_lines=3, fill="#333"))
    assay = model["assay"]
    if assay:
        parts.append(f'<text x="{x2}" y="{y + 148}" font-family="Arial" font-size="11" fill="#555">assay: {escape(assay[:70])}</text>')
    return parts


def _primitive_symbol(primitive: str, x: float, y: float) -> list[str]:
    key = str(primitive)
    parts = [f'<rect x="{x - 10}" y="{y - 8}" width="260" height="96" fill="#ffffff" stroke="#d7ded8" rx="4"/>']
    if key in {"decoherence_boundary_candidate", "sequence_regime_boundary_candidate"}:
        parts.extend(
            [
                f'<rect x="{x + 12}" y="{y + 22}" width="92" height="24" fill="#426b69" opacity="0.25"/>',
                f'<rect x="{x + 132}" y="{y + 22}" width="92" height="24" fill="#7b4b73" opacity="0.25"/>',
                f'<line x1="{x + 118}" y1="{y + 4}" x2="{x + 118}" y2="{y + 72}" stroke="#c6912f" stroke-width="3"/>',
                f'<circle cx="{x + 40}" cy="{y + 14}" r="3" fill="#426b69"/><circle cx="{x + 66}" cy="{y + 58}" r="3" fill="#426b69"/>',
                f'<circle cx="{x + 154}" cy="{y + 20}" r="2" fill="#7b4b73"/><circle cx="{x + 184}" cy="{y + 54}" r="2" fill="#7b4b73"/>',
                f'<path d="M {x + 52} {y + 35} C {x + 88} {y + 12}, {x + 102} {y + 20}, {x + 114} {y + 34}" fill="none" stroke="#65736f"/>',
                f'<path d="M {x + 122} {y + 36} C {x + 144} {y + 44}, {x + 168} {y + 44}, {x + 198} {y + 32}" fill="none" stroke="#65736f" stroke-dasharray="4 4"/>',
            ]
        )
    elif key == "criticality_tuner_candidate":
        parts.extend(
            [
                f'<line x1="{x + 24}" y1="{y + 66}" x2="{x + 220}" y2="{y + 66}" stroke="#65736f"/>',
                f'<line x1="{x + 24}" y1="{y + 66}" x2="{x + 24}" y2="{y + 12}" stroke="#65736f"/>',
                f'<path d="M {x + 26} {y + 62} C {x + 78} {y + 62}, {x + 102} {y + 58}, {x + 118} {y + 38} C {x + 134} {y + 16}, {x + 166} {y + 14}, {x + 218} {y + 14}" fill="none" stroke="#7b4b73" stroke-width="3"/>',
                f'<line x1="{x + 118}" y1="{y + 10}" x2="{x + 118}" y2="{y + 70}" stroke="#c6912f" stroke-dasharray="4 3"/>',
            ]
        )
    elif key == "replication_instability_candidate":
        parts.extend(
            [
                f'<line x1="{x + 26}" y1="{y + 40}" x2="{x + 218}" y2="{y + 40}" stroke="#426b69" stroke-width="3"/>',
                f'<path d="M {x + 92} {y + 40} L {x + 122} {y + 24} M {x + 92} {y + 40} L {x + 122} {y + 56}" stroke="#7b4b73" stroke-width="3" fill="none"/>',
                f'<rect x="{x + 132}" y="{y + 28}" width="30" height="24" fill="#fff5d9" stroke="#c6912f"/>',
                f'<path d="M {x + 146} {y + 28} C {x + 132} {y + 10}, {x + 168} {y + 10}, {x + 154} {y + 28}" fill="none" stroke="#9b5f3f"/>',
            ]
        )
    elif key == "resonant_pulse_decoder_candidate":
        parts.extend(
            [
                f'<path d="M {x + 24} {y + 50} L {x + 42} {y + 50} L {x + 42} {y + 22} L {x + 52} {y + 22} L {x + 52} {y + 50} L {x + 82} {y + 50} L {x + 82} {y + 22} L {x + 92} {y + 22} L {x + 92} {y + 50} L {x + 122} {y + 50} L {x + 122} {y + 22} L {x + 132} {y + 22} L {x + 132} {y + 50}" fill="none" stroke="#426b69" stroke-width="2"/>',
                f'<path d="M {x + 24} {y + 68} C {x + 60} {y + 68}, {x + 74} {y + 34}, {x + 104} {y + 34} C {x + 138} {y + 34}, {x + 154} {y + 68}, {x + 218} {y + 68}" fill="none" stroke="#7b4b73" stroke-width="3"/>',
            ]
        )
    elif key == "hysteresis_candidate":
        parts.extend(
            [
                f'<line x1="{x + 36}" y1="{y + 66}" x2="{x + 214}" y2="{y + 66}" stroke="#65736f"/>',
                f'<line x1="{x + 36}" y1="{y + 66}" x2="{x + 36}" y2="{y + 14}" stroke="#65736f"/>',
                f'<path d="M {x + 46} {y + 56} C {x + 86} {y + 10}, {x + 168} {y + 10}, {x + 204} {y + 56} C {x + 160} {y + 38}, {x + 96} {y + 38}, {x + 46} {y + 56}" fill="none" stroke="#7b4b73" stroke-width="3"/>',
            ]
        )
    elif key == "negative_space_element_candidate":
        parts.extend(
            [
                f'<rect x="{x + 24}" y="{y + 30}" width="58" height="20" fill="#426b69" opacity="0.35"/>',
                f'<rect x="{x + 162}" y="{y + 30}" width="58" height="20" fill="#426b69" opacity="0.35"/>',
                f'<rect x="{x + 92}" y="{y + 20}" width="58" height="40" fill="none" stroke="#c6912f" stroke-dasharray="4 4"/>',
                f'<text x="{x + 100}" y="{y + 45}" font-family="Arial" font-size="11" fill="#8a6a1f">void</text>',
            ]
        )
    elif key == "non_B_DNA_physical_susceptibility_candidate":
        parts.extend(
            [
                f'<line x1="{x + 26}" y1="{y + 54}" x2="{x + 214}" y2="{y + 54}" stroke="#426b69" stroke-width="2"/>',
                f'<path d="M {x + 74} {y + 54} C {x + 78} {y + 10}, {x + 128} {y + 10}, {x + 132} {y + 54}" fill="none" stroke="#7b4b73" stroke-width="3"/>',
                f'<circle cx="{x + 98}" cy="{y + 24}" r="5" fill="#c6912f"/><circle cx="{x + 112}" cy="{y + 24}" r="5" fill="#c6912f"/>',
                f'<path d="M {x + 154} {y + 18} L {x + 142} {y + 44} L {x + 162} {y + 44} L {x + 150} {y + 70}" fill="none" stroke="#9b5f3f" stroke-width="2"/>',
            ]
        )
    elif key == "fractal_scaffold_candidate":
        parts.extend(
            [
                f'<rect x="{x + 42}" y="{y + 22}" width="158" height="48" fill="none" stroke="#426b69" stroke-width="2"/>',
                f'<rect x="{x + 72}" y="{y + 32}" width="44" height="28" fill="none" stroke="#7b4b73" stroke-width="2"/>',
                f'<rect x="{x + 132}" y="{y + 32}" width="44" height="28" fill="none" stroke="#7b4b73" stroke-width="2"/>',
                f'<path d="M {x + 42} {y + 72} C {x + 82} {y + 92}, {x + 158} {y + 92}, {x + 200} {y + 72}" fill="none" stroke="#c6912f"/>',
            ]
        )
    elif key == "constraint_grammar_region_candidate":
        bases = ["A", "C", "G", "T", "A", "G"]
        for idx, base in enumerate(bases):
            cx = x + 42 + idx * 30
            parts.append(f'<circle cx="{cx}" cy="{y + 44}" r="12" fill="#426b69" opacity="0.25" stroke="#426b69"/>')
            parts.append(f'<text x="{cx - 4}" y="{y + 48}" font-family="Arial" font-size="11" fill="#333">{base}</text>')
        parts.append(f'<line x1="{x + 126}" y1="{y + 26}" x2="{x + 156}" y2="{y + 62}" stroke="#c6912f" stroke-width="3"/>')
    elif key == "TE_grammar_node_candidate":
        for idx, color in enumerate(["#426b69", "#7b4b73", "#9b5f3f", "#3f6b8f"]):
            parts.append(f'<rect x="{x + 28 + idx * 48}" y="{y + 30}" width="42" height="28" fill="{color}" opacity="0.35" stroke="{color}"/>')
        parts.append(f'<path d="M {x + 50} {y + 24} C {x + 88} {y + 8}, {x + 134} {y + 8}, {x + 188} {y + 24}" fill="none" stroke="#65736f" stroke-dasharray="4 4"/>')
    elif key == "possibility_gate_candidate":
        parts.extend(
            [
                f'<circle cx="{x + 46}" cy="{y + 44}" r="12" fill="#426b69" opacity="0.35"/>',
                f'<line x1="{x + 58}" y1="{y + 44}" x2="{x + 124}" y2="{y + 24}" stroke="#65736f"/>',
                f'<line x1="{x + 58}" y1="{y + 44}" x2="{x + 124}" y2="{y + 44}" stroke="#65736f"/>',
                f'<line x1="{x + 58}" y1="{y + 44}" x2="{x + 124}" y2="{y + 64}" stroke="#65736f"/>',
                f'<rect x="{x + 126}" y="{y + 18}" width="18" height="52" fill="#fff5d9" stroke="#c6912f"/>',
                f'<circle cx="{x + 184}" cy="{y + 24}" r="9" fill="#7b4b73" opacity="0.25"/>',
                f'<circle cx="{x + 184}" cy="{y + 44}" r="9" fill="#7b4b73" opacity="0.25"/>',
                f'<circle cx="{x + 184}" cy="{y + 64}" r="9" fill="#7b4b73" opacity="0.25"/>',
            ]
        )
    elif key == "chromatin_motion_oscillator_candidate":
        parts.extend(
            [
                f'<path d="M {x + 28} {y + 46} C {x + 58} {y + 16}, {x + 92} {y + 76}, {x + 124} {y + 46} C {x + 154} {y + 16}, {x + 188} {y + 76}, {x + 218} {y + 46}" fill="none" stroke="#7b4b73" stroke-width="3"/>',
                f'<circle cx="{x + 126}" cy="{y + 46}" r="8" fill="#426b69"/>',
                f'<circle cx="{x + 126}" cy="{y + 46}" r="22" fill="none" stroke="#426b69" stroke-dasharray="3 4"/>',
            ]
        )
    else:
        parts.extend(
            [
                f'<line x1="{x + 28}" y1="{y + 66}" x2="{x + 218}" y2="{y + 66}" stroke="#65736f"/>',
                f'<path d="M {x + 28} {y + 64} L {x + 76} {y + 58} L {x + 102} {y + 22} L {x + 128} {y + 58} L {x + 218} {y + 62}" fill="none" stroke="#7b4b73" stroke-width="3"/>',
                f'<circle cx="{x + 102}" cy="{y + 22}" r="5" fill="#c6912f"/>',
            ]
        )
    return parts


def _model_text(card: dict) -> dict[str, str]:
    bridge = card.get("mechanistic_bridge", {}) or {}
    observed = card.get("observed_feature_evidence", {}) or {}
    supporting = observed.get("supporting_features") or card.get("top_supporting_features") or []
    feature = bridge.get("measured_feature") or ", ".join(str(item) for item in supporting) or "sequence-derived proxy"
    behavior = bridge.get("proposed_dynamic_phenotype")
    if not behavior:
        hypothesis = card.get("primitive_hypothesis", {}) or {}
        behavior = hypothesis.get("hypothesis_statement") or card.get("predicted_hidden_property") or "candidate primitive behavior"
    return {
        "feature": str(feature),
        "behavior": str(behavior),
        "assay": str(card.get("recommended_primitive_assay", "") or ""),
    }


def _summary_rows(card: dict) -> str:
    rows = [
        ("Coordinates", card.get("coordinates", "")),
        ("Primitive class", card.get("primitive_class", "")),
        ("Confidence", _fmt_float(card.get("primitive_confidence", np.nan), 3)),
        ("Observed evidence", ", ".join((card.get("observed_feature_evidence", {}) or {}).get("supporting_features", []) or [])),
        ("Primitive hypothesis", (card.get("primitive_hypothesis", {}) or {}).get("hypothesis_statement", "")),
        ("Mechanistic bridge", _bridge_text(card)),
        ("Assay", card.get("recommended_primitive_assay", "")),
        ("Key test", card.get("key_interaction_test", "")),
        ("Caveat", card.get("interpretation_caveat", "")),
    ]
    return "\n".join(f"<tr><th>{escape(str(key))}</th><td>{escape(str(value))}</td></tr>" for key, value in rows if str(value))


def _bridge_text(card: dict) -> str:
    bridge = card.get("mechanistic_bridge", {}) or {}
    measured = bridge.get("measured_feature", "")
    behavior = bridge.get("proposed_dynamic_phenotype", "")
    return f"{measured} -> {behavior}" if measured or behavior else ""


def _visual_alt(card: dict, focus: dict, intervals: pd.DataFrame) -> str:
    unique = 0 if intervals is None or intervals.empty else intervals["primitive_class"].nunique()
    return f"{focus['name']} at {focus['chrom']}:{focus['start']}-{focus['end']}; {unique} primitive classes in local view; candidate behavior sketch for {card.get('primitive_class', 'candidate')}"


def _append_axis_ticks(parts: list[str], start: int, end: int, left: float, y: float, width: float) -> None:
    if end <= start:
        return
    ticks = np.linspace(start, end, num=4)
    for tick in ticks:
        x = _scale(float(tick), start, end, left, width)
        parts.append(f'<line x1="{x:.1f}" y1="{y - 5}" x2="{x:.1f}" y2="{y}" stroke="#65736f"/>')
        parts.append(f'<text x="{x - 20:.1f}" y="{y + 14}" font-family="Arial" font-size="10" fill="#555">{escape(_fmt_bp(float(tick)))}</text>')


def _scaled_interval(start: int, end: int, view: dict, left: float, width: float) -> tuple[float, float]:
    x1 = _scale(max(int(start), int(view["start"])), int(view["start"]), int(view["end"]), left, width)
    x2 = _scale(min(int(end), int(view["end"])), int(view["start"]), int(view["end"]), left, width)
    return x1, max(2.0, x2 - x1)


def _scale(value: float, start: float, end: float, left: float, width: float) -> float:
    if end <= start:
        return left
    return left + ((float(value) - float(start)) / (float(end) - float(start))) * width


def _wrapped_svg_text(text: str, x: float, y: float, max_chars: int, *, max_lines: int, fill: str) -> list[str]:
    words = str(text).split()
    lines: list[str] = []
    current = ""
    for word in words:
        if len(current) + len(word) + 1 <= max_chars:
            current = f"{current} {word}".strip()
        else:
            if current:
                lines.append(current)
            current = word
        if len(lines) >= max_lines:
            break
    if current and len(lines) < max_lines:
        lines.append(current)
    if len(lines) == max_lines and len(" ".join(words)) > len(" ".join(lines)):
        lines[-1] = lines[-1].rstrip(".") + "..."
    return [
        f'<text x="{x}" y="{y + idx * 15}" font-family="Arial" font-size="11" fill="{fill}">{escape(line)}</text>'
        for idx, line in enumerate(lines)
    ]


def _primitive_label(value: str, max_chars: int = 32) -> str:
    text = str(value).replace("_candidate_score", "").replace("_candidate", "").replace("_", " ")
    return text if len(text) <= max_chars else text[: max_chars - 3] + "..."


def _color_for_primitive(primitive: str) -> str:
    idx = sum(ord(ch) for ch in str(primitive)) % len(PALETTE)
    return PALETTE[idx]


def _interval_tooltip(row: object) -> str:
    pieces = [
        str(getattr(row, "region_id", "")),
        str(getattr(row, "primitive_class", "")),
        f"confidence {_fmt_float(getattr(row, 'primitive_confidence', np.nan), 3)}",
    ]
    if hasattr(row, "residual_zscore"):
        pieces.append(f"residual z {_fmt_float(getattr(row, 'residual_zscore', np.nan), 2)}")
    if hasattr(row, "matched_null_zscore"):
        pieces.append(f"matched-null z {_fmt_float(getattr(row, 'matched_null_zscore', np.nan), 2)}")
    if hasattr(row, "window_size"):
        pieces.append(f"window {getattr(row, 'window_size')} bp")
    return "; ".join(piece for piece in pieces if piece)


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    return (slug or "card")[:140]


def _unique_slug(slug: str, used: set[str]) -> str:
    base = slug
    idx = 2
    while slug in used:
        slug = f"{base}_{idx}"
        idx += 1
    used.add(slug)
    return slug


def _fmt_bp(value: float) -> str:
    number = float(value)
    if abs(number) >= 1_000_000:
        return f"{number / 1_000_000:.2f} Mb"
    if abs(number) >= 1_000:
        return f"{number / 1_000:.1f} kb"
    return f"{int(round(number))} bp"


def _fmt_float(value: object, digits: int) -> str:
    number = _finite_float(value, default=np.nan)
    if not np.isfinite(number):
        return "NA"
    return f"{number:.{digits}f}"


def _finite_float(value: object, *, default: float) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if np.isfinite(number) else default
