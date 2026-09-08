"""Locus-level candidate aggregation for overlapping multiscale windows."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from darkdna.io.bed import write_bed


LOCUS_COLUMNS = [
    "locus_id",
    "chrom",
    "start",
    "end",
    "primitive_class",
    "representative_region_id",
    "n_windows",
    "window_sizes",
    "scale_count",
    "scale_discovery_window_size",
    "scale_validation_window_sizes",
    "scale_validation_status",
    "support_score",
    "max_primitive_confidence",
    "max_residual_zscore",
    "max_matched_null_zscore",
    "raw_min_empirical_p_value",
    "locus_effective_test_count",
    "locus_empirical_p_value",
    "primitive_bh_q_value",
    "global_bh_q_value",
    "locus_significant_fdr_0_05",
    "locus_significant_fdr_0_10",
    "block_id",
    "chromosome_cv_fold",
    "member_region_ids",
    "artifact_risk_flags",
    "locus_level_note",
]


BOOTSTRAP_COLUMNS = [
    "primitive_class",
    "n_loci",
    "n_blocks",
    "block_size_bp",
    "mean_locus_support",
    "bootstrap_ci_low",
    "bootstrap_ci_high",
    "block_bootstrap_status",
]


def benjamini_hochberg_qvalues(p_values: pd.Series) -> pd.Series:
    """Return Benjamini-Hochberg q-values aligned to the input index."""

    p = pd.to_numeric(p_values, errors="coerce")
    q = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p.notna()
    if int(valid.sum()) == 0:
        return q
    ordered = p.loc[valid].sort_values(kind="mergesort")
    n = len(ordered)
    ranks = np.arange(1, n + 1, dtype=float)
    raw_q = ordered.to_numpy(dtype=float) * n / ranks
    monotone_q = np.minimum.accumulate(raw_q[::-1])[::-1]
    q.loc[ordered.index] = np.clip(monotone_q, 0.0, 1.0)
    return q


def _numeric(df: pd.DataFrame, column: str, default: float = math.nan) -> pd.Series:
    if column not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[column], errors="coerce")


def _candidate_windows(windows: pd.DataFrame, labels: pd.DataFrame, residuals: pd.DataFrame) -> pd.DataFrame:
    if windows is None or labels is None or windows.empty or labels.empty:
        return pd.DataFrame()
    coord_cols = [col for col in ["region_id", "chrom", "start", "end", "window_size", "artifact_risk_flags"] if col in windows.columns]
    if not {"region_id", "chrom", "start", "end"}.issubset(coord_cols):
        return pd.DataFrame()
    candidates = windows[coord_cols].merge(labels, on="region_id", how="inner")
    candidates = candidates[candidates["primitive_class"].fillna("").astype(str).ne("no_call")].copy()
    candidates = candidates[candidates["primitive_class"].fillna("").astype(str).ne("")]
    if candidates.empty:
        return candidates

    if residuals is not None and not residuals.empty and {"region_id", "primitive_score_name"}.issubset(candidates.columns):
        extra_cols = [
            col
            for col in [
                "region_id",
                "primitive",
                "observed_score",
                "predicted_classical_score",
                "residual_score",
                "classical_explanation_fraction",
                "model_method",
            ]
            if col in residuals.columns
        ]
        if {"region_id", "primitive"}.issubset(extra_cols):
            selected = residuals[extra_cols].copy()
            candidates = candidates.merge(
                selected,
                left_on=["region_id", "primitive_score_name"],
                right_on=["region_id", "primitive"],
                how="left",
            )

    candidates["start"] = pd.to_numeric(candidates["start"], errors="coerce").astype("Int64")
    candidates["end"] = pd.to_numeric(candidates["end"], errors="coerce").astype("Int64")
    candidates = candidates.dropna(subset=["start", "end"]).copy()
    candidates["start"] = candidates["start"].astype(int)
    candidates["end"] = candidates["end"].astype(int)
    if "window_size" not in candidates.columns:
        candidates["window_size"] = candidates["end"] - candidates["start"]
    candidates["window_size"] = pd.to_numeric(candidates["window_size"], errors="coerce").fillna(candidates["end"] - candidates["start"]).astype(int)
    residual_z = _numeric(candidates, "residual_zscore")
    matched_z = _numeric(candidates, "matched_null_zscore")
    candidates["support_score"] = residual_z.clip(lower=0.0).fillna(0.0) + 0.5 * matched_z.clip(lower=0.0).fillna(0.0)
    candidates["primitive_confidence"] = _numeric(candidates, "primitive_confidence")
    candidates["empirical_p_value"] = pd.to_numeric(candidates.get("empirical_p_value", np.nan), errors="coerce")
    return candidates.sort_values(["primitive_class", "chrom", "start", "end"], kind="mergesort")


def _window_sizes_text(group: pd.DataFrame) -> str:
    sizes = sorted({int(v) for v in pd.to_numeric(group["window_size"], errors="coerce").dropna()})
    return ";".join(str(v) for v in sizes)


def _flags_text(group: pd.DataFrame) -> str:
    flags: list[str] = []
    if "artifact_risk_flags" not in group.columns:
        return ""
    for value in group["artifact_risk_flags"].fillna("").astype(str):
        for flag in value.split(";"):
            flag = flag.strip()
            if flag and flag not in flags:
                flags.append(flag)
    return ";".join(flags)


def _member_region_ids(group: pd.DataFrame, limit: int = 80) -> str:
    ids = [str(region_id) for region_id in group["region_id"].dropna()]
    if len(ids) > limit:
        return ";".join(ids[:limit] + [f"...{len(ids) - limit}_more"])
    return ";".join(ids)


def _scale_status(group: pd.DataFrame, representative: pd.Series) -> tuple[int, str, str, str]:
    sizes = sorted({int(v) for v in pd.to_numeric(group["window_size"], errors="coerce").dropna()})
    discovery = int(representative.get("window_size", sizes[0] if sizes else 0) or 0)
    validation_sizes = [size for size in sizes if size != discovery]
    if not sizes:
        return 0, "", "", "unknown_scale"
    if validation_sizes:
        return len(sizes), str(discovery), ";".join(str(size) for size in validation_sizes), "cross_scale_supported"
    return len(sizes), str(discovery), "", "single_scale_only"


def _block_id(chrom: str, start: int, end: int, block_size: int) -> str:
    if block_size <= 0:
        return str(chrom)
    first = int(start) // block_size
    last = max(int(start), int(end) - 1) // block_size
    if first == last:
        return f"{chrom}:block_{first}"
    return ";".join(f"{chrom}:block_{idx}" for idx in range(first, last + 1))


def _locus_p_value(group: pd.DataFrame, scale_count: int) -> tuple[float, int, float]:
    p_values = pd.to_numeric(group.get("empirical_p_value", np.nan), errors="coerce").dropna()
    if p_values.empty:
        return math.nan, max(1, scale_count), math.nan
    raw_min = float(p_values.min())
    effective_tests = max(1, int(scale_count))
    return raw_min, effective_tests, float(min(1.0, raw_min * effective_tests))


def _build_locus_row(group: pd.DataFrame, block_size: int) -> dict:
    support = _numeric(group, "support_score")
    best_idx = support.idxmax()
    representative = group.loc[best_idx]
    start = int(group["start"].min())
    end = int(group["end"].max())
    chrom = str(representative["chrom"])
    primitive = str(representative["primitive_class"])
    scale_count, discovery_size, validation_sizes, scale_status = _scale_status(group, representative)
    raw_min_p, effective_tests, locus_p = _locus_p_value(group, scale_count)
    locus_id = f"{chrom}:{start}-{end}:{primitive}"
    return {
        "locus_id": locus_id,
        "chrom": chrom,
        "start": start,
        "end": end,
        "primitive_class": primitive,
        "representative_region_id": str(representative["region_id"]),
        "n_windows": int(len(group)),
        "window_sizes": _window_sizes_text(group),
        "scale_count": int(scale_count),
        "scale_discovery_window_size": discovery_size,
        "scale_validation_window_sizes": validation_sizes,
        "scale_validation_status": scale_status,
        "support_score": float(support.max()) if support.notna().any() else math.nan,
        "max_primitive_confidence": float(_numeric(group, "primitive_confidence").max()),
        "max_residual_zscore": float(_numeric(group, "residual_zscore").max()),
        "max_matched_null_zscore": float(_numeric(group, "matched_null_zscore").max()),
        "raw_min_empirical_p_value": raw_min_p,
        "locus_effective_test_count": int(effective_tests),
        "locus_empirical_p_value": locus_p,
        "block_id": _block_id(chrom, start, end, block_size),
        "chromosome_cv_fold": chrom,
        "member_region_ids": _member_region_ids(group),
        "artifact_risk_flags": _flags_text(group),
        "locus_level_note": "Merged overlapping candidate windows; count and FDR are locus-level, not window-level.",
    }


def merge_candidate_loci(
    windows: pd.DataFrame,
    labels: pd.DataFrame,
    residuals: pd.DataFrame,
    *,
    merge_gap: int = 0,
    block_size: int = 100_000,
) -> pd.DataFrame:
    """Merge overlapping candidate windows into countable locus-level intervals.

    Windows are merged within chromosome and primitive class. This deliberately
    avoids treating half-step or nested multiscale windows as independent
    candidate observations.
    """

    candidates = _candidate_windows(windows, labels, residuals)
    if candidates.empty:
        return pd.DataFrame(columns=LOCUS_COLUMNS)

    rows: list[dict] = []
    for (_, _), group in candidates.groupby(["primitive_class", "chrom"], sort=True):
        current: list[pd.Series] = []
        current_end: int | None = None
        for _, row in group.sort_values(["start", "end"], kind="mergesort").iterrows():
            starts_new = current_end is None or (int(row["start"]) > current_end + merge_gap if merge_gap > 0 else int(row["start"]) >= current_end)
            if starts_new and current:
                rows.append(_build_locus_row(pd.DataFrame(current), block_size))
                current = []
            current.append(row)
            current_end = max(int(current_end or row["end"]), int(row["end"]))
        if current:
            rows.append(_build_locus_row(pd.DataFrame(current), block_size))

    loci = pd.DataFrame(rows)
    if loci.empty:
        return pd.DataFrame(columns=LOCUS_COLUMNS)
    loci["global_bh_q_value"] = benjamini_hochberg_qvalues(loci["locus_empirical_p_value"])
    loci["primitive_bh_q_value"] = np.nan
    for _, idx in loci.groupby("primitive_class").groups.items():
        loci.loc[idx, "primitive_bh_q_value"] = benjamini_hochberg_qvalues(loci.loc[idx, "locus_empirical_p_value"])
    loci["locus_significant_fdr_0_05"] = pd.to_numeric(loci["global_bh_q_value"], errors="coerce") <= 0.05
    loci["locus_significant_fdr_0_10"] = pd.to_numeric(loci["global_bh_q_value"], errors="coerce") <= 0.10
    sort_q = pd.to_numeric(loci["global_bh_q_value"], errors="coerce").fillna(1.0)
    loci = loci.assign(_sort_q=sort_q).sort_values(
        ["_sort_q", "locus_empirical_p_value", "support_score", "chrom", "start"],
        ascending=[True, True, False, True, True],
        kind="mergesort",
    )
    loci = loci.drop(columns=["_sort_q"]).reset_index(drop=True)
    return loci[[col for col in LOCUS_COLUMNS if col in loci.columns]]


def block_bootstrap_locus_summary(
    loci: pd.DataFrame,
    *,
    block_size: int = 100_000,
    n_bootstrap: int = 500,
    seed: int = 13,
) -> pd.DataFrame:
    """Summarize support stability by resampling genomic blocks, not windows."""

    if loci is None or loci.empty:
        return pd.DataFrame(columns=BOOTSTRAP_COLUMNS)
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for primitive, group in loci.groupby("primitive_class", sort=True):
        support = pd.to_numeric(group["support_score"], errors="coerce").fillna(0.0)
        block_values = group.assign(_support=support).groupby("block_id")["_support"].mean()
        observed = float(support.mean()) if len(support) else 0.0
        if len(block_values) < 2:
            low = high = observed
            status = "insufficient_independent_blocks"
        else:
            samples = []
            values = block_values.to_numpy(dtype=float)
            for _ in range(max(1, int(n_bootstrap))):
                samples.append(float(rng.choice(values, size=len(values), replace=True).mean()))
            low, high = np.percentile(samples, [2.5, 97.5])
            status = "ok"
        rows.append(
            {
                "primitive_class": primitive,
                "n_loci": int(len(group)),
                "n_blocks": int(len(block_values)),
                "block_size_bp": int(block_size),
                "mean_locus_support": observed,
                "bootstrap_ci_low": float(low),
                "bootstrap_ci_high": float(high),
                "block_bootstrap_status": status,
            }
        )
    return pd.DataFrame(rows, columns=BOOTSTRAP_COLUMNS)


def write_candidate_locus_outputs(
    windows: pd.DataFrame,
    labels: pd.DataFrame,
    residuals: pd.DataFrame,
    outdir: str | Path,
    *,
    merge_gap: int = 0,
    block_size: int = 100_000,
) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    loci = merge_candidate_loci(windows, labels, residuals, merge_gap=merge_gap, block_size=block_size)
    summary = block_bootstrap_locus_summary(loci, block_size=block_size)

    paths = {
        "candidate_loci_parquet": out / "candidate_loci.parquet",
        "candidate_loci_tsv": out / "candidate_loci.tsv",
        "candidate_loci_bed": out / "candidate_loci.bed",
        "candidate_loci_block_bootstrap_parquet": out / "candidate_loci_block_bootstrap.parquet",
        "candidate_loci_block_bootstrap_tsv": out / "candidate_loci_block_bootstrap.tsv",
    }
    loci.to_parquet(paths["candidate_loci_parquet"], index=False)
    loci.to_csv(paths["candidate_loci_tsv"], sep="\t", index=False)
    summary.to_parquet(paths["candidate_loci_block_bootstrap_parquet"], index=False)
    summary.to_csv(paths["candidate_loci_block_bootstrap_tsv"], sep="\t", index=False)

    bed = loci.copy()
    if bed.empty:
        paths["candidate_loci_bed"].write_text("", encoding="utf-8")
    else:
        bed["name"] = bed["locus_id"]
        bed["score"] = (pd.to_numeric(bed["max_primitive_confidence"], errors="coerce").fillna(0.0).clip(0.0, 1.0) * 1000).round().astype(int)
        bed["strand"] = "."
        write_bed(bed, paths["candidate_loci_bed"], columns=["chrom", "start", "end", "name", "score", "strand"])
    return paths
