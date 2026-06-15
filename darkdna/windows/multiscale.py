"""Nested window relationships."""

from __future__ import annotations

import pandas as pd


def assign_multiscale_context(windows: pd.DataFrame) -> pd.DataFrame:
    """Assign parent and child region IDs for nested windows at the same locus."""

    if windows.empty:
        return windows.copy()
    out = windows.copy()
    out["parent_region_id"] = None
    out["child_region_ids"] = ""
    sizes = sorted(out["window_size"].dropna().unique())
    by_size = {size: out[out["window_size"] == size] for size in sizes}
    parent_for: dict[str, str | None] = {}
    children_for: dict[str, list[str]] = {str(r.region_id): [] for r in out.itertuples()}

    for size in sizes:
        bigger = [s for s in sizes if s > size]
        current = by_size[size]
        for row in current.itertuples():
            parent = None
            for candidate_size in bigger:
                candidates = by_size[candidate_size]
                same_chrom = candidates[candidates["chrom"] == row.chrom]
                hits = same_chrom[(same_chrom["start"] <= row.start) & (same_chrom["end"] >= row.end)]
                if not hits.empty:
                    parent = str(hits.iloc[0]["region_id"])
                    children_for.setdefault(parent, []).append(str(row.region_id))
                    break
            parent_for[str(row.region_id)] = parent

    out["parent_region_id"] = out["region_id"].map(parent_for)
    out["child_region_ids"] = out["region_id"].map(lambda rid: ",".join(children_for.get(str(rid), [])))
    return out


def scale_level_for_size(window_size: int, sizes: list[int]) -> int:
    return sorted(sizes).index(window_size)
