"""Nested window relationships."""

from __future__ import annotations

import numpy as np
import pandas as pd


def assign_multiscale_context(windows: pd.DataFrame) -> pd.DataFrame:
    """Assign parent and child region IDs for nested windows at the same locus."""

    if windows.empty:
        return windows.copy()
    out = windows.copy()
    out["parent_region_id"] = None
    out["child_region_ids"] = ""
    sizes = sorted(out["window_size"].dropna().unique())
    by_size_chrom: dict[object, dict[str, pd.DataFrame]] = {}
    for size in sizes:
        by_size_chrom[size] = {
            str(chrom): group.sort_values("start").reset_index(drop=True)
            for chrom, group in out[out["window_size"] == size].groupby("chrom", sort=False)
        }
    parent_for: dict[str, str | None] = {}
    children_for: dict[str, list[str]] = {str(r.region_id): [] for r in out.itertuples()}

    for size in sizes:
        bigger = [s for s in sizes if s > size]
        current_by_chrom = by_size_chrom[size]
        bigger_indexes = {}
        for candidate_size in bigger:
            bigger_indexes[candidate_size] = {}
            for chrom, candidates in by_size_chrom[candidate_size].items():
                bigger_indexes[candidate_size][chrom] = {
                    "df": candidates,
                    "starts": candidates["start"].to_numpy(dtype=int),
                    "ends": candidates["end"].to_numpy(dtype=int),
                }
        for current in current_by_chrom.values():
            for row in current.itertuples():
                parent = None
                for candidate_size in bigger:
                    indexed = bigger_indexes[candidate_size].get(str(row.chrom))
                    if indexed is None:
                        continue
                    starts = indexed["starts"]
                    ends = indexed["ends"]
                    pos = int(np.searchsorted(starts, int(row.start), side="right")) - 1
                    if pos >= 0 and int(ends[pos]) >= int(row.end):
                        candidates = indexed["df"]
                        parent = str(candidates.iloc[pos]["region_id"])
                        children_for.setdefault(parent, []).append(str(row.region_id))
                        break
                parent_for[str(row.region_id)] = parent

    out["parent_region_id"] = out["region_id"].map(parent_for)
    out["child_region_ids"] = out["region_id"].map(lambda rid: ",".join(children_for.get(str(rid), [])))
    return out


def scale_level_for_size(window_size: int, sizes: list[int]) -> int:
    return sorted(sizes).index(window_size)
