"""GFF3/GTF parsing with non-model genome defaults."""

from __future__ import annotations

from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd

from darkdna.utils.progress import ProgressReporter, progress_message


def parse_attributes(attr_text: str, fmt: str = "gff3") -> dict[str, str]:
    attrs: dict[str, str] = {}
    if not attr_text or attr_text == ".":
        return attrs
    if fmt == "gtf" or ' "' in attr_text:
        for token in attr_text.strip().rstrip(";").split(";"):
            token = token.strip()
            if not token:
                continue
            if " " in token:
                key, value = token.split(" ", 1)
                attrs[key] = value.strip().strip('"')
    else:
        for token in attr_text.split(";"):
            if not token:
                continue
            if "=" in token:
                key, value = token.split("=", 1)
                attrs[key] = unquote(value)
            else:
                attrs[token] = ""
    return attrs


def read_gff(path: str | Path | None, *, progress: bool = False, label: str = "gff") -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    rows = []
    reporter = ProgressReporter(label, total=None, min_interval=10.0) if progress else None
    if reporter:
        progress_message(label, f"reading {Path(path)}")
        reporter.start("parsing GFF/GTF records")
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                if reporter:
                    reporter.update(line_no, message=f"rows={len(rows)}")
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                if reporter:
                    reporter.update(line_no, message=f"rows={len(rows)}")
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = parts[:9]
            fmt = "gtf" if ' "' in attrs else "gff3"
            parsed = parse_attributes(attrs, fmt=fmt)
            rows.append(
                {
                    "chrom": chrom,
                    "source": source,
                    "feature": feature,
                    "start": int(start) - 1,
                    "end": int(end),
                    "score": score,
                    "strand": strand if strand in {"+", "-"} else ".",
                    "phase": phase,
                    "attributes": parsed,
                    "id": parsed.get("ID") or parsed.get("gene_id") or parsed.get("Name"),
                    "parent": parsed.get("Parent") or parsed.get("transcript_id"),
                    "name": parsed.get("Name") or parsed.get("gene_name") or parsed.get("gene_id") or parsed.get("ID"),
                    "biotype": parsed.get("biotype")
                    or parsed.get("gene_biotype")
                    or parsed.get("transcript_biotype")
                    or parsed.get("gbkey"),
                    "family": parsed.get("family") or parsed.get("Family") or parsed.get("Name"),
                    "class": parsed.get("class") or parsed.get("Class") or parsed.get("type"),
                    "superfamily": parsed.get("superfamily") or parsed.get("Superfamily"),
                    "divergence": parsed.get("divergence") or parsed.get("perc_div"),
                }
            )
            if reporter:
                reporter.update(line_no, message=f"rows={len(rows)}")
    if reporter:
        reporter.finish(f"rows={len(rows)}")
    return pd.DataFrame(rows)


def feature_table(path: str | Path | None, feature_names: set[str]) -> pd.DataFrame:
    df = read_gff(path)
    if df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "name", "strand", "feature"])
    mask = df["feature"].str.lower().isin({name.lower() for name in feature_names})
    out = df.loc[mask, ["chrom", "start", "end", "name", "strand", "feature", "biotype"]].copy()
    return out


def annotation_tables(path: str | Path | None, promoter_bp: int = 1000, *, progress: bool = False) -> dict[str, pd.DataFrame]:
    df = read_gff(path, progress=progress, label="annotation-gff")
    if df.empty:
        empty = pd.DataFrame(columns=["chrom", "start", "end", "name", "strand", "feature"])
        return {key: empty.copy() for key in ["genes", "exons", "utrs", "promoters", "introns", "tss"]}

    lower = df["feature"].str.lower()
    genes = df.loc[lower.isin({"gene", "pseudogene", "lnc_rna_gene"})].copy()
    exons = df.loc[lower.isin({"exon", "cds"})].copy()
    utrs = df.loc[lower.str.contains("utr", na=False)].copy()
    if progress:
        progress_message("annotation-gff", f"genes={len(genes):,} exons/cds={len(exons):,} utrs={len(utrs):,}")

    promoter_rows = []
    tss_rows = []
    promoter_reporter = ProgressReporter("annotation-promoters", total=len(genes), min_interval=10.0) if progress else None
    if promoter_reporter:
        promoter_reporter.start("building promoters/TSS")
    for idx, row in enumerate(genes.itertuples(), start=1):
        if row.strand == "-":
            tss = int(row.end)
            start, end = tss, tss + promoter_bp
        else:
            tss = int(row.start)
            start, end = max(0, tss - promoter_bp), tss
        name = row.name or row.id or f"{row.chrom}:{row.start}-{row.end}"
        promoter_rows.append({"chrom": row.chrom, "start": start, "end": end, "name": name, "strand": row.strand, "feature": "promoter"})
        tss_rows.append({"chrom": row.chrom, "start": tss, "end": tss + 1, "name": name, "strand": row.strand, "feature": "tss"})
        if promoter_reporter:
            promoter_reporter.update(idx, message=str(row.chrom))
    if promoter_reporter:
        promoter_reporter.finish()

    promoters = pd.DataFrame(promoter_rows)
    tss = pd.DataFrame(tss_rows)

    intron_rows = []
    exon_index = {}
    if not exons.empty:
        for chrom, group in exons.sort_values(["chrom", "start", "end"]).groupby("chrom", sort=False):
            exon_index[str(chrom)] = {
                "starts": group["start"].to_numpy(dtype=int),
                "ends": group["end"].to_numpy(dtype=int),
            }
    intron_reporter = ProgressReporter("annotation-introns", total=len(genes), min_interval=10.0) if progress else None
    if intron_reporter:
        intron_reporter.start("building intron intervals")
    for idx, gene in enumerate(genes.itertuples(), start=1):
        chrom_exons = exon_index.get(str(gene.chrom))
        if chrom_exons is None:
            gene_exon_intervals: list[tuple[int, int]] = []
        else:
            starts = chrom_exons["starts"]
            ends = chrom_exons["ends"]
            cutoff = int(np.searchsorted(starts, int(gene.end), side="left"))
            mask = ends[:cutoff] > int(gene.start)
            gene_exon_intervals = [
                (max(int(gene.start), int(exon_start)), min(int(gene.end), int(exon_end)))
                for exon_start, exon_end in zip(starts[:cutoff][mask], ends[:cutoff][mask])
                if int(exon_end) > int(gene.start) and int(exon_start) < int(gene.end)
            ]
        cursor = int(gene.start)
        for exon_start, exon_end in sorted(gene_exon_intervals):
            if exon_end <= cursor:
                continue
            if exon_start > cursor:
                intron_rows.append({"chrom": gene.chrom, "start": cursor, "end": exon_start, "name": gene.name, "strand": gene.strand, "feature": "intron"})
            cursor = max(cursor, exon_end)
        if cursor < int(gene.end):
            intron_rows.append({"chrom": gene.chrom, "start": cursor, "end": int(gene.end), "name": gene.name, "strand": gene.strand, "feature": "intron"})
        if intron_reporter:
            intron_reporter.update(idx, message=str(gene.chrom))
    if intron_reporter:
        intron_reporter.finish(f"introns={len(intron_rows)}")
    introns = pd.DataFrame(intron_rows)

    cols = ["chrom", "start", "end", "name", "strand", "feature"]
    return {
        "genes": genes.rename(columns={"feature": "feature"})[[c for c in cols if c in genes.columns]],
        "exons": exons[[c for c in cols if c in exons.columns]],
        "utrs": utrs[[c for c in cols if c in utrs.columns]],
        "promoters": promoters if not promoters.empty else pd.DataFrame(columns=cols),
        "introns": introns if not introns.empty else pd.DataFrame(columns=cols),
        "tss": tss if not tss.empty else pd.DataFrame(columns=cols),
    }


def read_te_annotation(path: str | Path | None, *, progress: bool = False) -> pd.DataFrame:
    if not path:
        return pd.DataFrame(columns=["chrom", "start", "end", "family", "class", "superfamily", "strand", "divergence"])
    p = Path(path)
    if p.suffix.lower() in {".gff", ".gff3", ".gtf"}:
        df = read_gff(path, progress=progress, label="te-gff")
        if df.empty:
            return pd.DataFrame(columns=["chrom", "start", "end", "family", "class", "superfamily", "strand", "divergence"])
        return df[["chrom", "start", "end", "family", "class", "superfamily", "strand", "divergence"]].copy()
    from .bed import read_bed

    bed = read_bed(path)
    if bed.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "family", "class", "superfamily", "strand", "divergence"])
    out = bed.copy()
    out["family"] = out.get("name", pd.Series([None] * len(out)))
    out["class"] = out.get("field_6", pd.Series([None] * len(out)))
    out["superfamily"] = out.get("field_7", pd.Series([None] * len(out)))
    out["divergence"] = out.get("field_8", pd.Series([None] * len(out)))
    return out[["chrom", "start", "end", "family", "class", "superfamily", "strand", "divergence"]]
