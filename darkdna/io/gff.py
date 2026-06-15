"""GFF3/GTF parsing with non-model genome defaults."""

from __future__ import annotations

from pathlib import Path
from urllib.parse import unquote

import pandas as pd


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


def read_gff(path: str | Path | None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    rows = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
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
    return pd.DataFrame(rows)


def feature_table(path: str | Path | None, feature_names: set[str]) -> pd.DataFrame:
    df = read_gff(path)
    if df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "name", "strand", "feature"])
    mask = df["feature"].str.lower().isin({name.lower() for name in feature_names})
    out = df.loc[mask, ["chrom", "start", "end", "name", "strand", "feature", "biotype"]].copy()
    return out


def annotation_tables(path: str | Path | None, promoter_bp: int = 1000) -> dict[str, pd.DataFrame]:
    df = read_gff(path)
    if df.empty:
        empty = pd.DataFrame(columns=["chrom", "start", "end", "name", "strand", "feature"])
        return {key: empty.copy() for key in ["genes", "exons", "utrs", "promoters", "introns", "tss"]}

    lower = df["feature"].str.lower()
    genes = df.loc[lower.isin({"gene", "pseudogene", "lnc_rna_gene"})].copy()
    exons = df.loc[lower.isin({"exon", "cds"})].copy()
    utrs = df.loc[lower.str.contains("utr", na=False)].copy()

    promoter_rows = []
    tss_rows = []
    for row in genes.itertuples():
        if row.strand == "-":
            tss = int(row.end)
            start, end = tss, tss + promoter_bp
        else:
            tss = int(row.start)
            start, end = max(0, tss - promoter_bp), tss
        name = row.name or row.id or f"{row.chrom}:{row.start}-{row.end}"
        promoter_rows.append({"chrom": row.chrom, "start": start, "end": end, "name": name, "strand": row.strand, "feature": "promoter"})
        tss_rows.append({"chrom": row.chrom, "start": tss, "end": tss + 1, "name": name, "strand": row.strand, "feature": "tss"})

    promoters = pd.DataFrame(promoter_rows)
    tss = pd.DataFrame(tss_rows)

    intron_rows = []
    for gene in genes.itertuples():
        gene_exons = exons[(exons["chrom"] == gene.chrom) & (exons["start"] >= gene.start) & (exons["end"] <= gene.end)]
        cursor = int(gene.start)
        for exon in gene_exons.sort_values("start").itertuples():
            if int(exon.start) > cursor:
                intron_rows.append({"chrom": gene.chrom, "start": cursor, "end": int(exon.start), "name": gene.name, "strand": gene.strand, "feature": "intron"})
            cursor = max(cursor, int(exon.end))
        if cursor < int(gene.end):
            intron_rows.append({"chrom": gene.chrom, "start": cursor, "end": int(gene.end), "name": gene.name, "strand": gene.strand, "feature": "intron"})
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


def read_te_annotation(path: str | Path | None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame(columns=["chrom", "start", "end", "family", "class", "superfamily", "strand", "divergence"])
    p = Path(path)
    if p.suffix.lower() in {".gff", ".gff3", ".gtf"}:
        df = read_gff(path)
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
