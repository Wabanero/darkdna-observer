"""Subset FASTA and GFF3 references for small integration fixtures."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from darkdna.io.fasta import read_fasta, write_chrom_sizes, write_fasta
from darkdna.io.gff import read_gff
from darkdna.toy_data import write_fai


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, span = region.split(":", 1)
    start_text, end_text = span.replace(",", "").split("-", 1)
    return chrom, int(start_text) - 1, int(end_text)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def maybe_decompress(path: Path, out: Path) -> Path:
    if path.suffix != ".gz":
        return path
    with gzip.open(path, "rt", encoding="utf-8") as src, out.open("w", encoding="utf-8") as dst:
        dst.write(src.read())
    return out


def subset_fasta(fasta: Path, out_fasta: Path, chrom: str, start: int, end: int, out_name: str | None = None) -> dict[str, str]:
    records = read_fasta(fasta)
    if chrom not in records:
        raise KeyError(f"Chromosome {chrom!r} not found in {fasta}")
    name = out_name or chrom
    subset = {name: records[chrom][start:end].upper()}
    write_fasta(subset, out_fasta)
    write_fai(subset, out_fasta)
    write_chrom_sizes({name: len(subset[name])}, out_fasta.parent / "chrom.sizes")
    return subset


def subset_gff3(gff3: Path, out_gff3: Path, chrom: str, start: int, end: int, out_chrom: str | None = None) -> None:
    df = read_gff(gff3)
    out_chrom = out_chrom or chrom
    if df.empty:
        out_gff3.write_text("##gff-version 3\n", encoding="utf-8")
        return
    subset = df[(df["chrom"].astype(str) == chrom) & (df["end"] > start) & (df["start"] < end)].copy()
    lines = ["##gff-version 3"]
    for row in subset.itertuples():
        new_start = max(0, int(row.start) - start) + 1
        new_end = min(end - start, int(row.end) - start)
        attrs = getattr(row, "attributes", {}) or {}
        attr_text = ";".join(f"{k}={v}" for k, v in attrs.items()) if attrs else f"ID={row.feature}_{new_start}_{new_end}"
        lines.append(
            "\t".join(
                [
                    out_chrom,
                    str(row.source),
                    str(row.feature),
                    str(new_start),
                    str(new_end),
                    str(row.score),
                    str(row.strand),
                    str(row.phase),
                    attr_text,
                ]
            )
        )
    out_gff3.write_text("\n".join(lines) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_checksums(outdir: Path) -> None:
    rows = []
    for path in sorted(outdir.glob("*")):
        if path.is_file() and path.name != "checksums.sha256":
            rows.append(f"{sha256(path)}  {path.name}")
    (outdir / "checksums.sha256").write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Subset FASTA/GFF3 into a small reference fixture.")
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--gff3", required=True, type=Path)
    parser.add_argument("--region", required=True, help="1-based region, for example Chr1:1-2000000")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--out-chrom", default=None)
    args = parser.parse_args()
    chrom, start, end = parse_region(args.region)
    args.out.mkdir(parents=True, exist_ok=True)
    subset_fasta(args.fasta, args.out / "genome.fa", chrom, start, end, args.out_chrom or chrom)
    subset_gff3(args.gff3, args.out / "genes.gff3", chrom, start, end, args.out_chrom or chrom)
    write_checksums(args.out)


if __name__ == "__main__":
    main()
