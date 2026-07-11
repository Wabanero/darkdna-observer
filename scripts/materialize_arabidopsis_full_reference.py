"""Materialize the downloaded Arabidopsis TAIR10 full reference.

This uses the already downloaded Ensembl Genomes source files from the small
fixture's ignored `_download` directory and writes a local full-genome reference
for full DarkDNA runs.
"""

from __future__ import annotations

import gzip
import hashlib
import shutil
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "data" / "reference" / "arabidopsis_TAIR10_chr1_2Mb" / "_download"
OUT_DIR = ROOT / "data" / "reference" / "arabidopsis_TAIR10_full"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_full_fasta(source: Path, destination: Path) -> list[tuple[str, int]]:
    sizes: list[tuple[str, int]] = []
    current: str | None = None
    length = 0
    with gzip.open(source, "rt", encoding="utf-8") as src, destination.open("w", encoding="utf-8", newline="\n") as dst:
        for line in src:
            if line.startswith(">"):
                if current is not None:
                    sizes.append((current, length))
                current = line[1:].split()[0]
                length = 0
                dst.write(line)
            else:
                seq = line.strip().upper()
                length += len(seq)
                dst.write(seq + "\n")
        if current is not None:
            sizes.append((current, length))
    return sizes


def write_checksums(outdir: Path) -> None:
    rows = []
    for path in sorted(outdir.glob("*")):
        if path.is_file() and path.name != "checksums.sha256":
            rows.append(f"{sha256(path)}  {path.name}")
    (outdir / "checksums.sha256").write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    fasta_source = SOURCE_DIR / "source.fa.gz"
    gff3_source = SOURCE_DIR / "source.gff3.gz"
    if not fasta_source.exists() or not gff3_source.exists():
        raise SystemExit(f"Missing downloaded source files in {SOURCE_DIR}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sizes = write_full_fasta(fasta_source, OUT_DIR / "genome.fa")

    with gzip.open(gff3_source, "rt", encoding="utf-8") as src, (OUT_DIR / "genes.gff3").open("w", encoding="utf-8", newline="\n") as dst:
        shutil.copyfileobj(src, dst)

    (OUT_DIR / "chrom.sizes").write_text(
        "".join(f"{chrom}\t{size}\n" for chrom, size in sizes),
        encoding="utf-8",
    )
    (OUT_DIR / "genome.fa.fai").write_text(
        "".join(f"{chrom}\t{size}\t0\t80\t81\n" for chrom, size in sizes),
        encoding="utf-8",
    )
    (OUT_DIR / "blacklist.bed").write_text("", encoding="utf-8")
    (OUT_DIR / "te_annotation.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    (OUT_DIR / "mappability.bedGraph").write_text(
        "".join(f"{chrom}\t0\t{size}\t1.0\n" for chrom, size in sizes),
        encoding="utf-8",
    )
    (OUT_DIR / "README.md").write_text(
        "# arabidopsis_TAIR10_full\n\n"
        "Arabidopsis thaliana TAIR10 full toplevel genome reference for local full-genome DarkDNA analysis.\n\n"
        f"- Prepared: {date.today().isoformat()}\n"
        "- Source FASTA: Ensembl Genomes Plants release 60 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz\n"
        "- Source GFF3: Ensembl Genomes Plants release 60 Arabidopsis_thaliana.TAIR10.60.gff3.gz\n"
        "- TE annotation: unavailable in downloaded source; empty placeholder written.\n"
        "- Mappability: placeholder 1.0 across all sequences.\n\n"
        "Full genomes should not be committed to this repository.\n",
        encoding="utf-8",
    )
    write_checksums(OUT_DIR)

    total_bp = sum(size for _, size in sizes)
    print(f"Wrote {len(sizes)} sequences and {total_bp} bp to {OUT_DIR}")
    for chrom, size in sizes:
        print(f"{chrom}\t{size}")


if __name__ == "__main__":
    main()
