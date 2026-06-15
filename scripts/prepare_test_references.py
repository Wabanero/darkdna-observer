"""Prepare optional small real reference fixtures for integration tests.

Normal pytest does not download data. Use --download explicitly, or provide
local source files through --source-dir with --no-download.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import urllib.request
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from darkdna.io.fasta import read_fasta, write_chrom_sizes, write_fasta
from darkdna.toy_data import write_fai
from scripts.subset_reference import parse_region, subset_gff3, subset_fasta, write_checksums


DATASETS = {
    "yeast_R64_chrI": {
        "region": "chrI:1-230218",
        "out_chrom": "chrI",
        "readme": "Saccharomyces cerevisiae S288C/R64 chromosome I technical integration fixture.",
        "fasta_candidates": ["genome.fa", "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"],
        "gff3_candidates": ["genes.gff3", "Saccharomyces_cerevisiae.R64-1-1.60.gff3"],
        "download_note": "Provide R64 FASTA and GFF3 from SGD/Ensembl/Fungi if automatic URLs fail.",
        "fasta_url": "",
        "gff3_url": "",
    },
    "arabidopsis_TAIR10_chr1_2Mb": {
        "region": "Chr1:1-2000000",
        "out_chrom": "Chr1",
        "readme": "Arabidopsis thaliana TAIR10 chromosome 1 first 2 Mb plant/non-model fixture.",
        "fasta_candidates": ["genome.fa", "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"],
        "gff3_candidates": ["genes.gff3", "Arabidopsis_thaliana.TAIR10.60.gff3"],
        "te_candidates": ["te_annotation.gff3", "TAIR10_TE.gff3"],
        "download_note": "Provide TAIR10 FASTA/GFF3/optional TE annotation from EnsemblPlants/Araport/TAIR if automatic URLs fail.",
        "fasta_url": "",
        "gff3_url": "",
    },
}


def find_first(source_dir: Path, names: list[str]) -> Path | None:
    for name in names:
        direct = source_dir / name
        gz = source_dir / f"{name}.gz"
        if direct.exists():
            return direct
        if gz.exists():
            return gz
    return None


def download_file(url: str, destination: Path) -> Path:
    if not url:
        raise RuntimeError("No download URL is configured for this dataset. Use --no-download --source-dir with local FASTA/GFF3 files.")
    print(f"Downloading {url} -> {destination}")
    urllib.request.urlretrieve(url, destination)
    return destination


def materialize_text_file(path: Path, out_path: Path) -> Path:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as src, out_path.open("w", encoding="utf-8") as dst:
            shutil.copyfileobj(src, dst)
        return out_path
    return path


def write_placeholder_tracks(outdir: Path, chrom: str, length: int) -> None:
    blacklist = outdir / "blacklist.bed"
    if not blacklist.exists():
        blacklist.write_text(f"{chrom}\t0\t1000\tsynthetic_test_blacklist\n", encoding="utf-8")
    mappability = outdir / "mappability.bedGraph"
    if not mappability.exists():
        mappability.write_text(f"{chrom}\t0\t{length}\t1.0\n", encoding="utf-8")


def prepare_dataset(args: argparse.Namespace) -> None:
    spec = DATASETS[args.dataset]
    out = args.out
    if out.exists() and args.force:
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    source = args.source_dir

    if args.download:
        source = out / "_download"
        source.mkdir(exist_ok=True)
        fasta_source = download_file(spec["fasta_url"], source / "source.fa.gz")
        gff3_source = download_file(spec["gff3_url"], source / "source.gff3.gz")
    else:
        if source is None:
            raise SystemExit("--no-download requires --source-dir containing local FASTA/GFF3 files.")
        fasta_source = find_first(source, spec["fasta_candidates"])
        gff3_source = find_first(source, spec["gff3_candidates"])
        if fasta_source is None or gff3_source is None:
            raise SystemExit(f"Missing local FASTA/GFF3 in {source}. {spec['download_note']}")

    fasta_plain = materialize_text_file(fasta_source, out / "_source.fa")
    gff3_plain = materialize_text_file(gff3_source, out / "_source.gff3")
    region = args.region or spec["region"]
    chrom, start, end = parse_region(region)
    out_chrom = spec["out_chrom"]
    subset = subset_fasta(fasta_plain, out / "genome.fa", chrom, start, end, out_chrom)
    subset_gff3(gff3_plain, out / "genes.gff3", chrom, start, end, out_chrom)

    if args.dataset == "arabidopsis_TAIR10_chr1_2Mb":
        te_source = find_first(source, spec.get("te_candidates", [])) if source else None
        if te_source:
            subset_gff3(materialize_text_file(te_source, out / "_source_te.gff3"), out / "te_annotation.gff3", chrom, start, end, out_chrom)
            te_note = "TE annotation was subset from local source."
        else:
            (out / "te_annotation.gff3").write_text("##gff-version 3\n", encoding="utf-8")
            te_note = "TE annotation unavailable; empty placeholder written for schema testing."
        write_placeholder_tracks(out, out_chrom, len(next(iter(subset.values()))))
    else:
        te_note = "No TE annotation required for yeast technical fixture."

    for temporary in ["_source.fa", "_source.gff3", "_source_te.gff3"]:
        path = out / temporary
        if path.exists():
            path.unlink()
    write_checksums(out)
    (out / "README.md").write_text(
        f"# {args.dataset}\n\n"
        f"{spec['readme']}\n\n"
        f"- Prepared: {date.today().isoformat()}\n"
        f"- Region: {region}\n"
        f"- Output chromosome name: {out_chrom}\n"
        f"- Download used: {args.download}\n"
        f"- Source directory: {source or 'download'}\n"
        f"- TE/mappability/blacklist note: {te_note}\n\n"
        "Full genomes should not be committed to this repository.\n",
        encoding="utf-8",
    )
    print(f"Prepared {args.dataset} in {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare small optional real reference fixtures.")
    parser.add_argument("--dataset", required=True, choices=sorted(DATASETS))
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--no-download", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--source-dir", type=Path)
    parser.add_argument("--region", default=None)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.download and args.no_download:
        raise SystemExit("Use only one of --download or --no-download.")
    if not args.download:
        args.no_download = True
    prepare_dataset(args)


if __name__ == "__main__":
    main()
