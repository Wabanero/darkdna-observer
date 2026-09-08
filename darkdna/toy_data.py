"""Deterministic Level A toy genome for tests and examples."""

from __future__ import annotations

import json
import random
from pathlib import Path

import pandas as pd

from darkdna.io.fasta import write_chrom_sizes, write_fasta
from darkdna.windows.make_windows import make_dark_windows, write_windows


def repeat_to_length(pattern: str, length: int) -> str:
    return (pattern * ((length // len(pattern)) + 1))[:length]


def random_dna(length: int, rng: random.Random, gc: float = 0.42) -> str:
    weights = [((1 - gc) / 2), (gc / 2), (gc / 2), ((1 - gc) / 2)]
    return "".join(rng.choices(["A", "C", "G", "T"], weights=weights, k=length))


def write_fai(records: dict[str, str], fasta_path: str | Path, fai_path: str | Path | None = None, line_bases: int = 80) -> Path:
    fasta = Path(fasta_path)
    fai = Path(fai_path) if fai_path else fasta.with_suffix(fasta.suffix + ".fai")
    offset = 0
    rows = []
    for name, seq in records.items():
        header = f">{name}\n"
        offset += len(header.encode("utf-8"))
        seq_offset = offset
        line_width = line_bases + 1
        full_lines, remainder = divmod(len(seq), line_bases)
        rows.append(f"{name}\t{len(seq)}\t{seq_offset}\t{line_bases}\t{line_width}")
        offset += full_lines * line_width + (remainder + 1 if remainder else 0)
    fai.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return fai


def build_toy_records(seed: int = 42) -> tuple[dict[str, str], list[dict[str, object]]]:
    rng = random.Random(seed)
    records = {
        "toy_chr1": list(random_dna(100_000, rng, gc=0.43)),
        "toy_chr2": list(random_dna(80_000, rng, gc=0.38)),
        "toy_scaffoldA": list(random_dna(50_000, rng, gc=0.47)),
    }
    candidates: list[dict[str, object]] = [
        {
            "region_label": "toy_fractal_scale",
            "chrom": "toy_chr1",
            "start": 10_000,
            "end": 12_000,
            "expected_candidate_type": "fractal_scaffold_candidate",
            "expected_high_score_family": "scale_fractal",
            "seq": repeat_to_length("ACGTACGTGGCC" + "ACGTGGCC" * 2 + "ATGC", 2_000),
        },
        {
            "region_label": "toy_constraint_grammar",
            "chrom": "toy_chr1",
            "start": 20_000,
            "end": 21_600,
            "expected_candidate_type": "constraint_grammar_region_candidate",
            "expected_high_score_family": "grammar",
            "seq": repeat_to_length("AAAACCCCCG" + "T" * 10 + "GGGGTTTTAA" + "A" * 10, 1_600),
        },
        {
            "region_label": "toy_g4_physical",
            "chrom": "toy_chr1",
            "start": 32_000,
            "end": 33_600,
            "expected_candidate_type": "non_B_DNA_physical_susceptibility_candidate",
            "expected_high_score_family": "physical_susceptibility",
            "seq": repeat_to_length("GGGAGGGCGGGTTGGGAA", 1_600),
        },
        {
            "region_label": "toy_replication_instability",
            "chrom": "toy_chr1",
            "start": 45_000,
            "end": 46_800,
            "expected_candidate_type": "replication_instability_candidate",
            "expected_high_score_family": "fork_texture_nonb",
            "seq": repeat_to_length("ATATATATCGCGCGATGCAT" + "GATC" + "ATGCAT", 1_800),
        },
        {
            "region_label": "toy_entropy_boundary",
            "chrom": "toy_chr2",
            "start": 8_000,
            "end": 10_000,
            "expected_candidate_type": "sequence_regime_boundary_candidate",
            "expected_high_score_family": "boundary",
            "seq": repeat_to_length("AAAAAAAAAATTTTTTTTTT", 1_000) + random_dna(1_000, rng, gc=0.62),
        },
        {
            "region_label": "toy_negative_space",
            "chrom": "toy_chr2",
            "start": 25_000,
            "end": 26_600,
            "expected_candidate_type": "negative_space_element_candidate",
            "expected_high_score_family": "negative_space",
            "seq": repeat_to_length("ATATTAATTAAATAAT", 1_600),
        },
        {
            "region_label": "toy_te_mosaic",
            "chrom": "toy_chr2",
            "start": 40_000,
            "end": 42_400,
            "expected_candidate_type": "TE_grammar_node_candidate",
            "expected_high_score_family": "TE_grammar",
            "seq": repeat_to_length("TGCATGCAAAAACCCCGGGGTTTT", 2_400),
        },
        {
            "region_label": "toy_resonant_periodic",
            "chrom": "toy_scaffoldA",
            "start": 5_000,
            "end": 6_800,
            "expected_candidate_type": "periodic_spacing_grammar_candidate",
            "expected_high_score_family": "periodicity_phase",
            "seq": repeat_to_length("G" + "A" * 9 + "C" + "T" * 9, 1_800),
        },
        {
            "region_label": "toy_hysteresis_proxy",
            "chrom": "toy_scaffoldA",
            "start": 18_000,
            "end": 20_200,
            "expected_candidate_type": "asymmetric_repeat_architecture_candidate",
            "expected_high_score_family": "asymmetry_repeat_nonb",
            "seq": repeat_to_length("GGGAGGG" * 4 + "ATATAT" * 3 + "CCCG" + "A" * 20, 2_200),
        },
    ]
    for candidate in candidates:
        chrom = str(candidate["chrom"])
        start = int(candidate["start"])
        end = int(candidate["end"])
        records[chrom][start:end] = list(str(candidate["seq"])[: end - start])
    return {chrom: "".join(seq) for chrom, seq in records.items()}, candidates


def write_expected_candidates(candidates: list[dict[str, object]], path: str | Path) -> None:
    columns = ["region_label", "chrom", "start", "end", "expected_candidate_type", "expected_high_score_family"]
    pd.DataFrame([{col: candidate[col] for col in columns} for candidate in candidates]).to_csv(path, sep="\t", index=False)


def make_toy_data(outdir: str | Path, seed: int = 42) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    records, candidates = build_toy_records(seed=seed)

    fasta = out / "toy_genome.fa"
    write_fasta(records, fasta)
    fai = write_fai(records, fasta)
    chrom_sizes = out / "toy.chrom.sizes"
    write_chrom_sizes({name: len(seq) for name, seq in records.items()}, chrom_sizes)

    gff3 = out / "toy_annotation.gff3"
    gff3.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "toy_chr1\ttoy\tgene\t1001\t7000\t.\t+\t.\tID=geneA;Name=ToyGeneA;biotype=protein_coding",
                "toy_chr1\ttoy\tmRNA\t1001\t7000\t.\t+\t.\tID=txA;Parent=geneA",
                "toy_chr1\ttoy\tfive_prime_UTR\t1001\t1200\t.\t+\t.\tID=utrA1;Parent=txA",
                "toy_chr1\ttoy\texon\t1201\t1600\t.\t+\t.\tID=exonA1;Parent=txA",
                "toy_chr1\ttoy\texon\t3000\t3500\t.\t+\t.\tID=exonA2;Parent=txA",
                "toy_chr1\ttoy\tthree_prime_UTR\t6800\t7000\t.\t+\t.\tID=utrA2;Parent=txA",
                "toy_chr1\ttoy\tgene\t31500\t35000\t.\t-\t.\tID=geneB;Name=ToyGeneB;biotype=lncRNA",
                "toy_chr1\ttoy\texon\t33000\t33400\t.\t-\t.\tID=exonB1;Parent=geneB",
                "toy_chr2\ttoy\tgene\t7800\t12000\t.\t+\t.\tID=geneC;Name=ToyGeneC;biotype=protein_coding",
                "toy_chr2\ttoy\texon\t7800\t8100\t.\t+\t.\tID=exonC1;Parent=geneC",
                "toy_chr2\ttoy\texon\t11200\t11600\t.\t+\t.\tID=exonC2;Parent=geneC",
                "toy_scaffoldA\ttoy\tgene\t17000\t22000\t.\t-\t.\tID=geneD;Name=ToyGeneD;biotype=protein_coding",
                "toy_scaffoldA\ttoy\texon\t20500\t21000\t.\t-\t.\tID=exonD1;Parent=geneD",
                "",
            ]
        ),
        encoding="utf-8",
    )

    te_gff3 = out / "toy_te_annotation.gff3"
    te_gff3.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "toy_chr2\ttoy\tLTR_retrotransposon\t40001\t41000\t.\t+\t.\tID=te1;Name=GypsyToy;Class=LTR;Family=GypsyToy;Superfamily=Gypsy;Strand=+;Divergence=0.12",
                "toy_chr2\ttoy\tDNA_transposon\t40600\t41600\t.\t-\t.\tID=te2;Name=hATToy;Class=DNA;Family=hATToy;Superfamily=hAT;Strand=-;Divergence=0.05",
                "toy_chr2\ttoy\tLINE_element\t41500\t42400\t.\t+\t.\tID=te3;Name=LINEToy;Class=LINE;Family=LINEToy;Superfamily=L1;Strand=+;Divergence=0.21",
                "scaffold_A\ttoy\tLTR_retrotransposon\t6601\t7250\t.\t+\t.\tID=te_scaffoldA;Name=GypsyScaffoldToy;Class=LTR;Family=GypsyScaffoldToy;Superfamily=Gypsy;Strand=+;Divergence=0.08",
                "",
            ]
        ),
        encoding="utf-8",
    )

    blacklist = out / "toy_blacklist.bed"
    blacklist.write_text(
        "toy_chr1\t72000\t72500\ttoy_blacklist_1\n"
        "toy_chr2\t60000\t60600\ttoy_blacklist_2\n"
        "toy_scaffoldA\t100\t900\ttoy_scaffold_edge_blacklist\n",
        encoding="utf-8",
    )
    mappability = out / "toy_mappability.bedGraph"
    mappability.write_text(
        "toy_chr1\t0\t70000\t1.0\n"
        "toy_chr1\t70000\t76000\t0.25\n"
        "toy_chr1\t76000\t100000\t1.0\n"
        "toy_chr2\t0\t80000\t1.0\n"
        "toy_scaffoldA\t0\t2000\t0.35\n"
        "toy_scaffoldA\t2000\t50000\t1.0\n"
        "contig_unplaced_01\t0\t5000\t0.25\n",
        encoding="utf-8",
    )

    expected_tsv = out / "expected_candidates.tsv"
    write_expected_candidates(candidates, expected_tsv)
    (out / "expected_candidates.json").write_text(json.dumps(candidates, indent=2, default=str), encoding="utf-8")
    (out / "README.md").write_text(
        "# DarkDNA-Observer Toy Data\n\n"
        "Deterministic Level A synthetic test data generated by `darkdna make-toy-data --seed 42`.\n"
        "The candidate intervals are artificial sequence signatures for schema and approximate-overlap tests; they are not biological claims.\n",
        encoding="utf-8",
    )

    # Backward-compatible aliases used by early MVP tests.
    write_fasta(records, out / "toy.fa")
    write_fai(records, out / "toy.fa")
    (out / "toy.gff3").write_text(gff3.read_text(encoding="utf-8"), encoding="utf-8")
    (out / "toy_te.gff3").write_text(te_gff3.read_text(encoding="utf-8"), encoding="utf-8")

    windows = make_dark_windows(
        fasta=fasta,
        window_sizes=[200, 1000],
        step_fraction=0.5,
        annotation_path=gff3,
        te_annotation_path=te_gff3,
        blacklist_path=blacklist,
        mappability_path=mappability,
        exclude_coding_exons=True,
    )
    paths = write_windows(windows, out)
    windows.to_csv(out / "toy_windows.tsv", sep="\t", index=False)
    return {
        "fasta": fasta,
        "chrom_sizes": chrom_sizes,
        "fai": fai,
        "annotation_gff3": gff3,
        "te_gff3": te_gff3,
        "blacklist": blacklist,
        "mappability": mappability,
        "expected_candidates": expected_tsv,
        "toy_windows": out / "toy_windows.tsv",
        **paths,
    }
