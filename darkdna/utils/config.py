"""Configuration loading and validation.

Pydantic is declared as a package dependency. The module also contains a small
fallback so tests and smoke paths can run in minimal environments before the
package is installed with all dependencies.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml


try:  # pragma: no cover - depends on environment.
    from pydantic import BaseModel, Field

    class ModeToggle(BaseModel):
        enabled: bool = False

    class AnalysisModesConfig(BaseModel):
        sequence_specific: ModeToggle = Field(default_factory=lambda: ModeToggle(enabled=True))
        sequence_indifferent: ModeToggle = Field(default_factory=ModeToggle)
        conformation: ModeToggle = Field(default_factory=ModeToggle)
        transcription_rna: ModeToggle = Field(default_factory=ModeToggle)
        dynamic: ModeToggle = Field(default_factory=ModeToggle)
        evolution: ModeToggle = Field(default_factory=ModeToggle)

    class ArchitectureAnalysisLevels(BaseModel):
        window: bool = True
        locus: bool = True
        repeat_family: bool = True
        chromosome: bool = False
        genome: bool = False

    class ArchitectureIntervalSources(BaseModel):
        candidate_loci: str | None = None
        repeats: str | None = None
        copy_number: str | None = None
        structural_variants: str | None = None
        presence_absence: str | None = None
        syntenic_intervals: str | None = None
        anchors: str | None = None
        chromatin_compartments: str | None = None
        heterochromatin: str | None = None
        replication_timing: str | None = None

    class ArchitectureNullConfig(BaseModel):
        matched_interval: bool = True
        equal_length_replacement: bool = True
        length_titration: bool = True
        copy_number_permutation: bool = True
        repeat_family_matched: bool = True
        genome_size_matched: bool = False
        phylogenetic: bool = False

    class SequenceIndifferentConfig(BaseModel):
        analysis_levels: ArchitectureAnalysisLevels = Field(default_factory=ArchitectureAnalysisLevels)
        interval_sources: ArchitectureIntervalSources = Field(default_factory=ArchitectureIntervalSources)
        phenotype_table: str | None = None
        null_models: ArchitectureNullConfig = Field(default_factory=ArchitectureNullConfig)
        transformation_replicates: int = 5
        kmer_size: int = 3

    class EvolutionConfig(BaseModel):
        pangenome_enabled: bool = False
        synteny_enabled: bool = False
        population_variation_enabled: bool = False
        evolutionary_null_enabled: bool = False
        phylogenetic_tree: str | None = None
        mutation_spectrum: str | None = None
        n_surrogates: int = 25

    class NullFrameworkConfig(BaseModel):
        block_size_bp: int = 100_000
        minimum_independent_blocks: int = 5
        n_controls: int | None = None
        agreement_z_threshold: float = 2.0

    class DefaultStateBenchmarkConfig(BaseModel):
        enabled: bool = False
        max_windows: int = 256
        local_block_size: int = 1_000
        kmer_size: int = 3
        foundation_model_paths: dict[str, str] = Field(default_factory=dict)

    class StatisticsConfig(BaseModel):
        cv_group_priority: list[str] = Field(
            default_factory=lambda: ["chromosome", "genomic_block", "haplotype", "strain"]
        )
        conformal_residuals: bool = True
        block_bootstrap: bool = True
        minimum_independent_blocks: int = 5

    class ObserverConfig(BaseModel):
        project_name: str = "darkdna-observer"
        random_seed: int = 13
        fasta: str | None = None
        chrom_sizes: str | None = None
        annotation: str | None = None
        blacklist: str | None = None
        te_annotation: str | None = None
        ccre: str | None = None
        enhancer: str | None = None
        promoter: str | None = None
        mappability: str | None = None
        assembly_gaps: str | None = None
        segmental_duplication: str | None = None
        centromere_telomere: str | None = None
        output_dir: str = "darkdna_run"
        window_sizes: list[int] = Field(default_factory=lambda: [200, 1000, 5000, 10000, 50000])
        step_fraction: float = 0.5
        exclude_coding_exons: bool = True
        promoter_bp: int = 1000
        n_null: int = 25
        enable_candidate_only_prompt1: bool = True
        enable_te_grammar: bool = True
        enable_gff3_parsing: bool = True
        plant_non_model_compatibility: bool = True
        artifact_thresholds: dict[str, float] = Field(default_factory=dict)
        primitive_thresholds: dict[str, float] = Field(default_factory=dict)
        analysis_modes: AnalysisModesConfig = Field(default_factory=AnalysisModesConfig)
        sequence_indifferent: SequenceIndifferentConfig = Field(default_factory=SequenceIndifferentConfig)
        evolution: EvolutionConfig = Field(default_factory=EvolutionConfig)
        null_models: NullFrameworkConfig = Field(default_factory=NullFrameworkConfig)
        default_state_benchmark: DefaultStateBenchmarkConfig = Field(default_factory=DefaultStateBenchmarkConfig)
        statistics: StatisticsConfig = Field(default_factory=StatisticsConfig)

        def to_dict(self) -> dict[str, Any]:
            return self.model_dump()

except Exception:  # pragma: no cover - covered where pydantic unavailable.

    @dataclass
    class ModeToggle:  # type: ignore[no-redef]
        enabled: bool = False

    @dataclass
    class AnalysisModesConfig:  # type: ignore[no-redef]
        sequence_specific: ModeToggle = field(default_factory=lambda: ModeToggle(enabled=True))
        sequence_indifferent: ModeToggle = field(default_factory=ModeToggle)
        conformation: ModeToggle = field(default_factory=ModeToggle)
        transcription_rna: ModeToggle = field(default_factory=ModeToggle)
        dynamic: ModeToggle = field(default_factory=ModeToggle)
        evolution: ModeToggle = field(default_factory=ModeToggle)

    @dataclass
    class ArchitectureAnalysisLevels:  # type: ignore[no-redef]
        window: bool = True
        locus: bool = True
        repeat_family: bool = True
        chromosome: bool = False
        genome: bool = False

    @dataclass
    class ArchitectureIntervalSources:  # type: ignore[no-redef]
        candidate_loci: str | None = None
        repeats: str | None = None
        copy_number: str | None = None
        structural_variants: str | None = None
        presence_absence: str | None = None
        syntenic_intervals: str | None = None
        anchors: str | None = None
        chromatin_compartments: str | None = None
        heterochromatin: str | None = None
        replication_timing: str | None = None

    @dataclass
    class ArchitectureNullConfig:  # type: ignore[no-redef]
        matched_interval: bool = True
        equal_length_replacement: bool = True
        length_titration: bool = True
        copy_number_permutation: bool = True
        repeat_family_matched: bool = True
        genome_size_matched: bool = False
        phylogenetic: bool = False

    @dataclass
    class SequenceIndifferentConfig:  # type: ignore[no-redef]
        analysis_levels: ArchitectureAnalysisLevels = field(default_factory=ArchitectureAnalysisLevels)
        interval_sources: ArchitectureIntervalSources = field(default_factory=ArchitectureIntervalSources)
        phenotype_table: str | None = None
        null_models: ArchitectureNullConfig = field(default_factory=ArchitectureNullConfig)
        transformation_replicates: int = 5
        kmer_size: int = 3

    @dataclass
    class EvolutionConfig:  # type: ignore[no-redef]
        pangenome_enabled: bool = False
        synteny_enabled: bool = False
        population_variation_enabled: bool = False
        evolutionary_null_enabled: bool = False
        phylogenetic_tree: str | None = None
        mutation_spectrum: str | None = None
        n_surrogates: int = 25

    @dataclass
    class NullFrameworkConfig:  # type: ignore[no-redef]
        block_size_bp: int = 100_000
        minimum_independent_blocks: int = 5
        n_controls: int | None = None
        agreement_z_threshold: float = 2.0

    @dataclass
    class DefaultStateBenchmarkConfig:  # type: ignore[no-redef]
        enabled: bool = False
        max_windows: int = 256
        local_block_size: int = 1_000
        kmer_size: int = 3
        foundation_model_paths: dict[str, str] = field(default_factory=dict)

    @dataclass
    class StatisticsConfig:  # type: ignore[no-redef]
        cv_group_priority: list[str] = field(
            default_factory=lambda: ["chromosome", "genomic_block", "haplotype", "strain"]
        )
        conformal_residuals: bool = True
        block_bootstrap: bool = True
        minimum_independent_blocks: int = 5

    @dataclass
    class ObserverConfig:  # type: ignore[no-redef]
        project_name: str = "darkdna-observer"
        random_seed: int = 13
        fasta: str | None = None
        chrom_sizes: str | None = None
        annotation: str | None = None
        blacklist: str | None = None
        te_annotation: str | None = None
        ccre: str | None = None
        enhancer: str | None = None
        promoter: str | None = None
        mappability: str | None = None
        assembly_gaps: str | None = None
        segmental_duplication: str | None = None
        centromere_telomere: str | None = None
        output_dir: str = "darkdna_run"
        window_sizes: list[int] = field(default_factory=lambda: [200, 1000, 5000, 10000, 50000])
        step_fraction: float = 0.5
        exclude_coding_exons: bool = True
        promoter_bp: int = 1000
        n_null: int = 25
        enable_candidate_only_prompt1: bool = True
        enable_te_grammar: bool = True
        enable_gff3_parsing: bool = True
        plant_non_model_compatibility: bool = True
        artifact_thresholds: dict[str, float] = field(default_factory=dict)
        primitive_thresholds: dict[str, float] = field(default_factory=dict)
        analysis_modes: AnalysisModesConfig = field(default_factory=AnalysisModesConfig)
        sequence_indifferent: SequenceIndifferentConfig = field(default_factory=SequenceIndifferentConfig)
        evolution: EvolutionConfig = field(default_factory=EvolutionConfig)
        null_models: NullFrameworkConfig = field(default_factory=NullFrameworkConfig)
        default_state_benchmark: DefaultStateBenchmarkConfig = field(default_factory=DefaultStateBenchmarkConfig)
        statistics: StatisticsConfig = field(default_factory=StatisticsConfig)

        @classmethod
        def model_validate(cls, data: dict[str, Any]) -> "ObserverConfig":
            allowed = cls().__dict__.keys()
            payload = {k: v for k, v in data.items() if k in allowed}
            modes = payload.get("analysis_modes")
            if isinstance(modes, dict):
                payload["analysis_modes"] = AnalysisModesConfig(
                    **{
                        key: ModeToggle(**value) if isinstance(value, dict) else value
                        for key, value in modes.items()
                    }
                )
            sequence_indifferent = payload.get("sequence_indifferent")
            if isinstance(sequence_indifferent, dict):
                sequence_indifferent = dict(sequence_indifferent)
                if isinstance(sequence_indifferent.get("analysis_levels"), dict):
                    sequence_indifferent["analysis_levels"] = ArchitectureAnalysisLevels(
                        **sequence_indifferent["analysis_levels"]
                    )
                if isinstance(sequence_indifferent.get("interval_sources"), dict):
                    sequence_indifferent["interval_sources"] = ArchitectureIntervalSources(
                        **sequence_indifferent["interval_sources"]
                    )
                if isinstance(sequence_indifferent.get("null_models"), dict):
                    sequence_indifferent["null_models"] = ArchitectureNullConfig(
                        **sequence_indifferent["null_models"]
                    )
                payload["sequence_indifferent"] = SequenceIndifferentConfig(**sequence_indifferent)
            for key, model in [
                ("evolution", EvolutionConfig),
                ("null_models", NullFrameworkConfig),
                ("default_state_benchmark", DefaultStateBenchmarkConfig),
                ("statistics", StatisticsConfig),
            ]:
                if isinstance(payload.get(key), dict):
                    payload[key] = model(**payload[key])
            cfg = cls(**payload)
            if not cfg.window_sizes or any(int(w) <= 0 for w in cfg.window_sizes):
                raise ValueError("window_sizes must contain positive integers")
            if not 0 < float(cfg.step_fraction) <= 1:
                raise ValueError("step_fraction must be in (0, 1]")
            return cfg

        def to_dict(self) -> dict[str, Any]:
            return asdict(self)


DEFAULT_ARTIFACT_THRESHOLDS = {
    "high_n_fraction": 0.2,
    "low_mappability": 0.5,
    "extreme_low_complexity": 0.55,
    "extreme_repeat_density": 0.65,
    "scaffold_edge_bp": 1000,
    "very_short_usable_fraction": 0.5,
}

DEFAULT_PRIMITIVE_THRESHOLDS = {
    "residual_zscore": 2.0,
    "matched_null_zscore": 2.0,
    "confidence_high": 0.75,
}


def load_config(path: str | Path | None = None, overrides: dict[str, Any] | None = None) -> ObserverConfig:
    data: dict[str, Any] = {}
    if path:
        text = Path(path).read_text(encoding="utf-8")
        data = yaml.safe_load(text) or {}
    if overrides:
        data.update({k: v for k, v in overrides.items() if v is not None})
    cfg = ObserverConfig.model_validate(data)
    cfg.artifact_thresholds = {**DEFAULT_ARTIFACT_THRESHOLDS, **(cfg.artifact_thresholds or {})}
    cfg.primitive_thresholds = {**DEFAULT_PRIMITIVE_THRESHOLDS, **(cfg.primitive_thresholds or {})}
    return cfg


def write_config_snapshot(config: ObserverConfig, path: str | Path) -> None:
    Path(path).write_text(yaml.safe_dump(config.to_dict(), sort_keys=False), encoding="utf-8")


def resolve_config_path(config_path: str | Path | None, value: str | Path | None) -> Path | None:
    if value is None or value == "":
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    cwd_path = Path.cwd() / path
    if cwd_path.exists():
        return cwd_path
    if config_path is not None:
        cfg_relative = Path(config_path).resolve().parent / path
        if cfg_relative.exists():
            return cfg_relative
    return cwd_path
