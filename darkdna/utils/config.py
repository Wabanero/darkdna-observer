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

        def to_dict(self) -> dict[str, Any]:
            return self.model_dump()

except Exception:  # pragma: no cover - covered where pydantic unavailable.

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

        @classmethod
        def model_validate(cls, data: dict[str, Any]) -> "ObserverConfig":
            allowed = cls().__dict__.keys()
            payload = {k: v for k, v in data.items() if k in allowed}
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
