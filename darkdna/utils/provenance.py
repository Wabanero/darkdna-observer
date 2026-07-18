"""Run provenance writers."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from darkdna import __version__

from .config import ObserverConfig, write_config_snapshot
from .optional_deps import write_optional_dependency_report


def file_checksum(path: str | Path, chunk_size: int = 1024 * 1024) -> str | None:
    p = Path(path)
    if not p.exists() or p.is_dir():
        return None
    digest = hashlib.sha256()
    with p.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect_versions() -> str:
    packages = ["numpy", "pandas", "scipy", "sklearn", "pyarrow", "jinja2", "networkx", "typer"]
    lines = [f"python={sys.version}", f"platform={platform.platform()}"]
    for package in packages:
        try:
            mod = __import__(package)
            lines.append(f"{package}={getattr(mod, '__version__', 'unknown')}")
        except Exception as exc:
            lines.append(f"{package}=unavailable ({exc})")
    return "\n".join(lines) + "\n"


def write_provenance(
    outdir: str | Path,
    command: str,
    config: ObserverConfig,
    input_paths: Iterable[str | Path | None] = (),
) -> None:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    checksums = {str(p): file_checksum(p) for p in input_paths if p}
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "darkdna_observer_version": __version__,
        "scientific_method_contract": "v2_phase1",
        "random_seed": config.random_seed,
        "input_file_checksums": checksums,
    }
    (out / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    write_config_snapshot(config, out / "config_snapshot.yaml")
    (out / "software_versions.txt").write_text(collect_versions(), encoding="utf-8")
    with (out / "command_log.txt").open("a", encoding="utf-8") as handle:
        handle.write(command + "\n")
    write_optional_dependency_report(out / "optional_dependency_report.json")


def shell_command_version(command: list[str]) -> str:
    try:
        return subprocess.check_output(command, text=True, stderr=subprocess.STDOUT, timeout=10)
    except Exception as exc:  # pragma: no cover - external command availability varies.
        return f"unavailable: {exc}"
