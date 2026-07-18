"""DNA-shape predictor adapters with explicit availability and provenance.

The historical dinucleotide table remains available only under ``legacy_*``
screening names.  It is never substituted for a validated DNAshapeR or
dnacurve value under the strong predictor-backed field names.
"""

from __future__ import annotations

import json
import math
from functools import lru_cache
from importlib import metadata
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from darkdna.utils.optional_deps import optional_import
from .sequence import clean_sequence


DINUC_BENDABILITY = {
    "AA": 0.92,
    "AT": 0.82,
    "TA": 1.05,
    "CA": 1.10,
    "GT": 1.10,
    "CT": 1.02,
    "GA": 1.02,
    "CG": 0.78,
    "GC": 1.20,
    "GG": 0.95,
    "CC": 0.95,
    "AC": 0.96,
    "AG": 0.98,
    "TC": 0.98,
    "TG": 0.96,
    "TT": 0.92,
}


def dinuc_values(seq: str, table: dict[str, float]) -> list[float]:
    seq = clean_sequence(seq)
    return [table[seq[index : index + 2]] for index in range(len(seq) - 1) if seq[index : index + 2] in table]


def _version(package: str) -> str:
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return "unavailable"


def _record(
    feature: str,
    *,
    predictor: str,
    predictor_version: str,
    raw_value: float | None = None,
    normalization: str = "none",
    availability_status: str = "unavailable",
    warning: str,
) -> dict[str, float | str]:
    return {
        "feature": feature,
        "predictor": predictor,
        "predictor_version": predictor_version,
        "raw_value": float(raw_value) if raw_value is not None and np.isfinite(raw_value) else None,
        "normalization": normalization,
        "availability_status": availability_status,
        "warning": warning,
    }


@lru_cache(maxsize=1)
def _predictor_availability() -> dict[str, str | bool]:
    """Discover local predictors once; no package is downloaded automatically."""

    rpy2, rpy2_warning = optional_import("rpy2")
    dna_shape_available = False
    dna_shape_reason = rpy2_warning or "rpy2 is unavailable."
    if rpy2 is not None:  # pragma: no cover - depends on local R installation.
        try:
            packages, _ = optional_import("rpy2.robjects.packages")
            dna_shape_available = bool(packages and packages.isinstalled("DNAshapeR"))
            dna_shape_reason = (
                "DNAshapeR is installed, but batch prediction requires a configured local adapter and is unavailable in this sequence-only call."
                if dna_shape_available
                else "rpy2 is installed but the local R package DNAshapeR is not installed."
            )
        except Exception as exc:
            dna_shape_reason = f"Could not inspect the local DNAshapeR installation: {exc}"
    _, dnacurve_warning = optional_import("dnacurve")
    return {
        "dnashaper_available": dna_shape_available,
        "dnashaper_reason": dna_shape_reason,
        "dnashaper_version": _version("rpy2") + "+local_R_package",
        "dnacurve_reason": dnacurve_warning or "dnacurve is installed but no validated aggregate adapter is configured for this call.",
        "dnacurve_version": _version("dnacurve"),
    }


def _run_dnashaper(seq: str) -> tuple[dict[str, float], str]:
    """Invoke a locally installed DNAshapeR predictor through rpy2."""

    try:  # pragma: no cover - depends on a local R + DNAshapeR installation.
        robjects, _ = optional_import("rpy2.robjects")
        packages, _ = optional_import("rpy2.robjects.packages")
        if robjects is None or packages is None:
            return {}, "rpy2 predictor modules are unavailable."
        dna_shape_r = packages.importr("DNAshapeR")
        with TemporaryDirectory(prefix="darkdna_dnashaper_") as temp_dir:
            fasta = Path(temp_dir) / "sequence.fa"
            fasta.write_text(f">sequence\n{seq}\n", encoding="ascii")
            predicted = dna_shape_r.getShape(str(fasta), shapeType="All")
        names = list(predicted.names or [])
        aliases = {
            "MGW": "minor_groove_width",
            "Roll": "roll",
            "ProT": "propeller_twist",
            "HelT": "helix_twist",
            "EP": "electrostatic_potential_proxy",
        }
        values: dict[str, float] = {}
        for r_name, output_name in aliases.items():
            if r_name not in names:
                continue
            raw = np.asarray(predicted.rx2(r_name), dtype=float).ravel()
            finite = raw[np.isfinite(raw)]
            if finite.size:
                values[output_name] = float(finite.mean())
        if not values:
            return {}, "DNAshapeR returned no finite recognized shape arrays."
        return values, "DNAshapeR local prediction completed; raw_value is the finite per-position mean."
    except Exception as exc:
        return {}, f"DNAshapeR local prediction failed: {exc}"


def _validated_shape_records(seq: str) -> list[dict[str, float | str]]:
    availability = _predictor_availability()
    predicted: dict[str, float] = {}
    prediction_reason = str(availability["dnashaper_reason"])
    if bool(availability["dnashaper_available"]):
        predicted, prediction_reason = _run_dnashaper(seq)
    shape_features = {
        "minor_groove_width": "DNAshapeR",
        "roll": "DNAshapeR",
        "propeller_twist": "DNAshapeR",
        "helix_twist": "DNAshapeR",
        "electrostatic_potential_proxy": "DNAshapeR_EP",
        "nucleosome_affinity": "user_supplied_validated_predictor",
        "curvature": "dnacurve",
        "bendability": "user_supplied_validated_predictor",
    }
    records: list[dict[str, float | str]] = []
    for feature, predictor in shape_features.items():
        if predictor.startswith("DNAshapeR"):
            reason = prediction_reason
            version = str(availability["dnashaper_version"])
            raw_value = predicted.get(feature, math.nan)
            status = "available" if np.isfinite(raw_value) else "unavailable"
        elif predictor == "dnacurve":
            reason = str(availability["dnacurve_reason"])
            version = str(availability["dnacurve_version"])
            raw_value = None
            status = "unavailable"
        else:
            reason = "No validated local predictor or model path was supplied."
            version = "unavailable"
            raw_value = None
            status = "unavailable"
        records.append(
            _record(
                feature,
                predictor=predictor,
                predictor_version=version,
                raw_value=raw_value,
                normalization="finite_per_position_mean" if status == "available" else "unavailable",
                availability_status=status,
                warning=reason,
            )
        )
    return records


def compute_physical_shape_features(seq: str) -> dict[str, float | str]:
    seq = clean_sequence(seq)
    legacy_values = dinuc_values(seq, DINUC_BENDABILITY)
    legacy_bendability = float(np.mean(legacy_values)) if legacy_values else math.nan
    legacy_curvature_variability = float(np.std(legacy_values)) if legacy_values else math.nan
    records = _validated_shape_records(seq)
    by_feature = {str(record["feature"]): record for record in records}
    warning = (
        "Deprecated compatibility alias from an unvalidated dinucleotide screening table. "
        "It is not a DNAshapeR/dnacurve prediction and must not be interpreted as measured DNA shape."
    )
    result: dict[str, float | str] = {
        "dna_shape_records_json": json.dumps(records, sort_keys=True),
        "dna_shape_availability_status": "available" if any(record["availability_status"] == "available" for record in records) else "unavailable",
        "dna_shape_method": "validated_predictor_only_no_silent_heuristic_substitution_v2",
        "legacy_dinucleotide_bendability_screen": legacy_bendability,
        "legacy_dinucleotide_curvature_variability_screen": legacy_curvature_variability,
        "legacy_dinucleotide_shape_warning": warning,
        # Compatibility aliases: explicit warnings accompany every value.
        "DNA_bendability_proxy": legacy_bendability,
        "DNA_curvature_proxy": legacy_curvature_variability,
        "DNA_stiffness_proxy": math.nan,
        "DNA_bendability_proxy_deprecation_warning": warning,
        "DNA_curvature_proxy_deprecation_warning": warning,
        "DNA_stiffness_proxy_deprecation_warning": (
            "Unavailable. The former stiffness=2-bendability transformation was scientifically unsupported and has been retired."
        ),
        "physical_shape_warning": "Validated shape values are NA when their local predictor is unavailable; legacy screens are separately named.",
    }
    output_names = {
        "minor_groove_width": "minor_groove_width",
        "roll": "roll",
        "propeller_twist": "propeller_twist",
        "helix_twist": "helix_twist",
        "electrostatic_potential_proxy": "electrostatic_potential_proxy",
        "nucleosome_affinity": "nucleosome_affinity",
        "curvature": "predicted_DNA_curvature",
        "bendability": "predicted_DNA_bendability",
    }
    for feature, output_name in output_names.items():
        record = by_feature[feature]
        result[output_name] = float(record["raw_value"]) if record["raw_value"] is not None else math.nan
        result[f"{output_name}_predictor"] = str(record["predictor"])
        result[f"{output_name}_predictor_version"] = str(record["predictor_version"])
        result[f"{output_name}_normalization"] = str(record["normalization"])
        result[f"{output_name}_availability_status"] = str(record["availability_status"])
        result[f"{output_name}_warning"] = str(record["warning"])
    # Historical *_proxy fields now point to the strong, unavailable values,
    # not to the legacy table.
    result["minor_groove_width_proxy"] = result["minor_groove_width"]
    result["helix_twist_proxy"] = result["helix_twist"]
    result["propeller_twist_proxy"] = result["propeller_twist"]
    result["roll_proxy"] = result["roll"]
    return result
