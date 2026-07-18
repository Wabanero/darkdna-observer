"""Mutation-spectrum records for evolutionary null simulation."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class MutationSpectrum:
    transition_rate: float = 1.0
    transversion_rate: float = 0.5
    cpg_deamination_multiplier: float = 4.0
    source: str = "generic_relative_rates"
    calibration_status: str = "generic_not_organism_specific"


def mutation_spectrum_from_table(table: pd.DataFrame | None) -> MutationSpectrum:
    if table is None or table.empty or not {"mutation_class", "relative_rate"}.issubset(table.columns):
        return MutationSpectrum()
    rates = dict(zip(table["mutation_class"].astype(str), pd.to_numeric(table["relative_rate"], errors="coerce")))
    return MutationSpectrum(
        transition_rate=float(rates.get("transition", 1.0)),
        transversion_rate=float(rates.get("transversion", 0.5)),
        cpg_deamination_multiplier=float(rates.get("CpG_deamination", 4.0)),
        source="user_supplied_mutation_spectrum",
        calibration_status="organism_conditioned_partial",
    )
