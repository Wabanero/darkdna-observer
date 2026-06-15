"""Physical susceptibility view composed from conservative sequence proxies."""

from __future__ import annotations

import numpy as np


def compute_physical_susceptibility_view(row: dict) -> dict[str, float]:
    g4 = float(row.get("G4_susceptibility_proxy", row.get("predicted_G_quadruplex_proxy_score", 0.0)) or 0.0)
    oxidation = float(row.get("G_tract_density", 0.0) or 0.0) + max(0.0, float(row.get("G_skew", 0.0) or 0.0))
    nonb = float(row.get("non_B_DNA_aggregate_score", 0.0) or 0.0)
    rloop = float(row.get("R_loop_forming_potential", 0.0) or 0.0)
    fork_texture = float(
        np.mean(
            [
                row.get("simple_repeat_fraction", 0.0) or 0.0,
                row.get("palindrome_density", 0.0) or 0.0,
                row.get("homopolymer_run_p95", 0.0) / max(1.0, row.get("length", 1.0) or 1.0),
                row.get("AT_skew", 0.0) or 0.0,
                row.get("spacing_periodicity_autocorrelation", 0.0) or 0.0,
            ]
        )
    )
    charge = float(np.mean([g4, oxidation, nonb]))
    aggregate = float(np.mean([g4, oxidation, nonb, rloop, fork_texture]))
    return {
        "physical_view_G4_susceptibility": g4,
        "oxidation_susceptibility": oxidation,
        "physical_view_non_B_DNA_propensity": nonb,
        "physical_view_R_loop_propensity": rloop,
        "replication_fork_texture": fork_texture,
        "charge_oxidation_susceptibility_score": charge,
        "fork_texture_score": fork_texture,
        "nonB_physical_susceptibility_score": aggregate,
    }
