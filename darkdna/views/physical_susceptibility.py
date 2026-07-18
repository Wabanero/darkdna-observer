"""Structure-specific physical-susceptibility screening views.

Unlike non-B structures are not averaged into a mechanistic score.  The legacy
aggregate is retained only as an explicitly deprecated priority-screen alias.
"""

from __future__ import annotations

import numpy as np


def compute_physical_susceptibility_view(row: dict) -> dict[str, float]:
    g4 = float(row.get("G4_sequence_potential", row.get("G4_susceptibility_proxy", row.get("predicted_G_quadruplex_proxy_score", 0.0))) or 0.0)
    i_motif = float(row.get("i_motif_sequence_potential", 0.0) or 0.0)
    z_dna = float(row.get("Z_DNA_sequence_potential", row.get("Z_DNA_propensity_proxy", 0.0)) or 0.0)
    triplex = float(row.get("triplex_H_DNA_sequence_potential", row.get("triplex_H_DNA_proxy", 0.0)) or 0.0)
    cruciform = float(row.get("cruciform_sequence_potential", row.get("cruciform_forming_potential", 0.0)) or 0.0)
    slipped = float(row.get("slipped_DNA_sequence_potential", 0.0) or 0.0)
    oxidation = float(row.get("G_tract_density", 0.0) or 0.0) + max(0.0, float(row.get("G_skew", 0.0) or 0.0))
    rloop = float(row.get("R_loop_susceptibility_sequence_potential", row.get("R_loop_forming_potential", 0.0)) or 0.0)
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
    charge = float(np.mean([g4, oxidation]))
    # A maximum is retained solely for compatibility ranking.  It means "at
    # least one Level-1 screen is high", not one shared formation mechanism.
    legacy_priority = float(max(g4, i_motif, z_dna, triplex, cruciform, slipped, rloop, fork_texture))
    return {
        "physical_view_G4_susceptibility": g4,
        "physical_view_i_motif_sequence_potential": i_motif,
        "physical_view_Z_DNA_sequence_potential": z_dna,
        "physical_view_triplex_H_DNA_sequence_potential": triplex,
        "physical_view_cruciform_sequence_potential": cruciform,
        "physical_view_slipped_DNA_sequence_potential": slipped,
        "oxidation_susceptibility": oxidation,
        "physical_view_R_loop_propensity": rloop,
        "replication_fork_texture": fork_texture,
        "charge_oxidation_susceptibility_score": charge,
        "fork_texture_score": fork_texture,
        "nonB_physical_susceptibility_score": legacy_priority,
        "nonB_physical_susceptibility_score_status": "deprecated_legacy_priority_screen",
        "nonB_physical_susceptibility_score_warning": (
            "Maximum of separate Level-1 structure screens retained for compatibility; it is not a mechanistic aggregate, formation probability, or observed structure."
        ),
    }
