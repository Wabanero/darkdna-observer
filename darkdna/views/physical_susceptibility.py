"""Structure-specific physical-susceptibility screening views.

Unlike non-B structures are not averaged into a mechanistic score.  The legacy
aggregate is retained only as an explicitly deprecated priority-screen alias.
Missing predictors stay NA and are omitted from means and maxima.
"""

from __future__ import annotations

import numpy as np

from darkdna.utils.stats import finite_max, finite_mean, optional_row_float


def compute_physical_susceptibility_view(row: dict) -> dict[str, float]:
    g4 = optional_row_float(
        row,
        "G4_sequence_potential",
        "G4_susceptibility_proxy",
        "predicted_G_quadruplex_proxy_score",
    )
    i_motif = optional_row_float(row, "i_motif_sequence_potential")
    z_dna = optional_row_float(row, "Z_DNA_sequence_potential", "Z_DNA_propensity_proxy")
    triplex = optional_row_float(row, "triplex_H_DNA_sequence_potential", "triplex_H_DNA_proxy")
    cruciform = optional_row_float(row, "cruciform_sequence_potential", "cruciform_forming_potential")
    slipped = optional_row_float(row, "slipped_DNA_sequence_potential")
    g_tract = optional_row_float(row, "G_tract_density")
    g_skew = optional_row_float(row, "G_skew")
    positive_skew = g_skew if np.isfinite(g_skew) and g_skew > 0 else float("nan")
    if np.isfinite(g_tract) and np.isfinite(positive_skew):
        oxidation = g_tract + positive_skew
    elif np.isfinite(g_tract):
        oxidation = g_tract
    else:
        oxidation = positive_skew
    rloop = optional_row_float(row, "R_loop_susceptibility_sequence_potential", "R_loop_forming_potential")
    homopolymer = optional_row_float(row, "homopolymer_run_p95")
    length = optional_row_float(row, "length")
    homopolymer_frac = homopolymer / length if np.isfinite(homopolymer) and np.isfinite(length) and length > 0 else float("nan")
    fork_texture = finite_mean(
        [
            optional_row_float(row, "simple_repeat_fraction"),
            optional_row_float(row, "palindrome_density"),
            homopolymer_frac,
            optional_row_float(row, "AT_skew"),
            optional_row_float(row, "spacing_periodicity_autocorrelation"),
        ]
    )
    charge = finite_mean([g4, oxidation])
    # A maximum is retained solely for compatibility ranking.  It means "at
    # least one Level-1 screen is high", not one shared formation mechanism.
    legacy_priority = finite_max(
        [g4, i_motif, z_dna, triplex, cruciform, slipped, rloop, fork_texture]
    )
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
