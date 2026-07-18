"""Schemas and scientific interpretation constants for Mode B."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


EVOLUTIONARY_CAVEAT = (
    "A causal or quantity-dependent effect does not establish that the element "
    "originated or is maintained by selection for that effect."
)

MODEL_CAVEAT = (
    "Transformation sensitivities are model-based perturbation evidence, not biological causality."
)


@dataclass(frozen=True)
class AvailabilityRecord:
    status: str
    value: float | None = None
    reason: str = ""
    source: str = ""
    limitations: tuple[str, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


MODE_B_CANDIDATE_LABELS = (
    "sequence_specific_architecture_candidate",
    "sequence_indifferent_architecture_candidate",
    "mixed_sequence_quantity_candidate",
    "length_constrained_interval_candidate",
    "copy_number_constraint_candidate",
    "spacing_constraint_candidate",
    "topological_spacer_candidate",
    "bulk_heterochromatin_candidate",
    "occupancy_dependent_region_candidate",
    "genome_quantity_difference_maker_candidate",
    "artifact_compatible_candidate",
    "unresolved_architecture_candidate",
)
