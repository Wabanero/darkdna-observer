"""Machine-readable severe-null registry with explicit availability checks."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable


@dataclass(frozen=True)
class NullModelSpec:
    null_model_id: str
    family: str
    execution_mode: str
    required_columns: tuple[str, ...] = ()
    required_input: str = ""
    description: str = ""


NULL_MODEL_SPECS = (
    NullModelSpec("same_length_gc_matched", "composition", "matched_table", ("gc_content",), description="Same-length controls matched on GC."),
    NullModelSpec("mononucleotide_preserving", "composition", "sequence_transform", required_input="sequence", description="Permutation preserving exact mononucleotide counts."),
    NullModelSpec("dinucleotide_preserving_shuffle", "composition", "sequence_transform", required_input="sequence", description="Eulerian sequence surrogate preserving dinucleotide transitions where feasible."),
    NullModelSpec("kmer_preserving_shuffle", "composition", "sequence_transform", required_input="sequence", description="Block surrogate preserving local k-mer content."),
    NullModelSpec("markov_chain_surrogate", "evolutionary_proxy", "sequence_transform", required_input="sequence", description="First-order Markov surrogate fitted to the supplied sequence."),
    NullModelSpec("local_genomic_matched", "genomic", "matched_table", ("chrom", "start"), description="Controls from independent genomic blocks, prioritising local context."),
    NullModelSpec("gene_tss_distance_matched", "genomic", "matched_table", ("distance_to_nearest_tss",), description="Controls matched on TSS distance."),
    NullModelSpec("te_family_matched", "repeat", "matched_table", ("TE_family",), description="Controls from the same TE family."),
    NullModelSpec("te_subfamily_matched", "repeat", "matched_table", ("TE_subfamily",), description="Controls from the same TE subfamily."),
    NullModelSpec("te_age_matched", "repeat", "matched_table", ("TE_age",), description="Controls matched on TE age or divergence."),
    NullModelSpec("repeat_density_matched", "repeat", "matched_table", ("local_TE_density",), description="Controls matched on local repeat burden."),
    NullModelSpec("chromatin_compartment_matched", "chromatin", "matched_table", ("chromatin_compartment",), description="Controls from the same chromatin compartment."),
    NullModelSpec("replication_timing_matched", "chromatin", "matched_table", ("replication_timing",), description="Controls matched on replication timing."),
    NullModelSpec("recombination_matched", "mutation", "matched_table", ("recombination_rate",), description="Controls matched on recombination context."),
    NullModelSpec("mutation_rate_matched", "mutation", "matched_table", ("mutation_rate",), description="Controls matched on regional mutation rate."),
    NullModelSpec("damage_environment_matched", "mutation", "matched_table", ("damage_environment",), description="Controls matched on measured damage environment."),
    NullModelSpec("mappability_matched", "technical", "matched_table", ("mappability",), description="Controls matched on mappability."),
    NullModelSpec("assembly_confidence_matched", "technical", "matched_table", ("assembly_confidence",), description="Controls matched on assembly confidence."),
    NullModelSpec("syntenic_ortholog", "comparative", "matched_table", ("synteny_group",), description="Syntenic ortholog controls across assemblies or species."),
    NullModelSpec("population_frequency_matched", "population", "matched_table", ("population_frequency",), description="Controls matched on population frequency."),
    NullModelSpec("copy_number_matched", "quantity", "matched_table", ("copy_number",), description="Controls matched on copy number."),
    NullModelSpec("presence_absence_matched", "quantity", "matched_table", ("presence_absence_frequency",), description="Controls matched on presence/absence frequency."),
    NullModelSpec("reversed_sequence", "orientation", "sequence_transform", required_input="sequence", description="Sequence reversal without complementation."),
    NullModelSpec("reverse_complement", "orientation", "sequence_transform", required_input="sequence", description="Reverse-complement orientation control."),
    NullModelSpec("synthetic_equal_composition", "composition", "sequence_transform", required_input="sequence", description="Independent synthetic sequence with the same base composition and length."),
    NullModelSpec("evolutionary_process_generated", "evolutionary", "generative", required_input="evolutionary_model", description="Context-conditioned neutral evolutionary-process surrogate."),
)


def null_model_registry() -> list[dict[str, object]]:
    return [asdict(spec) for spec in NULL_MODEL_SPECS]


def assess_null_availability(
    columns: Iterable[str] = (),
    *,
    sequence_available: bool = False,
    evolutionary_model_available: bool = False,
) -> list[dict[str, object]]:
    available_columns = set(columns)
    rows: list[dict[str, object]] = []
    for spec in NULL_MODEL_SPECS:
        missing_columns = [column for column in spec.required_columns if column not in available_columns]
        input_available = (
            not spec.required_input
            or (spec.required_input == "sequence" and sequence_available)
            or (spec.required_input == "evolutionary_model" and evolutionary_model_available)
        )
        available = not missing_columns and input_available
        if available:
            reason = "All declared inputs are available."
        elif missing_columns:
            reason = f"Missing required columns: {', '.join(missing_columns)}."
        else:
            reason = f"Missing required input: {spec.required_input}."
        rows.append(
            {
                **asdict(spec),
                "available": bool(available),
                "status": "available" if available else "unavailable",
                "reason": reason,
            }
        )
    return rows
