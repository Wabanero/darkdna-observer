"""Primitive-specific assay blueprints."""

from __future__ import annotations

from copy import deepcopy

from .ontology import get_primitive


KEY_INTERACTION = "effect = (Native_treatment - Native_control) - (ControlSequence_treatment - ControlSequence_control)"
TEMPORAL_INTERACTION = "Include sequence x treatment x time/history interaction."


ASSAY_BLUEPRINTS: dict[str, dict] = {
    "fractal_scaffold_candidate": {
        "assay": "Fractal Folding Assay",
        "sequence_controls": ["native", "GC-matched shuffled", "dinucleotide-preserved shuffled", "k-mer-preserved shuffled", "scale-shuffled", "repeat-scrambled"],
        "perturbations": ["salt", "Mg2+", "temperature shift", "crowding agent", "torsional stress"],
        "readouts": ["DNA compaction", "single-molecule extension", "AFM shape", "nucleosome assembly"],
        "classical_validation": ["tiling MPRA", "CRISPRi deletion", "ATAC/CUT&Tag before/after perturbation"],
    },
    "constraint_grammar_region_candidate": {
        "assay": "Grammar Scramble Tiling Assay",
        "sequence_controls": ["native", "motif-preserved spacing-randomized", "spacing-preserved motif-mutated", "orientation-inverted", "phase-shifted +5/+10 bp", "palindrome-broken", "forbidden-word inserted"],
        "perturbations": ["mechanical stretch", "osmotic stress", "heat pulse", "mild oxidative stress"],
        "readouts": ["DNA shape", "nucleosome positioning", "protein-free folding", "condensate partitioning"],
        "classical_validation": ["tiling MPRA", "synthetic reporter library"],
    },
    "quantum_susceptible_domain_candidate": {
        "assay": "Charge-Oxidation Susceptibility Assay",
        "sequence_controls": ["native", "G-tract disrupted", "G4-disrupting mutant", "GC-matched control", "reverse complement", "methylated native"],
        "perturbations": ["mild oxidative stress", "temperature shift", "altered ionic condition", "UV/photosensitized oxidation if appropriate"],
        "readouts": ["8-oxoG formation", "G4 folding/unfolding", "nanopore dwell-time anomaly", "charge-transfer proxy"],
        "classical_validation": ["G4-seq/G4-ChIP-seq overlap if available", "oxidative lesion mapping"],
    },
    "replication_instability_candidate": {
        "assay": "Fork-Texture Assay",
        "sequence_controls": ["native", "scrambled", "repeat-disrupted", "palindrome-broken", "nearby neutral genomic control"],
        "perturbations": ["mild replication stress", "low-dose aphidicolin/hydroxyurea if appropriate", "oxidative replication stress", "recovery"],
        "readouts": ["fork speed", "fork pausing", "origin firing", "replication timing", "fork asymmetry"],
        "classical_validation": ["replication timing tracks", "Repli-seq if available"],
    },
    "chromatin_motion_oscillator_candidate": {
        "assay": "Live-Locus Motion Assay",
        "sequence_controls": ["endogenous native locus", "matched genomic control", "safe-harbor native insertion", "phase-scrambled insertion"],
        "perturbations": ["osmotic stress", "mechanical stretch", "temperature pulse", "cell-state perturbation", "recovery"],
        "readouts": ["mean squared displacement", "confinement radius", "subdiffusion exponent", "locus-nuclear landmark distance"],
        "classical_validation": ["Hi-C/Micro-C", "ATAC/CUT&Tag context if available"],
    },
    "decoherence_boundary_candidate": {
        "assay": "Noise Propagation Barrier Assay",
        "sequence_controls": ["reporterA-native-reporterB", "reporterA-scrambled-reporterB", "reporterA-neutral-spacer-reporterB", "known insulator control"],
        "perturbations": ["transcriptional pulse", "mild stress", "oscillatory stimulus", "recovery"],
        "readouts": ["reporter variance", "covariance", "burst synchrony", "mutual information", "noise transfer coefficient"],
        "classical_validation": ["single-cell reporter", "Perturb-seq", "scATAC/scRNA variance readouts"],
    },
    "hysteresis_candidate": {
        "assay": "Hysteresis Sequence Assay",
        "sequence_controls": ["native", "scrambled", "repeat-disrupted", "palindrome-broken", "GC-matched control"],
        "perturbations": ["treatment A", "recovery", "treatment A again", "mild-to-strong stress", "strong-to-mild stress"],
        "readouts": ["structural state", "reporter dynamics", "single-cell variance", "locus motion", "replication timing memory"],
        "classical_validation": ["time-course ATAC/RNA after perturbation and recovery"],
        "temporal": True,
    },
    "resonant_pulse_decoder_candidate": {
        "assay": "Temporal Pulse Decoding Assay",
        "sequence_controls": ["native", "phase-scrambled", "periodicity-disrupted", "GC-matched", "motif-preserved spacing-randomized"],
        "perturbations": ["equal total dose continuous stimulus", "6x10 min pulses", "12x5 min pulses", "random pulses", "low/high frequency oscillatory stimulus"],
        "readouts": ["reporter burst timing", "structural switching", "locus motion", "noise propagation", "cell-state transition probability"],
        "classical_validation": ["live reporter time-course", "single-cell time-course"],
        "temporal": True,
    },
    "possibility_gate_candidate": {
        "assay": "Future-State Bias Assay",
        "sequence_controls": ["native", "scrambled", "deletion/no-sequence", "GC-matched", "neutral control"],
        "perturbations": ["weak differentiation stimulus", "subthreshold stress", "full treatment", "treatment followed by recovery"],
        "readouts": ["future-state probability", "transition time", "state stability", "state entropy"],
        "classical_validation": ["single-cell fate mapping", "pseudotime/state transition analysis"],
        "temporal": True,
    },
    "criticality_tuner_candidate": {
        "assay": "Transition-Threshold Assay",
        "sequence_controls": ["native", "scrambled", "phase-disrupted", "deletion/no-sequence", "neutral control"],
        "perturbations": ["dose gradient", "pulse-duration gradient"],
        "readouts": ["transition probability", "bifurcation point", "critical slowing down", "recovery rate"],
        "classical_validation": ["dose-response single-cell assay"],
        "temporal": True,
    },
    "negative_space_element_candidate": {
        "assay": "Negative-Space Rescue/Scramble Assay",
        "sequence_controls": ["native", "forbidden-word-inserted", "depleted-kmer-rescued", "GC-matched shuffled", "repeat-inserted"],
        "perturbations": ["baseline", "mild stress", "recovery"],
        "readouts": ["structural change", "reporter variance", "nucleosome positioning", "local chromatin accessibility", "sequence grammar disruption"],
        "classical_validation": ["tiling reporter with inserted/depleted tokens"],
    },
    "sequence_regime_boundary_candidate": {
        "assay": "Boundary Disruption Assay",
        "sequence_controls": ["native", "boundary-smoothed", "left-half duplicated", "right-half duplicated", "boundary-shuffled", "neutral spacer"],
        "perturbations": ["baseline", "stress pulse", "mechanical/osmotic perturbation"],
        "readouts": ["noise transfer", "nucleosome positioning", "reporter covariance", "chromatin accessibility boundary shift"],
        "classical_validation": ["ATAC/CUT&Tag boundary profiling", "reporter insulation assay"],
    },
    "TE_grammar_node_candidate": {
        "assay": "TE Grammar Reconstruction Assay",
        "sequence_controls": ["native TE mosaic", "TE-order scrambled", "TE-orientation inverted", "TE-boundary deleted", "family-matched random TE fragments"],
        "perturbations": ["baseline", "stress", "recovery", "cell-state perturbation"],
        "readouts": ["reporter response", "chromatin state", "noise", "folding", "accessibility", "TE-family-specific behavior"],
        "classical_validation": ["TE-overlap enrichment", "TE-derived cCRE comparison", "CRISPRi of TE-derived candidate"],
    },
    "unexplained_dark_anomaly_candidate": {
        "assay": "Matched-Null Prioritization Assay",
        "sequence_controls": ["native", "GC-matched shuffled", "dinucleotide-preserved shuffled", "nearby neutral genomic control"],
        "perturbations": ["baseline", "mild stress", "recovery"],
        "readouts": ["reporter response", "chromatin accessibility", "sequence-shape proxy", "single-cell variance"],
        "classical_validation": ["matched-null reanalysis", "artifact and mappability review"],
    },
}


def recommend_assay(primitive: str) -> dict:
    ontology = get_primitive(primitive)
    candidate_name = ontology.candidate_name
    blueprint = deepcopy(ASSAY_BLUEPRINTS.get(candidate_name, ASSAY_BLUEPRINTS["unexplained_dark_anomaly_candidate"]))
    blueprint["candidate_name"] = ontology.candidate_name
    blueprint["confirmed_name"] = ontology.confirmed_name
    blueprint["requires_dynamic_data"] = ontology.requires_dynamic_data
    blueprint["required_validation_data"] = ontology.required_input_level
    blueprint["prompt1_allowed_interpretation"] = ontology.prompt1_interpretation
    blueprint["suggested_prompt2_view"] = ontology.suggested_prompt2_view
    blueprint["key_interaction_test"] = KEY_INTERACTION
    if blueprint.get("temporal"):
        blueprint["temporal_interaction_test"] = TEMPORAL_INTERACTION
    blueprint.setdefault("allowed_interpretation", "A positive result supports a sequence-dependent, perturbation-sensitive hypothesis requiring replication and orthogonal validation.")
    blueprint.setdefault("forbidden_interpretation", "Do not claim confirmed function, quantum effects, teleology, or mechanism from this assay blueprint alone.")
    return blueprint
