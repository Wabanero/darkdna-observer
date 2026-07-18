"""Primitive-specific assay blueprints."""

from __future__ import annotations

from copy import deepcopy

from .ontology import get_primitive


KEY_INTERACTION = "effect = (Native_treatment - Native_control) - (ControlSequence_treatment - ControlSequence_control)"
TEMPORAL_INTERACTION = "Include sequence x treatment x time/history interaction."


def _bridge(
    measured_feature: str,
    proposed_dynamic_phenotype: str,
    candidate_processes: list[str],
    bridge_validation_steps: list[str],
    *,
    status: str = "hypothesized_bridge_requires_validation",
) -> dict:
    return {
        "measured_feature": measured_feature,
        "proposed_dynamic_phenotype": proposed_dynamic_phenotype,
        "candidate_intermediate_processes": candidate_processes,
        "required_bridge_evidence": bridge_validation_steps,
        "bridge_status": status,
        "assay_scope_if_bridge_missing": (
            "exploratory_only; this assay can generate bridge evidence but is not direct validation "
            "of the proposed primitive phenotype until the intermediate physical or molecular process is shown."
        ),
        "direct_primitive_validation_allowed": False,
    }


MECHANISTIC_BRIDGES: dict[str, dict] = {
    "fractal_scaffold_candidate": _bridge(
        "diagnostically valid multiscale texture, surrogate-standardized DFA, calibrated compression, or numeric-walk texture",
        "folding, compaction, or multiscale scaffold behavior",
        [
            "polymer physics or coarse-grained DNA simulations",
            "DNA shape or bendability models",
            "nucleosome occupancy or phasing predictions",
            "in vitro nucleosome assembly",
            "single-molecule force or extension response",
        ],
        [
            "show that the sequence-derived scale feature changes predicted polymer or shape behavior",
            "test length-matched controls",
            "test GC-, dinucleotide-, and k-mer-preserved controls",
            "test scale-shuffled controls before claiming a folding mechanism",
        ],
    ),
    "constraint_grammar_region_candidate": _bridge(
        "motif spacing, grammar entropy, transition surprise, or forbidden-word depletion",
        "grammar-dependent structural, binding, or context-response behavior",
        [
            "TF or nucleosome occupancy grammar",
            "DNA shape grammar",
            "protein cooperative binding",
            "local chromatin accessibility grammar",
        ],
        [
            "separate motif content from spacing/order effects",
            "use motif-preserved spacing-randomized controls",
            "use spacing-preserved motif-mutated controls",
            "show a molecular readout changes before interpreting any dynamic phenotype",
        ],
    ),
    "non_B_DNA_physical_susceptibility_candidate": _bridge(
        "G-richness, G4/non-B-DNA propensity, oxidation-prone context, or charge/oxidation proxy",
        "physical susceptibility under oxidation, ionic, thermal, or structural perturbation",
        [
            "G4 or non-B-DNA folding",
            "oxidative lesion formation",
            "sequence-dependent electronic or nanopore dwell-time proxy",
            "methylation or chromatin-state modulation of physical susceptibility",
        ],
        [
            "validate G4/non-B structure or oxidation sensitivity directly",
            "compare G-tract disrupted and GC-matched controls",
            "avoid claims of actual quantum effects from sequence proxies",
        ],
    ),
    "replication_instability_candidate": _bridge(
        "repeat, palindrome, skew, non-B propensity, or fork-texture proxy",
        "fork pausing, replication stress, origin firing, or timing instability",
        [
            "polymerase fork pausing at repeats/non-B structures",
            "R-loop or secondary-structure formation",
            "replication timing context",
            "repair or checkpoint sensitivity",
        ],
        [
            "measure fork speed, pausing, timing, or stress response",
            "test repeat-disrupted and palindrome-broken controls",
            "separate sequence instability from mappability or assembly artifacts",
        ],
    ),
    "chromatin_motion_oscillator_candidate": _bridge(
        "spacing periodicity, bendability proxy, entropy asymmetry, or phase-like sequence texture",
        "locus motion, confinement, or pulse-sensitive spatial dynamics",
        [
            "nucleosome phasing",
            "DNA stiffness or bendability",
            "loop extrusion or tethering context",
            "mechanosensitive chromatin response",
        ],
        [
            "predict or measure nucleosome phasing/shape first",
            "compare phase-scrambled insertions",
            "show live-locus motion changes after bridge-level sequence perturbation",
        ],
    ),
    "decoherence_boundary_candidate": _bridge(
        "entropy cliff, regime boundary, compression boundary, or local feature void",
        "noise barrier or variance-propagation boundary",
        [
            "nucleosome occupancy transition",
            "TF/cofactor occupancy discontinuity",
            "chromatin accessibility boundary",
            "reporter coupling or insulation geometry",
        ],
        [
            "measure a molecular boundary before interpreting noise insulation",
            "use boundary-smoothed and neutral-spacer controls",
            "show variance/covariance effects are sequence-dependent and not reporter-position artifacts",
        ],
    ),
    "hysteresis_candidate": _bridge(
        "asymmetry, nested repeats, non-B propensity, or recurrent-kmer orientation",
        "history-dependent response or recovery memory",
        [
            "secondary-structure metastability",
            "chromatin-state persistence",
            "replication or repair history dependence",
            "transcription-factor cooperative state retention",
        ],
        [
            "show a sequence-dependent molecular state can persist through recovery",
            "include first-pulse vs second-pulse comparisons",
            "use repeat-disrupted and palindrome-broken controls",
        ],
    ),
    "resonant_pulse_decoder_candidate": _bridge(
        "10 bp/147 bp periodicity, spacing Fourier power, autocorrelation, or phase grammar",
        "frequency- or pulse-specific response",
        [
            "nucleosome phasing",
            "TF cooperative binding at phased spacing",
            "mechanosensitivity",
            "replication dynamics",
            "chromatin looping or enhancer-promoter contact timing",
        ],
        [
            "show periodicity affects an intermediate molecular process",
            "compare phase-scrambled and periodicity-disrupted controls",
            "use equal-total-dose pulse controls before claiming pulse decoding",
        ],
    ),
    "possibility_gate_candidate": _bridge(
        "boundary-condition score, negative-space boundary, or forbidden-word depletion",
        "future-state probability, reachable-state bias, or transition gating",
        [
            "chromatin accessibility gating",
            "TF cooperative occupancy",
            "3D contact permissiveness",
            "noise/variance modulation before state transition",
        ],
        [
            "show sequence affects an intermediate chromatin or binding state",
            "measure baseline and weak-stimulus state probabilities",
            "avoid state-bias claims without perturbation and pseudotime/fate evidence",
        ],
    ),
    "criticality_tuner_candidate": _bridge(
        "entropy/compression boundary, sequence-regime transition, or threshold-like texture",
        "dose threshold, transition slope, or recovery-rate change",
        [
            "chromatin accessibility threshold",
            "cooperative TF binding",
            "nucleosome barrier strength",
            "stress-response or replication-response threshold",
        ],
        [
            "measure dose-response or pulse-duration curves",
            "test phase-disrupted/deletion controls",
            "demonstrate an intermediate molecular threshold before claiming criticality",
        ],
    ),
    "negative_space_element_candidate": _bridge(
        "depleted k-mers, motif desert, repeat/CpG/G-tract desert, or structured absence",
        "absence-dependent chromatin, binding, noise, or accessibility behavior",
        [
            "avoidance of nucleosome-disfavoring words",
            "absence of TF-binding or repeat-associated tokens",
            "DNA shape produced by missing motifs",
            "local chromatin accessibility or silence from structured absence",
        ],
        [
            "insert depleted/forbidden tokens and test rescue/disruption",
            "use GC-, length-, and k-mer-preserved controls",
            "show the absence itself changes a molecular readout before dynamic interpretation",
        ],
    ),
    "sequence_regime_boundary_candidate": _bridge(
        "left/right sequence-regime difference, entropy boundary, GC/CpG/repeat shift, or compression change",
        "boundary behavior, insulation proxy, or local readout discontinuity",
        [
            "nucleosome-position boundary",
            "DNA shape discontinuity",
            "TF/cofactor occupancy transition",
            "chromatin accessibility boundary",
        ],
        [
            "smooth or shuffle the boundary while preserving flanking content",
            "measure a molecular boundary before claiming regulatory boundary behavior",
            "include left-half/right-half and neutral-spacer controls",
        ],
    ),
    "TE_grammar_node_candidate": _bridge(
        "TE mosaic, TE boundary score, TE orientation entropy, or TE overlap grammar",
        "TE-order-dependent chromatin, binding, accessibility, or stress-response behavior",
        [
            "TE-derived TF binding",
            "TE boundary nucleosome positioning",
            "orientation/order-dependent repeat grammar",
            "stress-responsive TE chromatin context",
        ],
        [
            "separate simple TE overlap from TE order/orientation effects",
            "use TE-order scrambled and family-matched random TE controls",
            "show a molecular readout changes before TE exaptation claims",
        ],
    ),
    "unexplained_dark_anomaly_candidate": _bridge(
        "high residual sequence anomaly without a dominant primitive-specific feature",
        "unknown sequence-dependent behavior",
        [
            "unknown; requires feature audit",
            "matched-null review",
            "artifact, assembly, and mappability review",
        ],
        [
            "identify the measured feature class driving the anomaly",
            "build matched controls before wet-lab interpretation",
            "treat assays as exploratory until a bridge is specified",
        ],
        status="bridge_missing_feature_audit_required",
    ),
    "no_call": _bridge(
        "no candidate feature crossed threshold",
        "none",
        ["not applicable"],
        ["do not prioritize wet-lab validation until candidate thresholds are met"],
        status="not_applicable_no_candidate",
    ),
}


ASSAY_BLUEPRINTS: dict[str, dict] = {
    "fractal_scaffold_candidate": {
        "assay": "Fractal Folding Assay",
        "sequence_controls": ["native", "GC-matched shuffled", "dinucleotide-preserved shuffled", "k-mer-preserved shuffled", "scale-shuffled", "repeat-scrambled"],
        "perturbations": ["salt", "Mg2+", "temperature shift", "crowding agent", "torsional stress"],
        "readouts": ["DNA compaction", "single-molecule extension", "AFM shape", "nucleosome assembly"],
        "classical_validation": ["tiling MPRA", "CRISPRi deletion", "ATAC/CUT&Tag before/after perturbation"],
        "key_interaction_test": "folding_scale_effect = (Native_perturbed_compaction - Native_baseline_compaction) - (ScaleShuffled_perturbed_compaction - ScaleShuffled_baseline_compaction)",
    },
    "constraint_grammar_region_candidate": {
        "assay": "Grammar Scramble Tiling Assay",
        "sequence_controls": ["native", "motif-preserved spacing-randomized", "spacing-preserved motif-mutated", "orientation-inverted", "phase-shifted +5/+10 bp", "palindrome-broken", "forbidden-word inserted"],
        "perturbations": ["mechanical stretch", "osmotic stress", "heat pulse", "mild oxidative stress"],
        "readouts": ["DNA shape", "nucleosome positioning", "protein-free folding", "condensate partitioning"],
        "classical_validation": ["tiling MPRA", "synthetic reporter library"],
        "key_interaction_test": "grammar_effect = (Native_readout - MotifPreservedSpacingRandomized_readout) - (SpacingPreservedMotifMutated_readout - NeutralControl_readout)",
    },
    "non_B_DNA_physical_susceptibility_candidate": {
        "assay": "Charge-Oxidation Susceptibility Assay",
        "sequence_controls": ["native", "G-tract disrupted", "G4-disrupting mutant", "GC-matched control", "reverse complement", "methylated native"],
        "perturbations": ["mild oxidative stress", "temperature shift", "altered ionic condition", "UV/photosensitized oxidation if appropriate"],
        "readouts": ["8-oxoG formation", "G4 folding/unfolding", "nanopore dwell-time anomaly", "physical susceptibility proxy"],
        "classical_validation": ["G4-seq/G4-ChIP-seq overlap if available", "oxidative lesion mapping"],
        "key_interaction_test": "physical_susceptibility_effect = (Native_oxidation_or_G4_signal - G4Disrupted_signal) - (GCMatchedControl_oxidation_or_G4_signal - GCMatchedControl_baseline)",
    },
    "replication_instability_candidate": {
        "assay": "Fork-Texture Assay",
        "sequence_controls": ["native", "scrambled", "repeat-disrupted", "palindrome-broken", "nearby neutral genomic control"],
        "perturbations": ["mild replication stress", "low-dose aphidicolin/hydroxyurea if appropriate", "oxidative replication stress", "recovery"],
        "readouts": ["fork speed", "fork pausing", "origin firing", "replication timing", "fork asymmetry"],
        "classical_validation": ["replication timing tracks", "Repli-seq if available"],
        "key_interaction_test": "fork_instability_effect = (Native_stress_fork_pause - Native_baseline_fork_pause) - (RepeatDisrupted_stress_fork_pause - RepeatDisrupted_baseline_fork_pause)",
    },
    "chromatin_motion_oscillator_candidate": {
        "assay": "Live-Locus Motion Assay",
        "sequence_controls": ["endogenous native locus", "matched genomic control", "safe-harbor native insertion", "phase-scrambled insertion"],
        "perturbations": ["osmotic stress", "mechanical stretch", "temperature pulse", "cell-state perturbation", "recovery"],
        "readouts": ["mean squared displacement", "confinement radius", "subdiffusion exponent", "locus-nuclear landmark distance"],
        "classical_validation": ["Hi-C/Micro-C", "ATAC/CUT&Tag context if available"],
        "key_interaction_test": "motion_response_effect = (Native_locus_motion_after_pulse - Native_locus_motion_baseline) - (PhaseScrambled_locus_motion_after_pulse - PhaseScrambled_locus_motion_baseline)",
    },
    "decoherence_boundary_candidate": {
        "assay": "Noise Propagation Barrier Assay",
        "sequence_controls": ["reporterA-native-reporterB", "reporterA-scrambled-reporterB", "reporterA-neutral-spacer-reporterB", "known insulator control"],
        "perturbations": ["transcriptional pulse", "mild stress", "oscillatory stimulus", "recovery"],
        "readouts": ["reporter variance", "covariance", "burst synchrony", "mutual information", "noise transfer coefficient"],
        "classical_validation": ["single-cell reporter", "Perturb-seq", "scATAC/scRNA variance readouts"],
        "key_interaction_test": "noise_barrier_effect = (Native_flanking_reporter_covariance_after_pulse - Native_baseline_covariance) - (ScrambledSpacer_after_pulse_covariance - ScrambledSpacer_baseline_covariance)",
    },
    "hysteresis_candidate": {
        "assay": "Hysteresis Sequence Assay",
        "sequence_controls": ["native", "scrambled", "repeat-disrupted", "palindrome-broken", "GC-matched control"],
        "perturbations": ["treatment A", "recovery", "treatment A again", "mild-to-strong stress", "strong-to-mild stress"],
        "readouts": ["structural state", "reporter dynamics", "single-cell variance", "locus motion", "replication timing memory"],
        "classical_validation": ["time-course ATAC/RNA after perturbation and recovery"],
        "key_interaction_test": "hysteresis_effect = (Native_second_pulse_response - Native_first_pulse_response) - (Scrambled_second_pulse_response - Scrambled_first_pulse_response)",
        "temporal": True,
    },
    "resonant_pulse_decoder_candidate": {
        "assay": "Temporal Pulse Decoding Assay",
        "sequence_controls": ["native", "phase-scrambled", "periodicity-disrupted", "GC-matched", "motif-preserved spacing-randomized"],
        "perturbations": ["equal total dose continuous stimulus", "6x10 min pulses", "12x5 min pulses", "random pulses", "low/high frequency oscillatory stimulus"],
        "readouts": ["reporter burst timing", "structural switching", "locus motion", "noise propagation", "cell-state transition probability"],
        "classical_validation": ["live reporter time-course", "single-cell time-course"],
        "key_interaction_test": "pulse_decoding_effect = (Native_frequency_specific_response - Native_equal_dose_continuous_response) - (PeriodicityDisrupted_frequency_response - PeriodicityDisrupted_continuous_response)",
        "temporal": True,
    },
    "possibility_gate_candidate": {
        "assay": "Future-State Bias Assay",
        "sequence_controls": ["native", "scrambled", "deletion/no-sequence", "GC-matched", "neutral control"],
        "perturbations": ["weak differentiation stimulus", "subthreshold stress", "full treatment", "treatment followed by recovery"],
        "readouts": ["future-state probability", "transition time", "state stability", "state entropy"],
        "classical_validation": ["single-cell fate mapping", "pseudotime/state transition analysis"],
        "key_interaction_test": "state_bias_effect = (Native_future_state_probability_under_weak_stimulus - Native_baseline_probability) - (Deletion_or_scrambled_future_state_probability - Deletion_or_scrambled_baseline_probability)",
        "temporal": True,
    },
    "criticality_tuner_candidate": {
        "assay": "Transition-Threshold Assay",
        "sequence_controls": ["native", "scrambled", "phase-disrupted", "deletion/no-sequence", "neutral control"],
        "perturbations": ["dose gradient", "pulse-duration gradient"],
        "readouts": ["transition probability", "bifurcation point", "critical slowing down", "recovery rate"],
        "classical_validation": ["dose-response single-cell assay"],
        "key_interaction_test": "threshold_shift_effect = Native_transition_dose50_or_slope - Scrambled_or_deletion_transition_dose50_or_slope",
        "temporal": True,
    },
    "negative_space_element_candidate": {
        "assay": "Negative-Space Rescue/Scramble Assay",
        "sequence_controls": ["native", "forbidden-word-inserted", "depleted-kmer-rescued", "GC-matched shuffled", "repeat-inserted"],
        "perturbations": ["baseline", "mild stress", "recovery"],
        "readouts": ["structural change", "reporter variance", "nucleosome positioning", "local chromatin accessibility", "sequence grammar disruption"],
        "classical_validation": ["tiling reporter with inserted/depleted tokens"],
        "key_interaction_test": "negative_space_rescue_effect = (Native_stress_readout - Native_baseline_readout) - (DepletedKmerRescued_or_forbiddenWordInserted_stress_readout - Rescued_or_inserted_baseline_readout)",
    },
    "sequence_regime_boundary_candidate": {
        "assay": "Boundary Disruption Assay",
        "sequence_controls": ["native", "boundary-smoothed", "left-half duplicated", "right-half duplicated", "boundary-shuffled", "neutral spacer"],
        "perturbations": ["baseline", "stress pulse", "mechanical/osmotic perturbation"],
        "readouts": ["noise transfer", "nucleosome positioning", "reporter covariance", "chromatin accessibility boundary shift"],
        "classical_validation": ["ATAC/CUT&Tag boundary profiling", "reporter insulation assay"],
        "key_interaction_test": "boundary_effect = (Native_left_right_readout_discontinuity_after_pulse - Native_baseline_discontinuity) - (BoundarySmoothed_after_pulse_discontinuity - BoundarySmoothed_baseline_discontinuity)",
    },
    "TE_grammar_node_candidate": {
        "assay": "TE Grammar Reconstruction Assay",
        "sequence_controls": ["native TE mosaic", "TE-order scrambled", "TE-orientation inverted", "TE-boundary deleted", "family-matched random TE fragments"],
        "perturbations": ["baseline", "stress", "recovery", "cell-state perturbation"],
        "readouts": ["reporter response", "chromatin state", "noise", "folding", "accessibility", "TE-family-specific behavior"],
        "classical_validation": ["TE-overlap enrichment", "TE-derived cCRE comparison", "CRISPRi of TE-derived candidate"],
        "key_interaction_test": "TE_grammar_effect = (NativeTE_mosaic_response - TEOrderScrambled_response) - (FamilyMatchedRandomTE_response - NeutralControl_response)",
    },
    "unexplained_dark_anomaly_candidate": {
        "assay": "Matched-Null Prioritization Assay",
        "sequence_controls": ["native", "GC-matched shuffled", "dinucleotide-preserved shuffled", "nearby neutral genomic control"],
        "perturbations": ["baseline", "mild stress", "recovery"],
        "readouts": ["reporter response", "chromatin accessibility", "sequence-shape proxy", "single-cell variance"],
        "classical_validation": ["matched-null reanalysis", "artifact and mappability review"],
        "key_interaction_test": "dark_anomaly_effect = (Native_stress_or_context_readout - Native_baseline_readout) - (MatchedNullOrGCMatchedControl_stress_readout - MatchedNullOrGCMatchedControl_baseline_readout)",
    },
    "no_call": {
        "assay": "No primitive-specific assay recommended",
        "sequence_controls": ["re-run with stronger evidence", "review matched-null controls", "review artifact flags"],
        "perturbations": ["not recommended until candidate thresholds are met"],
        "readouts": ["residual z-score", "matched-null z-score", "artifact-risk review"],
        "classical_validation": ["inspect classical covariates and matched nulls before assay design"],
        "key_interaction_test": "no_call_review = below candidate threshold; do not prioritize a wet-lab interaction test without stronger residual or matched-null support",
    },
}


def recommend_assay(primitive: str) -> dict:
    if str(primitive) == "no_call":
        blueprint = deepcopy(ASSAY_BLUEPRINTS["no_call"])
        blueprint["candidate_name"] = "no_call"
        blueprint["confirmed_name"] = "no_call"
        blueprint["requires_dynamic_data"] = False
        blueprint["required_validation_data"] = "review_only"
        blueprint["prompt1_allowed_interpretation"] = "No sequence-derived primitive candidate was assigned."
        blueprint["suggested_prompt2_view"] = "none"
        blueprint["mechanistic_bridge"] = deepcopy(MECHANISTIC_BRIDGES["no_call"])
        blueprint.setdefault("allowed_interpretation", "No assay should be prioritized from this card unless future evidence crosses candidate thresholds.")
        blueprint.setdefault("forbidden_interpretation", "Do not interpret no_call regions as primitive candidates.")
        return blueprint
    ontology = get_primitive(primitive)
    candidate_name = ontology.candidate_name
    blueprint = deepcopy(ASSAY_BLUEPRINTS.get(candidate_name, ASSAY_BLUEPRINTS["unexplained_dark_anomaly_candidate"]))
    blueprint["candidate_name"] = ontology.candidate_name
    blueprint["confirmed_name"] = ontology.confirmed_name
    blueprint["requires_dynamic_data"] = ontology.requires_dynamic_data
    blueprint["required_validation_data"] = ontology.required_input_level
    blueprint["prompt1_allowed_interpretation"] = ontology.prompt1_interpretation
    blueprint["suggested_prompt2_view"] = ontology.suggested_prompt2_view
    blueprint["mechanistic_bridge"] = deepcopy(MECHANISTIC_BRIDGES.get(candidate_name, MECHANISTIC_BRIDGES["unexplained_dark_anomaly_candidate"]))
    blueprint.setdefault("key_interaction_test", KEY_INTERACTION)
    if blueprint.get("temporal"):
        blueprint["temporal_interaction_test"] = TEMPORAL_INTERACTION
    blueprint.setdefault("allowed_interpretation", "A positive result supports a sequence-dependent, perturbation-sensitive hypothesis requiring replication and orthogonal validation.")
    blueprint.setdefault("forbidden_interpretation", "Do not claim confirmed function, quantum effects, teleology, or mechanism from this assay blueprint alone.")
    return blueprint
