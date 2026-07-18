"""Structure-specific Level-1 non-B-DNA sequence-potential records.

Predicted motif potential, context-conditioned formation, and observed
structure are separate evidence levels.  This module has sequence alone and
therefore emits Level 1 only; it never averages unlike structures into a
mechanistic score.
"""

from __future__ import annotations

import json
import math
import re
from collections import Counter
from dataclasses import asdict, dataclass

import numpy as np

from .repeats import g_tract_density
from .sequence import clean_sequence


NONB_PREDICTOR_VERSION = "internal_motif_screen_v2"


@dataclass(frozen=True)
class NonBStructureRecord:
    region_id: str | None
    structure_type: str
    strand: str
    start_offset: int
    end_offset: int
    predictor: str
    predictor_version: str
    raw_score: float
    calibrated_score: float | None
    formation_context_status: str
    experimental_support: str
    prediction_agreement: str
    evidence_level: str = "level_1_sequence_potential"


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def longest_alternating(seq: str, groups: tuple[set[str], set[str]]) -> int:
    best = 0
    current = 0
    last_group = None
    for base in seq:
        group = 0 if base in groups[0] else 1 if base in groups[1] else None
        if group is None:
            current = 0
            last_group = None
        elif last_group is None or group != last_group:
            current += 1
            last_group = group
        else:
            current = 1
            last_group = group
        best = max(best, current)
    return best


def repeat_density(seq: str, mode: str = "direct", k: int = 6) -> float:
    seq = clean_sequence(seq)
    if len(seq) < 2 * k:
        return 0.0
    complement = str.maketrans("ACGT", "TGCA")
    counts = Counter(seq[index : index + k] for index in range(len(seq) - k + 1) if "N" not in seq[index : index + k])
    total = sum(counts.values())
    if total == 0:
        return 0.0
    hits = 0
    for token, count in counts.items():
        if mode == "direct" and count > 1:
            hits += count
        elif mode == "inverted" and counts.get(token.translate(complement)[::-1], 0) > 0:
            hits += count
        elif mode == "mirror" and counts.get(token[::-1], 0) > 0:
            hits += count
    return hits / total


def _record(
    structure_type: str,
    strand: str,
    start: int,
    end: int,
    raw_score: float,
    *,
    region_id: str | None,
    predictor: str,
) -> NonBStructureRecord:
    return NonBStructureRecord(
        region_id=region_id,
        structure_type=structure_type,
        strand=strand,
        start_offset=int(start),
        end_offset=int(end),
        predictor=predictor,
        predictor_version=NONB_PREDICTOR_VERSION,
        raw_score=float(raw_score),
        calibrated_score=None,
        formation_context_status="unavailable_requires_context_tracks_or_physical_conditions",
        experimental_support="unavailable_no_observed_structure_input",
        prediction_agreement="single_internal_screen_not_cross_predictor_validated",
    )


def _regex_records(
    seq: str,
    structure_type: str,
    pattern: str,
    *,
    region_id: str | None,
    predictor: str,
    both_strands: bool = True,
) -> list[NonBStructureRecord]:
    records: list[NonBStructureRecord] = []
    orientations = [("+", seq)]
    if both_strands:
        orientations.append(("-", reverse_complement(seq)))
    for strand, oriented in orientations:
        for match in re.finditer(pattern, oriented):
            if strand == "+":
                start, end = match.start(), match.end()
            else:
                start, end = len(seq) - match.end(), len(seq) - match.start()
            records.append(_record(structure_type, strand, start, end, match.end() - match.start(), region_id=region_id, predictor=predictor))
    return records


def _inverted_repeat_records(seq: str, *, region_id: str | None, maximum_records: int = 64) -> list[NonBStructureRecord]:
    records: list[NonBStructureRecord] = []
    for arm_length in range(6, 11):
        for start in range(0, max(0, len(seq) - 2 * arm_length - 3)):
            left = seq[start : start + arm_length]
            if "N" in left:
                continue
            reverse_arm = reverse_complement(left)
            minimum_right = start + arm_length + 3
            maximum_right = min(len(seq) - arm_length, start + arm_length + 20)
            right_start = seq.find(reverse_arm, minimum_right, maximum_right + arm_length)
            if right_start >= 0 and right_start <= maximum_right:
                right_end = right_start + arm_length
                loop_length = right_start - (start + arm_length)
                score = arm_length / max(1.0, loop_length)
                records.append(_record("hairpin", "+", start, right_end, score, region_id=region_id, predictor="inverted_repeat_hairpin_screen"))
                records.append(_record("cruciform", ".", start, right_end, score, region_id=region_id, predictor="palindromic_cruciform_screen"))
            if len(records) >= maximum_records:
                return records
    return records


def _rloop_records(seq: str, *, region_id: str | None, window: int = 30) -> list[NonBStructureRecord]:
    records: list[NonBStructureRecord] = []
    if len(seq) < window:
        return records
    for start in range(0, len(seq) - window + 1, max(1, window // 3)):
        segment = seq[start : start + window]
        g = segment.count("G")
        c = segment.count("C")
        g_fraction = g / window
        skew = (g - c) / max(1, g + c)
        if g_fraction >= 0.5 and skew >= 0.25:
            records.append(_record("R_loop_susceptibility", "+", start, start + window, g_fraction + skew, region_id=region_id, predictor="G_rich_non_template_strand_screen"))
    return records


def compute_nonb_structure_records(seq: str, region_id: str | None = None) -> list[dict[str, object]]:
    seq = clean_sequence(seq)
    records: list[NonBStructureRecord] = []
    records += _regex_records(
        seq,
        "G4",
        r"G{3,}[ACGTN]{1,7}G{3,}[ACGTN]{1,7}G{3,}[ACGTN]{1,7}G{3,}",
        region_id=region_id,
        predictor="four_G_tract_motif_screen",
    )
    records += _regex_records(
        seq,
        "i_motif",
        r"C{3,}[ACGTN]{1,7}C{3,}[ACGTN]{1,7}C{3,}[ACGTN]{1,7}C{3,}",
        region_id=region_id,
        predictor="four_C_tract_motif_screen",
    )
    records += _regex_records(seq, "Z_DNA", r"(?:(?:CG)|(?:GC)){4,}", region_id=region_id, predictor="alternating_GC_screen", both_strands=False)
    records += _regex_records(seq, "triplex_H_DNA", r"[AG]{12,}|[CT]{12,}", region_id=region_id, predictor="homopurine_homopyrimidine_screen")
    records += _regex_records(seq, "slipped_DNA", r"([ACGT]{2,6})\1{2,}", region_id=region_id, predictor="short_tandem_repeat_slippage_screen", both_strands=False)
    records += _regex_records(seq, "A_tract_curvature", r"A{4,}|T{4,}", region_id=region_id, predictor="A_tract_sequence_screen", both_strands=False)
    records += _inverted_repeat_records(seq, region_id=region_id)
    records += _rloop_records(seq, region_id=region_id)
    unique: dict[tuple[str, str, int, int], NonBStructureRecord] = {}
    for record in records:
        key = (record.structure_type, record.strand, record.start_offset, record.end_offset)
        unique[key] = record
    return [asdict(record) for record in sorted(unique.values(), key=lambda item: (item.start_offset, item.end_offset, item.structure_type, item.strand))]


def _maximum_record_score(records: list[dict[str, object]], structure_type: str) -> float:
    values = [float(record["raw_score"]) for record in records if record["structure_type"] == structure_type]
    return max(values) if values else 0.0


def compute_nonb_dna_features(seq: str) -> dict[str, float | int | str]:
    seq = clean_sequence(seq)
    length = max(1, len(seq))
    records = compute_nonb_structure_records(seq)
    types = sorted({str(record["structure_type"]) for record in records})
    g = seq.count("G")
    c = seq.count("C")
    g_skew = (g - c) / max(1, g + c)
    g4_proxy = g_tract_density(seq, min_run=3) * (1.0 + max(0.0, g_skew))
    z_proxy = longest_alternating(seq, ({"G", "C"}, {"A", "T"})) / length
    purine_pyrimidine = longest_alternating(seq, ({"A", "G"}, {"C", "T"})) / length
    phased_a = sum(1 for index in range(0, max(0, len(seq) - 4), 10) if seq[index : index + 4].count("A") >= 3) / max(1.0, length / 10)
    inverted = repeat_density(seq, "inverted")
    direct = repeat_density(seq, "direct")
    mirror = repeat_density(seq, "mirror")
    triplex = (purine_pyrimidine + max(0.0, g_skew) + seq.count("AGG") / length) / 3.0
    rloop = _maximum_record_score(records, "R_loop_susceptibility")
    deprecated_warning = (
        "non_B_DNA_aggregate_score has been retired because unlike structures cannot be averaged into one mechanistic score. "
        "Use nonb_structure_records_json and structure-specific fields."
    )
    return {
        "nonb_structure_records_json": json.dumps(records, sort_keys=True),
        "nonb_structure_record_count": int(len(records)),
        "nonb_structure_types_detected": ";".join(types),
        "nonb_evidence_level": "level_1_sequence_potential",
        "nonb_context_conditioned_formation_status": "unavailable_requires_context_tracks_or_physical_conditions",
        "nonb_observed_structure_status": "unavailable_no_experimental_structure_input",
        "G4_sequence_potential": float(g4_proxy),
        "i_motif_sequence_potential": float(_maximum_record_score(records, "i_motif") / length),
        "Z_DNA_sequence_potential": float(z_proxy),
        "triplex_H_DNA_sequence_potential": float(triplex),
        "cruciform_sequence_potential": float(_maximum_record_score(records, "cruciform")),
        "hairpin_sequence_potential": float(_maximum_record_score(records, "hairpin")),
        "slipped_DNA_sequence_potential": float(_maximum_record_score(records, "slipped_DNA") / length),
        "R_loop_susceptibility_sequence_potential": float(rloop),
        "A_tract_curvature_sequence_potential": float(phased_a),
        # Compatibility proxy fields remain Level-1 screens.
        "G4_susceptibility_proxy": float(g4_proxy),
        "Z_DNA_propensity_proxy": float(z_proxy),
        "A_phased_tract_score": float(phased_a),
        "inverted_repeat_density": float(inverted),
        "direct_repeat_density": float(direct),
        "mirror_repeat_density": float(mirror),
        "triplex_H_DNA_proxy": float(triplex),
        "cruciform_forming_potential": float(_maximum_record_score(records, "cruciform")),
        "R_loop_forming_potential": float(rloop),
        "non_B_DNA_aggregate_score": math.nan,
        "non_B_DNA_aggregate_score_status": "deprecated_unavailable",
        "non_B_DNA_aggregate_score_deprecation_warning": deprecated_warning,
    }
