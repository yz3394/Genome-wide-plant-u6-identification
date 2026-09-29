"""Strand-aware extraction and descriptive U6 promoter evidence.

Coordinates supplied to ``fetch`` are 1-based, closed genomic intervals. Motif
coordinates are transcript-oriented: the predicted first RNA base is +1 and
the immediately upstream base is -1 (there is no promoter position zero).
No function discovers a TSS, establishes transcription, or excludes a U6 locus.
Call these functions separately for each externally supported boundary hypothesis.
"""

from __future__ import annotations

import re
from typing import Callable


IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT",
    "M": "AC", "K": "GT", "S": "CG", "W": "AT", "B": "CGT",
    "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
USE_THREE_AT = "RTMCCACATCG"
USE_CLASSIC = "RTCCCACATCG"
TATA_REFERENCE = "TTTATATA"
MSP_CONSENSUS = "RGCCCR"  # Connelly et al. 1994, PMID 8065324.


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans(
        "ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN"
    ))[::-1]


def extract_upstream(
    fetch: Callable[[str, int, int], str], seqid: str, strand: str,
    tss1: int, length: int, seq_length: int,
) -> dict:
    """Extract the available upstream window, never pad an assembly edge.

    An out-of-reference boundary returns an explicit per-locus failure record;
    an ordinary truncated promoter is returned with its actual length. Malformed
    arguments or a failed FASTA fetch are errors, not missing biological evidence.
    """
    if strand not in {"+", "-"} or length <= 0 or seq_length <= 0:
        raise ValueError("Expected strand +/-, positive length and sequence length")
    result = {
        "seqid": seqid, "strand": strand, "tss1": tss1,
        "requested_length": length, "sequence": "", "start1": None,
        "end1": None, "actual_length": 0, "truncated": True, "N_count": 0,
        "ambiguous_base_count": 0, "first_relative_position": None,
        "boundary_interpretation": "externally_supplied_hypothesis_not_measured_TSS",
    }
    if not 1 <= tss1 <= seq_length:
        return {**result, "extraction_status": "boundary_outside_reference"}
    if strand == "+":
        start, end = max(1, tss1 - length), tss1 - 1
    else:
        start, end = tss1 + 1, min(seq_length, tss1 + length)
    if end < start:
        return {**result, "extraction_status": "no_upstream_sequence_at_reference_edge"}
    sequence = fetch(seqid, start, end).upper()
    if len(sequence) != end - start + 1:
        raise ValueError(f"FASTA fetch length mismatch: {seqid}:{start}-{end}")
    if strand == "-":
        sequence = reverse_complement(sequence)
    return {
        **result, "sequence": sequence, "start1": start, "end1": end,
        "actual_length": len(sequence), "truncated": len(sequence) < length,
        "N_count": sequence.count("N"),
        "ambiguous_base_count": sum(base not in "ACGT" for base in sequence),
        "first_relative_position": -len(sequence),
        "extraction_status": "truncated_reference_edge" if len(sequence) < length else "complete",
    }


def motif_mismatches(sequence: str, consensus: str) -> int:
    if len(sequence) != len(consensus):
        raise ValueError("Motif and consensus lengths differ")
    # An ambiguous assembly base is not a supported match to an IUPAC motif.
    return sum(base not in IUPAC[code] for base, code in zip(sequence.upper(), consensus))


def _pair_record(use: dict, tata: dict, pair_class: str) -> dict:
    edge = tata["start"] - use["end"]
    return {
        "use_sequence": use["sequence"], "use_start_rel_tss": use["start"],
        "use_end_rel_tss": use["end"],
        "use_mismatch_three_at_consensus": use["mm"],
        "use_mismatch_classic_consensus": use["classic_mm"],
        "tata_sequence": tata["sequence"], "tata_start_rel_tss": tata["start"],
        "tata_end_rel_tss": tata["end"], "tata_mismatch_at_reference": tata["mm"],
        "use_tata_start_separation": tata["start"] - use["start"],
        "use_tata_edge_distance": edge, "use_tata_intervening_bases": edge - 1,
        "use_tata_center_distance": (tata["start"] + tata["end"] - use["start"] - use["end"]) / 2,
        "motif_pair_class": pair_class,
        "motif_ambiguous_base_count": sum(base not in "ACGT" for base in use["sequence"] + tata["sequence"]),
    }


def scan_motifs(
    sequence: str, first_relative_position: int,
    profile: str = "legacy_arabidopsis", *, msp: bool = False,
) -> dict:
    """Report paired sequence/geometry evidence, not promoter activity.

    ``legacy_arabidopsis`` freezes the original two-species scoring and tie order.
    Its P3 'best' can be a high-mismatch algorithmic placeholder and is flagged.
    ``generic_provisional`` expands positional windows for exploratory review;
    it is NOT calibrated across plants and must not be used as a rejection rule.
    Both profiles currently use the same Arabidopsis-derived motifs. Multiple
    pairs with <=3 mismatches per motif are retained, with four spacing metrics.
    Optional exact RGCCCR MSP matches are descriptive only and never affect rank.
    """
    if profile not in {"legacy_arabidopsis", "generic_provisional"}:
        raise ValueError(f"Unknown motif profile: {profile}")
    sequence = sequence.upper()
    legacy = profile == "legacy_arabidopsis"
    use_window = (-90, -45) if legacy else (-120, -40)
    tata_window = (-42, -18) if legacy else (-55, -15)
    last_relative = first_relative_position + len(sequence) - 1
    complete = first_relative_position <= use_window[0] and last_relative >= tata_window[1] + len(TATA_REFERENCE) - 1
    use_hits, tata_hits = [], []
    for consensus, window, target in (
        (USE_THREE_AT, use_window, use_hits), (TATA_REFERENCE, tata_window, tata_hits),
    ):
        for index in range(len(sequence) - len(consensus) + 1):
            start = first_relative_position + index
            if not window[0] <= start <= window[1]:
                continue
            subsequence = sequence[index:index + len(consensus)]
            target.append({
                "sequence": subsequence, "start": start, "end": start + len(consensus) - 1,
                "mm": motif_mismatches(subsequence, consensus),
                "classic_mm": motif_mismatches(subsequence, USE_CLASSIC) if consensus == USE_THREE_AT else None,
            })
    ranked, reported = [], []
    for use in use_hits:
        for tata in tata_hits:
            if tata["start"] <= use["end"]:
                continue
            edge = tata["start"] - use["end"]
            strict = use["mm"] <= 1 and tata["mm"] <= 2 and -35 <= tata["start"] <= -25 and 24 <= edge <= 26
            relaxed = use["mm"] <= 2 and tata["mm"] <= 2 and 21 <= edge <= 29
            class_rank = 0 if strict else 1 if relaxed else 2
            pair_class = ("P1_canonical_geometry", "P2_near_canonical_geometry", "P3_uncertain")[class_rank]
            record = _pair_record(use, tata, pair_class)
            if not legacy:
                plausible = relaxed and record["motif_ambiguous_base_count"] == 0
                record["motif_pair_class"] = "provisional_type3_like_geometry" if plausible else "provisional_uncertain_geometry"
            key = (
                class_rank, use["mm"] + tata["mm"], abs(edge - 24.5),
                abs(use["start"] + 65.5) + abs(tata["start"] + 31),
            )
            if use["mm"] <= 3 and tata["mm"] <= 3:
                reported.append(record)
            # Keep P3 placeholders only for the historical reproducibility profile.
            if legacy or (use["mm"] <= 3 and tata["mm"] <= 3):
                ranked.append((key, record))
    best = dict(min(ranked, key=lambda item: item[0])[1]) if ranked else {}
    if best:
        best["motif_pair_selection_basis"] = (
            "canonical_geometry_then_minimum_mismatch"
            if legacy and best["motif_pair_class"] != "P3_uncertain"
            else "algorithm_minimum_mismatch_P3_not_a_functional_best_pair" if legacy
            else "provisional_Arabidopsis_motif_heuristic_not_cross_plant_calibrated"
        )
    reliable_use_pair = bool(best) and best["motif_ambiguous_base_count"] == 0 and (
        best["motif_pair_class"] in {"P1_canonical_geometry", "P2_near_canonical_geometry", "provisional_type3_like_geometry"}
    )
    msp_hits = []
    if msp:
        for index in range(len(sequence) - len(MSP_CONSENSUS) + 1):
            subsequence = sequence[index:index + len(MSP_CONSENSUS)]
            if motif_mismatches(subsequence, MSP_CONSENSUS):
                continue
            start = first_relative_position + index
            end = start + len(MSP_CONSENSUS) - 1
            if reliable_use_pair and end >= best["use_start_rel_tss"]:
                continue
            msp_hits.append({
                "sequence": subsequence, "start_rel_tss": start, "end_rel_tss": end,
                "consensus": MSP_CONSENSUS,
                "relation_to_USE": "upstream_of_selected_USE" if reliable_use_pair else "USE_relation_unresolved",
                "interpretation": "descriptive_exact_match_not_function_or_lineage_evidence",
            })
    return {
        "best": best, "pairs": reported, "msp_hits": msp_hits,
        "scan_status": "complete_window" if complete else "partial_window" if use_hits and tata_hits else "insufficient_window",
        "motif_evidence_status": "paired_heuristic_evidence" if reliable_use_pair else "uncertain_or_no_pair_not_U6_exclusion",
        "profile": profile,
        "calibration_status": "historical_Arabidopsis_heuristic_not_activity_calibrated" if legacy else "provisional_not_cross_plant_calibrated",
        "msp_scan_status": "descriptive_exact_consensus_only" if msp else "not_requested",
        "msp_search_range": "supplied_promoter_upstream_of_selected_USE" if reliable_use_pair else "supplied_promoter_USE_relation_unresolved",
        "use_start_scan_window": list(use_window), "tata_start_scan_window": list(tata_window),
        "sequence_ambiguous_base_count": sum(base not in "ACGT" for base in sequence),
        "absence_interpretation": "not_evidence_of_absent_promoter_or_pseudogene",
    }


def scan_termination(oriented_context: str, body_span: int) -> dict:
    """Describe T4+ tracts in context beginning at the candidate +1 base.

    Retains the old expected-start window (body_span-15 .. body_span+30).
    A body-end estimate is not an experimentally defined transcript 3' end.
    A hit elsewhere is NEVER returned as the chosen expected termination hit.
    'internal' means wholly before this terminal-proximal window, not merely
    overlapping the body; terminal oligo-T commonly overlaps a body estimate.
    """
    if body_span <= 0:
        raise ValueError("body_span must be positive")
    context = oriented_context.upper()
    hits = [
        {"start_rel_tss": match.start() + 1, "end_rel_tss": match.end(), "sequence": match.group(),
         "right_censored_by_context_end": match.end() == len(context)}
        for match in re.finditer(r"T{4,}", context)
    ]
    window_start, requested_end = max(1, body_span - 15), body_span + 30
    window_end = min(len(context), requested_end)
    expected = [hit for hit in hits if window_start <= hit["start_rel_tss"] <= window_end]
    internal = [hit for hit in hits if hit["end_rel_tss"] < window_start]
    chosen = min(expected, key=lambda hit: abs(hit["start_rel_tss"] - body_span + 3)) if expected else None
    fully_observed = len(context) >= requested_end
    ambiguous_count = sum(base not in "ACGT" for base in context[window_start - 1:window_end])
    return {
        "termination_t4plus_expected_window": "yes" if expected else "no" if fully_observed and not ambiguous_count else "not_fully_assessable",
        "termination_run_sequence": chosen["sequence"] if chosen else "",
        "termination_start_rel_tss": chosen["start_rel_tss"] if chosen else "",
        "termination_end_rel_tss": chosen["end_rel_tss"] if chosen else "",
        "termination_all_t4plus": ",".join(f"{hit['start_rel_tss']}:{hit['end_rel_tss']}:{hit['sequence']}" for hit in hits),
        "termination_scan_window_rel_tss": f"{window_start}..{window_end}" if window_end >= window_start else "unavailable",
        "termination_requested_window_rel_tss": f"{window_start}..{requested_end}",
        "termination_scan_status": "complete_window" if fully_observed else "partial_window" if window_end >= window_start else "unavailable_window",
        "termination_window_ambiguous_base_count": ambiguous_count,
        "termination_expected_hits": expected, "internal_t4plus_hits": internal,
        "body_overlapping_t4plus_hits": [hit for hit in hits if hit["start_rel_tss"] <= body_span],
        "terminal_window_overlapping_t4plus_hits": [hit for hit in hits if hit["end_rel_tss"] >= window_start and hit["start_rel_tss"] <= window_end],
        "all_t4plus_hits": hits,
        "termination_interpretation": "sequence_feature_not_measured_termination_or_expression",
    }

