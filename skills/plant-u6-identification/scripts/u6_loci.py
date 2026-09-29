#!/usr/bin/env python3
"""Evidence-bounded U6 locus curation; no promoter or activity inference.

The defaults reproduce the historical search profile, not universal plant
thresholds. Coordinates are one-based inclusive. Every strong BLAST decision
is made on ONE HSP, and every cluster retains an immutable seed anchor.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
import hashlib
import re


VERSION = "2.0.0"


@dataclass(frozen=True)
class CurationPolicy:
    name: str = "historical_80_80_80_single_HSP_fixed_anchor_v2"
    strong_identity: float = 80.0
    strong_coverage: float = 80.0
    strong_bitscore: float = 80.0
    moderate_identity: float = 70.0
    moderate_coverage: float = 60.0
    moderate_bitscore: float = 50.0
    minimum_overlap_nt: int = 15
    minimum_shorter_overlap: float = 0.50
    maximum_envelope_nt: int = 140
    duplicate_reciprocal_overlap: float = 0.90
    competition_reciprocal_overlap: float = 0.50
    historical_query_token: str = "X52528.1"


@dataclass
class Hit:
    raw_hit_id: str
    method: str
    query: str
    seqid: str
    start1: int
    end1: int
    strand: str
    source: str
    line: int
    pident: float | None = None
    qstart: int | None = None
    qend: int | None = None
    qlen: int | None = None
    qcov_reported: float | None = None
    qcov_exact: float | None = None
    bitscore: float | None = None
    gapopen: int | None = None
    alignment_length: int | None = None
    qseq: str | None = None
    sseq: str | None = None
    btop: str | None = None
    subject_length: int | None = None
    score: float | None = None
    evalue: float | None = None
    model_from: int | None = None
    model_to: int | None = None
    model_coverage: float | None = None
    trunc: str | None = None
    bias: float | None = None
    ga: float | None = None
    ga_pass: bool = False
    suppressed_by: str = ""
    locus_key: str = ""
    multi_anchor_conflict: bool = False

    @property
    def length(self):
        return self.end1 - self.start1 + 1


def read_cm_metadata(path):
    """Read the first CM header, never the appended HMM header or hardcoded GA."""
    wanted = {"ACC", "NAME", "CLEN", "GA", "W"}
    fields = {}
    with Path(path).open() as handle:
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "CM":
                break
            if parts[0] in wanted and parts[0] not in fields:
                fields[parts[0]] = parts[1].rstrip(";")
    missing = {"ACC", "CLEN", "GA"} - fields.keys()
    if missing:
        raise ValueError(f"CM header lacks {sorted(missing)}: {path}")
    return {"accession": fields["ACC"], "name": fields.get("NAME", ""),
            "clen": int(fields["CLEN"]), "ga": float(fields["GA"]),
            "window": int(fields["W"]) if "W" in fields else None,
            "source": str(path)}


def parse_blast(species, path):
    hits = []
    with Path(path).open() as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) not in {15, 18}:
                raise ValueError(f"Expected frozen 15-column or extended 18-column BLAST format: {path}:{line_no}")
            qs, qe, ss, se, qlen, slen = map(int, (f[6], f[7], f[8], f[9], f[12], f[13]))
            if not (1 <= qs <= qe <= qlen and 1 <= min(ss, se) <= max(ss, se) <= slen):
                raise ValueError(f"Invalid nucleotide coordinates: {path}:{line_no}")
            if len(f) == 18 and (len(f[15]) != int(f[3]) or len(f[16]) != int(f[3])
                                or len(f[15].replace("-", "")) != qe-qs+1
                                or len(f[16].replace("-", "")) != abs(se-ss)+1):
                raise ValueError(f"Aligned BLAST sequences disagree with coordinates: {path}:{line_no}")
            hits.append(Hit(
                f"{species}_BLAST_{len(hits)+1:06d}", "blastn", f[0], f[1],
                min(ss, se), max(ss, se), "+" if ss <= se else "-", str(path), line_no,
                pident=float(f[2]), qstart=qs, qend=qe, qlen=qlen,
                qcov_reported=float(f[14]), qcov_exact=100.0*(qe-qs+1)/qlen,
                bitscore=float(f[11]), gapopen=int(f[5]), alignment_length=int(f[3]),
                qseq=f[15] if len(f) == 18 else None,
                sseq=f[16] if len(f) == 18 else None,
                btop=f[17] if len(f) == 18 else None,
                subject_length=slen, evalue=float(f[10])))
    return hits


def parse_cm(species, path, model):
    hits = []
    with Path(path).open() as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            f = raw.split(maxsplit=17)
            if len(f) < 17 or f[4] != "cm":
                raise ValueError(f"Expected cmsearch CM tblout: {path}:{line_no}")
            if f[3] not in {model["accession"], "-"}:
                raise ValueError(f"Model/table accession mismatch: {path}:{line_no}")
            mf, mt, sf, st = map(int, f[5:9])
            if not (1 <= mf <= mt <= model["clen"] and min(sf, st) >= 1 and f[9] in "+-"):
                raise ValueError(f"Invalid CM coordinates: {path}:{line_no}")
            score = float(f[14])
            hits.append(Hit(
                f"{species}_{model['accession']}_{len(hits)+1:06d}", "cmsearch",
                model["accession"], f[0], min(sf, st), max(sf, st), f[9], str(path), line_no,
                score=score, evalue=float(f[15]), model_from=mf, model_to=mt,
                model_coverage=(mt-mf+1)/model["clen"], trunc=f[10], bias=float(f[13]),
                ga=model["ga"], ga_pass=score >= model["ga"]))
    return hits


def overlap_nt(a, b):
    if (a.seqid, a.strand) != (b.seqid, b.strand):
        return 0
    return max(0, min(a.end1, b.end1)-max(a.start1, b.start1)+1)


def single_hsp_support(hit, policy=CurationPolicy()):
    if hit.method != "blastn":
        return "none"
    if (hit.pident >= policy.strong_identity and hit.qcov_exact >= policy.strong_coverage
            and hit.bitscore >= policy.strong_bitscore):
        return "strong"
    if (hit.pident >= policy.moderate_identity and hit.qcov_exact >= policy.moderate_coverage
            and hit.bitscore >= policy.moderate_bitscore):
        return "moderate"
    return "weak"


def hit_rank(hit, policy=CurationPolicy()):
    support = {"none": 0, "weak": 1, "moderate": 2, "strong": 3}
    return (int(hit.ga_pass), support[single_hsp_support(hit, policy)],
            hit.score if hit.score is not None else hit.bitscore or 0,
            hit.qcov_exact or 0, hit.pident or 0, -hit.start1)


def suppress_duplicates(hits, policy=CurationPolicy()):
    """Suppress overlapping redundant calls, retaining all raw evidence."""
    groups = defaultdict(list)
    for hit in hits:
        groups[(hit.seqid, hit.strand, hit.method, hit.query)].append(hit)
    for group in groups.values():
        kept = []
        for hit in sorted(group, key=lambda h: hit_rank(h, policy), reverse=True):
            duplicate = next((other for other in kept if
                              overlap_nt(hit, other)/max(hit.length, other.length)
                              >= policy.duplicate_reciprocal_overlap), None)
            if duplicate:
                hit.suppressed_by = duplicate.raw_hit_id
            else:
                kept.append(hit)


def anchor_compatible(hit, anchor, members, policy=CurationPolicy()):
    overlap = overlap_nt(hit, anchor)
    if (overlap < policy.minimum_overlap_nt or
            overlap/min(hit.length, anchor.length) < policy.minimum_shorter_overlap):
        return False
    envelope = max([hit.end1] + [m.end1 for m in members]) - min(
        [hit.start1] + [m.start1 for m in members]) + 1
    # A long CM anchor is not deleted because of this clustering safeguard.
    return envelope <= max(policy.maximum_envelope_nt, anchor.length)


def cluster_hits(canonical_hits, atac_hits, policy=CurationPolicy()):
    active = [h for h in canonical_hits if not h.suppressed_by]
    anchors = sorted((h for h in active if h.query == "RF00026" and h.ga_pass),
                     key=lambda h: (h.seqid, h.strand, h.start1, h.end1))
    clusters = [{"anchor": h, "hits": [h], "atac": []} for h in anchors]
    anchor_ids = {h.raw_hit_id for h in anchors}
    others = sorted((h for h in active if h.raw_hit_id not in anchor_ids),
                    key=lambda h: hit_rank(h, policy), reverse=True)
    for hit in others:
        compatible = [c for c in clusters if anchor_compatible(hit, c["anchor"], c["hits"], policy)]
        if compatible:
            hit.multi_anchor_conflict = len(compatible) > 1
            chosen = max(compatible, key=lambda c: (overlap_nt(hit, c["anchor"]),
                                                   hit_rank(c["anchor"], policy)))
            chosen["hits"].append(hit)
        else:
            clusters.append({"anchor": hit, "hits": [hit], "atac": []})
    standalone = []
    for hit in atac_hits:
        if hit.suppressed_by:
            continue
        compatible = [c for c in clusters if
                      overlap_nt(hit, c["anchor"]) >= policy.minimum_overlap_nt and
                      overlap_nt(hit, c["anchor"])/max(hit.length, c["anchor"].length)
                      >= policy.competition_reciprocal_overlap]
        if compatible:
            hit.multi_anchor_conflict = len(compatible) > 1
            chosen = max(compatible, key=lambda c: overlap_nt(hit, c["anchor"]))
            chosen["atac"].append(hit)
        else:
            standalone.append(hit)
    return clusters, standalone


def summarize_cluster(species, cluster, policy=CurationPolicy()):
    hits, atac_hits, anchor = cluster["hits"], cluster["atac"], cluster["anchor"]
    blast = [h for h in hits if h.method == "blastn"]
    query_best = {}
    for h in blast:
        if h.query not in query_best or hit_rank(h, policy) > hit_rank(query_best[h.query], policy):
            query_best[h.query] = h
    strong = [q for q, h in query_best.items() if single_hsp_support(h, policy) == "strong"]
    rf = max((h for h in hits if h.query == "RF00026"), key=lambda h: h.score, default=None)
    atac = max(atac_hits, key=lambda h: h.score, default=None)
    rf_pass, atac_pass = bool(rf and rf.ga_pass), bool(atac and atac.ga_pass)
    if rf_pass and atac_pass:
        family, grade = "U6_U6atac_ambiguous", "B2_clan_ambiguous"
    elif rf_pass:
        family = "RF00026_supported"
        grade = "A1_concordant_U6" if len(strong) >= 2 else "A2_CM_led_U6"
    elif atac_pass:
        family, grade = "RF00619_supported", "U6atac_like"
    elif strong:
        family = "BLAST_supported_review"
        grade = ("B1_subGA_supported" if rf else "C1_BLAST_only") if len(strong) >= 2 else "fragment_or_weak"
    elif rf:
        family = "below_GA_review"
        grade = "C2_lowCM_only" if rf.model_coverage >= .70 and rf.trunc == "no" else "fragment_or_weak"
    else:
        family, grade = "weak_evidence_review", "fragment_or_weak"
    best = max(blast, key=lambda h: hit_rank(h, policy), default=None)
    # Coordinate choice is not a completeness decision. Preserve legacy body
    # selection where possible, but explicitly remove its old completeness claim.
    body = (rf if rf and rf.model_coverage >= .80 and rf.trunc == "no" else
            best if best and best.qcov_exact >= 90 else rf or best)
    historical = max((h for h in blast if policy.historical_query_token in h.query),
                     key=lambda h: hit_rank(h, policy), default=None)
    q1 = None
    anchor_status = "not_available"
    missing = None
    if historical:
        missing = historical.qstart - 1
        if missing == 0:
            q1 = historical.start1 if historical.strand == "+" else historical.end1
            anchor_status = "direct_query_position_1"
        elif historical.gapopen == 0:
            proposed = historical.start1 - missing if historical.strand == "+" else historical.end1 + missing
            if 1 <= proposed <= historical.subject_length:
                q1, anchor_status = proposed, "provisional_linear_5prime_extrapolation"
            else:
                anchor_status = "extrapolation_outside_reference"
        else:
            anchor_status = "gapped_HSP_requires_alignment_review"
    key = f"{species}:{anchor.seqid}:{anchor.strand}:{anchor.start1}-{anchor.end1}"
    for h in hits + atac_hits:
        h.locus_key = key
    return {
        "locus_key": key, "species_key": species, "seqid": body.seqid, "strand": body.strand,
        "anchor_start1": anchor.start1, "anchor_end1": anchor.end1,
        "anchor_source": anchor.query, "body_start1": body.start1, "body_end1": body.end1,
        "body_length": body.length, "body_boundary_source": body.query,
        "evidence_envelope_start1": min(h.start1 for h in hits),
        "evidence_envelope_end1": max(h.end1 for h in hits),
        "family_class": family, "legacy_grade": grade,
        "canonical_retained": rf_pass and not atac_pass,
        "strict_historical_query_pass": bool(historical and single_hsp_support(historical, policy) == "strong" and rf_pass and not atac_pass),
        "best_query": best.query if best else None,
        "best_query_5p_anchor": (best.start1 if best.strand == "+" else best.end1) if best and best.qstart == 1 else None,
        "best_pident": best.pident if best else None,
        "best_qcov_exact": best.qcov_exact if best else None,
        "best_bitscore": best.bitscore if best else None,
        "n_blast_queries": len(query_best), "n_strong_queries": len(strong),
        "strong_queries": sorted(strong), "best_q1_anchor": q1, "q1_missing_nt": missing,
        "q1_anchor_status": anchor_status, "q1_anchor_safe": missing == 0 if historical else False,
        "historical_hsp_id": historical.raw_hit_id if historical else None,
        "rf00026_score": rf.score if rf else None, "rf00026_ga_pass": rf_pass,
        "rf00026_evalue": rf.evalue if rf else None,
        "rf00026_model_coverage": rf.model_coverage if rf else None,
        "rf00026_trunc": rf.trunc if rf else None,
        "rf00619_score": atac.score if atac else None, "rf00619_ga_pass": atac_pass,
        "rf00619_evalue": atac.evalue if atac else None,
        "multi_anchor_conflict": any(h.multi_anchor_conflict for h in hits+atac_hits),
        "raw_hit_ids": [h.raw_hit_id for h in hits+atac_hits],
        "body_integrity": "unresolved", "integrity_reasons": ["structural_alignment_not_yet_examined"],
    }


def curate_species(species_key, blast_path, u6_tblout, atac_tblout, u6_model, atac_model,
                   *, policy=CurationPolicy()):
    models = [read_cm_metadata(u6_model), read_cm_metadata(atac_model)]
    if [m["accession"] for m in models] != ["RF00026", "RF00619"]:
        raise ValueError("Expected RF00026 and RF00619 models in that order")
    blast = parse_blast(species_key, blast_path)
    u6 = parse_cm(species_key, u6_tblout, models[0])
    atac = parse_cm(species_key, atac_tblout, models[1])
    all_hits = blast+u6+atac
    suppress_duplicates(all_hits, policy)
    clusters, standalone = cluster_hits(blast+u6, atac, policy)
    loci = [summarize_cluster(species_key, c, policy) for c in clusters]
    loci.sort(key=lambda r: (r["seqid"], r["body_start1"], r["body_end1"], r["strand"]))
    by_id = {h.raw_hit_id: h for h in all_hits}
    for hit in all_hits:
        if hit.suppressed_by:
            hit.locus_key = by_id[hit.suppressed_by].locus_key
    return {
        "loci": loci, "evidence": [asdict(h) for h in all_hits],
        "standalone_atac_hits": [asdict(h) for h in standalone],
        "diagnostics": {
            "version": VERSION, "policy": asdict(policy), "models": models,
            "raw_hit_counts": {"blastn": len(blast), "RF00026": len(u6), "RF00619": len(atac)},
            "locus_count": len(loci), "family_counts": dict(Counter(r["family_class"] for r in loci)),
            "legacy_grade_counts": dict(Counter(r["legacy_grade"] for r in loci)),
            "strict_historical_count": sum(r["strict_historical_query_pass"] for r in loci),
            "suppressed_hit_count": sum(bool(h.suppressed_by) for h in all_hits),
            "multi_anchor_hit_count": sum(h.multi_anchor_conflict for h in all_hits),
            "standalone_atac_GA_count": sum(h.ga_pass for h in standalone),
            "absence_interpretation": "not_reported_at_frozen_search_settings; not proof of absence",
        },
    }


def _wuss_pairs(structure):
    opening = {"<": ">", "(": ")", "[": "]", "{": "}"}
    opening.update({chr(n): chr(n+32) for n in range(ord("A"), ord("Z")+1)})
    closing = {v: k for k, v in opening.items()}
    stacks = defaultdict(list)
    pairs = []
    for pos, char in enumerate(structure):
        if char in opening:
            stacks[char].append(pos)
        elif char in closing:
            key = closing[char]
            if not stacks[key]:
                raise ValueError("Unbalanced Stockholm SS_cons closing symbol")
            pairs.append((stacks[key].pop(), pos))
    if any(stacks.values()):
        raise ValueError("Unbalanced Stockholm SS_cons opening symbol")
    return pairs


def parse_stockholm_metrics(path):
    """Parse interleaved cmalign Stockholm RF/SS_cons, ignoring per-sequence PP.

    Occupancy concerns RF match columns, NOT the model coordinate span.
    Pair compatibility is descriptive canonical/wobble compatibility, not a
    probabilistic folding score or proof of functional RNA structure.
    """
    sequences, annotations = defaultdict(str), defaultdict(str)
    ended = False
    with Path(path).open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line == "//":
                ended = True
                continue
            if ended and not line.startswith("#"):
                raise ValueError("Expected one Stockholm alignment, not concatenated alignments")
            if line.startswith("#=GC "):
                _, tag, value = line.split(maxsplit=2)
                annotations[tag] += value
            elif not line.startswith("#"):
                fields = line.split()
                if len(fields) != 2:
                    raise ValueError("Malformed Stockholm sequence line")
                sequences[fields[0]] += fields[1]
    if not ended:
        raise ValueError("Stockholm alignment lacks its terminating // marker")
    if not sequences or "RF" not in annotations:
        raise ValueError("Stockholm sequences and RF annotation are required")
    rf = annotations["RF"]
    if any(len(s) != len(rf) for s in sequences.values()):
        raise ValueError("Stockholm sequence/RF lengths differ")
    ss = annotations.get("SS_cons", "")
    if ss and len(ss) != len(rf):
        raise ValueError("Stockholm SS_cons/RF lengths differ")
    gaps = set(".-_~")
    match_cols = [i for i, c in enumerate(rf) if c not in gaps]
    if not match_cols:
        raise ValueError("Stockholm RF has no match columns")
    match_set = set(match_cols)
    pairs = [(a, b) for a, b in _wuss_pairs(ss) if a in match_set and b in match_set] if ss else []
    results = {}
    for name, sequence in sequences.items():
        occupied = [sequence[i] not in gaps for i in match_cols]
        if any(occupied):
            first = occupied.index(True)
            last = len(occupied)-1-occupied[::-1].index(True)
            interior = occupied[first:last+1]
        else:
            first, last, interior = len(occupied), -1, []
        longest = max((len(m.group(0)) for m in re.finditer("0+", "".join("1" if x else "0" for x in interior))), default=0)
        pair_observed = pair_evaluable = pair_compatible = 0
        for left, right in pairs:
            a, b = sequence[left].upper().replace("T", "U"), sequence[right].upper().replace("T", "U")
            if a not in gaps and b not in gaps:
                pair_observed += 1
                if a in "ACGU" and b in "ACGU":
                    pair_evaluable += 1
                    pair_compatible += a+b in {"AU", "UA", "GC", "CG", "GU", "UG"}
        ungapped = "".join(c for c in sequence if c not in gaps).upper().replace("U", "T")
        results[name] = {
            "aligned_sequence_length": len(ungapped), "rf_columns": len(match_cols),
            "aligned_sequence_sha256": hashlib.sha256(ungapped.encode()).hexdigest(),
            "occupied_rf_columns": sum(occupied), "rf_occupancy": sum(occupied)/len(occupied),
            "terminal_missing_5p": first, "terminal_missing_3p": len(occupied)-last-1 if any(occupied) else 0,
            "internal_missing_columns": sum(not x for x in interior), "longest_internal_deletion": longest,
            "inserted_nt": sum(sequence[i] not in gaps for i in range(len(rf)) if i not in match_set),
            "ambiguous_nt": sum(c not in "ACGT" for c in ungapped),
            "paired_columns": len(pairs)*2, "paired_positions": len(pairs),
            "paired_observed": pair_observed, "paired_evaluable": pair_evaluable,
            "paired_compatible": pair_compatible,
            "pair_compatibility": pair_compatible/pair_evaluable if pair_evaluable else None,
            "ss_cons_available": bool(ss), "alignment_scope": "submitted_sequence_only",
        }
    return results


def classify_integrity(sequence, metrics, cm_coverage=None, cm_trunc=None,
                       annotated_pseudogene=False):
    """Conservative review flags for submitted bodies; never declare function.

    Current numeric review cutoffs are a frozen pilot profile. They are not a
    universal evolutionary definition of intact plant U6 and require clade
    calibration before a new-species production run.
    """
    seq = sequence.upper().replace("U", "T")
    reasons = []
    required = {"rf_occupancy", "internal_missing_columns", "longest_internal_deletion",
                "terminal_missing_5p", "terminal_missing_3p", "aligned_sequence_length"}
    if not metrics or not required.issubset(metrics):
        status = "unresolved"
        reasons.append("structural_alignment_missing_or_incomplete")
    elif metrics["aligned_sequence_length"] != len(seq):
        status = "unresolved"
        reasons.append("alignment_sequence_length_mismatch")
    elif (metrics.get("aligned_sequence_sha256") is not None and
          metrics["aligned_sequence_sha256"] != hashlib.sha256(seq.encode()).hexdigest()):
        status = "unresolved"
        reasons.append("alignment_sequence_checksum_mismatch")
    elif any(c not in "ACGT" for c in seq):
        status = "unresolved"
        reasons.append("ambiguous_reference_bases")
    else:
        if metrics["internal_missing_columns"] > 5 or metrics["longest_internal_deletion"] > 3:
            reasons.append("substantial_internal_model_column_deletion")
        if metrics["terminal_missing_5p"] > 5 or metrics["terminal_missing_3p"] > 5:
            reasons.append("terminal_model_columns_missing_in_submitted_body")
        if metrics["rf_occupancy"] < .90:
            reasons.append("RF_match_column_occupancy_below_0.90")
        if not 80 <= len(seq) <= 130:
            reasons.append("length_outside_pilot_80_130nt_review_window")
        if (metrics.get("paired_evaluable", 0) >= 5 and
                metrics.get("pair_compatibility", 1) < .70):
            reasons.append("low_canonical_wobble_pair_compatibility")
        if "substantial_internal_model_column_deletion" in reasons:
            status = "structural_deletion_risk"
        elif "terminal_model_columns_missing_in_submitted_body" in reasons:
            status = "terminally_incomplete_submitted_body"
        elif reasons:
            status = "atypical_alignment_requires_review"
        else:
            status = "alignment_near_complete_not_functionally_validated"
    if annotated_pseudogene:
        reasons.append("annotated_pseudogene_requires_review")
    return {"body_integrity": status, "integrity_reasons": reasons,
            "annotated_pseudogene": bool(annotated_pseudogene),
            "cm_coordinate_coverage": cm_coverage, "cm_search_truncation": cm_trunc,
            "functional_status": "not_tested",
            "integrity_scope": "submitted_search_body; flanking_sequence_and_transcript_ends_not_validated"}

