#!/usr/bin/env python3
"""Portable, read-only-input U6 discovery and evidence reporting.

Run with --config CONFIG.json --out NEW_DIRECTORY. Inputs are explicit; no
downloads, dependency installation, expression inference, or old-project data
are required. Existing output directories are never reused or overwritten.
Only a fully hash-verified completed search may be reused with --reuse-search.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys
from urllib.parse import unquote

from u6_io import IndexedFasta, digest, read_fasta, revcomp, write_fasta, write_tsv
from u6_loci import curate_species, parse_stockholm_metrics, classify_integrity
from u6_promoters import extract_upstream, scan_motifs, scan_termination
from u6_context import GffContext
from search_u6 import BLAST_FIELDS, SCORING, build_plan, run_plan, capture_tool_versions, resolve_executable

SCRIPT_VERSION = "1.0.0"
CORE_FILES = ("run_u6.py", "search_u6.py", "u6_io.py", "u6_loci.py", "u6_promoters.py", "u6_context.py")
PATH_KEYS = ("genome", "query_panel", "historical_query", "u6_model", "atac_model", "gff")
LOCUS_FIELDS = ["locus_id", "species_key", "species", "assembly_id", "seqid", "strand", "body_start1", "body_end1", "body_length", "family_class", "canonical_retained", "strict_historical_query_pass"]
PROMOTER_FIELDS = ["locus_id", "species_key", "hypothesis", "is_primary_hypothesis", "seqid", "strand", "tss1", "requested_length", "sequence", "start1", "end1", "actual_length", "truncated", "N_count", "ambiguous_base_count", "first_relative_position", "boundary_interpretation", "extraction_status"]


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def overlap(a, b, c, d):
    return max(0, min(b, d) - max(a, c) + 1)


def map_query_position(qseq, sseq, qstart, sstart, qpos):
    """Map an actually aligned query residue; never extrapolate a missing q1."""
    if len(qseq) != len(sseq):
        raise ValueError("Alignment columns differ")
    q, s = qstart, sstart
    for qb, sb in zip(qseq, sseq):
        if qb != "-" and q == qpos:
            return s if sb != "-" else None
        q += qb != "-"
        s += sb != "-"
    return None


def audit_fasta(path, index_out):
    """Stream an uncompressed DNA FASTA and write a NEW audited FAI.

    No source-side index is trusted or written. Wrapped lines must have uniform
    width except the last line in each record; unsupported layouts fail loudly.
    """
    path, index_out = Path(path), Path(index_out)
    if index_out.exists():
        raise ValueError(f"Index output exists: {index_out}")
    entries, lengths, name, count, first_bases, first_width = [], {}, None, 0, None, None
    previous_size, offset, total_n, sha = None, None, 0, hashlib.sha256()

    def finish():
        if name is None:
            return
        if not count:
            raise ValueError(f"Empty FASTA record: {name}")
        lengths[name] = count
        entries.append((name, count, offset, first_bases, first_width))

    with path.open("rb") as handle:
        while True:
            position = handle.tell()
            line = handle.readline()
            if not line:
                break
            sha.update(line)
            if line.startswith(b">"):
                finish()
                header = line[1:].strip().split()
                if not header:
                    raise ValueError("Empty FASTA header")
                name = header[0].decode("utf-8")
                if name in lengths:
                    raise ValueError(f"Duplicate FASTA ID: {name}")
                count, first_bases, first_width, previous_size, offset = 0, None, None, None, None
                continue
            if name is None:
                raise ValueError(f"Expected uncompressed FASTA beginning with >: {path}")
            sequence = line.rstrip(b"\r\n").upper()
            if not sequence or set(sequence) - set(b"ACGTRYSWMKBDHVN"):
                raise ValueError(f"Blank/non-DNA sequence line in {name}; use ordinary uncompressed DNA FASTA")
            if previous_size is not None and previous_size != (first_bases, first_width):
                raise ValueError(f"Nonuniform FASTA wrapping before final line: {name}")
            if first_bases is None:
                first_bases, first_width, offset = len(sequence), len(line), position
            elif len(sequence) > first_bases:
                raise ValueError(f"FASTA final line longer than declared wrap width: {name}")
            previous_size = (len(sequence), len(line))
            count += len(sequence)
            total_n += sequence.count(b"N")
        finish()
    if not lengths:
        raise ValueError("FASTA has no records")
    with index_out.open("x") as handle:
        for entry in entries:
            handle.write("\t".join(map(str, entry)) + "\n")
    return {"path": str(path.resolve()), "sha256": sha.hexdigest(), "records": len(lengths),
            "bases": sum(lengths.values()), "N_bases": total_n, "lengths": lengths,
            "source_index_used": False, "generated_index": str(index_out.resolve()),
            "scope": "All supplied records included; nuclear identity and assembly choice are user-audited metadata"}


def audit_gff(path, lengths, excluded_seqids=()):
    """Reject mismatched/invalid included coordinates, recording explicit exclusions."""
    if path is None:
        return {"status": "not_provided", "included_features": None, "excluded_features": None,
                "context_interpretation": "unassessed, not conflict-free"}
    excluded = set(excluded_seqids)
    if excluded & lengths.keys():
        raise ValueError("excluded_gff_seqids may not hide annotations on included genome records")
    seen, omitted, kinds, declared_ranges = Counter(), Counter(), Counter(), []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region "):
                parts = line.split()
                seqid = unquote(parts[1])
                if seqid not in lengths and seqid not in excluded:
                    raise ValueError(f"GFF/FASTA seqid mismatch at line {number}: {seqid}")
                if seqid in lengths:
                    start, end = int(parts[2]), int(parts[3])
                    if not 1 <= start <= end <= lengths[seqid]:
                        raise ValueError(f"GFF sequence-region out of bounds: {seqid}")
                    declared_ranges.append({"seqid": seqid, "start1": start, "end1": end,
                                            "covers_whole_record": start == 1 and end == lengths[seqid]})
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GFF line {number}: expected nine columns")
            seqid = unquote(fields[0])
            if seqid in excluded:
                omitted[seqid] += 1
                continue
            if seqid not in lengths:
                raise ValueError(f"GFF/FASTA seqid mismatch at line {number}: {seqid}; explicit exclusions required")
            start, end = int(fields[3]), int(fields[4])
            if not 1 <= start <= end <= lengths[seqid] or fields[6] not in {"+", "-", ".", "?"}:
                raise ValueError(f"GFF coordinate/strand out of bounds at line {number}: {seqid}:{start}-{end}")
            seen[seqid] += 1
            kinds[fields[2]] += 1
    return {"status": "coordinate_checks_passed", "included_features": sum(seen.values()),
            "included_seqids": sorted(seen), "excluded_features": dict(omitted),
            "explicit_excluded_seqids": sorted(excluded), "feature_types": dict(kinds),
            "declared_sequence_regions": declared_ranges,
            "partial_sequence_region_note": "A valid declared subrange is permitted; it does not establish complete annotation coverage",
            "included_genome_records_without_features": sorted(lengths.keys() - seen.keys()),
            "context_interpretation": "Available classes are not guaranteed complete; missing class/seqid is NA"}


def load_config(path):
    source = Path(path).resolve(strict=True)
    c = json.loads(source.read_text())
    for key in ("key", "prefix", "species", "assembly", "genome", "query_panel", "historical_query", "u6_model", "atac_model"):
        if not c.get(key):
            raise ValueError(f"Missing required config key: {key}")
    for key in ("key", "prefix"):
        if not isinstance(c[key], str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", c[key]):
            raise ValueError(f"{key} must be a filename-safe identifier beginning with a letter")
    c.setdefault("gff", None)
    for key in PATH_KEYS:
        if c.get(key) is not None:
            value = Path(c[key])
            c[key] = str((source.parent / value).resolve(strict=True) if not value.is_absolute() else value.resolve(strict=True))
    for tool in ("blastn", "cmsearch", "cmalign"):
        value = c.get(tool, tool)
        if "/" in value and not Path(value).is_absolute():
            value = str(source.parent / value)
        c[tool] = resolve_executable(value)
    c.setdefault("threads", 2)
    c.setdefault("cm_min_score", 20.0)
    c.setdefault("excluded_gff_seqids", [])
    c.setdefault("activity_evidence", [])
    c.setdefault("reference_scope", "user_supplied_selected_reference")
    if not isinstance(c["threads"], int) or isinstance(c["threads"], bool) or c["threads"] < 1:
        raise ValueError("threads must be a positive integer")
    if not isinstance(c["excluded_gff_seqids"], list) or not all(isinstance(x, str) for x in c["excluded_gff_seqids"]):
        raise ValueError("excluded_gff_seqids must be an explicit list of sequence IDs")
    panel, historical = read_fasta(c["query_panel"]), read_fasta(c["historical_query"])
    if len(historical) != 1:
        raise ValueError("historical_query must have one historical X52528.1 transcript sequence")
    historical_id, historical_sequence = next(iter(historical.items()))
    if "X52528.1" not in historical_id or panel.get(historical_id) != historical_sequence:
        raise ValueError("Historical X52528.1 record must match the same ID and sequence in query_panel")
    if any(not seq or set(seq) - set("ACGTRYSWMKBDHVN") for seq in panel.values()):
        raise ValueError("Query FASTA must contain nonempty DNA records")
    if not isinstance(c["activity_evidence"], list):
        raise ValueError("activity_evidence must be a list of source-backed exact-body-coordinate records")
    for evidence in c["activity_evidence"]:
        if not all(k in evidence for k in ("seqid", "strand", "body_start1", "body_end1", "source", "evidence_type")):
            raise ValueError("activity_evidence needs seqid,strand,body_start1,body_end1,source,evidence_type")
        if evidence["strand"] not in {"+", "-"} or not evidence["source"] or not evidence["evidence_type"]:
            raise ValueError("Invalid activity_evidence strand/source/type")
        if evidence["evidence_type"] not in {"engineered_promoter_activity", "native_transcription", "genome_editing_with_promoter"}:
            raise ValueError("activity_evidence evidence_type must describe an experiment, not a prediction")
    return c


def search_args(c, out):
    return argparse.Namespace(genome=c["genome"], queries=c["query_panel"], u6_model=c["u6_model"],
                              atac_model=c["atac_model"], out=str(out), threads=c["threads"],
                              task="blastn-short", cm_min_score=c["cm_min_score"], dry_run=False,
                              blastn=c["blastn"], cmsearch=c["cmsearch"])


def validate_reused_search(search, c):
    """Verify complete ledger, every output hash and actual command parameters."""
    search = Path(search).resolve(strict=True)
    ledger = json.loads((search / "search_run.json").read_text())
    if ledger.get("status") != "completed" or ledger.get("dry_run"):
        raise ValueError("Reuse requires a completed non-dry-run search")
    if ledger.get("blast_scoring") != SCORING or ledger.get("blast_output_fields") != BLAST_FIELDS.split():
        raise ValueError("Search scoring/field schema mismatch; historical bit thresholds cannot be transferred")
    if ledger.get("cm_reporting_score") != c["cm_min_score"]:
        raise ValueError("Search CM reporting parameter mismatch")
    for old, current in (("genome", "genome"), ("queries", "query_panel"), ("u6_model", "u6_model"), ("atac_model", "atac_model")):
        if ledger.get("inputs", {}).get(old, {}).get("sha256") != digest(c[current]):
            raise ValueError(f"Search/extraction input mismatch: {old}")
    output_by_name = {}
    for record in ledger.get("outputs", []):
        name = Path(record["path"]).name
        if name in output_by_name:
            raise ValueError("Duplicate reused output manifest basename")
        file = search / name
        if not file.is_file() or file.stat().st_size != record["size_bytes"] or digest(file) != record["sha256"]:
            raise ValueError(f"Reused search output absent or changed: {name}")
        output_by_name[name] = record
    expected_names = {"u6_blastn_short_18col.tsv", "RF00026.tblout", "RF00026.report.txt", "RF00619.tblout", "RF00619.report.txt"}
    if not expected_names <= output_by_name.keys():
        raise ValueError("Reused search manifest is incomplete")
    commands = ledger.get("commands", [])
    if [entry.get("name") for entry in commands] != ["blastn", "RF00026", "RF00619"]:
        raise ValueError("Expected three independent whole-genome searches")
    # Rebuild, but do not execute or create, the approved command plan. Comparing
    # entire token lists prevents extra filtering/strand flags from being hidden
    # behind a superficially correct scoring dictionary.
    expected_plan = build_plan(search_args(c, search / ".unused_read_only_validation_plan"))

    def normalized_argv(command, plan):
        replacements = {record["path"]: "<input:" + key + ">" for key, record in plan["inputs"].items()}
        for record in command.get("expected_outputs", []) + command.get("optional_outputs", []):
            replacements[record] = "<output:" + Path(record).name + ">"
        argv = ["<tool:" + command["name"] + ">"] + [replacements.get(value, value) for value in command["argv"][1:]]
        for flag in ("-T", "--incT"):
            if flag in argv:
                argv[argv.index(flag) + 1] = str(float(argv[argv.index(flag) + 1]))
        return argv

    for actual, expected in zip(commands, expected_plan["commands"]):
        if normalized_argv(actual, ledger) != normalized_argv(expected, expected_plan):
            raise ValueError(f"Reused command differs from approved full argv: {actual['name']}")
        for kind in ("expected_outputs", "optional_outputs"):
            if [Path(p).name for p in actual.get(kind, [])] != [Path(p).name for p in expected.get(kind, [])]:
                raise ValueError("Reused command output contract differs")
    for command in commands:
        if command.get("exit_code") != 0 or not command.get("finished_utc"):
            raise ValueError("Search ledger has an unsuccessful or incomplete command")
        argv = command["argv"]
        required = ({"-task": "blastn-short", "-strand": "both", "-dust": "no", "-soft_masking": "false",
                     "-evalue": "1000", "-perc_identity": "60", "-qcov_hsp_perc": "50",
                     "-num_threads": str(c["threads"]), "-outfmt": "7 " + BLAST_FIELDS,
                     **{"-" + key: str(value) for key, value in SCORING.items()}}
                    if command["name"] == "blastn" else
                    {"--cpu": str(c["threads"]), "-T": str(c["cm_min_score"]), "--incT": str(c["cm_min_score"])})
        for flag, expected in required.items():
            if argv.count(flag) != 1 or argv.index(flag) + 1 >= len(argv):
                raise ValueError(f"Reused command missing/duplicate parameter: {flag}")
            actual = argv[argv.index(flag) + 1]
            if flag in {"-T", "--incT"}:
                equal = float(actual) == float(expected)
            else:
                equal = actual == expected
            if not equal:
                raise ValueError(f"Reused command parameter differs: {flag}={actual}")
        if command["name"] == "blastn":
            for flag, key in (("-query", "queries"), ("-subject", "genome")):
                if argv[argv.index(flag) + 1] != ledger["inputs"][key]["path"]:
                    raise ValueError(f"Reused {flag} path conflicts with input ledger")
            if int(argv[argv.index("-max_target_seqs") + 1]) < max(100000, ledger["inputs"]["genome"]["fasta_records"]):
                raise ValueError("Reused BLAST target cap may have omitted records")
        else:
            model_key = "u6_model" if command["name"] == "RF00026" else "atac_model"
            if argv[-2:] != [ledger["inputs"][model_key]["path"], ledger["inputs"]["genome"]["path"]]:
                raise ValueError("Reused CM did not scan the declared complete reference")
            if "--notextw" not in argv:
                raise ValueError("Reused CM command differs from frozen search profile")
        for expected in command.get("expected_outputs", []):
            if Path(expected).name not in output_by_name:
                raise ValueError("Reused command output absent from hash manifest")
    return ledger


def table(path, rows, empty_fields):
    write_tsv(path, rows, fields=None if rows else empty_fields)


def priority(row):
    """Descriptive review routing; never an activity score or family filter."""
    reasons = []
    core = row.get("core120_context_counts", {})
    if any(core.get(k, 0) for k in ("CDS", "exon", "five_prime_UTR", "three_prime_UTR", "repeat")):
        reasons.append("core_overlap_requires_interpretation")
    if row["body_integrity"] != "alignment_near_complete_not_functionally_validated":
        reasons.append("submitted_body_integrity_review")
    if row.get("multi_anchor_conflict"):
        reasons.append("multiple_anchor_evidence_review")
    if row.get("promoter500_extraction_status") != "complete":
        reasons.append("upstream_extraction_incomplete")
    if reasons:
        return "review_before_experimental_selection", reasons
    if row["activity_evidence"] != "not_tested":
        return "review_source_supported_candidate", ["provided_activity_source_not_independently_verified_by_pipeline"]
    return "candidate_for_comparison", ["activity_not_tested", "PCR_uniqueness_and_homeologs_not_assessed"]


def run(config_path, out, reuse_search=None):
    out = Path(out).absolute()
    if out.exists() or out.is_symlink():
        raise ValueError(f"Output exists; choose a new directory: {out}")
    c = load_config(config_path)
    out.mkdir(parents=True, exist_ok=False)
    for name in ("00_audit", "02_curation", "03_promoters", "04_delivery", "qa"):
        (out / name).mkdir()
    audit, cur, prom, delivery = (out / n for n in ("00_audit", "02_curation", "03_promoters", "04_delivery"))
    save_json(audit / "config.resolved.json", c)
    ledger = {"schema_version": SCRIPT_VERSION, "status": "running", "started_utc": now(),
              "config": "00_audit/config.resolved.json", "internal_qa": "not_run", "independent_qa": "not_run",
              "biological_function": "not_established_by_computation", "commands": []}
    save_json(out / "run.json", ledger)
    genome = None

    def execute(argv, name):
        entry = {"argv": list(map(str, argv)), "started_utc": now(), "log": "02_curation/" + name}
        ledger["commands"].append(entry)
        save_json(out / "run.json", ledger)
        with (cur / name).open("w") as log:
            result = subprocess.run(entry["argv"], stdout=log, stderr=subprocess.STDOUT, check=False)
        entry.update(exit_code=result.returncode, finished_utc=now())
        save_json(out / "run.json", ledger)
        if result.returncode:
            raise RuntimeError(f"Command failed: {entry['log']} (exit {result.returncode})")

    try:
        input_paths = [Path(config_path).resolve()] + [Path(c[k]) for k in PATH_KEYS if c.get(k)]
        input_paths += [Path(__file__).parent / name for name in CORE_FILES]
        manifest = [{"path": str(p), "bytes": p.stat().st_size, "sha256": digest(p), "md5": digest(p, "md5")} for p in dict.fromkeys(input_paths)]
        write_tsv(delivery / "input_manifest_with_md5.tsv", manifest)
        genome_audit = audit_fasta(c["genome"], audit / "genome.fai")
        if genome_audit["sha256"] != next(m["sha256"] for m in manifest if m["path"] == c["genome"]):
            raise ValueError("Genome changed while auditing")
        save_json(audit / "genome_audit.json", genome_audit)
        save_json(audit / "gff_audit.json", audit_gff(c["gff"], genome_audit["lengths"], c["excluded_gff_seqids"]))
        versions = capture_tool_versions(c["blastn"], c["cmsearch"])
        execute([c["cmalign"], "-h"], "cmalign_version.txt")
        save_json(audit / "versions.json", versions)
        if reuse_search:
            search = Path(reuse_search).resolve(strict=True)
            search_ledger = validate_reused_search(search, c)
            (out / "01_search").mkdir()
            save_json(out / "01_search/reuse_receipt.json", {"verified_utc": now(), "search": str(search),
                      "search_ledger_sha256": digest(search / "search_run.json"), "outputs": search_ledger["outputs"],
                      "policy": "completed commands; matching inputs, scoring and parameters; every output rehashed"})
        else:
            search = out / "01_search"
            search_ledger = build_plan(search_args(c, search))
            search_ledger["versions"] = versions
            run_plan(search_ledger)
        ledger["search_directory"] = str(search)
        ledger["search_reused"] = bool(reuse_search)
        data = curate_species(c["key"], search / "u6_blastn_short_18col.tsv", search / "RF00026.tblout", search / "RF00619.tblout", c["u6_model"], c["atac_model"])
        genome = IndexedFasta(c["genome"], audit / "genome.fai")
        rows, bodies, windows, coords = data["loci"], {}, {}, {}
        for i, row in enumerate(rows, 1):
            name = f"{c['prefix']}U6L{i:04d}"
            row.update(locus_id=name, assembly_id=c["assembly"], species=c["species"])
            seq = genome.fetch(row["seqid"], row["body_start1"], row["body_end1"])
            bodies[name] = seq if row["strand"] == "+" else revcomp(seq)
            row["body_sequence_sha256"] = hashlib.sha256(bodies[name].encode()).hexdigest()
            lo, hi = max(1, row["body_start1"] - 60), min(genome.lengths[row["seqid"]], row["body_end1"] + 60)
            seq = genome.fetch(row["seqid"], lo, hi)
            windows[name] = seq if row["strand"] == "+" else revcomp(seq)
            coords[name] = (lo, hi)
        write_fasta(cur / "all_search_bodies.fa", bodies)
        write_fasta(cur / "body_plus_60bp_windows.fa", windows)
        metrics, local_hits = {}, defaultdict(list)
        if rows:
            alignment = cur / "RF00026_bodies_global.sto"
            execute([c["cmalign"], "-g", "--dnaout", "--cpu", c["threads"], "--sfile", cur / "cmalign_scores.tsv", "-o", alignment, c["u6_model"], cur / "all_search_bodies.fa"], "cmalign.log")
            metrics = parse_stockholm_metrics(alignment)
            local = cur / "historical_query_window_blast.tsv"
            execute([c["blastn"], "-task", "blastn-short", "-query", c["historical_query"], "-subject", cur / "body_plus_60bp_windows.fa", "-strand", "plus", "-word_size", "7", "-reward", "1", "-penalty", "-3", "-gapopen", "5", "-gapextend", "2", "-dust", "no", "-soft_masking", "false", "-evalue", "1000", "-perc_identity", "60", "-qcov_hsp_perc", "50", "-max_target_seqs", str(max(100000, len(rows))), "-outfmt", "7 " + BLAST_FIELDS, "-out", local], "local_blast.log")
            by_id = {r["locus_id"]: r for r in rows}
            for line in local.read_text().splitlines():
                if not line or line.startswith("#"):
                    continue
                f = line.split("\t")
                name, qs, qe, ss, se, qlen = f[1], *map(int, (f[6], f[7], f[8], f[9], f[12]))
                row, (lo, hi) = by_id[name], coords[name]
                left, right = ((lo + ss - 1, lo + se - 1) if row["strand"] == "+" else (hi - se + 1, hi - ss + 1))
                if ss > se or overlap(left, right, row["body_start1"], row["body_end1"]) < .5 * min(se - ss + 1, row["body_length"]):
                    continue
                qcov = 100 * (qe - qs + 1) / qlen
                q1s = map_query_position(f[15], f[16], qs, ss, 1)
                q1 = None if q1s is None else lo + q1s - 1 if row["strand"] == "+" else hi - q1s + 1
                local_hits[name].append(dict(pident=float(f[2]), qcov_exact=qcov, bitscore=float(f[11]), qstart=qs, qend=qe, q1_genomic=q1, btop=f[17], qseq=f[15], sseq=f[16], strong=float(f[2]) >= 80 and qcov >= 80 and float(f[11]) >= 80))
        context = GffContext(c["gff"], genome.lengths)
        pseudogene_ids = {f.declared_id for f in context.features if f.declared_id and ("pseudogen" in f.kind.lower() or any("pseudogen" in b.lower() for b in f.own_biotypes) or any(v.lower() == "true" for v in f.attrs.get("pseudo", ())))}
        promoters, audits, pairs, boundaries, jobs, activity_matches = [], [], [], [], [], []
        for row in rows:
            name = row["locus_id"]
            self_ids = context.supporting_u6_ids(row["seqid"], row["body_start1"], row["body_end1"], row["strand"])
            pseudo = self_ids & pseudogene_ids
            row.update(classify_integrity(bodies[name], metrics.get(name), row["rf00026_model_coverage"], row["rf00026_trunc"], bool(pseudo)))
            row.update(metrics.get(name, {}))
            row.update(annotated_pseudogene=bool(pseudo) if self_ids else None, self_pseudogene_annotation_ids=sorted(pseudo), self_u6_annotation_ids=sorted(self_ids), pseudogene_annotation_scope="self_U6_annotation_only; absent_self_annotation_is_unassessed", homeolog_group=None, PCR_uniqueness=None)
            ctx = context.audit(row["seqid"], row["body_start1"], row["body_end1"], row["strand"], self_ids=self_ids)
            row.update(body_context_counts=ctx["counts_by_type"], body_context_status=ctx["status"])
            audits.append(dict(locus_id=name, region="body", **ctx))
            jobs.append(dict(job_id="B__" + name, seqid=row["seqid"], start1=row["body_start1"], end1=row["body_end1"], strand=row["strand"], sequence=bodies[name]))
            hsps = sorted(local_hits[name], key=lambda h: (h["strong"], h["bitscore"], h["qcov_exact"]), reverse=True)
            best = hsps[0] if hsps else None
            direct = best["q1_genomic"] if best and best["strong"] else None
            cm_boundary = row["body_start1"] if row["strand"] == "+" else row["body_end1"]
            tss = direct if direct is not None else cm_boundary
            primary = "aligned_historical_q1" if direct is not None else "search_body_5prime_hypothesis"
            plus1 = genome.fetch(row["seqid"], tss, tss)
            row.update(predicted_tss1=tss, boundary_confidence=primary, terminal_5prime_unobserved=direct is None, plus1_base=plus1 if row["strand"] == "+" else revcomp(plus1), boundary_interpretation="hypothesis_not_measured_native_TSS")
            hypotheses = {primary: tss, "search_body_5prime_hypothesis": cm_boundary}
            boundaries.append(dict(locus_id=name, primary_hypothesis=primary, primary_tss1=tss, hypotheses=hypotheses, local_best_hsp=best))
            matches = [e for e in c["activity_evidence"] if all(e[k] == row[k] for k in ("seqid", "strand", "body_start1", "body_end1"))]
            row["activity_evidence"] = "provided_source_support" if matches else "not_tested"
            row["activity_sources"] = matches
            activity_matches.extend(dict(locus_id=name, **e) for e in matches)
            if not row["canonical_retained"]:
                row["experimental_review_priority"] = "family_review_required"
                continue
            for label, boundary in hypotheses.items():
                for length in ((300, 500, 1000) if label == primary else (500,)):
                    region = extract_upstream(genome.fetch, row["seqid"], row["strand"], boundary, length, genome.lengths[row["seqid"]])
                    record = dict(locus_id=name, species_key=c["key"], hypothesis=label, is_primary_hypothesis=label == primary, **region)
                    promoters.append(record)
                    if region["actual_length"]:
                        jobs.append(dict(job_id=f"P__{name}__{label}__{length}", seqid=row["seqid"], start1=region["start1"], end1=region["end1"], strand=row["strand"], sequence=region["sequence"]))
                    if length != 500:
                        continue
                    if label == primary:
                        row["promoter500_extraction_status"] = region["extraction_status"]
                        row["promoter500_start1"], row["promoter500_end1"] = region["start1"], region["end1"]
                    if not region["actual_length"]:
                        continue
                    ctx = context.audit(row["seqid"], region["start1"], region["end1"], row["strand"], self_ids=self_ids)
                    audits.append(dict(locus_id=name, region="promoter_500", hypothesis=label, **ctx))
                    record.update(context_counts=ctx["counts_by_type"], context_status=ctx["status"])
                    for profile in ("legacy_arabidopsis", "generic_provisional"):
                        motif = scan_motifs(region["sequence"], region["first_relative_position"], profile, msp=True)
                        record[profile] = motif
                        pairs.extend(dict(locus_id=name, hypothesis=label, profile=profile, **p) for p in motif["pairs"])
                        if label == primary:
                            row[profile + "_class"] = motif["best"].get("motif_pair_class", "no_reportable_pair")
                            row[profile + "_best_pair"] = motif["best"]
                            row[profile + "_MSP_count"] = len(motif["msp_hits"])
                    if label == primary:
                        row.update(promoter500_context_counts=ctx["counts_by_type"], promoter500_context_status=ctx["status"])
                if label == primary:
                    core = extract_upstream(genome.fetch, row["seqid"], row["strand"], boundary, 120, genome.lengths[row["seqid"]])
                    if core["actual_length"]:
                        ctx = context.audit(row["seqid"], core["start1"], core["end1"], row["strand"], self_ids=self_ids)
                        row["core120_context_counts"] = ctx["counts_by_type"]
                        audits.append(dict(locus_id=name, region="core120", hypothesis=label, **ctx))
            if row["strand"] == "+":
                terminal = genome.fetch(row["seqid"], tss, min(genome.lengths[row["seqid"]], row["body_end1"] + 60))
                span = row["body_end1"] - tss + 1
            else:
                terminal = revcomp(genome.fetch(row["seqid"], max(1, row["body_start1"] - 60), tss))
                span = tss - row["body_start1"] + 1
            row.update(scan_termination(terminal, span))
            row["activity_geometry_interpretation"] = ("experiment_supported_noncanonical_geometry" if matches and row.get("generic_provisional_class") != "provisional_type3_like_geometry" else "source_support_and_geometry_are_separate" if matches else "activity_not_tested")
            row["experimental_review_priority"], row["review_priority_reasons"] = priority(row)
        selected = [r for r in rows if r["canonical_retained"]]
        table(cur / "all_search_loci.tsv", rows, LOCUS_FIELDS)
        table(cur / "raw_evidence.tsv", data["evidence"], ["raw_hit_id", "method", "query", "seqid", "start1", "end1", "strand"])
        table(cur / "U6atac_separate_evidence.tsv", data["standalone_atac_hits"], ["raw_hit_id", "seqid", "start1", "end1", "strand", "ga_pass"])
        table(cur / "boundary_hypotheses.tsv", boundaries, ["locus_id", "primary_hypothesis", "primary_tss1", "hypotheses", "local_best_hsp"])
        table(prom / "promoter_sequences.tsv", promoters, PROMOTER_FIELDS)
        table(prom / "annotation_context.tsv", audits, ["locus_id", "region", "hypothesis", "status", "counts_by_type"])
        table(prom / "all_reportable_motif_pairs.tsv", pairs, ["locus_id", "hypothesis", "profile", "motif_pair_class"])
        table(out / "qa/extraction_jobs.tsv", jobs, ["job_id", "seqid", "start1", "end1", "strand", "sequence"])
        save_json(out / "qa/activity_evidence_mapping.json", {"provided": c["activity_evidence"], "matched": activity_matches,
                  "interpretation": "Exact body-coordinate match; source claims supplied externally, not independently validated here"})
        table(delivery / f"{c['key']}_U6_candidates.tsv", selected, LOCUS_FIELDS)
        table(delivery / f"{c['key']}_U6_strict_AtU6_subset.tsv", [r for r in selected if r["strict_historical_query_pass"]], LOCUS_FIELDS)
        review = [dict(evidence_set="noncanonical_locus", **r) for r in rows if not r["canonical_retained"]]
        review += [dict(evidence_set="standalone_U6atac_hit", **h) for h in data["standalone_atac_hits"]]
        table(delivery / "U6atac_and_review.tsv", review, ["evidence_set", "seqid", "strand", "family_class"])
        write_fasta(delivery / f"{c['key']}_U6_body.fa", {r["locus_id"]: bodies[r["locus_id"]] for r in selected})
        with (delivery / f"{c['key']}_U6_loci.bed").open("w") as handle:
            for r in selected:
                handle.write(f"{r['seqid']}\t{r['body_start1']-1}\t{r['body_end1']}\t{r['locus_id']}\t0\t{r['strand']}\n")
        for length in (300, 500, 1000):
            write_fasta(delivery / f"{c['key']}_U6_promoter_{length}.fa", {p["locus_id"]: p["sequence"] for p in promoters if p["is_primary_hypothesis"] and p["requested_length"] == length})
        checks = {"all_bodies_aligned": set(bodies) == set(metrics), "unique_locus_ids": len(bodies) == len(rows),
                  "all_family_have_three_primary_promoter_records": sum(p["is_primary_hypothesis"] for p in promoters) == 3 * len(selected),
                  "inputs_and_code_preserved": all(digest(m["path"]) == m["sha256"] for m in manifest),
                  "all_body_lengths_match": all(len(bodies[r["locus_id"]]) == r["body_length"] for r in rows)}
        summary = {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks,
                   "species": c["species"], "assembly": c["assembly"], "discovery": data["diagnostics"],
                   "counts": {"all_search_loci": len(rows), "RF00026_supported": len(selected),
                              "strict_AtU6": sum(r["strict_historical_query_pass"] for r in selected),
                              "integrity": dict(Counter(r["body_integrity"] for r in selected)),
                              "generic_geometry": dict(Counter(r.get("generic_provisional_class", "unassessed") for r in selected)),
                              "extraction_jobs": len(jobs)},
                   "independent_qa": "not_run", "function_or_expression_validation": "not_performed"}
        save_json(out / "qa/summary.json", summary)
        if not all(checks.values()):
            raise RuntimeError("Internal QA failed; inspect qa/summary.json")
        write_reports(out, c, summary, ledger, search_ledger, versions)
        ledger.update(status="completed", internal_qa="PASS", finished_utc=now(), summary="qa/summary.json")
        save_json(out / "run.json", ledger)
        outputs = [{"path": str(p.relative_to(out)), "bytes": p.stat().st_size, "sha256": digest(p)} for p in sorted(out.rglob("*")) if p.is_file()]
        write_tsv(out / "output_manifest.tsv", outputs)
        return summary
    except (Exception, KeyboardInterrupt) as error:
        ledger.update(status="failed", error=str(error), finished_utc=now())
        save_json(out / "run.json", ledger)
        raise
    finally:
        if genome is not None:
            genome.close()


def write_reports(out, c, summary, ledger, search_ledger, versions):
    d, counts = out / "04_delivery", summary["counts"]
    report = f"""# U6 candidate report: {c['species']}

Assembly: {c['assembly']}. Reference scope: {c['reference_scope']}. All {search_ledger['inputs']['genome']['fasta_records']} supplied FASTA records were searched. Nuclear identity and reference selection are declared upstream, not inferred by file-format checks.

## Verification layers

- Computation and internal consistency: PASS. This is not independent extraction QA.
- Independent sequence/coordinate validation is a separate stage. This immutable report does not track its later status: use the separately linked independent_validation.json receipt. Run validate_outputs.py if that receipt is not yet available.
- Promoter function, expression, editing efficiency, native TSS, PCR uniqueness and homeolog assignment: not established by this computation.

## Results and scope

{counts['all_search_loci']} merged BLAST/RF00026 search loci include {counts['RF00026_supported']} RF00026-supported, non-U6atac-ambiguous family candidates. The strict historical AtU6 subset contains {counts['strict_AtU6']} candidates and is not the full family. Separate U6atac evidence and weak/ambiguous loci are retained in U6atac_and_review.tsv. Empty sets are explicit, not evidence of biological absence.

Separate RF00619 GA-supported hits: {summary['discovery']['standalone_atac_GA_count']}. All locus family classes: {json.dumps(summary['discovery']['family_counts'])}. Standalone hit counts are evidence records, not an independently curated U6atac gene catalogue.

Body integrity classes: {json.dumps(counts['integrity'], ensure_ascii=False)}.

Generic provisional geometry classes: {json.dumps(counts['generic_geometry'], ensure_ascii=False)}. These are Arabidopsis-derived heuristics, not calibrated plant-wide promoter-activity classes. Motifs never remove a family candidate. The known Kitaake OsU6a counterexample has USE/TATA-like sequences but noncanonical edge distance 20 (gap 19); source support and geometry remain separate.

## Interpretation boundaries

Candidate FASTA/BED files cover the whole family-supported set. 300/500/1000-bp promoter exports are transcript-oriented and exclude the predicted +1. Empty or shortened windows at assembly ends are retained and described in 03_promoters/promoter_sequences.tsv. A query-aligned q1 or search-body boundary is a hypothesis, not measured native TSS. Legacy provisional extrapolation diagnostics in raw locus evidence are not used as primary TSS coordinates.

Missing GFF, UTR, repeats or annotation categories are NA, never a clean/no-conflict claim. Self-U6 annotations are excluded using parent-aware provenance; foreign overlaps remain descriptive review evidence and do not by themselves invalidate family membership. Review-priority categories are not expression scores or an experiment-ready shortlist. Source-backed activity metadata, if supplied, is mapped by exact body coordinates, retained verbatim and not independently source-validated by this script.

Independent Infernal scans retain heuristic filters; full-reference scanning is not proof that every functional U6 was recovered. Refer to raw evidence, exact commands, input hashes and per-locus motif/context records when interpreting candidates.
"""
    (d / "U6_candidate_report.md").write_text(report)
    (d / "materials_and_methods.md").write_text(f"""# Materials and methods (computational draft)

The {c['species']} assembly {c['assembly']} was analyzed using all records in the supplied selected-reference FASTA (scope: {c['reference_scope']}). Nuclear-reference selection is an upstream source-audit decision and is not established by these file-format checks. FASTA sequence IDs, DNA alphabet and record lengths were checked, and a new index was generated outside the source directory. Supplied GFF3 sequence IDs and included-record coordinates were checked against the FASTA; explicitly excluded non-reference IDs were logged. Missing annotation classes were treated as unassessed.

The frozen nucleotide-query panel, including the historical AtU6 transcript identified by X52528.1, was independently searched with BLASTN-short. The explicit scoring profile was reward 1, mismatch penalty -3, gap opening 5, gap extension 2 and word size 7; dust and soft masking were disabled, the reporting E-value was 1000, and reporting identity/coverage limits were 60%/50%. The strict historical subset required one HSP meeting 80% identity, 80% exact query-coordinate coverage and 80 bits under this scoring profile. It was not used as the sole family-membership criterion.

RF00026 and RF00619 covariance models independently scanned every supplied FASTA record with Infernal at reporting score {c['cm_min_score']}; model-specific GA thresholds were read from the supplied files. Fixed-anchor locus curation separated U6, overlapping GA-supported U6atac competition and weaker evidence. All raw hits and ambiguous assignments were preserved. RF00026-supported loci without overlapping GA-supported RF00619 competition formed the family candidate set regardless of promoter motif classification.

Submitted body sequences were globally aligned to RF00026 with cmalign to assess model-position occupancy, missing termini, internal gaps and pair compatibility; these are sequence-review flags, not proven native pseudogene classifications. Actual aligned historical-query position 1 was used only when supported by a strong local HSP. Otherwise the oriented search-body 5-prime end remained an explicit boundary hypothesis. Strand-correct 300-, 500- and 1000-bp upstream windows excluded +1 and were truncated without padding at reference boundaries. USE/TATA pairs, four distance measures, optional descriptive MSP motifs and positional poly-T evidence were recorded using the frozen legacy_arabidopsis and generic_provisional profiles. GFF context was audited with self-U6 ancestry exclusion and missingness retained. Exact software version output, commands and input checksums accompany this draft.

This computation did not measure RNA expression, promoter activity, editing efficiency or native transcription start sites. Independent sequence extraction QA must be cited separately when completed.
""")
    (d / "results.md").write_text(f"""# Results (computational draft)

Independent nucleotide-homology and covariance-model searches of {c['species']} ({c['assembly']}) yielded {counts['all_search_loci']} merged search loci. Of these, {counts['RF00026_supported']} were supported by RF00026 without competing same-locus GA-supported RF00619 evidence and were retained as U6-family candidates. The strict historical AtU6 subset contained {counts['strict_AtU6']} loci. Weak, U6atac-like and ambiguous evidence was retained separately rather than incorporated into this family count.

Each family candidate received an oriented body sequence and three standard upstream extraction records. Submitted-body integrity categories were {json.dumps(counts['integrity'])}. Generic provisional motif geometry categories were {json.dumps(counts['generic_geometry'])}; they describe sequence arrangement under an Arabidopsis-derived profile and do not demonstrate promoter activity. Candidate-specific coordinates, extraction status, motifs and annotation context are provided in the accompanying tables. Functional activity and relative promoter strength remain untested by this analysis.
""")
    (d / "Methods_and_Results_current_stage.md").write_text((d / "materials_and_methods.md").read_text() + "\n" + (d / "results.md").read_text())
    lines = ["Exact executed argv (shell-display only; no shell was used):", ""]
    for command in search_ledger["commands"] + ledger["commands"]:
        lines.append(shlex.join(command["argv"]))
    lines += ["", "Search-stage versions (reuse retains original search versions):", json.dumps(search_ledger.get("versions", {}), indent=2),
              "", "Current postprocessing versions:", json.dumps(versions, indent=2), (out / "02_curation/cmalign_version.txt").read_text()]
    (d / "commands_and_versions.txt").write_text("\n".join(lines) + "\n")
    (out / "README_where_to_find.md").write_text(f"""# Where to find this run

Start with [U6 candidate report](04_delivery/U6_candidate_report.md). Computation/internal QA has passed; independent QA must be run separately and biological activity is not established.

- `04_delivery/{c['key']}_U6_candidates.tsv`: whole family-supported set, NOT the strict At subset or an experimental shortlist.
- `04_delivery/{c['key']}_U6_strict_AtU6_subset.tsv`: historical single-query subset only.
- `04_delivery/{c['key']}_U6_loci.bed`, `_U6_body.fa`, `_U6_promoter_300.fa`, `_500.fa`, `_1000.fa`: family-supported sequences/coordinates; promoter windows exclude +1.
- `04_delivery/U6atac_and_review.tsv`: weak, ambiguous and separate U6atac evidence.
- `02_curation/`: all search loci, all bodies, global alignment, raw evidence and boundary hypotheses.
- `03_promoters/`: every extraction record, motif pairs and raw parent-aware annotation overlap evidence.
- `00_audit/`: resolved config, independently built index, FASTA/GFF audits and software versions.
- `01_search/`: raw independent whole-genome searches or verified reuse receipt pointing to preserved originals.
- `qa/summary.json`: internal checks only; `run.json`: completion/failure ledger.
- `04_delivery/materials_and_methods.md`, `results.md`: editable evidence-bounded manuscript drafts.
- [Input manifest](04_delivery/input_manifest_with_md5.tsv), [commands and versions](04_delivery/commands_and_versions.txt), [output manifest](output_manifest.tsv): provenance and reproducibility. The output manifest is at the run root, not in 04_delivery.

No files in source reference directories were modified. A interrupted/failed run remains preserved; use a new output directory. Reuse only the complete search stage after exact input/parameter/output hash checks via --reuse-search. Completion of an old directory alone is insufficient.
""")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="JSON config; relative input paths are relative to this file")
    parser.add_argument("--out", required=True, help="New output directory; existing paths refused")
    parser.add_argument("--reuse-search", help="Completed search directory with matching ledger and output hashes")
    args = parser.parse_args(argv)
    try:
        summary = run(args.config, args.out, args.reuse_search)
        print(json.dumps(summary, indent=2))
        return 0
    except (OSError, ValueError, RuntimeError, KeyError, subprocess.SubprocessError) as error:
        print(f"U6 pipeline error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
