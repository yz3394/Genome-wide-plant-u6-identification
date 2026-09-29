#!/usr/bin/env python3
"""Independent BLASTN and U6/U6atac CM searches against a supplied genome.

The caller must freeze and audit the nuclear FASTA beforehand. This program
does not download, mask, rename, exclude scaffolds, infer expression or apply
candidate thresholds. Only the historical blastn-short scoring profile is
supported; 80-bit curation rules must not be transferred to another profile.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from u6_loci import read_cm_metadata


BLAST_FIELDS = (
    "qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
    "evalue bitscore qlen slen qcovhsp qseq sseq btop"
)
SCORING = {"reward": 1, "penalty": -3, "gapopen": 5, "gapextend": 2, "word_size": 7}


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def file_manifest(path, *, fasta=False):
    """Hash exact bytes; also count FASTA records without a second genome pass."""
    path = Path(path).resolve(strict=True)
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"Input must be a nonempty regular file: {path}")
    digest, records, previous, first = hashlib.sha256(), 0, b"", True
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
            if fasta:
                if first and not block.startswith(b">"):
                    raise ValueError(f"Expected uncompressed FASTA beginning with >: {path}")
                records += (previous + block).count(b"\n>") + int(first and block.startswith(b">"))
                previous = block[-1:]
            first = False
    result = {"path": str(path), "size_bytes": path.stat().st_size, "sha256": digest.hexdigest()}
    if fasta:
        result["fasta_records"] = records
    return result


def resolve_executable(value):
    executable = shutil.which(value)
    if not executable:
        raise ValueError(f"Executable not found or not executable: {value}; provide its path, no installation is attempted")
    return str(Path(executable).resolve())


def capture_tool_versions(blastn, cmsearch):
    records = {}
    for name, command in [("blastn", [blastn, "-version"]), ("cmsearch", [cmsearch, "-h"])]:
        started = utc_now()
        result = subprocess.run(command, capture_output=True, text=True, timeout=30, check=False)
        records[name] = {
            "argv": command, "started_utc": started, "finished_utc": utc_now(),
            "exit_code": result.returncode, "stdout": result.stdout, "stderr": result.stderr,
        }
        if result.returncode:
            raise ValueError(f"Cannot obtain {name} version: exit {result.returncode}: {result.stderr.strip()}")
    return records


def build_plan(args):
    out = Path(args.out).absolute()
    if out.exists() or out.is_symlink():
        raise ValueError(f"Output path already exists; choose a new directory: {out}")
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    if args.task != "blastn-short":
        raise ValueError("Only the frozen blastn-short profile is supported; other scores require separate validation")
    if not math.isfinite(args.cm_min_score):
        raise ValueError("--cm-min-score must be finite")
    inputs = {
        "genome": file_manifest(args.genome, fasta=True),
        "queries": file_manifest(args.queries, fasta=True),
        "u6_model": file_manifest(args.u6_model),
        "atac_model": file_manifest(args.atac_model),
    }
    models = {
        "RF00026": read_cm_metadata(inputs["u6_model"]["path"]),
        "RF00619": read_cm_metadata(inputs["atac_model"]["path"]),
    }
    for expected, metadata in models.items():
        if metadata["accession"] != expected:
            raise ValueError(f"Expected {expected}, found {metadata['accession']} in {metadata['source']}")
        # One search per model must not silently become a multi-model search.
        with Path(metadata["source"]).open() as handle:
            count = sum(line.startswith("INFERNAL1/") for line in handle)
        if count != 1:
            raise ValueError(f"Expected one Infernal CM per file, found {count}: {metadata['source']}")
        if args.cm_min_score > metadata["ga"]:
            raise ValueError(f"Reporting score {args.cm_min_score} exceeds {expected} GA {metadata['ga']}; this would hide GA-supported hits")
    blastn, cmsearch = resolve_executable(args.blastn), resolve_executable(args.cmsearch)
    genome, queries = inputs["genome"]["path"], inputs["queries"]["path"]
    blast_output = out / "u6_blastn_short_18col.tsv"
    # No subject scaffold is excluded solely by the default output-target cap.
    maximum_targets = max(100000, inputs["genome"]["fasta_records"])
    blast = [
        blastn, "-task", "blastn-short", "-query", queries, "-subject", genome,
        "-strand", "both", "-dust", "no", "-soft_masking", "false",
        "-evalue", "1000", "-perc_identity", "60", "-qcov_hsp_perc", "50",
        "-max_target_seqs", str(maximum_targets), "-num_threads", str(args.threads),
    ]
    for key, value in SCORING.items():
        blast += ["-" + key, str(value)]
    blast += ["-outfmt", "7 " + BLAST_FIELDS, "-out", str(blast_output)]
    commands = [{"name": "blastn", "argv": blast, "expected_outputs": [str(blast_output)]}]
    for accession, metadata in models.items():
        table, report, alignment = (out / f"{accession}.{suffix}" for suffix in ("tblout", "report.txt", "alignments.sto"))
        command = [
            cmsearch, "--cpu", str(args.threads), "--notextw", "-T", str(args.cm_min_score),
            "--incT", str(args.cm_min_score), "--tblout", str(table), "-o", str(report),
            "-A", str(alignment), metadata["source"], genome,
        ]
        # With zero significant hits Infernal may omit an alignment file.
        commands.append({"name": accession, "argv": command,
                         "expected_outputs": [str(table), str(report)],
                         "optional_outputs": [str(alignment)]})
    for command in commands:
        command["shell_display"] = shlex.join(command["argv"])
    return {
        "schema_version": "1.0", "created_utc": utc_now(), "dry_run": args.dry_run,
        "output_directory": str(out), "status": "planned", "inputs": inputs,
        "models": models, "commands": commands, "blast_scoring": SCORING.copy(),
        "blast_profile": "historical_blastn_short_explicit_1_minus3_gap5_2_word7",
        "scoring_validation": "BLAST 2.16.0+ local query-panel self-search pairwise footer confirmed matrix 1 -3 and gap 5/2 on 2026-09-28",
        "blast_output_fields": BLAST_FIELDS.split(), "cm_reporting_score": args.cm_min_score,
        "scope": "Each of RF00026, RF00619 and BLASTN independently scans every supplied FASTA record on both strands",
        "evidence_limits": [
            "Input FASTA nuclear content and assembly quality require a separate reference audit",
            "CM reporting cutoff is not GA; GA and CLEN are read from each frozen model",
            "Infernal heuristic filters remain enabled; absence is not an exhaustive negative proof",
            "No expression, promoter activity, pseudogene or candidate-grade decision is made here",
            "Historical bit-score thresholds are not validated for a different BLAST scoring profile",
        ],
        "tools": {"blastn": blastn, "cmsearch": cmsearch},
    }


def run_plan(plan):
    out = Path(plan["output_directory"])
    out.mkdir(parents=True, exist_ok=False)
    ledger_path = out / "search_run.json"
    plan["status"], plan["started_utc"] = "running", utc_now()

    def save():
        ledger_path.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + "\n")

    save()
    try:
        for command in plan["commands"]:
            command["started_utc"] = utc_now()
            stdout_path, stderr_path = out / f"{command['name']}.stdout.log", out / f"{command['name']}.stderr.log"
            command["stdout_path"], command["stderr_path"] = str(stdout_path), str(stderr_path)
            save()
            with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
                result = subprocess.run(command["argv"], stdout=stdout, stderr=stderr, check=False)
            command["exit_code"], command["finished_utc"] = result.returncode, utc_now()
            save()
            if result.returncode:
                raise RuntimeError(f"{command['name']} exited {result.returncode}; inspect {stderr_path}")
            missing = [path for path in command["expected_outputs"] if not Path(path).is_file()]
            if missing:
                raise RuntimeError(f"Command succeeded but expected outputs are missing: {missing}")
        plan["status"] = "completed"
        plan["outputs"] = []
        for path in sorted(out.iterdir()):
            if path != ledger_path and path.is_file():
                if path.stat().st_size:
                    plan["outputs"].append(file_manifest(path))
                else:
                    plan["outputs"].append({"path": str(path), "size_bytes": 0, "sha256": hashlib.sha256(b"").hexdigest()})
    except (Exception, KeyboardInterrupt) as error:
        plan["status"], plan["error"] = "failed", str(error)
        raise
    finally:
        plan["finished_utc"] = utc_now()
        save()
    return plan


def argument_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("genome", "queries", "u6-model", "atac-model", "out"):
        parser.add_argument("--" + option, required=True)
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--cmsearch", default="cmsearch")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--task", choices=["blastn-short"], default="blastn-short")
    parser.add_argument("--cm-min-score", type=float, default=20.0)
    parser.add_argument("--dry-run", action="store_true", help="Validate and print the plan without creating --out or running searches")
    return parser


def main(argv=None):
    args = argument_parser().parse_args(argv)
    try:
        plan = build_plan(args)
        plan["versions"] = capture_tool_versions(plan["tools"]["blastn"], plan["tools"]["cmsearch"])
        if not args.dry_run:
            run_plan(plan)
        print(json.dumps(plan, indent=2, ensure_ascii=False))
        return 0
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"U6 search error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

