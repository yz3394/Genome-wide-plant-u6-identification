#!/usr/bin/env python3
"""Independently audit a completed U6 run with read-only SeqKit extraction.

This standalone script imports no pipeline module. It rebuilds extraction jobs
from coordinates, checks exports and frozen inputs, and streams the reference
to SeqKit without creating or reading an index alongside the original genome.
PASS is computational consistency, not proof of U6/promoter biological activity.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def digest(path):
    hasher = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            hasher.update(block)
    return hasher.hexdigest()


def seq_digest(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def read_table(path):
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError('Missing TSV header: ' + str(path))
        return list(reader)


def read_fasta(path, bed_ids=None):
    sequences, current = {}, None
    with Path(path).open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith('>'):
                words = line[1:].split()
                if not words:
                    raise ValueError('Empty FASTA identifier: ' + str(path))
                if bed_ids is None:
                    current = words[0]
                else:
                    matches = [word for word in words if word in bed_ids]
                    if len(matches) != 1:
                        raise ValueError('No unique BED job ID in SeqKit header: ' + line)
                    current = matches[0]
                if current in sequences:
                    raise ValueError('Duplicate FASTA identifier: ' + current)
                sequences[current] = ''
            elif line:
                if current is None:
                    raise ValueError('Sequence before FASTA header: ' + str(path))
                sequences[current] += line.upper()
    return sequences


def genome_lengths(path):
    """Independently count bases without trusting a pipeline FASTA index."""
    lengths, current = {}, None
    with Path(path).open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith('>'):
                words = line[1:].split()
                if not words or words[0] in lengths:
                    raise ValueError('Empty or duplicate genome FASTA identifier')
                current = words[0]
                lengths[current] = 0
            elif line:
                if current is None or any(c.isspace() for c in line):
                    raise ValueError('Malformed genome FASTA sequence')
                lengths[current] += len(line)
    if not lengths or any(n == 0 for n in lengths.values()):
        raise ValueError('Empty genome FASTA or empty record')
    return lengths


def truth(value):
    normalized = str(value).lower()
    if normalized not in {'true', '1', 'yes', 'false', '0', 'no'}:
        raise ValueError('Invalid boolean value: ' + str(value))
    return normalized in {'true', '1', 'yes'}


def utc():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def write_table(path, rows, fields):
    with Path(path).open('w') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def validate(args):
    run = Path(args.run).resolve(strict=True)
    out = Path(args.out).absolute()
    if out.exists() or out.is_symlink():
        raise ValueError('Output already exists; use a fresh QA directory')
    executable = shutil.which(args.seqkit)
    if not executable:
        raise ValueError('SeqKit is unavailable; supply --seqkit; no installation attempted')
    executable = str(Path(executable).resolve())
    out.mkdir(parents=True, exist_ok=False)
    summary = dict(status='RUNNING', started_utc=utc(), run=str(run),
                   validator_sha256=digest(__file__), checks=[], commands=[],
                   scope='Independent coordinates, sequence exports, family bookkeeping and input preservation; not biological activity validation')
    before, index_before = {}, {}

    def check(name, condition, **detail):
        summary['checks'].append(dict(name=name, passed=bool(condition), **detail))

    def save():
        (out / 'independent_validation.json').write_text(json.dumps(summary, indent=2) + '\n')

    def protect(path):
        path = Path(path).resolve(strict=True)
        if str(path) not in before:
            before[str(path)] = digest(path)
        return path

    try:
        config_path = protect(run / '00_audit/config.resolved.json')
        config = json.loads(config_path.read_text())
        key = config['key']
        if not key or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in key):
            raise ValueError('Unsafe species key in config')
        genome = protect(config['genome'])
        summary['genome'] = str(genome)
        for suffix in ('.fai', '.seqkit.fai', '.gzi'):
            index = Path(str(genome) + suffix)
            index_before[str(index)] = digest(index) if index.is_file() else None
        ledger_path = protect(run / 'run.json')
        ledger = json.loads(ledger_path.read_text())
        check('pipeline_run_completed', ledger.get('status') == 'completed')
        manifest_path = protect(run / '04_delivery/input_manifest_with_md5.tsv')
        manifest = read_table(manifest_path)
        manifest_paths = []
        for row in manifest:
            path = protect(row['path'])
            manifest_paths.append(path)
            check('frozen_input_sha256:' + str(path), before[str(path)] == row['sha256'])
        check('manifest_has_one_genome_entry', manifest_paths.count(genome) == 1)
        check('manifest_paths_unique', len(manifest_paths) == len(set(manifest_paths)))
        lengths = genome_lengths(genome)
        paths = {
            'loci': run / '02_curation/all_search_loci.tsv',
            'bodies': run / '02_curation/all_search_bodies.fa',
            'promoters': run / '03_promoters/promoter_sequences.tsv',
            'jobs': run / 'qa/extraction_jobs.tsv',
            'candidates': run / '04_delivery' / f'{key}_U6_candidates.tsv',
            'family_bodies': run / '04_delivery' / f'{key}_U6_body.fa',
            'bed': run / '04_delivery' / f'{key}_U6_loci.bed',
        }
        paths.update({f'promoter_{n}': run / '04_delivery' / f'{key}_U6_promoter_{n}.fa' for n in (300, 500, 1000)})
        for path in paths.values():
            protect(path)
        loci = read_table(paths['loci'])
        promoters = read_table(paths['promoters'])
        bodies = read_fasta(paths['bodies'])
        by_id = {row['locus_id']: row for row in loci}
        check('unique_locus_ids', len(by_id) == len(loci))
        check('all_search_body_fasta_ids', set(bodies) == set(by_id))
        jobs = {}
        for name, row in by_id.items():
            sequence = bodies.get(name, '')
            start, end = int(row['body_start1']), int(row['body_end1'])
            check('body_geometry_and_digest:' + name,
                  1 <= start <= end <= lengths.get(row['seqid'], 0) and
                  len(sequence) == end - start + 1 == int(row['body_length']) and
                  seq_digest(sequence) == row['body_sequence_sha256'] and row['strand'] in {'+', '-'})
            if 'rf00026_ga_pass' in row and 'rf00619_ga_pass' in row:
                check('family_retention_boolean_consistency:' + name,
                      truth(row['canonical_retained']) == (truth(row['rf00026_ga_pass']) and not truth(row['rf00619_ga_pass'])))
            jobs['B__' + name] = dict(seqid=row['seqid'], start1=start, end1=end, strand=row['strand'], sequence=sequence)
        family = {name for name, row in by_id.items() if truth(row['canonical_retained'])}
        primary, promoter_keys, empty_count = {}, set(), 0
        for row in promoters:
            name, requested, tss = row['locus_id'], int(row['requested_length']), int(row['tss1'])
            pk = name, row['hypothesis'], requested
            if pk in promoter_keys:
                raise ValueError('Duplicate promoter key: ' + str(pk))
            promoter_keys.add(pk)
            check('promoter_locus_scope:' + str(pk), name in family and row['seqid'] == by_id[name]['seqid'] and row['strand'] == by_id[name]['strand'])
            if truth(row['is_primary_hypothesis']):
                if (name, requested) in primary:
                    raise ValueError('Duplicate primary promoter window')
                primary[name, requested] = row
            sequence, actual = row['sequence'].upper(), int(row['actual_length'])
            n = lengths.get(row['seqid'], 0)
            expected = ((max(1, tss - requested), tss - 1) if row['strand'] == '+'
                        else (tss + 1, min(n, tss + requested)))
            if not actual:
                empty_count += 1
                reason_valid = ((not 1 <= tss <= n and row['extraction_status'] == 'boundary_outside_reference') or
                                (1 <= tss <= n and expected[0] > expected[1] and row['extraction_status'] == 'no_upstream_sequence_at_reference_edge'))
                check('empty_promoter_geometry:' + str(pk), not sequence and reason_valid and truth(row['truncated']))
                continue
            start, end = int(row['start1']), int(row['end1'])
            check('promoter_geometry:' + str(pk),
                  requested > 0 and row['strand'] in {'+', '-'} and 1 <= tss <= n and
                  (start, end) == expected and 1 <= start <= end <= n and
                  len(sequence) == actual == end - start + 1 and actual <= requested and
                  int(row['first_relative_position']) == -actual and
                  truth(row['truncated']) == (actual < requested) and int(row['N_count']) == sequence.count('N'))
            jobs[f'P__{name}__{row["hypothesis"]}__{requested}'] = dict(seqid=row['seqid'], start1=start, end1=end, strand=row['strand'], sequence=sequence)
        supplied_jobs = read_table(paths['jobs'])
        observed_jobs = {row['job_id']: dict(seqid=row['seqid'], start1=int(row['start1']), end1=int(row['end1']), strand=row['strand'], sequence=row['sequence'].upper()) for row in supplied_jobs}
        check('driver_jobs_equal_independently_reconstructed_jobs', len(supplied_jobs) == len(observed_jobs) and observed_jobs == jobs,
              reconstructed_count=len(jobs), driver_count=len(supplied_jobs))
        candidates = read_table(paths['candidates'])
        check('candidate_table_ids', len(candidates) == len(family) and {row['locus_id'] for row in candidates} == family)
        check('candidate_shared_fields_match_all_loci', all(all(row[k] == by_id.get(row['locus_id'], {}).get(k) for k in row if k in by_id.get(row['locus_id'], {})) for row in candidates))
        check('family_body_fasta_exact', read_fasta(paths['family_bodies']) == {name: bodies[name] for name in family})
        beds = [line.rstrip().split('\t') for line in paths['bed'].read_text().splitlines() if line.strip()]
        expected_beds = {(r['seqid'], str(int(r['body_start1']) - 1), r['body_end1'], name, '0', r['strand']) for name, r in by_id.items() if name in family}
        check('family_BED_exact', len(beds) == len(family) and {tuple(bed) for bed in beds} == expected_beds)
        check('primary_promoters_exact_scope', set(primary) == {(name, length) for name in family for length in (300, 500, 1000)})
        for length in (300, 500, 1000):
            expected = {name: row['sequence'].upper() for (name, requested), row in primary.items() if requested == length}
            exported = read_fasta(paths[f'promoter_{length}'])
            check(f'promoter_FASTA_exact:{length}', exported == expected and set(exported) == family)
        for name in sorted(family):
            windows = [primary.get((name, length)) for length in (300, 500, 1000)]
            check('primary_windows_nested:' + name, all(windows) and
                  all(right['sequence'].upper().endswith(left['sequence'].upper()) for left, right in zip(windows, windows[1:])) and
                  len({row['tss1'] for row in windows}) == 1 and
                  all(int(row['tss1']) == int(by_id[name]['predicted_tss1']) for row in windows))
        bed = out / 'independent_intervals.bed'
        with bed.open('w') as handle:
            for name, job in sorted(jobs.items(), key=lambda item: (item[1]['seqid'], item[1]['start1'], item[0])):
                handle.write(f'{job["seqid"]}\t{job["start1"] - 1}\t{job["end1"]}\t{name}\t0\t{job["strand"]}\n')
        version = subprocess.run([executable, 'version'], text=True, capture_output=True, check=True)
        summary['seqkit_version'] = version.stdout.strip()
        independent = {}
        if jobs:
            extracted_path, stderr_path = out / 'independent_seqkit.fa', out / 'seqkit.stderr.log'
            argv = [executable, 'subseq', '--bed', str(bed), '-w', '0', '-j', '1', '-t', 'dna', '-']
            command = dict(argv=argv, stdin_source=str(genome), stdout=str(extracted_path), stderr=str(stderr_path), started_utc=utc())
            summary['commands'].append(command)
            save()
            with genome.open('rb') as stdin, extracted_path.open('wb') as stdout, stderr_path.open('wb') as stderr:
                completed = subprocess.run(argv, stdin=stdin, stdout=stdout, stderr=stderr, cwd=out, check=False)
            command.update(exit_code=completed.returncode, finished_utc=utc())
            if completed.returncode:
                raise RuntimeError('SeqKit extraction failed; inspect seqkit.stderr.log')
            independent = read_fasta(extracted_path, set(jobs))
        else:
            summary['sequence_extraction_status'] = 'no_search_loci_or_promoters; empty exports checked; SeqKit subseq not applicable'
        check('SeqKit_extracted_all_job_ids', set(independent) == set(jobs), expected=len(jobs), observed=len(independent))
        comparisons = []
        for name, job in sorted(jobs.items()):
            observed = independent.get(name, '')
            good = observed == job['sequence']
            comparisons.append(dict(job_id=name, expected_length=len(job['sequence']), observed_length=len(observed), expected_sha256=seq_digest(job['sequence']), observed_sha256=seq_digest(observed), passed=good))
            check('independent_sequence:' + name, good)
        write_table(out / 'sequence_comparisons.tsv', comparisons, ['job_id', 'expected_length', 'observed_length', 'expected_sha256', 'observed_sha256', 'passed'])
        summary['counts'] = dict(genome_records=len(lengths), genome_bases=sum(lengths.values()), all_search_loci=len(loci), family_candidates=len(family), promoter_table_rows=len(promoters), primary_promoter_rows=len(primary), explicitly_empty_promoters=empty_count, independent_extractions=len(independent), extraction_strands=dict(Counter(j['strand'] for j in jobs.values())))
    except Exception as error:
        summary['error'] = str(error)
        check('validator_completed_without_exception', False)
    finally:
        preservation = [dict(path=path, before_sha256=sha, after_sha256=digest(path) if Path(path).is_file() else 'MISSING') for path, sha in before.items()]
        for row in preservation:
            row['unchanged'] = row['before_sha256'] == row['after_sha256']
        write_table(out / 'source_preservation.tsv', preservation, ['path', 'before_sha256', 'after_sha256', 'unchanged'])
        check('all_read_sources_unchanged', all(row['unchanged'] for row in preservation))
        for path, old_hash in index_before.items():
            new_hash = digest(path) if Path(path).is_file() else None
            check('source_index_not_created_or_modified:' + path, new_hash == old_hash)
        summary['status'] = 'PASS' if all(c['passed'] for c in summary['checks']) else 'FAIL'
        summary['finished_utc'], summary['check_count'] = utc(), len(summary['checks'])
        save()
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--seqkit', default='seqkit')
    args = parser.parse_args()
    try:
        result = validate(args)
        print(json.dumps({k: v for k, v in result.items() if k not in {'checks', 'commands'}}, indent=2))
        return 0 if result['status'] == 'PASS' else 1
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print('Independent output QA failed: ' + str(error), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
