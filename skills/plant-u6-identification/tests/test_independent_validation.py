"""Real SeqKit checks on synthetic plus/minus/edge and corruption fixtures."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import random
import shutil
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'validate_outputs.py'
SPEC = importlib.util.spec_from_file_location('independent_output_validator', SCRIPT)
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)


def table(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def fasta(path, records):
    with path.open('w') as handle:
        for name, sequence in records.items():
            handle.write('>' + name + '\n' + sequence + '\n')


def seqsha(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def fixture(run, empty=False):
    """Fixture writer deliberately does not import any primary pipeline code."""
    run.mkdir()
    for name in ('00_audit', '02_curation', '03_promoters', '04_delivery', 'qa'):
        (run / name).mkdir()
    rng = random.Random(60821)
    sequence = ''.join(rng.choice('ACGT') for _ in range(2500))
    sequence = sequence[:750] + 'N' + sequence[751:]
    genome = run / 'genome.fa'
    fasta(genome, {'chr_test': sequence})
    (run / '00_audit/config.resolved.json').write_text(json.dumps({'key': 'test', 'genome': str(genome.resolve())}))
    (run / 'run.json').write_text(json.dumps({'status': 'completed', 'internal_qa': 'PASS'}))
    loci, bodies, promoters, jobs = [], {}, [], []
    targets = [] if empty else [('plus', 1600, '+'), ('minus', 900, '-'), ('edge', 1, '+'), ('clipped', 2300, '-')]
    for name, tss, strand in targets:
        start, end = (tss, tss + 102) if strand == '+' else (tss - 102, tss)
        body = sequence[start - 1:end]
        if strand == '-':
            body = reverse_complement(body)
        bodies[name] = body
        loci.append(dict(locus_id=name, seqid='chr_test', strand=strand,
                         body_start1=start, body_end1=end, body_length=103,
                         body_sequence_sha256=seqsha(body), canonical_retained=True,
                         rf00026_ga_pass=True, rf00619_ga_pass=False, predicted_tss1=tss))
        jobs.append(dict(job_id='B__' + name, seqid='chr_test', strand=strand,
                         start1=start, end1=end, sequence=body))
        for requested in (300, 500, 1000):
            lo, hi = (max(1, tss - requested), tss - 1) if strand == '+' else (tss + 1, min(2500, tss + requested))
            promoter = sequence[lo - 1:hi] if hi >= lo else ''
            if strand == '-':
                promoter = reverse_complement(promoter)
            status = ('complete' if len(promoter) == requested else 'truncated_reference_edge') if promoter else 'no_upstream_sequence_at_reference_edge'
            promoters.append(dict(locus_id=name, seqid='chr_test', strand=strand, tss1=tss,
                                  hypothesis='fixture', is_primary_hypothesis=True,
                                  requested_length=requested, actual_length=len(promoter),
                                  first_relative_position=-len(promoter) if promoter else 'NA',
                                  start1=lo if promoter else 'NA', end1=hi if promoter else 'NA',
                                  sequence=promoter, truncated=len(promoter) < requested,
                                  N_count=promoter.count('N'), extraction_status=status))
            if promoter:
                jobs.append(dict(job_id=f'P__{name}__fixture__{requested}', seqid='chr_test', strand=strand,
                                 start1=lo, end1=hi, sequence=promoter))
    locus_fields = ['locus_id', 'seqid', 'strand', 'body_start1', 'body_end1', 'body_length', 'body_sequence_sha256', 'canonical_retained', 'rf00026_ga_pass', 'rf00619_ga_pass', 'predicted_tss1']
    promoter_fields = ['locus_id', 'seqid', 'strand', 'tss1', 'hypothesis', 'is_primary_hypothesis', 'requested_length', 'actual_length', 'first_relative_position', 'start1', 'end1', 'sequence', 'truncated', 'N_count', 'extraction_status']
    table(run / '02_curation/all_search_loci.tsv', loci, locus_fields)
    table(run / '04_delivery/test_U6_candidates.tsv', loci, locus_fields)
    table(run / '03_promoters/promoter_sequences.tsv', promoters, promoter_fields)
    table(run / 'qa/extraction_jobs.tsv', jobs, ['job_id', 'seqid', 'strand', 'start1', 'end1', 'sequence'])
    table(run / '04_delivery/input_manifest_with_md5.tsv', [{'path': str(genome.resolve()), 'sha256': VALIDATOR.digest(genome)}])
    fasta(run / '02_curation/all_search_bodies.fa', bodies)
    fasta(run / '04_delivery/test_U6_body.fa', bodies)
    for requested in (300, 500, 1000):
        fasta(run / '04_delivery' / f'test_U6_promoter_{requested}.fa', {p['locus_id']: p['sequence'] for p in promoters if p['requested_length'] == requested})
    with (run / '04_delivery/test_U6_loci.bed').open('w') as handle:
        for row in loci:
            handle.write(f'chr_test\t{row["body_start1"] - 1}\t{row["body_end1"]}\t{row["locus_id"]}\t0\t{row["strand"]}\n')
    return genome


class TestIndependentValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.seqkit = shutil.which(os.environ.get('SEQKIT_BIN', 'seqkit'))
        if not cls.seqkit:
            raise unittest.SkipTest('Real SeqKit unavailable; set SEQKIT_BIN; extraction validation not executed')

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='u6-independent-test-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.run = self.root / 'run'

    def validate(self):
        return VALIDATOR.validate(argparse.Namespace(run=self.run, out=self.root / 'qa', seqkit=self.seqkit))

    def failed(self, result):
        self.assertEqual(result['status'], 'FAIL')
        return [r['name'] for r in result['checks'] if not r['passed']]

    def test_positive_both_strands_clipped_empty_and_input_indexes_untouched(self):
        genome = fixture(self.run)
        index = Path(str(genome) + '.fai')
        index.write_text('not a usable index; validator must not read this\n')
        old = index.read_bytes()
        result = self.validate()
        self.assertEqual(result['status'], 'PASS', result)
        self.assertEqual(result['counts']['family_candidates'], 4)
        self.assertEqual(result['counts']['independent_extractions'], 13)
        self.assertEqual(result['counts']['explicitly_empty_promoters'], 3)
        self.assertEqual(result['counts']['extraction_strands'], {'+': 5, '-': 8})
        self.assertEqual(index.read_bytes(), old)
        self.assertFalse(Path(str(genome) + '.seqkit.fai').exists())
        self.assertTrue(all(row['unchanged'] == 'True' for row in VALIDATOR.read_table(self.root / 'qa/source_preservation.tsv')))

    def test_empty_candidate_run_checked_without_empty_bed_seqkit(self):
        fixture(self.run, empty=True)
        result = self.validate()
        self.assertEqual(result['status'], 'PASS', result)
        self.assertEqual(result['counts']['family_candidates'], 0)
        self.assertEqual(result['commands'], [])
        self.assertIn('not applicable', result['sequence_extraction_status'])

    def test_corrupt_export_detected(self):
        fixture(self.run)
        path = self.run / '04_delivery/test_U6_body.fa'
        records = VALIDATOR.read_fasta(path)
        records['plus'] = ('A' if records['plus'][0] != 'A' else 'C') + records['plus'][1:]
        fasta(path, records)
        self.assertIn('family_body_fasta_exact', self.failed(self.validate()))

    def test_corrupt_driver_jobs_detected(self):
        fixture(self.run)
        path = self.run / 'qa/extraction_jobs.tsv'
        rows = VALIDATOR.read_table(path)
        rows[0]['start1'] = int(rows[0]['start1']) + 1
        table(path, rows)
        self.assertIn('driver_jobs_equal_independently_reconstructed_jobs', self.failed(self.validate()))

    def test_consistent_wrong_body_coordinates_fail_independent_sequence(self):
        fixture(self.run)
        for name in ('02_curation/all_search_loci.tsv', '04_delivery/test_U6_candidates.tsv'):
            path = self.run / name
            rows = VALIDATOR.read_table(path)
            rows[0]['body_start1'] = int(rows[0]['body_start1']) + 1
            rows[0]['body_end1'] = int(rows[0]['body_end1']) + 1
            table(path, rows)
        path = self.run / 'qa/extraction_jobs.tsv'
        rows = VALIDATOR.read_table(path)
        rows[0]['start1'] = int(rows[0]['start1']) + 1
        rows[0]['end1'] = int(rows[0]['end1']) + 1
        table(path, rows)
        path = self.run / '04_delivery/test_U6_loci.bed'
        beds = [line.split('\t') for line in path.read_text().splitlines()]
        beds[0][1], beds[0][2] = str(int(beds[0][1]) + 1), str(int(beds[0][2]) + 1)
        path.write_text('\n'.join('\t'.join(row) for row in beds) + '\n')
        failed = self.failed(self.validate())
        self.assertEqual(failed, ['independent_sequence:B__plus'])

    def test_modified_frozen_genome_detected(self):
        genome = fixture(self.run)
        with genome.open('a') as handle:
            handle.write('\n')
        self.assertTrue(any(n.startswith('frozen_input_sha256:') for n in self.failed(self.validate())))

    def test_discordant_family_boolean_detected(self):
        fixture(self.run)
        for name in ('02_curation/all_search_loci.tsv', '04_delivery/test_U6_candidates.tsv'):
            path = self.run / name
            rows = VALIDATOR.read_table(path)
            rows[0]['rf00619_ga_pass'] = True
            table(path, rows)
        self.assertIn('family_retention_boolean_consistency:plus', self.failed(self.validate()))

    def test_existing_output_refused(self):
        fixture(self.run)
        (self.root / 'qa').mkdir()
        with self.assertRaisesRegex(ValueError, 'fresh QA directory'):
            self.validate()


if __name__ == '__main__':
    unittest.main()
