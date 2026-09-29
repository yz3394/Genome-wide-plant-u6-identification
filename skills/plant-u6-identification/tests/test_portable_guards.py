"""Audit and recovery checks using a frozen, actually executed tiny search."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_u6 import audit_fasta, audit_gff, load_config, validate_reused_search
from u6_io import digest

FIXTURES = Path(__file__).parent / 'fixtures'


class PortableAuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def genome(self, text):
        path = self.root / 'genome.fa'
        path.write_text(text)
        return path

    def gff(self, text):
        path = self.root / 'models.gff3'
        path.write_text(text)
        return path

    def test_duplicate_fasta_ID_rejected_without_source_index_write(self):
        path = self.genome('>same\nACGT\n>same\nACGT\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate FASTA ID'):
            audit_fasta(path, self.root / 'new.fai')
        self.assertFalse(Path(str(path) + '.fai').exists())

    def test_audited_index_is_separate_and_source_hash_preserved(self):
        path = self.genome('>one\r\nACGT\r\nNN\r\n>two\r\nTGCA\r\n')
        before = digest(path)
        result = audit_fasta(path, self.root / 'new.fai')
        self.assertEqual(result['lengths'], {'one': 6, 'two': 4})
        self.assertEqual(result['N_bases'], 2)
        self.assertEqual(digest(path), before)
        self.assertFalse(Path(str(path) + '.fai').exists())

    def test_existing_generated_index_not_overwritten(self):
        index = self.root / 'new.fai'
        index.write_text('old index\n')
        with self.assertRaisesRegex(ValueError, 'Index output exists'):
            audit_fasta(self.genome('>one\nACGT\n'), index)
        self.assertEqual(index.read_text(), 'old index\n')

    def test_empty_record_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Empty FASTA record'):
            audit_fasta(self.genome('>empty\n'), self.root / 'new.fai')

    def test_nonuniform_wrapping_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Nonuniform FASTA wrapping'):
            audit_fasta(self.genome('>one\nACGT\nAC\nACGT\n'), self.root / 'new.fai')

    def test_invalid_DNA_rejected(self):
        with self.assertRaisesRegex(ValueError, 'non-DNA'):
            audit_fasta(self.genome('>protein\nPEPTIDE\n'), self.root / 'new.fai')

    def test_GFF_mismatched_record_fails(self):
        path = self.gff('wrong\tx\tgene\t1\t4\t.\t+\t.\tID=g\n')
        with self.assertRaisesRegex(ValueError, 'GFF/FASTA seqid mismatch'):
            audit_gff(path, {'right': 10})

    def test_GFF_out_of_bounds_fails(self):
        path = self.gff('right\tx\tgene\t1\t11\t.\t+\t.\tID=g\n')
        with self.assertRaisesRegex(ValueError, 'out of bounds'):
            audit_gff(path, {'right': 10})

    def test_GFF_sequence_region_mismatch_fails(self):
        path = self.gff('##gff-version 3\n##sequence-region right 1 11\n')
        with self.assertRaisesRegex(ValueError, 'sequence-region out of bounds'):
            audit_gff(path, {'right': 10})

    def test_valid_partial_GFF_sequence_region_is_not_rejected_or_full_coverage(self):
        path = self.gff('##gff-version 3\n##sequence-region right 20 80\n')
        result = audit_gff(path, {'right': 100})
        self.assertEqual(result['status'], 'coordinate_checks_passed')
        self.assertEqual(result['declared_sequence_regions'], [dict(seqid='right', start1=20, end1=80, covers_whole_record=False)])

    def test_GFF_explicit_exclusion_cannot_hide_included_records(self):
        path = self.gff('right\tx\tgene\t1\t10\t.\t+\t.\tID=g\n')
        with self.assertRaisesRegex(ValueError, 'may not hide'):
            audit_gff(path, {'right': 10}, ['right'])

    def test_GFF_explicit_nonincluded_exclusion_is_recorded(self):
        path = self.gff('plastid\tx\tgene\t1\t99\t.\t+\t.\tID=g\n')
        result = audit_gff(path, {'nuclear': 10}, ['plastid'])
        self.assertEqual(result['excluded_features'], {'plastid': 1})
        self.assertEqual(result['included_genome_records_without_features'], ['nuclear'])

    def test_no_GFF_is_unassessed_not_clean(self):
        result = audit_gff(None, {'one': 10})
        self.assertEqual(result['status'], 'not_provided')
        self.assertIsNone(result['included_features'])

    def config(self):
        value = json.loads((FIXTURES / 'mini_config.json').read_text())
        for key in ('genome', 'gff', 'query_panel', 'historical_query', 'u6_model', 'atac_model'):
            value[key] = str((FIXTURES / value[key]).resolve())
        for key in ('blastn', 'cmsearch', 'cmalign'):
            value[key] = sys.executable
        return value

    def load(self, value):
        path = self.root / 'config.json'
        path.write_text(json.dumps(value))
        return load_config(path)

    def test_missing_dependency_fails_without_install(self):
        value = self.config()
        value['cmsearch'] = 'plant_u6_no_such_dependency_20260928'
        with self.assertRaisesRegex(ValueError, 'no installation is attempted'):
            self.load(value)

    def test_nonexperimental_activity_class_rejected(self):
        value = self.config()
        value['activity_evidence'] = [dict(seqid='x', strand='+', body_start1=1, body_end1=102,
                                           source='prediction', evidence_type='motif_prediction')]
        with self.assertRaisesRegex(ValueError, 'must describe an experiment'):
            self.load(value)


class RealSearchRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.search = self.root / 'completed_search'
        shutil.copytree(FIXTURES / 'completed_search', self.search)
        self.ledger_path = self.search / 'search_run.json'
        self.ledger = json.loads(self.ledger_path.read_text())
        self.config = json.loads((FIXTURES / 'mini_config.json').read_text())
        for key in ('genome', 'gff', 'query_panel', 'historical_query', 'u6_model', 'atac_model'):
            self.config[key] = str((FIXTURES / self.config[key]).resolve())
        self.config.update(blastn=sys.executable, cmsearch=sys.executable, cmalign=sys.executable, cm_min_score=20.0)

    def mutate(self):
        self.ledger_path.write_text(json.dumps(self.ledger))

    def test_actual_complete_matching_search_reusable_without_modification(self):
        before = {p.name: digest(p) for p in self.search.iterdir() if p.is_file()}
        result = validate_reused_search(self.search, self.config)
        self.assertEqual(result['status'], 'completed')
        self.assertEqual(before, {p.name: digest(p) for p in self.search.iterdir() if p.is_file()})

    def test_running_or_failed_search_not_complete(self):
        for status in ('running', 'failed'):
            self.ledger['status'] = status
            self.mutate()
            with self.assertRaisesRegex(ValueError, 'completed non-dry-run'):
                validate_reused_search(self.search, self.config)

    def test_dry_run_cannot_be_reused(self):
        self.ledger['dry_run'] = True
        self.mutate()
        with self.assertRaisesRegex(ValueError, 'completed non-dry-run'):
            validate_reused_search(self.search, self.config)

    def test_fake_completed_marker_not_sufficient(self):
        self.ledger_path.write_text('{"status":"completed"}')
        with self.assertRaises(ValueError):
            validate_reused_search(self.search, self.config)

    def test_changed_input_hash_rejected(self):
        genome = self.root / 'changed.fa'
        genome.write_text(Path(self.config['genome']).read_text() + '\n')
        self.config['genome'] = str(genome)
        with self.assertRaisesRegex(ValueError, 'input mismatch: genome'):
            validate_reused_search(self.search, self.config)

    def test_missing_raw_output_rejected(self):
        (self.search / 'RF00619.tblout').unlink()
        with self.assertRaisesRegex(ValueError, 'absent or changed'):
            validate_reused_search(self.search, self.config)

    def test_changed_raw_output_rejected(self):
        path = self.search / 'RF00026.tblout'
        path.write_text(path.read_text() + '# unexpected change\n')
        with self.assertRaisesRegex(ValueError, 'absent or changed'):
            validate_reused_search(self.search, self.config)

    def test_incomplete_output_manifest_rejected(self):
        self.ledger['outputs'] = []
        self.mutate()
        with self.assertRaisesRegex(ValueError, 'manifest is incomplete'):
            validate_reused_search(self.search, self.config)

    def test_failed_command_under_completed_marker_rejected(self):
        self.ledger['commands'][1]['exit_code'] = 1
        self.mutate()
        with self.assertRaisesRegex(ValueError, 'unsuccessful or incomplete'):
            validate_reused_search(self.search, self.config)

    def test_unfinished_command_under_completed_marker_rejected(self):
        self.ledger['commands'][2].pop('finished_utc', None)
        self.mutate()
        with self.assertRaisesRegex(ValueError, 'unsuccessful or incomplete'):
            validate_reused_search(self.search, self.config)

    def test_extra_strand_filter_cannot_hide_in_matching_scoring_dict(self):
        self.ledger['commands'][1]['argv'].insert(1, '--toponly')
        self.mutate()
        with self.assertRaisesRegex(ValueError, 'full argv'):
            validate_reused_search(self.search, self.config)

    def test_extra_BLAST_hsp_limit_cannot_hide_in_matching_scoring_dict(self):
        self.ledger['commands'][0]['argv'].extend(['-max_hsps', '1'])
        self.mutate()
        with self.assertRaisesRegex(ValueError, 'full argv'):
            validate_reused_search(self.search, self.config)

    def test_extra_bottom_strand_filter_rejected(self):
        self.ledger['commands'][1]['argv'].insert(1, '--bottomonly')
        self.mutate()
        with self.assertRaisesRegex(ValueError, 'full argv'):
            validate_reused_search(self.search, self.config)

    def test_cutoff_above_GA_rejected_even_if_reporting_metadata_matches(self):
        self.config['cm_min_score'] = 50
        self.ledger['cm_reporting_score'] = 50
        self.mutate()
        with self.assertRaisesRegex(ValueError, 'hide GA-supported hits'):
            validate_reused_search(self.search, self.config)

    def test_changed_scoring_profile_not_reused(self):
        self.ledger['blast_scoring']['reward'] = 2
        self.mutate()
        with self.assertRaisesRegex(ValueError, 'scoring/field schema mismatch'):
            validate_reused_search(self.search, self.config)

    def test_changed_reporting_threshold_not_reused(self):
        self.config['cm_min_score'] = 21
        with self.assertRaisesRegex(ValueError, 'reporting parameter mismatch'):
            validate_reused_search(self.search, self.config)

    def test_changed_thread_parameter_requires_new_search_by_contract(self):
        self.config['threads'] = 3
        with self.assertRaisesRegex(ValueError, 'full argv'):
            validate_reused_search(self.search, self.config)


if __name__ == '__main__':
    unittest.main()
