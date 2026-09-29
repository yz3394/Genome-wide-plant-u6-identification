"""Source-backed regression cases, not independent biological accuracy tests."""
import hashlib
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from u6_io import read_fasta, revcomp
from u6_loci import classify_integrity, parse_blast, single_hsp_support
from u6_promoters import extract_upstream, scan_motifs
from u6_context import GffContext

FIXTURES = Path(__file__).parent / 'fixtures'


class LiteratureRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = json.loads((FIXTURES / 'known_cases.json').read_text())
        cls.records = {r['fixture_id']: r for r in cls.cases['derivative_intervals']}
        cls.genome = read_fasta(FIXTURES / 'mini_genome.fa')
        cls.published = read_fasta(FIXTURES / 'published_promoter_blocks.fa')

    def upstream(self, name, length=1000):
        row = self.records[name]
        tss = row.get('local_predicted_tss1', row['local_body_end1'] if row['strand'] == '-' else row['local_body_start1'])
        return extract_upstream(lambda seqid, start, end: self.genome[seqid][start-1:end],
                                name, row['strand'], tss, length, len(self.genome[name]))

    def test_OsU6a_source_block_contains_annotation_inferred_plus_one(self):
        source = self.published['OsU6a_Kitaake_Kim2019_FigS1']
        self.assertEqual((len(source), source[-1]), (442, 'G'))
        self.assertTrue(self.upstream('rice_OsU6a')['sequence'].endswith(source[:-1]))
        row = self.records['rice_OsU6a']
        plus_one = revcomp(self.genome['rice_OsU6a'][row['local_predicted_tss1']-1:row['local_predicted_tss1']])
        self.assertEqual(plus_one, source[-1])

    def test_OsU6a_detected_motifs_do_not_imply_canonical_geometry(self):
        promoter = self.upstream('rice_OsU6a', 500)
        for profile, expected in [('legacy_arabidopsis', 'P3_uncertain'),
                                  ('generic_provisional', 'provisional_uncertain_geometry')]:
            scan = scan_motifs(promoter['sequence'], promoter['first_relative_position'], profile)
            best = scan['best']
            self.assertEqual((best['use_sequence'], best['tata_sequence']), ('GTACCACCTCG', 'CTTATATG'))
            self.assertEqual((best['use_tata_edge_distance'], best['use_tata_intervening_bases']), (20, 19))
            self.assertEqual((best['use_start_rel_tss'], best['tata_start_rel_tss']), (-60, -30))
            self.assertEqual(best['motif_pair_class'], expected)
        self.assertEqual(self.cases['published_sequence_sources']['OsU6a'], 'https://doi.org/10.1186/s12284-019-0325-7')

    def test_grape_published_U6_sequences_equal_extracted_suffixes(self):
        for name, source in [('grape_VvU6_1', 'VvU6.1'), ('grape_VvU6_2', 'VvU6.2')]:
            self.assertTrue(self.upstream(name)['sequence'].endswith(self.published[source]))
        self.assertEqual((len(self.published['VvU6.1']), len(self.published['VvU6.2'])), (425, 591))

    def test_grape_U6_motifs_have_expected_geometry(self):
        for name in ('grape_VvU6_1', 'grape_VvU6_2'):
            promoter = self.upstream(name, 500)
            best = scan_motifs(promoter['sequence'], -500)['best']
            self.assertEqual(best['motif_pair_class'], 'P1_canonical_geometry')
            self.assertEqual((best['use_tata_edge_distance'], best['use_tata_intervening_bases']), (25, 24))

    def test_grape_U3_also_has_type3_like_motifs(self):
        for name in ('VvU3.1', 'VvU3.2'):
            seq = self.published[name]
            best = scan_motifs(seq, -len(seq), 'generic_provisional')['best']
            self.assertEqual(best['motif_pair_class'], 'provisional_type3_like_geometry')
        self.assertEqual(len(self.published['VvU3.2']), 511)
        # RNA family is NOT assigned by scan_motifs; U3 control label is independently sourced.
        self.assertTrue(self.records['grape_VvU3_1']['body_endpoints_are_synthetic'])

    def test_rice_0022_exact_coverage_changes_with_scoring_not_family_evidence(self):
        by_profile = {}
        for profile, path in [('frozen', 'rice_historical_blast.tsv'), ('alternative', 'rice_alternative_blast.tsv')]:
            hits = parse_blast(profile, FIXTURES / path)
            hist = [h for h in hits if 'X52528.1' in h.query]
            self.assertEqual(len(hist), 1)
            by_profile[profile] = hist[0]
        self.assertAlmostEqual(by_profile['frozen'].qcov_exact, 100*78/102)
        self.assertAlmostEqual(by_profile['alternative'].qcov_exact, 100*83/102)
        self.assertNotEqual(single_hsp_support(by_profile['frozen']), 'strong')
        self.assertLess(by_profile['frozen'].qcov_exact, 80)
        self.assertGreaterEqual(by_profile['alternative'].qcov_exact, 80)
        # No cross-scoring bitscore classification is used as a new biological threshold.

    def test_rounded_80_in_real_Vicia_hit_is_not_exact_80(self):
        hits = parse_blast('frozen', FIXTURES / 'rice_historical_blast.tsv')
        hit = next(h for h in hits if 'X04788.1' in h.query)
        self.assertEqual(hit.qcov_reported, 80)
        self.assertAlmostEqual(hit.qcov_exact, 100*78/98)
        self.assertNotEqual(single_hsp_support(hit), 'strong')

    def test_grape_partial_mapping_does_not_enter_strict_denominator(self):
        records = json.loads((FIXTURES / 'grape_annotation_controls.json').read_text())['records']
        u6 = [r for r in records if r['family'] == 'U6']
        self.assertEqual(len(u6), 15)
        strict = [r for r in u6 if r['seqid'] is not None]
        self.assertEqual(len(strict), 14)
        unresolved = next(r for r in u6 if r['control_id'] == 'ENSRNA049469705')
        self.assertIsNone(unresolved['seqid'])
        partial = self.cases['partial_mapping_hsps'][0]
        self.assertEqual(partial['qualifies'], 'False')
        self.assertEqual(partial['body_endpoint_bases_observed'], 'True')
        self.assertLess(float(partial['qcov_exact']), 99)
        self.assertEqual((partial['mapped_start1'], partial['mapped_end1']), ('3298', '3399'))

    def test_synthetic_short_body_is_not_complete_despite_CM_span(self):
        body = read_fasta(FIXTURES / 'synthetic_rice_U6_short_body.fa')['synthetic_rice_U6_3prime73']
        metrics = self.cases['synthetic_short_body_alignment']
        self.assertEqual(len(body), 73)
        self.assertLess(metrics['occupied_rf_columns'], metrics['rf_columns'])
        self.assertGreater(metrics['terminal_missing_5p'], 0)
        result = classify_integrity(body, metrics, cm_coverage=1, cm_trunc='no')
        self.assertEqual(result['body_integrity'], 'terminally_incomplete_submitted_body')
        self.assertEqual(result['functional_status'], 'not_tested')

    def test_fixture_sequences_match_source_hashes_and_preserve_orientation(self):
        self.assertEqual(len(self.genome), 8)
        for name, sequence in self.genome.items():
            row = self.records[name]
            self.assertEqual(hashlib.sha256(sequence.encode()).hexdigest(), row['sequence_sha256'])
            self.assertEqual(len(sequence), row['source_interval_end1'] - row['source_interval_start1'] + 1)

    def test_synthetic_annotation_self_exclusion_foreign_CDS_and_NA(self):
        context = GffContext(FIXTURES / 'mini_context.gff3', {k: len(v) for k, v in self.genome.items()})
        row = self.records['rice_CM_rescue']
        ids = context.supporting_u6_ids('rice_CM_rescue', row['local_body_start1'], row['local_body_end1'], '+')
        self.assertTrue(ids)
        body = context.audit('rice_CM_rescue', row['local_body_start1'], row['local_body_end1'], '+', self_ids=ids)
        self.assertEqual(body['counts_by_type']['ncRNA'], 0)
        self.assertIsNone(body['counts_by_type']['repeat'])
        upstream = context.audit('rice_CM_rescue', 881, 1000, '+', self_ids=ids)
        self.assertEqual(upstream['counts_by_type']['CDS'], 1)


if __name__ == '__main__':
    unittest.main()
