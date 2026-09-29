#!/usr/bin/env python3
"""Check real fixture cases after a complete portable run; no claimed sensitivity.

Run after run_u6.py using tests/fixtures/mini_config.json. These are known
development examples, not a blinded biological validation set. Also run the
independent validate_outputs.py against the actual exported sequences.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from u6_io import read_fasta
from u6_promoters import scan_motifs


def check(run, out):
    if out.exists() or out.is_symlink():
        raise FileExistsError('Choose a fresh fixture-check directory')
    delivery = run / '04_delivery'
    fixtures = Path(__file__).parent / 'fixtures'
    cases = json.loads((fixtures / 'known_cases.json').read_text())
    with (delivery / 'fixture_U6_candidates.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    bodies = read_fasta(delivery / 'fixture_U6_body.fa')
    promoters = {length: read_fasta(delivery / f'fixture_U6_promoter_{length}.fa') for length in (300, 500, 1000)}
    published = read_fasta(fixtures / 'published_promoter_blocks.fa')
    checks = []

    def expect(name, ok, detail=None):
        checks.append(dict(check=name, passed=bool(ok), detail=detail))

    by_case = {}
    for case in cases['derivative_intervals']:
        overlaps = [r for r in rows if r['seqid'] == case['fixture_id'] and r['strand'] == case['strand']
                    and min(int(r['body_end1']), case['local_body_end1']) >= max(int(r['body_start1']), case['local_body_start1'])]
        if case['family'] == 'U6':
            expect(case['fixture_id'] + ':family_retained', len(overlaps) == 1, [r['locus_id'] for r in overlaps])
            if len(overlaps) != 1:
                continue
            row = overlaps[0]
            by_case[case['fixture_id']] = row
            name = row['locus_id']
            expect(case['fixture_id'] + ':exact_source_body', hashlib.sha256(bodies[name].encode()).hexdigest() == case['expected_body_sha256'])
            expect(case['fixture_id'] + ':coordinates', (int(row['body_start1']), int(row['body_end1'])) == (case['local_body_start1'], case['local_body_end1']))
            expect(case['fixture_id'] + ':nested_promoters', promoters[1000][name].endswith(promoters[500][name]) and promoters[500][name].endswith(promoters[300][name]))
            expect(case['fixture_id'] + ':not_new_functional_claim', row['functional_status'] == 'not_tested')
        else:
            expect(case['fixture_id'] + ':not_canonical_U6', not overlaps)
    if 'rice_OsU6a' in by_case:
        row = by_case['rice_OsU6a']
        name = row['locus_id']
        expect('OsU6a:source_block_441plusG', promoters[500][name].endswith(published['OsU6a_Kitaake_Kim2019_FigS1'][:-1]))
        motif = scan_motifs(promoters[500][name], -len(promoters[500][name]), 'generic_provisional')['best']
        expect('OsU6a:both_motifs_detected', (motif.get('use_sequence'), motif.get('tata_sequence')) == ('GTACCACCTCG', 'CTTATATG'))
        expect('OsU6a:noncanonical_geometry_preserved', motif.get('use_tata_edge_distance') == 20 and motif.get('motif_pair_class') == 'provisional_uncertain_geometry')
    for fixture_id, source in [('grape_VvU6_1', 'VvU6.1'), ('grape_VvU6_2', 'VvU6.2')]:
        if fixture_id in by_case:
            row = by_case[fixture_id]
            expect(source + ':published_promoter_suffix', promoters[1000][row['locus_id']].endswith(published[source]))
    if 'rice_CM_rescue' in by_case:
        row = by_case['rice_CM_rescue']
        expect('OsU6L0022:CM_rescue_not_strict_At', row['strict_historical_query_pass'] == 'False')
    activity_mapping = json.loads((run / 'qa/activity_evidence_mapping.json').read_text())
    if activity_mapping['provided']:
        expect('source_activity:three_exact_matches', len(activity_mapping['provided']) == len(activity_mapping['matched']) == 3)
        for case_name in ('rice_OsU6a', 'grape_VvU6_1', 'grape_VvU6_2'):
            row = by_case[case_name]
            expect(case_name + ':source_activity_preserved', row['activity_evidence'] == 'provided_source_support')
        row = by_case['rice_OsU6a']
        expect('OsU6a:source_supported_noncanonical_not_P1',
               row['activity_geometry_interpretation'] == 'experiment_supported_noncanonical_geometry'
               and row['generic_provisional_class'] == 'provisional_uncertain_geometry'
               and row['legacy_arabidopsis_class'] == 'P3_uncertain')
        expect('OsU6a:source_caveat_retained', all('not native expression' in e['note'] for e in activity_mapping['matched']))
    out.mkdir(parents=True)
    summary = dict(status='PASS' if all(c['passed'] for c in checks) else 'FAIL', checks_run=len(checks),
                   checks=checks, candidate_count=len(rows), run=str(run.resolve()),
                   scope='Known compact-regression cases only; not whole-genome or independent biological validation')
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'checks'}))
    return 0 if summary['status'] == 'PASS' else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    raise SystemExit(check(args.run, args.out))
