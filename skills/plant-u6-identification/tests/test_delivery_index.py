"""Generated delivery links must resolve to the actual output locations."""
from pathlib import Path
import re
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_u6 import write_reports


class DeliveryIndexTests(unittest.TestCase):
    def test_generated_manifest_links_resolve(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            (out / '04_delivery').mkdir()
            (out / '02_curation').mkdir()
            (out / '02_curation/cmalign_version.txt').write_text('test version\n')
            (out / '04_delivery/input_manifest_with_md5.tsv').write_text('path\tmd5\n')
            (out / 'output_manifest.tsv').write_text('path\tsha256\n')
            config = dict(key='demo', species='Test reference', assembly='test-v1', reference_scope='diagnostic', cm_min_score=20)
            counts = dict(all_search_loci=0, RF00026_supported=0, strict_AtU6=0, integrity={}, generic_geometry={})
            summary = dict(counts=counts, discovery=dict(standalone_atac_GA_count=0, family_counts={}))
            ledger = dict(commands=[])
            search = dict(inputs=dict(genome=dict(fasta_records=1)), commands=[])
            write_reports(out, config, summary, ledger, search, {})
            targets = re.findall(r'\]\(([^)]+)\)', (out / 'README_where_to_find.md').read_text())
            self.assertIn('output_manifest.tsv', targets)
            self.assertIn('04_delivery/input_manifest_with_md5.tsv', targets)
            self.assertIn('04_delivery/commands_and_versions.txt', targets)
            self.assertTrue(all((out / target).is_file() for target in targets))


if __name__ == '__main__':
    unittest.main()
