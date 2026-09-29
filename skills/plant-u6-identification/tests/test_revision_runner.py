"""Pure-function and isolated-fixture checks for the revised runner/I/O."""

import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
from run_u6 import map_query_position, run
from u6_io import IndexedFasta, digest, revcomp, write_tsv


class AlignmentCoordinateTests(unittest.TestCase):
    def test_query_one_on_match(self):
        self.assertEqual(map_query_position("ACGT", "ACGT", 1, 61, 1), 61)

    def test_query_one_after_subject_insertion(self):
        self.assertEqual(map_query_position("--ACGT", "TTACGT", 1, 59, 1), 61)

    def test_deleted_query_residue_has_no_subject_coordinate(self):
        self.assertIsNone(map_query_position("ACGT", "-CGT", 1, 61, 1))
        self.assertEqual(map_query_position("ACGT", "-CGT", 1, 61, 2), 61)

    def test_missing_query_one_is_not_extrapolated(self):
        self.assertIsNone(map_query_position("CGT", "CGT", 2, 61, 1))

    def test_negative_genomic_projection_from_oriented_window(self):
        oriented_pos = map_query_position("ACGT", "ACGT", 1, 61, 1)
        window_high = 300
        self.assertEqual(window_high - oriented_pos + 1, 240)

    def test_column_count_mismatch_rejected(self):
        with self.assertRaises(ValueError):
            map_query_position("ACGT", "ACG", 1, 61, 1)


class BaselineAndPreservationTests(unittest.TestCase):
    def test_existing_output_refused_before_missing_config(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            output = root / "existing"
            output.mkdir()
            with self.assertRaisesRegex(ValueError, 'Output exists'):
                run(root / "missing.json", output)

    def test_dangling_output_link_refused(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            output = root / "existing"
            output.symlink_to(root / "absent")
            with self.assertRaisesRegex(ValueError, 'Output exists'):
                run(root / "missing.json", output)

    def test_existing_output_refused_before_other_input_access(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            config = root / "config.json"
            config.write_text("{}\n")
            output = root / "existing"
            output.mkdir()
            evidence = output / "keep.txt"
            evidence.write_text("Do not replace this pre-existing result.\n")
            before = digest(evidence)
            with self.assertRaisesRegex(ValueError, 'Output exists'):
                run(config, output)
            self.assertEqual(digest(evidence), before)
            self.assertEqual(list(output.iterdir()), [evidence])


class FastaAndTableTests(unittest.TestCase):
    def test_indexed_read_crosses_crlf_lines_and_preserves_source(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "genome.fa"
            path.write_bytes(b">c\r\nACGT\r\nATNN\r\n")
            index = Path(str(path) + ".fai")
            index.write_text("c\t8\t4\t4\t6\n")
            before = (digest(path), digest(index))
            reader = IndexedFasta(path)
            try:
                self.assertEqual(reader.fetch("c", 3, 7), "GTATN")
                with self.assertRaises(ValueError):
                    reader.fetch("c", 0, 3)
            finally:
                reader.close()
            self.assertEqual((digest(path), digest(index)), before)

    def test_all_standard_iupac_bases_supported(self):
        sequence = "ACGTRYSWKMBDHVN"
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "genome.fa"
            path.write_text(">c\n" + sequence + "\n")
            Path(str(path) + ".fai").write_text(f"c\t{len(sequence)}\t3\t{len(sequence)}\t{len(sequence)+1}\n")
            reader = IndexedFasta(path)
            try:
                self.assertEqual(reader.fetch("c", 1, len(sequence)), sequence)
            finally:
                reader.close()
        self.assertEqual(revcomp(revcomp(sequence)), sequence)

    def test_tsv_serializes_sets_deterministically(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "table.tsv"
            write_tsv(path, [{"locus": "one", "ids": {"b", "a"}}])
            with path.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(json.loads(rows[0]["ids"]), ["a", "b"])


if __name__ == "__main__":
    unittest.main()
