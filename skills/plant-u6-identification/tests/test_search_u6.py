"""No genome search is performed by these planning and guard tests."""

import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import search_u6


class SearchPlanningTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        (self.base / "genome.fa").write_text(">scaffold1\nACGT\n>scaffold2\nACGT\n")
        (self.base / "queries.fa").write_text(">q\nACGT\n")
        for accession, length, ga in [("RF00026", 105, 43.5), ("RF00619", 126, 31.25)]:
            (self.base / f"{accession}.cm").write_text(
                f"INFERNAL1/a\nNAME toy\nACC {accession}\nCLEN {length}\nW 199\nGA {ga}\nCM\n//\n")
        self.argv = [
            "--genome", str(self.base / "genome.fa"), "--queries", str(self.base / "queries.fa"),
            "--u6-model", str(self.base / "RF00026.cm"), "--atac-model", str(self.base / "RF00619.cm"),
            "--out", str(self.base / "search"), "--blastn", sys.executable, "--cmsearch", sys.executable,
        ]

    def plan(self, extras=()):
        return search_u6.build_plan(search_u6.argument_parser().parse_args(self.argv + list(extras)))

    def test_independent_three_searches_and_explicit_historical_scoring(self):
        plan = self.plan()
        self.assertEqual([c["name"] for c in plan["commands"]], ["blastn", "RF00026", "RF00619"])
        self.assertEqual(plan["blast_scoring"], {"reward": 1, "penalty": -3, "gapopen": 5, "gapextend": 2, "word_size": 7})
        self.assertEqual(len(plan["blast_output_fields"]), 18)
        self.assertEqual(plan["blast_output_fields"][-3:], ["qseq", "sseq", "btop"])
        self.assertEqual(plan["inputs"]["genome"]["fasta_records"], 2)
        self.assertEqual(plan["models"]["RF00026"]["ga"], 43.5)
        self.assertEqual(plan["models"]["RF00619"]["ga"], 31.25)
        for command in plan["commands"][1:]:
            self.assertEqual(command["argv"][-1], str((self.base / "genome.fa").resolve()))
            self.assertNotIn("--noali", command["argv"])
            self.assertIn("-A", command["argv"])
            self.assertIn("--incT", command["argv"])
        self.assertFalse((self.base / "search").exists())

    def test_dry_run_prints_plan_without_writing_output(self):
        stdout = io.StringIO()
        with patch("search_u6.capture_tool_versions", return_value={"mock_only": True}), contextlib.redirect_stdout(stdout):
            self.assertEqual(search_u6.main(self.argv + ["--dry-run"]), 0)
        result = json.loads(stdout.getvalue())
        self.assertTrue(result["dry_run"])
        self.assertEqual(result["status"], "planned")
        self.assertFalse((self.base / "search").exists())

    def test_overwrite_guard_including_dry_run(self):
        (self.base / "search").mkdir()
        sentinel = self.base / "search" / "keep.txt"
        sentinel.write_text("original")
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.plan(["--dry-run"])
        self.assertEqual(sentinel.read_text(), "original")

    def test_dangling_output_symlink_rejected(self):
        (self.base / "search").symlink_to(self.base / "missing")
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.plan()

    def test_missing_tool_is_explicit_and_no_install_attempted(self):
        with self.assertRaisesRegex(ValueError, "no installation is attempted"):
            self.plan(["--cmsearch", "absolutely_nonexistent_u6_tool_20260928"])

    def test_wrong_or_multiple_models_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Expected RF00026"):
            self.plan(["--u6-model", str(self.base / "RF00619.cm")])
        source = self.base / "RF00026.cm"
        source.write_text(source.read_text() * 2)
        with self.assertRaisesRegex(ValueError, "one Infernal CM"):
            self.plan()

    def test_retention_threshold_must_not_hide_ga_hits(self):
        with self.assertRaisesRegex(ValueError, "hide GA-supported hits"):
            self.plan(["--cm-min-score", "40"])

    def test_other_task_is_not_silently_supported(self):
        args = search_u6.argument_parser().parse_args(self.argv)
        args.task = "blastn"
        with self.assertRaisesRegex(ValueError, "Only the frozen"):
            search_u6.build_plan(args)


if __name__ == "__main__":
    unittest.main()

