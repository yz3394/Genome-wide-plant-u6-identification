"""Focused adversarial tests; fixtures are synthetic and not biological truths."""

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

MODULE_PATH = Path(__file__).parents[1] / "scripts/u6_loci.py"
SPEC = importlib.util.spec_from_file_location("u6_loci", MODULE_PATH)
m = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = m
SPEC.loader.exec_module(m)


def hit(name, start=100, end=204, query="q", method="blastn", strand="+", **kw):
    values = dict(raw_hit_id=name, method=method, query=query, seqid="Chr1", start1=start,
                  end1=end, strand=strand, source="synthetic", line=1)
    if method == "blastn":
        values.update(pident=95, qstart=1, qend=100, qlen=100, qcov_exact=100,
                      qcov_reported=100, bitscore=120, gapopen=0, subject_length=10000)
    else:
        values.update(score=60, ga=44, ga_pass=True, model_coverage=1, trunc="no")
    values.update(kw)
    return m.Hit(**values)


class LocusTests(unittest.TestCase):
    def test_single_hsp_cannot_combine_metrics(self):
        a = hit("a", qstart=1, qend=60, qcov_exact=60, pident=100, bitscore=120)
        b = hit("b", qstart=41, qend=100, qcov_exact=60, pident=90, bitscore=110)
        rf = hit("rf", query="RF00026", method="cmsearch")
        result = m.summarize_cluster("test", {"anchor": rf, "hits": [rf, a, b], "atac": []})
        self.assertEqual(result["n_strong_queries"], 0)
        self.assertEqual(result["legacy_grade"], "A2_CM_led_U6")

    def test_rounded_coverage_cannot_cross_cutoff(self):
        h = hit("h", qlen=102, qend=81, qcov_reported=80, qcov_exact=100*81/102)
        self.assertNotEqual(m.single_hsp_support(h), "strong")

    def test_fixed_anchor_stops_transitive_bridge(self):
        rf = hit("rf", query="RF00026", method="cmsearch")
        middle = hit("m", 150, 220, bitscore=150)
        end = hit("e", 205, 270, query="second")
        clusters, _ = m.cluster_hits([rf, middle, end], [])
        self.assertEqual(len(clusters), 2)
        self.assertNotIn(end, clusters[0]["hits"])

    def test_nonoverlapping_neighbors_not_merged(self):
        clusters, _ = m.cluster_hits([hit("a", 100, 150), hit("b", 151, 204, query="other")], [])
        self.assertEqual(len(clusters), 2)

    def test_duplicate_does_not_discard_strong_in_favor_of_score(self):
        strong = hit("strong", pident=81, bitscore=80)
        weak = hit("weak", pident=79, bitscore=100)
        m.suppress_duplicates([strong, weak])
        self.assertEqual(weak.suppressed_by, "strong")
        self.assertEqual(strong.suppressed_by, "")

    def test_competition_requires_same_strand_and_material_overlap(self):
        rf = hit("rf", query="RF00026", method="cmsearch")
        opposite = hit("opp", query="RF00619", method="cmsearch", strand="-")
        edge = hit("edge", 195, 300, query="RF00619", method="cmsearch")
        clusters, standalone = m.cluster_hits([rf], [opposite, edge])
        self.assertEqual(clusters[0]["atac"], [])
        self.assertEqual(len(standalone), 2)

    def test_both_ga_retained_as_unresolved_not_deleted(self):
        rf = hit("rf", query="RF00026", method="cmsearch")
        atac = hit("atac", 100, 225, query="RF00619", method="cmsearch")
        clusters, _ = m.cluster_hits([rf], [atac])
        result = m.summarize_cluster("test", clusters[0])
        self.assertEqual(result["family_class"], "U6_U6atac_ambiguous")
        self.assertFalse(result["canonical_retained"])

    def test_weak_atac_does_not_override_rf_ga(self):
        rf = hit("rf", query="RF00026", method="cmsearch")
        atac = hit("atac", query="RF00619", method="cmsearch", score=24, ga=33, ga_pass=False)
        result = m.summarize_cluster("test", {"anchor": rf, "hits": [rf], "atac": [atac]})
        self.assertEqual(result["family_class"], "RF00026_supported")

    def test_reverse_q1_extrapolation_and_gapped_rejection(self):
        h = hit("h", query="X52528.1", strand="-", qstart=5, qend=100, qcov_exact=96)
        result = m.summarize_cluster("test", {"anchor": h, "hits": [h], "atac": []})
        self.assertEqual(result["best_q1_anchor"], 208)
        self.assertFalse(result["q1_anchor_safe"])
        h.gapopen = 1
        result = m.summarize_cluster("test", {"anchor": h, "hits": [h], "atac": []})
        self.assertIsNone(result["best_q1_anchor"])

    def test_nonhistorical_query_not_silently_tss_anchor(self):
        h = hit("h", query="different_mature_U6")
        result = m.summarize_cluster("test", {"anchor": h, "hits": [h], "atac": []})
        self.assertIsNone(result["best_q1_anchor"])


class AlignmentTests(unittest.TestCase):
    def parse(self, content):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "test.sto"
            path.write_text(content)
            return m.parse_stockholm_metrics(path)

    def test_interleaved_stockholm_internal_deletion(self):
        metrics = self.parse("# STOCKHOLM 1.0\nseq A--CG\n#=GR seq PP 9..99\n#=GC RF xxxxx\n#=GC SS_cons <....\n\nseq tT\n#=GC RF .x\n#=GC SS_cons .>\n//\n")["seq"]
        self.assertEqual(metrics["rf_columns"], 6)
        self.assertEqual(metrics["occupied_rf_columns"], 4)
        self.assertEqual(metrics["longest_internal_deletion"], 2)
        self.assertEqual(metrics["inserted_nt"], 1)
        self.assertEqual(metrics["paired_compatible"], 1)

    def test_73nt_full_model_span_not_called_complete(self):
        metric = dict(rf_occupancy=73/105, internal_missing_columns=32,
                      longest_internal_deletion=25, terminal_missing_5p=0,
                      terminal_missing_3p=0, aligned_sequence_length=73)
        result = m.classify_integrity("A"*73, metric, cm_coverage=1, cm_trunc="no")
        self.assertEqual(result["body_integrity"], "structural_deletion_risk")
        self.assertEqual(result["functional_status"], "not_tested")

    def test_cm_coverage_alone_leaves_integrity_unresolved(self):
        result = m.classify_integrity("A"*105, None, cm_coverage=1, cm_trunc="no")
        self.assertEqual(result["body_integrity"], "unresolved")

    def test_pseudogene_annotation_independent_of_alignment(self):
        metric = dict(rf_occupancy=1, internal_missing_columns=0,
                      longest_internal_deletion=0, terminal_missing_5p=0,
                      terminal_missing_3p=0, aligned_sequence_length=105)
        result = m.classify_integrity("A"*105, metric, annotated_pseudogene=True)
        self.assertEqual(result["body_integrity"], "alignment_near_complete_not_functionally_validated")
        self.assertIn("annotated_pseudogene_requires_review", result["integrity_reasons"])

    def test_missing_rf_raises_not_false_complete(self):
        with self.assertRaises(ValueError):
            self.parse("# STOCKHOLM 1.0\nseq ACGT\n//\n")

    def test_partial_alignment_rejected(self):
        with self.assertRaises(ValueError):
            self.parse("# STOCKHOLM 1.0\nseq ACGT\n#=GC RF xxxx\n")

    def test_equal_length_wrong_sequence_cannot_be_used(self):
        metrics = self.parse("# STOCKHOLM 1.0\nseq ACGT\n#=GC RF xxxx\n//\n")["seq"]
        result = m.classify_integrity("TGCA", metrics)
        self.assertIn("alignment_sequence_checksum_mismatch", result["integrity_reasons"])

    def test_model_metadata_reads_frozen_header(self):
        skill = Path(__file__).parents[1]
        model = skill / "resources/models/RF00026_U6.cm"
        meta = m.read_cm_metadata(model)
        self.assertEqual((meta["accession"], meta["clen"], meta["ga"]), ("RF00026", 105, 44.0))


class FrozenEvidenceTests(unittest.TestCase):
    def test_extended_blast_retains_alignment_and_exact_coverage(self):
        row = "q\tChr1\t100\t4\t0\t0\t1\t4\t100\t103\t1e-4\t80\t5\t1000\t80\tACGT\tACGT\t4\n"
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "blast.tsv"
            path.write_text(row)
            hits = m.parse_blast("test", path)
        self.assertEqual((hits[0].qseq, hits[0].sseq, hits[0].btop), ("ACGT", "ACGT", "4"))
        self.assertEqual(hits[0].qcov_exact, 80.0)



if __name__ == "__main__":
    unittest.main()

