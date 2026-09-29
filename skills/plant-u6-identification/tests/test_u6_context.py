"""Small independent GFF3 fixtures exercise annotation evidence boundaries."""

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "u6_context.py"
SPEC = importlib.util.spec_from_file_location("u6_context", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)
GffContext = MODULE.GffContext


def row(kind, start, end, attrs, strand="+", seqid="chr1", source="fixture"):
    return f"{seqid}\t{source}\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n"


class GffContextTests(unittest.TestCase):
    def context(self, *rows, seq_lengths=None):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        path = Path(directory.name) / "test.gff3"
        path.write_text("##gff-version 3\n" + "".join(rows), encoding="utf-8")
        return GffContext(path, seq_lengths=seq_lengths)

    def test_multilevel_and_multiparent_graph_preserves_both_genes(self):
        context = self.context(
            row("gene", 1, 800, "ID=g1;gene_biotype=protein_coding"),
            row("gene", 5, 900, "ID=g2;gene_biotype=protein_coding"),
            row("primary_transcript", 10, 800, "ID=pre;Parent=g1"),
            row("mRNA", 10, 800, "ID=t1;Parent=pre"),
            row("mRNA", 10, 900, "ID=t2;Parent=g2"),
            row("exon", 100, 200, "ID=shared;Parent=t1,t2"),
        )
        audit = context.audit("chr1", 150, 160, "+")
        exon = next(x for x in audit["overlaps"] if x["feature_id"] == "shared")
        self.assertEqual(exon["gene_ids"], ["g1", "g2"])
        self.assertEqual(exon["parent_ids"], ["t1", "t2"])
        self.assertEqual(exon["biotype"], "protein_coding")
        self.assertEqual(exon["overlap_bp"], 11)
        self.assertEqual(exon["region_fraction"], 1.0)

    def test_self_u6_and_children_excluded_but_other_ncrna_retained(self):
        context = self.context(
            row("gene", 100, 201, "ID=gu;gene_biotype=snRNA"),
            row("snRNA", 100, 201, "ID=ru;Parent=gu;product=U6%20spliceosomal%20RNA"),
            row("exon", 100, 201, "ID=eu;Parent=ru"),
            row("snRNA", 190, 280, "ID=other;product=U5%20RNA"),
        )
        supported = context.supporting_u6_ids("chr1", 100, 201, "+")
        self.assertEqual(supported, {"gu", "ru"})
        audit = context.audit("chr1", 100, 201, "+", supported)
        self.assertEqual(audit["excluded_self_feature_count"], 3)
        self.assertEqual([x["feature_id"] for x in audit["overlaps"]], ["other"])
        self.assertEqual(audit["counts_by_type"]["exon"], 0)

    def test_containing_gene_is_not_automatically_self(self):
        context = self.context(
            row("gene", 1, 1000, "ID=broad;Name=U6;gene_biotype=snRNA"),
            row("snRNA", 100, 201, "ID=ru;Parent=broad;product=U6%20RNA"),
            row("exon", 100, 201, "ID=eu;Parent=ru"),
        )
        supported = context.supporting_u6_ids("chr1", 100, 201, "+")
        self.assertEqual(supported, {"ru"})
        audit = context.audit("chr1", 100, 201, "+", supported)
        self.assertEqual([x["feature_id"] for x in audit["overlaps"]], ["broad"])

    def test_gene_only_u6_label_and_u6atac_do_not_supply_self_ids(self):
        context = self.context(
            row("gene", 1, 1000, "ID=broad;Name=U6"),
            row("snRNA", 100, 201, "ID=atac;Name=U6atac"),
        )
        self.assertEqual(context.supporting_u6_ids("chr1", 100, 201, "+"), set())

    def test_same_coordinate_gene_with_foreign_sibling_not_excluded(self):
        context = self.context(
            row("gene", 100, 201, "ID=g;gene_biotype=snRNA"),
            row("snRNA", 100, 201, "ID=ru;Parent=g;Name=U6"),
            row("snRNA", 100, 201, "ID=r5;Parent=g;Name=U5"),
        )
        supported = context.supporting_u6_ids("chr1", 100, 201, "+")
        self.assertEqual(supported, {"ru"})
        audit = context.audit("chr1", 100, 201, "+", supported)
        self.assertEqual({x["feature_id"] for x in audit["overlaps"]}, {"g", "r5"})

    def test_whole_chromosome_region_is_not_a_foreign_feature(self):
        context = self.context(
            row("region", 1, 1000000, "ID=chr1;Is_circular=false"),
            row("gene", 100, 200, "ID=g"),
        )
        audit = context.audit("chr1", 500, 600, "+")
        self.assertEqual(audit["overlaps"], [])
        self.assertEqual(audit["coordinate_container_records_ignored"], 1)
        self.assertEqual(context._starts["chr1"], [100])

    def test_self_feature_shared_with_foreign_parent_is_retained(self):
        context = self.context(
            row("gene", 100, 201, "ID=gu;gene_biotype=snRNA"),
            row("gene", 100, 500, "ID=foreign;gene_biotype=protein_coding"),
            row("snRNA", 100, 201, "ID=ru;Parent=gu;Name=U6"),
            row("mRNA", 100, 500, "ID=tm;Parent=foreign"),
            row("exon", 100, 201, "ID=shared;Parent=ru,tm"),
        )
        audit = context.audit("chr1", 100, 201, "+", {"gu"})
        shared = next(x for x in audit["overlaps"] if x["feature_id"] == "shared")
        self.assertTrue(shared["shared_with_self"])
        self.assertEqual(shared["gene_ids"], ["foreign", "gu"])

    def test_reverse_cds_and_intron_only_gene_span_are_distinct(self):
        context = self.context(
            row("gene", 1, 1000, "ID=g;gene_biotype=protein_coding", "-"),
            row("mRNA", 1, 1000, "ID=t;Parent=g", "-"),
            row("exon", 100, 200, "ID=e;Parent=t", "-"),
            row("CDS", 120, 190, "ID=c;Parent=t", "-"),
        )
        intron = context.audit("chr1", 400, 500, "+")
        self.assertEqual(intron["counts_by_type"]["gene_span"], 1)
        self.assertEqual(intron["counts_by_type"]["CDS"], 0)
        self.assertEqual(intron["counts_by_type"]["exon"], 0)
        cds = context.audit("chr1", 150, 160, "+")
        self.assertEqual(cds["counts_by_type_and_strand"]["CDS"]["opposite"], 1)

    def test_lncrna_utr_repeat_and_pseudogene_independent_categories(self):
        context = self.context(
            row("gene", 100, 200, "ID=lg;gene_biotype=lncRNA"),
            row("lnc_RNA", 100, 200, "ID=lt;Parent=lg"),
            row("exon", 100, 200, "ID=le;Parent=lt"),
            row("five_prime_UTR", 110, 120, "ID=u5", "-"),
            row("three_prime_UTR", 130, 140, "ID=u3"),
            row("repeat_region", 145, 170, "ID=rep", ".", source="RepeatMasker"),
            row("pseudogene", 160, 300, "ID=pg;gene_biotype=processed_pseudogene"),
        )
        audit = context.audit("chr1", 100, 200, "+")
        self.assertEqual(audit["counts_by_type"]["ncRNA"], 3)
        self.assertEqual(audit["counts_by_type"]["five_prime_UTR"], 1)
        self.assertEqual(audit["counts_by_type"]["three_prime_UTR"], 1)
        self.assertEqual(audit["counts_by_type"]["repeat"], 1)
        self.assertEqual(audit["counts_by_type"]["pseudogene"], 1)
        repeat = next(x for x in audit["overlaps"] if x["feature_id"] == "rep")
        self.assertEqual(repeat["source"], "RepeatMasker")
        self.assertEqual(repeat["strand_relation"], "unknown")
        self.assertIsNone(audit["counts_by_type"]["CDS"])

    def test_missing_gff_and_missing_seqid_are_na(self):
        missing = GffContext(None).audit("chr1", 1, 100, "+")
        self.assertEqual(missing["reason"], "no_gff")
        self.assertTrue(all(v is None for v in missing["counts_by_type"].values()))
        context = self.context(row("gene", 1, 100, "ID=g"))
        mismatch = context.audit("Chr1", 1, 100, "+")
        self.assertEqual(mismatch["status"], "context_not_assessed")
        self.assertEqual(mismatch["reason"], "seqid_mismatch")
        self.assertTrue(all(v == "NA" for v in mismatch["annotation_type_availability"].values()))

    def test_missing_type_on_sequence_is_na_not_false_zero(self):
        context = self.context(
            row("gene", 1, 100, "ID=g1"),
            row("CDS", 1, 100, "ID=c2", seqid="chr2"),
        )
        audit = context.audit("chr1", 101, 200, "+")
        self.assertEqual(audit["status"], "assessed")
        self.assertEqual(audit["counts_by_type"]["gene_span"], 0)
        self.assertIsNone(audit["counts_by_type"]["CDS"])
        self.assertEqual(audit["annotation_type_availability"]["CDS"], "NA")

    def test_aliases_decoding_and_interval_boundaries(self):
        context = self.context(
            row("gene", 100, 201, "ID=g%2C1;Alias=old%2Cname"),
            row("snRNA", 100, 201, "ID=r;Parent=g%2C1;Name=U6"),
            row("exon", 201, 220, "ID=e"),
        )
        audit = context.audit("chr1", 201, 201, "+", {"old,name"})
        self.assertEqual(audit["resolved_self_ids"], ["g,1"])
        self.assertEqual(len(audit["overlaps"]), 1)
        self.assertEqual(audit["overlaps"][0]["overlap_bp"], 1)
        self.assertEqual(context.audit("chr1", 221, 222, "+")["overlaps"], [])

    def test_fasta_bounds_and_orphan_parent_are_not_silently_clean(self):
        context = self.context(row("exon", 10, 20, "ID=e;Parent=absent"), seq_lengths={"chr1": 100})
        hit = context.audit("chr1", 10, 20, "+")["overlaps"][0]
        self.assertEqual(hit["missing_parent_ids"], ["absent"])
        self.assertIsNone(hit["biotype"])
        self.assertEqual(context.audit("chr1", 90, 101, "+")["reason"], "invalid_region")
        invalid = self.context(row("gene", 1, 200, "ID=g"), seq_lengths={"chr1": 100})
        self.assertEqual(invalid.audit("chr1", 1, 20, "+")["reason"], "gff_coordinates_out_of_bounds")


if __name__ == "__main__":
    unittest.main()

