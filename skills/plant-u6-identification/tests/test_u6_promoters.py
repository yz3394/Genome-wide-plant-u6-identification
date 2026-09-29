"""Pure-function edge tests; private project regression data are not bundled."""

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from u6_promoters import extract_upstream, reverse_complement, scan_motifs, scan_termination


def promoter_with(*elements):
    sequence = list("C" * 150)
    for start, bases in elements:
        sequence[start + 150:start + 150 + len(bases)] = bases
    return "".join(sequence)


class ExtractionTests(unittest.TestCase):
    def setUp(self):
        self.sequence = "ACGTRYNACCTTAGGC"
        self.fetch = lambda seqid, start, end: self.sequence[start - 1:end]

    def test_plus_exact(self):
        result = extract_upstream(self.fetch, "chr", "+", 10, 5, 16)
        self.assertEqual((result["start1"], result["end1"]), (5, 9))
        self.assertEqual(result["sequence"], self.sequence[4:9])
        self.assertEqual(result["N_count"], 1)
        self.assertEqual(result["ambiguous_base_count"], 3)
        self.assertFalse(result["truncated"])

    def test_minus_is_oriented_and_truncated(self):
        result = extract_upstream(self.fetch, "chr", "-", 13, 8, 16)
        self.assertEqual(result["sequence"], reverse_complement(self.sequence[13:16]))
        self.assertEqual((result["start1"], result["end1"]), (14, 16))
        self.assertEqual(result["first_relative_position"], -3)
        self.assertTrue(result["truncated"])

    def test_empty_edges_and_invalid_boundary_do_not_fetch(self):
        def no_fetch(*args):
            self.fail("Empty intervals must not call FASTA fetch")
        for strand, tss in [("+", 1), ("-", 16)]:
            result = extract_upstream(no_fetch, "chr", strand, tss, 500, 16)
            self.assertEqual(result["actual_length"], 0)
            self.assertEqual(result["extraction_status"], "no_upstream_sequence_at_reference_edge")
        result = extract_upstream(no_fetch, "chr", "+", -8, 500, 16)
        self.assertEqual(result["extraction_status"], "boundary_outside_reference")

    def test_plus_edge_no_padding(self):
        result = extract_upstream(self.fetch, "chr", "+", 4, 500, 16)
        self.assertEqual(result["sequence"], "ACG")
        self.assertEqual(result["actual_length"], 3)
        self.assertEqual(result["first_relative_position"], -3)


class MotifTests(unittest.TestCase):
    def test_exact_pair_four_spacing_definitions(self):
        result = scan_motifs(promoter_with((-66, "GTCCCACATCG"), (-31, "TTTATATA")), -150)
        best = result["best"]
        self.assertEqual(best["motif_pair_class"], "P1_canonical_geometry")
        self.assertEqual(best["use_tata_start_separation"], 35)
        self.assertEqual(best["use_tata_edge_distance"], 25)
        self.assertEqual(best["use_tata_intervening_bases"], 24)
        self.assertEqual(best["use_tata_center_distance"], 33.5)
        self.assertGreater(len(result["pairs"]), 1)

    def test_noncanonical_pair_retained_without_family_rejection(self):
        result = scan_motifs(promoter_with((-103, "GTCCCACATCG"), (-49, "TTTATATA")), -150, "generic_provisional")
        self.assertTrue(result["pairs"])
        self.assertEqual(result["best"]["motif_pair_class"], "provisional_uncertain_geometry")
        self.assertEqual(result["calibration_status"], "provisional_not_cross_plant_calibrated")
        self.assertIn("not_evidence", result["absence_interpretation"])

    def test_at_rich_is_not_positive_without_use(self):
        result = scan_motifs("AT" * 75, -150, "generic_provisional")
        self.assertEqual(result["best"], {})
        self.assertEqual(result["motif_evidence_status"], "uncertain_or_no_pair_not_U6_exclusion")

    def test_incomplete_and_empty_windows(self):
        for sequence in ["", "AACG", "N" * 20]:
            result = scan_motifs(sequence, -len(sequence))
            self.assertEqual(result["scan_status"], "insufficient_window")
            self.assertEqual(result["best"], {})
        partial = scan_motifs("A" * 50, -50)
        self.assertEqual(partial["scan_status"], "partial_window")
        self.assertEqual(partial["best"]["motif_pair_class"], "P3_uncertain")

    def test_negative_strand_extraction_preserves_motif_orientation(self):
        promoter = promoter_with((-66, "GTCCCACATCG"), (-31, "TTTATATA"))
        reference = "G" + reverse_complement(promoter)
        fetch = lambda seqid, start, end: reference[start - 1:end]
        extracted = extract_upstream(fetch, "negative", "-", 1, 150, len(reference))
        self.assertEqual(extracted["sequence"], promoter)
        result = scan_motifs(extracted["sequence"], extracted["first_relative_position"])
        self.assertEqual(result["best"]["use_start_rel_tss"], -66)
        self.assertEqual(result["best"]["motif_pair_class"], "P1_canonical_geometry")

    def test_n_within_motif_is_not_supported_evidence(self):
        result = scan_motifs(promoter_with((-66, "GTCCCACNTCG"), (-31, "TTTATATA")), -150, "generic_provisional")
        self.assertNotEqual(result["motif_evidence_status"], "paired_heuristic_evidence")

    def test_optional_msp_has_no_effect_on_best_or_ranking(self):
        sequence = promoter_with((-130, "AGCCCG"), (-66, "GTCCCACATCG"), (-31, "TTTATATA"), (-10, "AGCCCG"))
        off, on = scan_motifs(sequence, -150), scan_motifs(sequence, -150, msp=True)
        self.assertEqual(off["best"], on["best"])
        self.assertEqual(off["pairs"], on["pairs"])
        self.assertEqual(len(on["msp_hits"]), 1)
        self.assertEqual(on["msp_hits"][0]["start_rel_tss"], -130)
        unknown = scan_motifs(promoter_with((-130, "AGCCCG")), -150, msp=True)
        self.assertEqual(unknown["msp_hits"][0]["relation_to_USE"], "USE_relation_unresolved")


class TerminationTests(unittest.TestCase):
    def test_internal_tract_not_substituted_for_expected_terminator(self):
        result = scan_termination("A" * 19 + "TTTT" + "A" * 117, 102)
        self.assertEqual(result["termination_t4plus_expected_window"], "no")
        self.assertEqual(result["termination_run_sequence"], "")
        self.assertEqual(len(result["internal_t4plus_hits"]), 1)

    def test_expected_and_internal_are_separate(self):
        result = scan_termination("A" * 19 + "TTTT" + "A" * 75 + "TTTTTT" + "A" * 36, 102)
        self.assertEqual(result["termination_start_rel_tss"], 99)
        self.assertEqual(len(result["internal_t4plus_hits"]), 1)
        self.assertEqual(len(result["body_overlapping_t4plus_hits"]), 2)

    def test_missing_window_is_unassessable_not_negative(self):
        result = scan_termination("A" * 40, 102)
        self.assertEqual(result["termination_t4plus_expected_window"], "not_fully_assessable")
        self.assertEqual(result["termination_scan_status"], "unavailable_window")
        self.assertEqual(result["termination_run_sequence"], "")

    def test_window_boundary_crossing_is_reported_but_not_relabelled(self):
        result = scan_termination("A" * 83 + "TTTTTT" + "A" * 60, 102)
        self.assertEqual(result["termination_t4plus_expected_window"], "no")
        self.assertEqual(len(result["terminal_window_overlapping_t4plus_hits"]), 1)
        self.assertEqual(result["termination_expected_hits"], [])

    def test_ambiguous_terminal_window_is_not_a_supported_absence(self):
        result = scan_termination("A" * 99 + "NNNN" + "A" * 40, 102)
        self.assertEqual(result["termination_scan_status"], "complete_window")
        self.assertEqual(result["termination_window_ambiguous_base_count"], 4)
        self.assertEqual(result["termination_t4plus_expected_window"], "not_fully_assessable")




if __name__ == "__main__":
    unittest.main()

