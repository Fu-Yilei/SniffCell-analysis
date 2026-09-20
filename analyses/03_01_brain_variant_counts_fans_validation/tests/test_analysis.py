"""Methodology regression tests. Fixture records are synthetic, not study calls."""
import csv
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS.parents[1] / "03_02_brain_tr_motif_composition" / "scripts"))
from rebuild_candidate_counts import category, number, shared_match, write_table
from summarize_fans_validation import summarize
from derive_tr_motif_composition import sequence_class, derive
from measure_fans_bam_evidence import indels, tr_values, sv_values


class MethodTests(unittest.TestCase):
    def test_sample_prefixed_group_and_ambiguous_group(self):
        self.assertEqual(category("pbmc", "donor.T-cell_NK-cell_B-cell"), "Lymphocyte group")
        with self.assertRaises(ValueError):
            category("brain", "Neuron+Oligodendrocyte")

    def test_nan_does_not_pass_threshold(self):
        for value in ["nan", "inf", ""]:
            with self.assertRaises(ValueError):
                number(value)

    def test_shared_sv_matching_boundaries(self):
        row = dict(id="candidate", sv_chr="test_contig", sv_pos="2000", sv_type="INS", sv_len="100")
        self.assertTrue(shared_match(row, [("other", "test_contig", 3000, "INS", 50)]))
        self.assertFalse(shared_match(row, [("other", "test_contig", 3001, "INS", 50)]))
        self.assertFalse(shared_match(row, [("other", "test_contig", 2000, "DEL", 100)]))
        self.assertFalse(shared_match(row, [("other", "test_contig", 2000, "INS", 49)]))

    def test_unique_locus_union_and_unscored(self):
        rows = [dict(variant_type="TR", tissue="brain", donor=d, chrom="test_contig",
                     start=s, end=s+10, observed_in_paired_FANS=status)
                for d, s, status in [("d1", 10, "Yes"), ("d2", 10, "No"), ("d1", 30, "unknown")]]
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            write_table(path / "input.tsv", list(rows[0]), rows)
            summarize(path / "input.tsv", path / "out")
            with (path / "out/fans_validation_denominators.tsv").open() as handle:
                result = list(csv.DictReader(handle, delimiter="\t"))
            record = next(r for r in result if r["tissue"] == "brain" and
                          r["variant_class"] == "TR" and r["counting_unit"] == "unique_loci")
            self.assertEqual((record["observed"], record["eligible"], record["unscored"]), ("1", "1", "1"))
            self.assertEqual(record["loci_with_mixed_yes_no"], "1")
            write_table(path / "input.tsv", list(rows[0]), rows + rows[:1])
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                summarize(path / "input.tsv", path / "out")

    def test_motif_thresholds_and_invalid_sequence(self):
        self.assertEqual(sequence_class("AAAAC"), "AT-rich")
        self.assertEqual(sequence_class("CCCCA"), "GC-rich")
        self.assertEqual(sequence_class("ATGC"), "Mixed")
        for motif in ["", "AN"]:
            with self.assertRaises(ValueError):
                sequence_class(motif)

    def test_motif_counts_donors_not_unique_loci(self):
        rows = [dict(variant_type="TR", tissue="brain", donor=d, chrom="test_contig",
                     start=10, end=20, motif="ATAT") for d in ["d1", "d2"]]
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            write_table(path / "input.tsv", list(rows[0]), rows)
            derive(path / "input.tsv", path / "out")
            with (path / "out/tr_motif_composition.tsv").open() as handle:
                output = list(csv.DictReader(handle, delimiter="\t"))
            row = next(r for r in output if r["tissue"] == "Brain" and r["category"] == "4")
            self.assertEqual(row["count"], "2")  # string length, not minimal period 2

    def test_cigar_reference_skip_does_not_shift_indel(self):
        class Read:
            reference_start = 0
            reference_end = 1000
            is_unmapped = is_secondary = is_supplementary = False
            cigartuples = [(0, 100), (3, 100), (1, 30), (0, 100), (2, 10)]
        self.assertEqual(list(indels(Read())), [(200, 30), (300, -10)])
        class Bam:
            references = ["test_contig"]
            def fetch(self, *args):
                return iter([Read()])
        self.assertEqual(tr_values(Bam(), "test_contig", 200, 250, 100), [20])
        self.assertEqual(sv_values(Bam(), "test_contig", 200, 60, 150), (1, 1))


if __name__ == "__main__":
    unittest.main()
