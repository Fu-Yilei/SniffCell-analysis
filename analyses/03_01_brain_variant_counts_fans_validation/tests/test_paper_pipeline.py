"""Synthetic-record tests; no real candidate calls are shipped as fixtures."""
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from paper_source import read_candidates, summarize, depth_counts, records_sha256
from derive_brain_tr_candidates import mechanistic_filter, reproducibility_class

gaa_path = SCRIPTS.parents[1] / "03_02_brain_tr_motif_composition/scripts/analyze_gaa_enrichment.py"
spec = importlib.util.spec_from_file_location("gaa_analysis", gaa_path)
gaa = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gaa)


def record(i=1, **updates):
    row = dict(variant_id=f"synthetic{i}", variant_type="TR", tissue="brain", donor="bcontrol1",
               chrom="synthetic_contig", start=str(100 * i), end=str(100 * i + 50),
               cell_type="Neuron", observed_in_paired_FANS="Yes", coverages_detected="full;60x")
    row.update(updates)
    return row


class PaperTests(unittest.TestCase):
    def parse(self, rows):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.json"
            h = list(rows[0])
            path.write_text(json.dumps(dict(values=[["Table title"], h, *[[r[k] for k in h] for r in rows], []])))
            return read_candidates(path)

    def test_header_blank_rows_whitespace_and_checksum(self):
        rows = self.parse([record(variant_id=" synthetic1 ")])
        self.assertEqual(rows[0]["variant_id"], "synthetic1")
        self.assertEqual(rows[0]["source_row"], 3)
        self.assertEqual(rows[0]["coverages_detected"], "60x;full")
        self.assertEqual(records_sha256(rows), records_sha256([dict(rows[0], source_row=900)]))
        self.assertNotEqual(records_sha256(rows), records_sha256([dict(rows[0], motif="GAA")]))

    def test_reject_duplicate_ids_and_unknown_validation(self):
        for rows in ([record(), record()], [record(observed_in_paired_FANS="Unscored")],
                     [record(coverages_detected="full;full")], [record(start="nan")]):
            with self.assertRaises(ValueError):
                self.parse(rows)

    def test_unique_loci_any_donor_not_sum_of_celltypes(self):
        rows = [record(1), record(2, donor="bcontrol2", start="100", end="150", cell_type="Oligodendrocyte"),
                record(3, variant_type="SV", observed_in_paired_FANS="No"),
                record(4, tissue="pbmc", donor="pbmc1", cell_type="Monocyte")]
        validation, composition, unique = summarize(self.parse(rows))
        self.assertEqual((unique["observed"], unique["total"]), (1, 1))
        self.assertEqual(validation[0]["total"], 3)
        self.assertEqual(sum(r["count"] for r in composition), 4)

    def test_depth_membership_and_zeros(self):
        rows = [record(i, donor=d, tissue=t, cell_type="Monocyte" if t == "pbmc" else "Neuron")
                for i, (d, t) in enumerate([( "pbmc1", "pbmc"), ("pbmc2", "pbmc"),
                    ("bcontrol1", "brain"), ("bcontrol2", "brain"), ("bmsa", "brain"), ("bparkinsons", "brain")], 1)]
        counts = depth_counts(self.parse(rows))
        self.assertEqual(len(counts), 18)
        self.assertEqual(sum(r["total"] for r in counts if r["coverage"] == "30x"), 0)
        self.assertEqual(sum(r["TR"] for r in counts if r["coverage"] == "full"), 6)

    def test_gaa_strand_rotation_long_motif_and_negative(self):
        for motif in ("GAA", "AGA", "AAG", "TTC", "CTT", "TCT", "AAAAGAAAAAAAAG"):
            self.assertTrue(gaa.is_gaa_element(motif), motif)
        for motif in ("A", "AG", "CAG", "CCG"):
            self.assertFalse(gaa.is_gaa_element(motif), motif)
        with self.assertRaises(ValueError):
            gaa.is_gaa_element("GAA,TTC")

    def test_fisher_known_results(self):
        self.assertAlmostEqual(gaa.fisher(18, 23, 13, 30)[2], .2588006906267615)
        self.assertAlmostEqual(gaa.fisher(18, 23, 7, 36)[2], .008267487570047)  # prose contingency, not inputs

    def test_catalog_tie_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "catalog.bed"
            path.write_text("synthetic_contig\t100\t150\tMOTIFS=GAA\nsynthetic_contig\t100\t150\tMOTIFS=CAG\n")
            with self.assertRaises(ValueError):
                gaa.match_catalog(self.parse([record()]), path)

    def test_historical_filter_boundaries(self):
        row = dict(support="3", delta="2000", reproducible="False", support_gc_median=".4",
                   support_entropy1_median=".95", support_top4mer_frac_median=".05", support_dinuc_purity_median=".1")
        self.assertTrue(mechanistic_filter(row)[0])
        self.assertFalse(mechanistic_filter(dict(row, delta="2001"))[0])
        self.assertTrue(mechanistic_filter(dict(row, delta="2001", reproducible="True"))[0])


if __name__ == "__main__":
    unittest.main()
