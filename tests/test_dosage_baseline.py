"""Plausibility of the chromosome baseline for the per-exon dosage hints (structural.dosage_baseline)."""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin"))
from brca12lib import structural  # noqa: E402


def exons(depth, chrom="chr17", n=10):
    return [{"gene": "BRCA1", "exon": i, "label": f"Exon {i}", "chrom": chrom, "meanDepth": depth,
             "dosageRatio": round(depth / 1.0, 2), "status": "ok"} for i in range(1, n + 1)]


class DosageBaseline(unittest.TestCase):
    def test_whole_genome_baseline_is_plausible(self):
        ok, info = structural.dosage_baseline(exons(30.0), {"chr17": 31.0, "chr13": 30.0})
        self.assertTrue(ok)

    def test_partial_input_baseline_is_rejected(self):
        # only the gene regions were read: the chromosome median is near zero
        ok, info = structural.dosage_baseline(exons(30.0), {"chr17": 1.0, "chr13": 0.9})
        self.assertFalse(ok)
        self.assertEqual(info["chromosome"], "chr17")
        self.assertAlmostEqual(info["medianExonDepth"], 30.0)

    def test_no_baseline(self):
        ok, _ = structural.dosage_baseline(exons(30.0), {})
        self.assertTrue(ok)


if __name__ == "__main__":
    unittest.main()
