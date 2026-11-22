"""Known weak spot of the method (BRCA1 3' UTR poly-T stretch): indels there carry a quality flag."""
import importlib.util
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
spec = importlib.util.spec_from_file_location("brca12_report", os.path.join(ROOT, "bin", "brca12_report.py"))
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)

from brca12lib import i18n  # noqa: E402


class WeakSpotTest(unittest.TestCase):
    def test_indels_in_the_poly_t_stretch_are_flagged(self):
        self.assertTrue(report.in_weak_spot({"chrom": "chr17", "pos": 43044803, "ref": "AT", "alt": "A"}))
        self.assertTrue(report.in_weak_spot({"chrom": "17", "pos": 43044810, "ref": "T", "alt": "TT"}))
        self.assertTrue(report.in_weak_spot({"chrom": "chr17", "pos": 43044817, "ref": "TTT", "alt": "T"}))

    def test_snvs_and_indels_elsewhere_are_not(self):
        self.assertFalse(report.in_weak_spot({"chrom": "chr17", "pos": 43044810, "ref": "T", "alt": "C"}))
        self.assertFalse(report.in_weak_spot({"chrom": "chr17", "pos": 43044900, "ref": "AT", "alt": "A"}))
        self.assertFalse(report.in_weak_spot({"chrom": "chr13", "pos": 43044810, "ref": "AT", "alt": "A"}))

    def test_flag_is_translated(self):
        self.assertTrue(i18n.tr(report.WEAK_SPOT_FLAG, "en").startswith("known weak spot: poly-T stretch"))


if __name__ == "__main__":
    unittest.main()
