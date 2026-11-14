"""GRCh37 coordinates and BIC legacy names (bin/brca12lib/legacy.py).

    python3 -m unittest discover -s tests
"""
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))

from brca12lib import legacy  # noqa: E402


class Grch37Test(unittest.TestCase):
    def test_brca1_founder_position(self):
        # BRCA1 c.68_69del: GRCh38 chr17:43124027 -> GRCh37 chr17:41276044 (VCF position)
        self.assertEqual(legacy.grch37("chr17", 43124027), 41276044)

    def test_brca2_locus_bounds(self):
        self.assertEqual(legacy.grch37("chr13", 32315077), 32889214)
        self.assertEqual(legacy.grch37("chr13", 32400268), 32974405)

    def test_outside_loci_is_not_converted(self):
        self.assertIsNone(legacy.grch37("chr17", 1000))
        self.assertIsNone(legacy.grch37("chr2", 43124027))


class BicTest(unittest.TestCase):
    def test_founder_variants(self):
        self.assertIn("185delAG", legacy.bic_names("BRCA1", "ENST00000357654.9:c.68_69del"))
        self.assertIn("5382insC", legacy.bic_names("BRCA1", "c.5266dup"))
        self.assertEqual(legacy.bic_names("BRCA2", "c.5946del"), ["6174delT"])

    def test_hgvs_position_comes_first(self):
        self.assertEqual(legacy.bic_names("BRCA1", "c.68_69del"), ["187delAG", "185delAG"])

    def test_substitutions(self):
        self.assertEqual(legacy.bic_name("BRCA1", "c.181T>G"), "300T>G")
        self.assertEqual(legacy.bic_name("BRCA2", "c.-26G>A"), "203G>A")     # no c.0 in 5'UTR

    def test_intron_names_follow_legacy_exons(self):
        self.assertEqual(legacy.bic_name("BRCA1", "c.-19-10T>C"), "IVS1-10T>C")
        self.assertEqual(legacy.bic_name("BRCA1", "c.-20+11C>T"), "IVS1+11C>T")
        self.assertEqual(legacy.bic_name("BRCA1", "c.135-3641C>T"), "IVS3-3641C>T")   # no legacy exon 4

    def test_long_deletion_uses_length(self):
        self.assertEqual(legacy.bic_name("BRCA1", "c.949_953del"), "1068del5")

    def test_unsupported_notation(self):
        self.assertEqual(legacy.bic_names("BRCA1", "c.100_102delinsT"), [])
        self.assertEqual(legacy.bic_names("TP53", "c.1A>G"), [])


if __name__ == "__main__":
    unittest.main()
