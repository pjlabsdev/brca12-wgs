"""Tests for the bundled ENIGMA v1.2 tables: PVS1/PM5_PTC (Table 4), PS3/BS3 (Table 9), PP4/BP5.

    python3 -m unittest discover -s tests
"""
import json
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))

from brca12lib import enigma, enigma_tables  # noqa: E402


def codes(gene, consequence, hgvsc, hgvsp="", exon=""):
    csq = {"Consequence": consequence, "HGVSc": f"ENST:{hgvsc}", "HGVSp": hgvsp, "EXON": exon}
    return [h["code"] for h in enigma.acmg_hints(gene, csq, None)]


class AssetsTest(unittest.TestCase):
    def test_source_and_rows(self):
        with open(os.path.join(ROOT, "assets", "enigma", "SOURCE.json")) as handle:
            meta = json.load(handle)
        self.assertEqual(meta["version"], "1.2")
        rows = {name: info["rows"] for name, info in meta["tracks"].items()}
        self.assertEqual(rows["BRCAsplicing"], 519)
        self.assertEqual(rows["BRCAfunctionalAssays"], 4729)
        self.assertEqual(rows["BRCAmfa"], 13182)
        self.assertTrue(enigma_tables.available())


class TableFourTest(unittest.TestCase):
    def test_ptc_in_brca1_exon_10_gets_pvs1_and_pm5_ptc(self):
        found = codes("BRCA1", "frameshift_variant", "c.3000del")
        self.assertIn("PVS1", found)
        self.assertIn("PM5_PTC_Strong", found)
        self.assertNotIn("PVS1 (Kandidat)", found)

    def test_ptc_after_nmd_boundary_is_not_pvs1(self):
        found = codes("BRCA1", "stop_gained", "c.5565G>A")       # c.5563-c.5571: PVS1_N/A, PM5_N/A
        self.assertIn("PVS1 nicht anwendbar", found)
        self.assertNotIn("PVS1", found)

    def test_brca2_k3326ter_keeps_known_exception(self):
        found = codes("BRCA2", "stop_gained", "c.9976A>T", "ENSP:p.Lys3326Ter", "27/27")
        self.assertIn("Bekannte Ausnahme", found)
        self.assertFalse(any(c == "PVS1" for c in found))

    def test_brca2_c1909_plus1_is_pvs1_in_v1_2(self):
        self.assertIn("PVS1", codes("BRCA2", "splice_donor_variant", "c.1909+1G>T"))

    def test_brca1_c4096_plus1_is_moderate(self):
        self.assertIn("PVS1_Moderate", codes("BRCA1", "splice_donor_variant", "c.4096+1G>T"))

    def test_splice_indel_falls_back_to_candidate(self):
        self.assertIn("PVS1 (Kandidat)", codes("BRCA1", "splice_donor_variant", "c.4096+1del"))

    def test_start_codon(self):
        found = codes("BRCA2", "start_lost", "c.1A>G")
        self.assertTrue(any(c.startswith("PVS1") for c in found))

    def test_exon_deletion_row(self):
        row = enigma_tables.exon_event_row("BRCA1", 10, "DEL")
        self.assertIsNotNone(row)
        self.assertTrue(row["code"].startswith("PVS1"))
        self.assertIsNotNone(enigma_tables.exon_event_row("BRCA2", 11, "DUP"))


class TableNineAndLikelihoodTest(unittest.TestCase):
    def test_functional_bs3(self):
        self.assertIn("BS3", codes("BRCA2", "missense_variant", "c.47A>G"))

    def test_functional_ps3(self):
        self.assertIn("PS3", codes("BRCA2", "start_lost", "c.3G>A"))

    def test_multifactorial_bp5(self):
        self.assertIn("BP5_VeryStrong", codes("BRCA1", "3_prime_UTR_variant", "c.*1382dup"))

    def test_unknown_variant_has_no_table_codes(self):
        found = codes("BRCA1", "missense_variant", "c.9999999A>G")
        self.assertFalse(any(c.startswith(("PS3", "BS3", "PP4", "BP5")) for c in found))

    def test_parse_code(self):
        self.assertEqual(enigma_tables.parse_code("PVS1_Strong (RNA)"), ("PVS1", "Strong", True))
        self.assertEqual(enigma_tables.parse_code("PM5_N/A"), ("PM5_PTC", "N/A", False))


if __name__ == "__main__":
    unittest.main()
