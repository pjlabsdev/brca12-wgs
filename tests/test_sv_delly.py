"""Structural variants from Delly (default SV caller) in the report parser (bin/brca12lib/structural.py)."""
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
from brca12lib import structural  # noqa: E402

REGIONS = [("chr17", 43044294, 43170327, "BRCA1"), ("chr13", 32315076, 32400268, "BRCA2")]
EXONS = [{"gene": "BRCA1", "label": "Exon 14 (Legacy 15)", "chrom": "chr17", "start": 43074331, "end": 43074521},
         {"gene": "BRCA1", "label": "Exon 15 (Legacy 16)", "chrom": "chr17", "start": 43070928, "end": 43071238}]
DELLY = (
    "##fileformat=VCFv4.2\n##source=DELLY\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA18949\n"
    "chr17\t43069972\tDEL00000001\tN\t<DEL>\t960\tPASS\tPRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv1.3.3;END=43075601;PE=13;SR=7;CT=3to5\tGT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV\t0/1:-50,0,-90:10000:PASS:900:820:880:1:30:13:25:7\n"
    "chr17\t43100000\tINV00000002\tN\t<INV>\t120\tPASS\tIMPRECISE;SVTYPE=INV;END=43105000;PE=5;CT=3to3\tGT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV\t0/1:-10,0,-20:60:PASS:0:0:0:2:20:5:0:0\n"
    "chr17\t43110000\tBND00000003\tN\tN]chr5:1000]\t300\tPASS\tPRECISE;SVTYPE=BND;CHR2=chr5;POS2=1000;PE=8;SR=4;CT=3to3\tGT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV\t0/1:-20,0,-40:90:PASS:0:0:0:2:25:8:20:4\n"
    "chr2\t5000\tDEL00000004\tN\t<DEL>\t500\tPASS\tPRECISE;SVTYPE=DEL;END=9000;PE=10\tGT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV\t0/1:-20,0,-40:90:PASS:0:0:0:1:25:10:0:0\n"
)


class DellyParser(unittest.TestCase):
    def parse(self, text):
        with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as handle:
            handle.write(text)
        try:
            return structural.structural_variants(handle.name, REGIONS, EXONS)
        finally:
            os.unlink(handle.name)

    def test_deletion_fields(self):
        svs = self.parse(DELLY)
        dele = next(s for s in svs if s["type"] == "DEL")
        self.assertEqual((dele["chrom"], dele["pos"], dele["end"], dele["gene"]), ("chr17", 43069972, 43075601, "BRCA1"))
        self.assertEqual(dele["caller"], "Delly")
        self.assertEqual(dele["pairedReads"], "30,13")
        self.assertEqual(dele["splitReads"], "25,7")
        self.assertEqual(dele["length"], 43075601 - 43069972)
        self.assertEqual(dele["genotype"], "0/1")
        self.assertEqual(set(dele["affectedExons"]), {"Exon 14 (Legacy 15)", "Exon 15 (Legacy 16)"})
        self.assertNotIn("unpräzise Bruchpunkte", dele["flags"])

    def test_inversion_imprecise_and_bnd(self):
        svs = self.parse(DELLY)
        inv = next(s for s in svs if s["type"] == "INV")
        self.assertIn("unpräzise Bruchpunkte", inv["flags"])
        bnd = next(s for s in svs if s["type"] == "BND")
        self.assertEqual((bnd["mateChrom"], bnd["matePos"]), ("chr5", 1000))
        self.assertEqual(bnd["pairedReads"], "25,8")

    def test_outside_genes_ignored(self):
        self.assertFalse(any(s["chrom"] == "chr2" for s in self.parse(DELLY)))

    def test_manta_fields_unchanged(self):
        manta = ("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                 "chr17\t43090000\tMantaDEL:1\tN\t<DEL>\t200\tPASS\tSVTYPE=DEL;END=43096000;SVLEN=-6000\tGT:PR:SR\t0/1:20,8:15,6\n")
        sv = self.parse(manta)[0]
        self.assertEqual((sv["caller"], sv["pairedReads"], sv["splitReads"], sv["length"]), ("Manta", "20,8", "15,6", 6000))


if __name__ == "__main__":
    unittest.main()


class FilterWindows(unittest.TestCase):
    def run_filter(self, records):
        import subprocess
        with tempfile.TemporaryDirectory() as d:
            bed = os.path.join(d, "sv.bed")
            with open(bed, "w") as handle:
                handle.write("chr13\t31315076\t33400268\tBRCA2\nchr17\t42044294\t44170327\tBRCA1\n")
            text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "".join(records)
            out = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "brca12_filter_sv.py"), "--windows", bed],
                                 input=text, capture_output=True, text=True, check=True).stdout
        return [l.split("\t")[2] for l in out.splitlines() if not l.startswith("#")]

    def test_window_and_length_rules(self):
        kept = self.run_filter([
            "chr17\t43070057\tDEL1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=43075601\n",       # kept
            "chr13\t32386567\tDEL2\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=32386598\n",       # 31 bp: small variant
            "chr13\t32383134\tINV3\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=100111305\n",      # END outside the windows
            "chr17\t43110000\tBND4\tN\tN]chr5:1000]\t.\tPASS\tSVTYPE=BND;CHR2=chr5\n",   # BND: partner anywhere
            "chr2\t5000\tDEL5\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=9000\n",                # outside
        ])
        self.assertEqual(kept, ["DEL1", "BND4"])
