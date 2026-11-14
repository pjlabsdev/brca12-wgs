"""Exon copy number against a reference cohort (bin/brca12lib/exon_cnv.py, bin/brca12_exon_cnv.py)."""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))
from brca12lib import exon_cnv  # noqa: E402

PANEL_ASSET = os.path.join(ROOT, "assets", "exon_cnv_panel_1kg30x.json")


def toy_panel():
    exons = {}
    for g, chrom, base in (("BRCA1", "chr17", 43000000), ("BRCA2", "chr13", 32000000)):
        for i in range(1, 9):
            exons[f"{g}_ex{i:02d}"] = {"chrom": chrom, "start": base + i * 5000, "end": base + i * 5000 + 150,
                                       "gene": g, "label": f"Exon {i}", "median": 0.0, "sd": 0.1}
    controls = {f"c{i}": {"chrom": "chr17", "start": 1000000 * i, "end": 1000000 * i + 4000} for i in range(1, 7)}
    return {"description": "toy", "size": 10, "pad": 200, "pChange": 0.05, "minQuality": 10,
            "exons": exons, "controls": controls}


def depths(panel, factor=None):
    out = {n: 30.0 for n in panel["controls"]}
    for n in panel["exons"]:
        out[n] = 30.0 * (factor or {}).get(n, 1.0)
    return out


class ExonCnvModel(unittest.TestCase):
    def test_no_call_in_normal_sample(self):
        p = toy_panel()
        r = exon_cnv.call(p, depths(p))
        self.assertTrue(r["enabled"])
        self.assertEqual(r["calls"], [])
        self.assertEqual(len(r["exons"]), 16)

    def test_two_exon_deletion(self):
        p = toy_panel()
        r = exon_cnv.call(p, depths(p, {"BRCA1_ex03": 0.5, "BRCA1_ex04": 0.5}))
        self.assertEqual(len(r["calls"]), 1)
        c = r["calls"][0]
        self.assertEqual((c["gene"], c["type"], c["copyNumber"]), ("BRCA1", "Deletion", 1))
        self.assertEqual(c["exons"], ["BRCA1_ex03", "BRCA1_ex04"])
        self.assertGreaterEqual(c["quality"], 10)

    def test_three_exon_duplication(self):
        p = toy_panel()
        r = exon_cnv.call(p, depths(p, {f"BRCA2_ex0{i}": 1.5 for i in (5, 6, 7)}))
        self.assertEqual([(c["gene"], c["type"], c["exons"]) for c in r["calls"]],
                         [("BRCA2", "Duplikation", ["BRCA2_ex05", "BRCA2_ex06", "BRCA2_ex07"])])

    def test_too_few_controls_disables(self):
        p = toy_panel()
        d = depths(p)
        for n in list(p["controls"])[:5]:
            d[n] = None
        r = exon_cnv.call(p, d)
        self.assertFalse(r["enabled"])

    def test_mean_depths_from_samtools_lines(self):
        rows = [("chr17", 99, 104, "a"), ("chr17", 200, 202, "b")]
        lines = [f"chr17\t{p}\t{10 if p <= 102 else 20}" for p in range(100, 106)] + ["chr17\t201\t7", "chr17\t202\t9"]
        d = exon_cnv.mean_depths(lines, rows)
        self.assertAlmostEqual(d["a"], (10 * 3 + 20 * 2) / 5)
        self.assertAlmostEqual(d["b"], 8.0)

    def test_build_panel_roundtrip(self):
        p = toy_panel()
        samples = [depths(p) for _ in range(5)]
        built = exon_cnv.build_panel(samples, {n: {k: e[k] for k in ("chrom", "start", "end", "gene", "label")}
                                               for n, e in p["exons"].items()}, p["controls"], 200, 0.05, 10)
        self.assertEqual(built["size"], 5)
        self.assertAlmostEqual(built["exons"]["BRCA1_ex01"]["median"], 0.0)
        self.assertAlmostEqual(built["exons"]["BRCA1_ex01"]["sd"], 0.07)


class ExonCnvHintsAndAlleles(unittest.TestCase):
    def test_hint_tier(self):
        p = toy_panel()
        p["minQuality"] = 100   # quality is capped at 99, so every segment becomes a hint
        r = exon_cnv.call(p, depths(p, {"BRCA1_ex03": 0.5, "BRCA1_ex04": 0.5}))
        self.assertEqual(r["calls"], [])
        self.assertEqual([h["exons"] for h in r["hints"]], [["BRCA1_ex03", "BRCA1_ex04"]])

    def test_allele_counts_filter_and_window(self):
        p = toy_panel()
        e = p["exons"]["BRCA1_ex02"]
        rec = lambda pos, gt, ad, ref="A", alt="G": {"chrom": "chr17", "pos": pos, "ref": ref, "alts": [alt],
                                                    "sample": {"GT": gt, "AD": ad}}
        records = [rec(e["start"] + 10, "0/1", "20,10"),        # counted
                   rec(e["start"] - 150, "0|1", "12,12"),       # inside +/- 200 bp, counted
                   rec(e["start"] + 10, "1/1", "0,30"),         # homozygous, skipped
                   rec(e["start"] + 10, "0/1", "29,1"),         # minor allele too weak, skipped
                   rec(e["start"] + 10, "0/1", "10,10", "AT", "A"),  # indel, skipped
                   rec(e["end"] + 500, "0/1", "10,10")]         # outside the window
        self.assertEqual(exon_cnv.allele_counts(records, p), {"BRCA1_ex02": [(20, 10), (12, 12)]})

    def test_allele_balance_supports_duplication(self):
        ab = exon_cnv.allele_loglik([(20, 10)])
        self.assertGreater(ab[3], ab[2])
        ab = exon_cnv.allele_loglik([(15, 15)])
        self.assertGreater(ab[2], ab[3])
        self.assertGreater(ab[2], ab[1] + 20)   # a heterozygous site argues strongly against a deletion


class ExonCnvAssetAndCli(unittest.TestCase):
    def test_panel_asset(self):
        p = exon_cnv.load_panel(PANEL_ASSET)
        self.assertEqual(p["size"], 80)
        self.assertEqual(len(p["exons"]), 50)
        self.assertEqual(len(p["controls"]), 8)
        self.assertEqual((p["pad"], p["pChange"], p["minQuality"]), (200, 0.05, 10))

    def test_cli_bed_and_call(self):
        p = toy_panel()
        with tempfile.TemporaryDirectory() as d:
            panel_path = os.path.join(d, "panel.json")
            json.dump(p, open(panel_path, "w"))
            bed = os.path.join(d, "regions.bed")
            script = os.path.join(ROOT, "bin", "brca12_exon_cnv.py")
            subprocess.run([sys.executable, script, "--panel", panel_path, "--write-bed", bed], check=True)
            rows = [l.split("\t") for l in open(bed)]
            self.assertEqual(len(rows), 16 + 6)
            lines = []
            for chrom, start, end, name in rows:
                depth = 15 if name.strip() in ("BRCA1_ex02", "BRCA1_ex03") else 30
                lines += [f"{chrom}\t{pos}\t{depth}" for pos in range(int(start) + 1, int(end) + 1)]
            out = os.path.join(d, "exon_cnv.json")
            subprocess.run([sys.executable, script, "--panel", panel_path, "--out", out],
                           input="\n".join(lines) + "\n", text=True, check=True)
            r = json.load(open(out))
            self.assertEqual([c["exons"] for c in r["calls"]], [["BRCA1_ex02", "BRCA1_ex03"]])

    def test_cli_disabled_panel(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "exon_cnv.json")
            bed = os.path.join(d, "r.bed")
            script = os.path.join(ROOT, "bin", "brca12_exon_cnv.py")
            subprocess.run([sys.executable, script, "--panel", os.path.join(ROOT, "assets", "NO_EXON_CNV_PANEL"),
                            "--write-bed", bed, "--out", out], check=True)
            self.assertEqual(os.path.getsize(bed), 0)
            self.assertFalse(json.load(open(out))["enabled"])


if __name__ == "__main__":
    unittest.main()
